//$Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ASV

#include "AliL3StandardIncludes.h"
#include "AliL3DataCompressor.h"
#include "AliL3FileHandler.h"
#include "AliL3Transform.h"
#include "AliL3SpacePointData.h"
#include "AliL3Compress.h"
#include "AliL3TrackArray.h"
#include "AliL3ModelTrack.h"
#include "AliL3Modeller.h"

#include <AliTPCParamSR.h>
#include <AliTPCDigitsArray.h>
#include <AliTPCClustersArray.h>
#include <AliTPCcluster.h>
#include <AliTPCClustersRow.h>
#include <AliSimDigits.h>
#include <AliTPC.h>
#include <AliTPCv2.h>
#include <AliRun.h>

#include <TFile.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TSystem.h>
//_____________________________________________________________
//
//  AliL3DataCompression
//
// Interface class; binary <-> AliROOT handling of TPC data compression classes.
//


ClassImp(AliL3DataCompressor)

AliL3DataCompressor::AliL3DataCompressor()
{
  fMemHandler=0;
  fMinSlice=0;
  fMaxSlice=0;
}

AliL3DataCompressor::AliL3DataCompressor(Char_t *path,Int_t minslice,Int_t maxslice)
{
  fMinSlice=minslice;
  fMaxSlice=maxslice;
  strcpy(fPath,path);
  fMemHandler = new AliL3MemHandler();
}

AliL3DataCompressor::~AliL3DataCompressor()
{
  if(fMemHandler)
    delete fMemHandler;
}

void AliL3DataCompressor::ProcessData(Char_t *trackpath,Int_t padoverlap,Int_t timeoverlap,Int_t padsearch,Int_t timesearch)
{
  //Find clusters based on the input tracks, and write them to file.
  
  for(Int_t slice=fMinSlice; slice<=fMaxSlice; slice++)
    {
      for(Int_t patch=0; patch<AliL3Transform::GetNPatches(); patch++)
	{
	  AliL3Modeller *modeller = new AliL3Modeller();

	  modeller->SetOverlap(padoverlap,timeoverlap);
	  modeller->SetSearchRange(padsearch,timesearch);
	  modeller->Init(slice,patch,trackpath,fPath,kTRUE,kTRUE);
	  modeller->FindClusters();
	  modeller->WriteRemaining();
	  
	  AliL3TrackArray *tracks = modeller->GetTracks();

	  AliL3Compress *comp = new AliL3Compress(slice,patch,fPath);
	  
	  comp->WriteFile(tracks);
	  
	  delete comp;
	  delete modeller;
	}
    }
  
}

void AliL3DataCompressor::CompressAndExpand(Int_t bitspad,Int_t bitstime,Int_t bitscharge,Int_t bitsshape)
{
  //Read tracks/clusters from file, compress data and uncompress it. Write compression rates to file.
  //Input parameters are number of bits to use to code the pad/time/charge/shape residuals.
  
  Char_t filename[1024];
  sprintf(filename,"%s/comprates.txt",fPath);
  FILE *file = fopen(filename,"w");
  for(Int_t slice=fMinSlice; slice<=fMaxSlice; slice++)
    {
      for(Int_t patch=0; patch < AliL3Transform::GetNPatches(); patch++)
	{
	  AliL3Compress *comp = new AliL3Compress(slice,patch,fPath);
	  comp->SetBitNumbers(bitspad,bitstime,bitscharge,bitsshape);

	  comp->CompressFile();
	  comp->ExpandFile();
	  comp->PrintCompRatio(file);
	  
	  delete comp;
	}
    }
  fclose(file);
}

void AliL3DataCompressor::WriteRemainingDigits()
{
  //Write the remaining digits resulting from ProcessData to a rootfile
  //which serves as an input to the offline cluster finder.
  
  Char_t filename[1024];
  sprintf(filename,"rm %s/comp/remains.root",fPath);
  gSystem->Exec(filename);
  sprintf(filename,"%s/comp/remains.root",fPath);
  
  cout<<"AliL3DataCompressor::WriteRemainingDigits : Removing old file : "<<filename<<endl;
  for(Int_t slice=fMinSlice; slice<=fMaxSlice; slice++)
    {
      for(Int_t patch=0; patch<AliL3Transform::GetNPatches(); patch++)
	{
	  AliL3Compress *comp = new AliL3Compress(slice,patch,fPath);
	  comp->WriteRootFile(filename);
	  delete comp;
	}
    }
  
}

void AliL3DataCompressor::FindOfflineClusters(Bool_t remains)
{
  //Code taken from macro AliTPCFindClusters.
  //remains = kTRUE : Find remaining clusters after comression
  //remains = kFALSE : Find offline clusters before compression.

  Char_t fname[1024];
  if(!remains)
    sprintf(fname,"%s/AliTPCclusters.root",fPath);
  else
    sprintf(fname,"%s/comp/AliTPCclusters_remains.root",fPath);
  
  TFile *out=TFile::Open(fname,"RECREATE");
    
  sprintf(fname,"%s/alirunfile.root",fPath);
  TFile *in = TFile::Open(fname);
  
  if (!(gAlice=(AliRun*)in->Get("gAlice"))) 
    {
      cerr<<"AliL3DataCompressor::FindOfflineClusters : gAlice have not been found on file "<<fname<<endl;
      return;
   }
  
  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC"); 
  Int_t ver = TPC->IsVersion(); 
  cerr<<"TPC version "<<ver<<" has been found !\n";
  
  AliTPCParamSR *dig=(AliTPCParamSR *)in->Get("75x40_100x60");
  if(dig){
    cerr<<"2 pad-length geom hits with 3 pad-lengths geom digits\n";
    delete dig;
    dig = new AliTPCParamSR();
  }
  else
    {
      dig=(AliTPCParamSR *)gDirectory->Get("75x40_100x60_150x60");
    }
  if (!dig) {cerr<<"TPC parameters have not been found !\n"; return;}
  
  if(!remains)
    sprintf(fname,"%s/digitfile.root",fPath);
  else
    sprintf(fname,"%s/comp/remains.root",fPath);
  TFile *dfile = TFile::Open(fname);
  
  cerr<<"Looking for clusters...\n";
  AliTPCv2 tpc; 
  tpc.SetParam(dig); 
  dfile->cd();  
  
  tpc.Digits2Clusters(out,0); //event 0
  
  
  delete gAlice; gAlice=0;
  out->Close();
  in->Close();
}

void AliL3DataCompressor::WriteRemainingClusters()
{
  //Write the output from the offline cluster finder into 
  //binary files.

  Char_t filename[1024];
  sprintf(filename,"%s/comp/AliTPCclusters_remains.root",fPath);
  AliL3FileHandler file;
  file.SetAliInput(filename);
  
  for(Int_t slice=fMinSlice; slice<=fMaxSlice; slice++)
    {
      for(Int_t patch=0; patch < AliL3Transform::GetNPatches(); patch++)
	{
	  file.Init(slice,patch);
	  sprintf(filename,"%s/comp/clusters_remains_%d_%d.raw",fPath,slice,patch);
	  file.SetBinaryOutput(filename);
	  file.AliPoints2Binary();
	  file.CloseBinaryOutput();
	}
    }
  file.CloseAliInput();
}

void AliL3DataCompressor::RestoreData()
{
  //Restore the uncompressed data together with the remaining clusters,
  //and write to a final cluster file which serves as an input to the
  //final offline tracker. 
  
  Char_t filename[1024];
  
  sprintf(filename,"%s/digitfile.root",fPath);
  TFile *rootfile = TFile::Open(filename);
  rootfile->cd();
  AliTPCParam *param = (AliTPCParam*)rootfile->Get(AliL3Transform::GetParamName());
  AliTPCDigitsArray *darray = new AliTPCDigitsArray();
  darray->Setup(param);
  darray->SetClass("AliSimDigits");
  sprintf(filename,"TreeD_%s_0",AliL3Transform::GetParamName());
  Bool_t ok = darray->ConnectTree(filename);
  if(!ok)
    {
      cerr<<"AliL3DataCompressor::RestoreData : Problems connecting tree"<<endl;
      return;
    }
  
  TDirectory *savedir = gDirectory;
  
  //Create new file for storing the final clusters:
  sprintf(filename,"%s/comp/AliTPCclusters.root",fPath);
  TFile *clusterfile = TFile::Open(filename,"RECREATE");
  clusterfile->cd();
  param->Write(param->GetTitle());
  
  AliTPCClustersArray carray;
  carray.Setup(param);
  carray.SetClusterType("AliTPCcluster");
  carray.MakeTree();

  //Now start the loop:
  for(Int_t slice=fMinSlice; slice<=fMaxSlice; slice++)
    {
      for(Int_t patch=0; patch < AliL3Transform::GetNPatches(); patch++)
	{
	  //Get the remaining clusters:
	  Int_t nclusters = FindRemaining(slice,patch);
	  UInt_t size=0;
	  AliL3SpacePointData *points = (AliL3SpacePointData*)fMemHandler->GetDataPointer(size);
	  Int_t counter=0;
	  cout<<"Found "<<nclusters<<" clusters in slice "<<slice<<" patch "<<patch<<endl;
	  
	  //Get the uncompressed data:
	  
	  AliL3Compress *comp = new AliL3Compress(slice,patch,fPath);

	  comp->ReadFile('u');

	  AliL3TrackArray *tracks = comp->GetTracks();
	  for(Int_t padrow=AliL3Transform::GetFirstRow(patch); padrow<=AliL3Transform::GetLastRow(patch); padrow++)
	    {
	      Int_t sec,row;
	      AliL3Transform::Slice2Sector(slice,padrow,sec,row);
	      AliTPCClustersRow *clrow=carray.CreateRow(sec,row);
	      AliSimDigits *digits = (AliSimDigits*)darray->LoadRow(sec,row);
	      
	      Float_t pad,time,xywidth,zwidth;
	      Int_t charge;
	      for(Int_t i=0; i<tracks->GetNTracks(); i++)
		{
		  AliL3ModelTrack *track = (AliL3ModelTrack*)tracks->GetCheckedTrack(i);
		  if(!track) continue;
		  if(!track->IsPresent(padrow)) continue;
		  track->GetPad(padrow,pad);
		  track->GetTime(padrow,time);
		  track->GetClusterCharge(padrow,charge);

		  Float_t xyz[3];
		  AliL3Transform::Raw2Local(xyz,sec,row,pad,time);
		  
		  
		  //Fill the remaining clusters first, if the order is such;
		  while(counter < nclusters && points[counter].fPadRow == padrow && points[counter].fY <= xyz[1])
		    {
		      if(points[counter].fY == xyz[1] && points[counter].fZ < xyz[2]) break;//same pad
		      AliTPCcluster *c = new AliTPCcluster();
		      c->SetY(points[counter].fY);
		      c->SetZ(points[counter].fZ);
		      c->SetQ(points[counter].fCharge);
		      
		      xywidth = points[counter].fXYErr * points[counter].fXYErr;
		      zwidth = points[counter].fZErr * points[counter].fZErr;
		      c->SetSigmaY2(xywidth);
		      c->SetSigmaZ2(zwidth);
		      
		      Float_t temp[3] = {points[counter].fX,points[counter].fY,points[counter].fZ};
		      AliL3Transform::Local2Raw(temp,sec,row);
		      Int_t tpad,ttime;
		      tpad = TMath::Nint(temp[1]);
		      ttime = TMath::Nint(temp[2]);
		      c->SetLabel(digits->GetTrackID(ttime,tpad,0),0);
		      c->SetLabel(digits->GetTrackID(ttime,tpad,1),1);
		      c->SetLabel(digits->GetTrackID(ttime,tpad,2),2);

		      clrow->InsertCluster(c);
		      delete c;
		      if(counter == nclusters - 1) break;
		      counter++;
		    }
		  
		  //Now, fill the track data:
		  track->GetXYWidth(padrow,xywidth);
		  track->GetZWidth(padrow,zwidth);

		  Int_t tpad,ttime;
		  tpad = TMath::Nint(pad);
		  ttime = TMath::Nint(time);
		  if(tpad < 0 || tpad >= AliL3Transform::GetNPads(padrow) ||
		     ttime < 0 || ttime >= AliL3Transform::GetNTimeBins())
		    {
		      cerr<<"AliL3DataCompressor::RestoreData : Wrong pad "<<tpad<<" or time "<<ttime<<endl;
		      track->Print();
		      return;
		    }
		  
		  xywidth = (xywidth+1./12)*pow(AliL3Transform::GetPadPitchWidth(patch),2);
		  zwidth = (zwidth+1./12)*pow(AliL3Transform::GetZWidth(),2);
		  
		  AliTPCcluster *c = new AliTPCcluster();
		  c->SetY(xyz[1]);
		  c->SetZ(xyz[2]);
		  c->SetSigmaY2(xywidth);
		  c->SetSigmaZ2(zwidth);
		  c->SetQ(charge);
		  
		  
		  c->SetLabel(digits->GetTrackID(ttime,tpad,0),0);
		  c->SetLabel(digits->GetTrackID(ttime,tpad,1),1);
		  c->SetLabel(digits->GetTrackID(ttime,tpad,2),2);
		  		  
		  clrow->InsertCluster(c);
		  delete c;
		}

	      //Fill the rest of the remaining clusters on this row:
	      while(counter < nclusters && points[counter].fPadRow == padrow)
		{
		  AliTPCcluster *c = new AliTPCcluster();
		  c->SetY(points[counter].fY);
		  c->SetZ(points[counter].fZ);
		  xywidth = points[counter].fXYErr * points[counter].fXYErr;
		  zwidth = points[counter].fZErr * points[counter].fZErr;
		  
		  Float_t temp[3] = {points[counter].fX,points[counter].fY,points[counter].fZ};
		  AliL3Transform::Local2Raw(temp,sec,row);

		  Int_t tpad,ttime;
		  tpad = TMath::Nint(temp[1]);
		  ttime = TMath::Nint(temp[2]);
		  c->SetSigmaY2(xywidth);
		  c->SetSigmaZ2(zwidth);
		  c->SetQ(points[counter].fCharge);
		  
		  c->SetLabel(digits->GetTrackID(ttime,tpad,0),0);
		  c->SetLabel(digits->GetTrackID(ttime,tpad,1),1);
		  c->SetLabel(digits->GetTrackID(ttime,tpad,2),2);
		  
		  clrow->InsertCluster(c);
		  delete c;
		  counter++;
		}

	      carray.StoreRow(sec,row);
	      carray.ClearRow(sec,row);
	      darray->ClearRow(sec,row);
	    }
	  delete comp;
	  fMemHandler->Free();
	  if(nclusters != counter)
	    cerr<<"AliL3DataCompressor::RestoreData : Mismatching counters : "<<nclusters<<" "<<counter<<endl;
	}
    }
  
  //Write the tree to file:
  sprintf(filename,"TreeC_TPC_%d",0);
  carray.GetTree()->SetName(filename);
  carray.GetTree()->Write();
  savedir->cd();
  delete darray;
  rootfile->Close();
  //  clusterfile->Close();
  
}


Int_t AliL3DataCompressor::FindRemaining(Int_t slice,Int_t patch)
{
  Char_t filename[1024];
  sprintf(filename,"%s/comp/clusters_remains_%d_%d.raw",fPath,slice,patch);
  fMemHandler->SetBinaryInput(filename);
  fMemHandler->Init(slice,patch);
  UInt_t pointsize = 10000*sizeof(AliL3SpacePointData);
  UInt_t npoints;
  AliL3SpacePointData *points = (AliL3SpacePointData*)fMemHandler->Allocate(pointsize);
  fMemHandler->Binary2Memory(npoints,points);
  
  return (Int_t)npoints;
}
