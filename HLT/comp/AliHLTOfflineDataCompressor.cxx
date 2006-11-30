// @(#) $Id$

//_____________________________________________________________
//
//  AliHLTOfflineDataCompression
//
// Class to compress data with offline tracks
// as seeds.
//
// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include <AliTPCParamSR.h>
#include <AliTPCClustersArray.h>
#include <AliTPCcluster.h>
#include <AliTPCClustersRow.h>
#include <AliTPC.h>
#include <AliTPCv2.h>
#include <AliTPCcluster.h>
#include <AliTPCtracker.h>
#include <AliTPCclusterMI.h>
#include <AliTPCtrackerMI.h>
#include <AliKalmanTrack.h>
#include <AliRun.h>

#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TH1F.h>

#include "AliHLTTransform.h"
#include "AliHLTModelTrack.h"
#include "AliHLTCompress.h"
#include "AliHLTTrackArray.h"
#include "bitio.h"
#include "AliHLTOfflineDataCompressor.h"

#if __GNUC__ == 3
using namespace std;
#endif

ClassImp(AliHLTOfflineDataCompressor)

AliHLTOfflineDataCompressor::AliHLTOfflineDataCompressor()
{ //constructor
  fMarian = kFALSE;
  fTracker=0;
}

AliHLTOfflineDataCompressor::AliHLTOfflineDataCompressor(Char_t *path,Bool_t keep,Bool_t writeshape,Bool_t MI) 
  : AliHLTDataCompressor(path,keep,writeshape)
{ //constructor
  fMarian = MI;
  fTracker=0;
}

AliHLTOfflineDataCompressor::~AliHLTOfflineDataCompressor()
{ //deconstructor
  if(fTracker)
    {
      fTracker->UnloadClusters();
      delete fTracker;
    }
}


void AliHLTOfflineDataCompressor::LoadData(Int_t event,Bool_t sp)
{
  //Take offline reconstructed tracks as an input.
  //In this case, no remaining clusters are written.
  
  if(fTracker)
    {
      fTracker->UnloadClusters();
      delete fTracker;
    }
  
  fSinglePatch = sp;
  fEvent = event;

  char filename[1024];
  //AliKalmanTrack::SetConvConst(1000/0.299792458/AliHLTTransform::GetSolenoidField());
  if(fMarian==kFALSE)
    sprintf(filename,"%s/offline/AliTPCclusters.root",fPath);
  else
    sprintf(filename,"%s/offline/AliTPCclustersMI.root",fPath);
  
  bool compressed=0;
  if(compressed)
    {
      cout<<"AliHLTOfflineDataCompressor::LoadData : Taking compressed offline files!!"<<endl;
      sprintf(filename,"%s/comp/offline/AliTPCclusters.root",fPath);
    }
  TFile *in = TFile::Open(filename);
  AliTPCParam *param=(AliTPCParam*)in->Get("75x40_100x60_150x60");

  if(fMarian==kFALSE)
    fTracker = new AliTPCtracker(param);
  else
    fTracker = new AliTPCtrackerMI(param);
  //fTracker->SetEventNumber(event);
#ifdef asvversion
  fTracker->LoadClusters();
#endif  
  if(fMarian==kTRUE)
    {
#ifdef asvversion
      AliTPCtrackerMI *mitracker = (AliTPCtrackerMI*)fTracker;
      mitracker->LoadInnerSectors();
      mitracker->LoadOuterSectors();
#endif
    }
  
  const Int_t kMAX=20000;
  Int_t nentr=0,i=0; TObjArray tarray(kMAX);
  if(fMarian==kFALSE)
    sprintf(filename,"%s/offline/AliTPCtracks.root",fPath);
  else
    sprintf(filename,"%s/offline/AliTPCtracksMI.root",fPath);
  
  if(compressed)
    sprintf(filename,"%s/comp/offline/AliTPCtracks.root",fPath);
  
  TFile *tf=TFile::Open(filename);
  
  char tname[100]; sprintf(tname,"TreeT_TPC_%d",event);
  TTree *tracktree=(TTree*)tf->Get(tname);
  
  TBranch *tbranch=tracktree->GetBranch("tracks");
  nentr=(Int_t)tracktree->GetEntries();
  AliTPCtrack *iotrack=0;

  for (i=0; i<nentr; i++) {
    iotrack=new AliTPCtrack;
    tbranch->SetAddress(&iotrack);
    tracktree->GetEvent(i);
    tarray.AddLast(iotrack);
  }   
  delete tracktree; 
  tf->Close();
  
  AliHLTTrackArray *comptracks = new AliHLTTrackArray("AliHLTModelTrack");
  cout<<"Loaded "<<nentr<<" offline tracks"<<endl;
  Int_t slice,padrow;
  Int_t totcounter=0;

  for(i=0; i<nentr; i++)
    {
      
      AliTPCtrack *track=(AliTPCtrack*)tarray.UncheckedAt(i);
      Int_t nhits = track->GetNumberOfClusters();
      Int_t idx = track->GetClusterIndex(nhits-1);
      Int_t sec=(idx&0xff000000)>>24, row=(idx&0x00ff0000)>>16;
      
      /*
	TPC sector numbering within the AliTPCtracker class:
	There are in total 18 inner sectors and 18 outer sectors.
	This means that one sector includes _both_ sides of the TPC.
	Example: sec=0 -> sector 0 and 18.
	         sec=18 -> sector 36 and 54
      */
      
      if(fMarian==kFALSE)
	if(sec >= 18)
	  sec += 18;
      
      AliHLTTransform::Sector2Slice(slice,padrow,sec,row);
      Double_t par[5],xk=AliHLTTransform::Row2X(padrow);
      track->PropagateTo(xk);
      track->GetExternalParameters(xk,par);
      Double_t psi = TMath::ASin(par[2]) + track->GetAlpha();
      if (psi<-TMath::Pi()) psi+=2*TMath::Pi();
      if (psi>=TMath::Pi()) psi-=2*TMath::Pi();
      Float_t pt1=TMath::Abs(par[4]);
      Int_t charge = 1;
      if(par[4] > 0)
	charge=-1;

      Float_t first[3];
      AliCluster *fcl=0;
      if(fMarian==kFALSE)
	fcl= fTracker->GetCluster(idx);
      else
	{
	  AliTPCtrackerMI *mitracker = (AliTPCtrackerMI*)fTracker;
	  fcl = mitracker->GetClusterMI(idx);
	}
      first[0] = xk;
      first[1] = par[0];
      first[2] = par[1];

      AliHLTTransform::Local2Global(first,slice);
      
      AliHLTModelTrack *outtrack = (AliHLTModelTrack*)comptracks->NextTrack();
      outtrack->SetNHits(nhits);
      outtrack->SetFirstPoint(first[0],first[1],first[2]);
      outtrack->SetPt(1./pt1);
      outtrack->SetPsi(psi);
      outtrack->SetTgl(par[3]);
      outtrack->SetCharge(charge);
      outtrack->CalculateHelix();
      outtrack->Init(0,-1);
      
      //cout<<"Loading track with "<<nhits<<" hits"<<endl;
      for(int j=nhits-1; j>=0; j--)
	{

	  Int_t index = track->GetClusterIndex(j);
	  if(index == 0)
	    continue;

	  Float_t xyz[3];
	  Int_t clustercharge =0;

	  //AliTPCcluster *cluster = (AliTPCcluster*)tracker->GetCluster(index);
	  AliCluster *cluster=0;
	  
	  if(fMarian==kFALSE)
	    cluster = fTracker->GetCluster(index);
	  else
	    {
	      AliTPCtrackerMI *mitracker = (AliTPCtrackerMI*)fTracker;
	      cluster = mitracker->GetClusterMI(index);
	    }

	  xyz[1] = cluster->GetY();
	  xyz[2] = cluster->GetZ();
	  if(fMarian==kFALSE)
	    {
	      AliTPCcluster *cl = (AliTPCcluster*)cluster;
	      clustercharge = (Int_t)cl->GetQ();
	    }
	  else
	    {
	      AliTPCclusterMI *cl = (AliTPCclusterMI*)cluster;
	      clustercharge = (Int_t)cl->GetQ();
	    }

	  cluster->Use();
	  
	  sec=(index&0xff000000)>>24; row=(index&0x00ff0000)>>16;
	  
	  if(fMarian==kFALSE)
	    {
	      if(sec >= 18)
		sec += 18;
	      
	      if(xyz[2] < 0)
		sec += 18;
	    }
	    
	  //cout<<"sector "<<sec<<" row "<<row<<endl;
	  if(!AliHLTTransform::Sector2Slice(slice,padrow,sec,row))
	    exit(5);
	  xyz[0] = AliHLTTransform::Row2X(padrow);
	  
	  //cout<<"Hit in slice "<<slice<<" padrow "<<padrow<<" index "<<index<<" y "<<cluster->GetY()<<" z "<<cluster->GetZ()<<endl;
	  AliHLTTransform::Local2Raw(xyz,sec,row);
	  //cout<<"slice "<<slice<<" padrow "<<padrow<<" pad "<<xyz[1]<<" time "<<xyz[2]<<endl;
	  
	  if(xyz[1] < -1 || xyz[1] > AliHLTTransform::GetNPads(padrow) ||
	     xyz[2] < -1 || xyz[2] > AliHLTTransform::GetNTimeBins())
	    {
	      cerr<<"AliHLTDataCompressor::FillOfflineData : Wrong time "<<xyz[2]<<" in slice "
		  <<slice<<" padrow "<<padrow<<endl;
	      cout<<"sector "<<sec<<" row "<<row<<endl;
	      //cout<<"Hit in slice "<<slice<<" padrow "<<padrow<<" y "<<cluster->GetY()<<" z "<<cluster->GetZ()<<endl;
	      cout<<"Track hit "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<endl;
	      exit(5);
	    }

	  Float_t angle = 0;
	  AliHLTTransform::Local2GlobalAngle(&angle,slice);
	  if(!outtrack->CalculateReferencePoint(angle,AliHLTTransform::Row2X(padrow)))
	    {
	      cerr<<"AliHLTDataCompressor::FillOfflineData : Error in crossing point calc on slice "
		  <<slice<<" row "<<padrow<<endl;
	      exit(5);
	    }

	  Float_t xyzcross[3] = {outtrack->GetPointX(),outtrack->GetPointY(),outtrack->GetPointZ()};
	  AliHLTTransform::Global2Raw(xyzcross,sec,row);
	  /*
	    if(fabs(xyzcross[1] - xyz[1]) > 10 ||
	    fabs(xyzcross[2] - xyz[2]) > 10)
	    {
	    cout<<"AliHLTDataCompressor::FillOfflineData : Wrong crossing slice "<<slice<<" padrow "
	    <<padrow<<" pad "<<xyz[1]<<" padhit "<<xyzcross[1]<<" time "<<xyz[2]<<" timehit "<<xyzcross[2]<<endl;
	    outtrack->Print();
	    exit(5);
	    }
	  */
	  //cout<<" crossing "<<xyzcross[0]<<" "<<xyzcross[1]<<" "<<xyzcross[2]<<endl;
	  outtrack->SetPadHit(padrow,xyzcross[1]);
	  outtrack->SetTimeHit(padrow,xyzcross[2]);
	  
	  if(fWriteClusterShape)
	    {
	      Float_t angle = outtrack->GetCrossingAngle(padrow,slice);
	      outtrack->SetCrossingAngleLUT(padrow,angle);
	      outtrack->CalculateClusterWidths(padrow,kTRUE);
	      Int_t patch = AliHLTTransform::GetPatch(padrow);
	      Float_t sigmaY2 = cluster->GetSigmaY2() / pow(AliHLTTransform::GetPadPitchWidth(patch),2);
	      Float_t sigmaZ2 = cluster->GetSigmaZ2() / pow(AliHLTTransform::GetZWidth(),2);
	      outtrack->SetCluster(padrow,xyz[1],xyz[2],clustercharge,sigmaY2,sigmaZ2,3);
	    }
	  else
	    outtrack->SetCluster(padrow,xyz[1],xyz[2],clustercharge,0,0,3);
	  totcounter++;
	  outtrack->GetClusterModel(padrow)->fSlice = slice;
	  fNusedClusters++;
	}
    }
  
  cout<<"AliHLTDataCompressor::FillOfflineData : Wrote "<<totcounter<<" clusters"<<endl;
  //Write tracks to file
  AliHLTCompress *comp = new AliHLTCompress(-1,-1,fPath,fWriteClusterShape,fEvent);
  comp->WriteFile(comptracks);
  delete comp;
  delete comptracks;
}

void AliHLTOfflineDataCompressor::WriteRemaining(Bool_t select)
{
  //Write remaining clusters (not assigned to any tracks) to file
  
  if(!fKeepRemaining)
    return;
  
  if(select)
    SelectRemainingClusters();
  
  Char_t filename[1024];
  
  if(!fSinglePatch)
    {
      cerr<<"AliHLTOfflineDataCompressor::WriteRemaining : You have to modify this function when not running singlepatch"<<endl;
      return;
    }
  
  ofstream idfile;
  if(fWriteIdsToFile)
    {
      sprintf(filename,"%s/comp/remains_ids.txt",fPath);
      idfile.open(filename);
    }
  
  cout<<"Writing remaining clusters "<<endl;
  Int_t nrows = AliHLTTransform::GetNRows(),sector,row,sec;
#ifdef asvversion
  AliTPCtracker *tracker = (AliTPCtracker*)fTracker;
#endif
  for(Int_t slice=0; slice<=35; slice++)
    {
      sprintf(filename,"%s/comp/remains_%d_%d_%d.raw",fPath,fEvent,slice,-1);
      BIT_FILE *output = OpenOutputBitFile(filename);
      if(!output)
	{
	  cerr<<"AliHLTOfflineDataCompressor::WriteRemaining : Cannot open file "<<filename<<endl;
	  exit(5);
	}
      
      //Write number of padrows with clusters
      OutputBits(output,nrows,8);
      
      for(Int_t padrow=0; padrow < nrows; padrow++)
	{
	  AliHLTTransform::Slice2Sector(slice,padrow,sector,row);
	  sec=sector;
	  
	  if(fMarian == kFALSE)
	    {
	      if(slice >= 18)
		sec -= 18;
	      if(sec >= 18)
		sec -= 18;
	    }
	  //cout<<"Getting clusters in sector "<<sec<<" row "<<row<<endl;
	  Int_t ncl = 0;
#ifdef asvversion
	  ncl = tracker->GetNClusters(sec,row);
#endif
	  	  
	  Int_t counter=0;
	  Int_t j;
	  for(j=0; j<ncl; j++)
	    {
	      AliTPCcluster *cluster = 0;
#ifdef asvversion
	      cluster=(AliTPCcluster*)tracker->GetCluster(sec,row,j);
#endif
	      if(cluster->GetZ() < 0 && slice < 18) continue;
	      if(cluster->GetZ() > 0 && slice > 17) continue;
	      if(cluster->IsUsed())
		continue;
	      counter++;
	    }
	  
	  OutputBits(output,padrow,8);//Write padrow #
	  OutputBits(output,counter,10);//Write number of clusters on this padrow

	  //cout<<"Found "<<counter<<" unused out of "<<ncl<<" clusters on slice "<<slice<<" padrow "<<padrow<<endl;
	  for(j=0; j<ncl; j++)
	    {
	      AliTPCcluster *cluster = 0;
#ifdef asvversion
	      cluster=(AliTPCcluster*)tracker->GetCluster(sec,row,j);
#endif
	      if(cluster->GetZ() < 0 && slice < 18) continue;
	      if(cluster->GetZ() > 0 && slice > 17) continue;
	      if(cluster->IsUsed())
		continue;
	      
	      if(fWriteIdsToFile)
		idfile << cluster->GetLabel(0)<<' ';
	      
	      Float_t xyz[3] = {AliHLTTransform::Row2X(padrow),cluster->GetY(),cluster->GetZ()};
	      AliHLTTransform::Local2Raw(xyz,sector,row);

	      Int_t patch = AliHLTTransform::GetPatch(padrow);
	      Float_t padw = cluster->GetSigmaY2()/pow(AliHLTTransform::GetPadPitchWidth(patch),2);
	      Float_t timew = cluster->GetSigmaZ2()/pow(AliHLTTransform::GetZWidth(),2);
	      
	      Int_t buff;
	      //Write pad
	      buff = (Int_t)rint(xyz[1]*10);
	      if(buff<0)
		{
		  cerr<<"AliHLTOfflineDataCompressor:WriteRemaining : Wrong pad value "<<buff<<endl;
		  exit(5);
		}
	      OutputBits(output,buff,11);

	      //Write time
	      buff = (Int_t)rint(xyz[2]*10);
	      if(buff<0)
		{
		  cerr<<"AliHLTOfflineDataCompressor:WriteRemaining : Wrong time value "<<buff<<endl;
		  exit(5);
		}
	      OutputBits(output,buff,13);

	      //Write widths
	      buff = (Int_t)rint(padw*100);
	      OutputBits(output,buff,8);
	      buff = (Int_t)rint(timew*100);
	      OutputBits(output,buff,8);
	      
	      //Write charge 
	      buff = (Int_t)cluster->GetQ();
	      if(buff >= 1<<14)
		buff = (1<<14)-1;
	      OutputBits(output,buff,14);
	      
	      fNunusedClusters++;
	    }
	}
      CloseOutputBitFile(output);
    }
  if(fWriteIdsToFile)
    {
      idfile << endl;
      idfile.close();
    }
}

void AliHLTOfflineDataCompressor::SelectRemainingClusters()
{
  //select the remaining clusters 
  //which were not compressed  
  cout<<"Cleaning up clusters"<<endl;
  Int_t nrows = AliHLTTransform::GetNRows();
  Int_t gap=(Int_t)(0.125*nrows), shift=(Int_t)(0.5*gap);
  
  Int_t sector,row,sec;

#ifdef asvversion
  AliTPCtracker *tracker = (AliTPCtracker*)fTracker;
#endif
  for(Int_t slice=0; slice<36; slice++)
    {
      for(Int_t padrow=0; padrow < nrows; padrow++)
	{
	  	  
	  AliHLTTransform::Slice2Sector(slice,padrow,sector,row);
	  sec=sector;
	  
	  if(fMarian == kFALSE)
	    {
	      if(slice >= 18)
		sec -= 18;
	      if(sec >= 18)
		sec -= 18;
	    }
	  Int_t ncl = 0;
#ifdef asvversion
	  ncl=tracker->GetNClusters(sec,row);
#endif
	  for(Int_t j=0; j<ncl; j++)
	    {
	      AliTPCcluster *cluster = 0;
#ifdef asvversion
	      cluster = (AliTPCcluster*)tracker->GetCluster(sec,row,j);
#endif
	      if(cluster->IsUsed())
		continue;

	      //Check the widths (errors) of the cluster, and remove big bastards:
	      Float_t xyw = cluster->GetSigmaY2() / pow(AliHLTTransform::GetPadPitchWidth(AliHLTTransform::GetPatch(padrow)),2);
	      Float_t zw  = cluster->GetSigmaZ2() / pow(AliHLTTransform::GetZWidth(),2);
	      if(xyw >= 2.55 || zw >= 2.55)//Because we use 1 byte to store
		{
		  cluster->Use();
		  continue;
		}
	      
	      //if(padrow >= nrows-1-gap-shift) continue;//save all the clusters in this region
	      
	      if(padrow == nrows - 1 || padrow == nrows - 1 - gap ||                 //First seeding
		 padrow == nrows - 1 - shift || padrow == nrows - 1 - gap - shift)   //Second seeding
		continue;
	      
	      if(cluster->GetZ() < 0 && slice < 18) continue;
	      if(cluster->GetZ() > 0 && slice > 17) continue;
	      if(cluster->IsUsed())
		continue;
	      
	      cluster->Use();
	    }
	}
      
    }
}
