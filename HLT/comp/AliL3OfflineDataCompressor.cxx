// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

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
#include <TH2F.h>

#include "AliL3Transform.h"
#include "AliL3ModelTrack.h"
#include "AliL3Compress.h"
#include "AliL3TrackArray.h"

#include "AliL3OfflineDataCompressor.h"

#if GCCVERSION == 3
using namespace std;
#endif

//_____________________________________________________________
//
//  AliL3OfflineDataCompression
//


ClassImp(AliL3OfflineDataCompressor)

AliL3OfflineDataCompressor::AliL3OfflineDataCompressor()
{
  fMarian = kFALSE;
  fTracker=0;
}

AliL3OfflineDataCompressor::AliL3OfflineDataCompressor(Char_t *path,Bool_t keep,Bool_t writeshape,Bool_t MI) 
  : AliL3DataCompressor(path,keep,writeshape)
{
  fMarian = MI;
  fTracker=0;
}

AliL3OfflineDataCompressor::~AliL3OfflineDataCompressor()
{
  if(fTracker)
    {
      fTracker->UnloadClusters();
      delete fTracker;
    }
}


void AliL3OfflineDataCompressor::LoadData(Int_t event,Bool_t sp)
{
  //Take offline reconstructed tracks as an input.
  //In this case, no remaining clusters are written.
  
  fSinglePatch = sp;
  
  char filename[1024];
  AliKalmanTrack::SetConvConst(1000/0.299792458/AliL3Transform::GetSolenoidField());
  if(fMarian==kFALSE)
    sprintf(filename,"%s/offline/AliTPCclusters.root",fPath);
  else
    sprintf(filename,"%s/offline/AliTPCclustersMI.root",fPath);
  
  TFile *in = TFile::Open(filename);
  AliTPCParam *param=(AliTPCParam*)in->Get("75x40_100x60_150x60");

  if(fMarian==kFALSE)
    fTracker = new AliTPCtracker(param);
  else
    fTracker = new AliTPCtrackerMI(param);
  fTracker->SetEventNumber(event);
#ifdef asvversion
  fTracker->LoadClusters();
#endif  
  if(fMarian==kTRUE)
    {
      AliTPCtrackerMI *mitracker = (AliTPCtrackerMI*)fTracker;
#ifdef asvversion
      mitracker->LoadInnerSectors();
      mitracker->LoadOuterSectors();
#endif
    }
  
  const Int_t MAX=20000;
  Int_t nentr=0,i=0; TObjArray tarray(MAX);
  if(fMarian==kFALSE)
    sprintf(filename,"%s/offline/AliTPCtracks.root",fPath);
  else
    sprintf(filename,"%s/offline/AliTPCtracksMI.root",fPath);
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
  
  AliL3TrackArray *comptracks = new AliL3TrackArray("AliL3ModelTrack");
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
      
      AliL3Transform::Sector2Slice(slice,padrow,sec,row);
      Double_t par[5],xk=AliL3Transform::Row2X(padrow);
      track->PropagateTo(xk);
      track->GetExternalParameters(xk,par);
      Double_t psi = TMath::ASin(par[2]) + track->GetAlpha();
      if (psi<-TMath::Pi()) psi+=2*TMath::Pi();
      if (psi>=TMath::Pi()) psi-=2*TMath::Pi();
      Float_t pt_1=TMath::Abs(par[4]);
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
      first[1] = fcl->GetY();
      first[2] = fcl->GetZ();

      AliL3Transform::Local2Global(first,slice);
      
      AliL3ModelTrack *outtrack = (AliL3ModelTrack*)comptracks->NextTrack();
      outtrack->SetNHits(nhits);
      outtrack->SetFirstPoint(first[0],first[1],first[2]);
      outtrack->SetPt(1/pt_1);
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
	  if(!AliL3Transform::Sector2Slice(slice,padrow,sec,row))
	    exit(5);
	  xyz[0] = AliL3Transform::Row2X(padrow);
	  
	  //cout<<"Hit in slice "<<slice<<" padrow "<<padrow<<" index "<<index<<" y "<<cluster->GetY()<<" z "<<cluster->GetZ()<<endl;
	  AliL3Transform::Local2Raw(xyz,sec,row);
	  //cout<<"slice "<<slice<<" padrow "<<padrow<<" pad "<<xyz[1]<<" time "<<xyz[2]<<endl;
	  
	  if(xyz[1] < -1 || xyz[1] > AliL3Transform::GetNPads(padrow) ||
	     xyz[2] < -1 || xyz[2] > AliL3Transform::GetNTimeBins())
	    {
	      cerr<<"AliL3DataCompressor::FillOfflineData : Wrong time "<<xyz[2]<<" in slice "
		  <<slice<<" padrow "<<padrow<<endl;
	      cout<<"sector "<<sec<<" row "<<row<<endl;
	      //cout<<"Hit in slice "<<slice<<" padrow "<<padrow<<" y "<<cluster->GetY()<<" z "<<cluster->GetZ()<<endl;
	      cout<<"Track hit "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<endl;
	      exit(5);
	    }
	  
	  Float_t angle = 0;
	  AliL3Transform::Local2GlobalAngle(&angle,slice);
	  if(!outtrack->CalculateReferencePoint(angle,AliL3Transform::Row2X(padrow)))
	    {
	      cerr<<"AliL3DataCompressor::FillOfflineData : Error in crossing point calc on slice "
		  <<slice<<" row "<<padrow<<endl;
	      exit(5);
	    }
	  Float_t xyz_cross[3] = {outtrack->GetPointX(),outtrack->GetPointY(),outtrack->GetPointZ()};
	  AliL3Transform::Global2Raw(xyz_cross,sec,row);
	  /*
	    if(fabs(xyz_cross[1] - xyz[1]) > 10 ||
	    fabs(xyz_cross[2] - xyz[2]) > 10)
	    {
	    cout<<"AliL3DataCompressor::FillOfflineData : Wrong crossing slice "<<slice<<" padrow "
	    <<padrow<<" pad "<<xyz[1]<<" padhit "<<xyz_cross[1]<<" time "<<xyz[2]<<" timehit "<<xyz_cross[2]<<endl;
	    outtrack->Print();
	    exit(5);
	    }
	  */
	  //cout<<" crossing "<<xyz_cross[0]<<" "<<xyz_cross[1]<<" "<<xyz_cross[2]<<endl;
	  outtrack->SetPadHit(padrow,xyz_cross[1]);
	  outtrack->SetTimeHit(padrow,xyz_cross[2]);
	  
	  if(fWriteClusterShape)
	    {
	      Float_t angle = outtrack->GetCrossingAngle(padrow,slice);
	      outtrack->SetCrossingAngleLUT(padrow,angle);
	      outtrack->CalculateClusterWidths(padrow,kTRUE);
	      Int_t patch = AliL3Transform::GetPatch(padrow);
	      Float_t sigmaY2 = cluster->GetSigmaY2() / pow(AliL3Transform::GetPadPitchWidth(patch),2);
	      Float_t sigmaZ2 = cluster->GetSigmaZ2() / pow(AliL3Transform::GetZWidth(),2);
	      outtrack->SetCluster(padrow,xyz[1],xyz[2],clustercharge,sigmaY2,sigmaZ2,3);
	    }
	  else
	    outtrack->SetCluster(padrow,xyz[1],xyz[2],clustercharge,0,0,3);
	  totcounter++;
	  outtrack->GetClusterModel(padrow)->fSlice = slice;
	  fNusedClusters++;
	}
    }
  
  cout<<"AliL3DataCompressor::FillOfflineData : Wrote "<<totcounter<<" clusters"<<endl;
  //Write tracks to file
  AliL3Compress *comp = new AliL3Compress(-1,-1,fPath,fWriteClusterShape,fEvent);
  comp->WriteFile(comptracks);
  delete comp;
  delete comptracks;

}

void AliL3OfflineDataCompressor::WriteRemaining(Bool_t select)
{
  //Write remaining clusters (not assigned to any tracks) to file
  
  if(!fKeepRemaining)
    return;
  
  if(select)
    SelectRemainingClusters();
  
  Char_t filename[1024];
  
  if(!fSinglePatch)
    {
      cerr<<"AliL3OfflineDataCompressor::WriteRemaining : You have to modify this function when not running singlepatch"<<endl;
      return;
    }
  
  ofstream idfile;
  if(fWriteIdsToFile)
    {
      sprintf(filename,"%s/comp/remains_ids.txt",fPath);
      idfile.open(filename);
    }
  
  cout<<"Writing remaining clusters "<<endl;
  Int_t nrows = AliL3Transform::GetNRows(),sector,row,sec;
  AliTPCtracker *tracker = (AliTPCtracker*)fTracker;
  for(Int_t slice=0; slice<=35; slice++)
    {
      sprintf(filename,"%s/comp/remains_%d_%d_%d.raw",fPath,fEvent,slice,-1);
      FILE *outfile = fopen(filename,"w");
      if(!outfile)
	{
	  cerr<<"AliL3OfflineDataCompressor::WriteRemaining : Cannot open file "<<filename<<endl;
	  exit(5);
	}
      
      for(Int_t padrow=0; padrow < nrows; padrow++)
	{
	  AliL3Transform::Slice2Sector(slice,padrow,sector,row);
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
	  tracker->GetNClusters(sec,row);
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

	  Int_t size = sizeof(AliL3RemainingRow) + counter*sizeof(AliL3RemainingCluster);
	  Byte_t *data = new Byte_t[size];
	  AliL3RemainingRow *tempPt = (AliL3RemainingRow*)data;
	  
	  tempPt->fPadRow = padrow;
	  tempPt->fNClusters = counter;
	  //cout<<"Found "<<counter<<" unused out of "<<ncl<<" clusters on slice "<<slice<<" padrow "<<padrow<<endl;
	  Int_t local_counter=0;
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
	      
	      if(local_counter > counter)
		{
		  cerr<<"AliL3OfflineDataCompressor::WriterRemaining : array out of range "<<local_counter<<" "<<counter<<endl;
		  return;
		}
	      Float_t xyz[3] = {AliL3Transform::Row2X(padrow),cluster->GetY(),cluster->GetZ()};
	      AliL3Transform::Local2Raw(xyz,sector,row);
	      //cout<<"y "<<cluster->GetY()<<" z "<<cluster->GetZ()<<" pad "<<xyz[1]<<" time "<<xyz[2]<<endl;
	      tempPt->fClusters[local_counter].fY = cluster->GetY();
	      tempPt->fClusters[local_counter].fZ = cluster->GetZ();
	      tempPt->fClusters[local_counter].fCharge = (UShort_t)cluster->GetQ();
	      tempPt->fClusters[local_counter].fSigmaY2 = cluster->GetSigmaY2();
	      tempPt->fClusters[local_counter].fSigmaZ2 = cluster->GetSigmaZ2();
	      local_counter++;
	      fNunusedClusters++;
	    }
	  
	  fwrite(tempPt,size,1,outfile);
	  delete [] data;
	}
      fclose(outfile);
    }
  if(fWriteIdsToFile)
    {
      idfile << endl;
      idfile.close();
    }
}

void AliL3OfflineDataCompressor::SelectRemainingClusters()
{
  
  cout<<"Cleaning up clusters"<<endl;
  Int_t nrows = AliL3Transform::GetNRows();
  Int_t gap=(Int_t)(0.125*nrows), shift=(Int_t)(0.5*gap);
  
  Int_t sector,row,sec;
  AliTPCtracker *tracker = (AliTPCtracker*)fTracker;
  for(Int_t slice=0; slice<36; slice++)
    {
      for(Int_t padrow=0; padrow < nrows; padrow++)
	{
	  if(padrow >= nrows-1-gap-shift) continue;//save all the clusters in this region
	  
	  AliL3Transform::Slice2Sector(slice,padrow,sector,row);
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
	      AliTPCcluster *cluster = 0; //todo consti (AliTPCcluster*)tracker->GetCluster(sec,row,j);
	      if(cluster->GetZ() < 0 && slice < 18) continue;
	      if(cluster->GetZ() > 0 && slice > 17) continue;
	      if(cluster->IsUsed())
		continue;
	      
	      cluster->Use();
	    }
	}
      
    }
}
