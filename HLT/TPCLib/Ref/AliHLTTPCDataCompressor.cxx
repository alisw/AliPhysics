// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTTPCStandardIncludes.h"

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCCompress.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCModelTrack.h"
#include "AliHLTTPCBenchmark.h"
#include "AliHLTTPCClusterFitter.h"

#ifdef use_aliroot
#include "AliHLTTPCFileHandler.h"
#include <AliTPCcluster.h>
#include <AliTPCParamSR.h>
#include <AliTPCDigitsArray.h>
#include <AliTPCClustersArray.h>
#include <AliTPCClustersRow.h>
#include <AliSimDigits.h>
#include <AliTPC.h>
#include <AliTPCv2.h>
#include <AliRun.h>
#endif

#ifdef use_root
#include <TFile.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TH2F.h>
#endif

#include "AliHLTTPCDataCompressor.h"

#if GCCVERSION == 3
using namespace std;
#endif

//_____________________________________________________________
//
//  AliHLTTPCDataCompression
//
// Interface class; binary <-> AliROOT handling of TPC data compression classes.
//


ClassImp(AliHLTTPCDataCompressor)

Int_t AliHLTTPCDataCompressor::fNumTimeBits = 12;
Int_t AliHLTTPCDataCompressor::fNumPadBits = 12;
Int_t AliHLTTPCDataCompressor::fNumChargeBits = 14;
Int_t AliHLTTPCDataCompressor::fNumShapeBits = 14;
Float_t AliHLTTPCDataCompressor::fXYResidualStep1 = 0.03;
Float_t AliHLTTPCDataCompressor::fXYResidualStep2 = 0.03;
Float_t AliHLTTPCDataCompressor::fXYResidualStep3 = 0.03;
Float_t AliHLTTPCDataCompressor::fZResidualStep1 = 0.05;
Float_t AliHLTTPCDataCompressor::fZResidualStep2 = 0.05;
Float_t AliHLTTPCDataCompressor::fZResidualStep3 = 0.05;
Float_t AliHLTTPCDataCompressor::fXYWidthStep = 0.005;
Float_t AliHLTTPCDataCompressor::fZWidthStep = 0.005;
Int_t AliHLTTPCDataCompressor::fClusterCharge = 100;

AliHLTTPCDataCompressor::AliHLTTPCDataCompressor()
{
  fBenchmark=0;
  fInputTracks=0;
  fKeepRemaining=kTRUE;
  fEvent=0;
  fWriteClusterShape=kFALSE;
  fOutputFile=0;
  fCompRatioFile=0;
  fNusedClusters=0;
  fNunusedClusters=0;
  memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
}

AliHLTTPCDataCompressor::AliHLTTPCDataCompressor(Char_t *path,Bool_t keep,Bool_t writeshape)
{
  strcpy(fPath,path);
  fBenchmark = new AliHLTTPCBenchmark();
  fInputTracks=0;
  fKeepRemaining=keep;
  fWriteClusterShape = writeshape;
  fEvent=0;
  fOutputFile=0;
  fNusedClusters=0;
  fNunusedClusters=0;
  memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
#ifdef use_root
  Char_t name[1024];
  sprintf(name,"rm -f %s/comp/*",path);//Clean the directory
  gSystem->Exec(name);
#endif
  OpenOutputFile();
}

AliHLTTPCDataCompressor::~AliHLTTPCDataCompressor()
{
  if(fInputTracks)
    delete fInputTracks;
  if(fBenchmark)
    delete fBenchmark;
  if(fClusters)
    {
      for(Int_t i=0; i<36; i++)
	for(Int_t j=0; j<6; j++)
	  if(fClusters[i][j])
	    delete fClusters[i][j];
    }
  CloseOutputFile();
}

void AliHLTTPCDataCompressor::DoBench(Char_t *fname)
{
  fBenchmark->Analyze(fname);
}

void AliHLTTPCDataCompressor::SetBitNumbers(Int_t pad,Int_t time,Int_t charge,Int_t shape)
{
  fNumPadBits = pad;
  fNumTimeBits = time;
  fNumChargeBits = charge;
  fNumShapeBits = shape;
}

void AliHLTTPCDataCompressor::SetTransverseResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width)
{
  fXYResidualStep1 = res1;
  fXYResidualStep2 = res2;
  fXYResidualStep3 = res3;
  fXYWidthStep = width;
}

void AliHLTTPCDataCompressor::SetLongitudinalResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width)
{
  fZResidualStep1 = res1;
  fZResidualStep2 = res2;
  fZResidualStep3 = res3;
  fZWidthStep = width;
}

const Float_t AliHLTTPCDataCompressor::GetXYResidualStep(Int_t row) 
{
  if(row < AliHLTTPCTransform::GetNRowLow())
    return fXYResidualStep1;
  else if(row < AliHLTTPCTransform::GetNRowLow() + AliHLTTPCTransform::GetNRowUp1())
    return fXYResidualStep2;
  else if(row < AliHLTTPCTransform::GetNRowLow() + AliHLTTPCTransform::GetNRowUp1() + AliHLTTPCTransform::GetNRowUp2())
    return fXYResidualStep3;
  else
    {
      cerr<<"AliHLTTPCDataCompressor::GetXYResidualStep : Wrong row number "<<row<<endl;
      return -1;
    }
}

const Float_t AliHLTTPCDataCompressor::GetZResidualStep(Int_t row) 
{
  if(row < AliHLTTPCTransform::GetNRowLow())
    return fZResidualStep1;
  else if(row < AliHLTTPCTransform::GetNRowLow() + AliHLTTPCTransform::GetNRowUp1())
    return fZResidualStep2;
  else if(row < AliHLTTPCTransform::GetNRowLow() + AliHLTTPCTransform::GetNRowUp1() + AliHLTTPCTransform::GetNRowUp2())
    return fZResidualStep3;
  else
    {
      cerr<<"AliHLTTPCDataCompressor::GetXYResidualStep : Wrong row number "<<row<<endl;
      return -1;
    }
}

void AliHLTTPCDataCompressor::OpenOutputFile()
{
#ifndef use_aliroot
   LOG(AliHLTTPCLog::kError,"AliHLTTPCDataCompressor::OpenOutputFile","Version")
     <<"You have to compile with use_aliroot flag in order to use this function"<<ENDLOG;
#else
  Char_t filename[1024];
  
  sprintf(filename,"%s/comp/comprates.txt",fPath);
  fCompRatioFile = new ofstream(filename);
  
  if(fOutputFile)
    if(fOutputFile->IsOpen())
      fOutputFile->Close();

  sprintf(filename,"%s/alirunfile.root",fPath);
  TFile *f = TFile::Open(filename);
  AliTPCParam *param = (AliTPCParam*)f->Get(AliHLTTPCTransform::GetParamName());
  sprintf(filename,"%s/comp/AliTPCclusters.root",fPath);
  fOutputFile = TFile::Open(filename,"RECREATE");
  param->Write(param->GetTitle());
  f->Close();
#endif
}

void AliHLTTPCDataCompressor::CloseOutputFile()
{
  if(fCompRatioFile)
    {
      fCompRatioFile->close();
      delete fCompRatioFile;
    }
  
  if(!fOutputFile)
    return;
#ifdef use_root
  if(!fOutputFile->IsOpen())
    return;
  fOutputFile->Close();
#else
  fclose(fOutputFile);
#endif
  fOutputFile=0;
}

void AliHLTTPCDataCompressor::LoadData(Int_t event,Bool_t sp)
{
  fSinglePatch=sp;
  fEvent=event;
  AliHLTTPCMemHandler *clusterfile[36][6];
  Char_t fname[1024];
  for(Int_t s=0; s<=35; s++)
    {
      for(Int_t p=0; p<6; p++)
	{
	  if(fClusters[s][p])
	    delete fClusters[s][p];
	  fClusters[s][p] = 0;
	  clusterfile[s][p] = new AliHLTTPCMemHandler();
	  if(fSinglePatch)
	    sprintf(fname,"%s/cf/points_%d_%d_%d.raw",fPath,fEvent,s,-1);
	  else
	    sprintf(fname,"%s/cf/points_%d_%d_%d.raw",fPath,fEvent,s,p);
	  clusterfile[s][p]->SetBinaryInput(fname);
	  
	  fClusters[s][p] = (AliHLTTPCSpacePointData*)clusterfile[s][p]->Allocate();
	  clusterfile[s][p]->Binary2Memory(fNcl[s][p],fClusters[s][p]);
	  clusterfile[s][p]->CloseBinaryInput();
	  
	  if(fSinglePatch)
	    break;
	}
    }
  
  sprintf(fname,"%s/cf/tracks_%d.raw",fPath,fEvent);
  AliHLTTPCMemHandler *tfile = new AliHLTTPCMemHandler();
  tfile->SetBinaryInput(fname);
  
  if(fInputTracks)
    delete fInputTracks;
  fInputTracks = new AliHLTTPCTrackArray();
  tfile->Binary2TrackArray(fInputTracks);
  tfile->CloseBinaryInput();
  delete tfile;
}

void AliHLTTPCDataCompressor::FillData(Int_t min_hits,Bool_t expand)
{
  
  //Fill the track data into track and cluster structures, and write to file.
  //Preparation for compressing it.
  
#if 0
  cout<<"Filling data; "<<fInputTracks->GetNTracks()<<" tracks"<<endl;
#endif
  AliHLTTPCTrackArray *comptracks = new AliHLTTPCTrackArray("AliHLTTPCModelTrack");
  fInputTracks->QSort();
  for(Int_t i=0; i<fInputTracks->GetNTracks(); i++)
    {
      AliHLTTPCTrack *intrack = fInputTracks->GetCheckedTrack(i);
      if(!intrack) continue;

      if(intrack->GetNHits()<min_hits) break;

      intrack->CalculateHelix();
      
      AliHLTTPCModelTrack *outtrack = (AliHLTTPCModelTrack*)comptracks->NextTrack();
      outtrack->SetNHits(intrack->GetNHits());
      outtrack->SetRowRange(intrack->GetFirstRow(),intrack->GetLastRow());
      outtrack->SetFirstPoint(intrack->GetFirstPointX(),intrack->GetFirstPointY(),intrack->GetFirstPointZ());
      outtrack->SetLastPoint(intrack->GetLastPointX(),intrack->GetLastPointY(),intrack->GetLastPointZ());
      outtrack->SetPt(intrack->GetPt());
      outtrack->SetPsi(intrack->GetPsi());
      outtrack->SetTgl(intrack->GetTgl());
      outtrack->SetCharge(intrack->GetCharge());
      outtrack->CalculateHelix();
      Int_t nhits = intrack->GetNHits();
      UInt_t *hitids = intrack->GetHitNumbers();
      Int_t origslice = (hitids[nhits-1]>>25)&0x7f;
      outtrack->Init(origslice,-1);
      for(Int_t j=nhits-1; j>=0; j--)
	{
	  UInt_t id=hitids[j];
	  Int_t slice = (id>>25)&0x7f;
	  Int_t patch = (id>>22)&0x7;
	  UInt_t pos = id&0x3fffff;	     

	  //UInt_t size;
	  AliHLTTPCSpacePointData *points = fClusters[slice][patch];//->GetDataPointer(size);
	  Float_t xyz[3] = {points[pos].fX,points[pos].fY,points[pos].fZ};
	  Int_t padrow = points[pos].fPadRow;

	  //Calculate the crossing point between track and padrow
	  Float_t angle = 0; //Perpendicular to padrow in local coordinates
	  AliHLTTPCTransform::Local2GlobalAngle(&angle,slice);
	  if(!intrack->CalculateReferencePoint(angle,AliHLTTPCTransform::Row2X(padrow)))
	    {
	      cerr<<"AliHLTTPCDataCompressor::FillData : Error in crossing point calc on slice "<<slice<<" row "<<padrow<<endl;
	      break;
	      //outtrack->Print(kFALSE);
	      //exit(5);
	    }
	  
	  Float_t xyz_cross[3] = {intrack->GetPointX(),intrack->GetPointY(),intrack->GetPointZ()};
	  
	  Int_t sector,row;
	  AliHLTTPCTransform::Slice2Sector(slice,padrow,sector,row);
	  AliHLTTPCTransform::Global2Raw(xyz_cross,sector,row);
	  AliHLTTPCTransform::Global2Raw(xyz,sector,row);
	  
	  outtrack->SetPadHit(padrow,xyz_cross[1]);
	  outtrack->SetTimeHit(padrow,xyz_cross[2]);

	  if(fWriteClusterShape)
	    {
	      Float_t angle = intrack->GetCrossingAngle(padrow,slice);
	      outtrack->SetCrossingAngleLUT(padrow,angle);
	      outtrack->CalculateClusterWidths(padrow,kTRUE);
	      Int_t patch = AliHLTTPCTransform::GetPatch(padrow);
	      Float_t sigmaY2 = points[pos].fSigmaY2 / pow(AliHLTTPCTransform::GetPadPitchWidth(patch),2);
	      Float_t sigmaZ2 = points[pos].fSigmaZ2 / pow(AliHLTTPCTransform::GetZWidth(),2);
	      outtrack->SetCluster(padrow,xyz[1],xyz[2],points[pos].fCharge,sigmaY2,sigmaZ2,3);
	    }
	  else
	    outtrack->SetCluster(padrow,xyz[1],xyz[2],points[pos].fCharge,0,0,3);
	  
	  //IMPORTANT: Set the slice in which cluster is, you need it in AliHLTTPCModelTrack::FillTrack!
	  outtrack->GetClusterModel(padrow)->fSlice=slice;
	  points[pos].fCharge = 0;//Mark this cluster as used.
	  fNusedClusters++;
	}
      if(!expand)
	outtrack->SetNClusters(AliHLTTPCTransform::GetNRows(-1));
    }
  
  if(expand)
    ExpandTrackData(comptracks);
  
#if 0
  cout<<"Writing "<<comptracks->GetNTracks()<<" tracks to file"<<endl;
#endif
  AliHLTTPCCompress *comp = new AliHLTTPCCompress(-1,-1,fPath,fWriteClusterShape,fEvent);
  comp->WriteFile(comptracks);
  delete comp;
  delete comptracks;
  
}

void AliHLTTPCDataCompressor::ExpandTrackData(AliHLTTPCTrackArray *tracks)
{
  //Loop over tracks and try to assign unused clusters.
  //Only clusters which are closer than the max. residual are taken.
  
#if 0
  cout<<"Expanding "<<tracks->GetNTracks()<<" tracks"<<endl;
#endif
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      if(track->GetNHits() == AliHLTTPCTransform::GetNRows()) continue;
      
      Int_t nhits = track->GetNHits();
      //cout<<"Expanding track with "<<nhits<<" clusters"<<endl;
      
      Int_t last_slice=-1;
      for(Int_t padrow=AliHLTTPCTransform::GetNRows()-1; padrow>=0; padrow--)
	{
	  if(track->IsPresent(padrow))
	    {
	      last_slice = track->GetClusterModel(padrow)->fSlice;
	      continue;
	    }
	  
	  if(last_slice < 0) //the outer cluster is missing, so skip it - it will be written anyhow.
	    continue;
	  
	  //Check the slice of the next padrow:
	  Int_t next_padrow = padrow-1;
	  Int_t next_slice = -1;
	  while(next_padrow >=0)
	    {
	      if(track->IsPresent(next_padrow))
		{
		  next_slice = track->GetClusterModel(next_padrow)->fSlice;
		  break;
		}
	      next_padrow--;
	    }
	  if(next_slice>=0)
	    if(next_slice != last_slice)//The track crosses a slice boundary here
	      continue;
	  
 	  //UInt_t size;
	  AliHLTTPCSpacePointData *points = fClusters[last_slice][0];//->GetDataPointer(size);
	  
	  Float_t angle = 0;
	  AliHLTTPCTransform::Local2GlobalAngle(&angle,last_slice);
	  if(!track->CalculateReferencePoint(angle,AliHLTTPCTransform::Row2X(padrow)))
	    continue;
	  Float_t xyz_cross[3] = {track->GetPointX(),track->GetPointY(),track->GetPointZ()};
	  AliHLTTPCTransform::Global2Local(xyz_cross,last_slice,kTRUE);
	  Float_t mindist = 123456789;
	  AliHLTTPCSpacePointData *closest=0;
	  for(UInt_t j=0; j<fNcl[last_slice][0]; j++)
	    {
	      if(points[j].fCharge == 0) continue;// || points[j].fPadRow != padrow) continue;
	      if(points[j].fPadRow < padrow) continue;
	      if(points[j].fPadRow > padrow) break;
	      Float_t xyz[3] = {points[j].fX,points[j].fY,points[j].fZ};
	      AliHLTTPCTransform::Global2Local(xyz,last_slice,kTRUE);
	      
	      //Check for overflow:
	      Int_t temp = (Int_t)rint((xyz_cross[1]-xyz[1])/GetXYResidualStep(padrow));
	      if( abs(temp) > 1<<(GetNPadBits()-1))
		continue;
	      
	      temp = (Int_t)rint((xyz_cross[2]-xyz[2])/GetZResidualStep(padrow));
	      if( abs(temp) > 1<<(GetNTimeBits()-1))
		continue;
	      
	      Float_t dist = sqrt( pow(xyz_cross[1]-xyz[1],2) + pow(xyz_cross[2]-xyz[2],2) );
	      if(dist < mindist)
		{
		  closest = &points[j];
		  mindist = dist;
		}
	    }
	  if(closest) //there was a cluster assigned
	    {
	      Int_t sector,row;
	      Float_t xyz[3] = {closest->fX,closest->fY,closest->fZ};
	      AliHLTTPCTransform::Slice2Sector(last_slice,padrow,sector,row);
	      AliHLTTPCTransform::Local2Raw(xyz_cross,sector,row);
	      AliHLTTPCTransform::Global2Raw(xyz,sector,row);
	      
	      track->SetPadHit(padrow,xyz_cross[1]);
	      track->SetTimeHit(padrow,xyz_cross[2]);
	      
	      if(fWriteClusterShape)
		{
		  Float_t angle = track->GetCrossingAngle(padrow,last_slice);
		  track->SetCrossingAngleLUT(padrow,angle);
		  track->CalculateClusterWidths(padrow,kTRUE);
		  Int_t patch = AliHLTTPCTransform::GetPatch(padrow);
		  Float_t sigmaY2 = closest->fSigmaY2 / pow(AliHLTTPCTransform::GetPadPitchWidth(patch),2);
		  Float_t sigmaZ2 = closest->fSigmaZ2 / pow(AliHLTTPCTransform::GetZWidth(),2);
		  track->SetCluster(padrow,xyz[1],xyz[2],closest->fCharge,sigmaY2,sigmaZ2,3);
		}
	      else
		track->SetCluster(padrow,xyz[1],xyz[2],closest->fCharge,0,0,3);
	      nhits++;
	      
	      //IMPORTANT: Set the slice in which cluster is, you need it in AliHLTTPCModelTrack::FillTrack!
	      track->GetClusterModel(padrow)->fSlice=last_slice;
	      closest->fCharge = 0;//Mark this cluster as used.
	    }
	}
      track->SetNClusters(AliHLTTPCTransform::GetNRows());
      //cout<<"Track was assigned "<<nhits<<" clusters"<<endl;
    }
  
}

void AliHLTTPCDataCompressor::WriteRemaining(Bool_t select)
{
  //Write remaining clusters (not assigned to any tracks) to file

  
  if(!fKeepRemaining)
    return;
  
  if(select)
    SelectRemainingClusters();
  
  Char_t filename[1024];
  
  if(!fSinglePatch)
    {
      cerr<<"AliHLTTPCCompressor::WriteRemaining : You have to modify this function when not running singlepatch"<<endl;
      return;
    }

#if 0
  cout<<"Writing remaining clusters "<<endl;
#endif
  Int_t nrows = AliHLTTPCTransform::GetNRows();
  Int_t *npoints = new Int_t[nrows];
  for(Int_t i=0; i<=35; i++)
    {
      for(Int_t patch=0; patch < 1; patch++)
	{
	  sprintf(filename,"%s/comp/remains_%d_%d_%d.raw",fPath,fEvent,i,-1);
	  FILE *outfile = fopen(filename,"w");
	  if(!outfile)
	    {
	      cerr<<"AliHLTTPCDataCompressor::WriteRemaining : Cannot open file "<<filename<<endl;
	      exit(5);
	    }
	  //UInt_t dummy;
	  AliHLTTPCSpacePointData *points = fClusters[i][patch];//->GetDataPointer(dummy);
	  
	  memset(npoints,0,nrows*sizeof(Int_t));
	  
	  for(UInt_t j=0; j<fNcl[i][patch]; j++)
	    {
	      if(points[j].fCharge == 0) continue; //has been used
	      npoints[points[j].fPadRow]++;
	    }
	  Int_t size =0;
	  Byte_t *data = 0;
	  AliHLTTPCRemainingRow *tempPt=0;
	  
	  Int_t last_row = -2;
	  Int_t localcounter=0;
	  
	  for(UInt_t j=0; j<fNcl[i][patch]; j++)
	    {
	      if(points[j].fCharge == 0) continue; //has been used
	      
	      Int_t padrow = points[j].fPadRow;
	      if(padrow != last_row)
		{
		  if(last_row != -2)
		    {
		      if(!tempPt)
			{
			  cerr<<"AliHLTTPCDataCompressor::WriteRemaining : Zero row pointer "<<endl;
			  exit(5);
			}
		      if(localcounter != tempPt->fNClusters)
			{
			  cerr<<"AliHLTTPCDataCompressor::WriteRemaining : Mismatching clustercounter "<<localcounter<<" "
			      <<(Int_t)tempPt->fNClusters<<endl;
			  exit(5);
			}
		      //cout<<"Writing row "<<(int)tempPt->fPadRow<<" with "<<(int)tempPt->fNClusters<<" clusters"<<endl;
		      fwrite(tempPt,size,1,outfile);
		    }
		  if(data)
		    delete [] data;
		  size = sizeof(AliHLTTPCRemainingRow) + npoints[padrow]*sizeof(AliHLTTPCRemainingCluster);
		  data = new Byte_t[size];
		  tempPt = (AliHLTTPCRemainingRow*)data;
		  
		  localcounter=0;
		  tempPt->fPadRow = padrow;
		  tempPt->fNClusters = npoints[padrow];
		  last_row = padrow;
		}
	      if(localcounter >= npoints[padrow])
		{
		  cerr<<"AliHLTTPCDataCompressor::WriteRemaining : Cluster counter out of range: "
		      <<localcounter<<" "<<npoints[padrow]<<endl;
		  exit(5);
		}
	      
	      Float_t xyz[3] = {points[j].fX,points[j].fY,points[j].fZ};
	      AliHLTTPCTransform::Global2Local(xyz,i,kTRUE);
	      
	      tempPt->fClusters[localcounter].fY = xyz[1];
	      tempPt->fClusters[localcounter].fZ = xyz[2];
	      tempPt->fClusters[localcounter].fCharge = points[j].fCharge;
	      tempPt->fClusters[localcounter].fSigmaY2 = points[j].fSigmaY2;
	      tempPt->fClusters[localcounter].fSigmaZ2 = points[j].fSigmaZ2;
	      localcounter++;
	      fNunusedClusters++;
	    }
	  
	  //Write the last row:
	  fwrite(tempPt,size,1,outfile);
	  if(data)
	    delete [] data;
	  fclose(outfile);
	}
    }
  delete [] npoints;
}

void AliHLTTPCDataCompressor::SelectRemainingClusters()
{
  //Select which remaining clusters to write in addition to the compressed data.
  //In particular one can here make sure that "important" clusters are not missed:
  //The offline track finder perform seed finding in the outer padrows;
  //the first seeding is using pair of points on outermost padrow and 
  //0.125*nrows more rows towards the vertex. The second seeding uses pair
  //of points on the outermost padrow-0.5*0.125*nrows and 0.125*nrows + 0.5*0.125*nrows
  //more rows towards the vertex. In order to evaluate the seeds, the track offline
  //track finder checks whether a certain amount of possible clusters (padrows) is 
  //attached to the track, and then the kalman filtering starts.
  //To ensure a minimal loss off efficiency, all clusters in this region should be
  //intact.....
  
#if 0
  cout<<"Cleaning up clusters"<<endl;
#endif
  Int_t nrows = AliHLTTPCTransform::GetNRows();
  Int_t gap=(Int_t)(0.125*nrows), shift=(Int_t)(0.5*gap);
  
  for(Int_t slice=0; slice<36; slice++)
    {
      //UInt_t dummy;
      AliHLTTPCSpacePointData *points = fClusters[slice][0];//->GetDataPointer(dummy);
      for(UInt_t i=0; i<fNcl[slice][0]; i++)
	{
	  if(points[i].fCharge == 0) continue; //Already removed
	  Int_t padrow = (Int_t)points[i].fPadRow;
	  
	  Float_t xyz[3] = {points[i].fX,points[i].fY,points[i].fZ};
	  Int_t sector,row;
	  AliHLTTPCTransform::Slice2Sector(slice,padrow,sector,row);
	  AliHLTTPCTransform::Global2Raw(xyz,sector,row);
	  
	  if(padrow >= nrows-1-gap-shift) continue;//save all the clusters in this region
	  
	  //if(padrow >= nrows-1-shift) continue;

	  //Save the clusters at the borders:
	  //if(xyz[1] < 3 || xyz[1] >= AliHLTTPCTransform::GetNPads(padrow)-4)
	  // continue;

	  //Save clusters on padrows used for offline seeding:
	  if(padrow == nrows - 1 || padrow == nrows - 1 - gap ||                 //First seeding
	     padrow == nrows - 1 - shift || padrow == nrows - 1 - gap - shift)   //Second seeding
	    continue;
	  
	  //Cluster did not meet any of the above criteria, so disregard it:
	  points[i].fCharge = 0;
	}
    }
  
}

void AliHLTTPCDataCompressor::CompressAndExpand()
{
  //Read tracks/clusters from file, compress data and uncompress it. Write compression rates to file.
#if 0
  cout<<"Compressing and expanding data"<<endl;
#endif
  AliHLTTPCCompress *comp = new AliHLTTPCCompress(-1,-1,fPath,fWriteClusterShape,fEvent);
  comp->CompressFile();
  comp->ExpandFile();
  comp->PrintCompRatio(fCompRatioFile);
  delete comp;
  
  //Write the ratio between used and unused clusters to comp file:
  ofstream &out = *fCompRatioFile;
  out<<fNusedClusters<<' '<<fNunusedClusters<<endl;
}


void AliHLTTPCDataCompressor::RestoreData(Bool_t remaining_only)
{
  //Restore the uncompressed data together with the remaining clusters,
  //and write to a final cluster file which serves as an input to the
  //final offline tracker.
  
#ifndef use_aliroot
   LOG(AliHLTTPCLog::kError,"AliHLTTPCDataCompressor::RestoreData","Version")
     <<"You have to compile with use_aliroot flag in order to use this function"<<ENDLOG;
#else

#if 0
  cout<<"Restoring data"<<endl;
#endif
  
  const Int_t maxpoints=500000;
  TempCluster **clusters = new TempCluster*[36];
  Int_t *ncl = new Int_t[36];
  for(Int_t i=0; i<36; i++)
    {
      ncl[i]=0;
      clusters[i] = new TempCluster[maxpoints];
    }
  
  if(!remaining_only)
    ReadUncompressedData(clusters,ncl,maxpoints);
  
  if(fKeepRemaining)
    ReadRemaining(clusters,ncl,maxpoints);
  
  Char_t filename[1024];
  sprintf(filename,"%s/digitfile.root",fPath);
  TFile *rootfile = TFile::Open(filename);
  rootfile->cd();
  AliTPCParam *param = (AliTPCParam*)rootfile->Get(AliHLTTPCTransform::GetParamName());

  AliTPCDigitsArray *darray = new AliTPCDigitsArray();
  darray->Setup(param);
  darray->SetClass("AliSimDigits");
  sprintf(filename,"TreeD_%s_%d",AliHLTTPCTransform::GetParamName(),fEvent);
  Bool_t ok = darray->ConnectTree(filename);
  if(!ok)
    {
      cerr<<"AliHLTTPCDataCompressor::RestoreData : Problems connecting tree"<<endl;
      return;
    }

  fOutputFile->cd();
    
  AliTPCClustersArray *carray = new AliTPCClustersArray();
  carray->Setup(param);
  carray->SetClusterType("AliTPCcluster");
  carray->MakeTree();
  
  Int_t totcounter=0;
  for(Int_t slice=0; slice<=35; slice++)
    {
      TempCluster **clPt = new TempCluster*[maxpoints];
#if 0
      cout<<"Sorting "<<ncl[slice]<<" clusters in slice "<<slice<<endl;
#endif
      for(Int_t i=0; i<ncl[slice]; i++)
	clPt[i] = &clusters[slice][i];
      
      QSort(clPt,0,ncl[slice]);
      
      //cout<<"padrow "<<clPt[i]->padrow<<" pad "<<clPt[i]->pad<<" time "<<clPt[i]->time<<endl;

      Int_t falseid=0;
      Int_t counter=0;
      for(Int_t padrow=AliHLTTPCTransform::GetFirstRow(-1); padrow<=AliHLTTPCTransform::GetLastRow(-1); padrow++)
	{
	  Int_t sec,row;
	  AliHLTTPCTransform::Slice2Sector(slice,padrow,sec,row);
	  AliTPCClustersRow *clrow=carray->CreateRow(sec,row);
	  AliSimDigits *digits = (AliSimDigits*)darray->LoadRow(sec,row);
	  digits->ExpandBuffer();
	  digits->ExpandTrackBuffer();
	  Int_t patch = AliHLTTPCTransform::GetPatch(padrow);
	  while(counter < ncl[slice] && clPt[counter]->padrow == padrow)
	    {
	      Float_t temp[3];
	      AliHLTTPCTransform::Raw2Local(temp,sec,row,clPt[counter]->pad,clPt[counter]->time);
	      
	      AliTPCcluster *c = new AliTPCcluster();
	      c->SetY(temp[1]);
	      c->SetZ(temp[2]);
	      c->SetQ(clPt[counter]->charge);
	      
	      c->SetSigmaY2(clPt[counter]->sigmaY2*pow(AliHLTTPCTransform::GetPadPitchWidth(patch),2));
	      c->SetSigmaZ2(clPt[counter]->sigmaZ2*pow(AliHLTTPCTransform::GetZWidth(),2));
	      Int_t pad = TMath::Nint(clPt[counter]->pad);
	      Int_t time = TMath::Nint(clPt[counter]->time);
	      
	      if(pad < 0)
		pad=0;
	      if(pad >= AliHLTTPCTransform::GetNPads(padrow))
		pad = AliHLTTPCTransform::GetNPads(padrow)-1;
	      if(time < 0 || time >= AliHLTTPCTransform::GetNTimeBins())
		cerr<<"row "<<padrow<<" pad "<<pad<<" time "<<time<<endl;
	      
	      for(Int_t lab=0; lab<3; lab++)
		{
		  Int_t label = digits->GetTrackIDFast(time,pad,lab);
		  if(label > 1)
		    c->SetLabel(label-2,lab);
		  else if(label==0)
		    c->SetLabel(-2,lab);
		  else
		    c->SetLabel(-1,lab);
		  if(lab==0 && c->GetLabel(0) < 0)
		    {
		      falseid++;
		      //AliHLTTPCTransform::Local2Global(temp,slice);
		      //cout<<"slice "<<slice<<" padrow "<<padrow<<" y "<<temp[1]<<" z "<<temp[2]<<" label "<<c->GetLabel(0)<<endl;
		    }
		}
	      //cout<<"row "<<padrow<<" pad "<<clPt[counter]->pad<<" time "<<clPt[counter]->time<<" sigmaY2 "<<c->GetSigmaY2()<<" sigmaZ2 "<<c->GetSigmaZ2()<<endl;
	      clrow->InsertCluster(c);
	      delete c;
	      counter++;
	      totcounter++;
	    }
	  carray->StoreRow(sec,row);
	  carray->ClearRow(sec,row);
	  darray->ClearRow(sec,row);
	}
      //cerr<<"Slice "<<slice<<" nclusters "<<counter<<" falseones "<<falseid<<endl;
      if(counter != ncl[slice])
	cerr<<"AliLDataCompressor::RestoreData : Mismatching cluster count :"<<counter<<" "<<ncl[slice]<<endl;
      delete [] clPt;
    }

#if 0
  cout<<"Writing "<<totcounter<<" clusters to rootfile "<<endl;
#endif

  sprintf(filename,"TreeC_TPC_%d",fEvent);
  carray->GetTree()->SetName(filename);
  carray->GetTree()->Write();
  delete carray;
  delete darray;
  rootfile->Close();
  
  for(Int_t i=0; i<36; i++)
    delete [] clusters[i];
  delete [] clusters;
  delete [] ncl;
#endif
}

void AliHLTTPCDataCompressor::ReadUncompressedData(TempCluster **clusters,Int_t *ncl,const Int_t maxpoints)
{

#if 0
  cout<<"Reading uncompressed tracks "<<endl;
#endif
  AliHLTTPCCompress *comp = new AliHLTTPCCompress(-1,-1,fPath,fWriteClusterShape,fEvent);
  
  if(!comp->ReadFile('u'))
    return;
  
  AliHLTTPCTrackArray *tracks = comp->GetTracks();
  
  Int_t charge;
  Float_t pad,time,sigmaY2,sigmaZ2;
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      for(Int_t padrow=0; padrow < AliHLTTPCTransform::GetNRows(-1); padrow++)
	{
	  if(!track->IsPresent(padrow)) continue;
	  track->GetPad(padrow,pad);
	  track->GetTime(padrow,time);
	  track->GetClusterCharge(padrow,charge);
	  track->GetXYWidth(padrow,sigmaY2);
	  track->GetZWidth(padrow,sigmaZ2);
	  Int_t slice = track->GetClusterModel(padrow)->fSlice;
	  /*
	    if(pad < -1 || pad > AliHLTTPCTransform::GetNPads(padrow) || time < -1 || time > AliHLTTPCTransform::GetNTimeBins())
	    {
	    cerr<<"AliHLTTPCDataCompressor::ReadUncompressData : Wrong pad "<<pad<<" or time "<<time<<" on row "<<padrow<<" track index "<<i<<endl;
	    track->Print();
	    exit(5);
	    }
	  */
	  if(ncl[slice] >= maxpoints)
	    {
	      cerr<<"AliHLTTPCDataCompressor::ReadUncompressedData : Too many clusters"<<endl;
	      exit(5);
	    }
	  clusters[slice][ncl[slice]].pad = pad;
	  clusters[slice][ncl[slice]].time = time;
	  clusters[slice][ncl[slice]].charge = charge;
	  clusters[slice][ncl[slice]].sigmaY2 = sigmaY2;
	  clusters[slice][ncl[slice]].sigmaZ2 = sigmaZ2;
	  clusters[slice][ncl[slice]].padrow = padrow;
	  //cout<<"row "<<padrow<<" pad "<<pad<<" time "<<time<<" charge "<<charge<<" sigmas "<<sigmaY2<<" "<<sigmaZ2<<endl;
	  ncl[slice]++;
	}
    }

  delete comp;
}

void AliHLTTPCDataCompressor::ReadRemaining(TempCluster **clusters,Int_t *ncl,const Int_t maxpoints)
{
  
  Char_t filename[1024];
#if 0
  cout<<"Reading remaining clusters "<<endl;
#endif
  AliHLTTPCMemHandler mem;
  
  for(Int_t slice=0; slice<=35; slice++)
    {
      for(Int_t p=0; p<1; p++)
	{
	  sprintf(filename,"%s/comp/remains_%d_%d_%d.raw",fPath,fEvent,slice,-1);
	  
	  mem.SetBinaryInput(filename);
	  AliHLTTPCRemainingRow *tempPt = (AliHLTTPCRemainingRow*)mem.Allocate();
	  
	  Int_t nrows=0;
	  FILE *infile = mem.GetFilePointer();
	  while(!feof(infile))
	    {
	      Byte_t *dPt = (Byte_t*)tempPt;
	      if(fread(tempPt,sizeof(AliHLTTPCRemainingRow),1,infile)!=1) break;
	      
	      dPt += sizeof(AliHLTTPCRemainingRow);
	      
	      Int_t size = sizeof(AliHLTTPCRemainingCluster)*tempPt->fNClusters;
	      
	      fread(dPt,size,1,infile);
	      dPt += size;
	      tempPt = (AliHLTTPCRemainingRow*)dPt;
	      nrows++;
	    }
	  
	  mem.CloseBinaryInput();
	  UInt_t dummy;
	  tempPt = (AliHLTTPCRemainingRow*)mem.GetDataPointer(dummy);
	  
	  for(Int_t i=0; i<nrows; i++)
	    {
	      AliHLTTPCRemainingCluster *points = tempPt->fClusters;
	      Int_t padrow = (Int_t)tempPt->fPadRow;
	      Int_t patch = AliHLTTPCTransform::GetPatch(padrow);
	      Int_t sector,row;
	      AliHLTTPCTransform::Slice2Sector(slice,padrow,sector,row);
	      //cout<<"Loading slice "<<slice<<" row "<<padrow<<" with "<<(Int_t)tempPt->fNClusters<<" clusters "<<endl;
	      for(Int_t j=0; j<tempPt->fNClusters; j++)
		{
		  
		  Float_t xyz[3] = {AliHLTTPCTransform::Row2X(padrow),points[j].fY,points[j].fZ};
		  
		  AliHLTTPCTransform::Local2Raw(xyz,sector,row);
		  
		  if(ncl[slice] >= maxpoints)
		    {
		      cerr<<"AliHLTTPCDataCompressor::ReadRemaining : Too many clusters"<<endl;
		      exit(5);
		    }
		  //cout<<"slice "<<slice<<" padrow "<<padrow<<" pad "<<xyz[1]<<" time "<<xyz[2]<<endl;
		  clusters[slice][ncl[slice]].pad = xyz[1];
		  clusters[slice][ncl[slice]].time = xyz[2];
		  clusters[slice][ncl[slice]].charge = points[j].fCharge;
		  clusters[slice][ncl[slice]].sigmaY2 = points[j].fSigmaY2/pow(AliHLTTPCTransform::GetPadPitchWidth(patch),2);
		  clusters[slice][ncl[slice]].sigmaZ2 = points[j].fSigmaZ2/pow(AliHLTTPCTransform::GetZWidth(),2);
		  clusters[slice][ncl[slice]].padrow = padrow;
		  ncl[slice]++;
		}
	      Byte_t *dPt = (Byte_t*)tempPt;
	      Int_t size = sizeof(AliHLTTPCRemainingRow) + tempPt->fNClusters*sizeof(AliHLTTPCRemainingCluster);
	      dPt += size;
	      tempPt = (AliHLTTPCRemainingRow*)dPt;
	    }
	  
	  mem.Free();
	}
    }
}

void AliHLTTPCDataCompressor::QSort(TempCluster **a, Int_t first, Int_t last)
{
  static TempCluster *tmp;
   static int i;           // "static" to save stack space
   int j;

   while (last - first > 1) {
      i = first;
      j = last;
      for (;;) {
	while (++i < last && Compare(a[i], a[first]) < 0)
	  ;
	while (--j > first && Compare(a[j], a[first]) > 0)
	  ;
         if (i >= j)
            break;

         tmp  = a[i];
         a[i] = a[j];
         a[j] = tmp;
      }
      if (j == first) {
         ++first;
         continue;
      }
      tmp = a[first];
      a[first] = a[j];
      a[j] = tmp;
      if (j - first < last - (j + 1)) {
         QSort(a, first, j);
         first = j + 1;   // QSort(j + 1, last);
      } else {
         QSort(a, j + 1, last);
         last = j;        // QSort(first, j);
      }
   }
}

Int_t AliHLTTPCDataCompressor::Compare(TempCluster *a,TempCluster *b)
{
  /*
  if(a->padrow < 0 || a->padrow > AliHLTTPCTransform::GetNRows(-1) ||
     b->padrow < 0 || b->padrow > AliHLTTPCTransform::GetNRows(-1))
    {
      cerr<<"AliHLTTPCCompressor::Compare : Wrong padrows "<<a->padrow<<" "<<b->padrow<<endl;
      exit(5);
    }
  else if(a->pad < 0 || a->pad > AliHLTTPCTransform::GetNPads(a->padrow) || 
	  b->pad < 0 || b->pad > AliHLTTPCTransform::GetNPads(b->padrow))
    {
      cerr<<"AliHLTTPCCompressor::Compare : Wrong pads "<<a->pad<<" "<<b->pad<<endl;
      exit(5);
    }
  else if(a->time < 0 || a->time > AliHLTTPCTransform::GetNTimeBins() || 
	  b->time < 0 || b->time > AliHLTTPCTransform::GetNTimeBins())
    {
      cerr<<"AliHLTTPCCompressor::Compare : Wrong timebins "<<a->time<<" "<<b->time<<endl;
      exit(5);
    }
  */
  if(a->padrow < b->padrow) return -1;
  if(a->padrow > b->padrow) return 1;

  if(rint(a->pad) == rint(b->pad) && rint(a->time) == rint(b->time)) return 0;
  
  if(rint(a->pad) < rint(b->pad)) return -1;
  if(rint(a->pad) == rint(b->pad) && rint(a->time) < rint(b->time)) return -1;
  
  return 1;
}

