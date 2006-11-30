// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group
//_____________________________________________________________
//
//  AliHLTDataCompression
//
// Interface class; binary <-> AliROOT handling of TPC data compression classes.
//


#include "AliHLTStandardIncludes.h"

#include "AliHLTLogging.h"
#include "AliHLTRootTypes.h"
#include "AliHLTTransform.h"
#include "AliHLTMemHandler.h"
#include "AliHLTSpacePointData.h"
#include "AliHLTCompressAC.h"
#include "AliHLTTrackArray.h"
#include "AliHLTModelTrack.h"
#include "AliHLTBenchmark.h"
#include "AliHLTClusterFitter.h"

#ifdef use_aliroot
#include "AliHLTFileHandler.h"
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
#include <TDirectory.h>
#include <TSystem.h>
#include <TH2F.h>
#endif

#include "AliHLTDataCompressorHelper.h"
#include "AliHLTDataCompressor.h"
#include <math.h>

#if __GNUC__ == 3
using namespace std;
#endif


ClassImp(AliHLTDataCompressor)

AliHLTDataCompressor::AliHLTDataCompressor()
{
  // default constructor
  fBenchmark=0;
  fInputTracks=0;
  fKeepRemaining=kTRUE;
  fNoCompression=kFALSE;
  fEvent=0;
  fWriteClusterShape=kFALSE;
  fOutputFile=0;
  fCompRatioFile=0;
  fNusedClusters=0;
  fNunusedClusters=0;
  memset(fClusters,0,36*6*sizeof(AliHLTSpacePointData*));
}

AliHLTDataCompressor::AliHLTDataCompressor(Char_t *path,Bool_t keep,Bool_t writeshape)
{
  // constructor
  strcpy(fPath,path);
  fBenchmark = new AliHLTBenchmark();
  fInputTracks=0;
  fKeepRemaining=keep;
  fWriteClusterShape = writeshape;
  fEvent=0;
  fOutputFile=0;
  fNusedClusters=0;
  fNunusedClusters=0;
  fNoCompression=kFALSE;
  memset(fClusters,0,36*6*sizeof(AliHLTSpacePointData*));
#ifdef use_root
  Char_t name[1024];
  sprintf(name,"rm -f %s/comp/*",path);//Clean the directory
  gSystem->Exec(name);
#endif
  OpenOutputFile();
}

AliHLTDataCompressor::~AliHLTDataCompressor()
{
  // destructor
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

void AliHLTDataCompressor::DoBench(Char_t *fname)
{
  // does benchmarking
  fBenchmark->Analyze(fname);
}

void AliHLTDataCompressor::OpenOutputFile()
{
  // opens the output file
#ifndef use_aliroot
   LOG(AliHLTLog::kError,"AliHLTDataCompressor::OpenOutputFile","Version")
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
  AliTPCParam *param = (AliTPCParam*)f->Get(AliHLTTransform::GetParamName());
  sprintf(filename,"%s/comp/AliTPCclusters.root",fPath);
  fOutputFile = TFile::Open(filename,"RECREATE");
  param->Write(param->GetTitle());
  f->Close();
#endif
}

void AliHLTDataCompressor::CloseOutputFile()
{
  // closes the output file
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

void AliHLTDataCompressor::LoadData(Int_t event,Bool_t sp)
{
  // Loads data
  fSinglePatch=sp;
  fEvent=event;
  AliHLTMemHandler *clusterfile[36][6];
  Char_t fname[1024];
  for(Int_t s=0; s<=35; s++)
    {
      for(Int_t p=0; p<6; p++)
	{
	  if(fClusters[s][p])
	    delete fClusters[s][p];
	  fClusters[s][p] = 0;
	  clusterfile[s][p] = new AliHLTMemHandler();
	  if(fSinglePatch)
	    sprintf(fname,"%s/cf/points_%d_%d_%d.raw",fPath,fEvent,s,-1);
	  else
	    sprintf(fname,"%s/cf/points_%d_%d_%d.raw",fPath,fEvent,s,p);
	  clusterfile[s][p]->SetBinaryInput(fname);
	  
	  fClusters[s][p] = (AliHLTSpacePointData*)clusterfile[s][p]->Allocate();
	  clusterfile[s][p]->Binary2Memory(fNcl[s][p],fClusters[s][p]);
	  clusterfile[s][p]->CloseBinaryInput();
	  
	  if(fSinglePatch)
	    break;
	}
    }
  
  //cout<<endl<<"Reading from offlinecf"<<endl<<endl;
  sprintf(fname,"%s/cf/tracks_%d.raw",fPath,fEvent);
  AliHLTMemHandler *tfile = new AliHLTMemHandler();
  tfile->SetBinaryInput(fname);
  
  if(fInputTracks)
    delete fInputTracks;
  fInputTracks = new AliHLTTrackArray();
  tfile->Binary2TrackArray(fInputTracks);
  tfile->CloseBinaryInput();
  delete tfile;
}

void AliHLTDataCompressor::FillData(Int_t minHits,Bool_t expand)
{
  
  //Fill the track data into track and cluster structures, and write to file.
  //Preparation for compressing it.
  
  cout<<"Filling data; "<<fInputTracks->GetNTracks()<<" tracks"<<endl;
  AliHLTTrackArray *comptracks = new AliHLTTrackArray("AliHLTModelTrack");
  fInputTracks->QSort();
  for(Int_t i=0; i<fInputTracks->GetNTracks(); i++)
    {
      AliHLTTrack *intrack = fInputTracks->GetCheckedTrack(i);
      if(!intrack) continue;

      if(intrack->GetNHits()<minHits) break;
      if(intrack->GetPt()<0.1) continue;
      
      intrack->CalculateHelix();
      
      AliHLTModelTrack *outtrack = (AliHLTModelTrack*)comptracks->NextTrack();
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
	  AliHLTSpacePointData *points = fClusters[slice][patch];//->GetDataPointer(size);
	  Float_t xyz[3] = {points[pos].fX,points[pos].fY,points[pos].fZ};
	  Int_t padrow = points[pos].fPadRow;

	  //Calculate the crossing point between track and padrow
	  Float_t angle = 0; //Perpendicular to padrow in local coordinates
	  AliHLTTransform::Local2GlobalAngle(&angle,slice);
	  if(!intrack->CalculateReferencePoint(angle,AliHLTTransform::Row2X(padrow)))
	    {
	      cerr<<"AliHLTDataCompressor::FillData : Error in crossing point calc on slice "<<slice<<" row "<<padrow<<endl;
	      break;
	      //outtrack->Print(kFALSE);
	      //exit(5);
	    }
	  
	  Float_t xyzCross[3] = {intrack->GetPointX(),intrack->GetPointY(),intrack->GetPointZ()};

	  Int_t sector,row;
	  AliHLTTransform::Slice2Sector(slice,padrow,sector,row);
	  AliHLTTransform::Global2Raw(xyzCross,sector,row);
	  AliHLTTransform::Global2Raw(xyz,sector,row);
	  
	  outtrack->SetPadHit(padrow,xyzCross[1]);
	  outtrack->SetTimeHit(padrow,xyzCross[2]);

	  outtrack->SetCrossingAngleLUT(padrow,intrack->GetCrossingAngle(padrow,slice));
	  outtrack->CalculateClusterWidths(padrow,kTRUE);

	  if(fWriteClusterShape)
	    {
	      Int_t patch = AliHLTTransform::GetPatch(padrow);
	      Float_t sigmaY2 = points[pos].fSigmaY2 / pow(AliHLTTransform::GetPadPitchWidth(patch),2);
	      Float_t sigmaZ2 = points[pos].fSigmaZ2 / pow(AliHLTTransform::GetZWidth(),2);
	      outtrack->SetCluster(padrow,xyz[1],xyz[2],points[pos].fCharge,sigmaY2,sigmaZ2,3);
	    }
	  else
	    outtrack->SetCluster(padrow,xyz[1],xyz[2],points[pos].fCharge,0,0,3);
	  
	  //IMPORTANT: Set the slice in which cluster is, you need it in AliHLTModelTrack::FillTrack!
	  outtrack->GetClusterModel(padrow)->fSlice=slice;
	  points[pos].fCharge = 0;//Mark this cluster as used.
	  fNusedClusters++;
	}
      if(!expand)
	outtrack->SetNClusters(AliHLTTransform::GetNRows(-1));
    }
  if(expand)
    ExpandTrackData(comptracks);
  
  cout<<"Writing "<<comptracks->GetNTracks()<<" tracks to file"<<endl;
  AliHLTCompress *comp = new AliHLTCompress(-1,-1,fPath,fWriteClusterShape,fEvent);
  comp->WriteFile(comptracks);
  delete comp;
  delete comptracks;
  
}

void AliHLTDataCompressor::ExpandTrackData(AliHLTTrackArray *tracks)
{
  //Loop over tracks and try to assign unused clusters.
  //Only clusters which are closer than the max. residual are taken.
  
  cout<<"Expanding "<<tracks->GetNTracks()<<" tracks"<<endl;
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliHLTModelTrack *track = (AliHLTModelTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      if(track->GetNHits() == AliHLTTransform::GetNRows()) continue;
      
      Int_t nhits = track->GetNHits();
      //cout<<"Expanding track with "<<nhits<<" clusters"<<endl;
      
      Int_t lastSlice=-1;
      for(Int_t padrow=AliHLTTransform::GetNRows()-1; padrow>=0; padrow--)
	{
	  if(track->IsPresent(padrow))
	    {
	      lastSlice = track->GetClusterModel(padrow)->fSlice;
	      continue;
	    }
	  
	  if(lastSlice < 0) //the outer cluster is missing, so skip it - it will be written anyhow.
	    continue;
	  
	  //Check the slice of the next padrow:
	  Int_t nextPadrow = padrow-1;
	  Int_t nextSlice = -1;
	  while(nextPadrow >=0)
	    {
	      if(track->IsPresent(nextPadrow))
		{
		  nextSlice = track->GetClusterModel(nextPadrow)->fSlice;
		  break;
		}
	      nextPadrow--;
	    }
	  if(nextSlice>=0)
	    if(nextSlice != lastSlice)//The track crosses a slice boundary here
	      continue;
	  
 	  //UInt_t size;
	  AliHLTSpacePointData *points = fClusters[lastSlice][0];//->GetDataPointer(size);
	  
	  Float_t angle = 0;
	  AliHLTTransform::Local2GlobalAngle(&angle,lastSlice);
	  if(!track->CalculateReferencePoint(angle,AliHLTTransform::Row2X(padrow)))
	    continue;
	  Float_t xyzCross[3] = {track->GetPointX(),track->GetPointY(),track->GetPointZ()};
	  AliHLTTransform::Global2LocHLT(xyzCross,lastSlice);
	  Float_t mindist = 123456789;
	  AliHLTSpacePointData *closest=0;
	  for(UInt_t j=0; j<fNcl[lastSlice][0]; j++)
	    {
	      if(points[j].fCharge == 0) continue;// || points[j].fPadRow != padrow) continue;
	      if(points[j].fPadRow < padrow) continue;
	      if(points[j].fPadRow > padrow) break;
	      Float_t xyz[3] = {points[j].fX,points[j].fY,points[j].fZ};
	      AliHLTTransform::Global2LocHLT(xyz,lastSlice);
	      
	      //Check for overflow:
	      Int_t temp = (Int_t)rint((xyzCross[1]-xyz[1])/AliHLTDataCompressorHelper::GetXYResidualStep(padrow));
	      if( abs(temp) > 1<<(AliHLTDataCompressorHelper::GetNPadBits()-1))
		continue;
	      
	      temp = (Int_t)rint((xyzCross[2]-xyz[2])/AliHLTDataCompressorHelper::GetZResidualStep(padrow));
	      if( abs(temp) > 1<<(AliHLTDataCompressorHelper::GetNTimeBits()-1))
		continue;
	      
	      Float_t dist = sqrt( pow(xyzCross[1]-xyz[1],2) + pow(xyzCross[2]-xyz[2],2) );
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
	      AliHLTTransform::Slice2Sector(lastSlice,padrow,sector,row);
	      AliHLTTransform::Local2Raw(xyzCross,sector,row);
	      AliHLTTransform::Global2Raw(xyz,sector,row);
	      
	      track->SetPadHit(padrow,xyzCross[1]);
	      track->SetTimeHit(padrow,xyzCross[2]);
	      
	      if(fWriteClusterShape)
		{
		  Float_t angle = track->GetCrossingAngle(padrow,lastSlice);
		  track->SetCrossingAngleLUT(padrow,angle);
		  track->CalculateClusterWidths(padrow,kTRUE);
		  Int_t patch = AliHLTTransform::GetPatch(padrow);
		  Float_t sigmaY2 = closest->fSigmaY2 / pow(AliHLTTransform::GetPadPitchWidth(patch),2);
		  Float_t sigmaZ2 = closest->fSigmaZ2 / pow(AliHLTTransform::GetZWidth(),2);
		  track->SetCluster(padrow,xyz[1],xyz[2],closest->fCharge,sigmaY2,sigmaZ2,3);
		}
	      else
		track->SetCluster(padrow,xyz[1],xyz[2],closest->fCharge,0,0,3);
	      nhits++;
	      
	      //IMPORTANT: Set the slice in which cluster is, you need it in AliHLTModelTrack::FillTrack!
	      track->GetClusterModel(padrow)->fSlice=lastSlice;
	      closest->fCharge = 0;//Mark this cluster as used.
	    }
	}
      track->SetNClusters(AliHLTTransform::GetNRows());
      //cout<<"Track was assigned "<<nhits<<" clusters"<<endl;
    }
  
}



void AliHLTDataCompressor::DetermineMinBits()
{
  //Make a pass through the modelled data (after FillData has been done) to determine
  //how many bits is needed to encode the residuals _without_ overflows.
  
  AliHLTCompress *comp = new AliHLTCompress(-1,-1,fPath,fWriteClusterShape,fEvent);
  comp->ReadFile('m');
  AliHLTTrackArray *tracks = comp->GetTracks();
  if(tracks->GetNTracks()==0)
    {
      delete comp;
      return;
    }
  
  Int_t maxtime=0,maxpad=0,maxsigma=0,maxcharge=0;
  Int_t dpad,dtime,charge,dsigmaY,dsigmaZ,npadbits,ntimebits,nchargebits,nshapebits=0;
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliHLTModelTrack *track = (AliHLTModelTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      for(Int_t padrow=0; padrow<AliHLTTransform::GetNRows(); padrow++)
	{
	  if(!track->IsPresent(padrow)) continue;
	  dpad = AliHLTDataCompressorHelper::Abs(AliHLTDataCompressorHelper::Nint(track->GetClusterModel(padrow)->fDPad));
	  dtime = AliHLTDataCompressorHelper::Abs(AliHLTDataCompressorHelper::Nint(track->GetClusterModel(padrow)->fDTime));
	  charge = AliHLTDataCompressorHelper::Abs((Int_t)track->GetClusterModel(padrow)->fDCharge);
	  dsigmaY = AliHLTDataCompressorHelper::Abs(AliHLTDataCompressorHelper::Nint(track->GetClusterModel(padrow)->fDSigmaY));
	  dsigmaZ = AliHLTDataCompressorHelper::Abs(AliHLTDataCompressorHelper::Nint(track->GetClusterModel(padrow)->fDSigmaZ));
	  if(dpad > maxpad)
	    maxpad=dpad;
	  if(dtime > maxtime)
	    maxtime=dtime;
	  if(charge > maxcharge)
	    maxcharge=charge;
	  if(dsigmaY > maxsigma)
	    maxsigma=dsigmaY;
	  if(dsigmaZ > maxsigma)
	    maxsigma=dsigmaZ;
	}
    }
  cout<<"maxpad "<<maxpad<<" maxtime "<<maxtime<<" maxcharge "<<maxcharge<<endl;
  npadbits = (Int_t)ceil(log(Double_t(maxpad))/log(2.)) + 1; //need 1 extra bit to encode the sign
  ntimebits = (Int_t)ceil(log(Double_t(maxtime))/log(2.)) + 1;
  nchargebits = (Int_t)ceil(log(Double_t(maxcharge))/log(2.)); //Store as a absolute value
  if(fWriteClusterShape)
    nshapebits = (Int_t)ceil(log(Double_t(maxsigma))/log(2.)) + 1;
  
  nchargebits = AliHLTDataCompressorHelper::GetNChargeBits();
  cout<<"Updating bitnumbers; pad "<<npadbits<<" time "<<ntimebits<<" charge "<<nchargebits<<" shape "<<nshapebits<<endl;
  AliHLTDataCompressorHelper::SetBitNumbers(npadbits,ntimebits,nchargebits,nshapebits);
}

void AliHLTDataCompressor::WriteRemaining(Bool_t select)
{
  //Write remaining clusters (not assigned to any tracks) to file

  
  if(!fKeepRemaining)
    return;
  
  if(select)
    SelectRemainingClusters();
  
  if(!fSinglePatch)
    {
      cerr<<"AliHLTCompressor::WriteRemaining : You have to modify this function when not running singlepatch"<<endl;
      return;
    }
  if(!fNoCompression)
    {
      cout<<"Compressing remaining clusters "<<endl;
      AliHLTCompress *comp = new AliHLTCompress(-1,-1,fPath,fWriteClusterShape,fEvent);
      comp->CompressRemaining(fClusters,fNcl);
      delete comp;
      return;
    }
  else
    {
      cout<<"Writing remaining clusters"<<endl;
      Int_t nrows = AliHLTTransform::GetNRows();
      Int_t *npoints = new Int_t[nrows];
      Char_t filename[1024];
      for(Int_t i=0; i<=35; i++)
	{
	  for(Int_t patch=0; patch < 1; patch++)
	    {
	      sprintf(filename,"%s/comp/remains_%d_%d_%d.raw",fPath,fEvent,i,-1);
	      FILE *outfile = fopen(filename,"w");
	      if(!outfile)
		{
		  cerr<<"AliHLTDataCompressor::WriteRemaining : Cannot open file "<<filename<<endl;
		  exit(5);
		}

	      AliHLTSpacePointData *points = fClusters[i][patch];
	  
	      memset(npoints,0,nrows*sizeof(Int_t));
	  
	      for(UInt_t j=0; j<fNcl[i][patch]; j++)
		{
		  if(points[j].fCharge == 0) continue; //has been used
		  npoints[points[j].fPadRow]++;
		}
	      Int_t size =0;
	      Byte_t *data = 0;
	      AliHLTRemainingRow *tempPt=0;
	  
	      Int_t lastRow = -2;
	      Int_t localcounter=0;
	  
	      for(UInt_t j=0; j<fNcl[i][patch]; j++)
		{
		  if(points[j].fCharge == 0) continue; //has been used
	      
		  Int_t padrow = points[j].fPadRow;
		  if(padrow != lastRow)
		    {
		      if(lastRow != -2)
			{
			  if(!tempPt)
			    {
			      cerr<<"AliHLTDataCompressor::WriteRemaining : Zero row pointer "<<endl;
			      exit(5);
			    }
			  if(localcounter != tempPt->fNClusters)
			    {
			      cerr<<"AliHLTDataCompressor::WriteRemaining : Mismatching clustercounter "<<localcounter<<" "
				  <<(Int_t)tempPt->fNClusters<<endl;
			      exit(5);
			    }
			  //cout<<"Writing row "<<(int)tempPt->fPadRow<<" with "<<(int)tempPt->fNClusters<<" clusters"<<endl;
			  fwrite(tempPt,size,1,outfile);
			}
		      if(data)
			delete [] data;
		      size = sizeof(AliHLTRemainingRow) + npoints[padrow]*sizeof(AliHLTRemainingCluster);
		      data = new Byte_t[size];
		      tempPt = (AliHLTRemainingRow*)data;
		  
		      localcounter=0;
		      tempPt->fPadRow = padrow;
		      tempPt->fNClusters = npoints[padrow];
		      lastRow = padrow;
		    }
		  if(localcounter >= npoints[padrow])
		    {
		      cerr<<"AliHLTDataCompressor::WriteRemaining : Cluster counter out of range: "
			  <<localcounter<<" "<<npoints[padrow]<<endl;
		      exit(5);
		    }
	      
		  Float_t xyz[3] = {points[j].fX,points[j].fY,points[j].fZ};
		  Int_t sector,row;
		  AliHLTTransform::Slice2Sector(i,padrow,sector,row);
		  AliHLTTransform::Global2Raw(xyz,sector,row);
		  
		  Float_t padw = points[j].fSigmaY2 / pow(AliHLTTransform::GetPadPitchWidth(AliHLTTransform::GetPatch(padrow)),2);
		  Float_t timew = points[j].fSigmaZ2 / pow(AliHLTTransform::GetZWidth(),2);
		  tempPt->fClusters[localcounter].fPad = xyz[1];
		  tempPt->fClusters[localcounter].fTime = xyz[2];
		  tempPt->fClusters[localcounter].fCharge = points[j].fCharge;
		  tempPt->fClusters[localcounter].fSigmaY2 = padw;
		  tempPt->fClusters[localcounter].fSigmaZ2 = timew;
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
}

void AliHLTDataCompressor::SelectRemainingClusters()
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
  
  cout<<"Cleaning up clusters"<<endl;
  Int_t nrows = AliHLTTransform::GetNRows();
  Int_t gap=(Int_t)(0.125*nrows), shift=(Int_t)(0.5*gap);
  
  for(Int_t slice=0; slice<36; slice++)
    {
      AliHLTSpacePointData *points = fClusters[slice][0];
      for(UInt_t i=0; i<fNcl[slice][0]; i++)
	{
	  if(points[i].fCharge == 0) continue; //Already removed
	  Int_t padrow = (Int_t)points[i].fPadRow;
	  
	  //Check the widths (errors) of the cluster, and remove big bastards:
	  Float_t padw = sqrt(points[i].fSigmaY2) / AliHLTTransform::GetPadPitchWidth(AliHLTTransform::GetPatch(padrow));
	  Float_t timew = sqrt(points[i].fSigmaZ2) / AliHLTTransform::GetZWidth();
	  if(padw >= 2.55 || timew >= 2.55)//Because we use 1 byte to store
	    {
	      points[i].fCharge = 0;
	      continue;
	    }

	  Float_t xyz[3] = {points[i].fX,points[i].fY,points[i].fZ};
	  Int_t sector,row;
	  AliHLTTransform::Slice2Sector(slice,padrow,sector,row);
	  AliHLTTransform::Global2Raw(xyz,sector,row);
	  
	  if(padrow >= nrows-1-gap-shift) continue;//save all the clusters in this region
	  
	  //if(padrow >= nrows-1-shift) continue;

	  //Save the clusters at the borders:
	  //if(xyz[1] < 3 || xyz[1] >= AliHLTTransform::GetNPads(padrow)-4)
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

void AliHLTDataCompressor::CompressAndExpand(Bool_t arithmeticCoding)
{
  //Read tracks/clusters from file, compress data and uncompress it. Write compression rates to file.
  if(fNoCompression)
    return;
  
  cout<<"Compressing and expanding data"<<endl;
  AliHLTCompress *comp = 0;
  if(arithmeticCoding)
    comp = new AliHLTCompressAC(-1,-1,fPath,fWriteClusterShape,fEvent);
  else
    comp = new AliHLTCompress(-1,-1,fPath,fWriteClusterShape,fEvent);
  comp->CompressFile();
  comp->ExpandFile();
  comp->PrintCompRatio(fCompRatioFile);
  delete comp;
  
  ofstream &out = *fCompRatioFile;
  out<<AliHLTDataCompressorHelper::GetNPadBits()<<' '<<AliHLTDataCompressorHelper::GetNTimeBits()<<' '
     <<AliHLTDataCompressorHelper::GetNChargeBits()<<' '<<AliHLTDataCompressorHelper::GetNShapeBits()<<' '
     <<AliHLTDataCompressorHelper::GetNPadBitsRemaining()<<' '<<AliHLTDataCompressorHelper::GetNTimeBitsRemaining()<<' '
     <<AliHLTDataCompressorHelper::GetNShapeBitsRemaining()<<endl;
  /*
  //Write the ratio between used and unused clusters to comp file:
  out<<fNusedClusters<<' '<<fNunusedClusters<<endl;
  */
}


void AliHLTDataCompressor::RestoreData(Bool_t remainingOnly)
{
  //Restore the uncompressed data together with the remaining clusters,
  //and write to a final cluster file which serves as an input to the
  //final offline tracker.
  
#ifndef use_aliroot
   LOG(AliHLTLog::kError,"AliHLTDataCompressor::RestoreData","Version")
     <<"You have to compile with use_aliroot flag in order to use this function"<<ENDLOG;
#else

  cout<<"Restoring data"<<endl;
  
  const Int_t kmaxpoints=500000;
  TempCluster **clusters = new TempCluster*[36];
  Int_t *ncl = new Int_t[36];
  for(Int_t i=0; i<36; i++)
    {
      ncl[i]=0;
      clusters[i] = new TempCluster[kmaxpoints];
    }
  
  if(!remainingOnly)
    ReadUncompressedData(clusters,ncl,kmaxpoints);
    
  if(fKeepRemaining)
    ReadRemaining(clusters,ncl,kmaxpoints);
  
  Char_t filename[1024];
  sprintf(filename,"%s/digitfile.root",fPath);
  TFile *rootfile = TFile::Open(filename);
  rootfile->cd();
  AliTPCParam *param = (AliTPCParam*)rootfile->Get(AliHLTTransform::GetParamName());

  AliTPCDigitsArray *darray = new AliTPCDigitsArray();
  darray->Setup(param);
  darray->SetClass("AliSimDigits");
  sprintf(filename,"TreeD_%s_%d",AliHLTTransform::GetParamName(),fEvent);
  Bool_t ok = darray->ConnectTree(filename);
  if(!ok)
    {
      cerr<<"AliHLTDataCompressor::RestoreData : Problems connecting tree"<<endl;
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
      TempCluster **clPt = new TempCluster*[kmaxpoints];
      cout<<"Sorting "<<ncl[slice]<<" clusters in slice "<<slice<<endl;
      for(Int_t i=0; i<ncl[slice]; i++)
	clPt[i] = &clusters[slice][i];
      
      if(fNusedClusters)
	QSort(clPt,0,ncl[slice]);
      
      //cout<<"padrow "<<clPt[i]->fPadrow<<" pad "<<clPt[i]->fPad<<" time "<<clPt[i]->fTime<<endl;

      Int_t falseid=0;
      Int_t counter=0;
      for(Int_t padrow=AliHLTTransform::GetFirstRow(-1); padrow<=AliHLTTransform::GetLastRow(-1); padrow++)
	{
	  Int_t sec,row;
	  AliHLTTransform::Slice2Sector(slice,padrow,sec,row);
	  AliTPCClustersRow *clrow=carray->CreateRow(sec,row);
	  AliSimDigits *digits = (AliSimDigits*)darray->LoadRow(sec,row);
	  digits->ExpandBuffer();
	  digits->ExpandTrackBuffer();
	  Int_t patch = AliHLTTransform::GetPatch(padrow);
	  while(counter < ncl[slice] && clPt[counter]->fPadrow == padrow)
	    {
	      Float_t temp[3];
	      AliHLTTransform::Raw2Local(temp,sec,row,clPt[counter]->fPad,clPt[counter]->fTime);
	      
	      AliTPCcluster *c = new AliTPCcluster();
	      c->SetY(temp[1]);
	      c->SetZ(temp[2]);
	      c->SetQ(clPt[counter]->fCharge);
	      
	      c->SetSigmaY2(clPt[counter]->fSigmaY2*pow(AliHLTTransform::GetPadPitchWidth(patch),2));
	      c->SetSigmaZ2(clPt[counter]->fSigmaZ2*pow(AliHLTTransform::GetZWidth(),2));
	      Int_t pad = AliHLTDataCompressorHelper::Nint(clPt[counter]->fPad);
	      Int_t time = AliHLTDataCompressorHelper::Nint(clPt[counter]->fTime);
	      
	      if(pad < 0)
		pad=0;
	      if(pad >= AliHLTTransform::GetNPads(padrow))
		pad = AliHLTTransform::GetNPads(padrow)-1;
	      if(time < 0 || time >= AliHLTTransform::GetNTimeBins())
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
		      //AliHLTTransform::Local2Global(temp,slice);
		      //cout<<"slice "<<slice<<" padrow "<<padrow<<" y "<<temp[1]<<" z "<<temp[2]<<" label "<<c->GetLabel(0)<<endl;
		    }
		}
	      //cout<<"row "<<padrow<<" pad "<<clPt[counter]->fPad<<" time "<<clPt[counter]->fTime<<" sigmaY2 "<<c->GetSigmaY2()<<" sigmaZ2 "<<c->GetSigmaZ2()<<endl;
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

  cout<<"Writing "<<totcounter<<" clusters to rootfile "<<endl;

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

void AliHLTDataCompressor::ReadUncompressedData(TempCluster **clusters,Int_t *ncl,const Int_t kmaxpoints)
{
  // Reads uncompressed data  
  AliHLTCompress *comp = new AliHLTCompress(-1,-1,fPath,fWriteClusterShape,fEvent);
  if(fNoCompression)
    {
      cout<<endl<<"Reading unmodified data, no compression has been done here!!!!"<<endl<<endl;
      comp->ReadFile('m');//Read the unmodified data (no compression has been done).
    }
  else
    {
      cout<<"Reading uncompressed tracks "<<endl;
      comp->ReadFile('u');
    }
  
  AliHLTTrackArray *tracks = comp->GetTracks();
  
  //Float_t totcounter=0,pcounter=0,tcounter=0;
  Int_t charge;
  Float_t pad,time,sigmaY2,sigmaZ2;
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliHLTModelTrack *track = (AliHLTModelTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      for(Int_t padrow=0; padrow < AliHLTTransform::GetNRows(-1); padrow++)
	{
	  if(!track->IsPresent(padrow)) continue;
	  track->GetPad(padrow,pad);
	  track->GetTime(padrow,time);
	  track->GetClusterCharge(padrow,charge);
	  track->GetSigmaY2(padrow,sigmaY2);
	  track->GetSigmaZ2(padrow,sigmaZ2);
	  Int_t slice = track->GetClusterModel(padrow)->fSlice;
	  /*
	    if(pad < -1 || pad > AliHLTTransform::GetNPads(padrow) || time < -1 || time > AliHLTTransform::GetNTimeBins())
	    {
	    cerr<<"AliHLTDataCompressor::ReadUncompressData : Wrong pad "<<pad<<" or time "<<time<<" on row "<<padrow<<" track index "<<i<<endl;
	    track->Print();
	    exit(5);
	    }
	  */
	  if(ncl[slice] >= kmaxpoints)
	    {
	      cerr<<"AliHLTDataCompressor::ReadUncompressedData : Too many clusters"<<endl;
	      exit(5);
	    }
	  clusters[slice][ncl[slice]].fPad = pad;
	  clusters[slice][ncl[slice]].fTime = time;
	  clusters[slice][ncl[slice]].fCharge = charge;
	  clusters[slice][ncl[slice]].fSigmaY2 = sigmaY2;
	  clusters[slice][ncl[slice]].fSigmaZ2 = sigmaZ2;
	  clusters[slice][ncl[slice]].fPadrow = padrow;
	  //cout<<"row "<<padrow<<" pad "<<pad<<" time "<<time<<" charge "<<charge<<" sigmas "<<sigmaY2<<" "<<sigmaZ2<<endl;
	  ncl[slice]++;
	}
    }
  delete comp;
}

void AliHLTDataCompressor::ReadRemaining(TempCluster **clusters,Int_t *ncl,const Int_t kmaxpoints)
{
  // reads remaining clusters  
  cout<<"Reading remaining clusters "<<endl;
  if(!fNoCompression)
    {
      AliHLTCompress *comp = new AliHLTCompress(-1,-1,fPath,fWriteClusterShape,fEvent);
      comp->ExpandRemaining(clusters,ncl,kmaxpoints);
      delete comp;
      return;
    }
  else
    {
      AliHLTMemHandler mem;
      Char_t filename[1024];
      for(Int_t slice=0; slice<=35; slice++)
	{
	  for(Int_t p=0; p<1; p++)
	    {
	      sprintf(filename,"%s/comp/remains_%d_%d_%d.raw",fPath,fEvent,slice,-1);
	  
	      mem.SetBinaryInput(filename);
	      AliHLTRemainingRow *tempPt = (AliHLTRemainingRow*)mem.Allocate();
	  
	      Int_t nrows=0;
	      FILE *infile = mem.GetFilePointer();
	      while(!feof(infile))
		{
		  Byte_t *dPt = (Byte_t*)tempPt;
		  if(fread(tempPt,sizeof(AliHLTRemainingRow),1,infile)!=1) break;
	      
		  dPt += sizeof(AliHLTRemainingRow);
	      
		  Int_t size = sizeof(AliHLTRemainingCluster)*tempPt->fNClusters;
	      
		  fread(dPt,size,1,infile);
		  dPt += size;
		  tempPt = (AliHLTRemainingRow*)dPt;
		  nrows++;
		}
	  
	      mem.CloseBinaryInput();
	      UInt_t dummy;
	      tempPt = (AliHLTRemainingRow*)mem.GetDataPointer(dummy);
	  
	      for(Int_t i=0; i<nrows; i++)
		{
		  AliHLTRemainingCluster *points = tempPt->fClusters;
		  Int_t padrow = (Int_t)tempPt->fPadRow;
		  //Int_t sector,row;
		  //AliHLTTransform::Slice2Sector(slice,padrow,sector,row);
		  //cout<<"Loading slice "<<slice<<" row "<<padrow<<" with "<<(Int_t)tempPt->fNClusters<<" clusters "<<endl;
		  for(Int_t j=0; j<tempPt->fNClusters; j++)
		    {
		  
		      //Float_t xyz[3] = {AliHLTTransform::Row2X(padrow),points[j].fY,points[j].fZ};
		      //AliHLTTransform::Local2Raw(xyz,sector,row);
		  
		      if(ncl[slice] >= kmaxpoints)
			{
			  cerr<<"AliHLTDataCompressor::ReadRemaining : Too many clusters"<<endl;
			  exit(5);
			}
		      //cout<<"slice "<<slice<<" padrow "<<padrow<<" pad "<<xyz[1]<<" time "<<xyz[2]<<endl;
		      clusters[slice][ncl[slice]].fPad = points[j].fPad;
		      clusters[slice][ncl[slice]].fTime = points[j].fTime;
		      clusters[slice][ncl[slice]].fCharge = points[j].fCharge;
		      clusters[slice][ncl[slice]].fSigmaY2 = points[j].fSigmaY2;
		      clusters[slice][ncl[slice]].fSigmaZ2 = points[j].fSigmaZ2;
		      clusters[slice][ncl[slice]].fPadrow = padrow;
		      //cout<<"padrow "<<padrow<<" pad "<<clusters[slice][ncl[slice]].fPad<<" time "<<clusters[slice][ncl[slice]].fTime<<" charge "<<clusters[slice][ncl[slice]].fCharge<<" widths "<<clusters[slice][ncl[slice]].fSigmaY2<<" "<<clusters[slice][ncl[slice]].fSigmaZ2<<endl;
		      ncl[slice]++;
		    }
		  Byte_t *dPt = (Byte_t*)tempPt;
		  Int_t size = sizeof(AliHLTRemainingRow) + tempPt->fNClusters*sizeof(AliHLTRemainingCluster);
		  dPt += size;
		  tempPt = (AliHLTRemainingRow*)dPt;
		}
	  
	      mem.Free();
	    }
	}
    }
}

void AliHLTDataCompressor::QSort(TempCluster **a, Int_t first, Int_t last)
{
  // Implementation of quick sort
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

Int_t AliHLTDataCompressor::Compare(TempCluster *a,TempCluster *b)
{
  // compares two clusters
  if(a->fPadrow < b->fPadrow) return -1;
  if(a->fPadrow > b->fPadrow) return 1;

  if(rint(a->fPad) == rint(b->fPad) && rint(a->fTime) == rint(b->fTime)) return 0;
  
  if(rint(a->fPad) < rint(b->fPad)) return -1;
  if(rint(a->fPad) == rint(b->fPad) && rint(a->fTime) < rint(b->fTime)) return -1;
  
  return 1;
}

