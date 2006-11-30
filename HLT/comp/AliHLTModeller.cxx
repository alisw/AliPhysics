// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group
//_____________________________________________________________
// AliHLTModeller
//
// Class for modeling TPC data.
// 
// This performs the cluster finding, based on track parameters.
// Basically it propagates the tracks to all padrows, and looks 
// for a corresponding cluster. For the moment only cog is calculated,
// and no deconvolution is done. 


#include "AliHLTStandardIncludes.h"

#include "AliHLTLogging.h"
#include "AliHLTModeller.h"
#include "AliHLTMemHandler.h"
#include "AliHLTTrackArray.h"
#include "AliHLTModelTrack.h"
#include "AliHLTDigitData.h"
#include "AliHLTTransform.h"
#include "AliHLTSpacePointData.h"

#ifdef use_aliroot
#include "AliHLTFileHandler.h"
#endif

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTModeller)

AliHLTModeller::AliHLTModeller()
{
  // default constructor
  fMemHandler=0;
  fTracks=0;
  fRow=0;
  fTrackThreshold=0;
  SetOverlap();
  SetTrackThreshold();
  SetSearchRange();
  SetMaxClusterRange(0,0);
  fDebug=kFALSE;
}


AliHLTModeller::~AliHLTModeller()
{
  // destructor
  if(fMemHandler)
    delete fMemHandler;
  if(fTracks)
    delete fTracks;
  if(fRow)
    delete [] fRow;
}

void AliHLTModeller::Init(Int_t slice,Int_t patch,Char_t *trackdata,Char_t *path,Bool_t houghtracks,Bool_t binary)
{
  // Initialization
  fSlice = slice;
  fPatch = patch;
  fHoughTracks=houghtracks;

  sprintf(fPath,"%s",path);
  
  fTracks = new AliHLTTrackArray("AliHLTModelTrack");
  
  Char_t fname[100];
  AliHLTMemHandler *file = new AliHLTMemHandler();
  if(!houghtracks)
    sprintf(fname,"%s/tracks_tr_%d_0.raw",trackdata,fSlice); //output tracks from the tracker (no merging)
  else 
    sprintf(fname,"%s/tracks_ho_%d.raw",trackdata,fSlice);
  //sprintf(fname,"%s/tracks_ho_%d_%d.raw",trackdata,fSlice,fPatch);
  if(!file->SetBinaryInput(fname))
    {
      cerr<<"AliHLTModeller::Init : Error opening trackfile: "<<fname<<endl;
      return;
    }
  file->Binary2TrackArray(fTracks);
  file->CloseBinaryInput();
  delete file;
  
  if(!houghtracks)
    fTracks->QSort();
  
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliHLTModelTrack *track = (AliHLTModelTrack*)fTracks->GetCheckedTrack(i);
      if(!track) continue;
      track->Init(fSlice,fPatch);

      //Only if the tracks has been merged across sector boundaries:
      //if(!houghtracks)
      //track->Rotate(fSlice,kTRUE); //!!!!!!!!!!!!!!!!!!!
      
      track->CalculateHelix();
    }    
  
  Int_t ntimes = AliHLTTransform::GetNTimeBins()+1;
  Int_t npads = AliHLTTransform::GetNPads(AliHLTTransform::GetLastRow(fPatch))+1;//Max num of pads.
  Int_t bounds = ntimes*npads;
  fRow = new Digit[bounds];
  
  
  UInt_t ndigits=0;
  AliHLTDigitRowData *digits=0;
#ifdef use_aliroot
  fMemHandler = new AliHLTFileHandler();
  fMemHandler->Init(slice,patch);
  if(binary == kFALSE)
    {
      sprintf(fname,"%s/digitfile.root",fPath);
      fMemHandler->SetAliInput(fname);
      digits = fMemHandler->AliAltroDigits2Memory(ndigits);
    }
  else
    {
      sprintf(fname,"%sdigits_%d_%d.raw",fPath,fSlice,fPatch);
      if(!fMemHandler->SetBinaryInput(fname))
	{
	  cerr<<"AliHLTModeller::Init : Error opening file "<<fname<<endl;
	  return;
	}
      digits=(AliHLTDigitRowData*)fMemHandler->CompBinary2Memory(ndigits);
    }
#else
  fMemHandler = new AliHLTMemHandler();
  fMemHandler->Init(slice,patch);
  if(binary == kFALSE)
    {
      cerr<<"AliHLTModeller::Init : Compile with AliROOT if you want rootfile as input"<<endl;
      return;
    }
  else
    {
      sprintf(fname,"%sdigits_%d_%d.raw",fPath,fSlice,fPatch);
      if(!fMemHandler->SetBinaryInput(fname))
	{
	  cerr<<"AliHLTModeller::Init : Error opening file "<<fname<<endl;
	  return;
	}
    }
  digits=(AliHLTDigitRowData*)fMemHandler->CompBinary2Memory(ndigits);
#endif
  
  SetInputData(digits);
}

void AliHLTModeller::FindClusters()
{
  // Finds clusters
  if(fDebug)
    cout<<"AliHLTModeller::FindClusters : Processing slice "<<fSlice<<" patch "<<fPatch<<endl;
  if(!fTracks)
    {
      cerr<<"AliHLTModeller::Process : No tracks"<<endl;
      return;
    }
  if(!fRowData)
    {
      cerr<<"AliHLTModeller::Process : No data "<<endl;
      return;
    }
  
  AliHLTDigitRowData *rowPt = fRowData;
  AliHLTDigitData *digPt=0;

  Int_t pad,time;
  Short_t charge;
  Cluster cluster;
  ClusterRegion region[200];
  
  for(Int_t i=AliHLTTransform::GetFirstRow(fPatch); i<=AliHLTTransform::GetLastRow(fPatch); i++)
    {
      if(i != (Int_t)rowPt->fRow)
	{
	  cerr<<"AliHLTModeller::FindClusters : Mismatching rownumbering "<<i<<" "<<rowPt->fRow<<endl;
	  return;
	}
      fCurrentPadRow = i;
      memset((void*)fRow,0,(AliHLTTransform::GetNTimeBins()+1)*(AliHLTTransform::GetNPads(i)+1)*sizeof(Digit));
      digPt = (AliHLTDigitData*)rowPt->fDigitData;
      //cout<<"Loading row "<<i<<" with "<<(Int_t)rowPt->fNDigit<<" digits"<<endl;
      for(UInt_t j=0; j<rowPt->fNDigit; j++)
	{
	  pad = digPt[j].fPad;
	  time = digPt[j].fTime;
	  charge = digPt[j].fCharge;
	  fRow[(AliHLTTransform::GetNTimeBins()+1)*pad + time].fCharge = charge;
	  fRow[(AliHLTTransform::GetNTimeBins()+1)*pad + time].fUsed = kFALSE;
	  //cout<<"Row "<<i<<" pad "<<pad<<" time "<<time<<" charge "<<charge<<endl;
	}
      
      for(Int_t k=0; k<fTracks->GetNTracks(); k++)
	{
	  AliHLTModelTrack *track = (AliHLTModelTrack*)fTracks->GetCheckedTrack(k);
	  if(!track) continue;
	  
	  if(track->GetPadHit(i)<0 || track->GetTimeHit(i)<0 || track->GetNOverlaps(i)>0)//track->GetOverlap(i)>=0)
	    {
	      //cout<<"Track "<<k<<" is empty on row "<<i<<" "<<track->GetPadHit(i)<<" "<<track->GetTimeHit(i)<<endl;
	      track->SetCluster(i,0,0,0,0,0,0); //The track has left the patch, or it is overlapping
	      continue;
	    }
	  
	  Int_t minpad,mintime,maxpad,maxtime;
	  minpad = mintime = 999;
	  maxpad = maxtime = 0;
	  
	  memset(&cluster,0,sizeof(Cluster));
	  LocateCluster(track,region,minpad,maxpad);//,mintime,maxtime);
	  if(maxpad - minpad + 1 > fMaxPads ||  // maxtime - mintime + 1 > fMaxTimebins ||
	     maxpad - minpad < 1)               //  || maxtime - mintime < 1)
	    {
	      //cout<<"Cluster not found on row "<<i<<" maxpad "<<maxpad<<" minpad "<<minpad<<" maxtime "<<maxtime<<" mintime "<<mintime
	      //  <<" padhit "<<track->GetPadHit(i)<<" timehit "<<track->GetTimeHit(i)<<endl;
		
	      track->SetCluster(i,0,0,0,0,0,0);
	      continue;
	    }
	  
	  Int_t npads=0;
	  for(pad=minpad; pad<=maxpad; pad++)
	    {
	      Int_t ntimes=0;
	      for(time=region[pad].fMintime; time<=region[pad].fMaxtime; time++)
		{
		  charge = fRow[(AliHLTTransform::GetNTimeBins()+1)*pad+time].fCharge;
		  if(!charge) continue;
		  if(fRow[(AliHLTTransform::GetNTimeBins()+1)*pad+time].fUsed == kTRUE)
		    continue;
		  ntimes++;
		  
		  //Update the cluster parameters with this timebin
		  cluster.fTime += time*charge;
		  cluster.fPad += pad*charge;
		  cluster.fCharge += charge;
		  cluster.fSigmaY2 += pad*pad*charge;
		  cluster.fSigmaZ2 += time*time*charge;
		  fRow[(AliHLTTransform::GetNTimeBins()+1)*pad+time].fUsed = kTRUE;
		}
	      if(ntimes)
		npads++;
	    }
	  FillCluster(track,&cluster,i,npads);
	}
      FillZeros(rowPt);
      fMemHandler->UpdateRowPointer(rowPt);
    }
  //cout<<"done processing"<<endl;
  

  //Debug:
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliHLTModelTrack *track = (AliHLTModelTrack*)fTracks->GetCheckedTrack(i);
      if(!track) continue;
      if(track->GetNClusters() != AliHLTTransform::GetNRows(fPatch))
	cerr<<endl<<"Mismatching hitcounts; nclusters: "<<track->GetNClusters()<<" nrows "<<AliHLTTransform::GetNRows(fPatch)<<endl<<endl;
    }
  
}


void AliHLTModeller::LocateCluster(AliHLTModelTrack *track,ClusterRegion *region,Int_t &padmin,Int_t &padmax)
{
  //Set the cluster range
  //This method searches for _all_ nonzeros timebins which are neigbours.
  //This makes it rather impractical when dealing with high occupancy,
  //because then you might have very large "cluster" areas from low
  //pt electrons/noise. 
  
  Int_t row=fCurrentPadRow,charge,prtmin=0,prtmax=999;
  Int_t hitpad = (Int_t)rint(track->GetPadHit(row));
  Int_t hittime = (Int_t)rint(track->GetTimeHit(row));
  Int_t tmin = hittime;
  Int_t tmax = tmin;
    
  Int_t clustercharge=0;
  Int_t pad=hitpad;
  Bool_t pm = kTRUE;
  Int_t npads=0,middlemax=tmax,middlemin=tmin;
  while(1)
    {
      Bool_t padpr=kFALSE;
      Int_t time = hittime;
      Bool_t tm = kTRUE;
      if(pad < 0)
	{
	  padmin = 0;
	  pad = hitpad+1;
	  pm = kFALSE;
	  prtmin = middlemin;
	  prtmax = middlemax;
	  continue;
	}
      else if(pad >= AliHLTTransform::GetNPads(row))
	{
	  padmax = AliHLTTransform::GetNPads(row)-1;
	  break;
	}
      
      tmin = 999;
      tmax = 0;
      //if(row==0)
      //cout<<"Starting to look in pad "<<pad<<" time "<<time<<endl;
      while(1)
	{
	  if(time < 0)
	    {
	      time = hittime+1;
	      tm = kFALSE;
	    }
	  else if(time >= AliHLTTransform::GetNTimeBins())
	    {
	      //timemax = AliHLTTransform::GetNTimeBins()-1;
	      break;
	    }
	  charge = fRow[(AliHLTTransform::GetNTimeBins()+1)*pad+time].fCharge;
	  //if(row==0)
	  //cout<<"charge "<<charge<<" at pad "<<pad<<" time "<<time<<endl;
	  if(charge>0)
	    {
	      clustercharge+=charge;
	      padpr = kTRUE;
	      if(time < tmin)
		tmin = time;
	      if(time > tmax)
		tmax = time;
	      if(tm)
		time--;
	      else
		time++;
	    }
	  else
	    {
	      if(tm)
		{
		  //if(abs(time - hittime) < fTimeSearch && padpr == kFALSE)//Keep looking
		  if(time > prtmin && npads!=0)
		    time--;
		  else
		    {
		      time = hittime+1;
		      tm=kFALSE;
		    }
		}
	      //else if(abs(time-hittime) < fTimeSearch && padpr == kFALSE)//Keep looking
	      else if(time < prtmax && npads != 0)
		time++;
	      else
		break;
	    }
	}
      if(npads==0)
	{
	  middlemax = tmax;
	  middlemin = tmin;
	}
      //      if(row==0)
      //cout<<"tmax "<<tmax<<" tmin "<<tmin<<" prtmin "<<prtmin<<" ptrmax "<<prtmax<<endl;
      
      if(padpr && tmax >= prtmin && tmin <= prtmax)//Sequence is overlapping with the previous
	{
	  //if(row==0)
	  //cout<<"Incrementing pad "<<endl;
	  npads++;
	  
	  region[pad].fMintime=tmin;
	  region[pad].fMaxtime=tmax;
	  
	  /*
	  if(tmin < timemin)
	    timemin=tmin;
	  if(tmax > timemax)
	    timemax=tmax;
	  */
	  if(pad < padmin)
	    padmin = pad;
	  if(pad > padmax)
	    padmax = pad;
	  if(pm)
	    pad--;
	  else
	    pad++;
	  
	  prtmin = tmin;
	  prtmax = tmax;
	}
      else
	{
	  if(pm)
	    {
	      if(abs(pad-hitpad)<fPadSearch && clustercharge == 0)
		pad--;
	      else
		{
		  //if(row==0)
		  //cout<<"Setting new pad "<<hitpad+1<<endl;
		  pad = hitpad+1;
		  pm = kFALSE;
		  prtmin = middlemin;
		  prtmax = middlemax;
		  continue;
		}
	    }
	  else 
	    {
	      if(abs(pad-hitpad)<fPadSearch && clustercharge==0)
		pad++;
	      else
		break;
	    }
	}
    }
  
}


void AliHLTModeller::FillCluster(AliHLTModelTrack *track,Cluster *cluster,Int_t row,Int_t npads)
{
  // Fill clusters
  if(cluster->fCharge==0)
    {
      track->SetCluster(row,0,0,0,0,0,0);
      return;
    }
  Float_t fcharge = (Float_t)cluster->fCharge;
  Float_t fpad = ((Float_t)cluster->fPad/fcharge);
  Float_t ftime = ((Float_t)cluster->fTime/fcharge);
  Float_t sigmaY2,sigmaZ2;
  CalcClusterWidth(cluster,sigmaY2,sigmaZ2);
  track->SetCluster(row,fpad,ftime,fcharge,sigmaY2,sigmaZ2,npads);
#ifdef do_mc
  Int_t trackID[3];
  GetTrackID((Int_t)rint(fpad),(Int_t)rint(ftime),trackID);
  track->SetClusterLabel(row,trackID);
#endif
}



void AliHLTModeller::FillZeros(AliHLTDigitRowData *rowPt,Bool_t reversesign)
{
  //Fill zero where data has been used.

  AliHLTDigitData *digPt = (AliHLTDigitData*)rowPt->fDigitData;
  for(UInt_t j=0; j<rowPt->fNDigit; j++)
    {
      Int_t pad = digPt[j].fPad;
      Int_t time = digPt[j].fTime;
      if(fRow[(AliHLTTransform::GetNTimeBins()+1)*pad+time].fUsed==kTRUE)
	{
	  if(reversesign)
	    {
	      if(digPt[j].fCharge < 1024)
		digPt[j].fCharge += 1024;
	    }
	  else
	    digPt[j].fCharge = 0;
	}
    }
}

void AliHLTModeller::WriteRemaining()
{
  //Write remaining (nonzero) digits to file.
  
  AliHLTDigitRowData *rowPt;
  rowPt = (AliHLTDigitRowData*)fRowData;
  Int_t digitcount=0;
  Int_t *ndigits=new Int_t[(AliHLTTransform::GetNRows(fPatch))];
  for(Int_t i=AliHLTTransform::GetFirstRow(fPatch); i<=AliHLTTransform::GetLastRow(fPatch); i++)
    {
      AliHLTDigitData *digPt = (AliHLTDigitData*)rowPt->fDigitData;
      ndigits[(i-AliHLTTransform::GetFirstRow(fPatch))]=0;
      for(UInt_t j=0; j<rowPt->fNDigit; j++)
	{
	  if(digPt[j].fCharge==0) continue;
	  digitcount++;
	  ndigits[(i-AliHLTTransform::GetFirstRow(fPatch))]++;
	}
      //cout<<"Difference "<<(int)ndigits[(i-AliHLTTransform::GetFirstRow(fPatch))]<<" "<<(int)rowPt->fNDigit<<endl;
      fMemHandler->UpdateRowPointer(rowPt);
    }
  
  Int_t size = digitcount*sizeof(AliHLTDigitData) + AliHLTTransform::GetNRows(fPatch)*sizeof(AliHLTDigitRowData);
  Byte_t *data = new Byte_t[size];
  memset(data,0,size);
  AliHLTDigitRowData *tempPt = (AliHLTDigitRowData*)data;
  rowPt = (AliHLTDigitRowData*)fRowData;
  
  for(Int_t i=AliHLTTransform::GetFirstRow(fPatch); i<=AliHLTTransform::GetLastRow(fPatch); i++)
    {
      Int_t localcount=0;
      tempPt->fRow = i;
      tempPt->fNDigit = ndigits[(i-AliHLTTransform::GetFirstRow(fPatch))];
      AliHLTDigitData *digPt = (AliHLTDigitData*)rowPt->fDigitData;
      for(UInt_t j=0; j<rowPt->fNDigit; j++)
	{
	  if(digPt[j].fCharge==0) continue;
	  if(localcount >= ndigits[(i-AliHLTTransform::GetFirstRow(fPatch))])
	    {
	      cerr<<"AliHLTModeller::WriteRemaining : Digitarray out of range!!"<<endl;
	      return;
	    }
	  tempPt->fDigitData[localcount].fCharge = digPt[j].fCharge;
	  tempPt->fDigitData[localcount].fPad = digPt[j].fPad;
	  tempPt->fDigitData[localcount].fTime = digPt[j].fTime;

	  localcount++;
	}
      if(ndigits[(i-AliHLTTransform::GetFirstRow(fPatch))] != localcount)
	{
	  cerr<<"AliHLTModeller::WriteRemaining : Mismatch in digitcount"<<endl;
	  return;
	}
      fMemHandler->UpdateRowPointer(rowPt);
      Byte_t *tmp = (Byte_t*)tempPt;
      Int_t size = sizeof(AliHLTDigitRowData) + ndigits[(i-AliHLTTransform::GetFirstRow(fPatch))]*sizeof(AliHLTDigitData);
      tmp += size;
      tempPt = (AliHLTDigitRowData*)tmp;
    }

  Char_t fname[100];
  AliHLTMemHandler *mem = new AliHLTMemHandler();
  sprintf(fname,"%s/comp/remains_%d_%d.raw",fPath,fSlice,fPatch);
  mem->Init(fSlice,fPatch);
  mem->SetBinaryOutput(fname);
  mem->Memory2CompBinary((UInt_t)AliHLTTransform::GetNRows(fPatch),(AliHLTDigitRowData*)data);
  mem->CloseBinaryOutput();
  delete mem;
  delete [] data;
  delete [] ndigits;
}

void AliHLTModeller::RemoveBadTracks()
{
  //Remove tracsk which should not be included in the compression scheme.

  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliHLTModelTrack *track = (AliHLTModelTrack*)fTracks->GetCheckedTrack(i);
      if(!track) continue;

      if(track->GetPt() < 0.08)
	{
	  fTracks->Remove(i);
	  continue;
	}

      if(!fHoughTracks)
	if(track->GetNHits() < fTrackThreshold)
	  fTracks->Remove(i);
    }
  fTracks->Compress();
  
}

void AliHLTModeller::CalculateCrossingPoints()
{
  // calculates crossing points
  if(fDebug)
    cout<<"Calculating crossing points on "<<fTracks->GetNTracks()<<" tracks"<<endl;
  if(!fTracks)
    {
      cerr<<"AliHLTModeller::CalculateCrossingPoints(): No tracks"<<endl;
      return;
    }
  Float_t hit[3];
  
  Int_t sector,row;
  for(Int_t i=AliHLTTransform::GetLastRow(fPatch); i>=AliHLTTransform::GetFirstRow(fPatch); i--)
    {
      for(Int_t j=0; j<fTracks->GetNTracks(); j++)
	{
	  AliHLTModelTrack *track = (AliHLTModelTrack*)fTracks->GetCheckedTrack(j);
	  if(!track) continue;

	  if(!track->GetCrossingPoint(i,hit)) 
	    {
	      //cerr<<"AliHLTModeller::CalculateCrossingPoints : Track "<<j<<" does not intersect row "<<i<<" :"<<endl<<
	      //	" pt "<<track->GetPt()<<
	      //	" tgl "<<track->GetTgl()<<" psi "<<track->GetPsi()<<" charge "<<track->GetCharge()<<endl;
	      //fTracks->Remove(j);
	      track->SetPadHit(i,-1);
	      track->SetTimeHit(i,-1);
	      continue;
	    }
	  //cout<<"X "<<hit[0]<<" Y "<<hit[1]<<" Z "<<hit[2]<<" tgl "<<track->GetTgl()<<endl;
	  
	  AliHLTTransform::Slice2Sector(fSlice,i,sector,row);
	  AliHLTTransform::Local2Raw(hit,sector,row);
	  //cout<<"Pad "<<hit[1]<<" time "<<hit[2]<<" in sector "<<sector<<" row "<<row<<endl;
	  if(hit[1]<0 || hit[1]>AliHLTTransform::GetNPads(i) ||
	     hit[2]<0 || hit[2]>AliHLTTransform::GetNTimeBins())
	    {//Track is leaving the patch, so flag the track hits (<0)
	      track->SetPadHit(i,-1);
	      track->SetTimeHit(i,-1);
	      continue;
	    }

	  track->SetPadHit(i,hit[1]);
	  track->SetTimeHit(i,hit[2]);
	  track->CalculateClusterWidths(i);

	  Double_t beta = track->GetCrossingAngle(i);
	  track->SetCrossingAngleLUT(i,beta);
	  
	  //if(hit[1]<0 || hit[2]>445)
	  //if(hit[2]<0 || hit[2]>445)
	  //cout<<"pad "<<hit[1]<<" time "<<hit[2]<<" pt "<<track->GetPt()<<" psi "<<track->GetPsi()<<" tgl "<<track->GetTgl()<<" firstpoint "<<track->GetFirstPointX()<<" "<<track->GetFirstPointY()<<" "<<track->GetFirstPointZ()<<endl;
	  //cout<<"Crossing pad "<<hit[1]<<" time "<<hit[2]<<endl;
	}
    }
  fTracks->Compress();
  if(fDebug)
    cout<<"And there are "<<fTracks->GetNTracks()<<" tracks remaining"<<endl;
}

void AliHLTModeller::CheckForOverlaps(Float_t dangle,Int_t *rowrange)
{
  //Flag the tracks that overlap
  
  if(fDebug)
    cout<<"Checking for overlaps on "<<fTracks->GetNTracks()<<endl;
  Int_t counter=0;
  
  for(Int_t k=AliHLTTransform::GetFirstRow(fPatch); k<=AliHLTTransform::GetLastRow(fPatch); k++)
    {
      if(rowrange)
	{
	  if(k < rowrange[0]) continue;
	  if(k > rowrange[1]) break;
	}
      for(Int_t i=0; i<fTracks->GetNTracks(); i++)
	{
	  AliHLTModelTrack *track1 = (AliHLTModelTrack*)fTracks->GetCheckedTrack(i);
	  if(!track1) continue;
	  if(track1->GetPadHit(k)<0 || track1->GetTimeHit(k)<0) continue;
	  
	  for(Int_t j=i+1; j<fTracks->GetNTracks(); j++)
	    {
	      AliHLTModelTrack *track2 = (AliHLTModelTrack*)fTracks->GetCheckedTrack(j);
	      if(!track2) continue;
	      if(track2->GetPadHit(k)<0 || track2->GetTimeHit(k)<0) continue;
	      
	      if(abs((Int_t)rint(track1->GetPadHit(k))-(Int_t)rint(track2->GetPadHit(k))) <= fPadOverlap &&
		 abs((Int_t)rint(track1->GetTimeHit(k))-(Int_t)rint(track2->GetTimeHit(k))) <= fTimeOverlap)
		{
		  if(dangle>0 && fabs(track1->GetCrossingAngleLUT(k) - track2->GetCrossingAngleLUT(k)) < dangle)
		    fTracks->Remove(j);
		  
		  //cout<<"row "<<k<<" "<<i<<" "<<j<<" "<<track1->GetPadHit(k)<<" "<<track2->GetPadHit(k)<<" "<<fabs(track1->GetCrossingAngleLUT(k) - track2->GetCrossingAngleLUT(k))<<endl;

		  else
		    track1->SetOverlap(k,j);
		  counter++;
		}
	    }
	}
    }
  fTracks->Compress();
  if(fDebug)
    cout<<"and there are "<<fTracks->GetNTracks()<<" track left"<<endl;
  //cout<<"found "<<counter<<" done"<<endl;
}


void AliHLTModeller::CalcClusterWidth(Cluster *cl,Float_t &sigmaY2,Float_t &sigmaZ2)
{
  // calculates cluster's width
  Float_t padw,timew;
  
  padw = AliHLTTransform::GetPadPitchWidth(fPatch);
  
  Float_t charge = (Float_t)cl->fCharge;
  Float_t pad = (Float_t)cl->fPad/charge;
  Float_t time = (Float_t)cl->fTime/charge;
  Float_t s2 = (Float_t)cl->fSigmaY2/charge - pad*pad;
  
  //Save the sigmas in pad and time:
  
  sigmaY2 = (s2);// + 1./12);//*padw*padw;
  
  /*Constants added by offline
    if(s2 != 0)
    {
    sigmaY2 = sigmaY2*0.108;
    if(fPatch<3)
    sigmaY2 = sigmaY2*2.07;
    }
  */

  s2 = (Float_t)cl->fSigmaZ2/charge - time*time;
  timew = AliHLTTransform::GetZWidth();
  sigmaZ2 = (s2);// +1./12);//*timew*timew;
  

  
  /*
    Constants added by offline
    if(s2 != 0)
    {
    sigmaZ2 = sigmaZ2*0.169;
    if(fPatch < 3)
    sigmaZ2 = sigmaZ2*1.77;
    }
  */
}

#ifdef do_mc
void AliHLTModeller::GetTrackID(Int_t pad,Int_t time,Int_t *trackID)
{
  // Gets track ID
  AliHLTDigitRowData *rowPt = (AliHLTDigitRowData*)fRowData;
  
  trackID[0]=trackID[1]=trackID[2]=-2;
  
  for(Int_t i=AliHLTTransform::GetFirstRow(fPatch); i<=AliHLTTransform::GetLastRow(fPatch); i++)
    {
      if(rowPt->fRow < (UInt_t)fCurrentPadRow)
	{
	  AliHLTMemHandler::UpdateRowPointer(rowPt);
	  continue;
	}
      AliHLTDigitData *digPt = (AliHLTDigitData*)rowPt->fDigitData;
      for(UInt_t j=0; j<rowPt->fNDigit; j++)
	{
	  Int_t cpad = digPt[j].fPad;
	  Int_t ctime = digPt[j].fTime;
	  if(cpad != pad) continue;
	  if(ctime != time) continue;
	  //if(cpad != pad && ctime != ctime) continue;
	  //cout<<"Reading row "<<fCurrentRow<<" pad "<<cpad<<" time "<<ctime<<" trackID "<<digPt[j].fTrackID[0]<<endl;
	  trackID[0] = digPt[j].fTrackID[0];
	  trackID[1] = digPt[j].fTrackID[1];
	  trackID[2] = digPt[j].fTrackID[2];
	  break;
	  //cout<<"Reading trackID "<<trackID[0]<<endl;
	}
      break;
    }
#else
  void AliHLTModeller::GetTrackID(Int_t /*pad*/,Int_t /*time*/,Int_t */*trackID*/)
{
  // Does nothing if do_mc undefined
  return;
#endif
}

