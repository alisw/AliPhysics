//$Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ASV

#include <stream.h>
#include <iostream.h>
#include <math.h>

#include "AliL3Modeller.h"
#include "AliL3MemHandler.h"
#include "AliL3TrackArray.h"
#include "AliL3ModelTrack.h"
#include "AliL3DigitData.h"
#include "AliL3Transform.h"

#include "AliL3Defs.h"

ClassImp(AliL3Modeller)

AliL3Modeller::AliL3Modeller()
{
  fMemHandler=0;
  fTracks=0;
  fTransform=0;
}


AliL3Modeller::~AliL3Modeller()
{
  if(fMemHandler)
    delete fMemHandler;
  if(fTracks)
    delete fTracks;
  if(fTransform)
    delete fTransform;

}

void AliL3Modeller::Init(Int_t slice,Int_t patch,Char_t *path)
{
  fSlice = slice;
  fPatch = patch;
  fPadOverlap=4;
  fTimeOverlap=4;
  fTransform = new AliL3Transform();
  fTracks = new AliL3TrackArray("AliL3ModelTrack");

  AliL3MemHandler *file = new AliL3MemHandler();
  if(!file->SetBinaryInput("tracks.raw"))
    {
      cerr<<"AliL3Modeller::Init : Error opening trackfile"<<endl;
      return;
    }
  file->Binary2TrackArray(fTracks);
  file->CloseBinaryInput();
  delete file;
  
  fTracks->QSort();
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliL3ModelTrack *track = (AliL3ModelTrack*)fTracks->GetCheckedTrack(i);
      if(!track) continue;
      track->Init(fSlice,fPatch);
      track->Rotate(fSlice,kTRUE);
      track->CalculateHelix();
    }    
  
  CalculateCrossingPoints();
  CheckForOverlaps();

  fMemHandler = new AliL3MemHandler();
  Char_t fname[100];
  sprintf(fname,"%sdigits_%d_%d.raw",path,fSlice,fPatch);
  if(!fMemHandler->SetBinaryInput(fname))
    {
      cerr<<"AliL3Modeller::Init : Error opening file "<<fname<<endl;
      return;
    }
  UInt_t ndigits;
  AliL3DigitRowData *digits=(AliL3DigitRowData*)fMemHandler->CompBinary2Memory(ndigits);
  
  SetInputData(digits);
}

void AliL3Modeller::Process()
{
  
  if(!fTracks)
    {
      cerr<<"AliL3Modeller::Process : No tracks"<<endl;
      return;
    }
  if(!fRowData)
    {
      cerr<<"AliL3Modeller::Process : No data "<<endl;
      return;
    }
  
  AliL3DigitRowData *rowPt = fRowData;
  AliL3DigitData *digPt=0;

  Int_t ntimes = fTransform->GetNTimeBins()+1;
  Int_t npads = fTransform->GetNPads(NRows[fPatch][1])+1;//Max num of pads.
  Digit *row = new Digit[(ntimes)*(npads)];
  
  Int_t seq_charge;
  Int_t pad,time;
  Short_t charge;
  Cluster cluster;

  for(Int_t i=NRows[fPatch][0]; i<=NRows[fPatch][1]; i++)
    {
      memset((void*)row,0,ntimes*npads*sizeof(Digit));
      digPt = (AliL3DigitData*)rowPt->fDigitData;
      for(UInt_t j=0; j<rowPt->fNDigit; j++)
	{
	  pad = digPt[j].fPad;
	  time = digPt[j].fTime;
	  charge = digPt[j].fCharge;
	  row[ntimes*pad+time].fCharge = charge;
	  row[ntimes*pad+time].fUsed = kFALSE;
	}
      
      for(Int_t k=0; k<fTracks->GetNTracks(); k++)
	{
	  AliL3ModelTrack *track = (AliL3ModelTrack*)fTracks->GetCheckedTrack(k);
	  if(!track) continue;
	  if(track->GetOverlap()>=0) continue;//Track is overlapping
	  if(track->GetPadHit(i)<0 || track->GetTimeHit(i)<0)
	    {
	      track->SetCluster(0,0,0,0,0); //The track has left the patch.
	      continue;
	    }
	  
	  Int_t hitpad = (Int_t)rint(track->GetPadHit(i));
	  Int_t hittime = (Int_t)rint(track->GetTimeHit(i));
	  //cout<<"Checking track with pad "<<hitpad<<" time "<<hittime<<endl;
	  pad = hitpad;
	  time = hittime;
	  Int_t padsign=-1;
	  Int_t timesign=-1;
	  
	  memset(&cluster,0,sizeof(Cluster));
	  
	  while(1)//Process this padrow
	    {
	      seq_charge=0;
	      timesign=-1;
	      time = hittime;
	      while(1) //Process sequence on this pad:
		{
		  charge = row[ntimes*pad+time].fCharge;
		  if(charge==0 && timesign==-1)
		    {time=hittime+1; timesign=1; continue;}
		  else if(charge==0 && timesign==1)
		    break;
		  //cout<<"Doing pad "<<pad<<" time "<<time<<" charge "<<charge<<endl;
		  
		  seq_charge += charge;
		  		  		  
		  cluster.fTime += time*charge;
		  cluster.fPad += pad*charge;
		  cluster.fCharge += charge;
		  cluster.fSigmaY2 += pad*pad*charge;
		  cluster.fSigmaZ2 += time*time*charge;
		  
		  row[ntimes*pad+time].fUsed = kTRUE;
		  time += timesign;
		}
	      //cout<<"Finished on pad "<<pad<<" and time "<<time<<endl;
	      if(seq_charge)
		pad += padsign;
	      else //Nothing more on this pad, goto next pad
		{
		  if(padsign==-1) 
		    {
		      if(cluster.fCharge==0 && abs(pad-hitpad) < fPadOverlap/2)
			{
			  pad--; //In this case, we haven't found anything yet, 
			}        //so we will try to expand our search within the natural boundaries.
		      else
			{
			  pad=hitpad+1; 
			  padsign=1; 
			}
		      continue;
		    }
		  else if(padsign==1 && cluster.fCharge==0 && abs(pad-hitpad) < fPadOverlap/2)
		    {
		      pad++;
		      continue;
		    }
		  else //Nothing more in this cluster
		    {
		      Float_t fcharge = (Float_t)cluster.fCharge;
		      Float_t fpad = ((Float_t)cluster.fPad/fcharge);
		      Float_t ftime = ((Float_t)cluster.fTime/fcharge);
		      Float_t sigmaY2,sigmaZ2;
		      CalcClusterWidth(&cluster,sigmaY2,sigmaZ2);
		      //cout<<"row "<<i<<" charge "<<fcharge<<endl;
		      track->SetCluster(fpad,ftime,fcharge,sigmaY2,sigmaZ2);
		      break;
		    } 
		}
	      // pad += padsign;
	    }
	}
      FillZeros(rowPt,row);
      fMemHandler->UpdateRowPointer(rowPt);
    }
  delete [] row;
  
}

void AliL3Modeller::FillZeros(AliL3DigitRowData *rowPt,Digit *row)
{
  //Fill zero where data has been used.
  
  Int_t ntimes = fTransform->GetNTimeBins()+1;
  AliL3DigitData *digPt = (AliL3DigitData*)rowPt->fDigitData;
  for(UInt_t j=0; j<rowPt->fNDigit; j++)
    {
      Int_t pad = digPt[j].fPad;
      Int_t time = digPt[j].fTime;
      if(row[ntimes*pad+time].fUsed==kTRUE)
	digPt[j].fCharge = 0;
    }
}

void AliL3Modeller::WriteRemaining(Char_t *output)
{
  //Write remaining (nonzero) digits to file.
  
  cout<<"Writing remaining data to file "<<output<<endl;
  AliL3DigitRowData *rowPt;
  rowPt = (AliL3DigitRowData*)fRowData;
  Int_t digitcount=0;
  Int_t ndigits[(NumRows[fPatch])];
  for(Int_t i=NRows[fPatch][0]; i<NRows[fPatch][1]; i++)
    {
      AliL3DigitData *digPt = (AliL3DigitData*)rowPt->fDigitData;
      ndigits[(i-NRows[fPatch][0])]=0;
      for(UInt_t j=0; j<rowPt->fNDigit; j++)
	{
	  if(digPt[j].fCharge==0) continue;
	  digitcount++;
	  ndigits[(i-NRows[fPatch][0])]++;
	}
      //cout<<"Difference "<<(int)ndigits[(i-NRows[fPatch][0])]<<" "<<(int)rowPt->fNDigit<<endl;
      fMemHandler->UpdateRowPointer(rowPt);
    }
  
  Int_t size = digitcount*sizeof(AliL3DigitData) + NumRows[fPatch]*sizeof(AliL3DigitRowData);
  Byte_t *data = new Byte_t[size];
  memset(data,0,size);
  AliL3DigitRowData *tempPt = (AliL3DigitRowData*)data;
  rowPt = (AliL3DigitRowData*)fRowData;
  
  for(Int_t i=NRows[fPatch][0]; i<NRows[fPatch][1]; i++)
    {
      Int_t localcount=0;
      tempPt->fRow = i;
      tempPt->fNDigit = ndigits[(i-NRows[fPatch][0])];
      AliL3DigitData *digPt = (AliL3DigitData*)rowPt->fDigitData;
      for(UInt_t j=0; j<rowPt->fNDigit; j++)
	{
	  if(digPt[j].fCharge==0) continue;
	  if(localcount >= ndigits[(i-NRows[fPatch][0])])
	    {
	      cerr<<"AliL3Modeller::WriteRemaining : Digitarray out of range!!"<<endl;
	      return;
	    }
	  tempPt->fDigitData[localcount].fCharge = digPt[j].fCharge;
	  tempPt->fDigitData[localcount].fPad = digPt[j].fPad;
	  tempPt->fDigitData[localcount].fTime = digPt[j].fTime;
	  localcount++;
	}
      if(ndigits[(i-NRows[fPatch][0])] != localcount)
	{
	  cerr<<"AliL3Modeller::WriteRemaining : Mismatch in digitcount"<<endl;
	  return;
	}
      fMemHandler->UpdateRowPointer(rowPt);
      Byte_t *tmp = (Byte_t*)tempPt;
      Int_t size = sizeof(AliL3DigitRowData) + ndigits[(i-NRows[fPatch][0])]*sizeof(AliL3DigitData);
      tmp += size;
      tempPt = (AliL3DigitRowData*)tmp;
    }

  AliL3MemHandler *mem = new AliL3MemHandler();
  mem->SetBinaryOutput(output);
  mem->Memory2CompBinary((UInt_t)NumRows[fPatch],(AliL3DigitRowData*)data);
  mem->CloseBinaryOutput();
  delete mem;
}


void AliL3Modeller::CalculateCrossingPoints()
{
  cout<<"Calculating crossing points on "<<fTracks->GetNTracks()<<" tracks"<<endl;
  if(!fTracks)
    {
      cerr<<"AliL3Modeller::CalculateCrossingPoints(): No tracks"<<endl;
      return;
    }
  Float_t hit[3];
  for(Int_t i=NRows[fPatch][1]; i>=NRows[fPatch][0]; i--)
    {
      for(Int_t j=0; j<fTracks->GetNTracks(); j++)
	{
	  AliL3ModelTrack *track = (AliL3ModelTrack*)fTracks->GetCheckedTrack(j);
	  if(!track) continue;
	  if(track->GetNHits() < 100)
	    fTracks->Remove(j);
	  if(!track->GetCrossingPoint(i,hit)) 
	    {
	      cerr<<"AliL3Modeller::CalculateCrossingPoints : Track does not intersect line "<<endl;
	      continue;
	    }
	  //cout<<" x "<<track->GetPointX()<<" y "<<track->GetPointY()<<" z "<<track->GetPointZ()<<endl;
	  
	  fTransform->Local2Raw(hit,fSlice,i);
	  if(hit[1]<0 || hit[1]>fTransform->GetNPads(i) ||
	     hit[2]<0 || hit[2]>fTransform->GetNTimeBins())
	    {//Track is leaving the patch, so flag the track hits (<0)
	      track->SetPadHit(i,-1);
	      track->SetTimeHit(i,-1);
	      continue;
	    }
	    
	  track->SetPadHit(i,hit[1]);
	  track->SetTimeHit(i,hit[2]);
	  
	  //if(hit[1]<0 || hit[2]>445)
	  //if(hit[2]<0 || hit[2]>445)
	  //cout<<"pad "<<hit[1]<<" time "<<hit[2]<<" pt "<<track->GetPt()<<" psi "<<track->GetPsi()<<" tgl "<<track->GetTgl()<<" firstpoint "<<track->GetFirstPointX()<<" "<<track->GetFirstPointY()<<" "<<track->GetFirstPointZ()<<endl;
	  //cout<<"Crossing pad "<<hit[1]<<" time "<<hit[2]<<endl;
	}
    }
  fTracks->Compress();
  //cout<<"And there are "<<fTracks->GetNTracks()<<" tracks remaining"<<endl;
}

void AliL3Modeller::CheckForOverlaps()
{
  //Flag the tracks that overlap
  
  cout<<"Checking for overlaps"<<endl;
  
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliL3ModelTrack *track1 = (AliL3ModelTrack*)fTracks->GetCheckedTrack(i);
      if(!track1) continue;
      for(Int_t j=i+1; j<fTracks->GetNTracks(); j++)
	{
	  AliL3ModelTrack *track2 = (AliL3ModelTrack*)fTracks->GetCheckedTrack(j);
	  if(!track2) continue;
	  for(Int_t k=NRows[fPatch][0]; k<NRows[fPatch][1]; k++)
	    {
	      if(track1->GetPadHit(k)<0 || track1->GetTimeHit(k)<0 ||
		 track2->GetPadHit(k)<0 || track2->GetTimeHit(k)<0)
		continue;
	      if(fabs(track1->GetPadHit(k)-track2->GetPadHit(k)) < fPadOverlap &&
		 fabs(track1->GetTimeHit(k)-track2->GetTimeHit(k)) < fTimeOverlap)
		{
		  track1->SetOverlap(j);
		  track2->SetOverlap(i);
		}
	    }
	}
    }
  
}


void AliL3Modeller::CalcClusterWidth(Cluster *cl,Float_t &sigmaY2,Float_t &sigmaZ2)
{
  
  Float_t padw,timew;
  if(fPatch < 3)
    padw = fTransform->GetPadPitchWidthLow();
  else
    padw = fTransform->GetPadPitchWidthUp();
  Float_t charge = (Float_t)cl->fCharge;
  Float_t pad = (Float_t)cl->fPad/charge;
  Float_t time = (Float_t)cl->fTime/charge;
  Float_t s2 = (Float_t)cl->fSigmaY2/charge - pad*pad;
  sigmaY2 = (s2 + 1./12)*padw*padw;

  if(s2 != 0)
    {
      sigmaY2 = sigmaY2*0.108;
      if(fPatch<3)
	sigmaY2 = sigmaY2*2.07;
    }
  
  s2 = (Float_t)cl->fSigmaZ2/charge - time*time;
  timew = fTransform->GetZWidth();
  sigmaZ2 = (s2 +1./12)*timew*timew;
  if(s2 != 0)
    {
      sigmaZ2 = sigmaZ2*0.169;
      if(fPatch < 3)
	sigmaZ2 = sigmaZ2*1.77;
    }
  
}
