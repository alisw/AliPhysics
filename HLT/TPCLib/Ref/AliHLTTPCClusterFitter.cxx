// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTTPCStandardIncludes.h"

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCClusterFitter.h"
#include "AliHLTTPCFitUtilities.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCModelTrack.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCHoughTrack.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCCompress.h"

#if GCCVERSION == 3
using namespace std;
#endif

//_____________________________________________________________
//
//  AliHLTTPCClusterFitter
//


ClassImp(AliHLTTPCClusterFitter)

Int_t AliHLTTPCClusterFitter::fBadFitError=0;
Int_t AliHLTTPCClusterFitter::fFitError=0;
Int_t AliHLTTPCClusterFitter::fResultError=0;
Int_t AliHLTTPCClusterFitter::fFitRangeError=0;

AliHLTTPCClusterFitter::AliHLTTPCClusterFitter()
{
  plane=0;
  fNmaxOverlaps = 3;
  fChiSqMax[0]=fChiSqMax[1]=12;
  fRowMin=-1;
  fRowMax=-1;
  fFitted=0;
  fFailed=0;
  fYInnerWidthFactor=1;
  fZInnerWidthFactor=1;
  fYOuterWidthFactor=1;
  fZOuterWidthFactor=1;
  fSeeds=0;
  fProcessTracks=0;
  fClusters=0;
  fNMaxClusters=0;
  fNClusters=0;
  fEvent=0;
}

AliHLTTPCClusterFitter::AliHLTTPCClusterFitter(Char_t *path)
{
  strcpy(fPath,path);
  plane=0;
  fNmaxOverlaps = 3;
  fChiSqMax[0]=fChiSqMax[1]=12;
  fRowMin=-1;
  fRowMax=-1;
  fFitted=0;
  fFailed=0;
  fYInnerWidthFactor=1;
  fZInnerWidthFactor=1;
  fYOuterWidthFactor=1;
  fZOuterWidthFactor=1;
  fSeeds=0;
  fProcessTracks=0;
  fNMaxClusters=100000;
  fClusters=0;
  fNClusters=0;
  fEvent=0;
}

AliHLTTPCClusterFitter::~AliHLTTPCClusterFitter()
{
  if(fSeeds)
    delete fSeeds;
  if(fClusters)
    delete [] fClusters;
}

void AliHLTTPCClusterFitter::Init(Int_t slice,Int_t patch,Int_t *rowrange,AliHLTTPCTrackArray *tracks)
{
  //Assuming tracklets found by the line transform

  fSlice=slice;
  fPatch=patch;
  
  if(rowrange[0] > AliHLTTPCTransform::GetLastRow(patch) || rowrange[1] < AliHLTTPCTransform::GetFirstRow(patch))
    cerr<<"AliHLTTPCClusterFitter::Init : Wrong rows "<<rowrange[0]<<" "<<rowrange[1]<<endl;
  fRowMin=rowrange[0];
  fRowMax=rowrange[1];

  if(fRowMin < 0)
    fRowMin = 0;
  if(fRowMax > AliHLTTPCTransform::GetLastRow(fPatch))
    fRowMax = AliHLTTPCTransform::GetLastRow(fPatch);
  
  fFitted=fFailed=0;
  
  Int_t ntimes = AliHLTTPCTransform::GetNTimeBins()+1;
  Int_t npads = AliHLTTPCTransform::GetNPads(AliHLTTPCTransform::GetLastRow(fPatch))+1;//Max num of pads.
  Int_t bounds = ntimes*npads;
  if(fRow)
    delete [] fRow;
  fRow = new Digit[bounds];
  if(fTracks)
    delete fTracks;
  
  fTracks = new AliHLTTPCTrackArray("AliHLTTPCModelTrack");
  
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliHLTTPCHoughTrack *track = (AliHLTTPCHoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      AliHLTTPCModelTrack *mtrack = (AliHLTTPCModelTrack*)fTracks->NextTrack();
      mtrack->Init(slice,patch);
      mtrack->SetTgl(track->GetTgl());
      mtrack->SetRowRange(rowrange[0],rowrange[1]);
      for(Int_t j=fRowMin; j<=fRowMax; j++)
	{
	  Float_t hit[3];
	  track->GetLineCrossingPoint(j,hit);
	  hit[0] += AliHLTTPCTransform::Row2X(track->GetFirstRow());
	  Float_t R = sqrt(hit[0]*hit[0] + hit[1]*hit[1]);
	  hit[2] = R*track->GetTgl();
	  Int_t se,ro;
	  AliHLTTPCTransform::Slice2Sector(slice,j,se,ro);
	  AliHLTTPCTransform::Local2Raw(hit,se,ro);
	  if(hit[1]<0 || hit[1]>=AliHLTTPCTransform::GetNPads(j) || hit[2]<0 || hit[2]>=AliHLTTPCTransform::GetNTimeBins())
	    {
	      mtrack->SetPadHit(j,-1);
	      mtrack->SetTimeHit(j,-1);
	      continue;
	    }
	  mtrack->SetPadHit(j,hit[1]);
	  mtrack->SetTimeHit(j,hit[2]);
	  mtrack->SetCrossingAngleLUT(j,fabs(track->GetPsiLine() - AliHLTTPCTransform::Pi()/2));
	  //if(mtrack->GetCrossingAngleLUT(j) > AliHLTTPCTransform::Deg2Rad(20))
	  //  cout<<"Angle "<<mtrack->GetCrossingAngleLUT(j)<<" psiline "<<track->GetPsiLine()*180/3.1415<<endl;
	  mtrack->CalculateClusterWidths(j);
	}
    }
  //  cout<<"Copied "<<fTracks->GetNTracks()<<" tracks "<<endl;
}

void AliHLTTPCClusterFitter::Init(Int_t slice,Int_t patch)
{
  fSlice=slice;
  fPatch=patch;
  
  fRowMin=AliHLTTPCTransform::GetFirstRow(patch);
  fRowMax=AliHLTTPCTransform::GetLastRow(patch);
  
  fFitted=fFailed=0;
  
  Int_t ntimes = AliHLTTPCTransform::GetNTimeBins()+1;
  Int_t npads = AliHLTTPCTransform::GetNPads(AliHLTTPCTransform::GetLastRow(fPatch))+1;//Max num of pads.
  Int_t bounds = ntimes*npads;
  if(fRow)
    delete [] fRow;
  fRow = new Digit[bounds];
  if(fTracks)
    delete fTracks;
  fTracks = new AliHLTTPCTrackArray("AliHLTTPCModelTrack");  

}

void AliHLTTPCClusterFitter::LoadLocalSegments()
{
  Char_t filename[1024];
  sprintf(filename,"%s/hough/tracks_ho_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);
  AliHLTTPCMemHandler mem;
  mem.SetBinaryInput(filename);
  mem.Binary2TrackArray(fTracks);
  mem.CloseBinaryInput();
  
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)fTracks->GetCheckedTrack(i);
      if(!track) continue;

      track->CalculateHelix();
      
      track->Init(fSlice,fPatch);

      for(Int_t j=fRowMin; j<=fRowMax; j++)
	{
	  //Calculate the crossing point between track and padrow
	  
	  Float_t xyz_cross[3];
	  if(!track->GetCrossingPoint(j,xyz_cross))
	    continue;
	  
	  Int_t sector,row;
	  AliHLTTPCTransform::Slice2Sector(fSlice,j,sector,row);
	  AliHLTTPCTransform::Local2Raw(xyz_cross,sector,row);
	  
	  if(xyz_cross[1] < 0 || xyz_cross[1] >= AliHLTTPCTransform::GetNPads(j) ||
	     xyz_cross[2] < 0 || xyz_cross[2] >= AliHLTTPCTransform::GetNTimeBins()) //track goes out of range
	    continue;
	  
	  track->SetPadHit(j,xyz_cross[1]);
	  track->SetTimeHit(j,xyz_cross[2]);

	  Float_t crossingangle = track->GetCrossingAngle(j);
	  track->SetCrossingAngleLUT(j,crossingangle);
	  track->CalculateClusterWidths(j);
	  track->GetClusterModel(j)->fSlice = fSlice;
	  
	}
    }
}

void AliHLTTPCClusterFitter::LoadSeeds(Int_t *rowrange,Bool_t offline,Int_t eventnr)
{
  //Function assumes _global_ tracks written to a single file.
  
#if 0
  cout<<"Loading the seeds"<<endl;
#endif
  Char_t fname[1024];
  fEvent = eventnr;
  
  if(offline)
    sprintf(fname,"%s/offline/tracks_%d.raw",fPath,fEvent);
  else
    sprintf(fname,"%s/hough/tracks_%d.raw",fPath,fEvent);
  
#if 0
  cout<<"AliHLTTPCClusterFitter::LoadSeeds : Loading input tracks from "<<fname<<endl;
#endif
  
  AliHLTTPCMemHandler tfile;
  tfile.SetBinaryInput(fname);
  
  if(fSeeds)
    delete fSeeds;
  fSeeds = new AliHLTTPCTrackArray("AliHLTTPCModelTrack");
  tfile.Binary2TrackArray(fSeeds);
  tfile.CloseBinaryInput();

  //if(!offline)
  //fSeeds->QSort();
  
  Int_t clustercount=0;
  for(Int_t i=0; i<fSeeds->GetNTracks(); i++)
    {
      AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)fSeeds->GetCheckedTrack(i);
      if(!track) continue;

      if(!offline)
	{
	  if(i==0) cerr<<"AliHLTTPCClusterFitter::LoadSeeds : Cutting on pt of 4GeV!!"<<endl;
	  if(track->GetPt() > 4.) 
	    {
	      fSeeds->Remove(i);
	      continue;
	    }
	}
      clustercount += track->GetNHits();
      track->CalculateHelix();
      
      Int_t nhits = track->GetNHits();
      UInt_t *hitids = track->GetHitNumbers();

      Int_t origslice = (hitids[nhits-1]>>25)&0x7f;//Slice of innermost point

      track->Init(origslice,-1);
      Int_t slice = origslice;
      
      //for(Int_t j=rowrange[1]; j>=rowrange[0]; j--)
      for(Int_t j=rowrange[0]; j<=rowrange[1]; j++)
	{
	  
	  //Calculate the crossing point between track and padrow
	  Float_t angle = 0; //Perpendicular to padrow in local coordinates
	  AliHLTTPCTransform::Local2GlobalAngle(&angle,slice);
	  if(!track->CalculateReferencePoint(angle,AliHLTTPCTransform::Row2X(j)))
	    {
	      //cerr<<"No crossing in slice "<<slice<<" padrow "<<j<<endl;
	      continue;
	      //track->Print();
	      //exit(5);
	    }
	  Float_t xyz_cross[3] = {track->GetPointX(),track->GetPointY(),track->GetPointZ()};
	  
	  Int_t sector,row;
	  AliHLTTPCTransform::Slice2Sector(slice,j,sector,row);
	  AliHLTTPCTransform::Global2Raw(xyz_cross,sector,row);
	  //cout<<"Examining slice "<<slice<<" row "<<j<<" pad "<<xyz_cross[1]<<" time "<<xyz_cross[2]<<endl;
	  if(xyz_cross[1] < 0 || xyz_cross[1] >= AliHLTTPCTransform::GetNPads(j)) //Track leaves the slice
	    {
	    newslice:
	      
	      Int_t tslice=slice;
	      Float_t lastcross=xyz_cross[1];
	      if(xyz_cross[1] > 0)
		{
		  if(slice == 17)
		    slice=0;
		  else if(slice == 35)
		    slice = 18;
		  else
		    slice += 1;
		}
	      else
		{
		  if(slice == 0)
		    slice = 17;
		  else if(slice==18)
		    slice = 35;
		  else
		    slice -= 1;
		}
	      if(slice < 0 || slice>35)
		{
		  cerr<<"Wrong slice "<<slice<<" on row "<<j<<endl;
		  exit(5);
		}
	      //cout<<"Track leaving, trying slice "<<slice<<endl;
	      angle=0;
	      AliHLTTPCTransform::Local2GlobalAngle(&angle,slice);
	      if(!track->CalculateReferencePoint(angle,AliHLTTPCTransform::Row2X(j)))
		{
		  cerr<<"No crossing in slice "<<slice<<" padrow "<<j<<endl;
		  continue;
		  //track->Print();
		  //exit(5);
		}
	      xyz_cross[0] = track->GetPointX();
	      xyz_cross[1] = track->GetPointY();
	      xyz_cross[2] = track->GetPointZ();
	      Int_t sector,row;
	      AliHLTTPCTransform::Slice2Sector(slice,j,sector,row);
	      AliHLTTPCTransform::Global2Raw(xyz_cross,sector,row);
	      if(xyz_cross[1] < 0 || xyz_cross[1] >= AliHLTTPCTransform::GetNPads(j)) //track is in the borderline
		{
		  if(xyz_cross[1] > 0 && lastcross > 0 || xyz_cross[1] < 0 && lastcross < 0)
		    goto newslice;
		  else
		    {
		      slice = tslice;//Track is on the border of two slices
		      continue;
		    }
		}
	    }
	  
	  if(xyz_cross[2] < 0 || xyz_cross[2] >= AliHLTTPCTransform::GetNTimeBins())//track goes out of range
	    continue;
	  
	  if(xyz_cross[1] < 0 || xyz_cross[1] >= AliHLTTPCTransform::GetNPads(j))
	    {
	      cerr<<"Slice "<<slice<<" padrow "<<j<<" pad "<<xyz_cross[1]<<" time "<<xyz_cross[2]<<endl;
	      track->Print();
	      exit(5);
	    }
	  
	  track->SetPadHit(j,xyz_cross[1]);
	  track->SetTimeHit(j,xyz_cross[2]);
	  angle=0;
	  AliHLTTPCTransform::Local2GlobalAngle(&angle,slice);
	  Float_t crossingangle = track->GetCrossingAngle(j,slice);
	  track->SetCrossingAngleLUT(j,crossingangle);
	  
	  track->CalculateClusterWidths(j);
	  
	  track->GetClusterModel(j)->fSlice = slice;
	  
	}
      memset(track->GetHitNumbers(),0,159*sizeof(UInt_t));//Reset the hitnumbers
      track->SetNHits(0);
    }
  fSeeds->Compress();
  
#if 0
  cout<<"Loaded "<<fSeeds->GetNTracks()<<" seeds and "<<clustercount<<" clusters"<<endl;
#endif
}

void AliHLTTPCClusterFitter::FindClusters()
{
  if(!fTracks)
    {
      cerr<<"AliHLTTPCClusterFitter::Process : No tracks"<<endl;
      return;
    }
  if(!fRowData)
    {
      cerr<<"AliHLTTPCClusterFitter::Process : No data "<<endl;
      return;
    }
  
  AliHLTTPCDigitRowData *rowPt = fRowData;
  AliHLTTPCDigitData *digPt=0;

  Int_t pad,time;
  Short_t charge;
  
  if(fRowMin < 0)
    {
      fRowMin = AliHLTTPCTransform::GetFirstRow(fPatch);
      fRowMax = AliHLTTPCTransform::GetLastRow(fPatch);
    }
  for(Int_t i=AliHLTTPCTransform::GetFirstRow(fPatch); i<=AliHLTTPCTransform::GetLastRow(fPatch); i++)
    {
      if((Int_t)rowPt->fRow < fRowMin)
	{
	  AliHLTTPCMemHandler::UpdateRowPointer(rowPt);
	  continue;
	}
      else if((Int_t)rowPt->fRow > fRowMax)
	break;
      else if((Int_t)rowPt->fRow != i)
	{
	  cerr<<"AliHLTTPCClusterFitter::FindClusters : Mismatching row numbering "<<i<<" "<<rowPt->fRow<<endl;
	  exit(5);
	}
      fCurrentPadRow = i;
      memset((void*)fRow,0,(AliHLTTPCTransform::GetNTimeBins()+1)*(AliHLTTPCTransform::GetNPads(i)+1)*sizeof(Digit));
      digPt = (AliHLTTPCDigitData*)rowPt->fDigitData;

      for(UInt_t j=0; j<rowPt->fNDigit; j++)
	{
	  pad = digPt[j].fPad;
	  time = digPt[j].fTime;
	  charge = digPt[j].fCharge;
	  fRow[(AliHLTTPCTransform::GetNTimeBins()+1)*pad+time].fCharge = charge;
	  fRow[(AliHLTTPCTransform::GetNTimeBins()+1)*pad+time].fUsed = kFALSE;
	  //cout<<"Row "<<i<<" pad "<<pad<<" time "<<time<<" charge "<<charge<<endl;
	}
      
      for(Int_t it=0; it<2; it++)
	{
	  if(it==0)
	    {
	      fProcessTracks = fSeeds;
	      fSeeding = kTRUE;
	    }
	  else
	    {
	      fProcessTracks = fTracks;
	      fSeeding = kFALSE;
	    }
	  if(!fProcessTracks)
	    continue;
	  
	  for(Int_t k=0; k<fProcessTracks->GetNTracks(); k++)
	    {
	      AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)fProcessTracks->GetCheckedTrack(k);
	      if(!track) continue;
	      
	      if(fSeeding)
		if(track->GetClusterModel(i)->fSlice != fSlice) continue;
	      
	      if(track->GetPadHit(i) < 0 || track->GetPadHit(i) > AliHLTTPCTransform::GetNPads(i)-1 ||
		 track->GetTimeHit(i) < 0 || track->GetTimeHit(i) > AliHLTTPCTransform::GetNTimeBins()-1)
		{
		  track->SetCluster(i,0,0,0,0,0,0);
		  continue;
		}
	      
	      if(CheckCluster(k) == kFALSE)
		fFailed++;
	    }
	}
      AliHLTTPCMemHandler::UpdateRowPointer(rowPt);
    }
  
  fSeeding = kTRUE;
  AddClusters();
  fSeeding = kFALSE;
  AddClusters();
    
#if 0
  cout<<"Fitted "<<fFitted<<" clusters, failed "<<fFailed<<endl;
  cout<<"Distribution:"<<endl;
  cout<<"Bad fit "<<fBadFitError<<endl;
  cout<<"Fit error "<<fFitError<<endl;
  cout<<"Result error "<<fResultError<<endl;
  cout<<"Fit range error "<<fFitRangeError<<endl;
#endif

}

Bool_t AliHLTTPCClusterFitter::CheckCluster(Int_t trackindex)
{
  //Check if this is a single or overlapping cluster
  
  AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)fProcessTracks->GetCheckedTrack(trackindex);
  
  Int_t row = fCurrentPadRow;
  
  if(track->IsSet(row)) //A fit has already be performed on this one
    return kTRUE;
  
  //Define the cluster region of this hit:
  Int_t padr[2]={999,-1};
  Int_t timer[2]={999,-1};
  
  if(!SetFitRange(track,padr,timer))
    {
      track->SetCluster(fCurrentPadRow,0,0,0,0,0,0);
      
      if(fDebug)
	  {
#if 0
	cout<<"Failed to fit cluster at row "<<row<<" pad "<<(Int_t)rint(track->GetPadHit(row))<<" time "
	    <<(Int_t)rint(track->GetTimeHit(row))<<" hitcharge "
	    <<fRow[(AliHLTTPCTransform::GetNTimeBins()+1)*(Int_t)rint(track->GetPadHit(row))+(Int_t)rint(track->GetTimeHit(row))].fCharge<<endl;
#endif
	  }
      fFitRangeError++;
      return kFALSE;
    }

  //Check if any other track contributes to this cluster:
  //This is done by checking if the tracks are overlapping within
  //the range defined by the track parameters
  
  for(Int_t t=trackindex+1; t<fProcessTracks->GetNTracks(); t++)
    {
      AliHLTTPCModelTrack *tr = (AliHLTTPCModelTrack*)fProcessTracks->GetCheckedTrack(t);
      if(!tr) continue;
      if(fSeeding)
	if(tr->GetClusterModel(row)->fSlice != fSlice) continue;

      Int_t xyw = (Int_t)ceil(sqrt(tr->GetParSigmaY2(row))*GetYWidthFactor()); 
      Int_t zw = (Int_t)ceil(sqrt(tr->GetParSigmaZ2(row))*GetZWidthFactor()); 
      
      if( 
	 (tr->GetPadHit(row) - xyw > padr[0] && tr->GetPadHit(row) - xyw < padr[1] &&
	  tr->GetTimeHit(row) - zw > timer[0] && tr->GetTimeHit(row) - zw < timer[1]) ||
	 
	 (tr->GetPadHit(row) + xyw > padr[0] && tr->GetPadHit(row) + xyw < padr[1] &&
	  tr->GetTimeHit(row) - zw > timer[0] && tr->GetTimeHit(row) - zw < timer[1]) ||
	 
	 (tr->GetPadHit(row) - xyw > padr[0] && tr->GetPadHit(row) - xyw < padr[1] &&
	  tr->GetTimeHit(row) + zw > timer[0] && tr->GetTimeHit(row) + zw < timer[1]) ||
	 
	 (tr->GetPadHit(row) + xyw > padr[0] && tr->GetPadHit(row) + xyw < padr[1] &&
	  tr->GetTimeHit(row) + zw > timer[0] && tr->GetTimeHit(row) + zw < timer[1]) 
	 )
	{
	  if(SetFitRange(tr,padr,timer)) //Expand the cluster fit range
	    track->SetOverlap(row,t);    //Set overlap
	}
    }

  if(fDebug)
      {
#if 0
    cout<<"Fitting cluster with "<<track->GetNOverlaps(fCurrentPadRow)<<" overlaps"<<endl;
#endif
      }
  FitClusters(track,padr,timer);
  return kTRUE;
}

Bool_t AliHLTTPCClusterFitter::SetFitRange(AliHLTTPCModelTrack *track,Int_t *padrange,Int_t *timerange)
{
  Int_t row = fCurrentPadRow;
  Int_t nt = AliHLTTPCTransform::GetNTimeBins()+1;
  
  Int_t nsearchbins=0;
  if(row < AliHLTTPCTransform::GetNRowLow())
    nsearchbins=25;
  else if(row < AliHLTTPCTransform::GetNRowLow() + AliHLTTPCTransform::GetNRowUp1())
    nsearchbins=49;
  else
    nsearchbins=49;
  
  /*
  Int_t padloop[49] = {0,0,0,-1,1,-1,1,-1,1,0,0,-1,1,-1,1 ,2,-2,2,-2,2,-2,2,-2,2,-2
		       ,0,1,2,3,3,3,3,3,3,3
		       ,2,1,0,-1,-2,-3
		       ,-3,-3,-3,-3,-3,-3,-1,-1};
  Int_t timeloop[49] = {0,1,-1,0,0,1,1,-1,-1,2,-2,2,2,-2,-2 ,0,0,1,1,-1,-1,2,2,-2,-2
                        ,-3,-3,-3,-3,-2,-1,0,1,2,3
			,3,3,3,3,3,3,2,1,0,-1,-2,-3,-3,-3};
  */
  
  Int_t padloop[49] = {0,0,0,-1,1,-1,1,-1,1,0,0,-1,1,-1,1 ,2,-2,2,-2,2,-2,2,-2,2,-2
		       ,-3,3,-3,3,-3,3,-3,3,-3,3
		       ,0,0,-1,-1,1,1,-2,-2,2,2,-3,-3,3,3};
  Int_t timeloop[49] = {0,1,-1,0,0,1,1,-1,-1,2,-2,2,2,-2,-2 ,0,0,1,1,-1,-1,2,2,-2,-2
			,0,0,-1,-1,1,1,-2,-2,2,2
			,-3,3,-3,3,-3,3,-3,3,-3,3,-3,3,-3,3};
  
  Int_t padhit = (Int_t)rint(track->GetPadHit(row));
  Int_t timehit = (Int_t)rint(track->GetTimeHit(row));
  Int_t padmax=-1;
  Int_t timemax=-1;

  for(Int_t index=0; index<nsearchbins; index++)
    {
      if(IsMaximum(padhit + padloop[index],timehit + timeloop[index])) 
	{
	  padmax = padhit + padloop[index];
	  timemax = timehit + timeloop[index];
	  break;
	}
    }


  //Define the cluster region of this hit:
  //The region we look for, is centered at the local maxima
  //and expanded around using the parametrized cluster width
  //according to track parameters.
  
  Int_t xyw = (Int_t)ceil(sqrt(track->GetParSigmaY2(row))*GetYWidthFactor());
  Int_t zw = (Int_t)ceil(sqrt(track->GetParSigmaZ2(row))*GetZWidthFactor());
  
  if(padmax>=0 && timemax>=0)
    {
      if(fDebug)
	{
#if 0
	  cout<<"Expanding cluster range using expected cluster widths: "<<xyw<<" "<<zw
	      <<" and setting local maxima pad "<<padmax<<" time "<<timemax<<endl;
	  if(xyw > 10 || zw > 10)
	    track->Print();
#endif
	}
      
      //Set the hit to the local maxima of the cluster.
      //Store the maxima in the cluster model structure,
      //-only temporary, it will be overwritten when calling SetCluster.
      
      track->GetClusterModel(row)->fDPad = padmax;
      track->GetClusterModel(row)->fDTime = timemax;

      for(Int_t i=padmax-xyw; i<=padmax+xyw; i++)
	{
	  for(Int_t j=timemax-zw; j<=timemax+zw; j++)
	    {
	      if(i<0 || i>=AliHLTTPCTransform::GetNPads(row) || j<0 || j>=AliHLTTPCTransform::GetNTimeBins()) continue;
	      if(fRow[nt*i+j].fCharge)
		{
		  if(i < padrange[0]) padrange[0]=i;
		  if(i > padrange[1]) padrange[1]=i;
		  if(j < timerange[0]) timerange[0]=j;
		  if(j > timerange[1]) timerange[1]=j;
		}
	    }
	}
      if(fDebug)
	  {
#if 0
	cout<<"New padrange "<<padrange[0]<<" "<<padrange[1]<<" "<<" time "<<timerange[0]<<" "<<timerange[1]<<endl;
#endif
	  }
      return kTRUE;
    }
  return kFALSE;
}

Bool_t AliHLTTPCClusterFitter::IsMaximum(Int_t pad,Int_t time)
{
  if(pad<0 || pad >= AliHLTTPCTransform::GetNPads(fCurrentPadRow) ||
     time<0 || time >= AliHLTTPCTransform::GetNTimeBins())
    return kFALSE;
  Int_t nt = AliHLTTPCTransform::GetNTimeBins()+1;
  if(fRow[nt*pad+time].fUsed == kTRUE) return kFALSE; //Peak has been assigned before
  Int_t charge = fRow[nt*pad+time].fCharge;
  if(charge == 1023 || charge==0) return kFALSE;
  
  //fRow[nt*pad+time].fUsed = kTRUE;
  //return kTRUE;

  if(charge < fRow[nt*(pad-1)+(time-1)].fCharge) return kFALSE;
  if(charge < fRow[nt*(pad)+(time-1)].fCharge) return kFALSE;
  if(charge < fRow[nt*(pad+1)+(time-1)].fCharge) return kFALSE;
  if(charge < fRow[nt*(pad-1)+(time)].fCharge) return kFALSE;
  if(charge < fRow[nt*(pad+1)+(time)].fCharge) return kFALSE;
  if(charge < fRow[nt*(pad-1)+(time+1)].fCharge) return kFALSE;
  if(charge < fRow[nt*(pad)+(time+1)].fCharge) return kFALSE;
  if(charge < fRow[nt*(pad+1)+(time+1)].fCharge) return kFALSE;
  fRow[nt*pad+time].fUsed = kTRUE;
  return kTRUE;
}

void AliHLTTPCClusterFitter::FitClusters(AliHLTTPCModelTrack *track,Int_t *padrange,Int_t *timerange)
{
  //Handle single and overlapping clusters
    
  //Check whether this cluster has been set before:
  
  Int_t size = FIT_PTS;
  Int_t max_tracks = FIT_MAXPAR/NUM_PARS;
  if(track->GetNOverlaps(fCurrentPadRow) > max_tracks)
    {
      cerr<<"AliHLTTPCClusterFitter::FitOverlappingClusters : Too many overlapping tracks"<<endl;
      return;
    }
  Int_t *overlaps = track->GetOverlaps(fCurrentPadRow);
  
  //Check if at least one cluster is not already fitted
  Bool_t all_fitted=kTRUE;
  
  Int_t k=-1;
  while(k < track->GetNOverlaps(fCurrentPadRow))
    {
      AliHLTTPCModelTrack *tr=0;
      if(k==-1)
 	tr = track;
      else
	tr = (AliHLTTPCModelTrack*)fProcessTracks->GetCheckedTrack(overlaps[k]);
      k++;
      if(!tr) continue;
      if(!tr->IsSet(fCurrentPadRow) && !tr->IsPresent(fCurrentPadRow))//cluster has not been set and is not present
	{
	  all_fitted = kFALSE;
	  break;
	}
    }
  if(all_fitted)
    {
      if(fDebug)
	  {
#if 0
	cout<<"But all the clusters were already fitted on row "<<fCurrentPadRow<<endl;
#endif
	  }
      return;
    }
  
  //Allocate fit parameters array; this is interface to the C code
  plane = new DPOINT[FIT_PTS];
  memset(plane,0,FIT_PTS*sizeof(DPOINT));

  Double_t x[FIT_PTS],y[FIT_PTS],s[FIT_PTS];
  
  //Fill the fit parameters:
  Double_t a[FIT_MAXPAR];
  Int_t lista[FIT_MAXPAR];
  Double_t dev[FIT_MAXPAR],chisq_f;
  
  Int_t fit_pars=0;
  
  Int_t n_overlaps=0;
  k=-1;
  
  //Fill the overlapping tracks:
  while(k < track->GetNOverlaps(fCurrentPadRow))
    {
      AliHLTTPCModelTrack *tr=0;
      if(k==-1)
	tr = track;
      else
	tr = (AliHLTTPCModelTrack*)fProcessTracks->GetCheckedTrack(overlaps[k]);
      k++;
      if(!tr) continue;
      
      if(tr->IsSet(fCurrentPadRow) && !tr->IsPresent(fCurrentPadRow)) continue;//Cluster fit failed before
      
      //Use the local maxima as the input to the fitting routine.
      //The local maxima is temporary stored in the cluster model:
      Int_t hitpad = (Int_t)rint(tr->GetClusterModel(fCurrentPadRow)->fDPad);  
      Int_t hittime = (Int_t)rint(tr->GetClusterModel(fCurrentPadRow)->fDTime);
      Int_t charge = fRow[(AliHLTTPCTransform::GetNTimeBins()+1)*hitpad + hittime].fCharge;
      
      if(fDebug)
	  {
#if 0
	cout<<"Fitting track cluster, pad "<<tr->GetPadHit(fCurrentPadRow)<<" time "
	    <<tr->GetTimeHit(fCurrentPadRow)<<" charge "<<charge<<" at local maxima in pad "<<hitpad
	    <<" time "<<hittime<<" xywidth "<<sqrt(tr->GetParSigmaY2(fCurrentPadRow))
	    <<" zwidth "<<sqrt(tr->GetParSigmaZ2(fCurrentPadRow))<<endl;
#endif
	  }
      
      if(charge==0)
	{
	  cerr<<"Charge still zero!"<<endl;
	  exit(5);
	}
            
      a[n_overlaps*NUM_PARS+2] = hitpad;
      a[n_overlaps*NUM_PARS+4] = hittime;
      
      if(!tr->IsSet(fCurrentPadRow)) //Cluster is not fitted before
	{
	  a[n_overlaps*NUM_PARS+1] = charge;
	  a[n_overlaps*NUM_PARS+3] = sqrt(tr->GetParSigmaY2(fCurrentPadRow)) * GetYWidthFactor();
	  a[n_overlaps*NUM_PARS+5] = sqrt(tr->GetParSigmaZ2(fCurrentPadRow)) * GetZWidthFactor();
	  a[n_overlaps*NUM_PARS+6] = sqrt(tr->GetParSigmaZ2(fCurrentPadRow)) * GetZWidthFactor();
	  lista[n_overlaps*NUM_PARS + 1] = 1;
	  lista[n_overlaps*NUM_PARS + 2] = 1;
	  lista[n_overlaps*NUM_PARS + 3] = 0;
	  lista[n_overlaps*NUM_PARS + 4] = 1;
	  lista[n_overlaps*NUM_PARS + 5] = 0;
	  lista[n_overlaps*NUM_PARS + 6] = 0;
	  fit_pars             += 3;
	}
      else  //Cluster was fitted before
	{
	  if(!tr->IsPresent(fCurrentPadRow))
	    {
	      cerr<<"AliHLTTPCClusterFitter::FindClusters : Cluster not present; there is a bug here"<<endl;
	      exit(5);
	    }
	  Int_t charge;
	  Float_t xywidth,zwidth,pad,time;
	  tr->GetPad(fCurrentPadRow,pad);
	  tr->GetTime(fCurrentPadRow,time);
	  tr->GetClusterCharge(fCurrentPadRow,charge);
	  xywidth = sqrt(tr->GetParSigmaY2(fCurrentPadRow));
	  zwidth = sqrt(tr->GetParSigmaZ2(fCurrentPadRow));
	  if(fDebug)
	      {
#if 0
	    cout<<"Cluster had been fitted before, pad "<<pad<<" time "<<time<<" charge "<<charge<<" width "<<xywidth<<" "<<zwidth<<endl;
#endif
	      }
	  
	  a[n_overlaps*NUM_PARS+2] = pad;
	  a[n_overlaps*NUM_PARS+4] = time;
	  a[n_overlaps*NUM_PARS+1] = charge;
	  a[n_overlaps*NUM_PARS+3] = sqrt(xywidth) * GetYWidthFactor();
	  a[n_overlaps*NUM_PARS+5] = sqrt(zwidth) * GetZWidthFactor();
	  a[n_overlaps*NUM_PARS+6] = sqrt(zwidth) * GetZWidthFactor();

	  lista[n_overlaps*NUM_PARS + 1] = 1;
	  lista[n_overlaps*NUM_PARS + 2] = 0;
	  lista[n_overlaps*NUM_PARS + 3] = 0;
	  lista[n_overlaps*NUM_PARS + 4] = 0;
	  lista[n_overlaps*NUM_PARS + 5] = 0;
	  lista[n_overlaps*NUM_PARS + 6] = 0;
	  fit_pars             += 1;
	}
      n_overlaps++;
    }
  
  if(n_overlaps==0) //No clusters here
    {
      delete [] plane;
      return;
    }

  Int_t pad_num=0;
  Int_t time_num_max=0;
  Int_t ndata=0;
  Int_t tot_charge=0;
  if(fDebug)
      {
#if 0
    cout<<"Padrange "<<padrange[0]<<" "<<padrange[1]<<" timerange "<<timerange[0]<<" "<<timerange[1]<<endl;
#endif
      }
  for(Int_t i=padrange[0]; i<=padrange[1]; i++)
    {
      Int_t max_charge = 0;
      Int_t time_num=0;
      for(Int_t j=timerange[0]; j<=timerange[1]; j++)
	{
	  Int_t charge = fRow[(AliHLTTPCTransform::GetNTimeBins()+1)*i + j].fCharge;
	  
	  if(charge <= 0) continue;

	  time_num++;
	  if(charge > max_charge)
	    {
	      max_charge = charge;
	      //time_num++;
	    }
	  if(fDebug)
	      {
#if 0
	    cout<<"Filling padrow "<<fCurrentPadRow<<" pad "<<i<<" time "<<j<<" charge "<<charge<<endl;
#endif
	      }
	  tot_charge += charge;
	  ndata++;
	  if(ndata >= size)
	    {
	      cerr<<"Too many points; row "<<fCurrentPadRow<<" padrange "<<padrange[0]<<" "<<padrange[1]<<" timerange "
		  <<timerange[0]<<" "<<timerange[1]<<endl;
	      exit(5);
	    }

	  plane[ndata].u = (Double_t)i;
	  plane[ndata].v = (Double_t)j;
	  x[ndata]=ndata;
	  y[ndata]=charge;
	  s[ndata]= 1 + sqrt((Double_t)charge);
	}
      if(max_charge) //there was charge on this pad
	pad_num++;
      if(time_num_max < time_num)
	time_num_max = time_num;
    }
  
  if(pad_num <= 1 || time_num_max <=1 || n_overlaps > fNmaxOverlaps || ndata <= fit_pars) //too few to do fit
    {
      SetClusterfitFalse(track);
      if(fDebug)
	  {
#if 0
	cout<<"Too few digits or too many overlaps: "<<pad_num<<" "<<time_num_max<<" "<<n_overlaps<<" ndata "<<ndata<<" fit_pars "<<fit_pars<<endl;
#endif
	  }
      delete [] plane;
      return;
    }

  
  Int_t npars = n_overlaps * NUM_PARS;
  if(fDebug)
      {
#if 0
    cout<<"Number of overlapping clusters "<<n_overlaps<<endl;
#endif
      }
  Int_t ret = lev_marq_fit( x, y, s, ndata, a, lista, dev, npars, &chisq_f, f2gauss5 );
  
  if(ret<0)
    {
      SetClusterfitFalse(track);
      fFailed++;
      fFitError++;
      delete [] plane;
      return;
      //exit(5);
    }

  chisq_f /= (ndata-fit_pars);
  if(fDebug)
      {
#if 0
    cout<<"Chisq "<<chisq_f<<endl;
#endif
      }
  
  Bool_t overlapping=kFALSE;
  if(track->GetNOverlaps(fCurrentPadRow) > 0)//There is a overlap
    overlapping=kTRUE;
  
  k=-1;
  n_overlaps=0;
  while(k < track->GetNOverlaps(fCurrentPadRow))
    {
      AliHLTTPCModelTrack *tr=0;
      if(k==-1)
	tr = track;
      else
	tr = (AliHLTTPCModelTrack*)fProcessTracks->GetCheckedTrack(overlaps[k]);
      k++;
      if(!tr) continue;
      if(!tr->IsPresent(fCurrentPadRow))
	{
	  if(tr->IsSet(fCurrentPadRow)) continue;//This cluster has been set before
	  
	  if(chisq_f < fChiSqMax[(Int_t)overlapping])//cluster fit is good enough
	    {
	      tot_charge = (Int_t)(a[n_overlaps*NUM_PARS+1] * a[n_overlaps*NUM_PARS+3] * a[n_overlaps*NUM_PARS+5]);
	      Float_t fpad = a[n_overlaps*NUM_PARS+2];
	      Float_t ftime = a[n_overlaps*NUM_PARS+4];
	      if(tot_charge < 0 || fpad < -1 || fpad > AliHLTTPCTransform::GetNPads(fCurrentPadRow) || 
		 ftime < -1 || ftime > AliHLTTPCTransform::GetNTimeBins())
		{
		  if(fDebug)
		      {
#if 0
		    cout<<"AliHLTTPCClusterFitter::Fatal result(s) in fit; in slice "<<fSlice<<" row "<<fCurrentPadRow
			<<"; pad "<<fpad<<" time "<<ftime<<" charge "<<tot_charge<<" xywidth "<<a[n_overlaps*NUM_PARS+3]
			<<" zwidth "<<a[n_overlaps*NUM_PARS+5]<<" peakcharge "<<a[n_overlaps*NUM_PARS+1]<<endl;
#endif
		      }
		  tr->SetCluster(fCurrentPadRow,0,0,0,0,0,0);
		  fFailed++;
		  fResultError++;
		  continue;
		}
	      
	      tr->SetCluster(fCurrentPadRow,fpad,ftime,tot_charge,0,0,pad_num);
	      if(fDebug)
		  {
#if 0
		cout<<"Setting cluster in pad "<<a[n_overlaps*NUM_PARS+2]<<" time "<<a[n_overlaps*NUM_PARS+4]<<" charge "<<tot_charge<<endl;
#endif
		  }
	      /*
	      //Set the digits to used:
	      for(Int_t i=padrange[0]; i<=padrange[1]; i++)
	      for(Int_t j=timerange[0]; j<=timerange[1]; j++)
	      fRow[(AliHLTTPCTransform::GetNTimeBins()+1)*i + j].fUsed = kTRUE;
	      */
	      fFitted++;
	    }
	  else //fit was too bad
	    {
	      if(fDebug)
		  {
#if 0
		cout<<"Cluster fit was too bad"<<endl;
#endif
		  }
	      tr->SetCluster(fCurrentPadRow,0,0,0,0,0,0);
	      fBadFitError++;
	      fFailed++;
	    }
	}
      n_overlaps++;
    }
  
  delete [] plane;
}

void AliHLTTPCClusterFitter::SetClusterfitFalse(AliHLTTPCModelTrack *track)
{
  //Cluster fit failed, so set the clusters to all the participating
  //tracks to zero.
  
  Int_t i=-1;
  Int_t *overlaps = track->GetOverlaps(fCurrentPadRow);
  while(i < track->GetNOverlaps(fCurrentPadRow))
    {
      AliHLTTPCModelTrack *tr=0;
      if(i==-1)
	tr = track;
      else
	tr = (AliHLTTPCModelTrack*)fProcessTracks->GetCheckedTrack(overlaps[i]);
      i++;
      if(!tr) continue;
      
      //Set the digit data to unused, so it can be fitted to another bastard:
      Int_t hitpad = (Int_t)rint(tr->GetClusterModel(fCurrentPadRow)->fDPad);  
      Int_t hittime = (Int_t)rint(tr->GetClusterModel(fCurrentPadRow)->fDTime);
      fRow[(AliHLTTPCTransform::GetNTimeBins()+1)*hitpad + hittime].fUsed = kFALSE;
      
      tr->SetCluster(fCurrentPadRow,0,0,0,0,0,0);
    }

}


void AliHLTTPCClusterFitter::AddClusters()
{
  if(!fClusters)
    {
      fClusters = new AliHLTTPCSpacePointData[fNMaxClusters];
      fNClusters=0;
    }
  
  if(fDebug)
      {
#if 0
    cout<<"Writing cluster in slice "<<fSlice<<" patch "<<fPatch<<endl;
#endif
      }
  
  AliHLTTPCTrackArray *tracks=0;
  if(fSeeding==kTRUE)
    tracks = fSeeds;
  else
    tracks = fTracks;
  
  if(!tracks)
    return;
  
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliHLTTPCModelTrack *tr = (AliHLTTPCModelTrack*)tracks->GetCheckedTrack(i);
      if(!tr) continue;
      
      UInt_t *hitids = tr->GetHitNumbers();
      Int_t nhits = tr->GetNHits();
      for(Int_t i=fRowMax; i>=fRowMin; i--)
	{
	  if(fSeeding)
	    if(tr->GetClusterModel(i)->fSlice != fSlice) continue;
	  if(!tr->IsPresent(i)) continue;
	  fCurrentPadRow = i;
	  Float_t pad,time,xywidth,zwidth;
	  Int_t charge;
	  tr->GetPad(i,pad);
	  tr->GetTime(i,time);
	  tr->GetClusterCharge(i,charge);

	  if(pad < -1 || pad >= AliHLTTPCTransform::GetNPads(i) || 
	     time < -1 || time >= AliHLTTPCTransform::GetNTimeBins())
	    {
	      continue;
#if 0
	      cout<<"slice "<<fSlice<<" row "<<i<<" pad "<<pad<<" time "<<time<<endl;
#endif
	      tr->Print();
	      exit(5);
	    }

	  tr->CalculateClusterWidths(i,kTRUE); //Parametrize errors
	  
	  tr->GetXYWidth(i,xywidth);
	  tr->GetZWidth(i,zwidth);
	  Float_t xyz[3];
	  Int_t sector,row;
	  AliHLTTPCTransform::Slice2Sector(fSlice,i,sector,row);
	  
	  AliHLTTPCTransform::Raw2Local(xyz,sector,row,pad,time);
	  
	  if(fNClusters >= fNMaxClusters)
	    {
	      cerr<<"AliHLTTPCClusterFitter::AddClusters : Too many clusters "<<fNClusters<<endl;
	      exit(5);
	    }
	  
	  fClusters[fNClusters].fX = xyz[0];
	  fClusters[fNClusters].fY = xyz[1];
	  fClusters[fNClusters].fZ = xyz[2];
	  fClusters[fNClusters].fCharge = charge;
	  fClusters[fNClusters].fPadRow = i;
	  Int_t pa = AliHLTTPCTransform::GetPatch(i);
	  if(xywidth==0 || zwidth==0)
	    cerr<<"AliHLTTPCClusterFitter::AddClusters : Cluster with zero width"<<endl;
	  if(xywidth>0)
	    fClusters[fNClusters].fSigmaY2 = xywidth*pow(AliHLTTPCTransform::GetPadPitchWidth(pa),2);
	  else
	    fClusters[fNClusters].fSigmaY2 = 1;
	  if(zwidth>0)
	    fClusters[fNClusters].fSigmaZ2 = zwidth*pow(AliHLTTPCTransform::GetZWidth(),2);
	  else
	    fClusters[fNClusters].fSigmaZ2 = 1;
	  Int_t pat=fPatch;
	  if(fPatch==-1)
	    pat=0;
	  fClusters[fNClusters].fID = fNClusters + ((fSlice&0x7f)<<25)+((pat&0x7)<<22);
	  
	  if(nhits >= AliHLTTPCTransform::GetNRows())
	    {
	      cerr<<"AliHLTTPCClusterFitter::AddClusters : Cluster counter of out range "<<nhits<<endl;
	      exit(5);
	    }
	  hitids[nhits++] = fClusters[fNClusters].fID;
	  
#ifdef do_mc
	  Int_t trackID[3];
	  Int_t fpad = (Int_t)rint(pad);
	  Int_t ftime = (Int_t)rint(time);
	  if(fpad < 0)
	    fpad=0;
	  if(fpad >= AliHLTTPCTransform::GetNPads(i))
	    fpad = AliHLTTPCTransform::GetNPads(i)-1;
	  if(ftime<0)
	    ftime=0;
	  if(ftime >= AliHLTTPCTransform::GetNTimeBins())
	    ftime = AliHLTTPCTransform::GetNTimeBins()-1;
	  GetTrackID(fpad,ftime,trackID);
	  fClusters[fNClusters].fTrackID[0] = trackID[0];
	  fClusters[fNClusters].fTrackID[1] = trackID[1];
	  fClusters[fNClusters].fTrackID[2] = trackID[2];
#endif  
	  //cout<<"Setting id "<<trackID[0]<<" on pad "<<pad<<" time "<<time<<" row "<<i<<endl;
	  fNClusters++;
	}
      
      //Copy back the number of assigned clusters
      tr->SetNHits(nhits);
      
    }
}

void AliHLTTPCClusterFitter::WriteTracks(Int_t min_hits)
{
  if(!fSeeds)
    return;
  
  AliHLTTPCTrackArray *fakes = new AliHLTTPCTrackArray();
  
  Int_t clustercount=0;
  for(Int_t i=0; i<fSeeds->GetNTracks(); i++)
    {
      AliHLTTPCModelTrack *tr = (AliHLTTPCModelTrack*)fSeeds->GetCheckedTrack(i);
      if(!tr) continue;
      if(tr->GetNHits() < min_hits)
	{
	  fakes->AddLast(tr);
	  fSeeds->Remove(i);
	}
      clustercount += tr->GetNHits();
    }
  
#if 0
  cout<<"Writing "<<clustercount<<" clusters"<<endl;
#endif
  fSeeds->Compress();
  AliHLTTPCMemHandler mem;
  Char_t filename[1024];
  sprintf(filename,"%s/fitter/tracks_%d.raw",fPath,fEvent);
  mem.SetBinaryOutput(filename);
  mem.TrackArray2Binary(fSeeds);
  mem.CloseBinaryOutput();
  
  //Write the fake tracks to its own file
  mem.Free();
  sprintf(filename,"%s/fitter/tracks_fakes_%d.raw",fPath,fEvent);
  mem.SetBinaryOutput(filename);
  mem.TrackArray2Binary(fakes);
  mem.CloseBinaryOutput();
  delete fakes;
}

void AliHLTTPCClusterFitter::WriteClusters(Bool_t global)
{
  AliHLTTPCMemHandler mem;
  if(fDebug)
      {
#if 0
    cout<<"Write "<<fNClusters<<" clusters to file"<<endl;
#endif
      }
  Char_t filename[1024];
  sprintf(filename,"%s/fitter/points_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);
  mem.SetBinaryOutput(filename);
  if(global == kTRUE)
    mem.Transform(fNClusters,fClusters,fSlice);
  mem.Memory2Binary(fNClusters,fClusters);
  mem.CloseBinaryOutput();
  mem.Free();
  
  delete [] fClusters;
  fClusters=0;
  fNClusters=0;
}

