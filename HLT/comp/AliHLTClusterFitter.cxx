// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group
/** /class AliHLTClusterFitter
//<pre>
//_____________________________________________________________
//
//  AliHLTClusterFitter
//
</pre>
*/


#include "AliHLTStandardIncludes.h"

#include "AliHLTClusterFitter.h"
#include "AliHLTFitUtilities.h"
#include "AliHLTDigitData.h"
#include "AliHLTModelTrack.h"
#include "AliHLTTrackArray.h"
#include "AliHLTMemHandler.h"
#include "AliHLTHoughTrack.h"
#include "AliHLTSpacePointData.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTClusterFitter)

Int_t AliHLTClusterFitter::fgBadFitError=0;
Int_t AliHLTClusterFitter::fgFitError=0;
Int_t AliHLTClusterFitter::fgResultError=0;
Int_t AliHLTClusterFitter::fgFitRangeError=0;

AliHLTClusterFitter::AliHLTClusterFitter()
{
  // default constructor
  plane=0;
  fNmaxOverlaps = 3;
  fChiSqMax[0]=fChiSqMax[1]=fChiSqMax[2]=12;
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

AliHLTClusterFitter::AliHLTClusterFitter(Char_t *path)
{
  // constructor
  strcpy(fPath,path);
  plane=0;
  fNmaxOverlaps = 3;
  fChiSqMax[0]=fChiSqMax[1]=fChiSqMax[2]=12;
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

AliHLTClusterFitter::~AliHLTClusterFitter()
{
  // destructor
  if(fSeeds)
    delete fSeeds;
  if(fClusters)
    delete [] fClusters;
}

void AliHLTClusterFitter::Init(Int_t slice,Int_t patch,Int_t *rowrange,AliHLTTrackArray *tracks)
{
  //Assuming tracklets found by the line transform

  fSlice=slice;
  fPatch=patch;
  
  if(rowrange[0] > AliHLTTransform::GetLastRow(patch) || rowrange[1] < AliHLTTransform::GetFirstRow(patch))
    cerr<<"AliHLTClusterFitter::Init : Wrong rows "<<rowrange[0]<<" "<<rowrange[1]<<endl;
  fRowMin=rowrange[0];
  fRowMax=rowrange[1];

  if(fRowMin < 0)
    fRowMin = 0;
  if(fRowMax > AliHLTTransform::GetLastRow(fPatch))
    fRowMax = AliHLTTransform::GetLastRow(fPatch);
  
  fFitted=fFailed=0;
  
  Int_t ntimes = AliHLTTransform::GetNTimeBins()+1;
  Int_t npads = AliHLTTransform::GetNPads(AliHLTTransform::GetLastRow(fPatch))+1;//Max num of pads.
  Int_t bounds = ntimes*npads;
  if(fRow)
    delete [] fRow;
  fRow = new Digit[bounds];
  if(fTracks)
    delete fTracks;
  
  fTracks = new AliHLTTrackArray("AliHLTModelTrack");
  
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliHLTHoughTrack *track = (AliHLTHoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      AliHLTModelTrack *mtrack = (AliHLTModelTrack*)fTracks->NextTrack();
      mtrack->Init(slice,patch);
      mtrack->SetTgl(track->GetTgl());
      mtrack->SetRowRange(rowrange[0],rowrange[1]);
      for(Int_t j=fRowMin; j<=fRowMax; j++)
	{
	  Float_t hit[3];
	  track->GetLineCrossingPoint(j,hit);
	  hit[0] += AliHLTTransform::Row2X(track->GetFirstRow());
	  Float_t r = sqrt(hit[0]*hit[0] + hit[1]*hit[1]);
	  hit[2] = r*track->GetTgl();
	  Int_t se,ro;
	  AliHLTTransform::Slice2Sector(slice,j,se,ro);
	  AliHLTTransform::Local2Raw(hit,se,ro);
	  if(hit[1]<0 || hit[1]>=AliHLTTransform::GetNPads(j) || hit[2]<0 || hit[2]>=AliHLTTransform::GetNTimeBins())
	    {
	      mtrack->SetPadHit(j,-1);
	      mtrack->SetTimeHit(j,-1);
	      continue;
	    }
	  mtrack->SetPadHit(j,hit[1]);
	  mtrack->SetTimeHit(j,hit[2]);
	  mtrack->SetCrossingAngleLUT(j,fabs(track->GetPsiLine() - AliHLTTransform::Pi()/2));
	  //if(mtrack->GetCrossingAngleLUT(j) > AliHLTTransform::Deg2Rad(20))
	  //  cout<<"Angle "<<mtrack->GetCrossingAngleLUT(j)<<" psiline "<<track->GetPsiLine()*180/3.1415<<endl;
	  mtrack->CalculateClusterWidths(j);
	}
    }
  //  cout<<"Copied "<<fTracks->GetNTracks()<<" tracks "<<endl;
}

void AliHLTClusterFitter::Init(Int_t slice,Int_t patch)
{
  // Initialization
  fSlice=slice;
  fPatch=patch;
  
  fRowMin=AliHLTTransform::GetFirstRow(patch);
  fRowMax=AliHLTTransform::GetLastRow(patch);
  
  fFitted=fFailed=0;
  
  Int_t ntimes = AliHLTTransform::GetNTimeBins()+1;
  Int_t npads = AliHLTTransform::GetNPads(AliHLTTransform::GetLastRow(fPatch))+1;//Max num of pads.
  Int_t bounds = ntimes*npads;
  if(fRow)
    delete [] fRow;
  fRow = new Digit[bounds];
  if(fTracks)
    delete fTracks;
  fTracks = new AliHLTTrackArray("AliHLTModelTrack");  

}

void AliHLTClusterFitter::LoadLocalSegments()
{
  // loads local segments
  Char_t filename[1024];
  sprintf(filename,"%s/hough/tracks_ho_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);
  AliHLTMemHandler mem;
  mem.SetBinaryInput(filename);
  mem.Binary2TrackArray(fTracks);
  mem.CloseBinaryInput();
  
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliHLTModelTrack *track = (AliHLTModelTrack*)fTracks->GetCheckedTrack(i);
      if(!track) continue;

      track->CalculateHelix();
      
      track->Init(fSlice,fPatch);

      for(Int_t j=fRowMin; j<=fRowMax; j++)
	{
	  //Calculate the crossing point between track and padrow
	  
	  Float_t xyzCross[3];
	  if(!track->GetCrossingPoint(j,xyzCross))
	    continue;
	  
	  Int_t sector,row;
	  AliHLTTransform::Slice2Sector(fSlice,j,sector,row);
	  AliHLTTransform::Local2Raw(xyzCross,sector,row);
	  
	  if(xyzCross[1] < 0 || xyzCross[1] >= AliHLTTransform::GetNPads(j) ||
	     xyzCross[2] < 0 || xyzCross[2] >= AliHLTTransform::GetNTimeBins()) //track goes out of range
	    continue;
	  
	  track->SetPadHit(j,xyzCross[1]);
	  track->SetTimeHit(j,xyzCross[2]);

	  Float_t crossingangle = track->GetCrossingAngle(j);
	  track->SetCrossingAngleLUT(j,crossingangle);
	  track->CalculateClusterWidths(j);
	  track->GetClusterModel(j)->fSlice = fSlice;
	  
	}
    }
}

void AliHLTClusterFitter::LoadSeeds(Int_t *rowrange,Bool_t offline,Int_t eventnr,Float_t zvertex)
{
  //Function assumes _global_ tracks written to a single file.
  cout<<"Loading the seeds"<<endl;
  Char_t fname[1024];
  fEvent = eventnr;
  
  if(offline)
    sprintf(fname,"%s/offline/tracks_%d.raw",fPath,fEvent);
  else
    sprintf(fname,"%s/hough/tracks_%d.raw",fPath,fEvent);
  
  cout<<"AliHLTClusterFitter::LoadSeeds : Loading input tracks from "<<fname<<endl;
  
  AliHLTMemHandler tfile;
  tfile.SetBinaryInput(fname);
  
  if(fSeeds)
    delete fSeeds;
  fSeeds = new AliHLTTrackArray("AliHLTModelTrack");

  tfile.Binary2TrackArray(fSeeds);
  tfile.CloseBinaryInput();

  //if(!offline)
  //fSeeds->QSort();
  
  Int_t clustercount=0;
  for(Int_t i=0; i<fSeeds->GetNTracks(); i++)
    {
      AliHLTModelTrack *track = (AliHLTModelTrack*)fSeeds->GetCheckedTrack(i);
      if(!track) continue;

      if(!offline)
	{
	  if(i==0) cerr<<"AliHLTClusterFitter::LoadSeeds : Cutting on pt of 4 GeV!!"<<endl;
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
      /*
      if(i==0) cerr<<"Cluster fitter only in HALF TPC!!!"<<endl;
      if(origslice > 17) 
	{
	  fSeeds->Remove(i);
	  continue;
	}
      */
      track->Init(origslice,-1);
      Int_t slice = origslice;
      
      //for(Int_t j=rowrange[1]; j>=rowrange[0]; j--)
      for(Int_t j=rowrange[0]; j<=rowrange[1]; j++)
	{
	  
	  //Calculate the crossing point between track and padrow
	  Float_t angle = 0; //Perpendicular to padrow in local coordinates
	  AliHLTTransform::Local2GlobalAngle(&angle,slice);
	  if(!track->CalculateReferencePoint(angle,AliHLTTransform::Row2X(j)))
	    {
	      //cerr<<"No crossing in slice "<<slice<<" padrow "<<j<<endl;
	      continue;
	      //track->Print();
	      //exit(5);
	    }
	  Float_t xyzCross[3] = {track->GetPointX(),track->GetPointY(),track->GetPointZ()};
	  xyzCross[2] += zvertex;
	  
	  Int_t sector,row;
	  AliHLTTransform::Slice2Sector(slice,j,sector,row);
	  AliHLTTransform::Global2Raw(xyzCross,sector,row);
	  //cout<<"Examining slice "<<slice<<" row "<<j<<" pad "<<xyzCross[1]<<" time "<<xyzCross[2]<<endl;
	  if(xyzCross[1] < 0 || xyzCross[1] >= AliHLTTransform::GetNPads(j)) //Track leaves the slice
	    {
	    newslice:
	      
	      Int_t tslice=slice;
	      Float_t lastcross=xyzCross[1];
	      if(xyzCross[1] > 0)
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
	      AliHLTTransform::Local2GlobalAngle(&angle,slice);
	      if(!track->CalculateReferencePoint(angle,AliHLTTransform::Row2X(j)))
		{
		  cerr<<"No crossing in slice "<<slice<<" padrow "<<j<<endl;
		  continue;
		  //track->Print();
		  //exit(5);
		}
	      xyzCross[0] = track->GetPointX();
	      xyzCross[1] = track->GetPointY();
	      xyzCross[2] = track->GetPointZ();
	      xyzCross[2] += zvertex;
	      Int_t sector,row;
	      AliHLTTransform::Slice2Sector(slice,j,sector,row);
	      AliHLTTransform::Global2Raw(xyzCross,sector,row);
	      if(xyzCross[1] < 0 || xyzCross[1] >= AliHLTTransform::GetNPads(j)) //track is in the borderline
		{
		  if(xyzCross[1] > 0 && lastcross > 0 || xyzCross[1] < 0 && lastcross < 0)
		    goto newslice;
		  else
		    {
		      slice = tslice;//Track is on the border of two slices
		      continue;
		    }
		}
	    }
	  
	  if(xyzCross[2] < 0 || xyzCross[2] >= AliHLTTransform::GetNTimeBins())//track goes out of range
	    continue;
	  
	  if(xyzCross[1] < 0 || xyzCross[1] >= AliHLTTransform::GetNPads(j))
	    {
	      cerr<<"Slice "<<slice<<" padrow "<<j<<" pad "<<xyzCross[1]<<" time "<<xyzCross[2]<<endl;
	      track->Print();
	      exit(5);
	    }
	  
	  track->SetPadHit(j,xyzCross[1]);
	  track->SetTimeHit(j,xyzCross[2]);
	  angle=0;
	  AliHLTTransform::Local2GlobalAngle(&angle,slice);
	  Float_t crossingangle = track->GetCrossingAngle(j,slice);
	  track->SetCrossingAngleLUT(j,crossingangle);
	  
	  track->CalculateClusterWidths(j);
	  
	  track->GetClusterModel(j)->fSlice = slice;
	  
	}
      memset(track->GetHitNumbers(),0,159*sizeof(UInt_t));//Reset the hitnumbers
      track->SetNHits(0);
    }
  fSeeds->Compress();
  
  cout<<"Loaded "<<fSeeds->GetNTracks()<<" seeds and "<<clustercount<<" clusters"<<endl;
}

void AliHLTClusterFitter::FindClusters()
{
  // finds clusters
  if(!fTracks)
    {
      cerr<<"AliHLTClusterFitter::Process : No tracks"<<endl;
      return;
    }
  if(!fRowData)
    {
      cerr<<"AliHLTClusterFitter::Process : No data "<<endl;
      return;
    }
  
  AliHLTDigitRowData *rowPt = fRowData;
  AliHLTDigitData *digPt=0;

  Int_t pad,time;
  Short_t charge;
  
  if(fRowMin < 0)
    {
      fRowMin = AliHLTTransform::GetFirstRow(fPatch);
      fRowMax = AliHLTTransform::GetLastRow(fPatch);
    }
  for(Int_t i=AliHLTTransform::GetFirstRow(fPatch); i<=AliHLTTransform::GetLastRow(fPatch); i++)
    {
      if((Int_t)rowPt->fRow < fRowMin)
	{
	  AliHLTMemHandler::UpdateRowPointer(rowPt);
	  continue;
	}
      else if((Int_t)rowPt->fRow > fRowMax)
	break;
      else if((Int_t)rowPt->fRow != i)
	{
	  cerr<<"AliHLTClusterFitter::FindClusters : Mismatching row numbering "<<i<<" "<<rowPt->fRow<<endl;
	  exit(5);
	}
      fCurrentPadRow = i;
      memset((void*)fRow,0,(AliHLTTransform::GetNTimeBins()+1)*(AliHLTTransform::GetNPads(i)+1)*sizeof(Digit));
      digPt = (AliHLTDigitData*)rowPt->fDigitData;

      for(UInt_t j=0; j<rowPt->fNDigit; j++)
	{
	  pad = digPt[j].fPad;
	  time = digPt[j].fTime;
	  charge = digPt[j].fCharge;
	  fRow[(AliHLTTransform::GetNTimeBins()+1)*pad+time].fCharge = charge;
	  fRow[(AliHLTTransform::GetNTimeBins()+1)*pad+time].fUsed = kFALSE;
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
	      AliHLTModelTrack *track = (AliHLTModelTrack*)fProcessTracks->GetCheckedTrack(k);
	      if(!track) continue;
	      
	      if(fSeeding)
		if(track->GetClusterModel(i)->fSlice != fSlice) continue;
	      
	      if(track->GetPadHit(i) < 0 || track->GetPadHit(i) > AliHLTTransform::GetNPads(i)-1 ||
		 track->GetTimeHit(i) < 0 || track->GetTimeHit(i) > AliHLTTransform::GetNTimeBins()-1)
		{
		  track->SetCluster(i,0,0,0,0,0,0);
		  continue;
		}
	      
	      if(CheckCluster(k) == kFALSE)
		fFailed++;
	    }
	}
      AliHLTMemHandler::UpdateRowPointer(rowPt);
    }
  
  fSeeding = kTRUE;
  AddClusters();
  fSeeding = kFALSE;
  AddClusters();
    
  cout<<"Fitted "<<fFitted<<" clusters, failed "<<fFailed<<endl;
  cout<<"Distribution:"<<endl;
  cout<<"Bad fit "<<fgBadFitError<<endl;
  cout<<"Fit error "<<fgFitError<<endl;
  cout<<"Result error "<<fgResultError<<endl;
  cout<<"Fit range error "<<fgFitRangeError<<endl;

}

Bool_t AliHLTClusterFitter::CheckCluster(Int_t trackindex)
{
  //Check if this is a single or overlapping cluster
  
  AliHLTModelTrack *track = (AliHLTModelTrack*)fProcessTracks->GetCheckedTrack(trackindex);
  
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
	cout<<"Failed to fit cluster at row "<<row<<" pad "<<(Int_t)rint(track->GetPadHit(row))<<" time "
	    <<(Int_t)rint(track->GetTimeHit(row))<<" hitcharge "
	    <<fRow[(AliHLTTransform::GetNTimeBins()+1)*(Int_t)rint(track->GetPadHit(row))+(Int_t)rint(track->GetTimeHit(row))].fCharge<<endl;
      fgFitRangeError++;
      return kFALSE;
    }

  //Check if any other track contributes to this cluster:

  for(Int_t t=trackindex+1; t<fProcessTracks->GetNTracks(); t++)
    {
      AliHLTModelTrack *tr = (AliHLTModelTrack*)fProcessTracks->GetCheckedTrack(t);
      if(!tr) continue;
      if(fSeeding)
	if(tr->GetClusterModel(row)->fSlice != fSlice) continue;

      if(tr->GetPadHit(row) > padr[0] && tr->GetPadHit(row) < padr[1] &&
	 tr->GetTimeHit(row) > timer[0] && tr->GetTimeHit(row) < timer[1])
	{
	  if(SetFitRange(tr,padr,timer))
	    track->SetOverlap(row,t);
	}
      /*
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
	    {
	      track->SetOverlap(row,t);    //Set overlap
	    }
	}
      */
    }

  if(fDebug)
    cout<<"Fitting cluster with "<<track->GetNOverlaps(fCurrentPadRow)<<" overlaps"<<endl;

  FitClusters(track,padr,timer);
  //CalculateWeightedMean(track,padr,timer);
  return kTRUE;
}

Bool_t AliHLTClusterFitter::SetFitRange(AliHLTModelTrack *track,Int_t *padrange,Int_t *timerange)
{
  // sets the fit range
  Int_t row = fCurrentPadRow;
  Int_t nt = AliHLTTransform::GetNTimeBins()+1;
  
  Int_t nsearchbins=0;
  if(row < AliHLTTransform::GetNRowLow())
    nsearchbins=25;
  else if(row < AliHLTTransform::GetNRowLow() + AliHLTTransform::GetNRowUp1())
    nsearchbins=25;
  else
    nsearchbins=25;
  
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
  
  Int_t xyw = 3;
  Int_t zw =  5;
  
  if(padmax>=0 && timemax>=0)
    {
      //Set the hit to the local maxima of the cluster.
      //Store the maxima in the cluster model structure,
      //-only temporary, it will be overwritten when calling SetCluster.
      
      track->GetClusterModel(row)->fDPad = padmax;
      track->GetClusterModel(row)->fDTime = timemax;
      
      Int_t i=padmax,j=timemax;
      for(Int_t pdir=-1; pdir<=1; pdir+=2)
	{
	  i=padmax;
	  while(abs(padmax-i) < xyw)
	    {
	      Bool_t chargeonpad=kFALSE;
	      for(Int_t tdir=-1; tdir<=1; tdir+=2)
		{
		  j=timemax;
		  while(abs(timemax-j) < zw)
		    {
		      if(i<0 || i>=AliHLTTransform::GetNPads(row) || j<0 || j>=AliHLTTransform::GetNTimeBins()) break;
		      if(fRow[nt*i+j].fCharge)
			{
			  if(i < padrange[0]) padrange[0]=i;
			  if(i > padrange[1]) padrange[1]=i;
			  if(j < timerange[0]) timerange[0]=j;
			  if(j > timerange[1]) timerange[1]=j;
			  chargeonpad=kTRUE;
			}
		      else
			break;
		      j+=tdir;
		    } 
		}
	      if(!chargeonpad)
		break;
	      i+=pdir;
	    }
	}
      /*
	for(Int_t i=padmax-xyw; i<=padmax+xyw; i++)
	{
	for(Int_t j=timemax-zw; j<=timemax+zw; j++)
	{
	if(i<0 || i>=AliHLTTransform::GetNPads(row) || j<0 || j>=AliHLTTransform::GetNTimeBins()) continue;
	if(fRow[nt*i+j].fCharge)
	{
	if(i < padrange[0]) padrange[0]=i;
	if(i > padrange[1]) padrange[1]=i;
	if(j < timerange[0]) timerange[0]=j;
	if(j > timerange[1]) timerange[1]=j;
	}
	}
	}
      */
      
      if(fDebug)
	cout<<"New padrange "<<padrange[0]<<" "<<padrange[1]<<" "<<" time "<<timerange[0]<<" "<<timerange[1]<<endl;
      return kTRUE;
    }
  return kFALSE;
}

Bool_t AliHLTClusterFitter::IsMaximum(Int_t pad,Int_t time)
{
  // checks the maximum
  if(pad<0 || pad >= AliHLTTransform::GetNPads(fCurrentPadRow) ||
     time<0 || time >= AliHLTTransform::GetNTimeBins())
    return kFALSE;
  Int_t nt = AliHLTTransform::GetNTimeBins()+1;
  if(fRow[nt*pad+time].fUsed == kTRUE) return kFALSE; //Peak has been assigned before
  Int_t charge = fRow[nt*pad+time].fCharge;
  if(charge == 1024 || charge==0) return kFALSE;
  //if(charge == 0) return kFALSE;
  
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

void AliHLTClusterFitter::CalculateWeightedMean(AliHLTModelTrack *track,Int_t *padrange,Int_t *timerange)
{
  // calculates weighted mean
  Float_t sum=0;
  Int_t npads=0;
  Float_t pad=0,time=0;
  Int_t nt = AliHLTTransform::GetNTimeBins()+1;
  for(Int_t i=padrange[0]; i<=padrange[1]; i++)
    {
      Int_t lsum=0;
      for(Int_t j=timerange[0]; j<=timerange[1]; j++)
	{
	  lsum += fRow[nt*(i-1)+(j-1)].fCharge;
	  time += j*fRow[nt*(i-1)+(j-1)].fCharge;
	}
      if(lsum)
	npads++;
      pad += i * lsum;
    }
  if(sum)
    {
      pad /= sum;
      time /= sum;
      track->SetCluster(fCurrentPadRow,pad,time,sum,0,0,npads);
    }
  else
    track->SetCluster(fCurrentPadRow,0,0,0,0,0,0);
}

void AliHLTClusterFitter::FitClusters(AliHLTModelTrack *track,Int_t *padrange,Int_t *timerange)
{
  //Handle single and overlapping clusters
  
  /*
  if( (Int_t)rint(track->GetClusterModel(fCurrentPadRow)->fDPad) <= 1 || 
      (Int_t)rint(track->GetClusterModel(fCurrentPadRow)->fDPad) >= AliHLTTransform::GetNPads(fCurrentPadRow)-2 ||
      (Int_t)rint(track->GetClusterModel(fCurrentPadRow)->fDTime) <= 1 ||
      (Int_t)rint(track->GetClusterModel(fCurrentPadRow)->fDTime) >= AliHLTTransform::GetNTimeBins()-2)
    {
      CalculateWeightedMean(track,padrange,timerange);
      return;
    }
  */
  Int_t size = FIT_PTS;
  Int_t maxTracks = FIT_MAXPAR/NUM_PARS;
  if(track->GetNOverlaps(fCurrentPadRow) > maxTracks)
    {
      cerr<<"AliHLTClusterFitter::FitOverlappingClusters : Too many overlapping tracks"<<endl;
      return;
    }
  Int_t *overlaps = track->GetOverlaps(fCurrentPadRow);
  
  //Check if at least one cluster is not already fitted
  Bool_t allFitted=kTRUE;
  
  Int_t k=-1;
  while(k < track->GetNOverlaps(fCurrentPadRow))
    {
      AliHLTModelTrack *tr=0;
      if(k==-1)
 	tr = track;
      else
	tr = (AliHLTModelTrack*)fProcessTracks->GetCheckedTrack(overlaps[k]);
      k++;
      if(!tr) continue;
      if(!tr->IsSet(fCurrentPadRow) && !tr->IsPresent(fCurrentPadRow))//cluster has not been set and is not present
	{
	  allFitted = kFALSE;
	  break;
	}
    }
  if(allFitted)
    {
      if(fDebug)
	cout<<"But all the clusters were already fitted on row "<<fCurrentPadRow<<endl;
      return;
    }
  
  //Allocate fit parameters array; this is interface to the C code
  plane = new DPOINT[FIT_PTS];
  memset(plane,0,FIT_PTS*sizeof(DPOINT));

  Double_t x[FIT_PTS],y[FIT_PTS],s[FIT_PTS];
  
  //Fill the fit parameters:
  Double_t a[FIT_MAXPAR];
  Int_t lista[FIT_MAXPAR];
  Double_t dev[FIT_MAXPAR],chisqF;
  
  Int_t fitPars=0;
  
  Int_t nOverlaps=0;
  k=-1;
  
  //Fill the overlapping tracks:
  while(k < track->GetNOverlaps(fCurrentPadRow))
    {
      AliHLTModelTrack *tr=0;
      if(k==-1)
	tr = track;
      else
	tr = (AliHLTModelTrack*)fProcessTracks->GetCheckedTrack(overlaps[k]);
      k++;
      if(!tr) continue;
      
      if(tr->IsSet(fCurrentPadRow) && !tr->IsPresent(fCurrentPadRow)) continue;//Cluster fit failed before
      
      //Use the local maxima as the input to the fitting routine.
      //The local maxima is temporary stored in the cluster model:
      Int_t hitpad = (Int_t)rint(tr->GetClusterModel(fCurrentPadRow)->fDPad);  
      Int_t hittime = (Int_t)rint(tr->GetClusterModel(fCurrentPadRow)->fDTime);
      Int_t charge = fRow[(AliHLTTransform::GetNTimeBins()+1)*hitpad + hittime].fCharge;
      
      if(fDebug)
	cout<<"Fitting track cluster, pad "<<tr->GetPadHit(fCurrentPadRow)<<" time "
	    <<tr->GetTimeHit(fCurrentPadRow)<<" charge "<<charge<<" at local maxima in pad "<<hitpad
	    <<" time "<<hittime<<" xywidth "<<sqrt(tr->GetParSigmaY2(fCurrentPadRow))
	    <<" zwidth "<<sqrt(tr->GetParSigmaZ2(fCurrentPadRow))<<endl;
      
      if(charge==0)
	{
	  cerr<<"Charge still zero!"<<endl;
	  exit(5);
	}
            
      a[nOverlaps*NUM_PARS+2] = hitpad;
      a[nOverlaps*NUM_PARS+4] = hittime;
      
      if(!tr->IsSet(fCurrentPadRow)) //Cluster is not fitted before
	{
	  a[nOverlaps*NUM_PARS+1] = charge;
	  a[nOverlaps*NUM_PARS+3] = sqrt(tr->GetParSigmaY2(fCurrentPadRow)) * GetYWidthFactor();
	  a[nOverlaps*NUM_PARS+5] = sqrt(tr->GetParSigmaZ2(fCurrentPadRow)) * GetZWidthFactor();
	  //a[nOverlaps*NUM_PARS+6] = sqrt(tr->GetParSigmaZ2(fCurrentPadRow)) * GetZWidthFactor();
	  lista[nOverlaps*NUM_PARS + 1] = 1;
	  lista[nOverlaps*NUM_PARS + 2] = 1;
	  lista[nOverlaps*NUM_PARS + 3] = 0;
	  lista[nOverlaps*NUM_PARS + 4] = 1;
	  lista[nOverlaps*NUM_PARS + 5] = 0;
	  //lista[nOverlaps*NUM_PARS + 6] = 0;
	  fitPars             += 3;          //<-------------------
	}
      else  //Cluster was fitted before
	{
	  if(!tr->IsPresent(fCurrentPadRow))
	    {
	      cerr<<"AliHLTClusterFitter::FindClusters : Cluster not present; there is a bug here"<<endl;
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
	    cout<<endl<<"Cluster had been fitted before, pad "<<pad<<" time "<<time<<" charge "<<charge<<" width "<<xywidth<<" "<<zwidth<<endl;
	  
	  a[nOverlaps*NUM_PARS+2] = pad;
	  a[nOverlaps*NUM_PARS+4] = time;
	  a[nOverlaps*NUM_PARS+1] = charge;
	  a[nOverlaps*NUM_PARS+3] = xywidth * GetYWidthFactor();
	  a[nOverlaps*NUM_PARS+5] = zwidth * GetZWidthFactor();
	  //a[nOverlaps*NUM_PARS+6] = zwidth * GetZWidthFactor();

	  lista[nOverlaps*NUM_PARS + 1] = 1;
	  lista[nOverlaps*NUM_PARS + 2] = 0;
	  lista[nOverlaps*NUM_PARS + 3] = 0;
	  lista[nOverlaps*NUM_PARS + 4] = 0;
	  lista[nOverlaps*NUM_PARS + 5] = 0;
	  //lista[nOverlaps*NUM_PARS + 6] = 0;
	  fitPars             += 1;
	}
      nOverlaps++;
    }
  
  if(nOverlaps==0) //No clusters here
    {
      delete [] plane;
      return;
    }

  Int_t padNum=0;
  Int_t timeNumMax=0;
  Int_t ndata=0;
  Int_t totCharge=0;
  if(fDebug)
    cout<<"Padrange "<<padrange[0]<<" "<<padrange[1]<<" timerange "<<timerange[0]<<" "<<timerange[1]<<endl;
  for(Int_t i=padrange[0]; i<=padrange[1]; i++)
    {
      Int_t maxCharge = 0;
      Int_t timeNum=0;
      for(Int_t j=timerange[0]; j<=timerange[1]; j++)
	{
	  Int_t charge = fRow[(AliHLTTransform::GetNTimeBins()+1)*i + j].fCharge;
	  
	  if(charge <= 0) continue;

	  timeNum++;
	  if(charge > maxCharge)
	    {
	      maxCharge = charge;
	      //timeNum++;
	    }
	  if(fDebug)
	    cout<<"Filling padrow "<<fCurrentPadRow<<" pad "<<i<<" time "<<j<<" charge "<<charge<<endl;
	  totCharge += charge;
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
      if(maxCharge) //there was charge on this pad
	padNum++;
      if(timeNumMax < timeNum)
	timeNumMax = timeNum;
    }

  if(padNum <= 1 || timeNumMax <=1 || nOverlaps > fNmaxOverlaps || ndata <= fitPars) //too few to do fit
    {
      SetClusterfitFalse(track);
      if(fDebug)
	cout<<"Too few digits or too many overlaps: "<<padNum<<" "<<timeNumMax<<" "<<nOverlaps<<" ndata "<<ndata<<" fitPars "<<fitPars<<endl;
      delete [] plane;
      return;
    }

  
  Int_t npars = nOverlaps * NUM_PARS;
  if(fDebug)
    cout<<"Number of overlapping clusters "<<nOverlaps<<endl;
  Int_t ret = lev_marq_fit( x, y, s, ndata, a, lista, dev, npars, &chisqF, f2gauss5 );
  
  if(ret<0)
    {
      SetClusterfitFalse(track);
      fFailed++;
      fgFitError++;
      delete [] plane;
      return;
      //exit(5);
    }

  chisqF /= (ndata-fitPars);
  if(fDebug)
    cout<<"Chisq "<<chisqF<<endl;
  
  Bool_t overlapping=kFALSE;
  if(track->GetNOverlaps(fCurrentPadRow) > 0)//There is a overlap
    overlapping=kTRUE;
  
  k=-1;
  nOverlaps=0;
  while(k < track->GetNOverlaps(fCurrentPadRow))
    {
      AliHLTModelTrack *tr=0;
      if(k==-1)
	tr = track;
      else
	tr = (AliHLTModelTrack*)fProcessTracks->GetCheckedTrack(overlaps[k]);
      k++;
      if(!tr) continue;
      if(!tr->IsPresent(fCurrentPadRow))
	{
	  if(tr->IsSet(fCurrentPadRow)) continue;//This cluster has been set before
	  
	  Int_t lpatch;
	  if(fCurrentPadRow < AliHLTTransform::GetNRowLow())
	    lpatch=0;
	  else if(fCurrentPadRow < AliHLTTransform::GetNRowLow() + AliHLTTransform::GetNRowUp1())
	    lpatch=1;
	  else 
	    lpatch=2;
	  
	  //if(chisqF < fChiSqMax[(Int_t)overlapping])//cluster fit is good enough
	  if(chisqF < fChiSqMax[lpatch])//cluster fit is good enough
	    {
	      totCharge = (Int_t)(2*AliHLTTransform::Pi() * a[nOverlaps*NUM_PARS+1] * a[nOverlaps*NUM_PARS+3] * a[nOverlaps*NUM_PARS+5]);
	      Float_t fpad = a[nOverlaps*NUM_PARS+2];
	      Float_t ftime = a[nOverlaps*NUM_PARS+4];
	      if(totCharge < 0 || fpad < -1 || fpad > AliHLTTransform::GetNPads(fCurrentPadRow) || 
		 ftime < -1 || ftime > AliHLTTransform::GetNTimeBins())
		{
		  if(fDebug)
		    cout<<"AliHLTClusterFitter::Fatal result(s) in fit; in slice "<<fSlice<<" row "<<fCurrentPadRow
			<<"; pad "<<fpad<<" time "<<ftime<<" charge "<<totCharge<<" xywidth "<<a[nOverlaps*NUM_PARS+3]
			<<" zwidth "<<a[nOverlaps*NUM_PARS+5]<<" peakcharge "<<a[nOverlaps*NUM_PARS+1]<<endl;
		  tr->SetCluster(fCurrentPadRow,0,0,0,0,0,0);
		  fFailed++;
		  fgResultError++;
		  continue;
		}
	      
	      tr->SetCluster(fCurrentPadRow,fpad,ftime,totCharge,0,0,padNum);
	      /*
		tr->SetCluster(fCurrentPadRow,fpad,ftime,totCharge,
		pow(a[nOverlaps*NUM_PARS+3],2),
		pow(a[nOverlaps*NUM_PARS+5],2),padNum);
	      */
	      if(fDebug)
		{
		  cout<<"Setting cluster in slice "<<fSlice<<" row "<<fCurrentPadRow<<" pad "<<a[nOverlaps*NUM_PARS+2]<<" time "<<a[nOverlaps*NUM_PARS+4]
		      <<" padwidth "<<a[nOverlaps*NUM_PARS+3]<<" timewidth "<<a[nOverlaps*NUM_PARS+5]
		      <<" charge "<<totCharge<<endl;
		}
	      /*
	      //Set the digits to used:
	      for(Int_t i=padrange[0]; i<=padrange[1]; i++)
	      for(Int_t j=timerange[0]; j<=timerange[1]; j++)
	      fRow[(AliHLTTransform::GetNTimeBins()+1)*i + j].fUsed = kTRUE;
	      */
	      fFitted++;
	    }
	  else //fit was too bad
	    {
	      if(fDebug)
		cout<<"Cluster fit was too bad"<<endl;
	      tr->SetCluster(fCurrentPadRow,0,0,0,0,0,0);
	      fgBadFitError++;
	      fFailed++;
	    }
	}
      nOverlaps++;
    }
  
  delete [] plane;
}

void AliHLTClusterFitter::SetClusterfitFalse(AliHLTModelTrack *track)
{
  //Cluster fit failed, so set the clusters to all the participating
  //tracks to zero.
  
  Int_t i=-1;
  Int_t *overlaps = track->GetOverlaps(fCurrentPadRow);
  while(i < track->GetNOverlaps(fCurrentPadRow))
    {
      AliHLTModelTrack *tr=0;
      if(i==-1)
	tr = track;
      else
	tr = (AliHLTModelTrack*)fProcessTracks->GetCheckedTrack(overlaps[i]);
      i++;
      if(!tr) continue;
      
      //Set the digit data to unused, so it can be fitted to another bastard:
      Int_t hitpad = (Int_t)rint(tr->GetClusterModel(fCurrentPadRow)->fDPad);  
      Int_t hittime = (Int_t)rint(tr->GetClusterModel(fCurrentPadRow)->fDTime);
      fRow[(AliHLTTransform::GetNTimeBins()+1)*hitpad + hittime].fUsed = kFALSE;
      
      tr->SetCluster(fCurrentPadRow,0,0,0,0,0,0);
    }

}


void AliHLTClusterFitter::AddClusters()
{
  // adds clusters
  if(!fClusters)
    {
      fClusters = new AliHLTSpacePointData[fNMaxClusters];
      fNClusters=0;
    }
  
  if(fDebug)
    cout<<"Writing cluster in slice "<<fSlice<<" patch "<<fPatch<<endl;
  
  AliHLTTrackArray *tracks=0;
  if(fSeeding==kTRUE)
    tracks = fSeeds;
  else
    tracks = fTracks;
  
  if(!tracks)
    return;
  
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliHLTModelTrack *tr = (AliHLTModelTrack*)tracks->GetCheckedTrack(i);
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

	  if(pad < -1 || pad >= AliHLTTransform::GetNPads(i) || 
	     time < -1 || time >= AliHLTTransform::GetNTimeBins())
	    {
	      continue;
//  	      cout<<"slice "<<fSlice<<" row "<<i<<" pad "<<pad<<" time "<<time<<endl;
//  	      tr->Print();
//  	      exit(5);
	    }
	  
	  tr->CalculateClusterWidths(i,kTRUE); //Parametrize errors
	  
	  tr->GetSigmaY2(i,xywidth);
	  tr->GetSigmaZ2(i,zwidth);
	  Float_t xyz[3];
	  Int_t sector,row;
	  AliHLTTransform::Slice2Sector(fSlice,i,sector,row);
	  
	  AliHLTTransform::Raw2Local(xyz,sector,row,pad,time);
	  
	  if(fNClusters >= fNMaxClusters)
	    {
	      cerr<<"AliHLTClusterFitter::AddClusters : Too many clusters "<<fNClusters<<endl;
	      exit(5);
	    }
	  
	  fClusters[fNClusters].fX = xyz[0];
	  fClusters[fNClusters].fY = xyz[1];
	  fClusters[fNClusters].fZ = xyz[2];
	  fClusters[fNClusters].fCharge = charge;
	  fClusters[fNClusters].fPadRow = i;
	  Int_t pa = AliHLTTransform::GetPatch(i);
	  if(xywidth==0 || zwidth==0)
	    cerr<<"AliHLTClusterFitter::AddClusters : Cluster with zero width"<<endl;
	  if(xywidth>0)
	    fClusters[fNClusters].fSigmaY2 = xywidth*pow(AliHLTTransform::GetPadPitchWidth(pa),2);
	  else
	    fClusters[fNClusters].fSigmaY2 = 1;
	  if(zwidth>0)
	    fClusters[fNClusters].fSigmaZ2 = zwidth*pow(AliHLTTransform::GetZWidth(),2);
	  else
	    fClusters[fNClusters].fSigmaZ2 = 1;
	  Int_t pat=fPatch;
	  if(fPatch==-1)
	    pat=0;
	  fClusters[fNClusters].fID = fNClusters + ((fSlice&0x7f)<<25)+((pat&0x7)<<22);
	  
	  if(nhits >= AliHLTTransform::GetNRows())
	    {
	      cerr<<"AliHLTClusterFitter::AddClusters : Cluster counter of out range "<<nhits<<endl;
	      exit(5);
	    }
	  hitids[nhits++] = fClusters[fNClusters].fID;
	  
#ifdef do_mc
	  Int_t trackID[3];
	  Int_t fpad = (Int_t)rint(pad);
	  Int_t ftime = (Int_t)rint(time);
	  if(fpad < 0)
	    fpad=0;
	  if(fpad >= AliHLTTransform::GetNPads(i))
	    fpad = AliHLTTransform::GetNPads(i)-1;
	  if(ftime<0)
	    ftime=0;
	  if(ftime >= AliHLTTransform::GetNTimeBins())
	    ftime = AliHLTTransform::GetNTimeBins()-1;
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

void AliHLTClusterFitter::WriteTracks(Int_t minHits)
{
  // writes tracks
  if(!fSeeds)
    return;
  
  AliHLTTrackArray *fakes = new AliHLTTrackArray();
  
  Int_t clustercount=0;
  for(Int_t i=0; i<fSeeds->GetNTracks(); i++)
    {
      AliHLTModelTrack *tr = (AliHLTModelTrack*)fSeeds->GetCheckedTrack(i);
      if(!tr) continue;
      if(tr->GetNHits() < minHits)
	{
	  fakes->AddLast(tr);
	  fSeeds->Remove(i);
	}
      clustercount += tr->GetNHits();
    }
  
  cout<<"Writing "<<clustercount<<" clusters"<<endl;
  fSeeds->Compress();
  AliHLTMemHandler mem;
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

void AliHLTClusterFitter::WriteClusters(Bool_t global)
{
  // writes clusters
  AliHLTMemHandler mem;
  if(fDebug)
    cout<<"Write "<<fNClusters<<" clusters to file"<<endl;
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
