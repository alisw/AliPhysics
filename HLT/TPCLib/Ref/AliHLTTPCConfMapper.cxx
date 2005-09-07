// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group
 
#include "AliHLTTPCStandardIncludes.h"
#include <sys/time.h>

#include "AliHLTTPCLogging.h" 
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCConfMapTrack.h"
#include "AliHLTTPCConfMapPoint.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCConfMapper.h"

/** \class AliHLTTPCConfMapper
<pre>
//_____________________________________________________________
// AliHLTTPCConfMapper
//
// Conformal mapping base class
//
</pre>
*/

ClassImp(AliHLTTPCConfMapper)

Double_t AliHLTTPCConfMapper::pi=3.14159265358979323846;
Double_t AliHLTTPCConfMapper::twopi=2*pi;
Double_t AliHLTTPCConfMapper::todeg=180./pi;

AliHLTTPCConfMapper::AliHLTTPCConfMapper()
{
  //Default constructor
  fVertex = NULL;
  fTrack = NULL;
  fHit = NULL;
  fVolume = NULL;
  fRow = NULL;
  fBench = (Bool_t)true;
  fVertexConstraint = (Bool_t)true;
  fParamSet[0]=0;
  fParamSet[1]=0;
}


AliHLTTPCConfMapper::~AliHLTTPCConfMapper()
{
  // Destructor.

  if(fVolume) {
    delete [] fVolume;
  }
  if(fRow) {
    delete [] fRow;
  }
  if(fHit) {
    delete [] fHit;
  }
  if(fTrack) {
    delete fTrack;
  }

}
 
void AliHLTTPCConfMapper::InitVolumes()
{
  //Data organization.
  //Allocate volumes, set conformal coordinates and pointers.
  
  //Should be done after setting the track parameters
  
  fNumRowSegmentPlusOne = AliHLTTPCTransform::GetNRows();//NumRows[0]; //Maximum 32.
  fNumPhiSegmentPlusOne = fNumPhiSegment+1;
  fNumEtaSegmentPlusOne = fNumEtaSegment+1;
  fNumPhiEtaSegmentPlusOne = fNumPhiSegmentPlusOne*fNumEtaSegmentPlusOne;
  fBounds = fNumRowSegmentPlusOne * fNumPhiSegmentPlusOne * fNumEtaSegmentPlusOne;
  
  //Allocate volumes:
  if(fVolume) delete [] fVolume;
  if(fRow) delete [] fRow;
  
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCConfMapper::InitVolumes","Memory")<<AliHLTTPCLog::kDec<<
    "Allocating "<<fBounds*sizeof(AliHLTTPCConfMapContainer)<<" Bytes to fVolume"<<ENDLOG;
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCConfMapper::InitVolumes","Memory")<<AliHLTTPCLog::kDec<<
    "Allocating "<<fNumRowSegmentPlusOne*sizeof(AliHLTTPCConfMapContainer)<<" Bytes to fRow"<<ENDLOG;
  
  fVolume = new AliHLTTPCConfMapContainer[fBounds];
  fRow = new AliHLTTPCConfMapContainer[fNumRowSegmentPlusOne];
  
  memset(fVolume,0,fBounds*sizeof(AliHLTTPCConfMapContainer));
  memset(fRow,0,fNumRowSegmentPlusOne*sizeof(AliHLTTPCConfMapContainer));
  
  //Int_t max_num_of_tracks = 2000;
  //Int_t max_num_of_hits = 120000;
  Int_t max_num_of_tracks = 3500;
  Int_t max_num_of_hits = 200000;
  
  if(fHit)
    delete [] fHit;
  if(fTrack)
    delete fTrack;
    
  fHit = new AliHLTTPCConfMapPoint[max_num_of_hits];
  fTrack = new AliHLTTPCTrackArray("AliHLTTPCConfMapTrack",max_num_of_tracks);
}

void AliHLTTPCConfMapper::InitSector(Int_t sector,Int_t *rowrange,Float_t *etarange)
{ //sector means slice here
  //Initialize tracker for tracking in a given sector.
  //Resets track and hit arrays.
  //Here it is also possible to specify a subsector, by defining
  //rowrange[0]=innermost row;
  //rowrange[1]=outermostrow;
  //Finally you can specify etaslices to save time (assuming a good seed from TRD...)
    
  //Define tracking area:
  if(rowrange)
    {
      fRowMin = rowrange[0];
      fRowMax = rowrange[1];
    }
  else //complete sector
    {
      fRowMin = 0;
      fRowMax = AliHLTTPCTransform::GetNRows() - 1;
    }
  if(etarange)
    {
      fEtaMin = etarange[0];
      fEtaMax = sector < 18 ? etarange[1] : -etarange[1];
    }
  else
    {
      fEtaMin = 0;
      fEtaMax = sector < 18 ? 0.9 : -0.9;
    }
  
  //Set the angles to sector 2:
  fPhiMin = -1.*10/todeg;//fParam->GetAngle(sector) - 10/todeg;
  fPhiMax = 10/todeg;//fParam->GetAngle(sector) + 10/todeg;

  nTracks=0;
  fMainVertexTracks = 0;
  fClustersUnused = 0;
  fEtaHitsOutOfRange=0;
  fPhiHitsOutOfRange=0;
  
  fNumRowSegment = fRowMax - fRowMin; //number of rows to be considered by tracker
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCConfMapper::InitSector","B-field")
    <<"Tracker initializing with a magnetic field of "<<AliHLTTPCTransform::GetBField()<<ENDLOG;
  
  fTrack->Reset();
}

Bool_t AliHLTTPCConfMapper::ReadHits(UInt_t count, AliHLTTPCSpacePointData* hits )
{
  Int_t nhit=(Int_t)count; 
  for(Int_t i=0;i<nhit;i++)
    {
      fHit[fClustersUnused].Reset();
      fHit[fClustersUnused].ReadHits(&(hits[i]));
      fClustersUnused++;    
    }
//   for (Int_t i=0;i<nhit;i++)
//     {
//       fHit[i].Reset();
//       fHit[i].ReadHits(&(hits[i]));
//     }
//   fClustersUnused += nhit;
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCConfMapper::ReadHits","#hits")
    <<AliHLTTPCLog::kDec<<"hit_counter: "<<nhit<<" count: "<<count<<ENDLOG;
  
  return true;
}

void AliHLTTPCConfMapper::SetPointers()
{
  //Check if there are not enough clusters to make a track in this sector
  //Can happen in pp events.

  if(fClustersUnused < fMinPoints[fVertexConstraint])
    return;
  
  //Reset detector volumes
  memset(fVolume,0,fBounds*sizeof(AliHLTTPCConfMapContainer));
  memset(fRow,0,fNumRowSegmentPlusOne*sizeof(AliHLTTPCConfMapContainer));
  
  Float_t phiSlice = (fPhiMax-fPhiMin)/fNumPhiSegment;
  Float_t etaSlice = (fEtaMax-fEtaMin)/fNumEtaSegment;

  Int_t volumeIndex;
  Int_t local_counter=0;
  for(Int_t j=0; j<fClustersUnused; j++)
    {
      //AliHLTTPCConfMapPoint *thisHit = (AliHLTTPCConfMapPoint*)fHit->At(j);
      AliHLTTPCConfMapPoint *thisHit = &(fHit[j]);

      thisHit->Setup(fVertex);
      
      Int_t localrow = thisHit->GetPadRow();
      
      if(localrow < fRowMin || localrow > fRowMax)
	continue;

      //Get indexes:
      thisHit->phiIndex=(Int_t)((thisHit->GetPhi()-fPhiMin)/phiSlice +1);
      
      if(thisHit->phiIndex<1 || thisHit->phiIndex>fNumPhiSegment)
	{
	  //cout << "Phiindex: " << thisHit->phiIndex << " " << thisHit->GetPhi() << endl;
	  fPhiHitsOutOfRange++;
	  continue;
	}
      
      thisHit->etaIndex=(Int_t)((thisHit->GetEta()-fEtaMin)/etaSlice + 1);
      if(thisHit->etaIndex<1 || thisHit->etaIndex>fNumEtaSegment)
	{
	  //cout << "Etaindex: " << thisHit->etaIndex << " " << thisHit->GetEta() << endl;
	  fEtaHitsOutOfRange++;
	  continue;
	}
      local_counter++;
      
      volumeIndex = (localrow-fRowMin)*fNumPhiEtaSegmentPlusOne+thisHit->phiIndex*fNumEtaSegmentPlusOne+thisHit->etaIndex;
      
      if(fVolume[volumeIndex].first == NULL)
	fVolume[volumeIndex].first = (void *)thisHit;
      else
 	((AliHLTTPCConfMapPoint *)fVolume[volumeIndex].last)->nextVolumeHit=thisHit;
      fVolume[volumeIndex].last = (void *)thisHit;
      
      
      //set row pointers
      if(fRow[(localrow-fRowMin)].first == NULL)
 	fRow[(localrow-fRowMin)].first = (void *)thisHit;
      else
 	((AliHLTTPCConfMapPoint *)(fRow[(localrow-fRowMin)].last))->nextRowHit = thisHit;
	fRow[(localrow-fRowMin)].last = (void *)thisHit;
    }
  
  if(fClustersUnused>0 && local_counter==0)
    LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapper::SetPointers","Parameters")
      <<AliHLTTPCLog::kDec<<"No points passed to track finder, hits out of range: "
      <<fEtaHitsOutOfRange+fPhiHitsOutOfRange<<ENDLOG;

  Int_t hits_accepted=fClustersUnused-(fEtaHitsOutOfRange+fPhiHitsOutOfRange);
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCConfMapper::SetPointers","Setup")
    <<"Setup finished, hits out of range: "<<fEtaHitsOutOfRange+fPhiHitsOutOfRange
    <<" hits accepted "<<hits_accepted<<ENDLOG;
}

void AliHLTTPCConfMapper::MainVertexTracking_a()
{
  //Tracking with vertex constraint.

  if(!fParamSet[(Int_t)kTRUE])
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapper::MainVertexTracking","Parameters")<<AliHLTTPCLog::kDec<<
	"Tracking parameters not set!"<<ENDLOG;
      return;
    }

  Double_t initCpuTime,cpuTime;
  initCpuTime = CpuTime();

  SetPointers();
  SetVertexConstraint(true);
  cpuTime = CpuTime() - initCpuTime;
  if(fBench)
    LOG(AliHLTTPCLog::kInformational,"AliHLTTPCConfMapper::MainVertexTracking_a","Timing")
      <<AliHLTTPCLog::kDec<<"Setup finished in "<<cpuTime*1000<<" ms"<<ENDLOG;
  
}

void AliHLTTPCConfMapper::MainVertexTracking_b()
{
  //Tracking with vertex constraint.

  if(!fParamSet[(Int_t)kTRUE])
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapper::MainVertexTracking","Parameters")<<AliHLTTPCLog::kDec<<
	"Tracking parameters not set!"<<ENDLOG;
      return;
    }
  Double_t initCpuTime,cpuTime;
  initCpuTime = CpuTime();
  
  ClusterLoop();
 
  cpuTime = CpuTime() - initCpuTime;
  if(fBench)
    LOG(AliHLTTPCLog::kInformational,"AliHLTTPCConfMapper::MainVertexTracking_b","Timing")
      <<AliHLTTPCLog::kDec<<"Main Tracking finished in "<<cpuTime*1000<<" ms"<<ENDLOG;
}

void AliHLTTPCConfMapper::MainVertexTracking()
{
  //Tracking with vertex constraint.

  if(!fParamSet[(Int_t)kTRUE])
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapper::MainVertexTracking","Parameters")<<AliHLTTPCLog::kDec<<
	"Tracking parameters not set!"<<ENDLOG;
      return;
    }

  Double_t initCpuTime,cpuTime;
  initCpuTime = CpuTime();
  
  SetPointers();
  SetVertexConstraint(true);
      
  ClusterLoop();

  cpuTime = CpuTime() - initCpuTime;
  if(fBench)
    LOG(AliHLTTPCLog::kInformational,"AliHLTTPCConfMapper::MainVertexTracking","Timing")<<AliHLTTPCLog::kDec<<
      "Tracking finished in "<<cpuTime*1000<<" ms"<<ENDLOG;
  
  return;
}

void AliHLTTPCConfMapper::NonVertexTracking()
{
  //Tracking with no vertex constraint. This should be called after doing MainVertexTracking,
  //in order to do tracking on the remaining clusters.
  //The conformal mapping is now done with respect to the first cluster
  //assosciated with this track.
  
  if(!fParamSet[(Int_t)kFALSE])
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapper::NonVertexTracking","Parameters")<<AliHLTTPCLog::kDec<<
	"Tracking parameters not set!"<<ENDLOG;
      return;
    }
  
  SetVertexConstraint(false);
  ClusterLoop();
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCConfMapper::NonVertexTracking","ntracks")<<AliHLTTPCLog::kDec<<
    "Number of nonvertex tracks found: "<<(nTracks-fMainVertexTracks)<<ENDLOG;
  return;
}

void AliHLTTPCConfMapper::MainVertexSettings(Int_t trackletlength, Int_t tracklength,
					 Int_t rowscopetracklet, Int_t rowscopetrack,
					 Double_t maxphi,Double_t maxeta)
{
  //Settings for main vertex tracking. The cuts are:
  //TrackletLength:      #hits on segment, before trying to build a track
  //TrackLength:         Minimum hits on a track
  //RowScopeTracklet:    Row search range for segments
  //RowScopeTrack:       Row search range for tracks
  
  SetTrackletLength(trackletlength,(Bool_t)true);
  SetRowScopeTracklet(rowscopetracklet, (Bool_t) true);
  SetRowScopeTrack(rowscopetrack, (Bool_t) true);
  SetMinPoints(tracklength,(Bool_t)true);
  fMaxPhi=maxphi;
  fMaxEta=maxeta;
  SetParamDone(kTRUE);
}

void AliHLTTPCConfMapper::NonVertexSettings(Int_t trackletlength, Int_t tracklength,
					Int_t rowscopetracklet, Int_t rowscopetrack)
{
  SetTrackletLength(trackletlength,(Bool_t)false);
  SetRowScopeTracklet(rowscopetracklet, (Bool_t)false);
  SetRowScopeTrack(rowscopetrack, (Bool_t)false);
  SetMinPoints(tracklength,(Bool_t)false);
  SetParamDone(kFALSE);
}

void AliHLTTPCConfMapper::SetTrackCuts(Double_t hitChi2Cut, Double_t goodHitChi2, Double_t trackChi2Cut,Int_t maxdist,Bool_t vertexconstraint)
{
  //Settings for tracks. The cuts are:
  //HitChi2Cut:     Maximum hit chi2
  //goodHitChi2:    Chi2 to stop look for next hit
  //trackChi2Cut:   Maximum track chi2
  //maxdist:        Maximum distance between two clusters when forming segments
  
  SetHitChi2Cut(hitChi2Cut,vertexconstraint);
  SetGoodHitChi2(goodHitChi2,vertexconstraint);
  SetTrackChi2Cut(trackChi2Cut,vertexconstraint);
  SetMaxDist(maxdist,vertexconstraint);
}

void AliHLTTPCConfMapper::SetTrackletCuts(Double_t maxangle,Double_t goodDist, Bool_t vertex_constraint)
{
  //Sets cuts of tracklets. Right now this is only:
  //maxangle:  Maximum angle when forming segments (if trackletlength > 2)
 
  fGoodDist=goodDist;
  SetMaxAngleTracklet(maxangle, vertex_constraint);
}

void AliHLTTPCConfMapper::ClusterLoop()
{
  //Loop over hits, starting at outermost padrow, and trying to build segments.
  
  //Check if there are not enough clusters to make a track in this sector
  //Can happen in pp events.
  if(fClustersUnused < fMinPoints[fVertexConstraint])
    return;
  
  Int_t row_segm,lastrow = fRowMin + fMinPoints[fVertexConstraint];
  AliHLTTPCConfMapPoint *hit;
  
  //Loop over rows, and try to create tracks from the hits.
  //Starts at the outermost row, and loops as long as a track can be build, due to length.
  
  for(row_segm = fRowMax; row_segm >= lastrow; row_segm--)
    {
      if(fRow[(row_segm-fRowMin)].first && ((AliHLTTPCConfMapPoint*)fRow[(row_segm-fRowMin)].first)->GetPadRow() < fRowMin + 1)
	break;

      for(hit = (AliHLTTPCConfMapPoint*)fRow[(row_segm-fRowMin)].first; hit!=0; hit=hit->nextRowHit)
	{
	  if(hit->GetUsage() == true)
	    continue;
	  else
	    CreateTrack(hit);
	}
    }
  
  return;
}


void AliHLTTPCConfMapper::CreateTrack(AliHLTTPCConfMapPoint *hit)
{
  //Tries to create a track from the initial hit given by ClusterLoop()

  AliHLTTPCConfMapPoint *closest_hit = NULL;
  AliHLTTPCConfMapTrack *track = NULL;
  
  Int_t point;
  Int_t tracks = nTracks;
  nTracks++;

  track = (AliHLTTPCConfMapTrack*)fTrack->NextTrack();

  //reset hit parameters:
  track->Reset();
  
  UInt_t *trackhitnumber = track->GetHitNumbers();
    
  //set conformal coordinates if we are looking for non vertex tracks
  if(!fVertexConstraint) 
    {
      hit->SetAllCoord(hit);
    }
  
  //fill fit parameters of initial track:
  track->UpdateParam(hit); //here the number of hits is incremented.
  trackhitnumber[track->GetNumberOfPoints()-1] = hit->GetHitNumber();
  
  Double_t dx,dy;
  //create tracklets:
  
  for(point=1; point<fTrackletLength[fVertexConstraint]; point++)
    {
      if((closest_hit = GetNextNeighbor(hit)))
	{//closest hit exist
	  
	  //   Calculate track length in sz plane
	  dx = ((AliHLTTPCConfMapPoint*)closest_hit)->GetX() - ((AliHLTTPCConfMapPoint*)hit)->GetX();
	  dy = ((AliHLTTPCConfMapPoint*)closest_hit)->GetY() - ((AliHLTTPCConfMapPoint*)hit)->GetY();
	  //track->fLength += (Double_t)sqrt ( dx * dx + dy * dy ) ;
	  Double_t length = track->GetLength()+(Double_t)sqrt ( dx * dx + dy * dy );
	  track->SetLength(length);

	  //closest_hit->SetS(track->fLength);
	  closest_hit->SetS(track->GetLength());

	  //update fit parameters
	  track->UpdateParam(closest_hit);
	  trackhitnumber[track->GetNumberOfPoints()-1] = closest_hit->GetHitNumber();
	
	  hit = closest_hit;
	}
      else
	{
	  //closest hit does not exist:
	  track->DeleteCandidate();
	  fTrack->RemoveLast();
	  nTracks--;
	  point = fTrackletLength[fVertexConstraint];
	}
    }
  
  //tracklet is long enough to be extended to a track
  if(track->GetNumberOfPoints() == fTrackletLength[fVertexConstraint])
    {
      
      track->SetProperties(true);
            
      if(TrackletAngle(track) > fMaxAngleTracklet[fVertexConstraint])
	{//proof if the first points seem to be a beginning of a track
	  track->SetProperties(false);
	  track->DeleteCandidate();
	  fTrack->RemoveLast();
	  nTracks--;
	}
      
      else//good tracklet ->proceed, follow the trackfit
	{
	  tracks++;
	  			  
	  //define variables to keep the total chi:
	  Double_t xyChi2 = track->fChiSq[0];
	  Double_t szChi2 = track->fChiSq[1];
	  
	  for(point = fTrackletLength[fVertexConstraint]; point <= fNumRowSegment; point++)
	    {
	      track->fChiSq[0] = fHitChi2Cut[fVertexConstraint];
	      closest_hit = GetNextNeighbor((AliHLTTPCConfMapPoint*)track->lastHit,track);
	      
	      if(closest_hit)
		{
		  
		  //keep total chi:
		  Double_t lxyChi2 = track->fChiSq[0]-track->fChiSq[1];
		  xyChi2 += lxyChi2;
		  closest_hit->xyChi2 = lxyChi2;
		  		  
		  //update track length:
		  //track->fLength = closest_hit->GetS();
		  track->SetLength(closest_hit->GetS());
		  szChi2 += track->fChiSq[1];
		  closest_hit->szChi2 = track->fChiSq[1];
		  
		  track->UpdateParam(closest_hit);
		  trackhitnumber[track->GetNumberOfPoints()-1] = closest_hit->GetHitNumber();
		  
		  //add closest hit to track
		  closest_hit->SetUsage(true);
		  closest_hit->SetTrackNumber(tracks-1);
		  
		}//closest_hit
	      
	      else
		{
		  //closest hit does not exist
		  point = fNumRowSegment; //continue with next hit in segment
		}//else
	      
	    }//create tracks
	  
	  //store track chi2:
	  track->fChiSq[0] = xyChi2;
	  track->fChiSq[1] = szChi2;
	  Double_t normalized_chi2 = (track->fChiSq[0]+track->fChiSq[1])/track->GetNumberOfPoints();
	  
	  //remove tracks with not enough points already now
	  if(track->GetNumberOfPoints() < fMinPoints[fVertexConstraint] || normalized_chi2 > fTrackChi2Cut[fVertexConstraint])
	    {
	      track->SetProperties(false);
	      nTracks--;
	      track->DeleteCandidate();
	      fTrack->RemoveLast();
	      tracks--;
	    }
	  
	  else
	    {
	      fClustersUnused -= track->GetNumberOfPoints();
	      track->ComesFromMainVertex(fVertexConstraint);
	      //mark track as main vertex track or not
	      track->SetSector(2); //only needed for testing purposes.
	      track->SetRowRange(fRowMin,fRowMax);

	      if(fVertexConstraint) 
		fMainVertexTracks++;
	    }
     
	}//good tracklet
      
    }
  
  return;
}

AliHLTTPCConfMapPoint *AliHLTTPCConfMapper::GetNextNeighbor(AliHLTTPCConfMapPoint *start_hit,
					  AliHLTTPCConfMapTrack *track)
{
  //When forming segments: Finds closest hit to input hit
  //When forming tracks: Find closest hit to track fit.
  
  Double_t dist,closest_dist = fMaxDist[fVertexConstraint];
  
  AliHLTTPCConfMapPoint *hit = NULL;
  AliHLTTPCConfMapPoint *closest_hit = NULL;
    
  Int_t sub_row_segm;
  Int_t sub_phi_segm;
  Int_t sub_eta_segm;
  Int_t volumeIndex;
  Int_t test_hit;

  Int_t max_row = start_hit->GetPadRow()-1;
  Int_t min_row;

  if(track) //finding hit close to trackfit
    {
      min_row = start_hit->GetPadRow()-fRowScopeTrack[fVertexConstraint];
    }
  else
    {
      min_row = start_hit->GetPadRow()-fRowScopeTracklet[fVertexConstraint];
    }

  //make a smart loop
  Int_t loop_eta[25] = {0,0,0,-1,-1,-1,1,1,1, 0,0,-1,-1,1,1,-2,-2,-2,-2,-2,2,2,2,2,2};
  Int_t loop_phi[25] = {0,-1,1,0,-1,1,0,-1,1, -2,2,-2,2,-2,2,-2,-1,0,1,2,-2,-1,0,1,2};
  
  if(min_row < fRowMin)
    min_row = fRowMin;
  if(max_row < fRowMin)
    return 0;  //reached the last padrow under consideration

  else
    {
      //loop over sub rows
      for(sub_row_segm=max_row; sub_row_segm>=min_row; sub_row_segm--)
	{
	  //loop over subsegments, in the order defined above.
	  for(Int_t i=0; i<9; i++)  
	    {
	      sub_phi_segm = start_hit->phiIndex + loop_phi[i];
	      
	      if(sub_phi_segm < 0 || sub_phi_segm >= fNumPhiSegment)
		continue;
	      /*
		if(sub_phi_segm<0)
		sub_phi_segm += fNumPhiSegment;
		
		else if(sub_phi_segm >=fNumPhiSegment)
		sub_phi_segm -= fNumPhiSegment;
	      */
	      //loop over sub eta segments
	      
	      sub_eta_segm = start_hit->etaIndex + loop_eta[i];
	      
	      if(sub_eta_segm < 0 || sub_eta_segm >=fNumEtaSegment)
		continue;//segment exceeds bounds->skip it
	      
	      //loop over hits in this sub segment:
	      volumeIndex= (sub_row_segm-fRowMin)*fNumPhiEtaSegmentPlusOne +
		sub_phi_segm*fNumEtaSegmentPlusOne + sub_eta_segm;
	      
	      if(volumeIndex<0)
		{//debugging
		  LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapper::GetNextNeighbor","Memory")<<AliHLTTPCLog::kDec<<
		    "VolumeIndex error "<<volumeIndex<<ENDLOG;
		}
	      
	      for(hit = (AliHLTTPCConfMapPoint*)fVolume[volumeIndex].first;
		  hit!=0; hit = hit->nextVolumeHit)
		{
		  
		  if(!hit->GetUsage())
		    {//hit was not used before
		      
		      //set conformal mapping if looking for nonvertex tracks:
		      if(!fVertexConstraint)
			{
			  hit->SetAllCoord(start_hit);
			}
		     
		      if(track)//track search - look for nearest neighbor to extrapolated track
			{
			  if(!VerifyRange(start_hit,hit))
			    continue;
			  			  
			  test_hit = EvaluateHit(start_hit,hit,track);
			  
			  if(test_hit == 0)//chi2 not good enough, keep looking
			    continue;
			  else if(test_hit==2)//chi2 good enough, return it
			    return hit;
			  else
			    closest_hit = hit;//chi2 acceptable, but keep looking
			  
			}//track search
		      
		      else //tracklet search, look for nearest neighbor
			{
			  
			  if((dist=CalcDistance(start_hit,hit)) < closest_dist)
			    {
			      if(!VerifyRange(start_hit,hit))
				continue;
			      closest_dist = dist;
			      closest_hit = hit;
			 
			      //if this hit is good enough, return it:
			      if(closest_dist < fGoodDist)
			        return closest_hit;
			    }
			  else
			    continue;//sub hit was farther away than a hit before
			  
			}//tracklet search
		      
		    }//hit not used before
		  
		  else continue; //sub hit was used before
		  
		}//loop over hits in sub segment
	     	      
	    }//loop over sub segments
	  	  
	}//loop over subrows
      
    }//else

  //closest hit found:
  if(closest_hit)// && closest_dist < mMaxDist)
    return closest_hit;
  else
    return 0;
}

Int_t AliHLTTPCConfMapper::EvaluateHit(AliHLTTPCConfMapPoint *start_hit,AliHLTTPCConfMapPoint *hit,AliHLTTPCConfMapTrack *track) 
{
  //Check if space point gives a fit with acceptable chi2.
  
  Double_t temp,dxy,lchi2,dx,dy,slocal,dsz,lszChi2;
  temp = (track->a2Xy*hit->GetXprime()-hit->GetYprime()+track->a1Xy);
  dxy = temp*temp/(track->a2Xy*track->a2Xy + 1.);
  
  //Calculate chi2
  lchi2 = (dxy*hit->GetXYWeight());
  
  if(lchi2 > track->fChiSq[0])//chi2 was worse than before.
    return 0;
    
  //calculate s and the distance hit-line
  dx = start_hit->GetX()-hit->GetX();
  dy = start_hit->GetY()-hit->GetY();
  //slocal = track->fLength+sqrt(dx*dx+dy*dy);
  slocal = track->GetLength()+sqrt(dx*dx+dy*dy);
  
  temp = (track->a2Sz*slocal-hit->GetZ()+track->a1Sz);
  dsz = temp*temp/(track->a2Sz*track->a2Sz+1);
  
  //calculate chi2
  lszChi2 = dsz*hit->GetZWeight();
  lchi2 += lszChi2;
  
    
  //check whether chi2 is better than previous one:
  if(lchi2 < track->fChiSq[0])
    {
      track->fChiSq[0] = lchi2;
      track->fChiSq[1] = lszChi2;
    
      hit->SetS(slocal);
  
      //if chi2 good enough, stop here:
      if(lchi2 < fGoodHitChi2[fVertexConstraint]) 
        return 2;
      
      return 1;
    }
  
  return 0;
  
}

Double_t AliHLTTPCConfMapper::CalcDistance(const AliHLTTPCConfMapPoint *hit1,const AliHLTTPCConfMapPoint *hit2) const
{
  //Return distance between two clusters, defined by Pablo
  
  Double_t phi_diff = fabs( hit1->GetPhi() - hit2->GetPhi() );
  if (phi_diff > pi) phi_diff = twopi - phi_diff;
  
  return todeg*fabs((Float_t)((hit1->GetPadRow() - hit2->GetPadRow()) * (phi_diff + fabs( hit1->GetEta() - hit2->GetEta()))));
}

Bool_t AliHLTTPCConfMapper::VerifyRange(const AliHLTTPCConfMapPoint *hit1,const AliHLTTPCConfMapPoint *hit2) const
{
  //Check if the hit are within reasonable range in phi and eta
  Double_t dphi,deta;//maxphi=0.1,maxeta=0.1;
  dphi = fabs(hit1->GetPhi() - hit2->GetPhi());
  if(dphi > pi) dphi = fabs(twopi - dphi);
  if(dphi > fMaxPhi) return false;
  
  deta = fabs(hit1->GetEta() - hit2->GetEta());
  if(deta > fMaxEta) return false;

  return true;

}

Double_t AliHLTTPCConfMapper::TrackletAngle(AliHLTTPCConfMapTrack *track,Int_t n) const
{
  // Returns the angle 'between' the last three points (started at point number n) on this track.
  
  if(n > track->GetNumberOfPoints())
    n = track->GetNumberOfPoints();
  
  if(n<3)
    return 0;
  
  Double_t x1[2];
  Double_t x2[2];
  Double_t x3[2];
  Double_t angle1,angle2;
  Int_t counter=0;
  for(track->StartLoop(); track->LoopDone(); track->GetNextHit())
    {
      AliHLTTPCConfMapPoint *p = (AliHLTTPCConfMapPoint*)track->currentHit;
      if( (n-1) == counter)
	{
	  x1[0] = p->GetX();
	  x1[1] = p->GetY();
	}
      else if( (n-2) == counter)
	{
	  x2[0] = p->GetX();
	  x2[1] = p->GetY();
	}
      else if( (n-3) == counter)
	{
	  x3[0] = p->GetX();
	  x3[1] = p->GetY();
	}
      counter++;
    }
  
  angle1 = atan2(x2[1]-x3[1],x2[0]-x3[0]);
  angle2 = atan2(x1[1]-x2[1],x1[0]-x2[0]);
  
  return fabs(angle1-angle2);
  
  /*
    Double_t x1[2];
  Double_t x2[2];
  Double_t angle1,angle2;
  TObjArray *hits = track->GetHits();
  
  if (n > track->GetNumberOfPoints()) {
    n = track->GetNumberOfPoints();
  }

  if (n<3) 
    return 0;
  

  x1[0] = ((AliHLTTPCConfMapPoint *)hits->At(n-2))->GetX() - ((AliHLTTPCConfMapPoint *)hits->At(n-3))->GetX();
  x1[1] = ((AliHLTTPCConfMapPoint *)hits->At(n-2))->GetY() - ((AliHLTTPCConfMapPoint *)hits->At(n-3))->GetY();

  x2[0] = ((AliHLTTPCConfMapPoint *)hits->At(n-1))->GetX() - ((AliHLTTPCConfMapPoint *)hits->At(n-2))->GetX();
  x2[1] = ((AliHLTTPCConfMapPoint *)hits->At(n-1))->GetY() - ((AliHLTTPCConfMapPoint *)hits->At(n-2))->GetY();
  
  angle1 = atan2(x1[1],x1[0]);
  angle2 = atan2(x2[1],x1[0]);
  return fabs(angle1-angle2);
  */
}

Int_t AliHLTTPCConfMapper::FillTracks()
{
  //Fill track parameters. Which basically means do a fit of helix in real space,
  //which should be done in order to get nice tracks.
  
  Int_t num_of_tracks = nTracks;
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCConfMapper::FillTracks","nTracks")<<AliHLTTPCLog::kDec<<
    "Number of found tracks: "<<nTracks<<ENDLOG;
  
  if(nTracks == 0)
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapper::FillTracks","nTracks")<<AliHLTTPCLog::kDec<<
	"No tracks found!!"<<ENDLOG;
      return 0;
    }
  
  //  fTrack->Sort();
  for(int i=0; i<num_of_tracks; i++)
    {
      AliHLTTPCConfMapTrack *track = (AliHLTTPCConfMapTrack*)fTrack->GetTrack(i);
      track->Fill(fVertex,fMaxDca);
      
    }
  return 1;

}

Double_t AliHLTTPCConfMapper::CpuTime()
{
  //Return the Cputime in seconds.
 struct timeval tv;
 gettimeofday( &tv, NULL );
 return tv.tv_sec+(((Double_t)tv.tv_usec)/1000000.);
 //return (Double_t)(clock()) / CLOCKS_PER_SEC;
}
