
//Author:        Anders Strand Vestbo
//Last Modified: 15.12.2000

#include <iostream.h>
#include <time.h>
#include <math.h>
#include "AliL3ConfMapper.h"

#include "AliL3Logging.h"
#include "AliL3Transform.h"
#include "AliL3Vertex.h"
#include "AliL3ConfMapTrack.h"
#include "AliL3ConfMapPoint.h"
#include "AliL3TrackArray.h"

//
//AliL3ConfMapper
//
//Conformal mapping base class

ClassImp(AliL3ConfMapper)

Double_t AliL3ConfMapper::pi=3.14159265358979323846;
Double_t AliL3ConfMapper::twopi=2*pi;
Double_t AliL3ConfMapper::todeg=180./pi;

AliL3ConfMapper::AliL3ConfMapper()
{
  //Default constructor
  fVertex = NULL;
  fTrack = NULL;
  fHit = NULL;
  fVolume = NULL;
  fRow = NULL;
  fBench = (Bool_t)true;
  fParamSet = (Bool_t)false;
  fVertexConstraint = (Bool_t)true;
}


AliL3ConfMapper::~AliL3ConfMapper()
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

void AliL3ConfMapper::InitSector(Int_t sector,Int_t *rowrange,Float_t *etarange)
{
  //Initialize tracker for tracking in a given sector.
  //Resets track and hit arrays.
  //Here it is also possible to specify a subsector, by defining
  //rowrange[0]=innermost row;
  //rowrange[1]=outermostrow;
  //Finally you can specify etaslices to save time (assuming a good seed from TRD...)
  
  if(fHit)
    {
      delete [] fHit;
    }
  
  if(fTrack) 
    {
      delete fTrack;
    }

  //Define tracking area:
  if(rowrange)
    {
      fRowMin = rowrange[0];
      fRowMax = rowrange[1];
    }
  else //complete sector
    {
      fRowMin = 0;
      fRowMax = 173;
    }
  if(etarange)
    {
      fEtaMin = etarange[0];
      fEtaMax = etarange[1];
    }
  else
    {
      fEtaMin = 0;
      fEtaMax = sector < 18 ? 0.9 : -0.9;
    }
  
  //Set the angles to sector 2:
  fPhiMin = -1.*10/todeg;//fParam->GetAngle(sector) - 10/todeg;
  fPhiMax = 10/todeg;//fParam->GetAngle(sector) + 10/todeg;

  //rotation angles for sector 2:
  //cos: 0.766044 sin: 0.642788
  
  Int_t max_num_of_tracks = 3000;
  Int_t max_num_of_hits = 90000;

  fHit = new AliL3ConfMapPoint[max_num_of_hits];
  fTrack = new AliL3TrackArray("AliL3ConfMapTrack",max_num_of_tracks);
  
  nTracks=0;
  fClustersUnused = 0;
  
  fNumRowSegment = fRowMax - fRowMin; //number of rows to be considered by tracker
 
}



Bool_t AliL3ConfMapper::ReadHits(UInt_t count, AliL3SpacePointData* hits )
{
  Int_t nhit=(Int_t)count; 
  for (Int_t i=0;i<nhit;i++)
    fHit[i].ReadHits(&(hits[i]));
  fClustersUnused += nhit;
  LOG(AliL3Log::kInformational,"AliL3ConfMapper::ReadHits","#hits")<<AliL3Log::kDec
  <<"hit_counter: "<<nhit<<" count: "<<count<<ENDLOG;

  return true;
}


void AliL3ConfMapper::SetPointers()
{
  //Data organization.
  //Allocate volumes, set conformal coordinates and pointers.
  fNumRowSegmentPlusOne = 174;//fNumRowSegment+1;
  fNumPhiSegmentPlusOne = fNumPhiSegment+1;
  fNumEtaSegmentPlusOne = fNumEtaSegment+1;
  fNumPhiEtaSegmentPlusOne = fNumPhiSegmentPlusOne*fNumEtaSegmentPlusOne;
  fBounds = fNumRowSegmentPlusOne * fNumPhiSegmentPlusOne * fNumEtaSegmentPlusOne;
  
  //Allocate volumes:
  if(fVolume) delete [] fVolume;
  if(fRow) delete [] fRow;
  
  LOG(AliL3Log::kInformational,"AliL3ConfMapper::SetPointers","Memory")<<AliL3Log::kDec<<
    "Allocating "<<fBounds*sizeof(AliL3ConfMapContainer)<<" Bytes to fVolume"<<ENDLOG;
  LOG(AliL3Log::kInformational,"AliL3ConfMapper::SetPointers","Memory")<<AliL3Log::kDec<<
    "Allocating "<<fNumRowSegmentPlusOne*sizeof(AliL3ConfMapContainer)<<" Bytes to fRow"<<ENDLOG;
  
  fVolume = new AliL3ConfMapContainer[fBounds];
  fRow = new AliL3ConfMapContainer[fNumRowSegmentPlusOne];

  //set volumes to zero:
  memset(fVolume,0,fBounds*sizeof(AliL3ConfMapContainer));
  memset(fRow,0,fNumRowSegmentPlusOne*sizeof(AliL3ConfMapContainer));
  
  Float_t phiSlice = (fPhiMax-fPhiMin)/fNumPhiSegment;
  Float_t etaSlice = (fEtaMax-fEtaMin)/fNumEtaSegment;
  
  Int_t volumeIndex;
  for(Int_t j=0; j<fClustersUnused; j++)
    {
      
      //AliL3ConfMapPoint *thisHit = (AliL3ConfMapPoint*)fHit->At(j);
      AliL3ConfMapPoint *thisHit = &(fHit[j]);

      thisHit->Setup(fVertex);
      
      Int_t localrow = thisHit->GetPadRow();
      
      //reset pointers:
      thisHit->nextVolumeHit=thisHit->nextRowHit=0;
      
      if(localrow < fRowMin || localrow > fRowMax)
	continue;
      

      //Get indexes:
      thisHit->phiIndex=(Int_t)((thisHit->GetPhi()-fPhiMin)/phiSlice +1);
      
      if(thisHit->phiIndex<1 || thisHit->phiIndex>fNumPhiSegment)
	{
	  fPhiHitsOutOfRange++;
	  continue;
	}
      
      thisHit->etaIndex=(Int_t)((thisHit->GetEta()-fEtaMin)/etaSlice + 1);
      if(thisHit->etaIndex<1 || thisHit->etaIndex>fNumEtaSegment)
	{
	  fEtaHitsOutOfRange++;
	  continue;
	}
                  
      //set volume pointers
      volumeIndex = localrow*fNumPhiEtaSegmentPlusOne+thisHit->phiIndex*fNumEtaSegmentPlusOne+thisHit->etaIndex;
      if(fVolume[volumeIndex].first == NULL)
	fVolume[volumeIndex].first = (void *)thisHit;
      else
	((AliL3ConfMapPoint *)fVolume[volumeIndex].last)->nextVolumeHit=thisHit;
      fVolume[volumeIndex].last = (void *)thisHit;
      
      
      //set row pointers
      if(fRow[localrow].first == NULL)
	fRow[localrow].first = (void *)thisHit;
      else
	((AliL3ConfMapPoint *)(fRow[localrow].last))->nextRowHit = thisHit;
      fRow[localrow].last = (void *)thisHit;
      
    }
  
}

void AliL3ConfMapper::MainVertexTracking_a()
{
  //Tracking with vertex constraint.

  if(!fParamSet)
    {
      LOG(AliL3Log::kError,"AliL3ConfMapper::MainVertexTracking","Parameters")<<AliL3Log::kDec<<
	"Tracking parameters not set!"<<ENDLOG;
      return;
    }

  
  SetPointers();
  SetVertexConstraint(true);
}

void AliL3ConfMapper::MainVertexTracking_b()
{
  //Tracking with vertex constraint.

  if(!fParamSet)
    {
      LOG(AliL3Log::kError,"AliL3ConfMapper::MainVertexTracking","Parameters")<<AliL3Log::kDec<<
	"Tracking parameters not set!"<<ENDLOG;
      return;
    }

  ClusterLoop();
}

void AliL3ConfMapper::MainVertexTracking()
{
  //Tracking with vertex constraint.

  if(!fParamSet)
    {
      LOG(AliL3Log::kError,"AliL3ConfMapper::MainVertexTracking","Parameters")<<AliL3Log::kDec<<
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
    LOG(AliL3Log::kInformational,"AliL3ConfMapper::MainVertexTracking","Timing")<<AliL3Log::kDec<<
      "Tracking finished in "<<cpuTime*1000<<" ms"<<ENDLOG;
    
  return;
}

void AliL3ConfMapper::NonVertexTracking()
{
  //Tracking with no vertex constraint. This should be called after doing MainVertexTracking,
  //in order to do tracking on the remaining clusters.
  //The conformal mapping is now done with respect to the first cluster
  //assosciated with this track.

  SetVertexConstraint(false);
  ClusterLoop();
  LOG(AliL3Log::kInformational,"AliL3ConfMapper::NonVertexTracking","ntracks")<<AliL3Log::kDec<<
    "Number of nonvertex tracks found: "<<nTracks-fMainVertexTracks<<ENDLOG;

  return;
}

void AliL3ConfMapper::MainVertexSettings(Int_t trackletlength, Int_t tracklength,
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
  
}

void AliL3ConfMapper::NonVertexSettings(Int_t trackletlength, Int_t tracklength,
					Int_t rowscopetracklet, Int_t rowscopetrack)
{
  SetTrackletLength(trackletlength,(Bool_t)false);
  SetRowScopeTracklet(rowscopetracklet, (Bool_t)false);
  SetRowScopeTrack(rowscopetrack, (Bool_t)false);
  SetMinPoints(tracklength,(Bool_t)false);
}

void AliL3ConfMapper::SetTrackCuts(Double_t hitChi2Cut, Double_t goodHitChi2, Int_t trackChi2Cut,Int_t maxdist)
{
  //Settings for tracks. The cuts are:
  //HitChi2Cut:     Maximum hit chi2
  //goodHitChi2:    Chi2 to stop look for next hit
  //trackChi2Cut:   Maximum track chi2
  //maxdist:        Maximum distance between two clusters when forming segments

  fHitChi2Cut = hitChi2Cut;
  fGoodHitChi2 = goodHitChi2;
  fTrackChi2Cut = trackChi2Cut;
  fMaxDist = maxdist;
}

void AliL3ConfMapper::SetTrackletCuts(Double_t maxangle,Double_t goodDist, Bool_t vertex_constraint)
{
  //Sets cuts of tracklets. Right now this is only:
  //maxangle:  Maximum angle when forming segments (if trackletlength > 2)
 
  fGoodDist=goodDist;
  SetMaxAngleTracklet(maxangle, vertex_constraint);
}

void AliL3ConfMapper::ClusterLoop()
{
  //Loop over hits, starting at outermost padrow, and trying to build segments.

  Int_t row_segm,lastrow = fRowMin + fMinPoints[fVertexConstraint];
  AliL3ConfMapPoint *hit;
  
  //Loop over rows, and try to create tracks from the hits.
  //Starts at the outermost row, and loops as long as a track can be build, due to length.
  
  for(row_segm = fRowMax; row_segm >= lastrow; row_segm--)
    {
      if(fRow[row_segm].first && ((AliL3ConfMapPoint*)fRow[row_segm].first)->GetPadRow() < fRowMin + 1)
	break;
      for(hit = (AliL3ConfMapPoint*)fRow[row_segm].first; hit!=0; hit=hit->nextRowHit)
	{
	  if(hit->GetUsage() == true)
	    continue;
	  else
	    CreateTrack(hit);
	}
    }
  
  return;
}


void AliL3ConfMapper::CreateTrack(AliL3ConfMapPoint *hit)
{
  //Tries to create a track from the initial hit given by ClusterLoop()

  AliL3ConfMapPoint *closest_hit = NULL;
  AliL3ConfMapTrack *track = NULL;
  
  Int_t point;
  Int_t tracks = nTracks;
  nTracks++;

  track = (AliL3ConfMapTrack*)fTrack->NextTrack();

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
	  dx = ((AliL3ConfMapPoint*)closest_hit)->GetX() - ((AliL3ConfMapPoint*)hit)->GetX();
	  dy = ((AliL3ConfMapPoint*)closest_hit)->GetY() - ((AliL3ConfMapPoint*)hit)->GetY();
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
	      track->fChiSq[0] = fHitChi2Cut;
	      closest_hit = GetNextNeighbor((AliL3ConfMapPoint*)track->lastHit,track);
	      
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
	  if(track->GetNumberOfPoints() < fMinPoints[fVertexConstraint] || normalized_chi2 > fTrackChi2Cut)
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

AliL3ConfMapPoint *AliL3ConfMapper::GetNextNeighbor(AliL3ConfMapPoint *start_hit,
					  AliL3ConfMapTrack *track)
{
  //When forming segments: Finds closest hit to input hit
  //When forming tracks: Find closest hit to track fit.
  
  Double_t dist,closest_dist = fMaxDist;
  
  AliL3ConfMapPoint *hit = NULL;
  AliL3ConfMapPoint *closest_hit = NULL;
    
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
  Int_t loop_eta[9] = {0,0,0,-1,-1,-1,1,1,1};
  Int_t loop_phi[9] = {0,-1,1,0,-1,1,0,-1,1};
  
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
	      
	      if(sub_phi_segm<0)
		sub_phi_segm += fNumPhiSegment;
	      
	      else if(sub_phi_segm >=fNumPhiSegment)
		sub_phi_segm -= fNumPhiSegment;
	      
	      //loop over sub eta segments
	      
	      sub_eta_segm = start_hit->etaIndex + loop_eta[i];
	      
	      if(sub_eta_segm < 0 || sub_eta_segm >=fNumEtaSegment)
		continue;//segment exceeds bounds->skip it
	      
	      //loop over hits in this sub segment:
	      volumeIndex= sub_row_segm*fNumPhiEtaSegmentPlusOne +
		sub_phi_segm*fNumEtaSegmentPlusOne + sub_eta_segm;
	      
	      if(volumeIndex<0)
		{//debugging
		  LOG(AliL3Log::kError,"AliL3ConfMapper::GetNextNeighbor","Memory")<<AliL3Log::kDec<<
		    "VolumeIndex error "<<volumeIndex<<ENDLOG;
		}
	      
	      for(hit = (AliL3ConfMapPoint*)fVolume[volumeIndex].first;
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

Int_t AliL3ConfMapper::EvaluateHit(AliL3ConfMapPoint *start_hit,AliL3ConfMapPoint *hit,AliL3ConfMapTrack *track) 
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
      if(lchi2 < fGoodHitChi2) 
        return 2;
      
      return 1;
    }
  
  return 0;
  
}

Double_t AliL3ConfMapper::CalcDistance(const AliL3ConfMapPoint *hit1,const AliL3ConfMapPoint *hit2) const
{
  //Return distance between two clusters, defined by Pablo
  
  Double_t phi_diff = fabs( hit1->GetPhi() - hit2->GetPhi() );
  if (phi_diff > pi) phi_diff = twopi - phi_diff;
  
  return todeg*fabs(hit1->GetPadRow() - hit2->GetPadRow()) * (phi_diff + fabs( hit1->GetEta() - hit2->GetEta() ));
}

Bool_t AliL3ConfMapper::VerifyRange(const AliL3ConfMapPoint *hit1,const AliL3ConfMapPoint *hit2) const
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

Double_t AliL3ConfMapper::TrackletAngle(const AliL3ConfMapTrack *track,Int_t n) const
{
  // Returns the angle 'between' the last three points (started at point number n) on this track.

  return 0;
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
  

  x1[0] = ((AliL3ConfMapPoint *)hits->At(n-2))->GetX() - ((AliL3ConfMapPoint *)hits->At(n-3))->GetX();
  x1[1] = ((AliL3ConfMapPoint *)hits->At(n-2))->GetY() - ((AliL3ConfMapPoint *)hits->At(n-3))->GetY();

  x2[0] = ((AliL3ConfMapPoint *)hits->At(n-1))->GetX() - ((AliL3ConfMapPoint *)hits->At(n-2))->GetX();
  x2[1] = ((AliL3ConfMapPoint *)hits->At(n-1))->GetY() - ((AliL3ConfMapPoint *)hits->At(n-2))->GetY();
  
  angle1 = atan2(x1[1],x1[0]);
  angle2 = atan2(x2[1],x1[0]);
  return fabs(angle1-angle2);
  */
}

Int_t AliL3ConfMapper::FillTracks()
{
  //Fill track parameters. Which basically means do a fit of helix in real space,
  //which should be done in order to get nice tracks.
  
  Int_t num_of_tracks = nTracks;
  LOG(AliL3Log::kInformational,"AliL3ConfMapper::FillTracks","nTracks")<<AliL3Log::kDec<<
    "Number of found tracks: "<<nTracks<<ENDLOG;
  
  if(nTracks == 0)
    {
      LOG(AliL3Log::kError,"AliL3ConfMapper::FillTracks","nTracks")<<AliL3Log::kDec<<
	"No tracks found!!"<<ENDLOG;
      return 0;
    }
  
  //  fTrack->Sort();
  for(int i=0; i<num_of_tracks; i++)
    {
      AliL3ConfMapTrack *track = (AliL3ConfMapTrack*)fTrack->GetTrack(i);
      track->Fill(fVertex,fMaxDca);
    }
  return 1;

}

Double_t AliL3ConfMapper::CpuTime()
{
  //Return the Cputime in seconds.

  return (Double_t)(clock()) / CLOCKS_PER_SEC;
}
