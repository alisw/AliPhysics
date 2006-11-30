// @(#) $Id$

/** \class AliHLTConfMapper
<pre>
//_____________________________________________________________
// AliHLTConfMapper
//
// Conformal mapping base class
//
// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
// Copyright &copy ALICE HLT Group
</pre>
*/

#include <sys/time.h>
 
#include "AliHLTStandardIncludes.h"
#include "AliHLTRootTypes.h"
#include "AliHLTSpacePointData.h"
#include "AliHLTLogging.h" 
#include "AliHLTVertex.h"
#include "AliHLTConfMapTrack.h"
#include "AliHLTConfMapPoint.h"
#include "AliHLTTrackArray.h"
#include "AliHLTTransform.h"
#include "AliHLTConfMapper.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTConfMapper)

AliHLTConfMapper::AliHLTConfMapper()
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


AliHLTConfMapper::~AliHLTConfMapper()
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
 
void AliHLTConfMapper::InitVolumes()
{
  //Data organization.
  //Allocate volumes, set conformal coordinates and pointers.
  
  //Should be done after setting the track parameters
  
  fNumRowSegmentPlusOne = AliHLTTransform::GetNRows();//NumRows[0]; //Maximum 32.
  fNumPhiSegmentPlusOne = fNumPhiSegment+1;
  fNumEtaSegmentPlusOne = fNumEtaSegment+1;
  fNumPhiEtaSegmentPlusOne = fNumPhiSegmentPlusOne*fNumEtaSegmentPlusOne;
  fBounds = fNumRowSegmentPlusOne * fNumPhiSegmentPlusOne * fNumEtaSegmentPlusOne;
  
  //Allocate volumes:
  if(fVolume) delete [] fVolume;
  if(fRow) delete [] fRow;
  
  LOG(AliHLTLog::kInformational,"AliHLTConfMapper::InitVolumes","Memory")<<AliHLTLog::kDec<<
    "Allocating "<<fBounds*sizeof(AliHLTConfMapContainer)<<" Bytes to fVolume"<<ENDLOG;
  LOG(AliHLTLog::kInformational,"AliHLTConfMapper::InitVolumes","Memory")<<AliHLTLog::kDec<<
    "Allocating "<<fNumRowSegmentPlusOne*sizeof(AliHLTConfMapContainer)<<" Bytes to fRow"<<ENDLOG;
  
  fVolume = new AliHLTConfMapContainer[fBounds];
  fRow = new AliHLTConfMapContainer[fNumRowSegmentPlusOne];
  
  memset(fVolume,0,fBounds*sizeof(AliHLTConfMapContainer));
  memset(fRow,0,fNumRowSegmentPlusOne*sizeof(AliHLTConfMapContainer));
  
  Int_t maxnumoftracks = 2000;
  Int_t maxnumofhits = 120000;
  
  if(fHit)
    delete [] fHit;
  if(fTrack)
    delete fTrack;
    
  fHit = new AliHLTConfMapPoint[maxnumofhits];
  fTrack = new AliHLTTrackArray("AliHLTConfMapTrack",maxnumoftracks);
}

void AliHLTConfMapper::InitSector(Int_t sector,Int_t *rowrange,Float_t *etarange)
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
      fRowMax = AliHLTTransform::GetNRows() - 1;
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
  fPhiMin = -10*AliHLTTransform::ToRad();//fParam->GetAngle(sector) - 10/todeg;
  fPhiMax =  10*AliHLTTransform::ToRad();//fParam->GetAngle(sector) + 10/todeg;

  fNTracks=0;
  fMainVertexTracks = 0;
  fClustersUnused = 0;
  fEtaHitsOutOfRange=0;
  fPhiHitsOutOfRange=0;
  
  fNumRowSegment = fRowMax - fRowMin; //number of rows to be considered by tracker
  LOG(AliHLTLog::kInformational,"AliHLTConfMapper::InitSector","B-field")
    <<"Tracker initializing with a magnetic field of "<<AliHLTTransform::GetBField()<<ENDLOG;
  
  fTrack->Reset();
}

Bool_t AliHLTConfMapper::ReadHits(UInt_t count, AliHLTSpacePointData* hits )
{
  //read hits
  Int_t nhit=(Int_t)count; 
  for (Int_t i=0;i<nhit;i++)
    {
      fHit[i].Reset();
      fHit[i].ReadHits(&(hits[i]));
    }
  fClustersUnused += nhit;
  LOG(AliHLTLog::kInformational,"AliHLTConfMapper::ReadHits","#hits")
    <<AliHLTLog::kDec<<"hit_counter: "<<nhit<<" count: "<<count<<ENDLOG;
  
  return true;
}

void AliHLTConfMapper::SetPointers()
{
  //Check if there are not enough clusters to make a track in this sector
  //Can happen in pp events.

  if(fClustersUnused < fMinPoints[fVertexConstraint])
    return;
  
  //Reset detector volumes
  memset(fVolume,0,fBounds*sizeof(AliHLTConfMapContainer));
  memset(fRow,0,fNumRowSegmentPlusOne*sizeof(AliHLTConfMapContainer));
  
  Float_t phiSlice = (fPhiMax-fPhiMin)/fNumPhiSegment;
  Float_t etaSlice = (fEtaMax-fEtaMin)/fNumEtaSegment;

  Int_t volumeIndex;
  Int_t localcounter=0;
  for(Int_t j=0; j<fClustersUnused; j++)
    {
      //AliHLTConfMapPoint *thisHit = (AliHLTConfMapPoint*)fHit->At(j);
      AliHLTConfMapPoint *thisHit = &(fHit[j]);

      thisHit->Setup(fVertex);
      
      Int_t localrow = thisHit->GetPadRow();
      
      if(localrow < fRowMin || localrow > fRowMax)
	continue;

      //Get indexes:
      thisHit->SetPhiIndex((Int_t)((thisHit->GetPhi()-fPhiMin)/phiSlice +1));
      
      if(thisHit->GetPhiIndex()<1 || thisHit->GetPhiIndex()>fNumPhiSegment)
	{
	  //cout << "Phiindex: " << thisHit->phiIndex << " " << thisHit->GetPhi() << endl;
	  fPhiHitsOutOfRange++;
	  continue;
	}
      
      thisHit->SetEtaIndex((Int_t)((thisHit->GetEta()-fEtaMin)/etaSlice + 1));
      if(thisHit->GetEtaIndex()<1 || thisHit->GetEtaIndex()>fNumEtaSegment)
	{
	  //cout << "Etaindex: " << thisHit->etaIndex << " " << thisHit->GetEta() << endl;
	  fEtaHitsOutOfRange++;
	  continue;
	}
      localcounter++;
      
      volumeIndex = (localrow-fRowMin)*fNumPhiEtaSegmentPlusOne + 
                    thisHit->GetPhiIndex()*fNumEtaSegmentPlusOne+thisHit->GetEtaIndex();
      
      if(fVolume[volumeIndex].first == NULL)
	fVolume[volumeIndex].first = (void *)thisHit;
      else
 	((AliHLTConfMapPoint *)fVolume[volumeIndex].last)->SetNextVolumeHit(thisHit);
      fVolume[volumeIndex].last = (void *)thisHit;
      
      
      //set row pointers
      if(fRow[(localrow-fRowMin)].first == NULL)
 	fRow[(localrow-fRowMin)].first = (void *)thisHit;
      else
 	((AliHLTConfMapPoint *)(fRow[(localrow-fRowMin)].last))->SetNextRowHit(thisHit);
	fRow[(localrow-fRowMin)].last = (void *)thisHit;
    }
  
  if(fClustersUnused>0 && localcounter==0)
    LOG(AliHLTLog::kError,"AliHLTConfMapper::SetPointers","Parameters")
      <<AliHLTLog::kDec<<"No points passed to track finder, hits out of range: "
      <<fEtaHitsOutOfRange+fPhiHitsOutOfRange<<ENDLOG;

  Int_t hits_accepted=fClustersUnused-(fEtaHitsOutOfRange+fPhiHitsOutOfRange);
  LOG(AliHLTLog::kInformational,"AliHLTConfMapper::SetPointers","Setup")
    <<"Setup finished, hits out of range: "<<fEtaHitsOutOfRange+fPhiHitsOutOfRange
    <<" hits accepted "<<hits_accepted<<ENDLOG;
}

void AliHLTConfMapper::MainVertexTracking_a()
{
  //Tracking with vertex constraint.

  if(!fParamSet[(Int_t)kTRUE])
    {
      LOG(AliHLTLog::kError,"AliHLTConfMapper::MainVertexTracking","Parameters")<<AliHLTLog::kDec<<
	"Tracking parameters not set!"<<ENDLOG;
      return;
    }

  Double_t initCpuTime,cpuTime;
  initCpuTime = CpuTime();

  SetPointers();
  SetVertexConstraint(true);
  cpuTime = CpuTime() - initCpuTime;
  if(fBench)
    LOG(AliHLTLog::kInformational,"AliHLTConfMapper::MainVertexTracking_a","Timing")
      <<AliHLTLog::kDec<<"Setup finished in "<<cpuTime*1000<<" ms"<<ENDLOG;
  
}

void AliHLTConfMapper::MainVertexTracking_b()
{
  //Tracking with vertex constraint.

  if(!fParamSet[(Int_t)kTRUE])
    {
      LOG(AliHLTLog::kError,"AliHLTConfMapper::MainVertexTracking","Parameters")<<AliHLTLog::kDec<<
	"Tracking parameters not set!"<<ENDLOG;
      return;
    }
  Double_t initCpuTime,cpuTime;
  initCpuTime = CpuTime();
  
  ClusterLoop();
 
  cpuTime = CpuTime() - initCpuTime;
  if(fBench)
    LOG(AliHLTLog::kInformational,"AliHLTConfMapper::MainVertexTracking_b","Timing")
      <<AliHLTLog::kDec<<"Main Tracking finished in "<<cpuTime*1000<<" ms"<<ENDLOG;
}

void AliHLTConfMapper::MainVertexTracking()
{
  //Tracking with vertex constraint.

  if(!fParamSet[(Int_t)kTRUE])
    {
      LOG(AliHLTLog::kError,"AliHLTConfMapper::MainVertexTracking","Parameters")<<AliHLTLog::kDec<<
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
    LOG(AliHLTLog::kInformational,"AliHLTConfMapper::MainVertexTracking","Timing")<<AliHLTLog::kDec<<
      "Tracking finished in "<<cpuTime*1000<<" ms"<<ENDLOG;
  
  return;
}

void AliHLTConfMapper::NonVertexTracking()
{
  //Tracking with no vertex constraint. This should be called after doing MainVertexTracking,
  //in order to do tracking on the remaining clusters.
  //The conformal mapping is now done with respect to the first cluster
  //assosciated with this track.
  
  if(!fParamSet[(Int_t)kFALSE])
    {
      LOG(AliHLTLog::kError,"AliHLTConfMapper::NonVertexTracking","Parameters")<<AliHLTLog::kDec<<
	"Tracking parameters not set!"<<ENDLOG;
      return;
    }
  
  SetVertexConstraint(false);
  ClusterLoop();
  LOG(AliHLTLog::kInformational,"AliHLTConfMapper::NonVertexTracking","ntracks")<<AliHLTLog::kDec<<
    "Number of nonvertex tracks found: "<<(fNTracks-fMainVertexTracks)<<ENDLOG;
  return;
}

void AliHLTConfMapper::MainVertexSettings(Int_t trackletlength, Int_t tracklength,
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

void AliHLTConfMapper::NonVertexSettings(Int_t trackletlength, Int_t tracklength,
					Int_t rowscopetracklet, Int_t rowscopetrack)
{
  //set parameters for non-vertex tracking
  SetTrackletLength(trackletlength,(Bool_t)false);
  SetRowScopeTracklet(rowscopetracklet, (Bool_t)false);
  SetRowScopeTrack(rowscopetrack, (Bool_t)false);
  SetMinPoints(tracklength,(Bool_t)false);
  SetParamDone(kFALSE);
}

void AliHLTConfMapper::SetTrackCuts(Double_t hitChi2Cut, Double_t goodHitChi2, Double_t trackChi2Cut,Int_t maxdist,Bool_t vertexconstraint)
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

void AliHLTConfMapper::SetTrackletCuts(Double_t maxangle,Double_t goodDist, Bool_t vc)
{
  //Sets cuts of tracklets. Right now this is only:
  //maxangle:  Maximum angle when forming segments (if trackletlength > 2)
 
  fGoodDist=goodDist;
  SetMaxAngleTracklet(maxangle, vc);
}

void AliHLTConfMapper::ClusterLoop()
{
  //Loop over hits, starting at outermost padrow, and trying to build segments.
  
  //Check if there are not enough clusters to make a track in this sector
  //Can happen in pp events.
  if(fClustersUnused < fMinPoints[fVertexConstraint])
    return;
  
  Int_t rowsegm,lastrow = fRowMin + fMinPoints[fVertexConstraint];
  AliHLTConfMapPoint *hit;
  
  //Loop over rows, and try to create tracks from the hits.
  //Starts at the outermost row, and loops as long as a track can be build, due to length.
  
  for(rowsegm = fRowMax; rowsegm >= lastrow; rowsegm--)
    {
      if(fRow[(rowsegm-fRowMin)].first && ((AliHLTConfMapPoint*)fRow[(rowsegm-fRowMin)].first)->GetPadRow() < fRowMin + 1)
	break;

      for(hit = (AliHLTConfMapPoint*)fRow[(rowsegm-fRowMin)].first; hit!=0; hit=hit->GetNextRowHit())
	{
	  if(hit->GetUsage() == true)
	    continue;
	  else
	    CreateTrack(hit);
	}
    }
  
  return;
}


void AliHLTConfMapper::CreateTrack(AliHLTConfMapPoint *hit)
{
  //Tries to create a track from the initial hit given by ClusterLoop()

  AliHLTConfMapPoint *closesthit = NULL;
  AliHLTConfMapTrack *track = NULL;
  
  Int_t point;
  Int_t tracks = fNTracks;
  fNTracks++;

  track = (AliHLTConfMapTrack*)fTrack->NextTrack();

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
      if((closesthit = GetNextNeighbor(hit)))
	{//closest hit exist
	  
	  //   Calculate track length in sz plane
	  dx = ((AliHLTConfMapPoint*)closesthit)->GetX() - ((AliHLTConfMapPoint*)hit)->GetX();
	  dy = ((AliHLTConfMapPoint*)closesthit)->GetY() - ((AliHLTConfMapPoint*)hit)->GetY();
	  //track->fLength += (Double_t)sqrt ( dx * dx + dy * dy ) ;
	  Double_t length = track->GetLength()+(Double_t)sqrt ( dx * dx + dy * dy );
	  track->SetLength(length);

	  //closesthit->SetS(track->fLength);
	  closesthit->SetS(track->GetLength());

	  //update fit parameters
	  track->UpdateParam(closesthit);
	  trackhitnumber[track->GetNumberOfPoints()-1] = closesthit->GetHitNumber();
	
	  hit = closesthit;
	}
      else
	{
	  //closest hit does not exist:
	  track->DeleteCandidate();
	  fTrack->RemoveLast();
	  fNTracks--;
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
	  fNTracks--;
	}
      
      else//good tracklet ->proceed, follow the trackfit
	{
	  tracks++;
	  			  
	  //define variables to keep the total chi:
	  Double_t xyChi2 = track->GetChiSq1();
	  Double_t szChi2 = track->GetChiSq2();
	  
	  for(point = fTrackletLength[fVertexConstraint]; point <= fNumRowSegment; point++)
	    {
	      track->SetChiSq1(fHitChi2Cut[fVertexConstraint]);
	      closesthit = GetNextNeighbor((AliHLTConfMapPoint*)track->GetLastHit(),track);
	      
	      if(closesthit)
		{
		  //keep total chi:
		  Double_t lxyChi2 = track->GetChiSq1()-track->GetChiSq2();
		  xyChi2 += lxyChi2;
		  closesthit->SetXYChi2(lxyChi2);
		  		  
		  //update track length:
		  track->SetLength(closesthit->GetS());
		  szChi2 += track->GetChiSq2();
		  closesthit->SetSZChi2(track->GetChiSq2());
		  
		  track->UpdateParam(closesthit);
		  trackhitnumber[track->GetNumberOfPoints()-1] = closesthit->GetHitNumber();
		  
		  //add closest hit to track
		  closesthit->SetUsage(true);
		  closesthit->SetTrackNumber(tracks-1);
		  
		}//closesthit
	      
	      else
		{
		  //closest hit does not exist
		  point = fNumRowSegment; //continue with next hit in segment
		}//else
	      
	    }//create tracks
	  
	  //store track chi2:
	  track->SetChiSq1(xyChi2);
	  track->SetChiSq2(szChi2);
	  Double_t normalizedchi2 = (track->GetChiSq1()+track->GetChiSq2())/track->GetNumberOfPoints();
	  
	  //remove tracks with not enough points already now
	  if(track->GetNumberOfPoints() < fMinPoints[fVertexConstraint] || normalizedchi2 > fTrackChi2Cut[fVertexConstraint])
	    {
	      track->SetProperties(false);
	      fNTracks--;
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

AliHLTConfMapPoint *AliHLTConfMapper::GetNextNeighbor(AliHLTConfMapPoint *starthit,
					  AliHLTConfMapTrack *track)
{
  //When forming segments: Finds closest hit to input hit
  //When forming tracks: Find closest hit to track fit.
  
  Double_t dist,closestdist = fMaxDist[fVertexConstraint];
  
  AliHLTConfMapPoint *hit = NULL;
  AliHLTConfMapPoint *closesthit = NULL;
    
  Int_t subrowsegm;
  Int_t subphisegm;
  Int_t subetasegm;
  Int_t volumeIndex;
  Int_t testhit;

  Int_t maxrow = starthit->GetPadRow()-1;
  Int_t minrow;

  if(track) //finding hit close to trackfit
    {
      minrow = starthit->GetPadRow()-fRowScopeTrack[fVertexConstraint];
    }
  else
    {
      minrow = starthit->GetPadRow()-fRowScopeTracklet[fVertexConstraint];
    }

  //make a smart loop
  Int_t loopeta[25] = {0,0,0,-1,-1,-1,1,1,1, 0,0,-1,-1,1,1,-2,-2,-2,-2,-2,2,2,2,2,2};
  Int_t loopphi[25] = {0,-1,1,0,-1,1,0,-1,1, -2,2,-2,2,-2,2,-2,-1,0,1,2,-2,-1,0,1,2};
  
  if(minrow < fRowMin)
    minrow = fRowMin;
  if(maxrow < fRowMin)
    return 0;  //reached the last padrow under consideration

  else
    {
      //loop over sub rows
      for(subrowsegm=maxrow; subrowsegm>=minrow; subrowsegm--)
	{
	  //loop over subsegments, in the order defined above.
	  for(Int_t i=0; i<9; i++)  
	    {
	      subphisegm = starthit->GetPhiIndex() + loopphi[i];
	      
	      if(subphisegm < 0 || subphisegm >= fNumPhiSegment)
		continue;
	      /*
		if(subphisegm<0)
		subphisegm += fNumPhiSegment;
		
		else if(subphisegm >=fNumPhiSegment)
		subphisegm -= fNumPhiSegment;
	      */
	      //loop over sub eta segments
	      
	      subetasegm = starthit->GetEtaIndex() + loopeta[i];
	      
	      if(subetasegm < 0 || subetasegm >=fNumEtaSegment)
		continue;//segment exceeds bounds->skip it
	      
	      //loop over hits in this sub segment:
	      volumeIndex=(subrowsegm-fRowMin)*fNumPhiEtaSegmentPlusOne +
		subphisegm*fNumEtaSegmentPlusOne + subetasegm;
	      
	      if(volumeIndex<0)
		{//debugging
		  LOG(AliHLTLog::kError,"AliHLTConfMapper::GetNextNeighbor","Memory")<<AliHLTLog::kDec<<
		    "VolumeIndex error "<<volumeIndex<<ENDLOG;
		}
	      
	      for(hit = (AliHLTConfMapPoint*)fVolume[volumeIndex].first;
		  hit!=0; hit = hit->GetNextVolumeHit())
		{
		  
		  if(!hit->GetUsage())
		    {//hit was not used before
		      
		      //set conformal mapping if looking for nonvertex tracks:
		      if(!fVertexConstraint)
			{
			  hit->SetAllCoord(starthit);
			}
		     
		      if(track)//track search - look for nearest neighbor to extrapolated track
			{
			  if(!VerifyRange(starthit,hit))
			    continue;
			  			  
			  testhit = EvaluateHit(starthit,hit,track);
			  
			  if(testhit == 0)//chi2 not good enough, keep looking
			    continue;
			  else if(testhit==2)//chi2 good enough, return it
			    return hit;
			  else
			    closesthit = hit;//chi2 acceptable, but keep looking
			  
			}//track search
		      
		      else //tracklet search, look for nearest neighbor
			{
			  
			  if((dist=CalcDistance(starthit,hit)) < closestdist)
			    {
			      if(!VerifyRange(starthit,hit))
				continue;
			      closestdist = dist;
			      closesthit = hit;
			 
			      //if this hit is good enough, return it:
			      if(closestdist < fGoodDist)
			        return closesthit;
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
  if(closesthit)// && closestdist < mMaxDist)
    return closesthit;
  else
    return 0;
}

Int_t AliHLTConfMapper::EvaluateHit(AliHLTConfMapPoint *starthit,AliHLTConfMapPoint *hit,AliHLTConfMapTrack *track) 
{
  //Check if space point gives a fit with acceptable chi2.
  
  Double_t temp,dxy,lchi2,dx,dy,slocal,dsz,lszChi2;
  temp = (track->GetA2Xy()*hit->GetXprime()-hit->GetYprime()+track->GetA1Xy());
  dxy = temp*temp/(track->GetA2Xy()*track->GetA2Xy() + 1.);
  
  //Calculate chi2
  lchi2 = (dxy*hit->GetXYWeight());
  
  if(lchi2 > track->GetChiSq1())//chi2 was worse than before.
    return 0;
    
  //calculate s and the distance hit-line
  dx = starthit->GetX()-hit->GetX();
  dy = starthit->GetY()-hit->GetY();
  //slocal = track->fLength+sqrt(dx*dx+dy*dy);
  slocal = track->GetLength()+sqrt(dx*dx+dy*dy);
  
  temp = (track->GetA2Sz()*slocal-hit->GetZ()+track->GetA1Sz());
  dsz = temp*temp/(track->GetA2Sz()*track->GetA2Sz()+1);
  
  //calculate chi2
  lszChi2 = dsz*hit->GetZWeight();
  lchi2 += lszChi2;
  
    
  //check whether chi2 is better than previous one:
  if(lchi2 < track->GetChiSq1())
    {
      track->SetChiSq1(lchi2);
      track->SetChiSq2(lszChi2);
    
      hit->SetS(slocal);
  
      //if chi2 good enough, stop here:
      if(lchi2 < fGoodHitChi2[fVertexConstraint]) 
        return 2;
      
      return 1;
    }
  
  return 0;
  
}

Double_t AliHLTConfMapper::CalcDistance(const AliHLTConfMapPoint *hit1,const AliHLTConfMapPoint *hit2) const
{
  //Return distance between two clusters, defined by Pablo
  
  Double_t phidiff = fabs( hit1->GetPhi() - hit2->GetPhi() );
  if (phidiff > AliHLTTransform::Pi()) phidiff = AliHLTTransform::TwoPi() - phidiff;
  
  return AliHLTTransform::ToDeg()*fabs((Float_t)((hit1->GetPadRow() - hit2->GetPadRow()) * 
         (phidiff + fabs( hit1->GetEta() - hit2->GetEta()))));
}

Bool_t AliHLTConfMapper::VerifyRange(const AliHLTConfMapPoint *hit1,const AliHLTConfMapPoint *hit2) const
{
  //Check if the hit are within reasonable range in phi and eta
  Double_t dphi,deta;//maxphi=0.1,maxeta=0.1;
  dphi = fabs(hit1->GetPhi() - hit2->GetPhi());
  if(dphi > AliHLTTransform::Pi()) dphi = fabs(AliHLTTransform::TwoPi() - dphi);
  if(dphi > fMaxPhi) return false;
  
  deta = fabs(hit1->GetEta() - hit2->GetEta());
  if(deta > fMaxEta) return false;

  return true;

}

Double_t AliHLTConfMapper::TrackletAngle(AliHLTConfMapTrack *track,Int_t n) const
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
      AliHLTConfMapPoint *p = (AliHLTConfMapPoint*)track->GetCurrentHit();
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
  

  x1[0] = ((AliHLTConfMapPoint *)hits->At(n-2))->GetX() - ((AliHLTConfMapPoint *)hits->At(n-3))->GetX();
  x1[1] = ((AliHLTConfMapPoint *)hits->At(n-2))->GetY() - ((AliHLTConfMapPoint *)hits->At(n-3))->GetY();

  x2[0] = ((AliHLTConfMapPoint *)hits->At(n-1))->GetX() - ((AliHLTConfMapPoint *)hits->At(n-2))->GetX();
  x2[1] = ((AliHLTConfMapPoint *)hits->At(n-1))->GetY() - ((AliHLTConfMapPoint *)hits->At(n-2))->GetY();
  
  angle1 = atan2(x1[1],x1[0]);
  angle2 = atan2(x2[1],x1[0]);
  return fabs(angle1-angle2);
  */
}

Int_t AliHLTConfMapper::FillTracks()
{
  //Fill track parameters. Which basically means do a fit of helix in real space,
  //which should be done in order to get nice tracks.
  
  Int_t numoftracks = fNTracks;
  if(fNTracks == 0)
    {
      LOG(AliHLTLog::kError,"AliHLTConfMapper::FillTracks","fNTracks")<<AliHLTLog::kDec<<
	"No tracks found!!"<<ENDLOG;
      return 0;
    }

  LOG(AliHLTLog::kInformational,"AliHLTConfMapper::FillTracks","fNTracks")<<AliHLTLog::kDec<<
    "Number of found tracks: "<<fNTracks<<ENDLOG;
  
  //  fTrack->Sort();
  for(Int_t i=0; i<numoftracks; i++)
    {
      AliHLTConfMapTrack *track = (AliHLTConfMapTrack*)fTrack->GetTrack(i);
      track->Fill(fVertex,fMaxDca);
    }
  return 1;
}

Double_t AliHLTConfMapper::CpuTime()
{
  //Return the Cputime in seconds.
 struct timeval tv;
 gettimeofday( &tv, NULL );
 return tv.tv_sec+(((Double_t)tv.tv_usec)/1000000.);
}
