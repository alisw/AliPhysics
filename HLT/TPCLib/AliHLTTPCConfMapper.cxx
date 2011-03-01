// @(#) $Id$
// Original: AliHLTConfMapper.cxx,v 1.26 2005/06/14 10:55:21 cvetan Exp $

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Anders Vestbo, maintained by
//*                  Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTPCConfMapper.cxx
    @author Anders Vestbo, Matthias Richter
    @date   Conformal mapping base class.
    @brief  
*/

#include <cassert>
#include <sys/time.h>
 
#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCLogging.h" 
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCConfMapTrack.h"
#include "AliHLTTPCConfMapPoint.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCConfMapper.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCConfMapper)

AliHLTTPCConfMapper::AliHLTTPCConfMapper()
  :
  fBench(kTRUE),
  fNTracks(0),
  fVertex(NULL),
  fVertexFinder(kFALSE),
  fHit(),
  fTrack(NULL),
  fMaxDca(0.0), // no clue whether this is reasonable, but at least better than without initialization
  fVolume(NULL),
  fRow(NULL),
  fNumRowSegment(0),
  fNumPhiSegment(0),
  fNumEtaSegment(0),
  fNumRowSegmentPlusOne(0),
  fNumPhiSegmentPlusOne(0),
  fNumEtaSegmentPlusOne(0),
  fNumPhiEtaSegmentPlusOne(0),
  fBounds(0),
  fPhiHitsOutOfRange(0),
  fEtaHitsOutOfRange(0),
  fPhiMin(0),
  fPhiMax(0),
  fEtaMin(0),
  fEtaMax(0),
  fRowMin(0),
  fRowMax(0),
  fVertexConstraint(kTRUE),
  fGoodDist(0.0),
  fMaxPhi(0.0),
  fMaxEta(0.0),
  fMainVertexTracks(0),
  fClustersUnused(0),
  fClusterCutZ(-1)
{
  //Default constructor
  memset(fParamSet, 0, sizeof(fParamSet));
  memset(fTrackletLength, 0, sizeof(fTrackletLength));
  memset(fRowScopeTracklet, 0, sizeof(fRowScopeTracklet));
  memset(fRowScopeTrack, 0, sizeof(fRowScopeTrack));
  memset(fMinPoints, 0, sizeof(fMinPoints));
  
  memset(fMaxAngleTracklet, 0, sizeof(fMaxAngleTracklet));
  memset(fMaxDist, 0, sizeof(fMaxDist));
  memset(fHitChi2Cut, 0, sizeof(fHitChi2Cut));
  memset(fGoodHitChi2, 0, sizeof(fGoodHitChi2));
  memset(fTrackChi2Cut, 0, sizeof(fTrackChi2Cut));
}

AliHLTTPCConfMapper::~AliHLTTPCConfMapper()
{
  // Destructor.
  if(fRow) {
    delete [] fRow;
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

  Reset();
  
  fTrack = new AliHLTTPCTrackArray("AliHLTTPCConfMapTrack",10);
}

void AliHLTTPCConfMapper::Reset()
{
  if(fVolume) delete [] fVolume;
  fVolume=NULL;
  if(fRow) delete [] fRow;
  fRow=NULL;
  
  fClustersUnused=0;
  fHit.clear();
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
  fPhiMin = -10*AliHLTTPCTransform::ToRad();//fParam->GetAngle(sector) - 10/todeg;
  fPhiMax =  10*AliHLTTPCTransform::ToRad();//fParam->GetAngle(sector) + 10/todeg;

  fNTracks=0;
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
  //read hits with ReadHitsChecked  
  return ReadHitsChecked(count,hits,0);
}

Bool_t AliHLTTPCConfMapper::ReadHitsChecked(UInt_t count, AliHLTTPCSpacePointData* hits, unsigned int sizeInByte )
{
  //read hits
  if(fClusterCutZ == -1){
    if (fHit.size()<fClustersUnused+count) fHit.resize(fClustersUnused+count);
    assert(fHit.size()>=fClustersUnused+count);
    for (Int_t i=0;(UInt_t)i<count;i++)
      {	
	AliHLTTPCSpacePointData *hit = &hits[i];
	if (sizeInByte>0 && ((AliHLTUInt8_t*)hit)+sizeof(AliHLTTPCSpacePointData)>((AliHLTUInt8_t*)hits)+sizeInByte) {
	  LOG(AliHLTTPCLog::kWarning,"AliHLTTPCConfMapper::ReadHits","")<<"Wrong size of data (" << sizeInByte << " byte), skipping array of AliHLTTPCSpacePointData" <<ENDLOG;;
	  break;
	}
	fHit[i+fClustersUnused].Reset();
	fHit[i+fClustersUnused].Read(hits[i]);
      }
    fClustersUnused += count;
  }
  else{
    //Skipping clusters with high Z. 
    UInt_t skipped=0;
    
    if (fHit.size()<fClustersUnused+count) fHit.resize(fClustersUnused+count);
    assert(fHit.size()>=fClustersUnused+count);
    for (Int_t i=0;(UInt_t)i<count;i++)
      {  
	AliHLTTPCSpacePointData *hit = &hits[i];
	if (sizeInByte>0 && ((AliHLTUInt8_t*)hit)+sizeof(AliHLTTPCSpacePointData)>((AliHLTUInt8_t*)hits)+sizeInByte) {
	  LOG(AliHLTTPCLog::kWarning,"AliHLTTPCConfMapper::ReadHits","")<<"Wrong size of data (" << sizeInByte << " byte), skipping array of AliHLTTPCSpacePointData" <<ENDLOG;;
	  break;
	}
	if(hits[i].fZ > fClusterCutZ || hits[i].fZ < -1*fClusterCutZ){
	  ++skipped;
	  continue;
	}
	fHit[i+fClustersUnused-skipped].Reset();
	fHit[i+fClustersUnused-skipped].Read(hits[i]);
      }
    fClustersUnused += count - skipped;
    fHit.resize(fClustersUnused);
  }

  LOG(AliHLTTPCLog::kDebug,"AliHLTTPCConfMapper::ReadHits","#hits")
    <<AliHLTTPCLog::kDec<<"#hits: "<<count<<" total: "<<fClustersUnused<<ENDLOG;
  
  return true;
}

void AliHLTTPCConfMapper::SetPointers()
{
  //Check if there are not enough clusters to make a track in this sector
  //Can happen in pp events.

  if(fClustersUnused < fMinPoints[fVertexConstraint])
    return;
  
  //Allocate detector volumes
  if (fVolume==NULL) {
    LOG(AliHLTTPCLog::kDebug,"AliHLTTPCConfMapper::InitVolumes","Memory")<<AliHLTTPCLog::kDec<<
      "Allocating "<<fBounds*sizeof(AliHLTTPCConfMapContainer)<<" Bytes to fVolume"<<ENDLOG;
    fVolume = new AliHLTTPCConfMapContainer[fBounds];
  }

  if (fRow==NULL) {
    LOG(AliHLTTPCLog::kDebug,"AliHLTTPCConfMapper::InitVolumes","Memory")<<AliHLTTPCLog::kDec<<
      "Allocating "<<fNumRowSegmentPlusOne*sizeof(AliHLTTPCConfMapContainer)<<" Bytes to fRow"<<ENDLOG;
    fRow = new AliHLTTPCConfMapContainer[fNumRowSegmentPlusOne];
  }
  
  memset(fVolume,0,fBounds*sizeof(AliHLTTPCConfMapContainer));
  memset(fRow,0,fNumRowSegmentPlusOne*sizeof(AliHLTTPCConfMapContainer));
  
  Float_t phiSlice = (fPhiMax-fPhiMin)/fNumPhiSegment;
  Float_t etaSlice = (fEtaMax-fEtaMin)/fNumEtaSegment;

  Int_t volumeIndex;
  Int_t localcounter=0;
  assert((int)fHit.size()>=fClustersUnused);
  for(Int_t j=0; j<fClustersUnused; j++)
    {
      AliHLTTPCConfMapPoint *thisHit = &(fHit[j]);

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
 	((AliHLTTPCConfMapPoint *)fVolume[volumeIndex].last)->SetNextVolumeHit(thisHit);
      fVolume[volumeIndex].last = (void *)thisHit;
      
      
      //set row pointers
      if(fRow[(localrow-fRowMin)].first == NULL)
 	fRow[(localrow-fRowMin)].first = (void *)thisHit;
      else
 	((AliHLTTPCConfMapPoint *)(fRow[(localrow-fRowMin)].last))->SetNextRowHit(thisHit);
	fRow[(localrow-fRowMin)].last = (void *)thisHit;
    }
  
  //If a cluster has an Eta outside the Eta or Phi range set in the Tracker, it will go in to
  //the if here. This has been seen for high Eta clusters most likely from signal from the gating grid.
  //These clusters are read in, but not used in the Tracking. 
#ifdef PACKAGE_STRING
  if(fClustersUnused>0 && localcounter==0)
    LOG(AliHLTTPCLog::kDebug,"AliHLTTPCConfMapper::SetPointers","Parameters")
      <<AliHLTTPCLog::kDec<<"No points passed to track finder, hits out of range: "
      <<fEtaHitsOutOfRange+fPhiHitsOutOfRange<<ENDLOG;

  Int_t hits_accepted=fClustersUnused-(fEtaHitsOutOfRange+fPhiHitsOutOfRange);
  LOG(AliHLTTPCLog::kDebug,"AliHLTTPCConfMapper::SetPointers","Setup")
    <<"Setup finished, hits out of range: "<<fEtaHitsOutOfRange+fPhiHitsOutOfRange
    <<" hits accepted "<<hits_accepted<<ENDLOG;
#endif //PACKAGE_STRING
}

void AliHLTTPCConfMapper::MainVertexTrackingA()
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
    LOG(AliHLTTPCLog::kBenchmark,"AliHLTTPCConfMapper::MainVertexTrackingA","Timing")
      <<AliHLTTPCLog::kDec<<"Setup finished in "<<cpuTime*1000<<" ms"<<ENDLOG;
  
}

void AliHLTTPCConfMapper::MainVertexTrackingB()
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
    LOG(AliHLTTPCLog::kBenchmark,"AliHLTTPCConfMapper::MainVertexTrackingB","Timing")
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
    LOG(AliHLTTPCLog::kBenchmark,"AliHLTTPCConfMapper::MainVertexTracking","Timing")<<AliHLTTPCLog::kDec<<
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
  
  SetPointers(); //To be able to do only nonvertextracking (more testing) 
  
  ClusterLoop();
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCConfMapper::NonVertexTracking","ntracks")<<AliHLTTPCLog::kDec<<
    "Number of nonvertex tracks found: "<<(fNTracks-fMainVertexTracks)<<ENDLOG;
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
  //set parameters for non-vertex tracking
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

void AliHLTTPCConfMapper::SetTrackletCuts(Double_t maxangle,Double_t goodDist, Bool_t vc)
{
  //Sets cuts of tracklets. Right now this is only:
  //maxangle:  Maximum angle when forming segments (if trackletlength > 2)
 
  fGoodDist=goodDist;
  SetMaxAngleTracklet(maxangle, vc);
}

void AliHLTTPCConfMapper::ClusterLoop()
{
  //Loop over hits, starting at outermost padrow, and trying to build segments.
  
  //Check if there are not enough clusters to make a track in this sector
  //Can happen in pp events.
  if(fClustersUnused < fMinPoints[fVertexConstraint])
    return;
  
  Int_t rowsegm,lastrow = fRowMin + fMinPoints[fVertexConstraint];
  AliHLTTPCConfMapPoint *hit;
  
  //Loop over rows, and try to create tracks from the hits.
  //Starts at the outermost row, and loops as long as a track can be build, due to length.
  
  for(rowsegm = fRowMax; rowsegm >= lastrow; rowsegm--)
    {
      if(fRow[(rowsegm-fRowMin)].first && ((AliHLTTPCConfMapPoint*)fRow[(rowsegm-fRowMin)].first)->GetPadRow() < fRowMin + 1)
	break;

      for(hit = (AliHLTTPCConfMapPoint*)fRow[(rowsegm-fRowMin)].first; hit!=0; hit=hit->GetNextRowHit())
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

  AliHLTTPCConfMapPoint *closesthit = NULL;
  AliHLTTPCConfMapTrack *track = NULL;
  
  Int_t point;
  Int_t tracks = fNTracks;
  fNTracks++;

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
      if((closesthit = GetNextNeighbor(hit)))
	{//closest hit exist
	  
	  //   Calculate track length in sz plane
	  dx = ((AliHLTTPCConfMapPoint*)closesthit)->GetX() - ((AliHLTTPCConfMapPoint*)hit)->GetX();
	  dy = ((AliHLTTPCConfMapPoint*)closesthit)->GetY() - ((AliHLTTPCConfMapPoint*)hit)->GetY();
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
	      closesthit = GetNextNeighbor((AliHLTTPCConfMapPoint*)track->GetLastHit(),track);
	      
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

AliHLTTPCConfMapPoint *AliHLTTPCConfMapper::GetNextNeighbor(AliHLTTPCConfMapPoint *starthit,
					  AliHLTTPCConfMapTrack *track)
{
  //When forming segments: Finds closest hit to input hit
  //When forming tracks: Find closest hit to track fit.
  
  Double_t dist,closestdist = fMaxDist[fVertexConstraint];
  
  AliHLTTPCConfMapPoint *hit = NULL;
  AliHLTTPCConfMapPoint *closesthit = NULL;
    
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
		  LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapper::GetNextNeighbor","Memory")<<AliHLTTPCLog::kDec<<
		    "VolumeIndex error "<<volumeIndex<<ENDLOG;
		}
	      
	      assert(fVolume!=NULL);
	      for(hit = (AliHLTTPCConfMapPoint*)fVolume[volumeIndex].first;
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
			  if (fVertexConstraint) {   
			    if(!VerifyRange(starthit,hit))
			      continue;
			  }
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
			      if (fVertexConstraint) {   
				if(!VerifyRange(starthit,hit))
				  continue;
			      }
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

Int_t AliHLTTPCConfMapper::EvaluateHit(AliHLTTPCConfMapPoint *starthit,AliHLTTPCConfMapPoint *hit,AliHLTTPCConfMapTrack *track) 
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

Double_t AliHLTTPCConfMapper::CalcDistance(const AliHLTTPCConfMapPoint *hit1,const AliHLTTPCConfMapPoint *hit2) const
{
  //Return distance between two clusters, defined by Pablo
  
  Double_t phidiff = fabs( hit1->GetPhi() - hit2->GetPhi() );
  if (phidiff > AliHLTTPCTransform::Pi()) phidiff = AliHLTTPCTransform::TwoPi() - phidiff;
  
  return AliHLTTPCTransform::ToDeg()*fabs((Float_t)((hit1->GetPadRow() - hit2->GetPadRow()) * 
         (phidiff + fabs( hit1->GetEta() - hit2->GetEta()))));
}

Bool_t AliHLTTPCConfMapper::VerifyRange(const AliHLTTPCConfMapPoint *hit1,const AliHLTTPCConfMapPoint *hit2) const
{
  //Check if the hit are within reasonable range in phi and eta
  Double_t dphi,deta;//maxphi=0.1,maxeta=0.1;
  dphi = fabs(hit1->GetPhi() - hit2->GetPhi());
  if(dphi > AliHLTTPCTransform::Pi()) dphi = fabs(AliHLTTPCTransform::TwoPi() - dphi);
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
  
  Double_t x1[2]={0,0};
  Double_t x2[2]={0,0};
  Double_t x3[2]={0,0};
  Double_t angle1,angle2;
  Int_t counter=0;
  for(track->StartLoop(); track->LoopDone(); track->GetNextHit())
    {
      AliHLTTPCConfMapPoint *p = (AliHLTTPCConfMapPoint*)track->GetCurrentHit();
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
  
  Int_t numoftracks = fNTracks;
  if(fNTracks == 0)
    {
      LOG(AliHLTTPCLog::kDebug,"AliHLTTPCConfMapper::FillTracks","fNTracks")<<AliHLTTPCLog::kDec<<
	"No tracks found!!"<<ENDLOG;
      return 0;
    }

  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCConfMapper::FillTracks","fNTracks")<<AliHLTTPCLog::kDec<<
    "Number of found tracks: "<<fNTracks<<ENDLOG;
  
  //  fTrack->Sort();
  for(Int_t i=0; i<numoftracks; i++)
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
}
