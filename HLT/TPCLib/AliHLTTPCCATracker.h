// @(#) $Id$
// Original: AliL3ConfMapper.h,v 1.11 2004/07/05 09:03:11 loizides 

#ifndef ALIHLTTPCCATRACKER_H
#define ALIHLTTPCCATRACKER_H

//
// CA Tracking class 
//
// Author: Ivan Kisel 
//*-- Copyright &copy ALICE HLT Group

#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"

#include "AliHLTTPCTrackSegmentData.h"

class AliHLTTPCConfMapPoint;
class AliHLTTPCConfMapTrack;
class AliHLTTPCVertex;
class AliHLTTPCTrackArray;
class AliHLTTPCSpacePointData;

class AliHLTTPCCATracker {

 public:

  AliHLTTPCCATracker();
  //  AliHLTTPCCATracker(AliTPCParam *param,AliHLTTPCVertex *vertex,Bool_t bench=(Bool_t)false);
  virtual ~AliHLTTPCCATracker();

  Bool_t ReadHits(UInt_t count, AliHLTTPCSpacePointData* hits );

  //getters
#if 0
  Int_t GetNumberOfTracks()    const {return fNTracks;}
  AliHLTTPCTrackArray *GetTracks() const {return fTrack;}
  Double_t GetMaxDca()         const {return fMaxDca;}
  AliHLTTPCVertex* GetVertex()     const {return fVertex;}
#endif


  //setters
#if 0
  void SetNSegments(Int_t f,Int_t g) {fNumPhiSegment=f,fNumEtaSegment=g;} //Set number of subvolumes (#segments in (phi,eta)
  void SetMaxDca(Double_t f) {fMaxDca = f;}

  //setter:
  void SetMinPoints(Int_t f,Bool_t vconstraint) {fMinPoints[(Int_t)vconstraint] = f; }  
  void SetVertexConstraint(Bool_t f) {fVertexConstraint =f;}
  
  void SetHitChi2Cut(Double_t f,Bool_t vert) {fHitChi2Cut[(Int_t)vert]=f;}
  void SetGoodHitChi2(Double_t f,Bool_t vert) {fGoodHitChi2[(Int_t)vert]=f;}
  void SetTrackChi2Cut(Double_t f,Bool_t vert) {fTrackChi2Cut[(Int_t)vert]=f;}
  void SetMaxDist(Int_t f,Bool_t vert) {fMaxDist[(Int_t)vert]=f;}
  void SetTrackletLength(Int_t f,Bool_t vert) {fTrackletLength[(Int_t)vert]=f;}
  void SetRowScopeTrack(Int_t f, Bool_t vc){fRowScopeTrack[(Int_t)vc] = f;}
  void SetRowScopeTracklet(Int_t f, Bool_t vc){fRowScopeTracklet[(Int_t)vc] = f;}
  void SetMaxAngleTracklet(Double_t f, Bool_t vc){fMaxAngleTracklet[(Int_t)vc] = f;}
#endif

  void CACreateHistos();
  void CAWriteHistos();

  void CAInitialize();
  void CAReadPatchHits(Int_t patch, UInt_t count, AliHLTTPCSpacePointData* hits );
  void CAFindPatchTracks(Int_t patch);
  void CAFindSliceTracks();


  // JMT 2006/11/13
  void SetOutPtr( AliHLTTPCTrackSegmentData* tr ){fOutputPtr = tr;}
  UInt_t GetOutputSize() { return fOutputSize; }
  UInt_t GetOutputNTracks() { return fOutputNTracks; }


 private:

  AliHLTTPCTrackSegmentData* fOutputPtr;
  UInt_t fOutputNTracks;
  UInt_t fOutputSize;

  struct CAHit{
    Double_t x, y, z;
    Double_t errx, erry, errz;
    Int_t index, counter;
  };
  
  std::vector<CAHit> vec_hits; 
  Int_t patch_first_hit_ind, patch_last_hit_ind; // indices of the first and the last hits in the current patch

  struct CATrack{
    Int_t patch, nhits, ndf, good, used, next;
    //parameters
    Double_t x, y, z, ty, tz;
    //cov matrix
    Double_t cov_y, cov_ty, cov_yty, cov_z, cov_tz, cov_ztz, chi2;
    //indices of hits
    std::vector<Int_t> vec_ihits;
  };

  std::vector<CATrack> vec_patch_tracks; 
  Int_t patch_first_track_ind, patch_last_track_ind; // indices of the first and the last tracks in the current patch
  std::vector<CATrack> vec_slice_tracks; 

  static bool compareCATracks(const CATrack &a, const CATrack &b){
    if (a.patch != b.patch)
      return (a.patch < b.patch);
    else
      return (a.nhits > b.nhits);
  }

  static bool compareCAHitsX(const Int_t &i, const Int_t &j){
    return (i < j);
  }

  struct AliHLTTPCConfMapContainer 
  {
    void *first; // first track
    void *last;  // last track
  };

  Bool_t fBench; //run-time measurements
  Int_t fNTracks; //number of tracks build.

  AliHLTTPCVertex *fVertex; //!
  Bool_t fParamSet[2];  //!
  Bool_t fVertexFinder; //Include vertexfinding or not 
                        //(latter case vertex=(0,0,0))

  AliHLTTPCConfMapPoint *fHit;  //!
  AliHLTTPCTrackArray *fTrack;  //!
  Double_t fMaxDca;      //cut value for momentum fit
  
  AliHLTTPCConfMapContainer *fVolume;  //!  Segment volume
  AliHLTTPCConfMapContainer *fRow;     //!  Row volume

   //Number of cells (segments)
  Int_t  fNumRowSegment;          // Total number of padrows
  Int_t  fNumPhiSegment;          // number of phi segments 
  Int_t  fNumEtaSegment;          // number of eta segments
  Int_t  fNumRowSegmentPlusOne;   // row+1
  Int_t  fNumPhiSegmentPlusOne;   // phi+1
  Int_t  fNumEtaSegmentPlusOne;   // eta+1
  Int_t  fNumPhiEtaSegmentPlusOne;// phieta+1
  Int_t  fBounds;                 // bounds
  Int_t  fPhiHitsOutOfRange;      // phi hits out of range
  Int_t  fEtaHitsOutOfRange;      // eta hits out of range

  //tracking range:
  Float_t fPhiMin; //MinPhi angle to consider
  Float_t fPhiMax; //MaxPhi angle to consider
  Float_t fEtaMin; //MinEta to consider
  Float_t fEtaMax; //MaxEta to consider
  Int_t fRowMin;   //Minimum row to consider
  Int_t fRowMax;   //Maximum row to consider

  Bool_t fVertexConstraint;       //vertex constraint (true or false)
  Int_t fTrackletLength[2];       //minimal length of tracks 
  Int_t fRowScopeTracklet[2];     //number of row segments to look for the next point of a tracklet
  Int_t fRowScopeTrack[2];        //number of row segments to look for the next point of a track
  Int_t fMinPoints[2];            //minimum number of points on one track
  
  // Cuts
  Double_t fMaxAngleTracklet[2];  //limit of angle between to pieces of a tracklet
  Int_t fMaxDist[2];              //maximum distance between two hits 
  Double_t fHitChi2Cut[2];        //Maximum hit chi2
  Double_t fGoodHitChi2[2];       //Chi2 to stop looking for next hit
  Double_t fTrackChi2Cut[2];      //Maximum track chi2
  Double_t fGoodDist;             //In segment building, distance consider good enough
  Double_t fMaxPhi;               //Maximum phi
  Double_t fMaxEta;               //Maximum eta

  // Tracking informtion
  Int_t fMainVertexTracks; //number of tracks coming from the main vertex
  Int_t fClustersUnused;   //number of unused clusters

  ClassDef(AliHLTTPCCATracker,1) //Base class for conformal mapping tracking
};

#endif
