// @(#) $Id$
// Original: AliHLTConfMapper.h,v 1.11 2004/07/05 09:03:11 loizides 

#ifndef ALIHLTTPCCONFMAPPER_H
#define ALIHLTTPCCONFMAPPER_H

//
//Conformal mapping base class
//
// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group


class AliHLTTPCConfMapPoint;
class AliHLTTPCConfMapTrack;
class AliHLTTPCVertex;
class AliHLTTPCTrackArray;
class AliHLTTPCSpacePointData;


class AliHLTTPCConfMapper {

 public:

  AliHLTTPCConfMapper();
  //  AliHLTTPCConfMapper(AliTPCParam *param,AliHLTTPCVertex *vertex,Bool_t bench=(Bool_t)false);
  virtual ~AliHLTTPCConfMapper();
  
  void InitVolumes();
  void InitSector(Int_t sector,Int_t *rowrange=0,Float_t *etarange=0);
  void SetVertex(AliHLTTPCVertex *vertex){fVertex = vertex;}
  void MainVertexTracking_a();
  void MainVertexTracking_b();
  void MainVertexTracking();
  void NonVertexTracking();
  void MainVertexSettings(Int_t trackletlength, Int_t tracklength, 
			  Int_t rowscopetracklet, Int_t rowscopetrack,Double_t maxphi=0.1,Double_t maxeta=0.1);
  void NonVertexSettings(Int_t trackletlength, Int_t tracklength, 
			 Int_t rowscopetracklet, Int_t rowscopetrack);
  Bool_t ReadHits(UInt_t count, AliHLTTPCSpacePointData* hits );
  void ClusterLoop();
  void CreateTrack(AliHLTTPCConfMapPoint *hit);
  AliHLTTPCConfMapPoint *GetNextNeighbor(AliHLTTPCConfMapPoint *start_hit,AliHLTTPCConfMapTrack *track=NULL);
  Int_t EvaluateHit(AliHLTTPCConfMapPoint *start_hit,AliHLTTPCConfMapPoint *hit,AliHLTTPCConfMapTrack *track);

  Double_t CalcDistance(const AliHLTTPCConfMapPoint *hit1,const AliHLTTPCConfMapPoint *hit2) const;
  Double_t TrackletAngle(AliHLTTPCConfMapTrack *track,Int_t n=3) const;
  Bool_t VerifyRange(const AliHLTTPCConfMapPoint *hit1,const AliHLTTPCConfMapPoint *hit2) const;
  Int_t FillTracks();
  
  //getters
  Int_t GetNumberOfTracks()    const {return fNTracks;}
  AliHLTTPCTrackArray *GetTracks() const {return fTrack;}
  Double_t GetMaxDca()         const {return fMaxDca;}
  AliHLTTPCVertex* GetVertex()     const {return fVertex;}
  Double_t CpuTime();

  //setters
  void SetTrackCuts(Double_t hitChi2Cut, Double_t goodHitChi2, Double_t trackChi2Cut, Int_t maxdist,Bool_t vertexconstraint); 
  void SetTrackletCuts(Double_t maxangle,Double_t goodDist,Bool_t vc);    //Set cut of tracklet for the given vertex_constraint
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

  void SetPointers();
  void SetParamDone(Bool_t vconstraint) {fParamSet[(Int_t)vconstraint] = kTRUE;}

 private:

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

  ClassDef(AliHLTTPCConfMapper,1) //Base class for conformal mapping tracking
};

#endif
