#ifndef ALIL3_ConfMapper
#define ALIL3_ConfMapper

#include "AliL3RootTypes.h"

#include "AliL3SpacePointData.h"

class AliL3ConfMapPoint;
class AliL3ConfMapTrack;
class AliL3Vertex;
class AliL3TrackArray;

struct AliL3ConfMapContainer 
{
  void *first;
  void *last;
};

class AliL3ConfMapper {

 private:

  Bool_t fBench; //run-time measurements)

  Int_t nTracks; //number of tracks build.

  //  AliTPCParam *fParam;   
  AliL3Vertex *fVertex; //!
  Bool_t fParamSet;
  Bool_t fVertexFinder;  //Include vertexfinding or not (latter case vertex=(0,0,0))

  AliL3ConfMapPoint *fHit;  //!
  AliL3TrackArray *fTrack;  //!
  Double_t fMaxDca;      //cut value for momentum fit
  
  AliL3ConfMapContainer *fVolume;  //!  Segment volume
  AliL3ConfMapContainer *fRow;     //!  Row volume

   //Number of cells (segments)
  Int_t  fNumRowSegment;          // Total number of padrows
  Int_t  fNumPhiSegment;          // number of phi segments 
  Int_t  fNumEtaSegment;          // number of eta segments
  Int_t  fNumRowSegmentPlusOne;
  Int_t  fNumPhiSegmentPlusOne;
  Int_t  fNumEtaSegmentPlusOne;
  Int_t  fNumPhiEtaSegmentPlusOne;
  Int_t  fBounds;         
  Int_t  fPhiHitsOutOfRange;
  Int_t  fEtaHitsOutOfRange;

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
  Int_t fMaxDist;                 //maximum distance between two hits 
  Double_t fHitChi2Cut;           //Maximum hit chi2
  Double_t fGoodHitChi2;          //Chi2 to stop looking for next hit
  Double_t fTrackChi2Cut;         //Maximum track chi2
  Double_t fGoodDist;             //In segment building, distance consider good enough
  Double_t fMaxPhi;
  Double_t fMaxEta;

  // Tracking informtion
  Int_t fMainVertexTracks; //number of tracks coming from the main vertex
  Int_t fClustersUnused;   //number of unused clusters

  static Double_t todeg;
  static Double_t pi;
  static Double_t twopi;

  //setter:
  void SetMinPoints(Int_t f,Bool_t vertex_constraint) {fMinPoints[(Int_t)vertex_constraint] = f; }  
  void SetVertexConstraint(Bool_t f) {fVertexConstraint =f;}
  
  void SetTrackletLength(Int_t f,Bool_t vert) {fTrackletLength[(Int_t)vert]=f;}
  void SetRowScopeTrack(Int_t f, Bool_t vertex_constraint)         { fRowScopeTrack[(Int_t)vertex_constraint] = f; }     // sets one row scope for tracks
  void SetRowScopeTracklet(Int_t f, Bool_t vertex_constraint)      { fRowScopeTracklet[(Int_t)vertex_constraint] = f; }  // sets one row scope for tracklets
  void SetMaxAngleTracklet(Double_t f, Bool_t vertex_constraint)   { fMaxAngleTracklet[(Int_t)vertex_constraint] = f; }  // sets one angle cut

  void SetPointers();
  Double_t CpuTime();

 public:

  AliL3ConfMapper();
  //  AliL3ConfMapper(AliTPCParam *param,AliL3Vertex *vertex,Bool_t bench=(Bool_t)false);
  virtual ~AliL3ConfMapper();
  
  void InitSector(Int_t sector,Int_t *rowrange=0,Float_t *etarange=0);
  void SetVertex(AliL3Vertex *vertex){fVertex = vertex;}
  void MainVertexTracking_a();
  void MainVertexTracking_b();
  void MainVertexTracking();
  void NonVertexTracking();
  void MainVertexSettings(Int_t trackletlength, Int_t tracklength, 
			  Int_t rowscopetracklet, Int_t rowscopetrack,Double_t maxphi=0.1,Double_t maxeta=0.1);
  void NonVertexSettings(Int_t trackletlength, Int_t tracklength, 
			 Int_t rowscopetracklet, Int_t rowscopetrack);
  Bool_t ReadHits(UInt_t count, AliL3SpacePointData* hits );
  void ClusterLoop();
  void CreateTrack(AliL3ConfMapPoint *hit);
  AliL3ConfMapPoint *GetNextNeighbor(AliL3ConfMapPoint *start_hit,AliL3ConfMapTrack *track=NULL);
  Int_t EvaluateHit(AliL3ConfMapPoint *start_hit,AliL3ConfMapPoint *hit,AliL3ConfMapTrack *track);

  Double_t CalcDistance(const AliL3ConfMapPoint *hit1,const AliL3ConfMapPoint *hit2) const;
  Double_t TrackletAngle(const AliL3ConfMapTrack *track,Int_t n=3) const;
  Bool_t VerifyRange(const AliL3ConfMapPoint *hit1,const AliL3ConfMapPoint *hit2) const;
  Int_t FillTracks();
  
  //getters
  Int_t GetNumberOfTracks() {return nTracks;}
  AliL3TrackArray *GetTracks() {return fTrack;}
  Double_t GetMaxDca()  const {return fMaxDca;}
  AliL3Vertex* GetVertex() const
	  {
	  return fVertex;
	  }

  //setters
  void SetTrackCuts(Double_t hitChi2Cut, Double_t goodHitChi2, Int_t trackChi2Cut, Int_t maxdist); 
  void SetTrackletCuts(Double_t maxangle,Double_t goodDist,Bool_t vertex_constraint);   //Set cut of tracklet for the given vertex_constraint
  void SetNSegments(Int_t f,Int_t g) {fNumPhiSegment=f,fNumEtaSegment=g;} //Set number of subvolumes (#segments in (phi,eta)
  void SetParamDone(Bool_t f) {fParamSet = (Bool_t)f;}
  void SetMaxDca(Double_t f) {fMaxDca = f;}

  ClassDef(AliL3ConfMapper,1) //Base class for conformal mapping tracking
};

#endif
