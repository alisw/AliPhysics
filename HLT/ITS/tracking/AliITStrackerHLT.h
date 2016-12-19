#ifndef ALIITSTRACKERHLT_H
#define ALIITSTRACKERHLT_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


class TTree;
class TTreeSRedirector;
class AliESDEvent;
class AliESDtrack;

class AliITSChannelStatus;
class AliITSDetTypeRec;
class AliITSRecoParam;
class AliHLTITSTrackPoint;
#include "AliHLTITSTrack.h" 
#include "AliHLTITSDetector.h"
#include "AliHLTITSLayer.h"

#include <TObjArray.h>

#include "AliITSRecPoint.h"
#include "AliTracker.h"
#include "AliHLTITSTrack.h"

//-------------------------------------------------------------------------
class AliITStrackerHLT : public AliTracker {
public:

  
  void StartLoadClusters( Int_t NClusters );
  void LoadCluster( const AliITSRecPoint &cluster);
  void Reconstruct( AliExternalTrackParam *tracksTPC, int *tracksTPCLab, int nTPCTracks );

  AliHLTITSTrack *Tracks() const { return fTracks;}
  Int_t NTracks() const { return fNTracks;}
  AliHLTITSTrack *ITSOutTracks() const { return fITSOutTracks;}
  Int_t NITSOutTracks() const { return fNITSOutTracks;}

  Bool_t TransportToX( AliExternalTrackParam *t, double x ) const;
  Bool_t TransportToPhiX( AliExternalTrackParam *t, double phi, double x ) const;

  AliITStrackerHLT();
  AliITStrackerHLT(const Char_t *geom);
  virtual ~AliITStrackerHLT();
  AliCluster *GetCluster(Int_t index) const;
  virtual Bool_t GetTrackPoint(Int_t index, AliTrackPoint& p) const;
  virtual Bool_t GetTrackPointTrackingError(Int_t index, 
			AliTrackPoint& p, const AliESDtrack *t);
  AliITSRecPoint *GetClusterLayer(Int_t layn, Int_t ncl) const
                  {return fLayers[layn].GetCluster(ncl);}
  Int_t GetNumberOfClustersLayer(Int_t layn) const 
                        {return fLayers[layn].GetNumberOfClusters();}
  Int_t LoadClusters(TTree *cf);
  void UnloadClusters();
  
  Int_t Clusters2Tracks(AliESDEvent *event);
  Int_t PropagateBack(AliESDEvent *event);
  Int_t RefitInward(AliESDEvent *event);
  

  static Int_t CorrectForTPCtoITSDeadZoneMaterial(AliHLTITSTrack *t);


  AliHLTITSLayer    & GetLayer(Int_t layer) const;
  AliHLTITSDetector & GetDetector(Int_t layer, Int_t n) const {return GetLayer(layer).GetDetector(n); }
 
  void FollowProlongationTree(AliHLTITSTrack * otrack);
  Int_t FitOutward(AliHLTITSTrack * track );

  void Init();

 // track point calculation

  Int_t GetTrackPoint( Int_t clusterIndex, AliHLTITSTrackPoint& p ) const ;

protected:

  const AliITSRecoParam *GetRecoParam() const { return fRecoParam; }
  Bool_t ComputeRoad(AliHLTITSTrack* track,Int_t ilayer,Int_t idet,Double_t &zmin,Double_t &zmax,Double_t &ymin,Double_t &ymax) const;
  
  
  void CookLabel(AliKalmanTrack *t,Float_t wrong) const;
  void CookLabel(AliHLTITSTrack *t,Float_t wrong) const;

  void BuildMaterialLUT(TString material);
  
  Int_t CorrectForPipeMaterial(AliHLTITSTrack *t, bool InwardDirection=1);
  Int_t CorrectForShieldMaterial(AliHLTITSTrack *t, Int_t  shieldindex, bool InwardDirection=1);
  Int_t CorrectForLayerMaterial(AliHLTITSTrack *t, Int_t layerindex, bool InwardDirection=1);
  void UpdateESDtrack(AliESDtrack *tESD,AliHLTITSTrack* track, ULong_t flags) const;
  
  Bool_t LocalModuleCoord(Int_t ilayer,Int_t idet,const AliHLTITSTrack *track,
			  Float_t &xloc,Float_t &zloc) const;

  static Bool_t CheckTrack( const AliExternalTrackParam *t );


  // 

  const AliITSRecoParam *fRecoParam; //! 

  AliHLTITSLayer* fLayers; //!
  
  Double_t fSPDdetzcentre[4];            // centres of SPD modules in z
  
  Int_t fUseTGeo;                        // use TGeo to get material budget

  Float_t fxOverX0Pipe;                  // material budget
  Float_t fxTimesRhoPipe;                // material budget
  Float_t fxOverX0Shield[2];             // material budget
  Float_t fxTimesRhoShield[2];           // material budget
  Float_t fxOverX0Layer[6];              // material budget
  Float_t fxTimesRhoLayer[6];            // material budget

  AliHLTITSTrack *fTracks; // array of its-updated tracks
  AliHLTITSTrack *fITSOutTracks; // array of tracks, fitted outward with ITS only
  int fNTracks;// n tracks
  int fNITSOutTracks;// n out tracks
  double fLoadTime; //
  double fRecoTime; //
  int fNEvents; //
  AliITSRecPoint *fClusters; //!
  int fNClusters; //

private:
  AliITStrackerHLT(const AliITStrackerHLT &tracker);
  AliITStrackerHLT & operator=(const AliITStrackerHLT &tracker);  
  ClassDef(AliITStrackerHLT,0)   //HLT ITS tracker
};




/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////





inline void AliITStrackerHLT::CookLabel(AliKalmanTrack *t,Float_t wrong) const {
  //--------------------------------------------------------------------
  //This function "cooks" a track label. If label<0, this track is fake.
  //--------------------------------------------------------------------
   Int_t tpcLabel=t->GetLabel();
   if (tpcLabel<0) return;
   AliTracker::CookLabel(t,wrong);
   if (tpcLabel!=TMath::Abs(t->GetLabel())){
     t->SetFakeRatio(1.);
   }
   if (tpcLabel !=t->GetLabel()) {
     t->SetLabel(-tpcLabel);      
   }
}



#endif
