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
#include "AliHLTITSTrack.h" 
#include "AliHLTITSDetector.h"
#include "AliHLTITSLayer.h"

#include <TObjArray.h>

#include "AliITSRecPoint.h"
#include "AliTracker.h"
#include "AliHLTITSTrack.h"
#include <vector>

//-------------------------------------------------------------------------
class AliITStrackerHLT : public AliTracker {
public:

  
  void StartLoadClusters( Int_t guessForNClusters=0 );
  void LoadCluster( const AliITSRecPoint &cluster);
  void Reconstruct( std::vector<AliExternalTrackParam> tracksTPC );

  std::vector< AliHLTITSTrack > &Tracks(){ return fTracks;}

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
  

  TTreeSRedirector *GetDebugStreamer() {return fDebugStreamer;}
  static Int_t CorrectForTPCtoITSDeadZoneMaterial(AliHLTITSTrack *t);


  AliHLTITSLayer    & GetLayer(Int_t layer) const;
  AliHLTITSDetector & GetDetector(Int_t layer, Int_t n) const {return GetLayer(layer).GetDetector(n); }
 
  void FollowProlongationTree(AliHLTITSTrack * otrack);

  void Init();


protected:

  const AliITSRecoParam *GetRecoParam(){ return fRecoParam; }
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



// method to be used for Plane Efficiency evaluation

  // 

  const AliITSRecoParam *fRecoParam;

  AliHLTITSLayer* fLayers; //!
  
  AliESDEvent  * fEsd;                   //! pointer to the ESD event
  Double_t fSPDdetzcentre[4];            // centres of SPD modules in z
  
  Int_t fUseTGeo;                        // use TGeo to get material budget

  Float_t fxOverX0Pipe;                  // material budget
  Float_t fxTimesRhoPipe;                // material budget
  Float_t fxOverX0Shield[2];             // material budget
  Float_t fxTimesRhoShield[2];           // material budget
  Float_t fxOverX0Layer[6];              // material budget
  Float_t fxTimesRhoLayer[6];            // material budget

  TTreeSRedirector *fDebugStreamer;      //!debug streamer
  AliITSChannelStatus *fITSChannelStatus;//! bitmaps with channel status for SPD and SDD
  std::vector< AliHLTITSTrack > fTracks;

  double fLoadTime;
  double fRecoTime;
  int fNEvents;
  std::vector<AliITSRecPoint> fClusters;

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
