//========================================================================  
// Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.  
// See cxx source for full Copyright notice                                
//========================================================================  
///                       
/// \class AliEMCALTracker 
/// \brief Steer EMCal-Track matching
///
/// Implementation of the track matching method between barrel tracks and
/// EMCAL clusters.
/// Besides algorithm implementation, some cuts are required to be set
/// in order to define, for each track, an acceptance window where clusters
/// are searched to find best match (if any).
/// The class accepts as input an ESD container, and works directly on it,
/// simply setting, for each of its tracks, the fEMCALindex flag, for each
/// track which is matched to a cluster.
/// In order to use method, one must launch PropagateBack().
///
/// \author: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
/// \author: Rong Rong Ma, Yale: Adapt to data analysis tools and other fixes.
/// \author: Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS, Run2 fixes
///
//=========================================================================

#ifndef ALIEMCALTRACKER_H
#define ALIEMCALTRACKER_H

#include "AliTracker.h"
#include <TMath.h>
#include <TVector3.h>
class TList;
class TTree;
class TObjArray;
class AliESDEvent;
class AliVCluster;
class AliESDCaloCluster;
class AliEMCALRecPoint;
class AliEMCALGeometry;

class AliEMCALTracker : public AliTracker 
{
 public:
  AliEMCALTracker();
  AliEMCALTracker(const AliEMCALTracker &t);
  AliEMCALTracker& operator=(const AliEMCALTracker &source);
  virtual ~AliEMCALTracker() {Clear();}
	
  virtual void        Clear(Option_t *option="ALL");

  virtual Int_t       Clusters2Tracks(AliESDEvent*) { return -1   ; }

          void        InitParameters();

  virtual Int_t       LoadClusters(TTree*);
          Int_t       LoadClusters(AliESDEvent* esd);
          Int_t       LoadTracks  (AliESDEvent* esd);
  virtual Int_t       PropagateBack(AliESDEvent* esd);
  virtual Int_t       RefitInward(AliESDEvent*)     { return -1   ; }
  virtual void        UnloadClusters();

  virtual AliCluster* GetCluster   (Int_t)    const { return NULL ; }

          void        SetCutEta    (Double_t value) { fCutEta     = value ; }
	  void        SetCutPhi    (Double_t value) { fCutPhi     = value ; }
	  void        SetCutPt     (Double_t value) { fCutPt      = value ; }
	  void        SetCutNITS   (Double_t value) { fCutNITS    = value ; }
	  void        SetCutNTPC   (Double_t value) { fCutNTPC    = value ; }
	  void        SetTrackInITS(Bool_t   value) { fTrackInITS = value ; }
	  void        SetStepLength(Float_t length) { fStep       = length; }

	  void        SetTrackCorrectionMode(Option_t *option);
	  void        SetEMCalSurfaceDistance(Double_t d) {fEMCalSurfaceDistance = d;}
	  void        SetGeometry  (AliEMCALGeometry *geom) { fGeom = geom ; }

  enum {	kUnmatched = -99999 };
	
  class  AliEMCALMatchCluster : public TObject
  {
   public:
    AliEMCALMatchCluster(Int_t ID, AliEMCALRecPoint *recPoint);
    AliEMCALMatchCluster(Int_t ID, AliESDCaloCluster *caloCluster);
    virtual ~AliEMCALMatchCluster() { }
    //----------------------------------------------------------------------------
    Int_t     Index() const {return fIndex;}
    Double_t  X() const {return fX;}
    Double_t  Y() const {return fY;} 
    Double_t  Z() const {return fZ;}
   private:
    Int_t     fIndex;  ///< Index of cluster in its native container (ESD or TClonesArray)
    Double_t  fX;      ///< Global X position
    Double_t  fY;      ///< Global Y position
    Double_t  fZ;      ///< Global Z position
  };
   
 private:
  Int_t  FindMatchedCluster(AliESDtrack *track);
  enum ETrackCorr { 
    kTrackCorrNone  = 0,        ///< Do not correct for energy loss
    kTrackCorrMMB   = 1,        ///< Use MeanMaterialBudget() function to evaluate correction
  };

  Double_t    fCutPt;           ///< Minimum pT cut on tracks
  Double_t    fCutNITS;         ///< Minimum number of track hits in the ITS
  Double_t    fCutNTPC;         ///< Minimum number of track hits in the TPC
	
  Float_t     fStep;            ///< Length of each step in propagation
  ETrackCorr  fTrackCorrMode;   ///< Material budget correction mode
  Float_t     fEMCalSurfaceDistance; ///< EMCal surface distance
  Double_t    fClusterWindow;   ///< Select clusters in the window to be matched to tracks
  Double_t    fCutEta;          ///< cut on eta difference
  Double_t    fCutPhi;          ///< cut on phi difference
  Bool_t      fITSTrackSA;      ///< If no TPC, use ITS Tracks	
  Bool_t      fTrackInITS;      ///< Requiere ITS refit with AliVTrack::kITSout	
	
  TObjArray  *fTracks;          //!<! Collection of tracks
  TObjArray  *fClusters;        //!<! Collection of EMCAL clusters (ESDCaloCluster or EMCALRecPoint)
	
  AliEMCALGeometry *fGeom;      //!<! EMCAL geometry
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALTracker, 8) ;
  /// \endcond
 
};
#endif
