#ifndef ALIEMCALRECOUTILSBASE_H
#define ALIEMCALRECOUTILSBASE_H

///////////////////////////////////////////////////////////////////////////////
///
/// \class AliEMCALRecoUtilsBase
/// \ingroup EMCALUtils
/// \brief Some utilities for track matching.
///
/// This is the base class to correct and select the clusters and cells at analysis level.
/// In this class only track-matching related methods.
/// Most of the corrections implemented in AliEMCALRecoUtils (in AliPhysics from v5-09-26)
///
/// \author:  Gustavo Conesa Balbastre, <Gustavo.Conesa.Balbastre@cern.ch>, LPSC- Grenoble 
/// \author:  Rongrong Ma, Yale. Track matching part
///
///////////////////////////////////////////////////////////////////////////////

// Root includes
#include <TNamed.h>

// AliRoot includes
class AliVCluster;

// EMCAL includes
class AliEMCALGeometry;
class AliExternalTrackParam;
class AliVTrack;

class AliEMCALRecoUtilsBase : public TNamed {
  
public:
  
  AliEMCALRecoUtilsBase();
  
  AliEMCALRecoUtilsBase(           const AliEMCALRecoUtilsBase&); 
  
  AliEMCALRecoUtilsBase& operator=(const AliEMCALRecoUtilsBase&); 
  
  virtual ~AliEMCALRecoUtilsBase() { ; }  
  
  // Track-cluster mathing extrapolation methods
  
  static Bool_t ExtrapolateTrackToEMCalSurface(AliVTrack *track, /*note, on success the call will change the track*/
                                               Double_t   emcalR=440, 
                                               Double_t   mass=0.1396,
                                               Double_t   step=20, 
                                               Double_t   minpT=0.35,
                                               Bool_t     useMassForTracking = kFALSE, 
                                               Bool_t     useDCA = kFALSE,
                                               Bool_t     useOuterParam = kFALSE);
 
  static Bool_t ExtrapolateTrackToEMCalSurface(AliExternalTrackParam *trkParam, 
                                               Double_t emcalR, 
                                               Double_t mass, 
                                               Double_t step, 
                                               Float_t &eta, 
                                               Float_t &phi, 
                                               Float_t &pt);
  
  static Bool_t ExtrapolateTrackToPosition(AliExternalTrackParam *trkParam, 
                                           const Float_t *clsPos, 
                                           Double_t mass, 
                                           Double_t step, 
                                           Float_t &tmpEta, 
                                           Float_t &tmpPhi);
  
  static Bool_t ExtrapolateTrackToCluster (AliExternalTrackParam *trkParam, 
                                           const AliVCluster *cluster, 
                                           Double_t mass, 
                                           Double_t step,
                                           Float_t &tmpEta, 
                                           Float_t &tmpPhi);
  
  /// Position depth enum list of possible particle types
  enum     ParticleType{kPhoton=0, kElectron=1, kHadron =2, kUnknown=-1};
  
  Float_t  GetDepth(Float_t eCluster, Int_t iParticle, Int_t iSM) const; 

  
  /// \cond CLASSIMP
  ClassDef(AliEMCALRecoUtilsBase, 1) ;
  /// \endcond

};

#endif // ALIEMCALRECOUTILSBASE_H


