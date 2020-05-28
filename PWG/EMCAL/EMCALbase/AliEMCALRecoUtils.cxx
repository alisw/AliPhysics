/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// ROOT includes
#include <TH2F.h>
#include <TArrayI.h>
#include <TArrayF.h>
#include <TObjArray.h>

// STEER includes
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliExternalTrackParam.h"
#include "AliESDfriendTrack.h"
#include "AliTrackerBase.h"
#include "AliMCEvent.h"

// EMCAL includes
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliTrackerBase.h"
#include "AliEMCALPIDUtils.h"

/// \cond CLASSIMP
ClassImp(AliEMCALRecoUtils) ;
/// \endcond

///
/// Constructor.
/// Initialize all constant values which have to be used
/// during Reco algorithm execution
///
//_____________________________________
AliEMCALRecoUtils::AliEMCALRecoUtils():
  fParticleType(0),                       fPosAlgo(0),                            
  fW0(0),                                 fShowerShapeCellLocationType(0),
  fNonLinearityFunction(0),               fNonLinearThreshold(0),                 fUseShaperNonlin(kFALSE),
  fSmearClusterEnergy(kFALSE),            fRandom(),
  fCellsRecalibrated(kFALSE),             fRecalibration(kFALSE),                 fUse1Drecalib(kFALSE),                  fEMCALRecalibrationFactors(),
  fCellsSingleChannelRecalibrated(kFALSE),fSingleChannelRecalibration(kFALSE),    fEMCALSingleChannelRecalibrationFactors(nullptr),
  fConstantTimeShift(0),                  fTimeRecalibration(kFALSE),             fEMCALTimeRecalibrationFactors(nullptr),       fLowGain(kFALSE),
  fUseL1PhaseInTimeRecalibration(kFALSE), fEMCALL1PhaseInTimeRecalibration(nullptr),
  fIsParRun(kFALSE),                      fCurrentParNumber(0),                   fGlobalEventID(),
  fDoUseMergedBC(kFALSE),
  fUseRunCorrectionFactors(kFALSE),       
  fRemoveBadChannels(kFALSE),             fRecalDistToBadChannels(kFALSE),        fEMCALBadChannelMap(nullptr),                  fUse1Dmap(kFALSE),
  fNCellsFromEMCALBorder(0),              fNoEMCALBorderAtEta0(kTRUE),
  fRejectExoticCluster(kFALSE),           fRejectExoticCells(kFALSE), 
  fExoticCellFraction(0),                 fExoticCellDiffTime(0),                 fExoticCellMinAmplitude(0),
  fPIDUtils(),                            fAODFilterMask(0),
  fAODHybridTracks(0),                    fAODTPCOnlyTracks(0),
  fMatchedTrackIndex(),                   fMatchedClusterIndex(), 
  fResidualEta(), fResidualPhi(),   fCutEtaPhiSum(kFALSE),                  fCutEtaPhiSeparate(kFALSE), 
  fCutR(0),                               fCutEta(0),                             fCutPhi(0),
  fClusterWindow(0),                      fMass(0),                           
  fStepSurface(0),                        fStepCluster(0),
  fITSTrackSA(kFALSE),                    fUseTrackDCA(kTRUE), // keep it active, but not working for old MC
  fUseOuterTrackParam(kFALSE),            fEMCalSurfaceDistance(440.),
  fTrackCutsType(0),                      fCutMinTrackPt(0),                      fCutMinNClusterTPC(0), 
  fCutMinNClusterITS(0),                  fCutMaxChi2PerClusterTPC(0),            fCutMaxChi2PerClusterITS(0),
  fCutRequireTPCRefit(kFALSE),            fCutRequireITSRefit(kFALSE),            fCutAcceptKinkDaughters(kFALSE),
  fCutMaxDCAToVertexXY(0),                fCutMaxDCAToVertexZ(0),                 fCutDCAToVertex2D(kFALSE),
  fCutRequireITSStandAlone(kFALSE),       fCutRequireITSpureSA(kFALSE),             
  fNMCGenerToAccept(0),                   fMCGenerToAcceptForTrack(1)
{
  // Init parameters
  InitParameters();
  
  for(Int_t j = 0; j <  5;    j++)  fMCGenerToAccept[j] =  "";
  
  fPIDUtils              = new AliEMCALPIDUtils();
  
  fBadStatusSelection[0] = kTRUE;
  fBadStatusSelection[1] = kTRUE;
  fBadStatusSelection[2] = kTRUE;
  fBadStatusSelection[3] = kTRUE;
}

//
// Copy constructor.
//
//______________________________________________________________________
AliEMCALRecoUtils::AliEMCALRecoUtils(const AliEMCALRecoUtils & reco) 
: AliEMCALRecoUtilsBase(reco), 
  fParticleType(reco.fParticleType),                         fPosAlgo(reco.fPosAlgo),     
  fW0(reco.fW0),                                             fShowerShapeCellLocationType(reco.fShowerShapeCellLocationType),
  fNonLinearityFunction(reco.fNonLinearityFunction),         fNonLinearThreshold(reco.fNonLinearThreshold),
  fUseShaperNonlin(reco.fUseShaperNonlin),
  fSmearClusterEnergy(reco.fSmearClusterEnergy),             fRandom(),
  fCellsRecalibrated(reco.fCellsRecalibrated),
  fRecalibration(reco.fRecalibration),                       fUse1Drecalib(reco.fUse1Drecalib),                   
  fEMCALRecalibrationFactors(NULL),
  fCellsSingleChannelRecalibrated(reco.fCellsRecalibrated),
  fSingleChannelRecalibration(reco.fRecalibration),          fEMCALSingleChannelRecalibrationFactors(NULL),
  fConstantTimeShift(reco.fConstantTimeShift),  
  fTimeRecalibration(reco.fTimeRecalibration),               fEMCALTimeRecalibrationFactors(NULL),
  fLowGain(reco.fLowGain),
  fUseL1PhaseInTimeRecalibration(reco.fUseL1PhaseInTimeRecalibration), 
  fEMCALL1PhaseInTimeRecalibration(reco.fEMCALL1PhaseInTimeRecalibration),
  fIsParRun(reco.fIsParRun),
  fCurrentParNumber(reco.fCurrentParNumber),
  fGlobalEventID(reco.fGlobalEventID),
  fDoUseMergedBC(reco.fDoUseMergedBC),
  fUseRunCorrectionFactors(reco.fUseRunCorrectionFactors),   
  fRemoveBadChannels(reco.fRemoveBadChannels),               fRecalDistToBadChannels(reco.fRecalDistToBadChannels),
  fEMCALBadChannelMap(NULL),                                 fUse1Dmap(reco.fUse1Dmap),
  fNCellsFromEMCALBorder(reco.fNCellsFromEMCALBorder),       fNoEMCALBorderAtEta0(reco.fNoEMCALBorderAtEta0),
  fRejectExoticCluster(reco.fRejectExoticCluster),           fRejectExoticCells(reco.fRejectExoticCells), 
  fExoticCellFraction(reco.fExoticCellFraction),             fExoticCellDiffTime(reco.fExoticCellDiffTime),               
  fExoticCellMinAmplitude(reco.fExoticCellMinAmplitude),
  fPIDUtils(reco.fPIDUtils),                                 fAODFilterMask(reco.fAODFilterMask),
  fAODHybridTracks(reco.fAODHybridTracks),                   fAODTPCOnlyTracks(reco.fAODTPCOnlyTracks),
  fMatchedTrackIndex(  reco.fMatchedTrackIndex),
  fMatchedClusterIndex(reco.fMatchedClusterIndex),
  fResidualEta(        reco.fResidualEta),
  fResidualPhi(        reco.fResidualPhi),
  fCutEtaPhiSum(reco.fCutEtaPhiSum),                         fCutEtaPhiSeparate(reco.fCutEtaPhiSeparate), 
  fCutR(reco.fCutR),        fCutEta(reco.fCutEta),           fCutPhi(reco.fCutPhi),
  fClusterWindow(reco.fClusterWindow),
  fMass(reco.fMass),        fStepSurface(reco.fStepSurface), fStepCluster(reco.fStepCluster),
  fITSTrackSA(reco.fITSTrackSA),                             fUseTrackDCA(reco.fUseTrackDCA),
  fUseOuterTrackParam(reco.fUseOuterTrackParam),             fEMCalSurfaceDistance(440.),
  fTrackCutsType(reco.fTrackCutsType),                       fCutMinTrackPt(reco.fCutMinTrackPt), 
  fCutMinNClusterTPC(reco.fCutMinNClusterTPC),               fCutMinNClusterITS(reco.fCutMinNClusterITS), 
  fCutMaxChi2PerClusterTPC(reco.fCutMaxChi2PerClusterTPC),   fCutMaxChi2PerClusterITS(reco.fCutMaxChi2PerClusterITS),
  fCutRequireTPCRefit(reco.fCutRequireTPCRefit),             fCutRequireITSRefit(reco.fCutRequireITSRefit),
  fCutAcceptKinkDaughters(reco.fCutAcceptKinkDaughters),     fCutMaxDCAToVertexXY(reco.fCutMaxDCAToVertexXY),    
  fCutMaxDCAToVertexZ(reco.fCutMaxDCAToVertexZ),             fCutDCAToVertex2D(reco.fCutDCAToVertex2D),
  fCutRequireITSStandAlone(reco.fCutRequireITSStandAlone),   fCutRequireITSpureSA(reco.fCutRequireITSpureSA),
  fNMCGenerToAccept(reco.fNMCGenerToAccept),                 fMCGenerToAcceptForTrack(reco.fMCGenerToAcceptForTrack)
{  
  for (Int_t i = 0; i < 15 ; i++) { fMisalRotShift[i]      = reco.fMisalRotShift[i]      ; 
                                    fMisalTransShift[i]    = reco.fMisalTransShift[i]    ; }
  for (Int_t i = 0; i < 10  ; i++){ fNonLinearityParams[i] = reco.fNonLinearityParams[i] ; }
  for (Int_t i = 0; i < 3  ; i++) { fSmearClusterParam[i]  = reco.fSmearClusterParam[i]  ; }
  for (Int_t j = 0; j < 5  ; j++) { fMCGenerToAccept[j]    = reco.fMCGenerToAccept[j]    ; }
  for (Int_t j = 0; j < 4  ; j++) { fBadStatusSelection[j] = reco.fBadStatusSelection[j] ; }

  if(reco.fEMCALBadChannelMap) {
    // Copy constructor - not taking ownership over calibration histograms
    fEMCALBadChannelMap = new TObjArray(reco.fEMCALBadChannelMap->GetEntries());
    fEMCALBadChannelMap->SetOwner(false);
    for(int ism = 0; ism < reco.fEMCALBadChannelMap->GetEntries(); ism++) fEMCALBadChannelMap->AddAt(reco.fEMCALBadChannelMap->At(ism), ism);
  }

  if(reco.fEMCALRecalibrationFactors) {
    // Copy constructor - not taking ownership over calibration histograms
    fEMCALRecalibrationFactors = new TObjArray(reco.fEMCALRecalibrationFactors->GetEntries());
    fEMCALRecalibrationFactors->SetOwner(false);
    for(int ism = 0; ism < reco.fEMCALRecalibrationFactors->GetEntries(); ism++) fEMCALRecalibrationFactors->AddAt(reco.fEMCALRecalibrationFactors->At(ism), ism);
  }

  if(reco.fEMCALSingleChannelRecalibrationFactors) {
    // Copy constructor - not taking ownership over calibration histograms
    fEMCALSingleChannelRecalibrationFactors = new TObjArray(reco.fEMCALSingleChannelRecalibrationFactors->GetEntries());
    fEMCALSingleChannelRecalibrationFactors->SetOwner(false);
    for(int ism = 0; ism < reco.fEMCALSingleChannelRecalibrationFactors->GetEntries(); ism++) fEMCALSingleChannelRecalibrationFactors->AddAt(reco.fEMCALSingleChannelRecalibrationFactors->At(ism), ism);
  }

  
  if(reco.fEMCALTimeRecalibrationFactors) {
    // Copy constructor - not taking ownership over calibration histograms
    fEMCALTimeRecalibrationFactors = new TObjArray(reco.fEMCALTimeRecalibrationFactors->GetEntries());
    fEMCALTimeRecalibrationFactors->SetOwner(false);
    for(int ism = 0; ism < reco.fEMCALTimeRecalibrationFactors->GetEntries(); ism++) fEMCALTimeRecalibrationFactors->AddAt(reco.fEMCALTimeRecalibrationFactors->At(ism), ism);
  }
}

///
/// Assignment operator.
///
//______________________________________________________________________
AliEMCALRecoUtils & AliEMCALRecoUtils::operator = (const AliEMCALRecoUtils & reco) 
{  
  if (this == &reco)return *this;
  ((TNamed *)this)->operator=(reco);
  
  for (Int_t i = 0; i < 15 ; i++) { fMisalTransShift[i]    = reco.fMisalTransShift[i]    ; 
    fMisalRotShift[i]      = reco.fMisalRotShift[i]      ; }
  for (Int_t i = 0; i < 10  ; i++) { fNonLinearityParams[i] = reco.fNonLinearityParams[i] ; }
  for (Int_t i = 0; i < 3  ; i++) { fSmearClusterParam[i]  = reco.fSmearClusterParam[i]  ; }   
  
  fParticleType              = reco.fParticleType;
  fPosAlgo                   = reco.fPosAlgo; 
  fW0                        = reco.fW0;
  fShowerShapeCellLocationType = reco.fShowerShapeCellLocationType;
  
  fNonLinearityFunction      = reco.fNonLinearityFunction;
  fNonLinearThreshold        = reco.fNonLinearThreshold;
  fUseShaperNonlin           = reco.fUseShaperNonlin;
  fSmearClusterEnergy        = reco.fSmearClusterEnergy;
  
  fCellsRecalibrated         = reco.fCellsRecalibrated;
  fRecalibration             = reco.fRecalibration;
  fUse1Drecalib              = reco.fUse1Drecalib;
  
  fConstantTimeShift         = reco.fConstantTimeShift;
  fTimeRecalibration         = reco.fTimeRecalibration;
  fLowGain                   = reco.fLowGain;

  fUseL1PhaseInTimeRecalibration   = reco.fUseL1PhaseInTimeRecalibration;
  fEMCALL1PhaseInTimeRecalibration = reco.fEMCALL1PhaseInTimeRecalibration;

  fIsParRun                  = reco.fIsParRun;
  fCurrentParNumber          = reco.fCurrentParNumber;
  fGlobalEventID             = reco.fGlobalEventID;
  
  fDoUseMergedBC             = reco.fDoUseMergedBC;
  
  fUseRunCorrectionFactors   = reco.fUseRunCorrectionFactors;
  
  fRemoveBadChannels         = reco.fRemoveBadChannels;
  fRecalDistToBadChannels    = reco.fRecalDistToBadChannels;
  fUse1Dmap                  = reco.fUse1Dmap;
  
  fNCellsFromEMCALBorder     = reco.fNCellsFromEMCALBorder;
  fNoEMCALBorderAtEta0       = reco.fNoEMCALBorderAtEta0;
  
  fRejectExoticCluster       = reco.fRejectExoticCluster;           
  fRejectExoticCells         = reco.fRejectExoticCells; 
  fExoticCellFraction        = reco.fExoticCellFraction;
  fExoticCellDiffTime        = reco.fExoticCellDiffTime;              
  fExoticCellMinAmplitude    = reco.fExoticCellMinAmplitude;
  
  fPIDUtils                  = reco.fPIDUtils;
  
  fAODFilterMask             = reco.fAODFilterMask;
  fAODHybridTracks           = reco.fAODHybridTracks;
  fAODTPCOnlyTracks          = reco.fAODTPCOnlyTracks;
  
  fCutEtaPhiSum              = reco.fCutEtaPhiSum;
  fCutEtaPhiSeparate         = reco.fCutEtaPhiSeparate;
  fCutR                      = reco.fCutR;
  fCutEta                    = reco.fCutEta;
  fCutPhi                    = reco.fCutPhi;
  fClusterWindow             = reco.fClusterWindow;
  fMass                      = reco.fMass;
  fStepSurface               = reco.fStepSurface;
  fStepCluster               = reco.fStepCluster;
  fITSTrackSA                = reco.fITSTrackSA;
  fUseTrackDCA               = reco.fUseTrackDCA;
  fUseOuterTrackParam        = reco.fUseOuterTrackParam;
  fEMCalSurfaceDistance      = reco.fEMCalSurfaceDistance;
  
  fTrackCutsType             = reco.fTrackCutsType;
  fCutMinTrackPt             = reco.fCutMinTrackPt;
  fCutMinNClusterTPC         = reco.fCutMinNClusterTPC;
  fCutMinNClusterITS         = reco.fCutMinNClusterITS; 
  fCutMaxChi2PerClusterTPC   = reco.fCutMaxChi2PerClusterTPC;
  fCutMaxChi2PerClusterITS   = reco.fCutMaxChi2PerClusterITS;
  fCutRequireTPCRefit        = reco.fCutRequireTPCRefit;
  fCutRequireITSRefit        = reco.fCutRequireITSRefit;
  fCutAcceptKinkDaughters    = reco.fCutAcceptKinkDaughters;
  fCutMaxDCAToVertexXY       = reco.fCutMaxDCAToVertexXY;
  fCutMaxDCAToVertexZ        = reco.fCutMaxDCAToVertexZ;
  fCutDCAToVertex2D          = reco.fCutDCAToVertex2D;
  fCutRequireITSStandAlone   = reco.fCutRequireITSStandAlone; 
  fCutRequireITSpureSA       = reco.fCutRequireITSpureSA;

  fNMCGenerToAccept          = reco.fNMCGenerToAccept;
  fMCGenerToAcceptForTrack   = reco.fMCGenerToAcceptForTrack;
  for (Int_t j = 0; j < 5  ; j++)  
    fMCGenerToAccept[j]     = reco.fMCGenerToAccept[j];

  //
  // Assign or copy construct the different TArrays
  //
  fResidualEta = reco.fResidualEta;
  fResidualPhi = reco.fResidualPhi;
  fMatchedTrackIndex = reco.fMatchedTrackIndex;
  fMatchedClusterIndex = reco.fMatchedClusterIndex;

  for (Int_t j = 0; j < 4  ; j++) 
   fBadStatusSelection[j] = reco.fBadStatusSelection[j] ; 
  
  if(fEMCALBadChannelMap) delete fEMCALBadChannelMap;
  if(reco.fEMCALBadChannelMap) {
    // Copy constructor - not taking ownership over calibration histograms
    fEMCALBadChannelMap = new TObjArray(reco.fEMCALBadChannelMap->GetEntries());
    fEMCALBadChannelMap->SetOwner(false);
    for(int ism = 0; ism < reco.fEMCALBadChannelMap->GetEntries(); ism++) fEMCALBadChannelMap->AddAt(reco.fEMCALBadChannelMap->At(ism), ism);
  }

  if(fEMCALRecalibrationFactors) delete fEMCALRecalibrationFactors;
  if(reco.fEMCALRecalibrationFactors) {
    // Copy constructor - not taking ownership over calibration histograms
    fEMCALRecalibrationFactors = new TObjArray(reco.fEMCALRecalibrationFactors->GetEntries());
    fEMCALRecalibrationFactors->SetOwner(false);
    for(int ism = 0; ism < reco.fEMCALRecalibrationFactors->GetEntries(); ism++) fEMCALRecalibrationFactors->AddAt(reco.fEMCALRecalibrationFactors->At(ism), ism);
  }

  if(fEMCALSingleChannelRecalibrationFactors) delete fEMCALSingleChannelRecalibrationFactors;
  if(reco.fEMCALSingleChannelRecalibrationFactors) {
    // Copy constructor - not taking ownership over calibration histograms
    fEMCALSingleChannelRecalibrationFactors = new TObjArray(reco.fEMCALSingleChannelRecalibrationFactors->GetEntries());
    fEMCALSingleChannelRecalibrationFactors->SetOwner(false);
    for(int ism = 0; ism < reco.fEMCALSingleChannelRecalibrationFactors->GetEntries(); ism++) fEMCALSingleChannelRecalibrationFactors->AddAt(reco.fEMCALSingleChannelRecalibrationFactors->At(ism), ism);
  }
  
  if(fEMCALTimeRecalibrationFactors) delete fEMCALTimeRecalibrationFactors;
  if(reco.fEMCALTimeRecalibrationFactors) {
    // Copy constructor - not taking ownership over calibration histograms
    fEMCALTimeRecalibrationFactors = new TObjArray(reco.fEMCALTimeRecalibrationFactors->GetEntries());
    fEMCALTimeRecalibrationFactors->SetOwner(false);
    for(int ism = 0; ism < reco.fEMCALTimeRecalibrationFactors->GetEntries(); ism++) fEMCALTimeRecalibrationFactors->AddAt(reco.fEMCALTimeRecalibrationFactors->At(ism), ism);
  }

  if(this->fEMCALL1PhaseInTimeRecalibration) delete fEMCALL1PhaseInTimeRecalibration;
  if(reco.fEMCALL1PhaseInTimeRecalibration) {
    // Copy constructor - not taking ownership over calibration histograms
    fEMCALL1PhaseInTimeRecalibration = new TObjArray(reco.fEMCALL1PhaseInTimeRecalibration->GetEntries());
    fEMCALL1PhaseInTimeRecalibration->SetOwner(false);
    for(int ism = 0; ism < reco.fEMCALL1PhaseInTimeRecalibration->GetEntries(); ism++) fEMCALL1PhaseInTimeRecalibration->AddAt(reco.fEMCALL1PhaseInTimeRecalibration->At(ism), ism);
  }

  return *this;
}

///
/// Destructor.
///
//_____________________________________
AliEMCALRecoUtils::~AliEMCALRecoUtils()
{  
  if (fEMCALRecalibrationFactors) 
  { 
    delete fEMCALRecalibrationFactors;
  }  

  if (fEMCALSingleChannelRecalibrationFactors) 
  { 
    delete fEMCALSingleChannelRecalibrationFactors;
  }  

  if (fEMCALTimeRecalibrationFactors) 
  { 
    delete fEMCALTimeRecalibrationFactors;
  }  

  if(fEMCALL1PhaseInTimeRecalibration) 
  {  
    delete fEMCALL1PhaseInTimeRecalibration;
  }

  if (fEMCALBadChannelMap) 
  { 
    delete fEMCALBadChannelMap;
  }
 
  delete fPIDUtils            ;

  InitTrackCuts();
}

///
/// Reject cell if acceptance criteria not passed (correct cell number, is it bad channel) 
/// and calibrate it in energy and time.
///
/// \param absID: absolute cell ID number
/// \param bc: bunch crossing number
/// \param amp: input cell energy amplitude, output calibrated amplitude
/// \param time: input cell time, output calibrated time
/// \param cells: list of cells
///
/// \return bool quality of cell, exists or not 
///
//_______________________________________________________________________________
Bool_t AliEMCALRecoUtils::AcceptCalibrateCell(Int_t absID, Int_t bc,
                                              Float_t  & amp,    Double_t & time, 
                                              AliVCaloCells* cells) 
{  
  AliEMCALGeometry* geom = AliEMCALGeometry::GetInstance();
  
  if(!geom)
  {
    AliError("No instance of the geometry is available");
    return kFALSE;
  }

  if ( absID < 0 || absID >= 24*48*geom->GetNumberOfSuperModules() ) 
    return kFALSE;
  
  Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1, status=0; 
  
  if (!geom->GetCellIndex(absID,imod,iTower,iIphi,iIeta)) 
  {
    // cell absID does not exist
    amp=0; time = 1.e9;
    return kFALSE; 
  }
  
  geom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);  

  // Do not include bad channels found in analysis,
  if ( IsBadChannelsRemovalSwitchedOn() )
  {
    Bool_t bad = kFALSE;

    if(fUse1Dmap)
      bad = GetEMCALChannelStatus1D(absID,status);
    else
      bad = GetEMCALChannelStatus(imod, ieta, iphi,status);
    
    if ( status > 0 )
      AliDebug(1,Form("Channel absId %d, status %d, set as bad %d",absID, status, bad));
    
    if ( bad ) return kFALSE;
  }
  Bool_t isLowGain = !(cells->GetCellHighGain(absID));//HG = false -> LG = true
  
  //Recalibrate energy
  amp  = cells->GetCellAmplitude(absID);
  if (!fCellsRecalibrated && IsRecalibrationOn()){
    if(fUse1Drecalib)
      amp *= GetEMCALChannelRecalibrationFactor1D(absID);
    else
      amp *= GetEMCALChannelRecalibrationFactor(imod,ieta,iphi);

    if(fUseShaperNonlin && isLowGain){
      amp = CorrectShaperNonLin(amp,fUse1Drecalib ? GetEMCALChannelRecalibrationFactor1D(absID) : GetEMCALChannelRecalibrationFactor(imod,ieta,iphi));
    }
  }
  // Recalibrate time
  time = cells->GetCellTime(absID);
  time-=fConstantTimeShift*1e-9; // only in case of old Run1 simulation

  RecalibrateCellTime(absID,bc,time,isLowGain);
  
  //Recalibrate time with L1 phase 
  RecalibrateCellTimeL1Phase(imod, bc, time, fCurrentParNumber);

  return kTRUE;
}

///
/// Given the list of AbsId cells of the cluster, get the maximum cell and 
/// check if there are fNCellsFromBorder from the calorimeter border.
///
/// \param geom: AliEMCALGeometry pointer
/// \param cluster: AliVCluster pointer
/// \param cells: list of cells
///
/// \return bool, true if cluster in defined fiducial region
//_____________________________________________________________________________
Bool_t AliEMCALRecoUtils::CheckCellFiducialRegion(const AliEMCALGeometry* geom, 
                                                  const AliVCluster* cluster, 
                                                  AliVCaloCells* cells) 
{  
  if (!cluster)
  {
    AliInfo("Cluster pointer null!");
    return kFALSE;
  }
  
  // If the distance to the border is 0 or negative just exit accept all clusters
  if (cells->GetType()==AliVCaloCells::kEMCALCell && fNCellsFromEMCALBorder <= 0 ) 
    return kTRUE;
  
  Int_t absIdMax  = -1, iSM =-1, ieta = -1, iphi = -1;
  Bool_t shared = kFALSE;
  GetMaxEnergyCell(geom, cells, cluster, absIdMax,  iSM, ieta, iphi, shared);
  
  AliDebug(2,Form("Cluster Max AbsId %d, Cell Energy %2.2f, Cluster Energy %2.2f, Ncells from border %d, EMCAL eta=0 %d\n", 
                  absIdMax, cells->GetCellAmplitude(absIdMax), cluster->E(), fNCellsFromEMCALBorder, fNoEMCALBorderAtEta0));
  
  if (absIdMax==-1) return kFALSE;
  
  // Check if the cell is close to the borders:
  Bool_t okrow = kFALSE;
  Bool_t okcol = kFALSE;
  
  if (iSM < 0 || iphi < 0 || ieta < 0 ) 
  {
    AliFatal(Form("Negative value for super module: %d, or cell ieta: %d, or cell iphi: %d, check EMCAL geometry name\n",
                  iSM,ieta,iphi));
    return kFALSE; // trick coverity
  }
  
  // Check rows/phi
  Int_t iPhiLast = 24;
   if      ( geom->GetSMType(iSM) == AliEMCALGeometry::kEMCAL_Half ) iPhiLast /= 2;
   else if ( geom->GetSMType(iSM) == AliEMCALGeometry::kEMCAL_3rd  ) iPhiLast /= 3;// 1/3 sm case
   else if ( geom->GetSMType(iSM) == AliEMCALGeometry::kDCAL_Ext   ) iPhiLast /= 3;// 1/3 sm case

  if(iphi >= fNCellsFromEMCALBorder && iphi < iPhiLast - fNCellsFromEMCALBorder) okrow = kTRUE; 

  // Check columns/eta
  Int_t iEtaLast = 48;
  if ( !fNoEMCALBorderAtEta0 || geom->GetSMType(iSM) == AliEMCALGeometry::kDCAL_Standard ) 
  {
    // consider inner border
    if ( geom->IsDCALSM(iSM) ) iEtaLast = iEtaLast*2/3;        
  
    if ( ieta  > fNCellsFromEMCALBorder && ieta < iEtaLast-fNCellsFromEMCALBorder ) okcol = kTRUE;  
  } 
  else 
  {
    if (iSM%2==0) 
    {
     if (ieta >= fNCellsFromEMCALBorder)           okcol = kTRUE;  
    } 
    else 
    {
     if (ieta <  iEtaLast-fNCellsFromEMCALBorder)  okcol = kTRUE; 
    }
  }//eta 0 not checked
  
  AliDebug(2,Form("EMCAL Cluster in %d cells fiducial volume: ieta %d, iphi %d, SM %d:  column? %d, row? %d\nq",
                  fNCellsFromEMCALBorder, ieta, iphi, iSM, okcol, okrow));
  
  if (okcol && okrow) 
  {
    return kTRUE;
  } 
  else
  {
    AliDebug(2,Form("Reject cluster in border, max cell : ieta %d, iphi %d, SM %d\n",ieta, iphi, iSM));
    
    return kFALSE;
  }
}  

///
/// Check that in the cluster cells, there is no bad channel of those stored 
/// in fEMCALBadChannelMap 
///
/// \param geom: AliEMCALGeometry pointer
/// \param cellList: list of cells absolute ID in cluster
/// \param nCells: number of cells in cluster
///
/// \return bool, true if cluster contains a bad channel
///
//_______________________________________________________________________________
Bool_t AliEMCALRecoUtils::ClusterContainsBadChannel(const AliEMCALGeometry* geom, 
                                                    const UShort_t* cellList, 
                                                    Int_t nCells)
{  
  if (!fRemoveBadChannels)  return kFALSE;
  if (!fEMCALBadChannelMap) return kFALSE;
  
  Int_t icol = -1;
  Int_t irow = -1;
  Int_t imod = -1;
  for (Int_t iCell = 0; iCell<nCells; iCell++) 
  {
    //Get the column and row
    Int_t iTower = -1, iIphi = -1, iIeta = -1; 
    geom->GetCellIndex(cellList[iCell],imod,iTower,iIphi,iIeta); 
  
    if (fEMCALBadChannelMap->GetEntries() <= imod) continue;
    
    geom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,irow,icol);      
    
    Int_t status = 0;

    if(fUse1Dmap){
      if (GetEMCALChannelStatus1D(cellList[iCell], status)) 
      {
        AliDebug(2,Form("Cluster with bad channel: ID %d, status %d\n",cellList[iCell], status));
        return kTRUE;
      }
    }else{
      if (GetEMCALChannelStatus(imod, icol, irow, status)) 
      {
        AliDebug(2,Form("Cluster with bad channel: SM %d, col %d, row %d, status %d\n",imod, icol, irow, status));
        return kTRUE;
      }
    }
  }// cell cluster loop

  return kFALSE;
}

///
/// Calculate the energy in the cross around the energy of a given cell.
/// Used in exotic clusters/cells rejection.
///
/// \param absID: controlled cell absolute ID number
/// \param tcell: time of cell under control
/// \param cells: full list of cells
/// \param bc: bunch crossing number
/// \param cellMinEn: add the cell energy if large enough (for high energy clusters)
/// \param useWeight: add the cell energy if w > 0
/// \param energy: cluster or cell max energy, used for weight calculation
///
/// \return float E_cross
///
//___________________________________________________________________________
Float_t AliEMCALRecoUtils::GetECross(Int_t absID, Double_t tcell,
                                     AliVCaloCells* cells, Int_t bc, 
                                     Float_t cellMinEn, Bool_t useWeight, Float_t energy )
{  
  AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance();
  
  if(!geom)
  {
    AliError("No instance of the geometry is available");
    return -1;
  }

  Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
  geom->GetCellIndex(absID,imod,iTower,iIphi,iIeta); 
  geom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);  
  
  // Get close cells index, energy and time, not in corners
  
  Int_t absID1 = -1;
  Int_t absID2 = -1;
  
  if ( iphi < AliEMCALGeoParams::fgkEMCALRows-1) absID1 = geom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
  if ( iphi > 0 )                                absID2 = geom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);
  
  // In case of cell in eta = 0 border, depending on SM shift the cross cell index
  
  Int_t absID3 = -1;
  Int_t absID4 = -1;
  
  if ( ieta == AliEMCALGeoParams::fgkEMCALCols-1 && !(imod%2) ) 
  {
    absID3 = geom-> GetAbsCellIdFromCellIndexes(imod+1, iphi, 0);
    absID4 = geom-> GetAbsCellIdFromCellIndexes(imod,   iphi, ieta-1); 
  } 
  else if ( ieta == 0 && imod%2 ) 
  {
    absID3 = geom-> GetAbsCellIdFromCellIndexes(imod,   iphi, ieta+1);
    absID4 = geom-> GetAbsCellIdFromCellIndexes(imod-1, iphi, AliEMCALGeoParams::fgkEMCALCols-1); 
  } 
  else 
  {
    if ( ieta < AliEMCALGeoParams::fgkEMCALCols-1 ) 
      absID3 = geom-> GetAbsCellIdFromCellIndexes(imod, iphi, ieta+1);
    if ( ieta > 0 )                                 
      absID4 = geom-> GetAbsCellIdFromCellIndexes(imod, iphi, ieta-1); 
  }
  
  //printf("IMOD %d, AbsId %d, a %d, b %d, c %d e %d \n",imod,absID,absID1,absID2,absID3,absID4);
  
  Float_t  ecell1  = 0, ecell2  = 0, ecell3  = 0, ecell4  = 0;
  Double_t tcell1  = 0, tcell2  = 0, tcell3  = 0, tcell4  = 0;

  AcceptCalibrateCell(absID1,bc, ecell1,tcell1,cells); 
  AcceptCalibrateCell(absID2,bc, ecell2,tcell2,cells); 
  AcceptCalibrateCell(absID3,bc, ecell3,tcell3,cells); 
  AcceptCalibrateCell(absID4,bc, ecell4,tcell4,cells); 
  
  if (TMath::Abs(tcell-tcell1)*1.e9 > fExoticCellDiffTime) ecell1 = 0 ;
  if (TMath::Abs(tcell-tcell2)*1.e9 > fExoticCellDiffTime) ecell2 = 0 ;
  if (TMath::Abs(tcell-tcell3)*1.e9 > fExoticCellDiffTime) ecell3 = 0 ;
  if (TMath::Abs(tcell-tcell4)*1.e9 > fExoticCellDiffTime) ecell4 = 0 ;
 
  Float_t w1 = 1, w2 = 1, w3 = 1, w4 = 1;
  if ( useWeight )
  {
    w1 = GetCellWeight(ecell1,energy);
    w2 = GetCellWeight(ecell2,energy);
    w3 = GetCellWeight(ecell3,energy);
    w4 = GetCellWeight(ecell4,energy);
  }
  
  if ( ecell1 < cellMinEn || w1 <= 0 ) ecell1 = 0 ;
  if ( ecell2 < cellMinEn || w2 <= 0 ) ecell2 = 0 ;
  if ( ecell3 < cellMinEn || w3 <= 0 ) ecell3 = 0 ;
  if ( ecell4 < cellMinEn || w4 <= 0 ) ecell4 = 0 ;
 
  return ecell1+ecell2+ecell3+ecell4;
}

//________________________________________________________________________________________
/// Check if 2 cells belong to the same T-Card
/// Only for EMCal.
///
///  \param absId1: Reference absId cell
///  \param absId2: Cross checked cell absId
///  \param rowDiff: Distance in rows
///  \param colDiff: Distance in columns
///  \return true if belong to same TCard
///
//________________________________________________________________________________________
Bool_t  AliEMCALRecoUtils::IsAbsIDsFromTCard(Int_t absId1, Int_t absId2, 
                                             Int_t & rowDiff, Int_t & colDiff) const
{  
  AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance();
  
  if ( !geom )
  {
    AliError("No instance of the geometry is available");
    return -1;
  }
  
  rowDiff = -100;
  colDiff = -100;
  
  if(absId1 == absId2) return kFALSE;
  
  // Check if in same SM, if not for sure not same TCard
  Int_t sm1 = geom->GetSuperModuleNumber(absId1);
  Int_t sm2 = geom->GetSuperModuleNumber(absId2);
  if ( sm1 != sm2 ) return kFALSE ;
  
  // Get the column and row of each absId
  Int_t iTower = -1, iIphi = -1, iIeta = -1;
  
  Int_t col1, row1;
  geom->GetCellIndex(absId1,sm1,iTower,iIphi,iIeta);
  geom->GetCellPhiEtaIndexInSModule(sm1,iTower,iIphi, iIeta,row1,col1);
  
  Int_t col2, row2;
  geom->GetCellIndex(absId2,sm2,iTower,iIphi,iIeta);
  geom->GetCellPhiEtaIndexInSModule(sm2,iTower,iIphi, iIeta,row2,col2);
  
  Int_t row0 = Int_t(row1-row1%8);
  Int_t col0 = Int_t(col1-col1%2);
  
  Int_t rowDiff0 = row2-row0;
  Int_t colDiff0 = col2-col0;
  
  rowDiff = row1-row2;
  colDiff = col1-col2;
  
  // TCard is made by 2x8 towers
  if ( colDiff0 >=0 && colDiff0 < 2 && rowDiff0 >=0 && rowDiff0 < 8 ) 
  {
    
    //    printf("\t absId (%d,%d), sm %d; col (%d,%d), colDiff %d; row (%d,%d),rowDiff %d\n",
    //           absId1 , absId2, sm1, 
    //           col1, col2, colDiff, 
    //           row1, row2, rowDiff);
    return kTRUE ;
  }
  else
    return kFALSE;
}

//________________________________________________________________________________________
/// Count the number of cells in same or different T-Card
/// as the highest energy cell in the cluster.
/// Only for EMCal.
///
///  \param clus: AliVCluster
///  \param absIdMax: Abs Id number of highest energy cell in the cluster
///  \param cells: AliVCaloCells
///  \param nDiff: number of cells in different T-Card
///  \param nSame: number of cells in same T-Card
///  \param eDiff: sum of energy of cells in different T-Card
///  \param eSame: sum of energy of cells in same T-Card
///  \param emin: apply a min energy cut on cells while counting
///
//________________________________________________________________________________________
void AliEMCALRecoUtils::GetEnergyAndNumberOfCellsInTCard
(AliVCluster* clus, Int_t absIdMax, AliVCaloCells* cells,
 Int_t   & nDiff, Int_t   & nSame, 
 Float_t & eDiff, Float_t & eSame, 
 Float_t   emin)
{
  Int_t nCaloCellsPerCluster = clus->GetNCells();
  
  nDiff = 0;
  nSame = 0;
  Int_t   absId   = -1;
  Float_t amp     = 0;
  Int_t   rowDiff = -100, colDiff = -100;
  
  // Loop on cluster cells count those in same or different T-Card
  // with respect highest energy cell.
  for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) 
  {
    absId = clus->GetCellsAbsId()[ipos];  
    
    amp   = cells->GetCellAmplitude(absId);
    
    if ( absId == absIdMax || amp < emin ) continue;
    
    if ( IsAbsIDsFromTCard(absIdMax,absId,rowDiff,colDiff) ) 
    { 
      nSame++; 
      eSame+=amp; 
    }
    else             
    { 
      nDiff++; 
      eDiff+= amp; 
    }
    
  } // cell cluster loop
  
}

///
/// Look to cell neighbourhood and reject if it seems exotic
/// Do before recalibrating the cells.
///
/// \param absID: controlled cell absolute ID number
/// \param cells: full list of cells
/// \param bc: bunch crossing number
///
/// \return bool true if cell is found exotic
///
//_____________________________________________________________________________________________
Bool_t AliEMCALRecoUtils::IsExoticCell(Int_t absID, AliVCaloCells* cells, Int_t bc)
{
  if (!fRejectExoticCells) return kFALSE;
  
  Float_t  ecell  = 0;
  Double_t tcell  = 0;
  Bool_t   accept = AcceptCalibrateCell(absID, bc, ecell ,tcell ,cells); 
  
  if (!accept) return kTRUE; // reject this cell
  
  if (ecell < fExoticCellMinAmplitude) return kFALSE; // do not reject low energy cells

  Float_t eCross = GetECross(absID,tcell,cells,bc);
  
  if (1-eCross/ecell > fExoticCellFraction) 
  {
    AliDebug(2,Form("AliEMCALRecoUtils::IsExoticCell() - EXOTIC CELL id %d, eCell %f, eCross %f, 1-eCross/eCell %f\n",
                    absID,ecell,eCross,1-eCross/ecell));
    return kTRUE;
  }

  return kFALSE;
}

///
/// Check if the cluster highest energy tower is exotic.
///
/// \param cluster: pointer to AliVCluster
/// \param cells: full list of cells
/// \param bc: bunch crossing number
///
/// \param cluster:  pointer to AliVCluster
///
//___________________________________________________________________
Bool_t AliEMCALRecoUtils::IsExoticCluster(const AliVCluster *cluster, 
                                          AliVCaloCells *cells, 
                                          Int_t bc) 
{  
  if (!cluster) 
  {
    AliInfo("Cluster pointer null!");
    return kFALSE;
  }
  
  if (!fRejectExoticCluster) return kFALSE;
  
  // Get highest energy tower
  AliEMCALGeometry* geom = AliEMCALGeometry::GetInstance();
  
  if(!geom)
  {
    AliError("No instance of the geometry is available");
    return kFALSE;
  }

  Int_t iSupMod = -1, absId = -1, ieta = -1, iphi = -1;
  Bool_t shared = kFALSE;
  GetMaxEnergyCell(geom, cells, cluster, absId, iSupMod, ieta, iphi, shared);
  
  return IsExoticCell(absId,cells,bc);
}

///
/// In case of MC analysis, smear energy to match resolution/calibration in real data 
/// (old, in principle not needed anymore).
///
/// \param cluster: pointer to AliVCluster
///
/// \return float with smeared energy
///
//_______________________________________________________________________
Float_t AliEMCALRecoUtils::SmearClusterEnergy(const AliVCluster* cluster) 
{  
  if (!cluster) 
  {
    AliInfo("Cluster pointer null!");
    return 0;
  }
  
  Float_t energy    = cluster->E() ;
  Float_t rdmEnergy = energy ;
  
  if (fSmearClusterEnergy) 
  {
    rdmEnergy = fRandom.Gaus(energy,fSmearClusterParam[0] * TMath::Sqrt(energy) +
                                    fSmearClusterParam[1] * energy +
                                    fSmearClusterParam[2] );
    AliDebug(2, Form("Energy: original %f, smeared %f\n", energy, rdmEnergy));
  }
  
  return rdmEnergy;
}

///
/// Correct cluster energy from non linearity functions, defined in enum NonlinearityFunctions
/// Recomended for data kBeamTestCorrectedv3 and for simulation kPi0MCv3
///
/// \param cluster: pointer to AliVCluster
///
/// \return float with corrected cluster energy
///
//____________________________________________________________________________
Float_t AliEMCALRecoUtils::CorrectClusterEnergyLinearity(AliVCluster* cluster)
{  
  if (!cluster) 
  {
    AliInfo("Cluster pointer null!");
    return 0;
  }
  
  Float_t energy = cluster->E();
  if (energy < 0.100) 
  {
    // Clusters with less than 100 MeV or negative are not possible
    AliInfo(Form("Too low cluster energy!, E = %f < 0.100 GeV",energy));
    return 0;
  }
  else if(energy > 300.)
  {
    AliInfo(Form("Too high cluster energy!, E = %f GeV, do not apply non linearity",energy));
    return energy;
  }
  
  switch (fNonLinearityFunction) 
  {
    case kPi0MC:
    {
      //Non-Linearity correction (from MC with function ([0]*exp(-[1]/E))+(([2]/([3]*2.*TMath::Pi())*exp(-(E-[4])^2/(2.*[3]^2)))))
      //fNonLinearityParams[0] = 1.014;
      //fNonLinearityParams[1] =-0.03329;
      //fNonLinearityParams[2] =-0.3853;
      //fNonLinearityParams[3] = 0.5423;
      //fNonLinearityParams[4] =-0.4335;
       energy *= (fNonLinearityParams[0]*exp(-fNonLinearityParams[1]/energy))+
                  ((fNonLinearityParams[2]/(fNonLinearityParams[3]*2.*TMath::Pi())*
                    exp(-(energy-fNonLinearityParams[4])*(energy-fNonLinearityParams[4])/(2.*fNonLinearityParams[3]*fNonLinearityParams[3]))));
      break;
    }
     
    case kPi0MCv2:
    {
      //Non-Linearity correction (from MC with function [0]/((x+[1])^[2]))+1;
      //fNonLinearityParams[0] = 3.11111e-02;
      //fNonLinearityParams[1] =-5.71666e-02; 
      //fNonLinearityParams[2] = 5.67995e-01;      
      
      energy *= fNonLinearityParams[0]/TMath::Power(energy+fNonLinearityParams[1],fNonLinearityParams[2])+1;
      break;
    }
    
    case kPi0MCv3:
    {
      //Same as beam test corrected, change parameters
      //fNonLinearityParams[0] =  9.81039e-01
      //fNonLinearityParams[1] =  1.13508e-01;
      //fNonLinearityParams[2] =  1.00173e+00; 
      //fNonLinearityParams[3] =  9.67998e-02;
      //fNonLinearityParams[4] =  2.19381e+02;
      //fNonLinearityParams[5] =  6.31604e+01;
      //fNonLinearityParams[6] =  1;
      energy *= fNonLinearityParams[6]/(fNonLinearityParams[0]*(1./(1.+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]))*1./(1.+fNonLinearityParams[3]*exp((energy-fNonLinearityParams[4])/fNonLinearityParams[5]))));
      
      break;
    }
      
      
    case kPi0GammaGamma:
    {
      //Non-Linearity correction (from Olga Data with function p0+p1*exp(-p2*E))
      //fNonLinearityParams[0] = 1.04;
      //fNonLinearityParams[1] = -0.1445;
      //fNonLinearityParams[2] = 1.046;
      energy /= (fNonLinearityParams[0]+fNonLinearityParams[1]*exp(-fNonLinearityParams[2]*energy)); //Olga function
      break;
    }
      
    case kPi0GammaConversion:
    {
      //Non-Linearity correction (Nicolas from Dimitri Data with function C*[1-a*exp(-b*E)])
      //fNonLinearityParams[0] = 0.139393/0.1349766;
      //fNonLinearityParams[1] = 0.0566186;
      //fNonLinearityParams[2] = 0.982133;
      energy /= fNonLinearityParams[0]*(1-fNonLinearityParams[1]*exp(-fNonLinearityParams[2]*energy));
      
      break;
    }
      
    case kBeamTest:
    {
      //From beam test, Alexei's results, for different ZS thresholds
      //                        th=30 MeV; th = 45 MeV; th = 75 MeV
      //fNonLinearityParams[0] = 1.007;      1.003;      1.002 
      //fNonLinearityParams[1] = 0.894;      0.719;      0.797 
      //fNonLinearityParams[2] = 0.246;      0.334;      0.358 
      //Rescale the param[0] with 1.03
      energy /= fNonLinearityParams[0]/(1+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]));
      
      break;
    }
      
    case kBeamTestCorrected:
    {
      // From beam test, corrected for material between beam and EMCAL
      //fNonLinearityParams[0] =  0.99078
      //fNonLinearityParams[1] =  0.161499;
      //fNonLinearityParams[2] =  0.655166; 
      //fNonLinearityParams[3] =  0.134101;
      //fNonLinearityParams[4] =  163.282;
      //fNonLinearityParams[5] =  23.6904;
      //fNonLinearityParams[6] =  0.978;
        energy *= fNonLinearityParams[6]/(fNonLinearityParams[0]*(1./(1.+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]))*1./(1.+fNonLinearityParams[3]*exp((energy-fNonLinearityParams[4])/fNonLinearityParams[5]))));

      break;
    }
     
    case kBeamTestCorrectedv2:
    {
      // From beam test, corrected for material between beam and EMCAL
      // Different parametrization to kBeamTestCorrected
      //fNonLinearityParams[0] =  0.983504;
      //fNonLinearityParams[1] =  0.210106;
      //fNonLinearityParams[2] =  0.897274;
      //fNonLinearityParams[3] =  0.0829064;
      //fNonLinearityParams[4] =  152.299;
      //fNonLinearityParams[5] =  31.5028;
      //fNonLinearityParams[6] =  0.968;
      energy *= fNonLinearityParams[6]/(fNonLinearityParams[0]*(1./(1.+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]))*1./(1.+fNonLinearityParams[3]*exp((energy-fNonLinearityParams[4])/fNonLinearityParams[5]))));
      
      break;
    }
      
    case kBeamTestCorrectedv3:
    {
      // Same function as kBeamTestCorrected, different default parametrization.
      //fNonLinearityParams[0] =  0.976941;
      //fNonLinearityParams[1] =  0.162310;
      //fNonLinearityParams[2] =  1.08689;
      //fNonLinearityParams[3] =  0.0819592;
      //fNonLinearityParams[4] =  152.338;
      //fNonLinearityParams[5] =  30.9594;
      //fNonLinearityParams[6] =  0.9615;
      energy *= fNonLinearityParams[6]/(fNonLinearityParams[0]*(1./(1.+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]))*1./(1.+fNonLinearityParams[3]*exp((energy-fNonLinearityParams[4])/fNonLinearityParams[5]))));
      
      break;
    }
      
    case kBeamTestCorrectedv4:
    {
      // New parametrization of kBeamTestCorrected, 
      // fitting new points for E>100 GeV.
      // I should have same performance as v3 in the low energies
      // See EMCal meeting 21/09/2018 slides
      // https://indico.cern.ch/event/759154/contributions/3148448/attachments/1721042/2778585/nonLinearityUpdate.pdf
      //  and jira ticket EMCAL-190
      // Not very smart copy pasting the same function for each new parametrization, need to think how to do it better.
      
//      fNonLinearityParams[0] = 0.9892;
//      fNonLinearityParams[1] = 0.1976;
//      fNonLinearityParams[2] = 0.865;
//      fNonLinearityParams[3] = 0.06775;
//      fNonLinearityParams[4] = 156.6;
//      fNonLinearityParams[5] = 47.18;
//      fNonLinearityParams[6] = 0.97;
      
      energy *= fNonLinearityParams[6]/(fNonLinearityParams[0]*(1./(1.+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]))*1./(1.+fNonLinearityParams[3]*exp((energy-fNonLinearityParams[4])/fNonLinearityParams[5]))));
      
      break;
    }
    case kBeamTestNS:
    {
      // New parametrization of testbeam data points, 
      // includes also points for E>100 GeV.
      // See EMCal meeting 07/12/2018 slides
      // https://indico.cern.ch/event/761682/contributions/3245317/attachments/1767706/2870846/2018_12_pp5TeV_NonlinearityStudies_update.pdf
      
//      fNonLinearityParams[0] = 0.986154;
//      fNonLinearityParams[1] = 0.214860;
//      fNonLinearityParams[2] = 0.717724;
//      fNonLinearityParams[3] = 0.069200;
//      fNonLinearityParams[4] = 155.497605;
//      fNonLinearityParams[5] = 48.868069;
//      fNonLinearityParams[6] = 0.972947;
      energy *= fNonLinearityParams[6]/(fNonLinearityParams[0]*(1./(1.+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]))*1./(1.+fNonLinearityParams[3]*exp((energy-fNonLinearityParams[4])/fNonLinearityParams[5]))));
      
      break;
    }
      
    case kSDMv5:
    {
      //Based on fit to the MC/data using kNoCorrection on the data - utilizes symmetric decay method and kPi0MCv5(MC) - 28 Oct 2013
      //fNonLinearityParams[0] =  1.0;
      //fNonLinearityParams[1] =  6.64778e-02;
      //fNonLinearityParams[2] =  1.570;
      //fNonLinearityParams[3] =  9.67998e-02;
      //fNonLinearityParams[4] =  2.19381e+02;
      //fNonLinearityParams[5] =  6.31604e+01;
      //fNonLinearityParams[6] =  1.01286;
      energy *= fNonLinearityParams[6]/(fNonLinearityParams[0]*(1./(1.+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]))*1./(1.+fNonLinearityParams[3]*exp((energy-fNonLinearityParams[4])/fNonLinearityParams[5])))) * (0.964 + exp(-3.132-0.435*energy*2.0));
      
      break;
    }
      
    case kPi0MCv5:
    {
      //Based on comparing MC truth information to the reconstructed energy of clusters. 
      //fNonLinearityParams[0] =  1.0;
      //fNonLinearityParams[1] =  6.64778e-02;
      //fNonLinearityParams[2] =  1.570;
      //fNonLinearityParams[3] =  9.67998e-02;
      //fNonLinearityParams[4] =  2.19381e+02;
      //fNonLinearityParams[5] =  6.31604e+01;
      //fNonLinearityParams[6] =  1.01286;
      energy *= fNonLinearityParams[6]/(fNonLinearityParams[0]*(1./(1.+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]))*1./(1.+fNonLinearityParams[3]*exp((energy-fNonLinearityParams[4])/fNonLinearityParams[5]))));
      
      break;
    }
      
    case kSDMv6:
    {
      //Based on fit to the MC/data using kNoCorrection on the data 
      //  - utilizes symmetric decay method and kPi0MCv6(MC) - 09 Dec 2014
      //  - parameters constrained by the test beam data as well
      // described in the note: https://aliceinfo.cern.ch/Notes/node/211 - Sec 3.1.2 (Test Beam Constrained SDM).
      //fNonLinearityParams[0] =  1.0;
      //fNonLinearityParams[1] =  0.237767;
      //fNonLinearityParams[2] =  0.651203;
      //fNonLinearityParams[3] =  0.183741;
      //fNonLinearityParams[4] =  155.427;
      //fNonLinearityParams[5] =  17.0335;
      //fNonLinearityParams[6] =  0.987054;
      energy *= fNonLinearityParams[6]/(fNonLinearityParams[0]*(1./(1.+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]))*1./(1.+fNonLinearityParams[3]*exp((energy-fNonLinearityParams[4])/fNonLinearityParams[5]))));
      
      break;
    }
      
    case kPi0MCv6:
    {
      //Based on comparing MC truth information to the reconstructed energy of clusters.
      // described in the note: https://aliceinfo.cern.ch/Notes/node/211 - Sec 3.1.2 (Test Beam Constrained SDM).
      //fNonLinearityParams[0] =  1.0;
      //fNonLinearityParams[1] =  0.0797873;
      //fNonLinearityParams[2] =  1.68322;
      //fNonLinearityParams[3] =  0.0806098;
      //fNonLinearityParams[4] =  244.586;
      //fNonLinearityParams[5] =  116.938;
      //fNonLinearityParams[6] =  1.00437;
      energy *= fNonLinearityParams[6]/(fNonLinearityParams[0]*(1./(1.+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]))*1./(1.+fNonLinearityParams[3]*exp((energy-fNonLinearityParams[4])/fNonLinearityParams[5]))));
      
      break;
    }
      
    case kPi0MCNS:
    {
      // New parametrization of testbeam MC points,
      // includes also points for E>100 GeV.
      // See EMCal meeting 07/12/2018 slides
      // https://indico.cern.ch/event/761682/contributions/3245317/attachments/1767706/2870846/2018_12_pp5TeV_NonlinearityStudies_update.pdf
      //fNonLinearityParams[0] =  1.009121;
      //fNonLinearityParams[1] =  0.083153;
      //fNonLinearityParams[2] =  1.444362;
      //fNonLinearityParams[3] =  0.100294;
      //fNonLinearityParams[4] =  416.897753;
      //fNonLinearityParams[5] =  324.246101;
      //fNonLinearityParams[6] =  1.004055;
      energy *= fNonLinearityParams[6]/(fNonLinearityParams[0]*(1./(1.+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]))*1./(1.+fNonLinearityParams[3]*exp((energy-fNonLinearityParams[4])/fNonLinearityParams[5]))));
      
      break;
    }

  case kPCMv1:
    {
      //based on symmetric decays of pi0 meson 
      // described in the note: https://aliceinfo.cern.ch/Notes/node/211 - Sec 3.1.2 (Test Beam Constrained SDM).
      // parameters vary from MC to MC
      //fNonLinearityParams[0] =   0.984876;
      //fNonLinearityParams[1] =  -9.999609;
      //fNonLinearityParams[2] =  -4.999891;
      //fNonLinearityParams[3] =  0.;
      //fNonLinearityParams[4] =  0.;
      //fNonLinearityParams[5] =  0.;
      //fNonLinearityParams[6] =  0.;
      energy /=  TMath::Power(fNonLinearityParams[0] + exp(fNonLinearityParams[1] + fNonLinearityParams[2]*energy),2);//result coming from calo-conv method
      
      break;
    }  

  case kPCMplusBTCv1:
    {
      //convolution of TestBeamCorrectedv3 with PCM method
      //Based on comparing MC truth information to the reconstructed energy of clusters.
      // described in the note: https://aliceinfo.cern.ch/Notes/node/211 - Sec 3.1.2 (Test Beam Constrained SDM).
      // parameters vary from MC to MC
      //fNonLinearityParams[0] =  0.976941;
      //fNonLinearityParams[1] =  0.162310;
      //fNonLinearityParams[2] =  1.08689;
      //fNonLinearityParams[3] =  0.0819592;
      //fNonLinearityParams[4] =  152.338;
      //fNonLinearityParams[5] =  30.9594;
      //fNonLinearityParams[6] =  0.9615;
      //fNonLinearityParams[7] =   0.984876;
      //fNonLinearityParams[8] =  -9.999609;
      //fNonLinearityParams[9] =  -4.999891;
      energy *= fNonLinearityParams[6]/(fNonLinearityParams[0]*(1./(1.+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]))*1./(1.+fNonLinearityParams[3]*exp((energy-fNonLinearityParams[4])/fNonLinearityParams[5]))));
      energy /= TMath::Power(fNonLinearityParams[7] + exp(fNonLinearityParams[8] + fNonLinearityParams[9]*energy),2);//result coming from calo-conv method
      
      break;
    }  

  case kPCMsysv1:
    {
      // Systematic variation of kPCMv1
      //Based on comparing MC truth information to the reconstructed energy of clusters.
      // described in the note: https://aliceinfo.cern.ch/Notes/node/211 - Sec 3.1.2 (Test Beam Constrained SDM).
      // parameters vary from MC to MC
      //fNonLinearityParams[0] =  0.0;
      //fNonLinearityParams[1] =  1.0;
      //fNonLinearityParams[2] =  1.0;
      //fNonLinearityParams[3] =  0.0;
      //fNonLinearityParams[4] =  1.0;
      //fNonLinearityParams[5] =  0.0;
      //fNonLinearityParams[6] =  0.0;
      energy /= TMath::Power( (fNonLinearityParams[0] + fNonLinearityParams[1] * TMath::Power(energy,fNonLinearityParams[2]) ) /
			      (fNonLinearityParams[3] + fNonLinearityParams[4] * TMath::Power(energy,fNonLinearityParams[5]) ) + fNonLinearityParams[6], 2);//result coming from calo-conv method
      
      break;
    }  
    case kTestBeamShaper:
    {
      // THIS PARAMETRIZATION HAS TO BE USED TOGETHER WITH THE SHAPER NONLINEARITY: 
      // Final parametrization of testbeam data points, 
      // includes also points for E>100 GeV and determined on shaper corrected data.

    //  fNonLinearityParams[0] = 1.91897;
    //  fNonLinearityParams[1] = 0.0264988;
    //  fNonLinearityParams[2] = 0.965663;
    //  fNonLinearityParams[3] = -187.501;
    //  fNonLinearityParams[4] = 2762.51;
      energy /= ( 1.0505 * (fNonLinearityParams[0] + fNonLinearityParams[1] * TMath::Log(energy) ) / ( 1 + ( fNonLinearityParams[2] * TMath::Exp( ( energy - fNonLinearityParams[3] ) / fNonLinearityParams[4] ) ) ) );
      
      break;
    }
    case kTestBeamFinalMC:
    {
      // Final parametrization of testbeam MC points, 
      // includes also points for E>100 GeV

    //  fNonLinearityParams[0] = 1.09357;
    //  fNonLinearityParams[1] = 0.0192266;
    //  fNonLinearityParams[2] = 0.291993;
    //  fNonLinearityParams[3] = 370.927;
    //  fNonLinearityParams[4] = 694.656;
      energy /= ( 1.00 * (fNonLinearityParams[0] + fNonLinearityParams[1] * TMath::Log(energy) ) / ( 1 + ( fNonLinearityParams[2] * TMath::Exp( ( energy - fNonLinearityParams[3] ) / fNonLinearityParams[4] ) ) ) );
      break;
    }
    
    case kNoCorrection:
      AliDebug(2,"No correction on the energy\n");
      break;
      
  }
  
  return energy;
}

///
/// Initialising non Linearity Parameters for the different
/// parametrizations available, defined in enum NonlinearityFunctions
///
//__________________________________________________
void AliEMCALRecoUtils::InitNonLinearityParam()
{
  if (fNonLinearityFunction == kPi0MC) {
    fNonLinearityParams[0] = 1.014;
    fNonLinearityParams[1] = -0.03329;
    fNonLinearityParams[2] = -0.3853;
    fNonLinearityParams[3] = 0.5423;
    fNonLinearityParams[4] = -0.4335;
  }
  
  if (fNonLinearityFunction == kPi0MCv2) {
    fNonLinearityParams[0] = 3.11111e-02;
    fNonLinearityParams[1] =-5.71666e-02; 
    fNonLinearityParams[2] = 5.67995e-01;      
  }
  
  if (fNonLinearityFunction == kPi0MCv3) {
    fNonLinearityParams[0] =  9.81039e-01;
    fNonLinearityParams[1] =  1.13508e-01;
    fNonLinearityParams[2] =  1.00173e+00; 
    fNonLinearityParams[3] =  9.67998e-02;
    fNonLinearityParams[4] =  2.19381e+02;
    fNonLinearityParams[5] =  6.31604e+01;
    fNonLinearityParams[6] =  1;
  }
  
  if (fNonLinearityFunction == kPi0GammaGamma) {
    fNonLinearityParams[0] = 1.04;
    fNonLinearityParams[1] = -0.1445;
    fNonLinearityParams[2] = 1.046;
  }  

  if (fNonLinearityFunction == kPi0GammaConversion) {
    fNonLinearityParams[0] = 0.139393;
    fNonLinearityParams[1] = 0.0566186;
    fNonLinearityParams[2] = 0.982133;
  }  

  if (fNonLinearityFunction == kBeamTest) {
    if (fNonLinearThreshold == 30) {
      fNonLinearityParams[0] = 1.007; 
      fNonLinearityParams[1] = 0.894; 
      fNonLinearityParams[2] = 0.246; 
    }
    if (fNonLinearThreshold == 45) {
      fNonLinearityParams[0] = 1.003; 
      fNonLinearityParams[1] = 0.719; 
      fNonLinearityParams[2] = 0.334; 
    }
    if (fNonLinearThreshold == 75) {
      fNonLinearityParams[0] = 1.002; 
      fNonLinearityParams[1] = 0.797; 
      fNonLinearityParams[2] = 0.358; 
    }
  }

  if (fNonLinearityFunction == kBeamTestCorrected) {
    fNonLinearityParams[0] =  0.99078;
    fNonLinearityParams[1] =  0.161499;
    fNonLinearityParams[2] =  0.655166; 
    fNonLinearityParams[3] =  0.134101;
    fNonLinearityParams[4] =  163.282;
    fNonLinearityParams[5] =  23.6904;
    fNonLinearityParams[6] =  0.978;
  }
  
  if (fNonLinearityFunction == kBeamTestCorrectedv2) {
    // Parameters until November 2015, use now kBeamTestCorrectedv3
    fNonLinearityParams[0] =  0.983504;
    fNonLinearityParams[1] =  0.210106;
    fNonLinearityParams[2] =  0.897274;
    fNonLinearityParams[3] =  0.0829064;
    fNonLinearityParams[4] =  152.299;
    fNonLinearityParams[5] =  31.5028;
    fNonLinearityParams[6] =  0.968;    
  }

  if (fNonLinearityFunction == kBeamTestCorrectedv3) {
    
    // New parametrization of kBeamTestCorrected
    // excluding point at 0.5 GeV from Beam Test Data
    // https://indico.cern.ch/event/438805/contribution/1/attachments/1145354/1641875/emcalPi027August2015.pdf
    
    fNonLinearityParams[0] =  0.976941;
    fNonLinearityParams[1] =  0.162310;
    fNonLinearityParams[2] =  1.08689;
    fNonLinearityParams[3] =  0.0819592;
    fNonLinearityParams[4] =  152.338;
    fNonLinearityParams[5] =  30.9594;
    fNonLinearityParams[6] =  0.9615;
  }

  if (fNonLinearityFunction == kBeamTestCorrectedv4) {
    
    // New parametrization of kBeamTestCorrected, 
    // fitting new points for E>100 GeV.
    // I should have same performance as v3 in the low energies
    // See EMCal meeting 21/09/2018 slides
    // https://indico.cern.ch/event/759154/contributions/3148448/attachments/1721042/2778585/nonLinearityUpdate.pdf
    //  and jira ticket EMCAL-190
    
    fNonLinearityParams[0] = 0.9892;
    fNonLinearityParams[1] = 0.1976;
    fNonLinearityParams[2] = 0.865;
    fNonLinearityParams[3] = 0.06775;
    fNonLinearityParams[4] = 156.6;
    fNonLinearityParams[5] = 47.18;
    fNonLinearityParams[6] = 0.97;
  }

  if (fNonLinearityFunction == kBeamTestNS) {
    
    // New parametrization of testbeam data points, 
    // includes also points for E>100 GeV.
    // See EMCal meeting 07/12/2018 slides
    // https://indico.cern.ch/event/761682/contributions/3245317/attachments/1767706/2870846/2018_12_pp5TeV_NonlinearityStudies_update.pdf
    
     fNonLinearityParams[0] = 0.986154;
     fNonLinearityParams[1] = 0.214860;
     fNonLinearityParams[2] = 0.717724;
     fNonLinearityParams[3] = 0.069200;
     fNonLinearityParams[4] = 155.497605;
     fNonLinearityParams[5] = 48.868069;
     fNonLinearityParams[6] = 0.972947;
  }
  
  if (fNonLinearityFunction == kSDMv5) {
    fNonLinearityParams[0] =  1.0;
    fNonLinearityParams[1] =  6.64778e-02;
    fNonLinearityParams[2] =  1.570;
    fNonLinearityParams[3] =  9.67998e-02;
    fNonLinearityParams[4] =  2.19381e+02;
    fNonLinearityParams[5] =  6.31604e+01;
    fNonLinearityParams[6] =  1.01286;
  }

  if (fNonLinearityFunction == kPi0MCv5) {
    fNonLinearityParams[0] =  1.0;
    fNonLinearityParams[1] =  6.64778e-02;
    fNonLinearityParams[2] =  1.570;
    fNonLinearityParams[3] =  9.67998e-02;
    fNonLinearityParams[4] =  2.19381e+02;
    fNonLinearityParams[5] =  6.31604e+01;
    fNonLinearityParams[6] =  1.01286;
  }

  if (fNonLinearityFunction == kSDMv6) {
    fNonLinearityParams[0] = 1.0;      
    fNonLinearityParams[1] = 0.237767; 
    fNonLinearityParams[2] = 0.651203; 
    fNonLinearityParams[3] = 0.183741; 
    fNonLinearityParams[4] = 155.427;  
    fNonLinearityParams[5] = 17.0335;  
    fNonLinearityParams[6] = 0.987054; 
  }

  if (fNonLinearityFunction == kPi0MCv6) {
    fNonLinearityParams[0] = 1.0;       
    fNonLinearityParams[1] = 0.0797873; 
    fNonLinearityParams[2] = 1.68322;   
    fNonLinearityParams[3] = 0.0806098; 
    fNonLinearityParams[4] = 244.586;   
    fNonLinearityParams[5] = 116.938;   
    fNonLinearityParams[6] = 1.00437;   
  }

  if (fNonLinearityFunction == kPi0MCNS) {
    fNonLinearityParams[0] =  1.009121;
    fNonLinearityParams[1] =  0.083153;
    fNonLinearityParams[2] =  1.444362;
    fNonLinearityParams[3] =  0.100294;
    fNonLinearityParams[4] =  416.897753;
    fNonLinearityParams[5] =  324.246101;
    fNonLinearityParams[6] =  1.004055;
  }

if (fNonLinearityFunction == kPCMv1) {
  //parameters change from MC production to MC production, they need to set for each period
    fNonLinearityParams[0] =  0.984876;
    fNonLinearityParams[1] = -9.999609;
    fNonLinearityParams[2] = -4.999891;
    fNonLinearityParams[3] = 0.;
    fNonLinearityParams[4] = 0.;
    fNonLinearityParams[5] = 0.;
    fNonLinearityParams[6] = 0.;
}

 if (fNonLinearityFunction == kPCMplusBTCv1) {
   // test beam corrected values convoluted with symmetric meson decays values
   // for test beam:
   // https://indico.cern.ch/event/438805/contribution/1/attachments/1145354/1641875/emcalPi027August2015.pdf
   // for PCM method:
   // https://aliceinfo.cern.ch/Notes/node/211
    fNonLinearityParams[0] =  0.976941;
    fNonLinearityParams[1] =  0.162310;
    fNonLinearityParams[2] =  1.08689;
    fNonLinearityParams[3] =  0.0819592;
    fNonLinearityParams[4] =  152.338;
    fNonLinearityParams[5] =  30.9594;
    fNonLinearityParams[6] =  0.9615;
    fNonLinearityParams[7] =   0.984876;
    fNonLinearityParams[8] =  -9.999609;
    fNonLinearityParams[9] =  -4.999891;
 }

 if (fNonLinearityFunction == kPCMsysv1) {
   //systematics for kPCMv1
   // for PCM method:
   // https://aliceinfo.cern.ch/Notes/node/211
   fNonLinearityParams[0] =  0.0;
   fNonLinearityParams[1] =  1.0;
   fNonLinearityParams[2] =  1.0;
   fNonLinearityParams[3] =  0.0;
   fNonLinearityParams[4] =  1.0;
   fNonLinearityParams[5] =  0.0;
   fNonLinearityParams[6] =  0.0;
 }

 if (fNonLinearityFunction == kTestBeamShaper) {
   fNonLinearityParams[0] =  1.91897;
   fNonLinearityParams[1] =  0.0264988;
   fNonLinearityParams[2] =  0.965663;
   fNonLinearityParams[3] =  -187.501;
   fNonLinearityParams[4] =  2762.51;
 }
 if (fNonLinearityFunction == kTestBeamFinalMC) {
   fNonLinearityParams[0] =  1.09357;
   fNonLinearityParams[1] =  0.0192266;
   fNonLinearityParams[2] =  0.291993;
   fNonLinearityParams[3] =  370.927;
   fNonLinearityParams[4] =  694.656;
 }
}

///
/// Set the type of channels to be declared as bad
/// \param all  : all cases are bad, default true
/// \param dead : dead channels are bad
/// \param hot  : hot channels are bad
/// \param warm : warm channels are bad
///
//____________________________________________________________________
void AliEMCALRecoUtils::SetEMCALBadChannelStatusSelection(Bool_t all, Bool_t dead, Bool_t hot, Bool_t warm)
{ 
  fBadStatusSelection[0] = all;  // declare all as bad if true, never mind the other settings
  fBadStatusSelection[1] = dead; 
  fBadStatusSelection[2] = hot; 
  fBadStatusSelection[3] = warm; 
}

///
/// \return declare channel as bad (true) or not good (false)
/// By default if status is not kAlive, all are declared bad,
/// but optionnaly 
///
/// \param iSM: supermodule number of channel
/// \param iCol: cell column in SM
/// \param iRow: cell row in SM
/// \param status: channel status
///
//____________________________________________________________________
Bool_t AliEMCALRecoUtils::GetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow, Int_t & status) const 
{ 
  if(!fUse1Dmap){
    if(fEMCALBadChannelMap) 
      status = (Int_t) ((TH2I*)fEMCALBadChannelMap->At(iSM))->GetBinContent(iCol,iRow); 
    else 
      status = 0; // Channel is ok by default
  }else{
    AliEMCALGeometry* geom = AliEMCALGeometry::GetInstance();
    Int_t CellID = geom->GetAbsCellIdFromCellIndexes(iSM, iRow, iCol);
    if(fEMCALBadChannelMap) 
      status = (Int_t) ((TH1C*)fEMCALBadChannelMap->At(0))->GetBinContent(CellID); 
    else 
      status = 0; // Channel is ok by default
  }

  if ( status == AliCaloCalibPedestal::kAlive ) 
  {
    return kFALSE; // Good channel
  }
  else
  {
    if      ( fBadStatusSelection[0]  == kTRUE ) 
    {
      return kTRUE; // consider bad hot, dead and warm
    }
    else
    {
      if      ( fBadStatusSelection[AliCaloCalibPedestal::kDead]    == kTRUE  && 
                status == AliCaloCalibPedestal::kDead    ) 
        return kTRUE; // consider bad dead
      else if ( fBadStatusSelection[AliCaloCalibPedestal::kHot]     == kTRUE  && 
                status == AliCaloCalibPedestal::kHot     ) 
        return kTRUE; // consider bad hot
      else if ( fBadStatusSelection[AliCaloCalibPedestal::kWarning] == kTRUE  && 
                status == AliCaloCalibPedestal::kWarning ) 
        return kTRUE; // consider bad warm 
    }
  }
  
  AliWarning(Form("Careful, bad channel selection not properly done: ism %d, icol %d, irow %d, status %d,\n"
                  " fBadAll %d, fBadHot %d, fBadWarm %d, fBadDead %d",
                  iSM, iCol, iRow, status,
                  fBadStatusSelection[0], fBadStatusSelection[1],
                  fBadStatusSelection[2], fBadStatusSelection[3]));
  
  return kFALSE; // if everything fails, accept it.
}

///
/// \return declare channel as bad (true) or not good (false)
/// By default if status is not kAlive, all are declared bad,
/// but optionnaly 
///
/// \param iCell: cell ID
/// \param status: channel status
///
//____________________________________________________________________
Bool_t AliEMCALRecoUtils::GetEMCALChannelStatus1D(Int_t iCell, Int_t & status) const 
{ 
  if(fEMCALBadChannelMap) 
    status = (Int_t) ((TH1C*)fEMCALBadChannelMap->At(0))->GetBinContent(iCell); 
  else 
    status = 0; // Channel is ok by default
  
  if ( status == AliCaloCalibPedestal::kAlive ) 
  {
    return kFALSE; // Good channel
  }
  else
  {
    if      ( fBadStatusSelection[0]  == kTRUE ) 
    {
      return kTRUE; // consider bad hot, dead and warm
    }
    else
    {
      if      ( fBadStatusSelection[AliCaloCalibPedestal::kDead]    == kTRUE  && 
                status == AliCaloCalibPedestal::kDead    ) 
        return kTRUE; // consider bad dead
      else if ( fBadStatusSelection[AliCaloCalibPedestal::kHot]     == kTRUE  && 
                status == AliCaloCalibPedestal::kHot     ) 
        return kTRUE; // consider bad hot
      else if ( fBadStatusSelection[AliCaloCalibPedestal::kWarning] == kTRUE  && 
                status == AliCaloCalibPedestal::kWarning ) 
        return kTRUE; // consider bad warm 
    }
  }
  
  AliWarning(Form("Careful, bad channel selection not properly done: icell %d, status %d,\n"
                  " fBadAll %d, fBadHot %d, fBadWarm %d, fBadDead %d",
                  iCell, status,
                  fBadStatusSelection[0], fBadStatusSelection[1],
                  fBadStatusSelection[2], fBadStatusSelection[3]));
  
  return kFALSE; // if everything fails, accept it.
}

///
/// For a given CaloCluster gets the absId of the cell 
/// with maximum energy deposit.
///
/// \param geom: AliEMCALGeometry pointer
/// \param cells: full list of cells
/// \param clu: pointer to AliVCluster
/// \param absId: absolute id number of cell with highest energy in cluster
/// \param iSupMod: supermodule number of cell with highest energy in cluster
/// \param ieta: column number of cell with highest energy in cluster
/// \param iphi: row number of cell with highest energy in cluster
/// \param shared: cluster is shared between 2 supermodules
///
//____________________________________________________________________
void AliEMCALRecoUtils::GetMaxEnergyCell(const AliEMCALGeometry *geom, 
                                         AliVCaloCells* cells, 
                                         const AliVCluster* clu, 
                                         Int_t  & absId,  
                                         Int_t  & iSupMod, 
                                         Int_t  & ieta, 
                                         Int_t  & iphi, 
                                         Bool_t & shared)
{  
  Double_t eMax        = -1.;
  Double_t eCell       = -1.;
  Float_t  fraction    = 1.;
  Float_t  recalFactor = 1.;
  Int_t    cellAbsId   = -1 ;

  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iSupMod0= -1;

  if (!clu) 
  {
    AliInfo("Cluster pointer null!");
    absId=-1; iSupMod0=-1, ieta = -1; iphi = -1; shared = -1;
    return;
  }
  
  for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++) 
  {
    cellAbsId = clu->GetCellAbsId(iDig);
    fraction  = clu->GetCellAmplitudeFraction(iDig);
    //printf("a Cell %d, id, %d, amp %f, fraction %f\n",iDig,cellAbsId,cells->GetCellAmplitude(cellAbsId),fraction);
    if (fraction < 1e-4) fraction = 1.; // in case unfolding is off
    
    geom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta); 
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
    
    if (iDig==0) 
    {
      iSupMod0=iSupMod;
    } else if (iSupMod0!=iSupMod)  
    {
      shared = kTRUE;
      //printf("AliEMCALRecoUtils::GetMaxEnergyCell() - SHARED CLUSTER\n");
    }

    if (!fCellsRecalibrated && IsRecalibrationOn()){ 
      if(fUse1Drecalib)
        recalFactor = GetEMCALChannelRecalibrationFactor1D(cellAbsId);
      else
        recalFactor = GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
    }
    
    eCell  = cells->GetCellAmplitude(cellAbsId)*fraction*recalFactor;
    //printf("b Cell %d, id, %d, amp %f, fraction %f\n",iDig,cellAbsId,eCell,fraction);
    if (eCell > eMax) 
    { 
      eMax  = eCell; 
      absId = cellAbsId;
      //printf("\t new max: cell %d, e %f, ecell %f\n",maxId, eMax,eCell);
    }
  }// cell loop
  
  //Get from the absid the supermodule, tower and eta/phi numbers
  geom->GetCellIndex(absId,iSupMod,iTower,iIphi,iIeta); 
  //Gives SuperModule and Tower numbers
  geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta); 
  //printf("Max id %d, iSM %d, col %d, row %d\n",absId,iSupMod,ieta,iphi);
  //printf("Max end---\n");
}

///
/// \return weight of cell for shower shape calculation
/// If fW0 parameter is negative, apply log weight without trimming.
///
/// \param eCell: cluster cell energy
/// \param eCluster: cluster Energy
/// \param iSM: supermodule number
///
//_________________________________________________________
Float_t  AliEMCALRecoUtils::GetCellWeight(Float_t eCell, Float_t eCluster) const 
{ 
  if (eCell > 0 && eCluster > 0) 
  {
   if ( fW0 > 0 ) return TMath::Max( 0., fW0 + TMath::Log( eCell / eCluster ) ) ;
   else           return TMath::Log( eCluster / eCell ) ; 
  }
  else                           
    return 0. ; 
}

///
/// Initialize data members with default values
///
//______________________________________
void AliEMCALRecoUtils::InitParameters()
{  
  fParticleType = kPhoton;
  fPosAlgo      = kUnchanged;
  fW0           = 4.5;
  
  fNonLinearityFunction = kNoCorrection;
  fNonLinearThreshold   = 30;
  
  fExoticCellFraction     = 0.97;
  fExoticCellDiffTime     = 1e6;
  fExoticCellMinAmplitude = 4.0;
  
  fAODFilterMask    = 128;
  fAODHybridTracks  = kFALSE;
  fAODTPCOnlyTracks = kTRUE;
  
  fCutEtaPhiSum      = kTRUE;
  fCutEtaPhiSeparate = kFALSE;
  
  fCutR   = 0.05; 
  fCutEta = 0.025; 
  fCutPhi = 0.05;
  
  fClusterWindow = 100;
  fMass          = 0.139;
  
  fStepSurface   = 20.;                      
  fStepCluster   = 5.;
  fTrackCutsType = kLooseCut;
  
  fCutMinTrackPt     = 0;
  fCutMinNClusterTPC = -1;
  fCutMinNClusterITS = -1;
  
  fCutMaxChi2PerClusterTPC  = 1e10;
  fCutMaxChi2PerClusterITS  = 1e10;
  
  fCutRequireTPCRefit     = kFALSE;
  fCutRequireITSRefit     = kFALSE;
  fCutAcceptKinkDaughters = kFALSE;
  
  fCutMaxDCAToVertexXY = 1e10;             
  fCutMaxDCAToVertexZ  = 1e10;              
  fCutDCAToVertex2D    = kFALSE;
  
  fCutRequireITSStandAlone = kFALSE; //MARCEL
  fCutRequireITSpureSA     = kFALSE; //Marcel
  
  // Misalignment matrices
  for (Int_t i = 0; i < 15 ; i++) 
  {
    fMisalTransShift[i] = 0.; 
    fMisalRotShift[i]   = 0.; 
  }
  
  // Non linearity
  for (Int_t i = 0; i < 10  ; i++) fNonLinearityParams[i] = 0.; 
  
  // For kBeamTestCorrectedv2 case, but default is no correction
  fNonLinearityParams[0] =  0.983504;
  fNonLinearityParams[1] =  0.210106;
  fNonLinearityParams[2] =  0.897274;
  fNonLinearityParams[3] =  0.0829064;
  fNonLinearityParams[4] =  152.299;
  fNonLinearityParams[5] =  31.5028;
  fNonLinearityParams[6] =  0.968;
  
  // Cluster energy smearing
  fSmearClusterEnergy   = kFALSE;
  fSmearClusterParam[0] = 0.07; // * sqrt E term
  fSmearClusterParam[1] = 0.00; // * E term
  fSmearClusterParam[2] = 0.00; // constant
}

///
/// Init EMCAL energy calibration factors container
///
//_____________________________________________________
void AliEMCALRecoUtils::InitEMCALRecalibrationFactors()
{
  AliDebug(2,"AliCalorimeterUtils::InitEMCALRecalibrationFactors()");
  
  // In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fEMCALRecalibrationFactors = new TObjArray(22);
  for (int i = 0; i < 22; i++) 
    fEMCALRecalibrationFactors->Add(new TH2F(Form("EMCALRecalFactors_SM%d",i),
                                             Form("EMCALRecalFactors_SM%d",i),  48, 0, 48, 24, 0, 24));
  //Init the histograms with 1
  for (Int_t sm = 0; sm < 22; sm++) 
  {
    for (Int_t i = 0; i < 48; i++) 
    {
      for (Int_t j = 0; j < 24; j++) 
      {
        SetEMCALChannelRecalibrationFactor(sm,i,j,1.);
      }
    }
  }
  
  fEMCALRecalibrationFactors->SetOwner(kTRUE);
  fEMCALRecalibrationFactors->Compress();
  
  // In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);    
}
///
/// Init 1D EMCAL energy calibration factors container
///
//_____________________________________________________
void AliEMCALRecoUtils::InitEMCALRecalibrationFactors1D()
{
  AliDebug(2,"AliCalorimeterUtils::InitEMCALRecalibrationFactors1D()");
  
  // In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fEMCALRecalibrationFactors = new TObjArray(1);

  fEMCALRecalibrationFactors->Add(new TH1S("EMCALRecalFactors","EMCALRecalFactors", 48*24*22,0.,48*24*22));

  //Init the histograms with 1
  for (UInt_t icell = 0; icell < 48*24*22; icell++) 
  {
    SetEMCALChannelRecalibrationFactor1D(icell,1.);
  }
  
  fEMCALRecalibrationFactors->SetOwner(kTRUE);
  fEMCALRecalibrationFactors->Compress();
  
  // In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);    
}

///
/// Init EMCAL single channel energy calibration factors container
///
//_____________________________________________________
void AliEMCALRecoUtils::InitEMCALSingleChannelRecalibrationFactors()
{
  AliDebug(2,"AliCalorimeterUtils::InitEMCALSingleChannelRecalibrationFactors()");
  
  // In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fEMCALSingleChannelRecalibrationFactors = new TObjArray(22);
  for (int i = 0; i < 22; i++) 
    fEMCALSingleChannelRecalibrationFactors->Add(new TH2F(Form("EMCALSCCalibMap_Mod%d",i),
							  Form("EMCALSCCalibMap_Mod%d",i),  48, 0, 48, 24, 0, 24));
  //Init the histograms with 1
  for (Int_t sm = 0; sm < 22; sm++) 
  {
    for (Int_t i = 0; i < 48; i++) 
    {
      for (Int_t j = 0; j < 24; j++) 
      {
        SetEMCALSingleChannelRecalibrationFactor(sm,i,j,1.);
      }
    }
  }
  
  fEMCALSingleChannelRecalibrationFactors->SetOwner(kTRUE);
  fEMCALSingleChannelRecalibrationFactors->Compress();
  
  // In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);    
}

///
/// Init EMCAL time calibration shifts container
///
//_________________________________________________________
void AliEMCALRecoUtils::InitEMCALTimeRecalibrationFactors()
{
  AliDebug(2,"AliCalorimeterUtils::InitEMCALRecalibrationFactors()");

  // In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  if(fDoUseMergedBC){

    if(fLowGain) fEMCALTimeRecalibrationFactors = new TObjArray(2);
    else fEMCALTimeRecalibrationFactors = new TObjArray(1);

    fEMCALTimeRecalibrationFactors->Add(new TH1S("hAllTimeAv",
                                                 "hAllTimeAv",  
                                                 48*24*22,0.,48*24*22)          );
    // Init the histograms with 0
    for (Int_t iCh = 0; iCh < 48*24*22; iCh++) 
      SetEMCALChannelTimeRecalibrationFactor(0,iCh,0.,kFALSE);

    if(fLowGain) {
      fEMCALTimeRecalibrationFactors->Add(new TH1F("hAllTimeAvLG",
                                                   "hAllTimeAvLG",  
                                                    48*24*22,0.,48*24*22)        );
      for (Int_t iCh = 0; iCh < 48*24*22; iCh++) 
        SetEMCALChannelTimeRecalibrationFactor(1,iCh,0.,kTRUE);
    }

  }else{
    if(fLowGain) fEMCALTimeRecalibrationFactors = new TObjArray(8);
    else fEMCALTimeRecalibrationFactors = new TObjArray(4);

    for (int i = 0; i < 4; i++) 
      fEMCALTimeRecalibrationFactors->Add(new TH1F(Form("hAllTimeAvBC%d",i),
                                                 Form("hAllTimeAvBC%d",i),  
                                                 48*24*22,0.,48*24*22)          );
    // Init the histograms with 0
    for (Int_t iBC = 0; iBC < 4; iBC++) 
    {
      for (Int_t iCh = 0; iCh < 48*24*22; iCh++) 
        SetEMCALChannelTimeRecalibrationFactor(iBC,iCh,0.,kFALSE);
    }

    if(fLowGain) {
      for (int iBC = 0; iBC < 4; iBC++) {
        fEMCALTimeRecalibrationFactors->Add(new TH1F(Form("hAllTimeAvLGBC%d",iBC),
                                                     Form("hAllTimeAvLGBC%d",iBC),  
                                                     48*24*22,0.,48*24*22)        );
        for (Int_t iCh = 0; iCh < 48*24*22; iCh++) 
          SetEMCALChannelTimeRecalibrationFactor(iBC,iCh,0.,kTRUE);
      }
    }

  }
  
  fEMCALTimeRecalibrationFactors->SetOwner(kTRUE);
  fEMCALTimeRecalibrationFactors->Compress();
  
  // In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);    
}

///
/// Init EMCAL bad channels map container
///
//____________________________________________________
void AliEMCALRecoUtils::InitEMCALBadChannelStatusMap()
{
  AliDebug(2,"AliEMCALRecoUtils::InitEMCALBadChannelStatusMap()");

  // In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fEMCALBadChannelMap = new TObjArray(22);
  //TH2F * hTemp = new  TH2I("EMCALBadChannelMap","EMCAL SuperModule bad channel map", 48, 0, 48, 24, 0, 24);
  
  for (int i = 0; i < 22; i++) 
    fEMCALBadChannelMap->Add(new TH2I(Form("EMCALBadChannelMap_Mod%d",i),Form("EMCALBadChannelMap_Mod%d",i), 48, 0, 48, 24, 0, 24));
  
  fEMCALBadChannelMap->SetOwner(kTRUE);
  fEMCALBadChannelMap->Compress();
  
  // In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);    
}

///
/// Init EMCAL bad channels 1 dimensional map container
///
//____________________________________________________
void AliEMCALRecoUtils::InitEMCALBadChannelStatusMap1D()
{
  AliDebug(2,"AliEMCALRecoUtils::InitEMCALBadChannelStatusMap1D()");

  fUse1Dmap = kTRUE;
  // In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fEMCALBadChannelMap = new TObjArray(1);
  
  fEMCALBadChannelMap->Add(new TH1C("EMCALBadChannelMap","EMCALBadChannelMap", 48*24*22,0.,48*24*22));
  
  fEMCALBadChannelMap->SetOwner(kTRUE);
  fEMCALBadChannelMap->Compress();
  
  // In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);    
}

///
/// Init EMCAL L1 phase shifts container
///
//___________________________________________________________
void AliEMCALRecoUtils::InitEMCALL1PhaseInTimeRecalibration()
{
  AliDebug(2,"AliEMCALRecoUtils::InitEMCALL1PhaseInTimeRecalibrationFactors()");
 
  // In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fEMCALL1PhaseInTimeRecalibration = new TObjArray(1);

  fEMCALL1PhaseInTimeRecalibration->Add(new TH1C("h0","EMCALL1phaseForSM", 22, 0, 22));
  
  for (Int_t i = 0; i < 22; i++) //loop over SMs, default value = 0
    SetEMCALL1PhaseInTimeRecalibrationForSM(i,0,0);
  
  fEMCALL1PhaseInTimeRecalibration->SetOwner(kTRUE);
  fEMCALL1PhaseInTimeRecalibration->Compress();
  
  // In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);    
}

///
/// Recalibrate the cluster energy and time, considering the recalibration map 
/// and the time and energy of the cells that compose the cluster.
///
/// \param geom: pointer to geometry
/// \param cluster: pointer to cluster
/// \param cells: list of cells
/// \param bc: bunch crossing number returned by esdevent->GetBunchCrossNumber()
///
//____________________________________________________________________________
void AliEMCALRecoUtils::RecalibrateClusterEnergy(const AliEMCALGeometry* geom, 
                                                 AliVCluster * cluster, 
                                                 AliVCaloCells * cells, 
                                                 Int_t bc)
{  
  if (!cluster) 
  {
    AliInfo("Cluster pointer null!");
    return;
  }  
  
  // Get the cluster number of cells and list of absId, check what kind of cluster do we have.
  UShort_t * index    = cluster->GetCellsAbsId() ;
  Double_t * fraction = cluster->GetCellsAmplitudeFraction() ;
  Int_t ncells = cluster->GetNCells();
  
  // Initialize some used variables
  Float_t energy = 0;
  Int_t   absId  =-1;
  Int_t   icol   =-1, irow =-1, imod=1;
  Float_t factor = 1, frac = 0;
  Float_t factorSC = 1;
  Int_t   absIdMax = -1;
  Float_t emax     = 0;
  
  // Loop on the cells, get the cell amplitude and recalibration factor, multiply and and to the new energy
  for (Int_t icell = 0; icell < ncells; icell++)
  {
    absId = index[icell];
    frac =  fraction[icell];
    if (frac < 1e-5) frac = 1; //in case of EMCAL, this is set as 0 since unfolding is off
    
    if (!fCellsRecalibrated && IsRecalibrationOn()) 
    {
      // Energy  
      Int_t iTower = -1, iIphi = -1, iIeta = -1; 
      geom->GetCellIndex(absId,imod,iTower,iIphi,iIeta); 
      if (fEMCALRecalibrationFactors->GetEntries() <= imod) 
        continue;
      geom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,irow,icol);

      if(fUse1Drecalib)
        factor = GetEMCALChannelRecalibrationFactor1D(absId);
      else
        factor = GetEMCALChannelRecalibrationFactor(imod,icol,irow);
      
      AliDebug(2,Form("AliEMCALRecoUtils::RecalibrateClusterEnergy - recalibrate cell: module %d, col %d, row %d, cell fraction %f,recalibration factor %f, cell energy %f\n",
                      imod,icol,irow,frac,factor,cells->GetCellAmplitude(absId)));
      
    }
    if (!fCellsSingleChannelRecalibrated && IsSingleChannelRecalibrationOn()) 
    {
      // Single Channel  
      Int_t iTower = -1, iIphi = -1, iIeta = -1; 
      geom->GetCellIndex(absId,imod,iTower,iIphi,iIeta); 
      if (fEMCALSingleChannelRecalibrationFactors->GetEntries() <= imod) 
        continue;
      geom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,irow,icol);      
      factorSC = GetEMCALSingleChannelRecalibrationFactor(imod,icol,irow);
      
      AliDebug(2,Form("AliEMCALRecoUtils::RecalibrateClusterEnergy - sinlge channel recalibrate cell: module %d, col %d, row %d, cell fraction %f,recalibration factor %f, cell energy %f\n",
                      imod,icol,irow,frac,factorSC,cells->GetCellAmplitude(absId)));
      
    }
    
    energy += cells->GetCellAmplitude(absId)*factor*factorSC*frac;
    
    if (emax < cells->GetCellAmplitude(absId)*factor*factorSC*frac) 
    {
      emax     = cells->GetCellAmplitude(absId)*factor*factorSC*frac;
      absIdMax = absId;
    }
  }
  
  AliDebug(2,Form("AliEMCALRecoUtils::RecalibrateClusterEnergy - Energy before %f, after %f \n",cluster->E(),energy));
  
  cluster->SetE(energy);
  
  // Recalculate time of cluster
  Double_t timeorg = cluster->GetTOF();
  Bool_t isLowGain = !(cells->GetCellHighGain(absIdMax));//HG = false -> LG = true

  Double_t time = cells->GetCellTime(absIdMax);
  if (!fCellsRecalibrated && IsTimeRecalibrationOn())
    RecalibrateCellTime(absIdMax,bc,time,isLowGain);
  time-=fConstantTimeShift*1e-9; // only in case of Run1 old simulations
  
  // Recalibrate time with L1 phase 
  if (!fCellsRecalibrated && IsL1PhaseInTimeRecalibrationOn())
    RecalibrateCellTimeL1Phase(imod, bc, time, fCurrentParNumber);
  
  cluster->SetTOF(time);
  
  AliDebug(2,Form("AliEMCALRecoUtils::RecalibrateClusterEnergy - Time before %f, after %f \n",timeorg,cluster->GetTOF()));
}

///
/// Recalibrate all the cells time and energy, considering the recalibration map and 
/// the energy and time of each cell.
///
/// \param cells: list of cells
/// \param bc: bunch crossing number returned by esdevent->GetBunchCrossNumber()
///
//_______________________________________________________________________
void AliEMCALRecoUtils::RecalibrateCells(AliVCaloCells * cells, Int_t bc)
{
  if (!IsRecalibrationOn() && !IsTimeRecalibrationOn() && !IsBadChannelsRemovalSwitchedOn()) 
    return;
  
  if (!cells) 
  {
    AliInfo("Cells pointer null!");
    return;
  }  
  
  Short_t  absId  =-1;
  Bool_t   accept = kFALSE;
  Float_t  ecell  = 0;
  Double_t tcell  = 0;
  Double_t ecellin = 0;
  Double_t tcellin = 0;
  Int_t  mclabel = -1;
  Double_t efrac = 0;
  
  Int_t nEMcell  = cells->GetNumberOfCells() ;  
  for (Int_t iCell = 0; iCell < nEMcell; iCell++) 
  { 
    cells->GetCell( iCell, absId, ecellin, tcellin, mclabel, efrac );
    
    accept = AcceptCalibrateCell(absId, bc, ecell ,tcell ,cells); 
    if (!accept)
    {
      ecell = 0;
      tcell = -1;
    }
    
    // Set new values
    cells->SetCell(iCell,absId,ecell, tcell, mclabel, efrac);
  }

  fCellsRecalibrated = kTRUE;
}

///
/// Recalibrate all the cells with energy>40 GeV for the shaper nonlinearity
///
/// \param Emeas: energy of cell
/// \param EcalibHG: HG energy calibration coefficient
///
//_______________________________________________________________________
Float_t AliEMCALRecoUtils::CorrectShaperNonLin(Float_t Emeas, Float_t EcalibHG)
{
  // The following conversion factor needs to be applied to go from energy to ADC
  // AliEMCALCalibData::fADCchannelRef = 0.0162;

  if(Emeas<40){
    return Emeas*16.3/16;
  }
  Float_t par[]={1, 29.8279, 0.607704, 0.00164896, -2.28595e-06, -8.54664e-10, 5.50191e-12, -3.28098e-15};
  Float_t x = par[0]*Emeas/EcalibHG/16/0.0162;

  Float_t res=par[1];
  res+=par[2]*x;
  res+=par[3]*x*x;
  res+=par[4]*x*x*x;
  res+=par[5]*x*x*x*x;
  res+=par[6]*x*x*x*x*x;
  res+=par[7]*x*x*x*x*x*x;

  return  EcalibHG*16.3*res*0.0162;
}

///
/// Recalibrate time of cell from AbsID number considering cell calibration map 
///
/// \param absId: cell absolute ID number
/// \param bc: bunch crossing number returned by esdevent->GetBunchCrossNumber()
/// \param celltime: cell time to be returned calibrated
/// \param isLGon: low gain time calibration on/off
///
//_______________________________________________________________________________________________________
void AliEMCALRecoUtils::RecalibrateCellTime(Int_t absId, Int_t bc, Double_t & celltime, Bool_t isLGon) const
{  
  if (!fCellsRecalibrated && IsTimeRecalibrationOn() && bc >= 0) {
    if(fLowGain)
      celltime -= GetEMCALChannelTimeRecalibrationFactor(bc%4,absId,isLGon)*1.e-9;
    else
      celltime -= GetEMCALChannelTimeRecalibrationFactor(bc%4,absId,kFALSE)*1.e-9;
  }
}

///
/// Recalibrate time of cell from SM number considering the L1 phase shift 
///
/// \param iSM: supermodule number
/// \param bc: bunch crossing number returned by esdevent->GetBunchCrossNumber()
/// \param celltime: cell time to be returned calibrated
/// \param par: Int_t, in case of PAR load another set of L1 shifts, 0-no or before PAR, 1-after 1st PAR etc
///
//_______________________________________________________________________________________________________
void AliEMCALRecoUtils::RecalibrateCellTimeL1Phase(Int_t iSM, Int_t bc, Double_t & celltime, Short_t par) const
{
  if (!fCellsRecalibrated && IsL1PhaseInTimeRecalibrationOn() && bc >= 0) 
  {
    bc=bc%4;

    Float_t offsetPerSM=0.;
    Int_t l1PhaseShift = GetEMCALL1PhaseInTimeRecalibrationForSM(iSM,par);
    Int_t l1Phase=l1PhaseShift & 3; //bit operation

    if(bc >= l1Phase)
      offsetPerSM = (bc - l1Phase)*25;
    else
      offsetPerSM = (bc - l1Phase + 4)*25;

    Int_t l1shiftOffset=l1PhaseShift>>2; //bit operation
    l1shiftOffset*=25;
   
    celltime -= offsetPerSM*1.e-9;
    celltime -= l1shiftOffset*1.e-9;
  }
}

/// 2 tasks:
///    * Recover cell MC labels from the original cluster from the fraction of 
///      deposited energy to pass them to the digitizer
///    * If added generators have to be removed, check the origin of the label 
///      and depending on the deposited energy, correct the amplitude. Quite crude.
///
/// This only works for MC productions done with aliroot release larger than v5-08
///
/// \param clus     Input AliVCaloCluster with the list of cell MC labels
/// \param mc       MC event pointer, to identify the generator
/// \param absID    ID of the cell
/// \param amp      amplitude of the cell, to be recalculated if extra generators are removed
/// \param labeArr  list of MC labels associated to the cell
/// \param eDepArr  list of MC energy depositions in the cell corresponding to each MC label
///
//______________________________________________________________________________
void AliEMCALRecoUtils::RecalculateCellLabelsRemoveAddedGenerator( Int_t absID, AliVCluster* clus, AliMCEvent* mc,
                                                                   Float_t & amp, TArrayI & labeArr, TArrayF & eDepArr) const
{
  TString genName ;
  Float_t eDepFrac[4];
  
  Float_t edepTotFrac = 1;
  Bool_t  found       = kFALSE;
  Float_t ampOrg      = amp;
  
  //
  // Get the energy deposition fraction from cluster.
  //
  for(Int_t icluscell = 0; icluscell < clus->GetNCells(); icluscell++ )
  {
    if(absID == clus->GetCellAbsId(icluscell)) 
    {
      clus->GetCellMCEdepFractionArray(icluscell,eDepFrac);
      
      found = kTRUE;
      
      break;
    }
  } 
  
  if ( !found )
  {
    AliWarning(Form("Cell abs ID %d NOT found in cluster",absID));
    return;
  } 
  
  //
  // Check if there is a particle contribution from a given generator name.
  // If it is not one of the selected generators, 
  // remove the constribution from the cell.
  //
  if ( mc && fNMCGenerToAccept > 0 )
  {
    //printf("Accept contribution from generator? \n");
    for(UInt_t imc = 0; imc < 4; imc++)
    {
      if ( eDepFrac[imc] > 0 && clus->GetNLabels() > imc )
      {
        mc->GetCocktailGenerator(clus->GetLabelAt(imc),genName);
                
        Bool_t generOK = kFALSE;
        for(Int_t ig = 0; ig < fNMCGenerToAccept; ig++) 
        {
          if ( genName.Contains(fMCGenerToAccept[ig]) ) generOK = kTRUE;
        }
        
        if ( !generOK )
        {
          amp-=ampOrg*eDepFrac[imc];
          
          edepTotFrac-=eDepFrac[imc];
          
          eDepFrac[imc] = 0;
        }
        
      } // eDep > 0      
    } // 4 possible loop
    
  } // accept at least one kind of generator
  
  //
  // Add MC label and Edep to corresponding array (to be used later in digits)
  //
  Int_t nLabels = 0;
  for(UInt_t imc = 0; imc < 4; imc++)
  {
    if ( eDepFrac[imc] > 0 && clus->GetNLabels() > imc && edepTotFrac > 0 )
    {
      nLabels++;
      
      labeArr.Set(nLabels);
      labeArr.AddAt(clus->GetLabelAt(imc), nLabels-1);
      
      eDepArr.Set(nLabels);
      eDepArr.AddAt( (eDepFrac[imc]/edepTotFrac) * amp, nLabels-1);
      // use as deposited energy a fraction of the simulated energy (smeared and with noise)
    }  // edep > 0
    
  } // mc cell label loop
  
//  Commented on 12/05/20, in embedding case, data input cells are removed
//  Also, it will remove MC noise
//  // If no label found, reject cell
//  // It can happen to have this case (4 MC labels per cell is not enough for some cases)
//  // better to remove. To be treated carefully.
//  //
//  if ( nLabels == 0 ) amp = 0;
  
}

///
/// For a given CaloCluster recalculates the position for a given set of misalignment shifts and puts it again in the CaloCluster.
///
/// \param geom: EMCal geometry pointer
/// \param cells: list of EMCal cells with signal
/// \param clu: EMCal cluster subject to position recalculation
///
//______________________________________________________________________________
void AliEMCALRecoUtils::RecalculateClusterPosition(const AliEMCALGeometry *geom, 
                                                   AliVCaloCells* cells, 
                                                   AliVCluster* clu)
{  
  if (!clu)
  {
    AliInfo("Cluster pointer null!");
    return;
  }
    
  if      (fPosAlgo==kPosTowerGlobal) RecalculateClusterPositionFromTowerGlobal( geom, cells, clu);
  else if (fPosAlgo==kPosTowerIndex)  RecalculateClusterPositionFromTowerIndex ( geom, cells, clu);
  else    AliDebug(2,"Algorithm to recalculate position not selected, do nothing.");
}  

///
/// For a given CaloCluster recalculates the position for a given set of misalignment shifts and puts it again in the CaloCluster.
/// The algorithm is a copy of what is done in AliEMCALRecPoint.
///
/// \param geom: EMCal geometry pointer
/// \param cells: list of EMCal cells with signal
/// \param clu: EMCal cluster subject to position recalculation
///
//_____________________________________________________________________________________________
void AliEMCALRecoUtils::RecalculateClusterPositionFromTowerGlobal(const AliEMCALGeometry *geom, 
                                                                  AliVCaloCells* cells, 
                                                                  AliVCluster* clu)
{  
  Double_t eCell       = 0.;
  Float_t  fraction    = 1.;
  Float_t  recalFactor = 1.;
  
  Int_t    absId   = -1;
  Int_t    iTower  = -1, iIphi  = -1, iIeta  = -1;
  Int_t    iSupModMax = -1, iSM=-1, iphi   = -1, ieta   = -1;
  Float_t  weight = 0.,  totalWeight=0.;
  Float_t  newPos[3] = {-1.,-1.,-1.};
  Double_t pLocal[3], pGlobal[3];
  Bool_t shared = kFALSE;

  Float_t  clEnergy = clu->E(); //Energy already recalibrated previously
  if (clEnergy <= 0) return;
  
  GetMaxEnergyCell(geom, cells, clu, absId,  iSupModMax, ieta, iphi,shared);
  Double_t depth = GetDepth(clEnergy,fParticleType,iSupModMax) ;
  
  //printf("** Cluster energy %f, ncells %d, depth %f\n",clEnergy,clu->GetNCells(),depth);
  
  for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++) 
  {
    absId = clu->GetCellAbsId(iDig);
    fraction  = clu->GetCellAmplitudeFraction(iDig);
    if (fraction < 1e-4) fraction = 1.; // in case unfolding is off
    
    if (!fCellsRecalibrated) 
    {
      geom->GetCellIndex(absId,iSM,iTower,iIphi,iIeta); 
      geom->GetCellPhiEtaIndexInSModule(iSM,iTower,iIphi, iIeta,iphi,ieta);      
      if (IsRecalibrationOn()) {
        if(fUse1Drecalib)
          recalFactor = GetEMCALChannelRecalibrationFactor1D(absId);
        else
          recalFactor = GetEMCALChannelRecalibrationFactor(iSM,ieta,iphi);
      }
    }
    
    eCell  = cells->GetCellAmplitude(absId)*fraction*recalFactor;
    
    weight = GetCellWeight(eCell,clEnergy);
    totalWeight += weight;
    
    geom->RelPosCellInSModule(absId,depth,pLocal[0],pLocal[1],pLocal[2]);
    //printf("pLocal (%f,%f,%f), SM %d, absId %d\n",pLocal[0],pLocal[1],pLocal[2],iSupModMax,absId);
    geom->GetGlobal(pLocal,pGlobal,iSupModMax);
    //printf("pLocal (%f,%f,%f)\n",pGlobal[0],pGlobal[1],pGlobal[2]);

    for (int i=0; i<3; i++ ) newPos[i] += (weight*pGlobal[i]);
  }// cell loop
  
  if (totalWeight>0) 
  {
    for (int i=0; i<3; i++ )    newPos[i] /= totalWeight;
  }
    
  //Float_t pos[]={0,0,0};
  //clu->GetPosition(pos);
  //printf("OldPos  : %2.3f,%2.3f,%2.3f\n",pos[0],pos[1],pos[2]);
  //printf("NewPos  : %2.3f,%2.3f,%2.3f\n",newPos[0],newPos[1],newPos[2]);
  
  if (iSupModMax > 1) { //sector 1
    newPos[0] +=fMisalTransShift[3];//-=3.093; 
    newPos[1] +=fMisalTransShift[4];//+=6.82;
    newPos[2] +=fMisalTransShift[5];//+=1.635;
    //printf("   +    : %2.3f,%2.3f,%2.3f\n",fMisalTransShift[3],fMisalTransShift[4],fMisalTransShift[5]);
  } else { //sector 0
    newPos[0] +=fMisalTransShift[0];//+=1.134;
    newPos[1] +=fMisalTransShift[1];//+=8.2;
    newPos[2] +=fMisalTransShift[2];//+=1.197;
    //printf("   +    : %2.3f,%2.3f,%2.3f\n",fMisalTransShift[0],fMisalTransShift[1],fMisalTransShift[2]);
  }
  //printf("NewPos : %2.3f,%2.3f,%2.3f\n",newPos[0],newPos[1],newPos[2]);

  clu->SetPosition(newPos);
}  

///
/// For a given CaloCluster recalculates the position for a given set of misalignment shifts and puts it again in the CaloCluster.
/// The algorithm works with the tower indeces, averages the indeces and from them it calculates the global position.
///
/// \param geom: EMCal geometry pointer
/// \param cells: list of EMCal cells with signal
/// \param clu: EMCal cluster subject to position recalculation
///
//____________________________________________________________________________________________
void AliEMCALRecoUtils::RecalculateClusterPositionFromTowerIndex(const AliEMCALGeometry *geom, 
                                                                 AliVCaloCells* cells, 
                                                                 AliVCluster* clu)
{
  Double_t eCell       = 1.;
  Float_t  fraction    = 1.;
  Float_t  recalFactor = 1.;
  
  Int_t absId   = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1, iIeta   = -1;
  Int_t iSupMod = -1, iSupModMax = -1;
  Int_t iphi = -1, ieta =-1;
  Bool_t shared = kFALSE;

  Float_t clEnergy = clu->E(); //Energy already recalibrated previously.
  
  if (clEnergy <= 0)
    return;
  GetMaxEnergyCell(geom, cells, clu, absId,  iSupModMax, ieta, iphi,shared);
  Float_t  depth = GetDepth(clEnergy,fParticleType,iSupMod) ;

  Float_t weight = 0., weightedCol = 0., weightedRow = 0., totalWeight=0.;
  Bool_t areInSameSM = kTRUE; //exclude clusters with cells in different SMs for now
  Int_t startingSM = -1;
  
  for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++) 
  {
    absId = clu->GetCellAbsId(iDig);
    fraction  = clu->GetCellAmplitudeFraction(iDig);
    if (fraction < 1e-4) fraction = 1.; // in case unfolding is off

    if      (iDig==0)  startingSM = iSupMod;
    else if (iSupMod != startingSM) areInSameSM = kFALSE;

    eCell  = cells->GetCellAmplitude(absId);
    
    geom->GetCellIndex(absId,iSupMod,iTower,iIphi,iIeta); 
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);    
    
    if (!fCellsRecalibrated)
    {
      if (IsRecalibrationOn()) {
        if(fUse1Drecalib)
          recalFactor = GetEMCALChannelRecalibrationFactor1D(absId);
        else
          recalFactor = GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
      }
    }
    
    eCell  = cells->GetCellAmplitude(absId)*fraction*recalFactor;
    
    weight = GetCellWeight(eCell,clEnergy);
    if (weight < 0) weight = 0;
    totalWeight += weight;
    weightedCol += ieta*weight;
    weightedRow += iphi*weight;
    
    //printf("Max cell? cell %d, amplitude org %f, fraction %f, recalibration %f, amplitude new %f \n",cellAbsId, cells->GetCellAmplitude(cellAbsId), fraction, recalFactor, eCell) ;
  }// cell loop
    
  Float_t xyzNew[]={-1.,-1.,-1.};
  if (areInSameSM == kTRUE) 
  {
    //printf("In Same SM\n");
    weightedCol = weightedCol/totalWeight;
    weightedRow = weightedRow/totalWeight;
    geom->RecalculateTowerPosition(weightedRow, weightedCol, iSupModMax, depth, fMisalTransShift, fMisalRotShift, xyzNew); 
  } 
  else 
  {
    //printf("In Different SM\n");
    geom->RecalculateTowerPosition(iphi,        ieta,        iSupModMax, depth, fMisalTransShift, fMisalRotShift, xyzNew); 
  }
  
  clu->SetPosition(xyzNew);
}

///
/// Re-evaluate distance to bad channel with updated bad map
///
/// \param geom: EMCal geometry pointer
/// \param cells: list of EMCal cells with signal
/// \param cluster: EMCal cluster subject to distance to bad channel recalculation
///
//___________________________________________________________________________________________
void AliEMCALRecoUtils::RecalculateClusterDistanceToBadChannel(const AliEMCALGeometry * geom, 
                                                               AliVCaloCells* cells, 
                                                               AliVCluster * cluster)
{             
  if (!fRecalDistToBadChannels) return;
  
  if (!cluster)
  {
    AliInfo("Cluster pointer null!");
    return;
  }  
  
  // Get channels map of the supermodule where the cluster is.
  Int_t absIdMax  = -1, iSupMod =-1, icolM = -1, irowM = -1;
  Bool_t shared = kFALSE;
  GetMaxEnergyCell(geom, cells, cluster, absIdMax,  iSupMod, icolM, irowM, shared);
  TH2D* hMap  = 0x0;
  TH1C* hMap1D = 0x0;

  if(!fUse1Dmap) 
    hMap  = (TH2D*)fEMCALBadChannelMap->At(iSupMod);
  else
    hMap1D  = (TH1C*)fEMCALBadChannelMap->At(0);
    
  Int_t dRrow, dRcol;  
  Float_t  minDist = 10000.;
  Float_t  dist    = 0.;
  
  // Loop on tower status map 
  for (Int_t irow = 0; irow < AliEMCALGeoParams::fgkEMCALRows; irow++)
  {
    for (Int_t icol = 0; icol < AliEMCALGeoParams::fgkEMCALCols; icol++)
    {
      // Check if tower is bad.
      Int_t status=0;
      if (fUse1Dmap)
        status = hMap1D->GetBinContent(geom->GetAbsCellIdFromCellIndexes(iSupMod, irow, icol));
      else
        status = hMap->GetBinContent(icol,irow);

      if(status==0) continue;
      
      //printf("AliEMCALRecoUtils::RecalculateDistanceToBadChannels() - \n \t Bad channel in SM %d, col %d, row %d, \n \t Cluster max in col %d, row %d\n",
      //       iSupMod,icol, irow, icolM,irowM);
      
      dRrow=TMath::Abs(irowM-irow);
      dRcol=TMath::Abs(icolM-icol);
      dist=TMath::Sqrt(dRrow*dRrow+dRcol*dRcol);
      if (dist < minDist)
      {
        //printf("MIN DISTANCE TO BAD %2.2f\n",dist);
        minDist = dist;
      }
    }
  }
  
  // In case the cluster is shared by 2 SuperModules, need to check the map of the second Super Module
  if (shared) 
  {
    TH2D* hMap2 = 0;
    TH1C* hMap1D2 = 0;
    Int_t iSupMod2 = -1;
    
    // The only possible combinations are (0,1), (2,3) ... (8,9)
    if (iSupMod%2) iSupMod2 = iSupMod-1;
    else           iSupMod2 = iSupMod+1;

    if(!fUse1Dmap) 
      hMap2  = (TH2D*)fEMCALBadChannelMap->At(iSupMod2);
    else
      hMap1D2  = (TH1C*)fEMCALBadChannelMap->At(0);
      
    
    // Loop on tower status map of second super module
    for (Int_t irow = 0; irow < AliEMCALGeoParams::fgkEMCALRows; irow++)
    {
      for (Int_t icol = 0; icol < AliEMCALGeoParams::fgkEMCALCols; icol++)
      {
        // Check if tower is bad.
        Int_t status=0;
        if (fUse1Dmap)
          status = hMap1D2->GetBinContent(geom->GetAbsCellIdFromCellIndexes(iSupMod2, irow, icol));
        else
          status = hMap2->GetBinContent(icol,irow);

        if(status==0) continue;
        
        //printf("AliEMCALRecoUtils::RecalculateDistanceToBadChannels(shared) - \n \t Bad channel in SM %d, col %d, row %d \n \t Cluster max in SM %d, col %d, row %d\n",
        //     iSupMod2,icol, irow,iSupMod,icolM,irowM);
        dRrow=TMath::Abs(irow-irowM);
        
        if (iSupMod%2) 
          dRcol=TMath::Abs(icol-(AliEMCALGeoParams::fgkEMCALCols+icolM));
        else 
          dRcol=TMath::Abs(AliEMCALGeoParams::fgkEMCALCols+icol-icolM);
        
        dist=TMath::Sqrt(dRrow*dRrow+dRcol*dRcol);
        if (dist < minDist) minDist = dist;        
      }
    }
  }// shared cluster in 2 SuperModules
  
  AliDebug(2,Form("Max cluster cell (SM,col,row)=(%d %d %d) - Distance to Bad Channel %2.2f",iSupMod, icolM, irowM, minDist));
  cluster->SetDistanceToBadChannel(minDist);
}

///
/// Re-evaluate identification parameters with bayesian
///
/// \param cluster: EMCal cluster subject to PID recalculation
///
///
//__________________________________________________________________
void AliEMCALRecoUtils::RecalculateClusterPID(AliVCluster * cluster)
{             
  if (!cluster)
  {
    AliInfo("Cluster pointer null!");
    return;
  }
  
  if (cluster->GetM02() != 0)
    fPIDUtils->ComputePID(cluster->E(),cluster->GetM02());
  
  Float_t pidlist[AliPID::kSPECIESCN+1];
  for (Int_t i = 0; i < AliPID::kSPECIESCN+1; i++) pidlist[i] = fPIDUtils->GetPIDFinal(i);
        
  cluster->SetPID(pidlist);
}

///
/// Calculates different types of shower shape parameters, dispersion, shower shape eigenvalues and other.
/// 
/// \param geom: EMCal geometry pointer
/// \param cells: list of EMCal cells with signal
/// \param cluster: EMCal cluster subject to shower shape recalculation
/// \param cellEcut: minimum cell energy to be considered in the shower shape recalculation
/// \param cellTimeCut: time window of cells to be considered in shower recalculation
/// \param bc: event bunch crossing number
/// \param enAfterCuts: cluster energy when applying the cell cuts cellEcut and cellTime cut
/// \param l0: main shower shape eigen value
/// \param l1: second eigenvalue of shower shape
/// \param disp: dispersion
/// \param dEta: dispersion in eta (cols) direction
/// \param dPhi: disperion in phi (rows) direction
/// \param sEta: shower shape in eta  (cols) direction
/// \param sPhi: shower shape in phi (rows) direction
/// \param sEtaPhi: shower shape on phi / eta directions term
///
///
//___________________________________________________________________________________________________________________
void AliEMCALRecoUtils::RecalculateClusterShowerShapeParametersWithCellCuts(const AliEMCALGeometry * geom, 
                                                                            AliVCaloCells* cells, AliVCluster * cluster,
                                                                            Float_t cellEcut, Float_t cellTimeCut, Int_t bc, 
                                                                            Float_t & enAfterCuts, Float_t & l0,   Float_t & l1,   
                                                                            Float_t & disp, Float_t & dEta, Float_t & dPhi,
                                                                            Float_t & sEta, Float_t & sPhi, Float_t & sEtaPhi)
{  
  if (!cluster) 
  {
    AliInfo("Cluster pointer null!");
    return;
  }
    
  Double_t eCell       = 0.;
  Double_t tCell       = 0.;
  Float_t  fraction    = 1.;
  Float_t  recalFactor = 1.;
  Bool_t   isLowGain   = kFALSE;

  Int_t    iSupMod = -1;
  Int_t    iTower  = -1;
  Int_t    iIphi   = -1;
  Int_t    iIeta   = -1;
  Int_t    iphi    = -1;
  Int_t    ieta    = -1;
  Double_t etai    = -1.;
  Double_t phii    = -1.;
  
  Int_t    nstat   = 0 ;
  Float_t  wtot    = 0.;
  Double_t w       = 0.;
  Double_t etaMean = 0.;
  Double_t phiMean = 0.;
  
  Double_t pGlobal[3];
  
  // Loop on cells, calculate the cluster energy, in case a cut on cell energy is added,
  // or the non linearity correction was applied
  // and to check if the cluster is between 2 SM in eta
  Int_t   iSM0   = -1;
  Bool_t  shared = kFALSE;
  Float_t energy = 0;
  
  enAfterCuts = 0;
  l0 = 0;  l1 = 0;   
  disp = 0; dEta = 0; dPhi = 0;
  sEta = 0; sPhi = 0; sEtaPhi = 0;
  
  for (Int_t iDigit=0; iDigit < cluster->GetNCells(); iDigit++)
  {
    // Get from the absid the supermodule, tower and eta/phi numbers
    
    Int_t absId = cluster->GetCellAbsId(iDigit);

    geom->GetCellIndex(absId,iSupMod,iTower,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
    
    // Check if there are cells of different SM
    if      (iDigit == 0   ) iSM0 = iSupMod;
    else if (iSupMod!= iSM0) shared = kTRUE;
    
    // Get the cell energy, if recalibration is on, apply factors
    fraction  = cluster->GetCellAmplitudeFraction(iDigit);
    if (fraction < 1e-4) fraction = 1.; // in case unfolding is off
    
    if (IsRecalibrationOn()){
      if(fUse1Drecalib)
        recalFactor = GetEMCALChannelRecalibrationFactor1D(absId);
      else
        recalFactor = GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
    }
    
    eCell  = cells->GetCellAmplitude(absId)*fraction*recalFactor;
    tCell  = cells->GetCellTime     (absId);
    isLowGain = !(cells->GetCellHighGain(absId));//HG = false -> LG = true

    RecalibrateCellTime(absId, bc, tCell,isLowGain);
    tCell*=1e9;
    tCell-=fConstantTimeShift;

    if(eCell > 0.05) // at least a minimum 50 MeV cut
      energy += eCell;

    if(eCell > cellEcut && TMath::Abs(tCell) < cellTimeCut)
      enAfterCuts += eCell;
  } // cell loop
    
  
  // Loop on cells to calculate weights and shower shape terms parameters
  for (Int_t iDigit=0; iDigit < cluster->GetNCells(); iDigit++) 
  {
    // Get from the absid the supermodule, tower and eta/phi numbers
    Int_t absId = cluster->GetCellAbsId(iDigit);

    geom->GetCellIndex(absId,iSupMod,iTower,iIphi,iIeta); 
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);        
    
    //Get the cell energy, if recalibration is on, apply factors
    fraction  = cluster->GetCellAmplitudeFraction(iDigit);
    if (fraction < 1e-4) fraction = 1.; // in case unfolding is off
    
    if (!fCellsRecalibrated && IsRecalibrationOn()){
        if(fUse1Drecalib)
          recalFactor = GetEMCALChannelRecalibrationFactor1D(absId);
        else
          recalFactor = GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
    }
    
    eCell  = cells->GetCellAmplitude(absId)*fraction*recalFactor;
    tCell  = cells->GetCellTime     (absId);
    isLowGain = !(cells->GetCellHighGain(absId));//HG = false -> LG = true

    RecalibrateCellTime(absId, bc, tCell,isLowGain);
    tCell*=1e9;
    tCell-=fConstantTimeShift;
    
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2
    // C Side impair SM, nSupMod%2=1; A side pair SM, nSupMod%2=0
    if (shared && iSupMod%2) ieta+=AliEMCALGeoParams::fgkEMCALCols;
    
    if ( energy > 0 && eCell > cellEcut && TMath::Abs(tCell) < cellTimeCut ) 
    {
      w  = GetCellWeight(eCell, energy);
      
      // Cell index
      if     ( fShowerShapeCellLocationType == 0 )
      {
        etai=(Double_t)ieta;
        phii=(Double_t)iphi;  
      }
      // Cell angle location
      else if( fShowerShapeCellLocationType == 1 )
      {
        geom->EtaPhiFromIndex(absId, etai, phii);
        etai *= TMath::RadToDeg(); // change units to degrees instead of radians
        phii *= TMath::RadToDeg(); // change units to degrees instead of radians       
      }
      else
      {
        geom->GetGlobal(absId,pGlobal);
        
        // Cell x-z location
        if( fShowerShapeCellLocationType == 2 )
        {
          etai = pGlobal[2];
          phii = pGlobal[0];
        }
        // Cell r-z location
        else
        {
          etai = pGlobal[2];
          phii = TMath::Sqrt(pGlobal[0]*pGlobal[0]+pGlobal[1]*pGlobal[1]); 
        }
      }
               
      if (w > 0.0) 
      {
        wtot += w ;
        nstat++;   
        
        // Shower shape
        sEta     += w * etai * etai ;
        etaMean  += w * etai ;
        sPhi     += w * phii * phii ;
        phiMean  += w * phii ; 
        sEtaPhi  += w * etai * phii ; 
      }
    } 
    else if(eCell > 0.05)
      AliDebug(2,Form("Wrong energy in cell %f and/or cluster %f\n", eCell, cluster->E()));
  } // cell loop
  
  //printf("sEta %f sPhi %f etaMean %f phiMean %f sEtaPhi %f wtot %f\n",sEta,sPhi,etaMean,phiMean,sEtaPhi, wtot);
  
  // Normalize to the weight  
  if (wtot > 0) 
  {
    etaMean /= wtot ;
    phiMean /= wtot ;
  } 
  else
    AliDebug(2,Form("Wrong weight %f\n", wtot));
  
  // Loop on cells to calculate dispersion  
  for (Int_t iDigit=0; iDigit < cluster->GetNCells(); iDigit++) 
  {
    // Get from the absid the supermodule, tower and eta/phi numbers
    Int_t absId = cluster->GetCellAbsId(iDigit);

    geom->GetCellIndex(absId,iSupMod,iTower,iIphi,iIeta); 
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
    
    //Get the cell energy, if recalibration is on, apply factors
    fraction  = cluster->GetCellAmplitudeFraction(iDigit);
    if (fraction < 1e-4) fraction = 1.; // in case unfolding is off
    if (IsRecalibrationOn()) {
      if(fUse1Drecalib)
        recalFactor = GetEMCALChannelRecalibrationFactor1D(absId);
      else
        recalFactor = GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
    }
    
    eCell  = cells->GetCellAmplitude(absId)*fraction*recalFactor;
    tCell  = cells->GetCellTime     (absId);
    isLowGain = !(cells->GetCellHighGain(absId));//HG = false -> LG = true

    RecalibrateCellTime(absId, bc, tCell,isLowGain);
    tCell*=1e9;
    tCell-=fConstantTimeShift;

    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2
    // C Side impair SM, nSupMod%2=1; A side pair SM, nSupMod%2=0
    if (shared && iSupMod%2) ieta+=AliEMCALGeoParams::fgkEMCALCols;
    
    if ( energy > 0 && eCell > cellEcut && TMath::Abs(tCell) < cellTimeCut ) 
    {
      w  = GetCellWeight(eCell,cluster->E());
      
      // Cell index
      if     ( fShowerShapeCellLocationType == 0 )
      {
        etai=(Double_t)ieta;
        phii=(Double_t)iphi;  
      }
      // Cell angle location
      else if( fShowerShapeCellLocationType == 1 )
      {
        geom->EtaPhiFromIndex(absId, etai, phii);
        etai *= TMath::RadToDeg(); // change units to degrees instead of radians
        phii *= TMath::RadToDeg(); // change units to degrees instead of radians       
      }
      else
      {
        geom->GetGlobal(absId,pGlobal);
        
        // Cell x-z location
        if( fShowerShapeCellLocationType == 2 )
        {
          etai = pGlobal[2];
          phii = pGlobal[0];
        }
        // Cell r-z location
        else
        {
          etai = pGlobal[2];
          phii = TMath::Sqrt(pGlobal[0]*pGlobal[0]+pGlobal[1]*pGlobal[1]); 
        }
      }
      
      if (w > 0.0) 
      { 
        disp +=  w *((etai-etaMean)*(etai-etaMean)+(phii-phiMean)*(phii-phiMean)); 
        dEta +=  w * (etai-etaMean)*(etai-etaMean) ; 
        dPhi +=  w * (phii-phiMean)*(phii-phiMean) ; 
      }
    } 
    else
      AliDebug(2,Form("Wrong energy in cell %f and/or cluster %f\n", eCell, cluster->E()));
  }// cell loop
  
  // Normalize to the weigth and set shower shape parameters
  if (wtot > 0 && nstat > 1)
  {
    disp    /= wtot ;
    dEta    /= wtot ;
    dPhi    /= wtot ;
    sEta    /= wtot ;
    sPhi    /= wtot ;
    sEtaPhi /= wtot ;
    
    sEta    -= etaMean * etaMean ;
    sPhi    -= phiMean * phiMean ;
    sEtaPhi -= etaMean * phiMean ;
    
    l0 = (0.5 * (sEta + sPhi) + TMath::Sqrt( 0.25 * (sEta - sPhi) * (sEta - sPhi) + sEtaPhi * sEtaPhi ));
    l1 = (0.5 * (sEta + sPhi) - TMath::Sqrt( 0.25 * (sEta - sPhi) * (sEta - sPhi) + sEtaPhi * sEtaPhi ));
    
    //printf("sEta %f sPhi %f etaMean %f phiMean %f sEtaPhi %f wtot %f l0 %f l1 %f\n",sEta,sPhi,etaMean,phiMean,sEtaPhi, wtot,l0,l1);

  } 
  else 
  {
    l0   = 0. ;
    l1   = 0. ;
    dEta = 0. ; dPhi = 0. ; disp    = 0. ;
    sEta = 0. ; sPhi = 0. ; sEtaPhi = 0. ;
  }  
}

///
/// Calculates different types of shower shape parameters, dispersion, shower shape eigenvalues and other.
/// Call to AliEMCALRecoUtils::RecalculateClusterShowerShapeParametersWithCellCuts 
/// with default cell cuts (50 MeV minimum cell energy and not cut on time)
/// 
/// \param geom: EMCal geometry pointer
/// \param cells: list of EMCal cells with signal
/// \param cluster: EMCal cluster subject to shower shape recalculation
/// \param l0: main shower shape eigen value
/// \param l1: second eigenvalue of shower shape
/// \param disp: dispersion
/// \param dEta: dispersion in eta (cols) direction
/// \param dPhi: disperion in phi (rows) direction
/// \param sEta: shower shape in eta  (cols) direction
/// \param sPhi: shower shape in phi (rows) direction
/// \param sEtaPhi: shower shape on phi / eta directions term
///
///
//___________________________________________________________________________________________________________________
void AliEMCALRecoUtils::RecalculateClusterShowerShapeParameters(const AliEMCALGeometry * geom, 
                                                                AliVCaloCells* cells, AliVCluster * cluster,
                                                                Float_t & l0,   Float_t & l1,   
                                                                Float_t & disp, Float_t & dEta, Float_t & dPhi,
                                                                Float_t & sEta, Float_t & sPhi, Float_t & sEtaPhi)
{
  Float_t newEnergy      = 0;
  Float_t cellEmin       = 0.05; // 50 MeV
  Float_t cellTimeWindow = 1000; // open cut
  Int_t   bc             = 0; 
  AliEMCALRecoUtils::RecalculateClusterShowerShapeParametersWithCellCuts(geom, cells, cluster,
                                                                         cellEmin, cellTimeWindow, bc,
                                                                         newEnergy, l0, l1, disp,
                                                                         dEta, dPhi, sEta, sPhi, sEtaPhi);
}

///
/// Calculates Dispersion and main axis and puts them into the cluster
/// Call to method RecalculateClusterShowerShapeParameters
///
/// \param geom: EMCal geometry pointer
/// \param cells: list of EMCal cells with signal
/// \param cluster: EMCal cluster subject to shower shape recalculation
///
//____________________________________________________________________________________________
void AliEMCALRecoUtils::RecalculateClusterShowerShapeParameters(const AliEMCALGeometry * geom, 
                                                                AliVCaloCells* cells, 
                                                                AliVCluster * cluster)
{  
  Float_t l0   = 0., l1   = 0.;
  Float_t disp = 0., dEta = 0., dPhi    = 0.; 
  Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;
  
  AliEMCALRecoUtils::RecalculateClusterShowerShapeParameters(geom,cells,cluster,l0,l1,disp,
                                                             dEta, dPhi, sEta, sPhi, sEtaPhi);
  
  cluster->SetM02(l0);
  cluster->SetM20(l1);
  if (disp > 0. ) cluster->SetDispersion(TMath::Sqrt(disp)) ;
  
} 

///
/// Calculates Dispersion and main axis and puts them into the cluster.
/// Possibility to restrict the cell Energy and time window in the calculations
/// 
/// \param geom: EMCal geometry pointer
/// \param cells: list of EMCal cells with signal
/// \param cluster: EMCal cluster subject to shower shape recalculation
/// \param cellEcut: minimum cell energy to be considered in the shower shape recalculation
/// \param cellTimeCut: time window of cells to be considered in shower recalculation
/// \param bc: event bunch crossing number
/// \param enAfterCuts: cluster energy when applying the cell cuts cellEcut and cellTime cut
///
//____________________________________________________________________________________________
void AliEMCALRecoUtils::RecalculateClusterShowerShapeParametersWithCellCuts(const AliEMCALGeometry * geom, 
                                                                            AliVCaloCells* cells, AliVCluster * cluster,
                                                                            Float_t cellEcut, Float_t cellTimeCut, Int_t bc,
                                                                            Float_t & enAfterCuts)
{  
  Float_t l0   = 0., l1   = 0.;
  Float_t disp = 0., dEta = 0., dPhi    = 0.; 
  Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;
  
  AliEMCALRecoUtils::RecalculateClusterShowerShapeParametersWithCellCuts(geom, cells, cluster,
                                                                         cellEcut, cellTimeCut, bc,
                                                                         enAfterCuts, l0, l1, disp,
                                                                         dEta, dPhi, sEta, sPhi, sEtaPhi);
  
  cluster->SetM02(l0);
  cluster->SetM20(l1);
  if (disp > 0. ) cluster->SetDispersion(TMath::Sqrt(disp)) ;
  
} 

///
/// Find the candidate cluster-track matchs.
///
/// This function should be called before the cluster loop.
/// Before call this function, please recalculate the cluster positions.
/// Given the input event, loop over all the tracks, select the closest cluster as matched with fCutR.
/// Store matched cluster indexes and residuals.
///
/// \param event: event pointer
/// \param clusterArr: list of clusters
/// \param geom: AliEMCALGeometry pointer
/// \param mc: AliMCEvent pointer
///
//____________________________________________________________________________
void AliEMCALRecoUtils::FindMatches(AliVEvent *event,
                                    TObjArray * clusterArr,  
                                    const AliEMCALGeometry *geom,
                                    AliMCEvent * mc)
{  
  fMatchedTrackIndex.Reset();
  fMatchedClusterIndex.Reset();
  fResidualPhi.Reset();
  fResidualEta.Reset();
  
  fMatchedTrackIndex.Set(1000);
  fMatchedClusterIndex.Set(1000);
  fResidualPhi.Set(1000);
  fResidualEta.Set(1000);
  
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*> (event);
  AliAODEvent* aodevent = dynamic_cast<AliAODEvent*> (event);
  
  // Init the magnetic field if not already on
  if (!TGeoGlobalMagField::Instance()->GetField()) 
  {
    if (!event->InitMagneticField())
    {
      AliInfo("Mag Field not initialized, null esd/aod evetn pointers");
    }
  } // Init mag field
  
  if (esdevent) 
  {
    UInt_t mask1 = esdevent->GetESDRun()->GetDetectorsInDAQ();
    UInt_t mask2 = esdevent->GetESDRun()->GetDetectorsInReco();
    Bool_t desc1 = (mask1 >> 3) & 0x1;
    Bool_t desc2 = (mask2 >> 3) & 0x1;
    if (desc1==0 || desc2==0) { 
      //       AliError(Form("TPC not in DAQ/RECO: %u (%u)/%u (%u)", 
      //       mask1, esdevent->GetESDRun()->GetDetectorsInReco(),
      //       mask2, esdevent->GetESDRun()->GetDetectorsInDAQ()));
      fITSTrackSA=kTRUE;
    }
  }
  
  TObjArray *clusterArray = 0x0;
  if (!clusterArr) 
  {
    clusterArray = new TObjArray(event->GetNumberOfCaloClusters());
    for (Int_t icl=0; icl<event->GetNumberOfCaloClusters(); icl++) 
    {
      AliVCluster *cluster = (AliVCluster*) event->GetCaloCluster(icl);
      if (!IsGoodCluster(cluster,geom,(AliVCaloCells*)event->GetEMCALCells())) 
        continue;
      clusterArray->AddAt(cluster,icl);
    }
  }
  
  Int_t    matched=0;
  Double_t cv[21];
  TString  genName;
  for (Int_t i=0; i<21;i++) cv[i]=0;
  for (Int_t itr=0; itr<event->GetNumberOfTracks(); itr++)
  {
    AliExternalTrackParam *trackParam = 0;
    Int_t mcLabel = -1;
    // If the input event is ESD, the starting point for extrapolation is TPCOut, if available, or TPCInner 
    AliESDtrack *esdTrack = 0;
    AliAODTrack *aodTrack = 0;
    if (esdevent) 
    {
      esdTrack = esdevent->GetTrack(itr);
      if (!esdTrack) continue;
      if (!IsAccepted(esdTrack)) continue;
      if (esdTrack->Pt()<fCutMinTrackPt) continue;
      
      if ( TMath::Abs(esdTrack->Eta()) > 0.9 ) continue;
      
      // Save some time and memory in case of no DCal present
      if( geom->GetNumberOfSuperModules() < 13 ) // Run1 10 (12, 2 not active but present)
      {
        Double_t phi = esdTrack->Phi()*TMath::RadToDeg();
        if ( phi <= 10 || phi >= 250 ) continue;
      }
      
      if (!fITSTrackSA) // if TPC Available
      {
        if ( fUseOuterTrackParam )
          trackParam =  const_cast<AliExternalTrackParam*>(esdTrack->GetOuterParam());  
        else
          trackParam =  const_cast<AliExternalTrackParam*>(esdTrack->GetInnerParam());  
      }
      else
        trackParam =  new AliExternalTrackParam(*esdTrack); // If ITS Track Standing alone		
      
      mcLabel = TMath::Abs(esdTrack->GetLabel());
    }
    
    // If the input event is AOD, the starting point for extrapolation is at vertex
    // AOD tracks are selected according to its filterbit.
    else if (aodevent) 
    {
      aodTrack = dynamic_cast<AliAODTrack*>(aodevent->GetTrack(itr));
      
      if (!aodTrack) AliFatal("Not a standard AOD");
     
      if (!aodTrack) continue;
      
      if (fAODTPCOnlyTracks) 
      { // Match with TPC only tracks, default from May 2013, before filter bit 32
        //printf("Match with TPC only tracks, accept? %d, test bit 128 <%d> \n", aodTrack->IsTPCOnly(), aodTrack->TestFilterMask(128));
        if (!aodTrack->IsTPCConstrained()) continue ;
      } 
      else if (fAODHybridTracks) 
      { // Match with hybrid tracks
        //printf("Match with Hybrid tracks, accept? %d \n", aodTrack->IsHybridGlobalConstrainedGlobal());
        if (!aodTrack->IsHybridGlobalConstrainedGlobal()) continue ;
      } else 
      { // Match with tracks on a mask
        //printf("Match with tracks having filter bit mask %d, accept? %d \n",fAODFilterMask,aodTrack->TestFilterMask(fAODFilterMask));
        if (!aodTrack->TestFilterMask(fAODFilterMask) ) continue; //Select AOD tracks
      }
      
      if (aodTrack->Pt() < fCutMinTrackPt) continue;
      
      if ( TMath::Abs(aodTrack->Eta()) > 0.9 ) continue;
      
      // Save some time and memory in case of no DCal present
      if( geom->GetNumberOfSuperModules() < 13 ) // Run1 10 (12, 2 not active but present)
      {
        Double_t phi = aodTrack->Phi()*TMath::RadToDeg();
        if ( phi <= 10 || phi >= 250 ) continue;
      }
      
      Double_t pos[3],mom[3];
      
      if ( fUseTrackDCA )
        aodTrack->GetXYZ(pos);
      else
        aodTrack->XvYvZv(pos);
        
      aodTrack->GetPxPyPz(mom);
      AliDebug(5,Form("aod track: i=%d | pos=(%5.4f,%5.4f,%5.4f) | mom=(%5.4f,%5.4f,%5.4f) | charge=%d\n",
                      itr,pos[0],pos[1],pos[2],mom[0],mom[1],mom[2],aodTrack->Charge()));
      
      trackParam= new AliExternalTrackParam(pos,mom,cv,aodTrack->Charge());
      
      mcLabel = TMath::Abs(aodTrack->GetLabel());
    }
    
    //Return if the input data is not "AOD" or "ESD"
    else
    {
      AliWarning("Wrong input data type! Should be \"AOD\" or \"ESD\" ");
      if (clusterArray) 
      {
        clusterArray->Clear();
        delete clusterArray;
      }
      return;
    }
    
    if (!trackParam) continue;
    
    //
    // Check if track comes from a particular MC generator, do not include it
    // if it is not a selected one
    //
    if( mc && fMCGenerToAcceptForTrack && fNMCGenerToAccept > 0 )
    {
      mc->GetCocktailGenerator(mcLabel,genName);
            
      Bool_t generOK = kFALSE;
      for(Int_t ig = 0; ig < fNMCGenerToAccept; ig++) 
      {
        if ( genName.Contains(fMCGenerToAccept[ig]) ) generOK = kTRUE;
      }

      if ( !generOK ) continue;
    }
    
    // Extrapolate the track to EMCal surface, see AliEMCALRecoUtilsBase
    AliExternalTrackParam emcalParam(*trackParam);
    Float_t eta, phi, pt;
    if (!ExtrapolateTrackToEMCalSurface(&emcalParam, fEMCalSurfaceDistance, fMass, fStepSurface, eta, phi, pt)) 
    {
      if (aodevent    && trackParam) delete trackParam;
      if (fITSTrackSA && trackParam) delete trackParam;
      continue;
    }
    
    if ( TMath::Abs(eta) > 0.75 ) 
    {
      if ( trackParam && (aodevent || fITSTrackSA) )   delete trackParam;
      continue;
    }
    
    // Save some time and memory in case of no DCal present
    if ( geom->GetNumberOfSuperModules() < 13 &&  // Run1 10 (12, 2 not active but present)
        ( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad())) 
    {
      if ( trackParam && (aodevent || fITSTrackSA) )   delete trackParam;
      continue;
    }
    
    //Find matched clusters
    Int_t index = -1;
    Float_t dEta = -999, dPhi = -999;
    if (!clusterArr) 
      index = FindMatchedClusterInClusterArr(&emcalParam, &emcalParam, clusterArray, dEta, dPhi);  
    else 
      index = FindMatchedClusterInClusterArr(&emcalParam, &emcalParam, clusterArr  , dEta, dPhi);  
    
    
    if (index>-1) 
    {
      fMatchedTrackIndex[itr] = matched;
      fMatchedClusterIndex[index] = matched;
      fResidualEta[dEta] = matched;
      fResidualPhi[dPhi] = matched;
      matched++;
    }
    
    if (aodevent && trackParam) delete trackParam;
    if (fITSTrackSA && trackParam) delete trackParam;
  }//track loop
  
  if (clusterArray) 
  {
    clusterArray->Clear();
    delete clusterArray;
  }
  
  AliDebug(2,Form("Number of matched pairs = %d !\n",matched));
  
  fMatchedTrackIndex.Set(matched);
  fMatchedClusterIndex.Set(matched);
  fResidualPhi.Set(matched);
  fResidualEta.Set(matched);
}

///
/// Find matched cluster in event. See Find MatchedClusterInClusterArr().
///
/// \param track: track pointer
/// \param event: event pointer
/// \param geom: AliEMCALGeometry pointer
/// \param dEta: found track-cluster match residual in eta direction
/// \param dPhi: found track-cluster match residual in phi direction
///
/// \return  the index of matched cluster to input track, returns -1 if no match is found
///
//________________________________________________________________________________
Int_t AliEMCALRecoUtils::FindMatchedClusterInEvent(const AliESDtrack *track, 
                                                   const AliVEvent *event, 
                                                   const AliEMCALGeometry *geom, 
                                                   Float_t &dEta, Float_t &dPhi)
{
  Int_t index = -1;
  
  if ( TMath::Abs(track->Eta()) > 0.9 ) return index;
  
  // Save some time and memory in case of no DCal present
  if( geom->GetNumberOfSuperModules() < 13 ) // Run1 10 (12, 2 not active but present)
  {
    Double_t phiV = track->Phi()*TMath::RadToDeg();
    if ( phiV <= 10 || phiV >= 250 ) return index;
  }
  
  AliExternalTrackParam *trackParam = 0;
  if (!fITSTrackSA) // If TPC
  {
    if ( fUseOuterTrackParam )
      trackParam = const_cast<AliExternalTrackParam*>(track->GetOuterParam());  
    else 
      trackParam = const_cast<AliExternalTrackParam*>(track->GetInnerParam());  
  }
  else
    trackParam = new AliExternalTrackParam(*track);
  
  if (!trackParam) return index;
  AliExternalTrackParam emcalParam(*trackParam);
  
  Float_t eta, phi, pt;  
  if (!AliEMCALRecoUtilsBase::ExtrapolateTrackToEMCalSurface(&emcalParam, fEMCalSurfaceDistance, fMass, fStepSurface, eta, phi, pt))	
  {
    if (fITSTrackSA) delete trackParam;
    return index;
  }
  
  if ( TMath::Abs(eta) > 0.75 ) 
  {
    if (fITSTrackSA) delete trackParam;
    return index;
  }
  
  // Save some time and memory in case of no DCal present
  if ( geom->GetNumberOfSuperModules() < 13 &&  // Run1 10 (12, 2 not active but present)
      ( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad())) 
  {
    if (fITSTrackSA) delete trackParam;
    return index;    
  }
  
  TObjArray *clusterArr = new TObjArray(event->GetNumberOfCaloClusters());
  
  for (Int_t icl=0; icl<event->GetNumberOfCaloClusters(); icl++)
  {
    AliVCluster *cluster = (AliVCluster*) event->GetCaloCluster(icl);
    if (!IsGoodCluster(cluster,geom,(AliVCaloCells*)event->GetEMCALCells())) continue;
    clusterArr->AddAt(cluster,icl);
  }
  
  index = FindMatchedClusterInClusterArr(&emcalParam, &emcalParam, clusterArr, dEta, dPhi);  
  clusterArr->Clear();
  delete clusterArr;
  if (fITSTrackSA) delete trackParam;
  
  return index;
}

///
/// Find matched cluster in input array of clusters.
/// 
/// \param emcalParam: emcal track parameters container?
/// \param trkParam: track parameters container?
/// \param clusterArr: input array of clusters
/// \param dEta: found track-cluster match residual in eta direction
/// \param dPhi: found track-cluster match residual in phi direction
///
/// \return  the index of matched cluster to input track.
//_______________________________________________________________________________________________
Int_t  AliEMCALRecoUtils::FindMatchedClusterInClusterArr(const AliExternalTrackParam *emcalParam, 
                                                         AliExternalTrackParam *trkParam, 
                                                         const TObjArray * clusterArr, 
                                                         Float_t &dEta, Float_t &dPhi)
{  
  dEta=-999, dPhi=-999;
  Float_t dRMax = fCutR, dEtaMax=fCutEta, dPhiMax=fCutPhi;
  Int_t index = -1;
  Float_t tmpEta=-999, tmpPhi=-999;

  Double_t exPos[3] = {0.,0.,0.};
  if (!emcalParam->GetXYZ(exPos)) return index;

  Float_t clsPos[3] = {0.,0.,0.};
  for (Int_t icl=0; icl<clusterArr->GetEntriesFast(); icl++)
  {
    AliVCluster *cluster = dynamic_cast<AliVCluster*> (clusterArr->At(icl)) ;
    
    if (!cluster || !cluster->IsEMCAL()) continue;
    
    cluster->GetPosition(clsPos);
    
    Double_t dR = TMath::Sqrt(TMath::Power(exPos[0]-clsPos[0],2)+TMath::Power(exPos[1]-clsPos[1],2)+TMath::Power(exPos[2]-clsPos[2],2));
    if (dR > fClusterWindow) continue;
    
    AliExternalTrackParam trkPamTmp (*trkParam);//Retrieve the starting point every time before the extrapolation
    
    if (!AliEMCALRecoUtilsBase::ExtrapolateTrackToCluster(&trkPamTmp, cluster, fMass, fStepCluster, tmpEta, tmpPhi)) continue;
    
    if (fCutEtaPhiSum) 
    {
      Float_t tmpR=TMath::Sqrt(tmpEta*tmpEta + tmpPhi*tmpPhi);
      if (tmpR<dRMax) 
      {
        dRMax=tmpR;
        dEtaMax=tmpEta;
        dPhiMax=tmpPhi;
        index=icl;
      }
    } 
    else if (fCutEtaPhiSeparate) 
    {
      if (TMath::Abs(tmpEta)<TMath::Abs(dEtaMax) && TMath::Abs(tmpPhi)<TMath::Abs(dPhiMax)) 
      {
        dEtaMax = tmpEta;
        dPhiMax = tmpPhi;
        index=icl;
      }
    } 
    else 
    {
      AliError("Please specify your cut criteria");
      AliError("To cut on sqrt(dEta^2+dPhi^2), use: SwitchOnCutEtaPhiSum()");
      AliError("To cut on dEta and dPhi separately, use: SwitchOnCutEtaPhiSeparate()");
      return index;
    }
  }

  dEta=dEtaMax;
  dPhi=dPhiMax;

  return index;
}

///
/// Return the residual by extrapolating a track param to a cluster.
/// Mass and step hypothesis are set via data members fStepCluster and fMass 
/// passed to the main method in AliEMCALRecoUtilsBase
///
/// \param trkParam: pointer to external track param
/// \param cluster: pointer to AliVCluster
/// \param tmpEta: residual eta
/// \param tmpPhi: residual phi
///
/// \return bool true if track could be extrapolated
///
//---------------------------------------------------------------------------------
Bool_t AliEMCALRecoUtils::ExtrapolateTrackToCluster(AliExternalTrackParam *trkParam, 
                                                    const AliVCluster *cluster, 
                                                    Float_t &tmpEta, 
                                                    Float_t &tmpPhi)
{
  return AliEMCALRecoUtilsBase::ExtrapolateTrackToCluster(trkParam, cluster, fMass, fStepCluster, tmpEta, tmpPhi);
}

///
/// Given a cluster index, get the residuals dEta and dPhi for this cluster 
/// to the closest track.
/// It works with ESDs and AODs.
///
/// \param clsIndex: cluster index as in AliESDEvent::GetCaloCluster(clsIndex)
/// \param dEta: residual eta
/// \param dPhi: residual phi
///
//_______________________________________________________________________
void AliEMCALRecoUtils::GetMatchedResiduals(Int_t clsIndex, 
                                            Float_t &dEta, Float_t &dPhi)
{
  if (FindMatchedPosForCluster(clsIndex) >= 999) 
  {
    AliDebug(2,"No matched tracks found!\n");
    dEta=999.;
    dPhi=999.;
    return;
  }
  
  dEta = fResidualEta[FindMatchedPosForCluster(clsIndex)];
  dPhi = fResidualPhi[FindMatchedPosForCluster(clsIndex)];
}

///
/// Given a track index, get the residuals dEta and dPhi for this track 
/// to the closest cluster.
/// It works with ESDs and AODs.
///
/// \param trkIndex: cluster index as in AliESDEvent::GetTrack(trkIndex)
/// \param dEta: residual eta
/// \param dPhi: residual phi
///
//______________________________________________________________________________________________
void AliEMCALRecoUtils::GetMatchedClusterResiduals(Int_t trkIndex, Float_t &dEta, Float_t &dPhi)
{
  if (FindMatchedPosForTrack(trkIndex) >= 999) 
  {
    AliDebug(2,"No matched cluster found!\n");
    dEta=999.;
    dPhi=999.;
    return;
  }
  
  dEta = fResidualEta[FindMatchedPosForTrack(trkIndex)];
  dPhi = fResidualPhi[FindMatchedPosForTrack(trkIndex)];
}

///
/// Given a cluster index , get the index of matched track to this cluster.
/// It works with ESDs and AODs.
///
/// \param clsIndex: cluster index as in AliESDEvent::GetCaloCluster(clsIndex)
///
//__________________________________________________________
Int_t AliEMCALRecoUtils::GetMatchedTrackIndex(Int_t clsIndex)
{
  if (IsClusterMatched(clsIndex))
    return fMatchedTrackIndex[FindMatchedPosForCluster(clsIndex)];
  else 
    return -1; 
}

///
/// Given a track index, get the index of matched cluster to this track.
/// It works with ESDs and AODs.
///
/// \param trkIndex: cluster index as in AliESDEvent::GetTrack(trkIndex)
///
//__________________________________________________________
Int_t AliEMCALRecoUtils::GetMatchedClusterIndex(Int_t trkIndex)
{  
  if (IsTrackMatched(trkIndex))
    return fMatchedClusterIndex[FindMatchedPosForTrack(trkIndex)];
  else 
    return -1; 
}

///
/// Given a cluster index, it returns if the cluster has a match.
///
/// \param clsIndex: cluster index as in AliESDEvent::GetCaloCluster(clsIndex)
///
/// \return bool true if cluster is matched
///
//______________________________________________________________
Bool_t AliEMCALRecoUtils::IsClusterMatched(Int_t clsIndex) const
{
  if (FindMatchedPosForCluster(clsIndex) < 999) 
    return kTRUE;
  else
    return kFALSE;
}

///
/// Given a track index, it returns if the track has a match.
///
/// \param trkIndex: track index as in AliESDEvent::GetTrack(trkIndex)
///
/// \return bool true if cluster is matched
///
//____________________________________________________________
Bool_t AliEMCALRecoUtils::IsTrackMatched(Int_t trkIndex) const 
{
  if (FindMatchedPosForTrack(trkIndex) < 999) 
    return kTRUE;
  else
    return kFALSE;
}

///
/// Given a cluster index, it returns the position of the match in the fMatchedClusterIndex array
///
/// \param clsIndex: cluster index as in AliESDEvent::GetCaloCluster(clsIndex)
///
/// \return cluster position index
///
//______________________________________________________________________
UInt_t AliEMCALRecoUtils::FindMatchedPosForCluster(Int_t clsIndex) const
{
  Float_t tmpR = fCutR;
  UInt_t pos = 999;
  
  for (Int_t i=0; i<fMatchedClusterIndex.GetSize(); i++) 
  {
    if (fMatchedClusterIndex[i]==clsIndex) 
    {
      Float_t r = TMath::Sqrt(fResidualEta[i]*fResidualEta[i] + fResidualPhi[i]*fResidualPhi[i]);
      if (r<tmpR) 
      {
        pos=i;
        tmpR=r;
        AliDebug(3,Form("Matched cluster index: index: %d, dEta: %2.4f, dPhi: %2.4f.\n",
                        fMatchedClusterIndex[i],fResidualEta[i],fResidualPhi[i]));
      }
    }
  }
  return pos;
}

///
/// Given a cluster index, it returns the position of the match in the fMatchedTrackIndex  array
///
/// \param trkIndex: cluster index as in AliESDEvent::GetCaloCluster(clsIndex)
///
/// \return track position index
///
//____________________________________________________________________
UInt_t AliEMCALRecoUtils::FindMatchedPosForTrack(Int_t trkIndex) const
{
  Float_t tmpR = fCutR;
  UInt_t pos = 999;
  
  for (Int_t i=0; i<fMatchedTrackIndex.GetSize(); i++) 
  {
    if (fMatchedTrackIndex[i]==trkIndex) 
    {
      Float_t r = TMath::Sqrt(fResidualEta[i]*fResidualEta[i] + fResidualPhi[i]*fResidualPhi[i]);
      if (r<tmpR) 
      {
        pos=i;
        tmpR=r;
        AliDebug(3,Form("Matched track index: index: %d, dEta: %2.4f, dPhi: %2.4f.\n",
                        fMatchedTrackIndex[i],fResidualEta[i],fResidualPhi[i]));
      }
    }
  }
  return pos;
}

///
/// Check if the cluster survives the following quality cuts:
///   * Cluster pointer is non null and is EMCal (or DCal)
///   * There is no bad channel cell inside
///   * The cluster is not too close to a fiducial border
///   * The cluster is not exotic
///
/// \param cluster: pointer to AliVCluster
/// \param geom: pointer to AliEMCALGeometry
/// \param cells: full list of cells
/// \param bc: event bunch crossing number
///
/// \return true if cluster passes all selection criteria
///
//__________________________________________________________________________
Bool_t AliEMCALRecoUtils::IsGoodCluster(AliVCluster *cluster, 
                                        const AliEMCALGeometry *geom, 
                                        AliVCaloCells* cells, Int_t bc)
{
  Bool_t isGood=kTRUE;

  if (!cluster || !cluster->IsEMCAL())              return kFALSE;
  if (ClusterContainsBadChannel(geom,cluster->GetCellsAbsId(),cluster->GetNCells())) return kFALSE;
  if (!CheckCellFiducialRegion(geom,cluster,cells)) return kFALSE;
  if (IsExoticCluster(cluster, cells,bc))           return kFALSE;

  return isGood;
}

///
/// Given a esd track, return whether the track survive all the cuts.
///
/// The different quality parameter are first.
/// retrieved from the track. then it is found out what cuts the
/// track did not survive and finally the cuts are imposed.
///
/// \param esdTrack: pointer to ESD track
///
/// \return true if track passes all selection criteria
//__________________________________________________________
Bool_t AliEMCALRecoUtils::IsAccepted(AliESDtrack *esdTrack)
{
  UInt_t status = esdTrack->GetStatus();

  Int_t nClustersITS = esdTrack->GetITSclusters(0);
  Int_t nClustersTPC = esdTrack->GetTPCclusters(0);

  Float_t chi2PerClusterITS = -1;
  Float_t chi2PerClusterTPC = -1;
  if (nClustersITS!=0)
    chi2PerClusterITS = esdTrack->GetITSchi2()/Float_t(nClustersITS);
  if (nClustersTPC!=0) 
    chi2PerClusterTPC = esdTrack->GetTPCchi2()/Float_t(nClustersTPC);

  //
  // DCA cuts
  // Only to be used for primary
  //
  if ( fTrackCutsType == kGlobalCut ) 
  {
    Float_t maxDCAToVertexXYPtDep = 0.0182 + 0.0350/TMath::Power(esdTrack->Pt(),1.01); 
    // This expression comes from AliESDtrackCuts::GetStandardITSTPCTrackCuts2010()
    
    //AliDebug(3,Form("Track pT = %f, DCAtoVertexXY = %f",esdTrack->Pt(),MaxDCAToVertexXYPtDep));
    
    SetMaxDCAToVertexXY(maxDCAToVertexXYPtDep); //Set pT dependent DCA cut to vertex in x-y plane
  } 
  else if( fTrackCutsType == kGlobalCut2011 )
  {
    Float_t maxDCAToVertexXYPtDep = 0.0105 + 0.0350/TMath::Power(esdTrack->Pt(),1.1); 
    // This expression comes from AliESDtrackCuts::GetStandardITSTPCTrackCuts2011()
    
     //AliDebug(3,Form("Track pT = %f, DCAtoVertexXY = %f",esdTrack->Pt(),MaxDCAToVertexXYPtDep));
    
    SetMaxDCAToVertexXY(maxDCAToVertexXYPtDep); //Set pT dependent DCA cut to vertex in x-y plane
  }

  Float_t b[2];
  Float_t bCov[3];
  esdTrack->GetImpactParameters(b,bCov);
  if (bCov[0]<=0 || bCov[2]<=0) 
  {
    AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }

  Float_t dcaToVertexXY = b[0];
  Float_t dcaToVertexZ = b[1];
  Float_t dcaToVertex = -1;

  if (fCutDCAToVertex2D)
    dcaToVertex = TMath::Sqrt(dcaToVertexXY*dcaToVertexXY / fCutMaxDCAToVertexXY/fCutMaxDCAToVertexXY + 
                              dcaToVertexZ*dcaToVertexZ   / fCutMaxDCAToVertexZ/fCutMaxDCAToVertexZ     );
  else
    dcaToVertex = TMath::Sqrt(dcaToVertexXY*dcaToVertexXY + dcaToVertexZ*dcaToVertexZ);
    
  // cut the track?
  
  Bool_t cuts[kNCuts];
  for (Int_t i=0; i<kNCuts; i++) cuts[i]=kFALSE;
  
  // track quality cuts
  if (fCutRequireTPCRefit && (status&AliESDtrack::kTPCrefit)==0)
    cuts[0]=kTRUE;
  if (fCutRequireITSRefit && (status&AliESDtrack::kITSrefit)==0)
    cuts[1]=kTRUE;
  if (nClustersTPC<fCutMinNClusterTPC)
    cuts[2]=kTRUE;
  if (nClustersITS<fCutMinNClusterITS) 
    cuts[3]=kTRUE;
  if (chi2PerClusterTPC>fCutMaxChi2PerClusterTPC) 
    cuts[4]=kTRUE; 
  if (chi2PerClusterITS>fCutMaxChi2PerClusterITS) 
    cuts[5]=kTRUE;  
  if (!fCutAcceptKinkDaughters && esdTrack->GetKinkIndex(0)>0)
    cuts[6]=kTRUE;
  if (fCutDCAToVertex2D && dcaToVertex > 1)
    cuts[7] = kTRUE;
  if (!fCutDCAToVertex2D && TMath::Abs(dcaToVertexXY) > fCutMaxDCAToVertexXY)
    cuts[8] = kTRUE;
  if (!fCutDCAToVertex2D && TMath::Abs(dcaToVertexZ)  > fCutMaxDCAToVertexZ)
    cuts[9] = kTRUE;

  if (fTrackCutsType == kGlobalCut || fTrackCutsType == kGlobalCut2011) 
  {
    //Require at least one SPD point + anything else in ITS
    if ( (esdTrack->HasPointOnITSLayer(0) || esdTrack->HasPointOnITSLayer(1)) == kFALSE)
      cuts[10] = kTRUE;
  }

  // ITS
  if (fCutRequireITSStandAlone || fCutRequireITSpureSA) 
  {
    if ((status & AliESDtrack::kITSin) == 0 || (status & AliESDtrack::kTPCin)) 
    {
      // TPC tracks
      cuts[11] = kTRUE; 
    } 
    else 
    {
      // ITS standalone tracks
      if (fCutRequireITSStandAlone && !fCutRequireITSpureSA) 
      {
        if (status & AliESDtrack::kITSpureSA)    cuts[11] = kTRUE;
      } 
      else if (fCutRequireITSpureSA) 
      {
        if (!(status & AliESDtrack::kITSpureSA)) cuts[11] = kTRUE;
      }
    }
  }
  
  Bool_t cut=kFALSE;
  for (Int_t i=0; i<kNCuts; i++)
    if (cuts[i]) { cut = kTRUE ; }

  // cut the track
  if (cut) 
    return kFALSE;
  else 
    return kTRUE;
}

///
/// Initialize the track cut criteria.
/// By default these cuts are set according to AliESDtrackCuts::GetStandardTPCOnlyTrackCuts().
/// Also, you can customize the cuts using the setters.
///
//_____________________________________
void AliEMCALRecoUtils::InitTrackCuts()
{  
  switch (fTrackCutsType)
  {
    case kTPCOnlyCut:
    {
      AliInfo(Form("Track cuts for matching: AliESDtrackCuts::GetStandardTPCOnlyTrackCuts()"));
      //TPC
      SetMinNClustersTPC(70);
      SetMaxChi2PerClusterTPC(4);
      SetAcceptKinkDaughters(kFALSE);
      SetRequireTPCRefit(kFALSE);
      
      //ITS
      SetRequireITSRefit(kFALSE);
      SetMaxDCAToVertexZ(3.2);
      SetMaxDCAToVertexXY(2.4);
      SetDCAToVertex2D(kTRUE);
      
      break;
    }
      
    case kGlobalCut:
    {
      AliInfo(Form("Track cuts for matching: AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE)"));
      //TPC
      SetMinNClustersTPC(70);
      SetMaxChi2PerClusterTPC(4);
      SetAcceptKinkDaughters(kFALSE);
      SetRequireTPCRefit(kTRUE);
      
      //ITS
      SetRequireITSRefit(kTRUE);
      SetMaxDCAToVertexZ(2);
      SetMaxDCAToVertexXY();
      SetDCAToVertex2D(kFALSE);
      
      break;
    }
      
    case kLooseCut:
    {
      AliInfo(Form("Track cuts for matching: Loose cut w/o DCA cut"));
      SetMinNClustersTPC(50);
      SetAcceptKinkDaughters(kTRUE);
      
      break;
    }
      
    case kITSStandAlone:
    {
      AliInfo(Form("Track cuts for matching: ITS Stand Alone tracks cut w/o DCA cut"));
      SetRequireITSRefit(kTRUE);
      SetRequireITSStandAlone(kTRUE);
      SetITSTrackSA(kTRUE);
      break;
    }
    
    case kGlobalCut2011:
    {
      AliInfo(Form("Track cuts for matching: AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE)"));
      //TPC
      SetMinNClustersTPC(50);
      SetMaxChi2PerClusterTPC(4);
      SetAcceptKinkDaughters(kFALSE);
      SetRequireTPCRefit(kTRUE);
      
      //ITS
      SetRequireITSRefit(kTRUE);
      SetMaxDCAToVertexZ(2);
      SetMaxDCAToVertexXY();
      SetDCAToVertex2D(kFALSE);
      
      break;
    }  
     
    case kLooseCutWithITSrefit:
    {
      AliInfo(Form("Track cuts for matching: Loose cut w/o DCA cut plus ITSrefit"));
      SetMinNClustersTPC(50);
      SetAcceptKinkDaughters(kTRUE);
      SetRequireITSRefit(kTRUE);

      break;
    }
  }
}

///
/// Check if tracks are matched to EMC clusters and set the matched EMCAL cluster index to ESD track. 
///
/// \param event: pointer to event
///
//________________________________________________________________________
void AliEMCALRecoUtils::SetClusterMatchedToTrack(const AliVEvent *event)
{
  Int_t nTracks = event->GetNumberOfTracks();
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) 
  {
    AliVTrack* track = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
    if (!track) 
    {
      AliWarning(Form("Could not receive track %d", iTrack));
      continue;
    }
    
    Int_t matchClusIndex = GetMatchedClusterIndex(iTrack);       
    track->SetEMCALcluster(matchClusIndex); //sets -1 if track not matched within residual
    
    /*the following can be done better if AliVTrack::SetStatus will be there. Patch pending with Andreas/Peter*/
    AliESDtrack* esdtrack = dynamic_cast<AliESDtrack*>(track);
    if (esdtrack) 
    { 
      if (matchClusIndex != -1) 
        esdtrack->SetStatus(AliESDtrack::kEMCALmatch);
      else
        esdtrack->ResetStatus(AliESDtrack::kEMCALmatch);
    } 
    else 
    {
      AliAODTrack* aodtrack = dynamic_cast<AliAODTrack*>(track);
      if (matchClusIndex != -1) 
        aodtrack->SetStatus(AliESDtrack::kEMCALmatch);
      else
        aodtrack->ResetStatus(AliESDtrack::kEMCALmatch);
    }
  }
  AliDebug(2,"Track matched to closest cluster");  
}

///
/// Checks if EMC clusters are matched to ESD track.
/// Adds track indexes of all the tracks matched to a cluster within residuals in AliVCluster.
///
/// \param event: pointer to event
///
//_________________________________________________________________________
void AliEMCALRecoUtils::SetTracksMatchedToCluster(const AliVEvent *event)
{
  for (Int_t iClus=0; iClus < event->GetNumberOfCaloClusters(); ++iClus) 
  {
    AliVCluster *cluster = event->GetCaloCluster(iClus);
    if (!cluster->IsEMCAL()) 
      continue;
    
    //
    // Remove old matches in cluster
    //
    if(cluster->GetNTracksMatched() > 0)
    {
      if(!strcmp("AliESDCaloCluster",Form("%s",cluster->ClassName())))
      {
        TArrayI arrayTrackMatched(0);
        ((AliESDCaloCluster*)cluster)->AddTracksMatched(arrayTrackMatched);
      }
      else
      {
        for(Int_t iTrack = 0; iTrack < cluster->GetNTracksMatched(); iTrack++)
        {
          ((AliAODCaloCluster*)cluster)->RemoveTrackMatched((TObject*)((AliAODCaloCluster*)cluster)->GetTrackMatched(iTrack));
        }
      }
    }
    
    //
    // Find new matches and put them in the cluster
    //
    Int_t nTracks = event->GetNumberOfTracks();
    TArrayI arrayTrackMatched(nTracks);
    
    // Get the closest track matched to the cluster
    Int_t nMatched = 0;
    Int_t matchTrackIndex = GetMatchedTrackIndex(iClus);
    if (matchTrackIndex != -1) 
    {
      arrayTrackMatched[nMatched] = matchTrackIndex;
      nMatched++;
    }
    
    // Get all other tracks matched to the cluster
    for (Int_t iTrk=0; iTrk<nTracks; ++iTrk) 
    {
      AliVTrack* track = dynamic_cast<AliVTrack*>(event->GetTrack(iTrk));
      
      if( !track ) continue;
      
      if ( iTrk == matchTrackIndex ) continue;
      
      if ( track->GetEMCALcluster() == iClus )
      {
        arrayTrackMatched[nMatched] = iTrk;
        ++nMatched;
      }
    }
    
    AliDebug(2,Form("cluster E %f, N matches %d, first match %d\n",cluster->E(),nMatched,arrayTrackMatched[0]));
    
    arrayTrackMatched.Set(nMatched);
    AliESDCaloCluster *esdcluster = dynamic_cast<AliESDCaloCluster*>(cluster);
    if (esdcluster) 
      esdcluster->AddTracksMatched(arrayTrackMatched);
    else if ( nMatched > 0 ) 
    {
      AliAODCaloCluster *aodcluster = dynamic_cast<AliAODCaloCluster*>(cluster);
      if (aodcluster)
      {
        aodcluster->AddTrackMatched(event->GetTrack(arrayTrackMatched.At(0)));
        //AliAODTrack *aodtrack=dynamic_cast<AliAODTrack*>(event->GetTrack(arrayTrackMatched.At(0)));
        //printf("Is the closest matching track with ID %d a 128? %d what's its full filter map? %u\n",aodtrack->GetID(), aodtrack->TestFilterBit(128),aodtrack->GetFilterMap());
        //printf("With specs: pt %.4f, eta %.4f, phi %.4f\n",aodtrack->Pt(),aodtrack->Eta(), aodtrack->Phi());
      }
    }
    
    Float_t eta= -999, phi = -999;
    if (matchTrackIndex != -1) 
      GetMatchedResiduals(iClus, eta, phi);
    
    cluster->SetTrackDistance(phi, eta);
  }
  
  AliDebug(2,"Cluster matched to tracks");  
}

void AliEMCALRecoUtils::SetEMCALChannelRecalibrationFactors(const TObjArray *map) { 
  if(fEMCALRecalibrationFactors) fEMCALRecalibrationFactors->Clear();
  else {
    fEMCALRecalibrationFactors = new TObjArray(map->GetEntries());
    fEMCALRecalibrationFactors->SetOwner(true);
  }
  if(!fEMCALRecalibrationFactors->IsOwner()){
    // Must claim ownership since the new objects are owend by this instance
    fEMCALRecalibrationFactors->SetOwner(kTRUE);
  }

  if(!fUse1Drecalib){
    for(int i = 0; i < map->GetEntries(); i++){
      TH2F *hist = dynamic_cast<TH2F *>(map->At(i));
      if(!hist) continue;
      this->SetEMCALChannelRecalibrationFactors(i, hist);
    }
  }else{
    TH1S *hist = dynamic_cast<TH1S *>(map->At(0));
    this->SetEMCALChannelRecalibrationFactors1D(hist);
  }

}

void AliEMCALRecoUtils::SetEMCALChannelRecalibrationFactors(Int_t iSM , const TH2F* h) { 
  if(!fEMCALRecalibrationFactors){
    fEMCALRecalibrationFactors = new TObjArray(iSM);
    fEMCALRecalibrationFactors->SetOwner(true);
  }
  if(fEMCALRecalibrationFactors->GetEntries() <= iSM) fEMCALRecalibrationFactors->Expand(iSM+1);
  if(fEMCALRecalibrationFactors->At(iSM)) fEMCALRecalibrationFactors->RemoveAt(iSM);
  TH2F *clone = new TH2F(*h);
  clone->SetDirectory(NULL);
  fEMCALRecalibrationFactors->AddAt(clone,iSM); 
}

void AliEMCALRecoUtils::SetEMCALChannelRecalibrationFactors1D(const TH1S* h) { 
  if(!fEMCALRecalibrationFactors){
    fEMCALRecalibrationFactors = new TObjArray(1);
    fEMCALRecalibrationFactors->SetOwner(true);
  }
  if(fEMCALRecalibrationFactors->At(0)) fEMCALRecalibrationFactors->RemoveAt(0);
  TH1S *clone = new TH1S(*h);
  clone->SetDirectory(NULL);
  fEMCALRecalibrationFactors->AddAt(clone,0); 
}

TH2F * AliEMCALRecoUtils::GetEMCALChannelRecalibrationFactors(Int_t iSM) const{

  if(!fUse1Drecalib){
    return (TH2F*)fEMCALRecalibrationFactors->At(iSM);
  }else{
    const  Double_t lowerLimit[]={0,1152,2304,3456,4608,5760,6912,8064,9216,10368,11520,11904,12288,13056,13824,14592,15360,16128,16896,17280};
    const  Double_t upperLimit[]={1151 ,2303 ,3455 ,4607 ,5759 ,6911 ,8063 ,9215 ,10367,11519,11903,12287,13055,13823,14591,15359,16127,16895,17279,17663};

    TH2F* hist = new TH2F(Form("EMCALRecalFactors_SM%d",iSM),Form("EMCALRecalFactors_SM%d",iSM), 48, 0, 48, 24, 0, 24);
    AliEMCALGeometry* geom = AliEMCALGeometry::GetInstance();
    Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 

    for(Int_t iCell=lowerLimit[iSM]; iCell<upperLimit[iSM]; iCell++){
      geom->GetCellIndex(iCell,imod,iTower,iIphi,iIeta); 
      geom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);
      hist->SetBinContent(ieta, iphi, ((TH1S*)fEMCALRecalibrationFactors->At(0))->GetBinContent(iCell)/1000);
    }
    return hist;
  }
}

/*
 Setting EMCAL and DCAL single channel calibration factors using a map
 */
void AliEMCALRecoUtils::SetEMCALSingleChannelRecalibrationFactors(const TObjArray *map) { 
  if(fEMCALSingleChannelRecalibrationFactors) fEMCALSingleChannelRecalibrationFactors->Clear();
  else {
    fEMCALSingleChannelRecalibrationFactors = new TObjArray(map->GetEntries());
    fEMCALSingleChannelRecalibrationFactors->SetOwner(true);
  }
  if(!fEMCALSingleChannelRecalibrationFactors->IsOwner()){
    // Must claim ownership since the new objects are owend by this instance
    fEMCALSingleChannelRecalibrationFactors->SetOwner(kTRUE);
  }
  for(int i = 0; i < map->GetEntries(); i++){
    TH2F *hist = dynamic_cast<TH2F *>(map->At(i));
    if(!hist) continue;
    this->SetEMCALSingleChannelRecalibrationFactors(i, hist);
  }
}

/*
 Setting EMCAL and DCAL single channel calibration factors using an SM by SM histogram
 */
void AliEMCALRecoUtils::SetEMCALSingleChannelRecalibrationFactors(Int_t iSM , const TH2F* h) { 
  if(!fEMCALSingleChannelRecalibrationFactors){
    fEMCALSingleChannelRecalibrationFactors = new TObjArray(iSM);
    fEMCALSingleChannelRecalibrationFactors->SetOwner(true);
  }
  if(fEMCALSingleChannelRecalibrationFactors->GetEntries() <= iSM) fEMCALSingleChannelRecalibrationFactors->Expand(iSM+1);
  if(fEMCALSingleChannelRecalibrationFactors->At(iSM)) fEMCALSingleChannelRecalibrationFactors->RemoveAt(iSM);
  TH2F *clone = new TH2F(*h);
  clone->SetDirectory(NULL);
  fEMCALSingleChannelRecalibrationFactors->AddAt(clone,iSM); 
}

void AliEMCALRecoUtils::SetEMCALChannelStatusMap(const TObjArray *map) { 
  if(fEMCALBadChannelMap) fEMCALBadChannelMap->Clear();
  else {
    fEMCALBadChannelMap = new TObjArray(map->GetEntries());
    fEMCALBadChannelMap->SetOwner(true);
  }
  if(!fEMCALBadChannelMap->IsOwner()) {
    // Must claim ownership since the new objects are owend by this instance
    fEMCALBadChannelMap->SetOwner(true);
  }

  if(!fUse1Dmap){
    for(int i = 0; i < map->GetEntries(); i++){
      TH2I *hist = dynamic_cast<TH2I *>(map->At(i));
      if(!hist) continue;
      this->SetEMCALChannelStatusMap(i, hist);
    }
  }else{
    TH1C *hist = dynamic_cast<TH1C *>(map->At(0));
    this->SetEMCALChannelStatusMap1D(hist);
  }

}

void AliEMCALRecoUtils::SetEMCALChannelStatusMap(Int_t iSM , const TH2I* h) {
  if(!fEMCALBadChannelMap){
    fEMCALBadChannelMap = new TObjArray(iSM);
    fEMCALBadChannelMap->SetOwner(true);
  }
  if(fEMCALBadChannelMap->GetEntries() <= iSM) fEMCALBadChannelMap->Expand(iSM+1);
  if(fEMCALBadChannelMap->At(iSM)) fEMCALBadChannelMap->RemoveAt(iSM);
  TH2I *clone = new TH2I(*h);
  clone->SetDirectory(NULL);
  fEMCALBadChannelMap->AddAt(clone,iSM); 
}

void AliEMCALRecoUtils::SetEMCALChannelStatusMap1D(const TH1C* h) {
  fUse1Dmap = kTRUE;
  if(!fEMCALBadChannelMap){
    fEMCALBadChannelMap = new TObjArray(1);
    fEMCALBadChannelMap->SetOwner(true);
  }
  if(fEMCALBadChannelMap->At(0)) fEMCALBadChannelMap->RemoveAt(0);
  TH1C *clone = new TH1C(*h);
  clone->SetDirectory(NULL);
  fEMCALBadChannelMap->AddAt(clone,0); 
}

TH2I* AliEMCALRecoUtils::GetEMCALChannelStatusMap(Int_t iSM) const{

  if(!fUse1Dmap){
    return (TH2I*)fEMCALBadChannelMap->At(iSM) ;
  }else{
    const  Double_t lowerLimit[]={0,1152,2304,3456,4608,5760,6912,8064,9216,10368,11520,11904,12288,13056,13824,14592,15360,16128,16896,17280};
    const  Double_t upperLimit[]={1151 ,2303 ,3455 ,4607 ,5759 ,6911 ,8063 ,9215 ,10367,11519,11903,12287,13055,13823,14591,15359,16127,16895,17279,17663};

    TH2I* hist = new TH2I(Form("EMCALBadChannelMap_Mod%d",iSM),Form("EMCALBadChannelMap_Mod%d",iSM), 48, 0, 48, 24, 0, 24);
    AliEMCALGeometry* geom = AliEMCALGeometry::GetInstance();
    Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 

    for(Int_t iCell=lowerLimit[iSM]; iCell<upperLimit[iSM]; iCell++){
      geom->GetCellIndex(iCell,imod,iTower,iIphi,iIeta); 
      geom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);
      hist->SetBinContent(ieta, iphi, (Int_t)((TH1C*)fEMCALBadChannelMap->At(0))->GetBinContent(iCell));
    }
    return hist;
  }
}

void  AliEMCALRecoUtils::SetEMCALChannelTimeRecalibrationFactors(const TObjArray *map) { 
  if(fEMCALTimeRecalibrationFactors) fEMCALTimeRecalibrationFactors->Clear();
  else {
    fEMCALTimeRecalibrationFactors = new TObjArray(map->GetEntries());
    fEMCALTimeRecalibrationFactors->SetOwner(true);
  }
  if(!fEMCALTimeRecalibrationFactors->IsOwner()) {
    // Must claim ownership since the new objects are owend by this instance
    fEMCALTimeRecalibrationFactors->SetOwner(kTRUE);
  }
  for(int i = 0; i < map->GetEntries(); i++){
    TH1F *hist = dynamic_cast<TH1F *>(map->At(i));
    if(!hist) continue;
    this->SetEMCALChannelTimeRecalibrationFactors(i, hist);
  }
}

void  AliEMCALRecoUtils::SetEMCALChannelTimeRecalibrationFactors(Int_t bc, const TH1* h){ 
  if(!fEMCALTimeRecalibrationFactors){
    fEMCALTimeRecalibrationFactors = new TObjArray(bc);
    fEMCALTimeRecalibrationFactors->SetOwner(true);
  }
  if(fEMCALTimeRecalibrationFactors->GetEntries() <= bc) fEMCALTimeRecalibrationFactors->Expand(bc+1);
  if(fEMCALTimeRecalibrationFactors->At(bc)) fEMCALTimeRecalibrationFactors->RemoveAt(bc);
  if(fDoUseMergedBC){
    TH1S *clone = new TH1S(*(TH1S*)h);
    clone->SetDirectory(NULL);
    fEMCALTimeRecalibrationFactors->AddAt(clone,bc); 
  }else{
    TH1F *clone = new TH1F(*(TH1F*)h);
    clone->SetDirectory(NULL);
    fEMCALTimeRecalibrationFactors->AddAt(clone,bc); 
  }
}

void AliEMCALRecoUtils::SetEMCALL1PhaseInTimeRecalibrationForAllSM(const TObjArray *map) { 
  if(fEMCALL1PhaseInTimeRecalibration) fEMCALL1PhaseInTimeRecalibration->Clear();
  else {
    fEMCALL1PhaseInTimeRecalibration = new TObjArray(map->GetEntries());
    fEMCALL1PhaseInTimeRecalibration->SetOwner(true);
  }
  if(!fEMCALL1PhaseInTimeRecalibration->IsOwner()) {
    // Must claim ownership since the new objects are owend by this instance
    fEMCALL1PhaseInTimeRecalibration->SetOwner(true);
  }
  for(int i = 0; i < map->GetEntries(); i++) {
    TH1C *hist = dynamic_cast<TH1C *>(map->At(i));
    if(!hist) continue;
    SetEMCALL1PhaseInTimeRecalibrationForAllSM(hist,i);    
  }
}

void AliEMCALRecoUtils::SetEMCALL1PhaseInTimeRecalibrationForAllSM(const TH1C* h, Short_t parNumber) { 
  if(!fEMCALL1PhaseInTimeRecalibration){
    fEMCALL1PhaseInTimeRecalibration = new TObjArray(1);
    fEMCALL1PhaseInTimeRecalibration->SetOwner(true);
  }
  if(fEMCALL1PhaseInTimeRecalibration->GetEntries()<parNumber+1) fEMCALL1PhaseInTimeRecalibration->Expand(parNumber+1); 
  TH1C *clone = new TH1C(*h);
  clone->SetDirectory(NULL);
  fEMCALL1PhaseInTimeRecalibration->AddAt(clone,parNumber); 
}

///
/// Print Parameters.
///
//___________________________________________________
void AliEMCALRecoUtils::Print(const Option_t *) const 
{ 
  printf("-------------------------------------------------------------------------------------------------------------------------------------- \n");
  printf("AliEMCALRecoUtils Settings: \n");
  printf("\tMisalignment shifts\n");
  for (Int_t i=0; i<5; i++) printf("\t\t sector %d, traslation (x,y,z)=(%f,%f,%f), rotation (x,y,z)=(%f,%f,%f)\n",i, 
                                  fMisalTransShift[i*3],fMisalTransShift[i*3+1],fMisalTransShift[i*3+2],
                                  fMisalRotShift[i*3],  fMisalRotShift[i*3+1],  fMisalRotShift[i*3+2]   );
  printf("\tNon linearity function %d, parameters:\n", fNonLinearityFunction);
  if (fNonLinearityFunction != 3) // print only if not kNoCorrection
    for (Int_t i=0; i<10; i++) printf("param[%d]=%f\n",i, fNonLinearityParams[i]);
  
  printf("\tPosition Recalculation option %d, Particle Type %d, fW0 %2.2f, Recalibrate Data %d \n",fPosAlgo,fParticleType,fW0, fRecalibration);

  printf("\tMatching criteria: ");
  if (fCutEtaPhiSum) {
    printf("\t\tsqrt(dEta^2+dPhi^2)<%4.3f\n",fCutR);
  } else if (fCutEtaPhiSeparate) {
    printf("\t\tdEta<%4.3f, dPhi<%4.3f\n",fCutEta,fCutPhi);
  } else {
    printf("\t\tError\n");
    printf("\t\tplease specify your cut criteria\n");
    printf("\t\tTo cut on sqrt(dEta^2+dPhi^2), use: SwitchOnCutEtaPhiSum()\n");
    printf("\t\tTo cut on dEta and dPhi separately, use: SwitchOnCutEtaPhiSeparate()\n");
  }

  printf("\tMass hypothesis = %2.3f [GeV/c^2], extrapolation step to surface = %2.2f[cm], step to cluster = %2.2f[cm]\n",fMass,fStepSurface, fStepCluster);
  printf("\tCluster selection window: dR < %2.0f\n",fClusterWindow);

  printf("\tTrack cuts: \n");
  printf("\t\tMinimum track pT: %1.2f\n",fCutMinTrackPt);
  printf("\t\tAOD track selection: tpc only %d, or hybrid %d, or mask: %d\n",fAODTPCOnlyTracks,fAODHybridTracks, fAODFilterMask);
  printf("\t\tTPCRefit = %d, ITSRefit = %d\n",fCutRequireTPCRefit,fCutRequireITSRefit);
  printf("\t\tAcceptKinks = %d\n",fCutAcceptKinkDaughters);
  printf("\t\tMinNCulsterTPC = %d, MinNClusterITS = %d\n",fCutMinNClusterTPC,fCutMinNClusterITS);
  printf("\t\tMaxChi2TPC = %2.2f, MaxChi2ITS = %2.2f\n",fCutMaxChi2PerClusterTPC,fCutMaxChi2PerClusterITS);
  printf("\t\tDCSToVertex2D = %d, MaxDCAToVertexXY = %2.2f, MaxDCAToVertexZ = %2.2f\n",fCutDCAToVertex2D,fCutMaxDCAToVertexXY,fCutMaxDCAToVertexZ);
  printf("-------------------------------------------------------------------------------------------------------------------------------------- \n");
}
