// $Id:
//
// Analysis task for neutral pions (into two gammas), and for direct photons by subtraction method
//
// Author: D. Herzig, based on code by C. Loizides and  B. Sahlmueller

#include "AliAnalysisTaskEMCALDirGamma.h"
#include <TAxis.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TProfile.h>
#include <TRegexp.h>
#include <TString.h>
#include <TVector2.h>
#include <TArray.h>
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliCDBManager.h"
#include "AliCentrality.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
//#include "AliEMCALRecPoint.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDCaloTrigger.h"
#include "AliESDEvent.h"
#include "AliESDUtils.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliEventplane.h"
#include "AliGeomManager.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliMagF.h"
#include "AliMultiplicity.h"
#include "AliStack.h"
#include "AliTrackerBase.h"
#include "AliTriggerAnalysis.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"

using std::cout;
using std::endl;
using std::max;

ClassImp(AliAnalysisTaskEMCALDirGamma)

//________________________________________________________________________
AliAnalysisTaskEMCALDirGamma::AliAnalysisTaskEMCALDirGamma()
: AliAnalysisTaskSE(),
fCentVar("V0M"),
fCentFrom(0),
fCentTo(100),
fVtxZMin(-10),
fVtxZMax(+10),
fUseQualFlag(1),
fClusName(),
fDoNtuple(0),
fDoConvAna(0),
fDoAfterburner(0),
fAsymMax1(0.3),
fAsymMax2(0.7),
fAsymMax3(1),
fNminCells(2),
fMinE(0.500),
//fM02(100),
fM02min(0.01),
fM02max(0.3),
fMinErat(0),
fMinEcc(-1),
fDoTrMtSmpl(1),
fDoManualRecal(0),
fGeoName("EMCAL_COMPLETEV1"),
fMinNClusPerTr(50),
fIsoDist(0.2),
fTrClassNames(""),
fTrCuts(0),
fPrimTrCuts(0),
fPrimTracksName(""),
fDoTrMatGeom(0),
fTrainMode(0),
fMarkCells(),
fMinL0Time(-1),
fMaxL0Time(1024),
fMcMode(-1),
fEmbedMode(0),
fGeom(0),
fReco(0),
fTrigName(),
fDoPSel(kFALSE),
fIsGeoMatsSet(0),
fRotateMixed(0),
fAddedSignal(0),
fNEvs(0),
fOutput(0),
fTrClassNamesArr(0),
fEsdEv(0),
fAodEv(0),
fPIDResponse(0),
fESDtrackCuts(0),
fRecPoints(0),
fDigits(0),
fEsdClusters(0),
fEsdCells(0),
fAodClusters(0),
fAodCells(0),
fPtRanges(0),
fSelTracks(0),
fSelPrimTracks(0),
fNtuple(0),
ftrcuts(0),
fHeader(0),
fPrimVert(0),
fSpdVert(0),
fTpcVert(0),
fClusters(0),
fTriggers(0),
fMcParts(0),
fHCuts(0x0),
fHVertexZ(0x0),
fHVertexZ2(0x0),
fHCent(0x0),
fHCentQual(0x0),
fHTclsBeforeCuts(0x0),
fHTclsAfterCuts(0x0),
fHClustNoEvt(0),
fHClustAccEvt(0),
fHClustEccentricity(0),
fHClustEtaPhi(0x0),
fHClustEtaPhiIsEMCal(0x0),         
fHEMCalModule0(0x0),   
fHClustEtaPhiRaw(0x0),
fHv0TrackEtaPhi(0x0),
fHv0ClustEtaPhi(0x0),
fHv0TrackEtaPhi2(0x0),
fHv0ClustEtaPhi2(0x0),
fHClustEnergyPt(0x0),
fHClustEnergySM(0x0),
fHClustEnergySigma(0x0),
fHClustSigmaSigma(0x0),
fHClustEtaM02(0x0),				
fHClustPhiM02(0x0),
fHv0TrackPtEMCal(0x0),
fHv0TrackPt(0x0),
fHClustETrackP(0x0),
fHClustEnergyRatioPhoton(0x0), 
fHClustEnergyRatioPion(0x0), 	
fHClustNCellEnergyRatio(0x0),
fHClustEnergyNCell(0x0),
fHClustEnergyNCellPion(0x0),  
fHClustEnergyNCellPhoton(0x0),
fHClustNCellM02Photon(0x0),
fHClustNCellM02Pion(0x0),
fHClustEnergyNCellRaw(0x0),
fHClustEnergyNCellPionRaw(0x0),
fHClustEnergyNCellPhotonRaw(0x0),
fHClustNCellM02PhotonRaw(0x0),
fHClustNCellM02PionRaw(0x0),
fHConvEnergyPt(0x0),
fHPrimTrackPt(0x0),
fHPrimTrackEta(0x0),
fHPrimTrackPhi(0x0),
fHMatchDr(0x0),
fHMatchDz(0x0),
fHMatchEp(0x0),
fHPionEtaPhi(0x0),
fHAddPionEtaPt(0x0),
fHAddPionEtaPtWgt(0x0),
fHPyPionEtaPt(0x0),
fHdr(0),
fHPionMggPt(0x0),
fHPionMggAsym(0x0),
fHPionPtAsym(0x0),
fHPionMggDgg(0x0),
fHPionEtaPhiConv(0x0),
fHPionMggPtConv(0x0),
fHPionMggAsymConv(0x0),
fHPionMggDggConv(0x0),
fHPionInvMassesConv(0x0),
fHPionInvMassesConvMix(0x0),
fHPionEtaPhiConvConv(0x0),
fHPionMggPtConvConv(0x0),
fHPionMggAsymConvConv(0x0),
fHPionMggDggConvConv(0x0),
fHPionInvMassesConvConv	(0x0),
fHPionInvMassesConvConvMix(0x0),
fHConversionPoint(0),
fHPionTruthPt(),
fHPionTruthPtIn(),
fHPionTruthPtAcc(),
fHPionTruthPtConvAcc(),
fHEtaTruthPt(),
fHEtaTruthPtIn(),
fHEtaTruthPtAcc(),
fHEtaTruthPtConvAcc(),
fHGamTruthPt(),
fHGamTruthPtIn(),
fHGamTruthPtAcc(),
fHPionTruthPtAdd(),
fHPionTruthPtInAdd(),
fHPionTruthPtAccAdd(),
fHPionTruthPtConvAccAdd(),
fHEtaTruthPtAdd(),
fHEtaTruthPtInAdd(),
fHEtaTruthPtAccAdd(),
fHEtaTruthPtConvAccAdd(),
fHNMothers(0x0),
//fHMixRotation(),
fHClustEnergyM02Gamma(0x0),
fHClustEnergyM02Pi0(0x0),
fHClustEnergyM02Pion(0x0),
fHClustEnergyM02Elektron(0x0),
fHClustEnergyM02PionTM(0x0),
fHClustEnergyM02AllTM(0x0),
fHClustEnergyM02GammaAll(0x0),  
fHClustEnergyM02All(0x0),
fHClustEnergyM02v0(0x0),
fHClustEnergyEPv0(0x0),
fHClustM02Gamma0TM(0x0),
fHClustM02Gamma1TM(0x0),
fHClustM02Gamma2TM(0x0),
fHClustM02Gamma3TM(0x0),
fHClustM02Gamma4TM(0x0),
fHClustM02Gamma5TM(0x0),
fHClustM02Gamma6TM(0x0),
fHClustM02Gamma7TM(0x0),
fHClustM02Pion0TM(0x0),
fHClustM02Pion1TM(0x0),
fHClustM02Pion2TM(0x0),
fHClustM02Pion3TM(0x0),
fHClustM02Pion4TM(0x0),
fHClustM02Pion5TM(0x0),
fHClustM02Pion6TM(0x0),
fHClustM02Pion7TM(0x0),
fHClustEnergyM02GammaRaw(0x0),    
fHClustEnergyM02Pi0Raw(0x0),      
fHClustEnergyM02PionRaw(0x0),      
fHClustEnergyM02ElektronRaw(0x0),
fHClustEnergyM02GammaAllRaw(0x0), 
fHClustM02M20NoGammaRaw(0x0),      
fHClustM02M20GammaAllRaw(0x0),     
fHClustEnergyM02AllRaw(0x0),	
fHClustEnergyM02GammaSmallCut(0x0),
fHClustEnergyM02PionSmallCut(0x0), 
fHClustEnergyM02GammaCell1(0x0),
fHClustEnergyM02PionCell1(0x0),
fHv0electrons(0x0),
fHPtSpecAll(0x0),
fHPtSpecGamma(0x0),
fHPtSpecPion(0x0),
fHPtSpecElectron(0x0),
fHPtSpecMyon(0x0),
fHPtSpecProton(0x0),
fHPtSpecNeutron(0x0),
fHPtSpecKaon(0x0),
fHPtSpecKaon0(0x0),
fHPtSpecNoGamma(0x0),
fHPtSpecCharged(0x0),
fHPtSpecElectronTM(0x0),
fHPtSpecPionTM(0x0),
fHPtSpecElectronNoTM(0x0),
fHPtSpecPionNoTM(0x0),
fHPtSpecGammaNoM02(0x0),
fHPtSpecPionNoM02(0x0),
fHPtSpecGammaM02(0x0),
fHPtSpecPionM02(0x0),

fHPtSpecAllRaw(0x0),
fHPtSpecGammaRaw(0x0),
fHPtSpecPionRaw(0x0),
fHPtSpecElectronRaw(0x0),
fHPtSpecMyonRaw(0x0),
fHPtSpecProtonRaw(0x0),
fHPtSpecNeutronRaw(0x0),
fHPtSpecKaonRaw(0x0),
fHPtSpecKaon0Raw(0x0),
fHPtSpecNoGammaRaw(0x0),
fHPtSpecChargedRaw(0x0),
fHPtSpecChargedTM(0x0),
fHPtSpecGammaTM(0x0),
fHPtSpecEffTM1(0x0),
fHPtSpecEffTM2(0x0),
fHPtSpecSysEnergy1(0x0),
fHPtSpecSysEnergy2(0x0),
fHPtSpecSysEnergy3(0x0),
fHPtSpecSysEnergy4(0x0),
fHPtSpecSysEnergy5(0x0),
fHPtSpecSysNcell1(0x0),
fHPtSpecSysNcell2(0x0),
fHPtSpecSysNcell3(0x0),
fHPtSpecSysBorder1(0x0),
fHPtSpecSysBorder2(0x0),
fHPtSpecSysBorder3(0x0),
fHPtSpecGammaCompare(0x0),
fHPtSpecPionCompare(0x0),
fHPtSpecCompare(0x0),
fHPtSpecEffParticle(0x0),
fHPtSpecEffPhoton(0x0),
fHPtSpecAccPhoton(0x0),
fHPtSpecEtaPhoton(0x0),
fHPtSpecPhiPhoton(0x0),
fHGenEtaPhi(0x0),
fHPtSpecEffCluster(0x0),
fHPtSpecEffNeutron(0x0),
fHPtSpecConversion(0x0),
fHPtSpecConversionNot(0x0),
fHPtSpecEffPhotonEta5(0x0),
fHPtSpecEffPhotonEta4(0x0),
fHPtSpecEffPhotonEta3(0x0),
fHPtSpecEffPhotonEta2(0x0),
fHPtSpecEffPhotonEta1(0x0),
fHPtSpecDecayPi0(0x0),
fHPtSpecDecayEta(0x0),
fHPtSpecDecayOmega(0x0),
fHPtSpecDecayEtap(0x0),
fHCutVariationM02Photon(0x0),
fHCutVariationM02Pion(0x0),
fHCutVariationM02PhotonTest(0x0),
fHCutVariationM02PionTest(0x0),
fHM02Photon(0x0),
fHM02Pion(0x0),
fHCutVariationPion(0x0),
fHCutVariationPionM021(0x0),
fHCutVariationPionM022(0x0),
fHCutVariationPionM023(0x0),
fHCutVariationPionM024(0x0),
fHCutVariationPionM025(0x0),
fHCutVariationPionM026(0x0),
fHCutVariationPionM027(0x0),
fHCutVariationPionEnergy1(0x0),
fHCutVariationPionEnergy2(0x0),
fHCutVariationPionEnergy3(0x0),
fHCutVariationPionEnergy4(0x0),
fHCutVariationPionEnergy5(0x0),
fHCutVariationPionEnergy6(0x0),
fHCutVariationPionEnergy7(0x0),
fHCutVariationPionNcell1(0x0),
fHCutVariationPionNcell2(0x0),
fHCutVariationPionNcell3(0x0),
fHCutVariationPionNcell4(0x0),
fHCutVariationPionNcell5(0x0),
fHCutVariationPionNcell6(0x0),
fHCutVariationPionNcell7(0x0),
fHCutVariationPhoton(0x0),
fHCutVariationPhotonM021(0x0),
fHCutVariationPhotonM022(0x0),
fHCutVariationPhotonM023(0x0),
fHCutVariationPhotonM024(0x0),
fHCutVariationPhotonM025(0x0),
fHCutVariationPhotonM026(0x0),
fHCutVariationPhotonM027(0x0),
fHCutVariationPhotonEnergy1(0x0),
fHCutVariationPhotonEnergy2(0x0),
fHCutVariationPhotonEnergy3(0x0),
fHCutVariationPhotonEnergy4(0x0),
fHCutVariationPhotonEnergy5(0x0),
fHCutVariationPhotonEnergy6(0x0),
fHCutVariationPhotonEnergy7(0x0),
fHCutVariationPhotonNcell1(0x0),
fHCutVariationPhotonNcell2(0x0),
fHCutVariationPhotonNcell3(0x0),
fHCutVariationPhotonNcell4(0x0),
fHCutVariationPhotonNcell5(0x0),
fHCutVariationPhotonNcell6(0x0),
fHCutVariationPhotonNcell7(0x0),
fHMixRotation(0x0),
fHCorrection(0x0),
fHPionSm(0x0),
fHParticles(0x0),
fHParticlesTM(0x0),
fHParticlesNoTM(0x0),
fHParticlesRaw(0x0),
fHParticlesEff(0x0),
fHParticlesCompare(0x0),
fHParticlesCompare2(0x0),
fHParticlesCompare3(0x0),
fHParticleR(0x0),
fHParticleRcut(0x0),
fHParticleRcutcut(0x0),
fHParticleRcutcutcut(0x0),
fHMotherR(0x0),
fHMotherR2(0x0),
fHClustNcell(0x0),
fHClustNcell1(0x0),
fHClustNcellPhoton(0x0),
fHClustNcellNoPhoton(0x0),
fHClustNcellPhotonCut(0x0),
fHClustNcellNoPhotonCut(0x0),
fHGammaMIP(0x0),
fHHadronMIP(0x0),
fHClustEP(0x0),
fHPtSpecElectronMerge(0x0),
fHClustElectronZR(0x0),
fHclusterTOFdistance(0x0),
fHistTOF(0x0),




ipymin(0),
ipymax(0),
ipi0min(0),
ipi0max(0),
ietamin(0),
ietamax(0),
  fReaderGammas(0),
  eventHeader(0),
  pythiaHeader(0),
  addedPi0Header(0),
  addedEtaHeader(0),
  fHMCpartfrac(0),
  fHECluEMC(0x0),
  fHECluEMCAddPi0(0x0),
  fHECluEMCAddEta(0x0),
//  fHRecTrue(),
//  fHRecTrueAddPi0(),
//  fHRecTrueAddEta(),
  fHECluEMCnofull(),
  fHECluEMCnofullAdd(),
  fHECluEMCelectron(),
  fHECluEMCpion(),
  fHECluEMCkaon(),
  fHECluEMCother(),
  fHECluEMCpi0single(),
//  fHCorrection(),
//  fHPionSm(),
  fHWgt(0),
 fV0cuts(0)
	
{
  // Constructor.
}

//________________________________________________________________________
AliAnalysisTaskEMCALDirGamma::AliAnalysisTaskEMCALDirGamma(const char *name)
: AliAnalysisTaskSE(name),
fCentVar("V0M"),
fCentFrom(0),
fCentTo(100),
fVtxZMin(-10),
fVtxZMax(+10),
fUseQualFlag(1),
fClusName(),
fDoNtuple(0),
fDoConvAna(0),
fDoAfterburner(0),
fAsymMax1(0.3),
fAsymMax2(0.7),
fAsymMax3(1),
fNminCells(2),
fMinE(0.500),
//fM02(100),
fM02min(0.01),
fM02max(0.3),
fMinErat(0),
fMinEcc(-1),
fDoTrMtSmpl(1),
fDoManualRecal(0),
fGeoName("EMCAL_COMPLETEV1"),
fMinNClusPerTr(50),
fIsoDist(0.2),
fTrClassNames(""),
fTrCuts(0),
fPrimTrCuts(0),
fPrimTracksName(""),
fDoTrMatGeom(0),
fTrainMode(0),
fMarkCells(),
fMinL0Time(-1),
fMaxL0Time(1024),
fMcMode(-1),
fEmbedMode(0),
fGeom(0),
fReco(0),
fTrigName(),
fDoPSel(kFALSE),
fIsGeoMatsSet(0),
fRotateMixed(0),
fAddedSignal(0),
fNEvs(0),
fOutput(0),
fTrClassNamesArr(0),
fEsdEv(0),
fAodEv(0),
fPIDResponse(0),
fESDtrackCuts(0),
fRecPoints(0),
fDigits(0),
fEsdClusters(0),
fEsdCells(0),
fAodClusters(0),
fAodCells(0),
fPtRanges(0),
fSelTracks(0),
fSelPrimTracks(0),
fNtuple(0),
ftrcuts(0),
fHeader(0),
fPrimVert(0),
fSpdVert(0),
fTpcVert(0),
fClusters(0),
fTriggers(0),
fMcParts(0),
fHCuts(0x0),
fHVertexZ(0x0),
fHVertexZ2(0x0),
fHCent(0x0),
fHCentQual(0x0),
fHTclsBeforeCuts(0x0),
fHTclsAfterCuts(0x0),
fHClustNoEvt(0),
fHClustAccEvt(0),
fHClustEccentricity(0),
fHClustEtaPhi(0x0),
fHClustEtaPhiIsEMCal(0x0),        
fHEMCalModule0(0x0),   
fHClustEtaPhiRaw(0x0),
fHv0TrackEtaPhi(0x0),
fHv0ClustEtaPhi(0x0),
fHv0TrackEtaPhi2(0x0),
fHv0ClustEtaPhi2(0x0),
fHClustEnergyPt(0x0),
fHClustEnergySM(0x0),
fHClustEnergySigma(0x0),
fHClustSigmaSigma(0x0),
fHClustEtaM02(0x0),				
fHClustPhiM02(0x0),
fHv0TrackPtEMCal(0x0),
fHv0TrackPt(0x0),
fHClustETrackP(0x0),
fHClustEnergyRatioPhoton(0x0), 
fHClustEnergyRatioPion(0x0), 	
fHClustNCellEnergyRatio(0x0),
fHClustEnergyNCell(0x0),
fHClustEnergyNCellPion(0x0),  
fHClustEnergyNCellPhoton(0x0),
fHClustNCellM02Photon(0x0),
fHClustNCellM02Pion(0x0),
fHClustEnergyNCellRaw(0x0),
fHClustEnergyNCellPionRaw(0x0),
fHClustEnergyNCellPhotonRaw(0x0),
fHClustNCellM02PhotonRaw(0x0),
fHClustNCellM02PionRaw(0x0),
fHConvEnergyPt(0x0),
fHPrimTrackPt(0x0),
fHPrimTrackEta(0x0),
fHPrimTrackPhi(0x0),
fHMatchDr(0x0),
fHMatchDz(0x0),
fHMatchEp(0x0),
fHPionEtaPhi(0x0),
fHAddPionEtaPt(0x0),
fHAddPionEtaPtWgt(0x0),
fHPyPionEtaPt(0x0),
fHdr(0),
fHPionMggPt(0x0),
fHPionMggAsym(0x0),
fHPionPtAsym(0x0),
fHPionMggDgg(0x0),
fHPionEtaPhiConv(0x0),
fHPionMggPtConv(0x0),
fHPionMggAsymConv(0x0),
fHPionMggDggConv(0x0),
fHPionInvMassesConv(0x0),
fHPionInvMassesConvMix(0x0),
fHPionEtaPhiConvConv(0x0),
fHPionMggPtConvConv(0x0),
fHPionMggAsymConvConv(0x0),
fHPionMggDggConvConv(0x0),
fHPionInvMassesConvConv	(0x0),
fHPionInvMassesConvConvMix(0x0),
fHConversionPoint(0),
fHPionTruthPt(),
fHPionTruthPtIn(),
fHPionTruthPtAcc(),
fHPionTruthPtConvAcc(),
fHEtaTruthPt(),
fHEtaTruthPtIn(),
fHEtaTruthPtAcc(),
fHEtaTruthPtConvAcc(),
fHGamTruthPt(),
fHGamTruthPtIn(),
fHGamTruthPtAcc(),
fHPionTruthPtAdd(),
fHPionTruthPtInAdd(),
fHPionTruthPtAccAdd(),
fHPionTruthPtConvAccAdd(),
fHEtaTruthPtAdd(),
fHEtaTruthPtInAdd(),
fHEtaTruthPtAccAdd(),
fHEtaTruthPtConvAccAdd(),
fHNMothers(0x0),
//fHMixRotation(),

fHClustEnergyM02Gamma(0x0),
fHClustEnergyM02Pi0(0x0),
fHClustEnergyM02Pion(0x0),
fHClustEnergyM02Elektron(0x0),
fHClustEnergyM02PionTM(0x0),
fHClustEnergyM02AllTM(0x0),
fHClustEnergyM02GammaAll(0x0),  
fHClustEnergyM02All(0x0),
fHClustEnergyM02v0(0x0),
fHClustEnergyEPv0(0x0),
fHClustM02Gamma0TM(0x0),
fHClustM02Gamma1TM(0x0),
fHClustM02Gamma2TM(0x0),
fHClustM02Gamma3TM(0x0),
fHClustM02Gamma4TM(0x0),
fHClustM02Gamma5TM(0x0),
fHClustM02Gamma6TM(0x0),
fHClustM02Gamma7TM(0x0),
fHClustM02Pion0TM(0x0),
fHClustM02Pion1TM(0x0),
fHClustM02Pion2TM(0x0),
fHClustM02Pion3TM(0x0),
fHClustM02Pion4TM(0x0),
fHClustM02Pion5TM(0x0),
fHClustM02Pion6TM(0x0),
fHClustM02Pion7TM(0x0),
fHClustEnergyM02GammaRaw(0x0),    
fHClustEnergyM02Pi0Raw(0x0),      
fHClustEnergyM02PionRaw(0x0),      
fHClustEnergyM02ElektronRaw(0x0),
fHClustEnergyM02GammaAllRaw(0x0), 
fHClustM02M20NoGammaRaw(0x0),      
fHClustM02M20GammaAllRaw(0x0),     
fHClustEnergyM02AllRaw(0x0),	
fHClustEnergyM02GammaSmallCut(0x0),
fHClustEnergyM02PionSmallCut(0x0), 
fHClustEnergyM02GammaCell1(0x0),
fHClustEnergyM02PionCell1(0x0),
fHv0electrons(0x0),
fHPtSpecAll(0x0),
fHPtSpecGamma(0x0),
fHPtSpecPion(0x0),
fHPtSpecElectron(0x0),
fHPtSpecMyon(0x0),
fHPtSpecProton(0x0),
fHPtSpecNeutron(0x0),
fHPtSpecKaon(0x0),
fHPtSpecKaon0(0x0),
fHPtSpecNoGamma(0x0),
fHPtSpecCharged(0x0),
fHPtSpecElectronTM(0x0),
fHPtSpecPionTM(0x0),
fHPtSpecElectronNoTM(0x0),
fHPtSpecPionNoTM(0x0),
fHPtSpecGammaNoM02(0x0),
fHPtSpecPionNoM02(0x0),
fHPtSpecGammaM02(0x0),
fHPtSpecPionM02(0x0),
fHPtSpecAllRaw(0x0),
fHPtSpecGammaRaw(0x0),
fHPtSpecPionRaw(0x0),
fHPtSpecElectronRaw(0x0),
fHPtSpecMyonRaw(0x0),
fHPtSpecProtonRaw(0x0),
fHPtSpecNeutronRaw(0x0),
fHPtSpecKaonRaw(0x0),
fHPtSpecKaon0Raw(0x0),
fHPtSpecNoGammaRaw(0x0),
fHPtSpecChargedRaw(0x0),
fHPtSpecChargedTM(0x0),
fHPtSpecGammaTM(0x0),
fHPtSpecEffTM1(0x0),
fHPtSpecEffTM2(0x0),
fHPtSpecSysEnergy1(0x0),
fHPtSpecSysEnergy2(0x0),
fHPtSpecSysEnergy3(0x0),
fHPtSpecSysEnergy4(0x0),
fHPtSpecSysEnergy5(0x0),
fHPtSpecSysNcell1(0x0),
fHPtSpecSysNcell2(0x0),
fHPtSpecSysNcell3(0x0),
fHPtSpecSysBorder1(0x0),
fHPtSpecSysBorder2(0x0),
fHPtSpecSysBorder3(0x0),
fHPtSpecGammaCompare(0x0),
fHPtSpecPionCompare(0x0),
fHPtSpecCompare(0x0),
fHPtSpecEffParticle(0x0),
fHPtSpecEffPhoton(0x0),
fHPtSpecAccPhoton(0x0),
fHPtSpecEtaPhoton(0x0),
fHPtSpecPhiPhoton(0x0),
fHGenEtaPhi(0x0),
fHPtSpecEffCluster(0x0),
fHPtSpecEffNeutron(0x0),
fHPtSpecConversion(0x0),
fHPtSpecConversionNot(0x0),
fHPtSpecEffPhotonEta5(0x0),
fHPtSpecEffPhotonEta4(0x0),
fHPtSpecEffPhotonEta3(0x0),
fHPtSpecEffPhotonEta2(0x0),
fHPtSpecEffPhotonEta1(0x0),
fHPtSpecDecayPi0(0x0),
fHPtSpecDecayEta(0x0),
fHPtSpecDecayOmega(0x0),
fHPtSpecDecayEtap(0x0),
fHCutVariationM02Photon(0x0),
fHCutVariationM02Pion(0x0),
fHCutVariationM02PhotonTest(0x0),
fHCutVariationM02PionTest(0x0),
fHM02Photon(0x0),
fHM02Pion(0x0),
fHCutVariationPion(0x0),
fHCutVariationPionM021(0x0),
fHCutVariationPionM022(0x0),
fHCutVariationPionM023(0x0),
fHCutVariationPionM024(0x0),
fHCutVariationPionM025(0x0),
fHCutVariationPionM026(0x0),
fHCutVariationPionM027(0x0),
fHCutVariationPionEnergy1(0x0),
fHCutVariationPionEnergy2(0x0),
fHCutVariationPionEnergy3(0x0),
fHCutVariationPionEnergy4(0x0),
fHCutVariationPionEnergy5(0x0),
fHCutVariationPionEnergy6(0x0),
fHCutVariationPionEnergy7(0x0),
fHCutVariationPionNcell1(0x0),
fHCutVariationPionNcell2(0x0),
fHCutVariationPionNcell3(0x0),
fHCutVariationPionNcell4(0x0),
fHCutVariationPionNcell5(0x0),
fHCutVariationPionNcell6(0x0),
fHCutVariationPionNcell7(0x0),
fHCutVariationPhoton(0x0),
fHCutVariationPhotonM021(0x0),
fHCutVariationPhotonM022(0x0),
fHCutVariationPhotonM023(0x0),
fHCutVariationPhotonM024(0x0),
fHCutVariationPhotonM025(0x0),
fHCutVariationPhotonM026(0x0),
fHCutVariationPhotonM027(0x0),
fHCutVariationPhotonEnergy1(0x0),
fHCutVariationPhotonEnergy2(0x0),
fHCutVariationPhotonEnergy3(0x0),
fHCutVariationPhotonEnergy4(0x0),
fHCutVariationPhotonEnergy5(0x0),
fHCutVariationPhotonEnergy6(0x0),
fHCutVariationPhotonEnergy7(0x0),
fHCutVariationPhotonNcell1(0x0),
fHCutVariationPhotonNcell2(0x0),
fHCutVariationPhotonNcell3(0x0),
fHCutVariationPhotonNcell4(0x0),
fHCutVariationPhotonNcell5(0x0),
fHCutVariationPhotonNcell6(0x0),
fHCutVariationPhotonNcell7(0x0),
fHMixRotation(0x0),
fHCorrection(0x0),
fHPionSm(0x0),
fHParticles(0x0),
fHParticlesTM(0x0),
fHParticlesNoTM(0x0),
fHParticlesRaw(0x0),
fHParticlesEff(0x0),
fHParticlesCompare(0x0),
fHParticlesCompare2(0x0),
fHParticlesCompare3(0x0),
fHParticleR(0x0),
fHParticleRcut(0x0),
fHParticleRcutcut(0x0),
fHParticleRcutcutcut(0x0),
fHMotherR(0x0),
fHMotherR2(0x0),
fHClustNcell(0x0),
fHClustNcell1(0x0),
fHClustNcellPhoton(0x0),
fHClustNcellNoPhoton(0x0),
fHClustNcellPhotonCut(0x0),
fHClustNcellNoPhotonCut(0x0),
fHGammaMIP(0x0),
fHHadronMIP(0x0),
fHClustEP(0x0),
fHPtSpecElectronMerge(0x0),
fHClustElectronZR(0x0),
fHclusterTOFdistance(0x0),
fHistTOF(0x0),


ipymin(0),
ipymax(0),
ipi0min(0),
ipi0max(0),
ietamin(0),
ietamax(0),
fReaderGammas(0),
eventHeader(0),
pythiaHeader(0),
addedPi0Header(0),
addedEtaHeader(0),
fHMCpartfrac(0),
fHECluEMC(0x0),
fHECluEMCAddPi0(0x0),
fHECluEMCAddEta(0x0),
//fHRecTrue(),
//fHRecTrueAddPi0(),
//fHRecTrueAddEta(),
fHECluEMCnofull(),
fHECluEMCnofullAdd(),
fHECluEMCelectron(),
fHECluEMCpion(),
fHECluEMCkaon(),
fHECluEMCother(),
fHECluEMCpi0single(),
//  fHCorrection(),
//  fHPionSm(),
  fHWgt(0),
 fV0cuts(0)
{
  // Constructor.
  DefineOutput(1, TList::Class());
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.,SPDVertex.,TPCVertex.,EMCALCells.,Tracks,EMCALTrigger.,SPDPileupVertices,TrkPileupVertices "
  "AOD:header,vertices,emcalCells,tracks";
}

//________________________________________________________________________
AliAnalysisTaskEMCALDirGamma::~AliAnalysisTaskEMCALDirGamma()
{
  // Destructor.
  
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput; 
    //fOutput = 0;
  }
  //delete fPtRanges; fPtRanges = 0;
  //fGeom = 0; // do not delete geometry when using instance
  //delete fReco;
  //fReco = 0;
  //delete fTrClassNamesArr;
  //delete fSelTracks;
  //delete fSelPrimTracks;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::UserCreateOutputObjects()
{
  // Create user objects here.
  
  cout << "AliAnalysisTaskEMCALDirGamma: Input settings" << endl;
  cout << " fCentVar:       " << fCentVar << endl;
  cout << " fCentFrom:      " << fCentFrom << endl;
  cout << " fCentTo:        " << fCentTo << endl;
  cout << " fVtxZMin:       " << fVtxZMin << endl;
  cout << " fVtxZMax:       " << fVtxZMax << endl;
  cout << " fUseQualFlag:   " << fUseQualFlag << endl;
  cout << " fClusName:      \"" << fClusName << "\"" << endl;
  cout << " fDoNtuple:      " << fDoNtuple << endl;
  cout << " fDoAfterburner: " << fDoAfterburner << endl;
  cout << " fAsymMax1:       " << fAsymMax1 << endl;
  cout << " fAsymMax2:       " << fAsymMax2 << endl;
  cout << " fAsymMax3:       " << fAsymMax3 << endl;
  cout << " fNminCells:     " << fNminCells << endl;
  cout << " fMinE:          " << fMinE << endl;
  cout << " fMinErat:       " << fMinErat << endl;
  cout << " fMinEcc:        " << fMinEcc << endl;
  cout << " fM02max:           " << fM02max << endl;
  cout << " fGeoName:       \"" << fGeoName << "\"" << endl;
  cout << " fMinNClusPerTr: " << fMinNClusPerTr << endl;
  cout << " fIsoDist:       " << fIsoDist << endl;
  cout << " fTrClassNames:  \"" << fTrClassNames << "\"" << endl;
  cout << " fTrCuts:        " << fTrCuts << endl;
  cout << " fPrimTrCuts:    " << fPrimTrCuts << endl;
  cout << " fDoTrMatGeom:   " << fDoTrMatGeom << endl;
  cout << " fTrainMode:     " << fTrainMode << endl;
  cout << " fMarkCells:     " << fMarkCells << endl;
  cout << " fMinL0Time:     " << fMinL0Time << endl;
  cout << " fMaxL0Time:     " << fMaxL0Time << endl;
  cout << " fMcMode:        " << fMcMode << endl;
  cout << " fEmbedMode:     " << fEmbedMode << endl;
  cout << " fGeom:          " << fGeom << endl;
  cout << " fReco:          " << fReco << endl;
  cout << " fTrigName:      " << fTrigName << endl;
  cout << " fDoPSel:        " << fDoPSel << endl;
	cout << " fRotateMixed:   " << fRotateMixed << endl;
  
  if (!fGeom)
    fGeom = AliEMCALGeometry::GetInstance(fGeoName);
  else {
    if (fGeom->GetMatrixForSuperModule(0))
      fIsGeoMatsSet = kTRUE;
  }
  if (!fReco)
    fReco = new AliEMCALRecoUtils();
  fTrClassNamesArr = fTrClassNames.Tokenize(" ");
  fOutput = new TList();
  fOutput->SetOwner();
  fSelTracks = new TObjArray;
  fSelPrimTracks = new TObjArray;
  if(fMcMode){
    if (TClass::GetClass("AliStaPart"))
      //TClass::GetClass("AliStaPart")->IgnoreTObjectStreamer();
      fMcParts = new TClonesArray("AliStaPart");
  }
  
  if (fDoNtuple) {
    TFile *f = OpenFile(1);
    TDirectory::TContext context(f);
    if (f) {
      f->SetCompressionLevel(2);
      fNtuple = new TTree(Form("tree%.0fto%.0f",fCentFrom,fCentTo), "StandaloneTree");
      fNtuple->SetDirectory(f);
      if (fTrainMode) {
        fNtuple->SetAutoFlush(-2*1024*1024);
        fNtuple->SetAutoSave(0);
      } else {
        fNtuple->SetAutoFlush(-32*1024*1024);
        fNtuple->SetAutoSave(0);
      }
      
      fHeader = new AliStaHeader;
      fNtuple->Branch("header", &fHeader, 16*1024, 99);
      fPrimVert = new AliStaVertex;
      fNtuple->Branch("primv", &fPrimVert, 16*1024, 99);
      fSpdVert = new AliStaVertex;
      fNtuple->Branch("spdv", &fSpdVert, 16*1024, 99);
      fTpcVert = new AliStaVertex;
      fNtuple->Branch("tpcv", &fTpcVert, 16*1024, 99);
      if (TClass::GetClass("AliStaCluster"))
        //TClass::GetClass("AliStaCluster")->IgnoreTObjectStreamer();
        fClusters = new TClonesArray("AliStaCluster");
      fNtuple->Branch("clusters", &fClusters, 8*16*1024, 99);
      if (TClass::GetClass("AliStaTrigger"))
        //TClass::GetClass("AliStaTrigger")->IgnoreTObjectStreamer();
        fTriggers = new TClonesArray("AliStaTrigger");
      fNtuple->Branch("l0prim", &fTriggers, 16*1024, 99);
      if (fMcMode||fEmbedMode) {
        fNtuple->Branch("mcparts", &fMcParts, 8*16*1024, 99);
      }
    }
  }
  
  AliEMCALEMCGeometry *emc = fGeom->GetEMCGeometry();
  Double_t phimin = emc->GetArm1PhiMin()*TMath::DegToRad();
  Double_t phimax = emc->GetArm1PhiMax()*TMath::DegToRad();
  Double_t etamin = emc->GetArm1EtaMin()*TMath::DegToRad();
  Double_t etamax = emc->GetArm1EtaMax()*TMath::DegToRad();
  
  cout << "phimin=" << phimin << ", phimax=" << phimax << endl;
  
  // histograms
  Bool_t th1 =   TH1::GetDefaultSumw2();
  TH1::SetDefaultSumw2(kTRUE);
  Bool_t th2 =   TH2::GetDefaultSumw2();
  TH2::SetDefaultSumw2(kTRUE);
  fHCuts = new TH1F("hEventCuts","",7,0.5,7.5);
  fHCuts->GetXaxis()->SetBinLabel(1,"All");
  fHCuts->GetXaxis()->SetBinLabel(2,"tracks > 1");   
  fHCuts->GetXaxis()->SetBinLabel(3,"AllDetectos");
  fHCuts->GetXaxis()->SetBinLabel(2,"PS");
  fHCuts->GetXaxis()->SetBinLabel(3,Form("%s: %.0f-%.0f",fCentVar.Data(),fCentFrom,fCentTo));
  fHCuts->GetXaxis()->SetBinLabel(4,"QFlag");
  fHCuts->GetXaxis()->SetBinLabel(5,Form("zvtx: %.0f-%.0f",fVtxZMin,fVtxZMax));
  fOutput->Add(fHCuts);
  fHVertexZ = new TH1F("hVertexZBeforeCut","",100,-25,25);
  fHVertexZ->SetXTitle("z (cm)");
  fOutput->Add(fHVertexZ);
  fHVertexZ2 = new TH1F("hVertexZAfterCut","",100,-25,25);
  fHVertexZ2->SetXTitle("z (cm)");
  fOutput->Add(fHVertexZ2);
  fHCent = new TH1F("hCentBeforeCut","",102,-1,101);
  fHCent->SetXTitle(fCentVar.Data());
  fOutput->Add(fHCent);
  fHCentQual = new TH1F("hCentAfterCut","",102,-1,101);
  fHCentQual->SetXTitle(fCentVar.Data());
  fOutput->Add(fHCentQual);
  fHTclsBeforeCuts = new TH1F("hTclsBeforeCuts","",fTrClassNamesArr->GetEntries(),0.5,0.5+fTrClassNamesArr->GetEntries());
  fHTclsAfterCuts = new TH1F("hTclsAfterCuts","",fTrClassNamesArr->GetEntries(),0.5,0.5+fTrClassNamesArr->GetEntries());
  for (Int_t i = 0; i<fTrClassNamesArr->GetEntries(); ++i) {
    const char *name = fTrClassNamesArr->At(i)->GetName();
    fHTclsBeforeCuts->GetXaxis()->SetBinLabel(1+i,name);
    fHTclsAfterCuts->GetXaxis()->SetBinLabel(1+i,name);
  }
  fOutput->Add(fHTclsBeforeCuts);
  fOutput->Add(fHTclsAfterCuts);
  
  // histograms for clusters
  Int_t nsm = fGeom->GetEMCGeometry()->GetNumberOfSuperModules();
  if (!fTrainMode) {
    fHClustNoEvt = new TH1F("hClustNoEvt","",2000,0,2000);
    fHClustNoEvt->SetXTitle("# Clusters");
    fOutput->Add(fHClustNoEvt);
    fHClustAccEvt = new TH1F("hClustAccEvt","",2000,0,2000);
    fHClustAccEvt->SetXTitle("# Clusters");
    fOutput->Add(fHClustAccEvt);
    fHClustEccentricity = new TH1F("hClustEccentricity","",500,-0.1,3);
    fHClustEccentricity->SetXTitle("#epsilon_{C}");
    fOutput->Add(fHClustEccentricity);
    fHClustEtaPhiRaw = new TH2F("hClustEtaPhiRaw","",160,-0.8,0.8,400,0,6);
    fHClustEtaPhiRaw->SetXTitle("#eta");
    fHClustEtaPhiRaw->SetYTitle("#varphi");
    fOutput->Add(fHClustEtaPhiRaw);
    fHClustEtaPhi = new TH2F("hClustEtaPhi","",400,-2,2,400,0,6);
    fHClustEtaPhi->SetXTitle("#eta");
    fHClustEtaPhi->SetYTitle("#varphi");
    fOutput->Add(fHClustEtaPhi);	
    fHClustEtaPhiIsEMCal = new TH2F("hClustEtaPhiIsEMCal","",320,-0.8,0.8,500,1.1,3.4);
    fHClustEtaPhiIsEMCal->SetXTitle("#eta");
    fHClustEtaPhiIsEMCal->SetYTitle("#varphi");
    fOutput->Add(fHClustEtaPhiIsEMCal);		
    fHEMCalModule0 = new TH2F("hEMCalModule0","",400,-2,2,400,0,6);
    fHEMCalModule0->SetXTitle("X");
    fHEMCalModule0->SetYTitle("Y");
    fOutput->Add(fHEMCalModule0);			
    fHClustEnergyPt = new TH2F("hClustEnergyPt","",250,0,50,250,0,50);
    fHClustEnergyPt->SetXTitle("E (GeV)");
    fHClustEnergyPt->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHClustEnergyPt);
    fHClustEnergySM = new TH2F("hClustEnergySM","",250,0,50,17,-0.5,16.5);
    fHClustEnergySM->SetXTitle("E (GeV)");
    fHClustEnergySM->SetYTitle("SM number");
    fOutput->Add(fHClustEnergySM);
    fHClustEnergySigma = new TH2F("hClustEnergySigma","",50,0,5,250,0,50);
    fHClustEnergySigma->SetXTitle("E_{C} * #sigma_{max} [GeV*cm]");
    fHClustEnergySigma->SetYTitle("E_{C} (GeV)");
    fOutput->Add(fHClustEnergySigma);
    fHClustSigmaSigma = new TH2F("hClustSigmaSigma","",100,0,10,50,0,1);
    fHClustSigmaSigma->SetXTitle("#lambda_{0} (cm)");
    fHClustSigmaSigma->SetYTitle("#sigma_{max} (cm)");
    fOutput->Add(fHClustSigmaSigma);
    fHClustNCellEnergyRatio = new TH2F("hClustNCellEnergyRatio","",27,-0.5,26.5,101,-0.05,1.05);
    fHClustNCellEnergyRatio->SetXTitle("N_{cells}");
    fHClustNCellEnergyRatio->SetYTitle("E^{max}_{cell}/E_{clus}");
    fOutput->Add(fHClustNCellEnergyRatio);
    fHClustEnergyNCell = new TH2F("hClustEnergyNCell","",200,0,100,50,0,50);
    fHClustEnergyNCell->SetXTitle("E_{clus}");
    fHClustEnergyNCell->SetYTitle("N_{cells}");
    fOutput->Add(fHClustEnergyNCell);
    fHClustEtaM02 = new TH2F("hClustEtaM02","",160,-0.8,0.8,200,0,3);
    fHClustEtaM02->SetXTitle("#eta");
    fHClustEtaM02->SetYTitle("M_{02}");
    fOutput->Add(fHClustEtaM02);
    fHClustPhiM02 = new TH2F("hClustPhiM02","",400,0,6,200,0,3);
    fHClustPhiM02->SetXTitle("#phi");
    fHClustPhiM02->SetYTitle("M_{02}");
    fOutput->Add(fHClustPhiM02);	
	
    fHPtSpecEtaPhoton = new TH1F("hPtSpecEtaPhoton","",400,-2,2);
    fHPtSpecEtaPhoton->SetXTitle("#eta");
    fHPtSpecEtaPhoton->SetYTitle("#varphi");
    fOutput->Add(fHPtSpecEtaPhoton);
    fHPtSpecPhiPhoton = new TH1F("hPtSpecPhiPhoton","",400,0,6);
    fHPtSpecPhiPhoton->SetXTitle("#eta");
    fHPtSpecPhiPhoton->SetYTitle("#varphi");
    fOutput->Add(fHPtSpecPhiPhoton);	
    fHGenEtaPhi = new TH2F("hGenEtaPhi","",400,-2,2,400,0,6);
    fHGenEtaPhi->SetXTitle("#eta");
    fHGenEtaPhi->SetYTitle("#varphi");
    fOutput->Add(fHGenEtaPhi);	
	
    fHClustEnergyRatioPhoton = new TH1F("hClustEnergyRatioPhoton","",101,-0.05,1.05);
    fHClustEnergyRatioPhoton->SetXTitle("E^{max}_{cell}/E_{clus}");
    fOutput->Add(fHClustEnergyRatioPhoton);
    fHClustEnergyRatioPion = new TH1F("hClustEnergyRatioPion","",101,-0.05,1.05);
    fHClustEnergyRatioPion->SetXTitle("E^{max}_{cell}/E_{clus}");
    fOutput->Add(fHClustEnergyRatioPion);	
  }
  
  if(fDoConvAna){
    // histogram for conversion photons
    fHConvEnergyPt = new TH2F("hConvEnergyPt","",250,0,50,250,0,50);
    fHConvEnergyPt->SetXTitle("E (GeV)");
    fHConvEnergyPt->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHConvEnergyPt);
  }
  
  // histograms for primary tracks
  fHPrimTrackPt = new TH1F("hPrimTrackPt",";p_{T} [GeV/c]",500,0,50);
  fOutput->Add(fHPrimTrackPt);
  fHPrimTrackEta = new TH1F("hPrimTrackEta",";#eta",40,-2,2);
  fOutput->Add(fHPrimTrackEta);
  fHPrimTrackPhi = new TH1F("hPrimTrackPhi",";#varPhi [rad]",63,0,6.3);
  fOutput->Add(fHPrimTrackPhi);
  
  // histograms for track matching
  if (fDoTrMatGeom) {
    fHMatchDr = new TH1F("hMatchDrDist",";dR (cm)",500,0,200);
    fOutput->Add(fHMatchDr);
    fHMatchDz = new TH1F("hMatchDzDist",";dZ (cm)",500,-100,100);
    fOutput->Add(fHMatchDz);
    fHMatchEp = new TH1F("hMatchEpDist",";E/p",100,0,10);
    fOutput->Add(fHMatchEp);
  }
  
  const Int_t nbins = 120;
  const Int_t ptmax = 30;
  
  const Int_t ebins = 250;
  const Int_t emax = 20;
  
  // d_r of pairs
  fHdr = new TH1F("hdr","",2000,0,2);
  fHdr->SetXTitle("pair d_r");
  fOutput->Add(fHdr);
  
  if (!fTrainMode) {
    
    const Int_t massbins = 170;
    const Double_t massmax = 0.85;

    
    // distribution of particle weights
    fHWgt = new TH1F("hWgt","hWgt",100,0,10);
    fOutput->Add(fHWgt);

		
      	fHClustEnergyM02Gamma = new TH2F("hClustEnergyM02Gamma","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02Gamma->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02Gamma->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02Gamma);
      	fHClustEnergyM02Pion = new TH2F("hClustEnergyM02Pion","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02Pion->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02Pion->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02Pion);	
      	fHClustEnergyM02Elektron = new TH2F("hClustEnergyM02Elektron","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02Elektron->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02Elektron->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02Elektron);			
      	fHClustEnergyM02GammaAll = new TH2F("hClustEnergyM02GammaAll","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02GammaAll->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02GammaAll->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02GammaAll);
      	fHClustEnergyM02All = new TH2F("hClustEnergyM02All","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02All->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02All->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02All);	
			
				
      	fHClustM02Gamma0TM = new TH1F("hClustM02Gamma0TM","",150,0,1);
      	fHClustM02Gamma0TM->SetXTitle("M_{02}");
		fHClustM02Gamma0TM->SetYTitle("Anzahl Teilchen");			
		fHClustM02Gamma0TM->Sumw2();
      	fOutput->Add(fHClustM02Gamma0TM);
      	fHClustM02Gamma1TM = new TH1F("hClustM02Gamma1TM","",150,0,1);
      	fHClustM02Gamma1TM->SetXTitle("M_{02}");
		fHClustM02Gamma1TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Gamma1TM->Sumw2();
      	fOutput->Add(fHClustM02Gamma1TM);
      	fHClustM02Gamma2TM = new TH1F("hClustM02Gamma2TM","",150,0,1);
      	fHClustM02Gamma2TM->SetXTitle("M_{02}");
		fHClustM02Gamma2TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Gamma2TM->Sumw2();
      	fOutput->Add(fHClustM02Gamma2TM);		
      	fHClustM02Gamma3TM = new TH1F("hClustM02Gamma3TM","",150,0,1);
      	fHClustM02Gamma3TM->SetXTitle("M_{02}");
		fHClustM02Gamma3TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Gamma3TM->Sumw2();	
      	fOutput->Add(fHClustM02Gamma3TM);	
      	fHClustM02Gamma4TM = new TH1F("hClustM02Gamma4TM","",150,0,1);
      	fHClustM02Gamma4TM->SetXTitle("M_{02}");
		fHClustM02Gamma4TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Gamma4TM->Sumw2();
      	fOutput->Add(fHClustM02Gamma4TM);
      	fHClustM02Gamma5TM = new TH1F("hClustM02Gamma5TM","",150,0,1);
      	fHClustM02Gamma5TM->SetXTitle("M_{02}");
		fHClustM02Gamma5TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Gamma5TM->Sumw2();
      	fOutput->Add(fHClustM02Gamma5TM);	
      	fHClustM02Gamma6TM = new TH1F("hClustM02Gamma6TM","",150,0,1);
      	fHClustM02Gamma6TM->SetXTitle("M_{02}");
		fHClustM02Gamma6TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Gamma6TM->Sumw2();
      	fOutput->Add(fHClustM02Gamma6TM);			
      	fHClustM02Gamma7TM = new TH1F("hClustM02Gamma7TM","",150,0,1);
      	fHClustM02Gamma7TM->SetXTitle("M_{02}");
		fHClustM02Gamma7TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Gamma7TM->Sumw2();
      	fOutput->Add(fHClustM02Gamma7TM);				
			
      	fHClustM02Pion0TM = new TH1F("hClustM02Pion0TM","",150,0,1);
      	fHClustM02Pion0TM->SetXTitle("M_{02}");
		fHClustM02Pion0TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Pion0TM->Sumw2();	
      	fOutput->Add(fHClustM02Pion0TM);				
      	fHClustM02Pion1TM = new TH1F("hClustM02Pion1TM","",150,0,1);
      	fHClustM02Pion1TM->SetXTitle("M_{02}");
		fHClustM02Pion1TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Pion1TM->Sumw2();
      	fOutput->Add(fHClustM02Pion1TM);
      	fHClustM02Pion2TM = new TH1F("hClustM02Pion2TM","",150,0,1);
      	fHClustM02Pion2TM->SetXTitle("M_{02}");
		fHClustM02Pion2TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Pion2TM->Sumw2();
      	fOutput->Add(fHClustM02Pion2TM);		
      	fHClustM02Pion3TM = new TH1F("hClustM02Pion3TM","",150,0,1);
      	fHClustM02Pion3TM->SetXTitle("M_{02}");
		fHClustM02Pion3TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Pion3TM->Sumw2();	
      	fOutput->Add(fHClustM02Pion3TM);			
      	fHClustM02Pion4TM = new TH1F("hClustM02Pion4TM","",150,0,1);
      	fHClustM02Pion4TM->SetXTitle("M_{02}");
		fHClustM02Pion4TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Pion4TM->Sumw2();	
      	fOutput->Add(fHClustM02Pion4TM);	
      	fHClustM02Pion5TM = new TH1F("hClustM02Pion5TM","",150,0,1);
      	fHClustM02Pion5TM->SetXTitle("M_{02}");
		fHClustM02Pion5TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Pion5TM->Sumw2();	
      	fOutput->Add(fHClustM02Pion5TM);	
      	fHClustM02Pion6TM = new TH1F("hClustM02Pion6TM","",150,0,1);
      	fHClustM02Pion6TM->SetXTitle("M_{02}");
		fHClustM02Pion6TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Pion6TM->Sumw2();	
      	fOutput->Add(fHClustM02Pion6TM);	
      	fHClustM02Pion7TM = new TH1F("hClustM02Pion7TM","",150,0,1);
      	fHClustM02Pion7TM->SetXTitle("M_{02}");
		fHClustM02Pion7TM->SetYTitle("Anzahl Teilchen");
		fHClustM02Pion7TM->Sumw2();	
      	fOutput->Add(fHClustM02Pion7TM);	
								
		
	    fHClustEnergyNCellPhoton = new TH2F("hClustEnergyNCellPhoton","",250,0,20,50,0,50);
	    fHClustEnergyNCellPhoton->SetXTitle("E_{clus}");
	    fHClustEnergyNCellPhoton->SetYTitle("N_{cells}");
	    fOutput->Add(fHClustEnergyNCellPhoton);
	    fHClustNCellM02Photon = new TH2F("hClustNCellM02Photon","",50,0,50,200,0,3);
	    fHClustNCellM02Photon->SetXTitle("N_{cells}");
	    fHClustNCellM02Photon->SetYTitle("M_{02}");
	    fOutput->Add(fHClustNCellM02Photon);		
	    fHClustEnergyNCellPion = new TH2F("hClustEnergyNCellPion","",250,0,20,50,0,50);
	    fHClustEnergyNCellPion->SetXTitle("E_{clus}");
	    fHClustEnergyNCellPion->SetYTitle("N_{cells}");
	    fOutput->Add(fHClustEnergyNCellPion);
	    fHClustNCellM02Pion = new TH2F("hClustNCellM02Pion","",50,0,50,200,0,3);
	    fHClustNCellM02Pion->SetXTitle("N_{cells}");
	    fHClustNCellM02Pion->SetYTitle("M_{02}");
	    fOutput->Add(fHClustNCellM02Pion);		
		
				
		  
	  	
      	fHClustEnergyM02GammaRaw = new TH2F("hClustEnergyM02GammaRaw","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02GammaRaw->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02GammaRaw->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02GammaRaw);
      	fHClustEnergyM02PionRaw = new TH2F("hClustEnergyM02PionRaw","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02PionRaw->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02PionRaw->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02PionRaw);	
      	fHClustEnergyM02ElektronRaw = new TH2F("hClustEnergyM02ElektronRaw","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02ElektronRaw->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02ElektronRaw->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02ElektronRaw);		
      	fHClustM02M20NoGammaRaw = new TH2F("hClustM02M20NoGammaRaw","M02 vs M20",200,0,3,200,0,3);
      	fHClustM02M20NoGammaRaw->SetXTitle("M_{20}");
      	fHClustM02M20NoGammaRaw->SetYTitle("M_{02}");
      	fOutput->Add(fHClustM02M20NoGammaRaw);			
      	fHClustEnergyM02GammaAllRaw = new TH2F("hClustEnergyM02GammaAllRaw","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02GammaAllRaw->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02GammaAllRaw->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02GammaAllRaw);
      	fHClustM02M20GammaAllRaw = new TH2F("hClustM02M20GammaAllRaw","M02 vs M20",200,0,3,200,0,3);
      	fHClustM02M20GammaAllRaw->SetXTitle("M_{20}");
      	fHClustM02M20GammaAllRaw->SetYTitle("M_{02}");
      	fOutput->Add(fHClustM02M20GammaAllRaw);			
      	fHClustEnergyM02AllRaw = new TH2F("hClustEnergyM02AllRaw","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02AllRaw->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02AllRaw->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02AllRaw);	
      	fHClustEnergyM02PionTM = new TH2F("hClustEnergyM02PionTM","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02PionTM->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02PionTM->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02PionTM);	
      	fHClustEnergyM02AllTM = new TH2F("hClustEnergyM02AllTM","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02AllTM->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02AllTM->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02AllTM);	
		
      	fHClustEnergyM02GammaSmallCut = new TH2F("hClustEnergyM02GammaSmallCut","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02GammaSmallCut->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02GammaSmallCut->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02GammaSmallCut);			
      	fHClustEnergyM02PionSmallCut = new TH2F("hClustEnergyM02PionSmallCut","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02PionSmallCut->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02PionSmallCut->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02PionSmallCut);		
		
      	fHClustEnergyM02GammaCell1 = new TH2F("hClustEnergyM02GammaCell1","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02GammaCell1->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02GammaCell1->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02GammaCell1);			
      	fHClustEnergyM02PionCell1 = new TH2F("hClustEnergyM02PionCell1","M02 vs cluster energy",ebins,0,emax,200,0,3);
      	fHClustEnergyM02PionCell1->SetXTitle("E_{C} (GeV)");
      	fHClustEnergyM02PionCell1->SetYTitle("M_{02}");
      	fOutput->Add(fHClustEnergyM02PionCell1);				
		
			

				
						
	    fHClustEnergyNCellPhotonRaw = new TH2F("hClustEnergyNCellPhotonRaw","",250,0,20,50,0,50);
	    fHClustEnergyNCellPhotonRaw->SetXTitle("E_{clus}");
	    fHClustEnergyNCellPhotonRaw->SetYTitle("N_{cells}");
	    fOutput->Add(fHClustEnergyNCellPhotonRaw);
	    fHClustNCellM02PhotonRaw = new TH2F("hClustNCellM02PhotonRaw","",50,0,50,200,0,3);
	    fHClustNCellM02PhotonRaw->SetXTitle("N_{cells}");
	    fHClustNCellM02PhotonRaw->SetYTitle("M_{02}");
	    fOutput->Add(fHClustNCellM02PhotonRaw);		
	    fHClustEnergyNCellPionRaw = new TH2F("hClustEnergyNCellPionRaw","",250,0,20,50,0,50);
	    fHClustEnergyNCellPionRaw->SetXTitle("E_{clus}");
	    fHClustEnergyNCellPionRaw->SetYTitle("N_{cells}");
	    fOutput->Add(fHClustEnergyNCellPionRaw);
	    fHClustNCellM02PionRaw = new TH2F("hClustNCellM02PionRaw","",50,0,50,200,0,3);
	    fHClustNCellM02PionRaw->SetXTitle("N_{cells}");
	    fHClustNCellM02PionRaw->SetYTitle("M_{02}");
	    fOutput->Add(fHClustNCellM02PionRaw);		
		
		
		
		
  		fHParticles = new TH1F("hParticles", "", 6000, -3000, 3000);	
      	fHParticles->SetXTitle("PDG Code");
      	fHParticles->SetYTitle("counts");	
      	fOutput->Add(fHParticles);		
  		fHParticlesRaw = new TH1F("hParticlesRaw", "", 6000, -3000, 3000);	
      	fHParticlesRaw->SetXTitle("PDG Code");
      	fHParticlesRaw->SetYTitle("counts");	
      	fOutput->Add(fHParticlesRaw);		
  		fHParticlesEff = new TH1F("hParticlesEff", "", 6000, -3000, 3000);	
      	fHParticlesEff->SetXTitle("PDG Code");
      	fHParticlesEff->SetYTitle("counts");	
      	fOutput->Add(fHParticlesEff);					
  		fHParticlesCompare = new TH1F("hParticlesCompare", "", 6000, -3000, 3000);	
      	fHParticlesCompare->SetXTitle("PDG Code");
      	fHParticlesCompare->SetYTitle("counts");	
      	fOutput->Add(fHParticlesCompare);					
  		fHParticlesCompare2 = new TH1F("hParticlesCompare2", "", 6000, -3000, 3000);	
      	fHParticlesCompare2->SetXTitle("PDG Code");
      	fHParticlesCompare2->SetYTitle("counts");	
      	fOutput->Add(fHParticlesCompare2);	
  		fHParticlesCompare3 = new TH1F("hParticlesCompare3", "", 6000, -3000, 3000);	
      	fHParticlesCompare3->SetXTitle("PDG Code");
      	fHParticlesCompare3->SetYTitle("counts");	
      	fOutput->Add(fHParticlesCompare3);			
  		fHParticlesTM = new TH1F("hParticlesTM", "", 6000, -3000, 3000);	
      	fHParticlesTM->SetXTitle("PDG Code");
      	fHParticlesTM->SetYTitle("counts");	
      	fOutput->Add(fHParticlesTM);			
  		fHParticlesNoTM = new TH1F("hParticlesNoTM", "", 6000, -3000, 3000);	
      	fHParticlesNoTM->SetXTitle("PDG Code");
      	fHParticlesNoTM->SetYTitle("counts");	
      	fOutput->Add(fHParticlesNoTM);					
		
			
  		fHParticleR = new TH1F("hParticleR", "", 600, 0, 500);	
      	fHParticleR->SetXTitle("R Vertex");
      	fHParticleR->SetYTitle("counts");	
      	fOutput->Add(fHParticleR);			
  		fHParticleRcut = new TH1F("hParticleRcut", "", 600, 0, 500);	
      	fHParticleRcut->SetXTitle("R Vertex");
      	fHParticleRcut->SetYTitle("counts");	
      	fOutput->Add(fHParticleRcut);
  		fHParticleRcutcut = new TH1F("hParticleRcutcut", "", 600, 0, 500);	
      	fHParticleRcutcut->SetXTitle("R Vertex");
      	fHParticleRcutcut->SetYTitle("counts");	
      	fOutput->Add(fHParticleRcutcut);		
  		fHParticleRcutcutcut = new TH1F("hParticleRcutcutcut", "", 600, 0, 500);	
      	fHParticleRcutcutcut->SetXTitle("R Vertex");
      	fHParticleRcutcutcut->SetYTitle("counts");	
      	fOutput->Add(fHParticleRcutcutcut);		
		
  		fHMotherR = new TH1F("hMotherR", "", 10, -5, 5);	
      	fHMotherR->SetXTitle("R Vertex");
      	fHMotherR->SetYTitle("counts");	
      	fOutput->Add(fHMotherR);			
  		fHMotherR2 = new TH1F("hMotherR2", "", 600, 0, 500);	
      	fHMotherR2->SetXTitle("R Vertex");
      	fHMotherR2->SetYTitle("counts");	
      	fOutput->Add(fHMotherR2);				
      	fHClustElectronZR = new TH2F("hClustElectronZR","M02 vs cluster energy",200,0,500,200,-400,400);
      	fHClustElectronZR->SetXTitle("R (cm)");
      	fHClustElectronZR->SetYTitle("Z (cm)");
      	fOutput->Add(fHClustElectronZR);		
		
		
		
		
		fHPtSpecAll = new TH1F("hPtSpecAll","", 50, 0.01, 14.9);
		fHPtSpecAll->SetXTitle("p_{T} (GeV/c)");
		fHPtSpecAll->SetYTitle("counts");	
		fHPtSpecAll->Sumw2();
		fOutput->Add(fHPtSpecAll);
    	fHPtSpecGamma = new TH1F("hPtSpecGamma","", 50, 0.01, 14.9);
    	fHPtSpecGamma->SetXTitle("p_{T} (GeV/c)");
    	fHPtSpecGamma->SetYTitle("counts");	
		fHPtSpecGamma->Sumw2();		
    	fOutput->Add(fHPtSpecGamma);			  
    	fHPtSpecPion = new TH1F("hPtSpecPion","", 50, 0.01, 14.9);
   		fHPtSpecPion->SetXTitle("p_{T} (GeV/c)");
    	fHPtSpecPion->SetYTitle("counts");	
		fHPtSpecPion->Sumw2();	
		fOutput->Add(fHPtSpecPion);	
   		fHPtSpecElectron = new TH1F("hPtSpecElectron","", 50, 0.01, 14.9);
   		fHPtSpecElectron->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecElectron->SetYTitle("counts");
		fHPtSpecElectron->Sumw2();	
   		fOutput->Add(fHPtSpecElectron);  		
   		fHPtSpecMyon = new TH1F("hPtSpecMyon","", 50, 0.01, 14.9);
   		fHPtSpecMyon->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecMyon->SetYTitle("counts");
		fHPtSpecMyon->Sumw2();	
   		fOutput->Add(fHPtSpecMyon);  			
   		fHPtSpecProton = new TH1F("hPtSpecProton","", 50, 0.01, 14.9);
   		fHPtSpecProton->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecProton->SetYTitle("counts");
		fHPtSpecProton->Sumw2();	
   		fOutput->Add(fHPtSpecProton);  
   		fHPtSpecNeutron = new TH1F("hPtSpecNeutron","", 50, 0.01, 14.9);
   		fHPtSpecNeutron->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecNeutron->SetYTitle("counts");
		fHPtSpecNeutron->Sumw2();	
   		fOutput->Add(fHPtSpecNeutron);  
   		fHPtSpecKaon = new TH1F("hPtSpecKaon","", 50, 0.01, 14.9);
   		fHPtSpecKaon->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecKaon->SetYTitle("counts");
		fHPtSpecKaon->Sumw2();	
   		fOutput->Add(fHPtSpecKaon);  
   		fHPtSpecKaon0 = new TH1F("hPtSpecKaon0","", 50, 0.01, 14.9);
   		fHPtSpecKaon0->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecKaon0->SetYTitle("counts");
		fHPtSpecKaon0->Sumw2();	
   		fOutput->Add(fHPtSpecKaon0);  				
   		fHPtSpecNoGamma = new TH1F("hPtSpecNoGamma","", 50, 0.01, 14.9);
   		fHPtSpecNoGamma->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecNoGamma->SetYTitle("counts");
		fHPtSpecNoGamma->Sumw2();	
   		fOutput->Add(fHPtSpecNoGamma); 
   		fHPtSpecCharged = new TH1F("hPtSpecCharged","", 50, 0.01, 14.9);
   		fHPtSpecCharged->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecCharged->SetYTitle("counts");
		fHPtSpecCharged->Sumw2();	
   		fOutput->Add(fHPtSpecCharged); 		
		

	
   		fHPtSpecAllRaw = new TH1F("hPtSpecAllRaw","", 50, 0.01, 14.9);
   		fHPtSpecAllRaw->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecAllRaw->SetYTitle("counts");	
		fHPtSpecAllRaw->Sumw2();			
   		fOutput->Add(fHPtSpecAllRaw);
   		fHPtSpecGammaRaw = new TH1F("hPtSpecGammaRaw","", 50, 0.01, 14.9);
   		fHPtSpecGammaRaw->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecGammaRaw->SetYTitle("counts");	
		fHPtSpecGammaRaw->Sumw2();
   		fOutput->Add(fHPtSpecGammaRaw);			  
   		fHPtSpecPionRaw = new TH1F("hPtSpecPionRaw","", 50, 0.01, 14.9);
   		fHPtSpecPionRaw->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecPionRaw->SetYTitle("counts");
		fHPtSpecPionRaw->Sumw2();	
   		fOutput->Add(fHPtSpecPionRaw);  
   		fHPtSpecElectronRaw = new TH1F("hPtSpecElectronRaw","", 50, 0.01, 14.9);
   		fHPtSpecElectronRaw->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecElectronRaw->SetYTitle("counts");
		fHPtSpecElectronRaw->Sumw2();	
   		fOutput->Add(fHPtSpecElectronRaw);  		
   		fHPtSpecMyonRaw = new TH1F("hPtSpecMyonRaw","", 50, 0.01, 14.9);
   		fHPtSpecMyonRaw->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecMyonRaw->SetYTitle("counts");
		fHPtSpecMyonRaw->Sumw2();	
   		fOutput->Add(fHPtSpecMyonRaw);  			
   		fHPtSpecProtonRaw = new TH1F("hPtSpecProtonRaw","", 50, 0.01, 14.9);
   		fHPtSpecProtonRaw->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecProtonRaw->SetYTitle("counts");
		fHPtSpecProtonRaw->Sumw2();	
   		fOutput->Add(fHPtSpecProtonRaw);  
   		fHPtSpecNeutronRaw = new TH1F("hPtSpecNeutronRaw","", 50, 0.01, 14.9);
   		fHPtSpecNeutronRaw->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecNeutronRaw->SetYTitle("counts");
		fHPtSpecNeutronRaw->Sumw2();	
   		fOutput->Add(fHPtSpecNeutronRaw);  
   		fHPtSpecKaonRaw = new TH1F("hPtSpecKaonRaw","", 50, 0.01, 14.9);
   		fHPtSpecKaonRaw->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecKaonRaw->SetYTitle("counts");
		fHPtSpecKaonRaw->Sumw2();	
   		fOutput->Add(fHPtSpecKaonRaw);  
   		fHPtSpecKaon0Raw = new TH1F("hPtSpecKaon0Raw","", 50, 0.01, 14.9);
   		fHPtSpecKaon0Raw->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecKaon0Raw->SetYTitle("counts");
		fHPtSpecKaon0Raw->Sumw2();	
   		fOutput->Add(fHPtSpecKaon0Raw);  				
   		fHPtSpecNoGammaRaw = new TH1F("hPtSpecNoGammaRaw","", 50, 0.01, 14.9);
   		fHPtSpecNoGammaRaw->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecNoGammaRaw->SetYTitle("counts");
		fHPtSpecNoGammaRaw->Sumw2();	
   		fOutput->Add(fHPtSpecNoGammaRaw);  		
   		fHPtSpecChargedRaw = new TH1F("hPtSpecCharged","", 50, 0.01, 14.9);
   		fHPtSpecChargedRaw->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecChargedRaw->SetYTitle("counts");
		fHPtSpecChargedRaw->Sumw2();	
   		fOutput->Add(fHPtSpecChargedRaw); 	
   		fHPtSpecChargedTM = new TH1F("hPtSpecChargedTM","", 50, 0.01, 14.9);
   		fHPtSpecChargedTM->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecChargedTM->SetYTitle("counts");
		fHPtSpecChargedTM->Sumw2();	
   		fOutput->Add(fHPtSpecChargedTM); 					
		
   		fHPtSpecElectronTM = new TH1F("hPtSpecElectronTM","", 50, 0.01, 14.9);
   		fHPtSpecElectronTM->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecElectronTM->SetYTitle("counts");
		fHPtSpecElectronTM->Sumw2();	
   		fOutput->Add(fHPtSpecElectronTM); 	
   		fHPtSpecPionTM = new TH1F("hPtSpecPionTM","", 50, 0.01, 14.9);
   		fHPtSpecPionTM->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecPionTM->SetYTitle("counts");
		fHPtSpecPionTM->Sumw2();	
   		fOutput->Add(fHPtSpecPionTM); 	
   		fHPtSpecElectronNoTM = new TH1F("hPtSpecElectronNoTM","", 50, 0.01, 14.9);
   		fHPtSpecElectronNoTM->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecElectronNoTM->SetYTitle("counts");
		fHPtSpecElectronNoTM->Sumw2();	
   		fOutput->Add(fHPtSpecElectronNoTM); 	
   		fHPtSpecPionNoTM = new TH1F("hPtSpecPionNoTM","", 50, 0.01, 14.9);
   		fHPtSpecPionNoTM->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecPionNoTM->SetYTitle("counts");
		fHPtSpecPionNoTM->Sumw2();	
   		fOutput->Add(fHPtSpecPionNoTM); 	
   		fHPtSpecPionNoM02 = new TH1F("hPtSpecPionNoM02","", 50, 0.01, 14.9);
   		fHPtSpecPionNoM02->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecPionNoM02->SetYTitle("counts");
		fHPtSpecPionNoM02->Sumw2();	
   		fOutput->Add(fHPtSpecPionNoM02); 			
   		fHPtSpecGammaNoM02 = new TH1F("hPtSpecGammaNoM02","", 50, 0.01, 14.9);
   		fHPtSpecGammaNoM02->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecGammaNoM02->SetYTitle("counts");
		fHPtSpecGammaNoM02->Sumw2();	
   		fOutput->Add(fHPtSpecGammaNoM02); 					
   		fHPtSpecPionM02 = new TH1F("hPtSpecPionM02","", 50, 0.01, 14.9);
   		fHPtSpecPionM02->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecPionM02->SetYTitle("counts");
		fHPtSpecPionM02->Sumw2();	
   		fOutput->Add(fHPtSpecPionM02); 			
   		fHPtSpecGammaM02 = new TH1F("hPtSpecGammaM02","", 50, 0.01, 14.9);
   		fHPtSpecGammaM02->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecGammaM02->SetYTitle("counts");
		fHPtSpecGammaM02->Sumw2();	
   		fOutput->Add(fHPtSpecGammaM02); 		
   		fHPtSpecGammaTM = new TH1F("hPtSpecGammaTM","", 50, 0.01, 14.9);
   		fHPtSpecGammaTM->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecGammaTM->SetYTitle("counts");
		fHPtSpecGammaTM->Sumw2();	
   		fOutput->Add(fHPtSpecGammaTM); 
				
		
   		fHPtSpecGammaCompare = new TH1F("hPtSpecGammaCompare","", 50, 0.01, 14.9);
   		fHPtSpecGammaCompare->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecGammaCompare->SetYTitle("counts");	
		fHPtSpecGammaCompare->Sumw2();
   		fOutput->Add(fHPtSpecGammaCompare);			  
   		fHPtSpecPionCompare = new TH1F("hPtSpecPionCompare","", 50, 0.01, 14.9);
   		fHPtSpecPionCompare->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecPionCompare->SetYTitle("counts");
		fHPtSpecPionCompare->Sumw2();	
   		fOutput->Add(fHPtSpecPionCompare);  
   		fHPtSpecCompare = new TH1F("hPtSpecCompare","", 50, 0.01, 14.9);
   		fHPtSpecCompare->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecCompare->SetYTitle("counts");
		fHPtSpecCompare->Sumw2();	
   		fOutput->Add(fHPtSpecCompare);  		
   		fHPtSpecElectronMerge = new TH1F("hPtSpecElectronMerge","", 50, 0.01, 14.9);
   		fHPtSpecElectronMerge->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecElectronMerge->SetYTitle("counts");
		fHPtSpecElectronMerge->Sumw2();	
   		fOutput->Add(fHPtSpecElectronMerge);  	
			
		
				
   		fHPtSpecEffPhoton = new TH1F("hPtSpecEffPhoton","", 50, 0.01, 14.9);
   		fHPtSpecEffPhoton->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecEffPhoton->SetYTitle("counts");
		fHPtSpecEffPhoton->Sumw2();	
   		fOutput->Add(fHPtSpecEffPhoton);
   		fHPtSpecAccPhoton = new TH1F("hPtSpecAccPhoton","", 50, 0.01, 14.9);
   		fHPtSpecAccPhoton->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecAccPhoton->SetYTitle("counts");
		fHPtSpecAccPhoton->Sumw2();	
   		fOutput->Add(fHPtSpecAccPhoton);		
   		fHPtSpecEffParticle = new TH1F("hPtSpecEffParticle","", 50, 0.01, 14.9);
   		fHPtSpecEffParticle->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecEffParticle->SetYTitle("counts");
		fHPtSpecEffParticle->Sumw2();	
   		fOutput->Add(fHPtSpecEffParticle);	
   		fHPtSpecEffCluster = new TH1F("hPtSpecEffCluster","", 50, 0.01, 14.9);
   		fHPtSpecEffCluster->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecEffCluster->SetYTitle("counts");
		fHPtSpecEffCluster->Sumw2();	
   		fOutput->Add(fHPtSpecEffCluster);
   		fHPtSpecConversion = new TH1F("hPtSpecConversion","", 50, 0.01, 14.9);
   		fHPtSpecConversion->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecConversion->SetYTitle("counts");
		fHPtSpecConversion->Sumw2();	
   		fOutput->Add(fHPtSpecConversion);	
   		fHPtSpecConversionNot = new TH1F("hPtSpecConversionNot","", 50, 0.01, 14.9);
   		fHPtSpecConversionNot->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecConversionNot->SetYTitle("counts");
		fHPtSpecConversionNot->Sumw2();	
   		fOutput->Add(fHPtSpecConversionNot);			
   		fHPtSpecEffPhotonEta5 = new TH1F("hPtSpecEffPhotonEta5","", 50, 0.01, 14.9);
   		fHPtSpecEffPhotonEta5->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecEffPhotonEta5->SetYTitle("counts");
		fHPtSpecEffPhotonEta5->Sumw2();	
   		fOutput->Add(fHPtSpecEffPhotonEta5);		
   		fHPtSpecEffPhotonEta4 = new TH1F("hPtSpecEffPhotonEta4","", 50, 0.01, 14.9);
   		fHPtSpecEffPhotonEta4->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecEffPhotonEta4->SetYTitle("counts");
		fHPtSpecEffPhotonEta4->Sumw2();	
   		fOutput->Add(fHPtSpecEffPhotonEta4);		
   		fHPtSpecEffPhotonEta3 = new TH1F("hPtSpecEffPhotonEta3","", 50, 0.01, 14.9);
   		fHPtSpecEffPhotonEta3->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecEffPhotonEta3->SetYTitle("counts");
		fHPtSpecEffPhotonEta3->Sumw2();	
   		fOutput->Add(fHPtSpecEffPhotonEta3);		
   		fHPtSpecEffPhotonEta2 = new TH1F("hPtSpecEffPhotonEta2","", 50, 0.01, 14.9);
   		fHPtSpecEffPhotonEta2->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecEffPhotonEta2->SetYTitle("counts");
		fHPtSpecEffPhotonEta2->Sumw2();	
   		fOutput->Add(fHPtSpecEffPhotonEta2);		
   		fHPtSpecEffPhotonEta1 = new TH1F("hPtSpecEffPhotonEta1","", 50, 0.01, 14.9);
   		fHPtSpecEffPhotonEta1->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecEffPhotonEta1->SetYTitle("counts");
		fHPtSpecEffPhotonEta1->Sumw2();	
   		fOutput->Add(fHPtSpecEffPhotonEta1);						

   		fHPtSpecDecayPi0 = new TH1F("hPtSpecDecayPi0","", 50, 0.01, 14.9);
   		fHPtSpecDecayPi0->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecDecayPi0->SetYTitle("counts");
		fHPtSpecDecayPi0->Sumw2();	
   		fOutput->Add(fHPtSpecDecayPi0);	
   		fHPtSpecDecayEta = new TH1F("hPtSpecDecayEta","", 50, 0.01, 14.9);
   		fHPtSpecDecayEta->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecDecayEta->SetYTitle("counts");
		fHPtSpecDecayEta->Sumw2();	
   		fOutput->Add(fHPtSpecDecayEta);	
   		fHPtSpecDecayOmega = new TH1F("hPtSpecDecayOmega","", 50, 0.01, 14.9);
   		fHPtSpecDecayOmega->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecDecayOmega->SetYTitle("counts");
		fHPtSpecDecayOmega->Sumw2();	
   		fOutput->Add(fHPtSpecDecayOmega);	
   		fHPtSpecDecayEtap = new TH1F("hPtSpecDecayEtap","", 50, 0.01, 14.9);
   		fHPtSpecDecayEtap->SetXTitle("p_{T} (GeV/c)");
   		fHPtSpecDecayEtap->SetYTitle("counts");
		fHPtSpecDecayEtap->Sumw2();	
   		fOutput->Add(fHPtSpecDecayEtap);	

//   		fHPtSpecM02Cut0Photon = new TH1F("hPtSpecM02Cut0Photon","", 50, 0.01, 14.9);
//   		fHPtSpecM02Cut0Photon->SetXTitle("p_{T} (GeV/c)");
//   		fHPtSpecM02Cut0Photon->SetYTitle("counts");	
//		fHPtSpecM02Cut0Photon->Sumw2();			
//   		fOutput->Add(fHPtSpecM02Cut0Photon);	
//   		fHPtSpecM02Cut1Pion = new TH1F("hPtSpecM02Cut1Pion","", 50, 0.01, 14.9);
//   		fHPtSpecM02Cut1Pion->SetXTitle("p_{T} (GeV/c)");
//   		fHPtSpecM02Cut1Pion->SetYTitle("counts");	
//		fHPtSpecM02Cut1Pion->Sumw2();			
//   		fOutput->Add(fHPtSpecM02Cut1Pion);
//   		fHPtSpecM02Cut1Photon = new TH1F("hPtSpecM02Cut1Photon","", 50, 0.01, 14.9);
//   		fHPtSpecM02Cut1Photon->SetXTitle("p_{T} (GeV/c)");
//   		fHPtSpecM02Cut1Photon->SetYTitle("counts");	
//		fHPtSpecM02Cut1Photon->Sumw2();			
//   		fOutput->Add(fHPtSpecM02Cut1Photon);	
   		//fHPtSpecM02Cut2Pion = new TH1F("hPtSpecM02Cut2Pion","", 50, 0.01, 14.9);
   		//fHPtSpecM02Cut2Pion->SetXTitle("p_{T} (GeV/c)");
   		//fHPtSpecM02Cut2Pion->SetYTitle("counts");	
		//fHPtSpecM02Cut2Pion->Sumw2();			
   		//fOutput->Add(fHPtSpecM02Cut2Pion);
   		//fHPtSpecM02Cut2Photon = new TH1F("hPtSpecM02Cut2Photon","", 50, 0.01, 14.9);
   		//fHPtSpecM02Cut2Photon->SetXTitle("p_{T} (GeV/c)");
   		//fHPtSpecM02Cut2Photon->SetYTitle("counts");	
		//fHPtSpecM02Cut2Photon->Sumw2();			
   		//fOutput->Add(fHPtSpecM02Cut2Photon);					
		
//   		fHPtSpecTrackCut0Pion = new TH1F("hPtSpecTrackCut0Pion","", 50, 0.01, 14.9);
//   		fHPtSpecTrackCut0Pion->SetXTitle("p_{T} (GeV/c)");
//   		fHPtSpecTrackCut0Pion->SetYTitle("counts");	
//		fHPtSpecTrackCut0Pion->Sumw2();			
//   		fOutput->Add(fHPtSpecTrackCut0Pion);
//   		fHPtSpecTrackCut0Photon = new TH1F("hPtSpecTrackCut0Photon","", 50, 0.01, 14.9);
//   		fHPtSpecTrackCut0Photon->SetXTitle("p_{T} (GeV/c)");
//   		fHPtSpecTrackCut0Photon->SetYTitle("counts");	
//		fHPtSpecTrackCut0Photon->Sumw2();			
//   		fOutput->Add(fHPtSpecTrackCut0Photon);	
//   		fHPtSpecTrackCut1Pion = new TH1F("hPtSpecTrackCut1Pion","", 50, 0.01, 14.9);
//   		fHPtSpecTrackCut1Pion->SetXTitle("p_{T} (GeV/c)");
//   		fHPtSpecTrackCut1Pion->SetYTitle("counts");	
//		fHPtSpecTrackCut1Pion->Sumw2();			
//   		fOutput->Add(fHPtSpecTrackCut1Pion);
//   		fHPtSpecTrackCut1Photon = new TH1F("hPtSpecTrackCut1Photon","", 50, 0.01, 14.9);
//   		fHPtSpecTrackCut1Photon->SetXTitle("p_{T} (GeV/c)");
//   		fHPtSpecTrackCut1Photon->SetYTitle("counts");	
//		fHPtSpecTrackCut1Photon->Sumw2();			
//   		fOutput->Add(fHPtSpecTrackCut1Photon);	
   		//fHPtSpecTrackCut2Pion = new TH1F("hPtSpecTrackCut2Pion","", 50, 0.01, 14.9);
   		//fHPtSpecTrackCut2Pion->SetXTitle("p_{T} (GeV/c)");
   		//fHPtSpecTrackCut2Pion->SetYTitle("counts");	
		//fHPtSpecTrackCut2Pion->Sumw2();			
   		//fOutput->Add(fHPtSpecTrackCut2Pion);
   		//fHPtSpecTrackCut2Photon = new TH1F("hPtSpecTrackCut2Photon","", 50, 0.01, 14.9);
   		//fHPtSpecTrackCut2Photon->SetXTitle("p_{T} (GeV/c)");
   		//fHPtSpecTrackCut2Photon->SetYTitle("counts");	
		//fHPtSpecTrackCut2Photon->Sumw2();			
   		//fOutput->Add(fHPtSpecTrackCut2Photon);	
		
  
    // histogram for conversion point
    fHConversionPoint = new TH2F("hConversionPoint","conversion point in xy",1000,-500,500,1000,-500,500);
    fHConversionPoint->SetXTitle("x");
    fHConversionPoint->SetYTitle("y");
    fOutput->Add(fHConversionPoint);
  }
  
  TH1::SetDefaultSumw2(th1);
  TH2::SetDefaultSumw2(th2);
  PostData(1, fOutput);
  
  //  const Int_t nbins = 20;
  //  Double_t xbins[nbins+1] = {0.,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12.5,15,20,25,50};
  //  const Int_t nbins = 300;
  //  const Int_t ptmax = 30;
  
  // MC histograms
  if(fMcMode){
    
    if(fAddedSignal){
      fHAddPionEtaPt = new TH2F("hAddPionEtaPt","",150,-2.5,2.5,100,0.,10.);
      fHAddPionEtaPt->SetXTitle("#eta_{#gamma#gamma}");
      fHAddPionEtaPt->SetYTitle("p_{T}");
      fOutput->Add(fHAddPionEtaPt);

      fHAddPionEtaPtWgt = new TH2F("hAddPionEtaPtWgt","",150,-2.5,2.5,100,0.,10.);
      fHAddPionEtaPtWgt->SetXTitle("#eta_{#gamma#gamma}");
      fHAddPionEtaPtWgt->SetYTitle("p_{T}");
      fOutput->Add(fHAddPionEtaPtWgt);
    }
    
    fHPyPionEtaPt = new TH2F("hPyPionEtaPt","",150,-2.5,2.5,100,0.,10.);
    fHPyPionEtaPt->SetXTitle("#eta_{#gamma#gamma}");
    fHPyPionEtaPt->SetYTitle("p_{T}");
    fOutput->Add(fHPyPionEtaPt);

    // pi0
    fHPionTruthPt = new TH1F("hPionTruthPt","pi0 truth pT from MC",nbins,0,ptmax);
    fHPionTruthPt->SetXTitle("p_{T}");
    fOutput->Add(fHPionTruthPt);
    
    fHPionTruthPtIn = new TH1F("hPionTruthPtIn","pi0 truth pT from MC within eta range",nbins,0,ptmax);
    fHPionTruthPtIn->SetXTitle("p_{T}");
    fOutput->Add(fHPionTruthPtIn);
    
    fHPionTruthPtAcc = new TH1F("hPionTruthPtAcc","pi0 truth pT from MC within acceptance",nbins,0,ptmax);
    fHPionTruthPtAcc->SetXTitle("p_{T}");
    fOutput->Add(fHPionTruthPtAcc);
    
    fHPionTruthPtConvAcc = new TH1F("hPionTruthPtConvAcc","pi0 truth pT from MC within acceptance (eta<0.8 for converted)",nbins,0,ptmax);
    fHPionTruthPtConvAcc->SetXTitle("p_{T}");
    fOutput->Add(fHPionTruthPtConvAcc);
    
    // eta
    fHEtaTruthPt = new TH1F("hEtaTruthPt","eta truth pT from MC",nbins,0,ptmax);
    fHEtaTruthPt->SetXTitle("p_{T}");
    fOutput->Add(fHEtaTruthPt);
    
    fHEtaTruthPtIn = new TH1F("hEtaTruthPtIn","eta truth pT from MC within eta range",nbins,0,ptmax);
    fHEtaTruthPtIn->SetXTitle("p_{T}");
    fOutput->Add(fHEtaTruthPtIn);
    
    fHEtaTruthPtAcc = new TH1F("hEtaTruthPtAcc","eta truth pT from MC within acceptance",nbins,0,ptmax);
    fHEtaTruthPtAcc->SetXTitle("p_{T}");
    fOutput->Add(fHEtaTruthPtAcc);
    
    fHEtaTruthPtConvAcc = new TH1F("hEtaTruthPtConvAcc","eta truth pT from MC within acceptance (eta<0.8 for converted)",nbins,0,ptmax);
    fHEtaTruthPtConvAcc->SetXTitle("p_{T}");
    fOutput->Add(fHEtaTruthPtConvAcc);
    
    // direct photons
    fHGamTruthPt = new TH1F("hGamTruthPt","gamma truth pT from MC",nbins,0,ptmax);
    fHGamTruthPt->SetXTitle("p_{T}");
    fOutput->Add(fHGamTruthPt);
    
    fHGamTruthPtIn = new TH1F("hGamTruthPtIn","gamma truth pT from MC within eta range",nbins,0,ptmax);
    fHGamTruthPtIn->SetXTitle("p_{T}");
    fOutput->Add(fHGamTruthPtIn);
    
    fHGamTruthPtAcc = new TH1F("hGamTruthPtAcc","gamma truth pT from MC within acceptance",nbins,0,ptmax);
    fHGamTruthPtAcc->SetXTitle("p_{T}");
    fOutput->Add(fHGamTruthPtAcc);
    
    // added signals
    // pi0
    fHPionTruthPtAdd = new TH1F("hPionTruthPtAdd","added pi0 truth pT from MC",nbins,0,ptmax);
    fHPionTruthPtAdd->SetXTitle("p_{T}");
    fHPionTruthPtAdd->Sumw2();
    fOutput->Add(fHPionTruthPtAdd);
    
    fHPionTruthPtInAdd = new TH1F("hPionTruthPtInAdd","added pi0 truth pT from MC within eta range",nbins,0,ptmax);
    fHPionTruthPtInAdd->SetXTitle("p_{T}");
    fHPionTruthPtInAdd->Sumw2();
    fOutput->Add(fHPionTruthPtInAdd);
    
    fHPionTruthPtAccAdd = new TH1F("hPionTruthPtAccAdd","added pi0 truth pT from MC within acceptance",nbins,0,ptmax);
    fHPionTruthPtAccAdd->SetXTitle("p_{T}");
    fHPionTruthPtAccAdd->Sumw2();
    fOutput->Add(fHPionTruthPtAccAdd);
    
    fHPionTruthPtConvAccAdd = new TH1F("hPionTruthPtConvAccAdd","added pi0 truth pT from MC within acceptance (eta<0.8 for converted)",nbins,0,ptmax);
    fHPionTruthPtConvAccAdd->SetXTitle("p_{T}");
    fOutput->Add(fHPionTruthPtConvAccAdd);
    
    // eta
    fHEtaTruthPtAdd = new TH1F("hEtaTruthPtAdd","added eta truth pT from MC",nbins,0,ptmax);
    fHEtaTruthPtAdd->SetXTitle("p_{T}");
    fHEtaTruthPtAdd->Sumw2();
    fOutput->Add(fHEtaTruthPtAdd);
    
    fHEtaTruthPtInAdd = new TH1F("hEtaTruthPtInAdd","added eta truth pT from MC within eta range",nbins,0,ptmax);
    fHEtaTruthPtInAdd->SetXTitle("p_{T}");
    fHEtaTruthPtInAdd->Sumw2();
    fOutput->Add(fHEtaTruthPtInAdd);
    
    fHEtaTruthPtAccAdd = new TH1F("hEtaTruthPtAccAdd","added eta truth pT from MC within acceptance",nbins,0,ptmax);
    fHEtaTruthPtAccAdd->SetXTitle("p_{T}");
    fHEtaTruthPtAccAdd->Sumw2();
    fOutput->Add(fHEtaTruthPtAccAdd);
    
    fHEtaTruthPtConvAccAdd = new TH1F("hEtaTruthPtConvAccAdd","added eta truth pT from MC within acceptance (eta<0.8 for converted)",nbins,0,ptmax);
    fHEtaTruthPtConvAccAdd->SetXTitle("p_{T}");
    fOutput->Add(fHEtaTruthPtConvAccAdd);
    
    // particle information
    fHMCpartfrac = new TH2F("hMCpartfrac","fraction of most energy MC particle vs. fraction of this particle to all",100,0,2,100,0.8,1.2);
    fHMCpartfrac->SetXTitle("E_{MC}/E_{Clu}");
    fHMCpartfrac->SetYTitle("E_{MC}^{highest}/E_{MC}^{all}");
    fOutput->Add(fHMCpartfrac);
    
    fHECluEMC = new TH2F("hECluEMC","energy of most contributing MC track vs. energy of cluster",200,0,10,200,0,10);
    fHECluEMC->SetXTitle("E_{MC}");
    fHECluEMC->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMC);
    
    fHECluEMCAddPi0 = new TH2F("hECluEMCAddPi0","energy of most contributing added pi0 MC track vs. energy of cluster",200,0,10,200,0,10);
    fHECluEMCAddPi0->SetXTitle("E_{MC}");
    fHECluEMCAddPi0->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCAddPi0);
    
    fHECluEMCAddEta = new TH2F("hECluEMCAddEta","energy of most contributing added eta MC track vs. energy of cluster",200,0,10,200,0,10);
    fHECluEMCAddEta->SetXTitle("E_{MC}");
    fHECluEMCAddEta->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCAddEta);

//    fHRecTrue = new TH2F("hRecTrue","truth energy/cluster energy vs. cluster energy",200,0,10,200,0,10);
//    fHRecTrue->SetXTitle("E_{clu}");
//    fHRecTrue->SetYTitle("E_{tru}/E_{clu}");
//    fOutput->Add(fHRecTrue);
//    
//    fHRecTrueAddPi0 = new TH2F("hRecTrueAddPi0","truth energy/cluster energy vs. cluster energy for added pi0 MC",200,0,10,200,0,10);
//    fHRecTrueAddPi0->SetXTitle("E_{clu}");
//    fHRecTrueAddPi0->SetYTitle("E_{tru}/E_{clu}");
//    fOutput->Add(fHRecTrueAddPi0);
//
//    fHRecTrueAddEta = new TH2F("hRecTrueAddEta","truth energy/cluster energy vs. cluster energy for added eta MC",200,0,10,200,0,10);
//    fHRecTrueAddEta->SetXTitle("E_{clu}");
//    fHRecTrueAddEta->SetYTitle("E_{tru}/E_{clu}");
//    fOutput->Add(fHRecTrueAddEta);
    
    // also for clusters with more than one contribution
    //    fHMCpartfracnofull = new TH2F("hMCpartfracnofull","fraction of most energy MC particle vs. fraction of this particle to all",100,0,2,100,0.8,1.2);
    //    fHMCpartfracnofull->SetXTitle("E_{MC}/E_{Clu}");
    //    fHMCpartfracnofull->SetYTitle("E_{MC}^{highest}/E_{MC}^{all}");
    //    fOutput->Add(fHMCpartfracnofull);
    
    fHECluEMCnofullAdd = new TH2F("hECluEMCnofullAdd","energy of most contributing added MC track vs. energy of cluster (more than one MC particle in cluster)",200,0,10,200,0,10);
    fHECluEMCnofullAdd->SetXTitle("E_{MC}");
    fHECluEMCnofullAdd->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCnofullAdd);
    
    fHECluEMCnofull = new TH2F("hECluEMCnofull","energy of most contributing MC track vs. energy of cluster (more than one MC particle in cluster)",200,0,10,200,0,10);
    fHECluEMCnofull->SetXTitle("E_{MC}");
    fHECluEMCnofull->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCnofull);
    
    fHECluEMCelectron = new TH2F("hECluEMCelectron","energy of most contributing MC track vs. energy of cluster for electrons",200,0,10,200,0,10);
    fHECluEMCelectron->SetXTitle("E_{MC}");
    fHECluEMCelectron->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCelectron);
    
    fHECluEMCpion = new TH2F("hECluEMCpion","energy of most contributing MC track vs. energy of cluster for pions",200,0,10,200,0,10);
    fHECluEMCpion->SetXTitle("E_{MC}");
    fHECluEMCpion->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCpion);
    
    fHECluEMCkaon = new TH2F("hECluEMCkaon","energy of most contributing MC track vs. energy of cluster for kaons",200,0,10,200,0,10);
    fHECluEMCkaon->SetXTitle("E_{MC}");
    fHECluEMCkaon->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCkaon);
    
    fHECluEMCother = new TH2F("hECluEMCother","energy of most contributing MC track vs. energy of cluster for others",200,0,10,200,0,10);
    fHECluEMCother->SetXTitle("E_{MC}");
    fHECluEMCother->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCother);
    
    fHECluEMCpi0single = new TH2F("hECluEMCpi0single","E pi0 in MC vs. E of merged cluster",200,5,15,200,5,15);
    fHECluEMCpi0single->SetXTitle("p_{T}^{MC}");
    fHECluEMCpi0single->SetYTitle("p_{T}^{Clu}");
    fOutput->Add(fHECluEMCpi0single);
    
    fHNMothers = new TH2F("hNMothers","fHNMothers",21,-0.5,20.5,60,0,30);
    fHNMothers->SetXTitle("mother iterations");
    fHNMothers->SetYTitle("Energy");
    fOutput->Add(fHNMothers);

    fHClustNcell = new TH1F("hClustNcell"," ",21,-0.5,20.5);
    fHClustNcell->SetXTitle("number of cells in a cluster");
    fHClustNcell->SetYTitle("counts");
    fOutput->Add(fHClustNcell);	
    fHClustNcell1 = new TH1F("hClustNcell1"," ",50,0,5);
    fHClustNcell1->SetXTitle("number of cells in a cluster");
    fHClustNcell1->SetYTitle("counts");
    fOutput->Add(fHClustNcell1);		
	

    fHClustNcellPhoton = new TH1F("hClustNcellPhoton"," ",21,-0.5,20.5);
    fHClustNcellPhoton->SetXTitle("number of cells in a cluster");
    fHClustNcellPhoton->SetYTitle("counts");
    fOutput->Add(fHClustNcellPhoton);	
	
    fHClustNcellNoPhoton = new TH1F("hClustNcellNoPhoton"," ",21,-0.5,20.5);
    fHClustNcellNoPhoton->SetXTitle("number of cells in a cluster");
    fHClustNcellNoPhoton->SetYTitle("counts");
    fOutput->Add(fHClustNcellNoPhoton);		
    fHClustNcellPhotonCut = new TH1F("hClustNcellPhotonCut"," ",21,-0.5,20.5);
    fHClustNcellPhotonCut->SetXTitle("number of cells in a cluster");
    fHClustNcellPhotonCut->SetYTitle("counts");
    fOutput->Add(fHClustNcellPhotonCut);	
    fHClustNcellNoPhotonCut = new TH1F("hClustNcellNoPhotonCut"," ",21,-0.5,20.5);
    fHClustNcellNoPhotonCut->SetXTitle("number of cells in a cluster");
    fHClustNcellNoPhotonCut->SetYTitle("counts");
    fOutput->Add(fHClustNcellNoPhotonCut);		
	
	fHClustEP = new TH1F("hClustEP", "", 100,0, 3);	
  	fHClustEP->SetXTitle("E/p");
  	fHClustEP->SetYTitle("counts");	
  	fOutput->Add(fHClustEP);	
	
    fHGammaMIP = new TH1F("hGammaMIP"," ",100,0,1);
    fHGammaMIP->SetXTitle("E_{C} (GeV)");
    fHGammaMIP->SetYTitle("counts");
    fOutput->Add(fHGammaMIP);					

    fHHadronMIP = new TH1F("hHadronMIP"," ",100,0,1);
    fHHadronMIP->SetXTitle("E_{C} (GeV)");
    fHHadronMIP->SetYTitle("counts");
    fOutput->Add(fHHadronMIP);	

    
  }
  
//	if(fRotateMixed){
//		// more histograms
//		fHMixRotation = new TH1F("hMixRotation","rotation angle of mixed events",100,-6.28,6.28);
//		fHMixRotation->SetXTitle("phi");
//		fOutput->Add(fHMixRotation);
//	}
//  fHCorrection = new TH1F("hCorrection","correction factor for single photon",100,0,2);
//  fHCorrection->SetXTitle("correction factor");
//  fOutput->Add(fHCorrection);

//  fHPionSm = new TH2F("hPionSm","mass shift due to energy scale",200,0,0.5,200,0,0.5);
//  fHPionSm->SetXTitle("original m_inv");
//  fHPionSm->SetYTitle("changed m_inv");
//  fOutput->Add(fHPionSm);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::UserExec(Option_t *)
{
  // Called for each event.
  
  if (!InputEvent())
    return;
  
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  
  
	AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (am->GetInputEventHandler());
	//if (!esdH)
	//{
	//	Printf("ERROR: Could not get ESDInputHandler");
	//}
  
	AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (am->GetInputEventHandler());
	if (!aodH && !esdH)
	{
		Printf("ERROR: Could not get AODInputHandler");
	}
  
	if(esdH)
    fEsdEv = dynamic_cast<AliESDEvent*> (esdH->GetEvent());
  
  else if(aodH) fAodEv = aodH->GetEvent();
  
  else{
    AliFatal("Neither ESD nor AOD event found");
    return;
  }
  
  // set all counters to nought
  ipymin=0;
  ipymax=0;
  ipi0min=0;
  ipi0max=0;
  ietamin=0;
  ietamax=0;
  
  if(fMcMode && fAddedSignal){
    // monte carlo headers from Evi
    AliMCEvent *mcEvent = MCEvent();
    if (!mcEvent){
      cout << "no MC event" << endl;
      return;
    }
    eventHeader = mcEvent->GenEventHeader();
    
    AliGenCocktailEventHeader *cocktail = dynamic_cast<AliGenCocktailEventHeader *>(eventHeader);
    
    if(!cocktail) return ;
    TList *genHeaders = cocktail->GetHeaders();
    
    Int_t nGenerators = genHeaders->GetEntries();
    Int_t pythiaLastP=0;
    //Int_t IsPi0Flat=0;
    //printf("N generators %d \n", nGenerators);
    //Int_t fNMCProducedMin = 0;
    //Int_t fNMCProducedMax = 0;
    
    // it seems gen1 and gen2 are pi0 and eta added signals, respectively, when looking at lhc12i3
    for(Int_t igen = 0; igen < nGenerators; igen++)
    {
      AliGenEventHeader* eventHeader2 = (AliGenEventHeader*)genHeaders->At(igen) ;
      TString name = eventHeader2->GetName();
      
      if (name.Contains("Pythia",TString::kIgnoreCase))
      {
        pythiaLastP=eventHeader2->NProduced();
        ipymax=eventHeader2->NProduced()-1;
        //Printf("pythia partcles = %d", pythiaLastP);
        pythiaHeader = (AliGenEventHeader*)genHeaders->At(igen);
      }
      //fNMCProducedMin = pythiaLastP;
      //fNMCProducedMax= fNMCProducedMin+eventHeader2->NProduced();
      
      //Printf("Generator %d: Class Name %s, Name %s, title %s \n, nProduced %d",igen, eventHeader2->ClassName(), name.Data(), eventHeader2->GetTitle(),fNMCProducedMax-fNMCProducedMin);
      
      //    if (name.Contains("Pi0Flat",TString::kIgnoreCase)) {IsPi0Flat=1; break;}
      if (igen==1) {
        addedPi0Header = (AliGenEventHeader*)genHeaders->At(igen);
        ipi0min=ipymax+1;
        ipi0max=ipi0min+addedPi0Header->NProduced()-1;
      }
      if (igen==2) {
        addedEtaHeader = (AliGenEventHeader*)genHeaders->At(igen);
        ietamin=ipi0max+1;
        ietamax=ietamin+addedEtaHeader->NProduced()-1;
        break;
      }
    }
    
    //Printf("min = %d, max = %d",fNMCProducedMin, fNMCProducedMax);
  }
  
  
  UInt_t offtrigger = 0;
  if (fEsdEv) {
    am->LoadBranch("AliESDRun.");
    am->LoadBranch("AliESDHeader.");
    UInt_t mask1 = fEsdEv->GetESDRun()->GetDetectorsInDAQ();
    UInt_t mask2 = fEsdEv->GetESDRun()->GetDetectorsInReco();
    Bool_t desc1 = (mask1 >> 18) & 0x1;
    Bool_t desc2 = (mask2 >> 18) & 0x1;
    if (desc1==0 || desc2==0) { //AliDAQ::OfflineModuleName(18)=="EMCAL"
      AliError(Form("EMCAL not in DAQ/RECO: %u (%u)/%u (%u)",
                    mask1, fEsdEv->GetESDRun()->GetDetectorsInReco(),
                    mask2, fEsdEv->GetESDRun()->GetDetectorsInDAQ()));
      return;
    }
  }
  
  if(fAodEv){
    am->LoadBranch("header");
    offtrigger =  ((AliVAODHeader*)fAodEv->GetHeader())->GetOfflineTrigger();
  }
  /*
   fEsdEv = dynamic_cast<AliESDEvent*>(InputEvent());
   UInt_t offtrigger = 0;
   if (fEsdEv) {
   am->LoadBranch("AliESDRun.");
   am->LoadBranch("AliESDHeader.");
   UInt_t mask1 = fEsdEv->GetESDRun()->GetDetectorsInDAQ();
   UInt_t mask2 = fEsdEv->GetESDRun()->GetDetectorsInReco();
   Bool_t desc1 = (mask1 >> 18) & 0x1;
   Bool_t desc2 = (mask2 >> 18) & 0x1;
   if (desc1==0 || desc2==0) { //AliDAQ::OfflineModuleName(18)=="EMCAL"
   AliError(Form("EMCAL not in DAQ/RECO: %u (%u)/%u (%u)",
   mask1, fEsdEv->GetESDRun()->GetDetectorsInReco(),
   mask2, fEsdEv->GetESDRun()->GetDetectorsInDAQ()));
   return;
   }
   offtrigger = ((AliInputEventHandler*)(am->GetInputEventHandler()))->IsEventSelected();
   }
   else {
   fAodEv = dynamic_cast<AliAODEvent*>(InputEvent());
   if (!fAodEv) {
   AliFatal("Neither ESD nor AOD event found");
   return;
   }
   am->LoadBranch("header");
   offtrigger =  fAodEv->GetHeader()->GetOfflineTrigger();
   }
   */
  if (!fMcMode && (offtrigger & AliVEvent::kFastOnly)) {
    AliWarning(Form("EMCAL not in fast only partition"));
    return;
  }
  
  
  if (fDoTrMatGeom && !AliGeomManager::GetGeometry()) { // get geometry
    AliWarning("Accessing geometry from OCDB, this is not very efficient!");
    AliCDBManager *cdb = AliCDBManager::Instance();
    if (!cdb->IsDefaultStorageSet())
      cdb->SetDefaultStorage("raw://");
    Int_t runno = InputEvent()->GetRunNumber();
    if (runno != cdb->GetRun())
      cdb->SetRun(runno);
    AliGeomManager::LoadGeometry();
  }
  
  if (!AliGeomManager::GetGeometry()&&!fIsGeoMatsSet) { // set misalignment matrices (stored in first event)
    Int_t nsm = fGeom->GetEMCGeometry()->GetNumberOfSuperModules();
    for (Int_t i=0; i<nsm; ++i) {
      const TGeoHMatrix *geom = 0;
      if (fEsdEv)
        geom = fEsdEv->GetESDRun()->GetEMCALMatrix(i);
      else{
        AliAODHeader * aodheader = dynamic_cast<AliAODHeader*>(fAodEv->GetHeader());
        if(!aodheader) AliFatal("Not a standard AOD");
        geom = aodheader->GetEMCALMatrix(i);
      }
      if (!geom)
        continue;
      geom->Print();
      fGeom->SetMisalMatrix(geom,i);
    }
    fIsGeoMatsSet = kTRUE;
  }
  
  if (!TGeoGlobalMagField::Instance()->GetField()) { // construct field map
    if (fEsdEv) {
      const AliESDRun *erun = fEsdEv->GetESDRun();
      erun->InitMagneticField();
    }
    else {
      Double_t pol = -1; //polarity
      Double_t be = -1;  //beam energy
      AliMagF::BeamType_t btype = AliMagF::kBeamTypepp;
      Int_t runno = fAodEv->GetRunNumber();
      if (runno>=136851 && runno<138275) {
        pol = -1;
        be = 2760;
        btype = AliMagF::kBeamTypeAA;
      } else if (runno>=138275 && runno<=139517) {
        pol = +1;
        be = 2760;
        btype = AliMagF::kBeamTypeAA;
      } else {
        AliError(Form("Do not know the bfield parameters for run %d! Using defaults!!!", runno));
      }
      TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", pol, pol, AliMagF::k5kG, btype, be));
    }
  }
  
  Int_t cut = 1;
  fHCuts->Fill(cut++);
  
  TString trgclasses;
  AliESDHeader *h = dynamic_cast<AliESDHeader*>(InputEvent()->GetHeader());
  if (h) {
    trgclasses = fEsdEv->GetFiredTriggerClasses();
  } else {
    AliAODHeader *h2 = dynamic_cast<AliAODHeader*>(InputEvent()->GetHeader());
    if (h2)
      trgclasses = h2->GetFiredTriggerClasses();
  }
  for (Int_t i = 0; i<fTrClassNamesArr->GetEntries(); ++i) {
    const char *name = fTrClassNamesArr->At(i)->GetName();
    TRegexp regexp(name);
    if (trgclasses.Contains(regexp))
      fHTclsBeforeCuts->Fill(1+i);
  }

//  
//  if(trgclasses.Contains("EMC") && trgclasses.Contains("INT")){
//    cout << "####################################################################### EMC AND INT SET #######################################################################" << endl;
//    cout << "trigger mask EMC: "<<  std::bitset<32>(fEsdEv->GetTriggerMask()) << " = " << fEsdEv->GetTriggerMask()<< endl;
//    cout << trgclasses << endl;
//  }
//
//  if(1){
//  if(trgclasses.Contains("EMC")){
//    cout << trgclasses << endl;
//    cout << "trigger mask EMC: "<<  std::bitset<32>(fEsdEv->GetTriggerMask()) << " = " << fEsdEv->GetTriggerMask()<< endl;
//    cout << "trigger bit 0 is " <<  std::bitset<32>(fEsdEv->GetTriggerMask()).test(0) << endl;
//  }
//  }
//  if(1){
//  if(trgclasses.Contains("INT")){
//    cout << trgclasses << endl;
//    cout << "trigger mask INT: "<<  std::bitset<32>(fEsdEv->GetTriggerMask()) << " = " << fEsdEv->GetTriggerMask()<< endl;
//    cout << "trigger bit 0 is " <<  std::bitset<32>(fEsdEv->GetTriggerMask()).test(0) << endl;
//  }
//  }

  if (fDoPSel && offtrigger==0)
    return;
  
  // cut on certain events
  if(fEsdEv){
    const AliESDVertex* vtxESD = fEsdEv->GetPrimaryVertexTracks();
    if( vtxESD->GetNContributors() < 1 )
      {
        return;
      }
	fHCuts->Fill(cut++);
	  
    TString trigClasses = fEsdEv->GetFiredTriggerClasses();
    // remove "fast cluster events": 
    if (trigClasses.Contains("FAST")  && !trigClasses.Contains("ALL"))
      return;
    if(!(fEsdEv->GetPrimaryVertex()->GetStatus()))   
      return;
  }
  else if(fAodEv){
    TString trigClasses = fAodEv->GetFiredTriggerClasses();
    // remove "fast cluster events": 
    if (trigClasses.Contains("FAST")  && !trigClasses.Contains("ALL"))
      return;
  }
  fHCuts->Fill(cut++);

  // in lhc11a, cut on even more events
  Int_t runnumber = InputEvent()->GetRunNumber();
  if ((runnumber>=144871) && (runnumber<=146860)) {
    
    AliVCaloCells *cells   = InputEvent()->GetEMCALCells();
    const Short_t nCells   = cells->GetNumberOfCells();
    
    if (InputEvent()->IsA()==AliESDEvent::Class()) AliAnalysisManager::GetAnalysisManager()->LoadBranch("EMCALCells.");
    
    AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!fInputHandler) return;
    
    // count cells above threshold
    Int_t nCellCount[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    for(Int_t iCell=0; iCell<nCells; ++iCell) {
      Short_t cellId = cells->GetCellNumber(iCell);
      Double_t cellE = cells->GetCellAmplitude(cellId);
      Int_t sm       = cellId / (24*48);
      if (cellE>0.1) ++nCellCount[sm];
    }
    
    Bool_t fIsLedEvent = kFALSE;
    if (nCellCount[4] > 100) {
      fIsLedEvent = kTRUE;
    } else {
      if ((runnumber>=146858) && (runnumber<=146860)) {
	if ((fInputHandler->IsEventSelected() & AliVEvent::kMB) && (nCellCount[3]>=21))
	  fIsLedEvent = kTRUE;
	else if ((fInputHandler->IsEventSelected() & AliVEvent::kEMC1) && (nCellCount[3]>=35))
	  fIsLedEvent = kTRUE;
      }
    }
    if (fIsLedEvent) {
      return;
    }
  }


  fHCuts->Fill(cut++);
  
  
  const AliCentrality *centP = InputEvent()->GetCentrality();
  Double_t cent = centP->GetCentralityPercentileUnchecked(fCentVar);
  fHCent->Fill(cent);
  if (cent<fCentFrom||cent>fCentTo)
    return;
  
  fHCuts->Fill(cut++);
  
  if (fUseQualFlag) {
    if (centP->GetQuality()>0)
      return;
  }
  
  fHCentQual->Fill(cent);
  fHCuts->Fill(cut++);
  
  if (fEsdEv) {
    am->LoadBranch("PrimaryVertex.");
    am->LoadBranch("SPDVertex.");
    am->LoadBranch("TPCVertex.");
  } else {
    fAodEv = dynamic_cast<AliAODEvent*>(InputEvent());
    am->LoadBranch("vertices");
    if (!fAodEv) return;
  }
  
  const AliVVertex *vertex = InputEvent()->GetPrimaryVertex();
  if (!vertex)
    return;
  
  fHVertexZ->Fill(vertex->GetZ());
  
  if(vertex->GetZ()<fVtxZMin||vertex->GetZ()>fVtxZMax)
    return;
  
  fHCuts->Fill(cut++);
  fHVertexZ2->Fill(vertex->GetZ());
  
  // count number of accepted events
  ++fNEvs;
  
  for (Int_t i = 0; i<fTrClassNamesArr->GetEntries(); ++i) {
    const char *name = fTrClassNamesArr->At(i)->GetName();
    TRegexp regexp(name);
    if (trgclasses.Contains(regexp))
      fHTclsAfterCuts->Fill(1+i);
  }
  
  fRecPoints   = 0; // will be set if fClusName is given and AliAnalysisTaskEMCALClusterizeFast is used
  fDigits      = 0; // will be set if fClusName is given and AliAnalysisTaskEMCALClusterizeFast is used
  fEsdClusters = 0; // will be set if ESD input used and if fRecPoints are not set or if clusters are attached
  fEsdCells    = 0; // will be set if ESD input used
  fAodClusters = 0; // will be set if AOD input used and if fRecPoints are not set or if clusters are attached
  //             or if fClusName is given and AliAnalysisTaskEMCALClusterizeFast in AOD output mode
  fAodCells    = 0; // will be set if AOD input used
  
  // deal with special output from AliAnalysisTaskEMCALClusterizeFast first
  Bool_t overwrite    = 0;
  Bool_t clusattached = 0;
  Bool_t recalibrated = 0;
  if (0 && !fClusName.IsNull()) {
    AliAnalysisTaskEMCALClusterizeFast *cltask = 0;
    TObjArray *ts = am->GetTasks();
    cltask = dynamic_cast<AliAnalysisTaskEMCALClusterizeFast*>(ts->FindObject(fClusName));
    if (cltask && cltask->GetClusters()) {
      fRecPoints   = cltask->GetClusters();
      fDigits      = cltask->GetDigits();
      clusattached = cltask->GetAttachClusters();
      overwrite    = cltask->GetOverwrite();
      if (cltask->GetCalibData()!=0)
        recalibrated = kTRUE;
    }
  }
  if (1 && !fClusName.IsNull()) {
    TList *l = 0;
    if (AODEvent())
      l = AODEvent()->GetList();
    else if (fAodEv)
      l = fAodEv->GetList();
    if (l) {
      fAodClusters = dynamic_cast<TClonesArray*>(l->FindObject(fClusName));
    }
  }
  
  if (fEsdEv) { // ESD input mode
    if (1 && (!fRecPoints||clusattached)) {
      if (!clusattached && !overwrite)
        am->LoadBranch("CaloClusters");
      TList *l = fEsdEv->GetList();
      if (clusattached) {
        fEsdClusters = dynamic_cast<TClonesArray*>(l->FindObject(fClusName));
      } else {
        fEsdClusters = dynamic_cast<TClonesArray*>(l->FindObject("CaloClusters"));
      }
    }
    if (1) {
      if (!recalibrated)
        am->LoadBranch("EMCALCells.");
      fEsdCells = fEsdEv->GetEMCALCells();
    }
  } else if (fAodEv) { // AOD input mode
    if (1 && (!fAodClusters || clusattached)) {
      if (!clusattached)
        am->LoadBranch("caloClusters");
      TList *l = fAodEv->GetList();
      if (l) {
        fAodClusters = dynamic_cast<TClonesArray*>(l->FindObject("caloClusters"));
      }
    }
    if (1) {
      if (!recalibrated)
        am->LoadBranch("emcalCells");
      fAodCells = fAodEv->GetEMCALCells();
    }
  } else {
    AliFatal("Impossible to not have either pointer to ESD or AOD event");
  }
  
  // need to get converted photons from something, maybe apply quality cuts, and then use them for invariant mass calculation
  // also need space in event mixing to store all those candidates
  //fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
  //if(!fV0Reader){printf("Error: No V0 Reader");} // GetV0Reader
  //Int_t eventQuality = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetEventQuality();
  
  //if(fV0Reader)
    //fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut // GetGoodGammas instead?
  //cout << "number of conversion photons: " << fV0Reader->GetNReconstructedGammas() << endl;
  
  //else
  fReaderGammas = NULL;
  
  if (1) {
    AliDebug(2,Form("fRecPoints   set: %p", fRecPoints));
    AliDebug(2,Form("fDigits      set: %p", fDigits));
    AliDebug(2,Form("fEsdClusters set: %p", fEsdClusters));
    AliDebug(2,Form("fEsdCells    set: %p", fEsdCells));
    AliDebug(2,Form("fAodClusters set: %p", fAodClusters));
    AliDebug(2,Form("fAodCells    set: %p", fAodCells));
  }
  
  if (fDoAfterburner)
    ClusterAfterburner();
  
  if (fMcMode)
  	CalcMcInfo();
  
  // mixed events
  
  // some info
  Int_t vtxClass = 1;
  Double_t vtxcuts = (fVtxZMax-fVtxZMin)/nZClass;
  if(vertex->GetZ()>fVtxZMin && vertex->GetZ()<=fVtxZMin+vtxcuts) vtxClass=0;
  else if(vertex->GetZ()>fVtxZMin+vtxcuts && vertex->GetZ()<=fVtxZMin+2*vtxcuts) vtxClass=1;
  else if(vertex->GetZ()>fVtxZMin+2*vtxcuts && vertex->GetZ()<=fVtxZMin+3*vtxcuts) vtxClass=2;
  else if(vertex->GetZ()>fVtxZMin+3*vtxcuts && vertex->GetZ()<=fVtxZMin+4*vtxcuts) vtxClass=3;
//  else if(vertex->GetZ()>fVtxZMin+4*vtxcuts && vertex->GetZ()<fVtxZMax) vtxClass=4;
  
  //    vtxClass = 0;
  
  Int_t MulClass = 0;
  
  GetMulClass(MulClass);
  
	Float_t phitrig = 0;
  Float_t thetatrig = 0;
  Double_t pt_max = 0;
  Int_t ptClass = 0;
  if (!fTrainMode) {
    pt_max = FillClusHists(phitrig, thetatrig);
    
    if(pt_max > 1 && pt_max<=3) ptClass = 1;
    else if(pt_max > 3 && pt_max<=6) ptClass = 2;
    else ptClass = 3;

    ptClass = 0;
    
    FillOtherHists();
  }
  FillMcHists();
  if(fDoNtuple)
    FillNtuple();
  
  if (fTrainMode) {
    fSelTracks->Clear();
    fSelPrimTracks->Clear();
    if (fMcParts)
      fMcParts->Clear();
    if (fTriggers)
      fTriggers->Clear();
    if (fClusters)
      fClusters->Clear();
  }
  
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::Terminate(Option_t *)
{
  // Terminate called at the end of analysis.
  
  if (fNtuple) {
    TFile *f = OpenFile(1);
    TDirectory::TContext context(f);
    if (f)
      fNtuple->Write();
  }
  
  AliInfo(Form("%s: Accepted %lld events          ", GetName(), fNEvs));
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::ClusterAfterburner()
{
  // Run custer reconstruction afterburner.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  
  if (!cells)
    return;
  
  Int_t ncells = cells->GetNumberOfCells();
  if (ncells<=0)
    return;
  
  Double_t cellMeanE = 0, cellSigE = 0;
  for (Int_t i = 0; i<ncells; ++i) {
    Double_t cellE = cells->GetAmplitude(i);
    cellMeanE += cellE;
    cellSigE += cellE*cellE;
  }
  cellMeanE /= ncells;
  cellSigE /= ncells;
  cellSigE -= cellMeanE*cellMeanE;
  if (cellSigE<0)
    cellSigE = 0;
  cellSigE = TMath::Sqrt(cellSigE / ncells);
  
  Double_t subE = cellMeanE - 7*cellSigE;
  if (subE<0)
    return;
  
  for (Short_t i = 0; i<ncells; ++i) {
    Short_t id=-1;
    Int_t mclabel = -1;
    Double_t amp=0,time=0, efrac = 0;
    if (!cells->GetCell(i, id, amp, time, mclabel, efrac))
      continue;
    amp -= cellMeanE;
    if (amp<0.001)
      amp = 0;
    cells->SetCell(i, id, amp, time);
  }
  
  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  if (!clusters)
    return;
  
  Int_t nclus = clusters->GetEntries();
  for (Int_t i = 0; i<nclus; ++i) {
    AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(i));
    if (!clus->IsEMCAL())
      continue;
    Int_t nc = clus->GetNCells();
    Double_t clusE = 0;
    UShort_t ids[100] = {0};
    Double_t fra[100] = {0};
    for (Int_t j = 0; j<nc; ++j) {
      Short_t id = TMath::Abs(clus->GetCellAbsId(j));
      Double_t cen = cells->GetCellAmplitude(id);
      clusE += cen;
      if (cen>0) {
        ids[nc] = id;
        ++nc;
      }
    }
    if (clusE<=0) {
      clusters->RemoveAt(i);
      continue;
    }
    
    for (Int_t j = 0; j<nc; ++j) {
      Short_t id = ids[j];
      Double_t cen = cells->GetCellAmplitude(id);
      fra[j] = cen/clusE;
    }
    clus->SetE(clusE);
    AliAODCaloCluster *aodclus = dynamic_cast<AliAODCaloCluster*>(clus);
    if (aodclus) {
      aodclus->Clear("");
      aodclus->SetNCells(nc);
      aodclus->SetCellsAmplitudeFraction(fra);
      aodclus->SetCellsAbsId(ids);
    }
  }
  clusters->Compress();
}


//________________________________________________________________________
Double_t AliAnalysisTaskEMCALDirGamma::FillClusHists(Float_t& max_phi, Float_t& max_theta)
{
   
    
  //cout << "filling cluster histograms ";
  
  Float_t momentum;
  Int_t PDG;
  Int_t tracksmatched;
  Int_t clusternumber = 0;
  Int_t energycut = 0;
  Int_t cellcut = 0;	
  Int_t bordercut = 0;	
  Int_t m02cut = 0;	
  Int_t trackcut = 0;		
  Float_t pos[3];

  bool bprint = 0;
  
  max_phi = 0;
  max_theta = 0;
	//Double_t max_pt = 0;
	// Fill histograms related to cluster properties.
  
  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  if (!clusters){
    cout << "no clusters node in event!" << endl;
    return 0;
  }
  
  Int_t nclus = clusters->GetEntries();
  
  //       cout << nclus << " clusters in event ";
  
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  fHClustNoEvt->Fill(nclus);
  
  int nclusters = 0;
  
  for(Int_t i = 0; i<nclus; ++i) {
    AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(i));
    if (!clus){continue;}
    // emcal cluster?
	//if (!clus->IsEMCAL()){continue;}
	
	
	Double_t maxAxis    = clus->GetTOF(); //sigma
	Double_t clusterEcc = clus->Chi2();   //eccentricity
	Double_t clusterDisp = clus->GetDispersion();
	fHClustEccentricity->Fill(clusterDisp);
	
	
	
    // here we only fill after fulfilling cuts)
    /*
     if(clusterVecCorr.Pt()>max_pt){
     max_phi = (float)clusterVecCorr.Phi();
     max_theta = (float)clusterVecCorr.Theta();
     max_pt = clusterVecCorr.Pt();
     }
     */
    // fill clusters into this event
    



    TLorentzVector clusterVec;
    clus->GetMomentum(clusterVec,vertex);

    if(bprint)
      clusterVec.Print();
    //  if(bDirGam){
      
    //}
	
	fHClustEtaPhiIsEMCal->Fill(clusterVec.Eta(),clusterVec.Phi());
    if(!clus->IsEMCAL()){continue;}






    // fill cluster histograms
    fHClustEtaPhiRaw->Fill(clusterVec.Eta(),clusterVec.Phi());
    fHClustEnergyPt->Fill(clusterVec.E(),clusterVec.Pt());
    fHClustEnergySM->Fill(clusterVec.E(),GetModuleNumber(clus));
    fHClustEnergySigma->Fill(clus->E()*maxAxis,clus->E());
    fHClustSigmaSigma->Fill(max(clus->GetM02(),clus->GetM20()),clus->E()*maxAxis);
    fHClustNCellEnergyRatio->Fill(clus->GetNCells(),GetMaxCellEnergy(clus)/clus->E());
    fHClustEnergyNCell->Fill(clus->E(),clus->GetNCells());
	fHPtSpecEffCluster->Fill(clusterVec.Pt());
	fHClustEtaM02->Fill(clusterVec.Eta(),clus->GetM02());
	fHClustPhiM02->Fill(clusterVec.Phi(),clus->GetM02());	
    nclusters++;

    if(!fMcMode){

      Double_t En = clus->E();
      // Jasons recalibration
      if(fDoManualRecal){
        Int_t fRecalibrator = 8;
        Double_t recalScale = PrivateEnergyRecal(clus->E(), fRecalibrator);
        En = clus->E()*recalScale;// TOTAL HACK - JJ

        clusterVec.SetPx(clusterVec.Px()*recalScale);// TOTAL HACK - JJ
        clusterVec.SetPy(clusterVec.Py()*recalScale);// TOTAL HACK - JJ
        clusterVec.SetPz(clusterVec.Pz()*recalScale);// TOTAL HACK - JJ
        clusterVec.SetE(clusterVec.E()*recalScale);// TOTAL HACK - JJ
        //      Double_t ecorr = 1;
        //Double_t ecorr = fcorrect->Eval(2.0*clusterVec1.E());
      }
      
      TLorentzVector clusterVecCorr1(clusterVec.Px(),clusterVec.Py(),clusterVec.Pz(),En);
	  
  	//BEFORE CUT
  	fHClustEnergyM02AllRaw->Fill(clus->E(),clus->GetM02());	
  	fHPtSpecAllRaw->Fill(clusterVec.Pt());




      // apply cluster cuts first
      if(clus->E()<fMinE){continue;}
      if(clus->GetNCells()<fNminCells){continue;}
      //if(GetMaxCellEnergy(clus)/clus->E()<fMinErat){continue;}
      //if(clus->Chi2()<fMinEcc){continue;}// eccentricity cut

  	//M02 CUT
      if(clus->GetM02()<fM02min){continue;}
    	if(clus->GetM02()<fM02min){continue;}	
  	if(clus->E() > 6 && clus->GetM02() > 0.40){continue;}		
  	if(clus->E() > 4 && clus->E() < 6 && clus->GetM02() > 0.365){continue;}		
  	if(clus->E() > 3 && clus->E() < 4 && clus->GetM02() > 0.365){continue;}	
  	if(clus->E() > 2 && clus->E() < 3 && clus->GetM02() > 0.385){continue;}	
  	if(clus->E() > 1.5 && clus->E() < 2 && clus->GetM02() > 0.32){continue;}	
  	if(clus->E() > 1 && clus->E() < 1.5 && clus->GetM02() > 0.32){continue;}		
  	if(clus->E() > 0.7 && clus->E() < 1 && clus->GetM02() > 0.35){continue;}		
  	if(clus->E() < 0.7 && clus->GetM02() > 0.35){continue;}		
 

  	//cut on 1 border bins of the EMCal
  	if(clusterVec.Phi() < 1.4215 || clusterVec.Phi() > 3.1145){continue;}	
  	if(clusterVec.Phi() > 1.7185 && clusterVec.Phi() < 1.7715){continue;}
  	if(clusterVec.Phi() > 2.0685 && clusterVec.Phi() < 2.1215){continue;}
  	if(clusterVec.Phi() > 2.4185 && clusterVec.Phi() < 2.4695){continue;}
  	if(clusterVec.Phi() > 2.7665 && clusterVec.Phi() < 2.8175){continue;}
  	if(clusterVec.Eta() < -0.653 || clusterVec.Eta() > 0.63691){continue;}
  	//if(clusterVec.Eta() > -0.01498 && clusterVec.Eta() < 0.01276){continue;}

  	//cut on 2 border bins of the EMCal
  	//if(clusterVec.Phi() < 1.435 || clusterVec.Phi() > 3.101){continue;}	
  	//if(clusterVec.Phi() > 1.705 && clusterVec.Phi() < 1.785){continue;}
  	//if(clusterVec.Phi() > 2.055 && clusterVec.Phi() < 2.135){continue;}
  	//if(clusterVec.Phi() > 2.405 && clusterVec.Phi() < 2.483){continue;}
  	//if(clusterVec.Phi() > 2.753 && clusterVec.Phi() < 2.831){continue;}
  	//if(clusterVec.Eta() < -0.63913 || clusterVec.Eta() > 0.63691){continue;}
  	//if(clusterVec.Eta() > -0.02885 && clusterVec.Eta() < 0.02663){continue;}

  	//TrackMatching
    if(clus->GetNTracksMatched() > 0){continue;}



  	//AFTER CUT
  	fHClustEnergyM02All->Fill(clus->E(),clus->GetM02());	
  	fHPtSpecAll->Fill(clusterVec.Pt());
  	fHClustEtaPhi->Fill(clusterVec.Eta(),clusterVec.Phi());
	  
	  
	  
	  
	  
	  
    }
    // go through MC information of clusters
    else{
      // MC labels
      int ilabel = -1;
      // highest contribution - this one might follow to the first mother
      ilabel = clus->GetLabel();
      // all contributors
      Int_t* mcarr =  clus->GetLabels();
      // number of contributors
      Int_t nl = clus->GetNLabels();
      
      //cout << clus->GetLabel() << endl;
      Bool_t bcl = 1;
      if(nl > 1){
        bcl = 0;
      }
      
      if(ilabel != -1){
        
        
        // get MC event
        AliMCEvent *mcEvent = MCEvent();
        if (!mcEvent){
          cout << "no MC event" << endl;
          return 0;
        }
        mcEvent->PreReadAll();
        // get MC particle for first label
        AliMCParticle *mcP = static_cast<AliMCParticle*>(mcEvent->GetTrack(ilabel));
        // for the case it did not work
        if (!mcP)
          continue;

	    const AliVVertex *evtVtx = mcEvent->GetPrimaryVertex();
	    //if(!evtVtx){return;}


        Int_t nTracksMC  = mcEvent->GetNumberOfTracks();
        Int_t nPTracksMC = mcEvent->GetNumberOfPrimaries();


		clus->GetPosition(pos);
		TVector3 vpos(pos[0],pos[1],pos[2]);


		Double_t dR = 1;			
		//dR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX()) * (mcP->Xv()-evtVtx->GetX()) + (mcP->Yv()-evtVtx->GetY()) * (mcP->Yv()-evtVtx->GetY()));

		/*
		//Test Cluster
		Int_t nSupMod, nModule, nIphi, nIeta, iphi, ieta; //var for tower-mapping
		//AliVCaloCells &cells= *(esd->GetEMCALCells());
	    AliVCaloCells *cells = fEsdCells;
	    if (!cells){continue;}
		
		Int_t testPhoton;
		if(mcP->PdgCode() == 22 && clus->GetM02() < 0.7 && clus->GetM02() > 0.5 && clus->E() > 0.7 && clus->E() < 1){
			if(testPhoton == 1){continue;}
			
			for (Int_t ic = 0; ic<clus->GetNCells(); ++ic) {
				Int_t absID   = TMath::Abs(clus->GetCellAbsId(ic));
				Float_t cellE = cells->GetCellAmplitude(absID); 
				
				fGeom->GetCellIndex(absID, nSupMod, nModule, nIphi, nIeta);
				fGeom->GetCellPhiEtaIndexInSModule(nSupMod, nModule, nIphi, nIeta, iphi, ieta); 	
					
				if(nSupMod==0){
					fHEMCalModule0->Fill(ieta,iphi);				
					testPhoton = 1;
				}
			}
		}
		
		Int_t testPion;
		if(mcP->PdgCode() == 211 && clus->GetM02() < 0.3 && clus->E() > 1 && clus->E() < 1.5){
			if(testPion == 1){continue;}
			
			for (Int_t ic = 0; ic<clus->GetNCells(); ++ic) {
				Int_t absID   = TMath::Abs(clus->GetCellAbsId(ic));
				Float_t cellE = cells->GetCellAmplitude(absID); 
				
				fGeom->GetCellIndex(absID, nSupMod, nModule, nIphi, nIeta);
				fGeom->GetCellPhiEtaIndexInSModule(nSupMod, nModule, nIphi, nIeta, iphi, ieta); 	
					
				if(nSupMod==0){
					fHEMCalModule0->Fill(ieta,iphi); 
					testPion = 1;
				}
			
			
			}	
		}	
		*/
		
		
		
		fHClustEP->Fill(clus->E()/clusterVec.Pt());
		
		
		//Histogramme mit clustern vor CUT
		
		//All Cluster
		fHClustEnergyM02AllRaw->Fill(clus->E(),clus->GetM02());	
		fHPtSpecAllRaw->Fill(clusterVec.Pt());
		fHParticlesRaw->Fill(mcP->PdgCode());
		fHClustNcell->Fill(clus->GetNCells());
		if(clus->GetNCells() == 1){fHClustNcell1->Fill(clusterVec.Pt());}
		
			

		//Photons
		if(mcP->PdgCode() == 22){
			fHPtSpecGammaRaw->Fill(clusterVec.Pt());
	        fHClustEnergyM02GammaAllRaw->Fill(clus->E(),clus->GetM02());
			if(clus->GetNCells() == 1){fHClustEnergyM02GammaCell1->Fill(clus->E(),clus->GetM02());}
			fHClustM02M20GammaAllRaw->Fill(clus->GetM20(),clus->GetM02());	
			fHClustEnergyNCellPhotonRaw->Fill(clus->E(),clus->GetNCells());
			fHClustNCellM02PhotonRaw->Fill(clus->GetNCells(),clus->GetM02());
			fHClustNcellPhoton->Fill(clus->GetNCells());
			fHClustEnergyRatioPhoton->Fill(GetMaxCellEnergy(clus)/clus->E());
			if(clus->E() < 1){fHGammaMIP->Fill(clus->E());}
		}	
		//All except photons
		if(mcP->PdgCode() != 22){
			fHPtSpecNoGammaRaw->Fill(clusterVec.Pt());
			fHClustM02M20NoGammaRaw->Fill(clus->GetM20(),clus->GetM02());	
			fHClustNcellNoPhoton->Fill(clus->GetNCells());			
		}
		//charged Hadrons (Pions, (Anti-) Protons,)
		if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){
			fHPtSpecPionRaw->Fill(clusterVec.Pt());
	        fHClustEnergyM02PionRaw->Fill(clus->E(),clus->GetM02());
			if(clus->GetNCells() == 1){fHClustEnergyM02PionCell1->Fill(clus->E(),clus->GetM02());}				
			fHClustEnergyNCellPionRaw->Fill(clus->E(),clus->GetNCells());
			fHClustNCellM02PionRaw->Fill(clus->GetNCells(),clus->GetM02());
			fHClustEnergyRatioPion->Fill(GetMaxCellEnergy(clus)/clus->E());
			if(clus->E() < 1){fHHadronMIP->Fill(clus->E());}
		}
		//Electrons + Positrons
		if(mcP->PdgCode() == 11 || mcP->PdgCode() == -11){				
			fHPtSpecElectronRaw->Fill(clusterVec.Pt());
			fHClustEnergyM02ElektronRaw->Fill(clus->E(),clus->GetM02());
		}	
		//charged particles
		if(mcP->Charge() > 0 || mcP->Charge() < 0){
			fHPtSpecChargedRaw->Fill(clusterVec.Pt());	
		}			
		//Myons
		if(mcP->PdgCode() == 13 || mcP->PdgCode() == -13){				
			fHPtSpecMyonRaw->Fill(clusterVec.Pt());
		}	
		//Anti- (Neutrons)
		if(mcP->PdgCode() == 2112 || mcP->PdgCode() == -2112){
			fHPtSpecNeutronRaw->Fill(clusterVec.Pt());
		}
		////Anti- (Protons)
		//if(mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){
		//	fHPtSpecProtonRaw->Fill(clusterVec.Pt());
		//}
		//charged Kaons
		if(mcP->PdgCode() == 321 || mcP->PdgCode() == -321 || 
		   mcP->PdgCode() == 323 || mcP->PdgCode() == -323 ){
			fHPtSpecKaonRaw->Fill(clusterVec.Pt());
		}	
		//neutral Kaons
		if(mcP->PdgCode() == 130 || mcP->PdgCode() == -130 || 
		   mcP->PdgCode() == 310 || mcP->PdgCode() == -310 ||
		   mcP->PdgCode() == 311 || mcP->PdgCode() == -311 ||
		   mcP->PdgCode() == 313 || mcP->PdgCode() == -313){
			fHPtSpecKaon0Raw->Fill(clusterVec.Pt());
		}	

		





		
		//check where electrons come from (creationvertex)
		dR = 1;			
		if(mcP->PdgCode() == 11 || mcP->PdgCode() == -11){
		
			dR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX()) * (mcP->Xv()-evtVtx->GetX()) + (mcP->Yv()-evtVtx->GetY()) * (mcP->Yv()-evtVtx->GetY()));
			
			fHParticlesCompare->Fill(mcP->PdgCode());
			fHParticleR->Fill(dR);		
			//fHClustElectronZR->Fill(dR,pos[2]);

						//emother = mcP->GetMother();
			AliMCParticle *mcMe = static_cast<AliMCParticle*>(mcEvent->GetTrack(mcP->GetMother()));
		
			dR = TMath::Sqrt((mcMe->Xv()-evtVtx->GetX()) * (mcMe->Xv()-evtVtx->GetX()) + (mcMe->Yv()-evtVtx->GetY()) * (mcMe->Yv()-evtVtx->GetY()));
			
			fHParticlesCompare2->Fill(mcMe->PdgCode());
			//fHMotherR->Fill(dR);		
			
			
			if(mcMe->PdgCode() == 22){
				AliMCParticle *mcMe2 = static_cast<AliMCParticle*>(mcEvent->GetTrack(mcMe->GetMother()));
			
			
				dR = TMath::Sqrt((mcMe2->Xv()-evtVtx->GetX()) * (mcMe2->Xv()-evtVtx->GetX()) + (mcMe2->Yv()-evtVtx->GetY()) * (mcMe2->Yv()-evtVtx->GetY()));
			
			
				fHParticlesCompare3->Fill(mcMe2->PdgCode());	
				fHMotherR2->Fill(dR);		
			}
			
		}
		

		
		
		
		
		
		  
		//if(dR < 400){
		//	
		//	fHPtSpecCompare->Fill(mcPe->Pt());
		//	
		//		
		//	//if(mcP->PdgCode() == 2112 || mcP->PdgCode() == -2112){
		//	//if(mcP->PdgCode() == 11 || mcP->PdgCode() == -11){
		//		fHParticlesCompare->Fill(mcPe->PdgCode());	
		//		//}
		//	
		//	//if(mcP->PdgCode() == 11 || mcP->PdgCode() == -11){
		//  		if(mcPe->PdgCode() == 22){fHPtSpecGammaCompare->Fill(mcPe->Pt());}
		//		//}
		//	
		//  	if(mcPe->PdgCode() == 211 || mcPe->PdgCode() == -211){
		//		fHPtSpecPionCompare->Fill(mcPe->Pt());
		//	}
		//	//if(mcPe->PdgCode() == 22){
		//    //    fHClustEnergyM02GammaAllRaw->Fill(clus->E(),clus->GetM02());				
		//	//	fHPtSpecGammaAllRaw->Fill(clusterVec.Pt());
		//	//	fHClustEnergyNCellPhoton->Fill(clus->E(),clus->GetNCells());
		//	//}			
		//}		
		//
		////if(emother < 0 && dR > 4.0){continue;}
        //
        //emother = mcPe->GetMother();
		//if(emother >= 0){epart = emother;}			  
		//  	  
		//	  
		//
		////}//while-Schleife
	
	
		
		
		
		//Apply Cluster Cuts
   	 	cellcut = 0;	
   		bordercut = 0;	
   		m02cut = 0;	
   		trackcut = 0;
   		energycut = 0;	
    
    
   		if(clus->E()<fMinE){energycut = 1;}	
    	if(clus->GetNCells()<fNminCells){cellcut = 1;}
    	//if(GetMaxCellEnergy(clus)/clus->E()<fMinErat){continue;}
    	//if(clus->Chi2()<fMinEcc){continue;}// eccentricity cut
    	
		//M02 CUT
    	if(clus->GetM02()<fM02min){m02cut = 1;}	
		if(clus->E() > 6 && clus->GetM02() > 0.40){m02cut = 1;}		
		if(clus->E() > 4 && clus->E() < 6 && clus->GetM02() > 0.365){m02cut = 1;}		
		if(clus->E() > 3 && clus->E() < 4 && clus->GetM02() > 0.365){m02cut = 1;}	
		if(clus->E() > 2 && clus->E() < 3 && clus->GetM02() > 0.385){m02cut = 1;}	
		if(clus->E() > 1.5 && clus->E() < 2 && clus->GetM02() > 0.32){m02cut = 1;}	
		if(clus->E() > 1 && clus->E() < 1.5 && clus->GetM02() > 0.32){m02cut = 1;}		
		if(clus->E() > 0.7 && clus->E() < 1 && clus->GetM02() > 0.35){m02cut = 1;}		
		if(clus->E() < 0.7 && clus->GetM02() > 0.35){m02cut = 1;}		
			
		//cut on 1 border bins of the EMCal
		if(clusterVec.Phi() < 1.4215 || clusterVec.Phi() > 3.1145){bordercut = 1;}	
		if(clusterVec.Phi() > 1.7185 && clusterVec.Phi() < 1.7715){bordercut = 1;}
		if(clusterVec.Phi() > 2.0685 && clusterVec.Phi() < 2.1215){bordercut = 1;}
		if(clusterVec.Phi() > 2.4185 && clusterVec.Phi() < 2.4695){bordercut = 1;}
		if(clusterVec.Phi() > 2.7665 && clusterVec.Phi() < 2.8175){bordercut = 1;}
		if(clusterVec.Eta() < -0.653 || clusterVec.Eta() > 0.63691){bordercut = 1;}
		//if(clusterVec.Eta() > -0.01498 && clusterVec.Eta() < 0.01276){bordercut = 1;} //Zwischen benachbarter Supermodule keine Lücke 
    	
		//cut on 2 border bins of the EMCal
		//if(clusterVec.Phi() < 1.435 || clusterVec.Phi() > 3.101){bordercut = 1;}	
		//if(clusterVec.Phi() > 1.705 && clusterVec.Phi() < 1.785){bordercut = 1;}
		//if(clusterVec.Phi() > 2.055 && clusterVec.Phi() < 2.135){bordercut = 1;}
		//if(clusterVec.Phi() > 2.405 && clusterVec.Phi() < 2.483){bordercut = 1;}
		//if(clusterVec.Phi() > 2.753 && clusterVec.Phi() < 2.831){bordercut = 1;}
		//if(clusterVec.Eta() < -0.63913 || clusterVec.Eta() > 0.63691){bordercut = 1;}
		//if(clusterVec.Eta() > -0.02885 && clusterVec.Eta() < 0.02663){bordercut = 1;}
		
    	//TrackMatching
		if(cellcut == 0 && bordercut == 0 && m02cut == 0 && energycut == 0){
			if(mcP->Charge() < 1 && mcP->Charge() > -1){fHPtSpecChargedTM->Fill(clusterVec.Pt());}
			if(mcP->PdgCode() == 22){fHPtSpecGammaTM->Fill(clusterVec.Pt());}	
		}
		fHMotherR->Fill(clus->GetNTracksMatched());
		tracksmatched = clus->GetNTracksMatched();
		if(tracksmatched > 0){
			trackcut = 1;
			fHPtSpecEffTM1->Fill(clusterVec.Pt());
			if(mcP->Charge() < 1 && mcP->Charge() > -1){fHPtSpecEffTM2->Fill(clusterVec.Pt());}
			//fHClustEP->Fill(clus->E()/clusterVec.Pt());
			fHParticlesTM->Fill(mcP->PdgCode());	
			if(mcP->PdgCode() == 11 || mcP->PdgCode() == -11){fHPtSpecElectronTM->Fill(clusterVec.Pt());}
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){fHPtSpecPionTM->Fill(clusterVec.Pt());}
		}
		else {
			fHClustEnergyM02AllTM->Fill(clus->E(),clus->GetM02());
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){
				fHClustEnergyM02PionTM->Fill(clus->E(),clus->GetM02());
				fHPtSpecPionNoTM->Fill(clusterVec.Pt());
			}
			if(mcP->PdgCode() == 11 || mcP->PdgCode() == -11){
				dR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX()) * (mcP->Xv()-evtVtx->GetX()) + (mcP->Yv()-evtVtx->GetY()) * (mcP->Yv()-evtVtx->GetY()));
				fHParticleRcut->Fill(dR);
				fHPtSpecElectronNoTM->Fill(clusterVec.Pt());
			}
			
			
		}
	
	
		//m02 studies after cell, border and energycut
	if(cellcut == 0 && bordercut == 0 && energycut == 0){
		if(mcP->PdgCode() == 22){fHClustEnergyM02GammaSmallCut->Fill(clus->E(),clus->GetM02());}	
		if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){fHClustEnergyM02PionSmallCut->Fill(clus->E(),clus->GetM02());}		
			
		if(clus->E() < 0.7){
			if(mcP->PdgCode() == 22){fHClustM02Gamma0TM->Fill(clus->GetM02());}
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){fHClustM02Pion0TM->Fill(clus->GetM02());}
		}
		if(clus->E() > 0.7 && clus->E() < 1){
			if(mcP->PdgCode() == 22){fHClustM02Gamma1TM->Fill(clus->GetM02());}
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){fHClustM02Pion1TM->Fill(clus->GetM02());}
		}
		if(clus->E() > 1 && clus->E() < 1.5){
			if(mcP->PdgCode() == 22){fHClustM02Gamma2TM->Fill(clus->GetM02());}
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){fHClustM02Pion2TM->Fill(clus->GetM02());}
		}	
		if(clus->E() > 1.5 && clus->E() < 2){
			if(mcP->PdgCode() == 22){fHClustM02Gamma3TM->Fill(clus->GetM02());}
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){fHClustM02Pion3TM->Fill(clus->GetM02());}
		}			
		if(clus->E() > 2 && clus->E() < 3){
			if(mcP->PdgCode() == 22){fHClustM02Gamma4TM->Fill(clus->GetM02());}
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){fHClustM02Pion4TM->Fill(clus->GetM02());}
		}				
		if(clus->E() > 3 && clus->E() < 4){
			if(mcP->PdgCode() == 22){fHClustM02Gamma5TM->Fill(clus->GetM02());}
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){fHClustM02Pion5TM->Fill(clus->GetM02());}
		}		
		if(clus->E() > 4 && clus->E() < 6){
			if(mcP->PdgCode() == 22){fHClustM02Gamma6TM->Fill(clus->GetM02());}
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){fHClustM02Pion6TM->Fill(clus->GetM02());}
		}		
		if(clus->E() > 6){
			if(mcP->PdgCode() == 22){fHClustM02Gamma7TM->Fill(clus->GetM02());}
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){fHClustM02Pion7TM->Fill(clus->GetM02());}
		}	
		
		
	}		
	
	
	

	//crosscheck trackmatching
		if(cellcut == 0 && bordercut == 0 && m02cut == 0 && energycut == 0){
			if(mcP->PdgCode() == 11 || mcP->PdgCode() == -11){		
				dR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX()) * (mcP->Xv()-evtVtx->GetX()) + (mcP->Yv()-evtVtx->GetY()) * (mcP->Yv()-evtVtx->GetY()));
				fHParticleRcutcut->Fill(dR);
				
			}
			if(mcP->PdgCode() == 22){fHPtSpecGammaM02->Fill(clusterVec.Pt());}
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){fHPtSpecPionM02->Fill(clusterVec.Pt());}
		}
		
		if(cellcut == 0 && bordercut == 0 && energycut == 0 && trackcut == 0){
			fHParticlesNoTM->Fill(mcP->PdgCode());
			if(mcP->PdgCode() == 22){fHPtSpecGammaNoM02->Fill(clusterVec.Pt());}
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){fHPtSpecPionNoM02->Fill(clusterVec.Pt());}
		}	
		
		
		
		
		
		
		
		
		
		
		
		//Histogramme mit clustern nach ALLEN CUTS	
		if(cellcut == 0 && bordercut == 0 && m02cut == 0 && trackcut == 0 && energycut == 0){
	
			dR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX()) * (mcP->Xv()-evtVtx->GetX()) + (mcP->Yv()-evtVtx->GetY()) * (mcP->Yv()-evtVtx->GetY()));
	
			//All Cluster
			fHClustEnergyM02All->Fill(clus->E(),clus->GetM02());	
			fHPtSpecAll->Fill(clusterVec.Pt());
			fHParticles->Fill(mcP->PdgCode());	
			fHClustEtaPhi->Fill(clusterVec.Eta(),clusterVec.Phi());

			
			//Photons
			if(mcP->PdgCode() == 22){
				fHPtSpecGamma->Fill(clusterVec.Pt());
	    	    fHClustEnergyM02GammaAll->Fill(clus->E(),clus->GetM02());	
				fHClustEnergyNCellPhoton->Fill(clus->E(),clus->GetNCells());
				fHClustNCellM02Photon->Fill(clus->GetNCells(),clus->GetM02());	
				fHClustNcellPhotonCut->Fill(clus->GetNCells());			
			}	
			//All except photons
			if(mcP->PdgCode() != 22){
				fHPtSpecNoGamma->Fill(clusterVec.Pt());
				fHClustNcellNoPhotonCut->Fill(clus->GetNCells());	
			}
			if(mcP->Charge() > 0 || mcP->Charge() < 0){
				fHPtSpecCharged->Fill(clusterVec.Pt());	
			}			
			
			//charged Pions
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){
				fHPtSpecPion->Fill(clusterVec.Pt());
	    	    fHClustEnergyM02Pion->Fill(clus->E(),clus->GetM02());	
				fHClustEnergyNCellPion->Fill(clus->E(),clus->GetNCells());
				fHClustNCellM02Pion->Fill(clus->GetNCells(),clus->GetM02());			
			}
			//Electrons + Positrons
			if(mcP->PdgCode() == 11 || mcP->PdgCode() == -11){				
				fHPtSpecElectron->Fill(clusterVec.Pt());
				fHClustEnergyM02ElektronRaw->Fill(clus->E(),clus->GetNCells());
				dR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX()) * (mcP->Xv()-evtVtx->GetX()) + (mcP->Yv()-evtVtx->GetY()) * (mcP->Yv()-evtVtx->GetY()));
				fHParticleRcutcutcut->Fill(dR);
			}	
			//Myons
			if(mcP->PdgCode() == 13 || mcP->PdgCode() == -13){				
				fHPtSpecMyon->Fill(clusterVec.Pt());
			}	
			//Anti- (Neutrons)
			if(mcP->PdgCode() == 2112 || mcP->PdgCode() == -2112){
				fHPtSpecNeutron->Fill(clusterVec.Pt());
			}
			//Anti- (Protons)
			if(mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){
				fHPtSpecProton->Fill(clusterVec.Pt());
			}
	
			//charged Kaons
			if(mcP->PdgCode() == 321 || mcP->PdgCode() == -321 || 
			   mcP->PdgCode() == 323 || mcP->PdgCode() == -323 ){
				fHPtSpecKaon->Fill(clusterVec.Pt());
			}	
			//neutral Kaons
			if(mcP->PdgCode() == 130 || mcP->PdgCode() == -130 || 
			   mcP->PdgCode() == 310 || mcP->PdgCode() == -310 ||
			   mcP->PdgCode() == 311 || mcP->PdgCode() == -311 ||
			   mcP->PdgCode() == 313 || mcP->PdgCode() == -313){
				fHPtSpecKaon0->Fill(clusterVec.Pt());
			}
			
			//check Ncell-Cut 
			if(mcP->PdgCode() == 22){fHClustNcellPhotonCut->Fill(clus->GetNCells());}
			if(mcP->PdgCode() == 211 || mcP->PdgCode() == -211 || mcP->PdgCode() == 2212 || mcP->PdgCode() == -2212){fHClustNcellNoPhotonCut->Fill(clus->GetNCells());}
	
			
		}
    
		
		
	
		
	
	
	
	/*
	
        // sum of energy of all MC particles in cluster
        Double_t esum = 0;
        for(int ip=0;ip<nl;ip++){
          Int_t entry = mcarr[ip];
          AliMCParticle *mcPart = static_cast<AliMCParticle*>(mcEvent->GetTrack(entry));
          esum += mcPart->E();
        }
        Double_t efrac = 0;
        
        Double_t mce = mcP->E();
        Double_t cle = clus->E();
        
        // energy fraction of "leading" particle
        if(esum!=0)
          efrac = mcP->E()/esum;
        
        //        fHRecTrue->SetXTtitle("E_{clu}");
        //        fHRecTrue->SetYTtitle("E_{tru}/E_{clu}");
        //cout << "filling histos" << endl;
        // if generator
        if(bGen){
          // if leading particle with high contribution and photon
          if(mcP->PdgCode() == 22 && efrac>=0.99){
            // lets fill histograms
            fHMCpartfrac->Fill(mce/cle,efrac);
            fHECluEMC->Fill(mce,cle);
//            fHRecTrue->Fill(cle,mce/cle);
          }
          // if photon
          else if(mcP->PdgCode() == 22 && efrac<0.99){
            fHMCpartfrac->Fill(mce/cle,efrac);
            fHECluEMCnofull->Fill(mce,cle);
          }
          // if electron with high contribution
          else if((mcP->PdgCode() == 11 || mcP->PdgCode() == -11) && efrac>=0.99){
            fHECluEMCelectron->Fill(mce,cle);
          }
          // if pion with high contribution
          else if((mcP->PdgCode() == 211 || mcP->PdgCode() == -211) && efrac>=0.99){
            fHECluEMCpion->Fill(mce,cle);
          }
          // if kaon with high contribution
          else if((mcP->PdgCode() == 321 || mcP->PdgCode() == -321) && efrac>=0.99){
            fHECluEMCkaon->Fill(mce,cle);
          }
          // if other particle with high contribution
          else if(efrac>=0.99){
            fHECluEMCother->Fill(mce,cle);
          }
          
          // should I look for two photon clusters? ingredients: 2 photons, from same mother.
          for(int ip=0;ip<nl-1;ip++){
            Int_t entry = mcarr[ip];
            AliMCParticle *mcPart = static_cast<AliMCParticle*>(mcEvent->GetTrack(entry));
            imother = mcPart->GetMother();
            for(int jp=ip;jp<nl;jp++){
              Int_t entry2 = mcarr[jp];
              AliMCParticle *mcPart2 = static_cast<AliMCParticle*>(mcEvent->GetTrack(entry2));
              Int_t jmother = mcPart2->GetMother();
              if(imother == jmother && imother >= 0){
                AliMCParticle* mcMother = static_cast<AliMCParticle*>(mcEvent->GetTrack(imother));
                if(mcMother->PdgCode() == 111){
                  Double_t mcmop = sqrt(mcMother->P()*mcMother->P() + 0.135*0.135);
                  Double_t clue = clus->E();
                  fHECluEMCpi0single->Fill(mcmop,clue);
                }
              }
            }
          } // end two photon clusters
        } // end if generator
        if(bAddPi0){
          // if leading particle with high contribution and photon
          if(mcP->PdgCode() == 22 && efrac>0.5){
            // lets fill histograms
            fHECluEMCAddPi0->Fill(mce,cle);
            //fHRecTrueAddPi0->Fill(cle,mce/cle);
          }
        }
        if(bAddEta){
          // if leading particle with high contribution and photon
          if(mcP->PdgCode() == 22 && efrac>0.5){
            // lets fill histograms
            fHECluEMCAddEta->Fill(mce,cle);
            //fHRecTrueAddEta->Fill(cle,mce/cle);
          }
        }
        if(bAddPi0 || bAddEta){
          // if photon
          if(mcP->PdgCode() == 22 && efrac<0.99){
            fHECluEMCnofullAdd->Fill(mce,cle);
          }
        } // end else
        

        // Jasons recalibration
        Double_t En = clus->E();
        if(fDoManualRecal){
          Int_t fRecalibrator = 9;
          Double_t recalScale = PrivateEnergyRecal(clus->E(), fRecalibrator);
          En = clus->E()*recalScale;// TOTAL HACK - JJ
        
          //clus->GetMomentum(clusterVec,vertex);
          clusterVec.SetPx(clusterVec.Px()*recalScale);// TOTAL HACK - JJ
          clusterVec.SetPy(clusterVec.Py()*recalScale);// TOTAL HACK - JJ
          clusterVec.SetPz(clusterVec.Pz()*recalScale);// TOTAL HACK - JJ
          clusterVec.SetE(clusterVec.E()*recalScale);// TOTAL HACK - JJ
          //      Double_t ecorr = 1;
          //Double_t ecorr = fcorrect->Eval(2.0*clusterVec1.E());
        }
        TLorentzVector clusterVecCorr1(clusterVec.Px(),clusterVec.Py(),clusterVec.Pz(),En);
        
        Bool_t bkeep = 1;
        // for added signals, if particle is an electron, see if it is the higher energetic partner from conversion
//        if(!bGen && (bAddPi0 || bAddEta)){
//          
//          if(mcP->PdgCode() == 11 || mcP->PdgCode() == -11){
//            imother = mcP->GetMother();
//            AliMCParticle *tmppart2 = static_cast<AliMCParticle*>(mcEvent->GetTrack(imother));
//            // check if mother is a photon
//            if(tmppart2->PdgCode() == 22){
//              Int_t d1 = tmppart2->GetDaughterFirst();
//              Int_t d2 = tmppart2->GetDaughterLast();
//              if (d1>0){
//                if (d2<0){
//                  d2=d1;
//                }
//                if(d2-d1 == 1){
//                  for (Int_t ida=d1;ida<=d2;++ida) {
//                    AliMCParticle *tmpdaughter = static_cast<AliMCParticle *>(mcEvent->GetTrack(ida));
//                    if(tmpdaughter->E() > mcP->E())
//                      bkeep = 0;
//                  }
//                }
//              }
//            }
//          }
//        }
        
        // keep only eta -> gammagamma
        if(bAddEta){
          Int_t d1 = McMo->GetDaughterFirst();
          Int_t d2 = McMo->GetDaughterLast();
          if (d1>0){
            if (d2<0){
              d2=d1;
            }
            if(d2-d1 != 1){
              bkeep = 0;
            }
            else{
              AliMCParticle *tmpdaughter1 = static_cast<AliMCParticle *>(mcEvent->GetTrack(d1));
              AliMCParticle *tmpdaughter2 = static_cast<AliMCParticle *>(mcEvent->GetTrack(d2));
              if(!(tmpdaughter1->PdgCode() == 22 && tmpdaughter2->PdgCode() == 22)){
                bkeep = 0;
              }
            }
          }
          
        }
        if(1){
          if(bGen){
            // use hittype to flag particle history
            // if it is not from a pi0, flag 100
            // from primary pi0, flag 101
            // from secondary pi0, flag 102
            // from K0, flag 103
            // from material, flag 104
            
          }
          else{
            if(bAddPi0){
              thisEvent.hit[nclusters-1].hittype=1;
            }
            else if(bAddEta){
              thisEvent.hit[nclusters-1].hittype=2;
            }
          }
	
		}*/
      } // end ilabel
    } // end fmcmode
	
  } // end cluster loop
  fHClustAccEvt->Fill(nclusters);
  return 1;


}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::CalcMcInfo()
{
  // Get Mc truth particle information.
  if (!fMcMode)
    return;
  
  if (!fMcParts)
    return;
  
  fMcParts->Clear();
  
  AliEMCALEMCGeometry *emc = fGeom->GetEMCGeometry();
  Double_t etamin = emc->GetArm1EtaMin();
  Double_t etamax = emc->GetArm1EtaMax();
  Double_t phimin = emc->GetArm1PhiMin()*TMath::DegToRad();
  Double_t phimax = emc->GetArm1PhiMax()*TMath::DegToRad();
  
  //cout<<"eta: "<<etamin<<" -- "<<etamax<<endl;
  //cout<<"phi: "<<phimin<<" -- "<<phimax<<endl;  
  
  //  if (fAodEv) {
  //    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  //    am->LoadBranch(AliAODMCParticle::StdBranchName());
  //    TClonesArray *tca = dynamic_cast<TClonesArray*>(fAodEv->FindListObject(AliAODMCParticle::StdBranchName()));
  //    if (!tca)
  //      return;
  //
  //    Int_t nents = tca->GetEntries();
  //    for(int it=0; it<nents; ++it) {
  //      AliAODMCParticle *part = static_cast<AliAODMCParticle*>(tca->At(it));
  //      part->Print();
  //
  //      // pion or eta meson or direct photon
  //      if(part->GetPdgCode() == 111) {
  //      } else if(part->GetPdgCode() == 221) {
  //      } else if(part->GetPdgCode() == 22 ) {
  //      }	else
  //        continue;
  //
  //      // primary particle
  //      Double_t dR = TMath::Sqrt((part->Xv()*part->Xv())+(part->Yv()*part->Yv()));
  //      if(dR > 1.0)
  //        continue;
  //
  //      // kinematic cuts
  //      Double_t pt = part->Pt() ;
  //      if (pt<0.5)
  //        continue;
  //      Double_t eta = part->Eta();
  //      if (eta<etamin||eta>etamax)
  //        continue;
  //      Double_t phi  = part->Phi();
  //      if (phi<phimin||phi>phimax)
  //        continue;
  //
  //      ProcessDaughters(part, it, tca);
  //      PrintDaughters(part, tca, 2);
  //    }
  //    return;
  //  }
  
  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent){
    cout << "no MC event" << endl;
    return;
  }
  
  const AliVVertex *evtVtx = mcEvent->GetPrimaryVertex();
  if (!evtVtx)
    return;
  
  mcEvent->PreReadAll();
  
  
  
//Efficiency
//count the number of photons that fly in EMCAL direction	

	Bool_t bprim = kFALSE;
	Int_t noTracks = mcEvent->GetNumberOfTracks();
	for (Int_t jTrack = 0; jTrack<noTracks; ++jTrack){
	AliMCParticle *mcPa = static_cast<AliMCParticle*>(mcEvent->GetTrack(jTrack));
  	if(!mcPa){continue;}
	
	//if(mcPa->IsPrimary()){cout<<" ---- PRIMARY"<<endl;}
	
	
 		Double_t dR = TMath::Sqrt((mcPa->Xv()-evtVtx->GetX())*(mcPa->Xv()-evtVtx->GetX()) +
 	                  (mcPa->Yv()-evtVtx->GetY())*(mcPa->Yv()-evtVtx->GetY()));

 		if(dR > 0.1){continue;}
	fHParticlesEff->Fill(mcPa->PdgCode());		
	
  	//only photons
  	if(mcPa->PdgCode() != 22){continue;}
	
	
	
	
	//bprim = mcPa->IsPrimary();
	//if(bprim){		
  		fHPtSpecEffParticle->Fill(mcPa->Pt());			
		//}

		
		//If photon isn't a primary, check if mother is pi0, eta, rho, ...
		//if(!bprim){
			Int_t jmother = mcPa->GetMother();
		AliMCParticle* mcMo = static_cast<AliMCParticle*>(mcEvent->GetTrack(jmother));
	  		if (!mcMo){continue;}
		
  		if(mcMo->PdgCode() == 111) {	
			fHPtSpecDecayPi0->Fill(mcPa->Pt());						//Pi0
  		} else if(mcMo->PdgCode() == 221 ) {
			fHPtSpecDecayEta->Fill(mcPa->Pt());						//Eta
  		} else if(mcMo->PdgCode() == 223 ) {
			fHPtSpecDecayOmega->Fill(mcPa->Pt());					//Omega
  		} else if(mcMo->PdgCode() == 333 ) {	//Phi	  
  		} else if(mcMo->PdgCode() == 113 ) {	//Rho
  		} else if(mcMo->PdgCode() == 331 ) {
			fHPtSpecDecayEtap->Fill(mcPa->Pt());					//Eta Prime		  	  
  		} else
    		continue;
		
		
		//}	
 	
 	
   		
 		fHPtSpecAccPhoton->Fill(mcPa->Pt());	
		fHPtSpecEtaPhoton->Fill(mcPa->Eta());
		fHPtSpecPhiPhoton->Fill(mcPa->Phi());
		fHGenEtaPhi->Fill(mcPa->Eta(),mcPa->Phi());
		
		
		
		if(mcPa->Eta() > -2.5 && mcPa->Eta() < 2.5){fHPtSpecEffPhotonEta5->Fill(mcPa->Pt());}	  		
		if(mcPa->Eta() > -2.0 && mcPa->Eta() < 2.0){fHPtSpecEffPhotonEta4->Fill(mcPa->Pt());}	  		
		if(mcPa->Eta() > -1.5 && mcPa->Eta() < 1.5){fHPtSpecEffPhotonEta3->Fill(mcPa->Pt());}	  			
		if(mcPa->Eta() > -1.0 && mcPa->Eta() < 1.0){fHPtSpecEffPhotonEta2->Fill(mcPa->Pt());}	  	
		if(mcPa->Eta() > -0.5 && mcPa->Eta() < 0.5){fHPtSpecEffPhotonEta1->Fill(mcPa->Pt());}	
		
	
	
	  	//Is EMCal?  
		if(mcPa->Phi() < phimin || mcPa->Phi() > phimax){continue;}
		if(mcPa->Eta() < etamin || mcPa->Eta() > etamax){continue;}	  	 // -0.7 - 0.7
		//cout<<"etamin: "<<phimin<<" -- etamax: "<<phimax<<endl;
	
	
		//Photonen auf EMCal 
		fHPtSpecEffPhoton->Fill(mcPa->Pt());	 	
		
	
	
		//Conversions
		//count photons, that fly to the EMCAl but convert before the detector
		Int_t jdaughter = mcPa->GetDaughterFirst();
		AliMCParticle* mcDa = static_cast<AliMCParticle*>(mcEvent->GetTrack(jdaughter));
		if(!mcDa){continue;}
	
	Double_t dRd = TMath::Sqrt((mcDa->Xv()-evtVtx->GetX())*(mcDa->Xv()-evtVtx->GetX()) +
  	                        (mcDa->Yv()-evtVtx->GetY())*(mcDa->Yv()-evtVtx->GetY()));
		//if(dRd > 420){continue;}
	//if(mcDa->PdgCode() != 11 && mcDa->PdgCode() != -11) {cout<<mcDa->PdgCode()<<endl;}
	
  	if(mcDa->PdgCode() == 11 || mcDa->PdgCode() == -11){
		if(dRd < 420){fHPtSpecConversion->Fill(mcPa->Pt());}
		if(dRd > 420){fHPtSpecConversionNot->Fill(mcPa->Pt());}
	}
			
	
	
		



}//end of track loop
  
  
	
	
	
	
  
  /*
  
    if(bGen && !bAddEta && !bAddPi0){
      // if both photons are on EMCal
      if(binp){
        if(mcP->PdgCode() == 111 && bacc){
          fHPionTruthPtAcc->Fill(clusterVec.Pt());
        }
        if(mcP->PdgCode() == 221 && bacc){
          fHEtaTruthPtAcc->Fill(clusterVec.Pt());
        }
      
        if(mcP->PdgCode() == 111 && baccconv){
          fHPionTruthPtConvAcc->Fill(clusterVec.Pt());
        }
        if(mcP->PdgCode() == 221 && baccconv){
          fHEtaTruthPtConvAcc->Fill(clusterVec.Pt());
        }
      }
    }
    
    
    if(bAddPi0 && !bAddEta && !bGen){
      // if both photons are on EMCal
      if(binp){
        if(mcP->PdgCode() == 111 && bacc){
          fHPionTruthPtAccAdd->Fill(clusterVec.Pt(),wgt2);
        }
      
        if(mcP->PdgCode() == 111 && baccconv){
          fHPionTruthPtConvAccAdd->Fill(clusterVec.Pt(),wgt2);
        }
      }
    }
    
    
    if(bAddEta && !bAddPi0 && !bGen){
      // if both photons are on EMCal
      if(binp){
        if(mcP->PdgCode() == 221 && bacc){
          fHEtaTruthPtAccAdd->Fill(clusterVec.Pt(),wgt2);
        }
        if(mcP->PdgCode() == 221 && baccconv){
          fHEtaTruthPtConvAccAdd->Fill(clusterVec.Pt(),wgt2);
        }
      }
      
    }
    ProcessDaughters(mcP, iTrack, mcEvent);
	*/
  //}
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::FillNtuple()
{
  // Fill ntuple.
  
  if (!fNtuple)
    return;
  
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (fAodEv) {
    AliAODHeader * aodheader = dynamic_cast<AliAODHeader*>(fAodEv->GetHeader());
    if(!aodheader) AliFatal("Not a standard AOD");
 
    fHeader->fRun            = fAodEv->GetRunNumber();
    fHeader->fOrbit          = aodheader->GetOrbitNumber(); 
    fHeader->fPeriod         = aodheader->GetPeriodNumber();
    fHeader->fBx             = aodheader->GetBunchCrossNumber();
    fHeader->fL0             = aodheader->GetL0TriggerInputs();
    fHeader->fL1             = aodheader->GetL1TriggerInputs();
    fHeader->fL2             = aodheader->GetL2TriggerInputs();
    fHeader->fTrClassMask    = aodheader->GetTriggerMask();
    fHeader->fTrCluster      = aodheader->GetTriggerCluster();
    fHeader->fOffTriggers    = aodheader->GetOfflineTrigger();
    fHeader->fFiredTriggers  = aodheader->GetFiredTriggerClasses();
  } else {
    fHeader->fRun            = fEsdEv->GetRunNumber();
    fHeader->fOrbit          = fEsdEv->GetHeader()->GetOrbitNumber();
    fHeader->fPeriod         = fEsdEv->GetESDRun()->GetPeriodNumber();
    fHeader->fBx             = fEsdEv->GetHeader()->GetBunchCrossNumber();
    fHeader->fL0             = fEsdEv->GetHeader()->GetL0TriggerInputs();
    fHeader->fL1             = fEsdEv->GetHeader()->GetL1TriggerInputs();
    fHeader->fL2             = fEsdEv->GetHeader()->GetL2TriggerInputs();
    fHeader->fTrClassMask    = fEsdEv->GetHeader()->GetTriggerMask();
    fHeader->fTrCluster      = fEsdEv->GetHeader()->GetTriggerCluster();
    fHeader->fOffTriggers    = ((AliInputEventHandler*)(am->GetInputEventHandler()))->IsEventSelected();
    fHeader->fFiredTriggers  = fEsdEv->GetFiredTriggerClasses();
    Float_t v0CorrR = 0;
    fHeader->fV0 = AliESDUtils::GetCorrV0(fEsdEv,v0CorrR);
    const AliMultiplicity *mult = fEsdEv->GetMultiplicity();
    if (mult)
      fHeader->fCl1 = mult->GetNumberOfITSClusters(1);
    fHeader->fTr = AliESDtrackCuts::GetReferenceMultiplicity(fEsdEv,1);
    AliTriggerAnalysis trAn; /// Trigger Analysis
    Bool_t v0B = trAn.IsOfflineTriggerFired(fEsdEv, AliTriggerAnalysis::kV0C);
    Bool_t v0A = trAn.IsOfflineTriggerFired(fEsdEv, AliTriggerAnalysis::kV0A);
    fHeader->fV0And        = v0A && v0B;
    fHeader->fIsHT         = (fHeader->fOffTriggers & AliVEvent::kEMC1) || (fHeader->fOffTriggers & AliVEvent::kEMC7);
    am->LoadBranch("SPDPileupVertices");
    am->LoadBranch("TrkPileupVertices");
    fHeader->fIsPileup     = fEsdEv->IsPileupFromSPD(3,0.8);
    fHeader->fIsPileup2    = fEsdEv->IsPileupFromSPD(3,0.4);
    fHeader->fIsPileup4    = fEsdEv->IsPileupFromSPD(3,0.2);
    fHeader->fIsPileup8    = fEsdEv->IsPileupFromSPD(3,0.1);
    fHeader->fNSpdVertices = fEsdEv->GetNumberOfPileupVerticesSPD();
    fHeader->fNTpcVertices = fEsdEv->GetNumberOfPileupVerticesTracks();
  }
  
  AliCentrality *cent = InputEvent()->GetCentrality();
  fHeader->fV0Cent    = cent->GetCentralityPercentileUnchecked("V0M");
  fHeader->fCl1Cent   = cent->GetCentralityPercentileUnchecked("CL1");
  fHeader->fTrCent    = cent->GetCentralityPercentileUnchecked("TRK");
  fHeader->fCqual     = cent->GetQuality();
  
  AliEventplane *ep = InputEvent()->GetEventplane();
  if (ep) {
    if (ep->GetQVector())
      fHeader->fPsi     = ep->GetQVector()->Phi()/2. ;
    else
      fHeader->fPsi = -1;
    if (ep->GetQsub1()&&ep->GetQsub2())
      fHeader->fPsiRes  = ep->GetQsub1()->Phi()/2.-ep->GetQsub2()->Phi()/2.;
    else
      fHeader->fPsiRes = 0;
  }
  
  Double_t val = 0;
  TString trgclasses(fHeader->fFiredTriggers);
  for (Int_t j = 0; j<fTrClassNamesArr->GetEntries(); ++j) {
    const char *name = fTrClassNamesArr->At(j)->GetName();
    TRegexp regexp(name);
    if (trgclasses.Contains(regexp))
      val += TMath::Power(2,j);
  }
  fHeader->fTcls = (UInt_t)val;
  
  fHeader->fNSelTr     = fSelTracks->GetEntries();
  fHeader->fNSelPrimTr = fSelPrimTracks->GetEntries();
  fHeader->fNSelPrimTr1   = 0;
  fHeader->fNSelPrimTr2   = 0;
  for(int iTracks=0; iTracks < fSelPrimTracks->GetEntries(); iTracks++){
    AliESDtrack *track = static_cast<AliESDtrack*>(fSelPrimTracks->At(iTracks));
    if(track->Pt()>1)
      ++fHeader->fNSelPrimTr1;
    if(track->Pt()>2)
      ++fHeader->fNSelPrimTr2;
  }
  
  fHeader->fNCells   = 0;
  fHeader->fNCells0  = 0;
  fHeader->fNCells01 = 0;
  fHeader->fNCells03 = 0;
  fHeader->fNCells1  = 0;
  fHeader->fNCells2  = 0;
  fHeader->fNCells5  = 0;
  fHeader->fMaxCellE = 0;
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  
  if (cells) {
    Int_t ncells = cells->GetNumberOfCells();
    for(Int_t j=0; j<ncells; ++j) {
      Double_t cellen = cells->GetAmplitude(j);
      if (cellen>0.045)
        ++fHeader->fNCells0;
      if (cellen>0.1)
        ++fHeader->fNCells01;
      if (cellen>0.3)
        ++fHeader->fNCells03;
      if (cellen>1)
        ++fHeader->fNCells1;
      if (cellen>2)
        ++fHeader->fNCells2;
      if (cellen>5)
        ++fHeader->fNCells5;
      if (cellen>fHeader->fMaxCellE)
        fHeader->fMaxCellE = cellen;
    }
    fHeader->fNCells = ncells;
  }
  
  fHeader->fNClus      = 0;
  fHeader->fNClus1     = 0;
  fHeader->fNClus2     = 0;
  fHeader->fNClus5     = 0;
  fHeader->fMaxClusE   = 0;
  
  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  
  if (clusters) {
    Int_t nclus = clusters->GetEntries();
    for(Int_t j=0; j<nclus; ++j) {
      AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(j));
      if (!clus->IsEMCAL())
        continue;
      Double_t clusen = clus->E();
      if (clusen>1)
        ++fHeader->fNClus1;
      if (clusen>2)
        ++fHeader->fNClus2;
      if (clusen>5)
        ++fHeader->fNClus5;
      if (clusen>fHeader->fMaxClusE)
        fHeader->fMaxClusE = clusen;
    }
    fHeader->fNClus = nclus;
  }
  
  fHeader->fMaxTrE     = 0;
  if (fTriggers) {
    Int_t ntrig = fTriggers->GetEntries();
    for (Int_t j = 0; j<ntrig; ++j) {
      AliStaTrigger *sta = static_cast<AliStaTrigger*>(fTriggers->At(j));
      if (!sta)
        continue;
      if (sta->fE>fHeader->fMaxTrE)
        fHeader->fMaxTrE = sta->fE;
    }
  }
  
  // count cells above 100 MeV on super modules
  fHeader->fNcSM0 = GetNCells(0, 0.100);
  fHeader->fNcSM1 = GetNCells(1, 0.100);
  fHeader->fNcSM2 = GetNCells(2, 0.100);
  fHeader->fNcSM3 = GetNCells(3, 0.100);
  fHeader->fNcSM4 = GetNCells(4, 0.100);
  fHeader->fNcSM5 = GetNCells(5, 0.100);
  fHeader->fNcSM6 = GetNCells(6, 0.100);
  fHeader->fNcSM7 = GetNCells(7, 0.100);
  fHeader->fNcSM8 = GetNCells(8, 0.100);
  fHeader->fNcSM9 = GetNCells(9, 0.100);
  
  if (fAodEv) {
    am->LoadBranch("vertices");
    AliAODVertex *pv = fAodEv->GetPrimaryVertex();
    FillVertex(fPrimVert, pv);
    AliAODVertex *sv = fAodEv->GetPrimaryVertexSPD();
    FillVertex(fSpdVert, sv);
  } else {
    am->LoadBranch("PrimaryVertex.");
    const AliESDVertex *pv = fEsdEv->GetPrimaryVertexTracks();
    FillVertex(fPrimVert, pv);
    am->LoadBranch("SPDVertex.");
    const AliESDVertex *sv = fEsdEv->GetPrimaryVertexSPD();
    FillVertex(fSpdVert, sv);
    am->LoadBranch("TPCVertex.");
    const AliESDVertex *tv = fEsdEv->GetPrimaryVertexTPC();
    FillVertex(fTpcVert, tv);
  }
  
  fNtuple->Fill();
}


//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::FillMcHists()
{
  // Fill additional MC information histograms.
  
  if (!fMcParts)
    return;
  
  // check if aod or esd mc mode and the fill histos
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::FillOtherHists()
{
  // Fill other histograms.
  
  for(int iTracks=0; iTracks < fSelPrimTracks->GetEntries(); ++iTracks){
    AliESDtrack *track = static_cast<AliESDtrack*>(fSelPrimTracks->At(iTracks));
    if(!track)
      continue;
    fHPrimTrackPt->Fill(track->Pt());
    fHPrimTrackEta->Fill(track->Eta());
    fHPrimTrackPhi->Fill(track->Phi());
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::FillTrackHists()
{
  // Fill track histograms.
  
  if (fSelPrimTracks) {
    for(int iTracks=0; iTracks < fSelPrimTracks->GetEntries(); iTracks++) {
      AliESDtrack *track = static_cast<AliESDtrack*>(fSelPrimTracks->At(iTracks));
      if(!track)
        continue;
      fHPrimTrackPt->Fill(track->Pt());
      fHPrimTrackEta->Fill(track->Eta());
      fHPrimTrackPhi->Fill(track->Phi());
    }
  }
}

//__________________________________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::FillVertex(AliStaVertex *v, const AliESDVertex *esdv)
{
  // Fill vertex from ESD vertex info.
  
  v->fVx   = esdv->GetX();
  v->fVy   = esdv->GetY();
  v->fVz   = esdv->GetZ();
  v->fVc   = esdv->GetNContributors();
  v->fDisp = esdv->GetDispersion();
  v->fZres = esdv->GetZRes();
  v->fChi2 = esdv->GetChi2();
  v->fSt   = esdv->GetStatus();
  v->fIs3D = esdv->IsFromVertexer3D();
  v->fIsZ  = esdv->IsFromVertexerZ();
}

//__________________________________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::FillVertex(AliStaVertex *v, const AliAODVertex *aodv)
{
  // Fill vertex from AOD vertex info.
  
  v->fVx   = aodv->GetX();
  v->fVy   = aodv->GetY();
  v->fVz   = aodv->GetZ();
  v->fVc   = aodv->GetNContributors();
  v->fChi2 = aodv->GetChi2();
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALDirGamma::GetCellIsolation(Double_t cEta, Double_t cPhi, Double_t radius) const
{
  // Compute isolation based on cell content.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  if (!cells)
    return 0;
  
  Double_t cellIsolation = 0;
  Double_t rad2 = radius*radius;
  Int_t ncells = cells->GetNumberOfCells();
  for (Int_t i = 0; i<ncells; ++i) {
    Int_t absID    = TMath::Abs(cells->GetCellNumber(i));
    Float_t eta=-1, phi=-1;
    fGeom->EtaPhiFromIndex(absID,eta,phi);
    Double_t phidiff = TVector2::Phi_mpi_pi(phi-cPhi);
    Double_t dist = (eta-cEta)*(eta-cEta)+phidiff*phidiff;
    if(dist>rad2)
      continue;
    Double_t cellE = cells->GetAmplitude(i);
    Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
    Double_t cellEt = cellE*sin(theta);
    cellIsolation += cellEt;
  }
  return cellIsolation;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALDirGamma::GetCellIsoNxM(Double_t cEta, Double_t cPhi, Int_t N, Int_t M) const
{
  // Compute isolation based on cell content, in a NxM rectangle.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  if (!cells)
    return 0;
  
  Double_t cellIsolation = 0;
  Int_t ncells = cells->GetNumberOfCells();
  for (Int_t i = 0; i<ncells; ++i) {
    Int_t absID    = TMath::Abs(cells->GetCellNumber(i));
    Float_t eta=-1, phi=-1;
    fGeom->EtaPhiFromIndex(absID,eta,phi);
    Double_t phidiff = TVector2::Phi_mpi_pi(phi-cPhi);
    Double_t etadiff = (eta-cEta)*(eta-cEta)+phidiff*phidiff;
    if(TMath::Abs(etadiff)/0.014>N)
      continue;
    if(TMath::Abs(phidiff)/0.014>M)
      continue;
    Double_t cellE = cells->GetAmplitude(i);
    Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
    Double_t cellEt = cellE*sin(theta);
    cellIsolation += cellEt;
  }
  return cellIsolation;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALDirGamma::GetCellEnergy(const AliVCluster *cluster) const
{
  // Get maximum energy of attached cell.
  
  Double_t ret = 0;
  Int_t ncells = cluster->GetNCells();
  if (fEsdCells) {
    for (Int_t i=0; i<ncells; i++) {
      Double_t e = fEsdCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
      ret += e;
    }
  } else {
    for (Int_t i=0; i<ncells; i++) {
      Double_t e = fAodCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
      ret += e;
    }
  }
  return ret;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALDirGamma::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
{
  // Get maximum energy of attached cell.
  
  id = -1;
  Double_t maxe = 0;
  Int_t ncells = cluster->GetNCells();
  if (fEsdCells) {
    for (Int_t i=0; i<ncells; i++) {
      Double_t e = fEsdCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
      if (e>maxe) {
        maxe = e;
        id   = cluster->GetCellAbsId(i);
      }
    }
  } else {
    for (Int_t i=0; i<ncells; i++) {
      Double_t e = fAodCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
      if (e>maxe)
        maxe = e;
      id   = cluster->GetCellAbsId(i);
    }
  }
  return maxe;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALDirGamma::GetSecondMaxCellEnergy(AliVCluster *clus, Short_t &id) const
{
  // Get second maximum cell.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  if (!cells)
    return -1;
  
  Double_t secondEmax=0, firstEmax=0;
  Double_t cellen;
  for(Int_t iCell=0;iCell<clus->GetNCells();iCell++){
    Int_t absId = clus->GetCellAbsId(iCell);
    cellen = cells->GetCellAmplitude(absId);
    if(cellen > firstEmax)
      firstEmax = cellen;
  }
  for(Int_t iCell=0;iCell<clus->GetNCells();iCell++){
    Int_t absId = clus->GetCellAbsId(iCell);
    cellen = cells->GetCellAmplitude(absId);
    if(cellen < firstEmax && cellen > secondEmax) {
      secondEmax = cellen;
      id = absId;
    }
  }
  return secondEmax;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::GetSigma(const AliVCluster *c, Double_t& sigmaMax, Double_t &sigmaMin) const
{
  // Calculate the (E) weighted variance along the longer (eigen) axis.
  
  sigmaMax = 0;          // cluster variance along its longer axis
  sigmaMin = 0;          // cluster variance along its shorter axis
  Double_t Ec  = c->E(); // cluster energy
  if(Ec<=0)
    return;
  Double_t Xc  = 0 ;     // cluster first moment along X
  Double_t Yc  = 0 ;     // cluster first moment along Y
  Double_t Sxx = 0 ;     // cluster second central moment along X (variance_X^2)
  Double_t Sxy = 0 ;     // cluster second central moment along Y (variance_Y^2)
  Double_t Syy = 0 ;     // cluster covariance^2
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  
  if (!cells)
    return;
  
  Int_t ncells = c->GetNCells();
  if (ncells==1)
    return;
  
  for(Int_t j=0; j<ncells; ++j) {
    Int_t id = TMath::Abs(c->GetCellAbsId(j));
    Double_t cellen = cells->GetCellAmplitude(id);
    TVector3 pos;
    fGeom->GetGlobal(id,pos);
    Xc  += cellen*pos.X();
    Yc  += cellen*pos.Y();
    Sxx += cellen*pos.X()*pos.X();
    Syy += cellen*pos.Y()*pos.Y();
    Sxy += cellen*pos.X()*pos.Y();
  }
  Xc  /= Ec;
  Yc  /= Ec;
  Sxx /= Ec;
  Syy /= Ec;
  Sxy /= Ec;
  Sxx -= Xc*Xc;
  Syy -= Yc*Yc;
  Sxy -= Xc*Yc;
  Sxx = TMath::Abs(Sxx);
  Syy = TMath::Abs(Syy);
  sigmaMax = (Sxx + Syy + TMath::Sqrt(TMath::Abs((Sxx-Syy)*(Sxx-Syy)+4.0*Sxy*Sxy)))/2.0;
  sigmaMax = TMath::Sqrt(TMath::Abs(sigmaMax));
  sigmaMin = TMath::Abs(Sxx + Syy - TMath::Sqrt(TMath::Abs((Sxx-Syy)*(Sxx-Syy)+4.0*Sxy*Sxy)))/2.0;
  sigmaMin = TMath::Sqrt(TMath::Abs(sigmaMin));
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::GetSigmaEtaEta(const AliVCluster *c, Double_t& sEtaEta, Double_t &sPhiPhi) const
{
  // Calculate the (E) weighted variance along the pseudorapidity.
  
  sEtaEta = 0;
  sPhiPhi = 0;
  
  Double_t Ec  = c->E(); // cluster energy
  if(Ec<=0)
    return;
  
  const Int_t ncells = c->GetNCells();
  
  Double_t EtaC    = 0;  // cluster first moment along eta
  Double_t PhiC    = 0;  // cluster first moment along phi
  Double_t Setaeta = 0;  // cluster second central moment along eta
  Double_t Sphiphi = 0;  // cluster second central moment along phi
  Double_t w[ncells];    // weight max(0,4.5*log(E_i/Ec))
  Double_t sumw = 0;
  Int_t id[ncells];
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  
  if (!cells)
    return;
  
  if (ncells==1)
    return;
  
  for(Int_t j=0; j<ncells; ++j) {
    id[j] = TMath::Abs(c->GetCellAbsId(j));
    Double_t cellen = cells->GetCellAmplitude(id[j]);
    w[j] = TMath::Max(0., 4.5+TMath::Log(cellen/Ec));
    TVector3 pos;
    fGeom->GetGlobal(id[j],pos);
    EtaC += w[j]*pos.Eta();
    PhiC += w[j]*pos.Phi();
    sumw += w[j];
  }
  EtaC /= sumw;
  PhiC /= sumw;
  
  for(Int_t j=0; j<ncells; ++j) {
    TVector3 pos;
    fGeom->GetGlobal(id[j],pos);
    Setaeta =  w[j]*(pos.Eta() - EtaC)*(pos.Eta() - EtaC);
    Sphiphi =  w[j]*(pos.Phi() - PhiC)*(pos.Phi() - PhiC);
  }
  Setaeta /= sumw;
  sEtaEta = TMath::Sqrt(Setaeta);
  Sphiphi /= sumw;
  sPhiPhi = TMath::Sqrt(Sphiphi);
}

//________________________________________________________________________
Int_t AliAnalysisTaskEMCALDirGamma::GetNCells(const AliVCluster *c, Double_t emin) const
{
  // Calculate number of attached cells above emin.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  if (!cells)
    return 0;
  
  Int_t n = 0;
  Int_t ncells = c->GetNCells();
  for(Int_t j=0; j<ncells; ++j) {
    Int_t id = TMath::Abs(c->GetCellAbsId(j));
    Double_t cellen = cells->GetCellAmplitude(id);
    if (cellen>=emin)
      ++n;
  }
  return n;
}

//________________________________________________________________________
Int_t AliAnalysisTaskEMCALDirGamma::GetNCells(Int_t sm, Double_t emin) const
{
  // Calculate number of cells per SM above emin.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  if (!cells)
    return 0;
  
  Int_t n = 0;
  Int_t ncells = cells->GetNumberOfCells();
  for(Int_t j=0; j<ncells; ++j) {
    Int_t id = TMath::Abs(cells->GetCellNumber(j));
    Double_t cellen = cells->GetCellAmplitude(id);
    if (cellen<emin)
      continue;
    Int_t fsm = fGeom->GetSuperModuleNumber(id);
    if (fsm != sm)
      continue;
    ++n;
  }
  return n;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALDirGamma::GetTrackIsolation(Double_t cEta, Double_t cPhi, Double_t radius, Double_t pt) const
{
  // Compute isolation based on tracks.
  
  Double_t trkIsolation = 0;
  Double_t rad2 = radius*radius;
  Int_t ntrks = fSelPrimTracks->GetEntries();
  for(Int_t j = 0; j<ntrks; ++j) {
    AliVTrack *track = static_cast<AliVTrack*>(fSelPrimTracks->At(j));
    if (!track)
      continue;
    if (track->Pt()<pt)
      continue;
    Float_t eta = track->Eta();
    Float_t phi = track->Phi();
    Double_t phidiff = TVector2::Phi_mpi_pi(phi-cPhi);
    Double_t dist = (eta-cEta)*(eta-cEta)+phidiff*phidiff;
    if(dist>rad2)
      continue;
    trkIsolation += track->Pt();
  }
  return trkIsolation;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALDirGamma::GetTrackIsoStrip(Double_t cEta, Double_t cPhi, Double_t dEta, Double_t dPhi, Double_t pt) const
{
  // Compute isolation based on tracks.
  
  Double_t trkIsolation = 0;
  Int_t ntrks = fSelPrimTracks->GetEntries();
  for(Int_t j = 0; j<ntrks; ++j) {
    AliVTrack *track = static_cast<AliVTrack*>(fSelPrimTracks->At(j));
    if (!track)
      continue;
    if (track->Pt()<pt)
      continue;
    Float_t eta = track->Eta();
    Float_t phi = track->Phi();
    Double_t phidiff = TVector2::Phi_mpi_pi(phi-cPhi);
    Double_t etadiff = (eta-cEta);
    if(TMath::Abs(etadiff)>dEta || TMath::Abs(phidiff)>dPhi)
      continue;
    trkIsolation += track->Pt();
  }
  return trkIsolation;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEMCALDirGamma::IsShared(const AliVCluster *c) const
{
  // Returns if cluster shared across super modules.
  
  if (!c)
    return 0;
  
  Int_t n = -1;
  Int_t ncells = c->GetNCells();
  for(Int_t j=0; j<ncells; ++j) {
    Int_t id = TMath::Abs(c->GetCellAbsId(j));
    Int_t got = id / (24*48);
    if (n==-1) {
      n = got;
      continue;
    }
    if (got!=n)
      return 1;
  }
  return 0;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEMCALDirGamma::IsIdPartOfCluster(const AliVCluster *c, Short_t id) const
{
  // Returns if id is part of cluster.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  if (!cells)
    return 0;
  
  Int_t ncells = c->GetNCells();
  for(Int_t j=0; j<ncells; ++j) {
    Int_t cid = TMath::Abs(c->GetCellAbsId(j));
    if (cid == id)
      return 1;
  }
  return 0;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::PrintDaughters(const AliVParticle *p, const TObjArray *arr, Int_t level) const
{
  // Print recursively daughter information.
  
  if (!p || !arr)
    return;
  
  const AliAODMCParticle *amc = dynamic_cast<const AliAODMCParticle*>(p);
  if (!amc)
    return;
  for (Int_t i=0; i<level; ++i) printf(" ");
  amc->Print();
  
  Int_t n = amc->GetNDaughters();
  for (Int_t i=0; i<n; ++i) {
    Int_t d = amc->GetDaughterLabel(i);
    const AliVParticle *dmc = static_cast<const AliVParticle*>(arr->At(d));
    PrintDaughters(dmc,arr,level+1);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::PrintDaughters(const AliMCParticle *p, const AliMCEvent *arr, Int_t level) const
{
  // Print recursively daughter information.
  
  if (!p || !arr)
    return;
  
  for (Int_t i=0; i<level; ++i) printf(" ");
  Int_t d1 = p->GetDaughterFirst();
  Int_t d2 = p->GetDaughterLast();
  printf("pid=%d: %.2f %.2f %.2f (%.2f %.2f %.2f); nd=%d,%d\n",
         p->PdgCode(),p->Px(),p->Py(),p->Pz(),p->Xv(),p->Yv(),p->Zv(),d1,d2);
  if (d1<0)
    return;
  if (d2<0)
    d2=d1;
  for (Int_t i=d1;i<=d2;++i) {
    const AliMCParticle *dmc = static_cast<const AliMCParticle *>(arr->GetTrack(i));
    PrintDaughters(dmc,arr,level+1);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::ProcessDaughters(AliVParticle *p, Int_t index, const TObjArray *arr)
{
  // Process and create daughters.
  
  if (!p || !arr)
    return;
  
  AliAODMCParticle *amc = dynamic_cast<AliAODMCParticle*>(p);
  if (!amc)
    return;
  
  //amc->Print();
  
  Int_t nparts = arr->GetEntries();
  Int_t nents  = fMcParts->GetEntries();
  
  AliStaPart *newp = static_cast<AliStaPart*>(fMcParts->New(nents));
  newp->fPt  = amc->Pt();
  newp->fEta = amc->Eta();
  newp->fPhi = amc->Phi();
  if (amc->Xv() != 0 || amc->Yv() != 0 || amc->Zv() != 0) {
    TVector3 vec(amc->Xv(),amc->Yv(),amc->Zv());
    newp->fVR = vec.Perp();
    newp->fVEta = vec.Eta();
    newp->fVPhi = vec.Phi();
  }
  newp->fPid  = amc->PdgCode();
  newp->fLab  = nents;
  Int_t moi = amc->GetMother();
  if (moi>=0&&moi<nparts) {
    const AliAODMCParticle *mmc = static_cast<const AliAODMCParticle*>(arr->At(moi));
    moi = mmc->GetUniqueID();
  }
  newp->fMo = moi;
  p->SetUniqueID(nents);
  
  // TODO: Determine which detector was hit
  //newp->fDet = ???
  
  Int_t n = amc->GetNDaughters();
  for (Int_t i=0; i<n; ++i) {
    Int_t d = amc->GetDaughterLabel(i);
    if (d<=index || d>=nparts)
      continue;
    AliVParticle *dmc = static_cast<AliVParticle*>(arr->At(d));
    ProcessDaughters(dmc,d,arr);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::ProcessDaughters(AliMCParticle *p, Int_t index, const AliMCEvent *arr)
{
  // Process and create daughters.
  
  if (!p || !arr)
    return;
  
  Int_t d1 = p->GetDaughterFirst();
  Int_t d2 = p->GetDaughterLast();
  if (0) {
    printf("%d pid=%d: %.3f %.3f %.3f (%.2f %.2f %.2f); nd=%d,%d, mo=%d\n",
           index,p->PdgCode(),p->Px(),p->Py(),p->Pz(),p->Xv(),p->Yv(),p->Zv(),d1,d2, p->GetMother());
  }
  Int_t nents  = fMcParts->GetEntries();
  
  AliStaPart *newp = static_cast<AliStaPart*>(fMcParts->New(nents));
  newp->fPt  = p->Pt();
  newp->fEta = p->Eta();
  newp->fPhi = p->Phi();
  if (p->Xv() != 0 || p->Yv() != 0 || p->Zv() != 0) {
    TVector3 vec(p->Xv(),p->Yv(),p->Zv());
    newp->fVR = vec.Perp();
    newp->fVEta = vec.Eta();
    newp->fVPhi = vec.Phi();
  }
  newp->fPid  = p->PdgCode();
  newp->fLab  = nents;
  Int_t moi = p->GetMother();
  if (moi>=0) {
    const AliMCParticle *mmc = static_cast<const AliMCParticle *>(arr->GetTrack(moi));
    moi = mmc->GetUniqueID();
  }
  newp->fMo = moi;
  p->SetUniqueID(nents);
  
  Int_t nref = p->GetNumberOfTrackReferences();
  if (nref>0) {
    AliTrackReference *ref = p->GetTrackReference(nref-1);
    if (ref) {
      newp->fDet = ref->DetectorId();
    }
  }
  
  if (d1<0)
    return;
  if (d2<0)
    d2=d1;
  for (Int_t i=d1;i<=d2;++i) {
    AliMCParticle *dmc = static_cast<AliMCParticle *>(arr->GetTrack(i));
    if (dmc->P()<0.01)
      continue;
    ProcessDaughters(dmc,i,arr);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::GetMulClass(Int_t& imcl)
{
  Int_t nclus = 0;
  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  
  if (clusters)
    nclus = clusters->GetEntries();
  
  //const int MultCut[8] = {5, 15, 30, 50, 80, 120, 300, 9999};
  const int MultCut[nMulClass] = {5, 12, 20, 9999};
  
  imcl=0;
  
  for (imcl=0; imcl<nMulClass; imcl++) {
    if (nclus < MultCut[imcl]) break;
  }
}


//_____________________________________________________________________
void AliAnalysisTaskEMCALDirGamma::SetDnDpT(Int_t i, Double_t par0, Double_t par1, Double_t par2, Double_t par3, Double_t par4)
{
  /*
  // case 1: modified hagedorn
  if(i==1){
    fPi0DnDpt = new TF1("fPi0DnDpt","[0]*pow([1]/(([1]*exp(-[3]*x)+x)),[2])",0.6,12);
    fPi0DnDpt->SetParameters(par0,par1,par2,par3);
  
    par4 += 0;
  }
  
  // get rid of the warnings
  else{
    par0+=0;
    par1+=0;
    par2+=0;
    par3+=0;
    par4+=0;
  }
  */
  return;
}

//_____________________________________________________________________
Double_t AliAnalysisTaskEMCALDirGamma::CalcWeight(Double_t pt, Double_t eta, Int_t i){
  Double_t weight = 1.;
  if(i==1){
    Double_t par0 = 2.52684e08;
    Double_t par1 = 0.730803;
    Double_t par2 = 5.32059;
    Double_t par3 = 0.548711;

    weight = (par0*pow(par1/((par1*exp(-par3*pt)+pt)),par2)) * (7.93979e-01 + eta*1.31404e-03 + eta*eta*4.96607e-01);
  }

  return weight;
}

//_____________________________________________________________________
Int_t AliAnalysisTaskEMCALDirGamma::GetModuleNumber(AliVCluster * cluster) const
{
  //Get the EMCAL/PHOS module number that corresponds to this cluster
  TLorentzVector lv;
  Double_t v[]={0.,0.,0.}; //not necessary to pass the real vertex.
  if(!cluster)
  {
    if(fDebug > 1) printf("AliCalorimeterUtils::GetModuleNumber() - NUL Cluster, please check!!!");
    return -1;
  }
  
  cluster->GetMomentum(lv,v);
  Float_t phi = lv.Phi();
  if(phi < 0) phi+=TMath::TwoPi();
  Int_t absId = -1;
  if(cluster->IsEMCAL()){
    fGeom->GetAbsCellIdFromEtaPhi(lv.Eta(),phi, absId);
    if(fDebug > 2)
      printf("AliCalorimeterUtils::GetModuleNumber() - EMCAL: cluster eta %f, phi %f, absid %d, SuperModule %d\n",
             lv.Eta(), phi*TMath::RadToDeg(),absId, fGeom->GetSuperModuleNumber(absId));
    return fGeom->GetSuperModuleNumber(absId) ;
  }//EMCAL
  else if(cluster->IsPHOS())
  {
    //    Int_t    relId[4];
    //    if ( cluster->GetNCells() > 0)
    //    {
    //      absId = cluster->GetCellAbsId(0);
    //      if(fDebug > 2)
    //        printf("AliCalorimeterUtils::GetModuleNumber() - PHOS: cluster eta %f, phi %f, e %f, absId %d\n",
    //               lv.Eta(), phi*TMath::RadToDeg(), lv.E(), absId);
    //    }
    //    else return -1;
    //
    //    if ( absId >= 0)
    //    {
    //      fPHOSGeo->AbsToRelNumbering(absId,relId);
    //      if(fDebug > 2)
    //        printf("AliCalorimeterUtils::GetModuleNumber() - PHOS: Module %d\n",relId[0]-1);
    //      return relId[0]-1;
    //    }
    //    else return -1;
    return -1;
  }//PHOS
  
  return -1;
}

// Jason's energy recalibration, use calib factor 8 for data and 9 for MC
//________________________________________________________________________
Double_t AliAnalysisTaskEMCALDirGamma::PrivateEnergyRecal(Double_t energy, Int_t iCalib){
  
  double recalibfactor = 0.0;
  
  if(iCalib==0){// no recalibration!
    recalibfactor = 1.0;
  }
  else if(iCalib==1){// just a scale factor:
    recalibfactor = 0.984;
  }
  else if(iCalib==2){// Symmetric Decay Fit - corrects data to uncorrected MC.
    Double_t p[3] = {0.96968, -2.68720, -0.831607};
    recalibfactor = p[0] + exp(p[1] + p[2]*energy*2.0);
  }
  else if(iCalib==3){// Jason's fit to the LHC12f1a MC single photons - 04 Aug 2013 (call it kPi0MCv4??)
    Double_t p[7] = {1.00000e+00, 3.04925e-02, 4.69043e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.00046e+00};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==4){// Jason's fit to the test beam data - 04 Aug 2013(call it kBTCv3??)
    Double_t p[7] = {9.78672e-01, 2.39745e-01, 6.41199e-01, 9.13538e-02, 1.46058e+02, 1.99469e+01, 9.72716e-01};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==5){// Based on kSDM/kTBCv3 (call it kPi0MCv4??)
    Double_t p[10] = {9.78672e-01, 2.39745e-01, 6.41199e-01, 9.13538e-02, 1.46058e+02, 1.99469e+01, 9.72716e-01, 0.96968, -2.68720, -0.831607};
    recalibfactor = ( (p[6]/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5]))))) ) / ( p[7] + exp(p[8] + p[9]*energy/2.0) );
  }
  else if(iCalib==6){// kBeamTestCorrectedv2 - in AliROOT!
    Double_t p[7] = {9.83504e-01, 2.10106e-01, 8.97274e-01, 8.29064e-02, 1.52299e+02, 3.15028e+01, 0.968};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==7){// kPi0MCv3 - in AliROOT!
    Double_t p[7] = {9.81039e-01, 1.13508e-01, 1.00173e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.0};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==8){// Jason's fit to the noNL MC/data- based on kSDM and kPi0MCv5 - 28 Oct 2013 (call it... ??)
    Double_t p[10] = {1.0, 6.64778e-02, 1.57000e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.01286, 0.964, -3.132, -0.435};
    //Double_t p[10] = {1.0, 6.64778e-02, 1.57000e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.01286, 0.96968, -2.68720, -0.831607};//same SDM piece as iCalib==2
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5]))))) * (p[7] + exp(p[8]+p[9]*energy*2.0));
  }
  else if(iCalib==9){// Jason's fit to the LHC12f1a/b MC single photons (above 400MeV), including conversions - 28 Oct 2013 (call it kPi0MCv5??)
    Double_t p[7] = {1.0, 6.64778e-02, 1.57000e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.01286};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==10){// Jason played with test beam data
    Double_t p[7] = {1.0, 0.237767, 0.651203, 0.183741, 155.427, 17.0335, 0.987054};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==11){// Jason played with test beam MC
    Double_t p[7] = {1.0, 0.0797873, 1.68322, 0.0806098, 244.586, 116.938, 1.00437};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  
  return recalibfactor;
}

/*
 //__________________________________________________________________________________________________
 void AliStaCluster::GetMom(TLorentzVector& p, Double_t *vertex)
 {
 // Calculate momentum.
 
 TVector3 pos;
 pos.SetPtEtaPhi(fR,fEta,fPhi);
 
 if(vertex){ //calculate direction relative to vertex
 pos -= vertex;
 }
 
 Double_t r = pos.Mag();
 p.SetPxPyPzE(fE*pos.x()/r, fE*pos.y()/r, fE*pos.z()/r, fE);
 }
 
 //__________________________________________________________________________________________________
 void AliStaCluster::GetMom(TLorentzVector& p, AliStaVertex *vertex)
 {
 // Calculate momentum.
 
 Double_t v[3] = {0,0,0};
 if (vertex) {
 v[0] = vertex->fVx;
 v[1] = vertex->fVy;
 v[2] = vertex->fVz;
 }
 GetMom(p, v);
 }
 */
