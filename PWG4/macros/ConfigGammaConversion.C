/** VERSION NUMBER 1.1 */

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Ana Marin, Kathrin Koch, Kenneth Aamodt
 * Contact: kenneth.aamodt@cern.ch
 * Version 1.1                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
const int c_array_size = 27; // RRnewTOF c_array_size was increased to 27 in order to put a TOF PID cut on electrons

class AliAnalysisDataContainer;
class AliGammaConversionHistograms;
class AliAnalysisTaskGammaConversion;
class AliAnalysisTaskGammaJet;

// set this to a number if you want to analyze a set number of files
// if it is 0 it will analyze the files listed in the data list
Int_t kGCnumberOfFilesToAnalyze=0;

Bool_t kGCrunNeutralMeson = kTRUE;
Bool_t kGCrunJet          = kFALSE;
Bool_t kGCrunChic         = kFALSE;
Bool_t kGCrunDalitz       = kFALSE;
Bool_t kGCrunCF           = kFALSE;
Bool_t kGCcalculateBackground = kTRUE;
Bool_t kGCdoNeutralMesonV0MCCheck =kFALSE;
Bool_t kGCrunOmegaMeson = kFALSE;
Bool_t kGCrunRES = kFALSE;
Bool_t kGCRecalculateV0ForGamma = kFALSE;
Bool_t kGCUseTrackMultiplicityForBG = kTRUE;
Bool_t kGCMoveParticlesAccordingToVertex = kFALSE;
Bool_t kGCLowMemoryConsumption = kFALSE;
Bool_t kGCApplyChi2Cut = kFALSE;

Bool_t kGCUseRotationMethodInBG=kTRUE;
Int_t kGCnDegreeRotationPMForBG=15;
Int_t kGCnumberOfRotationEventsForBG=15;
Bool_t kGCdoBGProbability=kFALSE;
//Svein 
Bool_t kGCRunGammaJetTask = kFALSE;
/** ---------------------------------- define cuts here ------------------------------------*/
TString kGCAnalysisCutSelectionId="900356204010033210220000000"; // do not change here, use -set-cut-selection in argument instead

Int_t kGCNEventsForBGCalculation=20;

Int_t kGCpidOfNegativeTrack=11;
Int_t kGCpidOfPositiveTrack=-11;

Double_t kGCmaxVertexZ   = 10.;
Double_t kGCmaxRCut   = 180.;
Double_t kGCminRCut   = 2.8;
Double_t kGCetaCut    = 0.9;
Double_t kGCptCut     = 0.02;
Double_t kGCsingleptCut = 0.02;
Double_t kGCmaxZCut     = 240.;
Double_t kGCminClsTPCCut= 0.;
Double_t kGCminClsTPCCutToF= 0.;
Double_t kGCchi2CutConversion   = 30.;
Double_t kGCchi2CutMeson   = 50.;
Double_t kGCalphaCutMeson   = 0.7;
Double_t kGCalphaMinCutMeson   = 0.0;
Double_t kGCrapidityCutMeson    = 0.9;

Double_t kGCLineCutZRSlope = tan(2*atan(exp(-kGCetaCut)));
Double_t kGCLineCutZValue = 7.;


Double_t kGCxVertexCut = 0.;
Double_t kGCyVertexCut = 0.;
Double_t kGCzVertexCut = 0.;

Double_t kGCsigmaCutGammaMass=0.0001;

Bool_t kGCuseImprovedVertex = kTRUE;

Bool_t kGCUseOnFlyV0Finder = kTRUE;

// define masses of different particles, this will be used by the KF particle
// together with the width to set mass constraints. Units in GeV.
Double_t kGCelectronMass = 0.00051099892;
Double_t kGCgammaMass    = 0.;
Double_t kGCpi0Mass      = 0.1349766;
Double_t kGCetaMass      = 0.54751;

// define the width constraint used by KF particle.
Double_t kGCgammaWidth = 0.01;
Double_t kGCpi0Width   = 0.01;
Double_t kGCetaWidth   = 0.01;

// define the probability of track being an electron
Double_t kGCprobElectron = 0.000;

Double_t kGCminOpeningAngleGhostCut = 0.005;

Bool_t kGCRemovePileUp = kFALSE;

Bool_t kGCSelectV0AND = kFALSE;
Bool_t kGCUseMultiplicity = kFALSE;
Int_t  kGCUseMultiplicityBin=0;

Int_t kGCIsHeavyIon = 0;
Int_t kGCUseCentrality = 0;
Int_t kGCUseCentralityBin = 0;
Int_t kGCUseCorrectedTPCClsInfo = 0;


/** ---------------------------------- define pi0 dalitz cuts here ------------------------------------*/

Bool_t kGCRunStandalone    = kTRUE;
Bool_t kGCUseBayesPID      = kFALSE;
Bool_t kGCUseTrackIndexCut = kTRUE;
Bool_t kGCUsePsiPairCut    = kTRUE;
Bool_t kGCUseMassCut       = kFALSE;
Bool_t kGCUseGammaCut      = kFALSE;
Bool_t kGCReadMagFieldSign = kTRUE;
Bool_t kGCUseAliKF         = kFALSE;

Double_t kGCPsiPairCut                 = 0.45;
Double_t kGCDeltaPhiCutMin             = 0.;
Double_t kGCDeltaPhiCutMax             = 0.12;
Double_t kGCMassCutMin                 = 0.;
Double_t kGCMassCutMax                 = 0.1;
Double_t kGCNSigmaBelowElecTPCbethe    = -2.;
Double_t kGCNSigmaAboveElecTPCbethe    = 3.;
Double_t kGCNSigmaAbovePionTPCbethe    = 2.;
Double_t kGCNSigmaAboveKaonTPCbethe    = 2.;
Double_t kGCNSigmaAboveProtonTPCbethe  = 2.;

Int_t kGCTrkSelectionCriteria          = 1;     // kITSsaTrack=0, kGlobalTrack=1, kITSsaGlobalTrack=2
// NOTE: for details in the track cuts and defined histograms see AddGammaConvDalitz.C

/** ----------------------------------end define cuts here----------------------------------*/

/** -------------------------------- Phi/R Mapping ---------------------------------------*/
Int_t kGCnPhiIndex = 8;
Int_t kGCnRIndex   = 14;

Double_t kGCminRadius   = 0.;
Double_t kGCmaxRadius   = 200.;
Double_t kGCminPhi      = -TMath::Pi();
Double_t kGCmaxPhi      = TMath::Pi();
/** ------------------------------- end Phi/R Mapping ------------------------------------*/

Bool_t kGCdoOwnXYZCalculation = kFALSE;
Bool_t kGCdoConstructGamma = kFALSE;

/** -------------AOD stuff ---------------------------------------------------------------*/
TString kGCDeltaAODFilename = "AliAODGammaConversion.root";  //If empty, writes to standard common aod file.
Bool_t kGCWriteAOD =kTRUE;    // Turn on AOD 
Bool_t kGCForceAOD = kFALSE;  // Call AliAnalysisManager::SetFillAOD(kTRUE) every ESD event. 

/** ------------------- define which histograms to plot here --------------------------------*/
/**   NB: to change the bin numbers, see below the histogram flags                           */

Bool_t kGCplotEventQuality              = kTRUE;


// MC GAMMA DECAY PLOTS
Bool_t kGCplotMCGammaFromDecay          = kTRUE;
// NEUTRAL MESON PLOTS
Bool_t kGCplotMCConversionR             = kTRUE;
Bool_t kGCplotMCConversionZR            = kTRUE;
Bool_t kGCplotMCConversionXY            = kTRUE;
Bool_t kGCplotMCConversionOpeningAngle  = kTRUE;
Bool_t kGCplotMCConvGammaEAsymmetryP    = kTRUE;
Bool_t kGCplotMCConvGammaPAsymmetryP    = kTRUE;


Bool_t kGCplotMCEEnergy  = kFALSE;
Bool_t kGCplotMCEPt      = kTRUE;
Bool_t kGCplotMCEEta     = kTRUE;
Bool_t kGCplotMCEPhi     = kTRUE;
Bool_t kGCplotMCENTPCClusters = kTRUE;
Bool_t kGCplotMCENITSClusters = kTRUE;

Bool_t kGCplotMCPEnergy  = kFALSE;
Bool_t kGCplotMCPPt      = kTRUE;
Bool_t kGCplotMCPEta     = kTRUE;
Bool_t kGCplotMCPPhi     = kTRUE;
Bool_t kGCplotMCPNTPCClusters = kTRUE;
Bool_t kGCplotMCPNITSClusters = kTRUE;

Bool_t kGCplotMCallGammaEnergy = kFALSE;
Bool_t kGCplotMCallGammaPt     = kTRUE;
Bool_t kGCplotMCallGammaEta    = kTRUE;
Bool_t kGCplotMCallGammaPhi    = kTRUE;
Bool_t kGCplotMCallGammaRapid  = kTRUE;


Bool_t kGCplotMCConvGammaEnergy  = kFALSE;
Bool_t kGCplotMCConvGammaPt      = kTRUE;
Bool_t kGCplotMCConvGammaEta     = kTRUE;
Bool_t kGCplotMCConvGammaPhi     = kTRUE;
Bool_t kGCplotMCConvGammaRapid   = kTRUE;
Bool_t kGCplotMCConvGammaPtvsEta = kTRUE;

Bool_t kGCplotMCallDirectGammaEnergy  = kFALSE;
Bool_t kGCplotMCallDirectGammaPt      = kTRUE;
Bool_t kGCplotMCallDirectGammaEta     = kTRUE;
Bool_t kGCplotMCallDirectGammaPhi     = kTRUE;
Bool_t kGCplotMCallDirectGammaRapid   = kTRUE;

Bool_t kGCplotMCConvDirectGammaEnergy  = kFALSE;
Bool_t kGCplotMCConvDirectGammaPt      = kTRUE;
Bool_t kGCplotMCConvDirectGammaEta     = kTRUE;
Bool_t kGCplotMCConvDirectGammaPhi     = kTRUE;
Bool_t kGCplotMCConvDirectGammaRapid   = kTRUE;

Bool_t kGCplotMCMotherEta					= kTRUE;
Bool_t kGCplotMCMotherRapid                                = kTRUE;
Bool_t kGCplotMCMotherPhi					= kTRUE;
Bool_t kGCplotMCMotherPt					= kTRUE;
Bool_t kGCplotMCMotherEnergy				= kFALSE;
Bool_t kGCplotMCMotherMass					= kTRUE;
Bool_t kGCplotMCMotherOpeningAngle				= kTRUE;
Bool_t kGCplotMCMotherR					= kTRUE;
Bool_t kGCplotMCMotherZR					= kFALSE;
Bool_t kGCplotMCMotherXY	       				= kFALSE;
Bool_t kGCplotMCMotherPtvsEtaWithinAcceptance              = kTRUE;
Bool_t kGCplotMCMotherPtvsRapidWithinAcceptance            = kTRUE;
Bool_t kGCplotMCMotherPtvsEtaConvGammaWithinAcceptance     = kTRUE;
Bool_t kGCplotMCMotherPtvsRapidConvGammaWithinAcceptance   = kTRUE;
Bool_t kGCplotMCMotherSpectra				= kTRUE;

Bool_t kGCplotMCPhysicalPrimaryChargedPt                = kTRUE;

Bool_t kGCplotMCPi0Eta			        	= kTRUE;
Bool_t kGCplotMCPi0Rapid                                   = kTRUE;
Bool_t kGCplotMCPi0PtvsRapid                               = kTRUE;
Bool_t kGCplotMCPi0Phi                                     = kTRUE;
Bool_t kGCplotMCPi0Pt                                      = kTRUE;
Bool_t kGCplotMCPi0PtFiducial                              = kTRUE;
Bool_t kGCplotMCPi0PtWithinAcceptanceFiducial              = kTRUE;
Bool_t kGCplotMCPi0PtConvGammaWithinAcceptanceFiducial     = kTRUE;
Bool_t kGCplotMCPi0OpeningPtConvGammaWithinAcceptance      = kTRUE;
Bool_t kGCplotMCPi0PtGammaPtConvGammaWithinAcceptance = kTRUE;
Bool_t kGCplotMCPi0Energy                                  = kFALSE;
Bool_t kGCplotMCPi0Mass                                    = kTRUE;
Bool_t kGCplotMCPi0Alpha                                   = kTRUE;
Bool_t kGCplotMCEtaAlpha                                   = kTRUE;
Bool_t kGCplotMCPi0OpeningAngle                            = kTRUE;
Bool_t kGCplotMCPi0R                                       = kTRUE;
Bool_t kGCplotMCPi0ZR                                      = kFALSE;
Bool_t kGCplotMCPi0XY                                      = kFALSE;
Bool_t kGCplotMCPi0PtvsEtaWithinAcceptance                 = kTRUE;
Bool_t kGCplotMCPi0PtvsRapidWithinAcceptance               = kTRUE;
Bool_t kGCplotMCPi0PtvsEtaConvGammaWithinAcceptance        = kTRUE;
Bool_t kGCplotMCPi0PtvsRapidConvGammaWithinAcceptance      = kTRUE;
Bool_t kGCplotMCPi0ZRConvGammaWithinAcceptance		= kTRUE;

Bool_t kGCplotMCPi0SecondaryEta                                = kTRUE;
Bool_t kGCplotMCPi0SecondaryRapid                              = kTRUE;
Bool_t kGCplotMCPi0SecondaryPhi                                = kTRUE;
Bool_t kGCplotMCPi0SecondaryPt                                 = kTRUE;
Bool_t kGCplotMCPi0SecondaryEnergy                             = kFALSE;
Bool_t kGCplotMCPi0SecondaryMass                               = kTRUE;
Bool_t kGCplotMCPi0SecondaryOpeningAngle                       = kTRUE;
Bool_t kGCplotMCPi0SecondaryR                                  = kTRUE;
Bool_t kGCplotMCPi0SecondaryZR                                 = kFALSE;
Bool_t kGCplotMCPi0SecondaryXY                                 = kFALSE;
Bool_t kGCplotMCPi0SecondaryPtvsEtaWithinAcceptance            = kTRUE;
Bool_t kGCplotMCPi0SecondaryPtvsRapidWithinAcceptance          = kTRUE;
Bool_t kGCplotMCPi0SecondaryPtvsEtaConvGammaWithinAcceptance   = kTRUE;
Bool_t kGCplotMCPi0SecondaryPtvsRapidConvGammaWithinAcceptance = kTRUE;

Bool_t kGCplotMCEtaEta                                = kTRUE;
Bool_t kGCplotMCEtaRapid                              = kTRUE;
Bool_t kGCplotMCEtaPtvsRapid                          = kTRUE;
Bool_t kGCplotMCEtaPhi                                = kTRUE;
Bool_t kGCplotMCEtaPt                                 = kTRUE;
Bool_t kGCplotMCEtaEnergy                             = kFALSE;
Bool_t kGCplotMCEtaMass                               = kTRUE;
Bool_t kGCplotMCEtaOpeningAngleGamma                  = kTRUE;
Bool_t kGCplotMCEtaR                                  = kTRUE;
Bool_t kGCplotMCEtaZR                                 = kFALSE;
Bool_t kGCplotMCEtaXY                                 = kFALSE;
Bool_t kGCplotMCEtaPtvsEtaWithinAcceptance		   = kTRUE;
Bool_t kGCplotMCEtaPtvsRapidWithinAcceptance	   = kTRUE;
Bool_t kGCplotMCEtaPtvsEtaConvGammaWithinAcceptance   = kTRUE;
Bool_t kGCplotMCEtaPtvsRapidConvGammaWithinAcceptance = kTRUE;
Bool_t kGCplotMCEtaOpeningPtConvGammaWithinAcceptance = kTRUE;
Bool_t kGCplotMCEtaPtGammaPtConvGammaWithinAcceptance = kTRUE;
Bool_t kGCplotMCEtaZRConvGammaWithinAcceptance = kTRUE;

// Histograms from esd tracks
Bool_t kGCplotESDConversionR                   = kTRUE;
Bool_t kGCplotESDConversionZR                  = kTRUE;
Bool_t kGCplotESDConversionXY                  = kTRUE;
Bool_t kGCplotESDConversionXYBeamPipe          = kTRUE;
Bool_t kGCplotESDConversionRPhiBeamPipe        = kTRUE;
Bool_t kGCplotESDConversionOpeningAngle        = kTRUE;
Bool_t kGCplotESDConvGammaCosPointingAngle     = kTRUE;
Bool_t kGCplotESDConvGammaDcaDaugthers         = kTRUE;
Bool_t kGCplotESDConvGammaNormDcaDistDaugthers = kTRUE;
Bool_t kGCplotESDConvGammaLikelihoodAP         = kTRUE;
Bool_t kGCplotESDConvGammaEAsymmetryP         = kTRUE;
Bool_t kGCplotESDConvGammaPAsymmetryP         = kTRUE;
Bool_t kGCplotESDConvGammaEdEdxP         = kTRUE;
Bool_t kGCplotESDConvGammaPdEdxP         = kTRUE;
Bool_t kGCplotESDConvGammaEeProbP         = kTRUE;
Bool_t kGCplotESDConvGammaPeProbP         = kTRUE;
Bool_t kGCplotESDConvGammaEmupiProbP         = kTRUE;
Bool_t kGCplotESDConvGammaPmupiProbP         = kTRUE;

Bool_t kGCplotESDConvGammaQtAlfa         = kTRUE;


Bool_t kGCplotESDEEnergy = kFALSE;
Bool_t kGCplotESDEPt     = kTRUE;
Bool_t kGCplotESDEEta    = kTRUE;
Bool_t kGCplotESDEPhi    = kTRUE;
Bool_t kGCplotESDENTPCClusters = kTRUE;
Bool_t kGCplotESDENITSClusters = kTRUE;
Bool_t kGCplotESDENTPCClustersToFP = kTRUE;
Bool_t kGCplotESDENTPCClustersToFR = kTRUE;
Bool_t kGCplotESDETPCchi2 = kTRUE;

Bool_t kGCplotESDPEnergy = kFALSE;
Bool_t kGCplotESDPPt     = kTRUE;
Bool_t kGCplotESDPEta    = kTRUE;
Bool_t kGCplotESDPPhi    = kTRUE;
Bool_t kGCplotESDPNTPCClusters = kTRUE; 
Bool_t kGCplotESDPNITSClusters = kTRUE;
Bool_t kGCplotESDPNTPCClustersToFP = kTRUE;
Bool_t kGCplotESDPNTPCClustersToFR = kTRUE;
Bool_t kGCplotESDPTPCchi2 = kTRUE;

Bool_t kGCplotESDConvGammaEnergy = kFALSE;
Bool_t kGCplotESDConvGammaPt     = kTRUE;
Bool_t kGCplotESDConvGammaEta    = kTRUE;
Bool_t kGCplotESDConvGammaPhi    = kTRUE;
Bool_t kGCplotESDConvGammaMass   = kTRUE;
Bool_t kGCplotESDConvGammaWidth  = kTRUE;
Bool_t kGCplotESDConvGammaChi2   = kTRUE;
Bool_t kGCplotESDConvGammaNDF    = kTRUE;
Bool_t kGCplotESDConvGammaRapid  = kTRUE;
Bool_t kGCplotESDConvGammaPtvsEta = kTRUE;
Bool_t kGCplotESDConvGammaPtvsChi2 = kTRUE;
Bool_t kGCplotESDConvGammaEtavsChi2 = kTRUE;
Bool_t kGCplotESDConvGammaPtvsTOF = kTRUE;


Bool_t kGCplotESDTrueDalitzContaminationR    = kTRUE;
Bool_t kGCplotESDTruePi0DalitzContaminationR = kTRUE;
Bool_t kGCplotESDTrueEtaDalitzContaminationR = kTRUE;
Bool_t kGCplotESDTrueCombinatorialContaminationR = kTRUE;
Bool_t kGCplotESDTrueCombinatorialElecContaminationR = kTRUE;
Bool_t kGCplotESDTrueHadronicContaminationR  = kTRUE;
Bool_t kGCplotESDTrueCombinatorialContaminationPt = kTRUE;
Bool_t kGCplotESDTrueCombinatorialElecContaminationPt = kTRUE;
Bool_t kGCplotESDTrueHadronicContaminationPt  = kTRUE;


Bool_t kGCplotESDTrueBackground              = kTRUE;
Bool_t kGCplotESDTrueConvGammaEnergy         = kFALSE;
Bool_t kGCplotESDTrueConvGammaPt             = kTRUE;
Bool_t kGCplotESDTrueConvGammaEta            = kTRUE;
Bool_t kGCplotESDTrueConvGammaPhi            = kTRUE;
Bool_t kGCplotESDTrueConvGammaMass           = kTRUE;
Bool_t kGCplotESDTrueConvGammaWidth          = kTRUE;
Bool_t kGCplotESDTrueConvGammaChi2           = kTRUE;
Bool_t kGCplotESDTrueConvGammaNDF            = kTRUE;
Bool_t kGCplotESDTrueConvGammaRapid          = kTRUE;
Bool_t kGCplotESDTrueConvGammaPtvsEta        = kTRUE;
Bool_t kGCplotESDTrueConversionR             = kTRUE;
Bool_t kGCplotESDTrueConversionZR            = kFALSE;
Bool_t kGCplotESDTrueConversionXY            = kFALSE;
Bool_t kGCplotESDTrueConversionOpeningAngle  = kTRUE;
Bool_t kGCplotESDTrueConvGammaCosPointingAngle     = kTRUE;
Bool_t kGCplotESDTrueConvGammaDcaDaugthers         = kTRUE;
Bool_t kGCplotESDTrueConvGammaNormDcaDistDaugthers = kTRUE;
Bool_t kGCplotESDTrueConvGammaLikelihoodAP         = kTRUE;
Bool_t kGCplotESDTrueConvGammaEAsymmetryP         = kTRUE;
Bool_t kGCplotESDTrueConvGammaPAsymmetryP         = kTRUE;
Bool_t kGCplotESDTrueConvGammaEdEdxP         = kTRUE;
Bool_t kGCplotESDTrueConvGammaPdEdxP         = kTRUE;
Bool_t kGCplotESDTrueConvGammaQtAlfa         = kTRUE;

Bool_t kGCplotESDTrueConvGammaPtvsChi2       = kTRUE;
Bool_t kGCplotESDTrueConvGammaEtavsChi2      = kTRUE;
Bool_t kGCplotESDTrueConvGammaMCPtEta        = kTRUE;
Bool_t kGCplotESDTrueConversionMCZR          = kFALSE;
Bool_t kGCplotESDTrueConversionMCXY          = kFALSE;

Bool_t kGCplotESDNoCutAllV0Pt                = kTRUE;
Bool_t kGCplotESDNoCutConvGammaEnergy         = kFALSE;
Bool_t kGCplotESDNoCutConvGammaPt             = kTRUE;
Bool_t kGCplotESDNoCutConvGammaEta            = kTRUE;
Bool_t kGCplotESDNoCutConvGammaPhi            = kTRUE;
Bool_t kGCplotESDNoCutConvGammaMass           = kTRUE;
Bool_t kGCplotESDNoCutConvGammaWidth          = kTRUE;
Bool_t kGCplotESDNoCutConvGammaChi2           = kTRUE;
Bool_t kGCplotESDNoCutConvGammaNDF            = kTRUE;
Bool_t kGCplotESDNoCutConvGammaRapid          = kTRUE;
Bool_t kGCplotESDNoCutConvGammaPtvsEta        = kTRUE;
Bool_t kGCplotESDNoCutConversionR             = kTRUE;
Bool_t kGCplotESDNoCutConversionZR            = kFALSE;
Bool_t kGCplotESDNoCutConversionXY            = kFALSE;
Bool_t kGCplotESDNoCutConversionOpeningAngle  = kTRUE;
Bool_t kGCplotESDNoCutConvGammaCosPointingAngle     = kTRUE;
Bool_t kGCplotESDNoCutConvGammaDcaDaugthers         = kTRUE;
Bool_t kGCplotESDNoCutConvGammaNormDcaDistDaugthers = kTRUE;
Bool_t kGCplotESDNoCutConvGammaLikelihoodAP         = kTRUE;

Bool_t kGCplotESDNoCutConvGammaEAsymmetryP         = kTRUE;
Bool_t kGCplotESDNoCutConvGammaPAsymmetryP         = kTRUE;
Bool_t kGCplotESDNoCutConvGammaEdEdxP         = kTRUE;
Bool_t kGCplotESDNoCutConvGammaPdEdxP         = kTRUE;
Bool_t kGCplotESDNoCutConvGammaPtvsChi2       = kTRUE;
Bool_t kGCplotESDNoCutConvGammaEtavsChi2      = kTRUE;
Bool_t kGCplotESDNoCutConvGammaMCPtEta        = kTRUE;
Bool_t kGCplotESDNoCutConversionMCZR          = kFALSE;
Bool_t kGCplotESDNoCutConversionMCXY          = kFALSE;

Bool_t kGCplotESDMotherChi2 = kTRUE;
Bool_t kGCplotESDMotherOpeningAngleGamma = kTRUE;
Bool_t kGCplotESDMotherEnergy            = kFALSE;
Bool_t kGCplotESDMotherPt                = kFALSE;
Bool_t kGCplotESDMotherEta               = kTRUE;
Bool_t kGCplotESDMotherPhi               = kFALSE;
Bool_t kGCplotESDMotherMass              = kFALSE;
Bool_t kGCplotESDMotherR                 = kFALSE;
Bool_t kGCplotESDMotherZR                = kFALSE;
Bool_t kGCplotESDMotherXY                = kFALSE;
Bool_t kGCplotESDMotherRapid             = kTRUE;

Bool_t kGCplotESDBackgroundOpeningAngleGamma = kTRUE;
Bool_t kGCplotESDBackgroundEnergy            = kFALSE;
Bool_t kGCplotESDBackgroundPt                = kFALSE;
Bool_t kGCplotESDBackgroundEta               = kFALSE;
Bool_t kGCplotESDBackgroundPhi               = kFALSE;
Bool_t kGCplotESDBackgroundMass              = kFALSE;
Bool_t kGCplotESDBackgroundR                 = kFALSE;
Bool_t kGCplotESDBackgroundZR                = kFALSE;
Bool_t kGCplotESDBackgroundXY                = kFALSE;
Bool_t kGCplotESDBackgroundRapid             = kFALSE;

Bool_t kGCplotMapping = kTRUE;       

Bool_t kGCplotResolutiondPt = kTRUE;
Bool_t kGCplotResolutiondR  = kTRUE;
Bool_t kGCplotResolutiondZ  = kTRUE;

Bool_t kGCplotResolutiondRAbs  = kTRUE;
Bool_t kGCplotResolutiondZAbs  = kTRUE;
Bool_t kGCplotResolutiondPhiAbs  = kTRUE;

Bool_t kGCplotResolutiondRdPt = kTRUE;

Bool_t kGCplotResolutionMCPt = kTRUE;
Bool_t kGCplotResolutionMCR  = kTRUE;
Bool_t kGCplotResolutionMCZ  = kTRUE;

Bool_t kGCplotResolutionESDPt = kTRUE;
Bool_t kGCplotResolutionESDR  = kTRUE;
Bool_t kGCplotResolutionESDZ  = kTRUE;

Bool_t kGCplotResolutionPtdPt = kTRUE;

Bool_t kGCplotESDNumberOfV0s          = kTRUE;
Bool_t kGCplotESDNumberOfSurvivingV0s = kTRUE;
Bool_t kGCplotESDNumberOfContributorsVtx = kTRUE;
Bool_t kGCplotESDNumberOfGoodESDTracks = kTRUE;

//  debug histograms
Bool_t kGCplotESDCutGetOnFly      = kTRUE;
Bool_t kGCplotESDCutNContributors = kTRUE;
Bool_t kGCplotESDCutLikeSign      = kTRUE;
Bool_t kGCplotESDCutRefit         = kTRUE;
Bool_t kGCplotESDCutKink          = kTRUE;
Bool_t kGCplotESDCutPIDProb       = kTRUE;
Bool_t kGCplotESDCutdedxSigmaElectronLine=kTRUE;
Bool_t kGCplotESDCutdedxSigmaPionLine=kTRUE;
Bool_t kGCplotESDCutPionRejectionLowP  =kTRUE;
Bool_t kGCplotESDCutProtonRejectionLowP=kTRUE;
Bool_t kGCplotESDCutKaonRejectionLowP  =kTRUE;
Bool_t kGCplotESDCutQtGammaSelection=kTRUE;
Bool_t kGCplotESDCutR             = kTRUE;
Bool_t kGCplotESDCutMinR          = kTRUE;
Bool_t kGCplotESDCutLine          = kTRUE;
Bool_t kGCplotESDCutZ             = kTRUE;
Bool_t kGCplotESDCutMinClsTPC     = kTRUE;
Bool_t kGCplotESDCutMinClsTPCToF  = kTRUE;
Bool_t kGCplotESDGoodV0s          = kTRUE;
Bool_t kGCplotESDAllV0s           = kTRUE;
Bool_t kGCplotESDAllV0sCurrentFinder = kTRUE;
Bool_t kGCplotESDAllV0sCurrentFinderQtAlfa = kTRUE;

Bool_t kGCplotESDCutNDF           = kTRUE;
Bool_t kGCplotESDCutChi2          = kTRUE;
Bool_t kGCplotESDCutEta           = kTRUE;
Bool_t kGCplotESDCutPt            = kTRUE;
Bool_t kGCplotESDCutSinglePt      = kTRUE;
Bool_t kGCplotESDTrueConvGammaTrackLength =kFALSE;
Bool_t kGCplotESDTrueConvGammaTrackLengthVSInvMass =kFALSE;

Bool_t kGCplotPi0Spectra = kTRUE;
Bool_t kGCplotEtaSpectra = kTRUE;
Bool_t kGCplotOmegaSpectra = kTRUE;

/////////////Chi_c Analysis//////////////////////////
Bool_t kGCplotStatsElectrons                                  = kTRUE;
Bool_t kGCplotRecENegJPsiPtDiff                               = kTRUE;
Bool_t kGCplotRecEPosJPsiPtDiff                               = kTRUE;
Bool_t kGCplotRecEPosENegR                                    = kTRUE;
Bool_t kGCplotRecEPosENegEta                                  = kTRUE;
Bool_t kGCplotESDInvMassePluseMinus                           = kTRUE;
Bool_t kGCplotESDInvMassGammaePluseMinusChiC                  = kTRUE;
Bool_t kGCplotESDInvMassGammaePluseMinusPi0                   = kTRUE;
Bool_t kGCplotESDElectronPosNegPt                             = kTRUE;
Bool_t kGCplotESDElectronPosNegEta                            = kTRUE;
Bool_t kGCplotESDElectronPosNegAngle                          = kTRUE;
Bool_t kGCplotMCElectronPosNegPt                              = kTRUE;
Bool_t kGCplotMCElectronPosNegEta                             = kTRUE;
Bool_t kGCplotMCElectronPosNegJPsiAngle                       = kTRUE;
Bool_t kGCplotESDElectronPosNegPi0Angle                       = kTRUE;
Bool_t kGCplotMCElectronPosNegPi0Angle                        = kTRUE;
Bool_t kGCplotTableElectrons                                  = kTRUE;
Bool_t kGCplotESDEPosBackground                               = kTRUE;
Bool_t kGCplotESDENegBackground                               = kTRUE;
Bool_t kGCplotESDEPosENegBackground                           = kTRUE;
Bool_t kGCplotESDEPosENegBackgroundCut                        = kTRUE;
Bool_t kGCplotESDePoseNegAngle                                = kTRUE;
Bool_t kGCplotESDEPosENegGammaBackgroundMX                    = kTRUE;
Bool_t kGCplotMCLabels                                        = kTRUE;
///////////////////////////////////////////////////////////////////

//---------------- Gamma Jet analysis ----------------------------
Bool_t kGCplotdPhiHdrGam            = kTRUE;
Bool_t kGCplotdPhiHdrGamIsolated    = kTRUE;
Bool_t kGCplotMinimumIsoDistance    = kTRUE;
Bool_t kGCplotFFzHdrGam             = kTRUE;
Bool_t kGCplotImbalanceHdrGam       = kTRUE;
//----------------------------------------------------------------


/** ----------------- end define which histograms to plot here -------------------------------*/



/** ----------- Define the binning for the different plot types here -------------------------*/
// Number of V0s-plot
Int_t kGCnXBinsNV0        = 1000;
Double_t kGCfirstXBinNV0  = -0.5;
Double_t kGClastXBinNV0   =  999.5;


// Number of ESD track-Plot
Int_t kGCnXBinsESDtrk      = 10000;
Double_t kGCfirstXBinESDtrk= -0.5;
Double_t kGClastXBinESDtrk = 9999.5;


//EventQuality-plot
Int_t kGCnXBinsEvtQ= 9;
Double_t kGCfirstXBinEvtQ=-1.5;
Double_t kGClastXBinEvtQ=7.5;

//R-plots
Int_t kGCnXBinsR = 400;
Double_t kGCfirstXBinR = 0.;
Double_t kGClastXBinR = 200.;

//ZR-plots
Int_t kGCnXBinsZR = 1000;
Double_t kGCfirstXBinZR = -250.;
Double_t kGClastXBinZR = 250.;
Int_t kGCnYBinsZR = 400;
Double_t kGCfirstYBinZR = 0.;
Double_t kGClastYBinZR = 200.;

//XY-plots
Int_t kGCnXBinsXY = 800;
Double_t kGCfirstXBinXY = -200.;
Double_t kGClastXBinXY = 200.;
Int_t kGCnYBinsXY = 800;
Double_t kGCfirstYBinXY = -200.;
Double_t kGClastYBinXY = 200.;

//XY-plots-BeamPipe
Int_t kGCnXBinsXYBP = 200;
Double_t kGCfirstXBinXYBP = -10.;
Double_t kGClastXBinXYBP = 10.;
Int_t kGCnYBinsXYBP = 200;
Double_t kGCfirstYBinXYBP = -10.;
Double_t kGClastYBinXYBP = 10.;

//Rphi-plots-BeamPipe
Int_t kGCnXBinsRPhiBP = 200;
Double_t kGCfirstXBinRPhiBP = -TMath::Pi();
Double_t kGClastXBinRPhiBP = TMath::Pi();
Int_t kGCnYBinsRPhiBP = 200;
Double_t kGCfirstYBinRPhiBP = 0.;
Double_t kGClastYBinRPhiBP = 10.;



//OpenAngle-plots
Int_t kGCnXBinsOpeningAngle = 400;
Double_t kGCfirstXBinOpeningAngle = 0.;
Double_t kGClastXBinOpeningAngle = TMath::Pi();

//CosPointingAngle-plots
Int_t kGCnXBinsCosPointingAngle = 400;
Double_t kGCfirstXBinCosPointingAngle = 0.99;
Double_t kGClastXBinCosPointingAngle = 1.01;

//DCA Daugthers-plots
Int_t kGCnXBinsDcaDaughters = 400;
Double_t kGCfirstXBinDcaDaughters= 0.;
Double_t kGClastXBinDcaDaughters = 5.;

//Norm DCA dist Daugthers-plots
Int_t kGCnXBinsNormDcaDistDaughters = 400;
Double_t kGCfirstXBinNormDcaDistDaughters= 0.;
Double_t kGClastXBinNormDcaDistDaughters = 10.;

//LikelihoodAP Plots
Int_t kGCnXBinsLikelihoodAP = 400;
Double_t kGCfirstXBinLikelihoodAP= 0.;
Double_t kGClastXBinLikelihoodAP = 2.;


//Energy-plots
Int_t kGCnXBinsEnergy = 200;
Double_t kGCfirstXBinEnergy = 0.;
Double_t kGClastXBinEnergy = 50.;

//P-plots
Int_t kGCnXBinsP = 200;
Double_t kGCfirstXBinP = 0.05;
Double_t kGClastXBinP = 50.;

//dEdx-plots
Int_t kGCnYBinsdEdx = 200;
Double_t kGCfirstYBindEdx = 0.;
Double_t kGClastYBindEdx = 200.;

//EProb-plots
Int_t kGCnYBinsEProb = 200;
Double_t kGCfirstYBinEProb = 0.;
Double_t kGClastYBinEProb = 1.;

//Qt-plots
Int_t kGCnYBinsQt = 250;
Double_t kGCfirstYBinQt = 0.;
Double_t kGClastYBinQt = 0.25;



//Asymmetry-plots
Int_t kGCnYBinsAsymmetry = 200;
Double_t kGCfirstYBinAsymmetry = 0.;
Double_t kGClastYBinAsymmetry = 1.;


//Pt-plots
Int_t kGCnXBinsPt = 500;
Double_t kGCfirstXBinPt = 0.;
Double_t kGClastXBinPt = 50.;

//TOF-plots RRnewTOF start /////////////////////////////////////////
Int_t kGCnXBinsTOFsignal = 500;
Double_t kGCfirstXBinTOFsignal = -10000.;
Double_t kGClastXBinTOFsignal = 40000.;
// end RRnewTOF ////////////////////////////////////////////////////

//Eta-plots
Int_t kGCnXBinsEta = 40;
Double_t kGCfirstXBinEta = -2.;
Double_t kGClastXBinEta = 2.;

//Rapidity
Int_t kGCnXBinsRapid = 200;
Double_t kGCfirstXBinRapid = -10.;
Double_t kGClastXBinRapid = 10.;

//Phi-plots
Int_t kGCnXBinsPhi = 72;
Double_t kGCfirstXBinPhi = -TMath::Pi();
Double_t kGClastXBinPhi = TMath::Pi();

//nTPCCluster-plots
Int_t kGCnXBinsNTPCClusters = 201;
Double_t kGCfirstXBinNTPCClusters = -0.5;
Double_t kGClastXBinNTPCClusters = 200.5;

//nITSCluster-plots
Int_t kGCnXBinsNITSClusters = 7;
Double_t kGCfirstXBinNITSClusters = -0.5;
Double_t kGClastXBinNITSClusters = 6.5;

//TPCcluster to Findable-plots
Int_t kGCnYBinsClsToF = 200;
Double_t kGCfirstYBinClsToF = 0.;
Double_t kGClastYBinClsToF = 2.0;

//TPCchi2 -plots
Int_t kGCnXBinsTPCchi2 = 100;
Double_t kGCfirstXBinTPCchi2 = 0.;
Double_t kGClastXBinTPCchi2 = 10;


//Mapping-plots
Int_t kGCnXBinsMapping = 800;
Double_t kGCfirstXBinMapping = -100.;
Double_t kGClastXBinMapping = 100.;
Int_t kGCnYBinsMapping = 40;
Double_t kGCfirstYBinMapping = -2;
Double_t kGClastYBinMapping = 2;

//ResolutionPlots
//RESdPt
Int_t kGCnXBinsResdPt=500;
Int_t kGCfirstXBinResdPt= 0;
Int_t kGClastXBinResdPt=100;
Int_t kGCnYBinsResdPt=200;
Int_t kGCfirstYBinResdPt= -10;
Int_t kGClastYBinResdPt=10;

//RESdR
Int_t kGCnXBinsResdR=400;
Int_t kGCfirstXBinResdR= 0;
Int_t kGClastXBinResdR=200;
Int_t kGCnYBinsResdR=100;
Int_t kGCfirstYBinResdR= -25;
Int_t kGClastYBinResdR=25;

//RESdZ
Int_t kGCnXBinsResdZ=80;
Int_t kGCfirstXBinResdZ= -20;
Int_t kGClastXBinResdZ=20;
Int_t kGCnYBinsResdZ=80;
Int_t kGCfirstYBinResdZ= -20;
Int_t kGClastYBinResdZ=20;

//RESdRdPt
Int_t kGCnYBinsResdRdPt=400;
Int_t kGCfirstYBinResdRdPt= -10;
Int_t kGClastYBinResdRdPt=10;

//RESMCPt
Int_t kGCnXBinsResPt=500;
Int_t kGCfirstXBinResPt= 0;
Int_t kGClastXBinResPt=100;

//RESMCR
Int_t kGCnXBinsResR=500;
Int_t kGCfirstXBinResR= 0;
Int_t kGClastXBinResR=250;

//RESMCZ
Int_t kGCnXBinsResZ=500;
Int_t kGCfirstXBinResZ= 0;
Int_t kGClastXBinResZ=250;

//GammaMass-plots
Int_t kGCnXBinsGammaMass = 4000;
Double_t kGCfirstXBinGammaMass = 0.;
Double_t kGClastXBinGammaMass = 1.;

//Pi0Mass-plots
Int_t kGCnXBinsPi0Mass = 1000;
Double_t kGCfirstXBinPi0Mass = 0.;
Double_t kGClastXBinPi0Mass = 1.;
Double_t kGCfirstXBinPi0Alpha = 0.;
Double_t kGClastXBinPi0Alpha = 1.;


//EtaMass-plots
Int_t kGCnXBinsEtaMass = 1000;
Double_t kGCfirstXBinEtaMass = 0.;
Double_t kGClastXBinEtaMass = 1.;

//GammaWidth-plots
Int_t kGCnXBinsGammaWidth = 100;
Double_t kGCfirstXBinGammaWidth = 0.;
Double_t kGClastXBinGammaWidth = 1.;

//GammaChi2-plots
Int_t kGCnXBinsMesonChi2 = 200;
Int_t kGCnXBinsGammaChi2 = 100;
Double_t kGCfirstXBinGammaChi2 = 0;
Double_t kGClastXBinGammaChi2 = 200.;

//GammaNDF-plots
Int_t kGCnXBinsGammaNDF = 10;
Double_t kGCfirstXBinGammaNDF = 0.;
Double_t kGClastXBinGammaNDF = 10.;

//Spectra-plots
Int_t kGCnXBinsSpectra = 500;
Double_t kGCfirstXBinSpectra = 0.;
Double_t kGClastXBinSpectra = 1.;
Int_t kGCnYBinsSpectra = 250;
Double_t kGCfirstYBinSpectra = 0.;
Double_t kGClastYBinSpectra = 25.;

Double_t kGCfirstXBinAlphaG = -1.;
Double_t kGCfirstXBinAlpha = 0.;
Double_t kGClastXBinAlpha = 1.;

//track length plots
Int_t kGCnXBinsTrackLength = 1000;
Double_t kGCfirstXBinTrackLength = 0;
Double_t kGClastXBinTrackLength = 500;

/////////Pi0 Dalitz decay AnalysisTask ///////////////////////////////////

Int_t kGCnXBinsDalitzMass = 4000;
Double_t kGCfirstXBinDalitzMass = 0.;
Double_t kGClastXBinDalitzMass = 4.;

Int_t kGCnXBinsPi0DalitzMass = 4000;
Double_t kGCfirstXBinPi0DalitzMass = 0.;
Double_t kGClastXBinPi0DalitzMass = 4.;

/////////Chic_Analysis///////////////////////////////////
Int_t kGCnXBinsEPt = 1000;
Double_t kGCfirstXBinEPt = 0.;
Double_t kGClastXBinJPsiPt  = 10;

Int_t kGCnXBinsJPsiMass = 1000;
Double_t kGCfirstXBinJPsiMass = 0.;
Double_t kGClastXBinJPsiMass = 10.;

Int_t kGCnXBinsChicMass = 1000;
Double_t kGCfirstXBinChicMass = 0.;
Double_t kGClastXBinChicMass  = 10.;

Int_t kGCnXBinsPi0Mass  = 1000;
Double_t kGCfirstXBinPi0Mass = 0.;
Double_t kGClastXBinPi0Mass  = 1.;

Int_t kGCnXBinsEPosNegPt = 1000;
Double_t kGCfirstXBinEPosNegPt = 0.;
Double_t kGClastXBinEPosNegPt  = 10.;

Int_t kGCnXBinsEPosNegEta = 200;
Double_t kGCfirstXBinEPosNegEta = -1.2;
Double_t kGClastXBinEPosNegEta  = 1.2;

Int_t kGCnXBinsEPosNegAngle = 200;
Double_t kGCfirstXBinEPosNegAngle = 0.;
Double_t kGClastXBinEPosNegAngle = TMath::Pi();

Int_t kGCnXBinsEBackground = 1000;
Double_t kGCfirstXBinEBackground = 0.;
Double_t kGClastXBinEBackground  = 10.;

Int_t kGCnXBinsEBackgroundCut = 100;
Double_t kGCfirstXBinEBackgroundCut = 0.;
Double_t kGClastXBinEBackgroundCut  = 0.015.;

Int_t kGCnXBinsMCLabels = 10;
Double_t kGCfirstXBinMCLabels = 0.;
Double_t kGClastXBinMCLabels  = 10.;

Int_t kGCnElementsElectronTable = 19;

//18 elements
const char * kGCelectronTable[] = {
  "Num. Events",  "MC e+ J/Psi |\\eta|<0.9","MC e- J/Psi |\\eta|<0.9","MC e+ e+ from J/Psi |\\eta|<0.9",
  "ESDtracks",    "Kink Cut",
  "Vertex Cut","TRDOut","TRDrefit","TPCrefit",
  "ITSrefit","TRDout+TPC+TPC+ITS+nsigma>3 Pass","pid!=0","ESDElec","ESD e+ JPsi",
  "ESD e- JPsi","ESD e+ e- JPSI","MC: gamma < 1.2","e+,e- < 0.9 g <1.2"
};


// for Gamma Jet analysis
Int_t kGCnXBinsdphiHdrGam = 100;
Double_t kGCfirstXBindphiHdrGam = -TMath::PiOver2();
Double_t kGClastXBindphiHdrGam = 3*TMath::PiOver2();

Int_t kGCnXBinsMinimumIsoDistance = 100;
Double_t kGCfirstXBinMinimumIsoDistance = 0.;
Double_t kGClastXBinMinimumIsoDistance = TMath::PiOver2();

Int_t kGCnXBinsFFzHdrGam = 100;
Double_t kGCfirstXBinFFzHdrGam = 0.;
Double_t kGClastXBinFFzHdrGam = 5;

Int_t kGCnXBinsImbalanceHdrGam = 100;
Double_t kGCfirstXBinImbalanceHdrGam = -5.;
Double_t kGClastXBinImbalanceHdrGam = 5.;
////////////////////////////////////////////////////////


/** ---------- end Define the binning for the different plot types here ----------------------*/


/************************************************************************************************
 *                                                                                              *
 *                                                                                              *
 *                     EVERYTHING BELOW IS FOR DEVELOPERS ONLY                                  *
 *                                                                                              *
 *                                                                                              *
 ************************************************************************************************/
TString kGCoutputFileName = "histogramsGammaConversion";
TString kGCoutputFileAppendix = "";
TString kGCdataList = "";
Bool_t kGCwriteNtuple = kFALSE;
// WE DOO NOT NEED TO CHANGE THIS (kGCusePWG4PartCorr) ANYMORE SINCE IT IS TAKEN CARE OF AUTOMATICALLY NOW
Bool_t kGCusePWG4PartCorr = kTRUE;

/** Flag to enable running on train  */
Bool_t kGCrunOnTrain = kFALSE;
Bool_t kGCrunOnGsiTrain = kFALSE;

/** ------------------------------ Monte Carlo flag -----------------------------------------*/
Bool_t kGCdoMCTruth = kTRUE;
/** ---------------------------- end Monte Carlo flag ---------------------------------------*/

/** ------------------------------ Selecting trigger CINT1B -----------------------------------*/
Bool_t kGCtriggerCINT1B = kFALSE;
/** ---------------------------- end Monte Carlo flag ---------------------------------------*/

/** ------------------------- Choose KFParticle OR ESDTrack  --------------------------------*/
Bool_t kGCuseKFParticle = kTRUE;
Bool_t kGCuseESDTrack   = kFALSE;
/** ----------------------- end Choose KFParticle OR ESDTrack  -----------------------------*/

/**------------------------------Flag to apply dEdx cut base on sigmas to electron line----------*/
Bool_t kGCdodEdxSigmaCut= kTRUE;
/**------------------------------end Flag to apply NsigmadEdx cut ----------*/

Double_t kGCPIDnSigmaAboveElectronLine=5;
Double_t kGCPIDnSigmaBelowElectronLine=-3;
Double_t kGCPIDnSigmaAbovePionLine=0;
Double_t kGCPIDMinPnSigmaAbovePionLine=1.;
Double_t kGCPIDMaxPnSigmaAbovePionLine=3.;

Bool_t    kGCuseTOFpid = kFALSE; // RRnewTOF start //////////
Double_t  kGCtofPIDnSigmaBelowElectronLine=-100;
Double_t  kGCtofPIDnSigmaAboveElectronLine=100; // RRnewTOF end

/**------- Flag to apply rejection at LowP of Kaons, protons , pions------------*/
Bool_t kGCdoPionRejectionLowP=kTRUE;
Bool_t kGCdoKaonRejectionLowP=kTRUE;
Bool_t kGCdoProtonRejectionLowP=kTRUE;
/**---------------------------------------*/

Double_t kGCPIDnSigmaAtLowPAroundKaonLine=0;
Double_t kGCPIDnSigmaAtLowPAroundProtonLine=0;
Double_t kGCPIDnSigmaAtLowPAroundPionLine=0;

Double_t kGCPIDMinPKaonRejectionLowP=1.5;
Double_t kGCPIDMinPProtonRejectionLowP=2.;
Double_t kGCPIDMinPPionRejectionLowP=0.5;


Bool_t kGCdoQtGammaSelection=kTRUE;
Double_t kGCQtMax=100.;
Bool_t kGCdoHighPtQtGammaSelection=kFALSE; // RRnew
Double_t kGCHighPtQtMax=100.;		   // RRnew
Double_t kGCPtBorderForQt=100.;	  	   // RRnew



Bool_t scanArguments(TString arguments){
  
  Bool_t iResult = kTRUE;
	
  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  cout<<"Arguments received: "<<allArgs.Data()<<endl;
	
  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
		
    for(int i=0; i<pTokens->GetEntries() && iResult==kTRUE; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
			
      if(argument.IsNull()) continue;
      // -- deconvolute-time option
      if(argument.CompareTo("-data-list") == 0){
	if((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	kGCdataList = ((TObjString*)pTokens->At(i))->GetString();
	if(kGCdataList.IsNull()){
	  cout<<"-data-list is NULL"<<endl;
	  iResult=kFALSE;
	}
	else{
	  cout<<"Data list is set to: "<<kGCdataList<<endl;
	}
      }
      else if(argument.CompareTo("-output-file-name") == 0){
	if((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	kGCoutputFileName = ((TObjString*)pTokens->At(i))->GetString();
	if(kGCoutputFileName.IsNull()){
	  cout<<"-output-file-name is NULL"<<endl;
	  iResult=kFALSE;
	}
	else{
	  cout<<"Setting output file name to: "<<kGCoutputFileName<<endl;
	}
      }
      else if (argument.CompareTo("-bg-off") == 0){
	kGCcalculateBackground =kFALSE;
      }
      else if (argument.CompareTo("-bg-prob-off") == 0){
	kGCdoBGProbability = kFALSE;
      }
      else if (argument.CompareTo("-bg-prob-on") == 0){
	kGCdoBGProbability = kTRUE;
      }
      else if (argument.CompareTo("-bg-rotation-off") == 0){
	kGCUseRotationMethodInBG = kFALSE;
      }
      else if (argument.CompareTo("-use-v0-multiplicity") == 0){
	kGCUseTrackMultiplicityForBG = kFALSE;
      }
      else if (argument.CompareTo("-apply-chi2-cut") == 0){
	kGCApplyChi2Cut = kTRUE;
      }
      else if(argument.CompareTo("-set-number-of-rotations") == 0){
	if((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	kGCnumberOfRotationEventsForBG = ((TObjString*)pTokens->At(i))->GetString().Atoi();
      }
      else if(argument.CompareTo("-set-number-of-degrees") == 0){
	if((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	kGCnDegreeRotationPMForBG = ((TObjString*)pTokens->At(i))->GetString().Atoi();
      }
      else if (argument.CompareTo("-low-memory") == 0){
	kGCLowMemoryConsumption	= kTRUE;
      }
      else if (argument.CompareTo("-move-bg-vertex") == 0){
	kGCMoveParticlesAccordingToVertex = kTRUE;
      }
      else if (argument.CompareTo("-check-neutralmeson-pi0s") == 0){
	kGCdoNeutralMesonV0MCCheck=kTRUE;
      }
      else if (argument.CompareTo("-use-offline-finder") == 0){
	kGCUseOnFlyV0Finder = kFALSE;
      }
      else if (argument.CompareTo("-write-ntuple") == 0){
	cout<<"Writing ntuple to file."<<endl;
	kGCwriteNtuple = kTRUE;
      }
      else if (argument.CompareTo("-run-on-train") == 0){
	cout<<"Running on train"<<endl;
	//kGCWriteStandardAOD=kTRUE;
	kGCDeltaAODFilename="";
	kGCrunOnTrain = kTRUE;
      }
      else if (argument.CompareTo("-run-on-gsi-train") == 0){
	cout<<"Running on gsi train"<<endl;
	//kGCWriteStandardAOD=kFALSE;
	kGCrunOnTrain = kTRUE;
	kGCrunOnGsiTrain = kTRUE;
      }
      else if (argument.CompareTo("-run-jet") == 0){
	cout<<"Running jet analysis"<<endl;
	kGCrunJet = kTRUE;
      }
      else if (argument.CompareTo("-run-gamma-jet-task") == 0){
	cout<<"Running gamma jet task"<<endl;
	kGCRunGammaJetTask= kTRUE;
      }
      else if (argument.CompareTo("-run-neutralmeson") == 0){
	cout<<"Running neutral meson analysis"<<endl;
	kGCrunNeutralMeson = kTRUE;
      }
      else if (argument.CompareTo("-run-neutral-meson") == 0){
	cout<<"Running neutral meson analysis"<<endl;
	kGCrunNeutralMeson = kTRUE;
      }
      else if (argument.CompareTo("-run-omega-meson") == 0){
	cout<<"Running omega meson analysis"<<endl;
	kGCrunOmegaMeson = kTRUE;
      }
      else if (argument.CompareTo("-run-chic") == 0){
	cout<<"Running Chi_c analysis"<<endl;
	kGCrunChic = kTRUE;
      }
      else if (argument.CompareTo("-run-dalitz") == 0){
	cout<<"Running Dalitz analysis"<<endl;
	kGCrunDalitz = kTRUE;
      }
      else if (argument.CompareTo("-run-cf") == 0){
	cout<<"Running CF"<<endl;
	kGCrunCF = kTRUE;
      }
      else if (argument.CompareTo("-run-resolution") == 0){
	cout<<"Running Resolution"<<endl;
	kGCrunRES = kTRUE;
      }
      else if (argument.CompareTo("-run-recalculateV0") == 0){
	cout<<"Running RecalculateV0ForGamma"<<endl;
	kGCRecalculateV0ForGamma=kTRUE;
      }
      else if (argument.CompareTo("-jet-off") == 0){
	cout<<"Skipping jet analysis"<<endl;
	kGCrunJet = kFALSE;
      }
      else if (argument.CompareTo("-neutralmeson-off") == 0){
	cout<<"Skipping neutral meson analysis"<<endl;
	kGCrunNeutralMeson = kFALSE;
      }
      else if (argument.CompareTo("-neutral-meson-off") == 0){
	cout<<"Skipping neutral meson analysis"<<endl;
	kGCrunNeutralMeson = kFALSE;
      }
      else if (argument.CompareTo("-chic-off") == 0){
	cout<<"Skipping Chi_c analysis"<<endl;
	kGCrunChic = kFALSE;
      }
      else if (argument.CompareTo("-mc-off") == 0){
	cout<<"Switching off kGCdoMCTruth"<<endl;
	kGCdoMCTruth = kFALSE;
      }
      else if (argument.CompareTo("-trigger-CINT1B") == 0){
        cout<<"Selecting ONLY kGCtriggerCINT1B"<<endl;
        kGCtriggerCINT1B = kTRUE;
      }
      else if (argument.CompareTo("-use-own-xyz") == 0){
	cout<<"Switching on use own xyz calculation"<<endl;
	kGCdoOwnXYZCalculation = kTRUE;
      }
      else if (argument.CompareTo("-use-ConstructGamma") == 0){
	cout<<"Switching on use ConstructGamma and OFF own xyz calculation"<<endl;
	kGCdoOwnXYZCalculation = kFALSE;
	kGCdoConstructGamma = kTRUE;
      }
      else if (argument.CompareTo("-no-aod") == 0){
	cout<<"Turning off AOD"<<endl;
	kGCWriteAOD = kFALSE;
      }
      else if (argument.CompareTo("-force-aod") == 0){
	cout<<"Turning on FillAOD = kTRUE every event; If running in common train this should probably not be done!!!"<<endl;
	kGCForceAOD = kTRUE;
      }
      else if (argument.CompareTo("-standard-aod") == 0){
	cout<<"Writing to standard AOD, will only work on train"<<endl;
	kGCDeltaAODFilename = "";
      }
      else if(argument.CompareTo("-append-to-output-file") == 0){
	if((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	kGCoutputFileAppendix = TString("_")+((TObjString*)pTokens->At(i))->GetString();
	if(kGCoutputFileAppendix.IsNull()){
	  cout<<"-appending-to-output-file is NULL"<<endl;
	  iResult=kFALSE;
	}
	else{
	  cout<<"Appending to the output file: "<<kGCoutputFileAppendix<<endl;
	}
      }
      else if(argument.CompareTo("-set-cut-selection") == 0){
	if((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	kGCAnalysisCutSelectionId = ((TObjString*)pTokens->At(i))->GetString();
	cout<<"The analysis cut selection is set to: "<<kGCAnalysisCutSelectionId.Data()<<endl;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    cout<<"Missing parameter for argument "<< argument.Data()<<endl;
    iResult=kFALSE;
  }
  return iResult;
}

void SetVersionLibrary(){
  // Check if the file $ALICE_ROOT/PWG4/GammaConv/AliAnalysisTaskGammaConversion.cxx exists.
  // If yes, we set kGCusePWG4PartCorr to false since we have a newer version
  // If no, kGCusePWG4PartCorr is true.
	
  TString file = gSystem->Getenv("ALICE_ROOT");
  file+="/PWG4/PartCorr/AliAnalysisTaskGammaConversion.cxx";
	
  ifstream stream;
  stream.open(file.Data());
	
  if(!stream){
    kGCusePWG4PartCorr=kFALSE;
  }
  else{
    kGCusePWG4PartCorr=kTRUE;
  }
  stream.close();
}


AliAnalysisTaskGammaConversion* ConfigGammaConversion(TString arguments, AliAnalysisDataContainer *cin_esd=NULL){

						      
	
  
  if(!scanArguments(arguments)){
    break;
  }
  	
  SetVersionLibrary(); // checks if PWG4GammaConv or PWG4PartCorr is used
	
  if(cin_esd == NULL && kGCrunOnTrain == kTRUE){
    cout<<"Error: kGCrunOnTrain flag is set to true but the input AliAnalysisDataContainer is NULL"<<endl;
    cout<<"       you must also supply the AliAnalysisDataContainer as an argument"<<endl;
    return;
  }
	
  if(cin_esd != NULL && kGCrunOnTrain == kFALSE){
    cout<<"Error: kGCrunOnTrain flag is set to false but the input AliAnalysisDataContainer is not null"<<endl;
    cout<<"       add -run-on-train to the arguments to turn switch kGCrunOnTrain to kTRUE"<<endl;
    return;
  }
  if(kGCrunOnTrain == kFALSE){
    if(kGCnumberOfFilesToAnalyze==0){
      ifstream dataInStream;
      dataInStream.open(kGCdataList.Data());
      if ( !dataInStream ){
	cout<<"Data list file does not exist: "<<kGCdataList.Data()<<endl;
	return 0;
      }
      string line;
      while ( !dataInStream.eof() )
	{
	  getline(dataInStream, line);
	  if(line.compare("") != 0){//checks if there is an empty line in the data list
	    kGCnumberOfFilesToAnalyze++;
	  }
	}
    }
    cout<<"Number Of files to analyze: "<<kGCnumberOfFilesToAnalyze<<endl;
		
    build();//build (if necessary) and load the libraries needed
    LoadLibraries();
    gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C"); // load the CreateChain macro
  }
		
  if(!kGCrunOnTrain){
    // for the train leave this to the steering macro
    AliLog::SetGlobalDebugLevel(0);
    AliLog::SetGlobalLogLevel(AliLog::kFatal);
  }
  // ------------------------------------------------------------------------
		
    // for CF
		
    //Container def.
    const Double_t ptmin = kGCfirstXBinPt;
    const Double_t ptmax = kGClastXBinPt;
    const Double_t etamin = kGCfirstXBinEta;
    const Double_t etamax = kGClastXBinEta;
    const Double_t massmin = kGCfirstXBinPi0Mass;
    const Double_t massmax = kGClastXBinPi0Mass;
		
		
    // sensitive variables
    UInt_t ipt = 0;
    UInt_t ieta = 1;
    UInt_t imass = 2;	
		
    //how many selection steps 
    UInt_t nstep = 20;
    const Int_t nvar = 3;

    Int_t kGCnXBinsPtCF=40;
    Int_t kGCnXBinsEtaCF=8;
    Int_t kGCnXBinsPi0MassCF=10;

    if(!kGCrunCF){
      nstep=1;
      kGCnXBinsPtCF=1;
      kGCnXBinsEtaCF=1;
      kGCnXBinsPi0MassCF=1;
    }
    const Int_t nbin0 = kGCnXBinsPtCF;  // do not use same variable for CF than for histogram
    const Int_t nbin1 = kGCnXBinsEtaCF; // do not use same variable for CF than for histogram
    const Int_t nbin2 = kGCnXBinsPi0MassCF; // do not use same variable for CF than for histogram	 	
		
    //arrays for the number of bins in each dimension
    Int_t iBin[nvar];
    iBin[0] = nbin0;
    iBin[1] = nbin1;
    iBin[2] = nbin2;	
		
    //arrays for lower bounds
    Double_t *binLim0 = new Double_t[nbin0+1];
    Double_t *binLim1 = new Double_t[nbin1+1];	
    Double_t *binLim2 = new Double_t[nbin2+1];	
		
    // values for lower bounds
    for(Int_t i = 0; i <= nbin0; i++) binLim0[i] = ptmin + (ptmax - ptmin)/nbin0*i;
    for(Int_t i = 0; i <= nbin1; i++) binLim1[i] = etamin + (etamax - etamin)/nbin1*i;
    for(Int_t i = 0; i <= nbin2; i++) binLim2[i] = massmin + (massmax - massmin)/nbin2*i;
		
    // create container
    AliCFContainer *container = new AliCFContainer("container","container for gammaconversion", nstep,nvar,iBin);
    container->SetBinLimits(ipt,binLim0);
    container->SetBinLimits(ieta,binLim1);
    container->SetBinLimits(imass,binLim2);	
		
    AliCFManager *man = new AliCFManager();
    man->SetParticleContainer(container);
    for(int i = 0;i<nstep;i++){
      TObjArray *tmp = new TObjArray(0);
      man->SetParticleCutsList(i,tmp) ;
    }

    // end ---------------------------------------------------------------------------
		
		

	
  AliGammaConversionHistograms* histograms = new AliGammaConversionHistograms();  
  AddHistograms(histograms);
  	
  // Create the Analysis manager
  AliAnalysisManager *mgr =NULL;
  if(kGCrunOnTrain == kFALSE){
    mgr  = new AliAnalysisManager("My Manager", "My Analysis");
  }
  else{
    mgr = AliAnalysisManager::GetAnalysisManager();
  }
	
  if (!mgr) {
    ::Error("ConfigGammaConversion", "No analysis manager to connect to.");
    return NULL;
  }
  if(kGCrunOnTrain == kTRUE){
    if (!mgr->GetInputEventHandler()) {
      ::Error("ConfigGammaConversion", "This task requires an input event handler");
      return NULL;
    }
  }
  AliESDInputHandler* inpHandler = NULL;
	
  if(kGCrunOnTrain == kFALSE){
    // Define Input Event Handler 
    inpHandler = new AliESDInputHandler();
  }
  // Define MC Truth Event Handler
  AliMCEventHandler* mcHandler = NULL;
  if(kGCdoMCTruth){
    if(kGCrunOnTrain == kFALSE){
      mcHandler = new AliMCEventHandler();
    }
    else{
      mcHandler = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
    }
    if (!mcHandler) {
      ::Error("", "No MC handler connected");
      return NULL;
    }
  }
	

	
  if(kGCrunOnTrain == kFALSE){
    mgr->SetInputEventHandler  (inpHandler);
    mgr->SetMCtruthEventHandler(mcHandler);
  }
  // Be sure you are told what you are doing
  // mgr->SetDebugLevel(10);
	
  // Declare Common Input Tchain
  AliAnalysisDataContainer *cinput1 = NULL;
  
  if(kGCusePWG4PartCorr){
  
    if(kGCrunOnTrain == kFALSE){
      cinput1 = mgr->CreateContainer("GammaConvChain",TChain::Class(),AliAnalysisManager::kInputContainer);
    
    } else{
      cinput1 = cin_esd;
    
    }
  }
  else{
      cinput1 = mgr->GetCommonInputContainer();
  }
	
  // Common Output Tree in common ??????default?????? output file
  // CKB kGCusePWG4PartCorr and writestandard are not mutually exclusive?
  AliAnalysisDataContainer *coutput1 = NULL;

  if(kGCWriteAOD) {
    if(kGCrunOnTrain && !(kGCrunOnGsiTrain)) {
      coutput1 = mgr->GetCommonOutputContainer();
    } else {
      coutput1 = mgr->CreateContainer("GammaConvTree",TTree::Class(),AliAnalysisManager::kOutputContainer, "default");  
    }
  }
	
  // Private output objects
  if(kGCoutputFileName.Contains(".root")){
    kGCoutputFileName.ReplaceAll(".root","");
  }
  if(kGCoutputFileAppendix.Contains(".root")){
    kGCoutputFileAppendix.ReplaceAll(".root","");
  }
  //TString fileOut = kGCoutputFileName + kGCoutputFileAppendix + ".root";

  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  cout<<"Analyis cut selection ID is: "<<kGCAnalysisCutSelectionId.Data()<<endl;
  //  outputfile += Form(":PWG4_GammaConversion_%s",kGCAnalysisCutSelectionId.Data());
  outputfile += Form(":PWG4_GammaConversion");

//   if(kGCrunNeutralMeson==kTRUE) outputfile +="1";  else outputfile +="0";

//   if(kGCrunJet==kTRUE) outputfile +="1"; else outputfile +="0";

//   if(kGCrunChic==kTRUE) outputfile +="1"; else outputfile +="0";

//   if(kGCrunCF==kTRUE) outputfile +="1"; else outputfile +="0";

//   if(kGCcalculateBackground==kTRUE) outputfile +="1"; else outputfile +="0";

//   if(kGCdoNeutralMesonV0MCCheck==kTRUE) outputfile +="1"; else outputfile +="0";

//   if(kGCrunOmegaMeson==kTRUE) outputfile +="1"; else outputfile +="0";

//   if(kGCrunRES==kTRUE) outputfile +="1"; else outputfile +="0";

  outputfile += Form("_%s",kGCAnalysisCutSelectionId.Data());

  cout<<"Ouput file::"<<  outputfile <<endl;
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("histogramsAliGammaConversion_%s",kGCAnalysisCutSelectionId.Data()), TList::Class(),AliAnalysisManager::kOutputContainer, outputfile);
  // for CF
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("GammaConvccontainer0_%s",kGCAnalysisCutSelectionId.Data()),AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile);
	
  //------------------------ END: Define input/output handlers ---------------------------------------------------
	
  //check for errors in the specified data
  if(kGCuseKFParticle == kTRUE && kGCuseESDTrack == kTRUE){
    //Print warning, cannot use both
    ::Error("ConfigGammaConversion","Both kGCuseKFParticle and kGCuseESDTracks can be true at the same time")
      }
  if(kGCuseKFParticle == kFALSE && kGCuseESDTrack == kFALSE){
    //Print warning, one have to be specified
    ::Error("ConfigGammaConversion","Both kGCuseKFParticle and kGCuseESDTracks can be false at the same time")
      }
	

  if(!SetAnalysisCutSelection(kGCAnalysisCutSelectionId)){
    cout<<"Error in analysis cut selection"<<endl;
    return 0;
  }
	
  //Create the V0Reader
  AliV0Reader * v0Reader = new AliV0Reader();
  if(kGCuseKFParticle){
    v0Reader->UseKFParticle();
  }
  else if(kGCuseESDTrack){
    v0Reader->UseESDTrack();
  }

  v0Reader->SetMaxVertexZ(kGCmaxVertexZ);
  v0Reader->SetNegativeTrackPID(kGCpidOfNegativeTrack);
  v0Reader->SetPositiveTrackPID(kGCpidOfPositiveTrack);
  v0Reader->SetMaxRCut(kGCmaxRCut);
  v0Reader->SetMinRCut(kGCminRCut);
  v0Reader->SetEtaCut(kGCetaCut);
  v0Reader->SetRapidityMesonCut(kGCrapidityCutMeson);
  v0Reader->SetPtCut(kGCptCut);
  v0Reader->SetSinglePtCut(kGCsingleptCut);
  v0Reader->SetLineCutZRSlope(kGCLineCutZRSlope);
  v0Reader->SetLineCutZValue(kGCLineCutZValue);	
  v0Reader->SetMaxZCut(kGCmaxZCut);	
  v0Reader->SetMinClsTPCCut(kGCminClsTPCCut);	
  v0Reader->SetMinClsTPCCutToF(kGCminClsTPCCutToF);	
  v0Reader->SetChi2CutConversion(kGCchi2CutConversion);
  v0Reader->SetChi2CutMeson(kGCchi2CutMeson);
  v0Reader->SetAlphaCutMeson(kGCalphaCutMeson);
  v0Reader->SetAlphaMinCutMeson(kGCalphaMinCutMeson);
  v0Reader->SetPIDProbability(kGCprobElectron);
  v0Reader->SetXVertexCut(kGCxVertexCut);
  v0Reader->SetYVertexCut(kGCyVertexCut);
  v0Reader->SetZVertexCut(kGCzVertexCut);
  v0Reader->SetSigmaMass(kGCsigmaCutGammaMass);
  v0Reader->SetUseImprovedVertex(kGCuseImprovedVertex);
  v0Reader->SetDoMCTruth(kGCdoMCTruth);
  v0Reader->SetUseOwnXYZCalculation(kGCdoOwnXYZCalculation);
  v0Reader->SetUseChargedTracksMultiplicityForBG(kGCUseTrackMultiplicityForBG);
  // for CF
  v0Reader->SetCFManager(man);
	
  // for dEdx N sigma Cut
  v0Reader->SetDodEdxSigmaCut(kGCdodEdxSigmaCut);
  v0Reader->SetPIDnSigmaAboveElectronLine(kGCPIDnSigmaAboveElectronLine);
  v0Reader->SetPIDnSigmaBelowElectronLine(kGCPIDnSigmaBelowElectronLine);
  v0Reader->SetPIDnSigmaAbovePionLine(kGCPIDnSigmaAbovePionLine);
  v0Reader->SetPIDMinPnSigmaAbovePionLine(kGCPIDMinPnSigmaAbovePionLine);
  v0Reader->SetPIDMaxPnSigmaAbovePionLine(kGCPIDMaxPnSigmaAbovePionLine);
  v0Reader->SetOnFlyFlag(kGCUseOnFlyV0Finder);
  v0Reader->SetCalculateBackground(kGCcalculateBackground);

  v0Reader->SetDoTOFsigmaCut(kGCuseTOFpid); // RRnewTOF
  v0Reader->SetTofPIDnSigmaAboveElectronLine(kGCtofPIDnSigmaAboveElectronLine); // RRnewTOF
  v0Reader->SetTofPIDnSigmaBelowElectronLine(kGCtofPIDnSigmaBelowElectronLine); // RRnewTOF

  // for the rejection at LowP based on Dedx signal

  v0Reader->SetDoKaonRejectionLowP(kGCdoKaonRejectionLowP);
  v0Reader->SetDoProtonRejectionLowP(kGCdoProtonRejectionLowP);
  v0Reader->SetDoPionRejectionLowP(kGCdoPionRejectionLowP);
  v0Reader->SetPIDnSigmaAtLowPAroundKaonLine(kGCPIDnSigmaAtLowPAroundKaonLine);
  v0Reader->SetPIDnSigmaAtLowPAroundProtonLine(kGCPIDnSigmaAtLowPAroundProtonLine);
  v0Reader->SetPIDnSigmaAtLowPAroundPionLine(kGCPIDnSigmaAtLowPAroundPionLine);
  v0Reader->SetPIDMinPKaonRejectionLowP(kGCPIDMinPKaonRejectionLowP);
  v0Reader->SetPIDMinPProtonRejectionLowP(kGCPIDMinPProtonRejectionLowP);
  v0Reader->SetPIDMinPPionRejectionLowP(kGCPIDMinPPionRejectionLowP);
  v0Reader->SetDoQtGammaSelection(kGCdoQtGammaSelection);
  v0Reader->SetDoHighPtQtGammaSelection(kGCdoHighPtQtGammaSelection); // RRnew
  v0Reader->SetQtMax(kGCQtMax);
  v0Reader->SetHighPtQtMax(kGCHighPtQtMax); // RRnew
  v0Reader->SetPtBorderForQt(kGCPtBorderForQt); // RRnew
  kGCNEventsForBGCalculation= kGCnumberOfRotationEventsForBG;
  cout<< "number of Events used for mixing::"<<kGCNEventsForBGCalculation<<endl;
  v0Reader->SetNEventsForBG(kGCNEventsForBGCalculation);

  // Create the GammaConversionTask


  AliAnalysisTaskGammaConversion *gammaconversion = 
    new AliAnalysisTaskGammaConversion(Form("GammaConversionTask_%s",kGCAnalysisCutSelectionId.Data()));

  cout<<"name of Task::"<< Form("GammaConversionTask_%s",kGCAnalysisCutSelectionId.Data())<< " "<<gammaconversion->GetName() <<endl;
  gammaconversion->SetDebugLevel(0);
	
  gammaconversion->SetWriteNtuple(kGCwriteNtuple);
	
  gammaconversion->SetV0Reader(v0Reader);
  gammaconversion->SetCalculateBackground(kGCcalculateBackground);
  gammaconversion->Init();
	
  gammaconversion->SetElectronMass(kGCelectronMass);
  gammaconversion->SetGammaMass(kGCgammaMass);
  gammaconversion->SetPi0Mass(kGCpi0Mass);
  gammaconversion->SetEtaMass(kGCetaMass);
	
  gammaconversion->SetGammaWidth(kGCgammaWidth);
  gammaconversion->SetPi0Width(kGCpi0Width);
  gammaconversion->SetEtaWidth(kGCetaWidth);
	
  gammaconversion->SetMinOpeningAngleGhostCut(kGCminOpeningAngleGhostCut);
  
	
  Double_t lowPtMapping=0.4;
  Double_t highPtMapping=1.5;
  gammaconversion->SetLowPtMapping(lowPtMapping);
  gammaconversion->SetHighPtMapping(highPtMapping);

  // define the width constraint used by KF particle.
  Double_t gammaWidth = 0.01;
  Double_t pi0Width   = 0.01;
  Double_t etaWidth   = 0.01;
	
  gammaconversion->SetHistograms(histograms);
  v0Reader->SetHistograms(histograms);// also give the pointer to the v0reader, for debugging cuts

  gammaconversion->SetTriggerCINT1B(kGCtriggerCINT1B);
  gammaconversion->SetDoMCTruth(kGCdoMCTruth);
	
  gammaconversion->SetDoNeutralMeson(kGCrunNeutralMeson);
  gammaconversion->SetDoNeutralMesonV0MCCheck(kGCdoNeutralMesonV0MCCheck);
  gammaconversion->SetDoJet(kGCrunJet);
  gammaconversion->SetDoChic(kGCrunChic);
  gammaconversion->SetDoOmegaMeson(kGCrunOmegaMeson);
  gammaconversion->SetRecalculateV0ForGamma(kGCRecalculateV0ForGamma);
  gammaconversion->SetUseChargedTracksMultiplicityForBG(kGCUseTrackMultiplicityForBG);
  gammaconversion->SetMoveParticleAccordingToVertex(kGCMoveParticlesAccordingToVertex);
  gammaconversion->SetApplyChi2Cut(kGCApplyChi2Cut);
  cout<<"NumberOfDegrees in rotatation method set to: "<<kGCnDegreeRotationPMForBG<<endl;
  gammaconversion->SetPMDegreesBG(kGCnDegreeRotationPMForBG);
  if(kGCUseRotationMethodInBG){
    cout<<"Using rotation method for bg."<<endl;
  }
  else{
    cout<<"Using mixed event for bg."<<endl;
  }
  if(kGCUseTrackMultiplicityForBG){
    cout<<"Using track multiplicity for bck"<<endl;
  }else{
    cout<<"Using V0 multiplicity for bck"<<endl;
  }
  gammaconversion->SetDoRotation(kGCUseRotationMethodInBG);
  cout<<"Using :"<<kGCnumberOfRotationEventsForBG<<" rotations in bg calculation"<<endl;
  gammaconversion->SetNumberOfRotationsBG(kGCnumberOfRotationEventsForBG);
  gammaconversion->SetCheckBGProbability(kGCdoBGProbability);

  gammaconversion->SetRemovePileUp(kGCRemovePileUp);
  gammaconversion->SetSelectV0AND(kGCSelectV0AND);
  gammaconversion->SetUseMultiplicity(kGCUseMultiplicity);
  gammaconversion->SetUseMultiplicityBin(kGCUseMultiplicityBin);

  v0Reader->SetIsHeavyIon(kGCIsHeavyIon);
  v0Reader->SetUseCorrectedTPCClsInfo(kGCUseCorrectedTPCClsInfo);
  gammaconversion->SetUseCentrality(kGCUseCentrality);
  if(kGCUseCentrality){
    gammaconversion->SetUseCentralityBin(kGCUseCentralityBin);
  }




  // for CF
  gammaconversion->SetCFManager(man);
  gammaconversion->SetDoCF(kGCrunCF);
  v0Reader->SetDoCF(kGCrunCF);

  // Add task to the manager 
  mgr->AddTask(gammaconversion);


  // Define Output Event Handler and add
  if(kGCWriteAOD){
    gammaconversion->SetForceAOD(kGCForceAOD);
    gammaconversion->SetAODBranchName(Form("GammaConv_%s", kGCAnalysisCutSelectionId.Data()));


    if( kGCrunOnTrain ) {
      
      AliAODHandler * aodHandler = dynamic_cast<AliAODHandler*>(mgr->GetOutputEventHandler());
      if(!aodHandler) {
	::Error("This task requires an AOD handler");
	return NULL;
      }
      
      gammaconversion->SetDeltaAODFileName(kGCDeltaAODFilename);
      
      if(kGCDeltaAODFilename.Length() > 0) {
	mgr->RegisterExtraFile(kGCDeltaAODFilename.Data());
      }
      
    } else {
      if(kGCDeltaAODFilename.Length() == 0 ) {
	cout << "Error:: Need a file name for the AOD"<<endl;
	return NULL;
      }
      AliAODHandler* aodHandler = new AliAODHandler();
      aodHandler->SetOutputFileName(kGCDeltaAODFilename);
      aodHandler->SetCreateNonStandardAOD();
      mgr->SetOutputEventHandler(aodHandler);	
    }
  }  else {
    gammaconversion->SetCreateAOD(kFALSE);
  }

  // Connect I/O to the task
  mgr->ConnectInput (gammaconversion, 0, cinput1);
  if(mgr->GetCommonOutputContainer())
    mgr->ConnectOutput(gammaconversion, 0, mgr->GetCommonOutputContainer());
  
  mgr->ConnectOutput(gammaconversion, 1, coutput2);
  mgr->ConnectOutput(gammaconversion, 2, coutput3);
  
  if(kGCRunGammaJetTask) {
    AliAnalysisTaskGammaJet * gammaJetTask = new AliAnalysisTaskGammaJet("GammaJetTask");
    if(kGCrunOnTrain) {
      gammaJetTask->SetDeltaAODFileName(kGCDeltaAODFileName);
    }
    
    mgr->ConnectInput(gammaJetTask, 0, cinput1);
    AliAnalysisDataContainer *coutputgj = mgr->CreateContainer("chistpt", TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResultsConvJets.root");
    mgr->ConnectOutput(gammaJetTask, 1, coutputgj);
  }


  if( kGCrunDalitz ){
   
   gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskGammaConvDalitz.C");
   AddTaskGammaConvDalitz( v0Reader, kGCcalculateBackground, kGCRunStandalone );

  }

  if(kGCrunOnTrain == kFALSE){
    if(kGCdataList.IsNull()){
      cout<<"Data list is not set, aborting."<<endl;
      return;
    }
    /*
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
      if(kGCdoMCTruth)physSelTask->GetPhysicsSelection()->SetAnalyzeMC();
      physSelTask->GetPhysicsSelection()->AddBackgroundIdentification(new AliBackgroundSelection());
      gammaconversion->SelectCollisionCandidates();	
    */

    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kGCdoMCTruth,kTRUE);
    gammaconversion->SelectCollisionCandidates(); 

    //    if(kGCrunOnTrain == kFALSE){
      TChain* chain= CreateESDChain(kGCdataList,kGCnumberOfFilesToAnalyze);
			
      mgr->InitAnalysis();
			
      mgr->PrintStatus();
			
      mgr->StartAnalysis("local",chain);
      //    }
  }
  return gammaconversion;
}

void LoadLibraries() {
	 	 
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom");
  
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG4GammaConv.so");
  
  //  gSystem->ChangeDirectory(pwd.Data());
  
}
void build() {
	
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom");
	
  TString pwd = gSystem->WorkingDirectory();
	
  ////
  //Setting up STEERBase.par//
  ////
  cout<<"compiling STEERBase"<<endl;
  setupPar("STEERBase");
  gSystem->Load("libSTEERBase.so");


   ////
  //Setting up ESD.par//
  ////
  cout<<"compiling ESD"<<endl;
  setupPar("ESD");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");

	
  ////
  //Setting up AOD.par//
  ////
  cout<<"compiling AOD"<<endl;
  setupPar("AOD");
  gSystem->Load("libAOD.so");
	
  ////
  //Setting up ANALYSIS.par//
  ////
  cout<<"compiling ANALYSIS"<<endl;
  setupPar("ANALYSIS");
  gSystem->Load("libANALYSIS.so");
	
  ////
  //Setting up ANALYSISalice.par//
  ////
  cout<<"compiling ANALYSISalice"<<endl;
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice.so");
	
  ////
  //Setting up CORRFW.par//
  ////
  cout<<"compiling CORRFW"<<endl;
  setupPar("CORRFW");
  gSystem->Load("libCORRFW.so");
	
  ////
  //Setting up PWG4GammaConv.par//
  ////
  cout<<"compiling PWG4GammaConv"<<endl;
  setupPar("PWG4GammaConv");
  gSystem->Load("libPWG4GammaConv.so");
  
  gSystem->ChangeDirectory(pwd.Data());
}

Int_t setupPar(const char* pararchivename) {
  ///////////////////
  // Setup PAR File//
  ///////////////////
  if (pararchivename) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchivename);
		
    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");
			
      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
	Error("runAnalysis","Cannot Build the PAR Archive! - Abort!");
	return -1;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*******************************\n");
      printf("*** Setup PAR archive       ***\n");
      printf("*******************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
		
    gSystem->ChangeDirectory("../");
  }                                                                                                                                               
  return 1;
}



void AddHistograms(AliGammaConversionHistograms *histograms){
  //---------------------------------------------- Jets ---------------------------------------------------------
  if(kGCrunJet == kTRUE){
    if (kGCplotdPhiHdrGam == kTRUE){
      histograms->AddHistogram("ESD_dphiHdrGam","ESD_dphiHdrGam", kGCnXBinsdphiHdrGam,kGCfirstXBindphiHdrGam,kGClastXBindphiHdrGam,"dphiHdrGam (rad)","Counts");
    }
		
    if (kGCplotdPhiHdrGamIsolated == kTRUE){
      histograms->AddHistogram("ESD_dphiHdrGamIsolated","ESD_dphiHdrGamIsolated",  kGCnXBinsdphiHdrGam,kGCfirstXBindphiHdrGam,kGClastXBindphiHdrGam,"dphiHdrGamIsolated (rad)","Counts");
    }
		
    if (kGCplotMinimumIsoDistance == kTRUE){
      histograms->AddHistogram("ESD_MinimumIsoDistance","ESD_MinimumIsoDistance", kGCnXBinsMinimumIsoDistance,kGCfirstXBinMinimumIsoDistance,kGClastXBinMinimumIsoDistance,"Minimum Iso Distance (rad)","Counts");
    }
		
    if (kGCplotFFzHdrGam == kTRUE){
      histograms->AddHistogram("ESD_FFzHdrGam","ESD_FFzHdrGam", kGCnXBinsFFzHdrGam, kGCfirstXBinFFzHdrGam,kGClastXBinFFzHdrGam,"FFz Hdr Gam","Counts");
    }
		
    if (kGCplotImbalanceHdrGam == kTRUE){
      histograms->AddHistogram("ESD_ImbalanceHdrGam","ESD_ImbalanceHdrGam", kGCnXBinsImbalanceHdrGam, kGCfirstXBinImbalanceHdrGam,kGClastXBinImbalanceHdrGam,"Imbalance Hdr Gam","Counts");
    }
  }//end if(kGCrunJet)
	
  //---------------------------------------------- Chi_c ---------------------------------------------------------
  if(kGCrunChic){
		
    if(kGCplotESDInvMassePluseMinus == kTRUE){histograms->AddHistogram("ESD_InvMass_ePluseMinus","",kGCnXBinsJPsiMass, kGCfirstXBinJPsiMass, kGClastXBinJPsiMass, "",
								       "");}
    if(kGCplotESDInvMassePluseMinus == kTRUE){histograms->AddHistogram("ESD_InvMass_ePluseMinusTest","",kGCnXBinsJPsiMass, kGCfirstXBinJPsiMass, kGClastXBinJPsiMass,
								       "","");}
    if(kGCplotESDInvMassePluseMinus == kTRUE){histograms->AddHistogram("ESD_InvMass_xPlusxMinus","",kGCnXBinsJPsiMass, kGCfirstXBinJPsiMass, kGClastXBinJPsiMass, "",
								       "");}
    if(kGCplotESDElectronPosNegPt == kTRUE){histograms->AddHistogram("ESD_ElectronPosNegPt","",kGCnXBinsEPosNegPt,kGCfirstXBinEPosNegPt,kGClastXBinEPosNegPt,"","");}
    if(kGCplotESDElectronPosNegEta == kTRUE){histograms->AddHistogram("ESD_ElectronPosNegEta","",kGCnXBinsEPosNegEta,kGCfirstXBinEPosNegEta,kGClastXBinEPosNegEta,"","
																		  ");}
		
    if(kGCplotESDElectronPosNegPt == kTRUE){histograms->AddHistogram("ESD_ElectronPosPt","",kGCnXBinsEPosNegPt,kGCfirstXBinEPosNegPt,kGClastXBinEPosNegPt,"","");}
    if(kGCplotESDElectronPosNegPt == kTRUE){histograms->AddHistogram("ESD_ElectronNegPt","",kGCnXBinsEPosNegPt,kGCfirstXBinEPosNegPt,kGClastXBinEPosNegPt,"","");}
		
    if(kGCplotESDElectronPosNegAngle == kTRUE){histograms->AddHistogram("ESD_ElectronPosNegJPsiAngle","",kGCnXBinsEPosNegAngle,kGCfirstXBinEPosNegAngle,kGClastXBinEPosNegAngle,"","");}

    if(kGCplotESDePoseNegAngle == kTRUE){histograms->AddHistogram("ESD_eNegePosAngleBeforeCut","",kGCnXBinsEPosNegAngle,kGCfirstXBinEPosNegAngle,kGClastXBinEPosNegAngle,"","");}
    if(kGCplotESDePoseNegAngle == kTRUE){histograms->AddHistogram("ESD_eNegePosAngleAfterCut","",kGCnXBinsEPosNegAngle,kGCfirstXBinEPosNegAngle,kGClastXBinEPosNegAngle,"","");}
    if(kGCplotESDInvMassGammaePluseMinusChiC == kTRUE) {histograms->AddHistogram("ESD_InvMass_GammaePluseMinusChiC","",kGCnXBinsChicMass,kGCfirstXBinChicMass,kGClastXBinChicMass,"","");}
    if(kGCplotESDInvMassGammaePluseMinusChiC == kTRUE) {histograms->AddHistogram("ESD_InvMass_GammaePluseMinusChiCDiff","",kGCnXBinsChicMass,kGCfirstXBinChicMass,kGClastXBinChicMass,"","");}
    if(kGCplotESDInvMassGammaePluseMinusPi0 == kTRUE) {histograms->AddHistogram("ESD_InvMass_GammaePluseMinusPi0","",kGCnXBinsPi0Mass,kGCfirstXBinPi0Mass,kGClastXBinPi0Mass,"","");}
    if(kGCplotESDElectronPosNegPi0Angle == kTRUE){histograms->AddHistogram("ESD_ElectronPosNegPi0Angle","",kGCnXBinsEPosNegAngle,kGCfirstXBinEPosNegAngle,kGClastXBinEPosNegAngle,"","");}
		
    if(kGCplotESDEPosBackground == kTRUE){histograms->AddHistogram("ESD_EPosBackground","",kGCnXBinsEBackground,kGCfirstXBinEBackground,kGClastXBinEBackground,"","");}
		
    if(kGCplotESDEPosBackground == kTRUE){histograms->AddHistogram("ESD_EPosENegNoJPsiBG","",kGCnXBinsEBackground,kGCfirstXBinEBackground,kGClastXBinEBackground,"","");}
		
		
    if(kGCplotESDENegBackground == kTRUE){histograms->AddHistogram("ESD_ENegBackground","",kGCnXBinsEBackground,kGCfirstXBinEBackground,kGClastXBinEBackground,"","");}
    if(kGCplotESDEPosENegBackground == kTRUE){histograms->AddHistogram("ESD_EPosENegBackground","",kGCnXBinsEBackground,kGCfirstXBinEBackground,kGClastXBinEBackground,"","");}
    if(kGCplotESDEPosENegBackgroundCut == kTRUE){histograms->AddHistogram("ESD_EPosENegBackgroundCut","",kGCnXBinsEBackgroundCut,kGCfirstXBinEBackgroundCut,kGClastXBinEBackgroundCut,"","");}
		
    if(kGCplotESDEPosENegGammaBackgroundMX == kTRUE){histograms->AddHistogram("ESD_EPosENegGammaBackgroundMX","",kGCnXBinsEBackground,kGCfirstXBinEBackground,kGClastXBinEBackground,"","");}
    if(kGCplotESDEPosENegGammaBackgroundMX == kTRUE){histograms->AddHistogram("ESD_EPosENegGammaBackgroundMXDiff","",kGCnXBinsEBackground,kGCfirstXBinEBackground,kGClastXBinEBackground,"","");}
		
    if(kGCplotTableElectrons == kTRUE){ histograms->AddTable("Table_Electrons","",kGCnElementsElectronTable,kGCelectronTable);}

    if(kGCdoMCTruth){
      if(kGCplotMCElectronPosNegPt == kTRUE){histograms->AddHistogram("MC_ElectronPosNegPt","",kGCnXBinsEPosNegPt,kGCfirstXBinEPosNegPt,kGClastXBinEPosNegPt,"","");}
      if(kGCplotMCElectronPosNegEta == kTRUE){histograms->AddHistogram("MC_ElectronPosNegEta","",kGCnXBinsEPosNegEta,kGCfirstXBinEPosNegEta,kGClastXBinEPosNegEta,"","");}
      if(kGCplotMCElectronPosNegJPsiAngle == kTRUE){histograms->AddHistogram("MC_ElectronPosNegJPsiAngle","",kGCnXBinsEPosNegAngle,kGCfirstXBinEPosNegAngle,kGClastXBinEPosNegAngle,"","");}
      if(kGCplotMCElectronPosNegPi0Angle == kTRUE){histograms->AddHistogram("MC_ElectronPosNegPi0Angle","",kGCnXBinsEPosNegAngle,kGCfirstXBinEPosNegAngle,kGClastXBinEPosNegAngle,"","");}
    }

  }// end kGCrunChic
	

  //---------------------------------------------- Gamma from Decay ----------------------------------------------------
  if(kGCplotMCGammaFromDecay == kTRUE && kGCdoMCTruth == kTRUE){
     histograms->AddHistogram("MC_DecayPi0Gamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");  
     histograms->AddHistogram("MC_DecayRho0Gamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
     histograms->AddHistogram("MC_DecayEtaGamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
     histograms->AddHistogram("MC_DecayOmegaGamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
     histograms->AddHistogram("MC_DecayK0sGamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
     histograms->AddHistogram("MC_DecayEtapGamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
     histograms->AddHistogram("MC_DecayPhiGamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
     histograms->AddHistogram("MC_DecayAllGamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
  }

  //---------------------------------------------- Neutral Meson ---------------------------------------------------------
  if(kGCrunNeutralMeson){
		
    // Histograms from esd tracks	
    if(kGCplotESDEEnergy == kTRUE){ histograms->AddHistogram("ESD_E_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
    if(kGCplotESDEPt == kTRUE){ histograms->AddHistogram("ESD_E_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
    if(kGCplotESDEEta == kTRUE){ histograms->AddHistogram("ESD_E_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
    if(kGCplotESDEPhi == kTRUE){ histograms->AddHistogram("ESD_E_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
    if(kGCplotESDENTPCClusters == kTRUE){ histograms->AddHistogram("ESD_E_nTPCClusters" ,"" , kGCnXBinsNTPCClusters, kGCfirstXBinNTPCClusters, kGClastXBinNTPCClusters, "", "");}
    if(kGCplotESDENITSClusters == kTRUE){ histograms->AddHistogram("ESD_E_nITSClusters" ,"" , kGCnXBinsNITSClusters, kGCfirstXBinNITSClusters, kGClastXBinNITSClusters, "", "");}
    if(kGCplotESDENTPCClustersToFP== kTRUE){ histograms->AddHistogram("ESD_E_nTPCClustersToFP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsClsToF, kGCfirstYBinClsToF, kGClastYBinClsToF,"", "",0);}
    if(kGCplotESDENTPCClustersToFR== kTRUE){ histograms->AddHistogram("ESD_E_nTPCClustersToFR" ,"" ,kGCnXBinsR, kGCfirstXBinR, kGClastXBinR,kGCnYBinsClsToF, kGCfirstYBinClsToF, kGClastYBinClsToF,"", "",0);}
    if(kGCplotESDENTPCClustersToFR== kTRUE){ histograms->AddHistogram("ESD_TrueConversion_E_nTPCClustersToFR" ,"" ,kGCnXBinsR, kGCfirstXBinR, kGClastXBinR,kGCnYBinsClsToF, kGCfirstYBinClsToF, kGClastYBinClsToF,"", "",0);}
    if(kGCplotESDETPCchi2 == kTRUE){ histograms->AddHistogram("ESD_E_TPCchi2" ,"" , kGCnXBinsTPCchi2, kGCfirstXBinTPCchi2, kGClastXBinTPCchi2, "", "");}

		
    if(kGCplotESDPEnergy == kTRUE){ histograms->AddHistogram("ESD_P_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
    if(kGCplotESDPPt == kTRUE){ histograms->AddHistogram("ESD_P_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
    if(kGCplotESDPEta == kTRUE){ histograms->AddHistogram("ESD_P_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
    if(kGCplotESDPPhi == kTRUE){ histograms->AddHistogram("ESD_P_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
    if(kGCplotESDPNTPCClusters == kTRUE){ histograms->AddHistogram("ESD_P_nTPCClusters" ,"" , kGCnXBinsNTPCClusters, kGCfirstXBinNTPCClusters, kGClastXBinNTPCClusters, "", "");}
    if(kGCplotESDPNITSClusters == kTRUE){ histograms->AddHistogram("ESD_P_nITSClusters" ,"" , kGCnXBinsNITSClusters, kGCfirstXBinNITSClusters, kGClastXBinNITSClusters, "", "");}
    if(kGCplotESDPNTPCClustersToFP== kTRUE){ histograms->AddHistogram("ESD_P_nTPCClustersToFP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsClsToF, kGCfirstYBinClsToF, kGClastYBinClsToF,"", "",0);}	
    if(kGCplotESDPNTPCClustersToFR== kTRUE){ histograms->AddHistogram("ESD_P_nTPCClustersToFR" ,"" ,kGCnXBinsR, kGCfirstXBinR, kGClastXBinR,kGCnYBinsClsToF, kGCfirstYBinClsToF, kGClastYBinClsToF,"", "",0);}
    if(kGCplotESDPNTPCClustersToFR== kTRUE){ histograms->AddHistogram("ESD_TrueConversion_P_nTPCClustersToFR" ,"" ,kGCnXBinsR, kGCfirstXBinR, kGClastXBinR,kGCnYBinsClsToF, kGCfirstYBinClsToF, kGClastYBinClsToF,"", "",0);}
    if(kGCplotESDPTPCchi2 == kTRUE){ histograms->AddHistogram("ESD_P_TPCchi2" ,"" , kGCnXBinsTPCchi2, kGCfirstXBinTPCchi2, kGClastXBinTPCchi2, "", "");}

    if(kGCplotESDConvGammaEnergy == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
    if(kGCplotESDConvGammaPt == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
    if(kGCplotESDConvGammaEta == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
    if(kGCplotESDConvGammaPhi == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
    if(kGCplotESDConvGammaMass == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Mass" ,"" ,  kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass, "", "");}
    if(kGCplotESDConvGammaWidth == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Width" ,"" , kGCnXBinsGammaWidth, kGCfirstXBinGammaWidth, kGClastXBinGammaWidth, "", "");}
    if(kGCplotESDConvGammaChi2 == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Chi2" ,"" , kGCnXBinsGammaChi2, kGCfirstXBinGammaChi2, kGClastXBinGammaChi2, "", "");}
    if(kGCplotESDConvGammaNDF == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_NDF" ,"" , kGCnXBinsGammaNDF, kGCfirstXBinGammaNDF, kGClastXBinGammaNDF, "", "");}
    if(kGCplotESDConvGammaRapid == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Rapid" ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
    if(kGCplotESDConvGammaPtvsEta == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Pt_Eta","", kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt,kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta,"","" );}
    if(kGCplotESDConvGammaPtvsChi2 == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Pt_Chi2" ,"" ,kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsGammaChi2, kGCfirstXBinGammaChi2, kGClastXBinGammaChi2, "", "");}
    if(kGCplotESDConvGammaEtavsChi2 == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Eta_Chi2" ,"" ,kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, kGCnXBinsGammaChi2, kGCfirstXBinGammaChi2, kGClastXBinGammaChi2, "", "");}
	histograms->AddHistogram("ESD_ConvGamma_EandP_P_dT" ,"" , kGCnXBinsP, kGCfirstXBinP, kGClastXBinP, kGCnXBinsTOFsignal, kGCfirstXBinTOFsignal, kGClastXBinTOFsignal, "", "",0); // RRnewTOF		
		
		
    if(kGCplotESDConversionR == kTRUE){ histograms->AddHistogram("ESD_Conversion_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
    if(kGCplotESDConversionZR == kTRUE){ histograms->AddHistogram("ESD_Conversion_ZR" ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");}
    if(kGCplotESDConversionXY == kTRUE){ histograms->AddHistogram("ESD_Conversion_XY" ,"" , kGCnXBinsXY, kGCfirstXBinXY, kGClastXBinXY, kGCnYBinsXY, kGCfirstYBinXY, kGClastYBinXY, "", "");}
    if(kGCplotESDConversionXYBeamPipe == kTRUE && kGCLowMemoryConsumption == kFALSE){ histograms->AddHistogram("ESD_Conversion_XY_BeamPipe" ,"" , kGCnXBinsXYBP, kGCfirstXBinXYBP, kGClastXBinXYBP, kGCnYBinsXYBP, kGCfirstYBinXYBP, kGClastYBinXYBP, "", "");}
    if(kGCplotESDConversionRPhiBeamPipe == kTRUE && kGCLowMemoryConsumption == kFALSE){ histograms->AddHistogram("ESD_Conversion_RPhi_BeamPipe" ,"" , kGCnXBinsRPhiBP, kGCfirstXBinRPhiBP, kGClastXBinRPhiBP, kGCnYBinsRPhiBP, kGCfirstYBinRPhiBP, kGClastYBinRPhiBP, "", "");}
    if(kGCplotESDConversionOpeningAngle == kTRUE){ histograms->AddHistogram("ESD_Conversion_OpeningAngle" ,"" , kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");}

    if(kGCplotESDConvGammaCosPointingAngle == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_CosPointingAngle" ,"" , kGCnXBinsCosPointingAngle, kGCfirstXBinCosPointingAngle, kGClastXBinCosPointingAngle, "", "");}
    if(kGCplotESDConvGammaDcaDaugthers == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_DcaDaughters" ,"" , kGCnXBinsDcaDaughters, kGCfirstXBinDcaDaughters, kGClastXBinDcaDaughters, "", "");}
    if(kGCplotESDConvGammaNormDcaDistDaugthers == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_NormDcaDistDaughters" ,"" , kGCnXBinsNormDcaDistDaughters, kGCfirstXBinNormDcaDistDaughters, kGClastXBinNormDcaDistDaughters, "", "");}
    if(kGCplotESDConvGammaLikelihoodAP == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_LikelihoodAP" ,"" , kGCnXBinsLikelihoodAP, kGCfirstXBinLikelihoodAP, kGClastXBinLikelihoodAP, "", "");}
    if(kGCplotESDConvGammaEAsymmetryP== kTRUE){ histograms->AddHistogram("ESD_ConvGamma_E_AsymmetryP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsAsymmetry, kGCfirstYBinAsymmetry, kGClastYBinAsymmetry,"", "");}
    if(kGCplotESDConvGammaPAsymmetryP== kTRUE){ histograms->AddHistogram("ESD_ConvGamma_P_AsymmetryP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsAsymmetry, kGCfirstYBinAsymmetry, kGClastYBinAsymmetry,"", "");}
    if(kGCplotESDConvGammaEdEdxP== kTRUE){ histograms->AddHistogram("ESD_ConvGamma_E_dEdxP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"", "",0);}
    if(kGCplotESDConvGammaPdEdxP== kTRUE){ histograms->AddHistogram("ESD_ConvGamma_P_dEdxP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"", "",0);}

    if(kGCplotESDConvGammaEeProbP== kTRUE){ histograms->AddHistogram("ESD_ConvGamma_E_EProbP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsEProb, kGCfirstYBinEProb, kGClastYBinEProb,"", "",0);}
    if(kGCplotESDConvGammaPeProbP== kTRUE){ histograms->AddHistogram("ESD_ConvGamma_P_EProbP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsEProb, kGCfirstYBinEProb, kGClastYBinEProb,"", "",0);}
    if(kGCplotESDConvGammaEmupiProbP== kTRUE){ histograms->AddHistogram("ESD_ConvGamma_E_mupiProbP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsEProb, kGCfirstYBinEProb, kGClastYBinEProb,"", "",0);}
    if(kGCplotESDConvGammaPmupiProbP== kTRUE){ histograms->AddHistogram("ESD_ConvGamma_P_mupiProbP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsEProb, kGCfirstYBinEProb, kGClastYBinEProb,"", "",0);}

    if(kGCplotESDConvGammaQtAlfa== kTRUE){ histograms->AddHistogram("ESD_ConvGamma_alfa_qt" ,"" ,kGCnXBinsP, kGCfirstXBinAlphaG, kGClastXBinAlpha,kGCnYBinsQt, kGCfirstYBinQt, kGClastYBinQt,"", "");}


		
    if(kGCplotESDTrueDalitzContaminationR == kTRUE){ histograms->AddHistogram("ESD_TrueDalitzContamination_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
    if(kGCplotESDTruePi0DalitzContaminationR == kTRUE){ histograms->AddHistogram("ESD_TrueConvDalitzPi0_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
    if(kGCplotESDTrueEtaDalitzContaminationR == kTRUE){ histograms->AddHistogram("ESD_TrueConvDalitzEta_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
    if(kGCplotESDTrueCombinatorialContaminationR == kTRUE){ histograms->AddHistogram("ESD_TrueConvCombinatorial_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
    if(kGCplotESDTrueCombinatorialElecContaminationR == kTRUE){ histograms->AddHistogram("ESD_TrueConvCombinatorialElec_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
    if(kGCplotESDTrueHadronicContaminationR == kTRUE){ histograms->AddHistogram("ESD_TrueConvHadronicBck_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
    if(kGCplotESDTrueCombinatorialContaminationPt == kTRUE){ histograms->AddHistogram("ESD_TrueConvCombinatorial_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
    if(kGCplotESDTrueCombinatorialElecContaminationPt == kTRUE){ histograms->AddHistogram("ESD_TrueConvCombinatorialElec_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
    if(kGCplotESDTrueHadronicContaminationPt == kTRUE){ histograms->AddHistogram("ESD_TrueConvHadronicBck_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}


    if(kGCplotESDTrueBackground){
       histograms->AddHistogram("ESD_TrueConvCombinatorialDaughter_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt,kGCnXBinsPt,kGCfirstXBinPt , kGClastXBinPt, "", "");
       histograms->AddHistogram("ESD_TrueConvHadronicBckDaughter_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt,kGCnXBinsPt,kGCfirstXBinPt , kGClastXBinPt, "", "");
       histograms->AddHistogram("ESD_TrueConvCombinatorialPiDaughter_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt,kGCnXBinsPt,kGCfirstXBinPt , kGClastXBinPt, "", "");
       histograms->AddHistogram("ESD_TrueConvCombinatorialPiPDaughter_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt,kGCnXBinsPt,kGCfirstXBinPt , kGClastXBinPt, "", "");

       histograms->AddHistogram("ESD_TrueConvCombinatorialPi_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
       histograms->AddHistogram("ESD_TrueConvCombinatorialPi_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");
       histograms->AddHistogram("ESD_TrueConvCombinatorialPiP_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
       histograms->AddHistogram("ESD_TrueConvCombinatorialPiP_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");
       histograms->AddHistogram("ESD_TrueConvCombinatorialElecPi_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
       histograms->AddHistogram("ESD_TrueConvCombinatorialElecPi_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");
       histograms->AddHistogram("ESD_TrueConvMeson_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
       histograms->AddHistogram("ESD_TrueConvMeson_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");
       histograms->AddHistogram("ESD_TrueConvLambda_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
       histograms->AddHistogram("ESD_TrueConvLambda_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");
    }

    histograms->AddHistogram("ESD_TrueConvDalitzPi0_SinglePos_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");
    histograms->AddHistogram("ESD_TrueConvDalitzPi0_SingleNeg_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");
    histograms->AddHistogram("ESD_TrueConvDalitzPi0_SinglePos_kFirst_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");
    histograms->AddHistogram("ESD_TrueConvDalitzPi0_SingleNeg_kFirst_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");
    histograms->AddHistogram("ESD_TrueConversion_SinglePos_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR,  kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt,"", "");
    histograms->AddHistogram("ESD_TrueConversion_SingleNeg_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR,  kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
    histograms->AddHistogram("ESD_TrueConversion_SinglePos_kFirst_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
    histograms->AddHistogram("ESD_TrueConversion_SingleNeg_kFirst_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");


    if(kGCplotESDTrueConvGammaEnergy == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
    if(kGCplotESDTrueConvGammaPt == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
    if(kGCplotESDTrueConvGammaEta == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
    if(kGCplotESDTrueConvGammaPhi == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
    if(kGCplotESDTrueConvGammaMass == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Mass" ,"" ,  kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass, "", "");}
    if(kGCplotESDTrueConvGammaWidth == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Width" ,"" , kGCnXBinsGammaWidth, kGCfirstXBinGammaWidth, kGClastXBinGammaWidth, "", "");}
    if(kGCplotESDTrueConvGammaChi2 == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Chi2" ,"" , kGCnXBinsGammaChi2, kGCfirstXBinGammaChi2, kGClastXBinGammaChi2, "", "");}
    if(kGCplotESDTrueConvGammaNDF == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_NDF" ,"" , kGCnXBinsGammaNDF, kGCfirstXBinGammaNDF, kGClastXBinGammaNDF, "", "");}
    if(kGCplotESDTrueConvGammaRapid == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Rapid" ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
    if(kGCplotESDTrueConvGammaPtvsEta == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Pt_Eta" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt,kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}

    if(kGCplotESDTrueConvGammaPtvsChi2 == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Pt_Chi2" ,"" ,kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsGammaChi2, kGCfirstXBinGammaChi2, kGClastXBinGammaChi2, "", "");}
    if(kGCplotESDTrueConvGammaEtavsChi2 == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Eta_Chi2" ,"" ,kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, kGCnXBinsGammaChi2, kGCfirstXBinGammaChi2, kGClastXBinGammaChi2, "", "");}
		
    if(kGCplotESDTrueConversionR == kTRUE){ histograms->AddHistogram("ESD_TrueConversion_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
    if(kGCplotESDTrueConversionZR == kTRUE){ histograms->AddHistogram("ESD_TrueConversion_ZR" ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");}
    if(kGCplotESDTrueConversionXY == kTRUE){ histograms->AddHistogram("ESD_TrueConversion_XY" ,"" , kGCnXBinsXY, kGCfirstXBinXY, kGClastXBinXY, kGCnYBinsXY, kGCfirstYBinXY, kGClastYBinXY, "", "");}
    if(kGCplotESDTrueConversionOpeningAngle == kTRUE){ histograms->AddHistogram("ESD_TrueConversion_OpeningAngle" ,"" , kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");}

    if(kGCplotESDTrueConvGammaCosPointingAngle == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_CosPointingAngle" ,"" , kGCnXBinsCosPointingAngle, kGCfirstXBinCosPointingAngle, kGClastXBinCosPointingAngle, "", "");}
    if(kGCplotESDTrueConvGammaDcaDaugthers == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_DcaDaughters" ,"" , kGCnXBinsDcaDaughters, kGCfirstXBinDcaDaughters, kGClastXBinDcaDaughters, "", "");}
    if(kGCplotESDTrueConvGammaNormDcaDistDaugthers == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_NormDcaDistDaughters" ,"" , kGCnXBinsNormDcaDistDaughters, kGCfirstXBinNormDcaDistDaughters, kGClastXBinNormDcaDistDaughters, "", "");}
    if(kGCplotESDTrueConvGammaLikelihoodAP == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_LikelihoodAP" ,"" , kGCnXBinsLikelihoodAP, kGCfirstXBinLikelihoodAP, kGClastXBinLikelihoodAP, "", "");}
    if(kGCplotESDTrueConvGammaEAsymmetryP== kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_E_AsymmetryP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsAsymmetry, kGCfirstYBinAsymmetry, kGClastYBinAsymmetry,"", "");}
    if(kGCplotESDTrueConvGammaPAsymmetryP== kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_P_AsymmetryP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsAsymmetry, kGCfirstYBinAsymmetry, kGClastYBinAsymmetry,"", "");}
    if(kGCplotESDTrueConvGammaEdEdxP== kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_E_dEdxP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"", "",0);}
    if(kGCplotESDTrueConvGammaPdEdxP== kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_P_dEdxP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"", "",0);}
    if(kGCplotESDTrueConvGammaQtAlfa== kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_alfa_qt" ,"" ,kGCnXBinsP, kGCfirstXBinAlphaG, kGClastXBinAlpha,kGCnYBinsQt, kGCfirstYBinQt, kGClastYBinQt,"", "");}
		
    if(kGCplotESDTrueConvGammaMCPtEta == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_MC_Pt_Eta" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
    if(kGCplotESDTrueConversionMCZR == kTRUE){ histograms->AddHistogram("ESD_TrueConversion_MC_ZR" ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");}
    if(kGCplotESDTrueConversionMCXY == kTRUE){ histograms->AddHistogram("ESD_TrueConversion_MC_XY" ,"" , kGCnXBinsXY, kGCfirstXBinXY, kGClastXBinXY, kGCnYBinsXY, kGCfirstYBinXY, kGClastYBinXY, "", "");}
	histograms->AddHistogram("ESD_TrueConvGamma_EandP_P_dT" ,"" , kGCnXBinsP, kGCfirstXBinP, kGClastXBinP, kGCnXBinsTOFsignal, kGCfirstXBinTOFsignal, kGClastXBinTOFsignal, "", "",0); // RRnewTOF
	histograms->AddHistogram("ESD_TrueConvCombinatorial_DaughtersNotElec_P_dT" ,"" , kGCnXBinsP, kGCfirstXBinP, kGClastXBinP, kGCnXBinsTOFsignal, kGCfirstXBinTOFsignal, kGClastXBinTOFsignal, "", "",0); // RRnewTOF
	histograms->AddHistogram("ESD_TrueConvHadronicBck_Daughters_P_dT" ,"" , kGCnXBinsP, kGCfirstXBinP, kGClastXBinP, kGCnXBinsTOFsignal, kGCfirstXBinTOFsignal, kGClastXBinTOFsignal, "", "",0); // RRnewTOF		

	histograms->AddHistogram("ESD_NoCutConvGamma_EandP_P_dT" ,"" , kGCnXBinsP, kGCfirstXBinP, kGClastXBinP, kGCnXBinsTOFsignal, kGCfirstXBinTOFsignal, kGClastXBinTOFsignal, "", "",0); // RRnewTOF		
    if(kGCplotESDNoCutAllV0Pt == kTRUE){ histograms->AddHistogram("ESD_NoCutAllV0_Pt" ,"" , kGCnXBinsPt,kGCfirstXBinPt , kGClastXBinPt, "", "");}		
    if(kGCplotESDNoCutConvGammaEnergy == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
    if(kGCplotESDNoCutConvGammaPt == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
    if(kGCplotESDNoCutConvGammaEta == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
    if(kGCplotESDNoCutConvGammaPhi == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
    if(kGCplotESDNoCutConvGammaMass == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Mass" ,"" ,  kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass, "", "");}
    if(kGCplotESDNoCutConvGammaWidth == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Width" ,"" , kGCnXBinsGammaWidth, kGCfirstXBinGammaWidth, kGClastXBinGammaWidth, "", "");}
    if(kGCplotESDNoCutConvGammaChi2 == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Chi2" ,"" , kGCnXBinsGammaChi2, kGCfirstXBinGammaChi2, kGClastXBinGammaChi2, "", "");}
    if(kGCplotESDNoCutConvGammaNDF == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_NDF" ,"" , kGCnXBinsGammaNDF, kGCfirstXBinGammaNDF, kGClastXBinGammaNDF, "", "");}
    if(kGCplotESDNoCutConvGammaRapid == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Rapid" ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
    if(kGCplotESDNoCutConvGammaPtvsEta == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Pt_Eta" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt,kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
    if(kGCplotESDNoCutConvGammaPtvsChi2 == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Pt_Chi2" ,"" ,kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsGammaChi2, kGCfirstXBinGammaChi2, kGClastXBinGammaChi2, "", "");}
    if(kGCplotESDNoCutConvGammaEtavsChi2 == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Eta_Chi2" ,"" ,kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, kGCnXBinsGammaChi2, kGCfirstXBinGammaChi2, kGClastXBinGammaChi2, "", "");}
		
    if(kGCplotESDNoCutConversionR == kTRUE){ histograms->AddHistogram("ESD_NoCutConversion_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
    if(kGCplotESDNoCutConversionZR == kTRUE){ histograms->AddHistogram("ESD_NoCutConversion_ZR" ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");}
    if(kGCplotESDNoCutConversionXY == kTRUE){ histograms->AddHistogram("ESD_NoCutConversion_XY" ,"" , kGCnXBinsXY, kGCfirstXBinXY, kGClastXBinXY, kGCnYBinsXY, kGCfirstYBinXY, kGClastYBinXY, "", "");}
    if(kGCplotESDNoCutConversionOpeningAngle == kTRUE){ histograms->AddHistogram("ESD_NoCutConversion_OpeningAngle" ,"" , kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");}

    if(kGCplotESDNoCutConvGammaCosPointingAngle == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_CosPointingAngle" ,"" , kGCnXBinsCosPointingAngle, kGCfirstXBinCosPointingAngle, kGClastXBinCosPointingAngle, "", "");}
    if(kGCplotESDNoCutConvGammaDcaDaugthers == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_DcaDaughters" ,"" , kGCnXBinsDcaDaughters, kGCfirstXBinDcaDaughters, kGClastXBinDcaDaughters, "", "");}
    if(kGCplotESDNoCutConvGammaNormDcaDistDaugthers == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_NormDcaDistDaughters" ,"" , kGCnXBinsNormDcaDistDaughters, kGCfirstXBinNormDcaDistDaughters, kGClastXBinNormDcaDistDaughters, "", "");}
    if(kGCplotESDNoCutConvGammaLikelihoodAP == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_LikelihoodAP" ,"" , kGCnXBinsLikelihoodAP, kGCfirstXBinLikelihoodAP, kGClastXBinLikelihoodAP, "", "");}
    if(kGCplotESDNoCutConvGammaEAsymmetryP== kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_E_AsymmetryP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsAsymmetry, kGCfirstYBinAsymmetry, kGClastYBinAsymmetry,"", "");}
    if(kGCplotESDNoCutConvGammaPAsymmetryP== kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_P_AsymmetryP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsAsymmetry, kGCfirstYBinAsymmetry, kGClastYBinAsymmetry,"", "");}


    if(kGCplotESDNoCutConvGammaEdEdxP== kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_E_dEdxP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"", "",0);}
    if(kGCplotESDNoCutConvGammaPdEdxP== kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_P_dEdxP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"", "",0);}

    if(kGCplotESDNoCutConvGammaMCPtEta == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_MC_Pt_Eta" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
    if(kGCplotESDNoCutConversionMCZR == kTRUE){ histograms->AddHistogram("ESD_NoCutConversion_MC_ZR" ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");}
    if(kGCplotESDNoCutConversionMCXY == kTRUE){ histograms->AddHistogram("ESD_NoCutConversion_MC_XY" ,"" , kGCnXBinsXY, kGCfirstXBinXY, kGClastXBinXY, kGCnYBinsXY, kGCfirstYBinXY, kGClastYBinXY, "", "");}
		
		
    if(kGCplotESDMotherChi2 == kTRUE){ histograms->AddHistogram("ESD_Mother_Chi2","" , kGCnXBinsMesonChi2, kGCfirstXBinGammaChi2, kGClastXBinGammaChi2, "", "");}
    if(kGCplotESDMotherOpeningAngleGamma == kTRUE){ histograms->AddHistogram("ESD_Mother_GammaDaughter_OpeningAngle" ,"" , kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");}
    if(kGCplotESDMotherEnergy == kTRUE){ histograms->AddHistogram("ESD_Mother_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
    if(kGCplotESDMotherPt == kTRUE){ histograms->AddHistogram("ESD_Mother_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
    if(kGCplotESDMotherEta == kTRUE){ histograms->AddHistogram("ESD_Mother_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
    if(kGCplotESDMotherPhi == kTRUE){ histograms->AddHistogram("ESD_Mother_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
    if(kGCplotESDMotherMass == kTRUE){ histograms->AddHistogram("ESD_Mother_Mass" ,"" , kGCnXBinsPi0Mass, kGCfirstXBinPi0Mass, kGClastXBinPi0Mass, "", "");}
    if(kGCplotESDMotherR == kTRUE){ histograms->AddHistogram("ESD_Mother_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
    if(kGCplotESDMotherZR == kTRUE){ histograms->AddHistogram("ESD_Mother_ZR" ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");}
    if(kGCplotESDMotherXY == kTRUE){ histograms->AddHistogram("ESD_Mother_XY" ,"" , kGCnXBinsXY, kGCfirstXBinXY, kGClastXBinXY, kGCnYBinsXY, kGCfirstYBinXY, kGClastYBinXY, "", "");}
    if(kGCplotESDMotherRapid == kTRUE){ histograms->AddHistogram("ESD_Mother_Rapid" ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
		
    for(Int_t z=0;z<8;z++){
      for(Int_t m=0;m<6;m++){
	if(kGCplotESDBackgroundOpeningAngleGamma == kTRUE){ histograms->AddHistogram(Form("%d%dESD_Background_GammaDaughter_OpeningAngle",z,m) ,"" , kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");}
	if(kGCplotESDBackgroundEnergy == kTRUE){ histograms->AddHistogram(Form("%d%dESD_Background_Energy",z,m) ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
	if(kGCplotESDBackgroundPt == kTRUE){ histograms->AddHistogram(Form("%d%dESD_Background_Pt",z,m) ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
	if(kGCplotESDBackgroundEta == kTRUE){ histograms->AddHistogram(Form("%d%dESD_Background_Eta",z,m) ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
	if(kGCplotESDBackgroundPhi == kTRUE){ histograms->AddHistogram(Form("%d%dESD_Background_Phi",z,m) ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
	if(kGCplotESDBackgroundMass == kTRUE){ histograms->AddHistogram(Form("%d%dESD_Background_Mass",z,m) ,"" , kGCnXBinsEtaMass, kGCfirstXBinEtaMass, kGClastXBinEtaMass, "", "");}
	if(kGCplotESDBackgroundR == kTRUE){ histograms->AddHistogram(Form("%d%dESD_Background_R",z,m) ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
	if(kGCplotESDBackgroundZR == kTRUE){ histograms->AddHistogram(Form("%d%dESD_Background_ZR",z,m) ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");}
	if(kGCplotESDBackgroundXY == kTRUE){ histograms->AddHistogram(Form("%d%dESD_Background_XY",z,m) ,"" , kGCnXBinsXY, kGCfirstXBinXY, kGClastXBinXY, kGCnYBinsXY, kGCfirstYBinXY, kGClastYBinXY, "", "");}
	if(kGCplotESDBackgroundRapid == kTRUE){ histograms->AddHistogram(Form("%d%dESD_Background_Rapid",z,m) ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
      }
    }

    if(kGCplotESDBackgroundOpeningAngleGamma == kTRUE){ histograms->AddHistogram("ESD_Background_GammaDaughter_OpeningAngle" ,"" , kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");}
    if(kGCplotESDBackgroundEnergy == kTRUE){ histograms->AddHistogram("ESD_Background_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
    if(kGCplotESDBackgroundPt == kTRUE){ histograms->AddHistogram("ESD_Background_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
    if(kGCplotESDBackgroundEta == kTRUE){ histograms->AddHistogram("ESD_Background_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
    if(kGCplotESDBackgroundPhi == kTRUE){ histograms->AddHistogram("ESD_Background_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
    if(kGCplotESDBackgroundMass == kTRUE){ histograms->AddHistogram("ESD_Background_Mass" ,"" , kGCnXBinsEtaMass, kGCfirstXBinEtaMass, kGClastXBinEtaMass, "", "");}
    if(kGCplotESDBackgroundR == kTRUE){ histograms->AddHistogram("ESD_Background_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
    if(kGCplotESDBackgroundZR == kTRUE){ histograms->AddHistogram("ESD_Background_ZR" ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");}
    if(kGCplotESDBackgroundXY == kTRUE){ histograms->AddHistogram("ESD_Background_XY" ,"" , kGCnXBinsXY, kGCfirstXBinXY, kGClastXBinXY, kGCnYBinsXY, kGCfirstYBinXY, kGClastYBinXY, "", "");}
    if(kGCplotESDBackgroundRapid == kTRUE){ histograms->AddHistogram("ESD_Background_Rapid" ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
		
		
    if(kGCplotMapping == kTRUE){
      histograms->InitializeMappingValues(kGCnPhiIndex,kGCnRIndex,kGCnXBinsMapping,kGCminRadius,kGCmaxRadius,kGCnYBinsMapping,kGCminPhi,kGCmaxPhi);
      histograms->AddMappingHistograms(kGCnPhiIndex,kGCnRIndex,kGCnXBinsMapping,kGCfirstXBinMapping,kGClastXBinMapping,kGCnYBinsMapping,kGCfirstYBinMapping,kGClastYBinMapping);
      //      histograms->AddMappingHistograms(kGCnPhiIndex,kGCnRIndex,kGCnXBinsMapping,kGCminRadius,kGCmaxRadius,kGCnYBinsMapping,kGCminPhi,kGCmaxPhi);
    }

    //
    //************************************* Defining Resolution histograms *******************************************************/
    //
    // written by Friederike Bock 
    // contact: Friederike.Bock@cern.ch
    //

    if(kGCrunRES == kTRUE){
	//------------------------------------------ Absolute Resolutions --------------------------------------------------------
    if(kGCplotResolutiondRAbs== kTRUE){
		histograms->AddHistogram("Resolution_dRAbs_VS_R","" ,kGCnXBinsResdR, kGCfirstXBinResdR, kGClastXBinResdR,kGCnYBinsResdR,kGCfirstYBinResdR, kGClastYBinResdR, "", "");}
    if(kGCplotResolutiondZAbs== kTRUE){
		histograms->AddHistogram("Resolution_dZAbs_VS_R","" ,kGCnXBinsResdR, kGCfirstXBinResdR, kGClastXBinResdR,kGCnYBinsResdR,kGCfirstYBinResdR, kGClastYBinResdR, "", "");}
    if(kGCplotResolutiondPhiAbs== kTRUE){
		histograms->AddHistogram("Resolution_dPhiAbs_VS_R","" ,kGCnXBinsResdR, kGCfirstXBinResdR, kGClastXBinResdR,kGCnYBinsResdR, -TMath::Pi()/30., TMath::Pi()/30., "", "");}

	//------------------------------------------ Relative Resolutions --------------------------------------------------------
    if(kGCplotResolutiondR == kTRUE){
		histograms->AddHistogram("Resolution_dR" ,"" , kGCnXBinsResdR, kGCfirstXBinResdR, kGClastXBinResdR, kGCnYBinsResdR, kGCfirstYBinResdR, kGClastYBinResdR, "", "");}
    if(kGCplotResolutiondZ == kTRUE){
		histograms->AddHistogram("Resolution_dZ" ,"" , kGCnXBinsResdZ, kGCfirstXBinResdZ, kGClastXBinResdZ, kGCnYBinsResdZ, kGCfirstYBinResdZ, kGClastYBinResdZ, "", "");}

	//------------------------------------------- Pt vs R ---------------------------------------------------------------------		
    if(kGCplotResolutiondRdPt == kTRUE){
		histograms->AddHistogram("Resolution_R_dPt" ,"" , kGCnXBinsResdR, kGCfirstXBinResdR, kGClastXBinResdR, kGCnYBinsResdRdPt, kGCfirstYBinResdRdPt, kGClastYBinResdRdPt, "", "");}


	// ------------------------------------------- Reconstruction Plots for Resolution ----------------------------------------		
    if(kGCplotResolutionMCPt == kTRUE){
		histograms->AddHistogram("Resolution_MC_Pt" ,"" , kGCnXBinsResPt, kGCfirstXBinResPt, kGClastXBinResPt,"","");}
    if(kGCplotResolutionMCR == kTRUE){
		histograms->AddHistogram("Resolution_MC_R" ,"" , kGCnXBinsResR, kGCfirstXBinResR, kGClastXBinResR,"","");}
    if(kGCplotResolutionMCZ == kTRUE){
		histograms->AddHistogram("Resolution_MC_Z" ,"" , kGCnXBinsResZ, kGCfirstXBinResZ, kGClastXBinResZ,"","");}
		
    if(kGCplotResolutionESDPt == kTRUE){
		histograms->AddHistogram("Resolution_ESD_Pt" ,"" , kGCnXBinsResPt, kGCfirstXBinResPt, kGClastXBinResPt,"","");}
    if(kGCplotResolutionESDR == kTRUE){
		histograms->AddHistogram("Resolution_ESD_R" ,"" , kGCnXBinsResR, kGCfirstXBinResR, kGClastXBinResR,"","");}
    if(kGCplotResolutionESDZ == kTRUE){
		histograms->AddHistogram("Resolution_ESD_Z" ,"" , kGCnXBinsResZ, kGCfirstXBinResZ, kGClastXBinResZ,"","");}

	// ------------------------------------------- Plots for specific Gamma Trigger Studies -----------------------------------	
    if(kGCplotResolutionPtdPt = kTRUE){
		// ::::::::::::::::::::::::::::::::::::::: histograms for gammas ::::::::::::::::::::::::::::::::::::::::::::::::::::::
		histograms->AddHistogram("Resolution_Gamma_dPt_Pt" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_Gamma_dPt_Phi" ,"" , kGCnYBinsResdR, -TMath::Pi(), TMath::Pi(), kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		
		// ::::::::::::::::::::::::::::::::::::::: histograms for electrons :::::::::::::::::::::::::::::::::::::::::::::::::::
		histograms->AddHistogram("Resolution_E_dPt_Pt" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_E_dPt_Pt_ITS0" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_E_dPt_Pt_ITS1" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_E_dPt_Pt_ITS2" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_E_dPt_Pt_ITS3" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_E_dPt_Pt_ITS4" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_E_dPt_Pt_ITS5" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_E_dPt_Pt_ITS6" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_E_dPt_Phi" ,"" , kGCnYBinsResdR, -TMath::Pi(), TMath::Pi(), kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_E_nTRDtracklets_ESDPt" ,"" ,kGCnXBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, 8, -0.5, 7.5,"", "");
		histograms->AddHistogram("Resolution_E_nTRDtracklets_MCPt","" ,kGCnXBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, 8, -0.5, 7.5,"", "");	
		histograms->AddHistogram("Resolution_E_nTRDclusters_ESDPt","",kGCnXBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, 201, -0.5, 200.5,"", "");
		histograms->AddHistogram("Resolution_E_nTRDclusters_MCPt","",kGCnXBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, 201, -0.5, 200.5,"", "");
		//		histograms->AddHistogram("Resolution_E_TRDsignal_ESDPt","", fV0Reader->GetNegativeTrackPt(), fV0Reader->GetNegativeESDTrack()->GetTRDsignal());
		
		// :::::::::::::::::::::::::::::::::::::::: histograms for positrons :::::::::::::::::::::::::::::::::::::::::::::::::::
		histograms->AddHistogram("Resolution_P_dPt_Pt" ,"" , kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_P_dPt_Pt_ITS0" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_P_dPt_Pt_ITS1" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_P_dPt_Pt_ITS2" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_P_dPt_Pt_ITS3" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_P_dPt_Pt_ITS4" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_P_dPt_Pt_ITS5" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_P_dPt_Pt_ITS6" ,"" ,kGCnYBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_P_dPt_Phi" ,"" , kGCnYBinsResdR, -TMath::Pi(), TMath::Pi(), kGCnYBinsResdPt, kGCfirstYBinResdPt, kGClastYBinResdPt, "", "");
		histograms->AddHistogram("Resolution_P_nTRDtracklets_ESDPt" ,"" ,kGCnXBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, 8, -0.5, 7.5,"", "");   
		histograms->AddHistogram("Resolution_P_nTRDtracklets_MCPt","", kGCnXBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, 8, -0.5, 7.5,"", "");
		histograms->AddHistogram("Resolution_P_nTRDclusters_ESDPt","",kGCnXBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, 201, -0.5, 200.5,"", "");
		histograms->AddHistogram("Resolution_P_nTRDclusters_MCPt","",kGCnXBinsResdPt, kGCfirstXBinResdPt, kGClastXBinResdPt, 201,-0.5, 200.5,"", "");
		//		histograms->AddHistogram("Resolution_P_TRDsignal_ESDPt", "",fV0Reader->GetPositiveTrackPt(), fV0Reader->GetPositiveESDTrack()->GetTRDsignal());
    } //end of specific trigger study resolution plots
    } //end if(kGCrunRES=true)
    
    // ___________________________________________________________________________________________________________________________________________________

    if(kGCplotEventQuality == kTRUE){histograms->AddHistogram("ESD_EventQuality","ESD_EventQuality",kGCnXBinsEvtQ,kGCfirstXBinEvtQ,kGClastXBinEvtQ,"","");}
    if(kGCplotESDNumberOfV0s == kTRUE){histograms->AddHistogram("ESD_NumberOfV0s","Number of v0s",kGCnXBinsNV0,kGCfirstXBinNV0 ,kGClastXBinNV0 ,"","");}
    if(kGCplotESDNumberOfSurvivingV0s == kTRUE){histograms->AddHistogram("ESD_NumberOfSurvivingV0s","Number of surviving v0s",kGCnXBinsNV0, kGCfirstXBinNV0 , kGClastXBinNV0,"","");}
    if(kGCplotESDNumberOfContributorsVtx == kTRUE){histograms->AddHistogram("ESD_NumberOfContributorsVtx","Number of contributors to vertex",kGCnXBinsESDtrk, kGCfirstXBinESDtrk, kGClastXBinESDtrk,"","");}
    if(kGCplotESDNumberOfGoodESDTracks == kTRUE){histograms->AddHistogram("ESD_NumberOfGoodESDTracks","Number of Good ESD tracks",kGCnXBinsESDtrk, kGCfirstXBinESDtrk, kGClastXBinESDtrk,"","");}
    if(kGCplotESDNumberOfGoodESDTracks == kTRUE){histograms->AddHistogram("ESD_NumberOfGoodESDTracksVtx","Number of Good ESD tracks",kGCnXBinsESDtrk, kGCfirstXBinESDtrk, kGClastXBinESDtrk,"","");}	
	
    //  debug histograms
    if(kGCplotESDCutGetOnFly == kTRUE){histograms->AddHistogram("ESD_CutGetOnFly_InvMass" ,"Not GetOnFly" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutNContributors == kTRUE){histograms->AddHistogram("ESD_CutNContributors_InvMass" ,"NContributors <= 0" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutLikeSign == kTRUE){histograms->AddHistogram("ESD_CutLikeSign_InvMass" ,"LikeSign" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutRefit == kTRUE){histograms->AddHistogram("ESD_CutRefit_InvMass" ,"No TPC refit" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutKink == kTRUE){histograms->AddHistogram("ESD_CutKink_InvMass" ,"Kinks" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutPIDProb == kTRUE){histograms->AddHistogram("ESD_CutPIDProb_InvMass" ,"wrong TPC PID" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}

    if(kGCplotESDCutdedxSigmaElectronLine == kTRUE){histograms->AddHistogram("ESD_CutdEdxSigmaElectronLine_InvMass" ,"dedx ElectronLine" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutdedxSigmaPionLine == kTRUE){histograms->AddHistogram("ESD_CutdEdxSigmaPionLine_InvMass" ,"dedx PionLine" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutPionRejectionLowP==kTRUE){histograms->AddHistogram("ESD_CutPionRejectionLowP_InvMass" ,"dedx PionRejection LowP" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutKaonRejectionLowP==kTRUE){histograms->AddHistogram("ESD_CutKaonRejectionLowP_InvMass" ,"dedx KaonRejection LowP" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutQtGammaSelection==kTRUE){histograms->AddHistogram("ESD_CutQt_InvMass","ESD_CutQt_InvMass",kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    histograms->AddHistogram("ESD_CutTOFsigmaElec_InvMass", "ESD_CutTOFsigmaElec_InvMass",kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"",""); // RRnewTOF

    if(kGCplotESDCutProtonRejectionLowP==kTRUE){histograms->AddHistogram("ESD_CutProtonRejectionLowP_InvMass" ,"dedx ProtonRejection LowP" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutR == kTRUE){histograms->AddHistogram("ESD_CutR_InvMass" ,"Above RMax" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutMinR == kTRUE){histograms->AddHistogram("ESD_CutMinR_InvMass" ,"Above RMax" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutNDF == kTRUE){histograms->AddHistogram("ESD_CutNDF_InvMass" ,"NDF <= 0" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutChi2 == kTRUE){histograms->AddHistogram("ESD_CutChi2_InvMass" ,"#chi^{2} > Max" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutEta == kTRUE){histograms->AddHistogram("ESD_CutEta_InvMass" ,"Above #eta max" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutSinglePt == kTRUE){histograms->AddHistogram("ESD_CutSinglePt_InvMass" ,"Below p_{t} min" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutPt == kTRUE){histograms->AddHistogram("ESD_CutPt_InvMass" ,"Below p_{t} min" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutLine == kTRUE){histograms->AddHistogram("ESD_CutLine_InvMass" ,"Out of reconstruction area" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutZ == kTRUE){histograms->AddHistogram("ESD_CutZ_InvMass" ,"Out of reconstruction area" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutMinClsTPC == kTRUE){histograms->AddHistogram("ESD_CutMinNClsTPC_InvMass" ,"Out of reconstruction area" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDCutMinClsTPCToF == kTRUE){histograms->AddHistogram("ESD_CutMinNClsTPCToF_InvMass" ,"Out of reconstruction area" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}

    if(kGCplotESDGoodV0s == kTRUE){histograms->AddHistogram("ESD_GoodV0s_InvMass" ,"Good V0s" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDAllV0s == kTRUE){histograms->AddHistogram("ESD_AllV0s_InvMass" ,"All V0s" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}
    if(kGCplotESDAllV0sCurrentFinder == kTRUE){histograms->AddHistogram("ESD_AllV0sCurrentFinder_InvMass" ,"All V0s Current Finder" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");}

    if(kGCplotESDAllV0sCurrentFinderQtAlfa== kTRUE){ histograms->AddHistogram("ESD_AllV0sCurrentFinder_alfa_qt" ,"" ,kGCnXBinsP, kGCfirstXBinAlphaG, kGClastXBinAlpha,kGCnYBinsQt, kGCfirstYBinQt, kGClastYBinQt,"", "");}
    if(kGCplotESDAllV0sCurrentFinderQtAlfa== kTRUE){ histograms->AddHistogram("ESD_AllV0sCurrentFinder_goodtracks_alfa_qt" ,"" ,kGCnXBinsP, kGCfirstXBinAlphaG, kGClastXBinAlpha,kGCnYBinsQt, kGCfirstYBinQt, kGClastYBinQt,"", "");}

    if(kGCplotESDTrueConvGammaTrackLength == kTRUE){histograms->AddHistogram("ESD_TrueConvGamma_TrackLength","Track length of TrueConvGamma",kGCnXBinsTrackLength,kGCfirstXBinTrackLength,kGClastXBinTrackLength,"","");}
    if(kGCplotESDTrueConvGammaTrackLengthVSInvMass == kTRUE){histograms->AddHistogram("ESD_TrueConvGamma_TrackLengthVSInvMass","Track length of TrueConvGamma vs Inv mass",kGCnXBinsTrackLength,kGCfirstXBinTrackLength,kGClastXBinTrackLength,kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt,"","");}
		
    if(kGCplotOmegaSpectra == kTRUE){
      histograms->AddHistogram("ESD_Omega_InvMass_vs_Pt" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
      histograms->AddHistogram("ESD_Omega_InvMass","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");

      histograms->AddHistogram("ESD_Omega_Bck_InvMass_vs_Pt" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
      histograms->AddHistogram("ESD_Omega_Bck_InvMass","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
      histograms->AddHistogram("ESD_OmegaPipPinPi0_InvMass_vs_Pt" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
      histograms->AddHistogram("ESD_OmegaPipPinPi0_InvMass","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");

    }

    if(kGCplotPi0Spectra == kTRUE){
      histograms->AddHistogram("ESD_Mother_alfa","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinAlpha, kGClastXBinAlpha,"#alpha","Counts");
      histograms->AddHistogram("ESD_Mother_alfa_Pi0","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinAlpha, kGClastXBinAlpha,"#alpha","Counts");
      histograms->AddHistogram("ESD_Mother_alfa_Eta","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinAlpha, kGClastXBinAlpha,"#alpha","Counts");
      histograms->AddHistogram("ESD_Background_alfa_Pi0","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinAlpha, kGClastXBinAlpha,"#alpha","Counts");
      histograms->AddHistogram("ESD_Background_alfa_Eta","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinAlpha, kGClastXBinAlpha,"#alpha","Counts");
      //      histograms->AddHistogram("ESD_Mother_InvMass_vs_Pt" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
      histograms->AddHistogram("ESD_Mother_InvMass_vs_Pt" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
      histograms->AddHistogram("ESD_Mother_InvMass_vs_Pt_alpha" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
      histograms->AddHistogram("ESD_Mother_InvMass_vs_E_alpha" ,"Invariant Mass vs E" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","E [GeV]");
      histograms->AddHistogram("ESD_Mother_InvMass","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
      histograms->AddHistogram("ESD_Mother_InvMass_1212","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
      histograms->AddHistogram("ESD_Mother_InvMass_0912","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
      histograms->AddHistogram("ESD_Mother_InvMass_0909","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
      histograms->AddHistogram("ESD_Mother_InvMass_vs_Pt1212" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
      histograms->AddHistogram("ESD_Mother_InvMass_vs_Pt0912" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
      histograms->AddHistogram("ESD_Mother_InvMass_vs_Pt0909" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");

      if(kGCLowMemoryConsumption == kFALSE){
	histograms->AddHistogram("ESD_Mother_InvMass_GammaConvPHOS","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
	histograms->AddHistogram("ESD_Mother_InvMass_vs_Pt_GammaConvPHOS" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
	histograms->AddHistogram("ESD_Mother_InvMass_GammaConvPHOS_OpanLow","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
	histograms->AddHistogram("ESD_Mother_InvMass_GammaConvPHOS_OpanHigh","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");

	histograms->AddHistogram("ESD_Mother_InvMass_GammaConvEMCAL","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
	histograms->AddHistogram("ESD_Mother_InvMass_GammaConvEMCAL_Bck","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
	
	histograms->AddHistogram("ESD_Mother_InvMass_vs_Pt_GammaConvEMCAL" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
	histograms->AddHistogram("ESD_Mother_InvMass_vs_Pt_GammaConvEMCAL_Bck" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
	
	histograms->AddHistogram("ESD_Mother_InvMass_GammaConvEMCAL_OpanLow","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
	histograms->AddHistogram("ESD_Mother_InvMass_GammaConvEMCAL_OpanHigh","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
      }

      //      if(kGCdoNeutralMesonV0MCCheck == kTRUE){
      histograms->AddHistogram("ESD_TrueBckGG_InvMass_vs_Pt","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Counts");
      histograms->AddHistogram("ESD_TrueBckCont_InvMass_vs_Pt","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Counts");
      histograms->AddHistogram("ESD_TruePi0Sec_InvMass_vs_Pt","Invariant mass vs Pt",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Counts");
      histograms->AddHistogram("ESD_TruePi0DalitzCont_InvMass_vs_Pt","Invariant mass vs Pt",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Counts");
	histograms->AddHistogram("ESD_TruePi0_InvMass","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
	histograms->AddHistogram("ESD_TruePi0_InvMass_1212","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
	histograms->AddHistogram("ESD_TruePi0_InvMass_0912","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
	histograms->AddHistogram("ESD_TruePi0_InvMass_0909","Invariant mass",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
	histograms->AddHistogram("ESD_TruePi0_OpeningAngle_1212" ,"" , kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");
	histograms->AddHistogram("ESD_TruePi0_OpeningAngle_0912" ,"" , kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");
	histograms->AddHistogram("ESD_TruePi0_OpeningAngle_0909" ,"" , kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");
	histograms->AddHistogram("ESD_TruePi0_InvMass_vs_Pt" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
histograms->AddHistogram("ESD_TruePi0_InvMass_vs_Pt_alpha" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");

 if(kGCLowMemoryConsumption == kFALSE){
   histograms->AddHistogram("ESD_TruePi0_InvMass_vs_Pt1212" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
   histograms->AddHistogram("ESD_TruePi0_InvMass_vs_Pt0912" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
   histograms->AddHistogram("ESD_TruePi0_InvMass_vs_Pt0909" ,"Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");

	//}

   histograms->AddHistogram("ESD_Mother_InvMass_vs_Pt_Fiducial" ,"Invariant Mass vs Pt |eta|<0.9" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
 }
      histograms->AddHistogram("ESD_Mother_InvMass_Fiducial","Invariant mass |eta|<0.9",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");


    }

    if(kGCplotPi0Spectra == kTRUE && kGCcalculateBackground == kTRUE){
      for(Int_t m=0;m<6;m++){
	histograms->AddHistogram(Form("%dESD_Mother_InvMass_vs_Pt",m) ,"Background Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
      }
      for(Int_t z=0;z<8;z++){
	for(Int_t m=0;m<6;m++){
	  histograms->AddHistogram(Form("%d%dESD_Mother_InvMass_vs_Pt",z,m) ,"ESD Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");

	  histograms->AddHistogram(Form("%d%dESD_Background_InvMass_vs_Pt",z,m) ,"Background Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");

	  
	  histograms->AddHistogram(Form("%d%dESD_Background_InvMass",z,m),"Invariant mass background",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass BG [GeV]","Counts");
	  histograms->AddHistogram(Form("%d%dESD_Mother_InvMass",z,m),"Invariant mass in z and m bins",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass BG [GeV]","Counts");


	  histograms->AddHistogram(Form("%d%dESD_Background_InvMassvsPtFid",z,m) ,"Background Invariant Mass vs Pt |eta|<0.9" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");

	 
	  histograms->AddHistogram(Form("%d%dESD_Background_InvMass_Fiducial",z,m),"Invariant mass background |eta|<0.9",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass BG [GeV]","Counts");
	}
      }
    
      histograms->AddHistogram("ESD_Background_InvMass_vs_Pt" ,"Background Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
      histograms->AddHistogram("ESD_Background_InvMass_vs_Pt_alpha" ,"Background Invariant Mass vs Pt" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");

     histograms->AddHistogram("ESD_Background_InvMass_vs_E_alpha" ,"Background Invariant Mass vs E" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","E [GeV]");

  histograms->AddHistogram("ESD_Background_InvMass","Invariant mass background",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass BG [GeV]","Counts");

      histograms->AddHistogram("ESD_Background_InvMass_vs_Pt_Fiducial" ,"Background Invariant Mass vs Pt |eta|<0.9" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
      histograms->AddHistogram("ESD_Background_InvMass_Fiducial","Invariant mass background |eta|<0.9",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass BG [GeV]","Counts");
    
    }
    
    if(kGCdoMCTruth){
      if(kGCplotMCConversionR == kTRUE){ histograms->AddHistogram("MC_Conversion_R","Radius of gamma conversion points",kGCnXBinsR, kGCfirstXBinR, kGClastXBinR,"counts","cm");}
      if(kGCplotMCConversionZR == kTRUE){ histograms->AddHistogram("MC_Conversion_ZR","Radius of gamma conversion points vs Z",kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "cm", "cm");}
      if(kGCplotMCConversionXY == kTRUE){ histograms->AddHistogram("MC_Conversion_XY","Gamma XY converison point.",kGCnXBinsXY, kGCfirstXBinXY, kGClastXBinXY, kGCnYBinsXY, kGCfirstYBinXY, kGClastYBinXY, "cm", "cm");}
      if(kGCplotMCConversionOpeningAngle == kTRUE){ histograms->AddHistogram("MC_Conversion_OpeningAngle","Opening angle of e+e- pairs from gamma conversion",kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "counts", "cm");}
      if(kGCplotMCConvGammaEAsymmetryP== kTRUE){ histograms->AddHistogram("MC_ConvGamma_E_AsymmetryP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsAsymmetry, kGCfirstYBinAsymmetry, kGClastYBinAsymmetry,"", "");}
      if(kGCplotMCConvGammaPAsymmetryP== kTRUE){ histograms->AddHistogram("MC_ConvGamma_P_AsymmetryP" ,"" ,kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsAsymmetry, kGCfirstYBinAsymmetry, kGClastYBinAsymmetry,"", "");}
		

      if(kGCplotMCEEnergy == kTRUE){ histograms->AddHistogram("MC_E_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
      if(kGCplotMCEPt == kTRUE){ histograms->AddHistogram("MC_E_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCEEta == kTRUE){ histograms->AddHistogram("MC_E_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCEPhi == kTRUE){ histograms->AddHistogram("MC_E_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
      if(kGCplotMCENTPCClusters == kTRUE){ histograms->AddHistogram("MC_E_nTPCClusters" ,"" , kGCnXBinsNTPCClusters, kGCfirstXBinNTPCClusters, kGClastXBinNTPCClusters, "", "");}
      if(kGCplotMCENITSClusters == kTRUE){ histograms->AddHistogram("MC_E_nITSClusters" ,"" , kGCnXBinsNITSClusters, kGCfirstXBinNITSClusters, kGClastXBinNITSClusters, "", "");}
		
      if(kGCplotMCPEnergy == kTRUE){ histograms->AddHistogram("MC_P_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
      if(kGCplotMCPPt == kTRUE){ histograms->AddHistogram("MC_P_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPEta == kTRUE){ histograms->AddHistogram("MC_P_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCPPhi == kTRUE){ histograms->AddHistogram("MC_P_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
      if(kGCplotMCPNTPCClusters == kTRUE){ histograms->AddHistogram("MC_P_nTPCClusters" ,"" , kGCnXBinsNTPCClusters, kGCfirstXBinNTPCClusters, kGClastXBinNTPCClusters, "", "");}
      if(kGCplotMCPNITSClusters == kTRUE){ histograms->AddHistogram("MC_P_nITSClusters" ,"" , kGCnXBinsNITSClusters, kGCfirstXBinNITSClusters, kGClastXBinNITSClusters, "", "");}
		
      if(kGCplotMCallGammaEnergy == kTRUE){ histograms->AddHistogram("MC_allGamma_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
      if(kGCplotMCallGammaPt == kTRUE){ histograms->AddHistogram("MC_allGamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCallGammaEta == kTRUE){ histograms->AddHistogram("MC_allGamma_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCallGammaPhi == kTRUE){ histograms->AddHistogram("MC_allGamma_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
      if(kGCplotMCallGammaRapid == kTRUE){ histograms->AddHistogram("MC_allGamma_Rapid" ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
		
      if(kGCplotMCConvGammaEnergy == kTRUE){ histograms->AddHistogram("MC_ConvGamma_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
      if(kGCplotMCConvGammaPt == kTRUE){ histograms->AddHistogram("MC_ConvGamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCConvGammaEta == kTRUE){ histograms->AddHistogram("MC_ConvGamma_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCConvGammaPhi == kTRUE){ histograms->AddHistogram("MC_ConvGamma_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
      if(kGCplotMCConvGammaRapid == kTRUE){ histograms->AddHistogram("MC_ConvGamma_Rapid" ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
      if(kGCplotMCConvGammaPtvsEta == kTRUE){ histograms->AddHistogram("MC_ConvGamma_Pt_Eta","", kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt,kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta,"","");}
		
      if(kGCplotMCallDirectGammaEnergy == kTRUE){ histograms->AddHistogram("MC_allDirectGamma_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
      if(kGCplotMCallDirectGammaPt == kTRUE){ histograms->AddHistogram("MC_allDirectGamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCallDirectGammaEta == kTRUE){ histograms->AddHistogram("MC_allDirectGamma_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCallDirectGammaPhi == kTRUE){ histograms->AddHistogram("MC_allDirectGamma_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
      if(kGCplotMCallDirectGammaRapid == kTRUE){ histograms->AddHistogram("MC_allDirectGamma_Rapid" ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
		
      if(kGCplotMCConvDirectGammaEnergy == kTRUE){ histograms->AddHistogram("MC_ConvDirectGamma_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
      if(kGCplotMCConvDirectGammaPt == kTRUE){ histograms->AddHistogram("MC_ConvDirectGamma_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCConvDirectGammaEta == kTRUE){ histograms->AddHistogram("MC_ConvDirectGamma_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCConvDirectGammaPhi == kTRUE){ histograms->AddHistogram("MC_ConvDirectGamma_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
      if(kGCplotMCConvDirectGammaRapid == kTRUE){ histograms->AddHistogram("MC_ConvDirectGamma_Rapid" ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
		
      if(kGCplotMCMotherEta == kTRUE){ histograms->AddHistogram("MC_Mother_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCMotherPhi == kTRUE){ histograms->AddHistogram("MC_Mother_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
      if(kGCplotMCMotherRapid == kTRUE){ histograms->AddHistogram("MC_Mother_Rapid" ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
      if(kGCplotMCMotherPt == kTRUE){ histograms->AddHistogram("MC_Mother_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCMotherEnergy == kTRUE){ histograms->AddHistogram("MC_Mother_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
      if(kGCplotMCMotherMass == kTRUE){ histograms->AddHistogram("MC_Mother_Mass" ,"" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass, "", "");}
      if(kGCplotMCMotherOpeningAngle == kTRUE){ histograms->AddHistogram("MC_Mother_GammaDaughter_OpeningAngle" ,"" , kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");}
      if(kGCplotMCMotherR == kTRUE){ histograms->AddHistogram("MC_Mother_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
      if(kGCplotMCMotherZR == kTRUE){ histograms->AddHistogram("MC_Mother_ZR" ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");}
      if(kGCplotMCMotherXY == kTRUE){ histograms->AddHistogram("MC_Mother_XY" ,"" , kGCnXBinsXY, kGCfirstXBinXY, kGClastXBinXY, kGCnYBinsXY, kGCfirstYBinXY, kGClastYBinXY, "", "");}
      if(kGCplotMCMotherPtvsEtaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Mother_Pt_Eta_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCMotherPtvsRapidWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Mother_Pt_Rapid_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
      if(kGCplotMCMotherPtvsEtaConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Mother_Pt_Eta_ConvGamma_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCMotherPtvsRapidConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Mother_Pt_Rapid_ConvGamma_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
		
      if(kGCplotMCMotherSpectra == kTRUE){ 
	histograms->AddHistogram("MC_Mother_InvMass_vs_Pt" ,"" ,kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra, kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra, "", "");
	histograms->AddHistogram("MC_Mother_InvMass_vs_Pt_withinAcceptance" ,"" ,kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra, kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra, "", "");
	histograms->AddHistogram("MC_Mother_InvMass_vs_Pt_ConvGamma_withinAcceptance" ,"" ,kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra, kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra, "", "");
      }
		
      if(kGCplotMCPhysicalPrimaryChargedPt == kTRUE){ histograms->AddHistogram("MC_PhysicalPrimaryCharged_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Eta == kTRUE){ histograms->AddHistogram("MC_Pi0_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}	
      if(kGCplotMCPi0Rapid == kTRUE){ histograms->AddHistogram("MC_Pi0_Rapid" ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}	
      if(kGCplotMCPi0PtvsRapid == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt_vs_Rapid" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
	if(kGCplotMCPi0PtvsRapid == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt_vs_Rapid_allDaughters" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");} // RR primary Pi0 debug
      if(kGCplotMCPi0Phi == kTRUE){ histograms->AddHistogram("MC_Pi0_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}

      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_SD_EvtQ1_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_SD_EvtQ2_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_SD_EvtQ3_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_SD_EvtQ4_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_SD_EvtQ5_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}

      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_DD_EvtQ1_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_DD_EvtQ2_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_DD_EvtQ3_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_DD_EvtQ4_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_DD_EvtQ5_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}


      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_ND_EvtQ1_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_ND_EvtQ2_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_ND_EvtQ3_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_ND_EvtQ4_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_ND_EvtQ5_Pi0_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}

      if(kGCplotMCPi0PtFiducial == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt_Fiducial" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0PtWithinAcceptanceFiducial == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt_withinAcceptance_Fiducial" ,"" , kGCnXBinsPt,kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0PtConvGammaWithinAcceptanceFiducial == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt_ConvGamma_withinAcceptance_Fiducial","" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0OpeningPtConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_ConvGamma_OpeningAngle_Pt","", kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");}
      if(kGCplotMCPi0PtGammaPtConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_ConvGamma_PtGamma_Pt","", kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0Energy == kTRUE){ histograms->AddHistogram("MC_Pi0_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
      if(kGCplotMCPi0Mass == kTRUE){ histograms->AddHistogram("MC_Pi0_Mass" ,"" , kGCnXBinsPi0Mass, kGCfirstXBinPi0Mass, kGClastXBinPi0Mass, "", "");}
      if(kGCplotMCPi0Alpha == kTRUE){ histograms->AddHistogram("MC_Pi0_alpha" ,"" , kGCnXBinsPi0Mass, kGCfirstXBinPi0Alpha, kGClastXBinPi0Alpha, "", "");}
      if(kGCplotMCEtaAlpha == kTRUE){ histograms->AddHistogram("MC_Eta_alpha" ,"" , kGCnXBinsPi0Mass, kGCfirstXBinPi0Alpha, kGClastXBinPi0Alpha, "", "");}
       if(kGCplotMCPi0OpeningAngle == kTRUE){ histograms->AddHistogram("MC_Pi0_GammaDaughter_OpeningAngle" ,"" , kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");}
      if(kGCplotMCPi0R == kTRUE){ histograms->AddHistogram("MC_Pi0_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
      if(kGCplotMCPi0ZR == kTRUE){ histograms->AddHistogram("MC_Pi0_ZR" ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");}
      if(kGCplotMCPi0XY == kTRUE){ histograms->AddHistogram("MC_Pi0_XY" ,"" , kGCnXBinsXY, kGCfirstXBinXY, kGClastXBinXY, kGCnYBinsXY, kGCfirstYBinXY, kGClastYBinXY, "", "");}
      if(kGCplotMCPi0PtvsEtaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt_Eta_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCPi0PtvsRapidWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt_Rapid_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
      if(kGCplotMCPi0PtvsEtaConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt_Eta_ConvGamma_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCPi0PtvsRapidConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt_Rapid_ConvGamma_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
      if(kGCplotMCPi0ZRConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_ZR_ConvGamma_withinAcceptance" ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");}
		
		
      if(kGCplotMCPi0SecondaryEta == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCPi0SecondaryRapid == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Rapid" ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
      if(kGCplotMCPi0SecondaryPhi == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
      if(kGCplotMCPi0SecondaryPt == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCPi0SecondaryEnergy == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
      if(kGCplotMCPi0SecondaryMass == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Mass" ,"" , kGCnXBinsPi0Mass, kGCfirstXBinPi0Mass, kGClastXBinPi0Mass, "", "");}
      if(kGCplotMCPi0SecondaryOpeningAngle == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_GammaDaughter_OpeningAngle" ,"" , kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");}
      if(kGCplotMCPi0SecondaryR == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
      if(kGCplotMCPi0SecondaryZR == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_ZR" ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");}
      if(kGCplotMCPi0SecondaryXY == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_XY" ,"" , kGCnXBinsXY, kGCfirstXBinXY, kGClastXBinXY, kGCnYBinsXY, kGCfirstYBinXY, kGClastYBinXY, "", "");}
      if(kGCplotMCPi0SecondaryPtvsEtaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Pt_Eta_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCPi0SecondaryPtvsRapidWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Pt_Rapid_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
      if(kGCplotMCPi0SecondaryPtvsEtaConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Pt_Eta_ConvGamma_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCPi0SecondaryPtvsRapidConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Pt_Rapid_ConvGamma_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
		
		
		
      if(kGCplotMCEtaEta == kTRUE){ histograms->AddHistogram("MC_Eta_Eta" ,"" , kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCEtaRapid == kTRUE){ histograms->AddHistogram("MC_Eta_Rapid" ,"" , kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
      if(kGCplotMCEtaPtvsRapid == kTRUE){ histograms->AddHistogram("MC_Eta_Pt_vs_Rapid" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
      if(kGCplotMCEtaPhi == kTRUE){ histograms->AddHistogram("MC_Eta_Phi" ,"" , kGCnXBinsPhi, kGCfirstXBinPhi, kGClastXBinPhi, "", "");}
      if(kGCplotMCEtaPt == kTRUE){ histograms->AddHistogram("MC_Eta_Pt" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
      if(kGCplotMCEtaEnergy == kTRUE){ histograms->AddHistogram("MC_Eta_Energy" ,"" , kGCnXBinsEnergy, kGCfirstXBinEnergy, kGClastXBinEnergy, "", "");}
      if(kGCplotMCEtaMass == kTRUE){ histograms->AddHistogram("MC_Eta_Mass" ,"" , kGCnXBinsEtaMass, kGCfirstXBinEtaMass, kGClastXBinEtaMass, "", "");}
      if(kGCplotMCEtaOpeningAngleGamma == kTRUE){ histograms->AddHistogram("MC_Eta_GammaDaughter_OpeningAngle" ,"" , kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");}
      if(kGCplotMCEtaR == kTRUE){ histograms->AddHistogram("MC_Eta_R" ,"" , kGCnXBinsR, kGCfirstXBinR, kGClastXBinR, "", "");}
      if(kGCplotMCEtaZR == kTRUE){ histograms->AddHistogram("MC_Eta_ZR" ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");}
      if(kGCplotMCEtaXY == kTRUE){ histograms->AddHistogram("MC_Eta_XY" ,"" , kGCnXBinsXY, kGCfirstXBinXY, kGClastXBinXY, kGCnYBinsXY, kGCfirstYBinXY, kGClastYBinXY, "", "");}
      if(kGCplotMCEtaPtvsEtaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Eta_Pt_Eta_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCEtaPtvsRapidWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Eta_Pt_Rapid_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
      if(kGCplotMCEtaPtvsEtaConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Eta_Pt_Eta_ConvGamma_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "", "");}
      if(kGCplotMCEtaPtvsRapidConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Eta_Pt_Rapid_ConvGamma_withinAcceptance" ,"" , kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, "", "");}
      if(kGCplotMCEtaZRConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Eta_ZR_ConvGamma_withinAcceptance" ,"" , kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "", "");} 
    }
    if(kGCplotMCEtaOpeningPtConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Eta_ConvGamma_OpeningAngle_Pt","", kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsOpeningAngle, kGCfirstXBinOpeningAngle, kGClastXBinOpeningAngle, "", "");}
    if(kGCplotMCEtaPtGammaPtConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Eta_ConvGamma_PtGamma_Pt","", kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");}
  }// end kGCrunNeutralMeson

  
  //---------------------------------------------------  2 gamma Background -------------------------------------------------------

  if(kGCcalculateBackground==kTRUE){
    histograms->AddHistogram("ESD_GlobalPrimaryVtxZ","Z primary vertex Global",300, -15., 15.,"Z[cm]","counts");
    histograms->AddHistogram("ESD_SPDPrimaryVtxZ","Z primary vertex SPD",300, -15., 15.,"Z[cm]","counts");
    histograms->AddHistogram("ESD_Z_distribution" ,"Z primary vertex" , 2000, -30, 30,"Z[cm]","counts");
    histograms->AddHistogram("ESD_multiplicity_distribution" ,"multiplicity distribution" , 200, 0, 200,"counts","Multiplicity");
    histograms->AddHistogram("ESD_ZvsMultiplicity" ,"Z vs Multiplicity" , 1000, -10, 10,200,0,200,"Z[cm]","Multiplicity");
  }
  
  if(kGCRecalculateV0ForGamma==kTRUE){
    histograms->AddHistogram("ESD_RecalculateV0_InvMass","", kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass,"","");
    histograms->AddHistogram("ESD_RecalculateV0_Pt","",kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt, "", "");
    histograms->AddHistogram("ESD_RecalculateV0_E_dEdxP","",kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"", "",0);
    histograms->AddHistogram("ESD_RecalculateV0_P_dEdxP","",kGCnXBinsP, kGCfirstXBinP, kGClastXBinP,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"", "",0);
    histograms->AddHistogram("ESD_RecalculateV0_cpa","",1000,-1.,1.,"","");
    histograms->AddHistogram("ESD_RecalculateV0_dca","",100,0.,10.,"","");
    histograms->AddHistogram("ESD_RecalculateV0_Rr","",250,0.,250.,"","");
    histograms->AddHistogram("ESD_RecalculateV0_normdistP","",250,0.,250.,"","");
    histograms->AddHistogram("ESD_RecalculateV0_normdistN","",250,0.,250.,"","");

    histograms->AddHistogram("ESD_RecalculateGG_InvMass","",kGCnXBinsSpectra,kGCfirstXBinSpectra, kGClastXBinSpectra,"InvMass [GeV]","Counts");
    histograms->AddHistogram("ESD_RecalculateGG_InvMass_vs_Pt","", kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,"InvMass [GeV]","Pt [GeV]");		

  }
  

}



Int_t SetAnalysisCutSelection(TString analysisCutSelection){
  Int_t iResult=0;
  
  // set the cuts depending on the Cut Selection Id
  // first number is dummy always set to 9 
  //  const char* cutSelection = analysisCutSelection.Data(); 
  if(analysisCutSelection.Length()!=c_array_size){
    cout<<"Cut selection has the wrong length!"<<endl;
    return 0;
  }

  char cutSelection[] = analysisCutSelection.Data();
  int array[c_array_size];
  const int N = sizeof(array) / sizeof(int);
  string2array( cutSelection, array );



  Int_t goodId=array[0];
  Int_t v0FinderType=array[1];
  Int_t eProbCut=array[2];
  Int_t ededxSigmaCut=array[3];
  Int_t pidedxSigmaCut=array[4];
  Int_t piMomdedxSigmaCut=array[5];
  Int_t chi2GammaCut=array[6];
  Int_t singlePtCut=array[7];
  Int_t clsTPCCut=array[8];
  Int_t etaCut=array[9];
  Int_t chi2MesonCut=array[10];
  Int_t LowPRejectionSigmaCut=array[11];
  Int_t QtMaxCut=array[12];
  Int_t piMaxMomdedxSigmaCut=array[13];
  Int_t alphaMesonCut=array[14];
  Int_t minRCut=array[15];
  Int_t RapidityMesonCut=array[16];
  Int_t BackgroundScheme=array[17];
  Int_t DegreesForRotationMethod=array[18];
  Int_t NumberOfRotations=array[19];
  Int_t removePileUp=array[20];
  Int_t selectV0AND=array[21];
  Int_t multiplicityBin=array[22];
  Int_t isHeavyIon=array[23];
  Int_t useCentrality=array[24];
  Int_t centralityBin=array[25];
  Int_t TOFelectronPID=array[26]; // RRnewTOF

  cout<<"CentralityBin::"<< centralityBin <<endl;
  cout<<"Use Centrality::"<< useCentrality <<endl;
  cout<<"Heavy Ion::"<< isHeavyIon<<endl;
  cout<<"Multiplicity Bin::"<< multiplicityBin<<endl;
  cout<<"Select V0AND::"<< selectV0AND<<endl;
  cout<<"Remove PileUp::"<< removePileUp<<endl;
  cout<<"NumberOfRotations::"<<NumberOfRotations<<endl;
  cout<<"DegreesForRotationMethod::"<<DegreesForRotationMethod<<endl;
  cout<<"BackgroundScheme::"<<BackgroundScheme<<endl;
  cout<<"RapidityMesonCut::"<<RapidityMesonCut<<endl;	
  cout<<"minRCut::"<<minRCut<<endl;
  cout<<"alphaMesonCut::"<<alphaMesonCut<<endl;	
  cout<<"piMaxMomdedxSigmaCut::"<<piMaxMomdedxSigmaCut<<endl;
  cout<<"QtMaxCut:"<<QtMaxCut<<endl;
  cout<<"LowPRejectionSigmaCut:"<<LowPRejectionSigmaCut<<endl;
  cout<<"chi2MesonCut: "<< chi2MesonCut<<endl;
  cout<<"etaCut: "<<etaCut<<endl;
  cout<<"clsTPCCut: "<<clsTPCCut<<endl;
  cout<<"singlePtCut: "<<singlePtCut<<endl;
  cout<<"chi2GammaCut: "<<chi2GammaCut<<endl;
  cout<<"piMomdedxSigmaCut: "<<piMomdedxSigmaCut<<endl;
  cout<<"pidedxSigmaCut: "<<pidedxSigmaCut <<endl;
  cout<<"ededxSigmaCut: "<<ededxSigmaCut <<endl;
  cout<<"eProbCut: "<< eProbCut<<endl;
  cout<<"v0FinderType: "<<v0FinderType <<endl;
  cout<<"goodId: "<<goodId <<endl;
  cout<<"TOFelectronPID: "<<TOFelectronPID<<endl; // RRnewTOF

  if(goodId !=9){
    cout<<"Analysis Cut Selection too short or does not start with 9"<<endl;
    return iResult;
  }

  switch (v0FinderType){
  case 0:  // on fly V0 finder
    kGCUseOnFlyV0Finder=kTRUE;
    break;
  case 1:  // offline V0 finder
    kGCUseOnFlyV0Finder=kFALSE;
    break;
  default:
    return iResult;
  }

  switch(eProbCut){
  case 0:
    kGCprobElectron=0;
    break;
  case 1:
    kGCprobElectron=0.1;
    break;
  case 2:
    kGCprobElectron=0.5;
    break;
  case 3:
    kGCprobElectron=0.7;
    break;
  default:
    return iResult;
  }

  switch(ededxSigmaCut){
  case 0: // -10,10
    kGCPIDnSigmaBelowElectronLine=-10;
    kGCPIDnSigmaAboveElectronLine=10;
    break;
  case 1: // -5,5 
    kGCPIDnSigmaBelowElectronLine=-5;
    kGCPIDnSigmaAboveElectronLine=5;
    break;
  case 2: // -3,5
    kGCPIDnSigmaBelowElectronLine=-3;
    kGCPIDnSigmaAboveElectronLine=5;
    break;
  case 3: // -4,5
    kGCPIDnSigmaBelowElectronLine=-4;
    kGCPIDnSigmaAboveElectronLine=5;
    break;
  default:
    return iResult;
  }
  
  switch(pidedxSigmaCut){
  case 0:  // -10
    kGCPIDnSigmaAbovePionLine=-10;
    break;
  case 1:   // 0
    kGCPIDnSigmaAbovePionLine=0;
    break;
  case 2:  // 1
    kGCPIDnSigmaAbovePionLine=1;
    break;
  case 3:  // 1
    kGCPIDnSigmaAbovePionLine=-1;
    break;
  case 4:  // 1
    kGCPIDnSigmaAbovePionLine=-1.5;
    break;
  case 5:  // 1
    kGCPIDnSigmaAbovePionLine=2.;
    break;

  default:
    return iResult;
  }
  
  switch(piMomdedxSigmaCut){
  case 0:  // 0.5 GeV
    kGCPIDMinPnSigmaAbovePionLine=0.5;
    break;
  case 1:  // 1. GeV
    kGCPIDMinPnSigmaAbovePionLine=1.;
    break;
  case 2:  // 1.5 GeV
    kGCPIDMinPnSigmaAbovePionLine=1.5;
    break;
  case 3:  // 20.0 GeV
    kGCPIDMinPnSigmaAbovePionLine=20.;
    break;
  case 4:  // 50.0 GeV
    kGCPIDMinPnSigmaAbovePionLine=50.;
    break;
  case 5:  // 0.3 GeV
    kGCPIDMinPnSigmaAbovePionLine=0.3;
    break;
  case 6:  // 0.25 GeV     // RRnew label was wrong in comment
    kGCPIDMinPnSigmaAbovePionLine=0.25;
    break;
  case 7:  // 0.4 GeV
    kGCPIDMinPnSigmaAbovePionLine=0.4;
    break;
  default:
    return iResult;
  }
  
  switch(chi2GammaCut){
  case 0: // 100
    kGCchi2CutConversion = 100.;
    break;
  case 1:  // 50
    kGCchi2CutConversion = 50.;
    break;
  case 2:  // 30
    kGCchi2CutConversion = 30.;
    break;
  case 3:
    kGCchi2CutConversion = 200.;
    break;
  case 4:
    kGCchi2CutConversion = 500.;
    break;
  case 5:
    kGCchi2CutConversion = 1000.;
    break;
  case 6:
    kGCchi2CutConversion = 5.;
    break;
  case 7:
    kGCchi2CutConversion = 10.;
    break;
  case 8:
    kGCchi2CutConversion = 20.;
    break;
  case 9:
    kGCchi2CutConversion = 15.;
    break;
  default:
    return iResult;
  }

  switch(singlePtCut){
  case 0: // 0.050 GeV
    kGCsingleptCut = 0.050;
    break;
  case 1:  // 0.100 GeV
    kGCsingleptCut = 0.100;
    break;
  case 2:  // 0.150 GeV
    kGCsingleptCut = 0.150;
    break;
  case 3:  // 0.200 GeV
    kGCsingleptCut = 0.200;
    break;
  case 4:  // 0.075 GeV
    kGCsingleptCut = 0.075;
    break;
  case 5:  // 0.125 GeV
    kGCsingleptCut = 0.125;
    break;
  default:
    return iResult;
 }

  switch(clsTPCCut){
  case 0: // 0 
    kGCminClsTPCCut= 0.;
    break;
  case 1:  // 70 
    kGCminClsTPCCut= 70.;
    break;
  case 2:  // 80 
    kGCminClsTPCCut= 80.;
    break;
  case 3:  // 100 
    kGCminClsTPCCut= 100.;
    break;
  case 4:  // 60% of findable clusters
    kGCminClsTPCCutToF= 0.6;
    kGCUseCorrectedTPCClsInfo=0;
    break;
  case 5:  // 0% of findable clusters
    kGCminClsTPCCutToF= 0.0;
    kGCUseCorrectedTPCClsInfo=1;
    break;
  case 6:  // 0% of findable clusters
    kGCminClsTPCCutToF= 0.7;
    kGCUseCorrectedTPCClsInfo=0;
    break;
  case 7:  // 0% of findable clusters
    kGCminClsTPCCutToF= 0.35;
    kGCUseCorrectedTPCClsInfo=0;
    break;
  case 8:
    kGCminClsTPCCutToF= 0.35;
    kGCUseCorrectedTPCClsInfo=1;
    break;
  case 9:
    kGCminClsTPCCutToF= 0.6;
    kGCUseCorrectedTPCClsInfo=1;
    break;
  default:
    return iResult;
  }

  switch(etaCut){
  case 0: // 0.9 
    kGCetaCut    = 0.9;
    kGCLineCutZRSlope = tan(2*atan(exp(-kGCetaCut)));
    break;
  case 1:  // 1.2
    kGCetaCut    = 1.2;
    kGCLineCutZRSlope = tan(2*atan(exp(-kGCetaCut)));
    break;
  case 2:  // 1.4
    kGCetaCut    = 1.4;
    kGCLineCutZRSlope = tan(2*atan(exp(-kGCetaCut)));
    break;
  default:
    return iResult;
  }

  switch(chi2MesonCut){
  case 0:  // 100.
    kGCchi2CutMeson = 100.;
    break;
  case 1:  // 50.
    kGCchi2CutMeson = 50.;
    break;
  case 2:  // 30.
    kGCchi2CutMeson = 30.;
    break;
  case 3:
    kGCchi2CutMeson = 200.;
    break;
  case 4:
    kGCchi2CutMeson = 500.;
    break;
  case 5:
    kGCchi2CutMeson = 1000.;
    break;
  default:
    return iResult;
  }


  switch(LowPRejectionSigmaCut){
  case 0:  // 
    kGCPIDnSigmaAtLowPAroundKaonLine=0;
    kGCPIDnSigmaAtLowPAroundProtonLine=0;
    kGCPIDnSigmaAtLowPAroundPionLine=0;
    break;
  case 1:  // 
    kGCPIDnSigmaAtLowPAroundKaonLine=0.5;
    kGCPIDnSigmaAtLowPAroundProtonLine=0.5;
    kGCPIDnSigmaAtLowPAroundPionLine=0.5;
    break;
  case 2:  // 
    kGCPIDnSigmaAtLowPAroundKaonLine=1;
    kGCPIDnSigmaAtLowPAroundProtonLine=1;
    kGCPIDnSigmaAtLowPAroundPionLine=1;
    break;
  case 3:  // 
    kGCPIDnSigmaAtLowPAroundKaonLine=2.;
    kGCPIDnSigmaAtLowPAroundProtonLine=2.;
    kGCPIDnSigmaAtLowPAroundPionLine=2.;
    break;
  case 4:  // 
    kGCPIDnSigmaAtLowPAroundKaonLine=0.;
    kGCPIDnSigmaAtLowPAroundProtonLine=0.;
    kGCPIDnSigmaAtLowPAroundPionLine=1;
    break;
  case 5:  // 
    kGCPIDnSigmaAtLowPAroundKaonLine=0.;
    kGCPIDnSigmaAtLowPAroundProtonLine=0.;
    kGCPIDnSigmaAtLowPAroundPionLine=1.5;
    break;
  case 6:  // 
    kGCPIDnSigmaAtLowPAroundKaonLine=0.;
    kGCPIDnSigmaAtLowPAroundProtonLine=0.;
    kGCPIDnSigmaAtLowPAroundPionLine=2.;
    break;
  default:
    return iResult;
  }
  switch(QtMaxCut){
  case 0: //
    kGCQtMax=1.;
    kGCdoHighPtQtGammaSelection=kFALSE; // RRnew
    kGCHighPtQtMax=100.;	        // RRnew
    kGCPtBorderForQt=100.;	        // RRnew
    break;
  case 1:
    kGCQtMax=0.1;
    kGCdoHighPtQtGammaSelection=kFALSE; // RRnew
    kGCHighPtQtMax=100.;	        // RRnew
    kGCPtBorderForQt=100.;	        // RRnew
    break;
  case 2:
    kGCQtMax=0.07;
    kGCdoHighPtQtGammaSelection=kFALSE; // RRnew
    kGCHighPtQtMax=100.;	        // RRnew
    kGCPtBorderForQt=100.;	        // RRnew
    break;
  case 3:
    kGCQtMax=0.05;
    kGCdoHighPtQtGammaSelection=kFALSE; // RRnew
    kGCHighPtQtMax=100.;	        // RRnew
    kGCPtBorderForQt=100.;	        // RRnew
    break;
  case 4:
    kGCQtMax=0.03;
    kGCdoHighPtQtGammaSelection=kFALSE; // RRnew
    kGCHighPtQtMax=100.;	        // RRnew
    kGCPtBorderForQt=100.;	        // RRnew
    break;
  case 5: // RR try to improve (get rid of) low InvMass peak in PbPb
    kGCQtMax=0.02;
    kGCdoHighPtQtGammaSelection=kFALSE; // RRnew
    kGCHighPtQtMax=100.;	        // RRnew
    kGCPtBorderForQt=100.;	        // RRnew
    break; // end RR ///////////////////////////////////////////////
  case 6:  // RRnew start: pT dependent qT cut
    kGCQtMax=0.02;
    kGCdoHighPtQtGammaSelection=kTRUE;
    kGCHighPtQtMax=0.06;
    kGCPtBorderForQt=2.5;
    break; // RRnew end ////////////////////////////////////////////
  case 7:
    kGCQtMax=0.15;
    kGCdoHighPtQtGammaSelection=kFALSE; // RRnew
    kGCHighPtQtMax=100.;	        // RRnew
    kGCPtBorderForQt=100.;	        // RRnew
    break;
  default:
    return iResult;
  }

  switch(piMaxMomdedxSigmaCut){
  case 0:  // 100. GeV
    kGCPIDMaxPnSigmaAbovePionLine=100.;
    break;
  case 1:  // 5. GeV
    kGCPIDMaxPnSigmaAbovePionLine=5.;
    break;
  case 2:  // 4. GeV
    kGCPIDMaxPnSigmaAbovePionLine=4.;
    break;
  case 3:  // 3.5 GeV
    kGCPIDMaxPnSigmaAbovePionLine=3.5;
    break;
  case 4:  // 3. GeV
    kGCPIDMaxPnSigmaAbovePionLine=3.;
    break;
  default:
    return iResult;
  }

	
  switch(alphaMesonCut){
  case 0:  // 0- 0.7
    kGCalphaMinCutMeson   = 0.0;
    kGCalphaCutMeson   = 0.7;
    break;
  case 1:  // 0-0.5
    kGCalphaMinCutMeson   = 0.0;
    kGCalphaCutMeson   = 0.5;
    break;
  case 2:  // 0.5-1 
    kGCalphaMinCutMeson   = 0.5;
    kGCalphaCutMeson   = 1.;
    break;
  case 3:  // 0.0-1 
    kGCalphaMinCutMeson   = 0.0;
    kGCalphaCutMeson   = 1.;
    break;
  default:
    return iResult;
  }
  
  switch(minRCut){
  case 0:
    kGCminRCut=0;
    kGCmaxRCut = 180.;
    break;
  case 1:
    kGCminRCut=2.8;
    kGCmaxRCut = 180.;
    break;
  case 2:
    kGCminRCut=5.;
    kGCmaxRCut = 180.;
    break;
  case 3:
    kGCmaxRCut = 70.;
    kGCminRCut = 10.;
    break;
  case 4:
    kGCmaxRCut = 70.;
    kGCminRCut = 5.;
    break;
  case 5:
    kGCmaxRCut = 180.;
    kGCminRCut = 10.;
    break;
    
  default:
    return iResult;
  }

  switch(RapidityMesonCut){
  case 0:  //
    kGCrapidityCutMeson   = 0.9;
    break;
  case 1:  // 
    kGCrapidityCutMeson   = 0.8;
    break;
  case 2:  // 
    kGCrapidityCutMeson   = 0.7;
    break;
    
  default:
    return iResult;
  }

  switch(BackgroundScheme){
  case 0: //Rotation
    kGCUseRotationMethodInBG=kTRUE;
    kGCdoBGProbability=kFALSE;
    break;
  case 1: // mixed event with track multiplicity
    kGCUseRotationMethodInBG=kFALSE;
    kGCUseTrackMultiplicityForBG=kFALSE;
    kGCdoBGProbability=kFALSE;
    break;
  case 2: // mixed event with V0 multiplicity
    kGCUseRotationMethodInBG=kFALSE;
    kGCUseTrackMultiplicityForBG=kTRUE;
    kGCdoBGProbability=kFALSE;
    break;
  case 3: //Rotation
    kGCUseRotationMethodInBG=kTRUE;
    kGCdoBGProbability=kTRUE;
    break;
  default:
    return iResult;
  }
  
  switch(DegreesForRotationMethod){
  case 0:
    kGCnDegreeRotationPMForBG = 5;
    break;
  case 1:
    kGCnDegreeRotationPMForBG = 10;
    break;
  case 2:
    kGCnDegreeRotationPMForBG = 15;
    break;
  case 3:
    kGCnDegreeRotationPMForBG = 20;
    break;

  default:
    return iResult;
  }

  switch(NumberOfRotations){
  case 0:
    kGCnumberOfRotationEventsForBG = 5;
    break;
  case 1:
    kGCnumberOfRotationEventsForBG = 10;
    break;
  case 2:
    kGCnumberOfRotationEventsForBG = 15;
    break;
  case 3:
    kGCnumberOfRotationEventsForBG = 20;
    break;
  case 4:
    kGCnumberOfRotationEventsForBG = 2;
    break;
  case 5:
    kGCnumberOfRotationEventsForBG = 50;
    break;
  case 6:
    kGCnumberOfRotationEventsForBG = 80;
    break;
  case 7:
    kGCnumberOfRotationEventsForBG = 100;
    break;
  default:
    return iResult;
  }

  switch(removePileUp){
  case 0:
    kGCRemovePileUp=kFALSE;
    break;
  case 1:
    kGCRemovePileUp=kTRUE;
    break;
  default:
    return iResult;
  }

  switch(selectV0AND){
  case 0:
    kGCSelectV0AND=kFALSE;
    break;
  case 1:
    kGCSelectV0AND=kTRUE;
    break;
  default:
    return iResult;
  }

  switch(multiplicityBin){
  case 0:
    kGCUseMultiplicity=kFALSE;
    kGCUseMultiplicityBin=0;
    break;
  case 1:
    kGCUseMultiplicity=kTRUE;
    kGCUseMultiplicityBin=1;
     break;
  case 2:
    kGCUseMultiplicity=kTRUE;
    kGCUseMultiplicityBin=2;
     break;
  case 3:
    kGCUseMultiplicity=kTRUE;
    kGCUseMultiplicityBin=3;
     break;
  case 4:
    kGCUseMultiplicity=kTRUE;
    kGCUseMultiplicityBin=4;
     break;
  case 5:
    kGCUseMultiplicity=kTRUE;
    kGCUseMultiplicityBin=5;
     break;
  default:
    return iResult;
  }

  switch(isHeavyIon){
  case 0:
    kGCIsHeavyIon=0;
    break;
  case 1:
    kGCIsHeavyIon=1;
    break;
   default:
    return iResult;
  }

  switch(useCentrality){
  case 0:
    kGCUseCentrality=0;
    break;
  case 1:
    kGCUseCentrality=1;//    GetCentralityClass10("V0M")
    break;
  case 2:
    kGCUseCentrality=2; //    GetCentralityClass10("CL1")
    break;
  case 3://////////////// RRnew start //////////////////////////////////////////////////////////////////////////////////////////////////////
    kGCUseCentrality=3;//    GetCentralityClass10("V0M") but with merged bins: see AliAnalysisTaskGammaConversion.cxx
    break;
  case 4:
    kGCUseCentrality=4;//    GetCentralityClass10("CL1") but with merged bins: see AliAnalysisTaskGammaConversion.cxx
    break;/////////////// RRnew end ////////////////////////////////////////////////////////////////////////////////////////////////////////
   default:
    return iResult;
  }

  switch(centralityBin){
  case 0:
    kGCUseCentralityBin=0;
    break;
  case 1:
    kGCUseCentralityBin=1;
    break;
  case 2:
    kGCUseCentralityBin=2;
    break;
  case 3:
    kGCUseCentralityBin=3;
    break;
  case 4:
    kGCUseCentralityBin=4;
    break;
  case 5:
    kGCUseCentralityBin=5;
    break;
  case 6:
    kGCUseCentralityBin=6;
    break;
  case 7:
    kGCUseCentralityBin=7;
    break;
  case 8:
    kGCUseCentralityBin=8;
    break;
  case 9:
    kGCUseCentralityBin=9;
    break;
   default:
    return iResult;
  }

  switch(TOFelectronPID){ // RRnewTOF start //////////////////////////////////////////////////////////////////////////
  case 0: // no cut
    kGCuseTOFpid = kFALSE;
    kGCtofPIDnSigmaBelowElectronLine=-100;
    kGCtofPIDnSigmaAboveElectronLine=100;
    break;
  case 1: // -7,7
    kGCuseTOFpid = kTRUE;
    kGCtofPIDnSigmaBelowElectronLine=-7;
    kGCtofPIDnSigmaAboveElectronLine=7;
    break;
  case 2: // -5,5 
    kGCuseTOFpid = kTRUE;
    kGCtofPIDnSigmaBelowElectronLine=-5;
    kGCtofPIDnSigmaAboveElectronLine=5;
    break;
  case 3: // -3,5
    kGCuseTOFpid = kTRUE;
    kGCtofPIDnSigmaBelowElectronLine=-3;
    kGCtofPIDnSigmaAboveElectronLine=5;
    break;
  case 4: // -2,3
    kGCuseTOFpid = kTRUE;
    kGCtofPIDnSigmaBelowElectronLine=-2;
    kGCtofPIDnSigmaAboveElectronLine=3;
    break;
  default:
    return iResult;
  } //////////////////////// RRnewTOF end //////////////////////////////////////////////////////////////////////////

  iResult = 1;
  return iResult;

}


void string2array(const std::string& number, int a[c_array_size]) 
{
    if (number.size() == c_array_size) {
#define ASSIGNARRAY(i)  a[i] = number[i] - '0'
        ASSIGNARRAY(0);
        ASSIGNARRAY(1);
        ASSIGNARRAY(2);
        ASSIGNARRAY(3);
        ASSIGNARRAY(4);
        ASSIGNARRAY(5);
        ASSIGNARRAY(6);
        ASSIGNARRAY(7);
        ASSIGNARRAY(8);
        ASSIGNARRAY(9);
        ASSIGNARRAY(10);
        ASSIGNARRAY(11);
        ASSIGNARRAY(12);
        ASSIGNARRAY(13);
        ASSIGNARRAY(14);
        ASSIGNARRAY(15);
        ASSIGNARRAY(16);
        ASSIGNARRAY(17);
        ASSIGNARRAY(18);
        ASSIGNARRAY(19);
        ASSIGNARRAY(20);
        ASSIGNARRAY(21);
        ASSIGNARRAY(22);
        ASSIGNARRAY(23);
        ASSIGNARRAY(24);
        ASSIGNARRAY(25);
        ASSIGNARRAY(26); // RRnewTOF
  }
}




