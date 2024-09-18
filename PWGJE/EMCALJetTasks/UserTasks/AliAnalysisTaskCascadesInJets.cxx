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
 
//-------------------------------------------------------------------------
// Task for Cascade analysis in charged jets with 
// the strange particles (instead of daughters) added to the  jet finder
// Author: Ekaterina Grecka (ermeeka@fjfi.cvut.cz)
//-------------------------------------------------------------------------
#include <vector>

#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THnSparse.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TVector3.h>
#include <TROOT.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliEventPoolManager.h"
#include "AliPIDResponse.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliEmcalParticle.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliEventCuts.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliParticleContainer.h"
#include "AliAODVertex.h"
#include "AliStack.h"
#include "AliEmcalJet.h"

#include "AliAnalysisTaskCascadesInJets.h"

ClassImp(AliAnalysisTaskCascadesInJets)

// upper edges of centrality bins
const Int_t AliAnalysisTaskCascadesInJets::fgkiCentBinRanges[] = {10}; // only central
// axis: pT of jets
const Double_t AliAnalysisTaskCascadesInJets::fgkdBinsPtJet[] = {0, 100};
const Int_t AliAnalysisTaskCascadesInJets::fgkiNBinsPtJet = sizeof(AliAnalysisTaskCascadesInJets::fgkdBinsPtJet) / sizeof(AliAnalysisTaskCascadesInJets::fgkdBinsPtJet[0]) - 1;
const Int_t AliAnalysisTaskCascadesInJets::fgkiNBinsPtJetInit = int(((AliAnalysisTaskCascadesInJets::fgkdBinsPtJet)[AliAnalysisTaskCascadesInJets::fgkiNBinsPtJet] - (AliAnalysisTaskCascadesInJets::fgkdBinsPtJet)[0]) / 10.); // bin width 10 GeV/c

// axis: pT of Cascade
const Double_t AliAnalysisTaskCascadesInJets::fgkdBinsPtCascade[] = {0, 18};
const Int_t AliAnalysisTaskCascadesInJets::fgkiNBinsPtCascade = sizeof(AliAnalysisTaskCascadesInJets::fgkdBinsPtCascade) / sizeof((AliAnalysisTaskCascadesInJets::fgkdBinsPtCascade)[0]) - 1;
const Int_t AliAnalysisTaskCascadesInJets::fgkiNBinsPtCascadeInit = int(((AliAnalysisTaskCascadesInJets::fgkdBinsPtCascade)[AliAnalysisTaskCascadesInJets::fgkiNBinsPtCascade] - (AliAnalysisTaskCascadesInJets::fgkdBinsPtCascade)[0]) / 0.1); // bin width 0.1 GeV/c
const Int_t AliAnalysisTaskCascadesInJets::fgkiNBinsPtCascadeInitInJet = int(((AliAnalysisTaskCascadesInJets::fgkdBinsPtCascade)[AliAnalysisTaskCascadesInJets::fgkiNBinsPtCascade] - (AliAnalysisTaskCascadesInJets::fgkdBinsPtCascade)[0]) / 0.5); // bin width 0.5 GeV/c

// axis: Xi invariant mass
const Int_t AliAnalysisTaskCascadesInJets::fgkiNBinsMassXi =200;
const Double_t AliAnalysisTaskCascadesInJets::fgkdMassXiMin = 1.3;// [GeV/c^2]
const Double_t AliAnalysisTaskCascadesInJets::fgkdMassXiMax = 1.35;// [GeV/c^2]
// axis: Omega invariant mass
const Int_t AliAnalysisTaskCascadesInJets::fgkiNBinsMassOmega =200;
const Double_t AliAnalysisTaskCascadesInJets::fgkdMassOmegaMin = 1.65; // [GeV/c^2]
const Double_t AliAnalysisTaskCascadesInJets::fgkdMassOmegaMax = 1.7; // [GeV/c^2]
// PDG codes of used particles
const Int_t AliAnalysisTaskCascadesInJets::iPdgCodePion = 211;
const Int_t AliAnalysisTaskCascadesInJets::iPdgCodeProton = 2212;
const Int_t AliAnalysisTaskCascadesInJets::iPdgCodeK0s = 310;
const Int_t AliAnalysisTaskCascadesInJets::iPdgCodeLambda = 3122;
const Int_t AliAnalysisTaskCascadesInJets::iPdgCodeKaon = 321;
const Int_t AliAnalysisTaskCascadesInJets::iPdgCodeXi = 3312;       //Int_t iPdgCodeXiPlus = -3312;
const Int_t AliAnalysisTaskCascadesInJets::iPdgCodeOmega = 3334;    //Int_t iPdgCodeOmegaPlus = -3334;
//Index codes to distinguish particles in the jet 
const Int_t AliAnalysisTaskCascadesInJets::iXiMinusId = 1000000;
const Int_t AliAnalysisTaskCascadesInJets::iXiPlusId = 2000000;
const Int_t AliAnalysisTaskCascadesInJets::iOmegaMinusId = 3000000;
const Int_t AliAnalysisTaskCascadesInJets::iOmegaPlusId = 4000000;

// Default constructor
AliAnalysisTaskCascadesInJets::AliAnalysisTaskCascadesInJets():
  AliAnalysisTaskEmcal(),
  fOutputListStd(0),
  fOutputListStdJets(0),
  fOutputListMC(0),
  fCascadeCandidateArray(0),
  fGenMCCascade(0),
  fNCand(0),
  
  fFastJetWrapper("AliEmcalJetTask","AliEmcalJetTask"),
  //fFastJetWrapperBG("AliEmcalJetTask","AliEmcalJetTask"),  
  fFastJetWrapperMCGen("AliEmcalJetTask","AliEmcalJetTask"), 

  fAODIn(0),
  fAODOut(0),
  fRandom(0),

  fbIsPbPb(0), 
  fbMCAnalysis(0),
  fsGeneratorName(""),   
  fbSignalInBG(0), 
  fdCutVertexZ(0),
  fdCutVertexR2(0),
  fdCutCentLow(0),
  fdCutCentHigh(0),
  fdCutDeltaZMax(0),
  fiNContribMin(0),
  fdCentrality(0),
  fbUseMultiplicity(0),
  fbUseIonutCut(0),
  fbTPCRefit(0),
  fbRejectKinks(0),
  fbFindableClusters(0),
  fdCutNCrossedRowsTPCMin(0),
  fdCutCrossedRowsOverFindMin(0),
  fdCutCrossedRowsOverFindMax(0),
  fdCutPtDaughterMin(0),
  fdCutDCAToPrimVtxMin(0),
  fdCutDCADaughtersMax(0),
  fdCutEtaDaughterMax(0),
  fdCutNSigmadEdxMax(0),
  fdPtProtonPIDMax(0),
  fbOnFly(0),
  fdCutCascadePtMin(1.),
  fdCutRadiusDecayMin(0),
  fdCutRadiusDecayMax(0),
  fdCutEtaCascadeMax(0),
  fdCutRapCascadeMax(0),

    
  fdGhostArea(0.005),
  fdRadius(0.4),
  fdMinJetArea(0.001),
  fdMinJetPt(5.0),
  fdJetPhiMin(-10),
  fdJetPhiMax(+10),
  fdJetEtaMin(-0.9),
  fdJetEtaMax(+0.9),
  fdJetTrackPtMin(0.15),
  fdJetTrackEtaMax(0.9),
  fdMaxEtaJetBG(0.7),
  fdBgRadius(0.4),
  fdCutPtJetMin(0),
  fdCutPtTrackJetMin(0),
  fdCutAreaPercJetMin(0.6),
  fdDistanceCascadeJetMax(0.4),
  fdDistPrimaryMax(0.01),
  
  fh1EventCounterCut(0),
  fh1EventCent(0),
  fh1EventCent2(0),
  fh2EventCentTracks(0),
  fh2EventCentMult(0),
  fh1NCascadesInJetStats(0),
  fh1NRndConeCent(0),
  fh1AreaExcluded(0),
  fh1MCStats(0)     
{	
  for(Int_t i = 0; i < fgkiNCategCascade; i++) {
    fh1CascadeInvMassXiMinusAll[i] = 0;
    fh1CascadeInvMassXiPlusAll[i] = 0;
    fh1CascadeInvMassOmegaMinusAll[i] = 0;
    fh1CascadeInvMassOmegaPlusAll[i] = 0;
  }	
  for(Int_t i = 0; i < fgkiNBinsCent; i++) {
    fh1EventCounterCutCent[i] = 0;
    fh1CascadeCounterCentXiMinus[i] = 0;
    fh1CascadeCandPerEventCentXiMinus[i] = 0;
    fh1CascadeInvMassXiMinusCent[i] = 0;
    fhnCascadeInclusiveXiMinus[i] = 0;
    fhnCascadeInvMassCutXiMinus[i] = 0;
    
    fh1CascadeCounterCentXiPlus[i] = 0;
    fh1CascadeCandPerEventCentXiPlus[i] = 0;
    fh1CascadeInvMassXiPlusCent[i] = 0;
    fhnCascadeInclusiveXiPlus[i] = 0;
    fhnCascadeInvMassCutXiPlus[i] = 0;

    fh1CascadeCounterCentOmegaMinus[i] = 0;
    fh1CascadeCandPerEventCentOmegaMinus[i] = 0;
    fh1CascadeInvMassOmegaMinusCent[i] = 0;
    fhnCascadeInclusiveOmegaMinus[i] = 0;
    fhnCascadeInvMassCutOmegaMinus[i] = 0;

    fh1CascadeCounterCentOmegaPlus[i] = 0;
    fh1CascadeCandPerEventCentOmegaPlus[i] = 0;
    fh1CascadeInvMassOmegaPlusCent[i] = 0;
    fhnCascadeInclusiveOmegaPlus[i] = 0;
    fhnCascadeInvMassCutOmegaPlus[i] = 0;

    fh1VtxZ[i] = 0; 
    fh2VtxXY[i] = 0;   
    
    fh1PtJet[i] = 0;
    fh1EtaJet[i] = 0;
    fh2EtaPtJet[i] = 0;
    fh1PhiJet[i] = 0;
    fh2PtJetPtTrackLeading[i] = 0;
    fh1NJetPerEvent[i] = 0;   
    fh1NCascadeJetPerEvent[i] = 0;
    fh2EtaPhiRndCone[i] = 0;

    //In jets 
    fhnCascadeInJetXiMinus[i] = 0;
    fhnCascadeInJetXiPlus[i] = 0;
    fhnCascadeInJetOmegaMinus[i] = 0;
    fhnCascadeInJetOmegaPlus[i] = 0;  

    //In cones 
    fhnCascadeInPerpXiMinus[i] = 0;
    fhnCascadeInRndXiMinus[i] = 0;
    fhnCascadeInMedXiMinus[i] = 0;
    fhnCascadeOutJetXiMinus[i] = 0;
    fhnCascadeNoJetXiMinus[i] = 0;

    fhnCascadeInPerpXiPlus[i] = 0;
    fhnCascadeInRndXiPlus[i] = 0;
    fhnCascadeInMedXiPlus[i] = 0;
    fhnCascadeOutJetXiPlus[i] = 0;
    fhnCascadeNoJetXiPlus[i] = 0;

    fhnCascadeInPerpOmegaMinus[i] = 0;
    fhnCascadeInRndOmegaMinus[i] = 0;
    fhnCascadeInMedOmegaMinus[i] = 0;
    fhnCascadeOutJetOmegaMinus[i] = 0;
    fhnCascadeNoJetOmegaMinus[i] = 0;
    
    fhnCascadeInPerpOmegaPlus[i] = 0;
    fhnCascadeInRndOmegaPlus[i] = 0;
    fhnCascadeInMedOmegaPlus[i] = 0;
    fhnCascadeOutJetOmegaPlus[i] = 0;
    fhnCascadeNoJetOmegaPlus[i] = 0;
    
    //MC 
    fh1CascadeXiMinusPtMCGen[i] = 0;
    fh2CascadeXiMinusPtMassMCRec[i] = 0;
    fh1CascadeXiMinusPtMCRecFalse[i] = 0;
    fh2CascadeXiMinusEtaPtMCGen[i] = 0;
    fh3CascadeXiMinusEtaPtMassMCRec[i] = 0;
    fh2CascadeXiMinusInJetPtMCGen[i] = 0;
    fh3CascadeXiMinusInJetPtMassMCRec[i] = 0;
    fh3CascadeXiMinusInJetEtaPtMCGen[i] = 0;
    fh4CascadeXiMinusInJetEtaPtMassMCRec[i] = 0;
    fh2CascadeXiMinusMCResolMPt[i] = 0;
    fh2CascadeXiMinusMCPtGenPtRec[i] = 0;
    
    fh1CascadeXiPlusPtMCGen[i] = 0;
    fh2CascadeXiPlusPtMassMCRec[i] = 0;
    fh1CascadeXiPlusPtMCRecFalse[i] = 0;
    fh2CascadeXiPlusEtaPtMCGen[i] = 0;
    fh3CascadeXiPlusEtaPtMassMCRec[i] = 0;
    fh2CascadeXiPlusInJetPtMCGen[i] = 0;
    fh3CascadeXiPlusInJetPtMassMCRec[i] = 0;
    fh3CascadeXiPlusInJetEtaPtMCGen[i] = 0;
    fh4CascadeXiPlusInJetEtaPtMassMCRec[i] = 0;
    fh2CascadeXiPlusMCResolMPt[i] = 0;
    fh2CascadeXiPlusMCPtGenPtRec[i] = 0;

    fh1CascadeOmegaMinusPtMCGen[i] = 0;
    fh2CascadeOmegaMinusPtMassMCRec[i] = 0;
    fh1CascadeOmegaMinusPtMCRecFalse[i] = 0;
    fh2CascadeOmegaMinusEtaPtMCGen[i] = 0;
    fh3CascadeOmegaMinusEtaPtMassMCRec[i] = 0;
    fh2CascadeOmegaMinusInJetPtMCGen[i] = 0;
    fh3CascadeOmegaMinusInJetPtMassMCRec[i] = 0;
    fh3CascadeOmegaMinusInJetEtaPtMCGen[i] = 0;
    fh4CascadeOmegaMinusInJetEtaPtMassMCRec[i] = 0;
    fh2CascadeOmegaMinusMCResolMPt[i] = 0;
    fh2CascadeOmegaMinusMCPtGenPtRec[i] = 0;

    fh1CascadeOmegaPlusPtMCGen[i] = 0;
    fh2CascadeOmegaPlusPtMassMCRec[i] = 0;
    fh1CascadeOmegaPlusPtMCRecFalse[i] = 0;
    fh2CascadeOmegaPlusEtaPtMCGen[i] = 0;
    fh3CascadeOmegaPlusEtaPtMassMCRec[i] = 0;
    fh2CascadeOmegaPlusInJetPtMCGen[i] = 0;
    fh3CascadeOmegaPlusInJetPtMassMCRec[i] = 0;
    fh3CascadeOmegaPlusInJetEtaPtMCGen[i] = 0;
    fh4CascadeOmegaPlusInJetEtaPtMassMCRec[i] = 0;
    fh2CascadeOmegaPlusMCResolMPt[i] = 0;
    fh2CascadeOmegaPlusMCPtGenPtRec[i] = 0;    

    // eta daughters
    fhnCascadeXiMinusInclDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeXiMinusInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeXiPlusInclDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeXiPlusInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeOmegaMinusInclDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeOmegaMinusInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeOmegaPlusInclDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeOmegaPlusInJetsDaughterEtaPtPtMCRec[i] = 0;
  }  
}


// Constructor
AliAnalysisTaskCascadesInJets::AliAnalysisTaskCascadesInJets(const char* name):
  AliAnalysisTaskEmcal(name, kTRUE),
  fOutputListStd(0),
  fOutputListStdJets(0),
  fOutputListMC(0),
  fCascadeCandidateArray(0),
  fGenMCCascade(0),                   //! contains MC generated Cascades

  fNCand(0),
  
  fFastJetWrapper(name,name),
  //fFastJetWrapperBG(name,name),
  fFastJetWrapperMCGen(name,name),
  
  fAODIn(0),
  fAODOut(0),
  fRandom(0),
      
  fbIsPbPb(0),  
  fbMCAnalysis(0),
  fsGeneratorName(""), 
  fbSignalInBG(0), 
  fdCutVertexZ(0),
  fdCutVertexR2(0),
  fdCutCentLow(0),
  fdCutCentHigh(0),
  fdCutDeltaZMax(0),
  fiNContribMin(0),
  fdCentrality(0),
  fbUseMultiplicity(0),
  fbUseIonutCut(0),
  fbTPCRefit(0),
  fbRejectKinks(0),
  fbFindableClusters(0),
  fdCutNCrossedRowsTPCMin(0),
  fdCutCrossedRowsOverFindMin(0),
  fdCutCrossedRowsOverFindMax(0),
  fdCutPtDaughterMin(0),
  fdCutDCAToPrimVtxMin(0),
  fdCutDCADaughtersMax(0),
  fdCutEtaDaughterMax(0),
  fdCutNSigmadEdxMax(0),
  fdPtProtonPIDMax(0),
  fbOnFly(0),
  fdCutCascadePtMin(1.),
  fdCutRadiusDecayMin(0),
  fdCutRadiusDecayMax(0),
  fdCutEtaCascadeMax(0),
  fdCutRapCascadeMax(0),

  fdGhostArea(0.005),
  fdRadius(0.4),
  fdMinJetArea(0.001),
  fdMinJetPt(5.0),
  fdJetPhiMin(-10),
  fdJetPhiMax(+10),
  fdJetEtaMin(-0.9),
  fdJetEtaMax(+0.9),
  fdJetTrackPtMin(0.15),
  fdJetTrackEtaMax(0.9),
  fdMaxEtaJetBG(0.7),
  fdBgRadius(0.4),
  fdCutPtJetMin(0),
  fdCutPtTrackJetMin(0),
  fdCutAreaPercJetMin(0.6),
  fdDistanceCascadeJetMax(0.4),
  fdDistPrimaryMax(0.01),
  
  fh1EventCounterCut(0),
  fh1EventCent(0),
  fh1EventCent2(0),
  fh2EventCentTracks(0),
  fh2EventCentMult(0),
  fh1NCascadesInJetStats(0),
  fh1NRndConeCent(0),
  fh1AreaExcluded(0),
  fh1MCStats(0)    
{
  for(Int_t i = 0; i < fgkiNCategCascade; i++) {
    fh1CascadeInvMassXiMinusAll[i] = 0;
    fh1CascadeInvMassXiPlusAll[i] = 0;
    fh1CascadeInvMassOmegaMinusAll[i] = 0;
    fh1CascadeInvMassOmegaPlusAll[i] = 0;
  }	
  for(Int_t i = 0; i < fgkiNBinsCent; i++) {
    fh1EventCounterCutCent[i] = 0;
    fh1CascadeCounterCentXiMinus[i] = 0;
    fh1CascadeCandPerEventCentXiMinus[i] = 0;
    fh1CascadeInvMassXiMinusCent[i] = 0;
    fhnCascadeInclusiveXiMinus[i] = 0;
    fhnCascadeInvMassCutXiMinus[i] = 0;
    
    fh1CascadeCounterCentXiPlus[i] = 0;
    fh1CascadeCandPerEventCentXiPlus[i] = 0;
    fh1CascadeInvMassXiPlusCent[i] = 0;
    fhnCascadeInclusiveXiPlus[i] = 0;
    fhnCascadeInvMassCutXiPlus[i] = 0;

    fh1CascadeCounterCentOmegaMinus[i] = 0;
    fh1CascadeCandPerEventCentOmegaMinus[i] = 0;
    fh1CascadeInvMassOmegaMinusCent[i] = 0;
    fhnCascadeInclusiveOmegaMinus[i] = 0;
    fhnCascadeInvMassCutOmegaMinus[i] = 0;

    fh1CascadeCounterCentOmegaPlus[i] = 0;
    fh1CascadeCandPerEventCentOmegaPlus[i] = 0;
    fh1CascadeInvMassOmegaPlusCent[i] = 0;
    fhnCascadeInclusiveOmegaPlus[i] = 0;
    fhnCascadeInvMassCutOmegaPlus[i] = 0;

    fh1VtxZ[i] = 0; 
    fh2VtxXY[i] = 0;   
    
    fh1PtJet[i] = 0;
    fh1EtaJet[i] = 0;
    fh2EtaPtJet[i] = 0;
    fh1PhiJet[i] = 0;
    fh2PtJetPtTrackLeading[i] = 0;
    fh1NJetPerEvent[i] = 0;   
    fh1NCascadeJetPerEvent[i] = 0;
    fh2EtaPhiRndCone[i] = 0;

    //In jets 
    fhnCascadeInJetXiMinus[i] = 0;
    fhnCascadeInJetXiPlus[i] = 0;
    fhnCascadeInJetOmegaMinus[i] = 0;
    fhnCascadeInJetOmegaPlus[i] = 0;  

    //In cones 
    fhnCascadeInPerpXiMinus[i] = 0;
    fhnCascadeInRndXiMinus[i] = 0;
    fhnCascadeInMedXiMinus[i] = 0;
    fhnCascadeOutJetXiMinus[i] = 0;
    fhnCascadeNoJetXiMinus[i] = 0;

    fhnCascadeInPerpXiPlus[i] = 0;
    fhnCascadeInRndXiPlus[i] = 0;
    fhnCascadeInMedXiPlus[i] = 0;
    fhnCascadeOutJetXiPlus[i] = 0;
    fhnCascadeNoJetXiPlus[i] = 0;

    fhnCascadeInPerpOmegaMinus[i] = 0;
    fhnCascadeInRndOmegaMinus[i] = 0;
    fhnCascadeInMedOmegaMinus[i] = 0;
    fhnCascadeOutJetOmegaMinus[i] = 0;
    fhnCascadeNoJetOmegaMinus[i] = 0;
    
    fhnCascadeInPerpOmegaPlus[i] = 0;
    fhnCascadeInRndOmegaPlus[i] = 0;
    fhnCascadeInMedOmegaPlus[i] = 0;
    fhnCascadeOutJetOmegaPlus[i] = 0;
    fhnCascadeNoJetOmegaPlus[i] = 0;
    
    //MC 
    fh1CascadeXiMinusPtMCGen[i] = 0;
    fh2CascadeXiMinusPtMassMCRec[i] = 0;
    fh1CascadeXiMinusPtMCRecFalse[i] = 0;
    fh2CascadeXiMinusEtaPtMCGen[i] = 0;
    fh3CascadeXiMinusEtaPtMassMCRec[i] = 0;
    fh2CascadeXiMinusInJetPtMCGen[i] = 0;
    fh3CascadeXiMinusInJetPtMassMCRec[i] = 0;
    fh3CascadeXiMinusInJetEtaPtMCGen[i] = 0;
    fh4CascadeXiMinusInJetEtaPtMassMCRec[i] = 0;
    fh2CascadeXiMinusMCResolMPt[i] = 0;
    fh2CascadeXiMinusMCPtGenPtRec[i] = 0;
    
    fh1CascadeXiPlusPtMCGen[i] = 0;
    fh2CascadeXiPlusPtMassMCRec[i] = 0;
    fh1CascadeXiPlusPtMCRecFalse[i] = 0;
    fh2CascadeXiPlusEtaPtMCGen[i] = 0;
    fh3CascadeXiPlusEtaPtMassMCRec[i] = 0;
    fh2CascadeXiPlusInJetPtMCGen[i] = 0;
    fh3CascadeXiPlusInJetPtMassMCRec[i] = 0;
    fh3CascadeXiPlusInJetEtaPtMCGen[i] = 0;
    fh4CascadeXiPlusInJetEtaPtMassMCRec[i] = 0;
    fh2CascadeXiPlusMCResolMPt[i] = 0;
    fh2CascadeXiPlusMCPtGenPtRec[i] = 0;

    fh1CascadeOmegaMinusPtMCGen[i] = 0;
    fh2CascadeOmegaMinusPtMassMCRec[i] = 0;
    fh1CascadeOmegaMinusPtMCRecFalse[i] = 0;
    fh2CascadeOmegaMinusEtaPtMCGen[i] = 0;
    fh3CascadeOmegaMinusEtaPtMassMCRec[i] = 0;
    fh2CascadeOmegaMinusInJetPtMCGen[i] = 0;
    fh3CascadeOmegaMinusInJetPtMassMCRec[i] = 0;
    fh3CascadeOmegaMinusInJetEtaPtMCGen[i] = 0;
    fh4CascadeOmegaMinusInJetEtaPtMassMCRec[i] = 0;
    fh2CascadeOmegaMinusMCResolMPt[i] = 0;
    fh2CascadeOmegaMinusMCPtGenPtRec[i] = 0;

    fh1CascadeOmegaPlusPtMCGen[i] = 0;
    fh2CascadeOmegaPlusPtMassMCRec[i] = 0;
    fh1CascadeOmegaPlusPtMCRecFalse[i] = 0;
    fh2CascadeOmegaPlusEtaPtMCGen[i] = 0;
    fh3CascadeOmegaPlusEtaPtMassMCRec[i] = 0;
    fh2CascadeOmegaPlusInJetPtMCGen[i] = 0;
    fh3CascadeOmegaPlusInJetPtMassMCRec[i] = 0;
    fh3CascadeOmegaPlusInJetEtaPtMCGen[i] = 0;
    fh4CascadeOmegaPlusInJetEtaPtMassMCRec[i] = 0;
    fh2CascadeOmegaPlusMCResolMPt[i] = 0;
    fh2CascadeOmegaPlusMCPtGenPtRec[i] = 0;    

    // eta daughters
    fhnCascadeXiMinusInclDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeXiMinusInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeXiPlusInclDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeXiPlusInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeOmegaMinusInclDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeOmegaMinusInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeOmegaPlusInclDaughterEtaPtPtMCRec[i] = 0;
    fhnCascadeOmegaPlusInJetsDaughterEtaPtPtMCRec[i] = 0;
  }  
	
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class()); //Histograms
  DefineOutput(2, TList::Class()); //Jet histograms
  DefineOutput(3, TList::Class()); //MC histograms
}

//Destructor
AliAnalysisTaskCascadesInJets::~AliAnalysisTaskCascadesInJets()
{
  delete fRandom;
  fRandom = 0;
  
  if (fCascadeCandidateArray) {
    delete fCascadeCandidateArray;
    fCascadeCandidateArray = 0;
  }
  if (fGenMCCascade) {
    delete fGenMCCascade;
    fGenMCCascade = 0;
  }
  if (fJets) {
    delete fJets;
    fJets = 0;
  }
}

void AliAnalysisTaskCascadesInJets::UserCreateOutputObjects()
{
  // Called once

  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  
  // Initialise random-number generator
  fRandom = new TRandom3(0);

  // Create histograms  
  fOutputListStd = new TList();
  fOutputListStd->SetOwner();
  fOutputListStdJets = new TList();
  fOutputListStdJets->SetOwner();
  fOutputListMC = new TList();
  fOutputListMC->SetOwner();
  //List of Cascade candidates 
  fCascadeCandidateArray = new TClonesArray("AliAODcascade", 10); 
  fCascadeCandidateArray->SetOwner();

  //List of Cascade candidates 
  fGenMCCascade = new TClonesArray("AliAODMCParticle", 10); 
  fGenMCCascade->SetOwner();

  
  //List of all jets  
  fJets = new TClonesArray("AliEmcalJet");
  fJets->SetOwner(); 	
  
  // event categories
  const Int_t iNCategEvent = 13;
  TString categEvent[iNCategEvent] = {
    "coll. candid.",      //0
    "AOD OK",             //1
    "pp SPD pile-up",     //2
    "PV contributors",    //3 
    "Ionut's cut",        //4
    "|z| PV",             //5
    "|#Deltaz| SPD PV",   //6 
    "r PV",               //7
    "centrality",         //8
    "with Cascade",       //9
    "with jets",          //10 
    "jet selection",      //11
    "jet with Cascade"    //12
  };

// labels for stages of Cascade selection
  TString categCascade[fgkiNCategCascade] = {
    "all",                  //0
    "min Cascade Pt",       //1
    "mass range",           //2
    "rec. method",          //3
    "tracks TPC",           //4
    "track pt",             //5     
    "DCA to PV",            //6
	  "DCA V0 daughters",     //7
	  "DCA bach to PV",       //8
	  "DCA V0 PV",            //9
	  "DCA bach to V0",       //10
	  "CPA V0",               //11
	  "CPA Casc",             //12
	  "volume V0",            //13
	  "volume Casc",          //14
	  "track #it{#eta}",      //15
	  "Cascade #it{y} & #it{#eta}", //16 
	  "lifetime",             //17
    "PID",                  //18
	  "inclusive",            //19
	  "in jet"                //20 
  };
  
  //labels for the in jet statistics
  const Int_t iNCategInJetStat = 11;
  TString categInJetStats[iNCategInJetStat] = {
    "N Jets",                    //0
    "Cascades in jets",          //1
    "XiMinus in jets",           //2
    "XiPlus in jets",            //3
    "OmegaMinus in jets",        //4
    "OmegaPlus in jets",         //5
    "More than 1 Cascade in jet", //6 
    "N jets before cuts(after bg sub)", //7
    "N jet after pt sel",           //8
    "N jet after area sel",         //9
    "N jet after lead track pt sel" //10
   }; 
 

//histograms 
  fh1EventCounterCut = new TH1D("fh1EventCounterCut", "Number of events after filtering;selection filter;counts", iNCategEvent, 0, iNCategEvent);
  for(Int_t i = 0; i < iNCategEvent; i++)
    fh1EventCounterCut->GetXaxis()->SetBinLabel(i + 1, categEvent[i].Data());
  fh1EventCent2 = new TH1D("fh1EventCent2", "Number of events vs centrality;centrality;counts", 100, 0, 100);
  fh2EventCentTracks = new TH2D("fh2EventCentTracks", "Number of tracks vs centrality;centrality;tracks;counts", 100, 0, 100, 150, 0, 15e3);
  fh2EventCentMult = new TH2D("fh2EventCentMult", "Ref. multiplicity vs centrality;centrality;multiplicity;counts", 100, 0, 100, 150, 0, 15e3);
  fh1EventCent = new TH1D("fh1EventCent", "Number of events in centrality bins;centrality;counts", fgkiNBinsCent, 0, fgkiNBinsCent);
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
    fh1EventCent->GetXaxis()->SetBinLabel(i + 1, GetCentBinLabel(i).Data());
  fh1CascadeCandPerEvent = new TH1D("fh1CascadeCandPerEvent", "Number of all Cascade candidates per event;candidates;events", 300, 0, 30000);
  fh1NRndConeCent = new TH1D("fh1NRndConeCent", "Number of rnd. cones in centrality bins;centrality;counts", fgkiNBinsCent, 0, fgkiNBinsCent);
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
    fh1NRndConeCent->GetXaxis()->SetBinLabel(i + 1, GetCentBinLabel(i).Data());
  fh1AreaExcluded = new TH1D("fh1AreaExcluded", "Area of excluded cones in centrality bins;centrality;area", fgkiNBinsCent, 0, fgkiNBinsCent);
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
  fh1AreaExcluded->GetXaxis()->SetBinLabel(i + 1, GetCentBinLabel(i).Data());
  
  fOutputListStd->Add(fh1EventCounterCut);
  fOutputListStd->Add(fh1EventCent);
  fOutputListStd->Add(fh1EventCent2);;
  fOutputListStd->Add(fh2EventCentTracks);
  fOutputListStd->Add(fh2EventCentMult);
  fOutputListStd->Add(fh1CascadeCandPerEvent);
  fOutputListStd->Add(fh1NRndConeCent);
  fOutputListStdJets->Add(fh1AreaExcluded);

  fh1NCascadesInJetStats = new TH1D("fh1NCascadesInJetStats", "Cascades in jets statistics", iNCategInJetStat, 0, iNCategInJetStat);   
  for(Int_t i = 0; i < iNCategInJetStat; i++)
  fh1NCascadesInJetStats->GetXaxis()->SetBinLabel(i + 1, categInJetStats[i].Data());
  fOutputListStdJets->Add(fh1NCascadesInJetStats);
  
  if(fbMCAnalysis) {
	//labels for theMC statistics
    const Int_t iNMCCategStat = 6;
    TString categMCStats[iNMCCategStat] = {
    "N rec Cascades",             //0
    "N rec jets",            //1
    "N rec Cascades in jets",     //2
    "N gen Cascades",             //3
    "N gen jets",            //4
    "N gen Cascades in jets",     //5
    };  
    fh1MCStats = new TH1D("fh1MCStats", "MC Cascades in jets statistics", iNMCCategStat, 0, iNMCCategStat);   
    for(Int_t i = 0; i < iNMCCategStat; i++)
      fh1MCStats->GetXaxis()->SetBinLabel(i + 1, categMCStats[i].Data());
    fOutputListMC->Add(fh1MCStats);
  }
  
  // pt binning for Cascade and jets
  Int_t iNBinsPtCascade = fgkiNBinsPtCascadeInit;
  Int_t iNBinsPtCascadeInJet = fgkiNBinsPtCascadeInitInJet;
  Double_t dPtCascadeMin = fgkdBinsPtCascade[0];
  Double_t dPtCascadeMax = fgkdBinsPtCascade[fgkiNBinsPtCascade];
  Int_t iNJetPtBins = fgkiNBinsPtJetInit;
  Double_t dJetPtMin = fgkdBinsPtJet[0];
  Double_t dJetPtMax = fgkdBinsPtJet[fgkiNBinsPtJet];
  
  Double_t dStepEtaCascade = 0.05;
  Double_t dRangeEtaCascadeMax = 0.8;
  const Int_t iNBinsEtaCascade = 2 * Int_t(dRangeEtaCascadeMax / dStepEtaCascade);
  //Binning parameters for invariant mass histograms
  const Int_t iNDimIncl = 3;
  Int_t binsXiIncl[iNDimIncl] = {fgkiNBinsMassXi, iNBinsPtCascade, iNBinsEtaCascade};
  Double_t xminXiIncl[iNDimIncl] = {fgkdMassXiMin, dPtCascadeMin, -dRangeEtaCascadeMax};
  Double_t xmaxXiIncl[iNDimIncl] = {fgkdMassXiMax, dPtCascadeMax, dRangeEtaCascadeMax};
  Int_t binsOmegaIncl[iNDimIncl] = {fgkiNBinsMassOmega, iNBinsPtCascade, iNBinsEtaCascade};
  Double_t xminOmegaIncl[iNDimIncl] = {fgkdMassOmegaMin, dPtCascadeMin, -dRangeEtaCascadeMax};
  Double_t xmaxOmegaIncl[iNDimIncl] = {fgkdMassOmegaMax, dPtCascadeMax, dRangeEtaCascadeMax};
  //Binning in jets
  const Int_t iNDimInJC = 4;
  Int_t binsXiInJC[iNDimInJC] = {fgkiNBinsMassXi, iNBinsPtCascadeInJet, iNBinsEtaCascade, iNJetPtBins};
  Double_t xminXiInJC[iNDimInJC] = {fgkdMassXiMin, dPtCascadeMin, -dRangeEtaCascadeMax, dJetPtMin};
  Double_t xmaxXiInJC[iNDimInJC] = {fgkdMassXiMax, dPtCascadeMax, dRangeEtaCascadeMax, dJetPtMax};
  Int_t binsOmegaInJC[iNDimInJC] = {fgkiNBinsMassOmega, iNBinsPtCascadeInJet, iNBinsEtaCascade, iNJetPtBins};
  Double_t xminOmegaInJC[iNDimInJC] = {fgkdMassOmegaMin, dPtCascadeMin, -dRangeEtaCascadeMax, dJetPtMin};
  Double_t xmaxOmegaInJC[iNDimInJC] = {fgkdMassOmegaMax, dPtCascadeMax, dRangeEtaCascadeMax, dJetPtMax};
 
  //Binning for MC
  // binning eff inclusive vs eta-pT
  Double_t dStepDeltaEta = 0.1;
  Double_t dRangeDeltaEtaMax = 0.5;
  const Int_t iNBinsDeltaEta = 2 * Int_t(dRangeDeltaEtaMax / dStepDeltaEta);
  Int_t binsEtaXi[3] = {fgkiNBinsMassXi, iNBinsPtCascade, iNBinsEtaCascade};
  Double_t xminEtaXi[3] = {fgkdMassXiMin, dPtCascadeMin, -dRangeEtaCascadeMax};
  Double_t xmaxEtaXi[3] = {fgkdMassXiMax, dPtCascadeMax, dRangeEtaCascadeMax};
  Int_t binsEtaOmega[3] = {fgkiNBinsMassOmega, iNBinsPtCascade, iNBinsEtaCascade};
  Double_t xminEtaOmega[3] = {fgkdMassOmegaMin, dPtCascadeMin, -dRangeEtaCascadeMax};
  Double_t xmaxEtaOmega[3] = {fgkdMassOmegaMax, dPtCascadeMax, dRangeEtaCascadeMax};
  // binning eff in jets vs eta-pT
  // associated
  Int_t binsEtaXiInRec[5] = {fgkiNBinsMassXi, iNBinsPtCascadeInJet, iNBinsEtaCascade, iNJetPtBins, iNBinsDeltaEta};
  Double_t xminEtaXiInRec[5] = {fgkdMassXiMin, dPtCascadeMin, -dRangeEtaCascadeMax, dJetPtMin, -dRangeDeltaEtaMax};
  Double_t xmaxEtaXiInRec[5] = {fgkdMassXiMax, dPtCascadeMax, dRangeEtaCascadeMax, dJetPtMax, dRangeDeltaEtaMax};
  Int_t binsEtaOmegaInRec[5] = {fgkiNBinsMassOmega, iNBinsPtCascadeInJet, iNBinsEtaCascade, iNJetPtBins, iNBinsDeltaEta};
  Double_t xminEtaOmegaInRec[5] = {fgkdMassOmegaMin, dPtCascadeMin, -dRangeEtaCascadeMax, dJetPtMin, -dRangeDeltaEtaMax};
  Double_t xmaxEtaOmegaInRec[5] = {fgkdMassOmegaMax, dPtCascadeMax, dRangeEtaCascadeMax, dJetPtMax, dRangeDeltaEtaMax};
  // generated
  Int_t binsEtaInGen[4] = {iNBinsPtCascadeInJet, iNBinsEtaCascade, iNJetPtBins, iNBinsDeltaEta};
  Double_t xminEtaInGen[4] = {dPtCascadeMin, -dRangeEtaCascadeMax, dJetPtMin, -dRangeDeltaEtaMax};
  Double_t xmaxEtaInGen[4] = {dPtCascadeMax, dRangeEtaCascadeMax, dJetPtMax, dRangeDeltaEtaMax};
  // daughter eta: charge-etaD-ptD-etaCascade-ptCascade-ptJet
  const Int_t iNDimEtaD = 6;
  Int_t binsEtaDaughter[iNDimEtaD] = {2, 20, iNBinsPtCascade, iNBinsEtaCascade, iNBinsPtCascade, iNJetPtBins};
  Double_t xminEtaDaughter[iNDimEtaD] = {0, -1, dPtCascadeMin, -dRangeEtaCascadeMax, dPtCascadeMin, dJetPtMin};
  Double_t xmaxEtaDaughter[iNDimEtaD] = {2, 1, dPtCascadeMax, dRangeEtaCascadeMax, dPtCascadeMax, dJetPtMax};
  // daughter track pt: pos-neg-ptCascade-ptJet-ptLead
  //const Int_t iNDimDaughter = 5;
  //Int_t binsDaughter[iNDimDaughter] = {80, 80, iNBinsPtCascade, iNJetPtBins, 80};
  //Double_t xminDaughter[iNDimDaughter] = {0, 0, dPtCascadeMin, dJetPtMin, 0};
  //Double_t xmaxDaughter[iNDimDaughter] = {20, 20, dPtCascadeMax, dJetPtMax, 20};
 
  for(Int_t i = 0; i < fgkiNBinsCent; i++) {
    fh1EventCounterCutCent[i] = new TH1D(Form("fh1EventCounterCutCent_%d", i), Form("Number of events after filtering, cent %s;selection filter;counts", GetCentBinLabel(i).Data()), iNCategEvent, 0, iNCategEvent);
    for(Int_t j = 0; j < iNCategEvent; j++)
      fh1EventCounterCutCent[i]->GetXaxis()->SetBinLabel(j + 1, categEvent[j].Data());
    fh1CascadeCandPerEventCentXiMinus[i] = new TH1D(Form("fh1CascadeCandPerEventCentXiMinus_%d", i), Form("Number of selected XiMinus candidates per event, cent %s;candidates;events", GetCentBinLabel(i).Data()), 100, 0, 200);
    fh1CascadeCandPerEventCentXiPlus[i] = new TH1D(Form("fh1CascadeCandPerEventCentXiPlus_%d", i), Form("Number of selected XiPlus candidates per event, cent %s;candidates;events", GetCentBinLabel(i).Data()), 100, 0, 200);
    fh1CascadeCandPerEventCentOmegaMinus[i] = new TH1D(Form("fh1CascadeCandPerEventCentOmegaMinus_%d", i), Form("Number of selected OmegaMinus candidates per event, cent %s;candidates;events", GetCentBinLabel(i).Data()), 100, 0, 200);
    fh1CascadeCandPerEventCentOmegaPlus[i] = new TH1D(Form("fh1CascadeCandPerEventCentOmegaPlus_%d", i), Form("Number of selected OmegaPlus candidates per event, cent %s;candidates;events", GetCentBinLabel(i).Data()), 100, 0, 200);
    fh1CascadeCounterCentXiMinus[i] = new TH1D(Form("fh1CascadeCounterCentXiMinus_%d", i), Form("Number of XiMinus candidates after cuts, cent %s;cut;counts", GetCentBinLabel(i).Data()), fgkiNCategCascade, 0, fgkiNCategCascade);
    fh1CascadeCounterCentXiPlus[i] = new TH1D(Form("fh1CascadeCounterCentXiPlus_%d", i), Form("Number of XiPlus candidates after cuts, cent %s;cut;counts", GetCentBinLabel(i).Data()), fgkiNCategCascade, 0, fgkiNCategCascade);
    fh1CascadeCounterCentOmegaMinus[i] = new TH1D(Form("fh1CascadeCounterCentOmegaMinus_%d", i), Form("Number of OmegaMinus candidates after cuts, cent %s;cut;counts", GetCentBinLabel(i).Data()), fgkiNCategCascade, 0, fgkiNCategCascade);
    fh1CascadeCounterCentOmegaPlus[i] = new TH1D(Form("fh1CascadeCounterCentOmegaPlus_%d", i), Form("Number of OmegaPlus candidates after cuts, cent %s;cut;counts", GetCentBinLabel(i).Data()), fgkiNCategCascade, 0, fgkiNCategCascade);

    for(Int_t j = 0; j < fgkiNCategCascade; j++) {	
      fh1CascadeCounterCentXiMinus[i]->GetXaxis()->SetBinLabel(j + 1, categCascade[j].Data());
      fh1CascadeCounterCentXiPlus[i]->GetXaxis()->SetBinLabel(j + 1, categCascade[j].Data());
      fh1CascadeCounterCentOmegaMinus[i]->GetXaxis()->SetBinLabel(j + 1, categCascade[j].Data());
      fh1CascadeCounterCentOmegaPlus[i]->GetXaxis()->SetBinLabel(j + 1, categCascade[j].Data());
    }
    fOutputListStd->Add(fh1EventCounterCutCent[i]);
    fOutputListStd->Add(fh1CascadeCandPerEventCentXiMinus[i]);
    fOutputListStd->Add(fh1CascadeCandPerEventCentXiPlus[i]);
    fOutputListStd->Add(fh1CascadeCandPerEventCentOmegaMinus[i]);
    fOutputListStd->Add(fh1CascadeCandPerEventCentOmegaPlus[i]);
    fOutputListStd->Add(fh1CascadeCounterCentXiMinus[i]);
    fOutputListStd->Add(fh1CascadeCounterCentXiPlus[i]);
    fOutputListStd->Add(fh1CascadeCounterCentOmegaMinus[i]);
    fOutputListStd->Add(fh1CascadeCounterCentOmegaPlus[i]);


    fh1CascadeInvMassXiMinusCent[i] = new TH1D(Form("fh1CascadeInvMassXiMinusCent_%d", i), Form("XiMinus: Cascade invariant mass, cent %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", GetCentBinLabel(i).Data()), fgkiNBinsMassXi, fgkdMassXiMin, fgkdMassXiMax);
    fh1CascadeInvMassXiPlusCent[i] = new TH1D(Form("fh1CascadeInvMassXiPlusCent_%d", i), Form("XiPlus: Cascade invariant mass, cent %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", GetCentBinLabel(i).Data()), fgkiNBinsMassXi, fgkdMassXiMin, fgkdMassXiMax);
    fh1CascadeInvMassOmegaMinusCent[i] = new TH1D(Form("fh1CascadeInvMassOmegaMinusCent_%d", i), Form("OmegaMinus: Cascade invariant mass, cent %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", GetCentBinLabel(i).Data()), fgkiNBinsMassOmega, fgkdMassOmegaMin, fgkdMassOmegaMax);
    fh1CascadeInvMassOmegaPlusCent[i] = new TH1D(Form("fh1CascadeInvMassOmegaPlusCent_%d", i), Form("OmegaPlus: Cascade invariant mass, cent %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", GetCentBinLabel(i).Data()), fgkiNBinsMassOmega, fgkdMassOmegaMin, fgkdMassOmegaMax);
    fOutputListStd->Add(fh1CascadeInvMassXiMinusCent[i]);
    fOutputListStd->Add(fh1CascadeInvMassXiPlusCent[i]);
    fOutputListStd->Add(fh1CascadeInvMassOmegaMinusCent[i]);
    fOutputListStd->Add(fh1CascadeInvMassOmegaPlusCent[i]);  
    // Inclusive
    fhnCascadeInclusiveXiMinus[i] = new THnSparseD(Form("fhnCascadeInclusiveXiMinus_C%d", i), "XiMinus: Cascade invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};counts", iNDimIncl, binsXiIncl, xminXiIncl, xmaxXiIncl);
    fhnCascadeInclusiveXiPlus[i] = new THnSparseD(Form("fhnCascadeInclusiveXiPlus_C%d", i), "XiPlus: Cascade invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};counts", iNDimIncl, binsXiIncl, xminXiIncl, xmaxXiIncl);
    fhnCascadeInclusiveOmegaMinus[i] = new THnSparseD(Form("fhnCascadeInclusiveOmegaMinus_C%d", i), "OmegaMinus: Cascade invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};counts", iNDimIncl, binsOmegaIncl, xminOmegaIncl, xmaxOmegaIncl);
    fhnCascadeInclusiveOmegaPlus[i] = new THnSparseD(Form("fhnCascadeInclusiveOmegaPlus_C%d", i), "OmegaPlus: Cascade invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};counts", iNDimIncl, binsOmegaIncl, xminOmegaIncl, xmaxOmegaIncl);
    fOutputListStd->Add(fhnCascadeInclusiveXiMinus[i]);
    fOutputListStd->Add(fhnCascadeInclusiveXiPlus[i]);
    fOutputListStd->Add(fhnCascadeInclusiveOmegaMinus[i]);
    fOutputListStd->Add(fhnCascadeInclusiveOmegaPlus[i]);    
    // After 3sigma invariant mass window cut
    fhnCascadeInvMassCutXiMinus[i] = new THnSparseD(Form("fhnCascadeInvMassCutXiMinus_C%d", i), "XiMinus after inv mass window cut: Cascade invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};counts", iNDimIncl, binsXiIncl, xminXiIncl, xmaxXiIncl);
    fhnCascadeInvMassCutXiPlus[i] = new THnSparseD(Form("fhnCascadeInvMassCutXiPlus_C%d", i), "XiPlus after inv mass window cut: Cascade invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};counts", iNDimIncl, binsXiIncl, xminXiIncl, xmaxXiIncl);
    fhnCascadeInvMassCutOmegaMinus[i] = new THnSparseD(Form("fhnCascadeInvMassCutOmegaMinus_C%d", i), "OmegaMinus after inv mass window cut: Cascade invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};counts", iNDimIncl, binsOmegaIncl, xminOmegaIncl, xmaxOmegaIncl);
    fhnCascadeInvMassCutOmegaPlus[i] = new THnSparseD(Form("fhnCascadeInvMassCutOmegaPlus_C%d", i), "OmegaPlus after inv mass window cut: Cascade invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};counts", iNDimIncl, binsOmegaIncl, xminOmegaIncl, xmaxOmegaIncl);
    fOutputListStd->Add(fhnCascadeInvMassCutXiMinus[i]);
    fOutputListStd->Add(fhnCascadeInvMassCutXiPlus[i]);
    fOutputListStd->Add(fhnCascadeInvMassCutOmegaMinus[i]);
    fOutputListStd->Add(fhnCascadeInvMassCutOmegaPlus[i]);

    // event histograms
    fh1VtxZ[i] = new TH1D(Form("fh1VtxZ_%d", i), Form("#it{z} coordinate of the primary vertex, cent: %s;#it{z} (cm)", GetCentBinLabel(i).Data()), 300, -15, 15);
    fOutputListStd->Add(fh1VtxZ[i]);
    fh2VtxXY[i] = new TH2D(Form("fh2VtxXY_%d", i), Form("#it{xy} coordinate of the primary vertex, cent: %s;#it{x} (cm);#it{y} (cm)", GetCentBinLabel(i).Data()), 200, -0.2, 0.2, 500, -0.5, 0.5);
    fOutputListStd->Add(fh2VtxXY[i]);
    
    //jet histograms
    fh1PtJet[i] = new TH1D(Form("fh1PtJet_%d", i), Form("Jet pt spectrum, cent: %s;#it{p}_{T} jet (GeV/#it{c})", GetCentBinLabel(i).Data()), 2 * iNJetPtBins, dJetPtMin, dJetPtMax);
    fh1EtaJet[i] = new TH1D(Form("fh1EtaJet_%d", i), Form("Jet eta spectrum, cent: %s;#it{#eta} jet", GetCentBinLabel(i).Data()), 80, -1., 1.);
    fh2EtaPtJet[i] = new TH2D(Form("fh2EtaPtJet_%d", i), Form("Jet eta vs pT spectrum, cent: %s;#it{#eta} jet;#it{p}_{T} jet (GeV/#it{c})", GetCentBinLabel(i).Data()), 80, -1., 1., 2 * iNJetPtBins, dJetPtMin, dJetPtMax);
    fh1PhiJet[i] = new TH1D(Form("fh1PhiJet_%d", i), Form("Jet phi spectrum, cent: %s;#it{#phi} jet", GetCentBinLabel(i).Data()), 90, 0., TMath::TwoPi());
    fh2PtJetPtTrackLeading[i] = new TH2D(Form("fh2PtJetPtTrackLeading_%d", i), Form("jet pt vs leading track pt, cent: %s;#it{p}_{T}^{jet} (GeV/#it{c});#it{p}_{T} leading track (GeV/#it{c})", GetCentBinLabel(i).Data()), 4 * iNJetPtBins, dJetPtMin, dJetPtMax, 200, 0., 20); 
    fh1NJetPerEvent[i] = new TH1D(Form("fh1NJetPerEvent_%d", i), Form("Number of selected jets per event, cent: %s;# jets;# events", GetCentBinLabel(i).Data()), 100, 0., 100.);
    fh1NCascadeJetPerEvent[i] = new TH1D(Form("fh1NCascadeJetPerEvent_%d", i), Form("Number of jets with Cascade per event, cent: %s;# jets;# events", GetCentBinLabel(i).Data()), 100, 0., 100.);
    fh2EtaPhiRndCone[i] = new TH2D(Form("fh2EtaPhiRndCone_%d", i), Form("Rnd. cones: eta vs phi, cent: %s;#it{#eta} cone;#it{#phi} cone", GetCentBinLabel(i).Data()), 80, -1., 1., 90, 0., TMath::TwoPi());

    // Cascades in jets
    fhnCascadeInJetXiMinus[i] = new THnSparseD(Form("fhnCascadeInJetXiMinus_%d", i), Form("XiMinus: Mass vs Pt in jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsXiInJC, xminXiInJC, xmaxXiInJC);
    fhnCascadeInJetXiPlus[i] = new THnSparseD(Form("fhnCascadeInJetXiPlus_%d", i), Form("XiPlus: Mass vs Pt in jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsXiInJC, xminXiInJC, xmaxXiInJC);
    fhnCascadeInJetOmegaMinus[i] = new THnSparseD(Form("fhnCascadeInJetOmegaMinus_%d", i), Form("OmegaMinus: Mass vs Pt in jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsOmegaInJC, xminOmegaInJC, xmaxOmegaInJC);
    fhnCascadeInJetOmegaPlus[i] = new THnSparseD(Form("fhnCascadeInJetOmegaPlus_%d", i), Form("OmegaPlus: Mass vs Pt in jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsOmegaInJC, xminOmegaInJC, xmaxOmegaInJC);
       
    fOutputListStdJets->Add(fh1PtJet[i]);
    fOutputListStdJets->Add(fh1EtaJet[i]);
    fOutputListStdJets->Add(fh2EtaPtJet[i]);
    fOutputListStdJets->Add(fh1PhiJet[i]);
    fOutputListStdJets->Add(fh2PtJetPtTrackLeading[i]);
    fOutputListStdJets->Add(fh1NJetPerEvent[i]);
    fOutputListStdJets->Add(fh1NCascadeJetPerEvent[i]);
    fOutputListStdJets->Add(fhnCascadeInJetXiMinus[i]);
    fOutputListStdJets->Add(fhnCascadeInJetXiPlus[i]);
    fOutputListStdJets->Add(fhnCascadeInJetOmegaMinus[i]);
    fOutputListStdJets->Add(fhnCascadeInJetOmegaPlus[i]);
    fOutputListStdJets->Add(fh2EtaPhiRndCone[i]);


    //Cascades in cones 
    fhnCascadeInPerpXiMinus[i] = new THnSparseD(Form("fhnCascadeInPerpXiMinus_%d", i), Form("XiMinus: Mass vs Pt in perp. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsXiInJC, xminXiInJC, xmaxXiInJC);
    fhnCascadeInRndXiMinus[i] = new THnSparseD(Form("fhnCascadeInRndXiMinus_%d", i), Form("XiMinus: Mass vs Pt in rnd. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsXiIncl, xminXiIncl, xmaxXiIncl);
    fhnCascadeInMedXiMinus[i] = new THnSparseD(Form("fhnCascadeInMedXiMinus_%d", i), Form("XiMinus: Mass vs Pt in med.-cl. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsXiIncl, xminXiIncl, xmaxXiIncl);
    fhnCascadeOutJetXiMinus[i] = new THnSparseD(Form("fhnCascadeOutJetXiMinus_%d", i), Form("XiMinus: Pt outside jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsXiIncl, xminXiIncl, xmaxXiIncl);
    fhnCascadeNoJetXiMinus[i] = new THnSparseD(Form("fhnCascadeNoJetXiMinus_%d", i), Form("XiMinus: Pt in jet-less events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsXiIncl, xminXiIncl, xmaxXiIncl);
    
    fOutputListStdJets->Add(fhnCascadeInPerpXiMinus[i]);
    fOutputListStdJets->Add(fhnCascadeInRndXiMinus[i]);
    fOutputListStdJets->Add(fhnCascadeInMedXiMinus[i]);
    fOutputListStdJets->Add(fhnCascadeOutJetXiMinus[i]);
    fOutputListStdJets->Add(fhnCascadeNoJetXiMinus[i]);

    fhnCascadeInPerpXiPlus[i] = new THnSparseD(Form("fhnCascadeInPerpXiPlus_%d", i), Form("XiPlus: Mass vs Pt in perp. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsXiInJC, xminXiInJC, xmaxXiInJC);
    fhnCascadeInRndXiPlus[i] = new THnSparseD(Form("fhnCascadeInRndXiPlus_%d", i), Form("XiPlus: Mass vs Pt in rnd. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsXiIncl, xminXiIncl, xmaxXiIncl);
    fhnCascadeInMedXiPlus[i] = new THnSparseD(Form("fhnCascadeInMedXiPlus_%d", i), Form("XiPlus: Mass vs Pt in med.-cl. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsXiIncl, xminXiIncl, xmaxXiIncl);
    fhnCascadeOutJetXiPlus[i] = new THnSparseD(Form("fhnCascadeOutJetXiPlus_%d", i), Form("XiPlus: Pt outside jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsXiIncl, xminXiIncl, xmaxXiIncl);
    fhnCascadeNoJetXiPlus[i] = new THnSparseD(Form("fhnCascadeNoJetXiPlus_%d", i), Form("XiPlus: Pt in jet-less events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsXiIncl, xminXiIncl, xmaxXiIncl);
    
    fOutputListStdJets->Add(fhnCascadeInPerpXiPlus[i]);
    fOutputListStdJets->Add(fhnCascadeInRndXiPlus[i]);
    fOutputListStdJets->Add(fhnCascadeInMedXiPlus[i]);
    fOutputListStdJets->Add(fhnCascadeOutJetXiPlus[i]);
    fOutputListStdJets->Add(fhnCascadeNoJetXiPlus[i]);

    fhnCascadeInPerpOmegaMinus[i] = new THnSparseD(Form("fhnCascadeInPerpOmegaMinus_%d", i), Form("OmegaMinus: Mass vs Pt in perp. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsOmegaInJC, xminOmegaInJC, xmaxOmegaInJC);
    fhnCascadeInRndOmegaMinus[i] = new THnSparseD(Form("fhnCascadeInRndOmegaMinus_%d", i), Form("OmegaMinus: Mass vs Pt in rnd. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsOmegaIncl, xminOmegaIncl, xmaxOmegaIncl);
    fhnCascadeInMedOmegaMinus[i] = new THnSparseD(Form("fhnCascadeInMedOmegaMinus_%d", i), Form("OmegaMinus: Mass vs Pt in med.-cl. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsOmegaIncl, xminOmegaIncl, xmaxOmegaIncl);
    fhnCascadeOutJetOmegaMinus[i] = new THnSparseD(Form("fhnCascadeOutJetOmegaMinus_%d", i), Form("OmegaMinus: Pt outside jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsOmegaIncl, xminOmegaIncl, xmaxOmegaIncl);
    fhnCascadeNoJetOmegaMinus[i] = new THnSparseD(Form("fhnCascadeNoJetOmegaMinus_%d", i), Form("OmegaMinus: Pt in jet-less events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsOmegaIncl, xminOmegaIncl, xmaxOmegaIncl);
    
    fOutputListStdJets->Add(fhnCascadeInPerpOmegaMinus[i]);
    fOutputListStdJets->Add(fhnCascadeInRndOmegaMinus[i]);
    fOutputListStdJets->Add(fhnCascadeInMedOmegaMinus[i]);
    fOutputListStdJets->Add(fhnCascadeOutJetOmegaMinus[i]);
    fOutputListStdJets->Add(fhnCascadeNoJetOmegaMinus[i]);

    fhnCascadeInPerpOmegaPlus[i] = new THnSparseD(Form("fhnCascadeInPerpOmegaPlus_%d", i), Form("OmegaPlus: Mass vs Pt in perp. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsOmegaInJC, xminOmegaInJC, xmaxOmegaInJC);
    fhnCascadeInRndOmegaPlus[i] = new THnSparseD(Form("fhnCascadeInRndOmegaPlus_%d", i), Form("OmegaPlus: Mass vs Pt in rnd. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsOmegaIncl, xminOmegaIncl, xmaxOmegaIncl);
    fhnCascadeInMedOmegaPlus[i] = new THnSparseD(Form("fhnCascadeInMedOmegaPlus_%d", i), Form("OmegaPlus: Mass vs Pt in med.-cl. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsOmegaIncl, xminOmegaIncl, xmaxOmegaIncl);
    fhnCascadeOutJetOmegaPlus[i] = new THnSparseD(Form("fhnCascadeOutJetOmegaPlus_%d", i), Form("OmegaPlus: Pt outside jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsOmegaIncl, xminOmegaIncl, xmaxOmegaIncl);
    fhnCascadeNoJetOmegaPlus[i] = new THnSparseD(Form("fhnCascadeNoJetOmegaPlus_%d", i), Form("OmegaPlus: Pt in jet-less events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{Cascade} (GeV/#it{c});#it{#eta}_{Cascade}", GetCentBinLabel(i).Data()), iNDimIncl, binsOmegaIncl, xminOmegaIncl, xmaxOmegaIncl);
    
    fOutputListStdJets->Add(fhnCascadeInPerpOmegaPlus[i]);
    fOutputListStdJets->Add(fhnCascadeInRndOmegaPlus[i]);
    fOutputListStdJets->Add(fhnCascadeInMedOmegaPlus[i]);
    fOutputListStdJets->Add(fhnCascadeOutJetOmegaPlus[i]);
    fOutputListStdJets->Add(fhnCascadeNoJetOmegaPlus[i]);

    if(fbMCAnalysis) {
      fh1CascadeXiMinusPtMCGen[i] = new TH1D(Form("fh1CascadeXiMinusPtMCGen_%d", i), Form("MC XiMinus generated: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fh2CascadeXiMinusPtMassMCRec[i] = new TH2D(Form("fh2CascadeXiMinusPtMassMCRec_%d", i), Form("MC XiMinus associated: pt-m spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{m}_{inv} (GeV/#it{c}^{2})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax, fgkiNBinsMassXi, fgkdMassXiMin, fgkdMassXiMax);
      fh1CascadeXiMinusPtMCRecFalse[i] = new TH1D(Form("fh1CascadeXiMinusPtMCRecFalse_%d", i), Form("MC XiMinus false: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh1CascadeXiMinusPtMCGen[i]);      
      fOutputListMC->Add(fh2CascadeXiMinusPtMassMCRec[i]);
      fOutputListMC->Add(fh1CascadeXiMinusPtMCRecFalse[i]);
      // inclusive pt-eta
      fh2CascadeXiMinusEtaPtMCGen[i] = new TH2D(Form("fh2CascadeXiMinusEtaPtMCGen_%d", i), Form("MC XiMinus generated: pt-eta spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax, iNBinsEtaCascade, -dRangeEtaCascadeMax, dRangeEtaCascadeMax);
      fh3CascadeXiMinusEtaPtMassMCRec[i] = new THnSparseD(Form("fh3CascadeXiMinusEtaPtMassMCRec_%d", i), Form("MC XiMinus associated: m-pt-eta spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), 3, binsEtaXi, xminEtaXi, xmaxEtaXi);
      fOutputListMC->Add(fh2CascadeXiMinusEtaPtMCGen[i]);
      fOutputListMC->Add(fh3CascadeXiMinusEtaPtMassMCRec[i]);
      // in jet pt
      fh2CascadeXiMinusInJetPtMCGen[i] = new TH2D(Form("fh2CascadeXiMinusInJetPtMCGen_%d", i), Form("MC XiMinus in jet generated: pt-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascadeInJet, dPtCascadeMin, dPtCascadeMax, iNJetPtBins, dJetPtMin, dJetPtMax);
      fh3CascadeXiMinusInJetPtMassMCRec[i] = new THnSparseD(Form("fh3CascadeXiMinusInJetPtMassMCRec_%d", i), Form("MC XiMinus in jet associated: m-pt-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsXiInJC, xminXiInJC, xmaxXiInJC); 
      fOutputListMC->Add(fh2CascadeXiMinusInJetPtMCGen[i]);      
      fOutputListMC->Add(fh3CascadeXiMinusInJetPtMassMCRec[i]);
      // in jet pt-eta
      fh3CascadeXiMinusInJetEtaPtMCGen[i] = new THnSparseD(Form("fh3CascadeXiMinusInJetEtaPtMCGen_%d", i), Form("MC XiMinus generated: pt-eta-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 4, binsEtaInGen, xminEtaInGen, xmaxEtaInGen);
      fh4CascadeXiMinusInJetEtaPtMassMCRec[i] = new THnSparseD(Form("fh4CascadeXiMinusInJetEtaPtMassMCRec_%d", i), Form("MC XiMinus associated: m-pt-eta-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 5, binsEtaXiInRec, xminEtaXiInRec, xmaxEtaXiInRec);
      fh2CascadeXiMinusMCResolMPt[i] = new TH2D(Form("fh2CascadeXiMinusMCResolMPt_%d", i), Form("MC XiMinus associated: #Delta#it{m} vs pt, cent %s;#Delta#it{m} = #it{m}_{inv} - #it{m}_{true} (GeV/#it{c}^{2});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), 100, -0.02, 0.02, iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fh2CascadeXiMinusMCPtGenPtRec[i] = new TH2D(Form("fh2CascadeXiMinusMCPtGenPtRec_%d", i), Form("MC XiMinus associated: pt gen vs pt rec, cent %s;#it{p}_{T}^{gen} (GeV/#it{c});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax, iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh3CascadeXiMinusInJetEtaPtMCGen[i]);
      fOutputListMC->Add(fh4CascadeXiMinusInJetEtaPtMassMCRec[i]);
      fOutputListMC->Add(fh2CascadeXiMinusMCResolMPt[i]);
      fOutputListMC->Add(fh2CascadeXiMinusMCPtGenPtRec[i]);

      // inclusive pt
      fh1CascadeXiPlusPtMCGen[i] = new TH1D(Form("fh1CascadeXiPlusPtMCGen_%d", i), Form("MC XiPlus generated: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh1CascadeXiPlusPtMCGen[i]);
      fh2CascadeXiPlusPtMassMCRec[i] = new TH2D(Form("fh2CascadeXiPlusPtMassMCRec_%d", i), Form("MC XiPlus associated: pt-m spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{m}_{inv} (GeV/#it{c}^{2})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax, fgkiNBinsMassXi, fgkdMassXiMin, fgkdMassXiMax);
      fOutputListMC->Add(fh2CascadeXiPlusPtMassMCRec[i]);
      fh1CascadeXiPlusPtMCRecFalse[i] = new TH1D(Form("fh1CascadeXiPlusPtMCRecFalse_%d", i), Form("MC XiPlus false: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh1CascadeXiPlusPtMCRecFalse[i]);
      // inclusive pt-eta
      fh2CascadeXiPlusEtaPtMCGen[i] = new TH2D(Form("fh2CascadeXiPlusEtaPtMCGen_%d", i), Form("MC XiPlus generated: pt-eta spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax, iNBinsEtaCascade, -dRangeEtaCascadeMax, dRangeEtaCascadeMax);
      fOutputListMC->Add(fh2CascadeXiPlusEtaPtMCGen[i]);
      fh3CascadeXiPlusEtaPtMassMCRec[i] = new THnSparseD(Form("fh3CascadeXiPlusEtaPtMassMCRec_%d", i), Form("MC XiPlus associated: m-pt-eta spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), 3, binsEtaXi, xminEtaXi, xmaxEtaXi);
      fOutputListMC->Add(fh3CascadeXiPlusEtaPtMassMCRec[i]);
      // in jet pt
      fh2CascadeXiPlusInJetPtMCGen[i] = new TH2D(Form("fh2CascadeXiPlusInJetPtMCGen_%d", i), Form("MC XiPlus in jet generated: pt-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascadeInJet, dPtCascadeMin, dPtCascadeMax, iNJetPtBins, dJetPtMin, dJetPtMax);
      fOutputListMC->Add(fh2CascadeXiPlusInJetPtMCGen[i]);
      fh3CascadeXiPlusInJetPtMassMCRec[i] = new THnSparseD(Form("fh3CascadeXiPlusInJetPtMassMCRec_%d", i), Form("MC XiPlus in jet associated: m-pt-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsXiInJC, xminXiInJC, xmaxXiInJC);
      fOutputListMC->Add(fh3CascadeXiPlusInJetPtMassMCRec[i]);
      // in jet pt-eta
      fh3CascadeXiPlusInJetEtaPtMCGen[i] = new THnSparseD(Form("fh3CascadeXiPlusInJetEtaPtMCGen_%d", i), Form("MC XiPlus generated: pt-eta-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 4, binsEtaInGen, xminEtaInGen, xmaxEtaInGen);
      fOutputListMC->Add(fh3CascadeXiPlusInJetEtaPtMCGen[i]);
      fh4CascadeXiPlusInJetEtaPtMassMCRec[i] = new THnSparseD(Form("fh4CascadeXiPlusInJetEtaPtMassMCRec_%d", i), Form("MC XiPlus associated: m-pt-eta-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 5, binsEtaXiInRec, xminEtaXiInRec, xmaxEtaXiInRec);
      fOutputListMC->Add(fh4CascadeXiPlusInJetEtaPtMassMCRec[i]);

      fh2CascadeXiPlusMCResolMPt[i] = new TH2D(Form("fh2CascadeXiPlusMCResolMPt_%d", i), Form("MC XiPlus associated: #Delta#it{m} vs pt, cent %s;#Delta#it{m} = #it{m}_{inv} - #it{m}_{true} (GeV/#it{c}^{2});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), 100, -0.02, 0.02, iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh2CascadeXiPlusMCResolMPt[i]);
      fh2CascadeXiPlusMCPtGenPtRec[i] = new TH2D(Form("fh2CascadeXiPlusMCPtGenPtRec_%d", i), Form("MC XiPlus associated: pt gen vs pt rec, cent %s;#it{p}_{T}^{gen} (GeV/#it{c});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax, iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh2CascadeXiPlusMCPtGenPtRec[i]);

      // inclusive pt
      fh1CascadeOmegaMinusPtMCGen[i] = new TH1D(Form("fh1CascadeOmegaMinusPtMCGen_%d", i), Form("MC OmegaMinus generated: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh1CascadeOmegaMinusPtMCGen[i]);
      fh2CascadeOmegaMinusPtMassMCRec[i] = new TH2D(Form("fh2CascadeOmegaMinusPtMassMCRec_%d", i), Form("MC OmegaMinus associated: pt-m spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{m}_{inv} (GeV/#it{c}^{2})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax, fgkiNBinsMassOmega, fgkdMassOmegaMin, fgkdMassOmegaMax);
      fOutputListMC->Add(fh2CascadeOmegaMinusPtMassMCRec[i]);
      fh1CascadeOmegaMinusPtMCRecFalse[i] = new TH1D(Form("fh1CascadeOmegaMinusPtMCRecFalse_%d", i), Form("MC OmegaMinus false: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh1CascadeOmegaMinusPtMCRecFalse[i]);
      // inclusive pt-eta
      fh2CascadeOmegaMinusEtaPtMCGen[i] = new TH2D(Form("fh2CascadeOmegaMinusEtaPtMCGen_%d", i), Form("MC OmegaMinus generated: pt-eta spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax, iNBinsEtaCascade, -dRangeEtaCascadeMax, dRangeEtaCascadeMax);
      fOutputListMC->Add(fh2CascadeOmegaMinusEtaPtMCGen[i]);
      fh3CascadeOmegaMinusEtaPtMassMCRec[i] = new THnSparseD(Form("fh3CascadeOmegaMinusEtaPtMassMCRec_%d", i), Form("MC OmegaMinus associated: m-pt-eta spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), 3, binsEtaOmega, xminEtaOmega, xmaxEtaOmega);
      fOutputListMC->Add(fh3CascadeOmegaMinusEtaPtMassMCRec[i]);
      // in jet pt
      fh2CascadeOmegaMinusInJetPtMCGen[i] = new TH2D(Form("fh2CascadeOmegaMinusInJetPtMCGen_%d", i), Form("MC OmegaMinus in jet generated: pt-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascadeInJet, dPtCascadeMin, dPtCascadeMax, iNJetPtBins, dJetPtMin, dJetPtMax);
      fOutputListMC->Add(fh2CascadeOmegaMinusInJetPtMCGen[i]);
      fh3CascadeOmegaMinusInJetPtMassMCRec[i] = new THnSparseD(Form("fh3CascadeOmegaMinusInJetPtMassMCRec_%d", i), Form("MC OmegaMinus in jet associated: m-pt-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsOmegaInJC, xminOmegaInJC, xmaxOmegaInJC);
      fOutputListMC->Add(fh3CascadeOmegaMinusInJetPtMassMCRec[i]);
      // in jet pt-eta
      fh3CascadeOmegaMinusInJetEtaPtMCGen[i] = new THnSparseD(Form("fh3CascadeOmegaMinusInJetEtaPtMCGen_%d", i), Form("MC OmegaMinus generated: pt-eta-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 4, binsEtaInGen, xminEtaInGen, xmaxEtaInGen);
      fOutputListMC->Add(fh3CascadeOmegaMinusInJetEtaPtMCGen[i]);
      fh4CascadeOmegaMinusInJetEtaPtMassMCRec[i] = new THnSparseD(Form("fh4CascadeOmegaMinusInJetEtaPtMassMCRec_%d", i), Form("MC OmegaMinus associated: m-pt-eta-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 5, binsEtaOmegaInRec, xminEtaOmegaInRec, xmaxEtaOmegaInRec);
      fOutputListMC->Add(fh4CascadeOmegaMinusInJetEtaPtMassMCRec[i]);

      fh2CascadeOmegaMinusMCResolMPt[i] = new TH2D(Form("fh2CascadeOmegaMinusMCResolMPt_%d", i), Form("MC OmegaMinus associated: #Delta#it{m} vs pt, cent %s;#Delta#it{m} = #it{m}_{inv} - #it{m}_{true} (GeV/#it{c}^{2});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), 100, -0.02, 0.02, iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh2CascadeOmegaMinusMCResolMPt[i]);
      fh2CascadeOmegaMinusMCPtGenPtRec[i] = new TH2D(Form("fh2CascadeOmegaMinusMCPtGenPtRec_%d", i), Form("MC OmegaMinus associated: pt gen vs pt rec, cent %s;#it{p}_{T}^{gen} (GeV/#it{c});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax, iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh2CascadeOmegaMinusMCPtGenPtRec[i]);
      
      // inclusive pt
      fh1CascadeOmegaPlusPtMCGen[i] = new TH1D(Form("fh1CascadeOmegaPlusPtMCGen_%d", i), Form("MC OmegaPlus generated: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh1CascadeOmegaPlusPtMCGen[i]);
      fh2CascadeOmegaPlusPtMassMCRec[i] = new TH2D(Form("fh2CascadeOmegaPlusPtMassMCRec_%d", i), Form("MC OmegaPlus associated: pt-m spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{m}_{inv} (GeV/#it{c}^{2})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax, fgkiNBinsMassOmega, fgkdMassOmegaMin, fgkdMassOmegaMax);
      fOutputListMC->Add(fh2CascadeOmegaPlusPtMassMCRec[i]);
      fh1CascadeOmegaPlusPtMCRecFalse[i] = new TH1D(Form("fh1CascadeOmegaPlusPtMCRecFalse_%d", i), Form("MC OmegaPlus false: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh1CascadeOmegaPlusPtMCRecFalse[i]);
      // inclusive pt-eta
      fh2CascadeOmegaPlusEtaPtMCGen[i] = new TH2D(Form("fh2CascadeOmegaPlusEtaPtMCGen_%d", i), Form("MC OmegaPlus generated: pt-eta spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax, iNBinsEtaCascade, -dRangeEtaCascadeMax, dRangeEtaCascadeMax);
      fOutputListMC->Add(fh2CascadeOmegaPlusEtaPtMCGen[i]);
      fh3CascadeOmegaPlusEtaPtMassMCRec[i] = new THnSparseD(Form("fh3CascadeOmegaPlusEtaPtMassMCRec_%d", i), Form("MC OmegaPlus associated: m-pt-eta spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), 3, binsEtaOmega, xminEtaOmega, xmaxEtaOmega);
      fOutputListMC->Add(fh3CascadeOmegaPlusEtaPtMassMCRec[i]);
      // in jet pt
      fh2CascadeOmegaPlusInJetPtMCGen[i] = new TH2D(Form("fh2CascadeOmegaPlusInJetPtMCGen_%d", i), Form("MC OmegaPlus in jet generated: pt-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascadeInJet, dPtCascadeMin, dPtCascadeMax, iNJetPtBins, dJetPtMin, dJetPtMax);
      fOutputListMC->Add(fh2CascadeOmegaPlusInJetPtMCGen[i]);
      fh3CascadeOmegaPlusInJetPtMassMCRec[i] = new THnSparseD(Form("fh3CascadeOmegaPlusInJetPtMassMCRec_%d", i), Form("MC OmegaPlus in jet associated: m-pt-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsOmegaInJC, xminOmegaInJC, xmaxOmegaInJC);
      fOutputListMC->Add(fh3CascadeOmegaPlusInJetPtMassMCRec[i]);
      // in jet pt-eta
      fh3CascadeOmegaPlusInJetEtaPtMCGen[i] = new THnSparseD(Form("fh3CascadeOmegaPlusInJetEtaPtMCGen_%d", i), Form("MC OmegaPlus generated: pt-eta-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 4, binsEtaInGen, xminEtaInGen, xmaxEtaInGen);
      fOutputListMC->Add(fh3CascadeOmegaPlusInJetEtaPtMCGen[i]);
      fh4CascadeOmegaPlusInJetEtaPtMassMCRec[i] = new THnSparseD(Form("fh4CascadeOmegaPlusInJetEtaPtMassMCRec_%d", i), Form("MC OmegaPlus associated: m-pt-eta-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 5, binsEtaOmegaInRec, xminEtaOmegaInRec, xmaxEtaOmegaInRec);
      fOutputListMC->Add(fh4CascadeOmegaPlusInJetEtaPtMassMCRec[i]);

      fh2CascadeOmegaPlusMCResolMPt[i] = new TH2D(Form("fh2CascadeOmegaPlusMCResolMPt_%d", i), Form("MC OmegaPlus associated: #Delta#it{m} vs pt, cent %s;#Delta#it{m} = #it{m}_{inv} - #it{m}_{true} (GeV/#it{c}^{2});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), 100, -0.02, 0.02, iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh2CascadeOmegaPlusMCResolMPt[i]);
      fh2CascadeOmegaPlusMCPtGenPtRec[i] = new TH2D(Form("fh2CascadeOmegaPlusMCPtGenPtRec_%d", i), Form("MC OmegaPlus associated: pt gen vs pt rec, cent %s;#it{p}_{T}^{gen} (GeV/#it{c});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax, iNBinsPtCascade, dPtCascadeMin, dPtCascadeMax);
      fOutputListMC->Add(fh2CascadeOmegaPlusMCPtGenPtRec[i]);

      // daughter eta
      fhnCascadeXiMinusInclDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnCascadeXiMinusInclDaughterEtaPtPtMCRec_%d", i), Form("MC K0S, inclusive, assoc., daughters: charge-etaD-ptD-etaCascade-ptCascade-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen Cascade;pT gen Cascade;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnCascadeXiMinusInJetsDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnCascadeXiMinusInJetsDaughterEtaPtPtMCRec_%d", i), Form("MC K0S, in JC, assoc., daughters: charge-etaD-ptD-etaCascade-ptCascade-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen Cascade;pT gen Cascade;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnCascadeXiPlusInclDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnCascadeXiPlusInclDaughterEtaPtPtMCRec_%d", i), Form("MC XiPlus, inclusive, assoc., daughters: charge-etaD-ptD-etaCascade-ptCascade-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen Cascade;pT gen Cascade;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnCascadeXiPlusInJetsDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnCascadeXiPlusInJetsDaughterEtaPtPtMCRec_%d", i), Form("MC XiPlus, in JC, assoc., daughters: charge-etaD-ptD-etaCascade-ptCascade-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen Cascade;pT gen Cascade;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnCascadeOmegaMinusInclDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnCascadeOmegaMinusInclDaughterEtaPtPtMCRec_%d", i), Form("MC OmegaMinus, inclusive, assoc., daughters: charge-etaD-ptD-etaCascade-ptCascade-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen Cascade;pT gen Cascade;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnCascadeOmegaMinusInJetsDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnCascadeOmegaMinusInJetsDaughterEtaPtPtMCRec_%d", i), Form("MC OmegaMinus, in JC, assoc., daughters: charge-etaD-ptD-etaCascade-ptCascade-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen Cascade;pT gen Cascade;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnCascadeOmegaPlusInclDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnCascadeOmegaPlusInclDaughterEtaPtPtMCRec_%d", i), Form("MC OmegaPlus, inclusive, assoc., daughters: charge-etaD-ptD-etaCascade-ptCascade-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen Cascade;pT gen Cascade;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnCascadeOmegaPlusInJetsDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnCascadeOmegaPlusInJetsDaughterEtaPtPtMCRec_%d", i), Form("MC OmegaPlus, in JC, assoc., daughters: charge-etaD-ptD-etaCascade-ptCascade-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen Cascade;pT gen Cascade;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);

      fOutputListMC->Add(fhnCascadeXiMinusInclDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnCascadeXiMinusInJetsDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnCascadeXiPlusInclDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnCascadeXiPlusInJetsDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnCascadeOmegaMinusInclDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnCascadeOmegaMinusInJetsDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnCascadeOmegaPlusInclDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnCascadeOmegaPlusInJetsDaughterEtaPtPtMCRec[i]);   
    }
    
  }
  for(Int_t i = 0; i < fgkiNCategCascade; i++) {
    fh1CascadeInvMassXiMinusAll[i] = new TH1D(Form("fh1CascadeInvMassXiMinusAll_%d", i), Form("XiMinus: Cascade invariant mass, %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", categCascade[i].Data()), fgkiNBinsMassXi, fgkdMassXiMin, fgkdMassXiMax);
    fh1CascadeInvMassXiPlusAll[i] = new TH1D(Form("fh1CascadeInvMassXiPlusAll_%d", i), Form("XiPlus: Cascade invariant mass, %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", categCascade[i].Data()), fgkiNBinsMassXi, fgkdMassXiMin, fgkdMassXiMax);
    fh1CascadeInvMassOmegaMinusAll[i] = new TH1D(Form("fh1CascadeInvMassOmegaMinusAll_%d", i), Form("OmegaMinus: Cascade invariant mass, %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", categCascade[i].Data()), fgkiNBinsMassOmega, fgkdMassOmegaMin, fgkdMassOmegaMax);
    fh1CascadeInvMassOmegaPlusAll[i] = new TH1D(Form("fh1CascadeInvMassOmegaPlusAll_%d", i), Form("OmegaPlus: Cascade invariant mass, %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", categCascade[i].Data()), fgkiNBinsMassOmega, fgkdMassOmegaMin, fgkdMassOmegaMax);
    fOutputListStd->Add(fh1CascadeInvMassXiMinusAll[i]);
    fOutputListStd->Add(fh1CascadeInvMassXiPlusAll[i]);
    fOutputListStd->Add(fh1CascadeInvMassOmegaMinusAll[i]);
    fOutputListStd->Add(fh1CascadeInvMassOmegaPlusAll[i]);
  }



  for(Int_t i = 0; i < fOutputListStd->GetEntries(); ++i) {
    TH1* h1 = dynamic_cast<TH1*>(fOutputListStd->At(i));
    if(h1) {
      h1->Sumw2();
      continue;
    }
    THnSparse* hn = dynamic_cast<THnSparse*>(fOutputListStd->At(i));
    if(hn) hn->Sumw2();
  } 
  for(Int_t i = 0; i < fOutputListStdJets->GetEntries(); ++i) {
    TH1* h1 = dynamic_cast<TH1*>(fOutputListStdJets->At(i));
    if(h1) {
      h1->Sumw2();
      continue;
    }
    THnSparse* hn = dynamic_cast<THnSparse*>(fOutputListStdJets->At(i));
    if(hn) hn->Sumw2();
  }  
  for(Int_t i = 0; i < fOutputListMC->GetEntries(); ++i) {
    TH1* h1 = dynamic_cast<TH1*>(fOutputListMC->At(i));
    if(h1) {
      h1->Sumw2();
      continue;
    }
    THnSparse* hn = dynamic_cast<THnSparse*>(fOutputListMC->At(i));
    if(hn) hn->Sumw2();
  }  


  PostData(1, fOutputListStd);
  PostData(2, fOutputListStdJets);
  PostData(3, fOutputListMC);
}

void AliAnalysisTaskCascadesInJets::ExecOnce()
{
  AliAnalysisTaskEmcal::ExecOnce();
  
  //Open AOD here? 

  // setup fj wrapper
  fFastJetWrapper.SetAreaType(fastjet::active_area_explicit_ghosts);
  fFastJetWrapper.SetGhostArea(fdGhostArea);
  fFastJetWrapper.SetR(fdRadius);
  fFastJetWrapper.SetAlgorithm(fastjet::antikt_algorithm); //rewrite with the right algoritm
  fFastJetWrapper.SetRecombScheme(fastjet::pt_scheme); //rewrite it
  fFastJetWrapper.SetMaxRap(1);
  
  // setup bg fj wrapper
  /*
  fFastJetWrapperBG.SetAreaType(fastjet::active_area_explicit_ghosts);
  fFastJetWrapperBG.SetGhostArea(fdGhostArea);
  fFastJetWrapperBG.SetR(fdRadius);
  fFastJetWrapperBG.SetAlgorithm(fastjet::kt_algorithm); //rewrite with the right algoritm
  fFastJetWrapperBG.SetRecombScheme(fastjet::pt_scheme); //rewrite it
  fFastJetWrapperBG.SetMaxRap(1);
  */
  
  // setup MC fj wrapper
  fFastJetWrapperMCGen.SetAreaType(fastjet::active_area_explicit_ghosts);
  fFastJetWrapperMCGen.SetGhostArea(fdGhostArea);
  fFastJetWrapperMCGen.SetR(fdRadius);
  fFastJetWrapperMCGen.SetAlgorithm(fastjet::antikt_algorithm); //rewrite with the right algoritm
  fFastJetWrapperMCGen.SetRecombScheme(fastjet::pt_scheme); //rewrite it
  fFastJetWrapperMCGen.SetMaxRap(1);
    
  selectorBG = !fastjet::SelectorNHardest(2) * fastjet::SelectorAbsEtaMax(fdMaxEtaJetBG) * fastjet::SelectorPtRange(fdJetTrackPtMin, 1000.0); //set the max eta cut on the estimator, then get rid of 2 highest pt jets
 
  printf("=======================================================\n");
  printf("%s: Configuration summary:\n", ClassName());
  printf("task name: %s\n", GetName());
  printf("-------------------------------------------------------\nEvent selection:\n");
  printf("collision system: %s\n", fbIsPbPb ? "Pb+Pb" : "p+p");
  printf("data type: %s\n", fbMCAnalysis ? "MC" : "real");
  if(fbMCAnalysis)
    printf("MC generator: %s\n", fsGeneratorName.Length() ? fsGeneratorName.Data() : "any");
  if(fbIsPbPb) {
    printf("centrality range: %g-%g %%\n", fdCutCentLow, fdCutCentHigh);
    printf("centrality estimator: %s\n", fbUseMultiplicity ? "AliMultSelection" : "AliCentrality");
  }
  if(fdCutCascadePtMin > 0.) printf("min Cascade pT [Gev/c]: %g\n", fdCutCascadePtMin);
  if(fdCutVertexZ > 0.) printf("max |z| of the prim vtx [cm]: %g\n", fdCutVertexZ);
  if(fdCutVertexR2 > 0.) printf("max r^2 of the prim vtx [cm^2]: %g\n", fdCutVertexR2);
  if(fiNContribMin > 0) printf("min number of prim vtx contributors: %d\n", fiNContribMin);
  if(fdCutDeltaZMax > 0.) printf("max |Delta z| between nominal prim vtx and SPD vtx [cm]: %g\n", fdCutDeltaZMax);
  if(fbUseIonutCut) printf("Ionut's cut\n");
  printf("-------------------------------------------------------\nCascade selection:\n");
  if(fbTPCRefit) printf("TPC refit for daughter tracks\n");
  if(fbRejectKinks) printf("reject kink-like production vertices of daughter tracks\n");
  if(fbFindableClusters) printf("require positive number of findable clusters\n");
  if(fdCutNCrossedRowsTPCMin > 0.) printf("min number of crossed TPC rows: %g\n", fdCutNCrossedRowsTPCMin);
  if(fdCutCrossedRowsOverFindMin > 0.) printf("min ratio crossed rows / findable clusters: %g\n", fdCutCrossedRowsOverFindMin);
  if(fdCutCrossedRowsOverFindMax > 0.) printf("max ratio crossed rows / findable clusters: %g\n", fdCutCrossedRowsOverFindMax);
  if(fdCutPtDaughterMin > 0.) printf("min pt of daughter tracks [GeV/c]: %g\n", fdCutPtDaughterMin);
  if(fdCutDCAToPrimVtxMin > 0.) printf("min DCA of daughters to the prim vtx [cm]: %g\n", fdCutDCAToPrimVtxMin);
  if(fdCutDCADaughtersMax > 0.) printf("max DCA between daughters [sigma of TPC tracking]: %g\n", fdCutDCADaughtersMax);
  if(fdCutEtaDaughterMax > 0.) printf("max |eta| of daughter tracks: %g\n", fdCutEtaDaughterMax);
  if(fdCutNSigmadEdxMax > 0. && (!fbIsPbPb || (fbIsPbPb && fdPtProtonPIDMax > 0.))) printf("max |Delta(dE/dx)| in the TPC [sigma dE/dx]: %g\n", fdCutNSigmadEdxMax);
  if(fdCutNSigmadEdxMax > 0. && fbIsPbPb && fdPtProtonPIDMax > 0.) printf("max pt of proton for applying PID cut [GeV/c]: %g\n", fdPtProtonPIDMax);
  printf("Cascade reconstruction method: %s\n", fbOnFly ? "on-the-fly" : "offline");
  if(fdCutCPACascadeMin > 0.) printf("min Cascade CPA: %g\n", fdCutCPACascadeMin);
  if(fdCutCPACascadeV0Min > 0.) printf("min daughter v0 CPA: %g\n", fdCutCPACascadeV0Min);
  if(fdCutRadiusDecayMin > 0. && fdCutRadiusDecayMax > 0.) printf("R of the decay vertex [cm]: %g-%g\n", fdCutRadiusDecayMin, fdCutRadiusDecayMax);
  if(fdCutEtaCascadeMax > 0.) printf("max |eta| of Cascade: %g\n", fdCutEtaCascadeMax);
  if(fdCutRapCascadeMax > 0.) printf("max |y| of Cascade: %g\n", fdCutRapCascadeMax);
  
  if(fdCutNTauXMax > 0.) printf("max proper lifetime, Xi [tau]: %g\n", fdCutNTauXMax);

  printf("Jet cuts:\n");   
  printf("Jet radius: %g\n", fdRadius); 
  printf("min jet area: %g\n", fdMinJetArea); 
  printf("min jet pt [GeV/c]: %g\n", fdMinJetPt); 
  printf("min jet phi: %g, max jet phi: %g\n", fdJetPhiMin, fdJetPhiMax); 
  printf("min jet eta: %g, max jet eta: %g\n", fdJetEtaMin, fdJetEtaMax); 
  printf("minimum jet constituent track pt [GeV/c]: %g\n", fdJetTrackPtMin);
  
  printf("min jet pt [GeV/c]: %g\n", fdCutPtJetMin);

  printf("=======================================================\n");
  
  //add configuration for jet analysis
    
}  

Bool_t AliAnalysisTaskCascadesInJets::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  //Move the analysis part from FillHistograms here? 

  return kTRUE; // If return kFALSE FillHistogram() will NOT be executed.
}

Bool_t AliAnalysisTaskCascadesInJets::FillHistograms()
{
  // Main loop, called for each event
  
  //open AOD 
  if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Start");

  //clear the TClonesArray from the previous event
  fCascadeCandidateArray->Clear();
  fGenMCCascade->Clear();
  fJets->Clear();
  fFastJetWrapper.Clear();
  //fFastJetWrapperBG.Clear();
  fFastJetWrapperMCGen.Clear();  
  
  std::vector<fastjet::PseudoJet> InputBgParticles;
  //InputBgParticles.clear();  //clear the particle vector for the bg estimator
  
  fh1EventCounterCut->Fill(0); // all available selected events (collision candidates)

  fAODIn = dynamic_cast<AliAODEvent*>(InputEvent()); // input AOD
  if(!fAODIn) {
    AliWarning("Input AOD not available from InputEvent() trying with AODEvent().");
    // Assume that the AOD is in the general output.
    fAODIn = AODEvent();
    if(!fAODIn) {
      AliError("No input AOD found!");
      return kFALSE;
    }
  }
  if(fDebug > 1) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Loading AOD OK");
  
  // PID Response Task object
  AliPIDResponse* fPIDResponse = 0;
  if(fdCutNSigmadEdxMax > 0. && (!fbIsPbPb || (fbIsPbPb && fdPtProtonPIDMax > 0.))) {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)mgr->GetInputEventHandler();
    fPIDResponse = inputHandler->GetPIDResponse();
    if(!fPIDResponse) {
      AliError("No PID response object found!");
      return kFALSE;
    }
  }

  // AOD files are OK
  fh1EventCounterCut->Fill(1);
	  
  //===== Event selection =====
  if(!IsSelectedForAnalysis()) {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Event rejected");
    return kFALSE;
  }
  if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Event accepted: cent. %g", fdCentrality));
  if(!fbIsPbPb)
    fdCentrality = 0.; // select the first bin for p+p data
  Int_t iCentIndex = GetCentralityBinIndex(fdCentrality); // get index of centrality bin
  if(iCentIndex < 0) {
    AliError(Form("Event is out of histogram range: cent: %g!", fdCentrality));
    return kFALSE;
  }
  fh1EventCounterCut->Fill(8); // selected events (centrality OK)
  fh1EventCounterCutCent[iCentIndex]->Fill(8);	
	
  //Get tracks from AOD	
  UInt_t iNTracks = fAODIn->GetNumberOfTracks(); // get number of tracks in event

  //if(fDebug > 2)
  printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("There are %d tracks in this event", iNTracks));
  TClonesArray* fTracks = fAODIn->GetTracks();
  if(!fTracks) 
	AliError("Could not get tracks from AOD.");
  AliAODHeader* header = dynamic_cast<AliAODHeader*>(fAODIn->GetHeader());
  UInt_t iNTracksRef = 0;
  if(!header)
    AliError("Cannot get AOD header");
  else
    iNTracksRef = header->GetRefMultiplicityComb08(); // get reference multiplicity

  //Get Cascade information from AOD 
  Int_t iNCascades = fAODIn->GetNumberOfCascades(); // get the number of Cascade candidates
  printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("There are %d Cascade candidates in the event", iNCascades));
  if(!iNCascades) {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "No Cascades found in event");
  } else {
    fh1EventCounterCut->Fill(9); // events with Cascades
    fh1EventCounterCutCent[iCentIndex]->Fill(9);
  }
  
  if(fbMCAnalysis) {
    fEventMC = MCEvent();
    if(!fEventMC) {
      AliError("No MC event found!");
      return kFALSE;
    } 
  }
  //===== Event is OK for the analysis =====
  //Fill event histograms
  fh1EventCent->Fill(iCentIndex);
  fh1EventCent2->Fill(fdCentrality);
  fh2EventCentTracks->Fill(fdCentrality, iNTracks);
  fh2EventCentMult->Fill(fdCentrality, iNTracksRef);
	
  //Prepare cuts for the Cascade analysis	
  AliAODcascade* Cascade = 0; // pointer to Cascade candidates
  TVector3 vecCascadeMomentum; // 3D vector of Cascade momentum

  Double_t dMassCascadeXi = 0; // invariant mass of the XiMinus candidate
  Double_t dMassCascadeOmega = 0; // invariant mass of the OmegaMinus candidate
  Int_t iNCascadeCandTot = 0; // counter of all Cascade candidates at the beginning
  Int_t iNCascadeCandXiMinus = 0; // counter of XiMinus candidates at the end
  Int_t iNCascadeCandXiPlus = 0; // counter of XiPlus candidates at the end
  Int_t iNCascadeCandOmegaMinus = 0; // counter of OmegaMinus candidates at the end
  Int_t iNCascadeCandOmegaPlus = 0; // counter of OmegaPlus candidates at the end

  // Values of Cascade reconstruction cuts:   
  // Daughter tracks
  Double_t dDCACascadeBachToPrimVtxMin = fdCutDCACascadeBachToPrimVtxMin; // 0.04; // [cm] min DCA of bachelor track to the prim vtx
  Double_t dDCACascadeV0ToPrimVtxMin = fdCutDCACascadeV0ToPrimVtxMin; // 0.06; // [cm] min DCA of V0 (from cascade decay) to the prim vtx
  Double_t dDCACascadeDaughtersToPrimVtxMin = fdCutDCACascadeDaughtersToPrimVtxMin; // 0.03; // [cm] min DCA of daughters (from secondary V0 decay) to the prim vtx  
  Double_t dDCACascadeV0DaughtersMax = fdCutDCACascadeV0DaughtersMax; // 1.5; // [sigma of TPC tracking] max DCA between Cascade V0 daughters  
  Double_t dDCACascadeBachToV0Max = fdCutDCACascadeBachToV0Max; // 1.3; // [cm] max DCA between bachelor track to V0     
  Double_t dCascadeEtaDaughterMax = 0.8; // max |pseudorapidity| of daughter tracks in cascade 
  // Cascade candidate
  Double_t dCPACascadeMin = fdCutCPACascadeMin;// 0.97; // min cosine of the pointing angle of the cascade
  Double_t dCPACascadeV0Min = fdCutCPACascadeV0Min;// 0.97; // min cosine of the pointing angle of the V0 in the cascade   
  Double_t dCascadeRadiusDecayMin =0.6; // [cm] min radial distance of the decay vertex of the cascade
  Double_t dCascadeV0RadiusDecayMin = 1.2; // [cm] min radial distance of the decay vertex of the V0 (from the cascade)  
  Double_t dCascadeRapMax = 0.75;
  Double_t dCascadeEtaMax = 0.75; // max |pseudorapidity| of Cascade
  // Selection of active cuts
  Bool_t bCutEtaV0Daughter = 1; // V0 daughter pseudorapidity
  Bool_t bCutEtaCascade = 1; // Cascade pseudorapidity

  Double_t dCTauXi = 4.917; // [cm] c tau of Xi
  Double_t dCTauOmega = 2.463 ; // [cm] c tau of Omega
  
  Bool_t bPrintCuts = 0; // print out which cuts are applied
  // Other cuts
  Double_t dNSigmaMassMax = 3.; // [sigma m] max difference between candidate mass and real particle mass (used only for mass peak method of signal extraction)
  // particle masses from PDG
  Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
  Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
  Double_t dMassPDGXiMinus = TDatabasePDG::Instance()->GetParticle(kXiMinus)->Mass();       //Double_t dMassPDGXiPlus = TDatabasePDG::Instance()->GetParticle(kXiPlus)->Mass();  
  Double_t dMassPDGOmegaMinus = TDatabasePDG::Instance()->GetParticle(kOmegaMinus)->Mass(); //Double_t dMassPDGXOmegaPlus = TDatabasePDG::Instance()->GetParticle(kOmegaPlus)->Mass();  

  fNCand = 0;
  Int_t iXiMinusId = 100000;
  Int_t iXiPlusId = 200000;
  Int_t iOmegaMinusId = 300000;
  Int_t iOmegaPlusId = 400000;

  std::vector<Int_t> ivecCascadeCandIndex;
  // Loading primary vertex info
  AliAODVertex* primVtx = fAODIn->GetPrimaryVertex(); // get the primary vertex
  Double_t dPrimVtxPos[3]; // primary vertex position {x,y,z}
  primVtx->GetXYZ(dPrimVtxPos);
  fh1VtxZ[iCentIndex]->Fill(dPrimVtxPos[2]);
  fh2VtxXY[iCentIndex]->Fill(dPrimVtxPos[0], dPrimVtxPos[1]);
  
  Double_t dCutEtaJetMax = fdCutEtaCascadeMax - fdDistanceCascadeJetMax; // max jet |pseudorapidity|, to make sure that V0s can appear in the entire jet area
  Double_t dRadiusExcludeCone = 2 * fdDistanceCascadeJetMax; // radius of cones around jets excluded for V0 outside jets

  Double_t dAreaPercJetMin =  fdCutAreaPercJetMin*TMath::Pi()*fdRadius*fdRadius;
//Loop over  Cascade candidates
//------------------------------------------------------------------------------------------

  //===== Start of loop over Cascade candidates =====
  if(fDebug > 0) printf(": Start of Cascade loop\n");
  for(Int_t iXi = 0; iXi < iNCascades; iXi++)
  {
    Cascade = fAODIn->GetCascade(iXi);   // get next candidate from the list in AOD 
	  if(!Cascade)
	    continue;
 
    iNCascadeCandTot++;

    Int_t iCutIndex = 0; // indicator of current selection step  
    // Initialization of status indicators
    Bool_t bIsCandidateXiMinus = kTRUE; // candidate for XiMinus
    Bool_t bIsCandidateXiPlus = kTRUE; // candidate for XiPlus
    Bool_t bIsCandidateOmegaMinus = kTRUE; // candidate for OmegaMinus
    Bool_t bIsCandidateOmegaPlus = kTRUE; // candidate for OmegaPlus
    Bool_t bIsInPeakXi = kFALSE; // candidate within the Xi mass peak
    Bool_t bIsInPeakOmega = kFALSE; // candidate within the Omega mass peak
         
    // Invariant mass calculation
    dMassCascadeXi = Cascade->MassXi();
    dMassCascadeOmega = Cascade->MassOmega();
    
    // 0
    // All Cascade candidates
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++;
 
    Double_t dPtCascade = TMath::Sqrt(Cascade->Pt2Xi()); // transverse momentum of Cascade
    vecCascadeMomentum = TVector3(Cascade->MomXiX(), Cascade->MomXiY(), Cascade->MomXiZ()); // set the vector of Cascade momentum
    
    //1
    //Cascade min pt 
    if(dPtCascade < fdCutCascadePtMin)
      continue;
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++;

    // Sigma of the mass peak window
    Double_t dMassPeakWindowXi = dNSigmaMassMax * MassPeakSigma(dPtCascade, 0); 
    Double_t dMassPeakWindowOmega = dNSigmaMassMax * MassPeakSigma(dPtCascade, 1);
    //Mean of the mass peak window
    Double_t dMassPeakWindowMeanXi = MassPeakMean(dPtCascade, 0);
    Double_t dMassPeakWindowMeanOmega = MassPeakMean(dPtCascade, 1);
    // Invariant mass peak selection
    if(fbSignalInBG == 0) { //Inv mass in signal region
      if(TMath::Abs(dMassCascadeXi - dMassPeakWindowMeanXi) < dMassPeakWindowXi)
        bIsInPeakXi = kTRUE;
      if(TMath::Abs(dMassCascadeOmega - dMassPeakWindowMeanOmega) < dMassPeakWindowOmega)
        bIsInPeakOmega = kTRUE;
    }  
    else if(fbSignalInBG == 1){ //Inv mass in BG region
      if( (dMassCascadeXi > dMassPeakWindowMeanXi + dMassPeakWindowXi) && (dMassCascadeXi < dMassPeakWindowMeanXi + 2 * dMassPeakWindowXi))
        bIsInPeakXi = kTRUE;
      if( (dMassCascadeOmega > dMassPeakWindowMeanOmega + dMassPeakWindowOmega) &&  (dMassCascadeOmega < dMassPeakWindowMeanOmega + 2 * dMassPeakWindowOmega))
        bIsInPeakOmega = kTRUE;
    }

    // Skip candidates outside the histogram range
    if((dMassCascadeXi < fgkdMassXiMin) || (dMassCascadeXi >= fgkdMassXiMax)) {
      bIsCandidateXiMinus = kFALSE; 
      bIsCandidateXiPlus  = kFALSE;
    }
    if((dMassCascadeOmega < fgkdMassOmegaMin) || (dMassCascadeOmega >= fgkdMassOmegaMax)) {
      bIsCandidateOmegaMinus = kFALSE; 
      bIsCandidateOmegaPlus  = kFALSE;
    }
    if(!bIsCandidateXiMinus && !bIsCandidateXiPlus && !bIsCandidateOmegaMinus && !bIsCandidateOmegaPlus )
      continue;
    //Charge check
    if(Cascade->ChargeXi() > 0) {
		bIsCandidateXiMinus = kFALSE;
		bIsCandidateOmegaMinus = kFALSE;
    }
    if(Cascade->ChargeXi() < 0) {
		bIsCandidateXiPlus = kFALSE;
		bIsCandidateOmegaPlus = kFALSE;
    } 
        
    // Retrieving all relevant properties of the Cascade candidate
    Bool_t bOnFlyStatus = Cascade->GetOnFlyStatus(); // online (on fly) reconstructed vs offline reconstructed 
    //Cascade vertex position 
    Double_t dXiVtxPos[3]; // Xi vertex position {x,y,z}
    dXiVtxPos[0] = Cascade->DecayVertexXiX();
    dXiVtxPos[1] = Cascade->DecayVertexXiY();
    dXiVtxPos[2] = Cascade->DecayVertexXiZ();    
    Double_t dRadiusDecayXi = TMath::Sqrt(dXiVtxPos[0] * dXiVtxPos[0] + dXiVtxPos[1] * dXiVtxPos[1]); // distance of the Cascade vertex from the z-axis    
    //V0 decay vertex position 
    Double_t dV0VtxPos[3]; // V0 vertex position {x,y,z}
    Cascade->GetSecondaryVtx(dV0VtxPos);
    Double_t dRadiusDecayV0 = TMath::Sqrt(dV0VtxPos[0] * dV0VtxPos[0] + dV0VtxPos[1] * dV0VtxPos[1]); // distance of the V0 vertex from the z-axis
	  // Tracks of the daughters variables 
    const AliAODTrack* trackBach = (AliAODTrack*)Cascade->GetDecayVertexXi()->GetDaughter(0); // bachelor daughter track    
    const AliAODTrack* trackPos  = (AliAODTrack*)Cascade->GetDaughter(0); // positive daughter track
    const AliAODTrack* trackNeg  = (AliAODTrack*)Cascade->GetDaughter(1); // negative daughter track  
    //Double_t dPtBach = trackBach->Pt(); 
    Double_t dPtDaughterPos = trackPos->Pt(); // transverse momentum of a daughter track
    Double_t dPtDaughterNeg = trackNeg->Pt();
    Double_t dEtaDaughterNeg = trackNeg->Eta(); // = Cascade->EtaProng(1), pseudorapidity of a daughter track
    Double_t dEtaDaughterPos = trackPos->Eta(); // = Cascade->EtaProng(0)      
    Double_t dNRowsPos = trackPos->GetTPCClusterInfo(2, 1); // crossed TPC pad rows of a daughter track
    Double_t dNRowsNeg = trackNeg->GetTPCClusterInfo(2, 1);
	  //DCA and CPA variables 
    Double_t dDCAToPrimVtxPos  = TMath::Abs(Cascade->DcaPosToPrimVertex()); // dca of a positive daughter to the primary vertex
    Double_t dDCAToPrimVtxNeg  = TMath::Abs(Cascade->DcaNegToPrimVertex()); // dca of a negative daughter to the primary vertex
    Double_t dDCAV0Daughters   = Cascade->DcaV0Daughters();                 // dca between V0 daughters
    Double_t dDCAV0ToPrimVtx   = Cascade->DcaV0ToPrimVertex();              // dca V0 to primary vertex 
    Double_t dDCABachToPrimVtx = Cascade->DcaBachToPrimVertex();            // dca bachelor to primary vertex 
    Double_t dDCABachToV0      = Cascade->DcaXiDaughters();                 // dca bachelor to V0 
    Double_t dCPAXi = Cascade->CosPointingAngleXi(dPrimVtxPos[0], dPrimVtxPos[1], dPrimVtxPos[2]); // cosine of the pointing angle of Xi          
    Double_t dCPAV0 = Cascade->CosPointingAngle(primVtx);  //??(dV0VtxPos); // cosine of the pointing angle
    Double_t dRapXi = Cascade->RapXi(); // rapidity calculated for Xi assumption
    Double_t dRapOmega = Cascade->RapOmega(); // rapidity calculated for Omega assumption
    Double_t dEtaCascade = Cascade->Eta(); // Cascade  pseudorapidity

    //Calculations for the proper lifetime      
    Double_t dXiDecayPath[3];
    for(Int_t iPos = 0; iPos < 3; iPos++)
      dXiDecayPath[iPos] = dXiVtxPos[iPos] - dPrimVtxPos[iPos]; // vector of the Xi path      
    Double_t dXiDecLen2D = TMath::Sqrt(dXiDecayPath[0] * dXiDecayPath[0] + dXiDecayPath[1] * dXiDecayPath[1]); // transverse path length R
    Double_t dROverPt = dXiDecLen2D / dPtCascade; // R/pT  
    Double_t dMROverPtXi = dMassPDGXiMinus * dROverPt; // m*R/pT
    Double_t dMROverPtOmega = dMassPDGOmegaMinus * dROverPt; // m*R/pT

    Double_t dNSigmaPosPion   = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPos, AliPID::kPion)) : 0.); // difference between measured and expected signal of the dE/dx in the TPC
    Double_t dNSigmaPosProton = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPos, AliPID::kProton)) : 0.);
    Double_t dNSigmaNegPion   = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNeg, AliPID::kPion)) : 0.);
    Double_t dNSigmaNegProton = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNeg, AliPID::kProton)) : 0.);
    Double_t dNSigmaBachPion  = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackBach, AliPID::kProton)) : 0.);
    Double_t dNSigmaBachKaon  = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackBach, AliPID::kKaon)) : 0.);

    AliAODVertex* prodVtxDaughterPos = (AliAODVertex*)(trackPos->GetProdVertex()); // production vertex of the positive daughter track
    Char_t cTypeVtxProdPos = prodVtxDaughterPos->GetType(); // type of the production vertex
    AliAODVertex* prodVtxDaughterNeg = (AliAODVertex*)(trackNeg->GetProdVertex()); // production vertex of the negative daughter track
    Char_t cTypeVtxProdNeg = prodVtxDaughterNeg->GetType(); // type of the production vertex

    Double_t dEnergy = 0;
     
    //===== Start of reconstruction cutting =====
    // 2
    // All Cascade candidates in mass range
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++;
    // Start of global cuts
    // 3
    // Reconstruction method
    if(bPrintCuts) printf("Rec: Applying cut: Reconstruction method: on-the-fly? %s\n", (fbOnFly ? "yes" : "no"));
    if(bOnFlyStatus != fbOnFly)
      continue;
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++;
    // 4
    // Tracks TPC OK
    if(bPrintCuts) printf("Rec: Applying cut: Correct charge of daughters\n");
    if(!trackNeg || !trackPos)
      continue;
    if(trackNeg->Charge() == trackPos->Charge())  // daughters have different charge?
      continue;
    if(trackNeg->Charge() != -1)  // daughters have expected charge?
      continue;
    if(trackPos->Charge() != 1)  // daughters have expected charge?
      continue;
      
    if(fbTPCRefit) {
      if(bPrintCuts) printf("Rec: Applying cut: TPC refit\n");
      if(!trackNeg->IsOn(AliAODTrack::kTPCrefit)) // TPC refit is ON?
        continue;
      if(!trackPos->IsOn(AliAODTrack::kTPCrefit))
        continue;
    }

    if(fbRejectKinks) {
      if(bPrintCuts) printf("Rec: Applying cut: Type of production vertex of daughter: No kinks\n");
      if(cTypeVtxProdNeg == AliAODVertex::kKink) // kink daughter rejection
        continue;
      if(cTypeVtxProdPos == AliAODVertex::kKink)
        continue;
    }

    if(fdCutNCrossedRowsTPCMin > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Number of TPC rows >= %g\n", fdCutNCrossedRowsTPCMin);
      if(dNRowsNeg < fdCutNCrossedRowsTPCMin) // Crossed TPC padrows
        continue;
      if(dNRowsPos < fdCutNCrossedRowsTPCMin)
        continue;
    }
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++;

    // 5
    // Daughters: transverse momentum cut
    if(fdCutPtDaughterMin > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: Daughter pt >= %g\n", fdCutPtDaughterMin);
      if((dPtDaughterNeg < fdCutPtDaughterMin) || (dPtDaughterPos < fdCutPtDaughterMin))
        continue;
      FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 6
    // Daughters: Impact parameter of V0 daughters to prim vtx
    if(bPrintCuts) printf("Rec: Applying cut: V0 daughters DCA to prim vtx: > %f\n", dDCACascadeDaughtersToPrimVtxMin);
    if((dDCAToPrimVtxNeg < dDCACascadeDaughtersToPrimVtxMin) || (dDCAToPrimVtxPos < dDCACascadeDaughtersToPrimVtxMin))
      continue;
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++;

    // 7
    // DCA between V0 daughters
    if(bPrintCuts) printf("Rec: Applying cut: DCA between V0 daughters: < %f\n", dDCACascadeV0DaughtersMax);
    if(dDCAV0Daughters > dDCACascadeV0DaughtersMax)
      continue;
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++;

    // 8
    // DCA: Bachelor to the primary vertex 
    if(bPrintCuts) printf("Rec: Applying cut: DCA bachelor to the primary vertex : < %f\n", dDCACascadeBachToPrimVtxMin);
    if(dDCABachToPrimVtx < dDCACascadeBachToPrimVtxMin)
      continue;
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++;

    // 9
    // DCA: V0 to the primary vertex
    if(bPrintCuts) printf("Rec: Applying cut: DCA V0 to the primary vertex: < %f\n", dDCACascadeV0ToPrimVtxMin);
    if(dDCAV0ToPrimVtx < dDCACascadeV0ToPrimVtxMin)
      continue;
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++;

    // 10
    // DCA: Bachelor to the V0
    if(bPrintCuts) printf("Rec: Applying cut: DCA Bachelor to the V0: < %f\n", dDCACascadeBachToV0Max);
    if(dDCABachToV0 > dDCACascadeBachToV0Max)
      continue;
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++;
    
    // 11
    // V0: Cosine of the pointing angle
    if(bPrintCuts) printf("Rec: Applying cut: V0 CPA: > %f\n", dCPACascadeV0Min);    
    if(dCPAV0 < dCPACascadeV0Min)
       continue;
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++;

    // 12
    // Xi: Cosine of the pointing angle
    if(bPrintCuts) printf("Rec: Applying cut: Cascade CPA: > %f\n", dCPACascadeMin);
    if(dCPAXi < dCPACascadeMin)
      continue;
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++; 

    // 13
    // V0: Fiducial volume
    if(bPrintCuts) printf("Rec: Applying cut: V0 Decay radius: > %f\n", dCascadeV0RadiusDecayMin);
    if(dRadiusDecayV0 < dCascadeV0RadiusDecayMin )
      continue;    
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++;

    // 14
    // Xi: Fiducial volume
    if(bPrintCuts) printf("Rec: Applying cut: Cascade Decay radius: > %f\n", dCascadeRadiusDecayMin);
    if(dRadiusDecayXi < dCascadeRadiusDecayMin )
      continue;
    FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    iCutIndex++; 

    // 15
    // Daughters: pseudorapidity cut
    if(bCutEtaV0Daughter) {
      if(bPrintCuts) printf("Rec: Applying cut: Daughter |eta|: < %f\n", dCascadeEtaDaughterMax);
      if((TMath::Abs(dEtaDaughterNeg) > dCascadeEtaDaughterMax) || (TMath::Abs(dEtaDaughterPos) > dCascadeEtaDaughterMax))
        continue;
      FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    }
    iCutIndex++;  
    // End of global cuts
     
    // Start of particle-dependent cuts
    // 16
    // Xi: rapidity cut & pseudorapidity cut
    if(fdCutRapCascadeMax > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: Cascade |y| < %g\n", fdCutRapCascadeMax);
      if(TMath::Abs(dRapXi) > fdCutRapCascadeMax) {  
        bIsCandidateXiMinus = kFALSE;
        bIsCandidateXiPlus  = kFALSE;
      }   
      if(TMath::Abs(dRapOmega) > fdCutRapCascadeMax) {  
        bIsCandidateOmegaMinus = kFALSE;
        bIsCandidateOmegaPlus  = kFALSE;
      }   
    }
    if(fdCutEtaCascadeMax) {
      if(bPrintCuts) printf("Rec: Applying cut: Cascade |eta|: < %f\n", fdCutEtaCascadeMax);
      if(TMath::Abs(dEtaCascade) > fdCutEtaCascadeMax) {
        bIsCandidateXiMinus = kFALSE;
        bIsCandidateXiPlus  = kFALSE;
        bIsCandidateOmegaMinus = kFALSE;
        bIsCandidateOmegaPlus  = kFALSE;
      }
      if(fdCutEtaCascadeMax || fdCutRapCascadeMax > 0.)
      FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 17
    // Lifetime cut
    if(fdCutNTauXMax > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: Cascade Proper lifetime: < %f\n", fdCutNTauXMax);
      if(dMROverPtXi > fdCutNTauXMax * dCTauXi) {
        bIsCandidateXiMinus = kFALSE;
        bIsCandidateXiPlus  = kFALSE;
      }
      if(dMROverPtOmega > fdCutNTauXMax * dCTauOmega) {
        bIsCandidateOmegaMinus = kFALSE;
        bIsCandidateOmegaPlus  = kFALSE;
      }
      FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    }
    iCutIndex++;
 
    // 18
    // Daughter PID
    if(fdCutNSigmadEdxMax > 0.) {
      if(fbIsPbPb && fdPtProtonPIDMax > 0.) { // Pb-Pb
        if(bPrintCuts) printf("Rec: Applying cut: Delta dE/dx (proton below %g GeV/c) < %g\n", fdPtProtonPIDMax, fdCutNSigmadEdxMax);
        if((dPtDaughterPos < fdPtProtonPIDMax) && (dNSigmaPosProton > fdCutNSigmadEdxMax)) // p+
          bIsCandidateXiMinus = kFALSE;
        if((dPtDaughterNeg < fdPtProtonPIDMax) && (dNSigmaNegProton > fdCutNSigmadEdxMax)) // p-
          bIsCandidateXiPlus = kFALSE;
      }
      else { // p-p
        if(bPrintCuts) printf("Rec: Applying cut: Delta dE/dx (both daughters): < %g\n", fdCutNSigmadEdxMax);
        if(dNSigmaBachPion > fdCutNSigmadEdxMax || dNSigmaPosProton > fdCutNSigmadEdxMax || dNSigmaNegPion > fdCutNSigmadEdxMax) // p+, pi-
          bIsCandidateXiMinus = kFALSE;
        if(dNSigmaBachPion > fdCutNSigmadEdxMax || dNSigmaNegProton > fdCutNSigmadEdxMax || dNSigmaPosPion > fdCutNSigmadEdxMax) //p-, pi+
          bIsCandidateXiPlus= kFALSE;
        //add for omega
      }   
      FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
    }
    iCutIndex++;          

    // End of particle-dependent cuts

    //===== End of reconstruction cutting =====

    if(!bIsCandidateXiMinus && !bIsCandidateXiPlus && !bIsCandidateOmegaMinus && !bIsCandidateOmegaPlus)
      continue;

    if(fbMCAnalysis) {
      AliEmcalJet *ZeroJet = 0;
      AssociateRecCascadeWithMC(Cascade, ZeroJet, bIsCandidateXiMinus, bIsCandidateXiPlus, bIsCandidateOmegaMinus, bIsCandidateOmegaPlus, iCentIndex);
    }

    Int_t uid = 0; //
    //===== Start of filling V0 spectra =====
    // iCutIndex = 19
    if(bIsCandidateXiMinus) {
      // 19
      //XiMinus candidates after cuts
      FillCandidates(dMassCascadeXi, dMassCascadeOmega, bIsCandidateXiMinus,  kFALSE,  kFALSE,  kFALSE, iCutIndex, iCentIndex);
      Double_t valueXiIncl[3] = {dMassCascadeXi, dPtCascade, dEtaCascade};
      fhnCascadeInclusiveXiMinus[iCentIndex]->Fill(valueXiIncl);
      fh1CascadeInvMassXiMinusCent[iCentIndex]->Fill(dMassCascadeXi);
      if(bIsInPeakXi){
        fhnCascadeInvMassCutXiMinus[iCentIndex]->Fill(valueXiIncl);
      }  
      uid = iXiMinusId + fNCand;
      dEnergy = Cascade->EXi();
      iNCascadeCandXiMinus++;
    }   
    if(bIsCandidateXiPlus) {
      // 19
      //XiPlus candidates after cuts
      FillCandidates(dMassCascadeXi, dMassCascadeOmega, kFALSE,  bIsCandidateXiPlus,  kFALSE,  kFALSE, iCutIndex, iCentIndex);
      Double_t valueXiIncl[3] = {dMassCascadeXi, dPtCascade, dEtaCascade};
      fhnCascadeInclusiveXiPlus[iCentIndex]->Fill(valueXiIncl);
      fh1CascadeInvMassXiPlusCent[iCentIndex]->Fill(dMassCascadeXi);
      if(bIsInPeakXi){
        fhnCascadeInvMassCutXiPlus[iCentIndex]->Fill(valueXiIncl);
      }  
      uid = iXiPlusId + fNCand;
      dEnergy = Cascade->EXi();
      iNCascadeCandXiPlus++;
    }  
    if(bIsCandidateOmegaMinus) {
      // 19
      //OmegaMinus candidates after cuts
      FillCandidates(dMassCascadeXi, dMassCascadeOmega, kFALSE,  kFALSE,  bIsCandidateOmegaMinus,  kFALSE, iCutIndex, iCentIndex);
      Double_t valueOmegaIncl[3] = {dMassCascadeOmega, dPtCascade, dEtaCascade};
      fhnCascadeInclusiveOmegaMinus[iCentIndex]->Fill(valueOmegaIncl);
      fh1CascadeInvMassOmegaMinusCent[iCentIndex]->Fill(dMassCascadeOmega);
      if(bIsInPeakOmega){
        fhnCascadeInvMassCutOmegaMinus[iCentIndex]->Fill(valueOmegaIncl);
      }  
      uid = iOmegaMinusId + fNCand;
      dEnergy = Cascade->EOmega();
      iNCascadeCandOmegaMinus++;
    }   
    if(bIsCandidateOmegaPlus) {
      // 19
      //OmegaPlus candidates after cuts
      FillCandidates(dMassCascadeOmega, dMassCascadeOmega, kFALSE,  kFALSE,  kFALSE,  bIsCandidateOmegaPlus, iCutIndex, iCentIndex);
      Double_t valueOmegaIncl[3] = {dMassCascadeOmega, dPtCascade, dEtaCascade};
      fhnCascadeInclusiveOmegaPlus[iCentIndex]->Fill(valueOmegaIncl);
      fh1CascadeInvMassOmegaPlusCent[iCentIndex]->Fill(dMassCascadeOmega);
      if(bIsInPeakOmega){
        fhnCascadeInvMassCutOmegaPlus[iCentIndex]->Fill(valueOmegaIncl);
      }  
      uid = iOmegaPlusId + fNCand;
      dEnergy = Cascade->EOmega();
      iNCascadeCandOmegaPlus++;
    }  

    if((bIsInPeakXi || bIsInPeakOmega) && uid > 0)  {  // to get rid of the situations when isinpeakcandidate is different from the iscandidate
        new ((*fCascadeCandidateArray)[fNCand]) AliAODcascade(*Cascade); //  
        ivecCascadeCandIndex.push_back(uid);
        //add the v0 vector to the fastjetwrapper
        fFastJetWrapper.AddInputVector(vecCascadeMomentum[0], vecCascadeMomentum[1], vecCascadeMomentum[2], dEnergy, uid);
        InputBgParticles.push_back(fastjet::PseudoJet(vecCascadeMomentum[0], vecCascadeMomentum[1], vecCascadeMomentum[2], dEnergy));
      
        fNCand++;

    }
    

    //===== End of filling Cascade spectra ===== 
  }
  //===== End of Cascade loop =====
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "End of Cascade loop");
  
  fh1CascadeCandPerEvent->Fill(iNCascadeCandTot);
  fh1CascadeCandPerEventCentXiMinus[iCentIndex]->Fill(iNCascadeCandXiMinus);
  fh1CascadeCandPerEventCentXiPlus[iCentIndex]->Fill(iNCascadeCandXiPlus);
  fh1CascadeCandPerEventCentOmegaMinus[iCentIndex]->Fill(iNCascadeCandOmegaMinus);
  fh1CascadeCandPerEventCentOmegaPlus[iCentIndex]->Fill(iNCascadeCandOmegaPlus);
  
  if(fCascadeCandidateArray->GetEntriesFast() > 0) {
    AddEventTracks(fCascadeCandidateArray, fTracks, InputBgParticles);
  }

  //Prepare for the jet comp for the UE subtraction 
  TClonesArray* arrayJetPerp = new TClonesArray("AliAODJet", 0); // object where the perp. cones are stored
  Int_t iNJetPerp = 0; // number of perpendicular cones
  //AliAODJet* jet = 0; // pointer to a jet
  AliAODJet* jetPerp = 0; // pointer to a perp. cone
  AliAODJet* jetRnd = 0; // pointer to a rand. cone
  TLorentzVector vecPerpPlus; // 4-momentum of perpendicular cone plus
  TLorentzVector vecPerpMinus; // 4-momentum of perpendicular cone minus

  //Background estimation 

  fastjet::JetMedianBackgroundEstimator bge;
  fastjet::Subtractor subtr(&bge);
  bge.set_selector(selectorBG); 
  fastjet::JetDefinition jetDefBG(fastjet::kt_algorithm, fdBgRadius, fastjet::pt_scheme); //define the kT jet finding which will do the average background estimation
  fastjet::AreaDefinition areaDefBG(fastjet::active_area_explicit_ghosts);
  fastjet::ClusterSequenceArea cluster_seq_BG(InputBgParticles, jetDefBG, areaDefBG);
  std::vector<fastjet::PseudoJet> jetsBG = sorted_by_pt(selectorBG(cluster_seq_BG.inclusive_jets())); //find the kT jets
  if (jetsBG.size() > 0) {
    bge.set_jets(jetsBG);  // give the kT jets to the background estimator
    bge.rho();
  }


  //Run fj wrapper and loop over fj jets
  fFastJetWrapper.Run(); 

  std::vector<fastjet::PseudoJet> vJetsIncl = fFastJetWrapper.GetInclusiveJets(); 
  if(vJetsIncl.size() > 0) 
	  fh1EventCounterCut->Fill(10);
  
  // sort jets according to jet pt
  static Int_t indexes[9999] = {-1};
  GetSortedArray(indexes, vJetsIncl); 
  AliDebug(1,Form("%i jets found", (Int_t)vJetsIncl.size()));
  
  Int_t iJetCount = 0;
  Int_t iCascadeJets = 0;
  fastjet::PseudoJet jetSub;

  for (UInt_t ijet = 0; ijet < vJetsIncl.size(); ++ijet) {
    Int_t ij = indexes[ijet];
    AliDebug(3,Form("Jet pt = %f, area = %f", vJetsIncl[ij].perp(), fFastJetWrapper.GetJetArea(ij)));
    if (vJetsIncl[ij].perp() < fdMinJetPt) continue;
    if (fFastJetWrapper.GetJetArea(ij) < fdMinJetArea) continue;
    if ((vJetsIncl[ij].eta() < fdJetEtaMin) || (vJetsIncl[ij].eta() > fdJetEtaMax) ||
        (vJetsIncl[ij].phi() < fdJetPhiMin) || (vJetsIncl[ij].phi() > fdJetPhiMax))
      continue;

    jetSub = subtr(vJetsIncl[ij]);
    if(jetSub == 0) 
      continue;  
    
    if(jetSub.perp() < fdCutPtJetMin) // selection of high-pt jets, needs to be applied on the pt after bg subtraction
      continue;
      
    fh1NCascadesInJetStats->Fill(7);  //N of jets before cuts (after bg sub)
    
    if(jetSub.perp() < fdCutPtJetMin) // selection of high-pt jets, needs to be applied on the pt after bg subtraction
      continue;
    fh1NCascadesInJetStats->Fill(8);
    
    if(jetSub.area() < dAreaPercJetMin) //selection of the jets with area bigger than the cut (cut*pi*R2)
      continue;
    fh1NCascadesInJetStats->Fill(9);
    
    std::vector<fastjet::PseudoJet> constituents(jetSub.constituents());  //(fFastJetWrapper.GetJetConstituents(ij));
    Int_t iNConstit = constituents.size();
    Double_t dMaxTrPt = 0;
    for(Int_t ic = 0; ic < iNConstit; ic++) {
		  Double_t dtrpt = constituents[ic].perp();
      if(dtrpt > dMaxTrPt) dMaxTrPt = dtrpt;      
    }  
    if(dMaxTrPt < fdCutPtTrackJetMin)             // selection of jets with high leading track pt
      continue;                                           
    fh1NCascadesInJetStats->Fill(10);

    fh2PtJetPtTrackLeading[iCentIndex]->Fill(jetSub.perp(), dMaxTrPt); // pt_jet vs pt of leading jet track   

    AliEmcalJet *jet = new ((*fJets)[iJetCount])
    		          AliEmcalJet(jetSub.perp(), jetSub.eta(), jetSub.phi(), jetSub.m());
    jet->SetLabel(ij);
    fastjet::PseudoJet area(jetSub.area_4vector()); //(fFastJetWrapper.GetJetAreaVector(ij));
    jet->SetArea(area.perp());
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaE(area.E());

    vecPerpPlus.SetPtEtaPhiM(jetSub.perp(), jetSub.eta(), jetSub.phi(), 0.);
    vecPerpMinus.SetPtEtaPhiM(jetSub.perp(), jetSub.eta(), jetSub.phi(), 0.);
    vecPerpPlus.RotateZ(TMath::Pi() / 2.); // rotate vector by +90 deg around z
    vecPerpMinus.RotateZ(-TMath::Pi() / 2.); // rotate vector by -90 deg around z
    new((*arrayJetPerp)[iNJetPerp++]) AliAODJet(vecPerpPlus); // write perp. cone to the array
    new((*arrayJetPerp)[iNJetPerp++]) AliAODJet(vecPerpMinus); // write perp. cone to the array
     

    fh1PtJet[iCentIndex]->Fill(jet->Pt()); // pt spectrum of selected jets
    fh1EtaJet[iCentIndex]->Fill(jet->Eta()); // eta spectrum of selected jets
    fh2EtaPtJet[iCentIndex]->Fill(jet->Eta(), jet->Pt()); // eta-pT spectrum of selected jets
    fh1PhiJet[iCentIndex]->Fill(jet->Phi()); // phi spectrum of selected jets
    Double_t dAreaExcluded = TMath::Pi() * dRadiusExcludeCone * dRadiusExcludeCone; // area of the cone
    dAreaExcluded -= AreaCircSegment(dRadiusExcludeCone, fdCutEtaCascadeMax - jet->Eta()); // positive eta overhang
    dAreaExcluded -= AreaCircSegment(dRadiusExcludeCone, fdCutEtaCascadeMax + jet->Eta()); // negative eta overhang
    fh1AreaExcluded->Fill(iCentIndex, dAreaExcluded);
    fh1NCascadesInJetStats->Fill(0);
    if(fbMCAnalysis) 
      fh1MCStats->Fill(1); //Reconstructed jets
    iJetCount++;
   //printf("JetPt: %f, sub jetPt: %f \n", vJetsIncl[ij].perp(), jetSub.perp());
  
    Int_t uid   = -1;
    Int_t ixm = 0;
    Int_t ixp = 0;
    Int_t iom = 0; 
    Int_t iop = 0;

    Int_t index = 0;
    AliAODcascade* jetcascade = 0;
    
    for(Int_t ic = 0; ic < iNConstit; ic++){
		  uid = constituents[ic].user_index();
		  
      if(uid > iXiMinusId && uid < iXiPlusId) { //XiMinus
			  index = uid-iXiMinusId; //will give the id of the Cascade particle 
		    jetcascade = (AliAODcascade*)fCascadeCandidateArray->At(index); 
		    Double_t valueXMInJet[4] = {jetcascade->MassXi(), TMath::Sqrt(jetcascade->Pt2Xi()), jetcascade->Eta(), jet->Pt()};
        fhnCascadeInJetXiMinus[iCentIndex]->Fill(valueXMInJet);
        fh1NCascadesInJetStats->Fill(2);
        if(fbMCAnalysis) {
          AssociateRecCascadeWithMC(jetcascade, jet, true, false, false, false, iCentIndex);
	      }  
		   ixm++;
		  }			  
		  if(uid > iXiPlusId && uid < iOmegaMinusId) {  //XiPlus
		    index = uid-iXiPlusId; 
		    jetcascade = (AliAODcascade*)fCascadeCandidateArray->At(index); 
		    Double_t valueXPInJet[4] = {jetcascade->MassXi(), TMath::Sqrt(jetcascade->Pt2Xi()), jetcascade->Eta(), jet->Pt()};
        fhnCascadeInJetXiPlus[iCentIndex]->Fill(valueXPInJet);
        fh1NCascadesInJetStats->Fill(3);
        if(fbMCAnalysis) {
          AssociateRecCascadeWithMC(jetcascade, jet, false, true, false, false, iCentIndex);
	      }
	    
        ixp++;
		  }
		  if(uid > iOmegaMinusId && uid < iOmegaPlusId) { //OmegaMinus
		    index = uid-iOmegaMinusId;  
		    jetcascade = (AliAODcascade*)fCascadeCandidateArray->At(index); 
		  	Double_t valueOMInJet[4] = {jetcascade->MassOmega(), TMath::Sqrt(jetcascade->Pt2Xi()), jetcascade->Eta(), jet->Pt()};
        fhnCascadeInJetOmegaMinus[iCentIndex]->Fill(valueOMInJet);
        fh1NCascadesInJetStats->Fill(4);
        if(fbMCAnalysis) {
          AssociateRecCascadeWithMC(jetcascade, jet, false, false, true, false, iCentIndex);  
	      }
        iom++;
		  }  
 		  if(uid > iOmegaPlusId) { //OmegaPlus
		    index = uid-iOmegaPlusId;  
		    jetcascade = (AliAODcascade*)fCascadeCandidateArray->At(index); 
		  	Double_t valueOPInJet[4] = {jetcascade->MassOmega(), TMath::Sqrt(jetcascade->Pt2Xi()), jetcascade->Eta(), jet->Pt()};
        fhnCascadeInJetOmegaPlus[iCentIndex]->Fill(valueOPInJet);
        fh1NCascadesInJetStats->Fill(5);
        if(fbMCAnalysis) {
          AssociateRecCascadeWithMC(jetcascade, jet, false, false, false, true, iCentIndex);  
	      }
        iop++;
		  } 
     } 
    Int_t isum = ixm + ixp + iom + iop;
    if(isum > 0 ) {
	    iCascadeJets++;
	    fh1NCascadesInJetStats->Fill(1);
	  }
	  if(isum > 1)
	    fh1NCascadesInJetStats->Fill(7);	  		  
  } 
  if(iJetCount != 0) 
	  fh1EventCounterCut->Fill(11);
  if(iCascadeJets > 0) 
    fh1EventCounterCut->Fill(12); // events with jets with V0s
  fh1NJetPerEvent[iCentIndex]->Fill(iJetCount);
  fh1NCascadeJetPerEvent[iCentIndex]->Fill(iCascadeJets);  
 
  if(iJetCount) {
    jetRnd = GetRandomCone( fJets, dCutEtaJetMax, 2 * fdDistanceCascadeJetMax); //max jet pseudorap (fdCutEtaCascadeMax - fdDistanceCascadeJetMax)
    if(jetRnd) {
      fh1NRndConeCent->Fill(iCentIndex);
      fh2EtaPhiRndCone[iCentIndex]->Fill(jetRnd->Eta(), jetRnd->Phi());
    }
  }
  //Find cascades in the perp and random jet cones 
  if(fCascadeCandidateArray->GetEntriesFast() > 0) {
    for(Int_t icasc = 0; icasc < fCascadeCandidateArray->GetEntriesFast(); icasc++){
      AliAODcascade* casccand = (AliAODcascade*)fCascadeCandidateArray->At(icasc);
      Int_t iind = ivecCascadeCandIndex[icasc] - icasc;  
      Double_t dCascadeMassXi = casccand->MassXi();
      Double_t dCascadeMassOmega = casccand->MassOmega();
      Double_t dPtv =  TMath::Sqrt(casccand->Pt2Xi());
      Double_t dEtav = casccand->Eta();

      //Selection of cascade in no-jet events
      if(!iJetCount) {
        if(iind == iXiMinusId) {
          Double_t valueXiMNoJet[3] = {dCascadeMassXi, dPtv,dEtav};
          fhnCascadeNoJetXiMinus[iCentIndex]->Fill(valueXiMNoJet);
        }
        if(iind == iXiPlusId) {
          Double_t valueXiPNoJet[3] = {dCascadeMassXi, dPtv,dEtav};
          fhnCascadeNoJetXiPlus[iCentIndex]->Fill(valueXiPNoJet);
        }
        if(iind == iOmegaMinusId) {
          Double_t valueOmMNoJet[3] = {dCascadeMassOmega, dPtv,dEtav};
          fhnCascadeNoJetOmegaMinus[iCentIndex]->Fill(valueOmMNoJet);
        }
        if(iind == iOmegaPlusId) {
          Double_t valueOmPNoJet[3] = {dCascadeMassOmega, dPtv,dEtav};
          fhnCascadeNoJetOmegaPlus[iCentIndex]->Fill(valueOmPNoJet);
        }        
      }
      // Selection of Cascades outside jet cones
      if(!OverlapWithJets(fJets, casccand, 2 * fdDistanceCascadeJetMax)) { 
        if(iind == iXiMinusId) {
          Double_t valueXiMOutJC[3] = {dCascadeMassXi, dPtv,dEtav};
          fhnCascadeOutJetXiMinus[iCentIndex]->Fill(valueXiMOutJC);
        }
        if(iind == iXiPlusId) {
          Double_t valueXiPOutJet[3] =  {dCascadeMassXi, dPtv,dEtav};
          fhnCascadeOutJetXiPlus[iCentIndex]->Fill(valueXiPOutJet);
        }
        if(iind == iOmegaMinusId) {
          Double_t valueOmMOutJet[3] =  {dCascadeMassOmega, dPtv,dEtav};
          fhnCascadeOutJetOmegaMinus[iCentIndex]->Fill(valueOmMOutJet);
        }
        if(iind == iOmegaPlusId) {
          Double_t valueOmPOutJet[3] =  {dCascadeMassOmega, dPtv,dEtav};
          fhnCascadeOutJetOmegaPlus[iCentIndex]->Fill(valueOmPOutJet);
        }
      }
      // Selection of Cascades in perp. cones 
      for(Int_t iJet = 0; iJet < iNJetPerp; iJet++) {
        jetPerp = (AliAODJet*)arrayJetPerp->At(iJet); // load a jet in the list   
        if(IsParticleInCone(casccand, jetPerp, fdDistanceCascadeJetMax)) { // Cascade in perp. cone
          if(iind == iXiMinusId) {
          Double_t valueXiMInPC[4] = {dCascadeMassXi, dPtv, dEtav, jetPerp->Pt()};
          fhnCascadeInPerpXiMinus[iCentIndex]->Fill(valueXiMInPC);
          }
          if(iind == iXiPlusId) {
          Double_t valueXiPInPC[4] = {dCascadeMassXi, dPtv, dEtav, jetPerp->Pt()};
          fhnCascadeInPerpXiPlus[iCentIndex]->Fill(valueXiPInPC);
          }
          if(iind == iOmegaMinusId) {
          Double_t valueOmMInPC[4] = {dCascadeMassOmega, dPtv, dEtav, jetPerp->Pt()};
          fhnCascadeInPerpOmegaMinus[iCentIndex]->Fill(valueOmMInPC);
          }
          if(iind == iOmegaPlusId) {
          Double_t valueOmPInPC[4] = {dCascadeMassOmega, dPtv, dEtav, jetPerp->Pt()};
          fhnCascadeInPerpOmegaPlus[iCentIndex]->Fill(valueOmPInPC);
          }
          break;
        }
      }
      // Selection of Cascades in random cones
      if(jetRnd) {
        if(IsParticleInCone(Cascade, jetRnd, fdDistanceCascadeJetMax)) {
          if(iind == iXiMinusId) {
            Double_t valueXiMInRnd[3] = {dCascadeMassXi, dPtv,dEtav};
            fhnCascadeInRndXiMinus[iCentIndex]->Fill(valueXiMInRnd);
          }
          if(iind == iXiPlusId) {
            Double_t valueXiPInRnd[3] = {dCascadeMassXi, dPtv,dEtav};
            fhnCascadeInRndXiPlus[iCentIndex]->Fill(valueXiPInRnd);
          }
          if(iind == iOmegaMinusId) {
            Double_t valueOmMInRnd[3] = {dCascadeMassOmega, dPtv,dEtav};
            fhnCascadeInRndOmegaMinus[iCentIndex]->Fill(valueOmMInRnd);
          }
          if(iind == iOmegaPlusId) {
            Double_t valueOmPInRnd[3] = {dCascadeMassOmega, dPtv,dEtav};
            fhnCascadeInRndOmegaPlus[iCentIndex]->Fill(valueOmPInRnd);
          }          
        }
      }
    }
  }

   // Spectra of generated particles
  if(fbMCAnalysis) {
    GeneratedMCParticles(fTracks, iCentIndex);
  }

  PostData(1, fOutputListStd);
  PostData(2, fOutputListStdJets);
  PostData(3, fOutputListMC);
    
  if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "End");

  return kFALSE; // Must be false to avoid calling PostData from AliAnalysisTaskEmcal. Otherwise, slot 1 is not stored.
}

void AliAnalysisTaskCascadesInJets::FillCandidates(Double_t mXi, Double_t mOmega, Bool_t isXiMinus, Bool_t isXiPlus, Bool_t isOmegaMinus, Bool_t isOmegaPlus, Int_t iCut, Int_t iCent)
{
  if(isXiMinus) {
    fh1CascadeCounterCentXiMinus[iCent]->Fill(iCut);
    fh1CascadeInvMassXiMinusAll[iCut]->Fill(mXi);
  }
  if(isXiPlus) {
    fh1CascadeCounterCentXiPlus[iCent]->Fill(iCut);
    fh1CascadeInvMassXiPlusAll[iCut]->Fill(mXi);
  } 
  if(isOmegaMinus) {
    fh1CascadeCounterCentOmegaMinus[iCent]->Fill(iCut);
    fh1CascadeInvMassOmegaMinusAll[iCut]->Fill(mOmega);
  }
  if(isOmegaPlus) {
    fh1CascadeCounterCentOmegaPlus[iCent]->Fill(iCut);
    fh1CascadeInvMassOmegaPlusAll[iCut]->Fill(mOmega);
  } 
}

Bool_t AliAnalysisTaskCascadesInJets::IsParticleInCone(const AliVParticle* part1, const AliVParticle* part2, Double_t dRMax) const
{
// decides whether a particle is inside a jet cone
  if(!part1 || !part2)
    return kFALSE;

  TVector3 vecMom2(part2->Px(), part2->Py(), part2->Pz());
  TVector3 vecMom1(part1->Px(), part1->Py(), part1->Pz());
  Double_t dR = vecMom2.DeltaR(vecMom1); // = sqrt(dEta*dEta+dPhi*dPhi)
  if(dR < dRMax) // momentum vectors of part1 and part2 are closer than dRMax
    return kTRUE;
  return kFALSE;
}
Bool_t AliAnalysisTaskCascadesInJets::OverlapWithJets(const TClonesArray* array, const AliVParticle* part, Double_t dDistance) const
{  //there is a problem with this function?? 
   // decides whether a cone overlaps with other jets
  if(!part) {
    AliError("No particle!");
    return kFALSE;
  }
  if(!array) {
    AliError("No jet array!");
    return kFALSE;
  }
  Int_t iNJets = array->GetEntriesFast();
  if(iNJets <= 0) {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Warning: No jets");
    return kFALSE;
  }
  AliVParticle* jet = 0;
  for(Int_t iJet = 0; iJet < iNJets; iJet++) {
    jet = (AliVParticle*)array->At(iJet);
    if(!jet) {
      AliError(Form("Failed to load jet %d/%d!", iJet, iNJets));
      continue;
    }
    if(IsParticleInCone(part, jet, dDistance))
      return kTRUE;
  }
  return kFALSE;
}
AliAODJet* AliAnalysisTaskCascadesInJets::GetRandomCone(const TClonesArray* array, Double_t dEtaConeMax, Double_t dDistance) const
{
  // generate a random cone which does not overlap with selected jets
  if(fDebug > 3) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Generating random cone...");
  TLorentzVector vecCone;
  AliAODJet* part = 0;
  Double_t dEta, dPhi;
  Int_t iNTrialsMax = 10;
  Bool_t bStatus = kFALSE;
  for(Int_t iTry = 0; iTry < iNTrialsMax; iTry++) {
    if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Try %d", iTry));
    dEta = dEtaConeMax * (2 * fRandom->Rndm() - 1.); // random eta in [-dEtaConeMax,+dEtaConeMax]
    dPhi = TMath::TwoPi() * fRandom->Rndm(); // random phi in [0,2*Pi]
    vecCone.SetPtEtaPhiM(1., dEta, dPhi, 0.);
    part = new AliAODJet(vecCone);
    if(!OverlapWithJets(array, part, dDistance)) {
      bStatus = kTRUE;
      if(fDebug > 1) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Random cone successfully generated");
      break;
    }
    else
      delete part;
  }
  if(!bStatus) {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Failed to find a random cone");
    part = 0;
  }
  return part;
}

Double_t AliAnalysisTaskCascadesInJets::AreaCircSegment(Double_t dRadius, Double_t dDistance) const
{
// calculate area of a circular segment defined by the circle radius and the (oriented) distance between the secant line and the circle centre
  Double_t dEpsilon = 1e-2;
  Double_t dR = dRadius;
  Double_t dD = dDistance;
  if(TMath::Abs(dR) < dEpsilon)
  {
    AliError(Form("Too small radius: %g < %g!", dR, dEpsilon));
    return 0.;
  }
  if(dD >= dR)
    return 0.;
  if(dD <= -dR)
    return TMath::Pi() * dR * dR;
  return dR * dR * TMath::ACos(dD / dR) - dD * TMath::Sqrt(dR * dR - dD * dD);
}

Bool_t AliAnalysisTaskCascadesInJets::IsSelectedForAnalysis() //Event selection
{
  if(!fbIsPbPb) {
    if(fAODIn->IsPileupFromSPD()) {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "SPD pile-up");
      return kFALSE;
    }
    fh1EventCounterCut->Fill(2); // not pile-up from SPD
  }
  AliAODVertex* vertex = fAODIn->GetPrimaryVertex();
  if(!vertex) {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "No vertex");
    return kFALSE;
  }
  if(fiNContribMin > 0) {
    if(vertex->GetNContributors() < fiNContribMin) {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Not enough contributors, %d", vertex->GetNContributors()));
      return kFALSE;
    }
    fh1EventCounterCut->Fill(3); // enough contributors
  }
  if(fbUseIonutCut) {
    if(!fEventCutsStrictAntipileup.AcceptEvent(fAODIn)) {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Ionut's cut");
      return kFALSE;
    }
    fh1EventCounterCut->Fill(4); // Ionut's cut
  }
  Double_t zVertex = vertex->GetZ();
  if(fdCutVertexZ > 0.) {
    if(TMath::Abs(zVertex) > fdCutVertexZ) {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Cut on z, %g", zVertex));
      return kFALSE;
    }
    fh1EventCounterCut->Fill(5); // PV z coordinate within range
  }
  if(fdCutDeltaZMax > 0.) { // cut on |delta z| between SPD vertex and nominal primary vertex
    AliAODVertex* vertexSPD = fAODIn->GetPrimaryVertexSPD();
    if(!vertexSPD) {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "No SPD vertex");
      return kFALSE;
    }
    Double_t zVertexSPD = vertexSPD->GetZ();
    if(TMath::Abs(zVertex - zVertexSPD) > fdCutDeltaZMax) {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Cut on Delta z = %g - %g = %g", zVertex, zVertexSPD, zVertex - zVertexSPD));
      return kFALSE;
    }
    fh1EventCounterCut->Fill(6); // delta z within range
  }
  if(fdCutVertexR2 > 0.) {
    Double_t xVertex = vertex->GetX();
    Double_t yVertex = vertex->GetY();
    Double_t radiusSq = yVertex * yVertex + xVertex * xVertex;
    if(radiusSq > fdCutVertexR2) {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Cut on r, %g", radiusSq));
      return kFALSE;
    }
    fh1EventCounterCut->Fill(7); // radius within range
  }
  if(fbIsPbPb) {
    if(fbUseMultiplicity) { // from https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentralityCodeSnippets
      AliMultSelection* MultSelection = 0x0;
      MultSelection = (AliMultSelection*)fAODIn->FindListObject("MultSelection");
      if(!MultSelection) {
        AliWarning("AliMultSelection object not found!");
        return kFALSE;
      }
      fdCentrality = MultSelection->GetMultiplicityPercentile("V0M");
    }
    else
      fdCentrality = ((AliVAODHeader*)fAODIn->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
    if(fdCentrality < 0) {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Negative centrality");
      return kFALSE;
    }
    if((fdCutCentHigh < 0) || (fdCutCentLow < 0) || (fdCutCentHigh > 100) || (fdCutCentLow > 100) || (fdCutCentLow > fdCutCentHigh)) {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Wrong centrality limits");
      return kFALSE;
    }
    if((fdCentrality < fdCutCentLow) || (fdCentrality > fdCutCentHigh)) {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Centrality cut, %g", fdCentrality));
      return kFALSE;
    }
  }
  return kTRUE;
}

Int_t AliAnalysisTaskCascadesInJets::GetCentralityBinIndex(Double_t centrality)
{
// returns index of the centrality bin corresponding to the provided value of centrality
  if(centrality < 0 || centrality > fgkiCentBinRanges[fgkiNBinsCent - 1])
    return -1;
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
  {
    if(centrality <= fgkiCentBinRanges[i])
      return i;
  }
  return -1;
}

Int_t AliAnalysisTaskCascadesInJets::GetCentralityBinEdge(Int_t index)
{
// returns the upper edge of the centrality bin corresponding to the provided value of index
  if(index < 0 || index >= fgkiNBinsCent)
    return -1;
  return fgkiCentBinRanges[index];
}

TString AliAnalysisTaskCascadesInJets::GetCentBinLabel(Int_t index)
{
// get string with centrality range for given bin
  TString lowerEdge = ((index == 0) ? "0" : Form("%d", GetCentralityBinEdge(index - 1)));
  TString upperEdge = Form("%d", GetCentralityBinEdge(index));
  return Form("%s-%s %%", lowerEdge.Data(), upperEdge.Data());
}

Bool_t AliAnalysisTaskCascadesInJets::GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const
{
  static Float_t pt[9999] = {0};

  const Int_t n = (Int_t)array.size();

  if (n < 1)
    return kFALSE;

  for (Int_t i = 0; i < n; i++)
    pt[i] = array[i].perp();

  TMath::Sort(n, pt, indexes);

  return kTRUE;
}


Bool_t AliAnalysisTaskCascadesInJets::IsFromGoodGenerator(Int_t index)
{
  if(!fEventMC) {
    AliError("No MC event!");
    return kFALSE;
  }
  if(fsGeneratorName.Length()) {
    TString sGenName = "";
    Bool_t bCocktailOK = fEventMC->GetCocktailGenerator(index, sGenName);
    if(!bCocktailOK || !sGenName.Contains(fsGeneratorName.Data()))
      return kFALSE;
  }
  return kTRUE;
}

void AliAnalysisTaskCascadesInJets::AddEventTracks(TClonesArray* coll, TClonesArray* tracks, std::vector<fastjet::PseudoJet>& VectorBgPart)
{ 
  // Add event tracks to a collection that already contains the Cascade candidates, excluding the daughters of the Cascade candidates

  if (!tracks) {
	printf("Track container was not found. Function AddEventTracks will not run! \n");
	return;
  }
  
  TObjArray allDaughters(10);
  allDaughters.SetOwner(kFALSE);

  TIter next(coll);
  AliAODcascade* cascadepart = 0;
  while ((cascadepart = static_cast<AliAODcascade*>(next()))) {
    AliDebug(2, Form("Found a Cascade candidtate with pT = %.3f, eta = %.3f, phi = %.3f \n", cascadepart->Pt(), cascadepart->Eta(), cascadepart->Phi()));
    if (cascadepart) AddDaughters(cascadepart, allDaughters);
  }

  AliVTrack* track = 0;
  Int_t n = coll->GetEntriesFast();
   
  TIter nexttr(tracks);
  Int_t numbtrack = 0;
  Int_t nadded = 0;
  Int_t nexcluded = 0;
  while ((track = static_cast<AliVTrack*>(nexttr()))) { 
	  numbtrack++;   
	  if (track->Pt() < fdJetTrackPtMin || TMath::Abs(track->Eta()) > fdJetTrackEtaMax) { //condition for the minimum track pt  and tracj eta
      nexcluded++;
      continue; 
    }  
    if (allDaughters.Remove(track) == 0) {
      //adding track to the fastjet
      fFastJetWrapper.AddInputVector(track->Px(), track->Py(), track->Pz(), track->E(), n);
      VectorBgPart.push_back(fastjet::PseudoJet(track->Px(), track->Py(), track->Pz(), track->E()));
      //InputBgParticles.push_back(fastjet::PseudoJet(track->Px(), track->Py(), track->Pz(), track->E()));
      n++;
      nadded++;
      AliDebug(2, Form("Track %d (pT = %.3f, eta = %.3f, phi = %.3f) is included", numbtrack, track->Pt(), track->Eta(), track->Phi()));
    }
    else {
		nexcluded++;
		AliDebug(2, Form("Track %d (pT = %.3f, eta = %.3f, phi = %.3f) is excluded", numbtrack, track->Pt(), track->Eta(), track->Phi()));
    }
  } 
  //printf("There were %d tracks, %d were added to the fj and %d excluded. \n", numbtrack, nadded, nexcluded);
}

Double_t AliAnalysisTaskCascadesInJets::AddDaughters(AliAODRecoDecay* cand, TObjArray& daughters)
{
  // Add all the dauthers of cand in an array. Follows all the decay cascades.

  Int_t n = cand->GetNDaughters();
  //printf("AddDaughters: the number of daughters is %d \n", n);

  Int_t ntot = 0;
  Double_t pt = 0;
  for (Int_t i = 0; i < n; i++) {
    AliVTrack* track = dynamic_cast<AliVTrack*>(cand->GetDaughter(i));
    if (!track) {
      printf("No daughter from the candidate in AddDaughters, continue \n ");		
     continue;
    }
    AliAODRecoDecay* cand2 = dynamic_cast<AliAODRecoDecay*>(track);

    if (cand2) {
      //printf("cand2 true (call adddaughter for cand2(has its own daughter)), Daughter pT = %.3f --> \n", track->Pt());
      pt += AddDaughters(cand2, daughters);
    }
    else {
      if (!track->InheritsFrom("AliAODTrack")) {
        printf("Warning: One of the daughters is not of type 'AliAODTrack' nor 'AliAODRecoDecay'.\n");
        continue;
      }
      //printf("cand2 false, will not have daughters, add to array, Daughter pT = %.3f\n", track->Pt());
      daughters.AddLast(track);
      pt += track->Pt();
      ntot++;
    }
  }
  //printf("Total pt of the daughters = %.3f \n", pt);

  return pt;
}  

void AliAnalysisTaskCascadesInJets::AddEventTracksMC(TClonesArray* coll, TClonesArray* tracks, std::vector<fastjet::PseudoJet>& VectorBgPartMC)
{ 
  if (!tracks) {
	printf("Track container was not found. Function AddEventTracks will not run! \n");
	return;
  }

  std::vector<Int_t>  vecDaughterLabels; //vector with labels of daughter particles

  TIter next(coll);
  AliAODMCParticle* cascadepart = 0;
  while ((cascadepart = static_cast<AliAODMCParticle*>(next()))) {
    AliDebug(2, Form("Found a MC generated Cascade with pT = %.3f, eta = %.3f, phi = %.3f \n", cascadepart->Pt(), cascadepart->Eta(), cascadepart->Phi())); 
    Int_t n = cascadepart->GetNDaughters();  //why 3 daughters occur? 
    cout << "Daughters: " << n << endl;
    for (Int_t i = 0; i < n; i++) {
      Int_t iDLabel = cascadepart->GetDaughterLabel(i);
      vecDaughterLabels.push_back(iDLabel);
    } 
  }
  
  AliAODTrack* track = 0;
  Int_t iN = coll->GetEntriesFast();
  //Int_t inlabels = int(vecDaughterLabels.size());
   
  TIter nexttr(tracks);
  Int_t numbtrack = 0;
  Int_t nadded = 0;
  Int_t nexcl = 0;
  while ((track = static_cast<AliAODTrack*>(nexttr()))) { //here add condition for the track pt fdJetTrackPtMin!!!
	numbtrack++;    
	if (track->Pt() < fdJetTrackPtMin || TMath::Abs(track->Eta()) > fdJetTrackEtaMax) //condition for the minimum track pt  and tracj eta
      continue; 

    Int_t iTrackLabel = track->GetLabel();
    if (std::find(vecDaughterLabels.begin(), vecDaughterLabels.end(), iTrackLabel) != vecDaughterLabels.end()) {
	    nexcl++;
	    continue;	 
    }
    //did not find track in the daughters, adding track to the fastjet
    fFastJetWrapperMCGen.AddInputVector(track->Px(), track->Py(), track->Pz(), track->E(), iN);
    VectorBgPartMC.push_back(fastjet::PseudoJet(track->Px(), track->Py(), track->Pz(), track->E()));
    iN++;
    nadded++;

  } 
  //printf("There were %i Cascades, %i daughters, %d tracks, %d were added to the fj, excluded : %i. \n", coll->GetEntriesFast(), inlabels, numbtrack, nadded, nexcl);
}


Bool_t AliAnalysisTaskCascadesInJets::AssociateRecCascadeWithMC( AliAODcascade* cascpart, AliEmcalJet *xjet, Bool_t bIsXiMinus, Bool_t bIsXiPlus, Bool_t bIsOmegaMinus, Bool_t bIsOmegaPlus, Int_t iCent)
{
  if(!cascpart) {
    printf("CANNOT FIND Cascade!!!\n");	
    return kFALSE;
  }  
  TClonesArray* MCPartArray = 0; // array particles in the MC event
  AliAODMCHeader* MCHeader = 0; // MC header
  Int_t iNTracksMC = 0; // number of MC tracks
  Double_t dPrimVtxMCX = 0., dPrimVtxMCY = 0., dPrimVtxMCZ = 0.; // position of the MC primary vertex
  
  // Simulation info
  MCPartArray = (TClonesArray*)fAODIn->FindListObject(AliAODMCParticle::StdBranchName());
  if(!MCPartArray) {
    printf("No MC array found!");
    return kFALSE;
  }
  if(fDebug > 2) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "MC array found");
  iNTracksMC = MCPartArray->GetEntriesFast();
  MCHeader = (AliAODMCHeader*)fAODIn->FindListObject(AliAODMCHeader::StdBranchName());
  if(!MCHeader) {
    printf("No MC header found!");
    return kFALSE;
  }
  // get position of the MC primary vertex
  dPrimVtxMCX = MCHeader->GetVtxX();
  dPrimVtxMCY = MCHeader->GetVtxY();
  dPrimVtxMCZ = MCHeader->GetVtxZ();

  Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
  Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
  Double_t dMassPDGXi = TDatabasePDG::Instance()->GetParticle(kXiMinus)->Mass();       
  Double_t dMassPDGOmega = TDatabasePDG::Instance()->GetParticle(kOmegaMinus)->Mass(); 



  const AliAODTrack* trB = (AliAODTrack*)cascpart->GetDecayVertexXi()->GetDaughter(0); // bachelor daughter track    
  const AliAODTrack* trP = (AliAODTrack*)cascpart->GetDaughter(0); // positive daughter track
  const AliAODTrack* trN = (AliAODTrack*)cascpart->GetDaughter(1); // negative daughter track

  Int_t iLabelBach = TMath::Abs(trB->GetLabel()); 
  Int_t iLabelPos = TMath::Abs(trP->GetLabel());
  Int_t iLabelNeg = TMath::Abs(trN->GetLabel());

  Double_t dMXi = cascpart->MassXi();
  Double_t dMOmega = cascpart->MassOmega();

  Double_t dptcasc = TMath::Sqrt(cascpart->Pt2Xi()); // transverse momentum of Cascade
  
  // Make sure MC daughters are in the array range
  if((iLabelNeg < 0) || (iLabelNeg >= iNTracksMC) || (iLabelPos < 0) || (iLabelPos >= iNTracksMC) || (iLabelBach < 0) || (iLabelBach >= iNTracksMC))
    return kFALSE;
  // Get MC particles corresponding to reconstructed daughter tracks
  AliAODMCParticle* particleMCDaughterBach = (AliAODMCParticle*)MCPartArray->At(iLabelBach);
  AliAODMCParticle* particleMCDaughterNeg = (AliAODMCParticle*)MCPartArray->At(iLabelNeg);
  AliAODMCParticle* particleMCDaughterPos = (AliAODMCParticle*)MCPartArray->At(iLabelPos);

  if(!particleMCDaughterNeg || !particleMCDaughterPos || !particleMCDaughterBach)
    return kFALSE;
  // Make sure MC daughter particles are not physical primary
  if((particleMCDaughterNeg->IsPhysicalPrimary()) || (particleMCDaughterPos->IsPhysicalPrimary()) || (particleMCDaughterBach->IsPhysicalPrimary()))
    return kFALSE;
  // Get identities of MC daughter particles
  Int_t iPdgCodeDaughterBach = particleMCDaughterBach->GetPdgCode();
  Int_t iPdgCodeDaughterPos = particleMCDaughterPos->GetPdgCode();
  Int_t iPdgCodeDaughterNeg = particleMCDaughterNeg->GetPdgCode();
  // Get index of the mother particle for each MC daughter particle
  Int_t iIndexMotherBach = particleMCDaughterBach->GetMother();
  Int_t iIndexMotherPos = particleMCDaughterPos->GetMother();
  Int_t iIndexMotherNeg = particleMCDaughterNeg->GetMother();
  if((iIndexMotherNeg < 0) || (iIndexMotherNeg >= iNTracksMC) || (iIndexMotherPos < 0) || (iIndexMotherPos >= iNTracksMC) || (iIndexMotherBach < 0) || (iIndexMotherBach >= iNTracksMC))
    return kFALSE;
  // Check whether MC daughter particles have the same mother
  if(iIndexMotherNeg != iIndexMotherPos)
    return kFALSE;
  // Get the MC mother particle of both MC daughter particles (->V0)
  AliAODMCParticle* particleMCMotherV0 = (AliAODMCParticle*)MCPartArray->At(iIndexMotherPos);
  if(!particleMCMotherV0)
    return kFALSE;
  // Get the MC mother particle of bachelor daughter particle (->Cascade)
  AliAODMCParticle* particleMCMother = (AliAODMCParticle*)MCPartArray->At(iIndexMotherBach);
  if(!particleMCMother)
    return kFALSE;  
  // Get identity of the MC mother particles
  Int_t iPdgCodeMotherV0 = particleMCMotherV0->GetPdgCode();
  Int_t iPdgCodeMother = particleMCMother->GetPdgCode();
  // Skip not interesting particles
  if((iPdgCodeMotherV0 != iPdgCodeK0s) && (TMath::Abs(iPdgCodeMotherV0) != iPdgCodeLambda))
    return kFALSE;
  if((TMath::Abs(iPdgCodeMother) != iPdgCodeXi) && (TMath::Abs(iPdgCodeMother) != iPdgCodeOmega))
    return kFALSE;
  // Check identity of the MC mother particle and the decay channel
  // Is MC mother particle XiM?
  Bool_t bCascadeMCIsXiMinus = ((iPdgCodeMother == iPdgCodeXi) && (iPdgCodeMotherV0 == iPdgCodeLambda) && (iPdgCodeDaughterBach == -iPdgCodePion) && (iPdgCodeDaughterPos == +iPdgCodeProton) && (iPdgCodeDaughterNeg == -iPdgCodePion));
  // Is MC mother particle XiP?
  Bool_t bCascadeMCIsXiPlus = ((iPdgCodeMother == -iPdgCodeXi) && (iPdgCodeMotherV0 == -iPdgCodeLambda) && (iPdgCodeDaughterBach == iPdgCodePion) && (iPdgCodeDaughterPos == +iPdgCodePion) && (iPdgCodeDaughterNeg == -iPdgCodeProton));
  // Is MC mother particle OmegaM?
  Bool_t bCascadeMCIsOmegaMinus= ((iPdgCodeMother == iPdgCodeOmega) && (iPdgCodeMotherV0 == iPdgCodeLambda) && (iPdgCodeDaughterBach == -iPdgCodeKaon) && (iPdgCodeDaughterPos == +iPdgCodeProton) && (iPdgCodeDaughterNeg == -iPdgCodePion));
  // Is MC mother particle OmegaP?
  Bool_t bCascadeMCIsOmegaPlus= ((iPdgCodeMother == -iPdgCodeOmega)  && (iPdgCodeMotherV0 == -iPdgCodeLambda) && (iPdgCodeDaughterBach == iPdgCodeKaon) && (iPdgCodeDaughterPos == +iPdgCodePion) && (iPdgCodeDaughterNeg == -iPdgCodeProton));

  Double_t dPtCascadeGen = particleMCMother->Pt();
  Double_t dRapCascadeGen = particleMCMother->Y();
  Double_t dEtaCascadeGen = particleMCMother->Eta();
  
  // Cascade pseudorapidity cut applied on generated particles
  if(fdCutEtaCascadeMax > 0.) {
    if((TMath::Abs(dEtaCascadeGen) > fdCutEtaCascadeMax))
      return kFALSE;
  }
  // Cascade rapidity cut applied on generated particles
  if(fdCutRapCascadeMax > 0.) {
    if((TMath::Abs(dRapCascadeGen) > fdCutRapCascadeMax))
      return kFALSE;
  }
  // Select only particles from a specific generator
  if(!IsFromGoodGenerator(iIndexMotherBach))
    return kFALSE;
  // Get the MC mother particle of the MC mother particle

  // Get the distance between production point of the MC mother particle and the primary vertex
  Double_t dx = dPrimVtxMCX - particleMCMother->Xv();
  Double_t dy = dPrimVtxMCY - particleMCMother->Yv();
  Double_t dz = dPrimVtxMCZ - particleMCMother->Zv();
  Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
  Bool_t bCascadeMCIsPrimaryDist = (dDistPrimary < fdDistPrimaryMax); // Is close enough to be considered primary-like?

  if(bCascadeMCIsPrimaryDist && (bCascadeMCIsXiMinus || bCascadeMCIsXiPlus || bCascadeMCIsOmegaMinus || bCascadeMCIsOmegaPlus)) {
    fh1MCStats->Fill(0); //Reconstructed Cascades
    if(xjet) fh1MCStats->Fill(2); // Reconstrucetd Cascades in jets
  }  
  // XiMinus
  if(bIsXiMinus) { // selected candidates with any mass
    if(bCascadeMCIsXiMinus && bCascadeMCIsPrimaryDist) {// well reconstructed candidates
      if(xjet == 0) {
        fh2CascadeXiMinusPtMassMCRec[iCent]->Fill(dPtCascadeGen, dMXi);
        Double_t valueEtaK[3] = {dMXi, dPtCascadeGen, dEtaCascadeGen};
        fh3CascadeXiMinusEtaPtMassMCRec[iCent]->Fill(valueEtaK);
        Double_t valueEtaDKNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaCascadeGen, dPtCascadeGen, 0};
        fhnCascadeXiMinusInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKNeg);
        Double_t valueEtaDKPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaCascadeGen, dPtCascadeGen, 0};
        fhnCascadeXiMinusInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKPos);
        fh2CascadeXiMinusMCResolMPt[iCent]->Fill(dMXi - dMassPDGXi, dptcasc);
        fh2CascadeXiMinusMCPtGenPtRec[iCent]->Fill(dPtCascadeGen, dptcasc);
      }
      //here add an array with the well reconstructed particles? 
      else { // true Cascade associated to a candidate in jet
        Double_t valueKInJCMC[4] = {dMXi, dPtCascadeGen, dEtaCascadeGen, xjet->Pt()};
        fh3CascadeXiMinusInJetPtMassMCRec[iCent]->Fill(valueKInJCMC);
        Double_t valueEtaKIn[5] = {dMXi, dPtCascadeGen, dEtaCascadeGen, xjet->Pt(), dEtaCascadeGen - xjet->Eta()};
        fh4CascadeXiMinusInJetEtaPtMassMCRec[iCent]->Fill(valueEtaKIn);

        Double_t valueEtaDKJCNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaCascadeGen, dPtCascadeGen, xjet->Pt()};
        fhnCascadeXiMinusInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKJCNeg);
        Double_t valueEtaDKJCPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaCascadeGen, dPtCascadeGen, xjet->Pt()};
        fhnCascadeXiMinusInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKJCPos);
      }
    }
    if(bCascadeMCIsXiMinus && !bCascadeMCIsPrimaryDist) { // not primary K0s
      fh1CascadeXiMinusPtMCRecFalse[iCent]->Fill(dPtCascadeGen);
    }
  }
  // XiPlus
  if(bIsXiPlus) { // selected candidates with any mass
    if(bCascadeMCIsXiPlus && bCascadeMCIsPrimaryDist) {// well reconstructed candidates
      if(xjet == 0) {
        fh2CascadeXiPlusPtMassMCRec[iCent]->Fill(dPtCascadeGen, dMXi);
        Double_t valueEtaK[3] = {dMXi, dPtCascadeGen, dEtaCascadeGen};
        fh3CascadeXiPlusEtaPtMassMCRec[iCent]->Fill(valueEtaK);
        Double_t valueEtaDKNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaCascadeGen, dPtCascadeGen, 0};
        fhnCascadeXiPlusInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKNeg);
        Double_t valueEtaDKPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaCascadeGen, dPtCascadeGen, 0};
        fhnCascadeXiPlusInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKPos);
        fh2CascadeXiPlusMCResolMPt[iCent]->Fill(dMXi - dMassPDGXi, dptcasc);
        fh2CascadeXiPlusMCPtGenPtRec[iCent]->Fill(dPtCascadeGen, dptcasc);
      }
      //here add an array with the well reconstructed particles? 
      else { // true Cascade associated to a candidate in jet
        Double_t valueKInJCMC[4] = {dMXi, dPtCascadeGen, dEtaCascadeGen, xjet->Pt()};
        fh3CascadeXiPlusInJetPtMassMCRec[iCent]->Fill(valueKInJCMC);
        Double_t valueEtaKIn[5] = {dMXi, dPtCascadeGen, dEtaCascadeGen, xjet->Pt(), dEtaCascadeGen - xjet->Eta()};
        fh4CascadeXiPlusInJetEtaPtMassMCRec[iCent]->Fill(valueEtaKIn);

        Double_t valueEtaDKJCNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaCascadeGen, dPtCascadeGen, xjet->Pt()};
        fhnCascadeXiPlusInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKJCNeg);
        Double_t valueEtaDKJCPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaCascadeGen, dPtCascadeGen, xjet->Pt()};
        fhnCascadeXiPlusInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKJCPos);
      }
    }
    if(bCascadeMCIsXiPlus && !bCascadeMCIsPrimaryDist) { // not primary K0s
      fh1CascadeXiPlusPtMCRecFalse[iCent]->Fill(dPtCascadeGen);
    }
  }
  // OmegaMinus
  if(bIsOmegaMinus) { // selected candidates with any mass
    if(bCascadeMCIsOmegaMinus && bCascadeMCIsPrimaryDist) {// well reconstructed candidates
      if(xjet == 0) {
        fh2CascadeOmegaMinusPtMassMCRec[iCent]->Fill(dPtCascadeGen, dMOmega);
        Double_t valueEtaK[3] = {dMOmega, dPtCascadeGen, dEtaCascadeGen};
        fh3CascadeOmegaMinusEtaPtMassMCRec[iCent]->Fill(valueEtaK);
        Double_t valueEtaDKNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaCascadeGen, dPtCascadeGen, 0};
        fhnCascadeOmegaMinusInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKNeg);
        Double_t valueEtaDKPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaCascadeGen, dPtCascadeGen, 0};
        fhnCascadeOmegaMinusInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKPos);
        fh2CascadeOmegaMinusMCResolMPt[iCent]->Fill(dMOmega - dMassPDGOmega, dptcasc);
        fh2CascadeOmegaMinusMCPtGenPtRec[iCent]->Fill(dPtCascadeGen, dptcasc);
      }
      //here add an array with the well reconstructed particles? 
      else { // true Cascade associated to a candidate in jet
        Double_t valueKInJCMC[4] = {dMOmega, dPtCascadeGen, dEtaCascadeGen, xjet->Pt()};
        fh3CascadeOmegaMinusInJetPtMassMCRec[iCent]->Fill(valueKInJCMC);
        Double_t valueEtaKIn[5] = {dMOmega, dPtCascadeGen, dEtaCascadeGen, xjet->Pt(), dEtaCascadeGen - xjet->Eta()};
        fh4CascadeOmegaMinusInJetEtaPtMassMCRec[iCent]->Fill(valueEtaKIn);

        Double_t valueEtaDKJCNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaCascadeGen, dPtCascadeGen, xjet->Pt()};
        fhnCascadeOmegaMinusInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKJCNeg);
        Double_t valueEtaDKJCPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaCascadeGen, dPtCascadeGen, xjet->Pt()};
        fhnCascadeOmegaMinusInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKJCPos);
      }
    }
    if(bCascadeMCIsOmegaMinus && !bCascadeMCIsPrimaryDist) { // not primary K0s
      fh1CascadeOmegaMinusPtMCRecFalse[iCent]->Fill(dPtCascadeGen);
    }
  }
  // OmegaPlus
  if(bIsOmegaPlus) { // selected candidates with any mass
    if(bCascadeMCIsOmegaPlus && bCascadeMCIsPrimaryDist) {// well reconstructed candidates
      if(xjet == 0) {
        fh2CascadeOmegaPlusPtMassMCRec[iCent]->Fill(dPtCascadeGen, dMOmega);
        Double_t valueEtaK[3] = {dMOmega, dPtCascadeGen, dEtaCascadeGen};
        fh3CascadeOmegaPlusEtaPtMassMCRec[iCent]->Fill(valueEtaK);
        Double_t valueEtaDKNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaCascadeGen, dPtCascadeGen, 0};
        fhnCascadeOmegaPlusInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKNeg);
        Double_t valueEtaDKPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaCascadeGen, dPtCascadeGen, 0};
        fhnCascadeOmegaPlusInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKPos);
        fh2CascadeOmegaPlusMCResolMPt[iCent]->Fill(dMOmega - dMassPDGOmega, dptcasc);
        fh2CascadeOmegaPlusMCPtGenPtRec[iCent]->Fill(dPtCascadeGen, dptcasc);
      }
      //here add an array with the well reconstructed particles? 
      else { // true Cascade associated to a candidate in jet
        Double_t valueKInJCMC[4] = {dMOmega, dPtCascadeGen, dEtaCascadeGen, xjet->Pt()};
        fh3CascadeOmegaPlusInJetPtMassMCRec[iCent]->Fill(valueKInJCMC);
        Double_t valueEtaKIn[5] = {dMOmega, dPtCascadeGen, dEtaCascadeGen, xjet->Pt(), dEtaCascadeGen - xjet->Eta()};
        fh4CascadeOmegaPlusInJetEtaPtMassMCRec[iCent]->Fill(valueEtaKIn);

        Double_t valueEtaDKJCNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaCascadeGen, dPtCascadeGen, xjet->Pt()};
        fhnCascadeOmegaPlusInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKJCNeg);
        Double_t valueEtaDKJCPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaCascadeGen, dPtCascadeGen, xjet->Pt()};
        fhnCascadeOmegaPlusInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKJCPos);
      }
    }
    if(bCascadeMCIsOmegaPlus && !bCascadeMCIsPrimaryDist) { // not primary K0s
      fh1CascadeOmegaPlusPtMCRecFalse[iCent]->Fill(dPtCascadeGen);
    }
  }
  return kTRUE;
}

Bool_t AliAnalysisTaskCascadesInJets::GeneratedMCParticles(TClonesArray* track, Int_t iCent)
{
  	
  std::vector<fastjet::PseudoJet> InputBgParticlesMC; 
  //InputBgParticles.clear(); //clear before using again fot the bg estim for the generated mc particles 
  
  TClonesArray* MCPartArray = 0; // array particles in the MC event
  AliAODMCHeader* MCHeader = 0; // MC header
  Int_t iNTracksMC = 0; // number of MC tracks
  Double_t dPrimVtxMCX = 0., dPrimVtxMCY = 0., dPrimVtxMCZ = 0.; // position of the MC primary vertex
  
  // Simulation info
  MCPartArray = (TClonesArray*)fAODIn->FindListObject(AliAODMCParticle::StdBranchName());
  if(!MCPartArray) {
    printf("No MC array found!");
    return kFALSE;
  }
  if(fDebug > 2) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "MC array found");
  iNTracksMC = MCPartArray->GetEntriesFast();
  MCHeader = (AliAODMCHeader*)fAODIn->FindListObject(AliAODMCHeader::StdBranchName());
  if(!MCHeader) {
    printf("No MC header found!");
    return kFALSE;
  }
  // get position of the MC primary vertex
  dPrimVtxMCX = MCHeader->GetVtxX();
  dPrimVtxMCY = MCHeader->GetVtxY();
  dPrimVtxMCZ = MCHeader->GetVtxZ();

  Int_t iNMCCand = 0;

  static Int_t indxs[9999] = {-1};	

  for(Int_t iPartMC = 0; iPartMC < iNTracksMC; iPartMC++) {
    // Get MC particle
    AliAODMCParticle* particleMC = (AliAODMCParticle*)MCPartArray->At(iPartMC);
    if(!particleMC)
      continue;

    //need only final state particles ?
    //if(!particleMC->IsPhysicalPrimary())
      //continue;
    // Get identity of MC particle
    Int_t iPdgCodeParticleMC = particleMC->GetPdgCode();

    //  skip non Cascade particles 
    if((TMath::Abs(iPdgCodeParticleMC) != iPdgCodeXi) && (TMath::Abs(iPdgCodeParticleMC) != iPdgCodeOmega)) {
      continue;
    }
    // Check identity of the MC Cascade particle
    // Is MC Cascade particle IsXiMinus?
    Bool_t bCascadeMCIsXiMinus = (iPdgCodeParticleMC == -iPdgCodeXi);
    // Is MC Cascade particle IsXiPlus?
    Bool_t bCascadeMCIsXiPlus = (iPdgCodeParticleMC == +iPdgCodeXi);
    // Is MC Cascade particle OmegaMinus?
    Bool_t bCascadeMCIsOmegaMinus = (iPdgCodeParticleMC == -iPdgCodeOmega);
    // Is MC Cascade particle OmegaPlus?
    Bool_t bCascadeMCIsOmegaPlus = (iPdgCodeParticleMC == -iPdgCodeOmega);

    Double_t dPtCascadeGen = particleMC->Pt();
    Double_t dRapCascadeGen = particleMC->Y();
    Double_t dEtaCascadeGen = particleMC->Eta();

    // Cascade pseudorapidity cut
    if(fdCutEtaCascadeMax > 0.) {
      if((TMath::Abs(dEtaCascadeGen) > fdCutEtaCascadeMax)) 
        continue;
    }
    // Cascade rapidity cut
    if(fdCutRapCascadeMax > 0.) {
      if((TMath::Abs(dRapCascadeGen) > fdCutRapCascadeMax)) 
        continue;
    }   

    // Get the distance between the production point of the MC Cascade particle and the primary vertex
    Double_t dx = dPrimVtxMCX - particleMC->Xv();
    Double_t dy = dPrimVtxMCY - particleMC->Yv();
    Double_t dz = dPrimVtxMCZ - particleMC->Zv();
    Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
    Bool_t bCascadeMCIsPrimaryDist = (dDistPrimary < fdDistPrimaryMax); // Is close enough to be considered primary-like?
    // Select only primary-like MC Cascade particles
    if(!bCascadeMCIsPrimaryDist)
      continue;

    // Select only particles from a specific generator
    if(!IsFromGoodGenerator(iPartMC))
      continue;
 

    Int_t ind = 0;
    fh1MCStats->Fill(3); //Generated Cascades
    // XiMinus
    if(bCascadeMCIsXiMinus) {// well reconstructed candidates
      fh1CascadeXiMinusPtMCGen[iCent]->Fill(dPtCascadeGen);
      fh2CascadeXiMinusEtaPtMCGen[iCent]->Fill(dPtCascadeGen, dEtaCascadeGen);
      ind = iXiMinusId + iNMCCand;
    }
    // XiPlus
    if(bCascadeMCIsXiPlus) {// well reconstructed candidates
      fh1CascadeXiPlusPtMCGen[iCent]->Fill(dPtCascadeGen);
      fh2CascadeXiPlusEtaPtMCGen[iCent]->Fill(dPtCascadeGen, dEtaCascadeGen);
      ind = iXiPlusId + iNMCCand;
    }
    // OmegaMinus
    if(bCascadeMCIsOmegaMinus) {// well reconstructed candidates
      fh1CascadeOmegaMinusPtMCGen[iCent]->Fill(dPtCascadeGen);
      fh2CascadeOmegaMinusEtaPtMCGen[iCent]->Fill(dPtCascadeGen, dEtaCascadeGen);
      ind = iOmegaMinusId + iNMCCand;
    }
    // OmegaPlus
    if(bCascadeMCIsOmegaPlus) {// well reconstructed candidates
      fh1CascadeOmegaPlusPtMCGen[iCent]->Fill(dPtCascadeGen);
      fh2CascadeOmegaPlusEtaPtMCGen[iCent]->Fill(dPtCascadeGen, dEtaCascadeGen);
      ind = iOmegaPlusId + iNMCCand;
    }
    
    new ((*fGenMCCascade)[iNMCCand]) AliAODMCParticle(*particleMC); //  
    //add the MC cascade vector to the fastjetwrapper with modified id
    fFastJetWrapperMCGen.AddInputVector(particleMC->Px(), particleMC->Py(), particleMC->Pz(), particleMC->E(), ind);
    InputBgParticlesMC.push_back(fastjet::PseudoJet(particleMC->Px(), particleMC->Py(), particleMC->Pz(), particleMC->E()));
    
    iNMCCand++;
  }


  if(fGenMCCascade->GetEntriesFast() > 0) {
    AddEventTracksMC(fGenMCCascade, track, InputBgParticlesMC);
  }

  Double_t dAreaPercJetMin =  fdCutAreaPercJetMin*TMath::Pi()*fdRadius*fdRadius;

  //Background estimation 

  fastjet::JetMedianBackgroundEstimator bgeMC;
  fastjet::Subtractor subtr(&bgeMC);
  bgeMC.set_selector(selectorBG); 
  fastjet::JetDefinition jetDefBG(fastjet::kt_algorithm, fdBgRadius, fastjet::pt_scheme); //define the kT jet finding which will do the average background estimation
  fastjet::AreaDefinition areaDefBG(fastjet::active_area_explicit_ghosts);
  fastjet::ClusterSequenceArea cluster_seq_BG(InputBgParticlesMC, jetDefBG, areaDefBG);
  std::vector<fastjet::PseudoJet> jetsBGMC = sorted_by_pt(selectorBG(cluster_seq_BG.inclusive_jets())); //find the kT jets
  if (jetsBGMC.size() > 0) {
    bgeMC.set_jets(jetsBGMC);  // give the kT jets to the background estimator
    bgeMC.rho();
  }
  
  // run fjw
  fFastJetWrapperMCGen.Run(); 

  std::vector<fastjet::PseudoJet> vJetsMC = fFastJetWrapperMCGen.GetInclusiveJets(); 
  
 
  // sort jets according to jet pt
  GetSortedArray(indxs, vJetsMC); 
  AliDebug(1,Form("%i jets found", (Int_t)vJetsMC.size()));

  //Int_t iJetCount = 0;
  //Int_t iCascadeJets = 0;
  fastjet::PseudoJet jetSubMC;

  for (UInt_t ijet = 0; ijet < vJetsMC.size(); ++ijet) {
    Int_t ij = indxs[ijet];
    AliDebug(3,Form("Jet pt = %f, area = %f", vJetsMC[ij].perp(), fFastJetWrapperMCGen.GetJetArea(ij)));
    if (vJetsMC[ij].perp() < fdMinJetPt) continue;
    if (fFastJetWrapperMCGen.GetJetArea(ij) < fdMinJetArea) continue;
    if ((vJetsMC[ij].eta() < fdJetEtaMin) || (vJetsMC[ij].eta() > fdJetEtaMax) ||
        (vJetsMC[ij].phi() < fdJetPhiMin) || (vJetsMC[ij].phi() > fdJetPhiMax))
      continue;

    jetSubMC = subtr(vJetsMC[ij]);

    if(jetSubMC == 0) 
      continue; 
    
    if(jetSubMC.perp() < fdCutPtJetMin) // selection of high-pt jets, needs to be applied on the pt after bg subtraction
      continue;
    if(jetSubMC.area() < dAreaPercJetMin) //selection of the jets with area bigger than the cut (cut*pi*R2)
      continue;

    std::vector<fastjet::PseudoJet> constits(jetSubMC.constituents());
    Int_t iNConst = constits.size();
    Double_t dMaxTrPt = 0;
    for(Int_t ic = 0; ic < iNConst; ic++) {
		  Double_t dtrpt = constits[ic].perp();
      if(dtrpt > dMaxTrPt) dMaxTrPt = dtrpt;  
    }  
    cout << endl;
    if(dMaxTrPt < fdCutPtTrackJetMin)             // selection of jets with high leading track pt
      continue;                                           
    //printf(" JetPt: %f, sub jetPt: %f \n", vJetsMC[ij].perp(), jetSubMC.perp());
  
    fh1MCStats->Fill(4); //Generated jets
          
    Int_t uid   = -1;
    Int_t ik = 0;
    Int_t il = 0;
    Int_t ial = 0; 
    Int_t index = 0;
    AliAODMCParticle* jetv0 = 0;

    for(Int_t ic = 0; ic < iNConst; ic++){
	  uid = constits[ic].user_index();
      if(uid > iXiMinusId)
      	  fh1MCStats->Fill(5); //Generated Cascades in jets
      //Fill histograms for each particle: 	  
      if(uid > iXiMinusId && uid < iXiPlusId) { //XiMinus
		    index = uid-iXiMinusId; //will give the id of the Cascade particle 
		    jetv0 = (AliAODMCParticle*)fGenMCCascade->At(index); 
        fh2CascadeXiMinusInJetPtMCGen[iCent]->Fill(jetv0->Pt(), jetSubMC.perp());
        Double_t valueEtaKInGen[4] = {jetv0->Pt(), jetv0->Eta(), jetSubMC.perp(), jetv0->Eta() - jetSubMC.eta()};
        fh3CascadeXiMinusInJetEtaPtMCGen[iCent]->Fill(valueEtaKInGen);
		    ik++;
	    }	
      if(uid > iXiPlusId && uid < iOmegaPlusId) { //XiPlus
		    index = uid-iXiPlusId; //will give the id of the Cascade particle 
		    jetv0 = (AliAODMCParticle*)fGenMCCascade->At(index); 
        fh2CascadeXiPlusInJetPtMCGen[iCent]->Fill(jetv0->Pt(), jetSubMC.perp());
        Double_t valueEtaKInGen[4] = {jetv0->Pt(), jetv0->Eta(), jetSubMC.perp(), jetv0->Eta() - jetSubMC.eta()};
        fh3CascadeXiPlusInJetEtaPtMCGen[iCent]->Fill(valueEtaKInGen);
		    ik++;
	    }			  
      if(uid > iOmegaMinusId && uid < iOmegaPlusId) { //OmegaMinus
		    index = uid-iOmegaMinusId; //will give the id of the Cascade particle 
		    jetv0 = (AliAODMCParticle*)fGenMCCascade->At(index); 
        fh2CascadeOmegaMinusInJetPtMCGen[iCent]->Fill(jetv0->Pt(), jetSubMC.perp());
        Double_t valueEtaKInGen[4] = {jetv0->Pt(), jetv0->Eta(), jetSubMC.perp(), jetv0->Eta() - jetSubMC.eta()};
        fh3CascadeOmegaMinusInJetEtaPtMCGen[iCent]->Fill(valueEtaKInGen);
		    ik++;
	    }	
	    if(uid > iOmegaPlusId) { //OmegaPlus
		    index = uid-iOmegaPlusId; //will give the id of the Cascade particle 
		    jetv0 = (AliAODMCParticle*)fGenMCCascade->At(index); 
        fh2CascadeOmegaPlusInJetPtMCGen[iCent]->Fill(jetv0->Pt(), jetSubMC.perp());
        Double_t valueEtaKInGen[4] = {jetv0->Pt(), jetv0->Eta(), jetSubMC.perp(), jetv0->Eta() - jetSubMC.eta()};
        fh3CascadeOmegaPlusInJetEtaPtMCGen[iCent]->Fill(valueEtaKInGen);
		    ik++;
	    }	
    }
    //Int_t isum = ik + il + ial;
	  //printf("Tehere were %i K0s, %i Ls and %i ALs, isum: %i \n", ik, il, ial, isum);	  		    
  }

return kTRUE;
}


Double_t AliAnalysisTaskCascadesInJets::MassPeakSigma(Double_t pt, Int_t particle)
{
// estimation of the sigma of the invariant-mass peak as a function of pT and particle type
  switch(particle) {
    case 0: //Xi
      return 0.00127749 + 0.000339224 * pt; 
      break; 
    case 1: //Omega
      return 0.0037; //!!!!
      break;  
    default:
      return 0;
      break;
  }
}

Double_t AliAnalysisTaskCascadesInJets::MassPeakMean(Double_t pt, Int_t particle)
{
// estimation of the sigma of the invariant-mass peak as a function of pT and particle type
  switch(particle) {
    case 0: //Xi
      return 1.32228; 
      break; 
    case 1: //Omega
      return 1.67245; //!!!!
      break;    
    default:
      return 0;
      break;
  }
}