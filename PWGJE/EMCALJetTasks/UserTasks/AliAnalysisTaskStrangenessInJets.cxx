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
// Task for V0 analysis in charged jets with 
// the strange particles (instead of daughters) added to the  jet finder
// Author: Ekaterina Grecka (ermeeka@fjfi.cvut.cz)
// Modification of the AlianalysisTaskV0sInJetsEmcal task (author Vit Kucera) 
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

#include "AliAnalysisTaskStrangenessInJets.h"

ClassImp(AliAnalysisTaskStrangenessInJets)

// upper edges of centrality bins
const Int_t AliAnalysisTaskStrangenessInJets::fgkiCentBinRanges[] = {10}; // only central
// axis: pT of jets
const Double_t AliAnalysisTaskStrangenessInJets::fgkdBinsPtJet[] = {0, 100};
const Int_t AliAnalysisTaskStrangenessInJets::fgkiNBinsPtJet = sizeof(AliAnalysisTaskStrangenessInJets::fgkdBinsPtJet) / sizeof(AliAnalysisTaskStrangenessInJets::fgkdBinsPtJet[0]) - 1;
const Int_t AliAnalysisTaskStrangenessInJets::fgkiNBinsPtJetInit = int(((AliAnalysisTaskStrangenessInJets::fgkdBinsPtJet)[AliAnalysisTaskStrangenessInJets::fgkiNBinsPtJet] - (AliAnalysisTaskStrangenessInJets::fgkdBinsPtJet)[0]) / 10.); // bin width 10 GeV/c

// axis: pT of V0
const Double_t AliAnalysisTaskStrangenessInJets::fgkdBinsPtV0[] = {0, 15};
const Int_t AliAnalysisTaskStrangenessInJets::fgkiNBinsPtV0 = sizeof(AliAnalysisTaskStrangenessInJets::fgkdBinsPtV0) / sizeof((AliAnalysisTaskStrangenessInJets::fgkdBinsPtV0)[0]) - 1;
const Int_t AliAnalysisTaskStrangenessInJets::fgkiNBinsPtV0Init = int(((AliAnalysisTaskStrangenessInJets::fgkdBinsPtV0)[AliAnalysisTaskStrangenessInJets::fgkiNBinsPtV0] - (AliAnalysisTaskStrangenessInJets::fgkdBinsPtV0)[0]) / 0.1); // bin width 0.1 GeV/c
const Int_t AliAnalysisTaskStrangenessInJets::fgkiNBinsPtV0InitInJet = int(((AliAnalysisTaskStrangenessInJets::fgkdBinsPtV0)[AliAnalysisTaskStrangenessInJets::fgkiNBinsPtV0] - (AliAnalysisTaskStrangenessInJets::fgkdBinsPtV0)[0]) / 0.5); // bin width 0.5 GeV/c

// axis: K0S invariant mass
const Int_t AliAnalysisTaskStrangenessInJets::fgkiNBinsMassK0s = 300;
const Double_t AliAnalysisTaskStrangenessInJets::fgkdMassK0sMin = 0.35; // [GeV/c^2]
const Double_t AliAnalysisTaskStrangenessInJets::fgkdMassK0sMax = 0.65; // [GeV/c^2]
// axis: Lambda invariant mass
const Int_t AliAnalysisTaskStrangenessInJets::fgkiNBinsMassLambda = 200;
const Double_t AliAnalysisTaskStrangenessInJets::fgkdMassLambdaMin = 1.05; // [GeV/c^2]
const Double_t AliAnalysisTaskStrangenessInJets::fgkdMassLambdaMax = 1.25; // [GeV/c^2]

// PDG codes of used particles
const Int_t AliAnalysisTaskStrangenessInJets::iPdgCodePion = 211;
const Int_t AliAnalysisTaskStrangenessInJets::iPdgCodeProton = 2212;
const Int_t AliAnalysisTaskStrangenessInJets::iPdgCodeK0s = 310;
const Int_t AliAnalysisTaskStrangenessInJets::iPdgCodeLambda = 3122;

//Index codes to distinguish particles in the jet 
const Int_t AliAnalysisTaskStrangenessInJets::iK0Id = 1000000;
const Int_t AliAnalysisTaskStrangenessInJets::iLambdaId = 2000000;
const Int_t AliAnalysisTaskStrangenessInJets::iALambdaId = 3000000;
const Int_t AliAnalysisTaskStrangenessInJets::iK0LId = 4000000;
const Int_t AliAnalysisTaskStrangenessInJets::iK0ALId = 5000000;
const Int_t AliAnalysisTaskStrangenessInJets::iXiId = 6000000;
const Int_t AliAnalysisTaskStrangenessInJets::iAXiId = 7000000;
const Int_t AliAnalysisTaskStrangenessInJets::iXi0Id = 8000000;
const Int_t AliAnalysisTaskStrangenessInJets::iAXi0Id = 9000000;

// Default constructor
AliAnalysisTaskStrangenessInJets::AliAnalysisTaskStrangenessInJets():
  AliAnalysisTaskEmcal(),
  fOutputListStd(0),
  fOutputListStdJets(0),
  fOutputListMC(0),
  fV0CandidateArray(0),
  fGenMCV0(0),
  fGenMCXis(0),
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
  fdNSigmas(11), 
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
  fdCutV0PtMin(1.),
  fdCutCPAKMin(0),
  fdCutCPALMin(0),
  fdCutRadiusDecayMin(0),
  fdCutRadiusDecayMax(0),
  fdCutEtaV0Max(0),
  fdCutRapV0Max(0),
  fdCutNTauKMax(0),
  fdCutNTauLMax(0),
  fbCutArmPod(0),
  fbCutCross(0),  
    
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
  fdDistanceV0JetMax(0.4),
  fdDistPrimaryMax(0.01),
  
  fh1EventCounterCut(0),
  fh1EventCent(0),
  fh1EventCent2(0),
  fh1EventCent2Jets(0),
  fh1EventCent2NoJets(0),
  fh2EventCentTracks(0),
  fh2EventCentMult(0),
  fh1NV0sInJetStats(0),
  fh1NRndConeCent(0),
  fh1AreaExcluded(0),
  fh1MCStats(0)     
{	
  for(Int_t i = 0; i < fgkiNCategV0; i++) {
    fh1V0InvMassK0sAll[i] = 0;
    fh1V0InvMassLambdaAll[i] = 0;
    fh1V0InvMassALambdaAll[i] = 0;
  }	
  for(Int_t i = 0; i < fgkiNBinsCent; i++) {
    fh1EventCounterCutCent[i] = 0;
    fh1V0CounterCentK0s[i] = 0;
    fh1V0CounterCentLambda[i] = 0;
    fh1V0CounterCentALambda[i] = 0;
    fh1V0CandPerEventCentK0s[i] = 0;
    fh1V0CandPerEventCentLambda[i] = 0;
    fh1V0CandPerEventCentALambda[i] = 0;
    fh1V0InvMassK0sCent[i] = 0;
    fh1V0InvMassLambdaCent[i] = 0;
    fh1V0InvMassALambdaCent[i] = 0;
    fhnV0InclusiveK0s[i] = 0;
    fhnV0InclusiveLambda[i] = 0;
    fhnV0InclusiveALambda[i] = 0;   
    fhnV0InvMassCutK0s[i] = 0;
    fhnV0InvMassCutLambda[i] = 0;
    fhnV0InvMassCutALambda[i] = 0;   
    fh1VtxZ[i] = 0; 
    fh2VtxXY[i] = 0;   
    
    fh1PtJet[i] = 0;
    fh1EtaJet[i] = 0;
    fh2EtaPtJet[i] = 0;
    fh1PhiJet[i] = 0;
    fh2PtJetPtTrackLeading[i] = 0;
    fh1NJetPerEvent[i] = 0;   
    fh1NV0JetPerEvent[i] = 0;
    fh2EtaPhiRndCone[i] = 0;

    //In jets 
    fhnV0InJetK0s[i] = 0;
    fhnV0InJetLambda[i] = 0;
    fhnV0InJetALambda[i] = 0;    

    //In cones 
    fhnV0InPerpK0s[i] = 0;
    fhnV0InRndK0s[i] = 0;
    fhnV0InMedK0s[i] = 0;
    fhnV0OutJetK0s[i] = 0;
    fhnV0NoJetK0s[i] = 0;

    fhnV0InPerpLambda[i] = 0;
    fhnV0InRndLambda[i] = 0;
    fhnV0InMedLambda[i] = 0;
    fhnV0OutJetLambda[i] = 0;
    fhnV0NoJetLambda[i] = 0;

    fhnV0InPerpALambda[i] = 0;
    fhnV0InRndALambda[i] = 0;
    fhnV0InMedALambda[i] = 0;
    fhnV0OutJetALambda[i] = 0;
    fhnV0NoJetALambda[i] = 0;

    //MC 
    fh1V0K0sPtMCGen[i] = 0;
    fh2V0K0sPtMassMCRec[i] = 0;
    fh1V0K0sPtMCRecFalse[i] = 0;
    fh2V0K0sEtaPtMCGen[i] = 0;
    fh3V0K0sEtaPtMassMCRec[i] = 0;
    fh2V0K0sInJetPtMCGen[i] = 0;
    fh3V0K0sInJetPtMassMCRec[i] = 0;
    fh3V0K0sInJetEtaPtMCGen[i] = 0;
    fh4V0K0sInJetEtaPtMassMCRec[i] = 0;
    fh2V0K0sMCResolMPt[i] = 0;
    fh2V0K0sMCPtGenPtRec[i] = 0;
    
    fh1V0LambdaPtMCGen[i] = 0;
    fh2V0LambdaPtMassMCRec[i] = 0;
    fh1V0LambdaPtMCRecFalse[i] = 0;
    fh2V0LambdaEtaPtMCGen[i] = 0;
    fh3V0LambdaEtaPtMassMCRec[i] = 0;
    fh2V0LambdaInJetPtMCGen[i] = 0;
    fh3V0LambdaInJetPtMassMCRec[i] = 0;
    fh3V0LambdaInJetEtaPtMCGen[i] = 0;
    fh4V0LambdaInJetEtaPtMassMCRec[i] = 0;
    fh2V0LambdaMCResolMPt[i] = 0;
    fh2V0LambdaMCPtGenPtRec[i] = 0;
    fhnV0LambdaInclMCFromXi[i] = 0;
    fhnV0LambdaInclMCFromXi0[i] = 0;
    fhnV0LambdaInJetsMCFromXi[i] = 0;
    fhnV0LambdaInJetsMCFromXi0[i] = 0;
    fh1V0XiPtMCGen[i] = 0;
    fh1V0Xi0PtMCGen[i] = 0;
    fh1V0XiInJetPtMCGen[i] = 0;
    fh1V0Xi0InJetPtMCGen[i] = 0;

    fh1V0ALambdaPtMCGen[i] = 0;
    fh2V0ALambdaPtMassMCRec[i] = 0;
    fh1V0ALambdaPtMCRecFalse[i] = 0;
    fh2V0ALambdaEtaPtMCGen[i] = 0;
    fh3V0ALambdaEtaPtMassMCRec[i] = 0;
    fh2V0ALambdaInJetPtMCGen[i] = 0;
    fh3V0ALambdaInJetPtMassMCRec[i] = 0;
    fh3V0ALambdaInJetEtaPtMCGen[i] = 0;
    fh4V0ALambdaInJetEtaPtMassMCRec[i] = 0;
    fh2V0ALambdaMCResolMPt[i] = 0;
    fh2V0ALambdaMCPtGenPtRec[i] = 0;
    fhnV0ALambdaInclMCFromAXi[i] = 0;
    fhnV0ALambdaInclMCFromAXi0[i] = 0;
    fhnV0ALambdaInJetsMCFromAXi[i] = 0;
    fhnV0ALambdaInJetsMCFromAXi0[i] = 0;
    fh1V0AXiPtMCGen[i] = 0;
    fh1V0AXi0PtMCGen[i] = 0;
    fh1V0AXiInJetPtMCGen[i] = 0;
    fh1V0AXi0InJetPtMCGen[i] = 0;
    // eta daughters
    fhnV0K0sInclDaughterEtaPtPtMCRec[i] = 0;
    fhnV0K0sInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnV0LambdaInclDaughterEtaPtPtMCRec[i] = 0;
    fhnV0LambdaInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnV0ALambdaInclDaughterEtaPtPtMCRec[i] = 0;
    fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[i] = 0; 
  }  
}


// Constructor
AliAnalysisTaskStrangenessInJets::AliAnalysisTaskStrangenessInJets(const char* name):
  AliAnalysisTaskEmcal(name, kTRUE),
  fOutputListStd(0),
  fOutputListStdJets(0),
  fOutputListMC(0),
  fV0CandidateArray(0),
  fGenMCV0(0),
  fGenMCXis(0),                  //! contains MC generated V0s

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
  fdNSigmas(11),
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
  fdCutV0PtMin(1.),
  fdCutCPAKMin(0),
  fdCutCPALMin(0),
  fdCutRadiusDecayMin(0),
  fdCutRadiusDecayMax(0),
  fdCutEtaV0Max(0),
  fdCutRapV0Max(0),
  fdCutNTauKMax(0),
  fdCutNTauLMax(0),
  fbCutArmPod(0),
  fbCutCross(0),  
  
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
  fdDistanceV0JetMax(0.4),
  fdDistPrimaryMax(0.01),
  
  fh1EventCounterCut(0),
  fh1EventCent(0),
  fh1EventCent2(0),
  fh1EventCent2Jets(0),
  fh1EventCent2NoJets(0),
  fh2EventCentTracks(0),
  fh2EventCentMult(0),
  fh1NV0sInJetStats(0),
  fh1NRndConeCent(0),
  fh1AreaExcluded(0),
  fh1MCStats(0)    
{
  for(Int_t i = 0; i < fgkiNCategV0; i++) {
    fh1V0InvMassK0sAll[i] = 0;
    fh1V0InvMassLambdaAll[i] = 0;
    fh1V0InvMassALambdaAll[i] = 0;
  }	
  for(Int_t i = 0; i < fgkiNBinsCent; i++) {
    fh1EventCounterCutCent[i] = 0;
    fh1V0CounterCentK0s[i] = 0;
    fh1V0CounterCentLambda[i] = 0;
    fh1V0CounterCentALambda[i] = 0;
    fh1V0CandPerEventCentK0s[i] = 0;
    fh1V0CandPerEventCentLambda[i] = 0;
    fh1V0CandPerEventCentALambda[i] = 0;
    fh1V0InvMassK0sCent[i] = 0;
    fh1V0InvMassLambdaCent[i] = 0;
    fh1V0InvMassALambdaCent[i] = 0;
    fhnV0InclusiveK0s[i] = 0;
    fhnV0InclusiveLambda[i] = 0;
    fhnV0InclusiveALambda[i] = 0;   
    fhnV0InvMassCutK0s[i] = 0;
    fhnV0InvMassCutLambda[i] = 0;
    fhnV0InvMassCutALambda[i] = 0;   
    fh1VtxZ[i] = 0;
    fh2VtxXY[i] = 0;  
    
    fh1PtJet[i] = 0;
    fh1EtaJet[i] = 0;
    fh2EtaPtJet[i] = 0;
    fh1PhiJet[i] = 0;
    fh2PtJetPtTrackLeading[i] = 0;
    fh1NJetPerEvent[i] = 0;  
    fh1NV0JetPerEvent[i] = 0;
    fh2EtaPhiRndCone[i] = 0;

    //In jets 
    fhnV0InJetK0s[i] = 0;
    fhnV0InJetLambda[i] = 0;
    fhnV0InJetALambda[i] = 0;  
    
    //In cones 
    fhnV0InPerpK0s[i] = 0;
    fhnV0InRndK0s[i] = 0;
    fhnV0InMedK0s[i] = 0;
    fhnV0OutJetK0s[i] = 0;
    fhnV0NoJetK0s[i] = 0;

    fhnV0InPerpLambda[i] = 0;
    fhnV0InRndLambda[i] = 0;
    fhnV0InMedLambda[i] = 0;
    fhnV0OutJetLambda[i] = 0;
    fhnV0NoJetLambda[i] = 0;

    fhnV0InPerpALambda[i] = 0;
    fhnV0InRndALambda[i] = 0;
    fhnV0InMedALambda[i] = 0;
    fhnV0OutJetALambda[i] = 0;
    fhnV0NoJetALambda[i] = 0;

    //MC 
    fh1V0K0sPtMCGen[i] = 0;
    fh2V0K0sPtMassMCRec[i] = 0;
    fh1V0K0sPtMCRecFalse[i] = 0;
    fh2V0K0sEtaPtMCGen[i] = 0;
    fh3V0K0sEtaPtMassMCRec[i] = 0;
    fh2V0K0sInJetPtMCGen[i] = 0;
    fh3V0K0sInJetPtMassMCRec[i] = 0;
    fh3V0K0sInJetEtaPtMCGen[i] = 0;
    fh4V0K0sInJetEtaPtMassMCRec[i] = 0;
    fh2V0K0sMCResolMPt[i] = 0;
    fh2V0K0sMCPtGenPtRec[i] = 0;
    
    fh1V0LambdaPtMCGen[i] = 0;
    fh2V0LambdaPtMassMCRec[i] = 0;
    fh1V0LambdaPtMCRecFalse[i] = 0;
    fh2V0LambdaEtaPtMCGen[i] = 0;
    fh3V0LambdaEtaPtMassMCRec[i] = 0;
    fh2V0LambdaInJetPtMCGen[i] = 0;
    fh3V0LambdaInJetPtMassMCRec[i] = 0;
    fh3V0LambdaInJetEtaPtMCGen[i] = 0;
    fh4V0LambdaInJetEtaPtMassMCRec[i] = 0;
    fh2V0LambdaMCResolMPt[i] = 0;
    fh2V0LambdaMCPtGenPtRec[i] = 0;
    fhnV0LambdaInclMCFromXi[i] = 0;
    fhnV0LambdaInclMCFromXi0[i] = 0;
    fhnV0LambdaInJetsMCFromXi[i] = 0;
    fhnV0LambdaInJetsMCFromXi0[i] = 0;
    fh1V0XiPtMCGen[i] = 0;
    fh1V0Xi0PtMCGen[i] = 0;
    fh1V0XiInJetPtMCGen[i] = 0;
    fh1V0Xi0InJetPtMCGen[i] = 0;

    fh1V0ALambdaPtMCGen[i] = 0;
    fh2V0ALambdaPtMassMCRec[i] = 0;
    fh1V0ALambdaPtMCRecFalse[i] = 0;
    fh2V0ALambdaEtaPtMCGen[i] = 0;
    fh3V0ALambdaEtaPtMassMCRec[i] = 0;
    fh2V0ALambdaInJetPtMCGen[i] = 0;
    fh3V0ALambdaInJetPtMassMCRec[i] = 0;
    fh3V0ALambdaInJetEtaPtMCGen[i] = 0;
    fh4V0ALambdaInJetEtaPtMassMCRec[i] = 0;
    fh2V0ALambdaMCResolMPt[i] = 0;
    fh2V0ALambdaMCPtGenPtRec[i] = 0;
    fhnV0ALambdaInclMCFromAXi[i] = 0;
    fhnV0ALambdaInclMCFromAXi0[i] = 0;
    fhnV0ALambdaInJetsMCFromAXi[i] = 0;
    fhnV0ALambdaInJetsMCFromAXi0[i] = 0;
    fh1V0AXiPtMCGen[i] = 0;
    fh1V0AXi0PtMCGen[i] = 0;
    fh1V0AXiInJetPtMCGen[i] = 0;
    fh1V0AXi0InJetPtMCGen[i] = 0;
    // eta daughters
    fhnV0K0sInclDaughterEtaPtPtMCRec[i] = 0;
    fhnV0K0sInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnV0LambdaInclDaughterEtaPtPtMCRec[i] = 0;
    fhnV0LambdaInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnV0ALambdaInclDaughterEtaPtPtMCRec[i] = 0;
    fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[i] = 0; 
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
AliAnalysisTaskStrangenessInJets::~AliAnalysisTaskStrangenessInJets()
{
  delete fRandom;
  fRandom = 0;
  
  if (fV0CandidateArray) {
    delete fV0CandidateArray;
    fV0CandidateArray = 0;
  }
  if (fGenMCV0) {
    delete fGenMCV0;
    fGenMCV0 = 0;
  }
  if (fGenMCXis) {
    delete fGenMCXis;
    fGenMCXis = 0;
  }
  if (fJets) {
    delete fJets;
    fJets = 0;
  }
}

void AliAnalysisTaskStrangenessInJets::UserCreateOutputObjects()
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
  //List of V0 candidates 
  fV0CandidateArray = new TClonesArray("AliAODv0", 10); 
  fV0CandidateArray->SetOwner();

  //List of V0 candidates 
  fGenMCV0 = new TClonesArray("AliAODMCParticle", 10); 
  fGenMCV0->SetOwner();

  fGenMCXis= new TClonesArray("AliAODMCParticle", 10); 
  fGenMCXis->SetOwner();

  
  //List of all jets  
  fJets = new TClonesArray("AliEmcalJet");
  fJets->SetOwner(); 	
  
  // event categories
  const Int_t iNCategEvent = 14;
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
    "with V^{0}",         //9
    "with selected  V^{0}", //10
    "all jets",           //11 
    "jet selection",      //12
    "jet with V0"         //13
  };

  // labels for stages of V0 selection
  TString categV0[fgkiNCategV0] = {
    "all",                  //0
    "min V0 Pt",            //1
    "mass range",           //2
    "rec. method",          //3
    "tracks TPC",           //4
    "track pt",             //5
    "DCA prim v",           //6
    "DCA daughters",        //7
    "CPA",                  //8  
    "volume",               //9    
    "track #it{#eta}",      //10
    "V0 #it{y} & #it{#eta}",//11
    "lifetime",             //12 
    "PID",                  //13
    "Arm.-Pod.",            //14
    "cross-cont.",          //15
    "inclusive",            //16
    "in jet"                //17
  };  
  
  //labels for the in jet statistics
  const Int_t iNCategInJetStat = 12;
  TString categInJetStats[iNCategInJetStat] = {
    "N Jets",                       //0
    "V0s in jets",                  //1
    "K0s in jets",                  //2
    "Lambdas in jets",              //3
    "ALambdas in jets",             //4
    "K0s/Lambda in jets",           //5
    "K0s/ALambda in jets",          //6
    "More than 1 V0 in jet",        //7 
    "N jets before cuts(after bg sub)", //8
    "N jet after pt sel",           //9
    "N jet after area sel",         //10
    "N jet after lead track pt sel" //11
   }; 
 

  //histograms 
  fh1EventCounterCut = new TH1D("fh1EventCounterCut", "Number of events after filtering;selection filter;counts", iNCategEvent, 0, iNCategEvent);
  for(Int_t i = 0; i < iNCategEvent; i++)
    fh1EventCounterCut->GetXaxis()->SetBinLabel(i + 1, categEvent[i].Data());
  fh1EventCent2 = new TH1D("fh1EventCent2", "Number of events vs centrality;centrality;counts", 100, 0, 100);
  fh1EventCent2Jets = new TH1D("fh1EventCent2Jets", "Number of sel.-jet events vs centrality;centrality;counts", 100, 0, 100);
  fh1EventCent2NoJets = new TH1D("fh1EventCent2NoJets", "Number of no-jet events vs centrality;centrality;counts", 100, 0, 100);
  fh2EventCentTracks = new TH2D("fh2EventCentTracks", "Number of tracks vs centrality;centrality;tracks;counts", 100, 0, 100, 150, 0, 15e3);
  fh2EventCentMult = new TH2D("fh2EventCentMult", "Ref. multiplicity vs centrality;centrality;multiplicity;counts", 100, 0, 100, 150, 0, 15e3);
  fh1EventCent = new TH1D("fh1EventCent", "Number of events in centrality bins;centrality;counts", fgkiNBinsCent, 0, fgkiNBinsCent);
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
    fh1EventCent->GetXaxis()->SetBinLabel(i + 1, GetCentBinLabel(i).Data());
  fh1V0CandPerEvent = new TH1D("fh1V0CandPerEvent", "Number of all V0 candidates per event;candidates;events", 300, 0, 30000);
  fh1NRndConeCent = new TH1D("fh1NRndConeCent", "Number of rnd. cones in centrality bins;centrality;counts", fgkiNBinsCent, 0, fgkiNBinsCent);
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
    fh1NRndConeCent->GetXaxis()->SetBinLabel(i + 1, GetCentBinLabel(i).Data());
  fh1AreaExcluded = new TH1D("fh1AreaExcluded", "Area of excluded cones in centrality bins;centrality;area", fgkiNBinsCent, 0, fgkiNBinsCent);
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
  fh1AreaExcluded->GetXaxis()->SetBinLabel(i + 1, GetCentBinLabel(i).Data());
  
  fOutputListStd->Add(fh1EventCounterCut);
  fOutputListStd->Add(fh1EventCent);
  fOutputListStd->Add(fh1EventCent2);
  fOutputListStd->Add(fh1EventCent2Jets);
  fOutputListStd->Add(fh1EventCent2NoJets);
  fOutputListStd->Add(fh2EventCentTracks);
  fOutputListStd->Add(fh2EventCentMult);
  fOutputListStd->Add(fh1V0CandPerEvent);
  fOutputListStd->Add(fh1NRndConeCent);
  fOutputListStdJets->Add(fh1AreaExcluded); 

  fh1NV0sInJetStats = new TH1D("fh1NV0sInJetStats", "V0s in jets statistics", iNCategInJetStat, 0, iNCategInJetStat);   
  for(Int_t i = 0; i < iNCategInJetStat; i++)
  fh1NV0sInJetStats->GetXaxis()->SetBinLabel(i + 1, categInJetStats[i].Data());
  fOutputListStdJets->Add(fh1NV0sInJetStats);
  
  if(fbMCAnalysis) {
	//labels for theMC statistics
    const Int_t iNMCCategStat = 6;
    TString categMCStats[iNMCCategStat] = {
    "N rec V0s",             //0
    "N rec jets",            //1
    "N rec V0s in jets",     //2
    "N gen V0s",             //3
    "N gen jets",            //4
    "N gen V0s in jets",     //5
    };  
    fh1MCStats = new TH1D("fh1MCStats", "MC V0s in jets statistics", iNMCCategStat, 0, iNMCCategStat);   
    for(Int_t i = 0; i < iNMCCategStat; i++)
      fh1MCStats->GetXaxis()->SetBinLabel(i + 1, categMCStats[i].Data());
    fOutputListMC->Add(fh1MCStats);
  }
  
  // pt binning for V0 and jets
  Int_t iNBinsPtV0 = fgkiNBinsPtV0Init;
  Int_t iNBinsPtV0InJet = fgkiNBinsPtV0InitInJet;
  Double_t dPtV0Min = fgkdBinsPtV0[0];
  Double_t dPtV0Max = fgkdBinsPtV0[fgkiNBinsPtV0];
  Int_t iNJetPtBins = fgkiNBinsPtJetInit;
  Double_t dJetPtMin = fgkdBinsPtJet[0];
  Double_t dJetPtMax = fgkdBinsPtJet[fgkiNBinsPtJet];
  
  Double_t dStepEtaV0 = 0.05;
  Double_t dRangeEtaV0Max = 0.8;
  const Int_t iNBinsEtaV0 = 2 * Int_t(dRangeEtaV0Max / dStepEtaV0);
  //Binning parameters for invariant mass histograms
  const Int_t iNDimIncl = 3;
  Int_t binsKIncl[iNDimIncl] = {fgkiNBinsMassK0s, iNBinsPtV0, iNBinsEtaV0};
  Double_t xminKIncl[iNDimIncl] = {fgkdMassK0sMin, dPtV0Min, -dRangeEtaV0Max};
  Double_t xmaxKIncl[iNDimIncl] = {fgkdMassK0sMax, dPtV0Max, dRangeEtaV0Max};
  Int_t binsLIncl[iNDimIncl] = {fgkiNBinsMassLambda, iNBinsPtV0, iNBinsEtaV0};
  Double_t xminLIncl[iNDimIncl] = {fgkdMassLambdaMin, dPtV0Min, -dRangeEtaV0Max};
  Double_t xmaxLIncl[iNDimIncl] = {fgkdMassLambdaMax, dPtV0Max, dRangeEtaV0Max};
  //Binning in jets
  const Int_t iNDimInJC = 4;
  Int_t binsKInJC[iNDimInJC] = {fgkiNBinsMassK0s, iNBinsPtV0InJet, iNBinsEtaV0, iNJetPtBins};
  Double_t xminKInJC[iNDimInJC] = {fgkdMassK0sMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin};
  Double_t xmaxKInJC[iNDimInJC] = {fgkdMassK0sMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax};
  Int_t binsLInJC[iNDimInJC] = {fgkiNBinsMassLambda, iNBinsPtV0InJet, iNBinsEtaV0, iNJetPtBins};
  Double_t xminLInJC[iNDimInJC] = {fgkdMassLambdaMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin};
  Double_t xmaxLInJC[iNDimInJC] = {fgkdMassLambdaMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax};
 
  //Binning for MC
  // binning eff inclusive vs eta-pT
  Double_t dStepDeltaEta = 0.1;
  Double_t dRangeDeltaEtaMax = 0.5;
  const Int_t iNBinsDeltaEta = 2 * Int_t(dRangeDeltaEtaMax / dStepDeltaEta);
  Int_t binsEtaK[3] = {fgkiNBinsMassK0s, iNBinsPtV0, iNBinsEtaV0};
  Double_t xminEtaK[3] = {fgkdMassK0sMin, dPtV0Min, -dRangeEtaV0Max};
  Double_t xmaxEtaK[3] = {fgkdMassK0sMax, dPtV0Max, dRangeEtaV0Max};
  Int_t binsEtaL[3] = {fgkiNBinsMassLambda, iNBinsPtV0, iNBinsEtaV0};
  Double_t xminEtaL[3] = {fgkdMassLambdaMin, dPtV0Min, -dRangeEtaV0Max};
  Double_t xmaxEtaL[3] = {fgkdMassLambdaMax, dPtV0Max, dRangeEtaV0Max};
  // binning eff in jets vs eta-pT
  // associated
  Int_t binsEtaKInRec[5] = {fgkiNBinsMassK0s, iNBinsPtV0InJet, iNBinsEtaV0, iNJetPtBins, iNBinsDeltaEta};
  Double_t xminEtaKInRec[5] = {fgkdMassK0sMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin, -dRangeDeltaEtaMax};
  Double_t xmaxEtaKInRec[5] = {fgkdMassK0sMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax, dRangeDeltaEtaMax};
  Int_t binsEtaLInRec[5] = {fgkiNBinsMassLambda, iNBinsPtV0InJet, iNBinsEtaV0, iNJetPtBins, iNBinsDeltaEta};
  Double_t xminEtaLInRec[5] = {fgkdMassLambdaMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin, -dRangeDeltaEtaMax};
  Double_t xmaxEtaLInRec[5] = {fgkdMassLambdaMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax, dRangeDeltaEtaMax};
  // generated
  Int_t binsEtaInGen[4] = {iNBinsPtV0InJet, iNBinsEtaV0, iNJetPtBins, iNBinsDeltaEta};
  Double_t xminEtaInGen[4] = {dPtV0Min, -dRangeEtaV0Max, dJetPtMin, -dRangeDeltaEtaMax};
  Double_t xmaxEtaInGen[4] = {dPtV0Max, dRangeEtaV0Max, dJetPtMax, dRangeDeltaEtaMax};
  // daughter eta: charge-etaD-ptD-etaV0-ptV0-ptJet
  const Int_t iNDimEtaD = 6;
  Int_t binsEtaDaughter[iNDimEtaD] = {2, 20, iNBinsPtV0, iNBinsEtaV0, iNBinsPtV0, iNJetPtBins};
  Double_t xminEtaDaughter[iNDimEtaD] = {0, -1, dPtV0Min, -dRangeEtaV0Max, dPtV0Min, dJetPtMin};
  Double_t xmaxEtaDaughter[iNDimEtaD] = {2, 1, dPtV0Max, dRangeEtaV0Max, dPtV0Max, dJetPtMax};
  // daughter track pt: pos-neg-ptV0-ptJet-ptLead
  //const Int_t iNDimDaughter = 5;
  //Int_t binsDaughter[iNDimDaughter] = {80, 80, iNBinsPtV0, iNJetPtBins, 80};
  //Double_t xminDaughter[iNDimDaughter] = {0, 0, dPtV0Min, dJetPtMin, 0};
  //Double_t xmaxDaughter[iNDimDaughter] = {20, 20, dPtV0Max, dJetPtMax, 20};
 
  for(Int_t i = 0; i < fgkiNBinsCent; i++) {
    fh1EventCounterCutCent[i] = new TH1D(Form("fh1EventCounterCutCent_%d", i), Form("Number of events after filtering, cent %s;selection filter;counts", GetCentBinLabel(i).Data()), iNCategEvent, 0, iNCategEvent);
    for(Int_t j = 0; j < iNCategEvent; j++)
      fh1EventCounterCutCent[i]->GetXaxis()->SetBinLabel(j + 1, categEvent[j].Data());
    fh1V0CandPerEventCentK0s[i] = new TH1D(Form("fh1V0CandPerEventCentK0s_%d", i), Form("Number of selected K0s candidates per event, cent %s;candidates;events", GetCentBinLabel(i).Data()), 100, 0, 200);
    fh1V0CandPerEventCentLambda[i] = new TH1D(Form("fh1V0CandPerEventCentLambda_%d", i), Form("Number of selected Lambda candidates per event, cent %s;candidates;events", GetCentBinLabel(i).Data()), 100, 0, 200);
    fh1V0CandPerEventCentALambda[i] = new TH1D(Form("fh1V0CandPerEventCentALambda_%d", i), Form("Number of selected ALambda candidates per event, cent %s;candidates;events", GetCentBinLabel(i).Data()), 100, 0, 200);
    fh1V0CounterCentK0s[i] = new TH1D(Form("fh1V0CounterCentK0s_%d", i), Form("Number of K0s candidates after cuts, cent %s;cut;counts", GetCentBinLabel(i).Data()), fgkiNCategV0, 0, fgkiNCategV0);
    fh1V0CounterCentLambda[i] = new TH1D(Form("fh1V0CounterCentLambda_%d", i), Form("Number of Lambda candidates after cuts, cent %s;cut;counts", GetCentBinLabel(i).Data()), fgkiNCategV0, 0, fgkiNCategV0);
    fh1V0CounterCentALambda[i] = new TH1D(Form("fh1V0CounterCentALambda_%d", i), Form("Number of ALambda candidates after cuts, cent %s;cut;counts", GetCentBinLabel(i).Data()), fgkiNCategV0, 0, fgkiNCategV0);

    for(Int_t j = 0; j < fgkiNCategV0; j++) {	
      fh1V0CounterCentK0s[i]->GetXaxis()->SetBinLabel(j + 1, categV0[j].Data());
      fh1V0CounterCentLambda[i]->GetXaxis()->SetBinLabel(j + 1, categV0[j].Data());
      fh1V0CounterCentALambda[i]->GetXaxis()->SetBinLabel(j + 1, categV0[j].Data());
    }
    fOutputListStd->Add(fh1EventCounterCutCent[i]);
    fOutputListStd->Add(fh1V0CandPerEventCentK0s[i]);
    fOutputListStd->Add(fh1V0CandPerEventCentLambda[i]);
    fOutputListStd->Add(fh1V0CandPerEventCentALambda[i]);
    fOutputListStd->Add(fh1V0CounterCentK0s[i]);
    fOutputListStd->Add(fh1V0CounterCentLambda[i]);
    fOutputListStd->Add(fh1V0CounterCentALambda[i]);

    fh1V0InvMassK0sCent[i] = new TH1D(Form("fh1V0InvMassK0sCent_%d", i), Form("K0s: V0 invariant mass, cent %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", GetCentBinLabel(i).Data()), fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax);
    fh1V0InvMassLambdaCent[i] = new TH1D(Form("fh1V0InvMassLambdaCent_%d", i), Form("Lambda: V0 invariant mass, cent %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", GetCentBinLabel(i).Data()), fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
    fh1V0InvMassALambdaCent[i] = new TH1D(Form("fh1V0InvMassALambdaCent_%d", i), Form("ALambda: V0 invariant mass, cent %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", GetCentBinLabel(i).Data()), fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
    fOutputListStd->Add(fh1V0InvMassK0sCent[i]);
    fOutputListStd->Add(fh1V0InvMassLambdaCent[i]);
    fOutputListStd->Add(fh1V0InvMassALambdaCent[i]);
    // Inclusive
    fhnV0InclusiveK0s[i] = new THnSparseD(Form("fhnV0InclusiveK0s_C%d", i), "K0s: V0 invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};counts", iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fhnV0InclusiveLambda[i] = new THnSparseD(Form("fhnV0InclusiveLambda_C%d", i), "Lambda: V0 invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};counts", iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fhnV0InclusiveALambda[i] = new THnSparseD(Form("fhnV0InclusiveALambda_C%d", i), "ALambda: V0 invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};counts", iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0InclusiveK0s[i]);
    fOutputListStd->Add(fhnV0InclusiveLambda[i]);
    fOutputListStd->Add(fhnV0InclusiveALambda[i]);
    // After 3sigma invariant mass window cut
    fhnV0InvMassCutK0s[i] = new THnSparseD(Form("fhnV0InvMassCutK0s_C%d", i), "K0s after inv mass window cut: V0 invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};counts", iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fhnV0InvMassCutLambda[i] = new THnSparseD(Form("fhnV0InvMassCutLambda_C%d", i), "Lambda after inv mass window cut: V0 invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};counts", iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fhnV0InvMassCutALambda[i] = new THnSparseD(Form("fhnV0InvMassCutALambda_C%d", i), "ALambda after inv mass window cut: V0 invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};counts", iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0InvMassCutK0s[i]);
    fOutputListStd->Add(fhnV0InvMassCutLambda[i]);
    fOutputListStd->Add(fhnV0InvMassCutALambda[i]);
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
    fh1NV0JetPerEvent[i] = new TH1D(Form("fh1NV0JetPerEvent_%d", i), Form("Number of jets with V0 per event, cent: %s;# jets;# events", GetCentBinLabel(i).Data()), 100, 0., 100.);
    fh2EtaPhiRndCone[i] = new TH2D(Form("fh2EtaPhiRndCone_%d", i), Form("Rnd. cones: eta vs phi, cent: %s;#it{#eta} cone;#it{#phi} cone", GetCentBinLabel(i).Data()), 80, -1., 1., 90, 0., TMath::TwoPi());

    // V0s in jets
    fhnV0InJetK0s[i] = new THnSparseD(Form("fhnV0InJetK0s_%d", i), Form("K0s: Mass vs Pt in jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsKInJC, xminKInJC, xmaxKInJC);
    fhnV0InJetLambda[i] = new THnSparseD(Form("fhnV0InJetLambda_%d", i), Form("Lambda: Mass vs Pt in jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
    fhnV0InJetALambda[i] = new THnSparseD(Form("fhnV0InJetALambda_%d", i), Form("ALambda: Mass vs Pt in jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
            
    fOutputListStdJets->Add(fh1PtJet[i]);
    fOutputListStdJets->Add(fh1EtaJet[i]);
    fOutputListStdJets->Add(fh2EtaPtJet[i]);
    fOutputListStdJets->Add(fh1PhiJet[i]);
    fOutputListStdJets->Add(fh2PtJetPtTrackLeading[i]);
    fOutputListStdJets->Add(fh1NJetPerEvent[i]);
    fOutputListStdJets->Add(fh1NV0JetPerEvent[i]);
    fOutputListStdJets->Add(fhnV0InJetK0s[i]);
    fOutputListStdJets->Add(fhnV0InJetLambda[i]);
    fOutputListStdJets->Add(fhnV0InJetALambda[i]);
    fOutputListStdJets->Add(fh2EtaPhiRndCone[i]);


    //V0s in cones 
    fhnV0InPerpK0s[i] = new THnSparseD(Form("fhnV0InPerpK0s_%d", i), Form("K0s: Mass vs Pt in perp. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsKInJC, xminKInJC, xmaxKInJC);
    fhnV0InRndK0s[i] = new THnSparseD(Form("fhnV0InRndK0s_%d", i), Form("K0s: Mass vs Pt in rnd. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fhnV0InMedK0s[i] = new THnSparseD(Form("fhnV0InMedK0s_%d", i), Form("K0s: Mass vs Pt in med.-cl. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fhnV0OutJetK0s[i] = new THnSparseD(Form("fhnV0OutJetK0s_%d", i), Form("K0s: Pt outside jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fhnV0NoJetK0s[i] = new THnSparseD(Form("fhnV0NoJetK0s_%d", i), Form("K0s: Pt in jet-less events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    
    fOutputListStdJets->Add(fhnV0InPerpK0s[i]);
    fOutputListStdJets->Add(fhnV0InRndK0s[i]);
    fOutputListStdJets->Add(fhnV0InMedK0s[i]);
    fOutputListStdJets->Add(fhnV0OutJetK0s[i]);
    fOutputListStdJets->Add(fhnV0NoJetK0s[i]);

    fhnV0InPerpLambda[i] = new THnSparseD(Form("fhnV0InPerpLambda_%d", i), Form("Lambda: Mass vs Pt in perp. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
    fhnV0InRndLambda[i] = new THnSparseD(Form("fhnV0InRndLambda_%d", i), Form("Lambda: Mass vs Pt in rnd. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fhnV0InMedLambda[i] = new THnSparseD(Form("fhnV0InMedLambda_%d", i), Form("Lambda: Mass vs Pt in med.-cl. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fhnV0OutJetLambda[i] = new THnSparseD(Form("fhnV0OutJetLambda_%d", i), Form("Lambda: Pt outside jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fhnV0NoJetLambda[i] = new THnSparseD(Form("fhnV0NoJetLambda_%d", i), Form("Lambda: Pt in jet-less events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    
    fOutputListStdJets->Add(fhnV0InPerpLambda[i]);
    fOutputListStdJets->Add(fhnV0InRndLambda[i]);
    fOutputListStdJets->Add(fhnV0InMedLambda[i]);
    fOutputListStdJets->Add(fhnV0OutJetLambda[i]);
    fOutputListStdJets->Add(fhnV0NoJetLambda[i]);

    fhnV0InPerpALambda[i] = new THnSparseD(Form("fhnV0InPerpALambda_%d", i), Form("ALambda: Mass vs Pt in perp. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
    fhnV0InRndALambda[i] = new THnSparseD(Form("fhnV0InRndALambda_%d", i), Form("ALambda: Mass vs Pt in rnd. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fhnV0InMedALambda[i] = new THnSparseD(Form("fhnV0InMedALambda_%d", i), Form("ALambda: Mass vs Pt in med.-cl. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fhnV0OutJetALambda[i] = new THnSparseD(Form("fhnV0OutJetALambda_%d", i), Form("ALambda: Pt outside jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fhnV0NoJetALambda[i] = new THnSparseD(Form("fhnV0NoJetALambda_%d", i), Form("ALambda: Pt in jet-less events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    
    fOutputListStdJets->Add(fhnV0InPerpALambda[i]);
    fOutputListStdJets->Add(fhnV0InRndALambda[i]);
    fOutputListStdJets->Add(fhnV0InMedALambda[i]);
    fOutputListStdJets->Add(fhnV0OutJetALambda[i]);
    fOutputListStdJets->Add(fhnV0NoJetALambda[i]);

    if(fbMCAnalysis) {
      // inclusive pt
      fh1V0K0sPtMCGen[i] = new TH1D(Form("fh1V0K0sPtMCGen_%d", i), Form("MC K0s generated: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max);
      fh2V0K0sPtMassMCRec[i] = new TH2D(Form("fh2V0K0sPtMassMCRec_%d", i), Form("MC K0s associated: pt-m spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{m}_{inv} (GeV/#it{c}^{2})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax);
      fh1V0K0sPtMCRecFalse[i] = new TH1D(Form("fh1V0K0sPtMCRecFalse_%d", i), Form("MC K0s false: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh1V0K0sPtMCGen[i]);      
      fOutputListMC->Add(fh2V0K0sPtMassMCRec[i]);
      fOutputListMC->Add(fh1V0K0sPtMCRecFalse[i]);
      // inclusive pt-eta
      fh2V0K0sEtaPtMCGen[i] = new TH2D(Form("fh2V0K0sEtaPtMCGen_%d", i), Form("MC K0s generated: pt-eta spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNBinsEtaV0, -dRangeEtaV0Max, dRangeEtaV0Max);
      fh3V0K0sEtaPtMassMCRec[i] = new THnSparseD(Form("fh3V0K0sEtaPtMassMCRec_%d", i), Form("MC K0s associated: m-pt-eta spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), 3, binsEtaK, xminEtaK, xmaxEtaK);
      fOutputListMC->Add(fh2V0K0sEtaPtMCGen[i]);
      fOutputListMC->Add(fh3V0K0sEtaPtMassMCRec[i]);
      // in jet pt
      fh2V0K0sInJetPtMCGen[i] = new TH2D(Form("fh2V0K0sInJetPtMCGen_%d", i), Form("MC K0s in jet generated: pt-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0InJet, dPtV0Min, dPtV0Max, iNJetPtBins, dJetPtMin, dJetPtMax);
      fh3V0K0sInJetPtMassMCRec[i] = new THnSparseD(Form("fh3V0K0sInJetPtMassMCRec_%d", i), Form("MC K0s in jet associated: m-pt-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsKInJC, xminKInJC, xmaxKInJC); 
      fOutputListMC->Add(fh2V0K0sInJetPtMCGen[i]);      
      fOutputListMC->Add(fh3V0K0sInJetPtMassMCRec[i]);
      // in jet pt-eta
      fh3V0K0sInJetEtaPtMCGen[i] = new THnSparseD(Form("fh3V0K0sInJetEtaPtMCGen_%d", i), Form("MC K0s generated: pt-eta-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 4, binsEtaInGen, xminEtaInGen, xmaxEtaInGen);
      fh4V0K0sInJetEtaPtMassMCRec[i] = new THnSparseD(Form("fh4V0K0sInJetEtaPtMassMCRec_%d", i), Form("MC K0s associated: m-pt-eta-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 5, binsEtaKInRec, xminEtaKInRec, xmaxEtaKInRec);
      fh2V0K0sMCResolMPt[i] = new TH2D(Form("fh2V0K0sMCResolMPt_%d", i), Form("MC K0s associated: #Delta#it{m} vs pt, cent %s;#Delta#it{m} = #it{m}_{inv} - #it{m}_{true} (GeV/#it{c}^{2});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), 100, -0.02, 0.02, iNBinsPtV0, dPtV0Min, dPtV0Max);
      fh2V0K0sMCPtGenPtRec[i] = new TH2D(Form("fh2V0K0sMCPtGenPtRec_%d", i), Form("MC K0s associated: pt gen vs pt rec, cent %s;#it{p}_{T}^{gen} (GeV/#it{c});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh3V0K0sInJetEtaPtMCGen[i]);
      fOutputListMC->Add(fh4V0K0sInJetEtaPtMassMCRec[i]);
      fOutputListMC->Add(fh2V0K0sMCResolMPt[i]);
      fOutputListMC->Add(fh2V0K0sMCPtGenPtRec[i]);

      // inclusive pt
      fh1V0LambdaPtMCGen[i] = new TH1D(Form("fh1V0LambdaPtMCGen_%d", i), Form("MC Lambda generated: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh1V0LambdaPtMCGen[i]);
      fh2V0LambdaPtMassMCRec[i] = new TH2D(Form("fh2V0LambdaPtMassMCRec_%d", i), Form("MC Lambda associated: pt-m spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{m}_{inv} (GeV/#it{c}^{2})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
      fOutputListMC->Add(fh2V0LambdaPtMassMCRec[i]);
      fh1V0LambdaPtMCRecFalse[i] = new TH1D(Form("fh1V0LambdaPtMCRecFalse_%d", i), Form("MC Lambda false: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh1V0LambdaPtMCRecFalse[i]);
      // inclusive pt-eta
      fh2V0LambdaEtaPtMCGen[i] = new TH2D(Form("fh2V0LambdaEtaPtMCGen_%d", i), Form("MC Lambda generated: pt-eta spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNBinsEtaV0, -dRangeEtaV0Max, dRangeEtaV0Max);
      fOutputListMC->Add(fh2V0LambdaEtaPtMCGen[i]);
      fh3V0LambdaEtaPtMassMCRec[i] = new THnSparseD(Form("fh3V0LambdaEtaPtMassMCRec_%d", i), Form("MC Lambda associated: m-pt-eta spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), 3, binsEtaL, xminEtaL, xmaxEtaL);
      fOutputListMC->Add(fh3V0LambdaEtaPtMassMCRec[i]);
      // in jet pt
      fh2V0LambdaInJetPtMCGen[i] = new TH2D(Form("fh2V0LambdaInJetPtMCGen_%d", i), Form("MC Lambda in jet generated: pt-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0InJet, dPtV0Min, dPtV0Max, iNJetPtBins, dJetPtMin, dJetPtMax);
      fOutputListMC->Add(fh2V0LambdaInJetPtMCGen[i]);
      fh3V0LambdaInJetPtMassMCRec[i] = new THnSparseD(Form("fh3V0LambdaInJetPtMassMCRec_%d", i), Form("MC Lambda in jet associated: m-pt-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
      fOutputListMC->Add(fh3V0LambdaInJetPtMassMCRec[i]);
      // in jet pt-eta
      fh3V0LambdaInJetEtaPtMCGen[i] = new THnSparseD(Form("fh3V0LambdaInJetEtaPtMCGen_%d", i), Form("MC Lambda generated: pt-eta-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 4, binsEtaInGen, xminEtaInGen, xmaxEtaInGen);
      fOutputListMC->Add(fh3V0LambdaInJetEtaPtMCGen[i]);
      fh4V0LambdaInJetEtaPtMassMCRec[i] = new THnSparseD(Form("fh4V0LambdaInJetEtaPtMassMCRec_%d", i), Form("MC Lambda associated: m-pt-eta-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 5, binsEtaLInRec, xminEtaLInRec, xmaxEtaLInRec);
      fOutputListMC->Add(fh4V0LambdaInJetEtaPtMassMCRec[i]);

      fh2V0LambdaMCResolMPt[i] = new TH2D(Form("fh2V0LambdaMCResolMPt_%d", i), Form("MC Lambda associated: #Delta#it{m} vs pt, cent %s;#Delta#it{m} = #it{m}_{inv} - #it{m}_{true} (GeV/#it{c}^{2});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), 100, -0.02, 0.02, iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh2V0LambdaMCResolMPt[i]);
      fh2V0LambdaMCPtGenPtRec[i] = new TH2D(Form("fh2V0LambdaMCPtGenPtRec_%d", i), Form("MC Lambda associated: pt gen vs pt rec, cent %s;#it{p}_{T}^{gen} (GeV/#it{c});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh2V0LambdaMCPtGenPtRec[i]);

      // inclusive pt
      fh1V0ALambdaPtMCGen[i] = new TH1D(Form("fh1V0ALambdaPtMCGen_%d", i), Form("MC ALambda generated: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh1V0ALambdaPtMCGen[i]);
      fh2V0ALambdaPtMassMCRec[i] = new TH2D(Form("fh2V0ALambdaPtMassMCRec_%d", i), Form("MC ALambda associated: pt-m spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{m}_{inv} (GeV/#it{c}^{2})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
      fOutputListMC->Add(fh2V0ALambdaPtMassMCRec[i]);
      fh1V0ALambdaPtMCRecFalse[i] = new TH1D(Form("fh1V0ALambdaPtMCRecFalse_%d", i), Form("MC ALambda false: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh1V0ALambdaPtMCRecFalse[i]);
      // inclusive pt-eta
      fh2V0ALambdaEtaPtMCGen[i] = new TH2D(Form("fh2V0ALambdaEtaPtMCGen_%d", i), Form("MC ALambda generated: pt-eta spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNBinsEtaV0, -dRangeEtaV0Max, dRangeEtaV0Max);
      fOutputListMC->Add(fh2V0ALambdaEtaPtMCGen[i]);
      fh3V0ALambdaEtaPtMassMCRec[i] = new THnSparseD(Form("fh3V0ALambdaEtaPtMassMCRec_%d", i), Form("MC ALambda associated: m-pt-eta spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), 3, binsEtaL, xminEtaL, xmaxEtaL);
      fOutputListMC->Add(fh3V0ALambdaEtaPtMassMCRec[i]);
      // in jet pt
      fh2V0ALambdaInJetPtMCGen[i] = new TH2D(Form("fh2V0ALambdaInJetPtMCGen_%d", i), Form("MC ALambda in jet generated: pt-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0InJet, dPtV0Min, dPtV0Max, iNJetPtBins, dJetPtMin, dJetPtMax);
      fOutputListMC->Add(fh2V0ALambdaInJetPtMCGen[i]);
      fh3V0ALambdaInJetPtMassMCRec[i] = new THnSparseD(Form("fh3V0ALambdaInJetPtMassMCRec_%d", i), Form("MC ALambda in jet associated: m-pt-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
      fOutputListMC->Add(fh3V0ALambdaInJetPtMassMCRec[i]);
      // in jet pt-eta
      fh3V0ALambdaInJetEtaPtMCGen[i] = new THnSparseD(Form("fh3V0ALambdaInJetEtaPtMCGen_%d", i), Form("MC ALambda generated: pt-eta-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 4, binsEtaInGen, xminEtaInGen, xmaxEtaInGen);
      fOutputListMC->Add(fh3V0ALambdaInJetEtaPtMCGen[i]);
      fh4V0ALambdaInJetEtaPtMassMCRec[i] = new THnSparseD(Form("fh4V0ALambdaInJetEtaPtMassMCRec_%d", i), Form("MC ALambda associated: m-pt-eta-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 5, binsEtaLInRec, xminEtaLInRec, xmaxEtaLInRec);
      fOutputListMC->Add(fh4V0ALambdaInJetEtaPtMassMCRec[i]);

      fh2V0ALambdaMCResolMPt[i] = new TH2D(Form("fh2V0ALambdaMCResolMPt_%d", i), Form("MC ALambda associated: #Delta#it{m} vs pt, cent %s;#Delta#it{m} = #it{m}_{inv} - #it{m}_{true} (GeV/#it{c}^{2});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), 100, -0.02, 0.02, iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh2V0ALambdaMCResolMPt[i]);
      fh2V0ALambdaMCPtGenPtRec[i] = new TH2D(Form("fh2V0ALambdaMCPtGenPtRec_%d", i), Form("MC ALambda associated: pt gen vs pt rec, cent %s;#it{p}_{T}^{gen} (GeV/#it{c});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh2V0ALambdaMCPtGenPtRec[i]);

      // daughter eta
      fhnV0K0sInclDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0K0sInclDaughterEtaPtPtMCRec_%d", i), Form("MC K0S, inclusive, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen V0;pT gen V0;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnV0K0sInJetsDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0K0sInJetsDaughterEtaPtPtMCRec_%d", i), Form("MC K0S, in JC, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen V0;pT gen V0;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnV0LambdaInclDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0LambdaInclDaughterEtaPtPtMCRec_%d", i), Form("MC Lambda, inclusive, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen V0;pT gen V0;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnV0LambdaInJetsDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0LambdaInJetsDaughterEtaPtPtMCRec_%d", i), Form("MC Lambda, in JC, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen V0;pT gen V0;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnV0ALambdaInclDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0ALambdaInclDaughterEtaPtPtMCRec_%d", i), Form("MC ALambda, inclusive, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen V0;pT gen V0;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0ALambdaInJetsDaughterEtaPtPtMCRec_%d", i), Form("MC ALambda, in JC, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen V0;pT gen V0;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);

      fOutputListMC->Add(fhnV0K0sInclDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnV0K0sInJetsDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnV0LambdaInclDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnV0LambdaInJetsDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnV0ALambdaInclDaughterEtaPtPtMCRec[i]);

      Int_t iNBinsPtXi = 80;
      Double_t dPtXiMin = 0;
      Double_t dPtXiMax = 18;
      const Int_t iNDimFD = 3;
      Int_t binsFD[iNDimFD] = {iNBinsPtV0, iNBinsPtXi, iNJetPtBins};
      Double_t xminFD[iNDimFD] = {dPtV0Min, dPtXiMin, dJetPtMin};
      Double_t xmaxFD[iNDimFD] = {dPtV0Max, dPtXiMax, dJetPtMax};
      
      fhnV0LambdaInclMCFromXi[i] = new THnSparseD(Form("fhnV0LambdaInclMCFromXi_%d", i), Form("MC Lambda associated, inclusive, from Xi: pt-pt, cent %s;#it{p}_{T}^{#Lambda,rec.} (GeV/#it{c});#it{p}_{T}^{#Xi,rec.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0LambdaInclMCFromXi[i]);
      fhnV0LambdaInclMCFromXi0[i] = new THnSparseD(Form("fhnV0LambdaInclMCFromXi0_%d", i), Form("MC Lambda associated, inclusive, from Xi0: pt-pt, cent %s;#it{p}_{T}^{#Lambda,rec.} (GeV/#it{c});#it{p}_{T}^{#Xi0,rec.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0LambdaInclMCFromXi0[i]);

      fhnV0ALambdaInclMCFromAXi[i] = new THnSparseD(Form("fhnV0ALambdaInclMCFromAXi_%d", i), Form("MC ALambda associated, inclusive, from AXi: pt-pt, cent %s;#it{p}_{T}^{#Lambda,rec.} (GeV/#it{c});#it{p}_{T}^{#AXi,rec.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0ALambdaInclMCFromAXi[i]);
      fhnV0ALambdaInclMCFromAXi0[i] = new THnSparseD(Form("fhnV0ALambdaInclMCFromAXi0_%d", i), Form("MC ALambda associated, inclusive, from AXi0: pt-pt, cent %s;#it{p}_{T}^{#Lambda,rec.} (GeV/#it{c});#it{p}_{T}^{#AXi0,rec.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0ALambdaInclMCFromAXi0[i]);
      
      fhnV0LambdaInJetsMCFromXi[i] = new THnSparseD(Form("fhnV0LambdaInJetsMCFromXi_%d", i), Form("MC Lambda associated, in JC, from Xi: pt-pt-ptJet, cent %s;#it{p}_{T}^{#Lambda,rec.} (GeV/#it{c});#it{p}_{T}^{#Xi,rec.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0LambdaInJetsMCFromXi[i]);
      fhnV0LambdaInJetsMCFromXi0[i] = new THnSparseD(Form("fhnV0LambdaInJetsMCFromXi0_%d", i), Form("MC Lambda associated, in JC, from Xi0: pt-pt-ptJet, cent %s;#it{p}_{T}^{#Lambda,rec.} (GeV/#it{c});#it{p}_{T}^{#Xi0,rec.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0LambdaInJetsMCFromXi0[i]);

      fhnV0ALambdaInJetsMCFromAXi[i] = new THnSparseD(Form("fhnV0ALambdaInJetsMCFromAXi_%d", i), Form("MC ALambda associated, in JC, from AXi: pt-pt-ptJet, cent %s;#it{p}_{T}^{#Lambda,rec.} (GeV/#it{c});#it{p}_{T}^{#AXi,rec.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0ALambdaInJetsMCFromAXi[i]);
      fhnV0ALambdaInJetsMCFromAXi0[i] = new THnSparseD(Form("fhnV0ALambdaInJetsMCFromAXi0_%d", i), Form("MC ALambda associated, in JC, from AXi0: pt-pt-ptJet, cent %s;#it{p}_{T}^{#Lambda,rec.} (GeV/#it{c});#it{p}_{T}^{#AXi0,rec.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0ALambdaInJetsMCFromAXi0[i]);

      fh1V0XiPtMCGen[i] = new TH1D(Form("fh1V0XiPtMCGen_%d", i), Form("MC Xi^{-} generated: Pt spectrum, cent %s;#it{p}_{T}^{#Xi^{-},gen.} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtXi, dPtXiMin, dPtXiMax);
      fh1V0AXiPtMCGen[i] = new TH1D(Form("fh1V0AXiPtMCGen_%d", i), Form("MC AXi^{-} generated: Pt spectrum, cent %s;#it{p}_{T}^{A#Xi^{-},gen.} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtXi, dPtXiMin, dPtXiMax);
      fh1V0Xi0PtMCGen[i] = new TH1D(Form("fh1V0Xi0PtMCGen_%d", i), Form("MC Xi^{0} generated: Pt spectrum, cent %s;#it{p}_{T}^{#Xi^{0},gen.} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtXi, dPtXiMin, dPtXiMax);
      fh1V0AXi0PtMCGen[i] = new TH1D(Form("fh1V0AXi0PtMCGen_%d", i), Form("MC AXi^{0} generated: Pt spectrum, cent %s;#it{p}_{T}^{A#Xi^{0},gen.} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtXi, dPtXiMin, dPtXiMax);
 
      fh1V0XiInJetPtMCGen[i] = new TH1D(Form("fh1V0XiInJetPtMCGen_%d", i), Form("MC Xi^{-} generated in jet: Pt spectrum, cent %s;#it{p}_{T}^{#Xi^{-},gen.} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtXi, dPtXiMin, dPtXiMax);
      fh1V0AXiInJetPtMCGen[i] = new TH1D(Form("fh1V0AXiInJetPtMCGen_%d", i), Form("MC AXi^{-} generated in jet: Pt spectrum, cent %s;#it{p}_{T}^{A#Xi^{-},gen.} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtXi, dPtXiMin, dPtXiMax);
      fh1V0Xi0InJetPtMCGen[i] = new TH1D(Form("fh1V0Xi0InJetPtMCGen_%d", i), Form("MC Xi^{0} generated in jet: Pt spectrum, cent %s;#it{p}_{T}^{#Xi^{0},gen.} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtXi, dPtXiMin, dPtXiMax);
      fh1V0AXi0InJetPtMCGen[i] = new TH1D(Form("fh1V0AXi0InJetPtMCGen_%d", i), Form("MC AXi^{0} generated in jet: Pt spectrum, cent %s;#it{p}_{T}^{A#Xi^{0},gen.} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtXi, dPtXiMin, dPtXiMax);

      fOutputListMC->Add(fh1V0XiPtMCGen[i]);
      fOutputListMC->Add(fh1V0AXiPtMCGen[i]);     
      fOutputListMC->Add(fh1V0Xi0PtMCGen[i]);
      fOutputListMC->Add(fh1V0AXi0PtMCGen[i]); 
      fOutputListMC->Add(fh1V0XiInJetPtMCGen[i]);
      fOutputListMC->Add(fh1V0AXiInJetPtMCGen[i]);     
      fOutputListMC->Add(fh1V0Xi0InJetPtMCGen[i]);
      fOutputListMC->Add(fh1V0AXi0InJetPtMCGen[i]);   
    }
    
  }
  for(Int_t i = 0; i < fgkiNCategV0; i++) {
    fh1V0InvMassK0sAll[i] = new TH1D(Form("fh1V0InvMassK0sAll_%d", i), Form("K0s: V0 invariant mass, %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", categV0[i].Data()), fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax);
    fh1V0InvMassLambdaAll[i] = new TH1D(Form("fh1V0InvMassLambdaAll_%d", i), Form("Lambda: V0 invariant mass, %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", categV0[i].Data()), fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
    fh1V0InvMassALambdaAll[i] = new TH1D(Form("fh1V0InvMassALambdaAll_%d", i), Form("ALambda: V0 invariant mass, %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", categV0[i].Data()), fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
    fOutputListStd->Add(fh1V0InvMassK0sAll[i]);
    fOutputListStd->Add(fh1V0InvMassLambdaAll[i]);
    fOutputListStd->Add(fh1V0InvMassALambdaAll[i]);
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

void AliAnalysisTaskStrangenessInJets::ExecOnce()
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
  if(fdCutV0PtMin > 0.) printf("min V0 pT [Gev/c]: %g\n", fdCutV0PtMin);
  if(fdCutVertexZ > 0.) printf("max |z| of the prim vtx [cm]: %g\n", fdCutVertexZ);
  if(fdCutVertexR2 > 0.) printf("max r^2 of the prim vtx [cm^2]: %g\n", fdCutVertexR2);
  if(fiNContribMin > 0) printf("min number of prim vtx contributors: %d\n", fiNContribMin);
  if(fdCutDeltaZMax > 0.) printf("max |Delta z| between nominal prim vtx and SPD vtx [cm]: %g\n", fdCutDeltaZMax);
  if(fbUseIonutCut) printf("Ionut's cut\n");
  printf("-------------------------------------------------------\nV0 selection:\n");
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
  printf("V0 reconstruction method: %s\n", fbOnFly ? "on-the-fly" : "offline");
  if(fdCutCPAKMin > 0.) printf("min CPA, K0S: %g\n", fdCutCPAKMin);
  if(fdCutCPALMin > 0.) printf("min CPA, (A)Lambda: %g\n", fdCutCPALMin);
  if(fdCutRadiusDecayMin > 0. && fdCutRadiusDecayMax > 0.) printf("R of the decay vertex [cm]: %g-%g\n", fdCutRadiusDecayMin, fdCutRadiusDecayMax);
  if(fdCutEtaV0Max > 0.) printf("max |eta| of V0: %g\n", fdCutEtaV0Max);
  if(fdCutRapV0Max > 0.) printf("max |y| of V0: %g\n", fdCutRapV0Max);
  if(fdCutNTauKMax > 0.) printf("max proper lifetime, K0S [tau]: %g\n", fdCutNTauKMax);
  if(fdCutNTauLMax > 0.) printf("max proper lifetime, (A)Lambda [tau]: %g\n", fdCutNTauLMax);
  
  //if(fdCutNTauXMax > 0.) printf("max proper lifetime, Xi [tau]: %g\n", fdCutNTauXMax);
  
  if(fbCutArmPod) printf("Armenteros-Podolanski cut for K0S\n");
  
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

Bool_t AliAnalysisTaskStrangenessInJets::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  //Move the analysis part from FillHistograms here? 

  return kTRUE; // If return kFALSE FillHistogram() will NOT be executed.
}

Bool_t AliAnalysisTaskStrangenessInJets::FillHistograms()
{
  // Main loop, called for each event
  
  //open AOD 
  if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Start");

  //clear the TClonesArray from the previous event
  fV0CandidateArray->Clear();
  fGenMCV0->Clear();
  fGenMCXis->Clear();
  fJets->Clear();
  fFastJetWrapper.Clear();
  //fFastJetWrapperBG.Clear();
  fFastJetWrapperMCGen.Clear();  
  
  std::vector<fastjet::PseudoJet> InputBgParticles;

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

  //Get V0 information from AOD 
  Int_t iNV0s = fAODIn->GetNumberOfV0s(); // get the number of V0 candidates
  //printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("There are %d V0 candidates in the event", iNV0s));
  if(!iNV0s) {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "No V0s found in event");
  } else {
    fh1EventCounterCut->Fill(9); // events with V0s
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
	
  //Prepare cuts for the V0 analysis	
  AliAODv0* v0 = 0; // pointer to V0 candidates
  TVector3 vecV0Momentum; // 3D vector of V0 momentum
  
  Double_t dMassV0K0s = 0; // invariant mass of the K0s candidate
  Double_t dMassV0Lambda = 0; // invariant mass of the Lambda candidate
  Double_t dMassV0ALambda = 0; // invariant mass of the Lambda candidate 
  Int_t iNV0CandTot = 0; // counter of all V0 candidates at the beginning
  Int_t iNSelV0 = 0;  //Counter of all selected V0 candidates
  Int_t iNV0CandK0s = 0; // counter of K0s candidates at the end
  Int_t iNV0CandLambda = 0; // counter of Lambda candidates at the end
  Int_t iNV0CandALambda = 0; // counter of Lambda candidates at the end
  Bool_t bPrintCuts = 0; // print out which cuts are applied
  // Other cuts
  Double_t dNSigmaMassMax = fdNSigmas; // [sigma m] N of multiple of sigma for in jet bg estimation (used only for mass peak method of signal extraction)
  //Double_t dDistPrimaryMax = 0.01; // [cm] max distance of production point to the primary vertex (criterion for choice of MC particles considered as primary)
  // Mean lifetime
  Double_t dCTauK0s = 2.6844; // [cm] c*tau of K0S
  Double_t dCTauLambda = 7.89; // [cm] c*tau of Lambda
  // particle masses from PDG
  Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
  Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();

  fNCand = 0;
  Int_t iK0Id = 100000;
  Int_t iLambdaId = 200000;
  Int_t iALambdaId = 300000;
  Int_t iK0LId = 400000;
  Int_t iK0ALId = 500000;
  std::vector<Int_t> ivecV0CandIndex;
  // Loading primary vertex info
  AliAODVertex* primVtx = fAODIn->GetPrimaryVertex(); // get the primary vertex
  Double_t dPrimVtxPos[3]; // primary vertex position {x,y,z}
  primVtx->GetXYZ(dPrimVtxPos);
  fh1VtxZ[iCentIndex]->Fill(dPrimVtxPos[2]);
  fh2VtxXY[iCentIndex]->Fill(dPrimVtxPos[0], dPrimVtxPos[1]);
  
  Double_t dCutEtaJetMax = fdCutEtaV0Max - fdDistanceV0JetMax; // max jet |pseudorapidity|, to make sure that V0s can appear in the entire jet area
  Double_t dRadiusExcludeCone = 2 * fdDistanceV0JetMax; // radius of cones around jets excluded for V0 outside jets

  Double_t dAreaPercJetMin =  fdCutAreaPercJetMin*TMath::Pi()*fdRadius*fdRadius;

  //===== Start of loop over V0 candidates =====
  if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Start of V0 loop");
  for(Int_t iV0 = 0; iV0 < iNV0s; iV0++) {
    v0 = fAODIn->GetV0(iV0); // get next candidate from the list in AOD
    if(!v0)
      continue;

    iNV0CandTot++;	
    
    Int_t iCutIndex = 0; // indicator of current selection step    
    // Initialization of status indicators
    Bool_t bIsCandidateK0s = kTRUE; // candidate for K0s
    Bool_t bIsCandidateLambda = kTRUE; // candidate for Lambda
    Bool_t bIsCandidateALambda = kTRUE; // candidate for anti-Lambda
    Bool_t bIsInPeakK0s = kFALSE; // candidate within the K0s mass peak
    Bool_t bIsInPeakLambda = kFALSE; // candidate within the Lambda mass peak
    Bool_t bIsInPeakALambda = kFALSE; // candidate within the anti-Lambda mass peak
    // Invariant mass calculation
    dMassV0K0s = v0->MassK0Short();
    dMassV0Lambda = v0->MassLambda();
    dMassV0ALambda = v0->MassAntiLambda();

    // 0
    // All V0 candidates
    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    Double_t dPtV0 = TMath::Sqrt(v0->Pt2V0()); // transverse momentum of V0
    vecV0Momentum = TVector3(v0->Px(), v0->Py(), v0->Pz()); // set the vector of V0 momentum
    
    //1
    //V0 min pt 
    if(dPtV0 < fdCutV0PtMin)
      continue;
    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // Sigma of the mass peak window
    Double_t dMassPeakWindowK0s = dNSigmaMassMax * MassPeakSigma(dPtV0, 0);
    Double_t dMassPeakWindowLambda = dNSigmaMassMax * MassPeakSigma(dPtV0, 1);
    //Mean of the mass peak window
    Double_t dMassPeakWindowMeanK0s = MassPeakMean(dPtV0, 0);
    Double_t dMassPeakWindowMeanLambda = MassPeakMean(dPtV0, 1);
    if(!fbIsPbPb) { // p-p
      dMassPeakWindowK0s = 0.010; // LF p-p
      dMassPeakWindowLambda = 0.005; // LF p-p
    }

    // Invariant mass peak selection
    if( (dMassV0K0s > dMassPeakWindowMeanK0s - dMassPeakWindowK0s) && (dMassV0K0s < dMassPeakWindowMeanK0s +  dMassPeakWindowK0s) )    //Signal+BG (-11sigma,11sigma)
      bIsInPeakK0s = kTRUE;
    if( (dMassV0Lambda > dMassPeakWindowMeanLambda - dMassPeakWindowLambda) && (dMassV0Lambda < dMassPeakWindowMeanLambda + dMassPeakWindowLambda))
      bIsInPeakLambda = kTRUE;
    if( (dMassV0ALambda > dMassPeakWindowMeanLambda - dMassPeakWindowLambda) && (dMassV0ALambda < dMassPeakWindowMeanLambda + dMassPeakWindowLambda))
      bIsInPeakALambda = kTRUE;
    /*if(fbSignalInBG == 0) { //Inv mass in signal region
      if(TMath::Abs(dMassV0K0s - dMassPeakWindowMeanK0s ) < dMassPeakWindowK0s)
        bIsInPeakK0s = kTRUE;
      if(TMath::Abs(dMassV0Lambda - dMassPeakWindowMeanLambda) < dMassPeakWindowLambda)
        bIsInPeakLambda = kTRUE;
      if(TMath::Abs(dMassV0ALambda - dMassPeakWindowMeanLambda) < dMassPeakWindowLambda)
        bIsInPeakALambda = kTRUE;
    }  
    else if(fbSignalInBG == 1){ //Inv mass in BG region
      if( (dMassV0K0s > dMassPeakWindowMeanK0s + dMassPeakWindowK0s) && (dMassV0K0s < dMassPeakWindowMeanK0s + 2 * dMassPeakWindowK0s) )
        bIsInPeakK0s = kTRUE;
      if( (dMassV0Lambda > dMassPeakWindowMeanLambda + dMassPeakWindowLambda) && (dMassV0Lambda < dMassPeakWindowMeanLambda + 2 * dMassPeakWindowLambda))
        bIsInPeakLambda = kTRUE;
      if( (dMassV0ALambda > dMassPeakWindowMeanLambda + dMassPeakWindowLambda) && (dMassV0ALambda < dMassPeakWindowMeanLambda + 2 * dMassPeakWindowLambda))
        bIsInPeakALambda = kTRUE;
    }*/

    // Skip candidates outside the histogram range
    if((dMassV0K0s < fgkdMassK0sMin) || (dMassV0K0s >= fgkdMassK0sMax))
      bIsCandidateK0s = kFALSE;
    if((dMassV0Lambda < fgkdMassLambdaMin) || (dMassV0Lambda >= fgkdMassLambdaMax))
      bIsCandidateLambda = kFALSE;
    if((dMassV0ALambda < fgkdMassLambdaMin) || (dMassV0ALambda >= fgkdMassLambdaMax))
      bIsCandidateALambda = kFALSE;
    if(!bIsCandidateK0s && !bIsCandidateLambda && !bIsCandidateALambda)
      continue;

    // Retrieving all relevant properties of the V0 candidate
    Bool_t bOnFlyStatus = v0->GetOnFlyStatus(); // online (on fly) reconstructed vs offline reconstructed
    const AliAODTrack* trackPos = (AliAODTrack*)v0->GetDaughter(0); // positive daughter track
    const AliAODTrack* trackNeg = (AliAODTrack*)v0->GetDaughter(1); // negative daughter track
    Double_t dPtDaughterPos = trackPos->Pt(); // transverse momentum of a daughter track calculated as if primary, != v0->PtProng(0)
    Double_t dPtDaughterNeg = trackNeg->Pt(); // != v0->PtProng(1)
    Double_t dNRowsPos = trackPos->GetTPCClusterInfo(2, 1); // crossed TPC pad rows of a daughter track
    Double_t dNRowsNeg = trackNeg->GetTPCClusterInfo(2, 1);
    Double_t dFindablePos = Double_t(trackPos->GetTPCNclsF()); // Findable clusters
    Double_t dFindableNeg = Double_t(trackNeg->GetTPCNclsF());
    Double_t dDCAToPrimVtxPos = TMath::Abs(v0->DcaPosToPrimVertex()); // dca of a daughter to the primary vertex
    Double_t dDCAToPrimVtxNeg = TMath::Abs(v0->DcaNegToPrimVertex());
    Double_t dDCADaughters = v0->DcaV0Daughters(); // dca between daughters
    Double_t dCPA = v0->CosPointingAngle(primVtx); // cosine of the pointing angle
      
    Double_t dSecVtxPos[3]; // V0 vertex position {x,y,z}
    // Double_t dSecVtxPos[3] = {v0->DecayVertexV0X(),v0->DecayVertexV0Y(),v0->DecayVertexV0Z()}; // V0 vertex position
    v0->GetSecondaryVtx(dSecVtxPos);
    Double_t dRadiusDecay = TMath::Sqrt(dSecVtxPos[0] * dSecVtxPos[0] + dSecVtxPos[1] * dSecVtxPos[1]); // distance of the V0 vertex from the z-axis
    Double_t dEtaDaughterPos = trackPos->Eta(); // pseudorapidity of a daughter track calculated as if primary, != v0->EtaProng(0)
    Double_t dEtaDaughterNeg = trackNeg->Eta(); // != v0->EtaProng(1);
    Double_t dRapK0s = v0->RapK0Short(); // rapidity calculated for K0s assumption
    Double_t dRapLambda = v0->RapLambda(); // rapidity calculated for Lambda assumption
    Double_t dEtaV0 = v0->Eta(); // V0 pseudorapidity
    //Double_t dPhiV0 = v0->Phi(); // V0 azimuth
    Double_t dDecayPath[3];
    for(Int_t iPos = 0; iPos < 3; iPos++)
      dDecayPath[iPos] = dSecVtxPos[iPos] - dPrimVtxPos[iPos]; // vector of the V0 path
    //Double_t dDecLen = TMath::Sqrt(dDecayPath[0] * dDecayPath[0] + dDecayPath[1] * dDecayPath[1] + dDecayPath[2] * dDecayPath[2]); // path length L
    Double_t dDecLen2D = TMath::Sqrt(dDecayPath[0] * dDecayPath[0] + dDecayPath[1] * dDecayPath[1]); // transverse path length R
    //Double_t dLOverP = dDecLen / v0->P(); // L/p
    Double_t dROverPt = dDecLen2D / dPtV0; // R/pT
    //Double_t dMLOverPK0s = dMassPDGK0s * dLOverP; // m*L/p = c*(proper lifetime)
    // Double_t dMLOverPLambda = dMassPDGLambda*dLOverP; // m*L/p
    Double_t dMROverPtK0s = dMassPDGK0s * dROverPt; // m*R/pT
    Double_t dMROverPtLambda = dMassPDGLambda * dROverPt; // m*R/pT
    Double_t dNSigmaPosPion   = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPos, AliPID::kPion)) : 0.); // difference between measured and expected signal of the dE/dx in the TPC
    Double_t dNSigmaPosProton = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPos, AliPID::kProton)) : 0.);
    Double_t dNSigmaNegPion   = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNeg, AliPID::kPion)) : 0.);
    Double_t dNSigmaNegProton = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNeg, AliPID::kProton)) : 0.);
    Double_t dAlpha = v0->AlphaV0(); // Armenteros-Podolanski alpha
    Double_t dPtArm = v0->PtArmV0(); // Armenteros-Podolanski pT
    AliAODVertex* prodVtxDaughterPos = (AliAODVertex*)(trackPos->GetProdVertex()); // production vertex of the positive daughter track
    Char_t cTypeVtxProdPos = prodVtxDaughterPos->GetType(); // type of the production vertex
    AliAODVertex* prodVtxDaughterNeg = (AliAODVertex*)(trackNeg->GetProdVertex()); // production vertex of the negative daughter track
    Char_t cTypeVtxProdNeg = prodVtxDaughterNeg->GetType(); // type of the production vertex
    
    Double_t dEnergy = 0;
   //===== Start of reconstruction cutting =====

    // 2
    // All V0 candidates in mass range
    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // Start of global cuts
    // 3
    // Reconstruction method
    if(bPrintCuts) printf("Rec: Applying cut: Reconstruction method: %s\n", (fbOnFly ? "on-the-fly" : "offline"));
    if(bOnFlyStatus != fbOnFly)
      continue;
    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 4
    // Tracks TPC OK
    if(bPrintCuts) printf("Rec: Applying cut: Correct charge of daughters\n");
    if(!trackNeg || !trackPos)
      continue;
    if(trackNeg->Charge() == trackPos->Charge()) // daughters have different charge?
      continue;
    if(trackNeg->Charge() != -1) // daughters have expected charge?
      continue;
    if(trackPos->Charge() != 1) // daughters have expected charge?
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

    if(fbFindableClusters) {
      if(bPrintCuts) printf("Rec: Applying cut: Positive number of findable clusters\n");
      if(dFindableNeg <= 0.)
        continue;
      if(dFindablePos <= 0.)
        continue;
    }

    if(fdCutNCrossedRowsTPCMin > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: Number of TPC rows >= %g\n", fdCutNCrossedRowsTPCMin);
      if(dNRowsNeg < fdCutNCrossedRowsTPCMin) // Crossed TPC padrows
        continue;
      if(dNRowsPos < fdCutNCrossedRowsTPCMin)
        continue;
    }

    if(fdCutCrossedRowsOverFindMin > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: rows/findable >= %g\n", fdCutCrossedRowsOverFindMin);
      if(dNRowsNeg / dFindableNeg < fdCutCrossedRowsOverFindMin)
        continue;
      if(dNRowsPos / dFindablePos < fdCutCrossedRowsOverFindMin)
        continue;
    }

    if(fdCutCrossedRowsOverFindMax > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: rows/findable <= %g\n", fdCutCrossedRowsOverFindMax);
      if(dNRowsNeg / dFindableNeg > fdCutCrossedRowsOverFindMax)
        continue;
      if(dNRowsPos / dFindablePos > fdCutCrossedRowsOverFindMax)
        continue;
    }

    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 5
    // Daughters: transverse momentum cut
    if(fdCutPtDaughterMin > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: Daughter pt >= %g\n", fdCutPtDaughterMin);
      if((dPtDaughterNeg < fdCutPtDaughterMin) || (dPtDaughterPos < fdCutPtDaughterMin))
        continue;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 6
    // Daughters: Impact parameter of daughters to prim vtx
    if(fdCutDCAToPrimVtxMin > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: Daughter DCA to prim vtx >= %g\n", fdCutDCAToPrimVtxMin);
      if((dDCAToPrimVtxNeg < fdCutDCAToPrimVtxMin) || (dDCAToPrimVtxPos < fdCutDCAToPrimVtxMin))
        continue;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 7
    // Daughters: DCA
    if(fdCutDCADaughtersMax > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: DCA between daughters <= %g\n", fdCutDCADaughtersMax);
      if(dDCADaughters > fdCutDCADaughtersMax)
        continue;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 8
    // V0: Cosine of the pointing angle
    if(fdCutCPAKMin > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: CPA >= %g (K)\n", fdCutCPAKMin);
      if(dCPA < fdCutCPAKMin)
        bIsCandidateK0s = kFALSE;
    }
    if(fdCutCPALMin > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: CPA >= %g (L, AL)\n", fdCutCPALMin);
      if(dCPA < fdCutCPALMin)
      {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
    }
    if(fdCutCPAKMin > 0. || fdCutCPALMin > 0.)
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 9
    // V0: Fiducial volume
    if(fdCutRadiusDecayMin > 0. && fdCutRadiusDecayMax > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: Decay radius >= %g, <= %g\n", fdCutRadiusDecayMin, fdCutRadiusDecayMax);
      if((dRadiusDecay < fdCutRadiusDecayMin) || (dRadiusDecay > fdCutRadiusDecayMax))
        continue;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 10
    // Daughters: pseudorapidity cut
    if(fdCutEtaDaughterMax > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: Daughter |eta| < %g\n", fdCutEtaDaughterMax);
      if((TMath::Abs(dEtaDaughterNeg) > fdCutEtaDaughterMax) || (TMath::Abs(dEtaDaughterPos) > fdCutEtaDaughterMax))
        continue;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;
    // End of global cuts

    // Start of particle-dependent cuts
    // 11
    // V0: pseudorapidity cut & rapidity cut
    if(fdCutEtaV0Max > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: V0 |eta| < %g\n", fdCutEtaV0Max);
      if(TMath::Abs(dEtaV0) > fdCutEtaV0Max) {
        bIsCandidateK0s = kFALSE;
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
    }
    if(fdCutRapV0Max > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: V0 |y| < %g\n", fdCutRapV0Max);
      if(TMath::Abs(dRapK0s) > fdCutRapV0Max)
        bIsCandidateK0s = kFALSE;
      if(TMath::Abs(dRapLambda) > fdCutRapV0Max) {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
    }
    if(fdCutEtaV0Max > 0. || fdCutRapV0Max > 0.)
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 12
    // Lifetime cut
    if(fdCutNTauKMax > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: Proper lifetime < %g (K)\n", fdCutNTauKMax);
      if(dMROverPtK0s > fdCutNTauKMax * dCTauK0s)
        bIsCandidateK0s = kFALSE;
    }
    if(fdCutNTauLMax > 0.) {
      if(bPrintCuts) printf("Rec: Applying cut: Proper lifetime < %g (L, AL)\n", fdCutNTauLMax);
      if(dMROverPtLambda > fdCutNTauLMax * dCTauLambda) {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
    }
    if(fdCutNTauKMax > 0. || fdCutNTauLMax > 0.)
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 13
    // Daughter PID
    if(fdCutNSigmadEdxMax > 0.) {
      if(fbIsPbPb && fdPtProtonPIDMax > 0.) { // Pb-Pb
        if(bPrintCuts) printf("Rec: Applying cut: Delta dE/dx (proton below %g GeV/c) < %g\n", fdPtProtonPIDMax, fdCutNSigmadEdxMax);
        if((dPtDaughterPos < fdPtProtonPIDMax) && (dNSigmaPosProton > fdCutNSigmadEdxMax)) // p+
          bIsCandidateLambda = kFALSE;
        if((dPtDaughterNeg < fdPtProtonPIDMax) && (dNSigmaNegProton > fdCutNSigmadEdxMax)) // p-
          bIsCandidateALambda = kFALSE;
      }
      else { // p-p
        if(bPrintCuts) printf("Rec: Applying cut: Delta dE/dx (both daughters): < %g\n", fdCutNSigmadEdxMax);
        if(dNSigmaPosPion > fdCutNSigmadEdxMax || dNSigmaNegPion > fdCutNSigmadEdxMax) // pi+, pi-
          bIsCandidateK0s = kFALSE;
        if(dNSigmaPosProton > fdCutNSigmadEdxMax || dNSigmaNegPion > fdCutNSigmadEdxMax) // p+, pi-
          bIsCandidateLambda = kFALSE;
        if(dNSigmaNegProton > fdCutNSigmadEdxMax || dNSigmaPosPion > fdCutNSigmadEdxMax) // p-, pi+
          bIsCandidateALambda = kFALSE;
      }
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 14
    // Armenteros-Podolanski cut
    if(fbCutArmPod) {
      if(bPrintCuts) printf("Rec: Applying cut: Armenteros-Podolanski (K0S) pT > %g * |alpha|\n", 0.2);
      if(dPtArm < TMath::Abs(0.2 * dAlpha))
        bIsCandidateK0s = kFALSE;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;
    if(fbCutCross) {
      if(bIsInPeakK0s) {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
      if(bIsInPeakLambda) {
        bIsCandidateK0s = kFALSE;
      }
      if(bIsInPeakALambda) {
        bIsCandidateK0s = kFALSE;
      }
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;
    // End of particle-dependent cuts

    //===== End of reconstruction cutting =====

    if(!bIsCandidateK0s && !bIsCandidateLambda && !bIsCandidateALambda)
      continue;

    if(fbMCAnalysis) {
      AliEmcalJet *ZeroJet = 0;
      AssociateRecV0withMC(v0, ZeroJet, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCentIndex);
    }
    
    Int_t uid = 0; //?
    //===== Start of filling V0 spectra =====
    // iCutIndex = 16
    if(bIsCandidateK0s) {
      // 16 K0s candidates after cuts
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, kFALSE, kFALSE, iCutIndex, iCentIndex);
      Double_t valueKIncl[3] = {dMassV0K0s, dPtV0, dEtaV0};
      fhnV0InclusiveK0s[iCentIndex]->Fill(valueKIncl);
      fh1V0InvMassK0sCent[iCentIndex]->Fill(dMassV0K0s);
      if(bIsInPeakK0s){
        fhnV0InvMassCutK0s[iCentIndex]->Fill(valueKIncl);
      }  
      uid = iK0Id + fNCand;
      dEnergy = v0->EK0Short();
      iNV0CandK0s++;
    }   

    if(bIsCandidateLambda) {
      // 16 Lambda candidates after cuts
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, bIsCandidateLambda, kFALSE, iCutIndex, iCentIndex);
      Double_t valueLIncl[3] = {dMassV0Lambda, dPtV0, dEtaV0};
      fhnV0InclusiveLambda[iCentIndex]->Fill(valueLIncl);
      fh1V0InvMassLambdaCent[iCentIndex]->Fill(dMassV0Lambda);
      if(bIsInPeakLambda){
        fhnV0InvMassCutLambda[iCentIndex]->Fill(valueLIncl);
        uid = iLambdaId + fNCand; 
        dEnergy = v0->ELambda();
      }  
      iNV0CandLambda++;
    }

    if(bIsCandidateALambda) {
      // 16 ALambda candidates after cuts
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, kFALSE, bIsCandidateALambda, iCutIndex, iCentIndex);
      Double_t valueALIncl[3] = {dMassV0ALambda, dPtV0, dEtaV0};
      fhnV0InclusiveALambda[iCentIndex]->Fill(valueALIncl);
      fh1V0InvMassALambdaCent[iCentIndex]->Fill(dMassV0ALambda);  
      if(bIsInPeakALambda){
        fhnV0InvMassCutALambda[iCentIndex]->Fill(valueALIncl);
        uid = iALambdaId + fNCand;
        dEnergy = v0->ELambda();
      } 
      iNV0CandALambda++;    
    }

    if(bIsCandidateK0s && bIsCandidateLambda) {
      if(bIsInPeakK0s){
        uid = iK0Id + fNCand;
        dEnergy = v0->EK0Short();
      }
      if(bIsInPeakLambda){
        uid = iLambdaId + fNCand; 
        dEnergy = v0->ELambda();
      }
      if(bIsInPeakK0s && bIsInPeakLambda)
        uid = iK0LId + fNCand;
    }	

    if(bIsCandidateK0s && bIsCandidateALambda) {
      if(bIsInPeakK0s){
        uid = iK0Id + fNCand;
        dEnergy = v0->EK0Short();
      }
      if(bIsInPeakALambda){
        uid = iALambdaId + fNCand; 
        dEnergy = v0->ELambda();
      }
      if(bIsInPeakK0s && bIsInPeakALambda)
        uid = iK0ALId + fNCand;
    }	

    if(bIsInPeakK0s || bIsInPeakLambda || bIsInPeakALambda){
      if (uid > 0) { // to get rid of the situations when isinpeakcandidate is different from the iscandidate
        new ((*fV0CandidateArray)[fNCand]) AliAODv0(*v0); //  
        ivecV0CandIndex.push_back(uid);
        //add the v0 vector to the fastjetwrapper
        fFastJetWrapper.AddInputVector(vecV0Momentum[0], vecV0Momentum[1], vecV0Momentum[2], dEnergy, uid);
        InputBgParticles.push_back(fastjet::PseudoJet(vecV0Momentum[0], vecV0Momentum[1], vecV0Momentum[2], dEnergy));
      
        fNCand++;

        iNSelV0++;  
      }  
    }

    //===== End of filling V0 spectra =====            
  }	 
  //===== End of V0 loop =====
  //printf("fjw inputs before add track: %i \n", fFastJetWrapper.GetInputVectors().size());
  if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "End of V0 loop");
  //printf("There are %i selected V0s in the event. \n", iNSelV0);   
  if(fV0CandidateArray->GetEntriesFast() > 0) {
    AddEventTracks(fV0CandidateArray, fTracks, InputBgParticles);
  }
  if(iNSelV0)
    fh1EventCounterCut->Fill(10);
  //printf("fjw inputs after add track: %i \n", fFastJetWrapper.GetInputVectors().size());
  fh1V0CandPerEvent->Fill(iNV0CandTot);
  fh1V0CandPerEventCentK0s[iCentIndex]->Fill(iNV0CandK0s);
  fh1V0CandPerEventCentLambda[iCentIndex]->Fill(iNV0CandLambda);
  fh1V0CandPerEventCentALambda[iCentIndex]->Fill(iNV0CandALambda);
  
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

  //if (fFastJetWrapper.GetInputVectors().size() == 0) 

  //Run fj wrapper and loop over fj jets
  fFastJetWrapper.Run(); 

  //printf("All found jets: %i \n", fFastJetWrapper.GetInclusiveJets().size());

  std::vector<fastjet::PseudoJet> vJetsIncl = fFastJetWrapper.GetInclusiveJets(); 
  if(vJetsIncl.size()) 
	  fh1EventCounterCut->Fill(11);
  
  // sort jets according to jet pt
  static Int_t indexes[9999] = {-1};
  GetSortedArray(indexes, vJetsIncl); 
  AliDebug(1,Form("%i jets found", (Int_t)vJetsIncl.size()));
  
  Int_t iJetCount = 0;
  Int_t iV0Jets = 0;
  fastjet::PseudoJet jetSub;

  for (UInt_t ijet = 0; ijet < vJetsIncl.size(); ++ijet) {
    Int_t ij = indexes[ijet];
    AliDebug(3,Form("Jet pt = %f, area = %f", vJetsIncl[ij].perp(), fFastJetWrapper.GetJetArea(ij)));
    if (vJetsIncl[ij].perp() < fdMinJetPt) continue;
    if (fFastJetWrapper.GetJetArea(ij) < fdMinJetArea) continue;
    if ((vJetsIncl[ij].eta() < fdJetEtaMin) || (vJetsIncl[ij].eta() > fdJetEtaMax) ||
        (vJetsIncl[ij].phi() < fdJetPhiMin) || (vJetsIncl[ij].phi() > fdJetPhiMax))
      continue;

    jetSub = subtr(vJetsIncl[ij]); //background subtraction
    if(jetSub == 0) 
      continue;
    fh1NV0sInJetStats->Fill(8);  //N of jets before cuts (after bg sub)
    
    if(jetSub.perp() < fdCutPtJetMin) // selection of high-pt jets, needs to be applied on the pt after bg subtraction
      continue;
    fh1NV0sInJetStats->Fill(9);
    
    if(jetSub.area() < dAreaPercJetMin) //selection of the jets with area bigger than the cut (cut*pi*R2)
      continue;
    fh1NV0sInJetStats->Fill(10);
    
    std::vector<fastjet::PseudoJet> constituents(jetSub.constituents());  //(fFastJetWrapper.GetJetConstituents(ij));
    Int_t iNConstit = constituents.size();
    Double_t dMaxTrPt = 0;
    for(Int_t ic = 0; ic < iNConstit; ic++) {
		  Double_t dtrpt = constituents[ic].perp();
      if(dtrpt > dMaxTrPt) dMaxTrPt = dtrpt;  
    }  
    if(dMaxTrPt < fdCutPtTrackJetMin)             // selection of jets with high leading track pt
      continue;                                           
    fh1NV0sInJetStats->Fill(11);

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
    dAreaExcluded -= AreaCircSegment(dRadiusExcludeCone, fdCutEtaV0Max - jet->Eta()); // positive eta overhang
    dAreaExcluded -= AreaCircSegment(dRadiusExcludeCone, fdCutEtaV0Max + jet->Eta()); // negative eta overhang
    fh1AreaExcluded->Fill(iCentIndex, dAreaExcluded);
    
    fh1NV0sInJetStats->Fill(0);
    if(fbMCAnalysis) 
      fh1MCStats->Fill(1); //Reconstructed jets
    iJetCount++;
   //printf("JetPt: %f, sub jetPt: %f \n", vJetsIncl[ij].perp(), jetSub.perp());
  
    Int_t uid   = -1;
    Int_t ik = 0;
    Int_t il = 0;
    Int_t ial = 0; 
    Int_t ikl = 0;
    Int_t ikal = 0;
    Int_t index = 0;
    AliAODv0* jetv0 = 0;
  
    for(Int_t ic = 0; ic < iNConstit; ic++){
		  uid = constituents[ic].user_index();
		  //if(constituents[ic].user_index() > 0) cout << ic << " Pt:"  << constituents[ic].perp() << " uid: " << constituents[ic].user_index() << "         "  ;
      if(uid > iK0Id && uid < iLambdaId) { //K0
			  index = uid-iK0Id; //will give the id of the V0 particle 
		    jetv0 = (AliAODv0*)fV0CandidateArray->At(index); 
		    Double_t valueKInJet[4] = {jetv0->MassK0Short(), TMath::Sqrt(jetv0->Pt2V0()), jetv0->Eta(), jet->Pt()};
        fhnV0InJetK0s[iCentIndex]->Fill(valueKInJet);
        fh1NV0sInJetStats->Fill(2);

        if(fbMCAnalysis) {
          AssociateRecV0withMC(jetv0, jet, true, false, false, iCentIndex);
	      }  
		   ik++;
		  }	
		  
		  if(uid > iLambdaId && uid < iALambdaId) {  //Lambda	
		    index = uid-iLambdaId; 
		    jetv0 = (AliAODv0*)fV0CandidateArray->At(index); 
		    Double_t valueLInJet[4] = {jetv0->MassLambda(), TMath::Sqrt(jetv0->Pt2V0()), jetv0->Eta(), jet->Pt()};
        fhnV0InJetLambda[iCentIndex]->Fill(valueLInJet);
        fh1NV0sInJetStats->Fill(3);
        if(fbMCAnalysis) {
         AssociateRecV0withMC(jetv0, jet, false, true, false, iCentIndex);
	      }   
        il++;
		  }
		  if(uid > iALambdaId && uid < iK0LId) { //ALambda
		    index = uid-iALambdaId;  
		    jetv0 = (AliAODv0*)fV0CandidateArray->At(index); 
		  	Double_t valueLInJet[4] = {jetv0->MassAntiLambda(), TMath::Sqrt(jetv0->Pt2V0()), jetv0->Eta(), jet->Pt()};
        fhnV0InJetALambda[iCentIndex]->Fill(valueLInJet);
        fh1NV0sInJetStats->Fill(4);
        if(fbMCAnalysis) {
         AssociateRecV0withMC(jetv0, jet, false, false, true, iCentIndex);  
	      }
        ial++;
		  }  
		  if(uid > iK0LId && uid < iK0ALId) { //K L correlation 	
		    index = uid-iK0LId; 
        jetv0 = (AliAODv0*)fV0CandidateArray->At(index); 
		    //Double_t valueKInJet[4] = {jetv0->MassK0Short(), TMath::Sqrt(jetv0->Pt2V0()), jetv0->Eta(), jet->Pt()};
        //fhnV0InJetK0s[iCentIndex]->Fill(valueKInJet);  		  
		    //Double_t valueLInJet[4] = {jetv0->MassLambda(), TMath::Sqrt(jetv0->Pt2V0()), jetv0->Eta(), jet->Pt()};
        //fhnV0InJetLambda[iCentIndex]->Fill(valueLInJet);
        fh1NV0sInJetStats->Fill(5);
        ikl++;
		  }
		  if(uid > iK0ALId) { //K AL correlation 	
		    index = uid-iK0ALId;  
		    jetv0 = (AliAODv0*)fV0CandidateArray->At(index); 
        //Double_t valueKInJet[4] = {jetv0->MassK0Short(), TMath::Sqrt(jetv0->Pt2V0()), jetv0->Eta(), jet->Pt()};
        //fhnV0InJetK0s[iCentIndex]->Fill(valueKInJet);
        //Double_t valueLInJet[4] = {jetv0->MassAntiLambda(), TMath::Sqrt(jetv0->Pt2V0()), jetv0->Eta(), jet->Pt()};
        //fhnV0InJetALambda[iCentIndex]->Fill(valueLInJet);
        fh1NV0sInJetStats->Fill(6);
		    ikal++;
		  }		 
     }
    Int_t isum = ik + il + ial + ikl + ikal;
    if(isum) {
	    iV0Jets++;
	    fh1NV0sInJetStats->Fill(1);
	  }
	  if(isum > 1)
	    fh1NV0sInJetStats->Fill(7);	  		  
  }
  //printf("All selected jets: %i \n",  iJetCount);
  if(iJetCount) {
    fh1EventCent2Jets->Fill(fdCentrality);
    fh1EventCounterCut->Fill(12);
  }  
  else 
    fh1EventCent2NoJets->Fill(fdCentrality);
  //printf("All selected jets with V0: %i \n",  iV0Jets);
  if(iV0Jets) 
    fh1EventCounterCut->Fill(13); // events with jets with V0s
  fh1NJetPerEvent[iCentIndex]->Fill(iJetCount);
  fh1NV0JetPerEvent[iCentIndex]->Fill(iV0Jets);  

  if(iJetCount) {
    jetRnd = GetRandomCone( fJets, dCutEtaJetMax, 2 * fdDistanceV0JetMax); //max jet pseudorap (fdCutEtaV0Max - fdDistanceV0JetMax)
    if(jetRnd) {
      fh1NRndConeCent->Fill(iCentIndex);
      fh2EtaPhiRndCone[iCentIndex]->Fill(jetRnd->Eta(), jetRnd->Phi());
    }
    /*if(fJetsBgCont) // 
    {
      jetMed = GetMedianCluster(fJetsBgCont, dCutEtaJetMax);
      if(jetMed)
      {
        fh1NMedConeCent->Fill(iCentIndex);
        fh2EtaPhiMedCone[iCentIndex]->Fill(jetMed->Eta(), jetMed->Phi());
      }
    }*/
  }
  //Find v0s in the perp and random jet cones 
  if(fV0CandidateArray->GetEntriesFast() > 0) {
    for(Int_t iv0 = 0; iv0 < fV0CandidateArray->GetEntriesFast(); iv0++){
      AliAODv0* v0cand = (AliAODv0*)fV0CandidateArray->At(iv0);
      Int_t iind = ivecV0CandIndex[iv0] - iv0;  
      Double_t dV0MassK = v0cand->MassK0Short();
      Double_t dV0MassL = v0cand->MassLambda();
      Double_t dV0MassAL = v0cand->MassAntiLambda();
      Double_t dPtv =  TMath::Sqrt(v0->Pt2V0());
      Double_t dEtav = v0cand->Eta();

      //Selection of v0 in no-jet events
      if(!iJetCount) {
        if(iind == iK0Id) {
          Double_t valueKNoJet[3] = {dV0MassK, dPtv,dEtav};
          fhnV0NoJetK0s[iCentIndex]->Fill(valueKNoJet);
        }
        if(iind == iLambdaId) {
          Double_t valueLNoJet[3] = {dV0MassL, dPtv,dEtav};
          fhnV0NoJetLambda[iCentIndex]->Fill(valueLNoJet);
        }
        if(iind == iALambdaId) {
          Double_t valueALNoJet[3] = {dV0MassAL, dPtv,dEtav};
          fhnV0NoJetALambda[iCentIndex]->Fill(valueALNoJet);
        }
      }
      // Selection of V0s outside jet cones
      if(!OverlapWithJets(fJets, v0cand, 2 * fdDistanceV0JetMax)) { 
        if(iind == iK0Id) {
          Double_t valueKOutJC[3] = {dV0MassK, dPtv,dEtav};
          fhnV0OutJetK0s[iCentIndex]->Fill(valueKOutJC);
        }
        if(iind == iLambdaId) {
          Double_t valueLOutJet[3] =  {dV0MassL, dPtv,dEtav};
          fhnV0OutJetLambda[iCentIndex]->Fill(valueLOutJet);
        }
        if(iind == iALambdaId) {
          Double_t valueALOutJet[3] =  {dV0MassAL, dPtv,dEtav};
          fhnV0OutJetALambda[iCentIndex]->Fill(valueALOutJet);
        }
      }
      // Selection of V0s in perp. cones 
      for(Int_t iJet = 0; iJet < iNJetPerp; iJet++) {
        jetPerp = (AliAODJet*)arrayJetPerp->At(iJet); // load a jet in the list   
        if(IsParticleInCone(v0cand, jetPerp, fdDistanceV0JetMax)) { // V0 in perp. cone
          if(iind == iK0Id) {
          Double_t valueKInPC[4] = {dV0MassK, dPtv, dEtav, jetPerp->Pt()};
          fhnV0InPerpK0s[iCentIndex]->Fill(valueKInPC);
          }
          if(iind == iLambdaId) {
          Double_t valueLInPC[4] = {dV0MassL, dPtv, dEtav, jetPerp->Pt()};
          fhnV0InPerpLambda[iCentIndex]->Fill(valueLInPC);
          }
          if(iind == iALambdaId) {
          Double_t valueALInPC[4] = {dV0MassAL, dPtv, dEtav, jetPerp->Pt()};
          fhnV0InPerpALambda[iCentIndex]->Fill(valueALInPC);
          }
          break;
        }
      }
      // Selection of V0s in random cones
      if(jetRnd) {
        if(IsParticleInCone(v0, jetRnd, fdDistanceV0JetMax)) {
          if(iind == iK0Id) {
            Double_t valueKInRnd[3] = {dV0MassK, dPtv,dEtav};
            fhnV0InRndK0s[iCentIndex]->Fill(valueKInRnd);
          }
          if(iind == iLambdaId) {
            Double_t valueLInRnd[3] = {dV0MassL, dPtv,dEtav};
            fhnV0InRndLambda[iCentIndex]->Fill(valueLInRnd);
          }
          if(iind == iALambdaId) {
            Double_t valueALInRnd[3] = {dV0MassAL, dPtv,dEtav};
            fhnV0InRndALambda[iCentIndex]->Fill(valueALInRnd);
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

void AliAnalysisTaskStrangenessInJets::FillCandidates(Double_t mK, Double_t mL, Double_t mAL, Bool_t isK, Bool_t isL, Bool_t isAL, Int_t iCut/*cut index*/, Int_t iCent/*cent index*/)
{
  if(isK) {
    fh1V0CounterCentK0s[iCent]->Fill(iCut);
    fh1V0InvMassK0sAll[iCut]->Fill(mK);
  }
  if(isL) {
    fh1V0CounterCentLambda[iCent]->Fill(iCut);
    fh1V0InvMassLambdaAll[iCut]->Fill(mL);
  }
  if(isAL) {
    fh1V0CounterCentALambda[iCent]->Fill(iCut);
    fh1V0InvMassALambdaAll[iCut]->Fill(mAL);
  }
}

Bool_t AliAnalysisTaskStrangenessInJets::IsParticleInCone(const AliVParticle* part1, const AliVParticle* part2, Double_t dRMax) const
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
Bool_t AliAnalysisTaskStrangenessInJets::OverlapWithJets(const TClonesArray* array, const AliVParticle* part, Double_t dDistance) const
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
AliAODJet* AliAnalysisTaskStrangenessInJets::GetRandomCone(const TClonesArray* array, Double_t dEtaConeMax, Double_t dDistance) const
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

Double_t AliAnalysisTaskStrangenessInJets::AreaCircSegment(Double_t dRadius, Double_t dDistance) const
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

Bool_t AliAnalysisTaskStrangenessInJets::IsSelectedForAnalysis() //Event selection
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

Int_t AliAnalysisTaskStrangenessInJets::GetCentralityBinIndex(Double_t centrality)
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

Int_t AliAnalysisTaskStrangenessInJets::GetCentralityBinEdge(Int_t index)
{
// returns the upper edge of the centrality bin corresponding to the provided value of index
  if(index < 0 || index >= fgkiNBinsCent)
    return -1;
  return fgkiCentBinRanges[index];
}

TString AliAnalysisTaskStrangenessInJets::GetCentBinLabel(Int_t index)
{
// get string with centrality range for given bin
  TString lowerEdge = ((index == 0) ? "0" : Form("%d", GetCentralityBinEdge(index - 1)));
  TString upperEdge = Form("%d", GetCentralityBinEdge(index));
  return Form("%s-%s %%", lowerEdge.Data(), upperEdge.Data());
}

Bool_t AliAnalysisTaskStrangenessInJets::GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const
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


Bool_t AliAnalysisTaskStrangenessInJets::IsFromGoodGenerator(Int_t index)
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

void AliAnalysisTaskStrangenessInJets::AddEventTracks(TClonesArray* coll, TClonesArray* tracks, std::vector<fastjet::PseudoJet>& VectorBgPart)
{ 
  // Add event tracks to a collection that already contains the V0 candidates, excluding the daughters of the V0 candidates

  if (!tracks) {
	printf("Track container was not found. Function AddEventTracks will not run! \n");
	return;
  }
  
  TObjArray allDaughters(10);
  allDaughters.SetOwner(kFALSE);

  TIter next(coll);
  AliAODv0* v0part = 0;
  while ((v0part = static_cast<AliAODv0*>(next()))) {
    AliDebug(2, Form("Found a V0 candidtate with pT = %.3f, eta = %.3f, phi = %.3f \n", v0part->Pt(), v0part->Eta(), v0part->Phi()));
    if (v0part) AddDaughters(v0part, allDaughters);
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

Double_t AliAnalysisTaskStrangenessInJets::AddDaughters(AliAODRecoDecay* cand, TObjArray& daughters)
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
void AliAnalysisTaskStrangenessInJets::AddEventTracksMC(TClonesArray* coll, TClonesArray* tracks, std::vector<fastjet::PseudoJet>& VectorBgPartMC, TClonesArray* GenXi)
{ 
  if (!tracks) {
	printf("Track container was not found. Function AddEventTracks will not run! \n");
	return;
  }

  std::vector<Int_t>  vecDaughterLabels; //vector with labels of daughter particles
  //std::vector<Int_t>  vecXiLabels;
  //std::vector<Int_t>  vecAXiLabels;
  //std::vector<Int_t>  vecXi0Labels;
  //std::vector<Int_t>  vecAXi0Labels;
  TIter next(coll);
  AliAODMCParticle* v0part = 0;
  while ((v0part = static_cast<AliAODMCParticle*>(next()))) {
    AliDebug(2, Form("Found a MC generated V0 with pT = %.3f, eta = %.3f, phi = %.3f \n", v0part->Pt(), v0part->Eta(), v0part->Phi())); 
    Int_t n = v0part->GetNDaughters();  //why 3 daughters occur? 
    for (Int_t i = 0; i < n; i++) {
      Int_t iDLabel = v0part->GetDaughterLabel(i);
      vecDaughterLabels.push_back(iDLabel);
    } 
  }

  //in case Xi, Xi0 add additional index to the iN (for FD plot)
  /*TIter nextXi(GenXi);
  AliAODMCParticle* xipart = 0;
  while ((xipart = static_cast<AliAODMCParticle*>(nextXi()))) {
    AliDebug(2, Form("Found a MC generated Xi(Xi0) with pT = %.3f, eta = %.3f, phi = %.3f \n", xipart->Pt(), xipart->Eta(), xipart->Phi())); 
    Int_t iPdgCodePartMC = xipart->GetPdgCode();
    Int_t iLabel = xipart->GetLabel();
    if(iPdgCodePartMC  == 3312)
      vecXiLabels.push_back(iLabel);
    else if(iPdgCodePartMC  == -3312)  
      vecAXiLabels.push_back(iLabel);
    else if(iPdgCodePartMC  == 3322)
      vecXi0Labels.push_back(iLabel);
    else if(iPdgCodePartMC  == -3322)  
      vecAXi0Labels.push_back(iLabel);      
  }*/

 
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
    /*Int_t iindex = 0;
    if (std::find(vecXiLabels.begin(), vecXiLabels.end(), iTrackLabel) != vecXiLabels.end()) 
	    iindex = iXiId; 
    if (std::find(vecAXiLabels.begin(), vecAXiLabels.end(), iTrackLabel) != vecAXiLabels.end()) 
	    iindex = iAXiId;
    if (std::find(vecXi0Labels.begin(), vecXi0Labels.end(), iTrackLabel) != vecXi0Labels.end()) 
	    iindex = iXi0Id; 
    if (std::find(vecAXi0Labels.begin(), vecAXi0Labels.end(), iTrackLabel) != vecAXi0Labels.end()) 
	    iindex = iAXi0Id; */

    //did not find track in the daughters, adding track to the fastjet
    fFastJetWrapperMCGen.AddInputVector(track->Px(), track->Py(), track->Pz(), track->E(), iN); //+ iindex);
    VectorBgPartMC.push_back(fastjet::PseudoJet(track->Px(), track->Py(), track->Pz(), track->E()));
    iN++;
    nadded++;
  } 
  //printf("There were %i V0s, %i daughters, %d tracks, %d were added to the fj, excluded : %i. \n", coll->GetEntriesFast(), inlabels, numbtrack, nadded, nexcl);
}


Bool_t AliAnalysisTaskStrangenessInJets::AssociateRecV0withMC( AliAODv0* vpart, AliEmcalJet *xjet, Bool_t bIsK, Bool_t bIsL, Bool_t bIsAL, Int_t iCent /*centrality index*/)
{
  if(!vpart) {
    printf("CANNOT FIND V0!!!\n");	
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

  const AliAODTrack* trP = (AliAODTrack*)vpart->GetDaughter(0); // positive daughter track
  const AliAODTrack* trN = (AliAODTrack*)vpart->GetDaughter(1); // negative daughter track
  Int_t iLabelPos = TMath::Abs(trP->GetLabel());
  Int_t iLabelNeg = TMath::Abs(trN->GetLabel());

  Double_t dMK0s = vpart->MassK0Short();
  Double_t dMLambda = vpart->MassLambda();
  Double_t dMALambda = vpart->MassAntiLambda();
  Double_t dptv0 = TMath::Sqrt(vpart->Pt2V0()); // transverse momentum of V0

  // Make sure MC daughters are in the array range
  if((iLabelNeg < 0) || (iLabelNeg >= iNTracksMC) || (iLabelPos < 0) || (iLabelPos >= iNTracksMC))
    return kFALSE;
  // Get MC particles corresponding to reconstructed daughter tracks
  AliAODMCParticle* particleMCDaughterNeg = (AliAODMCParticle*)MCPartArray->At(iLabelNeg);
  AliAODMCParticle* particleMCDaughterPos = (AliAODMCParticle*)MCPartArray->At(iLabelPos);
 
  if(!particleMCDaughterNeg || !particleMCDaughterPos)
    return kFALSE;
  // Make sure MC daughter particles are not physical primary
  if((particleMCDaughterNeg->IsPhysicalPrimary()) || (particleMCDaughterPos->IsPhysicalPrimary()))
    return kFALSE;
  // Get identities of MC daughter particles
  Int_t iPdgCodeDaughterPos = particleMCDaughterPos->GetPdgCode();
  Int_t iPdgCodeDaughterNeg = particleMCDaughterNeg->GetPdgCode();
  // Get index of the mother particle for each MC daughter particle
  Int_t iIndexMotherPos = particleMCDaughterPos->GetMother();
  Int_t iIndexMotherNeg = particleMCDaughterNeg->GetMother();
  if((iIndexMotherNeg < 0) || (iIndexMotherNeg >= iNTracksMC) || (iIndexMotherPos < 0) || (iIndexMotherPos >= iNTracksMC))
    return kFALSE;
  // Check whether MC daughter particles have the same mother
  if(iIndexMotherNeg != iIndexMotherPos)
    return kFALSE;
  // Get the MC mother particle of both MC daughter particles
  AliAODMCParticle* particleMCMother = (AliAODMCParticle*)MCPartArray->At(iIndexMotherPos);
  if(!particleMCMother)
    return kFALSE;
  // Get identity of the MC mother particle
  Int_t iPdgCodeMother = particleMCMother->GetPdgCode();
  // Skip not interesting particles
  if((iPdgCodeMother != iPdgCodeK0s) && (TMath::Abs(iPdgCodeMother) != iPdgCodeLambda))
    return kFALSE;
  
  // Check identity of the MC mother particle and the decay channel
  // Is MC mother particle K0S?
  Bool_t bV0MCIsK0s = ((iPdgCodeMother == iPdgCodeK0s) && (iPdgCodeDaughterPos == +iPdgCodePion) && (iPdgCodeDaughterNeg == -iPdgCodePion));
  // Is MC mother particle Lambda?
  Bool_t bV0MCIsLambda = ((iPdgCodeMother == +iPdgCodeLambda) && (iPdgCodeDaughterPos == +iPdgCodeProton) && (iPdgCodeDaughterNeg == -iPdgCodePion));
  // Is MC mother particle anti-Lambda?
  Bool_t bV0MCIsALambda = ((iPdgCodeMother == -iPdgCodeLambda) && (iPdgCodeDaughterPos == +iPdgCodePion) && (iPdgCodeDaughterNeg == -iPdgCodeProton));
  
  Double_t dPtV0Gen = particleMCMother->Pt();
  Double_t dRapV0Gen = particleMCMother->Y();
  Double_t dEtaV0Gen = particleMCMother->Eta();

  // V0 pseudorapidity cut applied on generated particles
  if(fdCutEtaV0Max > 0.) {
    if((TMath::Abs(dEtaV0Gen) > fdCutEtaV0Max))
      return kFALSE;
  }
  // V0 rapidity cut applied on generated particles
  if(fdCutRapV0Max > 0.) {
    if((TMath::Abs(dRapV0Gen) > fdCutRapV0Max))
      return kFALSE;
  }
  // Select only particles from a specific generator
  if(!IsFromGoodGenerator(iIndexMotherPos))
    return kFALSE;

  //For the feed-down correction:
  // Get the MC mother particle of the MC mother particle
  Int_t iIndexMotherOfMother = particleMCMother->GetMother();
  AliAODMCParticle* particleMCMotherOfMother = 0;
  if(iIndexMotherOfMother >= 0)
    particleMCMotherOfMother = (AliAODMCParticle*)MCPartArray->At(iIndexMotherOfMother);
  // Get identity of the MC mother particle of the MC mother particle if it exists
  Int_t iPdgCodeMotherOfMother = 0;
  if(particleMCMotherOfMother)
    iPdgCodeMotherOfMother = particleMCMotherOfMother->GetPdgCode();
  // Check if the MC mother particle of the MC mother particle is a Xi (3322 - Xi0, 3312 - Xi-)
  Bool_t bV0MCComesFromXi = ((particleMCMotherOfMother) && (iPdgCodeMotherOfMother == 3312)); // Is MC mother particle daughter of a Xi?
  Bool_t bV0MCComesFromAXi = ((particleMCMotherOfMother) && (iPdgCodeMotherOfMother == -3312)); // Is MC mother particle daughter of a anti-Xi?
  Bool_t bV0MCComesFromXi0 = ((particleMCMotherOfMother) && (iPdgCodeMotherOfMother == 3322)); // Is MC mother particle daughter of a Xi?
  Bool_t bV0MCComesFromAXi0 = ((particleMCMotherOfMother) && (iPdgCodeMotherOfMother == -3322)); // Is MC mother particle daughter of a anti-Xi?

  // Get the distance between production point of the MC mother particle and the primary vertex
  Double_t dx = dPrimVtxMCX - particleMCMother->Xv();
  Double_t dy = dPrimVtxMCY - particleMCMother->Yv();
  Double_t dz = dPrimVtxMCZ - particleMCMother->Zv();
  Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
  Bool_t bV0MCIsPrimaryDist = (dDistPrimary < fdDistPrimaryMax); // Is close enough to be considered primary-like?

  if(bV0MCIsPrimaryDist && (bV0MCIsK0s || bV0MCIsLambda || bV0MCIsALambda)) {
    fh1MCStats->Fill(0); //Reconstructed V0s
    if(xjet) fh1MCStats->Fill(2); // Reconstrucetd V0s in jets
  }  
  
  // K0s
  if(bIsK) { // selected candidates with any mass
    if(bV0MCIsK0s && bV0MCIsPrimaryDist) {// well reconstructed candidates
      if(xjet == 0) {
        fh2V0K0sPtMassMCRec[iCent]->Fill(dPtV0Gen, dMK0s);
        Double_t valueEtaK[3] = {dMK0s, dPtV0Gen, dEtaV0Gen};
        fh3V0K0sEtaPtMassMCRec[iCent]->Fill(valueEtaK);
        Double_t valueEtaDKNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, 0};
        fhnV0K0sInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKNeg);
        Double_t valueEtaDKPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, 0};
        fhnV0K0sInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKPos);
        fh2V0K0sMCResolMPt[iCent]->Fill(dMK0s - dMassPDGK0s, dptv0);
        fh2V0K0sMCPtGenPtRec[iCent]->Fill(dPtV0Gen, dptv0);
      }
      //here add an array with the well reconstructed particles? 
      else { // true V0 associated to a candidate in jet
        Double_t valueKInJCMC[4] = {dMK0s, dPtV0Gen, dEtaV0Gen, xjet->Pt()};
        fh3V0K0sInJetPtMassMCRec[iCent]->Fill(valueKInJCMC);
        Double_t valueEtaKIn[5] = {dMK0s, dPtV0Gen, dEtaV0Gen, xjet->Pt(), dEtaV0Gen - xjet->Eta()};
        fh4V0K0sInJetEtaPtMassMCRec[iCent]->Fill(valueEtaKIn);

        Double_t valueEtaDKJCNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, xjet->Pt()};
        fhnV0K0sInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKJCNeg);
        Double_t valueEtaDKJCPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, xjet->Pt()};
        fhnV0K0sInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDKJCPos);
      }
    }
    if(bV0MCIsK0s && !bV0MCIsPrimaryDist) { // not primary K0s
      fh1V0K0sPtMCRecFalse[iCent]->Fill(dPtV0Gen);
    }
  }
  // Lambda
  if(bIsL) { // selected candidates with any mass
    if(bV0MCIsLambda && bV0MCIsPrimaryDist) { // well reconstructed candidates
      if(xjet == 0) {
        fh2V0LambdaPtMassMCRec[iCent]->Fill(dPtV0Gen, dMLambda);
        Double_t valueEtaL[3] = {dMLambda, dPtV0Gen, dEtaV0Gen};
        fh3V0LambdaEtaPtMassMCRec[iCent]->Fill(valueEtaL);
        Double_t valueEtaDLNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, 0};
        fhnV0LambdaInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDLNeg);
        Double_t valueEtaDLPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, 0};
        fhnV0LambdaInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDLPos);
        fh2V0LambdaMCResolMPt[iCent]->Fill(dMLambda - dMassPDGLambda, dptv0);
        fh2V0LambdaMCPtGenPtRec[iCent]->Fill(dPtV0Gen, dptv0);
      }
      else { // true V0 associated to a reconstructed candidate in jet
        Double_t valueLInJCMC[4] = {dMLambda, dPtV0Gen, dEtaV0Gen, xjet->Pt()};
        fh3V0LambdaInJetPtMassMCRec[iCent]->Fill(valueLInJCMC);
        Double_t valueEtaLIn[5] = {dMLambda, dPtV0Gen, dEtaV0Gen, xjet->Pt(), dEtaV0Gen - xjet->Eta()};
        fh4V0LambdaInJetEtaPtMassMCRec[iCent]->Fill(valueEtaLIn);
        Double_t valueEtaDLJCNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, xjet->Pt()};
        fhnV0LambdaInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDLJCNeg);
        Double_t valueEtaDLJCPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, xjet->Pt()};
        fhnV0LambdaInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDLJCPos);
      }
    }
    // Fill the feed-down histograms
    if(bV0MCIsLambda && bV0MCComesFromXi) {
      if(xjet == 0) {  
        Double_t valueFDLIncl[3] = {dPtV0Gen, particleMCMotherOfMother->Pt(), 0.};
        fhnV0LambdaInclMCFromXi[iCent]->Fill(valueFDLIncl);
      }
      else {
        Double_t valueFDLInJets[3] = {dPtV0Gen, particleMCMotherOfMother->Pt(), xjet->Pt()};
        fhnV0LambdaInJetsMCFromXi[iCent]->Fill(valueFDLInJets);
      }
    }
    if(bV0MCIsLambda && bV0MCComesFromXi0) {
      if(xjet == 0) {  
        Double_t valueFDLIncl[3] = {dPtV0Gen, particleMCMotherOfMother->Pt(), 0.};
        fhnV0LambdaInclMCFromXi0[iCent]->Fill(valueFDLIncl);
      }
      else {
        Double_t valueFDLInJets[3] = {dPtV0Gen, particleMCMotherOfMother->Pt(), xjet->Pt()};
        fhnV0LambdaInJetsMCFromXi0[iCent]->Fill(valueFDLInJets);
      }
    }
  }
  // anti-Lambda
  if(bIsAL) {// selected candidates with any mass
    if(bV0MCIsALambda && bV0MCIsPrimaryDist) { // well reconstructed candidates
      if(xjet == 0) {    
        fh2V0ALambdaPtMassMCRec[iCent]->Fill(dPtV0Gen, dMALambda);
        Double_t valueEtaAL[3] = {dMALambda, dPtV0Gen, dEtaV0Gen};
        fh3V0ALambdaEtaPtMassMCRec[iCent]->Fill(valueEtaAL);
        Double_t valueEtaDALNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, 0};
        fhnV0ALambdaInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDALNeg);
        Double_t valueEtaDALPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, 0};
        fhnV0ALambdaInclDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDALPos);
        fh2V0ALambdaMCResolMPt[iCent]->Fill(dMALambda - dMassPDGLambda, dptv0);
        fh2V0ALambdaMCPtGenPtRec[iCent]->Fill(dPtV0Gen, dptv0);
      }
      else { // true V0 associated to a reconstructed candidate in jet
        Double_t valueALInJCMC[4] = {dMALambda, dPtV0Gen, dEtaV0Gen, xjet->Pt()};
        fh3V0ALambdaInJetPtMassMCRec[iCent]->Fill(valueALInJCMC);
        Double_t valueEtaALIn[5] = {dMALambda, dPtV0Gen, dEtaV0Gen, xjet->Pt(), dEtaV0Gen - xjet->Eta()};
        fh4V0ALambdaInJetEtaPtMassMCRec[iCent]->Fill(valueEtaALIn);
        Double_t valueEtaDALJCNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, xjet->Pt()};
        fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDALJCNeg);
        Double_t valueEtaDALJCPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, xjet->Pt()};
        fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[iCent]->Fill(valueEtaDALJCPos);
      }
    }
    // Fill the feed-down histograms
    if(bV0MCIsALambda && bV0MCComesFromAXi) {
      if(xjet == 0) {  
        Double_t valueFDALIncl[3] = {dPtV0Gen, particleMCMotherOfMother->Pt(), 0.};
        fhnV0ALambdaInclMCFromAXi[iCent]->Fill(valueFDALIncl);
      }
      else {
        Double_t valueFDLInJets[3] = {dPtV0Gen, particleMCMotherOfMother->Pt(), xjet->Pt()};
        fhnV0ALambdaInJetsMCFromAXi[iCent]->Fill(valueFDLInJets);
      }
    }
    if(bV0MCIsALambda && bV0MCComesFromAXi0) {
      if(xjet == 0) {  
        Double_t valueFDALIncl[3] = {dPtV0Gen, particleMCMotherOfMother->Pt(), 0.};
        fhnV0ALambdaInclMCFromAXi0[iCent]->Fill(valueFDALIncl);
      }
      else {
        Double_t valueFDALInJets[3] = {dPtV0Gen, particleMCMotherOfMother->Pt(), xjet->Pt()};
        fhnV0ALambdaInJetsMCFromAXi0[iCent]->Fill(valueFDALInJets);
      }
    }
  }

  return kTRUE;
}

Bool_t AliAnalysisTaskStrangenessInJets::GeneratedMCParticles(TClonesArray* track, Int_t iCent)
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
    if((TMath::Abs(particleMC->Y()) < 0.5) && IsFromGoodGenerator(iPartMC)) {
      if((iPdgCodeParticleMC == 3312)) {
        fh1V0XiPtMCGen[iCent]->Fill(particleMC->Pt());
        new ((*fGenMCXis)[iNMCCand]) AliAODMCParticle(*particleMC); //  
      }
      else if((iPdgCodeParticleMC == -3312)) {
        fh1V0AXiPtMCGen[iCent]->Fill(particleMC->Pt());
        new ((*fGenMCXis)[iNMCCand]) AliAODMCParticle(*particleMC); //  
      }
      else if((iPdgCodeParticleMC == 3322) ) {
        fh1V0Xi0PtMCGen[iCent]->Fill(particleMC->Pt());
        new ((*fGenMCXis)[iNMCCand]) AliAODMCParticle(*particleMC); //  
      }
      else if((iPdgCodeParticleMC == -3322)) {
        fh1V0AXi0PtMCGen[iCent]->Fill(particleMC->Pt());
        new ((*fGenMCXis)[iNMCCand]) AliAODMCParticle(*particleMC); //  
      }
    }
    //  skip non V0 particles 
    if((iPdgCodeParticleMC != iPdgCodeK0s) && (TMath::Abs(iPdgCodeParticleMC) != iPdgCodeLambda)) {// && (TMath::Abs(iPdgCodeParticleMC) != iPdgCodeXi))
      continue;
    }
    // Check identity of the MC V0 particle
    // Is MC V0 particle K0S?
    Bool_t bV0MCIsK0s = (iPdgCodeParticleMC == iPdgCodeK0s);
    // Is MC V0 particle Lambda?
    Bool_t bV0MCIsLambda = (iPdgCodeParticleMC == +iPdgCodeLambda);
    // Is MC V0 particle anti-Lambda?
    Bool_t bV0MCIsALambda = (iPdgCodeParticleMC == -iPdgCodeLambda);

    Double_t dPtV0Gen = particleMC->Pt();
    Double_t dRapV0Gen = particleMC->Y();
    Double_t dEtaV0Gen = particleMC->Eta();

    // V0 pseudorapidity cut
    if(fdCutEtaV0Max > 0.) {
      if((TMath::Abs(dEtaV0Gen) > fdCutEtaV0Max)) 
        continue;
    }
    // V0 rapidity cut
    if(fdCutRapV0Max > 0.) {
      if((TMath::Abs(dRapV0Gen) > fdCutRapV0Max)) 
        continue;
    }   

    // Get the distance between the production point of the MC V0 particle and the primary vertex
    Double_t dx = dPrimVtxMCX - particleMC->Xv();
    Double_t dy = dPrimVtxMCY - particleMC->Yv();
    Double_t dz = dPrimVtxMCZ - particleMC->Zv();
    Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
    Bool_t bV0MCIsPrimaryDist = (dDistPrimary < fdDistPrimaryMax); // Is close enough to be considered primary-like?
    // Select only primary-like MC V0 particles
    if(!bV0MCIsPrimaryDist)
      continue;

    // Select only particles from a specific generator
    if(!IsFromGoodGenerator(iPartMC))
      continue;
 

    Int_t ind = 0;
    fh1MCStats->Fill(3); //Generated V0s
    // K0s
    if(bV0MCIsK0s) {// well reconstructed candidates
      fh1V0K0sPtMCGen[iCent]->Fill(dPtV0Gen);
      fh2V0K0sEtaPtMCGen[iCent]->Fill(dPtV0Gen, dEtaV0Gen);
      ind = iK0Id + iNMCCand;
    }
    // Lambda
    if(bV0MCIsLambda) { // well reconstructed candidates
      fh1V0LambdaPtMCGen[iCent]->Fill(dPtV0Gen);
      fh2V0LambdaEtaPtMCGen[iCent]->Fill(dPtV0Gen, dEtaV0Gen);
      ind = iLambdaId + iNMCCand;
    }
    // anti-Lambda
    if(bV0MCIsALambda) { // well reconstructed candidates
      fh1V0ALambdaPtMCGen[iCent]->Fill(dPtV0Gen);
      fh2V0ALambdaEtaPtMCGen[iCent]->Fill(dPtV0Gen, dEtaV0Gen);
      ind = iALambdaId + iNMCCand;
    }
    
    new ((*fGenMCV0)[iNMCCand]) AliAODMCParticle(*particleMC); //  
    //add the MC v0 vector to the fastjetwrapper with modified id
    fFastJetWrapperMCGen.AddInputVector(particleMC->Px(), particleMC->Py(), particleMC->Pz(), particleMC->E(), ind);
    InputBgParticlesMC.push_back(fastjet::PseudoJet(particleMC->Px(), particleMC->Py(), particleMC->Pz(), particleMC->E()));
    
    iNMCCand++;
  }

  
  
  if(fGenMCV0->GetEntriesFast() > 0) {
    AddEventTracksMC(fGenMCV0, track, InputBgParticlesMC, fGenMCXis);
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
  //Int_t iV0Jets = 0;
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
      if(uid > iK0Id)
      	fh1MCStats->Fill(5); //Generated V0s in jets
      //Fill histograms for each particle: 	  
      if(uid > iK0Id && uid < iLambdaId) { //K0
		    index = uid-iK0Id; //will give the id of the V0 particle 
		    jetv0 = (AliAODMCParticle*)fGenMCV0->At(index); 
        fh2V0K0sInJetPtMCGen[iCent]->Fill(jetv0->Pt(), jetSubMC.perp());
        Double_t valueEtaKInGen[4] = {jetv0->Pt(), jetv0->Eta(), jetSubMC.perp(), jetv0->Eta() - jetSubMC.eta()};
        fh3V0K0sInJetEtaPtMCGen[iCent]->Fill(valueEtaKInGen);
		   ik++;
	    }	
		  
      if(uid > iLambdaId && uid < iALambdaId) {  //Lambda	
		    index = uid-iLambdaId; 
		    jetv0 = (AliAODMCParticle*)fGenMCV0->At(index); 
        fh2V0LambdaInJetPtMCGen[iCent]->Fill(jetv0->Pt(), jetSubMC.perp());
        Double_t valueEtaLInGen[4] = {jetv0->Pt(), jetv0->Eta(), jetSubMC.perp(), jetv0->Eta() - jetSubMC.eta()};
        fh3V0LambdaInJetEtaPtMCGen[iCent]->Fill(valueEtaLInGen);
        il++;
	    }
	
	    if(uid > iALambdaId && uid < iXiId) {//&& uid < iK0LId) { //ALambda
		    index = uid-iALambdaId; 
        cout << index << endl; 
		    jetv0 = (AliAODMCParticle*)fGenMCV0->At(index); 
        fh2V0ALambdaInJetPtMCGen[iCent]->Fill(jetv0->Pt(), jetSubMC.perp());
        Double_t valueEtaALInGen[4] = {jetv0->Pt(), jetv0->Eta(), jetSubMC.perp(), jetv0->Eta() - jetSubMC.eta()};
        fh3V0ALambdaInJetEtaPtMCGen[iCent]->Fill(valueEtaALInGen);
        ial++;
	    }  
      /*
      if(uid > iXiId && uid < iAXiId) { //for the Feeddown generated Xi plots
        index = uid-iXiId; 
		    jetv0 = (AliAODMCParticle*)fGenMCV0->At(index); 
        printf("There is a Xi in the jet. %d , pT= %f\n", uid, jetv0->Pt());
        fh1V0XiInJetPtMCGen[iCent]->Fill(jetv0->Pt());
      }
      if(uid > iAXiId && uid < iXi0Id) {
        index = uid-iAXiId; 
		    jetv0 = (AliAODMCParticle*)fGenMCV0->At(index); 
        printf("There is a AXi in the jet. %d , pT= %f\n", uid, jetv0->Pt());
        fh1V0AXiPtMCGen[iCent]->Fill(jetv0->Pt());
      }
      if(uid > iXi0Id && uid < iAXi0Id) {
        index = uid-iXi0Id; 
		    jetv0 = (AliAODMCParticle*)fGenMCV0->At(index); 
        printf("There is a Xi0 in the jet. %d , pT= %f\n", uid, jetv0->Pt());
        fh1V0Xi0PtMCGen[iCent]->Fill(jetv0->Pt());
      }
      if(uid > iAXi0Id) {
        index = uid-iAXi0Id; 
		    jetv0 = (AliAODMCParticle*)fGenMCV0->At(index); 
        printf("There is a A0Xi in the jet. %d , pT= %f\n", uid, jetv0->Pt());
        fh1V0AXi0PtMCGen[iCent]->Fill(jetv0->Pt());
      }*/
    //Int_t isum = ik + il + ial;
	  //printf("Tehere were %i K0s, %i Ls and %i ALs, isum: %i \n", ik, il, ial, isum);
    }	    
  } 

  return kTRUE;

}
Double_t AliAnalysisTaskStrangenessInJets::MassPeakSigma(Double_t pt, Int_t particle)
{
// estimation of the sigma of the invariant-mass peak as a function of pT and particle type
  switch(particle) {
    case 0: // K0S
      return 0.0037 + 0.000279 * pt + 0.000021 * pt * pt;//0.00362 + 0.000309 * pt + 0.000016 * pt * pt;// 0.00398 + 0.000103 * pt + 0.000042 * pt * pt;
      break;
    case 1: // Lambda
      return 0.0017 - 0.000118 * pt + 0.000029 * pt * pt;//0.00157 - 0.000026 * pt + 0.000017 * pt * pt; //0.00156 - 0.000021 * pt + 0.000016 * pt * pt;  //old
      break;   
    default:
      return 0;
      break;
  }
}

Double_t AliAnalysisTaskStrangenessInJets::MassPeakMean(Double_t pt, Int_t particle)
{
// estimation of the sigma of the invariant-mass peak as a function of pT and particle type
  switch(particle) {
    case 0: // K0S
      return 0.49874;  //old: 0.0044 + 0.0004 * (pt - 1.);
      break;
    case 1: // Lambda
      return 1.11615; //old: 0.0023 + 0.00034 * (pt - 1.);
      break; 
    default:
      return 0;
      break;
  }
}