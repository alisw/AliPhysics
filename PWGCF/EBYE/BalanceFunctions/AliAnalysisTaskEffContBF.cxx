#include "TChain.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include <AliPID.h>
#include <AliPIDCombined.h>
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h" 
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliCentrality.h"
#include "AliGenEventHeader.h"
#include "AliMultSelection.h"
#include "AliESDtrackCuts.h"

#include "AliLog.h"
#include "AliAnalysisTaskEffContBF.h"

// ---------------------------------------------------------------------
//
// Task for calculating the efficiency of the Balance Function 
// for single particles and pairs
// 
// ---------------------------------------------------------------------

ClassImp(AliAnalysisTaskEffContBF)

// ---------------------------------------------------------------------
AliAnalysisTaskEffContBF::AliAnalysisTaskEffContBF() : AliAnalysisTaskSE(),
    fAOD(0),
    fArrayMC(0),
    fQAList(0),
    fOutputList(0),
    fHistEventStats(0),
    fHistCentrality(0),
    fHistNMult(0),
    fHistVz(0),
    fHistDCA(0),
    fHistNSigmaTPCvsPtbeforePID(0),
    fHistNSigmaTPCvsPtafterPID(0),
    fHistContaminationSecondariesPlus(0),
    fHistContaminationSecondariesMinus(0), //
    fHistContaminationSecondariesMaterialPlus(0),
    fHistContaminationSecondariesMaterialMinus(0), //
    fHistContaminationSecondariesWeakDecPlus(0),
    fHistContaminationSecondariesWeakDecMinus(0), //
    fHistContaminationPrimariesPlus(0),
    fHistContaminationPrimariesMinus(0), //
    fHistGeneratedEtaPtPhiPlus(0),
    fHistSurvivedEtaPtPhiPlus(0),
    fHistGeneratedEtaPtPhiMinus(0),
    fHistSurvivedEtaPtPhiMinus(0),
    fHistGeneratedEtaPtPlusControl(0),
    fHistSurvivedEtaPtPlusControl(0),
    fHistGeneratedEtaPtMinusControl(0),
    fHistSurvivedEtaPtMinusControl(0),
    fHistGeneratedEtaPtPlusPlus(0),
    fHistSurvivedEtaPtPlusPlus(0),
    fHistGeneratedEtaPtMinusMinus(0),
    fHistSurvivedEtaPtMinusMinus(0),
    fHistGeneratedEtaPtPlusMinus(0),
    fHistSurvivedEtaPtPlusMinus(0),
    fHistGeneratedPhiEtaPlusPlus(0),
    fHistSurvivedPhiEtaPlusPlus(0),
    fHistGeneratedPhiEtaMinusMinus(0),
    fHistSurvivedPhiEtaMinusMinus(0),
    fHistGeneratedPhiEtaPlusMinus(0),
    fHistSurvivedPhiEtaPlusMinus(0),
    fUseCentrality(kFALSE),
    fCentralityEstimator("V0M"),
    fCentralityPercentileMin(0.0),
    fCentralityPercentileMax(5.0),
    fInjectedSignals(kFALSE),
    fRejectLabelAboveThreshold(kFALSE),
    fGenToBeKept("Hijing"),
    fRejectCheckGenName(kFALSE),
    fPIDResponse(0),
    fElectronRejection(kFALSE),
    fElectronOnlyRejection(kFALSE),
    fElectronRejectionNSigma(-1.),
    fElectronRejectionMinPt(0.),
    fElectronRejectionMaxPt(1000.),
    fExcludeElectronsInMC(kFALSE),
    fUseY(kFALSE),
    fPIDCombined(0),
    fUsePIDnSigmaComb(kTRUE),
    fBayesPIDThr(0.8),
    fUsePIDstrategy(kFALSE),
    fUsePIDFromPDG(kFALSE),
    fpartOfInterest(AliPID::kPion),
    fPDGCodeWanted(0),
    fMassParticleOfInterest(0.13957),
    fVxMax(3.0),
    fVyMax(3.0),
    fVzMax(10.),
    fDCAxyCut(-1),
    fDCAzCut(-1),
    fTPCchi2Cut(-1),
    fNClustersTPCCut(-1),
    fMinTPCCrossedRows(-1),
    fMinTPCRowsOverFindableCls(-1),
    fAODTrackCutBit(128),
    fMinNumberOfTPCClusters(80),
    fMaxChi2PerTPCCluster(4.0),
    fMaxDCAxy(3.0),
    fMaxDCAz(3.0),
    fMinPt(0.0),
    fMaxPt(20.0),
    fMinEta(-0.8),
    fMaxEta(0.8),
    fEtaRangeMin(0.0),
    fEtaRangeMax(1.6),
    fPtRangeMin(0.0),
    fPtRangeMax(20.0),
    fUseTOFBCPileUpCut(kFALSE),
    fUseTPCInOutRowsCut(kFALSE),
    fInRows(2),
    fOutRows(20), 
    fEtaBin(100), //=100 (BF) 16
    fdEtaBin(64), //=64 (BF)  16
    fPtBin(100), //=100 (BF)  36
    fHistSurvived4EtaPtPhiPlus(0),
    fHistSurvived8EtaPtPhiPlus(0),
    fESDtrackCuts(0x0),
    fUseRaaGeoCut(kFALSE),
    fDeadZoneWidth(3),
    fCutGeoNcrNclLength(130),
    fCutGeoNcrNclGeom1Pt(1.5),
    fCutGeoNcrNclFractionNcr(0.85),
    fCutGeoNcrNclFractionNcl(0.7),
    fHistPdgGen(0),
    fHistPdgSurv(0){
} 

//________________________________________________________________________
AliAnalysisTaskEffContBF::AliAnalysisTaskEffContBF(const char *name)
  : AliAnalysisTaskSE(name), 
    fAOD(0),
    fArrayMC(0), 
    fQAList(0), 
    fOutputList(0), 
    fHistEventStats(0), 
    fHistCentrality(0),
    fHistNMult(0), 
    fHistVz(0),
    fHistDCA(0), 
    fHistNSigmaTPCvsPtbeforePID(0),
    fHistNSigmaTPCvsPtafterPID(0),  
    fHistContaminationSecondariesPlus(0),
    fHistContaminationSecondariesMinus(0), //
    fHistContaminationSecondariesMaterialPlus(0),
    fHistContaminationSecondariesMaterialMinus(0), //
    fHistContaminationSecondariesWeakDecPlus(0),
    fHistContaminationSecondariesWeakDecMinus(0), //
    fHistContaminationPrimariesPlus(0),
    fHistContaminationPrimariesMinus(0), //
    fHistGeneratedEtaPtPhiPlus(0), 
    fHistSurvivedEtaPtPhiPlus(0),
    fHistGeneratedEtaPtPhiMinus(0),
    fHistSurvivedEtaPtPhiMinus(0),
    fHistGeneratedEtaPtPlusControl(0),
    fHistSurvivedEtaPtPlusControl(0),
    fHistGeneratedEtaPtMinusControl(0),
    fHistSurvivedEtaPtMinusControl(0),
    fHistGeneratedEtaPtPlusPlus(0),
    fHistSurvivedEtaPtPlusPlus(0),
    fHistGeneratedEtaPtMinusMinus(0),
    fHistSurvivedEtaPtMinusMinus(0),
    fHistGeneratedEtaPtPlusMinus(0),
    fHistSurvivedEtaPtPlusMinus(0),
    fHistGeneratedPhiEtaPlusPlus(0),
    fHistSurvivedPhiEtaPlusPlus(0),
    fHistGeneratedPhiEtaMinusMinus(0),
    fHistSurvivedPhiEtaMinusMinus(0),
    fHistGeneratedPhiEtaPlusMinus(0),
    fHistSurvivedPhiEtaPlusMinus(0),
    fUseCentrality(kFALSE),
    fCentralityEstimator("V0M"), 
    fCentralityPercentileMin(0.0), 
    fCentralityPercentileMax(5.0), 
    fInjectedSignals(kFALSE),
    fRejectLabelAboveThreshold(kFALSE),
    fGenToBeKept("Hijing"),
    fRejectCheckGenName(kFALSE),
    fPIDResponse(0),
    fElectronRejection(kFALSE),
    fElectronOnlyRejection(kFALSE),
    fElectronRejectionNSigma(-1.),
    fElectronRejectionMinPt(0.),
    fElectronRejectionMaxPt(1000.),
    fExcludeElectronsInMC(kFALSE),
    fUseY(kFALSE),
    fPIDCombined(0),
    fUsePIDnSigmaComb(kTRUE),
    fBayesPIDThr(0.8),
    fUsePIDstrategy(kFALSE),
    fUsePIDFromPDG(kFALSE),
    fpartOfInterest(AliPID::kPion),
    fPDGCodeWanted(0),
    fMassParticleOfInterest(0.13957),
    fVxMax(3.0), 
    fVyMax(3.0),
    fVzMax(10.), 
    fDCAxyCut(-1),
    fDCAzCut(-1),
    fTPCchi2Cut(-1),
    fNClustersTPCCut(-1),
    fMinTPCCrossedRows(-1),
    fMinTPCRowsOverFindableCls(-1),
    fAODTrackCutBit(128),
    fMinNumberOfTPCClusters(80),
    fMaxChi2PerTPCCluster(4.0),
    fMaxDCAxy(3.0),
    fMaxDCAz(3.0),
    fMinPt(0.0),
    fMaxPt(20.0),
    fMinEta(-0.8), 
    fMaxEta(0.8),
    fEtaRangeMin(0.0), 
    fEtaRangeMax(1.6),
    fPtRangeMin(0.0),
    fPtRangeMax(20.0),
    fUseTOFBCPileUpCut(kFALSE),
    fUseTPCInOutRowsCut(kFALSE),
    fInRows(2),
    fOutRows(20), 
    fEtaBin(100), //=100 (BF) 16
    fdEtaBin(64), //=64 (BF)  16
    fPtBin(100), //=100 (BF)  36
    fHistSurvived4EtaPtPhiPlus(0),
    fHistSurvived8EtaPtPhiPlus(0),
    fESDtrackCuts(0x0),
    fUseRaaGeoCut(kFALSE),
    fDeadZoneWidth(3),
    fCutGeoNcrNclLength(130),
    fCutGeoNcrNclGeom1Pt(1.5),
    fCutGeoNcrNclFractionNcr(0.85),
    fCutGeoNcrNclFractionNcl(0.7),
    fHistPdgGen(0),
    fHistPdgSurv(0)
   {   
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEffContBF::UserCreateOutputObjects() {
    
  // Create histograms
  // Called once
  fQAList = new TList();
  fQAList->SetName("QAList");
  fQAList->SetOwner();
  
  fOutputList = new TList();
  fOutputList->SetName("OutputList");
  fOutputList->SetOwner();
  
  //Event stats.
  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
                             "Event statistics;;N_{events}",
                             4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fQAList->Add(fHistEventStats);

  //====================================================//
  Int_t ptBin = 36;
  Int_t etaBin = 16;
  Int_t phiBin = 100;

  Double_t nArrayPt[37]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0};
  Double_t nArrayEta[17]={-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

  Double_t nArrayPhi[phiBin+1];
  for(Int_t iBin = 0; iBin <= phiBin; iBin++) 
    nArrayPhi[iBin] = iBin*TMath::TwoPi()/phiBin;

  Int_t detaBin = 16;
  Int_t dphiBin = 100;
  Double_t nArrayDPhi[dphiBin+1];
  for(Int_t iBin = 0; iBin <= dphiBin; iBin++) 
    nArrayDPhi[iBin] = iBin*TMath::TwoPi()/dphiBin;
  Double_t nArrayDEta[17]={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6};
  //====================================================//

  if (fUsePIDFromPDG){
    fPDGCodeWanted = AliPID::ParticleCode(fpartOfInterest);
    Printf("******************** fPDGCodeWanted =%d ******************", fPDGCodeWanted);
  }
  
  //====================================================//
  //AOD analysis
  fHistCentrality = new TH1F("fHistCentrality",";Centrality bin;Events",
			     1001,-0.5,100.5);
  fQAList->Add(fHistCentrality);
  
  //multiplicity (good MC tracks)
  TString histName;
  histName = "fHistNMult";
  fHistNMult = new TH1F(histName.Data(), 
			";N_{mult.}",
			200,0,20000);
  fQAList->Add(fHistNMult);

  //Vz addition+++++++++++++++++++++++++++++
  fHistVz = new TH1F("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Entries",100,-20.,20.);
  fQAList->Add(fHistVz);

  //Electron cuts -> PID QA
  fHistNSigmaTPCvsPtbeforePID = new TH2F ("NSigmaTPCvsPtbefore","NSigmaTPCvsPtbefore",200, 0, 20, 200, -10, 10); 
  fQAList->Add(fHistNSigmaTPCvsPtbeforePID);

  fHistNSigmaTPCvsPtafterPID = new TH2F ("NSigmaTPCvsPtafter","NSigmaTPCvsPtafter",200, 0, 20, 200, -10, 10); 
  fQAList->Add(fHistNSigmaTPCvsPtafterPID);

  fHistDCA  = new TH2F("fHistDCA","DCA (xy vs. z)",400,-5,5,400,-5,5); 
  fQAList->Add(fHistDCA);

  //Contamination for Secondaries
  fHistContaminationSecondariesPlus = new TH3D("fHistContaminationSecondariesPlus","Secondaries;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistContaminationSecondariesPlus);

  fHistContaminationSecondariesMinus = new TH3D("fHistContaminationSecondariesMinus","Secondaries;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistContaminationSecondariesMinus);

  fHistContaminationSecondariesMaterialPlus = new TH3D("fHistContaminationSecondariesMaterialPlus","Secondaries pos from Material ;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistContaminationSecondariesMaterialPlus);

  fHistContaminationSecondariesMaterialMinus = new TH3D("fHistContaminationSecondariesMaterialMinus","Secondaries neg from Material ;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistContaminationSecondariesMaterialMinus);
    
  fHistContaminationSecondariesWeakDecPlus = new TH3D("fHistContaminationSecondariesWeakDecPlus","Secondaries pos from Weak Decay ;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistContaminationSecondariesWeakDecPlus);

  fHistContaminationSecondariesWeakDecMinus = new TH3D("fHistContaminationSecondariesWeakDecMinus","Secondaries neg from Weak Decay ;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistContaminationSecondariesWeakDecMinus);

  //Contamination for Primaries
  fHistContaminationPrimariesPlus = new TH3D("fHistContaminationPrimariesPlus","Primaries;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistContaminationPrimariesPlus);

  fHistContaminationPrimariesMinus = new TH3D("fHistContaminationPrimariesMinus","Primaries;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistContaminationPrimariesMinus);
  
  //eta vs pt for MC positives
  fHistGeneratedEtaPtPhiPlus = new TH3D("fHistGeneratedEtaPtPhiPlus",
					"Generated positive primaries;#eta;p_{T} (GeV/c);#phi",
					etaBin,nArrayEta, ptBin, nArrayPt,phiBin, nArrayPhi);
  // fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax,fPhiBin,nArrayPhi);
  fOutputList->Add(fHistGeneratedEtaPtPhiPlus);
  fHistSurvivedEtaPtPhiPlus = new TH3D("fHistSurvivedEtaPtPhiPlus",
				       "Survived positive primaries;#eta;p_{T} (GeV/c);#phi",
				       etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistSurvivedEtaPtPhiPlus);
  
  //eta vs pt for MC negatives
  fHistGeneratedEtaPtPhiMinus = new TH3D("fHistGeneratedEtaPtPhiMinus",
					 "Generated positive primaries;#eta;p_{T} (GeV/c);#phi",
					 etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistGeneratedEtaPtPhiMinus);
  fHistSurvivedEtaPtPhiMinus = new TH3D("fHistSurvivedEtaPtPhiMinus",
					"Survived positive primaries;#eta;p_{T} (GeV/c);#phi",
					etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistSurvivedEtaPtPhiMinus);
 
  //eta vs pt for MC positives (control)
  fHistGeneratedEtaPtPlusControl = new TH2F("fHistGeneratedEtaPtPlusControl",
					    "Generated positive primaries;#eta;p_{T} (GeV/c)",
					    etaBin,nArrayEta,ptBin,nArrayPt);
  fOutputList->Add(fHistGeneratedEtaPtPlusControl);
  fHistSurvivedEtaPtPlusControl = new TH2F("fHistSurvivedEtaPtPlusControl",
					   "Survived positive primaries;#eta;p_{T} (GeV/c)",
					   etaBin,nArrayEta,ptBin,nArrayPt);
  fOutputList->Add(fHistSurvivedEtaPtPlusControl);
  
  //eta vs pt for MC negatives (control)
  fHistGeneratedEtaPtMinusControl = new TH2F("fHistGeneratedEtaPtMinusControl",
					     "Generated positive primaries;#eta;p_{T} (GeV/c)",
					     etaBin,nArrayEta,ptBin,nArrayPt);
  fOutputList->Add(fHistGeneratedEtaPtMinusControl);
  fHistSurvivedEtaPtMinusControl = new TH2F("fHistSurvivedEtaPtMinusControl",
					    "Survived positive primaries;#eta;p_{T} (GeV/c)",
					    etaBin,nArrayEta,ptBin,nArrayPt);
  fOutputList->Add(fHistSurvivedEtaPtMinusControl);
  
  //eta vs pt for MC ++
  fHistGeneratedEtaPtPlusPlus = new TH2F("fHistGeneratedEtaPtPlusPlus",
					 "Generated ++ primaries;#Delta#eta;p_{T} (GeV/c)",
					 detaBin,nArrayDEta,ptBin,nArrayPt);
  fOutputList->Add(fHistGeneratedEtaPtPlusPlus);
  fHistSurvivedEtaPtPlusPlus = new TH2F("fHistSurvivedEtaPtPlusPlus",
					"Survived ++ primaries;#Delta#eta;p_{T} (GeV/c)",
					detaBin,nArrayDEta,ptBin,nArrayPt);
  fOutputList->Add(fHistSurvivedEtaPtPlusPlus);
  
  //eta vs pt for MC --
  fHistGeneratedEtaPtMinusMinus = new TH2F("fHistGeneratedEtaPtMinusMinus",
					   "Generated -- primaries;#Delta#eta;p_{T} (GeV/c)",
					   detaBin,nArrayDEta,ptBin,nArrayPt);
  fOutputList->Add(fHistGeneratedEtaPtMinusMinus);
  fHistSurvivedEtaPtMinusMinus = new TH2F("fHistSurvivedEtaPtMinusMinus",
					  "Survived -- primaries;#Delta#eta;p_{T} (GeV/c)",
					  detaBin,nArrayDEta,ptBin,nArrayPt);
  fOutputList->Add(fHistSurvivedEtaPtMinusMinus);
 
  //eta vs pt for MC +-
  fHistGeneratedEtaPtPlusMinus = new TH2F("fHistGeneratedEtaPtPlusMinus",
					  "Generated +- primaries;#Delta#eta;p_{T} (GeV/c)",
					  detaBin,nArrayDEta,ptBin,nArrayPt);
  fOutputList->Add(fHistGeneratedEtaPtPlusMinus);
  fHistSurvivedEtaPtPlusMinus = new TH2F("fHistSurvivedEtaPtPlusMinus",
					 "Survived +- primaries;#Delta#eta;p_{T} (GeV/c)",
					 detaBin,nArrayDEta,ptBin,nArrayPt);
  fOutputList->Add(fHistSurvivedEtaPtPlusMinus);
 
  //=============================//
  //phi vs eta for MC ++
  fHistGeneratedPhiEtaPlusPlus = new TH2F("fHistGeneratedPhiEtaPlusPlus",
					  "Generated ++ primaries;#Delta#phi",
					  dphiBin,nArrayDPhi,detaBin,nArrayDEta);
  fOutputList->Add(fHistGeneratedPhiEtaPlusPlus);
  fHistSurvivedPhiEtaPlusPlus = new TH2F("fHistSurvivedPhiEtaPlusPlus",
					 "Survived ++ primaries;#Delta#phi;#Delta#eta",
					 dphiBin,nArrayDPhi,detaBin,nArrayDEta);
  fOutputList->Add(fHistSurvivedPhiEtaPlusPlus);
  
  //phi vs eta for MC --
  fHistGeneratedPhiEtaMinusMinus = new TH2F("fHistGeneratedPhiEtaMinusMinus",
					    "Generated -- primaries;#Delta#phi;#Delta#eta",
					    dphiBin,nArrayDPhi,detaBin,nArrayDEta);
  fOutputList->Add(fHistGeneratedPhiEtaMinusMinus);
  fHistSurvivedPhiEtaMinusMinus = new TH2F("fHistSurvivedPhiEtaMinusMinus",
					   "Survived -- primaries;#Delta#phi;#Delta#eta",
					   dphiBin,nArrayDPhi,detaBin,nArrayDEta);
  fOutputList->Add(fHistSurvivedPhiEtaMinusMinus);
  
  //phi vs eta for MC +-
  fHistGeneratedPhiEtaPlusMinus = new TH2F("fHistGeneratedPhiEtaPlusMinus",
					   "Generated +- primaries;#Delta#phi;#Delta#eta",
					   dphiBin,nArrayDPhi,detaBin,nArrayDEta);
  fOutputList->Add(fHistGeneratedPhiEtaPlusMinus);
  fHistSurvivedPhiEtaPlusMinus = new TH2F("fHistSurvivedPhiEtaPlusMinus",
					  "Survived +- primaries;#Delta#phi;#Delta#eta",
					  dphiBin,nArrayDPhi,detaBin,nArrayDEta);
  fOutputList->Add(fHistSurvivedPhiEtaPlusMinus);
  
  fHistSurvived4EtaPtPhiPlus = new TH3F("fHistSurvived4EtaPtPhiPlus",
                                        "Survived4 + primaries;#eta;p_{T} (GeV/c);#phi",
                                        etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistSurvived4EtaPtPhiPlus);
  fHistSurvived8EtaPtPhiPlus = new TH3F("fHistSurvived8EtaPtPhiPlus",
					"Survived8 + primaries;#eta;p_{T} (GeV/c);#phi",
					etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistSurvived8EtaPtPhiPlus);
    
  //check of pdg
  fHistPdgGen  = new TH1F("fHistPdgGen","Pdg code distribution;pdg code;Entries",6401,-3200.5,3200.5);
  fOutputList->Add(fHistPdgGen);
  fHistPdgSurv  = new TH1F("fHistPdgSurv","Pdg code distribution;pdg code;Entries",6401,-3200.5,3200.5);
  fOutputList->Add(fHistPdgSurv);
    
  //fQAList->Print();
  //fOutputList->Print(); 
  PostData(1, fQAList);
  PostData(2, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskEffContBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    printf("ERROR: fAOD not available\n");
    return;
  }

  fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));   
  if (!fArrayMC)  
    AliFatal("No array of MC particles found !!!"); // MW  no AliFatal use return values   
  
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    AliError("ERROR: Could not retrieve MC event");
    return;
  }

  // PID Response task active?
  fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
  if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");

  // ==============================================================================================
  // Copy from AliAnalysisTaskPhiCorrelations:
  // For productions with injected signals, figure out above which label to skip particles/tracks
  Int_t skipParticlesAbove = 0;
  if ((fInjectedSignals) && (fRejectLabelAboveThreshold))
  {
    AliGenEventHeader* eventHeader = 0;
    Int_t headers = 0;
    
    // AOD only
    AliAODMCHeader* header = (AliAODMCHeader*) fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!header){
      AliFatal("fInjectedSignals set but no MC header found");
      return;
    }
    
    headers = header->GetNCocktailHeaders();
    eventHeader = header->GetCocktailHeader(0);
    
    
    if (!eventHeader)
      {
	// We avoid AliFatal here, because the AOD productions sometimes have events where the MC header is missing 
	// (due to unreadable Kinematics) and we don't want to loose the whole job because of a few events
	AliError("First event header not found. Skipping this event.");
	return;
      }
    
    skipParticlesAbove = eventHeader->NProduced();
    AliInfo(Form("Injected signals in this event (%d headers). Keeping particles/tracks of %s. Will skip particles/tracks above %d.", headers, eventHeader->ClassName(), skipParticlesAbove)); 
  }
  // ==============================================================================================

  
  
  // arrays for 2 particle histograms
  Int_t nMCLabelCounter         = 0;
  const Int_t maxMCLabelCounter = 20000;
  
  Double_t eta[maxMCLabelCounter];
  Double_t vY[maxMCLabelCounter];
  Double_t pt[maxMCLabelCounter];
  Double_t phi[maxMCLabelCounter];
  Int_t level[maxMCLabelCounter];
  Int_t charge[maxMCLabelCounter];
  
  //AliInfo(Form("%d %d",mcEvent->GetNumberOfTracks(),fAOD->GetNumberOfTracks()));
  fHistEventStats->Fill(1); //all events

  AliAODHeader *headerAOD = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
    if (!headerAOD){
      AliFatal("AOD header found");
      return;
    }
    
    //Centrality stuff
    Double_t nCentrality = 0;
    
    if(fUseCentrality){
      if (fAOD->GetRunNumber()<244824) {
	
	AliCentrality *centrality = headerAOD->GetCentralityP();
	nCentrality =centrality->GetCentralityPercentile(fCentralityEstimator.Data());

        if (nCentrality!=-1){
	
	if(!centrality->IsEventInCentralityClass(fCentralityPercentileMin,
						 fCentralityPercentileMax,
						 fCentralityEstimator.Data()))
	  return;
        }

	else if (nCentrality==-1){

        AliMultSelection *multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
        if(!multSelection) {
          AliWarning("AliMultSelection object not found!");
        }
        
        nCentrality = multSelection->GetMultiplicityPercentile(fCentralityEstimator, kFALSE);
	if ((nCentrality < fCentralityPercentileMin) || (nCentrality >= fCentralityPercentileMax)) return;
	}	

      }
      
      else {
	AliMultSelection *multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
	if(!multSelection) {
	  AliWarning("AliMultSelection object not found!");
	}
	else nCentrality = multSelection->GetMultiplicityPercentile(fCentralityEstimator, kTRUE);

	if ((nCentrality < fCentralityPercentileMin) || (nCentrality >= fCentralityPercentileMax)) return;
      }
      
    fHistEventStats->Fill(2); //triggered + centrality
    fHistCentrality->Fill(nCentrality);
    }
      
    //Printf("Centrality selection: %lf - %lf",fCentralityPercentileMin,fCentralityPercentileMax);
    
    const AliAODVertex *vertex = fAOD->GetPrimaryVertex(); 
    if(vertex) {
      if(vertex->GetNContributors() > 0) {
	Double32_t fCov[6];    
	vertex->GetCovarianceMatrix(fCov);   
	if(fCov[5] != 0) {
	  fHistEventStats->Fill(3); //events with a proper vertex
	  if(TMath::Abs(vertex->GetX()) < fVxMax) {    // antes Xv
	    //Printf("X Vertex: %lf", vertex->GetX());
	    //Printf("Y Vertex: %lf", vertex->GetY());
	    if(TMath::Abs(vertex->GetY()) < fVyMax) {  // antes Yv
	      if(TMath::Abs(vertex->GetZ()) < fVzMax) {  // antes Zv
		//Printf("Z Vertex: %lf", vertex->GetZ());
		
		fHistEventStats->Fill(4); //analyzed events
		fHistVz->Fill(vertex->GetZ());

		if (fUseRaaGeoCut){
		  fESDtrackCuts = new AliESDtrackCuts();
		  fESDtrackCuts->SetCutGeoNcrNcl(fDeadZoneWidth, fCutGeoNcrNclLength, fCutGeoNcrNclGeom1Pt, fCutGeoNcrNclFractionNcr, fCutGeoNcrNclFractionNcl);
		}
		
		
	      //++++++++++++++++++CONTAMINATION++++++++++++++++++//
	      Int_t nGoodAODTracks = fAOD->GetNumberOfTracks();
	      Int_t nMCParticles = mcEvent->GetNumberOfTracks();
	      TArrayI labelMCArray(nMCParticles);
	      
	      for(Int_t jTracks = 0; jTracks < nGoodAODTracks; jTracks++) {
		AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(jTracks));
		if(!track) AliFatal("Not a standard AOD");
		if(!track) continue;
		
		if (!track->TestFilterBit(fAODTrackCutBit))
		  continue;
		
		//acceptance
		
		if (fUseY){
		  if(TMath::Abs(track->Y(fMassParticleOfInterest)) > fMaxEta)
		    continue;	  
		}
		else {
		  if(TMath::Abs(track->Eta()) > fMaxEta)
		    continue;
		}
		
		if((track->Pt() > fMaxPt)||(track->Pt() <  fMinPt))
		  continue;
	
		if( fTPCchi2Cut != -1 && track->Chi2perNDF() > fTPCchi2Cut){
		  continue;
     		}
      		
		if( fNClustersTPCCut != -1 && track->GetTPCNcls() < fNClustersTPCCut){
		  continue;
      		}
		
		if(fMinTPCCrossedRows != -1){
		  if ((Float_t)track->GetTPCNCrossedRows() < (120 - (5/(Float_t)track->Pt())) ){
		    continue;
		  }
		}
		
		if (fMinTPCRowsOverFindableCls != -1){
		  Float_t nTPCCrossedRowsOverFindCls = (((Float_t)track->GetTPCNCrossedRows())/((Float_t)track->GetTPCNclsF()));
		  if (nTPCCrossedRowsOverFindCls < fMinTPCRowsOverFindableCls){
		    continue;
		  }
		}
		
		Double_t pos[3];
      		Double_t v[3];
      		Float_t dcaXY = 0.;
      		Float_t dcaZ  = 0.;

        	if(fAODTrackCutBit == 128){
            	 dcaXY = track->DCA();
		 dcaZ  = track->ZAtDCA();
		}
        	else{
		  vertex->GetXYZ(v);
		  track->GetXYZ(pos);
		  dcaXY  = TMath::Sqrt((pos[0] - v[0])*(pos[0] - v[0]) + (pos[1] - v[1])*(pos[1] - v[1]));
		  dcaZ   = pos[2] - v[2];
        	}
		
		if( fDCAxyCut != -1 && fDCAzCut != -1){
		  if(TMath::Sqrt((dcaXY*dcaXY)/(fDCAxyCut*fDCAxyCut)+(dcaZ*dcaZ)/(fDCAzCut*fDCAzCut)) > 1 ){
		    continue;  // 2D cut
		  }
      		}

		Double_t phiRad = track->Phi();

		if(fUseRaaGeoCut){
		  if (!fESDtrackCuts->IsSelected(track))
		    continue;
		}

		if (fUseTOFBCPileUpCut) {
		  if (!track->GetTOFBunchCrossing()==0)
		    continue;
		}
		 
		if (fUseTPCInOutRowsCut) {
		  const TBits& bmap = track->GetTPCClusterMap();
		  // require at least 20 out of 25 and 3 out of 5 innermost rows
		  int nset25 = 0, nset5 = 0;
		  for (int i=0;i<25; i++) {
		    if (!bmap.TestBitNumber(i)) continue;
		    nset25++;
		    if (i<5) nset5++;
		  }
		  if((nset5<fInRows) || (nset25<fOutRows))
		    continue;
		}
		
		Int_t label = TMath::Abs(track->GetLabel());
		if(label > nMCParticles) continue;
		AliAODMCParticle *AODmcTrack = (AliAODMCParticle*) mcEvent->GetTrack(label); 
		Short_t gAODmcCharge = AODmcTrack->Charge();////
		//fHistContaminationPrimaries->Fill(track->Eta(),track->Pt(),phiDeg);
		//if (!(AODmcTrack->IsPhysicalPrimary())) {
		//fHistContaminationSecondaries->Fill(track->Eta(),track->Pt(),phiDeg);
		//}

		// ==============================================================================================
		// Partial copy from AliAnalyseLeadingTrackUE::RemoveInjectedSignals:
		// Skip tracks that come from injected signals
		if (fInjectedSignals){
		  AliAODMCParticle* mother = AODmcTrack;
		  
		  // find the primary mother (if not already physical primary)
		  while (!((AliAODMCParticle*)mother)->IsPhysicalPrimary())
		    {
		      if (((AliAODMCParticle*)mother)->GetMother() < 0)
			{
			  mother = 0;
			  break;
			}
		      
		      mother = (AliAODMCParticle*) fArrayMC->At(((AliAODMCParticle*)mother)->GetMother());
		      if (!mother)
			break;
		    }
		  
		  if (!mother)
		    {
		      AliError(Form("WARNING: No mother found for particle %d:", AODmcTrack->GetLabel()));
		      continue;
		    }
		  
		  if (fRejectLabelAboveThreshold)
		    if (mother->GetLabel() >= skipParticlesAbove)
		      {
			//AliInfo(Form("Remove particle %d (>= %d)",mother->GetLabel(),skipParticlesAbove));
			continue;
		      }
		  
		  if (fRejectCheckGenName){
		    TString generatorName;
		    Bool_t hasGenerator = mcEvent->GetCocktailGenerator(label,generatorName);
		    if((!hasGenerator) || (!generatorName.Contains(fGenToBeKept.Data())))
		      continue;
		    
		    //Printf("mother =%d, generatorName=%s", label, generatorName.Data()); 
		    
		    
		  }
		}
	      // ==============================================================================================


		//Exclude electrons with PDG
                if(fExcludeElectronsInMC){

                  if(TMath::Abs(AODmcTrack->GetPdgCode()) == 11) continue;
                
		}


		if (fUsePIDFromPDG || fUsePIDstrategy){
		  Int_t pdgcode = AODmcTrack->GetPdgCode();
		  if (TMath::Abs(pdgcode) != fPDGCodeWanted) continue;
		}
		
		if (AODmcTrack->IsPhysicalPrimary()) {
		  if(gAODmcCharge > 0){
            		if (fUseY)
			  fHistContaminationPrimariesPlus->Fill(track->Y(fMassParticleOfInterest),track->Pt(),phiRad);
            		else 
			  fHistContaminationPrimariesPlus->Fill(track->Eta(),track->Pt(),phiRad);
		  }
		  if(gAODmcCharge < 0){
		    if (fUseY)
		      fHistContaminationPrimariesMinus->Fill(track->Y(fMassParticleOfInterest),track->Pt(),phiRad);
		    else
		      fHistContaminationPrimariesMinus->Fill(track->Eta(),track->Pt(),phiRad);
		  }
		}
		
		else{
		  Bool_t isFromMaterial = kFALSE;
		  if (AODmcTrack->IsSecondaryFromMaterial()) isFromMaterial = kTRUE;
		  if(gAODmcCharge > 0){
		    if (fUseY)
		      fHistContaminationSecondariesPlus->Fill(track->Y(fMassParticleOfInterest),track->Pt(),phiRad);
		    else
		      fHistContaminationSecondariesPlus->Fill(track->Eta(),track->Pt(),phiRad);
		    if (isFromMaterial){
		      if (fUseY)
			fHistContaminationSecondariesMaterialPlus->Fill(track->Y(fMassParticleOfInterest),track->Pt(),phiRad);
		      else 
			fHistContaminationSecondariesMaterialPlus->Fill(track->Eta(),track->Pt(),phiRad);
		    }
		    else {
		      if (fUseY)
			fHistContaminationSecondariesWeakDecPlus->Fill(track->Y(fMassParticleOfInterest),track->Pt(),phiRad);
		      else 
			fHistContaminationSecondariesWeakDecPlus->Fill(track->Eta(),track->Pt(),phiRad);
		    }
		  }
		  if(gAODmcCharge < 0){
              		if (fUseY)
			  fHistContaminationSecondariesMinus->Fill(track->Y(fMassParticleOfInterest),track->Pt(),phiRad);
              		else
			  fHistContaminationSecondariesMinus->Fill(track->Eta(),track->Pt(),phiRad);
              	  	if (isFromMaterial){
			  if (fUseY)
			    fHistContaminationSecondariesMaterialMinus->Fill(track->Y(fMassParticleOfInterest),track->Pt(),phiRad);
			  else 
			    fHistContaminationSecondariesMaterialMinus->Fill(track->Eta(),track->Pt(),phiRad);
              	  	}
              	  	else {
			  if (fUseY)
			    fHistContaminationSecondariesWeakDecMinus->Fill(track->Y(fMassParticleOfInterest),track->Pt(),phiRad);
			  else
			    fHistContaminationSecondariesWeakDecMinus->Fill(track->Eta(),track->Pt(),phiRad);
                  	}
		  }
               }

		fHistDCA->Fill(dcaZ,dcaXY);

	      }//loop over tracks
	      
	      //++++++++++++++++++CONTAMINATION++++++++++++++++++//
	      
	      //++++++++++++++++++EFFICIENCY+++++++++++++++++++++//
	      for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
		AliAODMCParticle *mcTrack = (AliAODMCParticle*) mcEvent->GetTrack(iTracks); 
		if (!mcTrack) {
		  AliError(Form("ERROR: Could not receive track %d (mc loop)", iTracks));
		  continue;
		}
		
		//exclude particles generated out of the acceptance
		Double_t vz = mcTrack->Zv();
		if (TMath::Abs(vz) > 50.) continue;
		//acceptance
		
		if (fUseY){		  
		  if(TMath::Abs(mcTrack->Y()) > fMaxEta)
		    continue;  
		}
		else{
		  if(TMath::Abs(mcTrack->Eta()) > fMaxEta)
		    continue;
		}
		
		if((mcTrack->Pt() > fMaxPt)||(mcTrack->Pt() < fMinPt))
		  continue;
		
		if(!mcTrack->IsPhysicalPrimary()) continue;
		
		if (fInjectedSignals){
		  if (fRejectCheckGenName){
		    TString generatorName;
		    Bool_t hasGenerator = mcEvent->GetCocktailGenerator(iTracks,generatorName);
		    if((!hasGenerator) || (!generatorName.Contains(fGenToBeKept.Data())))
		      continue;
		    
		  }  
		}
		
		//Exclude electrons with PDG
                if(fExcludeElectronsInMC){

                  if(TMath::Abs(mcTrack->GetPdgCode()) == 11) continue;

                }

		if (fUsePIDFromPDG || fUsePIDstrategy ){

		  Int_t pdgcode = mcTrack->GetPdgCode();
		  if (TMath::Abs(pdgcode) != fPDGCodeWanted) continue;
		}
		
		fHistPdgGen->Fill(mcTrack->GetPdgCode());

		Short_t gMCCharge = mcTrack->Charge();
		Double_t phiRad = mcTrack->Phi();
		
          	if(gMCCharge > 0){
		  if (fUseY)
		    fHistGeneratedEtaPtPhiPlus->Fill(mcTrack->Y(),
						     mcTrack->Pt(),
						     phiRad);
		  else
		    fHistGeneratedEtaPtPhiPlus->Fill(mcTrack->Eta(),
						     mcTrack->Pt(),
						     phiRad);
		  
         	}	
          	else if(gMCCharge < 0){
		  if (fUseY)
		    fHistGeneratedEtaPtPhiMinus->Fill(mcTrack->Y(),
						      mcTrack->Pt(),
						      phiRad); 
		  else
		    fHistGeneratedEtaPtPhiMinus->Fill(mcTrack->Eta(),
						      mcTrack->Pt(),
						      phiRad);
		  
         	}
		      
		Bool_t labelTPC = kTRUE;
		if(labelTPC) {
		  labelMCArray.AddAt(iTracks,nMCLabelCounter);
		  if(nMCLabelCounter >= maxMCLabelCounter){
		    AliWarning(Form("MC Label Counter > Limit (%d) --> stop loop here",maxMCLabelCounter));
		    break;
		  }
		  //fill the arrays for 2 particle analysis
		  eta[nMCLabelCounter]    = mcTrack->Eta();
          	  vY[nMCLabelCounter]      = mcTrack->Y();
		  pt[nMCLabelCounter]     = mcTrack->Pt();
		  phi[nMCLabelCounter]    = mcTrack->Phi();
		  charge[nMCLabelCounter] = gMCCharge;
		  
		  level[nMCLabelCounter]  = 1;
		  nMCLabelCounter += 1;
		}  
	      }//loop over MC particles
	      fHistNMult->Fill(nMCLabelCounter);
	      
	      //AOD track loop
	      Int_t nGoodTracks = fAOD->GetNumberOfTracks();   
	      TArrayI labelArray(nGoodTracks);
	      Int_t labelCounter = 0;
	      
	      for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
              AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
              if(!trackAOD) continue;
		
              //track cuts
              if (!trackAOD->TestFilterBit(fAODTrackCutBit)) continue;

	      if(fUseRaaGeoCut){
		if (!fESDtrackCuts->IsSelected(trackAOD))
		  continue;
	      }

	      if (fUseTOFBCPileUpCut) {
		if (!trackAOD->GetTOFBunchCrossing()==0)
		  continue;
	      }
	      
	      if (fUseTPCInOutRowsCut) {
		const TBits& bmap = trackAOD->GetTPCClusterMap();
		// require at least 20 out of 25 and 3 out of 5 innermost rows
		int nset25 = 0, nset5 = 0;
		for (int i=0;i<25; i++) {
		  if (!bmap.TestBitNumber(i)) continue;
		  nset25++;
		  if (i<5) nset5++;
		}
		if((nset5<fInRows) || (nset25<fOutRows))
		  continue;
	      }
	     
              Int_t label = TMath::Abs(trackAOD->GetLabel());
              if(IsLabelUsed(labelArray,label)) continue;
              labelArray.AddAt(label,labelCounter);
              labelCounter += 1;

              if (fInjectedSignals){
                  if (fRejectCheckGenName){
                      TString generatorName;
                      Bool_t hasGenerator = mcEvent->GetCocktailGenerator(label,generatorName);
                      if((!hasGenerator) || (!generatorName.Contains(fGenToBeKept.Data()))) continue;
                  }
              }

              //Exclude electrons with PDG
              if(fExcludeElectronsInMC){

                  AliAODMCParticle *trackAODMC = (AliAODMCParticle*) mcEvent->GetTrack(label);
                  if(TMath::Abs(trackAODMC->GetPdgCode()) == 11) continue;

              } 
		

		if( fTPCchi2Cut != -1 && trackAOD->Chi2perNDF() > fTPCchi2Cut){
                  continue;
                }

                if( fNClustersTPCCut != -1 && trackAOD->GetTPCNcls() < fNClustersTPCCut){
                  continue;
                }

		if(fMinTPCCrossedRows != -1){
		  if ((Float_t)trackAOD->GetTPCNCrossedRows() < (120 - (5/(Float_t)trackAOD->Pt())) ){
		    continue;
		  }
		}
		
		if (fMinTPCRowsOverFindableCls != -1){
		  Float_t nTPCCrossedRowsOverFindCls = (((Float_t)trackAOD->GetTPCNCrossedRows())/((Float_t)trackAOD->GetTPCNclsF()));
		  if (nTPCCrossedRowsOverFindCls < fMinTPCRowsOverFindableCls){
		    continue;
		  }
		}
		
		Double_t pos[3];
                Double_t v[3];
                Float_t dcaXY = 0.;
                Float_t dcaZ  = 0.;
              
              	 if(fAODTrackCutBit == 128){
                 dcaXY = trackAOD->DCA();
                 dcaZ  = trackAOD->ZAtDCA();
                }
                else{
                 vertex->GetXYZ(v);
                 trackAOD->GetXYZ(pos);
                 dcaXY  = TMath::Sqrt((pos[0] - v[0])*(pos[0] - v[0]) + (pos[1] - v[1])*(pos[1] - v[1]));
                 dcaZ   = pos[2] - v[2];
                }

                if( fDCAxyCut != -1 && fDCAzCut != -1){
                  if(TMath::Sqrt((dcaXY*dcaXY)/(fDCAxyCut*fDCAxyCut)+(dcaZ*dcaZ)/(fDCAzCut*fDCAzCut)) > 1 ){
                    continue;  // 2D cut
 		  }
                }
		
              Int_t mcGoods = nMCLabelCounter;
              for (Int_t k = 0; k < mcGoods; k++) {
		  Int_t mcLabel = labelMCArray.At(k);
		  
                  if (mcLabel != TMath::Abs(label)) continue;
                  if(mcLabel != label) continue;
                  // if(label > trackAOD->GetLabel()) continue; // MODIFIED 11.01.2017 (take all labels for efficiency)
		  
                  if(fUsePIDFromPDG){
                      AliAODMCParticle *mcTracMatchedWithReco = (AliAODMCParticle*) mcEvent->GetTrack(label);
                      if (!mcTracMatchedWithReco) {
                          AliError(Form("ERROR: Could not receive track %d (match reco - gen)", label));
                          continue;
                      }
		    
                      Int_t pdgcode = mcTracMatchedWithReco->GetPdgCode();
                      if (TMath::Abs(pdgcode) != fPDGCodeWanted) continue;
                  }
		   
		   //acceptance
                   if (fUseY){
                      if(TMath::Abs(trackAOD->Y(fMassParticleOfInterest)) > fMaxEta)
			continue;
                  }
                  else{
		    if(TMath::Abs(trackAOD->Eta()) > fMaxEta)
		      continue;
                  }
                  
                  if((trackAOD->Pt() > fMinPt)&&(trackAOD->Pt() < fMaxPt)) {
		    level[k]  = 2;
                  }
		  else{ continue;}
		  
		   Short_t gCharge = trackAOD->Charge();
                  Double_t phiRad = trackAOD->Phi();
                  Double_t mom = trackAOD->P();

                  AliAODMCParticle *trackAODMCforpdg = (AliAODMCParticle*) mcEvent->GetTrack(label);
                  fHistPdgSurv->Fill(trackAODMCforpdg->GetPdgCode());
                  
		  
		  if(fElectronRejection) {
		    
		    // get the electron nsigma
		    Double_t nSigma = fPIDResponse->NumberOfSigmasTPC(trackAOD,(AliPID::EParticleType)AliPID::kElectron);
		    fHistNSigmaTPCvsPtbeforePID->Fill(trackAOD->Pt(),nSigma);		    
		    
		    // check only for given momentum range
		    if( trackAOD->Pt() > fElectronRejectionMinPt && trackAOD->Pt() < fElectronRejectionMaxPt ){
		      
		      //look only at electron nsigma
		      if(!fElectronOnlyRejection){
			
			//Make the decision based on the n-sigma of electrons
			if(TMath::Abs(nSigma) < fElectronRejectionNSigma) continue;
		      }
		      else{
			
			Double_t nSigmaPions   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackAOD,(AliPID::EParticleType)AliPID::kPion));
			Double_t nSigmaKaons   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackAOD,(AliPID::EParticleType)AliPID::kKaon));
			Double_t nSigmaProtons = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackAOD,(AliPID::EParticleType)AliPID::kProton));
			
			//Make the decision based on the n-sigma of electrons exclusively ( = track not in nsigma region for other species)
			if(TMath::Abs(nSigma) < fElectronRejectionNSigma
			   && nSigmaPions   > fElectronRejectionNSigma
			   && nSigmaKaons   > fElectronRejectionNSigma
			   && nSigmaProtons > fElectronRejectionNSigma ) continue;
		      }
		    }
		    
		    fHistNSigmaTPCvsPtafterPID->Fill(trackAOD->Pt(),nSigma);		    

		  }
                  
	
		   if(fUsePIDstrategy){


                      AliAODPid* pidObj = trackAOD->GetDetPid();
                
                      Double_t probTPC[AliPID::kSPECIES]={0.};
                      Double_t probTPCTOF[AliPID::kSPECIES]={0.};
                
                      fPIDCombined=new AliPIDCombined;
                      fPIDCombined->SetDefaultTPCPriors();
                      Bool_t ParticleFlag = kFALSE;
                      Double_t nSigmaTPCElectrons = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackAOD,(AliPID::EParticleType)AliPID::kElectron));
                      Double_t nSigmaTPCPions   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackAOD,(AliPID::EParticleType)AliPID::kPion));
                      Double_t nSigmaTPCKaons   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackAOD,(AliPID::EParticleType)AliPID::kKaon));
                      Double_t nSigmaTPCProtons = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackAOD,(AliPID::EParticleType)AliPID::kProton));
                
                      if(nSigmaTPCElectrons<3 && nSigmaTPCPions>3 && nSigmaTPCKaons>3 && nSigmaTPCProtons>3 ) continue; //electron rejection
                
                      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC); //firts check only TPC
                      UInt_t detUsed = fPIDCombined->ComputeProbabilities(trackAOD, fPIDResponse, probTPC);
                
                      if (detUsed  == (UInt_t)fPIDCombined->GetDetectorMask()){
                          if (fpartOfInterest==(AliPID::kPion)){ fHistNSigmaTPCvsPtbeforePID->Fill(trackAOD->Pt(),nSigmaTPCPions);}
                          if (fpartOfInterest==(AliPID::kKaon)) fHistNSigmaTPCvsPtbeforePID->Fill(trackAOD->Pt(),nSigmaTPCKaons);
                          if (fpartOfInterest==(AliPID::kProton)) fHistNSigmaTPCvsPtbeforePID->Fill(trackAOD->Pt(),nSigmaTPCProtons);
                          if(mom < 0.5){
                              if (fUsePIDnSigmaComb){
                                  if(fpartOfInterest==(AliPID::kPion) && nSigmaTPCPions<3){
                                      ParticleFlag = kTRUE;
                                  }else if(fpartOfInterest==(AliPID::kKaon) && nSigmaTPCKaons<3){
                                      ParticleFlag = kTRUE;
                                  }else if(fpartOfInterest==(AliPID::kProton) && nSigmaTPCProtons<3){
                                      ParticleFlag = kTRUE;
                                  }else{ParticleFlag = kFALSE;}
                              }else{
                                  if (probTPC[fpartOfInterest] > fBayesPIDThr) ParticleFlag = kTRUE;
                                  if (probTPC[fpartOfInterest] < fBayesPIDThr) ParticleFlag = kFALSE;
                              }
                          }
                      }
                      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF&AliPIDResponse::kDetTPC);
                      detUsed = fPIDCombined->ComputeProbabilities(trackAOD, fPIDResponse, probTPCTOF);
                
                      if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask()){
                          if(!pidObj || pidObj->GetTOFsignal() > 99999)  continue;
                    
                          Double_t nSigmaTOFElectrons = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackAOD,(AliPID::EParticleType)AliPID::kElectron));
                          Double_t nSigmaTOFPions   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackAOD,(AliPID::EParticleType)AliPID::kPion));
                          Double_t nSigmaTOFKaons   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackAOD,(AliPID::EParticleType)AliPID::kKaon));
                          Double_t nSigmaTOFProtons = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackAOD,(AliPID::EParticleType)AliPID::kProton));
                    
                          Double_t combSquaredSigmaElectrons = TMath::Sqrt((nSigmaTPCElectrons*nSigmaTPCElectrons) + (nSigmaTOFElectrons*nSigmaTOFElectrons));
                          Double_t combSquaredSigmaPions = TMath::Sqrt((nSigmaTPCPions*nSigmaTPCPions) + (nSigmaTOFPions*nSigmaTOFPions));
                          Double_t combSquaredSigmaKaons = TMath::Sqrt((nSigmaTPCKaons*nSigmaTPCKaons) + (nSigmaTOFKaons*nSigmaTOFKaons));
                          Double_t combSquaredSigmaProtons = TMath::Sqrt((nSigmaTPCProtons*nSigmaTPCProtons) + (nSigmaTOFProtons*nSigmaTOFProtons));

                          if(mom >= 0.5){
                              if (fUsePIDnSigmaComb){
                                  if(fpartOfInterest==(AliPID::kPion)){
                                      if ((mom <= 2.5) && (TMath::Abs(combSquaredSigmaPions)<3.)){ParticleFlag = kTRUE;}
                                      if ((mom>2.5) && (TMath::Abs(combSquaredSigmaPions)<2.)){ParticleFlag = kTRUE;}
                                      else{ParticleFlag = kFALSE;}
                                  }
                                  if(fpartOfInterest==(AliPID::kKaon)){
                                      if ((mom <= 2.)&&(TMath::Abs(combSquaredSigmaKaons)<2.5)){ParticleFlag = kTRUE;}
                                      if ((mom > 2. && mom<= 3.)&&(TMath::Abs(combSquaredSigmaKaons)<2.)){ParticleFlag = kTRUE;}
                                      if ((mom > 3.)&&(TMath::Abs(combSquaredSigmaKaons)<1.5)){ParticleFlag = kTRUE;}
                                      else{ParticleFlag = kFALSE;}
                                  }
                                  if(fpartOfInterest==(AliPID::kProton)){
                                      if ((mom <= 3.)&&(TMath::Abs(combSquaredSigmaProtons)<3.)){ParticleFlag = kTRUE;}
                                      if ((mom > 3.)&&(mom <= 5.)&&(TMath::Abs(combSquaredSigmaProtons)<2.)){ParticleFlag = kTRUE;}
                                      if ((mom > 5.)&&(TMath::Abs(combSquaredSigmaProtons)<1.5)){ParticleFlag = kTRUE;}
                                      else{ParticleFlag = kFALSE;}
                                  }
                            
                              }//(fUsePIDnSigmaComb)
                              else{
                                  if (probTPCTOF[fpartOfInterest] > fBayesPIDThr) ParticleFlag = kTRUE;
                                  if (probTPCTOF[fpartOfInterest] < fBayesPIDThr) ParticleFlag = kFALSE;
                              }//else
                          }//mom>0.5
                      }
                      if(ParticleFlag && gCharge > 0){
                          if (fUseY)
			    fHistSurvivedEtaPtPhiPlus->Fill(trackAOD->Y(fMassParticleOfInterest),trackAOD->Pt(),phiRad);
                          else
			    fHistSurvivedEtaPtPhiPlus->Fill(trackAOD->Eta(),trackAOD->Pt(),phiRad);
                      }
                      else if(ParticleFlag && gCharge < 0){
			if (fUseY)
			  fHistSurvivedEtaPtPhiMinus->Fill(trackAOD->Y(fMassParticleOfInterest),trackAOD->Pt(),phiRad);
			else
			 fHistSurvivedEtaPtPhiMinus->Fill(trackAOD->Eta(),trackAOD->Pt(),phiRad);
                      }
                  }
		  else {
		    if (gCharge > 0) {
		      if (fUseY)
			fHistSurvivedEtaPtPhiPlus->Fill(trackAOD->Y(fMassParticleOfInterest),trackAOD->Pt(),phiRad);
		      else
			fHistSurvivedEtaPtPhiPlus->Fill(trackAOD->Eta(),trackAOD->Pt(),phiRad);
		    }	
		    else if(gCharge < 0)  {
		      if (fUseY)
			fHistSurvivedEtaPtPhiMinus->Fill(trackAOD->Y(fMassParticleOfInterest),trackAOD->Pt(),phiRad);
		      else
			fHistSurvivedEtaPtPhiMinus->Fill(trackAOD->Eta(),trackAOD->Pt(),phiRad);
		    }	
		  }
		  
	       }//end of mcGoods
	      }//AOD track loop

	      if (fUseRaaGeoCut) delete fESDtrackCuts;
		 
	      labelMCArray.Reset();
	      labelArray.Reset();	       
	      
	      }//Vz cut
	    }//Vy cut
	  }//Vx cut
	}//Vz resolution
      }//number of contributors
    }//valid vertex  
    
    // Here comes the 2 particle analysis
  // loop over all good MC particles
    for (Int_t i = 0; i < nMCLabelCounter ; i++) {
    // control 1D histograms (charge might be different?)
    if(charge[i] > 0){
        if (fUseY){
            if(level[i] > 0) fHistGeneratedEtaPtPlusControl->Fill(vY[i],pt[i]);
            if(level[i] > 1) fHistSurvivedEtaPtPlusControl->Fill(vY[i],pt[i]);
        }
        else if (!fUseY){
            if(level[i] > 0) fHistGeneratedEtaPtPlusControl->Fill(eta[i],pt[i]);
            if(level[i] > 1) fHistSurvivedEtaPtPlusControl->Fill(eta[i],pt[i]);
        }
    }
    else if(charge[i] < 0){
        if (fUseY){
            if(level[i] > 0) fHistGeneratedEtaPtMinusControl->Fill(vY[i],pt[i]);
            if(level[i] > 1) fHistSurvivedEtaPtMinusControl->Fill(vY[i],pt[i]);
        }
        else if (!fUseY){
            if(level[i] > 0) fHistGeneratedEtaPtMinusControl->Fill(eta[i],pt[i]);
            if(level[i] > 1) fHistSurvivedEtaPtMinusControl->Fill(eta[i],pt[i]);
        }
    }
    
    
    for (Int_t j = i+1; j < nMCLabelCounter ; j++) {
      
      if(charge[i] > 0 && charge[j] > 0 ){
	if(level[i] > 0 && level[j] > 0) {
	  if (fUseY)
            fHistGeneratedEtaPtPlusPlus->Fill(TMath::Abs(vY[i]-vY[j]),pt[i]);
	  else
	    fHistGeneratedEtaPtPlusPlus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) <  TMath::Pi()){
	    if (fUseY)
	      fHistGeneratedPhiEtaPlusPlus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(vY[i]-vY[j]));
	    else
	      fHistGeneratedPhiEtaPlusPlus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	  }
	}
	if(level[i] > 1 && level[j] > 1) {
	  if (fUseY)
            fHistSurvivedEtaPtPlusPlus->Fill(TMath::Abs(vY[i]-vY[j]),pt[i]);
	  else
            fHistSurvivedEtaPtPlusPlus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) < TMath::Pi()){
	    if (fUseY)
	      fHistSurvivedPhiEtaPlusPlus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(vY[i]-vY[j]));
	    else
	      fHistSurvivedPhiEtaPlusPlus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
        }
	}
      }
      
      else if(charge[i] < 0 && charge[j] < 0 ){
	if(level[i] > 0 && level[j] > 0) {
	  if (fUseY)
            fHistGeneratedEtaPtMinusMinus->Fill(TMath::Abs(vY[i]-vY[j]),pt[i]);
	  else
            fHistGeneratedEtaPtMinusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) <  TMath::Pi()){
	    if (fUseY)
	      fHistGeneratedPhiEtaMinusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(vY[i]-vY[j]));
	    else
	      fHistGeneratedPhiEtaMinusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	  }
	}
	if(level[i] > 2 && level[j] > 1) {
	  if (fUseY)
            fHistSurvivedEtaPtMinusMinus->Fill(TMath::Abs(vY[i]-vY[j]),pt[i]);
	  else
            fHistSurvivedEtaPtMinusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) <  TMath::Pi()){
	    if (fUseY)
	      fHistSurvivedPhiEtaMinusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(vY[i]-vY[j]));
	    else
	      fHistSurvivedPhiEtaMinusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	  }
	}
      }
      
      else if((charge[i] > 0 && charge[j] < 0)||(charge[i] < 0 && charge[j] > 0)){
	if(level[i] > 0 && level[j] > 0) {
	  if (fUseY)
	    fHistGeneratedEtaPtPlusMinus->Fill(TMath::Abs(vY[i]-vY[j]),pt[i]);
	  else
	    fHistGeneratedEtaPtPlusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) <  TMath::Pi()){
	    if (fUseY)
	      fHistGeneratedPhiEtaPlusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(vY[i]-vY[j]));
	    else 
	      fHistGeneratedPhiEtaPlusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	  }
	}
	if(level[i] > 2 && level[j] > 1) {
	  if (fUseY)
            fHistSurvivedEtaPtPlusMinus->Fill(TMath::Abs(vY[i]-vY[j]),pt[i]);
	  else
            fHistSurvivedEtaPtPlusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) <  TMath::Pi()){
	    if (fUseY)
	      fHistSurvivedPhiEtaPlusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(vY[i]-vY[j]));
	    else
	      fHistSurvivedPhiEtaPlusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	  }
	}	
      }
    }
    }
}

//________________________________________________________________________
void AliAnalysisTaskEffContBF::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}

//____________________________________________________________________//
Bool_t AliAnalysisTaskEffContBF::IsLabelUsed(TArrayI labelArray, Int_t label) {
  //Checks if the label is used already
  Bool_t status = kFALSE;
  for(Int_t i = 0; i < labelArray.GetSize(); i++) {
    if(labelArray.At(i) == label)
      status = kTRUE;
  }

  return status;
}
