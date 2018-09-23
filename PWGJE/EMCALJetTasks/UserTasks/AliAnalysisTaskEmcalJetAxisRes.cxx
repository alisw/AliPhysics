// This Task is written based on Jet QG tagging analysis task (Author: D. Caffarri, L. Cunqueiro).
// Special thanks to the authors of the QG tagging analysis task.
//
// Author: R.Hosokawa 
//
//====================================================
// just testing 
// For now, most of features of the Jet QG tagging analysis task are kept. 
// Most of them should be purged in future.
//====================================================

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THnSparse.h>
#include <TTree.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"
#include "TVector2.h"
#include "AliEventCuts.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliEmcalPythiaInfo.h"
#include "TRandom3.h"

#include "AliAODEvent.h"
#include "AliAnalysisTaskEmcalJetAxisRes.h"

#include "AliEventPoolManager.h"
#include "AliBasicParticle.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEmcalJetAxisRes)

//________________________________________________________________________
AliAnalysisTaskEmcalJetAxisRes::AliAnalysisTaskEmcalJetAxisRes() : 
AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetAxisRes", kTRUE),
  fContainer(0),
  fMinFractionShared(0.5),
  fJetShapeType(kDetEmbPartPythia),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fSelectedShapes(0),
  fminpTTrig(20.),
  fmaxpTTrig(50.),
  fangWindowRecoil(0.6),
  fSemigoodCorrect(0),
  fHolePos(0),
  fHoleWidth(0), 
  fCentSelectOn(kTRUE),
  fCentMin(30),
  fCentMax(50),
  fOneConstSelectOn(kFALSE),
  fDerivSubtrOrder(0),
  fFlowQnVectorMgr(0x0),
  fEventCuts(),
  fIsEventSelDPG(kTRUE),
  fMinJetPtPart(10),
  fMinJetPtDet(0),
  fMinPtForJetRecalc(4.),
  fPoolMgr(),
  fEP(999),
  fEPName(""),
  fDivEvPlane(),
  fDivPt(),
  fh2ResponseUW(0x0),
  fh2ResponseW(0x0), 
  fPhiJetCorr6(0x0), 
  fPhiJetCorr7(0x0),
  fEtaJetCorr6(0x0),
  fEtaJetCorr7(0x0),
  fPtJetCorr(0x0),
  fPtJet(0x0),
  fhpTjetpT(0x0),
  fhPt(0x0),
  fhPhi(0x0),
  fHLundIterative(0x0),
  fNbOfConstvspT(0x0),
  fh2JetdEtadPhiEmbVsDet             (),
  fh2JetdEtadPhiEmbVsPart            (),
  fh2JetdEtadPhiDetVsPart            (),
  fh2JetdEtadPhiRecalcEmbVsRecalcDet (),
  fh2JetdEtadPhiRecalcEmbVsRecalcPart(),
  fh2JetdEtadPhiRecalcDetVsRecalcPart(),
  fh2JetdEtadPhiRecalcEmbVsDet       (),
  fh2JetdEtadPhiRecalcEmbVsPart      (),
  fh2JetdEtadPhiRecalcDetVsPart      (),
  fh2JetHadCorr(),
  fh2JetHadCorrNear(),
  fh2JetHadCorrNearSideBand(),
  fh2JetHadCorrMix(),
  fh2JetHadCorrMixNear(),
  fh2JetHadCorrMixNearSideBand(),
  fTreeObservableTagging(0)

{
  for(Int_t i=0;i<17;i++){
    fShapesVar[i]=0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetAxisRes::AliAnalysisTaskEmcalJetAxisRes(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fContainer(0),
  fMinFractionShared(0.5),
  fJetShapeType(kDetEmbPartPythia),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fPtThreshold(20.),
  fRMatching(0.2),
  fSelectedShapes(0), 
  fminpTTrig(20.),
  fmaxpTTrig(50.),
  fangWindowRecoil(0.6),
  fSemigoodCorrect(0),
  fHolePos(0),
  fHoleWidth(0), 
  fCentSelectOn(kTRUE),
  fCentMin(30),
  fCentMax(50),
  fOneConstSelectOn(kFALSE),
  fDerivSubtrOrder(0),
  fFlowQnVectorMgr(0x0),
  fEventCuts(),
  fMinJetPtPart(10.),
  fMinJetPtDet (0.),
  fIsEventSelDPG(kTRUE),
  fMinPtForJetRecalc(4.),
  fPoolMgr(),
  fEP(999),
  fEPName("VZEROC"),
  fDivEvPlane(),
  fDivPt(),
  fh2ResponseUW(0x0),
  fh2ResponseW(0x0),
  fPhiJetCorr6(0x0), 
  fPhiJetCorr7(0x0),
  fEtaJetCorr6(0x0),
  fEtaJetCorr7(0x0),
  fPtJetCorr(0x0),
  fPtJet(0x0),
  fhpTjetpT(0x0),
  fhPt(0x0),
  fhPhi(0x0),
  fHLundIterative(0x0),
  fNbOfConstvspT(0x0),
  fh2JetdEtadPhiEmbVsDet             (),
  fh2JetdEtadPhiEmbVsPart            (),
  fh2JetdEtadPhiDetVsPart            (),
  fh2JetdEtadPhiRecalcEmbVsRecalcDet (),
  fh2JetdEtadPhiRecalcEmbVsRecalcPart(),
  fh2JetdEtadPhiRecalcDetVsRecalcPart(),
  fh2JetdEtadPhiRecalcEmbVsDet       (),
  fh2JetdEtadPhiRecalcEmbVsPart      (),
  fh2JetdEtadPhiRecalcDetVsPart      (),
  fh2JetHadCorr(),
  fh2JetHadCorrNear(),
  fh2JetHadCorrNearSideBand(),
  fh2JetHadCorrMix(),
  fh2JetHadCorrMixNear(),
  fh2JetHadCorrMixNearSideBand(),
  fTreeObservableTagging(0)
  
{
  // Standard constructor.
  for(Int_t i=0;i<17;i++){
    fShapesVar[i]=0;
  }
  SetMakeGeneralHistograms(kTRUE);
  
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetAxisRes::~AliAnalysisTaskEmcalJetAxisRes()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetAxisRes::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  for(Int_t iep=0; iep<9; iep++)
    fDivEvPlane[iep] = (-TMath::Pi()/2.)+(iep*TMath::Pi()/8.);
  Double_t DivPt[6] = {0.15,0.7,1,2,3,4};
  for(Int_t ipt=0; ipt<6; ipt++)
    fDivPt[ipt] = DivPt[ipt];
 
  fh2ResponseUW= new TH2F("fh2ResponseUW", "fh2ResponseUW", 100, 0, 200,  100, 0, 200); 
  fOutput->Add(fh2ResponseUW);
  fh2ResponseW= new TH2F("fh2ResponseW", "fh2ResponseW", 100, 0, 200,  100, 0, 200);
  fOutput->Add(fh2ResponseW);
  fPhiJetCorr6= new TH2F("fPhiJetCorr6", "fPhiJetCorr6", 50, 0, 2*TMath::Pi(), 50, 0, 2*TMath::Pi());
  fOutput->Add(fPhiJetCorr6);
  fEtaJetCorr6= new TH2F("fEtaJetCorr6", "fEtaJetCorr6", 50, -1.5, 1.5, 50, -1.5, 1.5);
  fOutput->Add(fEtaJetCorr6);
  
  fPhiJetCorr7= new TH2F("fPhiJetCorr7", "fPhiJetCorr7", 50, 0, 2*TMath::Pi(), 50, 0, 2*TMath::Pi());
  fOutput->Add(fPhiJetCorr7);
  fEtaJetCorr7= new TH2F("fEtaJetCorr7", "fEtaJetCorr7", 50, -1.5, 1.5, 50, -1.5, 1.5);
  fOutput->Add(fEtaJetCorr7);
  
  fPtJetCorr= new TH2F("fPtJetCorr", "fPtJetCorr", 100, 0, 200,  100, 0, 200);
  fOutput->Add(fPtJetCorr);
  fPtJet= new TH1F("fPtJet", "fPtJet", 100, 0, 200);
  fOutput->Add(fPtJet);
  
  fhpTjetpT= new TH2F("fhpTjetpT", "fhpTjetpT", 200, 0, 200,  200, 0, 200);
  fOutput->Add(fhpTjetpT);
  fhPt= new TH1F("fhPt", "fhPt", 200, 0, 200);
  fOutput->Add(fhPt);
  fhPhi= new TH1F("fhPhi", "fhPhi", 100, -TMath::Pi(), TMath::Pi());
  fOutput->Add(fhPhi);

  for(Int_t iep=0; iep<9; iep++) {
    fh2JetdEtadPhiEmbVsDet             [iep] = new TH2D(Form("fh2JetdEtadPhiEmbVsDet_ep%d", iep),Form("fh2JetdEtadPhiEmbVsDet_ep%d ; #Delta#eta(#eta_{Emb}-#eta{Det}) ; #Delta#phi(#phi_{Emb}-#phi{Det})",  iep), 120, -0.529375, 0.520625, 120, -0.529375, 0.520625);
    fh2JetdEtadPhiEmbVsPart            [iep] = new TH2D(Form("fh2JetdEtadPhiEmbVsPart_ep%d", iep),Form("fh2JetdEtadPhiEmbVsPart_ep%d ; #Delta#eta(#eta_{Emb}-#eta{Part}); #Delta#phi(#phi_{Emb}-#phi{Part})", iep), 120, -0.529375, 0.520625, 120, -0.529375, 0.520625);
    
    fh2JetdEtadPhiRecalcEmbVsRecalcDet [iep] = new TH2D(Form("fh2JetdEtadPhiRecalcEmbVsRecalcDet_ep%d ", iep),Form("fh2JetdEtadPhiRecalcEmbVsRecalcDet_ep%d ; #Delta#eta(#eta_{Emb}-#eta{Det}) ; #Delta#phi(#phi_{Emb}-#phi{Det})", iep), 120, -0.529375, 0.520625, 120, -0.529375, 0.520625);
    fh2JetdEtadPhiRecalcEmbVsRecalcPart[iep] = new TH2D(Form("fh2JetdEtadPhiRecalcEmbVsRecalcPart_ep%d", iep),Form("fh2JetdEtadPhiRecalcEmbVsRecalcPart_ep%d; #Delta#eta(#eta_{Emb}-#eta{Part}); #Delta#phi(#phi_{Emb}-#phi{Part})", iep), 120, -0.529375, 0.520625, 120, -0.529375, 0.520625);
    
    fh2JetdEtadPhiRecalcEmbVsDet       [iep] = new TH2D(Form("fh2JetdEtadPhiRecalcEmbVsDet_ep%d", iep),Form("fh2JetdEtadPhiRecalcEmbVsDet_ep%d  ; #Delta#eta(#eta_{Emb}-#eta{Det}) ; #Delta#phi(#phi_{Emb}-#phi{Det})",  iep), 120, -0.529375, 0.520625, 120, -0.529375, 0.520625);
    fh2JetdEtadPhiRecalcEmbVsPart      [iep] = new TH2D(Form("fh2JetdEtadPhiRecalcEmbVsPart_ep%d", iep),Form("fh2JetdEtadPhiRecalcEmbVsPart_ep%d ; #Delta#eta(#eta_{Emb}-#eta{Part}); #Delta#phi(#phi_{Emb}-#phi{Part})", iep), 120, -0.529375, 0.520625, 120, -0.529375, 0.520625);
    

    fOutput->Add(fh2JetdEtadPhiEmbVsDet             [iep]);
    fOutput->Add(fh2JetdEtadPhiEmbVsPart            [iep]);
    fOutput->Add(fh2JetdEtadPhiRecalcEmbVsRecalcDet [iep]);
    fOutput->Add(fh2JetdEtadPhiRecalcEmbVsRecalcPart[iep]);
    fOutput->Add(fh2JetdEtadPhiRecalcEmbVsDet       [iep]);
    fOutput->Add(fh2JetdEtadPhiRecalcEmbVsPart      [iep]);
  }

  for(Int_t iep=0; iep<9; iep++) {
    for(Int_t ipt=0; ipt<6; ipt++) {
      fh2JetHadCorr[iep][ipt]                = new TH2D(Form("fh2JetHadCorr_ep%d_pt%d", iep, ipt), Form("fh2JetHadCorr_ep%d_pt%d; #Delta #eta; #Delta #phi", iep, ipt), 100, -1.6, 1.6, 100, -TMath::Pi(), TMath::Pi()); 
      fh2JetHadCorrNear[iep][ipt]            = new TH2D(Form("fh2JetHadCorrNear_ep%d_pt%d", iep, ipt), Form("fh2JetHadCorrNear_ep%d_pt%d; #Delta #eta; #Delta #phi", iep, ipt), 160, -0.8, 0.8, 160, -0.8, 0.8); 
      fh2JetHadCorrNearSideBand[iep][ipt]    = new TH2D(Form("fh2JetHadCorrNearSideBand_ep%d_pt%d", iep, ipt), Form("fh2JetHadCorrNearSideBand_ep%d_pt%d; #Delta #eta; #Delta #phi", iep, ipt), 160, -0.8, 0.8, 160, -0.8, 0.8);
      fh2JetHadCorrMix[iep][ipt]             = new TH2D(Form("fh2JetHadCorrMix_ep%d_pt%d", iep, ipt), Form("fh2JetHadCorrMix_ep%d_pt%d; #Delta #eta; #Delta #phi", iep, ipt), 100, -1.6, 1.6, 100, -TMath::Pi(), TMath::Pi()); 
      fh2JetHadCorrMixNear[iep][ipt]         = new TH2D(Form("fh2JetHadCorrMixNear_ep%d_pt%d", iep, ipt), Form("fh2JetHadCorrMixNear_ep%d_pt%d; #Delta #eta; #Delta #phi", iep, ipt), 160, -0.8, 0.8, 160, -0.8, 0.8); 
      fh2JetHadCorrMixNearSideBand[iep][ipt] = new TH2D(Form("fh2JetHadCorrMixNearSideBand_ep%d_pt%d", iep, ipt), Form("fh2JetHadCorrMixNearSideBand_ep%d_%d; #Delta #eta; #Delta #phi", iep, ipt), 160, -0.8, 0.8, 160, -0.8, 0.8);
    }
  }

  fh2JetdEtadPhiDetVsPart             = new TH2D("fh2JetdEtadPhiDetVsPart", "fh2JetdEtadPhiDetVsPart; #Delta#eta(#eta_{Det}-#eta{Part}); #Delta#phi(#phi_{Det}-#phi{Part})", 120, -0.529375, 0.520625, 120, -0.529375, 0.520625);
  fh2JetdEtadPhiRecalcDetVsRecalcPart = new TH2D("fh2JetdEtadPhiRecalcDetVsRecalcPart", "fh2JetdEtadPhiRecalcDetVsRecalcPart; #Delta#eta(#eta_{Det}-#eta{Part}); #Delta#phi(#phi_{Det}-#phi{Part})", 120, -0.529375, 0.520625, 120, -0.529375, 0.520625);
  fh2JetdEtadPhiRecalcDetVsPart       = new TH2D("fh2JetdEtadPhiRecalcDetVsPart", "fh2JetdEtadPhiRecalcDetVsPart; #Delta#eta(#eta_{Det}-#eta{Part}); #Delta#phi(#phi_{Det}-#phi{Part})", 120, -0.529375, 0.520625, 120, -0.529375, 0.520625);

  fOutput->Add(fh2JetdEtadPhiDetVsPart            );
  fOutput->Add(fh2JetdEtadPhiRecalcDetVsRecalcPart);    
  fOutput->Add(fh2JetdEtadPhiRecalcDetVsPart      );
  //log(1/theta),log(z*theta),jetpT,algo// 
  const Int_t dimSpec   = 5;
  const Int_t nBinsSpec[5]     = {50,50,10,3,10};
  const Double_t lowBinSpec[5] = {0.0,-10,  0,0,0};
  const Double_t hiBinSpec[5]  = {5.0,  0,200,3,10};
  fHLundIterative = new THnSparseF("fHLundIterative",
				   "LundIterativePlot [log(1/theta),log(z*theta),pTjet,algo]",
				   dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);
  fOutput->Add(fHLundIterative);  
  
  fNbOfConstvspT=new TH2F("fNbOfConstvspT", "fNbOfConstvspT", 100, 0, 100, 200, 0, 200);
  fOutput->Add(fNbOfConstvspT);
  
  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn)hn->Sumw2();
  }

 
  TH1::AddDirectory(oldStatus);
  const Int_t nVar = 12;
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeObservableTagging = new TTree(nameoutput, nameoutput);
  
 
  TString *fShapesVarNames = new TString [nVar];

  fShapesVarNames[0] = "partonCode"; 
  fShapesVarNames[1] = "ptJet"; 
  fShapesVarNames[2] = "ptDJet"; 
  fShapesVarNames[3] = "phiJet";
  // fShapesVarNames[4] = "nbOfConst";
  fShapesVarNames[4] = "angularity";
  //fShapesVarNames[5] = "circularity";
  fShapesVarNames[5] = "lesub";
  //fShapesVarNames[7] = "coronna";

  fShapesVarNames[6] = "ptJetMatch"; 
  fShapesVarNames[7] = "ptDJetMatch"; 
  fShapesVarNames[8] = "phiJetMatch";
  // fShapesVarNames[12] = "nbOfConstMatch";
  fShapesVarNames[9] = "angularityMatch";
  //fShapesVarNames[12] = "circularityMatch";
  fShapesVarNames[10] = "lesubMatch";
  //fShapesVarNames[14] = "coronnaMatch";
  fShapesVarNames[11]="weightPythia";
  //fShapesVarNames[14]="ntrksEvt";
  //fShapesVarNames[16]="rhoVal";
  //fShapesVarNames[17]="rhoMassVal";
  //fShapesVarNames[12]="ptUnsub";

  for(Int_t ivar=0; ivar < nVar; ivar++){
    cout<<"looping over variables"<<endl;
    fTreeObservableTagging->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/F", fShapesVarNames[ivar].Data()));}
 
  PostData(1,fOutput);
  PostData(2,fTreeObservableTagging);

  delete [] fShapesVarNames;

  if(fIsEventSelDPG) {
    fEventCuts.SetManualMode();
    fEventCuts.SetupLHC15o();
  }

  if(fEPName!="") {
    AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections *>
      (AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
    if (flowQnVectorTask != 0x0) {
      fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
    }
    else {
      AliFatal("Pointer to object of Flow Qn vector corrections framework is needed but the object is not present. ABORTING!!!");
    }
  }

  // Event Mixing
  Int_t poolSize = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
  // ZVertex
  Int_t nZVertexBins = 8;//16;//
  Double_t* zVertexBins = GenerateFixedBinArray(nZVertexBins, -8, 8);
  // Event activity (centrality of multiplicity)
  Int_t nEventActivityBins = 8;
  Double_t* eventActivityBins = 0;
  // +1 to accomodate the fact that we define bins rather than array entries.
  const Int_t kMixedEventMulitplictyBins=8;
  Double_t multiplicityBins[kMixedEventMulitplictyBins+1] = {0., 4., 9., 15., 25., 35., 55., 100., 500.};

  // Cannot use GetBeamType() since it is not available until UserExec()
  if (fForceBeamType != AliAnalysisTaskEmcal::kpp ) {   //all besides pp
    // Event Activity is centrality in AA, pA
    //nEventActivityBins = fNCentBinsMixedEvent;
    //nEventActivityBins = 10;
    nEventActivityBins = 20;
    eventActivityBins = GenerateFixedBinArray(nEventActivityBins, 0, 100);
  }
  else if (fForceBeamType == AliAnalysisTaskEmcal::kpp) { //for pp only
    // Event Activity is multiplicity in pp
    eventActivityBins = multiplicityBins;
  }
  Int_t fNMixingTracks = 50000;
  fPoolMgr = new AliEventPoolManager(poolSize, fNMixingTracks, nEventActivityBins, eventActivityBins, nZVertexBins, zVertexBins);
  fPoolMgr->SetTargetValues(fNMixingTracks, 1., 5);

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetAxisRes::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  if(fEPName!="") {
    fEP=999;
    const AliQnCorrectionsQnVector *QnVector = 0x0;
    if(fEPName=="VZEROC" || fEPName=="VZEROA")
      QnVector  = fFlowQnVectorMgr->GetDetectorQnVector(fEPName.Data(),"align","align");
    else
      QnVector  = fFlowQnVectorMgr->GetDetectorQnVector(fEPName.Data(),"latest","latest");

    if(QnVector)
      fEP = QnVector->EventPlane(2);
  }

  if(fIsEventSelDPG) {
    if(fRunNumber==244918 || fRunNumber==244975 || fRunNumber==244980 || fRunNumber==244982 || fRunNumber==244983 || fRunNumber==245064 || fRunNumber==245066 || fRunNumber==245068 || fRunNumber==246390 || fRunNumber==246391 || fRunNumber==246392 || fRunNumber==245148) {
      fEventCuts.fUseEstimatorsCorrelationCut = false;
      fEventCuts.fUseVariablesCorrelationCuts = false;
    }
    else if(fRunNumber==246428 || fRunNumber==246980) {
      fEventCuts.fUseEstimatorsCorrelationCut = true;
      fEventCuts.fUseVariablesCorrelationCuts = false;
    }
    else {
      fEventCuts.fUseEstimatorsCorrelationCut = true;
      fEventCuts.fUseVariablesCorrelationCuts = true;
    }
    //fEventCuts.fMinVtz=-1000.f; //For test                                                                                                                                                              
    //fEventCuts.fMaxVtz=1000.f;  //For test                                                                                                                                                              
    //fEventCuts.fRequireTrackVertex = false;                                                                                                                                                          
    if(fEventCuts.AcceptEvent((AliVEvent*)InputEvent()))
      return kTRUE;
    else
      return kFALSE;
  }
  else
    return kTRUE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetAxisRes::FillHistograms()
{
  // Fill histograms.
  //cout<<"base container"<<endl;
  AliEmcalJet* jet1 = NULL;
  AliJetContainer *jetCont = GetJetContainer(0);
  AliEmcalJet* jetrecalc=0x0;


  Float_t kWeight=1;
  if (fCentSelectOn)
    if ((fCent>fCentMax) || (fCent<fCentMin)) return 0;
  
  AliAODTrack *triggerHadron = 0x0;
  
  if (fJetSelection == kRecoil) {
    //Printf("Recoil jets!!!, fminpTTrig = %f, fmaxpTTrig = %f", fminpTTrig, fmaxpTTrig);
    Int_t triggerHadronLabel = SelectTrigger(fminpTTrig, fmaxpTTrig);
     
    
    if (triggerHadronLabel==-99999) {
      //Printf ("Trigger Hadron not found, return");
      return 0;}

    AliTrackContainer *PartCont =NULL;
    AliParticleContainer *PartContMC=NULL;

    if (fJetShapeSub==kConstSub){
      if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly) PartContMC = GetParticleContainer(1);
      else PartCont = GetTrackContainer(1);
    }
    else{
      if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly) PartContMC = GetParticleContainer(0);
      else PartCont = GetTrackContainer(0);
    }
    TClonesArray *TrackArray = NULL;
    TClonesArray *TrackArrayMC = NULL;
    if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly) TrackArrayMC = PartContMC->GetArray();
    else TrackArray = PartCont->GetArray();    
    if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly) triggerHadron = static_cast<AliAODTrack*>(TrackArrayMC->At(triggerHadronLabel));
    else triggerHadron = static_cast<AliAODTrack*>(TrackArray->At(triggerHadronLabel));




    
    if (!triggerHadron) {
      //Printf("No Trigger hadron with the found label!!");
      return 0;
    }

    if(fSemigoodCorrect){
      Double_t disthole=RelativePhi(triggerHadron->Phi(),fHolePos);
      if(TMath::Abs(disthole)+fHoleWidth>TMath::Pi()-fangWindowRecoil){
        return 0;}
    }
   
    fhPt->Fill(triggerHadron->Pt());

  }

  
  
  AliParticleContainer *partContAn = GetParticleContainer(0);
  TClonesArray *trackArrayAn = partContAn->GetArray();
  Int_t ntracksEvt = trackArrayAn->GetEntriesFast();


  // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
  TObjArray* tracksClone = 0;

  // event mixing
  
  // 1. First get an event pool corresponding in mult (cent) and
  //    zvertex to the current event. Once initialized, the pool
  //    should contain nMix (reduced) events. This routine does not
  //    pre-scan the chain. The first several events of every chain
  //    will be skipped until the needed pools are filled to the
  //    specified depth. If the pool categories are not too rare, this
  //    should not be a problem. If they are rare, you could lose
  //    statistics.
  
  // 2. Collect the whole pool's content of tracks into one TObjArray
  //    (bgTracks), which is effectively a single background super-event.
  
  // 3. The reduced and bgTracks arrays must both be passed into
  //    FillCorrelations(). Also nMix should be passed in, so a weight
  //    of 1./nMix can be applied.
  
  AliEventPool *pool = 0;
  if (fForceBeamType == kAA || fForceBeamType == kpA) {//everything but pp
    //pool = fPoolMgr->GetEventPool(fCent, zVertex);
    pool = fPoolMgr->GetEventPool(fCent, fVertex[2]);
  }
  else if (fForceBeamType == kpp) {//pp only
    //pool = fPoolMgr->GetEventPool(static_cast<Double_t>(tracks->GetNTracks()), zVertex);
    pool = fPoolMgr->GetEventPool(static_cast<Double_t>(partContAn->GetNAcceptedParticles()), fVertex[2]);
  }
  
  if (!pool){
    if (fForceBeamType == kAA || fForceBeamType == kpA) AliFatal(Form("No pool found for centrality = %f, zVertex = %f", fCent, fVertex[2]));
    else if (fForceBeamType == kpp) AliFatal(Form("No pool found for ntracks_pp = %d, zVertex = %f", partContAn->GetNAcceptedParticles(), fVertex[2]));
    //return kTRUE;
  }
  // The number of events in the pool
  Int_t nMix = pool->GetCurrentNEvents();

  
  Float_t rhoVal=0, rhoMassVal = 0.;
  if(jetCont) {

    jetCont->ResetCurrentID();
    if ((fJetShapeSub==kConstSub) || (fJetShapeSub==kDerivSub)){
      //rho                                                                                                   
      AliRhoParameter* rhoParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject("RhoSparseR020"));
      if (!rhoParam) {
	Printf("%s: Could not retrieve rho %s (some histograms will be filled with zero)!", GetName(), jetCont->GetRhoName().Data());
      } else rhoVal = rhoParam->GetVal();
      //rhom                                                                                                                                                          
      AliRhoParameter* rhomParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject("RhoMassSparseR020"));
      if (!rhomParam) {
	Printf("%s: Could not retrieve rho_m %s (some histograms will be filled with zero)!", GetName(), jetCont->GetRhoMassName().Data());	
      } else rhoMassVal = rhomParam->GetVal();
    }
    
    while((jet1 = jetCont->GetNextAcceptJet())) {
      if (!jet1) continue;
      AliEmcalJet* jet2 = 0x0;
      AliEmcalJet* jet3 = 0x0;
      fPtJet->Fill(jet1->Pt());
      AliEmcalJet *jetUS = NULL;
      Int_t ifound=0, jfound=0;
      Int_t ilab=-1, jlab=-1;
      
      if(fSemigoodCorrect && (fJetSelection != kRecoil)){
	Double_t disthole=RelativePhi(jet1->Phi(),fHolePos);
	if(TMath::Abs(disthole)<fHoleWidth){
	  continue;}
      } 
      
      Float_t dphiRecoil = 0.;
      if (fJetSelection == kRecoil){
        dphiRecoil = RelativePhi(triggerHadron->Phi(), jet1->Phi());
        if (TMath::Abs(dphiRecoil) < (TMath::Pi() - fangWindowRecoil)) {
	  // Printf("Recoil jets back to back not found! continuing");
          continue;
        }
        
        fhpTjetpT->Fill(triggerHadron->Pt(), jet1->Pt());
        //Printf(" ************ FILLING HISTOS****** shapeSub = %d, triggerHadron = %f, jet1 = %f", fJetShapeSub, triggerHadron->Pt(), jet1->Pt());
        fhPhi->Fill(RelativePhi(triggerHadron->Phi(), jet1->Phi()));
        
      }
      
      
      fShapesVar[0] = 0.;
      if(fJetShapeType == kDetEmbPartPythia){
        AliJetContainer *jetContTrue = GetJetContainer(1);
        AliJetContainer *jetContUS = GetJetContainer(2);
	AliEmcalJet* jetRecalcDetPythia=0x0;
	AliEmcalJet* jetRecalcPartPythia=0x0;

        if(fJetShapeSub==kConstSub){
	  for(Int_t i = 0; i<jetContUS->GetNJets(); i++) {
	    jetUS = jetContUS->GetJet(i);
            if(jetUS->GetLabel()==jet1->GetLabel()) {
              ifound++;
              if(ifound==1) ilab = i;
            }
          }
          if(ilab==-1) continue;
          jetUS=jetContUS->GetJet(ilab);
          jet2=jetUS->ClosestJet();
        }
        
        if(!(fJetShapeSub==kConstSub)) jet2 = jet1->ClosestJet();
        if (!jet2) {
          Printf("jet2 does not exist, returning");
          continue;
        }
        
        AliJetContainer *jetContPart=GetJetContainer(3);
        jet3=jet2->ClosestJet();
        
        if(!jet3){
          Printf("jet3 does not exist, returning");
          continue;
        }
        cout<<"jet 3 exists"<<jet3->Pt()<<endl;
	if(jet1->Pt()-GetRhoVal(0)*jet1->Area()<fPtThreshold) continue;
	if(jet2->Pt()<fMinJetPtDet)  continue;
	if(jet3->Pt()<fMinJetPtPart) continue;

	if(fMinPtForJetRecalc>0.) { // jet axis recalculation (pt center of constituents)
	  Double_t MSumRecalc=0.;   
	  Double_t PtSumRecalc=0.;  
	  Double_t PhiRecalc=0.;    
	  Double_t EtaRecalc=0.;    

	  // jet1(embeded jets)
	  for(Int_t icons=0; icons<jet1->GetNumberOfTracks(); icons++) {
	    AliVParticle *vp = static_cast<AliVParticle*>(jet1->TrackAt(icons, jetCont->GetParticleContainer()->GetArray()));
	    if(!vp) continue;
	    if(vp->Pt()<fMinPtForJetRecalc) continue;
	    MSumRecalc  += vp->M();
	    PtSumRecalc += vp->Pt();
	    PhiRecalc   += TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(vp->Phi()) - TVector2::Phi_0_2pi(jet1->Phi()))*vp->Pt();
	    EtaRecalc   += (vp->Eta()-jet1->Eta())*vp->Pt();
	  }
	  if(PtSumRecalc>0) {
	    PhiRecalc /= PtSumRecalc;
	    EtaRecalc /= PtSumRecalc;
	    PhiRecalc = jet1->Phi()+PhiRecalc;
	    EtaRecalc = jet1->Eta()+EtaRecalc;
	    jetrecalc  = new AliEmcalJet(PtSumRecalc, EtaRecalc, PhiRecalc, MSumRecalc);
	  }
	  

	  MSumRecalc=0.;   
	  PtSumRecalc=0.;  
	  PhiRecalc=0.;    
	  EtaRecalc=0.;    

	  // jet2 (detector level pythia jets)
	  for(Int_t icons=0; icons<jet2->GetNumberOfTracks(); icons++) {
	    AliVParticle *vp = static_cast<AliVParticle*>(jet2->TrackAt(icons, jetContTrue->GetParticleContainer()->GetArray()));
	    if(!vp) continue;
	    if(vp->Pt()<fMinPtForJetRecalc) continue;
	    MSumRecalc  += vp->M();
	    PtSumRecalc += vp->Pt();
	    PhiRecalc   += TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(vp->Phi()) - TVector2::Phi_0_2pi(jet2->Phi()))*vp->Pt();
	    EtaRecalc   += (vp->Eta()-jet2->Eta())*vp->Pt();
	  }
	  if(PtSumRecalc>0) {
	    PhiRecalc /= PtSumRecalc;
	    EtaRecalc /= PtSumRecalc;
	    PhiRecalc = jet2->Phi()+PhiRecalc;
	    EtaRecalc = jet2->Eta()+EtaRecalc;
	    jetRecalcDetPythia  = new AliEmcalJet(PtSumRecalc, EtaRecalc, PhiRecalc, MSumRecalc);
	  }

	  MSumRecalc=0.;   
	  PtSumRecalc=0.;  
	  PhiRecalc=0.;    
	  EtaRecalc=0.;

	  // jet3 (particle level pythia jets)
	  for(Int_t icons=0; icons<jet3->GetNumberOfTracks(); icons++) {
	    AliVParticle *vp = static_cast<AliVParticle*>(jet3->TrackAt(icons, jetContPart->GetParticleContainer()->GetArray()));
	    if(!vp) continue;
	    if(vp->Pt()<fMinPtForJetRecalc) continue;
	    MSumRecalc  += vp->M();
	    PtSumRecalc += vp->Pt();
	    PhiRecalc   += TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(vp->Phi()) - TVector2::Phi_0_2pi(jet3->Phi()))*vp->Pt();
	    EtaRecalc   += (vp->Eta()-jet3->Eta())*vp->Pt();
	  }
	  if(PtSumRecalc>0) {
	    PhiRecalc /= PtSumRecalc;
	    EtaRecalc /= PtSumRecalc;
	    PhiRecalc = jet3->Phi()+PhiRecalc;
	    EtaRecalc = jet3->Eta()+EtaRecalc;
	    jetRecalcPartPythia  = new AliEmcalJet(PtSumRecalc, EtaRecalc, PhiRecalc, MSumRecalc);
	  }
	}        
        
        fh2ResponseUW->Fill(jet1->Pt(),jet2->Pt());
        
        Double_t fraction=0;
        if(!(fJetShapeSub==kConstSub))  fraction = jetCont->GetFractionSharedPt(jet1);
        if(fJetShapeSub==kConstSub) fraction = jetContUS->GetFractionSharedPt(jetUS);
        //if (fraction > 0.1) cout<<"***** hey a jet matched with fraction"<<fraction<<"  "<<jet1->Pt()<<" "<<jet2->Pt()<<" "<<fCent<<endl;
        if(fraction<fMinFractionShared){ 
	  //InputEvent()->Print();
	  if(jetrecalc)
	    jetrecalc->Delete();
	  if(jetRecalcDetPythia) 
	    jetRecalcDetPythia->Delete();
	  if(jetRecalcPartPythia)
	    jetRecalcPartPythia->Delete();
	  continue;
	}
	Double_t Jet1AxisPhiWrtEvPlane=-999;
	Double_t Jet2AxisPhiWrtEvPlane=-999;
	Double_t Jet3AxisPhiWrtEvPlane=-999;
	Int_t    EPBinJet1=-1;
	Int_t    EPBinJet2=-1;
	Int_t    EPBinJet3=-1;

	Double_t Jet1RecalcAxisPhiWrtEvPlane=-999;
	Double_t Jet2RecalcAxisPhiWrtEvPlane=-999;
	Double_t Jet3RecalcAxisPhiWrtEvPlane=-999;
	Int_t    EPBinRecalcJet1=-1;
	Int_t    EPBinRecalcJet2=-1;
	Int_t    EPBinRecalcJet3=-1;

	if(TMath::Abs(TVector2::Phi_mpi_pi(jet1->Phi()) - fEP)>TMath::Pi()/2.)
	  Jet1AxisPhiWrtEvPlane = TVector2::Phi_mpi_pi(jet1->Phi()-TMath::Pi()) -fEP;
	else
	  Jet1AxisPhiWrtEvPlane = TVector2::Phi_mpi_pi(jet1->Phi()) - fEP;

	if(TMath::Abs(TVector2::Phi_mpi_pi(jet2->Phi()) - fEP)>TMath::Pi()/2.)
	  Jet2AxisPhiWrtEvPlane = TVector2::Phi_mpi_pi(jet2->Phi()-TMath::Pi()) -fEP;
	else
	  Jet2AxisPhiWrtEvPlane = TVector2::Phi_mpi_pi(jet2->Phi()) - fEP;

	if(TMath::Abs(TVector2::Phi_mpi_pi(jet3->Phi()) - fEP)>TMath::Pi()/2.)
	  Jet3AxisPhiWrtEvPlane = TVector2::Phi_mpi_pi(jet3->Phi()-TMath::Pi()) -fEP;
	else
	  Jet3AxisPhiWrtEvPlane = TVector2::Phi_mpi_pi(jet3->Phi()) - fEP;

	if(jetrecalc) {
	  if(TMath::Abs(TVector2::Phi_mpi_pi(jetrecalc->Phi()) - fEP)>TMath::Pi()/2.)
	    Jet1RecalcAxisPhiWrtEvPlane = TVector2::Phi_mpi_pi(jetrecalc->Phi()-TMath::Pi()) -fEP;
	  else
	    Jet1RecalcAxisPhiWrtEvPlane = TVector2::Phi_mpi_pi(jetrecalc->Phi()) - fEP;
	}
	
	if(jetRecalcDetPythia) {
	  if(TMath::Abs(TVector2::Phi_mpi_pi(jetRecalcDetPythia->Phi()) - fEP)>TMath::Pi()/2.)
	    Jet2RecalcAxisPhiWrtEvPlane = TVector2::Phi_mpi_pi(jetRecalcDetPythia->Phi()-TMath::Pi()) -fEP;
	  else
	    Jet2RecalcAxisPhiWrtEvPlane = TVector2::Phi_mpi_pi(jetRecalcDetPythia->Phi()) - fEP;
	}
	
	if(jetRecalcPartPythia) {
	  if(TMath::Abs(TVector2::Phi_mpi_pi(jetRecalcPartPythia->Phi()) - fEP)>TMath::Pi()/2.)
	    Jet3RecalcAxisPhiWrtEvPlane = TVector2::Phi_mpi_pi(jetRecalcPartPythia->Phi()-TMath::Pi()) -fEP;
	  else
	    Jet3RecalcAxisPhiWrtEvPlane = TVector2::Phi_mpi_pi(jetRecalcPartPythia->Phi()) - fEP;
	}

        for(Int_t icl=0; icl<8; icl++) {
	  if(fDivEvPlane[icl]<Jet1AxisPhiWrtEvPlane && Jet1AxisPhiWrtEvPlane<fDivEvPlane[icl+1])
	    EPBinJet1 = icl+1;
	  if(fDivEvPlane[icl]<Jet2AxisPhiWrtEvPlane && Jet2AxisPhiWrtEvPlane<fDivEvPlane[icl+1])
	    EPBinJet2 = icl+1;
          if(fDivEvPlane[icl]<Jet3AxisPhiWrtEvPlane && Jet3AxisPhiWrtEvPlane<fDivEvPlane[icl+1])
	    EPBinJet3 = icl+1;
	  if(jetrecalc) {
	    if(fDivEvPlane[icl]<Jet1RecalcAxisPhiWrtEvPlane && Jet1RecalcAxisPhiWrtEvPlane<fDivEvPlane[icl+1])
	      EPBinRecalcJet1 = icl+1;
	  }
	  if(jetRecalcDetPythia) {
	    if(fDivEvPlane[icl]<Jet2RecalcAxisPhiWrtEvPlane && Jet2RecalcAxisPhiWrtEvPlane<fDivEvPlane[icl+1])
	      EPBinRecalcJet2 = icl+1;
	  }
	  if(jetRecalcPartPythia) {
	    if(fDivEvPlane[icl]<Jet3RecalcAxisPhiWrtEvPlane && Jet3RecalcAxisPhiWrtEvPlane<fDivEvPlane[icl+1])
	      EPBinRecalcJet3 = icl+1;
	  }
	}

	partContAn->ResetCurrentID();
	Int_t PtBin=-1;
	Double_t DeltaEta=999;
	Double_t DeltaPhi=999;
	AliVParticle *part=0x0;
	while((part = partContAn->GetNextAcceptParticle())) {
	  if (!part) continue;

	  PtBin=-1;
	  DeltaEta=999;
	  DeltaPhi=999;	  
	  for(Int_t ipt=0; ipt<5; ipt++) {
	    if(fDivPt[ipt]<=part->Pt() && part->Pt()<fDivPt[ipt+1]) {
	      PtBin = ipt+1;
	      break;
	    }
	  }
	  if(jetrecalc) {
	    DeltaEta = part->Eta() - jetrecalc->Eta();
	    DeltaPhi =TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(part->Phi())-TVector2::Phi_0_2pi(jetrecalc->Phi()));
	  }
	  else if(!(fMinPtForJetRecalc>0.)) {
	    DeltaEta = part->Eta() - jet1->Eta();
	    DeltaPhi =TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(part->Phi())-TVector2::Phi_0_2pi(jet1->Phi()));
	  }

	  if(DeltaEta<999 && DeltaPhi<999) {
	    fh2JetHadCorr[0][0]               ->Fill(DeltaEta, DeltaPhi);
	    fh2JetHadCorrNear[0][0]           ->Fill(DeltaEta, DeltaPhi);
	    if(1.<DeltaEta && DeltaEta<1.5)
	      fh2JetHadCorrNearSideBand[0][0]->Fill(DeltaEta-1., DeltaPhi);
	    if(-1.5<DeltaEta && DeltaEta<-1.)
	      fh2JetHadCorrNearSideBand[0][0]->Fill(DeltaEta+1., DeltaPhi);
	 
	    if(fEP!=999) {   
	      if(jetrecalc) {
		fh2JetHadCorr[EPBinRecalcJet1][0]               ->Fill(DeltaEta, DeltaPhi);
		fh2JetHadCorrNear[EPBinRecalcJet1][0]           ->Fill(DeltaEta, DeltaPhi);
		if(1.<DeltaEta && DeltaEta<1.5)
		  fh2JetHadCorrNearSideBand[EPBinRecalcJet1][0]   ->Fill(DeltaEta-1., DeltaPhi);
		if(-1.5<DeltaEta && DeltaEta<-1.)
		  fh2JetHadCorrNearSideBand[EPBinRecalcJet1][0]   ->Fill(DeltaEta+1., DeltaPhi);
		if(PtBin>0) {
		  fh2JetHadCorr[EPBinRecalcJet1][PtBin]            ->Fill(DeltaEta, DeltaPhi);
		  fh2JetHadCorrNear[EPBinRecalcJet1][PtBin]        ->Fill(DeltaEta, DeltaPhi);
		  if(1.<DeltaEta && DeltaEta<1.5)
		    fh2JetHadCorrNearSideBand[EPBinRecalcJet1][PtBin]->Fill(DeltaEta-1., DeltaPhi);
		  if(-1.5<DeltaEta && DeltaEta<-1.)
		    fh2JetHadCorrNearSideBand[EPBinRecalcJet1][PtBin]->Fill(DeltaEta+1., DeltaPhi);
		}
	      }
	      else if(!(fMinPtForJetRecalc>0.)) {
		fh2JetHadCorr[EPBinJet1][0]               ->Fill(DeltaEta, DeltaPhi);
		fh2JetHadCorrNear[EPBinJet1][0]           ->Fill(DeltaEta, DeltaPhi);
		if(1.<DeltaEta && DeltaEta<1.5)
		  fh2JetHadCorrNearSideBand[EPBinJet1][0]   ->Fill(DeltaEta-1., DeltaPhi);
		if(-1.5<DeltaEta && DeltaEta<-1.)
		  fh2JetHadCorrNearSideBand[EPBinJet1][0]   ->Fill(DeltaEta+1., DeltaPhi);
		
		if(PtBin>0) {
		  fh2JetHadCorr[EPBinJet1][PtBin]            ->Fill(DeltaEta, DeltaPhi);
		  fh2JetHadCorrNear[EPBinJet1][PtBin]        ->Fill(DeltaEta, DeltaPhi);
		  if(1.<DeltaEta && DeltaEta<1.5)
		    fh2JetHadCorrNearSideBand[EPBinJet1][PtBin]->Fill(DeltaEta-1., DeltaPhi);
		  if(-1.5<DeltaEta && DeltaEta<-1.)
		    fh2JetHadCorrNearSideBand[EPBinJet1][PtBin]->Fill(DeltaEta+1., DeltaPhi);
		}
	      }   
	    }
	  }    
	}

 	if(pool->IsReady()) {
 	  // Fill mixed-event histos here  
 	  for (Int_t jMix=0; jMix < nMix; jMix++) {
 	    TObjArray* bgTracks = pool->GetEvent(jMix);
 	    for(Int_t ibg=0; ibg < bgTracks->GetEntries(); ibg++){		    
 	      AliBasicParticle *bgTrack = static_cast<AliBasicParticle*>(bgTracks->At(ibg));
 	      if(!bgTrack)
 		{
 		  AliError(Form("%s:Failed to retrieve tracks from mixed events", GetName()));
 		}

	      PtBin=-1;
	      DeltaEta=999;
	      DeltaPhi=999;	  
	      for(Int_t ipt=0; ipt<5; ipt++) {
		if(fDivPt[ipt]<=bgTrack->Pt() && bgTrack->Pt()<fDivPt[ipt+1]) {
		  PtBin = ipt+1;
		  break;
		}
	      }

	      if(jetrecalc) {
		DeltaEta = bgTrack->Eta()-jetrecalc->Phi();
		DeltaPhi = TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(bgTrack->Phi())-jetrecalc->Phi());
	      }
	      else if(!(fMinPtForJetRecalc>0.)) {
		DeltaEta = bgTrack->Eta()-jet1->Phi();
		DeltaPhi = TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(bgTrack->Phi())-jet1->Phi());
	      }

	      if(DeltaEta<999 && DeltaPhi<999) {
		fh2JetHadCorrMix[0][0]            ->Fill(DeltaEta, DeltaPhi);
		fh2JetHadCorrMixNear[0][0]        ->Fill(DeltaEta, DeltaPhi);
		if(1.<DeltaEta && DeltaEta<1.5)
		  fh2JetHadCorrMixNearSideBand[0][0]->Fill(DeltaEta-1., DeltaPhi);
		if(-1.5<DeltaEta && DeltaEta<-1.)
		  fh2JetHadCorrMixNearSideBand[0][0]->Fill(DeltaEta+1., DeltaPhi);
	      
		
		if(fEP!=999) {
		  if(jetrecalc) {
		    fh2JetHadCorrMix[EPBinRecalcJet1][0]               ->Fill(DeltaEta, DeltaPhi);
		    fh2JetHadCorrMixNear[EPBinRecalcJet1][0]           ->Fill(DeltaEta, DeltaPhi);
		    if(1.<DeltaEta && DeltaEta<1.5)
		    fh2JetHadCorrMixNearSideBand[EPBinRecalcJet1][0]   ->Fill(DeltaEta-1., DeltaPhi);
		    if(-1.5<DeltaEta && DeltaEta<-1.)
		      fh2JetHadCorrMixNearSideBand[EPBinRecalcJet1][0]   ->Fill(DeltaEta+1., DeltaPhi);
		    if(PtBin>0) {
		      fh2JetHadCorrMix[EPBinRecalcJet1][PtBin]            ->Fill(DeltaEta, DeltaPhi);
		      fh2JetHadCorrMixNear[EPBinRecalcJet1][PtBin]        ->Fill(DeltaEta, DeltaPhi);
		      if(1.<DeltaEta && DeltaEta<1.5)
			fh2JetHadCorrMixNearSideBand[EPBinRecalcJet1][PtBin]->Fill(DeltaEta-1., DeltaPhi);
		      if(-1.5<DeltaEta && DeltaEta<-1.)
			fh2JetHadCorrMixNearSideBand[EPBinRecalcJet1][PtBin]->Fill(DeltaEta+1., DeltaPhi);
		    }
		  }
		  else if(!(fMinPtForJetRecalc>0.)) {
		    fh2JetHadCorrMix[EPBinJet1][0]               ->Fill(DeltaEta, DeltaPhi);
		    fh2JetHadCorrMixNear[EPBinJet1][0]           ->Fill(DeltaEta, DeltaPhi);
		    if(1.<DeltaEta && DeltaEta<1.5)
		      fh2JetHadCorrMixNearSideBand[EPBinJet1][0]   ->Fill(DeltaEta-1., DeltaPhi);
		    if(-1.5<DeltaEta && DeltaEta<-1.)
		      fh2JetHadCorrMixNearSideBand[EPBinJet1][0]   ->Fill(DeltaEta+1., DeltaPhi);
		    
		    if(PtBin>0) {
		      fh2JetHadCorrMix[EPBinJet1][PtBin]            ->Fill(DeltaEta, DeltaPhi);
		      fh2JetHadCorrMixNear[EPBinJet1][PtBin]        ->Fill(DeltaEta, DeltaPhi);
		      if(1.<DeltaEta && DeltaEta<1.5)
			fh2JetHadCorrMixNearSideBand[EPBinJet1][PtBin]->Fill(DeltaEta-1., DeltaPhi);
		      if(-1.5<DeltaEta && DeltaEta<-1.)
			fh2JetHadCorrMixNearSideBand[EPBinJet1][PtBin]->Fill(DeltaEta+1., DeltaPhi);
		    }
		  }
		}
	      }
	    }
	  }
	}

	Double_t DeltaPhiEmbVsDet  = TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jet1->Phi())-TVector2::Phi_0_2pi(jet2->Phi()));
	Double_t DeltaEtaEmbVsDet  = jet1->Eta()-jet2->Eta();
	Double_t DeltaPhiEmbVsPart = TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jet1->Phi())-TVector2::Phi_0_2pi(jet3->Phi()));
	Double_t DeltaEtaEmbVsPart = jet1->Eta()-jet3->Eta();
	Double_t DeltaPhiDetVsPart = TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jet2->Phi())-TVector2::Phi_0_2pi(jet3->Phi()));
	Double_t DeltaEtaDetVsPart = jet2->Eta()-jet3->Eta();                                                                

	Double_t DeltaPhiRecalcEmbVsRecalcDet  = 999;
	Double_t DeltaEtaRecalcEmbVsRecalcDet  = 999;
	Double_t DeltaPhiRecalcEmbVsRecalcPart = 999;
	Double_t DeltaEtaRecalcEmbVsRecalcPart = 999;
	Double_t DeltaPhiRecalcDetVsRecalcPart = 999;
	Double_t DeltaEtaRecalcDetVsRecalcPart = 999;

	Double_t DeltaPhiRecalcEmbVsDet  = 999;
	Double_t DeltaEtaRecalcEmbVsDet  = 999;
	Double_t DeltaPhiRecalcEmbVsPart = 999;
	Double_t DeltaEtaRecalcEmbVsPart = 999;
	Double_t DeltaPhiRecalcDetVsPart = 999;
	Double_t DeltaEtaRecalcDetVsPart = 999;   

	if(jetrecalc) {
	  DeltaPhiRecalcEmbVsDet  = TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jetrecalc->Phi())-TVector2::Phi_0_2pi(jet2->Phi()));
	  DeltaEtaRecalcEmbVsDet  = jetrecalc->Eta()-jet2->Eta();                                                                
	  DeltaPhiRecalcEmbVsPart = TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jetrecalc->Phi())-TVector2::Phi_0_2pi(jet3->Phi()));
	  DeltaEtaRecalcEmbVsPart = jetrecalc->Eta()-jet3->Eta();                                                                
	  if(jetRecalcDetPythia) { 
	    DeltaPhiRecalcEmbVsRecalcDet  = TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jetrecalc->Phi())-TVector2::Phi_0_2pi(jetRecalcDetPythia->Phi()));
	    DeltaEtaRecalcEmbVsRecalcDet  = jetrecalc->Eta()-jetRecalcDetPythia->Eta();                                                                
	  }
	  if(jetRecalcPartPythia) {
	    DeltaPhiRecalcEmbVsRecalcPart = TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jetrecalc->Phi())-TVector2::Phi_0_2pi(jetRecalcPartPythia->Phi()));
	    DeltaEtaRecalcEmbVsRecalcPart = jetrecalc->Eta()-jetRecalcPartPythia->Eta();                                                                
	  }
	}

	if(jetRecalcDetPythia) {
	  DeltaPhiRecalcDetVsPart = TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jetRecalcDetPythia->Phi())-TVector2::Phi_0_2pi(jet3->Phi()));
	  DeltaEtaRecalcDetVsPart = jetRecalcDetPythia->Eta()-jet3->Eta();                                                                
	  if(jetRecalcPartPythia) {
	    DeltaPhiRecalcDetVsRecalcPart = TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jetRecalcDetPythia->Phi())-TVector2::Phi_0_2pi(jetRecalcPartPythia->Phi()));
	    DeltaEtaRecalcDetVsRecalcPart = jetRecalcDetPythia->Eta()-jetRecalcPartPythia->Eta();                                                                
	  }
	}

	fh2JetdEtadPhiEmbVsDet             [0]->Fill(DeltaEtaEmbVsDet , DeltaPhiEmbVsDet);
	fh2JetdEtadPhiEmbVsPart            [0]->Fill(DeltaEtaEmbVsPart, DeltaPhiEmbVsPart);
	fh2JetdEtadPhiDetVsPart               ->Fill(DeltaEtaDetVsPart, DeltaPhiDetVsPart);
	if(jetrecalc) {
	  fh2JetdEtadPhiRecalcEmbVsDet       [0]->Fill(DeltaEtaRecalcEmbVsDet, DeltaPhiRecalcEmbVsDet);
	  fh2JetdEtadPhiRecalcEmbVsPart      [0]->Fill(DeltaEtaRecalcEmbVsPart, DeltaPhiRecalcEmbVsPart);
	  if(jetRecalcDetPythia)
	    fh2JetdEtadPhiRecalcEmbVsRecalcDet [0]->Fill(DeltaEtaRecalcEmbVsRecalcDet, DeltaPhiRecalcEmbVsRecalcDet);
	  if(jetRecalcPartPythia)
	    fh2JetdEtadPhiRecalcEmbVsRecalcPart[0]->Fill(DeltaEtaRecalcEmbVsRecalcPart, DeltaPhiRecalcDetVsRecalcPart);
	}

	if(jetRecalcDetPythia) {
	  fh2JetdEtadPhiRecalcDetVsPart->Fill(DeltaEtaRecalcDetVsPart, DeltaPhiRecalcDetVsPart);
	  if(jetRecalcPartPythia) {
	    fh2JetdEtadPhiRecalcDetVsRecalcPart->Fill(DeltaEtaRecalcDetVsRecalcPart, DeltaPhiRecalcDetVsRecalcPart);
	  }

	  if(fEP!=999) {
	    fh2JetdEtadPhiEmbVsDet             [EPBinJet1]->Fill(DeltaEtaEmbVsDet , DeltaPhiEmbVsDet);
	    fh2JetdEtadPhiEmbVsPart            [EPBinJet1]->Fill(DeltaEtaEmbVsPart, DeltaPhiEmbVsPart);
	    if(jetrecalc) {
	      fh2JetdEtadPhiRecalcEmbVsDet       [EPBinRecalcJet1]->Fill(DeltaEtaRecalcEmbVsDet, DeltaPhiRecalcEmbVsDet);
	      fh2JetdEtadPhiRecalcEmbVsPart      [EPBinRecalcJet1]->Fill(DeltaEtaRecalcEmbVsPart, DeltaPhiRecalcEmbVsPart);
	      if(jetRecalcDetPythia)
		fh2JetdEtadPhiRecalcEmbVsRecalcDet [EPBinRecalcJet1]->Fill(DeltaEtaRecalcEmbVsRecalcDet, DeltaPhiRecalcEmbVsRecalcDet);
	      if(jetRecalcPartPythia)
		fh2JetdEtadPhiRecalcEmbVsRecalcPart[EPBinRecalcJet1]->Fill(DeltaEtaRecalcEmbVsRecalcPart, DeltaPhiRecalcDetVsRecalcPart);
	    }
	  }

	}
	if(jetrecalc)
	  jetrecalc->Delete();
	if(jetRecalcDetPythia)
	  jetRecalcDetPythia->Delete();
	if(jetRecalcPartPythia)
	  jetRecalcPartPythia->Delete();
      }
    
      
      

      if (fJetShapeType == kPythiaDef){
        
	AliJetContainer *jetContTrue = GetJetContainer(1);
	AliJetContainer *jetContUS = GetJetContainer(2);
	AliJetContainer *jetContPart = GetJetContainer(3);
        
	if(fJetShapeSub==kConstSub){
          
	  for(Int_t i = 0; i<jetContUS->GetNJets(); i++) {
	    jetUS = jetContUS->GetJet(i);
	    if(jetUS->GetLabel()==jet1->GetLabel()) {
	      ifound++;
	      if(ifound==1) ilab = i;
	    }
	  }
	  if(ilab==-1) continue;
	  jetUS=jetContUS->GetJet(ilab);
	  jet2=jetUS->ClosestJet();
        
	  if (!jet2) {
	    Printf("jet2 does not exist, returning");
	    continue;
	  }
          
	  for(Int_t j=0; j<jetContPart->GetNJets(); j++) {
            
	    jet3 = jetContPart->GetJet(j);
	    if(!jet3) continue;
	    if(jet3->GetLabel()==jet2->GetLabel()) {
	      jfound++;
	      if(jfound==1) jlab = j;
	    }
	  }
	  if(jlab==-1) continue;
	  jet3=jetContPart->GetJet(jlab);
	  if(!jet3){
	    Printf("jet3 does not exist, returning");
	    continue;
	  }
	}
	if(!(fJetShapeSub==kConstSub)) jet3 = jet1->ClosestJet();
	if (!jet3) {
	  Printf("jet3 does not exist, returning");
	  continue;
	}
        
      
	fh2ResponseUW->Fill(jet1->Pt(),jet3->Pt());
        
        
      }
      
      
      if (fJetShapeType == kGenOnTheFly){
	const AliEmcalPythiaInfo *partonsInfo = 0x0;
	partonsInfo = GetPythiaInfo();
	Double_t jp1=RelativePhi(jet1->Phi(),partonsInfo->GetPartonPhi6());
	Double_t detap1=(jet1->Eta())-(partonsInfo->GetPartonEta6());
	kWeight=partonsInfo->GetPythiaEventWeight();
	fh2ResponseW->Fill(jet1->Pt(),jet1->Pt(),kWeight);
        
	Float_t dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
	fEtaJetCorr6->Fill(jet1->Eta(), partonsInfo->GetPartonEta6());
	fPhiJetCorr6->Fill(jet1->Phi(), partonsInfo->GetPartonPhi6());
	if(dRp1 < fRMatching) {
	  fShapesVar[0] = partonsInfo->GetPartonFlag6();
	  fPtJetCorr ->Fill(partonsInfo->GetPartonPt6(), jet1->Pt());
	}
	else {
	  jp1=RelativePhi(jet1->Phi(),partonsInfo->GetPartonPhi7());
	  detap1=(jet1->Eta())-(partonsInfo->GetPartonEta7());
	  dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
	  fEtaJetCorr7->Fill(jet1->Eta(), partonsInfo->GetPartonEta7());
	  fPhiJetCorr7->Fill(jet1->Phi(), partonsInfo->GetPartonPhi7());
	  if(dRp1 < fRMatching) {
	    fShapesVar[0] = partonsInfo->GetPartonFlag7();
	    fPtJetCorr->Fill(partonsInfo->GetPartonPt7(), jet1->Pt());
	  }
	  else fShapesVar[0]=0;
	}
      }
    
    
    
      
      Double_t ptSubtracted = 0;
      if (fJetShapeSub==kConstSub) ptSubtracted= jet1->Pt();
      
      else if (fJetShapeSub==kDerivSub)  {
	ptSubtracted=jet1->Pt()-GetRhoVal(0)*jet1->Area();
      }
      
      else if (fJetShapeSub==kNoSub) {
	if ((fJetShapeType==kData) || (fJetShapeType==kDetEmbPartPythia)) ptSubtracted=jet1->Pt()-GetRhoVal(0)*jet1->Area();
	else if ((fJetShapeType==kPythiaDef) || (fJetShapeType==kMCTrue) || (fJetShapeType==kGenOnTheFly)) ptSubtracted= jet1->Pt();
      }

      //Printf("ptSubtracted=%f,fPtThreshold =%f ", ptSubtracted, fPtThreshold);
      if (ptSubtracted < fPtThreshold) continue;
      
      if (fOneConstSelectOn == kTRUE) fNbOfConstvspT->Fill(GetJetNumberOfConstituents(jet1,0), ptSubtracted);
      if ((fCentSelectOn == kFALSE) && (jet1->GetNumberOfTracks() <= 1)) continue;

  
      fShapesVar[1] = ptSubtracted;
      fShapesVar[2] = GetJetpTD(jet1,0);
      fShapesVar[3] =jet1->Phi();
      if(fJetShapeType==kData) fShapesVar[3]=RelativePhi(triggerHadron->Phi(), jet1->Phi());
      //GetJetMass(jet1,0);
      fShapesVar[4] = GetJetAngularity(jet1,0);
      //fShapesVar[5] = GetJetCircularity(jet1,0);
      fShapesVar[5] = GetJetLeSub(jet1,0);
      //fShapesVar[6] = GetJetCoronna(jet1,0);
      RecursiveParents(jet1,jetCont,0);
      RecursiveParents(jet1,jetCont,1);
      RecursiveParents(jet1,jetCont,2);
      
      Float_t ptMatch=0., ptDMatch=0., massMatch=0., constMatch=0.,angulMatch=0.,circMatch=0., lesubMatch=0., sigma2Match=0., coronnaMatch=0;
      Int_t kMatched = 0;

      if (fJetShapeType==kPythiaDef) {
	kMatched =1;
	if(fJetShapeSub==kConstSub) kMatched = 3;
        
	ptMatch=jet3->Pt();
	ptDMatch=GetJetpTD(jet3, kMatched);
	massMatch=jet3->Phi();
	// GetJetMass(jet3,kMatched);
	//constMatch=1.*GetJetNumberOfConstituents(jet2,kMatched);
	angulMatch=GetJetAngularity(jet3, kMatched);
	//circMatch=GetJetCircularity(jet3, kMatched);
	lesubMatch=GetJetLeSub(jet3, kMatched);
	//coronnaMatch=GetJetCoronna(jet3,kMatched); 
	//sigma2Match = GetSigma2(jet2, kMatched);
      }
      
      if (fJetShapeType==kDetEmbPartPythia) {
	if(fJetShapeSub==kConstSub) kMatched = 3;
	if(fJetShapeSub==kDerivSub) kMatched = 2;
	ptMatch=jet3->Pt();
	ptDMatch=GetJetpTD(jet3, kMatched);
	massMatch=jet3->Phi();
	//GetJetMass(jet3,kMatched);
	// constMatch=1.*GetJetNumberOfConstituents(jet3,kMatched);
	angulMatch=GetJetAngularity(jet3, kMatched);
	// circMatch=GetJetCircularity(jet3, kMatched);
	lesubMatch=GetJetLeSub(jet3, kMatched);
	//coronnaMatch = GetJetCoronna(jet3, kMatched);
        
      }


       
      if (fJetShapeType == kMCTrue || fJetShapeType == kData || fJetShapeType == kGenOnTheFly) {
	kMatched = 0;
	ptMatch=0.;
	ptDMatch=0.;
	massMatch=0.;
	//constMatch=0.;
	angulMatch=0.;
	// circMatch=0.;
	lesubMatch=0.;
	//coronnaMatch =0.;
        
      }
      
    

      fShapesVar[6] = ptMatch;
      fShapesVar[7] = ptDMatch;
      fShapesVar[8] = massMatch;
      fShapesVar[9] = angulMatch;
      //fShapesVar[12] = circMatch;
      fShapesVar[10] = lesubMatch;
      //  fShapesVar[14] = coronnaMatch;
      fShapesVar[11] = kWeight;
      //fShapesVar[16] = ntracksEvt;
      // fShapesVar[16] = rhoVal;
      //fShapesVar[17] = rhoMassVal;
      //fShapesVar[16] = jet1->Pt();


      fTreeObservableTagging->Fill();
      





    }
 
  }
  tracksClone = CloneAndReduceTrackList();	
  //update pool if jet in event or not
  pool->UpdatePool(tracksClone);

  if(jetrecalc)
    jetrecalc->Delete();
  return kTRUE;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetAxisRes::GetJetMass(AliEmcalJet *jet,Int_t jetContNb) {
  //calc subtracted jet mass
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtracted();
    else return jet->GetShapeProperties()->GetSecondOrderSubtracted();
  else 
    return jet->M();
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetAxisRes::Angularity(AliEmcalJet *jet, Int_t jetContNb){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0; 
  Double_t den=0.;
  Double_t num = 0.;
  AliVParticle *vp1 = 0x0;
  for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
      
    if (!vp1){
      Printf("AliVParticle associated to constituent not found");
      continue;
    }
      
    Double_t dphi = RelativePhi(vp1->Phi(),jet->Phi());
    Double_t dr2 = (vp1->Eta()-jet->Eta())*(vp1->Eta()-jet->Eta()) + dphi*dphi;
    Double_t dr = TMath::Sqrt(dr2);
    num=num+vp1->Pt()*dr;
    den=den+vp1->Pt();
  }
  return num/den;
} 

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetAxisRes::GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb){

  if((fJetShapeSub==kDerivSub) && (jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedAngularity();
    else return jet->GetShapeProperties()->GetSecondOrderSubtractedAngularity();
  else
    return Angularity(jet, jetContNb);
 
}

//____________________________________________________________________________

Float_t AliAnalysisTaskEmcalJetAxisRes::Coronna(AliEmcalJet *jet, Int_t jetContNb){

  AliTrackContainer *PartCont = NULL;
  AliParticleContainer *PartContMC = NULL;


  if (fJetShapeSub==kConstSub){
    if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly) PartContMC = GetParticleContainer(1);
    else PartCont = GetTrackContainer(1);
  }
  else{
    if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly) PartContMC = GetParticleContainer(0);
    else PartCont = GetTrackContainer(0);
  }

  TClonesArray *TracksArray = NULL;
  TClonesArray *TracksArrayMC = NULL;

  if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly) TracksArrayMC = PartContMC->GetArray();
  else TracksArray = PartCont->GetArray();
 
  if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly){
    if(!PartContMC || !TracksArrayMC) return -2;
  }
  else {
    if(!PartCont || !TracksArray) return -2;
  }


  AliAODTrack *Track = 0x0;
  Float_t sumpt=0;
  Int_t NTracks=0;
  if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly) NTracks = TracksArrayMC->GetEntriesFast();
  else NTracks = TracksArray->GetEntriesFast();

  for(Int_t i=0; i < NTracks; i++){
    if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly){
      if((Track = static_cast<AliAODTrack*>(PartContMC->GetAcceptParticle(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	Double_t dphi = RelativePhi(Track->Phi(),jet->Phi());
	Double_t dr2 = (Track->Eta()-jet->Eta())*(Track->Eta()-jet->Eta()) + dphi*dphi;
	Double_t dr = TMath::Sqrt(dr2);
	if((dr>=0.8) && (dr<1)) sumpt=sumpt+Track->Pt();
      }
    }
    else{ 
      if((Track = static_cast<AliAODTrack*>(PartCont->GetAcceptTrack(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	if (Track->Pt()<0.15) continue;
	Double_t dphi = RelativePhi(Track->Phi(),jet->Phi());
	Double_t dr2 = (Track->Eta()-jet->Eta())*(Track->Eta()-jet->Eta()) + dphi*dphi;
	Double_t dr = TMath::Sqrt(dr2);
	if((dr>=0.8) && (dr<1)) sumpt=sumpt+Track->Pt();

      }
    } 
  }
 

  
  return sumpt; 
 


  
} 

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetAxisRes::GetJetCoronna(AliEmcalJet *jet, Int_t jetContNb){

  if((fJetShapeSub==kDerivSub) && (jetContNb==0)) return -2;
  else
    return Coronna(jet, jetContNb);
 
}





//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetAxisRes::PTD(AliEmcalJet *jet, Int_t jetContNb){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0; 
  Double_t den=0.;
  Double_t num = 0.;
  AliVParticle *vp1 = 0x0;
  for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
      
    if (!vp1){
      Printf("AliVParticle associated to constituent not found");
      continue;
    }
      
    num=num+vp1->Pt()*vp1->Pt();
    den=den+vp1->Pt();
  }
  return TMath::Sqrt(num)/den;
} 

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetAxisRes::GetJetpTD(AliEmcalJet *jet, Int_t jetContNb){
  //calc subtracted jet mass
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedpTD();
    else return jet->GetShapeProperties()->GetSecondOrderSubtractedpTD();
  else
    return PTD(jet, jetContNb);
 
}

//_____________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetAxisRes::Circularity(AliEmcalJet *jet, Int_t jetContNb){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0; 
  Double_t mxx    = 0.;
  Double_t myy    = 0.;
  Double_t mxy    = 0.;
  int  nc     = 0;
  Double_t sump2  = 0.;
  Double_t pxjet=jet->Px();
  Double_t pyjet=jet->Py();
  Double_t pzjet=jet->Pz();
  
  
  //2 general normalized vectors perpendicular to the jet
  TVector3  ppJ1(pxjet, pyjet, pzjet);
  TVector3  ppJ3(- pxjet* pzjet, - pyjet * pzjet, pxjet * pxjet + pyjet * pyjet);
  ppJ3.SetMag(1.);
  TVector3  ppJ2(-pyjet, pxjet, 0);
  ppJ2.SetMag(1.);
  AliVParticle *vp1 = 0x0;
  for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));  
    
    if (!vp1){
      Printf("AliVParticle associated to constituent not found");
      continue;
    }
    
    TVector3 pp(vp1->Px(), vp1->Py(), vp1->Pz());
   
    //local frame
    TVector3 pLong = pp.Dot(ppJ1) / ppJ1.Mag2() * ppJ1;
    TVector3 pPerp = pp - pLong;
    //projection onto the two perpendicular vectors defined above
    
    Float_t ppjX = pPerp.Dot(ppJ2);
    Float_t ppjY = pPerp.Dot(ppJ3);
    Float_t ppjT = TMath::Sqrt(ppjX * ppjX + ppjY * ppjY);
    if(ppjT<=0) return 0;
    
    mxx += (ppjX * ppjX / ppjT);
    myy += (ppjY * ppjY / ppjT);
    mxy += (ppjX * ppjY / ppjT);
    nc++;
    sump2 += ppjT;}
  
  if(nc<2) return 0;
  if(sump2==0) return 0;
  // Sphericity Matrix
  Double_t ele[4] = {mxx / sump2, mxy / sump2, mxy / sump2, myy / sump2};
  TMatrixDSym m0(2,ele);
  
  // Find eigenvectors
  TMatrixDSymEigen m(m0);
  TVectorD eval(2);
  TMatrixD evecm = m.GetEigenVectors();
  eval  = m.GetEigenValues();
  // Largest eigenvector
  int jev = 0;
  //  cout<<eval[0]<<" "<<eval[1]<<endl;
  if (eval[0] < eval[1]) jev = 1;
  TVectorD evec0(2);
  // Principle axis
  evec0 = TMatrixDColumn(evecm, jev);
  Double_t compx=evec0[0];
  Double_t compy=evec0[1];
  TVector2 evec(compx, compy);
  Double_t circ=0;
  if(jev==1) circ=2*eval[0];
  if(jev==0) circ=2*eval[1];
  
  return circ;
  
  
  
}




//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetAxisRes::GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb){
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedCircularity();
    else return jet->GetShapeProperties()->GetSecondOrderSubtractedCircularity();
  else
    return Circularity(jet, jetContNb);
 
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetAxisRes::LeSub(AliEmcalJet *jet, Int_t jetContNb){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0;
  Double_t den=0.;
  Double_t num = 0.;
  AliVParticle *vp1 = 0x0;
  AliVParticle *vp2 = 0x0;
  std::vector<int> ordindex;
  ordindex=jet->GetPtSortedTrackConstituentIndexes(jetCont->GetParticleContainer()->GetArray());
  //Printf("Nbof const = %d", jet->GetNumberOfTracks());
  //Printf("ordindex[0] = %d, ordindex[1] = %d", ordindex[0], ordindex[1]);
  
  if(ordindex.size()<2) return -1;
  
  vp1 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[0], jetCont->GetParticleContainer()->GetArray()));
  if (!vp1){
    Printf("AliVParticle associated to Leading constituent not found");
    return -1;
  }
  
  vp2 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[1], jetCont->GetParticleContainer()->GetArray()));
  if (!vp2){
    Printf("AliVParticle associated to Subleading constituent not found");
    return -1;
  }
  
  
  num=vp1->Pt();
  den=vp2->Pt();
  //Printf("vp1->Pt() =%f, vp2->Pt() =%f", vp1->Pt(), vp2->Pt());
  
  return num-den;
} 

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetAxisRes::GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb) {
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedLeSub();
    else return jet->GetShapeProperties()->GetSecondOrderSubtractedLeSub();
  else
    return LeSub(jet, jetContNb);
 
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetAxisRes::GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb){
  //calc subtracted jet mass
  
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedConstituent();
    else return jet->GetShapeProperties()->GetSecondOrderSubtractedConstituent();
  else
    return jet->GetNumberOfTracks();
 
}
   
 
//______________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetAxisRes::Sigma2(AliEmcalJet *jet, Int_t jetContNb){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0; 
  Double_t mxx    = 0.;
  Double_t myy    = 0.;
  Double_t mxy    = 0.;
  int  nc     = 0;
  Double_t sump2  = 0.;
       
  AliVParticle *vp1 = 0x0;
  for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
       
    if (!vp1){
      Printf("AliVParticle associated to constituent not found");
      continue;
    }
       
    Double_t ppt=vp1->Pt();
    Double_t dphi = RelativePhi(vp1->Phi(),jet->Phi());
     
    Double_t deta = vp1->Eta()-jet->Eta();
    mxx += ppt*ppt*deta*deta;
    myy += ppt*ppt*dphi*dphi;
    mxy -= ppt*ppt*deta*TMath::Abs(dphi);
    nc++;
    sump2 += ppt*ppt;
       
  }  
  if(nc<2) return 0;
  if(sump2==0) return 0;
  // Sphericity Matrix
  Double_t ele[4] = {mxx , mxy , mxy , myy };
  TMatrixDSym m0(2,ele);
     
  // Find eigenvectors
  TMatrixDSymEigen m(m0);
  TVectorD eval(2);
  TMatrixD evecm = m.GetEigenVectors();
  eval  = m.GetEigenValues();
  // Largest eigenvector
  int jev = 0;
  //  cout<<eval[0]<<" "<<eval[1]<<endl;
  if (eval[0] < eval[1]) jev = 1;
  TVectorD evec0(2);
  // Principle axis
  evec0 = TMatrixDColumn(evecm, jev);
  Double_t compx=evec0[0];
  Double_t compy=evec0[1];
  TVector2 evec(compx, compy);
  Double_t sig=0;
  if(jev==1) sig=TMath::Sqrt(TMath::Abs(eval[0])/sump2);
  if(jev==0) sig=TMath::Sqrt(TMath::Abs(eval[1])/sump2);
     
  return sig;
     
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetAxisRes::GetSigma2(AliEmcalJet *jet, Int_t jetContNb){
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedSigma2();
    else return jet->GetShapeProperties()->GetSecondOrderSubtractedSigma2();
  else
    return Sigma2(jet, jetContNb);
 
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetAxisRes::SelectTrigger(Float_t minpT, Float_t maxpT){

  AliTrackContainer *PartCont = NULL;
  AliParticleContainer *PartContMC = NULL;


  if (fJetShapeSub==kConstSub){
    if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly) PartContMC = GetParticleContainer(1);
    else PartCont = GetTrackContainer(1);
  }
  else{
    if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly) PartContMC = GetParticleContainer(0);
    else PartCont = GetTrackContainer(0);
  }

  TClonesArray *TracksArray = NULL;
  TClonesArray *TracksArrayMC = NULL;

  if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly) TracksArrayMC = PartContMC->GetArray();
  else TracksArray = PartCont->GetArray();
 
  if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly){
    if(!PartContMC || !TracksArrayMC) return -99999;
  }
  else {
    if(!PartCont || !TracksArray) return -99999;
  }


  AliAODTrack *Track = 0x0;

 
  
  TList *trackList = new TList();
  Int_t triggers[100];
  for (Int_t iTrigger=0; iTrigger<100; iTrigger++) triggers[iTrigger] = 0;
  Int_t iTT = 0;
  Int_t NTracks=0;
  if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly) NTracks = TracksArrayMC->GetEntriesFast();
  else NTracks = TracksArray->GetEntriesFast();

  for(Int_t i=0; i < NTracks; i++){
    if (fJetShapeType == AliAnalysisTaskEmcalJetAxisRes::kGenOnTheFly){
      if((Track = static_cast<AliAODTrack*>(PartContMC->GetAcceptParticle(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	if (Track->Pt()<0.15) continue;
	if ((Track->Pt() >= minpT) && (Track->Pt()< maxpT)) {
	  triggers[iTT] = i;
	  iTT++;
	}
      }
    }
    else{ 
      if((Track = static_cast<AliAODTrack*>(PartCont->GetAcceptTrack(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	if (Track->Pt()<0.15) continue;
	if ((Track->Pt() >= minpT) && (Track->Pt()< maxpT)) {
	  triggers[iTT] = i;
	  iTT++;
	}
      }
    } 
  }
 

  if (iTT == 0) return -99999;
  Int_t nbRn = 0, index = 0 ; 
  TRandom3* random = new TRandom3(0); 
  nbRn = random->Integer(iTT);
  index = triggers[nbRn];
  //Printf("iTT Total= %d, nbRn = %d, Index = %d",iTT, nbRn, index );
  return index; 
  
}

//__________________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetAxisRes::RelativePhi(Double_t mphi,Double_t vphi){

  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
  return dphi;//dphi in [-Pi, Pi]
}


//_________________________________________________________________________
void AliAnalysisTaskEmcalJetAxisRes::RecursiveParents(AliEmcalJet *fJet,AliJetContainer *fJetCont, Int_t ReclusterAlgo){
 
  std::vector<fastjet::PseudoJet>  fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet  PseudoTracks;
  double xflagalgo=0; 
  AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();
  
  if (fTrackCont) for (Int_t i=0; i<fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk) continue; 
      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
      PseudoTracks.set_user_index(fJet->TrackAt(i)+100);
      fInputVectors.push_back(PseudoTracks);
     
    }
  fastjet::JetAlgorithm jetalgo(fastjet::antikt_algorithm);

  if(ReclusterAlgo==0){ xflagalgo=0.5;
    jetalgo=fastjet::kt_algorithm ;}
      
  if(ReclusterAlgo==1){ xflagalgo=1.5;
    jetalgo=fastjet::cambridge_algorithm;}
  if(ReclusterAlgo==2){ xflagalgo=2.5;
    jetalgo=fastjet::antikt_algorithm;} 
  
  fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 

  try {
    fastjet::ClusterSequence fClustSeqSA(fInputVectors, fJetDef);
    std::vector<fastjet::PseudoJet>   fOutputJets;
    fOutputJets.clear();
    fOutputJets=fClustSeqSA.inclusive_jets(0);
  
    fastjet::PseudoJet jj;
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;
    jj=fOutputJets[0];
    double ndepth=0;
    while(jj.has_parents(j1,j2)){
      ndepth=ndepth+1;
      if(j1.perp() < j2.perp()) swap(j1,j2);
      double delta_R=j1.delta_R(j2);
      double z=j2.perp()/(j1.perp()+j2.perp());
      double y =log(1.0/delta_R);
      double lnpt_rel=log(z*delta_R);
      Double_t LundEntries[5] = {y,lnpt_rel,fOutputJets[0].perp(),xflagalgo,ndepth};  
      fHLundIterative->Fill(LundEntries);
      jj=j1;} 




  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    //return -1;
  }




  return;

}







//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetAxisRes::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskEmcalJetAxisRes::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.

  // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
  // if (!fTreeObservableTagging){
  //   Printf("ERROR: fTreeObservableTagging not available"); 
  //   return;
  // }

}

/**
 * Clone tracks into a lighter object (AliBasicParticle) to store in the event pool. By using
 * lighter objects, it reduces the event pool size in memory. Adapted from CF event mixing
 * code PhiCorrelations.
 *
 * @return Array containing the lighter track objects.
 * Quoted from AliAnalysisTaskEmcalJetHCorrelations
 */
TObjArray* AliAnalysisTaskEmcalJetAxisRes::CloneAndReduceTrackList()
{
  // clones a track list by using AliBasicTrack which uses much less memory (used for event mixing)
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  // Loop over all tracks
  AliBasicParticle * clone = 0;
  //AliTrackContainer * tracks = GetTrackContainer("tracksForCorrelations");
  AliTrackContainer * tracks = GetTrackContainer(0);

  for (auto particle : tracks->accepted())
  {
    // Fill some QA information about the tracks
    //Int_t trackPtBin = GetTrackPtBin(particle->Pt());
    //if(trackPtBin > -1) fHistTrackEtaPhi[trackPtBin]->Fill(particle->Eta(),particle->Phi());

    // Create new particle
    clone = new AliBasicParticle(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge());
    // Set so that we can do comparisons using the IsEqual() function.
    clone ->SetUniqueID(particle->GetUniqueID());

    tracksClone->Add(clone);
  }

  return tracksClone;
}

