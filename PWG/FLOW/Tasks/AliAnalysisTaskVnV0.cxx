#include "AliAnalysisTaskVnV0.h"

// ROOT includes
#include <TMath.h>

// AliRoot includes
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliCentrality.h"
#include "AliVHeader.h"
#include "AliAODVZERO.h"
#include "TFile.h"
#include "AliOADBContainer.h"
#include "TH2F.h"
#include "TF1.h"
#include "AliGenHijingEventHeader.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "TChain.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"

// STL includes
//#include <iostream>
//using namespace std;

ClassImp(AliAnalysisTaskVnV0)
Bool_t AliAnalysisTaskVnV0::fgIsPsiComputed = kFALSE;
Float_t AliAnalysisTaskVnV0::fgPsi2v0a=999.;
Float_t AliAnalysisTaskVnV0::fgPsi2v0c=999.;
Float_t AliAnalysisTaskVnV0::fgPsi2tpc=999.;
Float_t AliAnalysisTaskVnV0::fgPsi3v0a=999.;
Float_t AliAnalysisTaskVnV0::fgPsi3v0c=999.;
Float_t AliAnalysisTaskVnV0::fgPsi3tpc=999.;
Float_t AliAnalysisTaskVnV0::fgPsi2v0aMC=999.;
Float_t AliAnalysisTaskVnV0::fgPsi2v0cMC=999.;
Float_t AliAnalysisTaskVnV0::fgPsi2tpcMC=999.;
Float_t AliAnalysisTaskVnV0::fgPsi3v0aMC=999.;
Float_t AliAnalysisTaskVnV0::fgPsi3v0cMC=999.;
Float_t AliAnalysisTaskVnV0::fgPsi3tpcMC=999.;

//_____________________________________________________________________________
AliAnalysisTaskVnV0::AliAnalysisTaskVnV0():
  AliAnalysisTaskSE(),
  fVtxCut(10.0),  // cut on |vertex| < fVtxCut
  fEtaCut(0.8),   // cut on |eta| < fEtaCut
  fMinPt(0.15),   // cut on pt > fMinPt
  fMinDistV0(0),
  fMaxDistV0(100),
  fV2(kTRUE),
  fV3(kTRUE),
  fIsMC(kFALSE),
  fQAsw(kFALSE),
  fRun(-1),
  fNcluster(70),
  fList(new TList()),
  fList2(new TList()),
  fList3(new TList()),
  fList4(new TList()),
  fMultV0(NULL),
  fV0Cpol(100),
  fV0Apol(100),
  fHResTPCv0A2(NULL),
  fHResTPCv0C2(NULL),
  fHResv0Cv0A2(NULL),
  fHResTPCv0A3(NULL),
  fHResTPCv0C3(NULL),
  fHResv0Cv0A3(NULL),
  fPhiRPv0A(NULL),
  fPhiRPv0C(NULL),
  fPhiRPv0Av3(NULL),
  fPhiRPv0Cv3(NULL),
  fQA(NULL),
  fQA2(NULL),
  fQAv3(NULL),
  fQA2v3(NULL),
  fPID(new AliFlowBayesianPID()),
  fTree(NULL),
  fCentrality(-1),
  evPlAngV0ACor2(0),
  evPlAngV0CCor2(0),
  evPlAng2(0),
  evPlAngV0ACor3(0),
  evPlAngV0CCor3(0),
  evPlAng3(0),
  fContAllChargesV0A(NULL),
  fContAllChargesV0C(NULL),
  fContAllChargesV0Av3(NULL),
  fContAllChargesV0Cv3(NULL),
  fContAllChargesMC(NULL),
  fHResMA2(NULL),
  fHResMC2(NULL),
  fHResAC2(NULL),
  fHResMA3(NULL),
  fHResMC3(NULL),
  fHResAC3(NULL),  
  fContAllChargesMCA(NULL),
  fContAllChargesMCC(NULL),
  fContAllChargesMCAv3(NULL),
  fContAllChargesMCCv3(NULL),
  fFillDCA(kFALSE),
  fContQApid(NULL),
  fModulationDEDx(kFALSE),
  fCutsDaughter(NULL)
{
  // Default constructor (should not be used)
  fList->SetName("resultsV2");
  fList2->SetName("resultsV3");
  fList3->SetName("resultsMC");
  fList4->SetName("QA");

  fList->SetOwner(kTRUE); 
  fList2->SetOwner(kTRUE); 
  fList3->SetOwner(kTRUE); 
  fList4->SetOwner(kTRUE); 

  fPID->SetNewTrackParam(); // Better tuning for TOF PID tracking effect in LHC10h
}

//______________________________________________________________________________
AliAnalysisTaskVnV0::AliAnalysisTaskVnV0(const char *name):
  AliAnalysisTaskSE(name),
  fVtxCut(10.0),  // cut on |vertex| < fVtxCut
  fEtaCut(0.8),   // cut on |eta| < fEtaCut
  fMinPt(0.15),   // cut on pt > fMinPt
  fMinDistV0(0),
  fMaxDistV0(100),
  fV2(kTRUE),
  fV3(kTRUE),
  fIsMC(kFALSE),
  fQAsw(kFALSE),
  fRun(-1),
  fNcluster(70),
  fList(new TList()),
  fList2(new TList()),
  fList3(new TList()),
  fList4(new TList()),
  fMultV0(NULL),
  fV0Cpol(100),
  fV0Apol(100),
  fHResTPCv0A2(NULL),
  fHResTPCv0C2(NULL),
  fHResv0Cv0A2(NULL),
  fHResTPCv0A3(NULL),
  fHResTPCv0C3(NULL),
  fHResv0Cv0A3(NULL),
  fPhiRPv0A(NULL),
  fPhiRPv0C(NULL),
  fPhiRPv0Av3(NULL),
  fPhiRPv0Cv3(NULL),
  fQA(NULL),
  fQA2(NULL),
  fQAv3(NULL),
  fQA2v3(NULL),
  fPID(new AliFlowBayesianPID()),
  fTree(NULL),
  fCentrality(-1),
  evPlAngV0ACor2(0),
  evPlAngV0CCor2(0),
  evPlAng2(0),
  evPlAngV0ACor3(0),
  evPlAngV0CCor3(0),
  evPlAng3(0),
  fContAllChargesV0A(NULL),
  fContAllChargesV0C(NULL),
  fContAllChargesV0Av3(NULL),
  fContAllChargesV0Cv3(NULL),
  fContAllChargesMC(NULL),
  fHResMA2(NULL),
  fHResMC2(NULL),
  fHResAC2(NULL),
  fHResMA3(NULL),
  fHResMC3(NULL),
  fHResAC3(NULL),  
  fContAllChargesMCA(NULL),
  fContAllChargesMCC(NULL),
  fContAllChargesMCAv3(NULL),
  fContAllChargesMCCv3(NULL),
  fFillDCA(kFALSE),
  fContQApid(NULL),
  fModulationDEDx(kFALSE),
  fCutsDaughter(NULL)
{

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());

  // Output slot #1 writes into a TTree
  fList->SetName("resultsV2");
  fList2->SetName("resultsV3");
  fList3->SetName("resultsMC");
  fList4->SetName("QA");

  fList->SetOwner(kTRUE); 
  fList2->SetOwner(kTRUE); 
  fList3->SetOwner(kTRUE); 
  fList4->SetOwner(kTRUE); 

  fPID->SetNewTrackParam(); // Better tuning for TOF PID tracking effect in LHC10h
}
//_____________________________________________________________________________
AliAnalysisTaskVnV0::~AliAnalysisTaskVnV0()
{

}

//______________________________________________________________________________
void AliAnalysisTaskVnV0::UserCreateOutputObjects()
{

  if(fIsMC) fPID->SetMC(kTRUE);


  // Tree for EP debug (comment the adding to v2 list id not needed)
  fTree = new TTree("tree","tree");
  fTree->Branch("evPlAngV0ACor2",&evPlAngV0ACor2,"evPlAngV0ACor2/F");
  fTree->Branch("evPlAngV0CCor2",&evPlAngV0CCor2,"evPlAngV0CCor2/F");
  fTree->Branch("evPlAng2",&evPlAng2,"evPlAng2/F"); 
  fTree->Branch("fCentrality",&fCentrality,"fCentrality/F"); 
  fTree->Branch("evPlAngV0ACor3",&evPlAngV0ACor3,"evPlAngV0ACor3/F");
  fTree->Branch("evPlAngV0CCor3",&evPlAngV0CCor3,"evPlAngV0CCor3/F");
  fTree->Branch("evPlAng3",&evPlAng3,"evPlAng3/F"); 
  

  // Container analyses (different steps mean different species)
  const Int_t nPtBinsTOF = 45;
  Double_t binsPtTOF[nPtBinsTOF+1] = {0., 0.05,  0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.25, 2.5, 2.75,3.0,3.25,3.5,3.75,4.0,4.5,5,5.5,6,6.5,7,8,9,10,12,15,20};
  const Int_t nCentrTOF = nCentrBin;
  const Int_t nPsiTOF = 10;  
  const Int_t nChargeBinsTOFres = 2; 
  const Int_t nCentrTOFres = nCentrBin;
  const Int_t nProbTOFres = 4;
  const Int_t nPsiTOFres = 10;
  const Int_t nMaskPID = 3;

  Int_t nDCABin = 1; // put to 1 not to store this info
  if(fFillDCA) nDCABin = 3;
  if(fIsMC && nDCABin>1)  nDCABin = 6;
  /*
    0 = DCAxy < 2.4 && all (or Physical primary if MC)
    1 = DCAxy > 2.4 && all (or Physical primary if MC)
    2 = DCAxy < 2.4 && not Physical Primary for MC
    3 = DCAxy > 2.4 && not Physical Primary for MC
  */
  
  Int_t binsTOF[6] = {nCentrTOFres,nChargeBinsTOFres,nProbTOFres,nPsiTOFres,nMaskPID,nDCABin};
  Int_t binsTOFmc[5] = {nCentrTOFres,nChargeBinsTOFres,1,nPsiTOFres,2};
  Int_t binsTOFmcPureMC[5] = {nCentrTOFres,nChargeBinsTOFres,1,nPsiTOFres,1};

  // v2 container
  fContAllChargesV0A = new AliFlowVZEROResults("v2A",6,binsTOF);
  fContAllChargesV0A->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
  fContAllChargesV0A->SetVarRange(1,-1.5,1.5);  // charge
  fContAllChargesV0A->SetVarRange(2,0.6,1.0001);// prob
  fContAllChargesV0A->SetVarRange(3,-TMath::Pi()/2,TMath::Pi()/2); // Psi
  fContAllChargesV0A->SetVarRange(4,-0.5,2.5); // pid mask
  fContAllChargesV0A->SetVarRange(5,-0.5,nDCABin-0.5); // DCA mask
  fContAllChargesV0A->SetVarName(0,"centrality");
  fContAllChargesV0A->SetVarName(1,"charge");
  fContAllChargesV0A->SetVarName(2,"prob");
  fContAllChargesV0A->SetVarName(3,"#Psi");
  fContAllChargesV0A->SetVarName(4,"PIDmask");
  fContAllChargesV0A->SetVarName(5,"DCAbin");
  if(fV2) fContAllChargesV0A->AddSpecies("all",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0A->AddSpecies("pi",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0A->AddSpecies("k",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0A->AddSpecies("pr",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0A->AddSpecies("e",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0A->AddSpecies("d",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0A->AddSpecies("t",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0A->AddSpecies("he3",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0A->AddSpecies("mu",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0A->AddSpecies("Ks",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0A->AddSpecies("Lambda",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0A->AddSpecies("pFromLambda",nPtBinsTOF,binsPtTOF);

  fContAllChargesV0C = new AliFlowVZEROResults("v2C",6,binsTOF);
  fContAllChargesV0C->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
  fContAllChargesV0C->SetVarRange(1,-1.5,1.5);  // charge
  fContAllChargesV0C->SetVarRange(2,0.6,1.0001);// prob
  fContAllChargesV0C->SetVarRange(3,-TMath::Pi()/2,TMath::Pi()/2); // Psi
  fContAllChargesV0C->SetVarRange(4,-0.5,2.5); // pid mask
  fContAllChargesV0C->SetVarRange(5,-0.5,nDCABin-0.5); // DCA mask
  fContAllChargesV0C->SetVarName(0,"centrality");
  fContAllChargesV0C->SetVarName(1,"charge");
  fContAllChargesV0C->SetVarName(2,"prob");
  fContAllChargesV0C->SetVarName(3,"#Psi");
  fContAllChargesV0C->SetVarName(4,"PIDmask");
  fContAllChargesV0C->SetVarName(5,"DCAbin");
  if(fV2) fContAllChargesV0C->AddSpecies("all",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0C->AddSpecies("pi",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0C->AddSpecies("k",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0C->AddSpecies("pr",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0C->AddSpecies("e",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0C->AddSpecies("d",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0C->AddSpecies("t",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0C->AddSpecies("he3",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0C->AddSpecies("mu",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0C->AddSpecies("Ks",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0C->AddSpecies("Lambda",nPtBinsTOF,binsPtTOF);
  if(fV2) fContAllChargesV0C->AddSpecies("pFromLambda",nPtBinsTOF,binsPtTOF);

  fList->Add(fContAllChargesV0A);
  fList->Add(fContAllChargesV0C);

  if(fIsMC && fV2){
    fContAllChargesMC = new AliFlowVZEROResults("v2mc",5,binsTOFmc);
    fContAllChargesMC->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
    fContAllChargesMC->SetVarRange(1,-1.5,1.5);  // charge
    fContAllChargesMC->SetVarRange(2,0.6,1.0001);// prob
    fContAllChargesMC->SetVarRange(3,-TMath::Pi()/2,TMath::Pi()/2); // Psi
    fContAllChargesMC->SetVarRange(4,-0.5,1.5); // pid mask
    fContAllChargesMC->SetVarName(0,"centrality");
    fContAllChargesMC->SetVarName(1,"charge");
    fContAllChargesMC->SetVarName(2,"prob");
    fContAllChargesMC->SetVarName(3,"#Psi");
    fContAllChargesMC->SetVarName(4,"PIDmask");
    fContAllChargesMC->AddSpecies("all",nPtBinsTOF,binsPtTOF);
    fContAllChargesMC->AddSpecies("pi",nPtBinsTOF,binsPtTOF);
    fContAllChargesMC->AddSpecies("k",nPtBinsTOF,binsPtTOF);
    fContAllChargesMC->AddSpecies("pr",nPtBinsTOF,binsPtTOF);
    fContAllChargesMC->AddSpecies("e",nPtBinsTOF,binsPtTOF);
    fContAllChargesMC->AddSpecies("mu",nPtBinsTOF,binsPtTOF);
    fList3->Add(fContAllChargesMC); 

    fContAllChargesMCA = new AliFlowVZEROResults("v2mcA",5,binsTOFmcPureMC);
    fContAllChargesMCA->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
    fContAllChargesMCA->SetVarRange(1,-1.5,1.5);  // charge
    fContAllChargesMCA->SetVarRange(2,0.6,1.0001);// prob
    fContAllChargesMCA->SetVarRange(3,-TMath::Pi()/2,TMath::Pi()/2); // Psi
    fContAllChargesMCA->SetVarRange(4,-0.5,1.5); // pid mask
    fContAllChargesMCA->SetVarName(0,"centrality");
    fContAllChargesMCA->SetVarName(1,"charge");
    fContAllChargesMCA->SetVarName(2,"prob");
    fContAllChargesMCA->SetVarName(3,"#Psi");
    fContAllChargesMCA->SetVarName(4,"PIDmask");
    fContAllChargesMCA->AddSpecies("all",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCA->AddSpecies("pi",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCA->AddSpecies("k",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCA->AddSpecies("pr",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCA->AddSpecies("e",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCA->AddSpecies("mu",nPtBinsTOF,binsPtTOF);
    fList3->Add(fContAllChargesMCA); 

    fContAllChargesMCC = new AliFlowVZEROResults("v2mcC",5,binsTOFmcPureMC);
    fContAllChargesMCC->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
    fContAllChargesMCC->SetVarRange(1,-1.5,1.5);  // charge
    fContAllChargesMCC->SetVarRange(2,0.6,1.0001);// prob
    fContAllChargesMCC->SetVarRange(3,-TMath::Pi()/2,TMath::Pi()/2); // Psi
    fContAllChargesMCC->SetVarRange(4,-0.5,1.5); // pid mask
    fContAllChargesMCC->SetVarName(0,"centrality");
    fContAllChargesMCC->SetVarName(1,"charge");
    fContAllChargesMCC->SetVarName(2,"prob");
    fContAllChargesMCC->SetVarName(3,"#Psi");
    fContAllChargesMCC->SetVarName(4,"PIDmask");
    fContAllChargesMCC->AddSpecies("all",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCC->AddSpecies("pi",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCC->AddSpecies("k",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCC->AddSpecies("pr",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCC->AddSpecies("e",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCC->AddSpecies("mu",nPtBinsTOF,binsPtTOF);
    fList3->Add(fContAllChargesMCC); 
  }

  // v3 container
  fContAllChargesV0Av3 = new AliFlowVZEROResults("v3A",6,binsTOF);
  fContAllChargesV0Av3->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
  fContAllChargesV0Av3->SetVarRange(1,-1.5,1.5);  // charge
  fContAllChargesV0Av3->SetVarRange(2,0.6,1.0001);// prob
  fContAllChargesV0Av3->SetVarRange(3,-TMath::Pi()/3,TMath::Pi()/3); // Psi
  fContAllChargesV0Av3->SetVarRange(5,-0.5,nDCABin-0.5); // DCA mask
  fContAllChargesV0Av3->SetVarRange(4,-0.5,2.5); // pid mask
  fContAllChargesV0Av3->SetVarName(0,"centrality");
  fContAllChargesV0Av3->SetVarName(1,"charge");
  fContAllChargesV0Av3->SetVarName(2,"prob");
  fContAllChargesV0Av3->SetVarName(3,"#Psi");
  fContAllChargesV0Av3->SetVarName(4,"PIDmask");
  fContAllChargesV0Av3->SetVarName(5,"DCAbin");
  if(fV3) fContAllChargesV0Av3->AddSpecies("all",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Av3->AddSpecies("pi",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Av3->AddSpecies("k",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Av3->AddSpecies("pr",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Av3->AddSpecies("e",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Av3->AddSpecies("d",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Av3->AddSpecies("t",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Av3->AddSpecies("he3",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Av3->AddSpecies("mu",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Av3->AddSpecies("Ks",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Av3->AddSpecies("Lambda",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Av3->AddSpecies("pFromLambda",nPtBinsTOF,binsPtTOF);

  fContAllChargesV0Cv3 = new AliFlowVZEROResults("v3C",6,binsTOF);
  fContAllChargesV0Cv3->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
  fContAllChargesV0Cv3->SetVarRange(1,-1.5,1.5);  // charge
  fContAllChargesV0Cv3->SetVarRange(2,0.6,1.0001);// prob
  fContAllChargesV0Cv3->SetVarRange(3,-TMath::Pi()/3,TMath::Pi()/3); // Psi
  fContAllChargesV0Cv3->SetVarRange(4,-0.5,2.5); // pid mask
  fContAllChargesV0Cv3->SetVarRange(5,-0.5,nDCABin-0.5); // DCA mask
  fContAllChargesV0Cv3->SetVarName(0,"centrality");
  fContAllChargesV0Cv3->SetVarName(1,"charge");
  fContAllChargesV0Cv3->SetVarName(2,"prob");
  fContAllChargesV0Cv3->SetVarName(3,"#Psi");
  fContAllChargesV0Cv3->SetVarName(4,"PIDmask");
  fContAllChargesV0Cv3->SetVarName(5,"DCAbin");
  if(fV3) fContAllChargesV0Cv3->AddSpecies("all",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Cv3->AddSpecies("pi",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Cv3->AddSpecies("k",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Cv3->AddSpecies("pr",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Cv3->AddSpecies("e",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Cv3->AddSpecies("d",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Cv3->AddSpecies("t",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Cv3->AddSpecies("he3",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Cv3->AddSpecies("mu",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Cv3->AddSpecies("Ks",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Cv3->AddSpecies("Lambda",nPtBinsTOF,binsPtTOF);
  if(fV3) fContAllChargesV0Cv3->AddSpecies("pFromLambda",nPtBinsTOF,binsPtTOF);

  fList2->Add(fContAllChargesV0Av3);
  fList2->Add(fContAllChargesV0Cv3);

  if(fIsMC && fV3){
    fContAllChargesMCAv3 = new AliFlowVZEROResults("v3mcA",5,binsTOFmcPureMC);
    fContAllChargesMCAv3->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
    fContAllChargesMCAv3->SetVarRange(1,-1.5,1.5);  // charge
    fContAllChargesMCAv3->SetVarRange(2,0.6,1.0001);// prob
    fContAllChargesMCAv3->SetVarRange(3,-TMath::Pi()/3,TMath::Pi()/3); // Psi
    fContAllChargesMCAv3->SetVarRange(4,-0.5,1.5); // pid mask
    fContAllChargesMCAv3->SetVarName(0,"centrality");
    fContAllChargesMCAv3->SetVarName(1,"charge");
    fContAllChargesMCAv3->SetVarName(2,"prob");
    fContAllChargesMCAv3->SetVarName(3,"#Psi");
    fContAllChargesMCAv3->SetVarName(4,"PIDmask");
    fContAllChargesMCAv3->AddSpecies("all",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCAv3->AddSpecies("pi",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCAv3->AddSpecies("k",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCAv3->AddSpecies("pr",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCAv3->AddSpecies("e",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCAv3->AddSpecies("mu",nPtBinsTOF,binsPtTOF);
    fList3->Add(fContAllChargesMCAv3); 

    fContAllChargesMCCv3 = new AliFlowVZEROResults("v3mcC",5,binsTOFmcPureMC);
    fContAllChargesMCCv3->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
    fContAllChargesMCCv3->SetVarRange(1,-1.5,1.5);  // charge
    fContAllChargesMCCv3->SetVarRange(2,0.6,1.0001);// prob
    fContAllChargesMCCv3->SetVarRange(3,-TMath::Pi()/3,TMath::Pi()/3); // Psi
    fContAllChargesMCCv3->SetVarRange(4,-0.5,1.5); // pid mask
    fContAllChargesMCCv3->SetVarName(0,"centrality");
    fContAllChargesMCCv3->SetVarName(1,"charge");
    fContAllChargesMCCv3->SetVarName(2,"prob");
    fContAllChargesMCCv3->SetVarName(3,"#Psi");
    fContAllChargesMCCv3->SetVarName(4,"PIDmask");
    fContAllChargesMCCv3->AddSpecies("all",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCCv3->AddSpecies("pi",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCCv3->AddSpecies("k",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCCv3->AddSpecies("pr",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCCv3->AddSpecies("e",nPtBinsTOF,binsPtTOF);
    fContAllChargesMCCv3->AddSpecies("mu",nPtBinsTOF,binsPtTOF);
    fList3->Add(fContAllChargesMCCv3); 
  }

  // TProfile for resolutions 3 subevents (V0A, V0C, TPC)
  // v2
  fHResTPCv0A2 = new TProfile("hResTPCv0A2","",nCentrBin,0,nCentrBin);
  fHResTPCv0C2 = new TProfile("hResTPCv0C2","",nCentrBin,0,nCentrBin);
  fHResv0Cv0A2 = new TProfile("hResv0Cv0A2","",nCentrBin,0,nCentrBin);

  fList->Add(fHResTPCv0A2);
  fList->Add(fHResTPCv0C2);
  fList->Add(fHResv0Cv0A2);

  // v3
  fHResTPCv0A3 = new TProfile("hResTPCv0A3","",nCentrBin,0,nCentrBin);
  fHResTPCv0C3 = new TProfile("hResTPCv0C3","",nCentrBin,0,nCentrBin);
  fHResv0Cv0A3 = new TProfile("hResv0Cv0A3","",nCentrBin,0,nCentrBin);

  fList2->Add(fHResTPCv0A3);
  fList2->Add(fHResTPCv0C3);
  fList2->Add(fHResv0Cv0A3);

  // MC as in the dataEP resolution (but using MC tracks)
  if(fIsMC && fV3){
    fHResMA2 = new TProfile("hResMA2","",nCentrBin,0,nCentrBin);
    fHResMC2 = new TProfile("hResMC2","",nCentrBin,0,nCentrBin);
    fHResAC2 = new TProfile("hResAC2","",nCentrBin,0,nCentrBin);
    fList3->Add(fHResMA2); 
    fList3->Add(fHResMC2); 
    fList3->Add(fHResAC2); 
  }
  if(fIsMC && fV3){
    fHResMA3 = new TProfile("hResMA3","",nCentrBin,0,nCentrBin);
    fHResMC3 = new TProfile("hResMC3","",nCentrBin,0,nCentrBin);
    fHResAC3 = new TProfile("hResAC3","",nCentrBin,0,nCentrBin);
    fList3->Add(fHResMA3); 
    fList3->Add(fHResMC3); 
    fList3->Add(fHResAC3); 
  }


  // V0A and V0C event plane distributions
  //v2 
  fPhiRPv0A = new TH2F("fPhiRPv0Av2","#phi distribution of EP VZERO-A;centrality;#phi (rad)",nCentrBin,0,nCentrBin,nPsiTOF,-TMath::Pi()/2,TMath::Pi()/2);
  fPhiRPv0C = new TH2F("fPhiRPv0Cv2","#phi distribution of EP VZERO-C;centrality;#phi (rad)",nCentrBin,0,nCentrBin,nPsiTOF,-TMath::Pi()/2,TMath::Pi()/2);

  //v3
  fPhiRPv0Av3 = new TH2F("fPhiRPv0Av3","#phi distribution of EP VZERO-A;centrality;#phi (rad)",nCentrBin,0,nCentrBin,nPsiTOF,-TMath::Pi()/3,TMath::Pi()/3);
  fPhiRPv0Cv3 = new TH2F("fPhiRPv0Cv3","#phi distribution of EP VZERO-C;centrality;#phi (rad)",nCentrBin,0,nCentrBin,nPsiTOF,-TMath::Pi()/3,TMath::Pi()/3);

  // QA container
  // v2
  const Int_t nDETsignal = 50;
  Double_t binDETsignal[nDETsignal+1];
  for(Int_t i=0;i<nDETsignal+1;i++){
    binDETsignal[i] = -5 + i*10. / nDETsignal;
  }
//   const Int_t nEta = 5;
//   Double_t binEta[nEta+1];
//   for(Int_t i=0;i<nEta+1;i++){
//     binEta[i] = -1 + i*2. / nEta;
//   }

  const Int_t nDeltaPhi = 5;
  const Int_t nDeltaPhiV3 = 7;

  Int_t binsQA[5] = {nCentrTOF,7,5,nDeltaPhi,2};
  Int_t binsQAv3[5] = {nCentrTOF,7,5,nDeltaPhiV3,2};


  fQA = new AliFlowVZEROQA("v2AQA",5,binsQA);
  fQA->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
  fQA->SetVarRange(1,0,7);  // pt
  fQA->SetVarRange(2,0.,1.0001);// prob
  fQA->SetVarRange(3,-TMath::Pi(),TMath::Pi()); // Psi
  fQA->SetVarRange(4,-0.5,1.5); // pid mask
  fQA->SetVarName(0,"centrality");
  fQA->SetVarName(1,"p_{t}");
  fQA->SetVarName(2,"prob");
  fQA->SetVarName(3,"#Psi");
  fQA->SetVarName(4,"PIDmask");
  if(fQAsw && fV2) fQA->AddSpecies("pi",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
  if(fQAsw && fV2) fQA->AddSpecies("k",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
  if(fQAsw && fV2) fQA->AddSpecies("pr",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   if(fQAsw && fV2) fQA->AddSpecies("e",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   fQA->AddSpecies("d",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   fQA->AddSpecies("t",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   fQA->AddSpecies("he3",nDETsignal,binDETsignal,nDETsignal,binDETsignal);

  fQA2 = new AliFlowVZEROQA("v2CQA",5,binsQA);
  fQA2->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
  fQA2->SetVarRange(1,0,7);  // pt
  fQA2->SetVarRange(2,0.,1.0001);// prob
  fQA2->SetVarRange(3,-TMath::Pi(),TMath::Pi()); // Psi
  fQA2->SetVarRange(4,-0.5,1.5); // pid mask
  fQA2->SetVarName(0,"centrality");
  fQA2->SetVarName(1,"p_{t}");
  fQA2->SetVarName(2,"prob");
  fQA2->SetVarName(3,"#Psi");
  fQA2->SetVarName(4,"PIDmask");
  if(fQAsw && fV2) fQA2->AddSpecies("pi",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
  if(fQAsw && fV2) fQA2->AddSpecies("k",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
  if(fQAsw && fV2) fQA2->AddSpecies("pr",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   if(fQAsw && fV2) fQA2->AddSpecies("e",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   fQA2->AddSpecies("d",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   fQA2->AddSpecies("t",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   fQA2->AddSpecies("he3",nDETsignal,binDETsignal,nDETsignal,binDETsignal);

  fQAv3 = new AliFlowVZEROQA("v3AQA",5,binsQAv3);
  fQAv3->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
  fQAv3->SetVarRange(1,0,7);  // pt
  fQAv3->SetVarRange(2,0.,1.0001);// prob
  fQAv3->SetVarRange(3,-TMath::Pi(),TMath::Pi()); // Psi
  fQAv3->SetVarRange(4,-0.5,1.5); // pid mask
  fQAv3->SetVarName(0,"centrality");
  fQAv3->SetVarName(1,"p_{t}");
  fQAv3->SetVarName(2,"prob");
  fQAv3->SetVarName(3,"#Psi");
  fQAv3->SetVarName(4,"PIDmask");
  if(fQAsw && fV3) fQAv3->AddSpecies("pi",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   if(fQAsw && fV3) fQAv3->AddSpecies("k",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   if(fQAsw && fV3) fQAv3->AddSpecies("pr",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   if(fQAsw && fV2) fQAv3->AddSpecies("e",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   fQAv3->AddSpecies("d",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   fQAv3->AddSpecies("t",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   fQAv3->AddSpecies("he3",nDETsignal,binDETsignal,nDETsignal,binDETsignal);

  fQA2v3 = new AliFlowVZEROQA("v3CQA",5,binsQAv3);
  fQA2v3->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
  fQA2v3->SetVarRange(1,0,7);  // pt
  fQA2v3->SetVarRange(2,0.,1.0001);// prob
  fQA2v3->SetVarRange(3,-TMath::Pi(),TMath::Pi()); // Psi
  fQA2v3->SetVarRange(4,-0.5,1.5); // pid mask
  fQA2v3->SetVarName(0,"centrality");
  fQA2v3->SetVarName(1,"p_{t}");
  fQA2v3->SetVarName(2,"prob");
  fQA2v3->SetVarName(3,"#Psi");
  fQA2v3->SetVarName(4,"PIDmask");
  if(fQAsw && fV3) fQA2v3->AddSpecies("pi",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   if(fQAsw && fV3) fQA2v3->AddSpecies("k",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   if(fQAsw && fV3) fQA2v3->AddSpecies("pr",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   if(fQAsw && fV2) fQA2v3->AddSpecies("e",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   fQA2v3->AddSpecies("d",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   fQA2v3->AddSpecies("t",nDETsignal,binDETsignal,nDETsignal,binDETsignal);
//   fQA2v3->AddSpecies("he3",nDETsignal,binDETsignal,nDETsignal,binDETsignal);

  fList->Add(fPhiRPv0A);
  fList->Add(fPhiRPv0C);

  if(fQAsw && fV2){
    fList4->Add(fQA);
    fList4->Add(fQA2);
  }

  fList2->Add(fPhiRPv0Av3);
  fList2->Add(fPhiRPv0Cv3);

  if(fQAsw && fV3){
   fList4->Add(fQAv3);
   fList4->Add(fQA2v3);
  }

  //  fList->Add(fTree); // comment if not needed

  const Int_t nDCA = 300;
  Double_t DCAbin[nDCA+1];
  for(Int_t i=0;i <= nDCA;i++){
    DCAbin[i] = -3 +i*6.0/nDCA;
  }
  
  char nameHistos[100];
  for(Int_t iC=0;iC < nCentrBin;iC++){
    snprintf(nameHistos,100,"fHdcaPtPiCent%i",iC);
    fHdcaPt[iC][0] = new TH2D(nameHistos,"DCA_{xy} for #pi;p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
    snprintf(nameHistos,100,"fHdcaPtKaCent%i",iC);
    fHdcaPt[iC][1] = new TH2D(nameHistos,"DCA_{xy} for K;p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
    snprintf(nameHistos,100,"fHdcaPtPrCent%i",iC);
    fHdcaPt[iC][2] = new TH2D(nameHistos,"DCA_{xy} for #bar{p};p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
    snprintf(nameHistos,100,"fHdcaPtElCent%i",iC);
    fHdcaPt[iC][3] = new TH2D(nameHistos,"DCA_{xy} for e;p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
    snprintf(nameHistos,100,"fHdcaPtDeCent%i",iC);
    fHdcaPt[iC][4] = new TH2D(nameHistos,"DCA_{xy} for #bar{d};p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
    snprintf(nameHistos,100,"fHdcaPtTrCent%i",iC);
    fHdcaPt[iC][5] = new TH2D(nameHistos,"DCA_{xy} for #bar{t};p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
    snprintf(nameHistos,100,"fHdcaPtHeCent%i",iC);
    fHdcaPt[iC][6] = new TH2D(nameHistos,"DCA_{xy} for #bar{^{3}He};p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
  }
  
  if(fFillDCA && fQAsw){
    for(Int_t i=0;i<7;i++)
      for(Int_t iC=0;iC < nCentrBin;iC++)
	fList4->Add(fHdcaPt[iC][i]);
  }
  if(fIsMC){
    for(Int_t iC=0;iC < nCentrBin;iC++){
      snprintf(nameHistos,100,"fHdcaPtPiSecCent%i",iC);
      fHdcaPtSec[iC][0] = new TH2D(nameHistos,"DCA_{xy} for secondary #pi;p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
      snprintf(nameHistos,100,"fHdcaPtKaSecCent%i",iC);
      fHdcaPtSec[iC][1] = new TH2D(nameHistos,"DCA_{xy} for secondary K;p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
      snprintf(nameHistos,100,"fHdcaPtPrSecCent%i",iC);
      fHdcaPtSec[iC][2] = new TH2D(nameHistos,"DCA_{xy} for secondary #bar{p};p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
      snprintf(nameHistos,100,"fHdcaPtElSecCent%i",iC);
      fHdcaPtSec[iC][3] = new TH2D(nameHistos,"DCA_{xy} for secondary e;p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
      snprintf(nameHistos,100,"fHdcaPtDeSecCent%i",iC);
      fHdcaPtSec[iC][4] = new TH2D(nameHistos,"DCA_{xy} for secondary #bar{d};p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
      snprintf(nameHistos,100,"fHdcaPtTrSecCent%i",iC);
      fHdcaPtSec[iC][5] = new TH2D(nameHistos,"DCA_{xy} for secondary #bar{t};p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
      snprintf(nameHistos,100,"fHdcaPtHeSecCent%i",iC);
      fHdcaPtSec[iC][6] = new TH2D(nameHistos,"DCA_{xy} for secondary #bar{^{3}He};p_{t} (GeV/c);DCA_{xy} (cm)",nPtBinsTOF,binsPtTOF,nDCA,DCAbin);
    }
    
    if(fFillDCA && fQAsw){
      for(Int_t i=0;i<7;i++)
	for(Int_t iC=0;iC < nCentrBin;iC++)
	  fList4->Add(fHdcaPtSec[iC][i]);
    }
  }
  
  // Add TProfile Extra QA
  const Int_t nBinQApid = 2;
  Int_t binQApid[nBinQApid] = {nCentrTOF,200};
  const Int_t nbinsigma = 100;
  Double_t nsigmaQA[nbinsigma+1];
  for(Int_t i=0;i<nbinsigma+1;i++){
    nsigmaQA[i] = -10 + 20.0*i/nbinsigma;
  }
  fContQApid = new AliFlowVZEROResults("qaPID",nBinQApid,binQApid);
  fContQApid->SetVarRange(0,-0.5,nCentrBin-0.5); // centrality
  fContQApid->SetVarRange(1,0,20);  // charge
  fContQApid->SetVarName(0,"centrality");
  fContQApid->SetVarName(1,"p_{t}");
  fContQApid->AddSpecies("piTPC",nbinsigma,nsigmaQA);
  fContQApid->AddSpecies("piTOF",nbinsigma,nsigmaQA);
  fContQApid->AddSpecies("kaTPC",nbinsigma,nsigmaQA);
  fContQApid->AddSpecies("kaTOF",nbinsigma,nsigmaQA);
  fContQApid->AddSpecies("prTPC",nbinsigma,nsigmaQA);
  fContQApid->AddSpecies("prTOF",nbinsigma,nsigmaQA);
  if(fV2) fList->Add(fContQApid);
  
  printf("Output creation ok!!\n\n\n\n");

  // Post output data.
  if(fV2) PostData(1, fList);
  if(fV3) PostData(2, fList2);
  if(fIsMC) PostData(3, fList3);
  if(fQAsw) PostData(4, fList4);

}

//______________________________________________________________________________
void AliAnalysisTaskVnV0::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
    
    fgIsPsiComputed = kFALSE;
    fgPsi2v0a=999.;
    fgPsi2v0c=999.;
    fgPsi2tpc=999.;
    fgPsi3v0a=999.;
    fgPsi3v0c=999.;
    fgPsi3tpc=999.;
    fgPsi2v0aMC=999.;
    fgPsi2v0cMC=999.;
    fgPsi2tpcMC=999.;
    fgPsi3v0aMC=999.;
    fgPsi3v0cMC=999.;
    fgPsi3tpcMC=999.;

    fOutputAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fOutputAOD){
	Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
	this->Dump();
	return;
    }
    
    Int_t run = fOutputAOD->GetRunNumber();

    if(run != fRun){
	// Load the calibrations run dependent
	OpenInfoCalbration(run);
	fRun=run;
    }

    Float_t zvtx = GetVertex(fOutputAOD);



    //Get the MC object
    if(fIsMC){
      AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(fOutputAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
      if (!mcHeader) {
	AliError("Could not find MC Header in AOD");
	return;
      }
    }

    /*
    AliMCEvent* mcEvent = MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
    
    Double_t gReactionPlane = -999., gImpactParameter = -999.;
    //Get the MC header
    AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());
    if (headerH) {
      //Printf("=====================================================");
      //Printf("Reaction plane angle: %lf",headerH->ReactionPlaneAngle());
      //Printf("=====================================================");
      gReactionPlane = headerH->ReactionPlaneAngle();
      gImpactParameter = headerH->ImpactParameter();
    }

*/

    //    printf("vertex = %f\n",zvtx);
    if (TMath::Abs(zvtx) < fVtxCut) {
      //Centrality
      Float_t v0Centr  = -10.;
      Float_t trkCentr  = -10.;
      AliCentrality *centrality = fOutputAOD->GetCentrality();
      if (centrality){
//	printf("v0centr = %f -- tpccnetr%f\n",centrality->GetCentralityPercentile("V0M"),centrality->GetCentralityPercentile("TRK"));
	v0Centr  = centrality->GetCentralityPercentile("V0M");
	trkCentr = centrality->GetCentralityPercentile("TRK"); 
      }

      if(TMath::Abs(v0Centr - trkCentr) < 5.0){ // consistency cut on centrality selection
	fPID->SetDetResponse(fOutputAOD, v0Centr); // Set the PID object for each event!!!!
	Analyze(fOutputAOD,v0Centr); // Do analysis!!!!

        fCentrality = v0Centr;
        if(fV2) fTree->Fill();
      }
    }
    
}

//________________________________________________________________________
void AliAnalysisTaskVnV0::Analyze(AliAODEvent* aodEvent, Float_t v0Centr)
{      
  Float_t mass[8] = {5.10998909999999971e-04, 1.05658000000000002e-01, 1.39570000000000000e-01, 4.93676999999999977e-01, 9.38271999999999995e-01,1.87783699999999998,2.81740199999999996,1.40805449999999999};
  
  // Event plane resolution for v2
  Float_t evPlRes[18] = {0.350582,0.505393,0.607845,0.632913,0.592230,0.502489,0.381717,0.249539,0.133180, // V0A vs. centrality
			 0.446480,0.612705,0.712222,0.736200,0.697907,0.610114,0.481009,0.327402,0.182277};// V0C vs. centrality
  
  Int_t iC = -1;    
  if (v0Centr >0 && v0Centr < 80){ // analysis only for 0-80% centrality classes

    fgIsPsiComputed = kTRUE;

    // centrality bins
    if(v0Centr < 5) iC = 0;
    else if(v0Centr < 10) iC = 1;
    else if(v0Centr < 20) iC = 2;
    else if(v0Centr < 30) iC = 3;
    else if(v0Centr < 40) iC = 4;
    else if(v0Centr < 50) iC = 5;
    else if(v0Centr < 60) iC = 6;
    else if(v0Centr < 70) iC = 7;
    else iC = 8;

    Int_t iCcal = iC;

    if(nCentrBin==16){
      iC = Int_t(v0Centr/5);
       if(iC >= nCentrBin) iC = nCentrBin-1;
    }
    
    //reset Q vector info	
    Double_t Qxa2 = 0, Qya2 = 0;
    Double_t Qxc2 = 0, Qyc2 = 0;
    Double_t Qxa3 = 0, Qya3 = 0;
    Double_t Qxc3 = 0, Qyc3 = 0;

    Int_t nAODTracks = aodEvent->GetNumberOfTracks();

    AliAODMCHeader *mcHeader = NULL;
    TClonesArray *mcArray = NULL;
    Float_t evplaneMC = 0;
    if(fIsMC){
      mcHeader = dynamic_cast<AliAODMCHeader*>(fOutputAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));

      if (mcHeader) {	
	evplaneMC = mcHeader->GetReactionPlaneAngle();
	if(evplaneMC > TMath::Pi()/2 && evplaneMC <=  TMath::Pi()*3/2) evplaneMC-=TMath::Pi(); 
	else if(evplaneMC > TMath::Pi()*3/2) evplaneMC-=2*TMath::Pi(); 
	mcArray = (TClonesArray*)fOutputAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
      
	if(mcArray){
	  Float_t QxMCv2[3] = {0,0,0};
	  Float_t QyMCv2[3] = {0,0,0};
	  Float_t QxMCv3[3] = {0,0,0};
	  Float_t QyMCv3[3] = {0,0,0};
	  Float_t EvPlaneMCV2[3] = {0,0,0};
	  Float_t EvPlaneMCV3[3] = {0,0,0};
	  Float_t etaMin[3] = {2.8,-3.6,-0.8}; // A-side, C-side M-barrel
	  Float_t etaMax[3] = {4.88,-1.8,0.8};

	  // analysis on MC tracks
	  Int_t nMCtrack = mcArray->GetEntries() ;

	  // EP computation with MC tracks
	  for(Int_t iT=0;iT < nMCtrack;iT++){
	    AliAODMCParticle *mctr = (AliAODMCParticle*) mcArray->At(iT);
	    if(!mctr || !(mctr->IsPrimary()) || !(mctr->Charge()) || mctr->Pt() < 0.2) continue;
	    
	    Float_t eta = mctr->Eta();
	    
	    for(Int_t iD=0;iD<3;iD++){
	      if(eta > etaMin[iD] && eta < etaMax[iD]){
		Float_t phi = mctr->Phi();
		QxMCv2[iD] += TMath::Cos(2*phi);
		QyMCv2[iD] += TMath::Sin(2*phi);
		QxMCv3[iD] += TMath::Cos(3*phi);
		QyMCv3[iD] += TMath::Sin(3*phi);
	      }
	    }
	  }

	  if(fV2){
	    EvPlaneMCV2[0] = TMath::ATan2(QyMCv2[0],QxMCv2[0])/2.;
	    EvPlaneMCV2[1] = TMath::ATan2(QyMCv2[1],QxMCv2[1])/2.;
	    EvPlaneMCV2[2] = TMath::ATan2(QyMCv2[2],QxMCv2[2])/2.;
	    fHResMA2->Fill(Double_t(iC), TMath::Cos(2*(EvPlaneMCV2[2]-EvPlaneMCV2[0])));
	    fHResMC2->Fill(Double_t(iC), TMath::Cos(2*(EvPlaneMCV2[2]-EvPlaneMCV2[1])));
	    fHResAC2->Fill(Double_t(iC), TMath::Cos(2*(EvPlaneMCV2[0]-EvPlaneMCV2[1])));
            fgPsi2v0aMC = EvPlaneMCV2[0];
            fgPsi2v0cMC = EvPlaneMCV2[1];
            fgPsi2tpcMC = EvPlaneMCV2[2];
	  }
	  if(fV3){
	    EvPlaneMCV3[0] = TMath::ATan2(QyMCv3[0],QxMCv3[0])/3.;
	    EvPlaneMCV3[1] = TMath::ATan2(QyMCv3[1],QxMCv3[1])/3.;
	    EvPlaneMCV3[2] = TMath::ATan2(QyMCv3[2],QxMCv3[2])/3.;
	    fHResMA3->Fill(Double_t(iC), TMath::Cos(3*(EvPlaneMCV3[2]-EvPlaneMCV3[0])));
	    fHResMC3->Fill(Double_t(iC), TMath::Cos(3*(EvPlaneMCV3[2]-EvPlaneMCV3[1])));
	    fHResAC3->Fill(Double_t(iC), TMath::Cos(3*(EvPlaneMCV3[0]-EvPlaneMCV3[1])));
            fgPsi3v0aMC = EvPlaneMCV3[0];
            fgPsi3v0cMC = EvPlaneMCV3[1];
            fgPsi3tpcMC = EvPlaneMCV3[2];
	  }

	  // flow A and C side
	  Float_t xMCepAv2[5] = {Float_t(iC),0/*charge*/,1,EvPlaneMCV2[0],1};
	  Float_t xMCepCv2[5] = {Float_t(iC),0/*charge*/,1,EvPlaneMCV2[1],1};
	  Float_t xMCepAv3[5] = {Float_t(iC),0/*charge*/,1,EvPlaneMCV3[0],1};
	  Float_t xMCepCv3[5] = {Float_t(iC),0/*charge*/,1,EvPlaneMCV3[1],1};
	  
	  for(Int_t iT=0;iT < nMCtrack;iT++){
	    AliAODMCParticle *mctr = (AliAODMCParticle*) mcArray->At(iT);
	    if(!mctr || !(mctr->IsPhysicalPrimary()) || !(mctr->Charge()) || TMath::Abs(mctr->Eta()) > 0.8 || mctr->Pt() < 0.2) continue;
	    Int_t iS = TMath::Abs(mctr->GetPdgCode());
	    Int_t charge = mctr->Charge();
	    Float_t pt = mctr->Pt();
	    Float_t phi = mctr->Phi();

	    if(charge > 0){
	      xMCepAv2[1] = 1;
	      xMCepCv2[1] = 1;
	      xMCepAv3[1] = 1;
	      xMCepCv3[1] = 1;
	    }
	    else{
	      xMCepAv2[1] = -1;
	      xMCepCv2[1] = -1;
	      xMCepAv3[1] = -1;
	      xMCepCv3[1] = -1;
	    }

	    fContAllChargesMCA->Fill(0,pt, TMath::Cos(2*(phi - EvPlaneMCV2[0])),xMCepAv2);
	    fContAllChargesMCC->Fill(0,pt, TMath::Cos(2*(phi - EvPlaneMCV2[1])),xMCepCv2);
	    fContAllChargesMCAv3->Fill(0,pt, TMath::Cos(3*(phi - EvPlaneMCV3[0])),xMCepAv3);
	    fContAllChargesMCCv3->Fill(0,pt, TMath::Cos(3*(phi - EvPlaneMCV3[1])),xMCepCv3);

	    if(iS==11){
	      fContAllChargesMCA->Fill(4,pt, TMath::Cos(2*(phi - EvPlaneMCV2[0])),xMCepAv2);
	      fContAllChargesMCC->Fill(4,pt, TMath::Cos(2*(phi - EvPlaneMCV2[1])),xMCepCv2);
	      fContAllChargesMCAv3->Fill(4,pt, TMath::Cos(3*(phi - EvPlaneMCV3[0])),xMCepAv3);
	      fContAllChargesMCCv3->Fill(4,pt, TMath::Cos(3*(phi - EvPlaneMCV3[1])),xMCepCv3);
	    }
	    else if(iS==13){
	      fContAllChargesMCA->Fill(5,pt, TMath::Cos(2*(phi - EvPlaneMCV2[0])),xMCepAv2);
	      fContAllChargesMCC->Fill(5,pt, TMath::Cos(2*(phi - EvPlaneMCV2[1])),xMCepCv2);
	      fContAllChargesMCAv3->Fill(5,pt, TMath::Cos(3*(phi - EvPlaneMCV3[0])),xMCepAv3);
	      fContAllChargesMCCv3->Fill(5,pt, TMath::Cos(3*(phi - EvPlaneMCV3[1])),xMCepCv3);
	    }
	    else if(iS==211){
	      fContAllChargesMCA->Fill(1,pt, TMath::Cos(2*(phi - EvPlaneMCV2[0])),xMCepAv2);
	      fContAllChargesMCC->Fill(1,pt, TMath::Cos(2*(phi - EvPlaneMCV2[1])),xMCepCv2);
	      fContAllChargesMCAv3->Fill(1,pt, TMath::Cos(3*(phi - EvPlaneMCV3[0])),xMCepAv3);
	      fContAllChargesMCCv3->Fill(1,pt, TMath::Cos(3*(phi - EvPlaneMCV3[1])),xMCepCv3);
	    }
	    else if(iS==321){
	      fContAllChargesMCA->Fill(2,pt, TMath::Cos(2*(phi - EvPlaneMCV2[0])),xMCepAv2);
	      fContAllChargesMCC->Fill(2,pt, TMath::Cos(2*(phi - EvPlaneMCV2[1])),xMCepCv2);
	      fContAllChargesMCAv3->Fill(2,pt, TMath::Cos(3*(phi - EvPlaneMCV3[0])),xMCepAv3);
	      fContAllChargesMCCv3->Fill(2,pt, TMath::Cos(3*(phi - EvPlaneMCV3[1])),xMCepCv3);
	    }
	    else if(iS==2212){
	      fContAllChargesMCA->Fill(3,pt, TMath::Cos(2*(phi - EvPlaneMCV2[0])),xMCepAv2);
	      fContAllChargesMCC->Fill(3,pt, TMath::Cos(2*(phi - EvPlaneMCV2[1])),xMCepCv2);
	      fContAllChargesMCAv3->Fill(3,pt, TMath::Cos(3*(phi - EvPlaneMCV3[0])),xMCepAv3);
	      fContAllChargesMCCv3->Fill(3,pt, TMath::Cos(3*(phi - EvPlaneMCV3[1])),xMCepCv3);
	    }
	  }
	}
      }
    }

    //V0 info    
    AliAODVZERO* aodV0 = aodEvent->GetVZEROData();

    for (Int_t iv0 = 0; iv0 < 64; iv0++) {
      Double_t phiV0 = TMath::PiOver4()*(0.5 + iv0 % 8);
      Float_t multv0 = aodV0->GetMultiplicity(iv0);
      if (iv0 < 32){ // V0C
	Qxc2 += TMath::Cos(2*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	Qyc2 += TMath::Sin(2*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	Qxc3 += TMath::Cos(3*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	Qyc3 += TMath::Sin(3*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
      } else {       // V0A
	Qxa2 += TMath::Cos(2*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
	Qya2 += TMath::Sin(2*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
	Qxa3 += TMath::Cos(3*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
	Qya3 += TMath::Sin(3*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
      }
    }

    //grab for each centrality the proper histo with the Qx and Qy to do the recentering
    Double_t Qxamean2 = fMeanQ[iCcal][1][0];
    Double_t Qxarms2  = fWidthQ[iCcal][1][0];
    Double_t Qyamean2 = fMeanQ[iCcal][1][1];
    Double_t Qyarms2  = fWidthQ[iCcal][1][1];
    Double_t Qxamean3 = fMeanQv3[iCcal][1][0];
    Double_t Qxarms3  = fWidthQv3[iCcal][1][0];
    Double_t Qyamean3 = fMeanQv3[iCcal][1][1];
    Double_t Qyarms3  = fWidthQv3[iCcal][1][1];
    
    Double_t Qxcmean2 = fMeanQ[iCcal][0][0];
    Double_t Qxcrms2  = fWidthQ[iCcal][0][0];
    Double_t Qycmean2 = fMeanQ[iCcal][0][1];
    Double_t Qycrms2  = fWidthQ[iCcal][0][1];	
    Double_t Qxcmean3 = fMeanQv3[iCcal][0][0];
    Double_t Qxcrms3  = fWidthQv3[iCcal][0][0];
    Double_t Qycmean3 = fMeanQv3[iCcal][0][1];
    Double_t Qycrms3  = fWidthQv3[iCcal][0][1];	
    
    Double_t QxaCor2 = (Qxa2 - Qxamean2)/Qxarms2;
    Double_t QyaCor2 = (Qya2 - Qyamean2)/Qyarms2;
    Double_t QxcCor2 = (Qxc2 - Qxcmean2)/Qxcrms2;
    Double_t QycCor2 = (Qyc2 - Qycmean2)/Qycrms2;
    Double_t QxaCor3 = (Qxa3 - Qxamean3)/Qxarms3;
    Double_t QyaCor3 = (Qya3 - Qyamean3)/Qyarms3;
    Double_t QxcCor3 = (Qxc3 - Qxcmean3)/Qxcrms3;
    Double_t QycCor3 = (Qyc3 - Qycmean3)/Qycrms3;
	
    evPlAngV0ACor2 = TMath::ATan2(QyaCor2, QxaCor2)/2.;
    evPlAngV0CCor2 = TMath::ATan2(QycCor2, QxcCor2)/2.;
    evPlAngV0ACor3 = TMath::ATan2(QyaCor3, QxaCor3)/3.;
    evPlAngV0CCor3 = TMath::ATan2(QycCor3, QxcCor3)/3.;

    fgPsi2v0a = evPlAngV0ACor2;
    fgPsi2v0c = evPlAngV0CCor2;
    fgPsi3v0a = evPlAngV0ACor3;
    fgPsi3v0c = evPlAngV0CCor3;
				 
    //loop track and get pid
    for(Int_t iT = 0; iT < nAODTracks; iT++) { // loop on the tracks
      AliAODTrack* aodTrack = aodEvent->GetTrack(iT);
	
      if (!aodTrack){
	aodTrack->Delete();
	continue;
      }
      
      Bool_t trkFlag = aodTrack->TestFilterBit(1); // TPC only tracks
      if(fFillDCA) trkFlag = aodTrack->TestFilterBit(4); // Global track, DCA loose cut

      if ((TMath::Abs(aodTrack->Eta()) > fEtaCut) || (aodTrack->Pt() < fMinPt) || (aodTrack->GetTPCNcls() < fNcluster) || !trkFlag){
	continue;
      }

      Double_t b[2] = {-99., -99.};
      Double_t bCov[3] = {-99., -99., -99.};
      if (!aodTrack->PropagateToDCA(fOutputAOD->GetPrimaryVertex(), fOutputAOD->GetMagneticField(), 100., b, bCov))
	continue;
	    
      if (!fFillDCA && ((TMath::Abs(b[0]) > 3.0) || (TMath::Abs(b[1]) > 2.4)))
	continue;
 	    
      if(fFillDCA && TMath::Abs(b[0]) > 3.0 && TMath::Abs(b[1]) > 3)
	continue;
      
      // re-map the container in an array to do the analysis for V0A and V0C within a loop
      Float_t evPlAngV0[2] = {evPlAngV0ACor2,evPlAngV0CCor2};
      AliFlowVZEROResults *contV0[2] = {fContAllChargesV0A,fContAllChargesV0C};
      AliFlowVZEROQA *QA[2] = {fQA,fQA2};

      Float_t evPlAngV0v3[2] = {evPlAngV0ACor3,evPlAngV0CCor3};
      AliFlowVZEROResults *contV0v3[2] = {fContAllChargesV0Av3,fContAllChargesV0Cv3};
      AliFlowVZEROQA *QAv3[2] = {fQAv3,fQA2v3};

      // Fill MC results
      if(fIsMC && mcArray){
	fPID->ComputeProb(aodTrack,fOutputAOD); // compute Bayesian probabilities
	Float_t tofMismProbMC = fPID->GetTOFMismProb(); // TOF mismatch probability requested to be lower than 50% for TOF analysis 

	Float_t xMC[5] = {Float_t(iC),Float_t(aodTrack->Charge()),1,evplaneMC,Float_t(fPID->GetCurrentMask(1)&&tofMismProbMC < 0.5)}; // to fill analysis v2 container

	Float_t v2mc = TMath::Cos(2*(aodTrack->Phi() - evplaneMC));

	fContAllChargesMC->Fill(0,aodTrack->Pt(),v2mc,xMC);
	
	Int_t iS = TMath::Abs(((AliAODMCParticle*)mcArray->At(TMath::Abs(aodTrack->GetLabel())))->GetPdgCode());
	if(iS==11){
	  fContAllChargesMC->Fill(4,aodTrack->Pt(),v2mc,xMC);
	}
	else if(iS==13){
	  fContAllChargesMC->Fill(5,aodTrack->Pt(),v2mc,xMC);	  
	}
	else if(iS==211){
	  fContAllChargesMC->Fill(1,aodTrack->Pt(),v2mc,xMC);
	}
	else if(iS==321){
	  fContAllChargesMC->Fill(2,aodTrack->Pt(),v2mc,xMC);
	}
	else if(iS==2212){
	  fContAllChargesMC->Fill(3,aodTrack->Pt(),v2mc,xMC);	  
	}
      }

      for(Int_t iV0=0;iV0<2;iV0++){ // loop on A and C side

	if(fModulationDEDx) fPID->SetPsiCorrectionDeDx(evPlAngV0[iV0],evPlRes[iV0*8+iC]); // set the PID dE/dx correction as a function of the v2-EP (resolution is needed)

	Float_t v2V0 = TMath::Cos(2*(aodTrack->Phi() - evPlAngV0[iV0]));
	Float_t v3V0 = TMath::Cos(3*(aodTrack->Phi() - evPlAngV0v3[iV0]));
	    
	fPID->ComputeProb(aodTrack,fOutputAOD); // compute Bayesian probabilities
	Float_t dedx = fPID->GetDeDx();//aodTrack->GetTPCsignal();
	Float_t *probRead = fPID->GetProb();
	Float_t prob[8] = {probRead[0],probRead[1],probRead[2],probRead[3],probRead[4],probRead[5],probRead[6],probRead[7]};
	Float_t tofMismProb = fPID->GetTOFMismProb(); // TOF mismatch probability requested to be lower than 50% for TOF analysis 
	Float_t x[6] = {Float_t(iC),Float_t(aodTrack->Charge()),1,evPlAngV0[iV0],Float_t(fPID->GetCurrentMask(1)&&tofMismProb < 0.5),0}; // to fill analysis v2 container
	Float_t x3[6] = {Float_t(iC),Float_t(aodTrack->Charge()),1,evPlAngV0v3[iV0],Float_t(fPID->GetCurrentMask(1)&&tofMismProb < 0.5),0}; // to fill analysis v3 container

	// in case fill DCA info
	if(fFillDCA){
	  if(TMath::Abs(b[0]) > 0.1){
	    x[5] = 1;
	    x3[5] = 1;
	  }
	  if(TMath::Abs(b[0]) > 0.3){
	    x[5] = 2;
	    x3[5] = 2;
	  }
	  if(fIsMC && mcArray){
	    if(!((AliAODMCParticle*)mcArray->At(TMath::Abs(aodTrack->GetLabel())))->IsPhysicalPrimary()){
	      x[5] += 3;
	      x3[5] += 3;
	    }
	  }
	}

	// Fill no PID
	if(fV2) contV0[iV0]->Fill(0,aodTrack->Pt(),v2V0,x);
	if(fV3) contV0v3[iV0]->Fill(0,aodTrack->Pt(),v3V0,x3);


	Double_t dedxExp[8];
	Float_t tof = -1;
	Double_t inttimes[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
	Double_t expTOFsigma[8] = {0.,0.,0.,0.,0.,0.,0.,0.};

	Float_t nsigmaTPC[8];
	Float_t nsigmaTOF[8];

	if(aodTrack->GetDetPid()){ // check the PID object is available
	  for(Int_t iS=0;iS < 8;iS++){
	    dedxExp[iS] = fPID->GetExpDeDx(aodTrack,iS);
	    nsigmaTPC[iS] = (dedx - dedxExp[iS])/(dedxExp[iS]*0.07);
	    //	    printf("TPC %i = %f (%f %f)\n",iS, nsigmaTPC[iS],dedx, dedxExp[iS]);
	  }
	  		
	  if(fPID->GetCurrentMask(1)){ // if TOF is present
	    Float_t ptrack = aodTrack->P();
	    tof = aodTrack->GetTOFsignal() - fPID->GetESDpid()->GetTOFResponse().GetStartTime(ptrack);
	    aodTrack->GetIntegratedTimes(inttimes);
	    
	    for(Int_t iS=5;iS < 8;iS++) // extra info for light nuclei
	      inttimes[iS] = inttimes[0] / ptrack * mass[iS] * TMath::Sqrt(1+ptrack*ptrack/mass[iS]/mass[iS]);
	    
	    for(Int_t iS=0;iS<8;iS++){
	      expTOFsigma[iS] = fPID->GetESDpid()->GetTOFResponse().GetExpectedSigma(ptrack, inttimes[iS], mass[iS]);
	      nsigmaTOF[iS] = (tof - inttimes[iS])/expTOFsigma[iS];
	      //	      printf("TOF %i = %f\n",iS, nsigmaTOF[iS]);
	    }
	  }
	}

	Float_t deltaPhiV0 = aodTrack->Phi() - evPlAngV0[iV0];
	if(deltaPhiV0 > TMath::Pi()) deltaPhiV0 -= 2*TMath::Pi();
	else if(deltaPhiV0 < -TMath::Pi()) deltaPhiV0 += 2*TMath::Pi();
	if(deltaPhiV0 > TMath::Pi()) deltaPhiV0 -= 2*TMath::Pi();
	else if(deltaPhiV0 < -TMath::Pi()) deltaPhiV0 += 2*TMath::Pi();
	
	Float_t deltaPhiV0v3 = aodTrack->Phi() - evPlAngV0v3[iV0];
	if(deltaPhiV0v3 > TMath::Pi()) deltaPhiV0v3 -= 2*TMath::Pi();
	else if(deltaPhiV0v3 < -TMath::Pi()) deltaPhiV0v3 += 2*TMath::Pi();
	if(deltaPhiV0v3 > TMath::Pi()) deltaPhiV0v3 -= 2*TMath::Pi();
	else if(deltaPhiV0v3 < -TMath::Pi()) deltaPhiV0v3 += 2*TMath::Pi();

	// variable to fill QA container
	Float_t xQA[5] = {Float_t(iC),Float_t(aodTrack->Pt()), 0.0,deltaPhiV0,x[4]}; // v2
	Float_t xQA3[5] = {Float_t(iC),Float_t(aodTrack->Pt()), 0.0,deltaPhiV0v3,x[4]}; // v3

	// extra QA TProfiles
	if(iV0==1 && aodTrack->Pt() < 20 && fPID->GetCurrentMask(0) && fPID->GetCurrentMask(1)){
	  Float_t xQApid[2] = {Float_t(iC),Float_t(aodTrack->Pt())};
	  fContQApid->Fill(0,nsigmaTPC[2],v2V0,xQApid); // v2 TPC (V0C) w.r.t pions
	  fContQApid->Fill(1,nsigmaTOF[2],v2V0,xQApid); // v2 TOF (V0C) w.r.t. pions
	  fContQApid->Fill(2,nsigmaTPC[3],v2V0,xQApid); // v2 TPC (V0C) w.r.t kaons
	  fContQApid->Fill(3,nsigmaTOF[3],v2V0,xQApid); // v2 TOF (V0C) w.r.t. kaons
	  fContQApid->Fill(4,nsigmaTPC[4],v2V0,xQApid); // v2 TPC (V0C) w.r.t protons
	  fContQApid->Fill(5,nsigmaTOF[4],v2V0,xQApid); // v2 TOF (V0C) w.r.t. protons
	}

	// QA fill
	if(!(fPID->GetCurrentMask(0)) || !aodTrack->GetDetPid() || dedx < 10. || aodTrack->Pt() < 0 || aodTrack->Pt() > 7){}
	else{
	  if(TMath::Abs(nsigmaTPC[2])<5 && (!(fPID->GetCurrentMask(1)) || (TMath::Abs(nsigmaTOF[2])<5))){ //pi
	    xQA[2] = prob[2];
	    xQA3[2] = xQA[2];
	    if(fV2) QA[iV0]->Fill(0,nsigmaTPC[2],nsigmaTOF[2],xQA);
 	    if(fV3) QAv3[iV0]->Fill(0,nsigmaTPC[2],nsigmaTOF[2],xQA3);
	  }
	  if(TMath::Abs(nsigmaTPC[3])<5 && (!(fPID->GetCurrentMask(1)) || (TMath::Abs(nsigmaTOF[3])<5))){ //K
	    xQA[2] = prob[3];
	    xQA3[2] = xQA[2];
	    if(fV2) QA[iV0]->Fill(1,nsigmaTPC[3],nsigmaTOF[3],xQA);
// 	    if(fV3) QAv3[iV0]->Fill(1,nsigmaTPC[3],nsigmaTOF[3],xQA3);	  
	  }
	  if(TMath::Abs(nsigmaTPC[4])<5 && (!(fPID->GetCurrentMask(1)) || (TMath::Abs(nsigmaTOF[4])<5))){//p
	    xQA[2] = prob[4];
	    xQA3[2] = xQA[2];
	    if(fV2) QA[iV0]->Fill(2,nsigmaTPC[4],nsigmaTOF[4],xQA);
// 	    if(fV3) QAv3[iV0]->Fill(2,nsigmaTPC[4],nsigmaTOF[4],xQA3);	  
	  }
	  if(TMath::Abs(nsigmaTPC[0])<5 && (!(fPID->GetCurrentMask(1)) || (TMath::Abs(nsigmaTOF[0])<5))){//e
	    xQA[2] = prob[0];
	    xQA3[2] = xQA[2];
// 	    if(fV2) QA[iV0]->Fill(3,nsigmaTPC[0],nsigmaTOF[0],xQA);
// 	    if(fV3) QAv3[iV0]->Fill(3,nsigmaTPC[0],nsigmaTOF[0],xQA3);	  
	  }
	  if(TMath::Abs(nsigmaTPC[5])<5 && (!(fPID->GetCurrentMask(1)) || (TMath::Abs(nsigmaTOF[5])<5))){//d
	    xQA[2] = prob[5];
	    xQA3[2] = xQA[2];
	    // 	  if(fV2) QA[iV0]->Fill(4,nsigmaTPC[5],nsigmaTOF[5],xQA);
	    // 	  if(fV3) QAv3[iV0]->Fill(4,nsigmaTPC[5],nsigmaTOF[5],xQA3);	  
	  }
	  if(TMath::Abs(nsigmaTPC[6])<5 && (!(fPID->GetCurrentMask(1)) || (TMath::Abs(nsigmaTOF[6])<5))){//t
	    xQA[2] = prob[6];
	    xQA3[2] = xQA[2];
	    // 	  if(fV2) QA[iV0]->Fill(5,nsigmaTPC[6],nsigmaTOF[6],xQA);
	    // 	  if(fV3) QAv3[iV0]->Fill(5,nsigmaTPC[6],nsigmaTOF[6],xQA3);	  
	  }
	  if(TMath::Abs(nsigmaTPC[7])<5 && (!(fPID->GetCurrentMask(1)) || (TMath::Abs(nsigmaTOF[7])<5))){//He3
	    xQA[2] = prob[7];
	    xQA3[2] = xQA[2];
	    // 	  if(fV2) QA[iV0]->Fill(6,nsigmaTPC[7],nsigmaTOF[7],xQA);
	    // 	  if(fV3) QAv3[iV0]->Fill(6,nsigmaTPC[7],nsigmaTOF[7],xQA3);	  
	  }
	}

	//pid selection
	if(!(fPID->GetCurrentMask(0)) || !aodTrack->GetDetPid()){} // TPC PID and PID object strictly required (very important!!!!)
	else if(prob[2] > 0.6){ // pi
	  x[2] = prob[2];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[2]) < 5){ // TPC 5 sigma extra cut to accept the track
	    if(fV2) contV0[iV0]->Fill(1,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(1,aodTrack->Pt(),v3V0,x3);
	    if(x[2] > 0.9 && x[5] < 3) fHdcaPt[iC][0]->Fill(aodTrack->Pt(),b[0]);
	    else if(x[2] > 0.9 && fIsMC) fHdcaPtSec[iC][0]->Fill(aodTrack->Pt(),b[0]);
	  }
	}
	else if(prob[3] > 0.6){ // K
	  x[2] = prob[3];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[3]) < 5){
	    if(fV2) contV0[iV0]->Fill(2,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(2,aodTrack->Pt(),v3V0,x3);
	    if(x[2] > 0.9 && x[5] < 3) fHdcaPt[iC][1]->Fill(aodTrack->Pt(),b[0]);
	    else if(x[2] > 0.9 && fIsMC) fHdcaPtSec[iC][1]->Fill(aodTrack->Pt(),b[0]);
	  }
	}
	else if(prob[4] > 0.6){ // p
	  x[2] = prob[4];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[4]) < 5){
	    if(fV2) contV0[iV0]->Fill(3,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(3,aodTrack->Pt(),v3V0,x3);
	    if(x[2] > 0.9 && x[5] < 3 && x[1] < 0) fHdcaPt[iC][2]->Fill(aodTrack->Pt(),b[0]);
	    else if(x[2] > 0.9 && fIsMC && x[1] < 0) fHdcaPtSec[iC][2]->Fill(aodTrack->Pt(),b[0]);
	  }
	}
	else if(prob[0] > 0.6){ // e
	  x[2] = prob[0];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[0]) < 5){
	    if(fV2) contV0[iV0]->Fill(4,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(4,aodTrack->Pt(),v3V0,x3);
	    if(x[2] > 0.9 && x[5] < 3) fHdcaPt[iC][3]->Fill(aodTrack->Pt(),b[0]);
	    else if(x[2] > 0.9 && fIsMC) fHdcaPtSec[iC][3]->Fill(aodTrack->Pt(),b[0]);
	  }
	}
	else if(prob[1] > 0.6){ // mu
	  x[2] = prob[1];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[1]) < 5){
	    if(fV2) contV0[iV0]->Fill(8,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(8,aodTrack->Pt(),v3V0,x3);
	  }
	}
	else if(prob[5] > 0.6){ // d
	  x[2] = prob[5];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[5]) < 5){
	    if(fV2) contV0[iV0]->Fill(5,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(5,aodTrack->Pt(),v3V0,x3);
	    if(x[2] > 0.9 && x[5] < 3 && x[1] < 0) fHdcaPt[iC][4]->Fill(aodTrack->Pt(),b[0]);
	    else if(x[2] > 0.9 && fIsMC && x[1] < 0) fHdcaPtSec[iC][4]->Fill(aodTrack->Pt(),b[0]);
	  }
	}
	else if(prob[6] > 0.6){ // t
	  x[2] = prob[6];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[6]) < 5){
	    if(fV2) contV0[iV0]->Fill(6,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(6,aodTrack->Pt(),v3V0,x3);
	    if(x[2] > 0.9 && x[5] < 3 && x[1] < 0) fHdcaPt[iC][5]->Fill(aodTrack->Pt(),b[0]);
	    else if(x[2] > 0.9 && fIsMC && x[1] < 0) fHdcaPtSec[iC][5]->Fill(aodTrack->Pt(),b[0]);
	  }
	}
	else if(prob[7] > 0.6){ // He3
	  x[2] = prob[7];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[7]) < 5){
	    if(fV2) contV0[iV0]->Fill(7,aodTrack->Pt()*2,v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(7,aodTrack->Pt()*2,v3V0,x3);
	    if(x[2] > 0.9 && x[5] < 3 && x[1] < 0) fHdcaPt[iC][6]->Fill(aodTrack->Pt(),b[0]);
	    else if(x[2] > 0.9 && fIsMC && x[1] < 0) fHdcaPtSec[iC][6]->Fill(aodTrack->Pt(),b[0]);
	  }
	}
	
	if(x[4] > 0.5){ // if TOF was present redo TPC stand alone PID to check the PID in the same acceptance (PID mask = 2)
	  fPID->ResetDetOR(1); // exclude TOF from PID
	  tofMismProb = 0;
	  
	  fPID->ComputeProb(aodTrack,fOutputAOD);
	  dedx = fPID->GetDeDx();//aodTrack->GetTPCsignal();
	  probRead = fPID->GetProb();
	  
	  fPID->SetDetOR(1); // include TOF for PID
	}
	Float_t probTPC[8] = {probRead[0],probRead[1],probRead[2],probRead[3],probRead[4],probRead[5],probRead[6],probRead[7]}; // TPC stand alone prbabilities

	//pid selection TPC S.A. with TOF matching
	x[4]*=2; // set the mask to 2 id TOF is present
	if(x[4]<1 || !(fPID->GetCurrentMask(0)) || !aodTrack->GetDetPid()){} // TPC PID S.A. PID in TOF acceptance
	else if(probTPC[2] > 0.6){ // pi
	  x[2] = probTPC[2];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[2]) < 5){
	    if(fV2) contV0[iV0]->Fill(1,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(1,aodTrack->Pt(),v3V0,x3);
	  }
	}
	else if(probTPC[3] > 0.6){ // K
	  x[2] = probTPC[3];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[3]) < 5){
	    if(fV2) contV0[iV0]->Fill(2,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(2,aodTrack->Pt(),v3V0,x3);
	  }
	}
	else if(probTPC[4] > 0.6){ // p
	  x[2] = probTPC[4];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[4]) < 5){
	    if(fV2) contV0[iV0]->Fill(3,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(3,aodTrack->Pt(),v3V0,x3);
	  }
	}
	else if(probTPC[0] > 0.6){ // e
	  x[2] = probTPC[0];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[0]) < 5){
	    if(fV2) contV0[iV0]->Fill(4,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(4,aodTrack->Pt(),v3V0,x3);
	  }
	}
	else if(probTPC[1] > 0.6){ // mu
	  x[2] = probTPC[1];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[1]) < 5){
	    if(fV2) contV0[iV0]->Fill(8,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(8,aodTrack->Pt(),v3V0,x3);
	  }
	}
	else if(probTPC[5] > 0.6){ // d
	  x[2] = probTPC[5];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[5]) < 5){
	    if(fV2) contV0[iV0]->Fill(5,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(5,aodTrack->Pt(),v3V0,x3);
	  }
	}
	else if(probTPC[6] > 0.6){ // t
	  x[2] = probTPC[6];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[6]) < 5){
	    if(fV2) contV0[iV0]->Fill(6,aodTrack->Pt(),v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(6,aodTrack->Pt(),v3V0,x3);
	  }
	}
	else if(probTPC[7] > 0.6){ // He3
	  x[2] = probTPC[7];
	  x3[2] = x[2];
	  if(TMath::Abs(nsigmaTPC[7]) < 5){
	    if(fV2) contV0[iV0]->Fill(7,aodTrack->Pt()*2,v2V0,x);
	    if(fV3) contV0v3[iV0]->Fill(7,aodTrack->Pt()*2,v3V0,x3);
	  }
	}
      } // end side loop
    } // end track loop

    // V0 loop
    Int_t nV0s = fOutputAOD->GetNumberOfV0s();
    AliAODv0 *myV0;
    //    Double_t dQT, dPT, dALPHA,
    Double_t dMASS=0.0;

    for (Int_t i=0; i!=nV0s; ++i) {
      myV0 = (AliAODv0*) fOutputAOD->GetV0(i);
      if(!myV0) continue;
      if(myV0->Pt()<0.1 || TMath::Abs(myV0->Eta()) > fEtaCut) continue; // skipping low momentum
      Int_t pass = PassesAODCuts(myV0,fOutputAOD,0);
      if(pass) {
	dMASS = myV0->MassK0Short();
	pass = 3;
      }
      else {
	pass = PassesAODCuts(myV0,fOutputAOD,1);
	if(pass) dMASS = myV0->MassLambda();
	if(pass==2) dMASS = myV0->MassAntiLambda();
      }
      if(pass){// 1 lambda, 2 antilambda, 3=K0s
	// dPT=myV0->Pt();
	// dQT=myV0->PtArmV0();
	// dALPHA=myV0->AlphaV0();

	Int_t iPos, iNeg;
	AliAODTrack *iT=(AliAODTrack*) myV0->GetDaughter(0);
	if(iT->Charge()>0) {
	  iPos = 0; iNeg = 1;
	} else {
	  iPos = 1; iNeg = 0;
	}
	iT=(AliAODTrack*) myV0->GetDaughter(iPos); // positive
	AliAODTrack *jT=(AliAODTrack*) myV0->GetDaughter(iNeg); // negative

	// re-map the container in an array to do the analysis for V0A and V0C within a loop
	Float_t evPlAngV0[2] = {evPlAngV0ACor2,evPlAngV0CCor2};
	AliFlowVZEROResults *contV0[2] = {fContAllChargesV0A,fContAllChargesV0C};
	
	Float_t evPlAngV0v3[2] = {evPlAngV0ACor3,evPlAngV0CCor3};
	AliFlowVZEROResults *contV0v3[2] = {fContAllChargesV0Av3,fContAllChargesV0Cv3};

	for(Int_t iV0=0;iV0<2;iV0++){ // loop on A and C side
	 
	  Float_t v2V0 = TMath::Cos(2*(myV0->Phi() - evPlAngV0[iV0]));
	  Float_t v3V0 = TMath::Cos(3*(myV0->Phi() - evPlAngV0v3[iV0]));
	  
	  Float_t x[6] = {Float_t(iC),1,1,evPlAngV0[iV0],1,0}; // to fill analysis v2 container
	  Float_t x3[6] = {Float_t(iC),1,1,evPlAngV0v3[iV0],1,0}; // to fill analysis v3 container
	  
	  Float_t decaylength = myV0->DecayLengthXY(fOutputAOD->GetPrimaryVertex());
	  //	  printf("decay length = %f\n",decaylength);

	  if(pass==2){ // anti-lambda charge = -1
	    x[1] = -1;
	    x3[1] = -1;
	  }

	  if(decaylength < fMinDistV0) pass = 0;	  
	  if(decaylength > fMaxDistV0) pass = 0;	  

	  Float_t nsigma = 0;
	  if(pass < 3)
	    nsigma = TMath::Abs(dMASS-1.116)/0.0016;
	  else if(pass == 3)
	    nsigma = TMath::Abs(dMASS-0.497)/0.005;

	  if(nsigma < 1)
	    x[2] = 0.95;
	  else if(nsigma < 2)
	    x[2] = 0.85;
	  else if(nsigma < 3)
	    x[2] = 0.75;
	  else if(nsigma < 4)
	    x[2] = 0.65;
	  else
	    x[2] = 0.5;
	  	    
	  x3[2] = x[2];

	  // Fill Container for lambda and Ks
	  if(fV2 && pass == 3 && x[2] > 0.6) contV0[iV0]->Fill(9,myV0->Pt(),v2V0,x);
	  if(fV3 && pass == 3 && x[2] > 0.6) contV0v3[iV0]->Fill(9,myV0->Pt(),v3V0,x3);
	  if(fV2 && pass < 3 && x[2] > 0.6) contV0[iV0]->Fill(10,myV0->Pt(),v2V0,x);
	  if(fV3 && pass < 3 && x[2] > 0.6) contV0v3[iV0]->Fill(10,myV0->Pt(),v3V0,x3);

	  if(pass < 3){ // lambda
	    AliAODTrack* aodTrack = iT;
	    if(pass==2) aodTrack=jT;

	    v2V0 = TMath::Cos(2*(aodTrack->Phi() - evPlAngV0[iV0]));
	    v3V0 = TMath::Cos(3*(aodTrack->Phi() - evPlAngV0v3[iV0]));

	    fPID->ComputeProb(aodTrack,fOutputAOD); // compute Bayesian probabilities
	    Float_t *probRead = fPID->GetProb();
	    Float_t prob[8] = {probRead[0],probRead[1],probRead[2],probRead[3],probRead[4],probRead[5],probRead[6],probRead[7]};
	    Float_t tofMismProb = fPID->GetTOFMismProb(); // TOF mismatch probability requested to be lower than 50% for TOF analysis 
	    
	    Float_t xdec[6] = {Float_t(iC),Float_t(aodTrack->Charge()),prob[4],evPlAngV0[iV0],Float_t(fPID->GetCurrentMask(1)&&tofMismProb < 0.5),0}; // to fill analysis v2 container
	    Float_t xdec3[6] = {Float_t(iC),Float_t(aodTrack->Charge()),prob[4],evPlAngV0v3[iV0],Float_t(fPID->GetCurrentMask(1)&&tofMismProb < 0.5),0}; // to fill analysis v3 container

	    // Fill Container for (anti)proton from lambda
	    if(nsigma < 2 && xdec[2] > 0.6){
	      if(fV2) contV0[iV0]->Fill(11,aodTrack->Pt(),v2V0,xdec);
	      if(fV3) contV0v3[iV0]->Fill(11,aodTrack->Pt(),v3V0,xdec3);
	    }
	  }
	}
	
      }
    } // end loop on V0


    // Fill EP distribution histograms
    if(fV2) fPhiRPv0A->Fill(iC,evPlAngV0ACor2);
    if(fV2) fPhiRPv0C->Fill(iC,evPlAngV0CCor2);
    
    if(fV3) fPhiRPv0Av3->Fill(iC,evPlAngV0ACor3);
    if(fV3) fPhiRPv0Cv3->Fill(iC,evPlAngV0CCor3);

    // TPC EP needed for resolution studies (TPC subevent)
    Double_t Qx2 = 0, Qy2 = 0;
    Double_t Qx3 = 0, Qy3 = 0;

    for(Int_t iT = 0; iT < nAODTracks; iT++) {
      
      AliAODTrack* aodTrack = aodEvent->GetTrack(iT);
      
      if (!aodTrack){
	aodTrack->Delete();
	continue;
      }
      
      Bool_t trkFlag = aodTrack->TestFilterBit(1);

      if ((TMath::Abs(aodTrack->Eta()) > 0.8) || (aodTrack->Pt() < 0.2) || (aodTrack->GetTPCNcls() < fNcluster)  || !trkFlag) 
	continue;
	
      Double_t b[2] = {-99., -99.};
      Double_t bCov[3] = {-99., -99., -99.};
      if (!aodTrack->PropagateToDCA(fOutputAOD->GetPrimaryVertex(), fOutputAOD->GetMagneticField(), 100., b, bCov))
	continue;
	    
      if ((TMath::Abs(b[0]) > 3.0) || (TMath::Abs(b[1]) > 2.4))
	continue;
      
      Qx2 += TMath::Cos(2*aodTrack->Phi()); 
      Qy2 += TMath::Sin(2*aodTrack->Phi());
      Qx3 += TMath::Cos(3*aodTrack->Phi()); 
      Qy3 += TMath::Sin(3*aodTrack->Phi());
      
    }
    
    evPlAng2 = TMath::ATan2(Qy2, Qx2)/2.;
    evPlAng3 = TMath::ATan2(Qy3, Qx3)/3.;

    fgPsi2tpc = evPlAng2;
    fgPsi3tpc = evPlAng3;

    // Fill histograms needed for resolution evaluation
    if(fV2) fHResTPCv0A2->Fill(Double_t(iC), TMath::Cos(2*(evPlAng2 - evPlAngV0ACor2)));
    if(fV2) fHResTPCv0C2->Fill(Double_t(iC), TMath::Cos(2*(evPlAng2 - evPlAngV0CCor2)));
    if(fV2) fHResv0Cv0A2->Fill(Double_t(iC), TMath::Cos(2*(evPlAngV0ACor2 - evPlAngV0CCor2)));
    
    if(fV3) fHResTPCv0A3->Fill(Double_t(iC), TMath::Cos(3*(evPlAng3 - evPlAngV0ACor3)));
    if(fV3) fHResTPCv0C3->Fill(Double_t(iC), TMath::Cos(3*(evPlAng3 - evPlAngV0CCor3)));
    if(fV3) fHResv0Cv0A3->Fill(Double_t(iC), TMath::Cos(3*(evPlAngV0ACor3 - evPlAngV0CCor3)));
  }
  
}

//_____________________________________________________________________________
Float_t AliAnalysisTaskVnV0::GetVertex(AliAODEvent* aod) const
{

  Float_t zvtx = -999;

  const AliAODVertex* vtxAOD = aod->GetPrimaryVertex();
  if (!vtxAOD)
    return zvtx;
  if(vtxAOD->GetNContributors()>0)
    zvtx = vtxAOD->GetZ();
  
  return zvtx;
}
//_____________________________________________________________________________
void AliAnalysisTaskVnV0::Terminate(Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
//_____________________________________________________________________________
void AliAnalysisTaskVnV0::OpenInfoCalbration(Int_t run){
    TString oadbfilename = "$ALICE_ROOT/OADB/PWGCF/VZERO/VZEROcalibEP.root";
    TFile *foadb = TFile::Open(oadbfilename.Data());

    if(!foadb){
	printf("OADB file %s cannot be opened\n",oadbfilename.Data());
	return;
    }

    AliOADBContainer *cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorr");
    if(!cont){
	printf("OADB object hMultV0BefCorr is not available in the file\n");
	return;	
    }

    if(!(cont->GetObject(run))){
	printf("OADB object hMultV0BefCorr is not available for run %i (used run 137366)\n",run);
	run = 137366;
    }
    fMultV0 = ((TH2F *) cont->GetObject(run))->ProfileX();

    TF1 *fpol0 = new TF1("fpol0","pol0"); 
    fMultV0->Fit(fpol0,"","",0,31);
    fV0Cpol = fpol0->GetParameter(0);
    fMultV0->Fit(fpol0,"","",32,64);
    fV0Apol = fpol0->GetParameter(0);

    for(Int_t iside=0;iside<2;iside++){
	for(Int_t icoord=0;icoord<2;icoord++){
	    for(Int_t i=0;i  < 9;i++){
		char namecont[100];
  		if(iside==0 && icoord==0)
		  snprintf(namecont,100,"hQxc2_%i",i);
		else if(iside==1 && icoord==0)
		  snprintf(namecont,100,"hQxa2_%i",i);
		else if(iside==0 && icoord==1)
		  snprintf(namecont,100,"hQyc2_%i",i);
		else if(iside==1 && icoord==1)
		  snprintf(namecont,100,"hQya2_%i",i);

		cont = (AliOADBContainer*) foadb->Get(namecont);
		if(!cont){
		    printf("OADB object %s is not available in the file\n",namecont);
		    return;	
		}
		
		if(!(cont->GetObject(run))){
		    printf("OADB object %s is not available for run %i (used run 137366)\n",namecont,run);
		    run = 137366;
		}
		fMeanQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
		fWidthQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();

		//for v3
		if(iside==0 && icoord==0)
		  snprintf(namecont,100,"hQxc3_%i",i);
		else if(iside==1 && icoord==0)
		  snprintf(namecont,100,"hQxa3_%i",i);
		else if(iside==0 && icoord==1)
		  snprintf(namecont,100,"hQyc3_%i",i);
		else if(iside==1 && icoord==1)
		  snprintf(namecont,100,"hQya3_%i",i);

		cont = (AliOADBContainer*) foadb->Get(namecont);
		if(!cont){
		    printf("OADB object %s is not available in the file\n",namecont);
		    return;	
		}
		
		if(!(cont->GetObject(run))){
		    printf("OADB object %s is not available for run %i (used run 137366)\n",namecont,run);
		    run = 137366;
		}
		fMeanQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
		fWidthQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();

     	    }
	}
    }
}
//=======================================================================
Int_t AliAnalysisTaskVnV0::PassesAODCuts(AliAODv0 *myV0, AliAODEvent *tAOD,Int_t specie)
{
  Int_t set = 2;
  Float_t fV0Cuts[9];
  // defines cuts to be used
  // fV0Cuts[9] dl dca ctp d0 d0d0 qt minEta maxEta PID
  switch(set) {
  case(0): // No cuts
    fV0Cuts[0] = -1e+6; fV0Cuts[1] = +1e+6; fV0Cuts[2] = -1e+6;
    fV0Cuts[3] = -1e+6; fV0Cuts[4] = +1e+6; fV0Cuts[5] = -1e+6;
    fV0Cuts[6] = -1e+6; fV0Cuts[7] = +1e+6; fV0Cuts[8] = 0;
    break;
  case(1): // Tight cuts
    fV0Cuts[0] = +0.5; fV0Cuts[1] = +0.5; fV0Cuts[2] = +0.998;
    fV0Cuts[3] = +0.1; fV0Cuts[4] = +0.0; fV0Cuts[5] = +0.105;
    fV0Cuts[6] = -0.8; fV0Cuts[7] = +0.8; fV0Cuts[8] = 0;
    break;
  case(2): // Tight cuts + PID
    fV0Cuts[0] = +0.5; fV0Cuts[1] = +0.5; fV0Cuts[2] = +0.998;
    fV0Cuts[3] = +0.1; fV0Cuts[4] = +0.0; fV0Cuts[5] = +0.105;
    fV0Cuts[6] = -0.8; fV0Cuts[7] = +0.8; fV0Cuts[8] = 1;
    break;
  case(3): // No cuts + PID
    fV0Cuts[0] = -1e+6; fV0Cuts[1] = +1e+6; fV0Cuts[2] = -1e+6;
    fV0Cuts[3] = -1e+6; fV0Cuts[4] = +1e+6; fV0Cuts[5] = -1e+6;
    fV0Cuts[6] = -1e+6; fV0Cuts[7] = +1e+6; fV0Cuts[8] = 1;
    break;
  }

  // daughter cuts
  if(! fCutsDaughter){
    fCutsDaughter = new AliESDtrackCuts(Form("daughter_cuts_%s","ESD") );
    fCutsDaughter->SetPtRange(0.2,10.0);
    fCutsDaughter->SetEtaRange(-0.8, 0.8 );
    fCutsDaughter->SetMinNClustersTPC(80);
    fCutsDaughter->SetMaxChi2PerClusterTPC(4.0);
    fCutsDaughter->SetRequireTPCRefit(kTRUE);
    fCutsDaughter->SetAcceptKinkDaughters(kFALSE);
  }

  if (myV0->GetOnFlyStatus() ) return 0;
  //the following is needed in order to evualuate track-quality
  AliAODTrack *iT, *jT;
  AliAODVertex *vV0s = myV0->GetSecondaryVtx();
  Double_t pos[3],cov[6];
  vV0s->GetXYZ(pos);
  vV0s->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);
  // TESTING CHARGE
  int iPos, iNeg;
  iT=(AliAODTrack*) myV0->GetDaughter(0);
  if(iT->Charge()>0) {
    iPos = 0; iNeg = 1;
  } else {
    iPos = 1; iNeg = 0;
  }
  // END OF TEST

  iT=(AliAODTrack*) myV0->GetDaughter(iPos); // positive
  AliESDtrack ieT( iT );
  ieT.SetTPCClusterMap( iT->GetTPCClusterMap() );
  ieT.SetTPCSharedMap( iT->GetTPCSharedMap() );
  ieT.SetTPCPointsF( iT->GetTPCNclsF() );
  ieT.RelateToVertex(&vESD, tAOD->GetMagneticField(), 100);
  if (!fCutsDaughter->IsSelected( &ieT ) ) return 0;

  jT=(AliAODTrack*) myV0->GetDaughter(iNeg); // negative
  AliESDtrack jeT( jT );
  jeT.SetTPCClusterMap( jT->GetTPCClusterMap() );
  jeT.SetTPCSharedMap( jT->GetTPCSharedMap() );
  jeT.SetTPCPointsF( jT->GetTPCNclsF() );
  jeT.RelateToVertex(&vESD, tAOD->GetMagneticField(), 100);
  if (!fCutsDaughter->IsSelected( &jeT ) ) return 0;

  Double_t pvertex[3];
  pvertex[0]=tAOD->GetPrimaryVertex()->GetX();
  pvertex[1]=tAOD->GetPrimaryVertex()->GetY();
  pvertex[2]=tAOD->GetPrimaryVertex()->GetZ();
  Double_t dDL=myV0->DecayLengthV0( pvertex );
  Double_t dDCA=myV0->DcaV0Daughters();
  Double_t dCTP=myV0->CosPointingAngle( pvertex );
  Double_t dD0P=ieT.GetD(pvertex[0],pvertex[1],tAOD->GetMagneticField());
  Double_t dD0M=jeT.GetD(pvertex[0],pvertex[1],tAOD->GetMagneticField());
  Double_t dD0D0=dD0P*dD0M;
  Double_t dQT=myV0->PtArmV0();
  Double_t dALPHA=myV0->AlphaV0(); // AlphaV0 -> AODRecoDecat::Alpha -> return 1.-2./(1.+QlProng(0)/QlProng(1));
  if(myV0->ChargeProng(iPos)<0) dALPHA = -dALPHA; // protects for a change in convention
//   Double_t dPT=myV0->Pt();
  Double_t dETA=myV0->Eta();
  Int_t passes = 1;
  if(dDL  <fV0Cuts[0]) passes = 0;
  if(dDCA >fV0Cuts[1]) passes = 0;
  if(dCTP <fV0Cuts[2]) passes = 0;
  if(TMath::Abs(dD0P) <fV0Cuts[3]) passes = 0;
  if(TMath::Abs(dD0M) <fV0Cuts[3]) passes = 0;
  if(dD0D0>fV0Cuts[4]) passes = 0;
  if(dETA <fV0Cuts[6]) passes = 0;
  if(dETA >fV0Cuts[7]) passes = 0;
  if(specie==0) if(dQT<fV0Cuts[5]) passes = 0;
  if(specie==1&&passes==1&&dALPHA<0) passes = 2; // antilambda
  if(passes&&fV0Cuts[8]) {

    Double_t dedxExp[8];
    fPID->ComputeProb(iT,tAOD); // compute Bayesian probabilities
    Float_t nsigmaTPC[8];
    if(iT->GetDetPid()){ // check the PID object is available
      for(Int_t iS=0;iS < 8;iS++){
	dedxExp[iS] = fPID->GetExpDeDx(iT,iS);
	nsigmaTPC[iS] = (fPID->GetDeDx() - dedxExp[iS])/(dedxExp[iS]*0.07);
      }
    }
    else{
      for(Int_t iS=0;iS < 8;iS++)
	nsigmaTPC[iS] = 10;
    }

    fPID->ComputeProb(jT,tAOD); // compute Bayesian probabilities
    Float_t nsigmaTPC2[8];
    if(jT->GetDetPid()){ // check the PID object is available
      for(Int_t iS=0;iS < 8;iS++){
	dedxExp[iS] = fPID->GetExpDeDx(jT,iS);
	nsigmaTPC2[iS] = (fPID->GetDeDx() - dedxExp[iS])/(dedxExp[iS]*0.07);
      }
    }
    else{
      for(Int_t iS=0;iS < 8;iS++)
	nsigmaTPC2[iS] = 10;
    }

    if(jT->GetTPCNcls() < fNcluster) passes = 0;
    else if(iT->GetTPCNcls() < fNcluster) passes = 0;

    switch(specie) {
    case 0: // K0 PID
      if( (jT->GetTPCmomentum()<15) &&
	  (TMath::Abs(nsigmaTPC2[2])>3.) )
	passes = 0;
      if( (iT->GetTPCmomentum()<15) &&
	  (TMath::Abs(nsigmaTPC[2])>3.) )
	passes = 0;
      break;
    case 1: // Lambda PID  i==pos j ==neg
      if(passes==1) {
	if( (iT->GetTPCmomentum()<15) &&
	    (TMath::Abs(nsigmaTPC[4])>3.) )
	  passes = 0;
	if( (jT->GetTPCmomentum()<15) &&
	    (TMath::Abs(nsigmaTPC2[2])>3.) )
	  passes = 0;
      }
      if(passes==2) {
	if( (iT->GetTPCmomentum()<15) &&
	    (TMath::Abs(nsigmaTPC[2])>3.) )
	  passes = 0;
	if( (jT->GetTPCmomentum()<15) &&
	    (TMath::Abs(nsigmaTPC2[4])>3.) )
	  passes = 0;
      }
      break;
    }
  }
  return passes;
}
