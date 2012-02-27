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

// STL includes
//#include <iostream>
//using namespace std;

ClassImp(AliAnalysisTaskVnV0)

/*
fMeanQ[nCentrBin][2][2];           // and recentering
  Float_t fWidthQ[nCentrBin][2][2];          // ...
  Float_t fMeanQv3[nCentrBin][2][2];         // also for v3
  Float_t fWidthQv3[nCentrBin][2][2];        // ...
*/

//_____________________________________________________________________________
AliAnalysisTaskVnV0::AliAnalysisTaskVnV0():
  AliAnalysisTaskSE(),
  fAOD(0),
  fVtxCut(10),    // cut on |vertex| < fVtxCut
  fEtaCut(0.8),   // cut on |eta| < fEtaCut
  fMinPt(0.15),   // cut on pt > fMinPt
  fRun(-1),
  fList(new TList()),
  fList2(new TList()),
  fMultV0(0),
  fV0Cpol(100),
  fV0Apol(100),
  fContAllChargesV0A(0),
  fContAllChargesV0C(0),
  fContAllChargesV0Av3(0),
  fContAllChargesV0Cv3(0),
  fHResTPCv0A2(0),
  fHResTPCv0C2(0),
  fHResv0Cv0A2(0),
  fHResTPCv0A3(0),
  fHResTPCv0C3(0),
  fHResv0Cv0A3(0),
  fPhiRPv0A(0),
  fPhiRPv0C(0),
  fPhiRPv0Av3(0),
  fPhiRPv0Cv3(0),
  fPhiTracks(0),
  fQA(0),
  fQA2(0),
  fQAv3(0),
  fQA2v3(0),
  fPID(new AliFlowBayesianPID()),
  fTree(0),
  fCentrality(0),
  evPlAngV0ACor2(0),
  evPlAngV0CCor2(0),
  evPlAng2(0),
  evPlAngV0ACor3(0),
  evPlAngV0CCor3(0),
  evPlAng3(0),
  fV2(kTRUE),
  fV3(kTRUE)
{
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());

  // Default constructor (should not be used)
  fList->SetName("resultsV2");
  fList2->SetName("resultsV3");

  fPID->SetNewTrackParam(); // Better tuning for TOF PID tracking effect in LHC10h

}

//______________________________________________________________________________
AliAnalysisTaskVnV0::AliAnalysisTaskVnV0(const char *name):
  AliAnalysisTaskSE(name),
  fAOD(0),
  fVtxCut(10),    // cut on |vertex| < fVtxCut
  fEtaCut(0.8),   // cut on |eta| < fEtaCut
  fMinPt(0.15),   // cut on pt > fMinPt
  fRun(-1),
  fList(new TList()),
  fList2(new TList()),
  fMultV0(0),
  fV0Cpol(100),
  fV0Apol(100),
  fContAllChargesV0A(0),
  fContAllChargesV0C(0),
  fContAllChargesV0Av3(0),
  fContAllChargesV0Cv3(0),
  fHResTPCv0A2(0),
  fHResTPCv0C2(0),
  fHResv0Cv0A2(0),
  fHResTPCv0A3(0),
  fHResTPCv0C3(0),
  fHResv0Cv0A3(0),
  fPhiRPv0A(0),
  fPhiRPv0C(0),
  fPhiRPv0Av3(0),
  fPhiRPv0Cv3(0),
  fPhiTracks(0),
  fQA(0),
  fQA2(0),
  fQAv3(0),
  fQA2v3(0),
  fPID(new AliFlowBayesianPID()),
  fTree(0),
  fCentrality(0),
  evPlAngV0ACor2(0),
  evPlAngV0CCor2(0),
  evPlAng2(0),
  evPlAngV0ACor3(0),
  evPlAngV0CCor3(0),
  evPlAng3(0),
  fV2(kTRUE),
  fV3(kTRUE)
{

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());

  // Output slot #1 writes into a TTree
  fList->SetName("resultsV2");
  fList2->SetName("resultsV3");

  fPID->SetNewTrackParam(); // Better tuning for TOF PID tracking effect in LHC10h
}

//_____________________________________________________________________________
AliAnalysisTaskVnV0::~AliAnalysisTaskVnV0()
{

}

//______________________________________________________________________________
void AliAnalysisTaskVnV0::UserCreateOutputObjects()
{

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
  const Int_t nChargeBinsTOF = 3;	  
  Double_t binChargeTOF[nChargeBinsTOF+1] = {-1.5,-0.5,0.5,1.5};
  const Int_t nv_2BinsTOF = 50;
  Double_t binV_2TOF[nv_2BinsTOF+1];
  for(Int_t i=0;i<nv_2BinsTOF+1;i++){
    binV_2TOF[i] = -1 + i*2./nv_2BinsTOF;
  }
  const Int_t nCentrTOF = 9;
  Double_t binCentrTOF[nCentrTOF+1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
  const Int_t nProbTOF = 4;
  Double_t binProbTOF[nProbTOF+1] = {0.6,0.8,0.9,0.95,1.00001};

  const Int_t nPsiTOF = 20;
  Double_t binPsiTOF[nPsiTOF+1];
  for(Int_t i=0;i<nPsiTOF+1;i++){
    binPsiTOF[i] = -TMath::Pi() + i*TMath::Pi()/nPsiTOF*2;
  }
  Double_t binPsiTOFv3[nPsiTOF+1];
  for(Int_t i=0;i<nPsiTOF+1;i++){
    binPsiTOFv3[i] = -TMath::Pi()/3*2 + i*TMath::Pi()/nPsiTOF*2/3*2;
  }
  
  const Int_t nMaskPID = 3;
  Double_t binMaskPID[nMaskPID+1] = {-0.5,0.5,1.5,2.5};

  Int_t binsTOF[7] = {nCentrTOF,nv_2BinsTOF,nChargeBinsTOF,nPtBinsTOF,nProbTOF,nPsiTOF,nMaskPID};

  // v2 container
  fContAllChargesV0A = new AliCFContainer("fContAllChargesV0Av2","centr:v2:charge:pt:prob:Psi:PIDmask",1+7,7,binsTOF);
  fContAllChargesV0C = new AliCFContainer("fContAllChargesV0Cv2","centr:v2:charge:pt:prob",1+7,7,binsTOF);
  fContAllChargesV0A->SetBinLimits(0,binCentrTOF);
  fContAllChargesV0A->SetBinLimits(1,binV_2TOF);
  fContAllChargesV0A->SetBinLimits(2,binChargeTOF);
  fContAllChargesV0A->SetBinLimits(3,binsPtTOF);
  fContAllChargesV0A->SetBinLimits(4,binProbTOF);
  fContAllChargesV0A->SetBinLimits(5,binPsiTOF);
  fContAllChargesV0A->SetBinLimits(6,binMaskPID);
  fContAllChargesV0C->SetBinLimits(0,binCentrTOF);
  fContAllChargesV0C->SetBinLimits(1,binV_2TOF);
  fContAllChargesV0C->SetBinLimits(2,binChargeTOF);
  fContAllChargesV0C->SetBinLimits(3,binsPtTOF);
  fContAllChargesV0C->SetBinLimits(4,binProbTOF);
  fContAllChargesV0C->SetBinLimits(5,binPsiTOF);
  fContAllChargesV0C->SetBinLimits(6,binMaskPID);

  fContAllChargesV0A->SetVarTitle(0,"centrality");
  fContAllChargesV0A->SetVarTitle(1,"v_{2}");
  fContAllChargesV0A->SetVarTitle(2,"charge");
  fContAllChargesV0A->SetVarTitle(3,"p_{t} (GeV/c)");
  fContAllChargesV0A->SetVarTitle(4,"Bayesian Probability");
  fContAllChargesV0A->SetVarTitle(5,"#Psi V0A");
  fContAllChargesV0A->SetVarTitle(6,"TOF PID");
  fContAllChargesV0C->SetVarTitle(0,"centrality");
  fContAllChargesV0C->SetVarTitle(1,"v_{2}");
  fContAllChargesV0C->SetVarTitle(2,"charge");
  fContAllChargesV0C->SetVarTitle(3,"p_{t} (GeV/c)");
  fContAllChargesV0C->SetVarTitle(4,"Bayesian Probability");
  fContAllChargesV0C->SetVarTitle(5,"#Psi V0C");
  fContAllChargesV0C->SetVarTitle(6,"TOF PID");

  fList->Add(fContAllChargesV0A);
  fList->Add(fContAllChargesV0C);

  // v3 container
  fContAllChargesV0Av3 = new AliCFContainer("fContAllChargesV0Av3","centr:v3:charge:pt:prob:Psi:PIDmask",1+7,7,binsTOF);
  fContAllChargesV0Cv3 = new AliCFContainer("fContAllChargesV0Cv3","centr:v3:charge:pt:prob:Psi:PIDmask",1+7,7,binsTOF);
  fContAllChargesV0Av3->SetBinLimits(0,binCentrTOF);
  fContAllChargesV0Av3->SetBinLimits(1,binV_2TOF);
  fContAllChargesV0Av3->SetBinLimits(2,binChargeTOF);
  fContAllChargesV0Av3->SetBinLimits(3,binsPtTOF);
  fContAllChargesV0Av3->SetBinLimits(4,binProbTOF);
  fContAllChargesV0Av3->SetBinLimits(5,binPsiTOFv3);
  fContAllChargesV0Av3->SetBinLimits(6,binMaskPID);
  fContAllChargesV0Cv3->SetBinLimits(0,binCentrTOF);
  fContAllChargesV0Cv3->SetBinLimits(1,binV_2TOF);
  fContAllChargesV0Cv3->SetBinLimits(2,binChargeTOF);
  fContAllChargesV0Cv3->SetBinLimits(3,binsPtTOF);
  fContAllChargesV0Cv3->SetBinLimits(4,binProbTOF);
  fContAllChargesV0Cv3->SetBinLimits(5,binPsiTOFv3);
  fContAllChargesV0Cv3->SetBinLimits(6,binMaskPID);

  fContAllChargesV0Av3->SetVarTitle(0,"centrality");
  fContAllChargesV0Av3->SetVarTitle(1,"v_{2}");
  fContAllChargesV0Av3->SetVarTitle(2,"charge");
  fContAllChargesV0Av3->SetVarTitle(3,"p_{t} (GeV/c)");
  fContAllChargesV0Av3->SetVarTitle(4,"Bayesian Probability");
  fContAllChargesV0Av3->SetVarTitle(5,"#Psi V0A");
  fContAllChargesV0Av3->SetVarTitle(6,"TOF PID");
  fContAllChargesV0Cv3->SetVarTitle(0,"centrality");
  fContAllChargesV0Cv3->SetVarTitle(1,"v_{3}");
  fContAllChargesV0Cv3->SetVarTitle(2,"charge");
  fContAllChargesV0Cv3->SetVarTitle(3,"p_{t} (GeV/c)");
  fContAllChargesV0Cv3->SetVarTitle(4,"Bayesian Probability");
  fContAllChargesV0Cv3->SetVarTitle(5,"#Psi V0C");
  fContAllChargesV0Cv3->SetVarTitle(6,"TOF PID");

  fList2->Add(fContAllChargesV0Av3);
  fList2->Add(fContAllChargesV0Cv3);

  // TProfile for resolutions 3 subevents (V0A, V0C, TPC)
  // v2
  fHResTPCv0A2 = new TProfile("hResTPCv0A2","",9,0,9);
  fHResTPCv0C2 = new TProfile("hResTPCv0C2","",9,0,9);
  fHResv0Cv0A2 = new TProfile("hResv0Cv0A2","",9,0,9);

  fList->Add(fHResTPCv0A2);
  fList->Add(fHResTPCv0C2);
  fList->Add(fHResv0Cv0A2);

  // v3
  fHResTPCv0A3 = new TProfile("hResTPCv0A3","",9,0,9);
  fHResTPCv0C3 = new TProfile("hResTPCv0C3","",9,0,9);
  fHResv0Cv0A3 = new TProfile("hResv0Cv0A3","",9,0,9);

  fList2->Add(fHResTPCv0A3);
  fList2->Add(fHResTPCv0C3);
  fList2->Add(fHResv0Cv0A3);

  // V0A and V0C event plane distributions
  //v2 
  fPhiRPv0A = new TH2F("fPhiRPv0Av2","#phi distribution of EP VZERO-A;centrality;#phi (rad)",9,0,9,nPsiTOF,-TMath::Pi(),TMath::Pi());
  fPhiRPv0C = new TH2F("fPhiRPv0Cv2","#phi distribution of EP VZERO-C;centrality;#phi (rad)",9,0,9,nPsiTOF,-TMath::Pi(),TMath::Pi());

  //v3
  fPhiRPv0Av3 = new TH2F("fPhiRPv0Av3","#phi distribution of EP VZERO-A;centrality;#phi (rad)",9,0,9,nPsiTOF,-TMath::Pi()/3*2,TMath::Pi()/3*2);
  fPhiRPv0Cv3 = new TH2F("fPhiRPv0Cv3","#phi distribution of EP VZERO-C;centrality;#phi (rad)",9,0,9,nPsiTOF,-TMath::Pi()/3*2,TMath::Pi()/3*2);

  // Track phi distribution container (only once) put in v2 output (only if needed)
  const Int_t nPhiBin = 20;
  Double_t binPhi[nPhiBin+1];
  for(Int_t i=0;i<nPhiBin+1;i++){
        binPhi[i] = i*TMath::Pi()*0.1;
  }  
  Int_t binsPhi[5] = {nCentrTOF,nChargeBinsTOF,nPtBinsTOF,nProbTOF,nPhiBin};
  fPhiTracks = new AliCFContainer("fPhiTracks","centr:charge:pt:prob:phi",1+7,5,binsPhi);
  fPhiTracks->SetBinLimits(0,binCentrTOF);
  fPhiTracks->SetBinLimits(1,binChargeTOF);
  fPhiTracks->SetBinLimits(2,binsPtTOF);
  fPhiTracks->SetBinLimits(3,binProbTOF);
  fPhiTracks->SetBinLimits(4,binPhi);
  fPhiTracks->SetVarTitle(0,"centrality");
  fPhiTracks->SetVarTitle(1,"charge");
  fPhiTracks->SetVarTitle(2,"p_{t} (GeV/c)");
  fPhiTracks->SetVarTitle(3,"Bayesian Probability");
  fPhiTracks->SetVarTitle(4,"#phi");

  // QA container
  // v2
  const Int_t nDETsignal = 50;
  Double_t binDETsignal[nDETsignal+1];
  for(Int_t i=0;i<nDETsignal+1;i++){
    binDETsignal[i] = -5 + i*10. / nDETsignal;
  }

  const Int_t nDeltaPhi = 5;
  Double_t binDeltaPhi[nDeltaPhi+1];
  for(Int_t i=0;i<nDeltaPhi+1;i++){
    binDeltaPhi[i] = -TMath::Pi() + i*2*TMath::Pi() / nDeltaPhi;
  }

  Int_t binsQA[7] = {nCentrTOF,nPtBinsTOF,nProbTOF,nDETsignal,nDETsignal,nDeltaPhi,nMaskPID};

  fQA  = new AliCFContainer("fQAv2","centr:pt:prob:TPCsig:TOFsig:DeltaPhi:maskPID",7,7,binsQA);
  fQA->SetBinLimits(0,binCentrTOF);
  fQA->SetBinLimits(1,binsPtTOF);
  fQA->SetBinLimits(2,binProbTOF);
  fQA->SetBinLimits(3,binDETsignal);
  fQA->SetBinLimits(4,binDETsignal);
  fQA->SetBinLimits(5,binDeltaPhi);
  fQA->SetBinLimits(6,binMaskPID);
  fQA->SetVarTitle(0,"centrality");
  fQA->SetVarTitle(1,"p_{t} (GeV/c)");
  fQA->SetVarTitle(2,"Bayesian Probability");
  fQA->SetVarTitle(3,"N_{#sigma}^{TPC}");
  fQA->SetVarTitle(4,"N_{#sigma}^{TOF}");
  fQA->SetVarTitle(5,"#Delta#phi (V0A)");
  fQA->SetVarTitle(6,"TOF PID");

  fQA2  = new AliCFContainer("fQA2v2","centr:pt:prob:TPCsig:TOFsig:DeltaPhi:maskPID",7,8,binsQA);
  fQA2->SetBinLimits(0,binCentrTOF);
  fQA2->SetBinLimits(1,binsPtTOF);
  fQA2->SetBinLimits(2,binProbTOF);
  fQA2->SetBinLimits(3,binDETsignal);
  fQA2->SetBinLimits(4,binDETsignal);
  fQA2->SetBinLimits(5,binDeltaPhi);
  fQA2->SetBinLimits(6,binMaskPID);
  fQA2->SetVarTitle(0,"centrality");
  fQA2->SetVarTitle(1,"p_{t} (GeV/c)");
  fQA2->SetVarTitle(2,"Bayesian Probability");
  fQA2->SetVarTitle(3,"N_{#sigma}^{TPC}");
  fQA2->SetVarTitle(4,"N_{#sigma}^{TOF}");
  fQA2->SetVarTitle(5,"#Delta#phi (V0C)");
  fQA2->SetVarTitle(6,"TOF PID");

  // v3
  const Int_t nDeltaPhiV3 = 7;
  Double_t binDeltaPhiV3[nDeltaPhiV3+1];
  for(Int_t i=0;i<nDeltaPhiV3+1;i++){
    binDeltaPhiV3[i] = -TMath::Pi() + i*2*TMath::Pi() / nDeltaPhiV3;
  }

  Int_t binsQAv3[7] = {nCentrTOF,nPtBinsTOF,nProbTOF,nDETsignal,nDETsignal,nDeltaPhiV3,nMaskPID};

  fQAv3  = new AliCFContainer("fQAv3","centr:pt:prob:TPCsig:TOFsig:DeltaPhi:maskPID",7,7,binsQAv3);
  fQAv3->SetBinLimits(0,binCentrTOF);
  fQAv3->SetBinLimits(1,binsPtTOF);
  fQAv3->SetBinLimits(2,binProbTOF);
  fQAv3->SetBinLimits(3,binDETsignal);
  fQAv3->SetBinLimits(4,binDETsignal);
  fQAv3->SetBinLimits(5,binDeltaPhiV3);
  fQAv3->SetBinLimits(6,binMaskPID);
  fQAv3->SetVarTitle(0,"centrality");
  fQAv3->SetVarTitle(1,"p_{t} (GeV/c)");
  fQAv3->SetVarTitle(2,"Bayesian Probability");
  fQAv3->SetVarTitle(3,"N_{#sigma}^{TPC}");
  fQAv3->SetVarTitle(4,"N_{#sigma}^{TOF}");
  fQAv3->SetVarTitle(5,"#Delta#phi (V0A)");
  fQAv3->SetVarTitle(6,"TOF PID");

  fQA2v3  = new AliCFContainer("fQA2v3","centr:pt:prob:TPCsig:TOFsig:DeltaPhi:maskPID",7,8,binsQA);
  fQA2v3->SetBinLimits(0,binCentrTOF);
  fQA2v3->SetBinLimits(1,binsPtTOF);
  fQA2v3->SetBinLimits(2,binProbTOF);
  fQA2v3->SetBinLimits(3,binDETsignal);
  fQA2v3->SetBinLimits(4,binDETsignal);
  fQA2v3->SetBinLimits(5,binDeltaPhiV3);
  fQA2v3->SetBinLimits(6,binMaskPID);
  fQA2v3->SetVarTitle(0,"centrality");
  fQA2v3->SetVarTitle(1,"p_{t} (GeV/c)");
  fQA2v3->SetVarTitle(2,"Bayesian Probability");
  fQA2v3->SetVarTitle(3,"N_{#sigma}^{TPC}");
  fQA2v3->SetVarTitle(4,"N_{#sigma}^{TOF}");
  fQA2v3->SetVarTitle(5,"#Delta#phi (V0C)");
  fQA2v3->SetVarTitle(6,"TOF PID");


  fList->Add(fPhiRPv0A);
  fList->Add(fPhiRPv0C);
  //  fList->Add(fPhiTracks); // comment if not needed
  fList->Add(fQA);
  fList->Add(fQA2);

  fList2->Add(fPhiRPv0Av3);
  fList2->Add(fPhiRPv0Cv3);
  fList2->Add(fQAv3);
  fList2->Add(fQA2v3);

  //  fList->Add(fTree); // comment if not needed

  printf("Output creation ok!!\n\n\n\n");

  // Post output data.
  if(fV2) PostData(1, fList);
  if(fV3) PostData(2, fList2);
}

//______________________________________________________________________________
void AliAnalysisTaskVnV0::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
	Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
	this->Dump();
	return;
    }
    
    Int_t run = fAOD->GetRunNumber();

    if(run != fRun){
	// Load the calibrations run dependent
	OpenInfoCalbration(run);
	fRun=run;
    }

    Float_t zvtx = GetVertex(fAOD);

    if (TMath::Abs(zvtx) < fVtxCut) {
      //Centrality
      Float_t v0Centr  = -10.;
      Float_t trkCentr  = -10.;
      AliCentrality *centrality = fAOD->GetCentrality();
      if (centrality){
	v0Centr  = centrality->GetCentralityPercentile("V0M");
	trkCentr = centrality->GetCentralityPercentile("TRK"); 
      }

      if(TMath::Abs(v0Centr - trkCentr) < 5.0){ // consistency cut on centrality selection
	fPID->SetDetResponse(fAOD, v0Centr); // Set the PID object for each event!!!!
	Analyze(fAOD,v0Centr); // Do analysis!!!!

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
    
    //reset Q vector info	
    Double_t Qxa2 = 0, Qya2 = 0;
    Double_t Qxc2 = 0, Qyc2 = 0;
    Double_t Qxa3 = 0, Qya3 = 0;
    Double_t Qxc3 = 0, Qyc3 = 0;
    
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
    Double_t Qxamean2 = fMeanQ[iC][1][0];
    Double_t Qxarms2  = fWidthQ[iC][1][0];
    Double_t Qyamean2 = fMeanQ[iC][1][1];
    Double_t Qyarms2  = fWidthQ[iC][1][1];
    Double_t Qxamean3 = fMeanQv3[iC][1][0];
    Double_t Qxarms3  = fWidthQv3[iC][1][0];
    Double_t Qyamean3 = fMeanQv3[iC][1][1];
    Double_t Qyarms3  = fWidthQv3[iC][1][1];
    
    Double_t Qxcmean2 = fMeanQ[iC][0][0];
    Double_t Qxcrms2  = fWidthQ[iC][0][0];
    Double_t Qycmean2 = fMeanQ[iC][0][1];
    Double_t Qycrms2  = fWidthQ[iC][0][1];	
    Double_t Qxcmean3 = fMeanQv3[iC][0][0];
    Double_t Qxcrms3  = fWidthQv3[iC][0][0];
    Double_t Qycmean3 = fMeanQv3[iC][0][1];
    Double_t Qycrms3  = fWidthQv3[iC][0][1];	
    
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
				 
    Int_t nAODTracks = aodEvent->GetNumberOfTracks();
    //loop track and get pid
    for(Int_t iT = 0; iT < nAODTracks; iT++) { // loop on the tracks
      AliAODTrack* aodTrack = aodEvent->GetTrack(iT);
	
      if (!aodTrack){
	aodTrack->Delete();
	continue;
      }
      
      Bool_t trkFlag = aodTrack->TestFilterBit(1); // TPC only tracks

      if ((TMath::Abs(aodTrack->Eta()) > fEtaCut) || (aodTrack->Pt() < fMinPt) || (aodTrack->GetTPCNcls() < 70) || !trkFlag){
	continue;
      }

      Double_t b[2] = {-99., -99.};
      Double_t bCov[3] = {-99., -99., -99.};
      if (!aodTrack->PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 100., b, bCov))
	continue;
	    
      if ((TMath::Abs(b[0]) > 3.0) || (TMath::Abs(b[1]) > 2.4))
	continue;
 	    
      // re-map the container in an array to do the analysis for V0A and V0C within a loop
      Float_t evPlAngV0[2] = {evPlAngV0ACor2,evPlAngV0CCor2};
      AliCFContainer *contV0[2] = {fContAllChargesV0A,fContAllChargesV0C};
      AliCFContainer *QA[2] = {fQA,fQA2};

      Float_t evPlAngV0v3[2] = {evPlAngV0ACor3,evPlAngV0CCor3};
      AliCFContainer *contV0v3[2] = {fContAllChargesV0Av3,fContAllChargesV0Cv3};
      AliCFContainer *QAv3[2] = {fQAv3,fQA2v3};

      for(Int_t iV0=0;iV0<2;iV0++){ // loop on A and C side

	fPID->SetPsiCorrectionDeDx(evPlAngV0[iV0],evPlRes[iV0*8+iC]); // set the PID dE/dx correction as a function of the v2-EP (resolution is needed)

	Double_t v2V0 = TMath::Cos(2*(aodTrack->Phi() - evPlAngV0[iV0]));
	Double_t v3V0 = TMath::Cos(3*(aodTrack->Phi() - evPlAngV0v3[iV0]));
	    
	fPID->ComputeProb(aodTrack); // compute Bayesian probabilities
	Float_t *probRead = fPID->GetProb();
	Float_t prob[8] = {probRead[0],probRead[1],probRead[2],probRead[3],probRead[4],probRead[5],probRead[6],probRead[7]};
	Float_t tofMismProb = fPID->GetTOFMismProb(); // TOF mismatch probability requested to be lower than 50% for TOF analysis 
	Double_t x[7] = {iC,v2V0,aodTrack->Charge(),aodTrack->Pt(),1,evPlAngV0[iV0],fPID->GetCurrentMask(1)&&tofMismProb < 0.5}; // to fill analysis v2 container
	Double_t x3[7] = {iC,v3V0,aodTrack->Charge(),aodTrack->Pt(),1,evPlAngV0v3[iV0],fPID->GetCurrentMask(1)&&tofMismProb < 0.5}; // to fill analysis v3 container

	Double_t phi[5] = {iC,aodTrack->Charge(),aodTrack->Pt(),1,aodTrack->Phi()}; // to fill track container 

	// Fill no PID
	if(iV0 && fV2) fPhiTracks->Fill(phi,0);
	if(fV2) contV0[iV0]->Fill(x,0);
	if(fV3) contV0v3[iV0]->Fill(x3,0);
		

	Float_t dedx = 0;	  
	Double_t dedxExp[8];
	Float_t tof = -1;
	Double_t inttimes[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
	Double_t expTOFsigma[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
	if(aodTrack->GetDetPid()){ // check the PID object is available
	  for(Int_t iS=0;iS < 8;iS++)
	    dedxExp[iS] = fPID->GetExpDeDx(aodTrack,iS);
	  
	  dedx = aodTrack->GetTPCsignal();
		
	  if(fPID->GetCurrentMask(1)){ // if TOF is present
	    Float_t ptrack = aodTrack->P();
	    tof = aodTrack->GetTOFsignal() - fPID->GetESDpid()->GetTOFResponse().GetStartTime(ptrack);
	    aodTrack->GetIntegratedTimes(inttimes);
	    
	    for(Int_t iS=5;iS < 8;iS++) // extra info for light nuclei
	      inttimes[iS] = inttimes[0] / ptrack * mass[iS] * TMath::Sqrt(1+ptrack*ptrack/mass[iS]/mass[iS]);
	    
	    for(Int_t iS=0;iS<8;iS++) expTOFsigma[iS] = fPID->GetESDpid()->GetTOFResponse().GetExpectedSigma(ptrack, inttimes[iS], mass[iS]);
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
	Double_t xQA[7] = {iC,aodTrack->Pt(), 0.0, 4.99, 4.99,deltaPhiV0,x[6]}; // v2
	Double_t xQA3[7] = {iC,aodTrack->Pt(), 0.0, 4.99, 4.99,deltaPhiV0v3,x[6]}; // v3

	//pid selection
	if(!(fPID->GetCurrentMask(0)) || !aodTrack->GetDetPid()){} // TPC PID and PID object strictly required (very important!!!!)
	else if(prob[2] > 0.6){ // pi
	  phi[3] = prob[2]; // set probability in the container variables
	  x[4] = prob[2];
	  xQA[2] = prob[2];
	  x3[4] = x[4];
	  xQA3[2] = xQA[2];
	  if(dedx > 10.){ // set TPC signal in the QA container variables
	    xQA[3] = (dedx - dedxExp[2])/(dedxExp[2]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(x[6] > 0.5){ // set TOF signal in the QA container variables if present
	    xQA[4] = (tof - inttimes[2])/expTOFsigma[2]; // TOF
	    xQA3[4] = xQA[4]; // TOF
	  }
	  if(TMath::Abs(xQA[3]) < 5){ // TPC 5 sigma extra cut to accept the track
	    if(iV0 && fV2) fPhiTracks->Fill(phi,1);
	    if(fV2) contV0[iV0]->Fill(x,1);
	    if(fV2) QA[iV0]->Fill(xQA,0);
	    if(fV3) contV0v3[iV0]->Fill(x3,1);
	    if(fV3) QAv3[iV0]->Fill(xQA3,0);
	  }
	}
	else if(prob[3] > 0.6){ // K
	  phi[3] = prob[3];
	  x[4] = prob[3];
	  xQA[2] = prob[3];
	  x3[4] = x[4];
	  xQA3[2] = xQA[2];
	  if(dedx > 10.){
	    xQA[3] = (dedx - dedxExp[3])/(dedxExp[3]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(x[6] > 0.5){
	    xQA[4] = (tof - inttimes[3])/expTOFsigma[3]; // TOF
	    xQA3[4] = xQA[4]; // TOF
	  }
	  if(TMath::Abs(xQA[3]) < 5){
	    if(iV0 && fV2) fPhiTracks->Fill(phi,2);
	    if(fV2) contV0[iV0]->Fill(x,2);
	    if(fV2) QA[iV0]->Fill(xQA,1);
	    if(fV3) contV0v3[iV0]->Fill(x3,2);
	    if(fV3) QAv3[iV0]->Fill(xQA3,1);
	  }
	}
	else if(prob[4] > 0.6){ // p
	  phi[3] = prob[4];
	  x[4] = prob[4];
	  xQA[2] = prob[4];
	  x3[4] = x[4];
	  xQA3[2] = xQA[2];
	  if(dedx > 10.){
	    xQA[3] = (dedx - dedxExp[4])/(dedxExp[4]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(x[6] > 0.5){
	    xQA[4] = (tof - inttimes[4])/expTOFsigma[4]; // TOF
	    xQA3[4] = xQA[4]; // TOF
	  }
	  if(TMath::Abs(xQA[3]) < 5){
	    if(iV0 && fV2) fPhiTracks->Fill(phi,3);
	    if(fV2) contV0[iV0]->Fill(x,3);
	    if(fV2) QA[iV0]->Fill(xQA,2);
	    if(fV3) contV0v3[iV0]->Fill(x3,3);
	    if(fV3) QAv3[iV0]->Fill(xQA3,2);
	  }
	}
	else if(prob[0] > 0.6){ // e
	  phi[3] = prob[0];
	  x[4] = prob[0];
	  xQA[2] = prob[0];
	  x3[4] = x[4];
	  xQA3[2] = xQA[2];
	  if(dedx > 10.){
	    xQA[3] = (dedx - dedxExp[0])/(dedxExp[0]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(x[6] > 0.5){
	    xQA[4] = (tof - inttimes[0])/expTOFsigma[0]; // TOF
	    xQA3[4] = xQA[4]; // TOF
	  }
	  if(TMath::Abs(xQA[3]) < 5){
	    if(iV0 && fV2) fPhiTracks->Fill(phi,4);
	    if(fV2) contV0[iV0]->Fill(x,4);
	    if(fV2) QA[iV0]->Fill(xQA,3);
	    if(fV3) contV0v3[iV0]->Fill(x3,4);
	    if(fV3) QAv3[iV0]->Fill(xQA3,3);
	  }
	}
	else if(prob[5] > 0.6){ // d
	  phi[3] = prob[5];
	  x[4] = prob[5];
	  xQA[2] = prob[5];
	  x3[4] = x[4];
	  xQA3[2] = xQA[2];
	  if(dedx > 10.){
	    xQA[3] = (dedx - dedxExp[5])/(dedxExp[5]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(x[6] > 0.5){
	    xQA[4] = (tof - inttimes[5])/expTOFsigma[5]; // TOF
	    xQA3[4] = xQA[4]; // TOF
	  }
	  if(TMath::Abs(xQA[3]) < 5){
	    if(iV0 && fV2) fPhiTracks->Fill(phi,5);
	    if(fV2) contV0[iV0]->Fill(x,5);
	    if(fV2) QA[iV0]->Fill(xQA,4);
	    if(fV3) contV0v3[iV0]->Fill(x3,5);
	    if(fV3) QAv3[iV0]->Fill(xQA3,4);
	  }
	}
	else if(prob[6] > 0.6){ // t
	  phi[3] = prob[6];
	  x[4] = prob[6];
	  xQA[2] = prob[6];
	  x3[4] = x[4];
	  xQA3[2] = xQA[2];
	  if(dedx > 10.){
	    xQA[3] = (dedx - dedxExp[6])/(dedxExp[6]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(x[6] > 0.5){
	    xQA[4] = (tof - inttimes[6])/expTOFsigma[6]; // TOF
	    xQA3[4] = xQA[4]; // TOF
	  }
	  if(TMath::Abs(xQA[3]) < 5){
	    if(iV0 && fV2) fPhiTracks->Fill(phi,6);
	    if(fV2) contV0[iV0]->Fill(x,6);
	    if(fV2) QA[iV0]->Fill(xQA,5);
	    if(fV3) contV0v3[iV0]->Fill(x3,6);
	    if(fV3) QAv3[iV0]->Fill(xQA3,5);
	  }
	}
	else if(prob[7] > 0.6){ // He3
	  phi[3] = prob[7];
	  phi[1] *= 2;
	  x[4] = prob[7];
	  x[3] *= 2;
	  xQA[2] = prob[7];
	  x3[3] = x[3];
	  x3[4] = x[4];
	  xQA3[2] = xQA[2];
	  if(dedx > 10.){
	    xQA[3] = (dedx - dedxExp[7])/(dedxExp[7]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(x[6] > 0.5){
	    xQA[4] = (tof - inttimes[7])/expTOFsigma[7]; // TOF
	    xQA3[4] = xQA[4]; // TOF
	  }
	  if(TMath::Abs(xQA[3]) < 5){
	    if(iV0 && fV2) fPhiTracks->Fill(phi,7);
	    if(fV2) contV0[iV0]->Fill(x,7);
	    if(fV2) QA[iV0]->Fill(xQA,6);
	    if(fV3) contV0v3[iV0]->Fill(x3,7);
	    if(fV3) QAv3[iV0]->Fill(xQA3,6);
	  }
	  x[3] *= 0.5;
	  x3[3] = x[3];
	}
	
	if(x[6] > 0.5){ // if TOF was present redo TPC stand alone PID to check the PID in the same acceptance (PID mask = 2)
	  fPID->ResetDetOR(1); // exclude TOF from PID
	  tofMismProb = 0;
	  
	  fPID->ComputeProb(aodTrack);
	  probRead = fPID->GetProb();
	  
	  fPID->SetDetOR(1); // include TOF for PID
	}
	Float_t probTPC[8] = {probRead[0],probRead[1],probRead[2],probRead[3],probRead[4],probRead[5],probRead[6],probRead[7]}; // TPC stand alone prbabilities

	//pid selection TPC S.A. with TOF matching
	x[6]*=2; // set the mask to 2 id TOF is present
	if(x[6]<1 || !(fPID->GetCurrentMask(0)) || !aodTrack->GetDetPid()){} // TPC PID S.A. PID in TOF acceptance
	else if(probTPC[2] > 0.6){ // pi
	  x[4] = probTPC[2];
	  x3[4] = x[4];
	  if(dedx > 10.){
	    xQA[3] = (dedx - dedxExp[2])/(dedxExp[2]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(TMath::Abs(xQA[3]) < 5){
	    if(fV2) contV0[iV0]->Fill(x,1);
	    if(fV3) contV0v3[iV0]->Fill(x3,1);
	  }
	}
	else if(probTPC[3] > 0.6){ // K
	  x[4] = probTPC[3];
	  x3[4] = x[4];
	  if(dedx > 10.){
	    xQA[3] = (dedx - dedxExp[3])/(dedxExp[3]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(TMath::Abs(xQA[3]) < 5){
	    if(fV2) contV0[iV0]->Fill(x,2);
	    if(fV3) contV0v3[iV0]->Fill(x3,2);
	  }
	}
	else if(probTPC[4] > 0.6){ // p
	  x[4] = probTPC[4];
	  x3[4] = x[4];
	  if(dedx > 10.){
	    xQA[3] = (dedx - dedxExp[4])/(dedxExp[4]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(TMath::Abs(xQA[3]) < 5){
	    if(fV2) contV0[iV0]->Fill(x,3);
	    if(fV3) contV0v3[iV0]->Fill(x3,3);
	  }
	}
	else if(probTPC[0] > 0.6){ // e
	  x[4] = probTPC[0];
	  x3[4] = x[4];
	  if(dedx > 10.){
	    xQA[3] = (dedx - dedxExp[0])/(dedxExp[0]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(TMath::Abs(xQA[3]) < 5){
	    if(fV2) contV0[iV0]->Fill(x,4);
	    if(fV3) contV0v3[iV0]->Fill(x3,4);
	  }
	}
	else if(probTPC[5] > 0.6){ // d
	  x[4] = probTPC[5];
	  x3[4] = x[4];
	  if(dedx > 10.){
	    xQA[3] = (dedx - dedxExp[5])/(dedxExp[5]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(TMath::Abs(xQA[3]) < 5){
	    if(fV2) contV0[iV0]->Fill(x,5);
	    if(fV3) contV0v3[iV0]->Fill(x3,5);
	  }
	}
	else if(probTPC[6] > 0.6){ // t
	  x[4] = probTPC[6];
	  x3[4] = x[4];
	  if(dedx > 10.){
	    xQA[3] = (dedx - dedxExp[6])/(dedxExp[6]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(TMath::Abs(xQA[3]) < 5){
	    if(fV2) contV0[iV0]->Fill(x,6);
	    if(fV3) contV0v3[iV0]->Fill(x3,6);
	  }
	}
	else if(probTPC[7] > 0.6){ // He3
	  x[4] = probTPC[7];
	  x3[4] = x[4];
	  x[3] *= 2;
	  x3[3] = x[3];
	  if(dedx > 10.){
	    xQA[3] = (dedx - dedxExp[7])/(dedxExp[7]*0.07); // TPC
	    xQA3[3] = xQA[3]; // TPC
	  }
	  if(TMath::Abs(xQA[3]) < 5){
	    if(fV2) contV0[iV0]->Fill(x,7);
	    if(fV3) contV0v3[iV0]->Fill(x3,7);
	  }
	  x[3] *= 0.5;
	  x3[3] = x[3];
	}
      } // end side loop
    } // end track loop

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

      if ((TMath::Abs(aodTrack->Eta()) > 0.8) || (aodTrack->Pt() < 0.2) || (aodTrack->GetTPCNcls() < 70)  || !trkFlag) 
	continue;
	
      Double_t b[2] = {-99., -99.};
      Double_t bCov[3] = {-99., -99., -99.};
      if (!aodTrack->PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 100., b, bCov))
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
	    for(Int_t i=0;i  < nCentrBin;i++){
		char namecont[100];
  		if(iside==0 && icoord==0)
		    sprintf(namecont,"hQxc2_%i",i);
		else if(iside==1 && icoord==0)
		    sprintf(namecont,"hQxa2_%i",i);
		else if(iside==0 && icoord==1)
		    sprintf(namecont,"hQyc2_%i",i);
		else if(iside==1 && icoord==1)
		    sprintf(namecont,"hQya2_%i",i);

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
		    sprintf(namecont,"hQxc3_%i",i);
		else if(iside==1 && icoord==0)
		    sprintf(namecont,"hQxa3_%i",i);
		else if(iside==0 && icoord==1)
		    sprintf(namecont,"hQyc3_%i",i);
		else if(iside==1 && icoord==1)
		    sprintf(namecont,"hQya3_%i",i);

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
