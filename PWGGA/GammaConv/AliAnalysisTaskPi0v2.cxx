#include <exception>
#include "TRandom3.h"
#include "TChain.h"
#include "AliLog.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskPi0v2.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliCentrality.h"
#include <iostream>

#include "TFile.h"
#include "AliOADBContainer.h"

// Author Daniel Lohner (Daniel.Lohner@cern.ch)

using namespace std;

ClassImp(AliAnalysisTaskPi0v2)


//________________________________________________________________________
    AliAnalysisTaskPi0v2::AliAnalysisTaskPi0v2(const char *name) : AliAnalysisTaskSE(name),
    fV0Reader(NULL),
    fConversionSelection(NULL),
    fConversionGammas(NULL),
    fNCentralityBins(1),
    fCentralityBins(NULL),
    fCentrality(-1),
    fCentralityBin(0),
    fNBinsPhi(6),
    fEP(NULL),
    fEtaMax(0.75),
    fEtaGap(1),
    fRPTPCEtaA(0),
    fRPTPCEtaC(0),
    fRPV0A(0),
    fRPV0C(0),
    fNCuts(0),
    fCutList(NULL),
    fConversionCuts(NULL),
    fMesonCutList(NULL),
    fMesonCuts(NULL),
    fRandomizer(NULL),
    fOutputList(NULL),
    fMesonPDGCode(kPi0),
    fInvMassRange(NULL),
    fDeltaPsiRP(0),
    fRunNumber(0),
    fRunIndex(0),
    fNEPMethods(knEPMethod),
    fFillQA(kTRUE),
   
    hNEvents(NULL),
    hRPTPC(NULL),
    hRPV0A(NULL),
    hRPV0C(NULL),
    hRPTPCAC(NULL),
    hRPV0ATPC(NULL),
    hRPV0CTPC(NULL),
    hRPV0AC(NULL),
    hCos2TPC(NULL),
    hCos2V0ATPC(NULL),
    hCos2V0CTPC(NULL),
    hCos2V0AC(NULL),
    hRPTPCEtaA(NULL),
    hRPTPCEtaC(NULL),
    hRPTPCEtaAC(NULL),
    hCos2TPCEta(NULL),
    hCos2TPCWeightedPhoton(NULL),
    hCos2TPCEtaWeightedPhoton(NULL),
    hCos2V0ATPCWeightedPhoton(NULL),
    hCos2V0CTPCWeightedPhoton(NULL),
    hCos2V0ACWeightedPhoton(NULL),
    hCos2TPCWeightedCharged(NULL),
    hCos2TPCEtaWeightedCharged(NULL),
    hCos2V0ATPCWeightedCharged(NULL),
    hCos2V0CTPCWeightedCharged(NULL),
    hCos2V0ACWeightedCharged(NULL),
    hCos2TPCWeightedV0Mult(NULL),
    hCos2TPCEtaWeightedV0Mult(NULL),
    hCos2V0ATPCWeightedV0Mult(NULL),
    hCos2V0CTPCWeightedV0Mult(NULL),
    hCos2V0ACWeightedV0Mult(NULL),
    hGammaMultCent(NULL),
    hGammaPhi(NULL),
    hMultChargedvsNGamma(NULL),
    hMultChargedvsVZERO(NULL),
    hMultChargedvsSPD(NULL),
    hGammadNdPhi(NULL),
    hGammaMultdPhiTRUE(NULL),
    hGammaMultdPhiRECOTRUE(NULL),
    hGammaMultTRUE(NULL),
    hGammaMultRECOTRUE(NULL),
    hGammaMultdPhi(NULL),
    hGammaMult(NULL),
    hGamma(NULL),
    hGammaFull(NULL),
    hCharged(NULL),
    hPi0(NULL),
    hPi0BG(NULL),

    fMultV0(0x0),
    fV0Cpol(0.),
    fV0Apol(0.),
    hEPVertex(NULL)

{
    fInvMassRange=new Double_t[2];
    fInvMassRange[0]=0.05;
    fInvMassRange[1]=0.3;

    fRandomizer=new TRandom3();
    fRandomizer->SetSeed(0);

    // Define input and output slots here
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskPi0v2::~AliAnalysisTaskPi0v2(){

    if(fRandomizer){
	delete fRandomizer;
	fRandomizer=NULL;
    }
    if(fCentralityBins){
	delete fCentralityBins;
	fCentralityBins=NULL;
    }
    if(fInvMassRange){
	delete fInvMassRange;
        fInvMassRange=NULL;
    }
    if(fCutList){
	delete fCutList;
	fCutList=NULL;
    }
    if(fMesonCutList){
        delete fMesonCutList;
        fMesonCutList=NULL;
    }

    if(fConversionSelection){
	delete fConversionSelection;
	fConversionSelection=NULL;
    }
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::UserCreateOutputObjects()
{
    OpenFile(1);

    // GetConversionCuts
    fConversionCuts=fV0Reader->GetConversionCuts();

    // Flags

    Bool_t IsMC=AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
    Bool_t IsHeavyIon=fConversionCuts->IsHeavyIon();

    if(!IsHeavyIon||IsMC)fNEPMethods=1;

    if(!fCentralityBins){
	fCentralityBins=new Double_t[fNCentralityBins+1];
	fCentralityBins[0]=-0.5;
	fCentralityBins[1]=0.5;
    }

    // Create histograms

    if(fOutputList != NULL){
	delete fOutputList;
	fOutputList = NULL;
    }
    if(fOutputList == NULL){
	fOutputList = new TList();
	fOutputList->SetOwner(kTRUE);
    }

    Int_t kGCnXBinsSpectra = Int_t((fInvMassRange[1]-fInvMassRange[0])*500);  //500 for range 0 - 1
    Double_t kGCfirstXBinSpectra = fInvMassRange[0];
    Double_t kGClastXBinSpectra = fInvMassRange[1];

    Int_t nbinspi0[knbinsPi0]={kGCnYBinsSpectra,kGCnXBinsSpectra,6,fNCentralityBins,fNEPMethods};
    Double_t minpi0[knbinsPi0]={kGCfirstYBinSpectra,kGCfirstXBinSpectra,0,-0.5,-0.5};
    Double_t maxpi0[knbinsPi0]={kGClastYBinSpectra,kGClastXBinSpectra,TMath::Pi()/2,fNCentralityBins-0.5,fNEPMethods-0.5};
    const char *binpi0[knbinsPi0]={"pt","mass","dPhi","centr","EPm"};

    Int_t nbinsg[knbinsGamma]={kGCnYBinsSpectra,6,fNCentralityBins,fNEPMethods};
    Double_t ming[knbinsGamma]={kGCfirstYBinSpectra,0,-0.5,-0.5};
    Double_t maxg[knbinsGamma]={kGClastYBinSpectra,TMath::Pi()/2,fNCentralityBins-0.5,fNEPMethods-0.5};
    const char *bingamma[knbinsGamma]={"pt","dPhi","centr","EPm"};

    // Define Binning

    if(!IsMC){

	hPi0=new THnSparseF*[fNCuts];
	hPi0BG=new THnSparseF*[fNCuts];
	hGamma=new THnSparseF*[fNCuts];

	// Photon Cuts
	for(Int_t iCut=0;iCut<fNCuts;iCut++){
	    AliConversionCuts *fPhotonv2Cuts=(AliConversionCuts*)fCutList->At(iCut);

	    TList *fCutOutputList=new TList();
	    fCutOutputList->SetName(fPhotonv2Cuts->GetCutNumber().Data());
	    fCutOutputList->SetOwner(kTRUE);
	    fOutputList->Add(fCutOutputList);

	    /*fPhotonv2Cuts->SetFillCutHistograms("",kFALSE);
	     TList *fCutQAList=fPhotonv2Cuts->GetCutHistograms();
	     if(fCutQAList){
	     fCutQAList->SetOwner(kTRUE);
	     fCutOutputList->Add(fCutQAList);
	     } */

	    hPi0[iCut]=new THnSparseF("Pi0_Sparse","Pi0_Sparse",knbinsPi0,nbinspi0,minpi0,maxpi0);
	    for(Int_t i=0;i<knbinsPi0;i++) hPi0[iCut]->GetAxis(i)->SetName(binpi0[i]);
	    hPi0[iCut]->Sumw2();
	    fCutOutputList->Add(hPi0[iCut]);

	    hPi0BG[iCut]=new THnSparseF("Pi0BG_Sparse","Pi0BG_Sparse",knbinsPi0,nbinspi0,minpi0,maxpi0);
	    for(Int_t i=0;i<knbinsPi0;i++) hPi0BG[iCut]->GetAxis(i)->SetName(binpi0[i]);
	    hPi0BG[iCut]->Sumw2();
	    fCutOutputList->Add(hPi0BG[iCut]);

	    // Gamma

	    hGamma[iCut]=new THnSparseF("Gamma_Sparse","Gamma_Sparse",knbinsGamma,nbinsg,ming,maxg);
	    for(Int_t i=0;i<knbinsGamma;i++) hGamma[iCut]->GetAxis(i)->SetName(bingamma[i]);
	    hGamma[iCut]->Sumw2();
	    fCutOutputList->Add(hGamma[iCut]);

        }
    }

    if(IsHeavyIon&&!IsMC){

	// RP Calculation
	TList *fRPList=new TList();
	fRPList->SetName("Event Plane");
	fRPList->SetOwner(kTRUE);
	fOutputList->Add(fRPList);

	hRPTPC=new TH2F("TPCAC" ,"TPC_AC" , fNCentralityBins,fCentralityBins, 180, 0, TMath::Pi());
	hRPTPC->Sumw2();
	fRPList->Add(hRPTPC);
	hRPTPCEtaA=new TH2F("TPCEtaA" ,"TPC_EtaA" , fNCentralityBins,fCentralityBins, 180, 0, TMath::Pi());
	hRPTPCEtaA->Sumw2();
	fRPList->Add(hRPTPCEtaA);
	hRPTPCEtaC=new TH2F("TPCEtaC" ,"TPC_EtaC" , fNCentralityBins,fCentralityBins, 180, 0, TMath::Pi());
	hRPTPCEtaC->Sumw2();
	fRPList->Add(hRPTPCEtaC);
	hRPV0A=new TH2F("V0A" ,"VZERO_A" , fNCentralityBins,fCentralityBins, 180, 0, TMath::Pi());
	hRPV0A->Sumw2();
	fRPList->Add(hRPV0A);
	hRPV0C=new TH2F("V0C" ,"VZERO_C" , fNCentralityBins,fCentralityBins, 180, 0, TMath::Pi());
	hRPV0C->Sumw2();
	fRPList->Add(hRPV0C);

	hRPTPCAC=new TH2F("TPCA_TPCC" ,"TPCA_TPCC" , 180, 0, TMath::Pi(), 180, 0, TMath::Pi());
	hRPTPCAC->Sumw2();
	fRPList->Add(hRPTPCAC);
	hRPV0ATPC=new TH2F("V0A_TPC" ,"V0A_TPC" , 180, 0, TMath::Pi(), 180, 0, TMath::Pi());
	hRPV0ATPC->Sumw2();
	fRPList->Add(hRPV0ATPC);
	hRPV0CTPC=new TH2F("V0C_TPC" ,"V0C_TPC" , 180, 0, TMath::Pi(), 180, 0, TMath::Pi());
	hRPV0CTPC->Sumw2();
	fRPList->Add(hRPV0CTPC);
	hRPV0AC=new TH2F("V0A_V0C" ,"V0A_V0C" , 180, 0, TMath::Pi(), 180, 0, TMath::Pi());
	hRPV0AC->Sumw2();
	fRPList->Add(hRPV0AC);
	hRPTPCEtaAC=new TH2F("TPCEtaA_TPCEtaC" ,"TPCEtaA_TPCEtaC" , 180, 0, TMath::Pi(), 180, 0, TMath::Pi());
	hRPTPCEtaAC->Sumw2();
	fRPList->Add(hRPTPCEtaAC);

	hCos2TPC=new TH2F("Cos2_TPCAC" ,"Cos2_TPCAC" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2TPC->Sumw2();
	fRPList->Add(hCos2TPC);
	hCos2TPCEta=new TH2F("Cos2_TPCEtaAC" ,"Cos2_TPCEtaAC" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2TPCEta->Sumw2();
	fRPList->Add(hCos2TPCEta);
	hCos2V0ATPC=new TH2F("Cos2_V0ATPC" ,"Cos2_V0ATPC" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2V0ATPC->Sumw2();
	fRPList->Add(hCos2V0ATPC);
	hCos2V0CTPC=new TH2F("Cos2_V0CTPC" ,"Cos2_V0CTPC" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2V0CTPC->Sumw2();
	fRPList->Add(hCos2V0CTPC);
	hCos2V0AC=new TH2F("Cos2_V0AC" ,"Cos2_V0AC" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2V0AC->Sumw2();
	fRPList->Add(hCos2V0AC);

	hCos2TPCWeightedPhoton=new TH2F("Cos2_TPCAC_WeightedPhoton" ,"Cos2_TPCAC_WeightedPhoton" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2TPCWeightedPhoton->Sumw2();
	fRPList->Add(hCos2TPCWeightedPhoton);
	hCos2TPCEtaWeightedPhoton=new TH2F("Cos2_TPCEtaAC_WeightedPhoton" ,"Cos2_TPCEtaAC_WeightedPhoton" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2TPCEtaWeightedPhoton->Sumw2();
	fRPList->Add(hCos2TPCEtaWeightedPhoton);
	hCos2V0ATPCWeightedPhoton=new TH2F("Cos2_V0ATPC_WeightedPhoton" ,"Cos2_V0ATPC_WeightedPhoton" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2V0ATPCWeightedPhoton->Sumw2();
	fRPList->Add(hCos2V0ATPCWeightedPhoton);
	hCos2V0CTPCWeightedPhoton=new TH2F("Cos2_V0CTPC_WeightedPhoton" ,"Cos2_V0CTPC_WeightedPhoton" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2V0CTPCWeightedPhoton->Sumw2();
	fRPList->Add(hCos2V0CTPCWeightedPhoton);
	hCos2V0ACWeightedPhoton=new TH2F("Cos2_V0AC_WeightedPhoton" ,"Cos2_V0AC_WeightedPhoton" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2V0ACWeightedPhoton->Sumw2();
	fRPList->Add(hCos2V0ACWeightedPhoton);

	hCos2TPCWeightedCharged=new TH2F("Cos2_TPCAC_WeightedCharged" ,"Cos2_TPCAC_WeightedCharged" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2TPCWeightedCharged->Sumw2();
	fRPList->Add(hCos2TPCWeightedCharged);
	hCos2TPCEtaWeightedCharged=new TH2F("Cos2_TPCEtaAC_WeightedCharged" ,"Cos2_TPCEtaAC_WeightedCharged" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2TPCEtaWeightedCharged->Sumw2();
	fRPList->Add(hCos2TPCEtaWeightedCharged);
	hCos2V0ATPCWeightedCharged=new TH2F("Cos2_V0ATPC_WeightedCharged" ,"Cos2_V0ATPC_WeightedCharged" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2V0ATPCWeightedCharged->Sumw2();
	fRPList->Add(hCos2V0ATPCWeightedCharged);
	hCos2V0CTPCWeightedCharged=new TH2F("Cos2_V0CTPC_WeightedCharged" ,"Cos2_V0CTPC_WeightedCharged" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2V0CTPCWeightedCharged->Sumw2();
	fRPList->Add(hCos2V0CTPCWeightedCharged);
	hCos2V0ACWeightedCharged=new TH2F("Cos2_V0AC_WeightedCharged" ,"Cos2_V0AC_WeightedCharged" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2V0ACWeightedCharged->Sumw2();
	fRPList->Add(hCos2V0ACWeightedCharged);

	hCos2TPCWeightedV0Mult=new TH2F("Cos2_TPCAC_WeightedV0Mult" ,"Cos2_TPCAC_WeightedV0Mult" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2TPCWeightedV0Mult->Sumw2();
	fRPList->Add(hCos2TPCWeightedV0Mult);
	hCos2TPCEtaWeightedV0Mult=new TH2F("Cos2_TPCEtaAC_WeightedV0Mult" ,"Cos2_TPCEtaAC_WeightedV0Mult" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2TPCEtaWeightedV0Mult->Sumw2();
	fRPList->Add(hCos2TPCEtaWeightedV0Mult);
	hCos2V0ATPCWeightedV0Mult=new TH2F("Cos2_V0ATPC_WeightedV0Mult" ,"Cos2_V0ATPC_WeightedV0Mult" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2V0ATPCWeightedV0Mult->Sumw2();
	fRPList->Add(hCos2V0ATPCWeightedV0Mult);
	hCos2V0CTPCWeightedV0Mult=new TH2F("Cos2_V0CTPC_WeightedV0Mult" ,"Cos2_V0CTPC_WeightedV0Mult" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2V0CTPCWeightedV0Mult->Sumw2();
	fRPList->Add(hCos2V0CTPCWeightedV0Mult);
	hCos2V0ACWeightedV0Mult=new TH2F("Cos2_V0AC_WeightedV0Mult" ,"Cos2_V0AC_WeightedV0Mult" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hCos2V0ACWeightedV0Mult->Sumw2();
	fRPList->Add(hCos2V0ACWeightedV0Mult);

        const Int_t nepbins=4;
	Int_t nbinsep[nepbins]={30,30,40,180};
	Double_t minep[nepbins]={-0.015,0.17,-10,0};
	Double_t maxep[nepbins]={0.015,0.20,10,TMath::Pi()};

	hEPVertex=new THnSparseF("EP_Vertex","EP_Vertex",nepbins,nbinsep,minep,maxep);
	fRPList->Add(hEPVertex);

	
    }

    TList *fPhotonQAList=new TList();
    fPhotonQAList->SetName("Gamma_QA");
    fPhotonQAList->SetOwner(kTRUE);
    fOutputList->Add(fPhotonQAList);

    if(fFillQA){
	// Gamma QA
	hGammaPhi=new TH2F*[fNCentralityBins];
	for(Int_t m=0;m<fNCentralityBins;m++){
	    hGammaPhi[m]=new TH2F(Form("%d_GammaPhi",m),"GammaPhi",kGCnYBinsSpectra,kGCfirstYBinSpectra,kGClastYBinSpectra,360,0,2*TMath::Pi());
	    hGammaPhi[m]->Sumw2();
	    fPhotonQAList->Add(hGammaPhi[m]);
	}

	hGammaMultCent=new TH2F("GammaMultvsCent","GammaMultvsCent",fNCentralityBins,fCentralityBins, 60,-0.5,59.5);
	hGammaMultCent->Sumw2();
	fPhotonQAList->Add(hGammaMultCent);

	hMultChargedvsSPD=new TH2F("Mult_ChargedvsSPD","Mult_ChargedvsSPD",250,0,2500, 250,0,5000);
	hMultChargedvsSPD->Sumw2();
	fPhotonQAList->Add(hMultChargedvsSPD);
	hMultChargedvsVZERO=new TH2F("Mult_ChargedvsVZERO","Mult_ChargedvsVZERO",250,0,2500, 200,0,20000);
	hMultChargedvsVZERO->Sumw2();
	fPhotonQAList->Add(hMultChargedvsVZERO);
	hMultChargedvsNGamma=new TH2F("Mult_ChargedvsNGamma","Mult_ChargedvsNGamma",250,0,2500,60,-0.5,59.5);
	hMultChargedvsNGamma->Sumw2();
	fPhotonQAList->Add(hMultChargedvsNGamma);

	Int_t nbinsgmult[knbinsGammaMult]={kGCnYBinsSpectra,400,fNCentralityBins};
	Double_t mingmult[knbinsGammaMult]={kGCfirstYBinSpectra,0,-0.5};
	Double_t maxgmult[knbinsGammaMult]={kGClastYBinSpectra,8000,fNCentralityBins-0.5};
	Double_t maxgmultdPhi[knbinsGammaMult]={kGClastYBinSpectra,2000,fNCentralityBins-0.5};
	const char *bingammamult[knbinsGammaMult]={"pt","gmult","centr"};

	if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){

	    hGammaMultdPhiTRUE=new THnSparseF("Gamma_MultdPhi_TRUE","Gamma_MultdPhi_TRUE",knbinsGammaMult,nbinsgmult,mingmult,maxgmultdPhi);
	    for(Int_t i=0;i<knbinsGammaMult;i++) hGammaMultdPhiTRUE->GetAxis(i)->SetName(bingammamult[i]);
	    hGammaMultdPhiTRUE->Sumw2();
	    fPhotonQAList->Add(hGammaMultdPhiTRUE);

	hGammaMultTRUE=new THnSparseF("Gamma_Mult_TRUE","Gamma_Mult_TRUE",knbinsGammaMult,nbinsgmult,mingmult,maxgmult);
	for(Int_t i=0;i<knbinsGammaMult;i++) hGammaMultTRUE->GetAxis(i)->SetName(bingammamult[i]);
	hGammaMultTRUE->Sumw2();
	fPhotonQAList->Add(hGammaMultTRUE);

	hGammaMultdPhiRECOTRUE=new THnSparseF("Gamma_MultdPhi_RECOTRUE","Gamma_MultdPhi_RECOTRUE",knbinsGammaMult,nbinsgmult,mingmult,maxgmultdPhi);
	for(Int_t i=0;i<knbinsGammaMult;i++) hGammaMultdPhiRECOTRUE->GetAxis(i)->SetName(bingammamult[i]);
	hGammaMultdPhiRECOTRUE->Sumw2();
	fPhotonQAList->Add(hGammaMultdPhiRECOTRUE);

	hGammaMultRECOTRUE=new THnSparseF("Gamma_Mult_RECOTRUE","Gamma_Mult_RECOTRUE",knbinsGammaMult,nbinsgmult,mingmult,maxgmult);
	for(Int_t i=0;i<knbinsGammaMult;i++) hGammaMultRECOTRUE->GetAxis(i)->SetName(bingammamult[i]);
	hGammaMultRECOTRUE->Sumw2();
	fPhotonQAList->Add(hGammaMultRECOTRUE);
	}

	hGammaMultdPhi=new THnSparseF*[fNEPMethods];
	hGammaMult=new THnSparseF*[fNEPMethods];

	hGammadNdPhi=new THnSparseF("Gamma_dNdPhi","Gamma_dNdPhi",knbinsGamma,nbinsg,ming,maxg);
	for(Int_t i=0;i<knbinsGamma;i++) hGammadNdPhi->GetAxis(i)->SetName(bingamma[i]);
	hGammadNdPhi->Sumw2();
	fPhotonQAList->Add(hGammadNdPhi);

	for(Int_t iEP=0;iEP<fNEPMethods;iEP++){
	    hGammaMultdPhi[iEP]=new THnSparseF(Form("Gamma_MultdPhi_%d",iEP),"Gamma_MultdPhi",knbinsGammaMult,nbinsgmult,mingmult,maxgmultdPhi);
	    for(Int_t i=0;i<knbinsGammaMult;i++) hGammaMultdPhi[iEP]->GetAxis(i)->SetName(bingammamult[i]);
	    hGammaMultdPhi[iEP]->Sumw2();
	    fPhotonQAList->Add(hGammaMultdPhi[iEP]);

	    hGammaMult[iEP]=new THnSparseF(Form("Gamma_Mult_%d",iEP),"Gamma_Mult",knbinsGammaMult,nbinsgmult,mingmult,maxgmult);
	    for(Int_t i=0;i<knbinsGammaMult;i++) hGammaMult[iEP]->GetAxis(i)->SetName(bingammamult[i]);
	    hGammaMult[iEP]->Sumw2();
	    fPhotonQAList->Add(hGammaMult[iEP]);
	}

	const Int_t knbinsCharged=3;
	Int_t nbinscharged[knbinsCharged]={6,fNCentralityBins,fNEPMethods};
	Double_t mincharged[knbinsCharged]={0,-0.5,-0.5};
	Double_t maxcharged[knbinsCharged]={TMath::Pi()/2,fNCentralityBins-0.5,fNEPMethods-0.5};
	hCharged=new THnSparseF("Charged","Charged",knbinsCharged,nbinscharged,mincharged,maxcharged);
	hCharged->GetAxis(0)->SetName("dPhi");
	hCharged->GetAxis(1)->SetName("centr");
	hCharged->GetAxis(2)->SetName("EPm");
	hCharged->Sumw2();
	fPhotonQAList->Add(hCharged);

	Int_t nbinsgfull[knbinsGamma]={kGCnYBinsSpectra,24,fNCentralityBins,fNEPMethods};
	Double_t mingfull[knbinsGamma]={kGCfirstYBinSpectra,0,-0.5,-0.5};
	Double_t maxgfull[knbinsGamma]={kGClastYBinSpectra,2*TMath::Pi(),fNCentralityBins-0.5,fNEPMethods-0.5};
	hGammaFull=new THnSparseF("Gamma_Sparse_Full","Gamma_Sparse_Full",knbinsGamma,nbinsgfull,mingfull,maxgfull);
	for(Int_t i=0;i<knbinsGamma;i++) hGammaFull->GetAxis(i)->SetName(bingamma[i]);
	hGammaFull->Sumw2();
	fPhotonQAList->Add(hGammaFull);

    }
    hNEvents=new TH1F("NEvents","NEvents",fNCentralityBins,fCentralityBins);
    fPhotonQAList->Add(hNEvents);
   

    // V0 Reader Cuts
    TList *fV0ReaderCuts=fConversionCuts->GetCutHistograms();
    fV0ReaderCuts->SetOwner(kTRUE);
    fOutputList->Add(fV0ReaderCuts);

    PostData(1, fOutputList);

}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::InitConversionSelection(){

    fConversionSelection=new AliConversionSelection*[fNCuts];
    for(Int_t iCut=0;iCut<fNCuts;iCut++){
	AliConversionCuts *fPhotonv2Cuts=(AliConversionCuts*)fCutList->At(iCut);
	if(iCut==0)fEtaMax=fPhotonv2Cuts->GetEtaCut();
        AliConversionMesonCuts *fMesonv2Cuts=(AliConversionMesonCuts*)fMesonCutList->At(iCut);
	fConversionSelection[iCut]=new AliConversionSelection(fPhotonv2Cuts,fMesonv2Cuts);
	fConversionSelection[iCut]->SetInvMassRange(fInvMassRange);
    }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskPi0v2::InitEvent(){

    if(!fV0Reader){AliError("Error: No V0 Reader and Pi0 Reconstructor");return kFALSE;}
    if(!fV0Reader->IsEventSelected())return kFALSE;
    fConversionGammas=fV0Reader->GetReconstructedGammas();

    if(!fConversionSelection){
	InitConversionSelection();
    }

    if(!SetCentrality()){return kFALSE;}

    if(fConversionCuts->IsHeavyIon()){

	if(fRunNumber!=fInputEvent->GetRunNumber()){
	    fRunNumber=fInputEvent->GetRunNumber();
	    fRunIndex=GetRunIndex(fRunNumber);
	    OpenInfoCalibration(fRunNumber); // Load Calibration of V0 Event Plane
	}
	if(fRunIndex<0)return kFALSE;

	fEP=fInputEvent->GetEventplane();
	if(!fEP)return kFALSE;

	fRPTPCEtaA=GetTPCSubEPEta(fEtaGap/2,1);
	fRPTPCEtaC=GetTPCSubEPEta(-1,-fEtaGap/2);

	// GetV0 Event Plane
	GetV0EP(fInputEvent);
    }
    return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::UserExec(Option_t *) 
{

    if(!InitEvent())return;

    // Process Cuts
    for(Int_t iCut=0;iCut<fNCuts;iCut++){

	if(fConversionSelection[iCut]->ProcessEvent(fConversionGammas,fInputEvent,fMCEvent)){    // May only be called once per event!

	    if(!fMCEvent){

		// Process EP methods
		for(Int_t iEP=0;iEP<fNEPMethods;iEP++){

		    ProcessPi0s(iCut,EEventPlaneMethod(iEP));

		    ProcessGammas(iCut,EEventPlaneMethod(iEP));
		}
	    }
	}
    }

    // Fill N Events
    hNEvents->Fill(fCentrality);

    // EventPlaneResolution
    ProcessEventPlane();

    // QA
    if(fFillQA)ProcessQA();

    PostData(1, fOutputList);

}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::ProcessPi0s(Int_t iCut,EEventPlaneMethod iEP){

    if(!fConversionSelection[iCut])return;

    Int_t photonMultiplicity=fConversionSelection[iCut]->GetMultiplicity(fInputEvent);
    if(photonMultiplicity==0)return;

    // Process Pi0s
    Float_t weight=1/Float_t(photonMultiplicity);

    for(Int_t ii=0;ii<fConversionSelection[iCut]->GetNumberOfPi0s();ii++){

	AliAODConversionMother *pi0cand=fConversionSelection[iCut]->GetPi0(ii);

	if(!pi0cand)continue;

	Double_t val[knbinsPi0];
	val[kPi0Pt]=pi0cand->Pt();
	val[kPi0Mass]=pi0cand->M();
	val[kPi0dPhi]=GetPi0PhiwrtRP(pi0cand,iEP);
	val[kPi0Cent]=fCentralityBin;
	val[kPi0EPM]=Int_t(iEP);

	hPi0[iCut]->Fill(val,weight);
    }

    // Pi0 BG
    for(Int_t ii=0;ii<fConversionSelection[iCut]->GetNumberOfBGs();ii++){

	AliAODConversionMother *pi0cand=fConversionSelection[iCut]->GetBG(ii);

	if(!pi0cand)continue;

	Double_t val[knbinsPi0];
	val[kPi0Pt]=pi0cand->Pt();
	val[kPi0Mass]=pi0cand->M();
	val[kPi0dPhi]=GetPi0PhiwrtRP(pi0cand,iEP);
	val[kPi0Cent]=fCentralityBin;
	val[kPi0EPM]=Int_t(iEP);

	hPi0BG[iCut]->Fill(val,pi0cand->GetWeight()*weight);
    }
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::ProcessGammas(Int_t iCut,EEventPlaneMethod iEP){

    if(!fConversionSelection[iCut])return;

    Int_t photonMultiplicity=fConversionSelection[iCut]->GetMultiplicity(fInputEvent);
    if(photonMultiplicity==0)return;

    Float_t weight=1/Float_t(photonMultiplicity);

    for(Int_t ii=0;ii<fConversionSelection[iCut]->GetNumberOfPhotons();ii++){

	AliAODConversionPhoton *gamma=fConversionSelection[iCut]->GetPhoton(ii);

	Double_t val[knbinsGamma];
	val[kGammaPt]=gamma->Pt();
	val[kGammadPhi]=GetPhotonPhiwrtRP(gamma,iEP);
	val[kGammaCent]=fCentralityBin;
	val[kGammaEPM]=Int_t(iEP);

	hGamma[iCut]->Fill(val,weight);

	if(iCut==0&&fFillQA){
	    hGammadNdPhi->Fill(val);

	    Double_t EPAngle=GetEventPlaneAngle(iEP,gamma->Eta(),gamma,NULL);
	    Double_t dPhi=gamma->Phi()-EPAngle;
	    if(dPhi>=(2*TMath::Pi()))dPhi-=2*TMath::Pi();
	    if(dPhi<0)dPhi+=2*TMath::Pi();
	    val[kGammadPhi]=dPhi;
	    hGammaFull->Fill(val);
	}
    }
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::ProcessQA(){

    AliStack *fMCStack=NULL;
    if(fMCEvent)fMCStack=fMCEvent->Stack();

    // Multiplicity

    Double_t multcharged=fConversionSelection[0]->GetNumberOfChargedTracks(fInputEvent);
    Double_t multVZERO=fConversionSelection[0]->GetVZEROMult(fInputEvent);
    Double_t multSPD=fConversionSelection[0]->GetSPDMult(fInputEvent);

    hMultChargedvsNGamma->Fill(multcharged,fConversionSelection[0]->GetNumberOfPhotons());
    hMultChargedvsVZERO->Fill(multcharged,multVZERO);
    hMultChargedvsSPD->Fill(multcharged,multSPD);

    // Efficiency Purity

    Int_t photonMultiplicity=fConversionSelection[0]->GetMultiplicity(fInputEvent);
    Float_t weight=1/Float_t(photonMultiplicity);

    Double_t valdPhi[knbinsGammaMult];
    Double_t val[knbinsGammaMult];

    Int_t dNdPhi[fNBinsPhi];
    Int_t ncharged;

    for(Int_t iEP=0;iEP<fNEPMethods;iEP++){
	GetChargeddNdPhi(&dNdPhi[0],ncharged,iEP);

	// Reco
	for(Int_t ii=0;ii<fConversionSelection[0]->GetNumberOfPhotons();ii++){
	    AliAODConversionPhoton *gamma=fConversionSelection[0]->GetPhoton(ii);
            val[0]=gamma->Pt();
	    val[1]=ncharged;
	    val[2]=fCentralityBin;

	    valdPhi[0]=gamma->Pt();
	    valdPhi[1]=dNdPhi[GetPhotonPhiBin(gamma,iEP)];
	    valdPhi[2]=fCentralityBin;
       
	    hGammaMult[iEP]->Fill(val,weight);
	    hGammaMultdPhi[iEP]->Fill(valdPhi,weight);

	    // Gamma Phi
	    hGammaPhi[fCentralityBin]->Fill(gamma->Pt(),gamma->Phi());
	    hGammaMultCent->Fill(fCentrality,Float_t(fConversionSelection[0]->GetNumberOfPhotons()));

	    if(fMCStack){
		if(gamma->IsTruePhoton(fMCStack)){
		    hGammaMultRECOTRUE->Fill(val,weight);
		    hGammaMultdPhiRECOTRUE->Fill(valdPhi,weight);
		}
	    }
	}

	// MC Truth
	if(fMCEvent){
	    for(Int_t i = 0; i < fMCStack->GetNprimary(); i++) {
		TParticle* particle = (TParticle *)fMCStack->Particle(i);
		if (!particle) continue;
		if(fConversionCuts->PhotonIsSelectedMC(particle,fMCStack)){
		    TParticle *daughter=(TParticle *)fMCStack->Particle(particle->GetDaughter(0));
		    if(daughter){
			val[0]=particle->Pt();
			val[1]=ncharged;
			val[2]=fCentralityBin;
       
			valdPhi[0]=particle->Pt();
			valdPhi[1]=dNdPhi[GetPhiBin(GetMCPhotonPhiwrtRP(particle,EEventPlaneMethod(iEP)))];
			valdPhi[2]=fCentralityBin;
       
			hGammaMultTRUE->Fill(val,weight);
			hGammaMultdPhiTRUE->Fill(valdPhi,weight);
		    }
		}
	    }
	}
    }
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::GetPhotondNdPhi(Int_t *dNdPhi,Int_t iEP,Int_t iCut){

    for(Int_t iPhi=0;iPhi<fNBinsPhi;iPhi++)dNdPhi[iPhi]=0;

    for(Int_t ii=0;ii<fConversionSelection[iCut]->GetNumberOfPhotons();ii++){
	AliAODConversionPhoton *gamma=fConversionSelection[iCut]->GetPhoton(ii);
	Int_t phibin=GetPhotonPhiBin(gamma,iEP);
	dNdPhi[phibin]++;
    }
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::GetChargeddNdPhi(Int_t *dNdPhi,Int_t &ntot,Int_t iEP){

    for(Int_t iPhi=0;iPhi<fNBinsPhi;iPhi++)dNdPhi[iPhi]=0;
    ntot=0;

    Double_t val[3];
    val[1]=fCentralityBin;
    val[2]=Int_t(iEP);
    
    for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
	AliVTrack* currentTrack = (AliVTrack*)fInputEvent->GetTrack(iTracks);
	if(!currentTrack) continue;
	if(TMath::Abs(currentTrack->Eta())>fEtaMax)continue;

	Double_t phiwrt=GetChargedPhiwrtRP(currentTrack,EEventPlaneMethod(iEP));
	Int_t phibin=GetPhiBin(phiwrt);
     
	val[0]=phiwrt;
	hCharged->Fill(val);

	dNdPhi[phibin]++;
	ntot++;
    }
}

//________________________________________________________________________
Int_t AliAnalysisTaskPi0v2::GetPhiBin(Double_t phiwrt){
    Double_t binrange=TMath::Pi()/(2*Double_t(fNBinsPhi));
    for(Int_t iPhi=0;iPhi<fNBinsPhi;iPhi++){
	if(phiwrt>=(binrange*iPhi)&&phiwrt<(binrange*(iPhi+1)))return iPhi;
    }
    return -1;
}

//________________________________________________________________________
Int_t AliAnalysisTaskPi0v2::GetPhotonPhiBin(AliAODConversionPhoton *gamma,Int_t iEP){
    Double_t phiwrt=GetPhotonPhiwrtRP(gamma,EEventPlaneMethod(iEP));
    return GetPhiBin(phiwrt);
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetPi0PhiwrtRP(AliAODConversionMother *pi0,EEventPlaneMethod iEP){

    AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(pi0->GetLabel1()));
    AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(pi0->GetLabel2()));

    Double_t EPAngle=GetEventPlaneAngle(iEP,pi0->Eta(),gamma0,gamma1);

    return TMath::ACos(TMath::Abs(TMath::Cos(pi0->Phi()-EPAngle)));
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetPhotonPhiwrtRP(AliAODConversionPhoton *gamma,EEventPlaneMethod iEP){

    Double_t EPAngle=GetEventPlaneAngle(iEP,gamma->Eta(),gamma,NULL);

    return TMath::ACos(TMath::Abs(TMath::Cos(gamma->Phi()-EPAngle)));
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetMCPhotonPhiwrtRP(TParticle *gamma,EEventPlaneMethod iEP){

    Double_t EPAngle=GetEventPlaneAngle(iEP,gamma->Eta(),NULL,NULL);

    return TMath::ACos(TMath::Abs(TMath::Cos(gamma->Phi()-EPAngle)));
}
//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetChargedPhiwrtRP(AliVTrack *track,EEventPlaneMethod iEP){

    Double_t EPAngle=GetEventPlaneAngle(iEP,track->Eta(),NULL,NULL);

    return TMath::ACos(TMath::Abs(TMath::Cos(track->Phi()-EPAngle)));
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::Terminate(Option_t *) 
{

}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::ProcessEventPlane()
{
    if(!fMCEvent&&fConversionCuts->IsHeavyIon()){

        Double_t val[4];
	val[0]=fInputEvent->GetPrimaryVertex()->GetX();
         val[1]=fInputEvent->GetPrimaryVertex()->GetY();
         val[2]=fInputEvent->GetPrimaryVertex()->GetZ();
	 val[3]=GetEventPlaneAngle(kTPC);
         hEPVertex->Fill(val);

	// TPC EP
	Double_t PsiRP1=fEP->GetQsub1()->Phi()/2;
	Double_t PsiRP2=fEP->GetQsub2()->Phi()/2;
	Double_t EPTPC=GetEventPlaneAngle(kTPC);

	hRPTPC->Fill(fCentrality,EPTPC);
	hRPTPCAC->Fill(PsiRP1,PsiRP2);
	hCos2TPC->Fill(fCentrality,TMath::Cos(2*(PsiRP1-PsiRP2)));

	// TPC Eta Gap
	hRPTPCEtaA->Fill(fCentrality,fRPTPCEtaA);
        hRPTPCEtaC->Fill(fCentrality,fRPTPCEtaC);
	hRPTPCEtaAC->Fill(fRPTPCEtaA,fRPTPCEtaC);
	hCos2TPCEta->Fill(fCentrality,TMath::Cos(2.*(fRPTPCEtaA-fRPTPCEtaC)));

	// V0
	hCos2V0ATPC->Fill(fCentrality,TMath::Cos(2*(EPTPC-fRPV0A)));
	hCos2V0CTPC->Fill(fCentrality,TMath::Cos(2*(EPTPC-fRPV0C)));
	hCos2V0AC->Fill(fCentrality,TMath::Cos(2*(fRPV0C-fRPV0A)));

	hRPV0A->Fill(fCentrality,fRPV0A);
	hRPV0C->Fill(fCentrality,fRPV0C);

	hRPV0ATPC->Fill(fRPV0A,EPTPC);
	hRPV0CTPC->Fill(fRPV0C,EPTPC);
	hRPV0AC->Fill(fRPV0A,fRPV0C);

        // Weight with Photon Mult
	Float_t weightphoton=Float_t(fConversionSelection[0]->GetNumberOfPhotons());

	hCos2TPCWeightedPhoton->Fill(fCentrality,TMath::Cos(2*(PsiRP1-PsiRP2)),weightphoton);
	hCos2TPCEtaWeightedPhoton->Fill(fCentrality,TMath::Cos(2.*(fRPTPCEtaA-fRPTPCEtaC)),weightphoton);
	hCos2V0ATPCWeightedPhoton->Fill(fCentrality,TMath::Cos(2*(EPTPC-fRPV0A)),weightphoton);
	hCos2V0CTPCWeightedPhoton->Fill(fCentrality,TMath::Cos(2*(EPTPC-fRPV0C)),weightphoton);
	hCos2V0ACWeightedPhoton->Fill(fCentrality,TMath::Cos(2*(fRPV0C-fRPV0A)),weightphoton);

	// Weight with charged Track Mult

	Float_t weightcharged=Float_t(fConversionSelection[0]->GetNumberOfChargedTracks(fInputEvent));

        hCos2TPCWeightedCharged->Fill(fCentrality,TMath::Cos(2*(PsiRP1-PsiRP2)),weightcharged);
	hCos2TPCEtaWeightedCharged->Fill(fCentrality,TMath::Cos(2.*(fRPTPCEtaA-fRPTPCEtaC)),weightcharged);
	hCos2V0ATPCWeightedCharged->Fill(fCentrality,TMath::Cos(2*(EPTPC-fRPV0A)),weightcharged);
	hCos2V0CTPCWeightedCharged->Fill(fCentrality,TMath::Cos(2*(EPTPC-fRPV0C)),weightcharged);
	hCos2V0ACWeightedCharged->Fill(fCentrality,TMath::Cos(2*(fRPV0C-fRPV0A)),weightcharged);

	// Weight with V0 mult

	Float_t weightv0Mult=Float_t(fConversionSelection[0]->GetVZEROMult(fInputEvent));

        hCos2TPCWeightedV0Mult->Fill(fCentrality,TMath::Cos(2*(PsiRP1-PsiRP2)),weightv0Mult);
	hCos2TPCEtaWeightedV0Mult->Fill(fCentrality,TMath::Cos(2.*(fRPTPCEtaA-fRPTPCEtaC)),weightv0Mult);
	hCos2V0ATPCWeightedV0Mult->Fill(fCentrality,TMath::Cos(2*(EPTPC-fRPV0A)),weightv0Mult);
	hCos2V0CTPCWeightedV0Mult->Fill(fCentrality,TMath::Cos(2*(EPTPC-fRPV0C)),weightv0Mult);
	hCos2V0ACWeightedV0Mult->Fill(fCentrality,TMath::Cos(2*(fRPV0C-fRPV0A)),weightv0Mult);

    }
}

//________________________________________________________________________
TVector2 AliAnalysisTaskPi0v2::GetEPContribution(AliAODConversionPhoton *gamma){
    TVector2 q;
    for(Int_t ii=0;ii<2;ii++){
	AliVTrack *fCurrentTrack=AliConversionCuts::GetTrack(fInputEvent,gamma->GetTrackLabel(ii));
	TVector2 qtrack(fEP->GetQContributionX(fCurrentTrack),fEP->GetQContributionY(fCurrentTrack));
	q+=qtrack;
    }
    return q;
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetCorrectedTPCEPAngle(AliAODConversionPhoton *gamma0,AliAODConversionPhoton *gamma1){
    // Correct Event Plane for Dilepton Tracks
    TVector2 q0(*fEP->GetQVector());
    if(gamma0)q0-=GetEPContribution(gamma0);
    if(gamma1)q0-=GetEPContribution(gamma1);
    Double_t EPangle=q0.Phi()/2;
    //EPangle=ApplyFlatteningTPC(EPangle,fCentrality);

    return EPangle;
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetTPCSubEPEta(Double_t etamin,Double_t etamax){
    TVector2 q;
    for(Int_t ii=0;ii<fInputEvent->GetNumberOfTracks();ii++){
	AliVTrack *fCurrentTrack=dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(ii));
	if(!fCurrentTrack)continue;
	if(fCurrentTrack->Eta()>=etamin&&fCurrentTrack->Eta()<=etamax){
	    TVector2 qtrack(fEP->GetQContributionX(fCurrentTrack),fEP->GetQContributionY(fCurrentTrack));
	    q+=qtrack;
	}
    }
    return q.Phi()/2;
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetEventPlaneAngle(EEventPlaneMethod EPmethod,Double_t eta,AliAODConversionPhoton *gamma0,AliAODConversionPhoton *gamma1)
{
    // If arguments are not null, the contribution of these photons is subtracted from the TPC EP

    if(fConversionCuts->IsHeavyIon()){
       
	// For MC select random EP angle in order to avoid correlations due to azimuth dependent efficiencies (ITS holes)
	if(fMCEvent){
	    return fRandomizer->Uniform(2*TMath::Pi());
	}

	switch(EPmethod){
	case kTPC:
	    return GetCorrectedTPCEPAngle(gamma0,gamma1);
	case kTPCEtaGap:
	    if(eta<0)return fRPTPCEtaA; // Use opposite EP
	    else return fRPTPCEtaC;
	case kV0A:
	    return fRPV0A;
	case kV0C:
	    return fRPV0C;
	default:
	    return 0;
	}
    }

    // NO EP in pp mode
    return 0;
}

///________________________________________________________________________
Bool_t AliAnalysisTaskPi0v2::SetCentrality(){

    // Set centrality bin for current event

    if(!fConversionCuts->IsHeavyIon()){
	fCentrality=0;
        fCentralityBin=0;
	return kTRUE;
    }

    fCentrality=fConversionCuts->GetCentrality(fInputEvent);

    if(fNCentralityBins>1){
	for(fCentralityBin=0;fCentralityBin<fNCentralityBins;fCentralityBin++){
	    if(fCentrality>=fCentralityBins[fCentralityBin]&&fCentrality<fCentralityBins[fCentralityBin+1])return kTRUE;
	}
	return kFALSE;
    }
    fCentralityBin=0;
    return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::SetCentralityBins(Double_t *bins,Int_t nbins)
{
    // Set Centrality bins for analysis

    fNCentralityBins=nbins;

    if(fCentralityBins){
	delete[] fCentralityBins;
	fCentralityBins=NULL;
    }

    fCentralityBins=new Double_t[fNCentralityBins+1];
    for(Int_t ii=0;ii<=fNCentralityBins;ii++){
	fCentralityBins[ii]=bins[ii];
    }
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::SetCuts(TString *cutarray,Int_t ncuts){

    if(fCutList){
	delete fCutList;
        fCutList=NULL;
    }

    fNCuts=ncuts;

    fCutList = new TList();
    fCutList->SetOwner(kTRUE);
    AliConversionCuts **analysisCuts = new AliConversionCuts*[fNCuts];
    for(Int_t i = 0; i<fNCuts; i++){
	analysisCuts[i] = new AliConversionCuts();
	analysisCuts[i]->InitializeCutsFromCutString(cutarray[i].Data());
	fCutList->Add(analysisCuts[i]);
    }
}

//___________________________________________________________________________
Int_t AliAnalysisTaskPi0v2::GetRunIndex(Int_t run){

  switch(run){
  case  139517 : return 137;
  case  139514 : return 136;
  case  139513 : return 135; 
  case  139511 : return 134; 
  case  139510 : return 133; 
  case  139507 : return 132; 
  case  139505 : return 131; 
  case  139504 : return 130; 
  case  139503 : return 129; 
  case  139470 : return 128; 
  case  139467 : return 127; 
  case  139466 : return 126; 
  case  139465 : return 125; 
  case  139440 : return 124; 
  case  139439 : return 123; 
  case  139438 : return 122; 
  case  139437 : return 121; 
  case  139360 : return 120; 
  case  139329 : return 119; 
  case  139328 : return 118; 
  case  139314 : return 117; 
  case  139311 : return 116; 
  case  139310 : return 115; 
  case  139309 : return 114; 
  case  139308 : return 113; 
  case  139173 : return 112; 
  case  139172 : return 111; 
  case  139110 : return 110; 
  case  139107 : return 109; 
  case  139105 : return 108; 
  case  139104 : return 107; 
  case  139042 : return 106; 
  case  139038 : return 105; 
  case  139037 : return 104; 
  case  139036 : return 103; 
  case  139029 : return 102; 
  case  139028 : return 101; 
  case  138983 : return 100; 
  case  138982 : return 99; 
  case  138980 : return 98; 
  case  138979 : return 97; 
  case  138978 : return 96; 
  case  138977 : return 95; 
  case  138976 : return 94; 
  case  138973 : return 93; 
  case  138972 : return 92; 
  case  138965 : return 91; 
  case  138924 : return 90; 
  case  138872 : return 89; 
  case  138871 : return 88; 
  case  138870 : return 87; 
  case  138837 : return 86; 
  case  138830 : return 85; 
  case  138828 : return 84; 
  case  138826 : return 83; 
  case  138796 : return 82; 
  case  138795 : return 81; 
  case  138742 : return 80; 
  case  138732 : return 79; 
  case  138730 : return 78; 
  case  138666 : return 77; 
  case  138662 : return 76; 
  case  138653 : return 75; 
  case  138652 : return 74; 
  case  138638 : return 73; 
  case  138624 : return 72; 
  case  138621 : return 71; 
  case  138583 : return 70; 
  case  138582 : return 69; 
  case  138579 : return 68; 
  case  138578 : return 67; 
  case  138534 : return 66; 
  case  138469 : return 65; 
  case  138442 : return 64; 
  case  138439 : return 63; 
  case  138438 : return 62; 
  case  138396 : return 61; 
  case  138364 : return 60; 
  case  138359 : return 59; 
  case  138275 : return 58; 
  case  138225 : return 57; 
  case  138201 : return 56; 
  case  138200 : return 55; 
  case  138197 : return 54; 
  case  138192 : return 53; 
  case  138190 : return 52; 
  case  138154 : return 51; 
  case  138153 : return 50; 
  case  138151 : return 49; 
  case  138150 : return 48; 
  case  138126 : return 47; 
  case  138125 : return 46; 
  case  137848 : return 45; 
  case  137847 : return 44; 
  case  137844 : return 43; 
  case  137843 : return 42; 
  case  137752 : return 41; 
  case  137751 : return 40; 
  case  137748 : return 39; 
  case  137724 : return 38; 
  case  137722 : return 37; 
  case  137718 : return 36; 
  case  137704 : return 35; 
  case  137693 : return 34; 
  case  137692 : return 33; 
  case  137691 : return 32; 
  case  137689 : return 31; 
  case  137686 : return 30; 
  case  137685 : return 29; 
  case  137639 : return 28; 
  case  137638 : return 27; 
  case  137608 : return 26; 
  case  137595 : return 25; 
  case  137549 : return 24; 
  case  137546 : return 23; 
  case  137544 : return 22; 
  case  137541 : return 21; 
  case  137539 : return 20; 
  case  137531 : return 19; 
  case  137530 : return 18; 
  case  137443 : return 17; 
  case  137441 : return 16; 
  case  137440 : return 15; 
  case  137439 : return 14; 
  case  137434 : return 13; 
  case  137432 : return 12; 
  case  137431 : return 11; 
  case  137430 : return 10; 
  case  137366 : return 9; 
  case  137243 : return 8; 
  case  137236 : return 7; 
  case  137235 : return 6; 
  case  137232 : return 5; 
  case  137231 : return 4; 
  case  137165 : return 3; 
  case  137162 : return 2; 
  case  137161 : return 1;
  default : return -1;
  }
}

//____________________________________________________________________________
void AliAnalysisTaskPi0v2::GetV0EP(AliVEvent * event){

    // Corrected VZERO EP (from AliAnalysisTaskPi0Flow)

    //VZERO data
    AliESDVZERO* esdV0 = (AliESDVZERO*)event->GetVZEROData();

    //reset Q vector info
    Double_t Qxa2 = 0, Qya2 = 0;
    Double_t Qxc2 = 0, Qyc2 = 0;

    for (Int_t iv0 = 0; iv0 < 64; iv0++) {
	Double_t phiV0 = TMath::PiOver4()*(0.5 + iv0 % 8);
	Float_t multv0 = esdV0->GetMultiplicity(iv0);
	if (iv0 < 32){ // V0C
	    Qxc2 += TMath::Cos(2*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	    Qyc2 += TMath::Sin(2*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	} else {       // V0A
	    Qxa2 += TMath::Cos(2*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
	    Qya2 += TMath::Sin(2*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
	}
    }

    Int_t iC = -1;
    // centrality bins
    if(fCentrality < 5) iC = 0;
    else if(fCentrality < 10) iC = 1;
    else if(fCentrality < 20) iC = 2;
    else if(fCentrality < 30) iC = 3;
    else if(fCentrality < 40) iC = 4;
    else if(fCentrality < 50) iC = 5;
    else if(fCentrality < 60) iC = 6;
    else if(fCentrality < 70) iC = 7;
    else iC = 8;

    //grab for each centrality the proper histo with the Qx and Qy to do the recentering
    Double_t Qxamean2 = fMeanQ[iC][1][0];
    Double_t Qxarms2  = fWidthQ[iC][1][0];
    Double_t Qyamean2 = fMeanQ[iC][1][1];
    Double_t Qyarms2  = fWidthQ[iC][1][1];
    
    Double_t Qxcmean2 = fMeanQ[iC][0][0];
    Double_t Qxcrms2  = fWidthQ[iC][0][0];
    Double_t Qycmean2 = fMeanQ[iC][0][1];
    Double_t Qycrms2  = fWidthQ[iC][0][1];

    Double_t QxaCor2 = (Qxa2 - Qxamean2)/Qxarms2;
    Double_t QyaCor2 = (Qya2 - Qyamean2)/Qyarms2;
    Double_t QxcCor2 = (Qxc2 - Qxcmean2)/Qxcrms2;
    Double_t QycCor2 = (Qyc2 - Qycmean2)/Qycrms2;

    fRPV0A = TMath::ATan2(QyaCor2, QxaCor2)/2.;
    fRPV0C = TMath::ATan2(QycCor2, QxcCor2)/2.;

    while(fRPV0A<0)fRPV0A+=TMath::Pi() ;
    while(fRPV0A>TMath::Pi())fRPV0A-=TMath::Pi() ;
    fRPV0A=ApplyFlatteningV0A(fRPV0A,fCentrality) ;
    while(fRPV0A<0)fRPV0A+=TMath::Pi() ;
    while(fRPV0A>TMath::Pi())fRPV0A-=TMath::Pi() ;

    while(fRPV0C<0)fRPV0C+=TMath::Pi() ;
    while(fRPV0C>TMath::Pi())fRPV0C-=TMath::Pi() ;
    fRPV0C=ApplyFlatteningV0C(fRPV0C,fCentrality) ;
    while(fRPV0C<0)fRPV0C+=TMath::Pi() ;
    while(fRPV0C>TMath::Pi())fRPV0C-=TMath::Pi() ;

}


//_____________________________________________________________________________
void AliAnalysisTaskPi0v2::OpenInfoCalibration(Int_t run){
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
    printf("Setting V0 calibration \n") ;
    fMultV0 = ((TH2F *) cont->GetObject(run))->ProfileX();

    TF1 *fpol0 = new TF1("fpol0","pol0"); 
    fMultV0->Fit(fpol0,"","",0,31);
    fV0Cpol = fpol0->GetParameter(0);
    fMultV0->Fit(fpol0,"","",32,64);
    fV0Apol = fpol0->GetParameter(0);

    for(Int_t iside=0;iside<2;iside++){
	for(Int_t icoord=0;icoord<2;icoord++){
	    for(Int_t i=0;i  < nCentrBinV0;i++){
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
		//		fMeanQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
		//		fWidthQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();

	    }
	}
    }
}

//____________________________________________________________________________
Double_t  AliAnalysisTaskPi0v2::ApplyFlatteningV0A(Double_t phi, Double_t c){
  //LHC10h

  Double_t v2c =0.;
  Double_t v2s =0; 

  if(fRunIndex==1){
      v2c=0.00271301 + (0.00340046)*c + (-0.000376369)*c*c + (8.56362e-06)*c*c*c + (-5.70955e-08)*c*c*c*c ;
      v2s=0.024278 + (0.000586593)*c + (-1.07748e-06)*c*c + (-5.01567e-07)*c*c*c + (3.65956e-09)*c*c*c*c ;
  }
  if(fRunIndex==2){
      v2c=0.111374 + (-0.0111331)*c + (0.000212342)*c*c + (-9.6774e-07)*c*c*c + (-2.54199e-09)*c*c*c*c ;
      v2s=-0.0252636 + (0.0152193)*c + (-0.000830793)*c*c + (1.5099e-05)*c*c*c + (-9.04687e-08)*c*c*c*c ;
  }
  if(fRunIndex==4){
      v2c=-0.0323376 + (0.0058672)*c + (-0.000257135)*c*c + (3.54353e-06)*c*c*c + (-1.50431e-08)*c*c*c*c ;
      v2s=0.0596758 + (-0.00520203)*c + (0.000178495)*c*c + (-2.94143e-06)*c*c*c + (1.55424e-08)*c*c*c*c ;
  }
  if(fRunIndex==5){
      v2c=-0.0415369 + (0.0098273)*c + (-0.000486907)*c*c + (8.00188e-06)*c*c*c + (-4.29107e-08)*c*c*c*c ;
      v2s=-0.0146049 + (0.00693627)*c + (-0.000348663)*c*c + (5.79331e-06)*c*c*c + (-3.33179e-08)*c*c*c*c ;
  }
  if(fRunIndex==6){
      v2c=0.251902 + (-0.0415419)*c + (0.0018553)*c*c + (-3.15181e-05)*c*c*c + (1.80073e-07)*c*c*c*c ;
      v2s=0.0771839 + (-0.00962468)*c + (0.000394453)*c*c + (-7.22485e-06)*c*c*c + (4.53562e-08)*c*c*c*c ;
  }
  if(fRunIndex==7){
      v2c=0.157038 + (-0.0204377)*c + (0.000847477)*c*c + (-1.4232e-05)*c*c*c + (8.04965e-08)*c*c*c*c ;
      v2s=0.0256337 + (-0.00242333)*c + (0.000140455)*c*c + (-3.5149e-06)*c*c*c + (2.44854e-08)*c*c*c*c ;
  }
  if(fRunIndex==8){
      v2c=-0.0916108 + (0.0146996)*c + (-0.000624971)*c*c + (9.63331e-06)*c*c*c + (-4.8574e-08)*c*c*c*c ;
      v2s=0.0332113 + (-0.00332475)*c + (0.000151857)*c*c + (-2.86187e-06)*c*c*c + (1.54796e-08)*c*c*c*c ;
  }
  if(fRunIndex==9){
      v2c=0.00925348 + (1.58854e-05)*c + (-3.85178e-05)*c*c + (6.75167e-07)*c*c*c + (-3.1502e-09)*c*c*c*c ;
      v2s=0.0249253 + (-0.00300902)*c + (0.000152653)*c*c + (-3.20923e-06)*c*c*c + (2.08818e-08)*c*c*c*c ;
  }
  if(fRunIndex==10){
      v2c=0.040903 + (0.00142711)*c + (-0.000264675)*c*c + (6.50624e-06)*c*c*c + (-4.46628e-08)*c*c*c*c ;
      v2s=0.0914196 + (-0.00957901)*c + (0.000408998)*c*c + (-7.49412e-06)*c*c*c + (4.63915e-08)*c*c*c*c ;
  }
  if(fRunIndex==11){
      v2c=-0.0130133 + (0.00296474)*c + (-0.000183518)*c*c + (3.65755e-06)*c*c*c + (-2.24484e-08)*c*c*c*c ;
      v2s=-0.00449031 + (0.00365474)*c + (-0.000226455)*c*c + (4.49792e-06)*c*c*c + (-2.96018e-08)*c*c*c*c ;
  }
  if(fRunIndex==12){
      v2c=-0.021987 + (0.00506675)*c + (-0.000301768)*c*c + (5.66088e-06)*c*c*c + (-3.21447e-08)*c*c*c*c ;
      v2s=0.0830858 + (-0.00871637)*c + (0.000299212)*c*c + (-4.10705e-06)*c*c*c + (1.76888e-08)*c*c*c*c ;
  }
  if(fRunIndex==13){
      v2c=0.0502952 + (-0.00152478)*c + (-0.000106289)*c*c + (3.39398e-06)*c*c*c + (-2.45015e-08)*c*c*c*c ;
      v2s=0.0639936 + (-0.00775191)*c + (0.000347449)*c*c + (-6.46454e-06)*c*c*c + (4.00439e-08)*c*c*c*c ;
  }
  if(fRunIndex==14){
      v2c=-0.116922 + (0.0122105)*c + (-0.000373923)*c*c + (4.08464e-06)*c*c*c + (-1.05659e-08)*c*c*c*c ;
      v2s=0.128821 + (-0.0204399)*c + (0.000831009)*c*c + (-1.33147e-05)*c*c*c + (7.24876e-08)*c*c*c*c ;
  }
  if(fRunIndex==15){
      v2c=0.0853741 + (-0.00924528)*c + (0.00034257)*c*c + (-5.69551e-06)*c*c*c + (3.44998e-08)*c*c*c*c ;
      v2s=0.0205492 + (-0.00052015)*c + (-5.42051e-05)*c*c + (1.82754e-06)*c*c*c + (-1.58393e-08)*c*c*c*c ;
  }
  if(fRunIndex==16){
      v2c=-0.0233639 + (0.00703612)*c + (-0.000403314)*c*c + (7.45693e-06)*c*c*c + (-4.4522e-08)*c*c*c*c ;
      v2s=-0.0712605 + (0.0124469)*c + (-0.000543644)*c*c + (8.22618e-06)*c*c*c + (-4.1073e-08)*c*c*c*c ;
  }
  if(fRunIndex==17){
      v2c=0.0113215 + (-0.0019788)*c + (9.81834e-05)*c*c + (-1.30009e-06)*c*c*c + (5.03998e-09)*c*c*c*c ;
      v2s=0.0459851 + (-0.00170929)*c + (-5.96459e-05)*c*c + (2.90146e-06)*c*c*c + (-2.85325e-08)*c*c*c*c ;
  }
  if(fRunIndex==18){
      v2c=-0.0924521 + (0.0162203)*c + (-0.000702717)*c*c + (1.0269e-05)*c*c*c + (-4.632e-08)*c*c*c*c ;
      v2s=-0.0456605 + (0.00706138)*c + (-0.000243957)*c*c + (2.7335e-06)*c*c*c + (-8.83296e-09)*c*c*c*c ;
  }
  if(fRunIndex==19){
      v2c=0.00150974 + (0.0113356)*c + (-0.000923111)*c*c + (2.0452e-05)*c*c*c + (-1.36949e-07)*c*c*c*c ;
      v2s=0.0673845 + (-0.0187379)*c + (0.00112044)*c*c + (-2.19997e-05)*c*c*c + (1.33115e-07)*c*c*c*c ;
  }
  if(fRunIndex==20){
      v2c=0.0223891 + (-0.000921021)*c + (-3.56276e-05)*c*c + (9.95622e-07)*c*c*c + (-5.54486e-09)*c*c*c*c ;
      v2s=0.0657305 + (-0.00354234)*c + (2.90428e-05)*c*c + (5.82994e-07)*c*c*c + (-7.53656e-09)*c*c*c*c ;
  }
  if(fRunIndex==21){
      v2c=0.00969434 + (-0.00195313)*c + (3.56409e-05)*c*c + (2.15514e-07)*c*c*c + (-5.5161e-09)*c*c*c*c ;
      v2s=0.000236584 + (0.00458234)*c + (-0.000306408)*c*c + (6.04476e-06)*c*c*c + (-3.79354e-08)*c*c*c*c ;
  }
  if(fRunIndex==22){
      v2c=0.0973536 + (-0.00907174)*c + (0.000239041)*c*c + (-2.64399e-06)*c*c*c + (1.12038e-08)*c*c*c*c ;
      v2s=0.0506578 + (-0.00256521)*c + (1.53982e-05)*c*c + (4.18193e-07)*c*c*c + (-4.99622e-09)*c*c*c*c ;
  }
  if(fRunIndex==23){
      v2c=0.120714 + (-0.0190965)*c + (0.000887955)*c*c + (-1.57144e-05)*c*c*c + (9.22823e-08)*c*c*c*c ;
      v2s=-0.00178837 + (0.00377706)*c + (-0.000358619)*c*c + (9.54949e-06)*c*c*c + (-7.39141e-08)*c*c*c*c ;
  }
  if(fRunIndex==24){
      v2c=0.0593053 + (-0.00687397)*c + (0.000200607)*c*c + (-2.34483e-06)*c*c*c + (9.60642e-09)*c*c*c*c ;
      v2s=0.0235174 + (0.00219788)*c + (-0.000225065)*c*c + (4.9594e-06)*c*c*c + (-3.2695e-08)*c*c*c*c ;
  }
  if(fRunIndex==25){
      v2c=0.0551353 + (-0.00340874)*c + (2.26718e-05)*c*c + (4.1675e-07)*c*c*c + (-3.37457e-09)*c*c*c*c ;
      v2s=0.0110038 + (-0.000180543)*c + (3.28645e-05)*c*c + (-1.08119e-06)*c*c*c + (6.7157e-09)*c*c*c*c ;
  }
  if(fRunIndex==26){
      v2c=-0.0138673 + (0.00475006)*c + (-0.00030584)*c*c + (6.0294e-06)*c*c*c + (-3.69906e-08)*c*c*c*c ;
      v2s=0.0395065 + (-0.00331262)*c + (9.40592e-05)*c*c + (-1.1541e-06)*c*c*c + (4.38329e-09)*c*c*c*c ;
  }
  if(fRunIndex==27){
      v2c=0.0547504 + (-0.00669067)*c + (0.000269785)*c*c + (-4.75217e-06)*c*c*c + (2.96254e-08)*c*c*c*c ;
      v2s=0.101508 + (-0.0107211)*c + (0.000385539)*c*c + (-6.28713e-06)*c*c*c + (3.75932e-08)*c*c*c*c ;
  }
  if(fRunIndex==28){
      v2c=-0.102173 + (0.0158262)*c + (-0.000706394)*c*c + (1.14263e-05)*c*c*c + (-6.14646e-08)*c*c*c*c ;
      v2s=0.02013 + (9.82197e-05)*c + (-1.71572e-05)*c*c + (-7.69646e-08)*c*c*c + (1.45113e-09)*c*c*c*c ;
  }
  if(fRunIndex==29){
      v2c=0.100539 + (-0.0115012)*c + (0.000295367)*c*c + (-2.57489e-06)*c*c*c + (5.41526e-09)*c*c*c*c ;
      v2s=0.138027 + (-0.00342397)*c + (-0.000392569)*c*c + (1.34626e-05)*c*c*c + (-1.06683e-07)*c*c*c*c ;
  }
  if(fRunIndex==30){
      v2c=0.0695326 + (-0.00752755)*c + (0.000237872)*c*c + (-3.31513e-06)*c*c*c + (1.69653e-08)*c*c*c*c ;
      v2s=0.0290133 + (-0.00073999)*c + (-3.35653e-05)*c*c + (7.61029e-07)*c*c*c + (-2.95792e-09)*c*c*c*c ;
  }
  if(fRunIndex==31){
      v2c=-0.103171 + (0.0134058)*c + (-0.000631602)*c*c + (1.10298e-05)*c*c*c + (-6.41177e-08)*c*c*c*c ;
      v2s=-0.137648 + (0.0111259)*c + (-0.000400035)*c*c + (5.90535e-06)*c*c*c + (-3.24676e-08)*c*c*c*c ;
  }
  if(fRunIndex==32){
      v2c=0.0522957 + (-0.00718492)*c + (0.000288586)*c*c + (-5.14534e-06)*c*c*c + (3.31303e-08)*c*c*c*c ;
      v2s=0.0731685 + (-0.00469449)*c + (9.57641e-05)*c*c + (-1.11891e-06)*c*c*c + (6.13984e-09)*c*c*c*c ;
  }
  if(fRunIndex==33){
      v2c=0.0882828 + (-0.00916109)*c + (0.00027821)*c*c + (-3.61312e-06)*c*c*c + (1.76418e-08)*c*c*c*c ;
      v2s=0.123302 + (-0.0156303)*c + (0.000623468)*c*c + (-9.83474e-06)*c*c*c + (5.24687e-08)*c*c*c*c ;
  }
  if(fRunIndex==34){
      v2c=0.0848636 + (-0.0142808)*c + (0.000723319)*c*c + (-1.32627e-05)*c*c*c + (7.77064e-08)*c*c*c*c ;
      v2s=-0.138483 + (0.0268553)*c + (-0.0013585)*c*c + (2.56034e-05)*c*c*c + (-1.63065e-07)*c*c*c*c ;
  }
  if(fRunIndex==35){
      v2c=-0.0134575 + (0.00267433)*c + (-0.000168784)*c*c + (3.24012e-06)*c*c*c + (-1.92766e-08)*c*c*c*c ;
      v2s=0.0596795 + (-0.00350808)*c + (2.97282e-05)*c*c + (3.52115e-07)*c*c*c + (-4.17462e-09)*c*c*c*c ;
  }
  if(fRunIndex==36){
      v2c=0.0784813 + (-0.00718389)*c + (0.000145303)*c*c + (-1.40121e-06)*c*c*c + (7.52775e-09)*c*c*c*c ;
      v2s=-0.0439181 + (0.0101363)*c + (-0.000494075)*c*c + (8.43356e-06)*c*c*c + (-4.68549e-08)*c*c*c*c ;
  }
  if(fRunIndex==37){
      v2c=-0.0118549 + (0.00328751)*c + (-0.000207156)*c*c + (3.92678e-06)*c*c*c + (-2.30059e-08)*c*c*c*c ;
      v2s=0.0344878 + (-0.000363669)*c + (-7.62481e-05)*c*c + (1.67135e-06)*c*c*c + (-9.64004e-09)*c*c*c*c ;
  }
  if(fRunIndex==38){
      v2c=-0.0254171 + (0.00297494)*c + (-0.000122361)*c*c + (1.78861e-06)*c*c*c + (-7.87609e-09)*c*c*c*c ;
      v2s=-0.0321399 + (0.00705636)*c + (-0.000325906)*c*c + (5.17669e-06)*c*c*c + (-2.86143e-08)*c*c*c*c ;
  }
  if(fRunIndex==39){
      v2c=-0.0368721 + (0.00253133)*c + (-0.000136793)*c*c + (2.51874e-06)*c*c*c + (-1.46055e-08)*c*c*c*c ;
      v2s=-0.0757763 + (0.005728)*c + (-0.000217489)*c*c + (3.11936e-06)*c*c*c + (-1.56919e-08)*c*c*c*c ;
  }
  if(fRunIndex==40){
      v2c=0.0166884 + (0.000509256)*c + (-0.000117911)*c*c + (2.65839e-06)*c*c*c + (-1.60875e-08)*c*c*c*c ;
      v2s=0.0336629 + (-0.00123023)*c + (-4.90123e-05)*c*c + (2.17862e-06)*c*c*c + (-1.98445e-08)*c*c*c*c ;
  }
  if(fRunIndex==41){
      v2c=-0.00887514 + (0.00491344)*c + (-0.000325782)*c*c + (6.54782e-06)*c*c*c + (-4.14143e-08)*c*c*c*c ;
      v2s=0.0720999 + (-0.00650524)*c + (0.000214752)*c*c + (-3.19681e-06)*c*c*c + (1.64573e-08)*c*c*c*c ;
  }
  if(fRunIndex==42){
      v2c=-0.0495903 + (0.00428229)*c + (-0.000181618)*c*c + (2.3913e-06)*c*c*c + (-8.55469e-09)*c*c*c*c ;
      v2s=0.0127982 + (-0.00690065)*c + (0.000342707)*c*c + (-6.40515e-06)*c*c*c + (3.83973e-08)*c*c*c*c ;
  }
  if(fRunIndex==43){
      v2c=0.0196783 + (-0.000221524)*c + (-7.09556e-05)*c*c + (1.76517e-06)*c*c*c + (-1.19612e-08)*c*c*c*c ;
      v2s=0.0350359 + (-0.0020458)*c + (5.18727e-05)*c*c + (-6.58401e-07)*c*c*c + (1.10936e-09)*c*c*c*c ;
  }
  if(fRunIndex==44){
      v2c=-0.0196041 + (0.00257682)*c + (-0.000174854)*c*c + (3.2168e-06)*c*c*c + (-1.85569e-08)*c*c*c*c ;
      v2s=-0.0870278 + (0.00634528)*c + (-0.000273588)*c*c + (4.77944e-06)*c*c*c + (-2.98422e-08)*c*c*c*c ;
  }
  if(fRunIndex==45){
      v2c=-0.0105827 + (0.00271714)*c + (-0.000152357)*c*c + (2.71899e-06)*c*c*c + (-1.5435e-08)*c*c*c*c ;
      v2s=0.135579 + (-0.0135952)*c + (0.000424184)*c*c + (-5.12316e-06)*c*c*c + (1.84497e-08)*c*c*c*c ;
  }
  if(fRunIndex==46){
      v2c=-0.017755 + (0.00411137)*c + (-0.000231724)*c*c + (4.17641e-06)*c*c*c + (-2.39164e-08)*c*c*c*c ;
      v2s=-0.0198628 + (0.0010032)*c + (-4.9708e-05)*c*c + (7.12455e-07)*c*c*c + (-4.22231e-09)*c*c*c*c ;
  }
  if(fRunIndex==47){
      v2c=-0.0111387 + (0.00258591)*c + (-0.000194287)*c*c + (3.84062e-06)*c*c*c + (-2.37192e-08)*c*c*c*c ;
      v2s=-0.0478616 + (0.00365349)*c + (-0.000152037)*c*c + (2.24323e-06)*c*c*c + (-1.23308e-08)*c*c*c*c ;
  }
  if(fRunIndex==48){
      v2c=0.070344 + (-0.0137332)*c + (0.000638536)*c*c + (-1.188e-05)*c*c*c + (7.4801e-08)*c*c*c*c ;
      v2s=-0.00208823 + (-0.00314157)*c + (9.9553e-05)*c*c + (-8.73891e-07)*c*c*c + (-1.72944e-09)*c*c*c*c ;
  }
  if(fRunIndex==51){
      v2c=0.017789 + (-0.00598224)*c + (0.000328477)*c*c + (-6.80787e-06)*c*c*c + (4.46905e-08)*c*c*c*c ;
      v2s=-0.0894734 + (0.007351)*c + (-0.000291085)*c*c + (4.38654e-06)*c*c*c + (-2.32261e-08)*c*c*c*c ;
  }
  if(fRunIndex==52){
      v2c=0.0555614 + (-0.00430576)*c + (0.00011519)*c*c + (-1.81584e-06)*c*c*c + (1.19214e-08)*c*c*c*c ;
      v2s=0.0901479 + (-0.0091505)*c + (0.000343757)*c*c + (-5.71204e-06)*c*c*c + (3.31217e-08)*c*c*c*c ;
  }
  if(fRunIndex==53){
      v2c=0.0181761 + (-0.00186378)*c + (4.60531e-05)*c*c + (-7.59581e-07)*c*c*c + (5.52465e-09)*c*c*c*c ;
      v2s=0.0456713 + (-0.00526813)*c + (0.000232978)*c*c + (-3.98724e-06)*c*c*c + (2.10384e-08)*c*c*c*c ;
  }
  if(fRunIndex==54){
      v2c=0.0456225 + (-0.00611173)*c + (0.000228414)*c*c + (-3.70527e-06)*c*c*c + (2.10774e-08)*c*c*c*c ;
      v2s=-0.00553874 + (0.00501799)*c + (-0.000331323)*c*c + (6.96978e-06)*c*c*c + (-4.72472e-08)*c*c*c*c ;
  }
  if(fRunIndex==55){
      v2c=-0.0578879 + (0.00664986)*c + (-0.000335443)*c*c + (6.0086e-06)*c*c*c + (-3.54293e-08)*c*c*c*c ;
      v2s=-0.0675269 + (0.00407581)*c + (-0.000148063)*c*c + (1.922e-06)*c*c*c + (-8.8635e-09)*c*c*c*c ;
  }
  if(fRunIndex==56){
      v2c=-0.000898903 + (0.00199222)*c + (-0.000207257)*c*c + (4.46148e-06)*c*c*c + (-2.71361e-08)*c*c*c*c ;
      v2s=-0.103794 + (0.0198667)*c + (-0.000879042)*c*c + (1.43177e-05)*c*c*c + (-8.09672e-08)*c*c*c*c ;
  }
  if(fRunIndex==57){
      v2c=0.00953087 + (0.000947327)*c + (-3.45822e-05)*c*c + (-4.22294e-07)*c*c*c + (9.08727e-09)*c*c*c*c ;
      v2s=0.0662314 + (-0.0076136)*c + (0.000288482)*c*c + (-4.82931e-06)*c*c*c + (2.67751e-08)*c*c*c*c ;
  }
  if(fRunIndex==58){
      v2c=0.00682728 + (-0.000293407)*c + (-2.88037e-05)*c*c + (6.63629e-07)*c*c*c + (-3.93121e-09)*c*c*c*c ;
      v2s=0.0085365 + (0.00138919)*c + (-8.15068e-05)*c*c + (1.26437e-06)*c*c*c + (-7.50322e-09)*c*c*c*c ;
  }
  if(fRunIndex==59){
      v2c=-0.0185226 + (0.00377881)*c + (-0.000226383)*c*c + (4.22897e-06)*c*c*c + (-2.51383e-08)*c*c*c*c ;
      v2s=-0.222605 + (0.0100699)*c + (-0.00029256)*c*c + (3.3855e-06)*c*c*c + (-1.33648e-08)*c*c*c*c ;
  }
  if(fRunIndex==60){
      v2c=-0.00111322 + (0.00291259)*c + (-0.000215769)*c*c + (4.13355e-06)*c*c*c + (-2.36883e-08)*c*c*c*c ;
      v2s=0.033686 + (-0.00124534)*c + (1.52616e-05)*c*c + (-2.99302e-07)*c*c*c + (1.05198e-09)*c*c*c*c ;
  }
  if(fRunIndex==61){
      v2c=0.00041873 + (-0.000231092)*c + (1.19226e-05)*c*c + (-6.29604e-07)*c*c*c + (7.26582e-09)*c*c*c*c ;
      v2s=0.00435438 + (0.0025265)*c + (-0.000153028)*c*c + (2.71169e-06)*c*c*c + (-1.67127e-08)*c*c*c*c ;
  }
  if(fRunIndex==62){
      v2c=-0.0573899 + (0.00972389)*c + (-0.000495772)*c*c + (9.0103e-06)*c*c*c + (-5.45756e-08)*c*c*c*c ;
      v2s=0.00702691 + (0.00698504)*c + (-0.000541929)*c*c + (1.19971e-05)*c*c*c + (-8.26942e-08)*c*c*c*c ;
  }
  if(fRunIndex==63){
      v2c=0.047222 + (-0.00604037)*c + (0.0002092)*c*c + (-3.07144e-06)*c*c*c + (1.68579e-08)*c*c*c*c ;
      v2s=0.0568815 + (-0.00607571)*c + (0.000244243)*c*c + (-4.22237e-06)*c*c*c + (2.38928e-08)*c*c*c*c ;
  }
  if(fRunIndex==64){
      v2c=0.0464686 + (-0.00463482)*c + (0.000146379)*c*c + (-2.23008e-06)*c*c*c + (1.29549e-08)*c*c*c*c ;
      v2s=-0.00571995 + (0.00273485)*c + (-0.000131945)*c*c + (2.08986e-06)*c*c*c + (-1.26717e-08)*c*c*c*c ;
  }
  if(fRunIndex==65){
      v2c=0.0544794 + (-0.00616918)*c + (0.000173808)*c*c + (-2.01338e-06)*c*c*c + (8.78123e-09)*c*c*c*c ;
      v2s=0.0557002 + (-0.00701287)*c + (0.000340418)*c*c + (-6.12754e-06)*c*c*c + (3.3384e-08)*c*c*c*c ;
  }
  if(fRunIndex==66){
      v2c=0.0111176 + (0.000262254)*c + (-7.20349e-05)*c*c + (1.38755e-06)*c*c*c + (-6.65016e-09)*c*c*c*c ;
      v2s=0.0644184 + (-0.00518515)*c + (0.000174227)*c*c + (-2.7347e-06)*c*c*c + (1.29547e-08)*c*c*c*c ;
  }
  if(fRunIndex==67){
      v2c=-0.0481037 + (0.00823686)*c + (-0.000418053)*c*c + (7.40186e-06)*c*c*c + (-4.19609e-08)*c*c*c*c ;
      v2s=0.0560171 + (-0.002556)*c + (2.54529e-05)*c*c + (-5.22269e-09)*c*c*c + (-2.19983e-09)*c*c*c*c ;
  }
  if(fRunIndex==68){
      v2c=-0.0139809 + (0.00274305)*c + (-0.00017744)*c*c + (3.20286e-06)*c*c*c + (-1.74137e-08)*c*c*c*c ;
      v2s=-0.218673 + (0.00764806)*c + (-0.000222739)*c*c + (3.14656e-06)*c*c*c + (-1.73548e-08)*c*c*c*c ;
  }
  if(fRunIndex==69){
      v2c=0.0736938 + (-0.0122004)*c + (0.000498227)*c*c + (-8.08779e-06)*c*c*c + (4.52749e-08)*c*c*c*c ;
      v2s=-0.236254 + (0.00979816)*c + (-0.000280698)*c*c + (3.4103e-06)*c*c*c + (-1.47299e-08)*c*c*c*c ;
  }
  if(fRunIndex==70){
      v2c=-0.0366545 + (0.00399604)*c + (-0.000196876)*c*c + (3.49686e-06)*c*c*c + (-2.01142e-08)*c*c*c*c ;
      v2s=0.0220911 + (-0.000529459)*c + (9.10408e-06)*c*c + (-3.57334e-07)*c*c*c + (1.87909e-09)*c*c*c*c ;
  }
  if(fRunIndex==71){
      v2c=-0.039698 + (0.00654615)*c + (-0.000323865)*c*c + (5.50998e-06)*c*c*c + (-3.00866e-08)*c*c*c*c ;
      v2s=0.0650078 + (-0.00502636)*c + (0.000128459)*c*c + (-1.3877e-06)*c*c*c + (3.21917e-09)*c*c*c*c ;
  }
  if(fRunIndex==72){
      v2c=0.0594227 + (-0.00872513)*c + (0.000368684)*c*c + (-6.34986e-06)*c*c*c + (3.79523e-08)*c*c*c*c ;
      v2s=0.050348 + (-0.00512087)*c + (0.000228643)*c*c + (-4.40797e-06)*c*c*c + (2.62482e-08)*c*c*c*c ;
  }
  if(fRunIndex==73){
      v2c=0.131281 + (-0.0162697)*c + (0.000584063)*c*c + (-8.64349e-06)*c*c*c + (4.61009e-08)*c*c*c*c ;
      v2s=0.0271721 + (-0.000487039)*c + (1.60624e-05)*c*c + (-5.38029e-07)*c*c*c + (2.08589e-09)*c*c*c*c ;
  }
  if(fRunIndex==74){
      v2c=-0.0128651 + (0.00946426)*c + (-0.00055547)*c*c + (9.37245e-06)*c*c*c + (-4.78224e-08)*c*c*c*c ;
      v2s=0.0903476 + (-0.0164972)*c + (0.000927481)*c*c + (-1.88459e-05)*c*c*c + (1.20249e-07)*c*c*c*c ;
  }
  if(fRunIndex==75){
      v2c=0.0204147 + (-0.00139139)*c + (-2.47732e-05)*c*c + (1.1214e-06)*c*c*c + (-7.22533e-09)*c*c*c*c ;
      v2s=0.0214247 + (-0.000147929)*c + (-7.58241e-06)*c*c + (-7.05565e-08)*c*c*c + (-4.86528e-10)*c*c*c*c ;
  }
  if(fRunIndex==76){
      v2c=0.0469095 + (-0.00541793)*c + (0.000133238)*c*c + (-1.14334e-06)*c*c*c + (3.4826e-09)*c*c*c*c ;
      v2s=0.0304763 + (-0.00213957)*c + (8.27588e-05)*c*c + (-1.34495e-06)*c*c*c + (4.81134e-09)*c*c*c*c ;
  }
  if(fRunIndex==77){
      v2c=0.113822 + (-0.0154286)*c + (0.000574533)*c*c + (-8.67874e-06)*c*c*c + (4.76727e-08)*c*c*c*c ;
      v2s=0.010623 + (0.00513467)*c + (-0.000372733)*c*c + (7.40692e-06)*c*c*c + (-4.64093e-08)*c*c*c*c ;
  }
  if(fRunIndex==78){
      v2c=0.020658 + (-0.00274725)*c + (5.8872e-05)*c*c + (-3.43999e-07)*c*c*c + (-1.04273e-10)*c*c*c*c ;
      v2s=0.00219379 + (0.000972033)*c + (-1.41816e-05)*c*c + (-4.41149e-07)*c*c*c + (3.81037e-09)*c*c*c*c ;
  }
  if(fRunIndex==79){
      v2c=0.227324 + (-0.0264915)*c + (0.000915298)*c*c + (-1.3343e-05)*c*c*c + (7.1044e-08)*c*c*c*c ;
      v2s=0.0474072 + (-0.00939515)*c + (0.000730108)*c*c + (-1.80811e-05)*c*c*c + (1.3003e-07)*c*c*c*c ;
  }
  if(fRunIndex==81){
      v2c=-0.0368028 + (0.00368927)*c + (-0.000228321)*c*c + (4.48174e-06)*c*c*c + (-2.75298e-08)*c*c*c*c ;
      v2s=-0.238807 + (0.0108667)*c + (-0.000350649)*c*c + (4.95568e-06)*c*c*c + (-2.53022e-08)*c*c*c*c ;
  }
  if(fRunIndex==82){
      v2c=-0.0214072 + (0.00119258)*c + (-7.80215e-05)*c*c + (1.37682e-06)*c*c*c + (-7.31103e-09)*c*c*c*c ;
      v2s=-0.233825 + (0.00735109)*c + (-0.000142737)*c*c + (9.15013e-07)*c*c*c + (-1.98463e-10)*c*c*c*c ;
  }
  if(fRunIndex==83){
      v2c=-0.0716708 + (0.00855289)*c + (-0.00037249)*c*c + (6.41638e-06)*c*c*c + (-3.66318e-08)*c*c*c*c ;
      v2s=-0.329702 + (0.033231)*c + (-0.00131012)*c*c + (2.1101e-05)*c*c*c + (-1.18905e-07)*c*c*c*c ;
  }
  if(fRunIndex==84){
      v2c=-0.0259022 + (0.00711773)*c + (-0.000434778)*c*c + (8.94774e-06)*c*c*c + (-5.84408e-08)*c*c*c*c ;
      v2s=-0.0936287 + (-0.00455282)*c + (0.000390982)*c*c + (-7.95224e-06)*c*c*c + (4.95612e-08)*c*c*c*c ;
  }
  if(fRunIndex==85){
      v2c=0.0263223 + (-0.00475593)*c + (9.63236e-05)*c*c + (-4.95311e-07)*c*c*c + (-5.47428e-10)*c*c*c*c ;
      v2s=-0.220593 + (0.00823892)*c + (-0.000248477)*c*c + (3.28082e-06)*c*c*c + (-1.45098e-08)*c*c*c*c ;
  }
  if(fRunIndex==86){
      v2c=-0.0157962 + (0.00408922)*c + (-0.000282087)*c*c + (5.81961e-06)*c*c*c + (-3.6481e-08)*c*c*c*c ;
      v2s=0.0340979 + (-0.00139471)*c + (6.48561e-06)*c*c + (3.02184e-08)*c*c*c + (-1.05879e-09)*c*c*c*c ;
  }
  if(fRunIndex==87){
      v2c=0.0577366 + (-0.00471394)*c + (0.000107329)*c*c + (-7.28927e-07)*c*c*c + (-1.51199e-09)*c*c*c*c ;
      v2s=0.0350171 + (-0.00402687)*c + (0.00016046)*c*c + (-2.39625e-06)*c*c*c + (9.98146e-09)*c*c*c*c ;
  }
  if(fRunIndex==88){
      v2c=-0.00547971 + (0.00161905)*c + (-0.00013831)*c*c + (2.94662e-06)*c*c*c + (-1.92504e-08)*c*c*c*c ;
      v2s=-0.0722251 + (0.0127634)*c + (-0.000612284)*c*c + (1.07553e-05)*c*c*c + (-6.34054e-08)*c*c*c*c ;
  }
  if(fRunIndex==89){
      v2c=0.165574 + (-0.0241913)*c + (0.00101221)*c*c + (-1.62279e-05)*c*c*c + (8.77252e-08)*c*c*c*c ;
      v2s=0.0692954 + (-0.0147219)*c + (0.00068782)*c*c + (-1.11091e-05)*c*c*c + (5.37295e-08)*c*c*c*c ;
  }
  if(fRunIndex==91){
      v2c=-0.0334666 + (0.00583866)*c + (-0.000329883)*c*c + (6.08932e-06)*c*c*c + (-3.59853e-08)*c*c*c*c ;
      v2s=-0.229775 + (0.010061)*c + (-0.00028144)*c*c + (3.3231e-06)*c*c*c + (-1.40799e-08)*c*c*c*c ;
  }
  if(fRunIndex==95){
      v2c=-0.0993203 + (0.0150547)*c + (-0.000713531)*c*c + (1.22733e-05)*c*c*c + (-6.98669e-08)*c*c*c*c ;
      v2s=-0.217746 + (0.0068375)*c + (-0.000171962)*c*c + (2.05506e-06)*c*c*c + (-8.78703e-09)*c*c*c*c ;
  }
  if(fRunIndex==96){
      v2c=-0.0398801 + (0.00406028)*c + (-0.000214066)*c*c + (3.69056e-06)*c*c*c + (-1.99339e-08)*c*c*c*c ;
      v2s=-0.240195 + (0.00782806)*c + (-0.000131365)*c*c + (3.18193e-07)*c*c*c + (4.88009e-09)*c*c*c*c ;
  }
  if(fRunIndex==97){
      v2c=-0.0325538 + (0.00250986)*c + (-0.000136762)*c*c + (2.35572e-06)*c*c*c + (-1.23413e-08)*c*c*c*c ;
      v2s=-0.240415 + (0.010435)*c + (-0.000291848)*c*c + (3.35357e-06)*c*c*c + (-1.35833e-08)*c*c*c*c ;
  }
  if(fRunIndex==98){
      v2c=-0.00692338 + (0.00379335)*c + (-0.000332833)*c*c + (7.35052e-06)*c*c*c + (-4.89956e-08)*c*c*c*c ;
      v2s=-0.263042 + (0.0134291)*c + (-0.000465334)*c*c + (7.09595e-06)*c*c*c + (-3.90332e-08)*c*c*c*c ;
  }
  if(fRunIndex==101){
      v2c=-0.049409 + (0.0073433)*c + (-0.000376933)*c*c + (7.39945e-06)*c*c*c + (-4.967e-08)*c*c*c*c ;
      v2s=0.00564511 + (-0.00101958)*c + (0.000121077)*c*c + (-3.44889e-06)*c*c*c + (2.42733e-08)*c*c*c*c ;
  }
  if(fRunIndex==102){
      v2c=-0.0393185 + (0.00610331)*c + (-0.000286099)*c*c + (4.61693e-06)*c*c*c + (-2.44104e-08)*c*c*c*c ;
      v2s=-0.00898987 + (0.00463766)*c + (-0.000269122)*c*c + (5.46456e-06)*c*c*c + (-3.85677e-08)*c*c*c*c ;
  }
  if(fRunIndex==103){
      v2c=-0.0171815 + (0.00271616)*c + (-9.49063e-05)*c*c + (5.68888e-07)*c*c*c + (2.77777e-09)*c*c*c*c ;
      v2s=-0.0348362 + (0.00661647)*c + (-0.000196625)*c*c + (1.60914e-06)*c*c*c + (-4.34769e-09)*c*c*c*c ;
  }
  if(fRunIndex==104){
      v2c=0.0453893 + (-0.00220981)*c + (-6.34642e-05)*c*c + (2.5268e-06)*c*c*c + (-1.81319e-08)*c*c*c*c ;
      v2s=0.0131555 + (0.00115036)*c + (-5.70975e-05)*c*c + (5.20224e-07)*c*c*c + (-3.11332e-09)*c*c*c*c ;
  }
  if(fRunIndex==105){
      v2c=0.0052505 + (0.000251094)*c + (-3.95143e-05)*c*c + (5.6354e-07)*c*c*c + (-1.86603e-09)*c*c*c*c ;
      v2s=0.079882 + (-0.00841369)*c + (0.000342924)*c*c + (-5.64554e-06)*c*c*c + (2.89511e-08)*c*c*c*c ;
  }
  if(fRunIndex==106){
      v2c=0.0561693 + (-0.00693389)*c + (0.000243013)*c*c + (-3.7816e-06)*c*c*c + (2.21666e-08)*c*c*c*c ;
      v2s=0.115847 + (-0.0145063)*c + (0.000597908)*c*c + (-9.37418e-06)*c*c*c + (4.57406e-08)*c*c*c*c ;
  }
  if(fRunIndex==107){
      v2c=-0.160892 + (0.0259075)*c + (-0.00120084)*c*c + (1.99814e-05)*c*c*c + (-1.09708e-07)*c*c*c*c ;
      v2s=0.0922001 + (-0.00861465)*c + (0.000196846)*c*c + (-5.34939e-07)*c*c*c + (-1.20442e-08)*c*c*c*c ;
  }
  if(fRunIndex==108){
      v2c=0.0325605 + (-0.00795036)*c + (0.000436987)*c*c + (-9.37794e-06)*c*c*c + (6.59912e-08)*c*c*c*c ;
      v2s=0.0894763 + (-0.00694006)*c + (0.000214928)*c*c + (-3.05745e-06)*c*c*c + (1.33679e-08)*c*c*c*c ;
  }
  if(fRunIndex==109){
      v2c=-8.30505e-05 + (0.00304264)*c + (-0.000228365)*c*c + (4.43511e-06)*c*c*c + (-2.5788e-08)*c*c*c*c ;
      v2s=0.0276356 + (-0.000486726)*c + (-1.12152e-05)*c*c + (1.56258e-07)*c*c*c + (-2.17826e-09)*c*c*c*c ;
  }
  if(fRunIndex==112){
      v2c=0.0273919 + (-0.00243837)*c + (2.41667e-05)*c*c + (3.06367e-07)*c*c*c + (-3.9028e-09)*c*c*c*c ;
      v2s=0.0119428 + (0.00161447)*c + (-0.00010526)*c*c + (1.76239e-06)*c*c*c + (-1.12844e-08)*c*c*c*c ;
  }
  if(fRunIndex==114){
      v2c=0.114287 + (-0.016125)*c + (0.000673002)*c*c + (-1.11635e-05)*c*c*c + (6.39533e-08)*c*c*c*c ;
      v2s=0.0353709 + (-0.00169691)*c + (-3.93183e-05)*c*c + (2.08815e-06)*c*c*c + (-2.03159e-08)*c*c*c*c ;
  }
  if(fRunIndex==115){
      v2c=-0.0109276 + (0.00378322)*c + (-0.000252651)*c*c + (4.599e-06)*c*c*c + (-2.36568e-08)*c*c*c*c ;
      v2s=0.0185448 + (0.00102852)*c + (-6.04966e-05)*c*c + (5.89933e-07)*c*c*c + (-2.5921e-09)*c*c*c*c ;
  }
  if(fRunIndex==116){
      v2c=-0.117563 + (0.00405449)*c + (-0.000112806)*c*c + (4.01903e-06)*c*c*c + (-4.25358e-08)*c*c*c*c ;
      v2s=0.0280797 + (-0.00528592)*c + (0.000375784)*c*c + (-9.96825e-06)*c*c*c + (7.72147e-08)*c*c*c*c ;
  }
  if(fRunIndex==117){
      v2c=0.0301869 + (-0.000751512)*c + (-9.43748e-05)*c*c + (2.60697e-06)*c*c*c + (-1.73948e-08)*c*c*c*c ;
      v2s=0.0594442 + (-0.00499442)*c + (0.000170617)*c*c + (-2.80903e-06)*c*c*c + (1.51075e-08)*c*c*c*c ;
  }
  if(fRunIndex==118){
      v2c=0.0118942 + (-0.000669595)*c + (-3.6634e-05)*c*c + (1.40794e-06)*c*c*c + (-1.15253e-08)*c*c*c*c ;
      v2s=-0.0174395 + (0.00460355)*c + (-0.00023192)*c*c + (3.83941e-06)*c*c*c + (-2.30761e-08)*c*c*c*c ;
  }
  if(fRunIndex==119){
      v2c=0.018229 + (-0.00144698)*c + (-1.07295e-06)*c*c + (6.22528e-07)*c*c*c + (-4.90504e-09)*c*c*c*c ;
      v2s=0.0122761 + (0.000985351)*c + (-6.02936e-05)*c*c + (9.41867e-07)*c*c*c + (-6.59756e-09)*c*c*c*c ;
  }
  if(fRunIndex==120){
      v2c=-0.0968065 + (0.0131861)*c + (-0.000585224)*c*c + (9.82175e-06)*c*c*c + (-5.51313e-08)*c*c*c*c ;
      v2s=0.0719334 + (-0.00752912)*c + (0.000362464)*c*c + (-6.65114e-06)*c*c*c + (3.70295e-08)*c*c*c*c ;
  }
  if(fRunIndex==121){
      v2c=0.0296329 + (-0.00448004)*c + (0.00010239)*c*c + (-2.40647e-07)*c*c*c + (-5.4555e-09)*c*c*c*c ;
      v2s=0.0774925 + (-0.00635491)*c + (0.0001939)*c*c + (-3.04755e-06)*c*c*c + (1.66055e-08)*c*c*c*c ;
  }
  if(fRunIndex==122){
      v2c=0.00956939 + (0.000460793)*c + (-0.000105912)*c*c + (2.38286e-06)*c*c*c + (-1.446e-08)*c*c*c*c ;
      v2s=0.0152677 + (0.00345263)*c + (-0.000262399)*c*c + (5.5895e-06)*c*c*c + (-3.84974e-08)*c*c*c*c ;
  }
  if(fRunIndex==123){
      v2c=0.0495618 + (-0.00059293)*c + (-0.00025227)*c*c + (6.61088e-06)*c*c*c + (-4.24063e-08)*c*c*c*c ;
      v2s=0.0862227 + (-0.0190161)*c + (0.0010647)*c*c + (-1.94957e-05)*c*c*c + (1.08665e-07)*c*c*c*c ;
  }
  if(fRunIndex==124){
      v2c=-0.0408078 + (0.00565108)*c + (-0.00027528)*c*c + (5.2165e-06)*c*c*c + (-3.3244e-08)*c*c*c*c ;
      v2s=0.00169862 + (0.00371982)*c + (-0.000169757)*c*c + (2.07347e-06)*c*c*c + (-8.30142e-09)*c*c*c*c ;
  }
  if(fRunIndex==125){
      v2c=0.0285065 + (-0.00348583)*c + (0.000106189)*c*c + (-1.56726e-06)*c*c*c + (9.40552e-09)*c*c*c*c ;
      v2s=0.0363645 + (-0.00252206)*c + (0.000122553)*c*c + (-2.74066e-06)*c*c*c + (1.66856e-08)*c*c*c*c ;
  }
  if(fRunIndex==126){
      v2c=-0.00864424 + (0.00157555)*c + (-0.000127367)*c*c + (2.47447e-06)*c*c*c + (-1.42302e-08)*c*c*c*c ;
      v2s=-0.213003 + (0.00643156)*c + (-0.000121219)*c*c + (6.06838e-07)*c*c*c + (1.48381e-09)*c*c*c*c ;
  }
  if(fRunIndex==127){
      v2c=-0.0336884 + (0.00299236)*c + (-0.000147557)*c*c + (2.47908e-06)*c*c*c + (-1.33495e-08)*c*c*c*c ;
      v2s=-0.271925 + (0.0158755)*c + (-0.000583802)*c*c + (9.14376e-06)*c*c*c + (-5.05302e-08)*c*c*c*c ;
  }
  if(fRunIndex==128){
      v2c=-0.0639711 + (0.00822881)*c + (-0.000428416)*c*c + (7.98075e-06)*c*c*c + (-4.81995e-08)*c*c*c*c ;
      v2s=-0.22634 + (0.00908632)*c + (-0.00029)*c*c + (4.14095e-06)*c*c*c + (-2.15823e-08)*c*c*c*c ;
  }
  if(fRunIndex==129){
      v2c=0.0652753 + (-0.00529737)*c + (0.000133883)*c*c + (-1.63339e-06)*c*c*c + (7.88294e-09)*c*c*c*c ;
      v2s=0.00357003 + (0.00503303)*c + (-0.000316531)*c*c + (6.45793e-06)*c*c*c + (-4.44467e-08)*c*c*c*c ;
  }
  if(fRunIndex==131){
      v2c=0.117255 + (-0.0149009)*c + (0.000560178)*c*c + (-8.75754e-06)*c*c*c + (4.85024e-08)*c*c*c*c ;
      v2s=-0.00180786 + (0.00224258)*c + (-9.84056e-05)*c*c + (1.39525e-06)*c*c*c + (-7.63512e-09)*c*c*c*c ;
  }
  if(fRunIndex==132){
      v2c=0.0100287 + (-0.00051338)*c + (-2.7702e-05)*c*c + (5.0205e-07)*c*c*c + (-1.12682e-09)*c*c*c*c ;
      v2s=0.00556895 + (0.00211082)*c + (-0.000125448)*c*c + (2.30409e-06)*c*c*c + (-1.63703e-08)*c*c*c*c ;
  }
  if(fRunIndex==133){
      v2c=-0.00755913 + (0.00327793)*c + (-0.00024545)*c*c + (5.03417e-06)*c*c*c + (-3.11617e-08)*c*c*c*c ;
      v2s=0.0398477 + (-0.000728372)*c + (-2.63935e-05)*c*c + (7.85276e-07)*c*c*c + (-7.72703e-09)*c*c*c*c ;
  }
  if(fRunIndex==134){
      v2c=0.00344056 + (0.000171769)*c + (-8.00887e-05)*c*c + (1.76696e-06)*c*c*c + (-1.08652e-08)*c*c*c*c ;
      v2s=-0.200474 + (0.00736532)*c + (-0.000221809)*c*c + (2.91763e-06)*c*c*c + (-1.40827e-08)*c*c*c*c ;
  }
  if(fRunIndex==135){
      v2c=-0.0251554 + (0.00227996)*c + (-0.000133757)*c*c + (2.50604e-06)*c*c*c + (-1.48528e-08)*c*c*c*c ;
      v2s=-0.229223 + (0.0091122)*c + (-0.000245035)*c*c + (2.64826e-06)*c*c*c + (-9.57268e-09)*c*c*c*c ;
  }
  if(fRunIndex==136){
      v2c=-0.0308547 + (0.00209237)*c + (-0.000183034)*c*c + (4.54164e-06)*c*c*c + (-3.26557e-08)*c*c*c*c ;
      v2s=-0.241702 + (0.00988355)*c + (-0.000263014)*c*c + (2.73376e-06)*c*c*c + (-9.79971e-09)*c*c*c*c ;
  }

/*  
  
  if(fRunIndex>=1 && fRunIndex<=45){ // 137161 -137848 //period 1
//    v2c = 1.133676e-02-1.759454e-04*x-6.217033e-05*x*x+1.532507e-06*x*x*x-9.665265e-09*x*x*x*x;
//    v2s = 2.173690e-02-8.410340e-04*x-1.391841e-05*x*x+5.003301e-07*x*x*x-4.175974e-09*x*x*x*x;
    v2c=2.295070e-02-1.050171e-03*x-3.432843e-05*x*x+1.156549e-06*x*x*x-7.781377e-09*x*x*x*x ;
    v2s=3.994886e-02-1.961957e-03*x+1.641958e-05*x*x+1.234025e-07*x*x*x-2.175698e-09*x*x*x*x ;
  }
  if(fRunIndex>=46 && fRunIndex<=58){ //138125-138275//period 2
    v2c =-1.394612e-03+1.333726e-03*x-1.115323e-04*x*x+2.103762e-06*x*x*x-1.186760e-08*x*x*x*x ;
    v2s =-7.788201e-03+7.178984e-04*x-4.553597e-05*x*x+6.941817e-07*x*x*x-4.465249e-09*x*x*x*x ;
  }
  if(fRunIndex>=59 && fRunIndex<=78){ //138359-138730//period 3
    v2c = 1.448809e-02+-7.894029e-04*x-2.693871e-05*x*x+7.944017e-07*x*x*x-4.272135e-09*x*x*x*x ;
    v2s = 3.442216e-03+-5.574560e-05*x-1.758343e-05*x*x+2.386547e-07*x*x*x-2.404210e-09*x*x*x*x ;
  }
  if(fRunIndex>=79){//period 4
    v2c =-1.553536e-03+7.761106e-04*x-9.529723e-05*x*x+2.026750e-06*x*x*x-1.217009e-08*x*x*x*x ;
    v2s =-7.787317e-02+3.843455e-03*x-1.221144e-04*x*x+1.478448e-06*x*x*x-7.442702e-09*x*x*x*x ;
  }
*/  
  
  return phi+v2c*TMath::Sin(2.*phi)-v2s*TMath::Cos(2.*phi) ;

}
//____________________________________________________________________________
Double_t  AliAnalysisTaskPi0v2::ApplyFlatteningV0C(Double_t phi, Double_t c){
  //LHC10h

    Double_t v2c =0.;
    Double_t v2s =0;

    if(fRunIndex==1){
	v2c=-0.038901 + (0.00539521)*c + (-0.000213271)*c*c + (3.01147e-06)*c*c*c + (-1.35977e-08)*c*c*c*c ;
	v2s=0.0693294 + (-0.00851319)*c + (0.000337408)*c*c + (-5.33097e-06)*c*c*c + (2.90883e-08)*c*c*c*c ;
    }
    if(fRunIndex==2){
	v2c=-0.0943695 + (0.0137829)*c + (-0.000569334)*c*c + (8.95522e-06)*c*c*c + (-4.78396e-08)*c*c*c*c ;
	v2s=0.00971525 + (-0.00181026)*c + (6.15798e-05)*c*c + (-8.51046e-07)*c*c*c + (4.57403e-09)*c*c*c*c ;
    }
    if(fRunIndex==4){
	v2c=-0.00585475 + (0.000307793)*c + (1.89109e-06)*c*c + (-3.14804e-07)*c*c*c + (3.54341e-09)*c*c*c*c ;
	v2s=0.00374568 + (-0.00123435)*c + (8.01941e-05)*c*c + (-1.70807e-06)*c*c*c + (1.14576e-08)*c*c*c*c ;
    }
    if(fRunIndex==5){
	v2c=-0.0301267 + (0.00206698)*c + (-4.82282e-05)*c*c + (4.39355e-07)*c*c*c + (-1.3747e-09)*c*c*c*c ;
	v2s=-0.0441947 + (0.00501774)*c + (-0.000181736)*c*c + (2.69235e-06)*c*c*c + (-1.37411e-08)*c*c*c*c ;
    }
    if(fRunIndex==6){
	v2c=0.190652 + (-0.0254912)*c + (0.00105707)*c*c + (-1.63872e-05)*c*c*c + (8.51113e-08)*c*c*c*c ;
	v2s=0.0803124 + (-0.0119724)*c + (0.000575997)*c*c + (-1.07892e-05)*c*c*c + (6.67214e-08)*c*c*c*c ;
    }
    if(fRunIndex==7){
	v2c=0.00875044 + (-0.00343147)*c + (0.0002198)*c*c + (-4.76292e-06)*c*c*c + (3.24844e-08)*c*c*c*c ;
	v2s=-0.081288 + (0.0146429)*c + (-0.000728327)*c*c + (1.32462e-05)*c*c*c + (-7.95362e-08)*c*c*c*c ;
    }
    if(fRunIndex==8){
	v2c=0.0381367 + (-0.0082048)*c + (0.000410914)*c*c + (-7.28851e-06)*c*c*c + (4.22008e-08)*c*c*c*c ;
	v2s=-0.023338 + (0.003381)*c + (-0.000167791)*c*c + (3.13998e-06)*c*c*c + (-1.97233e-08)*c*c*c*c ;
    }
    if(fRunIndex==9){
	v2c=-0.00734261 + (0.00121272)*c + (-5.31181e-05)*c*c + (7.88974e-07)*c*c*c + (-3.86356e-09)*c*c*c*c ;
	v2s=-0.0228333 + (0.00266052)*c + (-0.000112824)*c*c + (1.89515e-06)*c*c*c + (-1.07446e-08)*c*c*c*c ;
    }
    if(fRunIndex==10){
	v2c=-0.0456385 + (0.0102956)*c + (-0.000529689)*c*c + (9.62809e-06)*c*c*c + (-5.7301e-08)*c*c*c*c ;
	v2s=0.0142197 + (-0.0031677)*c + (0.000190389)*c*c + (-4.10911e-06)*c*c*c + (2.81791e-08)*c*c*c*c ;
    }
    if(fRunIndex==11){
	v2c=-0.00999249 + (0.00123939)*c + (-4.23895e-05)*c*c + (6.55415e-07)*c*c*c + (-3.98256e-09)*c*c*c*c ;
	v2s=-0.00809817 + (0.00199727)*c + (-0.000116718)*c*c + (2.33449e-06)*c*c*c + (-1.46407e-08)*c*c*c*c ;
    }
    if(fRunIndex==12){
	v2c=0.0183566 + (-0.0020805)*c + (6.79726e-05)*c*c + (-7.96708e-07)*c*c*c + (3.1214e-09)*c*c*c*c ;
	v2s=-0.00714715 + (0.000721418)*c + (-3.50355e-05)*c*c + (9.09996e-07)*c*c*c + (-7.92198e-09)*c*c*c*c ;
    }
    if(fRunIndex==13){
	v2c=0.0056915 + (0.000167211)*c + (-4.68018e-05)*c*c + (1.21891e-06)*c*c*c + (-8.70124e-09)*c*c*c*c ;
	v2s=0.00986223 + (-0.00247957)*c + (0.000133539)*c*c + (-2.51362e-06)*c*c*c + (1.5397e-08)*c*c*c*c ;
    }
    if(fRunIndex==14){
	v2c=0.0479623 + (-0.00948196)*c + (0.000551795)*c*c + (-1.05692e-05)*c*c*c + (6.38754e-08)*c*c*c*c ;
	v2s=-0.00438456 + (-0.00152976)*c + (2.00645e-05)*c*c + (9.91017e-07)*c*c*c + (-1.35776e-08)*c*c*c*c ;
    }
    if(fRunIndex==15){
	v2c=-0.00384088 + (0.00121719)*c + (-6.88665e-05)*c*c + (1.26138e-06)*c*c*c + (-7.13987e-09)*c*c*c*c ;
	v2s=0.0273017 + (-0.00314818)*c + (7.85031e-05)*c*c + (-9.07336e-08)*c*c*c + (-6.42421e-09)*c*c*c*c ;
    }
    if(fRunIndex==16){
	v2c=-0.0313175 + (0.00368025)*c + (-0.000126895)*c*c + (1.7443e-06)*c*c*c + (-9.12442e-09)*c*c*c*c ;
	v2s=-0.0812178 + (0.0117402)*c + (-0.000506287)*c*c + (8.36766e-06)*c*c*c + (-4.69866e-08)*c*c*c*c ;
    }
    if(fRunIndex==17){
	v2c=0.062893 + (-0.0128204)*c + (0.000591931)*c*c + (-1.00835e-05)*c*c*c + (5.72782e-08)*c*c*c*c ;
	v2s=-0.143105 + (0.0183921)*c + (-0.000787289)*c*c + (1.33237e-05)*c*c*c + (-7.71707e-08)*c*c*c*c ;
    }
    if(fRunIndex==18){
	v2c=-0.0029353 + (-0.0016097)*c + (0.000199618)*c*c + (-5.42008e-06)*c*c*c + (4.08897e-08)*c*c*c*c ;
	v2s=0.0482015 + (-0.00115233)*c + (3.96019e-05)*c*c + (-7.7193e-07)*c*c*c + (3.90004e-09)*c*c*c*c ;
    }
    if(fRunIndex==19){
	v2c=-0.0397728 + (0.0129154)*c + (-0.000703903)*c*c + (1.30972e-05)*c*c*c + (-7.93547e-08)*c*c*c*c ;
	v2s=0.181945 + (-0.0289807)*c + (0.00133685)*c*c + (-2.30866e-05)*c*c*c + (1.32746e-07)*c*c*c*c ;
    }
    if(fRunIndex==20){
	v2c=-0.00927858 + (0.0021033)*c + (-0.000113971)*c*c + (2.17508e-06)*c*c*c + (-1.36812e-08)*c*c*c*c ;
	v2s=0.0201645 + (-0.00331311)*c + (0.000164916)*c*c + (-3.08376e-06)*c*c*c + (1.92565e-08)*c*c*c*c ;
    }
    if(fRunIndex==21){
	v2c=-0.0254895 + (0.00388685)*c + (-0.000171429)*c*c + (2.87209e-06)*c*c*c + (-1.60625e-08)*c*c*c*c ;
	v2s=0.00237609 + (0.000664255)*c + (-4.93436e-05)*c*c + (8.86415e-07)*c*c*c + (-4.51913e-09)*c*c*c*c ;
    }
    if(fRunIndex==22){
	v2c=0.0308671 + (-0.00365824)*c + (0.000148483)*c*c + (-2.52946e-06)*c*c*c + (1.49859e-08)*c*c*c*c ;
	v2s=-0.00805023 + (0.00154452)*c + (-6.89678e-05)*c*c + (1.14317e-06)*c*c*c + (-5.98991e-09)*c*c*c*c ;
    }
    if(fRunIndex==23){
	v2c=0.122018 + (-0.0233973)*c + (0.00116385)*c*c + (-2.11232e-05)*c*c*c + (1.25836e-07)*c*c*c*c ;
	v2s=-0.0511879 + (0.00708815)*c + (-0.000227282)*c*c + (2.27542e-06)*c*c*c + (-4.83128e-09)*c*c*c*c ;
    }
    if(fRunIndex==24){
	v2c=0.00888562 + (-0.000839847)*c + (2.52839e-05)*c*c + (-2.28921e-07)*c*c*c + (1.01057e-11)*c*c*c*c ;
	v2s=0.00814357 + (-0.000788997)*c + (2.43948e-05)*c*c + (-2.99558e-07)*c*c*c + (1.21675e-09)*c*c*c*c ;
    }
    if(fRunIndex==25){
	v2c=0.00360833 + (0.000379977)*c + (-3.68193e-05)*c*c + (5.86861e-07)*c*c*c + (-2.01826e-09)*c*c*c*c ;
	v2s=-0.00642047 + (0.00060076)*c + (-2.42907e-05)*c*c + (4.12933e-07)*c*c*c + (-2.3557e-09)*c*c*c*c ;
    }
    if(fRunIndex==26){
	v2c=-0.0131971 + (0.00189761)*c + (-7.7551e-05)*c*c + (1.25884e-06)*c*c*c + (-7.38989e-09)*c*c*c*c ;
	v2s=0.00679891 + (-0.000614692)*c + (1.04559e-06)*c*c + (4.34463e-07)*c*c*c + (-4.75821e-09)*c*c*c*c ;
    }
    if(fRunIndex==27){
	v2c=-0.00887089 + (0.00126865)*c + (-7.15766e-05)*c*c + (1.43811e-06)*c*c*c + (-9.408e-09)*c*c*c*c ;
	v2s=-0.0302105 + (0.00349565)*c + (-0.000126108)*c*c + (1.835e-06)*c*c*c + (-9.27845e-09)*c*c*c*c ;
    }
    if(fRunIndex==28){
	v2c=-0.0187616 + (0.0021974)*c + (-5.11412e-05)*c*c + (2.94628e-07)*c*c*c + (7.4194e-10)*c*c*c*c ;
	v2s=-0.0119204 + (-0.00156879)*c + (0.000132746)*c*c + (-2.67371e-06)*c*c*c + (1.61097e-08)*c*c*c*c ;
    }
    if(fRunIndex==29){
	v2c=-0.121018 + (0.0228924)*c + (-0.00127093)*c*c + (2.48724e-05)*c*c*c + (-1.56121e-07)*c*c*c*c ;
	v2s=-0.0687197 + (0.015702)*c + (-0.000764074)*c*c + (1.38991e-05)*c*c*c + (-8.48812e-08)*c*c*c*c ;
    }
    if(fRunIndex==30){
	v2c=0.01582 + (-0.0012293)*c + (1.62768e-05)*c*c + (2.36172e-07)*c*c*c + (-3.67153e-09)*c*c*c*c ;
	v2s=-0.00875651 + (0.00135376)*c + (-5.16763e-05)*c*c + (6.17614e-07)*c*c*c + (-1.64099e-09)*c*c*c*c ;
    }
    if(fRunIndex==31){
	v2c=-0.044721 + (0.00756113)*c + (-0.000365759)*c*c + (6.15804e-06)*c*c*c + (-3.40075e-08)*c*c*c*c ;
	v2s=-0.0500828 + (-0.000695112)*c + (0.000130909)*c*c + (-3.38476e-06)*c*c*c + (2.41179e-08)*c*c*c*c ;
    }
    if(fRunIndex==32){
	v2c=0.0306236 + (-0.00443307)*c + (0.000196517)*c*c + (-3.48059e-06)*c*c*c + (2.10266e-08)*c*c*c*c ;
	v2s=0.00819214 + (-0.000560558)*c + (-5.03783e-06)*c*c + (4.37862e-07)*c*c*c + (-3.93834e-09)*c*c*c*c ;
    }
    if(fRunIndex==33){
	v2c=0.0120077 + (-0.00139976)*c + (5.85677e-05)*c*c + (-9.06209e-07)*c*c*c + (4.61403e-09)*c*c*c*c ;
	v2s=-0.0197677 + (0.00253589)*c + (-0.000101374)*c*c + (1.59347e-06)*c*c*c + (-8.65961e-09)*c*c*c*c ;
    }
    if(fRunIndex==34){
	v2c=0.0759546 + (-0.0131538)*c + (0.000577001)*c*c + (-9.39636e-06)*c*c*c + (5.16814e-08)*c*c*c*c ;
	v2s=-0.0727864 + (0.00956655)*c + (-0.000406191)*c*c + (7.50453e-06)*c*c*c + (-4.93819e-08)*c*c*c*c ;
    }
    if(fRunIndex==35){
	v2c=-0.00622564 + (0.00107647)*c + (-6.25644e-05)*c*c + (1.28553e-06)*c*c*c + (-8.46286e-09)*c*c*c*c ;
	v2s=0.0051704 + (-0.00165348)*c + (8.69531e-05)*c*c + (-1.69895e-06)*c*c*c + (1.09762e-08)*c*c*c*c ;
    }
    if(fRunIndex==36){
	v2c=0.0117872 + (-0.000807391)*c + (-1.70999e-05)*c*c + (1.19089e-06)*c*c*c + (-1.13278e-08)*c*c*c*c ;
	v2s=-0.000162495 + (1.04465e-05)*c + (-1.07156e-05)*c*c + (2.54308e-07)*c*c*c + (-7.04948e-10)*c*c*c*c ;
    }
    if(fRunIndex==37){
	v2c=-0.0119236 + (0.00290105)*c + (-0.000162086)*c*c + (3.22645e-06)*c*c*c + (-2.07075e-08)*c*c*c*c ;
	v2s=-0.0178535 + (0.00283024)*c + (-0.00012812)*c*c + (2.02404e-06)*c*c*c + (-1.0248e-08)*c*c*c*c ;
    }
    if(fRunIndex==38){
	v2c=0.023002 + (-0.00433112)*c + (0.000229879)*c*c + (-4.52792e-06)*c*c*c + (2.92715e-08)*c*c*c*c ;
	v2s=-0.0217933 + (0.00235291)*c + (-0.000108054)*c*c + (2.08937e-06)*c*c*c + (-1.3524e-08)*c*c*c*c ;
    }
    if(fRunIndex==39){
	v2c=0.00795585 + (-0.000372683)*c + (1.85243e-05)*c*c + (-4.28133e-07)*c*c*c + (3.02452e-09)*c*c*c*c ;
	v2s=-0.0849826 + (0.00529466)*c + (-0.000163148)*c*c + (2.1186e-06)*c*c*c + (-9.48208e-09)*c*c*c*c ;
    }
    if(fRunIndex==40){
	v2c=0.00653019 + (-0.00127651)*c + (7.59264e-05)*c*c + (-1.6363e-06)*c*c*c + (1.11375e-08)*c*c*c*c ;
	v2s=-0.0107686 + (0.000982587)*c + (-2.97061e-05)*c*c + (4.04795e-07)*c*c*c + (-1.91289e-09)*c*c*c*c ;
    }
    if(fRunIndex==41){
	v2c=0.012007 + (-0.00102703)*c + (3.48714e-05)*c*c + (-5.13011e-07)*c*c*c + (2.57827e-09)*c*c*c*c ;
	v2s=-0.00713083 + (0.00101952)*c + (-4.7675e-05)*c*c + (9.12723e-07)*c*c*c + (-5.99694e-09)*c*c*c*c ;
    }
    if(fRunIndex==42){
	v2c=-0.0486696 + (0.00723345)*c + (-0.000324232)*c*c + (5.25317e-06)*c*c*c + (-2.78981e-08)*c*c*c*c ;
	v2s=-0.0485877 + (-0.00064506)*c + (9.15858e-05)*c*c + (-1.97486e-06)*c*c*c + (1.28831e-08)*c*c*c*c ;
    }
    if(fRunIndex==43){
	v2c=0.0003541 + (-0.000236832)*c + (1.82448e-05)*c*c + (-4.20764e-07)*c*c*c + (2.96882e-09)*c*c*c*c ;
	v2s=0.00107616 + (-0.000462058)*c + (2.46761e-05)*c*c + (-4.33199e-07)*c*c*c + (2.47324e-09)*c*c*c*c ;
    }
    if(fRunIndex==44){
	v2c=-0.0778259 + (0.0128797)*c + (-0.000514371)*c*c + (7.70095e-06)*c*c*c + (-3.89026e-08)*c*c*c*c ;
	v2s=-0.106839 + (0.00670369)*c + (-0.000184721)*c*c + (1.64581e-06)*c*c*c + (-1.47049e-09)*c*c*c*c ;
    }
    if(fRunIndex==45){
	v2c=0.0513049 + (-0.00715356)*c + (0.000289106)*c*c + (-4.46938e-06)*c*c*c + (2.36599e-08)*c*c*c*c ;
	v2s=-0.0252464 + (0.00170124)*c + (-3.04376e-05)*c*c + (9.4935e-08)*c*c*c + (9.41807e-10)*c*c*c*c ;
    }
    if(fRunIndex==46){
	v2c=0.0317697 + (-0.00121243)*c + (5.94543e-05)*c*c + (-1.26341e-06)*c*c*c + (8.49329e-09)*c*c*c*c ;
	v2s=-0.092276 + (0.00299941)*c + (-5.24713e-05)*c*c + (3.07331e-07)*c*c*c + (7.44139e-10)*c*c*c*c ;
    }
    if(fRunIndex==47){
	v2c=0.0429796 + (-0.000861588)*c + (6.82123e-07)*c*c + (1.96388e-07)*c*c*c + (-1.56984e-09)*c*c*c*c ;
	v2s=-0.133047 + (0.00853153)*c + (-0.000287983)*c*c + (4.1155e-06)*c*c*c + (-2.00239e-08)*c*c*c*c ;
    }
    if(fRunIndex==48){
	v2c=0.0480979 + (-0.000354969)*c + (-6.45176e-05)*c*c + (2.34807e-06)*c*c*c + (-1.9069e-08)*c*c*c*c ;
	v2s=-0.138315 + (0.0126874)*c + (-0.000672052)*c*c + (1.34449e-05)*c*c*c + (-8.64756e-08)*c*c*c*c ;
    }
    if(fRunIndex==51){
	v2c=-0.0114422 + (0.00506157)*c + (-0.000227041)*c*c + (3.97457e-06)*c*c*c + (-2.35322e-08)*c*c*c*c ;
	v2s=-0.208135 + (0.0160488)*c + (-0.000536132)*c*c + (7.2714e-06)*c*c*c + (-3.28328e-08)*c*c*c*c ;
    }
    if(fRunIndex==52){
	v2c=0.00113321 + (0.000922798)*c + (-6.5971e-05)*c*c + (1.40416e-06)*c*c*c + (-9.29905e-09)*c*c*c*c ;
	v2s=0.0550392 + (-0.00725661)*c + (0.000298984)*c*c + (-4.85473e-06)*c*c*c + (2.72063e-08)*c*c*c*c ;
    }
    if(fRunIndex==53){
	v2c=-0.003508 + (6.54451e-06)*c + (1.28603e-05)*c*c + (-3.83674e-07)*c*c*c + (3.1053e-09)*c*c*c*c ;
	v2s=0.00714978 + (-0.00198147)*c + (0.000113909)*c*c + (-2.1949e-06)*c*c*c + (1.34933e-08)*c*c*c*c ;
    }
    if(fRunIndex==54){
	v2c=-0.0024598 + (5.87979e-05)*c + (9.36516e-06)*c*c + (-2.66045e-07)*c*c*c + (1.80016e-09)*c*c*c*c ;
	v2s=-0.00651135 + (0.00109883)*c + (-4.41501e-05)*c*c + (6.19398e-07)*c*c*c + (-2.61977e-09)*c*c*c*c ;
    }
    if(fRunIndex==55){
	v2c=0.0472936 + (-0.00156718)*c + (3.03713e-05)*c*c + (-1.97139e-07)*c*c*c + (-3.55205e-10)*c*c*c*c ;
	v2s=-0.155014 + (0.00717506)*c + (-0.000168137)*c*c + (1.55419e-06)*c*c*c + (-3.64579e-09)*c*c*c*c ;
    }
    if(fRunIndex==56){
	v2c=0.011328 + (0.000924574)*c + (-6.803e-05)*c*c + (9.76956e-07)*c*c*c + (-2.90494e-09)*c*c*c*c ;
	v2s=-0.0412661 + (0.0107218)*c + (-0.000533108)*c*c + (9.75563e-06)*c*c*c + (-5.95608e-08)*c*c*c*c ;
    }
    if(fRunIndex==57){
	v2c=0.00900452 + (0.000588343)*c + (-3.02208e-05)*c*c + (2.13239e-07)*c*c*c + (9.25727e-10)*c*c*c*c ;
	v2s=0.13086 + (-0.0204768)*c + (0.000928956)*c*c + (-1.63917e-05)*c*c*c + (9.8273e-08)*c*c*c*c ;
    }
    if(fRunIndex==58){
	v2c=-0.0143951 + (0.00247373)*c + (-0.000107909)*c*c + (1.69603e-06)*c*c*c + (-9.05555e-09)*c*c*c*c ;
	v2s=-0.0209397 + (0.00180446)*c + (-5.97034e-05)*c*c + (8.17388e-07)*c*c*c + (-3.77989e-09)*c*c*c*c ;
    }
    if(fRunIndex==59){
	v2c=-0.0619279 + (0.00393575)*c + (-0.00011867)*c*c + (1.45001e-06)*c*c*c + (-5.98792e-09)*c*c*c*c ;
	v2s=-0.248103 + (0.0158553)*c + (-0.000529631)*c*c + (7.46158e-06)*c*c*c + (-3.60091e-08)*c*c*c*c ;
    }
    if(fRunIndex==60){
	v2c=-0.00886115 + (0.00152741)*c + (-6.94062e-05)*c*c + (1.11727e-06)*c*c*c + (-5.81351e-09)*c*c*c*c ;
	v2s=-0.000911417 + (0.000122344)*c + (-1.24605e-06)*c*c + (-6.29249e-08)*c*c*c + (8.81971e-10)*c*c*c*c ;
    }
    if(fRunIndex==61){
	v2c=-0.00510166 + (0.00040596)*c + (-4.37341e-06)*c*c + (-1.15302e-07)*c*c*c + (1.53428e-09)*c*c*c*c ;
	v2s=-0.0170123 + (0.00220818)*c + (-9.50905e-05)*c*c + (1.59566e-06)*c*c*c + (-8.92053e-09)*c*c*c*c ;
    }
    if(fRunIndex==62){
	v2c=-0.00332885 + (-0.000430498)*c + (5.60341e-05)*c*c + (-1.47922e-06)*c*c*c + (1.09723e-08)*c*c*c*c ;
	v2s=-0.0519534 + (0.00892818)*c + (-0.000434243)*c*c + (7.8407e-06)*c*c*c + (-4.69363e-08)*c*c*c*c ;
    }
    if(fRunIndex==63){
	v2c=-0.00176529 + (0.00123097)*c + (-6.59677e-05)*c*c + (1.24591e-06)*c*c*c + (-7.90285e-09)*c*c*c*c ;
	v2s=-0.0181792 + (0.00222004)*c + (-7.96902e-05)*c*c + (1.13047e-06)*c*c*c + (-5.72272e-09)*c*c*c*c ;
    }
    if(fRunIndex==64){
	v2c=0.00546184 + (-0.000449463)*c + (1.50833e-05)*c*c + (-1.72501e-07)*c*c*c + (3.58774e-10)*c*c*c*c ;
	v2s=0.013284 + (-0.00260671)*c + (0.000134403)*c*c + (-2.461e-06)*c*c*c + (1.47694e-08)*c*c*c*c ;
    }
    if(fRunIndex==65){
	v2c=0.0108447 + (-0.0017699)*c + (4.62716e-05)*c*c + (-3.4647e-07)*c*c*c + (5.24536e-10)*c*c*c*c ;
	v2s=0.0302247 + (-0.00605136)*c + (0.000290676)*c*c + (-4.87706e-06)*c*c*c + (2.65104e-08)*c*c*c*c ;
    }
    if(fRunIndex==66){
	v2c=-0.0100092 + (0.00144457)*c + (-6.512e-05)*c*c + (1.05356e-06)*c*c*c + (-5.52231e-09)*c*c*c*c ;
	v2s=-0.00380228 + (0.000610003)*c + (-3.47943e-05)*c*c + (7.21188e-07)*c*c*c + (-4.76727e-09)*c*c*c*c ;
    }
    if(fRunIndex==67){
	v2c=-0.0229901 + (0.00400737)*c + (-0.000206024)*c*c + (3.93038e-06)*c*c*c + (-2.46388e-08)*c*c*c*c ;
	v2s=-0.00286041 + (0.000648607)*c + (-2.37397e-05)*c*c + (2.0372e-07)*c*c*c + (2.97679e-10)*c*c*c*c ;
    }
    if(fRunIndex==68){
	v2c=-0.0603647 + (0.0045685)*c + (-0.000171555)*c*c + (2.66377e-06)*c*c*c + (-1.45892e-08)*c*c*c*c ;
	v2s=-0.21447 + (0.00909524)*c + (-0.000265298)*c*c + (3.73262e-06)*c*c*c + (-1.84484e-08)*c*c*c*c ;
    }
    if(fRunIndex==69){
	v2c=-0.0207879 + (0.000139613)*c + (1.27802e-05)*c*c + (-6.40328e-07)*c*c*c + (6.35949e-09)*c*c*c*c ;
	v2s=-0.266459 + (0.0156802)*c + (-0.000481647)*c*c + (6.4597e-06)*c*c*c + (-3.047e-08)*c*c*c*c ;
    }
    if(fRunIndex==70){
	v2c=-0.0324843 + (0.00376852)*c + (-0.000136873)*c*c + (2.01389e-06)*c*c*c + (-1.03877e-08)*c*c*c*c ;
	v2s=-0.0178378 + (0.00145491)*c + (-3.55421e-05)*c*c + (2.25993e-07)*c*c*c + (7.06258e-10)*c*c*c*c ;
    }
    if(fRunIndex==71){
	v2c=-0.0138645 + (0.00167315)*c + (-6.68544e-05)*c*c + (1.07133e-06)*c*c*c + (-5.92659e-09)*c*c*c*c ;
	v2s=0.0359716 + (-0.00461076)*c + (0.000192868)*c*c + (-3.21409e-06)*c*c*c + (1.84468e-08)*c*c*c*c ;
    }
    if(fRunIndex==72){
	v2c=0.0111712 + (-0.00175546)*c + (6.87862e-05)*c*c + (-8.99029e-07)*c*c*c + (3.44492e-09)*c*c*c*c ;
	v2s=-0.00352823 + (0.000480673)*c + (-2.66091e-05)*c*c + (5.54922e-07)*c*c*c + (-3.84534e-09)*c*c*c*c ;
    }
    if(fRunIndex==73){
	v2c=-0.000939145 + (-0.000235945)*c + (2.51908e-05)*c*c + (-6.52435e-07)*c*c*c + (4.8218e-09)*c*c*c*c ;
	v2s=0.00666852 + (-0.000373105)*c + (-1.22838e-05)*c*c + (7.3621e-07)*c*c*c + (-7.2867e-09)*c*c*c*c ;
    }
    if(fRunIndex==74){
	v2c=0.0474235 + (-0.00890967)*c + (0.000427741)*c*c + (-7.38723e-06)*c*c*c + (4.17242e-08)*c*c*c*c ;
	v2s=0.0142517 + (-0.00326338)*c + (0.000142729)*c*c + (-2.27831e-06)*c*c*c + (1.14084e-08)*c*c*c*c ;
    }
    if(fRunIndex==75){
	v2c=0.0237808 + (-0.00310491)*c + (0.000134288)*c*c + (-2.27429e-06)*c*c*c + (1.32749e-08)*c*c*c*c ;
	v2s=-0.0258905 + (0.00417719)*c + (-0.000193901)*c*c + (3.41383e-06)*c*c*c + (-2.01695e-08)*c*c*c*c ;
    }
    if(fRunIndex==76){
	v2c=0.0149865 + (-0.00140739)*c + (4.24537e-05)*c*c + (-4.06273e-07)*c*c*c + (5.44496e-10)*c*c*c*c ;
	v2s=0.00820134 + (-0.000519432)*c + (6.82211e-06)*c*c + (1.06291e-07)*c*c*c + (-1.65344e-09)*c*c*c*c ;
    }
    if(fRunIndex==77){
	v2c=0.048123 + (-0.00608197)*c + (0.000246409)*c*c + (-4.05028e-06)*c*c*c + (2.33536e-08)*c*c*c*c ;
	v2s=-0.0214729 + (0.0029389)*c + (-0.000133826)*c*c + (2.29845e-06)*c*c*c + (-1.29253e-08)*c*c*c*c ;
    }
    if(fRunIndex==78){
	v2c=-0.0172924 + (0.00202216)*c + (-6.478e-05)*c*c + (8.08739e-07)*c*c*c + (-3.73211e-09)*c*c*c*c ;
	v2s=-0.0355638 + (0.00121518)*c + (3.97839e-05)*c*c + (-1.62346e-06)*c*c*c + (1.27105e-08)*c*c*c*c ;
    }
    if(fRunIndex==79){
	v2c=0.040762 + (-0.00311921)*c + (6.66868e-05)*c*c + (-1.18554e-06)*c*c*c + (1.07519e-08)*c*c*c*c ;
	v2s=-0.00558438 + (-0.00154821)*c + (0.000165732)*c*c + (-4.56377e-06)*c*c*c + (3.56988e-08)*c*c*c*c ;
    }
    if(fRunIndex==81){
	v2c=-0.0247995 + (0.00179444)*c + (-0.000102085)*c*c + (2.31093e-06)*c*c*c + (-1.61924e-08)*c*c*c*c ;
	v2s=-0.249009 + (0.0133297)*c + (-0.000415717)*c*c + (5.84494e-06)*c*c*c + (-2.88179e-08)*c*c*c*c ;
    }
    if(fRunIndex==82){
	v2c=-0.0249921 + (0.00199208)*c + (-7.86121e-05)*c*c + (1.40547e-06)*c*c*c + (-9.39651e-09)*c*c*c*c ;
	v2s=-0.254469 + (0.0150694)*c + (-0.000499848)*c*c + (7.08318e-06)*c*c*c + (-3.40806e-08)*c*c*c*c ;
    }
    if(fRunIndex==83){
	v2c=-0.0721 + (0.00951559)*c + (-0.000431762)*c*c + (8.4632e-06)*c*c*c + (-5.79485e-08)*c*c*c*c ;
	v2s=-0.00864625 + (-0.00706196)*c + (0.000313773)*c*c + (-5.10066e-06)*c*c*c + (2.9662e-08)*c*c*c*c ;
    }
    if(fRunIndex==84){
	v2c=-0.0198112 + (-0.0062788)*c + (0.000353589)*c*c + (-5.70525e-06)*c*c*c + (2.82617e-08)*c*c*c*c ;
	v2s=-0.0704647 + (0.0031015)*c + (-0.000244985)*c*c + (6.17258e-06)*c*c*c + (-4.43331e-08)*c*c*c*c ;
    }
    if(fRunIndex==85){
	v2c=-0.0161955 + (-0.0017361)*c + (0.000130612)*c*c + (-2.9548e-06)*c*c*c + (2.05684e-08)*c*c*c*c ;
	v2s=-0.204119 + (0.00714989)*c + (-0.000162549)*c*c + (1.87459e-06)*c*c*c + (-7.79143e-09)*c*c*c*c ;
    }
    if(fRunIndex==86){
	v2c=0.00995155 + (-0.00196822)*c + (0.00010617)*c*c + (-2.01071e-06)*c*c*c + (1.22728e-08)*c*c*c*c ;
	v2s=0.0148901 + (-0.002581)*c + (0.000132002)*c*c + (-2.44601e-06)*c*c*c + (1.48144e-08)*c*c*c*c ;
    }
    if(fRunIndex==87){
	v2c=-0.00411031 + (0.000510639)*c + (-2.27814e-05)*c*c + (6.41192e-07)*c*c*c + (-5.85467e-09)*c*c*c*c ;
	v2s=0.0121118 + (-0.00281)*c + (0.000143754)*c*c + (-2.62557e-06)*c*c*c + (1.61207e-08)*c*c*c*c ;
    }
    if(fRunIndex==88){
	v2c=-0.00570308 + (-0.000945116)*c + (7.42606e-05)*c*c + (-1.51212e-06)*c*c*c + (9.57993e-09)*c*c*c*c ;
	v2s=-0.0282666 + (0.00398475)*c + (-0.000165678)*c*c + (2.68542e-06)*c*c*c + (-1.4688e-08)*c*c*c*c ;
    }
    if(fRunIndex==89){
	v2c=0.0845785 + (-0.012611)*c + (0.000561197)*c*c + (-9.39333e-06)*c*c*c + (5.16139e-08)*c*c*c*c ;
	v2s=0.149689 + (-0.0193648)*c + (0.000771629)*c*c + (-1.21901e-05)*c*c*c + (6.6443e-08)*c*c*c*c ;
    }
    if(fRunIndex==91){
	v2c=-0.0289605 + (0.00127922)*c + (-4.47445e-06)*c*c + (-4.00293e-07)*c*c*c + (3.82966e-09)*c*c*c*c ;
	v2s=-0.241098 + (0.00906132)*c + (-0.000159552)*c*c + (6.47243e-07)*c*c*c + (4.61048e-09)*c*c*c*c ;
    }
    if(fRunIndex==95){
	v2c=0.0144072 + (-0.000966056)*c + (6.47254e-05)*c*c + (-1.84164e-06)*c*c*c + (1.51588e-08)*c*c*c*c ;
	v2s=-0.306678 + (0.0184746)*c + (-0.000540218)*c*c + (6.36681e-06)*c*c*c + (-2.37512e-08)*c*c*c*c ;
    }
    if(fRunIndex==96){
	v2c=0.00506182 + (-0.000632672)*c + (1.89752e-05)*c*c + (-2.14999e-07)*c*c*c + (3.12179e-10)*c*c*c*c ;
	v2s=-0.218693 + (0.00676786)*c + (-0.000115485)*c*c + (5.56737e-07)*c*c*c + (2.68891e-09)*c*c*c*c ;
    }
    if(fRunIndex==97){
	v2c=-0.0225062 + (0.00204445)*c + (-6.76028e-05)*c*c + (8.92976e-07)*c*c*c + (-4.00576e-09)*c*c*c*c ;
	v2s=-0.260286 + (0.0139325)*c + (-0.000425631)*c*c + (5.64038e-06)*c*c*c + (-2.51554e-08)*c*c*c*c ;
    }
    if(fRunIndex==98){
	v2c=-0.0261849 + (0.0063862)*c + (-0.000293076)*c*c + (4.37234e-06)*c*c*c + (-2.04899e-08)*c*c*c*c ;
	v2s=-0.17744 + (0.00236157)*c + (9.2782e-06)*c*c + (-5.63126e-07)*c*c*c + (4.9714e-09)*c*c*c*c ;
    }
    if(fRunIndex==101){
	v2c=0.101445 + (-0.015664)*c + (0.000747255)*c*c + (-1.36979e-05)*c*c*c + (8.36903e-08)*c*c*c*c ;
	v2s=0.0152975 + (-0.00329946)*c + (0.00018627)*c*c + (-3.88572e-06)*c*c*c + (2.64895e-08)*c*c*c*c ;
    }
    if(fRunIndex==102){
	v2c=-0.0344927 + (0.00470964)*c + (-0.000188941)*c*c + (3.0206e-06)*c*c*c + (-1.69735e-08)*c*c*c*c ;
	v2s=0.0142624 + (-0.00240623)*c + (9.16657e-05)*c*c + (-1.14663e-06)*c*c*c + (4.14195e-09)*c*c*c*c ;
    }
    if(fRunIndex==103){
	v2c=-0.0238316 + (0.00322675)*c + (-0.000123462)*c*c + (1.76428e-06)*c*c*c + (-8.66176e-09)*c*c*c*c ;
	v2s=0.0227184 + (-0.00299954)*c + (0.00014422)*c*c + (-2.81018e-06)*c*c*c + (1.82546e-08)*c*c*c*c ;
    }
    if(fRunIndex==104){
	v2c=0.0118014 + (-0.000998833)*c + (3.77271e-05)*c*c + (-6.19183e-07)*c*c*c + (3.49065e-09)*c*c*c*c ;
	v2s=0.000977264 + (0.00012542)*c + (1.40219e-08)*c*c + (-2.24472e-09)*c*c*c + (-3.39763e-10)*c*c*c*c ;
    }
    if(fRunIndex==105){
	v2c=0.0312332 + (-0.00348485)*c + (0.000108654)*c*c + (-1.28443e-06)*c*c*c + (4.80188e-09)*c*c*c*c ;
	v2s=-0.0344234 + (0.00443313)*c + (-0.000173309)*c*c + (2.7124e-06)*c*c*c + (-1.42518e-08)*c*c*c*c ;
    }
    if(fRunIndex==106){
	v2c=0.0369386 + (-0.00434063)*c + (0.00018894)*c*c + (-3.45688e-06)*c*c*c + (2.16617e-08)*c*c*c*c ;
	v2s=-0.00675251 + (0.000691284)*c + (-1.08033e-06)*c*c + (-5.57352e-07)*c*c*c + (6.28148e-09)*c*c*c*c ;
    }
    if(fRunIndex==107){
	v2c=-0.0441483 + (0.00803486)*c + (-0.00034871)*c*c + (4.76349e-06)*c*c*c + (-1.81393e-08)*c*c*c*c ;
	v2s=-0.141046 + (0.0165183)*c + (-0.000592546)*c*c + (8.956e-06)*c*c*c + (-4.86662e-08)*c*c*c*c ;
    }
    if(fRunIndex==108){
	v2c=-0.0531169 + (0.00637348)*c + (-0.000218334)*c*c + (2.66873e-06)*c*c*c + (-9.86859e-09)*c*c*c*c ;
	v2s=-0.00809622 + (-0.000552159)*c + (6.27984e-05)*c*c + (-1.33136e-06)*c*c*c + (7.45069e-09)*c*c*c*c ;
    }
    if(fRunIndex==109){
	v2c=-0.00189593 + (0.00153994)*c + (-9.92105e-05)*c*c + (1.98756e-06)*c*c*c + (-1.24677e-08)*c*c*c*c ;
	v2s=0.0252107 + (-0.00330651)*c + (0.000128362)*c*c + (-1.89951e-06)*c*c*c + (9.45477e-09)*c*c*c*c ;
    }
    if(fRunIndex==112){
	v2c=-0.00939028 + (0.00128226)*c + (-5.24125e-05)*c*c + (8.50382e-07)*c*c*c + (-4.65785e-09)*c*c*c*c ;
	v2s=0.00201633 + (-0.00116151)*c + (7.59727e-05)*c*c + (-1.59962e-06)*c*c*c + (1.04645e-08)*c*c*c*c ;
    }
    if(fRunIndex==114){
	v2c=-0.00783351 + (0.00122277)*c + (-3.37736e-05)*c*c + (3.12136e-07)*c*c*c + (-9.19359e-10)*c*c*c*c ;
	v2s=0.0103833 + (-0.00215968)*c + (9.96942e-05)*c*c + (-1.71974e-06)*c*c*c + (1.01511e-08)*c*c*c*c ;
    }
    if(fRunIndex==115){
	v2c=0.00588164 + (-0.00083683)*c + (5.13684e-05)*c*c + (-1.14177e-06)*c*c*c + (8.11899e-09)*c*c*c*c ;
	v2s=0.0037561 + (-0.000786163)*c + (3.88084e-05)*c*c + (-7.73678e-07)*c*c*c + (5.27396e-09)*c*c*c*c ;
    }
    if(fRunIndex==116){
	v2c=0.0610685 + (-0.00392239)*c + (6.23272e-05)*c*c + (9.84282e-07)*c*c*c + (-1.75418e-08)*c*c*c*c ;
	v2s=-0.106295 + (0.0124333)*c + (-0.000510299)*c*c + (8.51235e-06)*c*c*c + (-4.9787e-08)*c*c*c*c ;
    }
    if(fRunIndex==117){
	v2c=-0.0132264 + (0.0024143)*c + (-0.000106307)*c*c + (1.74632e-06)*c*c*c + (-9.79752e-09)*c*c*c*c ;
	v2s=-0.00788943 + (0.000443907)*c + (-1.83984e-05)*c*c + (3.81018e-07)*c*c*c + (-2.52604e-09)*c*c*c*c ;
    }
    if(fRunIndex==118){
	v2c=-0.045835 + (0.00849686)*c + (-0.000408675)*c*c + (7.3031e-06)*c*c*c + (-4.34655e-08)*c*c*c*c ;
	v2s=-0.0455199 + (0.00718851)*c + (-0.000340502)*c*c + (6.00457e-06)*c*c*c + (-3.51883e-08)*c*c*c*c ;
    }
    if(fRunIndex==119){
	v2c=-0.00627236 + (0.00146712)*c + (-9.05375e-05)*c*c + (1.96962e-06)*c*c*c + (-1.35867e-08)*c*c*c*c ;
	v2s=-0.00544096 + (0.000598445)*c + (-7.61032e-06)*c*c + (-2.09579e-07)*c*c*c + (3.00091e-09)*c*c*c*c ;
    }
    if(fRunIndex==120){
	v2c=-0.0376823 + (0.00376965)*c + (-8.01573e-05)*c*c + (-2.80792e-09)*c*c*c + (6.79508e-09)*c*c*c*c ;
	v2s=-0.0675803 + (0.00849852)*c + (-0.000298266)*c*c + (3.93178e-06)*c*c*c + (-1.68582e-08)*c*c*c*c ;
    }
    if(fRunIndex==121){
	v2c=0.0285153 + (-0.00340253)*c + (0.000130854)*c*c + (-2.02449e-06)*c*c*c + (1.09594e-08)*c*c*c*c ;
	v2s=0.0224991 + (-0.00524188)*c + (0.000269226)*c*c + (-4.93926e-06)*c*c*c + (2.97604e-08)*c*c*c*c ;
    }
    if(fRunIndex==122){
	v2c=0.00379321 + (0.000532446)*c + (-6.16649e-05)*c*c + (1.67786e-06)*c*c*c + (-1.29277e-08)*c*c*c*c ;
	v2s=-0.011158 + (0.00256124)*c + (-0.000150921)*c*c + (3.22488e-06)*c*c*c + (-2.16809e-08)*c*c*c*c ;
    }
    if(fRunIndex==123){
	v2c=-0.0607006 + (0.00831906)*c + (-0.000398705)*c*c + (6.66385e-06)*c*c*c + (-3.5222e-08)*c*c*c*c ;
	v2s=0.0389204 + (-0.00798361)*c + (0.000297135)*c*c + (-3.79751e-06)*c*c*c + (1.47914e-08)*c*c*c*c ;
    }
    if(fRunIndex==124){
	v2c=-0.0288952 + (0.00325728)*c + (-0.000104126)*c*c + (1.29062e-06)*c*c*c + (-6.03539e-09)*c*c*c*c ;
	v2s=0.0318884 + (-0.00301585)*c + (0.000131913)*c*c + (-2.5702e-06)*c*c*c + (1.71313e-08)*c*c*c*c ;
    }
    if(fRunIndex==125){
	v2c=-0.00413859 + (0.00024613)*c + (5.33902e-06)*c*c + (-3.86963e-07)*c*c*c + (3.89277e-09)*c*c*c*c ;
	v2s=-0.0492298 + (0.00624793)*c + (-0.000243844)*c*c + (3.69001e-06)*c*c*c + (-1.9156e-08)*c*c*c*c ;
    }
    if(fRunIndex==126){
	v2c=-0.018816 + (0.00134076)*c + (-4.50486e-05)*c*c + (5.84259e-07)*c*c*c + (-2.34235e-09)*c*c*c*c ;
	v2s=-0.295912 + (0.0159608)*c + (-0.000497569)*c*c + (6.92442e-06)*c*c*c + (-3.35159e-08)*c*c*c*c ;
    }
    if(fRunIndex==127){
	v2c=-0.039222 + (0.00400933)*c + (-0.0001341)*c*c + (1.65214e-06)*c*c*c + (-6.33829e-09)*c*c*c*c ;
	v2s=-0.290884 + (0.0147027)*c + (-0.000429021)*c*c + (5.64116e-06)*c*c*c + (-2.58784e-08)*c*c*c*c ;
    }
    if(fRunIndex==128){
	v2c=0.0168685 + (-0.00285775)*c + (0.000141099)*c*c + (-2.48846e-06)*c*c*c + (1.41845e-08)*c*c*c*c ;
	v2s=-0.278776 + (0.0113771)*c + (-0.000279052)*c*c + (3.16144e-06)*c*c*c + (-1.18545e-08)*c*c*c*c ;
    }
    if(fRunIndex==129){
	v2c=0.0490322 + (-0.00736827)*c + (0.000304529)*c*c + (-4.89472e-06)*c*c*c + (2.60746e-08)*c*c*c*c ;
	v2s=-0.0961198 + (0.0179264)*c + (-0.000914589)*c*c + (1.70145e-05)*c*c*c + (-1.03423e-07)*c*c*c*c ;
    }
    if(fRunIndex==131){
	v2c=0.075924 + (-0.0108106)*c + (0.000493471)*c*c + (-8.39074e-06)*c*c*c + (4.67109e-08)*c*c*c*c ;
	v2s=0.0700446 + (-0.0105027)*c + (0.00048267)*c*c + (-8.4422e-06)*c*c*c + (4.91292e-08)*c*c*c*c ;
    }
    if(fRunIndex==132){
	v2c=0.00752887 + (-0.00142672)*c + (7.37857e-05)*c*c + (-1.36894e-06)*c*c*c + (8.28054e-09)*c*c*c*c ;
	v2s=-0.0121719 + (0.00188256)*c + (-8.67104e-05)*c*c + (1.50514e-06)*c*c*c + (-8.79759e-09)*c*c*c*c ;
    }
    if(fRunIndex==133){
	v2c=0.00242068 + (-0.000274576)*c + (1.06631e-05)*c*c + (-1.45829e-07)*c*c*c + (7.18001e-10)*c*c*c*c ;
	v2s=0.0113182 + (-0.000529549)*c + (-7.75839e-06)*c*c + (4.86382e-07)*c*c*c + (-3.98634e-09)*c*c*c*c ;
    }
    if(fRunIndex==134){
	v2c=-0.0438721 + (0.00906564)*c + (-0.000407471)*c*c + (6.6254e-06)*c*c*c + (-3.61035e-08)*c*c*c*c ;
	v2s=-0.276264 + (0.0107322)*c + (-0.000281971)*c*c + (3.27634e-06)*c*c*c + (-1.17426e-08)*c*c*c*c ;
    }
    if(fRunIndex==135){
	v2c=0.0224469 + (-0.00324762)*c + (0.000165861)*c*c + (-3.16496e-06)*c*c*c + (2.011e-08)*c*c*c*c ;
	v2s=-0.308891 + (0.0136721)*c + (-0.00032343)*c*c + (2.93575e-06)*c*c*c + (-5.7637e-09)*c*c*c*c ;
    }
    if(fRunIndex==136){
	v2c=0.0279685 + (-0.00441)*c + (0.000294342)*c*c + (-6.56635e-06)*c*c*c + (4.53669e-08)*c*c*c*c ;
	v2s=-0.273287 + (0.00917778)*c + (-0.000144155)*c*c + (1.59881e-07)*c*c*c + (8.63874e-09)*c*c*c*c ;
    }

    /*
     if(fRunIndex>=1 && fRunIndex<=45){ // 137161 -137848 //period 1
     //    v2c =-3.523975e-03+8.132703e-04*x+-4.107765e-05*x*x+7.308763e-07*x*x*x+-4.285144e-09*x*x*x*x;
     //    v2s =-1.555676e-02+1.075371e-03*x+-3.597818e-05*x*x+4.927188e-07*x*x*x+-2.217080e-09*x*x*x*x;
     v2c=-2.007559e-05+3.143290e-04*x-2.065995e-05*x*x+4.071671e-07*x*x*x-2.562114e-09*x*x*x*x ;
     v2s=-4.920033e-03+6.807606e-04*x-3.065278e-05*x*x+5.132261e-07*x*x*x-2.767644e-09*x*x*x*x ;

  }
  if(fRunIndex>=46 && fRunIndex<=58){ //138125-138275//period 2
    v2c = 1.591471e-02+2.772445e-04*x-2.229627e-05*x*x+3.803341e-07*x*x*x-2.056551e-09*x*x*x*x;
    v2s =-5.800393e-02+2.007102e-03*x-3.292094e-05*x*x+7.266378e-08*x*x*x+1.828472e-09*x*x*x*x ;
  }
  if(fRunIndex>=59 && fRunIndex<=78){ //138359-138730//period 3
    v2c =-8.010228e-03+8.947184e-04*x-4.192609e-05*x*x+7.388210e-07*x*x*x-4.280541e-09*x*x*x*x ;
    v2s =-3.203888e-02+2.122650e-03*x-7.549389e-05*x*x+1.161380e-06*x*x*x-6.156621e-09*x*x*x*x ;
  }
  if(fRunIndex>=79){//period 4
    v2c =-1.025964e-02+1.128454e-03*x-4.603917e-05*x*x+7.386637e-07*x*x*x-4.035342e-09*x*x*x*x ;
    v2s =-9.746961e-02+4.672675e-03*x-1.311440e-04*x*x+1.650922e-06*x*x*x-7.083187e-09*x*x*x*x ;
  }
*/  
    return phi+v2c*TMath::Sin(2.*phi)-v2s*TMath::Cos(2.*phi) ;
}
//____________________________________________________________________________
Double_t  AliAnalysisTaskPi0v2::ApplyFlatteningTPC(Double_t phi, Double_t c){
  
 
  if(fRunIndex==39){ //137748
    Double_t v2c= 0.0222974-0.00132297*c+0.000204021*c*c-4.49827e-06*c*c*c+2.76986e-08*c*c*c*c ; 
    Double_t v2s=0.00615401-7.76462e-05*c ;   
    return phi+v2c*TMath::Sin(2.*phi)-v2s*TMath::Cos(2.*phi) ;
  }
  
  
  //Periods 1,2,3
  //fRunIndex - run index
  if(fRunIndex>=1 && fRunIndex<79){
    Double_t v2c=4.40516e-04+TMath::Exp(-4.71923-7.62089e-02*c) ;
    Double_t v2s=1.79859e-03+TMath::Exp(-4.99649-7.87523e-02*c) ;
   
    return phi+v2c*TMath::Sin(2.*phi)-v2s*TMath::Cos(2.*phi) ;
  }
  
  //period4
   //138826, 138828, 138830
  if(fRunIndex==83 ||fRunIndex==84 ||fRunIndex==85 ){
     Double_t v2c = -1.4e-03 ;
     Double_t v2s =  -5.60117e-02-TMath::Exp( -1.62827e+00-9.96071e-02*c);
    return phi + v2c*TMath::Sin(2.*phi) - v2s*TMath::Cos(2.*phi) ;
  }
  //Runs 138871,138872
  if(fRunIndex==87 ||fRunIndex==88 ){   
    Double_t v2c=-0.00518371 ;
    Double_t v2s= 0.00633324 ;
    return phi + v2c*TMath::Sin(2.*phi) - v2s*TMath::Cos(2.*phi) ;
  }   
  //run 139029
  if(fRunIndex==102 ){   
     Double_t v2c = -0.00354633 ;
     Double_t v2s =  0.00418512 ;
     return phi+v2c*TMath::Sin(2.*phi)-v2s*TMath::Cos(2.*phi) ;
  }  
  //run 139110  
/*
  if(fRunIndex==110 ){   
    Double_t v2c=4.78327e-04+TMath::Exp(-3.21625    -6.74306e-02*c) ;
    Double_t v2s=3.50176e-02+TMath::Exp(-9.73541e-01-9.19214e-02*c) ;
    return phi + v2c*TMath::Sin(2.*phi) - v2s*TMath::Cos(2.*phi) ;
  }
*/
  if(fRunIndex>=79 && fRunIndex<=139){
    Double_t v2c= 4.78327e-04 + TMath::Exp(-6.70587-2.33120e-02*c) ;
    Double_t v2s=-2.57731e-03 + TMath::Exp(-2.75493-1.05166e-01*c) ;
   
    return phi + v2c*TMath::Sin(2.*phi) - v2s*TMath::Cos(2.*phi) ;
  }
  return phi ;  
  
}  
