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

#include "AliEPSelectionTask.h"

// Author Daniel Lohner (Daniel.Lohner@cern.ch)

using namespace std;

ClassImp(AliAnalysisTaskPi0v2)

//________________________________________________________________________
    AliAnalysisTaskPi0v2::AliAnalysisTaskPi0v2(const char *name,Int_t harmonic) : AliAnalysisTaskSE(name),
    fV0Reader(NULL),
    fNCuts(0),
    fConversionSelection(NULL),
    fConversionGammas(NULL),
    fNCentralityBins(1),
    fCentralityBins(NULL),
    fCentrality(-1),
    fCentralityBin(0),
    fNBinsPhi(6),
    fEP(NULL),
    fUseTPCOnlyTracks(kTRUE),
    fEtaMax(0.75),
    fEtaGap(1),
    fRPTPCEtaA(0),
    fRPTPCEtaC(0),
    fRPV0A(0),
    fRPV0C(0),
    fRPTPC(0),
    fRPTPCEtaABF(0),
    fRPTPCEtaCBF(0),
    fRPV0ABF(0),
    fRPV0CBF(0),
    fRPTPCBF(0),
    fConversionCuts(NULL),
    fRandomizer(NULL),
    fOutputList(NULL),
    fMesonPDGCode(kPi0),
    fDeltaPsiRP(0),
    fRunNumber(0),
    fRunIndex(0),
    fNEPMethods(knEPMethod),
    fFillQA(kTRUE),
    fHarmonic(harmonic),
    fPsiMax(2*TMath::Pi()/Double_t(harmonic)),
    fPeriod("LHC10h"),
    fIsAOD(kFALSE),
    fSparseDist(NULL),
    fHruns(NULL),
    fDoEPFlattening(kTRUE),

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
    hCos2V0ATPCEtaA(NULL),
    hCos2V0ATPCEtaC(NULL),
    hCos2V0CTPCEtaA(NULL),
    hCos2V0CTPCEtaC(NULL),
    hCos2SumWeights(NULL),
    hEtaTPCEP(NULL),
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
    //hEPVertex(NULL)
    hEPQA(NULL)
{
    fInvMassRange[0]=0.05;
    fInvMassRange[1]=0.3;

    for(Int_t ii=0;ii<knEPMethod;ii++)fEPSelectionMask[ii]=1;

    fRandomizer=new TRandom3();
    fRandomizer->SetSeed(0);

    for(Int_t i = 0; i < 4; ++i) {
	fPhiDist[i] = 0;
    }

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

    if(fConversionSelection){
	for(Int_t ii=0;ii<fNCuts;ii++)delete fConversionSelection[ii];
	delete[] fConversionSelection;
	fConversionSelection=NULL;
    }

    if (fPeriod.CompareTo("LHC11h")==0){
	for(Int_t i = 0; i < 4; i++) {
	    if(fPhiDist[i]){
		delete fPhiDist[i];
		fPhiDist[i] = 0;
	    }
	}
	if(fHruns) delete fHruns;
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

    Int_t nbinspi0[knbinsPi0]={kGCnYBinsSpectra,kGCnXBinsSpectra,fNBinsPhi,fNCentralityBins,fNEPMethods};
    Double_t minpi0[knbinsPi0]={kGCfirstYBinSpectra,kGCfirstXBinSpectra,0,-0.5,-0.5};
    Double_t maxpi0[knbinsPi0]={kGClastYBinSpectra,kGClastXBinSpectra,fPsiMax/2.,fNCentralityBins-0.5,fNEPMethods-0.5};
    const char *binpi0[knbinsPi0]={"pt","mass","dPhi","centr","EPm"};

    Int_t nbinsg[knbinsGamma]={kGCnYBinsSpectra,fNBinsPhi,fNCentralityBins,fNEPMethods};
    Double_t ming[knbinsGamma]={kGCfirstYBinSpectra,0,-0.5,-0.5};
    Double_t maxg[knbinsGamma]={kGClastYBinSpectra,fPsiMax/2.,fNCentralityBins-0.5,fNEPMethods-0.5};
    const char *bingamma[knbinsGamma]={"pt","dPhi","centr","EPm"};

    // Define Binning

    if(!IsMC){

	hPi0=new THnSparseF*[fNCuts+1];
	hPi0BG=new THnSparseF*[fNCuts+1];
	hGamma=new THnSparseF*[fNCuts+1];

	// Photon Cuts
	for(Int_t iCut=0;iCut<fNCuts;iCut++){

	    TList *fCutOutputList=new TList();
	    fCutOutputList->SetName(fConversionSelection[iCut]->GetCutString().Data());
	    fCutOutputList->SetOwner(kTRUE);
	    fOutputList->Add(fCutOutputList);

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

	// no EP Flattening
	Int_t iCut=0;

	TList *fCutOutputList=new TList();
	fCutOutputList->SetName(Form("%s_BF",fConversionSelection[iCut]->GetCutString().Data()));
	fCutOutputList->SetOwner(kTRUE);
	fOutputList->Add(fCutOutputList);

	iCut=fNCuts;

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

    if(IsHeavyIon&&!IsMC){

	// RP Calculation
	TList *fRPList=new TList();
	fRPList->SetName("Event Plane");
	fRPList->SetOwner(kTRUE);
	fOutputList->Add(fRPList);

	hRPTPC=new TH2F("TPCAC" ,"TPC_AC" , fNCentralityBins,fCentralityBins, 180, 0, fPsiMax);
	hRPTPC->Sumw2();
	fRPList->Add(hRPTPC);
	hRPTPCEtaA=new TH2F("TPCEtaA" ,"TPC_EtaA" , fNCentralityBins,fCentralityBins, 180, 0, fPsiMax);
	hRPTPCEtaA->Sumw2();
	fRPList->Add(hRPTPCEtaA);
	hRPTPCEtaC=new TH2F("TPCEtaC" ,"TPC_EtaC" , fNCentralityBins,fCentralityBins, 180, 0, fPsiMax);
	hRPTPCEtaC->Sumw2();
	fRPList->Add(hRPTPCEtaC);
	hRPV0A=new TH2F("V0A" ,"VZERO_A" , fNCentralityBins,fCentralityBins, 180, 0, fPsiMax);
	hRPV0A->Sumw2();
	fRPList->Add(hRPV0A);
	hRPV0C=new TH2F("V0C" ,"VZERO_C" , fNCentralityBins,fCentralityBins, 180, 0, fPsiMax);
	hRPV0C->Sumw2();
	fRPList->Add(hRPV0C);

	hRPTPCAC=new TH2F("TPCA_TPCC" ,"TPCA_TPCC" , 180, 0, fPsiMax, 180, 0, fPsiMax);
	hRPTPCAC->Sumw2();
	fRPList->Add(hRPTPCAC);
	hRPV0ATPC=new TH2F("V0A_TPC" ,"V0A_TPC" , 180, 0, fPsiMax, 180, 0, fPsiMax);
	hRPV0ATPC->Sumw2();
	fRPList->Add(hRPV0ATPC);
	hRPV0CTPC=new TH2F("V0C_TPC" ,"V0C_TPC" , 180, 0, fPsiMax, 180, 0, fPsiMax);
	hRPV0CTPC->Sumw2();
	fRPList->Add(hRPV0CTPC);
	hRPV0AC=new TH2F("V0A_V0C" ,"V0A_V0C" , 180, 0, fPsiMax, 180, 0, fPsiMax);
	hRPV0AC->Sumw2();
	fRPList->Add(hRPV0AC);
	hRPTPCEtaAC=new TH2F("TPCEtaA_TPCEtaC" ,"TPCEtaA_TPCEtaC" , 180, 0, fPsiMax, 180, 0, fPsiMax);
	hRPTPCEtaAC->Sumw2();
	fRPList->Add(hRPTPCEtaAC);

        Int_t nsystep=4;// 3 different weights + unflattened EP

	hCos2TPC=new TH2F("Cos2_TPCAC" ,"Cos2_TPCAC" ,fNCentralityBins,fCentralityBins,nsystep+1,-0.5,nsystep+0.5);
	hCos2TPC->Sumw2();
	fRPList->Add(hCos2TPC);
	hCos2TPCEta=new TH2F("Cos2_TPCEtaAC" ,"Cos2_TPCEtaAC" ,fNCentralityBins,fCentralityBins,nsystep+1,-0.5,nsystep+0.5);
	hCos2TPCEta->Sumw2();
	fRPList->Add(hCos2TPCEta);
	hCos2V0ATPC=new TH2F("Cos2_V0ATPC" ,"Cos2_V0ATPC" ,fNCentralityBins,fCentralityBins,nsystep+1,-0.5,nsystep+0.5);
	hCos2V0ATPC->Sumw2();
	fRPList->Add(hCos2V0ATPC);
	hCos2V0CTPC=new TH2F("Cos2_V0CTPC" ,"Cos2_V0CTPC" ,fNCentralityBins,fCentralityBins,nsystep+1,-0.5,nsystep+0.5);
	hCos2V0CTPC->Sumw2();
	fRPList->Add(hCos2V0CTPC);
	hCos2V0AC=new TH2F("Cos2_V0AC" ,"Cos2_V0AC" ,fNCentralityBins,fCentralityBins,nsystep+1,-0.5,nsystep+0.5);
	hCos2V0AC->Sumw2();
	fRPList->Add(hCos2V0AC);
        hCos2V0ATPCEtaA=new TH2F("Cos2_V0ATPCEtaA" ,"Cos2_V0ATPCEtaA" ,fNCentralityBins,fCentralityBins,nsystep+1,-0.5,nsystep+0.5);
	hCos2V0ATPCEtaA->Sumw2();
	fRPList->Add(hCos2V0ATPCEtaA);
	hCos2V0ATPCEtaC=new TH2F("Cos2_V0ATPCEtaC" ,"Cos2_V0ATPCEtaC" ,fNCentralityBins,fCentralityBins,nsystep+1,-0.5,nsystep+0.5);
	hCos2V0ATPCEtaC->Sumw2();
	fRPList->Add(hCos2V0ATPCEtaC);
        hCos2V0CTPCEtaA=new TH2F("Cos2_V0CTPCEtaA" ,"Cos2_V0CTPCEtaA" ,fNCentralityBins,fCentralityBins,nsystep+1,-0.5,nsystep+0.5);
	hCos2V0CTPCEtaA->Sumw2();
	fRPList->Add(hCos2V0CTPCEtaA);
	hCos2V0CTPCEtaC=new TH2F("Cos2_V0CTPCEtaC" ,"Cos2_V0CTPCEtaC" ,fNCentralityBins,fCentralityBins,nsystep+1,-0.5,nsystep+0.5);
	hCos2V0CTPCEtaC->Sumw2();
	fRPList->Add(hCos2V0CTPCEtaC);
	hCos2SumWeights=new TH2F("Cos2_SumWeights" ,"Cos2_SumWeights" ,fNCentralityBins,fCentralityBins,nsystep+1,-0.5,nsystep+0.5);
	hCos2SumWeights->Sumw2();
	fRPList->Add(hCos2SumWeights);

	hEtaTPCEP=new TH2F("Eta_TPCEP" ,"Eta_TPCEP" ,fNCentralityBins,fCentralityBins,100,-1,1);
	hEtaTPCEP->Sumw2();
	fRPList->Add(hEtaTPCEP);


	/*const Int_t nepbins=4;
	 Int_t nbinsep[nepbins]={30,30,40,180};
	 Double_t minep[nepbins]={-0.015,0.17,-10,0};
	 Double_t maxep[nepbins]={0.015,0.20,10,fPsiMax};

	 hEPVertex=new THnSparseF("EP_Vertex","EP_Vertex",nepbins,nbinsep,minep,maxep);
	 fRPList->Add(hEPVertex);
	 */

	const Int_t nRuns=270;
	const Int_t nepbins=4;
	Int_t nbinsep[nepbins]={fNCentralityBins,180,nRuns,5};
	Double_t minep[nepbins]={-0.5,0,0.5,-0.5};
	Double_t maxep[nepbins]={fNCentralityBins-0.5,fPsiMax,Double_t(nRuns)+0.5,4.5};
	hEPQA=new THnSparseF("EP_QA","EP_QA",nepbins,nbinsep,minep,maxep);
	fRPList->Add(hEPQA);
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
	Int_t nbinscharged[knbinsCharged]={fNBinsPhi,fNCentralityBins,fNEPMethods};
	Double_t mincharged[knbinsCharged]={0,-0.5,-0.5};
	Double_t maxcharged[knbinsCharged]={fPsiMax/2,fNCentralityBins-0.5,fNEPMethods-0.5};
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
Bool_t AliAnalysisTaskPi0v2::InitEvent(){

    if(!fV0Reader){AliError("Error: No V0 Reader and Pi0 Reconstructor");return kFALSE;}
    if(!fV0Reader->IsEventSelected())return kFALSE;
    fConversionGammas=fV0Reader->GetReconstructedGammas();

    fIsAOD=(fInputEvent->IsA()==AliAODEvent::Class());

    if(!fConversionSelection){
	AliError("No Cut Selection");
	return kFALSE;
    }

    if(!SetCentrality()){return kFALSE;}

    if(fConversionCuts->IsHeavyIon()&&!fMCEvent){

	if(fRunNumber!=fInputEvent->GetRunNumber()){
	    fRunNumber=fInputEvent->GetRunNumber();
	    if (fRunNumber >= 136851 && fRunNumber <= 139515){fPeriod = "LHC10h";}
	    if (fRunNumber >= 166529 && fRunNumber <= 170593){fPeriod = "LHC11h";}
	    fRunIndex=GetRunIndex(fRunNumber);
	    LoadVZEROCalibration(fRunNumber); // Load Calibration for V0 Event Plane
            LoadTPCCalibration(fRunNumber); // Load Calibration for TPC Event Plane
	}
	if(fRunIndex<0){
	    AliInfo("Run not selected");
	    return kFALSE;
	}

	// TPC Event Plane
	if(!GetTPCEventPlane())return kFALSE;
	//fEP=fInputEvent->GetEventplane();
	if(!fEP)return kFALSE;
	fRPTPCBF=GetCorrectedTPCEPAngle(NULL,NULL,kFALSE);
	fRPTPC=GetEventPlaneAngle(kTPC);

	// TPC Eta Sub Events
	fRPTPCEtaABF=GetTPCSubEPEta(kEPTPCEtaA);
	fRPTPCEtaA=ApplyFlattening(fRPTPCEtaABF,kEPTPCEtaA);

	fRPTPCEtaCBF=GetTPCSubEPEta(kEPTPCEtaC);
	fRPTPCEtaC=ApplyFlattening(fRPTPCEtaCBF,kEPTPCEtaC);

	// GetV0 Event Plane
	GetV0EP(fInputEvent,fRPV0ABF,fRPV0CBF);
	fRPV0A=ApplyFlattening(fRPV0ABF,kEPV0A);
	fRPV0C=ApplyFlattening(fRPV0CBF,kEPV0C);

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

		    if(!fEPSelectionMask[iEP])continue; // dont fill THnSparse if not needed-> Save Memory

		    ProcessPi0s(iCut,EEventPlaneMethod(iEP));

		    ProcessGammas(iCut,EEventPlaneMethod(iEP));
		}
	    }

	    // QA
	    if(fFillQA&&iCut==0)ProcessQA();
	}
    }

    // Fill N Events
    hNEvents->Fill(fCentrality);

    // EventPlaneResolution
    ProcessEventPlane();

    PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::ProcessPi0s(Int_t iCut,EEventPlaneMethod iEP){

    if(!fConversionSelection[iCut])return;

    if(fConversionSelection[iCut]->GetNumberOfPhotons()==0)return;

    // Process Pi0s

    for(Int_t ii=0;ii<fConversionSelection[iCut]->GetNumberOfPi0s();ii++){

	AliAODConversionMother *pi0cand=fConversionSelection[iCut]->GetPi0(ii);

	if(!pi0cand)continue;

	Double_t val[knbinsPi0];
	val[kPi0Pt]=pi0cand->Pt();
	val[kPi0Mass]=pi0cand->M();
	val[kPi0dPhi]=GetPi0PhiwrtRP(pi0cand,iEP);
	val[kPi0Cent]=fCentralityBin;
	val[kPi0EPM]=Int_t(iEP);

	hPi0[iCut]->Fill(val);

	if(iCut==0){
	    // no flattening
	    val[kPi0dPhi]=GetPi0PhiwrtRP(pi0cand,iEP,kFALSE);
	    hPi0[fNCuts]->Fill(val);
	}
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

	hPi0BG[iCut]->Fill(val,pi0cand->GetWeight());

	if(iCut==0){
	    // no flattening
	    val[kPi0dPhi]=GetPi0PhiwrtRP(pi0cand,iEP,kFALSE);
	    hPi0BG[fNCuts]->Fill(val);
	}
    }
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::ProcessGammas(Int_t iCut,EEventPlaneMethod iEP){

    if(!fConversionSelection[iCut]){
	AliWarning("Conversion Selection does not exist");
	return;
    }

    for(Int_t ii=0;ii<fConversionSelection[iCut]->GetNumberOfPhotons();ii++){

	AliAODConversionPhoton *gamma=fConversionSelection[iCut]->GetPhoton(ii);

	Double_t val[knbinsGamma];
	val[kGammaPt]=gamma->Pt();
	val[kGammadPhi]=GetPhotonPhiwrtRP(gamma,iEP);
	val[kGammaCent]=fCentralityBin;
	val[kGammaEPM]=Int_t(iEP);

	hGamma[iCut]->Fill(val);

	if(iCut==0){
            // no flattening
	    val[kGammadPhi]=GetPhotonPhiwrtRP(gamma,iEP,kFALSE);
	    hGamma[fNCuts]->Fill(val);
	}

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

    if(!fConversionSelection[0]){
	AliWarning("Conversion Selection does not exist");
	return;
    }


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
       
	    hGammaMult[iEP]->Fill(val);
	    hGammaMultdPhi[iEP]->Fill(valdPhi);

	    // Gamma Phi
	    hGammaPhi[fCentralityBin]->Fill(gamma->Pt(),gamma->Phi());
	    hGammaMultCent->Fill(fCentrality,Float_t(fConversionSelection[0]->GetNumberOfPhotons()));

	    if(fMCStack){
		if(gamma->IsTruePhoton(fMCStack)){
		    hGammaMultRECOTRUE->Fill(val);
		    hGammaMultdPhiRECOTRUE->Fill(valdPhi);
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
       
			hGammaMultTRUE->Fill(val);
			hGammaMultdPhiTRUE->Fill(valdPhi);
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
    Double_t binrange=TMath::Pi()/(Double_t(fHarmonic*fNBinsPhi));
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
Double_t AliAnalysisTaskPi0v2::GetPi0PhiwrtRP(AliAODConversionMother *pi0,EEventPlaneMethod iEP,Bool_t bDoFlattening){

    AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(pi0->GetLabel1()));
    AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(pi0->GetLabel2()));

    Double_t EPAngle=GetEventPlaneAngle(iEP,pi0->Eta(),gamma0,gamma1,bDoFlattening);

    return GetPhiwrtRP(pi0->Phi()-EPAngle);
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetPhotonPhiwrtRP(AliAODConversionPhoton *gamma,EEventPlaneMethod iEP,Bool_t bDoFlattening){

    Double_t EPAngle=GetEventPlaneAngle(iEP,gamma->Eta(),gamma,NULL,bDoFlattening);

    return GetPhiwrtRP(gamma->Phi()-EPAngle);
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetMCPhotonPhiwrtRP(TParticle *gamma,EEventPlaneMethod iEP,Bool_t bDoFlattening){

    Double_t EPAngle=GetEventPlaneAngle(iEP,gamma->Eta(),NULL,NULL,bDoFlattening);

    return GetPhiwrtRP(gamma->Phi()-EPAngle);
}
//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetChargedPhiwrtRP(AliVTrack *track,EEventPlaneMethod iEP,Bool_t bDoFlattening){

    Double_t EPAngle=GetEventPlaneAngle(iEP,track->Eta(),NULL,NULL,bDoFlattening);

    return GetPhiwrtRP(track->Phi()-EPAngle);
}
//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetPhiwrtRP(Double_t dPhi){
    Double_t newdPhi=TMath::Abs(dPhi); // Cos is symmetric
    while(newdPhi>fPsiMax/2)newdPhi=TMath::Abs(newdPhi-fPsiMax);
    return newdPhi;
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::Terminate(Option_t *) 
{

}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::ProcessEventPlane()
{
    if(!fMCEvent&&fConversionCuts->IsHeavyIon()){

	/*  Double_t val[4];
	val[0]=fInputEvent->GetPrimaryVertex()->GetX();
	val[1]=fInputEvent->GetPrimaryVertex()->GetY();
	val[2]=fInputEvent->GetPrimaryVertex()->GetZ();
	val[3]=GetEventPlaneAngle(kTPC);
	hEPVertex->Fill(val);   */

	// Run by run monitoring (before flattening)

	Double_t val[4];
        val[0]=fCentralityBin;
	val[2]=fRunIndex;

	val[1]=fRPTPCBF;
	val[3]=kEPTPC;
	hEPQA->Fill(val);

	val[1]=fRPTPCEtaABF;
	val[3]=kEPTPCEtaA;
        hEPQA->Fill(val);
	val[1]=fRPTPCEtaCBF;
	val[3]=kEPTPCEtaC;
	hEPQA->Fill(val);

        val[1]=fRPV0ABF;
	val[3]=kEPV0A;
        hEPQA->Fill(val);
        val[1]=fRPV0CBF;
	val[3]=kEPV0C;
        hEPQA->Fill(val);

	// After Flattening

	// TPC EP
        Double_t PsiRP1BF=fEP->GetQsub1()->Phi()/Double_t(fHarmonic);
	Double_t PsiRP2BF=fEP->GetQsub2()->Phi()/Double_t(fHarmonic);
	Double_t PsiRP1=ApplyFlattening(PsiRP1BF,kEPTPC);
	Double_t PsiRP2=ApplyFlattening(PsiRP2BF,kEPTPC);

	hRPTPC->Fill(fCentrality,fRPTPC);
	hRPTPCAC->Fill(PsiRP1,PsiRP2);
       
	// TPC Eta Gap
	hRPTPCEtaA->Fill(fCentrality,fRPTPCEtaA);
        hRPTPCEtaC->Fill(fCentrality,fRPTPCEtaC);
	hRPTPCEtaAC->Fill(fRPTPCEtaA,fRPTPCEtaC);
     
	// V0
      
	hRPV0A->Fill(fCentrality,fRPV0A);
	hRPV0C->Fill(fCentrality,fRPV0C);

	hRPV0ATPC->Fill(fRPV0A,fRPTPC);
	hRPV0CTPC->Fill(fRPV0C,fRPTPC);
	hRPV0AC->Fill(fRPV0A,fRPV0C);

	Double_t cos2V0ATPC=TMath::Cos(Double_t(fHarmonic)*(fRPTPC-fRPV0A));
        Double_t cos2V0CTPC=TMath::Cos(Double_t(fHarmonic)*(fRPTPC-fRPV0C));
	Double_t cos2V0AV0C=TMath::Cos(Double_t(fHarmonic)*(fRPV0C-fRPV0A));
	Double_t cos2TPCEta=TMath::Cos(Double_t(fHarmonic)*(fRPTPCEtaA-fRPTPCEtaC));
        Double_t cos2TPC=TMath::Cos(Double_t(fHarmonic)*(PsiRP1-PsiRP2));
	Double_t cos2V0ATPCEtaA=TMath::Cos(Double_t(fHarmonic)*(fRPV0A-fRPTPCEtaA));
	Double_t cos2V0CTPCEtaA=TMath::Cos(Double_t(fHarmonic)*(fRPV0C-fRPTPCEtaA));
        Double_t cos2V0ATPCEtaC=TMath::Cos(Double_t(fHarmonic)*(fRPV0A-fRPTPCEtaC));
        Double_t cos2V0CTPCEtaC=TMath::Cos(Double_t(fHarmonic)*(fRPV0C-fRPTPCEtaC));

        const Int_t nfill=4;
	Double_t weight[nfill];
	weight[0]=1.;// Fill unweighted
	weight[1]=Float_t(fConversionSelection[0]->GetNumberOfPhotons());// Weight with Photon Mult
	weight[2]=Float_t(fConversionSelection[0]->GetNumberOfChargedTracks(fInputEvent)); // Weight with charged Track Mult
	weight[3]=Float_t(fConversionSelection[0]->GetVZEROMult(fInputEvent)); // Weight with V0 mult

	for(Int_t i=0;i<nfill;i++){

	    hCos2V0ATPC->Fill(fCentrality,i,cos2V0ATPC*weight[i]);
	    hCos2V0CTPC->Fill(fCentrality,i,cos2V0CTPC*weight[i]);
	    hCos2V0AC->Fill(fCentrality,i,cos2V0AV0C*weight[i]);
	    hCos2TPCEta->Fill(fCentrality,i,cos2TPCEta*weight[i]);
	    hCos2TPC->Fill(fCentrality,i,cos2TPC*weight[i]);
	    hCos2V0ATPCEtaA->Fill(fCentrality,i,cos2V0ATPCEtaA*weight[i]);
	    hCos2V0ATPCEtaC->Fill(fCentrality,i,cos2V0ATPCEtaC*weight[i]);
	    hCos2V0CTPCEtaA->Fill(fCentrality,i,cos2V0CTPCEtaA*weight[i]);
	    hCos2V0CTPCEtaC->Fill(fCentrality,i,cos2V0CTPCEtaC*weight[i]);

	    hCos2SumWeights->Fill(fCentrality,i,weight[i]);
	}

        // Fill Resolution before EP Flattening
	Double_t cos2V0ATPCBF=TMath::Cos(Double_t(fHarmonic)*(fRPTPCBF-fRPV0ABF));
        Double_t cos2V0CTPCBF=TMath::Cos(Double_t(fHarmonic)*(fRPTPCBF-fRPV0CBF));
	Double_t cos2V0AV0CBF=TMath::Cos(Double_t(fHarmonic)*(fRPV0CBF-fRPV0ABF));
	Double_t cos2TPCEtaBF=TMath::Cos(Double_t(fHarmonic)*(fRPTPCEtaABF-fRPTPCEtaCBF));
	Double_t cos2TPCBF=TMath::Cos(Double_t(fHarmonic)*(PsiRP1BF-PsiRP2BF));
	Double_t cos2V0ATPCEtaABF=TMath::Cos(Double_t(fHarmonic)*(fRPV0ABF-fRPTPCEtaABF));
	Double_t cos2V0CTPCEtaABF=TMath::Cos(Double_t(fHarmonic)*(fRPV0CBF-fRPTPCEtaABF));
	Double_t cos2V0ATPCEtaCBF=TMath::Cos(Double_t(fHarmonic)*(fRPV0ABF-fRPTPCEtaCBF));
	Double_t cos2V0CTPCEtaCBF=TMath::Cos(Double_t(fHarmonic)*(fRPV0CBF-fRPTPCEtaCBF));

	hCos2V0ATPC->Fill(fCentrality,nfill,cos2V0ATPCBF);
	hCos2V0CTPC->Fill(fCentrality,nfill,cos2V0CTPCBF);
	hCos2V0AC->Fill(fCentrality,nfill,cos2V0AV0CBF);
	hCos2TPCEta->Fill(fCentrality,nfill,cos2TPCEtaBF);
	hCos2TPC->Fill(fCentrality,nfill,cos2TPCBF);
	hCos2V0ATPCEtaA->Fill(fCentrality,nfill,cos2V0ATPCEtaABF);
	hCos2V0ATPCEtaC->Fill(fCentrality,nfill,cos2V0ATPCEtaCBF);
	hCos2V0CTPCEtaA->Fill(fCentrality,nfill,cos2V0CTPCEtaABF);
	hCos2V0CTPCEtaC->Fill(fCentrality,nfill,cos2V0CTPCEtaCBF);

	hCos2SumWeights->Fill(fCentrality,nfill);

    }
}

//________________________________________________________________________
TVector2 AliAnalysisTaskPi0v2::GetEPContribution(AliAODConversionPhoton *gamma){
    TVector2 q;
    for(Int_t ii=0;ii<2;ii++){
	AliVTrack *fCurrentTrack=AliConversionCuts::GetTrack(fInputEvent,gamma->GetTrackLabel(ii));
	TVector2 qtrack=GetContributionEP(fCurrentTrack);
	q+=qtrack;
    }
    return q;
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetCorrectedTPCEPAngle(AliAODConversionPhoton *gamma0,AliAODConversionPhoton *gamma1,Bool_t bDoFlattening){
    // Correct Event Plane for Dilepton Tracks
    TVector2 q0(*fEP->GetQVector());
    if(gamma0)q0-=GetEPContribution(gamma0);
    if(gamma1)q0-=GetEPContribution(gamma1);
    Double_t EPangle=GetPsiInRange(q0.Phi()/Double_t(fHarmonic));
    if(bDoFlattening)EPangle=ApplyFlattening(EPangle,kEPTPC);

    return EPangle;
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetTPCSubEPEta(EEventPlane ep){

    Double_t etamin,etamax;
    switch(ep){
    case kEPTPCEtaA:
	etamin=fEtaGap/2;
	etamax=1;
	break;
    case kEPTPCEtaC:
	etamin=-1;
	etamax=-fEtaGap/2;
        break;
    default:
	return 0;
    }

    TVector2 q;
    for(Int_t ii=0;ii<fInputEvent->GetNumberOfTracks();ii++){
	AliVTrack *fCurrentTrack=dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(ii));
	if(!fCurrentTrack)continue;
	if(fCurrentTrack->Eta()>=etamin&&fCurrentTrack->Eta()<=etamax){
	    TVector2 qtrack=GetContributionEP(fCurrentTrack);
	    q+=qtrack;
	}
    }

    Double_t phi=GetPsiInRange(q.Phi()/Double_t(fHarmonic));

    return phi;
}

//________________________________________________________________________
TVector2 AliAnalysisTaskPi0v2::GetContributionEP(AliVTrack *track){

    TVector2 q(0,0);

    TArrayF *fQContributionX=fEP->GetQContributionXArray();
    TArrayF *fQContributionY=fEP->GetQContributionYArray();

    Int_t trackID=track->GetID();

    if(fIsAOD){
	if((trackID>-1)&&fUseTPCOnlyTracks)return q;
	if((trackID<0)&&!fUseTPCOnlyTracks)return q;
	if (fUseTPCOnlyTracks) trackID = trackID*(-1) - 1;
    }

    q.Set(fQContributionX->GetAt(trackID),fQContributionY->GetAt(trackID));

    return q;
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetPsiInRange(Double_t phi){

    Double_t newphi=phi;
    while(newphi<0)newphi+=fPsiMax;
    while(newphi>fPsiMax)newphi-=fPsiMax;
    return newphi;
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetEventPlaneAngle(EEventPlaneMethod EPmethod,Double_t eta,AliAODConversionPhoton *gamma0,AliAODConversionPhoton *gamma1,Bool_t bDoFlattening)
{
    // If arguments are not null, the contribution of these photons is subtracted from the TPC EP

    if(fConversionCuts->IsHeavyIon()){

	// For MC select random EP angle in order to avoid correlations due to azimuth dependent efficiencies (ITS holes)
	if(fMCEvent){
	    return fRandomizer->Uniform(0,fPsiMax);
	}
	switch(EPmethod){
	case kTPC:
	    return GetCorrectedTPCEPAngle(gamma0,gamma1,bDoFlattening);
	case kTPCEtaGap:
	    // Use opposite EP
	    if(bDoFlattening){
		if(eta<0)return fRPTPCEtaA; 
		else return fRPTPCEtaC;
	    }
	    else{
		if(eta<0)return fRPTPCEtaABF;
		else return fRPTPCEtaCBF;
	    }
	case kV0A:
	    if(bDoFlattening)return fRPV0A;
	    else return fRPV0ABF;
	case kV0C:
	    if(bDoFlattening)return fRPV0C;
	    else return fRPV0CBF;
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
void AliAnalysisTaskPi0v2::SetCuts(AliConversionSelection **conversionselection,Int_t ncuts){

    if(fConversionSelection){
	for(Int_t ii=0;ii<fNCuts;ii++)delete fConversionSelection;
        delete[] fConversionSelection;
        fConversionSelection=NULL;
    }
    fNCuts=ncuts;
    fConversionSelection=new AliConversionSelection*[fNCuts];
    for(Int_t ii=0;ii<fNCuts;ii++)fConversionSelection[ii]=new AliConversionSelection(*conversionselection[ii]);
}

//___________________________________________________________________________
Int_t AliAnalysisTaskPi0v2::GetRunIndex(Int_t run){

    switch(run){

	//LHC11h (131 runs)
    case 167902 : return 140;
    case 167903 : return 141;
    case 167909 : return 142;
    case 167915 : return 143;
    case 167920 : return 144;
    case 167985 : return 145;
    case 167986 : return 146;
    case 167987 : return 147;
    case 167988 : return 148;
    case 168066 : return 149;
    case 168068 : return 150;
    case 168069 : return 151;
    case 168076 : return 152;
    case 168103 : return 153;
    case 168104 : return 154;
    case 168105 : return 155;
    case 168107 : return 156;
    case 168108 : return 157;
    case 168115 : return 158;
    case 168212 : return 159;
    case 168310 : return 160;
    case 168311 : return 161;
    case 168322 : return 162;
    case 168325 : return 163;
    case 168341 : return 164;
    case 168342 : return 165;
    case 168361 : return 166;
    case 168362 : return 167;
    case 168458 : return 168;
    case 168460 : return 169;
    case 168461 : return 170;
    case 168464 : return 171;
    case 168467 : return 172;
    case 168511 : return 173;
    case 168512 : return 174;
    case 168514 : return 175;
    case 168777 : return 176;
    case 168826 : return 177;
    case 168984 : return 178;
    case 168988 : return 179;
    case 168992 : return 180;
    case 169035 : return 181;
    case 169040 : return 182;
    case 169044 : return 183;
    case 169045 : return 184;
    case 169091 : return 185;
    case 169094 : return 186;
    case 169099 : return 187;
    case 169138 : return 188;
    case 169143 : return 189;
    case 169144 : return 190;
    case 169145 : return 191;
    case 169148 : return 192;
    case 169156 : return 193;
    case 169160 : return 194;
    case 169167 : return 195;
    case 169238 : return 196;
    case 169411 : return 197;
    case 169415 : return 198;
    case 169417 : return 199;
    case 169418 : return 200;
    case 169419 : return 201;
    case 169420 : return 202;
    case 169475 : return 203;
    case 169498 : return 204;
    case 169504 : return 205;
    case 169506 : return 206;
    case 169512 : return 207;
    case 169515 : return 208;
    case 169550 : return 209;
    case 169553 : return 210;
    case 169554 : return 211;
    case 169555 : return 212;
    case 169557 : return 213;
    case 169584 : return 214;
    case 169586 : return 215;
    case 169587 : return 216;
    case 169588 : return 217;
    case 169590 : return 218;
    case 169591 : return 219;
    case 169835 : return 220;
    case 169837 : return 221;
    case 169838 : return 222;
    case 169846 : return 223;
    case 169855 : return 224;
    case 169858 : return 225;
    case 169859 : return 226;
    case 169922 : return 227;
    case 169923 : return 228;
    case 169956 : return 229;
    case 169965 : return 230;
    case 169975 : return 231;
    case 169981 : return 232;
    case 170027 : return 233;
    case 170036 : return 234;
    case 170038 : return 235;
    case 170040 : return 236;
    case 170081 : return 237;
    case 170083 : return 238;
    case 170084 : return 239;
    case 170085 : return 240;
    case 170088 : return 241;
    case 170089 : return 242;
    case 170091 : return 243;
    case 170152 : return 244;
    case 170155 : return 245;
    case 170159 : return 246;
    case 170163 : return 247;
    case 170193 : return 248;
    case 170195 : return 249;
    case 170203 : return 250;
    case 170204 : return 251;
    case 170207 : return 252;
    case 170208 : return 253;
    case 170228 : return 254;
    case 170230 : return 255;
    case 170268 : return 256;
    case 170269 : return 257;
    case 170270 : return 258;
    case 170306 : return 259;
    case 170308 : return 260;
    case 170309 : return 261;
    case 170311 : return 262;
    case 170312 : return 263;
    case 170313 : return 264;
    case 170315 : return 265;
    case 170387 : return 266;
    case 170388 : return 267;
    case 170556 : return 268;
    case 170572 : return 269;
    case 170593 : return 270;

    //LHC10h (137 runs)
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
    // case  138579 : return 68;
    // case  138578 : return 67;
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

    // Default
    default : return -1;
    }
}

//____________________________________________________________________________
void AliAnalysisTaskPi0v2::GetV0EP(AliVEvent * event,Double_t &rpv0a,Double_t &rpv0c){

    // Corrected VZERO EP (from AliAnalysisTaskPi0Flow)

    //VZERO data
    AliESDVZERO* esdV0 = (AliESDVZERO*)event->GetVZEROData();

    //reset Q vector info
    Double_t Qxa2 = 0, Qya2 = 0;
    Double_t Qxc2 = 0, Qyc2 = 0;

    for (Int_t iv0 = 0; iv0 < 64; iv0++) {
	Double_t phiV0 = TMath::PiOver4()*(0.5 + iv0 % 8);
	Float_t multv0 = esdV0->GetMultiplicity(iv0);
	Double_t lqx=TMath::Cos(fHarmonic*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	Double_t lqy=TMath::Sin(fHarmonic*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	if (iv0 < 32){ // V0C
	    Qxc2 += lqx;
	    Qyc2 += lqy;
	} else {       // V0A
	    Qxa2 += lqx;
	    Qya2 += lqy;
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
    rpv0a = TMath::ATan2(QyaCor2, QxaCor2)/Double_t(fHarmonic);
    rpv0c = TMath::ATan2(QycCor2, QxcCor2)/Double_t(fHarmonic);

    //rpv0a = TMath::ATan2(Qya2, Qxa2)/Double_t(fHarmonic);
    //rpv0c = TMath::ATan2(Qyc2, Qxc2)/Double_t(fHarmonic);

   // cout<<"Compare v"<<fHarmonic<<" "<<rpv0a<<" "<<fInputEvent->GetEventplane()->GetEventplane("V0A",fInputEvent,fHarmonic)<<endl;

    rpv0a=GetPsiInRange(rpv0a);
    rpv0c=GetPsiInRange(rpv0c);

}

//_____________________________________________________________________________
void AliAnalysisTaskPi0v2::LoadVZEROCalibration(Int_t run){

    // VZERO Phi Weights and Recentering

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
    fMultV0->Fit(fpol0,"0","",0,31);
    fV0Cpol = fpol0->GetParameter(0);
    fMultV0->Fit(fpol0,"0","",32,64);
    fV0Apol = fpol0->GetParameter(0);

    for(Int_t iside=0;iside<2;iside++){
	for(Int_t icoord=0;icoord<2;icoord++){
	    for(Int_t i=0;i  < nCentrBinV0;i++){
		char namecont[100];

		if(iside==0 && icoord==0)
		    snprintf(namecont,100,"hQxc%i_%i",fHarmonic,i);
		else if(iside==1 && icoord==0)
		    snprintf(namecont,100,"hQxa%i_%i",fHarmonic,i);
		else if(iside==0 && icoord==1)
		    snprintf(namecont,100,"hQyc%i_%i",fHarmonic,i);
		else if(iside==1 && icoord==1)
		    snprintf(namecont,100,"hQya%i_%i",fHarmonic,i);

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
	    }
	}
    }
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0v2::LoadTPCCalibration(Int_t run){

    // TPC Event Plane Weights
    AliOADBContainer *fEPContainer=NULL;
    TString oadbfilename="";

    if (run >= 136851 && run <= 139515){
	// LHC10h

	if(fIsAOD){
	    oadbfilename = (Form("%s/COMMON/EVENTPLANE/data/epphidist.aod.root", AliAnalysisManager::GetOADBPath()));
	}
	else{
	    oadbfilename = (Form("%s/COMMON/EVENTPLANE/data/epphidist.root", AliAnalysisManager::GetOADBPath()));
	}

	TFile foadb(oadbfilename);
	if(!foadb.IsOpen()) AliFatal(Form("Cannot open OADB file %s", oadbfilename.Data()));

	AliInfo("Using Standard OADB");
	fEPContainer = (AliOADBContainer*) foadb.Get("epphidist");
	if (!fEPContainer) AliFatal("Cannot fetch OADB container for EP selection");
	foadb.Close();

	fPhiDist[0] = (TH1F*) fEPContainer->GetObject(fRunNumber, "Default");

	Bool_t emptybins;

	int iter = 0;
	while (iter<3){
	    emptybins = kFALSE;

	    for (int i=1; i<fPhiDist[0]->GetNbinsX(); i++){
		if (!((fPhiDist[0]->GetBinContent(i))>0)) {
		    emptybins = kTRUE;
		}
	    }
	    if (emptybins) {
		cout << "empty bins - rebinning!" << endl;
		fPhiDist[0]->Rebin();
		iter++;
	    }
	    else iter = 3;
	}

	if (emptybins) {
	    AliError("After Maximum of rebinning still empty Phi-bins!!!");
	}
    }

    if (run >= 166529 && run <= 170593){
        // LHC11h

	oadbfilename = (Form("%s/COMMON/EVENTPLANE/data/epphidist2011.root", AliAnalysisManager::GetOADBPath()));
	TFile *foadb = TFile::Open(oadbfilename);
	if(!foadb->IsOpen()) AliFatal(Form("Cannot open OADB file %s", oadbfilename.Data()));

	AliInfo("Using Standard OADB");
	fSparseDist = (THnSparse*) foadb->Get("Default");
	if (!fSparseDist) AliFatal("Cannot fetch OADB container for EP selection");
	foadb->Close();
	if(!fHruns){
	  fHruns = (TH1F*)fSparseDist->Projection(0); //projection on run axis;
           fHruns->SetName("runsHisto");
	}

	Int_t runbin=fHruns->FindBin(fRunNumber);
	if (fHruns->GetBinContent(runbin) > 1){
	    fSparseDist->GetAxis(0)->SetRange(runbin,runbin);
	}
	else if(fHruns->GetBinContent(runbin) < 2){
	    fSparseDist->GetAxis(0)->SetRange(1,2901); // not calibrated run, use integrated phi-weights
	    AliInfo("Using integrated Phi-weights for this run");
	}
	for (Int_t i = 0; i<4 ;i++)
	{
	    if(fPhiDist[i]){
		delete fPhiDist[i];
		fPhiDist[i] = 0x0;
	    }
	    if(i == 0){
		fSparseDist->GetAxis(1)->SetRange(1,1);  // neg charge
		fSparseDist->GetAxis(2)->SetRange(1,1);} // neg eta
	    if(i == 1){
		fSparseDist->GetAxis(1)->SetRange(2,2);  // pos charge
		fSparseDist->GetAxis(2)->SetRange(1,1);} // neg eta
	    if(i == 2){
		fSparseDist->GetAxis(1)->SetRange(1,1);  // neg charge
		fSparseDist->GetAxis(2)->SetRange(2,2);} // pos eta
	    if(i == 3){
		fSparseDist->GetAxis(1)->SetRange(2,2);  // pos charge
		fSparseDist->GetAxis(2)->SetRange(2,2);} // pos eta
	    fPhiDist[i] = (TH1F*)fSparseDist->Projection(3); // Projection on Phi
	    fPhiDist[i]->SetName(Form("phidist%d%d",i,fRunNumber));
	    fSparseDist->GetAxis(1)->SetRange(1,2); // reset axes
	    fSparseDist->GetAxis(2)->SetRange(1,2);
	}
	fSparseDist->GetAxis(0)->SetRange(1,2901);// reset run axis
    }

    if (!fPhiDist[0]) AliFatal(Form("Cannot find OADB phi distribution for run %d", run));

}

//_________________________________________________________________________
Int_t AliAnalysisTaskPi0v2::GetAODEPTrackFilterBit(){

    if(fUseTPCOnlyTracks){
	return 128;// TPC only with vertex constraint
    }
    return 1;// Use Global Tracks
}

//_________________________________________________________________________
TObjArray* AliAnalysisTaskPi0v2::GetEventPlaneTracks(Int_t &maxID)
{
    TObjArray *tracklist=NULL;

    AliESDtrackCuts *fEPESDtrackCuts=AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEPESDtrackCuts->SetPtRange(0.15,20.);
    fEPESDtrackCuts->SetEtaRange(-0.8,0.8);

    Int_t fAODfilterbit=GetAODEPTrackFilterBit();

    if(fInputEvent->IsA()==AliESDEvent::Class()){
	maxID=fInputEvent->GetNumberOfTracks();
	tracklist=fEPESDtrackCuts->GetAcceptedTracks((AliESDEvent*)fInputEvent,fUseTPCOnlyTracks);
    }
    if(fInputEvent->IsA()==AliAODEvent::Class()){

        // From AliEPSelectionTask
	tracklist = new TObjArray();

	AliAODTrack *tr = 0;
	Int_t maxid1 = 0;
	Int_t maxidtemp = -1;
	Float_t ptlow = 0;
	Float_t ptup = 0;
	Float_t etalow = 0;
	Float_t etaup = 0;
	fEPESDtrackCuts->GetPtRange(ptlow,ptup);
	fEPESDtrackCuts->GetEtaRange(etalow,etaup);
	Int_t ntpc = fEPESDtrackCuts->GetMinNClusterTPC();

	for (Int_t i = 0; i < fInputEvent->GetNumberOfTracks() ; i++){
	    tr = (AliAODTrack*)fInputEvent->GetTrack(i);
	    maxidtemp = tr->GetID();
	    if(maxidtemp < 0 && fAODfilterbit != 128) continue;// id<0 means filter bit 128
	    if(maxidtemp > -1 && fAODfilterbit == 128) continue;// id>01 means filter bit 1
	    if (fAODfilterbit == 128) maxidtemp = maxidtemp*(-1) - 1;
	    if (maxidtemp > maxid1) maxid1 = maxidtemp;
	    if(tr->TestFilterBit(fAODfilterbit) && tr->Pt() < ptup && tr->Pt() > ptlow && tr->Eta() < etaup && tr->Eta() > etalow && tr->GetTPCNcls() > ntpc){
		tracklist->Add(tr);
	    }
	}
	maxID = maxid1;
    }
    delete fEPESDtrackCuts;
    if(!tracklist)AliError("No tracklist");
    return tracklist;
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskPi0v2::GetTPCEventPlane(){

    if(fEP){
	delete fEP;
	fEP=NULL;
    }
    fEP=new AliEventplane();

    float mQx=0, mQy=0;
    float mQx1=0, mQy1=0, mQx2=0, mQy2=0;
    AliVTrack* track;
    Double_t weight;
    Int_t idtemp = -1;
    int trackcounter1=0, trackcounter2=0;

    Int_t maxID=0;

    TObjArray *tracklist=GetEventPlaneTracks(maxID);

    fEP->GetQContributionXArray()->Set(maxID);
    fEP->GetQContributionYArray()->Set(maxID);
    fEP->GetQContributionXArraysub1()->Set(maxID);
    fEP->GetQContributionYArraysub1()->Set(maxID);
    fEP->GetQContributionXArraysub2()->Set(maxID);
    fEP->GetQContributionYArraysub2()->Set(maxID);

    int nt = tracklist->GetEntries();
    for (int i=0; i<nt; i++){
	weight = 1;
	track = dynamic_cast<AliVTrack*> (tracklist->At(i));
	if (track) {
	    // Fill Eta Distribution
	    hEtaTPCEP->Fill(fCentrality,track->Eta());

	    weight=GetWeight(track);
	    idtemp = track->GetID();
	    // TPC only tracks have negative id ((-1)*IDESD - 1) in AOD
	    if (fIsAOD && (fUseTPCOnlyTracks)) idtemp = idtemp*(-1) - 1;

	    Double_t qx=weight*cos(Double_t(fHarmonic)*track->Phi());
	    Double_t qy=weight*sin(Double_t(fHarmonic)*track->Phi());
	    fEP->GetQContributionXArray()->AddAt(qx,idtemp);
	    fEP->GetQContributionYArray()->AddAt(qy,idtemp);

	    mQx += (qx);
	    mQy += (qy);

	    // This loop splits the track set into 2 random subsets
	    if( trackcounter1 < int(nt/2.) && trackcounter2 < int(nt/2.)){
		float random = fRandomizer->Rndm();
		if(random < .5){
		    mQx1 += (qx);
		    mQy1 += (qy);
		    fEP->GetQContributionXArraysub1()->AddAt(qx,idtemp);
		    fEP->GetQContributionYArraysub1()->AddAt(qy,idtemp);
		    trackcounter1++;
		    }
		else {
		    mQx2 += (qx);
		    mQy2 += (qy);
		    fEP->GetQContributionXArraysub2()->AddAt(qx,idtemp);
		    fEP->GetQContributionYArraysub2()->AddAt(qy,idtemp);
		    trackcounter2++;
		}
	    }
	    else{
		if( trackcounter1 >= int(nt/2.)){
		    mQx2 += (qx);
		    mQy2 += (qy);
		    fEP->GetQContributionXArraysub2()->AddAt(qx,idtemp);
		    fEP->GetQContributionYArraysub2()->AddAt(qy,idtemp);
		    trackcounter2++;
		}
		else {
		    mQx1 += (qx);
		    mQy1 += (qy);
		    fEP->GetQContributionXArraysub1()->AddAt(qx,idtemp);
		    fEP->GetQContributionYArraysub1()->AddAt(qy,idtemp);
		    trackcounter1++;
		}
	    }
	}
    }

    tracklist->Clear();
    delete tracklist;
    tracklist = NULL;

    TVector2 *mQ=new TVector2();
    mQ->Set(mQx,mQy);
    Double_t EPAngle=mQ->Phi()/Double_t(fHarmonic);
   
    TVector2 *fQsub1=new TVector2();
    TVector2 *fQsub2=new TVector2();
    fQsub1->Set(mQx1,mQy1);
    fQsub2->Set(mQx2,mQy2);

    fEP->SetQVector(mQ);
    fEP->SetEventplaneQ(EPAngle);
    fEP->SetQsub(fQsub1,fQsub2);
    fEP->SetQsubRes(fQsub1->Phi()/Double_t(fHarmonic) - fQsub2->Phi()/Double_t(fHarmonic));

    Int_t ntracks=trackcounter1+trackcounter2;

    if(ntracks<3)return kFALSE;// <3 -> no subevents
    return kTRUE;
}


//____________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::ApplyFlattening(Double_t phi,EEventPlane ep){

    Double_t c2=0;
    Double_t s2=0;
    Double_t c4=0;
    Double_t s4=0;

    // GetCorrection Parameter
    if(fHarmonic==2){
     
	if(ep==kEPTPC){

	    Double_t cc2[5]={0.00904396,0.00472483,0.00306154,0.00218462,0.00167447};
	    Double_t cs2[5]={0.00885519,0.00516223,0.00411065,0.00380145,0.00324424};
	    Double_t cc4[5]={-0.00110933,-0.00110521,-0.00124342,0.00104131,0.000651779};
	    Double_t cs4[5]={0.00163869,-0.00053565,0.000878745,-0.000563657,-0.000604021};

	    c2=cc2[fCentralityBin];
	    s2=cs2[fCentralityBin];
	    c4=cc4[fCentralityBin];
	    s4=cs4[fCentralityBin];
	}

	if(ep==kEPTPCEtaA){

	    Double_t cc2[5]={0.00529447,0.00278029,0.00315325,0.00173634,0.000763168};
	    Double_t cs2[5]={0.00314285,0.00170173,0.00263333,0.0018509,0.00223784};
	    Double_t cc4[5]={-0.000737254,-0.00037845,-0.000492715,0.000775897,0.000768656};
	    Double_t cs4[5]={0.000347583,3.79872e-05,0.000387037,-0.000186129,0.000432698};

	    c2=cc2[fCentralityBin];
	    s2=cs2[fCentralityBin];
	    c4=cc4[fCentralityBin];
	    s4=cs4[fCentralityBin];
	}

	if(ep==kEPTPCEtaC){

	    Double_t cc2[5]={-0.00562282,-0.00456735,-0.00306068,-0.0027173,-0.00172432};
	    Double_t cs2[5]={0.0101804,0.00430782,0.00394715,0.00350156,0.00302749};
	    Double_t cc4[5]={0.00150831,-0.00159271,-0.000964157,0.000525894,9.93172e-05};
	    Double_t cs4[5]={0.00119279,-4.74629e-05,0.000118845,0.000278554,3.20868e-05};

	    c2=cc2[fCentralityBin];
	    s2=cs2[fCentralityBin];
	    c4=cc4[fCentralityBin];
	    s4=cs4[fCentralityBin];
	}

	if(ep==kEPV0A){

	    Double_t cc2[5]={0.046427,0.0105401,-0.000152992,-0.00578274,-0.0108038};
	    Double_t cs2[5]={0.00551503,0.0158159,0.00965148,0.00135414,-0.00548846};
	    Double_t cc4[5]={0.00362833,0.00170777,0.000152998,0.00223823,0.00215164};
	    Double_t cs4[5]={0.00349056,0.00142802,0.00123298,0.00207995,0.00145625};

	    c2=cc2[fCentralityBin];
	    s2=cs2[fCentralityBin];
	    c4=cc4[fCentralityBin];
	    s4=cs4[fCentralityBin];
	}

	if(ep==kEPV0C){

	    Double_t cc2[5]={-0.00473277,-0.000371313,0.000857122,-1.54263e-05,-0.000686139};
	    Double_t cs2[5]={0.00408304,-0.00208615,-0.00149018,-0.000853616,-2.78855e-05};
	    Double_t cc4[5]={-0.00451741,-0.00399036,-0.00318784,-0.00186472,-0.00106299};
	    Double_t cs4[5]={0.00188045,-0.00713956,-0.00484254,-0.00448149,-0.00482164};

	    c2=cc2[fCentralityBin];
	    s2=cs2[fCentralityBin];
	    c4=cc4[fCentralityBin];
	    s4=cs4[fCentralityBin];
	}
    }

    if(fHarmonic==3){

	if(ep==kEPTPC){

	    Double_t cc2[5]={0.0116542,0.0103631,0.00897965,0.00707409,0.00605151};
	    Double_t cs2[5]={-0.0171191,-0.013024,-0.0114752,-0.0086613,-0.00706863};
	    Double_t cc4[5]={-0.000602948,0.00144836,-0.000193641,0.000108773,-0.000518333};
	    Double_t cs4[5]={-0.00164769,0.00134327,-0.00106369,7.96546e-06,-0.000261517};

	    c2=cc2[fCentralityBin];
	    s2=cs2[fCentralityBin];
	    c4=cc4[fCentralityBin];
	    s4=cs4[fCentralityBin];
	}

	if(ep==kEPTPCEtaA){

	    Double_t cc2[5]={0.000386277,0.000119225,0.00111969,0.000534801,0.000642703};
	    Double_t cs2[5]={-0.00581604,-0.00607255,-0.00443819,-0.00268834,-0.00299961};
	    Double_t cc4[5]={0.00051635,0.00036326,-0.000221272,4.66775e-05,-3.05784e-06};
	    Double_t cs4[5]={1.43285e-05,0.000514099,0.000619339,0.00106466,0.000344196};

	    c2=cc2[fCentralityBin];
	    s2=cs2[fCentralityBin];
	    c4=cc4[fCentralityBin];
	    s4=cs4[fCentralityBin];
	}

	if(ep==kEPTPCEtaC){

	    Double_t cc2[5]={0.0116475,0.0102385,0.00801121,0.00552336,0.00423273};
	    Double_t cs2[5]={-0.0112722,-0.00796059,-0.00683678,-0.00531097,-0.00430716};
	    Double_t cc4[5]={-0.000609051,1.36573e-08,-0.000464961,-0.000387943,-2.28363e-05};
	    Double_t cs4[5]={0.00125449,0.00168484,-0.000390491,-0.000219447,8.11997e-07};

	    c2=cc2[fCentralityBin];
	    s2=cs2[fCentralityBin];
	    c4=cc4[fCentralityBin];
	    s4=cs4[fCentralityBin];
	}

	if(ep==kEPV0A){

            Double_t cc2[5]={-0.0057427,-0.00482728,-0.00565919,-0.000717094,-0.00933233};
	    Double_t cs2[5]={0.0306554,-0.0144675,-0.0159243,-0.0120465,-0.00814124};
	    Double_t cc4[5]={-0.002868,0.00159533,0.00754171,0.00683898,0.00689441};
	    Double_t cs4[5]={0.00083196,0.00198133,4.68307e-05,-0.00018187,-0.0014258};

	    c2=cc2[fCentralityBin];
	    s2=cs2[fCentralityBin];
	    c4=cc4[fCentralityBin];
	    s4=cs4[fCentralityBin];
	}

	if(ep==kEPV0C){

	    Double_t cc2[5]={-0.00259141,-0.00115826,-0.000738658,-4.96667e-05,-0.000346694};
	    Double_t cs2[5]={-0.0111001,0.00258109,0.00110959,-0.000147296,-0.000199817};
	    Double_t cc4[5]={0.000968742,0.00157903,0.000206157,0.000444206,-0.00046573};
	    Double_t cs4[5]={-0.00307319,-0.0047952,-0.00412117,-0.00320344,-0.00386629};

	    c2=cc2[fCentralityBin];
	    s2=cs2[fCentralityBin];
	    c4=cc4[fCentralityBin];
	    s4=cs4[fCentralityBin];
	}
    }

    // Do correction
    Double_t newphi=phi;
    newphi+=1/Double_t(fHarmonic)*(2*c2*sin(Double_t(fHarmonic)*phi)-2*s2*cos(Double_t(fHarmonic)*phi)+c4*sin(2*Double_t(fHarmonic)*phi)-s4*cos(2*Double_t(fHarmonic)*phi));
    newphi=GetPsiInRange(newphi);

    return newphi;
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetWeight(TObject* track1)
{
    Double_t ptweight=1;
    AliVTrack* track = dynamic_cast<AliVTrack*>(track1);
    if (track) {
	if (track->Pt()<2) ptweight=track->Pt();
	else ptweight=2;
    }
    return ptweight*GetPhiWeight(track);
}

//________________________________________________________________________
Double_t AliAnalysisTaskPi0v2::GetPhiWeight(TObject* track1)
{
  Double_t phiweight=1;
  AliVTrack* track = dynamic_cast<AliVTrack*>(track1);

  TH1F *phiDist = 0x0;
  if(track) phiDist = SelectPhiDist(track);

  if (phiDist && track) {
      Double_t nParticles = phiDist->Integral();
      Double_t nPhibins = phiDist->GetNbinsX();

      Double_t Phi = track->Phi();

      while (Phi<0) Phi += TMath::TwoPi();
      while (Phi>TMath::TwoPi()) Phi -= TMath::TwoPi();

      Double_t PhiDistValue = phiDist->GetBinContent(1+TMath::FloorNint((track->Phi())*nPhibins/TMath::TwoPi()));

      if (PhiDistValue > 0) phiweight = nParticles/nPhibins/PhiDistValue;
  }
  return phiweight;
}

//_________________________________________________________________________
TH1F* AliAnalysisTaskPi0v2::SelectPhiDist(AliVTrack *track)
{ 
  if (fPeriod.CompareTo("LHC10h")==0) return fPhiDist[0];
  else if(fPeriod.CompareTo("LHC11h")==0)
    {
     if (track->Charge() < 0)
       {
        if(track->Eta() < 0.)       return fPhiDist[0];
        else if (track->Eta() > 0.) return fPhiDist[2];
       }
      else if (track->Charge() > 0)
       {
        if(track->Eta() < 0.)       return fPhiDist[1];
        else if (track->Eta() > 0.) return fPhiDist[3];
       }
       
    }
  return 0;
}

