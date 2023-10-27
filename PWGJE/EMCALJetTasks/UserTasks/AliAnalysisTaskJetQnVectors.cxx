/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************************************
// \class AliAnalysisTaskJetQnVectors
// \brief task used to load the Qn calibrations and get the calibrated Qn vectors for JE analyses
// \adapted from HF task
// \authors of this task
// C. Beattie, caitie.beattie@yale.edu
// M. Sas, mike.sas@cern.ch
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <TChain.h>

#include "AliAnalysisTaskJetQnVectors.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliAODVertex.h"
#include "AliDataFile.h"
#include "TGrid.h"

using namespace std;

ClassImp(AliAnalysisTaskJetQnVectors)

//________________________________________________________________________
AliAnalysisTaskJetQnVectors::AliAnalysisTaskJetQnVectors() :
    AliAnalysisTaskSE("JEQnVectors"),
    fOutputList(nullptr),
    fHistNEvents(nullptr),
    fHistCentrality(nullptr),
    fHistResolution_epV0AV0C_qV0M(nullptr),
    fHistResolution_epV0CTPC_qV0M(nullptr),
    fHistResolution_epV0ATPC_qV0M(nullptr),
    fHistResolution_epV0AV0C_qV0A(nullptr),
    fHistResolution_epV0CTPC_qV0A(nullptr),
    fHistResolution_epV0ATPC_qV0A(nullptr),
    fHistResolution_epV0AV0C_qV0C(nullptr),
    fHistResolution_epV0CTPC_qV0C(nullptr),
    fHistResolution_epV0ATPC_qV0C(nullptr),
    fHistResolution_epTPCpTPCn_qV0M(nullptr),
    fHistResolution_epV0MTPCp_qV0M(nullptr),
    fHistResolution_epV0MTPCn_qV0M(nullptr),
    fEnableTPCPhiVsCentrDistr(false),
    fEnableQvecTPCVsCentrDistr(false),
    fJEQnVecHandler1(nullptr),
    fJEQnVecHandler2(nullptr),
    fHarmonic(2),
    fCalibType(AliJEQnVectorHandler::kQnCalib),
    fNormMethod(AliJEQnVectorHandler::kQoverM),
    fOADBFileName1(""),
    fOADBFileName2(""),
    fAOD(nullptr),
    fPrevEventRun(-1),
    fTriggerClass(""),
    fTriggerMask(AliVEvent::kAny),
    fRejectTPCPileup(false),
    fq2V0M(0),
    fq2V0A(0),
    fq2V0C(0),
    fq2TPC(0),
    fEPangleFullTPC(0),
    fEPanglePosTPC(0),
    fEPangleNegTPC(0),
    fEPangleV0M(0),
    fEPangleV0A(0),
    fEPangleV0C(0),
    fEventCuts("")
{
    //
    // default constructor
    //
    for(int iDet=0; iDet<3; iDet++) {
        fHistEventPlaneTPC[iDet] = nullptr;
        fHistEventPlaneV0[iDet] = nullptr;
		fHistqnVsCentrTPC[iDet] = nullptr;
		fHistqnVsCentrV0[iDet] = nullptr;
        fSplineListqnPercTPC[iDet] = nullptr;
        fSplineListqnPercV0[iDet] = nullptr;
        fQvecTPCVsCentrDistr[iDet] = nullptr;
        if(iDet<2)
            fTPCPhiVsCentrDistr[iDet] = nullptr;
    }
}

//________________________________________________________________________
AliAnalysisTaskJetQnVectors::AliAnalysisTaskJetQnVectors(const char *name, int harmonic, int calibType, TString oadbFileName1, TString oadbFileName2) :
    AliAnalysisTaskSE(name),
    fOutputList(nullptr),
    fHistNEvents(nullptr),
    fHistCentrality(nullptr),
    fHistResolution_epV0AV0C_qV0M(nullptr),
    fHistResolution_epV0CTPC_qV0M(nullptr),
    fHistResolution_epV0ATPC_qV0M(nullptr),
    fHistResolution_epV0AV0C_qV0A(nullptr),
    fHistResolution_epV0CTPC_qV0A(nullptr),
    fHistResolution_epV0ATPC_qV0A(nullptr),
    fHistResolution_epV0AV0C_qV0C(nullptr),
    fHistResolution_epV0CTPC_qV0C(nullptr),
    fHistResolution_epV0ATPC_qV0C(nullptr),
    fHistResolution_epTPCpTPCn_qV0M(nullptr),
    fHistResolution_epV0MTPCp_qV0M(nullptr),
    fHistResolution_epV0MTPCn_qV0M(nullptr),
    fEnableTPCPhiVsCentrDistr(false),
    fEnableQvecTPCVsCentrDistr(false),
    fJEQnVecHandler1(nullptr),
    fJEQnVecHandler2(nullptr),
    fHarmonic(harmonic),
    fCalibType(calibType),
    fNormMethod(AliJEQnVectorHandler::kQoverSqrtM),
    fOADBFileName1(oadbFileName1),
    fOADBFileName2(oadbFileName2),
    fAOD(nullptr),
    fPrevEventRun(-1),
    fTriggerClass(""),
    fTriggerMask(AliVEvent::kAny),
    fRejectTPCPileup(false),
    fq2V0M(0),
    fq2V0A(0),
    fq2V0C(0),
    fq2TPC(0),
    fEPangleFullTPC(0),
    fEPanglePosTPC(0),
    fEPangleNegTPC(0),
    fEPangleV0M(0),
    fEPangleV0A(0),
    fEPangleV0C(0),
    fEventCuts("")
{
    //
    // standard constructor
    //
    for(int iDet=0; iDet<3; iDet++) {
        fHistEventPlaneTPC[iDet] = nullptr;
        fHistEventPlaneV0[iDet] = nullptr;
		fHistqnVsCentrTPC[iDet] = nullptr;
		fHistqnVsCentrV0[iDet] = nullptr;
        fSplineListqnPercTPC[iDet] = nullptr;
        fSplineListqnPercV0[iDet] = nullptr;
        fQvecTPCVsCentrDistr[iDet] = nullptr;
        if(iDet<2)
            fTPCPhiVsCentrDistr[iDet] = nullptr;
    }

    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TH2F::Class());
    DefineOutput(3, TH2F::Class());
    DefineOutput(4, TH3F::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetQnVectors::~AliAnalysisTaskJetQnVectors()
{
    //
    // standard destructor
    //
    if(fOutputList) {
		if(!fOutputList->IsOwner()) {
			delete fHistNEvents;
			delete fHistCentrality;
                        delete fHistResolution_epV0AV0C_qV0M;
                        delete fHistResolution_epV0CTPC_qV0M;
                        delete fHistResolution_epV0ATPC_qV0M;
                        delete fHistResolution_epV0AV0C_qV0A;
                        delete fHistResolution_epV0CTPC_qV0A;
                        delete fHistResolution_epV0ATPC_qV0A;
                        delete fHistResolution_epV0AV0C_qV0C;
                        delete fHistResolution_epV0CTPC_qV0C;
                        delete fHistResolution_epV0ATPC_qV0C;
                        delete fHistResolution_epTPCpTPCn_qV0M;
                        delete fHistResolution_epV0MTPCp_qV0M;
                        delete fHistResolution_epV0MTPCn_qV0M;

			for(int iDet=0; iDet<3; iDet++) {
				delete fHistEventPlaneTPC[iDet];
                                delete fHistEventPlaneV0[iDet];
				delete fHistqnVsCentrTPC[iDet];
      				delete fHistqnVsCentrV0[iDet];
				if(fSplineListqnPercTPC[iDet]) delete fSplineListqnPercTPC[iDet];
				if(fSplineListqnPercV0[iDet]) delete fSplineListqnPercV0[iDet];
			}
		}
		delete fOutputList;
    }
    if(fJEQnVecHandler1) delete fJEQnVecHandler1;
    if(fJEQnVecHandler2) delete fJEQnVecHandler2;
    for(int iDet=0; iDet<3; iDet++) {
        if(fQvecTPCVsCentrDistr[iDet]) delete fQvecTPCVsCentrDistr[iDet];
        if(iDet<2) {
            if(fTPCPhiVsCentrDistr[iDet]) delete fTPCPhiVsCentrDistr[iDet];
        }
    }
}

//________________________________________________________________________
void AliAnalysisTaskJetQnVectors::UserCreateOutputObjects()
{
    
    fOutputList = new TList();
    fOutputList->SetOwner(true);

    fHistNEvents = new TH1F("fHistNEvents","Number of processed events;;Number of events",5,0.5,5.5);
    fHistNEvents->Sumw2();
    fHistNEvents->SetMinimum(0);
    fHistNEvents->GetXaxis()->SetBinLabel(1,"Read from AOD");
    fHistNEvents->GetXaxis()->SetBinLabel(2,"Pass Phys. Sel. + Trig");
    fHistNEvents->GetXaxis()->SetBinLabel(3,"Centrality > 90%");
    fHistNEvents->GetXaxis()->SetBinLabel(4,"No vertex");
    fHistNEvents->GetXaxis()->SetBinLabel(5,"Fail Pileup");
    fOutputList->Add(fHistNEvents);

    fHistCentrality = new TH1F("fHistCentrality","Centrality distribution;;Number of events",100,0.,100.);
    fOutputList->Add(fHistCentrality);
   
    //V0M Resolution Hists
    fHistResolution_epV0AV0C_qV0M = new TH3F("fHistResolution_epV0AV0C_qV0M", "fHistResolution_epV0AV0C_qV0M", 200,-1.,1., 200,0.,20., 100,0.,100.);
    fOutputList->Add(fHistResolution_epV0AV0C_qV0M);
    fHistResolution_epV0CTPC_qV0M = new TH3F("fHistResolution_epV0CTPC_qV0M", "fHistResolution_epV0CTPC_qV0M", 200,-1.,1., 200,0.,20., 100,0.,100.);
    fOutputList->Add(fHistResolution_epV0CTPC_qV0M);
    fHistResolution_epV0ATPC_qV0M = new TH3F("fHistResolution_epV0ATPC_qV0M", "fHistResolution_epV0ATPC_qV0M", 200,-1.,1., 200,0.,20., 100,0.,100.);
    fOutputList->Add(fHistResolution_epV0ATPC_qV0M);
    fHistResolution_epTPCpTPCn_qV0M = new TH3F("fHistResolution_epTPCpTPCn_qV0M", "fHistResolution_epTPCpTPCn_qV0M", 200,-1.,1., 200,0.,20., 100,0.,100.);
    fOutputList->Add(fHistResolution_epTPCpTPCn_qV0M);
    fHistResolution_epV0MTPCp_qV0M = new TH3F("fHistResolution_epV0MTPCp_qV0M", "fHistResolution_epV0MTPCn_qV0M", 200,-1.,1., 200,0.,20., 100,0.,100.);
    fOutputList->Add(fHistResolution_epV0MTPCp_qV0M);
    fHistResolution_epV0MTPCn_qV0M = new TH3F("fHistResolution_epV0MTPCn_qV0M", "fHistResolution_epV0MTPCn_qV0M", 200,-1.,1., 200,0.,20., 100,0.,100.);
    fOutputList->Add(fHistResolution_epV0MTPCn_qV0M);

    //V0A Resolution Hists
    fHistResolution_epV0AV0C_qV0A = new TH3F("fHistResolution_epV0AV0C_qV0A", "fHistResolution_epV0AV0C_qV0A", 200,-1.,1., 200,0.,20., 100,0.,100.);
    fOutputList->Add(fHistResolution_epV0AV0C_qV0A);
    fHistResolution_epV0CTPC_qV0A = new TH3F("fHistResolution_epV0CTPC_qV0A", "fHistResolution_epV0CTPC_qV0A", 200,-1.,1., 200,0.,20., 100,0.,100.);
    fOutputList->Add(fHistResolution_epV0CTPC_qV0A);
    fHistResolution_epV0ATPC_qV0A = new TH3F("fHistResolution_epV0ATPC_qV0A", "fHistResolution_epV0ATPC_qV0A", 200,-1.,1., 200,0.,20., 100,0.,100.);
    fOutputList->Add(fHistResolution_epV0ATPC_qV0A);

    //V0C Resolution Hists
    fHistResolution_epV0AV0C_qV0C = new TH3F("fHistResolution_epV0AV0C_qV0C", "fHistResolution_epV0AV0C_qV0C", 200,-1.,1., 200,0.,20., 100,0.,100.);
    fOutputList->Add(fHistResolution_epV0AV0C_qV0C);
    fHistResolution_epV0CTPC_qV0C = new TH3F("fHistResolution_epV0CTPC_qV0C", "fHistResolution_epV0CTPC_qV0C", 200,-1.,1., 200,0.,20., 100,0.,100.);
    fOutputList->Add(fHistResolution_epV0CTPC_qV0C);
    fHistResolution_epV0ATPC_qV0C = new TH3F("fHistResolution_epV0ATPC_qV0C", "fHistResolution_epV0ATPC_qV0C", 200,-1.,1., 200,0.,20., 100,0.,100.);
    fOutputList->Add(fHistResolution_epV0ATPC_qV0C);

    TString TPCNames[3] = {"FullTPC","PosTPC","NegTPC"};
    TString V0Names[3] = {"FullV0","V0A","V0C"};
    for(int iDet=0; iDet<3; iDet++) {
        fHistEventPlaneTPC[iDet] = new TH1F(Form("fHistEventPlane%s",TPCNames[iDet].Data()),Form(";%s #psi_{%d};Entries",TPCNames[iDet].Data(),fHarmonic),200,0,TMath::Pi());
        fHistEventPlaneV0[iDet] = new TH1F(Form("fHistEventPlane%s",V0Names[iDet].Data()),Form(";%s #psi_{%d};Entries",V0Names[iDet].Data(),fHarmonic),200,0,TMath::Pi());
    	fHistqnVsCentrTPC[iDet] = new TH2F(Form("fHistqnVsCentr%s",TPCNames[iDet].Data()),Form(";#it{q}_{%d}^{%s};Entries",fHarmonic,TPCNames[iDet].Data()),100,0.,100.,2000,0.,20.);
		fHistqnVsCentrV0[iDet] = new TH2F(Form("fHistqnVsCentr%s",V0Names[iDet].Data()),Form(";#it{q}_{%d}^{%s};Entries",fHarmonic,V0Names[iDet].Data()),100,0.,100.,2000,0.,20.);
        fOutputList->Add(fHistEventPlaneTPC[iDet]);
        fOutputList->Add(fHistEventPlaneV0[iDet]);
        fOutputList->Add(fHistqnVsCentrTPC[iDet]);
        fOutputList->Add(fHistqnVsCentrV0[iDet]);
	}

    if(fCalibType==AliJEQnVectorHandler::kQnCalib) {
        if(fEnableTPCPhiVsCentrDistr) {
            OpenFile(2);
            fTPCPhiVsCentrDistr[0] = dynamic_cast<TH2F*>(fJEQnVecHandler1->GetPhiDistrHistosTPCPosEta()->Clone());
            fTPCPhiVsCentrDistr[0]->Add(dynamic_cast<TH2F*>(fJEQnVecHandler2->GetPhiDistrHistosTPCPosEta()->Clone())); 
            OpenFile(3);
            fTPCPhiVsCentrDistr[1] = dynamic_cast<TH2F*>(fJEQnVecHandler1->GetPhiDistrHistosTPCNegEta()->Clone());
            fTPCPhiVsCentrDistr[1]->Add(dynamic_cast<TH2F*>(fJEQnVecHandler2->GetPhiDistrHistosTPCNegEta()->Clone()));
        }
        if(fEnableQvecTPCVsCentrDistr) {
            OpenFile(4);
            fQvecTPCVsCentrDistr[0] = new TH2F("fQvecTPCVsCentrDistr",Form(";centrality (%%);#it{#Q_{%d}^{TPC}}",fHarmonic),10,0,100,100,0,5000);
            OpenFile(5);
            fQvecTPCVsCentrDistr[1] = new TH2F("fQvecTPCVsCentrDistrPosEta",Form(";centrality (%%);#it{#Q_{%d}^{TPCPosEta}}",fHarmonic),10,0,100,100,0,5000);
            OpenFile(6);
            fQvecTPCVsCentrDistr[2] = new TH2F("fQvecTPCVsCentrDistrNegEta",Form(";centrality (%%);#it{#Q_{%d}^{TPCNegEta}}",fHarmonic),10,0,100,100,0,5000);
        }
    }

    // post data
    PostData(1, fOutputList);
    if(fCalibType==AliJEQnVectorHandler::kQnCalib) {
        if(fEnableTPCPhiVsCentrDistr) {
            PostData(2, fTPCPhiVsCentrDistr[0]);
            PostData(3, fTPCPhiVsCentrDistr[1]);
        }
        if(fEnableQvecTPCVsCentrDistr) {
            PostData(4, fQvecTPCVsCentrDistr[0]);
            PostData(5, fQvecTPCVsCentrDistr[1]);
            PostData(6, fQvecTPCVsCentrDistr[2]);
        }
    }
}

//________________________________________________________________________
void AliAnalysisTaskJetQnVectors::UserExec(Option_t */*option*/)
{
    // main event loop
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAOD && AODEvent() && IsStandardAOD()) {
        // In case there is an AOD handler writing a standard AOD, use the AOD
        // event in memory rather than the input (ESD) event.
        fAOD = dynamic_cast<AliAODEvent*> (AODEvent());
    }
    if (!fAOD) {
        AliWarning("AliAnalysisTaskJetQnVectors::Exec(): bad AOD");
        return;
    }
    fHistNEvents->Fill(1);

    if(TMath::Abs(fAOD->GetMagneticField())<0.001) return;

    AliAODHandler* aodHandler = static_cast<AliAODHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    if(!aodHandler) {
        AliWarning("AliAnalysisTaskJetQnVectors::Exec(): No AliInputEventHandler!");
        return;
    }

    unsigned int maskPhysSel = aodHandler->IsEventSelected();
    TString firedTriggerClasses = fAOD->GetFiredTriggerClasses();
    if((fAOD->GetRunNumber()<136851 || fAOD->GetRunNumber()>139517)) {
        if(!(firedTriggerClasses.Contains(fTriggerClass.Data()))) return;
    }
    if((maskPhysSel & fTriggerMask)==0.) {
        return;
    }


    AliMultSelection *multSelection = dynamic_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
    if(!multSelection){
        AliWarning("AliMultSelection could not be found in the aod event list of objects");
        return;
    }
    fHistNEvents->Fill(2);

    float cent = multSelection->GetMultiplicityPercentile("V0M");
    if(cent>90) {
        fHistNEvents->Fill(3);
        return;
    }

    const AliAODVertex *vertex = fAOD->GetPrimaryVertex();
    if(!vertex || TMath::Abs(vertex->GetZ())>10. || vertex->GetNContributors()<=0) {
        fHistNEvents->Fill(4);
        return;
    }
    if (fRejectTPCPileup)   {
        fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,1);
        Bool_t acceptEventCuts = fEventCuts.AcceptEvent(InputEvent());
        if(!acceptEventCuts)   { 
           fHistNEvents->Fill(5);
            return;}
    }

    fHistCentrality->Fill(cent);

    //choose calibration file to use (run numbers specific to 2018 pass3)
    AliJEQnVectorHandler *fJEQnVecHandler;
    if (fAOD->GetRunNumber() <= 296623) fJEQnVecHandler = fJEQnVecHandler1;   //child1
    if (fAOD->GetRunNumber() >  296623) fJEQnVecHandler = fJEQnVecHandler2;   //child2

    fJEQnVecHandler->ResetAODEvent();
    fJEQnVecHandler->SetAODEvent(fAOD);
    fJEQnVecHandler->ComputeCalibratedQnVectorTPC();
    fJEQnVecHandler->ComputeCalibratedQnVectorV0();

    //fill histos with EP angle
    double PsinFullTPC = -1., PsinPosTPC = -1., PsinNegTPC = -1.;
    double PsinFullV0 = -1., PsinV0A = -1., PsinV0C = -1.;
    fJEQnVecHandler->GetEventPlaneAngleTPC(PsinFullTPC,PsinPosTPC,PsinNegTPC);
    fJEQnVecHandler->GetEventPlaneAngleV0(PsinFullV0,PsinV0A,PsinV0C);

    fEPangleFullTPC = PsinFullTPC;
    fEPanglePosTPC = PsinPosTPC;
    fEPangleNegTPC = PsinNegTPC;

    fEPangleV0M = PsinFullV0;
    fEPangleV0A = PsinV0A;
    fEPangleV0C = PsinV0C;

    fHistEventPlaneTPC[0]->Fill(PsinFullTPC);
    fHistEventPlaneTPC[1]->Fill(PsinPosTPC);
    fHistEventPlaneTPC[2]->Fill(PsinNegTPC);
    fHistEventPlaneV0[0]->Fill(PsinFullV0);
    fHistEventPlaneV0[1]->Fill(PsinV0A);
    fHistEventPlaneV0[2]->Fill(PsinV0C);

    //fill histos for q2 spline calibration
    double qnFullTPC = -1., qnPosTPC = -1., qnNegTPC = -1.;
    double qnFullV0 = -1., qnV0A = -1., qnV0C = -1.;
    fJEQnVecHandler->GetqnTPC(qnFullTPC,qnPosTPC,qnNegTPC);
    fJEQnVecHandler->GetqnV0(qnFullV0,qnV0A,qnV0C);
    fq2V0M = qnFullV0;
    fq2V0A = qnV0A;
    fq2V0C = qnV0C;
    fq2TPC = qnFullTPC;

	fHistqnVsCentrTPC[0]->Fill(cent,qnFullTPC);
	fHistqnVsCentrTPC[1]->Fill(cent,qnPosTPC);
	fHistqnVsCentrTPC[2]->Fill(cent,qnNegTPC);
	fHistqnVsCentrV0[0]->Fill(cent,qnFullV0);
	fHistqnVsCentrV0[1]->Fill(cent,qnV0A);
	fHistqnVsCentrV0[2]->Fill(cent,qnV0C);

    //fill histos with resolution
    double resolution1 = cos(2*(PsinV0C-PsinV0A));
    double resolution2 = cos(2*(PsinV0C-PsinFullTPC));
    double resolution3 = cos(2*(PsinV0A-PsinFullTPC));
    double resolution4 = cos(2*(PsinPosTPC-PsinNegTPC));
    double resolution5 = cos(2*(PsinFullV0-PsinPosTPC));
    double resolution6 = cos(2*(PsinFullV0-PsinNegTPC));

    fHistResolution_epV0AV0C_qV0M->Fill(resolution1, qnFullV0, cent);  //used to calculate reaction plane resolution
    fHistResolution_epV0CTPC_qV0M->Fill(resolution2, qnFullV0, cent);  //2nd argument gives detector being used for q2 calculation
    fHistResolution_epV0ATPC_qV0M->Fill(resolution3, qnFullV0, cent);
    fHistResolution_epTPCpTPCn_qV0M->Fill(resolution4, qnFullV0, cent);
    fHistResolution_epV0MTPCp_qV0M->Fill(resolution5, qnFullV0, cent);
    fHistResolution_epV0MTPCn_qV0M->Fill(resolution6, qnFullV0, cent);

    fHistResolution_epV0AV0C_qV0A->Fill(resolution1, qnV0A, cent); 
    fHistResolution_epV0CTPC_qV0A->Fill(resolution2, qnV0A, cent); 
    fHistResolution_epV0ATPC_qV0A->Fill(resolution3, qnV0A, cent);

    fHistResolution_epV0AV0C_qV0C->Fill(resolution1, qnV0C, cent);  
    fHistResolution_epV0CTPC_qV0C->Fill(resolution2, qnV0C, cent);  
    fHistResolution_epV0ATPC_qV0C->Fill(resolution3, qnV0C, cent);

    int runnumber = fAOD->GetRunNumber();

    if(fCalibType==AliJEQnVectorHandler::kQnCalib) {
        if(fEnableTPCPhiVsCentrDistr) {
            TH2F* hPhiPosEta = fJEQnVecHandler->GetPhiDistrHistosTPCPosEta();
            TH2F* hPhiNegEta = fJEQnVecHandler->GetPhiDistrHistosTPCNegEta();
            if(runnumber!=fPrevEventRun) {
                fTPCPhiVsCentrDistr[0]->SetName(Form("%s_%d",hPhiPosEta->GetName(),runnumber));
                fTPCPhiVsCentrDistr[1]->SetName(Form("%s_%d",hPhiNegEta->GetName(),runnumber));
            }

            if(hPhiPosEta)
                fTPCPhiVsCentrDistr[0]->Add((TH2F*)hPhiPosEta->Clone());
            if(hPhiNegEta)
                fTPCPhiVsCentrDistr[1]->Add((TH2F*)hPhiNegEta->Clone());

            PostData(2, fTPCPhiVsCentrDistr[0]);
            PostData(3, fTPCPhiVsCentrDistr[1]);
        }
        if(fEnableQvecTPCVsCentrDistr) {
            if(runnumber!=fPrevEventRun) {
                for(int iDet=0; iDet<3; iDet++)
                    fQvecTPCVsCentrDistr[iDet]->SetName(Form("%s_%d",fQvecTPCVsCentrDistr[iDet]->GetName(),runnumber));
            }

            double QnVecFullTPC[2], QnVecPosTPC[2], QnVecNegTPC[2];
            fJEQnVecHandler->GetUnNormQnVecTPC(QnVecFullTPC, QnVecPosTPC, QnVecNegTPC);
            fQvecTPCVsCentrDistr[0]->Fill(cent, TMath::Sqrt(QnVecFullTPC[0]*QnVecFullTPC[0]+QnVecFullTPC[1]*QnVecFullTPC[1]));
            fQvecTPCVsCentrDistr[1]->Fill(cent, TMath::Sqrt(QnVecPosTPC[0]*QnVecPosTPC[0]+QnVecPosTPC[1]*QnVecPosTPC[1]));
            fQvecTPCVsCentrDistr[2]->Fill(cent, TMath::Sqrt(QnVecNegTPC[0]*QnVecNegTPC[0]+QnVecNegTPC[1]*QnVecNegTPC[1]));

            PostData(4, fQvecTPCVsCentrDistr[0]);
            PostData(5, fQvecTPCVsCentrDistr[1]);
            PostData(6, fQvecTPCVsCentrDistr[2]);
        }
    }
    fPrevEventRun = runnumber;

    // Post output data
    PostData(1, fOutputList);
}

//
// Creates the QnVector Handlers. Needs to be run in add task so that calibration
// files are loaded at job creation.
//________________________________________________________________________
void AliAnalysisTaskJetQnVectors::CreateQnVectorHandlers() {

    fJEQnVecHandler1 = new AliJEQnVectorHandler(fCalibType,fNormMethod,fHarmonic,fOADBFileName1);
    fJEQnVecHandler2 = new AliJEQnVectorHandler(fCalibType,fNormMethod,fHarmonic,fOADBFileName2);
    fJEQnVecHandler1->EnablePhiDistrHistos();
    fJEQnVecHandler2->EnablePhiDistrHistos();

}

//________________________________________________________________________
TDirectoryFile* AliAnalysisTaskJetQnVectors::GetSplineForqnPercentileList(int det) const
{
    if(det<=kNegTPC) {
        return fSplineListqnPercTPC[det];
    }
    else if(det>=kFullV0 && det<=kV0C) {
        return fSplineListqnPercV0[det-3];
    }
    else {
        AliWarning("Spline List not found!");
        return nullptr;
    }
}

//________________________________________________________________________
void AliAnalysisTaskJetQnVectors::LoadSplinesForqnPercentile(TString splinesFilePath)
{
    // load splines for qn percentiles

    TString listnameTPC[3] = {"SplineListq2TPC", "SplineListq2TPCPosEta", "SplineListq2TPCNegEta"};
    TString listnameV0[3] = {"SplineListq2V0", "SplineListq2V0A", "SplineListq2V0C"};

    TString pathToFileCMVFNS = AliDataFile::GetFileName(splinesFilePath.Data());
    TString pathToFileLocal = splinesFilePath;

      TFile* splinesfile;
      if (splinesFilePath.BeginsWith("alien://") && !gGrid)
        {
          AliInfo("Trying to connect to AliEn ...");
          TGrid::Connect("alien://");
        }
      // Check access to CVMFS (will only be displayed locally)
      if (!pathToFileCMVFNS.IsNull())   splinesfile = TFile::Open(pathToFileCMVFNS.Data());
      if (pathToFileCMVFNS.IsNull())    splinesfile = TFile::Open(splinesFilePath.Data());
   
      if(!splinesfile)    AliFatal("File with splines for qn percentiles not found!");

    for(int iDet=0; iDet<3; iDet++) {
        fSplineListqnPercTPC[iDet] = (TDirectoryFile*)splinesfile->Get(listnameTPC[iDet].Data());
        if(!fSplineListqnPercTPC[iDet]) {
           std::cout << "Warning: TDirectoryFile with splines for qnTPC percentiles not found in the spline file!\n";
           continue;}
        fSplineListqnPercV0[iDet] = (TDirectoryFile*)splinesfile->Get(listnameV0[iDet].Data());
        if(!fSplineListqnPercV0[iDet])  { 
           std::cout << "Warning: TDirectoryFile with splines for qnV0 percentiles not found in the spline file!\n";
           continue;}

        //uncoment if using TLists
        /*fSplineListqnPercTPC[iDet]->SetOwner();
        fSplineListqnPercV0[iDet]->SetOwner();*/
    }
    splinesfile->Close();
}
