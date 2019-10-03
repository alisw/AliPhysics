/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************************************
// \class AliAnalysisTaskSEHFTenderQnVectors
// \brief task used to load the Qn calibrations and get the calibrated Qn vectors for HF analyses
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
// F. Catalano, fabio.catalano@cern.ch
// A. Dobrin, alexandru.dobrin@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// G. Luparello, grazia.luparello@cern.ch
// F. Prino, prino@to.infn.it
// A. Rossi, andrea.rossi@cern.ch
// S. Trogolo, stefano.trogolo@cern.ch
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <TChain.h>

#include "AliAnalysisTaskSEHFTenderQnVectors.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliAODVertex.h"

ClassImp(AliAnalysisTaskSEHFTenderQnVectors)

//________________________________________________________________________
AliAnalysisTaskSEHFTenderQnVectors::AliAnalysisTaskSEHFTenderQnVectors() :
    AliAnalysisTaskSE("HFTenderQnVectors"),
    fOutputList(nullptr),
    fHistNEvents(nullptr),
    fHistCentrality(nullptr),
    fEnableTPCPhiVsCentrDistr(false),
    fEnableQvecTPCVsCentrDistr(false),
    fHFQnVecHandler(nullptr),
    fHarmonic(2),
    fCalibType(AliHFQnVectorHandler::kQnCalib),
    fNormMethod(AliHFQnVectorHandler::kQoverM),
    fOADBFileName(""),
    fAOD(nullptr),
    fPrevEventRun(-1),
    fTriggerClass(""),
    fTriggerMask(AliVEvent::kAny)
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
AliAnalysisTaskSEHFTenderQnVectors::AliAnalysisTaskSEHFTenderQnVectors(const char *name, int harmonic, int calibType, TString oadbFileName) :
    AliAnalysisTaskSE(name),
    fOutputList(nullptr),
    fHistNEvents(nullptr),
    fHistCentrality(nullptr),
    fEnableTPCPhiVsCentrDistr(false),
    fEnableQvecTPCVsCentrDistr(false),
    fHFQnVecHandler(nullptr),
    fHarmonic(harmonic),
    fCalibType(calibType),
    fNormMethod(AliHFQnVectorHandler::kQoverSqrtM),
    fOADBFileName(oadbFileName),
    fAOD(nullptr),
    fPrevEventRun(-1),
    fTriggerClass(""),
    fTriggerMask(AliVEvent::kAny)
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
}

//________________________________________________________________________
AliAnalysisTaskSEHFTenderQnVectors::~AliAnalysisTaskSEHFTenderQnVectors()
{
    //
    // standard destructor
    //
    if(fOutputList) {
		if(!fOutputList->IsOwner()) {
			delete fHistNEvents;
			delete fHistCentrality;
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
    if(fHFQnVecHandler) delete fHFQnVecHandler;
    for(int iDet=0; iDet<3; iDet++) {
        if(fQvecTPCVsCentrDistr[iDet]) delete fQvecTPCVsCentrDistr[iDet];
        if(iDet<2) {
            if(fTPCPhiVsCentrDistr[iDet]) delete fTPCPhiVsCentrDistr[iDet];
        }
    }
}

//________________________________________________________________________
void AliAnalysisTaskSEHFTenderQnVectors::UserCreateOutputObjects()
{
    fHFQnVecHandler = new AliHFQnVectorHandler(fCalibType,fNormMethod,fHarmonic,fOADBFileName);
    fHFQnVecHandler->EnablePhiDistrHistos();

    fOutputList = new TList();
    fOutputList->SetOwner(true);

    fHistNEvents = new TH1F("fHistNEvents","Number of processed events;;Number of events",4,0.5,4.5);
    fHistNEvents->Sumw2();
    fHistNEvents->SetMinimum(0);
    fHistNEvents->GetXaxis()->SetBinLabel(1,"Read from AOD");
    fHistNEvents->GetXaxis()->SetBinLabel(2,"Pass Phys. Sel. + Trig");
    fHistNEvents->GetXaxis()->SetBinLabel(3,"Centrality > 90%");
    fHistNEvents->GetXaxis()->SetBinLabel(4,"No vertex");
    fOutputList->Add(fHistNEvents);

    fHistCentrality = new TH1F("fHistCentrality","Centrality distribution;;Number of events",100,0.,100.);
    fOutputList->Add(fHistCentrality);

    TString TPCNames[3] = {"FullTPC","PosTPC","NegTPC"};
    TString V0Names[3] = {"FullV0","V0A","V0C"};
    for(int iDet=0; iDet<3; iDet++) {
        fHistEventPlaneTPC[iDet] = new TH1F(Form("fHistEventPlane%s",TPCNames[iDet].Data()),Form(";%s #psi_{%d};Entries",TPCNames[iDet].Data(),fHarmonic),200,0,TMath::Pi());
        fHistEventPlaneV0[iDet] = new TH1F(Form("fHistEventPlane%s",V0Names[iDet].Data()),Form(";%s #psi_{%d};Entries",V0Names[iDet].Data(),fHarmonic),200,0,TMath::Pi());
    	fHistqnVsCentrTPC[iDet] = new TH2F(Form("fHistqnVsCentr%s",TPCNames[iDet].Data()),Form(";#it{q}_{%d}^{%s};Entries",fHarmonic,TPCNames[iDet].Data()),100,0.,100.,15000,0.,15.);
		fHistqnVsCentrV0[iDet] = new TH2F(Form("fHistqnVsCentr%s",V0Names[iDet].Data()),Form(";#it{q}_{%d}^{%s};Entries",fHarmonic,V0Names[iDet].Data()),100,0.,100.,15000,0.,15.);
        fOutputList->Add(fHistEventPlaneTPC[iDet]);
        fOutputList->Add(fHistEventPlaneV0[iDet]);
        fOutputList->Add(fHistqnVsCentrTPC[iDet]);
        fOutputList->Add(fHistqnVsCentrV0[iDet]);
	}

    if(fCalibType==AliHFQnVectorHandler::kQnCalib) {
        if(fEnableTPCPhiVsCentrDistr) {
            OpenFile(2);
            fTPCPhiVsCentrDistr[0] = dynamic_cast<TH2F*>(fHFQnVecHandler->GetPhiDistrHistosTPCPosEta()->Clone());
            OpenFile(3);
            fTPCPhiVsCentrDistr[1] = dynamic_cast<TH2F*>(fHFQnVecHandler->GetPhiDistrHistosTPCNegEta()->Clone());
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
    if(fCalibType==AliHFQnVectorHandler::kQnCalib) {
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
void AliAnalysisTaskSEHFTenderQnVectors::UserExec(Option_t */*option*/)
{
    // main event loop
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAOD && AODEvent() && IsStandardAOD()) {
        // In case there is an AOD handler writing a standard AOD, use the AOD
        // event in memory rather than the input (ESD) event.
        fAOD = dynamic_cast<AliAODEvent*> (AODEvent());
    }
    if (!fAOD) {
        AliWarning("AliAnalysisTaskSEHFTenderQnVectors::Exec(): bad AOD");
        return;
    }
    fHistNEvents->Fill(1);

    if(TMath::Abs(fAOD->GetMagneticField())<0.001) return;

    AliAODHandler* aodHandler = static_cast<AliAODHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    if(!aodHandler) {
        AliWarning("AliAnalysisTaskSEHFTenderQnVectors::Exec(): No AliInputEventHandler!");
        return;
    }

    unsigned int maskPhysSel = aodHandler->IsEventSelected();
    TString firedTriggerClasses = fAOD->GetFiredTriggerClasses();
    if((fAOD->GetRunNumber()<136851 || fAOD->GetRunNumber()>139517)) {
        if(!(firedTriggerClasses.Contains(fTriggerClass.Data()))) return;
    }
    if((maskPhysSel & fTriggerMask)==0) {
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
    fHistCentrality->Fill(cent);

    const AliAODVertex *vertex = fAOD->GetPrimaryVertex();
    if(!vertex || TMath::Abs(vertex->GetZ())>10. || vertex->GetNContributors()<=0) {
        fHistNEvents->Fill(4);
        return;
    }

    fHFQnVecHandler->ResetAODEvent();
    fHFQnVecHandler->SetAODEvent(fAOD);
    fHFQnVecHandler->ComputeCalibratedQnVectorTPC();
    fHFQnVecHandler->ComputeCalibratedQnVectorV0();

	//fill histos with EP angle
    double PsinFullTPC = -1., PsinPosTPC = -1., PsinNegTPC = -1.;
    double PsinFullV0 = -1., PsinV0A = -1., PsinV0C = -1.;
    fHFQnVecHandler->GetEventPlaneAngleTPC(PsinFullTPC,PsinPosTPC,PsinNegTPC);
    fHFQnVecHandler->GetEventPlaneAngleV0(PsinFullV0,PsinV0A,PsinV0C);

    fHistEventPlaneTPC[0]->Fill(PsinFullTPC);
    fHistEventPlaneTPC[1]->Fill(PsinPosTPC);
    fHistEventPlaneTPC[2]->Fill(PsinNegTPC);
    fHistEventPlaneV0[0]->Fill(PsinFullV0);
    fHistEventPlaneV0[1]->Fill(PsinV0A);
    fHistEventPlaneV0[2]->Fill(PsinV0C);

	//fill histos for q2 spline calibration
    double qnFullTPC = -1., qnPosTPC = -1., qnNegTPC = -1.;
    double qnFullV0 = -1., qnV0A = -1., qnV0C = -1.;
    fHFQnVecHandler->GetqnTPC(qnFullTPC,qnPosTPC,qnNegTPC);
    fHFQnVecHandler->GetqnV0(qnFullV0,qnV0A,qnV0C);

	fHistqnVsCentrTPC[0]->Fill(cent,qnFullTPC);
	fHistqnVsCentrTPC[1]->Fill(cent,qnPosTPC);
	fHistqnVsCentrTPC[2]->Fill(cent,qnNegTPC);
	fHistqnVsCentrV0[0]->Fill(cent,qnFullV0);
	fHistqnVsCentrV0[1]->Fill(cent,qnV0A);
	fHistqnVsCentrV0[2]->Fill(cent,qnV0C);

    int runnumber = fAOD->GetRunNumber();

    if(fCalibType==AliHFQnVectorHandler::kQnCalib) {
        if(fEnableTPCPhiVsCentrDistr) {
            TH2F* hPhiPosEta = fHFQnVecHandler->GetPhiDistrHistosTPCPosEta();
            TH2F* hPhiNegEta = fHFQnVecHandler->GetPhiDistrHistosTPCNegEta();
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
            fHFQnVecHandler->GetUnNormQnVecTPC(QnVecFullTPC, QnVecPosTPC, QnVecNegTPC);
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

//________________________________________________________________________
TList* AliAnalysisTaskSEHFTenderQnVectors::GetSplineForqnPercentileList(int det) const
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
void AliAnalysisTaskSEHFTenderQnVectors::LoadSplinesForqnPercentile(TString splinesfilepath)
{
    // load splines for qn percentiles

    TString listnameTPC[3] = {"SplineListq2TPC", "SplineListq2TPCPosEta", "SplineListq2TPCNegEta"};
    TString listnameV0[3] = {"SplineListq2V0", "SplineListq2V0A", "SplineListq2V0C"};

    TFile* splinesfile = TFile::Open(splinesfilepath.Data());
    if(!splinesfile)
        AliFatal("File with splines for qn percentiles not found!");

    for(int iDet=0; iDet<3; iDet++) {
        fSplineListqnPercTPC[iDet] = (TList*)splinesfile->Get(listnameTPC[iDet].Data());
        if(!fSplineListqnPercTPC[iDet])
            AliFatal("TList with splines for qnTPC percentiles not found in the spline file!");
        fSplineListqnPercV0[iDet] = (TList*)splinesfile->Get(listnameV0[iDet].Data());
        if(!fSplineListqnPercV0[iDet])
            AliFatal("TList with splines for qnV0 percentiles not found in the spline file!");

        fSplineListqnPercTPC[iDet]->SetOwner();
        fSplineListqnPercV0[iDet]->SetOwner();
    }
    splinesfile->Close();
}
