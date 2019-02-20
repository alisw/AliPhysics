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

#include <TMath.h>
#include <TChain.h>

#include "AliAnalysisTaskSEHFTenderQnVectors.h"
#include "AliAODHandler.h"
#include "AliMultSelection.h"

ClassImp(AliAnalysisTaskSEHFTenderQnVectors)

//________________________________________________________________________
AliAnalysisTaskSEHFTenderQnVectors::AliAnalysisTaskSEHFTenderQnVectors() :
    AliAnalysisTaskSE("HFTenderQnVectors"),
    fOutputList(nullptr),
    fHistNEvents(nullptr),
    fHistCentrality(nullptr),
    fHFQnVecHandler(nullptr),
    fHarmonic(2),
    fCalibType(AliHFQnVectorHandler::kQnCalib),
    fNormMethod(AliHFQnVectorHandler::kQoverM),
    fOADBFileName(""),
    fAOD(nullptr),
    fTriggerClass(""),
    fTriggerMask(AliVEvent::kAny)
{
    //
    // default constructor
    //
    for(int iDet=0; iDet<3; iDet++) {
        fHistEventPlaneTPC[iDet] = nullptr;
        fHistEventPlaneV0[iDet] = nullptr;
    }
}

//________________________________________________________________________
AliAnalysisTaskSEHFTenderQnVectors::AliAnalysisTaskSEHFTenderQnVectors(const char *name, int harmonic, int calibType, TString oadbFileName) :
    AliAnalysisTaskSE(name),
    fOutputList(nullptr),
    fHistNEvents(nullptr),
    fHistCentrality(nullptr),
    fHFQnVecHandler(nullptr),
    fHarmonic(harmonic),
    fCalibType(calibType),
    fNormMethod(AliHFQnVectorHandler::kQoverSqrtM),
    fOADBFileName(oadbFileName),
    fAOD(nullptr),
    fTriggerClass(""),
    fTriggerMask(AliVEvent::kAny)
{
    //
    // standard constructor
    //
    for(int iDet=0; iDet<3; iDet++) {
        fHistEventPlaneTPC[iDet] = nullptr;
        fHistEventPlaneV0[iDet] = nullptr;
    }
    
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskSEHFTenderQnVectors::~AliAnalysisTaskSEHFTenderQnVectors()
{
    //
    // standard destructor
    //
    if(fOutputList) {
        delete fHistNEvents;
        delete fHistCentrality;
        for(int iDet=0; iDet<3; iDet++) {
            delete fHistEventPlaneTPC[iDet];
            delete fHistEventPlaneV0[iDet];
        }
        delete fOutputList;
    }
    if(fHFQnVecHandler) delete fHFQnVecHandler;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFTenderQnVectors::UserCreateOutputObjects()
{
    fHFQnVecHandler = new AliHFQnVectorHandler(fCalibType,fNormMethod,fHarmonic,fOADBFileName);

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
        fHistEventPlaneTPC[iDet] = new TH1F(Form("fHistEventPlane%s",TPCNames[iDet].Data()),Form(";#psi_{%d};Entries",fHarmonic),200,0,TMath::Pi());
        fHistEventPlaneV0[iDet] = new TH1F(Form("fHistEventPlane%s",V0Names[iDet].Data()),Form(";#psi_{%d};Entries",fHarmonic),200,0,TMath::Pi());
        fOutputList->Add(fHistEventPlaneTPC[iDet]);
        fOutputList->Add(fHistEventPlaneV0[iDet]);
    }

    // post data
    PostData(1, fOutputList);
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
    
    // Post output data
    PostData(1, fOutputList);
}
