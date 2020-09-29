/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************************
// \class AliHFQnVectorHandler
// \brief helper class to handle the Qn-vectors computed with different calibrations
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
// F. Catalano, fabio.catalano@cern.ch
// A. Dobrin, alexandru.dobrin@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// G. Luparello, grazia.luparello@cern.ch
// F. Prino, prino@to.infn.it
// A. Rossi, andrea.rossi@cern.ch
// S. Trogolo, stefano.trogolo@cern.ch
///////////////////////////////////////////////////////////////////////////////////////

#include <TGrid.h>

#include "AliHFQnVectorHandler.h"
#include "AliAODVertex.h"
#include "AliOADBContainer.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliMultSelection.h"

/// \cond CLASSIMP
ClassImp(AliHFQnVectorHandler);
/// \endcond

//________________________________________________________________
AliHFQnVectorHandler::AliHFQnVectorHandler():
    fHarmonic(2),
    fNormMethod(kQoverM),
    fCalibType(kQnCalib),
    fEtaMinTPC(0.),
    fEtaMaxTPC(0.8),
    fPtMinTPC(0.2),
    fPtMaxTPC(5.),
    fMultFullTPC(0.),
    fMultPosTPC(0.),
    fMultNegTPC(0.),
    fMultFullV0(0.),
    fMultV0A(0.),
    fMultV0C(0.),
    fQnVecNormFullTPC(1.),
    fQnVecNormPosTPC(1.),
    fQnVecNormNegTPC(1.),
    fQnVecNormFullV0(1.),
    fQnVecNormV0A(1.),
    fQnVecNormV0C(1.),
    fUsedTrackPosIDs(),
    fUsedTrackNegIDs(),
    fFractionOfTracksForQnTPC(1.1),
    fAODEvent(nullptr),
    fV0(nullptr),
    fRun(0),
    fZvtx(-9999.),
    fCentrality(-9999.),
    fOADBFile(nullptr),
    fOADBFileName(""),
    fIsOADBFileOpen(false),
    fCalibObjRun(-9999),
    fHistMultV0(nullptr),
    fV0CalibZvtxDiff(true),
    fEnablePhiDistrHistos(false),
    fQnVectorTask(nullptr),
    fQnVectorMgr(nullptr)
{
    //
    // Default constructor
    //
    fUsedTrackPosIDs = TBits(1000);
    fUsedTrackNegIDs = TBits(1000);

    for(int iComp=0; iComp<2; iComp++) {
        fQnVecFullTPC[iComp]     = -9999.;
        fQnVecPosTPC[iComp]      = -9999.;
        fQnVecNegTPC[iComp]      = -9999.;
        fQnVecFullV0[iComp]      = -9999.;
        fQnVecV0A[iComp]         = -9999.;
        fQnVecV0C[iComp]         = -9999.;
    }

    for(int iZvtx = 0; iZvtx < 14; iZvtx++) {
        fQx2mV0A[iZvtx] = nullptr;
        fQy2mV0A[iZvtx] = nullptr;
        fQx2sV0A[iZvtx] = nullptr;
        fQy2sV0A[iZvtx] = nullptr;
        fQx2mV0C[iZvtx] = nullptr;
        fQy2mV0C[iZvtx] = nullptr;
        fQx2sV0C[iZvtx] = nullptr;
        fQy2sV0C[iZvtx] = nullptr;  
    }

    for(int iCent = 0; iCent < 9; iCent++) {
        fWeightsTPCPosEta[iCent] = nullptr;
        fWeightsTPCNegEta[iCent] = nullptr;
    }

    fPhiVsCentrTPC[0]=nullptr;
    fPhiVsCentrTPC[1]=nullptr;
}

//________________________________________________________________
AliHFQnVectorHandler::AliHFQnVectorHandler(int calibType, int normMeth, int harmonic, TString OADBfileName):
    fHarmonic(harmonic),
    fNormMethod(normMeth),
    fCalibType(calibType),
    fEtaMinTPC(0.),
    fEtaMaxTPC(0.8),
    fPtMinTPC(0.2),
    fPtMaxTPC(5.),
    fMultFullTPC(0.),
    fMultPosTPC(0.),
    fMultNegTPC(0.),
    fMultFullV0(0.),
    fMultV0A(0.),
    fMultV0C(0.),
    fQnVecNormFullTPC(1.),
    fQnVecNormPosTPC(1.),
    fQnVecNormNegTPC(1.),
    fQnVecNormFullV0(1.),
    fQnVecNormV0A(1.),
    fQnVecNormV0C(1.),
    fUsedTrackPosIDs(),
    fUsedTrackNegIDs(),
    fFractionOfTracksForQnTPC(1.1),
    fAODEvent(nullptr),
    fV0(nullptr),
    fRun(0),
    fZvtx(-9999.),
    fCentrality(-9999.),
    fOADBFile(nullptr),
    fOADBFileName(OADBfileName),
    fIsOADBFileOpen(false),
    fCalibObjRun(-9999),
    fHistMultV0(nullptr),
    fV0CalibZvtxDiff(true),
    fEnablePhiDistrHistos(false),
    fQnVectorTask(nullptr),
    fQnVectorMgr(nullptr)
{
    //
    // Standard constructor
    //
    fUsedTrackPosIDs = TBits(1000);
    fUsedTrackNegIDs = TBits(1000);

    for(int iComp=0; iComp<2; iComp++) {
        fQnVecFullTPC[iComp]     = -9999.;
        fQnVecPosTPC[iComp]      = -9999.;
        fQnVecNegTPC[iComp]      = -9999.;
        fQnVecFullV0[iComp]      = -9999.;
        fQnVecV0A[iComp]         = -9999.;
        fQnVecV0C[iComp]         = -9999.;
    }

    for(int iZvtx = 0; iZvtx < 14; iZvtx++) {
        fQx2mV0A[iZvtx] = nullptr;
        fQy2mV0A[iZvtx] = nullptr;
        fQx2sV0A[iZvtx] = nullptr;
        fQy2sV0A[iZvtx] = nullptr;
        fQx2mV0C[iZvtx] = nullptr;
        fQy2mV0C[iZvtx] = nullptr;
        fQx2sV0C[iZvtx] = nullptr;
        fQy2sV0C[iZvtx] = nullptr;  
    }

    for(int iCent = 0; iCent < 9; iCent++) {
        fWeightsTPCPosEta[iCent] = nullptr;
        fWeightsTPCNegEta[iCent] = nullptr;
    }

    fPhiVsCentrTPC[0]=nullptr;
    fPhiVsCentrTPC[1]=nullptr;
}

//________________________________________________________________
AliHFQnVectorHandler::~AliHFQnVectorHandler()
{
    //
    // Standard destructor
    //
    if(fPhiVsCentrTPC[0]) delete fPhiVsCentrTPC[0];
    if(fPhiVsCentrTPC[1]) delete fPhiVsCentrTPC[1];
}

//________________________________________________________________
bool AliHFQnVectorHandler::SetAODEvent(AliAODEvent* event)
{
    if(!event) {
        AliWarning("Event not found!");
        return false;
    }
    fRun = event->GetRunNumber();
    const AliAODVertex* trkVtx = dynamic_cast<AliAODVertex*>(event->GetPrimaryVertex());
    if (!trkVtx || trkVtx->GetNContributors()<=0) {
        AliWarning("Primary vertex not found!");
        return false;
    }
    else
        fZvtx = trkVtx->GetZ();

    AliMultSelection* MultSelection = dynamic_cast<AliMultSelection*>(event->FindListObject("MultSelection"));
    if(!MultSelection) {
        AliWarning("AliMultSelection object not found!");
        return false;
    }    
    else 
        fCentrality = MultSelection->GetMultiplicityPercentile("V0M");

    fV0 = dynamic_cast<AliAODVZERO*>(event->GetVZEROData());
    if(!fV0) {
        AliWarning("V0 info not found!");
        return false;
    }

    //only if everything ok get event
    fAODEvent = event;

    if(fCalibType!=kQnFrameworkCalib)
        return true;

    //get Qn-framework objects if needed
    fQnVectorTask = static_cast<AliAnalysisTaskFlowVectorCorrections*>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
    if (!fQnVectorTask) {
        AliWarning("Flow Qn vector corrections framework not found!\n");
        return false;
    }
    fQnVectorMgr = fQnVectorTask->GetAliQnCorrectionsManager();
    return true;
}

//________________________________________________________________
void AliHFQnVectorHandler::ResetAODEvent() 
{
    fAODEvent       = nullptr;
    fRun            = -9999;
    fCentrality     = -9999.;
    fZvtx           = -9999.;
    fV0             = nullptr;
    fQnVectorTask   = nullptr;
    fQnVectorMgr    = nullptr;

    //reset all Q-vectors
    for(int iComp=0; iComp<2; iComp++) {
        fQnVecFullTPC[iComp]     = -9999.;
        fQnVecPosTPC[iComp]      = -9999.;
        fQnVecNegTPC[iComp]      = -9999.;
        fQnVecFullV0[iComp]      = -9999.;
        fQnVecV0A[iComp]     = -9999.;
        fQnVecV0C[iComp]     = -9999.;
    }

    fMultFullTPC      = -9999.;
    fMultPosTPC       = -9999.;
    fMultNegTPC       = -9999.;
    fMultFullV0       = -9999.;
    fMultV0A          = -9999.;
    fMultV0C          = -9999.;
    fQnVecNormFullTPC = -9999.;
    fQnVecNormPosTPC  = -9999.;
    fQnVecNormNegTPC  = -9999.;
    fQnVecNormFullV0  = -9999.;
    fQnVecNormV0A     = -9999.;
    fQnVecNormV0C     = -9999.;

    fUsedTrackPosIDs.ResetAllBits();
    fUsedTrackNegIDs.ResetAllBits();
}

//________________________________________________________________
bool AliHFQnVectorHandler::ComputeCalibratedQnVectorTPC(bool forceCalc) 
{
    if(fQnVecFullTPC[0] == -9999. || forceCalc) {
        if(!fAODEvent) {
            AliError("AOD Event not set!");
            return false;
        }
        if(fCalibType==kQnCalib) {
            if(!fIsOADBFileOpen || fCalibObjRun!=fRun) {
                fIsOADBFileOpen = OpenInfoCalbration();
                if(!fIsOADBFileOpen)
                    return false;
                fCalibObjRun = fRun;
            }
            ComputeQvecTPC();
        }
        else if(fCalibType==kQnFrameworkCalib && fQnVectorTask && fQnVectorMgr) {
            ComputeQvecQnFrameworkTPC();
        }
    }

    return true;
}

//________________________________________________________________
bool AliHFQnVectorHandler::ComputeCalibratedQnVectorV0(bool forceCalc) 
{
    if(fQnVecFullV0[0] == -9999. || forceCalc) {
        if(!fAODEvent) {
            AliError("AOD Event not set!");
            return false;
        }
        if(fCalibType==kQnCalib) {
            if(!fIsOADBFileOpen || fCalibObjRun!=fRun) {
                fIsOADBFileOpen = OpenInfoCalbration();
                if(!fIsOADBFileOpen)
                    return false;
                fCalibObjRun = fRun;
            }
            ComputeQvecV0();
        }
        else if(fCalibType==kQnFrameworkCalib && fQnVectorTask && fQnVectorMgr) {
            ComputeQvecQnFrameworkV0();
        }
    }

    return true;
}

//________________________________________________________________
void AliHFQnVectorHandler::GetQnVecTPC(double QnVecFullTPC[2], double QnVecPosTPC[2], double QnVecNegTPC[2])
{
    switch(fNormMethod) {
        case 0: //kQoverQlength
            for(int iComp=0; iComp<2; iComp++) {
                fQnVecNormFullTPC>0 ? QnVecFullTPC[iComp] = fQnVecFullTPC[iComp] / fQnVecNormFullTPC : QnVecFullTPC[iComp] = fQnVecFullTPC[iComp];
                fQnVecNormPosTPC>0 ? QnVecPosTPC[iComp]  = fQnVecPosTPC[iComp] / fQnVecNormPosTPC : QnVecPosTPC[iComp]  = fQnVecPosTPC[iComp];
                fQnVecNormNegTPC>0 ? QnVecNegTPC[iComp]  = fQnVecNegTPC[iComp] / fQnVecNormNegTPC : QnVecNegTPC[iComp]  = fQnVecNegTPC[iComp];
            }
        break;
        case 1: //kQoverM
            for(int iComp=0; iComp<2; iComp++) {
                fMultFullTPC>0 ? QnVecFullTPC[iComp] = fQnVecFullTPC[iComp] / fMultFullTPC : QnVecFullTPC[iComp] = fQnVecFullTPC[iComp];
                fMultPosTPC>0 ? QnVecPosTPC[iComp]  = fQnVecPosTPC[iComp] / fMultPosTPC : QnVecPosTPC[iComp]  = fQnVecPosTPC[iComp];
                fMultNegTPC>0 ? QnVecNegTPC[iComp]  = fQnVecNegTPC[iComp] / fMultNegTPC : QnVecNegTPC[iComp]  = fQnVecNegTPC[iComp];
            }
        break;
        case 2: //kQoverSqrtM
            for(int iComp=0; iComp<2; iComp++) {
                fMultFullTPC>0 ? QnVecFullTPC[iComp] = fQnVecFullTPC[iComp] / TMath::Sqrt(fMultFullTPC) : QnVecFullTPC[iComp] = fQnVecFullTPC[iComp];
                fMultPosTPC>0 ? QnVecPosTPC[iComp]  = fQnVecPosTPC[iComp] / TMath::Sqrt(fMultPosTPC) : QnVecPosTPC[iComp]  = fQnVecPosTPC[iComp];
                fMultNegTPC>0 ? QnVecNegTPC[iComp]  = fQnVecNegTPC[iComp] / TMath::Sqrt(fMultNegTPC) : QnVecNegTPC[iComp]  = fQnVecNegTPC[iComp];
            }
        break;
        case 3: //kNone
            for(int iComp=0; iComp<2; iComp++) {
                QnVecFullTPC[iComp] = fQnVecFullTPC[iComp];
                QnVecPosTPC[iComp]  = fQnVecPosTPC[iComp];
                QnVecNegTPC[iComp]  = fQnVecNegTPC[iComp];
            }
        break;
    }
}

//________________________________________________________________
void AliHFQnVectorHandler::GetQnVecV0(double QnVecFullV0[2], double QnVecV0A[2], double QnVecV0C[2])
{
    switch(fNormMethod) {
        case 0: //kQoverQlength
            for(int iComp=0; iComp<2; iComp++) {
                fQnVecNormFullV0 > 0 ? QnVecFullV0[iComp] = fQnVecFullV0[iComp] / fQnVecNormFullV0 : QnVecFullV0[iComp] = fQnVecFullV0[iComp];
                fQnVecNormV0A > 0 ? QnVecV0A[iComp]  = fQnVecV0A[iComp] / fQnVecNormV0A : QnVecV0A[iComp] = fQnVecV0A[iComp];
                fQnVecNormV0C > 0 ? QnVecV0C[iComp]  = fQnVecV0C[iComp] / fQnVecNormV0C : QnVecV0C[iComp] = fQnVecV0C[iComp];
            }
        break;
        case 1: //kQoverM
            for(int iComp=0; iComp<2; iComp++) {
                fMultFullV0 > 0 ? QnVecFullV0[iComp] = fQnVecFullV0[iComp] / fMultFullV0 : QnVecFullV0[iComp] = fQnVecFullV0[iComp];
                fMultV0A > 0 ? QnVecV0A[iComp]  = fQnVecV0A[iComp] / fMultV0A : QnVecV0A[iComp] = fQnVecV0A[iComp];
                fMultV0C > 0 ? QnVecV0C[iComp]  = fQnVecV0C[iComp] / fMultV0C : QnVecV0C[iComp] = fQnVecV0C[iComp];
            }
        break;
        case 2: //kQoverSqrtM
            for(int iComp=0; iComp<2; iComp++) {
                fMultFullV0 > 0 ? QnVecFullV0[iComp] = fQnVecFullV0[iComp] / TMath::Sqrt(fMultFullV0) : QnVecFullV0[iComp] = fQnVecFullV0[iComp];
                fMultV0A > 0 ? QnVecV0A[iComp]  = fQnVecV0A[iComp] / TMath::Sqrt(fMultV0A) : QnVecV0A[iComp] = fQnVecV0A[iComp];
                fMultV0C > 0 ? QnVecV0C[iComp]  = fQnVecV0C[iComp] / TMath::Sqrt(fMultV0C) : QnVecV0C[iComp] = fQnVecV0C[iComp];
            }
        break;
        case 3: //kNone
            for(int iComp=0; iComp<2; iComp++) {
                QnVecFullV0[iComp] = fQnVecFullV0[iComp];
                QnVecV0A[iComp]  = fQnVecV0A[iComp];
                QnVecV0C[iComp]  = fQnVecV0C[iComp];
            }
        break;
    }
}

//________________________________________________________________
void AliHFQnVectorHandler::GetUnNormQnVecTPC(double QnVecFullTPC[2], double QnVecPosTPC[2], double QnVecNegTPC[2])
{
    for(int iComp=0; iComp<2; iComp++) {
        QnVecFullTPC[iComp] = fQnVecFullTPC[iComp];
        QnVecPosTPC[iComp]  = fQnVecPosTPC[iComp];
        QnVecNegTPC[iComp]  = fQnVecNegTPC[iComp];
    }
}

//________________________________________________________________
void AliHFQnVectorHandler::GetUnNormQnVecV0(double QnVecFullV0[2], double QnVecV0A[2], double QnVecV0C[2])
{
    for(int iComp=0; iComp<2; iComp++) {
        QnVecFullV0[iComp] = fQnVecFullV0[iComp];
        QnVecV0A[iComp]  = fQnVecV0A[iComp];
        QnVecV0C[iComp]  = fQnVecV0C[iComp];
    }
}

//________________________________________________________________
void AliHFQnVectorHandler::GetqnTPC(double &qnFullTPC, double &qnPosTPC, double &qnNegTPC) 
{
    qnFullTPC = fQnVecNormFullTPC / TMath::Sqrt(fMultFullTPC);
    qnPosTPC  = fQnVecNormPosTPC / TMath::Sqrt(fMultPosTPC);
    qnNegTPC  = fQnVecNormNegTPC / TMath::Sqrt(fMultNegTPC);
}
    
//________________________________________________________________
void AliHFQnVectorHandler::GetqnV0(double &qnFullV0, double &qnV0A, double &qnV0C)
{
    qnFullV0 = fQnVecNormFullV0 / TMath::Sqrt(fMultFullV0);
    qnV0A    = fQnVecNormV0A / TMath::Sqrt(fMultV0A);
    qnV0C    = fQnVecNormV0C / TMath::Sqrt(fMultV0C);
}

//________________________________________________________________
void AliHFQnVectorHandler::GetMultQnVecTPC(double &MultFullTPC, double &MultPosTPC, double &MultNegTPC)
{
    MultFullTPC = fMultFullTPC;
    MultPosTPC  = fMultPosTPC;
    MultNegTPC  = fMultNegTPC;
}

//________________________________________________________________
void AliHFQnVectorHandler::GetMultQnVecV0(double &MultFullV0, double &MultV0A, double &MultV0C)
{
    MultFullV0 = fMultFullV0;
    MultV0A  = fMultV0A;
    MultV0C  = fMultV0C;
}

//________________________________________________________________
void AliHFQnVectorHandler::GetEventPlaneAngleTPC(double &EvPlaneFullTPC, double &EvPlanePosTPC, double &EvPlaneNegTPC)
{
    EvPlaneFullTPC = ComputeEventPlaneAngle(fQnVecFullTPC[0],fQnVecFullTPC[1]);
    EvPlanePosTPC  = ComputeEventPlaneAngle(fQnVecPosTPC[0],fQnVecPosTPC[1]);
    EvPlaneNegTPC  = ComputeEventPlaneAngle(fQnVecNegTPC[0],fQnVecNegTPC[1]);
}

//________________________________________________________________
void AliHFQnVectorHandler::GetEventPlaneAngleV0(double &EvPlaneFullV0, double &EvPlaneV0A, double &EvPlaneV0C)
{
    EvPlaneFullV0 = ComputeEventPlaneAngle(fQnVecFullV0[0],fQnVecFullV0[1]);
    EvPlaneV0A    = ComputeEventPlaneAngle(fQnVecV0A[0],fQnVecV0A[1]);
    EvPlaneV0C    = ComputeEventPlaneAngle(fQnVecV0C[0],fQnVecV0C[1]);
}

//________________________________________________________________
void AliHFQnVectorHandler::RemoveTracksFromQnTPC(vector<AliAODTrack*> trToRem, double QnVecFullTPC[2], double QnVecPosTPC[2], double QnVecNegTPC[2], double &multFullTPC, double &multPosTPC, double &multNegTPC, bool getUnNormalised) 
{
    for(int iComp=0; iComp<2; iComp++) {
        QnVecFullTPC[iComp] = fQnVecFullTPC[iComp];
        QnVecPosTPC[iComp]  = fQnVecPosTPC[iComp];
        QnVecNegTPC[iComp]  = fQnVecNegTPC[iComp];
    }
    multFullTPC = fMultFullTPC;
    multPosTPC = fMultPosTPC;
    multNegTPC = fMultNegTPC;

    TBits posIDsAfterRem(fUsedTrackPosIDs);
    TBits negIDsAfterRem(fUsedTrackNegIDs);
    
    if(fCalibType==kQnCalib) {
        
        short centbin = GetCentBin();

        for(unsigned int iTr=0; iTr<trToRem.size(); iTr++) {
            bool isTrackUsed = false;
            short trID = trToRem[iTr]->GetID();
            if(trID>=0) 
                isTrackUsed = posIDsAfterRem.TestBitNumber(trID);
            else 
                isTrackUsed = negIDsAfterRem.TestBitNumber(TMath::Abs(trID));
            if(!isTrackUsed) continue; // --> track not used for Qn
            double eta=trToRem[iTr]->Eta();
            double phi=trToRem[iTr]->Phi();
            double weight = 1.;
            if(eta>0) {
                if(fWeightsTPCPosEta[centbin]) {
                    int phibin = fWeightsTPCPosEta[centbin]->GetXaxis()->FindBin(phi);
                    weight = 1./fWeightsTPCPosEta[centbin]->GetBinContent(phibin);
                }
                QnVecFullTPC[0] -= weight*TMath::Cos(fHarmonic*phi);
                QnVecFullTPC[1] -= weight*TMath::Sin(fHarmonic*phi);
                multFullTPC     -= weight;
                QnVecPosTPC[0]  -= weight*TMath::Cos(fHarmonic*phi);
                QnVecPosTPC[1]  -= weight*TMath::Sin(fHarmonic*phi);
                multPosTPC      -= weight;
            }
            else {
                if(fWeightsTPCNegEta[centbin]) {
                    int phibin = fWeightsTPCNegEta[centbin]->GetXaxis()->FindBin(phi);
                    weight = 1./fWeightsTPCNegEta[centbin]->GetBinContent(phibin);
                }
                QnVecFullTPC[0] -= weight*TMath::Cos(fHarmonic*phi);
                QnVecFullTPC[1] -= weight*TMath::Sin(fHarmonic*phi);
                multFullTPC     -= weight;
                QnVecNegTPC[0]  -= weight*TMath::Cos(fHarmonic*phi);
                QnVecNegTPC[1]  -= weight*TMath::Sin(fHarmonic*phi);
                multNegTPC      -= weight;
            }
            if(trID>=0)
                posIDsAfterRem.ResetBitNumber(trID);
            else 
                negIDsAfterRem.ResetBitNumber(TMath::Abs(trID));
        }
    }
    else if(fCalibType==kQnFrameworkCalib) {

        TString detTPCConfName[3] = {"TPC","TPCPosEta","TPCNegEta"};

        double qxToRem[3]    = {0.,0.,0.};
        double qyToRem[3]    = {0.,0.,0.};
        double multToRem[3]  = {0.,0.,0.};
        double M[3]          = {multFullTPC,multPosTPC,multNegTPC};
        double qxAfterRem[3] = {QnVecFullTPC[0],QnVecPosTPC[0],QnVecNegTPC[0]};
        double qyAfterRem[3] = {QnVecFullTPC[1],QnVecPosTPC[1],QnVecNegTPC[1]};
        
        for(unsigned int iTr=0; iTr<trToRem.size(); iTr++) {
            bool isTrackUsed = false;
            short trID = trToRem[iTr]->GetID();
            if(trID>=0) 
                isTrackUsed = posIDsAfterRem.TestBitNumber(trID);
            else 
                isTrackUsed = negIDsAfterRem.TestBitNumber(TMath::Abs(trID));
            if(!isTrackUsed) continue; // --> track not used for Qn   
            double eta=trToRem[iTr]->Eta();
            double phi=trToRem[iTr]->Phi();
            qxToRem[0]   += TMath::Cos(fHarmonic*phi);
            qyToRem[0]   += TMath::Sin(fHarmonic*phi);
            multToRem[0] += 1;
            if(eta>0) {
                qxToRem[1]   += TMath::Cos(fHarmonic*phi);
                qyToRem[1]   += TMath::Sin(fHarmonic*phi);
                multToRem[1] += 1;
            }
            else {
                qxToRem[2]   += TMath::Cos(fHarmonic*phi);
                qyToRem[2]   += TMath::Sin(fHarmonic*phi);
                multToRem[2] += 1;
            }
            if(trID>=0)
                posIDsAfterRem.ResetBitNumber(trID);
            else 
                negIDsAfterRem.ResetBitNumber(TMath::Abs(trID));
        }

        TList* qnlist = dynamic_cast<TList*>(fQnVectorMgr->GetQnVectorList());
        for(int iDet=0; iDet<3; iDet++) {
            TList *pQvecList = dynamic_cast<TList*>(qnlist->FindObject(detTPCConfName[iDet].Data()));
            if (pQvecList) {
                /* the detector is present */
                AliQnCorrectionsQnVector* theQnVectorUncorr = dynamic_cast<AliQnCorrectionsQnVector*>(pQvecList->FindObject("plain")); //raw step for TPC
                if (theQnVectorUncorr && theQnVectorUncorr->IsGoodQuality() && theQnVectorUncorr->GetN() > 0) { 

                    double qxPlain = theQnVectorUncorr->Qx(fHarmonic);
                    double qyPlain = theQnVectorUncorr->Qy(fHarmonic);

                    //track removal from rec step (if present)
                    AliQnCorrectionsQnVector* theQnVectorCorr = dynamic_cast<AliQnCorrectionsQnVector*>(pQvecList->FindObject("rec")); //rec step for TPC
                    double corrX = 0., corrY = 0., qxRec = -999., qyRec = -999.;
                    if (theQnVectorCorr && theQnVectorCorr->IsGoodQuality() && theQnVectorCorr->GetN() > 0) {
                        double qxRec = theQnVectorCorr->Qx(fHarmonic);
                        double qyRec = theQnVectorCorr->Qy(fHarmonic);
                        corrX = qxPlain - qxRec;
                        corrY = qyPlain - qyRec;
                    }
                    double qxRecSub = (qxPlain*M[iDet] - qxToRem[iDet]) / (M[iDet]-multToRem[iDet]) - corrX;
                    double qyRecSub = (qyPlain*M[iDet] - qyToRem[iDet]) / (M[iDet]-multToRem[iDet]) - corrY;

                    double qxTwist = -999., qyTwist = -999., Lbplus = 0., Lbminus = 0.;
                    //auto-correlation removal from twist step (if present)
                    theQnVectorCorr = dynamic_cast<AliQnCorrectionsQnVector*>(pQvecList->FindObject("twist")); //twist step for TPC
                    if (theQnVectorCorr && theQnVectorCorr->IsGoodQuality() && theQnVectorUncorr->GetN() > 0) {
                        qxTwist = theQnVectorCorr->Qx(fHarmonic);
                        qyTwist = theQnVectorCorr->Qy(fHarmonic);
                        Lbplus  = (qyRec - qyTwist)/qxTwist;
                        Lbminus = (qxRec - qxTwist)/qyTwist;
                        qxAfterRem[iDet] = (qxRecSub - Lbminus*qyRecSub)/(1-Lbminus*Lbplus);
                        qyAfterRem[iDet] = (qyRecSub - Lbplus*qxRecSub)/(1-Lbminus*Lbplus);
                    }
                    else {
                        qxAfterRem[iDet] = qxRecSub;
                        qyAfterRem[iDet] = qyRecSub;                        
                    }
                }
            }
            else {
                multToRem[iDet] = 0;
            }
        }
        QnVecFullTPC[0] = qxAfterRem[0];
        QnVecFullTPC[1] = qyAfterRem[0];
        QnVecPosTPC[0]  = qxAfterRem[1];
        QnVecPosTPC[1]  = qyAfterRem[1];
        QnVecNegTPC[0]  = qxAfterRem[2];
        QnVecNegTPC[1]  = qyAfterRem[2];
        multFullTPC     = M[0]-multToRem[0];
        multPosTPC      = M[1]-multToRem[1];
        multNegTPC      = M[2]-multToRem[2]; 
    }

    if(!getUnNormalised) {
        double QnVecNormFullTPC = 1., QnVecNormPosTPC = 1., QnVecNormNegTPC = 1.;
        switch(fNormMethod) {
            case 0: //kQoverQlength
                QnVecNormFullTPC = TMath::Sqrt(QnVecFullTPC[0]*QnVecFullTPC[0]+QnVecFullTPC[1]*QnVecFullTPC[1]);
                QnVecNormPosTPC  = TMath::Sqrt(QnVecPosTPC[0]*QnVecPosTPC[0]+QnVecPosTPC[1]*QnVecPosTPC[1]);
                QnVecNormNegTPC  = TMath::Sqrt(QnVecNegTPC[0]*QnVecNegTPC[0]+QnVecNegTPC[1]*QnVecNegTPC[1]);
                for(int iComp=0; iComp<2; iComp++) {
                    if(QnVecNormFullTPC > 0.) QnVecFullTPC[iComp] = QnVecFullTPC[iComp] / QnVecNormFullTPC;
                    if(QnVecNormPosTPC > 0.) QnVecPosTPC[iComp]   = QnVecPosTPC[iComp] / QnVecNormPosTPC;
                    if(QnVecNormNegTPC > 0.) QnVecNegTPC[iComp]   = QnVecNegTPC[iComp] / QnVecNormNegTPC;
                }
            break;
            case 1: //kQoverM
                for(int iComp=0; iComp<2; iComp++) {
                    if(multFullTPC > 0.) QnVecFullTPC[iComp] = QnVecFullTPC[iComp] / multFullTPC;
                    if(multPosTPC > 0.) QnVecPosTPC[iComp]   = QnVecPosTPC[iComp] / multPosTPC;
                    if(multNegTPC > 0.) QnVecNegTPC[iComp]   = QnVecNegTPC[iComp] / multNegTPC;
                }
            break;
            case 2: //kQoverSqrtM
                for(int iComp=0; iComp<2; iComp++) {
                    if(multFullTPC > 0.) QnVecFullTPC[iComp] = QnVecFullTPC[iComp] / TMath::Sqrt(multFullTPC);
                    if(multPosTPC > 0.)  QnVecPosTPC[iComp]  = QnVecPosTPC[iComp] / TMath::Sqrt(multPosTPC);
                    if(multNegTPC > 0.)  QnVecNegTPC[iComp]  = QnVecNegTPC[iComp] / TMath::Sqrt(multNegTPC);
                }
            break;
        }
    }
}

//________________________________________________________________
void AliHFQnVectorHandler::EnablePhiDistrHistos() 
{
    fEnablePhiDistrHistos=true;
    TString histonames[2] = {"TPCPosEta","TPCNegEta"};
    for(int iHisto=0; iHisto<2; iHisto++) {
        if(fPhiVsCentrTPC[iHisto]) {
            delete fPhiVsCentrTPC[iHisto];
            fPhiVsCentrTPC[iHisto]=nullptr;
        }
        fPhiVsCentrTPC[iHisto] = new TH2F(Form("fPhiVsCentrTPC%s",histonames[iHisto].Data()),";centrality (%);#varphi;Entries",10,0.,100.,180,0.,2*TMath::Pi());
    }
}

//________________________________________________________________
void AliHFQnVectorHandler::ComputeQvecQnFrameworkTPC() 
{
    if(!fQnVectorTask || !fQnVectorMgr) {
        AliWarning("Qn-framework task not found!");
        return;
    }
    
    TList* qnlist = dynamic_cast<TList*>(fQnVectorMgr->GetQnVectorList());
    if (!qnlist) return;  

    TString detTPCConfName[3] = {"TPC","TPCPosEta","TPCNegEta"};
    TString normMethodName = "";
    switch(fNormMethod) {
        case 0: //kQoverQlength
            normMethodName = "QoverQlength";
        break;
        case 1: //kQoverM
            normMethodName = "QoverM";
        break;
        case 2: //kQoverSqrtM
            normMethodName = "QoverSqrtM";
        break;
    }
    AliQnCorrectionsQnVector *theQnVector[3] = {nullptr,nullptr,nullptr};

    for(int iDet=0; iDet<3; iDet++) {
        TList *pQvecList = dynamic_cast<TList*>(qnlist->FindObject(Form("%s%s",detTPCConfName[iDet].Data(),normMethodName.Data())));
        if (!pQvecList) return;  

        /* the detector is present */
        theQnVector[iDet] = dynamic_cast<AliQnCorrectionsQnVector*>(pQvecList->First());
        if (!theQnVector[iDet] || !(theQnVector[iDet]->IsGoodQuality()) || !(theQnVector[iDet]->GetN() != 0)) {
            /* the Qn vector for the expected step was not there */
	        theQnVector[iDet] = dynamic_cast<AliQnCorrectionsQnVector*>(pQvecList->FindObject("plain"));
        }
    }

    //fill Q vectors
    if(theQnVector[0]) {
        fMultFullTPC = theQnVector[0]->GetSumOfWeights();
        switch(fNormMethod) {
            case 0: //kQoverQlength
                fQnVecFullTPC[0] = theQnVector[0]->Qx(fHarmonic)*theQnVector[0]->Length(fHarmonic);
                fQnVecFullTPC[1] = theQnVector[0]->Qy(fHarmonic)*theQnVector[0]->Length(fHarmonic);
            break;
            case 1: //kQoverM
                fQnVecFullTPC[0] = theQnVector[0]->Qx(fHarmonic)*fMultFullTPC;
                fQnVecFullTPC[1] = theQnVector[0]->Qy(fHarmonic)*fMultFullTPC;
            break;
            case 2: //kQoverSqrtM
                fQnVecFullTPC[0] = theQnVector[0]->Qx(fHarmonic)*TMath::Sqrt(fMultFullTPC);
                fQnVecFullTPC[1] = theQnVector[0]->Qy(fHarmonic)*TMath::Sqrt(fMultFullTPC);
            break;
        }
        fQnVecNormFullTPC = TMath::Sqrt(fQnVecFullTPC[0]*fQnVecFullTPC[0]+fQnVecFullTPC[1]*fQnVecFullTPC[1]);
    }
    if(theQnVector[1]) {
        fMultPosTPC = theQnVector[1]->GetSumOfWeights();
        switch(fNormMethod) {
            case 0: //kQoverQlength
                fQnVecPosTPC[0] = theQnVector[1]->Qx(fHarmonic)*theQnVector[1]->Length(fHarmonic);
                fQnVecPosTPC[1] = theQnVector[1]->Qy(fHarmonic)*theQnVector[1]->Length(fHarmonic);
            break;
            case 1: //kQoverM
                fQnVecPosTPC[0] = theQnVector[1]->Qx(fHarmonic)*fMultPosTPC;
                fQnVecPosTPC[1] = theQnVector[1]->Qy(fHarmonic)*fMultPosTPC;
            break;
            case 2: //kQoverSqrtM
                fQnVecPosTPC[0] = theQnVector[1]->Qx(fHarmonic)*TMath::Sqrt(fMultPosTPC);
                fQnVecPosTPC[1] = theQnVector[1]->Qy(fHarmonic)*TMath::Sqrt(fMultPosTPC);
            break;
        }
        fQnVecNormPosTPC  = TMath::Sqrt(fQnVecPosTPC[0]*fQnVecPosTPC[0]+fQnVecPosTPC[1]*fQnVecPosTPC[1]);
    }
    if(theQnVector[2]) {
        fMultNegTPC = theQnVector[2]->GetSumOfWeights();
        switch(fNormMethod) {
            case 0: //kQoverQlength
                fQnVecNegTPC[0] = theQnVector[2]->Qx(fHarmonic)*theQnVector[2]->Length(fHarmonic);
                fQnVecNegTPC[1] = theQnVector[2]->Qy(fHarmonic)*theQnVector[2]->Length(fHarmonic);
            break;
            case 1: //kQoverM
                fQnVecNegTPC[0] = theQnVector[2]->Qx(fHarmonic)*fMultNegTPC;
                fQnVecNegTPC[1] = theQnVector[2]->Qy(fHarmonic)*fMultNegTPC;
            break;
            case 2: //kQoverSqrtM
                fQnVecNegTPC[0] = theQnVector[2]->Qx(fHarmonic)*TMath::Sqrt(fMultNegTPC);
                fQnVecNegTPC[1] = theQnVector[2]->Qy(fHarmonic)*TMath::Sqrt(fMultNegTPC);
            break;
        }
        fQnVecNormNegTPC  = TMath::Sqrt(fQnVecNegTPC[0]*fQnVecNegTPC[0]+fQnVecNegTPC[1]*fQnVecNegTPC[1]);
    }

    //fill array of tracks needed if some tracks have to be removed
    fUsedTrackPosIDs.ResetAllBits();
    fUsedTrackNegIDs.ResetAllBits();
    
    int nTracks=fAODEvent->GetNumberOfTracks();
    for(int iTrack=0; iTrack<nTracks; iTrack++) {
        AliAODTrack* track=dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(iTrack));
        if(!track || !IsTrackSelected(track)) continue;
        if(track->GetID()>=0) fUsedTrackPosIDs.SetBitNumber(track->GetID());
        else fUsedTrackNegIDs.SetBitNumber(TMath::Abs(track->GetID()));
    }
}

//________________________________________________________________
void AliHFQnVectorHandler::ComputeQvecQnFrameworkV0() 
{
    if(!fQnVectorTask || !fQnVectorMgr) {
        AliWarning("Qn-framework task not found!");
        return;
    }
    
    TList* qnlist = dynamic_cast<TList*>(fQnVectorMgr->GetQnVectorList());
    if (!qnlist) return;  

    TString detV0ConfName[3] = {"VZERO","VZEROA","VZEROC"};
    TString normMethodName = "";
    switch(fNormMethod) {
        case 0: //kQoverQlength
            normMethodName = "QoverQlength";
        break;
        case 1: //kQoverM
            normMethodName = "QoverM";
        break;
        case 2: //kQoverSqrtM
            normMethodName = "QoverSqrtM";
        break;
    }
    AliQnCorrectionsQnVector *theQnVector[3] = {nullptr,nullptr,nullptr};

    for(int iDet=0; iDet<3; iDet++) {
        TList *pQvecList = dynamic_cast<TList*>(qnlist->FindObject(Form("%s%s",detV0ConfName[iDet].Data(),normMethodName.Data())));
        if (!pQvecList) return;  

        /* the detector is present */
        theQnVector[iDet] = dynamic_cast<AliQnCorrectionsQnVector*>(pQvecList->First());
        if (!theQnVector[iDet] || !(theQnVector[iDet]->IsGoodQuality()) || !(theQnVector[iDet]->GetN() != 0)) {
            /* the Qn vector for the expected step was not there */
	        theQnVector[iDet] = dynamic_cast<AliQnCorrectionsQnVector*>(pQvecList->FindObject("raw"));
        }
    }

    //fill Q vectors
    if(theQnVector[0]) {
        fMultFullV0 = theQnVector[0]->GetSumOfWeights();
        switch(fNormMethod) {
            case 0: //kQoverQlength
                fQnVecFullV0[0] = theQnVector[0]->Qx(fHarmonic)*theQnVector[0]->Length(fHarmonic);
                fQnVecFullV0[1] = theQnVector[0]->Qy(fHarmonic)*theQnVector[0]->Length(fHarmonic);
            break;
            case 1: //kQoverM
                fQnVecFullV0[0] = theQnVector[0]->Qx(fHarmonic)*fMultFullV0;
                fQnVecFullV0[1] = theQnVector[0]->Qy(fHarmonic)*fMultFullV0;
            break;
            case 2: //kQoverSqrtM
                fQnVecFullV0[0] = theQnVector[0]->Qx(fHarmonic)*TMath::Sqrt(fMultFullV0);
                fQnVecFullV0[1] = theQnVector[0]->Qy(fHarmonic)*TMath::Sqrt(fMultFullV0);
            break;
        }
        fQnVecNormFullV0 = TMath::Sqrt(fQnVecFullV0[0]*fQnVecFullV0[0]+fQnVecFullV0[1]*fQnVecFullV0[1]);
    }
    if(theQnVector[1]) {
        fMultV0A = theQnVector[1]->GetSumOfWeights();
        switch(fNormMethod) {
            case 0: //kQoverQlength
                fQnVecV0A[0] = theQnVector[1]->Qx(fHarmonic)*theQnVector[1]->Length(fHarmonic);
                fQnVecV0A[1] = theQnVector[1]->Qy(fHarmonic)*theQnVector[1]->Length(fHarmonic);
            break;
            case 1: //kQoverM
                fQnVecV0A[0] = theQnVector[1]->Qx(fHarmonic)*fMultV0A;
                fQnVecV0A[1] = theQnVector[1]->Qy(fHarmonic)*fMultV0A;
            break;
            case 2: //kQoverSqrtM
                fQnVecV0A[0] = theQnVector[1]->Qx(fHarmonic)*TMath::Sqrt(fMultV0A);
                fQnVecV0A[1] = theQnVector[1]->Qy(fHarmonic)*TMath::Sqrt(fMultV0A);
            break;
        }
        fQnVecNormV0A  = TMath::Sqrt(fQnVecV0A[0]*fQnVecV0A[0]+fQnVecV0A[1]*fQnVecV0A[1]);
    }
    if(theQnVector[2]) {
        fMultV0C = theQnVector[2]->GetSumOfWeights();
        switch(fNormMethod) {
            case 0: //kQoverQlength
                fQnVecV0C[0] = theQnVector[2]->Qx(fHarmonic)*theQnVector[2]->Length(fHarmonic);
                fQnVecV0C[1] = theQnVector[2]->Qy(fHarmonic)*theQnVector[2]->Length(fHarmonic);
            break;
            case 1: //kQoverM
                fQnVecV0C[0] = theQnVector[2]->Qx(fHarmonic)*fMultV0C;
                fQnVecV0C[1] = theQnVector[2]->Qy(fHarmonic)*fMultV0C;
            break;
            case 2: //kQoverSqrtM
                fQnVecV0C[0] = theQnVector[2]->Qx(fHarmonic)*TMath::Sqrt(fMultV0C);
                fQnVecV0C[1] = theQnVector[2]->Qy(fHarmonic)*TMath::Sqrt(fMultV0C);
            break;
        }
        fQnVecNormV0C  = TMath::Sqrt(fQnVecV0C[0]*fQnVecV0C[0]+fQnVecV0C[1]*fQnVecV0C[1]);
    }
}

//__________________________________________________________
void AliHFQnVectorHandler::ComputeQvecTPC() 
{
    short centbin = GetCentBin();

    //initialise Q vectors
    for(int iComp=0; iComp<2; iComp++) {
        fQnVecFullTPC[iComp] = 0.;
        fQnVecPosTPC[iComp]  = 0.;
        fQnVecNegTPC[iComp]  = 0.;
    }
    fMultPosTPC = 0., fMultFullTPC = 0., fMultNegTPC = 0.;
    fQnVecNormFullTPC = 1., fQnVecNormPosTPC  = 1., fQnVecNormNegTPC  = 1.;
    fUsedTrackPosIDs.ResetAllBits();
    fUsedTrackNegIDs.ResetAllBits();
    
    //reset phi distributions
    if(fEnablePhiDistrHistos) {
        fPhiVsCentrTPC[0]->Reset();
        fPhiVsCentrTPC[1]->Reset();
    }

    int nTracks=fAODEvent->GetNumberOfTracks();
    for(int iTrack=0; iTrack<nTracks; iTrack++) {
        AliAODTrack* track=dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(iTrack));
        double pt=track->Pt();
        double eta=track->Eta();
        if(!track || !IsTrackSelected(track)) continue;
        double phi=track->Phi();
        double pseudoRand = pt*1000.-(long)(pt*1000);
        if(pseudoRand>fFractionOfTracksForQnTPC) continue;
        if(track->GetID()>=0) fUsedTrackPosIDs.SetBitNumber(track->GetID());
        else fUsedTrackNegIDs.SetBitNumber(TMath::Abs(track->GetID()));
        double qx = TMath::Cos(fHarmonic*phi);
        double qy = TMath::Sin(fHarmonic*phi);
        double weight = 1.;
        if(eta>0) {
            if(fWeightsTPCPosEta[centbin]) {
                int phibin = fWeightsTPCPosEta[centbin]->GetXaxis()->FindBin(phi);
                weight = 1./fWeightsTPCPosEta[centbin]->GetBinContent(phibin);
            }
            fQnVecFullTPC[0] += weight*qx;
            fQnVecFullTPC[1] += weight*qy;
            fMultFullTPC     += weight;
            fQnVecPosTPC[0]  += weight*qx;
            fQnVecPosTPC[1]  += weight*qy;
            fMultPosTPC      += weight;
            if(fEnablePhiDistrHistos && fPhiVsCentrTPC[0])
                fPhiVsCentrTPC[0]->Fill(fCentrality,phi);
        }
        else {
            if(fWeightsTPCNegEta[centbin]) {
                int phibin = fWeightsTPCNegEta[centbin]->GetXaxis()->FindBin(phi);
                weight = 1./fWeightsTPCNegEta[centbin]->GetBinContent(phibin);
            }
            fQnVecFullTPC[0] += weight*qx;
            fQnVecFullTPC[1] += weight*qy;
            fMultFullTPC     += weight;
            fQnVecNegTPC[0]  += weight*qx;
            fQnVecNegTPC[1]  += weight*qy;
            fMultNegTPC      += weight;
            if(fEnablePhiDistrHistos && fPhiVsCentrTPC[1])
                fPhiVsCentrTPC[1]->Fill(fCentrality,phi);
        }
    }

    fQnVecNormFullTPC = TMath::Sqrt(fQnVecFullTPC[0]*fQnVecFullTPC[0]+fQnVecFullTPC[1]*fQnVecFullTPC[1]);
    fQnVecNormPosTPC  = TMath::Sqrt(fQnVecPosTPC[0]*fQnVecPosTPC[0]+fQnVecPosTPC[1]*fQnVecPosTPC[1]);
    fQnVecNormNegTPC  = TMath::Sqrt(fQnVecNegTPC[0]*fQnVecNegTPC[0]+fQnVecNegTPC[1]*fQnVecNegTPC[1]);
}

//__________________________________________________________
void AliHFQnVectorHandler::ComputeQvecV0()
{
    //initialise Q vectors
    for(int iComp=0; iComp<2; iComp++) {
        fQnVecFullV0[iComp] = 0.;
        fQnVecV0A[iComp]    = 0.;
        fQnVecV0C[iComp]    = 0.;
    }
    fMultFullV0 = 0., fMultV0A = 0., fMultV0C = 0.;
    fQnVecNormFullV0 = 1., fQnVecNormV0A  = 1., fQnVecNormV0C  = 1.;

    short zvtxbin = GetVertexZbin();

    for (int iCh = 0; iCh < 64; iCh++) {
        
        double phiCh = TMath::PiOver4()*(0.5 + iCh % 8);
        double multv0 = fV0->GetMultiplicity(iCh);
        
        if (iCh < 32) { // V0C side
            double multCorC = -10;
            
            if (iCh < 8)
                multCorC = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(1);
            else if (iCh >= 8 && iCh < 16)
                multCorC = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(9);
            else if (iCh >= 16 && iCh < 24)
                multCorC = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(17);
            else if (iCh >= 24 && iCh < 32)
                multCorC = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(25);
            
            if (multCorC < 0) {
                AliWarning("Problem with multiplicity in V0C");
                continue;
            }
            
            fQnVecV0C[0] += TMath::Cos(fHarmonic*phiCh) * multCorC;
            fQnVecV0C[1] += TMath::Sin(fHarmonic*phiCh) * multCorC;
            fQnVecFullV0[0] += TMath::Cos(fHarmonic*phiCh) * multCorC;
            fQnVecFullV0[1] += TMath::Sin(fHarmonic*phiCh) * multCorC;
 
            fMultV0C += multCorC;  
            fMultFullV0 += multCorC;  
        } 
        else { // V0A side
            double multCorA = -10;
            
            if (iCh >= 32 && iCh < 40)
                multCorA = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(33);
            else if (iCh >= 40 && iCh < 48)
                multCorA = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(41);
            else if (iCh >= 48 && iCh < 56)
                multCorA = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(49);
            else if (iCh >= 56 && iCh < 64)
                multCorA = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(57);
            
            if (multCorA < 0) {
                AliWarning("Problem with multiplicity in V0A");
                continue;
            }
            
            fQnVecV0A[0] += TMath::Cos(fHarmonic*phiCh) * multCorA;
            fQnVecV0A[1] += TMath::Sin(fHarmonic*phiCh) * multCorA;
            fQnVecFullV0[0] += TMath::Cos(fHarmonic*phiCh) * multCorA;
            fQnVecFullV0[1] += TMath::Sin(fHarmonic*phiCh) * multCorA;
 
            fMultV0A += multCorA;  
            fMultFullV0 += multCorA;              
        }
    }

    int iCentBin = static_cast<int>(fCentrality)+1;

    //only recentering and not width equalisation to preserve multiplicity dependence (needed for qn)
    fQnVecV0A[0] = (fQnVecV0A[0] - fQx2mV0A[zvtxbin]->GetBinContent(iCentBin));///fQx2sV0A[zvtxbin]->GetBinContent(iCentBin);
    fQnVecV0A[1] = (fQnVecV0A[1] - fQy2mV0A[zvtxbin]->GetBinContent(iCentBin));///fQy2sV0A[zvtxbin]->GetBinContent(iCentBin);
    fQnVecV0C[0] = (fQnVecV0C[0] - fQx2mV0C[zvtxbin]->GetBinContent(iCentBin));///fQx2sV0C[zvtxbin]->GetBinContent(iCentBin);   
    fQnVecV0C[1] = (fQnVecV0C[1] - fQy2mV0C[zvtxbin]->GetBinContent(iCentBin));///fQy2sV0C[zvtxbin]->GetBinContent(iCentBin);    

    fQnVecNormFullV0 = TMath::Sqrt(fQnVecFullV0[0]*fQnVecFullV0[0]+fQnVecFullV0[1]*fQnVecFullV0[1]);
    fQnVecNormV0A    = TMath::Sqrt(fQnVecV0A[0]*fQnVecV0A[0]+fQnVecV0A[1]*fQnVecV0A[1]);
    fQnVecNormV0C    = TMath::Sqrt(fQnVecV0C[0]*fQnVecV0C[0]+fQnVecV0C[1]*fQnVecV0C[1]);
}

//________________________________________________________________
bool AliHFQnVectorHandler::OpenInfoCalbration() 
{
    if(fOADBFile && fOADBFile->IsOpen()) {
        fOADBFile->Close();
        delete fOADBFile;
        fOADBFile = nullptr;
    }

    if (!gGrid) {
        TGrid::Connect("alien://");
    }
    fOADBFile = TFile::Open(fOADBFileName.Data());

    if(!fOADBFile) {
        AliWarning("OADB V0-TPC calibration file cannot be opened\n");
        return false;
    }
    
    //load V0 calibrations (mandatory)
    AliOADBContainer* cont = (AliOADBContainer*) fOADBFile->Get("hMultV0BefCorPfpx");
    if(!cont) {
        AliWarning("OADB object hMultV0BefCorPfpx is not available in the file\n");
        return false;
    }
    if(!(cont->GetObject(fRun))) {
        AliWarning(Form("OADB object hMultV0BefCorPfpx is not available for run %i\n", fRun));
        return false;
    }
    fHistMultV0 = ((TH1D*) cont->GetObject(fRun));
            
    for(int iZvtx = 0; iZvtx < 14; iZvtx++) {
        AliOADBContainer* contQx2am = (AliOADBContainer*) fOADBFile->Get(Form("fqxa%dm_%d", fHarmonic, iZvtx));
        if(!contQx2am) { //check if it is not Zvtx differential
            contQx2am = (AliOADBContainer*) fOADBFile->Get(Form("fqxa%dm", fHarmonic));
            if(contQx2am)
                fV0CalibZvtxDiff = false;
        }
        if(!contQx2am) {
            AliWarning(Form("OADB object fqxa%dm is not available in the file\n", fHarmonic));
            return false;
        }
        if(!(contQx2am->GetObject(fRun))) {
            AliWarning(Form("OADB object fqxa%dm is not available for run %i\n", fHarmonic, fRun));
            return false;
        }
        fQx2mV0A[iZvtx] = ((TH1D*) contQx2am->GetObject(fRun));
        
        AliOADBContainer* contQy2am = nullptr;
        if(fV0CalibZvtxDiff)
            contQy2am = (AliOADBContainer*) fOADBFile->Get(Form("fqya%dm_%d", fHarmonic, iZvtx));
        else
            contQy2am = (AliOADBContainer*) fOADBFile->Get(Form("fqya%dm", fHarmonic));
        if(!contQy2am) {
            AliWarning(Form("OADB object fqya%dm is not available in the file\n", fHarmonic));
            return false;
        }
        if(!(contQy2am->GetObject(fRun))) {
            AliWarning(Form("OADB object fqya%dm is not available for run %i\n", fHarmonic, fRun));
            return false;
        }
        fQy2mV0A[iZvtx] = ((TH1D*) contQy2am->GetObject(fRun));
        
        AliOADBContainer* contQx2as = nullptr;
        if(fV0CalibZvtxDiff)
            contQx2as = (AliOADBContainer*) fOADBFile->Get(Form("fqxa%ds_%d", fHarmonic, iZvtx));
        else
            contQx2as = (AliOADBContainer*) fOADBFile->Get(Form("fqxa%ds", fHarmonic));
        if(!contQx2as) {
            AliWarning(Form("OADB object fqxa%ds is not available in the file\n", fHarmonic));
            return false;
        }
        if(!(contQx2as->GetObject(fRun))) {
            AliWarning(Form("OADB object fqxa%ds is not available for run %i\n", fHarmonic, fRun));
            return false;
        }
        fQx2sV0A[iZvtx] = ((TH1D*) contQx2as->GetObject(fRun));
        
        AliOADBContainer* contQy2as = nullptr;
        if(fV0CalibZvtxDiff)
            contQy2as = (AliOADBContainer*) fOADBFile->Get(Form("fqya%ds_%d", fHarmonic, iZvtx));
        else
            contQy2as = (AliOADBContainer*) fOADBFile->Get(Form("fqya%ds", fHarmonic));
        if(!contQy2as) {
            AliWarning(Form("OADB object fqya%ds is not available in the file\n", fHarmonic));
            return false;
        }
        if(!(contQy2as->GetObject(fRun))) {
            AliWarning(Form("OADB object fqya%ds is not available for run %i\n", fHarmonic, fRun));
            return false;
        }
        fQy2sV0A[iZvtx] = ((TH1D*) contQy2as->GetObject(fRun));
        
        AliOADBContainer* contQx2cm = nullptr;
        if(fV0CalibZvtxDiff)
            contQx2cm = (AliOADBContainer*) fOADBFile->Get(Form("fqxc%dm_%d", fHarmonic, iZvtx));
        else
            contQx2cm = (AliOADBContainer*) fOADBFile->Get(Form("fqxc%dm", fHarmonic));
        if(!contQx2cm) {
            AliWarning(Form("OADB object fqxc%dm is not available in the file\n", fHarmonic));
            return false;
        }
        if(!(contQx2cm->GetObject(fRun))) {
            AliWarning(Form("OADB object fqxc%dm is not available for run %i\n", fHarmonic, fRun));
            return false;
        }
        fQx2mV0C[iZvtx] = ((TH1D*) contQx2cm->GetObject(fRun));
        
        AliOADBContainer* contQy2cm = nullptr;
        if(fV0CalibZvtxDiff)
            contQy2cm = (AliOADBContainer*) fOADBFile->Get(Form("fqyc%dm_%d", fHarmonic, iZvtx));
        else
            contQy2cm = (AliOADBContainer*) fOADBFile->Get(Form("fqyc%dm", fHarmonic));
        if(!contQy2cm) {
            AliWarning(Form("OADB object fqyc%dm is not available in the file\n", fHarmonic));
            return false;
        }
        if(!(contQy2cm->GetObject(fRun))) {
            AliWarning(Form("OADB object fqyc%dm is not available for run %i\n", fHarmonic, fRun));
            return false;
        }
        fQy2mV0C[iZvtx] = ((TH1D*) contQy2cm->GetObject(fRun));

        AliOADBContainer* contQx2cs = nullptr;
        if(fV0CalibZvtxDiff)
            contQx2cs = (AliOADBContainer*) fOADBFile->Get(Form("fqxc%ds_%d", fHarmonic, iZvtx));
        else
            contQx2cs = (AliOADBContainer*) fOADBFile->Get(Form("fqxc%ds", fHarmonic));
        if(!contQx2cs) {
            AliWarning(Form("OADB object fqxc%ds is not available in the file\n", fHarmonic));
            return false;
        }
        if(!(contQx2cs->GetObject(fRun))) {
            AliWarning(Form("OADB object fqxc%ds is not available for run %i\n", fHarmonic, fRun));
            return false;
        }
        fQx2sV0C[iZvtx] = ((TH1D*) contQx2cs->GetObject(fRun));
        
        AliOADBContainer* contQy2cs = nullptr;
        if(fV0CalibZvtxDiff)
            contQy2cs = (AliOADBContainer*) fOADBFile->Get(Form("fqyc%ds_%d", fHarmonic, iZvtx));
        else
            contQy2cs = (AliOADBContainer*) fOADBFile->Get(Form("fqyc%ds", fHarmonic));
        if(!contQy2cs) {
            AliWarning(Form("OADB object fqyc%ds is not available in the file\n", fHarmonic));
            return false;
        }
        if(!(contQy2cs->GetObject(fRun))) {
            AliWarning(Form("OADB object fqyc%ds is not available for run %i\n", fHarmonic, fRun));
            return false;
        }
        fQy2sV0C[iZvtx] = ((TH1D*) contQy2cs->GetObject(fRun));

        if(!fV0CalibZvtxDiff) //assign only first element of array if it is not Zvtx differential
            break;
    }

    //load TPC calibrations (not mandatory)
    for(int iCent = 0; iCent < 9; iCent++) {
        AliOADBContainer* contTPCposEta = (AliOADBContainer*) fOADBFile->Get(Form("fphidistr_poseta_%d_%d", iCent*10, (iCent+1)*10));
        if(!contTPCposEta) {
            AliWarning("OADB object fphidistr_poseta is not available in the file\n");
            fWeightsTPCPosEta[iCent] = nullptr;
        }
        else {
            if(!(contTPCposEta->GetObject(fRun))) {
                AliWarning(Form("OADB object fphidistr_poseta is not available for run %i\n", fRun));
                fWeightsTPCPosEta[iCent] = nullptr;
            }
            else {
                fWeightsTPCPosEta[iCent] = ((TH1D*) contTPCposEta->GetObject(fRun));   
            }
        }

        AliOADBContainer* contTPCnegEta = (AliOADBContainer*) fOADBFile->Get(Form("fphidistr_negeta_%d_%d", iCent*10, (iCent+1)*10));
        if(!contTPCnegEta) {
            AliWarning("OADB object fphidistr_negeta is not available in the file\n");
            fWeightsTPCNegEta[iCent] = nullptr;
            return true;
        }
        else {        
            if(!(contTPCnegEta->GetObject(fRun))) {
                AliWarning(Form("OADB object fphidistr_negeta is not available for run %i\n", fRun));
                fWeightsTPCNegEta[iCent] = nullptr;
            }
            else {
                fWeightsTPCNegEta[iCent] = ((TH1D*) contTPCnegEta->GetObject(fRun));   
            }
        }
    }

    return true;
}

//__________________________________________________________
short AliHFQnVectorHandler::GetVertexZbin() const
{
    if(!fV0CalibZvtxDiff)
        return 0; //if it is not Zvtx differential, always first bin

    short zvtxbin = -10;
    
    if (fZvtx >= -10. && fZvtx < -8.)
        zvtxbin = 0;
    else if (fZvtx >= -8. && fZvtx < -6.)
        zvtxbin = 1;
    else if (fZvtx >= -6. && fZvtx < -4.)
        zvtxbin = 2;
    else if (fZvtx >= -4. && fZvtx < -3.)
        zvtxbin = 3;
    else if (fZvtx >= -3. && fZvtx < -2.)
        zvtxbin = 4;
    else if (fZvtx >= -2. && fZvtx < -1.)
        zvtxbin = 5;
    else if (fZvtx >= -1. && fZvtx < 0)
        zvtxbin = 6;
    else if (fZvtx >= 0 && fZvtx < 1.)
        zvtxbin = 7;
    else if (fZvtx >= 1. && fZvtx < 2.)
        zvtxbin = 8;
    else if (fZvtx >= 2. && fZvtx < 3.)
        zvtxbin = 9;
    else if (fZvtx >= 3. && fZvtx < 4.)
        zvtxbin = 10;
    else if (fZvtx >= 4. && fZvtx < 6.)
        zvtxbin = 11;
    else if (fZvtx >= 6. && fZvtx < 8.)
        zvtxbin = 12;
    else if (fZvtx >= 8. && fZvtx <= 10.)
        zvtxbin = 13;
    
    return zvtxbin;
}

//__________________________________________________________
short AliHFQnVectorHandler::GetCentBin() const
{
    short centbin = -10;
    
    if (fCentrality >= 0. && fCentrality < 10.)
        centbin = 0;
    else if (fCentrality >= 10. && fCentrality < 20.)
        centbin = 1;
    else if (fCentrality >= 20. && fCentrality < 30.)
        centbin = 2;
    else if (fCentrality >= 30. && fCentrality < 40.)
        centbin = 3;
    else if (fCentrality >= 40. && fCentrality < 50.)
        centbin = 4;
    else if (fCentrality >= 50. && fCentrality < 60.)
        centbin = 5;
    else if (fCentrality >= 60. && fCentrality < 70.)
        centbin = 6;
    else if (fCentrality >= 70. && fCentrality < 80.)
        centbin = 7;
    else if (fCentrality >= 80. && fCentrality < 90.)
        centbin = 8;
    else if(fCentrality >= 90.)
        centbin = 8;

    return centbin;
}

//__________________________________________________________
bool AliHFQnVectorHandler::IsTrackSelected(AliAODTrack* track) {

    if(fRun>=244918 && fRun<=246994) {//PbPb2015
        if(fCalibType==kQnCalib) {
            if(!track->TestFilterBit(BIT(8)) && !track->TestFilterBit(BIT(9)))
                return false;
            if(track->Pt()<fPtMinTPC || track->Pt()>fPtMaxTPC || TMath::Abs(track->Eta())>fEtaMaxTPC || TMath::Abs(track->Eta())<fEtaMinTPC)
                return false;
        }
        else if(fCalibType==kQnFrameworkCalib) {
            if(!track->TestFilterBit(BIT(8)) && !track->TestFilterBit(BIT(9)))
                return false;
            if(track->Pt()<0.2 || track->Pt()>5 || TMath::Abs(track->Eta())>0.8 || TMath::Abs(track->Eta())<0.)
                return false;
        }
    }
    else if(fRun>=295581 && fRun<=297317) { //PbPb2018
        if(fCalibType==kQnCalib) {
            if(!track->TestFilterBit(BIT(8)) && !track->TestFilterBit(BIT(9)))
                return false;
            if(track->Pt()<fPtMinTPC || track->Pt()>fPtMaxTPC || TMath::Abs(track->Eta())>fEtaMaxTPC || TMath::Abs(track->Eta())<fEtaMinTPC)
                return false;
        }
        else if(fCalibType==kQnFrameworkCalib) {
            if(!track->TestFilterBit(BIT(8)) && !track->TestFilterBit(BIT(9)))
                return false;
            if(track->Pt()<0.2 || track->Pt()>5 || TMath::Abs(track->Eta())>0.8 || TMath::Abs(track->Eta())<0.)
                return false;
        }
    }
    else { //default
        if(fCalibType==kQnCalib) {
            if(!track->TestFilterBit(BIT(8)) && !track->TestFilterBit(BIT(9)))
                return false;
            if(track->Pt()<fPtMinTPC || track->Pt()>fPtMaxTPC || TMath::Abs(track->Eta())>fEtaMaxTPC || TMath::Abs(track->Eta())<fEtaMinTPC)
                return false;
        }
        else if(fCalibType==kQnFrameworkCalib) {
            if(!track->TestFilterBit(BIT(8)) && !track->TestFilterBit(BIT(9)))
                return false;
            if(track->Pt()<0.2 || track->Pt()>5 || TMath::Abs(track->Eta())>0.8 || TMath::Abs(track->Eta())<0.)
                return false;
        }
    }

    return true;
}