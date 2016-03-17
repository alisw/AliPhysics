/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, proviyaded that the above copyright notice appears in all *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purapose. It is         *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This analysis extracts information into TTree for pT-spectra of deuterons //
// Based on AliAnalysisDeuteronpA task of J. Anielski for deuteron analysis  //
// and AliAnalysisTaskExtractV0 by D. Chinellato for TTree interface         //
// L.Barnby October 2015                                                     //
///////////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "TH1F.h"
#include "TList.h"
#include "TChain.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliCentrality.h"
#include "AliESDUtils.h"
#include "AliMultiplicity.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliAnalysisUtils.h"

#include "AliAnalysisDeuteronTree.h"

ClassImp(AliAnalysisDeuteronTree)

//________________________________________________________________________
AliAnalysisDeuteronTree::AliAnalysisDeuteronTree()
: AliAnalysisTaskSE(), fListHist(0), fTree(0),fESDtrackCuts(0),fPIDResponse(0),
fMCtrue(0),
fRapCMSpA(0),
fUtils(0),
fTimeStamp(0),
fCentrality(0),
fZVertex(0),
fCL1(0),
fPtCor(0),
fPt(0),
fMom(0),
fRapd(0),
fPxd(0),
fPyd(0),
fPzd(0),
fNsigmaTPCd(0),
fNsigmaTOFd(0),
fDcaXYd(0),
fMcCode(0),
fNpion(0),
fhZVertex(0),
fhCentrality(0),
fhNsigmaTPCvsMom(0),
fhNsigmaTOFvsMom(0),
fMomCorrConstA(0),
fMomCorrConstB(0),
fMomCorrPower(0),
fMinPtCut(0),
fMinPtTOFCut(0),
fNsigmaTOFdMax(0),
fNsigmaTOFdMin(0),
fNsigmaTPCdMax(0),
fNsigmaTPCdMin(0),
fNsigmaTPCpiAbsCut(0)
{
    //Default constructor
}

//________________________________________________________________________
AliAnalysisDeuteronTree::AliAnalysisDeuteronTree(const char *name)
: AliAnalysisTaskSE(name), fListHist(0), fTree(0), fESDtrackCuts(0),fPIDResponse(0),
fMCtrue(0),
fRapCMSpA(0),
fUtils(0),
fTimeStamp(0),
fCentrality(0),
fZVertex(0),
fCL1(0),
fPtCor(0),
fPt(0),
fMom(0),
fRapd(0),
fPxd(0),
fPyd(0),
fPzd(0),
fNsigmaTPCd(0),
fNsigmaTOFd(0),
fDcaXYd(0),
fMcCode(0),
fNpion(0),
fhZVertex(0),
fhCentrality(0),
fhNsigmaTPCvsMom(0),
fhNsigmaTOFvsMom(0),
fMomCorrConstA(0),
fMomCorrConstB(0),
fMomCorrPower(0),
fMinPtCut(0),
fMinPtTOFCut(0),
fNsigmaTOFdMax(0),
fNsigmaTOFdMin(0),
fNsigmaTPCdMax(0),
fNsigmaTPCdMin(0),
fNsigmaTPCpiAbsCut(0)
{
    //Standard constructor
    fMCtrue = kFALSE;
    fRapCMSpA = kTRUE;

    // The array size would be better as a const Int_t type but this wouldn't compile
    for (Int_t t=0; t<400; t++) {
        fPx[t]=0;
        fPy[t]=0;
        fPz[t]=0;
        fNsigmaTPCpi[t]=0;
        fCharge[t]=0;
    }
    //fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
    Initialize();
    
    //AliInfo("About to define TChain input");
    //DefineInput(0, TChain::Class());
    
    AliInfo("About to define TList output");
    DefineOutput(1, TList::Class());

    AliInfo("About to define TTree output");
    DefineOutput(2, TTree::Class());
    
    AliInfo("Constructor Finished");
}

//________________________________________________________________________
void AliAnalysisDeuteronTree::Initialize()
{
    //
    // updating parameters in case of changes
    //
    AliInfo("Initialization started");
    
    fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,kTRUE);
    if (fESDtrackCuts == 0x0) {
        AliWarning("No ESDtrackCuts");
    }
    fESDtrackCuts->SetMaxDCAToVertexXY(3);
    fESDtrackCuts->SetMaxDCAToVertexZ(2);
    fESDtrackCuts->SetEtaRange(-0.9,0.9);
    //
    
    // Momentum correction
    fMomCorrConstA = 0.333303;
    fMomCorrConstB = 0.651111;
    fMomCorrPower = 5.27268;
    
    fMinPtCut = 0.5;
    fMinPtTOFCut = 1.0;
    fNsigmaTOFdMax = 5.0;
    fNsigmaTOFdMin = -5.0;
    fNsigmaTPCdMax = 5.0;
    fNsigmaTPCdMin = -5.0;

    fNsigmaTPCpiAbsCut = 3.0;

    AliInfo("Initialization complete");
}

//________________________________________________________________________
AliAnalysisDeuteronTree::~AliAnalysisDeuteronTree()
{
    if (fListHist) {
        delete fListHist;
        fListHist = 0x0;
    }
    if (fTree) {
        delete fTree;
        fTree = 0x0;
    }
    //cleanup esd track cuts object too...
    if (fESDtrackCuts) {
        delete fESDtrackCuts;
        fESDtrackCuts = 0x0;
    }
    if (fPIDResponse) {
        delete fPIDResponse;
        fPIDResponse = 0x0;
    }
    if (fUtils) {
        delete fUtils;
        fUtils = 0x0;
    }

}

//________________________________________________________________________
void AliAnalysisDeuteronTree::UserCreateOutputObjects()
{
    AliInfo("Creation started");
    
    // Proper handling of output objects ("magic" lines)
    // For file-resident Tree
    OpenFile(2);
    fTree = new TTree("fTree","dCandidates");

    fListHist = new TList();
    fListHist->SetOwner();

    // Get PID response object
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    if(!man)
        AliFatal("Could not find manager");
    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*> (man->GetInputEventHandler());
    if(!inputHandler)
        AliFatal("No input event handler");
        fPIDResponse = dynamic_cast<AliPIDResponse *>(inputHandler->GetPIDResponse());
    if (!fPIDResponse)
        AliFatal("PIDResponse object was not created"); // Escalated to fatal. Task is unusable without PID response.
    
    //Create analysis utils for event selection and pileup rejection
    fUtils = new AliAnalysisUtils();
    fUtils->SetCutOnZVertexSPD(kFALSE);

    
    // histograms
    fhZVertex = new TH1F("hZVertex","Vertex Z position;Z_{VTX} (cm);Counts/0.25 cm",100,-25,25);
    fListHist->Add(fhZVertex);
    fhCentrality = new TH1F("hCentrality",";Centrality(V0A) %;Counts/1%",100,0,100);
    fListHist->Add(fhCentrality);

    fhNsigmaTPCvsMom = new TH2F("hNsigmaTPCvsMom","n_{#sigma} TPC distribution;p (GeV/#it{c});n_{#sigma}^deuteron",300,0.0,6.0,200,-15.0,5);
    fListHist->Add(fhNsigmaTPCvsMom);

    fhNsigmaTOFvsMom = new TH2F("hNsigmaTOFvsMom","n_{#sigma} TOF distribution;p (GeV/#it{c});n_{#sigma}^deuteron",300,0.0,6.0,200,-15.0,5);
    fListHist->Add(fhNsigmaTOFvsMom);

    // fTree Branch definitions
    fTree->Branch("fCentrality",&fCentrality,"fCentrality/F");
    fTree->Branch("fCL1",&fCL1,"fCL1/F");
    fTree->Branch("fZVertex",&fZVertex,"fZVertex/F");
    fTree->Branch("fTimeStamp",&fTimeStamp,"fTimeStamp/F");
    fTree->Branch("fPtCor",&fPtCor,"fPtCor/F");
    fTree->Branch("fPt",&fPt,"fPt/F");
    fTree->Branch("fMom",&fMom,"fMom/F");
    fTree->Branch("fRapd",&fRapd,"fRapd/F");
    fTree->Branch("fPxd",&fPxd,"fPxd/F");
    fTree->Branch("fPyd",&fPyd,"fPyd/F");
    fTree->Branch("fPzd",&fPzd,"fPzd/F");
    fTree->Branch("fNsigmaTPCd",&fNsigmaTPCd,"fNsigmaTPCd/F");
    fTree->Branch("fNsigmaTOFd",&fNsigmaTOFd,"fNsigmaTOFd/F");
    fTree->Branch("fDcaXYd",&fDcaXYd,"fDcaXYd/F");
    fTree->Branch("fMcCode",&fMcCode,"fMcCode/I");
    fTree->Branch("fNpion",&fNpion,"fNpion/I");
    // Pion branches
    fTree->Branch("fPx",&fPx,"fPx[fNpion]/F");
    fTree->Branch("fPy",&fPy,"fPy[fNpion]/F");
    fTree->Branch("fPz",&fPz,"fPz[fNpion]/F");
    fTree->Branch("fNsigmaTPCpi",&fNsigmaTPCpi,"fNsigmaTPCpi[fNpion]/F");
    fTree->Branch("fCharge",&fCharge,"fCharge[fNpion]/I");

    PostData(1, fListHist);
    PostData(2, fTree);

    AliInfo("Creation finished");
}

//________________________________________________________________________
void AliAnalysisDeuteronTree::UserExec(Option_t *option){

    AliInfo("UserExec started");

    AliESDEvent *lESDevent = 0x0;
    lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!lESDevent) {
        AliWarning("ERROR: lESDevent not available \n");
        return;
    }
    
    
    AliCentrality *centrality = lESDevent->GetCentrality();
    
    // Later add code for accessing MC information
    // AliMCEventHandler, AliMCEvent, AliStack
    
    // Event checks
    //first event in chunck --> continue
    if(fUtils->IsFirstEventInChunk(lESDevent)) {
        PostData(1, fListHist);
        PostData(2, fTree);
        return;
    }
    //check for pileup
    if (fUtils->IsPileUpEvent(lESDevent)){
        PostData(1, fListHist);
        PostData(2, fTree);
       return;
    }
    //vertex cuts
    Bool_t isVtxOk = kTRUE;
    if(!fUtils->IsVertexSelected2013pA(lESDevent)) isVtxOk = kFALSE;

    if (!fESDtrackCuts) {
        AliFatal("ERROR: fESDtrackCuts not available"); // No sense to continue without track cuts
        return;
    }
    
    // Physics selection should be taken care of by analysis manager
    // and task->SetCollisionsCandidates(AliVEvent::kMB) in the macro
    
    const AliESDVertex *vertex = lESDevent->GetPrimaryVertex();
    if (vertex && isVtxOk) fhZVertex->Fill(vertex->GetZ());
   
    fCentrality = centrality->GetCentralityPercentile("V0A");
    fCL1 = centrality->GetCentralityPercentile("CL1");
    fZVertex = vertex->GetZ();
    fhCentrality->Fill(fCentrality);
    
    fTimeStamp = lESDevent->GetTimeStamp();
    
    Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
    Double_t lTOFNsigma, lTPCNsigma, lEnergyDeuteron, lPz;
    Double_t p[3], pd[3];
    Int_t lNpion;
    
    for (Int_t i=0;i<lESDevent->GetNumberOfTracks();++i) {
        //
        AliESDtrack *track = lESDevent->GetTrack(i);
        
        if (!track->GetInnerParam()) continue;
        
        Double_t ptot = track->GetInnerParam()->GetP(); // momentum for dEdx determination
        fMom=ptot*track->Charge(); // Shouldn't this be corrected for the different mass assumption as below
        // i.e find uncorrected pz, correct pt and then recombine to get corrected ptot?
        
        //
        // momentum correction for different mass assumption in tracking
        //
        fPt = track->Pt();
        fPtCor = track->Pt()/(1 - fMomCorrConstA/TMath::Power(track->Pt() + fMomCorrConstB, fMomCorrPower));
        
        track->GetImpactParameters(dca, cov);
        if (!fESDtrackCuts->AcceptTrack(track)) continue;
        if (track->Pt()< fMinPtCut) continue;
        
        track->GetPxPyPz(pd);
        fPxd=pd[0];
        fPyd=pd[1];
        fPzd=pd[2];
        
        AliPIDResponse::EDetPidStatus lETOFStatus, lETPCStatus;
        lETOFStatus = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, AliPID::kDeuteron, lTOFNsigma);
        fNsigmaTOFd = -999;
        if(lETOFStatus==AliPIDResponse::kDetPidOk) fNsigmaTOFd = lTOFNsigma;
        
        lETPCStatus = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kDeuteron, lTPCNsigma);
        fNsigmaTPCd = -999;
        if(lETPCStatus==AliPIDResponse::kDetPidOk) fNsigmaTPCd = lTPCNsigma;
        
        fDcaXYd = dca[0];
        
        lEnergyDeuteron = TMath::Sqrt(fPt*fPt + fPzd*fPzd +
                                      AliPID::ParticleMass(AliPID::kDeuteron)*AliPID::ParticleMass(AliPID::kDeuteron));
        fRapd = 0.5*TMath::Log((lEnergyDeuteron + fPzd)/(lEnergyDeuteron - fPzd));
        fMcCode = 0;
        
        fhNsigmaTPCvsMom->Fill(ptot,fNsigmaTPCd);
        fhNsigmaTOFvsMom->Fill(track->GetP(),fNsigmaTOFd);
        
        // Need some selections here to stop the tree getting too big (because of all the non-deuterons)
        // outside TPC dEd/dx band
        if (fNsigmaTPCd < fNsigmaTPCdMin || fNsigmaTPCd > fNsigmaTPCdMax) continue;
        // over TOF pt threshold and outside TOF band
        if ((fPt > fMinPtTOFCut && (fNsigmaTOFd < fNsigmaTOFdMin || fNsigmaTOFd > fNsigmaTOFdMax))) continue;

        // Only in this case loop to get pions
        lNpion = 0;
        for (Int_t j=0; j<lESDevent->GetNumberOfTracks();++j) {
            if (i!=j) { // The deuteron candidate should not be used as the pion
                AliESDtrack *pion = lESDevent->GetTrack(j);
                if (!fESDtrackCuts->AcceptTrack(pion)) continue;
                if (!pion->GetInnerParam()) continue;
                lETPCStatus = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, pion, AliPID::kPion, lTPCNsigma);
                if (TMath::Abs(lTPCNsigma)<fNsigmaTPCpiAbsCut) {
                    pion->GetPxPyPz(p);
                    fPx[lNpion] = p[0];
                    fPy[lNpion] = p[1];
                    fPz[lNpion] = p[2];
                    fNsigmaTPCpi[lNpion] = lTPCNsigma;
                    fCharge[lNpion] = pion->Charge();
                    lNpion++;
                }
            }

        }
        fNpion = lNpion;
        fTree->Fill();
        
        
    } // End track loop
    
    // Output data
    PostData(1, fListHist);
    PostData(2, fTree);
    return;
 
}

void AliAnalysisDeuteronTree::Terminate(Option_t *)
{
    // Draw control histograms (if any)
}