/**************************************************************************
 * Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/////////////////////////////////////////////////////////////
// Author: Tiantian Cheng (1,2)
// (1) Central China Normal University
// (2) GSI Helmholtz Centre for Heavy Ion Research
// E-mail: chengtiantian@mails.ccnu.edu.cn
/////////////////////////////////////////////////////////////

#include <TObjString.h>
#include <TSystem.h>
#include <iostream>
#include <iomanip>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TTree.h>
#include "TChain.h"
#include "AliRDHFCuts.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliExternalTrackParam.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSESemileptonicOmegac0KFP.h"
#include "AliTOFPIDResponse.h"
#include "AliAODPidHF.h"
#include "AliInputEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliCentrality.h"
#include "AliVertexerTracks.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"

// includes added to play with KFParticle
#ifndef HomogeneousField
#define HomogeneousField
#endif

using std::cout;
using std::endl;

class  AliAnalysisTaskSESemileptonicOmegac0KFP; // your analysis class

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSESemileptonicOmegac0KFP)  // classimp: necessary for root
/// \endcond

//------------------------------------------------------------------------------------------
AliAnalysisTaskSESemileptonicOmegac0KFP :: AliAnalysisTaskSESemileptonicOmegac0KFP():
AliAnalysisTaskSE(),
fAnalCuts(0),
fPID(0),
fpVtx(0),
fMCEvent(0),
fBzkG(0),
fOutputList(0),
fTree_Event(0),
fVar_Event(0),
fTree_Omegac0(0),
fVar_Omegac0(0),
fTree_Omegac0_QA(0),
fVar_Omegac0_QA(0),
fTree_Omegac0MCGen(0),
fVar_Omegac0MCGen(0),
fListCuts(0),
fUseMCInfo(kFALSE),
fCounter(0),
fHistEvents(0),
fHTrigger(0),
fHistoElectronTPCPID(0),
fHistoElectronTOFPID(0),
fHistoElectronTPCPIDSelTOF(0),
fHistoMassConversions(0),
fWriteOmegac0Tree(kFALSE),
fWriteOmegac0QATree(kFALSE),
fWriteOmegac0MCGenTree(kFALSE)

{
 //
 // Default Constructor
 //
}
//------------------------------------------------------------------------------------------------
AliAnalysisTaskSESemileptonicOmegac0KFP :: AliAnalysisTaskSESemileptonicOmegac0KFP(const char * name, AliRDHFCutsOmegactoeleOmegafromKFP*   analCuts):

AliAnalysisTaskSE(name),
fAnalCuts(analCuts),
fPID(0),
fpVtx(0),
fMCEvent(0),
fBzkG(0),
fOutputList(0),
fTree_Event(0),
fVar_Event(0),
fTree_Omegac0(0),
fVar_Omegac0(0),
fTree_Omegac0_QA(0),
fVar_Omegac0_QA(0),
fTree_Omegac0MCGen(0),
fVar_Omegac0MCGen(0),
fListCuts(0),
fUseMCInfo(kFALSE),
fCounter(0),
fHistEvents(0),
fHTrigger(0),
fHistoElectronTPCPID(0),
fHistoElectronTOFPID(0),
fHistoElectronTPCPIDSelTOF(0),
fHistoMassConversions(0),
fWriteOmegac0Tree(kFALSE),
fWriteOmegac0QATree(kFALSE),
fWriteOmegac0MCGenTree(kFALSE)

{
    //
    // Constructor. Initialization of Inputs and Outputs
    //
    
    DefineInput(0,  TChain::Class()); // define the input of the analysis: in this case
    DefineOutput(1, TList::Class());  // define the output of the analysis: in this case it is a list of histograms
    DefineOutput(2, AliNormalizationCounter::Class());
    DefineOutput(3,TList::Class()); // analysis histogram
    DefineOutput(4, TTree::Class()); // event
    DefineOutput(5, TTree::Class()); //Omegac0
    DefineOutput(6, TTree::Class()); //Omegac0 MC Gen
    DefineOutput(7, TTree::Class()); //Omegac0 QA Tree

}
//------------------------------------------------------------------------------------------------
AliAnalysisTaskSESemileptonicOmegac0KFP :: ~ AliAnalysisTaskSESemileptonicOmegac0KFP()
{
    
    // destructor
    if (fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
        fOutputList = 0;
    }
    
    if (fListCuts) {
        delete fListCuts;
        fListCuts = 0;
    }
    
    if (fAnalCuts) {
        delete fAnalCuts;
        fAnalCuts = 0;
    }
    
    if (fTree_Event) {
        delete fTree_Event;
        fTree_Event = 0;
    }
    
    if (fVar_Event) {
        delete fVar_Event;
        fVar_Event = 0;
    }
    
    if ( fTree_Omegac0) {
        delete  fTree_Omegac0;
         fTree_Omegac0 = 0;
    }
    
    if (fVar_Omegac0) {
        delete fVar_Omegac0;
        fVar_Omegac0 = 0;
    }
    
    if (fTree_Omegac0MCGen) {
        delete fTree_Omegac0MCGen;
        fTree_Omegac0MCGen = 0;
    }
    
    if (fVar_Omegac0MCGen) {
        delete fVar_Omegac0MCGen;
        fVar_Omegac0MCGen = 0;
    }
    
    if (fCounter) {
        delete fCounter;
        fCounter = 0;
    }
    
    if ( fTree_Omegac0_QA) {
        delete  fTree_Omegac0_QA;
         fTree_Omegac0_QA = 0;
    }
    
    if (fVar_Omegac0_QA) {
        delete fVar_Omegac0_QA;
        fVar_Omegac0_QA = 0;
    }
   
}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskSESemileptonicOmegac0KFP :: Init()
{
    // Initialization
    
    if(fDebug >1 ) AliInfo("Init");
    
    fListCuts = new TList();
    fListCuts -> SetOwner();
    fListCuts -> SetName("ListCuts");
    fListCuts -> Add(new AliRDHFCutsOmegactoeleOmegafromKFP(*fAnalCuts));
    PostData(1, fListCuts);
    
    return;
}
//----------------------------------------------------------------------------------------------
void AliAnalysisTaskSESemileptonicOmegac0KFP :: UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use
   
    // Counter for Normalization
    TString normName="NormalizationCounter";
    AliAnalysisDataContainer *cont = GetOutputSlot(2)->GetContainer();
    if(cont) normName = (TString)cont->GetName();
    fCounter = new AliNormalizationCounter(normName.Data());
    fCounter->Init();
    PostData(2, fCounter);
    
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    fOutputList->SetName("AnaHist");
    DefineAnaHist(); // define analysis histograms
    PostData(3,fOutputList);
    
    DefineEvent();
    PostData(4,fTree_Event);
    
    DefineTreeRecoOmegac0();
    PostData(5,fTree_Omegac0);
    
    DefineTreeMCGenOmegac0();
    PostData(6,fTree_Omegac0MCGen);
    
    DefineTreeRecoOmegac0_QA();
    PostData(7,fTree_Omegac0_QA); // used for QA check
    
    return;
    
    // fOutputList object. the manager will in the end take care of writing your output to file
    // so it needs to know what's in the output
}
//----------------------------------------------------------------------------------------------
void AliAnalysisTaskSESemileptonicOmegac0KFP :: UserExec(Option_t *)
{
    
    // user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you
    // have access to the current event.
    // once you return from the UserExec function, the manager will retrieve the next event from the chain

    
    if (!fInputEvent) { // if the event is empty (getting it failed) skip this event
        AliError("NO EVENT FOUND!");
        return;
    }
    AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
    fHistEvents->Fill(1);
    
    //--------------------------------------------------------------
    // First check if the event has magnetic field and proper vertex
    //--------------------------------------------------------------
    
    fBzkG = (Double_t)aodEvent->GetMagneticField();
    if (TMath::Abs(fBzkG)<0.001) return;
    // Setting magnetic field for KF vertexing
    KFParticle::SetField(fBzkG);
    
    fpVtx = (AliAODVertex*)aodEvent->GetPrimaryVertex();
    if (!fpVtx) return;
    fHistEvents->Fill(2);
    
    fCounter->StoreEvent(aodEvent,fAnalCuts,fUseMCInfo);
 
    //------------------------------------------------
    // MC analysis setting
    //------------------------------------------------
    
    TClonesArray *mcArray = 0;
    AliAODMCHeader *mcHeader = 0;
    
    if (fUseMCInfo) {
        fMCEvent = MCEvent(); // get the corresponding MC event fMCEvent
        if (!fMCEvent) {
            Printf("ERROR: Could not retrieve MC event");
            return;
        }
        
    // MC array need for matching
    mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle ::StdBranchName()));
    if (!mcArray){
        AliError("Could not find Monte-Carlo in AOD");
    }
     fHistEvents->Fill(7); // number of MC array exist
        
    // load MC header
    mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if ( !mcHeader ) {
        AliError("AliAnalysisTaskSESemileptonicOmegac0KFP::UserExec: MC header branch not found!\n");
        return;
    }
    fHistEvents->Fill(8); // number of MC header exist

    Double_t zMCvtx = mcHeader->GetVtxZ();
    if ( TMath::Abs(zMCvtx) > fAnalCuts->GetMaxVtxZ() ) {
        AliDebug(2,Form("Event rejected: fabs(zVtxMC)=%f > fAnalCuts->GetMaxVtxZ()=%f", zMCvtx, fAnalCuts->GetMaxVtxZ()));
        return;
    } else {
        fHistEvents->Fill(18);
    }
    if ((TMath::Abs(zMCvtx) < fAnalCuts->GetMaxVtxZ()) && (!fAnalCuts->IsEventRejectedDuePhysicsSelection()) && (!fAnalCuts->IsEventRejectedDueToTrigger())) {
        Bool_t selevt = MakeMCAnalysis(mcArray);
        if(!selevt) return;
    }
}
    //------------------------------------------------
    // Event selection
    //------------------------------------------------
  
    Bool_t IsTriggerNotOK = fAnalCuts->IsEventRejectedDueToTrigger();
    Bool_t IsPhysSelNotOK = fAnalCuts->IsEventRejectedDuePhysicsSelection();
    Bool_t IsNoVertex = fAnalCuts->IsEventRejectedDueToNotRecoVertex();
    if( !IsTriggerNotOK && !IsPhysSelNotOK && !IsNoVertex && fabs(fpVtx->GetZ())<fAnalCuts->GetMaxVtxZ() ) fHistEvents->Fill(3);
    
    Bool_t IsEventSelected = fAnalCuts->IsEventSelected(aodEvent);
    if(!IsEventSelected) {
        return;
    }
    fHistEvents->Fill(4);
    
    Bool_t IsMB = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB)==(AliVEvent::kMB);
    Bool_t IsSemi = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kSemiCentral)==(AliVEvent::kSemiCentral);
    Bool_t IsCent = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kCentral)==(AliVEvent::kCentral);
    Bool_t IsINT7 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kINT7)==(AliVEvent::kINT7);
    Bool_t IsEMC7 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kEMC7)==(AliVEvent::kEMC7);
    if(IsMB) fHTrigger->Fill(1);
    if(IsSemi) fHTrigger->Fill(2);
    if(IsCent) fHTrigger->Fill(3);
    if(IsINT7) fHTrigger->Fill(4);
    if(IsEMC7) fHTrigger->Fill(5);
    if(IsMB||IsSemi||IsCent) fHTrigger->Fill(7);
    if(IsINT7||IsEMC7) fHTrigger->Fill(8);
    if(IsMB&&IsSemi) fHTrigger->Fill(10);
    if(IsMB&&IsCent) fHTrigger->Fill(11);
    if(IsINT7&&IsEMC7) fHTrigger->Fill(12);
    
    // set primary vertex
    
    //------------------------------------------------
    // Check if the event has v0 candidate
    //------------------------------------------------
     Int_t num_v0 = aodEvent->GetNumberOfV0s();
     if (num_v0>0)
     fHistEvents->Fill(5);
    
    //------------------------------------------------
    // Check if the event has cascade candidate
    //------------------------------------------------
    Int_t num_casc = aodEvent->GetNumberOfCascades();
    fHistEvents->Fill(6);
    
    // set primary vertex
    KFPVertex pVertex;
    Double_t pos[3],cov[6];
    fpVtx->GetXYZ(pos);
    if ( fabs(pos[2])>10. ) return; // vertex cut on z-axis direction
    fpVtx->GetCovarianceMatrix(cov);
    
    pVertex.SetXYZ((Float_t)pos[0], (Float_t)pos[1], (Float_t)pos[2]);
    Float_t covF[6];
    for (Int_t i=0; i<6; i++) { covF[i] = (Float_t)cov[i]; }
    pVertex.SetCovarianceMatrix(covF);
    pVertex.SetChi2(fpVtx->GetChi2());
    pVertex.SetNDF(fpVtx->GetNDF());
    pVertex.SetNContributors(fpVtx->GetNContributors());

    KFParticle PV(pVertex);
    if(!fAnalCuts) return;
    
    FillEventROOTObjects();
    
    
    //------------------------------------------------
    // Main analysis done in this function
    //------------------------------------------------
   
    MakeAnaOmegacZeroFromCasc(aodEvent,mcArray, PV);
    
    PostData(2, fCounter);
    PostData(3, fOutputList);
    PostData(4, fTree_Event);  // stream the results the analysis of this event to the output manager which will take care of writing it to a file
    PostData(5, fTree_Omegac0); // used for Reco. level
    PostData(6, fTree_Omegac0MCGen); // used for Gen. level
    PostData(7, fTree_Omegac0_QA); // used for QA check

    return;
    
}
//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskSESemileptonicOmegac0KFP ::Terminate(Option_t *)
{
    
    // The Terminate() function is the last function to be called during
    // a query.
    
    AliAnalysisTaskSE::Terminate();
    
    fOutputList = dynamic_cast<TList*>(GetOutputData(3));
    if (!fOutputList){
        AliError("fOutputList not available");
        return;
    }
    
    // called at the END of the analysis (when all events are processed)

    return;
}
//---------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskSESemileptonicOmegac0KFP :: MakeMCAnalysis(TClonesArray *mcArray)
{
    // Analysis AliAODMCParticle
    
    Double_t Num_Omegac0 = 0;
    Int_t    nmcpart = mcArray -> GetEntriesFast();
    Int_t    NDaughters = 3;
    
    for (Int_t i=0; i<nmcpart; i++){
        AliAODMCParticle *mcpart = NULL;
        mcpart = (AliAODMCParticle*)mcArray->At(i);
        
        //     if (TMath::Abs(mcpart->GetPdgCode()) == 4332) cout << "GetNdaug  "<< mcpart->GetNDaughters() << "  and ndaught expected is " << NDaughters  << endl;  // check daughters
        
        if(TMath::Abs(mcpart->GetPdgCode())==4332 && mcpart->GetNDaughters()==NDaughters){  // 4332: Omegac0
            Int_t  CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcpart,kTRUE);
            Bool_t e_flag = kFALSE;
            Bool_t omega_flag = kFALSE;
            AliAODMCParticle *mcepart = NULL;
            AliAODMCParticle *mccascpart = NULL;
            for(Int_t idau = mcpart->GetDaughterFirst(); idau<mcpart->GetDaughterLast()+1; idau++){
                if(idau <0) break;
                AliAODMCParticle *mcdau = (AliAODMCParticle*)mcArray->At(idau);
                if(TMath::Abs(mcdau->GetPdgCode())==11){  // 11: electron
                    e_flag = kTRUE;
                    mcepart = mcdau;
                }
                if(TMath::Abs(mcdau->GetPdgCode())==3334){
                    omega_flag = kTRUE;
                    mccascpart = mcdau;
                }
            } // loop for daughters -> e and omega
        
          
            /*
            cout << "the cascade's daughters are:  ";
            AliAODMCParticle* mclam;
            for (Int_t iCascDau = mccascpart->GetDaughterFirst(); iCascDau<= mccascpart->GetDaughterLast(); iCascDau++ ) {
                AliAODMCParticle* dau = (AliAODMCParticle*)mcArray->At(iCascDau);
                //   cout << dau->GetPdgCode() << '\t';
                if (iCascDau == mccascpart->GetDaughterFirst() ) mclam = dau;
            }
            cout << endl << "and the first daughter's daughters are  ";
            for (Int_t ilamDau = mclam->GetDaughterFirst(); ilamDau<= mclam->GetDaughterLast(); ilamDau++ ) {
                AliAODMCParticle* dau = (AliAODMCParticle*)mcArray->At(ilamDau);
                //  cout << dau->GetPdgCode() << '\t';
            }
            cout << endl;
            */
              
            if(e_flag && omega_flag){
                
                AliAODMCParticle *mcdau_0 = (AliAODMCParticle*)mcArray->At(mcpart->GetDaughterFirst());
                Double_t MLOverP = sqrt( pow(mcpart->Xv() - mcdau_0->Xv(),2.) +  pow(mcpart->Yv() - mcdau_0->Yv(),2.) +  pow(mcpart->Zv() - mcdau_0->Zv(),2.)) * mcpart-> M() / mcpart->P()*1.e4;
                
                FillTreeGenOmegac0(mcpart,CheckOrigin, MLOverP);
            }
        } // Omegac0 4332
        
    }  // nmcpart
    
    return kTRUE;
    
//************************
}//MakeMCAnalysis
//---------------------------------------------------------------------------------------------
void AliAnalysisTaskSESemileptonicOmegac0KFP :: FillTreeGenOmegac0(AliAODMCParticle *mcpart, Int_t CheckOrigin, Double_t MLOverP)
{

    // Fill histograms or tree depending
    for (Int_t i=0; i<5; i++){
        fVar_Omegac0MCGen[i] = -9999.;
    }
    fVar_Omegac0MCGen[0] = mcpart->Y();
    fVar_Omegac0MCGen[1] = mcpart->Pt();
    fVar_Omegac0MCGen[2] = CheckOrigin;
    fVar_Omegac0MCGen[3] = mcpart->GetPdgCode();
    fVar_Omegac0MCGen[4] = MLOverP;
    
    if(fWriteOmegac0MCGenTree) fTree_Omegac0MCGen ->Fill();
    
    
} //FillTreeGenOmegac0
//--------------------------------------------------------------------------------------------
void AliAnalysisTaskSESemileptonicOmegac0KFP :: SelectTrack(const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks, Bool_t *seleFlags, TClonesArray *mcArray)
{
    // Select good tracks using fAnalCuts (AliRDHFCuts object) and return the array of their ids
    if(trkEntries==0) return;
    
    nSeleTrks =0;
    for(Int_t i=0;i<trkEntries;i++){
        seleFlags[i] = kFALSE;
        AliVTrack *track;
        track = (AliVTrack*)event -> GetTrack(i);
        
        Double_t convtest[21];
        if(!track->GetCovarianceXYZPxPyPz(convtest)) continue;
        if(!fAnalCuts) continue;
       
        AliAODTrack *aodt = (AliAODTrack*)track;
        if(fAnalCuts->GetProdUseAODFilterBit()){
            Int_t filterbit = fAnalCuts->GetProdAODFilterBit();
            if(filterbit==7){
                if(!aodt->TestFilterBit(BIT(filterbit))) continue;
            }else{
                if(!aodt->TestFilterMask(BIT(filterbit))) continue;
            }
        }
        
        AliAODTrack *aodtpid = 0;
        if(fAnalCuts->GetProdAODFilterBit()==4){
            aodtpid = aodt;
        }
        
        Double_t nsigma_tpcele = -9999.;
        Double_t nsigma_tofele = -9999.;
        Double_t nSigmaCombined_Ele = -9999.;
        
        if(fAnalCuts->GetIsUsePID()){
            nsigma_tpcele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(aodtpid,AliPID::kElectron);
            nsigma_tofele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(aodtpid,AliPID::kElectron);
        }
        if(fAnalCuts->SingleTrkCutsNoPID(aodt,aodtpid,fpVtx)){
            fHistoElectronTPCPID->Fill(aodt->Pt(),nsigma_tpcele);
            fHistoElectronTOFPID->Fill(aodt->Pt(),nsigma_tofele);
           
            nSigmaCombined_Ele = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigma_tpcele, nsigma_tofele);
           
            Double_t Ele_TPCPIDSelCombinedTPCTOF[3];
            Ele_TPCPIDSelCombinedTPCTOF[0] = nSigmaCombined_Ele;
            Ele_TPCPIDSelCombinedTPCTOF[1] = aodt->Pt();
            Ele_TPCPIDSelCombinedTPCTOF[2] = nsigma_tpcele;
            
            fHistoElectronTPCPIDSelTOF -> Fill(Ele_TPCPIDSelCombinedTPCTOF);
    
        }
        
        if(fAnalCuts->SingleTrkCuts(aodt,aodtpid,fpVtx)){
            seleFlags[i]=kTRUE;
            nSeleTrks++;
            
            nSigmaCombined_Ele = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigma_tpcele, nsigma_tofele);
            
            Double_t ElePID_TPC_TOF[4];
            ElePID_TPC_TOF[0] = nSigmaCombined_Ele;
            ElePID_TPC_TOF[1] = nsigma_tpcele;
            ElePID_TPC_TOF[2] = nsigma_tofele;
            ElePID_TPC_TOF[3] = aodt->Pt();
            
            fHistoElectronTPCTOFSelPID -> Fill(ElePID_TPC_TOF);
           
      }
    } // end loop on tracks
} // SelectTrack void()
//--------------------------------------------------------------------------------------------
void AliAnalysisTaskSESemileptonicOmegac0KFP :: SelectCascade(const AliVEvent *event,Int_t nCascs,Int_t &nSeleCasc, Bool_t *seleCascFlags, TClonesArray *mcArray)
{
    
    //----- Select good Casc Using fAnalCuts (AliRDHFCuts objects) and return the array of their ids
    
    Double_t primvert[3];
    fpVtx->GetXYZ(primvert);
    
    nSeleCasc = 0;
    for (Int_t iCasc = 0; iCasc<nCascs;iCasc++)
      {
          seleCascFlags[iCasc] = kFALSE;
          AliAODcascade *casc = ((AliAODEvent*)event)->GetCascade(iCasc);
          
          if(!fAnalCuts) continue;
          if(fAnalCuts->SingleCascadeCuts(casc,primvert)){
              seleCascFlags[iCasc] = kTRUE;
              nSeleCasc++;
          }
    }
    
}// SelectCascade()
//--------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskSESemileptonicOmegac0KFP ::PrefilterElectronULS(AliAODTrack *etrk, AliAODEvent *aodEvent, Double_t &minmass)
{
 
    Bool_t isconv = kFALSE;
    minmass = 9999.;
   
    Int_t trkid = etrk ->GetID();
    Double_t px1 = etrk->Px();
    Double_t py1 = etrk->Py();
    Double_t pz1 = etrk->Pz();
    Double_t E1 = sqrt(px1*px1+py1*py1+pz1*pz1+0.000511*0.000511);
    
    const Int_t nTrks = aodEvent -> GetNumberOfTracks();
    for(Int_t it=0;it < nTrks; it++){
        AliAODTrack *trk = (AliAODTrack*)aodEvent -> GetTrack(it);
        if(!trk) continue;
        if(!trk->TestFilterMask(BIT(fAnalCuts->GetProdAODFilterBit()))) continue;
        if(fabs(fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kElectron))>5) continue;
        Int_t trkid2 = trk ->GetID();
        if(fabs(trkid)==fabs(trkid2)) continue;
        if(trk->Charge()*etrk->Charge()>0) continue;
        
        Double_t px2 = trk->Px();
        Double_t py2 = trk->Py();
        Double_t pz2 = trk->Pz();
        Double_t E2 = sqrt(px2*px2+py2*py2+pz2*pz2+0.000511*0.000511);
        Double_t mass = sqrt(pow(E1+E2,2)-pow(px1+px2,2)-pow(py1+py2,2)-pow(pz1+pz2,2));
    
        if (mass<minmass)    minmass = mass;
     } // nTracks
    
    if(minmass<0.05) isconv = kTRUE;
    
    return isconv;
  
} //PrefilterElectronULS()
//--------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskSESemileptonicOmegac0KFP ::PrefilterElectronLS(AliAODTrack *etrk, AliAODEvent *aodEvent,  Double_t &minmass_ss)
{
 
    Bool_t isconv = kFALSE;
    minmass_ss = 9999.;
   
    Int_t trkid = etrk ->GetID();
    Double_t px1 = etrk->Px();
    Double_t py1 = etrk->Py();
    Double_t pz1 = etrk->Pz();
    Double_t E1 = sqrt(px1*px1+py1*py1+pz1*pz1+0.000511*0.000511);
    
    const Int_t nTrks = aodEvent -> GetNumberOfTracks();
    for(Int_t it=0;it < nTrks; it++){
        AliAODTrack *trk = (AliAODTrack*)aodEvent -> GetTrack(it);
        if(!trk) continue;
        if(!trk->TestFilterMask(BIT(fAnalCuts->GetProdAODFilterBit()))) continue;
        if(fabs(fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kElectron))>5) continue;
        Int_t trkid2 = trk ->GetID();
        if(fabs(trkid)==fabs(trkid2)) continue;
        if(trk->Charge()*etrk->Charge()<0) continue;
        
        Double_t px2 = trk->Px();
        Double_t py2 = trk->Py();
        Double_t pz2 = trk->Pz();
        Double_t E2 = sqrt(px2*px2+py2*py2+pz2*pz2+0.000511*0.000511);
        Double_t mass = sqrt(pow(E1+E2,2)-pow(px1+px2,2)-pow(py1+py2,2)-pow(pz1+pz2,2));
    
        if (mass<minmass_ss)    minmass_ss = mass;
     } // nTracks
    
    if(minmass_ss<0.05) isconv = kTRUE;
    
    return isconv;
   
} //PrefilterElectronLS()
//__________________________________________________________________________________________________________________

void AliAnalysisTaskSESemileptonicOmegac0KFP :: MakeAnaOmegacZeroFromCasc(AliAODEvent *aodEvent, TClonesArray *mcArray, KFParticle PV)
{
    // main analysis called from "UserExec"
    // set the magnetic field
    
    KFParticle :: SetField(fBzkG);
    Double_t covP[21], covN[21], covB[21];
    const Int_t NDaughters = 2;
    const Float_t massLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
    const Float_t massXi     = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
    const Float_t massK0S    = TDatabasePDG::Instance()->GetParticle(310)->Mass();
    const Float_t massOmega  = TDatabasePDG::Instance()->GetParticle(3334)->Mass();
    
    //========================== select good candidates for electron and cascade ==========================
    
    const Int_t nTracks = aodEvent -> GetNumberOfTracks();
    Bool_t seleTrkFlags[nTracks];
    Int_t  nSeleTrks=0;
    SelectTrack(aodEvent,nTracks,nSeleTrks,seleTrkFlags,mcArray);
    
    const Int_t nCascs = aodEvent -> GetNumberOfCascades();
    Bool_t  seleCascFlags[nCascs];
    Int_t   nSeleCasc=0;
    SelectCascade(aodEvent,nCascs,nSeleCasc,seleCascFlags,mcArray);
    
    Double_t mass; Double_t mass_ss; Double_t nSigmaCombined_Ele;
    for(Int_t itrk =0; itrk<nTracks; itrk++){
        if(!seleTrkFlags[itrk]) continue;
        AliAODTrack *trk = static_cast<AliAODTrack*>(aodEvent->GetTrack(itrk));
       
        Float_t nSigmaTOF_Ele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk, AliPID::kElectron);
        Float_t nSigmaTPC_Ele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk, AliPID::kElectron);  // need to know the data type conversion
        nSigmaCombined_Ele = AliVertexingHFUtils::CombineNsigmaTPCTOF(nSigmaTPC_Ele, nSigmaTOF_Ele);
        
        Double_t MassConversions[3];
        MassConversions[0] = nSigmaCombined_Ele;
        MassConversions[1] = mass;
        MassConversions[2] = mass_ss;
 
        fHistoMassConversions -> Fill(MassConversions);
       
        
    } // nTracks: e+/e-
        
    //--------------------- Omega candidates ------------------------------------
 
    for(Int_t iCasc =0; iCasc<nCascs; iCasc++){
        // cascade cut
        if(!seleCascFlags[iCasc]) continue;
        AliAODcascade *casc = aodEvent->GetCascade(iCasc);
        if(!casc) continue;
        
        AliAODTrack *ptrack = (AliAODTrack*)(casc->GetDaughter(0));
        AliAODTrack *ntrack = (AliAODTrack*)(casc->GetDaughter(1));
        AliAODTrack *btrack = (AliAODTrack*)(casc->GetDecayVertexXi()->GetDaughter(0));
        
        if( !ptrack || !ntrack || !btrack) continue;
        
        // define the charge of first daughter as positive, otherwise, define as the second one
        if(ptrack->Charge()<0){
            ptrack = (AliAODTrack*)(casc->GetDaughter(1));
            ntrack = (AliAODTrack*)(casc->GetDaughter(0));
        }
        
        if ( !ptrack->GetCovarianceXYZPxPyPz(covP) || !ntrack->GetCovarianceXYZPxPyPz(covN) || !btrack->GetCovarianceXYZPxPyPz(covB) ) continue;
       
        // check convariance matrix
        if( !AliVertexingHFUtils::CheckAODtrackCov(ptrack) || !AliVertexingHFUtils::CheckAODtrackCov(ntrack) || !AliVertexingHFUtils::CheckAODtrackCov(btrack)) continue;
   
        // using KFParticle to redefine daughters
        KFParticle kfpProton       = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ptrack,2212);
        KFParticle kfpPionMinus    = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ntrack,-211);
        KFParticle kfpAntiProton   = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ntrack,-2212);
        KFParticle kfpPionPlus     = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ptrack,211);
        
        KFParticle kfpElePlus      = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ptrack,-11);
        KFParticle kfpEleMinus     = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ntrack,11);
        
        //--- reconstrcut K0S
        KFParticle kfpK0Short;
        const KFParticle *vk0sDaughters[2] = {&kfpPionPlus, &kfpPionMinus};
        kfpK0Short.Construct(vk0sDaughters,NDaughters);
       
        //--- to reconstruct the V0 Lambda
       
        Int_t decaytype = -9999.;  // WS = 0; RS =1
        
        if( btrack->Charge() < 0){  // kaon-
            
            const KFParticle *vDaughters[2] = {&kfpProton, &kfpPionMinus};
            
            KFParticle kfpLambda;
            kfpLambda.Construct(vDaughters,NDaughters);
            Float_t massLambda_Rec, err_massLambda;
            kfpLambda.GetMass(massLambda_Rec, err_massLambda);
            
            //check rapidity of lambda
            if( TMath::Abs(kfpLambda.GetE()) <= TMath::Abs(kfpLambda.GetPz() )) continue;
               
            // chi2>0 && NDF>0 for selecting Lambda
            if ( (kfpLambda.GetNDF()<=0 || kfpLambda.GetChi2()<=0) ) continue;

            // check cov. of Lambda
            if( !AliVertexingHFUtils::CheckKFParticleCov(kfpLambda)) continue;
            
            // error_mass > 0 of Lambas
            if( err_massLambda <=0 ) continue;
              
            // Chi2geo cut of Lambda
            if( (kfpLambda.GetChi2()/kfpLambda.GetNDF()) >= fAnalCuts->GetKFPLam_Chi2geoMax()) continue; // defined
               
            //--------  calcualte l/Delatl for Lambda ------------------------------
               Double_t dx_Lambda = PV.GetX()-kfpLambda.GetX();
               Double_t dy_Lambda = PV.GetY()-kfpLambda.GetY();
               Double_t dz_Lambda = PV.GetZ()-kfpLambda.GetZ();
               Double_t l_Lambda = TMath::Sqrt(dx_Lambda*dx_Lambda + dy_Lambda*dy_Lambda + dz_Lambda*dz_Lambda);
               Double_t dl_Lambda = (PV.GetCovariance(0)+kfpLambda.GetCovariance(0))*dx_Lambda*dx_Lambda + (PV.GetCovariance(2)+kfpLambda.GetCovariance(2))*dy_Lambda*dy_Lambda + (PV.GetCovariance(5)+kfpLambda.GetCovariance(5))*dz_Lambda*dz_Lambda + 2*( (PV.GetCovariance(1)+kfpLambda.GetCovariance(1))*dx_Lambda*dy_Lambda + (PV.GetCovariance(3)+kfpLambda.GetCovariance(3))*dx_Lambda*dz_Lambda + (PV.GetCovariance(4)+kfpLambda.GetCovariance(4))*dy_Lambda*dz_Lambda );
               if ( fabs(l_Lambda)<1.e-8f ) l_Lambda = 1.e-8f;
               dl_Lambda = dl_Lambda<0. ? 1.e8f : sqrt(dl_Lambda)/l_Lambda;
               Double_t nErr_l_Lambda = l_Lambda/dl_Lambda;
            //---------------------------------------------------------------------
               
            // l/Delatl cut of Lambda
            if ( nErr_l_Lambda <= fAnalCuts->GetKFPLam_lDeltalMin() ) continue;
            // mass window cut of Lambda
            if(TMath::Abs(massLambda_Rec - massLambda) > (fAnalCuts->GetProdMassTolLambda()) ) continue;
        
            KFParticle kfpLambda_m = kfpLambda;
            kfpLambda_m.SetNonlinearMassConstraint(massLambda);
    
            if( !AliVertexingHFUtils::CheckKFParticleCov(kfpLambda_m) || TMath::Abs(kfpLambda.GetE()) <= TMath::Abs(kfpLambda.GetPz()) ) continue;
             
            KFParticle  kfpKaon = AliVertexingHFUtils::CreateKFParticleFromAODtrack(btrack, -321); // kaon-
      
            //--------------------- reconstruct OmegaMinus------------------------------------
    
            KFParticle kfpOmegaMinus;
            const KFParticle *vOmegaDs[2] = {&kfpKaon, &kfpLambda_m};
            kfpOmegaMinus.Construct(vOmegaDs,NDaughters);
        
            //check rapidity for Omega-
            if( TMath::Abs(kfpOmegaMinus.GetE()) <= TMath::Abs(kfpOmegaMinus.GetPz()) ) continue;
            //err_massOmega > 0
            Float_t massOmegaMinus_Rec, err_massOmegaMinus;
            kfpOmegaMinus.GetMass(massOmegaMinus_Rec,err_massOmegaMinus);
            if(err_massOmegaMinus <= 0 ) continue;
               
            //chi2 >0 && NDF>0
            if(kfpOmegaMinus.GetNDF() <=0 || kfpOmegaMinus.GetChi2() <=0) continue;
               
            //Chi2geo cut
            if(kfpOmegaMinus.GetChi2()/kfpOmegaMinus.GetNDF() >= fAnalCuts->GetKFPOmega_Chi2geoMax() ) continue;
            
            //check covariance matrix
            if( !AliVertexingHFUtils::CheckKFParticleCov(kfpOmegaMinus) )continue;
               
            //mass window cut of Omega-
            if(TMath::Abs(massOmegaMinus_Rec - massOmega) > (fAnalCuts->GetProdMassTolOmega() ) ) continue;
               
            KFParticle kfpOmegaMinus_m = kfpOmegaMinus;
            kfpOmegaMinus_m.SetNonlinearMassConstraint(massOmega); //  mass constraint on Omega
            
            if( !AliVertexingHFUtils::CheckKFParticleCov(kfpOmegaMinus_m) || TMath::Abs(kfpOmegaMinus_m.GetE()) <= TMath::Abs(kfpOmegaMinus_m.GetPz()) ) continue;
        
            
            //============================== reconstruct Omegac0 ( Omega+electron) ===============================
            
            for(Int_t itrkBE =0; itrkBE< nTracks; itrkBE++ ){  // loop for nTracks: e+/e-
                
                if(!seleTrkFlags[itrkBE]) continue;
                AliAODTrack *trkBE = (AliAODTrack*)aodEvent->GetTrack(itrkBE);
                
                Int_t pid = ptrack->GetID();
                Int_t nid = ntrack->GetID();
                Int_t bid = btrack->GetID();
                Int_t eid = trkBE->GetID();
                if ( (pid==eid)|| (nid==eid)|| (bid==eid) ) continue;
         
                KFParticle kfpBE;
                KFParticle kfpOmegac0;
                KFParticle kfpOmegac0_woMassConst;
               
                //=========== for electron charge > 0, make Right-Sign eOmega_KFP pairs =======================
                
                if(trkBE->Charge()>0 ) {
                   
                    decaytype = 1; // RS
                    kfpBE = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trkBE ,11);
           
                    //================= reconstrcut Omegac0_RS with Omega mass conostraint ===========
                    const KFParticle *vOmegacZeroDs[2] = {&kfpBE, &kfpOmegaMinus_m};
                    kfpOmegac0.Construct(vOmegacZeroDs,NDaughters);
             
                    // chi2>0 && NDF>0
                    if (kfpOmegac0.GetNDF()<= 0 || kfpOmegac0.GetChi2()<=0) continue;
                    if (kfpOmegac0.GetChi2()/kfpOmegac0.GetNDF() >= fAnalCuts->GetKFPOmegac0_Chi2geoMax()) continue;
                    //check rapidity of Omegac0
                    if (TMath::Abs(kfpOmegac0.GetE())<= TMath::Abs(kfpOmegac0.GetPz())) continue;
                    // err_massOmegac0 > 0
                    Float_t massOmegac0_Rec, err_massOmegac0;
                    kfpOmegac0.GetMass(massOmegac0_Rec, err_massOmegac0);
                    if(err_massOmegac0 <=0 ) continue;
                
                    //================ reconstrcut Omegac0_RS without Omega mass conostraint ============
                    const KFParticle *vOmegacZeroDs_woMassConst[2] = {&kfpBE, &kfpOmegaMinus};
                    kfpOmegac0_woMassConst.Construct(vOmegacZeroDs_woMassConst, NDaughters);
                
                    // chi2>0 && NDF>0
                    if (kfpOmegac0_woMassConst.GetNDF()<= 0 || kfpOmegac0_woMassConst.GetChi2()<=0) continue;
                    if (kfpOmegac0_woMassConst.GetChi2()/kfpOmegac0_woMassConst.GetNDF() >= fAnalCuts->GetKFPOmegac0_Chi2geoMax()) continue;
                    //check rapidity of Omegac0
                    if (TMath::Abs(kfpOmegac0_woMassConst.GetE())<= TMath::Abs(kfpOmegac0_woMassConst.GetPz())) continue;
                    // err_massOmegac0 > 0
                    Float_t massOmegac0_woMassConst_Rec, err_massOmegac0_woMassConst;
                    kfpOmegac0_woMassConst.GetMass(massOmegac0_woMassConst_Rec, err_massOmegac0_woMassConst);
                    if(err_massOmegac0_woMassConst <=0 ) continue;
                
                    if (fWriteOmegac0Tree){
                        Int_t lab_Omegac0 = -9999.;
                        if(fUseMCInfo){
                            lab_Omegac0 = MatchToMCOmegac0(ptrack,ntrack,btrack,trkBE,mcArray);
                         
                            if(lab_Omegac0>=0) FillTreeRecOmegac0FromCasc(kfpOmegac0,kfpOmegac0_woMassConst, trkBE, kfpBE, kfpOmegaMinus, kfpOmegaMinus_m, kfpKaon, btrack, casc, kfpK0Short,  kfpLambda, kfpLambda_m, ptrack, ntrack, PV, mcArray, aodEvent, lab_Omegac0, decaytype);
                      
                        } //fUseMCInfo
                        if(!fUseMCInfo){
                            FillTreeRecOmegac0FromCasc(kfpOmegac0,kfpOmegac0_woMassConst, trkBE, kfpBE, kfpOmegaMinus, kfpOmegaMinus_m, kfpKaon, btrack, casc, kfpK0Short,  kfpLambda, kfpLambda_m, ptrack, ntrack, PV, mcArray, aodEvent, lab_Omegac0, decaytype);
                        } // !fUseMCInfo
                    } //fWriteOmegac0Tree
                    kfpOmegac0.Clear();
                    kfpOmegac0_woMassConst.Clear();
                    kfpBE.Clear();
                } // for electron_charge > 0 loop
                
                
                //=========== for electron charge < 0, make Wrong-Sign eOmega_KFP pairs ================================
               
                if(trkBE->Charge()<0 ) {
                    
                    decaytype = 0; // WS
                    kfpBE = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trkBE, -11);
           
                    //================= reconstrcut Omegac0_WS with Omega mass conostraint ===========
                    const KFParticle *vOmegacZeroDs[2] = {&kfpBE, &kfpOmegaMinus_m};
                    kfpOmegac0.Construct(vOmegacZeroDs,NDaughters);
                      
                    // chi2>0 && NDF>0
                    if (kfpOmegac0.GetNDF()<= 0 || kfpOmegac0.GetChi2()<=0) continue;
                    if (kfpOmegac0.GetChi2()/kfpOmegac0.GetNDF() >= fAnalCuts->GetKFPOmegac0_Chi2geoMax()) continue;
                    //check rapidity of Omegac0
                    if (TMath::Abs(kfpOmegac0.GetE())<= TMath::Abs(kfpOmegac0.GetPz())) continue;
                    // err_massOmegac0 > 0
                    Float_t massOmegac0_Rec, err_massOmegac0;
                    kfpOmegac0.GetMass(massOmegac0_Rec, err_massOmegac0);
                    if(err_massOmegac0 <=0 ) continue;
                        
                    //================ reconstrcut Omegac0_WS without Omega mass conostraints ===========
                    const KFParticle *vOmegacZeroDs_woMassConst[2] = {&kfpBE, &kfpOmegaMinus};
                    kfpOmegac0_woMassConst.Construct(vOmegacZeroDs_woMassConst, NDaughters);
                         
                    // chi2>0 && NDF>0
                    if (kfpOmegac0_woMassConst.GetNDF()<= 0 || kfpOmegac0_woMassConst.GetChi2()<=0) continue;
                    if (kfpOmegac0_woMassConst.GetChi2()/kfpOmegac0_woMassConst.GetNDF() >= fAnalCuts->GetKFPOmegac0_Chi2geoMax()) continue;
                    //check rapidity of Omegac0
                    if (TMath::Abs(kfpOmegac0_woMassConst.GetE())<= TMath::Abs(kfpOmegac0_woMassConst.GetPz())) continue;
                    // err_massOmegac0 > 0
                    Float_t massOmegac0_woMassConst_Rec, err_massOmegac0_woMassConst;
                    kfpOmegac0_woMassConst.GetMass(massOmegac0_woMassConst_Rec, err_massOmegac0_woMassConst);
                    if(err_massOmegac0_woMassConst <=0 ) continue;
             
                    if (fWriteOmegac0Tree){
                        Int_t lab_Omegac0 = -9999.;
                        if(fUseMCInfo){
                            lab_Omegac0 = MatchToMCOmegac0(ptrack,ntrack,btrack,trkBE,mcArray);
                            
                            if(lab_Omegac0>=0)   FillTreeRecOmegac0FromCasc(kfpOmegac0,kfpOmegac0_woMassConst, trkBE, kfpBE, kfpOmegaMinus, kfpOmegaMinus_m,  kfpKaon, btrack, casc, kfpK0Short, kfpLambda, kfpLambda_m, ptrack, ntrack, PV, mcArray, aodEvent, lab_Omegac0, decaytype);
                        
                        } //fUseMCInfo
                        if(!fUseMCInfo){
                            FillTreeRecOmegac0FromCasc(kfpOmegac0,kfpOmegac0_woMassConst, trkBE, kfpBE, kfpOmegaMinus, kfpOmegaMinus_m, kfpKaon, btrack, casc, kfpK0Short, kfpLambda, kfpLambda_m, ptrack, ntrack, PV, mcArray, aodEvent, lab_Omegac0, decaytype);
                        } // !fUseMCInfo
                    } //fWriteOmegac0Tree
                    kfpOmegac0.Clear();
                    kfpOmegac0_woMassConst.Clear();
                    kfpBE.Clear();
                } // for electron_charge > 0 loop
            } // loop for nTracks: e+/e-
            kfpOmegaMinus_m.Clear();
            kfpOmegaMinus.Clear();
            kfpKaon.Clear();
            kfpLambda_m.Clear();
            kfpLambda.Clear();
    } // loop for btrack charge < 0 (kaon-)
        
    if ( btrack->Charge()>0 ) {  // kaon+
            
        const KFParticle *vAntiDaughters[2] = {&kfpPionPlus, &kfpAntiProton};
                
        KFParticle kfpAntiLambda;
        kfpAntiLambda.Construct(vAntiDaughters, NDaughters);
        Float_t massAntiLambda_Rec, err_massAntiLambda;
        kfpAntiLambda.GetMass(massAntiLambda_Rec, err_massAntiLambda);
                
        // check rapidity of Anti-Lambda
        if ( TMath::Abs(kfpAntiLambda.GetE())<=TMath::Abs(kfpAntiLambda.GetPz()) ) continue;

        // chi2>0 && NDF>0 for selecting Anti-Lambda
        if ( kfpAntiLambda.GetNDF()<=0 || kfpAntiLambda.GetChi2()<=0 ) continue;

        // check cov. of Anti-Lambda
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiLambda) ) continue;

        // err_mass>0 of Anti-Lambda
        if ( err_massAntiLambda<=0 ) continue;

        // Chi2geo cut of Anti-Lambda
        if ( (kfpAntiLambda.GetChi2()/kfpAntiLambda.GetNDF()) >= fAnalCuts->GetKFPLam_Chi2geoMax() ) continue;

        //---------------------------- calculate l/Δl for Anti-Lambda -----------------------------
        Double_t dx_AntiLambda = PV.GetX()-kfpAntiLambda.GetX();
        Double_t dy_AntiLambda = PV.GetY()-kfpAntiLambda.GetY();
        Double_t dz_AntiLambda = PV.GetZ()-kfpAntiLambda.GetZ();
        Double_t l_AntiLambda = TMath::Sqrt(dx_AntiLambda*dx_AntiLambda + dy_AntiLambda*dy_AntiLambda + dz_AntiLambda*dz_AntiLambda);
        Double_t dl_AntiLambda = (PV.GetCovariance(0)+kfpAntiLambda.GetCovariance(0))*dx_AntiLambda*dx_AntiLambda + (PV.GetCovariance(2)+kfpAntiLambda.GetCovariance(2))*dy_AntiLambda*dy_AntiLambda + (PV.GetCovariance(5)+kfpAntiLambda.GetCovariance(5))*dz_AntiLambda*dz_AntiLambda + 2*( (PV.GetCovariance(1)+kfpAntiLambda.GetCovariance(1))*dx_AntiLambda*dy_AntiLambda + (PV.GetCovariance(3)+kfpAntiLambda.GetCovariance(3))*dx_AntiLambda*dz_AntiLambda + (PV.GetCovariance(4)+kfpAntiLambda.GetCovariance(4))*dy_AntiLambda*dz_AntiLambda );
        if ( fabs(l_AntiLambda)<1.e-8f ) l_AntiLambda = 1.e-8f;
        dl_AntiLambda = dl_AntiLambda<0. ? 1.e8f : sqrt(dl_AntiLambda)/l_AntiLambda;
        Double_t nErr_l_AntiLambda = l_AntiLambda/dl_AntiLambda;
     //------------------------------------------------------------------------------------------
                
        // l/Deltal cut of Anti-Lambda
        if ( nErr_l_AntiLambda <= fAnalCuts->GetKFPLam_lDeltalMin() ) continue;

        // mass window cut of Anti-Lambda
        if ( TMath::Abs(massAntiLambda_Rec-massLambda) > (fAnalCuts->GetProdMassTolLambda()) ) continue;

        KFParticle kfpAntiLambda_m = kfpAntiLambda;
        kfpAntiLambda_m.SetNonlinearMassConstraint(massLambda);  // mass constraint on AntiLambda

        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiLambda_m) || TMath::Abs(kfpAntiLambda_m.GetE()) <= TMath::Abs(kfpAntiLambda_m.GetPz()) ) continue;
                
        KFParticle kfpKaon = AliVertexingHFUtils::CreateKFParticleFromAODtrack(btrack, 321); // kaon+
        
        //======================== reconstruct OmegaPlus ==============================
        
        KFParticle kfpOmegaPlus;
        const KFParticle *vOmegaDs[2] = {&kfpKaon, &kfpAntiLambda_m};
        kfpOmegaPlus.Construct(vOmegaDs, NDaughters);
                
        // check rapidity of Omega+
        if ( TMath::Abs(kfpOmegaPlus.GetE())<=TMath::Abs(kfpOmegaPlus.GetPz()) ) continue;

        // err_massXi > 0
        Float_t massOmegaPlus_Rec, err_massOmegaPlus;
        kfpOmegaPlus.GetMass(massOmegaPlus_Rec, err_massOmegaPlus);
        if ( err_massOmegaPlus<=0 ) continue;

        // chi2>0 && NDF>0
        if ( kfpOmegaPlus.GetNDF()<=0 || kfpOmegaPlus.GetChi2()<=0 ) continue;

        if ( kfpOmegaPlus.GetChi2()/kfpOmegaPlus.GetNDF() >= fAnalCuts->GetKFPOmega_Chi2geoMax() ) continue;

        // check covariance matrix
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpOmegaPlus) ) continue;
                
        // mass window cut of Omega
        if ( (TMath::Abs(massOmegaPlus_Rec-massOmega) > (fAnalCuts->GetProdMassTolOmega())) ) continue;
                
        // mass constraint on OmegaPlus
        KFParticle kfpOmegaPlus_m = kfpOmegaPlus;
        kfpOmegaPlus_m.SetNonlinearMassConstraint(massOmega);
                
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpOmegaPlus_m) || TMath::Abs(kfpOmegaPlus_m.GetE()) <= TMath::Abs(kfpOmegaPlus_m.GetPz()) ) continue;
        
        //================================ reconstruct AntiOmegac0  ==================================================
        
        for ( Int_t itrkBE = 0; itrkBE < nTracks; itrkBE++){  // loop for nTracks: e+/e-
          
            if(!seleTrkFlags[itrkBE]) continue;
            AliAODTrack *trkBE = (AliAODTrack*)aodEvent->GetTrack(itrkBE);
        
            Int_t pid = ptrack->GetID();
            Int_t nid = ntrack->GetID();
            Int_t bid = btrack->GetID();
            Int_t eid = trkBE->GetID();
            if ( (pid==eid)|| (nid==eid)|| (bid==eid) ) continue;
          
            KFParticle kfpBE;
            KFParticle kfpAntiOmegac0;
            KFParticle kfpAntiOmegac0_woMassConst;
            
            //=========== for electron charge < 0, to make Right-Sign eOmega_KFP pairs ==========================
           
            if (trkBE->Charge() <0 ){
                
                decaytype = 1; // RS
                kfpBE = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trkBE, -11);
            
                //============== reconstruct Anti-Omegac0_RS with Omega mass constraint =========
                const KFParticle *vOmegac0Ds[2] = {&kfpBE, &kfpOmegaPlus_m};
                kfpAntiOmegac0.Construct(vOmegac0Ds, NDaughters);
                // chi2>0 && NDF>0
                if ( kfpAntiOmegac0.GetNDF()<=0 || kfpAntiOmegac0.GetChi2()<=0 ) continue;
                // Prefilter
                if ( kfpAntiOmegac0.GetChi2()/kfpAntiOmegac0.GetNDF() >= fAnalCuts->GetKFPOmegac0_Chi2geoMax() ) continue;
                // check rapidity of Anti-Omegac0
                if ( TMath::Abs(kfpAntiOmegac0.GetE())<=TMath::Abs(kfpAntiOmegac0.GetPz()) ) continue;
                // check covariance matrix
                if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiOmegac0) ) continue;
                // err_massAntiXic0 > 0
                Float_t massAntiOmegac0_Rec, err_massAntiOmegac0;
                kfpAntiOmegac0.GetMass(massAntiOmegac0_Rec, err_massAntiOmegac0);
                if ( massAntiOmegac0_Rec<=0 ) continue;
              
                //============== reconstruct Anti-Omegac0_RS without Omega mass constraint =======
                const KFParticle *vOmegac0Ds_woMassConst[2] = {&kfpBE, &kfpOmegaPlus};
                kfpAntiOmegac0_woMassConst.Construct(vOmegac0Ds_woMassConst, NDaughters);
                // chi2>0 && NDF>0
                if ( kfpAntiOmegac0_woMassConst.GetNDF()<=0 || kfpAntiOmegac0_woMassConst.GetChi2()<=0 ) continue;
                // Prefilter
                if ( kfpAntiOmegac0_woMassConst.GetChi2()/kfpAntiOmegac0_woMassConst.GetNDF() >= fAnalCuts->GetKFPOmegac0_Chi2geoMax() ) continue;
                // check rapidity of Anti-Omegac0
                if ( TMath::Abs(kfpAntiOmegac0_woMassConst.GetE())<=TMath::Abs(kfpAntiOmegac0_woMassConst.GetPz()) ) continue;
                // check covariance matrix
                if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiOmegac0_woMassConst) ) continue;
                // err_massAntiXic0 > 0
                Float_t massAntiOmegac0_woMassConst_Rec, err_massAntiOmegac0_woMassConst;
                kfpAntiOmegac0_woMassConst.GetMass(massAntiOmegac0_woMassConst_Rec, err_massAntiOmegac0_woMassConst);
                if ( massAntiOmegac0_woMassConst_Rec<=0 ) continue;
            
                if (fWriteOmegac0Tree) {
                    Int_t lab_AntiOmegac0 = -9999.;
                    if (fUseMCInfo) {
                        lab_AntiOmegac0 = MatchToMCAntiOmegac0(ntrack, ptrack, btrack, trkBE, mcArray);
                
                    if (lab_AntiOmegac0>=0)
                        FillTreeRecOmegac0FromCasc(kfpAntiOmegac0, kfpAntiOmegac0_woMassConst, trkBE, kfpBE, kfpOmegaPlus, kfpOmegaPlus_m,  kfpKaon, btrack, casc, kfpK0Short, kfpAntiLambda, kfpAntiLambda_m, ntrack, ptrack, PV, mcArray,aodEvent, lab_AntiOmegac0, decaytype);
                    }
                    if (!fUseMCInfo) {
                        FillTreeRecOmegac0FromCasc(kfpAntiOmegac0, kfpAntiOmegac0_woMassConst, trkBE, kfpBE, kfpOmegaPlus, kfpOmegaPlus_m, kfpKaon, btrack, casc, kfpK0Short, kfpAntiLambda, kfpAntiLambda_m, ntrack, ptrack, PV, mcArray,aodEvent, lab_AntiOmegac0, decaytype);
                    }
                }  //fWriteOmegac0Tree
                kfpAntiOmegac0.Clear();
                kfpAntiOmegac0_woMassConst.Clear();
                kfpBE.Clear();
        } // electron charge < 0 loop
            
            //=========== for electron charge > 0, to make Wrong-Sign eOmega_KFP pairs ============================
         
            if (trkBE->Charge() >0 ){
                
                decaytype = 0; // WS
                kfpBE = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trkBE, 11);
            
                //================= reconstrcut AntiOmegac0_WS with Omega mass conostraint ===========
                const KFParticle *vOmegac0Ds[2] = {&kfpBE, &kfpOmegaPlus_m};
                kfpAntiOmegac0.Construct(vOmegac0Ds, NDaughters);
                // chi2>0 && NDF>0
                if ( kfpAntiOmegac0.GetNDF()<=0 || kfpAntiOmegac0.GetChi2()<=0 ) continue;
                // Prefilter
                if ( kfpAntiOmegac0.GetChi2()/kfpAntiOmegac0.GetNDF() >= fAnalCuts->GetKFPOmegac0_Chi2geoMax() ) continue;
                // check rapidity of Anti-Omegac0
                if ( TMath::Abs(kfpAntiOmegac0.GetE())<=TMath::Abs(kfpAntiOmegac0.GetPz()) ) continue;
                // check covariance matrix
                if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiOmegac0) ) continue;
                // err_massAntiXic0 > 0
                Float_t massAntiOmegac0_Rec, err_massAntiOmegac0;
                kfpAntiOmegac0.GetMass(massAntiOmegac0_Rec, err_massAntiOmegac0);
                if ( massAntiOmegac0_Rec<=0 ) continue;
            
                //============== reconstruct Anti-Omegac0_WS without Omega mass constraint ===========
                const KFParticle *vOmegac0Ds_woMassConst[2] = {&kfpBE, &kfpOmegaPlus};
                kfpAntiOmegac0_woMassConst.Construct(vOmegac0Ds_woMassConst, NDaughters);
                // chi2>0 && NDF>0
                if ( kfpAntiOmegac0_woMassConst.GetNDF()<=0 || kfpAntiOmegac0_woMassConst.GetChi2()<=0 ) continue;
                // Prefilter
                if ( kfpAntiOmegac0_woMassConst.GetChi2()/kfpAntiOmegac0_woMassConst.GetNDF() >= fAnalCuts->GetKFPOmegac0_Chi2geoMax() ) continue;
                // check rapidity of Anti-Omegac0
                if ( TMath::Abs(kfpAntiOmegac0_woMassConst.GetE())<=TMath::Abs(kfpAntiOmegac0_woMassConst.GetPz()) ) continue;
                // check covariance matrix
                if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiOmegac0_woMassConst) ) continue;
                // err_massAntiXic0 > 0
                Float_t massAntiOmegac0_woMassConst_Rec, err_massAntiOmegac0_woMassConst;
                kfpAntiOmegac0_woMassConst.GetMass(massAntiOmegac0_woMassConst_Rec, err_massAntiOmegac0_woMassConst);
                if ( massAntiOmegac0_woMassConst_Rec<=0 ) continue;
         
                if (fWriteOmegac0Tree) {
                    Int_t lab_AntiOmegac0 = -9999.;
                    if (fUseMCInfo) {
                        lab_AntiOmegac0 = MatchToMCAntiOmegac0(ntrack, ptrack, btrack, trkBE, mcArray);
                        
                    if (lab_AntiOmegac0>=0)
                        FillTreeRecOmegac0FromCasc(kfpAntiOmegac0, kfpAntiOmegac0_woMassConst, trkBE, kfpBE, kfpOmegaPlus, kfpOmegaPlus_m,  kfpKaon, btrack, casc, kfpK0Short, kfpAntiLambda, kfpAntiLambda_m, ntrack, ptrack, PV, mcArray,aodEvent, lab_AntiOmegac0, decaytype);
                    }
                    if (!fUseMCInfo) {
                        FillTreeRecOmegac0FromCasc(kfpAntiOmegac0, kfpAntiOmegac0_woMassConst, trkBE, kfpBE, kfpOmegaPlus, kfpOmegaPlus_m, kfpKaon, btrack, casc, kfpK0Short, kfpAntiLambda, kfpAntiLambda_m, ntrack, ptrack, PV, mcArray,aodEvent, lab_AntiOmegac0, decaytype);
                    }
                }  //fWriteOmegac0Tree
                kfpAntiOmegac0.Clear();
                kfpAntiOmegac0_woMassConst.Clear();
                kfpBE.Clear();
            } // electron charge > 0: KFP_eOmega_WS pairs
        } // loop for nTracks: e+/e-
        kfpOmegaPlus_m.Clear();
        kfpOmegaPlus.Clear();
        kfpKaon.Clear();
        kfpAntiLambda_m.Clear();
        kfpAntiLambda.Clear();
        } // loop for btrack charge > 0 (kaon+)
        kfpK0Short.Clear();
        kfpPionPlus.Clear();
        kfpAntiProton.Clear();
        kfpPionMinus.Clear();
        kfpProton.Clear();
    }// loop for cascade
        
    return;

}  //MakeAnaOmegacZeroFromCasc
//-----------------------------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskSESemileptonicOmegac0KFP :: MatchToMCOmegac0(AliAODTrack* trackProton, AliAODTrack *trackPionMinus, AliAODTrack *trackKaon, AliAODTrack *trackEletron, TClonesArray *mcArray )
{
    // check if all of the tracks is matched to a MC signal
    // if no, return -1;
    // if yes, return label(>=0) of the AliAODMCParticle
    
    Int_t labelProton = fabs(trackProton->GetLabel());
    if (labelProton<0) return -1;
    AliAODMCParticle *mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));

    Int_t labelPionMinus = fabs(trackPionMinus->GetLabel());
    if(labelPionMinus<0) return -1;
    AliAODMCParticle *mcPionMinus = static_cast<AliAODMCParticle*>(mcArray->At(labelPionMinus));
 
    Int_t labelKaon = fabs(trackKaon->GetLabel());
    if(labelKaon<0) return -1;
    AliAODMCParticle *mcKaon = static_cast<AliAODMCParticle*>(mcArray->At(labelKaon));

    Int_t labelElectron = fabs(trackEletron->GetLabel());
    if(labelElectron<0) return -1;
    AliAODMCParticle *mcElectron = static_cast<AliAODMCParticle*>(mcArray->At(labelElectron));

    if(mcProton->GetPdgCode() != 2122 || mcPionMinus->GetPdgCode() != -211 || mcKaon->GetPdgCode() != -321 ||  TMath::Abs(mcElectron->GetPdgCode()) != 11 ) return -1;
  
    //--------- check V0
    Int_t IndexMother[2];
    IndexMother[0] = mcProton -> GetMother();
    IndexMother[1] = mcPionMinus -> GetMother();
    if(IndexMother[0]<0 || IndexMother[1]<0) return -1; // check their mother -- Lambda exist
    if(IndexMother[0] != IndexMother[1]) return -1;// check they are from same mother
    
    AliAODMCParticle *mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
    if (mcMother-> GetPdgCode() != 3122 || mcMother->GetNDaughters()!=2 ) return -1; //check mother is Lambda

    //--------- check Cascade
    IndexMother[0] = mcMother->GetMother();
    IndexMother[1] = mcKaon->GetMother();
    if(IndexMother[0]<0 || IndexMother[1]<0) return -1; // check their mother -- Omega- exist
    if(IndexMother[0] != IndexMother[1]) return -1;// check they are from same mother
    
    mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
    if( mcMother -> GetPdgCode() != 3334 || mcMother ->GetNDaughters() != 2) return -1; // check the mother of Lambda and Kaon is Omega-
    
    //---------- check Omegac0
    IndexMother[0] = mcMother->GetMother();
    IndexMother[1] = mcElectron->GetMother();
    if(IndexMother[0]<0 || IndexMother[1]<0) return -1; // check their mother -- Omegac0 exist
    if(IndexMother[0] != IndexMother[1]) return -1;// check they are from same mother
    
    mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
    if ( mcMother->GetPdgCode() != 4332 || mcMother->GetNDaughters()!=3 ) return -1; // check mother is Omegac0
   
    Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcMother,kTRUE);
    return CheckOrigin;
    
}//MatchToMCOmegac0()
//-------------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskSESemileptonicOmegac0KFP :: MatchToMCAntiOmegac0(AliAODTrack* trackAntiProton, AliAODTrack *trackPionPlus, AliAODTrack *trackKaon, AliAODTrack *trackEletron, TClonesArray *mcArray )
{
    // check if all of the tracks is matched to a MC signal
    // if no, return -1;
    // iif yes, return label(>=0) of the AliAODMCParticle
    
    
    Int_t labelAntiProton = fabs(trackAntiProton->GetLabel());
    if (labelAntiProton<0) return -1;
    AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiProton));
    
    Int_t labelPionPlus  = fabs(trackPionPlus->GetLabel());
    if (labelPionPlus<0) return -1;
    AliAODMCParticle* mcPionPlus = static_cast<AliAODMCParticle*>(mcArray->At(labelPionPlus));
    
    Int_t labelKaon  = fabs(trackKaon->GetLabel());
    if (labelKaon<0) return -1;
    AliAODMCParticle* mcKaon = static_cast<AliAODMCParticle*>(mcArray->At(labelKaon));
   
    Int_t labelElectron  = fabs(trackEletron->GetLabel());
    if (labelElectron<0) return -1;
    AliAODMCParticle* mcElectron = static_cast<AliAODMCParticle*>(mcArray->At(labelElectron));
 
    if ( mcAntiProton->GetPdgCode() != -2212 || mcPionPlus->GetPdgCode() != 211 || mcKaon->GetPdgCode() != 321 || (TMath::Abs(mcElectron->GetPdgCode()) ) != 11) return -1; // check pdg
    
    //--------- check V0
    Int_t IndexMother[2];
    IndexMother[0] = mcAntiProton->GetMother();
    IndexMother[1] = mcPionPlus->GetMother();
    if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
    if ( IndexMother[0] != IndexMother[1] ) return -1; // check they are from same mother
    AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
    if ( mcMother->GetPdgCode() != -3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is anti-lambda

    //--------- check Cascade
    IndexMother[0] = mcMother->GetMother(); // mother of lambda
    IndexMother[1] = mcKaon->GetMother();
    if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
    if ( IndexMother[0] != IndexMother[1] ) return -1; // check they are from same mother
    mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
    if ( mcMother->GetPdgCode() != -3334 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Omega+
 
    //---------- check Omegac0
    IndexMother[0] = mcMother->GetMother(); // mother of Omega+
    IndexMother[1] = mcElectron->GetMother();
    if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
    if ( IndexMother[0] != IndexMother[1] ) return -1; // check thet are from same mother
    mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
    if ( mcMother->GetPdgCode() != -4332 || mcMother->GetNDaughters()!=3 ) return -1; // check mother is anti-Omegac0

    Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcMother,kTRUE);
    return CheckOrigin;
        
}//MatchToMCANtiOmegac0()
//----------------------------------------------------------------------------------------
void AliAnalysisTaskSESemileptonicOmegac0KFP :: DefineAnaHist()
{
    
    // example of histograms
    fHistEvents = new TH1F("fHistEvents", "fHistEvents", 18, 0.5, 18.5);
    fHistEvents->GetXaxis()->SetBinLabel(1,"Analyzed events");
    fHistEvents->GetXaxis()->SetBinLabel(2,"AliAODVertex exists");
    fHistEvents->GetXaxis()->SetBinLabel(3,"TriggerOK");
    fHistEvents->GetXaxis()->SetBinLabel(4,"IsEventSelected");
    fHistEvents->GetXaxis()->SetBinLabel(5,"V0 exists");
    fHistEvents->GetXaxis()->SetBinLabel(6,"Cascade exists");
    fHistEvents->GetXaxis()->SetBinLabel(7,"MCarray exists");
    fHistEvents->GetXaxis()->SetBinLabel(8,"MCheader exists");
    fHistEvents->GetXaxis()->SetBinLabel(9,"triggerClass!=CINT1");
    fHistEvents->GetXaxis()->SetBinLabel(10,"triggerMask!=kAnyINT");
    fHistEvents->GetXaxis()->SetBinLabel(11,"triggerMask!=kAny");
    fHistEvents->GetXaxis()->SetBinLabel(12,"vtxTitle.Contains(Z)");
    fHistEvents->GetXaxis()->SetBinLabel(13,"vtxTitle.Contains(3D)");
    fHistEvents->GetXaxis()->SetBinLabel(14,"vtxTitle.Doesn'tContain(Z-3D)");
    fHistEvents->GetXaxis()->SetBinLabel(15,Form("zVtx<=%2.0fcm", fAnalCuts->GetMaxVtxZ()));
    fHistEvents->GetXaxis()->SetBinLabel(16,"!IsEventSelected");
    fHistEvents->GetXaxis()->SetBinLabel(17,"triggerMask!=kAnyINT || triggerClass!=CINT1");
    fHistEvents->GetXaxis()->SetBinLabel(18,Form("zVtxMC<=%2.0fcm",fAnalCuts->GetMaxVtxZ()));
    fHistEvents->GetYaxis()->SetTitle("counts");
    
    fHTrigger = new TH1F("fHTrigger", "counter", 18, -0.5, 17.5);
    fHTrigger->SetStats(kTRUE);
    fHTrigger->GetXaxis()->SetBinLabel(1,"X1");
    fHTrigger->GetXaxis()->SetBinLabel(2,"kMB");
    fHTrigger->GetXaxis()->SetBinLabel(3,"kSemiCentral");
    fHTrigger->GetXaxis()->SetBinLabel(4,"kCentral");
    fHTrigger->GetXaxis()->SetBinLabel(5,"kINT7");
    fHTrigger->GetXaxis()->SetBinLabel(6,"kEMC7");
    //fHTrigger->GetXaxis()->SetBinLabel(7,"Space");
    fHTrigger->GetXaxis()->SetBinLabel(8,"kMB|kSemiCentral|kCentral");
    fHTrigger->GetXaxis()->SetBinLabel(9,"kINT7|kEMC7");
    fHTrigger->GetXaxis()->SetBinLabel(11,"kMB&kSemiCentral");
    fHTrigger->GetXaxis()->SetBinLabel(12,"kMB&kCentral");
    fHTrigger->GetXaxis()->SetBinLabel(13,"kINT7&kEMC7");
    
    fOutputList->Add(fHistEvents);
    fOutputList->Add(fHTrigger);
    

  // Define analysis histograms
    fHistoElectronTPCPID = new TH2F("fHistoElectronTPCPID","",50,0.,5.,50,-20.,20.);
    fOutputList -> Add(fHistoElectronTPCPID);
    fHistoElectronTOFPID = new TH2F("fHistoElectronTOFPID","",50,0.,5.,50,-20.,20.);
    fOutputList -> Add(fHistoElectronTOFPID);
    
    Int_t bins_ele[3]    = {1000, 10, 500};
    Double_t xmin_ele[3] = {0., 0., -10.};
    Double_t xmax_ele[3] = {10., 5., 10.};
    fHistoElectronTPCPIDSelTOF = new THnSparseF("fHistoElectronTPCPIDSelTOF","",3,bins_ele,xmin_ele,xmax_ele);
    fOutputList -> Add(fHistoElectronTPCPIDSelTOF);
    
    Int_t bins_TPCTOF[4]=    {1000,500,500,50};
    Double_t xmin_TPCTOF[4]={0.,-10.,-10.,0.};
    Double_t xmax_TPCTOF[4]={10.,10.,10.,5};
    fHistoElectronTPCTOFSelPID = new THnSparseF("fHistoElectronTPCTOFSelPID","",4,bins_TPCTOF,xmin_TPCTOF,xmax_TPCTOF);
    fOutputList -> Add(fHistoElectronTPCTOFSelPID);
    
    Int_t bins[3] = {1000, 500, 500};
    Double_t xmin[3] = {0.,0.,0.};
    Double_t xmax[3] = {10.,0.5,0.5};
    fHistoMassConversions = new THnSparseF("fHistoMassConversions","",3,bins,xmin,xmax);
    fOutputList->Add(fHistoMassConversions);
    
    
    return;
  
}// DefineAnaHist
//--------------------------------------------------------------------------------------------
void AliAnalysisTaskSESemileptonicOmegac0KFP :: DefineTreeMCGenOmegac0()
{

    const char* nameoutput = GetOutputSlot(6)->GetContainer()->GetName();
    fTree_Omegac0MCGen = new TTree(nameoutput,"Omegac0 MC varibales tree");
    Int_t nVar = 5;
    fVar_Omegac0MCGen = new Float_t[nVar];
    TString *fVarNames = new TString[nVar];
    
    fVarNames[0] = "rap_Omegac0";
    fVarNames[1] = "pT_Omegac0";
    fVarNames[2] = "CheckOrigin_SourceOmegac0";
    fVarNames[3] = "PDG_Omegac0";
    fVarNames[4] = "ct_Omegac0";
    
    for(Int_t ivar=0; ivar<nVar; ivar++){
        fTree_Omegac0MCGen->Branch(fVarNames[ivar].Data(),&fVar_Omegac0MCGen[ivar], Form("%s/F",fVarNames[ivar].Data()));
    }
    
    return;
    
}//DefineTreeMCGenOmegac0
//------------------------------------------------------------------------------------------------
void AliAnalysisTaskSESemileptonicOmegac0KFP::FillEventROOTObjects()
{
    for (Int_t i=0; i<7; i++) {
      fVar_Event[i] = 0.;
    }

    Double_t pos[3];
    fpVtx->GetXYZ(pos);

    fVar_Event[1] = pos[2];

    fTree_Event->Fill();

    return;

}  // FillEventROOTObjects
//------------------------------------------------------------------------------------------------
void AliAnalysisTaskSESemileptonicOmegac0KFP ::FillTreeRecOmegac0FromCasc(KFParticle kfpOmegac0, KFParticle kfpOmegac0_woMassConst, AliAODTrack *trackElectronFromOmegac0, KFParticle kfpBE, KFParticle kfpOmegaMinus, KFParticle kfpOmegaMinus_m,  KFParticle kfpKaon, AliAODTrack *trackKaonFromOmega, AliAODcascade *casc, KFParticle kfpK0Short,  KFParticle kfpLambda, KFParticle kfpLambda_m, AliAODTrack *trkProton, AliAODTrack *trkPion, KFParticle PV, TClonesArray *mcArray, AliAODEvent *aodEvent, Int_t lab_Omegac0, Int_t decaytype)
{
    
    for (Int_t i=0; i< 35 ; i++){
        fVar_Omegac0[i] = -9999.;
    }
    
    for (Int_t i =0; i< 33; i++){
        fVar_Omegac0_QA[i] = -9999.;   // QA check
    }

    //------------ ee pairs Conversion type ---------
    Double_t mass; Double_t mass_ss;
    Bool_t Convee_LS =  (PrefilterElectronLS(trackElectronFromOmegac0, aodEvent, mass_ss));
    Bool_t Convee_ULS = (PrefilterElectronULS(trackElectronFromOmegac0, aodEvent, mass));
    
    //------ Calculate EleOmegaOA ---------------
    Double_t pxe = trackElectronFromOmegac0->Px();
    Double_t pye = trackElectronFromOmegac0->Py();
    Double_t pze = trackElectronFromOmegac0->Pz();
    Double_t mome = sqrt(pxe*pxe+pye*pye+pze*pze);
    Double_t Ee = sqrt(mome*mome+0.000511*0.000511);

    Double_t pxOmega = casc->MomXiX();
    Double_t pyOmega = casc->MomXiY();
    Double_t pzOmega = casc->MomXiZ();
    
    Double_t momOmega = sqrt(pxOmega*pxOmega+pyOmega*pyOmega+pzOmega*pzOmega);
    Double_t EOmega = sqrt(momOmega*momOmega+1.67245*1.67245);
    
    Double_t cosoa = (pxe*pxOmega+pye*pyOmega+pze*pzOmega)/mome/momOmega; // CosOpeningAngle

    
    Float_t nSigmaTOF_EleFromOmegac0 = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trackElectronFromOmegac0, AliPID::kElectron);
    Float_t nSigmaTPC_EleFromOmegac0 = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trackElectronFromOmegac0, AliPID::kElectron);
    Float_t nSigmaTPC_PiFromLam = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkPion, AliPID::kPion);
    Float_t nSigmaTPC_PrFromLam = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkProton,AliPID::kProton);
    Float_t nSigmaTPC_KaonFromOmega = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trackKaonFromOmega,AliPID::kKaon);
    Float_t nSigmaTOF_KaonFromOmega = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trackKaonFromOmega,AliPID::kKaon);
    
    if ( fabs(nSigmaTPC_PiFromLam)>=4. || fabs(nSigmaTPC_PrFromLam)>=4. || fabs(nSigmaTPC_KaonFromOmega)>=4. ) return;

    KFParticle kfpOmegac0_PV = kfpOmegac0;
    kfpOmegac0_PV.SetProductionVertex(PV);
    
    KFParticle kfpOmegac0_woMassConst_PV = kfpOmegac0_woMassConst;
    kfpOmegac0_woMassConst_PV.SetProductionVertex(PV);
    
    KFParticle kfpOmegaMinus_Omegac0 = kfpOmegaMinus_m;
    kfpOmegaMinus_Omegac0.SetProductionVertex(kfpOmegac0);
    KFParticle kfpOmegaMinus_PV = kfpOmegaMinus_m;
    kfpOmegaMinus_PV.SetProductionVertex(PV);
    
    KFParticle kfpLambda_Omega = kfpLambda_m;
    kfpLambda_Omega.SetProductionVertex(kfpOmegaMinus);
    KFParticle kfpLambda_PV = kfpLambda_m;
    kfpLambda_PV.SetProductionVertex(PV);
    
    KFParticle kfpBE_Omegac0 = kfpBE;
    kfpBE_Omegac0.SetProductionVertex(kfpOmegac0);
    
    //----------------- calculate l/Δl for Lambda -----------------
    Double_t dx_Lambda = PV.GetX()-kfpLambda.GetX();
    Double_t dy_Lambda = PV.GetY()-kfpLambda.GetY();
    Double_t dz_Lambda = PV.GetZ()-kfpLambda.GetZ();
    Double_t l_Lambda = TMath::Sqrt(dx_Lambda*dx_Lambda + dy_Lambda*dy_Lambda + dz_Lambda*dz_Lambda);
    Double_t dl_Lambda = (PV.GetCovariance(0)+kfpLambda.GetCovariance(0))*dx_Lambda*dx_Lambda + (PV.GetCovariance(2)+kfpLambda.GetCovariance(2))*dy_Lambda*dy_Lambda + (PV.GetCovariance(5)+kfpLambda.GetCovariance(5))*dz_Lambda*dz_Lambda + 2*( (PV.GetCovariance(1)+kfpLambda.GetCovariance(1))*dx_Lambda*dy_Lambda + (PV.GetCovariance(3)+kfpLambda.GetCovariance(3))*dx_Lambda*dz_Lambda + (PV.GetCovariance(4)+kfpLambda.GetCovariance(4))*dy_Lambda*dz_Lambda );
    if ( fabs(l_Lambda)<1.e-8f ) l_Lambda = 1.e-8f;
    dl_Lambda = dl_Lambda<0. ? 1.e8f : sqrt(dl_Lambda)/l_Lambda;
    if ( dl_Lambda<=0 ) return;
 
    //----------------- calculate l/Δl for Omega- --------------------
    Double_t dx_Omega = PV.GetX()-kfpOmegaMinus.GetX();
    Double_t dy_Omega = PV.GetY()-kfpOmegaMinus.GetY();
    Double_t dz_Omega = PV.GetZ()-kfpOmegaMinus.GetZ();
    Double_t l_Omega  = TMath::Sqrt(dx_Omega*dx_Omega + dy_Omega*dy_Omega + dz_Omega*dz_Omega);
    Double_t dl_Omega = (PV.GetCovariance(0)+kfpOmegaMinus.GetCovariance(0))*dx_Omega*dx_Omega + (PV.GetCovariance(2)+kfpOmegaMinus.GetCovariance(2))*dy_Omega*dy_Omega + (PV.GetCovariance(5)+kfpOmegaMinus.GetCovariance(5))*dz_Omega*dz_Omega + 2*( (PV.GetCovariance(1)+kfpOmegaMinus.GetCovariance(1))*dx_Omega*dy_Omega + (PV.GetCovariance(3)+kfpOmegaMinus.GetCovariance(3))*dx_Omega*dz_Omega + (PV.GetCovariance(4)+kfpOmegaMinus.GetCovariance(4))*dy_Omega*dz_Omega );
    if ( fabs(l_Omega)<1.e-8f ) l_Omega = 1.e-8f;
    dl_Omega = dl_Omega<0. ? 1.e8f : sqrt(dl_Omega)/l_Omega;
    if ( dl_Omega<=0 ) return;

    //--------------------------------------------------------------------
    if ( kfpLambda_PV.GetChi2()/kfpLambda_PV.GetNDF() <= fAnalCuts->GetKFPLam_Chi2topoMin() ) return;
    if ( kfpOmegaMinus_PV.GetChi2()/kfpOmegaMinus_PV.GetNDF() >= fAnalCuts->GetKFPOmega_Chi2topoMax() ) return;
    if ( l_Omega/dl_Omega <= fAnalCuts->GetKFPOmega_lDeltalMin() ) return;
    
   
    fVar_Omegac0[0] = AliVertexingHFUtils::CombineNsigmaTPCTOF(nSigmaTPC_EleFromOmegac0, nSigmaTOF_EleFromOmegac0);  // for electron
    fVar_Omegac0[1] = AliVertexingHFUtils::CombineNsigmaTPCTOF(nSigmaTPC_KaonFromOmega,nSigmaTOF_KaonFromOmega); // for kaon
    fVar_Omegac0[2]  = nSigmaTPC_PiFromLam;  // pion
    fVar_Omegac0[3]  = nSigmaTPC_PrFromLam; // proton
    fVar_Omegac0[4]  = casc -> DcaXiDaughters();
    fVar_Omegac0[5]  = kfpLambda.GetChi2()/kfpLambda.GetNDF(); // chi2geo_Lam
    fVar_Omegac0[6]  = l_Lambda/dl_Lambda;
    fVar_Omegac0[7] = kfpOmegaMinus.GetChi2()/kfpOmegaMinus.GetNDF();
    fVar_Omegac0[8] = l_Omega/dl_Omega;
    fVar_Omegac0[9] = kfpOmegaMinus_PV.GetChi2()/kfpOmegaMinus_PV.GetNDF();
    
    Float_t DecayLxy_Lam, err_DecayLxy_Lam;
    kfpLambda_Omega.GetDecayLengthXY(DecayLxy_Lam, err_DecayLxy_Lam);
    fVar_Omegac0[10] = DecayLxy_Lam;
    
    Float_t DecayLxy_Omega, err_DecayLxy_Omega;
    kfpOmegaMinus_PV.GetDecayLengthXY(DecayLxy_Omega, err_DecayLxy_Omega);
    fVar_Omegac0[11] = DecayLxy_Omega;
    
    // calculate CosPointingAngle
    Double_t cosPA_v0toOmega = AliVertexingHFUtils::CosPointingAngleFromKF(kfpLambda_m, kfpOmegaMinus);
    Double_t cosPA_OmegaToPV = AliVertexingHFUtils::CosPointingAngleFromKF(kfpOmegaMinus_m, PV);
    fVar_Omegac0[12] = TMath::ACos(cosPA_v0toOmega); // PA_LamToOmega
    fVar_Omegac0[13] = TMath::ACos(cosPA_OmegaToPV); // PA_OmegaToPV
    
    Float_t mass_Lam, err_mass_Lam;
    kfpLambda.GetMass(mass_Lam, err_mass_Lam);
    fVar_Omegac0[14] = mass_Lam;

    Float_t mass_Omega, err_mass_Omega;
    kfpOmegaMinus.GetMass(mass_Omega, err_mass_Omega);
    fVar_Omegac0[15] = mass_Omega;

    fVar_Omegac0[16] = trackElectronFromOmegac0->Pt();
    
    fVar_Omegac0[17] = kfpOmegac0_PV.GetPt();
    if ( TMath::Abs(kfpOmegac0_PV.GetE())>TMath::Abs(kfpOmegac0_PV.GetPz()) ) {
      fVar_Omegac0[18] = kfpOmegac0_PV.GetRapidity();
    }
    
    const Float_t PDGmassOmegac0 = 2.6952;
    Float_t mass_Omegac0_PV, err_mass_Omegac0_PV;
    kfpOmegac0_PV.GetMass(mass_Omegac0_PV, err_mass_Omegac0_PV);
    fVar_Omegac0[19] = mass_Omegac0_PV;
 
    KFParticle kfpBE_PV = kfpBE;
    kfpBE_PV.SetProductionVertex(PV);
    fVar_Omegac0[20] = kfpBE_PV.GetChi2()/kfpBE_PV.GetNDF();  //chi2_prim of Electron to PV
    
    fVar_Omegac0[21] = kfpBE.GetDistanceFromVertexXY(PV); // DCA of electron in x-y to PV
    
    Float_t massK0S_Rec, err_massK0S;
    kfpK0Short.GetMass(massK0S_Rec, err_massK0S);
    fVar_Omegac0[22] = massK0S_Rec;
    
    fVar_Omegac0[23] = kfpLambda_Omega.GetChi2()/kfpLambda_Omega.GetNDF();
    fVar_Omegac0[24] = kfpOmegac0.GetChi2()/kfpOmegac0.GetNDF();
    fVar_Omegac0[25] = casc->DcaV0Daughters(); // DCA_LamDau
    fVar_Omegac0[26] = kfpKaon.GetDistanceFromParticle(kfpLambda_m); // DCA_OmegaDau_KF
    fVar_Omegac0[27] = kfpBE.GetDistanceFromParticle(kfpOmegaMinus_m); // DCA_Omegac0Dau_KF
    fVar_Omegac0[28] = kfpOmegaMinus_m.GetDistanceFromVertexXY(PV); // DCA of Omega in x-y plan to PV
    
    
    AliAODMCParticle *mcOmegac0 = 0;
    if (fUseMCInfo) {
        fVar_Omegac0[29] = lab_Omegac0;
        if (lab_Omegac0>=0) {
          Int_t labelEleFromOmegac0 = fabs(trackElectronFromOmegac0->GetLabel());

          AliAODMCParticle *mcEleFromOmegac0 = static_cast<AliAODMCParticle*>(mcArray->At(labelEleFromOmegac0));
          Int_t IndexOmegac0 = mcEleFromOmegac0->GetMother();
          mcOmegac0 = static_cast<AliAODMCParticle*>(mcArray->At(IndexOmegac0));
            fVar_Omegac0[30] = mcOmegac0->Pt();
            
        } // lab_Omegac0>=0
    } // fUseMCInfo
    
    KFParticle kfpPion_Rej;
    kfpPion_Rej = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackKaonFromOmega, -211); // -211 pion- for Omega analysis
    KFParticle kfpCasc_Rej;
    
    const KFParticle *vCasc_Rej_Ds[2] = {&kfpPion_Rej, &kfpLambda_m};
    kfpCasc_Rej.Construct(vCasc_Rej_Ds, 2);
    Float_t massCasc_Rej, err_massCasc_Rej;
    kfpCasc_Rej.GetMass(massCasc_Rej, err_massCasc_Rej); // massCasc_Rej
    fVar_Omegac0[31] = massCasc_Rej;
    fVar_Omegac0[32] = cosoa;
    fVar_Omegac0[33] = (Int_t) Convee_ULS + 2 *  (Int_t)Convee_LS;
    fVar_Omegac0[34] = decaytype;
   
 
    if (fWriteOmegac0Tree)
    fTree_Omegac0 -> Fill();
    
    //=================== QA tree ==============
    
    fVar_Omegac0_QA[0]  = trackElectronFromOmegac0 -> GetTPCsignalN();
    fVar_Omegac0_QA[1]  = trackElectronFromOmegac0 -> GetITSNcls();
    fVar_Omegac0_QA[2]  = trackElectronFromOmegac0 -> GetTPCNCrossedRows();
    fVar_Omegac0_QA[3]  = trackElectronFromOmegac0 -> GetTPCNclsF();
    fVar_Omegac0_QA[4]  = trackKaonFromOmega -> GetTPCNCrossedRows();
    fVar_Omegac0_QA[5]  = trackKaonFromOmega -> GetTPCNclsF();
    fVar_Omegac0_QA[6]  = trackKaonFromOmega -> GetTPCsignalN();
    fVar_Omegac0_QA[7]  = trkProton -> GetTPCNCrossedRows();
    fVar_Omegac0_QA[8]  = trkProton -> GetTPCNclsF();
    fVar_Omegac0_QA[9]  = trkPion -> GetTPCNCrossedRows();
    fVar_Omegac0_QA[10] = trkPion -> GetTPCNclsF();
    fVar_Omegac0_QA[11] = casc->DcaBachToPrimVertex();
    fVar_Omegac0_QA[12] = casc->DcaV0ToPrimVertex();
   
    Bool_t isparticle = kTRUE;
    if (casc -> ChargeXi() >0) isparticle = kFALSE;
    if (isparticle){
        fVar_Omegac0_QA[13] = casc->DcaPosToPrimVertex();
        fVar_Omegac0_QA[14] = casc->DcaNegToPrimVertex();
    }
    else {
        fVar_Omegac0_QA[14] = casc->DcaPosToPrimVertex();
        fVar_Omegac0_QA[13] = casc->DcaNegToPrimVertex();
    }
    
    fVar_Omegac0_QA[15] = casc -> DcaXiDaughters();
    fVar_Omegac0_QA[16] =  casc -> DcaV0Daughters();
    
    Double_t lPosOmega[3];
    lPosOmega[0] = casc->DecayVertexXiX();
    lPosOmega[1] = casc->DecayVertexXiY();
    lPosOmega[2] = casc->DecayVertexXiZ();
    Double_t decayvertOmega = TMath::Sqrt(lPosOmega[0]*lPosOmega[0]+lPosOmega[1]*lPosOmega[1]);
    
    Double_t lPosV0[3];
    lPosV0[0] = casc->DecayVertexV0X();
    lPosV0[1] = casc->DecayVertexV0Y();
    lPosV0[2] = casc->DecayVertexV0Z();
    Double_t decayvertV0 = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
   
    Double_t V0CosineOfPointingAngleOmega = casc->CosPointingAngle(lPosOmega);
    fpVtx = (AliAODVertex*)aodEvent->GetPrimaryVertex();
    Double_t primvert[3];
    fpVtx->GetXYZ(primvert);
    Double_t OmegaCosineOfPointingAngle = casc->CosPointingAngleXi(primvert[0],primvert[1],primvert[2]);
    
    fVar_Omegac0_QA[17] = decayvertOmega;  // CascDecayLength
    fVar_Omegac0_QA[18] = decayvertV0;     // V0DecayLength
    fVar_Omegac0_QA[19] = V0CosineOfPointingAngleOmega; // V0CosineOfPointingAngleOmega
    fVar_Omegac0_QA[20] = OmegaCosineOfPointingAngle; // OmegaCosineOfPointingAngle
    
    fVar_Omegac0_QA[21] = trackElectronFromOmegac0 -> Pt();
    fVar_Omegac0_QA[22] = trkPion -> Pt();
    fVar_Omegac0_QA[23] = trkProton -> Pt();
    fVar_Omegac0_QA[24] = trackKaonFromOmega -> Pt();
    
    //----- in the follwoing branches used for cross check with fTree_Omegac tree
    //----- with mass const.
    Float_t mass_Lam_const, err_mass_Lam_const;
    kfpLambda_m.GetMass(mass_Lam_const, err_mass_Lam_const);
    fVar_Omegac0_QA[25] = mass_Lam_const;  // mass of Lambda with mass const.
    
    Float_t mass_Omega_const, err_mass_Omega_const;
    kfpOmegaMinus_m.GetMass(mass_Omega_const, err_mass_Omega_const);
    fVar_Omegac0_QA[26] = mass_Omega_const; // Omega mass with mass const.
    
    
    fVar_Omegac0_QA[27] = kfpLambda_m.GetChi2()/kfpLambda_m.GetNDF(); // Chi2_MassConst of Lambda
    fVar_Omegac0_QA[28] = kfpOmegaMinus_m.GetChi2()/kfpOmegaMinus_m.GetNDF(); // Chi2_MassConst of Omega
    
    Float_t mass_Omegac0_woMassConst_PV, err_mass_Omegac0_woMassConst_PV;
    kfpOmegac0_woMassConst_PV.GetMass(mass_Omegac0_woMassConst_PV, err_mass_Omegac0_woMassConst_PV);
    fVar_Omegac0[29] = mass_Omegac0_woMassConst_PV; //  without Omega mass const.
    //----- without mass const.
    fVar_Omegac0_QA[30] = kfpOmegac0_woMassConst.GetPt(); // pt of EleOmega without Omega mass const.
    fVar_Omegac0_QA[31] = kfpOmegac0_woMassConst.GetRapidity(); // rap of EleOmega without Omega mass const.
    fVar_Omegac0_QA[32] = kfpOmegac0_woMassConst_PV.GetChi2()/kfpOmegac0_woMassConst_PV.GetNDF();  // without Omega mass const.
   
    if(fWriteOmegac0QATree) fTree_Omegac0_QA->Fill();

    return;

}//FillTreeRecOmegac0FromCasc
//____________________________________________________________________________
void AliAnalysisTaskSESemileptonicOmegac0KFP :: DefineEvent()
{
    // This is used to define tree variables
    const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
    fTree_Event = new TTree(nameoutput, "Event");
    
    Int_t nVar = 7;
    fVar_Event = new Float_t[nVar];
    TString *fVarNames = new TString[nVar];

    fVarNames[0]  = "centrality";
    fVarNames[1]  = "z_vtx_reco";
    fVarNames[2]  = "n_vtx_contributors";
    fVarNames[3]  = "n_tracks";
    fVarNames[4]  = "is_ev_rej";
    fVarNames[5]  = "run_number";
    fVarNames[6]  = "ev_id";

    for (Int_t ivar=0; ivar<nVar; ivar++) {
      fTree_Event->Branch(fVarNames[ivar].Data(), &fVar_Event[ivar], Form("%s/F", fVarNames[ivar].Data()));
    }

    return;

}  // DefineEvent()
//____________________________________________________________________________
void AliAnalysisTaskSESemileptonicOmegac0KFP :: DefineTreeRecoOmegac0()
{
 
    // This is to define tree variabels
    
    const char* nameoutput = GetOutputSlot(5)->GetContainer()->GetName();
    fTree_Omegac0 = new TTree(nameoutput, "Omegac0 variables tree");
    Int_t nVar = 35;
    fVar_Omegac0 = new Float_t[nVar];
    TString *fVarNames = new TString[nVar];
    
    fVarNames[0]  = "nSigmaCombined_Ele"; // nsigma_combined for electron
    fVarNames[1]  = "nSigmaCombined_Kaon"; // nsigma_combined for kaon
    fVarNames[2]  = "nSigmaTPC_PiFromLam"; // TPC nSigma for pion from Lambda
    fVarNames[3]  = "nSigmaTPC_PrFromLam"; // TPC nSigma for proton from Lambda
    fVarNames[4]  = "DCA_OmegaDau"; // DCA of Omega's daughters (calculated from AOD cascade)
    fVarNames[5]  = "chi2geo_Lam"; // chi2_geometry of Lambda (without Lambda mass const.)
    fVarNames[6]  = "ldl_Lam"; // l/dl of Lambda (without Lambda mass const. and without PV const.)
    fVarNames[7] = "chi2geo_Omega"; // chi2_geometry of Omega (with Lambda mass const.)
    fVarNames[8] = "ldl_Omega"; // l/dl of Omega (with Lambda mass const.)
    fVarNames[9] = "chi2topo_OmegaToPV"; // chi2_topo of Omega to PV  (with Omega mass const.)
    fVarNames[10] = "DecayLxy_Lambda"; // decay length of Lambda in x-y plane (with Lambda mass const.)
    fVarNames[11] = "DecayLxy_Omega";  // decay length of Omega in x-y plane (with Omega mass const.)
    fVarNames[12] = "PA_LamToOmega"; // pointing angle of Lmabda (pointing back to Omega) (with Lambda mass const.)
    fVarNames[13] = "PA_OmegaToPV"; // pointing angle of Omega (pointing back to PV) (with Omega mass const.)
    fVarNames[14] = "Mass_Lam"; // mass of Lambda (without Lambda mass const.)
    fVarNames[15] = "Mass_Omega"; // mass of Omega (without Omega mass const.)
    fVarNames[16] = "Pt_EleFromOmegac0"; // pt of electron from Omegac0
    fVarNames[17] = "Pt_EleOmega_KF"; // pt of EleOmega_pairs (with Omega mass const. and with PV const.)
    fVarNames[18] = "Rap_EleOmega_KF"; // rapidity of EleOmega_pairs  (with Omega mass const. and with PV const.)
    fVarNames[19] = "Mass_EleOmega_KF"; // mass of EleOmega_pairs  (with Omega mass const. and with PV const.)
    fVarNames[20] = "Chi2prim_EleFromOmegac0"; // chi2_topo of electron to PV
    fVarNames[21] = "DCAxy_EleFromOmegac0_KF"; // DCA of electrono coming from Omegac0 in x-y plane
    fVarNames[22] = "Mass_K0S"; // mass of Ks0
    fVarNames[23] = "Chi2topo_LamToOmega"; // chi2_topo of Lambda to Omega (with Lambda mass const. and PV const.)
    fVarNames[24] = "Chi2geo_Omegac0"; // chi2_geometry of Omegac0 (with Omega mass const.)
    fVarNames[25] = "DCA_LamDau"; // DCA of Lambda's daughters (calculated from AOD cascade)
    fVarNames[26] = "DCA_OmegaDau_KF"; // DCA of Omega's daughters (calculated from KF after Lambda mass constraint)
    fVarNames[27] = "DCA_Omegac0Dau_KF"; // DCA of Omegac0's daughters (calculated from KF after Omega mass constraint)
    fVarNames[28] = "DCAxy_OmegaToPV_KF"; // DCA of Omega to PV in x-y plane (with Omega mass constraint)
    fVarNames[29] = "Source_Omegac0"; // flag for Omegac0 MC (“4”:prompt, "5": feed-down, “<0”: background)
    fVarNames[30] = "Omegac0_Pt_MC"; // Omegac0 pt distribution in MC
    fVarNames[31] = "Mass_Xi"; // mass of Xi, (with Lambda mass const.) - - to reject
    fVarNames[32] = "EleOmegaOA";  // Calcualtion from AOD
    fVarNames[33] = "ConvType"; // ee pairs - - prefilter method
    fVarNames[34] = "DecayType";  // flags for WS and RS of EleOmega_pairs
    
    
    for (Int_t ivar = 0; ivar<nVar ; ivar++){
     
        fTree_Omegac0->Branch(fVarNames[ivar].Data(), &fVar_Omegac0[ivar], Form("%s/F", fVarNames[ivar].Data()));
    }
    
    return;
    
} // DefineTree_RecoOmegac0()
//____________________________________________________________________________
void AliAnalysisTaskSESemileptonicOmegac0KFP :: DefineTreeRecoOmegac0_QA()
{
    
    const char* nameoutput = GetOutputSlot(7)->GetContainer()->GetName();
    fTree_Omegac0_QA = new TTree(nameoutput, "QA of Omegac varibales tree");
    Int_t nVar = 33;
    fVar_Omegac0_QA = new Float_t[nVar];
    TString *fVarNames = new TString[nVar];
    
    fVarNames[0]  = "Ele_TPCPID";
    fVarNames[1]  = "Ele_ITS";
    fVarNames[2]  = "Ele_CrossedRows";
    fVarNames[3]  = "Ele_FindableRatio";
    fVarNames[4]  = "Kaon_CrossedRows";
    fVarNames[5]  = "Kaon_FindableRatio";
    fVarNames[6]  = "Kaon_TPCPID";
    fVarNames[7]  = "Proton_CrossedRows";
    fVarNames[8]  = "Proton_FindableRatio";
    fVarNames[9]  = "Pion_CrossedRows";
    fVarNames[10] = "Pion_FindableRatio";
    fVarNames[11] = "DCABachToPrimVertex";
    fVarNames[12] = "DcaV0ToPrimVertex";
    fVarNames[13] = "DCAV0PosToPrimVertex";
    fVarNames[14] = "DCAV0NegToPrimVertex";
    fVarNames[15] = "DCACascDaughters";
    fVarNames[16] = "DCAV0Daughters";
    fVarNames[17] = "CascDecayLength";
    fVarNames[18] = "V0DecayLength";
    fVarNames[19] = "V0CosineOfPoiningAngleOmega";
    fVarNames[20] = "OmegaCosineOfPoiningAngle";
    fVarNames[21] = "Ele_pt";
    fVarNames[22] = "Pion_pt";
    fVarNames[23] = "Proton_pt";
    fVarNames[24] = "Kaon_pt";
    fVarNames[25] = "Mass_Lambda"; // with lam. mass const.
    fVarNames[26] = "Mass_Omega"; // with Omega mass const.
    fVarNames[27] = "Chi2geo_Lam_WMC"; // with Lam mass const.
    fVarNames[28] = "Chi2geo_Omega_WMC"; // with Omega mass const.
    fVarNames[29] = "Mass_EleOmega_KFP_WoMC"; // without Omega mass const
    fVarNames[30] = "Pt_EleOmega_KFP_WoMC"; // without Omega mass constraint
    fVarNames[31] = "Rap_EleOmega_KFP_WoMC"; //without  Omega mass constraint
    fVarNames[32] = "Chi2topo_Omegac0_WoMC"; // topo_Omegac0 without Omega mass constraint
   
    
    for (Int_t ivar = 0; ivar<nVar; ivar++){
        fTree_Omegac0_QA->Branch(fVarNames[ivar].Data(), &fVar_Omegac0_QA[ivar], Form("%s/F", fVarNames[ivar].Data()));
    }
    
    return;
                
} //DefineTreeRecoOmegac0_QA
