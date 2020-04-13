
//
// Calculate flow in the forward and central regions using the Q cumulants method.
//
// Inputs:
//  - AliAODEvent
//
// Outputs:
//  - AnalysisResults.root or forward_flow.root
//
#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TList.h>
#include <THn.h>
#include "AliAODMCHeader.h"
#include "AliLog.h"
#include "AliForwardNUETask.h"
#include "AliForwardGenericFramework.h"
#include "AliForwardNUATask.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliForwardFlowUtil.h"

#include "AliVVZERO.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"

#include "AliESDEvent.h"
#include "AliVTrack.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"

#include "AliAnalysisFilter.h"
#include "AliMultSelection.h"
#include "AliMultiplicity.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "TMath.h"
#include "AliStack.h"

using namespace std;
ClassImp(AliForwardNUETask)
#if 0
; // For emacs
#endif

//_____________________________________________________________________
AliForwardNUETask::AliForwardNUETask() : AliAnalysisTaskSE(),
  fAOD(0),           // input event
  fOutputList(0),    // output list
  fEventList(0),
  fSettings(),
  fUtil(),
  nua_mode(kFALSE)
  {
  //
  //  Default constructor
  //
  }

  //_____________________________________________________________________
  AliForwardNUETask::AliForwardNUETask(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0),           // input event
  fOutputList(0),
  fEventList(0),
  fSettings(),
  fUtil(),
  nua_mode(kFALSE)
  {
  //
  //  Constructor
  //
  //  Parameters:
  //   name: Name of task
  //
  DefineInput(1, AliForwardTaskValidation::Class());

    DefineOutput(1, TList::Class());
  }

//_____________________________________________________________________
  void AliForwardNUETask::UserCreateOutputObjects()
  {
  //
  //  Create output objects
  //
    fOutputList = new TList();          // the final output list
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them if requested


    //..adding QA plots from AliEventCuts class

    fEventList = new TList();

    fEventList->Add(new TH1D("Centrality","Centrality",100,0,100));
    fEventList->Add(new TH1D("Vertex","Vertex",fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));
    fEventList->Add(new TH2F("hOutliers","Maximum #sigma from mean N_{ch} pr. bin",
     20, 0., 100., 500, 0., 5.)); //((fFlags & kMC) ? 15. : 5. // Sigma <M> histogram
    fEventList->Add(new TH1D("FMDHits","FMDHits",100,0,10));

    fEventList->SetName("EventInfo");

    fOutputList->Add(new TH2F("NUA_fmd","NUA_fmd", 200, -4.0, 6.0,fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));
    fOutputList->Add(new TH2F("NUA_fmd_prim","NUA_fmd_prim", 200, -4.0, 6.0,fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));

    fOutputList->Add(new TH2F("NUA_spd","NUA_spd", 400, -2.5, 2.5,fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));
    fOutputList->Add(new TH2F("NUA_spd_prim","NUA_spd_prim", 400, -2.5, 2.5,fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));

    // create hist for tpc (eta, pt, z, filterbit )
    Int_t dimensions = 4;
    Int_t bins[4] = {300, 100, fSettings.fNZvtxBins, 5} ;
    Double_t xmin[5] = {-1.2, 0.0, fSettings.fZVtxAcceptanceLowEdge, -0.5};
    Double_t xmax[5] = {1.2, 5.0, fSettings.fZVtxAcceptanceUpEdge, 4.5};

    fOutputList->Add(new THnD("NUA_tpc", "NUA_tpc", dimensions, bins, xmin, xmax)); //(eta, n)
    fOutputList->Add(new TH3F("NUA_tpc_prim","NUA_tpc_prim", 300, -1.2, 1.2, 100, 0.0, 5.0,fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));

    fOutputList->Add(fEventList);

    PostData(1, fOutputList);
  }


//_____________________________________________________________________
void AliForwardNUETask::UserExec(Option_t *)
{
  AliForwardTaskValidation* ev_val = dynamic_cast<AliForwardTaskValidation*>(this->GetInputData(1));

  //
  //  Analyses the event with use of the helper class AliForwardQCumulantRun2
  //
  //  Parameters:
  //   option: Not used
  //
  AliAODEvent* fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  //AliAODMCHeader* fAODMCHeader = static_cast<AliAODMCHeader*>(fAOD->FindListObject(AliAODMCHeader::StdBranchName()));

  AliMCEvent* fAODMC = this->MCEvent();
  //AliStack* stack = fAODMC->Stack();
  //..check if I have AOD

  if(!fAOD){
    Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    return;
  }

  //..AliEventCuts selection
  if (!fSettings.esd && !ev_val->IsValidEvent()){
    PostData(1, this->fOutputList);
    return;
  }

  // Get detector objects
  AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));

  //TList* eventList = static_cast<TList*>(fOutputList->FindObject("EventInfo"));
  TH2F* nua_fmd = static_cast<TH2F*>(fOutputList->FindObject("NUA_fmd"));
  THnD* nua_tpc = static_cast<THnD*>(fOutputList->FindObject("NUA_tpc"));
  TH2F* nua_spd = static_cast<TH2F*>(fOutputList->FindObject("NUA_spd"));
  TH2F* nua_fmd_prim = static_cast<TH2F*>(fOutputList->FindObject("NUA_fmd_prim"));
  TH3F* nua_tpc_prim = static_cast<TH3F*>(fOutputList->FindObject("NUA_tpc_prim"));
  TH2F* nua_spd_prim = static_cast<TH2F*>(fOutputList->FindObject("NUA_spd_prim"));
  //TH1D* fFMDHits = static_cast<TH1D*>(eventList->FindObject("FMDHits"));

  Int_t  iTracks(fAOD->GetNumberOfTracks());
  //AliMultSelection *MultSelection = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  //Double_t cent = MultSelection->GetMultiplicityPercentile("V0M"); //CL0

  TH2D& forwarddNdedp = aodfmult->GetHistogram();
  Double_t zvertex = fAOD->GetPrimaryVertex()->GetZ();

  // extra cut on the FMD
  //if (!fUtil.ExtraEventCutFMD(forwarddNdedp, cent, true)) useEvent = false;

    // loop for the SPD
    AliAODTracklets* aodTracklets = fAOD->GetTracklets();
    std::cout << "aodTracklets->GetNumberOfTracklets()" << aodTracklets->GetNumberOfTracklets() << std::endl;
    for (Int_t i = 0; i < aodTracklets->GetNumberOfTracklets(); i++) {
      std::cout << aodTracklets->GetEta(i)<<std::endl;
      nua_spd->Fill(aodTracklets->GetEta(i),zvertex,1);
    }

              //Double_t x[5] = {noSamples, zvertex, refEtaA, cent, fSettings.kW4Four};
              //x[4] = fSettings.kW2Two;
    // loop for the TPC
    for(Int_t i(0); i < iTracks; i++) {
      AliAODTrack* track = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(i));
      //if (track->Pt() >= 0.2 && track->Pt() <= 5){
        if (fabs(track->Eta()) > 1.1) continue;
        Double_t x[4] = {track->Eta(),track->Pt(),zvertex,0.0};

        if (track->TestFilterBit(fSettings.kTPCOnly)){
          nua_tpc->Fill(x,1);
        }
        if (track->TestFilterBit(fSettings.kHybrid)){
          x[3] = 1.0;
          nua_tpc->Fill(x,1);
        }
        if (track->TestFilterBit(fSettings.kGlobal)){
          x[3] = 2.0;
          nua_tpc->Fill(x,1);
        }
        if (track->TestFilterBit(fSettings.kGlobalComb)){
          x[3] = 3.0;
          nua_tpc->Fill(x,1);
        }
        if (track->TestFilterBit(fSettings.kGlobalLoose)){
          x[3] = 4.0;
          nua_tpc->Fill(x,1);
        }
      //}
    }

    // loop for the FMD
    Int_t phibins = forwarddNdedp.GetNbinsY();
    for (Int_t etaBin = 1; etaBin <= forwarddNdedp.GetNbinsX(); etaBin++) {
        Double_t eta = forwarddNdedp.GetXaxis()->GetBinCenter(etaBin);
        for (Int_t phiBin = 1; phiBin <= phibins; phiBin++) {

          if (forwarddNdedp.GetBinContent(etaBin, 0) == 0) break;

          Double_t weight = forwarddNdedp.GetBinContent(etaBin, phiBin);

          // If empty, do not fill hist
          if (weight == 0) continue;
          nua_fmd->Fill(eta,zvertex,weight);
        } // End of phi loop
      } // End of eta bin

    Int_t nTracksMC   = fAODMC->GetNumberOfTracks();

    for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
      AliAODMCParticle* p = static_cast< AliAODMCParticle* >(fAODMC->GetTrack(iTr));
      if (!p->IsPhysicalPrimary()) continue;
      if (p->Charge() == 0) continue;

      if ( (p->Eta() < -1.7 && p->Eta() > -3.4) || (p->Eta() > 1.7 && p->Eta() < 5.0) ){
        nua_fmd_prim->Fill(p->Eta(),zvertex);
      }
      //if (p->Pt()>=0.2 && p->Pt()<=5) {
        if (fabs(p->Eta()) < 1.1) {
          nua_tpc_prim->Fill(p->Eta(),p->Pt(),zvertex);
        }
      //}
      if (fabs(p->Eta()) < 2.5) {
        nua_spd_prim->Fill(p->Eta(),zvertex);
      }
    }

  PostData(1, fOutputList);
  return;
}


//_____________________________________________________________________
void AliForwardNUETask::Terminate(Option_t */*option*/)
{
  return;
}


//_____________________________________________________________________
//
//
// EOF
