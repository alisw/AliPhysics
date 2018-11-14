
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

#include "AliLog.h"
#include "AliForwardNUATask.h"
#include "AliForwardQCumulantRun2.h"
#include "AliForwardGenericFramework.h"
#include "AliForwardFlowRun2Task.h"

#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"

#include "AliForwardFlowUtil.h"

#include "AliVVZERO.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"

#include "AliESDEvent.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"

#include "AliAnalysisFilter.h"
#include "AliMultSelection.h"
#include "AliMultiplicity.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "TMath.h"

using namespace std;
ClassImp(AliForwardNUATask)
#if 0
; // For emacs
#endif

//_____________________________________________________________________
AliForwardNUATask::AliForwardNUATask() : AliAnalysisTaskSE(),
  fAOD(0),           // input event
  fOutputList(0),    // output list
  fEventList(0),
  centralDist(),
  forwardDist(),
  fSettings(),
  fUtil(),
  useEvent(kTRUE),
  fEventCuts()
  {
  //
  //  Default constructor
  //
  }

  //_____________________________________________________________________
  AliForwardNUATask::AliForwardNUATask(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0),           // input event
  fOutputList(0),
  fEventList(0),
  centralDist(),
  forwardDist(),
  fSettings(),
  fUtil(),
  useEvent(kTRUE),
  fEventCuts()
  {
  //
  //  Constructor
  //
  //  Parameters:
  //   name: Name of task
  //
    DefineOutput(1, TList::Class());
  }

//_____________________________________________________________________
  void AliForwardNUATask::UserCreateOutputObjects()
  {
  //
  //  Create output objects
  //
    fOutputList = new TList();          // the final output list
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them if requested

    //..adding QA plots from AliEventCuts class
    fEventCuts.AddQAplotsToList(fOutputList);

    fEventList = new TList();

    fEventList->Add(new TH1D("Centrality","Centrality",100,0,100));
    fEventList->Add(new TH1D("Vertex","Vertex",fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));

    fEventList->SetName("EventInfo");

    Double_t centralEta = (fSettings.useSPD ? 2.5 : 1.5);
    Int_t forwardBinsEta = (fSettings.use_primaries ? 400 : 200);
    Int_t forwardBinsPhi = (fSettings.use_primaries ? 400 : 20);

    fOutputList->Add(new TH3F("NUA_fwd","NUA_fwd", forwardBinsEta, -4.0, 6.0, forwardBinsPhi, 0., 2*TMath::Pi(),fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));

    fOutputList->Add(new TH3F("NUA_cen","NUA_cen", 400, -centralEta, centralEta, 400, 0., 2*TMath::Pi(),fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));

    fOutputList->Add(fEventList);

    PostData(1, fOutputList);
  }


//_____________________________________________________________________
void AliForwardNUATask::UserExec(Option_t *)
{
  //
  //  Analyses the event with use of the helper class AliForwardQCumulantRun2
  //
  //  Parameters:
  //   option: Not used
  //
  fUtil.fevent = fInputEvent;

  Double_t centralEta = (fSettings.useSPD ? 2.5 : 1.5);
  TH2D centralDist_tmp = TH2D("c","",400,-centralEta,centralEta,400,0,2*TMath::Pi());
  centralDist_tmp.SetDirectory(0);

  TH2D forwardTrRef  ("ft","",200,-4,6,20,0,TMath::TwoPi());
  TH2D forwardPrim  ("fp","",400,-4,6,400,0,TMath::TwoPi());
  forwardTrRef.SetDirectory(0);
  forwardPrim.SetDirectory(0);

  centralDist = &centralDist_tmp;
  centralDist->SetDirectory(0);

  if (!fSettings.mc) {

    AliAODEvent* aodevent = dynamic_cast<AliAODEvent*>(InputEvent());

    if(!aodevent)
      throw std::runtime_error("Not AOD as expected");

    //..AliEventCuts selection
    if (!fSettings.esd && !fEventCuts.AcceptEvent(fInputEvent)) {
      PostData(1, fOutputList);
      return;
    }

    AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(aodevent->FindListObject("Forward"));
    forwardDist = &aodfmult->GetHistogram();

    if (fSettings.useSPD) fUtil.FillFromTracklets(centralDist);
    else                  fUtil.FillFromTracks(centralDist, fSettings.tracktype);
  }
  else {
    AliMCEvent* mcevent = dynamic_cast<AliMCEvent*>(InputEvent());


    if(!mcevent)
      throw std::runtime_error("Not MC as expected");

    forwardDist = (fSettings.use_primaries ? &forwardPrim : &forwardTrRef);
    const AliAODVertex vertex = *((AliAODVertex*)(this->MCEvent()->GetPrimaryVertex()));
    if (!fSettings.use_primaries) fUtil.FillFromTrackrefs(centralDist, forwardDist, vertex);
    else                          fUtil.FillFromPrimaries(centralDist, forwardDist);
  }
  forwardDist->SetDirectory(0);

  Double_t zvertex = fUtil.GetZ();
  Double_t cent = fUtil.GetCentrality(fSettings.centrality_estimator);

  if (!fSettings.use_primaries){
    if (!fUtil.ExtraEventCutFMD(*forwardDist, cent, fSettings.mc)) {
      useEvent = false;
      static_cast<TH1D*>(fEventList->FindObject("EventCuts_FMD"))->Fill(1.0);
    }
  }

  TList* eventList = static_cast<TList*>(fOutputList->FindObject("EventInfo"));
  TH3F* nua_fmd = static_cast<TH3F*>(fOutputList->FindObject("NUA_fwd"));
  TH3F* nua_tpc = static_cast<TH3F*>(fOutputList->FindObject("NUA_cen"));

  if (useEvent) {
    // loop for the TPC
    for (Int_t etaBin = 1; etaBin <= centralDist->GetNbinsX(); etaBin++) {
      for (Int_t phiBin = 1; phiBin <= centralDist->GetNbinsX(); phiBin++) {
        nua_tpc->Fill(centralDist->GetBinCenter(etaBin),centralDist->GetBinCenter(phiBin),zvertex,centralDist->GetBinContent(etaBin,phiBin));
      }
    }

    // loop for the FMD
    for (Int_t etaBin = 1; etaBin <= centralDist->GetNbinsX(); etaBin++) {
      for (Int_t phiBin = 1; phiBin <= centralDist->GetNbinsX(); phiBin++) {
        if (forwardDist->GetBinContent(etaBin, 0) == 0) break;
        Double_t weight = forwardDist->GetBinContent(etaBin,phiBin);
        if (fSettings.nua_mode) weight = InterpolateWeight(*forwardDist,phiBin,etaBin,weight);
        if (weight == 0) continue;

        nua_fmd->Fill(forwardDist->GetBinCenter(etaBin),forwardDist->GetBinCenter(phiBin),zvertex,forwardDist->GetBinContent(etaBin,phiBin));
      }
    }
  }
  PostData(1, fOutputList);
  return;
}

Double_t AliForwardNUATask::InterpolateWeight(const TH2D& forwarddNdedp,Int_t phiBin, Int_t etaBin, Double_t weight)
{
  if ((phiBin == 17) && (etaBin > 125 && etaBin < 137)){

    //if (weight > 0 || forwarddNdedp.GetBinContent(etaBin, 18) > 0) return weight;
    if (weight > 0 || forwarddNdedp.GetBinContent(etaBin, 18) > 0) return weight;
    if (forwarddNdedp.GetBinContent(etaBin - 1, phiBin) > 0 ) return weight;
    if (forwarddNdedp.GetBinContent(etaBin - 2, phiBin) > 0 ) return weight;
    if (forwarddNdedp.GetBinContent(etaBin - 3, phiBin) > 0 ) return weight;

    // std::cout << "found hole, etaBin = " << etaBin << ", eta = " << eta << ", phiBin = " << phiBin << ", phi = " << phi << std::endl;
 //    if (detType == "forward") weight = 1.;
     Double_t up = forwarddNdedp.GetBinContent(etaBin, 19);
     Double_t low = forwarddNdedp.GetBinContent(etaBin, 16);
     weight = ((up+low/2)+low)/2;
    //std::cout << weight << std::endl;

   }
  if ((phiBin == 18) && (etaBin > 125 && etaBin < 137)){
    if (weight > 0 || forwarddNdedp.GetBinContent(etaBin, 17) > 0 ) return weight;
    if (forwarddNdedp.GetBinContent(etaBin - 1, phiBin) > 0 ) return weight;
    if (forwarddNdedp.GetBinContent(etaBin - 2, phiBin) > 0 ) return weight;
    if (forwarddNdedp.GetBinContent(etaBin - 3, phiBin) > 0 ) return weight;

    // std::cout << "found hole, etaBin = " << etaBin << ", eta = " << eta << ", phiBin = " << phiBin << ", phi = " << phi << std::endl;
     Double_t up = forwarddNdedp.GetBinContent(etaBin, 19);
     Double_t low = forwarddNdedp.GetBinContent(etaBin, 16);
     weight = ((up+low/2)+low)/2;
    //std::cout << weight << std::endl;
  }
  if ((phiBin == 14) && (etaBin > 168 && etaBin < 185)){
    if (weight > 0) return weight;
    if (forwarddNdedp.GetBinContent(etaBin + 1, phiBin) > 0 ) return weight;
    if (forwarddNdedp.GetBinContent(etaBin + 2, phiBin) > 0 ) return weight;
    if (forwarddNdedp.GetBinContent(etaBin + 3, phiBin) > 0 ) return weight;
    // std::cout << "found hole, etaBin = " << etaBin << ", eta = " << eta << ", phiBin = " << phiBin << ", phi = " << phi << std::endl;
    Double_t  up = forwarddNdedp.GetBinContent(etaBin, 15);
    Double_t  low = forwarddNdedp.GetBinContent(etaBin, 13);
     weight = (up+low/2);
    //std::cout << weight << std::endl;

   }
   return weight;
}


//_____________________________________________________________________
void AliForwardNUATask::Terminate(Option_t */*option*/)
{
  return;
}


//_____________________________________________________________________
//
//
// EOF
