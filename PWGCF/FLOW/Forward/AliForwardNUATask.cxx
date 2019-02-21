
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
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliLog.h"
#include "AliForwardNUATask.h"
#include "AliForwardQCumulantRun2.h"
#include "AliForwardGenericFramework.h"
#include "AliForwardFlowRun2Task.h"
#include "TObjectTable.h"
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
  refDist(),
  forwardDist(),
  fSettings(),
  fUtil(),
  useEvent(kTRUE)
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
  refDist(),
  forwardDist(),
  ev_val(),
  fSettings(),
  fUtil(),
  useEvent(kTRUE)
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
  void AliForwardNUATask::UserCreateOutputObjects()
  {

  //
  //  Create output objects
  //
    fOutputList = new TList();          // the final output list
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them if requested


    fEventList = new TList();
    //fEventList->Add(new TH1D("ImpactParam","ImpactParam",100,0,20));
    fEventList->Add(new TH1D("Centrality","Centrality",100,0,100));
    fEventList->Add(new TH1D("Vertex","Vertex",fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));
    fEventList->Add(new TH1D("EventCuts_FMD","EventCuts_FMD",3,0,3));
    fEventList->Add(new TH2D("hOutliers","Maximum #sigma from mean N_{ch} pr. bin",
       20, 0., 100., 500, 0., 5.)); //((fFlags & kMC) ? 15. : 5. // Sigma <M> histogram
    fEventList->SetName("EventInfo");

    fEventList->Add(new TH1F("dNdeta","dNdeta",100 /*fSettings.fNDiffEtaBins*/,fSettings.fEtaLowEdge,fSettings.fEtaUpEdge));


    Double_t centralEta = (fSettings.useSPD ? 2.5 : 1.5);
    Int_t forwardBinsEta = (fSettings.use_primaries ? 200 : 200);
    Int_t forwardBinsPhi = (fSettings.use_primaries ? 20 : 20);

    fOutputList->Add(new TH3D("NUA_fwd","NUA_fwd", forwardBinsEta, -4.0, 6.0, forwardBinsPhi, 0., 2*TMath::Pi(),fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));

    fOutputList->Add(new TH3D("NUA_cen","NUA_cen", 400, -centralEta, centralEta, 400, 0., 2*TMath::Pi(),fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));

    fOutputList->Add(fEventList);

    PostData(1, fOutputList);
  }


//_____________________________________________________________________
void AliForwardNUATask::UserExec(Option_t *)
{


 //gObjectTable->Print();

  //
  //  Analyses the event with use of the helper class AliForwardQCumulantRun2
  //
  //  Parameters:
  //   option: Not used
  //
  // Get the event validation object
  ev_val = dynamic_cast<AliForwardTaskValidation*>(this->GetInputData(1));
  if (!ev_val->IsValidEvent()){
     PostData(1, this->fOutputList);
    return;
  }
  if (!fSettings.esd) {
    AliAODEvent* aodevent = dynamic_cast<AliAODEvent*>(InputEvent());
    fUtil.fAODevent = aodevent;
    if(!aodevent) throw std::runtime_error("Not AOD as expected");
  }

  fUtil.fevent = fInputEvent;
  fUtil.fSettings = fSettings;
  if (fSettings.mc) fUtil.fMCevent = this->MCEvent();


  Double_t centralEta = (fSettings.useSPD ? 2.5 : 1.5);
  TH2D centralDist_tmp = TH2D("c","",400,-centralEta,centralEta,400,0,2*TMath::Pi());
  centralDist_tmp.SetDirectory(0);
  TH2D refDist_tmp = TH2D("c","",400,-centralEta,centralEta,400,0,2*TMath::Pi());
  refDist_tmp.SetDirectory(0);

  TH2D forwardTrRef  ("ft","",200,-4,6,20,0,TMath::TwoPi());
  TH2D forwardPrim  ("fp","",400,-4,6,400,0,TMath::TwoPi());
  forwardTrRef.SetDirectory(0);
  forwardPrim.SetDirectory(0);
  forwardDist = (fSettings.use_primaries_fwd ? &forwardPrim : &forwardTrRef);

  centralDist = &centralDist_tmp;
  centralDist->SetDirectory(0);
  refDist = &refDist_tmp;
  refDist->SetDirectory(0);

  TH1F* dNdeta = static_cast<TH1F*>(fEventList->FindObject("dNdeta"));
  dNdeta->SetDirectory(0);

  fUtil.dodNdeta = kTRUE;
  fUtil.dNdeta = dNdeta;
  fUtil.FillData(refDist,centralDist,forwardDist);

  forwardDist->SetDirectory(0);

  Double_t zvertex = fUtil.GetZ();
  Double_t cent = fUtil.GetCentrality(fSettings.centrality_estimator);

  TList* eventList = static_cast<TList*>(fOutputList->FindObject("EventInfo"));
  TH3D* nua_fmd = static_cast<TH3D*>(fOutputList->FindObject("NUA_fwd"));
  TH3D* nua_tpc = static_cast<TH3D*>(fOutputList->FindObject("NUA_cen"));
  nua_tpc->SetDirectory(0);
  nua_fmd->SetDirectory(0);

  if (useEvent) {
    static_cast<TH1D*>(fEventList->FindObject("Vertex"))->Fill(zvertex);
    static_cast<TH1D*>(fEventList->FindObject("Centrality"))->Fill(cent);

    // loop for the TPC

    for (Int_t etaBin = 1; etaBin <= centralDist->GetNbinsX(); etaBin++) {
      for (Int_t phiBin = 1; phiBin <= centralDist->GetNbinsY(); phiBin++) {
        if (centralDist->GetBinContent(etaBin,phiBin) == 0) continue;//TH2D
        //nua_tpc->SetBinContent(etaBin,phiBin,nua_tpc->GetZaxis()->FindBin(zvertex),centralDist->GetBinContent(etaBin,phiBin));
        nua_tpc->Fill(centralDist->GetXaxis()->GetBinCenter(etaBin),centralDist->GetYaxis()->GetBinCenter(phiBin),zvertex,centralDist->GetBinContent(etaBin,phiBin)); //TH3D
        //nua_tpc->AddBinContent(nua_fmd->GetBin(etaBin,phiBin,nua_tpc->GetZaxis()->FindBin(zvertex)),centralDist->GetBinContent(etaBin,phiBin));
      }
    }

    if (fSettings.makeFakeHoles) fUtil.MakeFakeHoles(*forwardDist);
    // loop for the FMD
    for (Int_t etaBin = 1; etaBin <= forwardDist->GetNbinsX(); etaBin++) {
      for (Int_t phiBin = 1; phiBin <= forwardDist->GetNbinsY(); phiBin++) {
        //if (forwardDist->GetBinContent(etaBin, 0) == 0) break;
        Double_t weight = forwardDist->GetBinContent(etaBin,phiBin);

        if (fSettings.nua_mode & fSettings.kInterpolate)
          weight = InterpolateWeight(*forwardDist,phiBin,etaBin,weight);

        if (weight == 0) continue;

        nua_fmd->Fill(forwardDist->GetXaxis()->GetBinCenter(etaBin),forwardDist->GetYaxis()->GetBinCenter(phiBin),zvertex,weight);
        //nua_fmd->SetBinContent(etaBin,phiBin,nua_fmd->GetZaxis()->FindBin(zvertex),weight);
        //nua_fmd->AddBinContent(nua_fmd->GetBin(etaBin,phiBin,nua_fmd->GetZaxis()->FindBin(zvertex)),weight);
      }
    }
  }
  PostData(1, fOutputList);
  return;
}



Double_t AliForwardNUATask::InterpolateWeight(const TH2D& forwarddNdedp,Int_t phiBin, Int_t etaBin, Double_t weight)
{
  if ((phiBin == 17) && (etaBin >= 125 && etaBin <= 137)){

    if (weight > 0 || forwarddNdedp.GetBinContent(etaBin, 18) > 0) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin - 1, phiBin) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin - 2, phiBin) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin - 3, phiBin) > 0 ) return weight;

    // std::cout << "found hole, etaBin = " << etaBin << ", eta = " << eta << ", phiBin = " << phiBin << ", phi = " << phi << std::endl;
 //    if (detType == "forward") weight = 1.;
     Double_t up = forwarddNdedp.GetBinContent(etaBin, 19);
     Double_t low = forwarddNdedp.GetBinContent(etaBin, 16);
     weight = ((up+low)/2+low)/2;
     return weight;
    //std::cout << weight << std::endl;

   }
  if ((phiBin == 18) && (etaBin >= 125 && etaBin <= 137)){
    if (weight > 0 || forwarddNdedp.GetBinContent(etaBin, 17) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin - 1, phiBin) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin - 2, phiBin) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin - 3, phiBin) > 0 ) return weight;

    // std::cout << "found hole, etaBin = " << etaBin << ", eta = " << eta << ", phiBin = " << phiBin << ", phi = " << phi << std::endl;
     Double_t up = forwarddNdedp.GetBinContent(etaBin, 19);
     Double_t low = forwarddNdedp.GetBinContent(etaBin, 16);
     weight = ((up+low)/2+up)/2;
    //std::cout << weight << std::endl;
  }
  if ((phiBin == 14) && (etaBin >= 168 && etaBin <= 185)){
    if (weight > 0) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin + 1, phiBin) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin + 2, phiBin) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin + 3, phiBin) > 0 ) return weight;
    // std::cout << "found hole, etaBin = " << etaBin << ", eta = " << eta << ", phiBin = " << phiBin << ", phi = " << phi << std::endl;
    Double_t  up = forwarddNdedp.GetBinContent(etaBin, 15);
    Double_t  low = forwarddNdedp.GetBinContent(etaBin, 13);
     weight = (up+low)/2;
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
