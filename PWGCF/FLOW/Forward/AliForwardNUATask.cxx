
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
  nua_cen(),
  nua_fmd(),
  dNdeta(),
  ev_val(),
  fSettings(),
  fUtil()
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
  nua_cen(),
  nua_fmd(),
  dNdeta(),
  ev_val(),
  fSettings(),
  fUtil()
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

    fEventList->Add(new TH1D("Vertex","Vertex",fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));
    fEventList->SetName("EventInfo");

    Int_t   centralEtaBins = (fSettings.useITS ? 200 : 300);
    Int_t   centralPhiBins = (fSettings.useITS ? 20 : 300);
    Double_t centralEtaMin = (fSettings.useSPD ? -2.5 : fSettings.useITS ? -4 : -1.2);
    Double_t centralEtaMax = (fSettings.useSPD ? 2.5 : fSettings.useITS ? 6 : 1.2);

    Int_t forwardBinsEta = (fSettings.use_primaries_fwd ? 200 : 200);
    Int_t forwardBinsPhi = (fSettings.use_primaries_fwd ? 20 : 20);

    nua_cen = new TH3D("NUA_cen","NUA_cen", centralEtaBins, centralEtaMin, centralEtaMax, centralPhiBins, 0., 2*TMath::Pi(),fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge);
    nua_fmd = new TH3D("NUA_fwd","NUA_fwd", forwardBinsEta, -4.0, 6.0, forwardBinsPhi, 0., 2*TMath::Pi(),fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge);
    fOutputList->Add(nua_fmd);

    fOutputList->Add(nua_cen);

    dNdeta = new TH1F("dNdeta","dNdeta",100 /*fSettings.fNDiffEtaBins*/,fSettings.fEtaLowEdge,fSettings.fEtaUpEdge);
    fOutputList->Add(dNdeta);
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
    fUtil.fAODevent =  dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fUtil.fAODevent) throw std::runtime_error("Not AOD as expected");
  }
  fSettings.nua_runnumber = fUtil.GetNUARunNumber(fInputEvent->GetRunNumber());

  fUtil.fevent = fInputEvent;
  fUtil.fSettings = fSettings;
  if (fSettings.mc) fUtil.fMCevent = this->MCEvent();

  Int_t   centralEtaBins = (fSettings.useITS ? 200 : 300);
  Int_t   centralPhiBins = (fSettings.useITS ? 20 : 300);
  Double_t centralEtaMin = (fSettings.useSPD ? -2.5 : fSettings.useITS ? -4 : -1.2);
  Double_t centralEtaMax = (fSettings.useSPD ? 2.5 : fSettings.useITS ? 6 : 1.2);

  TH2D centralDist_tmp = TH2D("c","",centralEtaBins,centralEtaMin,centralEtaMax,centralPhiBins,0,2*TMath::Pi());
  centralDist_tmp.SetDirectory(0);
  TH2D refDist_tmp = TH2D("c","",centralEtaBins,centralEtaMin,centralEtaMax,centralPhiBins,0,2*TMath::Pi());
  refDist_tmp.SetDirectory(0);

  TH2D forwardTrRef  ("ft","",200,-4,6,20,0,TMath::TwoPi());
  TH2D forwardPrim  ("fp","",200,-4,6,200,0,TMath::TwoPi());
  forwardTrRef.SetDirectory(0);
  forwardPrim.SetDirectory(0);
  forwardDist = (fSettings.use_primaries_fwd ? &forwardPrim : &forwardTrRef);
  forwardDist->SetDirectory(0);

  centralDist = &centralDist_tmp;
  centralDist->SetDirectory(0);
  refDist = &refDist_tmp;
  refDist->SetDirectory(0);

  fUtil.dodNdeta = kTRUE;
  fUtil.dNdeta = dNdeta;
  fUtil.FillData(refDist,centralDist,forwardDist);


  if (fUtil.dodNdeta){
    for (Int_t etaBin = 1; etaBin <= centralDist->GetNbinsX(); etaBin++) {
      for (Int_t phiBin = 1; phiBin <= centralDist->GetNbinsX(); phiBin++) {
        dNdeta->Fill(centralDist->GetXaxis()->GetBinCenter(etaBin),centralDist->GetBinContent(etaBin, phiBin));
      }
    }
    for (Int_t etaBin = 1; etaBin <= forwardDist->GetNbinsX(); etaBin++) {
      for (Int_t phiBin = 1; phiBin <= forwardDist->GetNbinsX(); phiBin++) {
        dNdeta->Fill(forwardDist->GetXaxis()->GetBinCenter(etaBin),forwardDist->GetBinContent(etaBin, phiBin));
      }
    }
  }

  Double_t zvertex = fUtil.GetZ();

  static_cast<TH1D*>(fEventList->FindObject("Vertex"))->Fill(zvertex);


    // loop for the TPC
    for (Int_t etaBin = 1; etaBin <= centralDist->GetNbinsX(); etaBin++) {
      for (Int_t phiBin = 1; phiBin <= centralDist->GetNbinsY(); phiBin++) {
        if (centralDist->GetBinContent(etaBin,phiBin) == 0) continue;//TH2D
        nua_cen->Fill(centralDist->GetXaxis()->GetBinCenter(etaBin),centralDist->GetYaxis()->GetBinCenter(phiBin),zvertex,centralDist->GetBinContent(etaBin,phiBin)); //TH3D
      }
    }

    if (fSettings.makeFakeHoles) fUtil.MakeFakeHoles(*forwardDist);
    // loop for the FMD
    for (Int_t etaBin = 1; etaBin <= forwardDist->GetNbinsX(); etaBin++) {
      for (Int_t phiBin = 1; phiBin <= forwardDist->GetNbinsY(); phiBin++) {
          Double_t eta = forwardDist->GetXaxis()->GetBinCenter(etaBin);
          Double_t phi = forwardDist->GetYaxis()->GetBinCenter(phiBin);
//        if (fSettings.mc & fSettings.esd){
//          if (!fUtil.FMDAcceptanceExistMC(eta,phi,zvertex)) continue;
//        }
        Double_t weight = forwardDist->GetBinContent(etaBin,phiBin);

        if (fSettings.nua_mode & fSettings.kInterpolate)
          weight = InterpolateWeight(*forwardDist,phiBin,etaBin,weight);

        if (weight == 0) continue;

        nua_fmd->Fill(eta,phi,zvertex,weight);
      }
    }

  PostData(1, fOutputList);
  return;
}



Double_t AliForwardNUATask::InterpolateWeight(TH2D& forwarddNdedp,Int_t phiBin, Int_t etaBin, Double_t weight)
{

  if ((phiBin == 17) && (etaBin >= 125 && etaBin <= 137)){

    //std::cout << "interpolating 1 " << std::endl;

    if (!(weight == 0 || forwarddNdedp.GetBinContent(etaBin, 18) == 0)) return weight;
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
  //std::cout << "weight = " << weight << std::endl;
    //std::cout << "interpolating 2 " << std::endl;

    if (!(weight == 0 || forwarddNdedp.GetBinContent(etaBin, 17) == 0 ) )return weight;
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
    //std::cout << "interpolating 3 " << std::endl;

    if (!(weight == 0)) return weight;

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





Double_t AliForwardNUATask::InterpolateWeight(TH2D*& forwarddNdedp,Int_t phiBin, Int_t etaBin, Double_t weight)
{
  if ((phiBin == 17) && (etaBin >= 125 && etaBin <= 137)){

    if (weight > 0 || forwarddNdedp->GetBinContent(etaBin, 18) > 0) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin - 1, phiBin) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin - 2, phiBin) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin - 3, phiBin) > 0 ) return weight;

    // std::cout << "found hole, etaBin = " << etaBin << ", eta = " << eta << ", phiBin = " << phiBin << ", phi = " << phi << std::endl;
 //    if (detType == "forward") weight = 1.;
     Double_t up = forwarddNdedp->GetBinContent(etaBin, 19);
     Double_t low = forwarddNdedp->GetBinContent(etaBin, 16);
     weight = ((up+low)/2+low)/2;
     return weight;
    //std::cout << weight << std::endl;

   }
  if ((phiBin == 18) && (etaBin >= 125 && etaBin <= 137)){
    if (weight > 0 || forwarddNdedp->GetBinContent(etaBin, 17) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin - 1, phiBin) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin - 2, phiBin) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin - 3, phiBin) > 0 ) return weight;

    // std::cout << "found hole, etaBin = " << etaBin << ", eta = " << eta << ", phiBin = " << phiBin << ", phi = " << phi << std::endl;
     Double_t up = forwarddNdedp->GetBinContent(etaBin, 19);
     Double_t low = forwarddNdedp->GetBinContent(etaBin, 16);
     weight = ((up+low)/2+up)/2;
    //std::cout << weight << std::endl;
  }
  if ((phiBin == 14) && (etaBin >= 168 && etaBin <= 185)){
    if (weight > 0) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin + 1, phiBin) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin + 2, phiBin) > 0 ) return weight;
    //if (forwarddNdedp.GetBinContent(etaBin + 3, phiBin) > 0 ) return weight;
    // std::cout << "found hole, etaBin = " << etaBin << ", eta = " << eta << ", phiBin = " << phiBin << ", phi = " << phi << std::endl;
    Double_t  up = forwarddNdedp->GetBinContent(etaBin, 15);
    Double_t  low = forwarddNdedp->GetBinContent(etaBin, 13);
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
