/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////
//
// Check basic detector results at ESD level
//   - Geometrical efficiency  
//   - Tracking efficiency  
//   - PID efficiency  
//   - Refit efficiency  
//
// Author
//   Alex Bercuci <A.Bercuci@gsi.de>
//
//////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH2I.h>
#include <TH2F.h>
#include <TH3S.h>
#include <TH3F.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TChain.h>
#include <TParticle.h>
#include <TTimeStamp.h>
#include <TRandom.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisCuts.h"
#include "AliESDEvent.h"
#include "AliESDkink.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDpid.h"

#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliTrackReference.h"
//#include "AliESDCentrality.h"
#include "AliMultiplicity.h"
#include "AliCFContainer.h"

#include "AliTRDcheckESD.h"
#include <iostream>
using namespace std;

ClassImp(AliTRDcheckESD)

const Float_t AliTRDcheckESD::fgkxTPC = 290.;
const Float_t AliTRDcheckESD::fgkxTOF = 365.;
const UChar_t AliTRDcheckESD::fgkNgraph[AliTRDcheckESD::kNrefs] ={
8, 4, 2, 20};
FILE* AliTRDcheckESD::fgFile = NULL;

const Float_t AliTRDcheckESD::fgkEvVertexZ = 15.;
const Int_t   AliTRDcheckESD::fgkEvVertexN = 1;
const Float_t AliTRDcheckESD::fgkTrkDCAxy  = 40.;
const Float_t AliTRDcheckESD::fgkTrkDCAz   = 15.;
const Int_t   AliTRDcheckESD::fgkNclTPC    = 100;
const Float_t AliTRDcheckESD::fgkPt        = 0.2;
const Float_t AliTRDcheckESD::fgkEta       = 0.9;
const Float_t AliTRDcheckESD::fgkQs        = 0.002;

//____________________________________________________________________
AliTRDcheckESD::AliTRDcheckESD():
  AliAnalysisTaskSE()
  ,fStatus(0)
  ,fNRefFigures(0)
  ,fESD(NULL)
  ,fMC(NULL)
  ,fESDpid(new AliESDpid)
  ,fHistos(NULL)
  ,fResults(NULL)
  ,fCfContainer(NULL)
  ,fReferenceTrackFilter(NULL)
{
  //
  // Default constructor
  //
  SetNameTitle("TRDcheckESD", "Check TRD @ ESD level");
  SetMC(kTRUE);
}

//____________________________________________________________________
AliTRDcheckESD::AliTRDcheckESD(char* name):
  AliAnalysisTaskSE(name)
  ,fStatus(0)
  ,fNRefFigures(0)
  ,fESD(NULL)
  ,fMC(NULL)
  ,fESDpid(new AliESDpid)
  ,fHistos(NULL)
  ,fResults(NULL)
  ,fCfContainer(NULL)
  ,fReferenceTrackFilter(NULL)
{
  //
  // Default constructor
  //
  SetMC(kTRUE);
  SetTitle("Check TRD @ ESD level");
  DefineOutput(1, TObjArray::Class());
}

//____________________________________________________________________
AliTRDcheckESD::~AliTRDcheckESD()
{
// Destructor
  if(fHistos && !(AliAnalysisManager::GetAnalysisManager() && AliAnalysisManager::GetAnalysisManager()->IsProofMode())){
    if(fHistos->IsOwner()) fHistos->Delete();
    delete fHistos;
    fHistos = NULL;
  }
  if(fResults){
    fResults->Delete();
    delete fResults;
  }
}

//____________________________________________________________________
void AliTRDcheckESD::UserCreateOutputObjects()
{	
  //
  // Create Output Containers (TObjectArray containing 1D histograms)
  //
  Histos();
  PostData(1, fHistos);
}

//____________________________________________________________________
void AliTRDcheckESD::MakeSummaryFromCF(Double_t* trendValues, Bool_t useIsolatedBC, Bool_t cutTOFbc){
  //
  // Draw summary plots for the ESDcheck task using the CF container
  //

  cout << "Make summary from CF" << endl;
  TCanvas *cOut=0x0;
  if(gROOT->FindObject("trackingSummary")) delete gROOT->FindObject("trackingSummary");
  cOut = new TCanvas("trackingSummary", "Tracking summary for the ESD task", 1600, 1200);
  cOut->cd();
  PlotTrackingSummaryFromCF(0, trendValues, useIsolatedBC, cutTOFbc);
  cOut->SaveAs("trackingSummary.gif");
  
  if(gROOT->FindObject("pidSummary")) delete gROOT->FindObject("pidSummary");
  cOut = new TCanvas("pidSummary", "PID summary for the ESD task", 1600, 1200);
  cOut->cd();
  //GetRefFigure(6);
  PlotPidSummaryFromCF(0, trendValues, useIsolatedBC, cutTOFbc);
  cOut->SaveAs("pidSummary.gif");

  if(gROOT->FindObject("centSummary")) delete gROOT->FindObject("centSummary");
  cOut = new TCanvas("centSummary", "Centrality summary for the ESD task", 1600, 1200);
  cOut->cd();
  //GetRefFigure(7);
  PlotCentSummaryFromCF(trendValues, useIsolatedBC, cutTOFbc);
  cOut->SaveAs("centSummary.gif");
    
}


//____________________________________________________________________
void AliTRDcheckESD::MakeSummary(Double_t* trendValues){
  //
  // Draw summary plots for the ESDcheck task
  //
  // Old method to draw summary pictures from histograms. Use the MakeSummaryFromCF() when CF container is present
  
  cout << "Make summary" << endl;
  TCanvas *cTracking=0x0;
  if(gROOT->FindObject("trackingSummary")) delete gROOT->FindObject("trackingSummary");
  cTracking = new TCanvas("trackingSummary", "Tracking summary for the ESD task", 1600, 1200);
  cTracking->cd();
  //GetRefFigure(5);
  if(PlotTrackingSummary(0, trendValues))
    cTracking->SaveAs("trackingSummary.gif");
  
  TCanvas* cPid = 0x0;
  if(gROOT->FindObject("pidSummary")) delete gROOT->FindObject("pidSummary");
  cPid = new TCanvas("pidSummary", "PID summary for the ESD task", 1600, 1200);
  cPid->cd();
  //GetRefFigure(6);
  if(PlotPidSummary(0, trendValues))
    cPid->SaveAs("pidSummary.gif");

  TCanvas* cCent=0x0;
  if(gROOT->FindObject("centSummary")) delete gROOT->FindObject("centSummary");
  cCent = new TCanvas("centSummary", "Centrality summary for the ESD task", 1600, 1200);
  cCent->cd();
  //GetRefFigure(7);
  if(PlotCentSummary(trendValues))
    cCent->SaveAs("centSummary.gif");
}

//____________________________________________________________________
Bool_t AliTRDcheckESD::GetRefFigure(Int_t ifig)
{
  //
  // Produce reference Plots during PostProcessing
  //
  if(ifig>=fNRefFigures){
    AliWarning(Form("Ref plot %d not available. Valid only up to %d", ifig, fNRefFigures));
    return kFALSE;
  }
  if(!gPad){
    AliWarning("Please provide a canvas to draw results.");
    return kFALSE;
  } else {
    gPad->SetLogx(0);gPad->SetLogy(0);
    gPad->SetMargin(0.125, 0.015, 0.1, 0.015);
  }

  const Char_t *title[20];
  TH1 *hF(NULL);
  if((hF=(TH1S*)gROOT->FindObject("hFcheckESD"))) delete hF;
  TLegend *leg(NULL);
  TList *l(NULL); TVirtualPad *pad(NULL);
  TGraphErrors *g(NULL);TGraphAsymmErrors *ga(NULL);
  TObjArray *arr(NULL);
  TLatex *lat=new TLatex();
  lat->SetTextSize(0.07);
  lat->SetTextColor(2);
  TLine line;
  TTimeStamp now;
  switch(ifig){
  case kNCl: // number of clusters/track
    if(!(arr = (TObjArray*)fResults->At(kNCl))) return kFALSE;

    leg = new TLegend(.83, .7, .99, .96);
    leg->SetHeader("Species");
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    for(Int_t ig(0); ig<fgkNgraph[kNCl-1]; ig++){
      if(!(g = (TGraphErrors*)arr->At(ig))) return kFALSE;
      if(!g->GetN()) continue;
      g->Draw(ig?"pc":"apc"); leg->AddEntry(g, g->GetTitle(), "pl");
      if(ig) continue;
      hF=g->GetHistogram();
      hF->SetXTitle("no of clusters");
      hF->SetYTitle("entries"); 
      hF->GetYaxis()->CenterTitle(1);
      hF->GetYaxis()->SetTitleOffset(1.2);
      hF->SetMinimum(5);
    }
    leg->Draw(); gPad->SetLogy();
    break;
  case kTRDstat: // Efficiency
    if(!(arr = (TObjArray*)fResults->At(kTRDstat))) return kFALSE;
    leg = new TLegend(.62, .77, .98, .98);
    leg->SetHeader("TRD Efficiency");
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    title[0] = "Geometrical (TRDin/TPCout)";
    title[1] = "Tracking (TRDout/TRDin)";
    title[2] = "PID (TRDpid/TRDin)";
    title[3] = "Refit (TRDrefit/TRDin)";
    hF = new TH1S("hFcheckESD", ";p [GeV/c];Efficiency", 10, 0.1, 10.);
    hF->SetMaximum(1.4);
    hF->GetXaxis()->SetMoreLogLabels();
    hF->GetYaxis()->CenterTitle(1);
    hF->Draw("p");
    for(Int_t ig(0); ig<fgkNgraph[kTRDstat-1]; ig++){
      if(!(g = (TGraphErrors*)arr->At(ig))) return kFALSE;
      g->Draw("pl"); leg->AddEntry(g, title[ig], "pl");
      //PutTrendValue(name[id], g->GetMean(2));
      //PutTrendValue(Form("%sRMS", name[id]), g->GetRMS(2));
    }
    leg->Draw(); gPad->SetLogx();
    break;
  case kTRDmom: // Energy loss
    if(!(arr = (TObjArray*)fResults->At(kTRDmom))) return kFALSE;
    leg = new TLegend(.65, .7, .95, .99);
    leg->SetHeader("Energy Loss");
    leg->SetBorderSize(1); leg->SetFillColor(0);
    title[0] = "Max & 90% quantile";
    title[1] = "Mean & 60% quantile";
    hF = new TH1S("hFcheckESD", ";layer;#Delta E", 6, -0.5, 5.5);
    hF->SetMaximum(1.3);hF->SetMinimum(-.3);
    hF->Draw("p");
    for(Int_t ig(0); ig<fgkNgraph[kTRDmom-1]; ig++){
      if(!(ga = (TGraphAsymmErrors*)arr->At(ig))) return kFALSE;
      ga->Draw("pl"); leg->AddEntry(ga, title[ig], "pl");
      //PutTrendValue(name[id], g->GetMean(2));
      //PutTrendValue(Form("%sRMS", name[id]), g->GetRMS(2));
    }
    leg->Draw();gPad->SetLogx(kFALSE);
    break;
  case kPtRes: // Pt resolution @ vertex
    if(!(arr = (TObjArray*)fResults->At(kPtRes))) return kFALSE;
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = ((TVirtualPad*)l->At(0)); pad->cd(); pad->SetLogx();
    pad->SetMargin(0.1, 0.022, 0.1, 0.023);
    hF = new TH1S("hFcheckESD", "ITS+TPC+TRD;p_{t} [GeV/c];#Delta p_{t} / p_{t} [%]", 10, 0.2, 10.);
    hF->SetMaximum(10.);hF->SetMinimum(-3.);
    hF->GetXaxis()->SetMoreLogLabels();
    hF->GetXaxis()->SetTitleOffset(1.2);
    hF->GetYaxis()->CenterTitle();
    hF->Draw("p");
    //for(Int_t ig(0); ig<fgkNgraph[kPtRes-1]/2; ig++){
    for(Int_t ig(2); ig<6; ig++){
      if(!(g = (TGraphErrors*)arr->At(ig))) continue;
      if(!g->GetN()) continue;
      g->Draw("pl");
      //PutTrendValue(name[id], g->GetMean(2));
      //PutTrendValue(Form("%sRMS", name[id]), g->GetRMS(2));
    }
    pad = ((TVirtualPad*)l->At(1)); pad->cd(); pad->SetLogx();
    pad->SetMargin(0.1, 0.22, 0.1, 0.023);
    hF = (TH1*)hF->Clone("hFcheckESD1");
    hF->SetTitle("ITS+TPC");
    hF->SetMaximum(10.);hF->SetMinimum(-3.);
    hF->Draw("p");
    leg = new TLegend(.78, .1, .99, .98);
    leg->SetHeader("P_{t} @ DCA");
    leg->SetBorderSize(1); leg->SetFillColor(0);
    leg->SetTextAlign(22);
    leg->SetTextFont(12);
    leg->SetTextSize(0.03813559);
    {
      Int_t nPlots(0);
      //for(Int_t ig(fgkNgraph[kPtRes-1]/2); ig<fgkNgraph[kPtRes-1]; ig++){
      for(Int_t ig(12); ig<16; ig++){
        if(!(g = (TGraphErrors*)arr->At(ig))) continue;
        if(!g->GetN()) continue;
        nPlots++;
        g->Draw("pl"); leg->AddEntry(g, g->GetTitle(), "pl");
        //PutTrendValue(name[id], g->GetMean(2));
        //PutTrendValue(Form("%sRMS", name[id]), g->GetRMS(2));
      }
      if(nPlots) leg->Draw();
    }
    break;
  case 5: // plot a 3x3 canvas with tracking related histograms
    PlotTrackingSummary(0);
    break;
    
  case 6: // plot a 3x3 canvas with PID related histograms
    PlotPidSummary(0);
    break;

  case 7: // plot a 3x3 canvas with centrality dependence histograms
    PlotCentSummary();
    break;

  }
  return kTRUE;
}

//____________________________________________________________________
void AliTRDcheckESD::UserExec(Option_t *){
  //
  // Run the Analysis
  //
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  fMC = MCEvent();

  if(!fESD){
    AliError("ESD event missing.");
    return;
  }
  
  // Get MC information if available
  AliStack * fStack = NULL;
  if(HasMC()){
    if(!fMC){ 
      AliWarning("MC event missing");
      SetMC(kFALSE);
    } else {
      if(!(fStack = fMC->Stack())){
        AliWarning("MC stack missing");
        SetMC(kFALSE);
      }
    }
  }
  TH1 *h(NULL);
  
  Double_t values[kNTrdCfVariables];      // array where the CF container variables are stored
  values[kEventVtxZ] = fESD->GetPrimaryVertex()->GetZv();
  values[kEventBC] = fESD->GetBunchCrossNumber();
  
  const AliMultiplicity* mult=fESD->GetMultiplicity();
  Double_t itsNTracklets = mult->GetNumberOfTracklets();
  if(itsNTracklets<1) return;
  Int_t multLimits[6] = {0, 700, 1400, 2100, 2800, 3500};
  Int_t centralityClass = 0;
  for(Int_t iCent=0; iCent<5; ++iCent) {
    if(itsNTracklets>=multLimits[iCent] && itsNTracklets<multLimits[iCent+1])
      centralityClass=iCent+1;
  }
  values[kEventMult] = itsNTracklets;
  if(centralityClass == 0) return;
      
  AliESDtrack *esdTrack(NULL);
  for(Int_t itrk = 0; itrk < fESD->GetNumberOfTracks(); itrk++){
    esdTrack = fESD->GetTrack(itrk);
    if(!fReferenceTrackFilter->IsSelected(esdTrack)) continue;
    
    ULong_t status = esdTrack->GetStatus(); //PrintStatus(status);
            
    // pid quality
    Bool_t kBarrel = Bool_t(status & AliESDtrack::kTRDin);

    // find position and momentum of the track at entrance in TRD
    Double_t localCoord[3] = {0., 0., 0.};
    Bool_t localCoordGood = esdTrack->GetXYZAt(298., fESD->GetMagneticField(), localCoord);
    Double_t localMom[3] = {0., 0., 0.};
    Bool_t localMomGood = esdTrack->GetPxPyPzAt(298., fESD->GetMagneticField(), localMom);
    //Double_t localPhi = (localMomGood ? TMath::ATan2(localMom[1], localMom[0]) : 0.0);
    Double_t localSagitaPhi = (localCoordGood ? TMath::ATan2(localCoord[1], localCoord[0]) : 0.0);

    values[kTrackTOFdeltaBC] = esdTrack->GetTOFDeltaBC();
    values[kTrackTOFBC]      = esdTrack->GetTOFBunchCrossing(fESD->GetMagneticField());
    Float_t dcaXY=0.0; Float_t dcaZ=0.0;
    esdTrack->GetImpactParameters(dcaXY, dcaZ);
    values[kTrackDCAxy]  = dcaXY;
    values[kTrackDCAz]   = dcaZ;
    values[kTrackCharge] = esdTrack->Charge();
    values[kTrackPhi]    = localSagitaPhi;
    values[kTrackEta]    = esdTrack->Eta();
    values[kTrackPt]     = esdTrack->Pt();
    values[kTrackP]      = esdTrack->P();
    values[kTrackTrdTracklets] = esdTrack->GetTRDntracklets();
    values[kTrackTrdClusters] = esdTrack->GetTRDncls();
    for(Int_t i=0; i<6; ++i) values[kTrackQtot+i] = 0.0;
        
    if(localCoordGood && localMomGood) fCfContainer->Fill(values, 0);      // fill the TPC reference step
            
    // TRD reference tracks
    if(esdTrack->GetTRDntracklets()>=1) {
      // (slicePH,sliceNo) distribution and Qtot from slices
      for(Int_t iPlane=0; iPlane<6; iPlane++) {
        Float_t qtot=esdTrack->GetTRDslice(iPlane, 0);
        for(Int_t iSlice=0; iSlice<8; iSlice++) {
	  if(esdTrack->GetTRDslice(iPlane, iSlice)>20.) {
	    h = (TH2F*)fHistos->At(kPHSlice); h->Fill(iSlice, esdTrack->GetTRDslice(iPlane, iSlice));
	    h = (TH2F*)fHistos->At(kPHSlice+centralityClass); h->Fill(iSlice, esdTrack->GetTRDslice(iPlane, iSlice));
	  }
	}
	values[kTrackQtot+iPlane] = fgkQs*qtot;
      }
            
      if(localCoordGood && localMomGood) {
        fCfContainer->Fill(values, 1);
        if(Bool_t(status & AliESDtrack::kTOFpid)) fCfContainer->Fill(values, 2);
      }
    }  // end if nTRDtrkl>=1
    
    // look at external track param
    const AliExternalTrackParam *op = esdTrack->GetOuterParam();
    const AliExternalTrackParam *ip = esdTrack->GetInnerParam();

    Double_t pt(0.), pt0(0.), ptTRD(0.); 
    // read MC info if available
    Bool_t kFOUND(kFALSE), kPhysPrim(kFALSE);
    AliMCParticle *mcParticle(NULL);
    if(HasMC()){
      AliTrackReference *ref(NULL); 
      Int_t fLabel(esdTrack->GetLabel());
      Int_t fIdx(TMath::Abs(fLabel));
      if(!fStack || fIdx > fStack->GetNtrack()) continue; 
      
      // read MC particle 
      if(!(mcParticle = (AliMCParticle*) fMC->GetTrack(fIdx))) {
        AliWarning(Form("MC particle missing. Label[ %d].", fLabel));
        continue;
      }
  
      pt   = esdTrack->Pt();
      pt0  = mcParticle->Pt();
      //Double_t eta0 = mcParticle->Eta();
      //Double_t phi0 = mcParticle->Phi();
      kPhysPrim = fMC->IsPhysicalPrimary(fIdx);

      // read track references
      Int_t nRefs = mcParticle->GetNumberOfTrackReferences();
      if(!nRefs){
        AliWarning(Form("No TR found for track @ Label[%d].", fLabel));
        continue;
      }
      Int_t iref = 0;
      while(iref<nRefs){
        ref = mcParticle->GetTrackReference(iref);
        if(ref->LocalX() > fgkxTPC) break;
        ref=NULL; iref++;
      }
      if(ref){ 
        if(ref->LocalX() > fgkxTOF){ // track skipping TRD fiducial volume
          ref = mcParticle->GetTrackReference(TMath::Max(iref-1, 0));
        }
      } else { // track stopped in TPC 
        ref = mcParticle->GetTrackReference(TMath::Max(iref-1, 0));
      }
      ptTRD = ref->Pt();kFOUND=kTRUE;
    } else { // use reconstructed values
      if(op){
        Double_t x(op->GetX());
        if(x<fgkxTOF && x>fgkxTPC){
          ptTRD=op->Pt();
          kFOUND=kTRUE;
        }
      }

      if(!kFOUND && ip){
        ptTRD=ip->Pt();
        kFOUND=kTRUE;
      }
    }     // end if(HasMC())

    if(kFOUND){
      h = (TH2I*)fHistos->At(kTRDstat);
      if(status & AliESDtrack::kTPCout) h->Fill(ptTRD, kTPCout);
      if(status & AliESDtrack::kTRDin) h->Fill(ptTRD, kTRDin);
      if(kBarrel && (status & AliESDtrack::kTRDout)) h->Fill(ptTRD, kTRDout);
      if(kBarrel && (status & AliESDtrack::kTRDpid)) h->Fill(ptTRD, kTRDpid);
      if(kBarrel && (status & AliESDtrack::kTRDrefit)) h->Fill(ptTRD, kTRDref);
    }
    Int_t idx(HasMC() ? Pdg2Idx(TMath::Abs(mcParticle->PdgCode())): 0)
         ,sgn(esdTrack->Charge()<0?0:1);
    if(kBarrel && kPhysPrim) {
      TH3 *h3 = (TH3S*)fHistos->At(kPtRes);
      Int_t offset = (status & AliESDtrack::kTRDrefit) ? 0 : 10; 
      h3->Fill(pt0, 1.e2*(pt/pt0-1.), 
        offset + 2*idx + sgn);
    }
    ((TH1*)fHistos->At(kNCl))->Fill(esdTrack->GetTRDncls(), 2*idx + sgn);
    if(ip){
      h = (TH2I*)fHistos->At(kTRDmom);
      Float_t pTRD(0.);
      for(Int_t ily=6; ily--;){
        if((pTRD=esdTrack->GetTRDmomentum(ily))<0.) continue;
        h->Fill(ip->GetP()-pTRD, ily);
      }
    }
  }  // end loop over tracks
  
  // fill the number of tracks histograms
  //h = (TH1I*)fHistos->At(kNTracksAcc);
  //h->Fill(nTracksAcc);
  //h = (TH1I*)fHistos->At(kNTracksTPC);
  //h->Fill(nTracksTPC);
  PostData(1, fHistos);
}

//____________________________________________________________________
TObjArray* AliTRDcheckESD::Histos()
{
// Retrieve histograms array if already build or build it

  if(fHistos) return fHistos;

  fHistos = new TObjArray(kNhistos+1);
  fHistos->SetOwner(kTRUE);

  TH1 *h = NULL;

  // clusters per track
  const Int_t kNpt(30);
  Float_t pt(0.2);
  Float_t binsPt[kNpt+1];
  for(Int_t i=0;i<kNpt+1; i++,pt+=(TMath::Exp(i*i*.001)-1.)) binsPt[i]=pt;
  if(!(h = (TH2I*)gROOT->FindObject("hNCl"))){
    h = new TH2I("hNCl", "Clusters per TRD track;N_{cl}^{TRD};SPECIES;entries", 60, 0., 180., 10, -0.5, 9.5);
    TAxis *ay(h->GetYaxis());
    ay->SetLabelOffset(0.015);
    for(Int_t i(0); i<AliPID::kSPECIES; i++){
      ay->SetBinLabel(2*i+1, Form("%s^{-}", AliPID::ParticleLatexName(i)));
      ay->SetBinLabel(2*i+2, Form("%s^{+}", AliPID::ParticleLatexName(i)));
    }
  } else h->Reset();
  fHistos->AddAt(h, kNCl); fNRefFigures++;

  // status bits histogram
  const Int_t kNbits(5);
  Float_t bits(.5);
  Float_t binsBits[kNbits+1];
  for(Int_t i=0; i<kNbits+1; i++,bits+=1.) binsBits[i]=bits;
  if(!(h = (TH2I*)gROOT->FindObject("hTRDstat"))){
    h = new TH2I("hTRDstat", "TRD status bits;p_{t} @ TRD [GeV/c];status;entries", kNpt, binsPt, kNbits, binsBits);
    TAxis *ay(h->GetYaxis());
    ay->SetBinLabel(1, "kTPCout");
    ay->SetBinLabel(2, "kTRDin");
    ay->SetBinLabel(3, "kTRDout");
    ay->SetBinLabel(4, "kTRDpid");
    ay->SetBinLabel(5, "kTRDrefit");
  } else h->Reset();
  fHistos->AddAt(h, kTRDstat);

  // energy loss
  if(!(h = (TH2I*)gROOT->FindObject("hTRDmom"))){
    h = new TH2I("hTRDmom", "TRD energy loss;p_{inner} - p_{ly} [GeV/c];ly;entries", 100, -1., 2., 6, -0.5, 5.5);
  } else h->Reset();
  fHistos->AddAt(h, kTRDmom);
  //if(!HasMC()) return fHistos;

  // pt resolution
  const Int_t kNdpt(100), kNspec(4*AliPID::kSPECIES);
  Float_t dpt(-3.), spec(-0.5);
  Float_t binsDPt[kNdpt+1], binsSpec[kNspec+1];
  for(Int_t i=0; i<kNdpt+1; i++,dpt+=6.e-2) binsDPt[i]=dpt;
  for(Int_t i=0; i<kNspec+1; i++,spec+=1.) binsSpec[i]=spec;
  if(!(h = (TH3S*)gROOT->FindObject("hPtRes"))){
    h = new TH3S("hPtRes", "P_{t} resolution @ DCA;p_{t}^{MC} [GeV/c];#Delta p_{t}/p_{t}^{MC} [%];SPECIES", kNpt, binsPt, kNdpt, binsDPt, kNspec, binsSpec);
    TAxis *az(h->GetZaxis());
    az->SetLabelOffset(0.015);
    for(Int_t i(0); i<AliPID::kSPECIES; i++){
      az->SetBinLabel(2*i+1, Form("%s^{-}", AliPID::ParticleLatexName(i)));
      az->SetBinLabel(2*i+2, Form("%s^{+}", AliPID::ParticleLatexName(i)));
      az->SetBinLabel(10+2*i+1, Form("%s^{-}", AliPID::ParticleLatexName(i)));
      az->SetBinLabel(10+2*i+2, Form("%s^{+}", AliPID::ParticleLatexName(i)));
    }
  } else h->Reset();
  fHistos->AddAt(h, kPtRes);

  // TPC event vertex distribution
  if(!(h = (TH1F*)gROOT->FindObject("hTPCVertex"))){
    h = new TH1F("hTPCVertex", "Event vertex Z coord. from TPC tracks", 100, -25., 25.);
  } else h->Reset();
  fHistos->AddAt(h, kTPCVertex);
  
  // Event vertex
  if(!(h = (TH1F*)gROOT->FindObject("hEventVertex"))){
    h = new TH1F("hEventVertex", "Event vertex Z coord.", 100, -25., 25.);
  } else h->Reset();
  fHistos->AddAt(h, kEventVertex);
  
  // Number of all tracks
  if(!(h = (TH1I*)gROOT->FindObject("hNTracksAll"))){
    h = new TH1I("hNTracksAll", "Number of tracks per event, event vertex cuts", 5000, 0, 5000);
  } else h->Reset();
  fHistos->AddAt(h, kNTracksAll);
  
  // Number of tracks in acceptance and DCA cut
  if(!(h = (TH1I*)gROOT->FindObject("hNTracksAcc"))){
    h = new TH1I("hNTracksAcc", Form("Number of tracks per event, |#eta|<%.1f, |DCAxy|<%.1f, |DCAz|<%.1f",
				     fgkEta, fgkTrkDCAxy, fgkTrkDCAz), 5000, 0, 5000);
  } else h->Reset();
  fHistos->AddAt(h, kNTracksAcc);
  
  // Number of tracks in TPC (Ncls>10)
  if(!(h = (TH1I*)gROOT->FindObject("hNTracksTPC"))){
    h = new TH1I("hNTracksTPC", Form("Number of tracks per event, |#eta|<%.1f, pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
				     fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 5000, 0, 5000);
  } else h->Reset();
  fHistos->AddAt(h, kNTracksTPC);
  
  // Distribution of DCA-xy
  if(!(h = (TH1F*)gROOT->FindObject("hDCAxy"))){
    h = new TH1F("hDCAxy", "Distribution of transverse DCA", 100, -100., 100.);
  } else h->Reset();
  fHistos->AddAt(h, kDCAxy);
  
  // Distribution of DCA-z
  if(!(h = (TH1F*)gROOT->FindObject("hDCAz"))){
    h = new TH1F("hDCAz", "Distribution of longitudinal DCA", 100, -100., 100.);
  } else h->Reset();
  fHistos->AddAt(h, kDCAz);
  
  Double_t binPtLimits[33] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
		 	     1.0, 1.1, 1.2, 1.3, 1.4, 
			     1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0,
			     3.4, 3.8, 4.2, 4.6, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  // Pt distributions
  if(!(h = (TH1F*)gROOT->FindObject("hPt1"))){
    h = new TH1F("hPt1", Form("dN/dpt, |#eta|<%.1f and pt>%.1f", fgkEta, fgkPt), 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kPt1);
  
  if(!(h = (TH1F*)gROOT->FindObject("hPt2"))){
    h = new TH1F("hPt2", Form("dN/dpt, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f",
			      fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz), 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kPt2);
  
  if(!(h = (TH1F*)gROOT->FindObject("hPt3pos"))){
    h = new TH1F("hPt3pos", Form("dN/dpt (positives), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kPt3pos);
  
  if(!(h = (TH1F*)gROOT->FindObject("hPt3neg"))){
    h = new TH1F("hPt3neg", Form("dN/dpt (negatives), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kPt3neg);
  
  if(!(h = (TH1F*)gROOT->FindObject("hPt4pos"))){
    h = new TH1F("hPt4pos", Form("dN/dpt (positives), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq 1",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kPt4pos);
  
  if(!(h = (TH1F*)gROOT->FindObject("hPt4neg"))){
    h = new TH1F("hPt4pos", Form("dN/dpt (negatives), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq 1",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kPt4neg);
  
  // theta distribution of TRD tracks
  if(!(h = (TH1F*)gROOT->FindObject("hTheta"))){
    h = new TH1F("hTheta", Form("dN/d#theta, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq 1",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 220,.5,2.7);
  } else h->Reset();
  fHistos->AddAt(h, kTheta);
  
  // phi distribution of TRD tracks
  if(!(h = (TH1F*)gROOT->FindObject("hPhi"))){
    h = new TH1F("hPhi", Form("dN/d#varphi, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq 1",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 157,0,6.28);
  } else h->Reset();
  fHistos->AddAt(h, kPhi);
  
  // TPC cluster distribution
  if(!(h = (TH1F*)gROOT->FindObject("hNTPCCl"))){
    h = new TH1I("hNTPCCl", Form("Number of TPC clusters/track, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz), 160, 0, 160);
  } else h->Reset();
  fHistos->AddAt(h, kNTPCCl);
  
  if(!(h = (TH1I*)gROOT->FindObject("hNTPCCl2"))){
    h = new TH1F("hNTPCCl2", Form("Number of TPC clusters/track, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, pt>1.0 GeV/c",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz), 160, 0, 160);
  } else h->Reset();
  fHistos->AddAt(h, kNTPCCl2);
  
  // dE/dx vs P for TPC reference tracks
  if(!(h = (TH2F*)gROOT->FindObject("hTPCDedx"))){
    h = new TH2F("hTPCDedx", Form("TPC dE/dx vs P, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, 0.1,10.1, 150, 0, 150.);
  } else h->Reset();
  fHistos->AddAt(h, kTPCDedx);
  
  // eta,phi distribution of TPC reference tracks
  if(!(h = (TH2F*)gROOT->FindObject("hEtaPhi"))){
    h = new TH2F("hEtaPhi", Form("TPC (#eta,#varphi), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 50, -1, 1, 157, 0, 6.28);
  } else h->Reset();
  fHistos->AddAt(h, kEtaPhi);
  
  // Nclusters vs eta distribution for TPC tracks
  if(!(h = (TH2F*)gROOT->FindObject("hEtaNclsTPC"))){
    h = new TH2F("hEtaNclsTPC", Form("TPC Nclusters vs. #eta, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz), 50, -1, 1, 160, 0, 160.);
  } else h->Reset();
  fHistos->AddAt(h, kEtaNclsTPC);
  
  // Nclusters vs phi distribution for TPC reference tracks
  if(!(h = (TH2F*)gROOT->FindObject("hPhiNclsTPC"))){
    h = new TH2F("hPhiNclsTPC", Form("TPC Nclusters vs. #varphi, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz), 157, 0, 6.28, 160, 0, 160.);
  } else h->Reset();
  fHistos->AddAt(h, kPhiNclsTPC);

  // SPD multiplicity distribution
  if(!(h = (TH1F*)gROOT->FindObject("hSPDMult"))){
    h = new TH1F("hSPDMult", "SPD multiplicity", 10000, -0.5, 9999.5);
  } else h->Reset();
  fHistos->AddAt(h, kSPDMult);

  // Ntracklets/track vs P for TRD reference tracks
  Double_t binsP[19] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0,
			2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 9.0, 12.0};
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH2F*)gROOT->FindObject(Form("hNTrackletsTRD_cent%d",iCent+1)))){
      h = new TH2F(Form("hNTrackletsTRD_cent%d",iCent+1), Form("TRD Ntracklets/track vs. P, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
							       iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 18, binsP, 7, -0.5, 6.5);
    } else h->Reset();
    fHistos->AddAt(h, kNTrackletsTRD+iCent);
  }
  
  // Nclusters/track vs P for TRD reference tracks
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH2F*)gROOT->FindObject(Form("hNClsTrackTRD_cent%d",iCent+1)))){
      h = new TH2F(Form("hNClsTrackTRD_cent%d",iCent+1), Form("TRD Nclusters/track vs. P, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
							      iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 18, binsP, 300, 0., 300.);
    } else h->Reset();
    fHistos->AddAt(h, kNClsTrackTRD+iCent);
  }  

  // <PH> vs slice number for TRD reference tracklets
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH2F*)gROOT->FindObject(Form("hPHSlice_cent%d",iCent+1)))){
      h = new TH2F(Form("hPHSlice_cent%d",iCent+1), Form("<PH> vs sliceNo, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
				    iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 8, -0.5, 7.5, 200, 0., 2000.);
    } else h->Reset();
    fHistos->AddAt(h, kPHSlice+iCent);
  }

  // <PH> vs slice number for TRD reference tracklets, from TPC pions
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH2F*)gROOT->FindObject(Form("hPHSliceTPCpions_cent%d",iCent+1)))){
      h = new TH2F(Form("hPHSliceTPCpions_cent%d",iCent+1), Form("<PH> vs sliceNo, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, TPC pions",
								 iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 8, -0.5, 7.5, 100, 0., 2000.);
    } else h->Reset();
    fHistos->AddAt(h, kPHSliceTPCpions+iCent);
  }

  // TPC dE/dx vs P for TRD reference tracks, pions
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH2F*)gROOT->FindObject(Form("hTPCdedxPions_cent%d",iCent+1)))){
      h = new TH2F(Form("hTPCdedxPions_cent%d",iCent+1), Form("TPC dE/dx vs P, TPC pions, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
					 iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, 0.1,10.1, 100, 0,100.);
    } else h->Reset();
    fHistos->AddAt(h, kTPCdedxPions+iCent);
  }

  // <PH> vs slice number for TRD reference tracklets, from TPC electrons
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH2F*)gROOT->FindObject(Form("hPHSliceTPCelectrons_cent%d",iCent+1)))){
      h = new TH2F(Form("hPHSliceTPCelectrons_cent%d",iCent+1), Form("<PH> vs sliceNo, centrality %d,|#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, TPC electrons",
						iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 8, -0.5, 7.5, 100, 0., 2000.);
    } else h->Reset();
    fHistos->AddAt(h, kPHSliceTPCelectrons+iCent);
  }

  // TPC dE/dx vs P for TRD reference tracks, electrons
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH2F*)gROOT->FindObject(Form("hTPCdedxElectrons_cent%d",iCent+1)))){
      h = new TH2F(Form("hTPCdedxElectrons_cent%d",iCent+1), Form("TPC dE/dx vs P, TPC electrons, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
					     iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, 0.1,10.1, 100, 0,100.);
    } else h->Reset();
    fHistos->AddAt(h, kTPCdedxElectrons+iCent);
  }

  // Qtot vs P for TRD reference tracklets
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH2F*)gROOT->FindObject(Form("hQtotP_cent%d",iCent+1)))){
      h = new TH2F(Form("hQtotP_cent%d",iCent+1), Form("Qtot(from slices) vs P, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
				  iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 18, binsP, 400, 0., 20);
    } else h->Reset();
    fHistos->AddAt(h, kQtotP+iCent);
  }
  
  // (X,Y,momentum) distribution after AliESDtrack::PropagateTo(r=300.)
  if(!(h = (TH3F*)gROOT->FindObject("hPropagXYvsP"))){
    h = new TH3F("hPropagXYvsP", Form("(x,y) vs P after AliESDtrack::PropagateTo(r=300.), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100,-500,500, 100,-500,500, 10, 0.,10.);
  } else h->Reset();
  fHistos->AddAt(h, kPropagXYvsP);
  
  // (R,Z,momentum) distribution after AliESDtrack::PropagateTo(r=300.)
  if(!(h = (TH3F*)gROOT->FindObject("hPropagRZvsP"))){
    h = new TH3F("hPropagRZvsP", Form("(r,z) vs P after AliESDtrack::PropagateTo(r=300.), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100,-350., 350., 100,0.,500., 10, 0.,10.);
  } else h->Reset();
  fHistos->AddAt(h, kPropagRZvsP);
  
  Double_t etaBinLimits[101];	
  for(Int_t i=0; i<101; i++) etaBinLimits[i] = -1.0 + i*2.0/100.;
  Double_t phiBinLimits[151];
  for(Int_t i=0; i<151; i++) phiBinLimits[i] = -1.1*TMath::Pi() + i*2.2*TMath::Pi()/150.;
  // (eta,detector phi,P) distribution of reference TPC positive tracks
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH3F*)gROOT->FindObject(Form("hTPCRefTracksPos_cent%d",iCent+1)))){
      h = new TH3F(Form("hTPCRefTracksPos_cent%d",iCent+1), Form("(#eta,detector #varphi,p) for TPC positive reference tracks, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
					    iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
    } else h->Reset();
    fHistos->AddAt(h, kTPCRefTracksPos+iCent);
  }
  
  // (eta,detector phi,P) distribution of reference TPC negative tracks
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH3F*)gROOT->FindObject(Form("hTPCRefTracksNeg_cent%d",iCent+1)))){
      h = new TH3F(Form("hTPCRefTracksNeg_cent%d",iCent+1), Form("(#eta,detector #varphi,p) for TPC negative reference tracks, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
					    iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
    } else h->Reset();
    fHistos->AddAt(h, kTPCRefTracksNeg+iCent);
  }
  
  // (eta,detector phi,P) distribution of reference TRD positive tracks
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH3F*)gROOT->FindObject(Form("hTRDRefTracksPos_cent%d",iCent+1)))){
      h = new TH3F(Form("hTRDRefTracksPos_cent%d",iCent+1), Form("(#eta,detector #varphi,p) for TRD positive reference tracks, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
					    iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
    } else h->Reset();
    fHistos->AddAt(h, kTRDRefTracksPos+iCent);
  }
  
  // (eta,detector phi,P) distribution of reference TRD negative tracks
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH3F*)gROOT->FindObject(Form("hTRDRefTracksNeg_cent%d",iCent+1)))){
      h = new TH3F(Form("hTRDRefTracksNeg_cent%d",iCent+1), Form("(#eta,detector #varphi,p) for TRD negative reference tracks, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
					    iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
    } else h->Reset();
    fHistos->AddAt(h, kTRDRefTracksNeg+iCent);
  }

  // (eta,detector phi,P) distribution of reference TRD positive tracks with 4 tracklets
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH3F*)gROOT->FindObject(Form("hTRDRefTracksPos4_cent%d",iCent+1)))){
      h = new TH3F(Form("hTRDRefTracksPos4_cent%d",iCent+1), Form("(#eta,detector #varphi,p) for TRD positive reference tracks, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
								 iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
    } else h->Reset();
    fHistos->AddAt(h, kTRDRefTracksPos4+iCent);
  }

  // (eta,detector phi,P) distribution of reference TRD positive tracks with 5 tracklets
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH3F*)gROOT->FindObject(Form("hTRDRefTracksPos5_cent%d",iCent+1)))){
      h = new TH3F(Form("hTRDRefTracksPos5_cent%d",iCent+1), Form("(#eta,detector #varphi,p) for TRD positive reference tracks, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
								  iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
    } else h->Reset();
    fHistos->AddAt(h, kTRDRefTracksPos5+iCent);
  }

  // (eta,detector phi,P) distribution of reference TRD positive tracks with 6 tracklets
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH3F*)gROOT->FindObject(Form("hTRDRefTracksPos6_cent%d",iCent+1)))){
      h = new TH3F(Form("hTRDRefTracksPos6_cent%d",iCent+1), Form("(#eta,detector #varphi,p) for TRD positive reference tracks, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
                                                                  iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
    } else h->Reset();
    fHistos->AddAt(h, kTRDRefTracksPos6+iCent);
  }
  
  // (eta,detector phi,P) distribution of reference TRD negative tracks with 4 tracklets
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH3F*)gROOT->FindObject(Form("hTRDRefTracksNeg4_cent%d",iCent+1)))){
      h = new TH3F(Form("hTRDRefTracksNeg4_cent%d",iCent+1), Form("(#eta,detector #varphi,p) for TRD negative reference tracks, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
								 iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
    } else h->Reset();
    fHistos->AddAt(h, kTRDRefTracksNeg4+iCent);
  }

  // (eta,detector phi,P) distribution of reference TRD negative tracks with 5 tracklets
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH3F*)gROOT->FindObject(Form("hTRDRefTracksNeg5_cent%d",iCent+1)))){
      h = new TH3F(Form("hTRDRefTracksNeg5_cent%d",iCent+1), Form("(#eta,detector #varphi,p) for TRD negative reference tracks, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
								  iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
    } else h->Reset();
    fHistos->AddAt(h, kTRDRefTracksNeg5+iCent);
  }

  // (eta,detector phi,P) distribution of reference TRD negative tracks with 6 tracklets
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TH3F*)gROOT->FindObject(Form("hTRDRefTracksNeg6_cent%d",iCent+1)))){
      h = new TH3F(Form("hTRDRefTracksNeg6_cent%d",iCent+1), Form("(#eta,detector #varphi,p) for TRD negative reference tracks, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
                                                                  iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
    } else h->Reset();
    fHistos->AddAt(h, kTRDRefTracksNeg6+iCent);
  }


  // (eta,detector phi) profile of average number of TRD tracklets/track
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TProfile2D*)gROOT->FindObject(Form("hTRDEtaPhiAvNtrkl_cent%d",iCent+1)))){
      h = new TProfile2D(Form("hTRDEtaPhiAvNtrkl_cent%d",iCent+1), Form("<Ntracklets/track> vs (#eta,detector #varphi) for TRD reference tracks, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
						   iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
    } else h->Reset();
    fHistos->AddAt(h, kTRDEtaPhiAvNtrkl+iCent);
  }

  // (eta,delta phi) profile of average number of TRD tracklets/track
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    if(!(h = (TProfile2D*)gROOT->FindObject(Form("hTRDEtaDeltaPhiAvNtrkl_cent%d",iCent+1)))){
      h = new TProfile2D(Form("hTRDEtaDeltaPhiAvNtrkl_cent%d",iCent+1), Form("<Ntracklets/track> vs (#eta, #Delta#varphi) for TRD reference tracks, centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
							iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, -1.0, 1.0, 50, -0.4*TMath::Pi(), 0.4*TMath::Pi());
    } else h->Reset();
    fHistos->AddAt(h, kTRDEtaDeltaPhiAvNtrkl+iCent);
  }  

  // (eta, detector phi) profile of average tracklet Qtot from slices
  for(Int_t iCent=0; iCent<=5; ++iCent) {
    for(Int_t iLayer=0;iLayer<6;iLayer++) {
      if(!(h = (TProfile2D*)gROOT->FindObject(Form("hTRDEtaPhiAvQtot_Layer%d_cent%d",iLayer,iCent+1)))) {
	h = new TProfile2D(Form("hTRDEtaPhiAvQtot_Layer%d_cent%d",iLayer,iCent+1),
			   Form("<Q_{tot}> vs (#eta, detector #varphi) for TRD reference tracks (layer %d), centrality %d, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
				iLayer, iCent+1, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
      } else h->Reset();
      fHistos->AddAt(h, kTRDEtaPhiAvQtot+iCent*6+iLayer);
    }
  }

  // create a CF container and add it to the list of histograms
  Int_t nbinsCf[kNTrdCfVariables];
  for(Int_t i=0;i<kNTrdCfVariables;++i) nbinsCf[i]=0;
  nbinsCf[kEventVtxZ]         =   12;
  nbinsCf[kEventMult]         =    5;
  nbinsCf[kEventBC]           = 3500;
  nbinsCf[kTrackTOFdeltaBC]   = 2001;
  nbinsCf[kTrackTOFBC]        =  201;
  nbinsCf[kTrackDCAxy]        =   26;
  nbinsCf[kTrackDCAz]         =  101;
  nbinsCf[kTrackCharge]       =    2;
  nbinsCf[kTrackPhi]          =  150;
  nbinsCf[kTrackEta]          =  100;
  nbinsCf[kTrackPt]           =   32;
  nbinsCf[kTrackP]            =   18;
  nbinsCf[kTrackTrdTracklets] =    7;
  nbinsCf[kTrackTrdClusters]  =  200;
  for(Int_t i=0;i<6;++i) nbinsCf[kTrackQtot+i] = 100;
  Double_t evVtxLims[2]      = {-12.,+12.};
  Double_t evMultLims[6]     = {0.0, 700., 1400., 2100., 2800., 3500.};
  Double_t evBCLims[2]       = {-0.5, +3499.5};
  Double_t trkTOFdeltaBClims[2] = {-1000.5, 1000.5};
  Double_t trkTOFBClims[2]   = {-100.5, 100.5};
  Double_t trkDCAxyLims[27]  = {-50.0, -40.0, -30.0, -25.0, -20.0, 
                                -15.0, -10.0,  -8.0,  -6.0,  -4.0, 
                                 -3.0,  -2.0,  -1.0,   0.0,  +1.0,
                                 +2.0,  +3.0,  +4.0,  +6.0,  +8.0,
                                +10.0, +15.0, +20.0, +25.0, +30.0,
                                +40.0, +50.0};    
  Double_t trkDCAzLims[2]    = {-50.5, +50.5};    
  Double_t trkChargeLims[2]  = {-1.5, +1.5};
  Double_t trkPhiLims[2]     = {-1.1*TMath::Pi(), +1.1*TMath::Pi()};
  Double_t trkEtaLims[2]     = {-1.0, +1.0};
  Double_t trkPtLims[33]     = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 
                                2.6, 2.8, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  Double_t trkPLims[19]      = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0,
                                2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 9.0, 12.0};
  Double_t trkTrdNLims[2]    = {-0.5, 6.5};
  Double_t trkTrdNclsLims[2] = {-0.5, 199.5};
  Double_t trkQtotLims[2]    = {0.0, 20.};
  fCfContainer = new AliCFContainer("TrdCfContainer", "TRD ESD CF container", 3, kNTrdCfVariables, nbinsCf);
  fCfContainer->SetBinLimits(kEventVtxZ, evVtxLims[0], evVtxLims[1]);
  fCfContainer->SetBinLimits(kEventMult, evMultLims);
  fCfContainer->SetBinLimits(kEventBC, evBCLims[0], evBCLims[1]);
  fCfContainer->SetBinLimits(kTrackTOFdeltaBC, trkTOFdeltaBClims[0], trkTOFdeltaBClims[1]);
  fCfContainer->SetBinLimits(kTrackTOFBC, trkTOFBClims[0], trkTOFBClims[1]);
  fCfContainer->SetBinLimits(kTrackDCAxy, trkDCAxyLims);
  fCfContainer->SetBinLimits(kTrackDCAz, trkDCAzLims[0], trkDCAzLims[1]);
  fCfContainer->SetBinLimits(kTrackCharge, trkChargeLims[0], trkChargeLims[1]);
  fCfContainer->SetBinLimits(kTrackPhi, trkPhiLims[0], trkPhiLims[1]);
  fCfContainer->SetBinLimits(kTrackEta, trkEtaLims[0], trkEtaLims[1]);
  fCfContainer->SetBinLimits(kTrackPt, trkPtLims);
  fCfContainer->SetBinLimits(kTrackP, trkPLims);
  fCfContainer->SetBinLimits(kTrackTrdTracklets, trkTrdNLims[0], trkTrdNLims[1]);
  fCfContainer->SetBinLimits(kTrackTrdClusters, trkTrdNclsLims[0], trkTrdNclsLims[1]);
  for(Int_t i=0; i<6; ++i) fCfContainer->SetBinLimits(kTrackQtot+i, trkQtotLims[0], trkQtotLims[1]);
  fCfContainer->SetVarTitle(kEventVtxZ, "vtxZ");
  fCfContainer->SetVarTitle(kEventMult, "multiplicity");
  fCfContainer->SetVarTitle(kEventBC, "BC");
  fCfContainer->SetVarTitle(kTrackTOFdeltaBC, "TOFdeltaBC");
  fCfContainer->SetVarTitle(kTrackTOFBC, "TOFBC");
  fCfContainer->SetVarTitle(kTrackDCAxy, "DCAxy");
  fCfContainer->SetVarTitle(kTrackDCAz, "DCAz");
  fCfContainer->SetVarTitle(kTrackCharge, "charge");
  fCfContainer->SetVarTitle(kTrackPhi, "phi");
  fCfContainer->SetVarTitle(kTrackEta, "eta");
  fCfContainer->SetVarTitle(kTrackPt, "pt");
  fCfContainer->SetVarTitle(kTrackP, "P");
  fCfContainer->SetVarTitle(kTrackTrdTracklets, "tracklets");
  fCfContainer->SetVarTitle(kTrackTrdClusters, "clusters");
  for(Int_t i=0; i<6; ++i) fCfContainer->SetVarTitle(kTrackQtot+i, Form("Qtot%d",i));
  
  fCfContainer->SetStepTitle(0, "TPC reference");
  fCfContainer->SetStepTitle(1, "TRD");
  fCfContainer->SetStepTitle(2, "TOF");
  
  //  fHistos->AddAt(fCfContainer, kNhistos);
  fHistos->AddAt(fCfContainer, 0);
  
  return fHistos;
}

//____________________________________________________________________
Bool_t AliTRDcheckESD::Load(const Char_t *file, const Char_t *dir, const Char_t *name)
{
// Load data from performance file

  if(!TFile::Open(file)){
    AliWarning(Form("Couldn't open file %s.", file));
    return kFALSE;
  }
  if(dir){
    if(!gFile->cd(dir)){
      AliWarning(Form("Couldn't cd to %s in %s.", dir, file));
      return kFALSE;
    }
  }
  TObjArray *o(NULL);
  const Char_t *tn=(name ? name : GetName());
  if(!(o = (TObjArray*)gDirectory->Get(tn))){
    AliWarning(Form("Missing histogram container %s.", tn));
    return kFALSE;
  }
  fHistos = (TObjArray*)o->Clone(GetName());
  fCfContainer = (AliCFContainer*)fHistos->At(0);
  gFile->Close();
  return kTRUE;
}

//_______________________________________________________
Bool_t AliTRDcheckESD::PutTrendValue(const Char_t *name, Double_t val)
{
// Dump trending value to default file

  if(!fgFile){
    fgFile = fopen("TRD.Performance.txt", "at");
  }
  fprintf(fgFile, "%s_%s %f\n", GetName(), name, val);
  return kTRUE;
}

//____________________________________________________________________
void AliTRDcheckESD::Terminate(Option_t *)
{
  // Steer post-processing 
  if(!fHistos){
    fHistos = dynamic_cast<TObjArray *>(GetOutputData(1));
    if(!fHistos){
      AliError("Histogram container not found in output");
      return;
    }
  }

  //  fNRefFigures = 15;
  //  return;

  const Char_t *name[kNrefs] = {
    "Ncl", "Eff", "Eloss", "PtResDCA"
  };

  TObjArray *arr(NULL); TGraph *g(NULL);
  if(!fResults){
    fResults = new TObjArray(kNrefs);
    fResults->SetOwner();
    fResults->SetName("results");
    for(Int_t iref(0); iref<kNrefs; iref++){
      fResults->AddAt(arr = new TObjArray(fgkNgraph[iref]), iref);
      arr->SetName(name[iref]);  arr->SetOwner();
      switch(iref+1){
      case kNCl:
	for(Int_t ig(0); ig<fgkNgraph[iref]; ig++){
          arr->AddAt(g = new TGraphErrors(), ig);
          g->SetLineColor(ig+1); 
          g->SetMarkerColor(ig+1); 
          g->SetMarkerStyle(ig+20); 
          g->SetName(Form("s%d", ig));
          switch(ig){
          case 0: g->SetTitle("ALL"); break;
          case 1: g->SetTitle("NEG"); break;
          case 2: g->SetTitle("POS"); break;
          default: g->SetTitle(AliPID::ParticleLatexName(ig-3)); break;
          };
        }
        break;
      case kTRDmom:
	for(Int_t ig(0); ig<fgkNgraph[iref]; ig++){
          arr->AddAt(g = new TGraphAsymmErrors(), ig);
          g->SetLineColor(ig+1); 
          g->SetMarkerColor(ig+1); 
          g->SetMarkerStyle(ig+20); 
        }
        break;
      case kPtRes:
	for(Int_t idx(0); idx<AliPID::kSPECIES; idx++){
          Int_t ig(2*idx);
          arr->AddAt(g = new TGraphErrors(), ig);
          g->SetLineColor(kRed-idx); 
          g->SetMarkerColor(kRed-idx); 
          g->SetMarkerStyle(20+idx); 
          g->SetNameTitle(Form("s%d", ig), Form("res %s", AliPID::ParticleLatexName(idx)));
          arr->AddAt(g = new TGraphErrors(), ig+1);
          g->SetLineColor(kBlue-idx); 
          g->SetMarkerColor(kBlue-idx); 
          g->SetMarkerStyle(20+idx); 
          g->SetNameTitle(Form("m%d", ig+1), Form("sys %s", AliPID::ParticleLatexName(idx)));

          ig+=10;
          arr->AddAt(g = new TGraphErrors(), ig);
          g->SetLineColor(kRed-idx); 
          g->SetMarkerColor(kRed-idx); 
          g->SetMarkerStyle(20+idx); 
          g->SetNameTitle(Form("s%d", ig), Form("sigma %s", AliPID::ParticleLatexName(idx)));
          arr->AddAt(g = new TGraphErrors(), ig+1);
          g->SetLineColor(kBlue-idx); 
          g->SetMarkerColor(kBlue-idx); 
          g->SetMarkerStyle(20+idx); 
          g->SetNameTitle(Form("m%d", ig+1), Form("mean %s", AliPID::ParticleLatexName(idx)));
        }
        break;
      default:
	for(Int_t ig(0); ig<fgkNgraph[iref]; ig++){
          arr->AddAt(g = new TGraphErrors(), ig);
          g->SetLineColor(ig+1); 
          g->SetMarkerColor(ig+1); 
          g->SetMarkerStyle(ig+20); 
        }
        break;
      }
    }
  }
  TH1 *h1[2] = {NULL, NULL};
  TH2I *h2(NULL);
  TAxis *ax(NULL);

  // No of clusters
  if(!(h2 = (TH2I*)fHistos->At(kNCl))) return;
  ax = h2->GetXaxis();
  arr = (TObjArray*)fResults->At(kNCl);
  
  // All tracks
  h1[0] = h2->ProjectionX("Ncl_px");
  TGraphErrors *ge=(TGraphErrors*)arr->At(0);
  for(Int_t ib=2; ib<=ax->GetNbins(); ib++){
    ge->SetPoint(ib-2, ax->GetBinCenter(ib), h1[0]->GetBinContent(ib));
  }
  
  // All charged tracks
  TH1 *hNclCh[2] = {(TH1D*)h1[0]->Clone("NEG"), (TH1D*)h1[0]->Clone("POS")};
  hNclCh[0]->Reset();hNclCh[1]->Reset();
  for(Int_t is(1); is<=AliPID::kSPECIES; is++){
    hNclCh[0]->Add(h2->ProjectionX("Ncl_px", 2*is-1, 2*is-1)); // neg
    hNclCh[1]->Add(h2->ProjectionX("Ncl_px", 2*is, 2*is));     // pos
  }
  if(Int_t(hNclCh[0]->GetEntries())){
    ge=(TGraphErrors*)arr->At(1);
    for(Int_t ib=2; ib<=ax->GetNbins(); ib++){
      ge->SetPoint(ib-2, ax->GetBinCenter(ib), hNclCh[0]->GetBinContent(ib));
    }
  }
  
  if(Int_t(hNclCh[1]->GetEntries())){
    ge=(TGraphErrors*)arr->At(2);
    for(Int_t ib=2; ib<=ax->GetNbins(); ib++){
      ge->SetPoint(ib-2, ax->GetBinCenter(ib), hNclCh[1]->GetBinContent(ib));
    }
  }
  // Species wise
  for(Int_t is(1); is<=AliPID::kSPECIES; is++){
    h1[0] = h2->ProjectionX("Ncl_px", 2*is-1, 2*is);
    if(!Int_t(h1[0]->GetEntries())) continue;
    ge=(TGraphErrors*)arr->At(2+is);
    for(Int_t ib=2; ib<=ax->GetNbins(); ib++){
      ge->SetPoint(ib-2, ax->GetBinCenter(ib), h1[0]->GetBinContent(ib));
    }
  }
  
  fNRefFigures = 1;

  // EFFICIENCY
  // geometrical efficiency
  if(!(h2 = (TH2I*)fHistos->At(kTRDstat))) return;
  arr = (TObjArray*)fResults->At(kTRDstat-1);
  h1[0] = h2->ProjectionX("checkESDx0", kTPCout, kTPCout);
  h1[1] = h2->ProjectionX("checkESDx1", kTRDin, kTRDin);
  Process(h1, (TGraphErrors*)arr->At(0));
  delete h1[0];delete h1[1];
  // tracking efficiency
  h1[0] = h2->ProjectionX("checkESDx0", kTRDin, kTRDin);
  h1[1] = h2->ProjectionX("checkESDx1", kTRDout, kTRDout);
  Process(h1, (TGraphErrors*)arr->At(1));
  delete h1[1];
  // PID efficiency
  h1[1] = h2->ProjectionX("checkESDx1", kTRDpid, kTRDpid);
  Process(h1, (TGraphErrors*)arr->At(2));
  delete h1[1];
  // Refit efficiency
  h1[1] = h2->ProjectionX("checkESDx1", kTRDref, kTRDref);
  Process(h1, (TGraphErrors*)arr->At(3));
  delete h1[1];
  fNRefFigures++;

  // ENERGY LOSS
  if(!(h2 = dynamic_cast<TH2I*>(fHistos->At(kTRDmom)))) return;
  arr = (TObjArray*)fResults->At(kTRDmom-1);
  TGraphAsymmErrors *g06 = (TGraphAsymmErrors*)arr->At(0), *g09 = (TGraphAsymmErrors*)arr->At(1);
  ax=h2->GetXaxis();
  const Int_t nq(4);
  const Double_t xq[nq] = {0.05, 0.2, 0.8, 0.95};
  Double_t yq[nq];
  for(Int_t ily=6; ily--;){
    h1[0] = h2->ProjectionX("checkESDp0", ily+1, ily+1);
    h1[0]->GetQuantiles(nq,yq,xq);
    g06->SetPoint(ily, Float_t(ily), ax->GetBinCenter(h1[0]->GetMaximumBin()));
    g06->SetPointError(ily, 0., 0., TMath::Abs(yq[0]), yq[3]);
    g09->SetPoint(ily, Float_t(ily), h1[0]->GetMean());
    g09->SetPointError(ily, 0., 0., TMath::Abs(yq[1]), yq[2]);

    //printf(" max[%f] mean[%f] q[%f %f %f %f]\n", ax->GetBinCenter(h1[0]->GetMaximumBin()), h1[0]->GetMean(), yq[0], yq[1], yq[2], yq[3]);
    delete h1[0];
  }
  fNRefFigures++;
//  if(!HasMC()) return;

  // Pt RESOLUTION @ DCA
  TH3S* h3(NULL); TGraphErrors *gg[2] = {NULL,NULL};
  if(!(h3 = dynamic_cast<TH3S*>(fHistos->At(kPtRes)))) return;
  arr = (TObjArray*)fResults->At(kPtRes-1);
  TAxis *az(h3->GetZaxis());
  for(Int_t i(0); i<AliPID::kSPECIES; i++){
    Int_t idx(2*i);
    az->SetRange(idx+1, idx+2); 
    gg[1] = (TGraphErrors*)arr->At(idx);
    gg[0] = (TGraphErrors*)arr->At(idx+1);
    Process2D((TH2*)h3->Project3D("yx"), gg);

    idx+=10;
    az->SetRange(idx+1, idx+2); 
    gg[1] = (TGraphErrors*)arr->At(idx);
    gg[0] = (TGraphErrors*)arr->At(idx+1);
    Process2D((TH2*)h3->Project3D("yx"), gg);
  }
  fNRefFigures++;
  
  fNRefFigures++;
  // 3x3 tracking summary canvases for every centrality class
  fNRefFigures++;
  // 3x3 PID summary canvases for every centrality class
  fNRefFigures++;
  // 3x3 for centrality dependent pictures
  fNRefFigures++;
  
  //DoTrending();  
}

//____________________________________________________________________
Int_t AliTRDcheckESD::Pdg2Idx(Int_t pdg) const
{
  //
  // Helper function converting PDG code into AliPID index
  //
  switch(pdg){
  case kElectron: 
  case kPositron: return AliPID::kElectron;  
  case kMuonPlus:
  case kMuonMinus: return AliPID::kMuon;  
  case kPiPlus: 
  case kPiMinus: return AliPID::kPion;  
  case kKPlus: 
  case kKMinus: return AliPID::kKaon;
  case kProton: 
  case kProtonBar: return AliPID::kProton;
  } 
  return -1;
}

//____________________________________________________________________
void AliTRDcheckESD::Process(TH1 **h1, TGraphErrors *g)
{
// Generic function to process one reference plot

  Int_t n1 = 0, n2 = 0, ip=0;
  Double_t eff = 0.;

  TAxis *ax = h1[0]->GetXaxis();
  for(Int_t ib=1; ib<=ax->GetNbins(); ib++){
    if(!(n1 = (Int_t)h1[0]->GetBinContent(ib))) continue;
    n2 = (Int_t)h1[1]->GetBinContent(ib);
    eff = n2/Float_t(n1);

    ip=g->GetN();
    g->SetPoint(ip, ax->GetBinCenter(ib), eff);
    g->SetPointError(ip, 0., n2 ? eff*TMath::Sqrt(1./n1+1./n2) : 0.);
  }
}  
//________________________________________________________
void AliTRDcheckESD::Process2D(TH2 * const h2, TGraphErrors **g)
{
  //
  // Do the processing
  //

  Int_t n = 0;
  if((n=g[0]->GetN())) for(;n--;) g[0]->RemovePoint(n);
  if((n=g[1]->GetN())) for(;n--;) g[1]->RemovePoint(n);
  TF1 f("fg", "gaus", -3.,3.);
  for(Int_t ibin = 1; ibin <= h2->GetNbinsX(); ibin++){
    Double_t x = h2->GetXaxis()->GetBinCenter(ibin);
    TH1D *h = h2->ProjectionY("py", ibin, ibin);
    if(h->GetEntries()<100) continue;
    //AdjustF1(h, f);

    h->Fit(&f, "QN");
    Int_t ip = g[0]->GetN();
    g[0]->SetPoint(ip, x, f.GetParameter(1));
    g[0]->SetPointError(ip, 0., f.GetParError(1));
    g[1]->SetPoint(ip, x, f.GetParameter(2));
    g[1]->SetPointError(ip, 0., f.GetParError(2));
  }
  return;
}
//____________________________________________________________________
void AliTRDcheckESD::PrintStatus(ULong_t status)
{
// Dump track status to stdout

  printf("ITS[i(%d) o(%d) r(%d)] TPC[i(%d) o(%d) r(%d) p(%d)] TRD[i(%d) o(%d) r(%d) p(%d) s(%d)] HMPID[o(%d) p(%d)]\n"
    ,Bool_t(status & AliESDtrack::kITSin)
    ,Bool_t(status & AliESDtrack::kITSout)
    ,Bool_t(status & AliESDtrack::kITSrefit)
    ,Bool_t(status & AliESDtrack::kTPCin)
    ,Bool_t(status & AliESDtrack::kTPCout)
    ,Bool_t(status & AliESDtrack::kTPCrefit)
    ,Bool_t(status & AliESDtrack::kTPCpid)
    ,Bool_t(status & AliESDtrack::kTRDin)
    ,Bool_t(status & AliESDtrack::kTRDout)
    ,Bool_t(status & AliESDtrack::kTRDrefit)
    ,Bool_t(status & AliESDtrack::kTRDpid)
    ,Bool_t(status & AliESDtrack::kTRDStop)
    ,Bool_t(status & AliESDtrack::kHMPIDout)
    ,Bool_t(status & AliESDtrack::kHMPIDpid)
  );
}

//____________________________________________________________________
TH1D* AliTRDcheckESD::Proj2D(TH2* hist) {
  //
  // project the PH vs Slice 2D-histo into a 1D histo
  //
  /*TH1D* hProjection = new TH1F("hProjection","", hist->GetXaxis()->GetXbins()->GetSize()-1, 
			       hist->GetXaxis()->GetXbins()->GetArray());*/
  TH1D* hProjection = (TH1D*)hist->ProjectionX(Form("hProjection_%f", gRandom->Rndm()));
  hProjection->Reset();
  //cout << "Proj2D: nbins = " << hist->GetXaxis()->GetXbins()->GetSize()-1 << endl;
  TF1* fitLandau = new TF1("landauFunc","landau",0.,2000.);
  TH1D *hD;
  for(Int_t iBin=1;iBin<=hist->GetXaxis()->GetNbins();iBin++) {
    if(gROOT->FindObject("projection"))
      delete gROOT->FindObject("projection");
    hD = (TH1D*)hist->ProjectionY("projection",iBin,iBin);
    hD->Rebin(4);
    if(hD->Integral()>10) {
      fitLandau->SetParameter(1, hD->GetBinCenter(hD->GetMaximumBin()));
      fitLandau->SetParLimits(1, 0.2*hD->GetBinCenter(hD->GetMaximumBin()), 3.0*hD->GetBinCenter(hD->GetMaximumBin()));
      fitLandau->SetParameter(0, 1000.);
      fitLandau->SetParLimits(0, 1., 10000000.);
      fitLandau->SetParameter(2, 0.5*hD->GetBinCenter(hD->GetMaximumBin()));
      fitLandau->SetParLimits(2, 0.01*hD->GetBinCenter(hD->GetMaximumBin()), 1.0*hD->GetBinCenter(hD->GetMaximumBin()));
      hD->Fit(fitLandau, "Q0", "", hD->GetXaxis()->GetXmin(), hD->GetXaxis()->GetXmax());
      hD->Fit(fitLandau, "Q0", "", hD->GetXaxis()->GetXmin(), hD->GetXaxis()->GetXmax());
      hProjection->SetBinContent(iBin, fitLandau->GetParameter(1));
      hProjection->SetBinError(iBin, fitLandau->GetParameter(2));
    }
    else{
      hProjection->SetBinContent(iBin, 0);
      hProjection->SetBinError(iBin, 0);
    }
  }
  return hProjection;
}

//____________________________________________________________________
TH2F* AliTRDcheckESD::Proj3D(TH3* hist, TH2* accMap, Int_t zbinLow, Int_t zbinHigh, Float_t &entries) {
  //
  //  Project a 3D histogram to a 2D histogram in the Z axis interval [zbinLow,zbinHigh] 
  //  Return the 2D histogram and also the number of entries into this projection (entries)

  Int_t nBinsX = hist->GetXaxis()->GetNbins();   // X and Y axis bins are assumed to be all equal
  Float_t minX = hist->GetXaxis()->GetXmin();
  Float_t maxX = hist->GetXaxis()->GetXmax();
  Int_t nBinsY = hist->GetYaxis()->GetNbins();
  Float_t minY = hist->GetYaxis()->GetXmin();
  Float_t maxY = hist->GetYaxis()->GetXmax();
  Int_t nBinsZ = hist->GetZaxis()->GetNbins();  // Z axis bins (pt) might have different widths

  TH2F* projHisto = (TH2F*)gROOT->FindObject("projHisto");
  if(projHisto) 
    projHisto->Reset();
  else
    projHisto = new TH2F("projHisto", "projection", nBinsX, minX, maxX, nBinsY, minY, maxY);

  for(Int_t iZ=1; iZ<=nBinsZ; iZ++) {
    if(iZ<zbinLow) continue;
    if(iZ>zbinHigh) continue;
    for(Int_t iX=1; iX<=nBinsX; iX++) {
      for(Int_t iY=1; iY<=nBinsY; iY++) {
        if(accMap) {
          if(accMap->GetBinContent(iX,iY)>0.1)
            projHisto->SetBinContent(iX, iY, projHisto->GetBinContent(iX, iY)+hist->GetBinContent(iX,iY,iZ));
        }
        else    // no acc. cut 
          projHisto->SetBinContent(iX, iY, projHisto->GetBinContent(iX, iY)+hist->GetBinContent(iX,iY,iZ));
        // count only the entries which are inside the acceptance map
        if(accMap) {
          if(accMap->GetBinContent(iX,iY)>0.1)
            entries+=hist->GetBinContent(iX,iY,iZ);
        }
        else    // no acc. cut
          entries+=hist->GetBinContent(iX,iY,iZ);
      }
    }
  }
  return projHisto;
}

//____________________________________________________________________
void AliTRDcheckESD::CheckActiveSM(TH1D* phiProj, Bool_t activeSM[18]) {
  //
  // Check the active super-modules
  //
  Double_t entries[18] = {0.0};
  Double_t smPhiLimits[19];
  for(Int_t ism=0; ism<=18; ++ism) smPhiLimits[ism] = -TMath::Pi() + (2.0*TMath::Pi()/18.0)*ism;
  for(Int_t phiBin=1; phiBin<=phiProj->GetXaxis()->GetNbins(); ++phiBin) {
    Double_t phi = phiProj->GetBinCenter(phiBin);
    Int_t sm = -1;
    for(Int_t ism=0; ism<18; ++ism) 
      if(phi>=smPhiLimits[ism] && phi<smPhiLimits[ism+1]) sm = ism;
    if(sm==-1) continue;
    entries[sm] += phiProj->GetBinContent(phiBin);
  }
  Double_t avEntries = Double_t(phiProj->Integral())/18.0;
  for(Int_t ism=0; ism<18; ++ism) 
    if(entries[ism]>0.5*avEntries) activeSM[ism] = kTRUE;
}

//____________________________________________________________________
TH1F* AliTRDcheckESD::EfficiencyTRD(TH3* tpc3D, TH3* trd3D, Bool_t useAcceptance) {
  //
  // Calculate the TRD-TPC matching efficiency as function of pt
  //
  
  if(!tpc3D || !trd3D) return NULL;
  Int_t nBinsZ = trd3D->GetZaxis()->GetNbins();
  // project everything on the eta-phi map to obtain an acceptance map
  Float_t nada = 0.;
  TH2F *trdAcc = (useAcceptance ? (TH2F*)Proj3D(trd3D, 0x0, 1, nBinsZ, nada)->Clone(Form("trdAcc%f", gRandom->Rndm())) : 0x0);
  TH1D *phiProj = (trdAcc ? trdAcc->ProjectionY(Form("phiProj%f", gRandom->Rndm())) : 0x0);
  
  // prepare the acceptance map
  Bool_t activeSM[18] = {kFALSE};
  Double_t smPhiLimits[19];
  for(Int_t ism=0; ism<=18; ++ism) smPhiLimits[ism] = -TMath::Pi() + (2.0*TMath::Pi()/18.0)*ism;
  if(phiProj) {
    CheckActiveSM(phiProj, activeSM);   // get the active SMs
    trdAcc->Reset();
    // Put 1 entry in every bin which belongs to an active SM
    for(Int_t iY=1; iY<=trdAcc->GetYaxis()->GetNbins(); ++iY) {
      Double_t phi = trdAcc->GetYaxis()->GetBinCenter(iY);
      Bool_t isActive = kFALSE;
      for(Int_t ism=0; ism<18; ++ism) {
        if(phi>=smPhiLimits[ism] && phi<smPhiLimits[ism+1] && activeSM[ism]) {
          isActive = kTRUE;
        }
      }
      if(!isActive) continue;
      for(Int_t iX=1; iX<=trdAcc->GetXaxis()->GetNbins(); ++iX) 
        if(trdAcc->GetXaxis()->GetBinCenter(iX)>=-0.85 && trdAcc->GetXaxis()->GetBinCenter(iX)<=0.85) trdAcc->SetBinContent(iX, iY, 1.0);
    }  // end for over Y(phi) bins
  }  // end if phiProj
    
  // get the bin limits from the Z axis of 3D histos
  Float_t *ptBinLimits = new Float_t[nBinsZ+1];
  for(Int_t i=1; i<=nBinsZ; i++) {
    ptBinLimits[i-1] = trd3D->GetZaxis()->GetBinLowEdge(i);
  }
  ptBinLimits[nBinsZ] = trd3D->GetZaxis()->GetBinUpEdge(nBinsZ);
  
  TH1F *efficiency = new TH1F(Form("eff%d", Int_t(1000000.0*gRandom->Rndm())), "TRD-TPC matching efficiency", nBinsZ, ptBinLimits);
  
  // loop over Z bins
  Bool_t effGood = kFALSE;
  for(Int_t i=1; i<=nBinsZ; i++) {
    Float_t tpcEntries = 0.0; Float_t trdEntries = 0.0;
    Proj3D(tpc3D, trdAcc, i, i, tpcEntries);
    Proj3D(trd3D, trdAcc, i, i, trdEntries);
    Float_t ratio = 0;
    if(tpcEntries>0) ratio = trdEntries/tpcEntries;
    Float_t error = 0;
    if(tpcEntries>0 && trdEntries>0 && (tpcEntries-trdEntries)>=0.0) 
      error = TMath::Sqrt(trdEntries*(tpcEntries-trdEntries)/tpcEntries/tpcEntries/tpcEntries);
    if(ratio>0.001) {
      efficiency->SetBinContent(i,ratio);
      efficiency->SetBinError(i,error);
      effGood = kTRUE;
    }
  }     // end loop over Z bins
  if(!effGood) return 0x0;
  
  return efficiency;
}

//__________________________________________________________________________________________________
void AliTRDcheckESD::PlotCentSummaryFromCF(Double_t* trendValues, Bool_t useIsolatedBC, Bool_t cutTOFbc) {
  //
  // Make the centrality summary figure from the CF container 
  // 
  if(!fCfContainer) return;
  
  TLatex* lat=new TLatex();
  lat->SetTextSize(0.06);
  lat->SetTextColor(2);

  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001); gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
  gPad->Divide(3,3,0.,0.);
  TList* l=gPad->GetListOfPrimitives();
  TVirtualPad* pad=0x0;
  
  fCfContainer->SetRangeUser(kTrackDCAxy, -0.999, +0.999);
  fCfContainer->SetRangeUser(kTrackDCAz, -3.0, +3.0);
  if(cutTOFbc) fCfContainer->SetRangeUser(kTrackTOFBC, 0.0, 0.0);
  
  // find all the isolated bunch crossings with entries
  Bool_t isIsolated[3500];
  TH1D* tempTH1D = (TH1D*)fCfContainer->Project(0, kEventBC);
  FindIsolatedBCs(tempTH1D, isIsolated); delete tempTH1D;
  Int_t nIsolatedBC = 0;
  for(Int_t ibc=0; ibc<3500; ++ibc) 
    if(isIsolated[ibc]) nIsolatedBC++;
    
  if(gROOT->FindObject("rangeEffPt")) delete gROOT->FindObject("rangeEffPt");
  TH2F* rangeEffPt=new TH2F("rangeEffPt", "",10,0.,10.,10,0.,1.3);
  rangeEffPt->SetStats(kFALSE);
  SetStyle(rangeEffPt->GetXaxis(), "p_{T} [GeV/c]", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeEffPt->GetYaxis(), "efficiency", 0.07, 0.8, kTRUE, 0.05);
  
  Int_t padsForEffs[5] = {0,3,6,1,4};
  for(Int_t iCent=1; iCent<6; ++iCent) {
    pad = ((TVirtualPad*)l->At(padsForEffs[iCent-1])); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02); pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    rangeEffPt->Draw();
    TLine line;
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.DrawLine(rangeEffPt->GetXaxis()->GetXmin(), 0.7, rangeEffPt->GetXaxis()->GetXmax(), 0.7);
    line.DrawLine(rangeEffPt->GetXaxis()->GetXmin(), 0.9, rangeEffPt->GetXaxis()->GetXmax(), 0.9);
    
    fCfContainer->SetRangeUser(kEventMult, Double_t(iCent), Double_t(iCent), kTRUE);
    fCfContainer->SetRangeUser(kTrackCharge, +1.0, +1.0);        // positive charges
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);
    
    TH3D* h3PosTPC = (TH3D*)fCfContainer->Project(0, kTrackEta, kTrackPhi, kTrackPt);
    if(h3PosTPC->GetEntries()<10) {
      delete h3PosTPC;
      continue;
    }
    TH3D* h3PosTRDall = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 4.0, 4.0);        // >= 4 TRD tracklets
    TH3D* h3PosTRDtrk4 = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 5.0, 5.0);        // >= 5 TRD tracklets
    TH3D* h3PosTRDtrk5 = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 6.0, 6.0);        // >= 6 TRD tracklets
    TH3D* h3PosTRDtrk6 = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);

    fCfContainer->SetRangeUser(kTrackCharge, -1.0, -1.0);        // negative charges
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);    // >= 0 TRD tracklets
    TH3D* h3NegTPC = (TH3D*)fCfContainer->Project(0, kTrackEta, kTrackPhi, kTrackPt);
    TH3D* h3NegTRDall = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 4.0, 4.0);        // 4 TRD tracklets
    TH3D* h3NegTRDtrk4 = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 5.0, 5.0);        // 5 TRD tracklets
    TH3D* h3NegTRDtrk5 = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 6.0, 6.0);        // 6 TRD tracklets
    TH3D* h3NegTRDtrk6 = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    
    TH1F* hEffPosAll = EfficiencyTRD(h3PosTPC, h3PosTRDall, kTRUE);
    TH1F* hEffPosTrk4 = EfficiencyTRD(h3PosTPC, h3PosTRDtrk4, kTRUE);
    TH1F* hEffPosTrk5 = EfficiencyTRD(h3PosTPC, h3PosTRDtrk5, kTRUE);
    TH1F* hEffPosTrk6 = EfficiencyTRD(h3PosTPC, h3PosTRDtrk6, kTRUE);
    TH1F* hEffNegAll = EfficiencyTRD(h3NegTPC, h3NegTRDall, kTRUE);
    TH1F* hEffNegTrk4 = EfficiencyTRD(h3NegTPC, h3NegTRDtrk4, kTRUE);
    TH1F* hEffNegTrk5 = EfficiencyTRD(h3NegTPC, h3NegTRDtrk5, kTRUE);
    TH1F* hEffNegTrk6 = EfficiencyTRD(h3NegTPC, h3NegTRDtrk6, kTRUE);
    
    delete h3PosTPC; delete h3NegTPC;
    delete h3PosTRDall; delete h3PosTRDtrk4; delete h3PosTRDtrk5; delete h3PosTRDtrk6;
    delete h3NegTRDall; delete h3NegTRDtrk4; delete h3NegTRDtrk5; delete h3NegTRDtrk6;
    
    // get matching efficiencies for isolated bunches
    TH3D* h3TPCrefPos_IsolatedBC=0x0; TH3D* h3TPCrefNeg_IsolatedBC=0x0;
    TH3D* h3TRDrefPosAll_IsolatedBC=0x0; TH3D* h3TRDrefPosTrk4_IsolatedBC=0x0; TH3D* h3TRDrefPosTrk5_IsolatedBC=0x0; TH3D* h3TRDrefPosTrk6_IsolatedBC=0x0;
    TH3D* h3TRDrefNegAll_IsolatedBC=0x0; TH3D* h3TRDrefNegTrk4_IsolatedBC=0x0; TH3D* h3TRDrefNegTrk5_IsolatedBC=0x0; TH3D* h3TRDrefNegTrk6_IsolatedBC=0x0;
    for(Int_t ibc=0; ibc<3500; ++ibc) {
      if(!isIsolated[ibc]) continue;
      fCfContainer->SetRangeUser(kEventBC, Double_t(ibc), Double_t(ibc));
    
      TH3D* tempTH3D;
      fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);
      fCfContainer->SetRangeUser(kTrackCharge, +1.0, +1.0);    // positive charges
      tempTH3D = (TH3D*)fCfContainer->Project(0, kTrackEta, kTrackPhi, kTrackPt);
      if(!h3TPCrefPos_IsolatedBC) h3TPCrefPos_IsolatedBC = tempTH3D;
      else {h3TPCrefPos_IsolatedBC->Add(tempTH3D); delete tempTH3D;};
      tempTH3D = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
      if(!h3TRDrefPosAll_IsolatedBC) h3TRDrefPosAll_IsolatedBC = tempTH3D;
      else {h3TRDrefPosAll_IsolatedBC->Add(tempTH3D); delete tempTH3D;}
      fCfContainer->SetRangeUser(kTrackTrdTracklets, 4.0, 4.0);
      tempTH3D = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
      if(!h3TRDrefPosTrk4_IsolatedBC) h3TRDrefPosTrk4_IsolatedBC = tempTH3D;
      else {h3TRDrefPosTrk4_IsolatedBC->Add(tempTH3D); delete tempTH3D;}
      fCfContainer->SetRangeUser(kTrackTrdTracklets, 5.0, 5.0);
      tempTH3D = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
      if(!h3TRDrefPosTrk5_IsolatedBC) h3TRDrefPosTrk5_IsolatedBC = tempTH3D;
      else {h3TRDrefPosTrk5_IsolatedBC->Add(tempTH3D); delete tempTH3D;}
      fCfContainer->SetRangeUser(kTrackTrdTracklets, 6.0, 6.0);
      tempTH3D = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
      if(!h3TRDrefPosTrk6_IsolatedBC) h3TRDrefPosTrk6_IsolatedBC = tempTH3D;
      else {h3TRDrefPosTrk6_IsolatedBC->Add(tempTH3D); delete tempTH3D;}
          
      fCfContainer->SetRangeUser(kTrackCharge, -1.0, -1.0);   // negative charges
      fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);
      tempTH3D = (TH3D*)fCfContainer->Project(0, kTrackEta, kTrackPhi, kTrackPt);
      if(!h3TPCrefNeg_IsolatedBC) h3TPCrefNeg_IsolatedBC = tempTH3D;
      else {h3TPCrefNeg_IsolatedBC->Add(tempTH3D); delete tempTH3D;}
      tempTH3D = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
      if(!h3TRDrefNegAll_IsolatedBC) h3TRDrefNegAll_IsolatedBC = tempTH3D;
      else {h3TRDrefNegAll_IsolatedBC->Add(tempTH3D); delete tempTH3D;}
      fCfContainer->SetRangeUser(kTrackTrdTracklets, 4.0, 4.0);
      tempTH3D = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
      if(!h3TRDrefNegTrk4_IsolatedBC) h3TRDrefNegTrk4_IsolatedBC = tempTH3D;
      else {h3TRDrefNegTrk4_IsolatedBC->Add(tempTH3D); delete tempTH3D;}
      fCfContainer->SetRangeUser(kTrackTrdTracklets, 5.0, 5.0);
      tempTH3D = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
      if(!h3TRDrefNegTrk5_IsolatedBC) h3TRDrefNegTrk5_IsolatedBC = tempTH3D;
      else {h3TRDrefNegTrk5_IsolatedBC->Add(tempTH3D); delete tempTH3D;}
      fCfContainer->SetRangeUser(kTrackTrdTracklets, 6.0, 6.0);
      tempTH3D = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
      if(!h3TRDrefNegTrk6_IsolatedBC) h3TRDrefNegTrk6_IsolatedBC = tempTH3D;
      else {h3TRDrefNegTrk6_IsolatedBC->Add(tempTH3D); delete tempTH3D;}
    }
    fCfContainer->SetRangeUser(kEventBC, 0.0, 3500.0);
    fCfContainer->SetRangeUser(kTrackCharge, -1.0, +1.0);
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);
    
    TH1F* hEffPosAll_IsolatedBC  = EfficiencyTRD(h3TPCrefPos_IsolatedBC, h3TRDrefPosAll_IsolatedBC, kTRUE);
    TH1F* hEffPosTrk4_IsolatedBC = EfficiencyTRD(h3TPCrefPos_IsolatedBC, h3TRDrefPosTrk4_IsolatedBC, kTRUE);
    TH1F* hEffPosTrk5_IsolatedBC = EfficiencyTRD(h3TPCrefPos_IsolatedBC, h3TRDrefPosTrk5_IsolatedBC, kTRUE);
    TH1F* hEffPosTrk6_IsolatedBC = EfficiencyTRD(h3TPCrefPos_IsolatedBC, h3TRDrefPosTrk6_IsolatedBC, kTRUE);
    TH1F* hEffNegAll_IsolatedBC  = EfficiencyTRD(h3TPCrefNeg_IsolatedBC, h3TRDrefNegAll_IsolatedBC, kTRUE);
    TH1F* hEffNegTrk4_IsolatedBC = EfficiencyTRD(h3TPCrefNeg_IsolatedBC, h3TRDrefNegTrk4_IsolatedBC, kTRUE);
    TH1F* hEffNegTrk5_IsolatedBC = EfficiencyTRD(h3TPCrefNeg_IsolatedBC, h3TRDrefNegTrk5_IsolatedBC, kTRUE);
    TH1F* hEffNegTrk6_IsolatedBC = EfficiencyTRD(h3TPCrefNeg_IsolatedBC, h3TRDrefNegTrk6_IsolatedBC, kTRUE);
    
    if(h3TPCrefPos_IsolatedBC) delete h3TPCrefPos_IsolatedBC;
    if(h3TPCrefNeg_IsolatedBC) delete h3TPCrefNeg_IsolatedBC;
    if(h3TRDrefPosAll_IsolatedBC) delete h3TRDrefPosAll_IsolatedBC;
    if(h3TRDrefPosTrk4_IsolatedBC) delete h3TRDrefPosTrk4_IsolatedBC;
    if(h3TRDrefPosTrk5_IsolatedBC) delete h3TRDrefPosTrk5_IsolatedBC;
    if(h3TRDrefPosTrk6_IsolatedBC) delete h3TRDrefPosTrk6_IsolatedBC;
    if(h3TRDrefNegAll_IsolatedBC) delete h3TRDrefNegAll_IsolatedBC;
    if(h3TRDrefNegTrk4_IsolatedBC) delete h3TRDrefNegTrk4_IsolatedBC;
    if(h3TRDrefNegTrk5_IsolatedBC) delete h3TRDrefNegTrk5_IsolatedBC;
    if(h3TRDrefNegTrk6_IsolatedBC) delete h3TRDrefNegTrk6_IsolatedBC;
    
    if(!useIsolatedBC) {
      SetStyle(hEffPosAll,  1, kRed, 1, 24, kRed, 1);
      SetStyle(hEffPosTrk4, 1, kRed, 1, 25, kRed, 1);
      SetStyle(hEffPosTrk5, 1, kRed, 1, 26, kRed, 1);
      SetStyle(hEffPosTrk6, 1, kRed, 1, 27, kRed, 1);
      SetStyle(hEffNegAll,  1, kBlue, 1, 24, kBlue, 1);
      SetStyle(hEffNegTrk4, 1, kBlue, 1, 25, kBlue, 1);
      SetStyle(hEffNegTrk5, 1, kBlue, 1, 26, kBlue, 1);
      SetStyle(hEffNegTrk6, 1, kBlue, 1, 27, kBlue, 1);
      hEffPosAll->Draw("same");
      hEffNegAll->Draw("same");
      hEffPosTrk4->Draw("same");
      hEffNegTrk4->Draw("same");
      hEffPosTrk5->Draw("same");
      hEffNegTrk5->Draw("same");
      hEffPosTrk6->Draw("same");
      hEffNegTrk6->Draw("same");
    }
    else {
      if(nIsolatedBC>0) { 
        SetStyle(hEffPosAll_IsolatedBC,  1, kRed, 1, 24, kRed, 1);
        SetStyle(hEffPosTrk4_IsolatedBC, 1, kRed, 1, 25, kRed, 1);
        SetStyle(hEffPosTrk5_IsolatedBC, 1, kRed, 1, 26, kRed, 1);
        SetStyle(hEffPosTrk6_IsolatedBC, 1, kRed, 1, 27, kRed, 1);
        SetStyle(hEffNegAll_IsolatedBC,  1, kBlue, 1, 24, kBlue, 1);
        SetStyle(hEffNegTrk4_IsolatedBC, 1, kBlue, 1, 25, kBlue, 1);
        SetStyle(hEffNegTrk5_IsolatedBC, 1, kBlue, 1, 26, kBlue, 1);
        SetStyle(hEffNegTrk6_IsolatedBC, 1, kBlue, 1, 27, kBlue, 1);
        hEffPosAll_IsolatedBC->Draw("same");
        hEffNegAll_IsolatedBC->Draw("same");
        hEffPosTrk4_IsolatedBC->Draw("same");
        hEffNegTrk4_IsolatedBC->Draw("same");
        hEffPosTrk5_IsolatedBC->Draw("same");
        hEffNegTrk5_IsolatedBC->Draw("same");
        hEffPosTrk6_IsolatedBC->Draw("same");
        hEffNegTrk6_IsolatedBC->Draw("same");
      }
    }    
        
    TLegend* leg=new TLegend(0.18, 0.7, 0.77, 0.89);
    if(iCent==1) {
      leg->SetFillColor(0);
      leg->SetNColumns(2);
      leg->SetMargin(0.1);
      leg->SetBorderSize(0);
      if(useIsolatedBC) {
        leg->SetHeader("Isolated bunch crossings");
        if(nIsolatedBC) {
          leg->AddEntry(hEffPosAll_IsolatedBC,  "pos. (#geq 1 tracklet)", "p");
          leg->AddEntry(hEffNegAll_IsolatedBC,  "neg. (#geq 1 tracklet)", "p");
          leg->AddEntry(hEffPosTrk4_IsolatedBC, "pos. (4 tracklets)", "p");
          leg->AddEntry(hEffNegTrk4_IsolatedBC, "neg. (4 tracklets)", "p");
          leg->AddEntry(hEffPosTrk5_IsolatedBC, "pos. (5 tracklets)", "p");
          leg->AddEntry(hEffNegTrk5_IsolatedBC, "neg. (5 tracklets)", "p");
          leg->AddEntry(hEffPosTrk6_IsolatedBC, "pos. (6 tracklets)", "p");     
          leg->AddEntry(hEffNegTrk6_IsolatedBC, "neg. (6 tracklets)", "p");
        }
      }
      else {
        leg->SetHeader("All bunch crossings");
        leg->AddEntry(hEffPosAll,  "pos. (#geq 1 tracklet)", "p");
        leg->AddEntry(hEffNegAll,  "neg. (#geq 1 tracklet)", "p");
        leg->AddEntry(hEffPosTrk4, "pos. (4 tracklets)", "p");
        leg->AddEntry(hEffNegTrk4, "neg. (4 tracklets)", "p");
        leg->AddEntry(hEffPosTrk5, "pos. (5 tracklets)", "p");
        leg->AddEntry(hEffNegTrk5, "neg. (5 tracklets)", "p");
        leg->AddEntry(hEffPosTrk6, "pos. (6 tracklets)", "p");     
        leg->AddEntry(hEffNegTrk6, "neg. (6 tracklets)", "p");
      }
      leg->Draw();
    }
    lat->DrawLatex(0.2, 1.32, Form("Centrality class %d", iCent));
  }   // end for loop over multiplicity classes
  
  // Reset the modified user ranges of the CF container
  fCfContainer->SetRangeUser(kEventMult, 0, 6, kTRUE);
  fCfContainer->SetRangeUser(kTrackCharge, -1.0, +1.0);
  fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);
   
  // Cluster distributions in all multiplicity classes
  pad = ((TVirtualPad*)l->At(2)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
  pad->SetTopMargin(0.02); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  if(gROOT->FindObject("rangeNcls")) delete gROOT->FindObject("rangeNcls");
  TH2F* rangeNcls = new TH2F("rangeNcls", "", 10, 0.0, 199.9, 10, 0.0, 1.199);
  SetStyle(rangeNcls->GetXaxis(), "# TRD clusters", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeNcls->GetYaxis(), "entries (a.u.)", 0.07, 0.8, kTRUE, 0.05);
  rangeNcls->SetStats(kFALSE);
  rangeNcls->Draw();
    
  TH1D* hNcls[6]={0x0};
  TLegend* legCls=new TLegend(0.7, 0.75, 0.97, 0.97);
  legCls->SetBorderSize(0);
  legCls->SetFillColor(0);
  legCls->SetMargin(0.15);
  for(Int_t iCent=0; iCent<6; ++iCent) {
    if(iCent>0)
      fCfContainer->SetRangeUser(kEventMult, Double_t(iCent), Double_t(iCent), kTRUE);
    if(!useIsolatedBC) hNcls[iCent] = (TH1D*)fCfContainer->Project(1, kTrackTrdClusters);
    else{
      for(Int_t ibc=0; ibc<3500; ++ibc) {
        if(!isIsolated[ibc]) continue;
        fCfContainer->SetRangeUser(kEventBC, Double_t(ibc), Double_t(ibc));
        tempTH1D = (TH1D*)fCfContainer->Project(1, kTrackTrdClusters);
        if(!hNcls[iCent]) hNcls[iCent] = tempTH1D;
        else {hNcls[iCent]->Add(tempTH1D); delete tempTH1D;}
      }
      fCfContainer->SetRangeUser(kEventBC, 0.0, 3500.0);
    }
    if(!hNcls[iCent]) continue;
    hNcls[iCent]->SetLineColor(iCent<4 ? iCent+1 : iCent+2);
    Double_t maximum = hNcls[iCent]->GetMaximum();
    if(maximum>1.0)
      hNcls[iCent]->Scale(1.0/maximum);
    hNcls[iCent]->SetStats(kFALSE);
    hNcls[iCent]->SetTitle("");
    hNcls[iCent]->SetLineWidth(2);
    
    if(hNcls[iCent]->Integral()>0.01) {
      hNcls[iCent]->Draw("same");
      legCls->AddEntry(hNcls[iCent], (iCent==0 ? "all centralities" : Form("centrality class %d", iCent)), "l");
    }
  }
  if(useIsolatedBC) legCls->SetHeader("Isolated bunch crossings");
  else legCls->SetHeader("All bunch crossings");
  legCls->Draw();
  fCfContainer->SetRangeUser(kEventMult, 0.0, 6.0, kTRUE);
  
  // Qtot vs P
  pad = ((TVirtualPad*)l->At(5)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
  pad->SetTopMargin(0.02); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  if(gROOT->FindObject("rangeQtot")) delete gROOT->FindObject("rangeQtot");
  TH2F* rangeQtot = new TH2F("rangeQtot", "", 10, 0.0, 9.999, 10, 0.0, 1.199);
  SetStyle(rangeQtot->GetXaxis(), "Q_{tot} (a.u.)", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeQtot->GetYaxis(), "entries (a.u.)", 0.07, 0.8, kTRUE, 0.05);
  rangeQtot->SetStats(kFALSE);
  rangeQtot->Draw();
  
  TH1D* hQtot[6]={0x0};
  TLegend* leg2=new TLegend(0.6, 0.7, 0.9, 0.97);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  for(Int_t iCent=0; iCent<6; ++iCent) {
    if(iCent>0)
      fCfContainer->SetRangeUser(kEventMult, Double_t(iCent), Double_t(iCent), kTRUE);
    if(useIsolatedBC) {
      for(Int_t ibc=0; ibc<3500; ++ibc) {
        if(!isIsolated[ibc]) continue;
        fCfContainer->SetRangeUser(kEventBC, Double_t(ibc), Double_t(ibc));
        for(Int_t il=0; il<6; ++il) {
          tempTH1D = (TH1D*)fCfContainer->Project(1, kTrackQtot+il);
          if(!hQtot[iCent] && il==0) hQtot[iCent] = tempTH1D;
          else {hQtot[iCent]->Add(tempTH1D); delete tempTH1D;}
        }
      }
      fCfContainer->SetRangeUser(kEventBC, 0.0, 3500.0);
    }
    else {
      for(Int_t il=0; il<6; ++il) {
        tempTH1D = (TH1D*)fCfContainer->Project(1, kTrackQtot+il);
        if(il==0) hQtot[iCent] = tempTH1D;
        else hQtot[iCent]->Add(tempTH1D);
      }
    } // end if(useIsolatedBC)
    
    if(!hQtot[iCent]) continue;
    hQtot[iCent]->SetBinContent(1, 0);
    Double_t maximum = hQtot[iCent]->GetMaximum();
    if(maximum>1.0)
      hQtot[iCent]->Scale(1.0/maximum);
    hQtot[iCent]->SetLineColor(iCent<4 ? iCent+1 : iCent+2);
    hQtot[iCent]->SetStats(kFALSE);
    hQtot[iCent]->SetTitle("");
    hQtot[iCent]->SetLineWidth(2);
    if(hQtot[iCent]->Integral()>0.01) {
      hQtot[iCent]->Draw(iCent==0 ? "" : "same");
      leg2->AddEntry(hQtot[iCent], (iCent==0 ? "all centralities" : Form("centrality class %d", iCent)), "l");
    }
  }
  if(useIsolatedBC) leg2->SetHeader("Isolated bunch crossings");
  else leg2->SetHeader("All bunch crossings");
  leg2->Draw();
  fCfContainer->SetRangeUser(kEventMult, 0.0, 5.0, kTRUE);
  if(cutTOFbc) fCfContainer->SetRangeUser(kTrackTOFBC, -1000.0, +1000.0);  // reset the cut on TOFbc
}


//_________________________________________________________________
void AliTRDcheckESD::PlotTrackingSummaryFromCF(Int_t centralityClass, Double_t* trendValues, Bool_t useIsolatedBC, Bool_t cutTOFbc) {

  if(!fCfContainer) return;
  
  TLatex *lat=new TLatex();
  lat->SetTextSize(0.06);
  lat->SetTextColor(2);
  
  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001);
  gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
  gPad->Divide(3,3,0.,0.);
  TList* l=gPad->GetListOfPrimitives();
  
  // eta-phi distr. for positive TPC tracks
  TVirtualPad* pad = ((TVirtualPad*)l->At(0)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
 
  fCfContainer->SetRangeUser(kTrackDCAxy, -0.999, +0.999);
  fCfContainer->SetRangeUser(kTrackDCAz, -3.0, +3.0);
  if(cutTOFbc) fCfContainer->SetRangeUser(kTrackTOFBC, 0.0, 0.0);
  
  // find all the isolated bunch crossings with entries
  Bool_t isIsolated[3500];
  TH1D* tempTH1D = (TH1D*)fCfContainer->Project(0, kEventBC);
  FindIsolatedBCs(tempTH1D, isIsolated); delete tempTH1D;
  Int_t nIsolatedBC = 0;
  for(Int_t ibc=0; ibc<3500; ++ibc) 
    if(isIsolated[ibc]) nIsolatedBC++;
  
  if(centralityClass>0) // select the multiplicity class
    fCfContainer->SetRangeUser(kEventMult, Double_t(centralityClass), Double_t(centralityClass), kTRUE);
  
  TH2D* hTPCrefPos = 0x0; TH2D* hTRDrefPos = 0x0; TH2D* hTOFrefPos = 0x0;
  TH2D* hTPCrefNeg = 0x0; TH2D* hTRDrefNeg = 0x0; TH2D* hTOFrefNeg = 0x0;
  if(!useIsolatedBC) {
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);
    fCfContainer->SetRangeUser(kTrackCharge, +1.0, +1.0);      // positive charges
    hTPCrefPos = (TH2D*)fCfContainer->Project(0, kTrackEta, kTrackPhi);
    hTRDrefPos = (TH2D*)fCfContainer->Project(1, kTrackEta, kTrackPhi);
    hTOFrefPos = (TH2D*)fCfContainer->Project(2, kTrackEta, kTrackPhi);
    fCfContainer->SetRangeUser(kTrackCharge, -1.0, -1.0);      // negative charges
    hTPCrefNeg = (TH2D*)fCfContainer->Project(0, kTrackEta, kTrackPhi);
    hTRDrefNeg = (TH2D*)fCfContainer->Project(1, kTrackEta, kTrackPhi);
    hTOFrefNeg = (TH2D*)fCfContainer->Project(2, kTrackEta, kTrackPhi);
    fCfContainer->SetRangeUser(kTrackCharge, -1.0, +1.0);      // reset charge cut
  }
  
  TH2D* hTPCrefPos_IsolatedBC=0x0; TH2D* hTPCrefNeg_IsolatedBC=0x0;
  TH2D* hTRDrefPos_IsolatedBC=0x0; TH2D* hTRDrefNeg_IsolatedBC=0x0;
  TH2D* hTOFrefPos_IsolatedBC=0x0; TH2D* hTOFrefNeg_IsolatedBC=0x0;
  if(useIsolatedBC) {
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);
    for(Int_t ibc=0; ibc<3500; ++ibc) {
      if(!isIsolated[ibc]) continue;
      fCfContainer->SetRangeUser(kEventBC, Double_t(ibc), Double_t(ibc));
      fCfContainer->SetRangeUser(kTrackCharge, +1.0, +1.0);      // positive charges
      TH2D* tempTH2D;
      tempTH2D = (TH2D*)fCfContainer->Project(0, kTrackEta, kTrackPhi);
      if(!hTPCrefPos_IsolatedBC) hTPCrefPos_IsolatedBC = tempTH2D;
      else {hTPCrefPos_IsolatedBC->Add(tempTH2D); delete tempTH2D;}
      tempTH2D = (TH2D*)fCfContainer->Project(1, kTrackEta, kTrackPhi);
      if(!hTRDrefPos_IsolatedBC) hTRDrefPos_IsolatedBC = tempTH2D;
      else {hTRDrefPos_IsolatedBC->Add(tempTH2D); delete tempTH2D;}
      tempTH2D = (TH2D*)fCfContainer->Project(2, kTrackEta, kTrackPhi);
      if(!hTOFrefPos_IsolatedBC) hTOFrefPos_IsolatedBC = tempTH2D;
      else {hTOFrefPos_IsolatedBC->Add(tempTH2D); delete tempTH2D;}
      fCfContainer->SetRangeUser(kTrackCharge, -1.0, -1.0);      // negative charges
      tempTH2D = (TH2D*)fCfContainer->Project(0, kTrackEta, kTrackPhi);
      if(!hTPCrefNeg_IsolatedBC) hTPCrefNeg_IsolatedBC = tempTH2D;
      else {hTPCrefNeg_IsolatedBC->Add(tempTH2D); delete tempTH2D;}
      tempTH2D = (TH2D*)fCfContainer->Project(1, kTrackEta, kTrackPhi);
      if(!hTRDrefNeg_IsolatedBC) hTRDrefNeg_IsolatedBC = tempTH2D;
      else {hTRDrefNeg_IsolatedBC->Add(tempTH2D); delete tempTH2D;}
      tempTH2D = (TH2D*)fCfContainer->Project(2, kTrackEta, kTrackPhi);
      if(!hTOFrefNeg_IsolatedBC) hTOFrefNeg_IsolatedBC = tempTH2D;
      else {hTOFrefNeg_IsolatedBC->Add(tempTH2D); delete tempTH2D;}
    }
    fCfContainer->SetRangeUser(kEventBC, 0.0, 3500.);   // reset the BC range
    fCfContainer->SetRangeUser(kTrackCharge, -1.0, +1.0);      // reset charge cut
  }
  
  if(gROOT->FindObject("rangeEtaPhi")) delete gROOT->FindObject("rangeEtaPhi");
  TH2F* rangeEtaPhi = new TH2F("rangeEtaPhi", "", 10, -0.99, +0.99, 10, -3.4, +3.4);
  SetStyle(rangeEtaPhi->GetXaxis(), "#eta", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeEtaPhi->GetYaxis(), "detector #varphi", 0.07, 0.8, kTRUE, 0.05);
  rangeEtaPhi->SetStats(kFALSE);  
  
  //----------------------------------------------
  // eta-phi efficiency for positive TRD tracks
  pad = ((TVirtualPad*)l->At(0)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  
  if(useIsolatedBC) {
    if(nIsolatedBC>0) {
      TH2D* hTRDeffPos_IsolatedBC = (TH2D*)hTRDrefPos_IsolatedBC->Clone("hTRDeffPos_IsolatedBC");
      hTRDeffPos_IsolatedBC->Reset();
      hTRDeffPos_IsolatedBC->SetStats(kFALSE);
      hTRDeffPos_IsolatedBC->Divide(hTRDrefPos_IsolatedBC, hTPCrefPos_IsolatedBC);  
      hTRDeffPos_IsolatedBC->SetMaximum(1.0);
      hTRDeffPos_IsolatedBC->Draw("samecolz");
    }
  }
  else {
    TH2D* hTRDeffPos = (TH2D*)hTRDrefPos->Clone("hTRDeffPos");
    hTRDeffPos->Reset();
    hTRDeffPos->SetStats(kFALSE);
    hTRDeffPos->Divide(hTRDrefPos, hTPCrefPos);
    hTRDeffPos->SetMaximum(1.0);
    hTRDeffPos->Draw("samecolz");
  }
  lat->DrawLatex(-0.9, 3.6, "TPC-TRD matching for positive tracks");
  DrawTRDGrid();
  
  //----------------------------------------------
  // eta-phi efficiency for negative TRD tracks
  pad = ((TVirtualPad*)l->At(3)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  
  if(useIsolatedBC) {
    if(nIsolatedBC>0) {
      TH2D* hTRDeffNeg_IsolatedBC = (TH2D*)hTRDrefNeg_IsolatedBC->Clone("hTRDeffNeg_IsolatedBC");
      hTRDeffNeg_IsolatedBC->Reset();
      hTRDeffNeg_IsolatedBC->SetStats(kFALSE);
      hTRDeffNeg_IsolatedBC->Divide(hTRDrefNeg_IsolatedBC, hTPCrefNeg_IsolatedBC);  
      hTRDeffNeg_IsolatedBC->SetMaximum(1.0);
      hTRDeffNeg_IsolatedBC->Draw("samecolz");
    }
  }
  else {
    TH2D* hTRDeffNeg = (TH2D*)hTRDrefNeg->Clone("hTRDeffNeg");
    hTRDeffNeg->Reset();
    hTRDeffNeg->SetStats(kFALSE);
    hTRDeffNeg->Divide(hTRDrefNeg, hTPCrefNeg);
    hTRDeffNeg->SetMaximum(1.0);
    hTRDeffNeg->Draw("samecolz");
  }
  lat->DrawLatex(-0.9, 3.6, "TPC-TRD matching for negative tracks");
  DrawTRDGrid();  
  
  //----------------------------------------------
  // eta-phi TRD-TOF matching efficiency for positive tracks
  pad = ((TVirtualPad*)l->At(1)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  
  if(useIsolatedBC) {
    if(nIsolatedBC>0) {
      TH2D* hTOFeffPos_IsolatedBC = (TH2D*)hTOFrefPos_IsolatedBC->Clone("hTOFeffPos_IsolatedBC");
      hTOFeffPos_IsolatedBC->Reset();
      hTOFeffPos_IsolatedBC->SetStats(kFALSE);
      hTOFeffPos_IsolatedBC->Divide(hTOFrefPos_IsolatedBC, hTRDrefPos_IsolatedBC);  
      hTOFeffPos_IsolatedBC->SetMaximum(1.0);
      hTOFeffPos_IsolatedBC->Draw("samecolz");
    }
  }
  else {
    TH2D* hTOFeffPos = (TH2D*)hTOFrefPos->Clone("hTOFeffPos");
    hTOFeffPos->Reset();
    hTOFeffPos->SetStats(kFALSE);
    hTOFeffPos->Divide(hTOFrefPos, hTRDrefPos);
    hTOFeffPos->SetMaximum(1.0);
    hTOFeffPos->Draw("samecolz");
  }
  lat->DrawLatex(-0.9, 3.6, "TRD-TOF matching for positive tracks");
  DrawTRDGrid();
  
  //----------------------------------------------
  // eta-phi TRD-TOF matching efficiency for negative tracks
  pad = ((TVirtualPad*)l->At(4)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  
  if(useIsolatedBC) {
    if(nIsolatedBC) {
      TH2D* hTOFeffNeg_IsolatedBC = (TH2D*)hTOFrefNeg_IsolatedBC->Clone("hTOFeffNeg_IsolatedBC");
      hTOFeffNeg_IsolatedBC->Reset();
      hTOFeffNeg_IsolatedBC->SetStats(kFALSE);
      hTOFeffNeg_IsolatedBC->Divide(hTOFrefNeg_IsolatedBC, hTRDrefNeg_IsolatedBC);  
      hTOFeffNeg_IsolatedBC->SetMaximum(1.0);
      hTOFeffNeg_IsolatedBC->Draw("samecolz");
    }
  }
  else {
    TH2D* hTOFeffNeg = (TH2D*)hTOFrefNeg->Clone("hTOFeffNeg");
    hTOFeffNeg->Reset();
    hTOFeffNeg->SetStats(kFALSE);
    hTOFeffNeg->Divide(hTOFrefNeg, hTRDrefNeg);
    hTOFeffNeg->SetMaximum(1.0);
    hTOFeffNeg->Draw("samecolz");
  }
  lat->DrawLatex(-0.9, 3.6, "TRD-TOF matching for negative tracks");
  DrawTRDGrid();
  
  if(hTRDrefPos) delete hTRDrefPos; if(hTPCrefPos) delete hTPCrefPos; if(hTOFrefPos) delete hTOFrefPos;
  if(hTRDrefNeg) delete hTRDrefNeg; if(hTPCrefNeg) delete hTPCrefNeg; if(hTOFrefNeg) delete hTOFrefNeg;
  if(hTRDrefPos_IsolatedBC) delete hTRDrefPos_IsolatedBC; if(hTPCrefPos_IsolatedBC) delete hTPCrefPos_IsolatedBC; if(hTOFrefPos_IsolatedBC) delete hTOFrefPos_IsolatedBC;
  if(hTRDrefNeg_IsolatedBC) delete hTRDrefNeg_IsolatedBC; if(hTPCrefNeg_IsolatedBC) delete hTPCrefNeg_IsolatedBC; if(hTOFrefNeg_IsolatedBC) delete hTOFrefNeg_IsolatedBC;
  
  fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);
  fCfContainer->SetRangeUser(kTrackCharge, +1.0, +1.0);    // positive charges
  TH3D* h3TPCrefPos = (TH3D*)fCfContainer->Project(0, kTrackEta, kTrackPhi, kTrackPt);
  TH3D* h3TRDrefPosAll = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
  TH3D* h3TOFrefPosAll = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
  fCfContainer->SetRangeUser(kTrackTrdTracklets, 4.0, 4.0);
  TH3D* h3TRDrefPosTrk4 = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
  TH3D* h3TOFrefPosTrk4 = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
  fCfContainer->SetRangeUser(kTrackTrdTracklets, 5.0, 5.0);
  TH3D* h3TRDrefPosTrk5 = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
  TH3D* h3TOFrefPosTrk5 = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
  fCfContainer->SetRangeUser(kTrackTrdTracklets, 6.0, 6.0);
  TH3D* h3TRDrefPosTrk6 = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
  TH3D* h3TOFrefPosTrk6 = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
  
  fCfContainer->SetRangeUser(kTrackCharge, -1.0, -1.0);   // negative charges
  fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);
  TH3D* h3TPCrefNeg = (TH3D*)fCfContainer->Project(0, kTrackEta, kTrackPhi, kTrackPt);
  TH3D* h3TRDrefNegAll = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
  TH3D* h3TOFrefNegAll = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
  fCfContainer->SetRangeUser(kTrackTrdTracklets, 4.0, 4.0);
  TH3D* h3TRDrefNegTrk4 = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
  TH3D* h3TOFrefNegTrk4 = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
  fCfContainer->SetRangeUser(kTrackTrdTracklets, 5.0, 5.0);
  TH3D* h3TRDrefNegTrk5 = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
  TH3D* h3TOFrefNegTrk5 = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
  fCfContainer->SetRangeUser(kTrackTrdTracklets, 6.0, 6.0);
  TH3D* h3TRDrefNegTrk6 = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
  TH3D* h3TOFrefNegTrk6 = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
  fCfContainer->SetRangeUser(kTrackCharge, -1.0, +1.0);
  fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);
  
  TH1F* hTRDEffPtPosAll = EfficiencyTRD(h3TPCrefPos, h3TRDrefPosAll, kTRUE);
  TH1F* hTRDEffPtNegAll = EfficiencyTRD(h3TPCrefNeg, h3TRDrefNegAll, kTRUE);
  TH1F* hTRDEffPtPosTrk4 = EfficiencyTRD(h3TPCrefPos, h3TRDrefPosTrk4, kTRUE);
  TH1F* hTRDEffPtNegTrk4 = EfficiencyTRD(h3TPCrefNeg, h3TRDrefNegTrk4, kTRUE);
  TH1F* hTRDEffPtPosTrk5 = EfficiencyTRD(h3TPCrefPos, h3TRDrefPosTrk5, kTRUE);
  TH1F* hTRDEffPtNegTrk5 = EfficiencyTRD(h3TPCrefNeg, h3TRDrefNegTrk5, kTRUE);
  TH1F* hTRDEffPtPosTrk6 = EfficiencyTRD(h3TPCrefPos, h3TRDrefPosTrk6, kTRUE);
  TH1F* hTRDEffPtNegTrk6 = EfficiencyTRD(h3TPCrefNeg, h3TRDrefNegTrk6, kTRUE);
  
  TH1F* hTOFEffPtPosAll = EfficiencyTRD(h3TRDrefPosAll, h3TOFrefPosAll, kFALSE);
  TH1F* hTOFEffPtNegAll = EfficiencyTRD(h3TRDrefNegAll, h3TOFrefNegAll, kFALSE);
  TH1F* hTOFEffPtPosTrk4 = EfficiencyTRD(h3TRDrefPosTrk4, h3TOFrefPosTrk4, kFALSE);
  TH1F* hTOFEffPtNegTrk4 = EfficiencyTRD(h3TRDrefNegTrk4, h3TOFrefNegTrk4, kFALSE);
  TH1F* hTOFEffPtPosTrk5 = EfficiencyTRD(h3TRDrefPosTrk5, h3TOFrefPosTrk5, kFALSE);
  TH1F* hTOFEffPtNegTrk5 = EfficiencyTRD(h3TRDrefNegTrk5, h3TOFrefNegTrk5, kFALSE);
  TH1F* hTOFEffPtPosTrk6 = EfficiencyTRD(h3TRDrefPosTrk6, h3TOFrefPosTrk6, kFALSE);
  TH1F* hTOFEffPtNegTrk6 = EfficiencyTRD(h3TRDrefNegTrk6, h3TOFrefNegTrk6, kFALSE);
  
  delete h3TPCrefPos; delete h3TPCrefNeg;
  delete h3TRDrefPosAll; delete h3TRDrefPosTrk4; delete h3TRDrefPosTrk5; delete h3TRDrefPosTrk6;
  delete h3TRDrefNegAll; delete h3TRDrefNegTrk4; delete h3TRDrefNegTrk5; delete h3TRDrefNegTrk6;
  delete h3TOFrefPosAll; delete h3TOFrefPosTrk4; delete h3TOFrefPosTrk5; delete h3TOFrefPosTrk6;
  delete h3TOFrefNegAll; delete h3TOFrefNegTrk4; delete h3TOFrefNegTrk5; delete h3TOFrefNegTrk6;
  
  TF1* funcConst = new TF1("constFunc", "[0]", 1.0, 3.0);
  if(trendValues && hTRDEffPtPosAll){ 
    if(hTRDEffPtPosAll->Integral()>0.1) {
      hTRDEffPtPosAll->Fit(funcConst, "Q0ME", "goff", 1.0, 3.0);
      trendValues[0] = funcConst->GetParameter(0);
      trendValues[1] = funcConst->GetParError(0);
    }
  }
  if(trendValues && hTRDEffPtNegAll){ 
    if(hTRDEffPtNegAll->Integral()>0.1) {
      hTRDEffPtNegAll->Fit(funcConst, "Q0ME", "goff", 1.0, 3.0);
      trendValues[2] = funcConst->GetParameter(0);
      trendValues[3] = funcConst->GetParError(0);
    }
  }
  if(trendValues && hTOFEffPtPosAll){ 
    if(hTOFEffPtPosAll->Integral()>0.1) {
      hTOFEffPtPosAll->Fit(funcConst, "Q0ME", "goff", 1.0, 3.0);
      trendValues[4] = funcConst->GetParameter(0);
      trendValues[5] = funcConst->GetParError(0);
    }
  }
  if(trendValues && hTOFEffPtNegAll){ 
    if(hTOFEffPtNegAll->Integral()>0.1) {
      hTOFEffPtNegAll->Fit(funcConst, "Q0ME", "goff", 1.0, 3.0);
      trendValues[6] = funcConst->GetParameter(0);
      trendValues[7] = funcConst->GetParError(0);
    }
  }  
  // get matching efficiencies for isolated bunches
  TH3D* h3TPCrefPos_IsolatedBC=0x0; TH3D* h3TPCrefNeg_IsolatedBC=0x0;
  TH3D* h3TRDrefPosAll_IsolatedBC=0x0; TH3D* h3TRDrefPosTrk4_IsolatedBC=0x0; TH3D* h3TRDrefPosTrk5_IsolatedBC=0x0; TH3D* h3TRDrefPosTrk6_IsolatedBC=0x0;
  TH3D* h3TRDrefNegAll_IsolatedBC=0x0; TH3D* h3TRDrefNegTrk4_IsolatedBC=0x0; TH3D* h3TRDrefNegTrk5_IsolatedBC=0x0; TH3D* h3TRDrefNegTrk6_IsolatedBC=0x0;
  TH3D* h3TOFrefPosAll_IsolatedBC=0x0; TH3D* h3TOFrefPosTrk4_IsolatedBC=0x0; TH3D* h3TOFrefPosTrk5_IsolatedBC=0x0; TH3D* h3TOFrefPosTrk6_IsolatedBC=0x0;
  TH3D* h3TOFrefNegAll_IsolatedBC=0x0; TH3D* h3TOFrefNegTrk4_IsolatedBC=0x0; TH3D* h3TOFrefNegTrk5_IsolatedBC=0x0; TH3D* h3TOFrefNegTrk6_IsolatedBC=0x0;
  for(Int_t ibc=0; ibc<3500; ++ibc) {
    if(!isIsolated[ibc]) continue;
    fCfContainer->SetRangeUser(kEventBC, Double_t(ibc), Double_t(ibc));
    
    TH3D* temp;
    fCfContainer->SetRangeUser(kTrackCharge, +1.0, +1.0);    // positive charges
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);
    temp = (TH3D*)fCfContainer->Project(0, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TPCrefPos_IsolatedBC) h3TPCrefPos_IsolatedBC = temp;
    else {h3TPCrefPos_IsolatedBC->Add(temp); delete temp;}
    temp = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TRDrefPosAll_IsolatedBC) h3TRDrefPosAll_IsolatedBC = temp;
    else {h3TRDrefPosAll_IsolatedBC->Add(temp); delete temp;}
    temp = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TOFrefPosAll_IsolatedBC) h3TOFrefPosAll_IsolatedBC = temp;
    else {h3TOFrefPosAll_IsolatedBC->Add(temp); delete temp;}
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 4.0, 4.0);
    temp = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TRDrefPosTrk4_IsolatedBC) h3TRDrefPosTrk4_IsolatedBC = temp;
    else {h3TRDrefPosTrk4_IsolatedBC->Add(temp); delete temp;}
    temp = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TOFrefPosTrk4_IsolatedBC) h3TOFrefPosTrk4_IsolatedBC = temp;
    else {h3TOFrefPosTrk4_IsolatedBC->Add(temp); delete temp;}
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 5.0, 5.0);
    temp = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TRDrefPosTrk5_IsolatedBC) h3TRDrefPosTrk5_IsolatedBC = temp;
    else {h3TRDrefPosTrk5_IsolatedBC->Add(temp); delete temp;}
    temp = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TOFrefPosTrk5_IsolatedBC) h3TOFrefPosTrk5_IsolatedBC = temp;
    else {h3TOFrefPosTrk5_IsolatedBC->Add(temp); delete temp;}
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 6.0, 6.0);
    temp = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TRDrefPosTrk6_IsolatedBC) h3TRDrefPosTrk6_IsolatedBC = temp;
    else {h3TRDrefPosTrk6_IsolatedBC->Add(temp); delete temp;}
    temp = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TOFrefPosTrk6_IsolatedBC) h3TOFrefPosTrk6_IsolatedBC = temp;
    else {h3TOFrefPosTrk6_IsolatedBC->Add(temp); delete temp;}  
    
    fCfContainer->SetRangeUser(kTrackCharge, -1.0, -1.0);   // negative charges
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);
    temp = (TH3D*)fCfContainer->Project(0, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TPCrefNeg_IsolatedBC) h3TPCrefNeg_IsolatedBC = temp;
    else {h3TPCrefNeg_IsolatedBC->Add(temp); delete temp;}
    temp = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TRDrefNegAll_IsolatedBC) h3TRDrefNegAll_IsolatedBC = temp;
    else {h3TRDrefNegAll_IsolatedBC->Add(temp); delete temp;}
    temp = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TOFrefNegAll_IsolatedBC) h3TOFrefNegAll_IsolatedBC = temp;
    else {h3TOFrefNegAll_IsolatedBC->Add(temp); delete temp;}
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 4.0, 4.0);
    temp = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TRDrefNegTrk4_IsolatedBC) h3TRDrefNegTrk4_IsolatedBC = temp;
    else {h3TRDrefNegTrk4_IsolatedBC->Add(temp); delete temp;}
    temp = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TOFrefNegTrk4_IsolatedBC) h3TOFrefNegTrk4_IsolatedBC = temp;
    else {h3TOFrefNegTrk4_IsolatedBC->Add(temp); delete temp;}
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 5.0, 5.0);
    temp = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TRDrefNegTrk5_IsolatedBC) h3TRDrefNegTrk5_IsolatedBC = temp;
    else {h3TRDrefNegTrk5_IsolatedBC->Add(temp); delete temp;}
    temp = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TOFrefNegTrk5_IsolatedBC) h3TOFrefNegTrk5_IsolatedBC = temp;
    else {h3TOFrefNegTrk5_IsolatedBC->Add(temp); delete temp;}
    fCfContainer->SetRangeUser(kTrackTrdTracklets, 6.0, 6.0);
    temp = (TH3D*)fCfContainer->Project(1, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TRDrefNegTrk6_IsolatedBC) h3TRDrefNegTrk6_IsolatedBC = temp;
    else {h3TRDrefNegTrk6_IsolatedBC->Add(temp); delete temp;}
    temp = (TH3D*)fCfContainer->Project(2, kTrackEta, kTrackPhi, kTrackPt);
    if(!h3TOFrefNegTrk6_IsolatedBC) h3TOFrefNegTrk6_IsolatedBC = temp;
    else {h3TOFrefNegTrk6_IsolatedBC->Add(temp); delete temp;}
  }
  fCfContainer->SetRangeUser(kEventBC, 0.0, 3500.0);
  fCfContainer->SetRangeUser(kTrackCharge, -1.0, +1.0);
  fCfContainer->SetRangeUser(kTrackTrdTracklets, 0.0, 6.0);
    
  TH1F* hTRDEffPtPosAll_IsolatedBC  = EfficiencyTRD(h3TPCrefPos_IsolatedBC, h3TRDrefPosAll_IsolatedBC, kTRUE);
  TH1F* hTRDEffPtNegAll_IsolatedBC  = EfficiencyTRD(h3TPCrefNeg_IsolatedBC, h3TRDrefNegAll_IsolatedBC, kTRUE);
  TH1F* hTRDEffPtPosTrk4_IsolatedBC = EfficiencyTRD(h3TPCrefPos_IsolatedBC, h3TRDrefPosTrk4_IsolatedBC, kTRUE);
  TH1F* hTRDEffPtNegTrk4_IsolatedBC = EfficiencyTRD(h3TPCrefNeg_IsolatedBC, h3TRDrefNegTrk4_IsolatedBC, kTRUE);
  TH1F* hTRDEffPtPosTrk5_IsolatedBC = EfficiencyTRD(h3TPCrefPos_IsolatedBC, h3TRDrefPosTrk5_IsolatedBC, kTRUE);
  TH1F* hTRDEffPtNegTrk5_IsolatedBC = EfficiencyTRD(h3TPCrefNeg_IsolatedBC, h3TRDrefNegTrk5_IsolatedBC, kTRUE);
  TH1F* hTRDEffPtPosTrk6_IsolatedBC = EfficiencyTRD(h3TPCrefPos_IsolatedBC, h3TRDrefPosTrk6_IsolatedBC, kTRUE);
  TH1F* hTRDEffPtNegTrk6_IsolatedBC = EfficiencyTRD(h3TPCrefNeg_IsolatedBC, h3TRDrefNegTrk6_IsolatedBC, kTRUE);
  
  TH1F* hTOFEffPtPosAll_IsolatedBC  = EfficiencyTRD(h3TRDrefPosAll_IsolatedBC, h3TOFrefPosAll_IsolatedBC, kFALSE);
  TH1F* hTOFEffPtNegAll_IsolatedBC  = EfficiencyTRD(h3TRDrefNegAll_IsolatedBC, h3TOFrefNegAll_IsolatedBC, kFALSE);
  TH1F* hTOFEffPtPosTrk4_IsolatedBC = EfficiencyTRD(h3TRDrefPosTrk4_IsolatedBC, h3TOFrefPosTrk4_IsolatedBC, kFALSE);
  TH1F* hTOFEffPtNegTrk4_IsolatedBC = EfficiencyTRD(h3TRDrefNegTrk4_IsolatedBC, h3TOFrefNegTrk4_IsolatedBC, kFALSE);
  TH1F* hTOFEffPtPosTrk5_IsolatedBC = EfficiencyTRD(h3TRDrefPosTrk5_IsolatedBC, h3TOFrefPosTrk5_IsolatedBC, kFALSE);
  TH1F* hTOFEffPtNegTrk5_IsolatedBC = EfficiencyTRD(h3TRDrefNegTrk5_IsolatedBC, h3TOFrefNegTrk5_IsolatedBC, kFALSE);
  TH1F* hTOFEffPtPosTrk6_IsolatedBC = EfficiencyTRD(h3TRDrefPosTrk6_IsolatedBC, h3TOFrefPosTrk6_IsolatedBC, kFALSE);
  TH1F* hTOFEffPtNegTrk6_IsolatedBC = EfficiencyTRD(h3TRDrefNegTrk6_IsolatedBC, h3TOFrefNegTrk6_IsolatedBC, kFALSE);
  
  if(h3TPCrefPos_IsolatedBC) delete h3TPCrefPos_IsolatedBC; if(h3TPCrefNeg_IsolatedBC) delete h3TPCrefNeg_IsolatedBC;
  if(h3TRDrefPosAll_IsolatedBC) delete h3TRDrefPosAll_IsolatedBC; if(h3TRDrefPosTrk4_IsolatedBC) delete h3TRDrefPosTrk4_IsolatedBC; 
  if(h3TRDrefPosTrk5_IsolatedBC) delete h3TRDrefPosTrk5_IsolatedBC; if(h3TRDrefPosTrk6_IsolatedBC) delete h3TRDrefPosTrk6_IsolatedBC;
  if(h3TRDrefNegAll_IsolatedBC) delete h3TRDrefNegAll_IsolatedBC; if(h3TRDrefNegTrk4_IsolatedBC) delete h3TRDrefNegTrk4_IsolatedBC; 
  if(h3TRDrefNegTrk5_IsolatedBC) delete h3TRDrefNegTrk5_IsolatedBC; if(h3TRDrefNegTrk6_IsolatedBC) delete h3TRDrefNegTrk6_IsolatedBC;
  if(h3TOFrefPosAll_IsolatedBC) delete h3TOFrefPosAll_IsolatedBC; if(h3TOFrefPosTrk4_IsolatedBC) delete h3TOFrefPosTrk4_IsolatedBC; 
  if(h3TOFrefPosTrk5_IsolatedBC) delete h3TOFrefPosTrk5_IsolatedBC; if(h3TOFrefPosTrk6_IsolatedBC) delete h3TOFrefPosTrk6_IsolatedBC;
  if(h3TOFrefNegAll_IsolatedBC) delete h3TOFrefNegAll_IsolatedBC; if(h3TOFrefNegTrk4_IsolatedBC) delete h3TOFrefNegTrk4_IsolatedBC; 
  if(h3TOFrefNegTrk5_IsolatedBC) delete h3TOFrefNegTrk5_IsolatedBC; if(h3TOFrefNegTrk6_IsolatedBC) delete h3TOFrefNegTrk6_IsolatedBC;
  
  if(trendValues &&  hTRDEffPtPosAll_IsolatedBC && hTRDEffPtPosAll_IsolatedBC->Integral()>0.1) {
    hTRDEffPtPosAll_IsolatedBC->Fit(funcConst, "Q0ME", "goff", 1.0, 3.0);
    trendValues[18] = funcConst->GetParameter(0);
    trendValues[19] = funcConst->GetParError(0);
  }
  if(trendValues &&  hTRDEffPtNegAll_IsolatedBC && hTRDEffPtNegAll_IsolatedBC->Integral()>0.1) {
    hTRDEffPtNegAll_IsolatedBC->Fit(funcConst, "Q0ME", "goff", 1.0, 3.0);
    trendValues[20] = funcConst->GetParameter(0);
    trendValues[21] = funcConst->GetParError(0);
  }
  if(trendValues &&  hTOFEffPtPosAll_IsolatedBC && hTOFEffPtPosAll_IsolatedBC->Integral()>0.1) {
    hTOFEffPtPosAll_IsolatedBC->Fit(funcConst, "Q0ME", "goff", 1.0, 3.0);
    trendValues[22] = funcConst->GetParameter(0);
    trendValues[23] = funcConst->GetParError(0);
  }
  if(trendValues &&  hTOFEffPtNegAll_IsolatedBC && hTOFEffPtNegAll_IsolatedBC->Integral()>0.1) {
    hTOFEffPtNegAll_IsolatedBC->Fit(funcConst, "Q0ME", "goff", 1.0, 3.0);
    trendValues[24] = funcConst->GetParameter(0);
    trendValues[25] = funcConst->GetParError(0);
  }
  
  //---------------------------------------------------------
  // TPC-TRD matching efficiency vs pt
  pad = ((TVirtualPad*)l->At(6)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  
  if(gROOT->FindObject("rangeEffPt2")) delete gROOT->FindObject("rangeEffPt2");
  TH2F* rangeEffPt=new TH2F("rangeEffPt2", "",10,0.,10.,10,0.,1.4);
  rangeEffPt->SetStats(kFALSE);
  SetStyle(rangeEffPt->GetXaxis(), "p_{T} [GeV/c]", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeEffPt->GetYaxis(), "efficiency", 0.07, 0.8, kTRUE, 0.05);
  rangeEffPt->Draw();
  lat->DrawLatex(0.2, 1.42, "TPC-TRD matching efficiency");
  //++++++++++++++++++
  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.DrawLine(rangeEffPt->GetXaxis()->GetXmin(), 0.7, rangeEffPt->GetXaxis()->GetXmax(), 0.7);
  line.DrawLine(rangeEffPt->GetXaxis()->GetXmin(), 0.9, rangeEffPt->GetXaxis()->GetXmax(), 0.9);
  TLegend* leg=new TLegend(0.2, 0.7, 0.7, 0.89);
  leg->SetNColumns(2);
  leg->SetMargin(0.15);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  
  if(!useIsolatedBC) {
    SetStyle(hTRDEffPtPosAll, 1, kRed, 1, 24, kRed, 1);
    SetStyle(hTRDEffPtNegAll, 1, kBlue, 1, 24, kBlue, 1);
    SetStyle(hTRDEffPtPosTrk4, 1, kRed, 1, 25, kRed, 1);
    SetStyle(hTRDEffPtNegTrk4, 1, kBlue, 1, 25, kBlue, 1);
    SetStyle(hTRDEffPtPosTrk5, 1, kRed, 1, 26, kRed, 1);
    SetStyle(hTRDEffPtNegTrk5, 1, kBlue, 1, 26, kBlue, 1);
    SetStyle(hTRDEffPtPosTrk6, 1, kRed, 1, 27, kRed, 1);
    SetStyle(hTRDEffPtNegTrk6, 1, kBlue, 1, 27, kBlue, 1);
    leg->SetHeader("All bunch crossings");
    hTRDEffPtPosAll->Draw("same"); leg->AddEntry(hTRDEffPtPosAll, "pos. (#geq 1 tracklet)", "p");
    hTRDEffPtNegAll->Draw("same"); leg->AddEntry(hTRDEffPtNegAll, "neg. (#geq 1 tracklet)", "p");
    hTRDEffPtPosTrk4->Draw("same"); leg->AddEntry(hTRDEffPtPosTrk4, "pos. (4 tracklets)", "p");
    hTRDEffPtNegTrk4->Draw("same"); leg->AddEntry(hTRDEffPtNegTrk4, "neg. (4 tracklets)", "p");
    hTRDEffPtPosTrk5->Draw("same"); leg->AddEntry(hTRDEffPtPosTrk5, "pos. (5 tracklets)", "p");
    hTRDEffPtNegTrk5->Draw("same"); leg->AddEntry(hTRDEffPtNegTrk5, "neg. (5 tracklets)", "p");
    hTRDEffPtPosTrk6->Draw("same"); leg->AddEntry(hTRDEffPtPosTrk6, "pos. (6 tracklets)", "p");
    hTRDEffPtNegTrk6->Draw("same"); leg->AddEntry(hTRDEffPtNegTrk6, "neg. (6 tracklets)", "p");
  }
  else {
    if(nIsolatedBC>0) {
      SetStyle(hTRDEffPtPosAll_IsolatedBC, 1, kRed, 1, 24, kRed, 1);
      SetStyle(hTRDEffPtNegAll_IsolatedBC, 1, kBlue, 1, 24, kBlue, 1);
      SetStyle(hTRDEffPtPosTrk4_IsolatedBC, 1, kRed, 1, 25, kRed, 1);
      SetStyle(hTRDEffPtNegTrk4_IsolatedBC, 1, kBlue, 1, 25, kBlue, 1);
      SetStyle(hTRDEffPtPosTrk5_IsolatedBC, 1, kRed, 1, 26, kRed, 1);
      SetStyle(hTRDEffPtNegTrk5_IsolatedBC, 1, kBlue, 1, 26, kBlue, 1);
      SetStyle(hTRDEffPtPosTrk6_IsolatedBC, 1, kRed, 1, 27, kRed, 1);
      SetStyle(hTRDEffPtNegTrk6_IsolatedBC, 1, kBlue, 1, 27, kBlue, 1);
      leg->SetHeader("Isolated bunch crossings");
      hTRDEffPtPosAll_IsolatedBC->Draw("same"); leg->AddEntry(hTRDEffPtPosAll_IsolatedBC, "pos. (#geq 1 tracklet)", "p");
      hTRDEffPtNegAll_IsolatedBC->Draw("same"); leg->AddEntry(hTRDEffPtNegAll_IsolatedBC, "neg. (#geq 1 tracklet)", "p");
      hTRDEffPtPosTrk4_IsolatedBC->Draw("same"); leg->AddEntry(hTRDEffPtPosTrk4_IsolatedBC, "pos. (4 tracklets)", "p");
      hTRDEffPtNegTrk4_IsolatedBC->Draw("same"); leg->AddEntry(hTRDEffPtNegTrk4_IsolatedBC, "neg. (4 tracklets)", "p");
      hTRDEffPtPosTrk5_IsolatedBC->Draw("same"); leg->AddEntry(hTRDEffPtPosTrk5_IsolatedBC, "pos. (5 tracklets)", "p");
      hTRDEffPtNegTrk5_IsolatedBC->Draw("same"); leg->AddEntry(hTRDEffPtNegTrk5_IsolatedBC, "neg. (5 tracklets)", "p");
      hTRDEffPtPosTrk6_IsolatedBC->Draw("same"); leg->AddEntry(hTRDEffPtPosTrk6_IsolatedBC, "pos. (6 tracklets)", "p");
      hTRDEffPtNegTrk6_IsolatedBC->Draw("same"); leg->AddEntry(hTRDEffPtNegTrk6_IsolatedBC, "neg. (6 tracklets)", "p");
    }
  }
  
  leg->Draw();
  
  //---------------------------------------------------------
  // TRD-TOF matching efficiency vs pt
  pad = ((TVirtualPad*)l->At(7)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  
  rangeEffPt->Draw();
  lat->DrawLatex(0.2, 1.42, "TRD-TOF matching efficiency");
  
  if(!useIsolatedBC) {
    SetStyle(hTOFEffPtPosAll, 1, kRed, 1, 24, kRed, 1);
    SetStyle(hTOFEffPtPosTrk4, 1, kRed, 1, 25, kRed, 1);
    SetStyle(hTOFEffPtPosTrk5, 1, kRed, 1, 26, kRed, 1);
    SetStyle(hTOFEffPtPosTrk6, 1, kRed, 1, 27, kRed, 1);
    SetStyle(hTOFEffPtNegAll, 1, kBlue, 1, 24, kBlue, 1);
    SetStyle(hTOFEffPtNegTrk4, 1, kBlue, 1, 25, kBlue, 1);
    SetStyle(hTOFEffPtNegTrk5, 1, kBlue, 1, 26, kBlue, 1);
    SetStyle(hTOFEffPtNegTrk6, 1, kBlue, 1, 27, kBlue, 1);
    hTOFEffPtPosAll->Draw("same"); 
    hTOFEffPtPosTrk4->Draw("same"); 
    hTOFEffPtPosTrk5->Draw("same"); 
    hTOFEffPtPosTrk6->Draw("same"); 
    hTOFEffPtNegAll->Draw("same"); 
    hTOFEffPtNegTrk4->Draw("same"); 
    hTOFEffPtNegTrk5->Draw("same"); 
    hTOFEffPtNegTrk6->Draw("same"); 
  }
  else {
    if(nIsolatedBC>0) {
      SetStyle(hTOFEffPtPosAll_IsolatedBC, 1, kRed, 1, 24, kRed, 1);
      SetStyle(hTOFEffPtPosTrk4_IsolatedBC, 1, kRed, 1, 25, kRed, 1);
      SetStyle(hTOFEffPtPosTrk5_IsolatedBC, 1, kRed, 1, 26, kRed, 1);
      SetStyle(hTOFEffPtPosTrk6_IsolatedBC, 1, kRed, 1, 27, kRed, 1);
      SetStyle(hTOFEffPtNegAll_IsolatedBC, 1, kBlue, 1, 24, kBlue, 1);
      SetStyle(hTOFEffPtNegTrk4_IsolatedBC, 1, kBlue, 1, 25, kBlue, 1);
      SetStyle(hTOFEffPtNegTrk5_IsolatedBC, 1, kBlue, 1, 26, kBlue, 1);
      SetStyle(hTOFEffPtNegTrk6_IsolatedBC, 1, kBlue, 1, 27, kBlue, 1);
      hTOFEffPtPosAll_IsolatedBC->Draw("same"); 
      hTOFEffPtPosTrk4_IsolatedBC->Draw("same"); 
      hTOFEffPtPosTrk5_IsolatedBC->Draw("same"); 
      hTOFEffPtPosTrk6_IsolatedBC->Draw("same"); 
      hTOFEffPtNegAll_IsolatedBC->Draw("same"); 
      hTOFEffPtNegTrk4_IsolatedBC->Draw("same"); 
      hTOFEffPtNegTrk5_IsolatedBC->Draw("same"); 
      hTOFEffPtNegTrk6_IsolatedBC->Draw("same"); 
    }
  } 
    
  //-----------------------------------------------------
  // <ntracklets> vs (phi,eta)
  pad = ((TVirtualPad*)l->At(2)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  
  rangeEtaPhi->Draw();
  lat->DrawLatex(-0.9, 3.6, "TRD <N_{tracklets}>");
  
  TH3D* hNtracklets = (TH3D*)fCfContainer->Project(1, kTrackPhi, kTrackEta, kTrackTrdTracklets);
  TH3D* hNtracklets_IsolatedBC=0x0;
  for(Int_t ibc=0; ibc<3500; ++ibc) {
    if(!isIsolated[ibc]) continue;
    fCfContainer->SetRangeUser(kEventBC, Double_t(ibc), Double_t(ibc));
    TH3D* temp = (TH3D*)fCfContainer->Project(1, kTrackPhi, kTrackEta, kTrackTrdTracklets);
    if(!hNtracklets_IsolatedBC) hNtracklets_IsolatedBC = temp;
    else {hNtracklets_IsolatedBC->Add(temp); delete temp;}
  }  
  fCfContainer->SetRangeUser(kEventBC, 0.0, 3500.0);   // reset the BC range
  
  TProfile2D* hNtrackletsProf=0x0;
  if(useIsolatedBC) { 
    if(nIsolatedBC>0)
      hNtrackletsProf = hNtracklets_IsolatedBC->Project3DProfile();
  }
  else hNtrackletsProf = hNtracklets->Project3DProfile();
  delete hNtracklets; if(hNtracklets_IsolatedBC) delete hNtracklets_IsolatedBC;
  if(hNtrackletsProf) {
    hNtrackletsProf->SetStats(kFALSE);
    hNtrackletsProf->SetMinimum(0.);
    hNtrackletsProf->SetMaximum(6.);
    hNtrackletsProf->Draw("samecolz");
    DrawTRDGrid();
  }
  
  // calculate the trend value for tracklets/track
  TH2D* hNtrackletsVsP = (TH2D*)fCfContainer->Project(1, kTrackP, kTrackTrdTracklets);
  if(trendValues &&  hNtrackletsVsP && hNtrackletsVsP->GetEntries()>0.1) {
    TProfile* hNtrackletsVsPprof = hNtrackletsVsP->ProfileX("hNtrackletsVsPprof");
    hNtrackletsVsPprof->Fit(funcConst, "QME0", "goff", 1.0, 3.0);
    trendValues[8] = funcConst->GetParameter(0);
    trendValues[9] = funcConst->GetParError(0);
    delete hNtrackletsVsP;
  }
  TH2D* hNtrackletsVsP_IsolatedBC=0x0;
  for(Int_t ibc=0; ibc<3500; ++ibc) {
    if(!isIsolated[ibc]) continue;
    fCfContainer->SetRangeUser(kEventBC, Double_t(ibc), Double_t(ibc));
    TH2D* temp = (TH2D*)fCfContainer->Project(1, kTrackP, kTrackTrdTracklets);
    if(!hNtrackletsVsP_IsolatedBC) hNtrackletsVsP_IsolatedBC = temp;
    else {hNtrackletsVsP_IsolatedBC->Add(temp); delete temp;}
  }
  fCfContainer->SetRangeUser(kEventBC, 0.0, 3500.0);   // reset the BC range
  if(trendValues &&  hNtrackletsVsP_IsolatedBC && hNtrackletsVsP_IsolatedBC->GetEntries()>0.1) {
    TProfile* hNtrackletsVsPprof_IsolatedBC = hNtrackletsVsP_IsolatedBC->ProfileX("hNtrackletsVsPprof_IsolatedBC");
    hNtrackletsVsPprof_IsolatedBC->Fit(funcConst, "QME0", "goff", 1.0, 3.0);
    trendValues[26] = funcConst->GetParameter(0);
    trendValues[27] = funcConst->GetParError(0);
    delete hNtrackletsVsPprof_IsolatedBC;
  }
      
  //--------------------------------------------------------------
  // Nclusters per TRD track vs momentum
  pad = ((TVirtualPad*)l->At(5)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.12);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  pad->SetLogz();
  
  if(gROOT->FindObject("rangeNclsP")) delete gROOT->FindObject("rangeNclsP");
  TH2F* rangeNclsP = new TH2F("rangeNclsP", "", 10, 0.0, 11.99, 10, 0.0, 199.0);
  SetStyle(rangeNclsP->GetXaxis(), "p [GeV/c]", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeNclsP->GetYaxis(), "#clusters", 0.07, 0.8, kTRUE, 0.05);
  rangeNclsP->SetStats(kFALSE);
  rangeNclsP->Draw();
  lat->DrawLatex(1.0, 205., "TRD Clusters / track");
  
  TH2D* hNclsVsP = (TH2D*)fCfContainer->Project(1, kTrackP, kTrackTrdClusters);
  TH2D* hNclsVsP_IsolatedBC=0x0;
  for(Int_t ibc=0; ibc<3500; ++ibc) {
    if(!isIsolated[ibc]) continue;
    fCfContainer->SetRangeUser(kEventBC, Double_t(ibc), Double_t(ibc));
    TH2D* temp = (TH2D*)fCfContainer->Project(1, kTrackP, kTrackTrdClusters);
    if(!hNclsVsP_IsolatedBC) hNclsVsP_IsolatedBC = temp;
    else {hNclsVsP_IsolatedBC->Add(temp); delete temp;}
  }
  fCfContainer->SetRangeUser(kEventBC, 0.0, 3500.0);   // reset the BC range
  
  if(useIsolatedBC && hNclsVsP_IsolatedBC) {
    hNclsVsP_IsolatedBC->SetStats(kFALSE);
    hNclsVsP_IsolatedBC->Draw("samecolz");
  }
  if(!useIsolatedBC && hNclsVsP) {
    hNclsVsP->SetStats(kFALSE);
    hNclsVsP->Draw("samecolz");
  }
    
  if(trendValues && hNclsVsP && hNclsVsP->GetEntries()>10) {
    TProfile* hNclsVsPprof = hNclsVsP->ProfileX("hNclsVsPprof");
    hNclsVsPprof->Fit(funcConst, "QME0", "goff", 1.0, 3.0);
    trendValues[10] = funcConst->GetParameter(0);
    trendValues[11] = funcConst->GetParError(0);
  }
  if(trendValues && hNclsVsP_IsolatedBC && hNclsVsP_IsolatedBC->GetEntries()>10) {
    TProfile* hNclsVsPprof_IsolatedBC = hNclsVsP_IsolatedBC->ProfileX("hNclsVsPprof_IsolateBC");
    hNclsVsPprof_IsolatedBC->Fit(funcConst, "QME0", "goff", 1.0, 3.0);
    trendValues[28] = funcConst->GetParameter(0);
    trendValues[29] = funcConst->GetParError(0);
  }
    
  //--------------------------------------------------------------
  // TRD-TPC and TOF-TRD matching efficiency vs bunch crossing
  pad = ((TVirtualPad*)l->At(8)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);

  TH1D* phiProj = (TH1D*)fCfContainer->Project(1, kTrackPhi);
  Double_t smPhiLimits[19];
  Bool_t activeSM[18] = {kFALSE};
  for(Int_t ism=0; ism<=18; ++ism) smPhiLimits[ism] = -TMath::Pi() + (2.0*TMath::Pi()/18.0)*ism;
  CheckActiveSM(phiProj, activeSM);
  for(Int_t ism=0; ism<18; ++ism) cout << "sm " << ism << " is active : " << (activeSM[ism] ? "yes" : "no") << endl;
  
  fCfContainer->SetRangeUser(kTrackPt, 1.01, 2.99);   // 1.0 < pt < 3.0 GeV/c
  TH2D* hTPCPhiBC = (TH2D*)fCfContainer->Project(0, kEventBC, kTrackPhi);
  TH2D* hTRDPhiBC = (TH2D*)fCfContainer->Project(1, kEventBC, kTrackPhi);
  TH2D* hTOFPhiBC = (TH2D*)fCfContainer->Project(2, kEventBC, kTrackPhi);
  TH1D* projectionBC = (TH1D*)fCfContainer->Project(0, kEventBC);
  fCfContainer->SetRangeUser(kTrackPt, 0.0, 100.0);   // reset the pt range
  TH1D* hTRDEffBC = new TH1D("hTRDEffBC", "", hTPCPhiBC->GetXaxis()->GetNbins(), hTPCPhiBC->GetXaxis()->GetXmin(), hTPCPhiBC->GetXaxis()->GetXmax());
  TH1D* hTOFEffBC = new TH1D("hTOFEffBC", "", hTPCPhiBC->GetXaxis()->GetNbins(), hTPCPhiBC->GetXaxis()->GetXmin(), hTPCPhiBC->GetXaxis()->GetXmax());
  
  for(Int_t bcBin=1; bcBin<=hTPCPhiBC->GetXaxis()->GetNbins(); ++bcBin) {
    if(projectionBC->GetBinContent(bcBin)<0.1) continue;
    Double_t tpcEntries = 0.0; Double_t trdEntries = 0.0; Double_t tofEntries = 0.0;
    for(Int_t phiBin=1; phiBin<=hTPCPhiBC->GetYaxis()->GetNbins(); ++phiBin) {
      Double_t phi = hTPCPhiBC->GetYaxis()->GetBinCenter(phiBin);
      for(Int_t ism=0; ism<18; ++ism) {
        if(phi>=smPhiLimits[ism] && phi<smPhiLimits[ism+1] && activeSM[ism]) {
          tpcEntries += hTPCPhiBC->GetBinContent(bcBin, phiBin);
          trdEntries += hTRDPhiBC->GetBinContent(bcBin, phiBin);
          tofEntries += hTOFPhiBC->GetBinContent(bcBin, phiBin);
        }
      }  // end loop over super-modules
    }  // end loop over phi bins
    hTRDEffBC->SetBinContent(bcBin, (tpcEntries>0.01 ? trdEntries/tpcEntries : 0.0));
    if(tpcEntries>0.01 && trdEntries>0.01 && (tpcEntries-trdEntries)>=0.01) 
    hTRDEffBC->SetBinError(bcBin, TMath::Sqrt(trdEntries*(tpcEntries-trdEntries)/tpcEntries/tpcEntries/tpcEntries));
    hTOFEffBC->SetBinContent(bcBin, (trdEntries>0.01 ? tofEntries/trdEntries : 0.0));
    if(trdEntries>0.01 && tofEntries>0.01 && (trdEntries-tofEntries)>=0.01) 
    hTOFEffBC->SetBinError(bcBin, TMath::Sqrt(tofEntries*(trdEntries-tofEntries)/trdEntries/trdEntries/trdEntries));    
  }  // end loop over BC bins
  delete hTPCPhiBC; delete hTRDPhiBC; delete hTOFPhiBC;
  delete projectionBC;
  
  if(gROOT->FindObject("rangeBC")) delete gROOT->FindObject("rangeBC");
  TH2F* rangeBC = new TH2F("rangeBC", "", 10, -0.5, 3499.5, 10, 0.0, 1.4);
  rangeBC->SetStats(kFALSE);
  SetStyle(rangeBC->GetXaxis(), "Bunch crossing", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeBC->GetYaxis(), "efficiency", 0.07, 0.8, kTRUE, 0.05);
  rangeBC->Draw();
  
  TLegend* legBC=new TLegend(0.8, 0.7, 0.95, 0.89);
  legBC->SetBorderSize(0);
  legBC->SetMargin(0.15);
  legBC->SetFillColor(0);
  if(hTRDEffBC) {
    hTRDEffBC->SetStats(kFALSE);
    SetStyle(hTRDEffBC, 1, kRed, 2, 24, kRed, 1); legBC->AddEntry(hTRDEffBC, "TPC-TRD", "p");
    SetStyle(hTOFEffBC, 1, kBlue, 2, 24, kBlue, 1); legBC->AddEntry(hTOFEffBC, "TRD-TOF", "p");
    hTRDEffBC->Draw("same");
    hTOFEffBC->Draw("same");
    legBC->Draw();
    lat->DrawLatex(200., 1.42, "Matching efficiency at 1<p_{T}<3 GeV/c");
  }
    
  // reset the user range on the event multiplicity
  fCfContainer->SetRangeUser(kEventMult, 0.0, 6.0, kTRUE);
  if(cutTOFbc) fCfContainer->SetRangeUser(kTrackTOFBC, -1000.0, +1000.0);  // reset the cut on TOFbc
  
  delete funcConst;
}


//_________________________________________________________________
void AliTRDcheckESD::PlotPidSummaryFromCF(Int_t centralityClass, Double_t* trendValues, Bool_t useIsolatedBC, Bool_t cutTOFbc) {

  if(!fCfContainer) return;
  
  TLatex *lat=new TLatex();
  lat->SetTextSize(0.07);
  lat->SetTextColor(2);
  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001);
  gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
  gPad->Divide(3,3,0.,0.);
  TList* l=gPad->GetListOfPrimitives();
  
  fCfContainer->SetRangeUser(kTrackDCAxy, -0.999, +0.999);
  fCfContainer->SetRangeUser(kTrackDCAz, -3.0, +3.0);
  if(cutTOFbc) fCfContainer->SetRangeUser(kTrackTOFBC, 0.0, 0.0);
  
  // find all the isolated bunch crossings with entries
  Bool_t isIsolated[3500];
  TH1D* tempTH1D = (TH1D*)fCfContainer->Project(0, kEventBC);
  FindIsolatedBCs(tempTH1D, isIsolated); delete tempTH1D;
  Int_t nIsolatedBC = 0;
  for(Int_t ibc=0; ibc<3500; ++ibc) 
    if(isIsolated[ibc]) nIsolatedBC++;
  
  if(centralityClass>0) // select the multiplicity class
    fCfContainer->SetRangeUser(kEventMult, Double_t(centralityClass), Double_t(centralityClass), kTRUE);
  
  if(gROOT->FindObject("rangeEtaPhi2")) delete gROOT->FindObject("rangeEtaPhi2");
  TH2F* rangeEtaPhi = new TH2F("rangeEtaPhi2", "", 10, -0.99, +0.99, 10, -3.4, +3.4);
  SetStyle(rangeEtaPhi->GetXaxis(), "#eta", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeEtaPhi->GetYaxis(), "detector #varphi", 0.07, 0.8, kTRUE, 0.05);
  rangeEtaPhi->SetStats(kFALSE);  
  
  // eta-phi distr. for <Qtot> in layer 0
  TVirtualPad* pad;
  TProfile2D* hProf2D;
  for(Int_t iLayer=0; iLayer<6; ++iLayer) {
    pad = ((TVirtualPad*)l->At((iLayer<3 ? iLayer*3 : (iLayer-3)*3+1))); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    rangeEtaPhi->Draw();
    
    TH3D* hQtotEtaPhi=0x0;
    if(useIsolatedBC) {
      for(Int_t ibc=0; ibc<3500; ++ibc) {
        if(!isIsolated[ibc]) continue;
        fCfContainer->SetRangeUser(kEventBC, Double_t(ibc), Double_t(ibc));
        TH3D* temp = (TH3D*)fCfContainer->Project(1, kTrackPhi, kTrackEta, kTrackQtot+iLayer);
        if(!hQtotEtaPhi) hQtotEtaPhi = temp;
        else {hQtotEtaPhi->Add(temp); delete temp;}
      }
      fCfContainer->SetRangeUser(kEventBC, 0.0, 3500.0);
    }
    else hQtotEtaPhi = (TH3D*)fCfContainer->Project(1, kTrackPhi, kTrackEta, kTrackQtot+iLayer);
    hProf2D = (hQtotEtaPhi ? hQtotEtaPhi->Project3DProfile() : 0x0);
    if(hQtotEtaPhi) delete hQtotEtaPhi;
    
    if(hProf2D) {
      hProf2D->SetStats(kFALSE);
      hProf2D->SetMinimum(0.);
      hProf2D->SetMaximum(4.);
      hProf2D->Draw("samecolz");
    }
    lat->DrawLatex(-0.9, 3.6, Form("TRD <Q_{tot}> Layer %d", iLayer));
    DrawTRDGrid();
  }
    
  // PH versus slice number
  pad = ((TVirtualPad*)l->At(2)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.03); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  
  if(gROOT->FindObject("rangePHslice")) delete gROOT->FindObject("rangePHslice");
  TH2F* rangePHslice=new TH2F("rangePHslice", "", 8, -0.5, 7.5, 10, 0.0, 2000.);
  rangePHslice->SetStats(kFALSE);
  SetStyle(rangePHslice->GetXaxis(), "slice", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangePHslice->GetYaxis(), "PH", 0.07, 0.8, kTRUE, 0.05);
  rangePHslice->Draw();
  
  TF1* funcPol1 = new TF1("funcPol1", "[0]+[1]*x", 2.9, 6.4);
  
  TH2F* h2F;
  TH1D* hF;
  if((h2F = dynamic_cast<TH2F*>(fHistos->At(kPHSlice+centralityClass)))) {
    hF = Proj2D(h2F);
    h2F->SetStats(kFALSE);
    h2F->Draw("samecolz");
    if(trendValues) {
      hF->Fit(funcPol1, "QME0", "goff", 2.9, 6.4);
      trendValues[12] = funcPol1->GetParameter(0);  // PH plateau
      trendValues[13] = funcPol1->GetParError(0);   // PH plateau
      trendValues[14] = funcPol1->GetParameter(1);  // PH slope
      trendValues[15] = funcPol1->GetParError(1);   // PH slope
    }
    hF->SetLineWidth(2);
    hF->SetLineStyle(2);
    hF->Draw("same");
  }
  delete funcPol1;
  
  // Qtot vs P
  pad = ((TVirtualPad*)l->At(5)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.03); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  pad->SetLogz();
  
  if(gROOT->FindObject("rangeQtotP")) delete gROOT->FindObject("rangeQtotP");
  TH2F* rangeQtotP = new TH2F("rangeQtotP", "", 10, 0.0, 11.99, 10, 0.0, 11.99);
  SetStyle(rangeQtotP->GetXaxis(), "P [GeV/c]", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeQtotP->GetYaxis(), "Q_{tot}", 0.07, 0.8, kTRUE, 0.05);
  rangeQtotP->SetStats(kFALSE);
  rangeQtotP->Draw();
  
  TH2D* hQtotP_IsolatedBC = 0x0;
  for(Int_t ibc=0; ibc<3500; ++ibc) {
    if(!isIsolated[ibc]) continue;
    fCfContainer->SetRangeUser(kEventBC, Double_t(ibc), Double_t(ibc));
    TH2D* temp = (TH2D*)fCfContainer->Project(1, kTrackP, kTrackQtot);
    if(!hQtotP_IsolatedBC) hQtotP_IsolatedBC = temp;
    else {hQtotP_IsolatedBC->Add(temp); delete temp;}
    for(Int_t il=1; il<6; ++il) {
      temp = (TH2D*)fCfContainer->Project(1, kTrackP, kTrackQtot+il);
      hQtotP_IsolatedBC->Add(temp); delete temp;
    }
  }
  fCfContainer->SetRangeUser(kEventBC, 0.0, 3500.0);
  
  TH2D* temp = (TH2D*)fCfContainer->Project(1, kTrackP, kTrackQtot);
  TH2D* hQtotP = temp;
  for(Int_t il=1; il<6; ++il) {
    temp = (TH2D*)fCfContainer->Project(1, kTrackP, kTrackQtot+il);
    hQtotP->Add(temp); delete temp;
  }
  
  if(hQtotP_IsolatedBC)
    for(Int_t i=1; i<=hQtotP_IsolatedBC->GetXaxis()->GetNbins(); ++i) 
      hQtotP_IsolatedBC->SetBinContent(i, 1, 0.0);  
  TH1D* hQtotProj_IsolatedBC = (hQtotP_IsolatedBC ? Proj2D(hQtotP_IsolatedBC) : 0x0);
  if(hQtotProj_IsolatedBC) SetStyle(hQtotProj_IsolatedBC, 2, kBlue, 2, 1, kBlue, 1);
  if(trendValues && hQtotProj_IsolatedBC && hQtotProj_IsolatedBC->GetEntries()>2) {
    trendValues[30] = hQtotProj_IsolatedBC->GetBinContent(hQtotProj_IsolatedBC->FindBin(1.0));   // Landau MPV at 1GeV/c
    trendValues[31] = hQtotProj_IsolatedBC->GetBinError(hQtotProj_IsolatedBC->FindBin(1.0));     // Landau width at 1 GeV/c
  }
  if(hQtotP)
    for(Int_t i=1; i<=hQtotP->GetXaxis()->GetNbins(); ++i) 
      hQtotP->SetBinContent(i, 1, 0.0);  
  TH1D* hQtotProj = (hQtotP ? Proj2D(hQtotP) : 0x0);
  if(hQtotProj) SetStyle(hQtotProj, 2, kBlue, 2, 1, kBlue, 1);
  if(trendValues && hQtotProj && hQtotProj->GetEntries()>2) {
    trendValues[16] = hQtotProj->GetBinContent(hQtotProj->FindBin(1.0));   // Landau MPV at 1GeV/c
    trendValues[17] = hQtotProj->GetBinError(hQtotProj->FindBin(1.0));     // Landau width at 1 GeV/c
  }
  if(useIsolatedBC && hQtotP_IsolatedBC) {
    hQtotP_IsolatedBC->SetStats(kFALSE);
    for(Int_t i=1; i<=hQtotP_IsolatedBC->GetXaxis()->GetNbins(); ++i) hQtotP_IsolatedBC->SetBinContent(i, 1, 0.0);  
    hQtotP_IsolatedBC->Draw("samecolz");
    hQtotProj_IsolatedBC->Draw("same");
  }
  if(!useIsolatedBC && hQtotP) {
    hQtotP->SetStats(kFALSE);
    hQtotP->Draw("samecolz");
    hQtotProj->Draw("same");
  }

  // reset the user range on the event multiplicity
  fCfContainer->SetRangeUser(kEventMult, 0.0, 6.0, kTRUE);
  if(cutTOFbc) fCfContainer->SetRangeUser(kTrackTOFBC, -1000.0, +1000.0);  // reset the cut on TOFbc
}


//_________________________________________________________________
Bool_t AliTRDcheckESD::PlotCentSummary(Double_t* trendValues) {

  Bool_t isGoodForSaving=kFALSE;
  
  TLatex* lat=new TLatex();
  lat->SetTextSize(0.06);
  lat->SetTextColor(2);

  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001);
  gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
  gPad->Divide(3,3,0.,0.);
  TList* l=gPad->GetListOfPrimitives();
  
  TPad* pad=0x0;  

  if(gROOT->FindObject("rangeEffPt")) delete gROOT->FindObject("rangeEffPt");
  TH2F* rangeEffPt=new TH2F("rangeEffPt", "",10,0.,10.,10,0.,1.4);
  rangeEffPt->SetStats(kFALSE);
  SetStyle(rangeEffPt->GetXaxis(), "p_{T} [GeV/c]", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeEffPt->GetYaxis(), "efficiency", 0.07, 0.8, kTRUE, 0.05);
  
  TH3F *h3(NULL), *h3p(NULL), *h3n(NULL);
  Int_t padsForEffs[5] = {0,3,6,1,4};
  for(Int_t iCent=1; iCent<6; ++iCent) {
    // TPC-TRD matching efficiencies
    pad = ((TPad*)l->At(padsForEffs[iCent-1])); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02); pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    
    if(!(h3p = dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksPos+iCent)))) continue;
    if(!(h3n = dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksNeg+iCent)))) continue;
    // =============================================
    if(!(h3 = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksPos+iCent)))) continue;
    TH1F* hFeffP = EfficiencyTRD(h3p, h3, kTRUE);
    //
    if(!(h3 = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksNeg+iCent)))) continue;
    TH1F* hFeffN = EfficiencyTRD(h3n, h3, kTRUE);
    // =============================================
    if(!(h3 = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksPos4+iCent)))) continue;
    TH1F* hFeffP4 = EfficiencyTRD(h3p, h3, kTRUE);
    //
    if(!(h3 = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksNeg4+iCent)))) continue;
    TH1F* hFeffN4 = EfficiencyTRD(h3n, h3, kTRUE);
    // =============================================
    if(!(h3 = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksPos5+iCent)))) continue;
    TH1F* hFeffP5 = EfficiencyTRD(h3p, h3, kTRUE);
    //
    if(!(h3 = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksNeg5+iCent)))) continue;
    TH1F* hFeffN5 = EfficiencyTRD(h3n, h3, kTRUE);
    // =============================================
    if(!(h3 = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksPos6+iCent)))) continue;
    TH1F* hFeffP6 = EfficiencyTRD(h3p, h3, kTRUE);
    //
    if(!(h3 = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksNeg6+iCent)))) continue;
    TH1F* hFeffN6 = EfficiencyTRD(h3n, h3, kTRUE);
  
    rangeEffPt->Draw();
    
    TLine line;
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.DrawLine(rangeEffPt->GetXaxis()->GetXmin(), 0.7, rangeEffPt->GetXaxis()->GetXmax(), 0.7);
    line.DrawLine(rangeEffPt->GetXaxis()->GetXmin(), 0.9, rangeEffPt->GetXaxis()->GetXmax(), 0.9);
    line.SetLineStyle(1);
    line.SetLineWidth(1);
    line.DrawLine(rangeEffPt->GetXaxis()->GetXmin(), 1.0, rangeEffPt->GetXaxis()->GetXmax(), 1.0);
    if(hFeffP) SetStyle(hFeffP, 1, kRed, 1, 24, kRed, 1);
    if(hFeffP4) SetStyle(hFeffP4, 1, kRed, 1, 25, kRed, 1);
    if(hFeffP5) SetStyle(hFeffP5, 1, kRed, 1, 26, kRed, 1);
    if(hFeffP6) SetStyle(hFeffP6, 1, kRed, 1, 27, kRed, 1);
    if(hFeffN) SetStyle(hFeffN, 1, kBlue, 1, 24, kBlue, 1);
    if(hFeffN4) SetStyle(hFeffN4, 1, kBlue, 1, 25, kBlue, 1);
    if(hFeffN5) SetStyle(hFeffN5, 1, kBlue, 1, 26, kBlue, 1);
    if(hFeffN6) SetStyle(hFeffN6, 1, kBlue, 1, 27, kBlue, 1);
    
    TLegend* leg=new TLegend(0.16, 0.7, 0.61, 0.89);
    leg->SetFillColor(0);
    leg->SetNColumns(2);
    leg->SetTextSize(0.039);
    leg->SetMargin(0.1);
    if(hFeffP && hFeffP->Integral()>0.001) {
      isGoodForSaving = kTRUE;
      hFeffP->Draw("same");
      leg->AddEntry(hFeffP, "pos. (#geq 1 trcklt)", "p");
    }
    if(hFeffN && hFeffN->Integral()>0.001) {
      isGoodForSaving = kTRUE;
      hFeffN->Draw("same");
      leg->AddEntry(hFeffN, "neg. (#geq 1 trcklt)", "p");
    }
    if(hFeffP4 && hFeffP4->Integral()>0.001) {
      isGoodForSaving = kTRUE;
      hFeffP4->Draw("same");
      leg->AddEntry(hFeffP4, "pos. (4 trcklts)", "p");
    }
    if(hFeffN4 && hFeffN4->Integral()>0.001) {
      isGoodForSaving = kTRUE;
      hFeffN4->Draw("same");
      leg->AddEntry(hFeffN4, "neg. (4 trcklts)", "p");
    }
    if(hFeffP5 && hFeffP5->Integral()>0.001) {
      isGoodForSaving = kTRUE;
      hFeffP5->Draw("same");
      leg->AddEntry(hFeffP5, "pos. (5 trcklts)", "p");
    }
    if(hFeffN5 && hFeffN5->Integral()>0.001) {
      isGoodForSaving = kTRUE;
      hFeffN5->Draw("same");
      leg->AddEntry(hFeffN5, "neg. (5 trcklts)", "p");
    }
    if(hFeffP6 && hFeffP6->Integral()>0.001) {
      isGoodForSaving = kTRUE;
      hFeffP6->Draw("same");
      leg->AddEntry(hFeffP6, "pos. (6 trcklts)", "p");
    }
    if(hFeffN6 && hFeffN6->Integral()>0.001) {
      isGoodForSaving = kTRUE;
      hFeffN6->Draw("same");
      leg->AddEntry(hFeffN6, "neg. (6 trklts)", "p");
    }
    
    if(isGoodForSaving) {
      if(iCent==1) leg->Draw();
      lat->DrawLatex(5.6, 1.3, Form("Centrality class %d", iCent));
      lat->DrawLatex(0.5, 1.42, "TPC-TRD matching efficiency");
    }
  }   // end loop over multiplicity intervals

  // Number of clusters per TRD track
  pad = ((TPad*)l->At(2)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
  pad->SetTopMargin(0.02); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  if(gROOT->FindObject("rangeNcls")) delete gROOT->FindObject("rangeNcls");
  TH2F* rangeNcls = new TH2F("rangeNcls", "", 10, 0.0, 199.9, 10, 0.0, 1.199);
  SetStyle(rangeNcls->GetXaxis(), "# TRD clusters", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeNcls->GetYaxis(), "entries (a.u.)", 0.07, 0.8, kTRUE, 0.05);
  rangeNcls->SetStats(kFALSE);
  rangeNcls->Draw();
  
  TH2F* h2F[6]; TH1D* proj[6];
  TLegend* leg=new TLegend(0.2, 0.7, 0.5, 0.95);
  leg->SetFillColor(0);
  Bool_t isGood=kFALSE;
  for(Int_t iCent=0; iCent<6; ++iCent) {
    h2F[iCent] = dynamic_cast<TH2F*>(fHistos->At(kNClsTrackTRD+iCent));
    proj[iCent] = (h2F[iCent] && h2F[iCent]->GetEntries()>10 ? h2F[iCent]->ProjectionY(Form("projCent%d",iCent)) : 0x0);
    if(proj[iCent]) {
      proj[iCent]->SetLineColor(iCent<4 ? iCent+1 : iCent+2);
      Double_t maximum = proj[iCent]->GetMaximum();
      if(maximum>1.0)
        proj[iCent]->Scale(1.0/maximum);
      proj[iCent]->SetStats(kFALSE);
      proj[iCent]->Draw("same");
      leg->AddEntry(proj[iCent], (iCent==0 ? "all centralities" : Form("centrality class %d", iCent)), "l");
      isGood = kTRUE;
    }
  }
  if(isGood) leg->Draw();
  isGoodForSaving = isGoodForSaving || isGood;
  
  // Qtot vs P
  pad = ((TPad*)l->At(5)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
  pad->SetTopMargin(0.02); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  if(gROOT->FindObject("rangeQtot")) delete gROOT->FindObject("rangeQtot");
  TH2F* rangeQtot = new TH2F("rangeQtot", "", 10, 0.0, 9.999, 10, 0.0, 1.199);
  SetStyle(rangeQtot->GetXaxis(), "Q_{tot} (a.u.)", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeQtot->GetYaxis(), "entries (a.u.)", 0.07, 0.8, kTRUE, 0.05);
  rangeQtot->SetStats(kFALSE);
  rangeQtot->Draw();
  
  TH1D* projQ[6];
  TLegend* leg2=new TLegend(0.6, 0.7, 0.9, 0.95);
  leg2->SetFillColor(0);
  isGood = kFALSE;
  for(Int_t iCent=0; iCent<6; ++iCent) {  
    h2F[iCent] = dynamic_cast<TH2F*>(fHistos->At(kQtotP+iCent));
    projQ[iCent] = (h2F[iCent] && h2F[iCent]->GetEntries()>10 ? h2F[iCent]->ProjectionY(Form("projQCent%d",iCent)) : 0x0);
    if(projQ[iCent]) {
      projQ[iCent]->SetLineColor(iCent<4 ? iCent+1 : iCent+2);
      Double_t maximum = projQ[iCent]->GetMaximum();
      if(maximum>1.0)
        projQ[iCent]->Scale(1.0/maximum);
      projQ[iCent]->SetStats(kFALSE);
      projQ[iCent]->Draw("same");
      leg2->AddEntry(projQ[iCent], (iCent==0 ? "all centralities" : Form("centrality class %d", iCent)), "l");
      isGood = kTRUE;
    }
  }
  if(isGood) leg2->Draw();
  isGoodForSaving = isGoodForSaving || isGood;
  return isGoodForSaving;
}


//_________________________________________________________________
Bool_t AliTRDcheckESD::PlotTrackingSummary(Int_t centralityClass, Double_t* trendValues) {

  Bool_t isGoodForSaving=kFALSE;
  
  TLatex *lat=new TLatex();
  lat->SetTextSize(0.07);
  lat->SetTextColor(2);
  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001);
  gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
  gPad->Divide(3,3,0.,0.);
  TList* l=gPad->GetListOfPrimitives();
  // eta-phi distr. for positive TPC tracks
  TVirtualPad* pad = ((TVirtualPad*)l->At(0)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  if(gROOT->FindObject("rangeEtaPhi")) delete gROOT->FindObject("rangeEtaPhi");
  TH2F* rangeEtaPhi = new TH2F("rangeEtaPhi", "", 10, -0.99, +0.99, 10, -3.4, +3.4);
  SetStyle(rangeEtaPhi->GetXaxis(), "#eta", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeEtaPhi->GetYaxis(), "detector #varphi", 0.07, 0.8, kTRUE, 0.05);
  rangeEtaPhi->SetStats(kFALSE);
  rangeEtaPhi->Draw();
  lat->DrawLatex(-0.9, 3.6, "TPC positive ref. tracks");
  
  TH3F* h3F = dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksPos+centralityClass));
  TH2F* h2FtpcP = 0x0;
  Float_t nada=0.0;
  if(h3F && h3F->GetEntries()>10) {
    h2FtpcP = (TH2F*)Proj3D((TH3F*)h3F, 0x0, 1, h3F->GetZaxis()->GetNbins(), nada)->Clone();
    h2FtpcP->SetStats(kFALSE);
    h2FtpcP->Draw("samecolz"); isGoodForSaving = kTRUE;
    isGoodForSaving = kTRUE;
  }
  //-----------------
  // eta-phi distr. for negative TPC tracks
  pad = ((TVirtualPad*)l->At(1)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  lat->DrawLatex(-0.9, 3.6, "TPC negative ref. tracks");
  
  h3F = dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksNeg+centralityClass));
  TH2F* h2FtpcN = 0x0;
  if(h3F && h3F->GetEntries()>10) {
    h2FtpcN = (TH2F*)Proj3D((TH3F*)h3F, 0x0, 1, h3F->GetZaxis()->GetNbins(), nada)->Clone();
    h2FtpcN->SetStats(kFALSE);
    h2FtpcN->Draw("samecolz"); 
    isGoodForSaving = kTRUE;
  }
  //----------------------------------------------
  // eta-phi distr. for positive TRD tracks
  pad = ((TVirtualPad*)l->At(3)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  lat->DrawLatex(-0.9, 3.6, "TRD positive ref. tracks");
  
  h3F = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksPos+centralityClass));
  TH2F* h2FtrdP = 0x0;
  if(h3F && h3F->GetEntries()>10) {
    h2FtrdP = (TH2F*)Proj3D((TH3F*)h3F, 0x0, 1, h3F->GetZaxis()->GetNbins(), nada)->Clone();
    h2FtrdP->SetStats(kFALSE);
    h2FtrdP->SetMaximum((h2FtpcP ? h2FtpcP->GetMaximum() : h2FtrdP->GetMaximum()));
    h2FtrdP->Draw("samecolz"); 
    isGoodForSaving=kTRUE;
  }
  //--------------------------------------------
  // eta-phi distr. for negative TRD tracks
  pad = ((TVirtualPad*)l->At(4)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  lat->DrawLatex(-0.9, 3.6, "TRD negative ref. tracks");
  
  h3F = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksNeg+centralityClass));
  TH2F* h2FtrdN = 0x0;
  if(h3F && h3F->GetEntries()>10) {
    h2FtrdN = (TH2F*)Proj3D((TH3F*)h3F, 0x0, 1, h3F->GetZaxis()->GetNbins(), nada)->Clone();
    h2FtrdN->SetStats(kFALSE);
    h2FtrdN->SetMaximum(h2FtpcN ? h2FtpcN->GetMaximum() : h2FtrdN->GetMaximum());
    h2FtrdN->Draw("samecolz"); 
    isGoodForSaving=kTRUE;
  }
  //----------------------------------------------
  // eta-phi efficiency for positive TRD tracks
  pad = ((TVirtualPad*)l->At(6)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  lat->DrawLatex(-0.9, 3.6, "Efficiency positive tracks");
  
  TH2F* h2Feff = (h2FtrdP ? (TH2F*)h2FtrdP->Clone("h2FeffPos") : 0x0);
  if(h2Feff) {
    h2Feff->Reset();
    h2Feff->SetStats(kFALSE);
    h2Feff->Divide(h2FtrdP, h2FtpcP);
    h2Feff->SetMaximum(1.0);
    if(h2Feff->GetEntries()>1) { h2Feff->Draw("samecolz"); isGoodForSaving=kTRUE; }
  }
  //-------------------------------------------------
  // eta-phi efficiency for negative TRD tracks
  pad = ((TVirtualPad*)l->At(7)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  lat->DrawLatex(-0.9, 3.6, "Efficiency negative tracks");
  
  h2Feff = (h2FtrdN ? (TH2F*)h2FtrdN->Clone("h2FeffNeg") : 0x0);
  if(h2Feff) {
    h2Feff->Reset();
    h2Feff->SetStats(kFALSE);
    h2Feff->Divide(h2FtrdN, h2FtpcN);
    h2Feff->SetMaximum(1.0);
    if(h2Feff->GetEntries()>0.1) { h2Feff->Draw("samecolz"); isGoodForSaving=kTRUE; }
  }
  //-----------------------------------------------------
  // <ntracklets> vs (phi,eta)
  pad = ((TVirtualPad*)l->At(2)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  lat->DrawLatex(-0.9, 3.6, "TRD <N_{tracklets}>");
  
  TProfile2D* hProf2D;
  if((hProf2D = dynamic_cast<TProfile2D*>(fHistos->At(kTRDEtaPhiAvNtrkl+centralityClass)))) {
    if(hProf2D->GetEntries()>10) {
      hProf2D->SetStats(kFALSE);
      hProf2D->SetMinimum(0.);
      hProf2D->SetMaximum(6.);
      if(hProf2D->GetEntries()>1) { hProf2D->Draw("samecolz"); isGoodForSaving = kTRUE; }
    }
  }
  //---------------------------------------------------------
  // TPC-TRD matching efficiency vs pt
  pad = ((TVirtualPad*)l->At(5)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  TH1F* hFeffP = EfficiencyTRD(dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksPos+centralityClass)),
              dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksPos+centralityClass)), kTRUE);
  TH1F* hFeffN = EfficiencyTRD(dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksNeg+centralityClass)),
              dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksNeg+centralityClass)), kTRUE);
  TH1F* hFeffP4 = EfficiencyTRD(dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksPos+centralityClass)),
        dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksPos4+centralityClass)), kTRUE);
  TH1F* hFeffN4 = EfficiencyTRD(dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksNeg+centralityClass)),
        dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksNeg4+centralityClass)), kTRUE);
  TH1F* hFeffP5 = EfficiencyTRD(dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksPos+centralityClass)),
        dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksPos5+centralityClass)), kTRUE);
  TH1F* hFeffN5 = EfficiencyTRD(dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksNeg+centralityClass)),
        dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksNeg5+centralityClass)), kTRUE);
  TH1F* hFeffP6 = EfficiencyTRD(dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksPos+centralityClass)),
        dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksPos6+centralityClass)), kTRUE);
  TH1F* hFeffN6 = EfficiencyTRD(dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksNeg+centralityClass)),
        dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksNeg6+centralityClass)), kTRUE);
  
  TF1* funcConst = new TF1("funcConst", "[0]", 1.0, 3.0);
  
  if(gROOT->FindObject("rangeEffPt2")) delete gROOT->FindObject("rangeEffPt2");
  TH2F* rangeEffPt2=new TH2F("rangeEffPt2", "",10,0.,10.,10,0.,1.4);
  rangeEffPt2->SetStats(kFALSE);
  SetStyle(rangeEffPt2->GetXaxis(), "p_{T} [GeV/c]", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeEffPt2->GetYaxis(), "efficiency", 0.07, 0.8, kTRUE, 0.05);
  rangeEffPt2->Draw();
  lat->DrawLatex(0.5, 1.42, "TRD-TPC matching efficiency");
  //++++++++++++++++++
  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.DrawLine(rangeEffPt2->GetXaxis()->GetXmin(), 0.7, rangeEffPt2->GetXaxis()->GetXmax(), 0.7);
  line.DrawLine(rangeEffPt2->GetXaxis()->GetXmin(), 0.9, rangeEffPt2->GetXaxis()->GetXmax(), 0.9);
  line.SetLineStyle(1);
  line.SetLineWidth(1);
  line.DrawLine(rangeEffPt2->GetXaxis()->GetXmin(), 1.0, rangeEffPt2->GetXaxis()->GetXmax(), 1.0);
  TLegend* leg=new TLegend(0.2, 0.7, 0.6, 0.89);
  leg->SetNColumns(2);
  leg->SetFillColor(0);
  if(hFeffP){
    hFeffP->SetMarkerStyle(24);
    hFeffP->SetMarkerColor(2);
    hFeffP->SetLineColor(2);
    if(trendValues && hFeffP->GetEntries()>1) {
      hFeffP->Fit(funcConst, "QME0", "goff", 1.0, 3.0);
      trendValues[0] = funcConst->GetParameter(0);
      trendValues[1] = funcConst->GetParError(0);
    }
    if(hFeffP->Integral()>0.001) {
      hFeffP->Draw("same"); 
      leg->AddEntry(hFeffP, "positives (#geq 1 tracklet)", "p");
    }
  }
  if(hFeffN){
    hFeffN->SetMarkerStyle(24);
    hFeffN->SetMarkerColor(4);
    hFeffN->SetLineColor(4);
    if(trendValues && hFeffN->GetEntries()>1) {
      hFeffN->Fit(funcConst, "QME0", "goff", 1.0, 3.0);
      trendValues[2] = funcConst->GetParameter(0);
      trendValues[3] = funcConst->GetParError(0);
    }
    if(hFeffN->Integral()>0.001) {
      hFeffN->Draw("same"); 
      leg->AddEntry(hFeffN, "negatives (#geq 1 tracklet)", "p");
    }
  }
  if(hFeffP4){
    hFeffP4->SetMarkerStyle(25);
    hFeffP4->SetMarkerColor(2);
    hFeffP4->SetLineColor(2);
    if(hFeffP4->Integral()>0.001) {
      hFeffP4->Draw("same"); 
      leg->AddEntry(hFeffP4, "positives (4 tracklets)", "p");
    }
  }
  if(hFeffN4){
    hFeffN4->SetMarkerStyle(25);
    hFeffN4->SetMarkerColor(4);
    hFeffN4->SetLineColor(4);
    if(hFeffN4->Integral()>0.001) {
      hFeffN4->Draw("same"); 
      leg->AddEntry(hFeffN4, "negatives (4 tracklets)", "p");
    }
  }
  if(hFeffP5){
    hFeffP5->SetMarkerStyle(26);
    hFeffP5->SetMarkerColor(2);
    hFeffP5->SetLineColor(2);
    if(hFeffP5->Integral()>0.001) {
      hFeffP5->Draw("same"); 
      leg->AddEntry(hFeffP5, "positives (5 tracklets)", "p");
    }
  }
  if(hFeffN5){
    hFeffN5->SetMarkerStyle(26);
    hFeffN5->SetMarkerColor(4);
    hFeffN5->SetLineColor(4);
    if(hFeffN5->Integral()>0.001) {
      hFeffN5->Draw("same"); 
      leg->AddEntry(hFeffN5, "negatives (5 tracklets)", "p");
    }
  }
  if(hFeffP6){
    hFeffP6->SetMarkerStyle(27);
    hFeffP6->SetMarkerColor(2);
    hFeffP6->SetLineColor(2);
    if(hFeffP6->Integral()>0.001) {
      hFeffP6->Draw("same"); 
      leg->AddEntry(hFeffP6, "positives (6 tracklets)", "p");
    }
  }
  if(hFeffN6){
    hFeffN6->SetMarkerStyle(27);
    hFeffN6->SetMarkerColor(4);
    hFeffN6->SetLineColor(4);
    if(hFeffN6->Integral()>0.001) {
      hFeffN6->Draw("same"); 
      leg->AddEntry(hFeffN6, "negatives (6 tracklets)", "p");
    }
  }
  leg->Draw();
  
  //--------------------------------------------------------------
  // Nclusters per TRD track
  pad = ((TVirtualPad*)l->At(8)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.12);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  pad->SetLogz();
  if(gROOT->FindObject("rangeNclsP")) delete gROOT->FindObject("rangeNclsP");
  TH2F* rangeNclsP = new TH2F("rangeNclsP", "", 10, 0.0, 11.99, 10, 0.0, 199.0);
  SetStyle(rangeNclsP->GetXaxis(), "p [GeV/c]", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeNclsP->GetYaxis(), "#clusters", 0.07, 0.8, kTRUE, 0.05);
  rangeNclsP->SetStats(kFALSE);
  rangeNclsP->Draw();
  lat->DrawLatex(1.0, 205., "TRD Clusters / track");
  
  TH2F* hNclsVsP=0x0;
  if((hNclsVsP = dynamic_cast<TH2F*>(fHistos->At(kNClsTrackTRD+centralityClass)))) {
    hNclsVsP->SetStats(kFALSE);
    if(hNclsVsP->GetEntries()>10) {
      hNclsVsP->Draw("samecolz"); isGoodForSaving=kTRUE;
      if(trendValues) {
        TProfile* h2FProf = hNclsVsP->ProfileX("nclsVsPprof");
        h2FProf->Fit(funcConst, "QME0", "goff", 1.0, 3.0);
        trendValues[4] = funcConst->GetParameter(0);
        trendValues[5] = funcConst->GetParError(0);
      }
    }
  }
  
  delete funcConst;
  return isGoodForSaving;
}


//_________________________________________________________________
Bool_t AliTRDcheckESD::PlotPidSummary(Int_t centralityClass, Double_t* trendValues) {

  Bool_t isGoodForSaving=kFALSE;
  
  TLatex *lat=new TLatex();
  lat->SetTextSize(0.07);
  lat->SetTextColor(2);
  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001);
  gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
  gPad->Divide(3,3,0.,0.);
  TList* l=gPad->GetListOfPrimitives();
  // eta-phi distr. for <Qtot> in layer 0
  TVirtualPad* pad;
  TProfile2D* hProf2D;
  if(gROOT->FindObject("rangeEtaPhi2")) delete gROOT->FindObject("rangeEtaPhi2");
  TH2F* rangeEtaPhi = new TH2F("rangeEtaPhi2", "", 10, -0.99, +0.99, 10, -3.4, +3.4);
  SetStyle(rangeEtaPhi->GetXaxis(), "#eta", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeEtaPhi->GetYaxis(), "detector #varphi", 0.07, 0.8, kTRUE, 0.05);
  rangeEtaPhi->SetStats(kFALSE);
  
  for(Int_t iLayer=0; iLayer<6; ++iLayer) {
    pad = ((TVirtualPad*)l->At((iLayer<3 ? iLayer*3 : (iLayer-3)*3+1))); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    rangeEtaPhi->Draw();
    lat->DrawLatex(-0.9, 3.6, Form("TRD <Q_{tot}> Layer %d", iLayer));
    
    if(!(hProf2D = dynamic_cast<TProfile2D*>(fHistos->At(kTRDEtaPhiAvQtot+6*centralityClass+iLayer)))) continue;
    if(hProf2D && hProf2D->GetEntries()>10) {
      hProf2D->SetStats(kFALSE);
      hProf2D->SetMinimum(0.);
      hProf2D->SetMaximum(4.);
      if(hProf2D->GetEntries()>10) { hProf2D->Draw("samecolz"); isGoodForSaving=kTRUE; }
    }
  }
    
  // PH versus slice number
  pad = ((TVirtualPad*)l->At(2)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  if(gROOT->FindObject("rangePHslice")) delete gROOT->FindObject("rangePHslice");
  TH2F* rangePHslice = new TH2F("rangePHslice", "", 10, -0.5, 7.5, 10, 0.0, 2000.0);
  SetStyle(rangePHslice->GetXaxis(), "slice", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangePHslice->GetYaxis(), "PH", 0.07, 0.8, kTRUE, 0.05);
  rangePHslice->SetStats(kFALSE);
  rangePHslice->Draw();
  
  TF1* funcPol1 = new TF1("funcPol1", "[0]+[1]*x", 2.9, 6.4);
  
  TH2F* h2F;
  TH1D* hF;
  if((h2F = dynamic_cast<TH2F*>(fHistos->At(kPHSlice+centralityClass)))) {
    if(h2F && h2F->GetEntries()>10) {
      hF = Proj2D(h2F);
      h2F->SetStats(kFALSE);
      h2F->Draw("samecolz");
      isGoodForSaving=kTRUE;
      if(trendValues) {
        hF->Fit(funcPol1, "QME0", "goff", 2.9, 6.4);
        trendValues[6] = funcPol1->GetParameter(0);
        trendValues[7] = funcPol1->GetParError(0);
        trendValues[8] = funcPol1->GetParameter(1);
        trendValues[9] = funcPol1->GetParError(1);
      }
      hF->SetLineWidth(2);
      hF->Draw("same");
    }
  }
  delete funcPol1;

  // Qtot vs P
  pad = ((TVirtualPad*)l->At(5)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  pad->SetLogz();
  if(gROOT->FindObject("rangeQtotP")) delete gROOT->FindObject("rangeQtotP");
  TH2F* rangeQtotP = new TH2F("rangeQtotP", "", 10, 0.0, 11.99, 10, 0.0, 11.99);
  SetStyle(rangeQtotP->GetXaxis(), "P [GeV/c]", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeQtotP->GetYaxis(), "Q_{tot}", 0.07, 0.8, kTRUE, 0.05);
  rangeQtotP->SetStats(kFALSE);
  rangeQtotP->Draw();
  
  if((h2F = dynamic_cast<TH2F*>(fHistos->At(kQtotP+centralityClass)))) {
    if(h2F && h2F->GetEntries()>10) {
      h2F->SetStats(kFALSE);
      h2F->Draw("samecolz");
      isGoodForSaving=kTRUE;
      hF = Proj2D(h2F);
      hF->SetLineWidth(2);
      hF->Draw("same");
      if(trendValues) {
        trendValues[10] = hF->GetBinContent(hF->FindBin(1.0));
        trendValues[11] = hF->GetBinError(hF->FindBin(1.0));
      }
    }
  }
  return isGoodForSaving;
}

//__________________________________________________________________________________________________
void AliTRDcheckESD::DrawTRDGrid() {
  //
  //   Draw a grid of lines showing the TRD supermodule and stack structure in (eta,phi) coordinates.
  //   The canvas on which to draw must already exist.
  //
  TLine line;
  line.SetLineColor(2);
  line.SetLineWidth(1);
  line.SetLineStyle(2);
  for(Int_t i=0; i<=9; ++i) {
    line.DrawLine(-1.0, 2.0*TMath::Pi()/18.0*i, +1.0, 2.0*TMath::Pi()/18.0*i);
    line.DrawLine(-1.0, -2.0*TMath::Pi()/18.0*i, +1.0, -2.0*TMath::Pi()/18.0*i);
  }
  line.DrawLine(-0.85, -3.2, -0.85, +3.2);
  line.DrawLine(-0.54, -3.2, -0.54, +3.2);
  line.DrawLine(-0.16, -3.2, -0.16, +3.2);
  line.DrawLine(+0.16, -3.2, +0.16, +3.2);
  line.DrawLine(+0.54, -3.2, +0.54, +3.2);
  line.DrawLine(+0.85, -3.2, +0.85, +3.2);
}

//_________________________________________________________________
void AliTRDcheckESD::SetStyle(TH1* hist, 
			      Int_t lineStyle, Int_t lineColor, Int_t lineWidth, 
			      Int_t markerStyle, Int_t markerColor, Int_t markerSize) {
  //
  // Set style settings for histograms
  //
  hist->SetLineStyle(lineStyle);
  hist->SetLineColor(lineColor);
  hist->SetLineWidth(lineWidth);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerColor(markerColor);
  hist->SetMarkerSize(markerSize);
}

//____________________________________________________________________
void AliTRDcheckESD::SetStyle(TAxis* axis, const Char_t* title, Float_t titleSize, Float_t titleOffset, Bool_t centerTitle, 
                              Float_t labelSize) {
  //
  // Set style settings for axes
  //
  axis->SetTitle(title);
  axis->SetTitleSize(titleSize);
  axis->SetTitleOffset(titleOffset); 
  axis->CenterTitle(centerTitle);
  axis->SetLabelSize(labelSize);
}

//____________________________________________________________________
void AliTRDcheckESD::FindIsolatedBCs(TH1D* bcHist, Bool_t isIsolated[3500]) {
  //
  // Find the isolated bunch crossings
  //
  Int_t isolationSize = 10;      // number of free bunches in both directions
  for(Int_t bcBin=1; bcBin<=bcHist->GetXaxis()->GetNbins(); ++bcBin) {
    Int_t bc = TMath::Nint(bcHist->GetBinCenter(bcBin));
    if(bc<-0.001 || bc>3499.01) {
      isIsolated[bc] = kFALSE;
      continue;
    }
    Double_t entries = bcHist->GetBinContent(bcBin);
    if(entries<0.001) {
      isIsolated[bc] = kFALSE;
      continue;     // no entries
    }
    
    // check isolation
    isIsolated[bc] = kTRUE;
    for(Int_t ibc = TMath::Max(1,bcBin-isolationSize); ibc <= TMath::Min(3499, bcBin+isolationSize); ++ibc) {
      if(ibc==bcBin) continue;
      if(bcHist->GetBinContent(ibc)>0.01) {
        isIsolated[bc] = kFALSE;
        break;
      }
    }
  }   // end loop over BC bins
  
  cout << "Isolated bunches: " << endl;
  for(Int_t ibc=0; ibc<3500; ++ibc) 
    if(isIsolated[ibc]) cout << "BC #" << ibc << endl; 
}
