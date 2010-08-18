#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>
#include <string>
#include <TCanvas.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "AliFMDAnalysisTaskBFCorrelation.h"
#include "AliAnalysisManager.h"
#include "AliESDFMD.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliESDVertex.h"
#include "TMath.h"
#include "AliFMDAnaParameters.h"
//#include "AliFMDGeometry.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliHeader.h"
//#include "TDatabasePDG.h"
//#include "TParticlePDG.h"
#include "AliESDInputHandler.h"

ClassImp(AliFMDAnalysisTaskBFCorrelation)

AliFMDAnalysisTaskBFCorrelation::AliFMDAnalysisTaskBFCorrelation()
: fDebug(0),
  fOutputList(0),
  fInputList(0),
  fInternalList(0),
  fVertexString(0x0),
  fStandalone(kTRUE),
  fEvent(0),
  fnBinsX(200),
  fXmin(-6),
  fXmax(6),
  fnBinsY(20),
  fYmin(0),
  fYmax(2 * TMath::Pi()),
  c(0),
  debug0(0),
  debug1(0)
{
  // Default constructor
  DefineInput (0, TList::Class());
  DefineOutput(0, TList::Class());
}

//_____________________________________________________________________
AliFMDAnalysisTaskBFCorrelation::AliFMDAnalysisTaskBFCorrelation(const char* name, Bool_t SE):
  AliAnalysisTask(name,name),
  fDebug(0),
  fOutputList(0),
  fInputList(0),
  fInternalList(0),
  fVertexString(0x0),
  fStandalone(kTRUE),
  fEvent(0),
  fnBinsX(200),
  fXmin(-6),
  fXmax(6),
  fnBinsY(20),
  fYmin(0),
  fYmax(2 * TMath::Pi()),
  c(0),
  debug0(0),
  debug1(0)
{
  fStandalone = SE;
  if(fStandalone) {
    DefineInput (0, TList::Class());
    DefineInput(1, TObjString::Class());
    DefineOutput(0, TList::Class());
  }
}

//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::CreateOutputObjects()
{
  // Setup the list for storing results, if it does not exist

  if(!fOutputList) {
    fOutputList = new TList();
    fOutputList->SetName("BackgroundCorrected");
  }

  // Setup the list for temporary storage during the analysis
  
  if(!fInternalList) {
    fInternalList = new TList();
    fInternalList->SetName("InternalBFList");
  }
  
  // Set up histograms for analysis and storage. 4 different binnings
  
  for (Int_t i = 1; i <= 4; i++) { 
    
    // Temporary histograms for storing hist per eta summed over phi
    
    TH1D *hESDMult = new TH1D(Form("hESDMult_binning%d", i), 
			      Form("Multiplicity per event vs eta ESD (%d bins)", 30*i), 
			      18*i, -9, 9);
    TH1D *hMCMult  = new TH1D(Form("hMCMult_binning%d", i), 
			      Form("Multiplicity per event vs eta MC-truth (%d bins)", 30*i), 
			      18*i, -9, 9);
    fInternalList->Add(hESDMult);
    fInternalList->Add(hMCMult);

    
    // Histograms for storing the 5 values to calculate b (ESD & MC)
    
    TH1D *pnFnB = new TH1D(Form("pnFnB_binning%d", i),
			   Form("Correlation ESD (%d bins)", i*9), 
			   9*i, 0, 9);
    TH1D *pnF2  = new TH1D(Form("pnF2_binning%d", i),
			   Form("Fwd^2 ESD (%d bins)", i*9),
			   9*i, 0, 9);
    TH1D *pnB2  = new TH1D(Form("pnB2_binning%d", i), 
			   Form("Bwd^2 ESD (%d bins)", i*9), 
			   9*i, 0, 9);
    TH1D *pnF   = new TH1D(Form("pnF_binning%d", i), 
			   Form("Fwd ESD (%d bins)", i*9), 
			   9*i, 0, 9);
    TH1D *pnB   = new TH1D(Form("pnB_binning%d", i),
			   Form("Bwd ESD (%d bins)", i*9),
			   9*i, 0, 9);    
    
    TH1D *pnFnB_MC = new TH1D(Form("pnFnB_MC_binning%d", i),
			      Form("Correlation MC (%d bins)", i*9), 
			      9*i, 0, 9);
    TH1D *pnF2_MC  = new TH1D(Form("pnF2_MC_binning%d", i),
			      Form("Fwd^2 MC (%d bins)", i*9), 
			      9*i, 0, 9);
    TH1D *pnB2_MC  = new TH1D(Form("pnB2_MC_binning%d", i),
			      Form("Bwd^2 MC (%d bins)", i*9),
			      9*i, 0, 9);
    TH1D *pnF_MC   = new TH1D(Form("pnF_MC_binning%d", i),
			      Form("Fwd MC (%d bins)", i*9), 
			      9*i, 0, 9);
    TH1D *pnB_MC   = new TH1D(Form("pnB_MC_binning%d", i),
			      Form("Backward MC (%d bins)" ,i*9), 
			      9*i, 0, 9);
    fOutputList->Add(pnFnB);
    fOutputList->Add(pnF2);
    fOutputList->Add(pnB2);
    fOutputList->Add(pnF);
    fOutputList->Add(pnB);
    
    fOutputList->Add(pnFnB_MC);
    fOutputList->Add(pnF2_MC);
    fOutputList->Add(pnB2_MC);
    fOutputList->Add(pnF_MC);
    fOutputList->Add(pnB_MC);
    
    // Debugging histograms : "dNdEta"
    
    //    TH1D *hESDMultvsEta = new TH1D(Form("hESDMultvsEta_binning%d" ,i), 
    //			   Form("Mult vs Eta (%d bins, ESD)", 18*i),
    //			   18*i, -9, 9);
    TH1D *hESDMultvsEta = new TH1D(Form("hESDMultvsEta_binning%d" ,i), 
				   Form("Mult vs Eta (%d bins, ESD)", i*18),
				   200, -6, 6);
    TH1D *hMCMultvsEta = new TH1D(Form("hMCMultvsEta_binning%d", i),
				  Form("Mult vs Eta (%d bins, MC)", 18*i), 
				  200, -6, 6);
    //TH1D *hMCMultvsEta = new TH1D(Form("hMCMultvsEta_binning%d", i),
    //			  Form("Mult vs Eta (%d bins, MC)", 18*i), 
    //			  18*i, -9, 9);
    fOutputList->Add(hESDMultvsEta);
    fOutputList->Add(hMCMultvsEta);

    TH1D *hESDMultvsEtaC = new TH1D(Form("hESDMultvsEtaC_binning%d" ,i), 
				   Form("Mult vs Eta C (%d bins, ESD)", 18*i),
				   18*i, -9, 9);
    TH1D *hMCMultvsEtaC  = new TH1D(Form("hMCMultvsEtaC_binning%d", i),
				  Form("Mult vs Eta C (%d bins, MC)", 18*i), 
				  18*i, -9, 9);
    fOutputList->Add(hESDMultvsEtaC);
    fOutputList->Add(hMCMultvsEtaC);

    // Debugging histograms : Different distributions

    TH2D *hESDBwdDist = new TH2D(Form("hESDBwdDist_binning%d", i),
				 Form("Distribution of Bwd values (%d bins)", i*9),
				 9*i, 0, 9, 20, 0, 20);
    
    TH2D *hESDBwd2Dist = new TH2D(Form("hESDBwd2Dist_binning%d", i),
				  Form("Distribution of Bwd^2 values (%d bins)", i*9),
				  9*i, 0, 9, 400, 0, 400);
    
    TH2D *hMCBwdDist = new TH2D(Form("hMCBwdDist_binning%d", i),
				Form("Distribution of Bwd values (%d bins)", i*9),
				9*i, 0, 9, 20, 0, 20);
    
    TH2D *hMCBwd2Dist = new TH2D(Form("hMCBwd2Dist_binning%d", i),
				 Form("Distribution of Bwd^2 values (%d bins)", i*9),
				 9*i, 0, 9, 400, 0, 400);
    fOutputList->Add(hESDBwdDist);
    fOutputList->Add(hESDBwd2Dist);
    fOutputList->Add(hMCBwdDist);
    fOutputList->Add(hMCBwd2Dist);
  }
}

//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::ConnectInputData(Option_t */*option*/)
{
  if(fStandalone) {
    fInputList   = (TList*)GetInputData(0);
    fVertexString = (TObjString*)GetInputData(1);
  }
}

//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::Exec(Option_t */*option*/)
{
  fEvent++;
  if (fEvent % 1000 == 0) 
    std::cout << "Event # " << fEvent << std::endl;
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  fVertexString = (TObjString*)fInputList->At(0);
  Int_t vtxbin   = fVertexString->GetString().Atoi();
  if(vtxbin !=4) return;
  SetBounds();
  CountESDHits();
  CalculateValues("ESD");
  
  if(pars->GetProcessPrimary())
    ProcessPrimary();
  
  if(fStandalone) {
    PostData(0, fOutputList); 
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::CountESDHits() {

  TH2F *dNdEtadPhiHist = (TH2F*)fInputList->FindObject("dNdetadphiHistogramSPDTrVtx");

  for (Int_t i = 1; i<=4; i++) {
    
    TH1D *hESDMult = (TH1D*)fInternalList->FindObject(Form("hESDMult_binning%d", i));
    hESDMult->Reset();
  }
  
  Int_t etamax = dNdEtadPhiHist->GetNbinsX();
  Int_t phimax = dNdEtadPhiHist->GetNbinsY();
  
  TH1D *hprod = dNdEtadPhiHist->ProjectionX();
  
  TH1D *hESDMultvsEta = (TH1D*)fOutputList->FindObject(Form("hESDMultvsEta_binning%d" ,4));
  hESDMultvsEta->Add(hprod);
  
  for (Int_t etabin = 1; etabin <= etamax; etabin++) {
    Float_t val = 0;
    Float_t eta = dNdEtadPhiHist->GetXaxis()->GetBinCenter(etabin);
    
    
    for (Int_t i = 1; i <= 4; i++) {
       TH1D *hESDMult = (TH1D*)fInternalList->FindObject(Form("hESDMult_binning%d", i));
       hESDMult->Fill(eta, hprod->GetBinContent(etabin));
    }
  }
  
  //hESDMult->Add(hprod);
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::CalculateValues(TString type) {

  type.ToUpper();

  const char *stype = type.Data(); 

  TString sinput;
  TString soutput;
  TString soutputC;
  TString snFnB;
  TString snF2;
  TString snB2;
  TString snF;
  TString snB;
 
  if (type == "ESD") {
    snFnB  = "pnFnB";
    snF2   = "pnF2";
    snB2   = "pnB2";
    snF    = "pnF";
    snB    = "pnB";
  } else {
    snFnB  = "pnFnB_MC";
    snF2   = "pnF2_MC";
    snB2   = "pnB2_MC";
    snF    = "pnF_MC";
    snB    = "pnB_MC";
  }

  for (Int_t i = 1; i <= 4; i++) {
    
    sinput   = Form("h%sMult_binning%d", stype, i);
    soutput  = Form("h%sMultvsEta_binning%d", stype, i);
    soutputC = Form("h%sMultvsEtaC_binning%d", stype, i);
    
    TH1D *input  = (TH1D*)fInternalList->FindObject(sinput); 
    TH1D *output = (TH1D*)fOutputList->FindObject(soutput);
    TH1D *outputC = (TH1D*)fOutputList->FindObject(soutputC);

    //output->Add(input);
    
    for (Int_t bin = 1; bin <= input->GetNbinsX()/2; bin++) {
      
      Double_t eta    = input->GetBinCenter(bin);
      Double_t bwd    = input->GetBinContent(bin);
      Double_t bwd2   = TMath::Power(bwd, 2);
      Double_t fwd    = input->GetBinContent(input->FindBin(TMath::Abs(eta)));
      Double_t fwd2   = TMath::Power(fwd, 2);
      Double_t fwdbwd = fwd*bwd; 

      snFnB += Form("_binning%d", i);
      snF2  += Form("_binning%d", i);
      snB2  += Form("_binning%d", i);
      snF   += Form("_binning%d", i);
      snB   += Form("_binning%d", i);
      
      TH1D *nFnB  = static_cast<TH1D*>(fOutputList->FindObject(snFnB)); 
      TH1D *nF2   = static_cast<TH1D*>(fOutputList->FindObject(snF2)); 
      TH1D *nB2   = static_cast<TH1D*>(fOutputList->FindObject(snB2)); 
      TH1D *nF    = static_cast<TH1D*>(fOutputList->FindObject(snF)); 
      TH1D *nB    = static_cast<TH1D*>(fOutputList->FindObject(snB)); 
      
      nFnB->Fill(TMath::Abs(eta), fwdbwd);
      nF2->Fill(TMath::Abs(eta),  fwd2);
      nB2->Fill(TMath::Abs(eta),  bwd2);
      nF->Fill(TMath::Abs(eta),   fwd);
      nB->Fill(TMath::Abs(eta),   bwd);
      
      outputC->Fill(TMath::Abs(eta), fwd);
      outputC->Fill(eta, bwd);

      if (type == "ESD") {
	
	snFnB.Remove(5,9);
	snF2.Remove(4,9);
	snB2.Remove(4,9);
	snF.Remove(3,9);
	snB.Remove(3,9);
	TH2D *hESDBwdDist = (TH2D*)fOutputList->FindObject(Form("hESDBwdDist_binning%d", i));
	TH2D *hESDBwd2Dist = (TH2D*)fOutputList->FindObject(Form("hESDBwd2Dist_binning%d", i));
	hESDBwdDist->Fill(TMath::Abs(eta), bwd);
	hESDBwd2Dist->Fill(TMath::Abs(eta), bwd2);

      } else {
	
	snFnB.Remove(8,9);
	snF2.Remove(7,9);
	snB2.Remove(7,9);
	snF.Remove(6,9);
	snB.Remove(6,9);
	TH2D *hMCBwdDist = (TH2D*)fOutputList->FindObject(Form("hMCBwdDist_binning%d", i));
	TH2D *hMCBwd2Dist = (TH2D*)fOutputList->FindObject(Form("hMCBwd2Dist_binning%d", i));
	hMCBwdDist->Fill(TMath::Abs(eta), bwd);
	hMCBwd2Dist->Fill(TMath::Abs(eta), bwd2);
      }
    }
  }
}

//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::SetBounds() {
  
  TH2D *hTemp = (TH2D*)fInputList->FindObject("multTrVtx_FMD1I_vtxbin0");
  
  fnBinsX = hTemp->GetNbinsX();
  fnBinsY = hTemp->GetNbinsY();
  fXmin   = hTemp->GetXaxis()->GetXmin();
  fXmax   = hTemp->GetXaxis()->GetXmax();
  fYmin   = hTemp->GetYaxis()->GetXmin();
  fYmax   = hTemp->GetYaxis()->GetXmax();
}  

//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::Terminate(Option_t */*option*/) {
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::ProcessPrimary() {
  
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  AliMCEvent* mcEvent = eventHandler->MCEvent();
  if(!mcEvent)
    return;
    
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  AliMCParticle* particle = 0;
  AliStack* stack = mcEvent->Stack();
  
  TH1F* hPrimary = (TH1F*)fOutputList->FindObject("hMultvsEta");
  AliHeader* header            = mcEvent->Header();
  AliGenEventHeader* genHeader = header->GenEventHeader();
  
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
  
  if (!pythiaGenHeader) {
    std::cout<<" no pythia header!"<<std::endl;
    return; 
  }

  Int_t pythiaType = pythiaGenHeader->ProcessType();
  
  if(pythiaType==92||pythiaType==93){
    std::cout<<"single diffractive"<<std::endl;
    return;
  }
  if(pythiaType==94){
    std::cout<<"double diffractive"<<std::endl;
    return;
  }


  TArrayF vertex;
  genHeader->PrimaryVertex(vertex);
  if(TMath::Abs(vertex.At(2)) > pars->GetVtxCutZ())
    return;
  //Double_t delta           = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
  //Double_t vertexBinDouble = (vertex.At(2) + pars->GetVtxCutZ()) / delta;
  //Int_t    vertexBin       = (Int_t)vertexBinDouble;
    
  //Bool_t firstTrack = kTRUE;
  
  // we loop over the primaries only unless we need the hits (diagnostics running slowly)
  Int_t nTracks = stack->GetNprimary();
  //  if(pars->GetProcessHits())
  //  nTracks = stack->GetNtrack();
  TH1D *hMCMultvsEta = (TH1D*)fOutputList->FindObject(Form("hMCMultvsEta_binning%d" ,4));  
  
  for (Int_t i = 1; i <= 4; i++) {
    
    TH1D *hMCMult = (TH1D*)fInternalList->FindObject(Form("hMCMult_binning%d", i));
    hMCMult->Reset();

    for(Int_t ii = 0 ;ii<nTracks;ii++) {
      particle = (AliMCParticle*) mcEvent->GetTrack(ii);
      if(!particle)
      continue;
      
      if(stack->IsPhysicalPrimary(ii) && particle->Charge() != 0) {
	if (i == 1) 
	  hPrimary->Fill(particle->Eta());
	if(i==4)
	  hMCMultvsEta->Fill(particle->Eta());
	hMCMult->Fill(particle->Eta());
      }
    }
  }
  CalculateValues("MC");
}
//_____________________________________________________________________
//
//
// EOF
