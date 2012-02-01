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
  fnBinsX(0),
  fXmin(0),
  fXmax(0),
  fnBinsY(0),
  fYmin(0),
  fYmax(0)
 
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
  fnBinsX(0),
  fXmin(0),
  fXmax(0),
  fnBinsY(0),
  fYmin(0),
  fYmax(0)
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
  
  // Nomenclature for naming histograms
  // SE    - Used in histogram names if used for containing Single Event information
  // Rebin - Used in histogram names if they have been rebinned for analysis purposes

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
  
  // Get the bounds for the input histogram.
  // This version is optimized for 200 bins in the eta range (-4, 6)

  AliFMDAnaParameters *pars = AliFMDAnaParameters::Instance();
  TH2F *hDefault = (TH2F*)pars->GetBackgroundCorrection(1, 'I', 1);

  fnBinsX = hDefault->GetNbinsX();
  fnBinsY = hDefault->GetNbinsY();
  fXmin   = hDefault->GetXaxis()->GetXmin();
  fXmax   = hDefault->GetXaxis()->GetXmax();
  fYmin   = hDefault->GetYaxis()->GetXmin();
  fYmax   = hDefault->GetYaxis()->GetXmax();
  
  // Histogram to contain an event of MC data. Same dimension as ESD event

  TH2F *hEtaPhiParticleMap = new TH2F("hEtaPhiParticleMap",
				      "Distributions of MC particles (Eta,Phi)",
				      fnBinsX, fXmin, fXmax,
				      fnBinsY, fYmin, fYmax);
  hEtaPhiParticleMap->Sumw2();

  fInternalList->Add(hEtaPhiParticleMap);

  // Create histograms with same binning as input for Response analysis
  // and control histograms. One temporary histogram is also created to
  // avoid the new and delete call. Must be reset after use !

  TH1D *hSEMultESD      = new TH1D("hSEMultESD", "Multiplicity", fnBinsX, fXmin, fXmax);
  TH1D *hSEMultMC       = new TH1D("hSEMultMC",  "Multiplicity", fnBinsX, fXmin, fXmax);

  TH1D *hTemp           = new TH1D("hTemp", "Temporary histogram", fnBinsX, fXmin, fXmax);

  fInternalList->Add(hSEMultESD);
  fInternalList->Add(hSEMultMC);
  fInternalList->Add(hTemp);

  TH1F *hMultvsEtaESD   = new TH1F("hMultvsEtaESD", "Multiplicity vs Eta (ESD)", fnBinsX, fXmin, fXmax);
  TH1F *hMultvsEtaMC    = new TH1F("hMultvsEtaMC",  "Multiplicity vs Eta (MC)",  fnBinsX, fXmin, fXmax);
  
  TH2F *hResponseMatrix = new TH2F("hResponseMatrix", "Response matrix", 151, -0.5, 150.5, 151, -0.5, 150.5);

  fOutputList->Add(hMultvsEtaESD);
  fOutputList->Add(hMultvsEtaMC);
  fOutputList->Add(hResponseMatrix);

  // Set up histograms for analysis and storage. 3 different binnings

  for (Int_t i = 1; i <= 4; i++) { 
    if (i == 3) continue;
    
    Int_t nBinsX  = fnBinsX/(5*i);

    // Histograms for the event-by-event analysis
    
    TH1D *hSERebinMultESD        = new TH1D(Form("hSERebinMultESD_binning%d", i), 
					    Form("Multiplicity per event vs eta ESD (%d bins)", nBinsX), 
					    nBinsX, fXmin, fXmax);
    
    TH1D *hSERebinMultMC         = new TH1D(Form("hSERebinMultMC_binning%d", i), 
					    Form("Multiplicity per event vs eta MC-truth (%d bins)", nBinsX), 
					    nBinsX, fXmin, fXmax);
    
    TH1D *hSERebinMultMirrorESD  = new TH1D(Form("hSERebinMultMirrorESD_binning%d", i), 
					    Form("Multiplicity per event vs eta ESD Mirrored (%d bins)", nBinsX), 
					    nBinsX, fXmin, fXmax);
    
    TH1D *hSERebinMultMirrorMC   = new TH1D(Form("hSERebinMultMirrorMC_binning%d", i), 
					    Form("Multiplicity per event vs eta MC-truth Mirrored (%d bins)", nBinsX), 
					    nBinsX, fXmin, fXmax);
    
    
    TH1D *hSERebinMultWESD       = new TH1D(Form("hSERebinMultWESD_binning%d", i), 
					    Form("Multiplicity per event vs eta ESD (%d bins)", nBinsX), 
					    nBinsX, fXmin, fXmax);
    
    TH1D *hSERebinMultWMC        = new TH1D(Form("hSERebinMultWMC_binning%d", i), 
					    Form("Multiplicity per event vs eta MC-truth (%d bins)", nBinsX), 
					    nBinsX, fXmin, fXmax);
    
    TH1D *hSERebinMultMirrorWESD = new TH1D(Form("hSERebinMultMirrorWESD_binning%d", i), 
					    Form("Multiplicity per event vs eta ESD Mirrored (%d bins)", nBinsX), 
					    nBinsX, fXmin, fXmax);
    
    TH1D *hSERebinMultMirrorWMC  = new TH1D(Form("hSERebinMultMirrorWMC_binning%d", i), 
					    Form("Multiplicity per event vs eta MC-truth Mirrored (%d bins)", nBinsX), 
					    nBinsX, fXmin, fXmax);
    
    fInternalList->Add(hSERebinMultESD);
    fInternalList->Add(hSERebinMultMC);
    fInternalList->Add(hSERebinMultMirrorESD);
    fInternalList->Add(hSERebinMultMirrorMC);

    fInternalList->Add(hSERebinMultWESD);
    fInternalList->Add(hSERebinMultWMC);
    fInternalList->Add(hSERebinMultMirrorWESD);
    fInternalList->Add(hSERebinMultMirrorWMC);
    
    // Histograms for storing the acummulated parameters.

    TH1F *hRebinnESD        = new TH1F(Form("hRebinnESD_binning%d", i),
				  Form("Counts ESD (%d bins)", nBinsX),
				  nBinsX, fXmin, fXmax);
    
    TH1F *hRebinnMirrorESD  = new TH1F(Form("hRebinnMirrorESD_binning%d", i),
				  Form("Counts ESD Mirrored (%d bins)", nBinsX),
				  nBinsX, fXmin, fXmax);
    
    TH1F *hRebinn2ESD       = new TH1F(Form("hRebinn2ESD_binning%d", i),
				  Form("Counts^2 ESD (%d bins)", nBinsX),
				  nBinsX, fXmin, fXmax);
    
    TH1F *hRebinn2MirrorESD = new TH1F(Form("hRebinn2MirrorESD_binning%d", i),
				  Form("Counts^2 ESD Mirrored (%d bins)", nBinsX),
				  nBinsX, fXmin, fXmax);
    
    TH1F *hRebinnfbESD      = new TH1F(Form("hRebinnfbESD_binning%d", i),
				  Form("Fwd*bwd ESD (%d bins)", nBinsX),
				  nBinsX, fXmin, fXmax);
    

    TH1F *hRebinnMC         = new TH1F(Form("hRebinnMC_binning%d", i),
				  Form("Counts MC (%d bins)", nBinsX),
				  nBinsX, fXmin, fXmax);
    
    TH1F *hRebinnMirrorMC   = new TH1F(Form("hRebinnMirrorMC_binning%d", i),
				  Form("Counts MC Mirrored (%d bins)", nBinsX),
				  nBinsX, fXmin, fXmax);
    
    TH1F *hRebinn2MC        = new TH1F(Form("hRebinn2MC_binning%d", i),
				  Form("Counts^2 MC (%d bins)", nBinsX),
				  nBinsX, fXmin, fXmax);
    
    TH1F *hRebinn2MirrorMC  = new TH1F(Form("hRebinn2MirrorMC_binning%d", i),
				  Form("Counts^2 MC Mirrored (%d bins)", nBinsX),
				  nBinsX, fXmin, fXmax);
    
    TH1F *hRebinnfbMC       = new TH1F(Form("hRebinnfbMC_binning%d", i),
				  Form("Fwd*bwd MC (%d bins)", nBinsX),
				  nBinsX, fXmin, fXmax);
    
    fOutputList->Add(hRebinnESD);
    fOutputList->Add(hRebinnMirrorESD);
    fOutputList->Add(hRebinn2ESD);
    fOutputList->Add(hRebinn2MirrorESD);
    fOutputList->Add(hRebinnfbESD);
    
    fOutputList->Add(hRebinnMC);
    fOutputList->Add(hRebinnMirrorMC);
    fOutputList->Add(hRebinn2MC);
    fOutputList->Add(hRebinn2MirrorMC);
    fOutputList->Add(hRebinnfbMC);

    // Histograms for storing the weights for the acummulated parameters.

    TH1F *hRebinnWESD        = new TH1F(Form("hRebinnWESD_binning%d", i),
				   Form("Counts ESD Weights (%d bins)", nBinsX),
				   nBinsX, fXmin, fXmax);
    
    TH1F *hRebinnMirrorWESD  = new TH1F(Form("hRebinnMirrorWESD_binning%d", i),
				   Form("Counts ESD Weights (%d bins)", nBinsX),
				   nBinsX, fXmin, fXmax);
    
    TH1F *hRebinn2WESD       = new TH1F(Form("hRebinn2WESD_binning%d", i),
				   Form("Counts^2 ESD Weights (%d bins)", nBinsX),
				   nBinsX, fXmin, fXmax);
    
    TH1F *hRebinn2MirrorWESD = new TH1F(Form("hRebinn2MirrorWESD_binning%d", i),
				   Form("Counts^2 ESD Weights (%d bins)", nBinsX),
				   nBinsX, fXmin, fXmax);
    
    TH1F *hRebinnfbWESD      = new TH1F(Form("hRebinnfbWESD_binning%d", i),
				   Form("Fwd*bwd ESD Weights (%d bins)", nBinsX),
				   nBinsX, fXmin, fXmax);
    
    
    TH1F *hRebinnWMC         = new TH1F(Form("hRebinnWMC_binning%d", i),
				   Form("Counts MC weights (%d bins)", nBinsX),
				   nBinsX, fXmin, fXmax);
    
    TH1F *hRebinnMirrorWMC   = new TH1F(Form("hRebinnMirrorWMC_binning%d", i),
				   Form("Counts MC weights (%d bins)", nBinsX),
				   nBinsX, fXmin, fXmax);
    
    TH1F *hRebinn2WMC        = new TH1F(Form("hRebinn2WMC_binning%d", i),
				   Form("Counts^2 MC weights (%d bins)", nBinsX),
				   nBinsX, fXmin, fXmax);

    TH1F *hRebinn2MirrorWMC  = new TH1F(Form("hRebinn2MirrorWMC_binning%d", i),
				   Form("Counts^2 MC weights (%d bins)", nBinsX),
				   nBinsX, fXmin, fXmax);
    
    TH1F *hRebinnfbWMC       = new TH1F(Form("hRebinnfbWMC_binning%d", i),
				   Form("Fwd*bwd MC weights (%d bins)", nBinsX),
				   nBinsX, fXmin, fXmax);

    fOutputList->Add(hRebinnWESD);
    fOutputList->Add(hRebinnMirrorWESD);
    fOutputList->Add(hRebinn2WESD);
    fOutputList->Add(hRebinn2MirrorWESD);
    fOutputList->Add(hRebinnfbWESD);
 
    fOutputList->Add(hRebinnWMC);
    fOutputList->Add(hRebinnMirrorWMC);
    fOutputList->Add(hRebinn2WMC);
    fOutputList->Add(hRebinn2MirrorWMC);
    fOutputList->Add(hRebinnfbWMC);
    
    // Histograms for the final result
    /*
    TH1D *hBFFcor = new TH1D(Form("hBFFcor_binning%d", i),
			     Form("Forward - backward correlations F (%d bins)", nBinsXBF),
			     nBinsXBF, 0, fXmax);
    hBFFcor->GetXaxis()->SetTitle("#eta");
    hBFFcor->GetXaxis()->SetTitle("b");
    TH1D *hBFBcor = new TH1D(Form("hBFBcor_binning%d", i),
			     Form("Forward - backward correlations B (%d bins)", nBinsXBF),
			     nBinsXBF, 0, fXmax);
    hBFBcor->GetXaxis()->SetTitle("#eta");
    hBFBcor->GetXaxis()->SetTitle("b");

    TH1D *hBFFcor_MC = new TH1D(Form("hBFFcor_MC_binning%d", i),
				Form("Forward - backward correlations F (%d bins) MC", nBinsXBF),
				nBinsXBF, 0, fXmax);
    hBFFcor_MC->GetXaxis()->SetTitle("#eta");
    hBFFcor_MC->GetXaxis()->SetTitle("b");
    TH1D *hBFBcor_MC = new TH1D(Form("hBFBcor_MC_binning%d", i),
				Form("Forward - backward correlations B (%d bins) MC", nBinsXBF),
				nBinsXBF, 0, fXmax);
    hBFBcor_MC->GetXaxis()->SetTitle("#eta");
    hBFBcor_MC->GetXaxis()->SetTitle("b");

    fOutputList->Add(hBFFcor);
    fOutputList->Add(hBFBcor);
    fOutputList->Add(hBFFcor_MC);
    fOutputList->Add(hBFBcor_MC);
    */
    
    // Temporary histogram to avoid new-delete

    TH1D* hRebinTemp = new TH1D(Form("hRebinTemp_binning%d", i),
				  Form("Temporary histogram (%d bins)", nBinsX),
				  nBinsX, fXmin, fXmax);

    fInternalList->Add(hRebinTemp);
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
void AliFMDAnalysisTaskBFCorrelation::Exec(Option_t */*option*/) {

  fEvent++;
  //if (fEvent % 1000 == 0) 
  //  std::cout << "Event # " << fEvent << std::endl;
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  fVertexString = (TObjString*)fInputList->At(0);
  
  //Int_t vtxbin = fVertexString->GetString().Atoi();
  //if (vtxbin != 5) return;

  ProjectAndMirror("ESD");
  CalculateValues("ESD");
  
  if(pars->GetProcessPrimary()) {
    
    ProcessPrimary();

    CreateResponseMatrix();
  }
  
  
  if(fStandalone) {
    PostData(0, fOutputList); 
  }
}

//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::ProjectAndMirror(TString sType) {
  
  sType.ToUpper();
  
  if (!sType.Contains("ESD") && !sType.Contains("MC")) {
    std::cout << "Wrong type specification for 'ProjectAndMirror'" << std::endl;
    return;
  }
  
  // Get Single Event histograms for storing hits without rebinning
  
  TH1D *hMult = dynamic_cast<TH1D*>(fInternalList->FindObject(Form("hSEMult%s", sType.Data())));
  if(!hMult) {
    AliWarning("no hist - returning"); 
    return; 
  }
  hMult->Reset();
  
  // Create generic names for retrieving histograms 
  
  TList *list = 0; // List for getting either ESD or MC "hit map"
  
  TString sEtaPhiMap;
  
  if (sType.Contains("ESD")) {
    
    list = fInputList;
    
    sEtaPhiMap = "dNdetadphiHistogramSPDTrVtx";
  }
  
  if (sType.Contains("MC")) {
    
    list = fInternalList;
    
    sEtaPhiMap = "hEtaPhiParticleMap";
  }
  
  TString sMult("hSERebinMult");
  TString sMultMirror("hSERebinMultMirror");
  TString sMultW("hSERebinMultW");
  TString sMultMirrorW("hSERebinMultMirrorW");
  
  sType += "_binning";
  
  sMult        += sType;
  sMultMirror  += sType;
  sMultW       += sType;
  sMultMirrorW += sType;
  
  // Get the 2D histogram containing the particles of the event being analyzed
  
  TH2F *hEtaPhiMap = dynamic_cast<TH2F*>(list->FindObject(sEtaPhiMap));
  
  TH1D *hProjection = hEtaPhiMap->ProjectionX("hTemporary");
  hMult->Add(hProjection);  
  
  // Loop over the 3 binnings
  
  for (Int_t i = 1; i<=4; i++) {
    if (i == 3) continue;
  
    // Project the 2D "hit map" onto the eta-axis and rebin
    
    TH1D *hProjRebin = hEtaPhiMap->ProjectionX("hProjRebin");
    hProjRebin->Rebin(i*5);    
   
    // Retrieve the histograms to store the Singe Event information and reset
    
    TH1D *hSEMult        = static_cast<TH1D*>(fInternalList->FindObject(Form("%s%d", sMult.Data(),        i)));
    hSEMult->Reset();
    
    TH1D *hSEMultMirror  = static_cast<TH1D*>(fInternalList->FindObject(Form("%s%d", sMultMirror.Data(),  i)));
    hSEMultMirror->Reset();
    
    TH1D *hSEMultW       = static_cast<TH1D*>(fInternalList->FindObject(Form("%s%d", sMultW.Data(),       i)));
    hSEMultW->Reset();
    
    TH1D *hSEMultMirrorW = static_cast<TH1D*>(fInternalList->FindObject(Form("%s%d", sMultMirrorW.Data(), i)));
    hSEMultMirrorW->Reset();
    
    // Fill the histograms with the Single Event information
    
    hSEMult->Add(hProjRebin);
    
    for (Int_t bin = 1; bin <= hSEMult->GetNbinsX(); bin++) {
      
      hSEMultMirror->SetBinContent(hSEMultMirror->FindBin(-hProjRebin->GetBinCenter(bin)), hProjRebin->GetBinContent(bin));
      hSEMultMirror->SetBinError(hSEMultMirror->FindBin(-hProjRebin->GetBinCenter(bin)), hProjRebin->GetBinError(bin));
      
      // hMultDist->Fill(bin, hMultTemp->GetBinContent(bin));
    }
    hProjRebin->Delete();
  }
}

//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::CalculateValues(TString sType) {
  
  sType.ToUpper();
  
  if (!sType.Contains("ESD") && !sType.Contains("MC")) {
    std::cout << "Wrong type specification for 'CalculateValues'" << std::endl;
    return;
  }
  
  TString sMult("hSERebinMult");
  TString sMultMirror("hSERebinMultMirror");
  //  TString_t *sMultW("hSEMultW");
  //  TString_t *sMultMirrorW("hSEMultMirrorW");
 
  TString sn("hRebinn");
  TString snMirror("hRebinnMirror");
  TString sn2("hRebinn2");
  TString sn2Mirror("hRebinn2Mirror");
  TString snfb("hRebinnfb");

  sType += "_binning";

  sMult       += sType;
  sMultMirror += sType;

  sn          += sType;
  snMirror    += sType;
  sn2         += sType;
  sn2Mirror   += sType;
  snfb        += sType;
  
  for (Int_t i = 1; i <= 4; i++) {
    if (i == 3) continue;
    
    TH1D *hSEMult          = (TH1D*)fInternalList->FindObject(Form("%s%d", sMult.Data(),       i)); 
    TH1D *hSEMultMirror    = (TH1D*)fInternalList->FindObject(Form("%s%d", sMultMirror.Data(), i));
    /*
    TH1D *hSEMultW         = (TH1D*)fInternalList->FindObject(Form("%s%d", cMultW,       cType, i)); 
    TH1D *hSEMultMirrorW   = (TH1D*)fInternalList->FindObject(Form("%s%d", cMultMirrorW, cType, i));
    */

    TH1F *hn               = (TH1F*)fOutputList->FindObject(Form("%s%d", sn.Data(),        i));
    TH1F *hnMirror         = (TH1F*)fOutputList->FindObject(Form("%s%d", snMirror.Data(),  i));
    TH1F *hn2              = (TH1F*)fOutputList->FindObject(Form("%s%d", sn2.Data(),       i));
    TH1F *hn2Mirror        = (TH1F*)fOutputList->FindObject(Form("%s%d", sn2Mirror.Data(), i));
    TH1F *hnfb             = (TH1F*)fOutputList->FindObject(Form("%s%d", snfb.Data(),      i));
    /*
    TH1F *hnW              = (TH1F*)fOutputList->FindObject(Form("%sW%s_binning%d", cn, cType, i));
    TH1F *hnMirrorW        = (TH1F*)fOutputList->FindObject(Form("%sW%s_binning%d", cnMirror, cType, i));
    TH1F *hn2W             = (TH1F*)fOutputList->FindObject(Form("%sW%s_binning%d", cn2, cType, i));
    TH1F *hn2MirrorW       = (TH1F*)fOutputList->FindObject(Form("%sW%s_binning%d", cn2Mirror, cType, i));
    TH1F *hnfbW            = (TH1F*)fOutputList->FindObject(Form("%sW%s_binning%d", cnfb, cType, i));
    */
    TH1D *hTemp            = (TH1D*)fInternalList->FindObject(Form("hRebinTemp_binning%d", i)); 

    hn->Add(hSEMult);

    hnMirror->Add(hSEMultMirror);

    hTemp->Reset();
    hTemp->Add(hSEMult);
    hTemp->Multiply(hSEMult);
    hn2->Add(hTemp);
        
    hTemp->Reset();
    hTemp->Add(hSEMultMirror);
    hTemp->Multiply(hSEMultMirror);
    hn2Mirror->Add(hTemp);

    hSEMultMirror->Multiply(hSEMult);
    hnfb->Add(hSEMultMirror);
  }
}

//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::MultiplicityVsEta(TString sType) {

  sType.ToUpper();
  
  if (!sType.Contains("ESD") && !sType.Contains("MC")) {
    std::cout << "Wrong type specification for 'MultiplicityVsEta'" << std::endl;
    return;
  }
}
  

//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::CreateResponseMatrix() {

  TH2F *hResponseMatrix = (TH2F*)fOutputList->FindObject("hResponseMatrix");

  TH1D *hSEMultESD = (TH1D*)fInternalList->FindObject("hSEMultESD");
  TH1D *hSEMultMC  = (TH1D*)fInternalList->FindObject("hSEMultMC");
  
  Float_t HitsESD = 0;
  Float_t HitsMC  = 0;
  
  for (Int_t bin = 1; bin<= fnBinsX; bin++) {
    
    if ((hSEMultMC->GetBinLowEdge(bin)  > -3.4) &&
	(hSEMultMC->GetBinLowEdge(bin+1) <  5.)) {

      HitsESD += hSEMultESD->GetBinContent(bin);
      HitsMC  += hSEMultMC->GetBinContent(bin);
    }
  }
  
  hResponseMatrix->Fill(HitsMC, HitsESD);
}

//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::Terminate(Option_t */*option*/) {
  /*
  std::cout << "Terminating !" << std::endl;

  TH1I *hnEvents   = (TH1I*)fOutputList->FindObject("nEvents");
  TH1I *hnMCEvents = (TH1I*)fOutputList->FindObject("nEvents");

  Int_t nEvents   = hnEvents->GetEntries();
  Int_t nMCEvents = hnMCEvents->GetEntries();

  for (Int_t i = 1; i <= 4; i++) {
    if (i == 3) continue;
    
    TH1F *hnFnB = (TH1F*)fOutputList->FindObject(Form("hnFnB_binning%d", i));
    TH1F *hnF2  = (TH1F*)fOutputList->FindObject(Form("hnF2_binning%d", i));
    TH1F *hnB2  = (TH1F*)fOutputList->FindObject(Form("hnB2_binning%d", i));
    TH1F *hnF   = (TH1F*)fOutputList->FindObject(Form("hnF_binning%d", i));
    TH1F *hnB   = (TH1F*)fOutputList->FindObject(Form("hnB_binning%d", i));

    TH1F *hnFnB_MC = (TH1F*)fOutputList->FindObject(Form("hnFnB_MC_binning%d", i));
    TH1F *hnF2_MC  = (TH1F*)fOutputList->FindObject(Form("hnF2_MC_binning%d", i));
    TH1F *hnB2_MC  = (TH1F*)fOutputList->FindObject(Form("hnB2_MC_binning%d", i));
    TH1F *hnF_MC   = (TH1F*)fOutputList->FindObject(Form("hnF_MC_binning%d", i));
    TH1F *hnB_MC   = (TH1F*)fOutputList->FindObject(Form("hnB_MC_binning%d", i));
    
    hnFnB->Scale(1./Float_t(nEvents));
    hnF2->Scale(1./Float_t(nEvents));
    hnB2->Scale(1./Float_t(nEvents));
    hnF->Scale(1./Float_t(nEvents));
    hnB->Scale(1./Float_t(nEvents));

    hnFnB_MC->Scale(1./Float_t(nMCEvents));
    hnF2_MC->Scale(1./Float_t(nMCEvents));
    hnB2_MC->Scale(1./Float_t(nMCEvents));
    hnF_MC->Scale(1./Float_t(nMCEvents));
    hnB_MC->Scale(1./Float_t(nMCEvents));

    for (Int_t bin = 1; bin <= hnFnB->GetNbinsX(); bin++) {
      
      Double_t nFnBav = hnFnB->GetBinContent(bin);
      Double_t nF2av  = hnF2->GetBinContent(bin);
      Double_t nB2av  = hnB2->GetBinContent(bin);
      Double_t nFav   = hnF->GetBinContent(bin);
      Double_t nBav   = hnB->GetBinContent(bin);

      Double_t nFnB_MCav = hnFnB_MC->GetBinContent(bin);
      Double_t nF2_MCav  = hnF2_MC->GetBinContent(bin);
      Double_t nB2_MCav  = hnB2_MC->GetBinContent(bin);
      Double_t nF_MCav   = hnF_MC->GetBinContent(bin);
      Double_t nB_MCav   = hnB_MC->GetBinContent(bin);

      Double_t bF = ((nF2av-nFav*nFav) == 0 ? 0. : (nFnBav-nFav*nBav)/(nF2av-nFav*nFav));
      Double_t bB = ((nF2av-nFav*nFav) == 0 ? 0. : (nFnBav-nFav*nBav)/(nB2av-nBav*nBav));

      Double_t bF_MC = ((nF2_MCav-nF_MCav*nF_MCav) == 0 ? 0. : 
			(nFnB_MCav-nF_MCav*nB_MCav)/(nF2_MCav-nF_MCav*nF_MCav));
      Double_t bB_MC = ((nF2_MCav-nF_MCav*nF_MCav) == 0 ? 0. : 
			(nFnB_MCav-nF_MCav*nB_MCav)/(nB2_MCav-nB_MCav*nB_MCav));
      
      TH1D *hBFFcor = (TH1D*)fOutputList->FindObject(Form("hBFFcor_binning%d", i));
      TH1D *hBFBcor = (TH1D*)fOutputList->FindObject(Form("hBFBcor_binning%d", i));
      TH1D *hBFFcor_MC = (TH1D*)fOutputList->FindObject(Form("hBFFcor_MC_binning%d", i));
      TH1D *hBFBcor_MC = (TH1D*)fOutputList->FindObject(Form("hBFBcor_MC_binning%d", i));

      hBFFcor->SetBinContent(bin, bF);
      hBFBcor->SetBinContent(bin, bB);
      hBFFcor_MC->SetBinContent(bin, bF_MC);
      hBFBcor_MC->SetBinContent(bin, bB_MC);
    }
    }*/
}
 
//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::ProcessPrimary() {
  
  AliMCEventHandler* eventHandler = 
    dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()
				      ->GetMCtruthEventHandler());
  if (!eventHandler) return;

  AliMCEvent* mcEvent = eventHandler->MCEvent();
  if(!mcEvent) return;
    
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  AliMCParticle* particle = 0;
  AliStack* stack = mcEvent->Stack();
  
  //TH1F* hPrimary = (TH1F*)fOutputList->FindObject("hMultvsEta");
  AliHeader* header            = mcEvent->Header();
  AliGenEventHeader* genHeader = header->GenEventHeader();
  /*
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
  */
  TArrayF vertex;
  genHeader->PrimaryVertex(vertex);
  if(TMath::Abs(vertex.At(2)) > pars->GetVtxCutZ())
    return;
  //  Double_t delta           = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
  //  Double_t vertexBinDouble = (vertex.At(2) + pars->GetVtxCutZ()) / delta;
  //  Int_t    vertexBin       = (Int_t)vertexBinDouble;
  
  //if (vertexBin != 5) return;
  
  //Bool_t firstTrack = kTRUE;
  
  // we loop over the primaries only unless we need the hits (diagnostics running slowly)
  Int_t nTracks = stack->GetNprimary();
  //  if(pars->GetProcessHits())
  //  nTracks = stack->GetNtrack();
  
  TH2F *hEtaPhiParticleMap = (TH2F*)fInternalList->FindObject("hEtaPhiParticleMap");
  hEtaPhiParticleMap->Reset();
  
  for(Int_t i= 0; i < nTracks; i++) {
    particle = (AliMCParticle*) mcEvent->GetTrack(i);
    if(!particle) continue;
    
    if((stack->IsPhysicalPrimary(i)) && (particle->Charge() != 0)) {
      hEtaPhiParticleMap->Fill(particle->Eta(), particle->Phi());
      //std::cout << "Hans er nederen" << std::endl;
    }
  }
  ProjectAndMirror("MC");
  CalculateValues("MC");
}

//_____________________________________________________________________
//
//
// EOF
