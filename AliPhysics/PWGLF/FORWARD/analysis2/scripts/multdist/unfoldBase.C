#ifndef __CINT__
#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TList.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <THistPainter.h>
#include <TObject.h>
#include <TMath.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TVirtualFitter.h>
//#include <AliUnfolding.h>
//#include <AliUnfolding.cxx>
//#include <AliPWG0Helper.h>
//#include <AliMultiplicityCorrection.h>
#include <TMinuit.h>
#include <TBox.h>
#include <TGaxis.h>
#include "TRandom.h"
#include "/home/caz/ALICE/AliRoot/PWGLF/FORWARD/analysis2/AliForwardMultiplicityDistribution.h"
#include "/home/caz/ALICE/AliRoot/PWGLF/FORWARD/analysis2/AliForwardMultiplicityDistribution.cxx"
#include <TROOT.h>
#include <TSystem.h>
#include "TH1D.h"
#include "/home/caz/Desktop/RooUnfold-1.0.3/src/RooUnfoldResponse.h"
#include "/home/caz/Desktop/RooUnfold-1.0.3/src/RooUnfoldBayes.h"
#include "unfoldChi2Method.C"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldBinByBin.h"
#else
class RooUnfoldResponse;
class RooUnfoldBayes;
class RooUnfold;
class TF1;
class TFile;
class TH1D;
class TH1;
class TH1F;
class TH2D;
class TGraphErrors;
class TGaxis;
class AliUnfolding;
#endif
using namespace std;
enum Method {  kBayes, kSvd };
enum { kBla = 0, 
       kPol0, 
       kPol1, 
       kLog, 
       kEntropy, 
       kCurvature, 
       kRatio };

void    createOneUnfold(TList* list, TFile* dataFile, TFile* out, Double_t l, Double_t h, Method method);
TList*  getList(TFile* f, const Char_t* parent, const Char_t* name);
const char* getPostfix(const TH1* h);
void    unfoldBayesSet(const TObjArray& a, TDirectory* dir, TH2D* response, TH1D* triggerBias);
void    unfoldChi2minSet(const TObjArray& a, TDirectory* dir, TH2D* response, TH1D* triggerBias, Int_t limit);
void    normaliseAndCorrect(TH1D* unfolded, TH1D* triggerBias);
void    normaliseAndCorrect(TH1D* unfolded, TH1D* triggerBias, Int_t limit);

TH1D*   unfoldBayes(TH1D* data, RooUnfoldResponse* response);
TH1D*   unfoldChi2min(TH1F* data, TH2D* response, TH1F* eff, Int_t regFunc, Float_t regWeight);
TH1D*   getTriggerBiasHistogram(TList* responseList);
void    doTriggerBiasCorrection(TH1D*& hist, TH1D* triggerBias);

// --- Steering routine ----------------------------------------------
void unfoldBase(const Char_t* outputFile="unfoldOutput.root", Method method= kBayes,
		  const Char_t* responseFileName="$ALICE_ROOT/PWGLF/FORWARD/analysis2/responseMatrices.root", 
	  const Char_t* dataFileName="$ALICE_ROOT/PWGLF/FORWARD/analysis2/forward_multiplicity.root"){

  /*#ifdef __CINT__
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWG0selectors");
  
  gROOT->Macro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/scripts/LoadLibs.C");
  gROOT->LoadMacro("/home/caz/ppData/FirstAnalysis/multDists/Unfold.C");
#endif
  */
  TFile* responseFile = TFile::Open(responseFileName,"READ");
  if (!responseFile) { 
    Error("createUnfoldingFile", "Couldn't open %s", responseFileName);
    return;
  }
  TList* list = static_cast<TList*>(responseFile->Get("ResponseMatrices"));
  TFile* dataFile = TFile::Open(dataFileName,"READ");
  if (!dataFile) { 
    Error("createUnfoldingFile", "Couldn't open %s", dataFileName);
    return;
  }
  TFile* out = TFile::Open(outputFile, "RECREATE");
  cout << " Loaded responseFile and dataFile" << endl;

  
  Double_t limits[] = { 2.4, 1, 0. };
  // Double_t limits[] = {5.1, 3.4, 3.0, 2.4, 2.0, 1.5, 1.0, 0.5, 0. };
  Double_t* limit = limits;
  while ((*limit) > 0.1){
    if((*limit) >5)
      createOneUnfold(list, dataFile, out, -3.4,(*limit), method);
    else{
      createOneUnfold(list, dataFile, out, -(*limit),(*limit), method);
    }
    limit++; 
  }
  
  /*
  limit = limits;
  while ((*limit) > 0.1) { 
    createOneUnfold(list, dataFile, out, 0,(*limit),method);
    createOneUnfold(list, dataFile, out, -(*limit), 0, method);
    limit++; }
  
  
  for (Double_t l = -3; l < -0.5; l += 0.5){ 
    createOneUnfold(list, dataFile, out, l, l+0.5, method);  }
  for (Double_t l = 0.5; l < 5; l += 0.5){ 
    createOneUnfold(list, dataFile, out, l, l+0.5, method);  }
  for (Double_t l = -3; l < 5; l += 0.2){ 
    createOneUnfold(list, dataFile, out, l, l+0.2, method);  }
  */
    
  out->Close();
  responseFile->Close();
  dataFile->Close();
}


// --- Unfold data for one eta bin -----------------------------------
void createOneUnfold(TList* list, TFile* dataFile, TFile* out, Double_t l, Double_t h, Method method)
{
  const char* name= AliForwardMultiplicityDistribution::FormBinName(l, h);
  std::cout << "Name=" << name << std::endl;
  TDirectory* dir = out->mkdir(name);
  dir->cd();
  
  TList* multList = getList(dataFile, "CentralSums", name);
  if (!multList) return;
  
  //_____if unfolding MC response files____________________
  //TList* multList = getList(dataFile, "ResponseMatrices", name);
  //if (!multList) return;
  //cout << multList << endl;
    
  Int_t maxMult;
  TH1D* tmp= (TH1D*)multList->FindObject("mult");
  for(Int_t n=5; n<=tmp->GetNbinsX();n++){
    if(tmp->GetBinContent(n)<1e-9){
      maxMult= n;
      break;
    }
  }

  
  TObjArray a;
  a.Add(multList->FindObject("mult"));
  a.Add(multList->FindObject("multPlusSys"));
  a.Add(multList->FindObject("multMinusSys"));
  
  TList* responseList = 0;
  responseList = static_cast<TList*>(list->FindObject(name));
    
  TH2D* response= (TH2D*)responseList->FindObject("responseMatrix");
  response->SetDirectory(0);
  
  TH1D* triggerBias = getTriggerBiasHistogram(responseList);
  
  switch (method) { 
  case kBayes: 
    unfoldBayesSet(a, dir, response, triggerBias);
    break;
  case kSvd: 
     unfoldChi2minSet(a, dir, response, triggerBias, maxMult+30);
     break;
  }
  
  
  TIter next(&a);
  TH1D* data = 0;
  while ((data = static_cast<TH1D*>(next()))) {
    normaliseAndCorrect(data,triggerBias);
  }
  
  dir->Add(multList->FindObject("multMC"));
  
  dir->Add(&a);
  dir->Write();
  
  //delete list;
  //delete responseList;
  
  //delete triggerBias;
  delete multList;
  delete dir;
}

// --- Get a list from a file ----------------------------------------
TList*
getList(TFile* f, const Char_t* parent, const Char_t* name)
{
  TList* list = static_cast<TList*>(f->Get(parent));
  if (!list) { 
    Error("getList", "Couldn't find list %s in %s", parent, f->GetName());
    return 0;
  }
  if (!name) return list;
  TList* ret = static_cast<TList*>(list->FindObject(name));
  if (!ret) { 
    Error("getList", "Couldn't find %s in %s", name, parent);
    return 0;
  }
  
  return ret;
 
}

// --- Get end of input histogram name -------------------------------
const Char_t* 
getPostfix(const TH1* h)
{
  static TString t;
  t = h->GetName();
  t.ReplaceAll("mult", "");
  return t.Data();
}

// --- Unfold a single hist using Bayes ------------------------------
void
unfoldBayesSet(const TObjArray& a, TDirectory* dir, TH2D* response, TH1D* triggerBias)
{
  TH1* projY= (TH1D*) response->ProjectionY("projY",1,response->GetNbinsX(),"");
  TH1* projX= (TH1D*) response->ProjectionX("projX",1,response->GetNbinsY(),"");
  projX->SetDirectory(0);
  projY->SetDirectory(0);
  
  TH2D* responseTranspose = (TH2D*)response->Clone("response");
  for(Int_t i=1;i<=response->GetNbinsX();i++){
    for(Int_t j=1;j<=response->GetNbinsY();j++){
      responseTranspose->SetBinContent(j,i, response->GetBinContent(i,j));
      responseTranspose->SetBinError(j,i, response->GetBinError(i,j));
    }
  }
  RooUnfoldResponse* responseObject = new RooUnfoldResponse(projY, projX, responseTranspose,"test", "bla");

  Int_t dC = gStyle->GetNumberOfColors() / a.GetEntriesFast(); 
  Int_t iC = 0;
  TIter next(&a);
  TH1D* h = 0;
  while ((h = static_cast<TH1D*>(next()))) {
    TH1D* unfolded = unfoldBayes(h, responseObject);
    normaliseAndCorrect(unfolded,triggerBias);
    unfolded->SetDirectory(dir);
    unfolded->SetMarkerColor(gStyle->GetColorPalette(iC += dC));
    unfolded->SetMarkerStyle(20);
  }
  
  delete projY;
  delete projX;
  delete response;
  delete responseTranspose;
  delete responseObject;
}

// --- Unfold a single hist using Jan-Fietes SVD (Chi2) method --------------------------
void
unfoldChi2minSet(const TObjArray& a, TDirectory* dir, TH2D* response, TH1D* triggerBias, Int_t limit)
{

  TIter next(&a);
  TH1D* data = 0;
  while ((data = static_cast<TH1D*>(next()))) {
    Int_t max = (limit < 0 ? data->GetNbinsX() : limit);
    TH1F* dataTmp     = new TH1F(data->GetName(),data->GetTitle(), max,-0.5, max-.5);
    dataTmp->SetDirectory(0);
    for(Int_t k=1;k<= max; k++){
      dataTmp->SetBinContent(k,data->GetBinContent(k));
      dataTmp->SetBinError(k, data->GetBinError(k));
    }
    
    TH2D* responseTmp = response; 
    TH1F* eff         = new TH1F("eff","unfolded",max, -.5, max-.5);
    if (limit > 0) {
      responseTmp = new TH2D("res","res", max, -.5, max-.5,max, -.5, max-.5);
      for(Int_t k=1;k<= max; k++){
	eff->SetBinContent(k,1);
	if (limit > 0) {
	  for(Int_t j=1;j<=max; j++){
	    responseTmp->SetBinContent(k,j,response->GetBinContent(k,j));
	    responseTmp->SetBinError(k,j, response->GetBinError(k,j));
	  }
	}
      }    
    }
    Int_t nR = 8;
    Int_t dC = gStyle->GetNumberOfColors() / (nR * a.GetEntriesFast()); 
    Int_t iC = 0;
    for(Int_t f = 3; f <= 3; f++){
      iC = 0;
      for(Float_t w = 1e-5; w <= 1e5 ;w *= 10){
	TH1D* unfolded = unfoldChi2min(dataTmp,responseTmp, eff, f, w);
	normaliseAndCorrect(unfolded,triggerBias, limit);
	unfolded->SetDirectory(dir);
	unfolded->SetMarkerColor(gStyle->GetColorPalette(iC += dC));
	unfolded->SetMarkerStyle(20);
      }
    }
    if (limit > 0) delete responseTmp;
    delete dataTmp;
    delete eff;
  }
}

// --- Normalise to 1 and correct for trigger ------------------------
void normaliseAndCorrect(TH1D* unfolded, TH1D* triggerBias)
{
  if (triggerBias) doTriggerBiasCorrection(unfolded, triggerBias);
  unfolded->SetDirectory(0);
  unfolded->Sumw2();
  unfolded->Scale(1/unfolded->Integral());
}

void normaliseAndCorrect(TH1D* unfolded, TH1D* triggerBias, Int_t limit)
{
  if (triggerBias) doTriggerBiasCorrection(unfolded, triggerBias);
  unfolded->SetDirectory(0);
  unfolded->Sumw2();
  unfolded->Scale(1/unfolded->Integral(1,limit));
}

// --- Use RooUnfold stuff -------------------------------------------
TH1D* unfoldBayes(TH1D* data, RooUnfoldResponse* response)
{
  
  
  RooUnfold* unfold = new RooUnfoldBayes(response, data, 20);
  TH1D* unfolded= (TH1D*) unfold->Hreco();
  
  TString t = TString::Format("unfolded_bayes_%s",getPostfix(data));
  
  unfolded->SetName(t.Data());
  delete unfold;

  return unfolded;
}

// --- Use Jan-Fiete routines ----------------------------------------
TH1D* unfoldChi2min(TH1F* data, TH2D* response, TH1F* eff, 
		    Int_t regFunc, Float_t regWeight)
{  
  TH1F* tmp = static_cast<TH1F*>(data->Clone("tmp"));
  tmp->Reset();
  UnfoldChi2Min(data, tmp, response, eff, regFunc, regWeight);  
  
  Int_t max = data->GetNbinsX();
  TH1D* unfolded= new TH1D("unfolded","unfolded", max, -.5, max);  
  for(Int_t n=1;n<=max; n++){
    unfolded->SetBinContent(n,tmp->GetBinContent(n));
    unfolded->SetBinError(n, tmp->GetBinError(n));
  }

  TString regString="0";
  switch (regFunc) { 
  case 1: regString="pol0";      break;
  case 2: regString="pol1";      break;
  case 3: regString="log";       break;
  case 4: regString="entropy";   break;
  case 5: regString="curvature"; break;
  default:
    Error("unfold", "Invalid regFunc type");
    return 0;
  }
  delete tmp; 

  Info("", "data name: %s, postfix: %s", data->GetName(), getPostfix(data));
  TString n = TString::Format("unfolded_chi2_%s_%s_%1.0e", getPostfix(data),regString.Data(),regWeight);
  n.ReplaceAll("-0", "m");
  n.ReplaceAll("+0", "p");
  //n.ReplaceAll(".", "d");
  unfolded->SetName(n);
  return unfolded;
}


// --- Find trigger bias histogram -----------------------------------
TH1D* getTriggerBiasHistogram(TList* responseList){
  TH1D* hMCNSD;
  TH1D* hMCESDNSD;
  TH1D* hESDNSD;
 
  hMCNSD = (TH1D*)responseList->FindObject("fMCNSD");
  hMCESDNSD = (TH1D*)responseList->FindObject("fMCESDNSD");
  hESDNSD = (TH1D*)responseList->FindObject("fESDNSD");
    
  hMCNSD->Sumw2();
  hESDNSD->Sumw2();
  hMCESDNSD->Sumw2();
  
  TH1D* corr     = new TH1D("corr", "corr", 50,-0.5,49.5);
  for(Int_t n=1;n<corr->GetNbinsX();n++){
    Double_t errorSquared=0;
    errorSquared= (1/hESDNSD->GetBinContent(n)+   1/hMCNSD->GetBinContent(n));   
    corr->SetBinContent(n, hESDNSD->GetBinContent(n)/hMCNSD->GetBinContent(n));
    corr->SetBinError(n,corr->GetBinContent(n)*TMath::Sqrt(errorSquared));
  }
  
  hMCESDNSD->Divide(hMCNSD);
  hESDNSD->Divide(hMCNSD);
   
  
  delete hMCNSD;
  delete hMCESDNSD;
  
  return corr;
    
}

// --- Apply trigger bias histogram ----------------------------------
void doTriggerBiasCorrection(TH1D*& hist, TH1D* triggerBias){
  for(Int_t i = 1; i<=35;i++){
    if(triggerBias->GetBinContent(i)>1e-5&&hist->GetBinContent(i)>0){
      hist->SetBinContent(i, hist->GetBinContent(i)/triggerBias->GetBinContent(i));
      hist->SetBinError(i,TMath::Sqrt(TMath::Power(hist->GetBinError(i)/hist->GetBinContent(i),2)+TMath::Power(triggerBias->GetBinError(i)/triggerBias->GetBinContent(i),2))*hist->GetBinContent(i));
      //cout << hist->GetBinContent(i) << " +-  " << hist->GetBinError(i) << endl; 
    }
    //if(triggerBias->GetBinContent(i)<1e-5)
    //  hist->SetBinContent(i, hist->GetBinContent(i)/triggerBias->GetBinContent(i+1));
    
  }
}

//
// EOF
//
