/**
 * @file   UnfoldMult.C
 * @author Valentina Zaccolo
 * @date   Tue Mar 05 12:56:56 2013
 * 
 * @brief  Bayesian Method Unfolding Script
 *
 * @ingroup pwglf_forward_scripts
 * @ingroup pwglf_forward_multdist
 */
#include <iostream>
#include <TList.h>
#include <TFile.h>
#include <TObject.h>
#include <TDirectory.h>
#include "TH1D.h"
#include "TH2D.h"
#include <TROOT.h>
#include <TSystem.h>
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include <TMath.h>
/** 
 * @defgroup pwglf_forward_scripts_unfoldmult UnfoldMult stuff 
 *
 * @ingroup pwglf_forward_multdist
 */
/** 
 * Get an object from a directory 
 * 
 * @param d     Directory 
 * @param name  Name of object 
 * 
 * @return Pointer to object or null if not found 
 *
 * @ingroup pwglf_forward_scripts_unfoldmult
 */
TObject* GetObject(TDirectory* d, const char* name)
{
  TObject* o = d->Get(name);
  if (!o) { 
    Error("GetCollection", "%s not found in %s", name, d->GetName());
    return 0;
  }
  return o;
}
/** 
 * Get an object from a collection
 * 
 * @param d     Collection 
 * @param name  Name of object
 * 
 * @return Pointer to object or null if not found 
 *
 * @ingroup pwglf_forward_scripts_unfoldmult
 */
TObject* GetObject(TCollection* d, const char* name)
{
  TObject* o = d->FindObject(name);
  if (!o) { 
    Error("GetCollection", "%s not found in %s", name, d->GetName());
    d->ls();
    return 0;
  }
  return o;
}
/** 
 * Get a collection from a directory
 * 
 * @param d     Parent directory 
 * @param name  Name of collection
 * 
 * @return Pointer to collection or null 
 *
 * @ingroup pwglf_forward_scripts_unfoldmult
 */
TCollection* GetCollection(TDirectory* d, const char* name)
{
  return static_cast<TCollection*>(GetObject(d, name));
}
/** 
 * Get a collection from a collection
 * 
 * @param d     Parent collection
 * @param name  Name of collection
 * 
 * @return Pointer to collection or null 
 *
 * @ingroup pwglf_forward_scripts_unfoldmult
 */
TCollection* GetCollection(TCollection* d, const char* name)
{
  return static_cast<TCollection*>(GetObject(d, name));
}
/** 
 * Get a 1D-histogram from a directory 
 * 
 * @param d     Parent directory 
 * @param name  Name of histogram 
 * 
 * @return Pointer to histogram or null 
 *
 * @ingroup pwglf_forward_scripts_unfoldmult
 */
TH1* GetH1(TDirectory* d, const char* name)
{
  return static_cast<TH1*>(GetObject(d, name));
}
/** 
 * Get a 1D-histogram from a collection 
 * 
 * @param d     Parent collectoin 
 * @param name  Name of histogram 
 * 
 * @return Pointer to histogram or null 
 *
 * @ingroup pwglf_forward_scripts_unfoldmult
 */
TH1* GetH1(TCollection* d, const char* name)
{
  return static_cast<TH1*>(GetObject(d, name));
}
/** 
 * Get a 2D-histogram from a directory 
 * 
 * @param d     Parent directory 
 * @param name  Name of histogram 
 * 
 * @return Pointer to histogram or null 
 */
TH2* GetH2(TDirectory* d, const char* name)
{
  return static_cast<TH2*>(GetObject(d, name));
}
/** 
 * Get a 2D-histogram from a collection
 * 
 * @param d     Parent collection
 * @param name  Name of histogram 
 * 
 * @return Pointer to histogram or null 
 *
 * @ingroup pwglf_forward_scripts_unfoldmult
 */
TH2* GetH2(TCollection* d, const char* name)
{
  return static_cast<TH2*>(GetObject(d, name));
}
/** 
 * Open a file 
 * 
 * @param filename      File name 
 * @param readNotWrite  If true, open read-only
 * 
 * @return Pointer to file, or null
 *
 * @ingroup pwglf_forward_scripts_unfoldmult
 */
TFile* OpenFile(const char* filename, Bool_t readNotWrite)
{
  TFile* ret = TFile::Open(filename, readNotWrite ? "READ" : "RECREATE");
  if (!ret) { 
    Error("OpenFile", "Failed to open the file \"%s\"", filename);
    return 0;
  }
  return ret;
}
/** 
 * A type definition to make it easier to switch to TDirectory structure 
 *
 * @ingroup pwglf_forward_scripts_unfoldmult
 */
typedef TCollection Dir;

/** 
 * Get histograms for @f$\eta@f$-bin from @a lowLim to @a highLim and
 * process it .
 * 
 * @param lowLim   Lower @f$\eta@f$ limit
 * @param highLim  Upper @f$\eta@f$ limit 
 * @param resp     Directory containing response matrices 
 * @param data     Directory containing data histogram 
 * @param out      Output directory  
 *
 * @ingroup pwglf_forward_scripts_unfoldmult
 */
void GetHistos(Double_t    lowLim, 
	       Double_t    highLim, 
	       Dir*        resp, 
	       Dir*        data, 
	       TDirectory* out);
/** 
 * Do the actual unfolding and corrections 
 * 
 * @param response   Response matrix 
 * @param measured   Measured distribution
 * @param mcHist     MC truth distribution
 * @param esdHist    MC distribution after ESD event selection
 *
 * @ingroup pwglf_forward_scripts_unfoldmult
 */
void ProcessUnfold(TH2* response, TH1* measured, TH1* mcHist, TH1* esdHist); 
/** 
 * Caclulate the trigger bias correction 
 * 
 * @param mcHist     MC truth distribution
 * @param esdHist    MC distribution after ESD event selection
 * 
 * @return  Histogram of trigger bias correction
 *
 * @ingroup pwglf_forward_scripts_unfoldmult
 */
TH1* GetTriggerBias(TH1* mcHist, TH1* esdHist);
/** 
 * Apply the trigger bias correction to the unfolded result
 * 
 * @param unfoldBefore Unfolded measured distribution
 * @param triggerBias  Trigger bias correction 
 * @param mcHist       MC truth distribution
 */
void DoCorrection(TH1* unfoldBefore, TH1* triggerBias, TH1* mcHist);


/** 
 * Main Loop to call eta range
 * 
 * @param respFile File containing response matrices 
 * @param dataFile File containing measured distributions
 * @param outFile  Output file 
 *
 * @ingroup pwglf_forward_scripts_unfoldmult
 */
void UnfoldMult(const char* respFile="forward_response.root", 
		const char* dataFile="forward_multiplicy.root", 
		const char* outFile="unfolded.root") 
{
  TFile* output = OpenFile(outFile,false);
  TFile* resp   = OpenFile(respFile, true);
  TFile* data   = OpenFile(dataFile, true);
  if (!output || !resp || !data) return;

  Dir* topResp  = GetCollection(resp, "ResponseMatricesSums");
  Dir* topData  = GetCollection(data, "MultSums");
  if (!topData || !topResp) return;
  
  Double_t limits[] = {5.1, 3.4, 3.0, 2.4, 2.0, 1.5, 1.0, 0.5, 0. };
  Double_t* limit = limits;

  while ((*limit) > 0.1) {
    if ((*limit) >5.) GetHistos(-3.4,      (*limit), topResp, topData, output);
    else              GetHistos(-(*limit), (*limit), topResp, topData, output);
    limit++; 
    output->cd();
  }
  output->Close();
  resp->Close();
  data->Close();
}

//____________________________________________________________________
void GetHistos(Double_t    lowLim, 
	       Double_t    highLim, 
	       Dir*        topResp, 
	       Dir*        topData,
	       TDirectory* out) 
{
  // Format the name of the bin 
  TString sLow( TString::Format("%+03d",int(10*lowLim)));
  TString sHigh(TString::Format("%+03d",int(10*highLim)));
  sLow .ReplaceAll("+", "plus");
  sLow .ReplaceAll("-", "minus");
  sHigh.ReplaceAll("+", "plus");
  sHigh.ReplaceAll("-", "minus");
  TString     name = TString::Format("%s_%s", sLow.Data(), sHigh.Data());
  TDirectory* dir = out->mkdir(name);
  dir->cd();

  // Get our collections 
  Dir* listResp = GetCollection(topResp, name);
  Dir* listData = GetCollection(topData,name);
  if(!listResp || !listData) return;

  // Get the MC histograms 
  TH2* response = GetH2(listResp, "responseMatrix");
  TH1* mcHist   = GetH1(listResp, "fMCNSD");
  TH1* esdHist  = GetH1(listResp, "fESDNSD");
  if (!response || !mcHist || !esdHist) return;

  // Get the data histogram 
  TH1* measured = GetH1(listData, "mult");
  if(!measured) return; 

  // Now do the unfolding 
  ProcessUnfold(response, measured, mcHist, esdHist);
} 



//____________________________________________________________________
void ProcessUnfold(TH2* response, TH1* measured, TH1* mcHist, TH1* esdHist) 
{
  Int_t nX   = response->GetNbinsX();
  Int_t nY   = response->GetNbinsY();
  TH1* projY = static_cast<TH1*>(response->ProjectionY("projY",1,nX,""));
  TH1* projX = static_cast<TH1*>(response->ProjectionX("projX",1,nY,""));
  projX->SetDirectory(0);
  projY->SetDirectory(0);
  
  TH2* responseTranspose = static_cast<TH2*>(response->Clone("response"));
  for(Int_t i = 1; i <= nX; i++) {
    for(Int_t j = 1; j <= nY; j++) {
      responseTranspose->SetBinContent(j,i, response->GetBinContent(i,j));
      responseTranspose->SetBinError(j,i, response->GetBinError(i,j));
    }
  }
  
  RooUnfoldResponse* responseObject = new RooUnfoldResponse(projY, projX, 
							    responseTranspose, 
							    "", "");
  RooUnfold*         unfold         = new RooUnfoldBayes(responseObject, 
							 measured, 5); 
  TH1*               unfolded       = static_cast<TH1*>(unfold->Hreco());
  TH1*               triggerBias    = GetTriggerBias(mcHist, esdHist);
  DoCorrection(unfolded, TriggerBias, mcHist);

  Double_t scale_unf = 1/unfolded->Integral();
  unfolded->Scale(scale_unf);
  unfolded->Write("Unfolded");

  Double_t scale_raw = 1/measured->Integral();
  measured->Scale(scale_raw);
  measured->Write("Raw");
  triggerBias->Write("TriggerBiasCorr");
}



//____________________________________________________________________
TH1* GetTriggerBias(TH1* mcHist, TH1* esdHist) 
{
  mcHist->Sumw2();
  esdHist->Sumw2();
  
  TH1* corr = new TH1D("corr", "corr", 40, -0.5, 39.5);  
  for (Int_t j=1; j<corr->GetNbinsX(); j++) {
    Double_t errSqr = 0.;
    Double_t cESD   = esdHist->GetBinContent(j);
    Double_t cMS    = mcHist->GetBinContent(j);
    Double_t c      = (cMC  == 0 ? 1 : 
		       cESD == 0 ? 1/cMC : cESD/cMC);
    corr->SetBinContent(j, c);
    corr->SetBinError(j, c*TMath::Sqrt(errSqr));
  }
  return corr;  
}

//____________________________________________________________________
void DoCorrection(TH1* unfoldBefore, TH1* triggerBias, TH1*) 
{
  for (Int_t k=1; k<=35; k++) {
    Double_t tb = triggerBias->GetBinContent(k);
    Double_t ub = unfoldBefore->GetBinContent(k);
    if (tb > 1e-5 && ub > 0) {
      unfoldBefore->SetBinContent(k, ub / tb);
      Double_t eub = UnfoldBefore->GetBinError(k);
      Double_t etb = TriggerBias->GetBinError(k);
      unfoldBefore->SetBinError(k, TMath::Sqrt(TMath::Power(eub/ub)+
					       TMath::Power(etb/tb))*ub/tb); 
    }
  }
}
//
// EOF
//



