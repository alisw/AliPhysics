#include "AliOADBPWG2Spectra.h"
#include "AliOADBContainer.h"
#include "TNamed.h"
#include "TString.h"
#include "TH1D.h"
#include "TObject.h"
#include "TList.h"
#include "AliAnalysisManager.h"
#include "TBrowser.h"
#include "AliLog.h"

ClassImp(AliOADBPWG2Spectra)

const char * AliOADBPWG2Spectra::fgkDetectorNames[] = {"ITS", "ITSTPC", "TPC", "TOF", "TOFTPC", "Dummy", "Dummy"};
const char * AliOADBPWG2Spectra::fgkPidTypeNames[]  = {"GaussFit", "NSigma", "Bayes", "Kinks"};
const char * AliOADBPWG2Spectra::fgkChargeTags[]    = {"Pos", "Neg"};
const char * AliOADBPWG2Spectra::fgkParticleNames[] = {"Pion", "Kaon", "Proton"};


AliOADBPWG2Spectra::AliOADBPWG2Spectra():
TNamed("Dummy", "OADB Object for PWG2 Spectra" ), fHistos(0)
{
  // ctor
  


}
AliOADBPWG2Spectra::AliOADBPWG2Spectra(const char* name) :
TNamed(name, "OADB Object for PWG2 Spectra" ), fHistos(0) 

{
  // ctor
  // name is appended to all histos (e.g. "Corrected")

  Init();

}

AliOADBPWG2Spectra::~AliOADBPWG2Spectra() {
  // dtor
  if(fHistos) delete fHistos;
}

void AliOADBPWG2Spectra::Init() {
  fHistos = new TList();
}


const char * AliOADBPWG2Spectra::GetOADBPWG2SpectraFileName()  {
  // get file name to the OADB
  static TString filename;
  filename.Form("%s/PWG2/SPECTRA/spectraResults.root", AliOADBContainer::GetOADBPath()); 
  return filename.Data();

}
const char * AliOADBPWG2Spectra::GetHistoName(Int_t det, Int_t pidType, Int_t part, 
						     Int_t charge, const char * centrTag, Int_t centrBin) {

  // Returns histogram name
  // h[Name]_[Detector(s)]_[PIDType]_[Particle]_[Pos|Neg]_[MultiplicityOrCentralityIndex]
  // where "name" is the name of this container

  
  static TString histoName;
  if (centrTag) {
    if(!strcmp(centrTag,"MB")){
      // don't put a index for MB spectra
      histoName.Form("h%s_%s_%s_%s_%s_%s", GetName(), fgkDetectorNames[det], fgkPidTypeNames[pidType], fgkParticleNames[part], fgkChargeTags[charge], centrTag);
    }
    else {
      histoName.Form("h%s_%s_%s_%s_%s_%s_%d", GetName(), fgkDetectorNames[det], fgkPidTypeNames[pidType], fgkParticleNames[part], fgkChargeTags[charge], centrTag, centrBin);
    }
  }
  else 
    histoName.Form("h%s_%s_%s_%s_%s",       GetName(), fgkDetectorNames[det], fgkPidTypeNames[pidType], fgkParticleNames[part], fgkChargeTags[charge]);

  return histoName.Data();
}

TH1D * AliOADBPWG2Spectra::GetHisto(Int_t det, Int_t pidType, Int_t part, 
				    Int_t charge, const char * centrTag, Int_t centrBin){

  // Get an histogram from the list
  const char * name = GetHistoName(det,pidType,part,charge,centrTag,centrBin);
  TH1D * h = (TH1D*) fHistos->FindObject(name);
  return h;

}

void  AliOADBPWG2Spectra::AddHisto(TH1D * h, Int_t det, Int_t pidType, Int_t part, 
				    Int_t charge, const char * centrTag, Int_t centrBin) {
  // Add a histogram to the list
  // Rename and rebinn it if necessary
  
  if(!h) {
    AliWarning("Empty pointer to histogram");
    return;
  }
  
  static TH1D * htest = BookHisto(kDetDummy, kGaussFit,kPion, kPos);
  const char * name = GetHistoName(det,pidType,part,charge,centrTag,centrBin);
  if(!CompareBinning(h,htest)){
    AliWarning("Histo have different binning! Rebinning to standard"){
      h = GetHistoStandardBinning(h,det,pidType,part,charge,centrTag,centrBin);
    }
  }
  if(!fHistos) {
    AliError("fHistos not allocated!!");
    return;
  }
    
  TH1D * hold = (TH1D*) fHistos->FindObject(name);
  if (hold) fHistos->Remove(hold);
  delete hold;
  if(strcmp(h->GetName(),name)){
    AliWarning(Form("Histogram names are not consinstent %s-%s, resetting", h->GetName(),name));
    h->SetName(name); 
  }
  fHistos->Add(h);


}

TH1D * AliOADBPWG2Spectra::BookHisto(Int_t det, Int_t pidType, Int_t part, 
				     Int_t charge, const char * centrTag, Int_t centrBin) {

  // book a histogram according to the template. All the histograms
  // should have the same binning (if possible/reasonable) to
  // facilitate the compoarison and the combination

  const Float_t templBins[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
  Int_t nbinsTempl=52;
  const char * name = GetHistoName(det,pidType,part,charge,centrTag,centrBin);
  TH1D * h = new TH1D(name, name, nbinsTempl, templBins);
  return h;

}

void AliOADBPWG2Spectra::Browse(TBrowser *b)
{
  // Browse this object.
   // If b=0, there is no Browse call TObject::Browse(0) instead.
   //         This means TObject::Inspect() will be invoked indirectly


  if (b) {
    b->Add(fHistos);        
  }     
   else
      TObject::Browse(b);
}

TH1D * AliOADBPWG2Spectra::GetHistoStandardBinning(const TH1D* h, Int_t det, Int_t pidType, Int_t part, 
						   Int_t charge, const char * centrTag, Int_t centrBin) {
  // Returns a histo with the standard binning and the same content of h
  // if the bins of h are not a subset of the standard binning, it crashes with a fatal error
  // under and overflows are ignored
  
  // 1. Create a histogram with the standard binning
  TH1D * hStd = BookHisto(det,  pidType,  part, charge, centrTag, centrBin);
  Int_t nBinsH1=hStd->GetNbinsX();
  Int_t nBinsH2=h->GetNbinsX();
  // Loop over standard bins, 
  for(Int_t iBin=1; iBin<=nBinsH1; iBin++){
    Float_t lowPtH1 =hStd->GetBinLowEdge(iBin);
    Float_t binWidH1=hStd->GetBinWidth(iBin);
    // Loop over H2 bins and find overlapping bins to H1
    for(Int_t jBin=1; jBin<=nBinsH2; jBin++){
      Float_t lowPtH2=h->GetBinLowEdge(jBin);
      Float_t binWidH2=h->GetBinWidth(jBin);
      if(TMath::Abs(lowPtH1-lowPtH2)<0.001 && TMath::Abs(binWidH2-binWidH1)<0.001){
	hStd->SetBinContent(iBin, h->GetBinContent(jBin));
	hStd->SetBinError  (iBin, h->GetBinError  (jBin));
	break;
      }
      if(TMath::Abs(lowPtH1-lowPtH2)<0.001){
	AliError(Form("Found partially overlapping bins! [(%f,%f)(%f,%f)]",lowPtH1,binWidH1,lowPtH2,binWidH2));
	continue;
      }
    }
  }
  return hStd;
}

Bool_t AliOADBPWG2Spectra::CompareBinning(TH1 * h1, TH1 * h2) {

  // returns true if h1 and h2 have the same binning
  Int_t nbins1 = h1->GetNbinsX();
  Int_t nbins2 = h2->GetNbinsX();
  
  if(nbins1 != nbins2) return kFALSE;
  
  for(Int_t ibin = 1; ibin <= nbins1; ibin++){
    if(TMath::Abs(h1->GetBinLowEdge(ibin) - h2->GetBinLowEdge(ibin))>0.001) return kFALSE;
    if(TMath::Abs(h1->GetBinWidth(ibin) - h2->GetBinWidth(ibin))>0.001) return kFALSE;
  }
  
  return kTRUE;
}
