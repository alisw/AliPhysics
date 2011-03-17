#include "AliOADBPWG2Spectra.h"
#include "TNamed.h"
#include "TString.h"
#include "TH1D.h"
#include "TObject.h"
#include "TList.h"
#include "AliAnalysisManager.h"
#include "TBrowser.h"
#include "AliLog.h"

ClassImp(AliOADBPWG2Spectra)

const char * AliOADBPWG2Spectra::fDetectorNames[] = {"ITS", "ITSTPC", "TPC", "TOF", "TOFTPC"};
const char * AliOADBPWG2Spectra::fPidTypeNames[]  = {"GaussFit", "NSigma", "Bayes", "Kinks"};
const char * AliOADBPWG2Spectra::fChargeTags[]    = {"Pos", "Neg"};
const char * AliOADBPWG2Spectra::fParticleNames[] = {"Pion", "Kaon", "Proton"};


AliOADBPWG2Spectra::AliOADBPWG2Spectra():
TNamed("Dummy", "OADB Object for PWG2 Spectra" ), fHistos(0)
{
  // ctor
  


}
AliOADBPWG2Spectra::AliOADBPWG2Spectra(char* name) :
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
  filename.Form("%s/PWG2/SPECTRA/spectraResults.root", AliAnalysisManager::GetOADBPath()); 
  return filename.Data();

}
const char * AliOADBPWG2Spectra::GetHistoName(EPWG2SpectraDetector det, EPWG2SpectraPIDType pidType, EPWG2SpectraParticle part, 
						     EPWG2SpectraCharge charge, const char * centrTag, Int_t centrBin) {

  // Returns histogram name
  // h[Name]_[Detector(s)]_[PIDType]_[Particle]_[Pos|Neg]_[MultiplicityOrCentralityIndex]
  // where "name" is the name of this container

  static TString histoName;
  if (centrTag)
    histoName.Form("h%s_%s_%s_%s_%s_%s_%d", GetName(), fDetectorNames[det], fPidTypeNames[pidType], fParticleNames[part], fChargeTags[charge], centrTag, centrBin);
  else 
    histoName.Form("h%s_%s_%s_%s_%s",       GetName(), fDetectorNames[det], fPidTypeNames[pidType], fParticleNames[part], fChargeTags[charge]);

  return histoName.Data();
}

TH1D * AliOADBPWG2Spectra::GetHisto(EPWG2SpectraDetector det, EPWG2SpectraPIDType pidType, EPWG2SpectraParticle part, 
				    EPWG2SpectraCharge charge, const char * centrTag, Int_t centrBin){

  // Get an histogram from the list
  const char * name = GetHistoName(det,pidType,part,charge,centrTag,centrBin);
  TH1D * h = (TH1D*) fHistos->FindObject(name);
  return h;

}

void  AliOADBPWG2Spectra::AddHisto(TH1D * h, EPWG2SpectraDetector det, EPWG2SpectraPIDType pidType, EPWG2SpectraParticle part, 
				    EPWG2SpectraCharge charge, const char * centrTag, Int_t centrBin) {
  // Add a histogram to the list
  // Rename it if necessary
  const char * name = GetHistoName(det,pidType,part,charge,centrTag,centrBin);
  TH1D * hold = (TH1D*) fHistos->FindObject(name);
  if (hold) fHistos->Remove(hold);
  delete hold;
  if(strcmp(h->GetName(),name)){
    AliError(Form("Histogram namws are not consinstent %s-%s, resetting", h->GetName(),name));
    h->SetName(name); 
  }
  fHistos->Add(h);


}

TH1D * AliOADBPWG2Spectra::BookHisto(EPWG2SpectraDetector det, EPWG2SpectraPIDType pidType, EPWG2SpectraParticle part, 
				     EPWG2SpectraCharge charge, const char * centrTag, Int_t centrBin) {

  // book a histogram according to the template. All the histograms
  // should have the same binning (if possible/reasonable) to
  // facilitate the compoarison and the combination

  const Float_t templBins[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2,2.2,2.4,2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
  Int_t nbinsTempl=48;
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
