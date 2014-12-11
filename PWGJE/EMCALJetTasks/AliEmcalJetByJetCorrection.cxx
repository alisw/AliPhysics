//
// Utility class to do jet-by-jet correction
// Based on templates of ptJet vs ptTrack vs r
//
// Author: M.Verweij

#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliVParticle.h"

#include "AliEmcalJetByJetCorrection.h"

ClassImp(AliEmcalJetByJetCorrection)

//__________________________________________________________________________
AliEmcalJetByJetCorrection::AliEmcalJetByJetCorrection() :
TNamed(),
  fh3JetPtDRTrackPt(0x0),
  fBinWidthJetPt(10.),
  fJetPtMin(0.),
  fJetPtMax(150.),
  fCollTemplates(),
  fInitialized(kFALSE),
  fEfficiencyFixed(1.),
  fhEfficiency(0)
{
  // Dummy constructor.
  fCollTemplates.SetOwner(kTRUE);
}

//__________________________________________________________________________
AliEmcalJetByJetCorrection::AliEmcalJetByJetCorrection(const char* name) :
  TNamed(name, name),
  fh3JetPtDRTrackPt(0x0),
  fBinWidthJetPt(10.),
  fJetPtMin(0.),
  fJetPtMax(150.),
  fCollTemplates(),
  fInitialized(kFALSE),
  fEfficiencyFixed(1.),
  fhEfficiency(0)
{
  // Default constructor.
  fCollTemplates.SetOwner(kTRUE);
}

//__________________________________________________________________________
AliEmcalJetByJetCorrection::AliEmcalJetByJetCorrection(const AliEmcalJetByJetCorrection &other) :
  TNamed(other),
  fh3JetPtDRTrackPt(other.fh3JetPtDRTrackPt),
  fBinWidthJetPt(other.fBinWidthJetPt),
  fJetPtMin(other.fJetPtMin),
  fJetPtMax(other.fJetPtMax),
  fCollTemplates(other.fCollTemplates),
  fInitialized(other.fInitialized),
  fEfficiencyFixed(other.fEfficiencyFixed),
  fhEfficiency(other.fhEfficiency)
{
  // Copy constructor.
}

//__________________________________________________________________________
AliEmcalJetByJetCorrection& AliEmcalJetByJetCorrection::operator=(const AliEmcalJetByJetCorrection &other)
{
  // Assignment
  if (&other == this) return *this;
  TNamed::operator=(other);
  fh3JetPtDRTrackPt = other.fh3JetPtDRTrackPt;
  fBinWidthJetPt     = other.fBinWidthJetPt;
  fJetPtMin          = other.fJetPtMin;
  fJetPtMax          = other.fJetPtMax;
  fCollTemplates     = other.fCollTemplates;
  fInitialized       = other.fInitialized;
  fEfficiencyFixed   = other.fEfficiencyFixed;
  fhEfficiency       = other.fhEfficiency;

  return *this;
}

//__________________________________________________________________________
AliEmcalJet* AliEmcalJetByJetCorrection::Eval(const AliEmcalJet *jet, TClonesArray *fTracks) {

  if(!fInitialized) {
    Printf("AliEmcalJetByJetCorrection %s not initialized",GetName());
    return NULL;
  }

  Int_t bin = GetJetPtBin(jet->Pt());
  if(bin<0 || bin>fCollTemplates.GetEntriesFast()) return NULL;

  TH2D *hTemplate = static_cast<TH2D*>(fCollTemplates.At(bin));
  
  Double_t meanPt = GetMeanPtConstituents(jet,fTracks);
  Double_t eff = GetEfficiency(meanPt);

  Int_t np = TMath::FloorNint((double)jet->GetNumberOfTracks() * (1./eff -1.));

  TLorentzVector corrVec; corrVec.SetPtEtaPhiM(jet->Pt(),jet->Eta(),jet->Phi(),jet->M());

  Double_t mass = 0.13957; //pion mass

  for(Int_t i = 0; i<np; i++) {
    Double_t r;
    Double_t pt;
    hTemplate->GetRandom2(r,pt);
    Double_t t = TMath::TwoPi()*gRandom->Uniform(1.);
    Double_t deta = r*TMath::Cos(t);
    Double_t dphi = r*TMath::Sin(t);
    TLorentzVector curVec; 
    curVec.SetPtEtaPhiM(pt,deta+jet->Eta(),dphi+jet->Phi(),mass);
    corrVec+=curVec;
  }

  AliEmcalJet *jetCorr = new AliEmcalJet(corrVec.Pt(),corrVec.Eta(),corrVec.Phi(),corrVec.M());

  return jetCorr;
}

//__________________________________________________________________________
void AliEmcalJetByJetCorrection::Init() {
  //Init templates

  if(!fh3JetPtDRTrackPt) {
    Printf("%s fh3JetPtDRTrackPt not known",GetName());
    fInitialized = kFALSE;
    return;
  }

  if(fJetPtMax>fh3JetPtDRTrackPt->GetXaxis()->GetXmax())
    fJetPtMax = fh3JetPtDRTrackPt->GetXaxis()->GetXmax();

  if(fJetPtMin<fh3JetPtDRTrackPt->GetXaxis()->GetXmin())
    fJetPtMin = fh3JetPtDRTrackPt->GetXaxis()->GetXmin();

  Double_t eps = 0.00001;
  Int_t counter = 0;
  for(Double_t ptmin = fJetPtMin; ptmin<fJetPtMax; ptmin+=fBinWidthJetPt) {
    Int_t binMin = fh3JetPtDRTrackPt->GetXaxis()->FindBin(ptmin+eps);
    Int_t binMax = fh3JetPtDRTrackPt->GetXaxis()->FindBin(ptmin+fBinWidthJetPt-eps);
    //    Printf("%d bins: %d - %d -> %f - %f",counter,binMin,binMax,fh3JetPtDRTrackPt->GetXaxis()->GetBinLowEdge(binMin),fh3JetPtDRTrackPt->GetXaxis()->GetBinUpEdge(binMax));

    fh3JetPtDRTrackPt->GetXaxis()->SetRange(binMin,binMax);
    TH2D *h2 = dynamic_cast<TH2D*>(fh3JetPtDRTrackPt->Project3D("zy"));
    h2->SetName(Form("hPtR_%.0f_%.0f",fh3JetPtDRTrackPt->GetXaxis()->GetBinLowEdge(binMin),fh3JetPtDRTrackPt->GetXaxis()->GetBinUpEdge(binMax)));
    fCollTemplates.Add(h2);
    counter++;
  }
  //  Int_t nt = TMath::FloorNint((fJetPtMax-fJetPtMin)/fBinWidthJetPt);
  //  Printf("nt: %d entries fCollTemplates: %d",nt,fCollTemplates.GetEntriesFast());

  fInitialized = kTRUE;

}

//__________________________________________________________________________
Int_t AliEmcalJetByJetCorrection::GetJetPtBin(const Double_t jetpt) const {

  if(jetpt<fJetPtMin || jetpt>=fJetPtMax)
    return -1;

  Int_t bin = TMath::FloorNint((jetpt - fJetPtMin)/fBinWidthJetPt);

  return bin;
}

//________________________________________________________________________
Double_t AliEmcalJetByJetCorrection::GetEfficiency(const Double_t pt) const {
  Double_t eff = 1.;
  if(fEfficiencyFixed<1.) return fEfficiencyFixed;
  else if(fhEfficiency) {
    Int_t bin = fhEfficiency->GetXaxis()->FindBin(pt);
    eff = fhEfficiency->GetBinContent(bin);
  }
  return eff;
}

//________________________________________________________________________
Double_t AliEmcalJetByJetCorrection::GetMeanPtConstituents(const AliEmcalJet *jet, TClonesArray *fTracks) const {

  if(!jet || !fTracks) return -1.;
  if(jet->GetNumberOfTracks()<1) return -1;

  AliVParticle *vp;
  Double_t sumPtCh = 0.;
  for(Int_t icc=0; icc<jet->GetNumberOfTracks(); icc++) {
    vp = static_cast<AliVParticle*>(jet->TrackAt(icc, fTracks));
    if(!vp) continue;
    sumPtCh+=vp->Pt();
  }
  Double_t meanpt = sumPtCh/(double)(jet->GetNumberOfTracks());
  return meanpt;
}
