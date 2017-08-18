//
// Utility class to do jet-by-jet correction
// Based on templates of ptJet vs ptTrack vs r
//
// Author: M.Verweij

#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliVParticle.h"
#include <TF1.h>
#include <TList.h>
#include <THnSparse.h>

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
  fhEfficiency(0),
  fhSmoothEfficiency(0),
  fCorrectpTtrack(0),
  fNpPoisson(0),
  fExternalNmissed(0),
  fRndm(0),
  fNMissedTracks(-1),
  fpAppliedEfficiency(0),
  fhNmissing(0),
  fListOfOutput(0)
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
  fhEfficiency(0),
  fhSmoothEfficiency(0),
  fCorrectpTtrack(0),
  fNpPoisson(0),
  fExternalNmissed(0),
  fRndm(0),
  fNMissedTracks(-1),
  fpAppliedEfficiency(0),
  fhNmissing(0),
  fListOfOutput(0)
{
  // Default constructor.
  fCollTemplates.SetOwner(kTRUE);

  const Int_t nBinPt = 118; //0-2: 20 bins; 2-100: 98 bins
  Double_t binLimitsPt[nBinPt+1];
  for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
    if(iPt<20){
      binLimitsPt[iPt] = 0. + (Double_t)iPt*0.15;
    } else {// 1.0
      binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 1.0;
    }
  }
  const Int_t nBinsPtJ  = 200;
  const Double_t minPtJ = -50.;
  const Double_t maxPtJ = 150.;

  fpAppliedEfficiency = new TProfile("fpAppliedEfficiency","fpAppliedEfficiency",nBinPt,binLimitsPt);
  
  const Int_t nvars = 5;
  Int_t nbins[nvars]  = {nBinsPtJ, 21 , 21 , 21 , 21};
  Double_t minbin[nvars] = {minPtJ  , 0. , 0. , 0. , 0.};
  Double_t maxbin[nvars] = {maxPtJ  , 20., 20., 20., 20.};
  TString title = "fhNmissing", nameh = title;
  TString axtitles[nvars] = {"#it{p}_{T,jet}", "N constituents added", "N_{constituents} #times (1/eff - 1)", "N_{truth}"};
  for(Int_t i = 0; i<nvars; i++){
     title+=axtitles[i];
  }
  fhNmissing = new THnSparseF(nameh.Data(), title.Data(), nvars, nbins, minbin, maxbin);
  
  fListOfOutput = new TList();
  fListOfOutput->SetName("JetByJetCorrectionOutput");
  fListOfOutput->SetOwner();
  fListOfOutput->Add(fpAppliedEfficiency);
  fListOfOutput->Add(fhNmissing);
 

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
  fhEfficiency(other.fhEfficiency),
  fhSmoothEfficiency(other.fhSmoothEfficiency),
  fCorrectpTtrack(other.fCorrectpTtrack),
  fpAppliedEfficiency(other.fpAppliedEfficiency),
  fRndm(other.fRndm),
  fhNmissing(other.fhNmissing)
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
  fhSmoothEfficiency = other.fhSmoothEfficiency;
  fCorrectpTtrack    = other.fCorrectpTtrack;
  fpAppliedEfficiency= other.fpAppliedEfficiency;
  fhNmissing          = other.fhNmissing;
  fRndm              = other.fRndm;
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
  fpAppliedEfficiency->Fill(meanPt,eff);

  Double_t fillarray[5]; //"#it{p}_{T,jet}", "N constituents added", "N_{constituents} #times (1/eff - 1)", "N_{truth}"
  fillarray[0] = jet->Pt();
  
  //np is the estimation of missed tracks
  Int_t np = TMath::FloorNint((double)jet->GetNumberOfTracks() * (1./eff -1.));
  fillarray[2] = np;
  
  if(fExternalNmissed) {
     np = fNMissedTracks; //if the number of missed tracks was calculated from external sources
     fillarray[3] = fNMissedTracks;
  }
  //npc is the number of added tracks
  Int_t npc=np; //take the particle missed as particle added
  if(fNpPoisson){
     npc=fRndm->Poisson(np); // smear the particle missed with a poissonian to get the number of added
  }
  fillarray[1] = npc;
  fhNmissing->Fill(fillarray);

  TLorentzVector corrVec; corrVec.SetPtEtaPhiM(jet->Pt(),jet->Eta(),jet->Phi(),jet->M());

  Double_t mass = 0.13957; //pion mass
  fArrayTrackCorr->Clear();
  for(Int_t i = 0; i<npc; i++) {
    Double_t r;
    Double_t pt;
    hTemplate->GetRandom2(r,pt);
    Double_t t = TMath::TwoPi()*gRandom->Uniform(1.);
    Double_t deta = r*TMath::Cos(t);
    Double_t dphi = r*TMath::Sin(t);
    TLorentzVector curVec; 
    curVec.SetPtEtaPhiM(pt,deta+jet->Eta(),dphi+jet->Phi(),mass);
    new ((*fArrayTrackCorr)[i]) TLorentzVector(curVec);
    corrVec+=curVec;
  }
  AliEmcalJet *jetCorr = new AliEmcalJet(corrVec.Pt(),corrVec.Eta(),corrVec.Phi(),corrVec.M());
  Printf("ERROR: Need to implemet a fix here to set the jet acceptance correctly -> not done in the constructor AliEmcalJet(pt,eta,phi,m) just used. See JIRA https://alice.its.cern.ch/jira/browse/ALPHY-65. Use something like jet->SetJetAcceptanceType(fJetTask->FindJetAcceptanceType(eta,phi,r))");
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
   if(!fhEfficiency && fEfficiencyFixed==1){
      Printf("%s fhEfficiency not known",GetName());
      fInitialized = kFALSE;
      return;
      
   }
   fRndm = new TRandom3(1234);
   if(fhEfficiency){
      fhSmoothEfficiency = (TH1D*) fh3JetPtDRTrackPt->ProjectionZ("fhSmoothEfficiency"); //copy the binning of the pTtrack axis
      Int_t nBinsEf = fhEfficiency->GetXaxis()->GetNbins();
      Int_t nBinsEfSmth = fhSmoothEfficiency->GetXaxis()->GetNbins();
      Double_t smallBinW = 0.9;
      Double_t pTFitRange[2]={6,100}; //reset in the loop, silent warning
      for(Int_t ibeff=0;ibeff<nBinsEf;ibeff++){
      	 Double_t bw = fhEfficiency->GetBinWidth(ibeff+1);
      	 if(bw>smallBinW){
      	    pTFitRange[0] = fhEfficiency->GetBinLowEdge(ibeff);     	 
      	    break;
      	 }
      }
      pTFitRange[1] = fhEfficiency->GetBinLowEdge(nBinsEf+1);
      
      TF1* fitfuncEff = new TF1("fitfuncEff", "[0]+x*[1]", pTFitRange[0], pTFitRange[1]);
      //fit function in the high pT region to smooth out fluctuations
      fhEfficiency->Fit("fitfuncEff", "R0");
      
      for(Int_t i=0;i<nBinsEfSmth;i++){
      	 Double_t binCentreT=fhSmoothEfficiency->GetBinCenter(i+1);
      	 Int_t bin = fhEfficiency->FindBin(binCentreT);
      	 //Double_t binEff = fhEfficiency->GetBinContent(bin);
      	 
      	 //fill histogram fhSmoothEfficiency by interpolation or function
      	 if(fhEfficiency->GetBinWidth(bin) > smallBinW){
      	    fhSmoothEfficiency->SetBinContent(i+1,fitfuncEff->Eval(binCentreT));      
      	 } else {
      	    Double_t effInterp = fhEfficiency->Interpolate(binCentreT);
      	    fhSmoothEfficiency ->SetBinContent(i+1,effInterp);
      	    fhSmoothEfficiency ->SetBinError(i+1,0.);
      	    
      	 }
      }
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
      
      fh3JetPtDRTrackPt->GetXaxis()->SetRange(binMin,binMax);
      Int_t nBinspTtr = fh3JetPtDRTrackPt->GetZaxis()->GetNbins();
      Int_t nBinsr = fh3JetPtDRTrackPt->GetYaxis()->GetNbins();
      TH2D *h2 = dynamic_cast<TH2D*>(fh3JetPtDRTrackPt->Project3D("zy"));
      if(h2){
      	 h2->SetName(Form("hPtR_%.0f_%.0f",fh3JetPtDRTrackPt->GetXaxis()->GetBinLowEdge(binMin),fh3JetPtDRTrackPt->GetXaxis()->GetBinUpEdge(binMax)));
      	 if(fCorrectpTtrack) {
      	    //apply efficiency correction to pTtrack
      	    for(Int_t ipTtr=0;ipTtr<nBinspTtr;ipTtr++){
      	       for(Int_t ir=0;ir<nBinsr;ir++){
      	       	  Double_t uncorr = h2->GetBinContent(ir+1,ipTtr+1);
      	       	  Double_t corr = uncorr/GetEfficiency(h2->GetYaxis()->GetBinCenter(ipTtr+1));
      	       	  h2->SetBinContent(ir+1, ipTtr+1, corr);
      	       }
      	       
      	    }
      	    
      	 }
      }
      fCollTemplates.Add(h2);
      fListOfOutput->Add(h2);
      counter++;
   }
   //  Int_t nt = TMath::FloorNint((fJetPtMax-fJetPtMin)/fBinWidthJetPt);
   //  Printf("nt: %d entries fCollTemplates: %d",nt,fCollTemplates.GetEntriesFast());
   
   fArrayTrackCorr = new TClonesArray("TLorentzVector", 20);
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
  else if(fhSmoothEfficiency) {
    Int_t bin = fhSmoothEfficiency->GetXaxis()->FindBin(pt);
    eff = fhSmoothEfficiency->GetBinContent(bin);
    //Printf("Efficiency (pT = %.2f, bin %d) = %.2f", pt, bin, eff);
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
