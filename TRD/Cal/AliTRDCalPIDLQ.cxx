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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Container for the distributions of dE/dx and the time bin of the          //
// max. cluster for electrons and pions                                      //
//                                                                           //
// Author:                                                                   //
//   Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de>                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TH1F.h>
#include <TFile.h>
#include <TROOT.h>

#include "AliLog.h"
#include "AliPID.h"

#include "AliTRDCalPIDLQ.h"

ClassImp(AliTRDCalPIDLQ)

Char_t* AliTRDCalPIDLQ::fpartName[AliPID::kSPECIES] = {"electron", "muon", "pion", "kaon", "proton"};
    
//_________________________________________________________________________
AliTRDCalPIDLQ::AliTRDCalPIDLQ()
  :TNamed()
  ,fNMom(0)
  ,fTrackMomentum(0)
  ,fMeanChargeRatio(0)
  ,fNbins(0)
  ,fBinSize(0)
  ,fHistdEdx(0)
  ,fHistTimeBin(0)
{
  //
  //  The Default constructor
  //
  
}

//_________________________________________________________________________
AliTRDCalPIDLQ::AliTRDCalPIDLQ(const Text_t *name, const Text_t *title) 
  :TNamed(name,title)
  ,fNMom(0)
  ,fTrackMomentum(0)
  ,fMeanChargeRatio(0)
  ,fNbins(0)
  ,fBinSize(0)
  ,fHistdEdx(0)
  ,fHistTimeBin(0)
{
  //
  //  The main constructor
  //
  
  Init();

}

//_____________________________________________________________________________
AliTRDCalPIDLQ::AliTRDCalPIDLQ(const AliTRDCalPIDLQ &c) 
  :TNamed(c)
  ,fNMom(c.fNMom)
  ,fTrackMomentum(0)
  ,fMeanChargeRatio(c.fMeanChargeRatio)
  ,fNbins(c.fNbins)
  ,fBinSize(c.fBinSize)
  ,fHistdEdx(0)
  ,fHistTimeBin(0)
{
  //
  // Copy constructor
  //

  AliTRDCalPIDLQ& target = (AliTRDCalPIDLQ &) c;
  
  target.fTrackMomentum = new Double_t[fNMom];
  for (Int_t i=0; i<fNMom; ++i) {
    target.fTrackMomentum[i] = fTrackMomentum[i];
  }
  if (fHistdEdx) {
    target.fHistdEdx = (TObjArray*) fHistdEdx->Clone();
  }  

  if (fHistTimeBin) {
    target.fHistTimeBin = (TObjArray*) fHistTimeBin->Clone();
  }

}

//_________________________________________________________________________
AliTRDCalPIDLQ::~AliTRDCalPIDLQ()
{
  //
  // Destructor
  //
  
  CleanUp();

}

//_________________________________________________________________________
void AliTRDCalPIDLQ::CleanUp()
{
  //
  // Delets all newly created objects
  //

  if (fHistdEdx) {
    delete fHistdEdx;
    fHistdEdx = 0;
  }
  
  if (fHistTimeBin) {
    delete fHistTimeBin;
    fHistTimeBin = 0;
  }

  if (fTrackMomentum) {
    delete[] fTrackMomentum;
    fTrackMomentum = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalPIDLQ &AliTRDCalPIDLQ::operator=(const AliTRDCalPIDLQ &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDCalPIDLQ &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDCalPIDLQ::Copy(TObject &c) const
{
  //
  // Copy function
  //

  AliTRDCalPIDLQ& target = (AliTRDCalPIDLQ &) c;
  
  target.CleanUp();
  
  target.fNMom            = fNMom;
  target.fNbins           = fNbins;
  target.fBinSize         = fBinSize;
  target.fMeanChargeRatio = fMeanChargeRatio;
  
  target.fTrackMomentum = new Double_t[fNMom];
  for (Int_t i=0; i<fNMom; ++i) {
    target.fTrackMomentum[i] = fTrackMomentum[i];
  }

  if (fHistdEdx) {
    target.fHistdEdx = (TObjArray*) fHistdEdx->Clone();
  }
  if (fHistTimeBin) {
    target.fHistTimeBin = (TObjArray*) fHistTimeBin->Clone();
  }

  TObject::Copy(c);

}

//_________________________________________________________________________
void AliTRDCalPIDLQ::Init()
{
  //
  // Initialization
  //

  const Int_t kNMom = 11;

  fNMom = kNMom;
  fTrackMomentum = new Double_t[fNMom];
  Double_t trackMomentum[kNMom] = {  0.6,  0.8,  1.0,  1.5,  2.0
                                  ,  3.0,  4.0,  5.0,  6.0,  8.0
                                  , 10.0 };
  for (Int_t imom = 0; imom < kNMom; imom++) {
    fTrackMomentum[imom] = trackMomentum[imom];
  }  

  fHistdEdx    = new TObjArray(AliPID::kSPECIES * fNMom);
  fHistdEdx->SetOwner();
  fHistTimeBin = new TObjArray(AliPID::kSPECIES * fNMom);
  fHistTimeBin->SetOwner();  
  
  // ADC Gain normalization
  fMeanChargeRatio = 1.0;
  
  // Number of bins and bin size
  fNbins   = 0;
  fBinSize = 0.0;

}

//_________________________________________________________________________
Bool_t AliTRDCalPIDLQ::ReadData(Char_t *responseFile)
{
  //
  // Read the TRD dEdx histograms.
  //

  // Read histogram Root file  
  TFile *histFile = new TFile(responseFile, "READ");
  if (!histFile || !histFile->IsOpen()) {
    AliError(Form("Opening TRD histgram file %s failed", responseFile));    
    return kFALSE;
  }
  gROOT->cd();

  // Read histograms
  Char_t text[200];
  for (Int_t particle = 0; particle < AliPID::kSPECIES; ++particle)
  {
    Char_t* particleKey = "";
    switch (particle)
    {
      case AliPID::kElectron: particleKey = "EL"; break;
      case AliPID::kPion: particleKey = "PI"; break;
      case AliPID::kMuon: particleKey = "MU"; break;
      case AliPID::kKaon: particleKey = "KA"; break;
      case AliPID::kProton: particleKey = "PR"; break;
    }
    
    for (Int_t imom = 0; imom < fNMom; imom++) 
    {
      sprintf(text, "h1dEdx%s%01d", particleKey, imom+1);
      TH1F* hist = (TH1F*)histFile->Get(text)->Clone();
      hist->Scale(1.0/hist->Integral());
      fHistdEdx->AddAt(hist, GetHistID(particle, imom));
  
      if (particle == AliPID::kElectron || particle == AliPID::kPion)
      {
        sprintf(text,"h1MaxTimBin%s%01d", particleKey, imom+1);
        TH1F* hist = (TH1F*)histFile->Get(text)->Clone();
        hist->Scale(1.0/hist->Integral());
        fHistTimeBin->AddAt(hist, GetHistID(particle,imom));
      }
    }
  }
  
  histFile->Close();
  delete histFile;
  
  // Number of bins and bin size
  TH1F* hist = (TH1F*) fHistdEdx->At(GetHistID(AliPID::kPion, 1));
  fNbins   = hist->GetNbinsX();
  fBinSize = hist->GetBinWidth(1);
  
  return kTRUE;

}

//_________________________________________________________________________
Double_t  AliTRDCalPIDLQ::GetMean(Int_t k, Int_t ip) const
{
  //
  // Gets mean of de/dx dist. of e
  //

  AliInfo(Form("Mean for particle = %s and momentum = %.2f is:\n"
              ,fpartName[k]
              ,fTrackMomentum[ip]));
  if (k < 0 || k > AliPID::kSPECIES) {
    return 0;
  }

  return ((TH1F*) fHistdEdx->At(GetHistID(k,ip)))->GetMean();

}

//_________________________________________________________________________
Double_t  AliTRDCalPIDLQ::GetNormalization(Int_t k, Int_t ip) const
{
  //
  // Gets Normalization of de/dx dist. of e
  //

  AliInfo(Form("Normalization for particle = %s and momentum = %.2f is:\n"
              ,fpartName[k]
              ,fTrackMomentum[ip]));
  if (k < 0 || k > AliPID::kSPECIES) {
    return 0;
  }
  
  return ((TH1F*) fHistdEdx->At(GetHistID(k,ip)))->Integral();

}

//_________________________________________________________________________
TH1F* AliTRDCalPIDLQ::GetHistogram(Int_t k, Int_t ip) const
{
  //
  // Returns one selected dEdx histogram
  //

  AliInfo(Form("Histogram for particle = %s and momentum = %.2f is:\n"
              ,fpartName[k]
              ,fTrackMomentum[ip]));
  if (k < 0 || k > AliPID::kSPECIES) {
    return 0;
  }
  
  return (TH1F*) fHistdEdx->At(GetHistID(k,ip));

}

//_________________________________________________________________________
TH1F* AliTRDCalPIDLQ::GetHistogramT(Int_t k, Int_t ip) const
{
  //
  // Returns one selected time bin max histogram
  //

  AliInfo(Form("Histogram for particle = %s and momentum = %.2f is:\n"
              ,fpartName[k]
              ,fTrackMomentum[ip]));
  if (k < 0 || k > AliPID::kSPECIES)
    return 0;
  
  return (TH1F*) fHistTimeBin->At(GetHistID(k,ip));

}

//_________________________________________________________________________
Double_t AliTRDCalPIDLQ::GetProbability(Int_t k, Double_t mom, Double_t dedx1) const
{
  //
  // Gets the Probability of having dedx at a given momentum (mom)
  // and particle type k (0 for e) and (2 for pi)
  // from the precalculated de/dx distributions 
  //
  
  Double_t dedx   = dedx1/fMeanChargeRatio;
  Int_t    iEnBin = ((Int_t) (dedx/fBinSize+1));
  if(iEnBin > fNbins) iEnBin = fNbins;

  if (k < 0 || k > AliPID::kSPECIES) {
    return 1;
  }
  
  TH1F* hist1 = 0;
  TH1F* hist2 = 0;
  Double_t mom1 = 0;
  Double_t mom2 = 0;
  
  // Lower limit
  if (mom<=fTrackMomentum[0])  {
    hist1 = (TH1F*) fHistdEdx->At(GetHistID(k,1));
    hist2 = (TH1F*) fHistdEdx->At(GetHistID(k,0));
    mom1 = fTrackMomentum[1];
    mom2 = fTrackMomentum[0];
  }
    
  // Upper Limit
  if(mom>=fTrackMomentum[fNMom-1]) {
    hist2 = (TH1F*) fHistdEdx->At(GetHistID(k,fNMom-1));
    hist1 = (TH1F*) fHistdEdx->At(GetHistID(k,fNMom-2));
    mom2 = fTrackMomentum[fNMom-1];
    mom1 = fTrackMomentum[fNMom-2];
  }
    
  // In the range
  for (Int_t ip=1; ip<fNMom; ip++) {
    if ((fTrackMomentum[ip-1]<= mom) && (mom<fTrackMomentum[ip])) {
      hist1 = (TH1F*) fHistdEdx->At(GetHistID(k,ip));
      hist2 = (TH1F*) fHistdEdx->At(GetHistID(k,ip-1));
      mom1 = fTrackMomentum[ip];
      mom2 = fTrackMomentum[ip-1];
    }
  }
  
  Double_t slop = (hist1->GetBinContent(iEnBin) - hist2->GetBinContent(iEnBin)) 
                / (mom1 - mom2);
  return hist2->GetBinContent(iEnBin) + slop * (mom - mom2);

}

//_________________________________________________________________________
Double_t AliTRDCalPIDLQ::GetProbabilityT(Int_t k, Double_t mom, Int_t timbin) const
{
  //
  // Gets the Probability of having timbin at a given momentum (mom)
  // and particle type k (0 for e) and (2 for pi)
  // from the precalculated timbin distributions 
  //
  
  if (timbin<=0) {
    return 0.0;
  }

  Int_t iTBin = timbin+1;
  
  // Everything which is not an electron counts as a pion for time bin max
  if (k != AliPID::kElectron) {
    k = AliPID::kPion;
  }

  if (mom<=fTrackMomentum[0]) {
    return ((TH1F*) fHistTimeBin->At(GetHistID(k,0)))->GetBinContent(iTBin);
  }
  if (mom>=fTrackMomentum[fNMom-1]) { 
    return ((TH1F*) fHistTimeBin->At(GetHistID(k,fNMom-1)))->GetBinContent(iTBin);
  }

  for (Int_t ip=1; ip<fNMom; ip++) {
    if ((fTrackMomentum[ip-1]<= mom) && (mom<fTrackMomentum[ip])) {
      Double_t slop = (((TH1F*) fHistTimeBin->At(GetHistID(k,ip)))->GetBinContent(iTBin) 
                     - ((TH1F*) fHistTimeBin->At(GetHistID(k,ip-1)))->GetBinContent(iTBin)) 
                    / (fTrackMomentum[ip] - fTrackMomentum[ip-1]);
      // Linear interpolation
      return ((TH1F*) fHistTimeBin->At(GetHistID(k,ip-1)))->GetBinContent(iTBin) 
              + slop * (mom - fTrackMomentum[ip-1]);
    }
  }
  
  return -1;

}
