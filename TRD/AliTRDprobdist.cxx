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

//-----------------------------------------------------------------
// Class for dE/dx and Time Bin of Max. Cluster for Electrons and 
// pions in TRD. 
// It is instantiated in class AliTRDpidESD for particle identification
// in TRD
// Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de>
//-----------------------------------------------------------------

#include "AliTRDprobdist.h"

ClassImp(AliTRDprobdist)

//_________________________________________________________________________
AliTRDprobdist::AliTRDprobdist(Char_t *responseFile)
{
  //
  //  The main constructor
  //
  Char_t *partName[5] = {"electron", "muon", "pion", "kaon", "proton"};
  for(Int_t i=0; i<5; i++) {
    fpartName[i] = partName[i]; 
  }
  // ADC Gain normalization
  fADCNorm=1.0;
  fNMom=kNMom;
  // Set track momenta for which response functions are available
  Double_t trackMomentum[kNMom]= {0.6, 0.8, 1, 1.5, 2, 3, 4, 5, 6, 8, 10};
  for(Int_t imom=0; imom<fNMom; imom++) {
    fTrackMomentum[imom] = trackMomentum[imom];
  }
  // Read the histogram file
  ReadData(responseFile);
  //
}

Bool_t AliTRDprobdist::ReadData(Char_t *responseFile)
{
  //
  // Read the TRD dEdx histograms.
  //
  // Read histogram Root file  
  TFile *histFile = new TFile(responseFile);
  if (!histFile || !histFile->IsOpen()) {
    Error("AliTRDprobdist", "opening TRD histgram file %s failed", responseFile);
    return kFALSE;
  }

  // Read histograms
  Char_t text[200];
  for (Int_t imom = 0; imom < kNMom; imom++) {
    sprintf(text,"h1dEdxEL%01d",imom+1);
    fh1dEdxEL[imom] = (TH1F*)histFile->Get(text);
    fh1dEdxEL[imom]->Scale(1.0/fh1dEdxEL[imom]->Integral());

    sprintf(text,"h1dEdxPI%01d",imom+1);
    fh1dEdxPI[imom] = (TH1F*)histFile->Get(text);
    fh1dEdxPI[imom]->Scale(1.0/fh1dEdxPI[imom]->Integral());

    sprintf(text,"h1dEdxMU%01d",imom+1);
    fh1dEdxMU[imom] = (TH1F*)histFile->Get(text);
    fh1dEdxMU[imom]->Scale(1.0/fh1dEdxMU[imom]->Integral());

    sprintf(text,"h1dEdxKA%01d",imom+1);
    fh1dEdxKA[imom] = (TH1F*)histFile->Get(text);
    fh1dEdxKA[imom]->Scale(1.0/fh1dEdxKA[imom]->Integral());

    sprintf(text,"h1dEdxPR%01d",imom+1);
    fh1dEdxPR[imom] = (TH1F*)histFile->Get(text);
    fh1dEdxPR[imom]->Scale(1.0/fh1dEdxPR[imom]->Integral());

    sprintf(text,"h1MaxTimBinEL%01d",imom+1);
    fh1MaxTimBinEL[imom] = (TH1F*)histFile->Get(text);
    fh1MaxTimBinEL[imom]->Scale(1.0/fh1MaxTimBinEL[imom]->Integral());

    sprintf(text,"h1MaxTimBinPI%01d",imom+1);
    fh1MaxTimBinPI[imom] = (TH1F*)histFile->Get(text);
    fh1MaxTimBinPI[imom]->Scale(1.0/fh1MaxTimBinPI[imom]->Integral());
  }
  // Number of bins and bin size
  fNbins = fh1dEdxPI[1]->GetNbinsX();
  fBinSize = fh1dEdxPI[1]->GetBinWidth(1);
  return kTRUE;
}

AliTRDprobdist::AliTRDprobdist(const AliTRDprobdist& pd):TNamed()
{
  //
  // Copy constructor.
  //
  for(Int_t i=0; i<5; i++) {
    fpartName[i] = pd.fpartName[i]; 
  }
  // ADC Gain normalization
  fADCNorm=pd.fADCNorm;
  fNMom=pd.fNMom;
  // Set track momenta for which response functions are available
  for(Int_t imom=0; imom<fNMom; imom++) {
    fTrackMomentum[imom] = pd.fTrackMomentum[imom];
  }

  fNbins=pd.fNbins;
  fBinSize=pd.fBinSize;  

  for (Int_t imom = 0; imom < kNMom; imom++) {
    fh1dEdxEL[imom] = pd.fh1dEdxEL[imom];
    fh1dEdxPI[imom] = pd.fh1dEdxPI[imom];
    fh1dEdxMU[imom] = pd.fh1dEdxMU[imom];
    fh1dEdxKA[imom] = pd.fh1dEdxKA[imom];
    fh1dEdxPR[imom] = pd.fh1dEdxPR[imom];
    fh1MaxTimBinEL[imom] = pd.fh1MaxTimBinEL[imom];
    fh1MaxTimBinPI[imom] = pd.fh1MaxTimBinPI[imom];
  }
}


AliTRDprobdist::~AliTRDprobdist()
{
  // Destructor

}


//_________________________________________________________________________
Double_t  AliTRDprobdist::GetMean(Int_t k, Int_t ip) const
{
  //
  // Gets mean of de/dx dist. of e
  printf("Mean for particle = %s and momentum = %.2f is:\n", fpartName[k], fTrackMomentum[ip]);
  if(k==0) return fh1dEdxEL[ip]->GetMean();
  if(k==1) return fh1dEdxMU[ip]->GetMean();
  if(k==2) return fh1dEdxPI[ip]->GetMean();
  if(k==3) return fh1dEdxKA[ip]->GetMean();
  if(k==4) return fh1dEdxPR[ip]->GetMean();
  return fh1dEdxPR[ip]->GetMean();
}

//_________________________________________________________________________
Double_t  AliTRDprobdist::GetNormalization(Int_t k, Int_t ip) const
{
  //
  // Gets Normalization of de/dx dist. of e

  printf("Normalization for particle = %s and momentum = %.2f is:\n",fpartName[k], fTrackMomentum[ip]);
  if(k==0) return fh1dEdxEL[ip]->Integral();
  if(k==1) return fh1dEdxMU[ip]->Integral();
  if(k==2) return fh1dEdxPI[ip]->Integral();
  if(k==3) return fh1dEdxKA[ip]->Integral();
  if(k==4) return fh1dEdxPR[ip]->Integral();
  return fh1dEdxPR[ip]->Integral();
}

TH1F* AliTRDprobdist::GetHistogram(Int_t k, Int_t ip) const
{
  //
  //
  printf("Histogram for particle = %s and momentum = %.2f is:\n", fpartName[k], fTrackMomentum[ip]);
  if(k==0) return fh1dEdxEL[ip];
  if(k==1) return fh1dEdxMU[ip];
  if(k==2) return fh1dEdxPI[ip];
  if(k==3) return fh1dEdxKA[ip];
  if(k==4) return fh1dEdxPR[ip];
  return fh1dEdxPR[ip];
}

//_________________________________________________________________________
Double_t AliTRDprobdist::GetProbability(Int_t k, Double_t mom, Double_t dedx1) const
{
  //
  // Gets the Probability of having dedx at a given momentum (mom)
  // and particle type k (0 for e) and (2 for pi)
  // from the precalculated de/dx distributions 
  Double_t probability = 1.0;
  Double_t dedx = dedx1/fADCNorm;
  Int_t iEnBin= ((Int_t) (dedx/fBinSize+1));
  if(iEnBin > fNbins) iEnBin = fNbins;

  ////Electron//////////////////////////
  if(k==0){    // electron 
    Double_t slop;
    // Lower limit
    if(mom<=fTrackMomentum[0]) {
      slop=(fh1dEdxEL[1]->GetBinContent(iEnBin)-fh1dEdxEL[0]->GetBinContent(iEnBin))/(fTrackMomentum[1] - fTrackMomentum[0]);
      probability= fh1dEdxEL[0]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[0]);
      return probability;
    }
    // Upper Limit
    if(mom>=fTrackMomentum[fNMom-1]) {
      slop=(fh1dEdxEL[fNMom-1]->GetBinContent(iEnBin)-fh1dEdxEL[fNMom-2]->GetBinContent(iEnBin))/(fTrackMomentum[fNMom-1] - fTrackMomentum[fNMom-2]);
      probability= fh1dEdxEL[fNMom-1]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[fNMom-1]);
      return probability;
    }
    // In the range
    for(Int_t ip=1; ip<fNMom; ip++) 
      if((fTrackMomentum[ip-1]<= mom) && (mom<fTrackMomentum[ip])) {
	slop=(fh1dEdxEL[ip]->GetBinContent(iEnBin)-fh1dEdxEL[ip-1]->GetBinContent(iEnBin))/(fTrackMomentum[ip] - fTrackMomentum[ip-1]);
	// Linear Interpolation
	probability= fh1dEdxEL[ip-1]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[ip-1]);
	return probability;
      }
  }

  ////Pion//////////////////////////
  if(k==2){    // Pion
    Double_t slop;
    // Lower limit
    if(mom<=fTrackMomentum[0]) {
      slop=(fh1dEdxPI[1]->GetBinContent(iEnBin)-fh1dEdxPI[0]->GetBinContent(iEnBin))/(fTrackMomentum[1] - fTrackMomentum[0]);
      probability= fh1dEdxPI[0]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[0]);
      return probability;
    }
    // Upper Limit
    if(mom>=fTrackMomentum[fNMom-1]) {
      slop=(fh1dEdxPI[fNMom-1]->GetBinContent(iEnBin)-fh1dEdxPI[fNMom-2]->GetBinContent(iEnBin))/(fTrackMomentum[fNMom-1] - fTrackMomentum[fNMom-2]);
      probability= fh1dEdxPI[fNMom-1]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[fNMom-1]);
      return probability;
    }
    // In the range
    for(Int_t ip=1; ip<fNMom; ip++) 
      if((fTrackMomentum[ip-1]<= mom) && (mom<fTrackMomentum[ip])) {
	slop=(fh1dEdxPI[ip]->GetBinContent(iEnBin)-fh1dEdxPI[ip-1]->GetBinContent(iEnBin))/(fTrackMomentum[ip] - fTrackMomentum[ip-1]);
	// Linear Interpolation
	probability= fh1dEdxPI[ip-1]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[ip-1]);
	return probability;
      }
  }

  ////Muon//////////////////////////
  if(k==1){    // Muon
    Double_t slop;
    // Lower limit
    if(mom<=fTrackMomentum[0]) {
      slop=(fh1dEdxMU[1]->GetBinContent(iEnBin)-fh1dEdxMU[0]->GetBinContent(iEnBin))/(fTrackMomentum[1] - fTrackMomentum[0]);
      probability= fh1dEdxMU[0]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[0]);
      return probability;
    }
    // Upper Limit
    if(mom>=fTrackMomentum[fNMom-1]) {
      slop=(fh1dEdxMU[fNMom-1]->GetBinContent(iEnBin)-fh1dEdxMU[fNMom-2]->GetBinContent(iEnBin))/(fTrackMomentum[fNMom-1] - fTrackMomentum[fNMom-2]);
      probability= fh1dEdxMU[fNMom-1]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[fNMom-1]);
      return probability;
    }
    // In the range
    for(Int_t ip=1; ip<fNMom; ip++) 
      if((fTrackMomentum[ip-1]<= mom) && (mom<fTrackMomentum[ip])) {
	slop=(fh1dEdxMU[ip]->GetBinContent(iEnBin)-fh1dEdxMU[ip-1]->GetBinContent(iEnBin))/(fTrackMomentum[ip] - fTrackMomentum[ip-1]);
	// Linear Interpolation
	probability= fh1dEdxMU[ip-1]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[ip-1]);
	return probability;
      }
  }

  ////Kaon//////////////////////////
  if(k==3){    // Kaon
    Double_t slop;
    // Lower limit
    if(mom<=fTrackMomentum[0]) {
      slop=(fh1dEdxKA[1]->GetBinContent(iEnBin)-fh1dEdxKA[0]->GetBinContent(iEnBin))/(fTrackMomentum[1] - fTrackMomentum[0]);
      probability= fh1dEdxKA[0]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[0]);
      return probability;
    }
    // Upper Limit
    if(mom>=fTrackMomentum[fNMom-1]) {
      slop=(fh1dEdxKA[fNMom-1]->GetBinContent(iEnBin)-fh1dEdxKA[fNMom-2]->GetBinContent(iEnBin))/(fTrackMomentum[fNMom-1] - fTrackMomentum[fNMom-2]);
      probability= fh1dEdxKA[fNMom-1]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[fNMom-1]);
      return probability;
    }
    // In the range
    for(Int_t ip=1; ip<fNMom; ip++) 
      if((fTrackMomentum[ip-1]<= mom) && (mom<fTrackMomentum[ip])) {
	slop=(fh1dEdxKA[ip]->GetBinContent(iEnBin)-fh1dEdxKA[ip-1]->GetBinContent(iEnBin))/(fTrackMomentum[ip] - fTrackMomentum[ip-1]);
	// Linear Interpolation
	probability= fh1dEdxKA[ip-1]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[ip-1]);
	return probability;
      }
  }

  ////Proton//////////////////////////
  if(k==4){    // Proton
    Double_t slop;
    // Lower limit
    if(mom<=fTrackMomentum[0]) {
      slop=(fh1dEdxPR[1]->GetBinContent(iEnBin)-fh1dEdxPR[0]->GetBinContent(iEnBin))/(fTrackMomentum[1] - fTrackMomentum[0]);
      probability= fh1dEdxPR[0]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[0]);
      return probability;
    }
    // Upper Limit
    if(mom>=fTrackMomentum[fNMom-1]) {
      slop=(fh1dEdxPR[fNMom-1]->GetBinContent(iEnBin)-fh1dEdxPR[fNMom-2]->GetBinContent(iEnBin))/(fTrackMomentum[fNMom-1] - fTrackMomentum[fNMom-2]);
      probability= fh1dEdxPR[fNMom-1]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[fNMom-1]);
      return probability;
    }
    // In the range
    for(Int_t ip=1; ip<fNMom; ip++) 
      if((fTrackMomentum[ip-1]<= mom) && (mom<fTrackMomentum[ip])) {
	slop=(fh1dEdxPR[ip]->GetBinContent(iEnBin)-fh1dEdxPR[ip-1]->GetBinContent(iEnBin))/(fTrackMomentum[ip] - fTrackMomentum[ip-1]);
	// Linear Interpolation
	probability= fh1dEdxPR[ip-1]->GetBinContent(iEnBin) + slop*(mom-fTrackMomentum[ip-1]);
	return probability;
      }
  }

  return probability;
}


//_________________________________________________________________________
Double_t AliTRDprobdist::GetProbabilityT(Int_t k, Double_t mom, Int_t timbin) const
{
  //
  // Gets the Probability of having timbin at a given momentum (mom)
  // and particle type k (0 for e) and (2 for pi)
  // from the precalculated timbin distributions 
  Double_t probabilityT = 1.0;
  if(timbin<=0) return 0.;
  Int_t iTBin=timbin+1;

  if(k==0){    // electron
    if(mom<=fTrackMomentum[0]) probabilityT = fh1MaxTimBinEL[0]->GetBinContent(iTBin);
    if(mom>=fTrackMomentum[fNMom-1]) probabilityT = fh1MaxTimBinEL[fNMom-1]->GetBinContent(iTBin);
  }
  if(k==1||k==2||k==3||k==4){    // pion
    if(mom<=fTrackMomentum[0]) probabilityT = fh1MaxTimBinPI[0]->GetBinContent(iTBin);
    if(mom>=fTrackMomentum[fNMom-1]) probabilityT = fh1MaxTimBinPI[fNMom-1]->GetBinContent(iTBin);
  }

  if(k==0)    // electron
  for(Int_t ip=1; ip<fNMom; ip++)
  if((fTrackMomentum[ip-1]<= mom) && (mom<fTrackMomentum[ip])) {
    Double_t slop=(fh1MaxTimBinEL[ip]->GetBinContent(iTBin)-fh1MaxTimBinEL[ip-1]->GetBinContent(iTBin))/(fTrackMomentum[ip] - fTrackMomentum[ip-1]);
    // Linear Interpolation
    probabilityT= fh1MaxTimBinEL[ip-1]->GetBinContent(iTBin) + slop*(mom-fTrackMomentum[ip-1]);
    return probabilityT;
  }

  if(k==1||k==2||k==3||k==4)   // pion and other particles
  for(Int_t ip=1; ip<fNMom; ip++)
  if((fTrackMomentum[ip-1]<= mom) && (mom<fTrackMomentum[ip])) {
    Double_t slop=(fh1MaxTimBinPI[ip]->GetBinContent(iTBin)-fh1MaxTimBinPI[ip-1]->GetBinContent(iTBin))/(fTrackMomentum[ip] - fTrackMomentum[ip-1]);
    // Linear Interpolation
    probabilityT= fh1MaxTimBinPI[ip-1]->GetBinContent(iTBin) + slop*(mom-fTrackMomentum[ip-1]);
    return probabilityT;
  }
  return probabilityT;
}


