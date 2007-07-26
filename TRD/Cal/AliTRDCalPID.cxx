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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Container for the distributions of dE/dx and the time bin of the     //
// max. cluster for electrons and pions                                 //
//                                                                      //
// Authors:                                                             //
//   Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de>               //
//   Alex Bercuci <a.bercuci@gsi.de>                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TROOT.h>

#include "AliLog.h"
#include "AliPID.h"
#include "AliESD.h"
#include "AliESDtrack.h"

#include "AliTRDCalPID.h"
#include "AliTRDcalibDB.h"

ClassImp(AliTRDCalPID)

Char_t* AliTRDCalPID::fpartName[AliPID::kSPECIES] = {"electron", "muon", "pion", "kaon", "proton"};

Char_t* AliTRDCalPID::fpartSymb[AliPID::kSPECIES] = {"EL", "MU", "PI", "KA", "PR"};

Float_t AliTRDCalPID::fTrackMomentum[kNMom] = {0.6,  0.8,  1.0,  1.5,  2.0,  3.0,  4.0,  5.0,  6.0,  8.0,  10.0};
  
Float_t AliTRDCalPID::fTrackSegLength[kNLength] = {3.7, 3.9, 4.2, 5.0};

    
//_________________________________________________________________________
AliTRDCalPID::AliTRDCalPID()
  :TNamed("pid", "PID for TRD")
  ,fMeanChargeRatio(0)
  ,fHistdEdx(0x0)
  ,fHistTimeBin(0x0)
{
  //
  //  The Default constructor
  //

  Init();

}

//_________________________________________________________________________
AliTRDCalPID::AliTRDCalPID(const Text_t *name, const Text_t *title) 
  :TNamed(name,title)
  ,fMeanChargeRatio(0)
  ,fHistdEdx(0x0)
  ,fHistTimeBin(0x0)
{
  //
  //  The main constructor
  //
  
  Init();

}

//_____________________________________________________________________________
AliTRDCalPID::AliTRDCalPID(const AliTRDCalPID &c) 
  :TNamed(c)
  ,fMeanChargeRatio(c.fMeanChargeRatio)
  ,fHistdEdx(0x0)
  ,fHistTimeBin(0x0)
{
  //
  // Copy constructor
  //

  if (this != &c) ((AliTRDCalPID &) c).Copy(*this);
  
}

//_________________________________________________________________________
AliTRDCalPID::~AliTRDCalPID()
{
  //
  // Destructor
  //
  
  CleanUp();

}

//_________________________________________________________________________
void AliTRDCalPID::CleanUp()
{
  //
  // Delets all newly created objects
  //

  if (fHistdEdx) {
    delete fHistdEdx;
    fHistdEdx = 0x0;
  }
  
  if (fHistTimeBin) {
    delete fHistTimeBin;
    fHistTimeBin = 0x0;
  }
}

//_____________________________________________________________________________
AliTRDCalPID &AliTRDCalPID::operator=(const AliTRDCalPID &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDCalPID &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDCalPID::Copy(TObject &c) const
{
  //
  // Copy function
  //

  AliTRDCalPID& target = (AliTRDCalPID &) c;
  
  target.CleanUp();
  
  target.fMeanChargeRatio = fMeanChargeRatio;

  if (fHistdEdx) {
    target.fHistdEdx = (TObjArray*) fHistdEdx->Clone();
  }
  if (fHistTimeBin) {
    target.fHistTimeBin = (TObjArray*) fHistTimeBin->Clone();
  }

  TObject::Copy(c);

}

//_________________________________________________________________________
void AliTRDCalPID::Init()
{
  //
  // Initialization
  //

  fHistdEdx    = new TObjArray(AliPID::kSPECIES * kNMom/* * kNLength*/);
  fHistdEdx->SetOwner();
  fHistTimeBin = new TObjArray(2 * kNMom);
  fHistTimeBin->SetOwner();  

	// Initialization of estimator at object instantiation because late
	// initialization in function GetProbability() is not working due to
	// constantness of this function. 
	// fEstimator = new AliTRDCalPIDRefMaker();
	
  // ADC Gain normalization
  fMeanChargeRatio = 1.0;
}

//_________________________________________________________________________
Bool_t AliTRDCalPID::LoadLQReferences(Char_t *refFile)
{
  //
  // Read the TRD dEdx histograms.
  //

	Int_t nTimeBins = 22;
	// Get number of time bins from CDB
	AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
	if(!calibration){
		AliWarning(Form("No AliTRDcalibDB available. Using %d time bins.", nTimeBins));
	}else{
		if(calibration->GetRun() > -1) nTimeBins = calibration->GetNumberOfTimeBins();
		else AliWarning(Form("Run number not set. Using %d time bins.", nTimeBins));
	}

	
  // Read histogram Root file  
  TFile *histFile = TFile::Open(refFile, "READ");
  if (!histFile || !histFile->IsOpen()) {
    AliError(Form("Opening TRD histgram file %s failed", refFile));
    return kFALSE;
  }
  gROOT->cd();

  // Read histograms
  for (Int_t iparticle = 0; iparticle < AliPID::kSPECIES; iparticle++){
    for (Int_t imom = 0; imom < kNMom; imom++){
			TH2D* hist = (TH2D*)histFile->Get(Form("h2dEdx%s%02d", fpartSymb[iparticle], imom/*, ilength*/))->Clone();
			hist->Scale(1./hist->Integral());
			fHistdEdx->AddAt(hist, GetHistID(iparticle, imom));

			if (iparticle != AliPID::kElectron && iparticle != AliPID::kPion) continue;

			TH1F* ht = (TH1F*)histFile->Get(Form("h1MaxTB%s%02d", fpartSymb[iparticle], imom))->Clone();
			if(ht->GetNbinsX() != nTimeBins) AliWarning(Form("The number of time bins %d defined in h1MaxTB%s%02d differs from calibration value of %d. This may lead to erroneous results.", ht->GetNbinsX(), fpartSymb[iparticle], imom, nTimeBins));
			ht->Scale(1./ht->Integral());
			fHistTimeBin->AddAt(ht, ((iparticle==AliPID::kElectron)?0:1)*kNMom + imom);
		}
  }
  
  histFile->Close();
  delete histFile;
  
  // Number of bins and bin size
  //TH1F* hist = (TH1F*) fHistdEdx->At(GetHistID(AliPID::kPion, 1));
  //fNbins   = hist->GetNbinsX();
  //fBinSize = hist->GetBinWidth(1);
  
  return kTRUE;

}

// //_________________________________________________________________________
// Double_t  AliTRDCalPID::GetMean(Int_t k, Int_t ip) const
// {
//   //
//   // Gets mean of de/dx dist. of e
//   //
// 
//   AliInfo(Form("Mean for particle = %s and momentum = %.2f is:\n"
//               ,fpartName[k]
//               ,fTrackMomentum[ip]));
//   if (k < 0 || k > AliPID::kSPECIES) {
//     return 0;
//   }
// 
//   return ((TH1F*) fHistdEdx->At(GetHistID(k,ip)))->GetMean();
// 
// }
// 
// //_________________________________________________________________________
// Double_t  AliTRDCalPID::GetNormalization(Int_t k, Int_t ip) const
// {
//   //
//   // Gets Normalization of de/dx dist. of e
//   //
// 
//   AliInfo(Form("Normalization for particle = %s and momentum = %.2f is:\n"
//               ,fpartName[k]
//               ,fTrackMomentum[ip]));
//   if (k < 0 || k > AliPID::kSPECIES) {
//     return 0;
//   }
//   
//   return ((TH1F*) fHistdEdx->At(GetHistID(k,ip)))->Integral();
// 
// }

//_________________________________________________________________________
TH1* AliTRDCalPID::GetHistogram(Int_t k, Int_t ip/*, Int_t il*/) const
{
  //
  // Returns one selected dEdx histogram
  //

  if (k < 0 || k >= AliPID::kSPECIES) return 0x0;
	if(ip<0 || ip>= kNMom ) return 0x0;

	AliInfo(Form("Retrive dEdx histogram for %s of %5.2f GeV/c", fpartName[k], fTrackMomentum[ip]));
  
  return (TH1*)fHistdEdx->At(GetHistID(k, ip));

}

//_________________________________________________________________________
TH1* AliTRDCalPID::GetHistogramT(Int_t k, Int_t ip) const
{
  //
  // Returns one selected time bin max histogram
  //

  if (k < 0 || k >= AliPID::kSPECIES) return 0x0;
	if(ip<0 || ip>= kNMom ) return 0x0;
	  
 	AliInfo(Form("Retrive MaxTB histogram for %s of %5.2f GeV/c", fpartName[k], fTrackMomentum[ip]));

	return (TH1*)fHistTimeBin->At(((k==AliPID::kElectron)?0:1)*kNMom+ip);
}



//_________________________________________________________________________
Double_t AliTRDCalPID::GetProbability(Int_t spec, Float_t mom, Float_t *dedx, Float_t length) const
{
  //
	// Core function of AliTRDCalPID class for calculating the
	// likelihood for species "spec" (see AliTRDtrack::kNspecie) of a
	// given momentum "mom" and a given dE/dx (slice "dedx") yield per
	// layer
  //

	if (spec < 0 || spec >= AliPID::kSPECIES) return 0.;
		
	//Double_t dedx   = dedx1/fMeanChargeRatio;
	
	// find the interval in momentum and track segment length which applies for this data
	Int_t ilength = 1;
  while(ilength<kNLength-1 && length>fTrackSegLength[ilength]){
		ilength++;
	}
	Int_t imom = 1;
  while(imom<kNMom-1 && mom>fTrackMomentum[imom]) imom++;
	
	Int_t nbinsx, nbinsy;
	TAxis *ax = 0x0, *ay = 0x0;
	Double_t LQ1, LQ2;
	Double_t mom1 = fTrackMomentum[imom-1], mom2 = fTrackMomentum[imom];
	TH2 *hist = 0x0;
	if(!(hist = (TH2D*)fHistdEdx->At(GetHistID(spec, imom-1/*, ilength*/)))){
		AliInfo(Form("Looking for spec(%d) mom(%f) Ex(%f) Ey(%f) length(%f)", spec, mom, dedx[0], dedx[1], length));
		AliError(Form("EHistogram id %d not found in DB.", GetHistID(spec, imom-1)));
		return 0.;
	}
	ax = hist->GetXaxis(); nbinsx = ax->GetNbins();
	ay = hist->GetYaxis(); nbinsy = ay->GetNbins();
	Float_t x = dedx[0]+dedx[1], y = dedx[2];
  Bool_t kX = (x < ax->GetBinUpEdge(nbinsx));
	Bool_t kY = (y < ay->GetBinUpEdge(nbinsy));
	if(kX)
		if(kY) LQ1 = hist->GetBinContent( hist->FindBin(x, y)); 
    //fEstimator->Estimate2D2(hist, x, y);
		else LQ1 = hist->GetBinContent(ax->FindBin(x), nbinsy);
	else
		if(kY) LQ1 = hist->GetBinContent(nbinsx, ay->FindBin(y));
	 	else LQ1 = hist->GetBinContent(nbinsx, nbinsy);


	if(!(hist = (TH2D*)fHistdEdx->At(GetHistID(spec, imom/*, ilength*/)))){
		AliInfo(Form("Looking for spec(%d) mom(%f) Ex(%f) Ey(%f) length(%f)", spec, mom, dedx[0], dedx[1], length));
		AliError(Form("EHistogram id %d not found in DB.", GetHistID(spec, imom)));
		return LQ1;
	}
	if(kX)
		if(kY) LQ2 = hist->GetBinContent( hist->FindBin(x, y)); 
    //fEstimator->Estimate2D2(hist, x, y);
		else LQ2 = hist->GetBinContent(ax->FindBin(x), nbinsy);
	else
		if(kY) LQ2 = hist->GetBinContent(nbinsx, ay->FindBin(y));
	 	else LQ2 = hist->GetBinContent(nbinsx, nbinsy);

	
	// return interpolation over momentum binning
  return LQ1 + (LQ2 - LQ1)*(mom - mom1)/(mom2 - mom1);

}

//_________________________________________________________________________
Double_t AliTRDCalPID::GetProbabilityT(Int_t spec, Double_t mom, Int_t timbin) const
{
  //
  // Gets the Probability of having timbin at a given momentum (mom)
  // and particle type (spec) (0 for e) and (2 for pi)
  // from the precalculated timbin distributions 
  //
  
	if (spec < 0 || spec >= AliPID::kSPECIES) return 0.;

  Int_t iTBin = timbin+1;
  
  // Everything which is not an electron counts as a pion for time bin max
  //if(spec != AliPID::kElectron) spec = AliPID::kPion;
  

  
	Int_t imom = 1;
  while(imom<kNMom-1 && mom>fTrackMomentum[imom]) imom++;

	Double_t mom1 = fTrackMomentum[imom-1], mom2 = fTrackMomentum[imom];
	TH1F *hist = 0x0;
	if(!(hist = (TH1F*) fHistTimeBin->At(((spec==AliPID::kElectron)?0:1)*kNMom+imom-1))){
		AliInfo(Form("Looking for spec(%d) mom(%f) timbin(%d)", spec, mom, timbin));
		AliError(Form("THistogram id %d not found in DB.", ((spec==AliPID::kElectron)?0:1)*kNMom+imom-1));
		return 0.;
	}
	Double_t LQ1 = hist->GetBinContent(iTBin);

	if(!(hist = (TH1F*) fHistTimeBin->At(((spec==AliPID::kElectron)?0:1)*kNMom+imom))){
		AliInfo(Form("Looking for spec(%d) mom(%f) timbin(%d)", spec, mom, timbin));
		AliError(Form("THistogram id %d not found in DB.", ((spec==AliPID::kElectron)?0:1)*kNMom+imom));
		return LQ1;
	}
	Double_t LQ2 = hist->GetBinContent(iTBin);

	// return interpolation over momentum binning
  return LQ1 + (LQ2 - LQ1)*(mom - mom1)/(mom2 - mom1);
}

