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

#include "AliTRDCalPIDLQ.h"
#include "AliTRDcalibDB.h"

ClassImp(AliTRDCalPIDLQ)

Float_t AliTRDCalPIDLQ::fTrackSegLength[kNLength] = { 3.7, 3.9, 4.2, 5.0 };

//_________________________________________________________________________
AliTRDCalPIDLQ::AliTRDCalPIDLQ()
  :AliTRDCalPID("pid", "LQ PID references for TRD")
{
  //
  //  The Default constructor
  //

  Init();

}

//_________________________________________________________________________
AliTRDCalPIDLQ::AliTRDCalPIDLQ(const Text_t *name, const Text_t *title)
  :AliTRDCalPID(name,title)
{
  //
  //  The main constructor
  //
  
  Init();

}

//_________________________________________________________________________
AliTRDCalPIDLQ::~AliTRDCalPIDLQ()
{
  //
  // Destructor
  //
  
}

//_________________________________________________________________________
Bool_t AliTRDCalPIDLQ::LoadReferences(Char_t *refFile)
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
			TH2D* hist = (TH2D*)histFile->Get(Form("h2dEdx%s%d", fPartSymb[iparticle], imom/*, ilength*/))->Clone();
			hist->Scale(1./hist->Integral());
			fModel->AddAt(hist, GetModelID(imom, iparticle, 0));

// 			if (iparticle != AliPID::kElectron && iparticle != AliPID::kPion) continue;
// 
// 			TH1F* ht = (TH1F*)histFile->Get(Form("h1MaxTB%s%02d", fPartSymb[iparticle], imom))->Clone();
// 			if(ht->GetNbinsX() != nTimeBins) AliWarning(Form("The number of time bins %d defined in h1MaxTB%s%02d differs from calibration value of %d. This may lead to erroneous results.", ht->GetNbinsX(), fPartSymb[iparticle], imom, nTimeBins));
// 			ht->Scale(1./ht->Integral());
// 			fHistTimeBin->AddAt(ht, ((iparticle==AliPID::kElectron)?0:1)*kNMom + imom);
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

//_________________________________________________________________________
TObject* AliTRDCalPIDLQ::GetModel(Int_t ip, Int_t iType, Int_t iplane) const          // iType not needed
{
  //
  // Returns one selected dEdx histogram
  //

  if (iType < 0 || iType >= AliPID::kSPECIES) return 0x0;
	if(ip<0 || ip>= kNMom ) return 0x0;

	AliInfo(Form("Retrive dEdx histogram for %s of %5.2f GeV/c", fPartName[iType], fTrackMomentum[ip]));
  
  return fModel->At(GetModelID(ip, iType, iplane));
}

//_________________________________________________________________________
Double_t AliTRDCalPIDLQ::GetProbability(Int_t spec, Float_t mom, Float_t *dedx, Float_t length, Int_t iplane) const
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
	if(!(hist = (TH2D*)fModel->At(GetModelID(imom-1, spec, iplane)))){
		AliInfo(Form("Looking for spec(%d) mom(%f) Ex(%f) Ey(%f) length(%f)", spec, mom, dedx[0], dedx[1], length));
		AliError(Form("EHistogram id %d not found in DB.", GetModelID(imom-1, spec, iplane)));
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


	if(!(hist = (TH2D*)fModel->At(GetModelID(imom, spec, iplane)))){
		AliInfo(Form("Looking for spec(%d) mom(%f) Ex(%f) Ey(%f) length(%f)", spec, mom, dedx[0], dedx[1], length));
		AliError(Form("EHistogram id %d not found in DB.", GetModelID(imom, spec, iplane)));
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
        if(mom < fTrackMomentum[0]) return LQ1;
        else if(mom > fTrackMomentum[kNMom-1]) return LQ2;
        else return LQ1 + (LQ2 - LQ1)*(mom - mom1)/(mom2 - mom1);

}

//_________________________________________________________________________
void AliTRDCalPIDLQ::Init()
{
  //
  // Initialization
  //

  fModel = new TObjArray(AliPID::kSPECIES  * kNMom);
  fModel -> SetOwner();

}

//_________________________________________________________________________
Int_t AliTRDCalPIDLQ::GetModelID(Int_t mom, Int_t spec, Int_t) const
{  
  // returns the ID of the LQ distribution (55 Histos, ID from 1 to 55)

	return spec * AliTRDCalPID::kNMom + mom;
}
