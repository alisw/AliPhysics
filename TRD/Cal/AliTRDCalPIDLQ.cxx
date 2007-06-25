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
#include <TTree.h>
#include <TROOT.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TPrincipal.h>

#include "AliLog.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliPID.h"
#include "AliESD.h"
#include "AliESDtrack.h"

#include "AliTRDCalPIDLQ.h"
#include "AliTRDCalPIDLQRef.h"
#include "AliTRDpidESD.h"
#include "AliTRDcalibDB.h"
#include "AliTRDtrack.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDCalPIDLQ)

Char_t* AliTRDCalPIDLQ::fpartName[AliPID::kSPECIES] = {"electron", "muon", "pion", "kaon", "proton"};
Char_t* AliTRDCalPIDLQ::fpartSymb[AliPID::kSPECIES] = {"EL", "MU", "PI", "KA", "PR"};
    
//_________________________________________________________________________
AliTRDCalPIDLQ::AliTRDCalPIDLQ()
  :TNamed("pid", "PID for TRD")
  ,fNMom(0)
  ,fNLength(0)
  ,fTrackMomentum(0x0)
  ,fTrackSegLength(0x0)
  ,fNTimeBins(0)
  ,fMeanChargeRatio(0)
  ,fNbins(0)
  ,fBinSize(0)
  ,fHistdEdx(0x0)
  ,fHistTimeBin(0x0)
  ,fEstimator(0x0)
{
  //
  //  The Default constructor
  //

  Init();

}

//_________________________________________________________________________
AliTRDCalPIDLQ::AliTRDCalPIDLQ(const Text_t *name, const Text_t *title) 
  :TNamed(name,title)
  ,fNMom(0)
  ,fNLength(0)
  ,fTrackMomentum(0x0)
  ,fTrackSegLength(0x0)
  ,fNTimeBins(0)
  ,fMeanChargeRatio(0)
  ,fNbins(0)
  ,fBinSize(0)
  ,fHistdEdx(0x0)
  ,fHistTimeBin(0x0)
  ,fEstimator(0x0)
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
  ,fNLength(c.fNLength)
  ,fTrackMomentum(0x0)
  ,fTrackSegLength(0x0)
  ,fNTimeBins(c.fNTimeBins)
  ,fMeanChargeRatio(c.fMeanChargeRatio)
  ,fNbins(c.fNbins)
  ,fBinSize(c.fBinSize)
  ,fHistdEdx(0x0)
  ,fHistTimeBin(0x0)
  ,fEstimator(0x0)
{
  //
  // Copy constructor
  //

  if (this != &c) ((AliTRDCalPIDLQ &) c).Copy(*this);
  
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
    fHistdEdx = 0x0;
  }
  
  if (fHistTimeBin) {
    delete fHistTimeBin;
    fHistTimeBin = 0x0;
  }

  if (fTrackMomentum) {
    delete[] fTrackMomentum;
    fTrackMomentum = 0x0;
  }

  if (fTrackSegLength) {
    delete[] fTrackSegLength;
    fTrackSegLength = 0x0;
  }

  if (fEstimator) {
    delete fEstimator;
    fEstimator = 0x0;
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
  
  target.fNbins           = fNbins;
  target.fBinSize         = fBinSize;
  target.fMeanChargeRatio = fMeanChargeRatio;
  target.fNTimeBins       = fNTimeBins;

  //target.fNMom            = fNMom;
  target.fTrackMomentum = new Double_t[fNMom];
  for (Int_t i=0; i<fNMom; ++i) {
    target.fTrackMomentum[i] = fTrackMomentum[i];
  }

  //target.fNLength         = fNLength;
  target.fTrackSegLength = new Double_t[fNLength];
  for (Int_t i=0; i<fNLength; ++i) {
    target.fTrackSegLength[i] = fTrackSegLength[i];
  }

  if (fHistdEdx) {
    target.fHistdEdx = (TObjArray*) fHistdEdx->Clone();
  }
  if (fHistTimeBin) {
    target.fHistTimeBin = (TObjArray*) fHistTimeBin->Clone();
  }

  target.fEstimator = new AliTRDCalPIDLQRef(*fEstimator);

  TObject::Copy(c);

}

//_________________________________________________________________________
void AliTRDCalPIDLQ::Init()
{
  //
  // Initialization
  //

  fNMom    = 11;
  fNLength =  4;

  fTrackMomentum = new Double_t[fNMom];
  fTrackMomentum[0] = 0.6;
  fTrackMomentum[1] = 0.8;
  fTrackMomentum[2] = 1.0;
  fTrackMomentum[3] = 1.5;
  fTrackMomentum[4] = 2.0;
  fTrackMomentum[5] = 3.0;
  fTrackMomentum[6] = 4.0;
  fTrackMomentum[7] = 5.0;
  fTrackMomentum[8] = 6.0;
  fTrackMomentum[9] = 8.0;
  fTrackMomentum[10] = 10.0;
  
  fTrackSegLength = new Double_t[fNLength];
  fTrackSegLength[0] = 3.7;
  fTrackSegLength[1] = 3.9;
  fTrackSegLength[2] = 4.2;
  fTrackSegLength[3] = 5.0;

  fHistdEdx    = new TObjArray(AliPID::kSPECIES * fNMom/* * fNLength*/);
  fHistdEdx->SetOwner();
  fHistTimeBin = new TObjArray(2 * fNMom);
  fHistTimeBin->SetOwner();  

	// Initialization of estimator at object instantiation because late
	// initialization in function GetProbability() is not working due to
	// constantness of this function. 
	fEstimator = new AliTRDCalPIDLQRef();

	// Number of Time bins
	AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
	if(!calibration){
	  AliWarning("No AliTRDcalibDB available. Using 22 time bins.");
	  fNTimeBins = 22;
	} else {
          if (calibration->GetRun() > -1) {
            fNTimeBins = calibration->GetNumberOfTimeBins();
	  }
          else {
  	    AliWarning("No run number set. Using 22 time bins.");
	    fNTimeBins = 22;
	  }
	}
	
  // ADC Gain normalization
  fMeanChargeRatio = 1.0;
  
  // Number of bins and bin size
  fNbins   = 0;
  fBinSize = 0.0;
}

//_________________________________________________________________________
Bool_t AliTRDCalPIDLQ::ReadReferences(Char_t *responseFile)
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
  for (Int_t iparticle = 0; iparticle < AliPID::kSPECIES; iparticle++){
    for (Int_t imom = 0; imom < fNMom; imom++){
			TH2D* hist = (TH2D*)histFile->Get(Form("h2dEdx%s%02d", fpartSymb[iparticle], imom/*, ilength*/))->Clone();
			hist->Scale(1./hist->Integral());
			fHistdEdx->AddAt(hist, GetHistID(iparticle, imom));

			if (iparticle != AliPID::kElectron && iparticle != AliPID::kPion) continue;

			TH1F* ht = (TH1F*)histFile->Get(Form("h1MaxTB%s%02d", fpartSymb[iparticle], imom))->Clone();
			if(ht->GetNbinsX() != fNTimeBins) AliWarning(Form("The number of time bins %d defined in h1MaxTB%s%02d differs from calibration value of %d. This may lead to erroneous results.", ht->GetNbinsX(), fpartSymb[iparticle], imom, fNTimeBins));
			ht->Scale(1./ht->Integral());
			fHistTimeBin->AddAt(ht, ((iparticle==AliPID::kElectron)?0:1)*fNMom + imom);
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
// Double_t  AliTRDCalPIDLQ::GetMean(Int_t k, Int_t ip) const
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
// Double_t  AliTRDCalPIDLQ::GetNormalization(Int_t k, Int_t ip) const
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
TH1* AliTRDCalPIDLQ::GetHistogram(Int_t k, Int_t ip/*, Int_t il*/) const
{
  //
  // Returns one selected dEdx histogram
  //

  if (k < 0 || k >= AliPID::kSPECIES) return 0x0;
	if(ip<0 || ip>= fNMom ) return 0x0;

	AliInfo(Form("Retrive dEdx histogram for %s of %5.2f GeV/c", fpartName[k], fTrackMomentum[ip]));
  
  return (TH1*)fHistdEdx->At(GetHistID(k, ip));

}

//_________________________________________________________________________
TH1* AliTRDCalPIDLQ::GetHistogramT(Int_t k, Int_t ip) const
{
  //
  // Returns one selected time bin max histogram
  //

  if (k < 0 || k >= AliPID::kSPECIES) return 0x0;
	if(ip<0 || ip>= fNMom ) return 0x0;
	  
 	AliInfo(Form("Retrive MaxTB histogram for %s of %5.2f GeV/c", fpartName[k], fTrackMomentum[ip]));

	return (TH1*)fHistTimeBin->At(((k==AliPID::kElectron)?0:1)*fNMom+ip);
}

//_________________________________________________________________________
Double_t AliTRDCalPIDLQ::GetProbability(Int_t spec, Float_t mom, Float_t *dedx, Float_t length) const
{
  //
	// Core function of AliTRDCalPIDLQ class for calculating the
	// likelihood for species "spec" (see AliTRDtrack::kNspecie) of a
	// given momentum "mom" and a given dE/dx (slice "dedx") yield per
	// layer
  //

	if (spec < 0 || spec >= AliPID::kSPECIES) return 0.;
		
	//Double_t dedx   = dedx1/fMeanChargeRatio;
	
	// find the interval in momentum and track segment length which applies for this data
	Int_t ilength = 1;
  while(ilength<fNLength-1 && length>fTrackSegLength[ilength]){
		ilength++;
	}
	Int_t imom = 1;
  while(imom<fNMom-1 && mom>fTrackMomentum[imom]) imom++;
	
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
Double_t AliTRDCalPIDLQ::GetProbabilityT(Int_t spec, Double_t mom, Int_t timbin) const
{
  //
  // Gets the Probability of having timbin at a given momentum (mom)
  // and particle type (spec) (0 for e) and (2 for pi)
  // from the precalculated timbin distributions 
  //
  
	if (spec < 0 || spec >= AliPID::kSPECIES) return 0.;
  if (timbin<0 || timbin >= fNTimeBins) return 0.;

  Int_t iTBin = timbin+1;
  
  // Everything which is not an electron counts as a pion for time bin max
  //if(spec != AliPID::kElectron) spec = AliPID::kPion;
  

  
	Int_t imom = 1;
  while(imom<fNMom-1 && mom>fTrackMomentum[imom]) imom++;

	Double_t mom1 = fTrackMomentum[imom-1], mom2 = fTrackMomentum[imom];
	TH1F *hist = 0x0;
	if(!(hist = (TH1F*) fHistTimeBin->At(((spec==AliPID::kElectron)?0:1)*fNMom+imom-1))){
		AliInfo(Form("Looking for spec(%d) mom(%f) timbin(%d)", spec, mom, timbin));
		AliError(Form("THistogram id %d not found in DB.", ((spec==AliPID::kElectron)?0:1)*fNMom+imom-1));
		return 0.;
	}
	Double_t LQ1 = hist->GetBinContent(iTBin);

	if(!(hist = (TH1F*) fHistTimeBin->At(((spec==AliPID::kElectron)?0:1)*fNMom+imom))){
		AliInfo(Form("Looking for spec(%d) mom(%f) timbin(%d)", spec, mom, timbin));
		AliError(Form("THistogram id %d not found in DB.", ((spec==AliPID::kElectron)?0:1)*fNMom+imom));
		return LQ1;
	}
	Double_t LQ2 = hist->GetBinContent(iTBin);

	// return interpolation over momentum binning
  return LQ1 + (LQ2 - LQ1)*(mom - mom1)/(mom2 - mom1);
}

//__________________________________________________________________
Bool_t AliTRDCalPIDLQ::WriteReferences(Char_t *File, Char_t *dir)
{
	// Build, Fill and write to file the histograms used for PID.
	// The simulations are looked in the
	// directories with the general form Form("p%3.1f", momentum)
	// starting from dir (default .). Here momentum belongs to the list
	// of known momentum to PID (see fTrackMomentum).
	// The output histograms are
	// written to the file "File" in cwd (default
	// TRDPIDHistograms.root). In order to build a DB entry
	// consider running $ALICE_ROOT/Cal/AliTRDCreateDummyCDB.C
	// 
	// Author:
	// Alex Bercuci (A.Bercuci@gsi.de)

	const Int_t nDirs = 1;
  Int_t partCode[AliPID::kSPECIES] =
    {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton};
	
	// minimal test of simulation location
	TFile *f = new TFile(Form("%s/p%3.1f/galice.root", dir, 2.));
	if(!f || f->IsZombie()){
		AliError(Form("Could not access file galice in directry \"%s/p%3.1f\". Please check the location and try again.", dir, 2.));
		return kFALSE;
	}
	f->Close(); delete f;
	f = new TFile(Form("%s/p%3.1f/AliESDs.root", dir, 2.));
	if(!f || f->IsZombie()){
		AliError(Form("Could not access file AliESDs in directry \"%s/p%3.1f\". Please check the location and try again.", dir, 2.));
		return kFALSE;
	}
	f->Close(); delete f;
	

	// Init statistics
	Int_t nPart[AliPID::kSPECIES], nTotPart;
	printf("P[GeV/c] ");
	for(Int_t ispec=0; ispec<AliPID::kSPECIES; ispec++) printf(" %s[%%] ", fpartSymb[ispec]);
	printf("\n-----------------------------------------------\n");
	
	

	// Build PID reference histograms and reference object
	const Int_t color[] = {4, 3, 2, 7, 6};
	for (Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
		if (ispec != AliPID::kElectron && ispec != AliPID::kPion) continue;
	
		h1MaxTB[(ispec>0)?1:0] = new TH1F(Form("h1%s", fpartSymb[ispec]), "", fNTimeBins, -.5, fNTimeBins-.5);
		h1MaxTB[(ispec>0)?1:0]->SetLineColor(color[ispec]);
  }
	AliTRDCalPIDLQRef ref;
	
	// momentum loop
	AliRunLoader *fRunLoader = 0x0;
	TFile *esdFile = 0x0;
	TTree *esdTree = 0x0;
	AliESD *esd = 0x0;
	AliESDtrack *esdTrack = 0x0;
	for (Int_t imom = 0; imom < fNMom; imom++) {
	//for (Int_t imom = 4; imom < 5; imom++) {
		ref.Reset();
		
		for(Int_t idir = 0; idir<nDirs; idir++){
			// open run loader and load gAlice, kinematics and header
			fRunLoader = AliRunLoader::Open(Form("%s/p%3.1f/galice.root", dir, fTrackMomentum[imom]));
			if (!fRunLoader) {
				AliError(Form("Getting run loader for momentum %3.1f failed.", fTrackMomentum[imom]));
				return kFALSE;
			}
			TString s; s.Form("%s/p%3.1f/", dir, fTrackMomentum[imom]);
			fRunLoader->SetDirName(s);
			fRunLoader->LoadgAlice();
			gAlice = fRunLoader->GetAliRun();
			if (!gAlice) {
				AliError(Form("galice object not found for momentum %3.1f.", fTrackMomentum[imom]));
				return kFALSE;
			}
			fRunLoader->LoadKinematics();
			fRunLoader->LoadHeader();
	
			// open the ESD file
			esdFile = TFile::Open(Form("%s/p%3.1f/AliESDs.root", dir, fTrackMomentum[imom]));
			if (!esdFile || esdFile->IsZombie()) {
				AliError(Form("Opening ESD file failed for momentum", fTrackMomentum[imom]));
				return kFALSE;
			}
			esd = new AliESD;
			esdTree = (TTree*)esdFile->Get("esdTree");
			if (!esdTree) {
				AliError(Form("ESD tree not found for momentum %3.1f.", fTrackMomentum[imom]));
				return kFALSE;
			}
			esdTree->SetBranchAddress("ESD", &esd);
			nTotPart = 0;
			for(Int_t ispec=0; ispec<AliPID::kSPECIES; ispec++) nPart[ispec] = 0;
	
			
			// Event loop
			for (Int_t iEvent = 0; iEvent < fRunLoader->GetNumberOfEvents(); iEvent++) {
				fRunLoader->GetEvent(iEvent);
	
				// read stack info
				AliStack* stack = gAlice->Stack();
				TArrayF vertex(3);
				fRunLoader->GetHeader()->GenEventHeader()->PrimaryVertex(vertex);
							
				// Load event summary data
				esdTree->GetEvent(iEvent);
				if (!esd) {
					AliWarning(Form("ESD object not found for event %d. [@ momentum %3.1f]", iEvent, imom));
					continue;
				}
	
				for(Int_t iTrack=0; iTrack<esd->GetNumberOfTracks(); iTrack++){
					esdTrack = esd->GetTrack(iTrack);
	
					if(!AliTRDpidESD::CheckTrack(esdTrack)) continue;
					//if((esdTrack->GetStatus() & AliESDtrack::kITSrefit) == 0) continue;
					//if(esdTrack->GetConstrainedChi2() > 1E9) continue;
					//if ((esdTrack->GetStatus() & AliESDtrack::kESDpid) == 0) continue;
					if (esdTrack->GetTRDsignal() == 0.) continue;
	
					// read MC info
					Int_t label = esdTrack->GetLabel();
					if(label<0) continue;
					if (label > stack->GetNtrack()) continue;     // background
					TParticle* particle = stack->Particle(label);
					if(!particle){
						AliWarning(Form("Retriving particle with index %d from AliStack failed. [@ event %d track %d]", label, iEvent, iTrack));
						continue;
					}
					if(particle->Pt() < 1.E-3) continue;
					//      if (TMath::Abs(particle->Eta()) > 0.3) continue;
					TVector3 dVertex(particle->Vx() - vertex[0],
										particle->Vy() - vertex[1],
										particle->Vz() - vertex[2]);
					if (dVertex.Mag() > 1.E-4){
						//AliInfo(Form("Particle with index %d generated too far from vertex. Skip from analysis. Details follows. [@ event %d track %d]", label, iEvent, iTrack));
						//particle->Print();
						continue;
					}
					Int_t iGen = -1;
					for (Int_t ispec=0; ispec<AliPID::kSPECIES; ispec++)
						if(TMath::Abs(particle->GetPdgCode()) == partCode[ispec]){
							iGen = ispec;
							break;
						}
					if(iGen<0) continue;
	
					nPart[iGen]++; nTotPart++;
					
					Float_t mom, length;
					Double_t dedx[AliTRDtrack::kNslice], dEdx;
					Int_t timebin;
					for (Int_t iPlane=0; iPlane<AliTRDgeometry::kNplan; iPlane++){
						// read data for track segment
						for(int iSlice=0; iSlice<AliTRDtrack::kNslice; iSlice++)
							dedx[iSlice] = esdTrack->GetTRDsignals(iPlane, iSlice);
						dEdx    = esdTrack->GetTRDsignals(iPlane, -1);
						timebin = esdTrack->GetTRDTimBin(iPlane);
			
						// check data
						if ((dEdx <=  0.) || (timebin <= -1.)) continue;
			
						// retrive kinematic info for this track segment
						// Temporary fix
						//if(!AliTRDpidESD::RecalculateTrackSegmentKine(esdTrack, iPlane, mom, length)) continue;
						mom = esdTrack->GetOuterParam()->GetP();

						// find segment length and momentum bin
						Int_t jmom = 1, refMom = -1;
						while(jmom<fNMom-1 && mom>fTrackMomentum[jmom]) jmom++;
						if(TMath::Abs(fTrackMomentum[jmom-1] - mom) < fTrackMomentum[jmom-1] * .2) refMom = jmom-1;
						else if(TMath::Abs(fTrackMomentum[jmom] - mom) < fTrackMomentum[jmom] * .2) refMom = jmom;
						if(refMom<0){
							AliInfo(Form("Momentum at plane %d entrance not in momentum window. [@ event %d track %d]", iPlane, iEvent, iTrack));
							continue;
						}
						/*while(jleng<fNLength-1 && length>fTrackSegLength[jleng]) jleng++;*/
						
						// this track segment has fulfilled all requierments
						//nPlanePID++;

						if(dedx[0] > 0. && dedx[1] > 0.){
							dedx[0] = log(dedx[0]); dedx[1] = log(dedx[1]);
							ref.GetPrincipal(iGen)->AddRow(dedx);
						}
						h1MaxTB[(iGen>0)?1:0]->Fill(timebin);
					} // end plane loop
				} // end track loop
			} // end events loop
			
			delete esd; esd = 0x0;
			esdFile->Close();
			delete esdFile; esdFile = 0x0;
	
			fRunLoader->UnloadHeader();
			fRunLoader->UnloadKinematics();
			delete fRunLoader; fRunLoader = 0x0;
		} // end directory loop
		
		// use data to prepare references
		ref.Prepare2DReferences();
		// save this dEdx references
		ref.SaveReferences(imom, File);
		// save MaxTB references
		SaveMaxTimeBin(imom, File);

			
		// print momentum statistics
		printf("  %3.1f  ", fTrackMomentum[imom]);
		for(Int_t ispec=0; ispec<AliPID::kSPECIES; ispec++) printf(" %5.2f ", 100.*nPart[ispec]/nTotPart);
		printf("\n");
	} // end momentum loop
	
	TFile *fSave = 0x0;
	TListIter it((TList*)gROOT->GetListOfFiles());
	while((fSave=(TFile*)it.Next()))
		if(strcmp(File, fSave->GetName())==0) break;

	fSave->cd();
	fSave->Close();
	delete fSave;

	return kTRUE;
}

//__________________________________________________________________
void	AliTRDCalPIDLQ::SaveMaxTimeBin(const Int_t mom, const char *fn)
{
  //
  // Save the histograms
  //

	TFile *fSave = 0x0;
	TListIter it((TList*)gROOT->GetListOfFiles());
	TDirectory *pwd = gDirectory;
	Bool_t kFOUND = kFALSE;
	while((fSave=(TFile*)it.Next()))
		if(strcmp(fn, fSave->GetName())==0){
			kFOUND = kTRUE;
			break;
		}
	if(!kFOUND) fSave = new TFile(fn, "RECREATE");
	fSave->cd();

	TH1 *h;
	h = (TH1F*)h1MaxTB[0]->Clone(Form("h1MaxTBEL%02d", mom));
	h->SetTitle(Form("Maximum Time Bin distribution for electrons @ %4.1f GeV", fTrackMomentum[mom]));
	h->GetXaxis()->SetTitle("time [100 ns]");
	h->GetYaxis()->SetTitle("Probability");
	h->Write();

	h = (TH1F*)h1MaxTB[1]->Clone(Form("h1MaxTBPI%02d", mom));
	h->SetTitle(Form("Maximum Time Bin distribution for pions @ %4.1f GeV", fTrackMomentum[mom]));
	h->GetXaxis()->SetTitle("time [100 ns]");
	h->GetYaxis()->SetTitle("Probability");
	h->Write();
	
	pwd->cd();
}

