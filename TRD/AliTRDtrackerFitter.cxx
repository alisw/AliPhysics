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
//  Track fitter                                                             //
//                                                                           //
//  Authors:                                                                 //
//    Alex Bercuci <A.Bercuci@gsi.de>                                        //
//    Markus Fasel <M.Fasel@gsi.de>                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TLinearFitter.h>
#include <TMath.h>
#include <TTreeStream.h>

#include "AliLog.h"

#include "AliTRDseed.h"
#include "AliTRDcalibDB.h"
#include "AliTRDcluster.h"
#include "AliTRDReconstructor.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackerFitter.h"

#define DEBUG

ClassImp(AliTRDtrackerFitter)

//_________________________________________________________________________
AliTRDtrackerFitter::AliTRDtrackerFitter()
  :TObject()
  ,fFitterTC(0x0)
  ,fFitterT2(0x0)
  ,fRieman1(0x0)
  ,fRieman2(0x0)
  ,fChi2TR(0)
  ,fChi2TC(0)
  ,fCR(0)
  ,fCC(0)
  ,fDca(0)
  ,fDzmf(0)
  ,fZmf(0)
  ,fNlayers(0)
  ,fDebugStream(0x0)
{
  //
  // Constructor
  //

	fRieman1   = new AliRieman(1000);
	fRieman2  = new AliRieman(1000);
	fFitterTC = new TLinearFitter(2,"hyp2");
	fFitterT2 = new TLinearFitter(4,"hyp4");
	fFitterTC->StoreData(kTRUE);
	fFitterT2->StoreData(kTRUE);
}

//_________________________________________________________________________
AliTRDtrackerFitter::AliTRDtrackerFitter(const AliTRDtrackerFitter &f)
  :TObject(f)
  ,fFitterTC(0x0)
  ,fFitterT2(0x0)
  ,fRieman1(0x0)
  ,fRieman2(0x0)
  ,fChi2TR(f.fChi2TR)
  ,fChi2TC(f.fChi2TC)
  ,fCR(f.fCR)
  ,fCC(f.fCC)
  ,fDca(f.fDca)
  ,fDzmf(f.fDzmf)
  ,fZmf(f.fZmf)
  ,fNlayers(f.fNlayers)
  ,fDebugStream(0x0)
{
  //
  // Copy Constructor (performs a deep copy) 
  //

	fRieman1  = new AliRieman(*f.fRieman1);
	fRieman2  = new AliRieman(*f.fRieman2);
	fFitterTC = new TLinearFitter(*f.fFitterTC);
	fFitterT2 = new TLinearFitter(*f.fFitterT2);
	fFitterTC->StoreData(kTRUE);
	fFitterT2->StoreData(kTRUE);
}

//_________________________________________________________________________
AliTRDtrackerFitter::~AliTRDtrackerFitter()
{
  //
  // Destructor
  //

	delete fRieman1;
	delete fRieman2;
	delete fFitterTC;
	delete fFitterT2;
}

//_________________________________________________________________________
AliTRDtrackerFitter &AliTRDtrackerFitter::operator=(const AliTRDtrackerFitter& fitter)
{
  //
  // Assignment operator using the copy Method of this class
  //
	if(this != &fitter){
		fitter.Copy(*this);
	}
	return *this;
}

//_________________________________________________________________________
void AliTRDtrackerFitter::Copy(TObject &f) const 
{
  //
  // Copy method. 
  // Performs a deep copy. Mainly used in the Assignment operator
  //

  AliTRDtrackerFitter &tf = (AliTRDtrackerFitter &) f;

	if(tf.fFitterTC)
		delete tf.fFitterTC;
	tf.fFitterTC = new TLinearFitter(*fFitterTC); 
	if(tf.fFitterT2)
		delete tf.fFitterT2;
	tf.fFitterT2 = new TLinearFitter(*fFitterT2);
	if(tf.fRieman1)
		delete tf.fRieman1;
	tf.fRieman1 = new AliRieman(*fRieman1);
	if(tf.fRieman2)
		delete fRieman2;
	tf.fRieman2 = new AliRieman(*fRieman2);
	tf.fChi2TR = fChi2TR;
	tf.fChi2TC = fChi2TC;
	tf.fCR = fCR;
	tf.fCC = fCC;
	tf.fDca =fDca;
	tf.fDzmf = fDzmf;
	tf.fZmf =fZmf;
	tf.fNlayers = fNlayers;
	tf.fDebugStream = 0x0;
	fFitterTC->StoreData(kTRUE);
	fFitterT2->StoreData(kTRUE);
}

//_________________________________________________________________________
void AliTRDtrackerFitter::FitRieman(AliTRDcluster **cl, Int_t nLayers)
{
  //
  // Performs a Rieman Fit including 4 clusters (one in each layer)
  //
	fRieman1->Reset();
	for(Int_t i = 0; i < nLayers; i++)
		fRieman1->AddPoint(cl[i]->GetX(), cl[i]->GetY(), cl[i]->GetZ(), 1, 10);
	fRieman1->Update();
}

//_________________________________________________________________________
void AliTRDtrackerFitter::FitRieman(AliTRDseedV1 *cseed, Int_t *planes)
{
  //
  // Performs a Rieman fit using four respectively six seeds
  // 2 times used: 1st with 4 parameters and Later in the full track
  // fitting Part with 6 parameters
  //

	Int_t lplanes[] = {0, 1, 2, 3, 4, 5}, *tplanes;
	Int_t nplanes;
	if(planes){
		nplanes = 4;
		tplanes = planes;
	} else {
		nplanes = 6;
		tplanes = &lplanes[0];
	}
	fRieman1->Reset();
	for(Int_t iLayer = 0; iLayer < nplanes; iLayer++){
		if(!cseed[tplanes[iLayer]].IsOK()) continue;
		fRieman1->AddPoint(cseed[tplanes[iLayer]].GetX0(), cseed[tplanes[iLayer]].GetYfitR(0), cseed[tplanes[iLayer]].GetZProb(), 1, 10);
	}
	fRieman1->Update();
}

//_________________________________________________________________________
Double_t AliTRDtrackerFitter::FitHyperplane(AliTRDseedV1 *cseed
                                          , Double_t chi2ZF
                                          , Double_t zFitter)
{
  //
  // performs a hyperplane fit
  // 
  // Two linear fitters:  one in which kz is fixed and one in which kz is a free parameter
  // checking if values are acceptable and improves the fitresults
  // Calculates curvature parameters and chisquares
  // Return: likelihood value as a measurement for the seedquality
  //

  //const Int_t kMaxTimebins = 100; // Buffer
	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	Int_t nTimeBins = cal->GetNumberOfTimeBins();
	Int_t npointsT = 0;
	fFitterTC->ClearPoints();
	fFitterT2->ClearPoints();
	
	Double_t fHl[6];
	Double_t xref2 = (cseed[2].GetX0() + cseed[3].GetX0())/2;
	fRieman2->Reset();
	for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
		if (!cseed[iLayer].IsOK()) continue;
		fHl[iLayer] = cseed[iLayer].GetTilt();

		for (Int_t itime = 0; itime < nTimeBins; itime++) {
			if (!cseed[iLayer].IsUsable(itime)) continue;

			// X relative to the middle chamber
			Double_t x  = cseed[iLayer].GetX(itime) + cseed[iLayer].GetX0() - xref2;
			Double_t y  = cseed[iLayer].GetY(itime);
			Double_t z  = cseed[iLayer].GetZ(itime);
			// ExB correction to the correction
			// Tilted fRieman1
			Double_t uvt[6];
			// Global x
			Double_t x2 = cseed[iLayer].GetX(itime) + cseed[iLayer].GetX0();
			Double_t t  = 1.0 / (x2*x2 + y*y);
			uvt[1] = t;                 // t
			uvt[0] = 2.0 * x2 * uvt[1]; // u
			uvt[2] = 2.0 * fHl[iLayer] * uvt[1];
			uvt[3] = 2.0 * fHl[iLayer] * x * uvt[1];
			uvt[4] = 2.0 * (y + fHl[iLayer]*z) * uvt[1];
			Double_t error = 2.0 * 0.2 * uvt[1];
			fFitterT2->AddPoint(uvt,uvt[4],error);

			//Constrained fRieman1
			z = cseed[iLayer].GetZ(itime);
			uvt[0] = 2.0 * x2 * t; // u
			uvt[1] = 2.0 * fHl[iLayer] * x2 * uvt[1];
			uvt[2] = 2.0 * (y + fHl[iLayer] * (z - zFitter)) * t;	// parameter that is most propably wrong
			fFitterTC->AddPoint(uvt,uvt[2],error);
			fRieman2->AddPoint(x2,y,z,1,10);
			npointsT++;
		}
	} // Loop: iLayer

	fRieman2->Update();
	fFitterTC->Eval();
	fFitterT2->Eval();
	Double_t rpolz0 = fFitterT2->GetParameter(3);
	Double_t rpolz1 = fFitterT2->GetParameter(4);

	// Linear fitter  - not possible to make boundaries
	// Do not accept non possible z and dzdx combinations
	Bool_t acceptablez = kTRUE;
	for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
		if(!cseed[iLayer].IsOK()) continue;
		Double_t zT2 = rpolz0 + rpolz1 * (cseed[iLayer].GetX0() - xref2);
		if (TMath::Abs(cseed[iLayer].GetZProb() - zT2) > cseed[iLayer].GetPadLength() * 0.5 + 1.0) acceptablez = kFALSE;
	}
	if (!acceptablez) {
		fFitterT2->FixParameter(3,fRieman1->GetZat(xref2));
		fFitterT2->FixParameter(4,fRieman1->GetDZat(xref2));
		fFitterT2->Eval();
		fFitterT2->ReleaseParameter(3);
		fFitterT2->ReleaseParameter(4);
		rpolz0 = fFitterT2->GetParameter(3);
		rpolz1 = fFitterT2->GetParameter(4);
	}

	fChi2TR = fFitterT2->GetChisquare() / Float_t(npointsT);
	fChi2TC = fFitterTC->GetChisquare() / Float_t(npointsT);
	Double_t polz1c = fFitterTC->GetParameter(2);
	Double_t polz0c = polz1c * xref2;
	Double_t aC     =  fFitterTC->GetParameter(0);
	Double_t bC     =  fFitterTC->GetParameter(1);
	fCC     =  aC / TMath::Sqrt(bC * bC + 1.0); // Curvature
	Double_t aR     =  fFitterT2->GetParameter(0);
	Double_t bR     =  fFitterT2->GetParameter(1);
	Double_t dR     =  fFitterT2->GetParameter(2);
	fCR     =  1.0 + bR*bR - dR*aR;
	fDca    =  0.0;
	if (fCR > 0.0) {
		fDca = -dR / (TMath::Sqrt(1.0 + bR*bR - dR*aR) + TMath::Sqrt(1.0 + bR*bR));
		fCR  =  aR / TMath::Sqrt(fCR);
	}

	
	Double_t chi2ZT2 = 0.0;
	Double_t chi2ZTC = 0.0;
	for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
		if(!cseed[iLayer].IsOK()) continue;
		Double_t zT2 = rpolz0 + rpolz1 * (cseed[iLayer].GetX0() - xref2);
		Double_t zTC = polz0c + polz1c * (cseed[iLayer].GetX0() - xref2);
		chi2ZT2 += TMath::Abs(cseed[iLayer].GetMeanz() - zT2);
		chi2ZTC += TMath::Abs(cseed[iLayer].GetMeanz() - zTC);
	}
	chi2ZT2 /= TMath::Max((fNlayers - 3.0),1.0);
	chi2ZTC /= TMath::Max((fNlayers - 3.0),1.0);
	 
	
	AliTRDseedV1::FitRiemanTilt(cseed, kTRUE);
	Float_t sumdaf = 0.0;
	for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
		if(!cseed[iLayer].IsOK()) continue;
		sumdaf += TMath::Abs((cseed[iLayer].GetYfit(1) - cseed[iLayer].GetYref(1))/ cseed[iLayer].GetSigmaY2());
	}
	sumdaf /= Float_t (fNlayers - 2.0);

	// Cook Likelihoods
	Double_t likezf     = TMath::Exp(-chi2ZF * 0.14);
	Double_t likechi2C  = TMath::Exp(-fChi2TC * 0.677);
	Double_t likechi2TR = TMath::Exp(-fChi2TR * 0.78);
	Double_t likeaf     = TMath::Exp(-sumdaf * 3.23);
	Double_t likef = likezf * likechi2TR * likeaf;

#ifdef DEBUG
	if(fDebugStream && AliTRDReconstructor::StreamLevel() >= 2){
		TTreeSRedirector &treeStreamer = *fDebugStream;
		treeStreamer << "FitHyperplane"
			<< "seed0.="     << &cseed[0]
			<< "seed1.="     << &cseed[1]
			<< "seed2.="     << &cseed[2]
			<< "seed3.="     << &cseed[3]
			<< "seed4.="     << &cseed[4]
			<< "seed5.="     << &cseed[5]
			<< "chi2TR="     << fChi2TR
			<< "chi2TC="     << fChi2TC
			<< "chi2ZT2="    << chi2ZT2
			<< "chi2ZTC="    << chi2ZTC
			<< "CC="         << fCC
			<< "CR="         << fCR
			<< "DR="         << dR
			<< "DCA="        << fDca
			<< "Polz0="      << polz0c
			<< "Polz1="      << polz1c
			<< "RPolz0="     << rpolz0
			<< "RPolz1="     << rpolz1
			<< "Likechi2C="  << likechi2C
			<< "Likechi2TR=" << likechi2TR
			<< "Likezf="     << likezf
			<< "Likeaf="     << likeaf
			<< "LikeF="      << likef
			<< "Rieman1.="   << fRieman1
			<< "Rieman2.="   << fRieman2
			<< "\n";
	}
#endif

	return likef;
}

//_________________________________________________________________________
void AliTRDtrackerFitter::GetHyperplaneFitChi2(Double_t *chisquares) const
{
  //
  // Getter method returns the chisquares of the hyperplane fit
  //

	chisquares[0] = fChi2TR;
	chisquares[1] = fChi2TC;
}

//_________________________________________________________________________
void AliTRDtrackerFitter::GetHyperplaneFitResults(Double_t *params) const
{
  //
  // Getter Method returning the Curvature parameters of the hyperplane fit
  //

	params[0] = fCC;
	params[1] = fCR;
	params[2] = fDca;
}

//_________________________________________________________________________
void AliTRDtrackerFitter::Reset()
{
  //
  // Resets the object
  //

	fRieman1->Reset();
	fRieman2->Reset();
	//fFitterTC->Clear();
	//fFitterT2->Clear();
}
