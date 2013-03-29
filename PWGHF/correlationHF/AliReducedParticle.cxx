/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
//
//
//             Base class for DStar - Hadron Correlations Analysis
//
//-----------------------------------------------------------------------
//          
//
//						   Author S.Bjelogrlic
//                         Utrecht University 
//                      sandro.bjelogrlic@cern.ch
//
//-----------------------------------------------------------------------

/* $Id$ */

#include "AliVParticle.h"
#include "AliReducedParticle.h"

AliReducedParticle::AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, Int_t mcLabel, Int_t trackid, Double_t impPar, Bool_t checkSoftPi) : 
fEta(eta), 
fPhi(phi), 
fpT(pt), 
fMcLabel(mcLabel), 
fid(trackid),
fImpPar(impPar),
fCheckSoftPi(checkSoftPi),
fCharge(0),
fInvMass(0),
fWeight(1.),
fPtBin(-1),
fOriginMother(-1)
{
	//
	// default constructor
	//
}

AliReducedParticle::AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, Int_t mcLabel, Int_t trackid, Double_t impPar, Bool_t checkSoftPi, Short_t charge) : 
fEta(eta), 
fPhi(phi), 
fpT(pt), 
fMcLabel(mcLabel), 
fid(trackid),
fImpPar(impPar),
fCheckSoftPi(checkSoftPi),
fCharge(charge),
fInvMass(0),
fWeight(1.),
fPtBin(-1),
fOriginMother(-1)
{
	//
	// default constructor
	//
}

AliReducedParticle::AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, Int_t mcLabel, Int_t trackid, Double_t impPar, Bool_t checkSoftPi, Short_t charge, Double_t weight) : 
fEta(eta), 
fPhi(phi), 
fpT(pt), 
fMcLabel(mcLabel), 
fid(trackid),
fImpPar(impPar),
fCheckSoftPi(checkSoftPi),
fCharge(charge),
fInvMass(0),
fWeight(weight),
fPtBin(-1),
fOriginMother(-1)
{
	//
	// default constructor
	//
}



AliReducedParticle::AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, Int_t McLabel) : 
fEta(eta), 
fPhi(phi), 
fpT(pt), 
fMcLabel(McLabel),
fid(0),
fImpPar(0.),
fCheckSoftPi(kFALSE),
fCharge(0),
fInvMass(0),
fWeight(1.),
fPtBin(-1),
fOriginMother(-1)
{
	//
	// default constructor
	//
}

AliReducedParticle::AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, int charge, int originmother) : 
fEta(eta), 
fPhi(phi), 
fpT(pt), 
fMcLabel(0),
fid(0),
fImpPar(0.),
fCheckSoftPi(kFALSE),
fCharge(charge),
fInvMass(0),
fWeight(1.),
fPtBin(-1),
fOriginMother(originmother)
{
	//
	// default constructor
	//
}

AliReducedParticle::AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, Double_t invmass, int ptbin, int originmother) : 
fEta(eta), 
fPhi(phi), 
fpT(pt), 
fMcLabel(0),
fid(0),
fImpPar(0.),
fCheckSoftPi(kFALSE),
fCharge(0),
fInvMass(invmass),
fWeight(1.),
fPtBin(ptbin),
fOriginMother(originmother)
{
	//
	// default constructor
	//
}


AliReducedParticle::~AliReducedParticle() {

	//
	// destructor
	//
	
	//Info("AliReducedParticle","Calling Destructor");
}
