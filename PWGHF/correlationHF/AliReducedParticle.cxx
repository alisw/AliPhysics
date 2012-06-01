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

AliReducedParticle::AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, Int_t McLabel, Int_t trackid) : 
fEta(eta), 
fPhi(phi), 
fpT(pt), 
fMcLabel(McLabel), 
fid(trackid)
{
	//
	// default constructor
	//
}

AliReducedParticle::~AliReducedParticle() {

	//
	// destructor
	//
	
	Info("AliReducedParticle","Calling Destructor");
}
