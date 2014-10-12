/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

// Comment describing what this class does needed!

// $Id: AliJPiZero.cxx,v 1.3 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliJPiZero.cxx
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.3 $
  \date $Date: 2008/05/08 13:44:45 $
*/
////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include  "AliJBaseTrack.h"
#include  "AliJPhoton.h"

#include "AliJPiZero.h"



AliJPiZero::AliJPiZero():
  AliJBaseTrack(),
  fAsymm(-999),
  fMassBin(-999),
  fIsGood(kFALSE),
  fConvPlaneDelta(-999),
  fPtBin(-1)
{
  //constructor
  fPhoton[0] = 0;
  fPhoton[1] = 0;
}
//______________________________________________________________________________
AliJPiZero::AliJPiZero(const AliJPiZero& a):
  AliJBaseTrack(a),
  fAsymm(a.fAsymm),
  fMassBin(a.fMassBin),
  fIsGood(a.fIsGood),
  fConvPlaneDelta(a.fConvPlaneDelta),
  fPtBin(a.fPtBin)
{
  //copy constructor
  fPhoton[0] = a.fPhoton[0];
  fPhoton[1] = a.fPhoton[1];
}

AliJPiZero& AliJPiZero::operator=(const AliJPiZero& pizero){
  //copy constructor
  AliJBaseTrack::operator=(pizero);
  fAsymm = pizero.fAsymm;
  fMassBin = pizero.fMassBin;
  fIsGood = pizero.fIsGood;
  fPhoton[0] = pizero.fPhoton[0];
  fPhoton[1] = pizero.fPhoton[1];
	fConvPlaneDelta = pizero.fConvPlaneDelta;
	fPtBin = pizero.fPtBin;
	
	return *this;
}


//______________________________________________________________________________

bool AliJPiZero::SetMass(AliJPhoton *g1, AliJPhoton *g2){
  //set pi0 four-momentum, assddymetry, inv. mass, pT, phi, theta, energy

  fAsymm    = TMath::Abs((g1->E()-g2->E())/(g1->E()+g2->E()));
  *(TLorentzVector*)this = *(TLorentzVector*)g1 + *(TLorentzVector*)g2;
  
  fPhoton[0] = g1;
  fPhoton[1] = g2;

//   Double_t pimass;
//   pimass = TMath::Sqrt( 2. * g1->E() * g2->E() *
//                          ( 1. - TMath::Cos( g1->Angle( g2->Vect() ))));

  return true;
}

double AliJPiZero::DeltaConversion(){
    // check pair plane vector delta from expected conversion pair plain
    
    TVector3 pplain, phot1, phot2, pplaint, cplain, zplain, convplane;
		
		// pair vector in the pair plane
		pplain = Vect();
// 		cout << "-------------------" << endl;
//  		pplain.Print();
		
		// perpendicular to the plane
		phot1 = fPhoton[0]->Vect();
		phot2 = fPhoton[1]->Vect();
		pplaint = phot1.Cross( phot2 );
// 		pplaint.Print();
		
		// get the pair plane vector
		cplain = pplain.Cross( pplaint ).Unit();
// 		cplain.Print();
		
		// z plane vector
		zplain.SetXYZ( 0, 0, 1 );
// 		zplain.Print();
		
		// expected conversion plane vector
		convplane = pplain.Cross( zplain ).Unit();
// 		convplane.Print();
		
		// get the plane vector cosine similarity
		// both vectors must be unit for this to work
		fConvPlaneDelta = TMath::Abs( TMath::ACos( cplain.Dot( convplane )));
// 		cout << fConvPlaneDelta << endl;
		
		// the result comes out in 0,pi range. Depending on how the z plane vector
		// was chosen (2 options). Merge the results into 0,1/2pi range
		if( fConvPlaneDelta > TMath::Pi() / 2. )
			fConvPlaneDelta = TMath::Pi() - fConvPlaneDelta;
		
    return fConvPlaneDelta;    
}



//______________________________________________________________________________

// double AliJPiZero::operator- (const AliJPiZero &pi0){
// 
//   //========================================
//   // retruns the closest distance between
//   // photon hits bolonging to two pi0
//   //========================================
// 
//   TLorentzVector  v;
//   double d = 1e3;
//   v = fPhoton1 - pi0.fPhoton1;
//   if(v.M() < d) d = v.M();
//   v = fPhoton2 - pi0.fPhoton2;
//   if(v.M() < d) d = v.M();
//   v = fPhoton1 - pi0.fPhoton2;
//   if(v.M() < d) d = v.M();
//   v = fPhoton2 - pi0.fPhoton2;
//   if(v.M() < d) d = v.M();
//   return d;
// }
// 
// //______________________________________________________________________________ 
// AliJPiZero& AliJPiZero::operator= (const AliJPiZero& piz){
//   //operator=
//   if(this != &piz){
//     AliJBaseTrack::operator=(piz);
//     fPhoton1 = piz.fPhoton1;
//     fPhoton2 = piz.fPhoton2;
//     PhotSum  = fPhoton1 + fPhoton2;
//     fAsymm   = piz.fAsymm;
//     fMassBin = piz.fMassBin;
//   }
// 
//   return *this;
// }



ClassImp(AliJPiZero)
