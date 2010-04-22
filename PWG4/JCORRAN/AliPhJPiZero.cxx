/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notifce   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id: AliPhJPiZero.cxx,v 1.3 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
//
//  \file AliPhJPiZero.cxx
//  \brief
//  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
//  \email: djkim@jyu.fi
//  \version $Revision: 1.3 $
//  \date $Date: 2008/05/08 13:44:45 $
//
//  Class encapsulates basic properties of  pi0
////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include  "AliPhJBaseTrack.h"
#include  "AliPhJPhoton.h"

#include "AliPhJPiZero.h"

ClassImp(AliPhJPiZero)

//___________________________________________________________
AliPhJPiZero::AliPhJPiZero():
  fV1(-999,-999,-999), 
  fV2(-999,-999,-999),
  fPi0P(-999,-999,-999),
  fPhotSum(0,0,0,0),
  fPhoton1(0,0,0,0),
  fPhoton2(0,0,0,0),
  fAsymm(-999),
  fPizM(-999),
  fMassBin(-999)
{
  //constructor
}
//______________________________________________________________________________
AliPhJPiZero::AliPhJPiZero(const AliPhJPiZero& a):
  AliPhJBaseTrack(a),
  fV1(a.fV1),
  fV2(a.fV2),
  fPi0P(a.fPi0P),
  fPhotSum(a.fPhotSum),
  fPhoton1(a.fPhoton1),
  fPhoton2(a.fPhoton2),
  fAsymm(a.fAsymm),
  fPizM(a.fPizM),
  fMassBin(a.fMassBin)
{
  //copy constructor
}

//______________________________________________________________________________

bool AliPhJPiZero::SetMass(AliPhJPhoton *g1, AliPhJPhoton *g2){
  //set pi0 four-momentum, assymetry, inv. mass, pT, phi, theta, energy
  fV1.SetXYZ(g1->GetX(), g1->GetY(), g1->GetZ() );  //Z is corrected in main  for vertex Z
  fPhoton1.SetPxPyPzE(g1->GetE()*sin(fV1.Theta())*cos(fV1.Phi()),
		     g1->GetE()*sin(fV1.Theta())*sin(fV1.Phi()),
		     g1->GetE()*cos(fV1.Theta()),
		     g1->GetE());
 
  fV2.SetXYZ(g2->GetX(), g2->GetY(), g2->GetZ() ); 
  fPhoton2.SetPxPyPzE(g2->GetE()*sin(fV2.Theta())*cos(fV2.Phi()),
		     g2->GetE()*sin(fV2.Theta())*sin(fV2.Phi()),
		     g2->GetE()*cos(fV2.Theta()),
		     g2->GetE());
  fPhotSum  = fPhoton1 + fPhoton2;
  fPi0P     = fPhotSum.Vect();
  fAsymm    = fabs((g1->GetE()-g2->GetE())/(g1->GetE()+g2->GetE()));
  fPizM     = fPhotSum.M();
  fBasePt   = fPhotSum.Pt();
  fBasePhi  = fPhotSum.Phi();
  fBaseTheta= fPhotSum.Theta();
  fBaseE    = fPhotSum.E();
  return true;
  
}

//______________________________________________________________________________

double AliPhJPiZero::operator- (const AliPhJPiZero &pi0){
  
    //========================================
    // retruns the closest distance between
    // photon hits bolonging to two pi0
    //========================================

  TVector3  v;
  double d = 1e3;
  v = fV1-pi0.fV1;
  if(v.Mag() < d) d = v.Mag();
  v = fV2-pi0.fV1;
  if(v.Mag() < d) d = v.Mag();
  v = fV1-pi0.fV2;
  if(v.Mag() < d) d = v.Mag();
  v = fV2-pi0.fV2;
  if(v.Mag() < d) d = v.Mag();
  return d;
}


//______________________________________________________________________________

AliPhJPiZero& AliPhJPiZero::operator= (const AliPhJPiZero& piz){
  //operator=
  if (this != &piz){
    AliPhJBaseTrack::operator=(piz);
    fV1      = piz.fV1;
    fV2      = piz.fV2;
    fPi0P    = piz.fPi0P;
    fPhotSum = piz.fPhotSum;
    fPhoton1 = piz.fPhoton1;
    fPhoton2 = piz.fPhoton2;
    fAsymm   = piz.fAsymm;
    fPizM    = piz.fPizM;
    fMassBin = piz.fMassBin;
  }

  return *this;
}


