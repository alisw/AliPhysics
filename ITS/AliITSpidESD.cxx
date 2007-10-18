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

//-----------------------------------------------------------------
//           Implementation of the ITS PID class
// Very naive one... Should be made better by the detector experts...
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include <TMath.h>

#include "AliITSpidESD.h"
#include "AliESDtrack.h"

ClassImp(AliITSpidESD)


//______________________________________________________________________
AliITSpidESD::AliITSpidESD():TObject(){
  //Default constructor
 
}

Double_t AliITSpidESD::Bethe(Double_t p,Double_t mass) {

  Double_t mom=p*1000;//MeV
  Double_t Mass=mass*1000;//Mev
  Float_t dens =2.33; //g cm-3
  Double_t K=0.307075;//MeVcm^2/g
  Double_t ZovA=0.49848;
  Double_t me=0.511;//MeV/c^2
  Double_t I=173./1000000.;//MeV
  Double_t En=TMath::Sqrt(mom*mom+Mass*Mass);//MeV
  Double_t gamma=En/Mass;
  Double_t beta=mom/En;
  Double_t Tmax=2*me*beta*beta*gamma*gamma/(1+2*gamma*me/Mass+(me/Mass)*(me/Mass));
  Double_t deltaover2=28.816*1e-6*TMath::Sqrt(dens*ZovA)+TMath::Log(beta*gamma)-0.5;
  Double_t FNor=0.009164; //normalizing to 1 at the minimum of ionization

  return K*ZovA*1/(beta*beta)*(0.5*TMath::Log(2*me*beta*beta*gamma*gamma*Tmax/(I*I))-beta*beta-deltaover2)*2.33*1000*0.03*FNor;


}
