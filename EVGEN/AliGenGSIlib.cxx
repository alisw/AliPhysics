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

/*
$Log$
Revision 1.7  2000/05/02 08:12:13  morsch
Coding rule violations corrected.

Revision 1.6  1999/09/29 09:24:14  fca
Introduction of the Copyright and cvs Log

*/

#include "AliGenGSIlib.h"
#include "AliRun.h"
#include "iostream.h"

ClassImp(AliGenGSIlib)

  Bool_t AliGenGSIlib::fgDebug =kFALSE;
//                      Upsilon
//
//
//                  pt-distribution
//____________________________________________________________
Double_t AliGenGSIlib::PtUpsilonRitman( Double_t *px, Double_t *dummy )
{
  // Upsilon pT
  /*   AliGenMUONlib parametrisation
  const Double_t kpt0 = 5.3;
  const Double_t kxn  = 2.5;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
  */
  if (fgDebug) cout<<"Ritman Pt paramtrisation\n";
  const Double_t kpt0 = 4.7;
  const Double_t kxn  = 3.5;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+((x*x)/(kpt0*kpt0));
  return x/TMath::Power(pass1,kxn);
   
}
//
//                    y-distribution
//
//____________________________________________________________
Double_t AliGenGSIlib::YUpsilonRitman(Double_t *py, Double_t *dummy)
{

    // Upsilon y

/*   original from AliGenMUON
    const Double_t ky0 = 3.;
    const Double_t kb=1.;
    Double_t yu;
    Double_t y=TMath::Abs(*py);
    //
    if (y < ky0)
    yu=kb;
    else
    yu=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);
    return yu;
    */
  if (fgDebug) cout<<"Ritman Y paramtrisation\n";
  return 0.003;  //GSI parametrisation 
}
//                 particle composition
//
Int_t AliGenGSIlib::IpUpsilonRitman()
{
// y composition
  if (fgDebug) cout<<"Ritman Ip paramtrisation\n";
  return 553;     
}


Double_t AliGenGSIlib::PtUpsilonKarel( Double_t *px, Double_t *dummy )
{
  // Upsilon pT
  //to implement
  if (fgDebug) cout<<"Karel Pt paramtrisation\n";

  return 0.;   
}
//
//                    y-distribution
//
//____________________________________________________________
Double_t AliGenGSIlib::YUpsilonKarel(Double_t *py, Double_t *dummy)
{
  
  // Upsilon y
//to implement
  if (fgDebug) cout<<"Karel Y paramtrisation\n";
  return 0.003;  //Karel parametrisation 
}

//                 particle composition
//

Int_t AliGenGSIlib::IpUpsilonKarel()
{
  // y composition//
  //to implement
  if (fgDebug) cout<<"Karel Ip paramtrisation\n";
  return 553;     
}



Double_t AliGenGSIlib::PtUpsilonMUON( Double_t *px, Double_t *dummy )
{
// Upsilon pT
  const Double_t kpt0 = 5.3;
  const Double_t kxn  = 2.5;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
//
//                    y-distribution
//
//____________________________________________________________
Double_t AliGenGSIlib::YUpsilonMUON(Double_t *py, Double_t *dummy)
{
// Upsilon y
  const Double_t ky0 = 3.;
  const Double_t kb=1.;
  Double_t yu;
  Double_t y=TMath::Abs(*py);
  //
  if (y < ky0)
    yu=kb;
  else
    yu=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);
  return yu;
}
//                 particle composition
//
Int_t AliGenGSIlib::IpUpsilonMUON()
{
// y composition
    return 553;
}


typedef Double_t (*GenFunc) (Double_t*,  Double_t*);

typedef Int_t (*GenFuncIp) ();



GenFunc AliGenGSIlib::GetPt(Param_t param, const char * tname)
{
// Return pointer to pT parameterisation
  GenFunc func=0;
  TString sname(tname);
  switch (param) 
    {
    case upsilon_p:
      if (sname=="MUON"){
	func= PtUpsilonMUON;
	break;
      }
      if (sname=="RITMAN"){
	func=PtUpsilonRitman;
	break;
      }
      if (sname=="KAREL"){
	func=PtUpsilonKarel;
	break;
      }
      break;
       func=0;
        printf("<AliGenGSIlib::GetPt> unknown parametrisation\n");
    default:
        func=0;
        printf("<AliGenGSIlib::GetPt> unknown parametrisation\n");
    }
    return func;
}



GenFunc AliGenGSIlib::GetY(Param_t param, const char * tname)
{
  // Return pointer to y- parameterisation
   GenFunc func=0;
    TString sname(tname);
    switch (param) 
    {
    case upsilon_p:
      if (sname=="MUON"){
	func= YUpsilonMUON;
	break;
      }
      if (sname=="RITMAN"){
	func=YUpsilonRitman;
	break;
      }
      if (sname=="KAREL"){
	func=YUpsilonKarel;
	break;
      }
      func=0;
      printf("<AliGenGSIlib::GetY> unknown parametrisation\n");
      break;
    default:
        func=0;
        printf("<AliGenGSIlib::GetY> unknown parametrisation\n");
    }
    return func;
}



GenFuncIp AliGenGSIlib::GetIp(Param_t param, const char * tname)
{
// Return pointer to particle type parameterisation
    GenFuncIp func=0;
    TString sname(tname);
    switch (param) 
    {
    case upsilon_p:
      if (sname=="MUON"){
	func=IpUpsilonMUON;
	break;
      }
      if (sname=="RITMAN"){
	func=IpUpsilonRitman;
	break;
      }
      if (sname=="KAREL"){
	func=IpUpsilonKarel;
	break;
      }
      func=0;
      printf("<AliGenGSIlib::GetIP> unknown parametrisation\n");
      break;
    default:
        func=0;
        printf("<AliGenGSIlib::GetIp> unknown parametrisation\n");
    }
    return func;
}













