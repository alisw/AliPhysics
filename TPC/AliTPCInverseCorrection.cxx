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

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliTPCInverseCorrection class                                              //
//                                                                            //
// This is a wrapper that inverts an AliTPCCorrection. This is done by        //
// swapping the CalculateCorrection and CalculateInverseCorrection functions. //
// The wrapped correction is supplied as a pointer and the class relies       //
// on the fact, that this pointer keeps pointing to the right object.         //
// However, the ownership is not changed, i.e. the wrapped correction         //
// will not be deleted when this correction is destructed.                    //
//                                                                            //
// date: 27/04/2010                                                           //
// Authors: Magnus Mager, Stefan Rossegger, Jim Thomas                        //
////////////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include "AliTPCInverseCorrection.h"

AliTPCInverseCorrection::AliTPCInverseCorrection()
  : fCorrection(0) {
  //
  // default constructor
  // (only meant for ROOT I/O)
  //
}

AliTPCInverseCorrection::AliTPCInverseCorrection(AliTPCCorrection *correction) 
  : fCorrection(correction) {
  //
  // Constructor that is creating the inverse of the supplied correction.
  // It automatically sets the name ("inv_[correction name]") and tile
  // ("Inverse of [correction title]").
  //
  TString name,title;
  name  ="inv_";
  name +=correction->GetName();
  title ="Inverse of ";
  title+=correction->GetTitle();
  SetName(name.Data());
  SetTitle(title.Data());
}

AliTPCInverseCorrection::~AliTPCInverseCorrection() {
  //
  // virtual destructor
  //
}

void AliTPCInverseCorrection::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // This is just calling the CalculateInverseCorrection of the wrapped
  // correction -- or puts dr=0 if the latter is 0.
  //
  if (fCorrection)
    fCorrection->GetDistortion(x,roc,dx);
  else
    for (Int_t j=0;j<3;++j) dx[j]=0.;
}

void AliTPCInverseCorrection::GetDistortion(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // This is just calling the CalculateCorrection of the wrapped
  // correction -- or puts dr=0 if the latter is 0.
  //
  if (fCorrection)
    fCorrection->GetCorrection(x,roc,dx);
  else
    for (Int_t j=0;j<3;++j) dx[j]=0.;
}

ClassImp(AliTPCInverseCorrection)
