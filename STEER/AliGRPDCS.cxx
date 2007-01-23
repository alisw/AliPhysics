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
//                AliGRPDCS class
//   This class deals with the DCS related info of the GRP
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

#include "AliGRPDCS.h"

#include "AliDCSValue.h"

//#include <TObjString.h>
class TObjString;

#include <TH1.h>
#include <TObjArray.h>

ClassImp(AliGRPDCS)

//_______________________________________________________________
AliGRPDCS::AliGRPDCS():
  TObject(), fDCSArray(new TObjArray()) {
  // default constructor
  
}

//_______________________________________________________________
AliGRPDCS::AliGRPDCS(TObjArray *dcsArray):
  TObject() {
  // constructor
  fDCSArray = new TObjArray();
  fDCSArray = dcsArray;
}

//___________________________________________________________________________
AliGRPDCS::AliGRPDCS(const AliGRPDCS& grpDcs):
  TObject(grpDcs) {
  //copy constructor

  if (grpDcs.fDCSArray) fDCSArray = new TObjArray();
}

//_______________________________________________________________
const char* AliGRPDCS::ProcessDCS(TH1 *h) {
  // process the dcs float values
  Float_t fFDCSArraySum = 0.0, fFDCSArrayMean = 0.0;
  for(Int_t i = 0; i < fDCSArray->GetEntries(); i++) {
    AliDCSValue *v = (AliDCSValue *)fDCSArray->At(i);
    h->Fill(v->GetFloat());
    fFDCSArraySum += v->GetFloat();
  }
  fFDCSArrayMean = fFDCSArraySum/fDCSArray->GetEntries();
  TString fDCSDataPointValue; fDCSDataPointValue += fFDCSArrayMean;

  return fDCSDataPointValue.Data();
}

