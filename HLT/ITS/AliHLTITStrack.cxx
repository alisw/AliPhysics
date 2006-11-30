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

//-------------------------------------------------------------------------
//                Implementation of the HLT ITS track class
//
//          Origin: Cvetan Cheshkov, CERN, Cvetan.Cheshkov@cern.ch
//-------------------------------------------------------------------------

#include <TMath.h>

#include "AliHLTITStrack.h"

ClassImp(AliHLTITStrack)

//____________________________________________________________________________
AliHLTITStrack::AliHLTITStrack()
              :AliITStrackV2()
{
  //------------------------------------------------------------------
  //Constructor
  //------------------------------------------------------------------
}

//____________________________________________________________________________
AliHLTITStrack::AliHLTITStrack(AliESDtrack& t)
              :AliITStrackV2(t)
{
  //------------------------------------------------------------------
  //Constructor
  //------------------------------------------------------------------
}

//____________________________________________________________________________
AliHLTITStrack::AliHLTITStrack(const AliHLTITStrack& t) 
              : AliITStrackV2(t)
{
  //------------------------------------------------------------------
  //Copy constructor
  //------------------------------------------------------------------
}

//_____________________________________________________________________________
Int_t AliHLTITStrack::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the their curvature
  //-----------------------------------------------------------------
  AliHLTITStrack *t=(AliHLTITStrack*)o;
  Double_t co=TMath::Abs(t->Get1Pt());
  Double_t c =TMath::Abs(Get1Pt());
  //  Double_t co=t->GetSigmaY2()*t->GetSigmaZ2();
  //  Double_t c =GetSigmaY2()*GetSigmaZ2();
  if (c>co) return 1;
  else if (c<co) return -1;
  return 0;
}
