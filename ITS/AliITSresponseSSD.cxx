/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
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


#include "AliITSresponseSSD.h"
//////////////////////////////////////////////////
//  Base response class for ITS                 //
//  It is used to set static data members       //
//  connected to parameters equal for all       //
//  the SSD modules                             //
//////////////////////////////////////////////////


const Float_t AliITSresponseSSD::fgkDiffCoeffDefault = 0.;
const TString AliITSresponseSSD::fgkOption1Default = "";
const TString AliITSresponseSSD::fgkOption2Default = "";
const Double_t AliITSresponseSSD::fgkfCouplingPR = 0.01;
const Double_t AliITSresponseSSD::fgkfCouplingPL = 0.01;
const Double_t AliITSresponseSSD::fgkfCouplingNR = 0.01;
const Double_t AliITSresponseSSD::fgkfCouplingNL = 0.01;
const Int_t AliITSresponseSSD::fgkZSThreshold = 3;


ClassImp(AliITSresponseSSD)

//______________________________________________________________________
AliITSresponseSSD::AliITSresponseSSD():AliITSresponse(),
fADCpereV(0),
fCouplingPR(0),
fCouplingPL(0),
fCouplingNR(0),
fCouplingNL(9),
fZSThreshold(0),
fOption1(),
fOption2(){
    // Default Constructor

  SetDiffCoeff(fgkDiffCoeffDefault,0.);
  SetParamOptions(fgkOption1Default.Data(),fgkOption2Default.Data());
  SetADCpereV();
  SetCouplings(fgkfCouplingPR,fgkfCouplingPL,fgkfCouplingNR,fgkfCouplingNL);
  SetZSThreshold(fgkZSThreshold);
}

//______________________________________________________________________
AliITSresponseSSD::AliITSresponseSSD(const AliITSresponseSSD &ob) : AliITSresponse(ob),
fADCpereV(ob.fADCpereV),
fCouplingPR(ob.fCouplingPR),
fCouplingPL(ob.fCouplingPL),
fCouplingNR(ob.fCouplingNR),
fCouplingNL(ob.fCouplingNL),
fZSThreshold(ob.fZSThreshold),
fOption1(ob.fOption1),
fOption2(ob.fOption2) {
  // Copy constructor

}

//______________________________________________________________________
AliITSresponseSSD& AliITSresponseSSD::operator=(const AliITSresponseSSD& ob){
  // Assignment operator
  this->~AliITSresponseSSD();
  new(this) AliITSresponseSSD(ob);
  return *this;
}

