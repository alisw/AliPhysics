/**************************************************************************
 * Author: Panos Christakoglou.                                           *
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
//                   AliLHCTagCuts class
//   This is the class to deal with the LHC tag level cuts
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

class AliLog;

#include "AliLHCTag.h"
#include "AliLHCTagCuts.h"

ClassImp(AliLHCTagCuts)


//___________________________________________________________________________
AliLHCTagCuts::AliLHCTagCuts() :
  TObject(),
  fLHCState(0),
  fLHCStateFlag(kFALSE),
  fLHCLuminosityMin(0),
  fLHCLuminosityMax(0),
  fLHCLuminosityFlag(kFALSE)
{
  //Default constructor which calls the Reset method.
  Reset();
}

//___________________________________________________________________________
AliLHCTagCuts::~AliLHCTagCuts() {  
  //Defaut destructor.
}

//___________________________________________________________________________
void AliLHCTagCuts::Reset() {
  //Sets dummy values to every private member.
  fLHCState = "init";
  fLHCStateFlag = kFALSE;
  fLHCLuminosityMin = -1.0;
  fLHCLuminosityMax = 0.0;
  fLHCLuminosityFlag = kFALSE;
}

//___________________________________________________________________________
Bool_t AliLHCTagCuts::IsAccepted(AliLHCTag *lhcTag) const {
  //Returns true if the event is accepted otherwise false.
  if(fLHCStateFlag)
    if((lhcTag->GetLHCState() != fLHCState))
      return kFALSE;
  if(fLHCLuminosityFlag)
    if((lhcTag->GetLuminosity() < fLHCLuminosityMin)||(lhcTag->GetLuminosity() > fLHCLuminosityMax))
      return kFALSE;
 
  return kTRUE;
}
