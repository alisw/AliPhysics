#include "AliEventCut.h"
//________________________________
///////////////////////////////////////////////////////////
//
// class AliRunAnalysis
//
//
//
//
///////////////////////////////////////////////////////////

#include <TObjArray.h>
//#include <TIter.h>

#include "AliBaseEventCut.h"

AliEventCut::AliEventCut():
 fBaseCuts(0x0)
{
//costructor

}
/*********************************************************/

AliEventCut::~AliEventCut()
{
//costructor
 delete fBaseCuts;
}

/*********************************************************/

Bool_t AliEventCut::Pass(AliESD* esd) const
{
  //returns kTRUE if rejected
  TIter iter(fBaseCuts);
  AliBaseEventCut* becut;
  while (( becut = (AliBaseEventCut*)iter() ))
   {
     if (becut->Pass(esd)) return kTRUE;
   }
  return kFALSE;
}
