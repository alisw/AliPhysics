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

Bool_t AliEventCut::Pass(AliAOD* aod) const
{
  //returns kTRUE if rejected
  if (aod == 0x0)
   {
     Error("Pass","Pointer to AOD is NULL. Not passed the cut");
     return kFALSE;
   }
   
  TIter iter(fBaseCuts);
  AliBaseEventCut* becut;
  while (( becut = (AliBaseEventCut*)iter() ))
   {
     if (becut->Pass(aod)) return kTRUE;
   }
  return kFALSE;
}
