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

ClassImp(AliEventCut)


AliEventCut::AliEventCut():
 fBaseCuts(10)
{
//costructor

}
/*********************************************************/
AliEventCut::AliEventCut(const AliEventCut& in):
 TObject(in),
 fBaseCuts(in.fBaseCuts)
{
  //cpy ctor
  fBaseCuts.SetOwner(kTRUE);
}
/*********************************************************/

AliEventCut::~AliEventCut()
{
//costructor
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
   
  TIter iter(&fBaseCuts);
  AliBaseEventCut* becut;
  while (( becut = (AliBaseEventCut*)iter() ))
   {
     if (becut->Pass(aod)) return kTRUE;
   }
  return kFALSE;
}

/*********************************************************/
/*********************************************************/
/*********************************************************/

ClassImp(AliEmptyEventCut)
