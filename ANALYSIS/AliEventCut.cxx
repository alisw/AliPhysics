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

#include "AliEventBaseCut.h"

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

Bool_t AliEventCut::Rejected(AliAOD* aod) const
{
  //returns kTRUE if rejected
  if (aod == 0x0)
   {
     Error("Pass","Pointer to AOD is NULL. Not passed the cut");
     return kFALSE;
   }
   
  TIter iter(&fBaseCuts);
  AliEventBaseCut* becut;
  while (( becut = (AliEventBaseCut*)iter() ))
   {
     if (becut->Rejected(aod)) return kTRUE;
   }
  return kFALSE;
}
/*********************************************************/
void AliEventCut::AddBasePartCut(AliEventBaseCut* ebcut)
{
//Adds a base cut
 if (ebcut == 0x0)
  {
    Error("AddBasePartCut","Pointer to base cut is NULL");
    return;
  }
 
 if (ebcut->GetProperty() != AliEventBaseCut::kNone)
  {
    if (FindCut(ebcut->GetProperty()))
     {
       Warning("AddBasePartCut","Cut with this property is already in the list of base cuts");
     }
  }  
  
 fBaseCuts.Add(ebcut->Clone());
 
}
/*********************************************************/

AliEventBaseCut* AliEventCut::FindCut(AliEventBaseCut::EEventCutProperty prop)
{
//Finds and returns pointer to the cut with given property
 Int_t n = fBaseCuts.GetEntries();
 for (Int_t i = 0; i<n; i++)
  {
    AliEventBaseCut* bcut = (AliEventBaseCut*)fBaseCuts.At(i);
    if (bcut->GetProperty() == prop)
       return bcut; //we found the cut we were searching for
  }

 return 0x0; //we did not found this cut

}
/*********************************************************/

void AliEventCut::SetNChargedRange(Int_t min,Int_t max,Double_t etamin,Double_t etamax)
{
 //Sets renge of number of charged particles
  AliNChargedCut* cut = dynamic_cast<AliNChargedCut*>(FindCut(AliEventBaseCut::kNChargedCut));
  if(cut) 
   { 
     cut->SetRange(min,max);
     cut->SetEtaRange(etamin,etamax);
   }  
  else fBaseCuts.Add(new AliNChargedCut(min,max,etamin,etamax));
}

/*********************************************************/
/*********************************************************/
/*********************************************************/

ClassImp(AliEventEmptyCut)
