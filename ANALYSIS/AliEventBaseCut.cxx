#include "AliEventBaseCut.h"
//________________________________
///////////////////////////////////////////////////////////
//
// class AliEventBaseCut
//
//
//
//
///////////////////////////////////////////////////////////

#include <AliAOD.h>
ClassImp(AliEventBaseCut)

AliEventBaseCut::AliEventBaseCut():
 fMin(0.0),
 fMax(0.0),
 fProperty(kNone)
{
//ctor  
}
/**********************************************************/

AliEventBaseCut::AliEventBaseCut(Double_t min, Double_t max, EEventCutProperty prop):
 fMin(min),
 fMax(max),
 fProperty(prop)
{
 //ctor
}
/**********************************************************/

Bool_t AliEventBaseCut::Rejected(AliAOD* aod) const
{
  if ( (GetValue(aod) < fMin) || (GetValue(aod) > fMax) ) return kTRUE;
  return kFALSE;
}
/**********************************************************/
/**********************************************************/
/**********************************************************/
ClassImp(AliPrimVertexXCut)

Double_t AliPrimVertexXCut::GetValue(AliAOD* aod) const
{
 //returns x coordinate of the primary vertex
  Double_t x = 0, y = 0, z = 0;
  if (aod) aod->GetPrimaryVertex(x,y,z);
  return x;
}
/**********************************************************/
/**********************************************************/
/**********************************************************/
ClassImp(AliPrimVertexYCut)

Double_t AliPrimVertexYCut::GetValue(AliAOD* aod) const
{
 //returns x coordinate of the primary vertex
  Double_t x = 0, y = 0, z = 0;
  if (aod) aod->GetPrimaryVertex(x,y,z);
  return y;
}
/**********************************************************/
/**********************************************************/
/**********************************************************/
ClassImp(AliPrimVertexZCut)

Double_t AliPrimVertexZCut::GetValue(AliAOD* aod) const
{
 //returns x coordinate of the primary vertex
  Double_t x = 0, y = 0, z = 0;
  if (aod) aod->GetPrimaryVertex(x,y,z);
  return z;
}

/**********************************************************/
/**********************************************************/
/**********************************************************/

ClassImp(AliNChargedCut)

Double_t AliNChargedCut::GetValue(AliAOD* aod) const
{
  //returns number of charged particles
  if (aod) 
   {
     return aod->GetNumberOfCharged(fEtaMin,fEtaMax);
   }  
  Error("GetValue","Pointer to AOD is NULL");
  return 0.0;
}
