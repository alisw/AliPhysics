#include "AliBaseEventCut.h"
//________________________________
///////////////////////////////////////////////////////////
//
// class AliBaseEventCut
//
//
//
//
///////////////////////////////////////////////////////////

#include <AliESDtrack.h>
#include <AliESD.h>
ClassImp(AliBaseEventCut)

AliBaseEventCut::AliBaseEventCut():
 fMin(0.0),
 fMax(0.0),
{
  
}
/**********************************************************/

Bool_t AliBaseEventCut::Pass(AliESD* esd) const
{
  if ( (GetValue() < fMin) || (GetValue() > fMax) ) return kTRUE;
  return kFALSE;
{
/**********************************************************/
