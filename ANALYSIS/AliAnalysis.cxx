#include "AliAnalysis.h"
//________________________________
///////////////////////////////////////////////////////////
//
// class AliAnalysis
//
// Base class for analysis.
// Each inheriting calss must define 3 methods:
//   - Init() : that is called before event processing
//   - ProcessEvent(AliESD*,AliStack*)
//   - 
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////

#include "AliEventCut.h"

ClassImp(AliAnalysis)

AliAnalysis::AliAnalysis():
 fEventCut(0x0),
 fCutOnSim(kTRUE),
 fCutOnRec(kTRUE)
{
 //ctor
}
/*********************************************************/

AliAnalysis::AliAnalysis(const char* name,const char* title):
 TTask(name,title),
 fEventCut(0x0),
 fCutOnSim(kTRUE),
 fCutOnRec(kTRUE)
{
 //ctor
}
/*********************************************************/

AliAnalysis::~AliAnalysis()
{
 //dtor
 delete fEventCut;
}
/*********************************************************/

void AliAnalysis::SetEventCut(AliEventCut* evcut)
{
//Sets event -  makes a private copy
  delete fEventCut;
  if (evcut) fEventCut = (AliEventCut*)evcut->Clone();
  else fEventCut = 0x0;
}
/*********************************************************/

Bool_t AliAnalysis::Pass(AliAOD* recevent, AliAOD* simevent)
{
  //checks the event cut
  if (fEventCut == 0x0) return kFALSE;
  
  if (fCutOnRec)
    if (fEventCut->Pass(recevent)) return kTRUE;
    
  if (fCutOnSim)
    if (fEventCut->Pass(simevent)) return kTRUE;
  
  return kFALSE;
}
