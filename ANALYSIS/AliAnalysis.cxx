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
#include "AliAODPairCut.h"

ClassImp(AliAnalysis)

AliAnalysis::AliAnalysis():
 fEventCut(0x0),
 fCutOnSim(kTRUE),
 fCutOnRec(kTRUE),
 fPairCut(new AliAODEmptyPairCut()),//empty cut - accepts all particles
 fkPass(&AliAnalysis::PassPartAndTrack), //by default perform cut on both track and particle pair 
 fkPass1(&AliAnalysis::PassPartAndTrack1), //used onluy by ProcessTracksAndParticles
 fkPass2(&AliAnalysis::PassPartAndTrack2),
 fkPassPairProp(&AliAnalysis::PassPairPropPartAndTrack)
{
 //ctor
}
/*********************************************************/

AliAnalysis::AliAnalysis(const char* name,const char* title):
 TTask(name,title),
 fEventCut(0x0),
 fCutOnSim(kTRUE),
 fCutOnRec(kTRUE),
 fPairCut(new AliAODEmptyPairCut()),//empty cut - accepts all particles
 fkPass(&AliAnalysis::PassPartAndTrack), //by default perform cut on both track and particle pair 
 fkPass1(&AliAnalysis::PassPartAndTrack1), //used onluy by ProcessTracksAndParticles
 fkPass2(&AliAnalysis::PassPartAndTrack2),
 fkPassPairProp(&AliAnalysis::PassPairPropPartAndTrack)
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
/*************************************************************************************/ 

void AliAnalysis::SetPairCut(AliAODPairCut* cut)
{
//Sets new Pair Cut. Old one is deleted
//Note that it is created new object instead of simple pointer set
//I do not want to have pointer
//to object created somewhere else
//because in that case I could not believe that
//it would always exist (sb could delete it)
//so we have always own copy

 if(!cut)
   {
     Error("AliHBTFunction::SetPairCut","argument is NULL");
     return;
   }
 delete fPairCut;
 fPairCut = (AliAODPairCut*)cut->Clone();

}
/******************************************************************/

void AliAnalysis::SetCutsOnSim()
{
 // -- aplies only to Process("TracksAndParticles")
 // (ProcessTracksAndParticles and ProcessTracksAndParticlesNonIdentAnal)
 // Only particles properties are checkes against cuts
  
  fCutOnRec = kFALSE;
  fCutOnSim = kTRUE;
  
  fkPass = &AliAnalysis::PassPart;
  fkPass1 = &AliAnalysis::PassPart1;
  fkPass2 = &AliAnalysis::PassPart2;
  fkPassPairProp = &AliAnalysis::PassPairPropPart;
  
}
/*************************************************************************************/ 

void AliAnalysis::SetCutsOnRec()
{
 // -- aplies only to Process("TracksAndParticles")
 // (ProcessTracksAndParticles and ProcessTracksAndParticlesNonIdentAnal)
 // Only tracks properties are checkes against cuts
  
  fCutOnRec = kTRUE;
  fCutOnSim = kFALSE;
 
  fkPass = &AliAnalysis::PassTrack;
  fkPass1 = &AliAnalysis::PassTrack1;
  fkPass2 = &AliAnalysis::PassTrack2;
  fkPassPairProp = &AliAnalysis::PassPairPropTrack;
 
}
/*************************************************************************************/ 

void AliAnalysis::SetCutsOnRecAndSim()
{
 // -- aplies only to Process("TracksAndParticles")
 // (ProcessTracksAndParticles and ProcessTracksAndParticlesNonIdentAnal)
 // Both, tracks and particles, properties are checked against cuts

  fCutOnRec = kTRUE;
  fCutOnSim = kTRUE;

  fkPass = &AliAnalysis::PassPartAndTrack;
  fkPass1 = &AliAnalysis::PassPartAndTrack1;
  fkPass2 = &AliAnalysis::PassPartAndTrack2;
  fkPassPairProp = &AliAnalysis::PassPairPropPartAndTrack;
}
