/*
 * **********************************************************
 * Virtual class for processing trees of AliReducedEventInfo
 * Authors: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no
 *                Jacobus Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
 * Creation date: 2015/10/01
 *********************************************************
 */

#include "AliReducedAnalysisTaskSE.h"

ClassImp(AliReducedAnalysisTaskSE);


//___________________________________________________________________________
AliReducedAnalysisTaskSE::AliReducedAnalysisTaskSE() :
  TObject(),
//  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fName(""),
  fTitle(""),
  fEvent(0x0)
{
  //
  // default constructor
  //
   //fHistosManager->SetUseDefaultVariableNames(kTRUE);
   //fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,AliReducedVarManager::fgVariableUnits);
  for(Int_t i=0; i<AliReducedVarManager::kNVars; ++i)
    fValues[i] = 0.0;
}


//___________________________________________________________________________
AliReducedAnalysisTaskSE::AliReducedAnalysisTaskSE(const Char_t* name, const Char_t* title) :
  TObject(),
  //fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fName(name),
  fTitle(title),
  fEvent(0x0)
{
  //
  // named constructor
  //
   //fHistosManager->SetUseDefaultVariableNames(kTRUE);
   //fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,AliReducedVarManager::fgVariableUnits);
  for(Int_t i=0; i<AliReducedVarManager::kNVars; ++i)
    fValues[i] = 0.0;
}


//___________________________________________________________________________
AliReducedAnalysisTaskSE::~AliReducedAnalysisTaskSE() 
{
  //
  // destructor
  //
  /*if(fEventCuts) {fEventCuts->Clear("C"); delete fEventCuts;}
  if(fTrackCuts) {fTrackCuts->Clear("C"); delete fTrackCuts;}
  if(fPairCuts) {fPairCuts->Clear("C"); delete fPairCuts;}*/
  //if(fHistosManager) delete fHistosManager;
}

//___________________________________________________________________________
void AliReducedAnalysisTaskSE::Init() {
   //
   // initialization (typically called in AliAnalysisTask::UserCreateOutputObjects())
   //
}

//___________________________________________________________________________
void AliReducedAnalysisTaskSE::Process() {
   //
   // process a given event (typically called in AliAnalysisTask::UserExec())
   //
}

//___________________________________________________________________________
void AliReducedAnalysisTaskSE::Finish() {
   //
   // finish, to be executed after all events were processed
   //
}
