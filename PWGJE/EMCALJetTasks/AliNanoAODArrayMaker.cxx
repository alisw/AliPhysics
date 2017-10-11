///
/// task to create arrays from NanoAODs which are used in the analysis
///
/// Author: M.Zimmermann
///

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliAODEvent.h"
#include "AliNanoAODTrack.h"
#include "AliAODTrack.h"

#include <TClonesArray.h>

#include "AliNanoAODArrayMaker.h"

ClassImp(AliNanoAODArrayMaker)

//________________________________________________________________________
AliNanoAODArrayMaker::AliNanoAODArrayMaker(const char *name) 
  : AliAnalysisTaskSE(name), fOutputArrayName(), fOutputArray(0), fOutputArrayPythiaName(), fPythiaArray(0), fOutputList(0x0)
{
  // Constructor

  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliNanoAODArrayMaker::UserCreateOutputObjects()
{
  // Create arrays and set names
  // Called once

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE); 

  fOutputArray = new TClonesArray("AliAODTrack");
  fOutputArray->SetName(fOutputArrayName.Data());
  Printf("%s \n", fOutputArrayName.Data());

  fDataArray = new TClonesArray("AliAODTrack");
  fDataArray->SetName(fOutputArrayDataName.Data());
  Printf("%s \n", fOutputArrayDataName.Data());

  fPythiaArray = new TClonesArray("AliAODTrack");
  fPythiaArray->SetName(fOutputArrayPythiaName.Data());
  Printf("%s \n", fOutputArrayPythiaName.Data());

  PostData(1,fOutputList);

}
//________________________________________________________________________
void AliNanoAODArrayMaker::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  InputEvent()->AddObject(fOutputArray);
  InputEvent()->AddObject(fPythiaArray);
  InputEvent()->AddObject(fDataArray);

  //find NanoAOD particle array
  TClonesArray *particleArray = static_cast<TClonesArray*> (InputEvent()->FindListObject("Nanotracks"));
  Int_t nTracks = particleArray->GetEntries();

  Int_t accTracks = 0;
  Int_t accTracksPythia = 0;
  Int_t accTracksData = 0;

  //get custom NanoAOD variables which had to be defined in the nanoAOD generation
  Int_t indexHybGlob = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstIsGlobalHybrid");
  Int_t indexIsPyth = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstIsPythiaTrack");

  fOutputArray->Clear("C");   
  fPythiaArray->Clear("C"); 
  fDataArray->Clear("C"); 

  //loop over particles in the event and add them to the correct arrays
  for(Int_t iPart=0; iPart<nTracks; iPart++){
   AliNanoAODTrack *nanoTrack = (AliNanoAODTrack*) particleArray->At(iPart);
    
    if (nanoTrack->GetVar(indexIsPyth)==1){
        new ((*fPythiaArray)[accTracksPythia]) AliAODTrack(*(GetAODTrack(nanoTrack)));
        new ((*fOutputArray)[accTracks]) AliAODTrack(*(GetAODTrack(nanoTrack)));
        accTracksPythia++;
    }else{
        new ((*fDataArray)[accTracksData]) AliAODTrack(*(GetAODTrack(nanoTrack,indexHybGlob)));
        new ((*fOutputArray)[accTracks]) AliAODTrack(*(GetAODTrack(nanoTrack,indexHybGlob)));
        accTracksData++;
    }
    accTracks++;
  }
}      

//________________________________________________________________________
void AliNanoAODArrayMaker::Terminate(Option_t *) 
{
  // Called once at the end of the query


}
//_____________________________________________________________________________________________________
AliAODTrack* AliNanoAODArrayMaker::GetAODTrack(AliNanoAODTrack* track, Int_t index)
{
  //create AOD track from NanoAOD track with the availbable information
  AliAODTrack* newTrack = new AliAODTrack();
  newTrack->SetPt(track->Pt());
  newTrack->SetTheta(2.*atan(exp(-track->Eta()))); // the same as in AliAnalysisTaskParticleRandomizer
  newTrack->SetPhi(track->Phi());
  newTrack->SetCharge(track->Charge());
  newTrack->SetLabel(track->GetLabel());
  if (track->GetVar(index) == 1) newTrack->SetIsHybridGlobalConstrainedGlobal();
  newTrack->SetFilterMap(track->GetFilterMap());

  return newTrack;
}
//_____________________________________________________________________________________________________
AliAODTrack* AliNanoAODArrayMaker::GetAODTrack(AliNanoAODTrack* track){

  //create AOD track from NanoAOD track with the availbable information
  AliAODTrack* newTrack = new AliAODTrack();
  newTrack->SetPt(track->Pt());
  newTrack->SetTheta(2.*atan(exp(-track->Eta()))); // the same as in AliAnalysisTaskParticleRandomizer
  newTrack->SetPhi(track->Phi());
  newTrack->SetCharge(track->Charge());
  newTrack->SetLabel(track->GetLabel());

  // Hybrid tracks (compatible with LHC11h)
  UInt_t filterMap = BIT(8) | BIT(9);
  newTrack->SetIsHybridGlobalConstrainedGlobal();
  newTrack->SetFilterMap(track->GetFilterMap());
  return newTrack;
}
