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
#include <TRandom3.h>

#include "AliNanoAODArrayMaker.h"

ClassImp(AliNanoAODArrayMaker)

//________________________________________________________________________
AliNanoAODArrayMaker::AliNanoAODArrayMaker(const char *name) 
  : AliAnalysisTaskSE(name), fOutputArrayName(), fOutputArray(0), fOutputArrayPythiaName(), fPythiaArray(0), fTrackEffPythia(1.0), fTrackEffData(1.0), fRandom(), fOutputList(0x0)
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

  fRandom = new TRandom3(0);

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
  if(fIsFirstLoop){
    InputEvent()->AddObject(fOutputArray);
    InputEvent()->AddObject(fPythiaArray);
    InputEvent()->AddObject(fDataArray);
    fIsFirstLoop = false;
  }

  //find NanoAOD particle array
  TClonesArray *particleArray = static_cast<TClonesArray*> (InputEvent()->FindListObject("Nanotracks"));
  Int_t nTracks = particleArray->GetEntries();

  Int_t accTracks = 0;
  Int_t accTracksPythia = 0;
  Int_t accTracksData = 0;

  //get custom NanoAOD variables which had to be defined in the nanoAOD generation
  Int_t indexHybGlob = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstIsGlobalHybrid");
  Int_t indexIsPyth = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstIsPythiaTrack");

  AliAODTrack* newTrack = new AliAODTrack();

  //Delte array members from the previous event
  fPythiaArray->Delete();
  fOutputArray->Delete();
  fDataArray->Delete();
  
  //loop over particles in the event and add them to the correct arrays
  for(Int_t iPart=0; iPart<nTracks; iPart++){
   AliNanoAODTrack *nanoTrack = (AliNanoAODTrack*) particleArray->At(iPart);
    
    if (nanoTrack->GetVar(indexIsPyth)==1){
      // Discard tracks due to lowered tracking efficiency
      if (fTrackEffPythia < 1.0 && fTrackEffPythia < fRandom->Rndm())
        continue;
      GetAODTrack(newTrack, nanoTrack);
      new ((*fPythiaArray)[accTracksPythia]) AliAODTrack(*newTrack);
      accTracksPythia++;
    }else{
      // Discard tracks due to lowered tracking efficiency
      if (fTrackEffData < 1.0 && fTrackEffData < fRandom->Rndm())
        continue;
      GetAODTrack(newTrack, nanoTrack,indexHybGlob);
      new ((*fDataArray)[accTracksData]) AliAODTrack(*newTrack);
      accTracksData++;
    }

    new ((*fOutputArray)[accTracks]) AliAODTrack(*newTrack);
    accTracks++;
  }
  
  delete newTrack;
}      

//________________________________________________________________________
void AliNanoAODArrayMaker::Terminate(Option_t *) 
{
  // Called once at the end of the query


}
//_____________________________________________________________________________________________________
void AliNanoAODArrayMaker::GetAODTrack(AliAODTrack* newTrack, AliNanoAODTrack* track, Int_t index)
{
  //create AOD track from NanoAOD track with the availbable information
  newTrack->SetPt(track->Pt());
  newTrack->SetTheta(2.*atan(exp(-track->Eta()))); // the same as in AliAnalysisTaskParticleRandomizer
  newTrack->SetPhi(track->Phi());
  newTrack->SetCharge(track->Charge());
  newTrack->SetLabel(track->GetLabel());
  if (index==-1 || track->GetVar(index) == 1) newTrack->SetIsHybridGlobalConstrainedGlobal();
  newTrack->SetFilterMap(track->GetFilterMap());

}

