/* $Id$ */


// Simple task to test the collision normalization. Can be called for
// both MC and data and shows how to fill the collision normalization
// class
//
// Author: Michele Floris
//         CERN


#include "AliCollisionNormalizationTask.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliHeader.h>

#include "AliCollisionNormalization.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

//#include "AliBackgroundSelection.h"
#include "AliMultiplicity.h"
#include "AliMCEvent.h"

ClassImp(AliCollisionNormalizationTask)

AliCollisionNormalizationTask::AliCollisionNormalizationTask() :
  AliAnalysisTaskSE("AliCollisionNormalizationTask"),
  fOutput(0),
  fIsMC(0),
  fCollisionNormalization(0)
{
  //
  // Default event handler
  //

  // Define input and output slots here
  DefineOutput(1, TList::Class());
  
}

AliCollisionNormalizationTask::AliCollisionNormalizationTask(const char* name) :
  AliAnalysisTaskSE(name),
  fOutput(0),
  fIsMC(0),
  fCollisionNormalization(new AliCollisionNormalization())
{
  //
  // Constructor. Initialization of pointers
  //

  // Define input and output slots here
  DefineOutput(1, TList::Class());
  
  //  AliLog::SetClassDebugLevel("AliCollisionNormalizationTask", AliLog::kWarning);
}

AliCollisionNormalizationTask::~AliCollisionNormalizationTask()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

void AliCollisionNormalizationTask::UserCreateOutputObjects()
{
  // create result objects and add to output list

  Printf("AliCollisionNormalizationTask::CreateOutputObjects");

  fOutput = new TList;
  fOutput->SetOwner();
  
  if (!fCollisionNormalization)
    fCollisionNormalization = new AliCollisionNormalization;
  
  fOutput->Add(fCollisionNormalization);
//   fOutput->Add(fCollisionNormalization->GetVzCorrZeroBin ());
//   fOutput->Add(fCollisionNormalization->GetVzMCGen       ());
//   fOutput->Add(fCollisionNormalization->GetVzMCRec       ());
//   fOutput->Add(fCollisionNormalization->GetVzMCTrg       ());
//   fOutput->Add(fCollisionNormalization->GetVzData        ());
//   fOutput->Add(fCollisionNormalization->GetNEvents       ());
//   fOutput->Add(fCollisionNormalization->GetStatBin0      ());

}

void AliCollisionNormalizationTask::UserExec(Option_t*)
{
  // process the event


  PostData(1, fOutput);

  // Get the ESD
  AliESDEvent * aESD = dynamic_cast<AliESDEvent*>(fInputEvent);
  if (strcmp(aESD->ClassName(),"AliESDEvent")) {
    AliFatal("Not processing ESDs");
  }
  // Get MC event, if needed
  AliMCEvent* mcEvent = fIsMC ? MCEvent() : 0;
  if (!mcEvent && fIsMC){
    AliFatal("Running on MC but no MC handler available");
  }

  // Physics selection. At least in the case of MC we cannot use
  // yourTask->SelectCollisionCandidates();, because we also need to
  // fill the "generated" histogram 
  // NB never call IsEventSelected more than once per event
  // (statistics histogram would be altered)

  Bool_t isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  // Get the Multiplicity cut
  const AliMultiplicity* mult = aESD->GetMultiplicity();
  if (!mult){
    AliError("Can't get mult object");
    return;
  }

  // Check if the event is "reconstructed" according to some quality
  // cuts. Tipically, we mean "it has a good-eonugh vertex"

  Int_t ntracklet = mult->GetNumberOfTracklets();
  const AliESDVertex * vtxESD = aESD->GetPrimaryVertexSPD();
  if (IsEventInBinZero()) {
    ntracklet = 0;
    vtxESD    = 0;
  }
  
  if (ntracklet > 0 && !vtxESD) {
    AliError("No vertex but reconstructed tracklets?");
  }

  // assign vz. For MC we use generated vz
  Float_t vz = 0;
  if (!fIsMC) vz = vtxESD ? vtxESD->GetZ() : 0; // FIXME : is zv used anywhere in Gen?
  else        vz = mcEvent->GetPrimaryVertex()->GetZ();

  if (fIsMC) {
    // Monte Carlo:  we fill 3 histos
    if (!isSelected || !vtxESD) ntracklet = 0; //If the event does not pass the physics selection or is not rec, it goes in the bin0
    fCollisionNormalization->FillVzMCGen(vz, ntracklet, mcEvent);      
    // If triggered == passing the physics selection
    if (isSelected) {
      fCollisionNormalization->FillVzMCTrg(vz, ntracklet, mcEvent);
      // If reconstructer == good enough vertex
      if (vtxESD) fCollisionNormalization->FillVzMCRec(vz, ntracklet, mcEvent);    
    }
  } else {
    if (isSelected) {
      // Passing the trigger
      fCollisionNormalization->FillVzData(vz,ntracklet);
    }
  }

}

void AliCollisionNormalizationTask::Terminate(Option_t *)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput)
    Printf("ERROR: fOutput not available");


}

Bool_t AliCollisionNormalizationTask::IsEventInBinZero() {

  // Returns true if an event is to be assigned to the zero bin
  //
  // You should have your own version of this method in the class: in
  // general, the definition of "reconstructed" event is subjective.

  Bool_t isZeroBin = kTRUE;
  const AliESDEvent* esd= dynamic_cast<AliESDEvent*>(fInputEvent);
  const AliMultiplicity* mult = esd->GetMultiplicity();
  if (!mult){
    Printf("AliAnalysisTaskBGvsTime::IsBinZero: Can't get mult object");
    return kFALSE;
  }
  Int_t ntracklet = mult->GetNumberOfTracklets();
  const AliESDVertex * vtxESD = esd->GetPrimaryVertexSPD();
  if(vtxESD) {
    // If there is a vertex from vertexer z with delta phi > 0.02 we
    // don't consider it rec (we keep the event in bin0). If quality
    // is good eneough we check the number of tracklets
    // if the vertex is more than 15 cm away, this is autamatically bin0
    if( TMath::Abs(vtxESD->GetZ()) <= 15 ) {
      if (vtxESD->IsFromVertexerZ()) {
	if (vtxESD->GetDispersion()<=0.02 ) {
	  if(ntracklet>0) isZeroBin = kFALSE;
	}
      } else if(ntracklet>0) isZeroBin = kFALSE; // if the event is not from Vz we chek the n of tracklets
    } 
  }
  return isZeroBin;

}
