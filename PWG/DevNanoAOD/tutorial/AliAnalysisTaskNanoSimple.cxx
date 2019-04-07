#include <TList.h>
#include <TMath.h>

#include "AliAnalysisTaskNanoSimple.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"

#include "AliNanoAODHeader.h"
#include "AliNanoAODTrack.h"

#include "AliLog.h"
#include "AliVParticle.h"
#include "AliVTrack.h"

ClassImp(AliAnalysisTaskNanoSimple)

//____________________________________________________________________
AliAnalysisTaskNanoSimple:: AliAnalysisTaskNanoSimple(const char* name):
AliAnalysisTaskSE(name),
// general configuration
fListOfHistos(0x0)
{
  // Default constructor

  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

AliAnalysisTaskNanoSimple::~AliAnalysisTaskNanoSimple() 
{ 
  // destructor
  
  if (fListOfHistos) 
    delete fListOfHistos;
}

//____________________________________________________________________
void  AliAnalysisTaskNanoSimple::UserCreateOutputObjects()
{
  // Create the output container
  
  // Initialize output list of containers
  if (fListOfHistos != NULL){
    delete fListOfHistos;
    fListOfHistos = NULL;
  }
  if (!fListOfHistos){
    fListOfHistos = new TList();
    fListOfHistos->SetOwner(kTRUE); 
  }

  // TODO add output histograms to this list. E.g.
  //fListOfHistos->Add(new TH2F("multVsMPI", ";n_{MPI};nch_alice;events", 20, -0.5, 19.5, 301, -0.5, 300.5));
  
  PostData(1, fListOfHistos);
}

//____________________________________________________________________
void  AliAnalysisTaskNanoSimple::UserExec(Option_t */*option*/)
{
  // exec (per event)

  if (!fInputEvent)
    return;
  
  const AliVVertex* vertex = fInputEvent->GetPrimaryVertex();
  if (vertex->GetNContributors() < 1)
    return;
  
  Printf("Vertex = %f", vertex->GetZ());
  
  // for custom variables, cast to nano AOD header
  AliNanoAODHeader* nanoHeader = dynamic_cast<AliNanoAODHeader*>(fInputEvent->GetHeader());
  const Int_t kRefMult = nanoHeader->GetVarIndex("MultSelection.RefMult08");
  if (kRefMult != -1)
    Printf("Ref Mult = %f", nanoHeader->GetVar(kRefMult));
  
  unsigned int nTracks = fInputEvent->GetNumberOfTracks();
  for (unsigned int i = 0; i < nTracks; i++) {
    AliVTrack* track = (AliVTrack*) fInputEvent->GetTrack(i);
    Printf("pt = %f   ITS cluster = %d %d %d %d %d %d", track->Pt(), track->HasPointOnITSLayer(0), track->HasPointOnITSLayer(1), 
           track->HasPointOnITSLayer(2), track->HasPointOnITSLayer(3), track->HasPointOnITSLayer(4), track->HasPointOnITSLayer(5));
    //Printf("TOF BC = %d", track->GetTOFBunchCrossing());
    
    // for custom variables, cast to nano AOD track
    AliNanoAODTrack* nanoTrack = dynamic_cast<AliNanoAODTrack*>(track);
    // NOTE Access to custom variables. Use static here for caching of index
    static Bool_t bPIDAvailable = AliNanoAODTrack::InitPIDIndex();
    static const Int_t kcstNSigmaTPCPr  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, AliPID::kProton);
    static const Int_t kcstNSigmaTOFPr  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, AliPID::kProton);
    if (nanoTrack && bPIDAvailable)
      Printf("TPC_sigma_proton = %f  hasTOF = %d  TOF_sigma_proton = %f", nanoTrack->GetVar(kcstNSigmaTPCPr), nanoTrack->HasTOFPID(), nanoTrack->GetVar(kcstNSigmaTOFPr));
  }
  
  // V0 access - as usual
  AliAODEvent* aod = dynamic_cast<AliAODEvent*> (fInputEvent);
  if (aod->GetV0s()) {
    for (int i = 0; i < aod->GetNumberOfV0s(); i++)
      Printf("V0 %d: dca = %f", i, aod->GetV0(i)->DcaV0ToPrimVertex());
  }

  // cascade access - as usual
  
  // TODO in current AliRoot tag (v5-09-46), GetCascades() produces a SEGV if not filled
  //if (aod->GetCascades()) {
  //  for (int i = 0; i < aod->GetNumberOfCascades(); i++)
  //    Printf("Cascade %d: xi mass = %f", i, aod->GetCascade(i)->MassXi());
  //}
}
