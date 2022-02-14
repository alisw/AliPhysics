#include <TList.h>
#include <TMath.h>

#include "AliAnalysisTaskNanoValidator.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"

#include "AliNanoAODHeader.h"
#include "AliNanoAODTrack.h"

#include "AliAODConversionPhoton.h"

#include "AliLog.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliMultSelection.h"

ClassImp(AliAnalysisTaskNanoValidator)

//____________________________________________________________________
AliAnalysisTaskNanoValidator:: AliAnalysisTaskNanoValidator(const char* name):
AliAnalysisTaskSE(name),
// general configuration
fListOfHistos(0x0)
{
  // Default constructor

  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

AliAnalysisTaskNanoValidator::~AliAnalysisTaskNanoValidator() 
{ 
  // destructor
  
  if (fListOfHistos) 
    delete fListOfHistos;
}

//____________________________________________________________________
void  AliAnalysisTaskNanoValidator::UserCreateOutputObjects()
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

  PostData(1, fListOfHistos);
}

//____________________________________________________________________
void  AliAnalysisTaskNanoValidator::UserExec(Option_t */*option*/)
{
  // exec (per event)
  // print all fields in ESD, AOD and NanoAOD to allow comparison

  if (!fInputEvent)
    return;
  
  // header validation
  
  Double_t CentrV0M = -1;
  Double_t CentrCL0 = -1;
  Double_t CentrCL1 = -1;
  Double_t CentrTRK = -1;
  UInt_t OfflineTrigger = 0;
  Double_t RefMult08 = -1;
  
  // same method
  UInt_t BunchCrossNumber = fInputEvent->GetHeader()->GetBunchCrossNumber();
  UInt_t OrbitNumber = fInputEvent->GetHeader()->GetOrbitNumber();
  UInt_t PeriodNumber = fInputEvent->GetHeader()->GetPeriodNumber();
  Int_t RunNumber = fInputEvent->GetRunNumber();
  Double_t MagField = fInputEvent->GetMagneticField();
  Int_t NumberOfESDTracks = fInputEvent->GetNumberOfESDTracks();
  Float_t T0Spread[4] = { -1 };
  for (int i=0; i<4; i++)
    T0Spread[i] = fInputEvent->GetT0spread(i);
  
  AliESDHeader* esdHeader = dynamic_cast<AliESDHeader*>(fInputEvent->GetHeader());
  AliAODHeader* aodHeader = dynamic_cast<AliAODHeader*>(fInputEvent->GetHeader());
  AliNanoAODHeader* nanoHeader = dynamic_cast<AliNanoAODHeader*>(fInputEvent->GetHeader());
  if (!esdHeader && !aodHeader && !nanoHeader)
    AliFatal("No valid header found");
  
  // centrality
  if (nanoHeader) {
    CentrV0M = nanoHeader->GetCentr("V0M");
    CentrCL0 = nanoHeader->GetCentr("CL0");
    CentrCL1 = nanoHeader->GetCentr("CL1");
    CentrTRK = nanoHeader->GetCentr("TRK");

    static const Int_t kRefMult = nanoHeader->GetVarIndex("MultSelection.RefMult08");
    RefMult08 = nanoHeader->GetVar(kRefMult);
  } else {
    AliMultSelection* MultSelection = dynamic_cast<AliMultSelection*> (fInputEvent->FindListObject("MultSelection"));
    if(!MultSelection)
      AliFatal("No multiplicity selection");
    
    CentrV0M = MultSelection->GetMultiplicityPercentile("V0M");
    CentrCL1 = MultSelection->GetMultiplicityPercentile("CL1");
    CentrCL0 = MultSelection->GetMultiplicityPercentile("CL0");
    CentrTRK = MultSelection->GetMultiplicityPercentile("TRK");  
    RefMult08 = MultSelection->GetMultiplicityPercentile("RefMult08");  
  }
  
  if (esdHeader) {
    OfflineTrigger = fInputHandler->IsEventSelected();
  } else if (aodHeader) {
    OfflineTrigger = aodHeader->GetOfflineTrigger();
  } else if (nanoHeader) {
    OfflineTrigger = nanoHeader->GetOfflineTrigger();
  }

  Printf("Base: %d \t %d \t %d \t %d \t %f \t %d \t %d", BunchCrossNumber, OrbitNumber, PeriodNumber, RunNumber, MagField, NumberOfESDTracks, OfflineTrigger);
  Printf("T0: %f \t %f \t %f \t %f", T0Spread[0], T0Spread[1], T0Spread[2], T0Spread[3]);
  Printf("Centrality: %f \t %f \t %f \t %f \t %f", CentrV0M, CentrCL0, CentrCL1, CentrTRK, RefMult08);
  Printf(" ");
  
  // track validation
  
  unsigned int nTracks = fInputEvent->GetNumberOfTracks();
  for (unsigned int i = 0; i < nTracks; i++) {
    AliVTrack* track = (AliVTrack*) fInputEvent->GetTrack(i);
    Double_t pos[3];
    track->GetXYZ(pos);
    
    Printf("1) %f %f %f %d %d %d %d %d %d %lld", track->Pt(), track->Phi(), track->Theta(), track->HasPointOnITSLayer(0), track->HasPointOnITSLayer(1), 
           track->HasPointOnITSLayer(2), track->HasPointOnITSLayer(3), track->HasPointOnITSLayer(4), track->HasPointOnITSLayer(5),
           track->GetStatus()
          );
    Printf("2) %f %f %f %d %d %d", pos[0], pos[1], pos[2], track->GetTPCncls(), track->GetID(), track->GetTPCNclsF()
          );
    Printf("3) %f %f %f %f %f %f %d %f %f", track->GetTrackPhiOnEMCal(), track->GetTrackEtaOnEMCal(), track->GetTrackPtOnEMCal(),
           track->GetITSsignal(), track->GetTPCsignal(), track->GetTPCsignalTunedOnData(), track->GetTPCsignalN(), track->GetTPCmomentum(), track->GetTPCTgl()
          );
    Printf("4) %f %f %f %f %f %f %f %d %d", track->GetTOFsignal(), track->GetIntegratedLength(), track->GetTOFsignalTunedOnData(), track->GetHMPIDsignal(), 
           track->GetHMPIDoccupancy(), track->GetTRDsignal(), track->GetTRDchi2(), track->GetNumberOfTRDslices(), track->GetTOFBunchCrossing()
          );
 
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*> (track);
    if (aodTrack) {
      Printf("5) %f %f %f %f %f %f %f %f %f %d %d %d", aodTrack->PxAtDCA(), aodTrack->PyAtDCA(), aodTrack->PzAtDCA(), aodTrack->Chi2perNDF(),
             aodTrack->XAtDCA(), aodTrack->YAtDCA(), aodTrack->ZAtDCA(), aodTrack->DCA(), aodTrack->GetRAtAbsorberEnd(), 
             aodTrack->GetFilterMap(), aodTrack->GetTPCNCrossedRows(), aodTrack->GetTPCnclsS()
      );
    }
    AliNanoAODTrack* nanoTrack = dynamic_cast<AliNanoAODTrack*> (track);
    if (nanoTrack) {
      Printf("5) %f %f %f %f %f %f %f %f %f %d %d %d", nanoTrack->PxAtDCA(), nanoTrack->PyAtDCA(), nanoTrack->PzAtDCA(), nanoTrack->Chi2perNDF(),
             nanoTrack->XAtDCA(), nanoTrack->YAtDCA(), nanoTrack->ZAtDCA(), nanoTrack->DCA(), nanoTrack->GetRAtAbsorberEnd(), 
             nanoTrack->GetFilterMap(), nanoTrack->GetTPCNCrossedRows(), nanoTrack->GetTPCnclsS()
      );
    }
    Printf(" ");
    
    // TODO missing: covmat
  }
}

 
