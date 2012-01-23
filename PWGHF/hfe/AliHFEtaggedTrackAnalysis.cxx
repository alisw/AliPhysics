/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Class AliHFEtaggedTrackAnalysis
// Analyses tracks with an apriori PID information (i.e. using the daugther
// tracks from well-identified decays of neutral charged particles). Tracks
// are processed in the Process function, where given tracks are filtered 
// via the track cuts and used for PID later. The plugin fills Correction
// Framework containers and additional PID QA containers
//
// Author:
//   Markus Fasel <M.Fasel@gsi.de>
//
#include "AliCFCutBase.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliExternalTrackParam.h"

#include "AliHFEcollection.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidTPC.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtaggedTrackAnalysis.h"
#include "AliHFEvarManager.h"

ClassImp(AliHFEtaggedTrackAnalysis)
//____________________________________________________________
AliHFEtaggedTrackAnalysis::AliHFEtaggedTrackAnalysis():
    TNamed()
  , fVarManager(NULL)
  , fContainer(NULL)
  , fPID(NULL)
  , fPIDqa(NULL)
  , fCuts(NULL)
  , fCFM(NULL)
  , fQAhistos(NULL)
  , fCentralityF(0)
  , fClean(kFALSE)
  , fMagneticField(0.0)
  , fVariablesTRD(kFALSE)
  , fIsPbPb(kFALSE)
{
  //
  // Dummy constructor
  //
}

//____________________________________________________________
AliHFEtaggedTrackAnalysis::AliHFEtaggedTrackAnalysis(const char *name):
    TNamed(name, "")
  , fVarManager(NULL)
  , fContainer(NULL)
  , fPID(NULL)
  , fPIDqa(NULL)
  , fCuts(NULL)
  , fCFM(NULL)
  , fQAhistos(NULL)
  , fCentralityF(0)
  , fClean(kFALSE)
  , fMagneticField(0.0)
  , fVariablesTRD(kFALSE)
  , fIsPbPb(kFALSE)
{
  //
  // Default constructor
  //
  fVarManager = new AliHFEvarManager("taggedTrackVarManager");
  //fVarManager->AddVariable("pt");
  //fVarManager->AddVariable("eta");
  //fVarManager->AddVariable("phi");
  //fVarManager->AddVariable("charge");
  //fVarManager->AddVariable("species");
  fPIDqa = new AliHFEpidQAmanager;
  fCFM = new AliCFManager;
  SetBit(kIsOwner, kTRUE);
}

//____________________________________________________________
AliHFEtaggedTrackAnalysis::AliHFEtaggedTrackAnalysis(const AliHFEtaggedTrackAnalysis &ref):
    TNamed(ref)
  , fVarManager(ref.fVarManager)
  , fContainer(NULL)
  , fPID(ref.fPID)
  , fPIDqa(ref.fPIDqa)
  , fCuts(ref.fCuts)
  , fCFM(ref.fCFM)
  , fQAhistos(ref.fQAhistos)
  , fCentralityF(ref.fCentralityF)
  , fClean(ref.fClean)
  , fMagneticField(ref.fMagneticField)
  , fVariablesTRD(ref.fVariablesTRD)
  , fIsPbPb(ref.fIsPbPb)
{
  //
  // Copy constructor
  //
  if(ref.fContainer){
    InitContainer();
  }
  SetBit(kIsOwner, kFALSE);
}

//____________________________________________________________
AliHFEtaggedTrackAnalysis &AliHFEtaggedTrackAnalysis::operator=(const AliHFEtaggedTrackAnalysis &ref){
  //
  // Assignment operator
  //
  TNamed::operator=(ref);
  if(&ref != this){
    fVarManager = ref.fVarManager;
    fPID = ref.fPID;
    fPIDqa = ref.fPIDqa;
    fCuts = ref.fCuts;
    fCFM = ref.fCFM;
    fQAhistos = ref.fQAhistos;
    fCentralityF = ref.fCentralityF;
    fClean = ref.fClean;
    fMagneticField = ref.fMagneticField;
    fVariablesTRD = ref.fVariablesTRD;
    fIsPbPb = ref.fIsPbPb;

    if(ref.fContainer) InitContainer();
   
    SetBit(kIsOwner, kFALSE); 
    SetBit(kIsOwnerCuts, kFALSE);
  }
  return *this;
}

//____________________________________________________________
AliHFEtaggedTrackAnalysis::~AliHFEtaggedTrackAnalysis(){
  //
  // Destructor
  //
  if(TestBit(kIsOwner)){
    if(fVarManager) delete fVarManager;
    if(fPIDqa) delete fPIDqa;
  }
  if(TestBit(kIsOwnerCuts)) delete fCuts;
  if(fContainer) delete fContainer;
}

//____________________________________________________________
void AliHFEtaggedTrackAnalysis::InitContainer(){
  //
  // Initialize output container
  //
  if(fContainer) return;
  Int_t nStepPID = 0;
  if(!fPID){
    AliError("No PID set - defining container without PID steps");
  } else {
    nStepPID = fPID->GetNumberOfPIDdetectors();
  }
  fContainer = new AliHFEcontainer("containerV0");
  fVarManager->DefineVariables(fContainer);
  fContainer->CreateContainer("taggedTrackContainerReco", "Container for Tagged Tracks", AliHFEcuts::kNcutStepsRecTrack + nStepPID);

  // Set the step titles
  for(Int_t istep = 0; istep < AliHFEcuts::kNcutStepsRecTrack; istep++)
    fContainer->SetStepTitle("taggedTrackContainerReco", AliHFEcuts::RecoCutName(istep), istep);
  for(Int_t ipid = 0; ipid < nStepPID; ipid++){
    fContainer->SetStepTitle("taggedTrackContainerReco", fPID->SortedDetectorName(ipid), ipid + AliHFEcuts::kNcutStepsRecTrack);
  }
  fCFM->SetParticleContainer(fContainer->GetCFContainer("taggedTrackContainerReco"));

  // temporarily special QA
  fQAhistos = new AliHFEcollection("taggedTrackQA", "Special QA for the TaggedTrackAnalysis");
  fQAhistos->CreateTH2F("TPCclusters2_1", "TPCclusterInfo for findable clusters for 2 neighbors", 30, 0.1, 10., 162, 0., 161.);
  fQAhistos->CreateTH2F("TPCclusters2_0", "TPCclusterInfo for the ratio for 2 neighbors", 30, 0.1, 10., 110, 0., 1.1);
  fQAhistos->CreateTH2F("TPCncls", "TPC number of clusters", 30, 0.1, 10., 162, 0., 161.);
  fQAhistos->CreateTH2F("TPCclr", "TPC cluster ratio", 30, 0.1, 10., 110, 0., 1.1);
  fQAhistos->BinLogAxis("TPCclusters2_1", 0);   // pt axis in logarithmic binning
  fQAhistos->BinLogAxis("TPCclusters2_0", 0);   // pt axis in logarithmic binning
  fQAhistos->BinLogAxis("TPCncls", 0);   // pt axis in logarithmic binning
  fQAhistos->BinLogAxis("TPCclr", 0);   // pt axis in logarithmic binning
}

//____________________________________________________________
void AliHFEtaggedTrackAnalysis::ProcessTrack(AliVParticle *track, Int_t abinitioPID){
  //
  // Filter tracks tagged by V0 PID class
  //
  //
  fVarManager->NewTrack(track, NULL, fCentralityF, abinitioPID, kTRUE);



  // Phi Angle
  if(fVariablesTRD) {
    AliESDtrack *esdtrackc = dynamic_cast<AliESDtrack *>(track);
    if(esdtrackc) {
      
      const AliExternalTrackParam *trueparam = NULL;
      if(esdtrackc->GetOuterParam()) {
	trueparam = esdtrackc->GetOuterParam();
	fVarManager->NewTrack((AliVParticle *)trueparam, NULL, fCentralityF, abinitioPID, kTRUE);
      }
      else return;
    }
  }
  

  // Try a loose cut to reject pion contamination
  if(fClean) {
    if(abinitioPID == AliPID::kElectron){
      AliHFEpidTPC *pidTPC = (AliHFEpidTPC *) fPID->GetDetPID(AliHFEpid::kTPCpid);
      if(pidTPC) {
	      Double_t numberOfSigmaTPC = pidTPC->GetPIDResponse()->NumberOfSigmasTPC(track,AliPID::kElectron);
	      if(numberOfSigmaTPC < -5) return;
      }
    }
  }
  // temporarily monitoring of the number of TPC clusters 
  AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
  if(esdtrack && abinitioPID == AliPID::kElectron){
    if((esdtrack->GetITSClusterMap() & (BIT(0) | BIT(1))) && (TMath::Abs(esdtrack->Eta()) < 0.8)){  // Only select quasi-primary tracks
      fQAhistos->Fill("TPCclusters2_1", track->Pt(), esdtrack->GetTPCClusterInfo(2,1));
      fQAhistos->Fill("TPCclusters2_0", track->Pt(), esdtrack->GetTPCNclsF() > 0 ? esdtrack->GetTPCClusterInfo(2,1)/esdtrack->GetTPCNclsF() : 0.);
      fQAhistos->Fill("TPCncls", track->Pt(), esdtrack->GetTPCNcls());
      fQAhistos->Fill("TPCclr", track->Pt(), esdtrack->GetTPCNclsF() > 0 ? static_cast<Double_t>(esdtrack->GetTPCNcls())/static_cast<Double_t>(esdtrack->GetTPCNclsF()) : 0.);
    }
   }
   
  Int_t offset = AliHFEcuts::kStepRecKineITSTPC;
  fVarManager->FillContainer(fCFM->GetParticleContainer(), 0); // Fill Container without filtering
  
  Bool_t survived = kTRUE;
  for(Int_t icut = AliHFEcuts::kStepRecKineITSTPC; icut <= AliHFEcuts::kStepHFEcutsTRD; icut++){
    AliDebug(2, Form("Checking cut %d for species %s", icut + AliHFEcuts::kNcutStepsMCTrack, AliPID::ParticleName(abinitioPID)));
    /*
      TObjArray *cutlist = fCFM->GetParticleCutsList(icut + AliHFEcuts::kNcutStepsMCTrack);
      if(!cutlist){
      AliDebug(2, Form("No cuts for step %d set", icut + AliHFEcuts::kNcutStepsMCTrack));
      } else {
      AliDebug(2, Form("Cut Collection %s", cutlist->GetName()));
      TIter cutiter(cutlist);
      AliCFCutBase *cut;
      while((cut = dynamic_cast<AliCFCutBase *>(cutiter()))){
      AliDebug(2, Form("Cut object %s, QA on? %s", cut->GetName(), cut->IsQAOn() ? "yes" : "no"));
      }
      }
    */
    //if(!fCFM->CheckParticleCuts(icut + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)){
    if(!fCuts->CheckParticleCuts(icut + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)){
      AliDebug(2, Form("Track didn' survive cut %d", icut + AliHFEcuts::kNcutStepsMCTrack));
      survived = kFALSE;
      break;
    }
    AliDebug(2, Form("Cut passed, filling container %d", icut - offset + 1));
    fVarManager->FillContainer(fCFM->GetParticleContainer(), icut - offset + 1);
  }
  
   if(survived){
     AliDebug(2, "Use track in the PID");
     // Apply PID
     AliHFEpidObject hfetrack;
     hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
     hfetrack.SetRecTrack(track);
     hfetrack.SetAbInitioPID(abinitioPID);
     hfetrack.SetCentrality(fCentralityF);
     if(fIsPbPb) hfetrack.SetPbPb();
     else hfetrack.SetPP();
     fPID->SetVarManager(fVarManager);
     fPID->IsSelected(&hfetrack, fContainer, "taggedTrackContainer", fPIDqa);
   }
}

//____________________________________________________________
void AliHFEtaggedTrackAnalysis::SetCuts(AliHFEcuts *cuts){
  //
  // Set HFE cuts to be used to filter the tagged tracks
  //
  if(!cuts){
    AliWarning("Nob cuts provided - Using standard cuts");
    fCuts = new AliHFEcuts("cutsTagged", "HFE Cuts for the V0 tagged tracks");
    fCuts->CreateStandardCuts();
    fCuts->SetQAOn();
    SetBit(kIsOwnerCuts);
  } else {
    AliDebug(1, "Setting single track cuts");
    fCuts = cuts;
  }
  const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack;
  printf("Setting Number of cut steps %d\n", kNcutSteps);
  fCFM->SetNStepParticle(kNcutSteps);
  for(Int_t istep = 0; istep < kNcutSteps; istep++)
    fCFM->SetParticleCutsList(istep, NULL);

  fCuts->Initialize(fCFM); 
}

//____________________________________________________________
void AliHFEtaggedTrackAnalysis::SetPID(AliHFEpid *pid){
  //
  // Set the PID and initialize the the QA manager
  //
  fPID = pid;
  fPIDqa->Initialize(fPID);
}

//____________________________________________________________
TList *AliHFEtaggedTrackAnalysis::GetPIDQA() const {
  //
  // return PID QA 
  //
  return fPIDqa->MakeList("PIDqa_taggedTracks");
}

//____________________________________________________________
TList *AliHFEtaggedTrackAnalysis::GetCutQA() const {
  //
  // return Cut QA
  //
  return fCuts->GetQAhistograms();
}


