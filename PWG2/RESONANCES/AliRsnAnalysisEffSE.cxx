//
// Class AliRsnAnalysisEffSE
//
// TODO
// TODO
//
// authors: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          Martin Vala (martin.vala@cern.ch)
//
#include <Riostream.h>
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"

#include "AliCFContainer.h"
#include "AliRsnCutManager.h"
#include "AliRsnValue.h"
#include "AliRsnAnalysisEffSE.h"
#include "AliRsnPairDef.h"
#include "AliRsnCutSet.h"

ClassImp(AliRsnAnalysisEffSE)

//_____________________________________________________________________________
AliRsnAnalysisEffSE::AliRsnAnalysisEffSE(const char *name) :
  AliRsnVAnalysisTaskSE(name),
  fUseITSSA(kTRUE),
  fUseGlobal(kTRUE),
  fStepListMC(0),
  fStepListESD(0),
  fAxisList("AliRsnValue", 0),
  fPairDefList(0),
  fContainerList(0x0),
  fOutList(0x0),
  fVar(0),
  fPair(),
  fEventCuts("eventCuts", AliRsnCut::kEvent)
{
//
// Default constructor.
//

  AliDebug(AliLog::kDebug+2,"<-");

  DefineOutput(2, TList::Class());

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
AliRsnAnalysisEffSE::AliRsnAnalysisEffSE(const AliRsnAnalysisEffSE& copy) :
  AliRsnVAnalysisTaskSE(copy),
  fUseITSSA(copy.fUseITSSA),
  fUseGlobal(copy.fUseGlobal),
  fStepListMC(copy.fStepListMC),
  fStepListESD(copy.fStepListESD),
  fAxisList(copy.fAxisList),
  fPairDefList(copy.fPairDefList),
  fContainerList(copy.fContainerList),
  fOutList(0x0),
  fVar(0),
  fPair(),
  fEventCuts(copy.fEventCuts)
{
//
// Copy constrtuctor
//
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::AddStepMC(AliRsnCutManager *mgr)
{
//
// Add a step on montecarlo
//

  fStepListMC.AddLast(mgr);
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::AddStepESD(AliRsnCutManager *mgr) 
{
//
// Add a step on ESD
//

  fStepListESD.AddLast(mgr);
}
    
//_____________________________________________________________________________
void AliRsnAnalysisEffSE::AddAxis(AliRsnValue *axis) 
{
//
// Add a new axis
//

  Int_t n = fAxisList.GetEntries();
  new (fAxisList[n]) AliRsnValue(*axis);
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::RsnUserCreateOutputObjects()
{
//
// Creation of output objects.
// These are created through the utility methods in the analysis manager,
// which produces a list of histograms for each specified set of pairs.
// Each of these lists is added to the main list of this task.
//

  AliDebug(AliLog::kDebug+2,"<-");

  // get number of steps and axes
  Int_t iaxis  = 0;
  Int_t nAxes  = fAxisList.GetEntries();
  Int_t nSteps = (Int_t)fStepListMC.GetEntries() + (Int_t)fStepListESD.GetEntries();

  if (!nSteps) {
    AliError("No steps defined");
    return;
  }
  if (!nAxes) {
    AliError("No axes defined");
    return;
  }

  // initialize variable list
  fVar.Set(nAxes);

  // retrieve number of bins for each axis
  Int_t *nBins = new Int_t[nAxes];
  for (iaxis = 0; iaxis < nAxes; iaxis++) 
  {
    AliRsnValue *fcnAxis = (AliRsnValue*)fAxisList.At(iaxis);
    nBins[iaxis] = fcnAxis->GetArray().GetSize() - 1;
  }

  // create ouput list of containers
  fContainerList = new TList;
  fContainerList->SetOwner();
  fContainerList->SetName(Form("%s_containers", GetName()));

  // initialize output list
  OpenFile(2);
  fOutList = new TList();
  fOutList->SetOwner();

  // create the containers
  Int_t i = 0, nDef = (Int_t)fPairDefList.GetEntries();
  for (i = 0; i < nDef; i++) 
  {
    AliRsnPairDef *def = (AliRsnPairDef*)fPairDefList[i];
    AliCFContainer *cont = new AliCFContainer(Form("%s", def->GetPairName()), "", nSteps, nAxes, nBins);
    // set the bin limits for each axis
    for (iaxis = 0; iaxis < nAxes; iaxis++) 
    {
      AliRsnValue *fcnAxis = (AliRsnValue*)fAxisList.At(iaxis);
      cont->SetBinLimits(iaxis, fcnAxis->GetArray().GetArray());
    }
    // add the container to output list
    fContainerList->Add(cont);
  }

  fOutList->Add(fContainerList);
  fOutList->Print();

  PostData(2, fOutList);
  
  // clear heap
  delete [] nBins;

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::RsnUserExec(Option_t*)
{
//
// Execution of the analysis task.
// Recovers the input event and processes it with all included pair objects.
// In this case, we NEED to have ESD and MC, otherwise we cannod do anything.
//ULTIMO UNDO
/*
  AliRsnDaughter trk;
  for (Int_t i = 0; i <= AliPID::kSPECIES; i++)
  {
    cout << AliPID::ParticleName((AliPID::EParticleType)i) << ": " << endl;
    for (Int_t m = 0; m < AliRsnDaughter::kMethods; m++)
    {
      cout << "-- method: " << AliRsnDaughter::MethodName((AliRsnDaughter::EPIDMethod)m) << endl;
      Char_t   sign[2] = {'+', '-'};
      for (Int_t s = 0; s < 2; s++)
      {
        TArrayI *a = fRsnPIDIndex.GetTracksArray((AliRsnDaughter::EPIDMethod)m, sign[s], (AliPID::EParticleType)i);
        Int_t n = a->GetSize();
        for (Int_t j = 0; j < n; j++)
        {
          Int_t k = a->At(j);
          cout << "-- -- track " << Form("%4d ", k) << ": ";
          fRsnEvent.SetDaughter(trk, k);
          cout << "charge = " << (trk.IsPos() ? "+ " : (trk.IsNeg() ? "- " : "0 "));
          cout << "truePID = " << Form("%10s ", AliPID::ParticleName(trk.PerfectPID()));
          cout << "realPID = " << Form("%10s ", AliPID::ParticleName(trk.RealisticPID()));
          cout << endl;
          cout << "-- -- -- weights (computed): ";
          for (Int_t q = 0; q < AliPID::kSPECIES; q++)
            cout << Form("%15.12f", trk.ComputedWeights()[q]) << ' ';
          cout << endl;
          cout << "-- -- -- weights (original): ";
          for (Int_t q = 0; q < AliPID::kSPECIES; q++)
            cout << Form("%15.12f", trk.GetRef()->PID()[q]) << ' ';
          cout << endl;
        }
      }
    }
  } return;
  */

  AliDebug(AliLog::kDebug+2,"<-");
  if (!fESDEvent || !fMCEvent) {
    AliError("This task can process only ESD + MC events");
    return;
  }
  fRsnEvent.SetRef(fESDEvent, fMCEvent);

  // if general event cuts are added to the task (recommended)
  // they are checked here on the RSN event interface and,
  // if the event does not pass them, it is skipped and ProcessInfo
  // is updated accordingly
  if (!fEventCuts.IsSelected(&fRsnEvent)) {
    fTaskInfo.SetEventUsed(kFALSE);
    return;
  }
  
  // if cuts are passed or not cuts were defined,
  // update the task info before processing the event
  fTaskInfo.SetEventUsed(kTRUE);

  // process first MC steps and then ESD steps
  AliRsnPairDef *pairDef = 0;
  TObjArrayIter iter(&fPairDefList);
  while ( (pairDef = (AliRsnPairDef*)iter.Next()) )
  {
    //ProcessEventMC(pairDef);
    //ProcessEventESD(pairDef);
    ProcessEvent(pairDef);
  }

  // Post the data
  PostData(2, fOutList);

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::ProcessEvent(AliRsnPairDef *pairDef)
{
//
// Process current event with the definitions of the specified step in MC list
// and store results in the container slot defined in second argument.
// It is associated with the AliCFContainer with the name of the pair.
//

  AliStack      *stack = fRsnEvent.GetRefMC()->Stack();
  AliESDEvent   *esd   = fRsnEvent.GetRefESD();
  AliMCEvent    *mc    = fRsnEvent.GetRefMC();
  AliMCParticle *mother;

  if (!pairDef) return;
  AliCFContainer *cont = (AliCFContainer*)fContainerList->FindObject(pairDef->GetPairName());
  if (!cont) return;
  
  // get informations from pairDef
  Int_t pdgM = 0, pdgD[2] = {0, 0};
  Short_t chargeD[2] = {0, 0};
  pdgM    = pairDef->GetMotherPDG();
  pdgD[0] = AliPID::ParticleCode(pairDef->GetPID(0));
  pdgD[1] = AliPID::ParticleCode(pairDef->GetPID(1));
  chargeD[0] = pairDef->GetChargeShort(0);
  chargeD[1] = pairDef->GetChargeShort(1);

  // other utility variables
  Int_t   label[2] = {-1, -1}, first, j, ipart;
  Short_t charge[2] = {0, 0};
  Short_t pairDefMatch[2] = {-1, -1};
  Int_t   esdIndex[2];
  TParticle *part[2] = {0, 0};

  // in this case, we first take the resonance from MC
  // and then we find its daughters and compute cuts on them
  for (ipart = 0; ipart < stack->GetNprimary(); ipart++) 
  {
    // take a track from the MC event
    mother = (AliMCParticle*) fMCEvent->GetTrack(ipart);
    
    // check that it is a binary decay and the PDG code matches
    if (mother->Particle()->GetNDaughters() != 2) continue;
    if (mother->Particle()->GetPdgCode() != pdgM) continue;

    // store the labels of the two daughters
    label[0] = mother->Particle()->GetFirstDaughter();
    label[1] = mother->Particle()->GetLastDaughter();
    
    // retrieve the particles and other stuff
    // check if they match the order in the pairDef
    for (j = 0; j < 2; j++)
    {
      if (label[j] < 0) continue;
      part[j]   = stack->Particle(label[j]);
      pdgD[j]   = TMath::Abs(part[j]->GetPdgCode());
      charge[j] = (Short_t)(part[j]->GetPDG()->Charge() / 3);
      if (pdgD[j] == AliPID::ParticleCode(pairDef->GetPID(0)) && charge[j] == pairDef->GetChargeShort(0))
        pairDefMatch[j] = 0;
      else if (pdgD[j] == AliPID::ParticleCode(pairDef->GetPID(1)) && charge[j] == pairDef->GetChargeShort(1))
        pairDefMatch[j] = 1;
      else
        pairDefMatch[j] = -1;
        
      // find corresponding ESD particle: first try rejecting fakes,
      // and in case of failure, try accepting fakes
      esdIndex[j] = FindESDtrack(label[j], esd, kTRUE);
      //TArrayI idx = FindESDtracks(label[j], esd);
      //for (Int_t kk = 0; kk < idx.GetSize(); kk++) cout << "DAUGHTER " << j << " --> FOUND INDEX: " << idx[kk] << endl;
      if (esdIndex[j] < 0) esdIndex[j] = FindESDtrack(label[j], esd, kFALSE);
      //cout << "DAUGHTER " << j << " SINGLE FOUND INDEX = " << esdIndex[j] << endl;
    }
    
    // since each candidate good resonance is taken once, we must check
    // that it matches the pair definition in any order, and reject in case
    // in none of them the pair is OK
    // anyway, we must associate the correct daughter to the correct data member
    if (pairDefMatch[0] == 0 && pairDefMatch[1] == 1)
    {
      // 1st track --> 1st member of PairDef
      fDaughter[0].SetRef(mc->GetTrack(label[0]));
      fDaughter[0].SetRefMC((AliMCParticle*)mc->GetTrack(label[0]));
      fDaughter[0].SetGood();
      // 2nd track --> 2nd member of PairDef
      fDaughter[1].SetRef(mc->GetTrack(label[1]));
      fDaughter[1].SetRefMC((AliMCParticle*)mc->GetTrack(label[1]));
      fDaughter[1].SetGood();
    }
    else if ((pairDefMatch[0] == 1 && pairDefMatch[1] == 0))
    {
      // 1st track --> 2nd member of PairDef
      fDaughter[0].SetRef(mc->GetTrack(label[1]));
      fDaughter[0].SetRefMC((AliMCParticle*)mc->GetTrack(label[1]));
      fDaughter[0].SetGood();
      // 2nd track --> 1st member of PairDef
      fDaughter[1].SetRef(mc->GetTrack(label[0]));
      fDaughter[1].SetRefMC((AliMCParticle*)mc->GetTrack(label[0]));
      fDaughter[1].SetGood();
    }
    else
    {
      fDaughter[0].SetBad();
      fDaughter[1].SetBad();
    }
    
    // reject the pair if the matching was unsuccessful
    if (!fDaughter[0].IsOK() || !fDaughter[1].IsOK()) continue;
    
    // first, we set the internal AliRsnMother object to
    // the MC particles and then operate the selections on MC
    fPair.SetDaughters(&fDaughter[0], pairDef->GetMass(0), &fDaughter[1], pairDef->GetMass(1));
    FillContainer(cont, &fStepListMC, pairDef, 0);
    
    // then, if both particles found a good match in the ESD
    // reassociate the ESD tracks to the pair and fill ESD containers
    if (esdIndex[0] < 0 || esdIndex[1] < 0) continue;
    if (pairDefMatch[0] == 0 && pairDefMatch[1] == 1)
    {
      // 1st track --> 1st member of PairDef
      fDaughter[0].SetRef(esd->GetTrack(esdIndex[0]));
      // 2nd track --> 2nd member of PairDef
      fDaughter[1].SetRef(esd->GetTrack(esdIndex[1]));
      //cout << "****** MATCHING SCHEME 1" << endl;
    }
    else if ((pairDefMatch[0] == 1 && pairDefMatch[1] == 0))
    {
      // 1st track --> 2nd member of PairDef
      fDaughter[0].SetRef(esd->GetTrack(esdIndex[1]));
      // 2nd track --> 1st member of PairDef
      fDaughter[1].SetRef(esd->GetTrack(esdIndex[0]));
      //cout << "****** MATCHING SCHEME 2" << endl;
    }
    //cout << "****** IDs = " << fDaughter[0].GetID() << ' ' << fDaughter[1].GetID() << endl;
    // here we must remember how many steps were already filled
    first = (Int_t)fStepListMC.GetEntries();
    FillContainer(cont, &fStepListESD, pairDef, first);
  }
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::ProcessEventMC(AliRsnPairDef */*pairDef*/)
{
//
// Process current event with the definitions of the specified step in MC list
// and store results in the container slot defined in second argument
//

  /*
  AliStack      *stack = fMCEvent->Stack();
  AliMCParticle *mother, *daughter;

  if (!pairDef) return;
  AliCFContainer *cont = (AliCFContainer*)fContainerList->FindObject(pairDef->GetPairName().Data());
  if (!cont) return;

  // other utility variables
  Int_t i[2], j, ipart;

  // in this case, we first take the resonance from MC
  // and then we find its daughters and compute cuts on them
  for (ipart = 0; ipart < stack->GetNprimary(); ipart++) 
  {
    mother = (AliMCParticle*) fMCEvent->GetTrack(ipart);
    if (mother->Particle()->GetNDaughters() != 2) continue;

    i[0] = mother->Particle()->GetFirstDaughter();
    i[1] = mother->Particle()->GetLastDaughter();

    for (j = 0; j < 2; j++) 
    {
      daughter = (AliMCParticle*) fMCEvent->GetTrack(i[j]);
      fDaughter[j].SetRef(daughter);
      fDaughter[j].SetRefMC(daughter);
      fDaughter[j].FindMotherPDG(stack);
    }

    if (fDaughter[0].ChargeC() != pairDef->GetCharge(0)) continue;
    if (fDaughter[1].ChargeC() != pairDef->GetCharge(1)) continue;
    if (fDaughter[0].PerfectPID() != pairDef->GetType(0)) continue;
    if (fDaughter[1].PerfectPID() != pairDef->GetType(1)) continue;

    fPair.SetDaughters(&fDaughter[0], &fDaughter[1]);

    // create pair
    FillContainer(cont, &fStepListMC, pairDef, 0);
  }
  */
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::ProcessEventESD(AliRsnPairDef */*pairDef*/)
{
//
// Process current event with the definitions of the specified step in ESD list
// and store results in the container slot defined in second argument
//

  /*
  Int_t i0, i1, first = (Int_t)fStepListMC.GetEntries();

  if (!pairDef) return;
  AliCFContainer *cont = (AliCFContainer*)fContainerList->FindObject(pairDef->GetPairName().Data());
  if (!cont) return;

  // external loop on tracks
  for (i0 = 0; i0 < a0->GetSize(); i0++) {
    // connect interface
    fRsnEvent.SetDaughter(fDaughter[0], a0->At(i0));
    if (!fDaughter[0].IsOK()) continue;
    fDaughter[0].SetRequiredPID(pairDef->GetType(0));

    // internal loop on tracks
    for (i1 = 0; i1 < a1->GetSize(); i1++) {
      // connect interface
      fRsnEvent.SetDaughter(fDaughter[1], a1->At(i1));
      if (!fDaughter[1].IsOK()) continue;
      fDaughter[1].SetRequiredPID(pairDef->GetType(1));
      // build pair
      fPair.SetPair(&fDaughter[0], &fDaughter[1]);
      if (TMath::Abs(fPair.CommonMother()) != pairDef->GetMotherPDG()) continue;
      // fill containers
      FillContainer(cont, &fStepListESD, pairDef, first);
    }
  }
  */
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::FillContainer(AliCFContainer *cont, const TObjArray *stepList, AliRsnPairDef *pd, Int_t firstOutStep)
{
//
// Fill the containers
//

  Int_t iaxis, nAxes  = fAxisList.GetEntries();
  Int_t istep, nSteps = stepList->GetEntries();
  
  // set daughters to pair
  fPair.SetDaughters(&fDaughter[0], pd->GetMass(0), &fDaughter[1], pd->GetMass(1));

  // compute values for all axes
  for (iaxis = 0; iaxis < nAxes; iaxis++) 
  {
    AliRsnValue *fcnAxis = (AliRsnValue*)fAxisList.At(iaxis);
    fVar[iaxis] = -1E10;
    if (fcnAxis->Eval(&fPair, pd, &fRsnEvent)) fVar[iaxis] = (Float_t)fcnAxis->GetValue();
  }

  // fill all steps
  for (istep = 0; istep < nSteps; istep++) 
  {
    AliRsnCutManager *cutMgr = (AliRsnCutManager*)stepList->At(istep);
    cutMgr->SetEvent(&fRsnEvent);
    if (!cutMgr->PassCommonDaughterCuts(&fDaughter[0])) break;
    if (!cutMgr->PassCommonDaughterCuts(&fDaughter[1])) break;
    if (!cutMgr->PassDaughter1Cuts(&fDaughter[0])) break;
    if (!cutMgr->PassDaughter2Cuts(&fDaughter[1])) break;
    if (!cutMgr->PassMotherCuts(&fPair)) break;
    //cout << "**************************************** FILLING STEP " << istep << endl;
    cont->Fill(fVar.GetArray(), istep + firstOutStep);
  }
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::RsnTerminate(Option_t*)
{
//
// Termination.
// Could be added some monitor histograms here.
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::AddPairDef(AliRsnPairDef* pairDef)
{
//
//  Adds pair definition
//
  fPairDefList.AddLast(pairDef);
}

//_____________________________________________________________________________
Int_t AliRsnAnalysisEffSE::FindESDtrack(Int_t label, AliESDEvent *esd, Bool_t rejectFakes)
{
//
// Finds in the ESD a track whose label corresponds to that in argument.
// When global tracks are enabled, tries first to find a global track 
// satisfying that requirement.
// If no global tracks are found, if ITS-SA are enable, tries to search among them
// otherwise return a negative number.
// If global tracks are disabled, search only among ITS SA
//

  Int_t   i = 0;
  Int_t   ntracks = esd->GetNumberOfTracks();
  ULong_t status;
  Bool_t  isTPC;
  Bool_t  isITSSA;
  
  // loop for global tracks
  if (fUseGlobal)
  {
    for (i = 0; i < ntracks; i++)
    {
      AliESDtrack *track = esd->GetTrack(i);
      status  = (ULong_t)track->GetStatus();
      isTPC   = ((status & AliESDtrack::kTPCin)  != 0);
      if (!isTPC) continue;
      
      // check that label match
      if (TMath::Abs(track->GetLabel()) != label) continue;
      
      // if required, reject fakes
      if (rejectFakes && track->GetLabel() < 0) continue;
      
      // if all checks are passed and we are searching among global
      // this means that thie track is a global one with the right label
      // then, the return value is set to this, and returned
      return i;
    }
  }
  
  // loop for ITS-SA tracks (this happens only if no global tracks were found
  // or searching among globals is disabled)
  if (fUseITSSA)
  {
    for (i = 0; i < ntracks; i++)
    {
      AliESDtrack *track = esd->GetTrack(i);
      status  = (ULong_t)track->GetStatus();
      isITSSA = ((status & AliESDtrack::kTPCin)  == 0 && (status & AliESDtrack::kITSrefit) != 0 && (status & AliESDtrack::kITSpureSA) == 0 && (status & AliESDtrack::kITSpid) != 0);
      if (!isITSSA) continue;
      
      // check that label match
      if (TMath::Abs(track->GetLabel()) != label) continue;
            
      // if required, reject fakes
      if (rejectFakes && track->GetLabel() < 0) continue;
      
      // if all checks are passed and we are searching among global
      // this means that thie track is a global one with the right label
      // then, the return value is set to this, and returned
      return i;
    }
  }
  
  // if we reach this point, no match were found
  return -1;
}

//_____________________________________________________________________________
TArrayI AliRsnAnalysisEffSE::FindESDtracks(Int_t label, AliESDEvent *esd)
{
//
// Finds in the ESD a track whose label corresponds to that in argument.
// When global tracks are enabled, tries first to find a global track 
// satisfying that requirement.
// If no global tracks are found, if ITS-SA are enable, tries to search among them
// otherwise return a negative number.
// If global tracks are disabled, search only among ITS SA
//

  Int_t   i = 0;
  Int_t   ntracks = esd->GetNumberOfTracks();
  ULong_t status;
  Bool_t  isTPC;
  Bool_t  isITSSA;
  TArrayI array(100);
  Int_t   nfound = 0;
  
  // loop for global tracks
  if (fUseGlobal)
  {
    for (i = 0; i < ntracks; i++)
    {
      AliESDtrack *track = esd->GetTrack(i);
      status  = (ULong_t)track->GetStatus();
      isTPC   = ((status & AliESDtrack::kTPCin)  != 0);
      if (!isTPC) continue;
      
      // check that label match
      if (TMath::Abs(track->GetLabel()) != label) continue;
      
      array[nfound++] = i;
    }
  }
  
  // loop for ITS-SA tracks (this happens only if no global tracks were found
  // or searching among globals is disabled)
  if (fUseITSSA)
  {
    for (i = 0; i < ntracks; i++)
    {
      AliESDtrack *track = esd->GetTrack(i);
      status  = (ULong_t)track->GetStatus();
      isITSSA = ((status & AliESDtrack::kTPCin)  == 0 && (status & AliESDtrack::kITSrefit) != 0 && (status & AliESDtrack::kITSpureSA) == 0 && (status & AliESDtrack::kITSpid) != 0);
      if (!isITSSA) continue;
      
      // check that label match
      if (TMath::Abs(track->GetLabel()) != label) continue;
            
      array[nfound++] = i;
    }
  }
  
  array.Set(nfound);
  return array;
}
