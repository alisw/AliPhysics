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

// Analysis task to produce trees of lightweight events
// evgeny.kryshen@cern.ch
// lmilano@cern.ch
// alice.ohlson@cern.ch
// igor.lakomov@cern.ch

// aliroot
#include "AliAnalysisTaskCFTree.h"
#include "AliCFParticle.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliCentrality.h"
#include "AliAnalysisFilter.h"
#include "AliVMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCParticle.h"
#include "AliCollisionGeometry.h"
#include "AliVVZERO.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliMuonTrackCuts.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"

#include "AliMultiplicity.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"

// root
#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TPythia6Decayer.h"
#include "TParticle.h"
#include "TBits.h"
ClassImp(AliAnalysisTaskCFTree)

//-----------------------------------------------------------------------------
AliAnalysisTaskCFTree::AliAnalysisTaskCFTree(const char* name) :
AliAnalysisTaskSE(name),
fTrackFilter(0x0),
fHybridConstrainedMask(0),
fTPConlyConstrainedMask(0),
fMuonTrackCuts(new AliMuonTrackCuts),
fUtils(0x0),
fitssatrackcuts(0x0),
fListOfHistos(0x0),
fEventStatistics(0x0),
fClassStatistics(0x0),
fV0chan(0x0),
fTree(0x0),
fTracks(0x0),
fTracklets(0x0),
fMuons(0x0),
fMcParticles(0x0),
fMcMuons(0x0),
fMuonOrigin(0x0),
fIs13TeV(kTRUE),
fClassesFired(0),
fField(0),
fCurrentRunNumber(0),
fCentrality(),
fVtxX(0),
fVtxY(0),
fVtxZ(0),
fVtxTPConly(0),
fVtxContributors(0),
fSPDVtxZ(0),
fSPDVtxZRes(0),
fSPDVtxContributors(0),
fPeriod(0),
fOrbit(),
fBc(),
fSelectMask(0),
fIsDiffractive(0),
fIsPileupSPD(0),
fIsPileupMV(0),
fIR1(0),
fIR2(0),
fonlineSPD(),
fofflineSPD(),
fV0ADecision(-999),
fV0CDecision(-999),
fIsIncomplete(0),
fNofITSClusters(),
fMultV0Aeq(0),
fMultV0Ceq(0),
fMultV0Meq(0),
fMultV0Aring(),
fMultV0Cring(),
fMultTKL(0),
fMultPercV0Aeq(-1.),
fMultPercV0Ceq(-1.),
fMultPercV0Meq(-1.),
fMultPercTKL(-1.),
fMultPercV0S(-1.),
fMultMeanV0A(-1.),
fMultMeanV0C(-1.),
fMultMeanV0M(-1.),
fMultMeanTKL(-1.),
fIsEventSel(-1),
fNchTPC(-1),
fNchTPCmc(-1),
fNchV0Amc(-1),
fNchV0Cmc(-1),
fNchCL1mc(-1),
fClassBit(0xffffffff),
fSelectBit(AliVEvent::kAny),
fZVertexCut(10.),
fTrackFilterBit(0xffffffff),
fTrackletEtaCut(2.0),
fTrackEtaCut(5.0),
fPtMin(0.15),
fSharedClusterCut(0.4),
fCrossedRowsCut(100),
fFoundFractionCut(0.8),
fDphiCut(1.e9),
fStoreTracks(0),
fStoreTracklets(0),
fStoreMuons(0),
fStoreMcTracks(0),
fStoreMcTracklets(0),
fStoreMcMuons(0),
fStoreTrackInfo(0),
fStorePidInfo(0),
fStoreAodDCAInfo(0),
fStoreMuonOrigin(0),
fApplyPhysicsSelectionCut(0),
fStoreOnlyEventsWithMuons(0),
fStoreCutBitsInTrackMask(0),
fDecayArray(0x0),
fDecayer(0x0),
fMapping(0x0)
{
  Info("AliAnalysisTaskCFTree","Calling Constructor");
  fMuonTrackCuts->SetCustomParamFromRun(197388,"muon_pass2");
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskCFTree::UserCreateOutputObjects(){
  fListOfHistos = new TList();
  fListOfHistos->SetOwner();
  fEventStatistics = new TH1I("fEventStatistics","",10,0,10);
  fClassStatistics = new TH1D("fClassStatistics","",16,0,16);

#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
  fEventStatistics->SetBit(TH1::kCanRebin);
#endif
  fListOfHistos->Add(fEventStatistics);
  fListOfHistos->Add(fClassStatistics);

  fV0chan = new TH1F("fV0chan","",64,-0.5,63.5);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
  fV0chan->SetBit(TH1::kCanRebin);
#endif
  fListOfHistos->Add(fV0chan);

  if (fStoreTracks)    fTracks      = new TClonesArray("AliCFParticle",2000);
  if (fStoreTracklets) fTracklets   = new TClonesArray("AliCFParticle",2000);
  if (fStoreMuons)     fMuons       = new TClonesArray("AliCFParticle",2000);
  if (fStoreMcTracks)  fMcParticles = new TClonesArray("AliCFParticle",2000);
  if (fStoreMuonOrigin)fMuonOrigin  = new TClonesArray("TObjString",2000);
  if (fStoreMcMuons)   fMcMuons     = new TClonesArray("AliCFParticle",2000);

  // create file-resident tree
  TDirectory *owd = gDirectory;
  OpenFile(2);
  fTree = new TTree("events","events");
  owd->cd();
  fTree->Branch("classes",&fClassesFired);
  fTree->Branch("cent",&fCentrality,"fCentrality[6]/F");
  fTree->Branch("vtxx",&fVtxX);
  fTree->Branch("vtxy",&fVtxY);
  fTree->Branch("vtxz",&fVtxZ);
  fTree->Branch("vtxTPConly",&fVtxTPConly);
  fTree->Branch("vtxContributors",&fVtxContributors);
  fTree->Branch("spdvtxz",&fSPDVtxZ);
  fTree->Branch("spdvtxzres",&fSPDVtxZRes);
  fTree->Branch("spdvtxContributors",&fSPDVtxContributors);
  fTree->Branch("field",&fField);
  fTree->Branch("run",&fCurrentRunNumber);
  fTree->Branch("period",&fPeriod);
  fTree->Branch("orbit",&fOrbit);
  fTree->Branch("bc",&fBc);
  fTree->Branch("mask",&fSelectMask);
  fTree->Branch("IsDiffractive",&fIsDiffractive);
  fTree->Branch("pileupspd",&fIsPileupSPD);
  fTree->Branch("pileupmv",&fIsPileupMV);
  fTree->Branch("IR1",&fIR1);
  fTree->Branch("IR2",&fIR2);
  fTree->Branch("onlineSPD",&fonlineSPD);
  fTree->Branch("offlineSPD",&fofflineSPD);
  fTree->Branch("V0ADecision",&fV0ADecision);
  fTree->Branch("V0CDecision",&fV0CDecision);
  fTree->Branch("IsIncomplete",&fIsIncomplete);
  fTree->Branch("nofITSClusters",&fNofITSClusters,"fNofITSClusters[6]/I");
  fTree->Branch("multV0Aeq",&fMultV0Aeq);
  fTree->Branch("multV0Ceq",&fMultV0Ceq);
  fTree->Branch("multV0Meq",&fMultV0Meq);
  fTree->Branch("multV0Aring",&fMultV0Aring,"fMultV0Aring[4]/F");
  fTree->Branch("multV0Cring",&fMultV0Cring,"fMultV0Cring[4]/F");
  fTree->Branch("multTKL",&fMultTKL);
  fTree->Branch("multpercV0Aeq",&fMultPercV0Aeq);
  fTree->Branch("multpercV0Ceq",&fMultPercV0Ceq);
  fTree->Branch("multpercV0Meq",&fMultPercV0Meq);
  fTree->Branch("multpercTKL",&fMultPercTKL);
  fTree->Branch("multpercV0S",&fMultPercV0S);
  fTree->Branch("multmeanV0A",&fMultMeanV0A);
  fTree->Branch("multmeanV0C",&fMultMeanV0C);
  fTree->Branch("multmeanV0M",&fMultMeanV0M);
  fTree->Branch("multmeanTKL",&fMultMeanTKL);
  fTree->Branch("iseventsel",&fIsEventSel);
  fTree->Branch("nchTPC",&fNchTPC);
  fTree->Branch("nchTPCmc",&fNchTPCmc);
  fTree->Branch("nchV0Amc",&fNchV0Amc);
  fTree->Branch("nchV0Cmc",&fNchV0Cmc);
  fTree->Branch("nchCL1mc",&fNchCL1mc);

  if (fTracks)      fTree->Branch("tracks",&fTracks);
  if (fTracklets)   fTree->Branch("tracklets",&fTracklets);
  if (fMuons)       fTree->Branch("muons",&fMuons);
  if (fMcParticles) fTree->Branch("mcparticles",&fMcParticles);
  if (fMuonOrigin)  fTree->Branch("muon_origin",&fMuonOrigin);
  if (fMcMuons)     fTree->Branch("mcmuons",&fMcMuons);

  Int_t iParameter=0; // Mapping Tracks
  for (Int_t i=0; i<10; i++) {
    if (fStoreTrackInfo) fMapping->MappingTracks()[i]=iParameter++; else fMapping->MappingTracks()[i]=-1;
  }
  for (Int_t i=10; i<13; i++) {
    if (fStorePidInfo) fMapping->MappingTracks()[i]=iParameter++; else fMapping->MappingTracks()[i]=-1;
  }
  for (Int_t i=13; i<14; i++) {
    if (fStoreMcTracks) fMapping->MappingTracks()[i]=iParameter++; else fMapping->MappingTracks()[i]=-1;
  }
  for (Int_t i=14; i<AliCFTreeMapping::kMappingTracks; i++) {
    if (fStoreAodDCAInfo) fMapping->MappingTracks()[i]=iParameter++; else fMapping->MappingTracks()[i]=-1;
  }

  iParameter=0; // Mapping Tracklets
  for (Int_t i=0; i<AliCFTreeMapping::kMappingTracklets; i++) {
      if (fStoreTracklets) fMapping->MappingTracklets()[i]=iParameter++; else fMapping->MappingTracklets()[i]=-1;
  }

  iParameter=0; // Mapping Muons
  for (Int_t i=0; i<4; i++) {
      if (fStoreMuons) fMapping->MappingMuons()[i]=iParameter++; else fMapping->MappingMuons()[i]=-1;
  }
  for (Int_t i=4; i<AliCFTreeMapping::kMappingMuons; i++) {
      if (fStoreMcMuons) fMapping->MappingMuons()[i]=iParameter++; else fMapping->MappingMuons()[i]=-1;
  }

  iParameter=0; // Mapping MCTracks
  for (Int_t i=0; i<AliCFTreeMapping::kMappingMCTracks; i++) {
    if (fStoreMcTracks) fMapping->MappingMCTracks()[i]=iParameter++; else fMapping->MappingMCTracks()[i]=-1;
  }

  fTree->GetUserInfo()->Add(fMapping);	//to retreive it afterwards one needs fTree->GetUserInfo()->At(0)

  fUtils = new AliAnalysisUtils();
  fUtils->SetUseSPDCutInMultBins(kTRUE);
  //  fUtils->SetMinPlpContribSPD(3);
  fUtils->SetMinPlpContribMV(3);

  fitssatrackcuts = AliAODITSsaTrackCuts::GetStandardAODITSsaTrackCuts2015();

  if(fMcMuons) {
    Int_t products[4] = {13,211,111,22};
    Int_t mult[4]     = {1 ,  1,  1, 1};
    Int_t npart =3;
    fDecayer = new TPythia6Decayer();
    fDecayer->ForceParticleDecay(321,products,mult,npart);
    fDecayer->ForceParticleDecay(130,products,mult,npart);
    fDecayer->ForceParticleDecay(211,13,1);
    fDecayArray = new TClonesArray("TParticle",10);
  }

  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskCFTree::UserExec(Option_t *){
  fEventStatistics->Fill("before cuts",1);

  if (fInputEvent) {
    fEventStatistics->Fill("after event check",1);

    //filter out incomplete events:
    if (((AliAODEvent*)fInputEvent)->IsIncompleteDAQ() && !fMCEvent) return;
    fIsIncomplete=0;
    //additional check on incomplete events (needed on old productions)
    if(fInputEvent->GetHeader()->GetL0TriggerInputs()==0 && !fMCEvent)fIsIncomplete=1;
    fEventStatistics->Fill("after incomplete event check",1);

    TString classes = fInputEvent->GetFiredTriggerClasses();
    //    Printf("An event: %s",classes.Data());
    fClassesFired = 0;
    if(fIs13TeV) {
      fClassesFired |= (classes.Contains("CINT7-B-")              << 0);
      fClassesFired |= (classes.Contains("CMSL7-B-")              << 1);
      fClassesFired |= (classes.Contains("CMSH7-B-")              << 2);
      fClassesFired |= (classes.Contains("CVHMV0M-B-")            << 3);
      fClassesFired |= (classes.Contains("CVHMSH2-B-")            << 4);
      fClassesFired |= (classes.Contains("CVHMV0MMSL-B-")         << 5);
      fClassesFired |= (classes.Contains("CVHMSH2MSL-B-")         << 6);
    }
    else {
      fClassesFired |= (classes.Contains("CINT7-S-") << 0);
      fClassesFired |= (classes.Contains("CINT7-B-") << 0);
      fClassesFired |= (classes.Contains("CSHM8-S-") << 1);
      fClassesFired |= (classes.Contains("CSHM8-B-") << 1);
      fClassesFired |= (classes.Contains("CMSL7-S-") << 2);
      fClassesFired |= (classes.Contains("CMSL7-B-") << 2);
      fClassesFired |= (classes.Contains("CMSH7-S-") << 3);
      fClassesFired |= (classes.Contains("CMSH7-B-") << 3);
      fClassesFired |= (classes.Contains("CMUL7-S-") << 4);
      fClassesFired |= (classes.Contains("CMUL7-B-") << 4);
      fClassesFired |= (classes.Contains("CMLL7-S-") << 5);
      fClassesFired |= (classes.Contains("CMLL7-B-") << 5);
      fClassesFired |= (classes.Contains("CINT8-S-") << 6);
      fClassesFired |= (classes.Contains("CSPI7-S-") << 7);
      fClassesFired |= (classes.Contains("CSPI8-S-") << 8);
      fClassesFired |= (classes.Contains("CMSL8-S-") << 9);
      fClassesFired |= (classes.Contains("CMSH8-S-") <<10);
      fClassesFired |= (classes.Contains("CMUL8-S-") <<11);
      fClassesFired |= (classes.Contains("CMLL8-S-") <<12);
      fClassesFired |= (classes.Contains("C0MUL-SA") <<13);
      fClassesFired |= (classes.Contains("C0MUL-SC") <<14);
      fClassesFired |= (classes.Contains("CINT1-B-") <<15);
      fClassesFired |= (classes.Contains("CINT1B-ABCE-") <<15);
      fClassesFired |= (classes.Contains("CSH1-B-") <<16);
    }
    fClassStatistics->Fill(fClassesFired+0.001);


    if (!fApplyPhysicsSelectionCut && !(fClassesFired & fClassBit) && !fMCEvent) return;
    fEventStatistics->Fill("after trigger check",1);

    //      Printf("an event: %s    (%d)",classes.Data(),fClassesFired);

    fSelectMask = fInputHandler->IsEventSelected();
    if (fApplyPhysicsSelectionCut && !(fSelectMask & fSelectBit)) return;
    fEventStatistics->Fill("after physics selection",1);

    //They should be equal to 1 for good events.
    AliVVZERO* vzero = fInputEvent->GetVZEROData();
    fV0ADecision = (Int_t)vzero->GetV0ADecision();
    fV0CDecision = (Int_t)vzero->GetV0CDecision();

    if(fMcParticles)
    {
      if(fMCEvent)
      {
        if(fMCEvent->GetNumberOfTracks()==0)
          return;
        fEventStatistics->Fill("generated mc particles > 0",1);
      }
    }

    fPeriod        = fInputEvent->GetPeriodNumber();
    fOrbit         = fInputEvent->GetOrbitNumber();
    fBc            = fInputEvent->GetBunchCrossNumber();
    fField         = fInputEvent->GetMagneticField();
    fCurrentRunNumber = fInputEvent->GetRunNumber() ;
    fCentrality[0] = fInputEvent->GetCentrality()->GetCentralityPercentile("V0M");
    fCentrality[1] = fInputEvent->GetCentrality()->GetCentralityPercentile("V0A");
    fCentrality[2] = fInputEvent->GetCentrality()->GetCentralityPercentile("V0C");
    fCentrality[3] = fInputEvent->GetCentrality()->GetCentralityPercentile("TKL");
    fCentrality[4] = fInputEvent->GetCentrality()->GetCentralityPercentile("ZNA");
    fCentrality[5] = fInputEvent->GetCentrality()->GetCentralityPercentile("ZNC");
    //pileup
    fIsPileupSPD   = fUtils->IsPileUpSPD(fInputEvent);
    fIsPileupMV    = fUtils->IsPileUpMV(fInputEvent);

    //out-of-bunch pileup
    fIR1 = fInputEvent->GetHeader()->GetIRInt1InteractionMap();// bit=90 correspond to current interaction IR1 contains V0 information (VIR)
    fIR2 = fInputEvent->GetHeader()->GetIRInt2InteractionMap();// bit=90 correspond to current interaction IR2 contains T0 information

    for (Int_t i=0;i<6;i++) fNofITSClusters[i] = fInputEvent->GetNumberOfITSClusters(i);
    //fMultTKL = fInputEvent->GetMultiplicity()->GetNumberOfTracklets();
    //fMultV0Aeq=0;   for (Int_t i=32;i<64;i++) fMultV0Aeq+=fInputEvent->GetVZEROEqMultiplicity(i);
    //fMultV0Ceq=0;   for (Int_t i=0; i<32;i++) fMultV0Ceq+=fInputEvent->GetVZEROEqMultiplicity(i);
    /*fMultPercV0Aeq = fUtils->GetMultiplicityPercentile(fInputEvent,"V0AEq");
    fMultPercV0Ceq = fUtils->GetMultiplicityPercentile(fInputEvent,"V0CEq");
    fMultPercV0Meq = fUtils->GetMultiplicityPercentile(fInputEvent,"V0MEq");
    fMultPercV0B = fUtils->GetMultiplicityPercentile(fInputEvent,"V0B");
    fMultPercV0Apartial = fUtils->GetMultiplicityPercentile(fInputEvent,"V0Apartial");
    fMultPercV0Cpartial = fUtils->GetMultiplicityPercentile(fInputEvent,"V0Cpartial");
    fMultPercV0S = fUtils->GetMultiplicityPercentile(fInputEvent,"V0S");
    fMultPercV0SB = fUtils->GetMultiplicityPercentile(fInputEvent,"V0SB");*/

    //new multiplicity calibration class -- AliMultSelection
    fMultV0Aeq = -1.;
    fMultV0Ceq = -1.;
    fMultV0Meq = -1.;
    fMultTKL = -1.;
    fMultPercV0Aeq = -1.;
    fMultPercV0Ceq = -1.;
    fMultPercV0Meq = -1.;
    fMultPercTKL = -1.;
    fMultMeanV0M = -1.;
    fMultMeanV0A = -1.;
    fMultMeanV0C = -1.;
    fMultMeanTKL = -1.;
    fIsEventSel = -1;
    AliMultSelection *MultSelection = (AliMultSelection*) fInputEvent-> FindListObject("MultSelection");
    if(MultSelection){
      fMultV0Aeq =  MultSelection->GetEstimator("V0A")->GetValue();
      fMultV0Ceq =  MultSelection->GetEstimator("V0C")->GetValue();
      fMultV0Meq =  MultSelection->GetEstimator("V0M")->GetValue();
      fMultTKL = MultSelection->GetEstimator("SPDTracklets")->GetValue();
      fMultPercV0Aeq = MultSelection->GetMultiplicityPercentile("V0A");
      fMultPercV0Ceq = MultSelection->GetMultiplicityPercentile("V0C");
      fMultPercV0Meq = MultSelection->GetMultiplicityPercentile("V0M");
      fMultPercTKL = MultSelection->GetMultiplicityPercentile("SPDTracklets");
      fMultMeanV0A =  MultSelection->GetEstimator("V0A")->GetMean();
      fMultMeanV0C =  MultSelection->GetEstimator("V0C")->GetMean();
      fMultMeanV0M =  MultSelection->GetEstimator("V0M")->GetMean();
      fMultMeanTKL =  MultSelection->GetEstimator("SPDTracklets")->GetMean();
      fIsEventSel = MultSelection->GetEvSelCode();
      fEventStatistics->Fill("MultSelection found",1);
    }
    else{
      //Printf("Didn't find MultSelection in run %i",fCurrentRunNumber);
      fMultTKL = fInputEvent->GetMultiplicity()->GetNumberOfTracklets();
      fMultV0Aeq=0;   for (Int_t i=32;i<64;i++) fMultV0Aeq+=fInputEvent->GetVZEROEqMultiplicity(i);
      fMultV0Ceq=0;   for (Int_t i=0; i<32;i++) fMultV0Ceq+=fInputEvent->GetVZEROEqMultiplicity(i);
    }

    //online/offline spd
    AliVMultiplicity* mult = fInputEvent->GetMultiplicity();
    fonlineSPD     = mult->GetFastOrFiredChipMap().CountBits(400);
    fofflineSPD    = mult->GetFiredChipMap().CountBits(400);

    for (Int_t iring=0;iring<4;iring++){
      fMultV0Aring[iring]=0;
      fMultV0Cring[iring]=0;
      for (Int_t isector=0;isector<8;isector++){
        fMultV0Aring[iring]+=vzero->GetMultiplicityV0A(iring*8+isector);
        fMultV0Cring[iring]+=vzero->GetMultiplicityV0C(iring*8+isector);
      }
    }
    for (Int_t i=0;i<64;i++) fV0chan->Fill(i,fInputEvent->GetVZEROEqMultiplicity(i));



    const AliVVertex* vertex  = fInputEvent->GetPrimaryVertex();
    fVtxX  = vertex->GetX();
    fVtxY  = vertex->GetY();
    fVtxZ  = vertex->GetZ();
    fVtxTPConly = TString(vertex->GetName()).CompareTo("PrimaryVertex") && TString(vertex->GetName()).CompareTo("SPDVertex");
    fVtxContributors = vertex->GetNContributors();
    if (TMath::Abs(fVtxZ) >= fZVertexCut)  return;
    fEventStatistics->Fill("after vertex cut",1);

    if(fMCEvent){//get process type
      fIsDiffractive=0;
      AliGenPythiaEventHeader * headPy  = 0;
      AliGenDPMjetEventHeader * headPho = 0;
      AliGenEventHeader * htmp = fMCEvent->GenEventHeader();
      if(!htmp) {
        AliError("Cannot Get MC Header!!");
        return;
      }
      if( TString(htmp->IsA()->GetName()) == "AliGenPythiaEventHeader") {
        headPy =  (AliGenPythiaEventHeader*) htmp;
      } else if (TString(htmp->IsA()->GetName()) == "AliGenDPMjetEventHeader") {
        headPho = (AliGenDPMjetEventHeader*) htmp;
      } else {
        AliWarning("Unknown header");
      }
      //  Check if the event is diffractive
      if(headPy){
        //        Printf("Process: %d", headPy->ProcessType());
        if(headPy->ProcessType() == 92 || headPy->ProcessType() == 93)fIsDiffractive=1; //single diffractive
        else if(headPy->ProcessType() == 94)fIsDiffractive=2;//double diffractive
      }
      if (headPho){
        if(headPho->ProcessType() == 5 || headPho->ProcessType() == 6 )fIsDiffractive=1; //single difractive
        else if(headPho->ProcessType() != 7)fIsDiffractive=2;//double diffractive
      }
    }

    AliAODEvent* aodvtx = dynamic_cast<AliAODEvent*> (fInputEvent);
    if(aodvtx){
      AliVVertex* spdVtx = aodvtx->GetPrimaryVertexSPD();
      fSPDVtxZ = spdVtx->GetZ();
      //TString vtxTyp = spdVtx->GetTitle();
      Double_t cov[6]={0};
      spdVtx->GetCovarianceMatrix(cov);
      fSPDVtxZRes = (Float_t)TMath::Sqrt(cov[5]);
      fSPDVtxContributors = spdVtx->GetNContributors();
    }

    fNchTPC = -1;
    fNchTPCmc = -1;
    fNchV0Amc = -1;
    fNchV0Cmc = -1;
    fNchCL1mc = -1;

    if (fTracks) {
      fTracks->Clear();
      Int_t countNch = 0;
      Int_t countNch09 = 0;
      Int_t countNch14 = 0;
      for (Int_t ipart=0;ipart<fInputEvent->GetNumberOfTracks();ipart++){
        AliVTrack* track = (AliVTrack*) fInputEvent->GetTrack(ipart);
        if (!track) continue;
        UInt_t mask = GetFilterMap(track);
        if ( track->InheritsFrom("AliAODTrack") ) {
	  fitssatrackcuts->ExtractAndSetPrimaryVertex(fInputEvent);
	}
        if(fTrackFilterBit==2){//we are saving only ITSsa, It's a muon_calo_production!
          if((mask & (1<<1)) && track->Pt()>0.2 && TMath::Abs(track->Eta())<0.9) countNch++; // with pt > 0.2 and within |eta|<0.9
        }else if((mask & (1<<0)) && track->Pt()>0.2 && TMath::Abs(track->Eta())<0.9) countNch++; // with pt > 0.2 and within |eta|<0.9

	if ( !(mask==0 && fTrackFilterBit==0) ){ //this would apply the fbit selection only to tracks with mask defined
	  if(!(mask & fTrackFilterBit))continue; //filter bit cut
	}
        else if ( track->InheritsFrom("AliAODTrack") && !(fitssatrackcuts->AcceptTrack((AliAODTrack*)track)) ) continue;

        if (track->InheritsFrom("AliAODTrack")) AddTrack(track,mask,0);
        else if (track->InheritsFrom("AliESDtrack")) {
          if (mask)                           AddTrack(track,mask,1);
          if (mask & fHybridConstrainedMask)  AddTrack(track,mask,2);
          if (mask & fTPConlyConstrainedMask) AddTrack(track,mask,3);
        }
      }
      fNchTPC = countNch;
    }

    if (fTracklets){
      fTracklets->Clear();
      AliVMultiplicity* mult = fInputEvent->GetMultiplicity();
      Int_t nTracklets = mult->GetNumberOfTracklets();
      for (Int_t i=0;i<nTracklets;i++){
        Float_t phi   = mult->GetPhi(i);
        Float_t eta   = -TMath::Log(TMath::Tan(mult->GetTheta(i)/2));
        Float_t dphi  = mult->GetDeltaPhi(i);
        if (TMath::Abs(dphi)>fDphiCut) continue;
        if (TMath::Abs(eta)>fTrackletEtaCut) continue;

        AliCFParticle* tracklet = new ((*fTracklets)[fTracklets->GetEntriesFast()]) AliCFParticle(dphi,eta,phi,0,0,fStoreMcTracklets?4:0);
        if (!fStoreMcTracklets || !fMCEvent) continue;
        Int_t label1 = mult->GetLabel(i,0);
        Int_t label2 = mult->GetLabel(i,1);
        if (label1!=label2){
          tracklet->SetMask(999);//999 means fake tracklet
          continue;
        }
        AliVParticle* particle = fMCEvent->GetTrack(label1);
        if (!particle) continue;
        Short_t charge  = particle->Charge();
        Float_t ptMC    = particle->Pt();
        Float_t etaMC   = particle->Eta();
        Float_t phiMC   = particle->Phi();
        Float_t pdg     = particle->PdgCode();
        Bool_t primary  = particle->InheritsFrom("AliAODMCParticle") ? ((AliAODMCParticle*) particle)->IsPhysicalPrimary() : fMCEvent->IsPhysicalPrimary(label1);
        tracklet->SetCharge(charge);
        tracklet->SetMask(primary);
        tracklet->SetAt(ptMC,0);
        tracklet->SetAt(etaMC,1);
        tracklet->SetAt(phiMC,2);
        tracklet->SetAt(pdg,3);
      }
    }

    AliAODEvent* aod = dynamic_cast<AliAODEvent*> (fInputEvent);
    if (fMuons && aod){ // aod only
      fMuons->Clear();
      if (fStoreMuonOrigin) fMuonOrigin->Clear();
      for (Int_t iTrack = 0; iTrack < aod->GetNumberOfTracks(); iTrack++) {
        AliAODTrack* track = (AliAODTrack*) aod->GetTrack(iTrack);
        if (!track->IsMuonTrack()) continue;
        Float_t pt     = track->Pt();
        Float_t eta    = track->Eta();
        Float_t phi    = track->Phi();
        Short_t charge = track->Charge();
        Float_t dca    = track->DCA();
        Float_t chi2   = track->Chi2perNDF();
        Float_t rabs   = track->GetRAtAbsorberEnd();
        Int_t   mask   = track->GetMatchTrigger();
        Bool_t pdca    = fMuonTrackCuts->IsSelected(track);
        if (rabs < 17.6 || rabs > 89.5) continue;
        if (eta < -4 || eta > -2.5) continue;
        AliCFParticle* part = new ((*fMuons)[fMuons->GetEntriesFast()]) AliCFParticle(pt,eta,phi,charge,mask,fStoreMcMuons?16:4);
        part->SetAt(dca,0);
        part->SetAt(chi2,1);
        part->SetAt(rabs,2);
        part->SetAt(pdca,3);
        if (!fStoreMcMuons || !fMCEvent) continue;
        Int_t label = TMath::Abs(track->GetLabel());
        AliVParticle* mcpart = fMCEvent->GetTrack(label);
        if (!mcpart) continue;
        Float_t mcpt  = mcpart->Pt();
        Float_t mceta = mcpart->Eta();
        Float_t mcphi = mcpart->Phi();
        Int_t   mcpdg = mcpart->PdgCode();
        part->SetAt(mcpt,4);
        part->SetAt(mceta,5);
        part->SetAt(mcphi,6);
        part->SetAt(mcpdg,7);
        TString origin = Form("%+i",mcpdg);

        Bool_t isPrimary = ((AliAODMCParticle*) mcpart)->IsPhysicalPrimary();
        label = mcpart->GetMother();

        while (!isPrimary && label>=0) {
          mcpart = (AliVParticle*) fMCEvent->GetTrack(label);
          label = mcpart->GetMother();
          isPrimary = ((AliAODMCParticle*) mcpart)->IsPhysicalPrimary();
          origin+=Form("<%+i",mcpart->PdgCode());
        }

        Float_t mcprimarypdg = mcpart->PdgCode();
        Float_t mcprimarypt  = mcpart->Pt();
        Float_t mcprimaryeta = mcpart->Eta();
        Float_t mcprimaryphi = mcpart->Phi();
        part->SetAt(mcprimarypt,  8);
        part->SetAt(mcprimaryeta, 9);
        part->SetAt(mcprimaryphi,10);
        part->SetAt(mcprimarypdg,11);

        while (label>=0 && label < fMCEvent->GetNumberOfPrimaries()) {
          mcpart = (AliVParticle*) fMCEvent->GetTrack(label);
          label = mcpart->GetMother();
          origin+=Form("<%+i",mcpart->PdgCode());
        }

        Float_t mcoriginpdg = mcpart->PdgCode();
        Float_t mcoriginpt  = mcpart->Pt();
        Float_t mcorigineta = mcpart->Eta();
        Float_t mcoriginphi = mcpart->Phi();
        part->SetAt(mcoriginpt, 12);
        part->SetAt(mcorigineta,13);
        part->SetAt(mcoriginphi,14);
        part->SetAt(mcoriginpdg,15);
        if (fStoreMuonOrigin) new ((*fMuonOrigin)[fMuonOrigin->GetEntriesFast()]) TObjString(origin);
      }
    }
  }



  if(fMcMuons) fMcMuons->Clear();
  if(fMcParticles) fMcParticles->Clear();

  if (fMCEvent) {
    Int_t countNchMc = 0, countNchV0AMc = 0, countNchV0CMc = 0, countNchCL1Mc = 0;
    TLorentzVector h;

    for (Int_t ipart=0;ipart<fMCEvent->GetNumberOfTracks();ipart++){
      AliVParticle* mcpart = fMCEvent->GetTrack(ipart);
      Bool_t isPrimary = mcpart->InheritsFrom("AliAODMCParticle") ? ((AliAODMCParticle*) mcpart)->IsPhysicalPrimary() : fMCEvent->IsPhysicalPrimary(ipart);
      Float_t pt     = mcpart->Pt();
      Float_t eta    = mcpart->Eta();
      Float_t phi    = mcpart->Phi();
      Char_t  charge = mcpart->Charge();
      Int_t   pdg    = mcpart->PdgCode();
      Int_t pdgabs   = TMath::Abs(pdg);

      if(charge != 0 && isPrimary && pt>0.2 && TMath::Abs(eta)<1.) countNchMc++;
      if(charge != 0 && isPrimary && pt>0.001 && pt<50 && TMath::Abs(eta)<1.4) countNchCL1Mc++;
      if(charge != 0 && isPrimary && pt>0.001 && pt<50 && eta>2.8 && eta<5.1) countNchV0AMc++;
      if(charge != 0 && isPrimary && pt>0.001 && pt<50 && eta>-3.7 && eta<-1.7) countNchV0CMc++;

      // keep track of which generator-level tracks to write out (so that they don't get written out multiple times for different reasons
      Bool_t writeOutMC = kFALSE; 

      // if we are doing a track analysis: write out generator-level tracks needed for efficiency calculation and MC closure test
      if(fMcParticles) {

        // write out all MC tracks in fiducial region (same range as tracks)
        if(pt >= fPtMin && TMath::Abs(eta) <= fTrackEtaCut) writeOutMC = kTRUE;

        // write out all generated tracks corresponding to reconstructed tracks, regardless of their pt and eta
        if(!writeOutMC && fTracks)
        {
          Int_t additionalParameterIndex = 0;
          if (fStoreTrackInfo) additionalParameterIndex+=10;
          if (fStorePidInfo)   additionalParameterIndex+=3;
          if (fStoreMcTracks)  additionalParameterIndex+=1;
          additionalParameterIndex -= 1;
          for(Int_t jpart = 0; jpart<fTracks->GetEntriesFast(); jpart++)
          {
            Int_t genIndexTest = (Int_t)(((AliCFParticle*)fTracks->UncheckedAt(jpart))->GetAt(additionalParameterIndex));
            if(ipart==TMath::Abs(genIndexTest)){writeOutMC=kTRUE; break;}
          }
        }
        //if(!writeOutMC) continue;

      }

      // if we are doing a muon analysis: write out all primary and decay muons within the muon arm acceptance
      if(fMcMuons) {
        if (isPrimary) {
          if (pdgabs==211 || pdgabs==321 || pdgabs==310 || pdgabs==130 || pdgabs==13) {
            h.SetPxPyPzE(mcpart->Px(),mcpart->Py(),mcpart->Pz(),mcpart->E());
            fDecayer->Decay(pdg,&h);
            fDecayer->ImportParticles(fDecayArray);
            for (Int_t k=0;k<fDecayArray->GetEntriesFast();k++){
              TParticle* mu = (TParticle*) fDecayArray->At(k);
              if (TMath::Abs(mu->GetPdgCode())==13) {
                Float_t mu_eta = mu->Eta();
                Float_t mu_pt  = mu->Pt();
                Float_t mu_phi = mu->Phi();
                Char_t mu_charge = mu->GetPdgCode()>0 ? -1: 1;
                if ((mu_eta >= -4.0 && mu_eta <= -2.5) || (mu_eta >= 2.5 && mu_eta <= 4.0)) {
                  if (mu_pt >= 0 && mu_pt <= 5) {
                    AliCFParticle* mupart = new ((*fMcMuons)[fMcMuons->GetEntriesFast()]) AliCFParticle(mu_pt,mu_eta,mu_phi,mu_charge,pdg);
                  }
                }
              }
            }
          }
          //****** if you want to write out any track information (in fMcParticles) then it can be done here by turning on the writeOutMc flag
          //if (charge==0) continue;
          //if (!(eta>-5 && eta<5)) continue;
          //      if ((eta>-3.7 && eta<-2.7) || (eta>2.8 && eta<3.9)) fMultGenV0S++;
          //AliCFParticle* part = new ((*fMcParticles)[fMcParticles->GetEntriesFast()]) AliCFParticle(pt,eta,phi,charge,pdg);
        }
      }

      if(writeOutMC) {
        Int_t nAdditionalParameters = 4;
        AliCFParticle* cftrack = new ((*fMcParticles)[fMcParticles->GetEntriesFast()]) AliCFParticle(pt,eta,phi,charge,pdg,nAdditionalParameters);
        cftrack->SetAt(ipart,0); // index in generated MC array (matched with label of reconstructed tracks)
        cftrack->SetAt(((AliAODMCParticle*)mcpart)->GetLabel(),1);
        cftrack->SetAt(((AliAODMCParticle*)mcpart)->IsPhysicalPrimary(),2);
        Int_t motherlabel = (Int_t)(((AliAODMCParticle*)mcpart)->GetMother());
        AliVParticle* mothertrack = 0;
        if(motherlabel>=0)
          mothertrack = (AliVParticle*)fMCEvent->GetTrack(motherlabel);
        if(mothertrack) cftrack->SetAt(mothertrack->PdgCode(),3);
        else cftrack->SetAt(-999,3);
      }

    }

    fNchTPCmc = countNchMc;
    fNchV0Amc = countNchV0AMc;
    fNchV0Cmc = countNchV0CMc;
    fNchCL1mc = countNchCL1Mc;
  }

  if (!fStoreOnlyEventsWithMuons) fTree->Fill();
  else { if (fMuons) if (fMuons->GetEntriesFast()>0) fTree->Fill(); }

  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
UInt_t AliAnalysisTaskCFTree::GetFilterMap(AliVTrack* track){
  UInt_t mask = 0;
  if (track->InheritsFrom("AliAODTrack")) {
    AliAODTrack* part = (AliAODTrack*) track;
    mask = part->GetFilterMap();
    if (fStoreCutBitsInTrackMask) {
      Double_t nCrossedRaws      = part->GetTPCNCrossedRows();
      Double_t nFindableClusters = part->GetTPCNclsF();
      Double_t nSharedClusters   = part->GetTPCnclsS();
      Double_t nClusters         = part->GetTPCncls();
      Bool_t itsRefit            = part->GetStatus() & AliVTrack::kITSrefit;
      if (nCrossedRaws/nFindableClusters > fFoundFractionCut) mask |= (1 << 26);
      if (nCrossedRaws>fCrossedRowsCut)                       mask |= (1 << 27);
      if (itsRefit)                                           mask |= (1 << 28);
      if (nSharedClusters/nClusters<=fSharedClusterCut)       mask |= (1 << 29);
      if (part->GetLabel()<0)                                 mask |= (1 << 31);
    }
  } else if (track->InheritsFrom("AliESDtrack")){
    AliESDtrack* part = (AliESDtrack*) track;
    if (!fTrackFilter) AliFatal("Track filter undefined");
    mask |= fTrackFilter->IsSelected(part);
  }

  return mask;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
AliCFParticle* AliAnalysisTaskCFTree::AddTrack(AliVTrack* track, UInt_t mask, UInt_t flag){

  // skip neutral mc trackicles
  Char_t charge = track->Charge();
  if (charge==0) return NULL;

  // set pt,eta,phi
  Float_t pt=0,eta=0,phi=0;
  if (flag==0 || flag==1){ // AOD or Global ESD tracks
    pt  = track->Pt();
    eta = track->Eta();
    phi = track->Phi();
    if  (flag==1) mask &= (~fHybridConstrainedMask) & (~fTPConlyConstrainedMask);
  }
  else if (flag==2) { // Hybrid constrained tracks (ESD)
    AliESDtrack* part = (AliESDtrack*) track;
    const AliExternalTrackParam* param = part->GetConstrainedParam();
    pt  = param->Pt();
    eta = param->Eta();
    phi = param->Phi();
    mask &= fHybridConstrainedMask;
  }
  else if (flag==3) { // TPC only constrained tracks (ESD)
    AliESDtrack* part = (AliESDtrack*) track;
    AliESDtrack tpcTrack;
    if (!part->FillTPCOnlyTrack(tpcTrack)) return NULL;
    AliExternalTrackParam param;
    const AliESDVertex* vtxSPD = ((AliESDEvent*) fInputEvent)->GetPrimaryVertexSPD();
    if (!tpcTrack.RelateToVertexTPC(vtxSPD,fField,1.e30,&param)) return NULL;
    pt  = param.Pt();
    eta = param.Eta();
    phi = param.Phi();
    mask &= fTPConlyConstrainedMask;
  }

  // kinematic cuts
  if (pt < fPtMin || TMath::Abs(eta) > fTrackEtaCut) return NULL;

  //cut on SPD hits (for filter bit zero)
  //if(!(((AliAODTrack*)track)->HasPointOnITSLayer(0) || ((AliAODTrack*)track)->HasPointOnITSLayer(1))) return NULL;

  Int_t nAdditionalParameters = 0;
  if (fStoreTrackInfo) nAdditionalParameters+=10;
  if (fStorePidInfo)   nAdditionalParameters+=3;
  if (fStoreMcTracks)  nAdditionalParameters+=1;
  AliCFParticle* cftrack = new ((*fTracks)[fTracks->GetEntriesFast()]) AliCFParticle(pt,eta,phi,charge,mask,nAdditionalParameters);
  Int_t iParameter=0;
  if (fStoreTrackInfo && flag==0){
    AliAODTrack* part = (AliAODTrack*) track;
    cftrack->SetBit(AliAODTrack::kIsDCA              ,part->TestBit(AliAODTrack::kIsDCA));
    cftrack->SetBit(AliAODTrack::kUsedForVtxFit      ,part->TestBit(AliAODTrack::kUsedForVtxFit));
    cftrack->SetBit(AliAODTrack::kUsedForPrimVtxFit  ,part->TestBit(AliAODTrack::kUsedForPrimVtxFit));
    cftrack->SetBit(AliAODTrack::kIsTPCConstrained   ,part->TestBit(AliAODTrack::kIsTPCConstrained));
    cftrack->SetBit(AliAODTrack::kIsHybridTPCCG      ,part->TestBit(AliAODTrack::kIsHybridTPCCG));
    cftrack->SetBit(AliAODTrack::kIsGlobalConstrained,part->TestBit(AliAODTrack::kIsGlobalConstrained));
    cftrack->SetBit(AliAODTrack::kIsHybridGCG        ,part->TestBit(AliAODTrack::kIsHybridGCG));
    ULong_t status = part->GetStatus();
    cftrack->SetAt(*((Float_t*) &status),iParameter++);
    cftrack->SetAt(part->GetITSClusterMap(),iParameter++);
    cftrack->SetAt(part->GetTPCNCrossedRows(),iParameter++);
    cftrack->SetAt(part->GetTPCNclsF(),iParameter++);
    cftrack->SetAt(part->GetTPCnclsS(),iParameter++);
    cftrack->SetAt(part->GetTPCncls(),iParameter++);
    cftrack->SetAt(part->Chi2perNDF(),iParameter++);
    Double_t xyz[3];
    part->GetXYZ(xyz);
    cftrack->SetAt(xyz[0],iParameter++);
    cftrack->SetAt(xyz[1],iParameter++);
    cftrack->SetAt(xyz[2],iParameter++);
  }

  if (fStorePidInfo){
    Float_t ncl  = track->GetTPCsignalN();
    Float_t dedx = track->GetTPCsignalTunedOnData(); if (dedx<=0) dedx = track->GetTPCsignal();
    Float_t beta = -99;
    if (track->GetStatus()&AliESDtrack::kTOFpid){
      Double_t tof[5];
      track->GetIntegratedTimes(tof);
      beta = tof[0]/track->GetTOFsignal();
    }
    cftrack->SetAt(ncl,iParameter++);
    cftrack->SetAt(dedx,iParameter++);
    cftrack->SetAt(beta,iParameter++);
  }

  if(fStoreMcTracks)
  {
    cftrack->SetAt(track->GetLabel(),iParameter++);
  }

  if(fStoreAodDCAInfo)
  {
    cftrack->SetAt(fitssatrackcuts->CalculateDCAXY((AliAODTrack*)track),iParameter++);
    cftrack->SetAt(fitssatrackcuts->CalculateDCAZ((AliAODTrack*)track),iParameter++);
  }

  return cftrack;
}
