#include "AliAnalysisTaskTOFSpectraPbPb.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliPhysicsSelection.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliTOFcalib.h"
#include "AliTOFT0maker.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "AliCDBManager.h"
#include "AliLog.h"
#include "AliESDtrack.h"
#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "AliTOFRunParams.h"
#include "AliAnalysisTrack.h"
#include "AliAnalysisParticle.h"
#include "AliAnalysisEvent.h"
#include "TClonesArray.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "TTree.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"

ClassImp(AliAnalysisTaskTOFSpectraPbPb)
  
//_______________________________________________________
  
AliAnalysisTaskTOFSpectraPbPb::AliAnalysisTaskTOFSpectraPbPb() :
  AliAnalysisTaskSE("TOFSpectraPbPb"),
  fInitFlag(kFALSE),
  fMCFlag(kFALSE),
  fMCTuneFlag(kFALSE),
  fPbPbFlag(kFALSE),
  fVertexSelectionFlag(kFALSE),
  fRunNumber(0),
  fESDEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fTrackCuts(NULL),
  fESDpid(new AliESDpid()),
  fIsCollisionCandidate(kFALSE),
  fIsEventSelected(0),
  fIsPileupFromSPD(kFALSE),
  fHasVertex(kFALSE),
  fVertexZ(0.),
  fMCTimeZero(0.),
  fCentrality(NULL),
  fAnalysisEvent(new AliAnalysisEvent()),
  fAnalysisTrackArray(new TClonesArray("AliAnalysisTrack")),
  fAnalysisTrack(new AliAnalysisTrack()),
  fAnalysisParticleArray(new TClonesArray("AliAnalysisParticle")),
  fAnalysisParticle(new AliAnalysisParticle()),
  fTOFcalib(new AliTOFcalib()),
  fTOFT0maker(new AliTOFT0maker(fESDpid)),
  fTimeResolution(80.),
  fVertexCut(10.0),
  fRapidityCut(1.0),
  fHistoList(new TList()),
  fMCHistoList(new TList())
{
  /* 
   * default constructor 
   */

}

//_______________________________________________________

AliAnalysisTaskTOFSpectraPbPb::~AliAnalysisTaskTOFSpectraPbPb()
{
  /*
   * default destructor
   */

  if (fTrackCuts) delete fTrackCuts;
  delete fESDpid;
  delete fTOFcalib;
  delete fTOFT0maker;
  delete fHistoList;
  delete fMCHistoList;
}

//_______________________________________________________

void
AliAnalysisTaskTOFSpectraPbPb::UserCreateOutputObjects()
{
  /*
   * user create output objects
   */

  /* output tree */
  OutputTree()->Branch("AnalysisEvent", "AliAnalysisEvent", &fAnalysisEvent);  
  OutputTree()->Branch("AnalysisTrack", "TClonesArray", &fAnalysisTrackArray);  
  if (fMCFlag)
    OutputTree()->Branch("AnalysisParticle", "TClonesArray", &fAnalysisParticleArray);  
}

//_______________________________________________________

Bool_t
AliAnalysisTaskTOFSpectraPbPb::InitRun()
{
  /*
   * init run
   */

  /* get ESD event */
  fESDEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESDEvent) {
    AliError("cannot get ESD event");
    return kFALSE;
  }
  /* get run number */
  Int_t runNb = fESDEvent->GetRunNumber();
  /* check run already initialized */
  if (fInitFlag && fRunNumber == runNb) return kTRUE;
  /* init cdb */
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(runNb);
  /* init TOF calib */
  if (!fTOFcalib->Init()) {
    AliError("cannot init TOF calib");
    return kFALSE;
  }
  AliInfo(Form("initialized for run %d", runNb));
  fInitFlag = kTRUE;
  fRunNumber = runNb;
  return kTRUE;
}

//_______________________________________________________

Bool_t
AliAnalysisTaskTOFSpectraPbPb::InitEvent()
{
  /*
   * init event
   */

  /* get ESD event */
  fESDEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESDEvent) return kFALSE;
  /* get MC event */
  if (fMCFlag) {
    fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fMCEvent) return kFALSE;
  }
  /* get stack */
  if (fMCFlag) {
    fMCStack = fMCEvent->Stack();
    if (!fMCStack) return kFALSE;
  }
  /* event selection */
  fIsCollisionCandidate = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAny);
  fIsEventSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  fIsPileupFromSPD = fESDEvent->IsPileupFromSPD();

  /* vertex selection */
  const AliESDVertex *vertex = fESDEvent->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() < 1) {
    vertex = fESDEvent->GetPrimaryVertexSPD();
    if (vertex->GetNContributors() < 1) fHasVertex = kFALSE;
    else fHasVertex = kTRUE;
    Double_t cov[6]={0};
    vertex->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (vertex->IsFromVertexerZ() && (zRes>0.25)) fHasVertex = kFALSE;
  }
  else fHasVertex = kTRUE;

  //  if (fVertexSelectionFlag && (!fHasVertex || TMath::Abs(vertex->GetZ()) > fVertexCut)) return kFALSE;
  
  fVertexZ = vertex->GetZ();
  /* centrality in PbPb */
  if (fPbPbFlag) {
    fCentrality = fESDEvent->GetCentrality();
  }
  printf("smack");
  /* calibrate ESD (also in MC to correct TExp) */
  fTOFcalib->CalibrateESD(fESDEvent);

  /* check MC flags */
  if (fMCFlag && fMCTuneFlag)
    fMCTimeZero = fTOFcalib->TuneForMC(fESDEvent, fTimeResolution);

#if 0 /* NEW METHOD */
  /* compute TOF-T0, fill ESD and make PID */
  fTOFT0maker->ComputeT0TOF(fESDEvent);
  fTOFT0maker->WriteInESD(fESDEvent);
  fESDpid->SetTOFResponse(fESDEvent, AliESDpid::kTOF_T0);
  fESDpid->MakePID(fESDEvent, kFALSE, 0.);
#endif

#if 1 /* OLD METHOD */
  /* compute and apply TOF-T0 */
  fTOFT0maker->ComputeT0TOF(fESDEvent);
  fESDpid->MakePID(fESDEvent, kFALSE, 0.);
#endif

  return kTRUE;
}

//_______________________________________________________

Bool_t
AliAnalysisTaskTOFSpectraPbPb::HasPrimaryDCA(AliESDtrack *track)
{
  /*
   * has primary DCA
   */
  
  // cut on transverse impact parameter
  Float_t d0z0[2],covd0z0[3];
  track->GetImpactParameters(d0z0, covd0z0);
  Float_t sigma= 0.0050 + 0.0060 / TMath::Power(track->Pt(), 0.9);
  Float_t d0max = 7. * sigma;
  //
  Float_t sigmaZ = 0.0146 + 0.0070 / TMath::Power(track->Pt(), 1.114758);
  if (track->Pt() > 1.) sigmaZ = 0.0216;
  Float_t d0maxZ = 5. * sigmaZ;
  //
  if(TMath::Abs(d0z0[0]) > d0max || TMath::Abs(d0z0[1]) > d0maxZ)
    return kFALSE;
  
  /* primary DCA ok */
  return kTRUE;
}

//_______________________________________________________

void
AliAnalysisTaskTOFSpectraPbPb::UserExec(Option_t *option)
{
  /*
   * user exec
   */

  /*** INITIALIZATION ***/

  /* unset fill AOD */
  ((AliAODHandler*)(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()))->SetFillAOD(kFALSE);

  /* init run */
  if (!InitRun()) return;
  /* init event */
  if (!InitEvent()) return;

  /* set fill AOD */
  ((AliAODHandler*)(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()))->SetFillAOD(kTRUE);

  /*** MC PRIMARY PARTICLES ***/

  Int_t mcmulti = 0;

  if (fMCFlag) {

    /* reset track array */
    fAnalysisParticleArray->Clear();
    
    /* loop over primary particles */
    Int_t nPrimaries = fMCStack->GetNprimary();
    TParticle *particle;
    TParticlePDG *particlePDG;
    /* loop over primary particles */
    for (Int_t ipart = 0; ipart < nPrimaries; ipart++) {
      /* check primary */
      if (!fMCStack->IsPhysicalPrimary(ipart)) continue;
      /* get particle */
      particle = fMCStack->Particle(ipart);
      if (!particle) continue;
      /* get particlePDG */
      particlePDG = particle->GetPDG();
      if (!particlePDG) continue;
      /* check charged */
      if (particlePDG->Charge() == 0.) continue;
      mcmulti++;
      /* check rapidity and pt cuts */
      if (TMath::Abs(particle->Y()) > fRapidityCut) continue;
      if (particle->Pt() < 0.15) continue;
      /* update and add analysis particle */
      fAnalysisParticle->Update(particle, ipart);
      new ((*fAnalysisParticleArray)[fAnalysisParticleArray->GetEntries()]) AliAnalysisParticle(*fAnalysisParticle);
    } /* end of loop over primary particles */

    
  }

  /*** GLOBAL EVENT INFORMATION ***/

  fAnalysisEvent->Reset();
  /* update global event info */
  fAnalysisEvent->SetIsCollisionCandidate(fIsCollisionCandidate);
  fAnalysisEvent->SetIsEventSelected(fIsEventSelected);
  fAnalysisEvent->SetIsPileupFromSPD(fIsPileupFromSPD);
  fAnalysisEvent->SetHasVertex(fHasVertex);
  fAnalysisEvent->SetVertexZ(fVertexZ);
  fAnalysisEvent->SetMCTimeZero(fMCTimeZero);
  /* update TOF event info */
  for (Int_t i = 0; i < 10; i++) {
    fAnalysisEvent->SetTimeZeroTOF(i, fESDpid->GetTOFResponse().GetT0bin(i));
    fAnalysisEvent->SetTimeZeroTOFSigma(i, fESDpid->GetTOFResponse().GetT0binRes(i));
  }
  /* update T0 event info */
  for (Int_t i = 0; i < 3; i++)
    fAnalysisEvent->SetTimeZeroT0(i, fESDEvent->GetT0TOF(i));

  /* update centrality info in PbPb */
  if (fPbPbFlag) {
    fAnalysisEvent->SetCentralityQuality(fCentrality->GetQuality());
    for (Int_t icent = 0; icent < AliAnalysisEvent::kNCentralityEstimators; icent++) 
      fAnalysisEvent->SetCentralityPercentile(icent, fCentrality->GetCentralityPercentileUnchecked(AliAnalysisEvent::fgkCentralityEstimatorName[icent]));
  }
  Int_t refmulti = AliESDtrackCuts::GetReferenceMultiplicity(fESDEvent, AliESDtrackCuts::kTrackletsITSTPC,0.5); 
  fAnalysisEvent->SetReferenceMultiplicity(refmulti);
  fAnalysisEvent->SetMCMultiplicity(mcmulti);

  /*** RECONSTRUCTED TRACKS ***/

  /* reset track array */
  fAnalysisTrackArray->Clear();

  /* loop over ESD tracks */
  Int_t nTracks = fESDEvent->GetNumberOfTracks();
  AliESDtrack *track;
  for (Int_t itrk = 0; itrk < nTracks; itrk++) {
    /* get track */
    track = fESDEvent->GetTrack(itrk);
    if (!track) continue;
    /* check accept track */
    if (!fTrackCuts->AcceptTrack(track)) continue;
    /* update and add analysis track */
    fAnalysisTrack->Update(track, fMCStack, fMCEvent);
    new ((*fAnalysisTrackArray)[fAnalysisTrackArray->GetEntries()]) AliAnalysisTrack(*fAnalysisTrack);

  } /* end of loop over ESD tracks */
  
}

