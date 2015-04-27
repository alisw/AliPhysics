#include "AliAnalysisTaskTPCTOFPID.h"
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
#include "AliAnalysisPIDTrack.h"
#include "AliAnalysisPIDParticle.h"
#include "AliAnalysisPIDEvent.h"
#include "TClonesArray.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "TTree.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliPIDResponse.h"
#include "AliPPVsMultUtils.h"
#include "TRandom.h"


ClassImp(AliAnalysisTaskTPCTOFPID)
  
//_______________________________________________________
  
AliAnalysisTaskTPCTOFPID::AliAnalysisTaskTPCTOFPID() :
  AliAnalysisTaskSE("AnalysisResults"),
  fESDtrackCuts(0),
  fInitFlag(kFALSE),
  fMCFlag(kFALSE),
  fMCTuneFlag(kFALSE),
  fPbPbFlag(kFALSE),
  fVertexSelectionFlag(kFALSE),
  fPrimaryDCASelectionFlag(kFALSE),
  fPIDTree(0),
  fPIDResponse(0),
  fMultiUtils(0),
  fRunNumber(0),
  fStartTime(0),
  fEndTime(0),
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
  fAnalysisEvent(new AliAnalysisPIDEvent()),
  fAnalysisTrackArray(new TClonesArray("AliAnalysisPIDTrack")),
  fAnalysisTrack(new AliAnalysisPIDTrack()),
  fAnalysisParticleArray(new TClonesArray("AliAnalysisPIDParticle")),
  fAnalysisParticle(new AliAnalysisPIDParticle()),
  fTOFcalib(new AliTOFcalib()),
  fTOFT0maker(new AliTOFT0maker(fESDpid)),
  fTimeResolution(80.),
  fVertexCut(10.0),
  fRapidityCut(1.0),
  V0MBinCount(0),
  fHistoList(new TList()),
  fMCHistoList(new TList())
{
  /* 
   * default constructor 
   */
  fTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE); //If not running, set to kFALSE;
  fTrackCuts->SetEtaRange(-0.9,0.9);
  //V0MBinCount = new TH1F("V0MBinCount","V0MBinCount",12,V0MBinsDefault);
  //printf("Jegele!\n\n\n\n\n\n");

}




AliAnalysisTaskTPCTOFPID::AliAnalysisTaskTPCTOFPID(Bool_t isMC) :
  AliAnalysisTaskSE("AnalysisResults"),
  fESDtrackCuts(0),
  fInitFlag(kFALSE),
  fMCFlag(kFALSE),
  fMCTuneFlag(kFALSE),
  fPbPbFlag(kFALSE),
  fVertexSelectionFlag(kFALSE),
  fPrimaryDCASelectionFlag(kFALSE),
  fPIDTree(0),
  fPIDResponse(0),
  fMultiUtils(0),
  fRunNumber(0),
  fStartTime(0),
  fEndTime(0),
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
  fAnalysisEvent(new AliAnalysisPIDEvent()),
  fAnalysisTrackArray(new TClonesArray("AliAnalysisPIDTrack")),
  fAnalysisTrack(new AliAnalysisPIDTrack()),
  fAnalysisParticleArray(new TClonesArray("AliAnalysisPIDParticle")),
  fAnalysisParticle(new AliAnalysisPIDParticle()),
  fTOFcalib(new AliTOFcalib()),
  fTOFT0maker(new AliTOFT0maker(fESDpid)),
  fTimeResolution(80.),
  fVertexCut(10.0),
  fRapidityCut(1.0),
  V0MBinCount(0),
  fHistoList(new TList()),
  fMCHistoList(new TList())
{
  /* 
   * default constructor 
   */
  fTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE); //If not running, set to false
  fTrackCuts->SetEtaRange(-0.9,0.9);
  fMCFlag = isMC;
  DefineOutput(1, TTree::Class());
}













//_______________________________________________________

AliAnalysisTaskTPCTOFPID::~AliAnalysisTaskTPCTOFPID()
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


//________________________________________________________________________

void
AliAnalysisTaskTPCTOFPID::UserCreateOutputObjects()
{
  /*
   * user create output objects
   */

  /* output tree */
  fPIDTree = new TTree("PIDTree","PIDTree");
  fPIDTree->Branch("AnalysisEvent", "AliAnalysisPIDEvent", &fAnalysisEvent);  
  fPIDTree->Branch("AnalysisTrack", "TClonesArray", &fAnalysisTrackArray);  
  if (fMCFlag)
    fPIDTree->Branch("AnalysisParticle", "TClonesArray", &fAnalysisParticleArray);


  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse(); 
  fMultiUtils = new AliPPVsMultUtils();
  Double_t V0MbinsDefault[13] = {0, 0.01, 0.1, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
  V0MBinCount = new TH1F("V0MBinCount","V0MBinCount",12,V0MbinsDefault);
  PostData(1,fPIDTree);
}

//_______________________________________________________





Bool_t
AliAnalysisTaskTPCTOFPID::InitRun()
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
AliAnalysisTaskTPCTOFPID::InitEvent()
{
  /*
   * init event
   */

  /* get ESD event */
  fESDEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  //printf("Checking if fESD...");
  if (!fESDEvent) return kFALSE;
  //printf("pass!\n");
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
    TString vtxTyp = vertex->GetTitle();
    Double_t cov[6]={0};
    vertex->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) fHasVertex = kFALSE;
  }
  else fHasVertex = kTRUE;

  //  if (fVertexSelectionFlag && (!fHasVertex || TMath::Abs(vertex->GetZ()) > fVertexCut)) return kFALSE;
  
  fVertexZ = vertex->GetZ();
  /* centrality in PbPb */
  if (fPbPbFlag) {
    fCentrality = fESDEvent->GetCentrality();
  }
  //printf("smack");
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
AliAnalysisTaskTPCTOFPID::HasPrimaryDCA(AliESDtrack *track)
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
AliAnalysisTaskTPCTOFPID::UserExec(Option_t *option)
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
  //TRandom *rn = new TRandom(0);
  //Float_t ScaleFactor[] = {0.,0.,1.000000, 1.474651, 1.886733, 2.249043, 2.592145, 3.144156, 3.997966, 5.010369, 6.822343, 12.020458};//{0.,0.,1.000000, 1.279444, 1.525135, 1.757802, 2.132135, 2.711126, 3.397664, 4.626411, 8.151391};
  Float_t V0MPercentile = fMultiUtils->GetMultiplicityPercentile(fESDEvent, "V0M");
  /*if(V0MPercentile<1) V0MBinCount->Fill(V0MPercentile); else {
    if(rn->Rndm()<ScaleFactor[V0MBinCount->GetXaxis()->FindBin(V0MPercentile)]/V0MBinCount->GetXaxis()->GetBinWidth(V0MBinCount->GetXaxis()->FindBin(V0MPercentile))) V0MBinCount->Fill(V0MPercentile); else
      return;
      };*/
  //if(!(V0MPercentile>0)) return;


  /*if(fMultiUtils->GetMultiplicityPercentile(fESDEvent, "V0M")>5) {
    TRandom *rn = new TRandom(0);
    if(rn->Rndm()>0.05) return;
    };*/
  /* set fill AOD */
  ((AliAODHandler*)(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()))->SetFillAOD(kTRUE);

  /*** MC PRIMARY PARTICLES ***/

  Int_t mcmulti = 0;
  // if(fMultiUtils->GetMultiplicityPercentile(fESDEvent, "V0M")>1) {
  //   TRandom *rn = new TRandom(0);
  //   if(rn->Rndm()>0.01) return;
  // };
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
      new ((*fAnalysisParticleArray)[fAnalysisParticleArray->GetEntries()]) AliAnalysisPIDParticle(*fAnalysisParticle);
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
    for (Int_t icent = 0; icent < AliAnalysisPIDEvent::kNCentralityEstimators; icent++) 
      fAnalysisEvent->SetCentralityPercentile(icent, fCentrality->GetCentralityPercentileUnchecked(AliAnalysisPIDEvent::fgkCentralityEstimatorName[icent]));
  }
  Int_t refmulti[6];
  for(Int_t i=0;i<6;i++) refmulti[i] = AliESDtrackCuts::GetReferenceMultiplicity(fESDEvent, AliESDtrackCuts::kTrackletsITSTPC,(i+3)*0.1); 
  fAnalysisEvent->SetReferenceMultiplicity(refmulti);
  fAnalysisEvent->SetMCMultiplicity(mcmulti);
  fAnalysisEvent->SetV0Mmultiplicity(fMultiUtils->GetMultiplicityPercentile(fESDEvent, "V0M"));
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
    //if (!fTrackCuts->AcceptTrack(track)) continue;
    
    /* update and add analysis track */
    fAnalysisTrack->Update(track, fMCStack, fMCEvent,fPIDResponse, fTrackCuts->AcceptTrack(track));
    new ((*fAnalysisTrackArray)[fAnalysisTrackArray->GetEntries()]) AliAnalysisPIDTrack(*fAnalysisTrack);

  } /* end of loop over ESD tracks */
  fPIDTree->Fill();

  PostData(1,fPIDTree);
}

void AliAnalysisTaskTPCTOFPID::Terminate(Option_t *) {
    PostData(1,fPIDTree);
  printf("Terminate!\n");
}

