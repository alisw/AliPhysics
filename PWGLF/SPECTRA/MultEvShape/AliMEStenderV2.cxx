#include <TROOT.h>
#include <TString.h>
#include <TH1.h>
#include <THnSparse.h>
#include <TTreeStream.h>
#include <TGeoGlobalMagField.h>

#include <AliAnalysisManager.h>
#include <AliGeomManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDtrackCuts.h>
#include <AliPIDCombined.h>
#include <AliESDInputHandler.h>
#include <AliESDInputHandlerRP.h>
#include <AliMCEventHandler.h>

#include <AliESDEvent.h>
#include <AliESDHeader.h>
#include <AliESDRun.h>
#include <AliESDtrack.h>
#include <AliStack.h>
#include <AliMCEvent.h>
#include <AliMCParticle.h>

#include <AliMultiplicity.h>
#include <AliCentrality.h>
#include <AliPID.h>
#include <AliPIDResponse.h>
#include <AliAnalysisUtils.h>
#include <AliPPVsMultUtils.h>
#include <AliPhysicsSelectionTask.h>
#include "AliMultSelection.h"

#include "AliMESbaseTask.h"
#include "AliMEStenderV2.h"
#include "AliMESeventInfo.h"
#include "AliMEStrackInfo.h"
#include "TTree.h"

using namespace std;

ClassImp(AliMEStenderV2::AliMESconfigTender)
    ClassImp(AliMEStenderV2)
 
    //________________________________________________________________________
    AliMEStenderV2::AliMESconfigTender::AliMESconfigTender() : TObject(), fTrackCuts(0), fEventCuts(0), fPIDpriors(0)
{
}

//________________________________________________________________________
void AliMEStenderV2::AliMESconfigTender::Print(Option_t *) const
{
  // Dump config info to stdout
  printf("MES TENDER CONFIGURATION\n   Event cuts : ");
  switch (fEventCuts)
  {
  case kNoEC:
    printf("No\n");
    break;
  case k7TeV:
    printf("Trigger[MB, HM], Vertex[Yes]\n");
    break;
  case k13TeV:
    printf("13TeV: Trigger[MB, HM], Vertex[Yes]\n");
    break;
  default:
    printf("Not defined [%d]\n", fEventCuts);
    break;
  }
  printf("   Track cuts : ");
  switch (fTrackCuts)
  {
  case kNoTC:
    printf("No\n");
    break;
  case kStandardITSTPCTrackCuts2010:
    printf("StandardITSTPCTrackCuts2010\n");
    break;
  case kStandardITSTPCTrackCuts2011:
    printf("StandardITSTPCTrackCuts2011\n");
    break;
  default:
    printf("Not defined [%d]\n", fTrackCuts);
    break;
  }
  printf("   PID priors : ");
  switch (fPIDpriors)
  {
  case kNoPP:
    printf("Flat Priors\n");
    break;
  case kTPC:
    printf("DefaultTPCPriors\n");
    break;
  case kIterative:
    printf("LHC10d Iterative Priors\n");
    break;
  default:
    printf("Not defined [%d]\n", fPIDpriors);
    break;
  }

  printf("\n");
}

//________________________________________________________________________
Bool_t AliMEStenderV2::ConfigTask(AliMESconfigTender::EMESconfigEventCuts ec,
                                  AliMESconfigTender::EMESconfigTrackCuts tc,
                                  AliMESconfigTender::EMESconfigPIDpriors pp)
{
  // configure tender for current run
  fConfig.fEventCuts = ec;
  fConfig.fTrackCuts = tc;
  fConfig.fPIDpriors = pp;
  fConfig.Print();
  return kTRUE;
}

//________________________________________________________________________
AliMEStenderV2::AliMEStenderV2()
    : AliAnalysisTaskSE(), fConfig(), fTrackFilter(NULL), fPIDcomb(NULL), fTracks(NULL), fEvInfo(NULL), fMCtracks(NULL), fMCevInfo(NULL), fTreeSRedirector(NULL), fTree(NULL), fUtils(NULL), fEventCutsQA(NULL)
{
  //
  // Constructor
  //
}

//________________________________________________________________________
AliMEStenderV2::AliMEStenderV2(const char *name)
    : AliAnalysisTaskSE(name), fConfig(), fTrackFilter(NULL), fPIDcomb(NULL), fTracks(NULL), fEvInfo(NULL), fMCtracks(NULL), fMCevInfo(NULL), fTreeSRedirector(NULL), fTree(NULL), fUtils(NULL), fEventCutsQA(NULL)
{ //
  // Constructor
  //
  DefineOutput(kQA, TList::Class());
  DefineOutput(kEventInfo + 1, AliMESeventInfo::Class());
  DefineOutput(kTracks + 1, TObjArray::Class());
  DefineOutput(kTree + 1, TTree::Class());
}

//________________________________________________________________________
void AliMEStenderV2::SetMCdata(Bool_t mc)
{
  // prepare task for MC processing
  SetBit(kMCdata, mc);
  if (mc)
  {
    DefineOutput(kMCeventInfo + 1, AliMESeventInfo::Class());
    DefineOutput(kMCtracks + 1, TObjArray::Class());
  }
}
//________________________________________________________________________
AliMEStenderV2::~AliMEStenderV2()
{
  //
  // Destructor
  //
  if (fTrackFilter)
    delete fTrackFilter;
  if (fPIDcomb)
    delete fPIDcomb;

  if (fEvInfo)
    delete fEvInfo;
  if (fTracks->GetEntries())
  {
    fTracks->Delete();
    delete fTracks;
  }
  if (fMCevInfo)
    delete fMCevInfo;
  if (fMCtracks)
  {
    fMCtracks->Delete();
    delete fMCtracks;
  }

  if (fUtils)
    delete fUtils;
}

//________________________________________________________________________
void AliMEStenderV2::UserCreateOutputObjects()
{

  // Build user objects
  BuildQAHistos();
  PostData(kQA, fHistosQA);

  OpenFile(1);
  fTreeSRedirector = new TTreeSRedirector();

  // ------- track cuts
  fTrackFilter = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts *lTrackCuts(NULL);
  switch (fConfig.fTrackCuts)
  {
  case AliMESconfigTender::kStandardITSTPCTrackCuts2010:
    lTrackCuts = new AliESDtrackCuts("std10TC", "Standard 2010");
    lTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE, 0);
    fTrackFilter->AddCuts(lTrackCuts);
    break;
  case AliMESconfigTender::kStandardITSTPCTrackCuts2011:
    lTrackCuts = new AliESDtrackCuts("std11TC", "Standard 2011");
    lTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 0); // kTRUE for primaries
    fTrackFilter->AddCuts(lTrackCuts);
    break;
  case AliMESconfigTender::kNoTC:
  default:
    AliDebug(2, "No track cuts selected");
    break;
  }

  fTree = ((*fTreeSRedirector) << "data/rec").GetTree();

  fTracks = new TObjArray(200);
  fTracks->SetOwner(kTRUE);
  fEvInfo = new AliMESeventInfo;
  PostData(kEventInfo + 1, fEvInfo);
  PostData(kTracks + 1, fTracks);
  PostData(kTree + 1, fTree);

  fUtils = new AliPPVsMultUtils();

  if (!HasMCdata())
    return;
  fMCtracks = new TObjArray(200);
  fMCtracks->SetOwner(kTRUE);
  fMCevInfo = new AliMESeventInfo;

  PostData(kMCeventInfo + 1, fMCevInfo);
  PostData(kMCtracks + 1, fMCtracks);
}

#include "AliGRPManager.h"
//________________________________________________________________________
void AliMEStenderV2::UserExec(Option_t * /*opt*/)
{

  AliVEvent *ev = InputEvent();
  if (!ev)
  {
    Error("UserExec", "Could not retrieve event");
    return;
  }

  AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESD)
  {
    AliError("ESD event not available");
    return;
  }

  //   AliMESeventInfo *fEvInfo   = dynamic_cast<AliMESeventInfo*>(GetOutputData(AliMESbaseTask::kEventInfo+1));
  if (!fEvInfo)
  {
    AliError("REC event info missing. Processing skipped");
    return;
  }
  //   TObjArray *fTracks   = dynamic_cast<TObjArray*>(GetOutputData(AliMESbaseTask::kTracks+1));
  if (!fTracks)
  {
    AliError("REC track array missing. Processing skipped");
    return;
  }
  fEvInfo->Clear("");
  fTracks->Delete();
  // check QA container
  if (!fHistosQA)
  {
    AliError("No QA container defined.");
    return;
  }

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
  AliPIDResponse *pidResponse = inputHandler->GetPIDResponse();
  if (!pidResponse)
    AliFatal("This Task needs the PID response attached to the inputHandler");

  // init magnetic field
  if (!TGeoGlobalMagField::Instance()->GetField() && !TGeoGlobalMagField::Instance()->IsLocked())
  {
    AliGRPManager grpManager;
    if (!grpManager.ReadGRPEntry())
      AliError("Cannot get GRP entry for magnetic field");
    if (!grpManager.SetMagField())
      AliError("Problem with magnetic field setup");
  }

  Bool_t incompleteDAQ = kFALSE; // if true should be rejected
  AliAnalysisUtils SPDclustersTrackletsPileup;

  switch (fConfig.fEventCuts)
  {
  case AliMESconfigTender::k7TeV:

    ((TH1 *)fHistosQA->At(kEfficiency))->Fill(0); // all events
    if (AliPPVsMultUtils::IsMinimumBias(fESD))
    {
      ((TH1 *)fHistosQA->At(kEfficiency))->Fill(1); // events after Physics Selection (for MB normalisation to INEL)
    }
    if ((!AliPPVsMultUtils::IsEventSelected(fESD, AliVEvent::kMB)) && (!AliPPVsMultUtils::IsEventSelected(fESD, AliVEvent::kHighMult)))
    {
      return;
    }
    ((TH1 *)fHistosQA->At(kEfficiency))->Fill(2); // analyzed events
    break;

  case AliMESconfigTender::k13TeV:

    incompleteDAQ = fESD->IsIncompleteDAQ();

    ((TH1 *)fHistosQA->At(kEfficiency))->Fill(0); // all events
    if (AliPPVsMultUtils::IsMinimumBias(fESD))
    {
      ((TH1 *)fHistosQA->At(kEfficiency))->Fill(1); // events after Physics Selection (for MB normalisation to INEL)
    }

    if (!(inputHandler->IsEventSelected()) ||
        incompleteDAQ ||
        SPDclustersTrackletsPileup.IsSPDClusterVsTrackletBG(fESD))
    {
      return;
    }
    ((TH1 *)fHistosQA->At(kEfficiency))->Fill(2); // analyzed events
    break;
  default:
    AliDebug(2, "No event cuts selected");
  }

  // TRIGGER SELECTION
  // MB & HM triggers
  Bool_t triggerMB = 0;
  Bool_t triggerHM = 0;
  switch (fConfig.fEventCuts)
  {
  case AliMESconfigTender::k7TeV:
    triggerMB = (inputHandler->IsEventSelected() & AliVEvent::kMB),
    triggerHM = (inputHandler->IsEventSelected() & AliVEvent::kHighMult);
    break;
  case AliMESconfigTender::k13TeV:
    triggerMB = (inputHandler->IsEventSelected() & AliVEvent::kINT7),   // default - for MB
        triggerHM = (inputHandler->IsEventSelected() & AliVEvent::kMB); // for crosschecks
    break;
  default:
    AliDebug(2, "No trigger selected");
  }

  if (!triggerHM && !triggerMB)
  {
    AliDebug(2, "Miss trigger");
    //     ((TH1*)fHistosQA->At(kEfficiency))->Fill(1);
  }
  if (triggerMB)
    fEvInfo->SetTriggerMB();
  if (triggerHM)
    fEvInfo->SetTriggerHM();

  // vertex selection
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() >= 1)
  {
    fEvInfo->SetVertex();
    fEvInfo->SetVertexGlob();
  }
  else
  { // try SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if (vertex->GetNContributors() >= 1)
      fEvInfo->SetVertex();
    else
    {
      AliDebug(2, "Miss vertex");
      //       ((TH1*)fHistosQA->At(kEfficiency))->Fill(2);
      //       return;
    }
  }
  /*
  if(!AliPPVsMultUtils::HasNoInconsistentSPDandTrackVertices(fESD)){
    ((TH1*)fHistosQA->At(kEfficiency))->Fill(2);
// 	  return;
  }
  if(!AliPPVsMultUtils::IsINELgtZERO(fESD)){
    ((TH1*)fHistosQA->At(kEfficiency))->Fill(2);
// 	  return;
  }
// 	((TH1*)fHistosQA->At(kEfficiency))->Fill(0);
*/

  fEvInfo->SetVertexZ(vertex->GetZ());

  // pile-up selection
  if (fESD->IsPileupFromSPDInMultBins())
    fEvInfo->SetPileUp();

  // multiplicity
  AliESDtrackCuts *tc(NULL);
  AliMultSelection *MultSelection = (AliMultSelection *)fESD->FindListObject("MultSelection");
  if ((tc = dynamic_cast<AliESDtrackCuts *>(fTrackFilter->GetCuts()->At(0) /*FindObject("std10TC")*/)))
  {
    // MakeMultiplicityESD(fESD, "Combined")
    fEvInfo->SetMultiplicity(AliMESeventInfo::kComb, tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8));
    if (tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1) >= 1)
    {
      fEvInfo->SetMultiplicity(AliMESeventInfo::kSPDtrk, tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 0.8));
    }
    // fEvInfo->SetMultiplicity(AliMESeventInfo::kComb, AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD));
    // MakeMultiplicityESD(fESD, "Global")
    fEvInfo->SetMultiplicity(AliMESeventInfo::kGlob08, tc->CountAcceptedTracks(fESD));
    // V0M
    // fEvInfo->SetMultiplicity(AliMESeventInfo::kV0M, fUtils->GetMultiplicityPercentile(fESD, "V0M"));
    fEvInfo->SetMultiplicity(AliMESeventInfo::kV0M, MultSelection->GetMultiplicityPercentile("V0M"));
    // Combined multiplicity for eta (-0.8,-0.4) & (0.4, 0.8)
    fEvInfo->SetMultiplicity(AliMESeventInfo::kComb0408, (tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8) - tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.4)));
  }
  else
  {
    AliWarning("No track cuts defined. No multiplicity computed for REC data.");
    fTrackFilter->GetCuts()->ls();
  }

  Double_t val[7] = {0.};
  THnSparse *H(NULL);
  H = (THnSparse *)fHistosQA->At(kTrkInfo);
  AliMEStrackInfo *tmes(NULL);

  // Int_t counter = 0;
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++)
  {

    AliESDtrack *track = fESD->GetTrack(iTracks);
    if (!track)
    {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    // printf("label[%d] ITS[%d] TPC[%d] TRD[%d]\n", track->GetLabel(),
    //  track->GetITSLabel(), track->GetTPCLabel(), track->GetTRDLabel());

    if (!(fTrackFilter->IsSelected(track)))
    {
      // printf("ESD track reject %d\n", iTracks);
      // track->Print("");
      continue;
    }
    tmes = new AliMEStrackInfo(track, pidResponse, fPIDcomb);

    // TOF matching (old version - Alex)
    // 	if(pidResponse->CheckPIDStatus(AliPIDResponse::kTOF, track) == AliPIDResponse::kDetPidOk) tmes->SetTOFmisProb(fPIDcomb->GetTOFmismatchProb());
    // 	printf("mismatch prob from AliMEStenderV2 = %g\n", fPIDcomb->GetTOFmismatchProb());

    // TOF matching (Cristi)
    Int_t fTOFout;
    if ((!(track->GetStatus() & AliESDtrack::kTOFout)) == 0)
    {
      fTOFout = 1;
    }
    else
    {
      fTOFout = 0;
    }
    Int_t ftime;
    if ((!(track->GetStatus() & AliESDtrack::kTIME)) == 0)
    {
      ftime = 1;
    }
    else
    {
      ftime = 0;
    }
    Double_t flength;
    flength = track->GetIntegratedLength();
    Double_t ftimetof;
    ftimetof = track->GetTOFsignal();
    Double_t inttime[5];
    track->GetIntegratedTimes(inttime); // Returns the array with integrated times for each particle hypothesis
    Double_t fexptimepi;
    Double_t fexptimeka;
    Double_t fexptimepr;
    fexptimepi = inttime[2];
    fexptimeka = inttime[3];
    fexptimepr = inttime[4];

    // Barbara
    if ((fTOFout == 1) && (ftime == 1) && (flength > 350) && (ftimetof > 10000) && (fexptimepi > 10000) && (fexptimeka > 10000) && (fexptimepr > 10000) && (ftimetof < 80000))
    {
      tmes->SetTOFmisProb(1);
    }
    else
    {
      tmes->SetTOFmisProb(0);
    }

    // counter ++;

    // fill tracks QA
    val[0] = tmes->Pt();
    val[1] = tmes->Eta();
    val[2] = tmes->Phi();
    val[3] = tmes->Rv();
    val[4] = tmes->Zv();
    const AliMEStrackInfo::AliMESpid *pid = tmes->GetPID();
    if (pid)
    {
      val[5] = pid->GetRaw(AliMEStrackInfo::kTPC);
      val[6] = pid->GetRaw(AliMEStrackInfo::kTOF);
    }
    if (H)
      H->Fill(val);
    fTracks->Add(tmes);

    // printf("tender: counter = %i eta = %f \t pT = %g\n", counter, val[1], val[0]);
  }

  // printf("tender: combined08 = %g \t old combined08 = %i \t global08 = %g\n", fEvInfo->GetMultiplicity(AliMESeventInfo::kComb), (tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8)), fEvInfo->GetMultiplicity(AliMESeventInfo::kGlob08));
  // printf("fTracks = %i\n\n", fTracks->GetEntries());
  // printf("event index = %i\n", fESD->GetEventNumberInFile());

  // leading particle
  // printf("\n\nNew event!\n");

  // printf("LP search\n");
  fEvInfo->FindLeadingParticle(fTracks);
  // shape
  fEvInfo->MakeDirectivity(fTracks);
  fEvInfo->MakeSphericity(fTracks);
  // fill event QA
  H = (THnSparse *)fHistosQA->At(kEvInfo);
  val[0] = fEvInfo->GetVertexZ();
  val[1] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);
  val[2] = fEvInfo->GetEventShape()->GetDirectivity(kTRUE);
  val[3] = fEvInfo->GetEventShape()->GetDirectivity(kFALSE);
  if (H)
    H->Fill(val);

  Double_t fMultComb08(0.), fMultSPDtrk08(0.), fV0M(0.), fSphericity(0.), fPhiLP(0.), fPtLP(0.), fEtaLP(0.);
  fMultComb08 = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);
  fMultSPDtrk08 = fEvInfo->GetMultiplicity(AliMESeventInfo::kSPDtrk);
  fV0M = fEvInfo->GetMultiplicity(AliMESeventInfo::kV0M);
  fSphericity = fEvInfo->GetEventShape()->GetSphericity();
  Double_t pxLP = fEvInfo->GetEventShape()->GetMomLeading(kTRUE);
  Double_t pyLP = fEvInfo->GetEventShape()->GetMomLeading(kFALSE);
  fPhiLP = TMath::ATan2(pyLP, pxLP);
  fPhiLP = (fPhiLP > 0) ? fPhiLP : (fPhiLP + TMath::TwoPi()); // if negative add 2*pi
  Int_t nTracks = fTracks->GetEntriesFast();
  Int_t run = fESD->GetRunNumber();
  Double_t dca[2] = {0.};
  Double_t fPt[nTracks]{}, fEta[nTracks]{}, fPhi[nTracks]{}, fCharge[nTracks]{}, fDeltaPhi[nTracks]{}, fDeltaEta[nTracks]{}, fDCAxy[nTracks]{}, fDCAz[nTracks]{};
  for (int i(0); i < fTracks->GetEntriesFast(); i++)
  {
    AliMEStrackInfo *t = (AliMEStrackInfo *)(*fTracks)[i];
    fPt[i] = t->Pt();
    fCharge[i] = t->Charge();
    fEta[i] = t->Eta();
    fPhi[i] = t->Phi();
    if (i == 0)
    {
      fPtLP = fPt[i];
      fEtaLP = fEta[i];
    }
    fDeltaPhi[i] = ComputeDeltaPhi(fPhi[i], fPhiLP);
    fDeltaEta[i] = fEtaLP - fEta[i];
    t->GetDCA(dca);
    fDCAxy[i] = dca[0];
    fDCAz[i] = dca[1];
  }
  if (!fTreeSRedirector)
    return;
  (*fTreeSRedirector) << "data/rec"
                      << "run=" << run
                      << "MultComb08=" << fMultComb08
                      << "MultSPDtrk08=" << fMultSPDtrk08
                      << "V0M=" << fV0M
                      << "Sphericity=" << fSphericity
                      << "nTracks=" << nTracks
                      << "Pt[nTracks]=" << fPt[nTracks]
                      << "Charge[nTracks]=" << fCharge[nTracks]
                      << "Eta[nTracks]=" << fEta[nTracks]
                      << "Phi[nTracks]=" << fPhi[nTracks]
                      << "PtLP=" << fPtLP
                      << "EtaLP=" << fEtaLP
                      << "PhiLP=" << fPhiLP
                      << "DeltaPhiLP[nTracks]=" << fDeltaPhi[nTracks]
                      << "DeltaEta[nTracks]=" << fDeltaEta[nTracks]
                      << "fDCAxy[nTracks]=" << fDCAxy[nTracks]
                      << "fDCAz[nTracks]=" << fDCAz[nTracks]
                      << "\n";

  if (!fEventCutsQA.AcceptEvent(ev))
  {
    PostData(kQA, fHistosQA);
    return;
  }

  PostData(kEventInfo + 1, fEvInfo);
  PostData(kTracks + 1, fTracks);
  PostData(kTree + 1, fTree);

  //____ _________________________________
  if (!HasMCdata())
    return;
  //   AliMESeventInfo *fMCevInfo   = dynamic_cast<AliMESeventInfo*>(GetOutputData(AliMESbaseTask::kMCeventInfo+1));
  if (!fMCevInfo)
  {
    AliError("MC event info missing. MC processing skipped");
    return;
  }
  //   TObjArray *fMCtracks         = dynamic_cast<TObjArray*>(GetOutputData(AliMESbaseTask::kMCtracks+1));
  if (!fMCtracks)
  {
    AliError("MC track array missing. MC processing skipped");
    return;
  }
  fMCevInfo->Clear("");
  fMCtracks->Delete();

  AliMCEvent *fMC = dynamic_cast<AliMCEvent *>(MCEvent());
  if (!fMC)
  {
    AliError("MC event not available.");
    return;
  }

  AliStack *fMCStack = fMC->Stack();
  if (!fMCStack)
  {
    AliError("MC stack not available.");
    return;
  }

  // multiplicity
  // multiplicity for eta (-0.8, 0.8)
  fMCevInfo->SetMultiplicity(AliMESeventInfo::kGlob08, MakeMultiplicityMC(fMC));
  // multiplicity for eta (-0.8,-0.4) & (0.4, 0.8)
  fMCevInfo->SetMultiplicity(AliMESeventInfo::kComb0408, MakeMultiplicity0408MC(fMC));
  // multiplicity for eta (-3.7,-1.7) & (2.8, 5.1)  -> V0M
  fMCevInfo->SetMultiplicity(AliMESeventInfo::kV0M, MakeMultiplicityV0MMC(fMC));

  memset(val, 0, 7 * sizeof(Double_t));
  H = (THnSparse *)fHistosQA->At(kMCtrkInfo);
  AliMCParticle *particle(NULL);
  AliMEStrackInfo *tmesRec(NULL);
  for (Int_t ipart = 0; ipart < fMC->GetNumberOfTracks(); ipart++)
  {

    if (!(particle = dynamic_cast<AliMCParticle *>(fMC->GetTrack(ipart))))
    {
      AliWarning("MC particle pointer is null !!!");
      continue;
    }
    if (particle->E() - TMath::Abs(particle->Pz()) < 0.)
    {
      AliWarning("pz > E !!");
      continue;
    }
    // Reject neutral particles (gamma, neutrons, pi0, etc) and quarks/gluons
    if (TMath::Abs(particle->Charge()) < 3)
    {
      // particle->Particle()->Print("");
      continue;
    }
    tmes = new AliMEStrackInfo(particle, fMCStack);
    // if (!(fTrackFilter->IsSelected(tmes)))
    // {
    //   printf("MC track reject %d\n", ipart);
    //   continue;
    // }
    fMCtracks->AddLast(tmes);
    // printf("accept[%d] -> %d\n", ipart, tmes->GetLabel());

    // fill tracks QA
    val[0] = tmes->Pt();
    val[1] = tmes->Eta();
    val[2] = tmes->Phi();
    val[3] = tmes->Rv();
    val[4] = tmes->Zv();
    const AliMEStrackInfo::AliMESpid *pid = tmes->GetPID();
    if (pid)
    {
      const Double_t *prob = pid->GetProb(AliMEStrackInfo::kITS);
      for (Int_t is(0); is < AliPID::kSPECIES; is++)
        if (prob[is] > 0.)
        {
          val[5] = is;
          break;
        }
    }
    val[6] = 0.;
    if (H)
      H->Fill(val);
  }
  //   printf("Found %d / %d ESD tracks MC tracks %d / %d\n",
  //  fTracks->GetEntriesFast(), fESD->GetNumberOfTracks(), fMCtracks->GetEntriesFast(), fMC->GetNumberOfTracks());

  //   Int_t nTracks = fTracks->GetEntriesFast(),
  Int_t nTracksMC = fMCtracks->GetEntriesFast();

  // leading particle
  // printf("generated LP search\n");
  fMCevInfo->FindLeadingParticle(fMCtracks);
  // shape
  fMCevInfo->MakeDirectivity(fMCtracks);
  fMCevInfo->MakeSphericity(fMCtracks);
  // fill event QA
  H = (THnSparse *)fHistosQA->At(kMCevInfo);
  val[0] = 0.;
  val[1] = fMCevInfo->GetMultiplicity(AliMESeventInfo::kGlob08);
  val[2] = fMCevInfo->GetEventShape()->GetDirectivity(kTRUE);
  val[3] = fMCevInfo->GetEventShape()->GetDirectivity(kFALSE);
  if (H)
    H->Fill(val);

  // Re-mapping ESD label to the MC label after both ESD and MC ?! track filtering
  // define matching with ESD track array
  std::vector<Int_t> alreadyMatched;
  for (Int_t ipart = 0; ipart < fMCtracks->GetEntries(); ipart++)
  {
    if (!(tmes = (AliMEStrackInfo *)fMCtracks->At(ipart)))
    {
      AliError("Missing MC track element !");
      continue;
    }
    // printf(" MC[%d] label[%d]\n", ipart, tmes->GetLabel());
    Bool_t found = false;
    for (Int_t iesd(0); iesd < fTracks->GetEntries(); iesd++)
    {
      if (!(tmesRec = (AliMEStrackInfo *)fTracks->At(iesd)))
      {
        AliError("Missing ESD track element !");
        continue;
      }
      if (std::find(alreadyMatched.begin(), alreadyMatched.end(), iesd) != alreadyMatched.end())
        continue;
      if (tmesRec->GetLabel() != tmes->GetLabel())
        continue;
      found = true;
      // printf("ESD[%d] -> label[%d]\n", iesd, tmesRec->GetLabel());
      alreadyMatched.push_back(iesd);
      tmesRec->SetLabel(ipart);
      tmes->SetLabel(iesd);
      break;
    }
    if (!found)
    {
      // printf("   no ESD track \n");
      tmes->SetLabel(-1);
    }
  }

  Double_t fMult08_MC(0.), fV0M_MC(0.), fSphericity_MC(0.), fPhiLP_MC(0.), fPtLP_MC(0.), fEtaLP_MC(0.);
  fMult08_MC = fMCevInfo->GetMultiplicity(AliMESeventInfo::kGlob08);
  fV0M_MC = fMCevInfo->GetMultiplicity(AliMESeventInfo::kV0M);
  fSphericity_MC = fMCevInfo->GetEventShape()->GetSphericity();
  Double_t pxLP_MC = fMCevInfo->GetEventShape()->GetMomLeading(kTRUE);
  Double_t pyLP_MC = fMCevInfo->GetEventShape()->GetMomLeading(kFALSE);
  fPhiLP_MC = TMath::ATan2(pyLP_MC, pxLP_MC);
  fPhiLP_MC = (fPhiLP_MC > 0) ? fPhiLP_MC : (fPhiLP_MC + TMath::TwoPi());
  Int_t nTracks_MC = fMCtracks->GetEntriesFast();
  Double_t fPt_MC[nTracks]{}, fEta_MC[nTracks]{}, fPhi_MC[nTracks]{}, fCharge_MC[nTracks]{}, fDeltaPhi_MC[nTracks]{}, fDeltaEta_MC[nTracks]{}, fPrimary_MC[nTracks]{}, fSecondary_MC[nTracks]{}, fMaterial_MC[nTracks]{};
  AliMEStrackInfo *t(NULL), *tMC(NULL);
  for (int i(0); i < fTracks->GetEntriesFast(); i++)
  {
    if (!(t = (AliMEStrackInfo *)(*fTracks)[i]))
      continue;
    // std::cout << "trk id " << i << " MC label " << t->GetLabel() << std::endl;
    if (t->GetLabel() >= nTracks_MC)
    {
      AliError(Form("MC label %d request outside range %d", t->GetLabel(), nTracks_MC));
      continue;
    }
    if (!(tMC = (AliMEStrackInfo *)fMCtracks->At(t->GetLabel())))
      continue;
    // std::cout << "MC trk id " << t->GetLabel() << " ESD label " << tMC->GetLabel() << std::endl;
    if (tMC->GetLabel() != i)
      AliError(Form("ESD label %d from MC track differ from ESD id %d", tMC->GetLabel(), i));
    fPt_MC[i] = tMC->Pt();
    fCharge_MC[i] = tMC->Charge();
    fEta_MC[i] = tMC->Eta();
    fPhi_MC[i] = tMC->Phi();
    if (i == 0)
    {
      fPtLP_MC = fPt_MC[i];
      fEtaLP_MC = fEta_MC[i];
    }
    fDeltaPhi_MC[i] = ComputeDeltaPhi(fPhi_MC[i], fPhiLP_MC);
    fDeltaEta_MC[i] = fEtaLP_MC - fEta_MC[i];
    if (tMC->HasOrigin(AliMEStrackInfo::kPrimary))
    {
      fPrimary_MC[i] = 1.;
      fSecondary_MC[i] = 0.;
      fMaterial_MC[i] = 0.;
    }
    else if (tMC->HasOrigin(AliMEStrackInfo::kSecondary))
    {
      fPrimary_MC[i] = 0.;
      fSecondary_MC[i] = 1.;
      fMaterial_MC[i] = 0.;
    }

    else if (tMC->HasOrigin(AliMEStrackInfo::kMaterial))
    {
      fPrimary_MC[i] = 0.;
      fSecondary_MC[i] = 0.;
      fMaterial_MC[i] = 1.;
    }
  }

  // if (DebugLevel() >= 1)
  // {
  //   (*DebugStream()) << "MCeventInfo"
  //                                    << "Mult08=" << fMult08_MC
  //                                    << "V0M_MC=" << fV0M_MC
  //                                    << "Sphericity_MC=" << fSphericity_MC
  //                                    << "nTracks=" << nTracks
  //                                    << "Pt_MC[nTracks]=" << fPt_MC[nTracks]
  //                                    << "Charge_MC[nTracks]=" << fCharge_MC[nTracks]
  //                                    << "Eta_MC[nTracks]=" << fEta_MC[nTracks]
  //                                    << "Phi_MC[nTracks]=" << fPhi_MC[nTracks]
  //                                    << "PtLP_MC=" << fPtLP_MC
  //                                    << "EtaLP_MC=" << fEtaLP_MC
  //                                    << "PhiLP_MC=" << fPhiLP_MC
  //                                    << "DeltaPhiLP_MC[nTracks]=" << fDeltaPhi_MC[nTracks]
  //                                    << "DeltaEta_MC[nTracks]=" << fDeltaEta_MC[nTracks]
  //                                    << "fPrimary_MC[nTracks]=" << fPrimary_MC[nTracks]
  //                                    << "fSecondary_MC[nTracks]=" << fSecondary_MC[nTracks]
  //                                    << "fMaterial_MC[nTracks]=" << fMaterial_MC[nTracks]
  //                                    << "\n";
  // }
  // cout << "Debug save " << fTracksIO->GetEntriesFast() << " MC " << fMCtracksIO->GetEntriesFast() << endl;
  // for (int i = 0; i < fTracksIO->GetEntries(); i++)
  // {
  //   std::cout << "REC: index" << i << "constructed at" << fTracksIO->ConstructedAt(i) << std::endl;
  //   // std::cout << "GEN: index" << i << "constructed at" << fMCtracksIO->ConstructedAt(i) << std::endl;
  // }
  // for (int i = 0; i < fMCtracksIO->GetEntries(); i++)
  // {
  //   std::cout << "GEN: index" << i << "constructed at" << fMCtracksIO->ConstructedAt(i) << std::endl;
  // }
  // printf("tracksIn %d tracksOut %d\n", fTracks->GetEntries(), fTracksIO->GetEntries());
  // printf("MCtracksIn %d MCtracksOut %d\n", fMCtracks->GetEntries(), fMCtracksIO->GetEntries());

  Double_t fPt_Gen[nTracks_MC]{}, fEta_Gen[nTracks_MC]{}, fPhi_Gen[nTracks_MC]{}, fCharge_Gen[nTracks_MC]{}, fDeltaPhi_Gen[nTracks_MC]{}, fDeltaEta_Gen[nTracks_MC]{}, fPrimary_Gen[nTracks_MC]{}, fSecondary_Gen[nTracks_MC]{}, fMaterial_Gen[nTracks_MC]{};
  for (int i(0); i < fMCtracks->GetEntriesFast(); i++)
  {
    tMC = (AliMEStrackInfo *)(*fMCtracks)[i];
    // std::cout << i << " ESD label " << tMC->GetLabel() << " reco tracks " << fTracksIO->GetEntries()  << std::endl;
    if (!tMC)
    {
      AliError(Form("Missing MC trk at %d", i));
      continue;
    }
    fPt_Gen[i] = tMC->Pt();
    fCharge_Gen[i] = tMC->Charge();
    fEta_Gen[i] = tMC->Eta();
    fPhi_Gen[i] = tMC->Phi();
    fDeltaPhi_Gen[i] = ComputeDeltaPhi(fPhi_Gen[i], fPhiLP_MC);
    fDeltaEta_Gen[i] = fEtaLP_MC - fEta_Gen[i];
    if (tMC->HasOrigin(AliMEStrackInfo::kPrimary))
    {
      fPrimary_Gen[i] = 1.;
      fSecondary_Gen[i] = 0.;
      fMaterial_Gen[i] = 0.;
    }
    else if (tMC->HasOrigin(AliMEStrackInfo::kSecondary))
    {
      fPrimary_Gen[i] = 0.;
      fSecondary_Gen[i] = 1.;
      fMaterial_Gen[i] = 0.;
    }

    else if (tMC->HasOrigin(AliMEStrackInfo::kMaterial))
    {
      fPrimary_Gen[i] = 0.;
      fSecondary_Gen[i] = 0.;
      fMaterial_Gen[i] = 1.;
    }
    // new ((*fMCGenTracksIO)[j++]) AliMEStrackInfo(*tMC);
  }
  // if (DebugLevel() >= 1)
  // {
  //   (*DebugStream()) << "GenTracksInfo"
  //                                    << "nTracks_MC=" << nTracks_MC
  //                                    << "Pt_Gen[nTracks_MC]=" << fPt_Gen[nTracks_MC]
  //                                    << "Charge_Gen[nTracks_MC]=" << fCharge_Gen[nTracks_MC]
  //                                    << "Eta_Gen[nTracks_MC]=" << fEta_Gen[nTracks_MC]
  //                                    << "Phi_Gen[nTracks_MC]=" << fPhi_Gen[nTracks_MC]
  //                                    << "DeltaPhiLP_Gen[nTracks_MC]=" << fDeltaPhi_Gen[nTracks_MC]
  //                                    << "DeltaEta_Gen[nTracks_MC]=" << fDeltaEta_Gen[nTracks_MC]
  //                                    << "fPrimary_Gen[nTracks_MC]=" << fPrimary_Gen[nTracks_MC]
  //                                    << "fSecondary_Gen[nTracks_MC]=" << fSecondary_Gen[nTracks_MC]
  //                                    << "fMaterial_Gen[nTracks_MC]=" << fMaterial_Gen[nTracks_MC]
  //                                    << "\n";
  // }
  // printf("MCtracksGenIn %d MCtracksGenOut %d\n", fMCtracks->GetEntries(), fMCGenTracksIO->GetEntries());

  Int_t nTracks_Miss = nTracks_MC - nTracks;
  Double_t fPt_Miss[nTracks_Miss]{}, fEta_Miss[nTracks_Miss]{}, fPhi_Miss[nTracks_Miss]{}, fCharge_Miss[nTracks_Miss]{}, fDeltaPhi_Miss[nTracks_Miss]{}, fDeltaEta_Miss[nTracks_Miss]{}, fPrimary_Miss[nTracks_Miss]{}, fSecondary_Miss[nTracks_Miss]{}, fMaterial_Miss[nTracks_Miss]{};
  for (int i(0); i < fMCtracks->GetEntriesFast(); i++)
  {
    tMC = (AliMEStrackInfo *)fMCtracks->At(i);
    if (tMC->GetLabel() < 0)
    {
      fPt_Miss[i] = tMC->Pt();
      fCharge_Miss[i] = tMC->Charge();
      fEta_Miss[i] = tMC->Eta();
      fPhi_Miss[i] = tMC->Phi();
      fDeltaPhi_Miss[i] = ComputeDeltaPhi(fPhi_Miss[i], fPhiLP_MC);
      fDeltaEta_Miss[i] = fEtaLP_MC - fEta_Miss[i];
      if (tMC->HasOrigin(AliMEStrackInfo::kPrimary))
      {
        fPrimary_Miss[i] = 1.;
        fSecondary_Miss[i] = 0.;
        fMaterial_Miss[i] = 0.;
      }
      else if (tMC->HasOrigin(AliMEStrackInfo::kSecondary))
      {
        fPrimary_Miss[i] = 0.;
        fSecondary_Miss[i] = 1.;
        fMaterial_Miss[i] = 0.;
      }

      else if (tMC->HasOrigin(AliMEStrackInfo::kMaterial))
      {
        fPrimary_Miss[i] = 0.;
        fSecondary_Miss[i] = 0.;
        fMaterial_Miss[i] = 1.;
      }
    }
    // new ((*fMCtracksMissIO)[j++]) AliMEStrackInfo(*tMC);
  }
  // if (DebugLevel() >= 1)
  // {
  //   (*DebugStream()) << "MissedTracksInfo"
  //                                    << "nTracks_Miss=" << nTracks_Miss
  //                                    << "Pt_Miss[nTracks_Miss]=" << fPt_Miss[nTracks_Miss]
  //                                    << "Charge_Miss[nTracks_Miss]=" << fCharge_Miss[nTracks_Miss]
  //                                    << "Eta_Miss[nTracks_Miss]=" << fEta_Miss[nTracks_Miss]
  //                                    << "Phi_Miss[nTracks_Miss]=" << fPhi_Miss[nTracks_Miss]
  //                                    << "DeltaPhiLP_Miss[nTracks_Miss]=" << fDeltaPhi_Miss[nTracks_Miss]
  //                                    << "DeltaEta_Miss[nTracks_Miss]=" << fDeltaEta_Miss[nTracks_Miss]
  //                                    << "fPrimary_Miss[nTracks_Miss]=" << fPrimary_Miss[nTracks_Miss]
  //                                    << "fSecondary_Miss[nTracks_Miss]=" << fSecondary_Miss[nTracks_Miss]
  //                                    << "fMaterial_Miss[nTracks_Miss]=" << fMaterial_Miss[nTracks_Miss]
  //                                    << "\n";
  // }
  // printf("MCtracks selected %d\n", fMCtracks->GetEntries());
  // printf("MCtracks matched %d\n", fMCtracksIO->GetEntries());
  // printf("MCtracks missed %d\n", fMCtracksMissIO->GetEntries());
  // printf("Closure %d\n", fMCtracks->GetEntries() - fMCtracksMissIO->GetEntries() - fMCtracksIO->GetEntries());

  AliDebug(2, Form("Tracks REC[%d] MC[%d]", fTracks->GetEntries(), fMCtracks ? fMCtracks->GetEntries() : 0));

  PostData(kQA, fHistosQA);
  PostData(kMCeventInfo + 1, fMCevInfo);
  PostData(kMCtracks + 1, fMCtracks);
}
//_____________________________________________________________________
void AliMEStenderV2::SetPriors()
{

  fPIDcomb = new AliPIDCombined();
  fPIDcomb->SetSelectedSpecies(AliPID::kSPECIES);

  switch (fConfig.fPIDpriors)
  {
  case AliMESconfigTender::kTPC:
    // default aliroot priors
    AliInfo("Setting default priors ...");
    fPIDcomb->SetDefaultTPCPriors();
    AliInfo("Done setting default priors.");
    break;
  case AliMESconfigTender::kIterative:
  { // data priors identified @ 15.04.2015 by Cristi for LHC10d
    AliInfo("Getting iterative data priors from file...");
    TDirectory *cwd(gDirectory);
    TFile *lPriors = TFile::Open("$ALICE_PHYSICS/PWGLF/SPECTRA/MultEvShape/priorsDist_data_LHC10d_newAliroot.root");
    if (lPriors->IsZombie())
    {
      AliError("Could not open the priors file");
      return;
    }
    const Char_t *pname[] = {"e", "e", "pi", "K", "p"};
    AliInfo("Setting iterative data priors ...");
    for (Int_t i = 0; i < 5; i++)
    {
      fPIDcomb->SetPriorDistribution(AliPID::EParticleType(i), (TH1F *)lPriors->Get(Form("priors_%s_final", pname[i])));
      fPIDcomb->GetPriorDistribution(AliPID::EParticleType(i))->SetDirectory(cwd);
    }
    lPriors->Close();
    delete lPriors;
    for (Int_t i = 0; i < 5; i++)
    {
      TH1 *h(fPIDcomb->GetPriorDistribution(AliPID::EParticleType(i)));
      printf("%s -> %s\n", h->GetName(), h->GetDirectory()->GetName());
    }
    fPIDcomb->SetEnablePriors(kTRUE);
    AliInfo("Done setting iterative data priors.");
    break;
  }
  case AliMESconfigTender::kNoPP:
  { // flat priors
    AliInfo("Setting flat priors ...");
    fPIDcomb->SetEnablePriors(kFALSE); // FLAT priors
    AliInfo("Done setting flat priors.");
    break;
  }
  default:
    AliWarning("No PID priors selected");
    break;
  }
}

//________________________________________________________________________
Bool_t AliMEStenderV2::PostProcess()
{
  return kTRUE;
}

//________________________________________________________
Bool_t AliMEStenderV2::BuildQAHistos()
{
  // Make QA sparse histos for
  // - task configuration
  // - efficiency
  // - event info
  // - track info

  if (fHistosQA)
  { // QA histos already created. Check consistency
    if (fHistosQA->GetEntries() == (HasMCdata() ? kNqa : kMCevInfo) && dynamic_cast<AliMEStenderV2::AliMESconfigTender *>(fHistosQA->At(kConfig)))
    {
      AliInfo("QA histos already created");
      return kTRUE;
    }
    else
    {
      AliError("QA histos wrongly defined.");
      fHistosQA->ls();
      return kFALSE;
    }
  }

  // build QA histos
  fHistosQA = new TList();
  fHistosQA->SetOwner(kTRUE);
  fHistosQA->AddAt(&fConfig, kConfig);

  fEventCutsQA.AddQAplotsToList(fHistosQA, kTRUE);

  TH1 *hEff = new TH1I("hEff", "Cuts efficiency;cut;entries", 3, -0.5, 2.5);
  TAxis *ax = hEff->GetXaxis();
  ax->SetBinLabel(1, "OK");
  ax->SetBinLabel(2, "Trigger");
  ax->SetBinLabel(3, "Vertex");
  //   ax->SetBinLabel(4, "PileUp");

  fHistosQA->AddAt(hEff, kEfficiency);

  TString st;
  THnSparseI *H(NULL);
  if (!(H = (THnSparseI *)gROOT->FindObject("EventInfo")))
  {
    const Int_t ndim(4);
    const Char_t *cldTitle[ndim] = {"Z_{vertex}", "Multiplicity Combined", "D_{+}", "D_{-}"};
    const Int_t cldNbins[ndim] = {100, 100, 20, 20};
    const Double_t cldMin[ndim] = {-15., 0., 0., 0.},
                   cldMax[ndim] = {15., 100., 1., 1.};
    st = "Event Info;";
    for (Int_t idim(0); idim < ndim; idim++)
    {
      st += cldTitle[idim];
      st += ";";
    }
    H = new THnSparseI("EventInfo", st.Data(), ndim, cldNbins, cldMin, cldMax);
  }
  else
    H->Reset();
  fHistosQA->AddAt(H, kEvInfo);

  if (!(H = (THnSparseI *)gROOT->FindObject("TrackInfo")))
  {
    const Int_t ndim(7);
    const Char_t *cldTitle[ndim] = {"p_{t}", "eta", "phi", "DCA_{r}", "DCA_{z}", "dEdx", "p_Beta"};
    const Int_t cldNbins[ndim] = {100, 20, 20, 61, 41, 100, 100};
    const Double_t cldMin[ndim] = {0., -1.0, 0., -3., -2.0, 0., 0.},
                   cldMax[ndim] = {10., 1., (2 * TMath::Pi()), 3., 2., 200., 1.1};
    st = "Track Info;";
    for (Int_t idim(0); idim < ndim; idim++)
    {
      st += cldTitle[idim];
      st += ";";
    }
    H = new THnSparseI("TrackInfo", st.Data(), ndim, cldNbins, cldMin, cldMax);
  }
  else
    H->Reset();
  fHistosQA->AddAt(H, kTrkInfo);

  if (HasMCdata())
  {
    if ((H = (THnSparseI *)fHistosQA->At(kEvInfo)))
    {
      H = (THnSparseI *)H->Clone("MCeventInfo");
      H->Reset();
      fHistosQA->AddAt(H, kMCevInfo);
    }
    else
      AliError("Missing EventInfo sparse for MC");
    if ((H = (THnSparseI *)fHistosQA->At(kTrkInfo)))
    {
      H = (THnSparseI *)H->Clone("MCtrackInfo");
      H->Reset();
      fHistosQA->AddAt(H, kMCtrkInfo);
    }
    else
      AliError("Missing TrackInfo sparse for MC");
  }
  AliInfo("Succesfully build QA");
  fHistosQA->ls();
  return kTRUE;
}

//________________________________________________________________________
Int_t AliMEStenderV2::MakeMultiplicityMC(AliMCEvent *const mc)
{
  AliStack *stack(NULL);
  if (!(stack = mc->Stack()))
    return -1;

  //     Int_t nPrim = stack->GetNprimary();
  Int_t charged(0);
  AliMCParticle *particle(NULL);
  for (Int_t ipart = 0; ipart < mc->GetNumberOfTracks(); ipart++)
  {
    if (!(particle = dynamic_cast<AliMCParticle *>(mc->GetTrack(ipart))))
      continue;

    if (particle->E() - TMath::Abs(particle->Pz()) < 0.)
    {
      printf(" - E - AliMEStenderV2::MakeMultiplicityMC : pz > E !!\n");
      continue;
    }

    if (!(stack->IsPhysicalPrimary(particle->GetLabel())))
      continue;

    //  ---------  Charged  ----------
    if (TMath::Abs(particle->Charge()) < 3)
      continue;

    if (TMath::Abs(particle->Eta()) > 0.8)
      continue;

    charged++;

  } // end track loop

  return charged;
}

Int_t AliMEStenderV2::MakeMultiplicity0408MC(AliMCEvent *const mc)
{
  AliStack *stack(NULL);
  if (!(stack = mc->Stack()))
    return -1;

  //     Int_t nPrim = stack->GetNprimary();
  Int_t charged(0);
  AliMCParticle *particle(NULL);
  for (Int_t ipart = 0; ipart < mc->GetNumberOfTracks(); ipart++)
  {
    if (!(particle = dynamic_cast<AliMCParticle *>(mc->GetTrack(ipart))))
      continue;

    if (particle->E() - TMath::Abs(particle->Pz()) < 0.)
    {
      printf(" - E - AliMEStenderV2::MakeMultiplicityMC : pz > E !!\n");
      continue;
    }

    if (!(stack->IsPhysicalPrimary(particle->GetLabel())))
      continue;

    //  ---------  Charged  ----------
    if (TMath::Abs(particle->Charge()) < 3)
      continue;

    if (TMath::Abs(particle->Eta()) > 0.8)
      continue;
    if (TMath::Abs(particle->Eta()) < 0.4)
      continue;

    charged++;

  } // end track loop

  return charged;
}

Int_t AliMEStenderV2::MakeMultiplicityV0MMC(AliMCEvent *const mc)
{
  AliStack *stack(NULL);
  if (!(stack = mc->Stack()))
    return -1;

  //     Int_t nPrim = stack->GetNprimary();
  Int_t charged(0);
  AliMCParticle *particle(NULL);
  for (Int_t ipart = 0; ipart < mc->GetNumberOfTracks(); ipart++)
  {
    if (!(particle = dynamic_cast<AliMCParticle *>(mc->GetTrack(ipart))))
      continue;

    if (particle->E() - TMath::Abs(particle->Pz()) < 0.)
    {
      printf(" - E - AliMEStenderV2::MakeMultiplicityMC : pz > E !!\n");
      continue;
    }

    if (!(stack->IsPhysicalPrimary(particle->GetLabel())))
      continue;

    //  ---------  Charged  ----------
    if (TMath::Abs(particle->Charge()) < 3)
      continue;

    if (particle->Eta() < 5.1 && particle->Eta() > 2.8)
      charged++; // V0A
    if (particle->Eta() < -1.7 && particle->Eta() > -3.7)
      charged++; // V0C

  } // end track loop

  return charged;
}
//_____________________________________________________________________
Double_t AliMEStenderV2::ComputeDeltaPhi(Double_t phi, Double_t phi_LP)
{
  Double_t delta_phi = phi - phi_LP;
  Double_t result = -9999.;

  if (TMath::Abs(delta_phi) > 0.0001)
  { // avoid LP
    if (delta_phi < -TMath::PiOver2())
    {
      result = delta_phi + TMath::TwoPi();
    }
    else if (delta_phi > (3 * TMath::PiOver2()))
    {
      result = delta_phi - TMath::TwoPi();
    }
    else
    {
      result = delta_phi;
    }
  }
  // printf("\n\ndelta_phi = %g\n", result);
  if (result == -9999)
  {
    return 0.5;
  }
  else if (-0.8 <= result && result <= 0.8)
  { // near peak
    // printf("returning 0.5\n");
    return 1.5;
  }
  else if (2.4 <= result && result <= 4.0)
  { // away peak
    // printf("returning 1.5\n");
    return 2.5;
  }
  else
  {
    // printf("returning 2.5\n");
    return 3.5;
  }
}
//_____________________________________________________________________________
void AliMEStenderV2::FinishTaskOutput()
{
  //
  // Called one at the end
  // locally on working node
  //
  if(fTreeSRedirector)
  {
    delete fTreeSRedirector;
    fTreeSRedirector = NULL;
  }

}