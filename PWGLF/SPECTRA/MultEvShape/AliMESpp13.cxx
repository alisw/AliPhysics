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
#include <AliAnalysisUtils.h>
#include <AliPPVsMultUtils.h>
#include <AliPhysicsSelectionTask.h>
#include "AliMultSelection.h"

#include "AliMESbaseTask.h"
#include "AliMEStender.h"
#include "AliMESppColTask.h"
#include "AliMESpidTask.h"
#include "AliMESpp13.h"
#include "AliMESeventInfo.h"
#include "AliMEStrackInfo.h"
#include "TTree.h"

using namespace std;

int eventsPassAllCutsMC(0);
int eventsPassSLCutsMC(0);
ClassImp(AliMESpp13)
    //________________________________________________________________________
    AliMESpp13::AliMESpp13()
    : AliAnalysisTaskSE(), fTrackFilter(NULL), fTracks(NULL), fEvInfo(NULL), fMCtracks(NULL), fMCevInfo(NULL), fTreeSRedirector(NULL), fTracksIO(NULL), fMCtracksIO(NULL), fMCGenTracksIO(NULL), fMCtracksMissIO(NULL), fTree(NULL), fMCGenTree(NULL), fMCMissTree(NULL), fUtils(NULL), fEventCutsQA(kFALSE)
{
  //
  // Constructor
  //
}

//________________________________________________________________________
AliMESpp13::AliMESpp13(const char *name)
    : AliAnalysisTaskSE(name), fTrackFilter(NULL), fTracks(NULL), fEvInfo(NULL), fMCtracks(NULL), fMCevInfo(NULL), fTreeSRedirector(NULL), fTracksIO(NULL), fMCtracksIO(NULL), fMCGenTracksIO(NULL), fMCtracksMissIO(NULL), fTree(NULL), fMCGenTree(NULL), fMCMissTree(NULL), fUtils(NULL), fEventCutsQA(kFALSE)
{ //
  // Constructor
  //
  DefineOutput(kQA, TList::Class());
  DefineOutput(kTree + 1, TTree::Class());
}

//________________________________________________________________________
void AliMESpp13::SetMCdata(Bool_t mc)
{
  // prepare task for MC processing
  SetBit(kMCdata, mc);
  if (mc)
  {
    DefineOutput(kMCGenTree + 1, TTree::Class());
    // DefineOutput(kMCMissTree + 1, TTree::Class());
  }
}
//________________________________________________________________________
AliMESpp13::~AliMESpp13()
{
  //
  // Destructor
  //
  if (fTrackFilter)
    delete fTrackFilter;

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

  if (fTracksIO)
    fTracksIO->Delete();

  if (fMCtracksIO)
    fMCtracksIO->Delete();

  if (fMCGenTracksIO)
    fMCGenTracksIO->Delete();

  if (fMCtracksMissIO)
    fMCtracksMissIO->Delete();

  // if (fTree)
  //   delete fTree;

  // if (fMCGenTree)
  //   delete fMCGenTree;

  // if (fMCMissTree)
  //   delete fMCMissTree;
}

//________________________________________________________________________
void AliMESpp13::UserCreateOutputObjects()
{

  // Build user objects
  BuildQAHistos();
  PostData(kQA, fHistosQA);

  OpenFile(1);
  TDirectory *savedir = gDirectory;
  fTreeSRedirector = new TTreeSRedirector();
  savedir->cd();

  // ------- track cuts
  fTrackFilter = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts *lTrackCuts(NULL);
  lTrackCuts = new AliESDtrackCuts("trkCuts", "Track Cuts");
  lTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE, 1); // kTRUE for primaries
  fTrackFilter->AddCuts(lTrackCuts);

  fTree = ((*fTreeSRedirector) << "ev").GetTree();

  fTracks = new TObjArray(200);
  fTracks->SetOwner(kTRUE);
  fEvInfo = new AliMESeventInfo;
  PostData(kTree + 1, fTree);

  fUtils = new AliPPVsMultUtils();

  if (!HasMCdata())
    return;
  fMCtracks = new TObjArray(200);
  fMCtracks->SetOwner(kTRUE);
  fMCevInfo = new AliMESeventInfo;

  fMCGenTree = ((*fTreeSRedirector) << "genTrk").GetTree();
  // fMCMissTree = ((*fTreeSRedirector) << "missedTrk").GetTree();

  PostData(kMCGenTree + 1, fMCGenTree);
  // PostData(kMCMissTree + 1, fMCMissTree);
}

#include "AliGRPManager.h"
//________________________________________________________________________
void AliMESpp13::UserExec(Option_t * /*opt*/)
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

  if (!fEvInfo)
  {
    AliError("REC event info missing. Processing skipped");
    return;
  }
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

  // init magnetic field
  if (!TGeoGlobalMagField::Instance()->GetField() && !TGeoGlobalMagField::Instance()->IsLocked())
  {
    AliGRPManager grpManager;
    if (!grpManager.ReadGRPEntry())
      AliError("Cannot get GRP entry for magnetic field");
    if (!grpManager.SetMagField())
      AliError("Problem with magnetic field setup");
  }

  AliAnalysisUtils analysisUtils;
  eventsPassSLCutsMC = 0;
  eventsPassAllCutsMC = 0;
  ((TH1 *)fHistosQA->At(kEfficiency))->Fill(0); // all events
  if (inputHandler->IsEventSelected())
  {
    ((TH1 *)fHistosQA->At(kEfficiency))->Fill(1); // physics selection
    {
      if (fEventCutsQA.PassedCut(AliEventCuts::kTrigger))
      {
        ((TH1 *)fHistosQA->At(kEfficiency))->Fill(2); // trigger selection
        if (fEventCutsQA.PassedCut(AliEventCuts::kDAQincomplete))
        {
          ((TH1 *)fHistosQA->At(kEfficiency))->Fill(3); // DAQ incomplete selection
          if (!analysisUtils.IsSPDClusterVsTrackletBG(fESD))
          {
            ((TH1 *)fHistosQA->At(kEfficiency))->Fill(4); // no background selection
            if (AliPPVsMultUtils::IsNotPileupSPDInMultBins(fESD))
            {
              ((TH1 *)fHistosQA->At(kEfficiency))->Fill(5); // reject pile up evts using SPD
              fEvInfo->SetPileUp();
              if (fEventCutsQA.PassedCut(AliEventCuts::kVertex))
              {
                ((TH1 *)fHistosQA->At(kEfficiency))->Fill(6); // vertex
                eventsPassSLCutsMC = 1;
                if (fEventCutsQA.PassedCut(AliEventCuts::kVertexPosition))
                {
                  ((TH1 *)fHistosQA->At(kEfficiency))->Fill(7); // vertex  position
                  eventsPassSLCutsMC = 1;
                  if (fEventCutsQA.PassedCut(AliEventCuts::kINELgt0))
                  {
                    ((TH1 *)fHistosQA->At(kEfficiency))->Fill(8); // INEL > 0
                    eventsPassAllCutsMC = 1;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if (!fEventCutsQA.AcceptEvent(fESD))
  {
    // PostData(kQA, fHistosQA);
    // eventsPassAllCutsMC = 0;
    return;
  }
  else
  {

    ((TH1 *)fHistosQA->At(kEfficiency))->Fill(9); // All cuts
    // eventsPassAllCutsMC = 1;
  }

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
    }
  }
  fEvInfo->SetVertexZ(vertex->GetZ());

  // multiplicity
  AliESDtrackCuts *tc(NULL);
  AliMultSelection *MultSelection = (AliMultSelection *)fESD->FindListObject("MultSelection");
  if ((tc = dynamic_cast<AliESDtrackCuts *>(fTrackFilter->GetCuts()->At(0))))
  {
    fEvInfo->SetMultiplicity(AliMESeventInfo::kComb, tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8));
    fEvInfo->SetMultiplicity(AliMESeventInfo::kSPDtrk, tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 0.8));
    fEvInfo->SetMultiplicity(AliMESeventInfo::kGlob08, tc->CountAcceptedTracks(fESD));
    fEvInfo->SetMultiplicity(AliMESeventInfo::kV0M, MultSelection->GetMultiplicityPercentile("V0M"));
    fEvInfo->SetMultiplicity(AliMESeventInfo::kComb0408, (tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8) - tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.4))); // Combined multiplicity for eta (-0.8,-0.4) & (0.4, 0.8)
  }
  else
  {
    AliWarning("No track cuts defined. No multiplicity computed for REC data.");
    fTrackFilter->GetCuts()->ls();
  }

  Double_t val[5] = {0.};
  THnSparse *H(NULL);
  H = (THnSparse *)fHistosQA->At(kTrkInfo);
  AliMESppColTask sort;
  AliMESpidTask deltaPhi;
  AliMEStrackInfo *tmes(NULL);
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++)
  {
    AliESDtrack *track = fESD->GetTrack(iTracks);
    if (!track)
    {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    if (!(fTrackFilter->IsSelected(track)) && !DebugLevel())
      continue;
    tmes = new AliMEStrackInfo(track);

    // fill tracks QA
    val[0] = tmes->Pt();
    val[1] = tmes->Eta();
    val[2] = tmes->Phi();
    val[3] = tmes->Rv();
    val[4] = tmes->Zv();
    if (H)
      H->Fill(val);
    fTracks->Add(tmes);
    // printf("tender: counter = %i eta = %f \t pT = %g\n", counter, val[1], val[0]);
  }

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

  Double_t fMultComb08(0.), fMultSPDtrk08(0.), fV0M(0.), fSphericity(0.);
  fMultComb08 = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);
  fMultSPDtrk08 = fEvInfo->GetMultiplicity(AliMESeventInfo::kSPDtrk);
  fV0M = fEvInfo->GetMultiplicity(AliMESeventInfo::kV0M);
  fSphericity = fEvInfo->GetEventShape()->GetSphericity();
  Int_t nTracks = fTracks->GetEntriesFast();
  Int_t run = fESD->GetRunNumber();
  Int_t event = fESD->GetEventNumberInFile();
  // if (!fTracksIO)
  //   fTracksIO = new TClonesArray("AliMEStrackInfo");
  Double_t dca[2] = {0.};
  Double_t fPhiLP(0.), fPtLP(0.), fEtaLP(0.), fPxLP(0.), fPyLP(0.);
  Double_t fPt(0.), fEta(0.), fPhi(0.), fCharge(0.), fDeltaPhi(0.), fDeltaEta(0.), fPx(0.), fPy(0.), fDCAxy(0.), fDCAz(0.), fPassDCA(0.);
  sort.QSortTracks(*fTracks, 0, nTracks);
  for (int i(0); i < fTracks->GetEntriesFast(); i++)
  {
    AliMEStrackInfo *t = (AliMEStrackInfo *)(*fTracks)[i];
    if (TMath::Abs(t->Eta()) >= 0.8 || t->Pt() <= 0.15)
      continue;
    // new ((*fTracksIO)[i]) AliMEStrackInfo(*t);
    fCharge = t->Charge();
    fEta = t->Eta();
    fPhi = t->Phi();
    fPt = t->Pt();
    fPx = t->Px();
    fPy = t->Py();
    // cout << "pT_Rec" << fPt << endl;
    if (i == 0)
    {
      fPtLP = fPt;
      fEtaLP = fEta;
      fPxLP = fPx;
      fPyLP = fPy;
      fPhiLP = TMath::ATan2(fPyLP, fPxLP);
      fPhiLP = (fPhiLP > 0) ? fPhiLP : (fPhiLP + TMath::TwoPi()); // if negative add 2*pi
    }
    // printf("pt[%i]=%f; pT(LP)=%f\n", i, fPt, fPtLP);
    fDeltaPhi = deltaPhi.ComputeDeltaPhi(fPhi, fPhiLP);
    if (i == 0)
    {
      fDeltaEta = -9999;
    }
    else
    {
      fDeltaEta = fEtaLP - fEta;
    }
    t->GetDCA(dca);
    fDCAxy = dca[0];
    fDCAz = dca[1];
    if (TMath::Abs(fDCAxy) < (0.0105 + 0.0350 / TMath::Power(fPt, 1.01)))
    {
      fPassDCA = 1.;
    }
    else
    {
      fPassDCA = 0.;
    }
    if (!fTreeSRedirector)
      return;
    if (!HasMCdata())
    {
      (*fTreeSRedirector) << "ev"
                          // << "run=" << run
                          << "MultComb08=" << fMultComb08
                          << "MultSPDtrk08=" << fMultSPDtrk08
                          // << "V0M=" << fV0M
                          << "Sphericity=" << fSphericity
                          << "Pt=" << fPt
                          // << "Eta=" << fEta
                          // << "Phi=" << fPhi.at(i)
                          // << "Charge=" << fCharge.at(i)
                          // << "DeltaPhi=" << fDeltaPhi.at(i)
                          // << "DeltaEta=" << fDeltaEta.at(i)
                          << "DCAxy=" << fDCAxy
                          << "PassDCA=" << fPassDCA
                          << "\n";
    }
  }
  // cout << "entries Rec " << fTracksIO->GetEntries() << endl;
  // for (int i(0); i < fTracksIO->GetEntries(); i++)
  // {
  //   AliMEStrackInfo *t = (AliMEStrackInfo *)(*fTracksIO)[i];
  //   // if (TMath::Abs(t->Eta()) > 0.8 && t->Pt() < 0.15)
  //   //   continue;
  //   fCharge.push_back(t->Charge());
  //   fEta.push_back(t->Eta());
  //   fPhi.push_back(t->Phi());
  //   fPt.push_back(t->Pt());
  //   fPx.push_back(t->Px());
  //   fPy.push_back(t->Py());
  //   // cout << "pT_Rec" << fPt.at(i) << endl;
  //   if (i == 0)
  //   {
  //     fPtLP = fPt.at(i);
  //     fEtaLP = fEta.at(i);
  //     fPxLP = fPx.at(i);
  //     fPyLP = fPy.at(i);
  //     fPhiLP = TMath::ATan2(fPyLP, fPxLP);
  //     fPhiLP = (fPhiLP > 0) ? fPhiLP : (fPhiLP + TMath::TwoPi()); // if negative add 2*pi
  //   }
  //   // printf("pt[%i]=%f; pT(LP)=%f\n", i, fPt, fPtLP);
  //   fDeltaPhi.push_back(deltaPhi.ComputeDeltaPhi(fPhi.at(i), fPhiLP));
  //   if (i == 0)
  //   {
  //     fDeltaEta.push_back(-9999);
  //   }
  //   else
  //   {
  //     fDeltaEta.push_back(fEtaLP - fEta.at(i));
  //   }
  //   t->GetDCA(dca);
  //   fDCAxy.push_back(dca[0]);
  //   fDCAz.push_back(dca[1]);
  //   if (TMath::Abs(fDCAxy.at(i)) < (0.0182 + 0.0350 / TMath::Power(fPt.at(i), 1.01)))
  //   {
  //     fPassDCA.push_back(1.);
  //   }
  //   else
  //   {
  //     fPassDCA.push_back(0.);
  //   }
  //   if (!fTreeSRedirector)
  //     return;
  //   if (!HasMCdata())
  //   {
  //     (*fTreeSRedirector) << "ev"
  //                         << "run=" << run
  //                         << "MultComb08=" << fMultComb08
  //                         << "MultSPDtrk08=" << fMultSPDtrk08
  //                         << "V0M=" << fV0M
  //                         << "Sphericity=" << fSphericity
  //                         << "Pt=" << fPt.at(i)
  //                         << "Eta=" << fEta.at(i)
  //                         // << "Phi=" << fPhi.at(i)
  //                         // << "Charge=" << fCharge.at(i)
  //                         // << "DeltaPhi=" << fDeltaPhi.at(i)
  //                         // << "DeltaEta=" << fDeltaEta.at(i)
  //                         << "DCAxy=" << fDCAxy.at(i)
  //                         << "PassDCA=" << fPassDCA.at(i)
  //                         << "\n";
  //   }
  // }
  // fTracksIO->Clear("C");
  PostData(kTree + 1, fTree);
  //____ _________________________________
  if (!HasMCdata())
    return;
  if (!fMCevInfo)
  {
    AliError("MC event info missing. MC processing skipped");
    return;
  }
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
  fMCevInfo->SetMultiplicity(AliMESeventInfo::kGlob08, AliMEStender::MakeMultiplicityMC(fMC));       // multiplicity for eta (-0.8, 0.8)
  fMCevInfo->SetMultiplicity(AliMESeventInfo::kComb0408, AliMEStender::MakeMultiplicity0408MC(fMC)); // multiplicity for eta (-0.8,-0.4) & (0.4, 0.8)
  fMCevInfo->SetMultiplicity(AliMESeventInfo::kV0M, AliMEStender::MakeMultiplicityV0MMC(fMC));       // multiplicity for eta (-3.7,-1.7) & (2.8, 5.1)  -> V0M

  memset(val, 0, 5 * sizeof(Double_t));
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
    fMCtracks->AddLast(tmes);
    // printf("accept[%d] -> %d\n", ipart, tmes->GetLabel());

    // fill tracks QA
    val[0] = tmes->Pt();
    val[1] = tmes->Eta();
    val[2] = tmes->Phi();
    val[3] = tmes->Rv();
    val[4] = tmes->Zv();
    if (H)
      H->Fill(val);
  }
  //   printf("Found %d / %d ESD tracks MC tracks %d / %d\n",
  //  fTracks->GetEntriesFast(), fESD->GetNumberOfTracks(), fMCtracks->GetEntriesFast(), fMC->GetNumberOfTracks());

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

  Double_t fMult08_MC(0.), fV0M_MC(0.), fSphericity_MC(0.);
  fMult08_MC = fMCevInfo->GetMultiplicity(AliMESeventInfo::kGlob08);
  fV0M_MC = fMCevInfo->GetMultiplicity(AliMESeventInfo::kV0M);
  fSphericity_MC = fMCevInfo->GetEventShape()->GetSphericity();
  Int_t nTracks_MC = fMCtracks->GetEntriesFast();
  Int_t nTracksMissed = nTracks_MC - nTracks;
  AliMEStrackInfo *t(NULL), *tMC(NULL);
  // if (!fTracksIO)
  //   fTracksIO = new TClonesArray("AliMEStrackInfo");
  // if (!fMCtracksIO)
  //   fMCtracksIO = new TClonesArray("AliMEStrackInfo");
  Double_t fPhiLP_MC(0.), fPtLP_MC(0.), fEtaLP_MC(0.), fPxLP_MC(0.), fPyLP_MC(0.);
  Double_t fPt_MC(0.), fEta_MC(0.), fPhi_MC(0.), fCharge_MC(0.), fDeltaPhi_MC(0.), fDeltaEta_MC(0.), fPrimary_MC(0.), fSecondary_MC(0.), fMaterial_MC(0.), fPx_MC(0.), fPy_MC(0.);
  for (int i(0) /*, j(0)*/; i < fTracks->GetEntriesFast(); i++)
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
    {
      AliError(Form("ESD label %d from MC track differ from ESD id %d", tMC->GetLabel(), i));
      continue;
    }
    fPt_MC = tMC->Pt();
    // cout << "pT_Gen matched" << fPt_MC.at(i) << endl;
    fCharge_MC = tMC->Charge();
    fEta_MC = tMC->Eta();
    fPhi_MC = tMC->Phi();
    fPx_MC = tMC->Px();
    fPy_MC = tMC->Py();
    if (i == 0)
    {
      fPtLP_MC = fPt_MC;
      fEtaLP_MC = fEta_MC;
      fPxLP_MC = fPx_MC;
      fPyLP_MC = fPy_MC;
      fPhiLP_MC = TMath::ATan2(fPyLP_MC, fPxLP_MC);
      fPhiLP_MC = (fPhiLP_MC > 0) ? fPhiLP_MC : (fPhiLP_MC + TMath::TwoPi());
    }
    // printf("pT(%i) = %f\n", i, fPt_MC.at(i));
    fDeltaPhi_MC = deltaPhi.ComputeDeltaPhi(fPhi_MC, fPhiLP_MC);
    if (i == 0)
    {
      fDeltaEta_MC = -9999;
    }
    else
    {
      fDeltaEta_MC = fEtaLP_MC - fEta_MC;
    }
    if (tMC->HasOrigin(AliMEStrackInfo::kPrimary))
    {
      fPrimary_MC = 1.;
      fSecondary_MC = 0.;
      fMaterial_MC = 0.;
    }
    else if (tMC->HasOrigin(AliMEStrackInfo::kSecondary))
    {
      fPrimary_MC = 0.;
      fSecondary_MC = 1.;
      fMaterial_MC = 0.;
    }
    else if (tMC->HasOrigin(AliMEStrackInfo::kMaterial))
    {
      fPrimary_MC = 0.;
      fSecondary_MC = 0.;
      fMaterial_MC = 1.;
    }
    if (!fTreeSRedirector)
      return;
    (*fTreeSRedirector) << "ev"
                        // << "run=" << run
                        << "MultComb08=" << fMultComb08
                        << "MultSPDtrk08=" << fMultSPDtrk08
                        // << "V0M=" << fV0M
                        << "Sphericity=" << fSphericity
                        << "Pt=" << fPt
                        // << "Eta=" << fEta
                        // << "Phi=" << fPhi.at(i)
                        // << "DeltaPhi=" << fDeltaPhi.at(i)
                        // << "DeltaEta=" << fDeltaEta.at(i)
                        // << "Charge=" << fCharge.at(i)
                        << "DCAxy=" << fDCAxy
                        << "PassDCA=" << fPassDCA
                        << "Mult08=" << fMult08_MC
                        // << "V0M_MC=" << fV0M_MC
                        << "Sphericity_MC=" << fSphericity_MC
                        // << "EventsPassSLCuts_MC=" << eventsPassSLCutsMC
                        // << "EventsPassAllCuts_MC=" << eventsPassAllCutsMC
                        << "Pt_MC=" << fPt_MC
                        // << "Phi_MC=" << fPhi_MC.at(i)
                        // << "Eta_MC=" << fEta_MC
                        // << "DeltaPhi_MC=" << fDeltaPhi_MC.at(i)
                        // << "DeltaEta_MC=" << fDeltaEta_MC.at(i)
                        // << "Charge_MC=" << fCharge_MC.at(i)
                        << "Primary_MC=" << fPrimary_MC
                        << "Secondary_MC=" << fSecondary_MC
                        << "Material_MC=" << fMaterial_MC
                        << "\n";

    // new ((*fTracksIO)[j]) AliMEStrackInfo(*t);
    // new ((*fMCtracksIO)[j]) AliMEStrackInfo(*tMC);
    // j++;
  }
  // cout << "!!!!!! fTracksIO entries  = " << fTracksIO->GetEntries() << endl;
  // cout << "!!!!!! fMCtracksIO entries  = " << fMCtracksIO->GetEntries() << endl;
  //

  // for (int i(0); i < fMCtracksIO->GetEntries(); i++)
  // {
  //   AliMEStrackInfo *tMC = (AliMEStrackInfo *)(*fMCtracksIO)[i];
  //   fPt_MC.push_back(tMC->Pt());
  //   // cout << "pT_Gen matched" << fPt_MC.at(i) << endl;
  //   fCharge_MC.push_back(tMC->Charge());
  //   fEta_MC.push_back(tMC->Eta());
  //   fPhi_MC.push_back(tMC->Phi());
  //   fPx_MC.push_back(tMC->Px());
  //   fPy_MC.push_back(tMC->Py());
  //   if (i == 0)
  //   {
  //     fPtLP_MC = fPt_MC.at(i);
  //     fEtaLP_MC = fEta_MC.at(i);
  //     fPxLP_MC = fPx_MC.at(i);
  //     fPyLP_MC = fPy_MC.at(i);
  //     fPhiLP_MC = TMath::ATan2(fPyLP_MC, fPxLP_MC);
  //     fPhiLP_MC = (fPhiLP_MC > 0) ? fPhiLP_MC : (fPhiLP_MC + TMath::TwoPi());
  //   }
  //   // printf("pT(%i) = %f\n", i, fPt_MC.at(i));
  //   fDeltaPhi_MC.push_back(deltaPhi.ComputeDeltaPhi(fPhi_MC.at(i), fPhiLP_MC));
  //   if (i == 0)
  //   {
  //     fDeltaEta_MC.push_back(-9999);
  //   }
  //   else
  //   {
  //     fDeltaEta_MC.push_back(fEtaLP_MC - fEta_MC.at(i));
  //   }
  //   if (tMC->HasOrigin(AliMEStrackInfo::kPrimary))
  //   {
  //     fPrimary_MC.push_back(1.);
  //     fSecondary_MC.push_back(0.);
  //     fMaterial_MC.push_back(0.);
  //   }
  //   else if (tMC->HasOrigin(AliMEStrackInfo::kSecondary))
  //   {
  //     fPrimary_MC.push_back(0.);
  //     fSecondary_MC.push_back(1.);
  //     fMaterial_MC.push_back(0.);
  //   }
  //   else if (tMC->HasOrigin(AliMEStrackInfo::kMaterial))
  //   {
  //     fPrimary_MC.push_back(0.);
  //     fSecondary_MC.push_back(0.);
  //     fMaterial_MC.push_back(1.);
  //   }
  //   if (fPrimary_MC.size() != fPt_MC.size() || fSecondary_MC.size() != fPt_MC.size() || fMaterial_MC.size() != fPt_MC.size())
  //     continue;
  //   if (!fTreeSRedirector)
  //     return;
  //   (*fTreeSRedirector) << "ev"
  //                       << "run=" << run
  //                       << "MultComb08=" << fMultComb08
  //                       << "MultSPDtrk08=" << fMultSPDtrk08
  //                       << "V0M=" << fV0M
  //                       << "Sphericity=" << fSphericity
  //                       << "Pt=" << fPt.at(i)
  //                       << "Eta=" << fEta.at(i)
  //                       // << "Phi=" << fPhi.at(i)
  //                       // << "DeltaPhi=" << fDeltaPhi.at(i)
  //                       // << "DeltaEta=" << fDeltaEta.at(i)
  //                       // << "Charge=" << fCharge.at(i)
  //                       // << "DCAxy=" << fDCAxy.at(i)
  //                       // << "PassDCA=" << fPassDCA.at(i)
  //                       << "Mult08=" << fMult08_MC
  //                       << "V0M_MC=" << fV0M_MC
  //                       << "Sphericity_MC=" << fSphericity_MC
  //                       << "EventsPassSLCuts_MC=" << eventsPassSLCutsMC
  //                       << "EventsPassAllCuts_MC=" << eventsPassAllCutsMC
  //                       << "Pt_MC=" << fPt_MC.at(i)
  //                       // << "Phi_MC=" << fPhi_MC.at(i)
  //                       << "Eta_MC=" << fEta_MC.at(i)
  //                       // << "DeltaPhi_MC=" << fDeltaPhi_MC.at(i)
  //                       // << "DeltaEta_MC=" << fDeltaEta_MC.at(i)
  //                       // << "Charge_MC=" << fCharge_MC.at(i)
  //                       << "Primary_MC=" << fPrimary_MC.at(i)
  //                       << "Secondary_MC=" << fSecondary_MC.at(i)
  //                       << "Material_MC=" << fMaterial_MC.at(i)
  //                       << "\n";
  // }

  // // cout << "Debug save " << fTracksIO->GetEntriesFast() << " MC " << fMCtracksIO->GetEntriesFast() << endl;
  // for (int i = 0; i < fTracksIO->GetEntries(); i++)
  // {
  //   // std::cout << "REC: index" << i << "constructed at" << fTracksIO->ConstructedAt(i) << std::endl;
  //   // std::cout << "GEN: index" << i << "constructed at" << fMCtracksIO->ConstructedAt(i) << std::endl;
  // }
  // for (int i = 0; i < fMCtracksIO->GetEntries(); i++)
  // {
  //   // std::cout << "GEN: index" << i << "constructed at" << fMCtracksIO->ConstructedAt(i) << std::endl;
  // }
  // printf("tracksIn %d tracksOut %d\n", fTracks->GetEntries(), fTracksIO->GetEntries());
  // printf("MCtracksIn %d MCtracksOut %d\n", fMCtracks->GetEntries(), fMCtracksIO->GetEntries());

  // if (!fMCGenTracksIO)
  //   fMCGenTracksIO = new TClonesArray("AliMEStrackInfo");
  sort.QSortTracks(*fMCtracks, 0, nTracks_MC);
  Double_t fPt_Gen(0.), fEta_Gen(0.), fPhi_Gen(0.), fEtaLP_Gen(0.), fPtLP_Gen(0.), fCharge_Gen(0.), fDeltaPhi_Gen(0.), fDeltaEta_Gen(0.);
  for (int i(0) /*, j(0)*/; i < nTracks_MC; i++)
  {
    tMC = (AliMEStrackInfo *)(*fMCtracks)[i];
    // std::cout << i << " ESD label " << tMC->GetLabel() << " reco tracks " << fTracksIO->GetEntries()  << std::endl;
    if (!tMC)
    {
      AliError(Form("Missing MC trk at %d", i));
      continue;
    }
    if (TMath::Abs(tMC->Eta()) > 0.8 || tMC->Pt() <= 0.15)
      continue;
    fPt_Gen = tMC->Pt();
    fCharge_Gen = tMC->Charge();
    fEta_Gen = tMC->Eta();
    fPhi_Gen = tMC->Phi();
    // fDeltaPhi_Gen = deltaPhi.ComputeDeltaPhi(fPhi_Gen, fPhiLP_MC);
    if (i == 0)
    {
      fDeltaEta_Gen = -9999;
      fEtaLP_Gen = fEta_Gen;
      fPtLP_Gen = fPt_Gen;
    }
    else
    {
      fDeltaEta_Gen = fEtaLP_Gen - fEta_Gen;
    }
    // if (!fTreeSRedirector)
    //   return;
    // (*fTreeSRedirector) << "genTrk"
    //                     // << "Mult08=" << fMult08_MC
    //                     // << "Sphericity_MC=" << fSphericity_MC
    //                     // << "Pt_Gen=" << fPt_Gen
    //                     // //                     << "Charge_Gen=" << fCharge_Gen
    //                     // << "Eta_Gen=" << fEta_Gen
    //                     // << "Phi_Gen=" << fPhi_Gen
    //                     // << "DeltaPhi_Gen=" << fDeltaPhi_Gen
    //                     // << "DeltaEta_Gen=" << fDeltaEta_Gen
    //                     << "\n";
    // new ((*fMCGenTracksIO)[j++]) AliMEStrackInfo(*tMC);
  }

  // for (int i(0); i < nTracks_MC; i++)
  // {
  //   tMC = (AliMEStrackInfo *)(*fMCGenTracksIO)[i];
  //   if (TMath::Abs(tMC->Eta()) > 0.8 && tMC->Pt() < 0.15)
  //     continue;
  //   fPt_Gen = tMC->Pt();
  //   fCharge_Gen = tMC->Charge();
  //   fEta_Gen = tMC->Eta();
  //   fPhi_Gen = tMC->Phi();
  //   // fDeltaPhi_Gen = deltaPhi.ComputeDeltaPhi(fPhi_Gen, fPhiLP_MC);
  //   if (i == 0)
  //   {
  //     fDeltaEta_Gen = -9999;
  //     fEtaLP_Gen = fEta_Gen;
  //     fPtLP_Gen = fPt_Gen;
  //   }
  //   else
  //   {
  //     fDeltaEta_Gen = fEtaLP_Gen - fEta_Gen;
  //   }
  //   if (!fTreeSRedirector)
  //     return;
  //   (*fTreeSRedirector) << "genTrk"
  //                       << "run=" << run
  //                       << "Mult08=" << fMult08_MC
  //                       << "Sphericity_MC=" << fSphericity_MC
  //                       << "Pt_Gen=" << fPt_Gen
  //                       //                     << "Charge_Gen=" << fCharge_Gen
  //                       << "Eta_Gen=" << fEta_Gen
  //                       << "Primary_MC=" << fPrimary_MC.at(i)
  //                       << "Secondary_MC=" << fSecondary_MC.at(i)
  //                       << "Material_MC=" << fMaterial_MC.at(i)
  //                       // << "Phi_Gen=" << fPhi_Gen
  //                       // << "DeltaPhi_Gen=" << fDeltaPhi_Gen
  //                       // << "DeltaEta_Gen=" << fDeltaEta_Gen
  //                       << "\n";
  // }
  // if (!fTreeSRedirector)
  //   return;
  // (*fTreeSRedirector) << "ev"
  // << "run=" << run
  //                     // << "nTracks_Gen=" << nTracks_MC;
  //                     // << "PtLP_Gen=" << fPtLP_Gen
  //                     // << "EtaLP_Gen=" << fEtaLP_Gen
  //                     // << "PhiLP_Gen=" << fPhiLP_Gen
  //                     << "\n";
  // printf("MCtracksGenIn %d MCtracksGenOut %d\n", fMCtracks->GetEntries(), fMCGenTracksIO->GetEntries());

  // if (!fMCtracksMissIO)
  //   fMCtracksMissIO = new TClonesArray("AliMEStrackInfo");
  // for (int i(0), j(0); i < fMCtracks->GetEntriesFast(); i++)
  // {
  //   tMC = (AliMEStrackInfo *)fMCtracks->At(i);
  //   if (tMC->GetLabel() < 0)
  //     new ((*fMCtracksMissIO)[j++]) AliMEStrackInfo(*tMC);
  // }
  // sort.QSortTracks(*fMCtracksMissIO, 0, nTracksMissed);
  // Double_t fPt_Miss(0.), fEta_Miss(0.), fPtLP_Miss(0.), fEtaLP_Miss(0.), fPhi_Miss(0.), fCharge_Miss(0.), fDeltaPhi_Miss(0.), fDeltaEta_Miss(0.);
  // for (int i(0); i < nTracksMissed; i++)
  // {
  //   tMC = (AliMEStrackInfo *)(*fMCtracksMissIO)[i];
  //   if (TMath::Abs(tMC->Eta()) > 0.8 && tMC->Pt() < 0.15)
  //     continue;
  //   fPt_Miss = tMC->Pt();
  //   fCharge_Miss = tMC->Charge();
  //   fEta_Miss = tMC->Eta();
  //   fPhi_Miss = tMC->Phi();
  //   fDeltaPhi_Miss = deltaPhi.ComputeDeltaPhi(fPhi_Miss, fPhiLP_MC);
  //   if (i == 0)
  //   {
  //     fDeltaEta_Miss = -9999;
  //     fPtLP_Miss = fPt_Miss;
  //     fEtaLP_Miss = fEta_Miss;
  //   }
  //   else
  //   {
  //     fDeltaEta_Miss = fEtaLP_Miss - fEta_Miss;
  //   }
  //   // if (!fTreeSRedirector)
  //   //   return;
  //   // (*fTreeSRedirector) << "missedTrk"
  //   // << "run=" << run
  //   //                     // << "Pt_Miss=" << fPt_Miss
  //   //                     // << "Charge_Miss=" << fCharge_Miss
  //   //                     // << "Eta_Miss=" << fEta_Miss
  //   //                     // << "Phi_Miss=" << fPhi_Miss
  //   //                     // << "DeltaPhi_Miss=" << fDeltaPhi_Miss
  //   //                     // << "DeltaEta_Miss=" << fDeltaEta_Miss
  //   //                     << "\n";
  // }
  // if (!fTreeSRedirector)
  //   return;
  // (*fTreeSRedirector) << "ev"
  // << "run=" << run
  //                     // << "nTracks_Miss=" << nTracksMissed
  //                     // << "PtLP_Miss=" << fPtLP_Miss
  //                     // << "EtaLP_Miss=" << fEtaLP_Miss
  //                     // << "PhiLP_Miss=" << fPhiLP_Miss
  //                     << "\n";

  // printf("MCtracks selected %d\n", fMCtracks->GetEntries());
  // printf("MCtracks matched %d\n", fMCtracksIO->GetEntries());
  // printf("MCtracks missed %d\n", fMCtracksMissIO->GetEntries());
  // printf("Closure %d\n", fMCtracks->GetEntries() - fMCtracksMissIO->GetEntries() - fMCtracksIO->GetEntries());

  // if ((fMCtracks->GetEntries() - fMCtracksMissIO->GetEntries() - fMCtracksIO->GetEntries()) != 0)
  // {
  //   AliError("the closure test for MC is not passed!");
  // }

  AliDebug(2, Form("Tracks REC[%d] MC[%d]", fTracks->GetEntries(), fMCtracks ? fMCtracks->GetEntries() : 0));

  // fTracksIO->Clear("C");
  // fMCtracksIO->Clear("C");
  // fMCGenTracksIO->Clear("C");
  // fMCtracksMissIO->Clear("C");

  PostData(kQA, fHistosQA);
  PostData(kTree + 1, fTree);
  PostData(kMCGenTree + 1, fMCGenTree);
  // PostData(kMCMissTree + 1, fMCMissTree);
}
//________________________________________________________
Bool_t AliMESpp13::BuildQAHistos()
{
  // Make QA sparse histos for
  // - efficiency
  // - event info
  // - track info

  // build QA histos
  fHistosQA = new TList();
  fHistosQA->SetOwner(kTRUE);

  fEventCutsQA.AddQAplotsToList(fHistosQA, kTRUE);

  TH1 *hEff = new TH1I("hEff", "Event selection;;Events", 10, -0.5, 9.5);
  TAxis *ax = hEff->GetXaxis();
  ax->SetBinLabel(1, "All events");
  ax->SetBinLabel(2, "Physics selection");
  ax->SetBinLabel(3, "Trigger (kINT7)");
  ax->SetBinLabel(4, "Complete DAQ");
  ax->SetBinLabel(5, "No background");
  ax->SetBinLabel(6, "No pile up");
  ax->SetBinLabel(7, "Good vertex");
  ax->SetBinLabel(8, "|Zvtx|<10 cm");
  ax->SetBinLabel(9, "INEL > 0");
  ax->SetBinLabel(10, "All cuts");
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
    const Int_t ndim(5);
    const Char_t *cldTitle[ndim] = {"p_{t}", "eta", "phi", "DCA_{r}", "DCA_{z}"};
    const Int_t cldNbins[ndim] = {100, 20, 20, 61, 41};
    const Double_t cldMin[ndim] = {0., -1.0, 0., -3., -2.0},
                   cldMax[ndim] = {10., 1., (2 * TMath::Pi()), 3., 2.};
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
//_____________________________________________________________________________
void AliMESpp13::FinishTaskOutput()
{
  //
  // Called one at the end
  // locally on working node
  //
  if (fTreeSRedirector)
  {
    delete fTreeSRedirector;
    fTreeSRedirector = NULL;
  }
}
