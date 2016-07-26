// ROOT
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TFormula.h"
#include "TRandom.h"
#include "TSpline.h"

// analysis framework
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliVMultiplicity.h"

// MC stuff
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliStack.h"

// ESD stuff
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"

// AOD stuff
#include "AliAODEvent.h"
#include "AliAODJet.h"
#include "AliAODTrack.h"

// EMCAL framework and jet tasks
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliEmcalJet.h"

#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisTaskJetsEvshape.h"

#include <iostream>
#include <cmath>

AliAnalysisTaskJetsEvshape::AliAnalysisTaskJetsEvshape(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fMCEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fRunNumber(-1),
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0),
  fThrust(0.),
  fSphericityT(0.),
  fUseMC(kTRUE),
  fEtaMaxForEvshape(.9),
  fOutputList(0x0),
  fShortTaskId("jets_evshape")
{
  // default ctor

  SetMakeGeneralHistograms(kTRUE);

  DefineOutput(kOutputEmcal, TList::Class());
  AliInfo(Form("creating output slot #%i", kOutputTask));
  DefineOutput(kOutputTask, TList::Class());
}

AliAnalysisTaskJetsEvshape::~AliAnalysisTaskJetsEvshape()
{
  // dtor
}

void AliAnalysisTaskJetsEvshape::UserCreateOutputObjects()
{
  // create user output objects
  AliInfo("creating output objects");

  // common EMCAL framework
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  PostData(kOutputEmcal, fOutput);
  fJetsCont           = GetJetContainer(0);
  printf("setup jet container: %p, from available:\n", fJetsCont);
  fJetCollArray.Print();
  printf("end\n");
  if(fJetsCont) { //get particles and clusters connected to jets
    fJetsCont->PrintCuts();
    fTracksCont       = fJetsCont->GetParticleContainer();
    fCaloClustersCont = fJetsCont->GetClusterContainer();
  } else {        //no jets, just analysis tracks and clusters
    fTracksCont       = GetParticleContainer(0);
    fCaloClustersCont = GetClusterContainer(0);
  }
  if(fTracksCont)
    fTracksCont->SetClassName("AliVParticle");
  if(fCaloClustersCont)
    fCaloClustersCont->SetClassName("AliVCluster");

  // setup list
  OpenFile(kOutputTask);
  fOutputList = new TList();
  fOutputList->SetOwner();

  // setup histograms
  TH1 *hist;
  TH1 *histStat = AddHistogram(ID(kHistStat), "event statistics;;counts",
                               kStatLast-1, .5, kStatLast-.5);
  histStat->GetXaxis()->SetBinLabel(ID(kStatSeen));
  histStat->GetXaxis()->SetBinLabel(ID(kStatTrg));
  histStat->GetXaxis()->SetBinLabel(ID(kStatUsed));
  histStat->GetXaxis()->SetBinLabel(ID(kStatEvCuts));

  AddHistogram(ID(kHistJetPt), "jet spectrum;p_{T}^{jet,ch} (GeV/#it{c});counts",
               40, 0., 40.);
  AddHistogram(ID(kHistMult), "multiplicity;multiplicity;counts",
               100, 0., 4000.);
  AddHistogram(ID(kHistJetPtVsMult), "jet p_{T} vs mult;multiplicity;p_{T}^{jet,ch} (GeV/#it{c})",
               100, 0., 4000.,
               40, 0., 40.);
  AddHistogram(ID(kHistJetPtMultEvshape), "jet spectrum;p_{T}^{jet,ch} (GeV/#it{c});multiplicity;S_{T}",
               40, 0., 40., 100, 0., 500., 50, 0., 1.);
  AddHistogram(ID(kHistThrust), "thrust;T;counts",
               50, 0., 1.);
  AddHistogram(ID(kHistSphericityT), "transverse sphericity;S_{T};counts",
               50, 0., 1.);
  AddHistogram(ID(kHistSphericityTVsMult), "transverse sphericity vs mult;multiplicity;S_{T};counts",
               100, 0., 4000.,
               50, 0., 1.);

  PostData(kOutputTask, fOutputList);
}

Bool_t AliAnalysisTaskJetsEvshape::Notify()
{
  // actions to be taken upon notification about input file change

  return AliAnalysisTaskEmcalJet::Notify();
}

void AliAnalysisTaskJetsEvshape::Terminate(const Option_t *option)
{
  // actions at task termination

  AliAnalysisTaskEmcalJet::Terminate(option);
}

void AliAnalysisTaskJetsEvshape::PrintTask(Option_t *option, Int_t indent) const
{
  AliAnalysisTaskEmcalJet::PrintTask(option, indent);

  std::cout << std::setw(indent) << " " << "nothing to say: " << std::endl;
}

Bool_t AliAnalysisTaskJetsEvshape::PrepareEvent()
{
  Bool_t eventGood = kTRUE;

  // setup pointers to input data (null if unavailable)
  // mcEvent:  MC input
  // esdEvent: ESD input
  // outEvent: AOD output
  // aodEvent: AOD input if available, otherwise AOD output

  fMCEvent   = this->MCEvent();
  fESDEvent  = dynamic_cast<AliESDEvent*>(this->InputEvent()); // could also be AOD input
  AliAODEvent* outEvent  = this->AODEvent();
  fAODEvent  = dynamic_cast<AliAODEvent*> (this->InputEvent());
  if (!fAODEvent)
    fAODEvent = outEvent;

  if ((fDebug > 0) && fESDEvent)
    printf("event: %s-%06i\n", CurrentFileName(), fESDEvent->GetEventNumberInFile());

  // check for run change
  if (fRunNumber != InputEvent()->GetRunNumber()) {
    fRunNumber = InputEvent()->GetRunNumber();
  }

  return eventGood;
}

Bool_t AliAnalysisTaskJetsEvshape::CleanUpEvent()
{
  return kTRUE;
}

void AliAnalysisTaskJetsEvshape::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
}

Bool_t AliAnalysisTaskJetsEvshape::Run()
{
  // record number of sampled events and detect trigger contributions
  FillH1(kHistStat, kStatSeen);

  // so far, no trigger selection, we accept all
  FillH1(kHistStat, kStatTrg);

  // prepare the event
  if (PrepareEvent()) {
    // InputEvent()->GetList()->ls();

    // multiplicity selection

    // event shape selection
    if (!CalculateThrust())
      return kFALSE;

    if (!CalculateSphericityT())
      return kFALSE;

    // here we have passed the event cuts
    FillH1(kHistStat, kStatEvCuts);

    return kTRUE;
  } else {
    return kFALSE;
  }
}

Bool_t AliAnalysisTaskJetsEvshape::FillHistograms()
{
  FillH1(kHistStat, kStatUsed);

  const AliStack *stack = fMCEvent ? fMCEvent->Stack() : 0x0;
  const AliVMultiplicity *mult = InputEvent()->GetMultiplicity();
  const Int_t nTracklets =
    stack ? stack->GetNprimary() :
    mult ? mult->GetNumberOfTracklets() :
    -1;
  FillH1(kHistMult, nTracklets);
  FillH1(kHistThrust, fThrust);
  FillH1(kHistSphericityT, fSphericityT);
  FillH2(kHistSphericityTVsMult, nTracklets, fSphericityT);

  if (fJetsCont) {
    // const Int_t nJets         = fJetsCont->GetNJets();
    // const Int_t nAcceptedJets = fJetsCont->GetNAcceptedJets();

    fJetsCont->ResetCurrentID();
    while (AliEmcalJet *jet = fJetsCont->GetNextAcceptJet()) {
      FillH1(kHistJetPt, jet->Pt());
      FillH2(kHistJetPtVsMult, nTracklets, jet->Pt());
      FillH3(kHistJetPtMultEvshape, jet->Pt(), nTracklets, fSphericityT);
    }
  }

  return kTRUE;
}

void AliAnalysisTaskJetsEvshape::UserExec(Option_t *option)
{
  // actual work

  AliAnalysisTaskEmcalJet::UserExec(option);

  CleanUpEvent();

  PostData(kOutputEmcal, fOutput);
  PostData(kOutputTask, fOutputList);
}

Bool_t AliAnalysisTaskJetsEvshape::CalculateThrust()
{
  const AliVEvent *event = fUseMC ? MCEvent() : InputEvent();

  // to be added: track cuts for event shape observables
  const Int_t nTracksAll = event ? event->GetNumberOfTracks() : 0;
  if (nTracksAll == 0)
    return kFALSE;

  std::vector<TVector3> p;
  for (Int_t iTrack = 0; iTrack < nTracksAll; ++iTrack) {
    AliVParticle *track = event->GetTrack(iTrack);
    if (fUseMC && !MCEvent()->IsPhysicalPrimary(iTrack))
      continue;
    if (TMath::Abs(track->Eta()) > fEtaMaxForEvshape)
      continue;
    Double_t pi[3];
    track->PxPyPz(pi);
    TVector3 p3(pi);
    p.push_back(p3);
  }

  // thrust calculation following H. Yamamoto, J. Comp. Physics 52, 1983
  const Int_t nTracks = p.size();
  Char_t *signset1 = new Char_t[nTracks];
  Char_t *signset2 = new Char_t[nTracks];
  Char_t *sign = signset1;
  Char_t *signmax = signset2;
  Double_t thrustmax = 0;
  Double_t norm = 0;

  // looping over pairs of tracks (i,j) w. j > i
  TVector3 Pn;
  TVector3 Pthrust;
  TVector3 Ptemp[4];
  for (Int_t i = 0; i < nTracks; i++) {
    norm += p[i].Mag();
    for (Int_t j = i+1; j < nTracks; j++) {
      Pn = p[i].Cross(p[j]);

      // calculate hypothetical thrust
      Double_t thrust;
      Pthrust.SetXYZ(0., 0., 0.);
      for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
        if ((iTrack == i) || (iTrack == j))
          continue;
        sign[iTrack] = (Pn.Dot(p[iTrack]) > 0) ? 1 : -1;
        if (sign[iTrack] > 0)
          Pthrust += p[iTrack];
        else
          Pthrust -= p[iTrack];
      }

      Ptemp[0] = Ptemp[1] = Ptemp[2] = Ptemp[3] = Pthrust;
      Ptemp[0] += p[i]; Ptemp[0] -= p[j];
      Ptemp[1] += p[i]; Ptemp[0] += p[j];
      Ptemp[2] -= p[i]; Ptemp[0] -= p[j];
      Ptemp[3] -= p[i]; Ptemp[0] += p[j];

      Double_t thrustmaxtemp = 0;
      Int_t kmax = -1;
      for (Int_t k = 0; k < 4; k++) {
        if (Ptemp[k].Mag() > thrustmaxtemp) {
          kmax = k;
          thrustmaxtemp = Ptemp[k].Mag();
        }
      }
      if (thrustmaxtemp > thrustmax) {
        Char_t *temp = signmax;
        signmax = sign;
        sign = temp;
        thrustmax = thrustmaxtemp;
      }
    }
  }

  if (norm > 0) {
    fThrust = thrustmax / norm;
    return kTRUE;
  } else {
    return kFALSE;
  }
}

Bool_t AliAnalysisTaskJetsEvshape::CalculateSphericityT()
{
  const AliVEvent *event = fUseMC ? MCEvent() : InputEvent();

  const Int_t nTracks = event->GetNumberOfTracks();
  if (nTracks == 0)
    return kFALSE;

  // calculate transverse sphericity
  Float_t smatrix[3] = { 0. };
  Float_t sumpt = 0.;
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    AliVParticle *part= event->GetTrack(iTrack);
    if (fUseMC && !MCEvent()->IsPhysicalPrimary(iTrack))
      continue;
    if (TMath::Abs(part->Eta()) > fEtaMaxForEvshape)
      continue;
    const Float_t pt = part->Pt();
    const Float_t px = part->Px();
    const Float_t py = TMath::Sqrt(pt * pt - px * px);
    smatrix[0] += px * px / pt;
    smatrix[1] += px * py / pt;
    smatrix[2] += py * py / pt;
    sumpt += part->Pt();
  }
  if (sumpt > 0) {
    smatrix[0] /= sumpt;
    smatrix[1] /= sumpt;
    smatrix[2] /= sumpt;
  }

  // calculate eigenvalues
  Float_t lambda1 = ((smatrix[0] + smatrix[2]) + TMath::Sqrt((smatrix[0] + smatrix[2])*(smatrix[0] + smatrix[2]) - 4 * (smatrix[0] * smatrix[2] - smatrix[1] * smatrix[1]))) / 2;
  Float_t lambda2 = ((smatrix[0] + smatrix[2]) - TMath::Sqrt((smatrix[0] + smatrix[2])*(smatrix[0] + smatrix[2]) - 4 * (smatrix[0] * smatrix[2] - smatrix[1] * smatrix[1]))) / 2;

  if (TMath::Abs(lambda1 + lambda2) > 1.e-9) {
    fSphericityT = 2. * TMath::Min(lambda1, lambda2) / (lambda1 + lambda2);
    return kTRUE;
  } else {
    return kFALSE;
  }
}

// ----- histogram management -----
TH1* AliAnalysisTaskJetsEvshape::AddHistogram(Hist_t hist, const char *hid, TString title,
                                              Int_t xbins, Float_t xmin, Float_t xmax,
                                              Int_t binType)
{
  TString hName;
  hName.Form("%s_%s", fShortTaskId, hid);
  hName.ToLower();
  TH1 *h = 0x0;
  if (binType == 0)
    h = new TH1I(hName.Data(), title,
                 xbins, xmin, xmax);
  else
    h = new TH1F(hName.Data(), title,
                 xbins, xmin, xmax);
  GetHistogram(hist) = h;
  fOutputList->Add(h);
  return h;
}

TH2* AliAnalysisTaskJetsEvshape::AddHistogram(Hist_t hist, const char *hid, TString title,
                                              Int_t xbins, Float_t xmin, Float_t xmax,
                                              Int_t ybins, Float_t ymin, Float_t ymax,
                                              Int_t binType)
{
  TString hName;
  hName.Form("%s_%s", fShortTaskId, hid);
  hName.ToLower();
  TH1 *h = 0x0;
  if (binType == 0)
    h = GetHistogram(hist) = new TH2I(hName.Data(), title,
                                      xbins, xmin, xmax,
                                      ybins, ymin, ymax);
  else
    h = GetHistogram(hist) = new TH2F(hName.Data(), title,
                                      xbins, xmin, xmax,
                                      ybins, ymin, ymax);
  fOutputList->Add(h);
  return (TH2*) h;
}

TH3* AliAnalysisTaskJetsEvshape::AddHistogram(Hist_t hist, const char *hid, TString title,
                                              Int_t xbins, Float_t xmin, Float_t xmax,
                                              Int_t ybins, Float_t ymin, Float_t ymax,
                                              Int_t zbins, Float_t zmin, Float_t zmax,
                                              Int_t binType)
{
  TString hName;
  hName.Form("%s_%s", fShortTaskId, hid);
  hName.ToLower();
  TH1 *h = 0x0;
  if (binType == 0)
    h = GetHistogram(hist) = new TH3I(hName.Data(), title,
                                      xbins, xmin, xmax,
                                      ybins, ymin, ymax,
                                      zbins, zmin, zmax);
  else
    h = GetHistogram(hist) = new TH3F(hName.Data(), title,
                                      xbins, xmin, xmax,
                                      ybins, ymin, ymax,
                                      zbins, zmin, zmax);
  fOutputList->Add(h);
  return (TH3*) h;
}
