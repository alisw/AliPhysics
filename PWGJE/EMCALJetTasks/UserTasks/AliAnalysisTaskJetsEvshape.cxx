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
  AddHistogram(ID(kHistMult), "tracklet multiplicity;N_{trkl};counts",
               100, 0., 400.);
  AddHistogram(ID(kHistJetPtVsMult), "jet p_{T} vs mult;multiplicity;p_{T}^{jet,ch} (GeV/#it{c})",
               100, 0., 400.,
               40, 0., 40.);

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
  return kTRUE;
}

Bool_t AliAnalysisTaskJetsEvshape::FillHistograms()
{
  FillH1(kHistStat, kStatUsed);

  const AliStack *stack = fMCEvent->Stack();
  const AliVMultiplicity *mult = InputEvent()->GetMultiplicity();
  const Int_t nTracklets =
    stack ? stack->GetNprimary() :
    mult ? mult->GetNumberOfTracklets() :
    -1;
  FillH1(kHistMult, nTracklets);

  if (fJetsCont) {
    // const Int_t nJets         = fJetsCont->GetNJets();
    // const Int_t nAcceptedJets = fJetsCont->GetNAcceptedJets();

    fJetsCont->ResetCurrentID();
    while (AliEmcalJet *jet = fJetsCont->GetNextAcceptJet()) {
      FillH1(kHistJetPt, jet->Pt());
      FillH1(kHistJetPtVsMult, nTracklets, jet->Pt());
    }
  }

  return kTRUE;
}

void AliAnalysisTaskJetsEvshape::UserExec(Option_t *option)
{
  // actual work

  AliAnalysisTaskEmcalJet::UserExec(option);

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

  // record number of sampled events and detect trigger contributions
  FillH1(kHistStat, kStatSeen);

  // so far, no trigger selection, we accept all
  FillH1(kHistStat, kStatTrg);

  // prepare the event
  // (make sure it is cleaned up in the end)
  if (PrepareEvent()) {
    // here we have passed the event cuts
    FillH1(kHistStat, kStatEvCuts);

    // multiplicity selection

    // event shape selection

    // InputEvent()->GetList()->ls();
  }

 stop:
  CleanUpEvent();

  PostData(kOutputEmcal, fOutput);
  PostData(kOutputTask, fOutputList);
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
