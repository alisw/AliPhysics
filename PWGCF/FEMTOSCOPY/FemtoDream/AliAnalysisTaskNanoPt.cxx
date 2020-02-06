#include "AliAnalysisTaskNanoPt.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliNanoAODTrack.h"
#include "AliPIDResponse.h"
ClassImp(AliAnalysisTaskNanoPt)

//--------------------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskNanoPt::AliAnalysisTaskNanoPt()
  : AliAnalysisTaskSE(),
    fInputEvent(nullptr),
    fEvent(nullptr),
    fEvtCuts(nullptr),
    fTrack(nullptr),
    fProtonTrack(nullptr),
    fAntiProtonTrack(nullptr),
    fDeuteronTrack(nullptr),
    fAntiDeuteronTrack(nullptr),
    fProtonTrackNoTOF(nullptr),
    fAntiProtonTrackNoTOF(nullptr),
    fDeuteronTrackNoTOF(nullptr),
    fAntiDeuteronTrackNoTOF(nullptr),
    fConfig(nullptr),
    fIsMC(false),
    fEvtList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fDeuteronList(nullptr),
    fDeuteronMCList(nullptr),
    fAntiDeuteronList(nullptr),
    fAntiDeuteronMCList(nullptr),
    fProtonNoTOFList(nullptr),
    fProtonMCNoTOFList(nullptr),
    fAntiProtonNoTOFList(nullptr),
    fAntiProtonMCNoTOFList(nullptr),
    fDeuteronNoTOFList(nullptr),
    fDeuteronMCNoTOFList(nullptr),
    fAntiDeuteronNoTOFList(nullptr),
    fAntiDeuteronMCNoTOFList(nullptr),
    fGTI(nullptr),
    fProtonRestMass(nullptr),
    fAntiProtonRestMass(nullptr),
    fDeuteronRestMass(nullptr),
    fAntiDeuteronRestMass(nullptr),
    fProtonRestMassNoTOF(nullptr),
    fAntiProtonRestMassNoTOF(nullptr),
    fDeuteronRestMassNoTOF(nullptr),
    fAntiDeuteronRestMassNoTOF(nullptr),
    fProtonRestMassMC(nullptr),
    fAntiProtonRestMassMC(nullptr),
    fDeuteronRestMassMC(nullptr),
    fAntiDeuteronRestMassMC(nullptr),
    fKaonRestMassMC(nullptr),
    fAntiKaonRestMassMC(nullptr),
    fDProtonRestMassMC(nullptr),
    fDKaonRestMassMC(nullptr),
    fAntiDProtonRestMassMC (nullptr),
    fAntiDKaonRestMassMC(nullptr),
    fPionRestMassMC(nullptr),
    fAntiPionRestMassMC(nullptr),
    fDPionRestMassMC(nullptr),
    fAntiDPionRestMassMC(nullptr),
    fProtonBackgroungMC(nullptr),
    fAntiProtonBackgroundMC(nullptr),
    fDeuteronBackgroundMC(nullptr),
    fAntiDeuteronBackgroundMC(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
    fTrackBufferSize(2500) {}
//-----------------------------------------------------------------------------------------------------------------------

AliAnalysisTaskNanoPt::AliAnalysisTaskNanoPt(
  const char *name, const bool isMC)
  : AliAnalysisTaskSE(name),
    fInputEvent(nullptr),
    fEvent(nullptr),
    fEvtCuts(nullptr),
    fTrack(nullptr),
    fProtonTrack(nullptr),
    fAntiProtonTrack(nullptr),
    fDeuteronTrack(nullptr),
    fAntiDeuteronTrack(nullptr),
    fProtonTrackNoTOF(nullptr),
    fAntiProtonTrackNoTOF(nullptr),
    fDeuteronTrackNoTOF(nullptr),
    fAntiDeuteronTrackNoTOF(nullptr),
    fConfig(nullptr),
    fIsMC(isMC),
    fEvtList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fDeuteronList(nullptr),
    fDeuteronMCList(nullptr),
    fAntiDeuteronList(nullptr),
    fAntiDeuteronMCList(nullptr),
    fProtonNoTOFList(nullptr),
    fProtonMCNoTOFList(nullptr),
    fAntiProtonNoTOFList(nullptr),
    fAntiProtonMCNoTOFList(nullptr),
    fDeuteronNoTOFList(nullptr),
    fDeuteronMCNoTOFList(nullptr),
    fAntiDeuteronNoTOFList(nullptr),
    fAntiDeuteronMCNoTOFList(nullptr),
    fGTI(nullptr),
    fProtonRestMass(nullptr),
    fAntiProtonRestMass(nullptr),
    fDeuteronRestMass(nullptr),
    fAntiDeuteronRestMass(nullptr),
    fProtonRestMassNoTOF(nullptr),
    fAntiProtonRestMassNoTOF(nullptr),
    fDeuteronRestMassNoTOF(nullptr),
    fAntiDeuteronRestMassNoTOF(nullptr),
    fProtonRestMassMC(nullptr),
    fAntiProtonRestMassMC(nullptr),
    fDeuteronRestMassMC(nullptr),
    fAntiDeuteronRestMassMC(nullptr),
    fKaonRestMassMC(nullptr),
    fAntiKaonRestMassMC(nullptr),
    fDProtonRestMassMC(nullptr),
    fDKaonRestMassMC(nullptr),
    fAntiDProtonRestMassMC (nullptr),
    fAntiDKaonRestMassMC(nullptr),
    fPionRestMassMC(nullptr),
    fAntiPionRestMassMC(nullptr),
    fDPionRestMassMC(nullptr),
    fAntiDPionRestMassMC(nullptr),
    fProtonBackgroungMC(nullptr),
    fAntiProtonBackgroundMC(nullptr),
    fDeuteronBackgroundMC(nullptr),
    fAntiDeuteronBackgroundMC(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
    fTrackBufferSize(2500) {
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Proton Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(4, TList::Class());  //Output for the Dueteron Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiDeuteron Cuts
  DefineOutput(6, TList::Class());  //Output for the ProtonNoTOF Cuts
  DefineOutput(7, TList::Class());  //Output for the AntiProtonNoTOF Cuts
  DefineOutput(8, TList::Class());  //Output for the DueteronNoTOF Cuts
  DefineOutput(9, TList::Class());  //Output for the AntiDeuteronNoTOF Cuts
  DefineOutput(10, TList::Class());  //Output for the Results
  DefineOutput(11, TList::Class());  //Output for the Results QA
  if (fIsMC) {
    DefineOutput(12, TList::Class());  //Output for the Proton MC
    DefineOutput(13, TList::Class());  //Output for the AntiProton MC
    DefineOutput(14, TList::Class());  //Output for the Deuteron MC
    DefineOutput(15, TList::Class());  //Output for the AntiDeuteron MC
  }
}

//---------------------------------------------------------------------------------------------------------------------------

AliAnalysisTaskNanoPt::~AliAnalysisTaskNanoPt() {
  delete fEvent;
  delete fTrack;
  delete fProtonTrack;
  delete fAntiProtonTrack;
  delete fDeuteronTrack;
  delete fAntiDeuteronTrack;
  delete fProtonTrackNoTOF;
  delete fAntiProtonTrackNoTOF;
  delete fDeuteronTrackNoTOF;
  delete fAntiDeuteronTrackNoTOF;
  delete fPairCleaner;
  delete fPartColl;
}
//-----------------------------------------------------------------------------------------------------------------

Float_t AliAnalysisTaskNanoPt::GetMass2sq(AliFemtoDreamTrack *track) const {
  Float_t p = track->GetP();
  Float_t mass2sq = -999;
  Float_t beta = track->GetbetaTOF();
  if (beta > 0) {
    mass2sq = ((1 / (beta * beta)) - 1) * (p * p);
  }
  return mass2sq;
}

//-------------------------------------------------------UserCreateOutPut----------------------------------------------------------------------------------

void AliAnalysisTaskNanoPt::UserCreateOutputObjects() {

  fGTI = new AliVTrack *[fTrackBufferSize];

  if (!fEvtCuts) {
    AliError("No Event cuts \n");
  } else {
    fEvtCuts->InitQA();
  }

  if (!fProtonTrack) {
    AliError("No Proton cuts \n");
  } else {
    fProtonTrack->Init();
    fProtonRestMass = new TH2F("fProtonRestMass", "Proton", 36, 0.5, 4.05, 400, 0.00, 3);
    fProtonRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fProtonRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fProtonList = fProtonTrack->GetQAHists();
    if (fIsMC) {
      fProtonMCList = fProtonTrack->GetMCQAHists();
    }
    fProtonList->Add(fProtonRestMass);
  }

  if (!fAntiProtonTrack) {
    AliError("No AntiProton cuts \n");
  } else {
    fAntiProtonTrack->Init();
    fAntiProtonRestMass = new TH2F("fAntiProtonRestMass", "AntiProton", 36, 0.5, 4.05, 400, 0.00, 3);
    fAntiProtonRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fAntiProtonRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiProtonList = fAntiProtonTrack->GetQAHists();
    if (fIsMC) {
      fAntiProtonMCList = fAntiProtonTrack->GetMCQAHists();
    }
    fAntiProtonList->Add(fAntiProtonRestMass);
  }

  if (!fDeuteronTrack) {
    AliError("No Proton cuts \n");
  } else {
    fDeuteronTrack->Init();
    fDeuteronRestMass = new TH2F("fDeuteronRestMass", "Deuteron", 72, 0.5, 8.05, 800, 0.00, 10.0);
    fDeuteronRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fDeuteronRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fDeuteronList = fDeuteronTrack->GetQAHists();
    if (fIsMC) {
      fDeuteronMCList = fDeuteronTrack->GetMCQAHists();
    }
    fDeuteronList->Add(fDeuteronRestMass);
  }

  if (!fAntiDeuteronTrack) {
    AliError("No Proton cuts \n");
  } else {
    fAntiDeuteronTrack->Init();
    fAntiDeuteronRestMass = new TH2F("fAntiDeuteronRestMass", "AntiDeuteron", 72, 0.5, 8.05, 800, 0.00, 10.0);
    fAntiDeuteronRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fAntiDeuteronRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiDeuteronList = fAntiDeuteronTrack->GetQAHists();
    if (fIsMC) {
      fAntiDeuteronMCList = fAntiDeuteronTrack->GetMCQAHists();
    }
    fAntiDeuteronList->Add(fAntiDeuteronRestMass);
  }
//--------------------------------------------------------------------------------------------------------------------
  if (!fProtonTrackNoTOF) {
    AliError("No ProtonNoTOF cuts \n");
  } else {
    fProtonTrackNoTOF->Init();
    fProtonRestMassNoTOF = new TH2F("fProtonRestMassNoTOF", "Proton NoTOF", 36, 0.5, 4.05, 400, 0.00, 3);
    fProtonRestMassNoTOF->GetXaxis()->SetTitle("pT(GeV)");
    fProtonRestMassNoTOF->GetYaxis()->SetTitle("m^2(Gev)^2");
    fProtonNoTOFList = fProtonTrackNoTOF->GetQAHists();
    fProtonNoTOFList->Add(fProtonRestMassNoTOF);
  }

  if (!fAntiProtonTrackNoTOF) {
    AliError("No AntiProtonNoTOF cuts \n");
  } else {
    fAntiProtonTrackNoTOF->Init();
    fAntiProtonRestMassNoTOF = new TH2F("fAntiProtonRestMassNoTOF", "AntiProtonNoTOF", 36, 0.5, 4.05, 400, 0.00, 3);
    fAntiProtonRestMassNoTOF->GetXaxis()->SetTitle("pT(GeV)");
    fAntiProtonRestMassNoTOF->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiProtonNoTOFList = fAntiProtonTrackNoTOF->GetQAHists();
    fAntiProtonNoTOFList->Add(fAntiProtonRestMassNoTOF);
  }

  if (!fDeuteronTrackNoTOF) {
    AliError("No DeuteronNoTOF cuts \n");
  } else {
    fDeuteronTrackNoTOF->Init();
    fDeuteronRestMassNoTOF = new TH2F("fDeuteronRestMassNoTOF", "DeuteronNoTOF", 72, 0.5, 8.05, 800, 0.00, 10.0);
    fDeuteronRestMassNoTOF->GetXaxis()->SetTitle("pT(GeV)");
    fDeuteronRestMassNoTOF->GetYaxis()->SetTitle("m^2(Gev)^2");
    fDeuteronNoTOFList = fDeuteronTrackNoTOF->GetQAHists();
    fDeuteronNoTOFList->Add(fDeuteronRestMassNoTOF);
  }

  if (!fAntiDeuteronTrackNoTOF) {
    AliError("No AntiDeuteronNoTOF cuts \n");
  } else {
    fAntiDeuteronTrackNoTOF->Init();
    fAntiDeuteronRestMassNoTOF = new TH2F("fAntiDeuteronRestMassNoTOF", "AntiDeuteronNoTOF", 72, 0.5, 8.05, 800, 0.00, 10.0);
    fAntiDeuteronRestMassNoTOF->GetXaxis()->SetTitle("pT(GeV)");
    fAntiDeuteronRestMassNoTOF->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiDeuteronNoTOFList = fAntiDeuteronTrackNoTOF->GetQAHists();
    fAntiDeuteronNoTOFList->Add(fAntiDeuteronRestMassNoTOF);

  }
//===================================================================== MC HistogramsfProtonRestMassMCfKaonRestMassMC
  if (fIsMC) {
    if (!fProtonTrackNoTOF ) {
      AliError("No Proton cuts \n");
    } else {
      fProtonRestMassMC = new TH2F("fProtonRestMassMC", "Proton", 36, 0.5, 4.05, 400, 0.0, 3);
      fProtonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fProtonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fKaonRestMassMC = new TH2F("fKaonRestMassMC", "kaon", 36, 0.5, 4.05, 400, 0.0, 3);
      fKaonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fKaonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fPionRestMassMC = new TH2F("fPionRestMassMC", "pion", 36, 0.5, 4.05, 400, 0.0, 3);
      fPionRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fPionRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fProtonBackgroungMC = new TH2F("fProtonBackgroungMC", "background", 36, 0.5, 4.05, 400, 0.0, 3);
      fProtonBackgroungMC->GetXaxis()->SetTitle("pT(GeV)");
      fProtonBackgroungMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fProtonMCList->Add(fProtonRestMassMC);
      fProtonMCList->Add(fKaonRestMassMC);
      fProtonMCList->Add(fPionRestMassMC);
      fProtonMCList->Add(fProtonBackgroungMC);
    }

    if (!fAntiProtonTrackNoTOF ) {
      AliError("No AntiProton cuts \n");
    } else {
      fAntiProtonRestMassMC = new TH2F("fAntiProtonRestMassMC", "AntiProton", 36, 0.5, 4.05, 400, 0.00, 3);
      fAntiProtonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiProtonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiKaonRestMassMC = new TH2F("fAntiKaonRestMassMC", "AntiProton", 36, 0.5, 4.05, 400, 0.00, 3);
      fAntiKaonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiKaonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiPionRestMassMC = new TH2F("fAntiPionRestMassMC", "AntiPion", 36, 0.5, 4.05, 400, 0.00, 3);
      fAntiPionRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiPionRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiProtonBackgroundMC = new TH2F("fAntiProtonBackgroundMC", "AntiProtonBackgroundMC", 36, 0.5, 4.05, 400, 0.00, 3);
      fAntiProtonBackgroundMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiProtonBackgroundMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiProtonMCList->Add(fAntiProtonRestMassMC);
      fAntiProtonMCList->Add(fAntiKaonRestMassMC);
      fAntiProtonMCList->Add(fAntiPionRestMassMC);
      fAntiProtonMCList->Add(fAntiProtonBackgroundMC);
    }

    if (!fDeuteronTrackNoTOF ) {
      AliError("No Proton cuts \n");
    } else {
      fDeuteronRestMassMC = new TH2F("fDeuteronRestMassMC", "Deuteron", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fDeuteronRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fDeuteronRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fDProtonRestMassMC = new TH2F("fDProtonRestMassMC", "Proton", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fDProtonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fDProtonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fDKaonRestMassMC = new TH2F("fDKaonRestMassMC", "Kaon",72, 0.5, 8.05, 800, 0.00, 10.0);
      fDKaonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fDKaonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fDPionRestMassMC = new TH2F("fDPionRestMassMC", "Pion", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fDPionRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fDPionRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fDeuteronBackgroundMC = new TH2F("fDeuteronBackgroundMC", "Pion", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fDeuteronBackgroundMC->GetXaxis()->SetTitle("pT(GeV)");
      fDeuteronBackgroundMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fDeuteronMCList->Add(fDeuteronRestMassMC);
      fDeuteronMCList->Add(fDProtonRestMassMC);
      fDeuteronMCList->Add(fDKaonRestMassMC);
      fDeuteronMCList->Add(fDPionRestMassMC);
      fDeuteronMCList->Add(fDeuteronBackgroundMC);
    }

    if (!fAntiDeuteronTrackNoTOF) {
      AliError("No Proton cuts \n");
    } else {
      fAntiDeuteronRestMassMC = new TH2F("fAntiDeuteronRestMassMC", "AntiDeuteron", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fAntiDeuteronRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiDeuteronRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiDProtonRestMassMC = new TH2F("fAntiDProtonRestMassMC", "AntiProton",72, 0.5, 8.05, 800, 0.00, 10.0);
      fAntiDProtonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiDProtonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiDKaonRestMassMC = new TH2F("fAntiDKaonRestMassMC", "AntiKaon", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fAntiDKaonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiDKaonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiDPionRestMassMC = new TH2F("fAntiDPionRestMassMC", "AntiPion", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fAntiDPionRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiDPionRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiDeuteronBackgroundMC = new TH2F("fAntiDeuteronBackgroundMC", "AntiDeuteronBackgroundMC", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fAntiDeuteronBackgroundMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiDeuteronBackgroundMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiDeuteronMCList->Add(fAntiDeuteronRestMassMC);
      fAntiDeuteronMCList->Add(fAntiDProtonRestMassMC);
      fAntiDeuteronMCList->Add(fAntiDKaonRestMassMC);
      fAntiDeuteronMCList->Add(fAntiDPionRestMassMC);
      fAntiDeuteronMCList->Add(fAntiDeuteronBackgroundMC);
    }

  }

  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
        fConfig->GetMinimalBookingME());
    fPairCleaner = new AliFemtoDreamPairCleaner(2, 0,
        fConfig->GetMinimalBookingME());
  }

  fEvent = new AliFemtoDreamEvent(false, true, GetCollisionCandidates(), false);
  fTrack = new AliFemtoDreamTrack();

  fTrack->SetUseMCInfo(fIsMC);

  if (!fEvtCuts->GetMinimalBooking()) {
    fEvtList = fEvtCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();

  }

  fResultsQA = new TList();
  fResultsQA->SetOwner();
  fResultsQA->SetName("ResultsQA");

  if (fConfig->GetUseEventMixing()) {
    fResults = fPartColl->GetHistList();
    if (!fConfig->GetMinimalBookingME()) {
      fResultsQA->Add(fPartColl->GetQAList());
      fResultsQA->Add(fPairCleaner->GetHistList());
    }
  } else {
    fResults = new TList();
    fResults->SetOwner();
    fResults->SetName("Results");
  }

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fDeuteronList);
  PostData(5, fAntiDeuteronList);
  PostData(6, fProtonNoTOFList);
  PostData(7, fAntiProtonNoTOFList);
  PostData(8, fDeuteronNoTOFList);
  PostData(9, fAntiDeuteronNoTOFList);
  PostData(10, fResults);
  PostData(11, fResultsQA);
  if (fProtonTrack->GetIsMonteCarlo()) {
    PostData(12, fProtonMCList);
  }
  if (fAntiProtonTrack->GetIsMonteCarlo()) {
    PostData(13, fAntiProtonMCList);
  }

  if (fDeuteronTrack->GetIsMonteCarlo()) {
    PostData(14, fDeuteronMCList);
  }
  if (fAntiDeuteronTrack->GetIsMonteCarlo()) {
    PostData(15, fAntiDeuteronMCList);
  }

}

//------------------------------------------UserExec()----------------------------------------------------------------------------


void AliAnalysisTaskNanoPt::UserExec(Option_t  *option ) {
  AliVEvent *fInputEvent = InputEvent();


  // PREAMBLE - CHECK EVERYTHING IS THERE
  if (!fInputEvent) {
    AliError("No Input event");
    return;
  }

  if (fIsMC && !fMCEvent) {
    AliError("No MC event");
    return;
  }


  if (!fEvtCuts) {
    AliError("Event Cuts missing");
    return;
  }

  if (!fProtonTrack || !fAntiProtonTrack) {
    AliError("Proton Cuts missing");
    return;
  }

  if (!fDeuteronTrack || !fAntiDeuteronTrack) {
    AliError("Deuteron Cuts missing");
    return;
  }


  if (!fProtonTrackNoTOF || !fAntiProtonTrackNoTOF) {
    AliError("ProtonNoTOF Cuts missing");
    return;
  }

  if (!fDeuteronTrackNoTOF || !fAntiDeuteronTrackNoTOF) {
    AliError("DeuteronNoTOF Cuts missing");
    return;
  }

  // EVENT SELECTION
  fEvent->SetEvent(fInputEvent);
  if (!fEvtCuts->isSelected(fEvent))
    return;

  // PROTON SELECTION
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }

    StoreGlobalTrackReference(track);
  }

  std::vector<AliFemtoDreamBasePart> Proton;
  std::vector<AliFemtoDreamBasePart> AntiProton;
  std::vector<AliFemtoDreamBasePart> Deuteron;
  std::vector<AliFemtoDreamBasePart> AntiDeuteron;
  std::vector<AliFemtoDreamBasePart> ProtonNoTOF;
  std::vector<AliFemtoDreamBasePart> AntiProtonNoTOF;
  std::vector<AliFemtoDreamBasePart> DeuteronNoTOF;
  std::vector<AliFemtoDreamBasePart> AntiDeuteronNoTOF;

  const int multiplicity = fEvent->GetMultiplicity();
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fTrack->SetTrack(track, fInputEvent, multiplicity);

    if (fProtonTrack->isSelected(fTrack)) {

      fProtonRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      Proton.push_back(*fTrack);
    }
    if (fAntiProtonTrack->isSelected(fTrack)) {
      fAntiProtonRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      AntiProton.push_back(*fTrack);
    }
    if (fDeuteronTrack->isSelected(fTrack)) {
      fDeuteronRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      Deuteron.push_back(*fTrack);
    }
    if (fAntiDeuteronTrack->isSelected(fTrack)) {
      fAntiDeuteronRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      AntiDeuteron.push_back(*fTrack);
    }

    if (fProtonTrackNoTOF->isSelected(fTrack)) {
      fProtonRestMassNoTOF->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      ProtonNoTOF.push_back(*fTrack);

      if (fIsMC) {

        if (fTrack->GetMCPDGCode() == 2212) {

          fProtonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else if (fTrack->GetMCPDGCode() == 321 ) {

          fKaonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else if (fTrack->GetMCPDGCode() == 211 ) {

          fPionRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else {
          fProtonBackgroungMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
        }
      }
    }
    if (fAntiProtonTrackNoTOF->isSelected(fTrack)) {
      fAntiProtonRestMassNoTOF->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      AntiProtonNoTOF.push_back(*fTrack);

      if (fIsMC) {
        if (fTrack->GetMCPDGCode() == -2212) {

          fAntiProtonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else if (fTrack->GetMCPDGCode() == -321 ) {

          fAntiKaonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else if (fTrack->GetMCPDGCode() == -211 ) {

          fAntiPionRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else {
          fAntiProtonBackgroundMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        }
      }

    }
    if (fDeuteronTrackNoTOF->isSelected(fTrack)) {

      fDeuteronRestMassNoTOF->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      DeuteronNoTOF.push_back(*fTrack);

      if (fIsMC) {

        if (fTrack->GetMCPDGCode() == 1000010020) {

          fDeuteronRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else if (fTrack->GetMCPDGCode() == 2212) {

          fDProtonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else if (fTrack->GetMCPDGCode() == 321 ) {

          fDKaonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else if (fTrack->GetMCPDGCode() == 211 ) {

          fDPionRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else {

          fDeuteronBackgroundMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
        }
      }

    }

    if (fAntiDeuteronTrackNoTOF->isSelected(fTrack)) {
      fAntiDeuteronRestMassNoTOF->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      AntiDeuteronNoTOF.push_back(*fTrack);

      if (fIsMC) {

        if (fTrack->GetMCPDGCode() == -1000010020) {

          fAntiDeuteronRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else if (fTrack->GetMCPDGCode() == -2212) {

          fAntiDProtonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else if (fTrack->GetMCPDGCode() == -321 ) {

          fAntiDKaonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else if (fTrack->GetMCPDGCode() == -211 ) {

          fAntiDPionRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

        } else {
          fAntiDeuteronBackgroundMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
        }
      }

    }

  }

  fPairCleaner->CleanTrackAndDecay(&Proton, &Deuteron, 0);
  fPairCleaner->CleanTrackAndDecay(&AntiProton, &AntiDeuteron, 1);
  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(Proton);
  fPairCleaner->StoreParticle(AntiProton);
  fPairCleaner->StoreParticle(Deuteron);
  fPairCleaner->StoreParticle(AntiDeuteron);
  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());

// flush the data
  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fDeuteronList);
  PostData(5, fAntiDeuteronList);
  PostData(6, fProtonNoTOFList);
  PostData(7, fAntiProtonNoTOFList);
  PostData(8, fDeuteronNoTOFList);
  PostData(9, fAntiDeuteronNoTOFList);
  PostData(10, fResults);
  PostData(11, fResultsQA);
//-----------------------------------------MCTracksStorage------------------------------------------------------------------------------
  if (fProtonTrack->GetIsMonteCarlo()) {
    PostData(12, fProtonMCList);
  }
  if (fAntiProtonTrack->GetIsMonteCarlo()) {
    PostData(13, fAntiProtonMCList);
  }

  if (fDeuteronTrack->GetIsMonteCarlo()) {
    PostData(14, fDeuteronMCList);
  }
  if (fAntiDeuteronTrack->GetIsMonteCarlo()) {
    PostData(15, fAntiDeuteronMCList);
  }
}


///------------------------------------------------------------------------
void AliAnalysisTaskNanoPt::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//-------------------------------------------------------------------------

void AliAnalysisTaskNanoPt::StoreGlobalTrackReference(AliVTrack * track) {
  // see AliFemtoDreamAnalysis for details
  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack *>(track);
  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }
  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }
  if (fGTI[trackID]) {

    if ((!nanoTrack->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    if (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap()
        || fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap(),
             nanoTrack->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}


