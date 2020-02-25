/*
 * AliAnalysisTaskFemtoDreamDeuteron.cxx
 *
 *  Created on: 21 Mar 2018
 *      Author: bernhardhohlweger
 */

#include "AliAnalysisTaskFemtoDreamDeuteron.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliVEvent.h"
ClassImp(AliAnalysisTaskFemtoDreamDeuteron)
AliAnalysisTaskFemtoDreamDeuteron::AliAnalysisTaskFemtoDreamDeuteron()
  : AliAnalysisTaskSE(),
    fIsMC(false),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsDeuteronDCA(nullptr),
    fTrackCutsDeuteronMass(nullptr),
    fTrackCutsAntiDeuteronDCA(nullptr),
    fTrackCutsAntiDeuteronMass(nullptr),
    fTrackCutsProtonDCA(nullptr),
    fTrackCutsProtonMass(nullptr),
    fTrackCutsAntiProtonDCA(nullptr),
    fTrackCutsAntiProtonMass(nullptr),
    fConfig(nullptr),
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
    fPairCleaner(nullptr),
    fPartColl(nullptr),
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
    fProtonBackgroundMC(nullptr),
    fAntiProtonBackgroundMC(nullptr),
    fDeuteronBackgroundMC(nullptr),
    fAntiDeuteronBackgroundMC(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
    fProtonProtonDump(nullptr),
    fProtonAntiProtonDump(nullptr),
    fProtonDeuteronDump(nullptr),
    fProtonAntiDeuteronDump(nullptr),
    fAntiProtonAntiProtonDump(nullptr),
    fAntiProtonDeuteronDump(nullptr),
    fAntiProtonAntiDeuteronDump(nullptr),
    fDeuteronAntiDeuteronDump(nullptr),
    fDumpster(nullptr),
    fUseDumpster(nullptr),
    fUseDumpsterRestPairs(nullptr),
    fTrackBufferSize() {
}

AliAnalysisTaskFemtoDreamDeuteron::AliAnalysisTaskFemtoDreamDeuteron(const char *name, bool isMC)
  : AliAnalysisTaskSE(name),
    fIsMC(isMC),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsDeuteronDCA(nullptr),
    fTrackCutsDeuteronMass(nullptr),
    fTrackCutsAntiDeuteronDCA(nullptr),
    fTrackCutsAntiDeuteronMass(nullptr),
    fTrackCutsProtonDCA(nullptr),
    fTrackCutsProtonMass(nullptr),
    fTrackCutsAntiProtonDCA(nullptr),
    fTrackCutsAntiProtonMass(nullptr),
    fConfig(nullptr),
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
    fPairCleaner(nullptr),
    fPartColl(nullptr),
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
    fProtonBackgroundMC(nullptr),
    fAntiProtonBackgroundMC(nullptr),
    fDeuteronBackgroundMC(nullptr),
    fAntiDeuteronBackgroundMC(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
    fProtonProtonDump(nullptr),
    fProtonAntiProtonDump(nullptr),
    fProtonDeuteronDump(nullptr),
    fProtonAntiDeuteronDump(nullptr),
    fAntiProtonAntiProtonDump(nullptr),
    fAntiProtonDeuteronDump(nullptr),
    fAntiProtonAntiDeuteronDump(nullptr),
    fDeuteronAntiDeuteronDump(nullptr),
    fDumpster(nullptr),
    fUseDumpster(nullptr),
    fUseDumpsterRestPairs(nullptr),
    fTrackBufferSize(2000) {
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
  DefineOutput(12, TList::Class());  //Output for the Dumpster
  if (fIsMC) {
    DefineOutput(13, TList::Class());  //Output for the Proton MC
    DefineOutput(14, TList::Class());  //Output for the AntiProton MC
    DefineOutput(15, TList::Class());  //Output for the Deuteron MC
    DefineOutput(16, TList::Class());  //Output for the AntiDeuteron MC
  }
}

AliAnalysisTaskFemtoDreamDeuteron::~AliAnalysisTaskFemtoDreamDeuteron() {
  delete fEvent;
  delete fTrack;
  delete fTrackCutsDeuteronDCA;
  delete fTrackCutsDeuteronMass;
  delete fTrackCutsAntiDeuteronDCA;
  delete fTrackCutsAntiDeuteronMass;
  delete fTrackCutsProtonDCA;
  delete fTrackCutsProtonMass;
  delete fTrackCutsAntiProtonDCA;
  delete fTrackCutsAntiProtonMass;
  delete fPairCleaner;
  delete fPartColl;

  if (fProtonProtonDump) {
    delete fProtonProtonDump;
  }
  if (fProtonAntiProtonDump) {
    delete fProtonAntiProtonDump;
  }
  if (fProtonDeuteronDump) {
    delete fProtonDeuteronDump;
  }
  if (fProtonAntiDeuteronDump) {
    delete fProtonAntiDeuteronDump;
  }
  if (fAntiProtonAntiProtonDump) {
    delete fAntiProtonAntiProtonDump;
  }
  if (fAntiProtonDeuteronDump) {
    delete fAntiProtonDeuteronDump;
  }
  if (fAntiProtonAntiDeuteronDump) {
    delete fAntiProtonAntiDeuteronDump;
  }
  if (fDeuteronAntiDeuteronDump) {
    delete fDeuteronAntiDeuteronDump;
  }
  if (fDumpster) {
    delete fDumpster;
  }
}

Float_t AliAnalysisTaskFemtoDreamDeuteron::GetMass2sq(AliFemtoDreamTrack *track) {
  Float_t p = track->GetP();
  Float_t mass2sq = -999;
  Float_t beta = track->GetbetaTOF();
  if (!(beta < 0)) {
    mass2sq = ((1 / (beta * beta)) - 1) * (p * p);
  }
  return mass2sq;
}

void AliAnalysisTaskFemtoDreamDeuteron::UserCreateOutputObjects() {

  fGTI = new AliAODTrack*[fTrackBufferSize];

  if (!fEventCuts) {
    AliError("No Event cuts \n");
  } else {
    fEventCuts->InitQA();
  }

  if (!fTrackCutsProtonDCA) {
    AliError("No Proton cuts \n");
  } else {
    fTrackCutsProtonDCA->Init();
    fProtonRestMass = new TH2F("fProtonRestMass", "Proton", 36, 0.5, 4.05, 400, 0.00, 3);
    fProtonRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fProtonRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fProtonList = fTrackCutsProtonDCA->GetQAHists();
    if (fIsMC) {
      fProtonMCList = fTrackCutsProtonDCA->GetMCQAHists();
    }
    fProtonList->Add(fProtonRestMass);
  }

  if (!fTrackCutsAntiProtonDCA) {
    AliError("No AntiProton cuts \n");
  } else {
    fTrackCutsAntiProtonDCA->Init();
    fAntiProtonRestMass = new TH2F("fAntiProtonRestMass", "AntiProton", 36, 0.5, 4.05, 400, 0.00, 3);
    fAntiProtonRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fAntiProtonRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiProtonList = fTrackCutsAntiProtonDCA->GetQAHists();
    if (fIsMC) {
      fAntiProtonMCList = fTrackCutsAntiProtonDCA->GetMCQAHists();
    }
    fAntiProtonList->Add(fAntiProtonRestMass);
  }

  if (!fTrackCutsDeuteronDCA) {
    AliError("No Proton cuts \n");
  } else {
    fTrackCutsDeuteronDCA->Init();
    fDeuteronRestMass = new TH2F("fDeuteronRestMass", "Deuteron", 72, 0.5, 8.05, 800, 0.00, 10.0);
    fDeuteronRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fDeuteronRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fDeuteronList = fTrackCutsDeuteronDCA->GetQAHists();
    if (fIsMC) {
      fDeuteronMCList = fTrackCutsDeuteronDCA->GetMCQAHists();
    }
    fDeuteronList->Add(fDeuteronRestMass);
  }

  if (!fTrackCutsAntiDeuteronDCA) {
    AliError("No Proton cuts \n");
  } else {
    fTrackCutsAntiDeuteronDCA->Init();
    fAntiDeuteronRestMass = new TH2F("fAntiDeuteronRestMass", "AntiDeuteron", 72, 0.5, 8.05, 800, 0.00, 10.0);
    fAntiDeuteronRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fAntiDeuteronRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiDeuteronList = fTrackCutsAntiDeuteronDCA->GetQAHists();
    if (fIsMC) {
      fAntiDeuteronMCList = fTrackCutsAntiDeuteronDCA->GetMCQAHists();
    }
    fAntiDeuteronList->Add(fAntiDeuteronRestMass);
  }
//--------------------------------------------------------------------------------------------------------------------
  if (!fTrackCutsProtonMass) {
    AliError("No ProtonNoTOF cuts \n");
  } else {
    fTrackCutsProtonMass->Init();
    fProtonRestMassNoTOF = new TH2F("fProtonRestMassNoTOF", "Proton NoTOF", 36, 0.5, 4.05, 400, 0.00, 3);
    fProtonRestMassNoTOF->GetXaxis()->SetTitle("pT(GeV)");
    fProtonRestMassNoTOF->GetYaxis()->SetTitle("m^2(Gev)^2");
    fProtonNoTOFList = fTrackCutsProtonMass->GetQAHists();
    fProtonNoTOFList->Add(fProtonRestMassNoTOF);
  }

  if (!fTrackCutsAntiProtonMass) {
    AliError("No AntiProtonNoTOF cuts \n");
  } else {
    fTrackCutsAntiProtonMass->Init();
    fAntiProtonRestMassNoTOF = new TH2F("fAntiProtonRestMassNoTOF", "AntiProtonNoTOF", 36, 0.5, 4.05, 400, 0.00, 3);
    fAntiProtonRestMassNoTOF->GetXaxis()->SetTitle("pT(GeV)");
    fAntiProtonRestMassNoTOF->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiProtonNoTOFList = fTrackCutsAntiProtonMass->GetQAHists();
    fAntiProtonNoTOFList->Add(fAntiProtonRestMassNoTOF);
  }

  if (!fTrackCutsDeuteronMass) {
    AliError("No DeuteronNoTOF cuts \n");
  } else {
    fTrackCutsDeuteronMass->Init();
    fDeuteronRestMassNoTOF = new TH2F("fDeuteronRestMassNoTOF", "DeuteronNoTOF", 72, 0.5, 8.05, 800, 0.00, 10.0);
    fDeuteronRestMassNoTOF->GetXaxis()->SetTitle("pT(GeV)");
    fDeuteronRestMassNoTOF->GetYaxis()->SetTitle("m^2(Gev)^2");
    fDeuteronNoTOFList = fTrackCutsDeuteronMass->GetQAHists();
    fDeuteronNoTOFList->Add(fDeuteronRestMassNoTOF);
  }

  if (!fTrackCutsAntiDeuteronMass) {
    AliError("No AntiDeuteronNoTOF cuts \n");
  } else {
    fTrackCutsAntiDeuteronMass->Init();
    fAntiDeuteronRestMassNoTOF = new TH2F("fAntiDeuteronRestMassNoTOF", "AntiDeuteronNoTOF", 72, 0.5, 8.05, 800, 0.00, 10.0);
    fAntiDeuteronRestMassNoTOF->GetXaxis()->SetTitle("pT(GeV)");
    fAntiDeuteronRestMassNoTOF->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiDeuteronNoTOFList = fTrackCutsAntiDeuteronMass->GetQAHists();
    fAntiDeuteronNoTOFList->Add(fAntiDeuteronRestMassNoTOF);
  }
//===================================================================== MC HistogramsfProtonRestMassMCfKaonRestMassMC
  if (fIsMC) {
    if (!fTrackCutsProtonMass) {
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
      fProtonBackgroundMC = new TH2F("fProtonBackgroundMC", "background", 36, 0.5, 4.05, 400, 0.0, 3);
      fProtonBackgroundMC->GetXaxis()->SetTitle("pT(GeV)");
      fProtonBackgroundMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fProtonMCList->Add(fProtonRestMassMC);
      fProtonMCList->Add(fKaonRestMassMC);
      fProtonMCList->Add(fPionRestMassMC);
      fProtonMCList->Add(fProtonBackgroundMC);
    }

    if (!fTrackCutsAntiProtonMass) {
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

    if (!fTrackCutsDeuteronMass) {
      AliError("No Proton cuts \n");
    } else {
      fDeuteronRestMassMC = new TH2F("fDeuteronRestMassMC", "Deuteron", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fDeuteronRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fDeuteronRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fDProtonRestMassMC = new TH2F("fDProtonRestMassMC", "Proton", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fDProtonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fDProtonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fDKaonRestMassMC = new TH2F("fDKaonRestMassMC", "Kaon", 72, 0.5, 8.05, 800, 0.00, 10.0);
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

    if (!fTrackCutsAntiDeuteronMass) {
      AliError("No Proton cuts \n");
    } else {
      fAntiDeuteronRestMassMC = new TH2F("fAntiDeuteronRestMassMC", "AntiDeuteron", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fAntiDeuteronRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiDeuteronRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiDProtonRestMassMC = new TH2F("fAntiDProtonRestMassMC", "AntiProton", 72, 0.5, 8.05, 800, 0.00, 10.0);
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

  fEvent = new AliFemtoDreamEvent(false, true, GetCollisionCandidates());
  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

  fDumpster = new TList();
  fDumpster->SetName("Dumpster");
  fDumpster->SetOwner(kTRUE);

  if (fUseDumpster) {
    fProtonDeuteronDump = new AliFemtoDreamDump("pd");
    fDumpster->Add(fProtonDeuteronDump->GetOutput());

    fAntiProtonAntiDeuteronDump = new AliFemtoDreamDump("ApAd");
    fDumpster->Add(fAntiProtonAntiDeuteronDump->GetOutput());
  }

  if (fUseDumpsterRestPairs) {

    fProtonProtonDump = new AliFemtoDreamDump("pp");
    fDumpster->Add(fProtonProtonDump->GetOutput());

    fProtonAntiProtonDump = new AliFemtoDreamDump("pAp");
    fDumpster->Add(fProtonAntiProtonDump->GetOutput());

    fProtonAntiDeuteronDump = new AliFemtoDreamDump("pAd");
    fDumpster->Add(fProtonAntiDeuteronDump->GetOutput());

    fAntiProtonAntiProtonDump = new AliFemtoDreamDump("ApAp");
    fDumpster->Add(fAntiProtonAntiProtonDump->GetOutput());

    fAntiProtonDeuteronDump = new AliFemtoDreamDump("Apd");
    fDumpster->Add(fAntiProtonDeuteronDump->GetOutput());

    fDeuteronAntiDeuteronDump = new AliFemtoDreamDump("dAd");
    fDumpster->Add(fDeuteronAntiDeuteronDump->GetOutput());
  }

  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
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
  PostData(12, fDumpster);

  if (fTrackCutsProtonDCA->GetIsMonteCarlo()) {
    PostData(13, fProtonMCList);
  }
  if (fTrackCutsAntiProtonDCA->GetIsMonteCarlo()) {
    PostData(14, fAntiProtonMCList);
  }

  if (fTrackCutsDeuteronDCA->GetIsMonteCarlo()) {
    PostData(15, fDeuteronMCList);
  }
  if (fTrackCutsAntiDeuteronDCA->GetIsMonteCarlo()) {
    PostData(16, fAntiDeuteronMCList);
  }
}

void AliAnalysisTaskFemtoDreamDeuteron::UserExec(Option_t *) {
  AliAODEvent *Event = static_cast<AliAODEvent*>(fInputEvent);

  if (!Event) {
    AliWarning("No Input Event");
  } else {

    fEvent->SetEvent(Event);

    if (fEventCuts->isSelected(fEvent)) {

      ResetGlobalTrackReference();

      for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
        AliAODTrack *track = static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard AOD");
          return;
        }
        StoreGlobalTrackReference(track);
      }

      fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

      static std::vector<AliFemtoDreamBasePart> DCADeuterons;
      DCADeuterons.clear();
      static std::vector<AliFemtoDreamBasePart> MassDeuterons;
      MassDeuterons.clear();
      static std::vector<AliFemtoDreamBasePart> DCAAntiDeuterons;
      DCAAntiDeuterons.clear();
      static std::vector<AliFemtoDreamBasePart> MassAntiDeuterons;
      MassAntiDeuterons.clear();

      static std::vector<AliFemtoDreamBasePart> DCAProtons;
      DCAProtons.clear();
      static std::vector<AliFemtoDreamBasePart> MassProtons;
      MassProtons.clear();
      static std::vector<AliFemtoDreamBasePart> DCAAntiProtons;
      DCAAntiProtons.clear();
      static std::vector<AliFemtoDreamBasePart> MassAntiProtons;
      MassAntiProtons.clear();
      //Now we loop over all the tracks in the reconstructed event.
      for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
        AliAODTrack *track = static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard AOD");
          return;
        }

        fTrack->SetTrack(track);

        if (fTrackCutsDeuteronDCA->isSelected(fTrack)) {
          DCADeuterons.push_back(*fTrack);
        }
        if (fTrackCutsDeuteronMass->isSelected(fTrack)) {
          MassDeuterons.push_back(*fTrack);
          fDeuteronRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
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

        if (fTrackCutsAntiDeuteronDCA->isSelected(fTrack)) {
          //.. we add it to our particle buffer
          DCAAntiDeuterons.push_back(*fTrack);
        }

        if (fTrackCutsAntiDeuteronMass->isSelected(fTrack)) {
          MassAntiDeuterons.push_back(*fTrack);
          fAntiDeuteronRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
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

        if (fTrackCutsProtonDCA->isSelected(fTrack)) {
          //.. we add it to our particle buffer
          fTrack->SetCPA(gRandom->Uniform());
          DCAProtons.push_back(*fTrack);
        }
        if (fTrackCutsProtonMass->isSelected(fTrack)) {
          MassProtons.push_back(*fTrack);
          fProtonRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
          if (fIsMC) {

            if (fTrack->GetMCPDGCode() == 2212) {

              fProtonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            } else if (fTrack->GetMCPDGCode() == 321 ) {

              fKaonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            } else if (fTrack->GetMCPDGCode() == 211 ) {

              fPionRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            } else {
              fProtonBackgroundMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
            }
          }
        }

        if (fTrackCutsAntiProtonDCA->isSelected(fTrack)) {
          //.. we add it to our particle buffer
          DCAAntiProtons.push_back(*fTrack);
        }

        if (fTrackCutsAntiProtonMass->isSelected(fTrack)) {
          MassAntiProtons.push_back(*fTrack);
          fAntiProtonRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
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
      }

      fPairCleaner->CleanTrackAndDecay(&DCAProtons, &DCADeuterons, 0);
      fPairCleaner->CleanTrackAndDecay(&DCAAntiProtons, &DCAAntiDeuterons, 1);
      fPairCleaner->ResetArray();
      fPairCleaner->StoreParticle(DCAProtons);
      fPairCleaner->StoreParticle(DCAAntiProtons);
      fPairCleaner->StoreParticle(DCADeuterons);
      fPairCleaner->StoreParticle(DCAAntiDeuterons);
      fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                          fEvent->GetRefMult08(), fEvent->GetV0MCentrality());
      void SetEvent(std::vector<AliFemtoDreamBasePart> &vec1,
                    std::vector<AliFemtoDreamBasePart> &vec2,
                    AliFemtoDreamEvent * evt, const int pdg1, const int pdg2);


      if (fUseDumpster) {
        if (fProtonDeuteronDump) {
          fProtonDeuteronDump->SetEvent(DCAProtons, DCADeuterons, fEvent, 2212, 1000010020);
        }

        if (fAntiProtonAntiDeuteronDump) {
          fAntiProtonAntiDeuteronDump->SetEvent(DCAAntiProtons, DCAAntiDeuterons, fEvent, -2212, -1000010020);
        }
      }
      if (fUseDumpsterRestPairs) {
        if (fProtonProtonDump) {
          fProtonProtonDump->SetEvent(DCAProtons, DCAProtons, fEvent, 2212, 2212);
        }
        if (fProtonAntiProtonDump) {
          fProtonAntiProtonDump->SetEvent(DCAProtons, DCAAntiProtons, fEvent, 2212, -2212);
        }
        if (fProtonAntiDeuteronDump) {
          fProtonAntiDeuteronDump->SetEvent(DCAProtons, DCAAntiDeuterons, fEvent, 2212, -1000010020);
        }
        if (fAntiProtonAntiProtonDump) {
          fAntiProtonAntiProtonDump->SetEvent(DCAAntiProtons, DCAAntiProtons, fEvent, -2212, -2212);
        }
        if (fAntiProtonDeuteronDump) {
          fAntiProtonDeuteronDump->SetEvent(DCAAntiProtons, DCADeuterons, fEvent, -2212, 1000010020);
        }
        if (fDeuteronAntiDeuteronDump) {
          fDeuteronAntiDeuteronDump->SetEvent(DCADeuterons, DCAAntiDeuterons, fEvent, 1000010020, -1000010020);
        }
      }
    }
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
  PostData(12, fDumpster);
//-----------------------------------------MCTracksStorage------------------------------------------------------------------------------
  if (fTrackCutsProtonDCA->GetIsMonteCarlo()) {
    PostData(13, fProtonMCList);
  }
  if (fTrackCutsAntiProtonDCA->GetIsMonteCarlo()) {
    PostData(14, fAntiProtonMCList);
  }
  if (fTrackCutsDeuteronDCA->GetIsMonteCarlo()) {
    PostData(15, fDeuteronMCList);
  }
  if (fTrackCutsAntiDeuteronDCA->GetIsMonteCarlo()) {
    PostData(16, fAntiDeuteronMCList);
  }
}

void AliAnalysisTaskFemtoDreamDeuteron::ResetGlobalTrackReference() {

  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

void AliAnalysisTaskFemtoDreamDeuteron::StoreGlobalTrackReference(AliAODTrack *track) {

  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }

  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n"
           , trackID, fTrackBufferSize);
    return;
  }

  if (fGTI[trackID]) {

    if ( (!track->GetFilterMap()) && (!track->GetTPCNcls()) ) {
      return;
    }

    if ( fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls()  ) {
      // If we come here, there's a problem
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}
