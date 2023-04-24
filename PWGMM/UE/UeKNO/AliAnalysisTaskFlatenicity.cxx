/****************************************************************************
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
 *                                                                        *
 * Author: Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)                *
 * Anatask to compute flatenicity (arXiv:2204.13733)                      *
 **************************************************************************/

class TTree;

class AliESDtrackCuts;
class AliESDAD; // AD

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliESDAD.h" //AD
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDUtils.h"
#include "AliESDVZERO.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliVZDC.h" //AD
#include <AliVAD.h>  //AD

#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultVariable.h"
#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "TParticle.h"
#include "TProfile.h"
#include "TVector3.h"
#include <AliAnalysisFilter.h>
#include <AliESDVertex.h>
#include <AliHeader.h>
#include <AliMultiplicity.h>
#include <Riostream.h>
#include <TBits.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TRandom.h>
#include <TTree.h>
using std::cout;
using std::endl;

#include "AliAnalysisTaskFlatenicity.h"

const Int_t nPtbins = 36;
Double_t Ptbins[nPtbins + 1] = {
    0.0, 0.1, 0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5,  0.6, 0.7, 0.8,
    0.9, 1.0, 1.25, 1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5,  5.0, 6.0, 7.0,
    8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 30.0, 40.0, 50.0};
const Int_t nCent = 9;
Double_t centClass[nCent + 1] = {0.0,  1.0,  5.0,  10.0, 20.0,
                                 30.0, 40.0, 50.0, 70.0, 100.0};

const Int_t nFlatClass = 8;
Double_t FlatClass[nFlatClass + 1] = {0.0,  0.01, 0.05, 0.10, 0.20,
                                      0.30, 0.40, 0.50, 1.0};

// const Int_t nCent = 2;
// Double_t centClass[nCent + 1] = {0.0, 50.0, 100.0};
const Int_t nDet = 4;
const Char_t *DetName[nDet] = {"ADC", "V0C", "V0A", "ADA"};
const Int_t nComb = 3;
const Char_t *CombName[nComb] = {"V0C_V0A", "ADC_ADA", "V0C_V0A_ADC_ADA"};

using namespace std; // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskFlatenicity) // classimp: necessary for root

    AliAnalysisTaskFlatenicity::AliAnalysisTaskFlatenicity()
    : AliAnalysisTaskSE(), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
      fUseMC(kFALSE), fIsCalib(kFALSE), fIsEqualALICE(kTRUE), fVtxz(-1),
      fParVtx(0x0), ParVtxNorm(-1), fV0Mindex(-1), fmultTPC(-1), fmultV0A(-1),
      fmultV0C(-1), fmultADA(-1), fmultADC(-1), fmultTPCmc(-1), fmultV0Amc(-1),
      fmultV0Cmc(-1), fmultADAmc(-1), fmultADCmc(-1), fDetFlat("V0"),
      fRemoveTrivialScaling(kFALSE), fnGen(-1), fnDetec(-1), fnRecon(-1),
      fPIDResponse(0x0), fTrackFilter(0x0), fOutputList(0), fEtaCut(0.8),
      fPtMin(0.5), ftrackmult08(0), fv0mpercentile(0), fFlat(-1), fFlatMC(-1),
      fMultSelection(0x0), hPtPrimIn(0), hPtPrimOut(0), hPtSecOut(0), hPtOut(0),
      hFlatV0vsFlatTPC(0), hFlatenicityBefore(0), hFlatenicity(0),
      hFlatenicityMC(0), hFlatResponse(0), hFlatVsPt(0), hFlatVsPtMC(0),
      hActivityV0DataSectBefore(0), hActivityV0DataSect(0), hV0vsVtxz(0),
      hActivityV0McSect(0), hFlatVsNchMC(0), hFlatVsV0M(0), hFlatMCVsV0M(0),
      hEtamc(0), hEtamcAlice(0), hCounter(0), hCountProduV0m(0),
      hCountAuthV0m(0), hCountProdu_FlatMC(0), hCountAuth_FlatMC(0),
      hMultMCmVsV0M(0), hMultMCaVsV0M(0), hMultMCcVsV0M(0), hMultmVsV0M(0),
      hMultmVsV0Malice(0), hMultaVsV0M(0), hMultcVsV0M(0), hV0MBadruns(0),
      hChgProdu_All_pt(0), hChgAuth_All_pt(0), hChgProdu_pt_V0(0),
      hChgAuth_pt_V0(0), hChgProdu_pt_Flat(0), hChgAuth_pt_Flat(0) {
  for (Int_t i_c = 0; i_c < nCent; ++i_c) {
    hFlatVsPtV0M[i_c] = 0;
  }
  for (Int_t i_c = 0; i_c < nCent; ++i_c) {
    hFlatVsPtV0MMC[i_c] = 0;
  }
  for (Int_t i_d = 0; i_d < nDet; ++i_d) {
    hComponentsMult[i_d] = 0;
  }
  for (Int_t i_c = 0; i_c < nComb; ++i_c) {
    hCombinedMult[i_c] = 0;
  }
  for (Int_t i_d = 0; i_d < nDet; ++i_d) {
    hComponentsMultmc[i_d] = 0;
  }
  for (Int_t i_c = 0; i_c < nComb; ++i_c) {
    hCombinedMultmc[i_c] = 0;
  }
  for (Int_t i_c = 0; i_c < nComb; ++i_c) {
    hRmCombinedMult[i_c] = 0;
  }
  for (Int_t i_c = 0; i_c < nCent; ++i_c) {
    hMultMCmVsFlat[i_c] = 0;
  }
  for (Int_t i_c = 0; i_c < nCent; ++i_c) {
    hMultmVsFlat[i_c] = 0;
  }
}

//_____________________________________________________________________________
AliAnalysisTaskFlatenicity::AliAnalysisTaskFlatenicity(const char *name)
    : AliAnalysisTaskSE(name), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
      fUseMC(kFALSE), fIsCalib(kFALSE), fIsEqualALICE(kTRUE), fVtxz(-1),
      fParVtx(0x0), ParVtxNorm(-1), fV0Mindex(-1), fmultTPC(-1), fmultV0A(-1),
      fmultV0C(-1), fmultADA(-1), fmultADC(-1), fmultTPCmc(-1), fmultV0Amc(-1),
      fmultV0Cmc(-1), fmultADAmc(-1), fmultADCmc(-1), fDetFlat("V0"),
      fRemoveTrivialScaling(kFALSE), fnGen(-1), fnDetec(-1), fnRecon(-1),
      fPIDResponse(0x0), fTrackFilter(0x0), fOutputList(0), fEtaCut(0.8),
      fPtMin(0.5), ftrackmult08(0), fv0mpercentile(0), fFlat(-1), fFlatMC(-1),
      fMultSelection(0x0), hPtPrimIn(0), hPtPrimOut(0), hPtSecOut(0), hPtOut(0),
      hFlatV0vsFlatTPC(0), hFlatenicityBefore(0), hFlatenicity(0),
      hFlatenicityMC(0), hFlatResponse(0), hFlatVsPt(0), hFlatVsPtMC(0),
      hActivityV0DataSectBefore(0), hActivityV0DataSect(0), hV0vsVtxz(0),
      hActivityV0McSect(0), hFlatVsNchMC(0), hFlatVsV0M(0), hFlatMCVsV0M(0),
      hEtamc(0), hEtamcAlice(0), hCounter(0), hCountProduV0m(0),
      hCountAuthV0m(0), hCountProdu_FlatMC(0), hCountAuth_FlatMC(0),
      hMultMCmVsV0M(0), hMultMCaVsV0M(0), hMultMCcVsV0M(0), hMultmVsV0M(0),
      hMultmVsV0Malice(0), hMultaVsV0M(0), hMultcVsV0M(0), hV0MBadruns(0),
      hChgProdu_All_pt(0), hChgAuth_All_pt(0), hChgProdu_pt_V0(0),
      hChgAuth_pt_V0(0), hChgProdu_pt_Flat(0), hChgAuth_pt_Flat(0)

{
  for (Int_t i_c = 0; i_c < nCent; ++i_c) {
    hFlatVsPtV0M[i_c] = 0;
  }
  for (Int_t i_c = 0; i_c < nCent; ++i_c) {
    hFlatVsPtV0MMC[i_c] = 0;
  }
  for (Int_t i_d = 0; i_d < nDet; ++i_d) {
    hComponentsMult[i_d] = 0;
  }
  for (Int_t i_c = 0; i_c < nComb; ++i_c) {
    hCombinedMult[i_c] = 0;
  }
  for (Int_t i_d = 0; i_d < nDet; ++i_d) {
    hComponentsMultmc[i_d] = 0;
  }
  for (Int_t i_c = 0; i_c < nComb; ++i_c) {
    hCombinedMultmc[i_c] = 0;
  }
  for (Int_t i_c = 0; i_c < nComb; ++i_c) {
    hRmCombinedMult[i_c] = 0;
  }
  for (Int_t i_c = 0; i_c < nCent; ++i_c) {
    hMultMCmVsFlat[i_c] = 0;
  }
  for (Int_t i_c = 0; i_c < nCent; ++i_c) {
    hMultmVsFlat[i_c] = 0;
  }

  DefineInput(0, TChain::Class()); // define the input of the analysis: in this
                                   // case you take a 'chain' of events
  DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this
                                   // case it's a list of histograms
}

//_____________________________________________________________________________
AliAnalysisTaskFlatenicity::~AliAnalysisTaskFlatenicity() {
  // destructor
  if (fOutputList) {
    delete fOutputList; // at the end of your task, it is deleted from memory by
                        // calling this function
    fOutputList = 0x0;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskFlatenicity::UserCreateOutputObjects() {

  // create track filters
  fTrackFilter = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts *fCuts = new AliESDtrackCuts();
  fCuts->SetAcceptKinkDaughters(kFALSE);
  fCuts->SetRequireTPCRefit(kTRUE);
  fCuts->SetRequireITSRefit(kTRUE);
  fCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  fCuts->SetDCAToVertex2D(kFALSE);
  fCuts->SetRequireSigmaToVertex(kFALSE);
  fCuts->SetEtaRange(-0.8, 0.8);
  fCuts->SetMinNCrossedRowsTPC(70);
  fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fCuts->SetMaxChi2PerClusterTPC(4);
  fCuts->SetMaxDCAToVertexZ(2);
  //  fCuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);
  //  fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  //  fCuts->SetMaxChi2PerClusterTPC(4);
  //  fCuts->SetMaxDCAToVertexZ(2);
  fCuts->SetMaxChi2PerClusterITS(36);
  fCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
  //  fCuts->SetMaxChi2PerClusterITS(36);
  fTrackFilter->AddCuts(fCuts);

  float min_flat = -0.01;
  float max_flat = 1.01;
  int nbins_flat = 1020;
  if (fRemoveTrivialScaling) {
    min_flat = -0.1;
    max_flat = 9.9;
    nbins_flat = 2000;
  }

  // create output objects
  fParVtx = new TF1("vtxpar", "pol2", -15, 15);
  fParVtx->SetParameters(89.8737, 0.127185, 0.00572492);
  ParVtxNorm = 89.943;

  OpenFile(1);
  fOutputList =
      new TList(); // this is a list which will contain all of your histograms
  fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all

  hFlatV0vsFlatTPC =
      new TH2D("hFlatV0vsFlatTPC", "counter", nbins_flat, min_flat, max_flat,
               nbins_flat, min_flat, max_flat);
  fOutputList->Add(hFlatV0vsFlatTPC);

  hFlatenicityBefore =
      new TH1D("hFlatenicityBefore", "counter", nbins_flat, min_flat, max_flat);
  fOutputList->Add(hFlatenicityBefore);

  hFlatenicity =
      new TH1D("hFlatenicity", "counter", nbins_flat, min_flat, max_flat);
  fOutputList->Add(hFlatenicity);

  hFlatVsPt =
      new TH2D("hFlatVsPt", "Measured; Flatenicity; #it{p}_{T} (GeV/#it{c})",
               nbins_flat, min_flat, max_flat, nPtbins, Ptbins);
  fOutputList->Add(hFlatVsPt);

  for (Int_t i_c = 0; i_c < nCent; ++i_c) {
    hFlatVsPtV0M[i_c] = new TH2D(
        Form("hFlatVsPtV0M_c%d", i_c),
        Form("Measured %1.0f-%1.0f%%V0M; Flatenicity; #it{p}_{T} (GeV/#it{c})",
             centClass[i_c], centClass[i_c + 1]),
        nbins_flat, min_flat, max_flat, nPtbins, Ptbins);
    fOutputList->Add(hFlatVsPtV0M[i_c]);
  }

  for (Int_t i_d = 0; i_d < nDet; ++i_d) {
    hComponentsMult[i_d] = new TH2D(Form("hAmpl_%s", DetName[i_d]), "", 5000,
                                    -0.5, 20000.0, 200, -0.5, 199.5);
    fOutputList->Add(hComponentsMult[i_d]);
  }
  for (Int_t i_c = 0; i_c < nComb; ++i_c) {
    hCombinedMult[i_c] = new TH2D(Form("hCombined_%s", CombName[i_c]), "", 500,
                                  -0.5, 499.5, 200, -0.5, 199.5);
    fOutputList->Add(hCombinedMult[i_c]);
  }

  if (fUseMC) {

    const int nEtaBinsAlice = 9;
    float EtaBinsAlice[nEtaBinsAlice + 1] = {-7.0, -4.9, -3.7, -1.7, 0.0,
                                             2.8,  4.5,  4.8,  6.3,  7.0};

    hEtamc = new TH1D("hEtamc", "", 140, -7.0, 7.0);
    fOutputList->Add(hEtamc);

    hEtamcAlice = new TH1D("hEtamcAlice", "", nEtaBinsAlice, EtaBinsAlice);
    fOutputList->Add(hEtamcAlice);

    hPtPrimIn =
        new TH1D("hPtPrimIn", "Prim In; #it{p}_{T} (GeV/#it{c}; counts)",
                 nPtbins, Ptbins);
    fOutputList->Add(hPtPrimIn);

    hPtPrimOut =
        new TH1D("hPtPrimOut", "Prim Out; #it{p}_{T} (GeV/#it{c}; counts)",
                 nPtbins, Ptbins);
    fOutputList->Add(hPtPrimOut);

    hPtSecOut =
        new TH1D("hPtSecOut", "Sec Out; #it{p}_{T} (GeV/#it{c}; counts)",
                 nPtbins, Ptbins);
    fOutputList->Add(hPtSecOut);

    hPtOut = new TH1D("hPtOut", "all Out; #it{p}_{T} (GeV/#it{c}; counts)",
                      nPtbins, Ptbins);
    fOutputList->Add(hPtOut);

    hFlatenicityMC =
        new TH1D("hFlatenicityMC", "counter", nbins_flat, min_flat, max_flat);
    fOutputList->Add(hFlatenicityMC);
    hFlatResponse =
        new TH2D("hFlatResponse", "; true flat; measured flat", nbins_flat,
                 min_flat, max_flat, nbins_flat, min_flat, max_flat);
    fOutputList->Add(hFlatResponse);
    hFlatVsPtMC =
        new TH2D("hFlatVsPtMC", "MC true; Flatenicity; #it{p}_{T} (GeV/#it{c})",
                 nbins_flat, min_flat, max_flat, nPtbins, Ptbins);
    fOutputList->Add(hFlatVsPtMC);

    hFlatVsNchMC = new TH2D("hFlatVsNchMC", "; true flat; true Nch", nbins_flat,
                            min_flat, max_flat, 100, -0.5, 99.5);
    fOutputList->Add(hFlatVsNchMC);

    for (Int_t i_d = 0; i_d < nDet; ++i_d) {
      hComponentsMultmc[i_d] = new TH2D(Form("hTrueMult_%s", DetName[i_d]), "",
                                        600, -0.5, 599.0, 200, -0.5, 199.5);
      fOutputList->Add(hComponentsMultmc[i_d]);
    }
    for (Int_t i_c = 0; i_c < nComb; ++i_c) {
      hCombinedMultmc[i_c] = new TH2D(Form("hTrueCombined_%s", CombName[i_c]),
                                      "", 500, -0.5, 499.5, 200, -0.5, 199.5);
      fOutputList->Add(hCombinedMultmc[i_c]);
    }
    for (Int_t i_c = 0; i_c < nComb; ++i_c) {
      hRmCombinedMult[i_c] =
          new TH2D(Form("hRmCombined_%s", CombName[i_c]),
                   "; measured combined mult.; true combined mult.", 500, -0.5,
                   499.5, 500, -0.5, 499.5);
      fOutputList->Add(hRmCombinedMult[i_c]);
    }
    for (Int_t i_c = 0; i_c < nCent; ++i_c) {
      hFlatVsPtV0MMC[i_c] =
          new TH2D(Form("hFlatVsPtV0MMC_c%d", i_c),
                   Form("Measured %1.0f-%1.0f%%V0M; true Flatenicity; "
                        "#it{p}_{T} (GeV/#it{c})",
                        centClass[i_c], centClass[i_c + 1]),
                   nbins_flat, min_flat, max_flat, nPtbins, Ptbins);
      fOutputList->Add(hFlatVsPtV0MMC[i_c]);
    }

    hFlatMCVsV0M = new TH2D("hFlatMCVsV0M", "", nCent, centClass, nbins_flat,
                            min_flat, max_flat);
    fOutputList->Add(hFlatMCVsV0M);

    for (Int_t i_c = 0; i_c < nCent; ++i_c) {
      hMultMCmVsFlat[i_c] = new TH2D(Form("hMultMCmVsFlat_c%d", i_c), "", 100,
                                     0.0, 1.0, 1000, -0.5, 999.5);
      fOutputList->Add(hMultMCmVsFlat[i_c]);
    }
    hMultMCmVsV0M =
        new TH2D("hMultMCmVsV0M", "", nCent, centClass, 1000, -0.5, 999.5);
    fOutputList->Add(hMultMCmVsV0M);

    hMultMCaVsV0M =
        new TH2D("hMultMCaVsV0M", "", nCent, centClass, 1000, -0.5, 999.5);
    fOutputList->Add(hMultMCaVsV0M);

    hMultMCcVsV0M =
        new TH2D("hMultMCcVsV0M", "", nCent, centClass, 1000, -0.5, 999.5);
    fOutputList->Add(hMultMCcVsV0M);

    //  hCountEvent = new TH1D("hCountEvent", "event counter", 2, 0.5, 2.5);
    hCountProduV0m = new TH1D("hCountProduV0m", "CountChgV0m_prod;Events;V0m",
                              nCent, centClass);
    fOutputList->Add(hCountProduV0m);

    hCountAuthV0m = new TH1D("hCountAuthV0m", "CountChgV0m_auth;Events;V0m",
                             nCent, centClass);
    fOutputList->Add(hCountAuthV0m);

    hCountProdu_FlatMC =
        new TH1D("hCountProdu_FlatMC", "CountChgFlat_prod;Events;Flatenicity",
                 nbins_flat, min_flat, max_flat);
    fOutputList->Add(hCountProdu_FlatMC);

    hCountAuth_FlatMC =
        new TH1D("hCountAuth_FlatMC", "CountChgFlat_auth;Events;Flatenicity",
                 nbins_flat, min_flat, max_flat);
    fOutputList->Add(hCountAuth_FlatMC);

    hChgProdu_All_pt =
        new TH1D("hChgProdu_All_pt",
                 "hChgProdu_All_pt;;#it{p}_{T} (GeV/#it{c})", nPtbins, Ptbins);
    fOutputList->Add(hChgProdu_All_pt);

    hChgAuth_All_pt =
        new TH1D("hChgAuth_All_pt", "hChgAuth_All_pt;;#it{p}_{T} (GeV/#it{c})",
                 nPtbins, Ptbins);
    fOutputList->Add(hChgAuth_All_pt);

    hChgProdu_pt_V0 =
        new TH2D("hChgProdu_pt_V0", "; V0M Percentile; #it{p}_{T} (GeV/#it{c})",
                 nCent, centClass, nPtbins, Ptbins);
    fOutputList->Add(hChgProdu_pt_V0);

    hChgAuth_pt_V0 =
        new TH2D("hChgAuth_pt_V0", "; V0M Percentile; #it{p}_{T} (GeV/#it{c})",
                 nCent, centClass, nPtbins, Ptbins);
    fOutputList->Add(hChgAuth_pt_V0);

    hChgAuth_pt_Flat =
        new TH2D("hChgAuth_pt_Flat", "; Flatenicity; #it{p}_{T} (GeV/#it{c})",
                 nbins_flat, min_flat, max_flat, nPtbins, Ptbins);
    fOutputList->Add(hChgAuth_pt_Flat);

    hChgProdu_pt_Flat =
        new TH2D("hChgProdu_pt_Flat", "; Flatenicity; #it{p}_{T} (GeV/#it{c})",
                 nbins_flat, min_flat, max_flat, nPtbins, Ptbins);
    fOutputList->Add(hChgProdu_pt_Flat);
  }

  hActivityV0DataSectBefore = new TProfile(
      "hActivityV0DataSectBefore",
      "rec; V0 sector; (before calib) #LTmultiplicity#GT", 64, -0.5, 63.5);
  fOutputList->Add(hActivityV0DataSectBefore);

  hActivityV0DataSect = new TProfile(
      "hActivityV0DataSect", "rec; V0 sector; (after calib) #LTmultiplicity#GT",
      64, -0.5, 63.5);
  fOutputList->Add(hActivityV0DataSect);

  hV0vsVtxz = new TProfile("hV0vsVtxz", ";total amplitude; vtx_z", 30, -15, 15);
  fOutputList->Add(hV0vsVtxz);

  if (fUseMC) {
    hActivityV0McSect =
        new TProfile("hActivityV0McSect", "true; V0 sector; #LTmultiplicity#GT",
                     64, -0.5, 63.5);
    fOutputList->Add(hActivityV0McSect);
  }

  hFlatVsV0M = new TH2D("hFlatVsV0M", "", nCent, centClass, nbins_flat,
                        min_flat, max_flat);
  fOutputList->Add(hFlatVsV0M);

  hCounter = new TH1D("hCounter", "counter", 15, -0.5, 14.5);
  fOutputList->Add(hCounter);

  for (Int_t i_c = 0; i_c < nCent; ++i_c) {
    hMultmVsFlat[i_c] = new TH2D(Form("hMultmVsFlat_c%d", i_c), "", 100, 0.0,
                                 1.0, 1000, -0.5, 999.5);
    fOutputList->Add(hMultmVsFlat[i_c]);
  }

  hMultmVsV0M =
      new TH2D("hMultmVsV0M", "", nCent, centClass, 1000, -0.5, 999.5);
  fOutputList->Add(hMultmVsV0M);

  hMultmVsV0Malice =
      new TH2D("hMultmVsV0Malice", "", nCent, centClass, 1000, -0.5, 999.5);
  fOutputList->Add(hMultmVsV0Malice);

  hMultaVsV0M =
      new TH2D("hMultaVsV0M", "", nCent, centClass, 1000, -0.5, 999.5);
  fOutputList->Add(hMultaVsV0M);

  hMultcVsV0M =
      new TH2D("hMultcVsV0M", "", nCent, centClass, 1000, -0.5, 999.5);
  fOutputList->Add(hMultcVsV0M);

  hV0MBadruns = new TH1D("hV0MBadruns", "", 1000, -0.5, 999.5);
  fOutputList->Add(hV0MBadruns);

  fEventCuts.AddQAplotsToList(fOutputList);
  PostData(1, fOutputList); // postdata will notify the analysis manager of
                            // changes / updates to the
}

//_____________________________________________________________________________
void AliAnalysisTaskFlatenicity::UserExec(Option_t *) {

  AliVEvent *event = InputEvent();
  if (!event) {
    Error("UserExec", "Could not retrieve event");
    return;
  }

  fESD = dynamic_cast<AliESDEvent *>(event);

  if (!fESD) {
    Printf("%s:%d ESDEvent not found in Input Manager", (char *)__FILE__,
           __LINE__);
    this->Dump();
    return;
  }

  hCounter->Fill(0.0);

  if (fUseMC) {
    //      E S D
    fMC = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fMC) {
      Printf("%s:%d MCEvent not found in Input Manager", (char *)__FILE__,
             __LINE__);
      this->Dump();
      return;
    }
    fMCStack = fMC->Stack();
  }

  AliHeader *headerMC;
  Bool_t isGoodVtxPosMC = kFALSE;

  if (fUseMC) {
    headerMC = fMC->Header();
    AliGenEventHeader *genHeader = headerMC->GenEventHeader();
    TArrayF vtxMC(3); // primary vertex  MC
    vtxMC[0] = 9999;
    vtxMC[1] = 9999;
    vtxMC[2] = 9999; // initialize with dummy
    if (genHeader)
      genHeader->PrimaryVertex(vtxMC);

    if (TMath::Abs(vtxMC[2]) <= 10)
      isGoodVtxPosMC = kTRUE;
  }

  fMultSelection = (AliMultSelection *)fESD->FindListObject("MultSelection");
  if (!fMultSelection)
    cout << "------- No AliMultSelection Object Found --------"
         << fMultSelection << endl;

  fv0mpercentile = -999;
  fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
  int v0multalice = fMultSelection->GetEstimator("V0M")->GetValue();

  fFlatAltMC = -1;
  if (fUseMC) {
    fFlatAltMC = GetFlatenicityMC();
    //    if (fFlatAltMC >= 0) {}
  }

  vector<Float_t> ptMC;
  vector<Int_t> idMC;
  if (isGoodVtxPosMC) {
    fnGen = FillMCarray(ptMC, idMC);
    GetMCchargedTrueDists(fnGen, ptMC, idMC);
  }

  vector<Float_t> ptRecon;
  vector<Int_t> idRecon;
  vector<Int_t> isprimRecon;
  fnRecon = FillArray(ptRecon, idRecon, isprimRecon);

  // Trigger selection
  UInt_t fSelectMask = fInputHandler->IsEventSelected();
  Bool_t isINT7selected = fSelectMask & AliVEvent::kINT7;
  if (!isINT7selected)
    return;
  hCounter->Fill(1.0);

  // Good events
  if (!fEventCuts.AcceptEvent(event)) {
    PostData(1, fOutputList);
    return;
  }

  hCounter->Fill(2.0);

  // Good vertex
  Bool_t hasRecVertex = kFALSE;
  hasRecVertex = HasRecVertex();
  if (!hasRecVertex)
    return;

  vector<Float_t> ptDetMC;
  vector<Int_t> idDetMC;
  if (isGoodVtxPosMC) {
    fnDetec = FillMCarray(ptDetMC, idDetMC);
    GetMCchargedDetDists(fnDetec, ptDetMC, idDetMC);
  }

  hCounter->Fill(3.0);

  // good multiplicity
  //  fMultSelection = (AliMultSelection
  //  *)fESD->FindListObject("MultSelection"); if (!fMultSelection)
  //    cout << "------- No AliMultSelection Object Found --------"
  //         << fMultSelection << endl;

  hCounter->Fill(4.0);
  // Multiplicity Estimation
  //  fv0mpercentile = -999;
  //  fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
  //  int v0multalice = fMultSelection->GetEstimator("V0M")->GetValue();
  hCounter->Fill(10.0);

  for (Int_t i_c = 0; i_c < nCent; ++i_c) {
    if (fv0mpercentile >= centClass[i_c] &&
        fv0mpercentile < centClass[i_c + 1]) {
      fV0Mindex = i_c;
    } else {
      continue;
    }
  }

  Double_t flatenicity_v0 = -1;
  Double_t flatenicity_tpc = GetFlatenicityTPC();
  if (fIsEqualALICE) {
    flatenicity_v0 = GetFlatenicityV0EqualALICE();
    ExtractMultiplicitiesEqualALICE();
  } else {
    flatenicity_v0 = GetFlatenicityV0();
    ExtractMultiplicities();
  }

  float com1mc = 0;
  float com2mc = 0;
  float com3mc = 0;
  if (fUseMC) {
    if (isGoodVtxPosMC) {
      ExtractMultiplicitiesMC();

      float activityMC[4] = {0, 0, 0, 0};
      activityMC[0] = fmultADCmc;
      activityMC[1] = fmultV0Cmc;
      activityMC[2] = fmultV0Amc;
      activityMC[3] = fmultADAmc;
      for (int i_a = 0; i_a < 4; ++i_a) {
        hComponentsMultmc[i_a]->Fill(activityMC[i_a], fmultTPCmc);
      }
      com1mc = fmultV0Amc + fmultV0Cmc;
      com2mc = fmultADAmc + fmultADCmc;
      com3mc = com1mc + com2mc;
      hCombinedMultmc[0]->Fill(com1mc, fmultTPCmc);
      hCombinedMultmc[1]->Fill(com2mc, fmultTPCmc);
      hCombinedMultmc[2]->Fill(com3mc, fmultTPCmc);
    }
  }

  // these values were obtained from LHC16l pass 2 (and MC)
  float avData[4] = {1819.91, 55.6384, 27.6564, 449.373};
  float avExpect[4] = {9.18629, 15.1672, 11.934, 7.47469};
  if (fUseMC) {
    avData[0] = 1380.66;
    avData[1] = 60.009;
    avData[2] = 35.1942;
    avData[3] = 202.348;
  }
  float weigths[4] = {1.0, 1.0, 1.0, 1.0};
  for (int i_a = 0; i_a < 4; ++i_a) {
    weigths[i_a] = avExpect[i_a] / avData[i_a];
  }

  float activity[4] = {0, 0, 0, 0};
  activity[0] = fmultADC;
  activity[1] = fmultV0C;
  activity[2] = fmultV0A;
  activity[3] = fmultADA;
  for (int i_a = 0; i_a < 4; ++i_a) {
    hComponentsMult[i_a]->Fill(activity[i_a], fmultTPC);
  }
  float com1 = weigths[1] * activity[1] + activity[2] * weigths[2];
  float com2 = weigths[0] * activity[0] + activity[3] * weigths[3];
  float com3 = com1 + com2;
  hCombinedMult[0]->Fill(com1, fmultTPC);
  hCombinedMult[1]->Fill(com2, fmultTPC);
  hCombinedMult[2]->Fill(com3, fmultTPC);

  if (fUseMC) {
    hRmCombinedMult[0]->Fill(com1, com1mc);
    hRmCombinedMult[1]->Fill(com2, com2mc);
    hRmCombinedMult[2]->Fill(com3, com3mc);
  }

  fFlat = flatenicity_v0; // default V0
  if (fDetFlat == "VO_TPC") {
    fFlat = (flatenicity_v0 + flatenicity_tpc) / 2.0;
  }
  if (fDetFlat == "TPC") {
    fFlat = flatenicity_tpc;
  }
  if (fDetFlat == "V0") {
    fFlat = flatenicity_v0;
  }

  fFlatMC = -1;
  if ((fUseMC) && (fmultV0Cmc) > 0 && (fmultV0Amc > 0)) {
    fFlatMC = GetFlatenicityMC();
    if (fFlatMC >= 0) {
      hFlatenicityMC->Fill(fFlatMC);
      hFlatResponse->Fill(fFlatMC, fFlat);
      hFlatMCVsV0M->Fill(fv0mpercentile, fFlatMC);
      if (fV0Mindex >= 0) {
        hMultMCmVsFlat[fV0Mindex]->Fill(fFlatMC, fmultV0Cmc + fmultV0Amc);
      }
      hMultMCmVsV0M->Fill(fv0mpercentile, fmultV0Cmc + fmultV0Amc);
      hMultMCcVsV0M->Fill(fv0mpercentile, fmultV0Cmc);
      hMultMCaVsV0M->Fill(fv0mpercentile, fmultV0Amc);

      MakeMCanalysis();
    }
  }

  hFlatV0vsFlatTPC->Fill(flatenicity_tpc, flatenicity_v0);

  if ((fFlat >= 0) && (fmultV0C) > 0 && (fmultV0A > 0)) {
    hFlatenicityBefore->Fill(fFlat);
    if (flatenicity_v0 < 0.9 && flatenicity_tpc < 0.9) {
      hFlatenicity->Fill(fFlat);
    }
    if (fV0Mindex >= 0) {
      if ((fV0Mindex == nCent - 1) && (v0multalice > 400)) {
        hV0MBadruns->Fill(v0multalice);
      } else {
        hFlatVsV0M->Fill(fv0mpercentile, fFlat);
        hMultmVsFlat[fV0Mindex]->Fill(fFlat, fmultV0C + fmultV0A);
        hMultmVsV0M->Fill(fv0mpercentile, fmultV0C + fmultV0A);
        hMultmVsV0Malice->Fill(fv0mpercentile, v0multalice);
        hMultcVsV0M->Fill(fv0mpercentile, fmultV0C);
        hMultaVsV0M->Fill(fv0mpercentile, fmultV0A);
        MakeDataanalysis();
      }
    }
  }

  PostData(1, fOutputList); // stream the result of this event to the output
                            // manager which will write it to a file
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicity::Terminate(Option_t *) {}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicity::MakeDataanalysis() {

  // rec
  Int_t nTracks = fESD->GetNumberOfTracks();
  for (Int_t iT = 0; iT < nTracks; ++iT) {

    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
        fESD->GetTrack(iT)); // get a track (type AliesdTrack)
    if (!esdtrack)
      continue;
    if (!fTrackFilter->IsSelected(esdtrack))
      continue;
    if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
      continue;
    if (esdtrack->Pt() < fPtMin)
      continue;
    hFlatVsPt->Fill(fFlat, esdtrack->Pt());
    hFlatVsPtV0M[fV0Mindex]->Fill(fFlat, esdtrack->Pt());
  }
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicity::MakeMCanalysis() {

  for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

    AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
    if (!particle)
      continue;
    if (!fMC->IsPhysicalPrimary(i))
      continue;
    if (TMath::Abs(particle->Eta()) > fEtaCut)
      continue;
    if (particle->Pt() < fPtMin)
      continue;
    if (TMath::Abs(particle->Charge()) < 0.1)
      continue;
    hFlatVsPtMC->Fill(fFlatMC, particle->Pt());
    hFlatVsPtV0MMC[fV0Mindex]->Fill(fFlatMC, particle->Pt());
    hPtPrimIn->Fill(particle->Pt());
  }

  // rec
  Int_t nTracks = fESD->GetNumberOfTracks();
  for (Int_t iT = 0; iT < nTracks; ++iT) {

    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
        fESD->GetTrack(iT)); // get a track (type AliesdTrack)
    if (!esdtrack)
      continue;
    if (!fTrackFilter->IsSelected(esdtrack))
      continue;
    if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
      continue;
    if (esdtrack->Pt() < fPtMin)
      continue;
    hPtOut->Fill(esdtrack->Pt());
    Int_t mcLabel = -1;
    mcLabel = TMath::Abs(esdtrack->GetLabel());
    if (fMC->IsPhysicalPrimary(mcLabel)) {
      hPtPrimOut->Fill(esdtrack->Pt());
    } else {
      hPtSecOut->Fill(esdtrack->Pt());
    }
  }
}

//______________________________________________________________________________
Double_t AliAnalysisTaskFlatenicity::GetFlatenicityTPC() {

  const int nRings2 = 4;
  const int nSectors2 = 8;
  const int nCells2 = nRings2 * nSectors2;
  float maxEta2[nRings2] = {-0.4, 0.0, +0.4, +0.8};
  float minEta2[nRings2] = {-0.8, -0.4, +0.0, +0.4};
  float maxPhi2[nSectors2] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  float minPhi2[nSectors2] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
  float RhoLattice2[nCells2];
  for (int iCh = 0; iCh < nCells2; iCh++) {
    RhoLattice2[iCh] = 0.0;
  }
  int mult_glob = 0;
  Int_t nTracks = fESD->GetNumberOfTracks();
  for (Int_t iT = 0; iT < nTracks; ++iT) {

    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
        fESD->GetTrack(iT)); // get a track (type AliesdTrack)
    if (!esdtrack)
      continue;
    if (!fTrackFilter->IsSelected(esdtrack))
      continue;
    float eta_a = esdtrack->Eta();
    float phi_a = esdtrack->Phi();

    if (TMath::Abs(eta_a) > fEtaCut)
      continue;
    if (esdtrack->Pt() < fPtMin)
      continue;
    int i_ch = 0;
    for (int ir = 0; ir < nRings2; ir++) {
      for (int is = 0; is < nSectors2; is++) {
        if (eta_a >= minEta2[ir] && eta_a < maxEta2[ir] &&
            phi_a >= minPhi2[is] * 2.0 * M_PI / (1.0 * nSectors2) &&
            phi_a < maxPhi2[is] * 2.0 * M_PI / (1.0 * nSectors2)) {
          RhoLattice2[i_ch]++;
          mult_glob++;
        }
        i_ch++;
      }
    }
  }

  double mRho_glob = 0;
  for (int iCell = 0; iCell < nCells2; ++iCell) {
    mRho_glob += 1.0 * RhoLattice2[iCell];
  }
  // average activity per cell
  mRho_glob /= (1.0 * nCells2);
  // get sigma
  double sRho_glob_tmp = 0;
  for (int iCell = 0; iCell < nCells2; ++iCell) {
    sRho_glob_tmp += TMath::Power(1.0 * RhoLattice2[iCell] - mRho_glob, 2);
  }
  sRho_glob_tmp /= (1.0 * nCells2 * nCells2);
  double sRho_glob = TMath::Sqrt(sRho_glob_tmp);
  float flatenicity_glob = 9999;
  if (mRho_glob > 0) {
    if (fRemoveTrivialScaling) {
      flatenicity_glob = TMath::Sqrt(mult_glob) * sRho_glob / mRho_glob;
    } else {
      flatenicity_glob = sRho_glob / mRho_glob;
    }
  }

  return flatenicity_glob;
}

//______________________________________________________________________________
Double_t AliAnalysisTaskFlatenicity::GetFlatenicityV0EqualALICE() {

  AliVVZERO *lVV0 = 0x0;
  AliVEvent *lVevent = 0x0;
  lVevent = dynamic_cast<AliVEvent *>(InputEvent());
  if (!lVevent) {
    AliWarning("ERROR: ESD / AOD event not available \n");
    return -1;
  }
  // Get VZERO Information for multiplicity later
  lVV0 = lVevent->GetVZEROData();
  if (!lVV0) {
    AliError("AliVVZERO not available");
    return 9999;
  }
  // Flatenicity calculation
  const Int_t nRings = 4;
  const Int_t nSectors = 8;
  Float_t minEtaV0C[nRings] = {-3.7, -3.2, -2.7, -2.2};
  Float_t maxEtaV0C[nRings] = {-3.2, -2.7, -2.2, -1.7};
  Float_t maxEtaV0A[nRings] = {5.1, 4.5, 3.9, 3.4};
  Float_t minEtaV0A[nRings] = {4.5, 3.9, 3.4, 2.8};
  // Grid
  const Int_t nCells = nRings * 2 * nSectors;
  float RhoLattice[nCells];
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    RhoLattice[iCh] = 0.0;
  }

  // before calibration
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    float mult = 0;
    // only corrected for vertex
    if (iCh < 32) { // V0C
      mult = AliESDUtils::GetCorrV0C(lVV0->GetMultiplicity(iCh), fVtxz);
    } else { // V0A
      mult = AliESDUtils::GetCorrV0A(lVV0->GetMultiplicity(iCh), fVtxz);
    }

    RhoLattice[iCh] = mult;
    hActivityV0DataSectBefore->Fill(iCh, lVV0->GetMultiplicity(iCh));
  }

  // Filling histos with mult info
  float total_v0_tmp = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    hActivityV0DataSect->Fill(iCh, RhoLattice[iCh]);
    total_v0_tmp += RhoLattice[iCh];
  }
  int total_v0 = total_v0_tmp;

  hV0vsVtxz->Fill(fVtxz, total_v0);

  Int_t nringA = 0;
  Int_t nringC = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    Float_t detaV0 = -1;
    // Float_t mult = lVV0->GetMultiplicity(iCh);
    if (iCh < 32) { // V0C
      if (iCh < 8) {
        nringC = 0;
      } else if (iCh >= 8 && iCh < 16) {
        nringC = 1;
      } else if (iCh >= 16 && iCh < 24) {
        nringC = 2;
      } else {
        nringC = 3;
      }
      detaV0 = maxEtaV0C[nringC] - minEtaV0C[nringC];
    } else { // V0A
      if (iCh < 40) {
        nringA = 0;
      } else if (iCh >= 40 && iCh < 48) {
        nringA = 1;
      } else if (iCh >= 48 && iCh < 56) {
        nringA = 2;
      } else {
        nringA = 3;
      }
      detaV0 = maxEtaV0A[nringA] - minEtaV0A[nringA];
    }
    // consider the different eta coverage
    RhoLattice[iCh] /= detaV0;
  }
  Float_t mRho = 0;
  Float_t flatenicity = -1;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    mRho += RhoLattice[iCh];
  }
  Float_t multiplicityV0M = mRho;
  // average activity per cell
  mRho /= (1.0 * nCells);
  // get sigma
  Double_t sRho_tmp = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
  }
  sRho_tmp /= (1.0 * nCells * nCells);
  Float_t sRho = TMath::Sqrt(sRho_tmp);
  if (mRho > 0) {
    if (fRemoveTrivialScaling) {
      flatenicity = TMath::Sqrt(multiplicityV0M) * sRho / mRho;
    } else {
      flatenicity = sRho / mRho;
    }
  } else {
    flatenicity = 9999;
  }
  return flatenicity;
}

//______________________________________________________________________________
Double_t AliAnalysisTaskFlatenicity::GetFlatenicityV0() {

  AliVVZERO *lVV0 = 0x0;
  AliVEvent *lVevent = 0x0;
  lVevent = dynamic_cast<AliVEvent *>(InputEvent());
  if (!lVevent) {
    AliWarning("ERROR: ESD / AOD event not available \n");
    return -1;
  }
  // Get VZERO Information for multiplicity later
  lVV0 = lVevent->GetVZEROData();
  if (!lVV0) {
    AliError("AliVVZERO not available");
    return 9999;
  }

  // Flatenicity calculation
  const Int_t nRings = 4;
  const Int_t nSectors = 8;
  Float_t minEtaV0C[nRings] = {-3.7, -3.2, -2.7, -2.2};
  Float_t maxEtaV0C[nRings] = {-3.2, -2.7, -2.2, -1.7};
  Float_t maxEtaV0A[nRings] = {5.1, 4.5, 3.9, 3.4};
  Float_t minEtaV0A[nRings] = {4.5, 3.9, 3.4, 2.8};
  // Grid
  const Int_t nCells = nRings * 2 * nSectors;
  float RhoLattice[nCells];
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    RhoLattice[iCh] = 0.0;
  }

  // before calibration
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    Float_t mult = lVV0->GetMultiplicity(iCh);
    RhoLattice[iCh] = mult;
    hActivityV0DataSectBefore->Fill(iCh, RhoLattice[iCh]);
  }
  // after calibration
  if (fIsCalib) {
    for (Int_t iCh = 0; iCh < nCells; iCh++) {
      RhoLattice[iCh] *= fParVtx->Eval(0.0) / fParVtx->Eval(fVtxz);
    }
  }

  // Filling histos with mult info
  float total_v0_tmp = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    hActivityV0DataSect->Fill(iCh, RhoLattice[iCh]);
    total_v0_tmp += RhoLattice[iCh];
  }
  int total_v0 = total_v0_tmp;
  hV0vsVtxz->Fill(fVtxz, total_v0);

  Int_t nringA = 0;
  Int_t nringC = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    Float_t detaV0 = -1;
    // Float_t mult = lVV0->GetMultiplicity(iCh);
    if (iCh < 32) { // V0C
      if (iCh < 8) {
        nringC = 0;
      } else if (iCh >= 8 && iCh < 16) {
        nringC = 1;
      } else if (iCh >= 16 && iCh < 24) {
        nringC = 2;
      } else {
        nringC = 3;
      }
      detaV0 = maxEtaV0C[nringC] - minEtaV0C[nringC];
    } else { // V0A
      if (iCh < 40) {
        nringA = 0;
      } else if (iCh >= 40 && iCh < 48) {
        nringA = 1;
      } else if (iCh >= 48 && iCh < 56) {
        nringA = 2;
      } else {
        nringA = 3;
      }
      detaV0 = maxEtaV0A[nringA] - minEtaV0A[nringA];
    }
    // consider the different eta coverage
    RhoLattice[iCh] /= detaV0;
  }
  Float_t mRho = 0;
  Float_t flatenicity = -1;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    mRho += RhoLattice[iCh];
  }
  Float_t multiplicityV0M = mRho;
  // average activity per cell
  mRho /= (1.0 * nCells);
  // get sigma
  Double_t sRho_tmp = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
  }
  sRho_tmp /= (1.0 * nCells * nCells);
  Float_t sRho = TMath::Sqrt(sRho_tmp);
  if (mRho > 0) {
    if (fRemoveTrivialScaling) {
      flatenicity = TMath::Sqrt(multiplicityV0M) * sRho / mRho;
    } else {
      flatenicity = sRho / mRho;
    }
  } else {
    flatenicity = 9999;
  }
  return flatenicity;
}

//______________________________________________________________________________
Double_t AliAnalysisTaskFlatenicity::GetFlatenicityMC() {

  // Flatenicity calculation
  const Int_t nRings = 8;
  Float_t maxEta[nRings] = {-3.2, -2.7, -2.2, -1.7, 5.1, 4.5, 3.9, 3.4};
  Float_t minEta[nRings] = {-3.7, -3.2, -2.7, -2.2, 4.5, 3.9, 3.4, 2.8};

  const Int_t nSectors = 8;
  Float_t PhiBins[nSectors + 1];
  Float_t deltaPhi = (2.0 * TMath::Pi()) / (1.0 * nSectors);
  for (int i_phi = 0; i_phi < nSectors + 1; ++i_phi) {
    PhiBins[i_phi] = 0;
    if (i_phi < nSectors) {
      PhiBins[i_phi] = i_phi * deltaPhi;
    } else {
      PhiBins[i_phi] = 2.0 * TMath::Pi();
    }
  }

  // Grid
  const Int_t nCells = nRings * nSectors;
  Float_t RhoLattice[nCells];
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    RhoLattice[iCh] = 0.0;
  }

  Int_t nMult = 0;
  for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

    AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
    if (!particle)
      continue;
    if (!fMC->IsPhysicalPrimary(i))
      continue;
    if (particle->Pt() <= 0.0)
      continue;
    if (TMath::Abs(particle->Charge()) < 0.1)
      continue;
    Double_t phi = particle->Phi();
    Double_t eta = particle->Eta();

    Int_t i_segment = 0;
    for (int i_eta = 0; i_eta < nRings; ++i_eta) {

      for (int i_phi = 0; i_phi < nSectors; ++i_phi) {

        if (eta >= minEta[i_eta] && eta < maxEta[i_eta] &&
            phi >= PhiBins[i_phi] && phi < PhiBins[i_phi + 1]) {
          nMult++;
          RhoLattice[i_segment] += 1.0;
        }
        i_segment++;
      }
    }
  }

  Int_t i_segment = 0;
  for (int i_eta = 0; i_eta < nRings; ++i_eta) {
    for (int i_phi = 0; i_phi < nSectors; ++i_phi) {
      Float_t deltaEta = TMath::Abs(maxEta[i_eta] - minEta[i_eta]);
      hActivityV0McSect->Fill(i_segment, RhoLattice[i_segment]);
      RhoLattice[i_segment] /= deltaEta;
      // Filling histos with mult info
      i_segment++;
    }
  }

  Float_t mRho = 0;
  Float_t flatenicity = -1;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    mRho += RhoLattice[iCh];
  }
  // average activity per cell
  mRho /= (1.0 * nCells);
  // get sigma
  Float_t sRho_tmp = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
  }
  sRho_tmp /= (1.0 * nCells * nCells);
  Float_t sRho = TMath::Sqrt(sRho_tmp);

  if (mRho > 0) {
    if (fRemoveTrivialScaling) {
      flatenicity = TMath::Sqrt(1.0 * nMult) * sRho / mRho;
    } else {
      flatenicity = sRho / mRho;
    }
  } else {
    sRho = 9999;
  }
  hFlatVsNchMC->Fill(flatenicity, nMult);
  return flatenicity;
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicity::ExtractMultiplicities() {

  fmultTPC = 0;
  Int_t nTracks = fESD->GetNumberOfTracks();
  for (Int_t iT = 0; iT < nTracks; ++iT) {

    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
        fESD->GetTrack(iT)); // get a track (type AliesdTrack)
    if (!esdtrack)
      continue;
    if (!fTrackFilter->IsSelected(esdtrack))
      continue;
    if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
      continue;
    if (esdtrack->Pt() < fPtMin)
      continue;
    fmultTPC++;
  }

  AliVVZERO *lVV0 = 0x0;
  AliVEvent *lVevent = 0x0;
  lVevent = dynamic_cast<AliVEvent *>(InputEvent());
  if (!lVevent) {
    AliWarning("ERROR: ESD / AOD event not available \n");
    return;
  }
  // Get VZERO Information for multiplicity later
  lVV0 = lVevent->GetVZEROData();
  if (!lVV0) {
    AliError("AliVVZERO not available");
    return;
  }

  const Int_t nChannels = 64;
  float fmultV0C_tmp = 0;
  float fmultV0A_tmp = 0;
  for (Int_t iCh = 0; iCh < nChannels; iCh++) {
    float mult = lVV0->GetMultiplicity(iCh);
    if (fIsCalib) {
      mult *= fParVtx->Eval(0.0) / fParVtx->Eval(fVtxz);
    }
    if (iCh < 32) { // V0C
      fmultV0C_tmp += mult;
    } else { // V0A
      fmultV0A_tmp += mult;
    }
  }

  fmultV0C = fmultV0C_tmp;
  fmultV0A = fmultV0A_tmp;

  AliVAD *lVAD = 0x0;
  lVAD = lVevent->GetADData();
  if (!lVAD) {
    AliError("AliVAD not available");
    return;
  }
  fmultADA = 0;
  fmultADC = 0;
  // Get Multiplicity info per AD 16 channel: C-side : 0-7, A-side 8-15
  for (Int_t i = 0; i < 8; i++) {
    fmultADA += lVAD->GetMultiplicityADA(i);
  }
  for (Int_t i = 8; i < 16; i++) {
    fmultADC += lVAD->GetMultiplicityADC(i - 8);
  }
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicity::ExtractMultiplicitiesEqualALICE() {

  fmultTPC = 0;
  Int_t nTracks = fESD->GetNumberOfTracks();
  for (Int_t iT = 0; iT < nTracks; ++iT) {

    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
        fESD->GetTrack(iT)); // get a track (type AliesdTrack)
    if (!esdtrack)
      continue;
    if (!fTrackFilter->IsSelected(esdtrack))
      continue;
    if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
      continue;
    if (esdtrack->Pt() < fPtMin)
      continue;
    fmultTPC++;
  }

  AliVVZERO *lVV0 = 0x0;
  AliVEvent *lVevent = 0x0;
  lVevent = dynamic_cast<AliVEvent *>(InputEvent());
  if (!lVevent) {
    AliWarning("ERROR: ESD / AOD event not available \n");
    return;
  }
  // Get VZERO Information for multiplicity later
  lVV0 = lVevent->GetVZEROData();
  if (!lVV0) {
    AliError("AliVVZERO not available");
    return;
  }
  float fmultV0A_tmp = 0;
  float fmultV0C_tmp = 0;
  fmultV0A_tmp = AliESDUtils::GetCorrV0A(lVV0->GetMTotV0A(), fVtxz);
  fmultV0C_tmp = AliESDUtils::GetCorrV0C(lVV0->GetMTotV0C(), fVtxz);

  fmultV0A = fmultV0A_tmp;
  fmultV0C = fmultV0C_tmp;

  AliVAD *lVAD = 0x0;
  lVAD = lVevent->GetADData();
  if (!lVAD) {
    AliError("AliVAD not available");
    return;
  }
  fmultADA = 0;
  fmultADC = 0;
  // Get Multiplicity info per AD 16 channel: C-side : 0-7, A-side 8-15
  for (Int_t i = 0; i < 8; i++) {
    fmultADA += lVAD->GetMultiplicityADA(i);
  }
  for (Int_t i = 8; i < 16; i++) {
    fmultADC += lVAD->GetMultiplicityADC(i - 8);
  }
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicity::ExtractMultiplicitiesMC() {

  fmultV0Amc = 0;
  fmultV0Cmc = 0;
  fmultADAmc = 0;
  fmultADCmc = 0;
  fmultTPCmc = 0;

  for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

    AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
    if (!particle)
      continue;
    if (!fMC->IsPhysicalPrimary(i))
      continue;
    if (particle->Pt() <= 0.0)
      continue;
    if (TMath::Abs(particle->Charge()) < 0.1)
      continue;
    Double_t eta_a = particle->Eta();
    hEtamc->Fill(eta_a);
    hEtamcAlice->Fill(eta_a);
    if (eta_a >= 2.8 && eta_a < 5.1) { // v0a acceptance (excluding first ring)
      fmultV0Amc++;
    }
    if (eta_a >= 4.8 && eta_a < 6.3) { // ada acceptance
      fmultADAmc++;
    }
    if (eta_a >= -3.7 && eta_a < -1.7) { // v0c
      fmultV0Cmc++;
    }
    if (eta_a >= -7.0 && eta_a < -4.9) { // adc
      fmultADCmc++;
    }
    if (TMath::Abs(eta_a) < 0.8) { // adc
      fmultTPCmc++;
    }
  }
}

//______________________________________________________________________
Bool_t AliAnalysisTaskFlatenicity::HasRecVertex() {

  float fMaxDeltaSpdTrackAbsolute = 0.5f;
  float fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
  float fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
  float fMaxResolutionSPDvertex = 0.25f;
  float fMaxDispersionSPDvertex = 1.e14f;

  Bool_t fRequireTrackVertex = true;
  unsigned long fFlag;
  fFlag = BIT(AliEventCuts::kNoCuts);

  const AliVVertex *vtTrc = fESD->GetPrimaryVertex();
  bool isTrackV = true;
  if (vtTrc->IsFromVertexer3D() || vtTrc->IsFromVertexerZ())
    isTrackV = false;
  const AliVVertex *vtSPD = fESD->GetPrimaryVertexSPD();

  if (vtSPD->GetNContributors() > 0)
    fFlag |= BIT(AliEventCuts::kVertexSPD);

  if (vtTrc->GetNContributors() > 1 && isTrackV)
    fFlag |= BIT(AliEventCuts::kVertexTracks);

  if (((fFlag & BIT(AliEventCuts::kVertexTracks)) || !fRequireTrackVertex) &&
      (fFlag & BIT(AliEventCuts::kVertexSPD)))
    fFlag |= BIT(AliEventCuts::kVertex);

  const AliVVertex *&vtx =
      bool(fFlag & BIT(AliEventCuts::kVertexTracks)) ? vtTrc : vtSPD;
  AliVVertex *fPrimaryVertex = const_cast<AliVVertex *>(vtx);
  if (!fPrimaryVertex)
    return kFALSE;

  /// Vertex quality cuts
  double covTrc[6], covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  double dz = bool(fFlag & AliEventCuts::kVertexSPD) &&
                      bool(fFlag & AliEventCuts::kVertexTracks)
                  ? vtTrc->GetZ() - vtSPD->GetZ()
                  : 0.; /// If one of the two vertices is not available this cut
                        /// is always passed.
  double errTot = TMath::Sqrt(covTrc[5] + covSPD[5]);
  double errTrc =
      bool(fFlag & AliEventCuts::kVertexTracks) ? TMath::Sqrt(covTrc[5]) : 1.;
  double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
  /// vertex dispersion for run1, only for ESD, AOD code to be added here
  const AliESDVertex *vtSPDESD = dynamic_cast<const AliESDVertex *>(vtSPD);
  double vtSPDdispersion = vtSPDESD ? vtSPDESD->GetDispersion() : 0;
  if ((TMath::Abs(dz) <= fMaxDeltaSpdTrackAbsolute &&
       nsigTot <= fMaxDeltaSpdTrackNsigmaSPD &&
       nsigTrc <=
           fMaxDeltaSpdTrackNsigmaTrack) && // discrepancy track-SPD vertex
      (!vtSPD->IsFromVertexerZ() ||
       TMath::Sqrt(covSPD[5]) <= fMaxResolutionSPDvertex) &&
      (!vtSPD->IsFromVertexerZ() ||
       vtSPDdispersion <= fMaxDispersionSPDvertex) /// vertex dispersion cut for
                                                   /// run1, only for ESD
      ) // quality cut on vertexer SPD z
    fFlag |= BIT(AliEventCuts::kVertexQuality);

  Bool_t hasVtx = (TESTBIT(fFlag, AliEventCuts::kVertex)) &&
                  (TESTBIT(fFlag, AliEventCuts::kVertexQuality));
  fVtxz = vtSPD->GetZ();
  return hasVtx;
}

//______________________________________________________________________
Int_t AliAnalysisTaskFlatenicity::FillMCarray(vector<Float_t> &ptArray,
                                              vector<Int_t> &idArray) {
  /*
     id 0: lambda, id 1: pion, 2: kaon, 3: proton, 4: sigma plus,
     5: sigma minus, 6: Omega, 7: Xi, 8: other charged
   */

  ptArray.clear();
  idArray.clear();

  Int_t nParticlesChg = 0;

  for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {
    AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
    if (!particle)
      continue;
    if (!fMC->IsPhysicalPrimary(i))
      continue;
    if (TMath::Abs(particle->Eta()) > fEtaCut)
      continue;
    if (particle->Pt() < fPtMin)
      continue;
    if (TMath::Abs(particle->Charge()) < 0.1)
      continue;
    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC))
      continue;

    Int_t idPart = -1;
    Int_t partPDG = TMath::Abs(particle->PdgCode());
    if (partPDG == 3122)
      idPart = 0; // lambda
    if (particle->Charge() != 0) {
      if (partPDG == 211)
        idPart = 1; // pions
      else if (partPDG == 321)
        idPart = 2; // kaons
      else if (partPDG == 2212)
        idPart = 3; // protons
      else if (partPDG == 3222)
        idPart = 4; // sigma plus
      else if (partPDG == 3112)
        idPart = 5; // sigma minus
      else if (partPDG == 3334)
        idPart = 6; // Omega
      else if (partPDG == 3312)
        idPart = 7; // Xi
      else
        idPart = 8; // rest of the charged particles
    }

    ptArray.push_back(particle->Pt());
    idArray.push_back(idPart);

    nParticlesChg++;
  }

  return nParticlesChg;
}

//______________________________________________________________________
Int_t AliAnalysisTaskFlatenicity::FillArray(vector<Float_t> &ptArray,
                                            vector<Int_t> &idArray,
                                            vector<Int_t> &isprimArray) {
  /*
     id 0: lambda, id 1: pion, 2: kaon, 3: proton, 4: sigma plus,
     5: sigma minus, 6: Omega, 7: Xi, 8: other charged
   */

  ptArray.clear();
  idArray.clear();
  isprimArray.clear();

  Int_t nParticlesChgRec = 0;

  Int_t nTracks = fESD->GetNumberOfTracks();

  for (Int_t iT = 0; iT < nTracks; ++iT) {
    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
        fESD->GetTrack(iT)); // get a track (type AliesdTrack)
    if (!esdtrack)
      continue;
    if (!fTrackFilter->IsSelected(esdtrack))
      continue;
    if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
      continue;
    if (esdtrack->Pt() < 0.)
      continue;

    //      AliESDtrack *newTrack = 0x0;
    //      newTrack = new AliESDtrack(*esdtrack);
    //      ptArray.push_back(newTrack->Pt());
    ptArray.push_back(esdtrack->Pt());

    Int_t isPrim = -1;
    Int_t mcLabel = -1;
    Int_t idTrack = -1;

    if (fUseMC) { // get label: 0: prim, 1: weak decays, 2: material
      mcLabel = TMath::Abs(esdtrack->GetLabel());
      TParticle *mcParticle = fMC->GetTrack(mcLabel)->Particle();
      if (!mcParticle) {
        printf("----ERROR: mcParticle not available------------------\n");
        continue;
      }
      if (fMC->IsPhysicalPrimary(mcLabel))
        isPrim = 0;

      Int_t partPDG_rec = TMath::Abs(mcParticle->GetPdgCode());
      if (partPDG_rec == 3122)
        idTrack = 0; // lambdas
      else if (partPDG_rec == 211)
        idTrack = 1; // pions
      else if (partPDG_rec == 321)
        idTrack = 2; // kaons
      else if (partPDG_rec == 2212)
        idTrack = 3; // protons
      else if (partPDG_rec == 3222)
        idTrack = 4; // sigma plus
      else if (partPDG_rec == 3112)
        idTrack = 5; // sigma minus
      else if (partPDG_rec == 3334)
        idTrack = 6; // Omega
      else if (partPDG_rec == 3312)
        idTrack = 7; // Xi
      else
        idTrack = 8; // rest of the charged particles
    }

    isprimArray.push_back(isPrim);
    idArray.push_back(idTrack);

    nParticlesChgRec++;
  }

  return nParticlesChgRec;
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicity::GetMCchargedTrueDists(
    Int_t multGen, const vector<Float_t> &ptGen, const vector<Int_t> &idGen) {

  if (multGen < 1)
    return;
  //! ap   The INEL>0 condition

  //  hCountEvent->Fill(1.0);
  hCountProduV0m->Fill(fv0mpercentile);
  hCountProdu_FlatMC->Fill(fFlatAltMC);

  for (Int_t i = 0; i < multGen; ++i) {
    hChgProdu_All_pt->Fill(ptGen[i]);
    hChgProdu_pt_V0->Fill(fv0mpercentile, ptGen[i]);
    hChgProdu_pt_Flat->Fill(fFlatAltMC, ptGen[i]);
  }
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicity::GetMCchargedDetDists(
    Int_t multRec, const vector<Float_t> &ptRec, const vector<Int_t> &idRec) {

  if (multRec < 1)
    return;
  //! ap   The INEL>0 condition
  hCountAuthV0m->Fill(fv0mpercentile);
  hCountAuth_FlatMC->Fill(fFlatAltMC);

  for (Int_t i = 0; i < multRec; ++i) {
    hChgAuth_All_pt->Fill(ptRec[i]);
    hChgAuth_pt_V0->Fill(fv0mpercentile, ptRec[i]);
    hChgAuth_pt_Flat->Fill(fFlatAltMC, ptRec[i]);
  }
}
