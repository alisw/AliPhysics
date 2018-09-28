/*
 * AliFemtoDreamCorrHists.cxx
 *
 *  Created on: Sep 12, 2017
 *      Author: gu74req
 */
#include <iostream>
#include <vector>
#include "AliFemtoDreamCorrHists.h"
#include "TMath.h"
#include "AliLog.h"
ClassImp(AliFemtoDreamCorrHists)
AliFemtoDreamCorrHists::AliFemtoDreamCorrHists()
    : fQA(nullptr),
      fResults(nullptr),
      fPairs(nullptr),
      fPairQA(nullptr),
      fMinimalBooking(false),
      fMomentumResolution(false),
      fPhiEtaPlots(false),
      fRelKThreshold(0.2),
      fSameEventDist(nullptr),
      fSameEventCommonAncestDist(nullptr),
      fSameEventNonCommonAncestDist(nullptr),
      fSameEventMultDist(nullptr),
      fSameEventCentDist(nullptr),
      fSameEventmTDist(nullptr),
      fSameEventkTDist(nullptr),
      fSameEventkTCentDist(nullptr),
      fPairCounterSE(nullptr),
      fMixedEventDist(nullptr),
      fMixedEventMultDist(nullptr),
      fMixedEventCentDist(nullptr),
      fMixedEventmTDist(nullptr),
      fMixedEventkTDist(nullptr),
      fMixedEventkTCentDist(nullptr),
      fPairCounterME(nullptr),
      fMomResolution(nullptr),
      fMomResolutionDist(nullptr),
      fRadiiEtaPhiSE(nullptr),
      fRadiiEtaPhiME(nullptr),
      fRadiiEtaPhiSEsmallK(nullptr),
      fRadiiEtaPhiMEsmallK(nullptr),
      fdEtadPhiSE(nullptr),
      fdEtadPhiME(nullptr),
      fEffMixingDepth(nullptr),
      fDoMultBinning(false),
      fDoCentBinning(false),
      fDokTBinning(false),
      fDomTBinning(false),
      fDokTCentralityBins(false),
      fDoMCCommonAncest(false),
      fdPhidEtaPlots(false),
      fCentBins(0) {
}

AliFemtoDreamCorrHists::AliFemtoDreamCorrHists(
    const AliFemtoDreamCorrHists& hists)
    : fQA(hists.fQA),
      fResults(hists.fResults),
      fPairs(hists.fPairs),
      fPairQA(hists.fPairQA),
      fMinimalBooking(hists.fMinimalBooking),
      fMomentumResolution(hists.fMomResolution),
      fPhiEtaPlots(hists.fPhiEtaPlots),
      fRelKThreshold(hists.fRelKThreshold),
      fSameEventDist(hists.fSameEventDist),
      fSameEventCommonAncestDist(hists.fSameEventCommonAncestDist),
      fSameEventNonCommonAncestDist(hists.fSameEventNonCommonAncestDist),
      fSameEventMultDist(hists.fSameEventMultDist),
      fSameEventCentDist(hists.fSameEventCentDist),
      fSameEventmTDist(hists.fSameEventmTDist),
      fSameEventkTDist(hists.fSameEventkTDist),
      fSameEventkTCentDist(hists.fSameEventkTCentDist),
      fPairCounterSE(hists.fPairCounterSE),
      fMixedEventDist(hists.fMixedEventDist),
      fMixedEventMultDist(hists.fMixedEventMultDist),
      fMixedEventCentDist(hists.fMixedEventCentDist),
      fMixedEventmTDist(hists.fMixedEventmTDist),
      fMixedEventkTDist(hists.fMixedEventkTDist),
      fMixedEventkTCentDist(hists.fMixedEventkTCentDist),
      fPairCounterME(hists.fPairCounterME),
      fMomResolution(hists.fMomResolution),
      fMomResolutionDist(hists.fMomResolutionDist),
      fRadiiEtaPhiSE(hists.fRadiiEtaPhiSE),
      fRadiiEtaPhiME(hists.fRadiiEtaPhiME),
      fRadiiEtaPhiSEsmallK(hists.fRadiiEtaPhiSEsmallK),
      fRadiiEtaPhiMEsmallK(hists.fRadiiEtaPhiMEsmallK),
      fdEtadPhiSE(hists.fdEtadPhiSE),
      fdEtadPhiME(hists.fdEtadPhiME),
      fEffMixingDepth(hists.fEffMixingDepth),
      fDoMultBinning(hists.fDoMultBinning),
      fDoCentBinning(hists.fDoCentBinning),
      fDokTBinning(hists.fDokTBinning),
      fDomTBinning(hists.fDomTBinning),
      fDokTCentralityBins(hists.fDokTCentralityBins),
      fDoMCCommonAncest(hists.fDoMCCommonAncest),
      fdPhidEtaPlots(hists.fdPhidEtaPlots),
      fCentBins(hists.fCentBins) {
}

AliFemtoDreamCorrHists::AliFemtoDreamCorrHists(AliFemtoDreamCollConfig *conf,
                                               bool MinimalBooking)
    : fQA(nullptr),
      fResults(nullptr),
      fPairs(nullptr),
      fPairQA(nullptr),
      fMinimalBooking(MinimalBooking),
      fMomentumResolution(false),
      fPhiEtaPlots(false),
      fRelKThreshold(0.2),
      fSameEventDist(nullptr),
      fSameEventCommonAncestDist(nullptr),
      fSameEventNonCommonAncestDist(nullptr),
      fSameEventMultDist(nullptr),
      fSameEventCentDist(nullptr),
      fSameEventmTDist(nullptr),
      fSameEventkTDist(nullptr),
      fSameEventkTCentDist(nullptr),
      fPairCounterSE(nullptr),
      fMixedEventDist(nullptr),
      fMixedEventMultDist(nullptr),
      fMixedEventCentDist(nullptr),
      fMixedEventmTDist(nullptr),
      fMixedEventkTDist(nullptr),
      fMixedEventkTCentDist(nullptr),
      fPairCounterME(nullptr),
      fMomResolution(nullptr),
      fMomResolutionDist(nullptr),
      fRadiiEtaPhiSE(nullptr),
      fRadiiEtaPhiME(nullptr),
      fRadiiEtaPhiSEsmallK(nullptr),
      fRadiiEtaPhiMEsmallK(nullptr),
      fdEtadPhiSE(nullptr),
      fdEtadPhiME(nullptr),
      fEffMixingDepth(nullptr),
      fDoMultBinning(false),
      fDoCentBinning(false),
      fDokTBinning(false),
      fDomTBinning(false),
      fDokTCentralityBins(false),
      fDoMCCommonAncest(false),
      fdPhidEtaPlots(false),
      fCentBins(0) {
//  fMinimalBooking=MinimalBooking;
  fMomentumResolution = conf->GetDoMomResolution();
  fDoMultBinning = conf->GetDoMultBinning();
  fDoCentBinning = conf->GetDoCentBinning();
  fDokTBinning = conf->GetDokTBinning();
  fDokTCentralityBins = conf->GetDokTCentralityBinning();
  fDomTBinning = conf->GetDomTBinning();
  fPhiEtaPlots = conf->GetDoPhiEtaBinning();
  fDoMCCommonAncest = conf->GetDoSECommonAncestor();
  fdPhidEtaPlots = conf->GetdPhidEtaPlots();
  if (fDokTCentralityBins && !fDokTBinning) {
    AliWarning(
        "Doing the Centrality binning without the kT Binning wont work!\n");
  }
  if (fDokTCentralityBins) {
    fCentBins = conf->GetCentBins();
    if (fCentBins.size() == 0) {
      AliWarning("Did you set the Centrality bins, since their size is 0?\n");
    }
  }
  int multbins = conf->GetNMultBins();
  int nParticles = conf->GetNParticles();
  const int nHists = conf->GetNParticleCombinations();
  std::vector<int> NBinsHist = conf->GetNBinsHist();
  std::vector<int>::iterator itNBins = NBinsHist.begin();
  std::vector<float> kRelMin = conf->GetMinKRel();
  std::vector<float>::iterator itKMin = kRelMin.begin();
  std::vector<float> kRelMax = conf->GetMaxKRel();
  std::vector<float>::iterator itKMax = kRelMax.begin();
  if (nHists != (int) NBinsHist.size() || nHists != (int) kRelMin.size()
      || nHists != (int) kRelMax.size()) {
    //Todo: Replace by AliFatal!
    std::cout << "Initialzing the bin sizes, something went horribly wrong!"
              << std::endl;
  }
  //The way the histograms are assigned later is going to be for example for
  //4 different particle species X1,X2,X3,X4:
  //    X1  X2  X3  X4
  //X1  1   2   3   4
  //X2      5   6   7
  //X3          8   9
  //X4              10<-----Number of the Histogram=Position in input vector
  //X1 corresponds the first particle array in your vector that you hand over
  //in the AliFemtoDreamPartCollection::SetEvent Method, X2 to the second and
  //so on.
  //in case we only book the most basic things we don't need this
  fResults = new TList();
  fResults->SetName("Results");
  fResults->SetOwner();

  if (!fMinimalBooking) {
    fQA = new TList();
    fQA->SetName("PairQA");
    fQA->SetOwner();

    fPairQA = new TList*[nHists];
    fPairCounterSE = new TH2F*[nHists];
    fPairCounterME = new TH2F*[nHists];
    fEffMixingDepth = new TH1F*[nHists];
    if (fMomentumResolution) {
      fMomResolution = new TH2F*[nHists];
      fMomResolutionDist = new TH2F*[nHists];
    } else {
      fMomResolution = nullptr;
      fMomResolutionDist = nullptr;
    }
    if (fPhiEtaPlots) {
      fRadiiEtaPhiSE = new TH2F***[nHists];
      fRadiiEtaPhiME = new TH2F***[nHists];
      fRadiiEtaPhiSEsmallK = new TH2F***[nHists];
      fRadiiEtaPhiMEsmallK = new TH2F***[nHists];
    } else {
      fRadiiEtaPhiSE = nullptr;
      fRadiiEtaPhiME = nullptr;
      fRadiiEtaPhiSEsmallK = nullptr;
      fRadiiEtaPhiMEsmallK = nullptr;
    }
  } else {
    fQA = nullptr;
    fPairQA = nullptr;
    fPairCounterSE = nullptr;
    fPairCounterME = nullptr;
    fEffMixingDepth = nullptr;
    fMomResolution = nullptr;
    fMomResolutionDist = nullptr;
    fRadiiEtaPhiSE = nullptr;
    fRadiiEtaPhiME = nullptr;
    fRadiiEtaPhiSEsmallK = nullptr;
    fRadiiEtaPhiMEsmallK = nullptr;
  }
  //we always want to do this, regardless of the booking type!
  fPairs = new TList*[nHists];
  fSameEventDist = new TH1F*[nHists];
  fMixedEventDist = new TH1F*[nHists];
  if (fDoMultBinning) {
    fSameEventMultDist = new TH2F*[nHists];
    fMixedEventMultDist = new TH2F*[nHists];
  } else {
    fSameEventMultDist = nullptr;
    fMixedEventMultDist = nullptr;
  }
  if (fDoCentBinning) {
    fSameEventCentDist = new TH2F*[nHists];
    fMixedEventCentDist = new TH2F*[nHists];
  } else {
    fSameEventCentDist = nullptr;
    fMixedEventCentDist = nullptr;
  }
  if (fDokTBinning) {
    fSameEventkTDist = new TH2F*[nHists];
    fMixedEventkTDist = new TH2F*[nHists];
  } else {
    fSameEventkTDist = nullptr;
    fMixedEventkTDist = nullptr;
  }
  if (fDokTCentralityBins) {
    fSameEventkTCentDist = new TH2F**[nHists];
    fMixedEventkTCentDist = new TH2F**[nHists];
  } else {
    fSameEventkTCentDist = nullptr;
    fMixedEventkTCentDist = nullptr;
  }
  if (fDomTBinning) {
    fSameEventmTDist = new TH2F*[nHists];
    fMixedEventmTDist = new TH2F*[nHists];
  } else {
    fSameEventmTDist = nullptr;
    fMixedEventmTDist = nullptr;
  }
  if (fDoMCCommonAncest) {
    fSameEventCommonAncestDist = new TH1F*[nHists];
    fSameEventNonCommonAncestDist = new TH1F*[nHists];
  } else {
    fSameEventCommonAncestDist = nullptr;
    fSameEventNonCommonAncestDist = nullptr;
  }
  if (fdPhidEtaPlots) {
    fdEtadPhiSE = new TH2F*[nHists];
    fdEtadPhiME = new TH2F*[nHists];
  } else {
    fdEtadPhiSE = nullptr;
    fdEtadPhiME = nullptr;
  }
  int Counter = 0;
  for (int iPar1 = 0; iPar1 < nParticles; ++iPar1) {
    for (int iPar2 = iPar1; iPar2 < nParticles; ++iPar2) {
      fPairs[Counter] = new TList();
      TString PairFolderName = Form("Particle%d_Particle%d", iPar1, iPar2);
      fPairs[Counter]->SetName(PairFolderName.Data());
      fPairs[Counter]->SetOwner();
      fResults->Add(fPairs[Counter]);

      TString SameEventName = Form("SEDist_Particle%d_Particle%d", iPar1,
                                   iPar2);
      fSameEventDist[Counter] = new TH1F(SameEventName.Data(),
                                         SameEventName.Data(), *itNBins,
                                         *itKMin, *itKMax);
      fSameEventDist[Counter]->Sumw2();
      fPairs[Counter]->Add(fSameEventDist[Counter]);

      TString MixedEventName = Form("MEDist_Particle%d_Particle%d", iPar1,
                                    iPar2);
      fMixedEventDist[Counter] = new TH1F(MixedEventName.Data(),
                                          MixedEventName.Data(), *itNBins,
                                          *itKMin, *itKMax);
      fMixedEventDist[Counter]->Sumw2();
      fPairs[Counter]->Add(fMixedEventDist[Counter]);

      if (fDoMultBinning) {
        TString SameMultEventName = Form("SEMultDist_Particle%d_Particle%d",
                                         iPar1, iPar2);
        fSameEventMultDist[Counter] = new TH2F(SameMultEventName.Data(),
                                               SameMultEventName.Data(),
                                               *itNBins, *itKMin, *itKMax,
                                               multbins, 1, multbins + 1);
        fSameEventMultDist[Counter]->Sumw2();
        fPairs[Counter]->Add(fSameEventMultDist[Counter]);

        TString MixedMultEventName = Form("MEMultDist_Particle%d_Particle%d",
                                          iPar1, iPar2);
        fMixedEventMultDist[Counter] = new TH2F(MixedMultEventName.Data(),
                                                MixedMultEventName.Data(),
                                                *itNBins, *itKMin, *itKMax,
                                                multbins, 1, multbins + 1);
        fMixedEventMultDist[Counter]->Sumw2();
        fPairs[Counter]->Add(fMixedEventMultDist[Counter]);

      }

      if (fDoCentBinning) {
        //        pp 13 TeV
        //        Mult (%)  Class η dNch/dη sys AN  Paper
        //        0 - 1 INEL>0  |η|<0.5 26.22 ±0.28 link  in preparation
        //        0 - 0.01      36.75 ±0.48
        //        0.01 - 0.05     32.43 ±0.37
        //        0.05 - 0.1      30.50 ±0.34
        //        0.1 - 1     25.59 ±0.27
        //        1 - 5     20.05 ±0.21
        //        5 - 10      16.18 ±0.17
        //        10 - 15     13.78 ±015
        //        15 - 20     12.02 ±0.13
        //        20 - 30     10.02 ±0.11
        //        30 - 40     7.93  ±0.08
        //        40 - 50     6.29  ±0.07
        //        50 - 70     4.45  ±0.05
        //        70 - 100      2.42  ±0.03
        Double_t centBins[14] = { 0, 0.01, 0.05, 0.1, 1, 5, 10, 15, 20, 30, 40,
            50, 70, 100 };
        TString SameCentEventName = Form("SECentDist_Particle%d_Particle%d",
                                         iPar1, iPar2);
        fSameEventCentDist[Counter] = new TH2F(SameCentEventName.Data(),
                                               SameCentEventName.Data(),
                                               *itNBins, *itKMin, *itKMax, 13,
                                               centBins);
        fSameEventCentDist[Counter]->Sumw2();
        fPairs[Counter]->Add(fSameEventCentDist[Counter]);

        TString MixedCentEventName = Form("MECentDist_Particle%d_Particle%d",
                                          iPar1, iPar2);
        fMixedEventCentDist[Counter] = new TH2F(MixedCentEventName.Data(),
                                                MixedCentEventName.Data(),
                                                *itNBins, *itKMin, *itKMax, 13,
                                                centBins);
        fMixedEventCentDist[Counter]->Sumw2();
        fPairs[Counter]->Add(fMixedEventCentDist[Counter]);
      }

      if (fDokTBinning) {
        TString SamekTEventName = Form("SEkTDist_Particle%d_Particle%d", iPar1,
                                       iPar2);
        fSameEventkTDist[Counter] = new TH2F(SamekTEventName.Data(),
                                             SamekTEventName.Data(), *itNBins,
                                             *itKMin, *itKMax, *itNBins / 10,
                                             *itKMin, *itKMax * 1.5);
        fSameEventkTDist[Counter]->Sumw2();
        fPairs[Counter]->Add(fSameEventkTDist[Counter]);

        TString MixedkTEventName = Form("MEkTDist_Particle%d_Particle%d", iPar1,
                                        iPar2);
        fMixedEventkTDist[Counter] = new TH2F(MixedkTEventName.Data(),
                                              MixedkTEventName.Data(), *itNBins,
                                              *itKMin, *itKMax, *itNBins / 10,
                                              *itKMin, *itKMax * 1.5);
        fMixedEventkTDist[Counter]->Sumw2();
        fPairs[Counter]->Add(fMixedEventkTDist[Counter]);

      }
      if (fDokTCentralityBins) {
        const int nCentBins = fCentBins.size();
        fSameEventkTCentDist[Counter] = new TH2F*[nCentBins];
        fMixedEventkTCentDist[Counter] = new TH2F*[nCentBins];
        for (int iCent = 0; iCent < nCentBins; ++iCent) {
          TString SamekTCentEventName = Form(
              "SEkTCentDist_Cent_%.0f_Particle%d_Particle%d", fCentBins[iCent],
              iPar1, iPar2);
          fSameEventkTCentDist[Counter][iCent] = new TH2F(
              SamekTCentEventName.Data(), SamekTCentEventName.Data(), *itNBins,
              *itKMin, *itKMax, *itNBins / 10, *itKMin, *itKMax * 1.5);
          fSameEventkTCentDist[Counter][iCent]->Sumw2();
          fPairs[Counter]->Add(fSameEventkTCentDist[Counter][iCent]);

          TString MixedkTCentEventName = Form(
              "MEkTCentDistCent%.0f_Particle%d_Particle%d", fCentBins[iCent],
              iPar1, iPar2);
          fMixedEventkTCentDist[Counter][iCent] = new TH2F(
              MixedkTCentEventName.Data(), MixedkTCentEventName.Data(), *itNBins,
              *itKMin, *itKMax, *itNBins / 10, *itKMin, *itKMax * 1.5);
          fMixedEventkTCentDist[Counter][iCent]->Sumw2();
          fPairs[Counter]->Add(fMixedEventkTCentDist[Counter][iCent]);
        }
      }

      if (fDomTBinning) {
        TString SamemTEventName = Form("SEmTDist_Particle%d_Particle%d", iPar1,
                                       iPar2);
        fSameEventmTDist[Counter] = new TH2F(SamemTEventName.Data(),
                                             SamemTEventName.Data(), *itNBins,
                                             *itKMin, *itKMax, *itNBins / 10,
                                             *itKMin, *itKMax * 1.5);
        fSameEventmTDist[Counter]->Sumw2();
        fPairs[Counter]->Add(fSameEventmTDist[Counter]);

        TString MixedmTEventName = Form("MEmTDist_Particle%d_Particle%d", iPar1,
                                        iPar2);
        fMixedEventmTDist[Counter] = new TH2F(MixedmTEventName.Data(),
                                              MixedmTEventName.Data(), *itNBins,
                                              *itKMin, *itKMax, *itNBins / 10,
                                              *itKMin, *itKMax * 1.5);
        fMixedEventmTDist[Counter]->Sumw2();
        fPairs[Counter]->Add(fMixedEventmTDist[Counter]);
      }

      if (fDoMCCommonAncest) {
        TString SameEventCommonAncestName = Form(
            "SECommonAncestDist_Particle%d_Particle%d", iPar1, iPar2);
        fSameEventCommonAncestDist[Counter] = new TH1F(
            SameEventCommonAncestName.Data(), SameEventCommonAncestName.Data(),
            *itNBins, *itKMin, *itKMax);
        fSameEventCommonAncestDist[Counter]->Sumw2();
        fPairs[Counter]->Add(fSameEventCommonAncestDist[Counter]);

        TString SameEventNonCommonAncestName = Form(
            "SENonCommonAncestDist_Particle%d_Particle%d", iPar1, iPar2);
        fSameEventNonCommonAncestDist[Counter] = new TH1F(
            SameEventNonCommonAncestName.Data(),
            SameEventNonCommonAncestName.Data(), *itNBins, *itKMin, *itKMax);
        fSameEventNonCommonAncestDist[Counter]->Sumw2();
        fPairs[Counter]->Add(fSameEventNonCommonAncestDist[Counter]);
      }

      if (fdPhidEtaPlots) {
        TString SameEventdPhidEtaName = Form(
            "SEdPhidEtaDist_Particle%d_Particle%d", iPar1, iPar2);
        fdEtadPhiSE[Counter] = new TH2F(SameEventdPhidEtaName.Data(),
                                        SameEventdPhidEtaName.Data(), 80, -2.,
                                        2., 84, -2 * TMath::Pi() / 3,
                                        2 * TMath::Pi());
        fdEtadPhiSE[Counter]->GetXaxis()->SetTitle("#Delta#eta");
        fdEtadPhiSE[Counter]->GetYaxis()->SetTitle("#Delta#phi");
        fdEtadPhiSE[Counter]->Sumw2();
        fPairs[Counter]->Add(fdEtadPhiSE[Counter]);

        TString MixedEventdPhidEtaName = Form(
            "MEdPhidEtaDist_Particle%d_Particle%d", iPar1, iPar2);
        fdEtadPhiME[Counter] = new TH2F(MixedEventdPhidEtaName.Data(),
                                        MixedEventdPhidEtaName.Data(), 80, -2.,
                                        2., 84, -2 * TMath::Pi() / 3,
                                        2 * TMath::Pi());
        fdEtadPhiME[Counter]->GetXaxis()->SetTitle("#Delta#eta");
        fdEtadPhiME[Counter]->GetYaxis()->SetTitle("#Delta#phi");
        fdEtadPhiME[Counter]->Sumw2();
        fPairs[Counter]->Add(fdEtadPhiME[Counter]);

      }

      if (!fMinimalBooking) {
        fPairQA[Counter] = new TList();
        TString PairQAName = Form("QA_Particle%d_Particle%d", iPar1, iPar2);
        fPairQA[Counter]->SetName(PairQAName.Data());
        fPairQA[Counter]->SetOwner();
        fQA->Add(fPairQA[Counter]);
        TString PairCounterSEName = Form("SEPairs_Particle%d_Particle%d", iPar1,
                                         iPar2);
        fPairCounterSE[Counter] = new TH2F(PairCounterSEName.Data(),
                                           PairCounterSEName.Data(), 20, 0, 20,
                                           20, 0, 20);
        fPairCounterSE[Counter]->Sumw2();
        fPairCounterSE[Counter]->GetXaxis()->SetTitle(
            Form("Particle%d", iPar1));
        fPairCounterSE[Counter]->GetYaxis()->SetTitle(
            Form("Particle%d", iPar2));
        fPairQA[Counter]->Add(fPairCounterSE[Counter]);

        TString PairCounterMEName = Form("MEPairs_Particle%d_Particle%d", iPar1,
                                         iPar2);
        fPairCounterME[Counter] = new TH2F(PairCounterMEName.Data(),
                                           PairCounterMEName.Data(), 20, 0, 20,
                                           20, 0, 20);
        fPairCounterME[Counter]->Sumw2();
        fPairCounterME[Counter]->GetXaxis()->SetTitle(
            Form("Particle%d", iPar1));
        fPairCounterME[Counter]->GetYaxis()->SetTitle(
            Form("Particle%d", iPar2));
        fPairQA[Counter]->Add(fPairCounterME[Counter]);

        TString EffMixingDepthName = Form(
            "EffMixingDepth_Particle%d_Particle%d", iPar1, iPar2);
        int MixingDepth = conf->GetMixingDepth();
        ++MixingDepth;  //+1 since counting starts at 0
        fEffMixingDepth[Counter] = new TH1F(EffMixingDepthName.Data(),
                                            EffMixingDepthName.Data(),
                                            MixingDepth, -0.5,
                                            MixingDepth - 0.5);
        fEffMixingDepth[Counter]->Sumw2();
        fEffMixingDepth[Counter]->GetXaxis()->SetTitle("MixingDepth");
        fPairQA[Counter]->Add(fEffMixingDepth[Counter]);

        if (fMomentumResolution) {
          //take twice the number of bins we use for the CF to be sure, the range is
          //hard coded. This assumed that the input is in GeV!
          //          *itNBins,*itKMin,*itKMax,
          double dNBin = (*itKMax - *itKMin) / (double) (*itNBins);
          dNBin /= 2;
          int nBims = (int) (1 / dNBin);
          TString MomResoName = Form("MomentumResolution_Particle%d_Particle%d",
                                     iPar1, iPar2);
          fMomResolution[Counter] = new TH2F(MomResoName.Data(),
                                             MomResoName.Data(), nBims, 0, 1,
                                             nBims, 0, 1);
          fMomResolution[Counter]->Sumw2();
          fMomResolution[Counter]->GetXaxis()->SetTitle("k_{Generated}");
          fMomResolution[Counter]->GetYaxis()->SetTitle("k_{Reco}");
          fPairQA[Counter]->Add(fMomResolution[Counter]);

          TString MomResoDistName = Form(
              "MomentumResolutionDist_Particle%d_Particle%d", iPar1, iPar2);
          fMomResolutionDist[Counter] = new TH2F(MomResoDistName.Data(),
                                                 MomResoDistName.Data(), 500,
                                                 -0.3, 0.3, nBims, 0, 1);
          fMomResolutionDist[Counter]->Sumw2();
          fMomResolutionDist[Counter]->GetXaxis()->SetTitle(
              "k_{Reco}-k_{Generated}");
          fMomResolutionDist[Counter]->GetYaxis()->SetTitle("k_{Generated}");
          fPairQA[Counter]->Add(fMomResolutionDist[Counter]);

        }
        if (fPhiEtaPlots) {
          fRadiiEtaPhiSE[Counter] = new TH2F**[3];
          fRadiiEtaPhiME[Counter] = new TH2F**[3];

          fRadiiEtaPhiSEsmallK[Counter] = new TH2F**[3];
          fRadiiEtaPhiMEsmallK[Counter] = new TH2F**[3];
          for (int iDaug = 0; iDaug < 3; ++iDaug) {
            const int nRad = conf->GetNRadii();
            fRadiiEtaPhiSE[Counter][iDaug] = new TH2F*[nRad];
            fRadiiEtaPhiME[Counter][iDaug] = new TH2F*[nRad];
            fRadiiEtaPhiSEsmallK[Counter][iDaug] = new TH2F*[nRad];
            fRadiiEtaPhiMEsmallK[Counter][iDaug] = new TH2F*[nRad];
            for (int iRad = 0; iRad < nRad; ++iRad) {
              TString RadNameSE = Form("SERad_%i_Particle%d_Particle%d_Daug%d",
                                       iRad, iPar1, iPar2, iDaug);
              TString RadNameME = Form("MERad_%i_Particle%d_Particle%d_Daug%d",
                                       iRad, iPar1, iPar2, iDaug);
              fRadiiEtaPhiSE[Counter][iDaug][iRad] = new TH2F(RadNameSE.Data(),
                                                              RadNameSE.Data(),
                                                              200, 0, 0.4, 200,
                                                              0, 0.4);
              fRadiiEtaPhiSE[Counter][iDaug][iRad]->GetXaxis()->SetTitle(
                  "#Delta#eta");
              fRadiiEtaPhiSE[Counter][iDaug][iRad]->GetYaxis()->SetTitle(
                  "#Delta#phi");
              fPairQA[Counter]->Add(fRadiiEtaPhiSE[Counter][iDaug][iRad]);
              fRadiiEtaPhiME[Counter][iDaug][iRad] = new TH2F(RadNameME.Data(),
                                                              RadNameME.Data(),
                                                              200, 0, 0.4, 200,
                                                              0, 0.4);
              fRadiiEtaPhiME[Counter][iDaug][iRad]->GetXaxis()->SetTitle(
                  "#Delta#eta");
              fRadiiEtaPhiME[Counter][iDaug][iRad]->GetYaxis()->SetTitle(
                  "#Delta#phi");
              fPairQA[Counter]->Add(fRadiiEtaPhiME[Counter][iDaug][iRad]);

              RadNameSE += "_smallK";
              RadNameME += "_smallK";

              fRadiiEtaPhiSEsmallK[Counter][iDaug][iRad] = new TH2F(
                  RadNameSE.Data(), RadNameSE.Data(), 200, 0, 0.4, 200, 0, 0.4);
              fRadiiEtaPhiSEsmallK[Counter][iDaug][iRad]->GetXaxis()->SetTitle(
                  "#Delta#eta");
              fRadiiEtaPhiSEsmallK[Counter][iDaug][iRad]->GetYaxis()->SetTitle(
                  "#Delta#phi");
              fPairQA[Counter]->Add(fRadiiEtaPhiSEsmallK[Counter][iDaug][iRad]);
              fRadiiEtaPhiMEsmallK[Counter][iDaug][iRad] = new TH2F(
                  RadNameME.Data(), RadNameME.Data(), 200, 0, 0.4, 200, 0, 0.4);
              fRadiiEtaPhiMEsmallK[Counter][iDaug][iRad]->GetXaxis()->SetTitle(
                  "#Delta#eta");
              fRadiiEtaPhiMEsmallK[Counter][iDaug][iRad]->GetYaxis()->SetTitle(
                  "#Delta#phi");
              fPairQA[Counter]->Add(fRadiiEtaPhiMEsmallK[Counter][iDaug][iRad]);
            }
          }
        }
      }

      ++Counter;
      ++itNBins;
      ++itKMin;
      ++itKMax;
    }
  }
}
AliFemtoDreamCorrHists &AliFemtoDreamCorrHists::operator=(
    const AliFemtoDreamCorrHists& hists) {
  if (this != &hists) {
    this->fQA = hists.fQA;
    this->fResults = hists.fResults;
    this->fPairs = hists.fPairs;
    this->fPairQA = hists.fPairQA;
    this->fMinimalBooking = hists.fMinimalBooking;
    this->fMomentumResolution = hists.fMomResolution;
    this->fPhiEtaPlots = hists.fPhiEtaPlots;
    this->fRelKThreshold = hists.fRelKThreshold;
    this->fSameEventDist = hists.fSameEventDist;
    this->fSameEventCommonAncestDist = hists.fSameEventCommonAncestDist;
    this->fSameEventNonCommonAncestDist = hists.fSameEventNonCommonAncestDist;
    this->fSameEventMultDist = hists.fSameEventMultDist;
    this->fSameEventCentDist = hists.fSameEventCentDist;
    this->fSameEventmTDist = hists.fSameEventmTDist;
    this->fSameEventkTDist = hists.fSameEventkTDist;
    this->fSameEventkTCentDist = hists.fSameEventkTCentDist;
    this->fPairCounterSE = hists.fPairCounterSE;
    this->fMixedEventDist = hists.fMixedEventDist;
    this->fMixedEventMultDist = hists.fMixedEventMultDist;
    this->fMixedEventCentDist = hists.fMixedEventCentDist;
    this->fMixedEventmTDist = hists.fMixedEventmTDist;
    this->fMixedEventkTDist = hists.fMixedEventkTDist;
    this->fMixedEventkTCentDist = hists.fMixedEventkTCentDist;
    this->fPairCounterME = hists.fPairCounterME;
    this->fMomResolution = hists.fMomResolution;
    this->fMomResolutionDist = hists.fMomResolutionDist;
    this->fRadiiEtaPhiSE = hists.fRadiiEtaPhiSE;
    this->fRadiiEtaPhiME = hists.fRadiiEtaPhiME;
    this->fRadiiEtaPhiSEsmallK = hists.fRadiiEtaPhiSEsmallK;
    this->fRadiiEtaPhiMEsmallK = hists.fRadiiEtaPhiMEsmallK;
    this->fdEtadPhiSE = hists.fdEtadPhiSE;
    this->fdEtadPhiME = hists.fdEtadPhiME;
    this->fEffMixingDepth = hists.fEffMixingDepth;
    this->fDoMultBinning = hists.fDoMultBinning;
    this->fDoCentBinning = hists.fDoCentBinning;
    this->fDokTBinning = hists.fDokTBinning;
    this->fDomTBinning = hists.fDomTBinning;
    this->fDokTCentralityBins = hists.fDokTCentralityBins;
    this->fDoMCCommonAncest = hists.fDoMCCommonAncest;
    this->fdPhidEtaPlots = hists.fdPhidEtaPlots;
    this->fCentBins = hists.fCentBins;
  }
  return *this;
}
AliFemtoDreamCorrHists::~AliFemtoDreamCorrHists() {
  if (fPairs) {
    delete[] fPairs;
    delete fPairs;
  }
  if (fSameEventDist) {
    delete[] fSameEventDist;
    delete fSameEventDist;
  }
  if (fSameEventMultDist) {
    delete[] fSameEventMultDist;
    delete fSameEventMultDist;
  }
  if (fSameEventmTDist) {
    delete[] fSameEventmTDist;
    delete fSameEventmTDist;
  }
  if (fSameEventkTDist) {
    delete[] fSameEventkTDist;
    delete fSameEventDist;
  }
  if (fSameEventkTCentDist) {
    delete[] fSameEventkTCentDist;
    delete fSameEventkTCentDist;
  }
  if (fMixedEventDist) {
    delete[] fMixedEventDist;
    delete fMixedEventDist;
  }
  if (fMixedEventMultDist) {
    delete[] fMixedEventMultDist;
    delete fMixedEventMultDist;
  }
  if (fMixedEventmTDist) {
    delete[] fMixedEventmTDist;
    delete fMixedEventmTDist;
  }
  if (fMixedEventkTDist) {
    delete[] fMixedEventkTDist;
    delete fMixedEventDist;
  }
  if (fMixedEventkTCentDist) {
    delete[] fMixedEventkTCentDist;
    delete fMixedEventkTCentDist;
  }
}

void AliFemtoDreamCorrHists::FillSameEventkTCentDist(int i, float kT,
                                                     float RelK, float cent) {
  int centBin = -1;
  if (0 < cent) {
    for (unsigned int iCent = 0; iCent < fCentBins.size(); ++iCent) {
      if (cent < fCentBins[iCent]) {
        centBin = iCent;
        break;
      }
    }
    if (centBin != -1) {
      fSameEventkTCentDist[i][centBin]->Fill(RelK, kT);
    }
  }
}

void AliFemtoDreamCorrHists::FillMixedEventkTCentDist(int i, float kT,
                                                      float RelK, float cent) {
  int centBin = -1;
  if (0 < cent) {
    for (unsigned int iCent = 0; iCent < fCentBins.size(); ++iCent) {
      if (cent < fCentBins[iCent]) {
        centBin = iCent;
        break;
      }
    }
    if (centBin != -1) {
      fMixedEventkTCentDist[i][centBin]->Fill(RelK, kT);
    }
  }
}
