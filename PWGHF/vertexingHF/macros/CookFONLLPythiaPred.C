// Macro to "cook" the FONLL + Pythia predictions of D-meson and Lc cross sections obtained from ComputeBtoDdecay.C .
//
// The BR and FF are set to the DPG values and can be changed. The output has the same structure and histogram names
// as the one used in D2H analyses.
//
// Author: F. Catalano, fabio.catalano@cern.ch
//

#include <string>
#include <array>
#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

enum { // options for branching ratios
  kBROriginal, // keep the values used in the input files
  kBRPDG, // use hard coded values from PDG (2018)
  kBRPDGmix // use hard coded values from PDG (2018) and the B0/B+- admixture
};

enum { // options for fragmentation fractions
  kFFOriginal, // keep the values used in the input files
  kFFppbar, // use hard coded values from ppbar, PDG (2018)
  kFFee  // use hard coded values from e+e-, PDG (2018)
};

enum { // options for re-weight
  kSimple, // estimate a global factor to scale the prediction (keep the input pT dependence)
  kAccurate // correct each b-hadron contribution indipendently (modify also the pT dependence)
};

void CookFONLLPythiaPred(std::string inFileNameMin = "DfromB_FONLLminPythia8_FFppbar_yDcut.root", // min FONLL predictions
                          std::string inFileNameCent = "DfromB_FONLLcentPythia8_FFppbar_yDcut.root", // central FONLL predictions
                          std::string inFileNameMax = "DfromB_FONLLmaxPythia8_FFppbar_yDcut.root",  // max FONLL predictions
                          std::string outFileName = "DmesonLcPredictions_502TeV_y05 _pythia8.root",
                          int brOpt = kBRPDG,
                          int ffOpt = kFFOriginal,
                          int wOpt = kSimple) {

  std::array<std::string, 3> inFileNames = {inFileNameMin, inFileNameCent, inFileNameMax};
  std::array<std::string, 3> edgeNames = {"min", "central", "max"};

  const int numMothers = 4;
  const int numDaughters = 6;
  std::array<std::string, numMothers> mothToDauHistos = {"hB0dau", "hBplusdau", "hBsdau", "hLbdau"};
  std::array<std::string, numDaughters> predFDHistos = {"hnonpromptD0pt", "hnonpromptDpluspt", "hnonpromptDspt",
                                                        "hnonpromptLcpt", "hnonpromptDstarpt", "hnonpromptLcpt"};
  std::array<std::string, numDaughters> predPromptHistos = {"hfonllPromptD0", "hfonllPromptDplus", "hfonllPromptDs",
                                                            "hfonllPromptLc", "hfonllPromptDstar", "hfonllPromptLc"};
  std::array<std::string, numDaughters> predTag = {"D0Kpi", "Dpluskpipi", "DsPhipitoKkpi", "Lcpkpi", "DstarD0pi", "LcK0sp"};
  std::array<std::string, numDaughters - 1> partTag = {"D0", "Dplus", "Ds", "Lc", "Dstar"};
  std::array<double, numMothers> origBFF = {0.34, 0.34, 0.101, 0.219}; // (B0, B+, Bs, Lb) FF used in the input predictions
  std::array<std::array<double, numMothers>, numDaughters> pdgBRfromB = {{{0.555, 0.876, 0.008, 0.},   // D0 and (BRfromB0, BRfromB+, BRfromBs, BRfromLb) from PDG (2018)
                                                                          {0.392, 0.124, 0., 0.},      // D+
                                                                          {0.117, 0.09, 0.93, 0.011},  // Ds
                                                                          {0.066, 0.049, 0., 0.333},   // Lc
                                                                          {0.23, 0.061, 0.003, 0.},    // D*+
                                                                          {0.066, 0.049, 0., 0.333}    // Lc
                                                                          }};
  // B0 and B+ BR from B0/B+- admixture
  std::array<std::array<double, numMothers>, numDaughters> pdgBRfromBmix = {{{0.624, 0.624, 0.008, 0.},   // D0 and (BRfromB0, BRfromB+, BRfromBs, BRfromLb) from PDG (2018) and B0/B+- admixture
                                                                             {0.241, 0.241, 0., 0.},      // D+
                                                                             {0.083, 0.083, 0.93, 0.011}, // Ds
                                                                             {0.036, 0.036, 0., 0.333},   // Lc
                                                                             {0.225, 0.225, 0.003, 0.},   // D*+
                                                                             {0.036, 0.036, 0., 0.333}    // Lc
                                                                             }};
  std::array<double, numMothers> ppbarBFF = {0.34, 0.34, 0.101, 0.219}; // (B0, B+, Bs, Lb) from PDG (2018)
  std::array<double, numMothers> eeBFF = {0.412, 0.412, 0.088, 0.089}; // (B0, B+, Bs, Lb) from PDG (2018)
  std::array<double, numDaughters> decayBR = {0.0389, 0.0898, 0.0227, 0.0623, 0.0263, 0.0158}; // (D0, D+, Ds, Lc->pKpi, D*+, LC->K0sp) from PDG (2018)
  std::array<std::array<double, numMothers>, numDaughters> origBR = {};

  // print fraction of b to b-hadrons FF times PDG BR 
  if(brOpt == kBRPDG || brOpt == kBRPDGmix) {
    std::cout<<"\nb to X factors with PDG BR and FF\n";
    for(int iDau = 0; iDau < numDaughters - 1; iDau++) {
        double frac = 0.;
        for(int iMother = 0; iMother < numMothers; iMother++) {
          double motherFF = 1.;
          double BRfromB = 1.;

          if(ffOpt == kFFOriginal)
            motherFF = origBFF[iMother];
          else if(ffOpt == kFFppbar)
            motherFF = ppbarBFF[iMother];
          else if(ffOpt == kFFee)
            motherFF = eeBFF[iMother];

          if(brOpt == kBRPDG)
            BRfromB = pdgBRfromB[iDau][iMother];
          else if(brOpt == kBRPDGmix)
            BRfromB = pdgBRfromBmix[iDau][iMother];

          frac += motherFF * BRfromB;
        }
        std::cout<<"Factor b to " + partTag[iDau] << ": "<<frac<<std::endl;
    }
    std::cout<<std::endl;
  }

  TFile outFile(outFileName.data(),"recreate");

  for(int iFile = 0; iFile < 3; iFile++) {
    TFile *inFile = TFile::Open(inFileNames[iFile].data());

    // get the original BR
    for(int iMother = 0; iMother < numMothers; iMother++) {
      TH1D *hMothToDau = (TH1D *)inFile->Get(mothToDauHistos[iMother].data());
      double totMothers = hMothToDau->GetBinContent(1);

      for(int iDau = 0; iDau < numDaughters - 1; iDau++) {
        double totDau = hMothToDau->GetBinContent(iDau + 2);
        origBR[iDau][iMother] = 1. * totDau / totMothers; 
      }
      // manually set origBR for the LcK0sp case, it is equal to the LcpKpi case
      origBR[5][iMother] = origBR[3][iMother];
    }

    // print original BR factors from file, only for central case
    if(brOpt == kBROriginal && edgeNames[iFile] == "central") {
      for(int iMother = 0; iMother < numMothers; iMother++) {
        std::cout<<"BR for "<<mothToDauHistos[iMother]<<" extracted from "<<inFileNames[iFile]<<" file\n";
        for(int iDau = 0; iDau < numDaughters - 1; iDau++) 
          std::cout<<"  to  " + partTag[iDau] << ": "<<origBR[iDau][iMother]<<std::endl;
      }
    }

    // print fraction of b to b-hadrons FF times BR extracted from file, only for central case
    if(brOpt == kBROriginal && edgeNames[iFile] == "central") {
      std::cout<<"b to X factors with BR extracted from "<<inFileNames[iFile]<<" file\n";
      for(int iDau = 0; iDau < numDaughters - 1; iDau++) {
          double frac = 0.;
          for(int iMother = 0; iMother < numMothers; iMother++) {
            double motherFF = 1.;
            if(ffOpt == kFFOriginal)
              motherFF = origBFF[iMother];
            else if(ffOpt == kFFppbar)
              motherFF = ppbarBFF[iMother];
            else if(ffOpt == kFFee)
              motherFF = eeBFF[iMother];
            frac += motherFF * origBR[iDau][iMother];
          }
          std::cout<<"Factor b to " + partTag[iDau] << ": "<<frac<<std::endl;
      }
      std::cout<<std::endl;
    }

    // get and correct the predictions
    for(int iDau = 0; iDau < numDaughters; iDau++) {
      TH1D *hDauFDPred = nullptr;

      // crude method, estimate a global correction factor
      if(wOpt == kSimple) {
        // calculate the final correction factor for non-prompt
        double newFrac = 0.;
        double oldFrac = 0.;
        for(int iMother = 0; iMother < numMothers; iMother++) {
          double motherFF = 1.;
          double BRfromB = 1.;

          if(ffOpt == kFFOriginal)
            motherFF = origBFF[iMother];
          else if(ffOpt == kFFppbar)
            motherFF = ppbarBFF[iMother];
          else if(ffOpt == kFFee)
            motherFF = eeBFF[iMother];

          if(brOpt == kBROriginal)
            BRfromB = origBR[iDau][iMother];
          else if(brOpt == kBRPDG)
            BRfromB = pdgBRfromB[iDau][iMother];
          else if(brOpt == kBRPDGmix)
            BRfromB = pdgBRfromBmix[iDau][iMother];

          newFrac += motherFF * BRfromB;
          oldFrac += origBFF[iMother] * origBR[iDau][iMother];
        }

        double corr = newFrac / oldFrac;
        hDauFDPred = (TH1D *)inFile->Get(predFDHistos[iDau].data());
        hDauFDPred->SetDirectory(0);
        hDauFDPred->Scale(corr);

        if(edgeNames[iFile] == "central" && iDau < (numDaughters - 1))
          std::cout<<"Prediction correction factor for " + partTag[iDau] << ", "<<edgeNames[iFile]<<" pred: "<<corr<<std::endl;
      }

      // more accurate method, re-weight each contribution individually
      if(wOpt == kAccurate) {
        TH1D *hTemp = (TH1D *)inFile->Get(predFDHistos[iDau].data());
        hTemp->SetDirectory(0);
        hDauFDPred = (TH1D *)hTemp->Clone();
        hDauFDPred->Reset();

        std::string hName = predFDHistos[iDau] + "ByOrigin";
        TH2D *hDauPredByOrigin = (TH2D *)inFile->Get(hName.data());
        hDauPredByOrigin->SetDirectory(0);

        for(int iMother = 0; iMother < numMothers; iMother++) {
          TH1D *hTemp = hDauPredByOrigin->ProjectionY(Form("hBpred%d", iMother), iMother + 1, iMother + 1);

          double motherFF = 1.;
          double BRfromB = 1.;

          if(ffOpt == kFFOriginal)
            motherFF = origBFF[iMother];
          else if(ffOpt == kFFppbar)
            motherFF = ppbarBFF[iMother];
          else if(ffOpt == kFFee)
            motherFF = eeBFF[iMother];

          if(brOpt == kBROriginal)
            BRfromB = origBR[iDau][iMother];
          else if(brOpt == kBRPDG)
            BRfromB = pdgBRfromB[iDau][iMother];
          else if(brOpt == kBRPDGmix)
            BRfromB = pdgBRfromBmix[iDau][iMother];

          hDauFDPred->Add(hTemp, motherFF * BRfromB / origBFF[iMother] / origBR[iDau][iMother]);
        }
      }

      // prompt predictions
      TH1D *hDauPromptPred = (TH1D *)inFile->Get(predPromptHistos[iDau].data());
      hDauPromptPred->SetDirectory(0);
      hDauPromptPred->Scale(decayBR[iDau] / 1.e-6);
      for(int iBin = 0; iBin <= hDauPromptPred->GetXaxis()->GetNbins(); iBin++)
        hDauPromptPred->SetBinError(iBin + 1, 1.e-18);
      std::string name = "h" + predTag[iDau] + "pred_" + edgeNames[iFile];
      std::string title = predTag[iDau] + " " + edgeNames[iFile] + " value prediction (with BR)";
      hDauPromptPred->SetName(name.data());
      hDauPromptPred->SetTitle(title.data());
      hDauPromptPred->GetXaxis()->SetTitle("p_{T} (GeV)");
      hDauPromptPred->GetYaxis()->SetTitle("d#sigma/dp_{T} x BR (pb/GeV)");
      outFile.cd();
      hDauPromptPred->Write();

      // non-prompt predictions
      hDauFDPred->Scale(decayBR[iDau] / 1.e-6);
      std::string nameFD = "h" + predTag[iDau] + "fromBpred_" + edgeNames[iFile] + "_corr";
      std::string titleFD = predTag[iDau] + " from B " + edgeNames[iFile] + " value prediction (with BR and B->D correction)";
      hDauFDPred->SetName(nameFD.data());
      hDauFDPred->SetTitle(titleFD.data());
      hDauFDPred->GetYaxis()->SetTitle("d#sigma/dp_{T} x BR (pb/GeV)");
      outFile.cd();
      hDauFDPred->Write();
    }
    inFile->Close();
  }
}