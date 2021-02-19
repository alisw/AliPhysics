// Macro to "cook" the FONLL + Pythia predictions of D-meson and Lc cross sections obtained from ComputeBtoDdecay.C .
//
// The BR and FF are set to the DPG values and can be changed. The output has the same structure and histogram names
// as the one used in D2H analyses.
//
// Author: F. Catalano, fabio.catalano@cern.ch
//         F. Grosa, fabrizio.grosa@cern.ch
//

#include <string>
#include <array>
#include <iostream>
#include <numeric>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

enum { // options for branching ratios
  kBROriginal, // keep the values used in the input files
  kBRPDG, // use hard coded values from PDG (2020)
  kBRPDGmix // use hard coded values from PDG (2020) and the B0/B+- admixture
};

enum { // options for fragmentation fractions
  kFFOriginal, // keep the values used in the input files (mandatory for simulations produced with LHCb pT-dependent FF)
  kFFppbar, // use hard coded values from ppbar, PDG (2020)
  kFFee  // use hard coded values from e+e-, PDG (2020)
};

enum { // options for re-weight
  kSimple, // estimate a global factor to scale the prediction (keep the input pT dependence)
  kAccurate // correct each b-hadron contribution indipendently (modify also the pT dependence)
};

void CookFONLLPythiaPred(std::string inFileNameMin = "DfromB_FONLLminPythia8_FFptDepmin_yDcut_pp13TeV_PDG2020.root", // min FONLL predictions
                         std::string inFileNameCent = "DfromB_FONLLcentPythia8_FFptDepcent_yDcut_pp13TeV_PDG2020.root", // central FONLL predictions
                         std::string inFileNameMax = "DfromB_FONLLmaxPythia8_FFptDepmax_yDcut_pp13TeV_PDG2020.root",  // max FONLL predictions
                         std::string outFileName = "DmesonLcPredictions_13TeV_y05_FFptDepLHCb_BRpythia8_PDG2020.root",
                         int brOpt = kBROriginal,
                         int ffOpt = kFFOriginal,
                         int wOpt = kAccurate) {

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
  std::array<std::string, numDaughters - 1> partTitle = {"D^{0}", "D^{+}", "D_{s}^{+}", "#Lambda_{c}^{+}", "D^{*+}"};
  std::array<std::string, numDaughters - 1> partBTitle = {"B^{0}", "B^{+}", "B_{s}^{0}", "#Lambda_{b}^{0}"};
  std::array<std::array<double, numMothers>, numDaughters> pdgBRfromB = {{{0.555, 0.876, 0.002, 0.},   // D0 and (BRfromB0, BRfromB+, BRfromBs, BRfromLb) from PDG (2020)
                                                                          {0.392, 0.124, 0., 0.},      // D+
                                                                          {0.117, 0.09, 0.93, 0.011},  // Ds
                                                                          {0.066, 0.049, 0., 0.352},   // Lc
                                                                          {0.233, 0.055, 0.003, 0.},    // D*+
                                                                          {0.066, 0.049, 0., 0.352}    // Lc
                                                                          }};
  // B0 and B+ BR from B0/B+- admixture
  std::array<std::array<double, numMothers>, numDaughters> pdgBRfromBmix = {{{0.624, 0.624, 0.002, 0.},   // D0 and (BRfromB0, BRfromB+, BRfromBs, BRfromLb) from PDG (2020) and B0/B+- admixture (2020)
                                                                             {0.241, 0.241, 0., 0.},      // D+
                                                                             {0.083, 0.083, 0.93, 0.011}, // Ds
                                                                             {0.036, 0.036, 0., 0.352},   // Lc
                                                                             {0.225, 0.225, 0.003, 0.},   // D*+
                                                                             {0.036, 0.036, 0., 0.352}    // Lc
                                                                             }};
  std::array<double, numMothers> ppbarBFF = {0.344, 0.344, 0.115, 0.198}; // (B0, B+, Bs, Lb) from PDG (2020)
  std::array<double, numMothers> eeBFF = {0.408, 0.408, 0.100, 0.084}; // (B0, B+, Bs, Lb) from PDG (2020)
  std::array<double, numDaughters> decayBR = {0.03950, 0.0938, 0.0224, 0.0628, 0.03950*0.677, 0.0159*0.692}; // (D0, D+, Ds, Lc->pKpi, D*+, LC->K0sp) from PDG (2020)
  std::array<std::array<double, numMothers>, numDaughters> origBR = {};
  std::array<double, numMothers> origBFF = {};

  TFile outFile(outFileName.data(),"recreate");

  for(int iFile = 0; iFile < 3; iFile++) {
    TFile *inFile = TFile::Open(inFileNames[iFile].data());

    // get the original FF
    TH1D *hFFBeauty = (TH1D *)inFile->Get("hbFragmFrac");
    for(int iMother = 0; iMother < numMothers; iMother++) {
      origBFF[iMother] = hFFBeauty->GetBinContent(iMother + 1);
    }

    double origBFFsum = std::accumulate(origBFF.begin(), origBFF.end(), 0.f);
    if(ffOpt != kFFOriginal && origBFFsum == 0) {
      std::cerr << "ERROR: only original FF can be set if pT-dependent LHCb parametrisation used in the simulation! Exit" << std::endl;
      return;
    }

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

    if(edgeNames[iFile] == "central") {
      if(brOpt == kBROriginal) {
        // print original BR factors from file, only for central case
        for(int iMother = 0; iMother < numMothers; iMother++) {
          std::cout<<"BR for "<<mothToDauHistos[iMother]<<" extracted from "<<inFileNames[iFile]<<" file\n";
          for(int iDau = 0; iDau < numDaughters - 1; iDau++)
            std::cout<<"  to  " + partTag[iDau] << ": "<<origBR[iDau][iMother]<<std::endl;
        }

        // print fraction of b to b-hadrons FF times BR extracted from file, only for central case
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
      else if(brOpt == kBRPDG || brOpt == kBRPDGmix) {
      // print fraction of b to b-hadrons FF times PDG BR
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
    }

    // get and correct the predictions
    for(int iDau = 0; iDau < numDaughters; iDau++) {
      TH1D *hDauFDPred = nullptr;
      TH1D *hDauFDPredFromBs = nullptr;
      TH1D *hDauFDPredFromNonStrangeB = nullptr;
      TH1D *hDauFDPredFromLb = nullptr;

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

          if(origBFF[iMother] > 0) {
            newFrac += motherFF * BRfromB;
            oldFrac += origBFF[iMother] * origBR[iDau][iMother];
          }
          else {// pt-dependent FF in original file --> always unmodified
            newFrac += BRfromB;
            oldFrac += origBR[iDau][iMother];
          }
        }

        double corr = 0.;
        if(oldFrac > 0)
          corr = newFrac / oldFrac;

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
        hDauFDPredFromBs = (TH1D *)hTemp->Clone();
        hDauFDPredFromBs->Reset();
        hDauFDPredFromLb = (TH1D *)hTemp->Clone();
        hDauFDPredFromLb->Reset();
        hDauFDPredFromNonStrangeB = (TH1D *)hTemp->Clone();
        hDauFDPredFromNonStrangeB->Reset();

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

          double corr = 0.;
          if(origBR[iDau][iMother] > 0 && origBFF[iMother] > 0)
            corr = motherFF * BRfromB / origBFF[iMother] / origBR[iDau][iMother];
          else if(origBR[iDau][iMother] > 0) // pt-dependent FF in original file --> always unmodified
            corr = BRfromB / origBR[iDau][iMother];
        
          hDauFDPred->Add(hTemp, corr);
          if(iMother < 2)
            hDauFDPredFromNonStrangeB->Add(hTemp, corr);
          else if(iMother == 2)
            hDauFDPredFromBs->Add(hTemp, corr);
          else
            hDauFDPredFromLb->Add(hTemp, corr);
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
      hDauFDPredFromNonStrangeB->Scale(decayBR[iDau] / 1.e-6);
      hDauFDPredFromBs->Scale(decayBR[iDau] / 1.e-6);
      hDauFDPredFromLb->Scale(decayBR[iDau] / 1.e-6);
      std::string nameFD = "h" + predTag[iDau] + "fromBpred_" + edgeNames[iFile] + "_corr";
      std::string titleFD = predTag[iDau] + " from B " + edgeNames[iFile] + " value prediction (with BR and B->D correction)";
      hDauFDPred->SetName(nameFD.data());
      hDauFDPred->SetTitle(titleFD.data());
      hDauFDPred->GetYaxis()->SetTitle("d#sigma/dp_{T} x BR (pb/GeV)");
      nameFD = "h" + predTag[iDau] + "fromB0Bpluspred_" + edgeNames[iFile] + "_corr";
      titleFD = predTag[iDau] + " from B0, B+ " + edgeNames[iFile] + " value prediction (with BR and B->D correction)";
      hDauFDPredFromNonStrangeB->SetName(nameFD.data());
      hDauFDPredFromNonStrangeB->SetTitle(titleFD.data());
      hDauFDPredFromNonStrangeB->GetYaxis()->SetTitle("d#sigma/dp_{T} x BR (pb/GeV)");
      nameFD = "h" + predTag[iDau] + "fromBspred_" + edgeNames[iFile] + "_corr";
      titleFD = predTag[iDau] + " from Bs " + edgeNames[iFile] + " value prediction (with BR and B->D correction)";
      hDauFDPredFromBs->SetName(nameFD.data());
      hDauFDPredFromBs->SetTitle(titleFD.data());
      hDauFDPredFromBs->GetYaxis()->SetTitle("d#sigma/dp_{T} x BR (pb/GeV)");
      nameFD = "h" + predTag[iDau] + "fromLbpred_" + edgeNames[iFile] + "_corr";
      titleFD = predTag[iDau] + " from Lb " + edgeNames[iFile] + " value prediction (with BR and B->D correction)";
      hDauFDPredFromLb->SetName(nameFD.data());
      hDauFDPredFromLb->SetTitle(titleFD.data());
      hDauFDPredFromLb->GetYaxis()->SetTitle("d#sigma/dp_{T} x BR (pb/GeV)");
      outFile.cd();
      hDauFDPred->Write();
      hDauFDPredFromNonStrangeB->Write();
      hDauFDPredFromBs->Write();
      hDauFDPredFromLb->Write();
    }
    inFile->Close();
  }

  // save also BRs and FFs used
  for(int iDau = 0; iDau < numDaughters-1; iDau++) {
    TH1D *hDauBR = new TH1D(Form("hBRHbto%s", partTag[iDau].data()), ";;BR", 4, 0.5, 4.5);
    for(int iMother = 0; iMother < numMothers; iMother++) {
      double BRfromB = 1.;
      if(brOpt == kBROriginal)
        BRfromB = origBR[iDau][iMother];
      else if(brOpt == kBRPDG)
        BRfromB = pdgBRfromB[iDau][iMother];
      else if(brOpt == kBRPDGmix)
        BRfromB = pdgBRfromBmix[iDau][iMother];
      hDauBR->GetXaxis()->SetBinLabel(iMother+1, Form("%s #rightarrow %s", partBTitle[iMother].data(), partTitle[iDau].data()));
      hDauBR->SetBinContent(iMother+1, BRfromB);
    }
    outFile.cd();
    hDauBR->Write();
  }

  TH1D *hMotherFF = new TH1D("hHbFF", ";;FF", 4, 0.5, 4.5);
  for(int iMother = 0; iMother < numMothers; iMother++) {
    double motherFF = 1.;
    if(ffOpt == kFFOriginal)
      motherFF = origBFF[iMother];
    else if(ffOpt == kFFppbar)
      motherFF = ppbarBFF[iMother];
    else if(ffOpt == kFFee)
      motherFF = eeBFF[iMother];
    hMotherFF->GetXaxis()->SetBinLabel(iMother+1, Form("b #rightarrow %s", partBTitle[iMother].data()));
    hMotherFF->SetBinContent(iMother+1, motherFF);
  }
  outFile.cd();
  hMotherFF->Write();
  outFile.Close();
}