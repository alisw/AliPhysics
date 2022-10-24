/**************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Joshua Koenig <joshua.konig@cern.ch>                                        *
 * Version 1.0                                                            *
 *                                                                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class to hanlde 2d and 4d response matrix
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////

#ifndef ALIRESPONSEMATRIXHELPER_h
#define ALIRESPONSEMATRIXHELPER_h

#include <array>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"

class MatrixHandler4D
{

 public:
  MatrixHandler4D() = default;
  MatrixHandler4D(std::vector<double> arrMesonX, std::vector<double> arrMesonY, std::vector<double> arrJetX, std::vector<double> arrJetY, bool useTHN = false);
  MatrixHandler4D(std::vector<double> arrMesonX, std::vector<double> arrMesonY, std::vector<double> arrJetX, std::vector<double> arrJetY, THnSparse *h = nullptr);
  MatrixHandler4D(const MatrixHandler4D&) = delete;            // copy ctor
  MatrixHandler4D(MatrixHandler4D&&) = delete;                 // move ctor
  MatrixHandler4D& operator=(const MatrixHandler4D&) = delete; // copy assignment
  MatrixHandler4D& operator=(MatrixHandler4D&&) = delete;      // move assignment
  virtual ~MatrixHandler4D();

  int getBinIndexMesonX(const double val) const;
  int getBinIndexMesonY(const double val) const;
  int getBinIndexJetX(const double val) const;
  int getBinIndexJetY(const double val) const;

  void Fill(double valJetX, double valJetY, double valMesonX, double valMesonY, double val = 1)
  {
    int binJetX = getBinIndexJetX(valJetX);
    if (binJetX < 0)
      return;
    int binJetY = getBinIndexJetY(valJetY);
    if (binJetY < 0)
      return;
    int binMesonX = getBinIndexMesonX(valMesonX);
    if (binMesonX < 0)
      return;
    int binMesonY = getBinIndexMesonY(valMesonY);
    if (binMesonY < 0)
      return;

    if (useTHNSparese) {
      std::array<double, 2> arrFill;
      arrFill[0] = binJetX * (vecBinsMesonX.size() - 1) + binMesonX + 1 - 0.5; // -0.5 due to the fact that this is not the bin but the bin center
      arrFill[1] = binJetY * (vecBinsMesonY.size() - 1) + binMesonY + 1 - 0.5;
      // hSparseResponse->AddBinContent(arrFill.data(), val); // what happens to the errors?
      hSparseResponse->Fill(arrFill.data(), val);
    } else {
      std::array<int, 2> arrFill;
      arrFill[0] = binJetX * (vecBinsMesonX.size() - 1) + binMesonX + 1; // + 1 due to underflow bin?
      arrFill[1] = binJetY * (vecBinsMesonY.size() - 1) + binMesonY + 1;
      h2d->Fill(h2d->GetXaxis()->GetBinCenter(arrFill[0]), h2d->GetYaxis()->GetBinCenter(arrFill[1]), val);
    }
  }
  void GetAxisBinning(std::vector<double>& vecXBins, std::vector<double>& vecYBins)
  {
    const int sizeX = (vecBinsJetX.size() - 1) * (vecBinsMesonX.size() - 1) + 1;
    const int sizeY = (vecBinsJetY.size() - 1) * (vecBinsMesonY.size() - 1) + 1;
    vecXBins.resize(sizeX);
    vecYBins.resize(sizeY);
    double rangeXBinsMesons = std::abs(vecBinsMesonX[0] - vecBinsMesonX.back());
    double rangeYBinsMesons = std::abs(vecBinsMesonY[0] - vecBinsMesonY.back());
    int counter = 0;
    for (unsigned int jx = 0; jx < vecBinsJetX.size() - 1; ++jx) {
      double rangeXBinsJet = std::abs(vecBinsJetX[jx] - vecBinsJetX[jx + 1]);
      double downscalingfacX = rangeXBinsJet / rangeXBinsMesons;
      for (unsigned int mx = 0; mx < vecBinsMesonX.size() - 1; ++mx) {
        vecXBins[counter] = vecBinsJetX[jx] + vecBinsMesonX[mx] * downscalingfacX;
        counter++;
      }
    }
    vecXBins[counter] = vecBinsJetX.back();
    for (unsigned int i = 0; i < vecXBins.size() - 1; ++i) {
      if (vecXBins[i] >= vecXBins[i + 1])
        printf("Binning is not in ascending order! FIX THIS!!\n");
    }
    counter = 0;
    for (unsigned int jy = 0; jy < vecBinsJetY.size() - 1; ++jy) {
      double rangeYBinsJet = std::abs(vecBinsJetY[jy] - vecBinsJetY[jy + 1]);
      double downscalingfacY = rangeYBinsJet / rangeYBinsMesons;
      for (unsigned int my = 0; my < vecBinsMesonY.size() - 1; ++my) {
        vecYBins[counter] = vecBinsJetY[jy] + vecBinsMesonY[my] * downscalingfacY;
        counter++;
      }
    }
    vecYBins[counter] = vecBinsJetY.back();
  }

  THnSparseF* GetTHnSparseClone(const char* name = "hSparseResponse_Clone")
  {
    return (THnSparseF*)hSparseResponse->Clone(name);
  }
  THnSparseF* GetTHnSparse(const char* name = "")
  {
    hSparseResponse->SetName(name);
    return hSparseResponse;
  }
  TH2F* GetTH2(const char* name = "hSparseResponse_Clone")
  {
    if (useTHNSparese) {
      if (!hSparseResponse) {
        printf("Attention! hSparseResponse does not exist yet!\n");
        return nullptr;
      }
      std::vector<double> vecXBins;
      std::vector<double> vecYBins;
      GetAxisBinning(vecXBins, vecYBins);

      if (h2d) {
        delete h2d;
      }
      h2d = new TH2F(name, name, vecXBins.size() - 1, vecXBins.data(), vecYBins.size() - 1, vecYBins.data());

      // transfer the values
      for (int x = 0; x < h2d->GetXaxis()->GetNbins(); ++x) {
        for (int y = 0; y < h2d->GetYaxis()->GetNbins(); ++y) {
          std::array<int, 2> arrBins = {x, y};
          double binCont = hSparseResponse->GetBinContent(arrBins.data());
          double binErr = hSparseResponse->GetBinError(arrBins.data());
          h2d->SetBinContent(x, y, binCont);
          h2d->SetBinError(x, y, binErr);
        }
      }
      return h2d;

    } else {
      if (!h2d) {
        printf("Attention! 2d histogram does not exist yet!\n");
        return nullptr;
      }
      h2d->SetName(name);
      return h2d;
    }
  }

  TH2F* GetResponseMatrix(int binX, int binY, const char* name = "dummy")
  {
    TString nameHist = Form("%s_%i_%i", name, binX, binY);
    TH2F* hMesonResp = new TH2F(nameHist, nameHist, vecBinsMesonX.size() - 1, vecBinsMesonX.data(), vecBinsMesonY.size() - 1, vecBinsMesonY.data());
    for (unsigned int x = 0; x < vecBinsMesonX.size(); ++x) {
      for (unsigned int y = 0; y < vecBinsMesonY.size(); ++y) {
        double binCont = 0;
        double binErr = 0;
        std::array<int, 2> arrBins = {static_cast<int>(x + (binX * (vecBinsMesonX.size() - 1)) + 1), static_cast<int>(y + (binY * (vecBinsMesonY.size() - 1)) + 1)};
        if (useTHNSparese) {
          binCont = hSparseResponse->GetBinContent(arrBins.data());
          binErr = hSparseResponse->GetBinError(hSparseResponse->GetBin(arrBins.data()));
        } else {
          binCont = h2d->GetBinContent(arrBins[0], arrBins[1]);
          binErr = h2d->GetBinError(arrBins[0], arrBins[1]);
        }
        hMesonResp->SetBinContent(x + 1, y + 1, binCont);
        hMesonResp->SetBinError(x + 1, y + 1, binErr);
      }
    }
    return hMesonResp;
  }

 private:
  bool useTHNSparese = false;
  int nBinsJet = 0;
  std::vector<double> vecBinsMesonX = {};
  std::vector<double> vecBinsMesonY = {};
  std::vector<double> vecBinsJetX = {};
  std::vector<double> vecBinsJetY = {};
  TH2F* h2d = nullptr;
  TH1F* h1dJet = nullptr;
  TH1F* h1dMeson = nullptr;
  THnSparseF* hSparseResponse = nullptr;

  ClassDef(MatrixHandler4D, 2)
};

#endif