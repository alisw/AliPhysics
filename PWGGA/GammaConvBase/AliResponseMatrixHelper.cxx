
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

#include "AliResponseMatrixHelper.h"

//____________________________________________________________________________________________________________________________
MatrixHandler4D::MatrixHandler4D(std::vector<double> arrMesonX, std::vector<double> arrMesonY, std::vector<double> arrJetX, std::vector<double> arrJetY, bool useTHN)
{
  vecBinsMesonX = arrMesonX;
  vecBinsMesonY = arrMesonY;
  vecBinsJetX = arrJetX;
  vecBinsJetY = arrJetY;
  useTHNSparese = useTHN;
  if (useTHNSparese) {
    // in case of thnsparse, just use equidistant binning
    const int nBinsX = (arrMesonX.size() - 1) * (arrJetX.size() - 1);
    const int nBinsY = (arrMesonY.size() - 1) * (arrJetY.size() - 1);
    std::array<int, 2> arrNBins = {nBinsX, nBinsY};
    std::array<double, 2> arrXBins = {0, 0};
    std::array<double, 2> arrYBins = {static_cast<double>(nBinsX), static_cast<double>(nBinsY)};

    hSparseResponse = new THnSparseF("ResponseMatrix_dyn", "ResponseMatrix_dyn", arrNBins.size(), arrNBins.data(), arrXBins.data(), arrYBins.data());

  } else {
    std::vector<double> vecXBins;
    std::vector<double> vecYBins;
    GetAxisBinning(vecXBins, vecYBins);
    if (h2d) {
      delete h2d;
    }
    h2d = new TH2F("ResponseMatrix_stat", "ResponseMatrix_stat", vecXBins.size() - 1, vecXBins.data(), vecYBins.size() - 1, vecYBins.data());
  }
}

//____________________________________________________________________________________________________________________________
MatrixHandler4D::MatrixHandler4D(std::vector<double> arrMesonX, std::vector<double> arrMesonY, std::vector<double> arrJetX, std::vector<double> arrJetY, THnSparse* h)
{
  vecBinsMesonX = arrMesonX;
  vecBinsMesonY = arrMesonY;
  vecBinsJetX = arrJetX;
  vecBinsJetY = arrJetY;
  useTHNSparese = true;
  hSparseResponse = (THnSparseF*)h->Clone();
}

//____________________________________________________________________________________________________________________________
MatrixHandler4D::MatrixHandler4D(std::vector<double> arrMesonX, std::vector<double> arrMesonY, std::vector<double> arrJetX, std::vector<double> arrJetY, TH2F* h)
{
  vecBinsMesonX = arrMesonX;
  vecBinsMesonY = arrMesonY;
  vecBinsJetX = arrJetX;
  vecBinsJetY = arrJetY;
  useTHNSparese = false;
  h2d = (TH2F*)h->Clone();
}

//____________________________________________________________________________________________________________________________
MatrixHandler4D::~MatrixHandler4D()
{
}

//____________________________________________________________________________________________________________________________
int MatrixHandler4D::getBinIndexMesonX(const double val) const
{
  for (unsigned int i = 0; i < vecBinsMesonX.size() - 1; ++i) {
    if (vecBinsMesonX[i] < val && vecBinsMesonX[i + 1] > val) {
      return i;
    }
  }
  return -1;
}

//____________________________________________________________________________________________________________________________
int MatrixHandler4D::getBinIndexMesonY(const double val) const
{
  for (unsigned int i = 0; i < vecBinsMesonY.size() - 1; ++i) {
    if (vecBinsMesonY[i] < val && vecBinsMesonY[i + 1] > val) {
      return i;
    }
  }
  return -1;
}

//____________________________________________________________________________________________________________________________
int MatrixHandler4D::getBinIndexJetX(const double val) const
{
  for (unsigned int i = 0; i < vecBinsJetX.size() - 1; ++i) {
    if (vecBinsJetX[i] < val && vecBinsJetX[i + 1] > val) {
      return i;
    }
  }
  return -1;
}

//____________________________________________________________________________________________________________________________
int MatrixHandler4D::getBinIndexJetY(const double val) const
{
  for (unsigned int i = 0; i < vecBinsJetY.size() - 1; ++i) {
    if (vecBinsJetY[i] < val && vecBinsJetY[i + 1] > val) {
      return i;
    }
  }
  return -1;
}

//____________________________________________________________________________________________________________________________
double MatrixHandler4D::getValueForBinIndexJetX(const int index) const
{
  const int indexJet = static_cast<int>(index / (vecBinsMesonX.size() - 1));
  return 0.5 * (vecBinsJetX[indexJet] + vecBinsJetX[indexJet + 1]);
}

//____________________________________________________________________________________________________________________________
double MatrixHandler4D::getValueForBinIndexMesonX(const int index) const
{
  const int indexMeson = index % (vecBinsMesonX.size() - 1);
  return 0.5 * (vecBinsMesonX[indexMeson] + vecBinsMesonX[indexMeson + 1]);
}

//____________________________________________________________________________________________________________________________
double MatrixHandler4D::getValueForBinIndexJetY(const int index) const
{
  const int indexJet = static_cast<int>(index / (vecBinsMesonY.size() - 1));
  return 0.5 * (vecBinsJetY[indexJet] + vecBinsJetY[indexJet + 1]);
}

//____________________________________________________________________________________________________________________________
double MatrixHandler4D::getValueForBinIndexMesonY(const int index) const
{
  const int indexMeson = index % (vecBinsMesonY.size() - 1);
  return 0.5 * (vecBinsMesonY[indexMeson] + vecBinsMesonY[indexMeson + 1]);
}

//____________________________________________________________________________________________________________________________
void MatrixHandler4D::Fill(double valJetX, double valJetY, double valMesonX, double valMesonY, double val)
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
    hSparseResponse->Fill(arrFill.data(), val);
  } else {
    std::array<int, 2> arrFill;
    arrFill[0] = binJetX * (vecBinsMesonX.size() - 1) + binMesonX + 1;
    arrFill[1] = binJetY * (vecBinsMesonY.size() - 1) + binMesonY + 1;
    h2d->Fill(h2d->GetXaxis()->GetBinCenter(arrFill[0]), h2d->GetYaxis()->GetBinCenter(arrFill[1]), val);
  }
}

//____________________________________________________________________________________________________________________________
void MatrixHandler4D::AddBinContent(double valJetX, double valJetY, double valMesonX, double valMesonY, double val, double err)
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
    hSparseResponse->Fill(arrFill.data(), val);
  } else {
    std::array<int, 2> arrFill;
    arrFill[0] = binJetX * (vecBinsMesonX.size() - 1) + binMesonX + 1;
    arrFill[1] = binJetY * (vecBinsMesonY.size() - 1) + binMesonY + 1;
    double currentBinCont = h2d->GetBinContent(arrFill[0], arrFill[1]);
    double currentBinErr = h2d->GetBinError(arrFill[0], arrFill[1]);
    h2d->SetBinContent(arrFill[0], arrFill[1], val + currentBinCont);
    h2d->SetBinError(arrFill[0], arrFill[1], std::sqrt(err * err + currentBinErr * currentBinErr));
  }
}

//____________________________________________________________________________________________________________________________
void MatrixHandler4D::GetAxisBinning(std::vector<double>& vecXBins, std::vector<double>& vecYBins)
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

//____________________________________________________________________________________________________________________________
THnSparseF* MatrixHandler4D::GetTHnSparseClone(const char* name)
{
  return (THnSparseF*)hSparseResponse->Clone(name);
}

//____________________________________________________________________________________________________________________________
THnSparseF* MatrixHandler4D::GetTHnSparse(const char* name)
{
  hSparseResponse->SetName(name);
  return hSparseResponse;
}

//____________________________________________________________________________________________________________________________
TH2F* MatrixHandler4D::GetTH2(const char* name)
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
        std::array<int, 2> arrBins = {x + 1, y + 1}; // underflow bin?
        double binCont = hSparseResponse->GetBinContent(arrBins.data());
        double binErr = hSparseResponse->GetBinError(arrBins.data());
        h2d->SetBinContent(x + 1, y + 1, binCont);
        h2d->SetBinError(x + 1, y + 1, binErr);
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

//____________________________________________________________________________________________________________________________
TH2F* MatrixHandler4D::GetResponseMatrix(int binX, int binY, const char* name)
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