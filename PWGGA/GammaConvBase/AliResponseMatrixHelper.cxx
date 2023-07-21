
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

// //____________________________________________________________________________________________________________________________
// MatrixHandler4D::~MatrixHandler4D()
// {
// }

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
      vecXBins[counter] = vecBinsJetX[jx] + (vecBinsMesonX[mx] - vecBinsMesonX[0]) * downscalingfacX;
      if (counter > 0) {
        if (vecXBins[counter] <= vecXBins[counter - 1])
          printf("Binning is not in ascending order! FIX THIS!!\n");
      }
      counter++;
    }
  }
  vecXBins[counter] = vecBinsJetX.back();

  counter = 0;
  for (unsigned int jy = 0; jy < vecBinsJetY.size() - 1; ++jy) {
    double rangeYBinsJet = std::abs(vecBinsJetY[jy] - vecBinsJetY[jy + 1]);
    double downscalingfacY = rangeYBinsJet / rangeYBinsMesons;
    for (unsigned int my = 0; my < vecBinsMesonY.size() - 1; ++my) {
      vecYBins[counter] = vecBinsJetY[jy] + (vecBinsMesonY[my] - vecBinsMesonY[0]) * downscalingfacY;
      if (counter > 0) {
        if (vecYBins[counter] <= vecYBins[counter - 1])
          printf("Binning is not in ascending order! FIX THIS!!\n");
      }
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

//____________________________________________________________________________________________________________________________
void MatrixHandler4D::WeightResponseMatrix(TF1* funcMeson, TF1* funcJet)
{

  for(unsigned int jx = 0; jx < vecBinsJetX.size()-1; ++jx){
    for(unsigned int jy = 0; jy < vecBinsJetY.size()-1; ++jy){
      for(unsigned int mx = 0; mx < vecBinsMesonX.size()-1; ++mx){
        for(unsigned int my = 0; my < vecBinsMesonY.size()-1; ++my){
          const double valJet = 0.5*(vecBinsJetY[jy] + vecBinsJetY[jy+1]);
          const double valMeson = 0.5*(vecBinsMesonY[my] + vecBinsMesonY[my+1]);
          const double weightMeson = (funcMeson == nullptr) ? 1 : funcMeson->Eval(valMeson);
          const double weightJet = (funcJet == nullptr) ? 1 : funcJet->Eval(valJet);
          const double weight = weightMeson*weightJet;
          
          if (useTHNSparese) {
            std::array<double, 2> arrFill;
            arrFill[0] = jx * (vecBinsMesonX.size() - 1) + mx + 1 - 0.5; // -0.5 due to the fact that this is not the bin but the bin center
            arrFill[1] = jy * (vecBinsMesonY.size() - 1) + my + 1 - 0.5;
            const double val = hSparseResponse->GetBinContent(hSparseResponse->GetBin(arrFill.data()));
            const double err = hSparseResponse->GetBinError(hSparseResponse->GetBin(arrFill.data()));
            hSparseResponse->SetBinContent(hSparseResponse->GetBin(arrFill.data()), val*weight);
            hSparseResponse->SetBinError(hSparseResponse->GetBin(arrFill.data()), err*weight);
          } else {
            std::array<int, 2> arrFill;
            arrFill[0] = jx * (vecBinsMesonX.size() - 1) + mx + 1;
            arrFill[1] = jy * (vecBinsMesonY.size() - 1) + my + 1;
            const double val = h2d->GetBinContent(arrFill[0], arrFill[1]);
            const double err = h2d->GetBinError(arrFill[0], arrFill[1]);
            h2d->SetBinContent(arrFill[0], arrFill[1], val*weight);
            h2d->SetBinError(arrFill[0], arrFill[1], err*weight);
          }
        }
      }
    }
  }
}

//____________________________________________________________________________________________________________________________
void MatrixHandler4D::WeightResponseMatrix(TF2* func)
{

  for(unsigned int jx = 0; jx < vecBinsJetX.size()-1; ++jx){
    for(unsigned int jy = 0; jy < vecBinsJetY.size()-1; ++jy){
      for(unsigned int mx = 0; mx < vecBinsMesonX.size()-1; ++mx){
        for(unsigned int my = 0; my < vecBinsMesonY.size()-1; ++my){
          const double valJet = 0.5*(vecBinsJetY[jy] + vecBinsJetY[jy+1]);
          const double valMeson = 0.5*(vecBinsMesonY[my] + vecBinsMesonY[my+1]);
          const double weight = func->Eval(valJet, valMeson);
          
          if (useTHNSparese) {
            std::array<double, 2> arrFill;
            arrFill[0] = jx * (vecBinsMesonX.size() - 1) + mx + 1 - 0.5; // -0.5 due to the fact that this is not the bin but the bin center
            arrFill[1] = jy * (vecBinsMesonY.size() - 1) + my + 1 - 0.5;
            const double val = hSparseResponse->GetBinContent(hSparseResponse->GetBin(arrFill.data()));
            const double err = hSparseResponse->GetBinError(hSparseResponse->GetBin(arrFill.data()));
            hSparseResponse->SetBinContent(hSparseResponse->GetBin(arrFill.data()), val*weight);
            hSparseResponse->SetBinError(hSparseResponse->GetBin(arrFill.data()), err*weight);
            
          } else {
            std::array<int, 2> arrFill;
            arrFill[0] = jx * (vecBinsMesonX.size() - 1) + mx + 1;
            arrFill[1] = jy * (vecBinsMesonY.size() - 1) + my + 1;
            const double val = h2d->GetBinContent(arrFill[0], arrFill[1]);
            const double err = h2d->GetBinError(arrFill[0], arrFill[1]);
            h2d->SetBinContent(arrFill[0], arrFill[1], val*weight);
            h2d->SetBinError(arrFill[0], arrFill[1], err*weight);
          }
            
        }
      }
    }
  }
}






//=============================================================//
//============ Implementation of N dimensional handler ========//
//=============================================================//
//____________________________________________________________________________________________________________________________
MatrixHandlerNDim::MatrixHandlerNDim(std::vector<std::vector<double>> arrBinsX, std::vector<std::vector<double>> arrBinsY, bool useTHN)
{
  if(arrBinsX.size() != arrBinsY.size()){
    printf("MatrixHandlerNDim::MatrixHandlerNDim: ERROR. Vectors of different size\n");
  }

  vecBinsX = arrBinsX;
  vecBinsY = arrBinsY;
  useTHNSparese = useTHN;
  if (useTHNSparese) {
    // in case of thnsparse, just use equidistant binning

    std::vector<double> vecXBins = MakeAxis1D(arrBinsX);
    std::vector<double> vecYBins = MakeAxis1D(arrBinsY);
    int nBinsX = vecXBins.size()-1;
    int nBinsY = vecYBins.size()-1;
    
    std::array<int, 2> arrNBins = {nBinsX, nBinsY};
    std::array<double, 2> arrXBins = {0, 0};
    std::array<double, 2> arrYBins = {static_cast<double>(nBinsX), static_cast<double>(nBinsY)};

    hSparseResponse = new THnSparseF("ResponseMatrix_dyn", "ResponseMatrix_dyn", arrNBins.size(), arrNBins.data(), arrXBins.data(), arrYBins.data());

  } else {
    std::vector<double> vecXBins = MakeAxis1D(arrBinsX);
    std::vector<double> vecYBins = MakeAxis1D(arrBinsY);
    if (h2d) {
      delete h2d;
    }
    h2d = new TH2F("ResponseMatrix_stat", "ResponseMatrix_stat", vecXBins.size() - 1, vecXBins.data(), vecYBins.size() - 1, vecYBins.data());
  }
}


//____________________________________________________________________________________________________________________________
MatrixHandlerNDim::MatrixHandlerNDim(std::vector<std::vector<double>> arrX, std::vector<std::vector<double>> arrY, THnSparse* h)
{
  vecBinsX = arrX;
  vecBinsY = arrY;
  useTHNSparese = true;
  hSparseResponse = (THnSparseF*)h->Clone();
}

//____________________________________________________________________________________________________________________________
MatrixHandlerNDim::MatrixHandlerNDim(std::vector<std::vector<double>> arrX, std::vector<std::vector<double>> arrY, TH2F* h)
{
  vecBinsX = arrX;
  vecBinsY = arrY;
  useTHNSparese = false;
  h2d = (TH2F*)h->Clone();
}


//____________________________________________________________________________________________________________________________
std::vector<double> MatrixHandlerNDim::MakeAxis1D(const std::vector<std::vector<double>> vecND){
  const unsigned int dim = vecND.size();

    std::vector<unsigned int> vecDim(dim, 0);

    unsigned int nBins = 1;
    for(const auto & v : vecND){
        nBins*=(v.size()-1);
    }
    std::vector<double> vecBins(nBins+1);

    for(unsigned int i = 0; i < nBins-1; ++i){
        if(i==0) vecBins[i] = GetBinLowEdge(vecDim, vecND);
        for(int d = dim-1; d >= 0; --d){
            if(vecDim[d] == vecND[d].size()-2){
                vecDim[d] = 0;
                continue;
            }
            if(vecDim[d] < vecND[d].size()-2) {
                vecDim[d]++;
                break;
            }            
        }
        vecBins[i+1] =  GetBinLowEdge(vecDim, vecND);        
    }
    vecBins.back()=vecND[0].back();
    // printf("nBins %d\n", nBins);
    // for(const auto & i : vecBins){
    //   printf("%f, ", i);
    // }
    // printf("\n");

    return vecBins;
}


//____________________________________________________________________________________________________________________________
double MatrixHandlerNDim::GetBinLowEdge(std::vector<unsigned int> vecDim, const std::vector<std::vector<double>> vecND){
    std::vector<std::vector<double>> vecNDTmp = vecND;
    std::vector<unsigned int> vecDimTmp = vecDim;

    const unsigned int nDim = vecND.size()-1;
    double scaleFac = 1;
    // shrink first vector into specified second one
    for(unsigned int d = 1; d <= nDim; ++d){
        double dist = std::abs(vecND[d-1][vecDim[d-1]] - vecND[d-1][vecDim[d-1]+1]);
        double SizeVec = std::abs(vecND[d][0] - vecND[d].back());
        scaleFac *= dist/SizeVec;
        // cout << "scaleFac " << scaleFac << "   dist " << dist << "  SizeVec " << SizeVec << endl;
        for(unsigned int i = 0; i < vecND[d].size(); ++i){
            // cout << scaleFac*(vecND[d][i] - vecND[d][0]) << ", ";
            vecNDTmp[d][i] = scaleFac*(vecND[d][i] - vecND[d][0]);
        }
        // cout << endl;
    }

    double binCont = 0;
    for(unsigned int d = 0; d <= nDim; ++d){
        if(d == 0) binCont+=vecND[d][vecDim[d]];
        else binCont+=(vecNDTmp[d][vecDim[d]] - vecNDTmp[d][0]);
    }

    // cout << "----> " << binCont << endl;
    return binCont;  
}

//____________________________________________________________________________________________________________________________
std::vector<double> MatrixHandlerNDim::getValueForBinIndex(const unsigned int index, bool isXAxis){
  std::vector<std::vector<double>> vecBins = isXAxis == true ? vecBinsX : vecBinsY;
  const unsigned int dim = vecBins.size();

  std::vector<unsigned int> vecDim(dim, 0);

  unsigned int nBins = 1;
  for(const auto & v : vecBins){
      nBins*=(v.size()-1);
  }

  for(unsigned int i = 0; i < nBins; ++i){
    if(i == index) break;
    for(int d = dim-1; d >= 0; --d){
      if(vecDim[d] == vecBins[d].size()-2){
          vecDim[d] = 0;
          continue;
      }
      if(vecDim[d] < vecBins[d].size()-2) {
          vecDim[d]++;
          break;
      }            
    }
  }

  // get the values
  std::vector<double> val(dim, 0);
  for(unsigned int i = 0; i < vecBins.size(); ++i){
    val[i] = vecBins[i][vecDim[i]] + 0.000001;
  }

  return val;

}

//____________________________________________________________________________________________________________________________
unsigned int MatrixHandlerNDim::getIndex(std::vector<double> vecVal, bool isXAxis){
  std::vector<std::vector<double>> vecBins = isXAxis == true ? vecBinsX : vecBinsY;
  const unsigned int dim = vecBins.size();

  std::vector<unsigned int> vecDim(dim, 0);

  unsigned int nBins = 1;
  for(const auto & v : vecBins){
      nBins*=(v.size()-1);
  }

  unsigned int index = 0;
  for(unsigned int i = 0; i < nBins; ++i){
    bool indexFound = true;
    index = i;
    for(int d = dim-1; d >= 0; --d){
      if(!(vecVal[d] > vecBins[d][vecDim[d]] && vecVal[d] < vecBins[d][vecDim[d]+1])){
        indexFound = false;
      }
    }
    if(indexFound) break;
    for(int d = dim-1; d >= 0; --d){
      if(vecDim[d] == vecBins[d].size()-2){
          vecDim[d] = 0;
          continue;
      }
      if(vecDim[d] < vecBins[d].size()-2) {
          vecDim[d]++;
          break;
      }            
    }
  }
  return index;
}

//____________________________________________________________________________________________________________________________
void MatrixHandlerNDim::Fill(std::vector<double> vecValX, std::vector<double> vecValY, double val)
{
  int indexX = getIndex(vecValX, true);
  
  int indexY = getIndex(vecValY, false);
  

  if (useTHNSparese) {
    std::array<double, 2> arrFill;
    arrFill[0] = indexX + 0.5;
    arrFill[1] = indexY + 0.5;
    hSparseResponse->Fill(arrFill.data(), val);
  } else {
    std::array<int, 2> arrFill;
    arrFill[0] = indexX;
    arrFill[1] = indexY;
    h2d->Fill(h2d->GetXaxis()->GetBinCenter(arrFill[0]), h2d->GetYaxis()->GetBinCenter(arrFill[1]), val);
  }
}

//____________________________________________________________________________________________________________________________
void MatrixHandlerNDim::AddBinContent(std::vector<double> vecValX, std::vector<double> vecValY, double val, double err)
{
  int indexX = getIndex(vecValX, true);
  
  int indexY = getIndex(vecValY, false);
  

  if (useTHNSparese) {
    std::array<double, 2> arrFill;
    arrFill[0] = indexX + 0.5;
    arrFill[1] = indexY + 0.5;
    hSparseResponse->Fill(arrFill.data(), val);
  } else {
    std::array<int, 2> arrFill;
    arrFill[0] = indexX;
    arrFill[1] = indexY;
    double currentBinCont = h2d->GetBinContent(arrFill[0], arrFill[1]);
    double currentBinErr = h2d->GetBinError(arrFill[0], arrFill[1]);
    h2d->SetBinContent(arrFill[0], arrFill[1], val + currentBinCont);
    h2d->SetBinError(arrFill[0], arrFill[1], std::sqrt(err * err + currentBinErr * currentBinErr));
  }
}

//____________________________________________________________________________________________________________________________
TH2F* MatrixHandlerNDim::GetTH2(const char* name)
{
  if (useTHNSparese) {
    if (!hSparseResponse) {
      printf("Attention! hSparseResponse does not exist yet!\n");
      return nullptr;
    }
    std::vector<double> vecXBins = MakeAxis1D(vecBinsX);
    std::vector<double> vecYBins = MakeAxis1D(vecBinsY);
    if (h2d) {
      delete h2d;
    }
    h2d = new TH2F(name, name, vecXBins.size() - 1, vecXBins.data(), vecYBins.size() - 1, vecYBins.data());

    // transfer the values
    for (int x = 0; x < h2d->GetXaxis()->GetNbins(); ++x) {
      for (int y = 0; y < h2d->GetYaxis()->GetNbins(); ++y) {
        std::array<int, 2> arrBins = {x + 1, y + 1};
        double binCont = hSparseResponse->GetBinContent(arrBins.data());
        double binErr = hSparseResponse->GetBinError(arrBins.data());
        // printf("x %d, y %d = %f\n", x, y, binCont);
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
TH2F* MatrixHandlerNDim::GetResponseMatrix(std::vector<double> binsX, std::vector<double> binsY, const char* name){
  
  TH2F* hResponse = new TH2F(name, name, vecBinsX.back().size()-1, vecBinsX.back().data(), vecBinsY.back().size()-1, vecBinsY.back().data());

  binsX.push_back(0.5*(vecBinsX.back()[0] + vecBinsX.back()[1]));
  binsY.push_back(0.5*(vecBinsY.back()[0] + vecBinsY.back()[1]));
  int indexX = getIndex(binsX, true);
  int indexY = getIndex(binsY, false);

  for(unsigned x = 0; x < vecBinsX.back().size()-1; ++x){
    for(unsigned y = 0; y < vecBinsY.back().size()-1; ++y){
      if(h2d){
        hResponse->SetBinContent(x+1, y+1, h2d->GetBinContent(indexX + static_cast<int>(x) + 1, indexY + static_cast<int>(y) + 1));
      } else {
        int tmp[2] = {static_cast<int>(indexX + x + 1), static_cast<int>(indexY + y + 1)};
        hResponse->SetBinContent(x+1, y+1, hSparseResponse->GetBinContent(tmp));
      }
    }
  }
  return hResponse;
}

// if index -1 is passed, this is considered to be the one taken for the axis
//____________________________________________________________________________________________________________________________
TH2F* MatrixHandlerNDim::GetResponseMatrix(std::vector<int> binsX, std::vector<int> binsY, const char* name){
  
  if(vecBinsX.size() != binsX.size() || vecBinsY.size() != binsY.size()){
    printf("MatrixHandlerNDim::GetResponseMatrix  The given indeces cannot be correct \n");
    return nullptr;
  }

  std::array<std::vector<double>, 2> arrBinning;
  for(unsigned int i = 0; i < binsX.size(); ++i){
    if(binsX[i] < 0){
      if(arrBinning[0].size() == 0) arrBinning[0] = vecBinsX[i];
      else arrBinning[1] = vecBinsX[i];
    }
  }
  for(unsigned int i = 0; i < binsY.size(); ++i){
    if(binsY[i] < 0){
      if(arrBinning[0].size() == 0) arrBinning[0] = vecBinsY[i];
      else arrBinning[1] = vecBinsY[i];
    }
  }

  const unsigned int dimX = vecBinsX.size();
  std::vector<unsigned int> vecDimX(dimX, 0);
  const unsigned int dimY = vecBinsY.size();
  std::vector<unsigned int> vecDimY(dimY, 0);

  unsigned int nBinsX = 1;
  for(const auto & v : vecBinsX){
      nBinsX*=(v.size()-1);
  }

  unsigned int nBinsY = 1;
  for(const auto & v : vecBinsY){
      nBinsY*=(v.size()-1);
  }
  printf("nBinsY %d\n", nBinsY );

  printf("name %s\n", name);
  TH2F* hResponse = new TH2F(name, name, arrBinning[0].size()-1, arrBinning[0].data(), arrBinning[1].size()-1, arrBinning[1].data());
  
  TH2F* hResponse_org = GetTH2();

  for(unsigned x = 0; x < nBinsX; ++x){

    std::vector<double> valX = getValueForBinIndex(x, true);

    bool acceptBinX = true;
    for(unsigned int i = 0; i < valX.size(); ++i){
      if(binsX[i] < 0){
        continue;
      } else if(valX[i] < vecBinsX[i][binsX[i]] || valX[i] > vecBinsX[i][binsX[i] + 1]){ // which bin??
        acceptBinX = false;
      }
    }

    for(unsigned y = 0; y < nBinsY; ++y){

      std::vector<double> valY = getValueForBinIndex(y, false);

      bool acceptBinY = true;
      for(unsigned int i = 0; i < valY.size(); ++i){
        // if(i==0){
        //   printf(" %f < valY[i] %f < %f\n", vecBinsY[i][binsY[i]], valY[i], vecBinsY[i][binsY[i]+1]);
        // }
        if(binsY[i] < 0){
          continue;
        } else if(valY[i] < vecBinsY[i][binsY[i]] || valY[i] > vecBinsY[i][binsY[i] + 1]){
          acceptBinY = false;
        }
      }

      // if(acceptBinY) printf("accept y bin\n");

      if(acceptBinX && acceptBinY){
        printf("all accepted\n");
        std::array<double, 2> val = {-1e12, -1e12};
        for(unsigned int i = 0; i < binsX.size(); ++i){
          if(binsX[i] < 0){
            if(val[0]  == -1e12){
              val[0] = valX[i];
            } else {
              val[1] = valX[i];
            }
          }
        }

        for(unsigned int i = 0; i < binsY.size(); ++i){
          if(binsY[i] < 0){
            if(val[0]  == -1e12){
              val[0] = valY[i];
            } else {
              val[1] = valY[i];
            }
          }
        }
        printf("val %f,   %f --> %f\n", val[0], val[1], hResponse_org->GetBinContent(x, y));   
        hResponse->SetBinContent(hResponse->GetXaxis()->FindBin(val[0]), hResponse->GetYaxis()->FindBin(val[1]), hResponse_org->GetBinContent(x+1, y+1));
        hResponse->SetBinError(hResponse->GetXaxis()->FindBin(val[0]), hResponse->GetYaxis()->FindBin(val[1]), hResponse_org->GetBinError(x+1, y+1));  
      
      }
      
    }
  }

  return hResponse;
}


//____________________________________________________________________________________________________________________________
THnSparseF* MatrixHandlerNDim::GetTHnSparseClone(const char* name)
{
  return (THnSparseF*)hSparseResponse->Clone(name);
}

//____________________________________________________________________________________________________________________________
THnSparseF* MatrixHandlerNDim::GetTHnSparse(const char* name)
{
  hSparseResponse->SetName(name);
  return hSparseResponse;
}

