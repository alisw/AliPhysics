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
#include "TF1.h"
#include "TF2.h"
#include "THnSparse.h"

class MatrixHandler4D
{

 public:
  MatrixHandler4D() = default;
  MatrixHandler4D(std::vector<double> arrMesonX, std::vector<double> arrMesonY, std::vector<double> arrJetX, std::vector<double> arrJetY, bool useTHN = false);
  MatrixHandler4D(std::vector<double> arrMesonX, std::vector<double> arrMesonY, std::vector<double> arrJetX, std::vector<double> arrJetY, THnSparse* h = nullptr);
  MatrixHandler4D(std::vector<double> arrMesonX, std::vector<double> arrMesonY, std::vector<double> arrJetX, std::vector<double> arrJetY, TH2F* h = nullptr);
  MatrixHandler4D(const MatrixHandler4D&) = delete;            // copy ctor
  MatrixHandler4D(MatrixHandler4D&&) = delete;                 // move ctor
  MatrixHandler4D& operator=(const MatrixHandler4D&) = delete; // copy assignment
  MatrixHandler4D& operator=(MatrixHandler4D&&) = delete;      // move assignment
  virtual ~MatrixHandler4D() = default;

  int getBinIndexMesonX(const double val) const;
  int getBinIndexMesonY(const double val) const;
  int getBinIndexJetX(const double val) const;
  int getBinIndexJetY(const double val) const;

  double getValueForBinIndexMesonX(const int index) const;
  double getValueForBinIndexJetX(const int index) const;
  double getValueForBinIndexMesonY(const int index) const;
  double getValueForBinIndexJetY(const int index) const;


  std::vector<double> getBinsMesonX() const { return vecBinsMesonX; }
  std::vector<double> getBinsMesonY() const { return vecBinsMesonY; }
  std::vector<double> getBinsJetX() const { return vecBinsJetX; }
  std::vector<double> getBinsJetY() const { return vecBinsJetY; }

  void Fill(double valJetX, double valJetY, double valMesonX, double valMesonY, double val = 1);

  void AddBinContent(double valJetX, double valJetY, double valMesonX, double valMesonY, double val = 1, double err = 1);
  void GetAxisBinning(std::vector<double>& vecXBins, std::vector<double>& vecYBins);

  THnSparseF* GetTHnSparseClone(const char* name = "hSparseResponse_Clone");
  THnSparseF* GetTHnSparse(const char* name = "");

  TH2F* GetTH2(const char* name = "hSparseResponse_Clone");
  TH2F* GetResponseMatrix(int binX, int binY, const char* name = "dummy");

  void WeightResponseMatrix(TF1* funcMeson = nullptr, TF1* funcJet = nullptr);
  void WeightResponseMatrix(TF2* func);



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

  ClassDef(MatrixHandler4D, 4)
};



class MatrixHandlerNDim
{

 public:
  MatrixHandlerNDim() = default;
  MatrixHandlerNDim(std::vector<std::vector<double>> arrBinsX, std::vector<std::vector<double>> arrBinsY, bool useTHN = false);
  MatrixHandlerNDim(std::vector<std::vector<double>> arrX, std::vector<std::vector<double>> arrY, THnSparse* h = nullptr);
  MatrixHandlerNDim(std::vector<std::vector<double>> arrX, std::vector<std::vector<double>> arrY, TH2F* h = nullptr);
  MatrixHandlerNDim(const MatrixHandlerNDim&) = delete;            // copy ctor
  MatrixHandlerNDim(MatrixHandlerNDim&&) = delete;                 // move ctor
  MatrixHandlerNDim& operator=(const MatrixHandlerNDim&) = delete; // copy assignment
  MatrixHandlerNDim& operator=(MatrixHandlerNDim&&) = delete;      // move assignment
  virtual ~MatrixHandlerNDim() = default;

  int getBinIndexMesonX(const double val) const;
  int getBinIndexMesonY(const double val) const;
  int getBinIndexJetX(const double val) const;
  int getBinIndexJetY(const double val) const;

  unsigned int getIndex(std::vector<double> vecVal, bool isXAxis);

  // double getValueForBinIndexMesonX(const int index) const;
  // double getValueForBinIndexJetX(const int index) const;
  // double getValueForBinIndexMesonY(const int index) const;
  // double getValueForBinIndexJetY(const int index) const;

  TH2F* GetResponseMatrix(std::vector<int> binsX, std::vector<int> binsY, const char* name);

  std::vector<double> getValueForBinIndex(const unsigned int index, bool isXAxis);

  void Fill(double valJetX, double valJetY, double valMesonX, double valMesonY, double val = 1);
  void Fill(std::vector<double> valRec, std::vector<double> valTrue, double val = 1);

  void AddBinContent(double valJetX, double valJetY, double valMesonX, double valMesonY, double val = 1, double err = 1);
  void AddBinContent(std::vector<double> valRec, std::vector<double> valTrue, double val = 1, double err = 1);

  double GetBinLowEdge(std::vector<unsigned int> vecDim, const std::vector<std::vector<double>> vecND);
  std::vector<double> MakeAxis1D(const std::vector<std::vector<double>> vecND);

  THnSparseF* GetTHnSparseClone(const char* name = "hSparseResponse_Clone");
  THnSparseF* GetTHnSparse(const char* name = "");

  TH2F* GetTH2(const char* name = "hSparseResponse_Clone");
  TH2F* GetResponseMatrix(std::vector<double> binsX, std::vector<double> binsY, const char* name = "dummy");

 private:
  bool useTHNSparese = false;
  int nBinsJet = 0;
  std::vector<std::vector<double>> vecBinsX = {};
  std::vector<std::vector<double>> vecBinsY = {};
  TH2F* h2d = nullptr;
  TH1F* h1dJet = nullptr;
  TH1F* h1dMeson = nullptr;
  THnSparseF* hSparseResponse = nullptr;

  ClassDef(MatrixHandlerNDim, 1)
};

#endif