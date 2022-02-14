#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <fstream>
#include <iostream> 
#include <string> 
#include <vector>

#include "AliExternalBDT.h"

#define DELTA 1.0e-6

int test_AliEsternalBDT(string path = "") {

  string tree_path, model_path;

  if (path == "") {
    tree_path  = "test_tree_pt8_12.root";
    model_path = "test_xgboost_pt8_12.model";
  } else {
    tree_path  = path + "/" + "test_tree_pt8_12.root";
    model_path = path + "/" + "test_xgboost_pt8_12.model";
  }

  fstream fAliExtBDT_Pred, fXGBoost_Pred;

  fAliExtBDT_Pred.open("ali_extbdt_pred.txt");

  TFile *fInput = new TFile(tree_path.data(), "READ");

  TTreeReader fReader("tree_real_data", fInput);

  TTreeReaderValue<float> fValueDeltaMass(fReader, "delta_mass_KK");
  TTreeReaderValue<float> fValueDLen(fReader, "d_len");
  TTreeReaderValue<float> fValueNormDLXY(fReader, "norm_dl_xy");
  TTreeReaderValue<float> fValueSigVert(fReader, "sig_vert");
  TTreeReaderValue<float> fValueCosPiKPhi(fReader, "cos_PiKPhi_3");
  TTreeReaderValue<float> fValueNormIP(fReader, "norm_IP");
  TTreeReaderValue<float> fValueSigCombK0(fReader, "sigComb_K_0");
  TTreeReaderValue<float> fValueSigCombK1(fReader, "sigComb_K_1");
  TTreeReaderValue<float> fValueSigCombK2(fReader, "sigComb_K_2");
  TTreeReaderValue<float> fValueSigCombPi0(fReader, "sigComb_Pi_0");
  TTreeReaderValue<float> fValueSigCombPi1(fReader, "sigComb_Pi_1");
  TTreeReaderValue<float> fValueSigCombPi2(fReader, "sigComb_Pi_2");

  AliExternalBDT *fBDT = new AliExternalBDT();

  if (!fBDT->LoadXGBoostModel(model_path.data())) {
    return 1;
  }

  while (fReader.Next()) {

    double features[12] = {*fValueDeltaMass,  *fValueDLen,       *fValueNormDLXY,
                           *fValueSigVert,    *fValueCosPiKPhi,  *fValueNormIP,
                           *fValueSigCombK0,  *fValueSigCombK1,  *fValueSigCombK2,
                           *fValueSigCombPi0, *fValueSigCombPi1, *fValueSigCombPi2};

    fAliExtBDT_Pred << Form("%.10f", fBDT->Predict(features, 12, true)) << std::endl;
  }
  fInput->Close();
  delete fBDT;

  fAliExtBDT_Pred.clear();
  fAliExtBDT_Pred.seekg(0, ios::beg);

  string aliextbdt_score, xgboost_score;

  fXGBoost_Pred.open("xgboost_pred.txt", ios::in);

  float scoreA, scoreX;

  while (getline(fAliExtBDT_Pred, aliextbdt_score)) {

    getline(fXGBoost_Pred, xgboost_score);

    scoreX = std::stod(xgboost_score);
    scoreA = std::stod(aliextbdt_score);

    if (abs(scoreX - scoreA) > DELTA) {
      fXGBoost_Pred.close();
      fAliExtBDT_Pred.close();

      std::cout << "TEST: Fail!" << std::endl;
      return 1;
    }
  }
  fXGBoost_Pred.close();
  fAliExtBDT_Pred.close();

  std::cout << "TEST: Success!" << std::endl;
  return 0;
}
