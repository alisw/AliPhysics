
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
MatrixHandler4D::MatrixHandler4D(std::vector<double> arrMesonX, std::vector<double> arrMesonY, std::vector<double> arrJetX, std::vector<double> arrJetY, bool useTHN){
    vecBinsMesonX = arrMesonX;
    vecBinsMesonY = arrMesonY;
    vecBinsJetX = arrJetX;
    vecBinsJetY = arrJetY;
    useTHNSparese = useTHN;
    const int nBinsX = (arrMesonX.size() - 1) * (arrJetX.size() - 1); // + 1;
    const int nBinsY = (arrMesonY.size() - 1) * (arrJetY.size() - 1); // + 1;
    if(useTHNSparese){
        // in case of thnsparse, just use equidistant binning
        std::array<int, 2> arrNBins = {nBinsX, nBinsY};
        std::array<double, 2> arrXBins = {0, 0};
        std::array<double, 2> arrYBins = {static_cast<double>(nBinsX), static_cast<double>(nBinsY)}; 

        hSparseResponse = new THnSparseF("ResponseMatrix_dyn", "ResponseMatrix_dyn", arrNBins.size(), arrNBins.data(), arrXBins.data(), arrYBins.data());

    } else {
        std::vector<double> vecXBins;
        std::vector<double> vecYBins;
        GetAxisBinning(vecXBins, vecYBins);
        if(h2d){
                delete h2d;
            }
        h2d = new TH2F("ResponseMatrix_stat", "ResponseMatrix_stat", vecXBins.size() - 1, vecXBins.data(), vecYBins.size() - 1, vecYBins.data() );
        
    }
}

//____________________________________________________________________________________________________________________________
MatrixHandler4D::MatrixHandler4D(std::vector<double> arrMesonX, std::vector<double> arrMesonY, std::vector<double> arrJetX, std::vector<double> arrJetY, THnSparse*  h){
    vecBinsMesonX = arrMesonX;
    vecBinsMesonY = arrMesonY;
    vecBinsJetX = arrJetX;
    vecBinsJetY = arrJetY;
    useTHNSparese = true;
    const int nBinsX = (arrMesonX.size() - 1) * (arrJetX.size() - 1); // + 1;
    const int nBinsY = (arrMesonY.size() - 1) * (arrJetY.size() - 1); // + 1;
    hSparseResponse = (THnSparseF*) h->Clone();
}

//____________________________________________________________________________________________________________________________
MatrixHandler4D::~MatrixHandler4D(){
//  if(h2d){
//      delete h2d;
//  }
//  if(hSparseResponse){
//      delete hSparseResponse;
//  }
}

//____________________________________________________________________________________________________________________________
int MatrixHandler4D::getBinIndexMesonX(const double val) const{
    for(unsigned int i = 0; i < vecBinsMesonX.size() - 1; ++i){
        if(vecBinsMesonX[i] < val && vecBinsMesonX[i+1] > val){
            return i;
        }
    }
    return -1;
}

//____________________________________________________________________________________________________________________________
int MatrixHandler4D::getBinIndexMesonY(const double val) const{
    for(unsigned int i = 0; i < vecBinsMesonY.size() - 1; ++i){
        if(vecBinsMesonY[i] < val && vecBinsMesonY[i+1] > val){
            return i;
        }
    }
    return -1;
}

//____________________________________________________________________________________________________________________________
int MatrixHandler4D::getBinIndexJetX(const double val) const{
    for(unsigned int i = 0; i < vecBinsJetX.size() - 1; ++i){
        if(vecBinsJetX[i] < val && vecBinsJetX[i+1] > val){
            return i;
        }
    }
    return -1;
}

//____________________________________________________________________________________________________________________________
int MatrixHandler4D::getBinIndexJetY(const double val) const{
    for(unsigned int i = 0; i < vecBinsJetY.size() - 1; ++i){
        if(vecBinsJetY[i] < val && vecBinsJetY[i+1] > val){
            return i;
        }
    }
    return -1;
}