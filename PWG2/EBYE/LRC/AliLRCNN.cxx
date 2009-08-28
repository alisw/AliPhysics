/**************************************************************************
 * Author: Panos Christakoglou.                                           *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//-------------------------------------------------------------------------
//    Description: 
//    This class include into LRC library for Long-Range Correlation analysis
//    it is the NN class
//    calculates NN correlations for abs and rel var
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------

/* $Id$ */

#include "AliLRCNN.h"

ClassImp(AliLRCNN) 


/******************************************************
 * AliLRCNN class
 ******************************************************/
AliLRCNN::AliLRCNN():AliLRCAnalysis(){
 //Empty constructor
}


AliLRCNN::~AliLRCNN() {
//Destructor   
}


AliLRCNN::AliLRCNN(char *name, TH2D* sourceHist):AliLRCAnalysis() {
//Make NN from 2d histogramm
    SetGraphics();
    CreateHist(name, (char*)"NN_abs", (char*)"NN_rel", (char*)"n_{F}", (char*)"<n_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<n_{B}>_{n_{F}}}{<n_{B}>}", sourceHist);
    SetErrors(sourceHist, name);
}

AliLRCNN::AliLRCNN(char *fileHistname, char *histname, char *profname):AliLRCAnalysis() {
//Make NN from 2d histogramm from root file
    SetGraphics();
    fileHist = new TFile(fileHistname);
    TH2D* sourceHist = (TH2D*) fileHist->Get(histname);
    CreateHist(profname, (char*)"NN_abs", (char*)"NN_rel", (char*)"n_{F}", (char*)"<n_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<n_{B}>_{n_{F}}}{<n_{B}>}", sourceHist);
    SetErrors(sourceHist, profname);
}

AliLRCNN::AliLRCNN(TList *LHist, char *histname, char *profname):AliLRCAnalysis() {
//Make NN from 2d histogramm from root file
    SetGraphics();
    TH2D* sourceHist = (TH2D*) LHist->FindObject(histname);
    CreateHist(profname, (char*)"NN_abs", (char*)"NN_rel", (char*)"n_{F}", (char*)"<n_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<n_{B}>_{n_{F}}}{<n_{B}>}", sourceHist);
    SetErrors(sourceHist, profname);
}


void AliLRCNN::MakeHistogramm(char *name, TH2D* sourceHist) {
//Make NN from 2d histogramm
    SetGraphics();
    CreateHist(name, (char*)"NN_abs", (char*)"NN_rel", (char*)"n_{F}", (char*)"<n_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<n_{B}>_{n_{F}}}{<n_{B}>}", sourceHist);
    SetErrors(sourceHist, name);
}

void AliLRCNN::MakeHistogramm(char *fileHistname, char *histname, char *profname) {
//Make NN from 2d histogramm from root file
    SetGraphics();
    fileHist = new TFile(fileHistname);
    TH2D* sourceHist = (TH2D*) fileHist->Get(histname);
    CreateHist(profname, (char*)"NN_abs", (char*)"NN_rel", (char*)"n_{F}", (char*)"<n_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<n_{B}>_{n_{F}}}{<n_{B}>}", sourceHist);
    SetErrors(sourceHist, profname);
}

void AliLRCNN::MakeHistogramm(TList *LHist, char *histname, char *profname) {
//Make NN from 2d histogramm from root file
    SetGraphics();
    TH2D* sourceHist = (TH2D*) LHist->FindObject(histname);
    CreateHist(profname, (char*)"NN_abs", (char*)"NN_rel", (char*)"n_{F}", (char*)"<n_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<n_{B}>_{n_{F}}}{<n_{B}>}", sourceHist);
    SetErrors(sourceHist, profname);
}


