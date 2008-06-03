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
//    it is the PtN class
//    calculates PtN correlations for abs and rel var
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------

/* $Id$ */

//-------------------------------------------------------------------------
//         LRC library for Long-Range Correlation analysis
//
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------

#include "AliLRCPtN.h"

ClassImp(AliLRCPtN) 

/******************************************************
 * AliLRCPtN class
 ******************************************************/



AliLRCPtN::AliLRCPtN():AliLRCAnalysis() {
//Empty constructor
}

AliLRCPtN::~AliLRCPtN() {
//Destructor
}

AliLRCPtN::AliLRCPtN(char *name, TH2D* sourceHist, double ptd, TH2D* nb):AliLRCAnalysis() {
//Make PtN form 2d histogram
    SetGraphics();
    CreateHist(name, (char*)"PtN_abs", (char*)"PtN_rel", (char*)"n_{F}", (char*)"Pt_{B}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<Pt_{B}>_{n_{F}}}{<Pt_{B}>}", sourceHist);
    SetErrors(sourceHist, name, ptd, nb);
}

AliLRCPtN::AliLRCPtN(char *fileHistname, char *histname, char *profname, double ptd, char *errhistname):AliLRCAnalysis() {
//Make PtN form 2d histogram from root file
    SetGraphics();
    fileHist = new TFile(fileHistname);
    TH2D* sourceHist = (TH2D*) fileHist->Get(histname);
    TH2D* nb = (TH2D*) fileHist->Get(errhistname);
    CreateHist(profname, (char*)"PtN_abs", (char*)"PtN_rel", (char*)"n_{F}", (char*)"Pt_{B}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<Pt_{B}>_{n_{F}}}{<Pt_{B}>}", sourceHist);
    SetErrors(sourceHist, profname, ptd, nb);
    
}

void AliLRCPtN::MakeHistogramm(char *name, TH2D* sourceHist, double ptd, TH2D* nb) {
//Make PtN form 2d histogram
    SetGraphics();
    CreateHist(name, (char*)"PtN_abs", (char*)"PtN_rel", (char*)"n_{F}", (char*)"Pt_{B}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<Pt_{B}>_{n_{F}}}{<Pt_{B}>}", sourceHist);
    SetErrors(sourceHist, name, ptd, nb);
}

void AliLRCPtN::MakeHistogramm(char *fileHistname, char *histname, char *profname, double ptd, char *errhistname) {
//Make PtN form 2d histogram from root file
    SetGraphics();
    fileHist = new TFile(fileHistname);
    TH2D* sourceHist = (TH2D*) fileHist->Get(histname);
    TH2D* nb = (TH2D*) fileHist->Get(errhistname);
    CreateHist(profname, (char*)"PtN_abs", (char*)"PtN_rel", (char*)"n_{F}", (char*)"Pt_{B}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<Pt_{B}>_{n_{F}}}{<Pt_{B}>}", sourceHist);
    SetErrors(sourceHist, profname, ptd, nb);
}
