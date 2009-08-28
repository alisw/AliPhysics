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
//    it is the PtPt class
//    calculates PtPt correlations for abs and rel var
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------

/* $Id$ */

//-------------------------------------------------------------------------
//         LRC library for Long-Range Correlation analysis
//
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------

#include "AliLRCPtPt.h"

ClassImp(AliLRCPtPt)

/******************************************************
 * AliLRCPtPt class
 ******************************************************/
AliLRCPtPt::AliLRCPtPt():AliLRCAnalysis() {
//Empty constructor
}

AliLRCPtPt::~AliLRCPtPt() {
//Destructor
}

AliLRCPtPt::AliLRCPtPt(char *name, TH2D* sourceHist, double ptd, TH2D* nb):AliLRCAnalysis() {
//Make PtPt form 2d histogram
    SetGraphics();
    CreateHist(name, (char*)"PtPt_abs", (char*)"PtPt_rel", (char*)"Pt_{F}", (char*)"<Pt_{B}>_{Pt_{F}}", (char*)"#frac{Pt_{F}}{<Pt_{F}>}", (char*)"#frac{<Pt_{B}>_{Pt_{F}}}{<Pt_{B}>}", sourceHist);
    SetErrors(sourceHist, name, ptd, nb);
}

AliLRCPtPt::AliLRCPtPt(char *name, TH2D* sourceHist, double ptd, TProfile* nb):AliLRCAnalysis() {
//Make PtPt form 2d histogram
    SetGraphics();
    CreateHist(name, (char*)"PtPt_abs", (char*)"PtPt_rel", (char*)"Pt_{F}", (char*)"<Pt_{B}>_{Pt_{F}}", (char*)"#frac{Pt_{F}}{<Pt_{F}>}", (char*)"#frac{<Pt_{B}>_{Pt_{F}}}{<Pt_{B}>}", sourceHist);
    SetErrors(sourceHist, name, ptd, nb);
}

AliLRCPtPt::AliLRCPtPt(char *fileHistname, char *histname, char *profname, double ptd, char *errhistname):AliLRCAnalysis() {
//Make PtPt form 2d histogram from root file
    SetGraphics();
    fileHist = new TFile(fileHistname);
    TH2D* sourceHist = (TH2D*) fileHist->Get(histname);
    CreateHist(profname, (char*)"PtPt_abs", (char*)"PtPt_rel", (char*)"Pt_{F}", (char*)"<Pt_{B}>_{Pt_{F}}", (char*)"#frac{Pt_{F}}{<Pt_{F}>}", (char*)"#frac{<Pt_{B}>_{Pt_{F}}}{<Pt_{B}>}", sourceHist);
	TProfile* nbP = (TProfile*) fileHist->Get(errhistname);
	SetErrors(sourceHist, profname, ptd, nbP);
}

AliLRCPtPt::AliLRCPtPt(TList *LHist, char *histname, char *profname, char *ptdname, char *errhistname):AliLRCAnalysis() {
//Make PtPt form 2d histogram from root file
    SetGraphics();
    TH2D* sourceHist = (TH2D*) LHist->FindObject(histname);
    CreateHist(profname, (char*)"PtPt_abs", (char*)"PtPt_rel", (char*)"Pt_{F}", (char*)"<Pt_{B}>_{Pt_{F}}", (char*)"#frac{Pt_{F}}{<Pt_{F}>}", (char*)"#frac{<Pt_{B}>_{Pt_{F}}}{<Pt_{B}>}", sourceHist);
	TProfile* nbP = (TProfile*) LHist->FindObject(errhistname);
	TProfile *dPtB = (TProfile*) LHist->FindObject(ptdname);
	double dptb=dPtB->GetBinError(1)*sqrt(dPtB->GetBinEntries(1));
	SetErrors(sourceHist, profname, dptb, nbP);
}

void AliLRCPtPt::MakeHistogramm(char *name, TH2D* sourceHist, double ptd, TH2D* nb) {
//Make PtPt form 2d histogram
    SetGraphics();
    CreateHist(name, (char*)"PtPt_abs", (char*)"PtPt_rel", (char*)"Pt_{F}", (char*)"<Pt_{B}>_{Pt_{F}}", (char*)"#frac{Pt_{F}}{<Pt_{F}>}", (char*)"#frac{<Pt_{B}>_{Pt_{F}}}{<Pt_{B}>}", sourceHist);
    SetErrors(sourceHist, name, ptd, nb);
}

void AliLRCPtPt::MakeHistogramm(char *name, TH2D* sourceHist, double ptd, TProfile* nb) {
//Make PtPt form 2d histogram
    SetGraphics();
    CreateHist(name, (char*)"PtPt_abs", (char*)"PtPt_rel", (char*)"Pt_{F}", (char*)"<Pt_{B}>_{Pt_{F}}", (char*)"#frac{Pt_{F}}{<Pt_{F}>}", (char*)"#frac{<Pt_{B}>_{Pt_{F}}}{<Pt_{B}>}", sourceHist);
    SetErrors(sourceHist, name, ptd, nb);
}


void AliLRCPtPt::MakeHistogramm(char *fileHistname, char *histname, char *profname, double ptd, char *errhistname) {
//Make PtPt form 2d histogram from root file
    SetGraphics();
    fileHist = new TFile(fileHistname);
    TH2D* sourceHist = (TH2D*) fileHist->Get(histname);
    CreateHist(profname, (char*)"PtPt_abs", (char*)"PtPt_rel", (char*)"Pt_{F}", (char*)"<Pt_{B}>_{Pt_{F}}", (char*)"#frac{Pt_{F}}{<Pt_{F}>}", (char*)"#frac{<Pt_{B}>_{Pt_{F}}}{<Pt_{B}>}", sourceHist);
	TProfile* nbP = (TProfile*) fileHist->Get(errhistname);
	SetErrors(sourceHist, profname, ptd, nbP);
}

void AliLRCPtPt::MakeHistogramm(TList *LHist, char *histname, char *profname, char *ptdname, char *errhistname) {
//Make PtPt form 2d histogram from root file
    SetGraphics();
    TH2D* sourceHist = (TH2D*) LHist->FindObject(histname);
    CreateHist(profname, (char*)"PtPt_abs", (char*)"PtPt_rel", (char*)"Pt_{F}", (char*)"<Pt_{B}>_{Pt_{F}}", (char*)"#frac{Pt_{F}}{<Pt_{F}>}", (char*)"#frac{<Pt_{B}>_{Pt_{F}}}{<Pt_{B}>}", sourceHist);
	TProfile* nbP = (TProfile*) LHist->FindObject(errhistname);
	TProfile *dPtB = (TProfile*) LHist->FindObject(ptdname);
	double dptb=dPtB->GetBinError(1)*sqrt(dPtB->GetBinEntries(1));
	SetErrors(sourceHist, profname, dptb, nbP);
}
