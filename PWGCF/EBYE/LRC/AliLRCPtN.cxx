/**************************************************************************
 * Author: Andrey Ivanov.                                                 *
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
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch,
//    Andrey Ivanov (SPbSU-CERN), Igor Altsebeev (SPbSU-CERN) 
//-------------------------------------------------------------------------

/* $Id$ */


#include "AliLRCPtN.h"
#include "TFile.h"
#include "TList.h"
#include "TProfile.h"
#include "math.h"

class TFile;
class TProfile;
class TList;
class TH2D;

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
    CreateHist(name, (char*)"PtN_abs", (char*)"PtN_rel", (char*)"n_{F}", (char*)"<Pt_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<Pt_{B}>_{n_{F}}}{<Pt_{B}>}", sourceHist);
    SetErrors(sourceHist, name, ptd, nb);
}

AliLRCPtN::AliLRCPtN(char *fileHistname, char *histname, char *profname, double ptd, char *errhistname):AliLRCAnalysis() {
//Make PtN form 2d histogram from root file
    SetGraphics();
    fileHist = new TFile(fileHistname);
    TH2D* sourceHist = (TH2D*) fileHist->Get(histname);
    
    CreateHist(profname, (char*)"PtN_abs", (char*)"PtN_rel", (char*)"n_{F}", (char*)"<Pt_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<Pt_{B}>_{n_{F}}}{<Pt_{B}>}", sourceHist);
	TProfile* nbP = (TProfile*) fileHist->Get(errhistname);
	SetErrors(sourceHist, profname, ptd, nbP);
	
}
//const
AliLRCPtN::AliLRCPtN(const TList *  const LHist, char *histname, char *profname, char *ptdname, char *errhistname):AliLRCAnalysis() {
//Make PtN form 2d histogram from root file
    SetGraphics();
    TH2D* sourceHist = (TH2D*) LHist->FindObject(histname);
    CreateHist(profname, (char*)"PtN_abs", (char*)"PtN_rel", (char*)"n_{F}", (char*)"<Pt_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<Pt_{B}>_{n_{F}}}{<Pt_{B}>}", sourceHist);
	TProfile* nbP = (TProfile*) LHist->FindObject(errhistname);
	TProfile *dPtB = (TProfile*) LHist->FindObject(ptdname);
	double dptb=dPtB->GetBinError(1)*sqrt(dPtB->GetBinEntries(1));
	SetErrors(sourceHist, profname, dptb, nbP);
}

AliLRCPtN::AliLRCPtN(char *name, TH2D* sourceHist, double ptd, TProfile* nb):AliLRCAnalysis() {
//Make PtN form 2d histogram
    SetGraphics();
    CreateHist(name, (char*)"PtN_abs", (char*)"PtN_rel", (char*)"n_{F}", (char*)"<Pt_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<Pt_{B}>_{n_{F}}}{<Pt_{B}>}", sourceHist);
    SetErrors(sourceHist, name, ptd, nb);
}

void AliLRCPtN::MakeHistogramm(char *name, TH2D* sourceHist, double ptd, TH2D* nb) {
//Make PtN form 2d histogram
    SetGraphics();
    CreateHist(name, (char*)"PtN_abs", (char*)"PtN_rel", (char*)"n_{F}", (char*)"<Pt_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<Pt_{B}>_{n_{F}}}{<Pt_{B}>}", sourceHist);
    SetErrors(sourceHist, name, ptd, nb);
}

void AliLRCPtN::MakeHistogramm(char *name, TH2D* sourceHist, double ptd, TProfile* nb) {
//Make PtN form 2d histogram
    SetGraphics();
    CreateHist(name, (char*)"PtN_abs", (char*)"PtN_rel", (char*)"n_{F}", (char*)"<Pt_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<Pt_{B}>_{n_{F}}}{<Pt_{B}>}", sourceHist);
    SetErrors(sourceHist, name, ptd, nb);
}

void AliLRCPtN::MakeHistogramm(char *fileHistname, char *histname, char *profname, double ptd, char *errhistname) {
//Make PtN form 2d histogram from root file
    SetGraphics();
    fileHist = new TFile(fileHistname);
    TH2D* sourceHist = (TH2D*) fileHist->Get(histname);
    CreateHist(profname, (char*)"PtN_abs", (char*)"PtN_rel", (char*)"n_{F}", (char*)"<Pt_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<Pt_{B}>_{n_{F}}}{<Pt_{B}>}", sourceHist);
	TProfile* nbP = (TProfile*) fileHist->Get(errhistname);
	SetErrors(sourceHist, profname, ptd, nbP);
}

void AliLRCPtN::MakeHistogramm(const TList * const LHist, char *histname, char *profname, char *ptdname, char *errhistname) {
//Make PtN form 2d histogram from root file
    SetGraphics();
    TH2D* sourceHist = (TH2D*) LHist->FindObject(histname);
    CreateHist(profname, (char*)"PtN_abs", (char*)"PtN_rel", (char*)"n_{F}", (char*)"<Pt_{B}>_{n_{F}}", (char*)"#frac{n_{F}}{<n_{F}>}", (char*)"#frac{<Pt_{B}>_{n_{F}}}{<Pt_{B}>}", sourceHist);
	TProfile* nbP = (TProfile*) LHist->FindObject(errhistname);
	TProfile *dPtB = (TProfile*) LHist->FindObject(ptdname);
	double dptb=dPtB->GetBinError(1)*sqrt(dPtB->GetBinEntries(1));
	SetErrors(sourceHist, profname, dptb, nbP);
}
