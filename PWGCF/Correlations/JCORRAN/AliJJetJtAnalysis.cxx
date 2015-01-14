/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

// Comment describing what this class does needed!

//===========================================================
// Dummy comment, should be replaced by a real one
// comment
// comment
// Simple class for the jt anlyais by Beomkyu Kim and Dongjo Kim
//===========================================================

#include <TRandom.h>
#include <TMath.h>
#include <TRegexp.h>
#include <TVector.h>
#include "AliJJet.h"
#include "AliJEfficiency.h"
#include "AliJJetJtAnalysis.h"
#include "AliJHistManager.h"
#include "TClonesArray.h"

AliJJetJtAnalysis::AliJJetJtAnalysis():
	fInputList(NULL)
	, fJetList(NULL)
	, fJetListOfList()
	//, fJetBgList(NULL)
	, fJetBgListOfList()
    , fJetTriggPtBorders(NULL)
    , fJetConstPtLowLimits(NULL)
    , fJetAssocPtBorders(NULL)
    , fDeltaRBorders(NULL)
    , nJetContainer(0)
	, fCard(NULL)
	, fJJetAnalysis(NULL)
    , fJetFinderName(0)
    , fConeSizes(0)
	, fEfficiency(0x0)
	, cBin(-1)
	, fcent(-999)
	, zBin(-1)
	, zVert(-999)
    , fTracks(NULL)
    , fHMG(NULL)
    , fJetFinderBin()
    , fJetTriggerBin()
    , fTrkPtBin()
    , fTrkLimPtBin()
    , fdRBin()
    , fiHist()
    , fhNumber()
    , fhKNumber()
    , fhJetPt()
    , fhJetPtBin()
    , fhZ()
    , fhZBin()
    , fhJt()
    , fhJtBin()
    , fhJtWeightBin()
    , fhLogJtWeightBin()
    , fhJtWithPtCutWeightBinBin()
    , fhLogJtWithPtCutWeightBinBin()
    , fhJtBinLimBin()
    , fhJtWeightBinLimBin()
    , fhLogJtWeightBinLimBin()
    , fhJetBgPt()
    , fhJetBgPtBin()
    , fhBgZ()
    , fhBgZBin()
    , fhBgJt()
    , fhBgJtBin()
    , fhBgJtWeightBin()
    , fhBgLogJtWeightBin()
    , fhBgJtWithPtCutWeightBinBin()
    , fhBgLogJtWithPtCutWeightBinBin()
    , fhBgJtWithPtCutWeightBinBinSmallerR()
    , fhBgLogJtWithPtCutWeightBinBinSmallerR()
    , fhBgJtWithPtCutWeightBinBinDiffR()
    , fhBgLogJtWithPtCutWeightBinBinDiffR()
    , fhBgJtBinLimBin()
    , fhBgJtWeightBinLimBin()
    , fhBgLogJtWeightBinLimBin()
    , fhdeltaE()
    , fhdeltaN()
    , fhFullJetEChJetBin()
    , fhFullChdRChJetBin()
    , fh2DFullEvsChEdN0()
    , fh2DFullEvsChEdNnot0()
{

}

AliJJetJtAnalysis::AliJJetJtAnalysis( AliJCard * card ):
	fInputList(NULL)
	, fJetList(NULL)
	, fJetListOfList()
	//, fJetBgList(NULL)
	, fJetBgListOfList()
    , fJetTriggPtBorders(NULL)
    , fJetConstPtLowLimits(NULL)
    , fJetAssocPtBorders(NULL)
    , fDeltaRBorders(NULL)
    , nJetContainer(0)
	, fCard(card)
	, fJJetAnalysis(NULL)
    , fJetFinderName(0)
    , fConeSizes(0)
	, fEfficiency(0x0)
	, cBin(-1)
	, fcent(-999)
	, zBin(-1)
	, zVert(-999)
    , fTracks(NULL)
    , fHMG(NULL)
    , fJetFinderBin()
    , fJetTriggerBin()
    , fTrkPtBin()
    , fTrkLimPtBin()
    , fdRBin()
    , fiHist()
    , fhNumber()
    , fhKNumber()
    , fhJetPt()
    , fhJetPtBin()
    , fhZ()
    , fhZBin()
    , fhJt()
    , fhJtBin()
    , fhJtWeightBin()
    , fhLogJtWeightBin()
    , fhJtWithPtCutWeightBinBin()
    , fhLogJtWithPtCutWeightBinBin()
    , fhJtBinLimBin()
    , fhJtWeightBinLimBin()
    , fhLogJtWeightBinLimBin()
    , fhJetBgPt()
    , fhJetBgPtBin()
    , fhBgZ()
    , fhBgZBin()
    , fhBgJt()
    , fhBgJtBin()
    , fhBgJtWeightBin()
    , fhBgLogJtWeightBin()
    , fhBgJtWithPtCutWeightBinBin()
    , fhBgLogJtWithPtCutWeightBinBin()
    , fhBgJtWithPtCutWeightBinBinSmallerR()
    , fhBgLogJtWithPtCutWeightBinBinSmallerR()
    , fhBgJtWithPtCutWeightBinBinDiffR()
    , fhBgLogJtWithPtCutWeightBinBinDiffR()
    , fhBgJtBinLimBin()
    , fhBgJtWeightBinLimBin()
    , fhBgLogJtWeightBinLimBin()
    , fhdeltaE()
    , fhdeltaN()
    , fhFullJetEChJetBin()
    , fhFullChdRChJetBin()
    , fh2DFullEvsChEdN0()
    , fh2DFullEvsChEdNnot0()
{

}

AliJJetJtAnalysis::AliJJetJtAnalysis(const AliJJetJtAnalysis& ap) :
	fInputList(ap.fInputList)
	, fJetList(ap.fJetList)
	, fJetListOfList(ap.fJetListOfList)
	//, fJetBgList(ap.fJetBgList)
	, fJetBgListOfList(ap.fJetBgListOfList)
    , fJetTriggPtBorders(ap.fJetTriggPtBorders)
    , fJetConstPtLowLimits(ap.fJetConstPtLowLimits)
    , fJetAssocPtBorders(ap.fJetAssocPtBorders)
    , fDeltaRBorders(ap.fDeltaRBorders)
    , nJetContainer(ap.nJetContainer)
	, fCard(ap.fCard)
	, fJJetAnalysis(ap.fJJetAnalysis)
    , fJetFinderName(ap.fJetFinderName)
    , fConeSizes(ap.fConeSizes)
	, fEfficiency(ap.fEfficiency)
	, cBin(-1)
	, fcent(-999)
	, zBin(-1)
	, zVert(-999)
    , fTracks(ap.fTracks)
    , fHMG(ap.fHMG)
    , fJetFinderBin(ap.fJetFinderBin)
    , fJetTriggerBin(ap.fJetTriggerBin)
    , fTrkPtBin(ap.fTrkPtBin)
    , fTrkLimPtBin(ap.fTrkLimPtBin)
    , fdRBin(ap.fdRBin)
    , fiHist(ap.fiHist)
    , fhNumber(ap.fhNumber)
    , fhKNumber(ap.fhKNumber)
    , fhJetPt(ap.fhJetPt)
    , fhJetPtBin(ap.fhJetPtBin)
    , fhZ(ap.fhZ)
    , fhZBin(ap.fhZBin)
    , fhJt(ap.fhJt)
    , fhJtBin(ap.fhJtBin)
    , fhJtWeightBin(ap.fhJtWeightBin)
    , fhLogJtWeightBin(ap.fhLogJtWeightBin)
    , fhJtWithPtCutWeightBinBin(ap.fhJtWithPtCutWeightBinBin)
    , fhLogJtWithPtCutWeightBinBin(ap.fhLogJtWithPtCutWeightBinBin)
    , fhJtBinLimBin(ap.fhJtBinLimBin)
    , fhJtWeightBinLimBin(ap.fhJtWeightBinLimBin)
    , fhLogJtWeightBinLimBin(ap.fhLogJtWeightBinLimBin)
    , fhJetBgPt(ap.fhJetBgPt)
    , fhJetBgPtBin(ap.fhJetBgPtBin)
    , fhBgZ(ap.fhBgZ)
    , fhBgZBin(ap.fhBgZBin)
    , fhBgJt(ap.fhBgJt)
    , fhBgJtBin(ap.fhBgJtBin)
    , fhBgJtWeightBin(ap.fhBgJtWeightBin)
    , fhBgLogJtWeightBin(ap.fhBgLogJtWeightBin)
    , fhBgJtWithPtCutWeightBinBin(ap.fhBgJtWithPtCutWeightBinBin)
    , fhBgLogJtWithPtCutWeightBinBin(ap.fhBgLogJtWithPtCutWeightBinBin)
    , fhBgJtWithPtCutWeightBinBinSmallerR(ap.fhBgJtWithPtCutWeightBinBinSmallerR)
    , fhBgLogJtWithPtCutWeightBinBinSmallerR(ap.fhBgLogJtWithPtCutWeightBinBinSmallerR)
    , fhBgJtWithPtCutWeightBinBinDiffR(ap.fhBgJtWithPtCutWeightBinBinDiffR)
    , fhBgLogJtWithPtCutWeightBinBinDiffR(ap.fhBgLogJtWithPtCutWeightBinBinDiffR)
    , fhBgJtBinLimBin(ap.fhBgJtBinLimBin)
    , fhBgJtWeightBinLimBin(ap.fhBgJtWeightBinLimBin)
    , fhBgLogJtWeightBinLimBin(ap.fhBgLogJtWeightBinLimBin)
    , fhdeltaE(ap.fhdeltaE)
    , fhdeltaN(ap.fhdeltaN)
    , fhFullJetEChJetBin(ap.fhFullJetEChJetBin)
    , fhFullChdRChJetBin(ap.fhFullChdRChJetBin)
    , fh2DFullEvsChEdN0(ap.fh2DFullEvsChEdN0)
    , fh2DFullEvsChEdNnot0(ap.fh2DFullEvsChEdNnot0)
{

}

AliJJetJtAnalysis& AliJJetJtAnalysis::operator = (const AliJJetJtAnalysis& ap)
{
	// assignment operator

	this->~AliJJetJtAnalysis();
	new(this) AliJJetJtAnalysis(ap);
	return *this;
}


AliJJetJtAnalysis::~AliJJetJtAnalysis(){


    delete fJJetAnalysis;
    fJetFinderName.clear();
    fConeSizes.clear();
    delete fEfficiency;
    delete fHMG;
   


}



void AliJJetJtAnalysis::UserCreateOutputObjects(){
	//fJetListOfList always point one address in the whole time of this analysis.
    //Thus mustn't be cleared in it's life.     
    //fJetListOfList.Clear();


    fJJetAnalysis = new AliJJetAnalysis();

    fJetTriggPtBorders = fCard->GetVector("JetTriggPtBorders");
    fJetConstPtLowLimits = fCard->GetVector("JetConstPtLowLimits");
    fJetAssocPtBorders = fCard->GetVector("JetAssocPtBorders");
    fDeltaRBorders = fCard->GetVector("DeltaRBorders");

	fEfficiency = new AliJEfficiency();
    // 0:NoEff, 1:Period 2:RunNum 3:Auto
	fEfficiency->SetMode( fCard->Get("EfficiencyMode") );
    // Efficiency root file location local or alien
	fEfficiency->SetDataPath("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data"); 

    TRegexp reg("R[0-9][0-9][0-9]");
    TRegexp reg2("[0-9][0-9][0-9]");

    //container name has information of cone size like **R040**
    //this cone size information will be pulled to a numerical variable
    nJetContainer = fJetFinderName.size();
    fJetBgListOfList.resize(nJetContainer, TClonesArray("AliJJet",100));
	for (int i=0; i<nJetContainer; i++){
        TString fullNameOfiJetContainer(fJetFinderName[i]);
        TString coneSizeName (fullNameOfiJetContainer(reg));
        TString coneSizeValue (coneSizeName(reg2));
        fConeSizes.push_back( (double) coneSizeValue.Atoi()/100.);
	}


    int NBINS=150;
    double LogBinsX[NBINS+1], LimL=0.1, LimH=150;
    double logBW = (log(LimH)-log(LimL))/NBINS;
    for(int ij=0;ij<=NBINS;ij++) LogBinsX[ij]=LimL*exp(ij*logBW);



    fHMG = new AliJHistManager( "AliJJetJtHistManager");
    fJetFinderBin .Set("JetFinderOrder","NFin","NFin:%d", AliJBin::kSingle).SetBin(nJetContainer);
    fJetTriggerBin .Set("JetTriggerBin","JetPt","p_{T,jet} : %.1f - %.1f").SetBin(fCard->GetVector("JetTriggPtBorders"));
    fTrkPtBin .Set("TrkPtBin","TrkPt","p_{T,constituent}:%.1f-%.1f").SetBin(fCard->GetVector("JetAssocPtBorders"));
    fTrkLimPtBin .Set("TrkLimitPtBin","TrkLimitPt","p_{T,Limit}<%.1f", AliJBin::kSingle).SetBin(fJetConstPtLowLimits->GetNoElements());
    fdRBin.Set("dRBin","dR","dR : %.1f - %.1f ").SetBin(fCard->GetVector("DeltaRBorders"));
    fiHist.Set("iHist","iHist","iHist : %d ", AliJBin::kSingle).SetBin(10);

    fhNumber
        << TH1D("hNumber","Number",6,0,6)
        << fJetFinderBin
        <<"END";
    fhKNumber
        << TH1D("hKNumber","KNumber",17,0,17)
        << fJetFinderBin
        <<"END";

    fhJetPt 
        << TH1D("JetPt","",NBINS, LogBinsX ) 
        << fJetFinderBin
        <<"END";
    fhJetPtBin 
        << TH1D("JetPtBin","",NBINS, LogBinsX ) 
        << fJetFinderBin << fJetTriggerBin
        <<"END";

    int NBINSZ=150;
    double LogBinsZ[NBINSZ+1], LimLZ=0.001, LimHZ=1.1;
    double logBWZ = (TMath::Log(LimHZ)-TMath::Log(LimLZ))/NBINSZ;
    for(int ij=0;ij<=NBINSZ;ij++) LogBinsZ[ij]=LimLZ*exp(ij*logBWZ);//

    fhZ 
        << TH1D("Z","",NBINSZ, LogBinsZ ) 
        << fJetFinderBin
        <<"END";
    fhZBin 
        << TH1D("ZBin","",NBINSZ, LogBinsZ ) 
        << fJetFinderBin << fJetTriggerBin
        <<"END";

    int NBINSJt=150;
    double LogBinsJt[NBINSJt+1], LimLJt=0.01, LimHJt=10;
    double logBWJt = (TMath::Log(LimHJt)-TMath::Log(LimLJt))/NBINSJt;
    for(int ij=0;ij<=NBINSJt;ij++) LogBinsJt[ij]=LimLJt*exp(ij*logBWJt);
    int NBINSJtW=150;
    double LimLJtW=TMath::Log(0.01), LimHJtW=TMath::Log(10);

    fhJt 
        << TH1D("Jt","",NBINSJt, LogBinsJt ) 
        << fJetFinderBin
        <<"END";
    fhJtBin 
        << TH1D("JtBin","",NBINSJt, LogBinsJt ) 
        << fJetFinderBin << fJetTriggerBin
        <<"END";
    fhJtWeightBin 
        << TH1D("JtWeightBin","",NBINSJt, LogBinsJt ) 
        << fJetFinderBin << fJetTriggerBin
        <<"END";
    fhLogJtWeightBin 
        << TH1D("LogJtWeightBin","",NBINSJtW, LimLJtW, LimHJtW ) 
        << fJetFinderBin << fJetTriggerBin
        <<"END";

    fhJtWithPtCutWeightBinBin 
        << TH1D("JtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
        << fJetFinderBin << fJetTriggerBin << fTrkPtBin
        <<"END";
    fhLogJtWithPtCutWeightBinBin 
        << TH1D("LogJtWeightBinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
        << fJetFinderBin << fJetTriggerBin << fTrkPtBin
        <<"END";

    fhJtBinLimBin 
        << TH1D("JtBinLimBin","",NBINSJt, LogBinsJt ) 
        << fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
        <<"END";
    fhJtWeightBinLimBin 
        << TH1D("JtWeightBinLimBin","",NBINSJt, LogBinsJt ) 
        << fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
        <<"END";
    fhLogJtWeightBinLimBin 
        << TH1D("LogJtWeightBinLimBin","",NBINSJtW, LimLJtW, LimHJtW ) 
        << fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
        <<"END";

    fhJetBgPt 
        << TH1D("JetBgPt","",NBINS, LogBinsX ) 
        << fJetFinderBin
        <<"END";
    fhJetBgPtBin 
        << TH1D("JetBgPtBin","",NBINS, LogBinsX ) 
        << fJetFinderBin << fJetTriggerBin
        <<"END";
    fhBgZ 
        << TH1D("BgZ","",NBINSZ, LogBinsZ ) 
        << fJetFinderBin
        <<"END";
    fhBgZBin 
        << TH1D("BgZBin","",NBINSZ, LogBinsZ ) 
        << fJetFinderBin << fJetTriggerBin
        <<"END";


    fhBgJt 
        << TH1D("BgJt","",NBINSJt, LogBinsJt ) 
        << fJetFinderBin
        <<"END";
    fhBgJtBin 
        << TH1D("BgJtBin","",NBINSJt, LogBinsJt ) 
        << fJetFinderBin << fJetTriggerBin
        <<"END";
    fhBgJtWeightBin 
        << TH1D("BgJtWeightBin","",NBINSJt, LogBinsJt ) 
        << fJetFinderBin << fJetTriggerBin
        <<"END";
    fhBgLogJtWeightBin 
        << TH1D("BgLogJtWeightBin","",NBINSJtW, LimLJtW, LimHJtW ) 
        << fJetFinderBin << fJetTriggerBin
        <<"END";

    fhBgJtWithPtCutWeightBinBin 
        << TH1D("BgJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
        << fJetFinderBin << fJetTriggerBin << fTrkPtBin
        <<"END";
    fhBgLogJtWithPtCutWeightBinBin 
        << TH1D("BgLogJtWeightBinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
        << fJetFinderBin << fJetTriggerBin << fTrkPtBin
        <<"END";

    fhBgJtWithPtCutWeightBinBinSmallerR
        << TH1D("BgJtWithPtCutWeightBinBinSmallerR","",NBINSJt, LogBinsJt ) 
        << fiHist << fJetTriggerBin << fTrkPtBin 
        <<"END";
    fhBgLogJtWithPtCutWeightBinBinSmallerR 
        << TH1D("BgLogJtWeightBinBinBinSmallerR","",NBINSJtW, LimLJtW, LimHJtW ) 
        << fiHist << fJetTriggerBin << fTrkPtBin 
        <<"END";

    fhBgJtWithPtCutWeightBinBinDiffR
        << TH1D("BgJtWithPtCutWeightBinBinDiffR","",NBINSJt, LogBinsJt ) 
        << fiHist << fJetTriggerBin << fTrkPtBin 
        <<"END";
    fhBgLogJtWithPtCutWeightBinBinDiffR 
        << TH1D("BgLogJtWeightBinBinBinDiffR","",NBINSJtW, LimLJtW, LimHJtW ) 
        << fiHist << fJetTriggerBin << fTrkPtBin 
        <<"END";

    fhBgJtBinLimBin 
        << TH1D("BgJtBinLimBin","",NBINSJt, LogBinsJt ) << fJetFinderBin 
        << fJetTriggerBin << fTrkLimPtBin
        <<"END";
    fhBgJtWeightBinLimBin 
        << TH1D("BgJtWeightBinLimBin","",NBINSJt, LogBinsJt ) 
        << fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
        <<"END";
    fhBgLogJtWeightBinLimBin 
        << TH1D("BgLogJtWeightBinLimBin","",NBINSJtW, LimLJtW, LimHJtW ) 
        << fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
        <<"END";

    
    int NBINSdeltaN=40;
    double LimLdeltaN=-19.5, LimHdeltaN=19.5;
   
    fhdeltaN
    << TH1D("hdeltaN","",NBINSdeltaN,LimLdeltaN,LimHdeltaN )
    << fJetTriggerBin << fdRBin <<"END";

    int NBINSdeltaE=400;
    double LimLdeltaE=-20, LimHdeltaE=20;

    fhdeltaE
    << TH1D("hdeltaE","",NBINSdeltaE,LimLdeltaE,LimHdeltaE )
    << fJetTriggerBin << fdRBin <<"END";
    fhFullJetEChJetBin 
        << TH1D("hFullJetEChJetBin","",NBINS, LogBinsX )  << fJetTriggerBin
        <<"END";

    int nDR = 1000;double xDR0= -10; double xDR1 = 10;
    fhFullChdRChJetBin 
        << TH1D("hFullChdRChJetBin","",nDR,xDR0,xDR1)  << fJetTriggerBin
        <<"END";
    fh2DFullEvsChEdN0
        << TH2D("h2DFullEvsChEdN0","",NBINS, LogBinsX, NBINS, LogBinsX )  
        <<"END";
    fh2DFullEvsChEdNnot0
        << TH2D("h2DFullEvsChEdNnot0","",NBINS, LogBinsX, NBINS, LogBinsX )  
        <<"END";

    fHMG->Print();
    fHMG->WriteConfig();




}

void AliJJetJtAnalysis::ClearBeforeEvent(){
    //fJetListOfList.Clear();


}

void AliJJetJtAnalysis::UserExec(){
    for( int i=0;i<fJetListOfList.GetEntries();i++ ){
        TObjArray * Jets = (TObjArray*) fJetListOfList[i];
        if(!Jets) {
            continue;
        }
        this->FillJtHistogram(Jets,i);
    }

    //The Function should be called after calling "FillJtHistogram"
    //FillBgJtWithSmallerR(Bg jet finder array, old R, new R )
    //fJetBgListOfList[i] where {i=0-5;full0.4,full0.5,full0.6,Ch0.4,Ch0.5,Ch0.6}
    //Caution!! these array number should be changed WHEN jet finders change
    this->FillBgJtWithSmallerR(fJetBgListOfList[1], 1, 0.4,0);
    this->FillBgJtWithSmallerR(fJetBgListOfList[2], 2, 0.4,1);
    this->FillBgJtWithSmallerR(fJetBgListOfList[2], 2, 0.5,2);
    this->FillBgJtWithSmallerR(fJetBgListOfList[4], 4, 0.4,3);
    this->FillBgJtWithSmallerR(fJetBgListOfList[5], 5, 0.4,4);
    this->FillBgJtWithSmallerR(fJetBgListOfList[5], 5, 0.5,5);
   
    //Fill jt with diff cone axes (old axis iContainer, new axis, iHist) 
    this->FillBgJtWithDiffAxes(1, 0,0);
    this->FillBgJtWithDiffAxes(2, 0,1);
    this->FillBgJtWithDiffAxes(2, 1,2);
    this->FillBgJtWithDiffAxes(4, 3,0);
    this->FillBgJtWithDiffAxes(5, 3,1);
    this->FillBgJtWithDiffAxes(5, 4,2);
    //End.

    int iS1 = 0; //full 0.4
    int iS2 = 3; //Ch   0.4
    TObjArray * jetfinder1 = (TObjArray*) fJetListOfList[iS1];
    TObjArray * jetfinder2 = (TObjArray*) fJetListOfList[iS2];
    AliJJet *jet1 = NULL;
    AliJJet *jet2 = NULL;
    double deltaeta; 
    int chEbin=-1, rbin=-1;
    int dN=-1000;
    double dE=-1000.;
    for (int ijet = 0; ijet<jetfinder1->GetEntriesFast(); ijet++){
        jet1 = dynamic_cast<AliJJet*>( jetfinder1->At(ijet) );
        if (!jet1) continue;
        for (int jjet = 0; jjet<jetfinder2->GetEntriesFast(); jjet++){
            jet2 = dynamic_cast<AliJJet*>( jetfinder2->At(jjet) );
            if (!jet2) continue;
            chEbin = GetBin(fJetTriggPtBorders,jet2->E());
            deltaeta = TMath::Abs(jet1->Eta()-jet2->Eta());
            rbin   = GetBin(fDeltaRBorders,deltaeta); 
            fJJetAnalysis->CompareTwoJets(jet1, jet2, dE, dN);
            if (chEbin < 0 || rbin < 0 ) continue;
            fhdeltaE[chEbin][rbin]->Fill(dE);
            fhdeltaN[chEbin][rbin]->Fill(dN);
            if (dN ==0) { 
                fhFullJetEChJetBin[chEbin]->Fill(jet1->E());
                fhFullChdRChJetBin[chEbin]->Fill(jet1->DeltaR(*jet2));
                fh2DFullEvsChEdN0->Fill(jet1->E(), jet2->E());
            } else { 
                fh2DFullEvsChEdNnot0->Fill(jet1->E(), jet2->E());
            }

        }
    }


}

void AliJJetJtAnalysis::WriteHistograms(){


    TDirectory * cwd = gDirectory;
    //const int nJetContainer = fJetListOfList.GetEntries();


    for (int i=0; i<nJetContainer; i++){
        TDirectory *nwd = gDirectory->mkdir(fJetFinderName[i]);
        //Under the folder name, save objects
        //nwd->cd();
        //cwd->cd();
    }


}



void AliJJetJtAnalysis::FillJtHistogram( TObjArray *Jets , int iContainer)
{	

    int iBin, iptaBin=0;
    int jBin=0;
    double pT = 0;
    double conPtMax =0;

    double z; double jt;
    double pta;
    //double Y , deltaY = 0;
    //double Phi, deltaPhi;
    //double deltaR= 0;
    //cout<<"histogram filling number of jets : "<<Jets->GetEntriesFast()<<endl;

    TLorentzVector  vOrtho;
    fJetBgListOfList[iContainer].Clear();
    TClonesArray & bgjets = fJetBgListOfList[iContainer];



    int k = 0;
    double deltaR = -1;
    double deltaEta = -999;
    double deltaPhi = -999;
    double effCorrection = -1;
    double thisConeSize = fConeSizes[iContainer] ;
    int iBgJet = 0;

    // iJet loop for an event
    for (int i = 0; i<Jets->GetEntries(); i++){
        AliJJet *jet = dynamic_cast<AliJJet*>( Jets->At(i) );
        pT = jet->Pt();
        if (pT<(*fJetTriggPtBorders)[1]) continue;
        iBin = GetBin(fJetTriggPtBorders,pT); // fill jetPt histos
        if( iBin < 0 ) continue;
        fhJetPt[iContainer]->Fill( pT );
        fhJetPtBin[iContainer][iBin]->Fill( pT );

        for (int icon = 0; icon<jet->GetConstituents()->GetEntries(); icon++){
            AliJBaseTrack *con = jet->GetConstituent(icon);
            if (con->Pt()>conPtMax) conPtMax = con->Pt();
        }

        for (int ii = fJetConstPtLowLimits->GetNoElements(); ii >= 1 ; ii--){   
            if (conPtMax > (*fJetConstPtLowLimits)[ii]) {               
                jBin = ii-1;                                               
                break;
            }
        }

        //iConstituent loop for the iJet
        //jt, z are calcualted and filled  
        for (int icon = 0; icon<jet->GetConstituents()->GetEntries(); icon++){
            AliJBaseTrack *constituent = jet->GetConstituent(icon);
            z = (constituent->Vect()*jet->Vect().Unit())/jet->P();
            pta = constituent->Pt();
            constituent->SetTrackEff( fEfficiency->GetCorrection( pta, 5, fcent) );
            effCorrection = 1.0/constituent->GetTrackEff();
            iptaBin = GetBin(fJetAssocPtBorders, pta);
            if( iptaBin < 0 ) continue;


            fhZ[iContainer]->Fill( z , effCorrection);
            fhZBin[iContainer][iBin]->Fill( z , effCorrection);
            jt = (constituent->Vect()-z*jet->Vect()).Mag();
            fhJt[iContainer]->Fill( jt , effCorrection);
            fhJtBin[iContainer][iBin]->Fill( jt , effCorrection);
            fhJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt * effCorrection );
            fhLogJtWeightBin[iContainer][iBin]
                ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );



            fhJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
                ->Fill( jt, 1.0/jt * effCorrection );
            fhLogJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
                ->Fill( TMath::Log(jt), 1.0/jt * effCorrection);

            for (int jj = 0; jj <= jBin ; jj++) {
                fhJtBinLimBin[iContainer][iBin][jj]->Fill( jt, effCorrection );
                fhJtWeightBinLimBin[iContainer][iBin][jj]
                    ->Fill( jt, 1.0/jt * effCorrection );
                fhLogJtWeightBinLimBin[iContainer][iBin][jj]
                    ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
            }

        }



        vOrtho.SetVect(jet->Vect().Orthogonal());
        vOrtho.SetE(jet->E());
        //if (Log) cout<<"Before R caluation, R = "<<TMath::Sqrt(it->area()/TMath::Pi())<<endl;
        //R_area = TMath::Sqrt(it->area()/TMath::Pi())*Rs[order]/R;
        //if (Log )cout<<"Rs[order] = "<<Rs[order]<<" R = "<<R<<" Bg R area : "<<R_area<<endl;

        //Background jet (iBgJet) will be produced. This background jet is orthogonal to the iJet.  
        //If there is another jJet, then iBgJet will be consecutevely moved not to 
        //have jJet in the cone size. 
        if (Jets->GetEntries()>1){
            fhNumber[iContainer]->Fill(3.5);
            for (int j = 0; j<Jets->GetEntries(); j++){
                if (i == j) continue;
                AliJJet *jet2 = dynamic_cast<AliJJet*>( Jets->At(j) );

                if (k>15) {
                    fhNumber[iContainer]->Fill(5.5);
                    break;
                }


                deltaEta = vOrtho.Eta() - jet2->Eta();
                deltaPhi = vOrtho.Phi() - jet2->Phi();
                deltaR   = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
                if ( deltaR < thisConeSize) {

                    vOrtho.Rotate(TMath::Pi()/8, jet->Vect());
                    j=0;
                    k++;
                    fhNumber[iContainer]->Fill(4.5);
                }

            }
            fhKNumber[iContainer]->Fill(k);
        }

        // Filling iBgJet,  Bgjt and Bgz
        // "k<16" means that we will select a iBgJet which hasn't moved 
        // more than 16 times by the process above
        double maxconpt = 0;
        if ( k<16 ){
            new (bgjets[iBgJet]) AliJJet(vOrtho.Px(),vOrtho.Py(), vOrtho.Pz(), vOrtho.E(), i ,0,0);
            AliJJet * jbg = (AliJJet*) fJetBgListOfList[iContainer][iBgJet];
            iBgJet++;

            pT = vOrtho.Pt(); 
            if (pT<(*fJetTriggPtBorders)[1]) continue;

            fhJetBgPt[iContainer]->Fill( pT );
            //bbfHistos[iContainer]->fhJetBgPtWeight->Fill( pT, 1./pT);
            iBin = GetBin(fJetTriggPtBorders, pT);
            if( iBin < 0 ) continue;
            fhJetBgPtBin[iContainer][iBin]->Fill( pT );

            
            for (int icon = 0; icon<fTracks->GetEntries(); icon++){
                AliJBaseTrack *track = dynamic_cast<AliJBaseTrack*>(fTracks->At(icon));
                if (!track) continue;
                deltaEta = vOrtho.Eta() - track->Eta();
                deltaPhi = vOrtho.Phi() - track->Phi();
                deltaR   = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
                if ( deltaR > thisConeSize) continue;
        
                jbg->AddConstituent(track);

                pta = track->Pt();
                if (pta > maxconpt) maxconpt = pta;
                track->SetTrackEff( fEfficiency->GetCorrection( pta, 5, fcent) );
                effCorrection = 1.0/track->GetTrackEff();
                iptaBin = GetBin(fJetAssocPtBorders, pta);
                if( iptaBin < 0 ) continue;


                z = (track->Vect()*vOrtho.Vect().Unit())/vOrtho.P();
                fhBgZ[iContainer]->Fill( z , effCorrection);
                fhBgZBin[iContainer][iBin]->Fill( z , effCorrection);

                jt = (track->Vect()-z*vOrtho.Vect()).Mag();
                fhBgJt[iContainer]->Fill( jt , effCorrection);
                fhBgJtBin[iContainer][iBin]->Fill( jt , effCorrection);
                fhBgJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt * effCorrection );
                fhBgLogJtWeightBin[iContainer][iBin]->Fill( TMath::Log(jt), 1.0/jt * effCorrection );

                if (iptaBin < 0) continue;
                fhBgJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
                    ->Fill( jt, 1.0/jt * effCorrection );
                fhBgLogJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
                    ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
            }
            
            for (int ii = fJetConstPtLowLimits->GetNoElements(); ii >= 1 ; ii--)   {     
                if (maxconpt > (*fJetConstPtLowLimits)[ii]) {   
                    jBin = ii-1;
                    break;
                }
            }
            for (int icon =0; icon<jbg->GetConstituents()->GetEntries();icon++){
                AliJBaseTrack *con = jbg->GetConstituent(icon);
                z = (con->Vect()*jbg->Vect().Unit())/jbg->P();
                pta = con->Pt();
                iptaBin = GetBin(fJetAssocPtBorders, pta);
                jt = (con->Vect()-z*jbg->Vect()).Mag();
                if( iptaBin < 0 ) continue;
                for (int jj = 0; jj <= jBin ; jj++) {
                    fhBgJtBinLimBin[iContainer][iBin][jj]
                        ->Fill( jt, effCorrection );
                    fhBgJtWeightBinLimBin[iContainer][iBin][jj]
                        ->Fill( jt, 1.0/jt * effCorrection );
                    fhBgLogJtWeightBinLimBin[iContainer][iBin][jj]
                        ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
                }
            }
            cout<<"check : trackN, contbgN = "<<  
            fTracks->GetEntries() <<" "<<
            jbg->GetConstituents()->GetEntries()<<endl;   
        }


    }
}

void AliJJetJtAnalysis::FillBgJtWithSmallerR(const TClonesArray &Jets,int iContainer, double nR, int iHist){
    double iBin = -1, iptaBin = -1;
    double pT=-1, z=-1,jt=-1,  pta=-1;
    double effCorrection = -1;
    for (int i = 0; i<Jets.GetEntries(); i++){
        AliJJet *jet = dynamic_cast<AliJJet*>( Jets.At(i) );
        pT = jet->Pt();
        if (pT<(*fJetTriggPtBorders)[1]) continue;
        iBin = GetBin(fJetTriggPtBorders,pT); // fill jetPt histos
        if( iBin < 0 ) continue;

        for (int icon = 0; icon<jet->GetConstituents()->GetEntries(); icon++){
            AliJBaseTrack *con = jet->GetConstituent(icon);
            z = (con->Vect()*jet->Vect().Unit())/jet->P();
            pta = con->Pt();
            con->SetTrackEff( fEfficiency->GetCorrection( pta, 5, fcent) );
            effCorrection = 1.0/con->GetTrackEff();
            iptaBin = GetBin(fJetAssocPtBorders, pta);
            if( iptaBin < 0 ) continue;
            if (jet->DeltaR(*con)>nR){
                cout<<" old Bg R  : "<<fConeSizes[iContainer]
                    <<" jet con R : "<<jet->DeltaR(*con)
                    <<endl;
                continue;

            }
            jt = (con->Vect()-z*jet->Vect()).Mag();
            fhBgJtWithPtCutWeightBinBinSmallerR[iHist][iBin][iptaBin]
                ->Fill( jt, 1.0/jt * effCorrection );
            fhBgLogJtWithPtCutWeightBinBinSmallerR[iHist][iBin][iptaBin]
                ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
        }


    }


}


void AliJJetJtAnalysis::FillBgJtWithDiffAxes (
      int iao
    , int ian
    , int iHist
){

    const TClonesArray &ao = fJetBgListOfList[iao];
    const TClonesArray &an = fJetBgListOfList[ian];

    double iBin = -1, iptaBin = -1;
    double pT=-1, z=-1,jt=-1,  pta=-1;
    double effCorrection = -1;
     
    for (int io = 0; io<ao.GetEntries(); io++){
        AliJJet *jo = dynamic_cast<AliJJet*>( ao.At(io) );
        for (int in = 0; in<an.GetEntries(); in++){
            AliJJet *jn = dynamic_cast<AliJJet*>( an.At(in) );
            if (jo->DeltaR(*jn) > fConeSizes[ian]) {
                continue;
            } else{
                cout << "iao ian new axis delta R "
                     << iao <<" "
                     << ian <<" "
                     << jo->DeltaR(*jn)
                     <<endl;
            }
            pT = jn->Pt();
            if (pT<(*fJetTriggPtBorders)[1]) continue;
            iBin = GetBin(fJetTriggPtBorders,pT); // fill jetPt histos
            if( iBin < 0 ) continue;

            for (int ic = 0; ic<jo->GetConstituents()->GetEntries(); ic++){
                AliJBaseTrack *con = jo->GetConstituent(ic);
                if (jn->DeltaR(*con) > fConeSizes[ian]) continue;
                z = (con->Vect()*jn->Vect().Unit())/jn->P();
                pta = con->Pt();
                con->SetTrackEff( fEfficiency->GetCorrection( pta, 5, fcent) );
                effCorrection = 1.0/con->GetTrackEff();
                iptaBin = GetBin(fJetAssocPtBorders, pta);
                if( iptaBin < 0 ) continue;
                jt = (con->Vect()-z*jn->Vect()).Mag();
                fhBgJtWithPtCutWeightBinBinDiffR[iHist][iBin][iptaBin]
                    ->Fill( jt, 1.0/jt * effCorrection );
                fhBgLogJtWithPtCutWeightBinBinDiffR[iHist][iBin][iptaBin]
                    ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
            }
        }
    }

}

