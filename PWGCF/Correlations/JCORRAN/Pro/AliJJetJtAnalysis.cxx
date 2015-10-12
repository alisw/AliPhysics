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
	, fJetListOfList() // Jet finders container
	//, fJetBgList(NULL) 
	, fJetBgListOfList() // Bg jet list of list container
	, fJetEtaCut() // Jet's eta cut
	, fJetTriggPtBorders(NULL)
	, fJetConstPtLowLimits(NULL)
	, fJetAssocPtBorders(NULL)
	, fDeltaRBorders(NULL)
	, nJetContainer(0) //number of Jets finders
	, fCard(NULL) // pointer to AliJCard
	, fJJetAnalysis(NULL) //pointer to AliJJetAnalysis
	, fJetFinderName(0) // string array for jet finder tags
	, fConeSizes(0)    //array of cone sizes of Jets finder
	, fEfficiency(0x0) //pointer to tracking efficiency
	, cBin(-1)
	, fcent(-999) // values for centrality
	, zBin(-1)
	, zVert(-999)
	, fTracks(NULL) //list of tracks
	, Nrandom(1)
	, fHMG(NULL)
	, fJetFinderBin() //histmanager axis
	, fJetTriggerBin() //histmanager axis
	, fTrkPtBin() //histmanager axis
	, fTrkLimPtBin() //histmanager axis
	, fdRBin() //histmanager axis
	, fiHist() //histmanager axis
	, fhNumber() // number of jets
	, fhKNumber() // number of axis rotation
	, fhJetPt() //Jet pt distribution
	, fhJetPtBin() //Jet pt as a fn of jet pt 
	, fhJetPtWeight() //Jet pt distribution
	, fhJetPtWeightBin() //Jet pt as a fn of jet pt 
	, fhZ() // z dist
	, fhZBin() // z dist as fn of jet pt
	, fhJt() // jt dist
	, fhJtBin() // jt as fn of jet pt
	, fhJtWeightBin() // jt in P.S as fn of jet pt
	, fhLogJtWeightBin() //log jt in P.S as fn of jet pt
	, fhLogJtWeight2Bin() //log jt in P.S as fn of jet pt
	, fhJtWithPtCutWeightBinBin() // jt in P.S as fns of jet pt, const pt
	, fhLogJtWithPtCutWeightBinBin()
	, fhLogJtWithPtCutWeight2BinBin()
	, fhJetConeTrkPt()
	, fhJetConeTrkPtBin()
	, fhJetConeTrkPtWeightBin()
	, fhJetConeZ()
	, fhJetConeZBin()
	, fhJetConeJt()
	, fhJetConeJtBin()
	, fhJetConeJtWeightBin()
	, fhJetConeLogJtWeightBin()
	, fhJetConeLogJtWeight2Bin()
	, fhJetConeJtWithPtCutWeightBinBin()
	, fhJetConeLogJtWithPtCutWeightBinBin()
	, fhJetConeLogJtWithPtCutWeight2BinBin()
	, fhJtBinLimBin() //jt with maximum const pt bin
	, fhJtWeightBinLimBin()
	, fhLogJtWeightBinLimBin()
	, fhLogJtWeight2BinLimBin()
	, fhJetBgPt()
	, fhJetBgPtBin()
	, fhBgZ() //Bg Z 
	, fhBgZBin() //Bg Z as fn of jet pt
	, fhBgJt() //Bg jt 
	, fhBgJtBin() //Bg jt as fn of jet pt
	, fhBgJtWeightBin() //Bg jt in PS as fn of jet pt
	, fhBgLogJtWeightBin() // log of Bg jt as fn of jet pt
, fhBgLogJtWeight2Bin() // log of Bg jt as fn of jet pt
	// Bg jt as fns of jet pt, track pt
, fhBgJtWithPtCutWeightBinBin() 
	// log of Bg jt in PS as fns of jet pt, track pt
	, fhBgLogJtWithPtCutWeightBinBin() 
, fhBgLogJtWithPtCutWeight2BinBin() 
	// Bg jt in PS as fns of jet pt, track pt
	// if R = 0.6, all bg jts are filled for all R<0.6
	, fhBgJtWithPtCutWeightBinBinSmallerR() 
	, fhBgLogJtWithPtCutWeightBinBinSmallerR()
, fhBgLogJtWithPtCutWeight2BinBinSmallerR()
	// Bg jt in PS as fns of jet pt, track pt
	// from new jet axis with present tracks
	, fhBgJtWithPtCutWeightBinBinDiffR()
	, fhBgLogJtWithPtCutWeightBinBinDiffR()
, fhBgLogJtWithPtCutWeight2BinBinDiffR()
	// Bg jt as a fn of jet pt 
	// with maximum constituent's pT bin in a jet
	, fhBgJtBinLimBin()
	, fhBgJtWeightBinLimBin()
	, fhBgLogJtWeightBinLimBin()
	, fhBgLogJtWeight2BinLimBin()
	, fhTrkPt()
	, fhTrkPtBin()
	, fhTrkPtWeightBin()
	, fhBgTrkPt()
	, fhBgTrkPtBin()
, fhBgTrkPtWeightBin()
	, fhBgTrkNumber()
	, fhBgTrkNumberBin()
	, fhBgRndmTrkPt()
	, fhBgRndmTrkNumber()
	, fhBgRndmZ()
	, fhBgRndmJt()
	, fhBgRndmLogJt()
	, fhBgRndmJtWithPtCutWeightBin()
	, fhBgRndmLogJtWithPtCutWeight2Bin()
	, fhBgRndmJtWithPtCutWeightBinBin()
	, fhBgRndmLogJtWithPtCutWeight2BinBin()

	// dE = E of charged jet - E of shared tracks of
	// charged and full jets
, fhdeltaE()
	// dN = N of charged jet - N of shared tracks of
	// charged and full jets
, fhdeltaN()
	// Fulljet E as a fn of charged jet E 
, fhFullJetEChJetBin()
	// dR between full and charged jet 
	// as a fn of charged jet E
, fhFullChdRChJetBin()
	// x-axis: fulljet E, y: ch jet E when dN == 0
, fh2DFullEvsChEdN0()
	// x-axis: fulljet E, y: ch jet E when dN != 0
, fh2DFullEvsChEdNnot0()
	//Jet eta and Phi 2D dist for each jet finder
	, fhJetEtaPhi()
	, fhTrackEtaPhi()
, frandom(0x0)
{

	frandom = new TRandom3(); // frandom generator for jt flow UE
	frandom->SetSeed(0);
}

AliJJetJtAnalysis::AliJJetJtAnalysis( AliJCard * card ):
	fInputList(NULL)
	, fJetList(NULL)
	, fJetListOfList()
	, fJetBgListOfList()
	, fJetEtaCut() 
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
	, Nrandom(1)
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
	, fhJetPtWeight() //Jet pt distribution
	, fhJetPtWeightBin() //Jet pt as a fn of jet pt 
	, fhZ()
	, fhZBin()
	, fhJt()
	, fhJtBin()
	, fhJtWeightBin()
	, fhLogJtWeightBin()
	, fhLogJtWeight2Bin()
	, fhJtWithPtCutWeightBinBin()
	, fhLogJtWithPtCutWeightBinBin()
	, fhLogJtWithPtCutWeight2BinBin()
	, fhJtBinLimBin()
	, fhJtWeightBinLimBin()
	, fhLogJtWeightBinLimBin()
	, fhLogJtWeight2BinLimBin()
	, fhJetConeTrkPt()
	, fhJetConeTrkPtBin()
	, fhJetConeTrkPtWeightBin()
	, fhJetConeZ()
	, fhJetConeZBin()
	, fhJetConeJt()
	, fhJetConeJtBin()
	, fhJetConeJtWeightBin()
	, fhJetConeLogJtWeightBin()
	, fhJetConeLogJtWeight2Bin()
	, fhJetConeJtWithPtCutWeightBinBin()
	, fhJetConeLogJtWithPtCutWeightBinBin()
	, fhJetConeLogJtWithPtCutWeight2BinBin()
	, fhJetBgPt()
	, fhJetBgPtBin()
	, fhBgZ()
	, fhBgZBin()
	, fhBgJt()
	, fhBgJtBin()
	, fhBgJtWeightBin()
	, fhBgLogJtWeightBin()
	, fhBgLogJtWeight2Bin()
	, fhBgJtWithPtCutWeightBinBin()
	, fhBgLogJtWithPtCutWeightBinBin()
	, fhBgLogJtWithPtCutWeight2BinBin()
	, fhBgJtWithPtCutWeightBinBinSmallerR()
	, fhBgLogJtWithPtCutWeightBinBinSmallerR()
	, fhBgLogJtWithPtCutWeight2BinBinSmallerR()
	, fhBgJtWithPtCutWeightBinBinDiffR()
	, fhBgLogJtWithPtCutWeightBinBinDiffR()
	, fhBgLogJtWithPtCutWeight2BinBinDiffR()
	, fhBgJtBinLimBin()
	, fhBgJtWeightBinLimBin()
	, fhBgLogJtWeightBinLimBin()
	, fhBgLogJtWeight2BinLimBin()
	, fhBgRndmTrkPt()
	, fhBgRndmTrkNumber()
	, fhBgRndmZ()
	, fhBgRndmJt()
	, fhBgRndmLogJt()
	, fhBgRndmJtWithPtCutWeightBin()
	, fhBgRndmLogJtWithPtCutWeight2Bin()
	, fhBgRndmJtWithPtCutWeightBinBin()
	, fhBgRndmLogJtWithPtCutWeight2BinBin()
	, fhTrkPt()
	, fhTrkPtBin()
	, fhTrkPtWeightBin()
	, fhBgTrkPt()
	, fhBgTrkPtBin()
	, fhBgTrkPtWeightBin()
	, fhBgTrkNumber()
	, fhBgTrkNumberBin()
	, fhdeltaE()
	, fhdeltaN()
	, fhFullJetEChJetBin()
	, fhFullChdRChJetBin()
	, fh2DFullEvsChEdN0()
	, fh2DFullEvsChEdNnot0()
	, fhJetEtaPhi()
	, fhTrackEtaPhi()
	, frandom(0x0)
{

	frandom = new TRandom3(); // frandom generator for jt flow UE
	frandom->SetSeed(0);
}

AliJJetJtAnalysis::AliJJetJtAnalysis(const AliJJetJtAnalysis& ap) :
	fInputList(ap.fInputList)
	, fJetList(ap.fJetList)
	, fJetListOfList(ap.fJetListOfList)
	//, fJetBgList(ap.fJetBgList)
	, fJetBgListOfList(ap.fJetBgListOfList)
	, fJetEtaCut(ap.fJetEtaCut) // Jet's eta cut
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
	, Nrandom(ap.Nrandom)
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
	, fhJetPtWeight(ap.fhJetPtWeight)
	, fhJetPtWeightBin(ap.fhJetPtWeightBin)
	, fhZ(ap.fhZ)
	, fhZBin(ap.fhZBin)
	, fhJt(ap.fhJt)
	, fhJtBin(ap.fhJtBin)
	, fhJtWeightBin(ap.fhJtWeightBin)
	, fhLogJtWeightBin(ap.fhLogJtWeightBin)
	, fhLogJtWeight2Bin(ap.fhLogJtWeight2Bin)
	, fhJtWithPtCutWeightBinBin(ap.fhJtWithPtCutWeightBinBin)
	, fhLogJtWithPtCutWeightBinBin(ap.fhLogJtWithPtCutWeightBinBin)
	, fhLogJtWithPtCutWeight2BinBin(ap.fhLogJtWithPtCutWeight2BinBin)
	, fhJtBinLimBin(ap.fhJtBinLimBin)
	, fhJtWeightBinLimBin(ap.fhJtWeightBinLimBin)
	, fhLogJtWeightBinLimBin(ap.fhLogJtWeightBinLimBin)
	, fhLogJtWeight2BinLimBin(ap.fhLogJtWeight2BinLimBin)
	, fhJetConeTrkPt(ap.fhJetConeTrkPt)
	, fhJetConeTrkPtBin(ap.fhJetConeTrkPtBin)
	, fhJetConeTrkPtWeightBin(ap.fhJetConeTrkPtWeightBin)
	, fhJetConeZ(ap.fhJetConeZ)
	, fhJetConeZBin(ap.fhJetConeZBin)
	, fhJetConeJt(ap.fhJetConeJt)
	, fhJetConeJtBin(ap.fhJetConeJtBin)
	, fhJetConeJtWeightBin(ap.fhJetConeJtWeightBin)
	, fhJetConeLogJtWeightBin(ap.fhJetConeLogJtWeightBin)
	, fhJetConeLogJtWeight2Bin(ap.fhJetConeLogJtWeight2Bin)
	, fhJetConeJtWithPtCutWeightBinBin(ap.fhJetConeJtWithPtCutWeightBinBin)
	, fhJetConeLogJtWithPtCutWeightBinBin(ap.fhJetConeLogJtWithPtCutWeightBinBin)
	, fhJetConeLogJtWithPtCutWeight2BinBin(ap.fhJetConeLogJtWithPtCutWeight2BinBin)
	, fhJetBgPt(ap.fhJetBgPt)
	, fhJetBgPtBin(ap.fhJetBgPtBin)
	, fhBgZ(ap.fhBgZ)
	, fhBgZBin(ap.fhBgZBin)
	, fhBgJt(ap.fhBgJt)
	, fhBgJtBin(ap.fhBgJtBin)
	, fhBgJtWeightBin(ap.fhBgJtWeightBin)
	, fhBgLogJtWeightBin(ap.fhBgLogJtWeightBin)
	, fhBgLogJtWeight2Bin(ap.fhBgLogJtWeight2Bin)
	, fhBgJtWithPtCutWeightBinBin(ap.fhBgJtWithPtCutWeightBinBin)
	, fhBgLogJtWithPtCutWeightBinBin(ap.fhBgLogJtWithPtCutWeightBinBin)
	, fhBgLogJtWithPtCutWeight2BinBin(ap.fhBgLogJtWithPtCutWeight2BinBin)
	, fhBgJtWithPtCutWeightBinBinSmallerR(ap.fhBgJtWithPtCutWeightBinBinSmallerR)
	, fhBgLogJtWithPtCutWeightBinBinSmallerR(ap.fhBgLogJtWithPtCutWeightBinBinSmallerR)
	, fhBgLogJtWithPtCutWeight2BinBinSmallerR(ap.fhBgLogJtWithPtCutWeight2BinBinSmallerR)
	, fhBgJtWithPtCutWeightBinBinDiffR(ap.fhBgJtWithPtCutWeightBinBinDiffR)
	, fhBgLogJtWithPtCutWeightBinBinDiffR(ap.fhBgLogJtWithPtCutWeightBinBinDiffR)
	, fhBgLogJtWithPtCutWeight2BinBinDiffR(ap.fhBgLogJtWithPtCutWeight2BinBinDiffR)
	, fhBgJtBinLimBin(ap.fhBgJtBinLimBin)
	, fhBgJtWeightBinLimBin(ap.fhBgJtWeightBinLimBin)
	, fhBgLogJtWeightBinLimBin(ap.fhBgLogJtWeightBinLimBin)
	, fhBgLogJtWeight2BinLimBin(ap.fhBgLogJtWeight2BinLimBin)
	, fhTrkPt(ap.fhTrkPt)
	, fhTrkPtBin(ap.fhTrkPtBin)
	, fhTrkPtWeightBin(ap.fhTrkPtWeightBin)
	, fhBgTrkPt(ap.fhBgTrkPt)
	, fhBgTrkPtBin(ap.fhBgTrkPtBin)
	, fhBgTrkPtWeightBin(ap.fhBgTrkPtWeightBin)
	, fhBgTrkNumber(ap.fhBgTrkNumber)
	, fhBgTrkNumberBin(ap.fhBgTrkNumberBin)
	, fhBgRndmTrkPt(ap.fhBgRndmTrkPt)
	, fhBgRndmTrkNumber(ap.fhBgRndmTrkNumber)
	, fhBgRndmZ(ap.fhBgRndmZ)
	, fhBgRndmJt(ap.fhBgRndmJt)
	, fhBgRndmLogJt(ap.fhBgRndmLogJt)
	, fhBgRndmJtWithPtCutWeightBin(ap.fhBgRndmJtWithPtCutWeightBin)
	, fhBgRndmLogJtWithPtCutWeight2Bin(ap.fhBgRndmLogJtWithPtCutWeight2Bin)
	, fhBgRndmJtWithPtCutWeightBinBin(ap.fhBgRndmJtWithPtCutWeightBinBin)
	, fhBgRndmLogJtWithPtCutWeight2BinBin(ap.fhBgRndmLogJtWithPtCutWeight2BinBin)
	, fhdeltaE(ap.fhdeltaE)
	, fhdeltaN(ap.fhdeltaN)
	, fhFullJetEChJetBin(ap.fhFullJetEChJetBin)
	, fhFullChdRChJetBin(ap.fhFullChdRChJetBin)
	, fh2DFullEvsChEdN0(ap.fh2DFullEvsChEdN0)
	, fh2DFullEvsChEdNnot0(ap.fh2DFullEvsChEdNnot0)
	, fhJetEtaPhi(ap.fhJetEtaPhi)
	, fhTrackEtaPhi(ap.fhTrackEtaPhi)
	, frandom(0x0)
{

	frandom = new TRandom3(); // frandom generator for jt flow UE
	frandom->SetSeed(0);
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
	fJetEtaCut = fCard-> Get("JetEtaCut");

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
		<< TH1D("hNumber","Number",20,0,20)
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


	fhJetPtWeight 
		<< TH1D("JetPtWeight","",NBINS, LogBinsX ) 
		<< fJetFinderBin
		<<"END";
	fhJetPtWeightBin 
		<< TH1D("JetPtWeightBin","",NBINS, LogBinsX ) 
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

	fhJetConeZ 
		<< TH1D("JetConeZ","",NBINSZ, LogBinsZ ) 
		<< fJetFinderBin
		<<"END";
	fhJetConeZBin 
		<< TH1D("JetConeZBin","",NBINSZ, LogBinsZ ) 
		<< fJetFinderBin << fJetTriggerBin
		<<"END";

	fhBgRndmZ 
		<< TH1D("BgRndmZ","",NBINSZ, LogBinsZ ) 
		<< fJetFinderBin
		<<"END";
	int NBINSJt=150;
	double LogBinsJt[NBINSJt+1], LimLJt=0.01, LimHJt=10;
	double logBWJt = (TMath::Log(LimHJt)-TMath::Log(LimLJt))/NBINSJt;
	for(int ij=0;ij<=NBINSJt;ij++) LogBinsJt[ij]=LimLJt*exp(ij*logBWJt);
	int NBINSJtW=150;
	double LimLJtW=TMath::Log(0.01), LimHJtW=TMath::Log(10);

	cout << "Debug 1" << endl;

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
	fhJetConeJt 
		<< TH1D("JetConeJt","",NBINSJt, LogBinsJt ) 
		<< fJetFinderBin
		<<"END";
	fhJetConeJtBin 
		<< TH1D("JetConeJtBin","",NBINSJt, LogBinsJt ) 
		<< fJetFinderBin << fJetTriggerBin
		<<"END";
	fhJetConeJtWeightBin 
		<< TH1D("JetConeJtWeightBin","",NBINSJt, LogBinsJt ) 
		<< fJetFinderBin << fJetTriggerBin
		<<"END";
	fhLogJtWeightBin 
		<< TH1D("LogJtWeightBin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin
		<<"END";
	fhJetConeLogJtWeightBin 
		<< TH1D("JetConeLogJtWeightBin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin
		<<"END";
	fhLogJtWeight2Bin 
		<< TH1D("LogJtWeight2Bin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin
		<<"END";
	fhJetConeLogJtWeight2Bin 
		<< TH1D("JetConeLogJtWeight2Bin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin
		<<"END";

	fhJtWithPtCutWeightBinBin 
		<< TH1D("JtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";
	fhJetConeJtWithPtCutWeightBinBin 
		<< TH1D("JetConeJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";
	fhLogJtWithPtCutWeightBinBin 
		<< TH1D("LogJtWithPtCutWeightBinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";
	fhLogJtWithPtCutWeight2BinBin 
		<< TH1D("LogJtWithPtCutWeight2BinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";

	fhJetConeLogJtWithPtCutWeightBinBin 
		<< TH1D("JetConeLogJtWithPtCutWeightBinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";
	fhJetConeLogJtWithPtCutWeight2BinBin 
		<< TH1D("JetConeLogJtWithPtCutWeight2BinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
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
	fhLogJtWeight2Bin 
		<< TH1D("LogJtWeight2Bin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin
		<<"END";

	fhJtWithPtCutWeightBinBin 
		<< TH1D("JtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";
	fhJetConeJtWithPtCutWeightBinBin 
		<< TH1D("JetConeJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";
	fhLogJtWithPtCutWeightBinBin 
		<< TH1D("LogJtWithPtCutWeightBinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";
	fhLogJtWithPtCutWeight2BinBin 
		<< TH1D("LogJtWithPtCutWeight2BinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";
	fhLogJtWeightBinLimBin 
		<< TH1D("LogJtWeightBinLimBin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
		<<"END";
	fhLogJtWeight2BinLimBin 
		<< TH1D("LogJtWeight2BinLimBin","",NBINSJtW, LimLJtW, LimHJtW ) 
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

	cout << "Debug 2" << endl;

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
	fhBgLogJtWeight2Bin 
		<< TH1D("BgLogJtWeight2Bin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin
		<<"END";

	fhBgJtWithPtCutWeightBinBin 
		<< TH1D("BgJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";
	fhBgLogJtWithPtCutWeightBinBin 
		<< TH1D("BgLogJtWithPtCutWeightBinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";
	fhBgLogJtWithPtCutWeight2BinBin 
		<< TH1D("BgLogJtWeight2BinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";

	fhBgJtWithPtCutWeightBinBinSmallerR
		<< TH1D("BgJtWithPtCutWeightBinBinSmallerR","",NBINSJt, LogBinsJt ) 
		<< fiHist << fJetTriggerBin << fTrkPtBin 
		<<"END";
	fhBgLogJtWithPtCutWeightBinBinSmallerR 
		<< TH1D("BgLogJtWithPtCutWeightBinBinSmallerR","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fiHist << fJetTriggerBin << fTrkPtBin 
		<<"END";
	fhBgLogJtWithPtCutWeight2BinBinSmallerR 
		<< TH1D("BgLogJtWithPtCutWeight2BinBinSmallerR","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fiHist << fJetTriggerBin << fTrkPtBin 
		<<"END";

	fhBgJtWithPtCutWeightBinBinDiffR
		<< TH1D("BgJtWithPtCutWeightBinBinDiffR","",NBINSJt, LogBinsJt ) 
		<< fiHist << fJetTriggerBin << fTrkPtBin 
		<<"END";
	fhBgLogJtWithPtCutWeightBinBinDiffR 
		<< TH1D("BgLogJtWithPtCutWeightBinBinDiffR","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fiHist << fJetTriggerBin << fTrkPtBin 
		<<"END";
	fhBgLogJtWithPtCutWeight2BinBinDiffR 
		<< TH1D("BgLogJtWithPtCutWeight2BinBinDiffR","",NBINSJtW, LimLJtW, LimHJtW ) 
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

	cout << "Debug 3" << endl;
	fhBgLogJtWeight2BinLimBin 
		<< TH1D("BgLogJtWeight2BinLimBin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
		<<"END";

	fhBgRndmJt 
		<< TH1D("BgRndmJt","",NBINSJt, LogBinsJt ) 
		<< fJetFinderBin
		<<"END";
	fhBgRndmLogJt 
		<< TH1D("BgRndmLogJt","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin
		<<"END";
	fhBgRndmJtWithPtCutWeightBin
		<< TH1D("BgRndmJtWithPtCutWeightBin","",NBINSJt, LogBinsJt ) 
		<< fJetFinderBin << fTrkPtBin
		<<"END";
	fhBgRndmLogJtWithPtCutWeight2Bin 
		<< TH1D("BgRndmLogJtWithPtCutWeight2Bin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fTrkPtBin
		<<"END";
	fhBgRndmJtWithPtCutWeightBinBin
		<< TH1D("BgRndmJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";
	fhBgRndmLogJtWithPtCutWeight2BinBin 
		<< TH1D("BgRndmLogJtWithPtCutWeight2BinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
		<< fJetFinderBin << fJetTriggerBin << fTrkPtBin
		<<"END";
	int NBINSPt=150;
	double LogBinsPt[NBINSPt+1], LimLPt=0.01, LimHPt=50;
	double logBWPt = (TMath::Log(LimHPt)-TMath::Log(LimLPt))/NBINSPt;
	for(int ij=0;ij<=NBINSPt;ij++) LogBinsPt[ij]=LimLPt*exp(ij*logBWPt);

	fhTrkPt
		<< TH1D("TrkPt","",NBINSPt,LogBinsPt)
		<< fJetFinderBin 
		<<"END";

	fhTrkPtBin
		<< TH1D("TrkPtBin","",NBINSPt,LogBinsPt)
		<< fJetFinderBin << fJetTriggerBin
		<<"END";

	fhTrkPtWeightBin
		<< TH1D("TrkPtWeightBin","",NBINSPt,LogBinsPt)
		<< fJetFinderBin << fJetTriggerBin
		<<"END";

	fhJetConeTrkPt
		<< TH1D("JetConeTrkPt","",NBINSPt,LogBinsPt)
		<< fJetFinderBin 
		<<"END";

	fhJetConeTrkPtBin
		<< TH1D("JetConeTrkPtBin","",NBINSPt,LogBinsPt)
		<< fJetFinderBin << fJetTriggerBin
		<<"END";

	fhJetConeTrkPtWeightBin
		<< TH1D("JetConeTrkPtWeightBin","",NBINSPt,LogBinsPt)
		<< fJetFinderBin << fJetTriggerBin
		<<"END";

	fhBgTrkPt
		<< TH1D("BgTrkPt","",NBINSPt,LogBinsPt)
		<< fJetFinderBin 
		<<"END";

	fhBgTrkPtBin
		<< TH1D("BgTrkPtBin","",NBINSPt,LogBinsPt)
		<< fJetFinderBin << fJetTriggerBin
		<<"END";

	fhBgTrkPtWeightBin
		<< TH1D("BgTrkPtWeightBin","",NBINSPt,LogBinsPt)
		<< fJetFinderBin << fJetTriggerBin
		<<"END";

	int NBINSNumber = 100;
	int LBinsNumber = 0;
	int HBinsNumber = 100;
	fhBgTrkNumber
		<< TH1D("BgTrkNumber","",NBINSNumber,LBinsNumber,HBinsNumber)
		<< fJetFinderBin 
		<<"END";
	fhBgTrkNumberBin
		<< TH1D("BgTrkNumberBin","",NBINSNumber,LBinsNumber,HBinsNumber)
		<< fJetFinderBin << fJetTriggerBin
		<<"END";


	fhBgRndmTrkPt
		<< TH1D("BgRndmTrkPt","",NBINSPt,LogBinsPt)
		<< fJetFinderBin 
		<<"END";
	fhBgRndmTrkNumber
		<< TH1D("BgRndmTrkNumber","",NBINSNumber,LBinsNumber,HBinsNumber)
		<< fJetFinderBin 
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

	fhJetEtaPhi
		<< TH2D("hJetEtaPhi","jet eta phi dist",100,-1.0,1.0,70,-1*TMath::Pi(),TMath::Pi())
		<< fJetFinderBin<<"END";

	fhTrackEtaPhi
		<< TH2D("hTrackEtaPhi","track eta phi dist",70,-1.,1.,70,-2*TMath::Pi(),2*TMath::Pi())
		<< fJetFinderBin<<"END";


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
	this->FillBgJtWithSmallerR(fJetBgListOfList[1], 0.4,0);
	this->FillBgJtWithSmallerR(fJetBgListOfList[2], 0.4,1);
	this->FillBgJtWithSmallerR(fJetBgListOfList[2], 0.5,2);
	this->FillBgJtWithSmallerR(fJetBgListOfList[4], 0.4,3);
	this->FillBgJtWithSmallerR(fJetBgListOfList[5], 0.4,4);
	this->FillBgJtWithSmallerR(fJetBgListOfList[5], 0.5,5);

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


// Fill jt to  histograms and calculate bg jt and fill it into histograms
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
	TLorentzVector  randomTrack;
	TLorentzVector summedJet;
	fJetBgListOfList[iContainer].Clear();
	TClonesArray & bgjets = fJetBgListOfList[iContainer];

	//cout << "REMOVE THIS " << endl << endl << endl;


	int k = 0;
	double deltaR = -1;
	double deltaEta = -999;
	double deltaPhi = -999;
	double effCorrection = -1;
	double thisConeSize = fConeSizes[iContainer] ;
	double etaRndm = 0;
	double phiRndm = 0;
	int iBgJet = 0;


	int doRndmBg = 0;
	// iJet loop for an event
	for (int i = 0; i<Jets->GetEntries(); i++){
		AliJJet *jet = dynamic_cast<AliJJet*>( Jets->At(i) );
		if (!jet) continue;
		if (TMath::Abs(jet->Eta()) > fJetEtaCut) continue;
		fhJetEtaPhi[iContainer]->Fill(jet->Eta(),jet->Phi());
		pT = jet->Pt();
		if (pT<(*fJetTriggPtBorders)[1]) continue;
		iBin = GetBin(fJetTriggPtBorders,pT); // fill jetPt histos
		if( iBin < 0 ) continue;
		doRndmBg++;
		fhJetPt[iContainer]->Fill( pT );
		fhJetPtBin[iContainer][iBin]->Fill( pT );
		fhJetPtWeight[iContainer]->Fill( pT , 1/pT);
		fhJetPtWeightBin[iContainer][iBin]->Fill( pT, 1/pT );

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

			fhTrkPt[iContainer]->Fill(pta,effCorrection);
			fhTrkPtBin[iContainer][iBin]->Fill(pta,effCorrection);
			fhTrkPtWeightBin[iContainer][iBin]->Fill(pta,effCorrection/pta);


			fhZ[iContainer]->Fill( z , effCorrection);
			fhZBin[iContainer][iBin]->Fill( z , effCorrection);
			jt = (constituent->Vect()-z*jet->Vect()).Mag();
			fhJt[iContainer]->Fill( jt , effCorrection);
			fhJtBin[iContainer][iBin]->Fill( jt , effCorrection);
			fhJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt * effCorrection );
			fhLogJtWeightBin[iContainer][iBin]
				->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
			fhLogJtWeight2Bin[iContainer][iBin]
				->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );



			fhJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
				->Fill( jt, 1.0/jt * effCorrection );
			fhLogJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
				->Fill( TMath::Log(jt), 1.0/jt * effCorrection);
			fhLogJtWithPtCutWeight2BinBin[iContainer][iBin][iptaBin]
				->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection);

			for (int jj = 0; jj <= jBin ; jj++) {
				fhJtBinLimBin[iContainer][iBin][jj]->Fill( jt, effCorrection );
				fhJtWeightBinLimBin[iContainer][iBin][jj]
					->Fill( jt, 1.0/jt * effCorrection );
				fhLogJtWeightBinLimBin[iContainer][iBin][jj]
					->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
				fhLogJtWeight2BinLimBin[iContainer][iBin][jj]
					->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
			}

		}



		//vOrtho.SetVect(jet->Vect().Orthogonal());
		vOrtho.SetVect(jet->Vect());
		vOrtho.SetE(jet->E());
		vOrtho.SetPhi(jet->Phi()+TMath::Pi()/2);

		//k=0; // count number of rotations of the orthogonal axis of a jet. 

		//Background jet (iBgJet) will be produced. This background jet is orthogonal to the iJet.  
		//If there is another jJet, then iBgJet will be consecutevely moved not to 
		//have jJet in the cone size. 
		int doBkg = 1;
		int counter = 0;
		if (Jets->GetEntries()>1){
			fhNumber[iContainer]->Fill(3.5);
			for (int j = 0; j<Jets->GetEntries(); j++){
				if (i == j) continue;
				AliJJet *jet2 = dynamic_cast<AliJJet*>( Jets->At(j) );
				if (!jet2) continue;

				/*if (k>15) {
				  fhNumber[iContainer]->Fill(5.5);
				  break;
				  }*/


				deltaR   = getDiffR(vOrtho.Phi(),jet2->Phi(),vOrtho.Eta(),jet2->Eta());
				if ( deltaR < 2*thisConeSize) {
					doBkg = 0;
					/*vOrtho.Rotate(TMath::Pi()/8, jet->Vect());
					  j=0;
					  k++;*/
					fhNumber[iContainer]->Fill(4.5);
					break;
				}

			}
			//fhKNumber[iContainer]->Fill(k);
		}

		if(doBkg){
			fhNumber[iContainer]->Fill(5.5); //TODO jet pt bins ->Fill(5.5+ iBin);
		}

		// Filling iBgJet,  Bgjt and Bgz
		double maxconpt = 0;
		AliJJet *jbg;
		if ( doBkg ){
			new (bgjets[iBgJet]) AliJJet(vOrtho.Px(),vOrtho.Py(), vOrtho.Pz(), vOrtho.E(), i ,0,0);
			//AliJJet * jbg = (AliJJet*) fJetBgListOfList[iContainer][iBgJet];
			jbg = (AliJJet*) fJetBgListOfList[iContainer][iBgJet];
			iBgJet++;

			pT = vOrtho.Pt(); 
			if (pT<(*fJetTriggPtBorders)[1]) continue;

			fhJetBgPt[iContainer]->Fill( pT );
			//bbfHistos[iContainer]->fhJetBgPtWeight->Fill( pT, 1./pT);
			iBin = GetBin(fJetTriggPtBorders, pT);
			if( iBin < 0 ) continue;
			fhJetBgPtBin[iContainer][iBin]->Fill( pT );
		}


		for (int icon = 0; icon<fTracks->GetEntries(); icon++){
			AliJBaseTrack *track = dynamic_cast<AliJBaseTrack*>(fTracks->At(icon));
			if (!track) continue;
			fhTrackEtaPhi[iContainer]->Fill(track->Eta(),track->Phi());
			pta = track->Pt();
			if (pta > maxconpt) maxconpt = pta;
			track->SetTrackEff( fEfficiency->GetCorrection( pta, 5, fcent) );
			effCorrection = 1.0/track->GetTrackEff();
			iptaBin = GetBin(fJetAssocPtBorders, pta);
			if( iptaBin < 0 ) continue;

			//Jet Cone Jt here
			deltaR   = getDiffR(jet->Phi(),track->Phi(),jet->Eta(),track->Eta());
			if ( deltaR < thisConeSize){ 
				fhJetConeTrkPt[iContainer]->Fill(pta,effCorrection);
				fhJetConeTrkPtBin[iContainer][iBin]->Fill(pta,effCorrection);
				fhJetConeTrkPtWeightBin[iContainer][iBin]->Fill(pta,effCorrection/pta);

				z = (track->Vect()*jet->Vect().Unit())/jet->P();
				fhJetConeZ[iContainer]->Fill( z , effCorrection);
				fhJetConeZBin[iContainer][iBin]->Fill( z , effCorrection);

				jt = (track->Vect()-z*jet->Vect()).Mag();
				fhJetConeJt[iContainer]->Fill( jt , effCorrection);
				fhJetConeJtBin[iContainer][iBin]->Fill( jt , effCorrection);
				fhJetConeJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt * effCorrection );
				fhJetConeLogJtWeightBin[iContainer][iBin]
					->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
				fhJetConeLogJtWeight2Bin[iContainer][iBin]
					->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );

				if (iptaBin < 0) continue;
				fhJetConeJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
					->Fill( jt, 1.0/jt * effCorrection );
				fhJetConeLogJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
					->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
				fhJetConeLogJtWithPtCutWeight2BinBin[iContainer][iBin][iptaBin]
					->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
			}


			if ( doBkg ){
				//Background jt
				deltaR   = getDiffR(vOrtho.Phi(),track->Phi(),vOrtho.Eta(),track->Eta());
				if ( deltaR < thisConeSize){ 
					counter++;
					jbg->AddConstituent(track);
					fhBgTrkPt[iContainer]->Fill(pta,effCorrection);
					fhBgTrkPtBin[iContainer][iBin]->Fill(pta,effCorrection);
					fhBgTrkPtWeightBin[iContainer][iBin]->Fill(pta,effCorrection/pta);
					if(moveJet){
						summedJet = track->GetLorentzVector() + vOrtho;
						z = (track->Vect()*summedJet.Vect().Unit())/summedJet.P();
						fhBgZ[iContainer]->Fill( z , effCorrection);
						fhBgZBin[iContainer][iBin]->Fill( z , effCorrection);

						jt = (track->Vect()-z*summedJet.Vect()).Mag();
						fhBgJt[iContainer]->Fill( jt , effCorrection);
						fhBgJtBin[iContainer][iBin]->Fill( jt , effCorrection);
						fhBgJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt * effCorrection );
						fhBgLogJtWeightBin[iContainer][iBin]
							->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
						fhBgLogJtWeight2Bin[iContainer][iBin]
							->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );

						if (iptaBin < 0) continue;
						fhBgJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
							->Fill( jt, 1.0/jt * effCorrection );
						fhBgLogJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
							->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
						fhBgLogJtWithPtCutWeight2BinBin[iContainer][iBin][iptaBin]
							->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
					}else{
						z = (track->Vect()*vOrtho.Vect().Unit())/vOrtho.P();
						fhBgZ[iContainer]->Fill( z , effCorrection);
						fhBgZBin[iContainer][iBin]->Fill( z , effCorrection);

						jt = (track->Vect()-z*vOrtho.Vect()).Mag();
						fhBgJt[iContainer]->Fill( jt , effCorrection);
						fhBgJtBin[iContainer][iBin]->Fill( jt , effCorrection);
						fhBgJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt * effCorrection );
						fhBgLogJtWeightBin[iContainer][iBin]
							->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
						fhBgLogJtWeight2Bin[iContainer][iBin]
							->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );

						if (iptaBin < 0) continue;
						fhBgJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
							->Fill( jt, 1.0/jt * effCorrection );
						fhBgLogJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
							->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
						fhBgLogJtWithPtCutWeight2BinBin[iContainer][iBin][iptaBin]
							->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
					}
				}
			}

		}
		if(doRndmBg){
			this->FillRandomBackground(jet->Pt(), jet->E(),  Jets, iContainer);
		}

		for (int ii = fJetConstPtLowLimits->GetNoElements(); ii >= 1 ; ii--)   {     
			if (maxconpt > (*fJetConstPtLowLimits)[ii]) {   
				jBin = ii-1;
				break;
			}
		}

		if(doBkg){
			fhBgTrkNumber[iContainer]->Fill(counter);
			fhBgTrkNumberBin[iContainer][iBin]->Fill(counter);
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
					fhBgLogJtWeight2BinLimBin[iContainer][iBin][jj]
						->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
				}
			}
		}
	}
	//doRndmBg = 1;
	/*if(doRndmBg){
	  this->FillRandomBackground(iContainer);
	  }*/
	/*for(int iBg = 0 ; iBg < doRndmBg ; iBg++){
	  this->FillRandomBackground(Jets, iContainer);
	  }*/
}

void AliJJetJtAnalysis::FillRandomBackground(TObjArray *Jets , int iContainer){
	int iBin, iptaBin=0;
	int jBin=0;
	double pT = 0;
	double conPtMax =0;
	iBin = GetBin(fJetTriggPtBorders,pT); // fill jetPt histos



	double z; double jt;
	double pta;
	//double Y , deltaY = 0;
	//double Phi, deltaPhi;
	//double deltaR= 0;
	//cout<<"histogram filling number of jets : "<<Jets->GetEntriesFast()<<endl;

	TLorentzVector  vRandomJet;
	TLorentzVector  randomTrack;

	int counter =0 ;
	int k = 0;
	double deltaR = -1;
	double deltaEta = -999;
	double deltaPhi = -999;
	double effCorrection = -1;
	double thisConeSize = fConeSizes[iContainer] ;
	double etaRndm = 0;
	double phiRndm = 0;
	etaRndm  = frandom->Uniform(-1.0+thisConeSize, 1.0-thisConeSize); //TODO Fix with realistic eta distribution
	phiRndm = kJPi * frandom->Uniform(-1, 1); //TODO realistic phi distribution
	vRandomJet.SetPtEtaPhiM(1,etaRndm,phiRndm,0);

	int countTrack = 1;
	fhNumber[iContainer]->Fill(6.5+iBin,Nrandom);
	for (int icon = 0; icon<fTracks->GetEntries(); icon++){
		countTrack = 1;
		AliJBaseTrack *track = dynamic_cast<AliJBaseTrack*>(fTracks->At(icon));
		if (!track) continue;
		if (Jets->GetEntries()>0){
			//fhNumber[iContainer]->Fill(3.5);
			for (int j = 0; j<Jets->GetEntries(); j++){
				AliJJet *jet = dynamic_cast<AliJJet*>( Jets->At(j) );
				if(jet->Pt() < 5) continue;
				if (!jet) continue;

				deltaR   = getDiffR(track->Phi(),jet->Phi(),track->Eta(),jet->Eta());
				if ( deltaR < thisConeSize) {
					countTrack = 0;
					break;
				}
			}
		}
		if(!countTrack) continue;
		//fhTrackEtaPhi[iContainer]->Fill(track->Eta(),track->Phi());
		pta = track->Pt();
		track->SetTrackEff( fEfficiency->GetCorrection( pta, 5, fcent) );
		effCorrection = 1.0/track->GetTrackEff();
		iptaBin = GetBin(fJetAssocPtBorders, pta);
		if( iptaBin < 0 ) continue;
		if(track->Eta()*track->Eta() < 1){
			//AliJBaseTrack *randomTrack = dynamic_cast<AliJBaseTrack*>(fTracks->At(icon));
			//etaRndm  = frandom->Uniform(-fmaxEtaRange, fmaxEtaRange);
			for(int irt = 0 ; irt < Nrandom ; irt++){
				etaRndm  = frandom->Uniform(-1, 1); //TODO Fix with realistic eta distribution
				phiRndm = kJPi * frandom->Uniform(-1, 1); //TODO realistic phi distribution
				randomTrack.SetVect(track->Vect());
				randomTrack.SetE(track->E());
				randomTrack.SetPtEtaPhiE(track->Pt(), etaRndm, phiRndm, track->E());
				//randomTrack.SetPtEtaPhiE(track->Pt(), track->Eta(), phiRndm, track->E());
				deltaR   = getDiffR(vRandomJet.Phi(),randomTrack.Phi(),vRandomJet.Eta(),randomTrack.Eta());
				if ( deltaR < thisConeSize){
					counter++;
					fhBgRndmTrkPt[iContainer]->Fill(pta,effCorrection);
					z = (randomTrack.Vect()*vRandomJet.Vect().Unit())/vRandomJet.P();
					fhBgRndmZ[iContainer]->Fill( z , effCorrection);

					jt = (randomTrack.Vect()-z*vRandomJet.Vect()).Mag();
					fhBgRndmJt[iContainer]->Fill( jt , 1.0/jt * effCorrection);
					fhBgRndmLogJt[iContainer]->Fill( jt , 1.0/jt/jt * effCorrection);

					if (iptaBin < 0) continue;
					fhBgRndmJtWithPtCutWeightBin[iContainer][iptaBin]
						->Fill( jt, 1.0/jt * effCorrection );
					fhBgRndmJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
						->Fill( jt, 1.0/jt * effCorrection );
					fhBgRndmLogJtWithPtCutWeight2Bin[iContainer][iptaBin]
						->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
					fhBgRndmLogJtWithPtCutWeight2BinBin[iContainer][iBin][iptaBin]
						->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
				}
			}
		}
	}
	fhBgRndmTrkNumber[iContainer]->Fill(counter/Nrandom);
}

void AliJJetJtAnalysis::FillRandomBackground(double jetpT, double jetE, TObjArray *Jets , int iContainer){
	int iBin, iptaBin=0;
	int jBin=0;
	iBin = GetBin(fJetTriggPtBorders,jetpT); // fill jetPt histos



	double z; double jt;
	double pta;
	//double Y , deltaY = 0;
	//double Phi, deltaPhi;
	//double deltaR= 0;
	//cout<<"histogram filling number of jets : "<<Jets->GetEntriesFast()<<endl;

	TLorentzVector  vRandomJet;
	TLorentzVector  randomTrack;
	TLorentzVector summedJet;

	int counter =0 ;
	int k = 0;
	double deltaR = -1;
	double deltaEta = -999;
	double deltaPhi = -999;
	double effCorrection = -1;
	double thisConeSize = fConeSizes[iContainer] ;
	double etaRndm = 0;
	double phiRndm = 0;
	etaRndm  = frandom->Uniform(-1.0+thisConeSize, 1.0-thisConeSize); //TODO Fix with realistic eta distribution
	phiRndm = kJPi * frandom->Uniform(-1, 1); //TODO realistic phi distribution
	vRandomJet.SetPtEtaPhiE(jetpT,etaRndm,phiRndm,jetE);

	int countTrack = 1;
	fhNumber[iContainer]->Fill(6.5+iBin,Nrandom);
	for (int icon = 0; icon<fTracks->GetEntries(); icon++){
		countTrack = 1;
		AliJBaseTrack *track = dynamic_cast<AliJBaseTrack*>(fTracks->At(icon));
		if (!track) continue;
		if (Jets->GetEntries()>0){
			//fhNumber[iContainer]->Fill(3.5);
			for (int j = 0; j<Jets->GetEntries(); j++){
				AliJJet *jet = dynamic_cast<AliJJet*>( Jets->At(j) );
				if(jet->Pt() < 5) continue;
				if (!jet) continue;

				deltaR   = getDiffR(track->Phi(),jet->Phi(),track->Eta(),jet->Eta());
				if ( deltaR < thisConeSize) {
					countTrack = 0;
					break;
				}
			}
		}
		if(!countTrack) continue;
		//fhTrackEtaPhi[iContainer]->Fill(track->Eta(),track->Phi());
		pta = track->Pt();
		track->SetTrackEff( fEfficiency->GetCorrection( pta, 5, fcent) );
		effCorrection = 1.0/track->GetTrackEff();
		iptaBin = GetBin(fJetAssocPtBorders, pta);
		if( iptaBin < 0 ) continue;
		if(track->Eta()*track->Eta() < 1){
			//AliJBaseTrack *randomTrack = dynamic_cast<AliJBaseTrack*>(fTracks->At(icon));
			//etaRndm  = frandom->Uniform(-fmaxEtaRange, fmaxEtaRange);
			for(int irt = 0 ; irt < Nrandom ; irt++){
				etaRndm  = frandom->Uniform(-1, 1); //TODO Fix with realistic eta distribution
				phiRndm = kJPi * frandom->Uniform(-1, 1); //TODO realistic phi distribution
				randomTrack.SetVect(track->Vect());
				randomTrack.SetE(track->E());
				randomTrack.SetPtEtaPhiE(track->Pt(), etaRndm, phiRndm, track->E());
				deltaR   = getDiffR(vRandomJet.Phi(),randomTrack.Phi(),vRandomJet.Eta(),randomTrack.Eta());
				if ( deltaR < thisConeSize){
					counter++;
					if(moveJet){
						summedJet = randomTrack + vRandomJet;
						fhBgRndmTrkPt[iContainer]->Fill(pta,effCorrection);
						z = (randomTrack.Vect()*summedJet.Vect().Unit())/summedJet.P();
						fhBgRndmZ[iContainer]->Fill( z , effCorrection);

						jt = (randomTrack.Vect()-z*summedJet.Vect()).Mag();
						fhBgRndmJt[iContainer]->Fill( jt , 1.0/jt * effCorrection);
						fhBgRndmLogJt[iContainer]->Fill( jt , 1.0/jt/jt * effCorrection);

						if (iptaBin < 0) continue;
						fhBgRndmJtWithPtCutWeightBin[iContainer][iptaBin]
							->Fill( jt, 1.0/jt * effCorrection );
						fhBgRndmJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
							->Fill( jt, 1.0/jt * effCorrection );
						fhBgRndmLogJtWithPtCutWeight2Bin[iContainer][iptaBin]
							->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
						fhBgRndmLogJtWithPtCutWeight2BinBin[iContainer][iBin][iptaBin]
							->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
					}
					else{
						//randomTrack.SetPtEtaPhiE(track->Pt(), track->Eta(), phiRndm, track->E());
						fhBgRndmTrkPt[iContainer]->Fill(pta,effCorrection);
						z = (randomTrack.Vect()*vRandomJet.Vect().Unit())/vRandomJet.P();
						fhBgRndmZ[iContainer]->Fill( z , effCorrection);

						jt = (randomTrack.Vect()-z*vRandomJet.Vect()).Mag();
						fhBgRndmJt[iContainer]->Fill( jt , 1.0/jt * effCorrection);
						fhBgRndmLogJt[iContainer]->Fill( jt , 1.0/jt/jt * effCorrection);

						if (iptaBin < 0) continue;
						fhBgRndmJtWithPtCutWeightBin[iContainer][iptaBin]
							->Fill( jt, 1.0/jt * effCorrection );
						fhBgRndmJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
							->Fill( jt, 1.0/jt * effCorrection );
						fhBgRndmLogJtWithPtCutWeight2Bin[iContainer][iptaBin]
							->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
						fhBgRndmLogJtWithPtCutWeight2BinBin[iContainer][iBin][iptaBin]
							->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
					}
				}
			}
		}
	}
	fhBgRndmTrkNumber[iContainer]->Fill(counter/Nrandom);
}
// new Bg jt will be filled with a new cone size nR at histograms in iHist
// cone size of Jets always shold be greather than nR 
// to calcualte new bg jt with smaller cone size from Jets constituents with larger cone size
void AliJJetJtAnalysis::FillBgJtWithSmallerR(const TClonesArray &Jets, double nR, int iHist){
	double iBin = -1, iptaBin = -1;
	double pT=-1, z=-1,jt=-1,  pta=-1;
	double effCorrection = -1;
	for (int i = 0; i<Jets.GetEntries(); i++){
		AliJJet *jet = dynamic_cast<AliJJet*>( Jets.At(i) );
		if (!jet) continue;
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
			if (jet->DeltaR(*con)>nR) continue;
			jt = (con->Vect()-z*jet->Vect()).Mag();
			fhBgJtWithPtCutWeightBinBinSmallerR[iHist][iBin][iptaBin]
				->Fill( jt, 1.0/jt * effCorrection );
			fhBgLogJtWithPtCutWeightBinBinSmallerR[iHist][iBin][iptaBin]
				->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
			fhBgLogJtWithPtCutWeight2BinBinSmallerR[iHist][iBin][iptaBin]
				->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
		}


	}


}

// bg jt is newly calculated by a jet axis of fJetBgListOfList[ian]
// with bg constituents from fJetBgListOfList[iao]
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
		if (!jo) continue;
		for (int in = 0; in<an.GetEntries(); in++){
			AliJJet *jn = dynamic_cast<AliJJet*>( an.At(in) );
			if (!jn) continue;
			if (jo->DeltaR(*jn) > fConeSizes[ian]) continue;
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
				fhBgLogJtWithPtCutWeight2BinBinDiffR[iHist][iBin][iptaBin]
					->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
			}
		}
	}

}

//Phi1 and Phi2 between 0 and 2 pi
Double_t AliJJetJtAnalysis::getDiffR(double phi1, double phi2, double eta1, double eta2){
	Double_t diffPhi = TMath::Abs(phi1-phi2);
	if(diffPhi > TMath::Pi()){
		diffPhi = 2*TMath::Pi() - diffPhi;
	}
	return TMath::Sqrt(TMath::Power(diffPhi,2)+TMath::Power(eta1-eta2,2));
}


