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

// Class for jT analysis in jets

//===========================================================
// Dummy comment, should be replaced by a real one
// comment
// comment
// jT Analysis class by Beomkyu Kim, Dongjo Kim and Tomas Snellman
//===========================================================

#include <TRandom.h>
#include <TMath.h>
#include <TRegexp.h>
#include <TVector.h>
#include "AliJJet.h"
#include "AliEmcalJet.h"
#include "AliJEfficiency.h"
#include "AliJJetJtAnalysis.h"
#include "AliJHistManager.h"
#include "AliJMCTrack.h"
#include "TClonesArray.h"
#include "AliAODMCParticle.h"

AliJJetJtAnalysis::AliJJetJtAnalysis():
  fInputList(NULL)
  , fJetList(NULL)
  , fJetListOfList() // Jet finders container
  , fMCJetList(NULL)
  , fMCJetListOfList() // Jet finders container
  //, fJetBgList(NULL) 
  , fJetBgListOfList() // Bg jet list of list container
  , fpythiaJets("AliJJet",10)
  , fJetEtaCut() // Jet's eta cut
  , fLeadingJets(0)
  , frandom(0x0)
  , fJetTriggPtBorders(NULL)
  , fJetConstPtLowLimits(NULL)
  , fJetAssocPtBorders(NULL)
  , fJetMultBorders(NULL)
  , fDeltaRBorders(NULL)
  , fJetLeadPtBorders(NULL)
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
  , fMCTracks(NULL) //list of tracks
  , fTrackJt(NULL)    //Store jT of tracks
  , fConstJt(NULL)    //Store jT of jet constituents
  , fTrackPt(NULL)    //Store pT of tracks
  , fConstPt(NULL)    //Store pT of jet constituents
  , fConstLabels(NULL)    //Store labels of jet constituents
  , fJetPt(NULL)    //Store pT of corresponding jets
  , fDiJetMjj(NULL) //Store di-jet invariant mass for each jet finder
  , fTrackFound(NULL) //Keep track of matched tracks
  , fConstFound(NULL) //Keep track of matched tracks
  , fBin2(NULL)    //Store iBin2 of corresponding tracks 
  , fBin3(NULL)    //Store iBin3 of corresponding tracks
  , fpta(NULL)    //Store pTa bin of corresponding tracks
  , fptt(NULL)    //Store pTt bin of corresponding tracks
  , Nrandom(1)
  , moveJet(1)
  , fDoMC(0)
  , fDoLog(0)
  , fHMG(NULL)
  , fHMGMC(NULL)
  , fJetFinderBin() //histmanager axis
  , fJetTriggerBin() //histmanager axis
  , fTrkPtBin() //histmanager axis
  , fTrkLimPtBin() //histmanager axis
  , fJetLeadPtBin() //histmanager axis
  , fJetMultBin() //histmanager axis
  , fdRBin() //histmanager axis
  , fiHist() //histmanager axis
  , fJetFinderBinMC() //histmanager axis
  , fJetTriggerBinMC() //histmanager axis
  , fTrkPtBinMC() //histmanager axis
  , fTrkLimPtBinMC() //histmanager axis
  , fJetLeadPtBinMC() //histmanager axis
  , fJetMultBinMC() //histmanager axis
  , fdRBinMC() //histmanager axis
  , fiHistMC() //histmanager axis
  , fBgTypeBin() //histmanager axis
  , fhNumber() // number of jets
  , fhJetPt() //Jet pt distribution
  , fhJetPtBin() //Jet pt as a fn of jet pt 
  , fhJetPtTrackCutBin()
  , fhJetPtMultiplicityCutBin()
  , fhJetPtWeight() //Jet pt distribution
  , fhJetPtWeightBin() //Jet pt as a fn of jet pt 
  , fhJetMultiplicityBin() //Number of tracks in jet
  , fhJetMultiplicityBin2() //Number of tracks with jT filled in jet
  , fhJetConeMultiplicityBin() //Number of tracks in jet cone
  , fhDiJetM() // Di-jet invariant mass
  , fhZ() // z dist
  , fhZBin() // z dist as fn of jet pt
  , fhRBin() // R dist as fn of jet pt
  , fhJt() // jt dist
  , fhJtBin() // jt as fn of jet pt
  , fhJtWeightBin() // jt in P.S as fn of jet pt
  , fhJtWeightBinTest() // jt in P.S as fn of jet pt
  , fhJtWeightBinTest2() // jt in P.S as fn of jet pt
  , fhLogJtWeightBin() //log jt in P.S as fn of jet pt
  , fhLogJtWeight2Bin() //log jt in P.S as fn of jet pt
  , fhJtWithPtCutWeightBinBin() // jt in P.S as fns of jet pt, const pt
  , fhJtBinLimBin() //jt with maximum const pt bin
  , fhJtWeightBinLimBin()
  , fhLeadingJt() // jt dist for leading track
  , fhLeadingJtBin() // jt as fn of jet pt for leading track
  , fhLeadingJtWeightBin() // jt in P.S as fn of jet pt for leading track
  , fhLeadingJtWithPtCutWeightBinBin() // jt in P.S as fns of jet pt, const pt for leading track
  , fhLogJtWithPtCutWeightBinBin()
  , fhLogJtWithPtCutWeight2BinBin()
  , fhLogJtWeightBinLimBin()
  , fhLogJtWeight2BinLimBin()
  , fhEventJtBin()
  , fhEventJtWeightBin()
  , fhEventJtWithPtCutWeightBinBin()
  , fhJetConeTrkPt()
  , fhJetConeTrkPtBin()
  , fhJetConeTrkPtWeightBin()
  , fhJetConeZ()
  , fhJetConeZBin()
  , fhJetConeJt()
  , fhJetConeJtBin()
  , fhJetConeJtWeightBin()
  , fhJetConeJtWeightWithTrackCutBinBin()
  , fhJetConeJtWeightWithMultiplicityCutBinBin()
  , fhJetConeLogJtWeightBin()
  , fhJetConeLogJtWeight2Bin()
  , fhJetConeJtWithPtCutWeightBinBin()
  , fhJetConeLogJtWithPtCutWeightBinBin()
  , fhJetConeLogJtWithPtCutWeight2BinBin()
  , fhJtUnfBg()
  , fhJtBinUnfBg()
  , fhJtWeightBinUnfBg()
  , fhLogJtWeightBinUnfBg()
  , fhLogJtWeight2BinUnfBg() 
  , fhJetConeJtUnfBg()
  , fhJetConeJtBinUnfBg()
  , fhJetConeJtWeightBinUnfBg()
  , fhJetConeLogJtWeightBinUnfBg()
  , fhJetConeLogJtWeight2BinUnfBg()
  , fhJetBgPt()
  , fhJetBgPtBin()
  , fhBgZ() //Bg Z 
  , fhBgZBin() //Bg Z as fn of jet pt
  , fhBgRBin() //Bg R as fn of jet pt
  , fhBgRndmRBin() //Rndm Bg R as fn of jet pt
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
  , fhLeadingTrkPtBin()
  , fhBgTrkPt()
  , fhBgTrkPtBin()
  , fhBgTrkPtWeightBin()
  , fhBgTrkNumber()
  , fhBgTrkNumberBin()
  , fhBgRndmTrkPt()
  , fhBgRndmZ()
  , fhBgRndmTrkNumber()
  , fhBgRndmJt()
  , fhBgRndmJtBin()
  , fhBgRndmJtWeight()
  , fhBgRndmJtWeightBin()
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

  , fhTrackJtCorrBin()
  //, fhTrackJtCorr2D()
  , fhConstJtCorrBin()
  , fhTrackJtCorrBinTest()
  , fhTrackJtCorrBinTest2()
  , fhJetPtCorr() //Correlation between true and reco mc jets
  , fhJetPtCorr2() 
  , fhJetPtCorrCoarse() //Correlation between true and reco mc jets
  , fhJetdR() //Jet distance distribution
  , fhJetdPt() //Jet momentum difference distribution
  , fhTrackMatchSuccess()
  , fhConstMatchSuccess()
  , fhTrackPtCorr() //Correlation between true and reco mc track pt
  , fhConstPtCorr() //Correlation between true and reco mc track pt for jet constituents
  , fhDiJetMjjCorr() //Correlation between true and reco mc di-jet invariant mass
  , fhTrackJtCorrBinPythia()
  , fhJetPtCorrPythia() //Correlation between true and reco mc jets
  , fhJetPtCorrPythiaCoarse() //Correlation between true and reco mc jets
  , fhTrackMatchSuccessPythia()
  , fhJetPtCorr2Pythia() 
  , fhTrackPtCorrPythia() //Correlation between true and reco mc track pt
  , fhJetdRPythia() //Jet distance distribution
  , fhJetdPtPythia() //Jet momentum difference distribution
  , fhJetPtPythia()
  , fhJetPtBinPythia()
  , fhJetMultiplicityBinPythia()
  , fhZPythia() // z dist
  , fhZBinPythia() // z dist as fn of jet pt
  , fhZBinPythiaRCut() // z dist as fn of jet pt where R < 0.4
  , fhRBinPythia() // R dist as fn of jet pt
  , fhJtPythia()
  , fhJtBinPythia()
  , fhJtBinPythiaRCut()
  , fhJtWeightBinPythia()
  , fhLogJtWeightBinPythia()
  , fhLogJtWeight2BinPythia()
  , fhJtWithPtCutWeightBinBinPythia()
  , fhLogJtWithPtCutWeightBinBinPythia()
  , fhLogJtWithPtCutWeight2BinBinPythia()
  , fhBgJtPythia()
  , fhBgJtPythiaTrackType()
  , fhBgJtBinPythia()
  , fhBgJtBinPythiaTypeBin()
  , fhBgJtWeightBinPythia()
  , fhBgLogJtWeightBinPythia()
  , fhBgLogJtWeight2BinPythia()
  , fhBgJtWithPtCutWeightBinBinPythia()
  , fhBgLogJtWithPtCutWeightBinBinPythia()
  , fhBgLogJtWithPtCutWeight2BinBinPythia()

{
  frandom = new TRandom3(); // frandom generator for jt flow UE
  frandom->SetSeed(0);
}

AliJJetJtAnalysis::AliJJetJtAnalysis( AliJCard * card ):
  fInputList(NULL)
  , fJetList(NULL)
  , fJetListOfList()
  , fMCJetList(NULL)
  , fMCJetListOfList() // Jet finders container
  , fJetBgListOfList()
  , fpythiaJets("AliJJet",10)
  , fJetEtaCut() 
  , fLeadingJets(0)
  , frandom(0x0)
  , fJetTriggPtBorders(NULL)
  , fJetConstPtLowLimits(NULL)
  , fJetAssocPtBorders(NULL)
  , fJetMultBorders(NULL)
  , fDeltaRBorders(NULL)
  , fJetLeadPtBorders(NULL)
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
  , fMCTracks(NULL)
  , fTrackJt(NULL)
  , fConstJt(NULL)
  , fTrackPt(NULL)
  , fConstPt(NULL)
  , fConstLabels(NULL)
  , fJetPt(NULL)
  , fDiJetMjj(NULL)
  , fTrackFound(NULL)
  , fConstFound(NULL)
  , Nrandom(1)
  , moveJet(1)
  , fDoMC(0)
  , fDoLog(0)
  , fHMG(NULL)
  , fHMGMC(NULL)
  , fJetFinderBin()
  , fJetTriggerBin()
  , fTrkPtBin()
  , fTrkLimPtBin()
  , fJetLeadPtBin() //histmanager axis
  , fJetMultBin() //histmanager axis
  , fdRBin()
  , fiHist()
  , fJetFinderBinMC() //histmanager axis
  , fJetTriggerBinMC() //histmanager axis
  , fTrkPtBinMC() //histmanager axis
  , fTrkLimPtBinMC() //histmanager axis
  , fJetLeadPtBinMC() //histmanager axsi
  , fJetMultBinMC() //histmanager axis
  , fdRBinMC() //histmanager axis
  , fiHistMC() //histmanager axis
  , fBgTypeBin()
  , fhNumber()
  , fhJetPt()
  , fhJetPtBin()
  , fhJetPtTrackCutBin()
  , fhJetPtMultiplicityCutBin()
  , fhJetPtWeight() //Jet pt distribution
  , fhJetPtWeightBin() //Jet pt as a fn of jet pt 
  , fhJetMultiplicityBin() //Number of tracks in jet
  , fhJetMultiplicityBin2() //Number of tracks with jT filled in jet
  , fhJetConeMultiplicityBin() //Number of tracks in jet cone
  , fhDiJetM() // Di-jet invariant mass
  , fhZ()
  , fhZBin()
  , fhJt()
  , fhJtBin()
  , fhJtWeightBin()
  , fhJtWeightBinTest()
  , fhJtWeightBinTest2()
  , fhLogJtWeightBin()
  , fhLogJtWeight2Bin()
  , fhJtWithPtCutWeightBinBin()
  , fhJtBinLimBin()
  , fhJtWeightBinLimBin()
  , fhLeadingJt() // jt dist for leading track
  , fhLeadingJtBin() // jt as fn of jet pt for leading track
  , fhLeadingJtWeightBin() // jt in P.S as fn of jet pt for leading track
  , fhLeadingJtWithPtCutWeightBinBin() // jt in P.S as fns of jet pt, const pt for leading track
  , fhLogJtWithPtCutWeightBinBin()
  , fhLogJtWithPtCutWeight2BinBin()
  , fhLogJtWeightBinLimBin()
  , fhLogJtWeight2BinLimBin()
  , fhEventJtBin()
  , fhEventJtWeightBin()
  , fhEventJtWithPtCutWeightBinBin()
  , fhJetConeTrkPt()
  , fhJetConeTrkPtBin()
  , fhJetConeTrkPtWeightBin()
  , fhJetConeZ()
  , fhJetConeZBin()
  , fhJetConeJt()
  , fhJetConeJtBin()
  , fhJetConeJtWeightBin()
  , fhJetConeJtWeightWithTrackCutBinBin()
  , fhJetConeJtWeightWithMultiplicityCutBinBin()
  , fhJetConeLogJtWeightBin()
  , fhJetConeLogJtWeight2Bin()
  , fhJetConeJtWithPtCutWeightBinBin()
  , fhJetConeLogJtWithPtCutWeightBinBin()
  , fhJetConeLogJtWithPtCutWeight2BinBin()
  , fhJtUnfBg()
  , fhJtBinUnfBg()
  , fhJtWeightBinUnfBg()
  , fhLogJtWeightBinUnfBg()
  , fhLogJtWeight2BinUnfBg()
  , fhJetConeJtUnfBg()
  , fhJetConeJtBinUnfBg()
  , fhJetConeJtWeightBinUnfBg()
  , fhJetConeLogJtWeightBinUnfBg()
  , fhJetConeLogJtWeight2BinUnfBg()
  , fhJetBgPt()
  , fhJetBgPtBin()
  , fhBgZ()
  , fhBgZBin()
  , fhBgRBin() //Bg R as fn of jet pt
  , fhBgRndmRBin() //Rndm Bg R as fn of jet pt
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
  , fhTrkPt()
  , fhTrkPtBin()
  , fhTrkPtWeightBin()
  , fhLeadingTrkPtBin()
  , fhBgTrkPt()
  , fhBgTrkPtBin()
  , fhBgTrkPtWeightBin()
  , fhBgTrkNumber()
  , fhBgTrkNumberBin()
  , fhBgRndmTrkPt()
  , fhBgRndmZ()
  , fhBgRndmTrkNumber()
  , fhBgRndmJt()
  , fhBgRndmJtBin()
  , fhBgRndmJtWeight()
  , fhBgRndmJtWeightBin()
  , fhBgRndmLogJt()
  , fhBgRndmJtWithPtCutWeightBin()
  , fhBgRndmLogJtWithPtCutWeight2Bin()
  , fhBgRndmJtWithPtCutWeightBinBin()
  , fhBgRndmLogJtWithPtCutWeight2BinBin()
  , fhdeltaE()
  , fhdeltaN()
  , fhFullJetEChJetBin()
  , fhFullChdRChJetBin()
  , fh2DFullEvsChEdN0()
  , fh2DFullEvsChEdNnot0()
  , fhJetEtaPhi()
  , fhTrackEtaPhi()
  // Bg jt as fns of jet pt, track pt
  // log of Bg jt in PS as fns of jet pt, track pt
  // Bg jt in PS as fns of jet pt, track pt
  // if R = 0.6, all bg jts are filled for all R<0.6
  // Bg jt in PS as fns of jet pt, track pt
  // from new jet axis with present tracks
  // Bg jt as a fn of jet pt 
  // with maximum constituent's pT bin in a jet

  // dE = E of charged jet - E of shared tracks of
  // charged and full jets
  // dN = N of charged jet - N of shared tracks of
  // charged and full jets
  // Fulljet E as a fn of charged jet E 
  // dR between full and charged jet 
  // as a fn of charged jet E
  // x-axis: fulljet E, y: ch jet E when dN == 0
  // x-axis: fulljet E, y: ch jet E when dN != 0
  //Jet eta and Phi 2D dist for each jet finder
  , fhTrackJtCorrBin()
  //, fhTrackJtCorr2D()
  , fhConstJtCorrBin()
  , fhTrackJtCorrBinTest()
  , fhTrackJtCorrBinTest2()
  , fhJetPtCorr()
  , fhJetPtCorr2() 
  , fhJetPtCorrCoarse()
  , fhJetdR() 
  , fhJetdPt() 
  , fhTrackMatchSuccess()
  , fhConstMatchSuccess()
  , fhTrackPtCorr() 
  , fhConstPtCorr() 
  , fhDiJetMjjCorr() //Correlation between true and reco mc di-jet invariant mass
  , fhTrackJtCorrBinPythia()
  , fhJetPtCorrPythia()
  , fhJetPtCorrPythiaCoarse()
  , fhTrackMatchSuccessPythia()
  , fhJetPtCorr2Pythia() 
  , fhTrackPtCorrPythia() 
  , fhJetdRPythia() 
  , fhJetdPtPythia() 
  , fhJetPtPythia()
  , fhJetPtBinPythia()
  , fhJetMultiplicityBinPythia()
  , fhZPythia() // z dist
  , fhZBinPythia() // z dist as fn of jet pt
  , fhZBinPythiaRCut() // z dist as fn of jet pt with R < 0.4
  , fhJtPythia()
  , fhJtBinPythia()
  , fhJtBinPythiaRCut()
  , fhJtWeightBinPythia()
  , fhLogJtWeightBinPythia()
  , fhLogJtWeight2BinPythia()
  , fhJtWithPtCutWeightBinBinPythia()
  , fhLogJtWithPtCutWeightBinBinPythia()
  , fhLogJtWithPtCutWeight2BinBinPythia()
  , fhBgJtPythia()
  , fhBgJtPythiaTrackType()
  , fhBgJtBinPythia()
  , fhBgJtBinPythiaTypeBin()
  , fhBgJtWeightBinPythia()
  , fhBgLogJtWeightBinPythia()
  , fhBgLogJtWeight2BinPythia()
  , fhBgJtWithPtCutWeightBinBinPythia()
  , fhBgLogJtWithPtCutWeightBinBinPythia()
  , fhBgLogJtWithPtCutWeight2BinBinPythia()
{

  frandom = new TRandom3(); // frandom generator for jt flow UE
  frandom->SetSeed(0);
}

AliJJetJtAnalysis::AliJJetJtAnalysis(const AliJJetJtAnalysis& ap) :
  fInputList(ap.fInputList)
  , fJetList(ap.fJetList)
  , fJetListOfList(ap.fJetListOfList)
  , fMCJetList(ap.fMCJetList)
  , fMCJetListOfList(ap.fMCJetListOfList)
  //, fJetBgList(ap.fJetBgList)
  , fJetBgListOfList(ap.fJetBgListOfList)
  , fpythiaJets(ap.fpythiaJets)
  , fJetEtaCut(ap.fJetEtaCut) // Jet's eta cut
  , fLeadingJets(ap.fLeadingJets)
  , frandom(0x0)
  , fJetTriggPtBorders(ap.fJetTriggPtBorders)
  , fJetConstPtLowLimits(ap.fJetConstPtLowLimits)
  , fJetAssocPtBorders(ap.fJetAssocPtBorders)
  , fJetMultBorders(ap.fJetMultBorders)
  , fDeltaRBorders(ap.fDeltaRBorders)
  , fJetLeadPtBorders(ap.fJetLeadPtBorders)
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
  , fMCTracks(ap.fMCTracks)
  , fTrackPt(ap.fTrackPt)
  , fConstPt(ap.fConstPt)
  , fConstLabels(ap.fConstLabels)
  , Nrandom(ap.Nrandom)
  , moveJet(ap.moveJet)
  , fDoMC(ap.fDoMC)
  , fDoLog(ap.fDoLog)
  , fHMG(ap.fHMG)
  , fHMGMC(ap.fHMGMC)
  , fJetFinderBin(ap.fJetFinderBin)
  , fJetTriggerBin(ap.fJetTriggerBin)
  , fTrkPtBin(ap.fTrkPtBin)
  , fTrkLimPtBin(ap.fTrkLimPtBin)
  , fJetLeadPtBin(ap.fJetLeadPtBin)
  , fJetMultBin(ap.fJetMultBin)
  , fdRBin(ap.fdRBin)
  , fiHist(ap.fiHist)
  , fJetFinderBinMC(ap.fJetFinderBinMC) 
  , fJetTriggerBinMC(ap.fJetTriggerBinMC)
  , fTrkPtBinMC(ap.fTrkPtBinMC) 
  , fTrkLimPtBinMC(ap.fTrkLimPtBinMC)
  , fJetLeadPtBinMC(ap.fJetLeadPtBinMC)
  , fJetMultBinMC(ap.fJetMultBinMC)
  , fdRBinMC(ap.fdRBinMC) 
  , fiHistMC(ap.fiHistMC) 
  , fBgTypeBin(ap.fBgTypeBin) 
  , fhNumber(ap.fhNumber)
  , fhJetPt(ap.fhJetPt)
  , fhJetPtBin(ap.fhJetPtBin)
  , fhJetPtTrackCutBin(ap.fhJetPtTrackCutBin)
  , fhJetPtMultiplicityCutBin(ap.fhJetPtMultiplicityCutBin)
  , fhJetPtWeight(ap.fhJetPtWeight)
  , fhJetPtWeightBin(ap.fhJetPtWeightBin)
  , fhJetMultiplicityBin(ap.fhJetMultiplicityBin) //Number of tracks in jet
  , fhJetMultiplicityBin2(ap.fhJetMultiplicityBin2) 
  , fhJetConeMultiplicityBin(ap.fhJetConeMultiplicityBin) 
  , fhDiJetM(ap.fhDiJetM)
  , fhZ(ap.fhZ)
  , fhZBin(ap.fhZBin)
  , fhJt(ap.fhJt)
  , fhJtBin(ap.fhJtBin)
  , fhJtWeightBin(ap.fhJtWeightBin)
  , fhJtWeightBinTest(ap.fhJtWeightBinTest)
  , fhJtWeightBinTest2(ap.fhJtWeightBinTest2)
  , fhLogJtWeightBin(ap.fhLogJtWeightBin)
  , fhLogJtWeight2Bin(ap.fhLogJtWeight2Bin)
  , fhJtWithPtCutWeightBinBin(ap.fhJtWithPtCutWeightBinBin)
  , fhJtBinLimBin(ap.fhJtBinLimBin)
  , fhJtWeightBinLimBin(ap.fhJtWeightBinLimBin)
  , fhLeadingJt(ap.fhLeadingJt)
  , fhLeadingJtBin(ap.fhLeadingJtBin)
  , fhLeadingJtWeightBin(ap.fhLeadingJtWeightBin)
  , fhLeadingJtWithPtCutWeightBinBin(ap.fhLeadingJtWithPtCutWeightBinBin)
  , fhLogJtWithPtCutWeightBinBin(ap.fhLogJtWithPtCutWeightBinBin)
  , fhLogJtWithPtCutWeight2BinBin(ap.fhLogJtWithPtCutWeight2BinBin)
  , fhLogJtWeightBinLimBin(ap.fhLogJtWeightBinLimBin)
  , fhLogJtWeight2BinLimBin(ap.fhLogJtWeight2BinLimBin)
  , fhEventJtBin(ap.fhEventJtBin)
  , fhEventJtWeightBin(ap.fhEventJtWeightBin)
  , fhEventJtWithPtCutWeightBinBin(ap.fhEventJtWithPtCutWeightBinBin)
  , fhJetConeTrkPt(ap.fhJetConeTrkPt)
  , fhJetConeTrkPtBin(ap.fhJetConeTrkPtBin)
  , fhJetConeTrkPtWeightBin(ap.fhJetConeTrkPtWeightBin)
  , fhJetConeZ(ap.fhJetConeZ)
  , fhJetConeZBin(ap.fhJetConeZBin)
  , fhJetConeJt(ap.fhJetConeJt)
  , fhJetConeJtBin(ap.fhJetConeJtBin)
  , fhJetConeJtWeightBin(ap.fhJetConeJtWeightBin)
  , fhJetConeJtWeightWithTrackCutBinBin(ap.fhJetConeJtWeightWithTrackCutBinBin)
  , fhJetConeJtWeightWithMultiplicityCutBinBin(ap.fhJetConeJtWeightWithMultiplicityCutBinBin)
  , fhJetConeLogJtWeightBin(ap.fhJetConeLogJtWeightBin)
  , fhJetConeLogJtWeight2Bin(ap.fhJetConeLogJtWeight2Bin)
  , fhJetConeJtWithPtCutWeightBinBin(ap.fhJetConeJtWithPtCutWeightBinBin)
  , fhJetConeLogJtWithPtCutWeightBinBin(ap.fhJetConeLogJtWithPtCutWeightBinBin)
  , fhJetConeLogJtWithPtCutWeight2BinBin(ap.fhJetConeLogJtWithPtCutWeight2BinBin)
  , fhJtUnfBg(ap.fhJtUnfBg)
  , fhJtBinUnfBg(ap.fhJtBinUnfBg)
  , fhJtWeightBinUnfBg(ap.fhJtWeightBinUnfBg)
  , fhLogJtWeightBinUnfBg(ap.fhLogJtWeightBinUnfBg)
  , fhLogJtWeight2BinUnfBg(ap.fhLogJtWeight2BinUnfBg)
  , fhJetConeJtUnfBg(ap.fhJetConeJtUnfBg)
  , fhJetConeJtBinUnfBg(ap.fhJetConeJtBinUnfBg)
  , fhJetConeJtWeightBinUnfBg(ap.fhJetConeJtWeightBinUnfBg)
  , fhJetConeLogJtWeightBinUnfBg(ap.fhJetConeLogJtWeightBinUnfBg)
  , fhJetConeLogJtWeight2BinUnfBg(ap.fhJetConeLogJtWeight2BinUnfBg)
  , fhJetBgPt(ap.fhJetBgPt)
  , fhJetBgPtBin(ap.fhJetBgPtBin)
  , fhBgZ(ap.fhBgZ)
  , fhBgZBin(ap.fhBgZBin)
  , fhBgRBin(ap.fhBgRBin) 
  , fhBgRndmRBin(ap.fhBgRndmRBin)
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
  , fhLeadingTrkPtBin(ap.fhLeadingTrkPtBin)
  , fhBgTrkPt(ap.fhBgTrkPt)
  , fhBgTrkPtBin(ap.fhBgTrkPtBin)
  , fhBgTrkPtWeightBin(ap.fhBgTrkPtWeightBin)
  , fhBgTrkNumber(ap.fhBgTrkNumber)
  , fhBgTrkNumberBin(ap.fhBgTrkNumberBin)
  , fhBgRndmTrkPt(ap.fhBgRndmTrkPt)
  , fhBgRndmZ(ap.fhBgRndmZ)
  , fhBgRndmTrkNumber(ap.fhBgRndmTrkNumber)
  , fhBgRndmJt(ap.fhBgRndmJt)
  , fhBgRndmJtBin(ap.fhBgRndmJtBin)
  , fhBgRndmJtWeight(ap.fhBgRndmJtWeight)
  , fhBgRndmJtWeightBin(ap.fhBgRndmJtWeightBin)
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
  , fhTrackJtCorrBin(ap.fhTrackJtCorrBin)
  //, fhTrackJtCorr2D(ap.fhTrackJtCorr2D)
  , fhConstJtCorrBin(ap.fhConstJtCorrBin)
  , fhTrackJtCorrBinTest(ap.fhTrackJtCorrBinTest)
  , fhTrackJtCorrBinTest2(ap.fhTrackJtCorrBinTest2)
  , fhJetPtCorr(ap.fhJetPtCorr)
  , fhJetPtCorr2(ap.fhJetPtCorr2) 
  , fhJetPtCorrCoarse(ap.fhJetPtCorrCoarse)
  , fhJetdR(ap.fhJetdR)
  , fhJetdPt(ap.fhJetdPt)
  , fhTrackMatchSuccess(ap.fhTrackMatchSuccess)
  , fhConstMatchSuccess(ap.fhConstMatchSuccess)
  , fhTrackPtCorr(ap.fhTrackPtCorr) 
  , fhConstPtCorr(ap.fhConstPtCorr)
  , fhDiJetMjjCorr(ap.fhDiJetMjjCorr) 
  , fhTrackJtCorrBinPythia(ap.fhTrackJtCorrBinPythia)
  , fhJetPtCorrPythia(ap.fhJetPtCorrPythia)
  , fhJetPtCorrPythiaCoarse(ap.fhJetPtCorrPythiaCoarse)
  , fhTrackMatchSuccessPythia(ap.fhTrackMatchSuccessPythia)
  , fhJetPtCorr2Pythia(ap.fhJetPtCorr2Pythia) 
  , fhTrackPtCorrPythia(ap.fhTrackPtCorrPythia) 
  , fhJetdRPythia(ap.fhJetdRPythia)
  , fhJetdPtPythia(ap.fhJetdPtPythia)
  , fhJetPtPythia(ap.fhJetPtPythia)
  , fhJetPtBinPythia(ap.fhJetPtBinPythia)
  , fhJetMultiplicityBinPythia(ap.fhJetMultiplicityBinPythia)
  , fhZPythia(ap.fhZPythia) 
  , fhZBinPythia(ap.fhZBinPythia) 
  , fhZBinPythiaRCut(ap.fhZBinPythiaRCut) 
  , fhJtPythia(ap.fhJtPythia)
  , fhJtBinPythia(ap.fhJtBinPythia)
  , fhJtBinPythiaRCut(ap.fhJtBinPythiaRCut)
  , fhJtWeightBinPythia(ap.fhJtWeightBinPythia)
  , fhLogJtWeightBinPythia(ap.fhLogJtWeightBinPythia)
  , fhLogJtWeight2BinPythia(ap.fhLogJtWeight2BinPythia)
  , fhJtWithPtCutWeightBinBinPythia(ap.fhJtWithPtCutWeightBinBinPythia)
  , fhLogJtWithPtCutWeightBinBinPythia(ap.fhLogJtWithPtCutWeightBinBinPythia)
  , fhLogJtWithPtCutWeight2BinBinPythia(ap.fhLogJtWithPtCutWeight2BinBinPythia)
  , fhBgJtPythia(ap.fhBgJtPythia)
  , fhBgJtPythiaTrackType(ap.fhBgJtPythiaTrackType)
  , fhBgJtBinPythia(ap.fhBgJtBinPythia)
  , fhBgJtBinPythiaTypeBin(ap.fhBgJtBinPythiaTypeBin)
  , fhBgJtWeightBinPythia(ap.fhBgJtWeightBinPythia)
  , fhBgLogJtWeightBinPythia(ap.fhBgLogJtWeightBinPythia)
  , fhBgLogJtWeight2BinPythia(ap.fhBgLogJtWeight2BinPythia)
  , fhBgJtWithPtCutWeightBinBinPythia(ap.fhBgJtWithPtCutWeightBinBinPythia)
  , fhBgLogJtWithPtCutWeightBinBinPythia(ap.fhBgLogJtWithPtCutWeightBinBinPythia)
  , fhBgLogJtWithPtCutWeight2BinBinPythia(ap.fhBgLogJtWithPtCutWeight2BinBinPythia)
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
  //fTrackJt.clear();
  //fTrackPt.clear();
  //fJetPt.clear();
  delete fEfficiency;
  delete fHMG;
  delete fHMGMC;



}



/// Create histograms
void AliJJetJtAnalysis::UserCreateOutputObjects(){
  //TODO Use .setWith
  cout << "AliJJetJetAnalysis UserCreateOutputObjects" << endl;
  //fJetListOfList always point one address in the whole time of this analysis.
  //Thus mustn't be cleared in it's life.     
  //fJetListOfList.Clear();


  fJJetAnalysis = new AliJJetAnalysis();

  fJetTriggPtBorders = fCard->GetVector("JetTriggPtBorders");
  fJetConstPtLowLimits = fCard->GetVector("JetConstPtLowLimits");
  fJetAssocPtBorders = fCard->GetVector("JetAssocPtBorders");
  fJetLeadPtBorders = fCard->GetVector("JetLeadPtBorders");
  fJetMultBorders = fCard->GetVector("JetMultBorders");
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


  if(fDoMC){
    CreateMCHistograms();
  }

  fHMG = new AliJHistManager( "AliJJetJtHistManager");
  fJetFinderBin .Set("JetFinderOrder","NFin","NFin:%d", AliJBin::kSingle).SetBin(nJetContainer);
  fJetTriggerBin .Set("JetTriggerBin","JetPt","p_{T,jet} : %.1f - %.1f").SetBin(fCard->GetVector("JetTriggPtBorders"));
  fTrkPtBin .Set("TrkPtBin","TrkPt","p_{T,constituent}:%.1f-%.1f").SetBin(fCard->GetVector("JetAssocPtBorders"));
  fTrkLimPtBin .Set("TrkLimitPtBin","TrkLimitPt","p_{T,Limit}<%.1f", AliJBin::kSingle).SetBin(fJetConstPtLowLimits->GetNoElements());
  fJetLeadPtBin.Set("LeadPtBin","LeadPt","p_{T,leading}:%.1f-%.1f").SetBin(fCard->GetVector("JetLeadPtBorders"));
  fJetMultBin.Set("MultBin","Multiplicity","N_{const.}:%.1f-%.1f").SetBin(fCard->GetVector("JetMultBorders"));
  fdRBin.Set("dRBin","dR","dR : %.1f - %.1f ").SetBin(fCard->GetVector("DeltaRBorders"));
  fiHist.Set("iHist","iHist","iHist : %d ", AliJBin::kSingle).SetBin(10);


  int NBINS=150;
  double LogBinsX[NBINS+1], LimL=0.1, LimH=500;
  double logBW = (log(LimH)-log(LimL))/NBINS;
  for(int ij=0;ij<=NBINS;ij++) LogBinsX[ij]=LimL*exp(ij*logBW);


  // Dijet Histograms - Oskari
  fhDiJetM
    << TH1D("DiJetM","Di-jet invariant mass",NBINS, LogBinsX ) 
    << fJetFinderBin
    <<"END";

  fhNumber
    << TH1D("hNumber","Number",20,0,20)
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

  fhJetPtMultiplicityCutBin 
    << TH1D("JetPtMultiplicityCutBin","",NBINS, LogBinsX ) 
    << fJetFinderBin << fJetTriggerBin << fJetMultBin
    <<"END";

  fhJetPtTrackCutBin
    << TH1D("JetPtTrackCutBin","",NBINS, LogBinsX ) 
    << fJetFinderBin << fJetTriggerBin << fJetLeadPtBin 
    <<"END";

  fhJetPtWeight 
    << TH1D("JetPtWeight","",NBINS, LogBinsX ) 
    << fJetFinderBin
    <<"END";
  fhJetPtWeightBin 
    << TH1D("JetPtWeightBin","",NBINS, LogBinsX ) 
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  fhJetMultiplicityBin 
    << TH1D("JetMultiplicityBin","Jet Mulplicity",50,0,50)
    << fJetFinderBin << fJetTriggerBin
    <<"END";

  fhJetMultiplicityBin2
    << TH1D("JetMultiplicityBin2","Jet Mulplicity",50,0,50)
    << fJetFinderBin << fJetTriggerBin
    <<"END";

  fhJetConeMultiplicityBin
    << TH1D("JetConeMultiplicityBin","Jet Cone Mulplicity",50,0,50)
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

  int NBINSR = 100;
  double LimLR = 0;
  double LimHR = 2;


  fhRBin
    << TH1D("RBin","",NBINSR,LimLR,LimHR)
    << fJetFinderBin << fJetTriggerBin
    << "END";

    
  fhBgRBin
    << TH1D("BgRBin","",NBINSR,LimLR,LimHR)
    << fJetFinderBin << fJetTriggerBin
    << "END";
  fhBgRndmRBin
    << TH1D("BgRndmRBin","",NBINSR,LimLR,LimHR)
    << fJetFinderBin << fJetTriggerBin
    << "END";

  int NBINSJt=150;
  double LogBinsJt[NBINSJt+1], LimLJt=0.01, LimHJt=10;
  double logBWJt = (TMath::Log(LimHJt)-TMath::Log(LimLJt))/NBINSJt;
  for(int ij=0;ij<=NBINSJt;ij++) LogBinsJt[ij]=LimLJt*exp(ij*logBWJt);
  int NBINSJtW=150;
  double LimLJtW=TMath::Log(0.01), LimHJtW=TMath::Log(10);

  fhEventJtBin 
    << TH1D("EventJtBin","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fJetTriggerBin
    <<"END";

  fhEventJtWeightBin 
    << TH1D("EventJtWeightBin","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fJetTriggerBin
    <<"END";

  fhEventJtWithPtCutWeightBinBin
    << TH1D("EventJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fJetTriggerBin << fTrkPtBin
    <<"END";

  fhLeadingJt 
    << TH1D("LeadingJt","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin
    <<"END";
  fhLeadingJtBin 
    << TH1D("LeadingJtBin","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  fhLeadingJtWeightBin 
    << TH1D("LeadingJtWeightBin","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  fhLeadingJtWithPtCutWeightBinBin
    << TH1D("LeadingJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fJetTriggerBin << fTrkPtBin
    <<"END";
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
  fhJtWeightBinTest
    << TH1D("JtWeightBinTest","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  fhJtWeightBinTest2
    << TH1D("JtWeightBinTest2","",NBINSJt, LogBinsJt ) 
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
  fhJetConeJtWeightWithTrackCutBinBin
    << TH1D("JetConeJtWeightWithTrackCutBinBin","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fJetTriggerBin << fJetLeadPtBin 
    <<"END";
  fhJetConeJtWeightWithMultiplicityCutBinBin
    << TH1D("JetConeJtWeightWithMultiplicityCutBinBin","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fJetTriggerBin << fJetMultBin 
    <<"END";

  fhJtWithPtCutWeightBinBin 
    << TH1D("JtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fJetTriggerBin << fTrkPtBin
    <<"END";
  fhJetConeJtWithPtCutWeightBinBin 
    << TH1D("JetConeJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
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

  fhBgJtWithPtCutWeightBinBin 
    << TH1D("BgJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fJetTriggerBin << fTrkPtBin
    <<"END";

  fhBgJtWithPtCutWeightBinBinSmallerR
    << TH1D("BgJtWithPtCutWeightBinBinSmallerR","",NBINSJt, LogBinsJt ) 
    << fiHist << fJetTriggerBin << fTrkPtBin 
    <<"END";

  fhBgJtWithPtCutWeightBinBinDiffR
    << TH1D("BgJtWithPtCutWeightBinBinDiffR","",NBINSJt, LogBinsJt ) 
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

  fhBgRndmJt 
    << TH1D("BgRndmJt","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin
    <<"END";
  fhBgRndmJtWeight
    << TH1D("BgRndmJtWeight","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin
    <<"END";
  fhBgRndmJtBin
    << TH1D("BgRndmJtBin","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  fhBgRndmJtWeightBin
    << TH1D("BgRndmJtWeightBin","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  fhBgRndmJtWithPtCutWeightBin
    << TH1D("BgRndmJtWithPtCutWeightBin","",NBINSJt, LogBinsJt ) 
    << fJetFinderBin << fTrkPtBin
    <<"END";
  fhBgRndmJtWithPtCutWeightBinBin
    << TH1D("BgRndmJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) 
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

  fhLeadingTrkPtBin
    << TH1D("LeadingTrkPtBin","",NBINSPt,LogBinsPt)
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

  if(fDoLog){
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
    fhLogJtWeightBinLimBin 
      << TH1D("LogJtWeightBinLimBin","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
      <<"END";
    fhLogJtWeight2BinLimBin 
      << TH1D("LogJtWeight2BinLimBin","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
      <<"END";
    fhBgLogJtWeightBin 
      << TH1D("BgLogJtWeightBin","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBin << fJetTriggerBin
      <<"END";
    fhBgLogJtWeight2Bin 
      << TH1D("BgLogJtWeight2Bin","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBin << fJetTriggerBin
      <<"END";
    fhBgLogJtWithPtCutWeightBinBin 
      << TH1D("BgLogJtWithPtCutWeightBinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBin << fJetTriggerBin << fTrkPtBin
      <<"END";
    fhBgLogJtWithPtCutWeight2BinBin 
      << TH1D("BgLogJtWithPtCutWeight2BinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBin << fJetTriggerBin << fTrkPtBin
      <<"END";
    fhBgLogJtWithPtCutWeightBinBinSmallerR 
      << TH1D("BgLogJtWithPtCutWeightBinBinSmallerR","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fiHist << fJetTriggerBin << fTrkPtBin 
      <<"END";
    fhBgLogJtWithPtCutWeight2BinBinSmallerR 
      << TH1D("BgLogJtWithPtCutWeight2BinBinSmallerR","",NBINSJtW, LimLJtW, LimHJtW ) 
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
    fhBgLogJtWeightBinLimBin 
      << TH1D("BgLogJtWeightBinLimBin","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
      <<"END";
    fhBgLogJtWeight2BinLimBin 
      << TH1D("BgLogJtWeight2BinLimBin","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
      <<"END";
    fhBgRndmLogJt 
      << TH1D("BgRndmLogJt","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBin
      <<"END";
    fhBgRndmLogJtWithPtCutWeight2Bin 
      << TH1D("BgRndmLogJtWithPtCutWeight2Bin","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBin << fTrkPtBin
      <<"END";
    fhBgRndmLogJtWithPtCutWeight2BinBin 
      << TH1D("BgRndmLogJtWithPtCutWeight2BinBin","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBin << fJetTriggerBin << fTrkPtBin
      <<"END";
  }


  fHMG->Print();
  fHMG->WriteConfig();

  cout << "UserCreateOutputObjects done " << endl;



}

/// Create Correlation and pythia histograms
void AliJJetJtAnalysis::CreateMCHistograms(){

  int NBINS=150;
  double LogBinsX[NBINS+1], LimL=0.1, LimH=500;
  double logBW = (log(LimH)-log(LimL))/NBINS;
  for(int ij=0;ij<=NBINS;ij++) LogBinsX[ij]=LimL*exp(ij*logBW);

  fHMGMC = new AliJHistManager( "AliJJetJtMCHistManager");
  fJetFinderBinMC .Set("JetFinderOrder","NFin","NFin:%d", AliJBin::kSingle).SetBin(nJetContainer);
  fJetTriggerBinMC .Set("JetTriggerBin","JetPt","p_{T,jet} : %.1f - %.1f").SetBin(fCard->GetVector("JetTriggPtBorders"));
  fTrkPtBinMC .Set("TrkPtBin","TrkPt","p_{T,constituent}:%.1f-%.1f").SetBin(fCard->GetVector("JetAssocPtBorders"));
  fTrkLimPtBinMC .Set("TrkLimitPtBin","TrkLimitPt","p_{T,Limit}<%.1f", AliJBin::kSingle).SetBin(fJetConstPtLowLimits->GetNoElements());
  fJetLeadPtBinMC.Set("LeadPtBin","LeadPt","p_{T,leading}:%.1f-%.1f").SetBin(fCard->GetVector("JetLeadPtBorders"));
  fJetMultBinMC.Set("MultBin","Multiplicity","N_{const.}:%.1f-%.1f").SetBin(fCard->GetVector("JetMultBorders"));
  fdRBinMC.Set("dRBin","dR","dR : %.1f - %.1f ").SetBin(fCard->GetVector("DeltaRBorders"));
  fiHistMC.Set("iHist","iHist","iHist : %d ", AliJBin::kSingle).SetBin(10);
  fBgTypeBin.Set("BgTypeBin","BgTypeBin","BgType: %d", AliJBin::kSingle).SetBin(6).Print();

  int NBINSC = fJetTriggPtBorders->GetNoElements();
  //cout << "Number of elements in fJetTriggPtBorders: " << NBINSC << endl;
  double BinsC[NBINSC];
  for(int i = 0 ; i < NBINSC ; i++){
    BinsC[i] = (*fJetTriggPtBorders)[i+1];
    //cout << "Bin " << i << " set to " << BinsC[i] << endl;
  }


  fhBgJtPythiaTrackType
    << TH1D("BgJtPythiaTrackType","0: No mother found, 1: Mother track is outside acceptance, 2: ijet == icon, 3: ijet < 0, 4: R > 0.5, 5: Origin pT was too small to be a jet",6,0,6)
    << "END";

  fhTrackMatchSuccess
    << TH1D("TrackMatchSuccess","",6,0,6)
    << fJetFinderBinMC << fJetTriggerBinMC
    << "END";

  fhConstMatchSuccess
    << TH1D("ConstMatchSuccess","",6,0,6)
    << fJetFinderBinMC << fJetTriggerBinMC
    << "END";

  fhTrackMatchSuccessPythia
    << TH1D("TrackMatchSuccessPythia","",6,0,6)
    << fJetTriggerBinMC
    << "END";

  fhJetPtCorr
    << TH2D("JetPtCorr","",NBINS,LogBinsX,NBINS,LogBinsX)
    << fJetFinderBinMC
    << "END";

  fhJetPtCorrCoarse
    << TH2D("JetPtCorrCoarse","",NBINSC-1,BinsC,NBINSC-1,BinsC)
    << fJetFinderBinMC
    << "END";

  fhJetPtCorrPythia
    << TH2D("JetPtCorrPythia","",NBINS,LogBinsX,NBINS,LogBinsX)
    << "END";

  fhJetPtCorrPythiaCoarse
    << TH2D("JetPtCorrPythiaCoarse","",NBINSC-1,BinsC,NBINSC-1,BinsC)
    << "END";

  fhJetPtCorr2
    << TH2D("JetPtCorr2","",NBINS,LogBinsX,NBINS,LogBinsX)
    << fJetFinderBinMC
    << "END";

  fhJetPtCorr3
    << TH2D("JetPtCorr3","",NBINS,LogBinsX,NBINS,LogBinsX)
    << fJetFinderBinMC
    << "END";

  fhJetPtCorr2Pythia
    << TH2D("JetPtCorr2Pythia","",NBINS,LogBinsX,NBINS,LogBinsX)
    << "END";

  fhJetdR
    << TH1D("hJetdR","",50,0,0.5)
    << fJetFinderBinMC << fJetTriggerBinMC
    << "END";

  fhJetdRPythia
    << TH1D("hJetdRPythia","",50,0,0.5)
    << fJetTriggerBinMC
    << "END";

  fhJetdPt
    << TH1D("hJetdPt","",50,0,0.5)
    << fJetFinderBinMC << fJetTriggerBinMC
    << "END";

  fhJetdPtPythia
    << TH1D("hJetdPtPythia","",50,0,0.5)
    << fJetTriggerBinMC
    << "END";

  fhDiJetMjjCorr
    << TH2D("JetDiJetMjjCorr","",NBINS,LogBinsX,NBINS,LogBinsX)
    << fJetFinderBinMC
    << "END";

  int NBINSZ=150;
  double LogBinsZ[NBINSZ+1], LimLZ=0.001, LimHZ=1.1;
  double logBWZ = (TMath::Log(LimHZ)-TMath::Log(LimLZ))/NBINSZ;
  for(int ij=0;ij<=NBINSZ;ij++) LogBinsZ[ij]=LimLZ*exp(ij*logBWZ);//


  fhZPythia
    << TH1D("ZPythia","",NBINSZ, LogBinsZ ) 
    <<"END";
  fhZBinPythia
    << TH1D("ZBinPythia","",NBINSZ, LogBinsZ ) 
    << fJetTriggerBinMC
    <<"END";

  fhZBinPythiaRCut
    << TH1D("ZBinPythiaRCut","",NBINSZ, LogBinsZ ) 
    << fJetTriggerBinMC
    <<"END";

  int NBINSR = 100;
  double LimLR = 0;
  double LimHR = 2;

  fhRBinPythia
    << TH1D("RBinPythia","",NBINSR,LimLR,LimHR)
    << fJetTriggerBinMC
    << "END";

  int NBINSJt=150;
  double LogBinsJt[NBINSJt+1], LimLJt=0.01, LimHJt=10;
  double logBWJt = (TMath::Log(LimHJt)-TMath::Log(LimLJt))/NBINSJt;
  for(int ij=0;ij<=NBINSJt;ij++) LogBinsJt[ij]=LimLJt*exp(ij*logBWJt);
  int NBINSJtW=150;
  double LimLJtW=TMath::Log(0.01), LimHJtW=TMath::Log(10);

  fhTrackJtCorrBin
    << TH2D("TrackJtCorrBin","",NBINSJt,LogBinsJt,NBINSJt,LogBinsJt)
    << fJetFinderBinMC << fJetTriggerBinMC
    << "END";

  /*fhTrackJtCorr2D
    << TH2D("TrackJtCorr2D","",NBINSJt,LogBinsJt,NBINSJt,LogBinsJt,NBINSC,BinsC)
    << fJetFinderBinMC << fJetTriggerBinMC
    << "END";*/

  fhConstJtCorrBin
    << TH2D("ConstJtCorrBin","",NBINSJt,LogBinsJt,NBINSJt,LogBinsJt)
    << fJetFinderBinMC << fJetTriggerBinMC
    << "END";

  fhTrackJtCorrBinTest
    << TH2D("TrackJtCorrBinTest","",NBINSJt,LogBinsJt,NBINSJt,LogBinsJt)
    << fJetFinderBinMC << fJetTriggerBinMC
    << "END";
  fhTrackJtCorrBinTest2
    << TH2D("TrackJtCorrBinTest2","",NBINSJt,LogBinsJt,NBINSJt,LogBinsJt)
    << fJetFinderBinMC << fJetTriggerBinMC
    << "END";

  fhTrackJtCorrBinPythia
    << TH2D("TrackJtCorrBinPythia","",NBINSJt,LogBinsJt,NBINSJt,LogBinsJt)
    << fJetTriggerBinMC
    << "END";

  fhJtUnfBg
    << TH1D("JtUnfBg","",NBINSJt, LogBinsJt ) 
    << fJetFinderBinMC
    <<"END";
  fhJtBinUnfBg
    << TH1D("JtBinUnfBg","",NBINSJt, LogBinsJt ) 
    << fJetFinderBinMC << fJetTriggerBinMC
    <<"END";
  fhJtWeightBinUnfBg
    << TH1D("JtWeightBinUnfBg","",NBINSJt, LogBinsJt ) 
    << fJetFinderBinMC << fJetTriggerBinMC
    <<"END";

  fhJetConeJtUnfBg
    << TH1D("JetConeJtUnfBg","",NBINSJt, LogBinsJt ) 
    << fJetFinderBinMC
    <<"END";
  fhJetConeJtBinUnfBg
    << TH1D("JetConeJtBinUnfBg","",NBINSJt, LogBinsJt ) 
    << fJetFinderBinMC << fJetTriggerBinMC
    <<"END";
  fhJetConeJtWeightBinUnfBg
    << TH1D("JetConeJtWeightBinUnfBg","",NBINSJt, LogBinsJt ) 
    << fJetFinderBinMC << fJetTriggerBinMC
    <<"END";



  int NBINSPt=150;
  double LogBinsPt[NBINSPt+1], LimLPt=0.01, LimHPt=50;
  double logBWPt = (TMath::Log(LimHPt)-TMath::Log(LimLPt))/NBINSPt;
  for(int ij=0;ij<=NBINSPt;ij++) LogBinsPt[ij]=LimLPt*exp(ij*logBWPt);

  fhTrackPtCorr
    << TH2D("TrackPtCorr","",NBINSPt,LogBinsPt,NBINSPt,LogBinsPt)
    << fJetFinderBinMC
    << "END";
  fhConstPtCorr
    << TH2D("ConstPtCorr","",NBINSPt,LogBinsPt,NBINSPt,LogBinsPt)
    << fJetFinderBinMC
    << "END";

  fhTrackPtCorrPythia
    << TH2D("TrackPtCorrPythia","",NBINSPt,LogBinsPt,NBINSPt,LogBinsPt)
    << "END";


  fhJetMultiplicityBinPythia
    << TH1D("JetMultiplicityBinPythia","Jet Mulplicity",50,0,50)
    << fJetTriggerBinMC
    <<"END";

  fhJetPtPythia
    << TH1D("JetPtPythia","",NBINS, LogBinsX ) 
    <<"END";

  fhJetPtBinPythia
    << TH1D("JetPtBinPythia","",NBINS, LogBinsX ) 
    << fJetTriggerBinMC
    <<"END";

  fhJtPythia
    << TH1D("JtPythia","",NBINSJt, LogBinsJt ) 
    <<"END";
  fhJtBinPythia
    << TH1D("JtBinPythia","",NBINSJt, LogBinsJt ) 
    << fJetTriggerBinMC
    <<"END";

  fhJtBinPythiaRCut
    << TH1D("JtBinPythiaRCut","",NBINSJt, LogBinsJt ) 
    << fJetTriggerBinMC
    <<"END";

  fhJtWeightBinPythia
    << TH1D("JtWeightBinPythia","",NBINSJt, LogBinsJt ) 
    << fJetTriggerBinMC
    <<"END";


  fhJtWithPtCutWeightBinBinPythia
    << TH1D("JtWithPtCutWeightBinBinPythia","",NBINSJt, LogBinsJt ) 
    << fJetTriggerBinMC << fTrkPtBinMC
    <<"END";
  fhBgJtPythia
    << TH1D("BgJtPythia","",NBINSJt, LogBinsJt ) 
    <<"END";

  fhBgJtBinPythia
    << TH1D("BgJtBinPythia","",NBINSJt, LogBinsJt ) 
    << fJetTriggerBinMC
    <<"END";
  fhBgJtBinPythiaTypeBin
    << TH1D("BgJtBinPythiaTypeBin","",NBINSJt, LogBinsJt ) 
    << fJetTriggerBinMC << fBgTypeBin
    <<"END";

  fhBgJtWeightBinPythia
    << TH1D("BgJtWeightBinPythia","",NBINSJt, LogBinsJt ) 
    << fJetTriggerBinMC
    <<"END";


  fhBgJtWithPtCutWeightBinBinPythia
    << TH1D("BgJtWithPtCutWeightBinBinPythia","",NBINSJt, LogBinsJt ) 
    << fJetTriggerBinMC << fTrkPtBinMC
    <<"END";

  if(fDoLog){
    fhLogJtWeightBinUnfBg
      << TH1D("LogJtWeightBinUnfBg","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBinMC << fJetTriggerBinMC
      <<"END";
    fhLogJtWeight2BinUnfBg 
      << TH1D("LogJtWeight2BinUnfBg","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBinMC << fJetTriggerBinMC
      <<"END";
    fhJetConeLogJtWeightBinUnfBg
      << TH1D("JetConeLogJtWeightBinUnfBg","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBinMC << fJetTriggerBinMC
      <<"END";
    fhJetConeLogJtWeight2BinUnfBg 
      << TH1D("JetConeLogJtWeight2BinUnfBg","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetFinderBinMC << fJetTriggerBinMC
      <<"END";
    fhLogJtWeightBinPythia
      << TH1D("LogJtWeightBinPythia","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetTriggerBinMC
      <<"END";
    fhLogJtWeight2BinPythia
      << TH1D("LogJtWeight2BinPythia","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetTriggerBinMC
      <<"END";
    fhLogJtWithPtCutWeightBinBinPythia
      << TH1D("LogJtWithPtCutWeightBinBinPythia","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetTriggerBinMC << fTrkPtBinMC
      <<"END";

    fhLogJtWithPtCutWeight2BinBinPythia
      << TH1D("LogJtWithPtCutWeight2BinBinPythia","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetTriggerBinMC << fTrkPtBinMC
      <<"END";
    fhBgLogJtWeightBinPythia
      << TH1D("BgLogJtWeightBinPythia","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetTriggerBinMC
      <<"END";
    fhBgLogJtWeight2BinPythia
      << TH1D("BgLogJtWeight2BinPythia","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetTriggerBinMC
      <<"END";
    fhBgLogJtWithPtCutWeightBinBinPythia
      << TH1D("BgLogJtWithPtCutWeightBinBinPythia","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetTriggerBinMC << fTrkPtBinMC
      <<"END";

    fhBgLogJtWithPtCutWeight2BinBinPythia
      << TH1D("BgLogJtWithPtCutWeight2BinBinPythia","",NBINSJtW, LimLJtW, LimHJtW ) 
      << fJetTriggerBinMC << fTrkPtBinMC
      <<"END";
  }

  fHMGMC->Print();
  fHMGMC->WriteConfig();
  TDirectory *thisDir = fHMGMC->GetDirectory();
  TDirectory *motherDir = thisDir->GetMotherDir();
  motherDir->cd();



}

void AliJJetJtAnalysis::ClearBeforeEvent(){
  //fJetListOfList.Clear();


}

///UserExec is executed for every event
void AliJJetJtAnalysis::UserExec(){
  //Initialize Vectors
  fTrackJt = new TVector(fTracks->GetEntries()); //Store jT data of reconstructed tracks
  fConstJt = new TVector(fTracks->GetEntries()); //Store jT data of reconstructed jet constituents
  fConstLabels = new TVector(fTracks->GetEntries()); //Store labels of jet constituents
  fTrackPt = new TVector(fTracks->GetEntries()); //pT data of reconstructed tracks
  fConstPt = new TVector(fTracks->GetEntries()); //pT data of reconstructed tracks
  fJetPt = new TVector(fTracks->GetEntries()); //parent Jet pT data of reconstructed tracks
  fTrackFound = new TVector(fTracks->GetEntries()); //Keep track of matched tracks
  fConstFound = new TVector(fTracks->GetEntries()); //Keep track of matched tracks
  fBin2 = new TVector(fTracks->GetEntries());  //Store iBin2 of corresponding tracks 
  fBin3 = new TVector(fTracks->GetEntries()); //Store iBin3 of corresponding tracks
  fpta = new TVector(fTracks->GetEntries());  //Store pTa bin of corresponding tracks
  fptt = new TVector(fTracks->GetEntries());  //Store pTt bin of corresponding tracks

  fDiJetMjj = new TVector(fJetListOfList.GetEntries()); //Di-jet invariant mass for each jet finder.

  int nR = 0;
  if(fDoMC){
    nR = fJetListOfList.GetEntries()/4;
  }else{
    nR = fJetListOfList.GetEntries()/2;
  }
  //Loop over jet finders (Full/charged jets R=0.4,0.5,0.6, Pythia/Reco)
  for( int i=0;i<fJetListOfList.GetEntries();i++ ){
    TObjArray * Jets = (TObjArray*) fJetListOfList[i];
    if(!Jets) {
      continue;
    }
    if(GetTrackOrMCParticle(i) == kJMCParticle){ //Check if jet finder is Detector or particle level
      this->FillJtHistogram(Jets,i,1); //Particle
    }else{
      this->FillJtHistogram(Jets,i,0); //Detector
    }

    if(fDoMC){
      //fTrackPt.clear();
      if(i+1>nR*2) continue;
      TObjArray * Jets = (TObjArray*) fJetListOfList[i];
      TObjArray * MCJets = (TObjArray*) fJetListOfList[i+nR*2];
      int doCor = 1;
      if(!MCJets) {
        doCor = 0;
      }
      if(!Jets) {
        doCor = 0;
      }
      if(doCor){
        this->FillCorrelation(Jets,MCJets,i,i+nR*2);  //Correlation histograms between reco and particle level
      }
      if(i == 0){
        this->FillPythia(Jets,i); //Fill histogram for jet jT using Pythia fragmentation information, including correlation between reco and pythia
      }
    }
  }

  //Fill correlation matrix for di-jet invariant mass between reco and particle level.
  if(fDoMC){
    for( int i=0;i<fJetListOfList.GetEntries()/2;i++ ){
      if((*fDiJetMjj)[i+nR*2] > 0 ) { //MC jets
        if((*fDiJetMjj)[i] > 0 ) {  //jets
          fhDiJetMjjCorr[i]->Fill((*fDiJetMjj)[i+nR*2],(*fDiJetMjj)[i]);
        } else {
          fhDiJetMjjCorr[i]->Fill((*fDiJetMjj)[i+nR*2],-1);
        }
      }
      /*if((*fDiJetMjj)[i+nR*2] == 0 && (*fDiJetMjj)[i] > 0){ //If event has a measured dijet mass but no true dijet mass
        cout << "Found measured Dijet M but no true M!!!!! " << endl << endl << endl;
        cout << "Measured M" << (*fDiJetMjj)[i] << endl;
        }*/
    }
  }



  //The Function should be called after calling "FillJtHistogram"
  //FillBgJtWithSmallerR(Bg jet finder array, old R, new R )
  //fJetBgListOfList[i] where {i=0-5;full0.4,full0.5,full0.6,Ch0.4,Ch0.5,Ch0.6}
  //Caution!! these array number should be changed WHEN jet finders change

  /*this->FillBgJtWithSmallerR(fJetBgListOfList[1], 0.4,0);
    this->FillBgJtWithSmallerR(fJetBgListOfList[2], 0.4,1);
    this->FillBgJtWithSmallerR(fJetBgListOfList[2], 0.5,2);
    this->FillBgJtWithSmallerR(fJetBgListOfList[4], 0.4,3);
    this->FillBgJtWithSmallerR(fJetBgListOfList[5], 0.4,4);
    this->FillBgJtWithSmallerR(fJetBgListOfList[5], 0.5,5);*/

  //Fill jt with diff cone axes (old axis iContainer, new axis, iHist) 
  /*this->FillBgJtWithDiffAxes(1, 0,0);
    this->FillBgJtWithDiffAxes(2, 0,1);
    this->FillBgJtWithDiffAxes(2, 1,2);
    this->FillBgJtWithDiffAxes(4, 3,0);
    this->FillBgJtWithDiffAxes(5, 3,1);
    this->FillBgJtWithDiffAxes(5, 4,2);*/
  //End.


  int iS1 = 0; //full 0.4
  int iS2 = nR; //Ch   0.4
  TObjArray * jetfinder1;
  TObjArray * jetfinder2;
  /*if(fDoMC){
    jetfinder1 = (TObjArray*) fMCJetListOfList[iS1];
    jetfinder2 = (TObjArray*) fMCJetListOfList[iS2];
    }else{*/
  jetfinder1 = (TObjArray*) fJetListOfList[iS1];
  jetfinder2 = (TObjArray*) fJetListOfList[iS2];
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

/// Does nothing
void AliJJetJtAnalysis::WriteHistograms(){

  //TDirectory * cwd = gDirectory;
  //const int nJetContainer = fJetListOfList.GetEntries();

  for (int i=0; i<nJetContainer; i++){
    TDirectory *nwd = gDirectory->mkdir(fJetFinderName[i]);
    //Under the folder name, save objects
    //nwd->cd();
    //cwd->cd();
  }

}


/// Fill jt to  histograms and calculate bg jt and fill it into histograms
/// \param Jets List of jets used
/// \param iContainer Index of jet finder
/// \param mc Whether or not this is particle level data
void AliJJetJtAnalysis::FillJtHistogram( TObjArray *Jets , int iContainer,int mc)
{	
  TClonesArray *trackArray;
  if(mc)
    trackArray = fMCTracks;
  else 
    trackArray = fTracks;

  int iBin, iptaBin=0;
  int jBin=0;
  int iBin2=0;
  int iBin3=0;
  double pT = 0;
  double conPtMax =0;


  double z; double jt;
  double pta;
  //double Y , deltaY = 0;
  //double deltaR= 0;
  //cout<<"histogram filling number of jets : "<<Jets->GetEntriesFast()<<endl;

  TLorentzVector  vOrtho;
  TLorentzVector  randomTrack;
  TLorentzVector summedJet;
  fJetBgListOfList[iContainer].Clear();
  TClonesArray & bgjets = fJetBgListOfList[iContainer];

  double deltaR = -1;
  double effCorrection = -1;
  double thisConeSize = fConeSizes[iContainer] ;
  double leadingTrackPt = 0;
  double leadingTrackJt = 0;
  double leadingTrackEff = 0;
  int leadingTrackIndex = 0;
  int iBgJet = 0;
  int nC = 0;

  //fTrackJt.reserve(fTracks->GetEntries());
  //fTrackPt.reserve(fTracks->GetEntries());
  //fJetPt.reserve(fTracks->GetEntries());

  int doRndmBg = 0;
  // iJet loop for an event
  if(!mc){
    iConst = 0;
    for (int icon = 0; icon<trackArray->GetEntries(); icon++){
      (*fTrackJt)[icon] = -1;
      (*fConstJt)[icon] = -1;
      (*fConstPt)[icon] = -1;
      (*fConstLabels)[icon] = -1;
      (*fConstFound)[icon] = -1;
      (*fJetPt)[icon] = -1;
      (*fTrackFound)[icon] = -1;
      AliJBaseTrack *track = dynamic_cast<AliJBaseTrack*>(trackArray->At(icon));
      if (!track) continue;
      (*fTrackPt)[icon] = track->Pt();
    }
  }

  if(Jets->GetEntries() == 0){
    //cout << "No Jets!" << endl;
    //cout << "iContainer: " << iContainer << " No jets" << endl;
    return;
  }

  // By Oskari Saarimäki
  // Taking the leading and subleading jet and calculating dijet mass.
  AliJJet *leadingjet;
  AliJJet *subleadingjet;
  TLorentzVector dijet;
  if (Jets->GetEntries()>1) { // Need at least two jets.
    leadingjet = dynamic_cast<AliJJet*>( Jets->At(0) ); // Highest p_T
    subleadingjet = dynamic_cast<AliJJet*>( Jets->At(1) ); // Second highest p_T
    dijet = leadingjet->GetLorentzVector() + subleadingjet->GetLorentzVector(); // Adding the jets together.
    if (leadingjet->Pt() > 20 && subleadingjet->Pt() > 20){ 
        (*fDiJetMjj)[iContainer] = dijet.M();
        //cout << "Dijet mass in an event: " << dijet.M() << ", leading pt: " << leadingjet->Pt() << ", subleading pt: " << subleadingjet->Pt() << endl;
        fhDiJetM[iContainer]->Fill((*fDiJetMjj)[iContainer]);
    } else { // For testing:
        //cout << "No dijet in an event!        leading pt: " << leadingjet->Pt() << ", subleading pt: " << subleadingjet->Pt() << endl;
    }
  }
  // End of dijet calculation.

  //Loop over jets
  //cout << "Container: " << iContainer << " Number of Jets: " << Jets->GetEntries() << endl;
  for (int i = 0; i<Jets->GetEntries(); i++){
    //cout << "Jet: " << (i+1) << "/" << Jets->GetEntries() << endl;
    AliJJet *jet = dynamic_cast<AliJJet*>( Jets->At(i) );
    if (!jet) continue;
    jet->ReSum(); 
    //if (TMath::Abs(jet->Eta()) > fJetEtaCut) {cout << "ijet: " << i << " Jet outside eta range, eta: " << jet->Eta() << " pT: " << jet->Pt() << endl; continue;}
    //if (TMath::Abs(jet->Eta()) > fJetEtaCut && i == 0) cout << "Leading jet skipped because outside eta range!, eta: " << jet->Eta() << endl;
    if (TMath::Abs(jet->Eta()) > fJetEtaCut) continue;
    //TODO Separate Eta cuts for charged and full jets
    //TODO Limit full jets to emcal acceptance Phi: from 1.855 to 2.685
    fhJetEtaPhi[iContainer]->Fill(jet->Eta(),jet->Phi());
    pT = jet->Pt();
    //cout << "ijet: " << i << " pT: " << pT  << endl;
    if (pT<(*fJetTriggPtBorders)[1]) continue;
    iBin = GetBin(fJetTriggPtBorders,pT); // fill jetPt histos
    if( iBin < 0 ) continue;
    doRndmBg++;
    int jetMult = jet->GetConstituents()->GetEntries();
    if(jetMult < 1){
      cout << "WARNING Jet multiplicity: " << jetMult << "!!!!!!" << endl;
      cout << "Skipping jet" << endl;
      continue;
    }
    if(i > 0 && fLeadingJets > 0) continue;
    fhJetPt[iContainer]->Fill( pT );
    fhJetPtBin[iContainer][iBin]->Fill( pT );
    fhJetPtWeight[iContainer]->Fill( pT , 1/pT);
    fhJetPtWeightBin[iContainer][iBin]->Fill( pT, 1/pT );
    fhJetMultiplicityBin[iContainer][iBin]->Fill(jetMult);
    leadingTrackPt = jet->LeadingParticlePt(); //FIXME? For MC tracks this is possibly a track with no charge
    fhLeadingTrkPtBin[iContainer][iBin]->Fill(leadingTrackPt);
    iBin2 = GetBin(fJetLeadPtBorders,leadingTrackPt);
    iBin3 = GetBin(fJetMultBorders,jetMult);
    if(iBin2 < 0 || iBin3 <0){
      cout << "Jet pT: " << pT << " Multiplicity: " << jetMult << " Leading track pT: " << leadingTrackPt << " Jet Finder: " << iContainer << endl;
    }
    if(iBin2 > -1){
      fhJetPtTrackCutBin[iContainer][iBin][iBin2]->Fill(pT);
    }
    if(iBin3 > -1){
      fhJetPtMultiplicityCutBin[iContainer][iBin][iBin3]->Fill(pT);
    }

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
    nC = 0;
    leadingTrackPt = 0;
    leadingTrackIndex = 0;
    leadingTrackEff = 0;
    leadingTrackJt = 0;
    for (int icon = 0; icon<jet->GetConstituents()->GetEntries(); icon++){
      AliJBaseTrack *constituent = jet->GetConstituent(icon);
      if(constituent->GetCharge() == 0 ) {
        //cout << "Track " << icon << "/" << jet->GetConstituents()->GetEntries() << " skipped because charge == 0, mc: " << mc << endl;  
        continue;
      }
      z = (constituent->Vect()*jet->Vect().Unit())/jet->P();
      pta = constituent->Pt();
      constituent->SetTrackEff( fEfficiency->GetCorrection( pta, 5, fcent) );
      if(mc){
        effCorrection = 1.0;
      }else{
        effCorrection = 1.0/constituent->GetTrackEff();
      }
      iptaBin = GetBin(fJetAssocPtBorders, pta);
      if( iptaBin < 0 ) {
        //cout << "Track " << icon << "/" << jet->GetConstituents()->GetEntries() << " skipped because iptaBin < 0, pta: " << pta << " mc: " << mc << endl;
        continue; 
      }

      fhTrkPt[iContainer]->Fill(pta,effCorrection);
      fhTrkPtBin[iContainer][iBin]->Fill(pta,effCorrection);
      fhTrkPtWeightBin[iContainer][iBin]->Fill(pta,effCorrection/pta);


      fhRBin[iContainer][iBin]->Fill(getDiffR(constituent->Phi(),jet->Phi(),constituent->Eta(),jet->Eta()),effCorrection);
      fhZ[iContainer]->Fill( z , effCorrection);
      fhZBin[iContainer][iBin]->Fill( z , effCorrection);
      jt = (constituent->Vect()-z*jet->Vect()).Mag();
      nC++;
      if(!mc){
        (*fConstJt)[iConst] = jt;
        (*fConstPt)[iConst] = pta;
        (*fConstLabels)[iConst] = constituent->GetLabel();
        iConst++;
      }
      fhJt[iContainer]->Fill( jt , effCorrection);
      fhJtBin[iContainer][iBin]->Fill( jt , effCorrection);
      fhJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt * effCorrection );
      if(fDoLog){
        fhLogJtWeightBin[iContainer][iBin]
          ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
        fhLogJtWeight2Bin[iContainer][iBin]
          ->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
      }



      //cout << "iContainer: " << iContainer << " iBin: " << iBin << " iptaBin: " << iptaBin << endl;
      fhJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
        ->Fill( jt, 1.0/jt * effCorrection );
      if(fDoLog){
        fhLogJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
          ->Fill( TMath::Log(jt), 1.0/jt * effCorrection);
        fhLogJtWithPtCutWeight2BinBin[iContainer][iBin][iptaBin]
          ->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection);
      }

      if(pta > leadingTrackPt){
        leadingTrackPt = pta;
        leadingTrackIndex = icon;
        leadingTrackJt = jt; //Save leading track jT for later
        leadingTrackEff = effCorrection;
      }


      for (int jj = 0; jj <= jBin ; jj++) {
        fhJtBinLimBin[iContainer][iBin][jj]->Fill( jt, effCorrection );
        fhJtWeightBinLimBin[iContainer][iBin][jj]
          ->Fill( jt, 1.0/jt * effCorrection );
        if(fDoLog){
          fhLogJtWeightBinLimBin[iContainer][iBin][jj]
            ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
          fhLogJtWeight2BinLimBin[iContainer][iBin][jj]
            ->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
        }
      }

    }
    if(leadingTrackPt > 0.1){
      iptaBin = GetBin(fJetAssocPtBorders, leadingTrackPt);
      if(iptaBin >= 0){
        fhLeadingJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
          ->Fill( leadingTrackJt, 1.0/leadingTrackJt * leadingTrackEff);
        fhLeadingJt[iContainer]->Fill(leadingTrackJt, leadingTrackEff);
        fhLeadingJtBin[iContainer][iBin]->Fill(leadingTrackJt, leadingTrackEff);
        fhLeadingJtWeightBin[iContainer][iBin]->Fill(leadingTrackJt, 1.0/leadingTrackJt * leadingTrackEff);
      }
    }

    fhJetMultiplicityBin2[iContainer][iBin]->Fill(nC);

    //vOrtho.SetVect(jet->Vect().Orthogonal());
    vOrtho.SetVect(jet->Vect());
    vOrtho.SetE(jet->E());
    vOrtho.SetPhi(jet->Phi()+TMath::Pi()/2);


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
        deltaR   = getDiffR(vOrtho.Phi(),jet2->Phi(),vOrtho.Eta(),jet2->Eta());
        if ( deltaR < 2*thisConeSize) {
          doBkg = 0;
          fhNumber[iContainer]->Fill(4.5);
          break;
        }

      }
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

    nC = 0;
    for (int icon = 0; icon<trackArray->GetEntries(); icon++){
      AliJBaseTrack *track = dynamic_cast<AliJBaseTrack*>(trackArray->At(icon));
      if (!track) continue;
      if(mc){
        if(!track->IsPrimary()) continue;
        if(track->GetCharge() == 0) continue;
        if(TMath::Abs(track->Eta()) > 1.0) continue;
      }else{
        if(!track->IsTrue(1)) continue;
      }
      fhTrackEtaPhi[iContainer]->Fill(track->Eta(),track->Phi());
      pta = track->Pt();
      if (pta > maxconpt) maxconpt = pta;
      track->SetTrackEff( fEfficiency->GetCorrection( pta, 5, fcent) );
      if(mc){
        effCorrection = 1.0;
      }
      else{
        effCorrection = 1.0/track->GetTrackEff();
      }
      iptaBin = GetBin(fJetAssocPtBorders, pta);
      if( iptaBin < 0 ) continue;


      z = (track->Vect()*jet->Vect().Unit())/jet->P();
      jt = (track->Vect()-z*jet->Vect()).Mag();
      //jT for all tracks in the event 

      deltaR   = getDiffR(jet->Phi(),track->Phi(),jet->Eta(),track->Eta());
      if(deltaR < TMath::Pi()/2){
        fhEventJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]->Fill(jt, 1.0/jt * effCorrection);
        fhEventJtWeightBin[iContainer][iBin]->Fill(jt, 1.0/jt * effCorrection); 
        fhEventJtBin[iContainer][iBin]->Fill(jt,effCorrection); 
      }
      //Jet Cone Jt here
      if ( deltaR < thisConeSize){ 
        fhJetConeTrkPt[iContainer]->Fill(pta,effCorrection);
        fhJetConeTrkPtBin[iContainer][iBin]->Fill(pta,effCorrection);
        fhJetConeTrkPtWeightBin[iContainer][iBin]->Fill(pta,effCorrection/pta);

        fhJetConeZ[iContainer]->Fill( z , effCorrection);
        fhJetConeZBin[iContainer][iBin]->Fill( z , effCorrection);
        nC++;
        if(!mc){
          (*fTrackJt)[icon] = jt;
          (*fJetPt)[icon] = jet->Pt();
        }
        fhJetConeJt[iContainer]->Fill( jt , effCorrection);
        fhJetConeJtBin[iContainer][iBin]->Fill( jt , effCorrection);
        fhJetConeJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt * effCorrection );
        if(iBin2 > -1){
          fhJetConeJtWeightWithTrackCutBinBin[iContainer][iBin][iBin2]->Fill( jt, 1.0/jt * effCorrection );
        }
        if(iBin3 > -1){
          fhJetConeJtWeightWithMultiplicityCutBinBin[iContainer][iBin][iBin3]->Fill( jt, 1.0/jt * effCorrection );
        }
        if(fDoLog){
          fhJetConeLogJtWeightBin[iContainer][iBin]
            ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
          fhJetConeLogJtWeight2Bin[iContainer][iBin]
            ->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
        }

        if (iptaBin < 0) continue;
        fhJetConeJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
          ->Fill( jt, 1.0/jt * effCorrection );
        if(fDoLog){
          fhJetConeLogJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
            ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
          fhJetConeLogJtWithPtCutWeight2BinBin[iContainer][iBin][iptaBin]
            ->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
        }
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
            fhBgRBin[iContainer][iBin]->Fill(getDiffR(summedJet.Phi(), track->Phi(), summedJet.Eta(), track->Eta()),effCorrection);
            z = (track->Vect()*summedJet.Vect().Unit())/summedJet.P();
            jt = (track->Vect()-z*summedJet.Vect()).Mag();
          }else{
            fhBgRBin[iContainer][iBin]->Fill(deltaR,effCorrection);
            z = (track->Vect()*vOrtho.Vect().Unit())/vOrtho.P();
            jt = (track->Vect()-z*vOrtho.Vect()).Mag();
          }
          fhBgZ[iContainer]->Fill( z , effCorrection);
          fhBgZBin[iContainer][iBin]->Fill( z , effCorrection);
          fhBgJt[iContainer]->Fill( jt , effCorrection);
          fhBgJtBin[iContainer][iBin]->Fill( jt , effCorrection);
          fhBgJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt * effCorrection );
          if(fDoLog){
            fhBgLogJtWeightBin[iContainer][iBin]
              ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
            fhBgLogJtWeight2Bin[iContainer][iBin]
              ->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
          }

          if (iptaBin < 0) continue;
          fhBgJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
            ->Fill( jt, 1.0/jt * effCorrection );
          if(fDoLog){
            fhBgLogJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
              ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
            fhBgLogJtWithPtCutWeight2BinBin[iContainer][iBin][iptaBin]
              ->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
          }
        }
      }
    }
    fhJetConeMultiplicityBin[iContainer][iBin]->Fill(nC);
    if(doRndmBg){
      this->FillRandomBackground(jet->Pt(), jet->E(),  Jets, iContainer, mc);
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
          if(fDoLog){
            fhBgLogJtWeightBinLimBin[iContainer][iBin][jj]
              ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
            fhBgLogJtWeight2BinLimBin[iContainer][iBin][jj]
              ->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
          }
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



///Fill jet jT with fragmentation information from pythia
///i.e. jT with respect to mother/grandmother/grand-grand... particle
///Fill also the correlation between reconstruced tracks/jets and pythia jet candidates
/// \param Jets List of jets used
/// \param iContainer Index of jet finder
void AliJJetJtAnalysis::FillPythia(TObjArray *Jets, int iContainer){
  int debug = 0;
  if(debug){
    cout << "FillPythia " << endl;
  }
  int iptaBin = 0;
  int iBin = 0;
  double pta = 0;
  double jt = 0;
  double z = 0;
  double deltaR = 0;
  //TVector *jetlist = new TVector(100);
  //TVector *bgtracks= new TVector(1000);
  int jetlist[100] = {0};
  int N_jets = 0;
  int bgtracks[1000] = {0};
  double bgtrackType[1000] = {0};
  int N_bg = 0;
  fpythiaJets.Clear();

  //Loop over Monte Carlo tracks
  for (int icon = 0; icon<fMCTracks->GetEntries(); icon++){
    //AliAODMCParticle *mcTrack = dynamic_cast<AliAODMCParticle*>(fMCTracks->At(icon));
    AliJMCTrack* mcTrack = dynamic_cast<AliJMCTrack*>(fMCTracks->At(icon));
    if (!mcTrack) continue;
    if (!mcTrack->IsPrimary()) continue; //Track must be primary
    if (TMath::Abs(mcTrack->Eta()) > 1.0) continue; //Eta cut
    if (mcTrack->GetCharge() == 0) continue; //Track must be charged
    pta = mcTrack->Pt();
    iptaBin = GetBin(fJetAssocPtBorders, pta);
    if( iptaBin < 0 ) {
      continue;
    }
    //cout << endl << "Track:" << endl;
    //mcTrack->Print();
    // cout << "pT: " << mcTrack->Pt() << "Eta: " << mcTrack->Eta() << "PDG code:" << mcTrack->GetPdgCode() << endl;
    AliJMCTrack *mother0 = dynamic_cast<AliJMCTrack*>(fMCTracks->At(mcTrack->GetMother(0))); //Get the mother particle
    if(!mother0) { //If no mother track is found this is background
      //cout << "Found bg track icon:" << icon << " N_bg: " << N_bg << endl;
      bgtracks[N_bg] = icon; //Add to list of bg tracks
      bgtrackType[N_bg] = 0.5;
      N_bg++;
      continue;
    };
    if(TMath::Abs(mother0->Eta()) > 0.5){ //If the mother track is outside acceptance this is background
      //cout << "Found bg track icon:" << icon << " N_bg: " << N_bg << endl;
      bgtracks[N_bg] = icon; 
      bgtrackType[N_bg] = 1.5;
      N_bg++;
      continue;
    }
    //cout << endl <<"Mother:" << endl;
    //cout << "pT: " << mother0->Pt() << "Eta: " << mother0->Eta() << "PDG code:" << mother0->GetPdgCode() << endl;
    int ijet = FindPythiaJet(iContainer,icon); //Look for the origin of this track
    if(ijet == icon){
      bgtracks[N_bg] = icon; //This was an original track -> Add to background
      bgtrackType[N_bg] = 2.5;
      N_bg++;
      continue; 
    }
    if(ijet < 0){
      if(debug) cout << "ijet: " << ijet << endl;
      //cout << "Found bg track icon:" << icon << " N_bg: " << N_bg << endl;
      bgtracks[N_bg] = icon;
      bgtrackType[N_bg] = 3.5;
      N_bg++;
      continue;
    }else{ //If the origin is found -> calculate jT
      AliJMCTrack *jetCandidate = dynamic_cast<AliJMCTrack*>(fMCTracks->At(ijet));
      double jetPt = jetCandidate->Pt();
      if (TMath::Abs(jetCandidate->Eta()) > fJetEtaCut) continue;
      if(jetPt > (*fJetTriggPtBorders)[1]){ //If the original particle is energetic enough
        iBin = GetBin(fJetTriggPtBorders,jetPt);
        if(iBin < 0) {  
          continue;
        }
        //cout << endl << "Parent Jet: " << endl;
        int found = 0;

        for(int i = 0 ; i < N_jets ; i++){
          if(jetlist[i] == ijet){
            //cout << "Jet already in the list" << endl;
            found = 1;
            AliJJet *j = (AliJJet*)fpythiaJets[i];
            j->AddConstituent(mcTrack);
            break;
          }
        }
        if(found == 0){
          //cout << "New jet" << endl;
          if(N_jets == 10){
            cout << "10 pythia jets already" << endl;
            continue;
          }
          jetlist[N_jets] = ijet;
          new (fpythiaJets[N_jets]) AliJJet(jetCandidate->Px(), jetCandidate->Py(), jetCandidate->Pz(), jetCandidate->E(), jetCandidate->GetLabel(),0,0);
          AliJJet *j = (AliJJet*)fpythiaJets[N_jets];
          j->AddConstituent(mcTrack);
          N_jets++;
        }
        deltaR   = getDiffR(jetCandidate->Phi(),mcTrack->Phi(),jetCandidate->Eta(),mcTrack->Eta());
        if(deltaR > 0.4){ //This could be background to another jet
          bgtracks[N_bg] = icon;
          bgtrackType[N_bg] = 4.5;
          N_bg++;
        }
      }else{ //Origin pT was too small for this to be considered a jet
        //cout << "Found bg track icon:" << icon << " N_bg: " << N_bg << endl;
        bgtracks[N_bg] = icon;  
        bgtrackType[N_bg] = 5.5;
        N_bg++;
        continue;
      }
    }
    for(int ij = 0 ; ij < N_jets ; ij++){
      int found = 0;
      AliJJet *jet= (AliJJet*)fpythiaJets[ij];
      AliJMCTrack *jetCandidate = dynamic_cast<AliJMCTrack*>(fMCTracks->At(jetlist[ij]));
      int jetMult = jet->GetConstituents()->GetEntries();
      if(jetMult < 5) continue;
      for(int ic = 0; ic < jetMult; ic++){
        AliJMCTrack *mcTrack = (AliJMCTrack*)jet->GetConstituent(ic);
        //iptaBin = GetBin(fJetAssocPtBorders, mcTrack->Pt());
        z = (mcTrack->Vect()*jetCandidate->Vect().Unit())/jetCandidate->P();
        fhZPythia->Fill(z);
        fhZBinPythia[iBin]->Fill(z);
        deltaR   = getDiffR(jetCandidate->Phi(),mcTrack->Phi(),jetCandidate->Eta(),mcTrack->Eta());
        fhRBinPythia[iBin]->Fill(deltaR); 
        jt = (mcTrack->Vect()-z*jetCandidate->Vect()).Mag();
        //if(debug) cout         << "Track index: " << icon << " Jet index: " << ijet << " Track pT: " << mcTrack->Pt() << " Eta: " << mcTrack->Eta() << " Phi: " << mcTrack->Phi() << " PDG code:" << mcTrack->GetPdgCode() << " z: " << z << " R: " << deltaR << " jT: " << jt << endl;
        if(z > 0.9 && debug) cout << endl << "Jet   pT: " << jetCandidate->Pt() << " Eta: " << jetCandidate->Eta() << " Phi: " << jetCandidate->Phi() <<  " PDG code:" << jetCandidate->GetPdgCode() << endl;

        fhJtPythia->Fill( jt);
        fhJtBinPythia[iBin]->Fill( jt);
        fhJtWeightBinPythia[iBin]->Fill( jt, 1.0/jt);
        fhJtWithPtCutWeightBinBinPythia[iBin][iptaBin]
          ->Fill( jt, 1.0/jt );
        if(fDoLog){
          fhLogJtWeightBinPythia[iBin]
            ->Fill( TMath::Log(jt), 1.0/jt );
          fhLogJtWeight2BinPythia[iBin]
            ->Fill( TMath::Log(jt), 1.0/jt/jt );
          fhLogJtWithPtCutWeightBinBinPythia[iBin][iptaBin]
            ->Fill( TMath::Log(jt), 1.0/jt );
          fhLogJtWithPtCutWeight2BinBinPythia[iBin][iptaBin]
            ->Fill( TMath::Log(jt), 1.0/jt/jt );
        }
        if(deltaR < 0.4){
          fhJtBinPythiaRCut[iBin]->Fill(jt);
          fhZBinPythiaRCut[iBin]->Fill(jt);
        }
        //Loop over reconstructed tracks to find corresponding tracks
        for(int icon2 = 0 ; icon2 < fTracks->GetEntries() ; icon2++){
          AliJBaseTrack * track = (AliJBaseTrack*) fTracks->At(icon2);
          if(!track) continue;
          if(!track->IsTrue(1)) continue;
          if(mcTrack && (TMath::Abs(track->GetLabel()) == TMath::Abs(mcTrack->GetLabel())) ){
            if((*fTrackJt)[icon2] < -0){
              if(debug) cout << "Reconstructed track jT not measured, MC jT was " << jt << endl;
              found = -4;
              fhTrackJtCorrBinPythia[iBin]->Fill(-1.0,jt);
              fhTrackPtCorrPythia->Fill((*fTrackPt)[icon2],mcTrack->Pt());
              fhJetPtCorr2Pythia->Fill(-1.0,jetCandidate->Pt());
              continue;
            }
            found++;
            if((*fTrackJt)[icon2] < 0){
              fhTrackJtCorrBinPythia[iBin]->Fill(-1,jt);
              if(debug) cout << "Reconstructed track jT not measured, MC jT was " << jt << endl;
            }else{
              fhTrackJtCorrBinPythia[iBin]->Fill((*fTrackJt)[icon2],jt);
              if(debug) cout << "Reconstructed track measured: " << (*fTrackJt)[icon2] << " , MC jT was " << jt << endl;
            }
            fhTrackPtCorrPythia->Fill((*fTrackPt)[icon2],mcTrack->Pt());
            fhJetPtCorr2Pythia->Fill((*fJetPt)[icon2],jetCandidate->Pt());
          }
        }
        if(found == 0){
          if(debug) cout << "FOUND NO MATCHING TRACK" << endl;
          fhTrackMatchSuccessPythia[iBin]->Fill(0.5);
          fhTrackJtCorrBinPythia[iBin]->Fill(-1.0,jt);
          fhJetPtCorr2Pythia->Fill(-1.0,jetCandidate->Pt());
          fhTrackPtCorrPythia->Fill(-1.0,mcTrack->Pt());
        }
        if(found > 2){
          fhTrackMatchSuccessPythia[iBin]->Fill(2.5);
          if(debug) cout << "FOUND MORE THAN 1 MATCHING TRACK" << endl;
        }
        if(found == 2){
          fhTrackMatchSuccessPythia[iBin]->Fill(4.5);
          if(debug) cout << "FOUND EXACTLY 2 MATCHING TRACKS" << endl;
        }
        if(found < 0){
          fhTrackMatchSuccessPythia[iBin]->Fill(3.5);
          if(debug) cout << "FOUND MATCHING TRACK WITHOUT JT" << endl;
        }
        if(found == 1){
          fhTrackMatchSuccessPythia[iBin]->Fill(1.5);
          if(debug) cout << "FOUND EXACTLY 1 MATCHING TRACK" << endl;
        }
      }
    }
  }
  if(debug) cout << "N_bg: " << N_bg << endl; //Number of background particles
  //Loop over background tracks
  for(int i = 0; i < N_bg; i++){
    //cout << "Bg track i: " << i << " icon: " << bgtracks[i] << endl;
    AliJMCTrack* mcTrack = dynamic_cast<AliJMCTrack*>(fMCTracks->At(bgtracks[i]));
    pta = mcTrack->Pt();
    iptaBin = GetBin(fJetAssocPtBorders, pta);
    if( iptaBin < 0 ) continue;
    //Loop over jets
    for(int ij = 0 ; ij < N_jets; ij++){
      if(ij == i) continue;
      AliJMCTrack *jetCandidate = dynamic_cast<AliJMCTrack*>(fMCTracks->At(jetlist[ij]));
      deltaR   = getDiffR(jetCandidate->Phi(),mcTrack->Phi(),jetCandidate->Eta(),mcTrack->Eta());
      //Check if the background track is close to the jet 
      if(deltaR < 0.4){
        z = (mcTrack->Vect()*jetCandidate->Vect().Unit())/jetCandidate->P();
        jt = (mcTrack->Vect()-z*jetCandidate->Vect()).Mag();
        //cout << endl << "Jet   pT: " << jetCandidate->Pt() << " Eta: " << jetCandidate->Eta() << " Phi: " << jetCandidate->Phi() <<  " PDG code:" << jetCandidate->GetPdgCode() << endl;
        if(debug) cout         << "Background Track index: " << i << " pT: " << mcTrack->Pt() << " Eta: " << mcTrack->Eta() << " Phi: " << mcTrack->Phi() << " PDG code:" << mcTrack->GetPdgCode() << " z: " << z << " R: " << deltaR << " jT: " << jt << endl;
        if(z > 0.99 && debug)  cout << endl << "Jet index: " << ij << " pT: " << jetCandidate->Pt() << " Eta: " << jetCandidate->Eta() << " Phi: " << jetCandidate->Phi() <<  " PDG code:" << jetCandidate->GetPdgCode() << endl;
        fhBgJtPythiaTrackType->Fill(bgtrackType[i]);
        fhBgJtBinPythiaTypeBin[iBin][bgtrackType[i]]->Fill(jt, 1.0/jt);
        fhBgJtPythia->Fill( jt);
        fhBgJtBinPythia[iBin]->Fill( jt);
        fhBgJtWeightBinPythia[iBin]->Fill( jt, 1.0/jt);
        fhBgJtWithPtCutWeightBinBinPythia[iBin][iptaBin]
          ->Fill( jt, 1.0/jt );
        if(fDoLog){
          fhBgLogJtWeightBinPythia[iBin]
            ->Fill( TMath::Log(jt), 1.0/jt );
          fhBgLogJtWeight2BinPythia[iBin]
            ->Fill( TMath::Log(jt), 1.0/jt/jt );
          fhBgLogJtWithPtCutWeightBinBinPythia[iBin][iptaBin]
            ->Fill( TMath::Log(jt), 1.0/jt );
          fhBgLogJtWithPtCutWeight2BinBinPythia[iBin][iptaBin]
            ->Fill( TMath::Log(jt), 1.0/jt/jt );
        }
      }
    }
  }



  int found = 0;
  double pT = 0;
  double pTmc = 0;
  //Loop over jet candidates to get correlation between them and reco jets
  //cout << "Number of pythia jet candidates: " << N_jets << endl;
  //TODO Move to jet loop?
  for(int ij = 0 ; ij < N_jets; ij++){
    AliJJet *pythiaJet = (AliJJet*)fpythiaJets[ij];
    AliJMCTrack *jetCandidate = dynamic_cast<AliJMCTrack*>(fMCTracks->At(jetlist[ij]));
    if (!jetCandidate) continue;
    if (TMath::Abs(jetCandidate->Eta()) > fJetEtaCut) continue;
    pTmc = jetCandidate->Pt();
    if (pTmc<(*fJetTriggPtBorders)[1]) continue;
    iBin = GetBin(fJetTriggPtBorders,pTmc);
    if( iBin < 0 ) continue;
    fhJetPtPythia->Fill(pTmc);
    fhJetPtBinPythia[iBin]->Fill(pTmc);
    int jetMult = pythiaJet->GetConstituents()->GetEntries();
    fhJetMultiplicityBinPythia[iBin]->Fill(jetMult);
    found = 0;
    //Look for matching reco jet
    for( int j = 0 ; j < Jets->GetEntries(); j++){
      AliJJet *jet = dynamic_cast<AliJJet*>( Jets->At(j));
      if(!jet) continue;
      if (TMath::Abs(jet->Eta()) > fJetEtaCut) continue;
      pT = jet->Pt();
      deltaR   = getDiffR(jet->Phi(),jetCandidate->Phi(),jet->Eta(),jetCandidate->Eta());
      if(deltaR < 0.5){
        fhJetPtCorrPythia->Fill(pTmc,pT);
        fhJetPtCorrPythiaCoarse->Fill(pTmc,pT);
        fhJetdRPythia[iBin]->Fill(deltaR);
        fhJetdPtPythia[iBin]->Fill(pTmc-pT);
        found = 1;
      }
      if(found == 0){
        fhJetPtCorrPythia->Fill(pTmc,-1);
        fhJetPtCorrPythiaCoarse->Fill(pTmc,-1);
      }
    }	
  }
}

/// Look for original particle of primary track
/// Returns -1 if particle is outside acceptance
/// Returns the track index if this is considered an original particle, i.e. mother outside acceptance or too far from this track
/// Otherwise return the index of mother particle
/// \param iContainer Index of jet finder, not needed!
/// \param itrack Index of track
/// \return Index of the mother particle


//FIXME iContainer not needed
int AliJJetJtAnalysis::FindPythiaJet(int iContainer, int itrack){
  //cout << "itrack: " << itrack << endl;
  AliJMCTrack *track= dynamic_cast<AliJMCTrack*>(fMCTracks->At(itrack));
  //if(!track) return -2;
  if(TMath::Abs(track->Eta()) > 1.0) return -1; //Track is outside acceptance -> return -1
  int iMother = track->GetMother(0); //Get the index of mother particle
  if(iMother == -1) return -1; //No mother was found, unlikely?
  int ijet = FindPythiaJet(iContainer,iMother); //Track is within acceptance -> move to the mother particle
  if(ijet == -1){ //If the mother particle was outside acceptance, this is considered the original particle
    return itrack;
  }else{ //If mother is found
    AliJMCTrack *jetCandidate = dynamic_cast<AliJMCTrack*>(fMCTracks->At(ijet)); 
    double deltaR   = getDiffR(jetCandidate->Phi(),track->Phi(),jetCandidate->Eta(),track->Eta()); //Check the separation of mother and daughter
    if(deltaR > 0.8){ //If the separation is too large this is considered an original track
      //cout << "deltaR: " << deltaR << endl;
      return itrack; 
    }
    return ijet; //Otherwise return the index of the mother particle
  }
}


/// Fill correlation histograms between particle and detector level for jT and jet pT
/// \param Jets List of jets to be used
/// \param MCJets List of corresponding MC jets
/// \param Index of detector level jet finder
/// \param Index of particle level jet finder
void AliJJetJtAnalysis::FillCorrelation(TObjArray *Jets, TObjArray *MCJets, int iContainer, int iContainerParticle){

  int iBin=0;
  double pT = 0;
  double pTmc = 0;
  double z; double jt;
  double pta;
  //double Y , deltaY = 0;
  //double deltaR= 0;
  //cout<<"histogram filling number of jets : "<<Jets->GetEntriesFast()<<endl;

  TLorentzVector  vOrtho;
  TLorentzVector  randomTrack;
  TLorentzVector summedJet;
  fJetBgListOfList[iContainer].Clear();

  //cout << "REMOVE THIS " << endl << endl << endl;


  double deltaR = -1;
  double thisConeSize = fConeSizes[iContainer] ;
  int found = 0;

  // Loop over Monte Carlo jets to get correlation between true and reco jet pT
  for (int i = 0; i<MCJets->GetEntries(); i++){
    if(i > 0 && fLeadingJets > 0) continue;
    AliJJet *mcjet = dynamic_cast<AliJJet*>( MCJets->At(i) );
    if (!mcjet) continue;
    if (TMath::Abs(mcjet->Eta()) > fJetEtaCut) continue;
    pTmc = mcjet->Pt();
    if (pTmc<(*fJetTriggPtBorders)[1]) continue;
    iBin = GetBin(fJetTriggPtBorders,pTmc); // fill jetPt histos	
    if( iBin < 0 ) continue;
    found = 0;
    for( int j = 0 ; j < Jets->GetEntries(); j++){
      AliJJet *jet = dynamic_cast<AliJJet*>( Jets->At(j));
      if(!jet) continue;
      if (TMath::Abs(jet->Eta()) > fJetEtaCut) continue;
      pT = jet->Pt();
      deltaR   = getDiffR(jet->Phi(),mcjet->Phi(),jet->Eta(),mcjet->Eta());
      if(deltaR < 0.5){
        fhJetPtCorr[iContainer]->Fill(pTmc,pT);
        fhJetPtCorrCoarse[iContainer]->Fill(pTmc,pT);
        fhJetdR[iContainer][iBin]->Fill(deltaR);
        fhJetdPt[iContainer][iBin]->Fill(pTmc-pT);
        found = 1;
      }
      if(found == 0){
        fhJetPtCorr[iContainer]->Fill(pTmc,-1);
        fhJetPtCorrCoarse[iContainer]->Fill(pTmc,-1);
      }
    }	
  }

  int iptaBin = 0;
  int iBinDet = 0;

  //Loop over Monte Carlo jets
  for (int i = 0 ; i<MCJets->GetEntries(); i++){
    AliJJet *mcjet = dynamic_cast<AliJJet*>( MCJets->At(i));
    if(!mcjet) continue;
    if (TMath::Abs(mcjet->Eta()) > fJetEtaCut) continue;
    pTmc = mcjet->Pt();
    if (pTmc<(*fJetTriggPtBorders)[1]) continue;
    iBin = GetBin(fJetTriggPtBorders,pTmc); // fill jetPt histos	
    if(iBin < 0) continue;
    for (int icon = 0; icon<mcjet->GetConstituents()->GetEntries(); icon++){
      AliJMCTrack *constituent = dynamic_cast<AliJMCTrack*>(mcjet->GetConstituent(icon));
      if (!constituent) continue;
      if (!constituent->IsPrimary()) continue;
      if (TMath::Abs(constituent->Eta()) > 1.0) continue;
      //if (constituent->GetCharge() == 0) continue;
      /*if (constituent->GetCharge() == 0) cout << "Constituent charge = 0" << endl;
        else cout << "Constituent charge not 0" << endl;*/

      z = (constituent->Vect()*mcjet->Vect().Unit())/mcjet->P();
      pta = constituent->Pt();
      //if(iptaBin < 0 ) continue;
      jt = (constituent->Vect()-z*mcjet->Vect()).Mag();
      iptaBin = GetBin(fJetAssocPtBorders, pta);
      int found = 0;

      for(int icon2 = 0 ; icon2 < iConst; icon2++){
        //iBinDet = GetBin(fJetTriggPtBorders,(*fJetPtConst)[icon2]);//TODO Create
        if(constituent && (TMath::Abs((*fConstLabels)[icon2]) == TMath::Abs(constituent->GetLabel())) ){ //If the track labels match
          //cout << "Found match" << endl;
          (*fConstFound)[icon2] = 1;
          found++;
          //fhConstJtCorrBin[iContainer][iBinDet]->Fill((*fConstJt)[icon2],jt); //FIXME
          fhConstJtCorrBin[iContainer][iBin]->Fill((*fConstJt)[icon2],jt);
          fhConstPtCorr[iContainer]->Fill((*fConstPt)[icon2],constituent->Pt());
          fhJetPtCorr3[iContainer]->Fill((*fJetPt)[icon2],mcjet->Pt());
        }
      }
      if(found == 0){ //Found no match
        //cout << "No Match found" << endl;
        fhConstMatchSuccess[iContainer][iBin]->Fill(0.5);
        fhConstJtCorrBin[iContainer][iBin]->Fill(-1.0,jt);
        fhJetPtCorr3[iContainer]->Fill(-1.0,mcjet->Pt());
        fhConstPtCorr[iContainer]->Fill(-1.0,constituent->Pt());
      }
      if(found > 2){ //Found more than 2 matching tracks
        fhConstMatchSuccess[iContainer][iBin]->Fill(2.5);
      }
      if(found == 2){ //Found exactly 2 matches
        fhConstMatchSuccess[iContainer][iBin]->Fill(4.5);
      }
      if(found == -1){ //Found match but jt was not measured
        fhConstMatchSuccess[iContainer][iBin]->Fill(3.5);
      }
      if(found == 1){ //One match was found
        fhConstMatchSuccess[iContainer][iBin]->Fill(1.5);
      }
    }


    for (int icon = 0; icon<fMCTracks->GetEntries(); icon++){
      AliJMCTrack *mcTrack = dynamic_cast<AliJMCTrack*>(fMCTracks->At(icon));
      if (!mcTrack) continue;
      if (!mcTrack->IsPrimary()) continue;
      if (TMath::Abs(mcTrack->Eta()) > 1.0) continue;
      if (mcTrack->GetCharge() == 0) continue;
      pta = mcTrack->Pt();
      iptaBin = GetBin(fJetAssocPtBorders, pta);
      if( iptaBin < 0 ) continue;

      double z = (mcTrack->Vect()*mcjet->Vect().Unit())/mcjet->P();
      //Jet Cone Jt here
      deltaR   = getDiffR(mcjet->Phi(),mcTrack->Phi(),mcjet->Eta(),mcTrack->Eta());
      double jt = 0;
      if ( deltaR < thisConeSize){ 
        jt = (mcTrack->Vect()-z*mcjet->Vect()).Mag();
        int found = 0;
        //AliJMCTrack *mcTrack = (AliJMCTrack*)fMCTracks->At(icon);
        if(!mcTrack || !mcTrack->IsPrimary()){
          cout << "Not primary " << endl;
          continue;
        }
        fhJtWeightBinTest[iContainerParticle][iBin]->Fill(jt,1.0/jt); //TODO Document
        //cout << " icon: " << icon << ", label: " << mcTrack->GetLabel() << " Eta: " << mcTrack->Eta() << " Phi: " << mcTrack->Phi() << ", pT: " << mcTrack->Pt() << ", jT: " << jt << " Charge: " << mcTrack->GetCharge()  << endl;


        for(int icon2 = 0 ; icon2 < fTracks->GetEntries() ; icon2++){
          AliJBaseTrack * track = (AliJBaseTrack*) fTracks->At(icon2);
          if(!track) continue;
          if(!track->IsTrue(1)) continue;
          iBinDet = GetBin(fJetTriggPtBorders,(*fJetPt)[icon2]); // fill jetPt histos	
          if(iBinDet < 0) continue;
          if(mcTrack && (TMath::Abs(track->GetLabel()) == TMath::Abs(mcTrack->GetLabel())) ){ //If the track labels match
            /*cout << "Found match: " << endl;
              cout << "Reco track index: " << icon2             << " mc Track index:  " << icon << endl;
              cout << "Reco track Label: " << track->GetLabel() << " mc Track Label:  " << mcTrack->GetLabel()  << endl;
              cout << "Reco track Eta:   " << track->Eta() <<      " mc Track Eta:    " << mcTrack->Eta()  << endl;
              cout << "Reco track Phi:   " << track->Phi() <<      " mc Track Eta:    " << mcTrack->Phi()  << endl;
              cout << "Reco track pT:    " << track->Pt() <<       " mc Track pT:     " << mcTrack->Pt()  << endl;
              */
            (*fTrackFound)[icon2] = 1;
            if((*fTrackJt)[icon2] < -0){ //Matching reco track was found but it has no jT value
              //cout << "Reconstructed track jT not measured, MC jT was " << jt << endl;
              found = -4;
              fhJtWeightBinTest2[iContainerParticle][iBin]->Fill(jt,1.0/jt); //TODO Document
              //cout << "fhJtWeightBinTest2 Filled, match found, but no jt";
              //cout << " icon2: " << icon2 << ", label: " << track->GetLabel() << " Eta: " << track->Eta() << " Phi: " << track->Phi() << ", pT: " << track->Pt() << " Charge: " << track->GetCharge() <<endl;

              fhTrackJtCorrBin[iContainer][iBin]->Fill(-1.0,jt);
              //fhTrackJtCorr2D[iContainer][iBin]->Fill(-1.0,jt,-1.0);
              fhTrackJtCorrBinTest[iContainer][iBin]->Fill(-1.0,jt);
              fhTrackPtCorr[iContainer]->Fill((*fTrackPt)[icon2],mcTrack->Pt());
              fhJetPtCorr2[iContainer]->Fill(-1.0,mcjet->Pt());
              continue;
            }
            found++;
            if((*fTrackJt)[icon2] < 0){
              fhTrackJtCorrBin[iContainer][iBin]->Fill(-1,jt);
              //fhTrackJtCorr2D[iContainer][iBin]->Fill(-1.0,jt,-1.0);
              fhTrackJtCorrBinTest[iContainer][iBin]->Fill(-1.0,jt);
              fhJtWeightBinTest2[iContainerParticle][iBin]->Fill(jt,1.0/jt); //TODO Document
              cout << "fhJtWeightBinTest2 Filled, Impossible!!!!" << endl;
              //fhJtWeightBinTest[iContainerParticle][iBin]->Fill(jt);
            }else{
              fhTrackJtCorrBin[iContainer][iBinDet]->Fill((*fTrackJt)[icon2],jt);
              //fhTrackJtCorr2D[iContainer][iBin]->Fill((*fTrackJt)[icon2],jt,(*fJetPt)[icon2]);
              fhTrackJtCorrBinTest[iContainer][iBinDet]->Fill((*fTrackJt)[icon2],jt);
              fhTrackJtCorrBinTest2[iContainer][iBinDet]->Fill((*fTrackJt)[icon2],jt);
              fhJtWeightBinTest[iContainer][iBinDet]->Fill((*fTrackJt)[icon2],1.0/(*fTrackJt)[icon2]); //TODO Document
              //fhJtWeightBinTest[iContainer+6][iBinDet]->Fill(jt,1.0/jt);
              fhJtWeightBinTest2[iContainer][iBinDet]->Fill((*fTrackJt)[icon2],1.0/(*fTrackJt)[icon2]); //TODO Document
              fhJtWeightBinTest2[iContainerParticle][iBinDet]->Fill(jt,1.0/jt); //TODO Document
            }
            fhTrackPtCorr[iContainer]->Fill((*fTrackPt)[icon2],mcTrack->Pt());
            fhJetPtCorr2[iContainer]->Fill((*fJetPt)[icon2],mcjet->Pt());
          }
        }
        if(found == 0){ //Found no match
          fhTrackMatchSuccess[iContainer][iBin]->Fill(0.5);
          fhTrackJtCorrBin[iContainer][iBin]->Fill(-1.0,jt);
          //fhTrackJtCorr2D[iContainer][iBin]->Fill(-1.0,jt,-1.0);
          fhTrackJtCorrBinTest[iContainer][iBin]->Fill(-1.0,jt);
          fhJtWeightBinTest2[iContainerParticle][iBin]->Fill(jt,1.0/jt); //TODO Document
          fhJetPtCorr2[iContainer]->Fill(-1.0,mcjet->Pt());
          fhTrackPtCorr[iContainer]->Fill(-1.0,mcTrack->Pt());
        }
        if(found > 2){ //Found more than 2 matching tracks
          fhTrackMatchSuccess[iContainer][iBin]->Fill(2.5);
        }
        if(found == 2){ //Found exactly 2 matches
          fhTrackMatchSuccess[iContainer][iBin]->Fill(4.5);
        }
        if(found == -1){ //Found match but jt was not measured
          fhTrackMatchSuccess[iContainer][iBin]->Fill(3.5);
        }
        if(found == 1){ //One match was found
          fhTrackMatchSuccess[iContainer][iBin]->Fill(1.5);
        }
      }
    }
  }
  double effCorrection;
  for(int icon2 = 0 ; icon2 < fTracks->GetEntries() ; icon2++){
    iBinDet = GetBin(fJetTriggPtBorders,(*fJetPt)[icon2]); // fill jetPt histos	
    if(iBinDet < 0) continue;
    AliJBaseTrack * track = (AliJBaseTrack*) fTracks->At(icon2);
    if(!track) continue;
    if(!track->IsTrue(1)) continue;
    if((*fTrackFound)[icon2] < 0 && (*fTrackJt)[icon2] > 0){
      fhTrackMatchSuccess[iContainer][iBinDet]->Fill(5.5); //Track with constructed j_T but no pythia counterpart
      //cout << "Constructed jT: " << (*fTrackJt)[icon2] << " with no pythia counterpart" << endl;
      pta = (*fTrackPt)[icon2];
      //track->SetTrackEff( fEfficiency->GetCorrection( pta, 5, fcent) );
      effCorrection = 1.0/fEfficiency->GetCorrection(pta,5,fcent);
      jt = (*fTrackJt)[icon2];
      fhJetConeJtUnfBg[iContainer]->Fill( jt , effCorrection);
      fhJetConeJtBinUnfBg[iContainer][iBinDet]->Fill( jt , effCorrection);
      fhJetConeJtWeightBinUnfBg[iContainer][iBinDet]->Fill( jt, 1.0/jt * effCorrection );
      if(fDoLog){
        fhJetConeLogJtWeightBinUnfBg[iContainer][iBinDet]
          ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
        fhJetConeLogJtWeight2BinUnfBg[iContainer][iBinDet]
          ->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
      }

      /*if (iptaBin < 0) continue;
        fhJetConeJtWithPtCutWeightBinBinUnfBg[iContainer][iBin][iptaBin]
        ->Fill( jt, 1.0/jt * effCorrection );
        fhJetConeLogJtWithPtCutWeightBinBinUnfBg[iContainer][iBin][iptaBin]
        ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
        fhJetConeLogJtWithPtCutWeight2BinBinUnfBg[iContainer][iBin][iptaBin]
        ->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );*/
    }
  }
}

/// Fill random background histograms
/// \param jetpT Jet transverse momentum
/// \param jetE Jet energy
/// \param Jets List of jets used
/// \param iContainer Index of jet finder
/// \param mc Whether or not this is particle level data

void AliJJetJtAnalysis::FillRandomBackground(double jetpT, double jetE, TObjArray *Jets , int iContainer, int mc){
  TClonesArray *trackArray;
  if(mc){
    trackArray = fMCTracks;
  }
  else 
    trackArray = fTracks;
  int iBin, iptaBin=0;
  iBin = GetBin(fJetTriggPtBorders,jetpT); // fill jetPt histos
  if(iBin < 0 ) return;



  double z; double jt;
  double pta;
  //double Y , deltaY = 0;
  //double deltaR= 0;
  //cout<<"histogram filling number of jets : "<<Jets->GetEntriesFast()<<endl;

  TLorentzVector  vRandomJet;
  TLorentzVector  randomTrack;
  TLorentzVector summedJet;

  int counter =0 ;
  double deltaR = -1;
  double effCorrection = -1;
  double thisConeSize = fConeSizes[iContainer] ;
  double etaRndm = 0;
  double phiRndm = 0;
  etaRndm  = frandom->Uniform(-1.0+thisConeSize, 1.0-thisConeSize); //TODO Fix with realistic eta distribution
  phiRndm = kJPi * frandom->Uniform(-1, 1); //TODO realistic phi distribution
  vRandomJet.SetPtEtaPhiE(jetpT,etaRndm,phiRndm,jetE);

  int countTrack = 1;
  fhNumber[iContainer]->Fill(6.5+iBin,Nrandom);
  for (int icon = 0; icon<trackArray->GetEntries(); icon++){
    //cout << "icon: " << icon << "/" << trackArray->GetEntries() << endl;
    countTrack = 1;
    AliJBaseTrack *track = dynamic_cast<AliJBaseTrack*>(trackArray->At(icon));
    if (!track) continue;
    if(!mc){
      if (!track->IsTrue(1)) continue;
    }
    if (Jets->GetEntries()>0){
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
    pta = track->Pt();
    if(!mc){
      track->SetTrackEff( fEfficiency->GetCorrection( pta, 5, fcent) );
      effCorrection = 1.0/track->GetTrackEff();
    }else{
      effCorrection = 1.0;
    }
    iptaBin = GetBin(fJetAssocPtBorders, pta);
    if( iptaBin < 0 ) continue;
    if(track->Eta()*track->Eta() < 1){
      //AliJBaseTrack *randomTrack = dynamic_cast<AliJBaseTrack*>(trackArray->At(icon));
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
            z = (randomTrack.Vect()*summedJet.Vect().Unit())/summedJet.P();
            fhBgRndmRBin[iContainer][iBin]->Fill(getDiffR(summedJet.Phi(), track->Phi(), summedJet.Eta(), track->Eta()),effCorrection);
            jt = (randomTrack.Vect()-z*summedJet.Vect()).Mag();
          }
          else{
            z = (randomTrack.Vect()*vRandomJet.Vect().Unit())/vRandomJet.P();
            jt = (randomTrack.Vect()-z*vRandomJet.Vect()).Mag();
            fhBgRndmRBin[iContainer][iBin]->Fill(deltaR,effCorrection);
          }
          fhBgRndmTrkPt[iContainer]->Fill(pta,effCorrection);
          fhBgRndmZ[iContainer]->Fill( z , effCorrection);
          fhBgRndmJt[iContainer]->Fill( jt , effCorrection);
          fhBgRndmJtWeight[iContainer]->Fill( jt , 1.0/jt* effCorrection);
          fhBgRndmJtBin[iContainer][iBin]->Fill(jt, effCorrection);
          fhBgRndmJtWeightBin[iContainer][iBin]->Fill(jt, 1.0/jt * effCorrection); 
          if(fDoLog) fhBgRndmLogJt[iContainer]->Fill( jt , 1.0/jt/jt * effCorrection);
          if (iptaBin < 0) continue;
          fhBgRndmJtWithPtCutWeightBin[iContainer][iptaBin]->Fill( jt, 1.0/jt * effCorrection );
          fhBgRndmJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]->Fill( jt, 1.0/jt * effCorrection );
          if(fDoLog){
            fhBgRndmLogJtWithPtCutWeight2Bin[iContainer][iptaBin]->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
            fhBgRndmLogJtWithPtCutWeight2BinBin[iContainer][iBin][iptaBin]->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
          }
        }
      }
    }
  }
  fhBgRndmTrkNumber[iContainer]->Fill(counter/Nrandom);
}

/// new Bg jt will be filled with a new cone size nR at histograms in iHist
/// cone size of Jets always shold be greather than nR 
/// to calcualte new bg jt with smaller cone size from Jets constituents with larger cone size
/// \param Jets List of jets
/// \param R new cone size
/// \param iHist Index of jet finder
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
      if(fDoLog){
        fhBgLogJtWithPtCutWeightBinBinSmallerR[iHist][iBin][iptaBin]
          ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
        fhBgLogJtWithPtCutWeight2BinBinSmallerR[iHist][iBin][iptaBin]
          ->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
      }
    }


  }


}

/// bg jt is newly calculated by a jet axis of fJetBgListOfList[ian]
/// with bg constituents from fJetBgListOfList[iao]
/// \param iao Comment needed
/// \param ian Comment needed
/// \param iHist Comment needed
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
        if(fDoMC){
          fhBgLogJtWithPtCutWeightBinBinDiffR[iHist][iBin][iptaBin]
            ->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
          fhBgLogJtWithPtCutWeight2BinBinDiffR[iHist][iBin][iptaBin]
            ->Fill( TMath::Log(jt), 1.0/jt/jt * effCorrection );
        }
      }
    }
  }

}

//Phi1 and Phi2 between 0 and 2 pi
///Get distance between two tracks or jet and track as \f$ \Delta R = \sqrt{\left( \Delta \phi \right)^2 + \left(\Delta \eta\right)^2}   \f$
/// 
/// 
Double_t AliJJetJtAnalysis::getDiffR(double phi1, double phi2, double eta1, double eta2){
  Double_t diffPhi = TMath::Abs(phi1-phi2);
  if(diffPhi > TMath::Pi()){
    diffPhi = 2*TMath::Pi() - diffPhi;
  }
  return TMath::Sqrt(TMath::Power(diffPhi,2)+TMath::Power(eta1-eta2,2));
}


