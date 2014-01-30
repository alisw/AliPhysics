/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//_________________________________________________________________________
//
// Split clusters with some criteria and calculate invariant mass
// to identify them as pi0 or conversion
//
//
//-- Author: Gustavo Conesa (LPSC-Grenoble)  
//_________________________________________________________________________

//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TList.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TH2F.h>
#include <TDatabasePDG.h>

// --- Analysis system --- 
#include "AliAnaInsideClusterInvariantMass.h" 
#include "AliCaloTrackReader.h"
#include "AliMCAnalysisUtils.h"
#include "AliStack.h"
#include "AliFiducialCut.h"
#include "TParticle.h"
#include "AliVCluster.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliEMCALGeoParams.h"

// --- Detectors ---
//#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"

ClassImp(AliAnaInsideClusterInvariantMass)
  
//__________________________________________________________________
AliAnaInsideClusterInvariantMass::AliAnaInsideClusterInvariantMass() : 
  AliAnaCaloTrackCorrBaseClass(),            fCalorimeter(""),
  fMinNCells(0),                             fMinBadDist(0),
  fHistoECut(0),                             fCheckSplitDistToBad(0),                   fFillAngleHisto(kFALSE),
  fFillTMHisto(kFALSE),                      fFillTMResidualHisto(kFALSE),              fFillSSExtraHisto(kFALSE),
  fFillMCHisto(kFALSE),                      fFillSSWeightHisto(kFALSE),                fFillEbinHisto(0),
  fFillMCOverlapHisto(0),                    fFillNCellHisto(0),                        fFillIdConvHisto(0),
  fFillIdEtaHisto(0),                        fFillHighMultHisto(0),
  fFillArmenterosHisto(0),                   fFillThetaStarHisto(0),
  fSSWeightN(0),                             fSSECellCutN(0),                           fWSimu(0),
  fhMassAsyCutNLocMax1(0),                   fhMassAsyCutNLocMax2(0),                   fhMassAsyCutNLocMaxN(0),
  fhM02AsyCutNLocMax1(0),                    fhM02AsyCutNLocMax2(0),                    fhM02AsyCutNLocMaxN(0),
  fhMassM02CutNLocMax1(0),                   fhMassM02CutNLocMax2(0),                   fhMassM02CutNLocMaxN(0),
  fhAsymM02CutNLocMax1(0),                   fhAsymM02CutNLocMax2(0),                   fhAsymM02CutNLocMaxN(0),
  fhMassEnCutNLocMax1(0),                    fhMassEnCutNLocMax2(0),                    fhMassEnCutNLocMaxN(0),
  fhM02EnCutNLocMax1(0),                     fhM02EnCutNLocMax2(0),                     fhM02EnCutNLocMaxN(0),
  fhAsymEnCutNLocMax1(0),                    fhAsymEnCutNLocMax2(0),                    fhAsymEnCutNLocMaxN(0),
  fhSplitEFracEnCutNLocMax1(0),              fhSplitEFracEnCutNLocMax2(0),              fhSplitEFracEnCutNLocMaxN(0),
  fhMassSplitECutNLocMax1(0),                fhMassSplitECutNLocMax2(0),                fhMassSplitECutNLocMaxN(0),
  fhMCGenFracAfterCutsNLocMax1MCPi0(0),      fhMCGenFracAfterCutsNLocMax2MCPi0(0),      fhMCGenFracAfterCutsNLocMaxNMCPi0(0),
  fhMCGenSplitEFracAfterCutsNLocMax1MCPi0(0),fhMCGenSplitEFracAfterCutsNLocMax2MCPi0(0),fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0(0),
  fhNCellMassEHighNLocMax1MCPi0(0),          fhNCellM02EHighNLocMax1MCPi0(0),
  fhNCellMassELowNLocMax1MCPi0(0),           fhNCellM02ELowNLocMax1MCPi0(0),
  fhNCellMassEHighNLocMax2MCPi0(0),          fhNCellM02EHighNLocMax2MCPi0(0),
  fhNCellMassELowNLocMax2MCPi0(0),           fhNCellM02ELowNLocMax2MCPi0(0),
  fhNCellMassEHighNLocMaxNMCPi0(0),          fhNCellM02EHighNLocMaxNMCPi0(0),
  fhNCellMassELowNLocMaxNMCPi0(0),           fhNCellM02ELowNLocMaxNMCPi0(0),
  fhAnglePairPrimPi0RecoNLocMax1(0),         fhAnglePairPrimPi0RecoNLocMax2(0),         fhAnglePairPrimPi0RecoNLocMaxN(0),
  fhAnglePairPrimPi0vsRecoNLocMax1(0),       fhAnglePairPrimPi0vsRecoNLocMax2(0),       fhAnglePairPrimPi0vsRecoNLocMaxN(0),
  fhAnglePairPrimPi0OverM02NLocMax1(0),      fhAnglePairPrimPi0OverM02NLocMax2(0),      fhAnglePairPrimPi0OverM02NLocMaxN(0),
  fhCentralityPi0NLocMax1(0),                fhCentralityEtaNLocMax1(0),
  fhCentralityPi0NLocMax2(0),                fhCentralityEtaNLocMax2(0),
  fhCentralityPi0NLocMaxN(0),                fhCentralityEtaNLocMaxN(0),
  fhEventPlanePi0NLocMax1(0),                fhEventPlaneEtaNLocMax1(0),
  fhEventPlanePi0NLocMax2(0),                fhEventPlaneEtaNLocMax2(0),
  fhEventPlanePi0NLocMaxN(0),                fhEventPlaneEtaNLocMaxN(0),
  fhClusterEtaPhiNLocMax1(0),                fhClusterEtaPhiNLocMax2(0),                fhClusterEtaPhiNLocMaxN(0),
  fhPi0EtaPhiNLocMax1(0),                    fhPi0EtaPhiNLocMax2(0),                    fhPi0EtaPhiNLocMaxN(0),
  fhEtaEtaPhiNLocMax1(0),                    fhEtaEtaPhiNLocMax2(0),                    fhEtaEtaPhiNLocMaxN(0),
  fhPi0EPairDiffTimeNLM1(0),                 fhPi0EPairDiffTimeNLM2(0),                 fhPi0EPairDiffTimeNLMN(0),
  fhEtaEPairDiffTimeNLM1(0),                 fhEtaEPairDiffTimeNLM2(0),                 fhEtaEPairDiffTimeNLMN(0),
  fhMCPi0HighNLMPair(0),                     fhMCPi0LowNLMPair(0),
  fhMCPi0AnyNLMPair(0),                      fhMCPi0NoneNLMPair(0),
  fhMCPi0HighNLMPairNoMCMatch(0),            fhMCPi0LowNLMPairNoMCMatch(0),
  fhMCPi0AnyNLMPairNoMCMatch(0),             fhMCPi0NoneNLMPairNoMCMatch(0),
  fhMCPi0HighNLMPairOverlap(0),              fhMCPi0LowNLMPairOverlap(0),
  fhMCPi0AnyNLMPairOverlap(0),               fhMCPi0NoneNLMPairOverlap(0),
  fhMCPi0HighNLMPairNoMCMatchOverlap(0),     fhMCPi0LowNLMPairNoMCMatchOverlap(0),
  fhMCPi0AnyNLMPairNoMCMatchOverlap(0),      fhMCPi0NoneNLMPairNoMCMatchOverlap(0),
  fhMCPi0DecayPhotonHitHighLM(0),            fhMCPi0DecayPhotonAdjHighLM(0),
  fhMCPi0DecayPhotonHitOtherLM(0),           fhMCPi0DecayPhotonAdjOtherLM(0),
  fhMCPi0DecayPhotonAdjacent(0),             fhMCPi0DecayPhotonHitNoLM(0),
  fhMCPi0DecayPhotonHitHighLMOverlap(0),     fhMCPi0DecayPhotonAdjHighLMOverlap(0),
  fhMCPi0DecayPhotonHitOtherLMOverlap(0),    fhMCPi0DecayPhotonAdjOtherLMOverlap(0),
  fhMCPi0DecayPhotonAdjacentOverlap(0),      fhMCPi0DecayPhotonHitNoLMOverlap(0),
  fhMCEOverlapType(0),                       fhMCEOverlapTypeMatch(0)
{
  //default ctor
  
  // Init array of histograms
  for(Int_t i = 0; i < 7; i++)
  {
    for(Int_t j = 0; j < 2; j++)
    {
      fhMassNLocMax1[i][j]  = 0;
      fhMassNLocMax2[i][j]  = 0;
      fhMassNLocMaxN[i][j]  = 0;
      fhMassSplitENLocMax1[i][j]  = 0;
      fhMassSplitENLocMax2[i][j]  = 0;
      fhMassSplitENLocMaxN[i][j]  = 0;
      fhNLocMax[i][j]       = 0;
      fhNLocMaxM02Cut[i][j] = 0;
      fhSplitClusterENLocMax   [i][j] = 0;
      fhSplitClusterEPi0NLocMax[i][j] = 0;
      fhM02NLocMax1[i][j]   = 0;
      fhM02NLocMax2[i][j]   = 0;
      fhM02NLocMaxN[i][j]   = 0;
      fhNCellNLocMax1[i][j] = 0;
      fhNCellNLocMax2[i][j] = 0;
      fhNCellNLocMaxN[i][j] = 0;
      fhM02Pi0NLocMax1[i][j] = 0;
      fhM02EtaNLocMax1[i][j] = 0;
      fhM02ConNLocMax1[i][j] = 0;
      fhM02Pi0NLocMax2[i][j] = 0;
      fhM02EtaNLocMax2[i][j] = 0;
      fhM02ConNLocMax2[i][j] = 0;
      fhM02Pi0NLocMaxN[i][j] = 0;
      fhM02EtaNLocMaxN[i][j] = 0;
      fhM02ConNLocMaxN[i][j] = 0;
      
      fhMassPi0NLocMax1[i][j] = 0;
      fhMassEtaNLocMax1[i][j] = 0;
      fhMassConNLocMax1[i][j] = 0;
      fhMassPi0NLocMax2[i][j] = 0;
      fhMassEtaNLocMax2[i][j] = 0;
      fhMassConNLocMax2[i][j] = 0;
      fhMassPi0NLocMaxN[i][j] = 0;
      fhMassEtaNLocMaxN[i][j] = 0;
      fhMassConNLocMaxN[i][j] = 0;
      
      fhNCellPi0NLocMax1[i][j] = 0;
      fhNCellEtaNLocMax1[i][j] = 0;
      fhNCellPi0NLocMax2[i][j] = 0;
      fhNCellEtaNLocMax2[i][j] = 0;
      fhNCellPi0NLocMaxN[i][j] = 0;
      fhNCellEtaNLocMaxN[i][j] = 0;
      
      fhAsyPi0NLocMax1[i][j] = 0;
      fhAsyEtaNLocMax1[i][j] = 0;
      fhAsyConNLocMax1[i][j] = 0;
      fhAsyPi0NLocMax2[i][j] = 0;
      fhAsyEtaNLocMax2[i][j] = 0;
      fhAsyConNLocMax2[i][j] = 0;
      fhAsyPi0NLocMaxN[i][j] = 0;
      fhAsyEtaNLocMaxN[i][j] = 0;
      fhAsyConNLocMaxN[i][j] = 0;      
      
      fhMassM02NLocMax1[i][j]= 0;
      fhMassM02NLocMax2[i][j]= 0;
      fhMassM02NLocMaxN[i][j]= 0;
      
      fhMassSplitEPi0NLocMax1[i][j]  = 0;
      fhMassSplitEPi0NLocMax2[i][j]  = 0;
      fhMassSplitEPi0NLocMaxN[i][j]  = 0;

      fhMassSplitEAfterCutsNLocMax1[i][j]  = 0;
      fhMassSplitEAfterCutsNLocMax2[i][j]  = 0;
      fhMassSplitEAfterCutsNLocMaxN[i][j]  = 0;

      
      fhMassDispEtaNLocMax1[i][j]= 0;
      fhMassDispEtaNLocMax2[i][j]= 0;
      fhMassDispEtaNLocMaxN[i][j]= 0;      
      fhMassDispPhiNLocMax1[i][j]= 0;
      fhMassDispPhiNLocMax2[i][j]= 0;
      fhMassDispPhiNLocMaxN[i][j]= 0;      
      fhMassDispAsyNLocMax1[i][j]= 0;
      fhMassDispAsyNLocMax2[i][j]= 0;
      fhMassDispAsyNLocMaxN[i][j]= 0;      
      
      fhSplitEFractionNLocMax1[i][j]=0;
      fhSplitEFractionNLocMax2[i][j]=0;
      fhSplitEFractionNLocMaxN[i][j]=0;
      
      fhAnglePairNLocMax1         [i][j] = 0;
      fhAnglePairNLocMax2         [i][j] = 0;
      fhAnglePairNLocMaxN         [i][j] = 0;

      fhAnglePairAfterCutsNLocMax1[i][j] = 0;
      fhAnglePairAfterCutsNLocMax2[i][j] = 0;
      fhAnglePairAfterCutsNLocMaxN[i][j] = 0;

      fhAnglePairPi0NLocMax1      [i][j] = 0;
      fhAnglePairPi0NLocMax2      [i][j] = 0;
      fhAnglePairPi0NLocMaxN      [i][j] = 0;
      
      fhAnglePairMassNLocMax1     [i][j] = 0;
      fhAnglePairMassNLocMax2     [i][j] = 0;
      fhAnglePairMassNLocMaxN     [i][j] = 0;

      fhAnglePairM02NLocMax1      [i][j] = 0;
      fhAnglePairM02NLocMax2      [i][j] = 0;
      fhAnglePairM02NLocMaxN      [i][j] = 0;
      
      fhAnglePairOverM02NLocMax1  [i][j] = 0;
      fhAnglePairOverM02NLocMax2  [i][j] = 0;
      fhAnglePairOverM02NLocMaxN  [i][j] = 0;
      
      fhAnglePairOverM02NLocMax1Overlap0[i][j] = 0;
      fhAnglePairOverM02NLocMax2Overlap0[i][j] = 0;
      fhAnglePairOverM02NLocMaxNOverlap0[i][j] = 0;
      
      fhCosThStarNLocMax1         [i][j] = 0;
      fhCosThStarNLocMax2         [i][j] = 0;
      fhCosThStarNLocMaxN         [i][j] = 0;
      
      fhCosThStarAfterCutsNLocMax1[i][j] = 0;
      fhCosThStarAfterCutsNLocMax2[i][j] = 0;
      fhCosThStarAfterCutsNLocMaxN[i][j] = 0;
      
      fhCosThStarPi0NLocMax1      [i][j] = 0;
      fhCosThStarPi0NLocMax2      [i][j] = 0;
      fhCosThStarPi0NLocMaxN      [i][j] = 0;
      
      fhMCGenFracNLocMax1[i][j]= 0;
      fhMCGenFracNLocMax2[i][j]= 0;
      fhMCGenFracNLocMaxN[i][j]= 0;

      fhMCGenFracNLocMax1NoOverlap[i][j]= 0;
      fhMCGenFracNLocMax2NoOverlap[i][j]= 0;
      fhMCGenFracNLocMaxNNoOverlap[i][j]= 0;
      
      fhMCGenSplitEFracNLocMax1[i][j]= 0;
      fhMCGenSplitEFracNLocMax2[i][j]= 0;
      fhMCGenSplitEFracNLocMaxN[i][j]= 0;    

      fhMCGenSplitEFracNLocMax1NoOverlap[i][j]= 0;
      fhMCGenSplitEFracNLocMax2NoOverlap[i][j]= 0;
      fhMCGenSplitEFracNLocMaxNNoOverlap[i][j]= 0;
      
      fhMCGenEFracvsSplitEFracNLocMax1[i][j]= 0;
      fhMCGenEFracvsSplitEFracNLocMax2[i][j]= 0;
      fhMCGenEFracvsSplitEFracNLocMaxN[i][j]= 0;    
      
      fhMCGenEvsSplitENLocMax1[i][j]= 0;
      fhMCGenEvsSplitENLocMax2[i][j]= 0;
      fhMCGenEvsSplitENLocMaxN[i][j]= 0;     
      
      fhAsymNLocMax1 [i][j] = 0;
      fhAsymNLocMax2 [i][j] = 0;
      fhAsymNLocMaxN [i][j] = 0;
      
      fhMassAfterCutsNLocMax1[i][j] = 0;
      fhMassAfterCutsNLocMax2[i][j] = 0;
      fhMassAfterCutsNLocMaxN[i][j] = 0;

      
      fhSplitEFractionAfterCutsNLocMax1[i][j] = 0 ;
      fhSplitEFractionAfterCutsNLocMax2[i][j] = 0 ;
      fhSplitEFractionAfterCutsNLocMaxN[i][j] = 0 ;
    }
   
    for(Int_t jj = 0; jj < 4; jj++)
    {
      fhM02MCGenFracNLocMax1Ebin[i][jj] = 0;
      fhM02MCGenFracNLocMax2Ebin[i][jj] = 0;
      fhM02MCGenFracNLocMaxNEbin[i][jj] = 0;
      
      fhMassMCGenFracNLocMax1Ebin[i][jj]= 0;
      fhMassMCGenFracNLocMax2Ebin[i][jj]= 0;
      fhMassMCGenFracNLocMaxNEbin[i][jj]= 0;
      
      fhMCGenFracNLocMaxEbin[i][jj]       = 0;
      fhMCGenFracNLocMaxEbinMatched[i][jj]= 0;
      
      fhMassSplitEFractionNLocMax1Ebin[i][jj] = 0;
      fhMassSplitEFractionNLocMax2Ebin[i][jj] = 0;
      fhMassSplitEFractionNLocMaxNEbin[i][jj] = 0;
    }
    
    fhTrackMatchedDEtaNLocMax1[i] = 0;
    fhTrackMatchedDPhiNLocMax1[i] = 0;
    fhTrackMatchedDEtaNLocMax2[i] = 0;
    fhTrackMatchedDPhiNLocMax2[i] = 0; 
    fhTrackMatchedDEtaNLocMaxN[i] = 0; 
    fhTrackMatchedDPhiNLocMaxN[i] = 0; 

    fhTrackMatchedDEtaNLocMax1Pos[i] = 0;
    fhTrackMatchedDPhiNLocMax1Pos[i] = 0;
    fhTrackMatchedDEtaNLocMax2Pos[i] = 0;
    fhTrackMatchedDPhiNLocMax2Pos[i] = 0;
    fhTrackMatchedDEtaNLocMaxNPos[i] = 0;
    fhTrackMatchedDPhiNLocMaxNPos[i] = 0;

    fhTrackMatchedDEtaNLocMax1Neg[i] = 0;
    fhTrackMatchedDPhiNLocMax1Neg[i] = 0;
    fhTrackMatchedDEtaNLocMax2Neg[i] = 0;
    fhTrackMatchedDPhiNLocMax2Neg[i] = 0;
    fhTrackMatchedDEtaNLocMaxNNeg[i] = 0;
    fhTrackMatchedDPhiNLocMaxNNeg[i] = 0;
    
    for(Int_t nlm = 0; nlm < 3; nlm++)
    {
      fhMCEM02Overlap0     [nlm][i] = 0;
      fhMCEM02Overlap1     [nlm][i] = 0;
      fhMCEM02OverlapN     [nlm][i] = 0;
      fhMCEM02Overlap0Match[nlm][i] = 0;
      fhMCEM02Overlap1Match[nlm][i] = 0;
      fhMCEM02OverlapNMatch[nlm][i] = 0;
      
      fhMCEMassOverlap0     [nlm][i] = 0;
      fhMCEMassOverlap1     [nlm][i] = 0;
      fhMCEMassOverlapN     [nlm][i] = 0;
      fhMCEMassOverlap0Match[nlm][i] = 0;
      fhMCEMassOverlap1Match[nlm][i] = 0;
      fhMCEMassOverlapNMatch[nlm][i] = 0;

      fhMCEAsymOverlap0     [nlm][i] = 0;
      fhMCEAsymOverlap1     [nlm][i] = 0;
      fhMCEAsymOverlapN     [nlm][i] = 0;
      fhMCEAsymOverlap0Match[nlm][i] = 0;
      fhMCEAsymOverlap1Match[nlm][i] = 0;
      fhMCEAsymOverlapNMatch[nlm][i] = 0;

      fhMCENCellOverlap0     [nlm][i] = 0;
      fhMCENCellOverlap1     [nlm][i] = 0;
      fhMCENCellOverlapN     [nlm][i] = 0;
      fhMCENCellOverlap0Match[nlm][i] = 0;
      fhMCENCellOverlap1Match[nlm][i] = 0;
      fhMCENCellOverlapNMatch[nlm][i] = 0;
      
      fhMCEEpriOverlap0     [nlm][i] = 0;
      fhMCEEpriOverlap1     [nlm][i] = 0;
      fhMCEEpriOverlapN     [nlm][i] = 0;
      fhMCEEpriOverlap0Match[nlm][i] = 0;
      fhMCEEpriOverlap1Match[nlm][i] = 0;
      fhMCEEpriOverlapNMatch[nlm][i] = 0;

      fhMCEEpriOverlap0IdPi0[nlm][i] = 0;
      fhMCEEpriOverlap1IdPi0[nlm][i] = 0;
      fhMCEEpriOverlapNIdPi0[nlm][i] = 0;
      
      fhMCESplitEFracOverlap0     [nlm][i] = 0;
      fhMCESplitEFracOverlap1     [nlm][i] = 0;
      fhMCESplitEFracOverlapN     [nlm][i] = 0;
      fhMCESplitEFracOverlap0Match[nlm][i] = 0;
      fhMCESplitEFracOverlap1Match[nlm][i] = 0;
      fhMCESplitEFracOverlapNMatch[nlm][i] = 0;
      
      fhMCENOverlaps       [nlm][i] = 0;
      fhMCENOverlapsMatch  [nlm][i] = 0;
            
      if(i > 3) continue ;
      
      fhMCPi0MassM02Overlap0     [nlm][i] = 0;
      fhMCPi0MassM02Overlap1     [nlm][i] = 0;
      fhMCPi0MassM02OverlapN     [nlm][i] = 0;
      fhMCPi0MassM02Overlap0Match[nlm][i] = 0;
      fhMCPi0MassM02Overlap1Match[nlm][i] = 0;
      fhMCPi0MassM02OverlapNMatch[nlm][i] = 0;
    }
  }
   
  for(Int_t i = 0; i < 2; i++)
  {
    fhSplitEFractionvsAsyNLocMax1[i] = 0;
    fhSplitEFractionvsAsyNLocMax2[i] = 0; 
    fhSplitEFractionvsAsyNLocMaxN[i] = 0;    
  }
  
  for(Int_t i = 0; i < 4; i++)
  {
    fhMassM02NLocMax1Ebin[i] = 0 ;
    fhMassM02NLocMax2Ebin[i] = 0 ;
    fhMassM02NLocMaxNEbin[i] = 0 ;

    fhMassAsyNLocMax1Ebin[i] = 0 ;
    fhMassAsyNLocMax2Ebin[i] = 0 ;
    fhMassAsyNLocMaxNEbin[i] = 0 ;

    fhAsyMCGenRecoNLocMax1EbinPi0[i] = 0 ;
    fhAsyMCGenRecoNLocMax2EbinPi0[i] = 0 ;
    fhAsyMCGenRecoNLocMaxNEbinPi0[i] = 0 ;
    
    fhMassDispEtaNLocMax1Ebin[i] = 0 ;
    fhMassDispEtaNLocMax2Ebin[i] = 0 ;
    fhMassDispEtaNLocMaxNEbin[i] = 0 ;
    
    fhMassDispPhiNLocMax1Ebin[i] = 0 ;
    fhMassDispPhiNLocMax2Ebin[i] = 0 ;
    fhMassDispPhiNLocMaxNEbin[i] = 0 ;    
    
    fhMassDispAsyNLocMax1Ebin[i] = 0 ;
    fhMassDispAsyNLocMax2Ebin[i] = 0 ;
    fhMassDispAsyNLocMaxNEbin[i] = 0 ;    

    fhMCAsymM02NLocMax1MCPi0Ebin[i] = 0 ;
    fhMCAsymM02NLocMax2MCPi0Ebin[i] = 0 ;
    fhMCAsymM02NLocMaxNMCPi0Ebin[i] = 0 ;
  }
  
  for(Int_t nlm = 0; nlm < 3; nlm++)
  {
    
    fhMCPi0DecayPhotonHitHighLMDiffELM1 [nlm] = 0 ;
    fhMCPi0DecayPhotonAdjHighLMDiffELM1 [nlm] = 0 ;           
    fhMCPi0DecayPhotonHitOtherLMDiffELM1[nlm] = 0 ;           
    fhMCPi0DecayPhotonAdjOtherLMDiffELM1[nlm] = 0 ;            
    
    fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1 [nlm] = 0 ;     
    fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1 [nlm] = 0 ;      
    fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM1[nlm] = 0 ;     
    fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM1[nlm] = 0 ;     
    
    fhMCPi0DecayPhotonHitHighLMDiffELM2 [nlm] = 0 ;           
    fhMCPi0DecayPhotonAdjHighLMDiffELM2 [nlm] = 0 ;            
    fhMCPi0DecayPhotonHitOtherLMDiffELM2[nlm] = 0 ;            
    fhMCPi0DecayPhotonAdjOtherLMDiffELM2[nlm] = 0 ;         
    
    fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2 [nlm] = 0 ;    
    fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2 [nlm] = 0 ;      
    fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM2[nlm] = 0 ;     
    fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM2[nlm] = 0 ;
    
    fhMCPi0DecayPhotonHitHighLMDiffELM1vsELM1 [nlm] = 0 ;
    fhMCPi0DecayPhotonAdjHighLMDiffELM1vsELM1 [nlm] = 0 ;
    fhMCPi0DecayPhotonHitOtherLMDiffELM1vsELM1[nlm] = 0 ;
    fhMCPi0DecayPhotonAdjOtherLMDiffELM1vsELM1[nlm] = 0 ;
    
    fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1vsELM1 [nlm] = 0 ;
    fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1vsELM1 [nlm] = 0 ;
    fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM1vsELM1[nlm] = 0 ;
    fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM1vsELM1[nlm] = 0 ;
    
    fhMCPi0DecayPhotonHitHighLMDiffELM2vsELM2 [nlm] = 0 ;
    fhMCPi0DecayPhotonAdjHighLMDiffELM2vsELM2 [nlm] = 0 ;
    fhMCPi0DecayPhotonHitOtherLMDiffELM2vsELM2[nlm] = 0 ;
    fhMCPi0DecayPhotonAdjOtherLMDiffELM2vsELM2[nlm] = 0 ;
    
    fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2vsELM2 [nlm] = 0 ;
    fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2vsELM2 [nlm] = 0 ;
    fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM2vsELM2[nlm] = 0 ;
    fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM2vsELM2[nlm] = 0 ;
    
    fhMCPi0DecayPhotonHitHighLMMass [nlm] = 0 ;
    fhMCPi0DecayPhotonAdjHighLMMass [nlm] = 0 ;                 
    fhMCPi0DecayPhotonHitOtherLMMass[nlm] = 0 ;              
    fhMCPi0DecayPhotonAdjOtherLMMass[nlm] = 0 ;               
    fhMCPi0DecayPhotonAdjacentMass  [nlm] = 0 ;                  
    fhMCPi0DecayPhotonHitNoLMMass   [nlm] = 0 ;                  
    
    fhMCPi0DecayPhotonHitHighLMOverlapMass [nlm] = 0 ;
    fhMCPi0DecayPhotonAdjHighLMOverlapMass [nlm] = 0 ;          
    fhMCPi0DecayPhotonHitOtherLMOverlapMass[nlm] = 0 ;        
    fhMCPi0DecayPhotonAdjOtherLMOverlapMass[nlm] = 0 ;        
    fhMCPi0DecayPhotonAdjacentOverlapMass  [nlm] = 0 ;
    fhMCPi0DecayPhotonHitNoLMOverlapMass   [nlm] = 0 ;          
    
    fhPi0CellE       [nlm] = 0 ;
    fhPi0CellEFrac   [nlm] = 0 ;
    fhPi0CellLogEFrac[nlm] = 0 ;
    
    fhPi0CellEMaxEMax2Frac   [nlm] = 0 ;
    fhPi0CellEMaxClusterFrac [nlm] = 0 ;
    fhPi0CellEMax2ClusterFrac[nlm] = 0 ;

    fhPi0CellEMaxFrac [nlm] = 0 ;
    fhPi0CellEMax2Frac[nlm] = 0 ;
    
    for(Int_t i = 0; i < 10; i++)
    {
      fhM02WeightPi0  [nlm][i] = 0;
      fhM02ECellCutPi0[nlm][i] = 0;
    }
    
    fhMassBadDistClose[nlm] = 0;
    fhM02BadDistClose [nlm] = 0;
    fhMassOnBorder    [nlm] = 0;
    fhM02OnBorder     [nlm] = 0;
    
    fhAsyMCGenRecoDiffMCPi0    [nlm] = 0;
    fhAsyMCGenRecoDiffMCPi0Conv[nlm] = 0;

  }
  
  for(Int_t i = 0; i < 7; i++)
  {
    for(Int_t j = 0; j < 4; j++)
    {
      
      fhArmNLocMax1[i][j]  = 0;
      fhArmNLocMax2[i][j]  = 0;
      fhArmNLocMaxN[i][j]  = 0;
      
      fhArmPi0NLocMax1[i][j] = 0;
      fhArmPi0NLocMax2[i][j] = 0;
      fhArmPi0NLocMaxN[i][j] = 0;
      
      fhArmAfterCutsNLocMax1[i][j] = 0;
      fhArmAfterCutsNLocMax2[i][j] = 0;
      fhArmAfterCutsNLocMaxN[i][j] = 0;
      
    }
  }
  
  InitParameters();

}

//___________________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::CheckLocalMaximaMCOrigin(AliVCluster* cluster, Int_t mcindex, Int_t noverlaps,
                                                                Float_t e1,    Float_t e2,    Float_t mass)
                                                                //Float_t m02,
                                                                //TLorentzVector l1, TLorentzVector l2)
{
  // Check origin NLM tower of the cluster, when MC gives merged pi0
  
  if(mcindex != kmcPi0 && mcindex != kmcPi0Conv) return;

  const UInt_t nc = cluster->GetNCells();
  Int_t   list[nc];
  Float_t elist[nc];
  Int_t nMax = GetCaloUtils()->GetNumberOfLocalMaxima(cluster, GetEMCALCells(),list, elist);
  
  
  //// PRINTS /////
  
  //if(mcindex==kmcPi0)     printf("** Normal Pi0 **\n");
  //if(mcindex==kmcPi0Conv) printf("** Converted Pi0 **\n");

//  if((mass < 0.06 || mass > 1.8) && mcindex==kmcPi0 && noverlaps == 0)
//  {
//     printf("** N max %d - Overlaps = %d **, mass %2.2f, m02 %2.2f, Cl1(E,eta,phi)=(%2.2f,%2.2f,%2.2f),Cl2(E,eta,phi)=(%2.2f,%2.2f,%2.2f), mass(1,2) %2.2f \n",
//            nMax, noverlaps,mass,m02,
//            l1.E(),l1.Eta(),l1.Phi()*TMath::RadToDeg(),
//            l2.E(),l2.Eta(),l2.Phi()*TMath::RadToDeg(), (l1+l2).M());
//    
//    // Study the mothers of cluster
//    printf("Cluster MC labels %d \n", cluster->GetNLabels());
//    for (UInt_t ilab = 0; ilab < cluster->GetNLabels(); ilab++ )
//    {
//      Int_t mclabel = cluster->GetLabels()[ilab];
//      
//      Bool_t  mOK = 0;
//      Int_t   mpdg = -999999;
//      Int_t   mstatus = -1;
//      Int_t   grandLabel = -1;
//      TLorentzVector mother = GetMCAnalysisUtils()->GetMother(mclabel,GetReader(),mpdg,mstatus,mOK,grandLabel);
//      
//      printf("******** mother %d : Label %d, pdg %d; status %d, E %2.2f, Eta %2.2f, Phi %2.2f, ok %d, mother label %d\n",
//             ilab, mclabel, mpdg, mstatus,mother.E(), mother.Eta(),mother.Phi()*TMath::RadToDeg(),mOK,grandLabel);
//      
//      if( ( mpdg == 22 || TMath::Abs(mpdg)==11 ) && grandLabel >=0 )
//      {
//        while( ( mpdg == 22 || TMath::Abs(mpdg)==11 ) && grandLabel >=0 )
//        {
//          Int_t newLabel = -1;
//          TLorentzVector grandmother = GetMCAnalysisUtils()->GetMother(grandLabel,GetReader(),mpdg,mstatus,mOK,newLabel);
//          printf("\t grandmother %d : Label %d, pdg %d; status %d, E %2.2f, Eta %2.2f, Phi %2.2f, ok %d, mother label %d\n",
//                 ilab, grandLabel, mpdg, mstatus,grandmother.E(), grandmother.Eta(), grandmother.Phi()*TMath::RadToDeg(),mOK,newLabel);
//          grandLabel = newLabel;
//          
//        }
//      }
//    }
//    
//    printf("Cells in cluster %d\n",cluster->GetNCells() );
//    for(Int_t icell = 0; icell < cluster->GetNCells(); icell++)
//    {
//      Int_t absIdCell = cluster->GetCellAbsId(icell);
//      Int_t mcLabel   = GetEMCALCells()->GetCellMCLabel(absIdCell);
//      GetReader()->RemapMCLabelForAODs(mcLabel);
//      Int_t ietac=-1; Int_t iphic = 0; Int_t rcuc = 0;
//      Int_t smc = GetModuleNumberCellIndexes(absIdCell,fCalorimeter, ietac, iphic, rcuc);
//
//      printf(" \t cell i %d, abs %d, amp %2.3f, mclabel %d, (sm,ieta,iphi)=(%d,%d,%d)\n",icell,absIdCell,GetEMCALCells()->GetCellAmplitude(absIdCell),mcLabel,smc,ietac,iphic);
//    }
//  }
  //// PRINTS /////
  
  
  //If only one maxima, consider all the towers in the cluster
  if(nMax==1)
  {
    for (UInt_t icell = 0; icell < nc; icell++ )
    {
      list [icell] = cluster->GetCellAbsId(icell);
      elist[icell] = GetEMCALCells()->GetCellAmplitude(list[icell]);
    }
  }
  
  Int_t nmaxima = nMax;
  if(nMax==1) nmaxima = nc ;
  
  //Find highest energy Local Maxima Towers
  Int_t   imax  = -1;
  Int_t   imax2 = -1;
  Float_t emax  = -1;
  Float_t emax2 = -1;
  for(Int_t i = 0; i < nmaxima; i++)
  {
    //printf("i %d: AbsId %d; E %2.3f\n",i,list[i],elist[i]);
    if(elist[i] > emax)
    {
      imax = i;
      emax = elist[i];
    }
  }
  //Find second highest
  for(Int_t i = 0; i < nmaxima; i++)
  {
    if(i==imax) continue;
    
    //printf("j %d: AbsId %d; E %2.3f\n",i,list[i],elist[i]);
    
    
    if(elist[i] > emax2)
    {
      imax2 = i;
      emax2 = elist[i];
    }
  }
  
// if((mass < 0.06 || mass > 1.8) && mcindex==kmcPi0 && noverlaps == 0)
//    printf("Local maxima: a) index %d, absId %d; b) index %d, absId %d\n",imax, list[imax], imax2, list[imax2]);
  
  //---------------------------------------------------------
  //---------------------------------------------------------
  // Compare ancestors of all local maxima at cell MC level
  //---------------------------------------------------------
  //---------------------------------------------------------
  
  // Check that the highest mc label and the max cluster label are the same
  Int_t mcLabelMax = -1 ;
  if(imax >=0 )
  {
    mcLabelMax = GetEMCALCells()->GetCellMCLabel(list[imax]);
    GetReader()->RemapMCLabelForAODs(mcLabelMax);
  }
  
  Int_t mcLabelMax2 = -1 ;
  if(imax >=0 )
  {
    mcLabelMax2 = GetEMCALCells()->GetCellMCLabel(list[imax2]);
    GetReader()->RemapMCLabelForAODs(mcLabelMax2);
  }
  
  Int_t mcLabelclusterMax = cluster->GetLabels()[0];
  Bool_t matchHighLMAndHighMC = kFALSE;
  
  //printf("MC label: LM1 %d, LM2 %d, cluster %d\n",mcLabelMax,mcLabelMax2,mcLabelclusterMax);
  
  if(mcLabelclusterMax == mcLabelMax && mcLabelclusterMax >= 0)
  {
    matchHighLMAndHighMC = kTRUE;
    //printf("\t *** MATCH cluster and LM maximum ***\n");
  }
  else
  {
     //printf("\t *** NO MATCH cluster and LM maximum, check second ***\n");
    if(mcLabelclusterMax == mcLabelMax2 && mcLabelclusterMax >= 0)
    {
      //printf("\t \t *** MATCH cluster and 2nd LM maximum ***\n");
      matchHighLMAndHighMC = kTRUE;
    }
    else
    {
      //printf("\t \t *** NO MATCH***\n");
      matchHighLMAndHighMC = kFALSE;
    }
  }
  
  // Compare the common ancestors of the 2 highest energy local maxima
  Int_t ancPDG = 0, ancStatus = -1;
  TLorentzVector momentum; TVector3 prodVertex;
  Int_t ancLabel = 0;
  Bool_t high = kFALSE;
  Bool_t low  = kFALSE;

//  // print maxima origin
//  for(Int_t i = 0; i < nMax; i++)
//  {
//    Int_t mcLabel1 = GetEMCALCells()->GetCellMCLabel(list[i]);
//    GetReader()->RemapMCLabelForAODs(mcLabel1);
//    
//    Bool_t ok  =kFALSE,gok = kFALSE;
//    Int_t pdg    = -22222, status   = -1;
//    Int_t gpdg   = -22222, gstatus  = -1;
//    Int_t ggpdg  = -22222, ggstatus = -1;
//    Int_t gLabel = -1, ggLabel = -1;
//    TLorentzVector primary   =GetMCAnalysisUtils()->GetMother     (mcLabel1,GetReader(),  pdg,  status, ok);
//    TLorentzVector gprimary  =GetMCAnalysisUtils()->GetGrandMother(mcLabel1,GetReader(), gpdg, gstatus,gok, gLabel,ggLabel);
//    TLorentzVector ggprimary =GetMCAnalysisUtils()->GetMother(ggLabel  ,GetReader(),ggpdg,ggstatus,gok);
//    printf("Max index %d; mother: Label %d; PDG %d; E %2.2f - grand mother label %d; PDG %d; E %2.2f- great grand mother label %d; PDG %d; E %2.2f\n",
//           i,mcLabel1,pdg,primary.E(), gLabel,gpdg,gprimary.E(), ggLabel,ggpdg,ggprimary.E());
//  }

  for(Int_t i = 0; i < nmaxima-1; i++)
  {
    Int_t mcLabel1 = GetEMCALCells()->GetCellMCLabel(list[i]);
    GetReader()->RemapMCLabelForAODs(mcLabel1);
 
    for(Int_t j = i+1; j < nmaxima; j++)
    {
      Int_t mcLabel2 = GetEMCALCells()->GetCellMCLabel(list[j]);
      GetReader()->RemapMCLabelForAODs(mcLabel2);
      
      if(mcLabel1 < 0 || mcLabel2 < 0 )
      {
        //printf("\t i %d label %d - j %d label %d; skip!\n",i,mcLabel1,j,mcLabel2);
        continue;
      }
      
      ancLabel = GetMCAnalysisUtils()->CheckCommonAncestor(mcLabel1,mcLabel2,
                                                           GetReader(),ancPDG,ancStatus,momentum,prodVertex);
      if(ancPDG==111)
      {
        if((i==imax && j==imax2) ||  (j==imax && i==imax2))
          high = kTRUE;
        else
          low = kTRUE;
      }
      else if(ancPDG==22 || TMath::Abs(ancPDG)==11)
      {
        // If both bits are set, it could be that one of the maxima had a conversion
        // reset the bit in this case
        if(high && low)
        {
          //printf("\t Reset low bit\n");
          low = kFALSE;
        }
      }
     
      Bool_t ok  =kFALSE;
      Int_t pdg = -22222, status = -1;
      TLorentzVector primary  =GetMCAnalysisUtils()->GetMother(ancLabel,GetReader(), pdg, status, ok);
      //printf("\t i %d label %d - j %d label %d; ancestor label %d, PDG %d-%d; E %2.2f; high %d, any %d \n",i,mcLabel1,j,mcLabel2, ancLabel, ancPDG,pdg, primary.E(), high, low);

    }
  }
  
  Float_t en = cluster->E();
  
// if((mass < 0.06 || mass > 1.8) && mcindex==kmcPi0 && noverlaps == 0)
//    printf("Cell MC match: nMax %d; Match MC? %d; high %d; low %d\n",nMax,matchHighLMAndHighMC,high,low);
  
  if(!noverlaps)
  {
    if(matchHighLMAndHighMC)
    {
      if     (high && !low)  fhMCPi0HighNLMPair->Fill(en,nMax);
      else if(low  && !high) fhMCPi0LowNLMPair ->Fill(en,nMax);
      else if(low  &&  high) fhMCPi0AnyNLMPair ->Fill(en,nMax);
      else                   fhMCPi0NoneNLMPair->Fill(en,nMax);
    }
    else
    {
      if     (high && !low)  fhMCPi0HighNLMPairNoMCMatch->Fill(en,nMax);
      else if(low  && !high) fhMCPi0LowNLMPairNoMCMatch ->Fill(en,nMax);
      else if(low  &&  high) fhMCPi0AnyNLMPairNoMCMatch ->Fill(en,nMax);
      else                   fhMCPi0NoneNLMPairNoMCMatch->Fill(en,nMax);
    }
  }
  else
  {
    if(matchHighLMAndHighMC)
    {
      if     (high && !low)  fhMCPi0HighNLMPairOverlap->Fill(en,nMax);
      else if(low  && !high) fhMCPi0LowNLMPairOverlap->Fill(en,nMax);
      else if(low  &&  high) fhMCPi0AnyNLMPairOverlap->Fill(en,nMax);
      else                   fhMCPi0NoneNLMPairOverlap->Fill(en,nMax);
    }
    else
    {
      if     (high && !low)  fhMCPi0HighNLMPairNoMCMatchOverlap->Fill(en,nMax);
      else if(low  && !high) fhMCPi0LowNLMPairNoMCMatchOverlap->Fill(en,nMax);
      else if(low  &&  high) fhMCPi0AnyNLMPairNoMCMatchOverlap->Fill(en,nMax);
      else                   fhMCPi0NoneNLMPairNoMCMatchOverlap->Fill(en,nMax);
    }  
  }
  
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // Compare MC decay photon projection to cell location and Local Maxima
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  
  // Get the mother pi0
  
  Bool_t ok     = kFALSE;
  Int_t pdg    = -22222, status   = -1;
  Int_t gLabel = -1;
  
  Int_t label = cluster->GetLabel();
  TLorentzVector pi0Kine;
    
  while( pdg!=111 && label>=0 )
  {
    pi0Kine = GetMCAnalysisUtils()->GetGrandMother(label,GetReader(),pdg,status,ok, label,gLabel);
  }
  
  if(pdg!=111 || label < 0)
  {
    Info("CheckLocalMaximaMCOrigin","Mother Pi0 not found!\n");
    return;
  }
  
  Int_t nDaugthers = GetMCAnalysisUtils()->GetNDaughters(label,GetReader(),ok);
  
  if(nDaugthers != 2)
  {
    Info("CheckLocalMaximaMCOrigin","N daughters %d !=2!\n",nDaugthers);
    return;
  }
  
  // Get daughter photon kinematics
  Int_t pdg0 = -22222, status0   = -1; Int_t label0 = -1;
  TLorentzVector photon0Kine = GetMCAnalysisUtils()->GetDaughter(0,label,GetReader(),pdg0,status0,ok,label0);
  Int_t pdg1 = -22222, status1   = -1; Int_t label1 = -1;
  TLorentzVector photon1Kine = GetMCAnalysisUtils()->GetDaughter(1,label,GetReader(),pdg1,status1,ok,label1);

  if(pdg1!=22 || pdg0 != 22)
  {
    Info("CheckLocalMaximaMCOrigin","Wrong daughters PDG: photon0 %d - photon1 %d\n",pdg0,pdg1);
    return;
  }
  
  // In what cells did the photons hit
  Float_t eta0 = photon0Kine.Eta();
  Float_t eta1 = photon1Kine.Eta();
  
  Float_t phi0 = photon0Kine.Phi();
  Float_t phi1 = photon1Kine.Phi();

// if((mass < 0.06 || mass > 1.8) && mcindex==kmcPi0 && noverlaps == 0)
//  {
//    printf("MC pi0 label %d E  %2.2f, eta %2.2f, phi %2.2f, mass (ph1, ph2) %2.2f: \n \t photon0 label %d E %2.2f, eta %2.2f, phi %2.2f \n \t photon1 label %d E %2.2f eta %2.2f, phi %2.2f\n",
//           label , pi0Kine.E()    , pi0Kine.Eta(),pi0Kine.Phi()*TMath::RadToDeg(), (photon0Kine+photon1Kine).M(),
//           label0, photon0Kine.E(),          eta0,         phi0*TMath::RadToDeg(),
//           label1, photon1Kine.E(),          eta1,         phi1*TMath::RadToDeg());
//    
//    TLorentzVector momclus;
//    cluster->GetMomentum(momclus,GetVertex(0));
//    printf("Cluster E %2.2F eta %2.2f, phi %2.2f, dist to bad %2.2f\n",momclus.E(),momclus.Eta(),momclus.Phi()*TMath::RadToDeg(), cluster->GetDistanceToBadChannel());
//  }
  
  if(phi0 < 0 ) phi0+=TMath::TwoPi();
  if(phi1 < 0 ) phi1+=TMath::TwoPi();
  
  Int_t absId0=-1, absId1=-1;
  GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(eta0, phi0, absId0);
  GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(eta1, phi1, absId1);
  
  if(absId1 < 0 || absId1 < 0)
  {
    //printf("AliAnaInsideClusterInvariantMass::CheckLocalMaximaMCOrigin(() -  Photon hit AbsId: photon0 %d - photon1 %d\n",absId0,absId1);
    return;
  }
  
  //-----------------------------------------------
  // Check that the 2 photons hit the Local Maxima
  //-----------------------------------------------
  
  
// if((mass < 0.06 || mass > 1.8) && mcindex==kmcPi0 && noverlaps == 0)
//  {
//    printf("Photons AbsId (%d,%d); Local Maxima AbsId(%d,%d)\n",absId0,absId1,list[imax],list[imax2]);
//    printf("Photon1 (eta,phi)=(%f,%f); Photon2 (eta,phi)=(%f,%f);\n",eta0,phi0*TMath::RadToDeg(),eta1,phi1*TMath::RadToDeg());
//
//    Int_t ieta0=-1; Int_t iphi0 = 0; Int_t rcu0 = 0;
//    Int_t sm0 = GetModuleNumberCellIndexes(absId0,fCalorimeter, ieta0, iphi0, rcu0);
//    Int_t ieta1=-1; Int_t iphi1 = 0; Int_t rcu1 = 0;
//    Int_t sm1 = GetModuleNumberCellIndexes(absId1,fCalorimeter, ieta1, iphi1, rcu1);
//    
//    printf("Photon1 (id,sm,eta,phi)=(%d,%d,%d,%d), Photon2 (id,sm,eta,phi)=(%d,%d,%d,%d)\n",
//           absId0,sm0,ieta0,iphi0,absId1,sm1,ieta1,iphi1);
//    
//    Int_t ietam0=-1; Int_t iphim0 = 0; Int_t rcum0 = 0; Int_t smm0 = -1 ;
//    if(imax  >= 0) smm0 = GetModuleNumberCellIndexes(list[imax] ,fCalorimeter, ietam0, iphim0, rcum0);
//    Int_t ietam1=-1; Int_t iphim1 = 0; Int_t rcum1 = 0; Int_t smm1 = -1 ;
//    if(imax2 >= 0) smm1 = GetModuleNumberCellIndexes(list[imax2],fCalorimeter, ietam1, iphim1, rcum1);
//    
//    printf("Max (id, sm,eta,phi)=(%d,%d,%d,%d), Max2 (id, sm,eta,phi)=(%d,%d,%d,%d)\n",
//           list[imax],smm0,ietam0,iphim0,list[imax2],smm1,ietam1,iphim1);
//  }

  Int_t inlm = nMax-1;
  if(inlm > 2) inlm = 2;
  
  Bool_t match0  = kFALSE;
  Bool_t match1  = kFALSE;
  Int_t imatch0  = -1;
  Int_t imatch1  = -1;
  if(imax >= 0 && imax2 >=0 && absId0 > 0 && absId1 > 0 )
  {
    if     (absId0 == list[imax] ) { match0 = kTRUE ; imatch0 = imax  ; }
    else if(absId0 == list[imax2]) { match0 = kTRUE ; imatch0 = imax2 ; }
    
    if     (absId1 == list[imax] ) { match1 = kTRUE ; imatch1 = imax  ; }
    else if(absId1 == list[imax2]) { match1 = kTRUE ; imatch1 = imax2 ; }
  }
  
  //printf("primary imatch0 %d, imatch1 %d\n",imatch0,imatch1);

  // If one or the 2 not matched, check with the other MC labels
  // only in case there was a conversion
  
  Int_t   absId0second  = -1;
  Int_t   absId1second  = -1;
  Int_t   secLabel0     = -1;
  Int_t   secLabel1     = -1;
  Int_t   mcLabel0      = -1;
  Int_t   mcLabel1      = -1;
  Bool_t  secOK         = 0;
  Int_t   secpdg        = -999999;
  Int_t   secstatus     = -1;
  Int_t   secgrandLabel = -1;

  if(match0) { secLabel0 = label0 ; mcLabel0 = label0 ; }
  if(match1) { secLabel1 = label1 ; mcLabel1 = label1 ; }
  
  if((!match0 || !match1) && mcindex == kmcPi0Conv)
  {
    for (UInt_t ilab = 0; ilab < cluster->GetNLabels(); ilab++ )
    {
      Int_t mclabel = cluster->GetLabels()[ilab];
      
      //printf("Check label %d - %d\n",ilab,mclabel);
      
      if(mclabel == label0 || mclabel == label1)
      {
        //printf("continue: secLabel %d, label0 %d, label1 %d\n",mclabel,label0,label1);
        if(mclabel == label0 && secLabel0 < 0) { secLabel0 = label0 ; mcLabel0 = label0 ; }
        if(mclabel == label1 && secLabel1 < 0) { secLabel1 = label1 ; mcLabel1 = label1 ; }
        continue ;
      }
      
      //printf("Before while: secLabel0 %d, secLabel1 %d\n",secLabel0,secLabel1);
      
      // match mc label and parent photon
      Int_t tmplabel   = mclabel;
      while((secLabel0 < 0 || secLabel1 < 0) && tmplabel > 0 )
      {
        TLorentzVector mother = GetMCAnalysisUtils()->GetMother(tmplabel,GetReader(),secpdg,secstatus,secOK,secgrandLabel);
        
        //printf("\t \t while secLabel %d, mom %d, granmom %d\n",mclabel,tmplabel,secgrandLabel);
        
        if((secgrandLabel == label0) || (secgrandLabel == label1 ))
        {
          //printf("mcMatch! grand label %d, secLabel %d\n",secgrandLabel, mclabel);
          if(!match0 && mcLabel1 != secgrandLabel) { secLabel0 = mclabel; mcLabel0 = secgrandLabel; }
          if(!match1 && mcLabel0 != secgrandLabel) { secLabel1 = mclabel; mcLabel1 = secgrandLabel; }
        }
        
        //printf("\t GrandMother %d, secLabel0 %d, secLabel1 %d \n",secgrandLabel, secLabel0,secLabel1);

        tmplabel = secgrandLabel;
      }
    }
    
    // Get the position of the found secondaries mother
    if(!match0 && secLabel0 > 0)
    {
      TLorentzVector mother = GetMCAnalysisUtils()->GetMother(secLabel0,GetReader(),secpdg,secstatus,secOK,secgrandLabel);
      
      //Float_t eta = mother.Eta();
      //Float_t phi = mother.Phi();
      //if(phi < 0 ) phi+=TMath::TwoPi();
      //GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(eta, phi, absId0second);
      
      //printf("Secondary MC0 label %d, absId %d E %2.2F eta %2.2f, phi %f\n", secLabel0,absId0second, mother.E(),mother.Eta(),mother.Phi()*TMath::RadToDeg());
      
      if(absId0second == list[imax] ) { match0 = kTRUE ; imatch0 = imax  ; }
      if(absId0second == list[imax2]) { match0 = kTRUE ; imatch0 = imax2 ; }
    }

    if(!match1 && secLabel1 > 0)
    {
      TLorentzVector mother = GetMCAnalysisUtils()->GetMother(secLabel1,GetReader(),secpdg,secstatus,secOK,secgrandLabel);
      
      //Float_t eta = mother.Eta();
      //Float_t phi = mother.Phi();
      //if(phi < 0 ) phi+=TMath::TwoPi();
      //GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(eta, phi, absId1second);
      
      //printf("Secondary MC1 label %d absId %d E %2.2F eta %2.2f, phi %f\n",secLabel1, absId1second, mother.E(),mother.Eta(),mother.Phi()*TMath::RadToDeg());
      
      if(absId1second == list[imax] ) { match1 = kTRUE ; imatch1 = imax  ; }
      if(absId1second == list[imax2]) { match1 = kTRUE ; imatch1 = imax2 ; }
    }

    //printf("secondary label mc0 %d, mc1 %d, imatch0 %d, imatch1 %d\n",secLabel0,secLabel1,imatch0,imatch1);
    
  }
    
  //printf("imatch0 %d, imatch1 %d\n",imatch0,imatch1);
  if( match0 && match1 )
  {
//   if((mass < 0.06 || mass > 1.8) && mcindex==kmcPi0 && noverlaps == 0)
//      printf("a) Both Photons hit local maxima \n");
    
    if(!noverlaps)
    {
      fhMCPi0DecayPhotonHitHighLM          ->Fill(en,nMax);
      fhMCPi0DecayPhotonHitHighLMMass[inlm]->Fill(en,mass);
      if(match0 && imatch0 == imax)
      {
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonHitHighLMDiffELM1[inlm]->Fill(en,(e1-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonHitHighLMDiffELM2[inlm]->Fill(en,(e2-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonHitHighLMDiffELM1vsELM1[inlm]->Fill(e1,(e1-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonHitHighLMDiffELM2vsELM2[inlm]->Fill(e2,(e2-photon1Kine.E())/photon1Kine.E());
      }
      else
      {
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonHitHighLMDiffELM1[inlm]->Fill(en,(e1-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonHitHighLMDiffELM2[inlm]->Fill(en,(e2-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonHitHighLMDiffELM1vsELM1[inlm]->Fill(e1,(e1-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonHitHighLMDiffELM2vsELM2[inlm]->Fill(e2,(e2-photon0Kine.E())/photon0Kine.E());
      }
    }
    else
    {
      fhMCPi0DecayPhotonHitHighLMOverlap          ->Fill(en,nMax);
      fhMCPi0DecayPhotonHitHighLMOverlapMass[inlm]->Fill(en,mass);
      if(match0 && imatch0 == imax )
      {
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1[inlm]->Fill(en,(e1-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2[inlm]->Fill(en,(e2-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1vsELM1[inlm]->Fill(e1,(e1-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2vsELM2[inlm]->Fill(e2,(e2-photon1Kine.E())/photon1Kine.E());
      }
      else
      {
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1[inlm]->Fill(en,(e1-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2[inlm]->Fill(en,(e2-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1vsELM1[inlm]->Fill(e1,(e1-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2vsELM2[inlm]->Fill(e2,(e2-photon0Kine.E())/photon0Kine.E());
      }

    }
    
    return ;
  }
  
  //printf("Any match? photon0 %d, photon1 %d\n",match0,match1);
  //if(!match0 && !match1) printf("WARNING, LM not matched to any photon decay!\n");
  
  //---------------------------------------------
  // Check the adjacent cells to the local maxima
  //---------------------------------------------
  
  if(!match0)
  {
    if(imatch1!=imax  && GetCaloUtils()->AreNeighbours(fCalorimeter,absId0,list[imax]))   { match0 = kTRUE; imatch0 = imax  ; }
    //printf("imax - match0? (%d-%d)=%d, (%d-%d)=%d\n",ieta0,ietam0,ieta0-ietam0, iphi0,iphim0,iphi0-iphim0);
    if(imatch1!=imax2 && GetCaloUtils()->AreNeighbours(fCalorimeter,absId0,list[imax2]) ) { match0 = kTRUE; imatch0 = imax2 ; }
    //printf("imax2 - match0? (%d-%d)=%d, (%d-%d)=%d\n",ieta0,ietam1,ieta0-ietam1, iphi0,iphim1,iphi0-iphim1);
  }
  
  if(!match1)
  {
    if(imatch0!=imax  && GetCaloUtils()->AreNeighbours(fCalorimeter,absId1,list[imax]) ) { match1 = kTRUE; imatch1 = imax  ; }
    //printf("imax - match1? (%d-%d)=%d, (%d-%d)=%d\n",ieta1,ietam0,ieta1-ietam0, iphi1,iphim0,iphi1-iphim0);
  
    if(imatch0!=imax2 && GetCaloUtils()->AreNeighbours(fCalorimeter,absId1,list[imax2])) { match1 = kTRUE; imatch1 = imax2 ; }
    //printf("imax2 - match1? (%d-%d)=%d, (%d-%d)=%d\n",ieta1,ietam1,ieta1-ietam1, iphi1,iphim1,iphi1-iphim1);
  }
    
  //printf("Local Maxima: adjacent0 %d,adjacent1 %d \n",match0,match1);
  
  if(match0 && match1)
  {
//   if((mass < 0.06 || mass > 1.8) && mcindex==kmcPi0 && noverlaps == 0)
//      printf("b) Both Photons hit local maxima or cell adjacent or 2 cells adjacent \n");
    
    if(!noverlaps)
    {
      fhMCPi0DecayPhotonAdjHighLM          ->Fill(en,nMax);
      fhMCPi0DecayPhotonAdjHighLMMass[inlm]->Fill(en,mass);

      if(match0 && imatch0 == imax)
      {
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMDiffELM1[inlm]->Fill(en,(e1-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMDiffELM2[inlm]->Fill(en,(e2-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMDiffELM1vsELM1[inlm]->Fill(e1,(e1-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMDiffELM2vsELM2[inlm]->Fill(e2,(e2-photon1Kine.E())/photon1Kine.E());
      }
      else
      {
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMDiffELM1[inlm]->Fill(en,(e1-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMDiffELM2[inlm]->Fill(en,(e2-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMDiffELM1vsELM1[inlm]->Fill(e1,(e1-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMDiffELM2vsELM2[inlm]->Fill(e2,(e2-photon0Kine.E())/photon0Kine.E());
      }
    }
    else
    {
      fhMCPi0DecayPhotonAdjHighLMOverlap          ->Fill(en,nMax);
      fhMCPi0DecayPhotonAdjHighLMOverlapMass[inlm]->Fill(en,mass);
      if(match0 && imatch0 == imax)
      {
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1[inlm]->Fill(en,(e1-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2[inlm]->Fill(en,(e2-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1vsELM1[inlm]->Fill(e1,(e1-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2vsELM2[inlm]->Fill(e2,(e2-photon1Kine.E())/photon1Kine.E());
      }
      else
      {
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1[inlm]->Fill(en,(e1-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2[inlm]->Fill(en,(e2-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1vsELM1[inlm]->Fill(e1,(e1-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2vsELM2[inlm]->Fill(e2,(e2-photon0Kine.E())/photon0Kine.E());
      }
    }
    
    return;
  }
  
  // Decay photon cells are adjacent?
  
  if( (match0 || match1) && GetCaloUtils()->AreNeighbours(fCalorimeter,absId0,absId1) )
  {
//   if((mass < 0.06 || mass > 1.8) && mcindex==kmcPi0 && noverlaps == 0)
//      printf("c) Both Photons hit a local maxima and in adjacent cells \n");
    if(!noverlaps)
    {
      fhMCPi0DecayPhotonAdjacent          ->Fill(en,nMax);
      fhMCPi0DecayPhotonAdjacentMass[inlm]->Fill(en,mass);
    }
    else
    {
      fhMCPi0DecayPhotonAdjacentOverlap          ->Fill(en,nMax);
      fhMCPi0DecayPhotonAdjacentOverlapMass[inlm]->Fill(en,mass);
    }
    
    return;
  }
  
  //--------------------
  // Other Local maxima
  //--------------------
  
  Bool_t matchMCHitOtherLM = kFALSE;
  if(!match1)
  {
    for(Int_t i = 0; i < nmaxima; i++)
    {
      if(imax!=i && imax2!=i && absId1 == list[i]) { match1 = kTRUE; matchMCHitOtherLM = kTRUE; }
    }
  }
  
  if(!match0)
  {
    for(Int_t i = 0; i < nmaxima; i++)
    {
      if(imax!=i && imax2!=i && absId0 == list[i]) { match0 = kTRUE; matchMCHitOtherLM = kTRUE; }
    }
  }
  
  if(matchMCHitOtherLM)
  {
//   if((mass < 0.06 || mass > 1.8) && mcindex==kmcPi0 && noverlaps == 0)
//      printf("d) One Photon hits a local maxima, the other another not high \n");
    
    if(!noverlaps)
    {
      fhMCPi0DecayPhotonHitOtherLM          ->Fill(en,nMax);
      fhMCPi0DecayPhotonHitOtherLMMass[inlm]->Fill(en,mass);
      if(match0 && imatch0 == imax)
      {
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonHitOtherLMDiffELM1[inlm]->Fill(en,(e1-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonHitOtherLMDiffELM2[inlm]->Fill(en,(e2-photon1Kine.E())/photon1Kine.E());
      }
      else
      {
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonHitOtherLMDiffELM1[inlm]->Fill(en,(e1-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonHitOtherLMDiffELM2[inlm]->Fill(en,(e2-photon0Kine.E())/photon0Kine.E());
      }
    }
    else
    {
      fhMCPi0DecayPhotonHitOtherLMOverlap   ->Fill(en,nMax);
      fhMCPi0DecayPhotonHitOtherLMMass[inlm]->Fill(en,mass);
      if(match0 && imatch0 == imax)
      {
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM1[inlm]->Fill(en,(e1-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM2[inlm]->Fill(en,(e2-photon1Kine.E())/photon1Kine.E());
      }
      else
      {
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM1[inlm]->Fill(en,(e1-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM2[inlm]->Fill(en,(e2-photon0Kine.E())/photon0Kine.E());
      }
    }
    
    return ;
  }
  
  // Adjacent to other maxima
  
  Bool_t adjacentOther1 = kFALSE;
  if(match0)
  {
    for(Int_t i = 0; i < nmaxima; i++)
    {
      Int_t ieta=-1; Int_t iphi = 0; Int_t rcu = 0;
      GetModuleNumberCellIndexes(list[i] ,fCalorimeter, ieta, iphi, rcu);
      
      //printf(" Other Max (eta,phi)=(%d,%d)\n",ieta,iphi);
      
      if(GetCaloUtils()->AreNeighbours(fCalorimeter,absId1,list[i]) ) adjacentOther1 = kTRUE;
      
      //printf("Other Maxima: adjacentOther1 %d\n",adjacentOther1);
    }
  }
  
  Bool_t adjacentOther0 = kFALSE;
  if(match1)
  {
    for(Int_t i = 0; i < nmaxima; i++)
    {
      Int_t ieta=-1; Int_t iphi = 0; Int_t rcu = 0;
      GetModuleNumberCellIndexes(list[i] ,fCalorimeter, ieta, iphi, rcu);
      
      //printf(" Other Max (eta,phi)=(%d,%d)\n",ieta,iphi);
      
      if(GetCaloUtils()->AreNeighbours(fCalorimeter,absId0,list[i]) ) adjacentOther0 = kTRUE;
      
      //printf("Other Maxima: adjacentOther0 %d\n",adjacentOther0);
    }
  }
  
  if((match0 && adjacentOther1) || (match1 && adjacentOther0))
  {
    
//   if((mass < 0.06 || mass > 1.8) && mcindex==kmcPi0 && noverlaps == 0)
//      printf("e) One Photon hits a local maxima, the other another not high, adjacent \n");
    
    if(!noverlaps)
    {
      fhMCPi0DecayPhotonAdjOtherLM       ->Fill(en,nMax);
      fhMCPi0DecayPhotonAdjOtherLMMass[inlm]->Fill(en,mass);
      if(match0 && imatch0 == imax)
      {
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonAdjOtherLMDiffELM1[inlm]->Fill(en,(e1-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonAdjOtherLMDiffELM2[inlm]->Fill(en,(e2-photon1Kine.E())/photon1Kine.E());
      }
      else
      {
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonAdjOtherLMDiffELM1[inlm]->Fill(en,(e1-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonAdjOtherLMDiffELM2[inlm]->Fill(en,(e2-photon0Kine.E())/photon0Kine.E());
      }
    }
    else
    {
      fhMCPi0DecayPhotonAdjOtherLMOverlap          ->Fill(en,nMax);
      fhMCPi0DecayPhotonAdjOtherLMOverlapMass[inlm]->Fill(en,mass);
      if(match0 && imatch0 == imax)
      {
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM1[inlm]->Fill(en,(e1-photon0Kine.E())/photon0Kine.E());
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM2[inlm]->Fill(en,(e2-photon1Kine.E())/photon1Kine.E());
      }
      else
      {
        if(photon1Kine.E()>0)fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM1[inlm]->Fill(en,(e1-photon1Kine.E())/photon1Kine.E());
        if(photon0Kine.E()>0)fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM2[inlm]->Fill(en,(e2-photon0Kine.E())/photon0Kine.E());
      }
    }
    
    return;
  }
  
// if((mass < 0.06 || mass > 1.8) && mcindex==kmcPi0 && noverlaps == 0)
//    printf("f) No hit found \n");
  if(!noverlaps)
  {
    fhMCPi0DecayPhotonHitNoLM          ->Fill(en,nMax);
    fhMCPi0DecayPhotonHitNoLMMass[inlm]->Fill(en,mass);
  }
  else
  {
    fhMCPi0DecayPhotonHitNoLMOverlap          ->Fill(en,nMax);
    fhMCPi0DecayPhotonHitNoLMOverlapMass[inlm]->Fill(en,mass);
  }
  
}

//___________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillAngleHistograms(Int_t   nMax,      Bool_t  matched, Int_t mcIndex,
                                                           Float_t en,        Float_t e1,      Float_t e2,
                                                           Float_t angle,     Float_t mass,
                                                           Float_t anglePrim, Float_t m02,
                                                           Float_t asym,      Int_t   pid,     Int_t noverlaps)
{
  // Fill histograms related to opening angle
  
  Bool_t m02OK = GetCaloPID()->IsInPi0M02Range(en,m02,nMax);
  Bool_t asyOK = GetCaloPID()->IsInPi0SplitAsymmetryRange(en,asym,nMax);
  Bool_t m02On = GetCaloPID()->IsSplitShowerShapeCutOn();
  Bool_t asyOn = GetCaloPID()->IsSplitAsymmetryCutOn();

  Bool_t eCutOK= kFALSE;
  Int_t inlm = nMax-1;
  if(inlm > 2 ) inlm = 2;
  Float_t ensubcut = GetCaloPID()->GetSubClusterEnergyMinimum(inlm);
  if     (ensubcut > 0.1 && ensubcut < e1 && ensubcut < e2 ) eCutOK = kTRUE;
  else if(ensubcut < 0.1)                                    eCutOK = kTRUE;

  if     (nMax==1)
  {
    fhAnglePairNLocMax1[0][matched]->Fill(en,angle);
    
    if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
      fhAnglePairAfterCutsNLocMax1[0][matched]->Fill(en,angle);
    if(pid==AliCaloPID::kPi0)
      fhAnglePairPi0NLocMax1[0][matched]->Fill(en,angle);
    
    if(m02 > 0)
    {
      fhAnglePairOverM02NLocMax1[0][matched]->Fill(en,angle/m02);
      if(noverlaps == 0) fhAnglePairOverM02NLocMax1Overlap0[0][matched]->Fill(en,angle/m02);
    }
    
    if( en > 15 )
    {
      fhAnglePairMassNLocMax1[0][matched]->Fill(mass,angle);
      fhAnglePairM02NLocMax1 [0][matched]->Fill(m02 ,angle);
    }
  }
  else if(nMax==2)
  {
    fhAnglePairNLocMax2[0][matched]->Fill(en,angle);
    
    if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
      fhAnglePairAfterCutsNLocMax2[0][matched]->Fill(en,angle);
    if(pid==AliCaloPID::kPi0)
      fhAnglePairPi0NLocMax2[0][matched]->Fill(en,angle);
    
    if(m02 > 0)
    {
      fhAnglePairOverM02NLocMax2[0][matched]->Fill(en,angle/m02);
      if(noverlaps == 0) fhAnglePairOverM02NLocMax2Overlap0[0][matched]->Fill(angle/m02,en);
    }
    
    if( en > fHistoECut )
    {
      fhAnglePairMassNLocMax2[0][matched]->Fill(mass,angle);
      fhAnglePairM02NLocMax2 [0][matched]->Fill(m02,angle);
    }
  }
  else if(nMax >2)
  {
    fhAnglePairNLocMaxN[0][matched]->Fill(en,angle);
    
    if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
      fhAnglePairAfterCutsNLocMaxN[0][matched]->Fill(en,angle);
    if(pid==AliCaloPID::kPi0)
      fhAnglePairPi0NLocMaxN[0][matched]->Fill(en,angle);
    
    if(m02 > 0)
    {
      fhAnglePairOverM02NLocMaxN[0][matched]->Fill(en,angle/m02);
      if(noverlaps == 0) fhAnglePairOverM02NLocMaxNOverlap0[0][matched]->Fill(angle/m02,en);
    }
    
    if( en > fHistoECut )
    {
      fhAnglePairMassNLocMaxN[0][matched]->Fill(mass,angle);
      fhAnglePairM02NLocMaxN [0][matched]->Fill(m02,angle);
    }
  }
  
  if(IsDataMC() && mcIndex >  0 && mcIndex < 7)
  {
    if     (nMax==1)
    {
      fhAnglePairNLocMax1[mcIndex][matched]->Fill(en,angle);
      if( en > 15 )
      {
        fhAnglePairMassNLocMax1[mcIndex][matched]->Fill(mass,angle);
        fhAnglePairM02NLocMax1 [mcIndex][matched]->Fill(m02,angle);
      }
      if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
        fhAnglePairAfterCutsNLocMax1[mcIndex][matched]->Fill(en,angle);
      if(pid==AliCaloPID::kPi0)
         fhAnglePairPi0NLocMax1[mcIndex][matched]->Fill(en,angle);
      
      if(m02 > 0)
      {
        fhAnglePairOverM02NLocMax1[mcIndex][matched]->Fill(en,angle/m02);
        if(noverlaps == 0) fhAnglePairOverM02NLocMax1Overlap0[mcIndex][matched]->Fill(angle/m02,en);
      }
      
      if((mcIndex == kmcPi0 || mcIndex == kmcPi0Conv) && !matched && anglePrim > 0)
      {
        fhAnglePairPrimPi0RecoNLocMax1->Fill(en,angle/anglePrim);
        if(m02>0)fhAnglePairPrimPi0OverM02NLocMax1->Fill(en,anglePrim/m02);
        if(en > 15) fhAnglePairPrimPi0vsRecoNLocMax1->Fill(anglePrim,angle);

      }
    }
    else if(nMax==2)
    {
      fhAnglePairNLocMax2[mcIndex][matched]->Fill(en,angle);
      if( en > fHistoECut )
      {
        fhAnglePairMassNLocMax2[mcIndex][matched]->Fill(mass,angle);
        fhAnglePairM02NLocMax2 [mcIndex][matched]->Fill(m02 ,angle);
      }
      
      if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
        fhAnglePairAfterCutsNLocMax2[mcIndex][matched]->Fill(en,angle);
      if(pid==AliCaloPID::kPi0)
        fhAnglePairPi0NLocMax2[mcIndex][matched]->Fill(en,angle);
      
      if(m02 > 0)
      {
        fhAnglePairOverM02NLocMax2[mcIndex][matched]->Fill(en,angle/m02);
        if(noverlaps == 0) fhAnglePairOverM02NLocMax2Overlap0[mcIndex][matched]->Fill(angle/m02,en);
      }
      
      if((mcIndex == kmcPi0 || mcIndex == kmcPi0Conv) && !matched && anglePrim > 0)
      {
        fhAnglePairPrimPi0RecoNLocMax2->Fill(en,angle/anglePrim);
        if(m02>0)fhAnglePairPrimPi0OverM02NLocMax2->Fill(en,anglePrim/m02);
        if(en > 10) fhAnglePairPrimPi0vsRecoNLocMax2->Fill(anglePrim,angle);
      }
    }
    else if(nMax >2)
    {
      fhAnglePairNLocMaxN[mcIndex][matched]->Fill(en,angle);
      if( en > fHistoECut )
      {
        fhAnglePairMassNLocMaxN[mcIndex][matched]->Fill(mass,angle);
        fhAnglePairM02NLocMaxN [mcIndex][matched]->Fill(m02 ,angle);
      }
      if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
        fhAnglePairAfterCutsNLocMaxN[mcIndex][matched]->Fill(en,angle);
      if(pid==AliCaloPID::kPi0)
        fhAnglePairPi0NLocMaxN[mcIndex][matched]->Fill(en,angle);
      
      if(m02 > 0)
      {
        fhAnglePairOverM02NLocMaxN[mcIndex][matched]->Fill(en,angle/m02);
        if(noverlaps == 0) fhAnglePairOverM02NLocMaxNOverlap0[mcIndex][matched]->Fill(angle/m02,en);
      }
      
      if((mcIndex == kmcPi0 || mcIndex == kmcPi0Conv) && !matched && anglePrim > 0)
      {
        fhAnglePairPrimPi0RecoNLocMaxN->Fill(en,angle/anglePrim);
        if(m02>0)fhAnglePairPrimPi0OverM02NLocMaxN->Fill(en,anglePrim/m02);
        if(en > 10) fhAnglePairPrimPi0vsRecoNLocMaxN->Fill(anglePrim,angle);
      }
    }
    
  }
  
}

//______________________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillArmenterosHistograms(Int_t nMax, Int_t ebin, Int_t mcIndex,
                                                                Float_t en, TLorentzVector g1, TLorentzVector g2,
                                                                Float_t m02, Int_t pid)
{
  // Fill Armeteros type histograms
  
  // Get pTArm and AlphaArm
  TLorentzVector pi0 = g1+g2;
  Float_t momentumSquaredMother = pi0.P()*pi0.P();
  Float_t momentumDaughter1AlongMother = 0.;
  Float_t momentumDaughter2AlongMother = 0.;

  if (momentumSquaredMother > 0.)
  {
    momentumDaughter1AlongMother = (g1.Px()*pi0.Px() + g1.Py()*pi0.Py()+ g1.Pz()*pi0.Pz()) / sqrt(momentumSquaredMother);
    momentumDaughter2AlongMother = (g2.Px()*pi0.Px() + g2.Py()*pi0.Py()+ g2.Pz()*pi0.Pz()) / sqrt(momentumSquaredMother);
  }

  Float_t momentumSquaredDaughter1 = g1.P()*g1.P();
  Float_t ptArmSquared = momentumSquaredDaughter1 - momentumDaughter1AlongMother*momentumDaughter1AlongMother;
  
  Float_t pTArm = 0.;
  if (ptArmSquared > 0.)
    pTArm = sqrt(ptArmSquared);
  
  Float_t alphaArm = 0.;
  if(momentumDaughter1AlongMother +momentumDaughter2AlongMother > 0)
    alphaArm = (momentumDaughter1AlongMother -momentumDaughter2AlongMother) / (momentumDaughter1AlongMother + momentumDaughter2AlongMother);
  
  Float_t asym = TMath::Abs( g1.Energy()-g2.Energy() )/( g1.Energy()+g2.Energy() ) ;
  
   if(GetDebug() > 2 ) Info("FillArmenterosHistograms()","E %f, alphaArm %f, pTArm %f\n",en,alphaArm,pTArm);
  
  Bool_t m02OK = GetCaloPID()->IsInPi0M02Range(en,m02,nMax);
  Bool_t asyOK = GetCaloPID()->IsInPi0SplitAsymmetryRange(en,asym,nMax);
  Bool_t m02On = GetCaloPID()->IsSplitShowerShapeCutOn();
  Bool_t asyOn = GetCaloPID()->IsSplitAsymmetryCutOn();
  
  Bool_t eCutOK= kFALSE;
  Int_t inlm = nMax-1;
  if(inlm > 2 ) inlm = 2;
  Float_t ensubcut = GetCaloPID()->GetSubClusterEnergyMinimum(inlm);
  if     (ensubcut > 0.1 && ensubcut < g1.E() && ensubcut < g2.E() ) eCutOK = kTRUE;
  else if(ensubcut < 0.1)                                            eCutOK = kTRUE;

  
  if     (nMax==1)
  {
    fhArmNLocMax1[0][ebin]->Fill(alphaArm,pTArm);
    if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
      fhArmAfterCutsNLocMax1[0][ebin]->Fill(alphaArm,pTArm);
    if(pid==AliCaloPID::kPi0)
      fhArmPi0NLocMax1[0][ebin]->Fill(alphaArm,pTArm);
  }
  else if(nMax==2)
  {
    fhArmNLocMax2[0][ebin]->Fill(alphaArm,pTArm);
    if((m02OK && asyOK) && (asyOn || m02On))
      fhArmAfterCutsNLocMax2[0][ebin]->Fill(alphaArm,pTArm);
    if(pid==AliCaloPID::kPi0)
      fhArmPi0NLocMax2[0][ebin]->Fill(alphaArm,pTArm);
  }
  else if(nMax >2)
  {
    fhArmNLocMaxN[0][ebin]->Fill(alphaArm,pTArm);
    if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
      fhArmAfterCutsNLocMaxN[0][ebin]->Fill(alphaArm,pTArm);
    if(pid==AliCaloPID::kPi0)
      fhArmPi0NLocMaxN[0][ebin]->Fill(alphaArm,pTArm);
  }

  if(IsDataMC() && mcIndex >  0 && mcIndex < 7)
  {
    if     (nMax==1)
    {
      fhArmNLocMax1[mcIndex][ebin]->Fill(alphaArm,pTArm);
      if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
        fhArmAfterCutsNLocMax1[mcIndex][ebin]->Fill(alphaArm,pTArm);
      if(pid==AliCaloPID::kPi0)
        fhArmPi0NLocMax1[mcIndex][ebin]->Fill(alphaArm,pTArm);
    }
    else if(nMax==2)
    {
      fhArmNLocMax2[mcIndex][ebin]->Fill(alphaArm,pTArm);
      if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
        fhArmAfterCutsNLocMax2[mcIndex][ebin]->Fill(alphaArm,pTArm);
      if(pid==AliCaloPID::kPi0)
        fhArmPi0NLocMax2[mcIndex][ebin]->Fill(alphaArm,pTArm);
    }
    else if(nMax >2)
    {
      fhArmNLocMaxN[mcIndex][ebin]->Fill(alphaArm,pTArm);
      if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
        fhArmAfterCutsNLocMaxN[mcIndex][ebin]->Fill(alphaArm,pTArm);
      if(pid==AliCaloPID::kPi0)
        fhArmPi0NLocMaxN[mcIndex][ebin]->Fill(alphaArm,pTArm);
    }  
  
  }
  
}

//______________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillThetaStarHistograms(Int_t nMax, Bool_t matched, Int_t mcIndex,
                                                               Float_t en, TLorentzVector g1, TLorentzVector g2,
                                                               Float_t m02, Int_t pid)
{
  // Fill cos Theta^star histograms

  
  // Get cos Theta^star
  TLorentzVector pi0 = g1+g2;
  TLorentzVector g1Boost = g1;
  g1Boost.Boost(-pi0.BoostVector());
  Float_t  cosThStar=TMath::Cos(g1Boost.Vect().Angle(pi0.Vect()));
  
  Float_t asym = TMath::Abs( g1.Energy()-g2.Energy() )/( g1.Energy()+g2.Energy() ) ;

  Bool_t m02OK = GetCaloPID()->IsInPi0M02Range(en,m02,nMax);
  Bool_t asyOK = GetCaloPID()->IsInPi0SplitAsymmetryRange(en,asym,nMax);
  Bool_t m02On = GetCaloPID()->IsSplitShowerShapeCutOn();
  Bool_t asyOn = GetCaloPID()->IsSplitAsymmetryCutOn();
  
  Bool_t eCutOK= kFALSE;
  Int_t inlm = nMax-1;
  if(inlm > 2 ) inlm = 2;
  Float_t ensubcut = GetCaloPID()->GetSubClusterEnergyMinimum(inlm);
  if     (ensubcut > 0.1 && ensubcut < g1.E() && ensubcut < g2.E() ) eCutOK = kTRUE;
  else if(ensubcut < 0.1)                                            eCutOK = kTRUE;

  //printf("Reco cos %f, asy %f\n",cosThStar,asym);
  
  if     (nMax==1)
  {
    fhCosThStarNLocMax1[0][matched]->Fill(en,cosThStar);
    
    if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
      fhCosThStarAfterCutsNLocMax1[0][matched]->Fill(en,cosThStar);
    if(pid==AliCaloPID::kPi0)
      fhCosThStarPi0NLocMax1[0][matched]->Fill(en,cosThStar);
  }
  else if(nMax==2)
  {
    fhCosThStarNLocMax2[0][matched]->Fill(en,cosThStar);
    
    if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
      fhCosThStarAfterCutsNLocMax2[0][matched]->Fill(en,cosThStar);
    if(pid==AliCaloPID::kPi0)
      fhCosThStarPi0NLocMax2[0][matched]->Fill(en,cosThStar);
  }
  else if(nMax >2)
  {
    fhCosThStarNLocMaxN[0][matched]->Fill(en,cosThStar);
    
    if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
      fhCosThStarAfterCutsNLocMaxN[0][matched]->Fill(en,cosThStar);
    if(pid==AliCaloPID::kPi0)
      fhCosThStarPi0NLocMaxN[0][matched]->Fill(en,cosThStar);
  }
  
  if(IsDataMC() && mcIndex >  0 && mcIndex < 7)
  {
    if     (nMax==1)
    {
      fhCosThStarNLocMax1[mcIndex][matched]->Fill(en,cosThStar);
      
      if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
        fhCosThStarAfterCutsNLocMax1[mcIndex][matched]->Fill(en,cosThStar);
      if(pid==AliCaloPID::kPi0)
        fhCosThStarPi0NLocMax1[mcIndex][matched]->Fill(en,cosThStar);
    }
    else if(nMax==2)
    {
      fhCosThStarNLocMax2[mcIndex][matched]->Fill(en,cosThStar);
      
      if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
        fhCosThStarAfterCutsNLocMax2[mcIndex][matched]->Fill(en,cosThStar);
      if(pid==AliCaloPID::kPi0)
        fhCosThStarPi0NLocMax2[mcIndex][matched]->Fill(en,cosThStar);
    }
    else if(nMax >2)
    {
      fhCosThStarNLocMaxN[mcIndex][matched]->Fill(en,cosThStar);
      
      if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
        fhCosThStarAfterCutsNLocMaxN[mcIndex][matched]->Fill(en,cosThStar);
      if(pid==AliCaloPID::kPi0)
        fhCosThStarPi0NLocMaxN[mcIndex][matched]->Fill(en,cosThStar);
    }
    
  }

}

//__________________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillEBinHistograms(Int_t   ebin     , Int_t   nMax, Int_t mcindex,
                                                          Float_t splitFrac, Float_t mass, Float_t asym, Float_t l0)
{
  // Fill some histograms integrating in few energy bins
    
  if     (nMax==1)
  {
    fhMassSplitEFractionNLocMax1Ebin[0][ebin]->Fill(splitFrac,  mass);
    if(IsDataMC() && mcindex > 0 && mcindex < 7)fhMassSplitEFractionNLocMax1Ebin[mcindex][ebin]->Fill(splitFrac,  mass);
    
    fhMassM02NLocMax1Ebin    [ebin]->Fill(l0  ,  mass );
    fhMassAsyNLocMax1Ebin    [ebin]->Fill(asym,  mass );
  }
  else if(nMax==2)
  {
    fhMassSplitEFractionNLocMax2Ebin[0][ebin]->Fill(splitFrac,  mass);
    if(IsDataMC() && mcindex > 0 && mcindex < 7)fhMassSplitEFractionNLocMax2Ebin[mcindex][ebin]->Fill(splitFrac,  mass);
    
    fhMassM02NLocMax2Ebin    [ebin]->Fill(l0  ,  mass );
    fhMassAsyNLocMax2Ebin    [ebin]->Fill(asym,  mass );
  }
  else if(nMax > 2 )
  {
    fhMassSplitEFractionNLocMaxNEbin[0][ebin]->Fill(splitFrac,  mass);
    if(IsDataMC() && mcindex > 0 && mcindex < 7)fhMassSplitEFractionNLocMaxNEbin[mcindex][ebin]->Fill(splitFrac,  mass);
    
    fhMassM02NLocMaxNEbin    [ebin]->Fill(l0  ,  mass );
    fhMassAsyNLocMaxNEbin    [ebin]->Fill(asym,  mass );
  }
  
}

//________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillHistograms1(Float_t en,     Float_t e1,     Float_t e2,
                                                       Int_t nMax,     Float_t mass,   Float_t l0,
                                                       Float_t eta,    Float_t phi,
                                                       Bool_t matched, Int_t mcindex)
{
  // Fill histograms for clusters before any selection after spliting
  
  Float_t splitFrac = (e1+e2)/en;
  
  Float_t asym = -10;
  if(e1+e2>0) asym = (e1-e2)/(e1+e2);
  
  fhNLocMax   [0][matched]->Fill(en,nMax);
  fhLM1NLocMax[0][matched]->Fill(e1,nMax);
  fhLM2NLocMax[0][matched]->Fill(e2,nMax);
  fhSplitClusterENLocMax[0][matched]->Fill(e1,nMax);
  fhSplitClusterENLocMax[0][matched]->Fill(e2,nMax);
  
  if(IsDataMC() && mcindex > 0 && mcindex < 7)
  {
    fhNLocMax   [mcindex][matched]->Fill(en,nMax);
    fhLM1NLocMax[mcindex][matched]->Fill(e1,nMax);
    fhLM2NLocMax[mcindex][matched]->Fill(e2,nMax);
    fhSplitClusterENLocMax[mcindex][matched]->Fill(e1,nMax);
    fhSplitClusterENLocMax[mcindex][matched]->Fill(e2,nMax);
  }
  
  if     ( nMax == 1  )
  {
    fhM02NLocMax1[0][matched]->Fill(en,l0) ;
    fhSplitEFractionNLocMax1[0][matched]->Fill(en,splitFrac) ;
    
    if(IsDataMC() && mcindex > 0 && mcindex < 7)
    {
      fhM02NLocMax1[mcindex][matched]->Fill(en,l0) ;
      fhSplitEFractionNLocMax1[mcindex][matched]->Fill(en,splitFrac) ;
    }
    
    if(en > fHistoECut)
    {
      fhMassM02NLocMax1[0][matched]->Fill(l0, mass);
      if( IsDataMC() && mcindex > 0 && mcindex < 7 ) fhMassM02NLocMax1[mcindex][matched]->Fill(l0, mass);
      
      fhSplitEFractionvsAsyNLocMax1[matched]->Fill(asym,splitFrac) ;
      if(!matched)fhClusterEtaPhiNLocMax1->Fill(eta,phi);
    }
  }
  else if( nMax == 2  )
  {
    fhM02NLocMax2[0][matched]->Fill(en,l0) ;
    fhSplitEFractionNLocMax2[0][matched]->Fill(en,splitFrac) ;
    
    if(IsDataMC() && mcindex > 0 && mcindex < 7)
    {
      fhM02NLocMax2[mcindex][matched]->Fill(en,l0) ;
      fhSplitEFractionNLocMax2[mcindex][matched]->Fill(en,splitFrac) ;
    }
    
    if(en > fHistoECut)
    {
      fhMassM02NLocMax2[0][matched]->Fill(l0,  mass );
      if( IsDataMC() && mcindex > 0 && mcindex < 7 ) fhMassM02NLocMax2[mcindex][matched]->Fill(l0,mass);
      
      fhSplitEFractionvsAsyNLocMax2[matched]->Fill(asym,splitFrac) ;
      if(!matched)fhClusterEtaPhiNLocMax2->Fill(eta,phi);
    }
  }
  else if( nMax >= 3  )
  {
    fhM02NLocMaxN[0][matched]->Fill(en,l0) ;
    fhSplitEFractionNLocMaxN[0][matched]->Fill(en,splitFrac) ;
    
    if(IsDataMC() && mcindex > 0 && mcindex < 7)
    {
      fhM02NLocMaxN[mcindex][matched]->Fill(en,l0) ;
      fhSplitEFractionNLocMaxN[mcindex][matched]->Fill(en,splitFrac) ;
    }
    
    if(en > fHistoECut)
    {
      
      fhMassM02NLocMaxN[0][matched]->Fill(l0,mass);
      if( IsDataMC() && mcindex > 0 && mcindex < 7 ) fhMassM02NLocMaxN[mcindex][matched]->Fill(l0,mass);
      
      fhSplitEFractionvsAsyNLocMaxN[matched]->Fill(asym,splitFrac) ;
      if(!matched)fhClusterEtaPhiNLocMaxN->Fill(eta,phi);
    }
  }
  
  
}

//________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillHistograms2(Float_t en,     Float_t eprim,
                                                       Float_t e1,     Float_t e2,
                                                       Int_t nMax,     Float_t mass,   Float_t l0,
                                                       Bool_t matched, Int_t mcindex)
{
  // Fill histograms for clusters passing the first M02 selection
  
  Float_t efrac      = eprim/en;
  Float_t efracSplit = 0;
  if(e1+e2 > 0) efracSplit = eprim/(e1+e2);
  
  Float_t splitFrac = (e1+e2)/en;
  
  Float_t asym = -10;
  if(e1+e2>0) asym = (e1-e2)/(e1+e2);
  
  Int_t inlm = nMax-1;
  if(inlm > 2) inlm = 2;
  Float_t splitFracMin = GetCaloPID()->GetSplitEnergyFractionMinimum(inlm) ;
  
  Bool_t m02OK = GetCaloPID()->IsInPi0M02Range(en,l0,nMax);
  Bool_t asyOK = GetCaloPID()->IsInPi0SplitAsymmetryRange(en,asym,nMax);
  Bool_t m02On = GetCaloPID()->IsSplitShowerShapeCutOn();
  Bool_t asyOn = GetCaloPID()->IsSplitAsymmetryCutOn();
  
  Bool_t eCutOK = kFALSE;
  Float_t ensubcut = GetCaloPID()->GetSubClusterEnergyMinimum(inlm);
  if     (ensubcut > 0.1 && ensubcut < e1 && ensubcut < e2 ) eCutOK = kTRUE;
  else if(ensubcut < 0.1)                                    eCutOK = kTRUE;

  //printf("splitFracMin %f, val %f, m02ok %d, asyok %d, m02On %d, asyOn %d, ecutOK %d\n",splitFracMin,splitFrac,m02OK,asyOK,m02On,asyOn,eCutOK);
  
  if(m02On && m02OK)
  {
    fhNLocMaxM02Cut   [0][matched]->Fill(en,nMax);
    fhLM1NLocMaxM02Cut[0][matched]->Fill(e1,nMax);
    fhLM2NLocMaxM02Cut[0][matched]->Fill(e2,nMax);
    if(IsDataMC() && mcindex > 0 && mcindex < 7)
    {
      fhNLocMaxM02Cut   [mcindex][matched]->Fill(en,nMax);
      fhLM1NLocMaxM02Cut[mcindex][matched]->Fill(e1,nMax);
      fhLM2NLocMaxM02Cut[mcindex][matched]->Fill(e2,nMax);
    }
  }
  
  if     (nMax==1)
  {
    fhMassNLocMax1[0][matched]->Fill(en,mass );
    fhAsymNLocMax1[0][matched]->Fill(en,asym );
    fhMassSplitENLocMax1[0][matched]->Fill(e1+e2,mass);
    
    // Effect of cuts in mass histograms

    if(!matched && asyOK && asyOn )
    {
      fhMassAsyCutNLocMax1->Fill(en,mass);
      fhM02AsyCutNLocMax1 ->Fill(en,l0 );
    }
    
    if(!matched && m02OK && m02On )
    {
      fhMassM02CutNLocMax1->Fill(en,mass);
      fhAsymM02CutNLocMax1->Fill(en,asym );
      if(splitFrac > splitFracMin && fhMassSplitECutNLocMax1) fhMassSplitECutNLocMax1->Fill(en,mass );
    } 
    
    if(!matched && eCutOK && ensubcut > 0.1)
    {
      fhMassEnCutNLocMax1->Fill(en,mass );
      fhM02EnCutNLocMax1 ->Fill(en,l0   );
      fhAsymEnCutNLocMax1->Fill(en,asym );
      fhSplitEFracEnCutNLocMax1->Fill(en,splitFrac );
    }
    
    if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
    {
      fhSplitEFractionAfterCutsNLocMax1[0][matched]->Fill(en,splitFrac);
      if(splitFrac > splitFracMin)
      {
        fhMassAfterCutsNLocMax1[0][matched]->Fill(en,mass);
        fhMassSplitEAfterCutsNLocMax1[0][matched]->Fill(e1+e2,mass);
      }
      if(!matched && IsDataMC() && fFillMCHisto && mcindex==kmcPi0)
      {
        fhMCGenFracAfterCutsNLocMax1MCPi0      ->Fill(en   ,  efrac     );
        fhMCGenSplitEFracAfterCutsNLocMax1MCPi0->Fill(en   ,  efracSplit);
      }
    }
  }
  else if( nMax == 2 )
  {
    fhMassNLocMax2[0][matched]->Fill(en,mass );
    fhAsymNLocMax2[0][matched]->Fill(en,asym );
    fhMassSplitENLocMax2[0][matched]->Fill(e1+e2,mass);

    // Effect of cuts in mass histograms
    
    if(!matched && asyOK && asyOn )
    {
      fhMassAsyCutNLocMax2->Fill(en,mass);
      fhM02AsyCutNLocMax2 ->Fill(en,l0 );
    }
    
    if(!matched && m02OK && m02On )
    {
      fhMassM02CutNLocMax2->Fill(en,mass);
      fhAsymM02CutNLocMax2->Fill(en,asym );
      if(splitFrac > splitFracMin && fhMassSplitECutNLocMax2) fhMassSplitECutNLocMax2->Fill(en,mass );
    } 
    
    if(!matched && eCutOK && ensubcut > 0.1)
    {
      fhMassEnCutNLocMax2->Fill(en,mass );
      fhM02EnCutNLocMax2 ->Fill(en,l0   );
      fhAsymEnCutNLocMax2->Fill(en,asym );
      fhSplitEFracEnCutNLocMax2->Fill(en,splitFrac );
    }
    
    if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
    {
      fhSplitEFractionAfterCutsNLocMax2[0][matched]->Fill(en,splitFrac);
      if(splitFrac > splitFracMin)
      {
        fhMassAfterCutsNLocMax2[0][matched]->Fill(en,mass);
        fhMassSplitEAfterCutsNLocMax2[0][matched]->Fill(e1+e2,mass);
      }
      
      if(!matched && IsDataMC() && fFillMCHisto && mcindex==kmcPi0)
      {
        fhMCGenFracAfterCutsNLocMax2MCPi0      ->Fill(en   ,  efrac     );
        fhMCGenSplitEFracAfterCutsNLocMax2MCPi0->Fill(en   ,  efracSplit);
      }
    }
  }
  else if( nMax > 2 )
  {
    fhMassNLocMaxN[0][matched]->Fill(en,mass);
    fhAsymNLocMaxN[0][matched]->Fill(en,asym);
    fhMassSplitENLocMaxN[0][matched]->Fill(e1+e2,mass);

    // Effect of cuts in mass histograms
    
    if(!matched && asyOK && asyOn )
    {
      fhMassAsyCutNLocMaxN->Fill(en,mass);
      fhM02AsyCutNLocMaxN ->Fill(en,l0 );
    }
    
    if(!matched && m02OK && m02On )
    {
      fhMassM02CutNLocMaxN->Fill(en,mass);
      fhAsymM02CutNLocMaxN->Fill(en,asym );
      if(splitFrac > splitFracMin && fhMassSplitECutNLocMaxN) fhMassSplitECutNLocMaxN->Fill(en,mass );
    } 
    
    if(!matched && eCutOK && ensubcut > 0.1 )
    {
      fhMassEnCutNLocMaxN->Fill(en,mass );
      fhM02EnCutNLocMaxN ->Fill(en,l0   );
      fhAsymEnCutNLocMaxN->Fill(en,asym );
      fhSplitEFracEnCutNLocMaxN->Fill(en,splitFrac );
    }
    
    if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
    {
      fhSplitEFractionAfterCutsNLocMaxN[0][matched]->Fill(en,splitFrac);
      if(splitFrac > splitFracMin)
      {
        fhMassAfterCutsNLocMaxN[0][matched]->Fill(en,mass);
        fhMassSplitEAfterCutsNLocMaxN[0][matched]->Fill(e1+e2,mass);
      }
      
      if(!matched && IsDataMC() && fFillMCHisto && mcindex==kmcPi0)
      {
        fhMCGenFracAfterCutsNLocMaxNMCPi0      ->Fill(en   ,  efrac     );
        fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0->Fill(en   ,  efracSplit);
      }
    }
  }
  
  if(IsDataMC() && mcindex > 0 && mcindex < 7)
  {
    if     (nMax==1)
    {
      fhMassNLocMax1[mcindex][matched]->Fill(en,mass);
      fhAsymNLocMax1[mcindex][matched]->Fill(en,asym);
      fhMassSplitENLocMax1[mcindex][matched]->Fill(e1+e2,mass);

      if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
      {
        fhSplitEFractionAfterCutsNLocMax1[mcindex][matched]->Fill(en,splitFrac);
        if(splitFrac > splitFracMin)
        {
          fhMassAfterCutsNLocMax1[mcindex][matched]->Fill(en,mass);
          fhMassSplitEAfterCutsNLocMax1[mcindex][matched]->Fill(e1+e2,mass);
        }
      }
    }
    else if(nMax==2)
    {
      fhMassNLocMax2[mcindex][matched]->Fill(en,mass);
      fhAsymNLocMax2[mcindex][matched]->Fill(en,asym);
      fhMassSplitENLocMax2[mcindex][matched]->Fill(e1+e2,mass);

      if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
      {
        fhSplitEFractionAfterCutsNLocMax2[mcindex][matched]->Fill(en,splitFrac);
        if(splitFrac > splitFracMin)
        {
          fhMassAfterCutsNLocMax2[mcindex][matched]->Fill(en,mass);
          fhMassSplitEAfterCutsNLocMax2[mcindex][matched]->Fill(e1+e2,mass);
        }

      }
    }
    else if(nMax >2)
    {
      fhMassNLocMaxN[mcindex][matched]->Fill(en,mass);
      fhAsymNLocMaxN[mcindex][matched]->Fill(en,asym);
      fhMassSplitENLocMaxN[mcindex][matched]->Fill(e1+e2,mass);

      if((m02OK && asyOK) && (asyOn || m02On) && eCutOK)
      {
        fhSplitEFractionAfterCutsNLocMaxN[mcindex][matched]->Fill(en,splitFrac);
        if(splitFrac > splitFracMin)
        {
          fhMassAfterCutsNLocMaxN[mcindex][matched]->Fill(en,mass);
          fhMassSplitEAfterCutsNLocMaxN[mcindex][matched]->Fill(e1+e2,mass);
        }
      }
    }
  }//Work with MC truth
}


//_________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillIdPi0Histograms(Float_t en,     Float_t e1,     Float_t e2,
                                                           Int_t nc,       Int_t nMax,     Float_t t12diff,
                                                           Float_t mass,   Float_t l0,
                                                           Float_t eta,    Float_t phi,
                                                           Bool_t matched, Int_t mcindex)
{
  // Fill histograms for clusters passing the pi0 selection
  
  Float_t asym = -10;
  if(e1+e2>0) asym = (e1-e2)/(e1+e2);
  
  fhNLocMaxIdPi0   [0][matched]->Fill(en,nMax);
  fhLM1NLocMaxIdPi0[0][matched]->Fill(e1,nMax);
  fhLM2NLocMaxIdPi0[0][matched]->Fill(e2,nMax);
  
  fhSplitClusterEPi0NLocMax[0][matched]->Fill(e1,nMax);
  fhSplitClusterEPi0NLocMax[0][matched]->Fill(e2,nMax);
  
  if(IsDataMC() && mcindex > 0 && mcindex < 7)
  {
    fhSplitClusterEPi0NLocMax[mcindex][matched]->Fill(e1,nMax);
    fhSplitClusterEPi0NLocMax[mcindex][matched]->Fill(e2,nMax);
  }
  
  if     (nMax==1)
  {
    fhM02Pi0NLocMax1 [0][matched]->Fill(en,l0);
    fhMassPi0NLocMax1[0][matched]->Fill(en,mass);
    fhAsyPi0NLocMax1 [0][matched]->Fill(en,asym);
    fhMassSplitEPi0NLocMax1[0][matched]->Fill(e1+e2,mass);
    if(fFillNCellHisto) fhNCellPi0NLocMax1[0][matched]->Fill(en,nc);
    
    if(!matched)
    {
      if(fFillHighMultHisto)
      {
        fhCentralityPi0NLocMax1->Fill(en,GetEventCentrality()) ;
        fhEventPlanePi0NLocMax1->Fill(en,GetEventPlaneAngle()) ;
      }
      if(en > fHistoECut)fhPi0EtaPhiNLocMax1->Fill(eta,phi);
      fhPi0EPairDiffTimeNLM1->Fill(e1+e2,t12diff);
    }
  }
  else if(nMax==2)
  {
    fhM02Pi0NLocMax2 [0][matched]->Fill(en,l0);
    fhMassPi0NLocMax2[0][matched]->Fill(en,mass);
    fhAsyPi0NLocMax2 [0][matched]->Fill(en,asym);
    fhMassSplitEPi0NLocMax2[0][matched]->Fill(e1+e2,mass);
    if(fFillNCellHisto) fhNCellPi0NLocMax2[0][matched]->Fill(en,nc);
    
    if(!matched)
    {
      if(fFillHighMultHisto)
      {
        fhCentralityPi0NLocMax2->Fill(en,GetEventCentrality()) ;
        fhEventPlanePi0NLocMax2->Fill(en,GetEventPlaneAngle()) ;
      }
      if(en > fHistoECut)fhPi0EtaPhiNLocMax2->Fill(eta,phi);
      fhPi0EPairDiffTimeNLM2->Fill(e1+e2,t12diff);
    }
  }
  else if(nMax >2)
  {
    fhM02Pi0NLocMaxN [0][matched]->Fill(en,l0);
    fhMassPi0NLocMaxN[0][matched]->Fill(en,mass);
    fhAsyPi0NLocMaxN [0][matched]->Fill(en,asym);
    fhMassSplitEPi0NLocMaxN[0][matched]->Fill(e1+e2,mass);
    if(fFillNCellHisto) fhNCellPi0NLocMaxN[0][matched]->Fill(en,nc);
    
    if(!matched)
    {
      if(fFillHighMultHisto)
      {
        fhCentralityPi0NLocMaxN->Fill(en,GetEventCentrality()) ;
        fhEventPlanePi0NLocMaxN->Fill(en,GetEventPlaneAngle()) ;
      }
      if(en > fHistoECut)fhPi0EtaPhiNLocMaxN->Fill(eta,phi);
      fhPi0EPairDiffTimeNLMN->Fill(e1+e2,t12diff);
    }
  }
  
  if(IsDataMC() && mcindex > 0 && mcindex < 7)
  {
    fhNLocMaxIdPi0   [mcindex][matched]->Fill(en,nMax);
    fhLM1NLocMaxIdPi0[mcindex][matched]->Fill(e1,nMax);
    fhLM2NLocMaxIdPi0[mcindex][matched]->Fill(e2,nMax);
    
    if     (nMax==1)
    {
      fhM02Pi0NLocMax1 [mcindex][matched]->Fill(en,l0);
      fhMassPi0NLocMax1[mcindex][matched]->Fill(en,mass);
      fhAsyPi0NLocMax1 [mcindex][matched]->Fill(en,asym);
      fhMassSplitEPi0NLocMax1[mcindex][matched]->Fill(e1+e2,mass);
      if(fFillNCellHisto) fhNCellPi0NLocMax1[mcindex][matched]->Fill(en,nc);
      
    }
    else if(nMax==2)
    {
      fhM02Pi0NLocMax2 [mcindex][matched]->Fill(en,l0);
      fhMassPi0NLocMax2[mcindex][matched]->Fill(en,mass);
      fhAsyPi0NLocMax2 [mcindex][matched]->Fill(en,asym);
      fhMassSplitEPi0NLocMax2[mcindex][matched]->Fill(e1+e2,mass);
      if(fFillNCellHisto) fhNCellPi0NLocMax2[mcindex][matched]->Fill(en,nc);
    }
    else if(nMax >2)
    {
      fhM02Pi0NLocMaxN [mcindex][matched]->Fill(en,l0);
      fhMassPi0NLocMaxN[mcindex][matched]->Fill(en,mass);
      fhAsyPi0NLocMaxN [mcindex][matched]->Fill(en,asym);
      fhMassSplitEPi0NLocMaxN[mcindex][matched]->Fill(e1+e2,mass);
      if(fFillNCellHisto) fhNCellPi0NLocMaxN[mcindex][matched]->Fill(en,nc);
    }
  }//Work with MC truth
}

//______________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillIdEtaHistograms(Float_t en,     Float_t e1,  Float_t e2,
                                                           Int_t nc,       Int_t nMax,  Float_t t12diff,
                                                           Float_t mass,   Float_t l0,
                                                           Float_t eta,    Float_t phi,
                                                           Bool_t matched, Int_t mcindex)
{
  // Fill histograms for clusters passing the eta selection
  
  Float_t asym = -10;
  if(e1+e2>0) asym = (e1-e2)/(e1+e2);
  
  if     (nMax==1)
  {
    fhM02EtaNLocMax1 [0][matched]->Fill(en,l0);
    fhMassEtaNLocMax1[0][matched]->Fill(en,mass);
    fhAsyEtaNLocMax1 [0][matched]->Fill(en,asym);
    if(fFillNCellHisto) fhNCellEtaNLocMax1[0][matched]->Fill(en,nc);
    
    if(!matched)
    {
      if(fFillHighMultHisto)
      {
        fhCentralityEtaNLocMax1->Fill(en,GetEventCentrality()) ;
        fhEventPlaneEtaNLocMax1->Fill(en,GetEventPlaneAngle()) ;
      }
      if(en > fHistoECut)fhEtaEtaPhiNLocMax1->Fill(eta,phi);
      fhEtaEPairDiffTimeNLM1->Fill(e1+e2,t12diff);
    }
  }
  else if(nMax==2)
  {
    fhM02EtaNLocMax2 [0][matched]->Fill(en,l0);
    fhMassEtaNLocMax2[0][matched]->Fill(en,mass);
    fhAsyEtaNLocMax2 [0][matched]->Fill(en,asym);
    if(fFillNCellHisto) fhNCellEtaNLocMax2[0][matched]->Fill(en,nc);
    
    if(!matched)
    {
      if(fFillHighMultHisto)
      {
        fhCentralityEtaNLocMax2->Fill(en,GetEventCentrality()) ;
        fhEventPlaneEtaNLocMax2->Fill(en,GetEventPlaneAngle()) ;
      }
      if(en > fHistoECut)fhEtaEtaPhiNLocMax2->Fill(eta,phi);
      fhEtaEPairDiffTimeNLM2->Fill(e1+e2,t12diff);
    }
  }
  else if(nMax >2)
  {
    fhM02EtaNLocMaxN [0][matched]->Fill(en,l0);
    fhMassEtaNLocMaxN[0][matched]->Fill(en,mass);
    fhAsyEtaNLocMaxN [0][matched]->Fill(en,asym);
    if(fFillNCellHisto) fhNCellEtaNLocMaxN[0][matched]->Fill(en,nc);
    
    if(!matched)
    {
      if(fFillHighMultHisto)
      {
        fhCentralityEtaNLocMaxN->Fill(en,GetEventCentrality()) ;
        fhEventPlaneEtaNLocMaxN->Fill(en,GetEventPlaneAngle()) ;
      }
      if(en > fHistoECut)fhEtaEtaPhiNLocMaxN->Fill(eta,phi);
      fhEtaEPairDiffTimeNLMN->Fill(e1+e2,t12diff);
    }
  }
  
  if(IsDataMC() && mcindex > 0 && mcindex < 7)
  {
    if     (nMax==1)
    {
      fhM02EtaNLocMax1[mcindex][matched]->Fill(en,l0);
      fhMassEtaNLocMax1[mcindex][matched]->Fill(en,mass);
      fhAsyEtaNLocMax1[mcindex][matched]->Fill(en,asym);
      if(fFillNCellHisto) fhNCellEtaNLocMax1[mcindex][matched]->Fill(en,nc);
    }
    else if(nMax==2)
    {
      fhM02EtaNLocMax2 [mcindex][matched]->Fill(en,l0);
      fhMassEtaNLocMax2[mcindex][matched]->Fill(en,mass);
      fhAsyEtaNLocMax2 [mcindex][matched]->Fill(en,asym);
      if(fFillNCellHisto) fhNCellEtaNLocMax2[mcindex][matched]->Fill(en,nc);
      
    }
    else if(nMax >2)
    {
      fhM02EtaNLocMaxN[mcindex][matched]->Fill(en,l0);
      fhMassEtaNLocMaxN[mcindex][matched]->Fill(en,mass);
      fhAsyEtaNLocMaxN[mcindex][matched]->Fill(en,asym);
      if(fFillNCellHisto) fhNCellEtaNLocMaxN[mcindex][matched]->Fill(en,nc);
    }
  }//Work with MC truth
}


//__________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillIdConvHistograms(Float_t en,    Int_t nMax, Float_t asym,
                                                            Float_t mass,   Float_t l0,
                                                            Bool_t matched, Int_t mcindex)
{
  // Fill histograms for clusters passing the photon selection
  
  if     (nMax==1)
  {
    fhM02ConNLocMax1 [0][matched]->Fill(en,l0);
    fhMassConNLocMax1[0][matched]->Fill(en,mass);
    fhAsyConNLocMax1 [0][matched]->Fill(en,asym);
  }
  else if(nMax==2)
  {
    fhM02ConNLocMax2 [0][matched]->Fill(en,l0);
    fhMassConNLocMax2[0][matched]->Fill(en,mass);
    fhAsyConNLocMax2 [0][matched]->Fill(en,asym);
  }
  else if(nMax >2)
  {
    fhM02ConNLocMaxN [0][matched]->Fill(en,l0);
    fhMassConNLocMaxN[0][matched]->Fill(en,mass);
    fhAsyConNLocMaxN [0][matched]->Fill(en,asym);
  }
  
  if(IsDataMC() && mcindex > 0 && mcindex < 7)
  {
    if     (nMax==1)
    {
      fhM02ConNLocMax1 [mcindex][matched]->Fill(en,l0);
      fhMassConNLocMax1[mcindex][matched]->Fill(en,mass);
      fhAsyConNLocMax1 [mcindex][matched]->Fill(en,asym);
    }
    else if(nMax==2)
    {
      fhM02ConNLocMax2 [mcindex][matched]->Fill(en,l0);
      fhMassConNLocMax2[mcindex][matched]->Fill(en,mass);
      fhAsyConNLocMax2 [mcindex][matched]->Fill(en,asym);
    }
    else if(nMax >2)
    {
      fhM02ConNLocMaxN [mcindex][matched]->Fill(en,l0);
      fhMassConNLocMaxN[mcindex][matched]->Fill(en,mass);
      fhAsyConNLocMaxN [mcindex][matched]->Fill(en,asym);
    }
    
  }//Work with MC truth
}

//_______________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillMCHistograms(Float_t en,        Float_t e1  , Float_t e2,
                                                        Int_t ebin,        Int_t mcindex,Int_t noverlaps,
                                                        Float_t l0,        Float_t mass,
                                                        Int_t nMax,        Bool_t  matched,
                                                        Float_t splitFrac, Float_t asym,
                                                        Float_t eprim,     Float_t asymGen)
{
  // Fill histograms needing some MC input
    
  Float_t efrac      = eprim/en;
  Float_t efracSplit = 0;
  if(e1+e2 > 0) efracSplit = eprim/(e1+e2);
  Float_t asymDiff = TMath::Abs(asym) - TMath::Abs(asymGen);
  
  //printf("e1 %2.2f, e2 %2.2f, eprim %2.2f, ereco %2.2f, esplit/ereco %2.2f, egen/ereco %2.2f, egen/esplit %2.2f\n",
  //       e1,e2,eprim,en,splitFrac,efrac,efracSplit);
  
  if(ebin >= 0 && fFillEbinHisto)
  {
    if( !matched ) fhMCGenFracNLocMaxEbin       [mcindex][ebin]->Fill(efrac,nMax);
    else           fhMCGenFracNLocMaxEbinMatched[mcindex][ebin]->Fill(efrac,nMax);
  }

  if     (nMax==1)
  {
    fhMCGenFracNLocMax1      [mcindex][matched]->Fill(en     ,  efrac );
    fhMCGenSplitEFracNLocMax1[mcindex][matched]->Fill(en     ,  efracSplit );
    fhMCGenEvsSplitENLocMax1 [mcindex][matched]->Fill(eprim  ,  e1+e2);
    if(asym > 0 && !matched)
    {
      if      (mcindex==kmcPi0)    fhAsyMCGenRecoDiffMCPi0[0]    ->Fill(en, asymDiff );
      else  if(mcindex==kmcPi0Conv)fhAsyMCGenRecoDiffMCPi0Conv[0]->Fill(en, asymDiff );
    }

    if(noverlaps==0)
    {
      fhMCGenFracNLocMax1NoOverlap      [mcindex][matched]->Fill(en ,  efrac );
      fhMCGenSplitEFracNLocMax1NoOverlap[mcindex][matched]->Fill(en ,  efracSplit );
    }
    
    if( en > fHistoECut )
    {
      fhMCGenEFracvsSplitEFracNLocMax1[mcindex][matched]->Fill(efrac,splitFrac );
      
      if(!matched && ebin >= 0 && fFillEbinHisto)
      {
        fhM02MCGenFracNLocMax1Ebin [mcindex][ebin]->Fill(efrac  ,  l0    );
        fhMassMCGenFracNLocMax1Ebin[mcindex][ebin]->Fill(efrac  ,  mass  );
        
        if(mcindex==kmcPi0 || mcindex==kmcPi0Conv)
        {
          fhMCAsymM02NLocMax1MCPi0Ebin [ebin]->Fill(l0  ,  asymGen );
          fhAsyMCGenRecoNLocMax1EbinPi0[ebin]->Fill(asym,  asymGen );
        }
      }
    }
  }
  else if(nMax==2)
  {
    fhMCGenFracNLocMax2      [mcindex][matched]->Fill(en     ,  efrac );
    fhMCGenSplitEFracNLocMax2[mcindex][matched]->Fill(en     ,  efracSplit );
    fhMCGenEvsSplitENLocMax2 [mcindex][matched]->Fill(eprim  ,  e1+e2);

    if(asym > 0 && !matched)
    {
     if      (mcindex==kmcPi0)    fhAsyMCGenRecoDiffMCPi0[1]    ->Fill(en, asymDiff );
     else  if(mcindex==kmcPi0Conv)fhAsyMCGenRecoDiffMCPi0Conv[1]->Fill(en, asymDiff );
    }
    
    if(noverlaps==0)
    {
      fhMCGenFracNLocMax2NoOverlap      [mcindex][matched]->Fill(en ,  efrac );
      fhMCGenSplitEFracNLocMax2NoOverlap[mcindex][matched]->Fill(en ,  efracSplit );
    }
    
    if( en > fHistoECut )
    {
      fhMCGenEFracvsSplitEFracNLocMax2[mcindex][matched]->Fill(efrac,splitFrac );
      
      if(!matched && ebin >= 0 && fFillEbinHisto)
      {
        fhM02MCGenFracNLocMax2Ebin [mcindex][ebin]->Fill(efrac  ,  l0    );
        fhMassMCGenFracNLocMax2Ebin[mcindex][ebin]->Fill(efrac  ,  mass  );
        if(mcindex==kmcPi0 || mcindex==kmcPi0Conv)
        {
          fhMCAsymM02NLocMax2MCPi0Ebin [ebin]->Fill(l0  ,  asymGen );
          fhAsyMCGenRecoNLocMax2EbinPi0[ebin]->Fill(asym,  asymGen );
        }
      }
    }

  }
  else if(nMax > 2 )
  {
    fhMCGenFracNLocMaxN      [mcindex][matched]->Fill(en     ,  efrac );
    fhMCGenSplitEFracNLocMaxN[mcindex][matched]->Fill(en     ,  efracSplit );
    fhMCGenEvsSplitENLocMaxN [mcindex][matched]->Fill(eprim  ,  e1+e2);
    if(asym > 0 && !matched)
    {
      if      (mcindex==kmcPi0)    fhAsyMCGenRecoDiffMCPi0[2]    ->Fill(en, asymDiff );
      else  if(mcindex==kmcPi0Conv)fhAsyMCGenRecoDiffMCPi0Conv[2]->Fill(en, asymDiff );
    }

    if(noverlaps==0)
    {
      fhMCGenFracNLocMaxNNoOverlap      [mcindex][matched]->Fill(en ,  efrac );
      fhMCGenSplitEFracNLocMaxNNoOverlap[mcindex][matched]->Fill(en ,  efracSplit );
    }
    
    if( en > fHistoECut )
    {
      fhMCGenEFracvsSplitEFracNLocMaxN[mcindex][matched]->Fill(efrac,splitFrac );
      
      if(!matched && ebin >= 0 && fFillEbinHisto)
      {
        fhM02MCGenFracNLocMaxNEbin [mcindex][ebin]->Fill(efrac  ,  l0    );
        fhMassMCGenFracNLocMaxNEbin[mcindex][ebin]->Fill(efrac  ,  mass  );
        
        if(mcindex==kmcPi0 || mcindex==kmcPi0Conv)
        {
          fhMCAsymM02NLocMaxNMCPi0Ebin [ebin]->Fill(l0  ,  asymGen );
          fhAsyMCGenRecoNLocMaxNEbinPi0[ebin]->Fill(asym,  asymGen );
        }
      }
    }
  }
}

//__________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillMCOverlapHistograms(Float_t en,      Float_t enprim,
                                                               Int_t   nc,      Float_t mass,    Float_t l0,
                                                               Float_t asym,    Float_t splitFrac,
                                                               Int_t   inlm,    Int_t ebin, Bool_t matched,
                                                               Int_t   mcindex, Int_t noverlaps)
{
  // Fill histograms for MC Overlaps
  
  //printf("en %f,mass %f,l0 %f,inlm %d,ebin %d,matched %d,mcindex %d,noverlaps %d \n",en,mass,l0,inlm,ebin,matched,mcindex,noverlaps);
    
  //printf("AliAnaInsideClusterInvariantMass::FillMCOverlapHistograms - NLM bin=%d, mcIndex %d, n Overlaps %d\n",inlm,mcindex,noverlaps);
  
  if(!matched)
  {
    fhMCENOverlaps[inlm][mcindex]->Fill(en,noverlaps);
    
    if     (noverlaps == 0)
    {
      fhMCEM02Overlap0  [inlm][mcindex]->Fill(en, l0);
      fhMCEMassOverlap0 [inlm][mcindex]->Fill(en, mass);
      fhMCEEpriOverlap0 [inlm][mcindex]->Fill(en, enprim);
      fhMCEAsymOverlap0 [inlm][mcindex]->Fill(en, TMath::Abs(asym));
      if(fFillNCellHisto) fhMCENCellOverlap0[inlm][mcindex]->Fill(en, nc);
      fhMCESplitEFracOverlap0[inlm][mcindex]->Fill(en, splitFrac);
      if((mcindex==kmcPi0 || mcindex == kmcPi0Conv) && ebin >=0) fhMCPi0MassM02Overlap0[inlm][ebin]->Fill(l0,mass);
    }
    else if(noverlaps == 1)
    {
      fhMCEM02Overlap1  [inlm][mcindex]->Fill(en, l0);
      fhMCEMassOverlap1 [inlm][mcindex]->Fill(en, mass);
      fhMCEEpriOverlap1 [inlm][mcindex]->Fill(en, enprim);
      fhMCEAsymOverlap1 [inlm][mcindex]->Fill(en, TMath::Abs(asym));
      if(fFillNCellHisto) fhMCENCellOverlap1[inlm][mcindex]->Fill(en, nc);
      fhMCESplitEFracOverlap1[inlm][mcindex]->Fill(en, splitFrac);
      if((mcindex==kmcPi0 || mcindex == kmcPi0Conv) && ebin >=0) fhMCPi0MassM02Overlap1[inlm][ebin]->Fill(l0,mass);
    }
    else if(noverlaps  > 1)
    {
      fhMCEM02OverlapN  [inlm][mcindex]->Fill(en, l0);
      fhMCEMassOverlapN [inlm][mcindex]->Fill(en, mass);
      fhMCEEpriOverlapN [inlm][mcindex]->Fill(en, enprim);
      fhMCEAsymOverlapN [inlm][mcindex]->Fill(en, TMath::Abs(asym));
      if(fFillNCellHisto) fhMCENCellOverlapN[inlm][mcindex]->Fill(en, nc);
      fhMCESplitEFracOverlapN[inlm][mcindex]->Fill(en, splitFrac);
      if((mcindex==kmcPi0 || mcindex == kmcPi0Conv) && ebin >=0) fhMCPi0MassM02OverlapN[inlm][ebin]->Fill(l0,mass);
    }
    else
      Info("FillMCOverlapHistograms","n overlaps = %d!!", noverlaps);
  }
  else if(fFillTMHisto)
  {
    fhMCENOverlapsMatch[inlm][mcindex]->Fill(en,noverlaps);
    
    if     (noverlaps == 0)
    {
      fhMCEM02Overlap0Match  [inlm][mcindex]->Fill(en, l0);
      fhMCEMassOverlap0Match [inlm][mcindex]->Fill(en, mass);
      fhMCEEpriOverlap0Match [inlm][mcindex]->Fill(en, enprim);
      fhMCEAsymOverlap0Match [inlm][mcindex]->Fill(en, TMath::Abs(asym));
      if(fFillNCellHisto) fhMCENCellOverlap0Match[inlm][mcindex]->Fill(en, nc);
      fhMCESplitEFracOverlap0Match[inlm][mcindex]->Fill(en, splitFrac);
      if((mcindex==kmcPi0 || mcindex == kmcPi0Conv) && ebin >=0) fhMCPi0MassM02Overlap0Match[inlm][ebin]->Fill(l0,mass);
    }
    else if(noverlaps == 1)
    {
      fhMCEM02Overlap1Match  [inlm][mcindex]->Fill(en, l0);
      fhMCEMassOverlap1Match [inlm][mcindex]->Fill(en, mass);
      fhMCEEpriOverlap1Match [inlm][mcindex]->Fill(en, enprim);
      fhMCEAsymOverlap1Match [inlm][mcindex]->Fill(en, TMath::Abs(asym));
      if(fFillNCellHisto) fhMCENCellOverlap1Match[inlm][mcindex]->Fill(en, nc);
      fhMCESplitEFracOverlap1Match[inlm][mcindex]->Fill(en, splitFrac);
      if((mcindex==kmcPi0 || mcindex == kmcPi0Conv) && ebin >=0) fhMCPi0MassM02Overlap1Match[inlm][ebin]->Fill(l0,mass);
    }
    else if(noverlaps  > 1)
    {
      fhMCEM02OverlapNMatch  [inlm][mcindex]->Fill(en, l0);
      fhMCEMassOverlapNMatch [inlm][mcindex]->Fill(en, mass);
      fhMCEEpriOverlapNMatch [inlm][mcindex]->Fill(en, enprim);
      fhMCEAsymOverlapNMatch [inlm][mcindex]->Fill(en, TMath::Abs(asym));
      if(fFillNCellHisto) fhMCENCellOverlapNMatch[inlm][mcindex]->Fill(en, nc);
      fhMCESplitEFracOverlapN[inlm][mcindex]->Fill(en, splitFrac);
      if((mcindex==kmcPi0 || mcindex == kmcPi0Conv) && ebin >=0) fhMCPi0MassM02OverlapNMatch[inlm][ebin]->Fill(l0,mass);
    }
    else
        Info("FillMCOverlapHistograms()","n overlaps in matched = %d!!", noverlaps);
  }
}


//_____________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillNCellHistograms(Int_t   ncells,  Float_t energy, Int_t nMax,
                                                           Bool_t  matched, Int_t mcindex,
                                                           Float_t mass   , Float_t l0)

{
  // Fill optional histograms with more SS parameters
    
  if     (nMax==1)
  {
    fhNCellNLocMax1[0][matched]->Fill(energy,ncells) ;
    if(mcindex > 0 )  fhNCellNLocMax1[mcindex][matched]->Fill(energy,ncells) ;
    
    if (mcindex==kmcPi0 && !matched)
    {
      if( energy > fHistoECut)
      {
        fhNCellMassEHighNLocMax1MCPi0->Fill(ncells,mass);
        fhNCellM02EHighNLocMax1MCPi0 ->Fill(ncells,l0);
      }
      else
      {
        fhNCellMassELowNLocMax1MCPi0->Fill(ncells,mass);
        fhNCellM02ELowNLocMax1MCPi0 ->Fill(ncells,l0);
      }
    }
  }
  else if( nMax == 2  )
  {
    fhNCellNLocMax2[0][matched]->Fill(energy,ncells) ;
    if(mcindex > 0 )  fhNCellNLocMax2[mcindex][matched]->Fill(energy,ncells) ;
    
    
    if (mcindex==kmcPi0 && !matched)
    {
      if( energy > fHistoECut)
      {
        fhNCellMassEHighNLocMax2MCPi0->Fill(ncells,mass);
        fhNCellM02EHighNLocMax2MCPi0 ->Fill(ncells,l0);
      }
      else
      {
        fhNCellMassELowNLocMax2MCPi0->Fill(ncells,mass);
        fhNCellM02ELowNLocMax2MCPi0 ->Fill(ncells,l0);
      }
    }
  }
  else if( nMax >= 3  )
  {
    fhNCellNLocMaxN[0][matched]->Fill(energy,ncells) ;
    if(mcindex > 0 )  fhNCellNLocMaxN[mcindex][matched]->Fill(energy,ncells) ;
    
    if (mcindex==kmcPi0 && !matched)
    {
      if( energy > fHistoECut)
      {
        fhNCellMassEHighNLocMaxNMCPi0->Fill(ncells,mass);
        fhNCellM02EHighNLocMaxNMCPi0 ->Fill(ncells,l0);
      }
      else
      {
        fhNCellMassELowNLocMaxNMCPi0->Fill(ncells,mass);
        fhNCellM02ELowNLocMaxNMCPi0 ->Fill(ncells,l0);
      }
    }
  }
}

//_____________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillSSExtraHistograms(AliVCluster  *cluster, Int_t nMax,
                                                             Bool_t  matched, Int_t mcindex,
                                                             Float_t mass   , Int_t ebin)
{
  // Fill optional histograms with more SS parameters
    
  Float_t en = cluster->E();
  
  // Get more Shower Shape parameters
  Float_t ll0  = 0., ll1  = 0.;
  Float_t disp= 0., dispEta = 0., dispPhi    = 0.;
  Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;
  
  GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), GetReader()->GetInputEvent()->GetEMCALCells(), cluster,
                                                                               ll0, ll1, disp, dispEta, dispPhi, sEta, sPhi, sEtaPhi);
  
  Float_t dispAsy = -1;
  if(dispEta+dispPhi >0 ) dispAsy = (dispPhi-dispEta) / (dispPhi+dispEta);
  
  if     (nMax==1)
  {
    if( en > fHistoECut )
    {
      fhMassDispEtaNLocMax1[0][matched]->Fill(dispEta,  mass );
      fhMassDispPhiNLocMax1[0][matched]->Fill(dispPhi,  mass );
      fhMassDispAsyNLocMax1[0][matched]->Fill(dispAsy,  mass );
      
      if(IsDataMC() && mcindex > 0 && mcindex < 7)
      {
        fhMassDispEtaNLocMax1[mcindex][matched]->Fill(dispEta,  mass );
        fhMassDispPhiNLocMax1[mcindex][matched]->Fill(dispPhi,  mass );
        fhMassDispAsyNLocMax1[mcindex][matched]->Fill(dispAsy,  mass );
      }
    }
    
    if(!matched && ebin >= 0 && fFillEbinHisto)
    {
      fhMassDispEtaNLocMax1Ebin[ebin]->Fill(dispEta,  mass );
      fhMassDispPhiNLocMax1Ebin[ebin]->Fill(dispPhi,  mass );
      fhMassDispAsyNLocMax1Ebin[ebin]->Fill(dispAsy,  mass );
    }
  }
  else if( nMax == 2  )
  {
    if( en > fHistoECut )
    {
      fhMassDispEtaNLocMax2[0][matched]->Fill(dispEta,  mass );
      fhMassDispPhiNLocMax2[0][matched]->Fill(dispPhi,  mass );
      fhMassDispAsyNLocMax2[0][matched]->Fill(dispAsy,  mass );
      
      if(IsDataMC() && mcindex > 0 && mcindex < 7)
      {
        fhMassDispEtaNLocMax2[mcindex][matched]->Fill(dispEta,  mass );
        fhMassDispPhiNLocMax2[mcindex][matched]->Fill(dispPhi,  mass );
        fhMassDispAsyNLocMax2[mcindex][matched]->Fill(dispAsy,  mass );
      }
    }
    
    if(!matched && ebin >= 0 && fFillEbinHisto)
    {
      fhMassDispEtaNLocMax2Ebin[ebin]->Fill(dispEta,  mass );
      fhMassDispPhiNLocMax2Ebin[ebin]->Fill(dispPhi,  mass );
      fhMassDispAsyNLocMax2Ebin[ebin]->Fill(dispAsy,  mass );
    }
    
  }
  else if( nMax >= 3  )
  {
    if( en > fHistoECut )
    {
      fhMassDispEtaNLocMaxN[0][matched]->Fill(dispEta,  mass );
      fhMassDispPhiNLocMaxN[0][matched]->Fill(dispPhi,  mass );
      fhMassDispAsyNLocMaxN[0][matched]->Fill(dispAsy,  mass );
      
      if(IsDataMC() && mcindex > 0 && mcindex < 7)
      {
        fhMassDispEtaNLocMaxN[mcindex][matched]->Fill(dispEta,  mass );
        fhMassDispPhiNLocMaxN[mcindex][matched]->Fill(dispPhi,  mass );
        fhMassDispAsyNLocMaxN[mcindex][matched]->Fill(dispAsy,  mass );
      }
    }
    
    if(!matched && ebin >= 0 && fFillEbinHisto)
    {
      fhMassDispEtaNLocMaxNEbin[ebin]->Fill(dispEta,  mass );
      fhMassDispPhiNLocMaxNEbin[ebin]->Fill(dispPhi,  mass );
      fhMassDispAsyNLocMaxNEbin[ebin]->Fill(dispAsy,  mass );
    }

  }
  
}

//__________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillSSWeightHistograms(AliVCluster *clus,  Int_t nlm,
                                                              Int_t absId1, Int_t absId2)
{
  // Calculate weights and fill histograms
    
  AliVCaloCells* cells = 0;
  if(fCalorimeter == "EMCAL") cells = GetEMCALCells();
  else                        cells = GetPHOSCells();
  
  // First recalculate energy in case non linearity was applied
  Float_t  energy = 0;
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++)
  {
    
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp,fCalorimeter, id);
    
    energy    += amp;
      
  } // energy loop
  
  if(energy <=0 )
  {
    Info("WeightHistograms()","Wrong calculated energy %f\n",energy);
    return;
  }
  
  //Get amplitude of  main local maxima, recalibrate if needed
  Float_t amp1 = cells->GetCellAmplitude(absId1);
  GetCaloUtils()->RecalibrateCellAmplitude(amp1,fCalorimeter, absId1);
  Float_t amp2 = cells->GetCellAmplitude(absId2);
  GetCaloUtils()->RecalibrateCellAmplitude(amp2,fCalorimeter, absId2);

  if(amp1 < amp2)        Info("FillSSWeightHistograms","Bad local maxima E ordering : id1 E %f, id2 E %f\n ",amp1,amp2);
  if(amp1==0 || amp2==0) Info("FillSSWeightHistograms","Null E local maxima : id1 E %f, id2 E %f\n "        ,amp1,amp2);
  
  if(amp1>0)fhPi0CellEMaxEMax2Frac   [nlm]->Fill(energy,amp2/amp1);
  fhPi0CellEMaxClusterFrac [nlm]->Fill(energy,amp1/energy);
  fhPi0CellEMax2ClusterFrac[nlm]->Fill(energy,amp2/energy);
  
  //Get the ratio and log ratio to all cells in cluster
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++)
  {
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp,fCalorimeter, id);
    
    if(amp > 0)fhPi0CellE       [nlm]->Fill(energy,amp);
    fhPi0CellEFrac   [nlm]->Fill(energy,amp/energy);
    fhPi0CellLogEFrac[nlm]->Fill(energy,TMath::Log(amp/energy));
    
    if     (id!=absId1 && id!=absId2)
    {
      if(amp1>0)fhPi0CellEMaxFrac [nlm]->Fill(energy,amp/amp1);
      if(amp2>0)fhPi0CellEMax2Frac[nlm]->Fill(energy,amp/amp2);
    }

  }
  
  //Recalculate shower shape for different W0
  if(fCalorimeter=="EMCAL")
  {
    Float_t l0org = clus->GetM02();
    Float_t l1org = clus->GetM20();
    Float_t dorg  = clus->GetDispersion();
    Float_t w0org =  GetCaloUtils()->GetEMCALRecoUtils()->GetW0();
    
    for(Int_t iw = 0; iw < fSSWeightN; iw++)
    {
      GetCaloUtils()->GetEMCALRecoUtils()->SetW0(fSSWeight[iw]);
      //GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), cells, clus);
      
      Float_t l0   = 0., l1   = 0.;
      Float_t disp = 0., dEta = 0., dPhi    = 0.;
      Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;
      
      RecalculateClusterShowerShapeParametersWithCellCut(GetEMCALGeometry(), cells, clus,l0,l1,disp,
                                                         dEta, dPhi, sEta, sPhi, sEtaPhi,0);

      
      fhM02WeightPi0[nlm][iw]->Fill(energy,clus->GetM02());
      
    } // w0 loop
    
    // Set the original values back
    clus->SetM02(l0org);
    clus->SetM20(l1org);
    clus->SetDispersion(dorg);
    GetCaloUtils()->GetEMCALRecoUtils()->SetW0(w0org);

    for(Int_t iec = 0; iec < fSSECellCutN; iec++)
    {
      Float_t l0   = 0., l1   = 0.;
      Float_t disp = 0., dEta = 0., dPhi    = 0.;
      Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;
      
      RecalculateClusterShowerShapeParametersWithCellCut(GetEMCALGeometry(), cells, clus,l0,l1,disp,
                                                         dEta, dPhi, sEta, sPhi, sEtaPhi,fSSECellCut[iec]);
      
      //printf("E %f, l0 org %f, l0 new %f, slope %f\n",clus->E(),l0org,l0,fSSECellCut[iec]);
      fhM02ECellCutPi0[nlm][iec]->Fill(energy,l0);
      
    } // w0 loop
  
  }// EMCAL
}

//____________________________________________________________________________________________
void  AliAnaInsideClusterInvariantMass::FillTrackMatchingHistograms(AliVCluster * cluster,
                                                                    Int_t nMax, Int_t mcindex)
{
  // Fill histograms related to track matching
    
  Float_t dZ  = cluster->GetTrackDz();
  Float_t dR  = cluster->GetTrackDx();
  Float_t en  = cluster->E();
  
  if(cluster->IsEMCAL() && GetCaloUtils()->IsRecalculationOfClusterTrackMatchingOn())
  {
    dR = 2000., dZ = 2000.;
    GetCaloUtils()->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dZ,dR);
  }
  
  //printf("Pi0EbE: dPhi %f, dEta %f\n",dR,dZ);
  
  if(TMath::Abs(dR) < 999)
  {
    if     ( nMax == 1  ) { fhTrackMatchedDEtaNLocMax1[0]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax1[0]->Fill(en,dR); }
    else if( nMax == 2  ) { fhTrackMatchedDEtaNLocMax2[0]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax2[0]->Fill(en,dR); }
    else if( nMax >= 3  ) { fhTrackMatchedDEtaNLocMaxN[0]->Fill(en,dZ); fhTrackMatchedDPhiNLocMaxN[0]->Fill(en,dR); }
    
    if(IsDataMC() && mcindex > 0 && mcindex < 7)
    {
      if     ( nMax == 1  ) { fhTrackMatchedDEtaNLocMax1[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax1[mcindex]->Fill(en,dR); }
      else if( nMax == 2  ) { fhTrackMatchedDEtaNLocMax2[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax2[mcindex]->Fill(en,dR); }
      else if( nMax >= 3  ) { fhTrackMatchedDEtaNLocMaxN[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiNLocMaxN[mcindex]->Fill(en,dR); }
    }
    
    AliVTrack *track = GetCaloUtils()->GetMatchedTrack(cluster, GetReader()->GetInputEvent());
    
    Bool_t positive = kFALSE;
    if(track) positive = (track->Charge()>0);

    if(track)
    {
      if(positive)
      {
        if     ( nMax == 1  ) { fhTrackMatchedDEtaNLocMax1Pos[0]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax1Pos[0]->Fill(en,dR); }
        else if( nMax == 2  ) { fhTrackMatchedDEtaNLocMax2Pos[0]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax2Pos[0]->Fill(en,dR); }
        else if( nMax >= 3  ) { fhTrackMatchedDEtaNLocMaxNPos[0]->Fill(en,dZ); fhTrackMatchedDPhiNLocMaxNPos[0]->Fill(en,dR); }
        
        if(IsDataMC() && mcindex > 0 && mcindex < 7)
        {
          if     ( nMax == 1  ) { fhTrackMatchedDEtaNLocMax1Pos[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax1Pos[mcindex]->Fill(en,dR); }
          else if( nMax == 2  ) { fhTrackMatchedDEtaNLocMax2Pos[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax2Pos[mcindex]->Fill(en,dR); }
          else if( nMax >= 3  ) { fhTrackMatchedDEtaNLocMaxNPos[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiNLocMaxNPos[mcindex]->Fill(en,dR); }
        }
      }
      else
      {
        if     ( nMax == 1  ) { fhTrackMatchedDEtaNLocMax1Neg[0]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax1Neg[0]->Fill(en,dR); }
        else if( nMax == 2  ) { fhTrackMatchedDEtaNLocMax2Neg[0]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax2Neg[0]->Fill(en,dR); }
        else if( nMax >= 3  ) { fhTrackMatchedDEtaNLocMaxNNeg[0]->Fill(en,dZ); fhTrackMatchedDPhiNLocMaxNNeg[0]->Fill(en,dR); }
        
        if(IsDataMC() && mcindex > 0 && mcindex < 7)
        {
          if     ( nMax == 1  ) { fhTrackMatchedDEtaNLocMax1Neg[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax1Neg[mcindex]->Fill(en,dR); }
          else if( nMax == 2  ) { fhTrackMatchedDEtaNLocMax2Neg[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax2Neg[mcindex]->Fill(en,dR); }
          else if( nMax >= 3  ) { fhTrackMatchedDEtaNLocMaxNNeg[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiNLocMaxNNeg[mcindex]->Fill(en,dR); }
        }
      }
      
    }// track exists
    
  }
}

//_______________________________________________________________
TObjString *  AliAnaInsideClusterInvariantMass::GetAnalysisCuts()
{	
	//Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaInsideClusterInvariantMass ---\n") ;
  parList+=onePar ;	
  
  snprintf(onePar,buffersize,"Calorimeter: %s\n",        fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fNLocMaxCutE =%2.2f \n",    GetCaloUtils()->GetLocalMaximaCutE()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fNLocMaxCutEDiff =%2.2f \n",GetCaloUtils()->GetLocalMaximaCutEDiff()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinNCells =%d \n",        fMinNCells) ;
  parList+=onePar ;    
  snprintf(onePar,buffersize,"fMinBadDist =%1.1f \n",    fMinBadDist) ;
  parList+=onePar ;  
  if(fFillSSWeightHisto)
  {
    snprintf(onePar,buffersize," N w %d - N e cut %d \n",fSSWeightN,fSSECellCutN);
    parList+=onePar ;
  }
  
  return new TObjString(parList) ;
  
}

//________________________________________________________________
TList * AliAnaInsideClusterInvariantMass::GetCreateOutputObjects()
{
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("InsideClusterHistos") ;
  
  Int_t nptbins  = GetHistogramRanges()->GetHistoPtBins();           Float_t ptmax  = GetHistogramRanges()->GetHistoPtMax();           Float_t ptmin  = GetHistogramRanges()->GetHistoPtMin();
  Int_t ssbins   = GetHistogramRanges()->GetHistoShowerShapeBins();  Float_t ssmax  = GetHistogramRanges()->GetHistoShowerShapeMax();  Float_t ssmin  = GetHistogramRanges()->GetHistoShowerShapeMin();
  Int_t mbins    = GetHistogramRanges()->GetHistoMassBins();         Float_t mmax   = GetHistogramRanges()->GetHistoMassMax();         Float_t mmin   = GetHistogramRanges()->GetHistoMassMin();
  Int_t ncbins   = GetHistogramRanges()->GetHistoNClusterCellBins(); Int_t   ncmax  = GetHistogramRanges()->GetHistoNClusterCellMax(); Int_t   ncmin  = GetHistogramRanges()->GetHistoNClusterCellMin(); 
  Int_t nphibins = GetHistogramRanges()->GetHistoPhiBins();          Float_t phimax = GetHistogramRanges()->GetHistoPhiMax();          Float_t phimin = GetHistogramRanges()->GetHistoPhiMin();
  Int_t netabins = GetHistogramRanges()->GetHistoEtaBins();          Float_t etamax = GetHistogramRanges()->GetHistoEtaMax();          Float_t etamin = GetHistogramRanges()->GetHistoEtaMin();
  
  Int_t   nresetabins = GetHistogramRanges()->GetHistoTrackResidualEtaBins();          
  Float_t resetamax   = GetHistogramRanges()->GetHistoTrackResidualEtaMax();          
  Float_t resetamin   = GetHistogramRanges()->GetHistoTrackResidualEtaMin();
  Int_t   nresphibins = GetHistogramRanges()->GetHistoTrackResidualPhiBins();          
  Float_t resphimax   = GetHistogramRanges()->GetHistoTrackResidualPhiMax();          
  Float_t resphimin   = GetHistogramRanges()->GetHistoTrackResidualPhiMin();  
  
  Bool_t m02On   = GetCaloPID()->IsSplitShowerShapeCutOn();
  Bool_t asyOn   = GetCaloPID()->IsSplitAsymmetryCutOn();
  Bool_t splitOn = kFALSE;
  if(GetCaloPID()->GetSplitEnergyFractionMinimum(0) > 0 ||
     GetCaloPID()->GetSplitEnergyFractionMinimum(1) > 0 ||
     GetCaloPID()->GetSplitEnergyFractionMinimum(2) > 0) splitOn = kTRUE;
  
  TString ptype[] ={"","#gamma","#gamma->e^{#pm}","#pi^{0}","#pi^{0} (#gamma->e^{#pm})","#eta", "hadron"};
  TString pname[] ={"","Photon","Conversion",     "Pi0",    "Pi0Conv",                  "Eta","Hadron"};
  TString snlm [] = {"1","2","N"};

  TString sEBin[] = {"8 < E < 12 GeV","12 < E < 16 GeV", "16 < E < 20 GeV", "E > 20 GeV" };

  Int_t n = 1;
  
  if(IsDataMC()) n = 7;
  
  Int_t nMaxBins = 10;
  
  TString sMatched[] = {"","Matched"};
  
  Int_t nMatched = 2;
  if(!fFillTMHisto) nMatched = 1;
  
  if(fCheckSplitDistToBad)
  {
    for(Int_t inlm = 0; inlm < 3; inlm++)
    {
      fhMassBadDistClose[inlm]  = new TH2F(Form("hMassBadDistCloseNLocMax%s",snlm[inlm].Data()),
                                           Form("Invariant mass of splitted cluster with NLM=%d vs E, 2nd LM close to bad channel",inlm),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassBadDistClose[inlm]->SetYTitle("M (GeV/c^{2})");
      fhMassBadDistClose[inlm]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassBadDistClose[inlm]) ;
      
      fhM02BadDistClose[inlm]  = new TH2F(Form("hM02BadDistCloseNLocMax%s",snlm[inlm].Data()),
                                          Form("#lambda_{0}^{2} for cluster with NLM=%d vs E, 2nd LM close to bad channel",inlm),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhM02BadDistClose[inlm]->SetYTitle("#lambda_{0}^{2}");
      fhM02BadDistClose[inlm]->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02BadDistClose[inlm]) ;
      
      fhMassOnBorder[inlm]  = new TH2F(Form("hMassOnBorderNLocMax%s",snlm[inlm].Data()),
                                       Form("Invariant mass of splitted cluster with NLM=%d vs E, 2nd LM close to border",inlm),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassOnBorder[inlm]->SetYTitle("M (GeV/c^{2})");
      fhMassOnBorder[inlm]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassOnBorder[inlm]) ;
      
      fhM02OnBorder[inlm]  = new TH2F(Form("hM02OnBorderNLocMax%s",snlm[inlm].Data()),
                                      Form("#lambda_{0}^{2} for cluster with NLM=%d vs E, 2nd LM close to border",inlm),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhM02OnBorder[inlm]->SetYTitle("#lambda_{0}^{2}");
      fhM02OnBorder[inlm]->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02OnBorder[inlm]) ;
      
    }
  }
  
  for(Int_t i = 0; i < n; i++)
  {
    for(Int_t j = 0; j < nMatched; j++)
    {
      
      fhNLocMax[i][j]     = new TH2F(Form("hNLocMax%s%s",pname[i].Data(),sMatched[j].Data()),
                                     Form("Number of local maxima in cluster %s %s",ptype[i].Data(),sMatched[j].Data()),
                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
      fhNLocMax[i][j]   ->SetYTitle("N maxima");
      fhNLocMax[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhNLocMax[i][j]) ;

      fhLM1NLocMax[i][j]     = new TH2F(Form("hLM1NLocMax%s%s",pname[i].Data(),sMatched[j].Data()),
                                     Form("Number of local maxima in cluster for split cluster 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
      fhLM1NLocMax[i][j]   ->SetYTitle("N maxima");
      fhLM1NLocMax[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhLM1NLocMax[i][j]) ;

      fhLM2NLocMax[i][j]     = new TH2F(Form("hLM2NLocMax%s%s",pname[i].Data(),sMatched[j].Data()),
                                        Form("Number of local maxima in cluster for split cluster 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
      fhLM2NLocMax[i][j]   ->SetYTitle("N maxima");
      fhLM2NLocMax[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhLM2NLocMax[i][j]) ;
      
      if(m02On)
      {
        fhNLocMaxM02Cut[i][j] = new TH2F(Form("hNLocMaxM02Cut%s%s",pname[i].Data(),sMatched[j].Data()),
                                         Form("Number of local maxima in cluster %s %s, M02 cut",ptype[i].Data(),sMatched[j].Data()),
                                         nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
        fhNLocMaxM02Cut[i][j]->SetYTitle("N maxima");
        fhNLocMaxM02Cut[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhNLocMaxM02Cut[i][j]) ;
        
        fhLM1NLocMaxM02Cut[i][j]     = new TH2F(Form("hLM1NLocMaxM02Cut%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Number of local maxima in cluster for split cluster 1 %s %s, M02 cut",ptype[i].Data(),sMatched[j].Data()),
                                          nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
        fhLM1NLocMaxM02Cut[i][j]   ->SetYTitle("N maxima");
        fhLM1NLocMaxM02Cut[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhLM1NLocMaxM02Cut[i][j]) ;
        
        fhLM2NLocMaxM02Cut[i][j]     = new TH2F(Form("hLM2NLocMaxM02Cut%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Number of local maxima in cluster for split cluster 2 %s %s, M02 cut",ptype[i].Data(),sMatched[j].Data()),
                                          nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
        fhLM2NLocMaxM02Cut[i][j]   ->SetYTitle("N maxima");
        fhLM2NLocMaxM02Cut[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhLM2NLocMaxM02Cut[i][j]) ;

      }
      

      fhNLocMaxIdPi0[i][j]     = new TH2F(Form("hNLocMaxIdPi0%s%s",pname[i].Data(),sMatched[j].Data()),
                                     Form("Number of local maxima in pi0 ID cluster %s %s",ptype[i].Data(),sMatched[j].Data()),
                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
      fhNLocMaxIdPi0[i][j]   ->SetYTitle("N maxima");
      fhNLocMaxIdPi0[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhNLocMaxIdPi0[i][j]) ;

      
      fhLM1NLocMaxIdPi0[i][j]     = new TH2F(Form("hLM1NLocMaxIdPi0%s%s",pname[i].Data(),sMatched[j].Data()),
                                        Form("Number of local maxima in cluster for split cluster 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
      fhLM1NLocMaxIdPi0[i][j]   ->SetYTitle("N maxima");
      fhLM1NLocMaxIdPi0[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhLM1NLocMaxIdPi0[i][j]) ;
      
      fhLM2NLocMaxIdPi0[i][j]     = new TH2F(Form("hLM2NLocMaxIdPi0%s%s",pname[i].Data(),sMatched[j].Data()),
                                        Form("Number of local maxima in cluster for split cluster 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
      fhLM2NLocMaxIdPi0[i][j]   ->SetYTitle("N maxima");
      fhLM2NLocMaxIdPi0[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhLM2NLocMaxIdPi0[i][j]) ;
      

      
      
      fhSplitClusterENLocMax[i][j]     = new TH2F(Form("hSplitEClusterNLocMax%s%s",pname[i].Data(),sMatched[j].Data()),
                                                  Form("Number of local maxima vs E of split clusters %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                  nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
      fhSplitClusterENLocMax[i][j]   ->SetYTitle("N maxima");
      fhSplitClusterENLocMax[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhSplitClusterENLocMax[i][j]) ;
      
      
      fhSplitClusterEPi0NLocMax[i][j]     = new TH2F(Form("hSplitEClusterPi0NLocMax%s%s",pname[i].Data(),sMatched[j].Data()),
                                                     Form("Number of local maxima vs E of split clusters, id as pi0, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
      fhSplitClusterEPi0NLocMax[i][j]   ->SetYTitle("N maxima");
      fhSplitClusterEPi0NLocMax[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhSplitClusterEPi0NLocMax[i][j]) ;

      if(fFillNCellHisto)
      {
        fhNCellNLocMax1[i][j]  = new TH2F(Form("hNCellNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("n cells vs E for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                          nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
        fhNCellNLocMax1[i][j] ->SetYTitle("N cells");
        fhNCellNLocMax1[i][j] ->SetXTitle("E (GeV)");
        outputContainer->Add(fhNCellNLocMax1[i][j]) ;
        
        fhNCellNLocMax2[i][j]     = new TH2F(Form("hNCellNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("n cells vs E for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                             nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
        fhNCellNLocMax2[i][j]   ->SetYTitle("N cells");
        fhNCellNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhNCellNLocMax2[i][j]) ;
        
        
        fhNCellNLocMaxN[i][j]     = new TH2F(Form("hNCellNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("n cells vs E for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                             nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
        fhNCellNLocMaxN[i][j]   ->SetYTitle("N cells");
        fhNCellNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhNCellNLocMaxN[i][j]) ;
      }

      fhMassNLocMax1[i][j]  = new TH2F(Form("hMassNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of splitted cluster with NLM=1 vs E, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassNLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassNLocMax1[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassNLocMax1[i][j]) ;   
      
      fhMassNLocMax2[i][j]  = new TH2F(Form("hMassNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of splitted cluster with NLM=2 vs E, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassNLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassNLocMax2[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassNLocMax2[i][j]) ;   
      
      fhMassNLocMaxN[i][j]  = new TH2F(Form("hMassNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of splitted cluster with NLM>2 vs E, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassNLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassNLocMaxN[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassNLocMaxN[i][j]) ;

      fhMassSplitENLocMax1[i][j]  = new TH2F(Form("hMassSplitENLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of splitted cluster with NLM=1 vs E1+E2, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassSplitENLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassSplitENLocMax1[i][j]->SetXTitle("E1+E2 (GeV)");
      outputContainer->Add(fhMassSplitENLocMax1[i][j]) ;
      
      fhMassSplitENLocMax2[i][j]  = new TH2F(Form("hMassSplitENLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of splitted cluster with NLM=2 vs E1+E2, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassSplitENLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassSplitENLocMax2[i][j]->SetXTitle("E1+E2 (GeV)");
      outputContainer->Add(fhMassSplitENLocMax2[i][j]) ;
      
      fhMassSplitENLocMaxN[i][j]  = new TH2F(Form("hMassSplitENLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of splitted cluster with NLM>2 vs E1+E2, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassSplitENLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassSplitENLocMaxN[i][j]->SetXTitle("E1+E2 (GeV)");
      outputContainer->Add(fhMassSplitENLocMaxN[i][j]) ;
      
      fhM02NLocMax1[i][j]     = new TH2F(Form("hM02NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                         Form("#lambda_{0}^{2} vs E for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhM02NLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02NLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02NLocMax1[i][j]) ;
      
      fhM02NLocMax2[i][j]     = new TH2F(Form("hM02NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                         Form("#lambda_{0}^{2} vs E for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhM02NLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02NLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02NLocMax2[i][j]) ;
      
      fhM02NLocMaxN[i][j]    = new TH2F(Form("hM02NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                        Form("#lambda_{0}^{2} vs E for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhM02NLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02NLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02NLocMaxN[i][j]) ;
      
      fhAsymNLocMax1[i][j]  = new TH2F(Form("hAsymNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Asymmetry of NLM=1  vs cluster Energy, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,200,-1,1);
      fhAsymNLocMax1[i][j]->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
      fhAsymNLocMax1[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsymNLocMax1[i][j]) ;
      
      fhAsymNLocMax2[i][j]  = new TH2F(Form("hAsymNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Asymmetry of NLM=2  vs cluster Energy, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,200,-1,1);
      fhAsymNLocMax2[i][j]->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
      fhAsymNLocMax2[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsymNLocMax2[i][j]) ;
      
      fhAsymNLocMaxN[i][j]  = new TH2F(Form("hAsymNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Asymmetry of NLM>2  vs cluster Energy, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,200,-1,1);
      fhAsymNLocMaxN[i][j]->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
      fhAsymNLocMaxN[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsymNLocMaxN[i][j]) ;
      
      fhSplitEFractionNLocMax1[i][j]     = new TH2F(Form("hSplitEFractionNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                    Form("(E1+E2)/E_{cluster} vs E_{cluster} for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                    nptbins,ptmin,ptmax,120,0,1.2);
      fhSplitEFractionNLocMax1[i][j]   ->SetXTitle("E_{cluster} (GeV)");
      fhSplitEFractionNLocMax1[i][j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
      outputContainer->Add(fhSplitEFractionNLocMax1[i][j]) ;
      
      fhSplitEFractionNLocMax2[i][j]     = new TH2F(Form("hSplitEFractionNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                    Form("(E1+E2)/E_{cluster} vs E_{cluster} for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                    nptbins,ptmin,ptmax,120,0,1.2);
      fhSplitEFractionNLocMax2[i][j]   ->SetXTitle("E_{cluster} (GeV)");
      fhSplitEFractionNLocMax2[i][j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
      outputContainer->Add(fhSplitEFractionNLocMax2[i][j]) ;
      
      fhSplitEFractionNLocMaxN[i][j]    = new TH2F(Form("hSplitEFractionNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                   Form("(E1+E2)/E_{cluster} vs E_{cluster} for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                   nptbins,ptmin,ptmax,120,0,1.2);
      fhSplitEFractionNLocMaxN[i][j]   ->SetXTitle("E_{cluster} (GeV)");
      fhSplitEFractionNLocMaxN[i][j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
      outputContainer->Add(fhSplitEFractionNLocMaxN[i][j]) ;
      
      if(i==0 && j==0 )
      {
        if(m02On)
        {
          fhMassM02CutNLocMax1  = new TH2F("hMassM02CutNLocMax1","Invariant mass of splitted cluster with NLM=1 vs E, M02 cut, no TM",
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax);
          fhMassM02CutNLocMax1->SetYTitle("M (GeV/c^{2})");
          fhMassM02CutNLocMax1->SetXTitle("E (GeV)");
          outputContainer->Add(fhMassM02CutNLocMax1) ;
          
          fhMassM02CutNLocMax2  = new TH2F("hMassM02CutNLocMax2","Invariant mass of splitted cluster with NLM=2 vs E, M02 cut, no TM",
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax);
          fhMassM02CutNLocMax2->SetYTitle("M (GeV/c^{2})");
          fhMassM02CutNLocMax2->SetXTitle("E (GeV)");
          outputContainer->Add(fhMassM02CutNLocMax2) ;
          
          fhMassM02CutNLocMaxN  = new TH2F("hMassM02CutNLocMaxN","Invariant mass of splitted cluster with NLM>2 vs E, M02 cut, no TM",
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax);
          fhMassM02CutNLocMaxN->SetYTitle("M (GeV/c^{2})");
          fhMassM02CutNLocMaxN->SetXTitle("E (GeV)");
          outputContainer->Add(fhMassM02CutNLocMaxN) ;
          
          fhAsymM02CutNLocMax1  = new TH2F("hAsymM02CutNLocMax1","Asymmetry of NLM=1  vs cluster Energy, M02Cut, no TM", nptbins,ptmin,ptmax,200,-1,1);
          fhAsymM02CutNLocMax1->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
          fhAsymM02CutNLocMax1->SetXTitle("E (GeV)");
          outputContainer->Add(fhAsymM02CutNLocMax1) ;
          
          fhAsymM02CutNLocMax2  = new TH2F("hAsymM02CutNLocMax2","Asymmetry of NLM=2  vs cluster Energy, M02Cut, no TM", nptbins,ptmin,ptmax,200,-1,1);
          fhAsymM02CutNLocMax2->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
          fhAsymM02CutNLocMax2->SetXTitle("E (GeV)");
          outputContainer->Add(fhAsymM02CutNLocMax2) ;
          
          fhAsymM02CutNLocMaxN  = new TH2F("hAsymM02CutNLocMaxN","Asymmetry of NLM>2  vs cluster Energy, M02Cut, no TM", nptbins,ptmin,ptmax,200,-1,1);
          fhAsymM02CutNLocMaxN->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
          fhAsymM02CutNLocMaxN->SetXTitle("E (GeV)");
          outputContainer->Add(fhAsymM02CutNLocMaxN) ;
          
          if(splitOn)
          {
            fhMassSplitECutNLocMax1  = new TH2F("hMassSplitECutNLocMax1","Invariant mass of splitted cluster with NLM=1 vs E, (E1+E2)/E cut, M02 cut, no TM",
                                                nptbins,ptmin,ptmax,mbins,mmin,mmax);
            fhMassSplitECutNLocMax1->SetYTitle("M (GeV/c^{2})");
            fhMassSplitECutNLocMax1->SetXTitle("E (GeV)");
            outputContainer->Add(fhMassSplitECutNLocMax1) ;
            
            fhMassSplitECutNLocMax2  = new TH2F("hMassSplitECutNLocMax2","Invariant mass of splitted cluster with NLM=2 vs E, (E1+E2)/E cut, M02 cut, no TM",
                                                nptbins,ptmin,ptmax,mbins,mmin,mmax);
            fhMassSplitECutNLocMax2->SetYTitle("M (GeV/c^{2})");
            fhMassSplitECutNLocMax2->SetXTitle("E (GeV)");
            outputContainer->Add(fhMassSplitECutNLocMax2) ;
            
            fhMassSplitECutNLocMaxN  = new TH2F("hMassSplitECutNLocMaxN","Invariant mass of splitted cluster with NLM>2 vs E, (E1+E2)/E cut, M02 cut, no TM",
                                                nptbins,ptmin,ptmax,mbins,mmin,mmax);
            fhMassSplitECutNLocMaxN->SetYTitle("M (GeV/c^{2})");
            fhMassSplitECutNLocMaxN->SetXTitle("E (GeV)");
            outputContainer->Add(fhMassSplitECutNLocMaxN) ;
          }
        }//m02on
        
        if(asyOn)
        {
          fhMassAsyCutNLocMax1  = new TH2F("hMassAsyCutNLocMax1","Invariant mass of splitted cluster with NLM=1 vs E, Asy cut, no TM",
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax);
          fhMassAsyCutNLocMax1->SetYTitle("M (GeV/c^{2})");
          fhMassAsyCutNLocMax1->SetXTitle("E (GeV)");
          outputContainer->Add(fhMassAsyCutNLocMax1) ;
          
          fhMassAsyCutNLocMax2  = new TH2F("hMassAsyCutNLocMax2","Invariant mass of splitted cluster with NLM=2 vs E, Asy cut, no TM",
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax);
          fhMassAsyCutNLocMax2->SetYTitle("M (GeV/c^{2})");
          fhMassAsyCutNLocMax2->SetXTitle("E (GeV)");
          outputContainer->Add(fhMassAsyCutNLocMax2) ;
          
          fhMassAsyCutNLocMaxN  = new TH2F("hMassAsyCutNLocMaxN","Invariant mass of splitted cluster with NLM>2 vs E, Asy cut, no TM",
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax);
          fhMassAsyCutNLocMaxN->SetYTitle("M (GeV/c^{2})");
          fhMassAsyCutNLocMaxN->SetXTitle("E (GeV)");
          outputContainer->Add(fhMassAsyCutNLocMaxN) ;
          
          fhM02AsyCutNLocMax1  = new TH2F("hM02AsyCutNLocMax1","#lambda_{0}^{2} of NLM=1  vs cluster Energy, AsyCut, no TM",
                                          nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
          fhM02AsyCutNLocMax1->SetYTitle("#lambda_{0}^{2}");
          fhM02AsyCutNLocMax1->SetXTitle("E (GeV)");
          outputContainer->Add(fhM02AsyCutNLocMax1) ;
          
          fhM02AsyCutNLocMax2  = new TH2F("hM02AsyCutNLocMax2","#lambda_{0}^{2} of NLM=2  vs cluster Energy, AsyCut, no TM",
                                          nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
          fhM02AsyCutNLocMax2->SetYTitle("#lambda_{0}^{2}");
          fhM02AsyCutNLocMax2->SetXTitle("E (GeV)");
          outputContainer->Add(fhM02AsyCutNLocMax2) ;
          
          fhM02AsyCutNLocMaxN  = new TH2F("hM02AsyCutNLocMaxN","#lambda_{0}^{2} of NLM>2  vs cluster Energy, AsyCut, no TM",
                                          nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
          fhM02AsyCutNLocMaxN->SetYTitle("#lambda_{0}^{2}");
          fhM02AsyCutNLocMaxN->SetXTitle("E (GeV)");
          outputContainer->Add(fhM02AsyCutNLocMaxN) ;
        }
        
        if(GetCaloPID()->GetSubClusterEnergyMinimum(0) > 0.1)
        {
          fhMassEnCutNLocMax1  = new TH2F("hMassEnCutNLocMax1",Form("Invariant mass of splitted cluster with NLM=1 vs E, E > %1.1f GeV, no TM",GetCaloPID()->GetSubClusterEnergyMinimum(0)),
                                          nptbins,ptmin,ptmax,mbins,mmin,mmax);
          fhMassEnCutNLocMax1->SetYTitle("M (GeV/c^{2})");
          fhMassEnCutNLocMax1->SetXTitle("E (GeV)");
          outputContainer->Add(fhMassEnCutNLocMax1) ;
          
          fhMassEnCutNLocMax2  = new TH2F("hMassEnCutNLocMax2",Form("Invariant mass of splitted cluster with NLM=2 vs E, E > %1.1f GeV, no TM",GetCaloPID()->GetSubClusterEnergyMinimum(1)),
                                          nptbins,ptmin,ptmax,mbins,mmin,mmax);
          fhMassEnCutNLocMax2->SetYTitle("M (GeV/c^{2})");
          fhMassEnCutNLocMax2->SetXTitle("E (GeV)");
          outputContainer->Add(fhMassEnCutNLocMax2) ;
          
          fhMassEnCutNLocMaxN  = new TH2F("hMassEnCutNLocMaxN",Form("Invariant mass of splitted cluster with NLM>2 vs E, E > %1.1f GeV, no TM",GetCaloPID()->GetSubClusterEnergyMinimum(2)),
                                          nptbins,ptmin,ptmax,mbins,mmin,mmax);
          fhMassEnCutNLocMaxN->SetYTitle("M (GeV/c^{2})");
          fhMassEnCutNLocMaxN->SetXTitle("E (GeV)");
          outputContainer->Add(fhMassEnCutNLocMaxN) ;
          
          fhM02EnCutNLocMax1  = new TH2F("hM02EnCutNLocMax1",Form("#lambda_{0}^{2} of NLM=1  vs cluster Energy, E > %1.1f GeV, no TM",GetCaloPID()->GetSubClusterEnergyMinimum(0)),
                                         nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
          fhM02EnCutNLocMax1->SetYTitle("#lambda_{0}^{2}");
          fhM02EnCutNLocMax1->SetXTitle("E (GeV)");
          outputContainer->Add(fhM02EnCutNLocMax1) ;
          
          fhM02EnCutNLocMax2  = new TH2F("hM02EnCutNLocMax2",Form("#lambda_{0}^{2} of NLM=2  vs cluster Energy, E > %1.1f GeV, no TM",GetCaloPID()->GetSubClusterEnergyMinimum(1)),
                                         nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
          fhM02EnCutNLocMax2->SetYTitle("#lambda_{0}^{2}");
          fhM02EnCutNLocMax2->SetXTitle("E (GeV)");
          outputContainer->Add(fhM02EnCutNLocMax2) ;
          
          fhM02EnCutNLocMaxN  = new TH2F("hM02EnCutNLocMaxN",Form("#lambda_{0}^{2} of NLM>2  vs cluster Energy, E > %1.1f GeV, no TM",GetCaloPID()->GetSubClusterEnergyMinimum(2)),
                                         nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
          fhM02EnCutNLocMaxN->SetYTitle("#lambda_{0}^{2}");
          fhM02EnCutNLocMaxN->SetXTitle("E (GeV)");
          outputContainer->Add(fhM02EnCutNLocMaxN) ;
          
          fhAsymEnCutNLocMax1  = new TH2F("hAsymEnCutNLocMax1",Form("Asymmetry of NLM=1  vs cluster Energy, E > %1.1f GeV, no TM",GetCaloPID()->GetSubClusterEnergyMinimum(0))
                                          , nptbins,ptmin,ptmax,200,-1,1);
          fhAsymEnCutNLocMax1->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
          fhAsymEnCutNLocMax1->SetXTitle("E (GeV)");
          outputContainer->Add(fhAsymEnCutNLocMax1) ;
          
          fhAsymEnCutNLocMax2  = new TH2F("hAsymEnCutNLocMax2",Form("Asymmetry of NLM=2  vs cluster Energy, E > %1.1f GeV, no TM",GetCaloPID()->GetSubClusterEnergyMinimum(1))
                                          , nptbins,ptmin,ptmax,200,-1,1);
          fhAsymEnCutNLocMax2->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
          fhAsymEnCutNLocMax2->SetXTitle("E (GeV)");
          outputContainer->Add(fhAsymEnCutNLocMax2) ;
          
          fhAsymEnCutNLocMaxN  = new TH2F("hAsymEnCutNLocMaxN",Form("Asymmetry of NLM>2  vs cluster Energy, E > %1.1f GeV, no TM",GetCaloPID()->GetSubClusterEnergyMinimum(2))
                                          , nptbins,ptmin,ptmax,200,-1,1);
          fhAsymEnCutNLocMaxN->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
          fhAsymEnCutNLocMaxN->SetXTitle("E (GeV)");
          outputContainer->Add(fhAsymEnCutNLocMaxN) ;

          fhSplitEFracEnCutNLocMax1  = new TH2F("hSplitEFracEnCutNLocMax1",Form("SplitEFracmetry of NLM=1  vs cluster Energy, E > %1.1f GeV, no TM",GetCaloPID()->GetSubClusterEnergyMinimum(0))
                                          , nptbins,ptmin,ptmax,120,0,1.2);
          fhSplitEFracEnCutNLocMax1->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
          fhSplitEFracEnCutNLocMax1->SetXTitle("E (GeV)");
          outputContainer->Add(fhSplitEFracEnCutNLocMax1) ;
          
          fhSplitEFracEnCutNLocMax2  = new TH2F("hSplitEFracEnCutNLocMax2",Form("SplitEFracmetry of NLM=2  vs cluster Energy, E > %1.1f GeV, no TM",GetCaloPID()->GetSubClusterEnergyMinimum(1))
                                          , nptbins,ptmin,ptmax,120,0,1.2);
          fhSplitEFracEnCutNLocMax2->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
          fhSplitEFracEnCutNLocMax2->SetXTitle("E (GeV)");
          outputContainer->Add(fhSplitEFracEnCutNLocMax2) ;
          
          fhSplitEFracEnCutNLocMaxN  = new TH2F("hSplitEFracEnCutNLocMaxN",Form("SplitEFracmetry of NLM>2  vs cluster Energy, E > %1.1f GeV, no TM",GetCaloPID()->GetSubClusterEnergyMinimum(2))
                                          , nptbins,ptmin,ptmax,120,0,1.2);
          fhSplitEFracEnCutNLocMaxN->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
          fhSplitEFracEnCutNLocMaxN->SetXTitle("E (GeV)");
          outputContainer->Add(fhSplitEFracEnCutNLocMaxN) ;
        }
        
      } // no MC
      
      if(asyOn || m02On )
      {
        fhMassAfterCutsNLocMax1[i][j]     = new TH2F(Form("hMassAfterCutsNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                     Form("Mass vs E, %s %s, for NLM = 1, M02 and asy cut",ptype[i].Data(),sMatched[j].Data()),
                                                     nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassAfterCutsNLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
        fhMassAfterCutsNLocMax1[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassAfterCutsNLocMax1[i][j]) ;
        
        fhMassAfterCutsNLocMax2[i][j]     = new TH2F(Form("hMassAfterCutsNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                     Form("Mass vs E, %s %s, for NLM = 2, M02 and asy cut",ptype[i].Data(),sMatched[j].Data()),
                                                     nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassAfterCutsNLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
        fhMassAfterCutsNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassAfterCutsNLocMax2[i][j]) ;
        
        fhMassAfterCutsNLocMaxN[i][j]     = new TH2F(Form("hMassAfterCutsNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                     Form("Mass vs E, %s %s, for NLM > 2, M02 and asy cut",ptype[i].Data(),sMatched[j].Data()),
                                                     nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassAfterCutsNLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
        fhMassAfterCutsNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassAfterCutsNLocMaxN[i][j]) ;

        fhMassSplitEAfterCutsNLocMax1[i][j]     = new TH2F(Form("hMassSplitEAfterCutsNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                     Form("Mass vs E1+E2, %s %s, for NLM = 1, M02 and asy cut",ptype[i].Data(),sMatched[j].Data()),
                                                     nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassSplitEAfterCutsNLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
        fhMassSplitEAfterCutsNLocMax1[i][j]   ->SetXTitle("E1+E2 (GeV)");
        outputContainer->Add(fhMassSplitEAfterCutsNLocMax1[i][j]) ;
        
        fhMassSplitEAfterCutsNLocMax2[i][j]     = new TH2F(Form("hMassSplitEAfterCutsNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                     Form("Mass vs E1+E2, %s %s, for NLM = 2, M02 and asy cut",ptype[i].Data(),sMatched[j].Data()),
                                                     nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassSplitEAfterCutsNLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
        fhMassSplitEAfterCutsNLocMax2[i][j]   ->SetXTitle("E1+E2 (GeV)");
        outputContainer->Add(fhMassSplitEAfterCutsNLocMax2[i][j]) ;
        
        fhMassSplitEAfterCutsNLocMaxN[i][j]     = new TH2F(Form("hMassSplitEAfterCutsNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                     Form("Mass vs E1+E2, %s %s, for NLM > 2, M02 and asy cut",ptype[i].Data(),sMatched[j].Data()),
                                                     nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassSplitEAfterCutsNLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
        fhMassSplitEAfterCutsNLocMaxN[i][j]   ->SetXTitle("E1+E2 (GeV)");
        outputContainer->Add(fhMassSplitEAfterCutsNLocMaxN[i][j]) ;

        
        fhSplitEFractionAfterCutsNLocMax1[i][j]     = new TH2F(Form("hSplitEFractionAfterCutsNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                               Form("(E1+E2)/E_{cluster} vs E_{cluster} for N max  = 1, M02 and Asy cut on, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                               nptbins,ptmin,ptmax,120,0,1.2);
        fhSplitEFractionAfterCutsNLocMax1[i][j]   ->SetXTitle("E_{cluster} (GeV)");
        fhSplitEFractionAfterCutsNLocMax1[i][j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
        outputContainer->Add(fhSplitEFractionAfterCutsNLocMax1[i][j]) ;
        
        fhSplitEFractionAfterCutsNLocMax2[i][j]     = new TH2F(Form("hSplitEFractionAfterCutsNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                               Form("(E1+E2)/E_{cluster} vs E_{cluster} for N max  = 2, M02 and Asy cut on, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                               nptbins,ptmin,ptmax,120,0,1.2);
        fhSplitEFractionAfterCutsNLocMax2[i][j]   ->SetXTitle("E_{cluster} (GeV)");
        fhSplitEFractionAfterCutsNLocMax2[i][j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
        outputContainer->Add(fhSplitEFractionAfterCutsNLocMax2[i][j]) ;
        
        fhSplitEFractionAfterCutsNLocMaxN[i][j]    = new TH2F(Form("hSplitEFractionAfterCutsNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                              Form("(E1+E2)/E_{cluster} vs E_{cluster} for N max  > 2, M02 and Asy cut on, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                              nptbins,ptmin,ptmax,120,0,1.2);
        fhSplitEFractionAfterCutsNLocMaxN[i][j]   ->SetXTitle("E_{cluster} (GeV)");
        fhSplitEFractionAfterCutsNLocMaxN[i][j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
        outputContainer->Add(fhSplitEFractionAfterCutsNLocMaxN[i][j]) ;
      }
      
      fhMassM02NLocMax1[i][j]  = new TH2F(Form("hMassM02NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Invariant mass of splitted cluster with NLM=1, #lambda_{0}^{2}, E > 12 GeV,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                          ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMax1[i][j]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMax1[i][j]) ;   
      
      fhMassM02NLocMax2[i][j]  = new TH2F(Form("hMassM02NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Invariant mass of splitted cluster with NLM=2, #lambda_{0}^{2}, E > 12 GeV, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                          ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMax2[i][j]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMax2[i][j]) ;   
      
      fhMassM02NLocMaxN[i][j]  = new TH2F(Form("hMassM02NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Invariant mass of splitted cluster with NLM>2, vs #lambda_{0}^{2}, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                          ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMaxN[i][j]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMaxN[i][j]) ;   
      
      if(fFillSSExtraHisto)
      {
        fhMassDispEtaNLocMax1[i][j]  = new TH2F(Form("hMassDispEtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of splitted cluster with NLM=1, #sigma_{#eta #eta}^{2}, E > 12 GeV,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispEtaNLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispEtaNLocMax1[i][j]->SetXTitle("#sigma_{#eta #eta}^{2}");
        outputContainer->Add(fhMassDispEtaNLocMax1[i][j]) ;   
        
        fhMassDispEtaNLocMax2[i][j]  = new TH2F(Form("hMassDispEtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of splitted cluster with NLM=2 #sigma_{#eta #eta}^{2}, E > 12 GeV, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispEtaNLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispEtaNLocMax2[i][j]->SetXTitle("#sigma_{#eta #eta}^{2}");
        outputContainer->Add(fhMassDispEtaNLocMax2[i][j]) ;   
        
        fhMassDispEtaNLocMaxN[i][j]  = new TH2F(Form("hMassDispEtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of splitted cluster with NLM>2, #sigma_{#eta #eta}^{2}, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispEtaNLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispEtaNLocMaxN[i][j]->SetXTitle("#sigma_{#eta #eta}^{2}");
        outputContainer->Add(fhMassDispEtaNLocMaxN[i][j]) ;   
        
        fhMassDispPhiNLocMax1[i][j]  = new TH2F(Form("hMassDispPhiNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of 2 highest energy cells #sigma_{#phi #phi}^{2}, E > 12 GeV,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispPhiNLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispPhiNLocMax1[i][j]->SetXTitle("#sigma_{#phi #phi}^{2}");
        outputContainer->Add(fhMassDispPhiNLocMax1[i][j]) ;   
        
        fhMassDispPhiNLocMax2[i][j]  = new TH2F(Form("hMassDispPhiNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of 2 local maxima cells #sigma_{#phi #phi}^{2}, E > 12 GeV, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispPhiNLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispPhiNLocMax2[i][j]->SetXTitle("#sigma_{#phi #phi}^{2}");
        outputContainer->Add(fhMassDispPhiNLocMax2[i][j]) ;   
        
        fhMassDispPhiNLocMaxN[i][j]  = new TH2F(Form("hMassDispPhiNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of N>2 local maxima cells vs #sigma_{#phi #phi}^{2}, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispPhiNLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispPhiNLocMaxN[i][j]->SetXTitle("#sigma_{#phi #phi}^{2}");
        outputContainer->Add(fhMassDispPhiNLocMaxN[i][j]) ;   
        
        fhMassDispAsyNLocMax1[i][j]  = new TH2F(Form("hMassDispAsyNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of 2 highest energy cells A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), E > 12 GeV,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                                200,-1,1,mbins,mmin,mmax); 
        fhMassDispAsyNLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispAsyNLocMax1[i][j]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
        outputContainer->Add(fhMassDispAsyNLocMax1[i][j]) ;   
        
        fhMassDispAsyNLocMax2[i][j]  = new TH2F(Form("hMassDispAsyNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of 2 local maxima cells A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), E > 12 GeV, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                200,-1,1,mbins,mmin,mmax); 
        fhMassDispAsyNLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispAsyNLocMax2[i][j]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
        outputContainer->Add(fhMassDispAsyNLocMax2[i][j]) ;   
        
        fhMassDispAsyNLocMaxN[i][j]  = new TH2F(Form("hMassDispAsyNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of N>2 local maxima cells vsA = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                200,-1,1,mbins,mmin,mmax); 
        fhMassDispAsyNLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispAsyNLocMaxN[i][j]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
        outputContainer->Add(fhMassDispAsyNLocMaxN[i][j]) ;   
      }
      
      
      if(i > 0 && fFillMCHisto) // skip first entry in array, general case not filled
      {
        fhMCGenFracNLocMax1[i][j]     = new TH2F(Form("hMCGenFracNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                 Form("#lambda_{0}^{2} vs E for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                 nptbins,ptmin,ptmax,200,0,2); 
        fhMCGenFracNLocMax1[i][j]   ->SetYTitle("E_{gen} / E_{reco}");
        fhMCGenFracNLocMax1[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenFracNLocMax1[i][j]) ; 
        
        fhMCGenFracNLocMax2[i][j]     = new TH2F(Form("hMCGenFracNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                 Form("#lambda_{0}^{2} vs E for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                 nptbins,ptmin,ptmax,200,0,2); 
        fhMCGenFracNLocMax2[i][j]   ->SetYTitle("E_{gen} / E_{reco}");
        fhMCGenFracNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenFracNLocMax2[i][j]) ; 
        
        fhMCGenFracNLocMaxN[i][j]    = new TH2F(Form("hMCGenFracNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("#lambda_{0}^{2} vs E for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                nptbins,ptmin,ptmax,200,0,2); 
        fhMCGenFracNLocMaxN[i][j]   ->SetYTitle("E_{gen} / E_{reco}");
        fhMCGenFracNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenFracNLocMaxN[i][j]) ; 
        
        fhMCGenFracNLocMax1NoOverlap[i][j]     = new TH2F(Form("hMCGenFracNoOverlapNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                          Form("#lambda_{0}^{2} vs E for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                          nptbins,ptmin,ptmax,200,0,2);
        fhMCGenFracNLocMax1NoOverlap[i][j]   ->SetYTitle("E_{gen} / E_{reco}");
        fhMCGenFracNLocMax1NoOverlap[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenFracNLocMax1NoOverlap[i][j]) ;
        
        fhMCGenFracNLocMax2NoOverlap[i][j]     = new TH2F(Form("hMCGenFracNoOverlapNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                          Form("#lambda_{0}^{2} vs E for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                          nptbins,ptmin,ptmax,200,0,2);
        fhMCGenFracNLocMax2NoOverlap[i][j]   ->SetYTitle("E_{gen} / E_{reco}");
        fhMCGenFracNLocMax2NoOverlap[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenFracNLocMax2NoOverlap[i][j]) ;
        
        fhMCGenFracNLocMaxNNoOverlap[i][j]    = new TH2F(Form("hMCGenFracNoOverlapNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                         Form("#lambda_{0}^{2} vs E for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                         nptbins,ptmin,ptmax,200,0,2);
        fhMCGenFracNLocMaxNNoOverlap[i][j]   ->SetYTitle("E_{gen} / E_{reco}");
        fhMCGenFracNLocMaxNNoOverlap[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenFracNLocMaxNNoOverlap[i][j]) ;
        
        
        fhMCGenSplitEFracNLocMax1[i][j]     = new TH2F(Form("hMCGenSplitEFracNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                       Form("E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                       nptbins,ptmin,ptmax,200,0,2); 
        fhMCGenSplitEFracNLocMax1[i][j]   ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
        fhMCGenSplitEFracNLocMax1[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenSplitEFracNLocMax1[i][j]) ; 
        
        fhMCGenSplitEFracNLocMax2[i][j]     = new TH2F(Form("hMCGenSplitEFracNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                       Form("E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                       nptbins,ptmin,ptmax,200,0,2); 
        fhMCGenSplitEFracNLocMax2[i][j]   ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
        fhMCGenSplitEFracNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenSplitEFracNLocMax2[i][j]) ; 
        
        fhMCGenSplitEFracNLocMaxN[i][j]    = new TH2F(Form("hMCGenSplitEFracNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                      Form("E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                      nptbins,ptmin,ptmax,200,0,2); 
        fhMCGenSplitEFracNLocMaxN[i][j]   ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
        fhMCGenSplitEFracNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenSplitEFracNLocMaxN[i][j]) ; 
        
        fhMCGenSplitEFracNLocMax1NoOverlap[i][j]     = new TH2F(Form("hMCGenSplitEFracNoOverlapNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                                Form("E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                                nptbins,ptmin,ptmax,200,0,2);
        fhMCGenSplitEFracNLocMax1NoOverlap[i][j]   ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
        fhMCGenSplitEFracNLocMax1NoOverlap[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenSplitEFracNLocMax1NoOverlap[i][j]) ;
        
        fhMCGenSplitEFracNLocMax2NoOverlap[i][j]     = new TH2F(Form("hMCGenSplitEFracNoOverlapNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                                Form("E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                                nptbins,ptmin,ptmax,200,0,2);
        fhMCGenSplitEFracNLocMax2NoOverlap[i][j]   ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
        fhMCGenSplitEFracNLocMax2NoOverlap[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenSplitEFracNLocMax2NoOverlap[i][j]) ;
        
        fhMCGenSplitEFracNLocMaxNNoOverlap[i][j]    = new TH2F(Form("hMCGenSplitEFracNoOverlapNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                               Form("E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                               nptbins,ptmin,ptmax,200,0,2);
        fhMCGenSplitEFracNLocMaxNNoOverlap[i][j]   ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
        fhMCGenSplitEFracNLocMaxNNoOverlap[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenSplitEFracNLocMaxNNoOverlap[i][j]) ;
        
        fhMCGenEFracvsSplitEFracNLocMax1[i][j]     = new TH2F(Form("hMCGenEFracvsSplitEFracNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                              Form("(E_{1 split}+E_{2 split})/E_{reco} vs E_{gen} / E_{reco} for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                              200,0,2,200,0,2); 
        fhMCGenEFracvsSplitEFracNLocMax1[i][j]   ->SetYTitle("(E_{1 split}+E_{2 split})/E_{reco}");
        fhMCGenEFracvsSplitEFracNLocMax1[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
        outputContainer->Add(fhMCGenEFracvsSplitEFracNLocMax1[i][j]) ; 
        
        fhMCGenEFracvsSplitEFracNLocMax2[i][j]     = new TH2F(Form("hMCGenEFracvsSplitEFracNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                              Form("(E_{1 split}+E_{2 split})/E_{reco} vs E_{gen} / E_{reco} for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                              200,0,2,200,0,2); 
        fhMCGenEFracvsSplitEFracNLocMax2[i][j]   ->SetYTitle("(E_{1 split}+E_{2 split})/E_{reco}");
        fhMCGenEFracvsSplitEFracNLocMax2[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
        outputContainer->Add(fhMCGenEFracvsSplitEFracNLocMax2[i][j]) ; 
        
        
        fhMCGenEFracvsSplitEFracNLocMaxN[i][j]    = new TH2F(Form("hMCGenEFracvsSplitEFracNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                             Form("(E_{1 split}+E_{2 split})/E_{reco} vs E_{gen} / E_{reco} for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                             200,0,2,200,0,2); 
        fhMCGenEFracvsSplitEFracNLocMaxN[i][j]   ->SetYTitle("(E_{1 split}+E_{2 split})/E_{reco}");
        fhMCGenEFracvsSplitEFracNLocMaxN[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
        outputContainer->Add(fhMCGenEFracvsSplitEFracNLocMaxN[i][j]) ; 
        
        
        fhMCGenEvsSplitENLocMax1[i][j]     = new TH2F(Form("hMCGenEvsSplitENLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                      Form("E_{1 split}+E_{2 split} vs E_{gen} for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                      nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
        fhMCGenEvsSplitENLocMax1[i][j]   ->SetYTitle("E_{1 split}+E_{2 split} (GeV)");
        fhMCGenEvsSplitENLocMax1[i][j]   ->SetXTitle("E_{gen} (GeV)");
        outputContainer->Add(fhMCGenEvsSplitENLocMax1[i][j]) ; 
        
        fhMCGenEvsSplitENLocMax2[i][j]     = new TH2F(Form("hMCGenEvsSplitENLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                      Form("E_{1 split}+E_{2 split} vs E_{gen} for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                      nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
        fhMCGenEvsSplitENLocMax2[i][j]   ->SetYTitle("E_{1 split}+E_{2 split} (GeV)");
        fhMCGenEvsSplitENLocMax2[i][j]   ->SetXTitle("E_{gen} (GeV)");
        outputContainer->Add(fhMCGenEvsSplitENLocMax2[i][j]) ; 
        
        
        fhMCGenEvsSplitENLocMaxN[i][j]    = new TH2F(Form("hMCGenEvsSplitENLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                     Form("E_{1 split}+E_{2 split} vs E_{gen} for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                     nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
        fhMCGenEvsSplitENLocMaxN[i][j]   ->SetYTitle("E_{1 split}+E_{2 split} (GeV)");
        fhMCGenEvsSplitENLocMaxN[i][j]   ->SetXTitle("E_{gen} (GeV)");
        outputContainer->Add(fhMCGenEvsSplitENLocMaxN[i][j]) ; 
      }
          
      // Histograms after cluster identification
      
      
      // Pi0 //
      
      fhM02Pi0NLocMax1[i][j]     = new TH2F(Form("hM02Pi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("#lambda_{0}^{2} vs E, %s, for NLM = 1",ptype[i].Data()),
                                            nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhM02Pi0NLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0NLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0NLocMax1[i][j]) ;
      
      fhM02Pi0NLocMax2[i][j]     = new TH2F(Form("hM02Pi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("#lambda_{0}^{2} vs E, %s, for NLM = 2",ptype[i].Data()),
                                            nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02Pi0NLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0NLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0NLocMax2[i][j]) ; 
      
      fhM02Pi0NLocMaxN[i][j]     = new TH2F(Form("hM02Pi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("#lambda_{0}^{2} vs E, %s, for NLM > 2",ptype[i].Data()),
                                            nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02Pi0NLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0NLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0NLocMaxN[i][j]) ; 
      
      fhMassPi0NLocMax1[i][j]     = new TH2F(Form("hMassPi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("Mass vs E, %s, for NLM = 1",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassPi0NLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassPi0NLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassPi0NLocMax1[i][j]) ; 
      
      fhMassPi0NLocMax2[i][j]     = new TH2F(Form("hMassPi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("Mass vs E , %s, for NLM = 2",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassPi0NLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassPi0NLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassPi0NLocMax2[i][j]) ; 
      
      fhMassPi0NLocMaxN[i][j]     = new TH2F(Form("hMassPi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("Mass vs E, %s, for NLM > 2",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassPi0NLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassPi0NLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassPi0NLocMaxN[i][j]) ; 

      fhMassSplitEPi0NLocMax1[i][j]     = new TH2F(Form("hMassSplitEPi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("Mass vs E1+E2, %s, for NLM = 1",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassSplitEPi0NLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassSplitEPi0NLocMax1[i][j]   ->SetXTitle("E1+E2 (GeV)");
      outputContainer->Add(fhMassSplitEPi0NLocMax1[i][j]) ;
      
      fhMassSplitEPi0NLocMax2[i][j]     = new TH2F(Form("hMassSplitEPi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("Mass vs E1+E2 , %s, for NLM = 2",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassSplitEPi0NLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassSplitEPi0NLocMax2[i][j]   ->SetXTitle("E1+E2 (GeV)");
      outputContainer->Add(fhMassSplitEPi0NLocMax2[i][j]) ;
      
      fhMassSplitEPi0NLocMaxN[i][j]     = new TH2F(Form("hMassSplitEPi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("Mass vs E1+E2, %s, for NLM > 2",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassSplitEPi0NLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassSplitEPi0NLocMaxN[i][j]   ->SetXTitle("E1+E2 (GeV)");
      outputContainer->Add(fhMassSplitEPi0NLocMaxN[i][j]) ;
      
      fhAsyPi0NLocMax1[i][j]     = new TH2F(Form("hAsyPi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E, %s, for NLM = 1",ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyPi0NLocMax1[i][j]   ->SetYTitle("Asymmetry");
      fhAsyPi0NLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyPi0NLocMax1[i][j]) ; 
      
      fhAsyPi0NLocMax2[i][j]     = new TH2F(Form("hAsyPi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E, %s, for NLM = 2",ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyPi0NLocMax2[i][j]   ->SetYTitle("Asymmetry");
      fhAsyPi0NLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyPi0NLocMax2[i][j]) ; 
      
      fhAsyPi0NLocMaxN[i][j]     = new TH2F(Form("hAsyPi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E, %s, for NLM > 2",ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyPi0NLocMaxN[i][j]   ->SetYTitle("Asymmetry");
      fhAsyPi0NLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyPi0NLocMaxN[i][j]) ; 
      
      if(fFillNCellHisto)
      {
        fhNCellPi0NLocMax1[i][j]     = new TH2F(Form("hNCellPi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("n cells vs E, %s, for NLM = 1",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
        fhNCellPi0NLocMax1[i][j]   ->SetYTitle("n cells");
        fhNCellPi0NLocMax1[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhNCellPi0NLocMax1[i][j]) ;
        
        fhNCellPi0NLocMax2[i][j]     = new TH2F(Form("hNCellPi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("n cells vs E, %s, for NLM = 2",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
        fhNCellPi0NLocMax2[i][j]   ->SetYTitle("n cells");
        fhNCellPi0NLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhNCellPi0NLocMax2[i][j]) ;
        
        fhNCellPi0NLocMaxN[i][j]     = new TH2F(Form("hNCellPi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("n cells vs E, %s, for NLM > 2",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
        fhNCellPi0NLocMaxN[i][j]   ->SetYTitle("n cells");
        fhNCellPi0NLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhNCellPi0NLocMaxN[i][j]) ;
      }
      
      // Eta
    
      if(fFillIdEtaHisto)
      {
        fhM02EtaNLocMax1[i][j]     = new TH2F(Form("hM02EtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("#lambda_{0}^{2} vs E, %s, for NLM = 1",ptype[i].Data()),
                                              nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhM02EtaNLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
        fhM02EtaNLocMax1[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhM02EtaNLocMax1[i][j]) ;
        
        
        fhM02EtaNLocMax2[i][j]     = new TH2F(Form("hM02EtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("#lambda_{0}^{2} vs E, %s, for NLM = 2",ptype[i].Data()),
                                              nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhM02EtaNLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
        fhM02EtaNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhM02EtaNLocMax2[i][j]) ;
        
        fhM02EtaNLocMaxN[i][j]     = new TH2F(Form("hM02EtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("#lambda_{0}^{2} vs E, %s, for NLM > 2",ptype[i].Data()),
                                              nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhM02EtaNLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
        fhM02EtaNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhM02EtaNLocMaxN[i][j]) ;
        
        fhMassEtaNLocMax1[i][j]     = new TH2F(Form("hMassEtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                               Form("Mass vs E, %s, for NLM = 1",ptype[i].Data()),
                                               nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassEtaNLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
        fhMassEtaNLocMax1[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassEtaNLocMax1[i][j]) ;
        
        fhMassEtaNLocMax2[i][j]     = new TH2F(Form("hMassEtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                               Form("Mass vs E, %s, for NLM = 2",ptype[i].Data()),
                                               nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassEtaNLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
        fhMassEtaNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassEtaNLocMax2[i][j]) ;
        
        fhMassEtaNLocMaxN[i][j]     = new TH2F(Form("hMassEtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                               Form("Mass vs E, %s, for NLM > 2",ptype[i].Data()),
                                               nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassEtaNLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
        fhMassEtaNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassEtaNLocMaxN[i][j]) ;
        
        fhAsyEtaNLocMax1[i][j]     = new TH2F(Form("hAsyEtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("Asymmetry vs E, %s, for NLM = 1",ptype[i].Data()),
                                              nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhAsyEtaNLocMax1[i][j]   ->SetYTitle("Asymmetry");
        fhAsyEtaNLocMax1[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhAsyEtaNLocMax1[i][j]) ;
        
        fhAsyEtaNLocMax2[i][j]     = new TH2F(Form("hAsyEtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("Asymmetry vs E, %s, for NLM = 2",ptype[i].Data()),
                                              nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhAsyEtaNLocMax2[i][j]   ->SetYTitle("Asymmetry");
        fhAsyEtaNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhAsyEtaNLocMax2[i][j]) ;
        
        fhAsyEtaNLocMaxN[i][j]     = new TH2F(Form("hAsyEtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("Asymmetry vs E, %s, for NLM > 2",ptype[i].Data()),
                                              nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhAsyEtaNLocMaxN[i][j]   ->SetYTitle("Asymmetry");
        fhAsyEtaNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhAsyEtaNLocMaxN[i][j]) ;
        
        if(fFillNCellHisto)
        {
          fhNCellEtaNLocMax1[i][j]     = new TH2F(Form("hNCellEtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                  Form("n cells vs E, %s, for NLM = 1",ptype[i].Data()),
                                                  nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
          fhNCellEtaNLocMax1[i][j]   ->SetYTitle("n cells");
          fhNCellEtaNLocMax1[i][j]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhNCellEtaNLocMax1[i][j]) ;
          
          fhNCellEtaNLocMax2[i][j]     = new TH2F(Form("hNCellEtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                  Form("n cells vs E, %s, for NLM = 2",ptype[i].Data()),
                                                  nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
          fhNCellEtaNLocMax2[i][j]   ->SetYTitle("n cells");
          fhNCellEtaNLocMax2[i][j]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhNCellEtaNLocMax2[i][j]) ;
          
          fhNCellEtaNLocMaxN[i][j]     = new TH2F(Form("hNCellEtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                  Form("n cells vs E, %s, for NLM > 2",ptype[i].Data()),
                                                  nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
          fhNCellEtaNLocMaxN[i][j]   ->SetYTitle("n cells");
          fhNCellEtaNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhNCellEtaNLocMaxN[i][j]) ;
          
        }
      }
      
      if(fFillIdConvHisto)
      {
        fhM02ConNLocMax1[i][j]    = new TH2F(Form("hM02ConNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("#lambda_{0}^{2} vs E, %s, for NLM = 1",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhM02ConNLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
        fhM02ConNLocMax1[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhM02ConNLocMax1[i][j]) ;
        
        fhM02ConNLocMax2[i][j]    = new TH2F(Form("hM02ConNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("#lambda_{0}^{2} vs E, %s, for NLM = 2",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhM02ConNLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
        fhM02ConNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhM02ConNLocMax2[i][j]) ;
        
        fhM02ConNLocMaxN[i][j]    = new TH2F(Form("hM02ConNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("#lambda_{0}^{2} vs E, %s, for NLM > 2",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhM02ConNLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
        fhM02ConNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhM02ConNLocMaxN[i][j]) ;
        
        
        fhMassConNLocMax1[i][j]    = new TH2F(Form("hMassConNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("Mass vs E, %s, for NLM = 1",ptype[i].Data()),
                                              nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassConNLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
        fhMassConNLocMax1[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassConNLocMax1[i][j]) ;
        
        fhMassConNLocMax2[i][j]    = new TH2F(Form("hMassConNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("Mass vs E, %s, for NLM = 2",ptype[i].Data()),
                                              nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassConNLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
        fhMassConNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassConNLocMax2[i][j]) ;
        
        fhMassConNLocMaxN[i][j]    = new TH2F(Form("hMassConNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("Mass vs E, %s, for NLM > 2",ptype[i].Data()),
                                              nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassConNLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
        fhMassConNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassConNLocMaxN[i][j]) ;
        
        fhAsyConNLocMax1[i][j]    = new TH2F(Form("hAsyConNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("Asymmetry vs E, %s, for NLM = 1",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhAsyConNLocMax1[i][j]   ->SetYTitle("Asymmetry");
        fhAsyConNLocMax1[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhAsyConNLocMax1[i][j]) ;
        
        fhAsyConNLocMax2[i][j]    = new TH2F(Form("hAsyConNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("Asymmetry vs E, %s, for NLM = 2",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhAsyConNLocMax2[i][j]   ->SetYTitle("Asymmetry");
        fhAsyConNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhAsyConNLocMax2[i][j]) ;
        
        fhAsyConNLocMaxN[i][j]    = new TH2F(Form("hAsyConNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("Asymmetry vs E, %s, for NLM > 2",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhAsyConNLocMaxN[i][j]   ->SetYTitle("Asymmetry");
        fhAsyConNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhAsyConNLocMaxN[i][j]) ;
        
      }
      
    } // matched, not matched
    
    if(fFillEbinHisto)
    {
      for(Int_t j = 0; j < 4; j++)
      {
        
        fhMassSplitEFractionNLocMax1Ebin[i][j]  = new TH2F(Form("hMassSplitEFractionNLocMax1%sEbin%d",pname[i].Data(),j),
                                                           Form("Invariant mass of 2 highest energy cells vs (E1+E2)/Ecluster, %s, %s",ptype[i].Data(),sEBin[j].Data()),
                                                           120,0,1.2,mbins,mmin,mmax);
        fhMassSplitEFractionNLocMax1Ebin[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassSplitEFractionNLocMax1Ebin[i][j]->SetXTitle("(E_{split1}+E_{split2})/E_{cluster}");
        outputContainer->Add(fhMassSplitEFractionNLocMax1Ebin[i][j]) ;
        
        fhMassSplitEFractionNLocMax2Ebin[i][j]  = new TH2F(Form("hMassSplitEFractionNLocMax2%sEbin%d",pname[i].Data(),j),
                                                           Form("Invariant mass of 2 local maxima cells vs (E1+E2)/Ecluster, %s, %s",ptype[i].Data(),sEBin[j].Data()),
                                                           120,0,1.2,mbins,mmin,mmax);
        fhMassSplitEFractionNLocMax2Ebin[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassSplitEFractionNLocMax2Ebin[i][j]->SetXTitle("(E_{split1}+E_{split2})/E_{cluster}");
        outputContainer->Add(fhMassSplitEFractionNLocMax2Ebin[i][j]) ;
        
        fhMassSplitEFractionNLocMaxNEbin[i][j]  = new TH2F(Form("hMassSplitEFractionNLocMaxN%sEbin%d",pname[i].Data(),j),
                                                           Form("Invariant mass of N>2 local maxima cells vs (E1+E2)/Ecluster, %s, %s",ptype[i].Data(),sEBin[j].Data()),
                                                           120,0,1.2,mbins,mmin,mmax);
        fhMassSplitEFractionNLocMaxNEbin[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassSplitEFractionNLocMaxNEbin[i][j]->SetXTitle("(E_{split1}+E_{split2})/E_{cluster}");
        outputContainer->Add(fhMassSplitEFractionNLocMaxNEbin[i][j]) ;
        
        if(i>0 && fFillMCHisto) // skip first entry in array, general case not filled
        {
          fhMCGenFracNLocMaxEbin[i][j]  = new TH2F(Form("hMCGenFracNLocMax%sEbin%d",pname[i].Data(),j),
                                                   Form("NLM vs E, %s, %s",ptype[i].Data(),sEBin[j].Data()),
                                                   200,0,2,nMaxBins,0,nMaxBins);
          fhMCGenFracNLocMaxEbin[i][j]->SetYTitle("NLM");
          fhMCGenFracNLocMaxEbin[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMCGenFracNLocMaxEbin[i][j]) ;
          
          fhMCGenFracNLocMaxEbinMatched[i][j]  = new TH2F(Form("hMCGenFracNLocMax%sEbin%dMatched",pname[i].Data(),j),
                                                          Form("NLM vs E, %s, %s, matched to a track",ptype[i].Data(),sEBin[j].Data()),
                                                          200,0,2,nMaxBins,0,nMaxBins);
          fhMCGenFracNLocMaxEbinMatched[i][j]->SetYTitle("NLM");
          fhMCGenFracNLocMaxEbinMatched[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMCGenFracNLocMaxEbinMatched[i][j]) ;
          
          fhMassMCGenFracNLocMax1Ebin[i][j]  = new TH2F(Form("hMassMCGenFracNLocMax1%sEbin%d",pname[i].Data(),j),
                                                        Form("Invariant mass of 2 highest energy cells vs E, %s, %s",ptype[i].Data(),sEBin[j].Data()),
                                                        200,0,2,mbins,mmin,mmax);
          fhMassMCGenFracNLocMax1Ebin[i][j]->SetYTitle("M (GeV/c^{2})");
          fhMassMCGenFracNLocMax1Ebin[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMassMCGenFracNLocMax1Ebin[i][j]) ;
          
          fhMassMCGenFracNLocMax2Ebin[i][j]  = new TH2F(Form("hMassMCGenFracNLocMax2%sEbin%d",pname[i].Data(),j),
                                                        Form("Invariant mass of 2 local maxima cells vs E, %s, %s",ptype[i].Data(),sEBin[j].Data()),
                                                        200,0,2,mbins,mmin,mmax);
          fhMassMCGenFracNLocMax2Ebin[i][j]->SetYTitle("M (GeV/c^{2})");
          fhMassMCGenFracNLocMax2Ebin[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMassMCGenFracNLocMax2Ebin[i][j]) ;
          
          fhMassMCGenFracNLocMaxNEbin[i][j]  = new TH2F(Form("hMassMCGenFracNLocMaxN%sEbin%d",pname[i].Data(),j),
                                                        Form("Invariant mass of N>2 local maxima cells vs E, %s, %s",ptype[i].Data(),sEBin[j].Data()),
                                                        200,0,2,mbins,mmin,mmax);
          fhMassMCGenFracNLocMaxNEbin[i][j]->SetYTitle("M (GeV/c^{2})");
          fhMassMCGenFracNLocMaxNEbin[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMassMCGenFracNLocMaxNEbin[i][j]) ;
          
          fhM02MCGenFracNLocMax1Ebin[i][j]     = new TH2F(Form("hM02MCGenFracNLocMax1%sEbin%d",pname[i].Data(),j),
                                                          Form("#lambda_{0}^{2} vs E for N max  = 1 %s, %s",ptype[i].Data(),sEBin[j].Data()),
                                                          200,0,2,ssbins,ssmin,ssmax);
          fhM02MCGenFracNLocMax1Ebin[i][j]   ->SetYTitle("#lambda_{0}^{2}");
          fhM02MCGenFracNLocMax1Ebin[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhM02MCGenFracNLocMax1Ebin[i][j]) ;
          
          fhM02MCGenFracNLocMax2Ebin[i][j]     = new TH2F(Form("hM02MCGenFracNLocMax2%sEbin%d",pname[i].Data(),j),
                                                          Form("#lambda_{0}^{2} vs E for N max  = 2 %s, %s",ptype[i].Data(),sEBin[j].Data()),
                                                          200,0,2,ssbins,ssmin,ssmax);
          fhM02MCGenFracNLocMax2Ebin[i][j]   ->SetYTitle("#lambda_{0}^{2}");
          fhM02MCGenFracNLocMax2Ebin[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhM02MCGenFracNLocMax2Ebin[i][j]) ;
          
          fhM02MCGenFracNLocMaxNEbin[i][j]    = new TH2F(Form("hM02MCGenFracNLocMaxN%sEbin%d",pname[i].Data(),j),
                                                         Form("#lambda_{0}^{2} vs E for N max  > 2 %s, %s",ptype[i].Data(),sEBin[j].Data()),
                                                         200,0,2,ssbins,ssmin,ssmax);
          fhM02MCGenFracNLocMaxNEbin[i][j]   ->SetYTitle("#lambda_{0}^{2}");
          fhM02MCGenFracNLocMaxNEbin[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhM02MCGenFracNLocMaxNEbin[i][j]) ;
        }
      }
    }
  } // MC particle list
  
  if(fFillHighMultHisto)
  {
    // E vs centrality

    fhCentralityPi0NLocMax1  = new TH2F("hCentralityPi0NLocMax1",
                                        "E vs Centrality, selected pi0 cluster with NLM=1",
                                        nptbins,ptmin,ptmax,100,0,100);
    fhCentralityPi0NLocMax1->SetYTitle("Centrality");
    fhCentralityPi0NLocMax1->SetXTitle("E (GeV)");
    outputContainer->Add(fhCentralityPi0NLocMax1) ;
    
    fhCentralityPi0NLocMax2  = new TH2F("hCentralityPi0NLocMax2",
                                        "E vs Centrality, selected pi0 cluster with NLM=2",
                                        nptbins,ptmin,ptmax,100,0,100);
    fhCentralityPi0NLocMax2->SetYTitle("Centrality");
    fhCentralityPi0NLocMax2->SetXTitle("E (GeV)");
    outputContainer->Add(fhCentralityPi0NLocMax2) ;
    
    fhCentralityPi0NLocMaxN  = new TH2F("hCentralityPi0NLocMaxN",
                                        "E vs Centrality, selected pi0 cluster with NLM>1",
                                        nptbins,ptmin,ptmax,100,0,100);
    fhCentralityPi0NLocMaxN->SetYTitle("Centrality");
    fhCentralityPi0NLocMaxN->SetXTitle("E (GeV)");
    outputContainer->Add(fhCentralityPi0NLocMaxN) ;
    
    if(fFillIdEtaHisto)
    {
      fhCentralityEtaNLocMax1  = new TH2F("hCentralityEtaNLocMax1",
                                          "E vs Centrality, selected pi0 cluster with NLM=1",
                                          nptbins,ptmin,ptmax,100,0,100);
      fhCentralityEtaNLocMax1->SetYTitle("Centrality");
      fhCentralityEtaNLocMax1->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityEtaNLocMax1) ;
      
      fhCentralityEtaNLocMax2  = new TH2F("hCentralityEtaNLocMax2",
                                          "E vs Centrality, selected pi0 cluster with NLM=2",
                                          nptbins,ptmin,ptmax,100,0,100);
      fhCentralityEtaNLocMax2->SetYTitle("Centrality");
      fhCentralityEtaNLocMax2->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityEtaNLocMax2) ;
      
      fhCentralityEtaNLocMaxN  = new TH2F("hCentralityEtaNLocMaxN",
                                          "E vs Centrality, selected pi0 cluster with NLM>1",
                                          nptbins,ptmin,ptmax,100,0,100);
      fhCentralityEtaNLocMaxN->SetYTitle("Centrality");
      fhCentralityEtaNLocMaxN->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityEtaNLocMaxN) ;
    }
    
    // E vs Event plane angle
    
    fhEventPlanePi0NLocMax1  = new TH2F("hEventPlanePi0NLocMax1","E vs Event Plane Angle, selected pi0 cluster with NLM=1",
                                        nptbins,ptmin,ptmax,100,0,TMath::Pi());
    fhEventPlanePi0NLocMax1->SetYTitle("Event Plane Angle (rad)");
    fhEventPlanePi0NLocMax1->SetXTitle("E (GeV)");
    outputContainer->Add(fhEventPlanePi0NLocMax1) ;
    
    fhEventPlanePi0NLocMax2  = new TH2F("hEventPlanePi0NLocMax2","E vs Event Plane Angle, selected pi0 cluster with NLM=2",
                                        nptbins,ptmin,ptmax,100,0,TMath::Pi());
    fhEventPlanePi0NLocMax2->SetYTitle("Event Plane Angle (rad)");
    fhEventPlanePi0NLocMax2->SetXTitle("E (GeV)");
    outputContainer->Add(fhEventPlanePi0NLocMax2) ;
    
    fhEventPlanePi0NLocMaxN  = new TH2F("hEventPlanePi0NLocMaxN","E vs Event Plane Angle, selected pi0 cluster with NLM>1",
                                        nptbins,ptmin,ptmax,100,0,TMath::Pi());
    fhEventPlanePi0NLocMaxN->SetYTitle("Event Plane Angle (rad)");
    fhEventPlanePi0NLocMaxN->SetXTitle("E (GeV)");
    outputContainer->Add(fhEventPlanePi0NLocMaxN) ;
    
    if(fFillIdEtaHisto)
    {
      fhEventPlaneEtaNLocMax1  = new TH2F("hEventPlaneEtaNLocMax1","E vs Event Plane Angle, selected pi0 cluster with NLM=1",
                                          nptbins,ptmin,ptmax,100,0,TMath::Pi());
      fhEventPlaneEtaNLocMax1->SetYTitle("Event Plane Angle (rad)");
      fhEventPlaneEtaNLocMax1->SetXTitle("E (GeV)");
      outputContainer->Add(fhEventPlaneEtaNLocMax1) ;
      
      fhEventPlaneEtaNLocMax2  = new TH2F("hEventPlaneEtaNLocMax2","E vs Event Plane Angle, selected pi0 cluster with NLM=2",
                                          nptbins,ptmin,ptmax,100,0,TMath::Pi());
      fhEventPlaneEtaNLocMax2->SetYTitle("Event Plane Angle (rad)");
      fhEventPlaneEtaNLocMax2->SetXTitle("E (GeV)");
      outputContainer->Add(fhEventPlaneEtaNLocMax2) ;
      
      fhEventPlaneEtaNLocMaxN  = new TH2F("hEventPlaneEtaNLocMaxN","E vs Event Plane Angle, selected pi0 cluster with NLM>1",
                                          nptbins,ptmin,ptmax,100,0,TMath::Pi());
      fhEventPlaneEtaNLocMaxN->SetYTitle("Event Plane Angle (rad)");
      fhEventPlaneEtaNLocMaxN->SetXTitle("E (GeV)");
      outputContainer->Add(fhEventPlaneEtaNLocMaxN) ;
    }
  }
  
  if(fFillEbinHisto)
  {
    for(Int_t i = 0; i < 4; i++)
    {
      fhMassM02NLocMax1Ebin[i]  = new TH2F(Form("hMassM02NLocMax1Ebin%d",i),
                                           Form("Invariant mass of split clusters vs #lambda_{0}^{2}, NLM=1, %s",sEBin[i].Data()),
                                           ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMax1Ebin[i]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMax1Ebin[i]) ;   
      
      fhMassM02NLocMax2Ebin[i]  = new TH2F(Form("hMassM02NLocMax2Ebin%d",i),
                                           Form("Invariant mass of split clusters vs #lambda_{0}^{2}, NLM=2, %s",sEBin[i].Data()),
                                           ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMax2Ebin[i]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMax2Ebin[i]) ;   
      
      fhMassM02NLocMaxNEbin[i]  = new TH2F(Form("hMassM02NLocMaxNEbin%d",i),
                                           Form("Invariant mass of split clusters vs vs #lambda_{0}^{2}, NLM>2, %s",sEBin[i].Data()),
                                           ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMaxNEbin[i]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMaxNEbin[i]) ; 
      
      
      fhMassAsyNLocMax1Ebin[i]  = new TH2F(Form("hMassAsyNLocMax1Ebin%d",i),
                                           Form("Invariant mass of split clusters vs split asymmetry, NLM=1, %s",sEBin[i].Data()),
                                           200,-1,1,mbins,mmin,mmax);
      fhMassAsyNLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassAsyNLocMax1Ebin[i]->SetXTitle("asymmetry");
      outputContainer->Add(fhMassAsyNLocMax1Ebin[i]) ;
      
      fhMassAsyNLocMax2Ebin[i]  = new TH2F(Form("hMassAsyNLocMax2Ebin%d",i),
                                           Form("Invariant mass of split clusters vs split asymmetry, NLM=2, %s",sEBin[i].Data()),
                                           200,-1,1,mbins,mmin,mmax);
      fhMassAsyNLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassAsyNLocMax2Ebin[i]->SetXTitle("asymmetry");
      outputContainer->Add(fhMassAsyNLocMax2Ebin[i]) ;
      
      fhMassAsyNLocMaxNEbin[i]  = new TH2F(Form("hMassAsyNLocMaxNEbin%d",i),
                                           Form("Invariant mass of split clusters vs split asymmetry, NLM>2, %s",sEBin[i].Data()),
                                           200,-1,1,mbins,mmin,mmax);
      fhMassAsyNLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassAsyNLocMaxNEbin[i]->SetXTitle("asymmetry");
      outputContainer->Add(fhMassAsyNLocMaxNEbin[i]) ;
      
      
      if(IsDataMC() && fFillMCHisto)
      {
        fhMCAsymM02NLocMax1MCPi0Ebin[i]  = new TH2F(Form("hMCAsymM02NLocMax1MCPi0Ebin%d",i),
                                                    Form("Asymmetry of MC #pi^{0} vs #lambda_{0}^{2}, NLM=1, %s",sEBin[i].Data()),
                                                    ssbins,ssmin,ssmax,100,0,1);
        fhMCAsymM02NLocMax1MCPi0Ebin[i]->SetYTitle("Decay asymmetry");
        fhMCAsymM02NLocMax1MCPi0Ebin[i]->SetXTitle("#lambda_{0}^{2}");
        outputContainer->Add(fhMCAsymM02NLocMax1MCPi0Ebin[i]) ;
        
        fhMCAsymM02NLocMax2MCPi0Ebin[i]  = new TH2F(Form("hMCAsymM02NLocMax2MCPi0Ebin%d",i),
                                                    Form("Asymmetry of MC #pi^{0} vs #lambda_{0}^{2}, NLM=2, %s",sEBin[i].Data()),
                                                    ssbins,ssmin,ssmax,100,0,1);
        fhMCAsymM02NLocMax2MCPi0Ebin[i]->SetYTitle("Decay asymmetry");
        fhMCAsymM02NLocMax2MCPi0Ebin[i]->SetXTitle("#lambda_{0}^{2}");
        outputContainer->Add(fhMCAsymM02NLocMax2MCPi0Ebin[i]) ;
        
        fhMCAsymM02NLocMaxNMCPi0Ebin[i]  = new TH2F(Form("hMCAsymM02NLocMaxNMCPi0Ebin%d",i),
                                                    Form("Asymmetry of MC #pi^{0} vs #lambda_{0}^{2}, NLM>2, %s",sEBin[i].Data()),
                                                    ssbins,ssmin,ssmax,100,0,1);
        fhMCAsymM02NLocMaxNMCPi0Ebin[i]->SetYTitle("Decay asymmetry");
        fhMCAsymM02NLocMaxNMCPi0Ebin[i]->SetXTitle("#lambda_{0}^{2}");
        outputContainer->Add(fhMCAsymM02NLocMaxNMCPi0Ebin[i]) ;    
        
        
        fhAsyMCGenRecoNLocMax1EbinPi0[i]  = new TH2F(Form("hAsyMCGenRecoNLocMax1Ebin%dPi0",i),
                                                     Form("Generated vs reconstructed asymmetry of split clusters from pi0, NLM=1, %s",sEBin[i].Data()),
                                                     200,-1,1,200,-1,1);
        fhAsyMCGenRecoNLocMax1EbinPi0[i]->SetYTitle("M (GeV/c^{2})");
        fhAsyMCGenRecoNLocMax1EbinPi0[i]->SetXTitle("asymmetry");
        outputContainer->Add(fhAsyMCGenRecoNLocMax1EbinPi0[i]) ;
        
        fhAsyMCGenRecoNLocMax2EbinPi0[i]  = new TH2F(Form("hAsyMCGenRecoNLocMax2Ebin%dPi0",i),
                                                     Form("Generated vs reconstructed asymmetry of split clusters from pi0, NLM=2, %s",sEBin[i].Data()),
                                                     200,-1,1,200,-1,1);
        fhAsyMCGenRecoNLocMax2EbinPi0[i]->SetYTitle("M (GeV/c^{2})");
        fhAsyMCGenRecoNLocMax2EbinPi0[i]->SetXTitle("asymmetry");
        outputContainer->Add(fhAsyMCGenRecoNLocMax2EbinPi0[i]) ;
        
        fhAsyMCGenRecoNLocMaxNEbinPi0[i]  = new TH2F(Form("hAsyMCGenRecoNLocMaxNEbin%dPi0",i),
                                                     Form("Generated vs reconstructed asymmetry of split clusters from pi0, NLM>2, %s",sEBin[i].Data()),
                                                     200,-1,1,200,-1,1);
        fhAsyMCGenRecoNLocMaxNEbinPi0[i]->SetYTitle("M (GeV/c^{2})");
        fhAsyMCGenRecoNLocMaxNEbinPi0[i]->SetXTitle("asymmetry");
        outputContainer->Add(fhAsyMCGenRecoNLocMaxNEbinPi0[i]) ;
      }
      
      if(fFillSSExtraHisto)
      {
        fhMassDispEtaNLocMax1Ebin[i]  = new TH2F(Form("hMassDispEtaNLocMax1Ebin%d",i),
                                                 Form("Invariant mass of 2 highest energy cells #sigma_{#eta #eta}^{2}, %s",sEBin[i].Data()),
                                                 ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispEtaNLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispEtaNLocMax1Ebin[i]->SetXTitle("#sigma_{#eta #eta}^{2}");
        outputContainer->Add(fhMassDispEtaNLocMax1Ebin[i]) ;   
        
        fhMassDispEtaNLocMax2Ebin[i]  = new TH2F(Form("hMassDispEtaNLocMax2Ebin%d",i),
                                                 Form("Invariant mass of 2 local maxima cells #sigma_{#eta #eta}^{2}, %s",sEBin[i].Data()),
                                                 ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispEtaNLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispEtaNLocMax2Ebin[i]->SetXTitle("#sigma_{#eta #eta}^{2}");
        outputContainer->Add(fhMassDispEtaNLocMax2Ebin[i]) ;   
        
        fhMassDispEtaNLocMaxNEbin[i]  = new TH2F(Form("hMassDispEtaNLocMaxNEbin%d",i),
                                                 Form("Invariant mass of N>2 local maxima cells vs #sigma_{#eta #eta}^{2}, %s",sEBin[i].Data()),
                                                 ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispEtaNLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispEtaNLocMaxNEbin[i]->SetXTitle("#sigma_{#eta #eta}^{2}");
        outputContainer->Add(fhMassDispEtaNLocMaxNEbin[i]) ;   
        
        fhMassDispPhiNLocMax1Ebin[i]  = new TH2F(Form("hMassDispPhiNLocMax1Ebin%d",i),
                                                 Form("Invariant mass of 2 highest energy cells #sigma_{#phi #phi}^{2}, %s",sEBin[i].Data()),
                                                 ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispPhiNLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispPhiNLocMax1Ebin[i]->SetXTitle("#sigma_{#phi #phi}^{2}");
        outputContainer->Add(fhMassDispPhiNLocMax1Ebin[i]) ;   
        
        fhMassDispPhiNLocMax2Ebin[i]  = new TH2F(Form("hMassDispPhiNLocMax2Ebin%d",i),
                                                 Form("Invariant mass of 2 local maxima cells #sigma_{#phi #phi}^{2}, %s",sEBin[i].Data()),
                                                 ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispPhiNLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispPhiNLocMax2Ebin[i]->SetXTitle("#sigma_{#phi #phi}^{2}");
        outputContainer->Add(fhMassDispPhiNLocMax2Ebin[i]) ;   
        
        fhMassDispPhiNLocMaxNEbin[i]  = new TH2F(Form("hMassDispPhiNLocMaxNEbin%d",i),
                                                 Form("Invariant mass of N>2 local maxima cells vs #sigma_{#phi #phi}^{2}, %s",sEBin[i].Data()),
                                                 ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispPhiNLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispPhiNLocMaxNEbin[i]->SetXTitle("#sigma_{#phi #phi}^{2}");
        outputContainer->Add(fhMassDispPhiNLocMaxNEbin[i]) ;   
        
        fhMassDispAsyNLocMax1Ebin[i]  = new TH2F(Form("hMassDispAsyNLocMax1Ebin%d",i),
                                                 Form("Invariant mass of 2 highest energy cells A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), %s",sEBin[i].Data()),
                                                 200,-1,1,mbins,mmin,mmax); 
        fhMassDispAsyNLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispAsyNLocMax1Ebin[i]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
        outputContainer->Add(fhMassDispAsyNLocMax1Ebin[i]) ;   
        
        fhMassDispAsyNLocMax2Ebin[i]  = new TH2F(Form("hMassDispAsyNLocMax2Ebin%d",i),
                                                 Form("Invariant mass of 2 local maxima cells A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), %s",sEBin[i].Data()),
                                                 200,-1,1,mbins,mmin,mmax); 
        fhMassDispAsyNLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispAsyNLocMax2Ebin[i]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
        outputContainer->Add(fhMassDispAsyNLocMax2Ebin[i]) ;   
        
        fhMassDispAsyNLocMaxNEbin[i]  = new TH2F(Form("hMassDispAsyNLocMaxNEbin%d",i),
                                                 Form("Invariant mass of N>2 local maxima cells vs A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), %s",sEBin[i].Data()),
                                                 200,-1,1,mbins,mmin,mmax); 
        fhMassDispAsyNLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispAsyNLocMaxNEbin[i]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
        outputContainer->Add(fhMassDispAsyNLocMaxNEbin[i]) ;   
      }
    }
  }
  
  if(IsDataMC() && fFillMCHisto && (asyOn || m02On))
  {
    fhMCGenSplitEFracAfterCutsNLocMax1MCPi0     = new TH2F("hMCGenSplitEFracAfterCutsNLocMax1MCPi0",
                                                           "E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  = 1 MC Pi0, after M02 and Asym cut",
                                                           nptbins,ptmin,ptmax,200,0,2);
    fhMCGenSplitEFracAfterCutsNLocMax1MCPi0   ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
    fhMCGenSplitEFracAfterCutsNLocMax1MCPi0   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCGenSplitEFracAfterCutsNLocMax1MCPi0) ;
    
    fhMCGenSplitEFracAfterCutsNLocMax2MCPi0    = new TH2F("hMCGenSplitEFracAfterCutsNLocMax2MCPi0",
                                                          "E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  = 2 MC Pi0, after M02 and Asym cut",
                                                          nptbins,ptmin,ptmax,200,0,2);
    fhMCGenSplitEFracAfterCutsNLocMax2MCPi0  ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
    fhMCGenSplitEFracAfterCutsNLocMax2MCPi0  ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCGenSplitEFracAfterCutsNLocMax2MCPi0) ;
    
    
    fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0    = new TH2F("hMCGenSplitEFracAfterCutsNLocMaxNMCPi0",
                                                          "E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  > 2 MC Pi0, after M02 and Asym cut",
                                                          nptbins,ptmin,ptmax,200,0,2);
    fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0  ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
    fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0  ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0) ;
    
    fhMCGenFracAfterCutsNLocMax1MCPi0     = new TH2F("hMCGenFracAfterCutsNLocMax1MCPi0",
                                                     "E_{gen} / E_{reco} vs E_{reco} for N max  = 1 MC Pi0, after M02 and Asym cut",
                                                     nptbins,ptmin,ptmax,200,0,2);
    fhMCGenFracAfterCutsNLocMax1MCPi0   ->SetYTitle("E_{gen} / E_{reco}");
    fhMCGenFracAfterCutsNLocMax1MCPi0   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCGenFracAfterCutsNLocMax1MCPi0) ;
    
    fhMCGenFracAfterCutsNLocMax2MCPi0    = new TH2F("hMCGenFracAfterCutsNLocMax2MCPi0",
                                                    " E_{gen} / E_{reco} vs E_{reco} for N max  = 2 MC Pi0, after M02 and Asym cut",
                                                    nptbins,ptmin,ptmax,200,0,2);
    fhMCGenFracAfterCutsNLocMax2MCPi0   ->SetYTitle("E_{gen} / E_{reco}");
    fhMCGenFracAfterCutsNLocMax2MCPi0   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCGenFracAfterCutsNLocMax2MCPi0) ;
    
    
    fhMCGenFracAfterCutsNLocMaxNMCPi0   = new TH2F("hMCGenFracAfterCutsNLocMaxNMCPi0",
                                                   " E_{gen} / E_{reco}  vs E_{reco} for N max  > 2 MC Pi0, after M02 and Asym cut",
                                                   nptbins,ptmin,ptmax,200,0,2);
    fhMCGenFracAfterCutsNLocMaxNMCPi0   ->SetYTitle("E_{gen} / E_{reco}");
    fhMCGenFracAfterCutsNLocMaxNMCPi0   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCGenFracAfterCutsNLocMaxNMCPi0) ;
    
  }
  
  if(fFillTMResidualHisto && fFillTMHisto)
  {
    for(Int_t i = 0; i < n; i++)
    {  
      
      fhTrackMatchedDEtaNLocMax1[i]  = new TH2F
      (Form("hTrackMatchedDEtaNLocMax1%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, 1 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaNLocMax1[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaNLocMax1[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNLocMax1[i]  = new TH2F
      (Form("hTrackMatchedDPhiNLocMax1%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, 1 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiNLocMax1[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiNLocMax1[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaNLocMax1[i]) ; 
      outputContainer->Add(fhTrackMatchedDPhiNLocMax1[i]) ;
      
      fhTrackMatchedDEtaNLocMax2[i]  = new TH2F
      (Form("hTrackMatchedDEtaNLocMax2%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, 2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaNLocMax2[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaNLocMax2[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNLocMax2[i]  = new TH2F
      (Form("hTrackMatchedDPhiNLocMax2%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, 2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiNLocMax2[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiNLocMax2[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaNLocMax2[i]) ; 
      outputContainer->Add(fhTrackMatchedDPhiNLocMax2[i]) ;
      
      fhTrackMatchedDEtaNLocMaxN[i]  = new TH2F
      (Form("hTrackMatchedDEtaNLocMaxN%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, N>2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaNLocMaxN[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaNLocMaxN[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNLocMaxN[i]  = new TH2F
      (Form("hTrackMatchedDPhiNLocMaxN%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, N>2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiNLocMaxN[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiNLocMaxN[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaNLocMaxN[i]) ; 
      outputContainer->Add(fhTrackMatchedDPhiNLocMaxN[i]) ;
      
      fhTrackMatchedDEtaNLocMax1Pos[i]  = new TH2F
      (Form("hTrackMatchedDEtaNLocMax1Pos%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, 1 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEtaNLocMax1Pos[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaNLocMax1Pos[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNLocMax1Pos[i]  = new TH2F
      (Form("hTrackMatchedDPhiNLocMax1Pos%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, 1 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhiNLocMax1Pos[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiNLocMax1Pos[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaNLocMax1Pos[i]) ;
      outputContainer->Add(fhTrackMatchedDPhiNLocMax1Pos[i]) ;
      
      fhTrackMatchedDEtaNLocMax2Pos[i]  = new TH2F
      (Form("hTrackMatchedDEtaNLocMax2Pos%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, 2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEtaNLocMax2Pos[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaNLocMax2Pos[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNLocMax2Pos[i]  = new TH2F
      (Form("hTrackMatchedDPhiNLocMax2Pos%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, 2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhiNLocMax2Pos[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiNLocMax2Pos[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaNLocMax2Pos[i]) ;
      outputContainer->Add(fhTrackMatchedDPhiNLocMax2Pos[i]) ;
      
      fhTrackMatchedDEtaNLocMaxNPos[i]  = new TH2F
      (Form("hTrackMatchedDEtaNLocMaxNPos%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, N>2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEtaNLocMaxNPos[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaNLocMaxNPos[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNLocMaxNPos[i]  = new TH2F
      (Form("hTrackMatchedDPhiNLocMaxNPos%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, N>2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhiNLocMaxNPos[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiNLocMaxNPos[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaNLocMaxNPos[i]) ;
      outputContainer->Add(fhTrackMatchedDPhiNLocMaxNPos[i]) ;
      
      fhTrackMatchedDEtaNLocMax1Neg[i]  = new TH2F
      (Form("hTrackMatchedDEtaNLocMax1Neg%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, 1 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEtaNLocMax1Neg[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaNLocMax1Neg[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNLocMax1Neg[i]  = new TH2F
      (Form("hTrackMatchedDPhiNLocMax1Neg%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, 1 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhiNLocMax1Neg[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiNLocMax1Neg[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaNLocMax1Neg[i]) ;
      outputContainer->Add(fhTrackMatchedDPhiNLocMax1Neg[i]) ;
      
      fhTrackMatchedDEtaNLocMax2Neg[i]  = new TH2F
      (Form("hTrackMatchedDEtaNLocMax2Neg%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, 2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEtaNLocMax2Neg[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaNLocMax2Neg[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNLocMax2Neg[i]  = new TH2F
      (Form("hTrackMatchedDPhiNLocMax2Neg%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, 2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhiNLocMax2Neg[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiNLocMax2Neg[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaNLocMax2Neg[i]) ;
      outputContainer->Add(fhTrackMatchedDPhiNLocMax2Neg[i]) ;
      
      fhTrackMatchedDEtaNLocMaxNNeg[i]  = new TH2F
      (Form("hTrackMatchedDEtaNLocMaxNNeg%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, N>2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEtaNLocMaxNNeg[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaNLocMaxNNeg[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNLocMaxNNeg[i]  = new TH2F
      (Form("hTrackMatchedDPhiNLocMaxNNeg%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, N>2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhiNLocMaxNNeg[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiNLocMaxNNeg[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaNLocMaxNNeg[i]) ;
      outputContainer->Add(fhTrackMatchedDPhiNLocMaxNNeg[i]) ;
      
    }
  }
  
  if(fFillAngleHisto)
  {
    for(Int_t i = 0; i < n; i++)
    {
      for(Int_t j = 0; j < nMatched; j++)
      {
        
        fhAnglePairNLocMax1[i][j]  = new TH2F(Form("hAnglePairNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Opening angle split sub-clusters of cluster NLM=1 vs pair Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                           nptbins,ptmin,ptmax,200,0,0.2);
        fhAnglePairNLocMax1[i][j]->SetYTitle("#alpha (rad)");
        fhAnglePairNLocMax1[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhAnglePairNLocMax1[i][j]) ;
        
        fhAnglePairNLocMax2[i][j]  = new TH2F(Form("hAnglePairNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Opening angle split sub-clusters of cluster NLM=2 cells vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                           nptbins,ptmin,ptmax,200,0,0.2);
        fhAnglePairNLocMax2[i][j]->SetYTitle("#alpha (rad)");
        fhAnglePairNLocMax2[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhAnglePairNLocMax2[i][j]) ;
        
        fhAnglePairNLocMaxN[i][j]  = new TH2F(Form("hAnglePairNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Opening angle split sub-clusters of cluster NLM>2 vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                           nptbins,ptmin,ptmax,200,0,0.2);
        fhAnglePairNLocMaxN[i][j]->SetYTitle("#alpha (rad)");
        fhAnglePairNLocMaxN[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhAnglePairNLocMaxN[i][j]) ;
        
        if(asyOn || m02On)
        {
          fhAnglePairAfterCutsNLocMax1[i][j]  = new TH2F(Form("hAnglePairAfterCutsNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Opening angle split sub-clusters of cluster NLM=1, after cuts, vs pair Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                nptbins,ptmin,ptmax,200,0,0.2);
          fhAnglePairAfterCutsNLocMax1[i][j]->SetYTitle("#alpha (rad)");
          fhAnglePairAfterCutsNLocMax1[i][j]->SetXTitle("E (GeV)");
          outputContainer->Add(fhAnglePairAfterCutsNLocMax1[i][j]) ;
          
          fhAnglePairAfterCutsNLocMax2[i][j]  = new TH2F(Form("hAnglePairAfterCutsNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Opening angle split sub-clusters of cluster, after cuts, NLM=2 cells vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                nptbins,ptmin,ptmax,200,0,0.2);
          fhAnglePairAfterCutsNLocMax2[i][j]->SetYTitle("#alpha (rad)");
          fhAnglePairAfterCutsNLocMax2[i][j]->SetXTitle("E (GeV)");
          outputContainer->Add(fhAnglePairAfterCutsNLocMax2[i][j]) ;
          
          fhAnglePairAfterCutsNLocMaxN[i][j]  = new TH2F(Form("hAnglePairAfterCutsNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Opening angle split sub-clusters of cluster, after cuts, NLM>2 vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                nptbins,ptmin,ptmax,200,0,0.2);
          fhAnglePairAfterCutsNLocMaxN[i][j]->SetYTitle("#alpha (rad)");
          fhAnglePairAfterCutsNLocMaxN[i][j]->SetXTitle("E (GeV)");
          outputContainer->Add(fhAnglePairAfterCutsNLocMaxN[i][j]) ;

        }
        
        fhAnglePairPi0NLocMax1[i][j]  = new TH2F(Form("hAnglePairPi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("Opening angle split sub-clusters of cluster, Pi0 ID, NLM=1 vs pair Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                              nptbins,ptmin,ptmax,200,0,0.2);
        fhAnglePairPi0NLocMax1[i][j]->SetYTitle("#alpha (rad)");
        fhAnglePairPi0NLocMax1[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhAnglePairPi0NLocMax1[i][j]) ;
        
        fhAnglePairPi0NLocMax2[i][j]  = new TH2F(Form("hAnglePairPi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("Opening angle split sub-clusters of cluster, Pi0 ID, NLM=2 cells vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                              nptbins,ptmin,ptmax,200,0,0.2);
        fhAnglePairPi0NLocMax2[i][j]->SetYTitle("#alpha (rad)");
        fhAnglePairPi0NLocMax2[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhAnglePairPi0NLocMax2[i][j]) ;
        
        fhAnglePairPi0NLocMaxN[i][j]  = new TH2F(Form("hAnglePairPi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("Opening angle split sub-clusters of cluster, Pi0 ID, NLM>2 vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                              nptbins,ptmin,ptmax,200,0,0.2);
        fhAnglePairPi0NLocMaxN[i][j]->SetYTitle("#alpha (rad)");
        fhAnglePairPi0NLocMaxN[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhAnglePairPi0NLocMaxN[i][j]) ;
        
        fhAnglePairMassNLocMax1[i][j]  = new TH2F(Form("hAnglePairMassNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                               Form("Opening angle split sub-clusters of cluster, Pi0 ID, NLM=1 vs Mass for E > 12 GeV, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                               mbins,mmin,mmax,200,0,0.2);
        fhAnglePairMassNLocMax1[i][j]->SetXTitle("M (GeV/c^{2})");
        fhAnglePairMassNLocMax1[i][j]->SetYTitle("#alpha (rad)");
        outputContainer->Add(fhAnglePairMassNLocMax1[i][j]) ;
        
        
        fhAnglePairMassNLocMax2[i][j]  = new TH2F(Form("hAnglePairMassNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                               Form("Opening angle split sub-clusters of cluster NLM=2 vs Mass for E > 12 GeV, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                               mbins,mmin,mmax,200,0,0.2);
        fhAnglePairMassNLocMax2[i][j]->SetXTitle("M (GeV/c^{2})");
        fhAnglePairMassNLocMax2[i][j]->SetYTitle("#alpha (rad)");
        outputContainer->Add(fhAnglePairMassNLocMax2[i][j]) ;
        
        fhAnglePairMassNLocMaxN[i][j]  = new TH2F(Form("hAnglePairMassNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                               Form("Opening angle split sub-clusters of cluster NLM>2 vs Mass for E > 12 GeV, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                               mbins,mmin,mmax,200,0,0.2);
        fhAnglePairMassNLocMaxN[i][j]->SetXTitle("M (GeV/c^{2})");
        fhAnglePairMassNLocMaxN[i][j]->SetYTitle("#alpha (rad)");
        outputContainer->Add(fhAnglePairMassNLocMaxN[i][j]) ;  

        
        fhAnglePairM02NLocMax1[i][j]  = new TH2F(Form("hAnglePairM02NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                  Form("Opening angle split sub-clusters of cluster, Pi0 ID, NLM=1 vs Mass for E > 15 GeV, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                  ssbins,ssmin,ssmax,200,0,0.2);
        fhAnglePairM02NLocMax1[i][j]->SetXTitle("#lambda_{0}^{2}");
        fhAnglePairM02NLocMax1[i][j]->SetYTitle("#alpha (rad)");
        outputContainer->Add(fhAnglePairM02NLocMax1[i][j]) ;
        
        
        fhAnglePairM02NLocMax2[i][j]  = new TH2F(Form("hAnglePairM02NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                  Form("Opening angle split sub-clusters of cluster NLM=2 vs Mass for E > 12 GeV, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                  ssbins,ssmin,ssmax,200,0,0.2);
        fhAnglePairM02NLocMax2[i][j]->SetXTitle("#lambda_{0}^{2}");
        fhAnglePairM02NLocMax2[i][j]->SetYTitle("#alpha (rad)");
        outputContainer->Add(fhAnglePairM02NLocMax2[i][j]) ;
        
        fhAnglePairM02NLocMaxN[i][j]  = new TH2F(Form("hAnglePairM02NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                  Form("Opening angle split sub-clusters of cluster NLM>2 vs Mass for E > 12 GeV, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                  ssbins,ssmin,ssmax,200,0,0.2);
        fhAnglePairM02NLocMaxN[i][j]->SetXTitle("#lambda_{0}^{2}");
        fhAnglePairM02NLocMaxN[i][j]->SetYTitle("#alpha (rad)");
        outputContainer->Add(fhAnglePairM02NLocMaxN[i][j]) ;

        fhAnglePairOverM02NLocMax1[i][j]  = new TH2F(Form("hAnglePairOverM02NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("Opening angle split sub-clusters  / M02 of cluster NLM=1 vs pair Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                              nptbins,ptmin,ptmax,200,0,0.2);
        fhAnglePairOverM02NLocMax1[i][j]->SetYTitle("#alpha / #lambda_{0}^{2}");
        fhAnglePairOverM02NLocMax1[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhAnglePairOverM02NLocMax1[i][j]) ;
        
        fhAnglePairOverM02NLocMax2[i][j]  = new TH2F(Form("hAnglePairOverM02NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("Opening angle split sub-clusters  / M02 of cluster NLM=2 cells vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                              nptbins,ptmin,ptmax,200,0,0.2);
        fhAnglePairOverM02NLocMax2[i][j]->SetYTitle("#alpha / #lambda_{0}^{2}");
        fhAnglePairOverM02NLocMax2[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhAnglePairOverM02NLocMax2[i][j]) ;
        
        fhAnglePairOverM02NLocMaxN[i][j]  = new TH2F(Form("hAnglePairOverM02NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("Opening angle split sub-clusters  / M02  of cluster NLM>2 vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                              nptbins,ptmin,ptmax,200,0,0.2);
        fhAnglePairOverM02NLocMaxN[i][j]->SetYTitle("#alpha / #lambda_{0}^{2}");
        fhAnglePairOverM02NLocMaxN[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhAnglePairOverM02NLocMaxN[i][j]) ;

        
        fhAnglePairOverM02NLocMax1Overlap0[i][j]  = new TH2F(Form("hAnglePairOverM02NLocMax1Overlap0%s%s",pname[i].Data(),sMatched[j].Data()),
                                                     Form("Opening angle split sub-clusters  / M02  of cluster NLM=1 vs pair Energy, no overlap, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                     nptbins,ptmin,ptmax,200,0,0.2);
        fhAnglePairOverM02NLocMax1Overlap0[i][j]->SetYTitle("#alpha / #lambda_{0}^{2}");
        fhAnglePairOverM02NLocMax1Overlap0[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhAnglePairOverM02NLocMax1Overlap0[i][j]) ;
        
        fhAnglePairOverM02NLocMax2Overlap0[i][j]  = new TH2F(Form("hAnglePairOverM02NLocMax2Overlap0%s%s",pname[i].Data(),sMatched[j].Data()),
                                                     Form("Opening angle split sub-clusters / M02 of cluster NLM=2 cells vs Energy, no overlap, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                     nptbins,ptmin,ptmax,200,0,0.2);
        fhAnglePairOverM02NLocMax2Overlap0[i][j]->SetYTitle("#alpha / #lambda_{0}^{2}");
        fhAnglePairOverM02NLocMax2Overlap0[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhAnglePairOverM02NLocMax2Overlap0[i][j]) ;
        
        fhAnglePairOverM02NLocMaxNOverlap0[i][j]  = new TH2F(Form("hAnglePairOverM02NLocMaxNOverlap0%s%s",pname[i].Data(),sMatched[j].Data()),
                                                     Form("Opening angle split sub-clusters  / M02  of cluster NLM>2 vs Energy, no overlap, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                     nptbins,ptmin,ptmax,200,0,0.2);
        fhAnglePairOverM02NLocMaxNOverlap0[i][j]->SetYTitle("#alpha / #lambda_{0}^{2}");
        fhAnglePairOverM02NLocMaxNOverlap0[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhAnglePairOverM02NLocMaxNOverlap0[i][j]) ;

        
      }
    }
    
    if(IsDataMC())
    {
      fhAnglePairPrimPi0RecoNLocMax1  = new TH2F("fhAnglePairPrimPi0RecoNLocMax1",
                                                 "Opening angle split neutral sub-clusters reconstructed / generated #pi^{0} vs pair Energy, NLM=1",
                                                 nptbins,ptmin,ptmax,300,0,3);
      fhAnglePairPrimPi0RecoNLocMax1->SetYTitle("#alpha_{reco} / #alpha_{gen}");
      fhAnglePairPrimPi0RecoNLocMax1->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairPrimPi0RecoNLocMax1) ;

      fhAnglePairPrimPi0RecoNLocMax2  = new TH2F("fhAnglePairPrimPi0RecoNLocMax2",
                                                 "Opening angle split neutral sub-clusters reconstructed / generated #pi^{0} vs pair Energy, NLM=2",
                                                 nptbins,ptmin,ptmax,300,0,3);
      fhAnglePairPrimPi0RecoNLocMax2->SetYTitle("#alpha_{reco} / #alpha_{gen}");
      fhAnglePairPrimPi0RecoNLocMax2->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairPrimPi0RecoNLocMax2) ;

      fhAnglePairPrimPi0RecoNLocMaxN  = new TH2F("fhAnglePairPrimPi0RecoNLocMaxN",
                                                 "Opening angle split neutral sub-clusters reconstructed / generated #pi^{0} vs pair Energy, NLM>2",
                                                 nptbins,ptmin,ptmax,300,0,3);
      fhAnglePairPrimPi0RecoNLocMaxN->SetYTitle("#alpha_{reco} / #alpha_{gen}");
      fhAnglePairPrimPi0RecoNLocMaxN->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairPrimPi0RecoNLocMaxN) ;

      
      fhAnglePairPrimPi0vsRecoNLocMax1  = new TH2F("fhAnglePairPrimPi0vsRecoNLocMax1",
                                                   "Opening angle split neutral sub-clusters reconstructed vs generated #pi^{0} for E > 15 GeV, NLM=1",
                                                   200,0,0.2,200,0,0.2);
      fhAnglePairPrimPi0vsRecoNLocMax1->SetYTitle("#alpha_{reco} (rad)");
      fhAnglePairPrimPi0vsRecoNLocMax1->SetXTitle("#alpha_{gen} (rad)");
      outputContainer->Add(fhAnglePairPrimPi0vsRecoNLocMax1) ;

      fhAnglePairPrimPi0vsRecoNLocMax2  = new TH2F("fhAnglePairPrimPi0vsRecoNLocMax2",
                                                   "Opening angle split neutral sub-clusters reconstructed vs generated #pi^{0} for E > 10 GeV, NLM=2",
                                                   200,0,0.2,200,0,0.2);
      fhAnglePairPrimPi0vsRecoNLocMax2->SetYTitle("#alpha_{reco} (rad)");
      fhAnglePairPrimPi0vsRecoNLocMax2->SetXTitle("#alpha_{gen} (rad)");
      outputContainer->Add(fhAnglePairPrimPi0vsRecoNLocMax2) ;

      fhAnglePairPrimPi0vsRecoNLocMaxN  = new TH2F("fhAnglePairPrimPi0vsRecoNLocMaxN",
                                                   "Opening angle split neutral sub-clusters reconstructed vs generated #pi^{0} for E > 10 GeV, NLM=2",
                                                   200,0,0.2,200,0,0.2);
      fhAnglePairPrimPi0vsRecoNLocMaxN->SetYTitle("#alpha_{reco} (rad)");
      fhAnglePairPrimPi0vsRecoNLocMaxN->SetXTitle("#alpha_{gen} (rad)");
      outputContainer->Add(fhAnglePairPrimPi0vsRecoNLocMaxN) ;
      
      
      fhAnglePairPrimPi0OverM02NLocMax1  = new TH2F("fhAnglePairPrimPi0OverM02NLocMax1",
                                                 "Primary Opening angle split neutral sub-clusters reconstructed / Over vs cluster Energy, NLM=1",
                                                 nptbins,ptmin,ptmax,200,0,0.2);
      fhAnglePairPrimPi0OverM02NLocMax1->SetYTitle("#alpha_{gen} / #lambda_{0}^{2}");
      fhAnglePairPrimPi0OverM02NLocMax1->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairPrimPi0OverM02NLocMax1) ;
      
      fhAnglePairPrimPi0OverM02NLocMax2  = new TH2F("fhAnglePairPrimPi0OverM02NLocMax2",
                                                 "Primary Opening angle split neutral sub-clusters reconstructed / Over vs cluster Energy, NLM=2",
                                                 nptbins,ptmin,ptmax,200,0,0.2);
      fhAnglePairPrimPi0OverM02NLocMax2->SetYTitle("#alpha_{gen} / #lambda_{0}^{2}");
      fhAnglePairPrimPi0OverM02NLocMax2->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairPrimPi0OverM02NLocMax2) ;
      
      fhAnglePairPrimPi0OverM02NLocMaxN  = new TH2F("fhAnglePairPrimPi0OverM02NLocMaxN",
                                                 "Primary Opening angle split neutral sub-clusters reconstructed / M02 vs cluster Energy, NLM>2",
                                                 nptbins,ptmin,ptmax,200,0,0.2);
      fhAnglePairPrimPi0OverM02NLocMaxN->SetYTitle("#alpha_{gen} / #lambda_{0}^{2}");
      fhAnglePairPrimPi0OverM02NLocMaxN->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairPrimPi0OverM02NLocMaxN) ;

    }
  }
 
  // Same as asymmetry ...
  if(fFillThetaStarHisto)
  {
    for(Int_t i = 0; i < n; i++)
    {
      for(Int_t j = 0; j < nMatched; j++)
      {
        
        fhCosThStarNLocMax1[i][j]  = new TH2F(Form("hCosThStarNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("cos(#theta^{*}) split sub-clusters of cluster NLM=1 vs pair Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                              nptbins,ptmin,ptmax,200,-1,1);
        fhCosThStarNLocMax1[i][j]->SetYTitle("cos(#theta^{*})");
        fhCosThStarNLocMax1[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhCosThStarNLocMax1[i][j]) ;
        
        fhCosThStarNLocMax2[i][j]  = new TH2F(Form("hCosThStarNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("cos(#theta^{*}) split sub-clusters of cluster NLM=2 cells vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                              nptbins,ptmin,ptmax,200,-1,1);
        fhCosThStarNLocMax2[i][j]->SetYTitle("cos(#theta^{*})");
        fhCosThStarNLocMax2[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhCosThStarNLocMax2[i][j]) ;
        
        fhCosThStarNLocMaxN[i][j]  = new TH2F(Form("hCosThStarNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                              Form("cos(#theta^{*}) split sub-clusters of cluster NLM>2 vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                              nptbins,ptmin,ptmax,200,-1,1);
        fhCosThStarNLocMaxN[i][j]->SetYTitle("cos(#theta^{*})");
        fhCosThStarNLocMaxN[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhCosThStarNLocMaxN[i][j]) ;
        
        if(asyOn || m02On)
        {
          fhCosThStarAfterCutsNLocMax1[i][j]  = new TH2F(Form("hCosThStarAfterCutsNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                         Form("cos(#theta^{*}) split sub-clusters of cluster NLM=1, after cuts, vs pair Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                         nptbins,ptmin,ptmax,200,-1,1);
          fhCosThStarAfterCutsNLocMax1[i][j]->SetYTitle("cos(#theta^{*})");
          fhCosThStarAfterCutsNLocMax1[i][j]->SetXTitle("E (GeV)");
          outputContainer->Add(fhCosThStarAfterCutsNLocMax1[i][j]) ;
          
          fhCosThStarAfterCutsNLocMax2[i][j]  = new TH2F(Form("hCosThStarAfterCutsNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                         Form("cos(#theta^{*}) split sub-clusters of cluster, after cuts, NLM=2 cells vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                         nptbins,ptmin,ptmax,200,-1,1);
          fhCosThStarAfterCutsNLocMax2[i][j]->SetYTitle("cos(#theta^{*})");
          fhCosThStarAfterCutsNLocMax2[i][j]->SetXTitle("E (GeV)");
          outputContainer->Add(fhCosThStarAfterCutsNLocMax2[i][j]) ;
          
          fhCosThStarAfterCutsNLocMaxN[i][j]  = new TH2F(Form("hCosThStarAfterCutsNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                         Form("cos(#theta^{*}) split sub-clusters of cluster, after cuts, NLM>2 vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                         nptbins,ptmin,ptmax,200,-1,1);
          fhCosThStarAfterCutsNLocMaxN[i][j]->SetYTitle("cos(#theta^{*})");
          fhCosThStarAfterCutsNLocMaxN[i][j]->SetXTitle("E (GeV)");
          outputContainer->Add(fhCosThStarAfterCutsNLocMaxN[i][j]) ;
          
        }
        
        fhCosThStarPi0NLocMax1[i][j]  = new TH2F(Form("hCosThStarPi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                 Form("cos(#theta^{*}) split sub-clusters of cluster, Pi0 ID, NLM=1 vs pair Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                 nptbins,ptmin,ptmax,200,-1,1);
        fhCosThStarPi0NLocMax1[i][j]->SetYTitle("cos(#theta^{*})");
        fhCosThStarPi0NLocMax1[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhCosThStarPi0NLocMax1[i][j]) ;
        
        fhCosThStarPi0NLocMax2[i][j]  = new TH2F(Form("hCosThStarPi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                 Form("cos(#theta^{*}) split sub-clusters of cluster, Pi0 ID, NLM=2 cells vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                 nptbins,ptmin,ptmax,200,-1,1);
        fhCosThStarPi0NLocMax2[i][j]->SetYTitle("cos(#theta^{*})");
        fhCosThStarPi0NLocMax2[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhCosThStarPi0NLocMax2[i][j]) ;
        
        fhCosThStarPi0NLocMaxN[i][j]  = new TH2F(Form("hCosThStarPi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                 Form("cos(#theta^{*}) split sub-clusters of cluster, Pi0 ID, NLM>2 vs Energy, %s, %s",ptype[i].Data(),sMatched[j].Data()),
                                                 nptbins,ptmin,ptmax,200,-1,1);
        fhCosThStarPi0NLocMaxN[i][j]->SetYTitle("cos(#theta^{*})");
        fhCosThStarPi0NLocMaxN[i][j]->SetXTitle("E (GeV)");
        outputContainer->Add(fhCosThStarPi0NLocMaxN[i][j]) ;
        
      }
    }
  }

  
  for(Int_t j = 0; j < nMatched; j++)
  {
    fhSplitEFractionvsAsyNLocMax1[j]     = new TH2F(Form("hSplitEFractionvsAsyNLocMax1%s",sMatched[j].Data()),
                                                    Form("(E1+E2)/E_{cluster} vs (E_{split1}-E_{split2})/(E_{split1}+E_{split2}) for N max  = 1, E>12, %s",sMatched[j].Data()),
                                                    100,-1,1,120,0,1.2); 
    fhSplitEFractionvsAsyNLocMax1[j]   ->SetXTitle("(E_{split1}-E_{split2})/(E_{split1}+E_{split2})");
    fhSplitEFractionvsAsyNLocMax1[j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
    outputContainer->Add(fhSplitEFractionvsAsyNLocMax1[j]) ; 
    
    fhSplitEFractionvsAsyNLocMax2[j]     = new TH2F(Form("hSplitEFractionvsAsyNLocMax2%s",sMatched[j].Data()),
                                                    Form("(E1+E2)/E_{cluster} vs (E_{split1}-E_{split2})/(E_{split1}+E_{split2}) for N max  = 2,E>12, %s",sMatched[j].Data()),
                                                    100,-1,1,120,0,1.2); 
    fhSplitEFractionvsAsyNLocMax2[j]   ->SetXTitle("(E_{split1}-E_{split2})/(E_{split1}+E_{split2})");
    fhSplitEFractionvsAsyNLocMax2[j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
    outputContainer->Add(fhSplitEFractionvsAsyNLocMax2[j]) ; 
    
    fhSplitEFractionvsAsyNLocMaxN[j]    = new TH2F(Form("hSplitEFractionvsAsyNLocMaxN%s",sMatched[j].Data()),
                                                   Form("(E1+E2)/E_{cluster} vs (E_{split1}-E_{split2})/(E_{split1}+E_{split2}) for N max  > 2, E>12, %s",sMatched[j].Data()),
                                                   100,-1,1,120,0,1.2); 
    fhSplitEFractionvsAsyNLocMaxN[j]   ->SetXTitle("(E_{split1}-E_{split2})/(E_{split1}+E_{split2})");
    fhSplitEFractionvsAsyNLocMaxN[j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
    outputContainer->Add(fhSplitEFractionvsAsyNLocMaxN[j]) ; 
  }
  
  
  fhClusterEtaPhiNLocMax1  = new TH2F
  ("hClusterEtaPhiNLocMax1","Neutral Clusters with E > 8 GeV, NLM = 1: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhClusterEtaPhiNLocMax1->SetYTitle("#phi (rad)");
  fhClusterEtaPhiNLocMax1->SetXTitle("#eta");
  outputContainer->Add(fhClusterEtaPhiNLocMax1) ;
  
  fhClusterEtaPhiNLocMax2  = new TH2F
  ("hClusterEtaPhiNLocMax2","Neutral Clusters with E > 8 GeV, NLM = 2: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhClusterEtaPhiNLocMax2->SetYTitle("#phi (rad)");
  fhClusterEtaPhiNLocMax2->SetXTitle("#eta");
  outputContainer->Add(fhClusterEtaPhiNLocMax2) ;
  
  fhClusterEtaPhiNLocMaxN  = new TH2F
  ("hClusterEtaPhiNLocMaxN","Neutral Clusters with E > 8 GeV, NLM > 2: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhClusterEtaPhiNLocMaxN->SetYTitle("#phi (rad)");
  fhClusterEtaPhiNLocMaxN->SetXTitle("#eta");
  outputContainer->Add(fhClusterEtaPhiNLocMaxN) ;
  
  fhPi0EtaPhiNLocMax1  = new TH2F
  ("hPi0EtaPhiNLocMax1","Selected #pi^{0}'s with E > 8 GeV, NLM = 1: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhPi0EtaPhiNLocMax1->SetYTitle("#phi (rad)");
  fhPi0EtaPhiNLocMax1->SetXTitle("#eta");
  outputContainer->Add(fhPi0EtaPhiNLocMax1) ;
  
  fhPi0EtaPhiNLocMax2  = new TH2F
  ("hPi0EtaPhiNLocMax2","Selected #pi^{0}'s with E > 8 GeV, NLM = 2: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhPi0EtaPhiNLocMax2->SetYTitle("#phi (rad)");
  fhPi0EtaPhiNLocMax2->SetXTitle("#eta");
  outputContainer->Add(fhPi0EtaPhiNLocMax2) ;
  
  fhPi0EtaPhiNLocMaxN  = new TH2F
  ("hPi0EtaPhiNLocMaxN","Selected #pi^{0}'s with E > 8 GeV, NLM > 2: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhPi0EtaPhiNLocMaxN->SetYTitle("#phi (rad)");
  fhPi0EtaPhiNLocMaxN->SetXTitle("#eta");
  outputContainer->Add(fhPi0EtaPhiNLocMaxN) ;
  
  if(fFillIdEtaHisto)
  {
    fhEtaEtaPhiNLocMax1  = new TH2F
    ("hEtaEtaPhiNLocMax1","Selected #eta's with E > 8 GeV, NLM = 1: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
    fhEtaEtaPhiNLocMax1->SetYTitle("#phi (rad)");
    fhEtaEtaPhiNLocMax1->SetXTitle("#eta");
    outputContainer->Add(fhEtaEtaPhiNLocMax1) ;
    
    fhEtaEtaPhiNLocMax2  = new TH2F
    ("hEtaEtaPhiNLocMax2","Selected #eta's with E > 8 GeV, NLM = 2: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
    fhEtaEtaPhiNLocMax2->SetYTitle("#phi (rad)");
    fhEtaEtaPhiNLocMax2->SetXTitle("#eta");
    outputContainer->Add(fhEtaEtaPhiNLocMax2) ;
    
    fhEtaEtaPhiNLocMaxN  = new TH2F
    ("hEtaEtaPhiNLocMaxN","Selected #eta's with E > 8 GeV, NLM > 2: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
    fhEtaEtaPhiNLocMaxN->SetYTitle("#phi (rad)");
    fhEtaEtaPhiNLocMaxN->SetXTitle("#eta");
    outputContainer->Add(fhEtaEtaPhiNLocMaxN) ;
  }
  
  if(fFillSSWeightHisto)
  {
    for(Int_t nlm = 0; nlm < 3; nlm++)
    {
      fhPi0CellE[nlm]  = new TH2F(Form("hPi0CellENLocMax%s",snlm[nlm].Data()),
                                  Form("Selected #pi^{0}'s, NLM = %s: cluster E vs cell E",snlm[nlm].Data()),
                                  nptbins,ptmin,ptmax, nptbins,ptmin,ptmax);
      fhPi0CellE[nlm]->SetYTitle("E_{cell}");
      fhPi0CellE[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellE[nlm]) ;
      
      fhPi0CellEFrac[nlm]  = new TH2F(Form("hPi0CellEFracNLocMax%s",snlm[nlm].Data()),
                                      Form("Selected #pi^{0}'s, NLM = %s: cluster E vs cell E / cluster E",snlm[nlm].Data()),
                                      nptbins,ptmin,ptmax, 100,0,1);
      fhPi0CellEFrac[nlm]->SetYTitle("E_{cell} / E_{cluster}");
      fhPi0CellEFrac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellEFrac[nlm]) ;
      
      fhPi0CellLogEFrac[nlm]  = new TH2F(Form("hPi0CellLogEFracNLocMax%s",snlm[nlm].Data()),
                                         Form("Selected #pi^{0}'s, NLM = %s: cluster E vs Log(cell E / cluster E)",snlm[nlm].Data()),
                                         nptbins,ptmin,ptmax, 100,-10,0);
      fhPi0CellLogEFrac[nlm]->SetYTitle("Log(E_{cell} / E_{cluster})");
      fhPi0CellLogEFrac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellLogEFrac[nlm]) ;
      
      
      fhPi0CellEMaxEMax2Frac[nlm]  = new TH2F(Form("hPi0CellEMaxEMax2FracNLocMax%s",snlm[nlm].Data()),
                                              Form("Selected #pi^{0}'s, NLM = %s: cluster E vs 2nd loc. max. E / 1st loc. max.  E",snlm[nlm].Data()),
                                              nptbins,ptmin,ptmax, 100,0,1);
      fhPi0CellEMaxEMax2Frac[nlm]->SetYTitle("E_{Loc Max 2} / E_{Loc Max 1}");
      fhPi0CellEMaxEMax2Frac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellEMaxEMax2Frac[nlm]) ;
      
      fhPi0CellEMaxClusterFrac[nlm]  = new TH2F(Form("hPi0CellEMaxClusterFracNLocMax%s",snlm[nlm].Data()),
                                                Form("Selected #pi^{0}'s, NLM = %s: cluster E vs 1st loc. max. E / E cluster",snlm[nlm].Data()),
                                                nptbins,ptmin,ptmax, 100,0,1);
      fhPi0CellEMaxClusterFrac[nlm]->SetYTitle("E_{Loc Max 1} / E_{cluster}");
      fhPi0CellEMaxClusterFrac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellEMaxClusterFrac[nlm]) ;
      
      fhPi0CellEMax2ClusterFrac[nlm]  = new TH2F(Form("hPi0CellEMax2ClusterFracNLocMax%s",snlm[nlm].Data()),
                                                 Form("Selected #pi^{0}'s, NLM = %s: cluster E vs 2nd loc. max. E / E cluster",snlm[nlm].Data()),
                                                 nptbins,ptmin,ptmax, 100,0,1);
      fhPi0CellEMax2ClusterFrac[nlm]->SetYTitle("E_{Loc Max 2} / E_{cluster}");
      fhPi0CellEMax2ClusterFrac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellEMax2ClusterFrac[nlm]) ;
      
      fhPi0CellEMaxFrac[nlm]  = new TH2F(Form("hPi0CellEMaxFracNLocMax%s",snlm[nlm].Data()),
                                         Form("Selected #pi^{0}'s, NLM = %s: cluster E vs 1st loc. max. E / E cell i",snlm[nlm].Data()),
                                         nptbins,ptmin,ptmax, 100,0,1);
      fhPi0CellEMaxFrac[nlm]->SetYTitle("E_{Loc Max 1} / E_{cell i}");
      fhPi0CellEMaxFrac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellEMaxFrac[nlm]) ;
      
      fhPi0CellEMax2Frac[nlm]  = new TH2F(Form("hPi0CellEMax2FracNLocMax%s",snlm[nlm].Data()),
                                          Form("Selected #pi^{0}'s, NLM = %s: cluster E vs 2nd loc. max. E / E cell i",snlm[nlm].Data()),
                                          nptbins,ptmin,ptmax, 200,0,2);
      fhPi0CellEMax2Frac[nlm]->SetYTitle("E_{Loc Max 2} / E_{cell i}");
      fhPi0CellEMax2Frac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellEMax2Frac[nlm]) ;
      
      
      for(Int_t i = 0; i < fSSWeightN; i++)
      {
        fhM02WeightPi0[nlm][i]     = new TH2F(Form("hM02Pi0NLocMax%s_W%d",snlm[nlm].Data(),i),
                                              Form("#lambda_{0}^{2} vs E, with W0 = %2.2f, for NLM = %s", fSSWeight[i], snlm[nlm].Data()),
                                              nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhM02WeightPi0[nlm][i]   ->SetYTitle("#lambda_{0}^{2}");
        fhM02WeightPi0[nlm][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhM02WeightPi0[nlm][i]) ;
      }
      
      for(Int_t i = 0; i < fSSECellCutN; i++)
      {
        fhM02ECellCutPi0[nlm][i]     = new TH2F(Form("hM02Pi0NLocMax%s_Ecell%d",snlm[nlm].Data(),i),
                                                Form("#lambda_{0}^{2} vs E, with Ecell > %2.2f, for NLM = %s", fSSECellCut[i], snlm[nlm].Data()),
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhM02ECellCutPi0[nlm][i]   ->SetYTitle("#lambda_{0}^{2}");
        fhM02ECellCutPi0[nlm][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhM02ECellCutPi0[nlm][i]) ;
      }
      
    }
  }
  
  Int_t tdbins   = GetHistogramRanges()->GetHistoDiffTimeBins() ;    Float_t tdmax  = GetHistogramRanges()->GetHistoDiffTimeMax();     Float_t tdmin  = GetHistogramRanges()->GetHistoDiffTimeMin();
  
  fhPi0EPairDiffTimeNLM1 = new TH2F("hPi0EPairDiffTimeNLocMax1","cluster pair time difference vs E, selected #pi, NLM=1",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
  fhPi0EPairDiffTimeNLM1->SetXTitle("E_{pair} (GeV)");
  fhPi0EPairDiffTimeNLM1->SetYTitle("#Delta t (ns)");
  outputContainer->Add(fhPi0EPairDiffTimeNLM1);
  
  fhPi0EPairDiffTimeNLM2 = new TH2F("hPi0EPairDiffTimeNLocMax2","cluster pair time difference vs E, selected #pi, NLM=2",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
  fhPi0EPairDiffTimeNLM2->SetXTitle("E_{pair} (GeV)");
  fhPi0EPairDiffTimeNLM2->SetYTitle("#Delta t (ns)");
  outputContainer->Add(fhPi0EPairDiffTimeNLM2);
  
  fhPi0EPairDiffTimeNLMN = new TH2F("hPi0EPairDiffTimeNLocMaxN","cluster pair time difference vs E, selected #pi, NLM>2",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
  fhPi0EPairDiffTimeNLMN->SetXTitle("E_{pair} (GeV)");
  fhPi0EPairDiffTimeNLMN->SetYTitle("#Delta t (ns)");
  outputContainer->Add(fhPi0EPairDiffTimeNLMN);
  
  if(fFillIdEtaHisto)
  {
    fhEtaEPairDiffTimeNLM1 = new TH2F("hEtaEPairDiffTimeNLocMax1","cluster pair time difference vs E, selected #eta, NLM=1",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
    fhEtaEPairDiffTimeNLM1->SetXTitle("E_{pair} (GeV)");
    fhEtaEPairDiffTimeNLM1->SetYTitle("#Delta t (ns)");
    outputContainer->Add(fhEtaEPairDiffTimeNLM1);
    
    fhEtaEPairDiffTimeNLM2 = new TH2F("hEtaEPairDiffTimeNLocMax2","cluster pair time difference vs E, selected #eta, NLM=2",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
    fhEtaEPairDiffTimeNLM2->SetXTitle("E_{pair} (GeV)");
    fhEtaEPairDiffTimeNLM2->SetYTitle("#Delta t (ns)");
    outputContainer->Add(fhEtaEPairDiffTimeNLM2);
    
    fhEtaEPairDiffTimeNLMN = new TH2F("hEtaEPairDiffTimeNLocMaxN","cluster pair time difference vs E, selected #eta, NLM>2",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
    fhEtaEPairDiffTimeNLMN->SetXTitle("E_{pair} (GeV)");
    fhEtaEPairDiffTimeNLMN->SetYTitle("#Delta t (ns)");
    outputContainer->Add(fhEtaEPairDiffTimeNLMN);
  }
  
  if(fFillNCellHisto && IsDataMC())
  {
    
    fhNCellMassEHighNLocMax1MCPi0 = new TH2F("hNCellMassEHighNLocMax1MCPi0","n cells vs mass for MC pi0, high energy, NLM=1",ncbins,ncmin,ncmax,mbins,mmin,mmax);
    fhNCellMassEHighNLocMax1MCPi0->SetYTitle("M (GeV/c^{2})");
    fhNCellMassEHighNLocMax1MCPi0->SetXTitle("n cells");
    outputContainer->Add(fhNCellMassEHighNLocMax1MCPi0) ;
    
    fhNCellMassELowNLocMax1MCPi0 = new TH2F("hNCellMassELowNLocMax1MCPi0","n cells vs mass for MC pi0, low energy, NLM=1",ncbins,ncmin,ncmax,mbins,mmin,mmax);
    fhNCellMassELowNLocMax1MCPi0->SetYTitle("M (GeV/c^{2})");
    fhNCellMassELowNLocMax1MCPi0->SetXTitle("n cells");
    outputContainer->Add(fhNCellMassELowNLocMax1MCPi0) ;
    
    fhNCellM02EHighNLocMax1MCPi0 = new TH2F("hNCellM02EHighNLocMax1MCPi0","n cells vs #lambda_{0}^{2} for MC pi0, high energy, NLM=1",ncbins,ncmin,ncmax,ssbins,ssmin,ssmax);
    fhNCellM02EHighNLocMax1MCPi0->SetYTitle("#lambda_{0}^{2}");
    fhNCellM02EHighNLocMax1MCPi0->SetXTitle("n cells");
    outputContainer->Add(fhNCellM02EHighNLocMax1MCPi0) ;
    
    fhNCellM02ELowNLocMax1MCPi0 = new TH2F("hNCellM02ELowNLocMax1MCPi0","n cells vs #lambda_{0}^{2} for MC pi0, low energy, NLM=1",ncbins,ncmin,ncmax,ssbins,ssmin,ssmax);
    fhNCellM02ELowNLocMax1MCPi0->SetYTitle("#lambda_{0}^{2}");
    fhNCellM02ELowNLocMax1MCPi0->SetXTitle("n cells");
    outputContainer->Add(fhNCellM02ELowNLocMax1MCPi0) ;
    
    fhNCellMassEHighNLocMax2MCPi0 = new TH2F("hNCellMassEHighNLocMax2MCPi0","n cells vs mass for MC pi0, high energy, NLM=2",ncbins,ncmin,ncmax,mbins,mmin,mmax);
    fhNCellMassEHighNLocMax2MCPi0->SetYTitle("M (GeV/c^{2})");
    fhNCellMassEHighNLocMax2MCPi0->SetXTitle("n cells");
    outputContainer->Add(fhNCellMassEHighNLocMax2MCPi0) ;
    
    fhNCellMassELowNLocMax2MCPi0 = new TH2F("hNCellMassELowNLocMax2MCPi0","n cells vs mass for MC pi0, low energy, NLM=2",ncbins,ncmin,ncmax,mbins,mmin,mmax);
    fhNCellMassELowNLocMax2MCPi0->SetYTitle("M (GeV/c^{2})");
    fhNCellMassELowNLocMax2MCPi0->SetXTitle("n cells");
    outputContainer->Add(fhNCellMassELowNLocMax2MCPi0) ;
    
    fhNCellM02EHighNLocMax2MCPi0 = new TH2F("hNCellM02EHighNLocMax2MCPi0","n cells vs #lambda_{0}^{2} for MC pi0, high energy, NLM=2",ncbins,ncmin,ncmax,ssbins,ssmin,ssmax);
    fhNCellM02EHighNLocMax2MCPi0->SetYTitle("#lambda_{0}^{2}");
    fhNCellM02EHighNLocMax2MCPi0->SetXTitle("n cells");
    outputContainer->Add(fhNCellM02EHighNLocMax2MCPi0) ;
    
    fhNCellM02ELowNLocMax2MCPi0 = new TH2F("hNCellM02ELowNLocMax2MCPi0","n cells vs #lambda_{0}^{2} for MC pi0, low energy, NLM=2",ncbins,ncmin,ncmax,ssbins,ssmin,ssmax);
    fhNCellM02ELowNLocMax2MCPi0->SetYTitle("#lambda_{0}^{2}");
    fhNCellM02ELowNLocMax2MCPi0->SetXTitle("n cells");
    outputContainer->Add(fhNCellM02ELowNLocMax2MCPi0) ;
    
    fhNCellMassEHighNLocMaxNMCPi0 = new TH2F("hNCellMassEHighNLocMaxNMCPi0","n cells vs mass for MC pi0, high energy, NLM>2",ncbins,ncmin,ncmax,mbins,mmin,mmax);
    fhNCellMassEHighNLocMaxNMCPi0->SetYTitle("M (GeV/c^{2})");
    fhNCellMassEHighNLocMaxNMCPi0->SetXTitle("n cells");
    outputContainer->Add(fhNCellMassEHighNLocMaxNMCPi0) ;
    
    fhNCellMassELowNLocMaxNMCPi0 = new TH2F("hNCellMassELowNLocMaxNMCPi0","n cells vs mass for MC pi0, low energy, NLM>2",ncbins,ncmin,ncmax,mbins,mmin,mmax);
    fhNCellMassELowNLocMaxNMCPi0->SetYTitle("M (GeV/c^{2})");
    fhNCellMassELowNLocMaxNMCPi0->SetXTitle("n cells");
    outputContainer->Add(fhNCellMassELowNLocMaxNMCPi0) ;
    
    fhNCellM02EHighNLocMaxNMCPi0 = new TH2F("hNCellM02EHighNLocMaxNMCPi0","n cells vs #lambda_{0}^{2} for MC pi0, high energy, NLM>2",ncbins,ncmin,ncmax,ssbins,ssmin,ssmax);
    fhNCellM02EHighNLocMaxNMCPi0->SetYTitle("#lambda_{0}^{2}");
    fhNCellM02EHighNLocMaxNMCPi0->SetXTitle("n cells");
    outputContainer->Add(fhNCellM02EHighNLocMaxNMCPi0) ;
    
    fhNCellM02ELowNLocMaxNMCPi0 = new TH2F("hNCellM02ELowNLocMaxNMCPi0","n cells vs #lambda_{0}^{2} for MC pi0, low energy, NLM>2",ncbins,ncmin,ncmax,ssbins,ssmin,ssmax);
    fhNCellM02ELowNLocMaxNMCPi0->SetYTitle("#lambda_{0}^{2}");
    fhNCellM02ELowNLocMaxNMCPi0->SetXTitle("n cells");
    outputContainer->Add(fhNCellM02ELowNLocMaxNMCPi0) ;
    
  }
  
  if(IsDataMC() && fFillMCOverlapHisto)
  {
    for(Int_t i = 1; i < n; i++)
    {
      for(Int_t j = 0; j < 3; j++)
      {
        fhMCENOverlaps[j][i]     = new TH2F(Form("hMCENOverlapsNLocMax%s%s",snlm[j].Data(),pname[i].Data()),
                                            Form("# overlaps vs E for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,10,0,10);
        fhMCENOverlaps[j][i]   ->SetYTitle("# overlaps");
        fhMCENOverlaps[j][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCENOverlaps[j][i]) ;
        
        fhMCEM02Overlap0[j][i]     = new TH2F(Form("hMCEM02Overlap0NLocMax%s%s",snlm[j].Data(),pname[i].Data()),
                                              Form("Overlap 0, #lambda_{0}^{2} vs E for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                              nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhMCEM02Overlap0[j][i]   ->SetYTitle("#lambda_{0}^{2}");
        fhMCEM02Overlap0[j][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCEM02Overlap0[j][i]) ;
        
        fhMCEM02Overlap1[j][i]     = new TH2F(Form("hMCEM02Overlap1NLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                              Form("Overlap 1, #lambda_{0}^{2} vs E for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                              nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhMCEM02Overlap1[j][i]   ->SetYTitle("#lambda_{0}^{2}");
        fhMCEM02Overlap1[j][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCEM02Overlap1[j][i]) ;
        
        fhMCEM02OverlapN[j][i]     = new TH2F(Form("hMCEM02OverlapNNLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                              Form("Overlap N, #lambda_{0}^{2} vs E for NLM=%s %s",snlm[j].Data(),ptype[i].Data()),
                                              nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhMCEM02OverlapN[j][i]   ->SetYTitle("#lambda_{0}^{2}");
        fhMCEM02OverlapN[j][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCEM02OverlapN[j][i]) ;
        
        fhMCEMassOverlap0[j][i]     = new TH2F(Form("hMCEMassOverlap0NLocMax%s%s",snlm[j].Data(),pname[i].Data()),
                                               Form("Overlap 0, Mass vs E for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                               nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMCEMassOverlap0[j][i]   ->SetYTitle("Mass (GeV/c^{2}");
        fhMCEMassOverlap0[j][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCEMassOverlap0[j][i]) ;
        
        fhMCEMassOverlap1[j][i]     = new TH2F(Form("hMCEMassOverlap1NLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                               Form("Overalap 1, Mass vs E for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                               nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMCEMassOverlap1[j][i]   ->SetYTitle("Mass (GeV/c^{2}");
        fhMCEMassOverlap1[j][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCEMassOverlap1[j][i]) ;
        
        fhMCEMassOverlapN[j][i]     = new TH2F(Form("hMCEMassOverlapNNLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                               Form("Overlap N, Mass vs E for NLM=%s %s",snlm[j].Data(),ptype[i].Data()),
                                               nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMCEMassOverlapN[j][i]   ->SetYTitle("Mass (GeV/c^{2})");
        fhMCEMassOverlapN[j][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCEMassOverlapN[j][i]) ;
        
        fhMCEAsymOverlap0[j][i]     = new TH2F(Form("hMCEAsymOverlap0NLocMax%s%s",snlm[j].Data(),pname[i].Data()),
                                               Form("Overlap 0, Asymmetry vs E for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                               nptbins,ptmin,ptmax,100,0,1);
        fhMCEAsymOverlap0[j][i]   ->SetYTitle("|A|");
        fhMCEAsymOverlap0[j][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCEAsymOverlap0[j][i]) ;
        
        fhMCEAsymOverlap1[j][i]     = new TH2F(Form("hMCEAsymOverlap1NLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                               Form("Overalap 1, Asymmetry vs E for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                               nptbins,ptmin,ptmax,100,0,1);
        fhMCEAsymOverlap1[j][i]   ->SetYTitle("|A|");
        fhMCEAsymOverlap1[j][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCEAsymOverlap1[j][i]) ;
        
        fhMCEAsymOverlapN[j][i]     = new TH2F(Form("hMCEAsymOverlapNNLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                               Form("Overlap N, Asymmetry vs E for NLM=%s %s",snlm[j].Data(),ptype[i].Data()),
                                               nptbins,ptmin,ptmax,100,0,1);
        fhMCEAsymOverlapN[j][i]   ->SetYTitle("|A|");
        fhMCEAsymOverlapN[j][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCEAsymOverlapN[j][i]) ;

        
        if(fFillNCellHisto)
        {
          fhMCENCellOverlap0[j][i]     = new TH2F(Form("hMCENCellOverlap0NLocMax%s%s",snlm[j].Data(),pname[i].Data()),
                                                  Form("Overlap 0, n cells vs E for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                                  nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
          fhMCENCellOverlap0[j][i]   ->SetYTitle("n cells");
          fhMCENCellOverlap0[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCENCellOverlap0[j][i]) ;
          
          fhMCENCellOverlap1[j][i]     = new TH2F(Form("hMCENCellOverlap1NLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                                  Form("Overalap 1, n cells vs E for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                                  nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
          fhMCENCellOverlap1[j][i]   ->SetYTitle("n cells");
          fhMCENCellOverlap1[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCENCellOverlap1[j][i]) ;
          
          fhMCENCellOverlapN[j][i]     = new TH2F(Form("hMCENCellOverlapNNLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                                  Form("Overlap N, n cells vs E for NLM=%s %s",snlm[j].Data(),ptype[i].Data()),
                                                  nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
          fhMCENCellOverlapN[j][i]   ->SetYTitle("n cells");
          fhMCENCellOverlapN[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCENCellOverlapN[j][i]) ;
        }
        
        fhMCEEpriOverlap0[j][i]     = new TH2F(Form("hMCEEpriOverlap0NLocMax%s%s",snlm[j].Data(),pname[i].Data()),
                                               Form("Overlap 0, E reco vs E prim for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                               nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCEEpriOverlap0[j][i]   ->SetYTitle("E_{gen} (GeV)");
        fhMCEEpriOverlap0[j][i]   ->SetXTitle("E_{reco} (GeV)");
        outputContainer->Add(fhMCEEpriOverlap0[j][i]) ;
        
        fhMCEEpriOverlap1[j][i]     = new TH2F(Form("hMCEEpriOverlap1NLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                               Form("Overalap 1, E reco vs E prim for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                               nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCEEpriOverlap1[j][i]   ->SetYTitle("E_{gen} (GeV)");
        fhMCEEpriOverlap1[j][i]   ->SetXTitle("E_{reco} (GeV)");
        outputContainer->Add(fhMCEEpriOverlap1[j][i]) ;
        
        fhMCEEpriOverlapN[j][i]     = new TH2F(Form("hMCEEpriOverlapNNLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                               Form("Overlap N, E reco vs E prim for NLM=%s %s",snlm[j].Data(),ptype[i].Data()),
                                               nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCEEpriOverlapN[j][i]   ->SetYTitle("E_{gen} (GeV)");
        fhMCEEpriOverlapN[j][i]   ->SetXTitle("E_{reco} (GeV)");
        outputContainer->Add(fhMCEEpriOverlapN[j][i]) ;
        
        
        fhMCEEpriOverlap0IdPi0[j][i]     = new TH2F(Form("hMCEEpriOverlap0IdPi0NLocMax%s%s",snlm[j].Data(),pname[i].Data()),
                                               Form("Overlap 0, E reco vs E prim for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                               nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCEEpriOverlap0IdPi0[j][i]   ->SetYTitle("E_{gen} (GeV)");
        fhMCEEpriOverlap0IdPi0[j][i]   ->SetXTitle("E_{reco} (GeV)");
        outputContainer->Add(fhMCEEpriOverlap0IdPi0[j][i]) ;
        
        fhMCEEpriOverlap1IdPi0[j][i]     = new TH2F(Form("hMCEEpriOverlap1IdPi0NLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                               Form("Overalap 1, E reco vs E prim for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                               nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCEEpriOverlap1IdPi0[j][i]   ->SetYTitle("E_{gen} (GeV)");
        fhMCEEpriOverlap1IdPi0[j][i]   ->SetXTitle("E_{reco} (GeV)");
        outputContainer->Add(fhMCEEpriOverlap1IdPi0[j][i]) ;
        
        fhMCEEpriOverlapNIdPi0[j][i]     = new TH2F(Form("hMCEEpriOverlapNIdPi0NLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                               Form("Overlap N, E reco vs E prim for NLM=%s %s",snlm[j].Data(),ptype[i].Data()),
                                               nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCEEpriOverlapNIdPi0[j][i]   ->SetYTitle("E_{gen} (GeV)");
        fhMCEEpriOverlapNIdPi0[j][i]   ->SetXTitle("E_{reco} (GeV)");
        outputContainer->Add(fhMCEEpriOverlapNIdPi0[j][i]) ;

        
        fhMCESplitEFracOverlap0[j][i]     = new TH2F(Form("hMCESplitEFracOverlap0NLocMax%s%s",snlm[j].Data(),pname[i].Data()),
                                                     Form("Overlap 0, SplitEFrac vs E for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                                     nptbins,ptmin,ptmax,120,0,1.2);
        fhMCESplitEFracOverlap0[j][i]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
        fhMCESplitEFracOverlap0[j][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCESplitEFracOverlap0[j][i]) ;
        
        fhMCESplitEFracOverlap1[j][i]     = new TH2F(Form("hMCESplitEFracOverlap1NLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                                     Form("Overalap 1, SplitEFrac vs E for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                                     nptbins,ptmin,ptmax,120,0,1.2);
        fhMCESplitEFracOverlap1[j][i]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
        fhMCESplitEFracOverlap1[j][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCESplitEFracOverlap1[j][i]) ;
        
        fhMCESplitEFracOverlapN[j][i]     = new TH2F(Form("hMCESplitEFracOverlapNNLocMax%s%s",snlm[j].Data(), pname[i].Data()),
                                                     Form("Overlap N, SplitEFrac vs E for NLM=%s %s",snlm[j].Data(),ptype[i].Data()),
                                                     nptbins,ptmin,ptmax,120,0,1.2);
        fhMCESplitEFracOverlapN[j][i]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
        fhMCESplitEFracOverlapN[j][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCESplitEFracOverlapN[j][i]) ;
        
        if(i < 5)
        {
          fhMCPi0MassM02Overlap0[j][i-1]  = new TH2F(Form("hMCPi0MassM02Overlap0NLocMax%sEbin%d",snlm[j].Data(),i-1),
                                                     Form("Overlap 0, Mass vs #lambda_{0}^{2}, NLM=%s, %s",snlm[j].Data(),sEBin[i-1].Data()),
                                                     ssbins,ssmin,ssmax,mbins,mmin,mmax);
          fhMCPi0MassM02Overlap0[j][i-1]->SetYTitle("M (GeV/c^{2})");
          fhMCPi0MassM02Overlap0[j][i-1]->SetXTitle("#lambda_{0}^{2}");
          outputContainer->Add(fhMCPi0MassM02Overlap0[j][i-1]) ;
          
          fhMCPi0MassM02Overlap1[j][i-1]  = new TH2F(Form("hMCPi0MassM02Overlap1NLocMax%sEbin%d",snlm[j].Data(),i-1),
                                                     Form("Overlap 1, Mass vs #lambda_{0}^{2}, NLM=%s, %s",snlm[j].Data(),sEBin[i-1].Data()),
                                                     ssbins,ssmin,ssmax,mbins,mmin,mmax);
          fhMCPi0MassM02Overlap1[j][i-1]->SetYTitle("M (GeV/c^{2})");
          fhMCPi0MassM02Overlap1[j][i-1]->SetXTitle("#lambda_{0}^{2}");
          outputContainer->Add(fhMCPi0MassM02Overlap1[j][i-1]) ;
          
          fhMCPi0MassM02OverlapN[j][i-1]  = new TH2F(Form("hMCPi0MassM02OverlapNNLocMax%sEbin%d",snlm[j].Data(),i-1),
                                                     Form("Overlap N, Mass vs #lambda_{0}^{2}, NLM=%s, %s",snlm[j].Data(),sEBin[i-1].Data()),
                                                     ssbins,ssmin,ssmax,mbins,mmin,mmax);
          fhMCPi0MassM02OverlapN[j][i-1]->SetYTitle("M (GeV/c^{2})");
          fhMCPi0MassM02OverlapN[j][i-1]->SetXTitle("#lambda_{0}^{2}");
          outputContainer->Add(fhMCPi0MassM02OverlapN[j][i-1]) ;
        }
        
        if(fFillTMHisto)
        {
          fhMCENOverlapsMatch[j][i]     = new TH2F(Form("hMCENOverlapsNLocMax%s%sMatched",snlm[j].Data(),pname[i].Data()),
                                                   Form("# overlaps vs E for NLM=%s, %s",snlm[j].Data(),ptype[i].Data()),
                                                   nptbins,ptmin,ptmax,10,0,10);
          fhMCENOverlapsMatch[j][i]   ->SetYTitle("# overlaps");
          fhMCENOverlapsMatch[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCENOverlapsMatch[j][i]) ;
          
          fhMCEM02Overlap0Match[j][i]     = new TH2F(Form("hMCEM02Overlap0NLocMax%s%sMatched",snlm[j].Data(),pname[i].Data()),
                                                     Form("#lambda_{0}^{2} vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhMCEM02Overlap0Match[j][i]   ->SetYTitle("#lambda_{0}^{2}");
          fhMCEM02Overlap0Match[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCEM02Overlap0Match[j][i]) ;
          
          fhMCEM02Overlap1Match[j][i]     = new TH2F(Form("hMCEM02Overlap1NLocMax%s%sMatched",snlm[j].Data(), pname[i].Data()),
                                                     Form("#lambda_{0}^{2} vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhMCEM02Overlap1Match[j][i]   ->SetYTitle("#lambda_{0}^{2}");
          fhMCEM02Overlap1Match[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCEM02Overlap1Match[j][i]) ;
          
          fhMCEM02OverlapNMatch[j][i]     = new TH2F(Form("hMCEM02OverlapNNLocMax%s%sMatched",snlm[j].Data(), pname[i].Data()),
                                                     Form("#lambda_{0}^{2} vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhMCEM02OverlapNMatch[j][i]   ->SetYTitle("#lambda_{0}^{2}");
          fhMCEM02OverlapNMatch[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCEM02OverlapNMatch[j][i]) ;
          
          fhMCEMassOverlap0Match[j][i]     = new TH2F(Form("hMCEMassOverlap0NLocMax%s%sMatched",snlm[j].Data(),pname[i].Data()),
                                                      Form("Mass vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                      nptbins,ptmin,ptmax,mbins,mmin,mmax);
          fhMCEMassOverlap0Match[j][i]   ->SetYTitle("Mass (GeV/c^{2}");
          fhMCEMassOverlap0Match[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCEMassOverlap0Match[j][i]) ;
          
          fhMCEMassOverlap1Match[j][i]     = new TH2F(Form("hMCEMassOverlap1NLocMax%s%sMatched",snlm[j].Data(), pname[i].Data()),
                                                      Form("Mass vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                      nptbins,ptmin,ptmax,mbins,mmin,mmax);
          fhMCEMassOverlap1Match[j][i]   ->SetYTitle("Mass (GeV/c^{2}");
          fhMCEMassOverlap1Match[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCEMassOverlap1Match[j][i]) ;
          
          fhMCEMassOverlapNMatch[j][i]     = new TH2F(Form("hMCEMassOverlapNNLocMax%s%sMatched",snlm[j].Data(), pname[i].Data()),
                                                      Form("Mass vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                      nptbins,ptmin,ptmax,mbins,mmin,mmax);
          fhMCEMassOverlapNMatch[j][i]   ->SetYTitle("Mass (GeV/c^{2}");
          fhMCEMassOverlapNMatch[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCEMassOverlapNMatch[j][i]) ;
          
          
          fhMCEAsymOverlap0Match[j][i]     = new TH2F(Form("hMCEAsymOverlap0NLocMax%s%sMatched",snlm[j].Data(),pname[i].Data()),
                                                      Form("Overlap 0, Asymmetry vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                      nptbins,ptmin,ptmax,100,0,1);
          fhMCEAsymOverlap0Match[j][i]   ->SetYTitle("|A|");
          fhMCEAsymOverlap0Match[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCEAsymOverlap0Match[j][i]) ;
          
          fhMCEAsymOverlap1Match[j][i]     = new TH2F(Form("hMCEAsymOverlap1NLocMax%s%sMatched",snlm[j].Data(), pname[i].Data()),
                                                      Form("Overalap 1, Asymmetry vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                      nptbins,ptmin,ptmax,100,0,1);
          fhMCEAsymOverlap1Match[j][i]   ->SetYTitle("|A|");
          fhMCEAsymOverlap1Match[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCEAsymOverlap1Match[j][i]) ;
          
          fhMCEAsymOverlapNMatch[j][i]     = new TH2F(Form("hMCEAsymOverlapNNLocMax%s%sMatched",snlm[j].Data(), pname[i].Data()),
                                                      Form("Overlap N, Asymmetry vs E for NLM=%s %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                      nptbins,ptmin,ptmax,100,0,1);
          fhMCEAsymOverlapNMatch[j][i]   ->SetYTitle("|A|");
          fhMCEAsymOverlapNMatch[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCEAsymOverlapNMatch[j][i]) ;

          
          fhMCENCellOverlap0Match[j][i]     = new TH2F(Form("hMCENCellOverlap0NLocMax%s%sMatched",snlm[j].Data(),pname[i].Data()),
                                                      Form("Overlap 0, n cells vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                      nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
          fhMCENCellOverlap0Match[j][i]   ->SetYTitle("n cells");
          fhMCENCellOverlap0Match[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCENCellOverlap0Match[j][i]) ;
          
          fhMCENCellOverlap1Match[j][i]     = new TH2F(Form("hMCENCellOverlap1NLocMax%s%sMatched",snlm[j].Data(), pname[i].Data()),
                                                      Form("Overalap 1, n cell vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                      nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
          fhMCENCellOverlap1Match[j][i]   ->SetYTitle("n cells");
          fhMCENCellOverlap1Match[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCENCellOverlap1Match[j][i]) ;
          
          fhMCENCellOverlapNMatch[j][i]     = new TH2F(Form("hMCENCellOverlapNNLocMax%s%sMatched",snlm[j].Data(), pname[i].Data()),
                                                      Form("Overlap N, n cell vs E for NLM=%s %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                      nptbins,ptmin,ptmax,ncbins,ncmin,ncmax);
          fhMCENCellOverlapNMatch[j][i]   ->SetYTitle("n cells");
          fhMCENCellOverlapNMatch[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCENCellOverlapNMatch[j][i]) ;

          fhMCEEpriOverlap0Match[j][i]     = new TH2F(Form("hMCEEpriOverlap0NLocMax%s%sMatched",snlm[j].Data(),pname[i].Data()),
                                                      Form("Overlap 0, Asymmetry vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                      nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
          fhMCEEpriOverlap0Match[j][i]   ->SetYTitle("E_{gen} (GeV)");
          fhMCEEpriOverlap0Match[j][i]   ->SetXTitle("E_{reco} (GeV)");
          outputContainer->Add(fhMCEEpriOverlap0Match[j][i]) ;
          
          fhMCEEpriOverlap1Match[j][i]     = new TH2F(Form("hMCEEpriOverlap1NLocMax%s%sMatched",snlm[j].Data(), pname[i].Data()),
                                                      Form("Overalap 1, Asymmetry vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                      nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
          fhMCEEpriOverlap1Match[j][i]   ->SetYTitle("E_{gen} (GeV)");
          fhMCEEpriOverlap1Match[j][i]   ->SetXTitle("E_{reco} (GeV)");
          outputContainer->Add(fhMCEEpriOverlap1Match[j][i]) ;
          
          fhMCEEpriOverlapNMatch[j][i]     = new TH2F(Form("hMCEEpriOverlapNNLocMax%s%sMatched",snlm[j].Data(), pname[i].Data()),
                                                      Form("Overlap N, Asymmetry vs E for NLM=%s %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                      nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
          fhMCEEpriOverlapNMatch[j][i]   ->SetYTitle("E_{gen} (GeV)");
          fhMCEEpriOverlapNMatch[j][i]   ->SetXTitle("E_{reco} (GeV)");
          outputContainer->Add(fhMCEEpriOverlapNMatch[j][i]) ;
          
          
          fhMCESplitEFracOverlap0Match[j][i]     = new TH2F(Form("hMCESplitEFracOverlap0NLocMax%s%sMatched",snlm[j].Data(),pname[i].Data()),
                                                            Form("SplitEFrac vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                            nptbins,ptmin,ptmax,120,0,1.2);
          fhMCESplitEFracOverlap0Match[j][i]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
          fhMCESplitEFracOverlap0Match[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCESplitEFracOverlap0Match[j][i]) ;
          
          fhMCESplitEFracOverlap1Match[j][i]     = new TH2F(Form("hMCESplitEFracOverlap1NLocMax%s%sMatched",snlm[j].Data(), pname[i].Data()),
                                                            Form("SplitEFrac vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                            nptbins,ptmin,ptmax,120,0,1.2);
          fhMCESplitEFracOverlap1Match[j][i]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
          fhMCESplitEFracOverlap1Match[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCESplitEFracOverlap1Match[j][i]) ;
          
          fhMCESplitEFracOverlapNMatch[j][i]     = new TH2F(Form("hMCESplitEFracOverlapNNLocMax%s%sMatched",snlm[j].Data(), pname[i].Data()),
                                                            Form("SplitEFrac vs E for NLM=%s, %s, Track Matched",snlm[j].Data(),ptype[i].Data()),
                                                            nptbins,ptmin,ptmax,120,0,1.2);
          fhMCESplitEFracOverlapNMatch[j][i]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
          fhMCESplitEFracOverlapNMatch[j][i]   ->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCESplitEFracOverlapNMatch[j][i]) ;
          
          
          if(i < 5)
          {
            fhMCPi0MassM02Overlap0Match[j][i-1]  = new TH2F(Form("hMCPi0MassM02Overlap0NLocMax%sEbin%dMatched",snlm[j].Data(),i-1),
                                                            Form("Overlap 0, Mass vs #lambda_{0}^{2}, NLM=%s, %s, Track Matched",snlm[j].Data(),sEBin[i-1].Data()),
                                                            ssbins,ssmin,ssmax,mbins,mmin,mmax);
            fhMCPi0MassM02Overlap0Match[j][i-1]->SetYTitle("M (GeV/c^{2})");
            fhMCPi0MassM02Overlap0Match[j][i-1]->SetXTitle("#lambda_{0}^{2}");
            outputContainer->Add(fhMCPi0MassM02Overlap0Match[j][i-1]) ;
            
            fhMCPi0MassM02Overlap1Match[j][i-1]  = new TH2F(Form("hMCPi0MassM02Overlap1NLocMax%sEbin%dMatched",snlm[j].Data(),i-1),
                                                            Form("Overlap 1, Mass vs #lambda_{0}^{2}, NLM=%s, %s, Track Matched",snlm[j].Data(),sEBin[i-1].Data()),
                                                            ssbins,ssmin,ssmax,mbins,mmin,mmax);
            fhMCPi0MassM02Overlap1Match[j][i-1]->SetYTitle("M (GeV/c^{2})");
            fhMCPi0MassM02Overlap1Match[j][i-1]->SetXTitle("#lambda_{0}^{2}");
            outputContainer->Add(fhMCPi0MassM02Overlap1Match[j][i-1]) ;
            
            fhMCPi0MassM02OverlapNMatch[j][i-1]  = new TH2F(Form("hMCPi0MassM02OverlapNNLocMax%sEbin%dMatched",snlm[j].Data(),i-1),
                                                            Form("Overlap N, Mass vs #lambda_{0}^{2}, NLM=%s, %s, Track Matched",snlm[j].Data(),sEBin[i-1].Data()),
                                                            ssbins,ssmin,ssmax,mbins,mmin,mmax);
            fhMCPi0MassM02OverlapNMatch[j][i-1]->SetYTitle("M (GeV/c^{2})");
            fhMCPi0MassM02OverlapNMatch[j][i-1]->SetXTitle("#lambda_{0}^{2}");
            outputContainer->Add(fhMCPi0MassM02OverlapNMatch[j][i-1]) ;
            
          }
          
        }
      }
    }
    
    fhMCPi0HighNLMPair    = new TH2F("hMCPi0HighNLMPair","NLM vs E for merged pi0 cluster, high energy NLM pair are decays",
                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0HighNLMPair   ->SetYTitle("N maxima");
    fhMCPi0HighNLMPair   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0HighNLMPair) ;
    
    fhMCPi0LowNLMPair     = new TH2F("hMCPi0LowNLMPair","NLM vs E for merged pi0 cluster, lower energy NLM pair are decays",
                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0LowNLMPair   ->SetYTitle("N maxima");
    fhMCPi0LowNLMPair   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0LowNLMPair) ;
    
    fhMCPi0AnyNLMPair     = new TH2F("hMCPi0AnyNLMPair","NLM vs E for merged pi0 cluster, both high and other energy NLM pair are decays",
                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0AnyNLMPair   ->SetYTitle("N maxima");
    fhMCPi0AnyNLMPair   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0AnyNLMPair) ;
    
    fhMCPi0NoneNLMPair     = new TH2F("hMCPi0NoneNLMPair","NLM vs E for merged pi0 cluster, no NLM pair are decays",
                                      nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0NoneNLMPair   ->SetYTitle("N maxima");
    fhMCPi0NoneNLMPair   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0NoneNLMPair) ;
    
    
    fhMCPi0HighNLMPairNoMCMatch    = new TH2F("hMCPi0HighNLMPairNoMCMatch","NLM vs E for merged pi0 cluster, high energy NLM pair are decays",
                                              nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0HighNLMPairNoMCMatch   ->SetYTitle("N maxima");
    fhMCPi0HighNLMPairNoMCMatch   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0HighNLMPairNoMCMatch) ;
    
    fhMCPi0LowNLMPairNoMCMatch     = new TH2F("hMCPi0LowNLMPairNoMCMatch","NLM vs E for merged pi0 cluster, lower energy NLM pair are decays",
                                              nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0LowNLMPairNoMCMatch   ->SetYTitle("N maxima");
    fhMCPi0LowNLMPairNoMCMatch   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0LowNLMPairNoMCMatch) ;
    
    fhMCPi0AnyNLMPairNoMCMatch     = new TH2F("hMCPi0AnyNLMPairNoMCMatch","NLM vs E for merged pi0 cluster, both high and other energy NLM pair are decays",
                                              nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0AnyNLMPairNoMCMatch   ->SetYTitle("N maxima");
    fhMCPi0AnyNLMPairNoMCMatch   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0AnyNLMPairNoMCMatch) ;
    
    fhMCPi0NoneNLMPairNoMCMatch     = new TH2F("hMCPi0NoneNLMPairNoMCMatch","NLM vs E for merged pi0 cluster, no NLM pair are decays",
                                               nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0NoneNLMPairNoMCMatch   ->SetYTitle("N maxima");
    fhMCPi0NoneNLMPairNoMCMatch   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0NoneNLMPairNoMCMatch) ;
    
    
    fhMCPi0HighNLMPairOverlap    = new TH2F("hMCPi0HighNLMPairOverlap","NLM vs E for merged pi0 cluster, high energy NLM pair are decays",
                                            nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0HighNLMPairOverlap   ->SetYTitle("N maxima");
    fhMCPi0HighNLMPairOverlap   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0HighNLMPairOverlap) ;
    
    fhMCPi0LowNLMPairOverlap     = new TH2F("hMCPi0LowNLMPairOverlap","NLM vs E for merged pi0 cluster, lower energy NLM pair are decays",
                                            nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0LowNLMPairOverlap   ->SetYTitle("N maxima");
    fhMCPi0LowNLMPairOverlap   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0LowNLMPairOverlap) ;
    
    fhMCPi0AnyNLMPairOverlap     = new TH2F("hMCPi0AnyNLMPairOverlap","NLM vs E for merged pi0 cluster, both high and other energy NLM pair are decays",
                                            nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0AnyNLMPairOverlap   ->SetYTitle("N maxima");
    fhMCPi0AnyNLMPairOverlap   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0AnyNLMPairOverlap) ;
    
    fhMCPi0NoneNLMPairOverlap     = new TH2F("hMCPi0NoneNLMPairOverlap","NLM vs E for merged pi0 cluster, no NLM pair are decays",
                                             nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0NoneNLMPairOverlap   ->SetYTitle("N maxima");
    fhMCPi0NoneNLMPairOverlap   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0NoneNLMPairOverlap) ;
    
    fhMCPi0HighNLMPairNoMCMatchOverlap    = new TH2F("hMCPi0HighNLMPairNoMCMatchOverlap","NLM vs E for merged pi0 cluster, high energy NLM pair are decays",
                                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0HighNLMPairNoMCMatchOverlap   ->SetYTitle("N maxima");
    fhMCPi0HighNLMPairNoMCMatchOverlap   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0HighNLMPairNoMCMatchOverlap) ;
    
    fhMCPi0LowNLMPairNoMCMatchOverlap     = new TH2F("hMCPi0LowNLMPairNoMCMatchOverlap","NLM vs E for merged pi0 cluster, lower energy NLM pair are decays",
                                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0LowNLMPairNoMCMatchOverlap   ->SetYTitle("N maxima");
    fhMCPi0LowNLMPairNoMCMatchOverlap   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0LowNLMPairNoMCMatchOverlap) ;
    
    fhMCPi0AnyNLMPairNoMCMatchOverlap     = new TH2F("hMCPi0AnyNLMPairNoMCMatchOverlap","NLM vs E for merged pi0 cluster, both high and other energy NLM pair are decays",
                                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0AnyNLMPairNoMCMatchOverlap   ->SetYTitle("N maxima");
    fhMCPi0AnyNLMPairNoMCMatchOverlap   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0AnyNLMPairNoMCMatchOverlap) ;
    
    fhMCPi0NoneNLMPairNoMCMatchOverlap     = new TH2F("hMCPi0NoneNLMPairNoMCMatchOverlap","NLM vs E for merged pi0 cluster, no NLM pair are decays",
                                                      nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0NoneNLMPairNoMCMatchOverlap   ->SetYTitle("N maxima");
    fhMCPi0NoneNLMPairNoMCMatchOverlap   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0NoneNLMPairNoMCMatchOverlap) ;
  
    
    fhMCPi0DecayPhotonHitHighLM    = new TH2F("hMCPi0DecayPhotonHitHighLM","NLM vs E for merged pi0 cluster, decay photon hit High Local Maxima",
                                              nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0DecayPhotonHitHighLM    ->SetYTitle("N maxima");
    fhMCPi0DecayPhotonHitHighLM    ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0DecayPhotonHitHighLM ) ;
    
    fhMCPi0DecayPhotonAdjHighLM    = new TH2F("hMCPi0DecayPhotonAdjHighLM","NLM vs E for merged pi0 cluster, decay photon hit cells adjacent to High Local Maxima",
                                              nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0DecayPhotonAdjHighLM    ->SetYTitle("N maxima");
    fhMCPi0DecayPhotonAdjHighLM    ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0DecayPhotonAdjHighLM ) ;
    
    fhMCPi0DecayPhotonHitOtherLM    = new TH2F("hMCPi0DecayPhotonHitOtherLM","NLM vs E for merged pi0 cluster, decay photon hit Other Local Maxima",
                                               nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0DecayPhotonHitOtherLM    ->SetYTitle("N maxima");
    fhMCPi0DecayPhotonHitOtherLM    ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0DecayPhotonHitOtherLM ) ;
    
    fhMCPi0DecayPhotonAdjOtherLM    = new TH2F("hMCPi0DecayPhotonAdjOtherLM","NLM vs E for merged pi0 cluster, decay photon hit cells adjacent to Other Local Maxima",
                                               nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0DecayPhotonAdjOtherLM    ->SetYTitle("N maxima");
    fhMCPi0DecayPhotonAdjOtherLM    ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0DecayPhotonAdjOtherLM ) ;
    
    fhMCPi0DecayPhotonAdjacent   = new TH2F("hMCPi0DecayPhotonAdjacent","NLM vs E for merged pi0 cluster, decay photon hit adjacent cells",
                                            nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0DecayPhotonAdjacent    ->SetYTitle("N maxima");
    fhMCPi0DecayPhotonAdjacent    ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0DecayPhotonAdjacent ) ;
    
    fhMCPi0DecayPhotonHitNoLM    = new TH2F("hMCPi0DecayPhotonHitNoLM","NLM vs E for merged pi0 cluster, decay photon do not hit Local Maxima",
                                            nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0DecayPhotonHitNoLM    ->SetYTitle("N maxima");
    fhMCPi0DecayPhotonHitNoLM    ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0DecayPhotonHitNoLM ) ;
    
    
    fhMCPi0DecayPhotonHitHighLMOverlap    = new TH2F("hMCPi0DecayPhotonHitHighLMOverlap","NLM vs E for merged pi0 cluster, decay photon hit High Local Maxima, there was an overlap",
                                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0DecayPhotonHitHighLMOverlap    ->SetYTitle("N maxima");
    fhMCPi0DecayPhotonHitHighLMOverlap    ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0DecayPhotonHitHighLMOverlap ) ;
    
    fhMCPi0DecayPhotonAdjHighLMOverlap    = new TH2F("hMCPi0DecayPhotonAdjHighLMOverlap","NLM vs E for merged pi0 cluster, decay photon hit cells adjacent to High Local Maxima, there was an overlap",
                                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0DecayPhotonAdjHighLMOverlap    ->SetYTitle("N maxima");
    fhMCPi0DecayPhotonAdjHighLMOverlap    ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0DecayPhotonAdjHighLMOverlap ) ;
    
    fhMCPi0DecayPhotonHitOtherLMOverlap    = new TH2F("hMCPi0DecayPhotonHitOtherLMOverlap","NLM vs E for merged pi0 cluster, decay photon hit Other Local Maxima, there was an overlap",
                                                      nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0DecayPhotonHitOtherLMOverlap    ->SetYTitle("N maxima");
    fhMCPi0DecayPhotonHitOtherLMOverlap    ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0DecayPhotonHitOtherLMOverlap ) ;
    
    fhMCPi0DecayPhotonAdjOtherLMOverlap    = new TH2F("hMCPi0DecayPhotonAdjOtherLMOverlap","NLM vs E for merged pi0 cluster, decay photon hit cells adjacent to Other Local Maxima, there was an overlap",
                                                      nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0DecayPhotonAdjOtherLMOverlap    ->SetYTitle("N maxima");
    fhMCPi0DecayPhotonAdjOtherLMOverlap    ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0DecayPhotonAdjOtherLMOverlap ) ;
    
    fhMCPi0DecayPhotonAdjacentOverlap   = new TH2F("hMCPi0DecayPhotonAdjacentOverlap","NLM vs E for merged pi0 cluster, decay photon hit adjacent cells, there was an overlap",
                                                   nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0DecayPhotonAdjacentOverlap    ->SetYTitle("N maxima");
    fhMCPi0DecayPhotonAdjacentOverlap    ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0DecayPhotonAdjacentOverlap ) ;
    
    fhMCPi0DecayPhotonHitNoLMOverlap    = new TH2F("hMCPi0DecayPhotonHitNoLMOverlap","NLM vs E for merged pi0 cluster, decay photon do not hit Local Maxima, there was an overlap",
                                                   nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins);
    fhMCPi0DecayPhotonHitNoLMOverlap    ->SetYTitle("N maxima");
    fhMCPi0DecayPhotonHitNoLMOverlap    ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCPi0DecayPhotonHitNoLMOverlap ) ;
    
    
    for(Int_t nlm = 0; nlm < 3; nlm++)
    {
      fhMCPi0DecayPhotonHitHighLMMass[nlm]  = new TH2F(Form("hMCPi0DecayPhotonHitHighLM%sMass",snlm[nlm].Data()),
                                                       Form("Mass vs E for merged pi0 cluster, NLM=%s, decay photon hit High Local Maxima",snlm[nlm].Data()),
                                                       nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMCPi0DecayPhotonHitHighLMMass[nlm]  ->SetYTitle("Mass (MeV/c^{2})");
      fhMCPi0DecayPhotonHitHighLMMass[nlm]  ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitHighLMMass[nlm] ) ;
      
      fhMCPi0DecayPhotonAdjHighLMMass[nlm]  = new TH2F(Form("hMCPi0DecayPhotonAdjHighLM%sMass",snlm[nlm].Data()),
                                                       Form("Mass vs E for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to High Local Maxima",snlm[nlm].Data()),
                                                       nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMCPi0DecayPhotonAdjHighLMMass[nlm]  ->SetYTitle("Mass (MeV/c^{2})");
      fhMCPi0DecayPhotonAdjHighLMMass[nlm]  ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjHighLMMass[nlm] ) ;
      
      fhMCPi0DecayPhotonHitOtherLMMass[nlm]  = new TH2F(Form("hMCPi0DecayPhotonHitOtherLM%sMass",snlm[nlm].Data()),
                                                        Form("Mass vs E for merged pi0 cluster, NLM=%s, decay photon hit Other Local Maxima",snlm[nlm].Data()),
                                                        nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMCPi0DecayPhotonHitOtherLMMass[nlm]  ->SetYTitle("Mass (MeV/c^{2})");
      fhMCPi0DecayPhotonHitOtherLMMass[nlm]  ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitOtherLMMass[nlm] ) ;
      
      fhMCPi0DecayPhotonAdjOtherLMMass[nlm]  = new TH2F(Form("hMCPi0DecayPhotonAdjOtherLM%sMass",snlm[nlm].Data()),
                                                        Form("Mass vs E for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to Other Local Maxima",snlm[nlm].Data()),
                                                        nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMCPi0DecayPhotonAdjOtherLMMass[nlm]  ->SetYTitle("Mass (MeV/c^{2})");
      fhMCPi0DecayPhotonAdjOtherLMMass[nlm]  ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjOtherLMMass[nlm] ) ;
      
      fhMCPi0DecayPhotonAdjacentMass[nlm]   = new TH2F(Form("hMCPi0DecayPhotonAdjacentLM%sMass",snlm[nlm].Data()),
                                                       Form("Mass vs E for merged pi0 cluster, NLM=%s, decay photon hit adjacent cells",snlm[nlm].Data()),
                                                       nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMCPi0DecayPhotonAdjacentMass[nlm]    ->SetYTitle("Mass (MeV/c^{2})");
      fhMCPi0DecayPhotonAdjacentMass[nlm]    ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjacentMass[nlm] ) ;
      
      fhMCPi0DecayPhotonHitNoLMMass[nlm]  = new TH2F(Form("hMCPi0DecayPhotonHitNoLM%sMass",snlm[nlm].Data()),
                                                     Form("Mass vs E for merged pi0 cluster, NLM=%s, decay photon do not hit Local Maxima",snlm[nlm].Data()),
                                                     nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMCPi0DecayPhotonHitNoLMMass[nlm]  ->SetYTitle("Mass (MeV/c^{2})");
      fhMCPi0DecayPhotonHitNoLMMass[nlm]    ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitNoLMMass[nlm] ) ;
      
      
      fhMCPi0DecayPhotonHitHighLMOverlapMass[nlm]   = new TH2F(Form("hMCPi0DecayPhotonHitHighLM%sOverlapMass",snlm[nlm].Data()),
                                                               Form("Mass vs E for merged pi0 cluster, NLM=%s, decay photon hit High Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                               nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMCPi0DecayPhotonHitHighLMOverlapMass[nlm]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMCPi0DecayPhotonHitHighLMOverlapMass[nlm]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitHighLMOverlapMass[nlm]) ;
      
      fhMCPi0DecayPhotonAdjHighLMOverlapMass[nlm]   = new TH2F(Form("hMCPi0DecayPhotonAdjHighLM%sOverlapMass",snlm[nlm].Data()),
                                                               Form("Mass vs E for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to High Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                               nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMCPi0DecayPhotonAdjHighLMOverlapMass[nlm]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMCPi0DecayPhotonAdjHighLMOverlapMass[nlm]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjHighLMOverlapMass[nlm]) ;
      
      fhMCPi0DecayPhotonHitOtherLMOverlapMass[nlm]   = new TH2F(Form("hMCPi0DecayPhotonHitOtherLM%sOverlapMass",snlm[nlm].Data()),
                                                                Form("Mass vs E for merged pi0 cluster, NLM=%s, decay photon hit Other Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMCPi0DecayPhotonHitOtherLMOverlapMass[nlm]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMCPi0DecayPhotonHitOtherLMOverlapMass[nlm]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitOtherLMOverlapMass[nlm]) ;
      
      fhMCPi0DecayPhotonAdjOtherLMOverlapMass[nlm]   = new TH2F(Form("hMCPi0DecayPhotonAdjOtherLM%sOverlapMass",snlm[nlm].Data()),
                                                                Form("Mass vs E for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to Other Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMCPi0DecayPhotonAdjOtherLMOverlapMass[nlm]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMCPi0DecayPhotonAdjOtherLMOverlapMass[nlm]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjOtherLMOverlapMass[nlm]) ;
      
      fhMCPi0DecayPhotonAdjacentOverlapMass[nlm]  = new TH2F(Form("hMCPi0DecayPhotonAdjacentLM%sOverlapMass",snlm[nlm].Data()),
                                                             Form("Mass vs E for merged pi0 cluster, NLM=%s, decay photon hit adjacent cells, there was an overlap",snlm[nlm].Data()),
                                                             nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMCPi0DecayPhotonAdjacentOverlapMass[nlm]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMCPi0DecayPhotonAdjacentOverlapMass[nlm]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjacentOverlapMass[nlm]) ;
      
      fhMCPi0DecayPhotonHitNoLMOverlapMass[nlm]   = new TH2F(Form("hMCPi0DecayPhotonHitNoLM%sOverlapMass",snlm[nlm].Data()),
                                                             Form("Mass vs E for merged pi0 cluster, NLM=%s, decay photon do not hit Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                             nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMCPi0DecayPhotonHitNoLMOverlapMass[nlm]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMCPi0DecayPhotonHitNoLMOverlapMass[nlm]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitNoLMOverlapMass[nlm]) ;
      
      fhMCPi0DecayPhotonHitHighLMDiffELM1[nlm]  = new TH2F(Form("hMCPi0DecayPhotonHitHighLM%sDiffELM1",snlm[nlm].Data()),
                                                           Form("(E_{reco}-E_{gen})/E_{gen} vs E pi0 for merged pi0 cluster, NLM=%s, decay photon hit High Local Maxima",snlm[nlm].Data()),
                                                           nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitHighLMDiffELM1[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitHighLMDiffELM1[nlm]  ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitHighLMDiffELM1[nlm] ) ;
      
      fhMCPi0DecayPhotonAdjHighLMDiffELM1[nlm]  = new TH2F(Form("hMCPi0DecayPhotonAdjHighLM%sDiffELM1",snlm[nlm].Data()),
                                                           Form("(E_{reco}-E_{gen})/E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to High Local Maxima",snlm[nlm].Data()),
                                                           nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjHighLMDiffELM1[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjHighLMDiffELM1[nlm]  ->SetXTitle("E pi0 (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjHighLMDiffELM1[nlm] ) ;
      
      fhMCPi0DecayPhotonHitOtherLMDiffELM1[nlm]  = new TH2F(Form("hMCPi0DecayPhotonHitOtherLM%sDiffELM1",snlm[nlm].Data()),
                                                            Form("(E_{reco}-E_{gen})/E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit Other Local Maxima",snlm[nlm].Data()),
                                                            nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitOtherLMDiffELM1[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitOtherLMDiffELM1[nlm]  ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitOtherLMDiffELM1[nlm] ) ;
      
      fhMCPi0DecayPhotonAdjOtherLMDiffELM1[nlm]  = new TH2F(Form("hMCPi0DecayPhotonAdjOtherLM%sDiffELM1",snlm[nlm].Data()),
                                                            Form("(E_{reco}-E_{gen})/E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to Other Local Maxima",snlm[nlm].Data()),
                                                            nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjOtherLMDiffELM1[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjOtherLMDiffELM1[nlm]  ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjOtherLMDiffELM1[nlm] ) ;
      
      fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1[nlm]   = new TH2F(Form("hMCPi0DecayPhotonHitHighLM%sOverlapDiffELM1",snlm[nlm].Data()),
                                                                   Form("(E_{reco}-E_{gen})/E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit High Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                   nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1[nlm]   ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1[nlm]) ;
      
      fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1[nlm]   = new TH2F(Form("hMCPi0DecayPhotonAdjHighLM%sOverlapDiffELM1",snlm[nlm].Data()),
                                                                   Form("(E_{reco}-E_{gen})/E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to High Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                   nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1[nlm]   ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1[nlm]) ;
      
      fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM1[nlm]   = new TH2F(Form("hMCPi0DecayPhotonHitOtherLM%sOverlapDiffELM1",snlm[nlm].Data()),
                                                                    Form("(E_{reco}-E_{gen})/E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit Other Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                    nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM1[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM1[nlm]   ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM1[nlm]) ;
      
      fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM1[nlm]   = new TH2F(Form("hMCPi0DecayPhotonAdjOtherLM%sOverlapDiffELM1",snlm[nlm].Data()),
                                                                    Form("E_{reco}-E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to Other Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                    nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM1[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM1[nlm]   ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM1[nlm]) ;
      
      fhMCPi0DecayPhotonHitHighLMDiffELM2[nlm]  = new TH2F(Form("hMCPi0DecayPhotonHitHighLM%sDiffELM2",snlm[nlm].Data()),
                                                           Form("(E_{reco}-E_{gen})/E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit High Local Maxima",snlm[nlm].Data()),
                                                           nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitHighLMDiffELM2[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitHighLMDiffELM2[nlm]  ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitHighLMDiffELM2[nlm] ) ;
      
      fhMCPi0DecayPhotonAdjHighLMDiffELM2[nlm]  = new TH2F(Form("hMCPi0DecayPhotonAdjHighLM%sDiffELM2",snlm[nlm].Data()),
                                                           Form("(E_{reco}-E_{gen})/E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to High Local Maxima",snlm[nlm].Data()),
                                                           nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjHighLMDiffELM2[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjHighLMDiffELM2[nlm]  ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjHighLMDiffELM2[nlm] ) ;
      
      fhMCPi0DecayPhotonHitOtherLMDiffELM2[nlm]  = new TH2F(Form("hMCPi0DecayPhotonHitOtherLM%sDiffELM2",snlm[nlm].Data()),
                                                            Form("E_{reco}-E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit Other Local Maxima",snlm[nlm].Data()),
                                                            nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitOtherLMDiffELM2[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitOtherLMDiffELM2[nlm]  ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitOtherLMDiffELM2[nlm] ) ;
      
      fhMCPi0DecayPhotonAdjOtherLMDiffELM2[nlm]  = new TH2F(Form("hMCPi0DecayPhotonAdjOtherLM%sDiffELM2",snlm[nlm].Data()),
                                                            Form("E_{reco}-E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to Other Local Maxima",snlm[nlm].Data()),
                                                            nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjOtherLMDiffELM2[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjOtherLMDiffELM2[nlm]  ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjOtherLMDiffELM2[nlm] ) ;
      
      fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2[nlm]   = new TH2F(Form("hMCPi0DecayPhotonHitHighLM%sOverlapDiffELM2",snlm[nlm].Data()),
                                                                   Form("(E_{reco}-E_{gen})/E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit High Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                   nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2[nlm]   ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2[nlm]) ;
      
      fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2[nlm]   = new TH2F(Form("hMCPi0DecayPhotonAdjHighLM%sOverlapDiffELM2",snlm[nlm].Data()),
                                                                   Form("E_{reco}-E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to High Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                   nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2[nlm]   ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2[nlm]) ;
      
      fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM2[nlm]   = new TH2F(Form("hMCPi0DecayPhotonHitOtherLM%sOverlapDiffELM2",snlm[nlm].Data()),
                                                                    Form("(E_{reco}-E_{gen})/E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit Other Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                    nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM2[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM2[nlm]   ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM2[nlm]) ;
      
      fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM2[nlm]   = new TH2F(Form("hMCPi0DecayPhotonAdjOtherLM%sOverlapDiffELM2",snlm[nlm].Data()),
                                                                    Form("(E_{reco}-E_{gen})/E_{gen} vs E_{pi0 reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to Other Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                    nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM2[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM2[nlm]   ->SetXTitle("E_{pi0 reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM2[nlm]) ;
      
      
      fhMCPi0DecayPhotonHitHighLMDiffELM1vsELM1[nlm]  = new TH2F(Form("hMCPi0DecayPhotonHitHighLM%sDiffELM1vsELM1",snlm[nlm].Data()),
                                                                 Form("(E_{reco}-E_{gen})/E_{gen} vs E for merged pi0 cluster, NLM=%s, decay photon hit High Local Maxima",snlm[nlm].Data()),
                                                                 nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitHighLMDiffELM1vsELM1[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitHighLMDiffELM1vsELM1[nlm]  ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitHighLMDiffELM1vsELM1[nlm] ) ;
      
      fhMCPi0DecayPhotonAdjHighLMDiffELM1vsELM1[nlm]  = new TH2F(Form("hMCPi0DecayPhotonAdjHighLM%sDiffELM1vsELM1",snlm[nlm].Data()),
                                                                 Form("(E_{reco}-E_{gen})/E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to High Local Maxima",snlm[nlm].Data()),
                                                                 nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjHighLMDiffELM1vsELM1[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjHighLMDiffELM1vsELM1[nlm]  ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjHighLMDiffELM1vsELM1[nlm] ) ;
      
      fhMCPi0DecayPhotonHitOtherLMDiffELM1vsELM1[nlm]  = new TH2F(Form("hMCPi0DecayPhotonHitOtherLM%sDiffELM1vsELM1",snlm[nlm].Data()),
                                                                  Form("(E_{reco}-E_{gen})/E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit Other Local Maxima",snlm[nlm].Data()),
                                                                  nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitOtherLMDiffELM1vsELM1[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitOtherLMDiffELM1vsELM1[nlm]  ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitOtherLMDiffELM1vsELM1[nlm] ) ;
      
      fhMCPi0DecayPhotonAdjOtherLMDiffELM1vsELM1[nlm]  = new TH2F(Form("hMCPi0DecayPhotonAdjOtherLM%sDiffELM1vsELM1",snlm[nlm].Data()),
                                                                  Form("(E_{reco}-E_{gen})/E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to Other Local Maxima",snlm[nlm].Data()),
                                                                  nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjOtherLMDiffELM1vsELM1[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjOtherLMDiffELM1vsELM1[nlm]  ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjOtherLMDiffELM1vsELM1[nlm] ) ;
      
      fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1vsELM1[nlm]   = new TH2F(Form("hMCPi0DecayPhotonHitHighLM%sOverlapDiffELM1vsELM1",snlm[nlm].Data()),
                                                                         Form("(E_{reco}-E_{gen})/E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit High Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                         nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1vsELM1[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1vsELM1[nlm]   ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1vsELM1[nlm]) ;
      
      fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1vsELM1[nlm]   = new TH2F(Form("hMCPi0DecayPhotonAdjHighLM%sOverlapDiffELM1vsELM1",snlm[nlm].Data()),
                                                                         Form("(E_{reco}-E_{gen})/E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to High Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                         nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1vsELM1[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1vsELM1[nlm]   ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1vsELM1[nlm]) ;
      
      fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM1vsELM1[nlm]   = new TH2F(Form("hMCPi0DecayPhotonHitOtherLM%sOverlapDiffELM1vsELM1",snlm[nlm].Data()),
                                                                          Form("(E_{reco}-E_{gen})/E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit Other Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                          nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM1vsELM1[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM1vsELM1[nlm]   ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM1vsELM1[nlm]) ;
      
      fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM1vsELM1[nlm]   = new TH2F(Form("hMCPi0DecayPhotonAdjOtherLM%sOverlapDiffELM1vsELM1",snlm[nlm].Data()),
                                                                          Form("E_{reco}-E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to Other Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                          nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM1vsELM1[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM1vsELM1[nlm]   ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM1vsELM1[nlm]) ;
      
      fhMCPi0DecayPhotonHitHighLMDiffELM2vsELM2[nlm]  = new TH2F(Form("hMCPi0DecayPhotonHitHighLM%sDiffELM2vsELM2",snlm[nlm].Data()),
                                                                 Form("(E_{reco}-E_{gen})/E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit High Local Maxima",snlm[nlm].Data()),
                                                                 nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitHighLMDiffELM2vsELM2[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitHighLMDiffELM2vsELM2[nlm]  ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitHighLMDiffELM2vsELM2[nlm] ) ;
      
      fhMCPi0DecayPhotonAdjHighLMDiffELM2vsELM2[nlm]  = new TH2F(Form("hMCPi0DecayPhotonAdjHighLM%sDiffELM2vsELM2",snlm[nlm].Data()),
                                                                 Form("(E_{reco}-E_{gen})/E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to High Local Maxima",snlm[nlm].Data()),
                                                                 nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjHighLMDiffELM2vsELM2[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjHighLMDiffELM2vsELM2[nlm]  ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjHighLMDiffELM2vsELM2[nlm] ) ;
      
      fhMCPi0DecayPhotonHitOtherLMDiffELM2vsELM2[nlm]  = new TH2F(Form("hMCPi0DecayPhotonHitOtherLM%sDiffELM2vsELM2",snlm[nlm].Data()),
                                                                  Form("E_{reco}-E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit Other Local Maxima",snlm[nlm].Data()),
                                                                  nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitOtherLMDiffELM2vsELM2[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitOtherLMDiffELM2vsELM2[nlm]  ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitOtherLMDiffELM2vsELM2[nlm] ) ;
      
      fhMCPi0DecayPhotonAdjOtherLMDiffELM2vsELM2[nlm]  = new TH2F(Form("hMCPi0DecayPhotonAdjOtherLM%sDiffELM2vsELM2",snlm[nlm].Data()),
                                                                  Form("E_{reco}-E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to Other Local Maxima",snlm[nlm].Data()),
                                                                  nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjOtherLMDiffELM2vsELM2[nlm]  ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjOtherLMDiffELM2vsELM2[nlm]  ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjOtherLMDiffELM2vsELM2[nlm] ) ;
      
      fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2vsELM2[nlm]   = new TH2F(Form("hMCPi0DecayPhotonHitHighLM%sOverlapDiffELM2vsELM2",snlm[nlm].Data()),
                                                                         Form("(E_{reco}-E_{gen})/E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit High Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                         nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2vsELM2[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2vsELM2[nlm]   ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2vsELM2[nlm]) ;
      
      fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2vsELM2[nlm]   = new TH2F(Form("hMCPi0DecayPhotonAdjHighLM%sOverlapDiffELM2vsELM2",snlm[nlm].Data()),
                                                                         Form("E_{reco}-E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to High Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                         nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2vsELM2[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2vsELM2[nlm]   ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2vsELM2[nlm]) ;
      
      fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM2vsELM2[nlm]   = new TH2F(Form("hMCPi0DecayPhotonHitOtherLM%sOverlapDiffELM2vsELM2",snlm[nlm].Data()),
                                                                          Form("(E_{reco}-E_{gen})/E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit Other Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                          nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM2vsELM2[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM2vsELM2[nlm]   ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM2vsELM2[nlm]) ;
      
      fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM2vsELM2[nlm]   = new TH2F(Form("hMCPi0DecayPhotonAdjOtherLM%sOverlapDiffELM2vsELM2",snlm[nlm].Data()),
                                                                          Form("(E_{reco}-E_{gen})/E_{gen} vs E_{reco} for merged pi0 cluster, NLM=%s, decay photon hit cells adjacent to Other Local Maxima, there was an overlap",snlm[nlm].Data()),
                                                                          nptbins,ptmin,ptmax,200,-2,2);//nptbins,-ptmax/4,ptmax/4);
      fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM2vsELM2[nlm]   ->SetYTitle("(E_{reco}-E_{gen})/E_{gen}");
      fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM2vsELM2[nlm]   ->SetXTitle("E_{reco} (GeV)");
      outputContainer->Add(fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM2vsELM2[nlm]) ;
      
    }
    
    fhMCEOverlapType = new TH2F("hMCEOverlapType","Kind of overlap particle, neutral clusters",
                                nptbins,ptmin,ptmax,5,0,5);
    //fhMCEOverlapType   ->SetYTitle("Overlap Type");
    fhMCEOverlapType->GetYaxis()->SetBinLabel(1 ,"#gamma");
    fhMCEOverlapType->GetYaxis()->SetBinLabel(2 ,"e^{#pm}");
    fhMCEOverlapType->GetYaxis()->SetBinLabel(3 ,"hadron^{#pm}");
    fhMCEOverlapType->GetYaxis()->SetBinLabel(4 ,"hadron^{0}");
    fhMCEOverlapType->GetYaxis()->SetBinLabel(5 ,"??");
    fhMCEOverlapType->SetXTitle("Cluster E (GeV)");
    outputContainer->Add(fhMCEOverlapType) ;
    
    fhMCEOverlapTypeMatch = new TH2F("hMCEOverlapTypeMatched","Kind of overlap particle, charged clusters",
                                     nptbins,ptmin,ptmax,5,0,5);
    //fhMCEOverlapTypeMatch   ->SetYTitle("Overlap Type");
    fhMCEOverlapTypeMatch->GetYaxis()->SetBinLabel(1 ,"#gamma");
    fhMCEOverlapTypeMatch->GetYaxis()->SetBinLabel(2 ,"e^{#pm}");
    fhMCEOverlapTypeMatch->GetYaxis()->SetBinLabel(3 ,"hadron^{#pm}");
    fhMCEOverlapTypeMatch->GetYaxis()->SetBinLabel(4 ,"hadron^{0}");
    fhMCEOverlapTypeMatch->GetYaxis()->SetBinLabel(5 ,"??");
    fhMCEOverlapTypeMatch->SetXTitle("Cluster E (GeV)");
    outputContainer->Add(fhMCEOverlapTypeMatch) ;
    
  }// MC analysis, check overlaps
  
  
  if(IsDataMC())
  {
    for(Int_t inlm = 0; inlm < 3; inlm++)
    {
      fhAsyMCGenRecoDiffMCPi0[inlm]  = new TH2F(Form("hAsyMCGenRecoFracNLocMax%sMCPi0",snlm[inlm].Data()),
                                                Form("Reconstructed - Generated asymmetry with NLM=%d vs E, MC Pi0",inlm),
                                                nptbins,ptmin,ptmax,200,-1,1);
      fhAsyMCGenRecoDiffMCPi0[inlm]->SetYTitle("A_{reco} - A_{gen}");
      fhAsyMCGenRecoDiffMCPi0[inlm]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyMCGenRecoDiffMCPi0[inlm]) ;
      
      fhAsyMCGenRecoDiffMCPi0Conv[inlm]  = new TH2F(Form("hAsyMCGenRecoFracNLocMax%sMCPi0Conv",snlm[inlm].Data()),
                                                Form("Reconstructed - Generated asymmetry  with NLM=%d vs E, MC Pi0Conv",inlm),
                                                nptbins,ptmin,ptmax,200,-1,1);
      fhAsyMCGenRecoDiffMCPi0Conv[inlm]->SetYTitle("A_{reco} - A_{gen}");
      fhAsyMCGenRecoDiffMCPi0Conv[inlm]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyMCGenRecoDiffMCPi0Conv[inlm]) ;
    }
  }
  
  if(fFillArmenterosHisto)
  {
    Int_t narmbins = 400;
    Float_t armmin = 0;
    Float_t armmax = 0.4;
    
    for(Int_t i = 0; i < 7; i++) // MC bin
    {
      for(Int_t j = 0; j < 4; j++) // E bin
      {
        fhArmNLocMax1[i][j]  = new TH2F(Form("hArmNLocMax1EBin%d%s",j,pname[i].Data()),
                                         Form("Armenteros of splitted cluster with NLM=1, %s, %s",sEBin[j].Data(),ptype[i].Data()),
                                         200, -1, 1, narmbins,armmin,armmax);
        fhArmNLocMax1[i][j]->SetYTitle("p_{T}^{Arm}");
        fhArmNLocMax1[i][j]->SetXTitle("#alpha^{Arm}");
        outputContainer->Add(fhArmNLocMax1[i][j]) ;
        
        fhArmNLocMax2[i][j]  = new TH2F(Form("hArmNLocMax2EBin%d%s",j,pname[i].Data()),
                                         Form("Armenteros of splitted cluster with NLM=2, %s, %s",sEBin[j].Data(),ptype[i].Data()),
                                         200, -1, 1, narmbins,armmin,armmax);
        fhArmNLocMax2[i][j]->SetYTitle("p_{T}^{Arm}");
        fhArmNLocMax2[i][j]->SetXTitle("#alpha^{Arm}");
        outputContainer->Add(fhArmNLocMax2[i][j]) ;
        
        fhArmNLocMaxN[i][j]  = new TH2F(Form("hArmNLocMaxNEBin%d%s",j,pname[i].Data()),
                                         Form("Armenteros of splitted cluster with NLM>2, %s, %s",sEBin[j].Data(),ptype[i].Data()),
                                         200, -1, 1, narmbins,armmin,armmax);
        fhArmNLocMaxN[i][j]->SetYTitle("p_{T}^{Arm}");
        fhArmNLocMaxN[i][j]->SetXTitle("#alpha^{Arm}");
        outputContainer->Add(fhArmNLocMaxN[i][j]) ;
        
        if(asyOn || m02On)
        {
          fhArmAfterCutsNLocMax1[i][j]  = new TH2F(Form("hArmAfterCutsNLocMax1EBin%d%s",j,pname[i].Data()),
                                          Form("Armenteros of splitted cluster with NLM=1, %s, %s",sEBin[j].Data(),ptype[i].Data()),
                                          200, -1, 1, narmbins,armmin,armmax);
          fhArmAfterCutsNLocMax1[i][j]->SetYTitle("p_{T}^{Arm}");
          fhArmAfterCutsNLocMax1[i][j]->SetXTitle("#alpha^{Arm}");
          outputContainer->Add(fhArmAfterCutsNLocMax1[i][j]) ;
          
          fhArmAfterCutsNLocMax2[i][j]  = new TH2F(Form("hArmAfterCutsNLocMax2EBin%d%s",j,pname[i].Data()),
                                          Form("Armenteros of splitted cluster with NLM=2, %s, %s",sEBin[j].Data(),ptype[i].Data()),
                                          200, -1, 1, narmbins,armmin,armmax);
          fhArmAfterCutsNLocMax2[i][j]->SetYTitle("p_{T}^{Arm}");
          fhArmAfterCutsNLocMax2[i][j]->SetXTitle("#alpha^{Arm}");
          outputContainer->Add(fhArmAfterCutsNLocMax2[i][j]) ;
          
          fhArmAfterCutsNLocMaxN[i][j]  = new TH2F(Form("hArmAfterCutsNLocMaxNEBin%d%s",j,pname[i].Data()),
                                          Form("Armenteros of splitted cluster with NLM>2, %s, %s",sEBin[j].Data(),ptype[i].Data()),
                                          200, -1, 1, narmbins,armmin,armmax);
          fhArmAfterCutsNLocMaxN[i][j]->SetYTitle("p_{T}^{Arm}");
          fhArmAfterCutsNLocMaxN[i][j]->SetXTitle("#alpha^{Arm}");
          outputContainer->Add(fhArmAfterCutsNLocMaxN[i][j]) ;
        }

        fhArmPi0NLocMax1[i][j]  = new TH2F(Form("hArmPi0NLocMax1EBin%d%s",j,pname[i].Data()),
                                                 Form("Armenteros of splitted cluster with NLM=1, %s, %s",sEBin[j].Data(),ptype[i].Data()),
                                                 200, -1, 1, narmbins,armmin,armmax);
        fhArmPi0NLocMax1[i][j]->SetYTitle("p_{T}^{Arm}");
        fhArmPi0NLocMax1[i][j]->SetXTitle("#alpha^{Arm}");
        outputContainer->Add(fhArmPi0NLocMax1[i][j]) ;
        
        fhArmPi0NLocMax2[i][j]  = new TH2F(Form("hArmPi0NLocMax2EBin%d%s",j,pname[i].Data()),
                                                 Form("Armenteros of splitted cluster with NLM=2, %s, %s",sEBin[j].Data(),ptype[i].Data()),
                                                 200, -1, 1, narmbins,armmin,armmax);
        fhArmPi0NLocMax2[i][j]->SetYTitle("p_{T}^{Arm}");
        fhArmPi0NLocMax2[i][j]->SetXTitle("#alpha^{Arm}");
        outputContainer->Add(fhArmPi0NLocMax2[i][j]) ;
        
        fhArmPi0NLocMaxN[i][j]  = new TH2F(Form("hArmPi0NLocMaxNEBin%d%s",j,pname[i].Data()),
                                                 Form("Armenteros of splitted cluster with NLM>2, %s, %s",sEBin[j].Data(),ptype[i].Data()),
                                                 200, -1, 1, narmbins,armmin,armmax);
        fhArmPi0NLocMaxN[i][j]->SetYTitle("p_{T}^{Arm}");
        fhArmPi0NLocMaxN[i][j]->SetXTitle("#alpha^{Arm}");
        outputContainer->Add(fhArmPi0NLocMaxN[i][j]) ;
 
      }
    }
  }
  
  return outputContainer ;
  
}

//_____________________________________________________________________________
void AliAnaInsideClusterInvariantMass::GetMCIndex(AliVCluster* cluster,
                                                  Int_t & mcindex, Int_t & tag)
{
  
  // Assign mc index depending on MC bit set, to be used in histograms arrays
    
  tag	= GetMCAnalysisUtils()->CheckOrigin(cluster->GetLabels(),cluster->GetNLabels(), GetReader());
  
  if      ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) &&
           !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = kmcPi0;
  else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)  )      mcindex = kmcPi0Conv;
  else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)  )      mcindex = kmcEta;
  else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) &&
           !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = kmcPhoton;
  else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) &&
            GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = kmcConversion;
  else if (!GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron)  ) mcindex = kmcHadron;
 
  //printf("MC index %d\n",mcindex);
  
}

//____________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::GetMCPrimaryKine(AliVCluster* cluster, Int_t mcindex,
                                                        Int_t mctag, Bool_t matched,
                                                        Float_t & eprim, Float_t & asymGen,
                                                        Float_t & angleGen, Int_t & noverlaps )
{
  // Check origin of the candidates, get primary kinematics if overlapped meson decay
  
  Bool_t ok      = kFALSE;
  Int_t  mcLabel = cluster->GetLabel();
  
  TLorentzVector primary = GetMCAnalysisUtils()->GetMother(mcLabel,GetReader(),ok);
  eprim = primary.E();
  
  Int_t mesonLabel = -1;
  
  if(mcindex == kmcPi0 || mcindex == kmcEta || mcindex == kmcPi0Conv)
  {
    if(mcindex == kmcPi0 || mcindex == kmcPi0Conv)
    {
      GetMCAnalysisUtils()->GetMCDecayAsymmetryAngleForPDG(mcLabel,111,GetReader(),asymGen,angleGen,ok);
      asymGen = TMath::Abs(asymGen);
      TLorentzVector grandmom = GetMCAnalysisUtils()->GetMotherWithPDG(mcLabel,111,GetReader(),ok,mesonLabel);
      if(grandmom.E() > 0 && ok) eprim =  grandmom.E();
    }
    else
    {
      GetMCAnalysisUtils()->GetMCDecayAsymmetryAngleForPDG(mcLabel,221,GetReader(),asymGen,angleGen,ok);
      asymGen = TMath::Abs(asymGen);
      TLorentzVector grandmom = GetMCAnalysisUtils()->GetMotherWithPDG(mcLabel,221,GetReader(),ok,mesonLabel);
      if(grandmom.E() > 0 && ok) eprim =  grandmom.E();
    }
  }
  
  if(!fFillMCOverlapHisto) return;
    
  const UInt_t nlabels = cluster->GetNLabels();
  Int_t overpdg[nlabels];
  noverlaps = GetMCAnalysisUtils()->GetNOverlaps(cluster->GetLabels(), nlabels,mctag,mesonLabel,GetReader(),overpdg);

  for(Int_t iover = 0; iover < noverlaps; iover++)
  {
    Float_t histobin = -1;
    Int_t   mpdg     = overpdg[iover];
    
    if     (mpdg==22)      histobin = 0.5;
    else if(TMath::Abs(mpdg)==11) histobin = 1.5;
    else if(mpdg==-999999) histobin = 4.5;
    else
    {
      Double_t charge = TDatabasePDG::Instance()->GetParticle(mpdg)->Charge();
      if(TMath::Abs(charge) > 0 ) histobin = 2.5;
      else                        histobin = 3.5;
      //printf("charge %f\n",charge);
    }
    
    //printf("\t pdg = %d, histobin %2.1f\n",mpdg,histobin);
    if(histobin > 0)
    {
      if(matched)fhMCEOverlapType     ->Fill(cluster->E(),histobin);
      else       fhMCEOverlapTypeMatch->Fill(cluster->E(),histobin);
    }
  }
}

//___________________________________________
void AliAnaInsideClusterInvariantMass::Init()
{
  //Init
  //Do some checks
  
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn() && NewOutputAOD())
  {
    AliFatal("!!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD())
  {
    AliFatal("!!STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
  }
  
  if( GetReader()->GetDataType() == AliCaloTrackReader::kMC )
  {
    AliFatal("!!STOP: You want to use pure MC data!!\n");
  }
  
}

//_____________________________________________________
void AliAnaInsideClusterInvariantMass::InitParameters()
{
  //Initialize the parameters of the analysis.  
  AddToHistogramsName("AnaPi0InsideClusterInvariantMass_");
  
  fCalorimeter = "EMCAL" ;

  fMinNCells   = 4 ;
  fMinBadDist  = 2 ;
  
  fHistoECut   = 8 ;
  
  fSSWeightN   = 5;
  fSSWeight[0] = 4.6;  fSSWeight[1] = 4.7; fSSWeight[2] = 4.8; fSSWeight[3] = 4.9; fSSWeight[4] = 5.0;
  fSSWeight[5] = 5.1;  fSSWeight[6] = 5.2; fSSWeight[7] = 5.3; fSSWeight[8] = 5.4; fSSWeight[9] = 5.5;
  
  fSSECellCutN   = 10;
  fSSECellCut[0] = 0.16;  fSSECellCut[1] = 0.18; fSSECellCut[2] = 0.2; fSSECellCut[3] = 0.22; fSSECellCut[4] = 0.24;
  fSSECellCut[5] = 0.26;  fSSECellCut[6] = 0.28; fSSECellCut[7] = 0.3; fSSECellCut[8] = 0.32; fSSECellCut[9] = 0.34;

}


//__________________________________________________________________
void  AliAnaInsideClusterInvariantMass::MakeAnalysisFillHistograms() 
{
  //Search for pi0 in fCalorimeter with shower shape analysis 
  
  TObjArray * pl       = 0x0; 
  AliVCaloCells* cells = 0x0;

  //Select the Calorimeter of the photon
  if(fCalorimeter == "PHOS")
  {
    pl    = GetPHOSClusters();
    cells = GetPHOSCells();
  }
  else if (fCalorimeter == "EMCAL")
  {
    pl    = GetEMCALClusters();
    cells = GetEMCALCells();
  }
  
  if(!pl || !cells) 
  {
    Info("MakeAnalysisFillHistograms","TObjArray with %s clusters is NULL!\n",fCalorimeter.Data());
    return;
  }  
  
	if(fCalorimeter == "PHOS") return; // Not implemented for PHOS yet

  for(Int_t icluster = 0; icluster < pl->GetEntriesFast(); icluster++)
  {
    AliVCluster * cluster = (AliVCluster*) (pl->At(icluster));	

    //-------------------------------------------
    // Get cluster parameters, do some rejection
    //-------------------------------------------
    
    Float_t en = cluster->E();
    Float_t l0 = cluster->GetM02();
    Int_t   nc = cluster->GetNCells();
    Float_t bd = cluster->GetDistanceToBadChannel() ; 
    
    //If too small or big E or low number of cells, or close to a bad channel skip it
    
    if( en < GetMinEnergy() || en > GetMaxEnergy() || nc < fMinNCells || bd < fMinBadDist) continue ;
    
    //printf("en %2.3f  GetMinEnergy() %2.3f, GetMaxEnergy() %2.3f , nc %d fMinNCells %d\n",en, GetMinEnergy(),GetMaxEnergy(),nc ,fMinNCells);
    
    // Track-cluster matching
    
    Bool_t  matched   = IsTrackMatched(cluster,GetReader()->GetInputEvent());
    if(!fFillTMHisto && matched) continue ;

    // Get cluster angles
    
    TLorentzVector lv;
    cluster->GetMomentum(lv, GetVertex(0));
    Float_t eta = lv.Eta();
    Float_t phi = lv.Phi();
    if(phi<0) phi=+TMath::TwoPi();
    
    //printf("en %2.2f, GetMinEnergy() %2.2f, GetMaxEnergy() %2.2f, nc %d, fMinNCells %d,  bd %2.2f, fMinBadDist %2.2f\n",
    //       en,GetMinEnergy(), GetMaxEnergy(), nc, fMinNCells, bd, fMinBadDist);
    
    // Get PID, N local maximum, *** split cluster ***
    
    Int_t    nMax = 0;
    Double_t mass = 0., angle = 0.;
    TLorentzVector    lv1, lv2;
    Int_t    absId1   =-1; Int_t   absId2   =-1;
    Float_t  distbad1 =-1; Float_t distbad2 =-1;
    Bool_t   fidcut1  = 0; Bool_t  fidcut2  = 0;
    
    Int_t pidTag = GetCaloPID()->GetIdentifiedParticleTypeFromClusterSplitting(cluster,cells,GetCaloUtils(),
                                                                               GetVertex(0), nMax, mass, angle,
                                                                               lv1,lv2,absId1,absId2,
                                                                               distbad1,distbad2,fidcut1,fidcut2);
    if (nMax <= 0) 
    {
      if(GetDebug() > 0 )
        Info("MakeAnalysisFillHistograms","No local maximum found! It did not pass CaloPID selection criteria \n");
      
      continue;
    }
    
    // Set some index for array histograms
    
    Int_t inlm = -1;
    if     (nMax == 1) inlm = 0;
    else if(nMax == 2) inlm = 1;
    else if(nMax >  2) inlm = 2;
    else Info("MakeAnalysisFillHistograms","Wrong N local maximum -> %d, n cells in cluster %d \n",nMax,nc);

    // Skip events where one of the new clusters (lowest energy) is close to an EMCal border or a bad channel
    if( (fCheckSplitDistToBad) &&
        (!fidcut2 || !fidcut1 || distbad1 < fMinBadDist || distbad2 < fMinBadDist))
    {
      if(GetDebug() > 1)
        Info("MakeAnalysisFillHistograms","Dist to bad channel cl1 %f, cl2 %f; fid cl1 %d, cl2 %d \n",
                                 distbad1,distbad2, fidcut1,fidcut2);
      
      if(distbad1 < fMinBadDist || distbad2 < fMinBadDist)
      {
        fhMassBadDistClose[inlm]->Fill(en,mass);
        fhM02BadDistClose [inlm]->Fill(en,l0  );
      }
      
      if(!fidcut1 || !fidcut2)
      {
        fhMassOnBorder[inlm]->Fill(en,mass);
        fhM02OnBorder [inlm]->Fill(en,l0  );
      }
      
      continue ;
    }

    // Get sub-cluster parameters
    
    Float_t e1 = lv1.Energy();
    Float_t e2 = lv2.Energy();
    
    Double_t tof1  = cells->GetCellTime(absId1);
    GetCaloUtils()->RecalibrateCellTime(tof1, fCalorimeter, absId1,GetReader()->GetInputEvent()->GetBunchCrossNumber());
    tof1*=1.e9;
    
    Double_t tof2  = cells->GetCellTime(absId2);
    GetCaloUtils()->RecalibrateCellTime(tof2, fCalorimeter, absId2,GetReader()->GetInputEvent()->GetBunchCrossNumber());
    tof2*=1.e9;
    
    Double_t t12diff = tof1-tof2;
    
    Float_t splitFrac = (e1+e2)/en;

    Float_t asym = -10;
    if(e1+e2>0) asym = (e1-e2)/(e1+e2);
    
    //
    
    Int_t ebin = -1;
    if(en > 8  && en <= 12) ebin = 0;
    if(en > 12 && en <= 16) ebin = 1;
    if(en > 16 && en <= 20) ebin = 2;
    if(en > 20)             ebin = 3;
    
    // MC data histograms and some related calculations
    // mc tag, n overlaps, asym of generated mesons
    
    Int_t   mcindex   = -1;
    Int_t   mctag     = -1;
    Float_t eprim     = -1;
    Float_t asymGen   = -2;
    Float_t angleGen  =  2000;
    Int_t   noverlaps =  0;
    
    if(IsDataMC())
    {
      // MC indexes

      GetMCIndex(cluster,mcindex,mctag);
      
      // MC primary kine, generation fractions

      GetMCPrimaryKine(cluster,mcindex,mctag,matched,eprim,asymGen,angleGen,noverlaps);
        
      // For cluster with MC pi0 and more than 1 maxima
      
     }
    
    //
    
    FillHistograms1(en, e1, e2, nMax, mass, l0, eta, phi, matched, mcindex);

    //
    
    if(fFillNCellHisto)
      FillNCellHistograms(nc,en, nMax,matched, mcindex,mass,l0);
    
    //
    
    if(fFillSSExtraHisto)
      FillSSExtraHistograms(cluster, nMax, matched,mcindex,mass,ebin)  ;
    
    //

    if(!matched && ebin >= 0 && fFillEbinHisto)
      FillEBinHistograms(ebin,nMax,mcindex,splitFrac,mass,asym,l0);

    //
    
    if(fFillAngleHisto)
      FillAngleHistograms(nMax,matched,mcindex,en,e1,e2,angle,mass,angleGen,l0, asym,pidTag,noverlaps);

    if(fFillArmenterosHisto && ebin >= 0)
      FillArmenterosHistograms(nMax, ebin, mcindex, en, lv1, lv2, l0, pidTag);

    if(fFillThetaStarHisto)
      FillThetaStarHistograms(nMax,matched,mcindex, en, lv1, lv2, l0, pidTag);

    
    //---------------------------------------------------------------------
    // From here start applying some cuts
    //---------------------------------------------------------------------
    
    // If set, check just on MC clusters when SS cut is applied
    if( IsDataMC() && mcindex > 0 && mcindex < 7 && GetCaloPID()->IsInM02Range(l0) )
    {
      // For cluster with MC pi0 and more than 1 maxima
      
      if(fFillMCOverlapHisto)
        CheckLocalMaximaMCOrigin(cluster, mcindex,noverlaps,e1,e2,mass);
      //l0, l1, l2);
      
      //
      
      if(fFillMCHisto)
        FillMCHistograms(en,e1,e2,ebin,mcindex,noverlaps,l0,mass,
                         nMax,matched,splitFrac, asym, eprim,asymGen);
      
      //
      
      if(fFillMCOverlapHisto)
        FillMCOverlapHistograms(en,eprim,nc,mass,l0,asym,splitFrac,inlm,ebin,matched,mcindex,noverlaps);
      
    }
    
    // Fill few histograms, some still without cuts
    FillHistograms2(en, eprim, e1, e2, nMax, mass, l0, matched, mcindex);
    
    if(pidTag==AliCaloPID::kPi0)
    {
      FillIdPi0Histograms(en, e1, e2, nc, nMax, t12diff, mass, l0, eta, phi, matched, mcindex);
      
      if(fFillSSWeightHisto)
        FillSSWeightHistograms(cluster, inlm, absId1, absId2);
      
      if(fFillTMHisto && fFillTMResidualHisto)
        FillTrackMatchingHistograms(cluster,nMax,mcindex);
      
      if(fFillMCOverlapHisto && IsDataMC() && mcindex > 0 && mcindex < 7 && !matched)
      {
        if     (noverlaps == 0) fhMCEEpriOverlap0IdPi0 [inlm][mcindex]->Fill(en, eprim);
        else if(noverlaps == 1) fhMCEEpriOverlap1IdPi0 [inlm][mcindex]->Fill(en, eprim);
        else if(noverlaps  > 1) fhMCEEpriOverlapNIdPi0 [inlm][mcindex]->Fill(en, eprim);
      }
    }
    else if(fFillIdEtaHisto && pidTag==AliCaloPID::kEta)
    {
      FillIdEtaHistograms(en, e1, e2, nc, nMax, t12diff, mass, l0, eta, phi, matched, mcindex);
    }
    else if(fFillIdConvHisto && pidTag==AliCaloPID::kPhoton)
    {
      FillIdConvHistograms(en, nMax, asym, mass, l0, matched, mcindex);
    }
    
  }//loop
  
  if(GetDebug() > 1) Info("MakeAnalysisFillHistograms","END \n");

}

//______________________________________________________________________
void AliAnaInsideClusterInvariantMass::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print("");
  printf("Calorimeter     =     %s\n",  fCalorimeter.Data()) ;
  if(GetCaloUtils()) printf("Loc. Max. E > %2.2f\n",       GetCaloUtils()->GetLocalMaximaCutE());
  if(GetCaloUtils()) printf("Loc. Max. E Diff > %2.2f\n",  GetCaloUtils()->GetLocalMaximaCutEDiff());
  printf("Min. N Cells =%d \n",         fMinNCells) ;
  printf("Min. Dist. to Bad =%1.1f \n", fMinBadDist) ;
  if(fFillSSWeightHisto) printf(" N w %d - N e cut %d \n",fSSWeightN,fSSECellCutN);

  printf("    \n") ;
  
} 

//___________________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::RecalculateClusterShowerShapeParametersWithCellCut(const AliEMCALGeometry * geom,
                                                                                          AliVCaloCells* cells,
                                                                                          AliVCluster * cluster,
                                                                                          Float_t & l0,   Float_t & l1,
                                                                                          Float_t & disp, Float_t & dEta, Float_t & dPhi,
                                                                                          Float_t & sEta, Float_t & sPhi, Float_t & sEtaPhi,
                                                                                          Float_t eCellMin)
{
  // Calculates new center of gravity in the local EMCAL-module coordinates
  // and tranfers into global ALICE coordinates
  // Calculates Dispersion and main axis
  
  if(!cluster)
  {
    AliInfo("Cluster pointer null!");
    return;
  }
  
  Double_t eCell       = 0.;
  Float_t  fraction    = 1.;
  Float_t  recalFactor = 1.;
  
  Int_t    iSupMod = -1;
  Int_t    iTower  = -1;
  Int_t    iIphi   = -1;
  Int_t    iIeta   = -1;
  Int_t    iphi    = -1;
  Int_t    ieta    = -1;
  Double_t etai    = -1.;
  Double_t phii    = -1.;
  
  Int_t    nstat   = 0 ;
  Float_t  wtot    = 0.;
  Double_t w       = 0.;
  Double_t etaMean = 0.;
  Double_t phiMean = 0.;
    
  //Loop on cells, calculate the cluster energy, in case a cut on cell energy is added
  // and to check if the cluster is between 2 SM in eta
  Int_t   iSM0   = -1;
  Bool_t  shared = kFALSE;
  Float_t energy = 0;

  for(Int_t iDigit=0; iDigit < cluster->GetNCells(); iDigit++)
  {
    //Get from the absid the supermodule, tower and eta/phi numbers
    geom->GetCellIndex(cluster->GetCellAbsId(iDigit),iSupMod,iTower,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
    
    //Check if there are cells of different SM
    if     (iDigit == 0   ) iSM0 = iSupMod;
    else if(iSupMod!= iSM0) shared = kTRUE;
    
    //Get the cell energy, if recalibration is on, apply factors
    fraction  = cluster->GetCellAmplitudeFraction(iDigit);
    if(fraction < 1e-4) fraction = 1.; // in case unfolding is off
    
    if(GetCaloUtils()->GetEMCALRecoUtils()->IsRecalibrationOn())
    {
      recalFactor = GetCaloUtils()->GetEMCALRecoUtils()->GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
    }
    
    eCell  = cells->GetCellAmplitude(cluster->GetCellAbsId(iDigit))*fraction*recalFactor;
    
    if(eCell > eCellMin) energy += eCell;
    
  }//cell loop
  
  //Loop on cells, get weighted parameters
  for(Int_t iDigit=0; iDigit < cluster->GetNCells(); iDigit++)
  {
    //Get from the absid the supermodule, tower and eta/phi numbers
    geom->GetCellIndex(cluster->GetCellAbsId(iDigit),iSupMod,iTower,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
    
    //Get the cell energy, if recalibration is on, apply factors
    fraction  = cluster->GetCellAmplitudeFraction(iDigit);
    if(fraction < 1e-4) fraction = 1.; // in case unfolding is off
    
    if(GetCaloUtils()->GetEMCALRecoUtils()->IsRecalibrationOn())
    {
      recalFactor = GetCaloUtils()->GetEMCALRecoUtils()->GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
    }
    
    eCell  = cells->GetCellAmplitude(cluster->GetCellAbsId(iDigit))*fraction*recalFactor;
    
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2
    // C Side impair SM, nSupMod%2=1; A side pair SM, nSupMod%2=0
    if(shared && iSupMod%2) ieta+=AliEMCALGeoParams::fgkEMCALCols;
    
    if(energy > 0 && eCell > eCellMin)
    {
      w  = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(eCell,energy);

      //correct weight, ONLY in simulation
      w *= (1 - fWSimu * w );

      etai=(Double_t)ieta;
      phii=(Double_t)iphi;
      
      if(w > 0.0)
      {
        wtot += w ;
        nstat++;
        //Shower shape
        sEta     += w * etai * etai ;
        etaMean  += w * etai ;
        sPhi     += w * phii * phii ;
        phiMean  += w * phii ;
        sEtaPhi  += w * etai * phii ;
      }
    }
    else if(energy == 0 || (eCellMin <0.01 && eCell == 0)) AliError(Form("Wrong energy %f and/or amplitude %f\n", eCell, energy));
    
  }//cell loop
  
  //Normalize to the weight
  if (wtot > 0)
  {
    etaMean /= wtot ;
    phiMean /= wtot ;
  }
  else
    AliError(Form("Wrong weight %f\n", wtot));
  
  //Calculate dispersion
  for(Int_t iDigit=0; iDigit < cluster->GetNCells(); iDigit++)
  {
    //Get from the absid the supermodule, tower and eta/phi numbers
    geom->GetCellIndex(cluster->GetCellAbsId(iDigit),iSupMod,iTower,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
    
    //Get the cell energy, if recalibration is on, apply factors
    fraction  = cluster->GetCellAmplitudeFraction(iDigit);
    if(fraction < 1e-4) fraction = 1.; // in case unfolding is off
    if (GetCaloUtils()->GetEMCALRecoUtils()->IsRecalibrationOn())
    {
      recalFactor = GetCaloUtils()->GetEMCALRecoUtils()->GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
    }
    
    eCell  = cells->GetCellAmplitude(cluster->GetCellAbsId(iDigit))*fraction*recalFactor;
    
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2
    // C Side impair SM, nSupMod%2=1; A side pair SM, nSupMod%2=0
    if(shared && iSupMod%2) ieta+=AliEMCALGeoParams::fgkEMCALCols;
    
    if(energy > 0 && eCell > eCellMin)
    {
      w  = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(eCell,energy);
      
      //correct weight, ONLY in simulation
      w *= (1 - fWSimu * w );

      etai=(Double_t)ieta;
      phii=(Double_t)iphi;
      if(w > 0.0)
      {
        disp +=  w *((etai-etaMean)*(etai-etaMean)+(phii-phiMean)*(phii-phiMean));
        dEta +=  w * (etai-etaMean)*(etai-etaMean) ;
        dPhi +=  w * (phii-phiMean)*(phii-phiMean) ;
      }
    }
    else if(energy == 0 || (eCellMin <0.01 && eCell == 0)) AliError(Form("Wrong energy %f and/or amplitude %f\n", eCell, energy));
  }// cell loop
  
  //Normalize to the weigth and set shower shape parameters
  if (wtot > 0 && nstat > 1)
  {
    disp    /= wtot ;
    dEta    /= wtot ;
    dPhi    /= wtot ;
    sEta    /= wtot ;
    sPhi    /= wtot ;
    sEtaPhi /= wtot ;
    
    sEta    -= etaMean * etaMean ;
    sPhi    -= phiMean * phiMean ;
    sEtaPhi -= etaMean * phiMean ;
    
    l0 = (0.5 * (sEta + sPhi) + TMath::Sqrt( 0.25 * (sEta - sPhi) * (sEta - sPhi) + sEtaPhi * sEtaPhi ));
    l1 = (0.5 * (sEta + sPhi) - TMath::Sqrt( 0.25 * (sEta - sPhi) * (sEta - sPhi) + sEtaPhi * sEtaPhi ));
  }
  else
  {
    l0   = 0. ;
    l1   = 0. ;
    dEta = 0. ; dPhi = 0. ; disp    = 0. ;
    sEta = 0. ; sPhi = 0. ; sEtaPhi = 0. ;
  }
  
}


