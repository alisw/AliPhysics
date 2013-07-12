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
  AliAnaCaloTrackCorrBaseClass(),  
  fCalorimeter(""), 
  fM02MaxCut(0),    fM02MinCut(0),       
  fMinNCells(0),    fMinBadDist(0),
  fHistoECut(0),
  fFillAngleHisto(kFALSE),
  fFillTMHisto(kFALSE),
  fFillTMResidualHisto(kFALSE),
  fFillSSExtraHisto(kFALSE),
  fFillMCHisto(kFALSE),
  fFillSSWeightHisto(kFALSE),
  fFillEbinHisto(0),
  fFillMCOverlapHisto(0),
  fSSWeightN(0),              fSSECellCutN(0),
  fWSimu(0),
  fhMassM02CutNLocMax1(0),    fhMassM02CutNLocMax2(0),    fhMassM02CutNLocMaxN(0),
  fhAsymM02CutNLocMax1(0),    fhAsymM02CutNLocMax2(0),    fhAsymM02CutNLocMaxN(0),
  fhMassSplitECutNLocMax1(0), fhMassSplitECutNLocMax2(0), fhMassSplitECutNLocMaxN(0),
  fhMCGenFracAfterCutsNLocMax1MCPi0(0),
  fhMCGenFracAfterCutsNLocMax2MCPi0(0),
  fhMCGenFracAfterCutsNLocMaxNMCPi0(0),
  fhMCGenSplitEFracAfterCutsNLocMax1MCPi0(0),
  fhMCGenSplitEFracAfterCutsNLocMax2MCPi0(0),
  fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0(0),
  fhEventPlanePi0NLocMax1(0), fhEventPlaneEtaNLocMax1(0),
  fhEventPlanePi0NLocMax2(0), fhEventPlaneEtaNLocMax2(0),
  fhEventPlanePi0NLocMaxN(0), fhEventPlaneEtaNLocMaxN(0),
  fhClusterEtaPhiNLocMax1(0), fhClusterEtaPhiNLocMax2(0),  fhClusterEtaPhiNLocMaxN(0),
  fhPi0EtaPhiNLocMax1(0),     fhPi0EtaPhiNLocMax2(0),      fhPi0EtaPhiNLocMaxN(0),
  fhEtaEtaPhiNLocMax1(0),     fhEtaEtaPhiNLocMax2(0),      fhEtaEtaPhiNLocMaxN(0),
  fhPi0EPairDiffTimeNLM1(0),  fhPi0EPairDiffTimeNLM2(0),   fhPi0EPairDiffTimeNLMN(0),
  fhEtaEPairDiffTimeNLM1(0),  fhEtaEPairDiffTimeNLM2(0),   fhEtaEPairDiffTimeNLMN(0),
  fhMCPi0HighNLMPair(0),      fhMCPi0LowNLMPair(0),
  fhMCPi0AnyNLMPair(0),       fhMCPi0NoneNLMPair(0),
  fhMCPi0HighNLMPairNoMCMatch(0), fhMCPi0LowNLMPairNoMCMatch(0),
  fhMCPi0AnyNLMPairNoMCMatch(0),  fhMCPi0NoneNLMPairNoMCMatch(0),
  fhMCEOverlapType(0),            fhMCEOverlapTypeMatch(0)
{
  //default ctor
  
  // Init array of histograms
  for(Int_t i = 0; i < 8; i++)
  {
    for(Int_t j = 0; j < 2; j++)
    {
      fhMassNLocMax1[i][j]  = 0;
      fhMassNLocMax2[i][j]  = 0;
      fhMassNLocMaxN[i][j]  = 0;
      fhNLocMax[i][j]       = 0;
      fhNLocMaxM02Cut[i][j] = 0;
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
      
      fhMCGenFracNLocMax1[i][j]= 0;
      fhMCGenFracNLocMax2[i][j]= 0;
      fhMCGenFracNLocMaxN[i][j]= 0;
      
      fhMCGenSplitEFracNLocMax1[i][j]= 0;
      fhMCGenSplitEFracNLocMax2[i][j]= 0;
      fhMCGenSplitEFracNLocMaxN[i][j]= 0;    
      
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
      
      fhCentralityPi0NLocMax1[i][j] = 0 ;
      fhCentralityEtaNLocMax1[i][j] = 0 ;

      fhCentralityPi0NLocMax2[i][j] = 0 ;
      fhCentralityEtaNLocMax2[i][j] = 0 ;

      fhCentralityPi0NLocMaxN[i][j] = 0 ;
      fhCentralityEtaNLocMaxN[i][j] = 0 ;      
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
    fhAnglePairNLocMax1    [i] = 0;
    fhAnglePairNLocMax2    [i] = 0;
    fhAnglePairNLocMaxN    [i] = 0;
    fhAnglePairMassNLocMax1[i] = 0;
    fhAnglePairMassNLocMax2[i] = 0;
    fhAnglePairMassNLocMaxN[i] = 0;
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
  }
  
  InitParameters();

}

//_______________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::CheckLocalMaximaMCOrigin(AliVCluster* cluster, const Int_t mcindex)
{
  // Check origin NLM tower of the cluster, when MC gives merged pi0
  
  if(!fFillMCOverlapHisto) return;

  if(!IsDataMC()) return;
  
  if(mcindex != kmcPi0 && mcindex != kmcPi0Conv) return;
  
  const UInt_t nc = cluster->GetNCells();
  Int_t   list[nc];
  Float_t elist[nc];
  Int_t nMax = GetCaloUtils()->GetNumberOfLocalMaxima(cluster, GetEMCALCells(),list, elist);
    
  if(nMax < 2) return;
  
//  printf("AliAnaInsideClusterInvariantMass::CheckLocalMaximaMCOrigin() - Cluster E %2.2f; NLM = %d, cluster MC labels:\n",cluster->E(),nMax);
//  
//  for (UInt_t ilab = 0; ilab < cluster->GetNLabels(); ilab++ )
//  {
//    Bool_t ok  =kFALSE,gok = kFALSE;
//    Int_t pdg    = -22222, status   = -1;
//    Int_t gpdg   = -22222, gstatus  = -1;
//    Int_t ggpdg  = -22222, ggstatus = -1;
//    Int_t gLabel = -1, ggLabel = -1;
//
//    Int_t label = cluster->GetLabels()[ilab];
//    TLorentzVector primary   =GetMCAnalysisUtils()->GetMother     (label,GetReader(),  pdg,  status, ok);
//    TLorentzVector gprimary  =GetMCAnalysisUtils()->GetGrandMother(label,GetReader(), gpdg, gstatus,gok, gLabel,ggLabel);
//    TLorentzVector ggprimary =GetMCAnalysisUtils()->GetMother(ggLabel  ,GetReader(),ggpdg,ggstatus,gok);
//    printf("\t %d; mother: Label %d; PDG %d; E %2.2f - grand mother label %d; PDG %d; E %2.2f- great grand mother label %d; PDG %d; E %2.2f\n",
//           ilab,label,pdg,primary.E(), gLabel,gpdg,gprimary.E(), ggLabel,ggpdg,ggprimary.E());
//
//  }
//  
//  printf("AliAnaInsideClusterInvariantMass::CheckLocalMaximaMCOrigin() - Cluster Cells MC labels:\n");
//  
  for (UInt_t icell = 0; icell < nc; icell++ )
  {
//    Bool_t ok  =kFALSE,gok = kFALSE;
//    Int_t pdg    = -22222, status   = -1;
//    Int_t gpdg   = -22222, gstatus  = -1;
//    Int_t ggpdg  = -22222, ggstatus = -1;
//    Int_t gLabel = -1, ggLabel = -1;
//    
//    Int_t absId = cluster->GetCellAbsId(icell);
//    
//    printf("cell abs Id %d, amplitude %f\n",absId,GetEMCALCells()->GetCellAmplitude(absId));
//    
//    Int_t label = GetEMCALCells()->GetCellMCLabel(absId);
//    TLorentzVector primary   =GetMCAnalysisUtils()->GetMother     (label,GetReader(),  pdg,  status, ok);
//    TLorentzVector gprimary  =GetMCAnalysisUtils()->GetGrandMother(label,GetReader(), gpdg, gstatus,gok, gLabel,ggLabel);
//    TLorentzVector ggprimary =GetMCAnalysisUtils()->GetMother(ggLabel  ,GetReader(),ggpdg,ggstatus,gok);
//    printf(" %d; mother: Label %d; PDG %d; E %2.2f - grand mother label %d; PDG %d; E %2.2f- great grand mother label %d; PDG %d; E %2.2f\n",
//           icell,label,pdg,primary.E(), gLabel,gpdg,gprimary.E(), ggLabel,ggpdg,ggprimary.E());
    
  }
  
  //Find highest energy Local Maxima Towers
  Int_t   imax  = -1;
  Int_t   imax2 = -1;
  Float_t emax  = -1;
  Float_t emax2 = -1;
  for(Int_t i = 0; i < nMax; i++)
  {
    //printf("i %d: AbsId %d; E %2.3f\n",i,list[i],elist[i]);
    if(elist[i] > emax)
    {
      imax = i;
      emax = elist[i];
    }
  }
  //Find second highest
  for(Int_t i = 0; i < nMax; i++)
  {
    if(i==imax) continue;
    if(elist[i] > emax2)
    {
      imax2 = i;
      emax2 = elist[i];
    }
  }
  
  //printf("Highest : %d and %d\n",imax,imax2);
  
  // Check that the highest mc label and the max cluster label are the same
  Int_t mcLabelMax = GetEMCALCells()->GetCellMCLabel(list[imax]);
  GetReader()->RemapMCLabelForAODs(mcLabelMax);
  Int_t mcLabelMax2 = GetEMCALCells()->GetCellMCLabel(list[imax2]);
  GetReader()->RemapMCLabelForAODs(mcLabelMax2);
  
  Int_t mcLabelclusterMax = cluster->GetLabels()[0];
  Bool_t matchHighLMAndHighMC = kFALSE;
  
  if(mcLabelclusterMax==mcLabelMax)
  {
    matchHighLMAndHighMC = kTRUE;
    //printf("*** MATCH cluster and LM maximum ***\n");
  }
  else
  {
     //printf("*** NO MATCH cluster and LM maximum, check second ***\n");
    if(mcLabelclusterMax==mcLabelMax2)
    {
      //printf("\t *** MATCH cluster and 2nd LM maximum ***\n");
      matchHighLMAndHighMC = kTRUE;
    }
    else
    {
      //printf("\t *** NO MATCH***\n");
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
  
  // Compare ancestors of all local maxima
  for(Int_t i = 0; i < nMax-1; i++)
  {
    Int_t mcLabel1 = GetEMCALCells()->GetCellMCLabel(list[i]);
    GetReader()->RemapMCLabelForAODs(mcLabel1);
 
    for(Int_t j = i+1; j < nMax; j++)
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
  
  //printf("Match MC? %d; high %d; low %d\n",matchHighLMAndHighMC,high,low);
  
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

//___________________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillAngleHistograms(const Int_t nMax, const Bool_t matched,
                                                           const Float_t en, const Float_t angle, const Float_t mass)
{
  // Fill histograms related to opening angle
  
  if(!fFillAngleHisto) return;
  
  if     (nMax==1)
  {
    fhAnglePairNLocMax1[matched]->Fill(en,angle);
    if( en > fHistoECut ) fhAnglePairMassNLocMax1[matched]->Fill(mass,angle);
  }
  else if(nMax==2)
  {
    fhAnglePairNLocMax2[matched]->Fill(en,angle);
    if( en > fHistoECut ) fhAnglePairMassNLocMax2[matched]->Fill(mass,angle);
  }
  else if(nMax >2)
  {
    fhAnglePairNLocMaxN[matched]->Fill(en,angle);
    if( en > fHistoECut ) fhAnglePairMassNLocMaxN[matched]->Fill(mass,angle);
  }
  
}

//__________________________________________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillEBinHistograms(const Int_t   ebin     , const Int_t   nMax, const Int_t mcindex,
                                                          const Float_t splitFrac, const Float_t mass, const Float_t asym, const Float_t l0)
{
  // Fill some histograms integrating in few energy bins
  
  if(ebin < 0 || !fFillEbinHisto) return ;
  
  if     (nMax==1)
  {
    fhMassSplitEFractionNLocMax1Ebin[0][ebin]->Fill(splitFrac,  mass);
    if(IsDataMC())fhMassSplitEFractionNLocMax1Ebin[mcindex][ebin]->Fill(splitFrac,  mass);
    
    fhMassM02NLocMax1Ebin    [ebin]->Fill(l0  ,  mass );
    fhMassAsyNLocMax1Ebin    [ebin]->Fill(asym,  mass );
  }
  else if(nMax==2)
  {
    fhMassSplitEFractionNLocMax2Ebin[0][ebin]->Fill(splitFrac,  mass);
    if(IsDataMC())fhMassSplitEFractionNLocMax2Ebin[mcindex][ebin]->Fill(splitFrac,  mass);
    
    fhMassM02NLocMax2Ebin    [ebin]->Fill(l0  ,  mass );
    fhMassAsyNLocMax2Ebin    [ebin]->Fill(asym,  mass );
  }
  else if(nMax > 2 )
  {
    fhMassSplitEFractionNLocMaxNEbin[0][ebin]->Fill(splitFrac,  mass);
    if(IsDataMC())fhMassSplitEFractionNLocMaxNEbin[mcindex][ebin]->Fill(splitFrac,  mass);
    
    fhMassM02NLocMaxNEbin    [ebin]->Fill(l0  ,  mass );
    fhMassAsyNLocMaxNEbin    [ebin]->Fill(asym,  mass );
  }
  
}

//_____________________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillMCHistograms(const Float_t en,        const Float_t e1  , const Float_t e2,
                                                        const Int_t ebin,        const Int_t mcindex,
                                                        const Float_t l0,        const Float_t mass,
                                                        const Int_t nMax,        const Bool_t  matched,
                                                        const Float_t splitFrac, const Float_t asym,
                                                        const Float_t eprim,     const Float_t asymGen)
{
  // Fill histograms needing some MC input
  
  if(!IsDataMC() || !fFillMCHisto) return;
  
  Float_t efrac      = eprim/en;
  Float_t efracSplit = 0;
  if(e1+e2 > 0) efracSplit = eprim/(e1+e2);

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
    
    if( en > fHistoECut )
    {
      fhMCGenEFracvsSplitEFracNLocMax1[mcindex][matched]->Fill(efrac,splitFrac );
      
      if(!matched && ebin >= 0 && fFillEbinHisto)
      {
        fhM02MCGenFracNLocMax1Ebin [mcindex][ebin]->Fill(efrac  ,  l0    );
        fhMassMCGenFracNLocMax1Ebin[mcindex][ebin]->Fill(efrac  ,  mass  );
        
        fhMCAsymM02NLocMax1MCPi0Ebin [ebin]->Fill(l0  ,  asymGen );
        fhAsyMCGenRecoNLocMax1EbinPi0[ebin]->Fill(asym,  asymGen );
      }
    }
  }
  else if(nMax==2)
  {
    fhMCGenFracNLocMax2      [mcindex][matched]->Fill(en     ,  efrac );
    fhMCGenSplitEFracNLocMax2[mcindex][matched]->Fill(en     ,  efracSplit );
    fhMCGenEvsSplitENLocMax2 [mcindex][matched]->Fill(eprim  ,  e1+e2);
    
    if( en > fHistoECut )
    {
      fhMCGenEFracvsSplitEFracNLocMax2[mcindex][matched]->Fill(efrac,splitFrac );
      
      if(!matched && ebin >= 0 && fFillEbinHisto)
      {
        fhM02MCGenFracNLocMax2Ebin [mcindex][ebin]->Fill(efrac  ,  l0    );
        fhMassMCGenFracNLocMax2Ebin[mcindex][ebin]->Fill(efrac  ,  mass  );
        
        fhMCAsymM02NLocMax2MCPi0Ebin [ebin]->Fill(l0  ,  asymGen );
        fhAsyMCGenRecoNLocMax2EbinPi0[ebin]->Fill(asym,  asymGen );
      }
    }

  }
  else if(nMax > 2 )
  {
    fhMCGenFracNLocMaxN      [mcindex][matched]->Fill(en     ,  efrac );
    fhMCGenSplitEFracNLocMaxN[mcindex][matched]->Fill(en     ,  efracSplit );
    fhMCGenEvsSplitENLocMaxN [mcindex][matched]->Fill(eprim  ,  e1+e2);
    
    if( en > fHistoECut )
    {
      fhMCGenEFracvsSplitEFracNLocMaxN[mcindex][matched]->Fill(efrac,splitFrac );
      
      if(!matched && ebin >= 0 && fFillEbinHisto)
      {
        fhM02MCGenFracNLocMaxNEbin [mcindex][ebin]->Fill(efrac  ,  l0    );
        fhMassMCGenFracNLocMaxNEbin[mcindex][ebin]->Fill(efrac  ,  mass  );
        
        fhMCAsymM02NLocMaxNMCPi0Ebin [ebin]->Fill(l0  ,  asymGen );
        fhAsyMCGenRecoNLocMaxNEbinPi0[ebin]->Fill(asym,  asymGen );
      }
    }
  }
}

//_________________________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillMCOverlapHistograms(const Float_t en,      const Float_t mass, const Float_t l0,
                                                               const Int_t   inlm,    const Int_t ebin, const Bool_t matched,
                                                               const Int_t   mcindex, const Int_t noverlaps)
{
  
  // Fill histograms for MC Overlaps
  
  //printf("en %f,mass %f,l0 %f,inlm %d,ebin %d,matched %d,mcindex %d,noverlaps %d \n",en,mass,l0,inlm,ebin,matched,mcindex,noverlaps);
  
  if(!fFillMCOverlapHisto || !IsDataMC()) return;
  
  //printf("AliAnaInsideClusterInvariantMass::FillMCOverlapHistograms - NLM bin=%d, mcIndex %d, n Overlaps %d\n",inlm,mcindex,noverlaps);
  
  if(!matched)
  {
    fhMCENOverlaps[inlm][mcindex]->Fill(en,noverlaps);
    
    if     (noverlaps == 0)
    {
      fhMCEM02Overlap0[inlm][mcindex]->Fill(en, l0);
      fhMCEMassOverlap0[inlm][mcindex]->Fill(en, mass);
      if((mcindex==kmcPi0 || mcindex == kmcPi0Conv) && ebin >=0) fhMCPi0MassM02Overlap0[inlm][ebin]->Fill(l0,mass);
    }
    else if(noverlaps == 1)
    {
      fhMCEM02Overlap1[inlm][mcindex]->Fill(en, l0);
      fhMCEMassOverlap1[inlm][mcindex]->Fill(en, mass);
      if((mcindex==kmcPi0 || mcindex == kmcPi0Conv) && ebin >=0) fhMCPi0MassM02Overlap1[inlm][ebin]->Fill(l0,mass);
    }
    else if(noverlaps  > 1)
    {
      fhMCEM02OverlapN[inlm][mcindex]->Fill(en, l0);
      fhMCEMassOverlapN[inlm][mcindex]->Fill(en, mass);
      if((mcindex==kmcPi0 || mcindex == kmcPi0Conv) && ebin >=0) fhMCPi0MassM02OverlapN[inlm][ebin]->Fill(l0,mass);
    }
    else
      printf("AliAnaInsideClusterInvariantMass::FillMCOverlapHistograms() - n overlaps = %d!!", noverlaps);
  }
  else if(fFillTMHisto)
  {
    fhMCENOverlapsMatch[inlm][mcindex]->Fill(en,noverlaps);
    
    if     (noverlaps == 0)
    {
      fhMCEM02Overlap0Match[inlm][mcindex]->Fill(en, l0);
      fhMCEMassOverlap0Match[inlm][mcindex]->Fill(en, mass);
      if((mcindex==kmcPi0 || mcindex == kmcPi0Conv) && ebin >=0) fhMCPi0MassM02Overlap0Match[inlm][ebin]->Fill(l0,mass);
    }
    else if(noverlaps == 1)
    {
      fhMCEM02Overlap1Match[inlm][mcindex]->Fill(en, l0);
      fhMCEMassOverlap1Match[inlm][mcindex]->Fill(en, mass);
      if((mcindex==kmcPi0 || mcindex == kmcPi0Conv) && ebin >=0) fhMCPi0MassM02Overlap1Match[inlm][ebin]->Fill(l0,mass);
    }
    else if(noverlaps  > 1)
    {
      fhMCEM02OverlapNMatch[inlm][mcindex]->Fill(en, l0);
      fhMCEMassOverlapNMatch[inlm][mcindex]->Fill(en, mass);
      if((mcindex==kmcPi0 || mcindex == kmcPi0Conv) && ebin >=0) fhMCPi0MassM02OverlapNMatch[inlm][ebin]->Fill(l0,mass);
    }
    else
        printf("AliAnaInsideClusterInvariantMass::FillMCOverlapHistograms() - n overlaps in matched = %d!!", noverlaps);
  }
}

//______________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillSSExtraHistograms(AliVCluster  *cluster, const Int_t nMax,
                                                             const Bool_t  matched, const Int_t mcindex,
                                                             const Float_t mass   , const Int_t ebin)
{
  // Fill optional histograms with more SS parameters
  
  if(!fFillSSExtraHisto) return ;
  
  Float_t en = cluster->E();
  Float_t nc = cluster->GetNCells();
  
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
    fhNCellNLocMax1[0][matched]->Fill(en,nc) ;
    if(mcindex > 0 )  fhNCellNLocMax1[mcindex][matched]->Fill(en,nc) ;
    
    if( en > fHistoECut )
    {
      fhMassDispEtaNLocMax1[0][matched]->Fill(dispEta,  mass );
      fhMassDispPhiNLocMax1[0][matched]->Fill(dispPhi,  mass );
      fhMassDispAsyNLocMax1[0][matched]->Fill(dispAsy,  mass );
      
      if(IsDataMC())
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
    fhNCellNLocMax2[0][matched]->Fill(en,nc) ;
    if(mcindex > 0 )  fhNCellNLocMax2[mcindex][matched]->Fill(en,nc) ;
    
    if( en > fHistoECut )
    {
      fhMassDispEtaNLocMax2[0][matched]->Fill(dispEta,  mass );
      fhMassDispPhiNLocMax2[0][matched]->Fill(dispPhi,  mass );
      fhMassDispAsyNLocMax2[0][matched]->Fill(dispAsy,  mass );
      
      if(IsDataMC())
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
    fhNCellNLocMaxN[0][matched]->Fill(en,nc) ;
    if(mcindex > 0 )  fhNCellNLocMaxN[mcindex][matched]->Fill(en,nc) ;
    
    if( en > fHistoECut )
    {
      fhMassDispEtaNLocMaxN[0][matched]->Fill(dispEta,  mass );
      fhMassDispPhiNLocMaxN[0][matched]->Fill(dispPhi,  mass );
      fhMassDispAsyNLocMaxN[0][matched]->Fill(dispAsy,  mass );
      
      if(IsDataMC())
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

//__________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillSSWeightHistograms(AliVCluster *clus,  const Int_t nlm,
                                                              const Int_t absId1, const Int_t absId2)
{
  // Calculate weights and fill histograms
  
  if(!fFillSSWeightHisto) return;
  
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
    printf("AliAnaInsideClusterInvatiantMass::WeightHistograms()- Wrong calculated energy %f\n",energy);
    return;
  }
  
  //Get amplitude of  main local maxima, recalibrate if needed
  Float_t amp1 = cells->GetCellAmplitude(absId1);
  GetCaloUtils()->RecalibrateCellAmplitude(amp1,fCalorimeter, absId1);
  Float_t amp2 = cells->GetCellAmplitude(absId2);
  GetCaloUtils()->RecalibrateCellAmplitude(amp2,fCalorimeter, absId2);

  if(amp1 < amp2)        printf("Bad local maxima E ordering : id1 E %f, id2 E %f\n ",amp1,amp2);
  if(amp1==0 || amp2==0) printf("Null E local maxima : id1 E %f, id2 E %f\n "        ,amp1,amp2);
  
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

//________________________________________________________________________________________
void  AliAnaInsideClusterInvariantMass::FillTrackMatchingHistograms(AliVCluster * cluster,
                                                                    const Int_t nMax,
                                                                    const Int_t mcindex)
{
  // Fill histograms related to track matching
  
  if(!fFillTMResidualHisto) return;
  
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
    
    if(IsDataMC())
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
        
        if(IsDataMC())
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
        
        if(IsDataMC())
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
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaInsideClusterInvariantMass ---\n") ;
  parList+=onePar ;	
  
  snprintf(onePar,buffersize,"Calorimeter: %s\n",        fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fNLocMaxCutE =%2.2f \n",    GetCaloUtils()->GetLocalMaximaCutE()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fNLocMaxCutEDiff =%2.2f \n",GetCaloUtils()->GetLocalMaximaCutEDiff()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"%2.2f< M02 < %2.2f \n",    fM02MinCut, fM02MaxCut) ;
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
  
  TString ptype[] ={"","#gamma","#gamma->e^{#pm}","#pi^{0}","#eta","e^{#pm}", "hadron","#pi^{0} (#gamma->e^{#pm})"}; 
  TString pname[] ={"","Photon","Conversion",     "Pi0",    "Eta", "Electron","Hadron","Pi0Conv"};
  
  Int_t n = 1;
  
  if(IsDataMC()) n = 8;
  
  Int_t nMaxBins = 10;
  
  TString sMatched[] = {"","Matched"};
  
  Int_t nMatched = 2;
  if(!fFillTMHisto) nMatched = 1;
  
  for(Int_t i = 0; i < n; i++)
  {
    for(Int_t j = 0; j < nMatched; j++)
    {
      
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
     
      if(i==0 && j==0)
      {
        fhMassSplitECutNLocMax1  = new TH2F("hMassSplitECutNLocMax1","Invariant mass of splitted cluster with NLM=1 vs E, (E1+E2)/E cut, no TM",
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassSplitECutNLocMax1->SetYTitle("M (GeV/c^{2})");
        fhMassSplitECutNLocMax1->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassSplitECutNLocMax1) ;
        
        fhMassSplitECutNLocMax2  = new TH2F("hMassSplitECutNLocMax2","Invariant mass of splitted cluster with NLM=2 vs E, (E1+E2)/E cut, no TM",
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassSplitECutNLocMax2->SetYTitle("M (GeV/c^{2})");
        fhMassSplitECutNLocMax2->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassSplitECutNLocMax2) ;
        
        fhMassSplitECutNLocMaxN  = new TH2F("hMassSplitECutNLocMaxN","Invariant mass of splitted cluster with NLM>2 vs E, (E1+E2)/E cut, no TM",
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassSplitECutNLocMaxN->SetYTitle("M (GeV/c^{2})");
        fhMassSplitECutNLocMaxN->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassSplitECutNLocMaxN) ;
        
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
        
      }
      
      fhMassAfterCutsNLocMax1[i][j]     = new TH2F(Form("hMassAfterCutsNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                   Form("Mass vs E, %s %s, for N Local max = 1, M02 and asy cut",ptype[i].Data(),sMatched[j].Data()),
                                                   nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassAfterCutsNLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassAfterCutsNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassAfterCutsNLocMax1[i][j]) ;
      
      fhMassAfterCutsNLocMax2[i][j]     = new TH2F(Form("hMassAfterCutsNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                   Form("Mass vs E, %s %s, for N Local max = 2, M02 and asy cut",ptype[i].Data(),sMatched[j].Data()),
                                                   nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassAfterCutsNLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassAfterCutsNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassAfterCutsNLocMax2[i][j]) ;
      
      
      fhMassAfterCutsNLocMaxN[i][j]     = new TH2F(Form("hMassAfterCutsNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                   Form("Mass vs E, %s %s, for N Local max > 2, M02 and asy cut",ptype[i].Data(),sMatched[j].Data()),
                                                   nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassAfterCutsNLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassAfterCutsNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassAfterCutsNLocMaxN[i][j]) ;
            
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
      
      if(i==0 && j==0)
      {
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
      }
      
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
      
      fhNLocMax[i][j]     = new TH2F(Form("hNLocMax%s%s",pname[i].Data(),sMatched[j].Data()),
                                     Form("Number of local maxima in cluster %s %s",ptype[i].Data(),sMatched[j].Data()),
                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins); 
      fhNLocMax[i][j]   ->SetYTitle("N maxima");
      fhNLocMax[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhNLocMax[i][j]) ; 
            
      fhNLocMaxM02Cut[i][j] = new TH2F(Form("hNLocMaxM02Cut%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Number of local maxima in cluster %s for %2.2f < M02 < %2.2f",ptype[i].Data(),fM02MinCut,fM02MaxCut),
                                       nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins); 
      fhNLocMaxM02Cut[i][j]->SetYTitle("N maxima");
      fhNLocMaxM02Cut[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhNLocMaxM02Cut[i][j]) ; 
      
      
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
      
      if(fFillSSExtraHisto)
      {
        fhNCellNLocMax1[i][j]  = new TH2F(Form("hNCellNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                          nptbins,ptmin,ptmax,ncbins,ncmin,ncmax); 
        fhNCellNLocMax1[i][j] ->SetYTitle("N cells");
        fhNCellNLocMax1[i][j] ->SetXTitle("E (GeV)");
        outputContainer->Add(fhNCellNLocMax1[i][j]) ; 
        
        fhNCellNLocMax2[i][j]     = new TH2F(Form("hNCellNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("#lambda_{0}^{2} vs E for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                             nptbins,ptmin,ptmax,ncbins,ncmin,ncmax); 
        fhNCellNLocMax2[i][j]   ->SetYTitle("N cells");
        fhNCellNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhNCellNLocMax2[i][j]) ; 
        
        
        fhNCellNLocMaxN[i][j]     = new TH2F(Form("hNCellNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("#lambda_{0}^{2} vs E for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                             nptbins,ptmin,ptmax,ncbins,ncmin,ncmax); 
        fhNCellNLocMaxN[i][j]   ->SetYTitle("N cells");
        fhNCellNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhNCellNLocMaxN[i][j]) ;
      }
      
      
      // E vs centrality
      
      fhCentralityPi0NLocMax1[i][j]  = new TH2F(Form("hCentralityPi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("E vs Centrality, selected pi0 cluster with NLM=1, %s",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,100,0,100);
      fhCentralityPi0NLocMax1[i][j]->SetYTitle("Centrality");
      fhCentralityPi0NLocMax1[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityPi0NLocMax1[i][j]) ;
      
      fhCentralityPi0NLocMax2[i][j]  = new TH2F(Form("hCentralityPi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("E vs Centrality, selected pi0 cluster with NLM=2, %s",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,100,0,100);
      fhCentralityPi0NLocMax2[i][j]->SetYTitle("Centrality");
      fhCentralityPi0NLocMax2[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityPi0NLocMax2[i][j]) ;
      
      fhCentralityPi0NLocMaxN[i][j]  = new TH2F(Form("hCentralityPi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("E vs Centrality, selected pi0 cluster with NLM>1, %s",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,100,0,100);
      fhCentralityPi0NLocMaxN[i][j]->SetYTitle("Centrality");
      fhCentralityPi0NLocMaxN[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityPi0NLocMaxN[i][j]) ;
      
      fhCentralityEtaNLocMax1[i][j]  = new TH2F(Form("hCentralityEtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("E vs Centrality, selected pi0 cluster with NLM=1, %s",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,100,0,100);
      fhCentralityEtaNLocMax1[i][j]->SetYTitle("Centrality");
      fhCentralityEtaNLocMax1[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityEtaNLocMax1[i][j]) ;
      
      fhCentralityEtaNLocMax2[i][j]  = new TH2F(Form("hCentralityEtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("E vs Centrality, selected pi0 cluster with NLM=2, %s",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,100,0,100);
      fhCentralityEtaNLocMax2[i][j]->SetYTitle("Centrality");
      fhCentralityEtaNLocMax2[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityEtaNLocMax2[i][j]) ;
      
      fhCentralityEtaNLocMaxN[i][j]  = new TH2F(Form("hCentralityEtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("E vs Centrality, selected pi0 cluster with NLM>1, %s",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,100,0,100);
      fhCentralityEtaNLocMaxN[i][j]->SetYTitle("Centrality");
      fhCentralityEtaNLocMaxN[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityEtaNLocMaxN[i][j]) ;
      
      
      fhM02Pi0NLocMax1[i][j]     = new TH2F(Form("hM02Pi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max = 1",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhM02Pi0NLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0NLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0NLocMax1[i][j]) ;
      
      fhM02EtaNLocMax1[i][j]     = new TH2F(Form("hM02EtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 1",
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhM02EtaNLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02EtaNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02EtaNLocMax1[i][j]) ;
      
      fhM02ConNLocMax1[i][j]    = new TH2F(Form("hM02ConNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 1",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02ConNLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02ConNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02ConNLocMax1[i][j]) ; 
      
      fhM02Pi0NLocMax2[i][j]     = new TH2F(Form("hM02Pi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max = 2",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02Pi0NLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0NLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0NLocMax2[i][j]) ; 
      
      fhM02EtaNLocMax2[i][j]     = new TH2F(Form("hM02EtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 2",
                                               GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02EtaNLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02EtaNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02EtaNLocMax2[i][j]) ;
      
      fhM02ConNLocMax2[i][j]    = new TH2F(Form("hM02ConNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 2",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02ConNLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02ConNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02ConNLocMax2[i][j]) ; 
      
      fhM02Pi0NLocMaxN[i][j]     = new TH2F(Form("hM02Pi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max > 2",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02Pi0NLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0NLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0NLocMaxN[i][j]) ; 
      
      fhM02EtaNLocMaxN[i][j]     = new TH2F(Form("hM02EtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max > 2",
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02EtaNLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02EtaNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02EtaNLocMaxN[i][j]) ; 
      
      fhM02ConNLocMaxN[i][j]    = new TH2F(Form("hM02ConNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f], %s, for N Local max > 2",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02ConNLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02ConNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02ConNLocMaxN[i][j]) ;
            
      
      fhMassPi0NLocMax1[i][j]     = new TH2F(Form("hMassPi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max = 1",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassPi0NLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassPi0NLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassPi0NLocMax1[i][j]) ; 

      
      fhMassEtaNLocMax1[i][j]     = new TH2F(Form("hMassEtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 1",
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassEtaNLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassEtaNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassEtaNLocMax1[i][j]) ; 
      
      fhMassConNLocMax1[i][j]    = new TH2F(Form("hMassConNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 1",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassConNLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassConNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassConNLocMax1[i][j]) ; 
      
      fhMassPi0NLocMax2[i][j]     = new TH2F(Form("hMassPi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max = 2",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassPi0NLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassPi0NLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassPi0NLocMax2[i][j]) ; 
      
      
      fhMassEtaNLocMax2[i][j]     = new TH2F(Form("hMassEtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 2",
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassEtaNLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassEtaNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassEtaNLocMax2[i][j]) ; 
      
      fhMassConNLocMax2[i][j]    = new TH2F(Form("hMassConNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 2",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassConNLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassConNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassConNLocMax2[i][j]) ; 
      
      fhMassPi0NLocMaxN[i][j]     = new TH2F(Form("hMassPi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max > 2",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassPi0NLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassPi0NLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassPi0NLocMaxN[i][j]) ; 
      
      fhMassEtaNLocMaxN[i][j]     = new TH2F(Form("hMassEtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max > 2", 
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassEtaNLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassEtaNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassEtaNLocMaxN[i][j]) ; 
      
      fhMassConNLocMaxN[i][j]    = new TH2F(Form("hMassConNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Mass vs E for mass range [%2.2f-%2.2f], %s, for N Local max > 2",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassConNLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassConNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassConNLocMaxN[i][j]) ; 
      
      
      fhAsyPi0NLocMax1[i][j]     = new TH2F(Form("hAsyPi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max = 1",
                                                 GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyPi0NLocMax1[i][j]   ->SetYTitle("Asymmetry");
      fhAsyPi0NLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyPi0NLocMax1[i][j]) ; 
      
      fhAsyEtaNLocMax1[i][j]     = new TH2F(Form("hAsyEtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 1",
                                                 GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyEtaNLocMax1[i][j]   ->SetYTitle("Asymmetry");
      fhAsyEtaNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyEtaNLocMax1[i][j]) ; 
      
      fhAsyConNLocMax1[i][j]    = new TH2F(Form("hAsyConNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 1",
                                                GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyConNLocMax1[i][j]   ->SetYTitle("Asymmetry");
      fhAsyConNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyConNLocMax1[i][j]) ; 
      
      fhAsyPi0NLocMax2[i][j]     = new TH2F(Form("hAsyPi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max = 2",
                                                 GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyPi0NLocMax2[i][j]   ->SetYTitle("Asymmetry");
      fhAsyPi0NLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyPi0NLocMax2[i][j]) ; 
      
      fhAsyEtaNLocMax2[i][j]     = new TH2F(Form("hAsyEtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 2",
                                                 GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyEtaNLocMax2[i][j]   ->SetYTitle("Asymmetry");
      fhAsyEtaNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyEtaNLocMax2[i][j]) ; 
      
      fhAsyConNLocMax2[i][j]    = new TH2F(Form("hAsyConNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 2",
                                                GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyConNLocMax2[i][j]   ->SetYTitle("Asymmetry");
      fhAsyConNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyConNLocMax2[i][j]) ; 
      
      fhAsyPi0NLocMaxN[i][j]     = new TH2F(Form("hAsyPi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max > 2",
                                                 GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyPi0NLocMaxN[i][j]   ->SetYTitle("Asymmetry");
      fhAsyPi0NLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyPi0NLocMaxN[i][j]) ; 
      
      fhAsyEtaNLocMaxN[i][j]     = new TH2F(Form("hAsyEtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max > 2",
                                                 GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyEtaNLocMaxN[i][j]   ->SetYTitle("Asymmetry");
      fhAsyEtaNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyEtaNLocMaxN[i][j]) ; 
      
      fhAsyConNLocMaxN[i][j]    = new TH2F(Form("hAsyConNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Asymmetry vs E for mass range [%2.2f-%2.2f], %s, for N Local max > 2",
                                                GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyConNLocMaxN[i][j]   ->SetYTitle("Asymmetry");
      fhAsyConNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyConNLocMaxN[i][j]) ; 
      
    } // matched, not matched
    
    if(fFillEbinHisto)
    {
      for(Int_t j = 0; j < 4; j++)
      {
        
        fhMassSplitEFractionNLocMax1Ebin[i][j]  = new TH2F(Form("hMassSplitEFractionNLocMax1%sEbin%d",pname[i].Data(),j),
                                                           Form("Invariant mass of 2 highest energy cells vs (E1+E2)/Ecluster, %s, E bin %d",ptype[i].Data(),j),
                                                           120,0,1.2,mbins,mmin,mmax);
        fhMassSplitEFractionNLocMax1Ebin[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassSplitEFractionNLocMax1Ebin[i][j]->SetXTitle("(E_{split1}+E_{split2})/E_{cluster}");
        outputContainer->Add(fhMassSplitEFractionNLocMax1Ebin[i][j]) ;
        
        fhMassSplitEFractionNLocMax2Ebin[i][j]  = new TH2F(Form("hMassSplitEFractionNLocMax2%sEbin%d",pname[i].Data(),j),
                                                           Form("Invariant mass of 2 local maxima cells vs (E1+E2)/Ecluster, %s, E bin %d",ptype[i].Data(),j),
                                                           120,0,1.2,mbins,mmin,mmax);
        fhMassSplitEFractionNLocMax2Ebin[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassSplitEFractionNLocMax2Ebin[i][j]->SetXTitle("(E_{split1}+E_{split2})/E_{cluster}");
        outputContainer->Add(fhMassSplitEFractionNLocMax2Ebin[i][j]) ;
        
        fhMassSplitEFractionNLocMaxNEbin[i][j]  = new TH2F(Form("hMassSplitEFractionNLocMaxN%sEbin%d",pname[i].Data(),j),
                                                           Form("Invariant mass of N>2 local maxima cells vs (E1+E2)/Ecluster, %s, E bin %d",ptype[i].Data(),j),
                                                           120,0,1.2,mbins,mmin,mmax);
        fhMassSplitEFractionNLocMaxNEbin[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassSplitEFractionNLocMaxNEbin[i][j]->SetXTitle("(E_{split1}+E_{split2})/E_{cluster}");
        outputContainer->Add(fhMassSplitEFractionNLocMaxNEbin[i][j]) ;
        
        if(i>0 && fFillMCHisto) // skip first entry in array, general case not filled
        {
          fhMCGenFracNLocMaxEbin[i][j]  = new TH2F(Form("hMCGenFracNLocMax%sEbin%d",pname[i].Data(),j),
                                                   Form("NLM vs E, %s, E bin %d",ptype[i].Data(),j),
                                                   200,0,2,nMaxBins,0,nMaxBins);
          fhMCGenFracNLocMaxEbin[i][j]->SetYTitle("NLM");
          fhMCGenFracNLocMaxEbin[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMCGenFracNLocMaxEbin[i][j]) ;
          
          fhMCGenFracNLocMaxEbinMatched[i][j]  = new TH2F(Form("hMCGenFracNLocMax%sEbin%dMatched",pname[i].Data(),j),
                                                          Form("NLM vs E, %s, E bin %d, matched to a track",ptype[i].Data(),j),
                                                          200,0,2,nMaxBins,0,nMaxBins);
          fhMCGenFracNLocMaxEbinMatched[i][j]->SetYTitle("NLM");
          fhMCGenFracNLocMaxEbinMatched[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMCGenFracNLocMaxEbinMatched[i][j]) ;
          
          fhMassMCGenFracNLocMax1Ebin[i][j]  = new TH2F(Form("hMassMCGenFracNLocMax1%sEbin%d",pname[i].Data(),j),
                                                        Form("Invariant mass of 2 highest energy cells vs E, %s, E bin %d",ptype[i].Data(),j),
                                                        200,0,2,mbins,mmin,mmax);
          fhMassMCGenFracNLocMax1Ebin[i][j]->SetYTitle("M (GeV/c^{2})");
          fhMassMCGenFracNLocMax1Ebin[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMassMCGenFracNLocMax1Ebin[i][j]) ;
          
          fhMassMCGenFracNLocMax2Ebin[i][j]  = new TH2F(Form("hMassMCGenFracNLocMax2%sEbin%d",pname[i].Data(),j),
                                                        Form("Invariant mass of 2 local maxima cells vs E, %s, E bin %d",ptype[i].Data(),j),
                                                        200,0,2,mbins,mmin,mmax);
          fhMassMCGenFracNLocMax2Ebin[i][j]->SetYTitle("M (GeV/c^{2})");
          fhMassMCGenFracNLocMax2Ebin[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMassMCGenFracNLocMax2Ebin[i][j]) ;
          
          fhMassMCGenFracNLocMaxNEbin[i][j]  = new TH2F(Form("hMassMCGenFracNLocMaxN%sEbin%d",pname[i].Data(),j),
                                                        Form("Invariant mass of N>2 local maxima cells vs E, %s, E bin %d",ptype[i].Data(),j),
                                                        200,0,2,mbins,mmin,mmax);
          fhMassMCGenFracNLocMaxNEbin[i][j]->SetYTitle("M (GeV/c^{2})");
          fhMassMCGenFracNLocMaxNEbin[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMassMCGenFracNLocMaxNEbin[i][j]) ;
          
          fhM02MCGenFracNLocMax1Ebin[i][j]     = new TH2F(Form("hM02MCGenFracNLocMax1%sEbin%d",pname[i].Data(),j),
                                                          Form("#lambda_{0}^{2} vs E for N max  = 1 %s, E bin %d",ptype[i].Data(), j),
                                                          200,0,2,ssbins,ssmin,ssmax);
          fhM02MCGenFracNLocMax1Ebin[i][j]   ->SetYTitle("#lambda_{0}^{2}");
          fhM02MCGenFracNLocMax1Ebin[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhM02MCGenFracNLocMax1Ebin[i][j]) ;
          
          fhM02MCGenFracNLocMax2Ebin[i][j]     = new TH2F(Form("hM02MCGenFracNLocMax2%sEbin%d",pname[i].Data(),j),
                                                          Form("#lambda_{0}^{2} vs E for N max  = 2 %s, E bin %d",ptype[i].Data(),j),
                                                          200,0,2,ssbins,ssmin,ssmax);
          fhM02MCGenFracNLocMax2Ebin[i][j]   ->SetYTitle("#lambda_{0}^{2}");
          fhM02MCGenFracNLocMax2Ebin[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhM02MCGenFracNLocMax2Ebin[i][j]) ;
          
          fhM02MCGenFracNLocMaxNEbin[i][j]    = new TH2F(Form("hM02MCGenFracNLocMaxN%sEbin%d",pname[i].Data(),j),
                                                         Form("#lambda_{0}^{2} vs E for N max  > 2 %s, E bin %d",ptype[i].Data(),j),
                                                         200,0,2,ssbins,ssmin,ssmax);
          fhM02MCGenFracNLocMaxNEbin[i][j]   ->SetYTitle("#lambda_{0}^{2}");
          fhM02MCGenFracNLocMaxNEbin[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhM02MCGenFracNLocMaxNEbin[i][j]) ;
        }
      }
    }
  } // MC particle list
 
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
  
  if(fFillEbinHisto)
  {
    for(Int_t i = 0; i < 4; i++)
    {
      fhMassM02NLocMax1Ebin[i]  = new TH2F(Form("hMassM02NLocMax1Ebin%d",i),
                                           Form("Invariant mass of split clusters vs #lambda_{0}^{2}, NLM=1, E bin %d",i),
                                           ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMax1Ebin[i]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMax1Ebin[i]) ;   
      
      fhMassM02NLocMax2Ebin[i]  = new TH2F(Form("hMassM02NLocMax2Ebin%d",i),
                                           Form("Invariant mass of split clusters vs #lambda_{0}^{2}, NLM=2, E bin %d",i),
                                           ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMax2Ebin[i]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMax2Ebin[i]) ;   
      
      fhMassM02NLocMaxNEbin[i]  = new TH2F(Form("hMassM02NLocMaxNEbin%d",i),
                                           Form("Invariant mass of split clusters vs vs #lambda_{0}^{2}, NLM>2, E bin %d",i),
                                           ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMaxNEbin[i]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMaxNEbin[i]) ; 
      
      
      fhMassAsyNLocMax1Ebin[i]  = new TH2F(Form("hMassAsyNLocMax1Ebin%d",i),
                                           Form("Invariant mass of split clusters vs split asymmetry, NLM=1, E bin %d",i),
                                           200,-1,1,mbins,mmin,mmax);
      fhMassAsyNLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassAsyNLocMax1Ebin[i]->SetXTitle("asymmetry");
      outputContainer->Add(fhMassAsyNLocMax1Ebin[i]) ;
      
      fhMassAsyNLocMax2Ebin[i]  = new TH2F(Form("hMassAsyNLocMax2Ebin%d",i),
                                           Form("Invariant mass of split clusters vs split asymmetry, NLM=2, E bin %d",i),
                                           200,-1,1,mbins,mmin,mmax);
      fhMassAsyNLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassAsyNLocMax2Ebin[i]->SetXTitle("asymmetry");
      outputContainer->Add(fhMassAsyNLocMax2Ebin[i]) ;
      
      fhMassAsyNLocMaxNEbin[i]  = new TH2F(Form("hMassAsyNLocMaxNEbin%d",i),
                                           Form("Invariant mass of split clusters vs split asymmetry, NLM>2, E bin %d",i),
                                           200,-1,1,mbins,mmin,mmax);
      fhMassAsyNLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassAsyNLocMaxNEbin[i]->SetXTitle("asymmetry");
      outputContainer->Add(fhMassAsyNLocMaxNEbin[i]) ;
      
      
      if(IsDataMC() && fFillMCHisto)
      {
        fhMCAsymM02NLocMax1MCPi0Ebin[i]  = new TH2F(Form("hMCAsymM02NLocMax1MCPi0Ebin%d",i),
                                                    Form("Asymmetry of MC #pi^{0} vs #lambda_{0}^{2}, NLM=1, E bin %d",i),
                                                    ssbins,ssmin,ssmax,100,0,1);
        fhMCAsymM02NLocMax1MCPi0Ebin[i]->SetYTitle("Decay asymmetry");
        fhMCAsymM02NLocMax1MCPi0Ebin[i]->SetXTitle("#lambda_{0}^{2}");
        outputContainer->Add(fhMCAsymM02NLocMax1MCPi0Ebin[i]) ;
        
        fhMCAsymM02NLocMax2MCPi0Ebin[i]  = new TH2F(Form("hMCAsymM02NLocMax2MCPi0Ebin%d",i),
                                                    Form("Asymmetry of MC #pi^{0} vs #lambda_{0}^{2}, NLM=2, E bin %d",i),
                                                    ssbins,ssmin,ssmax,100,0,1);
        fhMCAsymM02NLocMax2MCPi0Ebin[i]->SetYTitle("Decay asymmetry");
        fhMCAsymM02NLocMax2MCPi0Ebin[i]->SetXTitle("#lambda_{0}^{2}");
        outputContainer->Add(fhMCAsymM02NLocMax2MCPi0Ebin[i]) ;
        
        fhMCAsymM02NLocMaxNMCPi0Ebin[i]  = new TH2F(Form("hMCAsymM02NLocMaxNMCPi0Ebin%d",i),
                                                    Form("Asymmetry of MC #pi^{0} vs #lambda_{0}^{2}, NLM>2, E bin %d",i),
                                                    ssbins,ssmin,ssmax,100,0,1);
        fhMCAsymM02NLocMaxNMCPi0Ebin[i]->SetYTitle("Decay asymmetry");
        fhMCAsymM02NLocMaxNMCPi0Ebin[i]->SetXTitle("#lambda_{0}^{2}");
        outputContainer->Add(fhMCAsymM02NLocMaxNMCPi0Ebin[i]) ;    
        
        
        fhAsyMCGenRecoNLocMax1EbinPi0[i]  = new TH2F(Form("hAsyMCGenRecoNLocMax1Ebin%dPi0",i),
                                                     Form("Generated vs reconstructed asymmetry of split clusters from pi0, NLM=1, E bin %d",i),
                                                     200,-1,1,200,-1,1);
        fhAsyMCGenRecoNLocMax1EbinPi0[i]->SetYTitle("M (GeV/c^{2})");
        fhAsyMCGenRecoNLocMax1EbinPi0[i]->SetXTitle("asymmetry");
        outputContainer->Add(fhAsyMCGenRecoNLocMax1EbinPi0[i]) ;
        
        fhAsyMCGenRecoNLocMax2EbinPi0[i]  = new TH2F(Form("hAsyMCGenRecoNLocMax2Ebin%dPi0",i),
                                                     Form("Generated vs reconstructed asymmetry of split clusters from pi0, NLM=2, E bin %d",i),
                                                     200,-1,1,200,-1,1);
        fhAsyMCGenRecoNLocMax2EbinPi0[i]->SetYTitle("M (GeV/c^{2})");
        fhAsyMCGenRecoNLocMax2EbinPi0[i]->SetXTitle("asymmetry");
        outputContainer->Add(fhAsyMCGenRecoNLocMax2EbinPi0[i]) ;
        
        fhAsyMCGenRecoNLocMaxNEbinPi0[i]  = new TH2F(Form("hAsyMCGenRecoNLocMaxNEbin%dPi0",i),
                                                     Form("Generated vs reconstructed asymmetry of split clusters from pi0, NLM>2, E bin %d",i),
                                                     200,-1,1,200,-1,1);
        fhAsyMCGenRecoNLocMaxNEbinPi0[i]->SetYTitle("M (GeV/c^{2})");
        fhAsyMCGenRecoNLocMaxNEbinPi0[i]->SetXTitle("asymmetry");
        outputContainer->Add(fhAsyMCGenRecoNLocMaxNEbinPi0[i]) ;
      }
      
      if(fFillSSExtraHisto)
      {
        fhMassDispEtaNLocMax1Ebin[i]  = new TH2F(Form("hMassDispEtaNLocMax1Ebin%d",i),
                                                 Form("Invariant mass of 2 highest energy cells #sigma_{#eta #eta}^{2}, E bin %d",i),
                                                 ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispEtaNLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispEtaNLocMax1Ebin[i]->SetXTitle("#sigma_{#eta #eta}^{2}");
        outputContainer->Add(fhMassDispEtaNLocMax1Ebin[i]) ;   
        
        fhMassDispEtaNLocMax2Ebin[i]  = new TH2F(Form("hMassDispEtaNLocMax2Ebin%d",i),
                                                 Form("Invariant mass of 2 local maxima cells #sigma_{#eta #eta}^{2}, E bin %d",i),
                                                 ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispEtaNLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispEtaNLocMax2Ebin[i]->SetXTitle("#sigma_{#eta #eta}^{2}");
        outputContainer->Add(fhMassDispEtaNLocMax2Ebin[i]) ;   
        
        fhMassDispEtaNLocMaxNEbin[i]  = new TH2F(Form("hMassDispEtaNLocMaxNEbin%d",i),
                                                 Form("Invariant mass of N>2 local maxima cells vs #sigma_{#eta #eta}^{2}, E bin %d",i),
                                                 ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispEtaNLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispEtaNLocMaxNEbin[i]->SetXTitle("#sigma_{#eta #eta}^{2}");
        outputContainer->Add(fhMassDispEtaNLocMaxNEbin[i]) ;   
        
        fhMassDispPhiNLocMax1Ebin[i]  = new TH2F(Form("hMassDispPhiNLocMax1Ebin%d",i),
                                                 Form("Invariant mass of 2 highest energy cells #sigma_{#phi #phi}^{2}, E bin %d",i),
                                                 ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispPhiNLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispPhiNLocMax1Ebin[i]->SetXTitle("#sigma_{#phi #phi}^{2}");
        outputContainer->Add(fhMassDispPhiNLocMax1Ebin[i]) ;   
        
        fhMassDispPhiNLocMax2Ebin[i]  = new TH2F(Form("hMassDispPhiNLocMax2Ebin%d",i),
                                                 Form("Invariant mass of 2 local maxima cells #sigma_{#phi #phi}^{2}, E bin %d",i),
                                                 ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispPhiNLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispPhiNLocMax2Ebin[i]->SetXTitle("#sigma_{#phi #phi}^{2}");
        outputContainer->Add(fhMassDispPhiNLocMax2Ebin[i]) ;   
        
        fhMassDispPhiNLocMaxNEbin[i]  = new TH2F(Form("hMassDispPhiNLocMaxNEbin%d",i),
                                                 Form("Invariant mass of N>2 local maxima cells vs #sigma_{#phi #phi}^{2}, E bin %d",i),
                                                 ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispPhiNLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispPhiNLocMaxNEbin[i]->SetXTitle("#sigma_{#phi #phi}^{2}");
        outputContainer->Add(fhMassDispPhiNLocMaxNEbin[i]) ;   
        
        fhMassDispAsyNLocMax1Ebin[i]  = new TH2F(Form("hMassDispAsyNLocMax1Ebin%d",i),
                                                 Form("Invariant mass of 2 highest energy cells A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), E bin %d",i),
                                                 200,-1,1,mbins,mmin,mmax); 
        fhMassDispAsyNLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispAsyNLocMax1Ebin[i]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
        outputContainer->Add(fhMassDispAsyNLocMax1Ebin[i]) ;   
        
        fhMassDispAsyNLocMax2Ebin[i]  = new TH2F(Form("hMassDispAsyNLocMax2Ebin%d",i),
                                                 Form("Invariant mass of 2 local maxima cells A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), E bin %d",i),
                                                 200,-1,1,mbins,mmin,mmax); 
        fhMassDispAsyNLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispAsyNLocMax2Ebin[i]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
        outputContainer->Add(fhMassDispAsyNLocMax2Ebin[i]) ;   
        
        fhMassDispAsyNLocMaxNEbin[i]  = new TH2F(Form("hMassDispAsyNLocMaxNEbin%d",i),
                                                 Form("Invariant mass of N>2 local maxima cells vs A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), E bin %d",i),
                                                 200,-1,1,mbins,mmin,mmax); 
        fhMassDispAsyNLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
        fhMassDispAsyNLocMaxNEbin[i]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
        outputContainer->Add(fhMassDispAsyNLocMaxNEbin[i]) ;   
      }
    }
  }
    
  if(IsDataMC() && fFillMCHisto)
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
    for(Int_t j = 0; j < nMatched; j++)
    {  
      
      fhAnglePairNLocMax1[j]  = new TH2F(Form("hAnglePairNLocMax1%s",sMatched[j].Data()),
                                        Form("Opening angle of 2 highest energy cells vs pair Energy, %s",sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,200,0,0.2); 
      fhAnglePairNLocMax1[j]->SetYTitle("#alpha (rad)");
      fhAnglePairNLocMax1[j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairNLocMax1[j]) ;   
      
      fhAnglePairNLocMax2[j]  = new TH2F(Form("hAnglePairNLocMax2%s",sMatched[j].Data()),
                                        Form("Opening angle of 2 local maxima cells vs Energy, %s",sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,200,0,0.2); 
      fhAnglePairNLocMax2[j]->SetYTitle("#alpha (rad)");
      fhAnglePairNLocMax2[j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairNLocMax2[j]) ;   
      
      fhAnglePairNLocMaxN[j]  = new TH2F(Form("hAnglePairNLocMaxN%s",sMatched[j].Data()),
                                        Form("Opening angle of N>2 local maxima cells vs Energy, %s",sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,200,0,0.2); 
      fhAnglePairNLocMaxN[j]->SetYTitle("#alpha (rad)");
      fhAnglePairNLocMaxN[j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairNLocMaxN[j]) ;   
      
      fhAnglePairMassNLocMax1[j]  = new TH2F(Form("hAnglePairMassNLocMax1%s",sMatched[j].Data()),
                                            Form("Opening angle of 2 highest energy cells vs Mass for E > 12 GeV, %s",sMatched[j].Data()),
                                            mbins,mmin,mmax,200,0,0.2); 
      fhAnglePairMassNLocMax1[j]->SetXTitle("M (GeV/c^{2})");
      fhAnglePairMassNLocMax1[j]->SetYTitle("#alpha (rad)");
      outputContainer->Add(fhAnglePairMassNLocMax1[j]) ;   
      
      fhAnglePairMassNLocMax2[j]  = new TH2F(Form("hAnglePairMassNLocMax2%s",sMatched[j].Data()),
                                            Form("Opening angle of 2 local maxima cells vs Mass for E > 12 GeV, %s",sMatched[j].Data()),
                                            mbins,mmin,mmax,200,0,0.2); 
      fhAnglePairMassNLocMax2[j]->SetXTitle("M (GeV/c^{2})");
      fhAnglePairMassNLocMax2[j]->SetYTitle("#alpha (rad)");
      outputContainer->Add(fhAnglePairMassNLocMax2[j]) ;   
      
      fhAnglePairMassNLocMaxN[j]  = new TH2F(Form("hAnglePairMassNLocMaxN%s",sMatched[j].Data()),
                                            Form("Opening angle of N>2 local maxima cells vs Mass for E > 12 GeV, %s",sMatched[j].Data()),
                                            mbins,mmin,mmax,200,0,0.2); 
      fhAnglePairMassNLocMaxN[j]->SetXTitle("M (GeV/c^{2})");
      fhAnglePairMassNLocMaxN[j]->SetYTitle("#alpha (rad)");
      outputContainer->Add(fhAnglePairMassNLocMaxN[j]) ;  
      
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
  
  TString snlm[] = {"1","2","N"};

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
                                              Form("#lambda_{0}^{2} vs E, with W0 = %2.2f, for N Local max = %s", fSSWeight[i], snlm[nlm].Data()),
                                              nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhM02WeightPi0[nlm][i]   ->SetYTitle("#lambda_{0}^{2}");
        fhM02WeightPi0[nlm][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhM02WeightPi0[nlm][i]) ;
      }
      
      for(Int_t i = 0; i < fSSECellCutN; i++)
      {
        fhM02ECellCutPi0[nlm][i]     = new TH2F(Form("hM02Pi0NLocMax%s_Ecell%d",snlm[nlm].Data(),i),
                                              Form("#lambda_{0}^{2} vs E, with Ecell > %2.2f, for N Local max = %s", fSSECellCut[i], snlm[nlm].Data()),
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
        
        
        if(i < 5)
        {
          fhMCPi0MassM02Overlap0[j][i-1]  = new TH2F(Form("hMCPi0MassM02Overlap0NLocMax%sEbin%d",snlm[j].Data(),i-1),
                                                   Form("Overlap 0, Mass vs #lambda_{0}^{2}, NLM=%s, E bin %d",snlm[j].Data(),i-1),
                                                   ssbins,ssmin,ssmax,mbins,mmin,mmax);
          fhMCPi0MassM02Overlap0[j][i-1]->SetYTitle("M (GeV/c^{2})");
          fhMCPi0MassM02Overlap0[j][i-1]->SetXTitle("#lambda_{0}^{2}");
          outputContainer->Add(fhMCPi0MassM02Overlap0[j][i-1]) ;
          
          fhMCPi0MassM02Overlap1[j][i-1]  = new TH2F(Form("hMCPi0MassM02Overlap1NLocMax%sEbin%d",snlm[j].Data(),i-1),
                                                   Form("Overlap 1, Mass vs #lambda_{0}^{2}, NLM=%s, E bin %d",snlm[j].Data(),i-1),
                                                   ssbins,ssmin,ssmax,mbins,mmin,mmax);
          fhMCPi0MassM02Overlap1[j][i-1]->SetYTitle("M (GeV/c^{2})");
          fhMCPi0MassM02Overlap1[j][i-1]->SetXTitle("#lambda_{0}^{2}");
          outputContainer->Add(fhMCPi0MassM02Overlap1[j][i-1]) ;
          
          fhMCPi0MassM02OverlapN[j][i-1]  = new TH2F(Form("hMCPi0MassM02OverlapNNLocMax%sEbin%d",snlm[j].Data(),i-1),
                                                   Form("Overlap N, Mass vs #lambda_{0}^{2}, NLM=%s, E bin %d",snlm[j].Data(),i-1),
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
          
          if(i < 5)
          {
            fhMCPi0MassM02Overlap0Match[j][i-1]  = new TH2F(Form("hMCPi0MassM02Overlap0NLocMax%sEbin%dMatched",snlm[j].Data(),i-1),
                                                     Form("Overlap 0, Mass vs #lambda_{0}^{2}, NLM=%s, E bin %d, Track Matched",snlm[j].Data(),i-1),
                                                     ssbins,ssmin,ssmax,mbins,mmin,mmax);
            fhMCPi0MassM02Overlap0Match[j][i-1]->SetYTitle("M (GeV/c^{2})");
            fhMCPi0MassM02Overlap0Match[j][i-1]->SetXTitle("#lambda_{0}^{2}");
            outputContainer->Add(fhMCPi0MassM02Overlap0Match[j][i-1]) ;
            
            fhMCPi0MassM02Overlap1Match[j][i-1]  = new TH2F(Form("hMCPi0MassM02Overlap1NLocMax%sEbin%dMatched",snlm[j].Data(),i-1),
                                                     Form("Overlap 1, Mass vs #lambda_{0}^{2}, NLM=%s, E bin %d, Track Matched",snlm[j].Data(),i-1),
                                                     ssbins,ssmin,ssmax,mbins,mmin,mmax);
            fhMCPi0MassM02Overlap1Match[j][i-1]->SetYTitle("M (GeV/c^{2})");
            fhMCPi0MassM02Overlap1Match[j][i-1]->SetXTitle("#lambda_{0}^{2}");
            outputContainer->Add(fhMCPi0MassM02Overlap1Match[j][i-1]) ;
            
            fhMCPi0MassM02OverlapNMatch[j][i-1]  = new TH2F(Form("hMCPi0MassM02OverlapNNLocMax%sEbin%dMatched",snlm[j].Data(),i-1),
                                                     Form("Overlap N, Mass vs #lambda_{0}^{2}, NLM=%s, E bin %d, Track Matched",snlm[j].Data(),i-1),
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

    
    fhMCEOverlapType = new TH2F("hMCEOverlapType","Kind of overlap particle, neutral clusters",
                                nptbins,ptmin,ptmax,5,0,5);
    //fhMCEOverlapType   ->SetYTitle("Overlap Type");
    fhMCEOverlapType->GetYaxis()->SetBinLabel(1 ,"#gamma");
    fhMCEOverlapType->GetYaxis()->SetBinLabel(2 ,"e^{#pm}");
    fhMCEOverlapType->GetYaxis()->SetBinLabel(3 ,"hadron^{#pm}");
    fhMCEOverlapType->GetYaxis()->SetBinLabel(4 ,"hadron^{0}");
    fhMCEOverlapType->GetYaxis()->SetBinLabel(5 ,"??");
    fhMCEOverlapType   ->SetXTitle("Cluster E (GeV)");
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
  
  return outputContainer ;
  
}

//_____________________________________________________________________________
void AliAnaInsideClusterInvariantMass::GetMCIndex(AliVCluster* cluster, Int_t & mcindex)
{
  
  // Assign mc index depending on MC bit set, to be used in histograms arrays
  
  if(!IsDataMC()) return;
  
  Int_t tag	= GetMCAnalysisUtils()->CheckOrigin(cluster->GetLabels(),cluster->GetNLabels(), GetReader());
  
  if      ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) &&
           !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = kmcPi0;
  else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)  )      mcindex = kmcPi0Conv;
  else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)  )      mcindex = kmcEta;
  else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) &&
           !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = kmcPhoton;
  else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) &&
            GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = kmcConversion;
  else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron))   mcindex = kmcElectron;
  else                                                                                mcindex = kmcHadron;
  
}


//______________________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::GetMCPrimaryKine(AliVCluster* cluster, const Int_t mcindex, const Bool_t matched,
                                                        Float_t & eprim, Float_t & asymGen, Int_t & noverlaps )
{
  // Check origin of the candidates, get primary kinematics if overlapped meson decay
  
  if(!IsDataMC()) return ;
  
  Bool_t ok      = kFALSE;
  Int_t  mcLabel = cluster->GetLabel();
  
  TLorentzVector primary = GetMCAnalysisUtils()->GetMother(mcLabel,GetReader(),ok);
  eprim = primary.E();
  
  if(mcindex == kmcPi0 || mcindex == kmcEta || mcindex == kmcPi0Conv)
  {
    if(mcindex == kmcPi0 || mcindex == kmcPi0Conv)
    {
      asymGen = TMath::Abs(GetMCAnalysisUtils()->GetMCDecayAsymmetryForPDG(mcLabel,111,GetReader(),ok));
      TLorentzVector grandmom = GetMCAnalysisUtils()->GetMotherWithPDG(mcLabel,111,GetReader(),ok);
      if(grandmom.E() > 0 && ok) eprim =  grandmom.E();
    }
    else
    {
      asymGen = TMath::Abs(GetMCAnalysisUtils()->GetMCDecayAsymmetryForPDG(mcLabel,221,GetReader(),ok));
      TLorentzVector grandmom = GetMCAnalysisUtils()->GetMotherWithPDG(mcLabel,221,GetReader(),ok);
      if(grandmom.E() > 0 && ok) eprim =  grandmom.E();
    }
  }
  
  if(!fFillMCOverlapHisto) return;
  
  // Compare the primary depositing more energy with the rest,
  // if no photon/electron (conversion) or neutral meson as comon ancestor, consider it as other particle contributing
  Int_t ancPDG = 0, ancStatus = -1;
  TLorentzVector momentum; TVector3 prodVertex;
  Int_t ancLabel = 0;
  noverlaps = 0;
  for (UInt_t ilab = 1; ilab < cluster->GetNLabels(); ilab++ )
  {
    ancLabel = GetMCAnalysisUtils()->CheckCommonAncestor(cluster->GetLabels()[0],cluster->GetLabels()[ilab],
                                                         GetReader(),ancPDG,ancStatus,momentum,prodVertex);
    
    if(ancPDG!=22 && TMath::Abs(ancPDG)!=11 && ancPDG!=111 && ancPDG!=221)
    {
      noverlaps++;
      
      // What is the origin of the overlap?
      Bool_t  mOK = 0,      gOK = 0;
      Int_t   mpdg = -999999,  gpdg = -1;
      Int_t   mstatus = -1, gstatus = -1;
      Int_t gLabel = -1, ggLabel = -1;
      TLorentzVector mother      = GetMCAnalysisUtils()->GetMother     (cluster->GetLabels()[ilab],GetReader(),mpdg,mstatus,mOK);
      TLorentzVector grandmother = GetMCAnalysisUtils()->GetGrandMother(cluster->GetLabels()[ilab],GetReader(),gpdg,gstatus,gOK, gLabel,ggLabel);
      
      //printf("Overlap!, mother pdg %d; grand mother pdg %d",mpdg,gpdg);
      
      if( ( mpdg == 22 || TMath::Abs(mpdg==11) ) &&
          ( gpdg == 22 || TMath::Abs(gpdg==11) ) &&
            gLabel >=0 )
      {
        Int_t label = gLabel;
        while( ( gpdg == 22 || TMath::Abs(gpdg==11) ) && gLabel >=0 )
        {
          mpdg=gpdg;
          grandmother = GetMCAnalysisUtils()->GetGrandMother(label,GetReader(),gpdg,gstatus,ok, gLabel,ggLabel);
          label=gLabel;
        }
      }
         
      //printf("; Final PDG %d\n",mpdg);
      Float_t histobin = -1;
      if     (mpdg==22)      histobin = 0.5;
      else if(TMath::Abs(mpdg)==11) histobin = 1.5;
      else if(mpdg==-999999) histobin = 4.5;
      else {
        Double_t charge = TDatabasePDG::Instance()->GetParticle(mpdg)->Charge();
        if(TMath::Abs(charge) > 0 ) histobin = 2.5;
        else                        histobin = 3.5;
        //printf("charge %f\n",charge);
      }
      
      //printf("pdg = %d, histobin %2.1f\n",mpdg,histobin);
      if(histobin > 0)
      {
        if(matched)fhMCEOverlapType     ->Fill(cluster->E(),histobin);
        else       fhMCEOverlapTypeMatch->Fill(cluster->E(),histobin);
      }
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
    printf("AliAnaInsideClusterInvariantMass::Init() - !!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD())
  {
    printf("AliAnaInsideClusterInvariantMass::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
  if( GetReader()->GetDataType() == AliCaloTrackReader::kMC )
  {
    printf("AliAnaInsideClusterInvariantMass::Init() - !!STOP: You want to use pure MC data!!\n");
    abort();
    
  }
  
}

//_____________________________________________________
void AliAnaInsideClusterInvariantMass::InitParameters()
{
  //Initialize the parameters of the analysis.  
  AddToHistogramsName("AnaPi0InsideClusterInvariantMass_");
  
  fCalorimeter = "EMCAL" ;

  fM02MinCut   = 0.26 ;
  fM02MaxCut   = 10 ;
  
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
    TLorentzVector    l1, l2;
    Int_t    absId1 = -1; Int_t absId2 = -1;
    
    Int_t pidTag = GetCaloPID()->GetIdentifiedParticleTypeFromClusterSplitting(cluster,cells,GetCaloUtils(),
                                                                               GetVertex(0), nMax, mass, angle,
                                                                               l1,l2,absId1,absId2);
    if (nMax <= 0) 
    {
      if(GetDebug() > 0 )
        printf("AliAnaInsideClusterInvariantMass::MakeAnalysisFillHistograms() - No local maximum found! It did not pass CaloPID selection criteria \n");
      
      return;
    }
    
    // Set some index for array histograms
    
    Int_t inlm = -1;
    if     (nMax == 1) inlm = 0;
    else if(nMax == 2) inlm = 1;
    else if(nMax >  2) inlm = 2;
    else printf("Wrong N local maximum -> %d, n cells in cluster %d \n",nMax,nc);

    // Get sub-cluster parameters
    
    Float_t e1 = l1.Energy();
    Float_t e2 = l2.Energy();
    
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
    
    // MC indexes
    
    Int_t mcindex   = -1;
    //Int_t mcLabel   = cluster->GetLabel();
    GetMCIndex(cluster,mcindex);

    // MC primary kine, generation fractions
    
    Float_t eprim   = -1;
    Float_t asymGen = -2;
    Int_t noverlaps =  0;
    GetMCPrimaryKine(cluster,mcindex,matched,eprim,asymGen,noverlaps);
    
    //
    
    // For cluster with MC pi0 and more than 1 maxima
    CheckLocalMaximaMCOrigin(cluster, mcindex);
    
    //----------------
    // Fill histograms
    //----------------
    
    fhNLocMax[0][matched]->Fill(en,nMax);
    if(IsDataMC())
      fhNLocMax[mcindex][matched]->Fill(en,nMax);
    
    if     ( nMax == 1  )
    { 
      fhM02NLocMax1[0][matched]->Fill(en,l0) ; 
      fhSplitEFractionNLocMax1[0][matched]->Fill(en,splitFrac) ;
      
      if(IsDataMC())
      {
        fhM02NLocMax1[mcindex][matched]->Fill(en,l0) ;
        fhSplitEFractionNLocMax1[mcindex][matched]->Fill(en,splitFrac) ;
      }

      if(en > fHistoECut)
      {
        fhMassM02NLocMax1[0][matched]->Fill(l0, mass);
        if( IsDataMC() ) fhMassM02NLocMax1[mcindex][matched]->Fill(l0, mass);

        fhSplitEFractionvsAsyNLocMax1[matched]->Fill(asym,splitFrac) ;
        if(!matched)fhClusterEtaPhiNLocMax1->Fill(eta,phi);
      }
      

    }
    else if( nMax == 2  ) 
    { 
      fhM02NLocMax2[0][matched]->Fill(en,l0) ; 
      fhSplitEFractionNLocMax2[0][matched]->Fill(en,splitFrac) ;
      
      if(IsDataMC())
      {
        fhM02NLocMax2[mcindex][matched]->Fill(en,l0) ;
        fhSplitEFractionNLocMax2[mcindex][matched]->Fill(en,splitFrac) ;
      }
      
      if(en > fHistoECut)
      {
        fhMassM02NLocMax2[0][matched]->Fill(l0,  mass );
        if( IsDataMC() ) fhMassM02NLocMax2[mcindex][matched]->Fill(l0,mass);
        
        fhSplitEFractionvsAsyNLocMax2[matched]->Fill(asym,splitFrac) ;
        if(!matched)fhClusterEtaPhiNLocMax2->Fill(eta,phi);
      }
    }
    else if( nMax >= 3  )
    { 
      fhM02NLocMaxN[0][matched]->Fill(en,l0) ; 
      fhSplitEFractionNLocMaxN[0][matched]->Fill(en,splitFrac) ;
      
      if(IsDataMC())
      {
        fhM02NLocMaxN[mcindex][matched]->Fill(en,l0) ;
        fhSplitEFractionNLocMaxN[mcindex][matched]->Fill(en,splitFrac) ;
      }
      
      if(en > fHistoECut)
      {
        
        fhMassM02NLocMaxN[0][matched]->Fill(l0,mass);
        if( IsDataMC() ) fhMassM02NLocMaxN[mcindex][matched]->Fill(l0,mass);
        
        fhSplitEFractionvsAsyNLocMaxN[matched]->Fill(asym,splitFrac) ;
        if(!matched)fhClusterEtaPhiNLocMaxN->Fill(eta,phi);
      }
    }
    
    //
    
    FillTrackMatchingHistograms(cluster,nMax,mcindex);
    
    //
      
    FillSSExtraHistograms(cluster, nMax, matched,mcindex,mass,ebin)  ;
    
    //
    
    FillMCHistograms(en,e1,e2,ebin,mcindex,l0,mass,
                     nMax,matched,splitFrac, asym, eprim,asymGen);
    
    //
    
    FillMCOverlapHistograms(en,mass,l0,inlm,ebin,matched,mcindex,noverlaps);

    //
    
    if(!matched) FillEBinHistograms(ebin,nMax,mcindex,splitFrac,mass,asym,l0);
    
    //---------------------------------------------------------------------
    // From here only if M02 is large but not too large, fill histograms 
    //---------------------------------------------------------------------
    
    if( l0 < fM02MinCut || l0 > fM02MaxCut ) continue ;
    
    FillAngleHistograms(matched,nMax,en,angle,mass);
    
    Bool_t m02OK = GetCaloPID()->IsInPi0M02Range(en,l0,nMax);
    Bool_t asyOK = GetCaloPID()->IsInPi0SplitAsymmetryRange(en,asym,nMax);
    
    Float_t efrac      = eprim/en;
    Float_t efracSplit = 0;
    if(e1+e2 > 0) efracSplit = eprim/(e1+e2);
    
    Float_t cent = GetEventCentrality();
    Float_t evp  = GetEventPlaneAngle();
    
    fhNLocMaxM02Cut[0][matched]->Fill(en,nMax);
    if(IsDataMC()) fhNLocMaxM02Cut[mcindex][matched]->Fill(en,nMax);
    
    Float_t splitFracMin = GetCaloPID()->GetSplitEnergyFractionMinimum(inlm) ;
    
    if     (nMax==1) 
    {
      fhMassNLocMax1[0][matched]->Fill(en,mass );
      fhAsymNLocMax1[0][matched]->Fill(en,asym );
      
      // Effect of cuts in mass histograms

      if(!matched)
      {
        if(m02OK)
        {
          fhMassM02CutNLocMax1->Fill(en,mass);
          fhAsymM02CutNLocMax1->Fill(en,asym );
          if(splitFrac > splitFracMin) fhMassSplitECutNLocMax1->Fill(en,mass );
        } // m02
      } // split frac
      
      if(m02OK && asyOK)
      {
        fhSplitEFractionAfterCutsNLocMax1[0][matched]->Fill(en,splitFrac);
        if(splitFrac > splitFracMin) fhMassAfterCutsNLocMax1[0][matched]->Fill(en,mass);
        
        if(!matched && IsDataMC() && fFillMCHisto && mcindex==kmcPi0)
        {
          fhMCGenFracAfterCutsNLocMax1MCPi0      ->Fill(en   ,  efrac     );
          fhMCGenSplitEFracAfterCutsNLocMax1MCPi0->Fill(en   ,  efracSplit);
        }
      }
      
            
      if     (pidTag==AliCaloPID::kPhoton)
      { fhM02ConNLocMax1 [0][matched]->Fill(en,l0);
        fhMassConNLocMax1[0][matched]->Fill(en,mass);
        fhAsyConNLocMax1 [0][matched]->Fill(en,asym);
      }
      else if(pidTag==AliCaloPID::kPi0   )
      {
        fhM02Pi0NLocMax1[0][matched]->Fill(en,l0);
        fhMassPi0NLocMax1[0][matched]->Fill(en,mass);
        fhAsyPi0NLocMax1[0][matched]->Fill(en,asym);
        fhCentralityPi0NLocMax1[0][matched]->Fill(en,cent) ;
        if(!matched)
        {
          fhEventPlanePi0NLocMax1->Fill(en,evp) ;
          if(en > fHistoECut)fhPi0EtaPhiNLocMax1->Fill(eta,phi);
          FillSSWeightHistograms(cluster, 0, absId1, absId2);
          fhPi0EPairDiffTimeNLM1->Fill(e1+e2,t12diff);
        }
      }
      else if(pidTag==AliCaloPID::kEta)
      {
        fhM02EtaNLocMax1 [0][matched]->Fill(en,l0);
        fhMassEtaNLocMax1[0][matched]->Fill(en,mass);
        fhAsyEtaNLocMax1 [0][matched]->Fill(en,asym);
        fhCentralityEtaNLocMax1[0][matched]->Fill(en,cent) ;
        if(!matched)
        {
          fhEventPlaneEtaNLocMax1->Fill(en,evp) ;
          if(en > fHistoECut)fhEtaEtaPhiNLocMax1->Fill(eta,phi);
          fhEtaEPairDiffTimeNLM1->Fill(e1+e2,t12diff);
        }
      }
      
    }
    else if(nMax==2) 
    {
      fhMassNLocMax2[0][matched]->Fill(en,mass );
      fhAsymNLocMax2[0][matched]->Fill(en,asym );
      
      // Effect of cuts in mass histograms
      if(!matched)
      {
        if(m02OK)
        {
          fhMassM02CutNLocMax2->Fill(en,mass);
          fhAsymM02CutNLocMax2->Fill(en,asym );
          if(splitFrac > splitFracMin) fhMassSplitECutNLocMax2->Fill(en,mass );
        } // m02
      } // split frac
      
      if(m02OK && asyOK)
      {
        fhSplitEFractionAfterCutsNLocMax2[0][matched]->Fill(en,splitFrac);
        if(splitFrac >splitFracMin) fhMassAfterCutsNLocMax2[0][matched]->Fill(en,mass);
        
        if(!matched && IsDataMC() && fFillMCHisto && mcindex==kmcPi0)
        {
          
          fhMCGenFracAfterCutsNLocMax2MCPi0      ->Fill(en   ,  efrac     );
          fhMCGenSplitEFracAfterCutsNLocMax2MCPi0->Fill(en   ,  efracSplit);
        }
      }
            
      if     (pidTag==AliCaloPID::kPhoton)
      {
        fhM02ConNLocMax2 [0][matched]->Fill(en,l0);
        fhMassConNLocMax2[0][matched]->Fill(en,mass);
        fhAsyConNLocMax2 [0][matched]->Fill(en,asym);
      }
      else if(pidTag==AliCaloPID::kPi0   )
      {
        fhM02Pi0NLocMax2 [0][matched]->Fill(en,l0);
        fhMassPi0NLocMax2[0][matched]->Fill(en,mass);
        fhAsyPi0NLocMax2 [0][matched]->Fill(en,asym);
        fhCentralityPi0NLocMax2[0][matched]->Fill(en,cent) ;
        if(!matched)
        {
          fhEventPlanePi0NLocMax2->Fill(en,evp) ;
          if(en > fHistoECut)fhPi0EtaPhiNLocMax2->Fill(eta,phi);
          FillSSWeightHistograms(cluster, 1, absId1, absId2);
          fhPi0EPairDiffTimeNLM2->Fill(e1+e2,t12diff);
        }
      }
      else if(pidTag==AliCaloPID::kEta)
      {
        fhM02EtaNLocMax2 [0][matched]->Fill(en,l0);
        fhMassEtaNLocMax2[0][matched]->Fill(en,mass);
        fhAsyEtaNLocMax2 [0][matched]->Fill(en,asym);
        fhCentralityEtaNLocMax2[0][matched]->Fill(en,cent) ;
        if(!matched)
        {
          fhEventPlaneEtaNLocMax2->Fill(en,evp) ;
          if(en > fHistoECut)fhEtaEtaPhiNLocMax2->Fill(eta,phi);
          fhEtaEPairDiffTimeNLM2->Fill(e1+e2,t12diff);
        }
      }
      
    }
    else if(nMax >2) 
    {
      fhMassNLocMaxN[0][matched]->Fill(en,mass);
      fhAsymNLocMaxN[0][matched]->Fill(en,asym);
      
      // Effect of cuts in mass histograms
      if(!matched)
      {
        if(m02OK)
        {
          fhMassM02CutNLocMaxN->Fill(en,mass);
          fhAsymM02CutNLocMaxN->Fill(en,asym );
          if(splitFrac > splitFracMin)fhMassSplitECutNLocMaxN->Fill(en,mass );
        } // m02
      } // split frac
      
      if(m02OK && asyOK)
      {
        fhSplitEFractionAfterCutsNLocMaxN[0][matched]->Fill(en,splitFrac);
        if(splitFrac > splitFracMin) fhMassAfterCutsNLocMaxN[0][matched]->Fill(en,mass);
        
        if(!matched && IsDataMC() && fFillMCHisto && mcindex==kmcPi0)
        {
          
          fhMCGenFracAfterCutsNLocMaxNMCPi0      ->Fill(en   ,  efrac     );
          fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0->Fill(en   ,  efracSplit);
        }
      }
      
      if     (pidTag==AliCaloPID::kPhoton)
      {
        fhM02ConNLocMaxN [0][matched]->Fill(en,l0);
        fhMassConNLocMaxN[0][matched]->Fill(en,mass);
        fhAsyConNLocMaxN [0][matched]->Fill(en,asym); }
      else if(pidTag==AliCaloPID::kPi0   )
      {
        fhM02Pi0NLocMaxN [0][matched]->Fill(en,l0);
        fhMassPi0NLocMaxN[0][matched]->Fill(en,mass);
        fhAsyPi0NLocMaxN [0][matched]->Fill(en,asym);
        fhCentralityPi0NLocMaxN[0][matched]->Fill(en,cent) ;
        if(!matched)
        {
          fhEventPlanePi0NLocMaxN->Fill(en,evp) ;
          if(en > fHistoECut)fhPi0EtaPhiNLocMaxN->Fill(eta,phi);
          FillSSWeightHistograms(cluster, 2,  absId1, absId2);
          fhPi0EPairDiffTimeNLMN->Fill(e1+e2,t12diff);
        }
      }
      else if(pidTag==AliCaloPID::kEta)
      {
        fhM02EtaNLocMaxN [0][matched]->Fill(en,l0);
        fhMassEtaNLocMaxN[0][matched]->Fill(en,mass);
        fhAsyEtaNLocMaxN [0][matched]->Fill(en,asym);
        fhCentralityEtaNLocMaxN[0][matched]->Fill(en,cent) ;
        if(!matched)
        {
          fhEventPlaneEtaNLocMaxN->Fill(en,evp) ;
          if(en > fHistoECut)fhEtaEtaPhiNLocMaxN->Fill(eta,phi);
          fhEtaEPairDiffTimeNLMN->Fill(e1+e2,t12diff);
        }
      }
    }
    
    if(IsDataMC())
    {
      if     (nMax==1) 
      { 
        fhMassNLocMax1[mcindex][matched]->Fill(en,mass);
        fhAsymNLocMax1[mcindex][matched]->Fill(en,asym);
        
        if(asyOK && m02OK)
        {
          fhSplitEFractionAfterCutsNLocMax1[mcindex][matched]->Fill(en,splitFrac);
          if(splitFrac > splitFracMin)
            fhMassAfterCutsNLocMax1[mcindex][matched]->Fill(en,mass);
        }

        if     (pidTag==AliCaloPID::kPhoton) { fhM02ConNLocMax1[mcindex][matched]->Fill(en,l0); fhMassConNLocMax1[mcindex][matched]->Fill(en,mass); fhAsyConNLocMax1[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kPi0   ) { fhM02Pi0NLocMax1[mcindex][matched]->Fill(en,l0); fhMassPi0NLocMax1[mcindex][matched]->Fill(en,mass); fhAsyPi0NLocMax1[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kEta   ) { fhM02EtaNLocMax1[mcindex][matched]->Fill(en,l0); fhMassEtaNLocMax1[mcindex][matched]->Fill(en,mass); fhAsyEtaNLocMax1[mcindex][matched]->Fill(en,asym); }
        
        if     (pidTag==AliCaloPID::kPi0) fhCentralityPi0NLocMax1[mcindex][matched]->Fill(en,cent) ;
        else if(pidTag==AliCaloPID::kEta) fhCentralityEtaNLocMax1[mcindex][matched]->Fill(en,cent) ;
      }
      else if(nMax==2) 
      {
        fhMassNLocMax2[mcindex][matched]->Fill(en,mass);
        fhAsymNLocMax2[mcindex][matched]->Fill(en,asym);
        
        if(asyOK && m02OK)
        {
          fhSplitEFractionAfterCutsNLocMax2[mcindex][matched]->Fill(en,splitFrac);
          if(splitFrac >splitFracMin)
            fhMassAfterCutsNLocMax2[mcindex][matched]->Fill(en,mass);
        }
        
        if     (pidTag==AliCaloPID::kPhoton) { fhM02ConNLocMax2[mcindex][matched]->Fill(en,l0); fhMassConNLocMax2[mcindex][matched]->Fill(en,mass); fhAsyConNLocMax2[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kPi0   ) { fhM02Pi0NLocMax2[mcindex][matched]->Fill(en,l0); fhMassPi0NLocMax2[mcindex][matched]->Fill(en,mass); fhAsyPi0NLocMax2[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kEta   ) { fhM02EtaNLocMax2[mcindex][matched]->Fill(en,l0); fhMassEtaNLocMax2[mcindex][matched]->Fill(en,mass); fhAsyEtaNLocMax2[mcindex][matched]->Fill(en,asym); } 
       
        if     (pidTag==AliCaloPID::kPi0) fhCentralityPi0NLocMax2[mcindex][matched]->Fill(en,cent) ;
        else if(pidTag==AliCaloPID::kEta) fhCentralityEtaNLocMax2[mcindex][matched]->Fill(en,cent) ;
        
      }
      else if(nMax >2) 
      {
        fhMassNLocMaxN[mcindex][matched]->Fill(en,mass);
        fhAsymNLocMaxN[mcindex][matched]->Fill(en,asym);
        
        if(asyOK && m02OK)
        {
          fhSplitEFractionAfterCutsNLocMaxN[mcindex][matched]->Fill(en,splitFrac);
          if(splitFrac > splitFracMin )
            fhMassAfterCutsNLocMaxN[mcindex][matched]->Fill(en,mass);
        }
        
        if     (pidTag==AliCaloPID::kPhoton) { fhM02ConNLocMaxN[mcindex][matched]->Fill(en,l0); fhMassConNLocMaxN[mcindex][matched]->Fill(en,mass); fhAsyConNLocMaxN[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kPi0   ) { fhM02Pi0NLocMaxN[mcindex][matched]->Fill(en,l0); fhMassPi0NLocMaxN[mcindex][matched]->Fill(en,mass); fhAsyPi0NLocMaxN[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kEta   ) { fhM02EtaNLocMaxN[mcindex][matched]->Fill(en,l0); fhMassEtaNLocMaxN[mcindex][matched]->Fill(en,mass); fhAsyEtaNLocMaxN[mcindex][matched]->Fill(en,asym); }
        
        if     (pidTag==AliCaloPID::kPi0) fhCentralityPi0NLocMaxN[mcindex][matched]->Fill(en,cent) ;
        else if(pidTag==AliCaloPID::kEta) fhCentralityEtaNLocMaxN[mcindex][matched]->Fill(en,cent) ;
      }
      
    }//Work with MC truth first
    
  }//loop
  
  if(GetDebug() > 1) printf("AliAnaInsideClusterInvariantMass::MakeAnalysisFillHistograms() - END \n");  

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
  printf("%2.2f < lambda_0^2 <%2.2f \n",fM02MinCut,fM02MaxCut);
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


