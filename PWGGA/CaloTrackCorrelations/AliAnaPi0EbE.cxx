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

// --- ROOT system ---
#include <TList.h>
#include <TClonesArray.h>
#include <TObjString.h>

// --- Analysis system ---
#include "AliAnaPi0EbE.h"
#include "AliCaloTrackReader.h"
#include "AliIsolationCut.h"
#include "AliNeutralMesonSelection.h"
#include "AliCaloPID.h"
#include "AliMCAnalysisUtils.h"
#include "AliStack.h"
#include "AliFiducialCut.h"
#include "TParticle.h"
#include "AliVCluster.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

/// \cond CLASSIMP
ClassImp(AliAnaPi0EbE) ;
/// \endcond

//____________________________
/// Default constructor.
/// Initialize parameters.
//____________________________
AliAnaPi0EbE::AliAnaPi0EbE() :
AliAnaCaloTrackCorrBaseClass(),
fAnaType(kIMCalo),
fNLMCutMin(-1),                     fNLMCutMax(10),
fFillAllNLMHistograms(0),           fFillTMHisto(0),
fFillSelectClHisto(0),              
fM02MaxCutForIM(0),                 fM02MinCutForIM(0),
fSelectIsolatedDecay(kFALSE),       fSelectPairInIsoCone(0),           
fR(0),                              fIsoCandMinPt(0),
fMinDist(0.),    fMinDist2(0.),     fMinDist3(0.),
fTimeCutMin(-10000),                fTimeCutMax(10000),
fRejectTrackMatch(kTRUE),           fCheckSplitDistToBad(0),      
fFillWeightHistograms(kFALSE),      fFillOnlySimpleSSHisto(1),
fFillEMCALBCHistograms(0),
fInputAODGammaConvName(""),
fMomentum(),  fMomentum1(),  fMomentum2(),
fMomentum12(),fPrimaryMom(), fGrandMotherMom(),
// Histograms
fhPt(0),                            fhE(0),
fhPtEta(0),                         fhPtPhi(0),                         fhEtaPhi(0),
fhEtaPhiEMCALBC0(0),                fhEtaPhiEMCALBC1(0),                fhEtaPhiEMCALBCN(0),
fhTimeTriggerEMCALBC0UMReMatchOpenTime(0),
fhTimeTriggerEMCALBC0UMReMatchCheckNeigh(0),
fhTimeTriggerEMCALBC0UMReMatchBoth(0),
fhPtCentrality(),                   fhPtEventPlane(0),                  fhMCPtCentrality(),
fhPtReject(0),                      fhEReject(0),
fhPtEtaReject(0),                   fhPtPhiReject(0),                   fhEtaPhiReject(0),
fhMass(0),                          fhMassPt(0),                        
fhMassPtMaxPair(0),                 fhMassPtMinPair(0),                 fhMassSplitPt(0),
fhSelectedMass(0),                  fhSelectedMassPt(0),                fhSelectedMassSplitPt(0),
fhMassPtIsoRCut(0),
fhMassNoOverlap(0),                 fhMassPtNoOverlap(0),               fhMassSplitPtNoOverlap(0),
fhSelectedMassNoOverlap(0),         fhSelectedMassPtNoOverlap(0),       fhSelectedMassSplitPtNoOverlap(0),
fhMCPi0PtRecoPtPrim(0),                       fhMCEtaPtRecoPtPrim(0),
fhMCPi0PtRecoPtPrimNoOverlap(0),              fhMCEtaPtRecoPtPrimNoOverlap(0),
fhMCPi0SplitPtRecoPtPrim(0),                  fhMCEtaSplitPtRecoPtPrim(0),
fhMCPi0SplitPtRecoPtPrimNoOverlap(0),         fhMCEtaSplitPtRecoPtPrimNoOverlap(0),
fhMCPi0SelectedPtRecoPtPrim(0),               fhMCEtaSelectedPtRecoPtPrim(0),
fhMCPi0SelectedPtRecoPtPrimNoOverlap(0),      fhMCEtaSelectedPtRecoPtPrimNoOverlap(0),
fhMCPi0SelectedSplitPtRecoPtPrim(0),          fhMCEtaSelectedSplitPtRecoPtPrim(0),
fhMCPi0SelectedSplitPtRecoPtPrimNoOverlap(0), fhMCEtaSelectedSplitPtRecoPtPrimNoOverlap(0),
fhAsymmetry(0),                     fhSelectedAsymmetry(0),
fhSplitE(0),                        fhSplitPt(0),
fhSplitPtEta(0),                    fhSplitPtPhi(0),
fhNLocMaxSplitPt(0),
fhPtDecay(0),

// Shower shape histos
fhPtDispersion(0),                  fhPtLambda0(0),                     fhPtLambda0NoSplitCut(0),
fhPtLambda1(0),                     fhPtLambda0NoTRD(0),                fhPtLambda0FracMaxCellCut(0),
fhPtFracMaxCell(0),                 fhPtFracMaxCellNoTRD(0),
fhPtNCells(0),                      fhPtTime(0),                        fhEPairDiffTime(0),
fhPtDispEta(0),                     fhPtDispPhi(0),
fhPtSumEta(0),                      fhPtSumPhi(0),                      fhPtSumEtaPhi(0),
fhPtDispEtaPhiDiff(0),              fhPtSphericity(0),

// MC histos
fhMCPtDecayLostPairPi0(0),          fhMCPtDecayLostPairEta(0),
fhMCE(),                            fhMCPt(),
fhMCPtPhi(),                        fhMCPtEta(),
fhMCEReject(),                      fhMCPtReject(),
fhMCPi0PtGenRecoFraction(0),        fhMCEtaPtGenRecoFraction(0),
fhMCPi0DecayPt(0),                  fhMCPi0DecayPtFraction(0),
fhMCEtaDecayPt(0),                  fhMCEtaDecayPtFraction(0),
fhMCOtherDecayPt(0),
fhMassPairMCPi0(0),                 fhMassPairMCEta(0),
fhAnglePairMCPi0(0),                fhAnglePairMCEta(0),
fhMCPi0PtOrigin(0x0),               fhMCEtaPtOrigin(0x0),
fhMCNotResonancePi0PtOrigin(0),     fhMCPi0PtStatus(0x0),
fhMCPi0ProdVertex(0),               fhMCEtaProdVertex(0),

// Weight studies
fhECellClusterRatio(0),             fhECellClusterLogRatio(0),
fhEMaxCellClusterRatio(0),          fhEMaxCellClusterLogRatio(0),
fhTrackMatchedDEta(0),              fhTrackMatchedDPhi(0),              fhTrackMatchedDEtaDPhi(0),
fhTrackMatchedDEtaPos(0),           fhTrackMatchedDPhiPos(0),           fhTrackMatchedDEtaDPhiPos(0),
fhTrackMatchedDEtaNeg(0),           fhTrackMatchedDPhiNeg(0),           fhTrackMatchedDEtaDPhiNeg(0),
fhTrackMatchedMCParticlePt(0),
fhTrackMatchedMCParticleDEta(0),    fhTrackMatchedMCParticleDPhi(0),
fhdEdx(0),                          fhEOverP(0),                        fhEOverPNoTRD(0),
// Number of local maxima in cluster
fhNLocMaxPt(0),                     fhNLocMaxPtReject(0),
// PileUp
fhTimePtNoCut(0),                   fhTimePtSPD(0),                     fhTimePtSPDMulti(0),
fhTimeNPileUpVertSPD(0),            fhTimeNPileUpVertTrack(0),
fhTimeNPileUpVertContributors(0),
fhTimePileUpMainVertexZDistance(0), fhTimePileUpMainVertexZDiamond(0),
fhPtNPileUpSPDVtx(0),               fhPtNPileUpTrkVtx(0),
fhPtNPileUpSPDVtxTimeCut(0),        fhPtNPileUpTrkVtxTimeCut(0),
fhPtNPileUpSPDVtxTimeCut2(0),       fhPtNPileUpTrkVtxTimeCut2(0)
{
  for(Int_t i = 0; i < fgkNmcTypes; i++)
  {
    fhMCE              [i] = 0;
    fhMCPt             [i] = 0;
    fhMCPtPhi          [i] = 0;
    fhMCPtEta          [i] = 0;
    fhMCPtCentrality   [i] = 0;
    
    fhMCSplitE         [i] = 0;
    fhMCSplitPt        [i] = 0;
    fhMCSplitPtPhi     [i] = 0;
    fhMCSplitPtEta     [i] = 0;
    
    fhMCNLocMaxPt      [i] = 0;
    fhMCNLocMaxSplitPt [i] = 0;
    fhMCNLocMaxPtReject[i] = 0;
    
    fhMCPtDecay         [i] = 0;
    fhMCPtLambda0       [i] = 0;
    fhMCPtLambda0NoTRD  [i] = 0;
    fhMCPtLambda0FracMaxCellCut[i]= 0;
    fhMCPtFracMaxCell   [i] = 0;
    fhMCPtLambda1       [i] = 0;
    fhMCPtDispersion    [i] = 0;
    
    fhMCPtDispEta       [i] = 0;
    fhMCPtDispPhi       [i] = 0;
    fhMCPtSumEtaPhi     [i] = 0;
    fhMCPtDispEtaPhiDiff[i] = 0;
    fhMCPtSphericity    [i] = 0;
    fhMCPtAsymmetry     [i] = 0;
    
    fhMCMassPt             [i]=0;
    fhMCMassSplitPt        [i]=0;
    fhMCSelectedMassPt     [i]=0;
    fhMCSelectedMassSplitPt[i]=0;
    
    fhMCMassPtNoOverlap             [i]=0;
    fhMCMassSplitPtNoOverlap        [i]=0;
    fhMCSelectedMassPtNoOverlap     [i]=0;
    fhMCSelectedMassSplitPtNoOverlap[i]=0;
    
    for(Int_t j = 0; j < 7; j++)
    {
      fhMCLambda0DispEta    [j][i] = 0;
      fhMCLambda0DispPhi    [j][i] = 0;
      fhMCDispEtaDispPhi    [j][i] = 0;
      fhMCAsymmetryLambda0  [j][i] = 0;
      fhMCAsymmetryDispEta  [j][i] = 0;
      fhMCAsymmetryDispPhi  [j][i] = 0;
    }
  }
  
  for(Int_t j = 0; j < 7; j++)
  {
    fhLambda0DispEta    [j] = 0;
    fhLambda0DispPhi    [j] = 0;
    fhDispEtaDispPhi    [j] = 0;
    fhAsymmetryLambda0  [j] = 0;
    fhAsymmetryDispEta  [j] = 0;
    fhAsymmetryDispPhi  [j] = 0;
    
    fhPtPileUp       [j] = 0;
  }
  
  for(Int_t i = 0; i < 3; i++)
  {
    fhPtLambda0LocMax       [i] = 0;
    fhPtLambda1LocMax       [i] = 0;
    fhPtDispersionLocMax    [i] = 0;
    fhPtDispEtaLocMax       [i] = 0;
    fhPtDispPhiLocMax       [i] = 0;
    fhPtSumEtaPhiLocMax     [i] = 0;
    fhPtDispEtaPhiDiffLocMax[i] = 0;
    fhPtSphericityLocMax    [i] = 0;
    fhPtAsymmetryLocMax     [i] = 0;
    fhMassPtLocMax          [i] = 0;
    fhSelectedMassPtLocMax  [i] = 0;
    for(Int_t ipart = 0; ipart < fgkNmcTypes; ipart++)
    {
      fhMCPtLambda0LocMax     [ipart][i] = 0;
      fhMCSelectedMassPtLocMax[ipart][i] = 0;
    }
    
    fhMCPi0PtRecoPtPrimLocMax             [i] = 0;
    fhMCEtaPtRecoPtPrimLocMax             [i] = 0;
    fhMCPi0SplitPtRecoPtPrimLocMax        [i] = 0;
    fhMCEtaSplitPtRecoPtPrimLocMax        [i] = 0;
    
    fhMCPi0SelectedPtRecoPtPrimLocMax     [i] = 0;
    fhMCEtaSelectedPtRecoPtPrimLocMax     [i] = 0;
    fhMCPi0SelectedSplitPtRecoPtPrimLocMax[i] = 0;
    fhMCEtaSelectedSplitPtRecoPtPrimLocMax[i] = 0;
    
  }
  
  // Weight studies
  for(Int_t i =0; i < 14; i++){
    fhLambda0ForW0[i] = 0;
    //fhLambda1ForW0[i] = 0;
    if(i<8)fhMassPairLocMax[i] = 0;
  }
  
  for(Int_t i = 0; i < 11; i++)
  {
    fhEtaPhiTriggerEMCALBC       [i] = 0 ;
    fhTimeTriggerEMCALBC         [i] = 0 ;
    fhTimeTriggerEMCALBCPileUpSPD[i] = 0 ;
    
    fhEtaPhiTriggerEMCALBCUM     [i] = 0 ;
    fhTimeTriggerEMCALBCUM       [i] = 0 ;
    
  }
  
  for(Int_t iSM = 0; iSM < 22; iSM++)
  {
    fhNLocMaxPtSM[iSM] = 0;
    for(Int_t inlm = 0; inlm < 3; inlm++)
    {
      fhSelectedMassPtLocMaxSM    [inlm][iSM] = 0;
      fhSelectedLambda0PtLocMaxSM [inlm][iSM] = 0;
    }
  }
  
  // Initialize parameters
  InitParameters();
  
}

//______________________________________________________________________________________________
/// EMCal trigger cluster Bunch Crossing time peak studies.
/// Fill timing histograms for different event timing peaks,
/// main BC (0), next BC (1), or any other (N).
/// Definitions of BC in AliCaloTrackReader (check).
//______________________________________________________________________________________________
void AliAnaPi0EbE::FillEMCALBCHistograms(Float_t energy, Float_t eta, Float_t phi, Float_t time)
{
  Int_t id = GetReader()->GetTriggerClusterId();
  if( id < 0 ) return;
  
  Int_t bc = GetReader()->GetTriggerClusterBC();
  if(TMath::Abs(bc) >= 6) AliInfo(Form("Trigger BC not expected = %d",bc));
  
  if(phi < 0) phi+=TMath::TwoPi();
  
  if(energy > 2)
  {
    Double_t timeUS = TMath::Abs(time);
    
    if      (timeUS < 25) fhEtaPhiEMCALBC0->Fill(eta, phi, GetEventWeight());
    else if (timeUS < 75) fhEtaPhiEMCALBC1->Fill(eta, phi, GetEventWeight());
    else                  fhEtaPhiEMCALBCN->Fill(eta, phi, GetEventWeight());
  }
  
  if(TMath::Abs(bc) >= 6) return ;
  
  if(GetReader()->IsBadCellTriggerEvent() || GetReader()->IsExoticEvent()) return ;
  
  if(GetReader()->IsTriggerMatched())
  {
    if(energy > 2)
        fhEtaPhiTriggerEMCALBC[bc+5]->Fill(eta, phi, GetEventWeight());
      
    fhTimeTriggerEMCALBC[bc+5]->Fill(energy, time, GetEventWeight());
      
    if(GetReader()->IsPileUpFromSPD())
        fhTimeTriggerEMCALBCPileUpSPD[bc+5]->Fill(energy, time, GetEventWeight());
  }
  else
  {
    if(energy > 2)
        fhEtaPhiTriggerEMCALBCUM[bc+5]->Fill(eta, phi, GetEventWeight());
      
    fhTimeTriggerEMCALBCUM[bc+5]->Fill(energy, time, GetEventWeight());
    
    if(bc==0)
    {
      if(GetReader()->IsTriggerMatchedOpenCuts(0))
          fhTimeTriggerEMCALBC0UMReMatchOpenTime   ->Fill(energy, time, GetEventWeight());
        
      if(GetReader()->IsTriggerMatchedOpenCuts(1))
          fhTimeTriggerEMCALBC0UMReMatchCheckNeigh ->Fill(energy, time, GetEventWeight());
        
      if(GetReader()->IsTriggerMatchedOpenCuts(2))
          fhTimeTriggerEMCALBC0UMReMatchBoth       ->Fill(energy, time, GetEventWeight());
    }
  }
}

//___________________________________________________________________________________
/// Fill some histograms to understand pile-up, depending on SPD/EMCal pile-up
/// definitions. EMCal definition depends on timing of clusters in event.
//___________________________________________________________________________________
void AliAnaPi0EbE::FillPileUpHistograms(Float_t pt, Float_t time, AliVCluster * calo)
{
  //printf("E %f, time %f\n",energy,time);
  AliVEvent * event = GetReader()->GetInputEvent();
  
  fhTimePtNoCut->Fill(pt, time, GetEventWeight());
    
  if(GetReader()->IsPileUpFromSPD())             { fhPtPileUp[0]->Fill(pt, GetEventWeight()); fhTimePtSPD->Fill(pt, time, GetEventWeight()); }
  if(GetReader()->IsPileUpFromEMCal())             fhPtPileUp[1]->Fill(pt, GetEventWeight());
  if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtPileUp[2]->Fill(pt, GetEventWeight());
  if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtPileUp[3]->Fill(pt, GetEventWeight());
  if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtPileUp[4]->Fill(pt, GetEventWeight());
  if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtPileUp[5]->Fill(pt, GetEventWeight());
  if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtPileUp[6]->Fill(pt, GetEventWeight());
  
  if(event->IsPileupFromSPDInMultBins()) fhTimePtSPDMulti->Fill(pt, time, GetEventWeight());
  
  // cells in cluster
  
  AliVCaloCells* cells = 0;
  if(GetCalorimeter() == kEMCAL) cells = GetEMCALCells();
  else                           cells = GetPHOSCells();
  
  Float_t maxCellFraction = 0.;
  Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(cells,calo,maxCellFraction);
  
  Double_t tmax  = cells->GetCellTime(absIdMax);
  GetCaloUtils()->RecalibrateCellTime(tmax, GetCalorimeter(), absIdMax,GetReader()->GetInputEvent()->GetBunchCrossNumber());
  tmax*=1.e9;
  
  // Loop on cells inside cluster, max cell must be over 100 MeV and time in BC=0
  if(cells->GetCellAmplitude(absIdMax) > 0.1 && TMath::Abs(tmax) < 30)
  {
    for (Int_t ipos = 0; ipos < calo->GetNCells(); ipos++)
    {
      Int_t absId  = calo->GetCellsAbsId()[ipos];
      
      if( absId == absIdMax ) continue ;
      
      Double_t timecell  = cells->GetCellTime(absId);
      Float_t  amp       = cells->GetCellAmplitude(absId);
      Int_t    bc        = GetReader()->GetInputEvent()->GetBunchCrossNumber();
      GetCaloUtils()->GetEMCALRecoUtils()->AcceptCalibrateCell(absId,bc,amp,timecell,cells);
      timecell*=1e9;
      
      Float_t diff = (tmax-timecell);
      
      if( cells->GetCellAmplitude(absIdMax) < 0.1 ) continue ;
      
      if(GetReader()->IsPileUpFromSPD())
      {
        fhPtCellTimePileUp[0]->Fill(pt, timecell, GetEventWeight());
        fhPtTimeDiffPileUp[0]->Fill(pt, diff    , GetEventWeight());
      }
      
      if(GetReader()->IsPileUpFromEMCal())
      {
        fhPtCellTimePileUp[1]->Fill(pt, timecell, GetEventWeight());
        fhPtTimeDiffPileUp[1]->Fill(pt, diff    , GetEventWeight());
      }
      
      if(GetReader()->IsPileUpFromSPDOrEMCal())
      {
        fhPtCellTimePileUp[2]->Fill(pt, timecell, GetEventWeight());
        fhPtTimeDiffPileUp[2]->Fill(pt, diff    , GetEventWeight());
      }
      
      if(GetReader()->IsPileUpFromSPDAndEMCal())
      {
        fhPtCellTimePileUp[3]->Fill(pt, timecell, GetEventWeight());
        fhPtTimeDiffPileUp[3]->Fill(pt, diff    , GetEventWeight());
      }
      
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())
      {
        fhPtCellTimePileUp[4]->Fill(pt, timecell, GetEventWeight());
        fhPtTimeDiffPileUp[4]->Fill(pt, diff    , GetEventWeight());
      }
      
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())
      {
        fhPtCellTimePileUp[5]->Fill(pt, timecell, GetEventWeight());
        fhPtTimeDiffPileUp[5]->Fill(pt, diff    , GetEventWeight());
      }
      
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal())
      {
        fhPtCellTimePileUp[6]->Fill(pt, timecell, GetEventWeight());
        fhPtTimeDiffPileUp[6]->Fill(pt, diff    , GetEventWeight());
      }
    }//loop
  }
  
  if(pt < 8) return; // Fill time figures for high energy clusters not too close to trigger threshold
  
  AliESDEvent* esdEv = dynamic_cast<AliESDEvent*> (event);
  AliAODEvent* aodEv = dynamic_cast<AliAODEvent*> (event);
  
  // N pile up vertices
  Int_t nVtxSPD = -1;
  Int_t nVtxTrk = -1;
  
  if      (esdEv)
  {
    nVtxSPD = esdEv->GetNumberOfPileupVerticesSPD();
    nVtxTrk = esdEv->GetNumberOfPileupVerticesTracks();
    
  } // ESD
  else if (aodEv)
  {
    nVtxSPD = aodEv->GetNumberOfPileupVerticesSPD();
    nVtxTrk = aodEv->GetNumberOfPileupVerticesTracks();
  } // AOD
  
  fhTimeNPileUpVertSPD  ->Fill(time, nVtxSPD, GetEventWeight());
  fhTimeNPileUpVertTrack->Fill(time, nVtxTrk, GetEventWeight());
  
  fhPtNPileUpSPDVtx->Fill(pt, nVtxSPD, GetEventWeight());
  fhPtNPileUpTrkVtx->Fill(pt, nVtxTrk, GetEventWeight());
  
  if(TMath::Abs(time) < 25)
  {
    fhPtNPileUpSPDVtxTimeCut ->Fill(pt, nVtxSPD, GetEventWeight());
    fhPtNPileUpTrkVtxTimeCut ->Fill(pt, nVtxTrk, GetEventWeight());
  }
  
  if(time < 75 && time > -25)
  {
    fhPtNPileUpSPDVtxTimeCut2->Fill(pt, nVtxSPD, GetEventWeight());
    fhPtNPileUpTrkVtxTimeCut2->Fill(pt, nVtxTrk, GetEventWeight());
  }
  
  //printf("Is SPD %d, Is SPD Multi %d, n spd %d, n track %d\n",
  //       GetReader()->IsPileUpFromSPD(),event->IsPileupFromSPDInMultBins(),nVtxSPD,nVtxTracks);
  
  Int_t ncont = -1;
  Float_t z1 = -1, z2 = -1;
  Float_t diamZ = -1;
  for(Int_t iVert=0; iVert<nVtxSPD;iVert++)
  {
    if      (esdEv)
    {
      const AliESDVertex* pv=esdEv->GetPileupVertexSPD(iVert);
      ncont=pv->GetNContributors();
      z1 = esdEv->GetPrimaryVertexSPD()->GetZ();
      z2 = pv->GetZ();
      diamZ = esdEv->GetDiamondZ();
    }//ESD
    else if (aodEv)
    {
      AliAODVertex *pv=aodEv->GetVertex(iVert);
      if(pv->GetType()!=AliAODVertex::kPileupSPD) continue;
      ncont=pv->GetNContributors();
      z1=aodEv->GetPrimaryVertexSPD()->GetZ();
      z2=pv->GetZ();
      diamZ = aodEv->GetDiamondZ();
    } // AOD
    
    Double_t distZ  = TMath::Abs(z2-z1);
    diamZ  = TMath::Abs(z2-diamZ);
    
    fhTimeNPileUpVertContributors  ->Fill(time, ncont, GetEventWeight());
    fhTimePileUpMainVertexZDistance->Fill(time, distZ, GetEventWeight());
    fhTimePileUpMainVertexZDiamond ->Fill(time, diamZ, GetEventWeight());
    
  } // vertex loop
}

//__________________________________________________________________________
/// Fill histograms for clusters that do not pass the identification
/// (kSSCalo case only).
//__________________________________________________________________________
void AliAnaPi0EbE::FillRejectedClusterHistograms(Int_t mctag, Int_t nMaxima)
{
  Float_t ener  = fMomentum.E();
  Float_t pt    = fMomentum.Pt();
    
  Float_t phi   = fMomentum.Phi();
  if(phi < 0) phi+=TMath::TwoPi();
  
  Float_t eta = fMomentum.Eta();
  
  fhPtReject     ->Fill(pt  , GetEventWeight());
  fhEReject      ->Fill(ener, GetEventWeight());
  
  fhPtEtaReject  ->Fill(pt , eta, GetEventWeight());
  fhPtPhiReject  ->Fill(pt , phi, GetEventWeight());
  fhEtaPhiReject ->Fill(eta, phi, GetEventWeight());
    
  fhNLocMaxPtReject->Fill(pt, nMaxima, GetEventWeight());
  
  if(IsDataMC())
  {
    Int_t mcIndex = GetMCIndex(mctag);
    fhMCEReject  [mcIndex] ->Fill(ener, GetEventWeight());
    fhMCPtReject [mcIndex] ->Fill(pt  , GetEventWeight());
      
    if(fFillAllNLMHistograms) fhMCNLocMaxPtReject[mcIndex]->Fill(pt, nMaxima, GetEventWeight());
  }
}

//_______________________________________________________________________________________________
/// Fill shower shape, timing and other histograms for selected clusters.
//_______________________________________________________________________________________________
void AliAnaPi0EbE::FillSelectedClusterHistograms(AliVCluster* cluster, Float_t pt, Int_t nMaxima,
                                                 Int_t tag, Float_t asy)
{
  Float_t ener = cluster->E();
  Float_t disp = cluster->GetDispersion()*cluster->GetDispersion();
  Float_t l0   = cluster->GetM02();
  Float_t l1   = cluster->GetM20();
  Int_t   nSM  = GetModuleNumber(cluster);
  
  Int_t ptbin = -1;
  if      (pt < 2 ) ptbin = 0;
  else if (pt < 4 ) ptbin = 1;
  else if (pt < 6 ) ptbin = 2;
  else if (pt < 10) ptbin = 3;
  else if (pt < 15) ptbin = 4;
  else if (pt < 20) ptbin = 5;
  else              ptbin = 6;
  
  Int_t indexMax = -1;
  if     (nMaxima==1) indexMax = 0 ;
  else if(nMaxima==2) indexMax = 1 ;
  else                indexMax = 2 ;
  
  FillWeightHistograms(cluster);
  
  fhPtLambda0->Fill(pt, l0, GetEventWeight());
  fhPtLambda1->Fill(pt, l1, GetEventWeight());
  
  fhNLocMaxPt->Fill(pt, nMaxima, GetEventWeight());
  
  if(fFillAllNLMHistograms)
  {
    if(nSM < fNModules && nSM >=0)
      fhNLocMaxPtSM[nSM]->Fill(pt, nMaxima, GetEventWeight());
    
    fhPtLambda0LocMax   [indexMax]->Fill(pt, l0, GetEventWeight());
    fhPtLambda1LocMax   [indexMax]->Fill(pt, l1, GetEventWeight());
  }
  
  Float_t ll0  = 0., ll1  = 0.;
  Float_t dispp= 0., dEta = 0., dPhi    = 0.;
  Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;
  AliVCaloCells * cell = 0x0;
  Float_t maxCellFraction = 0;
  
  if(GetCalorimeter() == kEMCAL && !fFillOnlySimpleSSHisto)
  {
    cell = GetEMCALCells();
    
    GetCaloUtils()->GetMaxEnergyCell(cell, cluster, maxCellFraction);
    fhPtFracMaxCell->Fill(pt, maxCellFraction, GetEventWeight());
    
    if(maxCellFraction < 0.5)
      fhPtLambda0FracMaxCellCut->Fill(pt, l0, GetEventWeight());
    
    GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(),
                                                                                 cell, cluster,
                                                                                 ll0, ll1, dispp, dEta, dPhi,
                                                                                 sEta, sPhi, sEtaPhi);
    fhPtDispersion    -> Fill(pt, disp     , GetEventWeight());
    fhPtDispEta       -> Fill(pt, dEta     , GetEventWeight());
    fhPtDispPhi       -> Fill(pt, dPhi     , GetEventWeight());
    fhPtSumEta        -> Fill(pt, sEta     , GetEventWeight());
    fhPtSumPhi        -> Fill(pt, sPhi     , GetEventWeight());
    fhPtSumEtaPhi     -> Fill(pt, sEtaPhi  , GetEventWeight());
    fhPtDispEtaPhiDiff-> Fill(pt, dPhi-dEta, GetEventWeight());
    if(dEta+dPhi>0)
        fhPtSphericity-> Fill(pt, (dPhi-dEta)/(dEta+dPhi), GetEventWeight());
    
    fhDispEtaDispPhi[ptbin]->Fill(dEta, dPhi, GetEventWeight());
    fhLambda0DispEta[ptbin]->Fill(l0  , dEta, GetEventWeight());
    fhLambda0DispPhi[ptbin]->Fill(l0  , dPhi, GetEventWeight());
    
    if (fAnaType==kSSCalo)
    {
      // Asymmetry histograms
      fhAsymmetryLambda0[ptbin]->Fill(l0  , asy, GetEventWeight());
      fhAsymmetryDispEta[ptbin]->Fill(dEta, asy, GetEventWeight());
      fhAsymmetryDispPhi[ptbin]->Fill(dPhi, asy, GetEventWeight());
    }
    
    if(fFillAllNLMHistograms)
    {
      fhPtDispersionLocMax    [indexMax]->Fill(pt, disp     , GetEventWeight());
      fhPtDispEtaLocMax       [indexMax]->Fill(pt, dEta     , GetEventWeight());
      fhPtDispPhiLocMax       [indexMax]->Fill(pt, dPhi     , GetEventWeight());
      fhPtSumEtaPhiLocMax     [indexMax]->Fill(pt, sEtaPhi  , GetEventWeight());
      fhPtDispEtaPhiDiffLocMax[indexMax]->Fill(pt, dPhi-dEta, GetEventWeight());
      if(dEta+dPhi>0)
          fhPtSphericityLocMax[indexMax]->Fill(pt, (dPhi-dEta)/(dEta+dPhi), GetEventWeight());
      if(fAnaType==kSSCalo)
          fhPtAsymmetryLocMax [indexMax]->Fill(pt, asy, GetEventWeight());
    }
  }
  
  if(GetCalorimeter()==kEMCAL &&  GetFirstSMCoveredByTRD() >= 0 &&
     GetModuleNumber(cluster) < GetFirstSMCoveredByTRD() )
  {
    fhPtLambda0NoTRD    ->Fill(pt, l0, GetEventWeight());
    if(!fFillOnlySimpleSSHisto)
      fhPtFracMaxCellNoTRD->Fill(pt, maxCellFraction, GetEventWeight());
  }
  
  fhPtTime  ->Fill(pt, cluster->GetTOF()*1.e9, GetEventWeight());
  fhPtNCells->Fill(pt, cluster->GetNCells()  , GetEventWeight());
  
  // Fill Track matching control histograms
  if(fFillTMHisto)
  {
    Float_t dZ  = cluster->GetTrackDz();
    Float_t dR  = cluster->GetTrackDx();
    
//    if(cluster->IsEMCAL() && GetCaloUtils()->IsRecalculationOfClusterTrackMatchingOn())
//    {
//      dR = 2000., dZ = 2000.;
//      GetCaloUtils()->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dZ,dR);
//    }
    //printf("Pi0EbE: dPhi %f, dEta %f\n",dR,dZ);
    
    AliVTrack *track = GetCaloUtils()->GetMatchedTrack(cluster, GetReader()->GetInputEvent());
    
    Bool_t positive = kFALSE;
    if(track) positive = (track->Charge()>0);
    
    if(fhTrackMatchedDEta && TMath::Abs(dR) < 999)
    {
      fhTrackMatchedDEta->Fill(pt, dZ, GetEventWeight());
      fhTrackMatchedDPhi->Fill(pt, dR, GetEventWeight());
      if(ener > 0.5)
          fhTrackMatchedDEtaDPhi->Fill(dZ, dR, GetEventWeight());
      
      if(track)
      {
        if(positive)
        {
          fhTrackMatchedDEtaPos->Fill(pt, dZ, GetEventWeight());
          fhTrackMatchedDPhiPos->Fill(pt, dR, GetEventWeight());
          if(ener > 0.5)
              fhTrackMatchedDEtaDPhiPos->Fill(dZ, dR, GetEventWeight());
        }
        else
        {
          fhTrackMatchedDEtaNeg->Fill(pt, dZ, GetEventWeight());
          fhTrackMatchedDPhiNeg->Fill(pt, dR, GetEventWeight());
          if(ener > 0.5)
              fhTrackMatchedDEtaDPhiNeg->Fill(dZ, dR, GetEventWeight());
        }
      }
    }
      
    // Check dEdx and E/p of matched clusters
    
    if(TMath::Abs(dZ) < 0.05 && TMath::Abs(dR) < 0.05)
    {
      if(track)
      {
        Float_t dEdx = track->GetTPCsignal();
        fhdEdx->Fill(pt, dEdx, GetEventWeight());
        
        Float_t eOverp = cluster->E()/track->P();
        fhEOverP->Fill(pt, eOverp, GetEventWeight());
        
        // Change nSM for year > 2011 (< 4 in 2012-13, none after)
        if(GetCalorimeter()==kEMCAL &&  GetFirstSMCoveredByTRD() >= 0 &&
           GetModuleNumber(cluster) < GetFirstSMCoveredByTRD() )
          fhEOverPNoTRD->Fill(pt, eOverp, GetEventWeight());
        
      }
      //else
      //  printf("AliAnaPi0EbE::FillSelectedClusterHistograms() - Residual OK but (dR, dZ)= (%2.4f,%2.4f) no track associated WHAT? \n", dR,dZ);
      
      if(IsDataMC())
      {
        Float_t mctag = -1;
        if  ( !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)  )
        {
          if       ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)      ||
                     GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)       ) mctag =  2.5 ;
          else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton)    ) mctag =  0.5 ;
          else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron)  ) mctag =  1.5 ;
          else                                                                                 mctag =  3.5 ;
          
        }
        else
        {
          if       ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)      ||
                     GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)       ) mctag =  6.5 ;
          else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton)    ) mctag =  4.5 ;
          else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron)  ) mctag =  5.5 ;
          else                                                                                 mctag =  7.5 ;
        }
        
        fhTrackMatchedMCParticlePt  ->Fill(pt, mctag, GetEventWeight());
        fhTrackMatchedMCParticleDEta->Fill(dZ, mctag, GetEventWeight());
        fhTrackMatchedMCParticleDPhi->Fill(dR, mctag, GetEventWeight());
        
      }  // MC
    }
  } // Track matching histograms
  
  if(IsDataMC())
  {
    Int_t mcIndex = GetMCIndex(tag);
    
    fhMCPtLambda0[mcIndex]    ->Fill(pt, l0, GetEventWeight());
    fhMCPtLambda1[mcIndex]    ->Fill(pt, l1, GetEventWeight());
    if(fFillAllNLMHistograms)
        fhMCPtLambda0LocMax[mcIndex][indexMax]->Fill(pt, l0, GetEventWeight());
    
    if(GetCalorimeter()==kEMCAL &&
       GetFirstSMCoveredByTRD() >= 0 &&
       GetModuleNumber(cluster) <  GetFirstSMCoveredByTRD() )
      fhMCPtLambda0NoTRD[mcIndex]->Fill(pt, l0, GetEventWeight());
    
    if(GetCalorimeter() == kEMCAL && !fFillOnlySimpleSSHisto)
    {
      if(maxCellFraction < 0.5)
        fhMCPtLambda0FracMaxCellCut[mcIndex]->Fill(pt, l0, GetEventWeight());
      
      fhMCPtFracMaxCell    [mcIndex]->Fill(pt, maxCellFraction, GetEventWeight());
        
      fhMCPtDispersion     [mcIndex]->Fill(pt, disp     , GetEventWeight());
      fhMCPtDispEta        [mcIndex]->Fill(pt, dEta     , GetEventWeight());
      fhMCPtDispPhi        [mcIndex]->Fill(pt, dPhi     , GetEventWeight());
      fhMCPtSumEtaPhi      [mcIndex]->Fill(pt, sEtaPhi  , GetEventWeight());
      fhMCPtDispEtaPhiDiff [mcIndex]->Fill(pt, dPhi-dEta, GetEventWeight());
      if(dEta+dPhi > 0)
          fhMCPtSphericity[mcIndex]-> Fill(pt, (dPhi-dEta)/(dEta+dPhi), GetEventWeight());
      
      if (fAnaType==kSSCalo)
      {
        fhMCAsymmetryLambda0[ptbin][mcIndex]->Fill(l0  , asy, GetEventWeight());
        fhMCAsymmetryDispEta[ptbin][mcIndex]->Fill(dEta, asy, GetEventWeight());
        fhMCAsymmetryDispPhi[ptbin][mcIndex]->Fill(dPhi, asy, GetEventWeight());
      }
      
      fhMCDispEtaDispPhi[ptbin][mcIndex]->Fill(dEta, dPhi, GetEventWeight());
      fhMCLambda0DispEta[ptbin][mcIndex]->Fill(l0  , dEta, GetEventWeight());
      fhMCLambda0DispPhi[ptbin][mcIndex]->Fill(l0  , dPhi, GetEventWeight());
    } // only SS simple?
  } // MC
}

//________________________________________________________
/// Calculate cluster energy weights and fill histograms.
//________________________________________________________
void AliAnaPi0EbE::FillWeightHistograms(AliVCluster *clus)
{
  if(!fFillWeightHistograms || GetMixedEvent()) return;
  
  AliVCaloCells* cells = 0;
  if(GetCalorimeter() == kEMCAL) cells = GetEMCALCells();
  else                        cells = GetPHOSCells();
  
  // First recalculate energy in case non linearity was applied
  Float_t  energy = 0;
  Float_t  ampMax = 0;
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++)
  {
    
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    // Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp,GetCalorimeter(), id);
    
    energy    += amp;
    
    if(amp> ampMax)
      ampMax = amp;
    
  } // energy loop
  
  if(energy <=0 )
  {
    AliInfo(Form("Wrong calculated energy %f",energy));
    return;
  }
  
  fhEMaxCellClusterRatio   ->Fill(energy, ampMax/energy            , GetEventWeight());
  fhEMaxCellClusterLogRatio->Fill(energy, TMath::Log(ampMax/energy), GetEventWeight());
  
  // Get the ratio and log ratio to all cells in cluster
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++)
  {
    Int_t id = clus->GetCellsAbsId()[ipos];
    
    // Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp,GetCalorimeter(), id);
    
    fhECellClusterRatio   ->Fill(energy, amp/energy            , GetEventWeight());
    fhECellClusterLogRatio->Fill(energy, TMath::Log(amp/energy), GetEventWeight());
  }
  
  // Recalculate shower shape for different W0
  if ( GetCalorimeter() == kEMCAL )
  {
    Float_t l0org = clus->GetM02();
    Float_t l1org = clus->GetM20();
    Float_t dorg  = clus->GetDispersion();
    
    for(Int_t iw = 0; iw < 14; iw++)
    {
      GetCaloUtils()->GetEMCALRecoUtils()->SetW0(1+iw*0.5);
      GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), cells, clus);
      
      fhLambda0ForW0[iw]->Fill(energy, clus->GetM02(), GetEventWeight());
//    fhLambda1ForW0[iw]->Fill(energy, clus->GetM20(), GetEventWeight());
      
    } // w0 loop
    
    // Set the original values back
    clus->SetM02(l0org);
    clus->SetM20(l1org);
    clus->SetDispersion(dorg);
    
  } // EMCAL
}

//__________________________________________
/// Save parameters used for analysis.
//__________________________________________
TObjString * AliAnaPi0EbE::GetAnalysisCuts()
{
  TString parList ; // this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaPi0EbE ---\n") ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fAnaType=%d (selection type) \n",fAnaType) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Calorimeter: %s;",GetCalorimeterString().Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Local maxima in cluster: %d < nlm < %d;",fNLMCutMin,fNLMCutMax) ;
  parList+=onePar ;
  
  if(fAnaType == kSSCalo)
  {
    snprintf(onePar,buffersize,"E cut: %2.2f<E<%2.2f;",GetMinEnergy(),GetMaxEnergy()) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"N cell cut: N > %d;",GetCaloPID()->GetClusterSplittingMinNCells()) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"Min Dist to Bad channel: fMinDist =%2.2f; fMinDist2=%2.2f, fMinDist3=%2.2f;",fMinDist, fMinDist2,fMinDist3) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"Min E cut for NLM cases: 1) %2.2f; 2) %2.2f; 3) %2.2f;",fNLMECutMin[0],fNLMECutMin[1],fNLMECutMin[2]) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"Reject Matched tracks?: %d;",fRejectTrackMatch) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"Reject split cluster close to border or bad?: %d;",fCheckSplitDistToBad) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"Time cut: %2.2f<t<%2.2f;",fTimeCutMin,fTimeCutMax) ;
    parList+=onePar ;
    
    // Get parameters set in PID class.
    parList += GetCaloPID()->GetPIDParametersList() ;
  }
  else if(fAnaType == kIMCalo || fAnaType == kIMCaloTracks)
  {
    snprintf(onePar,buffersize,"Select %s;", (GetNeutralMesonSelection()->GetParticle()).Data()) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"Mass cut: %2.2f<M<%2.2f;",GetNeutralMesonSelection()->GetInvMassMinCut() ,GetNeutralMesonSelection()->GetInvMassMaxCut()) ;
    parList+=onePar ;
  }
  else if(fAnaType == kIMCaloTracks)
  {
    snprintf(onePar,buffersize,"Photon Conv Array: %s;",fInputAODGammaConvName.Data()) ;
    parList+=onePar ;
  }
  else if(fAnaType == kIMCalo)
  {
    snprintf(onePar,buffersize,"Time Diff: %2.2f;",GetPairTimeCut()) ;
    parList+=onePar ;
  }
  
  // Get parameters set in base class.
  //parList += GetBaseParametersList() ;
  
  return new TObjString(parList) ;
}

//_____________________________________________
// Create histograms to be saved in output file and
// store them in outputContainer.
//_____________________________________________
TList *  AliAnaPi0EbE::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ;
  outputContainer->SetName("Pi0EbEHistos") ;
  
  Int_t nptbins  = GetHistogramRanges()->GetHistoPtBins();           Float_t ptmax  = GetHistogramRanges()->GetHistoPtMax();           Float_t ptmin  = GetHistogramRanges()->GetHistoPtMin();
  Int_t nphibins = GetHistogramRanges()->GetHistoPhiBins();          Float_t phimax = GetHistogramRanges()->GetHistoPhiMax();          Float_t phimin = GetHistogramRanges()->GetHistoPhiMin();
  Int_t netabins = GetHistogramRanges()->GetHistoEtaBins();          Float_t etamax = GetHistogramRanges()->GetHistoEtaMax();          Float_t etamin = GetHistogramRanges()->GetHistoEtaMin();
  Int_t ssbins   = GetHistogramRanges()->GetHistoShowerShapeBins();  Float_t ssmax  = GetHistogramRanges()->GetHistoShowerShapeMax();  Float_t ssmin  = GetHistogramRanges()->GetHistoShowerShapeMin();
  Int_t tdbins   = GetHistogramRanges()->GetHistoDiffTimeBins() ;    Float_t tdmax  = GetHistogramRanges()->GetHistoDiffTimeMax();     Float_t tdmin  = GetHistogramRanges()->GetHistoDiffTimeMin();
  Int_t tbins    = GetHistogramRanges()->GetHistoTimeBins() ;        Float_t tmax   = GetHistogramRanges()->GetHistoTimeMax();         Float_t tmin   = GetHistogramRanges()->GetHistoTimeMin();
  Int_t nbins    = GetHistogramRanges()->GetHistoNClusterCellBins(); Int_t   nmax   = GetHistogramRanges()->GetHistoNClusterCellMax(); Int_t   nmin   = GetHistogramRanges()->GetHistoNClusterCellMin();
  
  Int_t   nmassbins   = GetHistogramRanges()->GetHistoMassBins();
  Float_t massmin     = GetHistogramRanges()->GetHistoMassMin();
  Float_t massmax     = GetHistogramRanges()->GetHistoMassMax();
  
  Int_t   nresetabins = GetHistogramRanges()->GetHistoTrackResidualEtaBins();
  Float_t resetamax   = GetHistogramRanges()->GetHistoTrackResidualEtaMax();
  Float_t resetamin   = GetHistogramRanges()->GetHistoTrackResidualEtaMin();
  Int_t   nresphibins = GetHistogramRanges()->GetHistoTrackResidualPhiBins();
  Float_t resphimax   = GetHistogramRanges()->GetHistoTrackResidualPhiMax();
  Float_t resphimin   = GetHistogramRanges()->GetHistoTrackResidualPhiMin();
  
  Int_t   ndedxbins   = GetHistogramRanges()->GetHistodEdxBins();
  Float_t dedxmax     = GetHistogramRanges()->GetHistodEdxMax();
  Float_t dedxmin     = GetHistogramRanges()->GetHistodEdxMin();
  Int_t   nPoverEbins = GetHistogramRanges()->GetHistoPOverEBins();
  Float_t pOverEmax   = GetHistogramRanges()->GetHistoPOverEMax();
  Float_t pOverEmin   = GetHistogramRanges()->GetHistoPOverEMin();
  
  Int_t   ntimptbins  = GetHistogramRanges()->GetHistoTimeBins();
  Float_t timemax     = GetHistogramRanges()->GetHistoTimeMax();
  Float_t timemin     = GetHistogramRanges()->GetHistoTimeMin();
  
  // Init the number of modules, set in the class AliCalorimeterUtils
  //
  InitCaloParameters(); // See AliCaloTrackCorrBaseClass
  
  TString nlm[]   = {"1 Local Maxima","2 Local Maxima", "NLM > 2"};
  
  TString ptype [] = {"#pi^{0}", "#eta", "#gamma (direct)","#gamma (#pi^{0})", "#gamma (#eta)", "#gamma (other)",  "e^{#pm}"  , "hadron/other combinations"};
  TString pname [] = {"Pi0"    , "Eta" , "Photon"         ,"Pi0Decay"        , "EtaDecay"     , "OtherDecay"    ,   "Electron", "Hadron"};
  
  Int_t   bin[]   = {0,2,4,6,10,15,20,100}; // energy bins
  
  fhPt  = new TH1F("hPt","Number of identified  #pi^{0} (#eta) decay",nptbins,ptmin,ptmax);
  fhPt->SetYTitle("#it{N}");
  fhPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPt) ;
  
  fhE  = new TH1F("hE","Number of identified  #pi^{0} (#eta) decay pairs",nptbins,ptmin,ptmax);
  fhE->SetYTitle("#it{N}");
  fhE->SetXTitle("#it{E} (GeV)");
  outputContainer->Add(fhE) ;
  
  fhPtPhi  = new TH2F
  ("hPtPhi","Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #phi",nptbins,ptmin,ptmax, nphibins,phimin,phimax);
  fhPtPhi->SetYTitle("#phi (rad)");
  fhPtPhi->SetXTitle("#it{E} (GeV)");
  outputContainer->Add(fhPtPhi) ;
  
  fhPtEta  = new TH2F
  ("hPtEta","Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
  fhPtEta->SetYTitle("#eta");
  fhPtEta->SetXTitle("#it{E} (GeV)");
  outputContainer->Add(fhPtEta) ;
  
  fhEtaPhi  = new TH2F
  ("hEtaPhi","Selected #pi^{0} (#eta) pairs: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhi->SetYTitle("#phi (rad)");
  fhEtaPhi->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhi) ;
  
  if(GetCalorimeter()==kEMCAL && fFillEMCALBCHistograms)
  {
    fhEtaPhiEMCALBC0  = new TH2F
    ("hEtaPhiEMCALBC0","cluster, #it{E} > 2 GeV, #eta vs #phi, for clusters with |#it{t}| < 25 ns, EMCAL-BC=0",netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiEMCALBC0->SetYTitle("#phi (rad)");
    fhEtaPhiEMCALBC0->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiEMCALBC0) ;
    
    fhEtaPhiEMCALBC1  = new TH2F
    ("hEtaPhiEMCALBC1","cluster, #it{E} > 2 GeV, #eta vs #phi, for clusters with 25 < |#it{t}| < 75 ns, EMCAL-BC=1",netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiEMCALBC1->SetYTitle("#phi (rad)");
    fhEtaPhiEMCALBC1->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiEMCALBC1) ;
    
    fhEtaPhiEMCALBCN  = new TH2F
    ("hEtaPhiEMCALBCN","cluster, #it{E} > 2 GeV, #eta vs #phi, for clusters with |#it{t}| > 75 ns, EMCAL-BC>1",netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiEMCALBCN->SetYTitle("#phi (rad)");
    fhEtaPhiEMCALBCN->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiEMCALBCN) ;
    
    for(Int_t i = 0; i < 11; i++)
    {
      fhEtaPhiTriggerEMCALBC[i] = new TH2F
      (Form("hEtaPhiTriggerEMCALBC%d",i-5),
       Form("meson #it{E} > 2 GeV, #eta vs #phi, Trigger EMCAL-BC=%d",i-5),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhEtaPhiTriggerEMCALBC[i]->SetYTitle("#phi (rad)");
      fhEtaPhiTriggerEMCALBC[i]->SetXTitle("#eta");
      outputContainer->Add(fhEtaPhiTriggerEMCALBC[i]) ;
      
      fhTimeTriggerEMCALBC[i] = new TH2F
      (Form("hTimeTriggerEMCALBC%d",i-5),
       Form("meson #it{t} vs #it{E}, Trigger EMCAL-BC=%d",i-5),
       nptbins,ptmin,ptmax, ntimptbins,timemin,timemax);
      fhTimeTriggerEMCALBC[i]->SetXTitle("#it{E} (GeV)");
      fhTimeTriggerEMCALBC[i]->SetYTitle("#it{t} (ns)");
      outputContainer->Add(fhTimeTriggerEMCALBC[i]);
      
      fhTimeTriggerEMCALBCPileUpSPD[i] = new TH2F
      (Form("hTimeTriggerEMCALBC%dPileUpSPD",i-5),
       Form("meson #it{t} vs #it{E}, Trigger EMCAL-BC=%d",i-5),
       nptbins,ptmin,ptmax, ntimptbins,timemin,timemax);
      fhTimeTriggerEMCALBCPileUpSPD[i]->SetXTitle("#it{E} (GeV)");
      fhTimeTriggerEMCALBCPileUpSPD[i]->SetYTitle("#it{t} (ns)");
      outputContainer->Add(fhTimeTriggerEMCALBCPileUpSPD[i]);
      
      fhEtaPhiTriggerEMCALBCUM[i] = new TH2F
      (Form("hEtaPhiTriggerEMCALBC%d_UnMatch",i-5),
       Form("meson #it{E} > 2 GeV, #eta vs #phi, unmatched trigger EMCAL-BC=%d",i-5),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhEtaPhiTriggerEMCALBCUM[i]->SetYTitle("#phi (rad)");
      fhEtaPhiTriggerEMCALBCUM[i]->SetXTitle("#eta");
      outputContainer->Add(fhEtaPhiTriggerEMCALBCUM[i]) ;
      
      fhTimeTriggerEMCALBCUM[i] = new TH2F
      (Form("hTimeTriggerEMCALBC%d_UnMatch",i-5),
       Form("meson #it{t} vs #it{E}, unmatched trigger EMCAL-BC=%d",i-5),
       nptbins,ptmin,ptmax, ntimptbins,timemin,timemax);
      fhTimeTriggerEMCALBCUM[i]->SetXTitle("#it{E} (GeV)");
      fhTimeTriggerEMCALBCUM[i]->SetYTitle("#it{t} (ns)");
      outputContainer->Add(fhTimeTriggerEMCALBCUM[i]);
      
    }
    
    fhTimeTriggerEMCALBC0UMReMatchOpenTime = new TH2F("hTimeTriggerBC0_UnMatch_ReMatch_OpenTime",
                                                      "cluster #it{t} vs #it{E} of clusters, no match, rematch open time",
                                                      nptbins,ptmin,ptmax, ntimptbins,timemin,timemax);
    fhTimeTriggerEMCALBC0UMReMatchOpenTime->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBC0UMReMatchOpenTime->SetYTitle("#it{t} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBC0UMReMatchOpenTime);
    
    
    fhTimeTriggerEMCALBC0UMReMatchCheckNeigh = new TH2F("hTimeTriggerBC0_UnMatch_ReMatch_CheckNeighbours",
                                                        "cluster #it{t} vs #it{E} of clusters, no match, rematch with neigbour parches",
                                                        nptbins,ptmin,ptmax, ntimptbins,timemin,timemax);
    fhTimeTriggerEMCALBC0UMReMatchCheckNeigh->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBC0UMReMatchCheckNeigh->SetYTitle("#it{t} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBC0UMReMatchCheckNeigh);
    
    fhTimeTriggerEMCALBC0UMReMatchBoth = new TH2F("hTimeTriggerBC0_UnMatch_ReMatch_Both",
                                                  "cluster #it{t} vs #it{E} of clusters, no match, rematch open time and neigbour",
                                                  nptbins,ptmin,ptmax, ntimptbins,timemin,timemax);
    fhTimeTriggerEMCALBC0UMReMatchBoth->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBC0UMReMatchBoth->SetYTitle("#it{t} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBC0UMReMatchBoth);
    
  }
  
  if(IsHighMultiplicityAnalysisOn())
  {
    fhPtCentrality  = new TH2F("hPtCentrality","centrality vs #it{p}_{T}",nptbins,ptmin,ptmax, 100,0,100);
    fhPtCentrality->SetYTitle("centrality");
    fhPtCentrality->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtCentrality) ;
    
    fhPtEventPlane  = new TH2F("hPtEventPlane","event plane angle vs #it{p}_{T}",nptbins,ptmin,ptmax, 100,0,TMath::Pi());
    fhPtEventPlane->SetYTitle("Event plane angle (rad)");
    fhPtEventPlane->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtEventPlane) ;
  }
  
  if(fAnaType == kSSCalo)
  {
    fhPtReject  = new TH1F("hPtReject","Number of rejected as #pi^{0} (#eta) decay",nptbins,ptmin,ptmax);
    fhPtReject->SetYTitle("#it{N}");
    fhPtReject->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtReject) ;
    
    fhEReject  = new TH1F("hEReject","Number of rejected as  #pi^{0} (#eta) decay pairs",nptbins,ptmin,ptmax);
    fhEReject->SetYTitle("#it{N}");
    fhEReject->SetXTitle("#it{E} (GeV)");
    outputContainer->Add(fhEReject) ;
    
    fhPtPhiReject  = new TH2F
    ("hPtPhiReject","Rejected #pi^{0} (#eta) cluster: #it{p}_{T} vs #phi",nptbins,ptmin,ptmax, nphibins,phimin,phimax);
    fhPtPhiReject->SetYTitle("#phi (rad)");
    fhPtPhiReject->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtPhiReject) ;
    
    fhPtEtaReject  = new TH2F
    ("hPtEtaReject","Rejected #pi^{0} (#eta) cluster: #it{p}_{T} vs #eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhPtEtaReject->SetYTitle("#eta");
    fhPtEtaReject->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtEtaReject) ;
    
    fhEtaPhiReject  = new TH2F
    ("hEtaPhiReject","Rejected #pi^{0} (#eta) cluster: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
    fhEtaPhiReject->SetYTitle("#phi (rad)");
    fhEtaPhiReject->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiReject) ;
    
    fhNLocMaxPtReject = new TH2F("hNLocMaxPtReject","Number of local maxima in cluster, rejected clusters",
                                 nptbins,ptmin,ptmax,20,0,20);
    fhNLocMaxPtReject ->SetYTitle("N maxima");
    fhNLocMaxPtReject ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhNLocMaxPtReject) ;
  }
  
  fhMass  = new TH2F
  ("hMass","all pairs #it{M}: #it{E} vs #it{M}",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
  fhMass->SetYTitle("#it{M} (GeV/#it{c}^{2})");
  fhMass->SetXTitle("#it{E} (GeV)");
  outputContainer->Add(fhMass) ;
  
  fhSelectedMass  = new TH2F
  ("hSelectedMass","Selected #pi^{0} (#eta) pairs #it{M}: E vs #it{M}",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
  fhSelectedMass->SetYTitle("#it{M} (GeV/#it{c}^{2})");
  fhSelectedMass->SetXTitle("#it{E} (GeV)");
  outputContainer->Add(fhSelectedMass) ;
  
  fhMassPt  = new TH2F
  ("hMassPt","all pairs #it{M}: #it{p}_{T} vs #it{M}",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
  fhMassPt->SetYTitle("#it{M} (GeV/#it{c}^{2})");
  fhMassPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhMassPt) ;
  
  fhMassPtMaxPair  = new TH2F
  ("hMassPtMaxPair","all pairs #it{M}: #it{p}_{T}^{max} vs #it{M}",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
  fhMassPtMaxPair->SetYTitle("#it{M} (GeV/#it{c}^{2})");
  fhMassPtMaxPair->SetXTitle("#it{p}_{T}^{max} (GeV/#it{c})");
  outputContainer->Add(fhMassPtMaxPair) ;
  
  fhMassPtMinPair  = new TH2F
  ("hMassPtMinPair","all pairs #it{M}: #it{p}_{T}^{min} vs #it{M}",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
  fhMassPtMinPair->SetYTitle("#it{M} (GeV/#it{c}^{2})");
  fhMassPtMinPair->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhMassPtMinPair) ;
  
  fhSelectedMassPt  = new TH2F
  ("hSelectedMassPt","Selected #pi^{0} (#eta) pairs #it{M}: #it{p}_{T} vs #it{M}",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
  fhSelectedMassPt->SetYTitle("#it{M} (GeV/#it{c}^{2})");
  fhSelectedMassPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhSelectedMassPt) ;
  
  if(fAnaType != kSSCalo && fSelectPairInIsoCone)
  {
    fhMassPtIsoRCut  = new TH2F
    ("hMassPtIsoRCut",Form("#it{M}: #it{p}_{T} vs #it{M}, for R = %1.1f, #it{p}_{T,1} < %2.2f",fR,fIsoCandMinPt),
     nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
    fhMassPtIsoRCut->SetYTitle("#it{M} (GeV/#it{c}^{2})");
    fhMassPtIsoRCut->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhMassPtIsoRCut) ;
  }
  
  if(fAnaType == kSSCalo)
  {
    fhPtLambda0NoSplitCut  = new TH2F
    ("hPtLambda0NoSplitCut","all clusters: #it{p}_{T} vs #lambda_{0}^{2}",nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
    fhPtLambda0NoSplitCut->SetYTitle("#lambda_{0}^{2}");
    fhPtLambda0NoSplitCut->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtLambda0NoSplitCut) ;
    
    for(Int_t inlm = 0; inlm < 3; inlm++)
    {
      fhMassPtLocMax[inlm]  = new TH2F
      (Form("hMassPtNLocMax%d",inlm+1),Form("all pairs #it{M}: #it{p}_{T} vs #it{M} and NLM=%s",nlm[inlm].Data()),nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
      fhMassPtLocMax[inlm]->SetYTitle("#it{M} (GeV/#it{c}^{2})");
      fhMassPtLocMax[inlm]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMassPtLocMax[inlm]) ;
      
      fhSelectedMassPtLocMax[inlm]  = new TH2F
      (Form("hSelectedMassPtLocMax%d",inlm+1),Form("Selected #pi^{0} (#eta) pairs #it{M}: #it{p}_{T} vs #it{M}, NLM=%s",nlm[inlm].Data()),nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
      fhSelectedMassPtLocMax[inlm]->SetYTitle("#it{M} (GeV/#it{c}^{2})");
      fhSelectedMassPtLocMax[inlm]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhSelectedMassPtLocMax[inlm]) ;
      
      if(fFillAllNLMHistograms)
      {
        for(Int_t iSM = 0; iSM < fNModules; iSM++)
        {
          if(iSM < fFirstModule || iSM > fLastModule) continue;

          fhSelectedMassPtLocMaxSM[inlm][iSM]  = new TH2F
          (Form("hSelectedMassPtLocMax%d_SM%d",inlm+1,iSM),Form("Selected #pi^{0} (#eta) pairs #it{M}: #it{p}_{T} vs #it{M}, NLM=%s for SM=%d",nlm[inlm].Data(),iSM),nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
          fhSelectedMassPtLocMaxSM[inlm][iSM]->SetYTitle("#it{M} (GeV/#it{c}^{2})");
          fhSelectedMassPtLocMaxSM[inlm][iSM]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhSelectedMassPtLocMaxSM[inlm][iSM]) ;
          
          fhSelectedLambda0PtLocMaxSM[inlm][iSM]  = new TH2F
          (Form("hSelectedLambda0PtLocMax%d_SM%d",inlm+1,iSM),Form("Selected #pi^{0} (#eta) pairs #lambda_{0}^{2}: #it{p}_{T} vs #it{M}, NLM=%s for SM=%d",nlm[inlm].Data(),iSM),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhSelectedLambda0PtLocMaxSM[inlm][iSM]->SetYTitle("#lambda_{0}^{2}");
          fhSelectedLambda0PtLocMaxSM[inlm][iSM]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhSelectedLambda0PtLocMaxSM[inlm][iSM]) ;
        }
      }
      
      if(IsDataMC())
      {
        for(Int_t ipart = 0; ipart < fgkNmcTypes; ipart++)
        {
          fhMCSelectedMassPtLocMax[ipart][inlm]  = new TH2F
          (Form("hSelectedMassPtLocMax%d_MC%s",inlm+1,pname[ipart].Data()),
           Form("Selected #pi^{0} (#eta) pairs #it{M}: #it{p}_{T} vs #it{M}, NLM=%s, from MC %s",nlm[inlm].Data(),ptype[ipart].Data()),
           nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
          fhMCSelectedMassPtLocMax[ipart][inlm]->SetYTitle("#it{M} (GeV/#it{c}^{2})");
          fhMCSelectedMassPtLocMax[ipart][inlm]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhMCSelectedMassPtLocMax[ipart][inlm]) ;
        }
      }
    }
    
    if(IsDataMC())
    {
      fhMassNoOverlap  = new TH2F
      ("hMassNoOverlap","all pairs #it{M}: #it{E} vs #it{M}, no overlap",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
      fhMassNoOverlap->SetYTitle("#it{M} (GeV/#it{c}^{2})");
      fhMassNoOverlap->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhMassNoOverlap) ;
      
      fhSelectedMassNoOverlap  = new TH2F
      ("hSelectedMassNoOverlap","Selected #pi^{0} (#eta) pairs #it{M}: #it{E} vs #it{M}, no overlap",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
      fhSelectedMassNoOverlap->SetYTitle("#it{M} (GeV/#it{c}^{2})");
      fhSelectedMassNoOverlap->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhSelectedMassNoOverlap) ;
      
      fhMassPtNoOverlap  = new TH2F
      ("hMassPtNoOverlap","all pairs #it{M}: #it{p}_{T} vs #it{M}, no overlap",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
      fhMassPtNoOverlap->SetYTitle("#it{M} (GeV/#it{c}^{2})");
      fhMassPtNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMassPtNoOverlap) ;
      
      fhSelectedMassPtNoOverlap  = new TH2F
      ("hSelectedMassPtNoOverlap","Selected #pi^{0} (#eta) pairs #it{M}: #it{p}_{T} vs #it{M}, no overlap",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
      fhSelectedMassPtNoOverlap->SetYTitle("#it{M} (GeV/#it{c}^{2})");
      fhSelectedMassPtNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhSelectedMassPtNoOverlap) ;
    }
  }
  
  if(fAnaType != kSSCalo)
  {
    fhPtDecay  = new TH1F("hPtDecay","Selected #pi^{0} (#eta) decay photons",nptbins,ptmin,ptmax);
    fhPtDecay->SetYTitle("#it{N}");
    fhPtDecay->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtDecay) ;
    
    if(IsDataMC())
    {
      fhMCPtDecayLostPairPi0  = new TH1F("hPtDecay_MCPi0DecayLostPair","Selected  #pi^{0} (#eta) decay photons, from MC #gamma #pi^{0} decay, companion lost",
                                         nptbins,ptmin,ptmax);
      fhMCPtDecayLostPairPi0->SetYTitle("#it{N}");
      fhMCPtDecayLostPairPi0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCPtDecayLostPairPi0) ;
      
      fhMCPtDecayLostPairEta  = new TH1F("hPtDecay_MCEtaDecayLostPair","Selected  #pi^{0} (#eta) decay photons, from MC #gamma #eta decay, companion lost",
                                         nptbins,ptmin,ptmax);
      fhMCPtDecayLostPairEta->SetYTitle("#it{N}");
      fhMCPtDecayLostPairEta->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCPtDecayLostPairEta) ;
      
      for(Int_t ipart = 0; ipart < fgkNmcTypes; ipart++)
      {
        fhMCPtDecay[ipart]  = new TH1F(Form("hPtDecay_MC%s",pname[ipart].Data()),
                                       Form("Selected  #pi^{0} (#eta) decay photons, from MC %s",ptype[ipart].Data()),
                                       nptbins,ptmin,ptmax);
        fhMCPtDecay[ipart]->SetYTitle("#it{N}");
        fhMCPtDecay[ipart]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCPtDecay[ipart]) ;
      }
    }
  }
  
  ////////
  
  if( fFillSelectClHisto )
  {
    fhPtLambda0  = new TH2F
    ("hPtLambda0","Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #lambda_{0}",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
    fhPtLambda0->SetYTitle("#lambda_{0}^{2}");
    fhPtLambda0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtLambda0) ;
    
    fhPtLambda1  = new TH2F
    ("hPtLambda1","Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #lambda_{1}",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
    fhPtLambda1->SetYTitle("#lambda_{1}^{2}");
    fhPtLambda1->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtLambda1) ;
    
    if(GetCalorimeter()==kEMCAL &&  GetFirstSMCoveredByTRD() >=0 )
    {
      fhPtLambda0NoTRD  = new TH2F
      ("hPtLambda0NoTRD","Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #lambda_{0}, not behind TRD",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhPtLambda0NoTRD->SetYTitle("#lambda_{0}^{2}");
      fhPtLambda0NoTRD->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtLambda0NoTRD) ;
      
      if(!fFillOnlySimpleSSHisto)
      {
        fhPtFracMaxCellNoTRD  = new TH2F
        ("hPtFracMaxCellNoTRD","Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #lambda_{0}, Max cell fraction of energy, not behind TRD",nptbins,ptmin,ptmax,100,0,1);
        fhPtFracMaxCellNoTRD->SetYTitle("Fraction");
        fhPtFracMaxCellNoTRD->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtFracMaxCellNoTRD) ;
      }
    }
    
    if(!fFillOnlySimpleSSHisto)
    {
      fhPtDispersion  = new TH2F
      ("hPtDispersion","Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs dispersion",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhPtDispersion->SetYTitle("D^{2}");
      fhPtDispersion->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtDispersion) ;
      
      fhPtLambda0FracMaxCellCut  = new TH2F
      ("hPtLambda0FracMaxCellCut","Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #lambda_{0}, Max cell fraction of energy < 0.5",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhPtLambda0FracMaxCellCut->SetYTitle("#lambda_{0}^{2}");
      fhPtLambda0FracMaxCellCut->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtLambda0FracMaxCellCut) ;
      
      fhPtFracMaxCell  = new TH2F
      ("hPtFracMaxCell","Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #lambda_{0}, Max cell fraction of energy",nptbins,ptmin,ptmax,100,0,1);
      fhPtFracMaxCell->SetYTitle("Fraction");
      fhPtFracMaxCell->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtFracMaxCell) ;
      
      fhPtDispEta  = new TH2F ("hPtDispEta","#sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i} - <#eta>)^{2}/ #Sigma w_{i} vs #it{p}_{T}",  nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
      fhPtDispEta->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPtDispEta->SetYTitle("#sigma^{2}_{#eta #eta}");
      outputContainer->Add(fhPtDispEta);
      
      fhPtDispPhi  = new TH2F ("hPtDispPhi","#sigma^{2}_{#phi #phi} = #Sigma w_{i}(#phi_{i} - <#phi>)^{2} / #Sigma w_{i} vs #it{p}_{T}",  nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
      fhPtDispPhi->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPtDispPhi->SetYTitle("#sigma^{2}_{#phi #phi}");
      outputContainer->Add(fhPtDispPhi);
      
      fhPtSumEta  = new TH2F ("hPtSumEta","#sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i})^{2} / #Sigma w_{i} - <#eta>^{2} vs #it{p}_{T}",  nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
      fhPtSumEta->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPtSumEta->SetYTitle("#delta^{2}_{#eta #eta}");
      outputContainer->Add(fhPtSumEta);
      
      fhPtSumPhi  = new TH2F ("hPtSumPhi","#sigma^{2}_{#phi #phi} = #Sigma w_{i}(#phi_{i})^{2}/ #Sigma w_{i} - <#phi>^{2} vs #it{p}_{T}",
                              nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
      fhPtSumPhi->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPtSumPhi->SetYTitle("#delta^{2}_{#phi #phi}");
      outputContainer->Add(fhPtSumPhi);
      
      fhPtSumEtaPhi  = new TH2F ("hPtSumEtaPhi","#delta^{2}_{#eta #phi} = #Sigma w_{i}(#phi_{i} #eta_{i} ) / #Sigma w_{i} - <#phi><#eta> vs #it{p}_{T}",
                                 nptbins,ptmin,ptmax, 2*ssbins,-ssmax,ssmax);
      fhPtSumEtaPhi->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPtSumEtaPhi->SetYTitle("#delta^{2}_{#eta #phi}");
      outputContainer->Add(fhPtSumEtaPhi);
      
      fhPtDispEtaPhiDiff  = new TH2F ("hPtDispEtaPhiDiff","#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta} vs #it{p}_{T}",
                                      nptbins,ptmin,ptmax,200, -10,10);
      fhPtDispEtaPhiDiff->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPtDispEtaPhiDiff->SetYTitle("#sigma^{2}_{#phi #phi}-#sigma^{2}_{#eta #eta}");
      outputContainer->Add(fhPtDispEtaPhiDiff);
      
      fhPtSphericity  = new TH2F ("hPtSphericity","(#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi}) vs #it{p}_{T} (GeV/#it{c})",
                                  nptbins,ptmin,ptmax, 200, -1,1);
      fhPtSphericity->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPtSphericity->SetYTitle("s = (#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi})");
      outputContainer->Add(fhPtSphericity);
      
      for(Int_t i = 0; i < 7; i++)
      {
        fhDispEtaDispPhi[i] = new TH2F (Form("hDispEtaDispPhi_EBin%d",i),Form("#sigma^{2}_{#phi #phi} vs #sigma^{2}_{#eta #eta} for %d < E < %d GeV",bin[i],bin[i+1]),
                                        ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
        fhDispEtaDispPhi[i]->SetXTitle("#sigma^{2}_{#eta #eta}");
        fhDispEtaDispPhi[i]->SetYTitle("#sigma^{2}_{#phi #phi}");
        outputContainer->Add(fhDispEtaDispPhi[i]);
        
        fhLambda0DispEta[i] = new TH2F (Form("hLambda0DispEta_EBin%d",i),Form("#lambda^{2}_{0} vs #sigma^{2}_{#eta #eta} for %d < E < %d GeV",bin[i],bin[i+1]),
                                        ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
        fhLambda0DispEta[i]->SetXTitle("#lambda^{2}_{0}");
        fhLambda0DispEta[i]->SetYTitle("#sigma^{2}_{#eta #eta}");
        outputContainer->Add(fhLambda0DispEta[i]);
        
        fhLambda0DispPhi[i] = new TH2F (Form("hLambda0DispPhi_EBin%d",i),Form("#lambda^{2}_{0}} vs #sigma^{2}_{#phi #phi} for %d < E < %d GeV",bin[i],bin[i+1]),
                                        ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
        fhLambda0DispPhi[i]->SetXTitle("#lambda^{2}_{0}");
        fhLambda0DispPhi[i]->SetYTitle("#sigma^{2}_{#phi #phi}");
        outputContainer->Add(fhLambda0DispPhi[i]);
      }
    }
    
    fhNLocMaxPt = new TH2F("hNLocMaxPt","Number of local maxima in cluster, selected clusters",
                           nptbins,ptmin,ptmax,20,0,20);
    fhNLocMaxPt ->SetYTitle("N maxima");
    fhNLocMaxPt ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhNLocMaxPt) ;
    
    if(fFillAllNLMHistograms)
    {
      for(Int_t iSM = 0; iSM < fNModules; iSM++)
      {
        if(iSM < fFirstModule || iSM > fLastModule) continue;

        fhNLocMaxPtSM[iSM] = new TH2F(Form("hNLocMaxPt_SM%d",iSM),Form("Number of local maxima in cluster, selected clusters in SM %d",iSM),
                                      nptbins,ptmin,ptmax,20,0,20);
        fhNLocMaxPtSM[iSM] ->SetYTitle("N maxima");
        fhNLocMaxPtSM[iSM] ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhNLocMaxPtSM[iSM]) ;
      }
      
      for (Int_t i = 0; i < 3; i++)
      {
        fhPtLambda0LocMax[i]  = new TH2F(Form("hPtLambda0LocMax%d",i+1),
                                         Form("Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #lambda_{0}, NLM=%s",nlm[i].Data()),
                                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhPtLambda0LocMax[i]->SetYTitle("#lambda_{0}^{2}");
        fhPtLambda0LocMax[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtLambda0LocMax[i]) ;
        
        if(IsDataMC())
        {
          for(Int_t ipart = 0; ipart < fgkNmcTypes; ipart++)
          {
            fhMCPtLambda0LocMax[ipart][i]  = new TH2F
            (Form("hPtLambda0LocMax%d_MC%s",i+1,pname[ipart].Data()),
             Form("Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #lambda_{0}, NLM=%s, MC %s",nlm[i].Data(),ptype[ipart].Data()),
             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhMCPtLambda0LocMax[ipart][i]->SetYTitle("#lambda_{0}^{2}");
            fhMCPtLambda0LocMax[ipart][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            outputContainer->Add(fhMCPtLambda0LocMax[ipart][i]) ;
          }
        }
        
        fhPtLambda1LocMax[i]  = new TH2F(Form("hPtLambda1LocMax%d",i+1),
                                         Form("Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #lambda_{1}, %s",nlm[i].Data()),
                                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhPtLambda1LocMax[i]->SetYTitle("#lambda_{1}^{2}");
        fhPtLambda1LocMax[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtLambda1LocMax[i]) ;
        
        if(GetCalorimeter() == kEMCAL && !fFillOnlySimpleSSHisto)
        {
          fhPtDispersionLocMax[i]  = new TH2F(Form("hPtDispersionLocMax%d",i+1),
                                              Form("Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs dispersion^{2}, %s",nlm[i].Data()),
                                              nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtDispersionLocMax[i]->SetYTitle("dispersion^{2}");
          fhPtDispersionLocMax[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtDispersionLocMax[i]) ;
          
          fhPtDispEtaLocMax[i]  = new TH2F(Form("hPtDispEtaLocMax%d",i+1),
                                           Form("Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #sigma_{#eta #eta}, %s",nlm[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtDispEtaLocMax[i]->SetYTitle("#sigma_{#eta #eta}");
          fhPtDispEtaLocMax[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtDispEtaLocMax[i]) ;
          
          fhPtDispPhiLocMax[i]  = new TH2F(Form("hPtDispPhiLocMax%d",i+1),
                                           Form("Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #sigma_{#phi #phi}, %s",nlm[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtDispPhiLocMax[i]->SetYTitle("#sigma_{#phi #phi}");
          fhPtDispPhiLocMax[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtDispPhiLocMax[i]) ;
          
          fhPtSumEtaPhiLocMax[i]  = new TH2F(Form("hPtSumEtaPhiLocMax%d",i+1),
                                             Form("Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #sigma_{#eta #phi}, %s",nlm[i].Data()),
                                             nptbins,ptmin,ptmax,2*ssbins,-ssmax,ssmax);
          fhPtSumEtaPhiLocMax[i]->SetYTitle("#sigma_{#eta #phi}");
          fhPtSumEtaPhiLocMax[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtSumEtaPhiLocMax[i]) ;
          
          fhPtDispEtaPhiDiffLocMax[i]  = new TH2F(Form("hPtDispEtaPhiDiffLocMax%d",i+1),
                                                  Form("Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #sigma_{#phi #phi} - #sigma_{#eta #eta}, %s",nlm[i].Data()),
                                                  nptbins,ptmin,ptmax,200, -10,10);
          fhPtDispEtaPhiDiffLocMax[i]->SetYTitle("#sigma_{#phi #phi} - #sigma_{#eta #eta}");
          fhPtDispEtaPhiDiffLocMax[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtDispEtaPhiDiffLocMax[i]) ;
          
          fhPtSphericityLocMax[i]  = new TH2F(Form("hPtSphericityLocMax%d",i+1),
                                              Form("Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #sigma_{#phi #phi} - #sigma_{#eta #eta} / (#sigma_{#phi #phi} + #sigma_{#eta #eta}), %s",nlm[i].Data()),
                                              nptbins,ptmin,ptmax,200, -1,1);
          fhPtSphericityLocMax[i]->SetYTitle("#sigma_{#phi #phi} - #sigma_{#eta #eta} / (#sigma_{#phi #phi} + #sigma_{#eta #eta})");
          fhPtSphericityLocMax[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtSphericityLocMax[i]) ;
        }
      }
    } // all NLM histos
    
    fhPtNCells  = new TH2F ("hPtNCells","N cells in cluster vs E ", nptbins,ptmin,ptmax, nbins,nmin,nmax);
    fhPtNCells->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhPtNCells->SetYTitle("# of cells in cluster");
    outputContainer->Add(fhPtNCells);
    
    fhPtTime = new TH2F("hPtTime","cluster time vs pair E",nptbins,ptmin,ptmax, tbins,tmin,tmax);
    fhPtTime->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhPtTime->SetYTitle("t (ns)");
    outputContainer->Add(fhPtTime);
    
  }
  
  if(fAnaType != kIMCaloTracks)
  {
    fhEPairDiffTime = new TH2F("hEPairDiffTime","cluster pair time difference vs E",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
    fhEPairDiffTime->SetXTitle("#it{E}_{pair} (GeV)");
    fhEPairDiffTime->SetYTitle("#Delta t (ns)");
    outputContainer->Add(fhEPairDiffTime);
  }
  
  if(fAnaType == kIMCalo)
  {
    TString combiName [] = {"1LocMax","2LocMax","NLocMax","1LocMax2LocMax","1LocMaxNLocMax","2LocMaxNLocMax","1LocMaxSSBad","NLocMaxSSGood"};
    TString combiTitle[] = {"1 Local Maxima in both clusters","2 Local Maxima in both clusters","more than 2 Local Maxima in both clusters",
      "1 Local Maxima paired with 2 Local Maxima","1 Local Maxima paired with more than 2 Local Maxima",
      "2 Local Maxima paired with more than 2 Local Maxima",
      "1 Local Maxima paired with #lambda_{0}^{2}>0.3","N Local Maxima paired with 0.1<#lambda_{0}^{2}<0.3"};
    
    if(fFillAllNLMHistograms)
    {
      for (Int_t i = 0; i < 8 ; i++)
      {
        if (fAnaType == kIMCaloTracks && i > 2 ) continue ;
        
        fhMassPairLocMax[i]  = new TH2F
        (Form("MassPairLocMax%s",combiName[i].Data()),
         Form("#it{M} for decay #gamma pair vs #it{E}_{pair}, origin #pi^{0}, %s", combiTitle[i].Data()),
         nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhMassPairLocMax[i]->SetYTitle("#it{M} (GeV/#it{c}^{2})");
        fhMassPairLocMax[i]->SetXTitle("#it{E}_{pair} (GeV)");
        outputContainer->Add(fhMassPairLocMax[i]) ;
      }
    }
  }
  
  if(fFillTMHisto)
  {
    fhTrackMatchedDEta  = new TH2F
    ("hTrackMatchedDEta",
     "d#eta of cluster-track vs cluster #it{p}_{T}",
     nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
    fhTrackMatchedDEta->SetYTitle("d#eta");
    fhTrackMatchedDEta->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    
    fhTrackMatchedDPhi  = new TH2F
    ("hTrackMatchedDPhi",
     "d#phi of cluster-track vs cluster #it{p}_{T}",
     nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDPhi->SetYTitle("d#phi (rad)");
    fhTrackMatchedDPhi->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    
    fhTrackMatchedDEtaDPhi  = new TH2F
    ("hTrackMatchedDEtaDPhi",
     "d#eta vs d#phi of cluster-track",
     nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDEtaDPhi->SetYTitle("d#phi (rad)");
    fhTrackMatchedDEtaDPhi->SetXTitle("d#eta");
    
    outputContainer->Add(fhTrackMatchedDEta) ;
    outputContainer->Add(fhTrackMatchedDPhi) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhi) ;
    
    fhTrackMatchedDEtaPos  = new TH2F
    ("hTrackMatchedDEtaPos",
     "d#eta of cluster-track vs cluster #it{p}_{T}",
     nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
    fhTrackMatchedDEtaPos->SetYTitle("d#eta");
    fhTrackMatchedDEtaPos->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    
    fhTrackMatchedDPhiPos  = new TH2F
    ("hTrackMatchedDPhiPos",
     "d#phi of cluster-track vs cluster #it{p}_{T}",
     nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDPhiPos->SetYTitle("d#phi (rad)");
    fhTrackMatchedDPhiPos->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    
    fhTrackMatchedDEtaDPhiPos  = new TH2F
    ("hTrackMatchedDEtaDPhiPos",
     "d#eta vs d#phi of cluster-track",
     nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDEtaDPhiPos->SetYTitle("d#phi (rad)");
    fhTrackMatchedDEtaDPhiPos->SetXTitle("d#eta");
    
    outputContainer->Add(fhTrackMatchedDEtaPos) ;
    outputContainer->Add(fhTrackMatchedDPhiPos) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhiPos) ;
    
    fhTrackMatchedDEtaNeg  = new TH2F
    ("hTrackMatchedDEtaNeg",
     "d#eta of cluster-track vs cluster #it{p}_{T}",
     nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
    fhTrackMatchedDEtaNeg->SetYTitle("d#eta");
    fhTrackMatchedDEtaNeg->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    
    fhTrackMatchedDPhiNeg  = new TH2F
    ("hTrackMatchedDPhiNeg",
     "d#phi of cluster-track vs cluster #it{p}_{T}",
     nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDPhiNeg->SetYTitle("d#phi (rad)");
    fhTrackMatchedDPhiNeg->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    
    fhTrackMatchedDEtaDPhiNeg  = new TH2F
    ("hTrackMatchedDEtaDPhiNeg",
     "d#eta vs d#phi of cluster-track",
     nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDEtaDPhiNeg->SetYTitle("d#phi (rad)");
    fhTrackMatchedDEtaDPhiNeg->SetXTitle("d#eta");
    
    outputContainer->Add(fhTrackMatchedDEtaNeg) ;
    outputContainer->Add(fhTrackMatchedDPhiNeg) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhiNeg) ;
    
    fhdEdx  = new TH2F ("hdEdx","matched track <dE/dx> vs cluster #it{p}_{T}", nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
    fhdEdx->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhdEdx->SetYTitle("<#it{dE}/#it{dx}>");
    outputContainer->Add(fhdEdx);
    
    fhEOverP  = new TH2F ("hEOverP","matched track E/p vs cluster #it{p}_{T}", nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
    fhEOverP->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhEOverP->SetYTitle("#it{E}/#it{p}");
    outputContainer->Add(fhEOverP);
    
    if(GetCalorimeter()==kEMCAL &&  GetFirstSMCoveredByTRD() >=0)
    {
      fhEOverPNoTRD  = new TH2F ("hEOverPNoTRD","matched track E/p vs cluster E, SM not behind TRD ", nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
      fhEOverPNoTRD->SetXTitle("#it{E} (GeV)");
      fhEOverPNoTRD->SetYTitle("#it{E}/#it{p}");
      outputContainer->Add(fhEOverPNoTRD);
    }
    
    if(IsDataMC() && fFillTMHisto)
    {
      fhTrackMatchedMCParticlePt  = new TH2F
      ("hTrackMatchedMCParticlePt",
       "Origin of particle vs energy",
       nptbins,ptmin,ptmax,8,0,8);
      fhTrackMatchedMCParticlePt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      //fhTrackMatchedMCParticlePt->SetYTitle("Particle type");
      
      fhTrackMatchedMCParticlePt->GetYaxis()->SetBinLabel(1 ,"Photon");
      fhTrackMatchedMCParticlePt->GetYaxis()->SetBinLabel(2 ,"Electron");
      fhTrackMatchedMCParticlePt->GetYaxis()->SetBinLabel(3 ,"Meson Merged");
      fhTrackMatchedMCParticlePt->GetYaxis()->SetBinLabel(4 ,"Rest");
      fhTrackMatchedMCParticlePt->GetYaxis()->SetBinLabel(5 ,"Conv. Photon");
      fhTrackMatchedMCParticlePt->GetYaxis()->SetBinLabel(6 ,"Conv. Electron");
      fhTrackMatchedMCParticlePt->GetYaxis()->SetBinLabel(7 ,"Conv. Merged");
      fhTrackMatchedMCParticlePt->GetYaxis()->SetBinLabel(8 ,"Conv. Rest");
      
      outputContainer->Add(fhTrackMatchedMCParticlePt);
      
      fhTrackMatchedMCParticleDEta  = new TH2F
      ("hTrackMatchedMCParticleDEta",
       "Origin of particle vs #eta residual",
       nresetabins,resetamin,resetamax,8,0,8);
      fhTrackMatchedMCParticleDEta->SetXTitle("#Delta #eta");
      //fhTrackMatchedMCParticleDEta->SetYTitle("Particle type");
      
      fhTrackMatchedMCParticleDEta->GetYaxis()->SetBinLabel(1 ,"Photon");
      fhTrackMatchedMCParticleDEta->GetYaxis()->SetBinLabel(2 ,"Electron");
      fhTrackMatchedMCParticleDEta->GetYaxis()->SetBinLabel(3 ,"Meson Merged");
      fhTrackMatchedMCParticleDEta->GetYaxis()->SetBinLabel(4 ,"Rest");
      fhTrackMatchedMCParticleDEta->GetYaxis()->SetBinLabel(5 ,"Conv. Photon");
      fhTrackMatchedMCParticleDEta->GetYaxis()->SetBinLabel(6 ,"Conv. Electron");
      fhTrackMatchedMCParticleDEta->GetYaxis()->SetBinLabel(7 ,"Conv. Merged");
      fhTrackMatchedMCParticleDEta->GetYaxis()->SetBinLabel(8 ,"Conv. Rest");
      
      outputContainer->Add(fhTrackMatchedMCParticleDEta);
      
      fhTrackMatchedMCParticleDPhi  = new TH2F
      ("hTrackMatchedMCParticleDPhi",
       "Origin of particle vs #phi residual",
       nresphibins,resphimin,resphimax,8,0,8);
      fhTrackMatchedMCParticleDPhi->SetXTitle("#Delta #phi");
      //fhTrackMatchedMCParticleDPhi->SetYTitle("Particle type");
      
      fhTrackMatchedMCParticleDPhi->GetYaxis()->SetBinLabel(1 ,"Photon");
      fhTrackMatchedMCParticleDPhi->GetYaxis()->SetBinLabel(2 ,"Electron");
      fhTrackMatchedMCParticleDPhi->GetYaxis()->SetBinLabel(3 ,"Meson Merged");
      fhTrackMatchedMCParticleDPhi->GetYaxis()->SetBinLabel(4 ,"Rest");
      fhTrackMatchedMCParticleDPhi->GetYaxis()->SetBinLabel(5 ,"Conv. Photon");
      fhTrackMatchedMCParticleDPhi->GetYaxis()->SetBinLabel(6 ,"Conv. Electron");
      fhTrackMatchedMCParticleDPhi->GetYaxis()->SetBinLabel(7 ,"Conv. Merged");
      fhTrackMatchedMCParticleDPhi->GetYaxis()->SetBinLabel(8 ,"Conv. Rest");
      
      outputContainer->Add(fhTrackMatchedMCParticleDPhi);
    }
  }
  
  if(fFillWeightHistograms)
  {
    fhECellClusterRatio  = new TH2F ("hECellClusterRatio"," cell energy / cluster energy vs cluster energy, for selected decay photons from neutral meson",
                                     nptbins,ptmin,ptmax, 100,0,1.);
    fhECellClusterRatio->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhECellClusterRatio->SetYTitle("#it{E}_{cell i}/#it{E}_{cluster}");
    outputContainer->Add(fhECellClusterRatio);
    
    fhECellClusterLogRatio  = new TH2F ("hECellClusterLogRatio"," Log(cell energy / cluster energy) vs cluster energy, for selected decay photons from neutral meson",
                                        nptbins,ptmin,ptmax, 100,-10,0);
    fhECellClusterLogRatio->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhECellClusterLogRatio->SetYTitle("Log (#it{E}_{max cell}/#it{E}_{cluster})");
    outputContainer->Add(fhECellClusterLogRatio);
    
    fhEMaxCellClusterRatio  = new TH2F ("hEMaxCellClusterRatio"," max cell energy / cluster energy vs cluster energy, for selected decay photons from neutral meson",
                                        nptbins,ptmin,ptmax, 100,0,1.);
    fhEMaxCellClusterRatio->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhEMaxCellClusterRatio->SetYTitle("#it{E}_{max cell}/#it{E}_{cluster}");
    outputContainer->Add(fhEMaxCellClusterRatio);
    
    fhEMaxCellClusterLogRatio  = new TH2F ("hEMaxCellClusterLogRatio"," Log(max cell energy / cluster energy) vs cluster energy, for selected decay photons from neutral meson",
                                           nptbins,ptmin,ptmax, 100,-10,0);
    fhEMaxCellClusterLogRatio->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhEMaxCellClusterLogRatio->SetYTitle("Log (#it{E}_{max cell}/#it{E}_{cluster})");
    outputContainer->Add(fhEMaxCellClusterLogRatio);
    
    for(Int_t iw = 0; iw < 14; iw++)
    {
      fhLambda0ForW0[iw]  = new TH2F (Form("hLambda0ForW0%d",iw),Form("shower shape, #lambda^{2}_{0} vs E, w0 = %1.1f, for selected decay photons from neutral meson",1+0.5*iw),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhLambda0ForW0[iw]->SetXTitle("#it{E}_{cluster}");
      fhLambda0ForW0[iw]->SetYTitle("#lambda^{2}_{0}");
      outputContainer->Add(fhLambda0ForW0[iw]);
      
      //      fhLambda1ForW0[iw]  = new TH2F (Form("hLambda1ForW0%d",iw),Form("shower shape, #lambda^{2}_{1} vs E, w0 = %1.1f, for selected decay photons from neutral meson",0.5+0.5*iw),
      //                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      //      fhLambda1ForW0[iw]->SetXTitle("#it{E}_{cluster}");
      //      fhLambda1ForW0[iw]->SetYTitle("#lambda^{2}_{1}");
      //      outputContainer->Add(fhLambda1ForW0[iw]);
    }
  }
  
  if(IsDataMC())
  {
    // Origin
    
    fhMCPi0PtOrigin     = new TH2F("hMCPi0PtOrigin","Reconstructed pair from generated #pi^{0} #it{p}_{T} vs origin",nptbins,ptmin,ptmax,11,0,11) ;
    fhMCPi0PtOrigin->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhMCPi0PtOrigin->SetYTitle("Origin");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(1 ,"Status 21");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(2 ,"Quark");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(3 ,"qq Resonances");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(4 ,"Resonances");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(5 ,"#rho");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(6 ,"#omega");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(7 ,"K");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(8 ,"Other");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(9 ,"#eta");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(10 ,"#eta prime");
    outputContainer->Add(fhMCPi0PtOrigin) ;
    
    fhMCNotResonancePi0PtOrigin     = new TH2F("hMCNotResonancePi0PtOrigin","Reconstructed pair from generated #pi^{0} #it{p}_{T} vs origin",nptbins,ptmin,ptmax,11,0,11) ;
    fhMCNotResonancePi0PtOrigin->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhMCNotResonancePi0PtOrigin->SetYTitle("Origin");
    fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(1 ,"Status 21");
    fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(2 ,"Quark");
    fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(3 ,"qq Resonances");
    fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(4 ,"Resonances");
    fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(5 ,"#rho");
    fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(6 ,"#omega");
    fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(7 ,"K");
    fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(8 ,"Other");
    fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(9 ,"#eta");
    fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(10 ,"#eta prime");
    outputContainer->Add(fhMCNotResonancePi0PtOrigin) ;
    
    fhMCPi0PtStatus     = new TH2F("hMCPi0PtStatus","Reconstructed pair from generated #pi^{0} #it{p}_{T} vs status",nptbins,ptmin,ptmax,101,-50,50) ;
    fhMCPi0PtStatus->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhMCPi0PtStatus->SetYTitle("Status");
    outputContainer->Add(fhMCPi0PtStatus) ;
    
    fhMCEtaPtOrigin     = new TH2F("hMCEtaPtOrigin","Reconstructed pair from generated #pi^{0} #it{p}_{T} vs origin",nptbins,ptmin,ptmax,7,0,7) ;
    fhMCEtaPtOrigin->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhMCEtaPtOrigin->SetYTitle("Origin");
    fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(1 ,"Status 21");
    fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(2 ,"Quark");
    fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(3 ,"qq Resonances");
    fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(4 ,"Resonances");
    fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(5 ,"Other");
    fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(6 ,"#eta prime");
    outputContainer->Add(fhMCEtaPtOrigin) ;
    
    fhMCPi0ProdVertex = new TH2F("hMCPi0ProdVertex","Selected reco pair from generated #pi^{0} #it{p}_{T} vs production vertex",200,ptmin,20+ptmin,5000,0,500) ;
    fhMCPi0ProdVertex->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhMCPi0ProdVertex->SetYTitle("#it{R} (cm)");
    outputContainer->Add(fhMCPi0ProdVertex) ;
    
    fhMCEtaProdVertex = new TH2F("hMCEtaProdVertex","Selected reco pair from generated #eta #it{p}_{T} vs production vertex",200,ptmin,20+ptmin,5000,0,500) ;
    fhMCEtaProdVertex->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhMCEtaProdVertex->SetYTitle("#it{R} (cm)");
    outputContainer->Add(fhMCEtaProdVertex) ;
    
    if(GetReader()->GetDataType() != AliCaloTrackReader::kMC && fAnaType==kSSCalo)
    {
      fhMCPi0PtGenRecoFraction = new TH2F("hMCPi0PtGenRecoFraction","Number of clusters from #pi^{0} (2 #gamma) identified as #pi^{0} (#eta), #it{p}_{T} versus E primary #pi^{0} / E reco",
                                          nptbins,ptmin,ptmax,200,0,2);
      fhMCPi0PtGenRecoFraction->SetXTitle("#it{p}^{rec}_{T} (GeV/#it{c})");
      fhMCPi0PtGenRecoFraction->SetYTitle("#it{E}^{#pi^{0} mother} / #it{E}^{rec}");
      outputContainer->Add(fhMCPi0PtGenRecoFraction) ;
      
      fhMCEtaPtGenRecoFraction = new TH2F("hMCEtaPtGenRecoFraction","Number of clusters from #eta (2 #gamma) identified as #pi^{0} (#eta),#it{p}_{T} versus E primary #eta / E reco",
                                          nptbins,ptmin,ptmax,200,0,2);
      fhMCEtaPtGenRecoFraction->SetXTitle("#it{p}^{rec}_{T} (GeV/#it{c})");
      fhMCEtaPtGenRecoFraction->SetYTitle("#it{E}^{ #eta mother} / #it{E}^{rec}");
      outputContainer->Add(fhMCEtaPtGenRecoFraction) ;
      
      fhMCPi0DecayPt = new TH1F("hMCPi0DecayPt","Number of #gamma from #pi^{0} decay  identified as #pi^{0} (#eta)",nptbins,ptmin,ptmax);
      fhMCPi0DecayPt->SetYTitle("#it{N}");
      fhMCPi0DecayPt->SetXTitle("#it{p}^{rec}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCPi0DecayPt) ;
      
      fhMCPi0DecayPtFraction = new TH2F("hMCPi0DecayPtFraction","Number of #gamma from #pi^{0} decay  identified as #pi^{0} (#eta), #it{p}_{T} versus E primary  #gamma / #it{E} primary #pi^{0}",
                                        nptbins,ptmin,ptmax,100,0,1);
      fhMCPi0DecayPtFraction->SetXTitle("p^{rec}_{T} (GeV/#it{c})");
      fhMCPi0DecayPtFraction->SetYTitle("E^{gen} / E^{gen-mother}");
      outputContainer->Add(fhMCPi0DecayPtFraction) ;
      
      fhMCEtaDecayPt = new TH1F("hMCEtaDecayPt","Number of #gamma from #eta decay  identified as #pi^{0} (#eta)",nptbins,ptmin,ptmax);
      fhMCEtaDecayPt->SetYTitle("#it{N}");
      fhMCEtaDecayPt->SetXTitle("#it{p}^{rec}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCEtaDecayPt) ;
      
      fhMCEtaDecayPtFraction = new TH2F("hMCEtaDecayPtFraction","Number of #gamma from #eta decay  identified as #pi^{0} (#eta), #it{p}_{T} versus E primary  #gamma / E primary #eta",
                                        nptbins,ptmin,ptmax,100,0,1);
      fhMCEtaDecayPtFraction->SetXTitle("#it{p}^{rec}_{T} (GeV/#it{c})");
      fhMCEtaDecayPtFraction->SetYTitle("#it{E}^{gen} / #it{E}^{gen-mother}");
      outputContainer->Add(fhMCEtaDecayPtFraction) ;
      
      fhMCOtherDecayPt = new TH1F("hMCOtherDecayPt","Number of #gamma decay (not #eta or #pi^{0})  identified as #pi^{0} (#eta)",nptbins,ptmin,ptmax);
      fhMCOtherDecayPt->SetYTitle("#it{N}");
      fhMCOtherDecayPt->SetXTitle("#it{p}^{rec}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCOtherDecayPt) ;
    }
    
    if(fAnaType!=kSSCalo)
    {
      fhAnglePairMCPi0  = new TH2F
      ("AnglePairMCPi0",
       "Angle between decay #gamma pair vs #it{E}_{pair}, origin #pi^{0}",nptbins,ptmin,ptmax,250,0,0.5);
      fhAnglePairMCPi0->SetYTitle("#alpha (rad)");
      fhAnglePairMCPi0->SetXTitle("#it{E}_{pair} (GeV)");
      outputContainer->Add(fhAnglePairMCPi0) ;
      
      fhAnglePairMCEta  = new TH2F
      ("AnglePairMCEta",
       "Angle between decay #gamma pair vs #it{E}_{pair}, origin #eta",nptbins,ptmin,ptmax,250,0,0.5);
      fhAnglePairMCEta->SetYTitle("#alpha (rad)");
      fhAnglePairMCEta->SetXTitle("#it{E}_{pair} (GeV)");
      outputContainer->Add(fhAnglePairMCEta) ;
      
      fhMassPairMCPi0  = new TH2F
      ("MassPairMCPi0",
       "#it{M} for decay #gamma pair vs #it{E}_{pair}, origin #pi^{0}",nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
      fhMassPairMCPi0->SetYTitle("#it{M} (GeV/#it{c}^{2})");
      fhMassPairMCPi0->SetXTitle("#it{E}_{pair} (GeV)");
      outputContainer->Add(fhMassPairMCPi0) ;
      
      fhMassPairMCEta  = new TH2F
      ("MassPairMCEta",
       "#it{M} for decay #gamma pair vs #it{E}_{pair}, origin #eta",nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
      fhMassPairMCEta->SetYTitle("#it{M} (GeV/#it{c}^{2})");
      fhMassPairMCEta->SetXTitle("#it{E}_{pair} (GeV)");
      outputContainer->Add(fhMassPairMCEta) ;
    }
    
    Int_t ntypes = fgkNmcTypes;
    if(fAnaType != kSSCalo) ntypes = 2;
    
    for(Int_t i = 0; i < ntypes; i++)
    {
      fhMCE[i]  = new TH1F
      (Form("hE_MC%s",pname[i].Data()),
       Form("Identified as #pi^{0} (#eta), cluster from %s",
            ptype[i].Data()),
       nptbins,ptmin,ptmax);
      fhMCE[i]->SetYTitle("#it{N}");
      fhMCE[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhMCE[i]) ;
      
      fhMCPt[i]  = new TH1F
      (Form("hPt_MC%s",pname[i].Data()),
       Form("Identified as #pi^{0} (#eta), cluster from %s",
            ptype[i].Data()),
       nptbins,ptmin,ptmax);
      fhMCPt[i]->SetYTitle("#it{N}");
      fhMCPt[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCPt[i]) ;
      
      if(IsHighMultiplicityAnalysisOn())
      {
        fhMCPtCentrality[i]  = new TH2F
        (Form("hPtCentrality_MC%s",pname[i].Data()),
         Form("Identified as #pi^{0} (#eta), cluster from %s",
              ptype[i].Data()),
         nptbins,ptmin,ptmax, 100,0,100);
        fhMCPtCentrality[i]->SetYTitle("centrality");
        fhMCPtCentrality[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCPtCentrality[i]) ;
      }
      
      if(fAnaType == kSSCalo)
      {
        fhMCNLocMaxPt[i] = new TH2F
        (Form("hNLocMaxPt_MC%s",pname[i].Data()),
         Form("cluster from %s, #it{p}_{T} of cluster vs NLM, accepted",ptype[i].Data()),
         nptbins,ptmin,ptmax,20,0,20);
        fhMCNLocMaxPt[i] ->SetYTitle("#it{NLM}");
        fhMCNLocMaxPt[i] ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCNLocMaxPt[i]) ;
        
        fhMCNLocMaxPtReject[i] = new TH2F
        (Form("hNLocMaxPtReject_MC%s",pname[i].Data()),
         Form("cluster from %s, #it{p}_{T} of cluster vs NLM, rejected",ptype[i].Data()),
         nptbins,ptmin,ptmax,20,0,20);
        fhMCNLocMaxPtReject[i] ->SetYTitle("#it{NLM}");
        fhMCNLocMaxPtReject[i] ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCNLocMaxPtReject[i]) ;
        
        fhMCEReject[i]  = new TH1F
        (Form("hEReject_MC%s",pname[i].Data()),
         Form("Rejected as #pi^{0} (#eta), cluster from %s",
              ptype[i].Data()),
         nptbins,ptmin,ptmax);
        fhMCEReject[i]->SetYTitle("#it{N}");
        fhMCEReject[i]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCEReject[i]) ;
        
        fhMCPtReject[i]  = new TH1F
        (Form("hPtReject_MC%s",pname[i].Data()),
         Form("Rejected as #pi^{0} (#eta), cluster from %s",
              ptype[i].Data()),
         nptbins,ptmin,ptmax);
        fhMCPtReject[i]->SetYTitle("#it{N}");
        fhMCPtReject[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCPtReject[i]) ;
      }
      
      fhMCPtPhi[i]  = new TH2F
      (Form("hPtPhi_MC%s",pname[i].Data()),
       Form("Identified as #pi^{0} (#eta), cluster from %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nphibins,phimin,phimax);
      fhMCPtPhi[i]->SetYTitle("#phi");
      fhMCPtPhi[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCPtPhi[i]) ;
      
      fhMCPtEta[i]  = new TH2F
      (Form("hPtEta_MC%s",pname[i].Data()),
       Form("Identified as #pi^{0} (#eta), cluster from %s",
            ptype[i].Data()),nptbins,ptmin,ptmax,netabins,etamin,etamax);
      fhMCPtEta[i]->SetYTitle("#eta");
      fhMCPtEta[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCPtEta[i]) ;
      
      fhMCMassPt[i]  = new TH2F
      (Form("hMassPt_MC%s",pname[i].Data()),
       Form("all pairs #it{M}: #it{p}_{T} vs #it{M} from %s",ptype[i].Data()),
       nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
      fhMCMassPt[i]->SetYTitle("#it{M} (GeV/#it{c}^{2})");
      fhMCMassPt[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCMassPt[i]) ;
      
      fhMCSelectedMassPt[i]  = new TH2F
      (Form("hSelectedMassPt_MC%s",pname[i].Data()),
       Form("Selected #pi^{0} (#eta) pairs #it{M}: #it{p}_{T} vs #it{M} from %s",ptype[i].Data()),
       nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
      fhMCSelectedMassPt[i]->SetYTitle("#it{M} (GeV/#it{c}^{2})");
      fhMCSelectedMassPt[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCSelectedMassPt[i]) ;
      
      if(fAnaType == kSSCalo)
      {
        fhMCMassPtNoOverlap[i]  = new TH2F
        (Form("hMassPtNoOverlap_MC%s",pname[i].Data()),
         Form("all pairs #it{M}: #it{p}_{T} vs #it{M} from %s, no overlap",ptype[i].Data()),
         nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
        fhMCMassPt[i]->SetYTitle("#it{M} (GeV/#it{c}^{2})");
        fhMCMassPt[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCMassPtNoOverlap[i]) ;
        
        fhMCSelectedMassPtNoOverlap[i]  = new TH2F
        (Form("hSelectedMassPtNoOverlap_MC%s",pname[i].Data()),
         Form("Selected #pi^{0} (#eta) pairs #it{M}: #it{p}_{T} vs #it{M} from %s, no overlap",ptype[i].Data()),
         nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
        fhMCSelectedMassPtNoOverlap[i]->SetYTitle("#it{M} (GeV/#it{c}^{2})");
        fhMCSelectedMassPtNoOverlap[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCSelectedMassPtNoOverlap[i]) ;
      }
      
      if( fFillSelectClHisto )
      {
        fhMCPtLambda0[i]  = new TH2F(Form("hELambda0_MC%s",pname[i].Data()),
                                     Form("Selected pair, cluster from %s : #it{p}_{T} vs #lambda_{0}^{2}",ptype[i].Data()),
                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhMCPtLambda0[i]->SetYTitle("#lambda_{0}^{2}");
        fhMCPtLambda0[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCPtLambda0[i]) ;
        
        fhMCPtLambda1[i]  = new TH2F(Form("hELambda1_MC%s",pname[i].Data()),
                                     Form("Selected pair, cluster from %s : #it{p}_{T} vs #lambda_{1}^{2}",ptype[i].Data()),
                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhMCPtLambda1[i]->SetYTitle("#lambda_{1}^{2}");
        fhMCPtLambda1[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCPtLambda1[i]) ;
        
        if(GetCalorimeter()==kEMCAL &&  GetFirstSMCoveredByTRD() >= 0)
        {
          fhMCPtLambda0NoTRD[i]  = new TH2F(Form("hELambda0NoTRD_MC%s",pname[i].Data()),
                                            Form("Selected pair, cluster from %s : #it{p}_{T} vs #lambda_{0}^{2}, NoTRD",ptype[i].Data()),
                                            nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhMCPtLambda0NoTRD[i]->SetYTitle("#lambda_{0}^{2}");
          fhMCPtLambda0NoTRD[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhMCPtLambda0NoTRD[i]) ;
        }
        
        if(!fFillOnlySimpleSSHisto)
        {
          fhMCPtDispersion[i]  = new TH2F(Form("hEDispersion_MC%s",pname[i].Data()),
                                          Form("Selected pair, cluster from %s : #it{p}_{T} vs dispersion^{2}",ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhMCPtDispersion[i]->SetYTitle("#it{D}^{2}");
          fhMCPtDispersion[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhMCPtDispersion[i]) ;
          
          fhMCPtDispEta[i]  = new TH2F (Form("hPtDispEta_MC%s",pname[i].Data()),
                                        Form("cluster from %s : #sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i} - <#eta>)^{2}/ #Sigma w_{i} vs #it{p}_{T}",ptype[i].Data()),
                                        nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
          fhMCPtDispEta[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhMCPtDispEta[i]->SetYTitle("#sigma^{2}_{#eta #eta}");
          outputContainer->Add(fhMCPtDispEta[i]);
          
          fhMCPtDispPhi[i]  = new TH2F (Form("hPtDispPhi_MC%s",pname[i].Data()),
                                        Form("cluster from %s : #sigma^{2}_{#phi #phi} = #Sigma w_{i}(#phi_{i} - <#phi>)^{2} / #Sigma w_{i} vs #it{p}_{T}",ptype[i].Data()),
                                        nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
          fhMCPtDispPhi[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhMCPtDispPhi[i]->SetYTitle("#sigma^{2}_{#phi #phi}");
          outputContainer->Add(fhMCPtDispPhi[i]);
          
          fhMCPtSumEtaPhi[i]  = new TH2F (Form("hPtSumEtaPhi_MC%s",pname[i].Data()),
                                          Form("cluster from %s : #delta^{2}_{#eta #phi} = #Sigma w_{i}(#phi_{i} #eta_{i} ) / #Sigma w_{i} - <#phi><#eta> vs #it{p}_{T}",ptype[i].Data()),
                                          nptbins,ptmin,ptmax, 2*ssbins,-ssmax,ssmax);
          fhMCPtSumEtaPhi[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhMCPtSumEtaPhi[i]->SetYTitle("#delta^{2}_{#eta #phi}");
          outputContainer->Add(fhMCPtSumEtaPhi[i]);
          
          fhMCPtDispEtaPhiDiff[i]  = new TH2F (Form("hPtDispEtaPhiDiff_MC%s",pname[i].Data()),
                                               Form("cluster from %s : #sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta} vs #it{p}_{T}",ptype[i].Data()),
                                               nptbins,ptmin,ptmax,200,-10,10);
          fhMCPtDispEtaPhiDiff[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhMCPtDispEtaPhiDiff[i]->SetYTitle("#sigma^{2}_{#phi #phi}-#sigma^{2}_{#eta #eta}");
          outputContainer->Add(fhMCPtDispEtaPhiDiff[i]);
          
          fhMCPtSphericity[i]  = new TH2F (Form("hPtSphericity_MC%s",pname[i].Data()),
                                           Form("cluster from %s : (#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi}) vs E",ptype[i].Data()),
                                           nptbins,ptmin,ptmax, 200,-1,1);
          fhMCPtSphericity[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhMCPtSphericity[i]->SetYTitle("#it{s} = (#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi})");
          outputContainer->Add(fhMCPtSphericity[i]);
          
          for(Int_t ie = 0; ie < 7; ie++)
          {
            fhMCDispEtaDispPhi[ie][i] = new TH2F (Form("hMCDispEtaDispPhi_EBin%d_MC%s",ie,pname[i].Data()),
                                                  Form("cluster from %s : #sigma^{2}_{#phi #phi} vs #sigma^{2}_{#eta #eta} for %d < E < %d GeV",ptype[i].Data(),bin[ie],bin[ie+1]),
                                                  ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
            fhMCDispEtaDispPhi[ie][i]->SetXTitle("#sigma^{2}_{#eta #eta}");
            fhMCDispEtaDispPhi[ie][i]->SetYTitle("#sigma^{2}_{#phi #phi}");
            outputContainer->Add(fhMCDispEtaDispPhi[ie][i]);
            
            fhMCLambda0DispEta[ie][i] = new TH2F (Form("hMCLambda0DispEta_EBin%d_MC%s",ie,pname[i].Data()),
                                                  Form("cluster from %s : #lambda^{2}_{0} vs #sigma^{2}_{#eta #eta} for %d < E < %d GeV",ptype[i].Data(),bin[ie],bin[ie+1]),
                                                  ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
            fhMCLambda0DispEta[ie][i]->SetXTitle("#lambda^{2}_{0}");
            fhMCLambda0DispEta[ie][i]->SetYTitle("#sigma^{2}_{#phi #phi}");
            outputContainer->Add(fhMCLambda0DispEta[ie][i]);
            
            fhMCLambda0DispPhi[ie][i] = new TH2F (Form("hMCLambda0DispPhi_EBin%d_MC%s",ie,pname[i].Data()),
                                                  Form("cluster from %s :#lambda^{2}_{0} vs #sigma^{2}_{#phi #phi} for %d < E < %d GeV",ptype[i].Data(),bin[ie],bin[ie+1]),
                                                  ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
            fhMCLambda0DispPhi[ie][i]->SetXTitle("#lambda^{2}_{0}");
            fhMCLambda0DispPhi[ie][i]->SetYTitle("#sigma^{2}_{#phi #phi}");
            outputContainer->Add(fhMCLambda0DispPhi[ie][i]);
          }
          
          fhMCPtLambda0FracMaxCellCut[i]  = new TH2F(Form("hELambda0FracMaxCellCut_MC%s",pname[i].Data()),
                                                     Form("Selected pair, cluster from %s : #it{p}_{T} vs #lambda_{0}^{2}, Max cell fraction of energy < 0.5 ",ptype[i].Data()),
                                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhMCPtLambda0FracMaxCellCut[i]->SetYTitle("#lambda_{0}^{2}");
          fhMCPtLambda0FracMaxCellCut[i]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhMCPtLambda0FracMaxCellCut[i]) ;
          
          fhMCPtFracMaxCell[i]  = new TH2F(Form("hEFracMaxCell_MC%s",pname[i].Data()),
                                           Form("Selected pair, cluster from %s : #it{p}_{T} vs Max cell fraction of energy",ptype[i].Data()),
                                           nptbins,ptmin,ptmax,100,0,1);
          fhMCPtFracMaxCell[i]->SetYTitle("#it{Fraction}");
          fhMCPtFracMaxCell[i]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhMCPtFracMaxCell[i]) ;
        }
      }
    }// MC particle loop
  }//Histos with MC
  
  if(fAnaType==kSSCalo)
  {
    fhAsymmetry  = new TH2F ("hAsymmetry","#it{A} = ( #it{E}_{1} - #it{E}_{2} ) / ( #it{E}_{1} + #it{E}_{2} ) vs #it{E}",
                             nptbins,ptmin,ptmax, 200, -1,1);
    fhAsymmetry->SetXTitle("#it{E} (GeV)");
    fhAsymmetry->SetYTitle("#it{A} = ( #it{E}_{1} - #it{E}_{2} ) / ( #it{E}_{1} + #it{E}_{2} )");
    outputContainer->Add(fhAsymmetry);
    
    fhSelectedAsymmetry  = new TH2F ("hSelectedAsymmetry","#it{A} = ( #it{E}_{1} - #it{E}_{2} ) / ( #it{E}_{1} + #it{E}_{2} ) vs #it{E}",
                                     nptbins,ptmin,ptmax, 200, -1,1);
    fhSelectedAsymmetry->SetXTitle("#it{E} (GeV)");
    fhSelectedAsymmetry->SetYTitle("#it{A} = ( #it{E}_{1} - #it{E}_{2} ) / ( #it{E}_{1} + #it{E}_{2} )");
    outputContainer->Add(fhSelectedAsymmetry);
    
    fhSplitE  = new TH1F
    ("hSplitE","Selected #pi^{0} (#eta) pairs energy sum of split sub-clusters",nptbins,ptmin,ptmax);
    fhSplitE->SetYTitle("counts");
    fhSplitE->SetXTitle("#it{E} (GeV)");
    outputContainer->Add(fhSplitE) ;
    
    fhSplitPt  = new TH1F
    ("hSplitPt","Selected #pi^{0} (#eta) pairs #it{p}_{T} sum of split sub-clusters",nptbins,ptmin,ptmax);
    fhSplitPt->SetYTitle("counts");
    fhSplitPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhSplitPt) ;
    
    
    fhSplitPtPhi  = new TH2F
    ("hSplitPtPhi","Selected #pi^{0} (#eta) pairs: sum split sub-cluster #it{p}_{T} vs #phi",nptbins,ptmin,ptmax, nphibins,phimin,phimax);
    fhSplitPtPhi->SetYTitle("#phi (rad)");
    fhSplitPtPhi->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhSplitPtPhi) ;
    
    fhSplitPtEta  = new TH2F
    ("hSplitPtEta","Selected #pi^{0} (#eta) pairs: sum split sub-cluster #it{p}_{T} vs #eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhSplitPtEta->SetYTitle("#eta");
    fhSplitPtEta->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhSplitPtEta) ;
    
    
    fhNLocMaxSplitPt = new TH2F("hNLocMaxSplitPt","Number of local maxima in cluster",
                                nptbins,ptmin,ptmax,20,0,20);
    fhNLocMaxSplitPt ->SetYTitle("#it{NLM}");
    fhNLocMaxSplitPt ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhNLocMaxSplitPt) ;
    
    
    fhMassSplitPt  = new TH2F
    ("hMassSplitPt","all pairs #it{M}: sum split sub-cluster #it{p}_{T} vs #it{M}",
     nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
    fhMassSplitPt->SetYTitle("#it{M} (GeV/#it{c}^{2})");
    fhMassSplitPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhMassSplitPt) ;
    
    fhSelectedMassSplitPt  = new TH2F
    ("hSelectedMassSplitPt","Selected #pi^{0} (#eta) pairs #it{M}: sum split sub-cluster #it{p}_{T} vs #it{M}",
     nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
    fhSelectedMassSplitPt->SetYTitle("#it{M} (GeV/#it{c}^{2})");
    fhSelectedMassSplitPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhSelectedMassSplitPt) ;
    
    if(IsDataMC())
    {
      fhMassSplitPtNoOverlap  = new TH2F
      ("hMassSplitPtNoOverlap","all pairs #it{M}: sum split sub-cluster #it{p}_{T} vs #it{M}, no overlap",
       nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
      fhMassSplitPtNoOverlap->SetYTitle("#it{M} (GeV/#it{c}^{2})");
      fhMassSplitPtNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMassSplitPtNoOverlap) ;
      
      fhSelectedMassSplitPtNoOverlap  = new TH2F
      ("hSelectedMassSplitPtNoOverlap","Selected #pi^{0} (#eta) pairs #it{M}: sum split sub-cluster #it{p}_{T} vs #it{M}, no overlap",
       nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
      fhSelectedMassSplitPtNoOverlap->SetYTitle("#it{M} (GeV/#it{c}^{2})");
      fhSelectedMassSplitPtNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhSelectedMassSplitPtNoOverlap) ;
      
      
      fhMCPi0PtRecoPtPrim  = new TH2F
      ("hMCPi0PtRecoPtPrim","#it{p}_{T,reco} vs #it{p}_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0PtRecoPtPrim ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCPi0PtRecoPtPrim ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCPi0PtRecoPtPrim ) ;
      
      fhMCPi0PtRecoPtPrimNoOverlap  = new TH2F
      ("hMCPi0PtRecoPtPrimNoOverlap","#it{p}_{T,reco} vs #it{p}_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0PtRecoPtPrimNoOverlap ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCPi0PtRecoPtPrimNoOverlap ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCPi0PtRecoPtPrimNoOverlap ) ;
      
      fhMCPi0SelectedPtRecoPtPrim  = new TH2F
      ("hMCPi0SelectedPtRecoPtPrim","#it{p}_{T,reco} vs #it{p}_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0SelectedPtRecoPtPrim ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCPi0SelectedPtRecoPtPrim ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCPi0SelectedPtRecoPtPrim ) ;
      
      fhMCPi0SelectedPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCPi0SelectedPtRecoPtPrimNoOverlap","#it{p}_{T,reco} vs #it{p}_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0SelectedPtRecoPtPrimNoOverlap ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCPi0SelectedPtRecoPtPrimNoOverlap ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCPi0SelectedPtRecoPtPrimNoOverlap ) ;
      
      
      fhMCPi0SplitPtRecoPtPrim  = new TH2F
      ("hMCPi0SplitPtRecoPtPrim","#it{p}_{T,reco} (split sum) vs #it{p}_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0SplitPtRecoPtPrim ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCPi0SplitPtRecoPtPrim ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCPi0SplitPtRecoPtPrim ) ;
      
      fhMCPi0SplitPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCPi0SplitPtRecoPtPrimNoOverlap","#it{p}_{T,reco} (split sum) vs #it{p}_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0SplitPtRecoPtPrimNoOverlap ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCPi0SplitPtRecoPtPrimNoOverlap ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCPi0SplitPtRecoPtPrimNoOverlap ) ;
      
      fhMCPi0SelectedSplitPtRecoPtPrim  = new TH2F
      ("hMCPi0SelectedSplitPtRecoPtPrim","#it{p}_{T,reco} (split sum) vs #it{p}_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0SelectedSplitPtRecoPtPrim ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCPi0SelectedSplitPtRecoPtPrim ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCPi0SelectedSplitPtRecoPtPrim ) ;
      
      fhMCPi0SelectedSplitPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCPi0SelectedSplitPtRecoPtPrimNoOverlap","#it{p}_{T,reco} (split sum) vs #it{p}_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0SelectedSplitPtRecoPtPrimNoOverlap ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCPi0SelectedSplitPtRecoPtPrimNoOverlap ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCPi0SelectedSplitPtRecoPtPrimNoOverlap ) ;
      
      fhMCEtaPtRecoPtPrim  = new TH2F
      ("hMCEtaPtRecoPtPrim","#it{p}_{T,reco} vs #it{p}_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaPtRecoPtPrim ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCEtaPtRecoPtPrim ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCEtaPtRecoPtPrim ) ;
      
      fhMCEtaPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCEtaPtRecoPtPrimNoOverlap","#it{p}_{T,reco} vs #it{p}_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaPtRecoPtPrimNoOverlap ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCEtaPtRecoPtPrimNoOverlap ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCEtaPtRecoPtPrimNoOverlap ) ;
      
      fhMCEtaSelectedPtRecoPtPrim  = new TH2F
      ("hMCEtaSelectedPtRecoPtPrim","#it{p}_{T,reco} vs #it{p}_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaSelectedPtRecoPtPrim ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCEtaSelectedPtRecoPtPrim ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCEtaSelectedPtRecoPtPrim ) ;
      
      fhMCEtaSelectedPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCEtaSelectedPtRecoPtPrimNoOverlap","#it{p}_{T,reco} vs #it{p}_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaSelectedPtRecoPtPrimNoOverlap ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCEtaSelectedPtRecoPtPrimNoOverlap ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCEtaSelectedPtRecoPtPrimNoOverlap ) ;
      
      
      fhMCEtaSplitPtRecoPtPrim  = new TH2F
      ("hMCEtaSplitPtRecoPtPrim","#it{p}_{T,reco} (split sum) vs #it{p}_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaSplitPtRecoPtPrim ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCEtaSplitPtRecoPtPrim ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCEtaSplitPtRecoPtPrim ) ;
      
      fhMCEtaSplitPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCEtaSplitPtRecoPtPrimNoOverlap","#it{p}_{T,reco} (split sum) vs #it{p}_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaSplitPtRecoPtPrimNoOverlap ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCEtaSplitPtRecoPtPrimNoOverlap ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCEtaSplitPtRecoPtPrimNoOverlap ) ;
      
      fhMCEtaSelectedSplitPtRecoPtPrim  = new TH2F
      ("hMCEtaSelectedSplitPtRecoPtPrim","#it{p}_{T,reco} (split sum) vs #it{p}_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaSelectedSplitPtRecoPtPrim ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCEtaSelectedSplitPtRecoPtPrim ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCEtaSelectedSplitPtRecoPtPrim ) ;
      
      fhMCEtaSelectedSplitPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCEtaSelectedSplitPtRecoPtPrimNoOverlap","#it{p}_{T,reco} (split sum) vs #it{p}_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaSelectedSplitPtRecoPtPrimNoOverlap ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
      fhMCEtaSelectedSplitPtRecoPtPrimNoOverlap ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
      outputContainer->Add(fhMCEtaSelectedSplitPtRecoPtPrimNoOverlap ) ;
      
      
      for(Int_t inlm = 0; inlm < 3; inlm++)
      {
        fhMCPi0PtRecoPtPrimLocMax[inlm]  = new TH2F
        (Form("hMCPi0PtRecoPtPrimLocMax%d",inlm+1),Form("#it{p}_{T,reco} vs #it{p}_{T,gen}, %s",nlm[inlm].Data()),
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCPi0PtRecoPtPrimLocMax[inlm] ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
        fhMCPi0PtRecoPtPrimLocMax[inlm] ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
        outputContainer->Add(fhMCPi0PtRecoPtPrimLocMax[inlm] ) ;
        
        fhMCPi0SelectedPtRecoPtPrimLocMax[inlm]  = new TH2F
        (Form("hMCPi0SelectedPtRecoPtPrimLocMax%d",inlm+1),Form("#it{p}_{T,reco} vs #it{p}_{T,gen}, %s",nlm[inlm].Data()),
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCPi0SelectedPtRecoPtPrimLocMax[inlm] ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
        fhMCPi0SelectedPtRecoPtPrimLocMax[inlm] ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
        outputContainer->Add(fhMCPi0SelectedPtRecoPtPrimLocMax[inlm] ) ;
        
        fhMCPi0SplitPtRecoPtPrimLocMax[inlm]  = new TH2F
        (Form("hMCPi0SplitPtRecoPtPrimLocMax%d",inlm+1),Form("#it{p}_{T,reco} (split sum) vs #it{p}_{T,gen}, %s",nlm[inlm].Data()),
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCPi0SplitPtRecoPtPrimLocMax[inlm] ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
        fhMCPi0SplitPtRecoPtPrimLocMax[inlm] ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
        outputContainer->Add(fhMCPi0SplitPtRecoPtPrimLocMax[inlm] ) ;
        
        fhMCPi0SelectedSplitPtRecoPtPrimLocMax[inlm]  = new TH2F
        (Form("hMCPi0SelectedSplitPtRecoPtPrimLocMax%d",inlm+1),Form("#it{p}_{T,reco} (split sum) vs #it{p}_{T,gen}, %s",nlm[inlm].Data()),
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCPi0SelectedSplitPtRecoPtPrimLocMax[inlm] ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
        fhMCPi0SelectedSplitPtRecoPtPrimLocMax[inlm] ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
        outputContainer->Add(fhMCPi0SelectedSplitPtRecoPtPrimLocMax[inlm] ) ;
        
        fhMCEtaPtRecoPtPrimLocMax[inlm]  = new TH2F
        (Form("hMCEtaPtRecoPtPrimLocMax%d",inlm+1),Form("#it{p}_{T,reco} vs #it{p}_{T,gen}, %s",nlm[inlm].Data()),
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCEtaPtRecoPtPrimLocMax[inlm] ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
        fhMCEtaPtRecoPtPrimLocMax[inlm] ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
        outputContainer->Add(fhMCEtaPtRecoPtPrimLocMax[inlm] ) ;
        
        fhMCEtaSelectedPtRecoPtPrimLocMax[inlm]  = new TH2F
        (Form("hMCEtaSelectedPtRecoPtPrimLocMax%d",inlm+1),Form("#it{p}_{T,reco} vs #it{p}_{T,gen}, %s",nlm[inlm].Data()),
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCEtaSelectedPtRecoPtPrimLocMax[inlm] ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
        fhMCEtaSelectedPtRecoPtPrimLocMax[inlm] ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
        outputContainer->Add(fhMCEtaSelectedPtRecoPtPrimLocMax[inlm] ) ;
        
        fhMCEtaSplitPtRecoPtPrimLocMax[inlm]  = new TH2F
        (Form("hMCEtaSplitPtRecoPtPrimLocMax%d",inlm+1),Form("#it{p}_{T,reco} (split sum) vs #it{p}_{T,gen}, %s",nlm[inlm].Data()),
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCEtaSplitPtRecoPtPrimLocMax[inlm] ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
        fhMCEtaSplitPtRecoPtPrimLocMax[inlm] ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
        outputContainer->Add(fhMCEtaSplitPtRecoPtPrimLocMax[inlm] ) ;
        
        fhMCEtaSelectedSplitPtRecoPtPrimLocMax[inlm]  = new TH2F
        (Form("hMCEtaSelectedSplitPtRecoPtPrimLocMax%d",inlm+1),Form("#it{p}_{T,reco} (split sum) vs #it{p}_{T,gen}, %s",nlm[inlm].Data()),
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhMCEtaSelectedSplitPtRecoPtPrimLocMax[inlm] ->SetYTitle("#it{p}_{T,gen} (GeV/#it{c})");
        fhMCEtaSelectedSplitPtRecoPtPrimLocMax[inlm] ->SetXTitle("#it{p}_{T,reco} (GeV/#it{c})");
        outputContainer->Add(fhMCEtaSelectedSplitPtRecoPtPrimLocMax[inlm] ) ;
      }
      
      for(Int_t i = 0; i < fgkNmcTypes; i++)
      {
        fhMCPtAsymmetry[i]  = new TH2F (Form("hEAsymmetry_MC%s",pname[i].Data()),
                                        Form("cluster from %s : #it{A} = ( #it{E}_{1} - #it{E}_{2} ) / ( #it{E}_{1} + #it{E}_{2} ) vs #it{E}",ptype[i].Data()),
                                        nptbins,ptmin,ptmax, 200,-1,1);
        fhMCPtAsymmetry[i]->SetXTitle("#it{E} (GeV)");
        fhMCPtAsymmetry[i]->SetYTitle("#it{A} = ( #it{E}_{1} - #it{E}_{2} ) / ( #it{E}_{1} + #it{E}_{2} )");
        outputContainer->Add(fhMCPtAsymmetry[i]);
        
        fhMCSplitE[i]  = new TH1F
        (Form("hSplitE_MC%s",pname[i].Data()),
         Form("cluster from %s, energy sum of split sub-clusters",ptype[i].Data()),
         nptbins,ptmin,ptmax);
        fhMCSplitE[i]->SetYTitle("counts");
        fhMCSplitE[i]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCSplitE[i]) ;
        
        fhMCSplitPt[i]  = new TH1F
        (Form("hSplitPt_MC%s",pname[i].Data()),
         Form("cluster from %s, #it{p}_{T} sum of split sub-clusters",ptype[i].Data()),
         nptbins,ptmin,ptmax);
        fhMCSplitPt[i]->SetYTitle("counts");
        fhMCSplitPt[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCSplitPt[i]) ;
        
        
        fhMCSplitPtPhi[i]  = new TH2F
        (Form("hSplitPtPhi_MC%s",pname[i].Data()),
         Form("Identified as #pi^{0} (#eta), cluster from %s",ptype[i].Data()),
         nptbins,ptmin,ptmax,nphibins,phimin,phimax);
        fhMCSplitPtPhi[i]->SetYTitle("#phi");
        fhMCSplitPtPhi[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCSplitPtPhi[i]) ;
        
        fhMCSplitPtEta[i]  = new TH2F
        (Form("hSplitPtEta_MC%s",pname[i].Data()),
         Form("Identified as #pi^{0} (#eta), cluster from %s",
              ptype[i].Data()),nptbins,ptmin,ptmax,netabins,etamin,etamax);
        fhMCSplitPtEta[i]->SetYTitle("#eta");
        fhMCSplitPtEta[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCSplitPtEta[i]) ;
        
        
        fhMCNLocMaxSplitPt[i] = new TH2F
        (Form("hNLocMaxSplitPt_MC%s",pname[i].Data()),
         Form("cluster from %s, #it{p}_{T} sum of split sub-clusters, for NLM",ptype[i].Data()),
         nptbins,ptmin,ptmax,20,0,20);
        fhMCNLocMaxSplitPt[i] ->SetYTitle("#it{NLM}");
        fhMCNLocMaxSplitPt[i] ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCNLocMaxSplitPt[i]) ;
        
        fhMCMassSplitPt[i]  = new TH2F
        (Form("hMassSplitPt_MC%s",pname[i].Data()),
         Form("all pairs #it{M}: split #it{p}_{T} vs #it{M} from %s",ptype[i].Data()),
         nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
        fhMCMassSplitPt[i]->SetYTitle("#it{M} (GeV/#it{c}^{2})");
        fhMCMassSplitPt[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCMassSplitPt[i]) ;
        
        fhMCSelectedMassSplitPt[i]  = new TH2F
        (Form("hSelectedMassSplitPt_MC%s",pname[i].Data()),
         Form("Selected #pi^{0} (#eta) pairs #it{M}: split #it{p}_{T} vs #it{M} from %s",ptype[i].Data()),
         nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
        fhMCSelectedMassSplitPt[i]->SetYTitle("#it{M} (GeV/#it{c}^{2})");
        fhMCSelectedMassSplitPt[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCSelectedMassSplitPt[i]) ;
        
        fhMCMassSplitPtNoOverlap[i]  = new TH2F
        (Form("hMassSplitPtNoOverlap_MC%s",pname[i].Data()),
         Form("all pairs #it{M}: split #it{p}_{T} vs #it{M} from %s, no overlap",ptype[i].Data()),
         nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
        fhMCMassSplitPtNoOverlap[i]->SetYTitle("#it{M} (GeV/#it{c}^{2})");
        fhMCMassSplitPtNoOverlap[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCMassSplitPtNoOverlap[i]) ;
        
        fhMCSelectedMassSplitPtNoOverlap[i]  = new TH2F
        (Form("hSelectedMassSplitPtNoOverlap_MC%s",pname[i].Data()),
         Form("Selected #pi^{0} (#eta) pairs #it{M}: split #it{p}_{T} vs #it{M} from %s, no overlap",ptype[i].Data()),
         nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
        fhMCSelectedMassSplitPtNoOverlap[i]->SetYTitle("#it{M} (GeV/#it{c}^{2})");
        fhMCSelectedMassSplitPtNoOverlap[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCSelectedMassSplitPtNoOverlap[i]) ;
      }
    }
  }
  
  if(fAnaType==kSSCalo && fFillSelectClHisto && !fFillOnlySimpleSSHisto )
  {
    for(Int_t i = 0; i< 3; i++)
    {
      fhPtAsymmetryLocMax[i]  = new TH2F(Form("hEAsymmetryLocMax%d",i+1),
                                         Form("Selected #pi^{0} (#eta) pairs: #it{p}_{T} vs #it{A} = ( #it{E}_{1} - #it{E}_{2} ) / ( #it{E}_{1} + #it{E}_{2} ), %s",nlm[i].Data()),
                                         nptbins,ptmin,ptmax,200, -1,1);
      fhPtAsymmetryLocMax[i]->SetYTitle("#it{A} = ( #it{E}_{1} - #it{E}_{2} ) / ( #it{E}_{1} + #it{E}_{2} )");
      fhPtAsymmetryLocMax[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtAsymmetryLocMax[i]) ;
    }
    
    for(Int_t ie = 0; ie < 7; ie++)
    {
      
      fhAsymmetryLambda0[ie] = new TH2F (Form("hAsymmetryLambda0_EBin%d",ie),
                                         Form("#lambda_{0}^{2} vs A for %d < #it{E} < %d GeV",bin[ie],bin[ie+1]),
                                         ssbins,ssmin,ssmax , 200,-1,1);
      fhAsymmetryLambda0[ie]->SetXTitle("#lambda_{0}^{2}");
      fhAsymmetryLambda0[ie]->SetYTitle("A = ( E1 - E2 ) / ( E1 + E2 )");
      outputContainer->Add(fhAsymmetryLambda0[ie]);
      
      fhAsymmetryDispEta[ie] = new TH2F (Form("hAsymmetryDispEta_EBin%d",ie),
                                         Form("#sigma^{2}_{#eta #eta} vs #it{A} for %d < #it{E} < %d GeV",bin[ie],bin[ie+1]),
                                         ssbins,ssmin,ssmax , 200,-1,1);
      fhAsymmetryDispEta[ie]->SetXTitle("#sigma^{2}_{#eta #eta}");
      fhAsymmetryDispEta[ie]->SetYTitle("#it{A} = ( #it{E}_{1} - #it{E}_{2} ) / ( #it{E}_{1} + #it{E}_{2} )");
      outputContainer->Add(fhAsymmetryDispEta[ie]);
      
      fhAsymmetryDispPhi[ie] = new TH2F (Form("hAsymmetryDispPhi_EBin%d",ie),
                                         Form("#sigma^{2}_{#phi #phi} vs #it{A} for %d < #it{E} < %d GeV",bin[ie],bin[ie+1]),
                                         ssbins,ssmin,ssmax , 200,-1,1);
      fhAsymmetryDispPhi[ie]->SetXTitle("#sigma^{2}_{#phi #phi}");
      fhAsymmetryDispPhi[ie]->SetYTitle("#it{A} = ( #it{E}_{1} - #it{E}_{2} ) / ( #it{E}_{1} + #it{E}_{2} )");
      outputContainer->Add(fhAsymmetryDispPhi[ie]);
    }
    
    
    if(IsDataMC())
    {
      for(Int_t i = 0; i < fgkNmcTypes; i++)
      {
        for(Int_t ie = 0; ie < 7; ie++)
        {
          fhMCAsymmetryLambda0[ie][i] = new TH2F (Form("hMCAsymmetryLambda0_EBin%d_MC%s",ie,pname[i].Data()),
                                                  Form("cluster from %s : #lambda_{0}^{2} vs A for %d < #it{E} < %d GeV",ptype[i].Data(),bin[ie],bin[ie+1]),
                                                  ssbins,ssmin,ssmax , 200,-1,1);
          fhMCAsymmetryLambda0[ie][i]->SetXTitle("#lambda_{0}^{2}");
          fhMCAsymmetryLambda0[ie][i]->SetYTitle("#it{A} = ( #it{E}_{1} - #it{E}_{2} ) / ( #it{E}_{1} + #it{E}_{2} )");
          outputContainer->Add(fhMCAsymmetryLambda0[ie][i]);
          
          fhMCAsymmetryDispEta[ie][i] = new TH2F (Form("hMCAsymmetryDispEta_EBin%d_MC%s",ie,pname[i].Data()),
                                                  Form("cluster from %s : #sigma^{2}_{#eta #eta} vs #it{A} for %d < #it{E} < %d GeV",ptype[i].Data(),bin[ie],bin[ie+1]),
                                                  ssbins,ssmin,ssmax , 200,-1,1);
          fhMCAsymmetryDispEta[ie][i]->SetXTitle("#sigma^{2}_{#eta #eta}");
          fhMCAsymmetryDispEta[ie][i]->SetYTitle("#it{A} = ( #it{E}_{1} - #it{E}_{2} ) / ( #it{E}_{1} + #it{E}_{2} )");
          outputContainer->Add(fhMCAsymmetryDispEta[ie][i]);
          
          fhMCAsymmetryDispPhi[ie][i] = new TH2F (Form("hMCAsymmetryDispPhi_EBin%d_MC%s",ie,pname[i].Data()),
                                                  Form("cluster from %s : #sigma^{2}_{#phi #phi} vs #it{A} for %d < #it{E} < %d GeV",ptype[i].Data(),bin[ie],bin[ie+1]),
                                                  ssbins,ssmin,ssmax , 200,-1,1);
          fhMCAsymmetryDispPhi[ie][i]->SetXTitle("#sigma^{2}_{#phi #phi}");
          fhMCAsymmetryDispPhi[ie][i]->SetYTitle("#it{A} = ( #it{E}_{1} - #it{E}_{2} ) / ( #it{E}_{1} + #it{E}_{2} )");
          outputContainer->Add(fhMCAsymmetryDispPhi[ie][i]);
        }
      }
    }
  }
  
  if(IsPileUpAnalysisOn())
  {
    
    TString pileUpName[] = {"SPD","EMCAL","SPDOrEMCAL","SPDAndEMCAL","SPDAndNotEMCAL","EMCALAndNotSPD","NotSPDAndNotEMCAL"} ;
    
    for(Int_t i = 0 ; i < 7 ; i++)
    {
      fhPtPileUp[i]  = new TH1F(Form("hPtPileUp%s",pileUpName[i].Data()),
                                Form("Selected #pi^{0} (#eta) #it{p}_{T} distribution, %s Pile-Up event",pileUpName[i].Data()), nptbins,ptmin,ptmax);
      fhPtPileUp[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPileUp[i]);
      
      fhPtCellTimePileUp[i]  = new TH2F(Form("hPtCellTimePileUp%s",pileUpName[i].Data()),
                                        Form("Pt vs cell time in cluster, %s Pile-Up event",pileUpName[i].Data()),
                                        nptbins,ptmin,ptmax,ntimptbins,timemin,timemax);
      fhPtCellTimePileUp[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPtCellTimePileUp[i]->SetYTitle("#it{t}_{cell} (ns)");
      outputContainer->Add(fhPtCellTimePileUp[i]);
      
      fhPtTimeDiffPileUp[i]  = new TH2F(Form("hPtTimeDiffPileUp%s",pileUpName[i].Data()),
                                        Form("Pt vs t_{max}-t_{cell} in cluster, %s Pile-Up event",pileUpName[i].Data()),
                                        nptbins,ptmin,ptmax,400,-200,200);
      fhPtTimeDiffPileUp[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPtTimeDiffPileUp[i]->SetYTitle("#it{t}_{max}-#it{t}_{cell} (ns)");
      outputContainer->Add(fhPtTimeDiffPileUp[i]);
    }
    
    fhTimePtNoCut  = new TH2F ("hTimePt_NoCut","#it{t} of cluster vs #it{E} of clusters, no cut", nptbins,ptmin,ptmax, ntimptbins,timemin,timemax);
    fhTimePtNoCut->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhTimePtNoCut->SetYTitle("#it{t} (ns)");
    outputContainer->Add(fhTimePtNoCut);
    
    fhTimePtSPD  = new TH2F ("hTimePt_SPD","#it{t} of cluster vs #it{E} of clusters, SPD cut", nptbins,ptmin,ptmax, ntimptbins,timemin,timemax);
    fhTimePtSPD->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhTimePtSPD->SetYTitle("#it{t} (ns)");
    outputContainer->Add(fhTimePtSPD);
    
    fhTimePtSPDMulti  = new TH2F ("hTimePt_SPDMulti","time of cluster vs #it{E} of clusters, SPD multi cut", nptbins,ptmin,ptmax, ntimptbins,timemin,timemax);
    fhTimePtSPDMulti->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhTimePtSPDMulti->SetYTitle("#it{t} (ns)");
    outputContainer->Add(fhTimePtSPDMulti);
    
    fhTimeNPileUpVertSPD  = new TH2F ("hTime_NPileUpVertSPD","#it{t} of cluster vs #it{N} pile-up SPD vertex", ntimptbins,timemin,timemax,50,0,50);
    fhTimeNPileUpVertSPD->SetYTitle("# vertex ");
    fhTimeNPileUpVertSPD->SetXTitle("#it{t} (ns)");
    outputContainer->Add(fhTimeNPileUpVertSPD);
    
    fhTimeNPileUpVertTrack  = new TH2F ("hTime_NPileUpVertTracks","#it{t} of cluster vs #it{N} pile-up Tracks vertex", ntimptbins,timemin,timemax, 50,0,50 );
    fhTimeNPileUpVertTrack->SetYTitle("# vertex ");
    fhTimeNPileUpVertTrack->SetXTitle("#it{t} (ns)");
    outputContainer->Add(fhTimeNPileUpVertTrack);
    
    fhTimeNPileUpVertContributors  = new TH2F ("hTime_NPileUpVertContributors","#it{t} of cluster vs #it{N} constributors to pile-up SPD vertex", ntimptbins,timemin,timemax,50,0,50);
    fhTimeNPileUpVertContributors->SetYTitle("# vertex ");
    fhTimeNPileUpVertContributors->SetXTitle("#it{t} (ns)");
    outputContainer->Add(fhTimeNPileUpVertContributors);
    
    fhTimePileUpMainVertexZDistance  = new TH2F ("hTime_PileUpMainVertexZDistance","#it{t} of cluster vs distance in #it{Z} pile-up SPD vertex - main SPD vertex",ntimptbins,timemin,timemax,100,0,50);
    fhTimePileUpMainVertexZDistance->SetYTitle("distance #it{Z} (cm) ");
    fhTimePileUpMainVertexZDistance->SetXTitle("#it{t} (ns)");
    outputContainer->Add(fhTimePileUpMainVertexZDistance);
    
    fhTimePileUpMainVertexZDiamond  = new TH2F ("hTime_PileUpMainVertexZDiamond","#it{t} of cluster vs distance in #it{Z} pile-up SPD vertex - z diamond",ntimptbins,timemin,timemax,100,0,50);
    fhTimePileUpMainVertexZDiamond->SetYTitle("diamond distance #it{Z} (cm) ");
    fhTimePileUpMainVertexZDiamond->SetXTitle("#it{t} (ns)");
    outputContainer->Add(fhTimePileUpMainVertexZDiamond);
    
    fhPtNPileUpSPDVtx  = new TH2F ("hPt_NPileUpVertSPD","#it{p}_{T} of cluster vs #it{N} pile-up SPD vertex",
                                   nptbins,ptmin,ptmax,20,0,20);
    fhPtNPileUpSPDVtx->SetYTitle("# vertex ");
    fhPtNPileUpSPDVtx->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtNPileUpSPDVtx);
    
    fhPtNPileUpTrkVtx  = new TH2F ("hPt_NPileUpVertTracks","#it{p}_{T} of cluster vs #it{N} pile-up Tracks vertex",
                                   nptbins,ptmin,ptmax, 20,0,20 );
    fhPtNPileUpTrkVtx->SetYTitle("# vertex ");
    fhPtNPileUpTrkVtx->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtNPileUpTrkVtx);
    
    fhPtNPileUpSPDVtxTimeCut  = new TH2F ("hPt_NPileUpVertSPD_TimeCut","#it{p}_{T} of cluster vs N pile-up SPD vertex, |tof| < 25 ns",
                                          nptbins,ptmin,ptmax,20,0,20);
    fhPtNPileUpSPDVtxTimeCut->SetYTitle("# vertex ");
    fhPtNPileUpSPDVtxTimeCut->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtNPileUpSPDVtxTimeCut);
    
    fhPtNPileUpTrkVtxTimeCut  = new TH2F ("hPt_NPileUpVertTracks_TimeCut","#it{p}_{T} of cluster vs N pile-up Tracks vertex, |tof| < 25 ns",
                                          nptbins,ptmin,ptmax, 20,0,20 );
    fhPtNPileUpTrkVtxTimeCut->SetYTitle("# vertex ");
    fhPtNPileUpTrkVtxTimeCut->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtNPileUpTrkVtxTimeCut);
    
    fhPtNPileUpSPDVtxTimeCut2  = new TH2F ("hPt_NPileUpVertSPD_TimeCut2","#it{p}_{T} of cluster vs N pile-up SPD vertex, -25 < tof < 75 ns",
                                           nptbins,ptmin,ptmax,20,0,20);
    fhPtNPileUpSPDVtxTimeCut2->SetYTitle("# vertex ");
    fhPtNPileUpSPDVtxTimeCut2->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtNPileUpSPDVtxTimeCut2);
    
    fhPtNPileUpTrkVtxTimeCut2  = new TH2F ("hPt_NPileUpVertTracks_TimeCut2","#it{p}_{T} of cluster vs N pile-up Tracks vertex, -25 < tof < 75 ns",
                                           nptbins,ptmin,ptmax, 20,0,20 );
    fhPtNPileUpTrkVtxTimeCut2->SetYTitle("# vertex ");
    fhPtNPileUpTrkVtxTimeCut2->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtNPileUpTrkVtxTimeCut2);
    
  }
  
  // Keep neutral meson selection histograms if requiered
  // Setting done in AliNeutralMesonSelection
  
  if(fAnaType!=kSSCalo && GetNeutralMesonSelection())
  {
    TList * nmsHistos = GetNeutralMesonSelection()->GetCreateOutputObjects() ;
    
    if(GetNeutralMesonSelection()->AreNeutralMesonSelectionHistosKept())
      for(Int_t i = 0; i < nmsHistos->GetEntries(); i++) outputContainer->Add(nmsHistos->At(i)) ;
    
    delete nmsHistos;
  }
  
  return outputContainer ;
}

//_____________________________________________
/// Assign MC index depending on MC bit set.
/// Index used when filling array of MC dependent histograms.
//_____________________________________________
Int_t AliAnaPi0EbE::GetMCIndex(const Int_t tag)
{
  if       ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)  )
  {
    return kmcPi0 ;
  }//pi0
  else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)  )
  {
    return kmcEta ;
  }//eta
  else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt) ||
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation) ||
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR))
  {
    return kmcPhoton ;
  }//direct photon
  else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) &&
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay) )
  {
    return kmcPi0Decay ;
  }//decay photon from pi0
  else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) &&
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) )
  {
    return kmcEtaDecay ;
  }//decay photon from eta
  else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) )
  {
    return kmcOtherDecay ;
  }//decay photon from other than eta or pi0
  else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron))
  {
    return kmcElectron ;
  }//electron
  else
  {
    return kmcHadron ;
  }//other particles
}

//__________________________________________________________________
/// Check the labels of pair in case mother was same pi0 or eta
/// Set the new AOD accordingly. Only valid for kIMCalo and kIMCaloTracks
//__________________________________________________________________
void AliAnaPi0EbE::HasPairSameMCMother(Int_t label1 , Int_t label2,
                                       Int_t tag1   , Int_t tag2,
                                       Int_t & label, Int_t & tag)
{
  if(label1 < 0 || label2 < 0 || label1 == label2) return ;
  
  AliDebug(0,Form("Origin of: photon1 %d; photon2 %d",tag1, tag2));
  
  if( (GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCPi0Decay) &&
       GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCPi0Decay)    ) ||
      (GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCEtaDecay) &&
       GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCEtaDecay)    )
     )
  {
    Int_t pdg1    = -1;//, pdg2    = -1;
    Int_t ndaugh1 = -1, ndaugh2 = -1;
    //Check if pi0/eta mother is the same
    if(GetReader()->ReadStack())
    {
      if ( label1 >= 0 )
      {
        TParticle * mother1 = GetMCStack()->Particle(label1);//photon in kine tree
        label1 = mother1->GetFirstMother();
        mother1 = GetMCStack()->Particle(label1);//pi0
        pdg1=mother1->GetPdgCode();
        ndaugh1 = mother1->GetNDaughters();
      }
      if ( label2 >= 0 )
      {
        TParticle * mother2 = GetMCStack()->Particle(label2);//photon in kine tree
        label2 = mother2->GetFirstMother();
        mother2 = GetMCStack()->Particle(label2);//pi0
        //pdg2=mother2->GetPdgCode();
        ndaugh2 = mother2->GetNDaughters();
      }
    } // STACK
    else if(GetReader()->ReadAODMCParticles())
    {//&& (input > -1)){
      if ( label1 >= 0 )
      {
        AliAODMCParticle * mother1 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles())->At(label1);//photon in kine tree
        label1 = mother1->GetMother();
        mother1 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles())->At(label1);//pi0
        pdg1=mother1->GetPdgCode();
        ndaugh1 = mother1->GetNDaughters();
      }
      if ( label2 >= 0 )
      {
        AliAODMCParticle * mother2 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles())->At(label2);//photon in kine tree
        label2 = mother2->GetMother();
        mother2 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles())->At(label2);//pi0
        //pdg2=mother2->GetPdgCode();
        ndaugh2 = mother2->GetNDaughters();
      }
    } // AOD
    
//  printf("mother1 %d, mother2 %d\n",label1,label2);
    if ( label1 == label2 && label1 >= 0  && ndaugh1 == ndaugh2 && ndaugh1 == 2 )
    {
      label = label1;
      
      Double_t angle = fMomentum2.Angle(fMomentum1.Vect());
      Double_t mass  = (fMomentum1+fMomentum2).M();
      Double_t epair = (fMomentum1+fMomentum2).E();
      
      if ( pdg1 == 111 )
      {
        fhMassPairMCPi0 ->Fill(epair, mass , GetEventWeight());
        fhAnglePairMCPi0->Fill(epair, angle, GetEventWeight());
          
        GetMCAnalysisUtils()->SetTagBit(tag, AliMCAnalysisUtils::kMCPi0);
          
//      printf("Real pi0 pair: pt %f, mass %f\n",(fMomentum1+fMomentum2).Pt(),mass);
//      printf(" Lab1 %d (%d), lab2 %d (%d), pdg1 %d, pdg2 %d, Is In calo %d, %d, Is lost %d, %d\n",
//             label1,photon1->GetLabel(),label2,photon2->GetLabel(), pdg1, pdg2,
//             GetMCAnalysisUtils()->CheckTagBit(photon1->GetTag(),AliMCAnalysisUtils::kMCDecayPairInCalo),
//             GetMCAnalysisUtils()->CheckTagBit(photon2->GetTag(),AliMCAnalysisUtils::kMCDecayPairInCalo),
//             GetMCAnalysisUtils()->CheckTagBit(photon1->GetTag(),AliMCAnalysisUtils::kMCDecayPairLost),
//             GetMCAnalysisUtils()->CheckTagBit(photon2->GetTag(),AliMCAnalysisUtils::kMCDecayPairLost));
      }
      else  if ( pdg1 == 221 )
      {
        fhMassPairMCEta ->Fill(epair, mass , GetEventWeight());
        fhAnglePairMCEta->Fill(epair, angle, GetEventWeight());
          
        GetMCAnalysisUtils()->SetTagBit(tag, AliMCAnalysisUtils::kMCEta);
          
//      printf("Real eta pair\n");
      }
      
    } // same label
  } // both from eta or pi0 decay
}

//_______________________
/// Do some checks to run the analysis.
/// If no calorimeter list of clusters provided, abort.
//_______________________
void AliAnaPi0EbE::Init()
{
  if       ( GetCalorimeter() == kPHOS  && !GetReader()->IsPHOSSwitchedOn()  && NewOutputAOD() )
    AliFatal("STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!");
  else  if ( GetCalorimeter() == kEMCAL && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD() )
    AliFatal("STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!");
}

//_________________________________
/// Initialize the parameters of the analysis with default values.
//_________________________________
void AliAnaPi0EbE::InitParameters()
{
  AddToHistogramsName("AnaPi0EbE_");
  
  fInputAODGammaConvName = "PhotonsCTS" ;
  fAnaType = kIMCalo ;
  fMinDist  = 2.;
  fMinDist2 = 4.;
  fMinDist3 = 5.;
  
  fNLMECutMin[0] = 10.;
  fNLMECutMin[1] = 6. ;
  fNLMECutMin[2] = 6. ;
  
  fIsoCandMinPt = 8;
  fR            = 0.4;
  
  fM02MinCutForIM = -1 ;
  fM02MaxCutForIM = -1 ;
}

//_______________________________________
/// Do analysis and fill AOD branch with selected pi0/eta pairs/clusters.
/// Execute here the requested analysis type: kIMCalo, kSSCalo or kIMCaloTracks.
//_______________________________________
void  AliAnaPi0EbE::MakeAnalysisFillAOD()
{
  AliDebug(1,Form("Start analysis type %d",fAnaType));
  
  switch(fAnaType)
  {
    case kIMCalo:
      MakeInvMassInCalorimeter();
      break;
      
    case kSSCalo:
      MakeShowerShapeIdentification();
      break;
      
    case kIMCaloTracks:
      MakeInvMassInCalorimeterAndCTS();
      break;
  }
  
  AliDebug(1,"End");
}

//____________________________________________
/// Read photon list from AOD, produced in class AliAnaPhoton.
/// Check if 2 photons have the mass of the pi0/eta, depending on the
/// settings of AliNeutralMesonSelection. If the pair is selected create the
/// particle AOD object with the pair and fill the AOD branch.
/// Fill also some control histograms.
//____________________________________________
void  AliAnaPi0EbE::MakeInvMassInCalorimeter()
{
  if(!GetInputAODBranch())
  {
    AliFatal(Form("No input calo photons in AOD with name branch < %s >, STOP",GetInputAODName().Data()));
    return; // coverity
  }
  
  // Create obj array of clusters in case it is requested
  TObjArray *clusters = 0;
  if( ( fFillSelectClHisto|| IsPileUpAnalysisOn() )&&
     GetReader()->GetDataType()!=AliCaloTrackReader::kMC )
  {
    if     (GetCalorimeter()==kEMCAL) clusters = GetEMCALClusters();
    else if(GetCalorimeter()==kPHOS)  clusters = GetPHOSClusters() ;
  }
  
  AliVCluster *cluster1 = 0;
  AliVCluster *cluster2 = 0;
  
  Int_t nphoton = GetInputAODBranch()->GetEntriesFast();
  for( Int_t iphoton = 0; iphoton < nphoton-1; iphoton++ )
  {
    AliAODPWG4Particle * photon1 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(iphoton));
    
    // Vertex cut in case of mixed events
    // Not really in use or tested recently, remove ...
    Int_t evtIndex1 = 0 ;
    if(GetMixedEvent())
    {
      evtIndex1 = GetMixedEvent()->EventIndexForCaloCluster(photon1->GetCaloLabel(0)) ;
      if(TMath::Abs(GetVertex(evtIndex1)[2]) > GetZvertexCut()) continue ;  //vertex cut
    }
    
    // Get original cluster if requested
    Int_t iclus = -1;
    if( ( fFillSelectClHisto || IsPileUpAnalysisOn() ) &&
       GetReader()->GetDataType()!=AliCaloTrackReader::kMC)
    {
      cluster1 = FindCluster(clusters,photon1->GetCaloLabel(0),iclus);
      
      if(!cluster1) AliWarning("First cluster not found");
    }
    
    // Get kinematics and other parameters
    fMomentum1 = *(photon1->Momentum());
    
    Float_t       e1 = photon1->E();
    Float_t     tof1 = photon1->GetTime();
    Int_t   nMaxima1 = photon1->GetNLM();
    Int_t       lab1 = photon1->GetLabel();
    Int_t       tag1 = photon1->GetTag();
    Bool_t isolated1 = ((AliAODPWG4ParticleCorrelation*) photon1)->IsIsolated();
    
    for( Int_t jphoton = iphoton+1; jphoton < nphoton; jphoton++ )
    {
      AliAODPWG4Particle * photon2 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(jphoton));
      
      // Do analysis only when one of the decays is isolated
      // Run AliAnaParticleIsolation before
      if(fSelectIsolatedDecay)
      {
        Bool_t isolated2 = ((AliAODPWG4ParticleCorrelation*) photon2)->IsIsolated();
        if(!isolated1 && !isolated2) continue;
      }
      
      // Vertex cut in case of mixed events
      // Not really in use or tested recently, remove ...
      Int_t evtIndex2 = 0 ;
      if(GetMixedEvent())
      {
        evtIndex2 = GetMixedEvent()->EventIndexForCaloCluster(photon2->GetCaloLabel(0)) ;
        
        if(evtIndex1 == evtIndex2)
          continue ;
        
        if(TMath::Abs(GetVertex(evtIndex2)[2]) > GetZvertexCut()) continue ;  //vertex cut
      }
      
      // Get kinematics and other parameters
      fMomentum2 = *(photon2->Momentum());
      
      Float_t     e2 = photon2->E();
      Float_t   tof2 = photon2->GetTime();
      Int_t nMaxima2 = photon2->GetNLM();
      
      //
      // Skip pairs of EMCal and DCal clusters
      Double_t angleOp = fMomentum1.Angle(fMomentum2.Vect());
      
      if(angleOp > DegToRad(100.)) continue; 
      
      //
      // Select clusters with good time window difference
      Double_t t12diff = tof1-tof2;
        
      fhEPairDiffTime->Fill(e1+e2, t12diff, GetEventWeight());
        
      if(TMath::Abs(t12diff) > GetPairTimeCut()) continue;
      
      //
      // Play with the MC stack if available
      //
      Int_t mcIndex = kmcHadron;
      Int_t tag     = 0;
      Int_t label   =-1;
      if( IsDataMC() )
      {
        HasPairSameMCMother(lab1, photon2->GetLabel(),
                            tag1, photon2->GetTag(),
                            label, tag) ;
        mcIndex = GetMCIndex(tag);
      }
      
      // Check the invariant mass for different selection on the local maxima
      
      fMomentum = fMomentum1+fMomentum2;
      
      Double_t mass  = fMomentum.M();
      Double_t epair = fMomentum.E();
      Float_t ptpair = fMomentum.Pt();
      
      //
      // NLM dependent histograms
      //
      if(fFillAllNLMHistograms)
      {
        if(nMaxima1==nMaxima2)
        {
          if     (nMaxima1==1) fhMassPairLocMax[0]->Fill(epair, mass, GetEventWeight());
          else if(nMaxima1==2) fhMassPairLocMax[1]->Fill(epair, mass, GetEventWeight());
          else                 fhMassPairLocMax[2]->Fill(epair, mass, GetEventWeight());
        }
        else if(nMaxima1==1 || nMaxima2==1)
        {
          if  (nMaxima1==2 || nMaxima2==2) fhMassPairLocMax[3]->Fill(epair, mass, GetEventWeight());
          else                             fhMassPairLocMax[4]->Fill(epair, mass, GetEventWeight());
        }
        else
          fhMassPairLocMax[5]->Fill(epair, mass, GetEventWeight());
        
        // Combinations with SS axis cut and NLM cut
        if(nMaxima1 == 1 && photon2->GetM02() > 0.3)
            fhMassPairLocMax[6]->Fill(epair, mass, GetEventWeight());
          
        if(nMaxima2 == 1 && photon1->GetM02() > 0.3)
            fhMassPairLocMax[6]->Fill(epair, mass, GetEventWeight());
          
        if(nMaxima1 >  1 && photon2->GetM02() < 0.3 && photon2->GetM02() > 0.1 )
            fhMassPairLocMax[7]->Fill(epair, mass, GetEventWeight());
          
        if(nMaxima2 >  1 && photon1->GetM02() < 0.3 && photon1->GetM02() > 0.1 )
            fhMassPairLocMax[7]->Fill(epair, mass, GetEventWeight());
      }
      
      //
      // Skip pairs with too few or too many  NLM
      //
      if((nMaxima1 < fNLMCutMin || nMaxima1 > fNLMCutMax) || (nMaxima2 < fNLMCutMin || nMaxima2 > fNLMCutMax))
      {
        AliDebug(1,Form("NLM out of range: cluster1 %d, cluster2 %d",nMaxima1, nMaxima2));
        continue ;
      }
      
      //
      // Skip pairs with too large/small SS
      //
      if(fM02MaxCutForIM > 0 && fM02MinCutForIM > 0)  // 0.1 < M02 < 0.35
      {
        Float_t m02Ph1 = photon1->GetM02();
        Float_t m02Ph2 = photon2->GetM02();

        if(m02Ph1 > fM02MaxCutForIM || m02Ph1 < fM02MinCutForIM) continue ;
        if(m02Ph2 > fM02MaxCutForIM || m02Ph2 < fM02MinCutForIM) continue ;
      }

      //
      // Mass of all pairs before selection
      //
      fhMass  ->Fill( epair, mass, GetEventWeight());
      fhMassPt->Fill(ptpair, mass, GetEventWeight());
        
      if(photon1->Pt() >= photon2->Pt())
      {
        fhMassPtMaxPair->Fill(photon1->Pt(), mass, GetEventWeight());
        fhMassPtMinPair->Fill(photon2->Pt(), mass, GetEventWeight());
      }
      else
      {
        fhMassPtMaxPair->Fill(photon2->Pt(), mass, GetEventWeight());
        fhMassPtMinPair->Fill(photon1->Pt(), mass, GetEventWeight());
      }
      
      if(IsDataMC() && mcIndex < 2) fhMCMassPt[mcIndex]->Fill(ptpair, mass, GetEventWeight());
      
      if(fSelectPairInIsoCone && fMomentum1.Pt() > fIsoCandMinPt)
      {
        Double_t radius   = GetIsolationCut()->Radius(fMomentum1.Eta(),fMomentum1.Phi(),fMomentum2.Eta(),fMomentum2.Phi());
        
        if(radius < fR) fhMassPtIsoRCut->Fill(ptpair, mass, GetEventWeight());
                    
//        printf("pT pair (%f, %f), opening angle %f, radius %f; fR %1.1f, MinPt1 %2.2f\n",fMomentum1.Pt(),fMomentum2.Pt(),angleOp,radius,fR,fIsoCandMinPt);
      }
      
      //
      // Select good pair (good phi, pt cuts, aperture and invariant mass)
      //
      if(!GetNeutralMesonSelection()->SelectPair(fMomentum1, fMomentum2,GetCalorimeter())) continue;
      
      AliDebug(1,Form("Selected gamma pair: pt %f, phi %f, eta%f",
                      fMomentum.Pt(), fMomentum.Phi()*TMath::RadToDeg(), fMomentum.Eta()));
      
      //
      // Tag both photons as decay if not done before
      // set the corresponding bit for pi0 or eta or "side" case
      //
      Int_t bit1 = photon1->DecayTag();
      if( bit1 < 0 ) bit1 = 0 ;
      if( !GetNeutralMesonSelection()->CheckDecayBit(bit1) )
      {
        AliDebug(1,Form("pT1 %2.2f; bit requested %d; decay bit1: In %d",fMomentum1.Pt(), GetNeutralMesonSelection()->GetDecayBit(), bit1));
        
        GetNeutralMesonSelection()->SetDecayBit(bit1);
        photon1->SetDecayTag(bit1);
        
        AliDebug(1,Form("\t Out %d", bit1));
        
        fhPtDecay->Fill(photon1->Pt(), GetEventWeight());
        
        //
        // Fill some histograms about shower shape
        //
        if(fFillSelectClHisto && cluster1 && GetReader()->GetDataType()!=AliCaloTrackReader::kMC)
        {
          FillSelectedClusterHistograms(cluster1, fMomentum1.Pt(), nMaxima1, photon1->GetTag());
        }
        
        //
        // MC origin, lost pair
        //
        if(IsDataMC())
        {
          Int_t mcIndex1 = GetMCIndex(photon1->GetTag());
          fhMCPtDecay[mcIndex1]->Fill(photon1->Pt(), GetEventWeight());
            
          if(GetMCAnalysisUtils()->CheckTagBit(photon1->GetTag(), AliMCAnalysisUtils::kMCDecayPairLost))
          {
            if     ( mcIndex1 == kmcPi0Decay ) fhMCPtDecayLostPairPi0->Fill(photon1->Pt(), GetEventWeight());
            else if( mcIndex1 == kmcEtaDecay ) fhMCPtDecayLostPairEta->Fill(photon1->Pt(), GetEventWeight());
          }
        }
      }
      
      //
      // Histogram for cluster2 tagged as decay
      //
      Int_t bit2 = photon2->DecayTag();
      if( bit2 < 0 ) bit2 = 0 ;
      if( !GetNeutralMesonSelection()->CheckDecayBit(bit2) )
      {
        AliDebug(1,Form("pT2 %2.2f; bit requested %d; decay bit2: In %d",fMomentum2.Pt(), GetNeutralMesonSelection()->GetDecayBit(), bit2));
        
        GetNeutralMesonSelection()->SetDecayBit(bit2);
        photon2->SetDecayTag(bit2);
        
        AliDebug(1,Form("\t Out %d", bit2));
        
        fhPtDecay->Fill(photon2->Pt(), GetEventWeight());
        
        // Fill some histograms about shower shape
        if(fFillSelectClHisto && GetReader()->GetDataType()!=AliCaloTrackReader::kMC)
        {
          // Get original cluster, to recover some information
          
          Int_t iclus2 = -1;
          cluster2 = FindCluster(clusters,photon2->GetCaloLabel(0),iclus2,iclus+1);
          // start new loop from iclus1+1 to gain some time
          
          if(!cluster2) AliWarning("Second cluster not found");
          else          FillSelectedClusterHistograms(cluster2, fMomentum2.Pt(), nMaxima2, photon2->GetTag());
        }
        
        //
        // MC origin, lost pair
        //
        if(IsDataMC())
        {
          Int_t mcIndex2 = GetMCIndex(photon2->GetTag());
          fhMCPtDecay[mcIndex2]->Fill(photon2->Pt(), GetEventWeight());
            
          if(GetMCAnalysisUtils()->CheckTagBit(photon2->GetTag(), AliMCAnalysisUtils::kMCDecayPairLost))
          {
            if     ( mcIndex2 == kmcPi0Decay ) fhMCPtDecayLostPairPi0->Fill(photon2->Pt(), GetEventWeight());
            else if( mcIndex2 == kmcEtaDecay ) fhMCPtDecayLostPairEta->Fill(photon2->Pt(), GetEventWeight());
          }
        }
      }
      
      //
      // Mass of selected pairs
      //
      fhSelectedMass  ->Fill( epair, mass, GetEventWeight());
      fhSelectedMassPt->Fill(ptpair, mass, GetEventWeight());
      if(IsDataMC()  && mcIndex < 2) fhMCSelectedMassPt[mcIndex]->Fill(ptpair, mass, GetEventWeight());
      
      // Fill histograms to undertand pile-up before other cuts applied
      // Remember to relax time cuts in the reader
      if(cluster1 && IsPileUpAnalysisOn()) FillPileUpHistograms(ptpair,( photon1->GetTime() + photon2->GetTime() ) / 2, cluster1);
      
      //
      // Create AOD for analysis
      //
      AliAODPWG4Particle pi0 = AliAODPWG4Particle(fMomentum);
      
      if     ( (GetNeutralMesonSelection()->GetParticle()).Contains("Pi0") ) pi0.SetIdentifiedParticleType(AliCaloPID::kPi0);
      else if( (GetNeutralMesonSelection()->GetParticle()).Contains("Eta") ) pi0.SetIdentifiedParticleType(AliCaloPID::kEta);
      else
      {
        AliWarning("Particle type declared in AliNeutralMeson not correct, do not add");
        return ;
      }
      pi0.SetDetectorTag(photon1->GetDetectorTag());
      
      // MC
      pi0.SetLabel(label);
      pi0.SetTag(tag);
      
      //Set the indeces of the original caloclusters
      pi0.SetCaloLabel(photon1->GetCaloLabel(0), photon2->GetCaloLabel(0));
      //pi0.SetInputFileIndex(input);
      
      AddAODParticle(pi0);
      
    }//2n photon loop
    
  }//1st photon loop
  
  AliDebug(1,"End fill AODs");
}

//__________________________________________________
/// Read photon list from AOD, produced in class AliAnaPhoton and AliGammaConversion
/// Check if 2 photons from each list have the mass of the pi0/eta, depending on the
/// settings of AliNeutralMesonSelection. If the pair is selected create the
/// particle AOD object with the pair and fill the AOD branch.
/// Fill also some control histograms.
//__________________________________________________
void  AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS()
{
  // Check calorimeter input
  if(!GetInputAODBranch())
  {
    AliFatal(Form("No input calo photons in AOD branch with name < %s > , STOP",GetInputAODName().Data()));
  }
  
  // Get the array with conversion photons
  TClonesArray * inputAODGammaConv = (TClonesArray *) GetReader()->GetOutputEvent()->FindListObject(fInputAODGammaConvName);
  if(!inputAODGammaConv)
  {
    inputAODGammaConv = (TClonesArray *) GetReader()->GetInputEvent()->FindListObject(fInputAODGammaConvName);
    
    if(!inputAODGammaConv)
    {
      AliFatal(Form("No input gamma conversions in AOD branch with name < %s >",fInputAODGammaConvName.Data()));
      return; // coverity
    }
  }
  
  // Get shower shape information of clusters
  TObjArray *clusters = 0;
  if( fFillSelectClHisto || IsPileUpAnalysisOn() )
  {
    if     (GetCalorimeter() == kEMCAL) clusters = GetEMCALClusters() ;
    else if(GetCalorimeter() == kPHOS ) clusters = GetPHOSClusters () ;
  }
  
  AliVCluster * cluster = 0;
  
  Int_t nCTS  = inputAODGammaConv->GetEntriesFast();
  Int_t nCalo = GetInputAODBranch()->GetEntriesFast();
  if(nCTS<=0 || nCalo <=0)
  {
    AliDebug(1,Form("nCalo %d, nCTS %d, cannot loop",nCalo,nCTS));
    return;
  }
  
  AliDebug(1,Form("Number of conversion photons %d and number of clusters %d",nCTS,nCalo));
  
  // Do the loop, first calo, second CTS
  for(Int_t iphoton = 0; iphoton < GetInputAODBranch()->GetEntriesFast(); iphoton++)
  {
    AliAODPWG4Particle * photon1 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(iphoton));
    fMomentum1 = *(photon1->Momentum());
    
    // Do analysis only when one of the decays is isolated
    // Run AliAnaParticleIsolation before
    if(fSelectIsolatedDecay)
    {
      Bool_t isolated1 = ((AliAODPWG4ParticleCorrelation*) photon1)->IsIsolated();
      if(!isolated1) continue;
    }
    
    // Get original cluster, to recover some information
    Int_t iclus = -1;
    if( fFillSelectClHisto || IsPileUpAnalysisOn() )
    {
      cluster = FindCluster(clusters,photon1->GetCaloLabel(0),iclus);
      
      if(!cluster) AliWarning("Cluster not found");
    }
    
    // Skip cluster with too large/small SS
    if(fM02MaxCutForIM > 0 && fM02MinCutForIM > 0)  // 0.1 < M02 < 0.35
    {
      Float_t m02 = photon1->GetM02();
      
      if(m02 > fM02MaxCutForIM || m02 < fM02MinCutForIM) continue ;
    }
    
    for(Int_t jphoton = 0; jphoton < nCTS; jphoton++)
    {
      AliAODPWG4Particle * photon2 =  (AliAODPWG4Particle*) (inputAODGammaConv->At(jphoton));
      
      Int_t evtIndex = 0;
      if(GetMixedEvent())
      {
        evtIndex = GetMixedEvent()->EventIndexForCaloCluster(photon2->GetCaloLabel(0)) ;
        if(TMath::Abs(GetVertex(evtIndex)[2]) > GetZvertexCut()) continue ;  //vertex cut
      }
      
      fMomentum2 = *(photon2->Momentum());
      
      fMomentum = fMomentum1+fMomentum2;
      
      Double_t mass  = fMomentum.M();
      Double_t epair = fMomentum.E();
      Float_t ptpair = fMomentum.Pt();
      
      Int_t nMaxima = photon1->GetNLM();
      if(fFillAllNLMHistograms)
      {
        if     (nMaxima==1) fhMassPairLocMax[0]->Fill(epair, mass, GetEventWeight());
        else if(nMaxima==2) fhMassPairLocMax[1]->Fill(epair, mass, GetEventWeight());
        else                fhMassPairLocMax[2]->Fill(epair, mass, GetEventWeight());
      }
      
      if(nMaxima < fNLMCutMin || nMaxima > fNLMCutMax)
      {
        AliDebug(1,Form("NLM %d out of range",nMaxima));
        continue ;
      }
            
      // Play with the MC stack if available
      Int_t mcIndex = kmcHadron;
      Int_t tag     = 0;
      Int_t label   =-1;
      if(IsDataMC())
      {
        Int_t	label2 = photon2->GetLabel();
        if(label2 >= 0 )photon2->SetTag(GetMCAnalysisUtils()->CheckOrigin(label2, GetReader(),kCTS));
        
        HasPairSameMCMother(photon1->GetLabel(), photon2->GetLabel(),
                            photon1->GetTag()  , photon2->GetTag(),
                            label, tag) ;
        mcIndex = GetMCIndex(tag);
      }
      
      // Mass of selected pairs
      fhMass  ->Fill( epair, mass, GetEventWeight());
      fhMassPt->Fill(ptpair, mass, GetEventWeight());
  
      if(photon1->Pt() >= photon2->Pt())
      {
        fhMassPtMaxPair->Fill(photon1->Pt(), mass, GetEventWeight());
        fhMassPtMinPair->Fill(photon2->Pt(), mass, GetEventWeight());
      }
      else
      {
        fhMassPtMaxPair->Fill(photon2->Pt(), mass, GetEventWeight());
        fhMassPtMinPair->Fill(photon1->Pt(), mass, GetEventWeight());
      }
      
      if(IsDataMC() && mcIndex < 2 ) fhMCMassPt[mcIndex]->Fill(ptpair, mass, GetEventWeight());
      
      //
      // Select good pair (good phi, pt cuts, aperture and invariant mass)
      //
      if(!GetNeutralMesonSelection()->SelectPair(fMomentum1, fMomentum2,GetCalorimeter())) continue ;
      
      AliDebug(1,Form("Selected gamma pair: pt %f, phi %f, eta%f",fMomentum.Pt(), fMomentum.Phi()*TMath::RadToDeg(), fMomentum.Eta()));
      
      //
      // Tag both photons as decay if not done before
      // set the corresponding bit for pi0 or eta or "side" case
      //
      Int_t bit1 = photon1->DecayTag();
      if( bit1 < 0 ) bit1 = 0 ;
      if( !GetNeutralMesonSelection()->CheckDecayBit(bit1) )
      {
        GetNeutralMesonSelection()->SetDecayBit(bit1);
        photon1->SetDecayTag(bit1);
        
        fhPtDecay->Fill(photon1->Pt(), GetEventWeight());
        
        // Fill some histograms about shower shape
        if(fFillSelectClHisto && cluster && GetReader()->GetDataType()!=AliCaloTrackReader::kMC)
          FillSelectedClusterHistograms(cluster, fMomentum1.Pt(), nMaxima, photon1->GetTag());
        
        if(IsDataMC())
        {
          Int_t mcIndex1 = GetMCIndex(photon1->GetTag());
          fhMCPtDecay[mcIndex1]->Fill(photon1->Pt(), GetEventWeight());
            
          if(GetMCAnalysisUtils()->CheckTagBit(photon1->GetTag(), AliMCAnalysisUtils::kMCDecayPairLost))
          {
            if     ( mcIndex1 == kmcPi0Decay ) fhMCPtDecayLostPairPi0->Fill(photon1->Pt(), GetEventWeight());
            else if( mcIndex1 == kmcEtaDecay ) fhMCPtDecayLostPairEta->Fill(photon1->Pt(), GetEventWeight());
          }
        }
      }
      
      Int_t bit2 = photon2->DecayTag();
      if( bit2 < 0 ) bit2 = 0 ;
      if( !GetNeutralMesonSelection()->CheckDecayBit(bit2) )
      {
        GetNeutralMesonSelection()->SetDecayBit(bit2);
        photon2->SetDecayTag(bit2);
      }
      
      //
      // Mass of selected pairs
      //
      fhSelectedMass  ->Fill( epair, mass, GetEventWeight());
      fhSelectedMassPt->Fill(ptpair, mass, GetEventWeight());
      if(IsDataMC()  && mcIndex < 2) fhMCSelectedMassPt[mcIndex]->Fill(ptpair, mass, GetEventWeight());
      
      // Fill histograms to undertand pile-up before other cuts applied
      // Remember to relax time cuts in the reader
      if( cluster && IsPileUpAnalysisOn() ) FillPileUpHistograms(fMomentum.Pt(),photon1->GetTime(),cluster);
      
      //
      // Create AOD for analysis
      //
      AliAODPWG4Particle pi0 = AliAODPWG4Particle(fMomentum);
      
      if     ( (GetNeutralMesonSelection()->GetParticle()).Contains("Pi0") ) pi0.SetIdentifiedParticleType(AliCaloPID::kPi0);
      else if( (GetNeutralMesonSelection()->GetParticle()).Contains("Eta") ) pi0.SetIdentifiedParticleType(AliCaloPID::kEta);
      else
      {
        AliWarning("Particle type declared in AliNeutralMeson not correct, do not add");
        return ;
      }
      
      pi0.SetDetectorTag(photon1->GetDetectorTag());
      
      // MC
      pi0.SetLabel(label);
      pi0.SetTag(tag);
      
      // Set the indeces of the original tracks or caloclusters
      pi0.SetCaloLabel (photon1->GetCaloLabel(0) , -1);
      pi0.SetTrackLabel(photon2->GetTrackLabel(0), photon2->GetTrackLabel(1));
      //pi0.SetInputFileIndex(input);
      
      AddAODParticle(pi0);
      
    }//2n photon loop
    
  }//1st photon loop
  
  AliDebug(1,"End fill AODs");
}

//_________________________________________________
/// Read the EMCal cluster list (not implemented for PHOS)
/// select the clusters and tag them as pi0/eta (or conversion) depending on:
///  * Track matching rejection.
///  * Fiducial acceptance, distance to bad channels.
///  * Minimum number of cells.
///  * Number of local maxima range.
///  * Shower shape.
///  * sub-cluster split mass and asymmetry.
/// The cluster splitting and selection is done
/// via the utils in AliCaloPID and AliCalorimeterUtils.
/// Fill also some control histograms.
//_________________________________________________
void  AliAnaPi0EbE::MakeShowerShapeIdentification()
{
  TObjArray * pl        = 0x0;
  AliVCaloCells * cells = 0x0;
  
  // Select the Calorimeter of the photon
  if      (GetCalorimeter() == kEMCAL )
  {
    pl    = GetEMCALClusters();
    cells = GetEMCALCells();
  }
  else if (GetCalorimeter() == kPHOS)
  {
    AliFatal("kSSCalo case not implememted for PHOS");
    return; // for coverity
    
    //pl    = GetPHOSClusters();
    //cells = GetPHOSCells();
  }
  
  if(!pl)
  {
    AliInfo(Form("TObjArray with %s clusters is NULL!",GetCalorimeterString().Data()));
    return;
  }
  
  for(Int_t icalo = 0; icalo < pl->GetEntriesFast(); icalo++)
  {
    AliVCluster * calo = (AliVCluster*) (pl->At(icalo));
    
    Int_t evtIndex = 0 ;
    if (GetMixedEvent())
    {
      evtIndex=GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ;
      if(TMath::Abs(GetVertex(evtIndex)[2]) > GetZvertexCut()) continue ;  //vertex cut
    }
    
    // Get Momentum vector,
    Double_t vertex[]={0,0,0};
    if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    {
      calo->GetMomentum(fMomentum,GetVertex(evtIndex)) ;
    } // Assume that come from vertex in straight line
    else
    {
      calo->GetMomentum(fMomentum,vertex) ;
    }
    
    // If too small or big pt, skip it
    if(fMomentum.E() < GetMinEnergy() || fMomentum.E() > GetMaxEnergy() ) continue ;
    
    // Check acceptance selection
    if(IsFiducialCutOn())
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),GetCalorimeter()) ;
      if(! in ) continue ;
    }
    
    AliDebug(1,Form("Min pt cut and fiducial cut passed: pt %3.2f, phi %2.2f, eta %1.2f",fMomentum.Pt(),fMomentum.Phi(),fMomentum.Eta()));
    
    // Play with the MC stack if available
    // Check origin of the candidates
    Int_t tag	= 0 ;
    if(IsDataMC())
    {
      tag = GetMCAnalysisUtils()->CheckOrigin(calo->GetLabels(),calo->GetNLabels(),GetReader(),GetCalorimeter());
      //GetMCAnalysisUtils()->CheckMultipleOrigin(calo->GetLabels(),calo->GetNLabels(), GetReader(), aodpi0.GetInputFileIndex(), tag);
      AliDebug(1,Form("Origin of candidate %d",tag));
    }
    
    //Int_t nMaxima = GetCaloUtils()->GetNumberOfLocalMaxima(calo, cells); // NLM
    
    // Check Distance to Bad channel, set bit.
    Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
    if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
    if(distBad < fMinDist)
    {
      // In bad channel (PHOS cristal size 2.2x2.2 cm)
      // FillRejectedClusterHistograms(tag,nMaxima);
      continue ;
    }
    
    AliDebug(1,Form("Bad channel cut passed %4.2f",distBad));
    
    // If too low number of cells, skip it
    if ( calo->GetNCells() < GetCaloPID()->GetClusterSplittingMinNCells())
    {
      //FillRejectedClusterHistograms(tag,nMaxima);
      continue ;
    }
    
    AliDebug(1,Form("N cells cut passed %d > %d",calo->GetNCells(), GetCaloPID()->GetClusterSplittingMinNCells()));
    
    //.......................................
    // TOF cut, BE CAREFUL WITH THIS CUT
    Double_t tof = calo->GetTOF()*1e9;
    if(tof < fTimeCutMin || tof > fTimeCutMax)
    {
      //FillRejectedClusterHistograms(tag,nMaxima);
      continue ;
    }
    
    // Check PID
    // PID selection or bit setting
    Int_t    nMaxima  = 0;
    Double_t mass     = 0, angle    = 0;
    Int_t    absId1   =-1, absId2   =-1;
    Float_t  distbad1 =-1, distbad2 =-1;
    Bool_t   fidcut1  = 0, fidcut2  = 0;
    
    Int_t idPartType = GetCaloPID()->GetIdentifiedParticleTypeFromClusterSplitting(calo,cells,GetCaloUtils(),
                                                                                   GetVertex(evtIndex),nMaxima,
                                                                                   mass,angle,
                                                                                   fMomentum1,fMomentum2,
                                                                                   absId1,absId2,
                                                                                   distbad1,distbad2,
                                                                                   fidcut1,fidcut2) ;
    
    
    AliDebug(1,Form("PDG of identified particle %d",idPartType));
    
    // Skip events where one of the new clusters (lowest energy) is close to an EMCal border or a bad channel
    if( (fCheckSplitDistToBad) &&
       (!fidcut2 || !fidcut1 || distbad1 < fMinDist || distbad2 < fMinDist))
    {
      AliDebug(1,Form("Dist to bad channel cl %f, cl1 %f, cl2 %f; fid cl1 %d, cl2 %d",
                      calo->GetDistanceToBadChannel(),distbad1,distbad2, fidcut1,fidcut2));
      
      //FillRejectedClusterHistograms(tag,nMaxima);
      continue ;
    }
    
    // Skip events with too few or too many  NLM
    if(nMaxima < fNLMCutMin || nMaxima > fNLMCutMax)
    {
      //FillRejectedClusterHistograms(tag,nMaxima);
      continue ;
    }
    
    AliDebug(1,Form("NLM %d accepted",nMaxima));
    
    // Skip matched clusters with tracks
    if(fRejectTrackMatch && IsTrackMatched(calo, GetReader()->GetInputEvent()))
    {
      FillRejectedClusterHistograms(tag,nMaxima);
      continue ;
    }
    
    Float_t l0 = calo->GetM02();
    Float_t l1 = calo->GetM20();
    Float_t e1 = fMomentum1.Energy();
    Float_t e2 = fMomentum2.Energy();
    fMomentum12 = fMomentum1+fMomentum2;
    Float_t ptSplit = fMomentum12.Pt();
    Float_t  eSplit = e1+e2;
    
    // Mass of all clusters
    fhMass  ->Fill(fMomentum.E() , mass, GetEventWeight());
    fhMassPt->Fill(fMomentum.Pt(), mass, GetEventWeight());
      
    if(fMomentum1.Pt() >= fMomentum2.Pt())
    {
      fhMassPtMaxPair->Fill(fMomentum1.Pt(), mass, GetEventWeight());
      fhMassPtMinPair->Fill(fMomentum2.Pt(), mass, GetEventWeight());
    }
    else
    {
      fhMassPtMaxPair->Fill(fMomentum2.Pt(), mass, GetEventWeight());
      fhMassPtMinPair->Fill(fMomentum1.Pt(), mass, GetEventWeight());
    }
      
    fhMassSplitPt->Fill(ptSplit, mass, GetEventWeight());
    fhPtLambda0NoSplitCut->Fill(fMomentum.Pt(), l0, GetEventWeight());
    
    // Asymmetry of all clusters
    Float_t asy =-10;
    if(e1+e2 > 0) asy = (e1-e2) / (e1+e2);
    fhAsymmetry->Fill(fMomentum.E(), asy, GetEventWeight());
    
    // Divide NLM in 3 cases, 1 local maxima, 2 local maxima, more than 2 local maxima
    Int_t indexMax = -1;
    if     (nMaxima==1) indexMax = 0 ;
    else if(nMaxima==2) indexMax = 1 ;
    else                indexMax = 2 ;
      
    fhMassPtLocMax[indexMax]->Fill(fMomentum.Pt(), mass, GetEventWeight());
    
    Int_t   mcIndex    = -1;
    Int_t   noverlaps  =  0;
    Float_t ptprim     =  0;
    Int_t   mesonLabel = -1;
    
    if(IsDataMC())
    {
      mcIndex = GetMCIndex(tag);
      
      Bool_t ok      = kFALSE;
      Int_t  mcLabel = calo->GetLabel();
      
      fPrimaryMom = GetMCAnalysisUtils()->GetMother(mcLabel,GetReader(),ok);
      
      
      if(mcIndex == kmcPi0 || mcIndex == kmcEta)
      {
        if(mcIndex == kmcPi0)
        {
          fGrandMotherMom = GetMCAnalysisUtils()->GetMotherWithPDG(mcLabel,111,GetReader(),ok,mesonLabel);
          if(fGrandMotherMom.E() > 0 && ok) ptprim =  fGrandMotherMom.Pt();
        }
        else
        {
          fGrandMotherMom = GetMCAnalysisUtils()->GetMotherWithPDG(mcLabel,221,GetReader(),ok,mesonLabel);
          if(fGrandMotherMom.E() > 0 && ok) ptprim = fGrandMotherMom.Pt();
        }
      }
      
      const UInt_t nlabels = calo->GetNLabels();
      Int_t overpdg[nlabels];
      Int_t overlab[nlabels];
      noverlaps = GetMCAnalysisUtils()->GetNOverlaps(calo->GetLabels(), nlabels, tag, mesonLabel, GetReader(), overpdg, overlab);
      
      fhMCMassPt     [mcIndex]->Fill(fMomentum.Pt(), mass, GetEventWeight());
      fhMCMassSplitPt[mcIndex]->Fill(ptSplit       , mass, GetEventWeight());
        
      if(mcIndex==kmcPi0)
      {
        fhMCPi0PtRecoPtPrim                     ->Fill(fMomentum.Pt(), ptprim, GetEventWeight());
        fhMCPi0SplitPtRecoPtPrim                ->Fill(ptSplit       , ptprim, GetEventWeight());
        fhMCPi0PtRecoPtPrimLocMax     [indexMax]->Fill(fMomentum.Pt(), ptprim, GetEventWeight());
        fhMCPi0SplitPtRecoPtPrimLocMax[indexMax]->Fill(ptSplit       , ptprim, GetEventWeight());
        
      }
      else if(mcIndex==kmcEta)
      {
        fhMCEtaPtRecoPtPrim                     ->Fill(fMomentum.Pt(), ptprim, GetEventWeight());
        fhMCEtaSplitPtRecoPtPrim                ->Fill(ptSplit       , ptprim, GetEventWeight());
        fhMCEtaPtRecoPtPrimLocMax     [indexMax]->Fill(fMomentum.Pt(), ptprim, GetEventWeight());
        fhMCEtaSplitPtRecoPtPrimLocMax[indexMax]->Fill(ptSplit       , ptprim, GetEventWeight());
      }
      
      if(noverlaps==0)
      {
        if(mcIndex==kmcPi0)
        {
          fhMCPi0PtRecoPtPrimNoOverlap     ->Fill(fMomentum.Pt(), ptprim, GetEventWeight());
          fhMCPi0SplitPtRecoPtPrimNoOverlap->Fill(ptSplit       , ptprim, GetEventWeight());
        }
        else if(mcIndex==kmcEta)
        {
          fhMCEtaPtRecoPtPrimNoOverlap     ->Fill(fMomentum.Pt(), ptprim, GetEventWeight());
          fhMCEtaSplitPtRecoPtPrimNoOverlap->Fill(ptSplit       , ptprim, GetEventWeight());
        }
        
        fhMassNoOverlap       ->Fill(fMomentum.E() , mass, GetEventWeight());
        fhMassPtNoOverlap     ->Fill(fMomentum.Pt(), mass, GetEventWeight());
        fhMassSplitPtNoOverlap->Fill(ptSplit       , mass, GetEventWeight());
        
        fhMCMassPtNoOverlap     [mcIndex]->Fill(fMomentum.Pt(), mass, GetEventWeight());
        fhMCMassSplitPtNoOverlap[mcIndex]->Fill(ptSplit       , mass, GetEventWeight());
      }
      
      fhMCPtAsymmetry[mcIndex]->Fill(fMomentum.Pt(), asy, GetEventWeight());
    }
    
    // If cluster does not pass pid, not pi0/eta, skip it.
    if     (GetOutputAODName().Contains("Pi0") && idPartType != AliCaloPID::kPi0)
    {
      AliDebug(1,"Cluster is not Pi0");
      FillRejectedClusterHistograms(tag,nMaxima);
      continue ;
    }
    
    else if(GetOutputAODName().Contains("Eta") && idPartType != AliCaloPID::kEta)
    {
      AliDebug(1,"Cluster is not Eta");
      FillRejectedClusterHistograms(tag,nMaxima);
      continue ;
    }
    
    AliDebug(1,Form("Pi0/Eta selection cuts passed: pT %3.2f, pdg %d",fMomentum.Pt(), idPartType));
    
    // Mass and asymmetry of selected pairs
    fhSelectedAsymmetry  ->Fill(fMomentum.E() , asy , GetEventWeight());
    fhSelectedMass       ->Fill(fMomentum.E() , mass, GetEventWeight());
    fhSelectedMassPt     ->Fill(fMomentum.Pt(), mass, GetEventWeight());
    fhSelectedMassSplitPt->Fill(ptSplit       , mass, GetEventWeight());
      
    fhSelectedMassPtLocMax[indexMax]->Fill(fMomentum.Pt(), mass, GetEventWeight());
    
    Int_t   nSM  = GetModuleNumber(calo);
    if(nSM < fNModules && nSM >=0 && fFillAllNLMHistograms)
    {
      fhSelectedMassPtLocMaxSM   [indexMax][nSM]->Fill(fMomentum.Pt(), mass, GetEventWeight());
      fhSelectedLambda0PtLocMaxSM[indexMax][nSM]->Fill(fMomentum.Pt(), l0  , GetEventWeight());
    }
    
    if(IsDataMC())
    {
      if(mcIndex==kmcPi0)
      {
        fhMCPi0SelectedPtRecoPtPrim                     ->Fill(fMomentum.Pt(), ptprim, GetEventWeight());
        fhMCPi0SelectedSplitPtRecoPtPrim                ->Fill(ptSplit       , ptprim, GetEventWeight());
        fhMCPi0SelectedPtRecoPtPrimLocMax     [indexMax]->Fill(fMomentum.Pt(), ptprim, GetEventWeight());
        fhMCPi0SelectedSplitPtRecoPtPrimLocMax[indexMax]->Fill(ptSplit       , ptprim, GetEventWeight());
      }
      else if(mcIndex==kmcEta)
      {
        fhMCEtaSelectedPtRecoPtPrim                     ->Fill(fMomentum.Pt(), ptprim, GetEventWeight());
        fhMCEtaSelectedSplitPtRecoPtPrim                ->Fill(ptSplit       , ptprim, GetEventWeight());
        fhMCEtaSelectedPtRecoPtPrimLocMax     [indexMax]->Fill(fMomentum.Pt(), ptprim, GetEventWeight());
        fhMCEtaSelectedSplitPtRecoPtPrimLocMax[indexMax]->Fill(ptSplit       , ptprim, GetEventWeight());
      }
      
      if(noverlaps==0)
      {
        fhSelectedMassNoOverlap       ->Fill(fMomentum.E() , mass, GetEventWeight());
        fhSelectedMassPtNoOverlap     ->Fill(fMomentum.Pt(), mass, GetEventWeight());
        fhSelectedMassSplitPtNoOverlap->Fill(ptSplit       , mass, GetEventWeight());
        
        if(mcIndex==kmcPi0)
        {
          fhMCPi0SelectedPtRecoPtPrimNoOverlap     ->Fill(fMomentum.Pt(), ptprim, GetEventWeight());
          fhMCPi0SelectedSplitPtRecoPtPrimNoOverlap->Fill(ptSplit       , ptprim, GetEventWeight());
        }
        else if(mcIndex==kmcEta)
        {
          fhMCEtaSelectedPtRecoPtPrimNoOverlap     ->Fill(fMomentum.Pt(), ptprim, GetEventWeight());
          fhMCEtaSelectedSplitPtRecoPtPrimNoOverlap->Fill(ptSplit       , ptprim, GetEventWeight());
        }
      }
    }
    
    Float_t phi = fMomentum.Phi();
    if(phi<0) phi+=TMath::TwoPi();
      
    fhSplitE        ->Fill( eSplit,                  GetEventWeight());
    fhSplitPt       ->Fill(ptSplit,                  GetEventWeight());
    fhSplitPtPhi    ->Fill(ptSplit, phi            , GetEventWeight());
    fhSplitPtEta    ->Fill(ptSplit, fMomentum.Eta(), GetEventWeight());
    fhNLocMaxSplitPt->Fill(ptSplit, nMaxima        , GetEventWeight());
    
    // Check split-clusters with good time window difference
    Double_t tof1  = cells->GetCellTime(absId1);
    GetCaloUtils()->RecalibrateCellTime(tof1, GetCalorimeter(), absId1,GetReader()->GetInputEvent()->GetBunchCrossNumber());
    tof1*=1.e9;
    
    Double_t tof2  = cells->GetCellTime(absId2);
    GetCaloUtils()->RecalibrateCellTime(tof2, GetCalorimeter(), absId2,GetReader()->GetInputEvent()->GetBunchCrossNumber());
    tof2*=1.e9;
    
    Double_t t12diff = tof1-tof2;
    fhEPairDiffTime->Fill(e1+e2, t12diff, GetEventWeight());
    
    if(IsDataMC())
    {
      fhMCSplitE        [mcIndex]->Fill( eSplit,                         GetEventWeight());
      fhMCSplitPt       [mcIndex]->Fill(ptSplit,                         GetEventWeight());
      fhMCSplitPtPhi    [mcIndex]->Fill(ptSplit       , phi            , GetEventWeight());
      fhMCSplitPtEta    [mcIndex]->Fill(ptSplit       , fMomentum.Eta(), GetEventWeight());
      fhMCNLocMaxSplitPt[mcIndex]->Fill(ptSplit       , nMaxima        , GetEventWeight());
      fhMCNLocMaxPt     [mcIndex]->Fill(fMomentum.Pt(), nMaxima        , GetEventWeight());
      
      fhMCSelectedMassPt      [mcIndex]->Fill(fMomentum.Pt(), mass, GetEventWeight());
      fhMCSelectedMassSplitPt [mcIndex]->Fill(ptSplit       , mass, GetEventWeight());
      fhMCSelectedMassPtLocMax[mcIndex][indexMax]->Fill(fMomentum.Pt(), mass, GetEventWeight());
      
      if(noverlaps==0)
      {
        fhMCSelectedMassPtNoOverlap     [mcIndex]->Fill(fMomentum.Pt(), mass, GetEventWeight());
        fhMCSelectedMassSplitPtNoOverlap[mcIndex]->Fill(ptSplit       , mass, GetEventWeight());
      }
    }
    
    // Remove clusters with NLM=x depeding on a minimim energy cut
    if(nMaxima == 1 && fNLMECutMin[0] > fMomentum.E()) continue;
    if(nMaxima == 2 && fNLMECutMin[1] > fMomentum.E()) continue;
    if(nMaxima >  2 && fNLMECutMin[2] > fMomentum.E()) continue;
    
    // Fill some histograms about shower shape
    if(fFillSelectClHisto && GetReader()->GetDataType()!=AliCaloTrackReader::kMC)
    {
      FillSelectedClusterHistograms(calo, fMomentum.Pt(), nMaxima, tag, asy);
    }
    
    // Fill histograms to undertand pile-up before other cuts applied
    // Remember to relax time cuts in the reader
    Double_t tofcluster   = calo->GetTOF()*1e9;
    
    if ( IsPileUpAnalysisOn() ) FillPileUpHistograms(fMomentum.Pt(),tofcluster,calo);
    
    if(fFillEMCALBCHistograms && GetCalorimeter()==kEMCAL)
      FillEMCALBCHistograms(fMomentum.E(), fMomentum.Eta(), fMomentum.Phi(), tofcluster);
    
    //------------------------
    // Create AOD for analysis
    
    AliAODPWG4Particle aodpi0 = AliAODPWG4Particle(fMomentum);
    aodpi0.SetLabel(mesonLabel);
    
    // Set the indeces of the original caloclusters
    aodpi0.SetCaloLabel(calo->GetID(),-1);
    aodpi0.SetDetectorTag(GetCalorimeter());
    
    if     (distBad > fMinDist3) aodpi0.SetDistToBad(2) ;
    else if(distBad > fMinDist2) aodpi0.SetDistToBad(1) ;
    else                         aodpi0.SetDistToBad(0) ;
    
    // Check if cluster is pi0 via cluster splitting
    aodpi0.SetIdentifiedParticleType(idPartType);
    
    aodpi0.SetM02(l0);
    aodpi0.SetM20(l1);
    aodpi0.SetNLM(nMaxima);
    aodpi0.SetTime(tofcluster);
    aodpi0.SetNCells(calo->GetNCells());
    aodpi0.SetSModNumber(nSM);
    
    aodpi0.SetTag(tag);
    
    // Add AOD with pi0 object to aod branch
    AddAODParticle(aodpi0);
    
  } // loop
  
  AliDebug(1,"End fill AODs");
}

//______________________________________________
/// Fill histograms for selected pairs/single clusters
/// Mostly common histograms for the 3 implemented methods.
//______________________________________________
void  AliAnaPi0EbE::MakeAnalysisFillHistograms()
{
  if(!GetOutputAODBranch())
  {
    AliFatal(Form("No output pi0 in AOD branch with name < %s >,STOP",GetOutputAODName().Data()));
    return;
  }
  
  // Loop on stored AOD pi0
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  
  AliDebug(1,Form("AOD branch entries %d", naod));
  
  Float_t cen = GetEventCentrality();
  Float_t ep  = GetEventPlaneAngle();
  
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliAODPWG4Particle* pi0 =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = pi0->GetIdentifiedParticleType();
    
    if( ( pdg != AliCaloPID::kPi0 && pdg != AliCaloPID::kEta ) ) continue;
    
    // Fill pi0 histograms
    Float_t ener  = pi0->E();
    Float_t pt    = pi0->Pt();
    Float_t phi   = pi0->Phi();
    if(phi < 0) phi+=TMath::TwoPi();
    Float_t eta = pi0->Eta();
    
    fhPt     ->Fill(pt  , GetEventWeight());
    fhE      ->Fill(ener, GetEventWeight());
    
    fhPtEta  ->Fill(pt  , eta, GetEventWeight());
    fhPtPhi  ->Fill(pt  , phi, GetEventWeight());
    fhEtaPhi ->Fill(eta , phi, GetEventWeight());
    
    if(IsHighMultiplicityAnalysisOn())
    {
      fhPtCentrality ->Fill(pt, cen, GetEventWeight()) ;
      fhPtEventPlane ->Fill(pt, ep , GetEventWeight()) ;
    }
    
    if(IsDataMC())
    {
      Int_t tag     = pi0->GetTag();
      Int_t label   = pi0->GetLabel();
      Int_t mcIndex = GetMCIndex(tag);
      
      if(fAnaType != kSSCalo && mcIndex > 1) continue;
      
      fhMCE    [mcIndex] ->Fill(ener,      GetEventWeight());
      fhMCPt   [mcIndex] ->Fill(pt  ,      GetEventWeight());
      fhMCPtPhi[mcIndex] ->Fill(pt  , phi, GetEventWeight());
      fhMCPtEta[mcIndex] ->Fill(pt  , eta, GetEventWeight());
      
      if ( IsHighMultiplicityAnalysisOn() ) fhMCPtCentrality[mcIndex]->Fill(pt, cen, GetEventWeight());
      
      if((mcIndex==kmcPi0Decay || mcIndex==kmcEtaDecay ||
          mcIndex==kmcPi0      || mcIndex==kmcEta         ) &&
         fAnaType==kSSCalo)
      {
        Float_t efracMC   = 0;
        Int_t   momlabel  = -1;
        Bool_t  ok        = kFALSE;
        
        fPrimaryMom = GetMCAnalysisUtils()->GetMother(label,GetReader(),ok);
        if(!ok) continue;
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0))
        {
          fGrandMotherMom = GetMCAnalysisUtils()->GetMotherWithPDG(label,111,GetReader(),ok,momlabel);
          if(fGrandMotherMom.E() > 0 && ok)
          {
            efracMC =  fGrandMotherMom.E()/ener;
            fhMCPi0PtGenRecoFraction ->Fill(pt, efracMC, GetEventWeight());
          }
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))
        {
          fhMCPi0DecayPt->Fill(pt, GetEventWeight());
          fGrandMotherMom = GetMCAnalysisUtils()->GetMotherWithPDG(label,111,GetReader(),ok,momlabel);
            
          if(fGrandMotherMom.E() > 0 && ok)
          {
            efracMC =  fPrimaryMom.E()/fGrandMotherMom.E();
            fhMCPi0DecayPtFraction ->Fill(pt, efracMC, GetEventWeight());
          }
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta))
        {
          fGrandMotherMom = GetMCAnalysisUtils()->GetMotherWithPDG(label,221,GetReader(),ok,momlabel);
          if(fGrandMotherMom.E() > 0 && ok)
          {
            efracMC =  fGrandMotherMom.E()/ener;
            fhMCEtaPtGenRecoFraction ->Fill(pt, efracMC, GetEventWeight());
          }
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay))
        {
          fhMCEtaDecayPt->Fill(pt, GetEventWeight());
          fGrandMotherMom = GetMCAnalysisUtils()->GetMotherWithPDG(label,221,GetReader(),ok,momlabel);
            
          if(fGrandMotherMom.E() > 0 && ok)
          {
            efracMC =  fPrimaryMom.E()/fGrandMotherMom.E();
            fhMCEtaDecayPtFraction ->Fill(pt, efracMC, GetEventWeight());
          }
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
        {
          fhMCOtherDecayPt->Fill(pt, GetEventWeight());
        }
      }
      
      if( mcIndex==kmcPi0 || mcIndex==kmcEta )
      {
        Float_t prodR     = -1;
        Int_t   momindex  = -1;
        Int_t   mompdg    = -1;
        Int_t   momstatus = -1;
        Int_t status      = -1;
        if(GetReader()->ReadStack())
        {
          TParticle* ancestor = GetMCStack()->Particle(label);
          status = ancestor->GetStatusCode();
          momindex  = ancestor->GetFirstMother();
          if(momindex < 0) return;
          TParticle* mother = GetMCStack()->Particle(momindex);
          mompdg    = TMath::Abs(mother->GetPdgCode());
          momstatus = mother->GetStatusCode();
          prodR = mother->R();
        }
        else
        {
          TClonesArray * mcparticles = GetReader()->GetAODMCParticles();
          AliAODMCParticle* ancestor = (AliAODMCParticle *) mcparticles->At(label);
          status    = ancestor->GetStatus();
          momindex  = ancestor->GetMother();
          
          if(momindex < 0) return;
            
          AliAODMCParticle* mother = (AliAODMCParticle *) mcparticles->At(momindex);
          mompdg    = TMath::Abs(mother->GetPdgCode());
          momstatus = mother->GetStatus();
          prodR     = TMath::Sqrt(mother->Xv()*mother->Xv()+mother->Yv()*mother->Yv());
        }
        
        if( mcIndex==kmcPi0 )
        {
          fhMCPi0ProdVertex->Fill(pt, prodR , GetEventWeight());
          fhMCPi0PtStatus  ->Fill(pt, status, GetEventWeight());
          
          if     (momstatus  == 21) fhMCPi0PtOrigin->Fill(pt, 0.5, GetEventWeight());//parton
          else if(mompdg     < 22 ) fhMCPi0PtOrigin->Fill(pt, 1.5, GetEventWeight());//quark
          else if(mompdg     > 2100  && mompdg   < 2210) fhMCPi0PtOrigin->Fill(pt, 2.5, GetEventWeight());// resonances
          else if(mompdg    == 221) fhMCPi0PtOrigin->Fill(pt, 8.5, GetEventWeight());//eta
          else if(mompdg    == 331) fhMCPi0PtOrigin->Fill(pt, 9.5, GetEventWeight());//eta prime
          else if(mompdg    == 213) fhMCPi0PtOrigin->Fill(pt, 4.5, GetEventWeight());//rho
          else if(mompdg    == 223) fhMCPi0PtOrigin->Fill(pt, 5.5, GetEventWeight());//omega
          else if(mompdg    >= 310   && mompdg    <= 323) fhMCPi0PtOrigin->Fill(pt, 6.5, GetEventWeight());//k0S, k+-,k*
          else if(mompdg    == 130) fhMCPi0PtOrigin->Fill(pt, 6.5, GetEventWeight());//k0L
          else if(momstatus == 11 || momstatus  == 12   ) fhMCPi0PtOrigin->Fill(pt, 3.5, GetEventWeight());//resonances
          else                      fhMCPi0PtOrigin->Fill(pt, 7.5, GetEventWeight());//other?
          
          if(status != 11)
          {
            if     (momstatus  == 21) fhMCNotResonancePi0PtOrigin->Fill(pt, 0.5, GetEventWeight());//parton
            else if(mompdg     < 22 ) fhMCNotResonancePi0PtOrigin->Fill(pt, 1.5, GetEventWeight());//quark
            else if(mompdg     > 2100  && mompdg   < 2210) fhMCNotResonancePi0PtOrigin->Fill(pt, 2.5, GetEventWeight());// resonances
            else if(mompdg    == 221) fhMCNotResonancePi0PtOrigin->Fill(pt, 8.5, GetEventWeight());//eta
            else if(mompdg    == 331) fhMCNotResonancePi0PtOrigin->Fill(pt, 9.5, GetEventWeight());//eta prime
            else if(mompdg    == 213) fhMCNotResonancePi0PtOrigin->Fill(pt, 4.5, GetEventWeight());//rho
            else if(mompdg    == 223) fhMCNotResonancePi0PtOrigin->Fill(pt, 5.5, GetEventWeight());//omega
            else if(mompdg    >= 310   && mompdg    <= 323) fhMCNotResonancePi0PtOrigin->Fill(pt, 6.5, GetEventWeight());//k0S, k+-,k*
            else if(mompdg    == 130) fhMCNotResonancePi0PtOrigin->Fill(pt, 6.5, GetEventWeight());//k0L
            else if(momstatus == 11 || momstatus  == 12   ) fhMCNotResonancePi0PtOrigin->Fill(pt, 3.5, GetEventWeight());//resonances
            else                      fhMCNotResonancePi0PtOrigin->Fill(pt, 7.5, GetEventWeight());//other?
          }
        }
        else if (mcIndex==kmcEta )
        {
          fhMCEtaProdVertex->Fill(pt, prodR, GetEventWeight());
          
          if     (momstatus  == 21) fhMCEtaPtOrigin->Fill(pt, 0.5, GetEventWeight());//parton
          else if(mompdg     < 22 ) fhMCEtaPtOrigin->Fill(pt, 1.5, GetEventWeight());//quark
          else if(mompdg     > 2100  && mompdg   < 2210) fhMCEtaPtOrigin->Fill(pt, 2.5, GetEventWeight());// resonances
          else if(mompdg    == 221) fhMCEtaPtOrigin->Fill(pt, 8.5, GetEventWeight());//eta
          else if(mompdg    == 331) fhMCEtaPtOrigin->Fill(pt, 9.5, GetEventWeight());//eta prime
          else if(mompdg    == 213) fhMCEtaPtOrigin->Fill(pt, 4.5, GetEventWeight());//rho
          else if(mompdg    == 223) fhMCEtaPtOrigin->Fill(pt, 5.5, GetEventWeight());//omega
          else if(mompdg    >= 310   && mompdg    <= 323) fhMCEtaPtOrigin->Fill(pt, 6.5, GetEventWeight());//k0S, k+-,k*
          else if(mompdg    == 130) fhMCEtaPtOrigin->Fill(pt, 6.5, GetEventWeight());//k0L
          else if(momstatus == 11 || momstatus  == 12   ) fhMCEtaPtOrigin->Fill(pt, 3.5, GetEventWeight());//resonances
          else                      fhMCEtaPtOrigin->Fill(pt, 7.5, GetEventWeight());//other?
        }
      }
    } // Histograms with MC
  }// aod loop
  
  AliDebug(1,"End");
}

//__________________________________________________
/// Print some relevant parameters set for the analysis.
//__________________________________________________
void AliAnaPi0EbE::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
 
  AliAnaCaloTrackCorrBaseClass::Print("");
    
  printf("Analysis Type = %d \n",  fAnaType) ;
    
  if(fAnaType == kSSCalo)
  {
    printf("Calorimeter            =     %s\n", GetCalorimeterString().Data()) ;
    printf("Min Distance to Bad Channel   = %2.1f\n",fMinDist);
    printf("Min Distance to Bad Channel 2 = %2.1f\n",fMinDist2);
    printf("Min Distance to Bad Channel 3 = %2.1f\n",fMinDist3);
  }
    
  printf("    \n") ;
}





