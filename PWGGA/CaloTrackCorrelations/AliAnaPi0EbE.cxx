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
// Class for the analysis of high pT pi0 event by event
// Pi0/Eta identified by one of the following:
//  -Invariant mass of 2 cluster in calorimeter
//  -Shower shape analysis in calorimeter
//  -Invariant mass of one cluster in calorimeter and one photon reconstructed in CTS
//
// -- Author: Gustavo Conesa (LNF-INFN) &  Raphaelle Ichou (SUBATECH)
//////////////////////////////////////////////////////////////////////////////


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

ClassImp(AliAnaPi0EbE)

//____________________________
AliAnaPi0EbE::AliAnaPi0EbE() :
AliAnaCaloTrackCorrBaseClass(),     fAnaType(kIMCalo),                  fCalorimeter(""),
fMinDist(0.),fMinDist2(0.),         fMinDist3(0.),
fNLMCutMin(-1),                     fNLMCutMax(10),
fTimeCutMin(-10000),                fTimeCutMax(10000),
fRejectTrackMatch(kTRUE),
fFillPileUpHistograms(0),
fFillWeightHistograms(kFALSE),      fFillTMHisto(0),
fFillSelectClHisto(0),              fFillOnlySimpleSSHisto(1),          fFillEMCALBCHistograms(0),
fInputAODGammaConvName(""),
fCheckSplitDistToBad(0),
// Histograms
fhPt(0),                            fhE(0),
fhEEta(0),                          fhEPhi(0),
fhPtEta(0),                         fhPtPhi(0),                         fhEtaPhi(0),
fhEtaPhiEMCALBC0(0),                fhEtaPhiEMCALBC1(0),                fhEtaPhiEMCALBCN(0),
fhTimeTriggerEMCALBC0UMReMatchOpenTime(0),
fhTimeTriggerEMCALBC0UMReMatchCheckNeigh(0),
fhTimeTriggerEMCALBC0UMReMatchBoth(0),
fhPtCentrality(),                   fhPtEventPlane(0),
fhPtReject(0),                      fhEReject(0),
fhEEtaReject(0),                    fhEPhiReject(0),                    fhEtaPhiReject(0),
fhMass(0),                          fhMassPt(0),                        fhMassSplitPt(0),
fhSelectedMass(0),                  fhSelectedMassPt(0),                fhSelectedMassSplitPt(0),
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
fhPtDecay(0),                       fhEDecay(0),
// Shower shape histos
fhEDispersion(0),                   fhELambda0(0),                      fhELambda1(0),
fhELambda0NoTRD(0),                 fhELambda0FracMaxCellCut(0),
fhEFracMaxCell(0),                  fhEFracMaxCellNoTRD(0),
fhENCells(0),                       fhETime(0),                         fhEPairDiffTime(0),
fhDispEtaE(0),                      fhDispPhiE(0),
fhSumEtaE(0),                       fhSumPhiE(0),                       fhSumEtaPhiE(0),
fhDispEtaPhiDiffE(0),               fhSphericityE(0),

// MC histos
fhMCE(),                            fhMCPt(),
fhMCPhi(),                          fhMCEta(),
fhMCEReject(),                      fhMCPtReject(),
fhMCPtCentrality(),
fhMCPi0PtGenRecoFraction(0),        fhMCEtaPtGenRecoFraction(0),
fhMCPi0DecayPt(0),                  fhMCPi0DecayPtFraction(0),
fhMCEtaDecayPt(0),                  fhMCEtaDecayPtFraction(0),
fhMCOtherDecayPt(0),
fhMassPairMCPi0(0),                 fhMassPairMCEta(0),
fhAnglePairMCPi0(0),                fhAnglePairMCEta(0),
// Weight studies
fhECellClusterRatio(0),             fhECellClusterLogRatio(0),
fhEMaxCellClusterRatio(0),          fhEMaxCellClusterLogRatio(0),
fhTrackMatchedDEta(0),              fhTrackMatchedDPhi(0),              fhTrackMatchedDEtaDPhi(0),
fhTrackMatchedDEtaPos(0),           fhTrackMatchedDPhiPos(0),           fhTrackMatchedDEtaDPhiPos(0),
fhTrackMatchedDEtaNeg(0),           fhTrackMatchedDPhiNeg(0),           fhTrackMatchedDEtaDPhiNeg(0),
fhTrackMatchedMCParticleE(0),
fhTrackMatchedMCParticleDEta(0),    fhTrackMatchedMCParticleDPhi(0),
fhdEdx(0),                          fhEOverP(0),                        fhEOverPNoTRD(0),
// Number of local maxima in cluster
fhNLocMaxE(0),                      fhNLocMaxPt(0),                     fhNLocMaxPtReject(0),
// PileUp
fhTimePtNoCut(0),                   fhTimePtSPD(0),                     fhTimePtSPDMulti(0),
fhTimeNPileUpVertSPD(0),            fhTimeNPileUpVertTrack(0),
fhTimeNPileUpVertContributors(0),
fhTimePileUpMainVertexZDistance(0), fhTimePileUpMainVertexZDiamond(0),
fhPtNPileUpSPDVtx(0),               fhPtNPileUpTrkVtx(0),
fhPtNPileUpSPDVtxTimeCut(0),        fhPtNPileUpTrkVtxTimeCut(0),
fhPtNPileUpSPDVtxTimeCut2(0),       fhPtNPileUpTrkVtxTimeCut2(0)
{
  //default ctor
  
  for(Int_t i = 0; i < 6; i++)
  {
    fhMCE              [i] = 0;
    fhMCPt             [i] = 0;
    fhMCPhi            [i] = 0;
    fhMCEta            [i] = 0;
    fhMCPtCentrality   [i] = 0;
    
    fhMCSplitE         [i] = 0;
    fhMCSplitPt        [i] = 0;
    fhMCSplitPtPhi     [i] = 0;
    fhMCSplitPtEta     [i] = 0;
    
    fhMCNLocMaxPt      [i] = 0;
    fhMCNLocMaxSplitPt [i] = 0;
    fhMCNLocMaxPtReject[i] = 0;
    
    fhEMCLambda0       [i] = 0;
    fhEMCLambda0NoTRD  [i] = 0;
    fhEMCLambda0FracMaxCellCut[i]= 0;
    fhEMCFracMaxCell   [i] = 0;
    fhEMCLambda1       [i] = 0;
    fhEMCDispersion    [i] = 0;
    
    fhMCEDispEta       [i] = 0;
    fhMCEDispPhi       [i] = 0;
    fhMCESumEtaPhi     [i] = 0;
    fhMCEDispEtaPhiDiff[i] = 0;
    fhMCESphericity    [i] = 0;
    fhMCEAsymmetry     [i] = 0;
    
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
    fhELambda0LocMax       [i] = 0;
    fhELambda1LocMax       [i] = 0;
    fhEDispersionLocMax    [i] = 0;
    fhEDispEtaLocMax       [i] = 0;
    fhEDispPhiLocMax       [i] = 0;
    fhESumEtaPhiLocMax     [i] = 0;
    fhEDispEtaPhiDiffLocMax[i] = 0;
    fhESphericityLocMax    [i] = 0;
    fhEAsymmetryLocMax     [i] = 0;
  }
  
  //Weight studies
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
  
  //Initialize parameters
  InitParameters();
  
}

//___________________________________________________________________________________
void AliAnaPi0EbE::FillPileUpHistograms(Float_t pt, Float_t time, AliVCluster * calo)
{
  // Fill some histograms to understand pile-up
  if(!fFillPileUpHistograms) return;
  
  //printf("E %f, time %f\n",energy,time);
  AliVEvent * event = GetReader()->GetInputEvent();
  
  fhTimePtNoCut->Fill(pt,time);
  if(GetReader()->IsPileUpFromSPD())     
  
  if(GetReader()->IsPileUpFromSPD())             { fhPtPileUp[0]->Fill(pt); fhTimePtSPD     ->Fill(pt,time); }
  if(GetReader()->IsPileUpFromEMCal())             fhPtPileUp[1]->Fill(pt);
  if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtPileUp[2]->Fill(pt);
  if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtPileUp[3]->Fill(pt);
  if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtPileUp[4]->Fill(pt);
  if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtPileUp[5]->Fill(pt);
  if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtPileUp[6]->Fill(pt);
  
  if(event->IsPileupFromSPDInMultBins()) fhTimePtSPDMulti->Fill(pt,time);
  
  // cells in cluster
  
  AliVCaloCells* cells = 0;
  if(fCalorimeter == "EMCAL") cells = GetEMCALCells();
  else                        cells = GetPHOSCells();

  Float_t maxCellFraction = 0.;
  Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(cells,calo,maxCellFraction);
  
  Double_t tmax  = cells->GetCellTime(absIdMax);
  GetCaloUtils()->RecalibrateCellTime(tmax, fCalorimeter, absIdMax,GetReader()->GetInputEvent()->GetBunchCrossNumber());
  tmax*=1.e9;
    
  //Loop on cells inside cluster, max cell must be over 100 MeV and time in BC=0
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
        fhPtCellTimePileUp[0]->Fill(pt, timecell);
        fhPtTimeDiffPileUp[0]->Fill(pt, diff);
       }
      
      if(GetReader()->IsPileUpFromEMCal())
      {
        fhPtCellTimePileUp[1]->Fill(pt, timecell);
        fhPtTimeDiffPileUp[1]->Fill(pt, diff);
      }
      
      if(GetReader()->IsPileUpFromSPDOrEMCal())
      {
        fhPtCellTimePileUp[2]->Fill(pt, timecell);
        fhPtTimeDiffPileUp[2]->Fill(pt, diff);
      }
      
      if(GetReader()->IsPileUpFromSPDAndEMCal())
      {
        fhPtCellTimePileUp[3]->Fill(pt, timecell);
        fhPtTimeDiffPileUp[3]->Fill(pt, diff);
      }
      
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())
      {
        fhPtCellTimePileUp[4]->Fill(pt, timecell);
        fhPtTimeDiffPileUp[4]->Fill(pt, diff);
      }
      
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())
      {
        fhPtCellTimePileUp[5]->Fill(pt, timecell);
        fhPtTimeDiffPileUp[5]->Fill(pt, diff);
      }
      
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal())
      {
        fhPtCellTimePileUp[6]->Fill(pt, timecell);
        fhPtTimeDiffPileUp[6]->Fill(pt, diff);
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
    
  }//ESD
  else if (aodEv)
  {
    nVtxSPD = aodEv->GetNumberOfPileupVerticesSPD();
    nVtxTrk = aodEv->GetNumberOfPileupVerticesTracks();
  }//AOD
  
  fhTimeNPileUpVertSPD  ->Fill(time,nVtxSPD);
  fhTimeNPileUpVertTrack->Fill(time,nVtxTrk);
  
	fhPtNPileUpSPDVtx->Fill(pt,nVtxSPD);
	fhPtNPileUpTrkVtx->Fill(pt,nVtxTrk);
	
	if(TMath::Abs(time) < 25)
	{
		fhPtNPileUpSPDVtxTimeCut ->Fill(pt,nVtxSPD);
		fhPtNPileUpTrkVtxTimeCut ->Fill(pt,nVtxTrk);
  }
  
  if(time < 75 && time > -25)
  {
    fhPtNPileUpSPDVtxTimeCut2->Fill(pt,nVtxSPD);
    fhPtNPileUpTrkVtxTimeCut2->Fill(pt,nVtxTrk);
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
    }// AOD
    
    Double_t distZ  = TMath::Abs(z2-z1);
    diamZ  = TMath::Abs(z2-diamZ);
    
    fhTimeNPileUpVertContributors  ->Fill(time,ncont);
    fhTimePileUpMainVertexZDistance->Fill(time,distZ);
    fhTimePileUpMainVertexZDiamond ->Fill(time,diamZ);
    
  }// vertex loop
}


//______________________________________________________________________________________________
void AliAnaPi0EbE::FillRejectedClusterHistograms(TLorentzVector mom, Int_t mctag, Int_t nMaxima)
{
  // Fill histograms that do not pass the identification (SS case only)
  
  Float_t ener  = mom.E();
  Float_t pt    = mom.Pt();
  Float_t phi   = mom.Phi();
  if(phi < 0) phi+=TMath::TwoPi();
  Float_t eta = mom.Eta();
  
  fhPtReject     ->Fill(pt);
  fhEReject      ->Fill(ener);
  
  fhEEtaReject   ->Fill(ener,eta);
  fhEPhiReject   ->Fill(ener,phi);
  fhEtaPhiReject ->Fill(eta,phi);
  
  fhNLocMaxPtReject->Fill(pt,nMaxima);

  if(IsDataMC())
  {
    Int_t mcIndex = GetMCIndex(mctag);
    fhMCEReject  [mcIndex] ->Fill(ener);
    fhMCPtReject [mcIndex] ->Fill(pt);
    fhMCNLocMaxPtReject[mcIndex]->Fill(pt,nMaxima);
  }
}

//___________________________________________________________________________________
void AliAnaPi0EbE::FillSelectedClusterHistograms(AliVCluster* cluster, Int_t nMaxima,
                                                 Int_t tag, Float_t asy)
{
  // Fill shower shape, timing and other histograms for selected clusters from decay
  
  Float_t e    = cluster->E();
  Float_t disp = cluster->GetDispersion()*cluster->GetDispersion();
  Float_t l0   = cluster->GetM02();
  Float_t l1   = cluster->GetM20();
  Int_t   nSM  = GetModuleNumber(cluster);
  
  Int_t ebin = -1;
  if      (e < 2 ) ebin = 0;
  else if (e < 4 ) ebin = 1;
  else if (e < 6 ) ebin = 2;
  else if (e < 10) ebin = 3;
  else if (e < 15) ebin = 4;
  else if (e < 20) ebin = 5;
  else             ebin = 6;
  
  Int_t indexMax = -1;
  if     (nMaxima==1) indexMax = 0 ;
  else if(nMaxima==2) indexMax = 1 ;
  else                indexMax = 2 ;
  
  
  AliVCaloCells * cell = 0x0;
  if(fCalorimeter == "PHOS")
    cell = GetPHOSCells();
  else
    cell = GetEMCALCells();
  
  Float_t maxCellFraction = 0;
  GetCaloUtils()->GetMaxEnergyCell(cell, cluster, maxCellFraction);
  fhEFracMaxCell->Fill(e,maxCellFraction);
  
  FillWeightHistograms(cluster);
  
  fhEDispersion->Fill(e, disp);
  fhELambda0   ->Fill(e, l0  );
  fhELambda1   ->Fill(e, l1  );
  
  Float_t ll0  = 0., ll1  = 0.;
  Float_t dispp= 0., dEta = 0., dPhi    = 0.;
  Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;
  if(fCalorimeter == "EMCAL" && !fFillOnlySimpleSSHisto)
  {
    GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), GetReader()->GetInputEvent()->GetEMCALCells(), cluster,
                                                                                 ll0, ll1, dispp, dEta, dPhi, sEta, sPhi, sEtaPhi);
    
    fhDispEtaE        -> Fill(e,dEta);
    fhDispPhiE        -> Fill(e,dPhi);
    fhSumEtaE         -> Fill(e,sEta);
    fhSumPhiE         -> Fill(e,sPhi);
    fhSumEtaPhiE      -> Fill(e,sEtaPhi);
    fhDispEtaPhiDiffE -> Fill(e,dPhi-dEta);
    if(dEta+dPhi>0)fhSphericityE -> Fill(e,(dPhi-dEta)/(dEta+dPhi));
    
    fhDispEtaDispPhi[ebin]->Fill(dEta,dPhi);
    fhLambda0DispEta[ebin]->Fill(l0  ,dEta);
    fhLambda0DispPhi[ebin]->Fill(l0  ,dPhi);
    
    if (fAnaType==kSSCalo)
    {
      // Asymmetry histograms
      fhAsymmetryLambda0[ebin]->Fill(l0 ,asy);
      fhAsymmetryDispEta[ebin]->Fill(dEta,asy);
      fhAsymmetryDispPhi[ebin]->Fill(dPhi,asy);
    }
  }
  
  fhNLocMaxE ->Fill(e ,nMaxima);
  
  fhELambda0LocMax   [indexMax]->Fill(e,l0);
  fhELambda1LocMax   [indexMax]->Fill(e,l1);
  fhEDispersionLocMax[indexMax]->Fill(e,disp);
  
  if(fCalorimeter=="EMCAL" && !fFillOnlySimpleSSHisto)
  {
    fhEDispEtaLocMax       [indexMax]-> Fill(e,dEta);
    fhEDispPhiLocMax       [indexMax]-> Fill(e,dPhi);
    fhESumEtaPhiLocMax     [indexMax]-> Fill(e,sEtaPhi);
    fhEDispEtaPhiDiffLocMax[indexMax]-> Fill(e,dPhi-dEta);
    if(dEta+dPhi>0)       fhESphericityLocMax[indexMax]->Fill(e,(dPhi-dEta)/(dEta+dPhi));
    if(fAnaType==kSSCalo) fhEAsymmetryLocMax [indexMax]->Fill(e  ,asy);
    
  }
  
  if(fCalorimeter=="EMCAL" && nSM < 6)
  {
    fhELambda0NoTRD->Fill(e, l0  );
    fhEFracMaxCellNoTRD->Fill(e,maxCellFraction);
  }
  
  if(maxCellFraction < 0.5)
    fhELambda0FracMaxCellCut->Fill(e, l0  );
  
  fhETime  ->Fill(e, cluster->GetTOF()*1.e9);
  fhENCells->Fill(e, cluster->GetNCells());
  
  // Fill Track matching control histograms
  if(fFillTMHisto)
  {
    Float_t dZ  = cluster->GetTrackDz();
    Float_t dR  = cluster->GetTrackDx();
    
    if(cluster->IsEMCAL() && GetCaloUtils()->IsRecalculationOfClusterTrackMatchingOn())
    {
      dR = 2000., dZ = 2000.;
      GetCaloUtils()->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dZ,dR);
    }
    //printf("Pi0EbE: dPhi %f, dEta %f\n",dR,dZ);
    
    AliVTrack *track = GetCaloUtils()->GetMatchedTrack(cluster, GetReader()->GetInputEvent());
    
    Bool_t positive = kFALSE;
    if(track) positive = (track->Charge()>0);
    
    if(fhTrackMatchedDEta && TMath::Abs(dR) < 999)
    {
      fhTrackMatchedDEta->Fill(e,dZ);
      fhTrackMatchedDPhi->Fill(e,dR);
      if(e > 0.5) fhTrackMatchedDEtaDPhi->Fill(dZ,dR);
      
      if(track)
      {
        if(positive)
        {
          fhTrackMatchedDEtaPos->Fill(cluster->E(),dZ);
          fhTrackMatchedDPhiPos->Fill(cluster->E(),dR);
          if(cluster->E() > 0.5) fhTrackMatchedDEtaDPhiPos->Fill(dZ,dR);
        }
        else
        {
          fhTrackMatchedDEtaNeg->Fill(cluster->E(),dZ);
          fhTrackMatchedDPhiNeg->Fill(cluster->E(),dR);
          if(cluster->E() > 0.5) fhTrackMatchedDEtaDPhiNeg->Fill(dZ,dR);
        }
    }
    }
    // Check dEdx and E/p of matched clusters
    
    if(TMath::Abs(dZ) < 0.05 && TMath::Abs(dR) < 0.05)
    {      
      if(track)
      {
        Float_t dEdx = track->GetTPCsignal();
        fhdEdx->Fill(e, dEdx);
        
        Float_t eOverp = e/track->P();
        fhEOverP->Fill(e,  eOverp);
        
        if(fCalorimeter=="EMCAL" && nSM < 6) fhEOverPNoTRD->Fill(e,  eOverp);
        
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
        
        fhTrackMatchedMCParticleE   ->Fill(e , mctag);
        fhTrackMatchedMCParticleDEta->Fill(dZ, mctag);
        fhTrackMatchedMCParticleDPhi->Fill(dR, mctag);
        
      }  // MC
    }
  }// Track matching histograms
  
  if(IsDataMC())
  {
    Int_t mcIndex = GetMCIndex(tag);
    
    fhEMCLambda0[mcIndex]    ->Fill(e, l0);
    fhEMCLambda1[mcIndex]    ->Fill(e, l1);
    fhEMCDispersion[mcIndex] ->Fill(e, disp);
    fhEMCFracMaxCell[mcIndex]->Fill(e,maxCellFraction);
    
    if(fCalorimeter=="EMCAL" && nSM < 6)
      fhEMCLambda0NoTRD[mcIndex]->Fill(e, l0  );
    
    if(maxCellFraction < 0.5)
      fhEMCLambda0FracMaxCellCut[mcIndex]->Fill(e, l0  );
    
    if(fCalorimeter == "EMCAL" && !fFillOnlySimpleSSHisto)
    {
      fhMCEDispEta        [mcIndex]-> Fill(e,dEta);
      fhMCEDispPhi        [mcIndex]-> Fill(e,dPhi);
      fhMCESumEtaPhi      [mcIndex]-> Fill(e,sEtaPhi);
      fhMCEDispEtaPhiDiff [mcIndex]-> Fill(e,dPhi-dEta);
      if(dEta+dPhi>0)fhMCESphericity[mcIndex]-> Fill(e,(dPhi-dEta)/(dEta+dPhi));
      
      if (fAnaType==kSSCalo)
      {
        fhMCAsymmetryLambda0[ebin][mcIndex]->Fill(l0 ,asy);
        fhMCAsymmetryDispEta[ebin][mcIndex]->Fill(dEta,asy);
        fhMCAsymmetryDispPhi[ebin][mcIndex]->Fill(dPhi,asy);
      }
      
      fhMCDispEtaDispPhi[ebin][mcIndex]->Fill(dEta,dPhi);
      fhMCLambda0DispEta[ebin][mcIndex]->Fill(l0  ,dEta);
      fhMCLambda0DispPhi[ebin][mcIndex]->Fill(l0  ,dPhi);
      
    }
    
  }//MC
  
}

//________________________________________________________
void AliAnaPi0EbE::FillWeightHistograms(AliVCluster *clus)
{
  // Calculate weights and fill histograms
  
  if(!fFillWeightHistograms || GetMixedEvent()) return;
  
  AliVCaloCells* cells = 0;
  if(fCalorimeter == "EMCAL") cells = GetEMCALCells();
  else                        cells = GetPHOSCells();
  
  // First recalculate energy in case non linearity was applied
  Float_t  energy = 0;
  Float_t  ampMax = 0;
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++)
  {
    
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp,fCalorimeter, id);
    
    energy    += amp;
    
    if(amp> ampMax)
      ampMax = amp;
    
  } // energy loop
  
  if(energy <=0 )
  {
    printf("AliAnaPi0EbE::WeightHistograms()- Wrong calculated energy %f\n",energy);
    return;
  }
  
  fhEMaxCellClusterRatio   ->Fill(energy,ampMax/energy);
  fhEMaxCellClusterLogRatio->Fill(energy,TMath::Log(ampMax/energy));
  
  //Get the ratio and log ratio to all cells in cluster
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++)
  {
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp,fCalorimeter, id);
    
    fhECellClusterRatio   ->Fill(energy,amp/energy);
    fhECellClusterLogRatio->Fill(energy,TMath::Log(amp/energy));
  }
  
  //Recalculate shower shape for different W0
  if(fCalorimeter=="EMCAL"){
    
    Float_t l0org = clus->GetM02();
    Float_t l1org = clus->GetM20();
    Float_t dorg  = clus->GetDispersion();
    
    for(Int_t iw = 0; iw < 14; iw++)
    {
      GetCaloUtils()->GetEMCALRecoUtils()->SetW0(1+iw*0.5);
      GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), cells, clus);
      
      fhLambda0ForW0[iw]->Fill(energy,clus->GetM02());
      //fhLambda1ForW0[iw]->Fill(energy,clus->GetM20());
      
    } // w0 loop
    
    // Set the original values back
    clus->SetM02(l0org);
    clus->SetM20(l1org);
    clus->SetDispersion(dorg);
    
  }// EMCAL
}

//__________________________________________
TObjString * AliAnaPi0EbE::GetAnalysisCuts()
{
	//Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaPi0EbE ---\n") ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fAnaType=%d (Pi0 selection type) \n",fAnaType) ;
  parList+=onePar ;
  
  if(fAnaType == kSSCalo)
  {
    snprintf(onePar,buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fMinDist =%2.2f (Minimal distance to bad channel to accept cluster) \n",fMinDist) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fMinDist2=%2.2f (Cuts on Minimal distance to study acceptance evaluation) \n",fMinDist2) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fMinDist3=%2.2f (One more cut on distance used for acceptance-efficiency study) \n",fMinDist3) ;
    parList+=onePar ;
  }
  
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in PID class.
  if(fAnaType == kSSCalo) parList += GetCaloPID()->GetPIDParametersList() ;
  
  return new TObjString(parList) ;
}

//_____________________________________________
TList *  AliAnaPi0EbE::GetCreateOutputObjects()
{
  // Create histograms to be saved in output file and
  // store them in outputContainer
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
  
  Int_t   ntimebins= GetHistogramRanges()->GetHistoTimeBins();
  Float_t timemax  = GetHistogramRanges()->GetHistoTimeMax();
  Float_t timemin  = GetHistogramRanges()->GetHistoTimeMin();
  
  TString nlm[]   ={"1 Local Maxima","2 Local Maxima", "NLM > 2"};
  TString ptype[] ={"#gamma","#gamma->e^{#pm}","#pi^{0}","#eta","e^{#pm}", "hadron"};
  TString pname[] ={"Photon","Conversion",     "Pi0",    "Eta", "Electron","Hadron"};
  Int_t   bin[]   = {0,2,4,6,10,15,20,100}; // energy bins
  
  fhPt  = new TH1F("hPt","Number of identified  #pi^{0} (#eta) decay",nptbins,ptmin,ptmax);
  fhPt->SetYTitle("N");
  fhPt->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPt) ;
  
  fhE  = new TH1F("hE","Number of identified  #pi^{0} (#eta) decay pairs",nptbins,ptmin,ptmax);
  fhE->SetYTitle("N");
  fhE->SetXTitle("E (GeV)");
  outputContainer->Add(fhE) ;
  
  fhEPhi  = new TH2F
  ("hEPhi","Selected #pi^{0} (#eta) pairs: E vs #phi",nptbins,ptmin,ptmax, nphibins,phimin,phimax);
  fhEPhi->SetYTitle("#phi (rad)");
  fhEPhi->SetXTitle("E (GeV)");
  outputContainer->Add(fhEPhi) ;
  
  fhEEta  = new TH2F
  ("hEEta","Selected #pi^{0} (#eta) pairs: E vs #eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
  fhEEta->SetYTitle("#eta");
  fhEEta->SetXTitle("E (GeV)");
  outputContainer->Add(fhEEta) ;
  
  fhPtPhi  = new TH2F
  ("hPtPhi","Selected #pi^{0} (#eta) pairs: p_{T} vs #phi",nptbins,ptmin,ptmax, nphibins,phimin,phimax);
  fhPtPhi->SetYTitle("#phi (rad)");
  fhPtPhi->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtPhi) ;
  
  fhPtEta  = new TH2F
  ("hPtEta","Selected #pi^{0} (#eta) pairs: p_{T} vs #eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
  fhPtEta->SetYTitle("#eta");
  fhPtEta->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtEta) ;
  
  fhEtaPhi  = new TH2F
  ("hEtaPhi","Selected #pi^{0} (#eta) pairs: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhi->SetYTitle("#phi (rad)");
  fhEtaPhi->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhi) ;
  
  if(fCalorimeter=="EMCAL" && fFillEMCALBCHistograms)
  {
    fhEtaPhiEMCALBC0  = new TH2F
    ("hEtaPhiEMCALBC0","cluster,E > 2 GeV, #eta vs #phi, for clusters with |time| < 25 ns, EMCAL-BC=0",netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiEMCALBC0->SetYTitle("#phi (rad)");
    fhEtaPhiEMCALBC0->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiEMCALBC0) ;
    
    fhEtaPhiEMCALBC1  = new TH2F
    ("hEtaPhiEMCALBC1","cluster,E > 2 GeV, #eta vs #phi, for clusters with 25 < |time| < 75 ns, EMCAL-BC=1",netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiEMCALBC1->SetYTitle("#phi (rad)");
    fhEtaPhiEMCALBC1->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiEMCALBC1) ;
    
    fhEtaPhiEMCALBCN  = new TH2F
    ("hEtaPhiEMCALBCN","cluster,E > 2 GeV, #eta vs #phi, for clusters with |time| > 75 ns, EMCAL-BC>1",netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiEMCALBCN->SetYTitle("#phi (rad)");
    fhEtaPhiEMCALBCN->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiEMCALBCN) ;
    
    for(Int_t i = 0; i < 11; i++)
    {
      fhEtaPhiTriggerEMCALBC[i] = new TH2F
      (Form("hEtaPhiTriggerEMCALBC%d",i-5),
       Form("meson E > 2 GeV, #eta vs #phi, Trigger EMCAL-BC=%d",i-5),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhEtaPhiTriggerEMCALBC[i]->SetYTitle("#phi (rad)");
      fhEtaPhiTriggerEMCALBC[i]->SetXTitle("#eta");
      outputContainer->Add(fhEtaPhiTriggerEMCALBC[i]) ;
      
      fhTimeTriggerEMCALBC[i] = new TH2F
      (Form("hTimeTriggerEMCALBC%d",i-5),
       Form("meson time vs E, Trigger EMCAL-BC=%d",i-5),
       nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
      fhTimeTriggerEMCALBC[i]->SetXTitle("E (GeV)");
      fhTimeTriggerEMCALBC[i]->SetYTitle("time (ns)");
      outputContainer->Add(fhTimeTriggerEMCALBC[i]);
      
      fhTimeTriggerEMCALBCPileUpSPD[i] = new TH2F
      (Form("hTimeTriggerEMCALBC%dPileUpSPD",i-5),
       Form("meson time vs E, Trigger EMCAL-BC=%d",i-5),
       nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
      fhTimeTriggerEMCALBCPileUpSPD[i]->SetXTitle("E (GeV)");
      fhTimeTriggerEMCALBCPileUpSPD[i]->SetYTitle("time (ns)");
      outputContainer->Add(fhTimeTriggerEMCALBCPileUpSPD[i]);
      
      fhEtaPhiTriggerEMCALBCUM[i] = new TH2F
      (Form("hEtaPhiTriggerEMCALBC%d_UnMatch",i-5),
       Form("meson E > 2 GeV, #eta vs #phi, unmatched trigger EMCAL-BC=%d",i-5),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhEtaPhiTriggerEMCALBCUM[i]->SetYTitle("#phi (rad)");
      fhEtaPhiTriggerEMCALBCUM[i]->SetXTitle("#eta");
      outputContainer->Add(fhEtaPhiTriggerEMCALBCUM[i]) ;
      
      fhTimeTriggerEMCALBCUM[i] = new TH2F
      (Form("hTimeTriggerEMCALBC%d_UnMatch",i-5),
       Form("meson time vs E, unmatched trigger EMCAL-BC=%d",i-5),
       nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
      fhTimeTriggerEMCALBCUM[i]->SetXTitle("E (GeV)");
      fhTimeTriggerEMCALBCUM[i]->SetYTitle("time (ns)");
      outputContainer->Add(fhTimeTriggerEMCALBCUM[i]);
      
    }
    
    fhTimeTriggerEMCALBC0UMReMatchOpenTime = new TH2F("hTimeTriggerBC0_UnMatch_ReMatch_OpenTime",
                                                      "cluster time vs E of clusters, no match, rematch open time",
                                                      nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBC0UMReMatchOpenTime->SetXTitle("E (GeV)");
    fhTimeTriggerEMCALBC0UMReMatchOpenTime->SetYTitle("time (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBC0UMReMatchOpenTime);
    
    
    fhTimeTriggerEMCALBC0UMReMatchCheckNeigh = new TH2F("hTimeTriggerBC0_UnMatch_ReMatch_CheckNeighbours",
                                                        "cluster time vs E of clusters, no match, rematch with neigbour parches",
                                                        nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBC0UMReMatchCheckNeigh->SetXTitle("E (GeV)");
    fhTimeTriggerEMCALBC0UMReMatchCheckNeigh->SetYTitle("time (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBC0UMReMatchCheckNeigh);
    
    fhTimeTriggerEMCALBC0UMReMatchBoth = new TH2F("hTimeTriggerBC0_UnMatch_ReMatch_Both",
                                                  "cluster time vs E of clusters, no match, rematch open time and neigbour",
                                                  nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBC0UMReMatchBoth->SetXTitle("E (GeV)");
    fhTimeTriggerEMCALBC0UMReMatchBoth->SetYTitle("time (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBC0UMReMatchBoth);
    
  }
  
  fhPtCentrality  = new TH2F("hPtCentrality","centrality vs p_{T}",nptbins,ptmin,ptmax, 100,0,100);
  fhPtCentrality->SetYTitle("centrality");
  fhPtCentrality->SetXTitle("p_{T}(GeV/c)");
  outputContainer->Add(fhPtCentrality) ;
  
  fhPtEventPlane  = new TH2F("hPtEventPlane","event plane angle vs p_{T}",nptbins,ptmin,ptmax, 100,0,TMath::Pi());
  fhPtEventPlane->SetYTitle("Event plane angle (rad)");
  fhPtEventPlane->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtEventPlane) ;
  
  if(fAnaType == kSSCalo)
  {
    fhPtReject  = new TH1F("hPtReject","Number of rejected as #pi^{0} (#eta) decay",nptbins,ptmin,ptmax);
    fhPtReject->SetYTitle("N");
    fhPtReject->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtReject) ;
    
    fhEReject  = new TH1F("hEReject","Number of rejected as  #pi^{0} (#eta) decay pairs",nptbins,ptmin,ptmax);
    fhEReject->SetYTitle("N");
    fhEReject->SetXTitle("E (GeV)");
    outputContainer->Add(fhEReject) ;
    
    fhEPhiReject  = new TH2F
    ("hEPhiReject","Rejected #pi^{0} (#eta) cluster: E vs #phi",nptbins,ptmin,ptmax, nphibins,phimin,phimax);
    fhEPhiReject->SetYTitle("#phi (rad)");
    fhEPhiReject->SetXTitle("E (GeV)");
    outputContainer->Add(fhEPhiReject) ;
    
    fhEEtaReject  = new TH2F
    ("hEEtaReject","Rejected #pi^{0} (#eta) cluster: E vs #eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhEEtaReject->SetYTitle("#eta");
    fhEEtaReject->SetXTitle("E (GeV)");
    outputContainer->Add(fhEEtaReject) ;
    
    fhEtaPhiReject  = new TH2F
    ("hEtaPhiReject","Rejected #pi^{0} (#eta) cluster: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
    fhEtaPhiReject->SetYTitle("#phi (rad)");
    fhEtaPhiReject->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiReject) ;
  }
  
  fhMass  = new TH2F
  ("hMass","all pairs mass: E vs mass",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
  fhMass->SetYTitle("mass (GeV/c^{2})");
  fhMass->SetXTitle("E (GeV)");
  outputContainer->Add(fhMass) ;
  
  fhSelectedMass  = new TH2F
  ("hSelectedMass","Selected #pi^{0} (#eta) pairs mass: E vs mass",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
  fhSelectedMass->SetYTitle("mass (GeV/c^{2})");
  fhSelectedMass->SetXTitle("E (GeV)");
  outputContainer->Add(fhSelectedMass) ;
  
  fhMassPt  = new TH2F
  ("hMassPt","all pairs mass: p_{T} vs mass",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
  fhMassPt->SetYTitle("mass (GeV/c^{2})");
  fhMassPt->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhMassPt) ;
  
  fhSelectedMassPt  = new TH2F
  ("hSelectedMassPt","Selected #pi^{0} (#eta) pairs mass: p_{T} vs mass",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
  fhSelectedMassPt->SetYTitle("mass (GeV/c^{2})");
  fhSelectedMassPt->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhSelectedMassPt) ;

  if(IsDataMC() && fAnaType == kSSCalo)
  {
    fhMassNoOverlap  = new TH2F
    ("hMassNoOverlap","all pairs mass: E vs mass, no overlap",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
    fhMassNoOverlap->SetYTitle("mass (GeV/c^{2})");
    fhMassNoOverlap->SetXTitle("E (GeV)");
    outputContainer->Add(fhMassNoOverlap) ;
    
    fhSelectedMassNoOverlap  = new TH2F
    ("hSelectedMassNoOverlap","Selected #pi^{0} (#eta) pairs mass: E vs mass, no overlap",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
    fhSelectedMassNoOverlap->SetYTitle("mass (GeV/c^{2})");
    fhSelectedMassNoOverlap->SetXTitle("E (GeV)");
    outputContainer->Add(fhSelectedMassNoOverlap) ;
    
    fhMassPtNoOverlap  = new TH2F
    ("hMassPtNoOverlap","all pairs mass: p_{T} vs mass, no overlap",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
    fhMassPtNoOverlap->SetYTitle("mass (GeV/c^{2})");
    fhMassPtNoOverlap->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMassPtNoOverlap) ;
    
    fhSelectedMassPtNoOverlap  = new TH2F
    ("hSelectedMassPtNoOverlap","Selected #pi^{0} (#eta) pairs mass: p_{T} vs mass, no overlap",nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
    fhSelectedMassPtNoOverlap->SetYTitle("mass (GeV/c^{2})");
    fhSelectedMassPtNoOverlap->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhSelectedMassPtNoOverlap) ;
  }
  
  if(fAnaType != kSSCalo)
  {
    fhPtDecay  = new TH1F("hPtDecay","Number of identified  #pi^{0} (#eta) decay photons",nptbins,ptmin,ptmax);
    fhPtDecay->SetYTitle("N");
    fhPtDecay->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtDecay) ;
    
    fhEDecay  = new TH1F("hEDecay","Number of identified  #pi^{0} (#eta) decay photons",nptbins,ptmin,ptmax);
    fhEDecay->SetYTitle("N");
    fhEDecay->SetXTitle("E (GeV)");
    outputContainer->Add(fhEDecay) ;
  }
  
  ////////
  
  if( fFillSelectClHisto )
  {
    
    fhEDispersion  = new TH2F
    ("hEDispersion","Selected #pi^{0} (#eta) pairs: E vs dispersion",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
    fhEDispersion->SetYTitle("D^{2}");
    fhEDispersion->SetXTitle("E (GeV)");
    outputContainer->Add(fhEDispersion) ;
    
    fhELambda0  = new TH2F
    ("hELambda0","Selected #pi^{0} (#eta) pairs: E vs #lambda_{0}",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
    fhELambda0->SetYTitle("#lambda_{0}^{2}");
    fhELambda0->SetXTitle("E (GeV)");
    outputContainer->Add(fhELambda0) ;
    
    fhELambda1  = new TH2F
    ("hELambda1","Selected #pi^{0} (#eta) pairs: E vs #lambda_{1}",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
    fhELambda1->SetYTitle("#lambda_{1}^{2}");
    fhELambda1->SetXTitle("E (GeV)");
    outputContainer->Add(fhELambda1) ;
    
    fhELambda0FracMaxCellCut  = new TH2F
    ("hELambda0FracMaxCellCut","Selected #pi^{0} (#eta) pairs: E vs #lambda_{0}, Max cell fraction of energy < 0.5",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
    fhELambda0FracMaxCellCut->SetYTitle("#lambda_{0}^{2}");
    fhELambda0FracMaxCellCut->SetXTitle("E (GeV)");
    outputContainer->Add(fhELambda0FracMaxCellCut) ;
    
    fhEFracMaxCell  = new TH2F
    ("hEFracMaxCell","Selected #pi^{0} (#eta) pairs: E vs #lambda_{0}, Max cell fraction of energy",nptbins,ptmin,ptmax,100,0,1);
    fhEFracMaxCell->SetYTitle("Fraction");
    fhEFracMaxCell->SetXTitle("E (GeV)");
    outputContainer->Add(fhEFracMaxCell) ;
    
    if(fCalorimeter=="EMCAL")
    {
      fhELambda0NoTRD  = new TH2F
      ("hELambda0NoTRD","Selected #pi^{0} (#eta) pairs: E vs #lambda_{0}, not behind TRD",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhELambda0NoTRD->SetYTitle("#lambda_{0}^{2}");
      fhELambda0NoTRD->SetXTitle("E (GeV)");
      outputContainer->Add(fhELambda0NoTRD) ;
      
      fhEFracMaxCellNoTRD  = new TH2F
      ("hEFracMaxCellNoTRD","Selected #pi^{0} (#eta) pairs: E vs #lambda_{0}, Max cell fraction of energy, not behind TRD",nptbins,ptmin,ptmax,100,0,1);
      fhEFracMaxCellNoTRD->SetYTitle("Fraction");
      fhEFracMaxCellNoTRD->SetXTitle("E (GeV)");
      outputContainer->Add(fhEFracMaxCellNoTRD) ;
      
      if(!fFillOnlySimpleSSHisto)
      {
        fhDispEtaE  = new TH2F ("hDispEtaE","#sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i} - <#eta>)^{2}/ #Sigma w_{i} vs E",  nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
        fhDispEtaE->SetXTitle("E (GeV)");
        fhDispEtaE->SetYTitle("#sigma^{2}_{#eta #eta}");
        outputContainer->Add(fhDispEtaE);
        
        fhDispPhiE  = new TH2F ("hDispPhiE","#sigma^{2}_{#phi #phi} = #Sigma w_{i}(#phi_{i} - <#phi>)^{2} / #Sigma w_{i} vs E",  nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
        fhDispPhiE->SetXTitle("E (GeV)");
        fhDispPhiE->SetYTitle("#sigma^{2}_{#phi #phi}");
        outputContainer->Add(fhDispPhiE);
        
        fhSumEtaE  = new TH2F ("hSumEtaE","#sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i})^{2} / #Sigma w_{i} - <#eta>^{2} vs E",  nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
        fhSumEtaE->SetXTitle("E (GeV)");
        fhSumEtaE->SetYTitle("#delta^{2}_{#eta #eta}");
        outputContainer->Add(fhSumEtaE);
        
        fhSumPhiE  = new TH2F ("hSumPhiE","#sigma^{2}_{#phi #phi} = #Sigma w_{i}(#phi_{i})^{2}/ #Sigma w_{i} - <#phi>^{2} vs E",
                               nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
        fhSumPhiE->SetXTitle("E (GeV)");
        fhSumPhiE->SetYTitle("#delta^{2}_{#phi #phi}");
        outputContainer->Add(fhSumPhiE);
        
        fhSumEtaPhiE  = new TH2F ("hSumEtaPhiE","#delta^{2}_{#eta #phi} = #Sigma w_{i}(#phi_{i} #eta_{i} ) / #Sigma w_{i} - <#phi><#eta> vs E",
                                  nptbins,ptmin,ptmax, 2*ssbins,-ssmax,ssmax);
        fhSumEtaPhiE->SetXTitle("E (GeV)");
        fhSumEtaPhiE->SetYTitle("#delta^{2}_{#eta #phi}");
        outputContainer->Add(fhSumEtaPhiE);
        
        fhDispEtaPhiDiffE  = new TH2F ("hDispEtaPhiDiffE","#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta} vs E",
                                       nptbins,ptmin,ptmax,200, -10,10);
        fhDispEtaPhiDiffE->SetXTitle("E (GeV)");
        fhDispEtaPhiDiffE->SetYTitle("#sigma^{2}_{#phi #phi}-#sigma^{2}_{#eta #eta}");
        outputContainer->Add(fhDispEtaPhiDiffE);
        
        fhSphericityE  = new TH2F ("hSphericityE","(#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi}) vs E",
                                   nptbins,ptmin,ptmax, 200, -1,1);
        fhSphericityE->SetXTitle("E (GeV)");
        fhSphericityE->SetYTitle("s = (#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi})");
        outputContainer->Add(fhSphericityE);
        
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
    }
    
    fhNLocMaxE = new TH2F("hNLocMaxE","Number of local maxima in cluster",
                          nptbins,ptmin,ptmax,10,0,10);
    fhNLocMaxE ->SetYTitle("N maxima");
    fhNLocMaxE ->SetXTitle("E (GeV)");
    outputContainer->Add(fhNLocMaxE) ;
    
    if(fAnaType == kSSCalo)
    {
      fhNLocMaxPt = new TH2F("hNLocMaxPt","Number of local maxima in cluster, selected clusters",
                             nptbins,ptmin,ptmax,10,0,10);
      fhNLocMaxPt ->SetYTitle("N maxima");
      fhNLocMaxPt ->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhNLocMaxPt) ;

      fhNLocMaxPtReject = new TH2F("hNLocMaxPtReject","Number of local maxima in cluster, rejected clusters",
                             nptbins,ptmin,ptmax,10,0,10);
      fhNLocMaxPtReject ->SetYTitle("N maxima");
      fhNLocMaxPtReject ->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhNLocMaxPtReject) ;
    }
    
    for (Int_t i = 0; i < 3; i++)
    {
      fhELambda0LocMax[i]  = new TH2F(Form("hELambda0LocMax%d",i+1),
                                      Form("Selected #pi^{0} (#eta) pairs: E vs #lambda_{0}, %s",nlm[i].Data()),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhELambda0LocMax[i]->SetYTitle("#lambda_{0}^{2}");
      fhELambda0LocMax[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhELambda0LocMax[i]) ;
      
      fhELambda1LocMax[i]  = new TH2F(Form("hELambda1LocMax%d",i+1),
                                      Form("Selected #pi^{0} (#eta) pairs: E vs #lambda_{1}, %s",nlm[i].Data()),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhELambda1LocMax[i]->SetYTitle("#lambda_{1}^{2}");
      fhELambda1LocMax[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhELambda1LocMax[i]) ;
      
      fhEDispersionLocMax[i]  = new TH2F(Form("hEDispersionLocMax%d",i+1),
                                         Form("Selected #pi^{0} (#eta) pairs: E vs dispersion^{2}, %s",nlm[i].Data()),
                                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhEDispersionLocMax[i]->SetYTitle("dispersion^{2}");
      fhEDispersionLocMax[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhEDispersionLocMax[i]) ;
      
      if(fCalorimeter == "EMCAL" && !fFillOnlySimpleSSHisto)
      {
        fhEDispEtaLocMax[i]  = new TH2F(Form("hEDispEtaLocMax%d",i+1),
                                        Form("Selected #pi^{0} (#eta) pairs: E vs #sigma_{#eta #eta}, %s",nlm[i].Data()),
                                        nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhEDispEtaLocMax[i]->SetYTitle("#sigma_{#eta #eta}");
        fhEDispEtaLocMax[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhEDispEtaLocMax[i]) ;
        
        fhEDispPhiLocMax[i]  = new TH2F(Form("hEDispPhiLocMax%d",i+1),
                                        Form("Selected #pi^{0} (#eta) pairs: E vs #sigma_{#phi #phi}, %s",nlm[i].Data()),
                                        nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhEDispPhiLocMax[i]->SetYTitle("#sigma_{#phi #phi}");
        fhEDispPhiLocMax[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhEDispPhiLocMax[i]) ;
        
        fhESumEtaPhiLocMax[i]  = new TH2F(Form("hESumEtaPhiLocMax%d",i+1),
                                          Form("Selected #pi^{0} (#eta) pairs: E vs #sigma_{#eta #phi}, %s",nlm[i].Data()),
                                          nptbins,ptmin,ptmax,2*ssbins,-ssmax,ssmax);
        fhESumEtaPhiLocMax[i]->SetYTitle("#sigma_{#eta #phi}");
        fhESumEtaPhiLocMax[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhESumEtaPhiLocMax[i]) ;
        
        fhEDispEtaPhiDiffLocMax[i]  = new TH2F(Form("hEDispEtaPhiDiffLocMax%d",i+1),
                                               Form("Selected #pi^{0} (#eta) pairs: E vs #sigma_{#phi #phi} - #sigma_{#eta #eta}, %s",nlm[i].Data()),
                                               nptbins,ptmin,ptmax,200, -10,10);
        fhEDispEtaPhiDiffLocMax[i]->SetYTitle("#sigma_{#phi #phi} - #sigma_{#eta #eta}");
        fhEDispEtaPhiDiffLocMax[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhEDispEtaPhiDiffLocMax[i]) ;
        
        fhESphericityLocMax[i]  = new TH2F(Form("hESphericityLocMax%d",i+1),
                                           Form("Selected #pi^{0} (#eta) pairs: E vs #sigma_{#phi #phi} - #sigma_{#eta #eta} / (#sigma_{#phi #phi} + #sigma_{#eta #eta}), %s",nlm[i].Data()),
                                           nptbins,ptmin,ptmax,200, -1,1);
        fhESphericityLocMax[i]->SetYTitle("#sigma_{#phi #phi} - #sigma_{#eta #eta} / (#sigma_{#phi #phi} + #sigma_{#eta #eta})");
        fhESphericityLocMax[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhESphericityLocMax[i]) ;
      }
      
    }
    
    fhENCells  = new TH2F ("hENCells","N cells in cluster vs E ", nptbins,ptmin,ptmax, nbins,nmin,nmax);
    fhENCells->SetXTitle("E (GeV)");
    fhENCells->SetYTitle("# of cells in cluster");
    outputContainer->Add(fhENCells);
    
    fhETime = new TH2F("hETime","cluster time vs pair E",nptbins,ptmin,ptmax, tbins,tmin,tmax);
    fhETime->SetXTitle("E (GeV)");
    fhETime->SetYTitle("t (ns)");
    outputContainer->Add(fhETime);
    
  }
  
  
  fhEPairDiffTime = new TH2F("hEPairDiffTime","cluster pair time difference vs E",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
  fhEPairDiffTime->SetXTitle("E_{pair} (GeV)");
  fhEPairDiffTime->SetYTitle("#Delta t (ns)");
  outputContainer->Add(fhEPairDiffTime);
  
  if(fAnaType == kIMCalo)
  {
    TString combiName [] = {"1LocMax","2LocMax","NLocMax","1LocMax2LocMax","1LocMaxNLocMax","2LocMaxNLocMax","1LocMaxSSBad","NLocMaxSSGood"};
    TString combiTitle[] = {"1 Local Maxima in both clusters","2 Local Maxima in both clusters","more than 2 Local Maxima in both clusters",
      "1 Local Maxima paired with 2 Local Maxima","1 Local Maxima paired with more than 2 Local Maxima",
      "2 Local Maxima paired with more than 2 Local Maxima",
      "1 Local Maxima paired with #lambda_{0}^{2}>0.3","N Local Maxima paired with 0.1<#lambda_{0}^{2}<0.3"};
    
    for (Int_t i = 0; i < 8 ; i++)
    {
      
      if (fAnaType == kIMCaloTracks && i > 2 ) continue ;
      
      fhMassPairLocMax[i]  = new TH2F
      (Form("MassPairLocMax%s",combiName[i].Data()),
       Form("Mass for decay #gamma pair vs E_{pair}, origin #pi^{0}, %s", combiTitle[i].Data()),
       nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
      fhMassPairLocMax[i]->SetYTitle("Mass (MeV/c^{2})");
      fhMassPairLocMax[i]->SetXTitle("E_{pair} (GeV)");
      outputContainer->Add(fhMassPairLocMax[i]) ;
    }
  }
  
  if(fFillTMHisto)
  {
    fhTrackMatchedDEta  = new TH2F
    ("hTrackMatchedDEta",
     "d#eta of cluster-track vs cluster energy",
     nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
    fhTrackMatchedDEta->SetYTitle("d#eta");
    fhTrackMatchedDEta->SetXTitle("E_{cluster} (GeV)");
    
    fhTrackMatchedDPhi  = new TH2F
    ("hTrackMatchedDPhi",
     "d#phi of cluster-track vs cluster energy",
     nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDPhi->SetYTitle("d#phi (rad)");
    fhTrackMatchedDPhi->SetXTitle("E_{cluster} (GeV)");
    
    fhTrackMatchedDEtaDPhi  = new TH2F
    ("hTrackMatchedDEtaDPhi",
     "d#eta vs d#phi of cluster-track vs cluster energy",
     nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDEtaDPhi->SetYTitle("d#phi (rad)");
    fhTrackMatchedDEtaDPhi->SetXTitle("d#eta");
    
    outputContainer->Add(fhTrackMatchedDEta) ;
    outputContainer->Add(fhTrackMatchedDPhi) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhi) ;

    fhTrackMatchedDEtaPos  = new TH2F
    ("hTrackMatchedDEtaPos",
     "d#eta of cluster-track vs cluster energy",
     nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
    fhTrackMatchedDEtaPos->SetYTitle("d#eta");
    fhTrackMatchedDEtaPos->SetXTitle("E_{cluster} (GeV)");
    
    fhTrackMatchedDPhiPos  = new TH2F
    ("hTrackMatchedDPhiPos",
     "d#phi of cluster-track vs cluster energy",
     nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDPhiPos->SetYTitle("d#phi (rad)");
    fhTrackMatchedDPhiPos->SetXTitle("E_{cluster} (GeV)");
    
    fhTrackMatchedDEtaDPhiPos  = new TH2F
    ("hTrackMatchedDEtaDPhiPos",
     "d#eta vs d#phi of cluster-track vs cluster energy",
     nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDEtaDPhiPos->SetYTitle("d#phi (rad)");
    fhTrackMatchedDEtaDPhiPos->SetXTitle("d#eta");
    
    outputContainer->Add(fhTrackMatchedDEtaPos) ;
    outputContainer->Add(fhTrackMatchedDPhiPos) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhiPos) ;

    fhTrackMatchedDEtaNeg  = new TH2F
    ("hTrackMatchedDEtaNeg",
     "d#eta of cluster-track vs cluster energy",
     nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
    fhTrackMatchedDEtaNeg->SetYTitle("d#eta");
    fhTrackMatchedDEtaNeg->SetXTitle("E_{cluster} (GeV)");
    
    fhTrackMatchedDPhiNeg  = new TH2F
    ("hTrackMatchedDPhiNeg",
     "d#phi of cluster-track vs cluster energy",
     nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDPhiNeg->SetYTitle("d#phi (rad)");
    fhTrackMatchedDPhiNeg->SetXTitle("E_{cluster} (GeV)");
    
    fhTrackMatchedDEtaDPhiNeg  = new TH2F
    ("hTrackMatchedDEtaDPhiNeg",
     "d#eta vs d#phi of cluster-track vs cluster energy",
     nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDEtaDPhiNeg->SetYTitle("d#phi (rad)");
    fhTrackMatchedDEtaDPhiNeg->SetXTitle("d#eta");
    
    outputContainer->Add(fhTrackMatchedDEtaNeg) ;
    outputContainer->Add(fhTrackMatchedDPhiNeg) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhiNeg) ;
    
    fhdEdx  = new TH2F ("hdEdx","matched track <dE/dx> vs cluster E ", nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
    fhdEdx->SetXTitle("E (GeV)");
    fhdEdx->SetYTitle("<dE/dx>");
    outputContainer->Add(fhdEdx);
    
    fhEOverP  = new TH2F ("hEOverP","matched track E/p vs cluster E ", nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
    fhEOverP->SetXTitle("E (GeV)");
    fhEOverP->SetYTitle("E/p");
    outputContainer->Add(fhEOverP);
    
    if(fCalorimeter=="EMCAL")
    {
      fhEOverPNoTRD  = new TH2F ("hEOverPNoTRD","matched track E/p vs cluster E, SM not behind TRD ", nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
      fhEOverPNoTRD->SetXTitle("E (GeV)");
      fhEOverPNoTRD->SetYTitle("E/p");
      outputContainer->Add(fhEOverPNoTRD);
    }
    
    if(IsDataMC() && fFillTMHisto)
    {
      fhTrackMatchedMCParticleE  = new TH2F
      ("hTrackMatchedMCParticleE",
       "Origin of particle vs energy",
       nptbins,ptmin,ptmax,8,0,8);
      fhTrackMatchedMCParticleE->SetXTitle("E (GeV)");
      //fhTrackMatchedMCParticleE->SetYTitle("Particle type");
      
      fhTrackMatchedMCParticleE->GetYaxis()->SetBinLabel(1 ,"Photon");
      fhTrackMatchedMCParticleE->GetYaxis()->SetBinLabel(2 ,"Electron");
      fhTrackMatchedMCParticleE->GetYaxis()->SetBinLabel(3 ,"Meson Merged");
      fhTrackMatchedMCParticleE->GetYaxis()->SetBinLabel(4 ,"Rest");
      fhTrackMatchedMCParticleE->GetYaxis()->SetBinLabel(5 ,"Conv. Photon");
      fhTrackMatchedMCParticleE->GetYaxis()->SetBinLabel(6 ,"Conv. Electron");
      fhTrackMatchedMCParticleE->GetYaxis()->SetBinLabel(7 ,"Conv. Merged");
      fhTrackMatchedMCParticleE->GetYaxis()->SetBinLabel(8 ,"Conv. Rest");
      
      outputContainer->Add(fhTrackMatchedMCParticleE);
      
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
    fhECellClusterRatio->SetXTitle("E_{cluster} (GeV) ");
    fhECellClusterRatio->SetYTitle("E_{cell i}/E_{cluster}");
    outputContainer->Add(fhECellClusterRatio);
    
    fhECellClusterLogRatio  = new TH2F ("hECellClusterLogRatio"," Log(cell energy / cluster energy) vs cluster energy, for selected decay photons from neutral meson",
                                        nptbins,ptmin,ptmax, 100,-10,0);
    fhECellClusterLogRatio->SetXTitle("E_{cluster} (GeV) ");
    fhECellClusterLogRatio->SetYTitle("Log (E_{max cell}/E_{cluster})");
    outputContainer->Add(fhECellClusterLogRatio);
    
    fhEMaxCellClusterRatio  = new TH2F ("hEMaxCellClusterRatio"," max cell energy / cluster energy vs cluster energy, for selected decay photons from neutral meson",
                                        nptbins,ptmin,ptmax, 100,0,1.);
    fhEMaxCellClusterRatio->SetXTitle("E_{cluster} (GeV) ");
    fhEMaxCellClusterRatio->SetYTitle("E_{max cell}/E_{cluster}");
    outputContainer->Add(fhEMaxCellClusterRatio);
    
    fhEMaxCellClusterLogRatio  = new TH2F ("hEMaxCellClusterLogRatio"," Log(max cell energy / cluster energy) vs cluster energy, for selected decay photons from neutral meson",
                                           nptbins,ptmin,ptmax, 100,-10,0);
    fhEMaxCellClusterLogRatio->SetXTitle("E_{cluster} (GeV) ");
    fhEMaxCellClusterLogRatio->SetYTitle("Log (E_{max cell}/E_{cluster})");
    outputContainer->Add(fhEMaxCellClusterLogRatio);
    
    for(Int_t iw = 0; iw < 14; iw++)
    {
      fhLambda0ForW0[iw]  = new TH2F (Form("hLambda0ForW0%d",iw),Form("shower shape, #lambda^{2}_{0} vs E, w0 = %1.1f, for selected decay photons from neutral meson",1+0.5*iw),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhLambda0ForW0[iw]->SetXTitle("E_{cluster}");
      fhLambda0ForW0[iw]->SetYTitle("#lambda^{2}_{0}");
      outputContainer->Add(fhLambda0ForW0[iw]);
      
      //      fhLambda1ForW0[iw]  = new TH2F (Form("hLambda1ForW0%d",iw),Form("shower shape, #lambda^{2}_{1} vs E, w0 = %1.1f, for selected decay photons from neutral meson",0.5+0.5*iw),
      //                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      //      fhLambda1ForW0[iw]->SetXTitle("E_{cluster}");
      //      fhLambda1ForW0[iw]->SetYTitle("#lambda^{2}_{1}");
      //      outputContainer->Add(fhLambda1ForW0[iw]);
      
    }
  }
  
  if(IsDataMC())
  {
    if(GetReader()->GetDataType() != AliCaloTrackReader::kMC && fAnaType==kSSCalo)
    {
      fhMCPi0PtGenRecoFraction = new TH2F("hMCPi0PtGenRecoFraction","Number of clusters from #pi^{0} (2 #gamma) identified as #pi^{0} (#eta), pT versus E primary #pi^{0} / E reco",
                                          nptbins,ptmin,ptmax,200,0,2);
      fhMCPi0PtGenRecoFraction->SetXTitle("p^{rec}_{T} (GeV/c)");
      fhMCPi0PtGenRecoFraction->SetYTitle("E^{ #pi^{0} mother} / E^{rec}");
      outputContainer->Add(fhMCPi0PtGenRecoFraction) ;
      
      fhMCEtaPtGenRecoFraction = new TH2F("hMCEtaPtGenRecoFraction","Number of clusters from #eta (2 #gamma) identified as #pi^{0} (#eta),pT versus E primary #eta / E reco",
                                          nptbins,ptmin,ptmax,200,0,2);
      fhMCEtaPtGenRecoFraction->SetXTitle("p^{rec}_{T} (GeV/c)");
      fhMCEtaPtGenRecoFraction->SetYTitle("E^{ #eta mother} / E^{rec}");
      outputContainer->Add(fhMCEtaPtGenRecoFraction) ;
      
      fhMCPi0DecayPt = new TH1F("hMCPi0DecayPt","Number of #gamma from #pi^{0} decay  identified as #pi^{0} (#eta)",nptbins,ptmin,ptmax);
      fhMCPi0DecayPt->SetYTitle("N");
      fhMCPi0DecayPt->SetXTitle("p^{rec}_{T} (GeV/c)");
      outputContainer->Add(fhMCPi0DecayPt) ;
      
      fhMCPi0DecayPtFraction = new TH2F("hMCPi0DecayPtFraction","Number of #gamma from #pi^{0} decay  identified as #pi^{0} (#eta), pT versus E primary  #gamma / E primary #pi^{0}",
                                        nptbins,ptmin,ptmax,100,0,1);
      fhMCPi0DecayPtFraction->SetXTitle("p^{rec}_{T} (GeV/c)");
      fhMCPi0DecayPtFraction->SetYTitle("E^{gen} / E^{gen-mother}");
      outputContainer->Add(fhMCPi0DecayPtFraction) ;
      
      fhMCEtaDecayPt = new TH1F("hMCEtaDecayPt","Number of #gamma from #eta decay  identified as #pi^{0} (#eta)",nptbins,ptmin,ptmax);
      fhMCEtaDecayPt->SetYTitle("N");
      fhMCEtaDecayPt->SetXTitle("p^{rec}_{T} (GeV/c)");
      outputContainer->Add(fhMCEtaDecayPt) ;
      
      fhMCEtaDecayPtFraction = new TH2F("hMCEtaDecayPtFraction","Number of #gamma from #eta decay  identified as #pi^{0} (#eta), pT versus E primary  #gamma / E primary #eta",
                                        nptbins,ptmin,ptmax,100,0,1);
      fhMCEtaDecayPtFraction->SetXTitle("p^{rec}_{T} (GeV/c)");
      fhMCEtaDecayPtFraction->SetYTitle("E^{gen} / E^{gen-mother}");
      outputContainer->Add(fhMCEtaDecayPtFraction) ;
      
      fhMCOtherDecayPt = new TH1F("hMCOtherDecayPt","Number of #gamma decay (not #eta or #pi^{0})  identified as #pi^{0} (#eta)",nptbins,ptmin,ptmax);
      fhMCOtherDecayPt->SetYTitle("N");
      fhMCOtherDecayPt->SetXTitle("p^{rec}_{T} (GeV/c)");
      outputContainer->Add(fhMCOtherDecayPt) ;
      
    }
    
    if((GetReader()->GetDataType() == AliCaloTrackReader::kMC && fAnaType!=kSSCalo) ||
       GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    {
      
      fhAnglePairMCPi0  = new TH2F
      ("AnglePairMCPi0",
       "Angle between decay #gamma pair vs E_{pair}, origin #pi^{0}",nptbins,ptmin,ptmax,250,0,0.5);
      fhAnglePairMCPi0->SetYTitle("#alpha (rad)");
      fhAnglePairMCPi0->SetXTitle("E_{pair} (GeV)");
      outputContainer->Add(fhAnglePairMCPi0) ;
      
      if (fAnaType!= kSSCalo)
      {
        fhAnglePairMCEta  = new TH2F
        ("AnglePairMCEta",
         "Angle between decay #gamma pair vs E_{pair}, origin #eta",nptbins,ptmin,ptmax,250,0,0.5);
        fhAnglePairMCEta->SetYTitle("#alpha (rad)");
        fhAnglePairMCEta->SetXTitle("E_{pair} (GeV)");
        outputContainer->Add(fhAnglePairMCEta) ;
        
        fhMassPairMCPi0  = new TH2F
        ("MassPairMCPi0",
         "Mass for decay #gamma pair vs E_{pair}, origin #pi^{0}",nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhMassPairMCPi0->SetYTitle("Mass (MeV/c^{2})");
        fhMassPairMCPi0->SetXTitle("E_{pair} (GeV)");
        outputContainer->Add(fhMassPairMCPi0) ;
        
        fhMassPairMCEta  = new TH2F
        ("MassPairMCEta",
         "Mass for decay #gamma pair vs E_{pair}, origin #eta",nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhMassPairMCEta->SetYTitle("Mass (MeV/c^{2})");
        fhMassPairMCEta->SetXTitle("E_{pair} (GeV)");
        outputContainer->Add(fhMassPairMCEta) ;
      }
      
      for(Int_t i = 0; i < 6; i++)
      {
        
        fhMCE[i]  = new TH1F
        (Form("hE_MC%s",pname[i].Data()),
         Form("Identified as #pi^{0} (#eta), cluster from %s",
              ptype[i].Data()),
         nptbins,ptmin,ptmax);
        fhMCE[i]->SetYTitle("N");
        fhMCE[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCE[i]) ;
        
        fhMCPt[i]  = new TH1F
        (Form("hPt_MC%s",pname[i].Data()),
         Form("Identified as #pi^{0} (#eta), cluster from %s",
              ptype[i].Data()),
         nptbins,ptmin,ptmax);
        fhMCPt[i]->SetYTitle("N");
        fhMCPt[i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCPt[i]) ;
        
        fhMCPtCentrality[i]  = new TH2F
        (Form("hPtCentrality_MC%s",pname[i].Data()),
         Form("Identified as #pi^{0} (#eta), cluster from %s",
              ptype[i].Data()),
         nptbins,ptmin,ptmax, 100,0,100);
        fhMCPtCentrality[i]->SetYTitle("centrality");
        fhMCPtCentrality[i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCPtCentrality[i]) ;
        
        if(fAnaType == kSSCalo)
        {
          fhMCNLocMaxPt[i] = new TH2F
          (Form("hNLocMaxPt_MC%s",pname[i].Data()),
           Form("cluster from %s, pT of cluster vs NLM, accepted",ptype[i].Data()),
           nptbins,ptmin,ptmax,10,0,10);
          fhMCNLocMaxPt[i] ->SetYTitle("N maxima");
          fhMCNLocMaxPt[i] ->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhMCNLocMaxPt[i]) ;
 
          fhMCNLocMaxPtReject[i] = new TH2F
          (Form("hNLocMaxPtReject_MC%s",pname[i].Data()),
           Form("cluster from %s, pT of cluster vs NLM, rejected",ptype[i].Data()),
           nptbins,ptmin,ptmax,10,0,10);
          fhMCNLocMaxPtReject[i] ->SetYTitle("N maxima");
          fhMCNLocMaxPtReject[i] ->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhMCNLocMaxPtReject[i]) ;
          
          fhMCEReject[i]  = new TH1F
          (Form("hEReject_MC%s",pname[i].Data()),
           Form("Rejected as #pi^{0} (#eta), cluster from %s",
                ptype[i].Data()),
           nptbins,ptmin,ptmax);
          fhMCEReject[i]->SetYTitle("N");
          fhMCEReject[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCEReject[i]) ;
          
          fhMCPtReject[i]  = new TH1F
          (Form("hPtReject_MC%s",pname[i].Data()),
           Form("Rejected as #pi^{0} (#eta), cluster from %s",
                ptype[i].Data()),
           nptbins,ptmin,ptmax);
          fhMCPtReject[i]->SetYTitle("N");
          fhMCPtReject[i]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhMCPtReject[i]) ;
        }
        
        fhMCPhi[i]  = new TH2F
        (Form("hPhi_MC%s",pname[i].Data()),
         Form("Identified as #pi^{0} (#eta), cluster from %s",ptype[i].Data()),
         nptbins,ptmin,ptmax,nphibins,phimin,phimax);
        fhMCPhi[i]->SetYTitle("#phi");
        fhMCPhi[i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCPhi[i]) ;
        
        fhMCEta[i]  = new TH2F
        (Form("hEta_MC%s",pname[i].Data()),
         Form("Identified as #pi^{0} (#eta), cluster from %s",
              ptype[i].Data()),nptbins,ptmin,ptmax,netabins,etamin,etamax);
        fhMCEta[i]->SetYTitle("#eta");
        fhMCEta[i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCEta[i]) ;
        
        fhMCMassPt[i]  = new TH2F
        (Form("hMassPt_MC%s",pname[i].Data()),
         Form("all pairs mass: p_{T} vs massfrom %s",ptype[i].Data()),
         nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
        fhMCMassPt[i]->SetYTitle("mass (GeV/c^{2})");
        fhMCMassPt[i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCMassPt[i]) ;
        
        fhMCSelectedMassPt[i]  = new TH2F
        (Form("hSelectedMassPt_MC%s",pname[i].Data()),
         Form("Selected #pi^{0} (#eta) pairs mass: p_{T} vs massfrom %s",ptype[i].Data()),
         nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
        fhMCSelectedMassPt[i]->SetYTitle("mass (GeV/c^{2})");
        fhMCSelectedMassPt[i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCSelectedMassPt[i]) ;
        
        if(fAnaType == kSSCalo)
        {
          fhMCMassPtNoOverlap[i]  = new TH2F
          (Form("hMassPtNoOverlap_MC%s",pname[i].Data()),
           Form("all pairs mass: p_{T} vs massfrom %s, no overlap",ptype[i].Data()),
           nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
          fhMCMassPt[i]->SetYTitle("mass (GeV/c^{2})");
          fhMCMassPt[i]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhMCMassPtNoOverlap[i]) ;
          
          fhMCSelectedMassPtNoOverlap[i]  = new TH2F
          (Form("hSelectedMassPtNoOverlap_MC%s",pname[i].Data()),
           Form("Selected #pi^{0} (#eta) pairs mass: p_{T} vs massfrom %s, no overlap",ptype[i].Data()),
           nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
          fhMCSelectedMassPtNoOverlap[i]->SetYTitle("mass (GeV/c^{2})");
          fhMCSelectedMassPtNoOverlap[i]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhMCSelectedMassPtNoOverlap[i]) ;
        }
        
        if( fFillSelectClHisto )
        {
          fhEMCLambda0[i]  = new TH2F(Form("hELambda0_MC%s",pname[i].Data()),
                                      Form("Selected pair, cluster from %s : E vs #lambda_{0}^{2}",ptype[i].Data()),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhEMCLambda0[i]->SetYTitle("#lambda_{0}^{2}");
          fhEMCLambda0[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCLambda0[i]) ;
          
          fhEMCLambda1[i]  = new TH2F(Form("hELambda1_MC%s",pname[i].Data()),
                                      Form("Selected pair, cluster from %s : E vs #lambda_{1}^{2}",ptype[i].Data()),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhEMCLambda1[i]->SetYTitle("#lambda_{1}^{2}");
          fhEMCLambda1[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCLambda1[i]) ;
          
          fhEMCDispersion[i]  = new TH2F(Form("hEDispersion_MC%s",pname[i].Data()),
                                         Form("Selected pair, cluster from %s : E vs dispersion^{2}",ptype[i].Data()),
                                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhEMCDispersion[i]->SetYTitle("D^{2}");
          fhEMCDispersion[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCDispersion[i]) ;
          
          if(fCalorimeter=="EMCAL")
          {
            fhEMCLambda0NoTRD[i]  = new TH2F(Form("hELambda0NoTRD_MC%s",pname[i].Data()),
                                             Form("Selected pair, cluster from %s : E vs #lambda_{0}^{2}, NoTRD",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhEMCLambda0NoTRD[i]->SetYTitle("#lambda_{0}^{2}");
            fhEMCLambda0NoTRD[i]->SetXTitle("E (GeV)");
            outputContainer->Add(fhEMCLambda0NoTRD[i]) ;
            
            if(!fFillOnlySimpleSSHisto)
            {
              fhMCEDispEta[i]  = new TH2F (Form("hEDispEtaE_MC%s",pname[i].Data()),
                                           Form("cluster from %s : #sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i} - <#eta>)^{2}/ #Sigma w_{i} vs E",ptype[i].Data()),
                                           nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
              fhMCEDispEta[i]->SetXTitle("E (GeV)");
              fhMCEDispEta[i]->SetYTitle("#sigma^{2}_{#eta #eta}");
              outputContainer->Add(fhMCEDispEta[i]);
              
              fhMCEDispPhi[i]  = new TH2F (Form("hEDispPhiE_MC%s",pname[i].Data()),
                                           Form("cluster from %s : #sigma^{2}_{#phi #phi} = #Sigma w_{i}(#phi_{i} - <#phi>)^{2} / #Sigma w_{i} vs E",ptype[i].Data()),
                                           nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
              fhMCEDispPhi[i]->SetXTitle("E (GeV)");
              fhMCEDispPhi[i]->SetYTitle("#sigma^{2}_{#phi #phi}");
              outputContainer->Add(fhMCEDispPhi[i]);
              
              fhMCESumEtaPhi[i]  = new TH2F (Form("hESumEtaPhiE_MC%s",pname[i].Data()),
                                             Form("cluster from %s : #delta^{2}_{#eta #phi} = #Sigma w_{i}(#phi_{i} #eta_{i} ) / #Sigma w_{i} - <#phi><#eta> vs E",ptype[i].Data()),
                                             nptbins,ptmin,ptmax, 2*ssbins,-ssmax,ssmax);
              fhMCESumEtaPhi[i]->SetXTitle("E (GeV)");
              fhMCESumEtaPhi[i]->SetYTitle("#delta^{2}_{#eta #phi}");
              outputContainer->Add(fhMCESumEtaPhi[i]);
              
              fhMCEDispEtaPhiDiff[i]  = new TH2F (Form("hEDispEtaPhiDiffE_MC%s",pname[i].Data()),
                                                  Form("cluster from %s : #sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta} vs E",ptype[i].Data()),
                                                  nptbins,ptmin,ptmax,200,-10,10);
              fhMCEDispEtaPhiDiff[i]->SetXTitle("E (GeV)");
              fhMCEDispEtaPhiDiff[i]->SetYTitle("#sigma^{2}_{#phi #phi}-#sigma^{2}_{#eta #eta}");
              outputContainer->Add(fhMCEDispEtaPhiDiff[i]);
              
              fhMCESphericity[i]  = new TH2F (Form("hESphericity_MC%s",pname[i].Data()),
                                              Form("cluster from %s : (#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi}) vs E",ptype[i].Data()),
                                              nptbins,ptmin,ptmax, 200,-1,1);
              fhMCESphericity[i]->SetXTitle("E (GeV)");
              fhMCESphericity[i]->SetYTitle("s = (#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi})");
              outputContainer->Add(fhMCESphericity[i]);
              
              for(Int_t ie = 0; ie < 7; ie++)
              {
                fhMCDispEtaDispPhi[ie][i] = new TH2F (Form("hMCDispEtaDispPhi_EBin%d_MC%s",ie,pname[i].Data()),
                                                      Form("cluster from %s : #sigma^{2}_{#phi #phi} vs #sigma^{2}_{#eta #eta} for %d < E < %d GeV",pname[i].Data(),bin[ie],bin[ie+1]),
                                                      ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
                fhMCDispEtaDispPhi[ie][i]->SetXTitle("#sigma^{2}_{#eta #eta}");
                fhMCDispEtaDispPhi[ie][i]->SetYTitle("#sigma^{2}_{#phi #phi}");
                outputContainer->Add(fhMCDispEtaDispPhi[ie][i]);
                
                fhMCLambda0DispEta[ie][i] = new TH2F (Form("hMCLambda0DispEta_EBin%d_MC%s",ie,pname[i].Data()),
                                                      Form("cluster from %s : #lambda^{2}_{0} vs #sigma^{2}_{#eta #eta} for %d < E < %d GeV",pname[i].Data(),bin[ie],bin[ie+1]),
                                                      ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
                fhMCLambda0DispEta[ie][i]->SetXTitle("#lambda^{2}_{0}");
                fhMCLambda0DispEta[ie][i]->SetYTitle("#sigma^{2}_{#phi #phi}");
                outputContainer->Add(fhMCLambda0DispEta[ie][i]);
                
                fhMCLambda0DispPhi[ie][i] = new TH2F (Form("hMCLambda0DispPhi_EBin%d_MC%s",ie,pname[i].Data()),
                                                      Form("cluster from %s :#lambda^{2}_{0} vs #sigma^{2}_{#phi #phi} for %d < E < %d GeV",pname[i].Data(),bin[ie],bin[ie+1]),
                                                      ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
                fhMCLambda0DispPhi[ie][i]->SetXTitle("#lambda^{2}_{0}");
                fhMCLambda0DispPhi[ie][i]->SetYTitle("#sigma^{2}_{#phi #phi}");
                outputContainer->Add(fhMCLambda0DispPhi[ie][i]);
                
              }
            }
          }
          
          fhEMCLambda0FracMaxCellCut[i]  = new TH2F(Form("hELambda0FracMaxCellCut_MC%s",pname[i].Data()),
                                                    Form("Selected pair, cluster from %s : E vs #lambda_{0}^{2}, Max cell fraction of energy < 0.5 ",ptype[i].Data()),
                                                    nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhEMCLambda0FracMaxCellCut[i]->SetYTitle("#lambda_{0}^{2}");
          fhEMCLambda0FracMaxCellCut[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCLambda0FracMaxCellCut[i]) ;
          
          fhEMCFracMaxCell[i]  = new TH2F(Form("hEFracMaxCell_MC%s",pname[i].Data()),
                                          Form("Selected pair, cluster from %s : E vs Max cell fraction of energy",ptype[i].Data()),
                                          nptbins,ptmin,ptmax,100,0,1);
          fhEMCFracMaxCell[i]->SetYTitle("Fraction");
          fhEMCFracMaxCell[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCFracMaxCell[i]) ;
          
        }//
      } // shower shape histo
      
    } //Not MC reader
  }//Histos with MC
  
  if(fAnaType==kSSCalo)
  {
    fhAsymmetry  = new TH2F ("hAsymmetry","A = ( E1 - E2 ) / ( E1 + E2 ) vs E",
                             nptbins,ptmin,ptmax, 200, -1,1);
    fhAsymmetry->SetXTitle("E (GeV)");
    fhAsymmetry->SetYTitle("A = ( E1 - E2 ) / ( E1 + E2 )");
    outputContainer->Add(fhAsymmetry);
    
    fhSelectedAsymmetry  = new TH2F ("hSelectedAsymmetry","A = ( E1 - E2 ) / ( E1 + E2 ) vs E",
                                     nptbins,ptmin,ptmax, 200, -1,1);
    fhSelectedAsymmetry->SetXTitle("E (GeV)");
    fhSelectedAsymmetry->SetYTitle("A = ( E1 - E2 ) / ( E1 + E2 )");
    outputContainer->Add(fhSelectedAsymmetry);
    
    fhSplitE  = new TH1F
    ("hSplitE","Selected #pi^{0} (#eta) pairs energy sum of split sub-clusters",nptbins,ptmin,ptmax);
    fhSplitE->SetYTitle("counts");
    fhSplitE->SetXTitle("E (GeV)");
    outputContainer->Add(fhSplitE) ;
    
    fhSplitPt  = new TH1F
    ("hSplitPt","Selected #pi^{0} (#eta) pairs pT sum of split sub-clusters",nptbins,ptmin,ptmax);
    fhSplitPt->SetYTitle("counts");
    fhSplitPt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhSplitPt) ;
    
    
    fhSplitPtPhi  = new TH2F
    ("hSplitPtPhi","Selected #pi^{0} (#eta) pairs: sum split sub-cluster p_{T} vs #phi",nptbins,ptmin,ptmax, nphibins,phimin,phimax);
    fhSplitPtPhi->SetYTitle("#phi (rad)");
    fhSplitPtPhi->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhSplitPtPhi) ;
    
    fhSplitPtEta  = new TH2F
    ("hSplitPtEta","Selected #pi^{0} (#eta) pairs: sum split sub-cluster p_{T} vs #eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhSplitPtEta->SetYTitle("#eta");
    fhSplitPtEta->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhSplitPtEta) ;
    
    
    fhNLocMaxSplitPt = new TH2F("hNLocMaxSplitPt","Number of local maxima in cluster",
                                nptbins,ptmin,ptmax,10,0,10);
    fhNLocMaxSplitPt ->SetYTitle("N maxima");
    fhNLocMaxSplitPt ->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhNLocMaxSplitPt) ;
    
    
    fhMassSplitPt  = new TH2F
    ("hMassSplitPt","all pairs mass: sum split sub-cluster p_{T} vs mass",
     nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
    fhMassSplitPt->SetYTitle("mass (GeV/c^{2})");
    fhMassSplitPt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMassSplitPt) ;
    
    fhSelectedMassSplitPt  = new TH2F
    ("hSelectedMassSplitPt","Selected #pi^{0} (#eta) pairs mass: sum split sub-cluster p_{T} vs mass",
     nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
    fhSelectedMassSplitPt->SetYTitle("mass (GeV/c^{2})");
    fhSelectedMassSplitPt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhSelectedMassSplitPt) ;
    
    if(IsDataMC())
    {
      fhMassSplitPtNoOverlap  = new TH2F
      ("hMassSplitPtNoOverlap","all pairs mass: sum split sub-cluster p_{T} vs mass, no overlap",
       nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
      fhMassSplitPtNoOverlap->SetYTitle("mass (GeV/c^{2})");
      fhMassSplitPtNoOverlap->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhMassSplitPtNoOverlap) ;
      
      fhSelectedMassSplitPtNoOverlap  = new TH2F
      ("hSelectedMassSplitPtNoOverlap","Selected #pi^{0} (#eta) pairs mass: sum split sub-cluster p_{T} vs mass, no overlap",
       nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
      fhSelectedMassSplitPtNoOverlap->SetYTitle("mass (GeV/c^{2})");
      fhSelectedMassSplitPtNoOverlap->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhSelectedMassSplitPtNoOverlap) ;

      
      fhMCPi0PtRecoPtPrim  = new TH2F
      ("hMCPi0PtRecoPtPrim","p_{T,reco} vs p_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0PtRecoPtPrim ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCPi0PtRecoPtPrim ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCPi0PtRecoPtPrim ) ;
      
      fhMCPi0PtRecoPtPrimNoOverlap  = new TH2F
      ("hMCPi0PtRecoPtPrimNoOverlap","p_{T,reco} vs p_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0PtRecoPtPrimNoOverlap ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCPi0PtRecoPtPrimNoOverlap ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCPi0PtRecoPtPrimNoOverlap ) ;
      
      fhMCPi0SelectedPtRecoPtPrim  = new TH2F
      ("hMCPi0SelectedPtRecoPtPrim","p_{T,reco} vs p_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0SelectedPtRecoPtPrim ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCPi0SelectedPtRecoPtPrim ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCPi0SelectedPtRecoPtPrim ) ;
      
      fhMCPi0SelectedPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCPi0SelectedPtRecoPtPrimNoOverlap","p_{T,reco} vs p_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0SelectedPtRecoPtPrimNoOverlap ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCPi0SelectedPtRecoPtPrimNoOverlap ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCPi0SelectedPtRecoPtPrimNoOverlap ) ;

      
      fhMCPi0SplitPtRecoPtPrim  = new TH2F
      ("hMCPi0SplitPtRecoPtPrim","p_{T,reco} (split sum) vs p_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0SplitPtRecoPtPrim ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCPi0SplitPtRecoPtPrim ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCPi0SplitPtRecoPtPrim ) ;
      
      fhMCPi0SplitPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCPi0SplitPtRecoPtPrimNoOverlap","p_{T,reco} (split sum) vs p_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0SplitPtRecoPtPrimNoOverlap ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCPi0SplitPtRecoPtPrimNoOverlap ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCPi0SplitPtRecoPtPrimNoOverlap ) ;
      
      fhMCPi0SelectedSplitPtRecoPtPrim  = new TH2F
      ("hMCPi0SelectedSplitPtRecoPtPrim","p_{T,reco} (split sum) vs p_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0SelectedSplitPtRecoPtPrim ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCPi0SelectedSplitPtRecoPtPrim ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCPi0SelectedSplitPtRecoPtPrim ) ;
      
      fhMCPi0SelectedSplitPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCPi0SelectedSplitPtRecoPtPrimNoOverlap","p_{T,reco} (split sum) vs p_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCPi0SelectedSplitPtRecoPtPrimNoOverlap ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCPi0SelectedSplitPtRecoPtPrimNoOverlap ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCPi0SelectedSplitPtRecoPtPrimNoOverlap ) ;

      fhMCEtaPtRecoPtPrim  = new TH2F
      ("hMCEtaPtRecoPtPrim","p_{T,reco} vs p_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaPtRecoPtPrim ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCEtaPtRecoPtPrim ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCEtaPtRecoPtPrim ) ;
      
      fhMCEtaPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCEtaPtRecoPtPrimNoOverlap","p_{T,reco} vs p_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaPtRecoPtPrimNoOverlap ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCEtaPtRecoPtPrimNoOverlap ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCEtaPtRecoPtPrimNoOverlap ) ;
      
      fhMCEtaSelectedPtRecoPtPrim  = new TH2F
      ("hMCEtaSelectedPtRecoPtPrim","p_{T,reco} vs p_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaSelectedPtRecoPtPrim ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCEtaSelectedPtRecoPtPrim ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCEtaSelectedPtRecoPtPrim ) ;
      
      fhMCEtaSelectedPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCEtaSelectedPtRecoPtPrimNoOverlap","p_{T,reco} vs p_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaSelectedPtRecoPtPrimNoOverlap ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCEtaSelectedPtRecoPtPrimNoOverlap ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCEtaSelectedPtRecoPtPrimNoOverlap ) ;
      
      
      fhMCEtaSplitPtRecoPtPrim  = new TH2F
      ("hMCEtaSplitPtRecoPtPrim","p_{T,reco} (split sum) vs p_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaSplitPtRecoPtPrim ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCEtaSplitPtRecoPtPrim ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCEtaSplitPtRecoPtPrim ) ;
      
      fhMCEtaSplitPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCEtaSplitPtRecoPtPrimNoOverlap","p_{T,reco} (split sum) vs p_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaSplitPtRecoPtPrimNoOverlap ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCEtaSplitPtRecoPtPrimNoOverlap ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCEtaSplitPtRecoPtPrimNoOverlap ) ;
      
      fhMCEtaSelectedSplitPtRecoPtPrim  = new TH2F
      ("hMCEtaSelectedSplitPtRecoPtPrim","p_{T,reco} (split sum) vs p_{T,gen}",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaSelectedSplitPtRecoPtPrim ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCEtaSelectedSplitPtRecoPtPrim ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCEtaSelectedSplitPtRecoPtPrim ) ;
      
      fhMCEtaSelectedSplitPtRecoPtPrimNoOverlap  = new TH2F
      ("hMCEtaSelectedSplitPtRecoPtPrimNoOverlap","p_{T,reco} (split sum) vs p_{T,gen}, no overlap",
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMCEtaSelectedSplitPtRecoPtPrimNoOverlap ->SetYTitle("p_{T,gen} (GeV/c)");
      fhMCEtaSelectedSplitPtRecoPtPrimNoOverlap ->SetXTitle("p_{T,reco} (GeV/c)");
      outputContainer->Add(fhMCEtaSelectedSplitPtRecoPtPrimNoOverlap ) ;
      
      for(Int_t i = 0; i< 6; i++)
      {
        fhMCEAsymmetry[i]  = new TH2F (Form("hEAsymmetry_MC%s",pname[i].Data()),
                                       Form("cluster from %s : A = ( E1 - E2 ) / ( E1 + E2 ) vs E",ptype[i].Data()),
                                       nptbins,ptmin,ptmax, 200,-1,1);
        fhMCEAsymmetry[i]->SetXTitle("E (GeV)");
        fhMCEAsymmetry[i]->SetYTitle("A = ( E1 - E2 ) / ( E1 + E2 )");
        outputContainer->Add(fhMCEAsymmetry[i]);
        
        fhMCSplitE[i]  = new TH1F
        (Form("hSplitE_MC%s",pname[i].Data()),
         Form("cluster from %s, energy sum of split sub-clusters",ptype[i].Data()),
         nptbins,ptmin,ptmax);
        fhMCSplitE[i]->SetYTitle("counts");
        fhMCSplitE[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCSplitE[i]) ;
        
        fhMCSplitPt[i]  = new TH1F
        (Form("hSplitPt_MC%s",pname[i].Data()),
         Form("cluster from %s, pT sum of split sub-clusters",ptype[i].Data()),
         nptbins,ptmin,ptmax);
        fhMCSplitPt[i]->SetYTitle("counts");
        fhMCSplitPt[i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCSplitPt[i]) ;
        
        
        fhMCSplitPtPhi[i]  = new TH2F
        (Form("hSplitPtPhi_MC%s",pname[i].Data()),
         Form("Identified as #pi^{0} (#eta), cluster from %s",ptype[i].Data()),
         nptbins,ptmin,ptmax,nphibins,phimin,phimax);
        fhMCSplitPtPhi[i]->SetYTitle("#phi");
        fhMCSplitPtPhi[i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCSplitPtPhi[i]) ;
        
        fhMCSplitPtEta[i]  = new TH2F
        (Form("hSplitPtEta_MC%s",pname[i].Data()),
         Form("Identified as #pi^{0} (#eta), cluster from %s",
              ptype[i].Data()),nptbins,ptmin,ptmax,netabins,etamin,etamax);
        fhMCSplitPtEta[i]->SetYTitle("#eta");
        fhMCSplitPtEta[i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCSplitPtEta[i]) ;
        
        
        fhMCNLocMaxSplitPt[i] = new TH2F
        (Form("hNLocMaxSplitPt_MC%s",pname[i].Data()),
         Form("cluster from %s, pT sum of split sub-clusters, for NLM",ptype[i].Data()),
         nptbins,ptmin,ptmax,10,0,10);
        fhMCNLocMaxSplitPt[i] ->SetYTitle("N maxima");
        fhMCNLocMaxSplitPt[i] ->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCNLocMaxSplitPt[i]) ;
        
        fhMCMassSplitPt[i]  = new TH2F
        (Form("hMassSplitPt_MC%s",pname[i].Data()),
         Form("all pairs mass: split p_{T} vs mass from %s",ptype[i].Data()),
         nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
        fhMCMassSplitPt[i]->SetYTitle("mass (GeV/c^{2})");
        fhMCMassSplitPt[i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCMassSplitPt[i]) ;
        
        fhMCSelectedMassSplitPt[i]  = new TH2F
        (Form("hSelectedMassSplitPt_MC%s",pname[i].Data()),
         Form("Selected #pi^{0} (#eta) pairs mass: split p_{T} vs mass from %s",ptype[i].Data()),
         nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
        fhMCSelectedMassSplitPt[i]->SetYTitle("mass (GeV/c^{2})");
        fhMCSelectedMassSplitPt[i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCSelectedMassSplitPt[i]) ;

        fhMCMassSplitPtNoOverlap[i]  = new TH2F
        (Form("hMassSplitPtNoOverlap_MC%s",pname[i].Data()),
         Form("all pairs mass: split p_{T} vs mass from %s, no overlap",ptype[i].Data()),
         nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
        fhMCMassSplitPtNoOverlap[i]->SetYTitle("mass (GeV/c^{2})");
        fhMCMassSplitPtNoOverlap[i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCMassSplitPtNoOverlap[i]) ;
        
        fhMCSelectedMassSplitPtNoOverlap[i]  = new TH2F
        (Form("hSelectedMassSplitPtNoOverlap_MC%s",pname[i].Data()),
         Form("Selected #pi^{0} (#eta) pairs mass: split p_{T} vs mass from %s, no overlap",ptype[i].Data()),
         nptbins,ptmin,ptmax, nmassbins,massmin,massmax);
        fhMCSelectedMassSplitPtNoOverlap[i]->SetYTitle("mass (GeV/c^{2})");
        fhMCSelectedMassSplitPtNoOverlap[i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCSelectedMassSplitPtNoOverlap[i]) ;
      }
    }
  }
  
  if(fAnaType==kSSCalo && fFillSelectClHisto && !fFillOnlySimpleSSHisto )
  {
    
    
    for(Int_t i = 0; i< 3; i++)
    {
      fhEAsymmetryLocMax[i]  = new TH2F(Form("hEAsymmetryLocMax%d",i+1),
                                        Form("Selected #pi^{0} (#eta) pairs: E vs A = ( E1 - E2 ) / ( E1 + E2 ), %s",nlm[i].Data()),
                                        nptbins,ptmin,ptmax,200, -1,1);
      fhEAsymmetryLocMax[i]->SetYTitle("A = ( E1 - E2 ) / ( E1 + E2 )");
      fhEAsymmetryLocMax[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhEAsymmetryLocMax[i]) ;
    }
    
    for(Int_t ie = 0; ie< 7; ie++)
    {
      
      fhAsymmetryLambda0[ie] = new TH2F (Form("hAsymmetryLambda0_EBin%d",ie),
                                         Form("#lambda_{0}^{2} vs A for %d < E < %d GeV",bin[ie],bin[ie+1]),
                                         ssbins,ssmin,ssmax , 200,-1,1);
      fhAsymmetryLambda0[ie]->SetXTitle("#lambda_{0}^{2}");
      fhAsymmetryLambda0[ie]->SetYTitle("A = ( E1 - E2 ) / ( E1 + E2 )");
      outputContainer->Add(fhAsymmetryLambda0[ie]);
      
      fhAsymmetryDispEta[ie] = new TH2F (Form("hAsymmetryDispEta_EBin%d",ie),
                                         Form("#sigma^{2}_{#eta #eta} vs A for %d < E < %d GeV",bin[ie],bin[ie+1]),
                                         ssbins,ssmin,ssmax , 200,-1,1);
      fhAsymmetryDispEta[ie]->SetXTitle("#sigma^{2}_{#eta #eta}");
      fhAsymmetryDispEta[ie]->SetYTitle("A = ( E1 - E2 ) / ( E1 + E2 )");
      outputContainer->Add(fhAsymmetryDispEta[ie]);
      
      fhAsymmetryDispPhi[ie] = new TH2F (Form("hAsymmetryDispPhi_EBin%d",ie),
                                         Form("#sigma^{2}_{#phi #phi} vs A for %d < E < %d GeV",bin[ie],bin[ie+1]),
                                         ssbins,ssmin,ssmax , 200,-1,1);
      fhAsymmetryDispPhi[ie]->SetXTitle("#sigma^{2}_{#phi #phi}");
      fhAsymmetryDispPhi[ie]->SetYTitle("A = ( E1 - E2 ) / ( E1 + E2 )");
      outputContainer->Add(fhAsymmetryDispPhi[ie]);
    }
    
    
    if(IsDataMC())
    {
      for(Int_t i = 0; i< 6; i++)
      {
        for(Int_t ie = 0; ie < 7; ie++)
        {
          fhMCAsymmetryLambda0[ie][i] = new TH2F (Form("hMCAsymmetryLambda0_EBin%d_MC%s",ie,pname[i].Data()),
                                                  Form("cluster from %s : #lambda_{0}^{2} vs A for %d < E < %d GeV",pname[i].Data(),bin[ie],bin[ie+1]),
                                                  ssbins,ssmin,ssmax , 200,-1,1);
          fhMCAsymmetryLambda0[ie][i]->SetXTitle("#lambda_{0}^{2}");
          fhMCAsymmetryLambda0[ie][i]->SetYTitle("A = ( E1 - E2 ) / ( E1 + E2 )");
          outputContainer->Add(fhMCAsymmetryLambda0[ie][i]);
          
          fhMCAsymmetryDispEta[ie][i] = new TH2F (Form("hMCAsymmetryDispEta_EBin%d_MC%s",ie,pname[i].Data()),
                                                  Form("cluster from %s : #sigma^{2}_{#eta #eta} vs A for %d < E < %d GeV",pname[i].Data(),bin[ie],bin[ie+1]),
                                                  ssbins,ssmin,ssmax , 200,-1,1);
          fhMCAsymmetryDispEta[ie][i]->SetXTitle("#sigma^{2}_{#eta #eta}");
          fhMCAsymmetryDispEta[ie][i]->SetYTitle("A = ( E1 - E2 ) / ( E1 + E2 )");
          outputContainer->Add(fhMCAsymmetryDispEta[ie][i]);
          
          fhMCAsymmetryDispPhi[ie][i] = new TH2F (Form("hMCAsymmetryDispPhi_EBin%d_MC%s",ie,pname[i].Data()),
                                                  Form("cluster from %s : #sigma^{2}_{#phi #phi} vs A for %d < E < %d GeV",pname[i].Data(),bin[ie],bin[ie+1]),
                                                  ssbins,ssmin,ssmax , 200,-1,1);
          fhMCAsymmetryDispPhi[ie][i]->SetXTitle("#sigma^{2}_{#phi #phi}");
          fhMCAsymmetryDispPhi[ie][i]->SetYTitle("A = ( E1 - E2 ) / ( E1 + E2 )");
          outputContainer->Add(fhMCAsymmetryDispPhi[ie][i]);
        }
      }
    }
  }
  
  if(fFillPileUpHistograms)
  {
    
    TString pileUpName[] = {"SPD","EMCAL","SPDOrEMCAL","SPDAndEMCAL","SPDAndNotEMCAL","EMCALAndNotSPD","NotSPDAndNotEMCAL"} ;
    
    for(Int_t i = 0 ; i < 7 ; i++)
    {
      fhPtPileUp[i]  = new TH1F(Form("hPtPileUp%s",pileUpName[i].Data()),
                                   Form("Selected #pi^{0} (#eta) p_{T} distribution, %s Pile-Up event",pileUpName[i].Data()), nptbins,ptmin,ptmax);
      fhPtPileUp[i]->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtPileUp[i]);
      
      fhPtCellTimePileUp[i]  = new TH2F(Form("hPtCellTimePileUp%s",pileUpName[i].Data()),
                                             Form("Pt vs cell time in cluster, %s Pile-Up event",pileUpName[i].Data()),
                                             nptbins,ptmin,ptmax,ntimebins,timemin,timemax);
      fhPtCellTimePileUp[i]->SetXTitle("p_{T} (GeV/c)");
      fhPtCellTimePileUp[i]->SetYTitle("t_{cell} (ns)");
      outputContainer->Add(fhPtCellTimePileUp[i]);
      
      fhPtTimeDiffPileUp[i]  = new TH2F(Form("hPtTimeDiffPileUp%s",pileUpName[i].Data()),
                                             Form("Pt vs t_{max}-t_{cell} in cluster, %s Pile-Up event",pileUpName[i].Data()),
                                             nptbins,ptmin,ptmax,400,-200,200);
      fhPtTimeDiffPileUp[i]->SetXTitle("p_{T} (GeV/c");
      fhPtTimeDiffPileUp[i]->SetYTitle("t_{max}-t_{cell} (ns)");
      outputContainer->Add(fhPtTimeDiffPileUp[i]);

    }
    
    fhTimePtNoCut  = new TH2F ("hTimePt_NoCut","time of cluster vs E of clusters, no cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimePtNoCut->SetXTitle("p_{T} (GeV/c)");
    fhTimePtNoCut->SetYTitle("time (ns)");
    outputContainer->Add(fhTimePtNoCut);
    
    fhTimePtSPD  = new TH2F ("hTimePt_SPD","time of cluster vs E of clusters, SPD cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimePtSPD->SetXTitle("p_{T} (GeV/c)");
    fhTimePtSPD->SetYTitle("time (ns)");
    outputContainer->Add(fhTimePtSPD);
    
    fhTimePtSPDMulti  = new TH2F ("hTimePt_SPDMulti","time of cluster vs E of clusters, SPD multi cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimePtSPDMulti->SetXTitle("p_{T} (GeV/c)");
    fhTimePtSPDMulti->SetYTitle("time (ns)");
    outputContainer->Add(fhTimePtSPDMulti);
    
    fhTimeNPileUpVertSPD  = new TH2F ("hTime_NPileUpVertSPD","time of cluster vs N pile-up SPD vertex", ntimebins,timemin,timemax,50,0,50);
    fhTimeNPileUpVertSPD->SetYTitle("# vertex ");
    fhTimeNPileUpVertSPD->SetXTitle("time (ns)");
    outputContainer->Add(fhTimeNPileUpVertSPD);
    
    fhTimeNPileUpVertTrack  = new TH2F ("hTime_NPileUpVertTracks","time of cluster vs N pile-up Tracks vertex", ntimebins,timemin,timemax, 50,0,50 );
    fhTimeNPileUpVertTrack->SetYTitle("# vertex ");
    fhTimeNPileUpVertTrack->SetXTitle("time (ns)");
    outputContainer->Add(fhTimeNPileUpVertTrack);
    
    fhTimeNPileUpVertContributors  = new TH2F ("hTime_NPileUpVertContributors","time of cluster vs N constributors to pile-up SPD vertex", ntimebins,timemin,timemax,50,0,50);
    fhTimeNPileUpVertContributors->SetYTitle("# vertex ");
    fhTimeNPileUpVertContributors->SetXTitle("time (ns)");
    outputContainer->Add(fhTimeNPileUpVertContributors);
    
    fhTimePileUpMainVertexZDistance  = new TH2F ("hTime_PileUpMainVertexZDistance","time of cluster vs distance in Z pile-up SPD vertex - main SPD vertex",ntimebins,timemin,timemax,100,0,50);
    fhTimePileUpMainVertexZDistance->SetYTitle("distance Z (cm) ");
    fhTimePileUpMainVertexZDistance->SetXTitle("time (ns)");
    outputContainer->Add(fhTimePileUpMainVertexZDistance);
    
    fhTimePileUpMainVertexZDiamond  = new TH2F ("hTime_PileUpMainVertexZDiamond","time of cluster vs distance in Z pile-up SPD vertex - z diamond",ntimebins,timemin,timemax,100,0,50);
    fhTimePileUpMainVertexZDiamond->SetYTitle("diamond distance Z (cm) ");
    fhTimePileUpMainVertexZDiamond->SetXTitle("time (ns)");
    outputContainer->Add(fhTimePileUpMainVertexZDiamond);
		
		fhPtNPileUpSPDVtx  = new TH2F ("hPt_NPileUpVertSPD","pT of cluster vs N pile-up SPD vertex",
																	 nptbins,ptmin,ptmax,20,0,20);
		fhPtNPileUpSPDVtx->SetYTitle("# vertex ");
		fhPtNPileUpSPDVtx->SetXTitle("p_{T} (GeV/c)");
		outputContainer->Add(fhPtNPileUpSPDVtx);
	  
		fhPtNPileUpTrkVtx  = new TH2F ("hPt_NPileUpVertTracks","pT of cluster vs N pile-up Tracks vertex",
																	 nptbins,ptmin,ptmax, 20,0,20 );
		fhPtNPileUpTrkVtx->SetYTitle("# vertex ");
		fhPtNPileUpTrkVtx->SetXTitle("p_{T} (GeV/c)");
		outputContainer->Add(fhPtNPileUpTrkVtx);
		
		fhPtNPileUpSPDVtxTimeCut  = new TH2F ("hPt_NPileUpVertSPD_TimeCut","pT of cluster vs N pile-up SPD vertex, |tof| < 25 ns",
                                          nptbins,ptmin,ptmax,20,0,20);
		fhPtNPileUpSPDVtxTimeCut->SetYTitle("# vertex ");
		fhPtNPileUpSPDVtxTimeCut->SetXTitle("p_{T} (GeV/c)");
		outputContainer->Add(fhPtNPileUpSPDVtxTimeCut);
	  
		fhPtNPileUpTrkVtxTimeCut  = new TH2F ("hPt_NPileUpVertTracks_TimeCut","pT of cluster vs N pile-up Tracks vertex, |tof| < 25 ns",
																					nptbins,ptmin,ptmax, 20,0,20 );
		fhPtNPileUpTrkVtxTimeCut->SetYTitle("# vertex ");
		fhPtNPileUpTrkVtxTimeCut->SetXTitle("p_{T} (GeV/c)");
		outputContainer->Add(fhPtNPileUpTrkVtxTimeCut);
    
    fhPtNPileUpSPDVtxTimeCut2  = new TH2F ("hPt_NPileUpVertSPD_TimeCut2","pT of cluster vs N pile-up SPD vertex, -25 < tof < 75 ns",
                                           nptbins,ptmin,ptmax,20,0,20);
		fhPtNPileUpSPDVtxTimeCut2->SetYTitle("# vertex ");
		fhPtNPileUpSPDVtxTimeCut2->SetXTitle("p_{T} (GeV/c)");
		outputContainer->Add(fhPtNPileUpSPDVtxTimeCut2);
	  
		fhPtNPileUpTrkVtxTimeCut2  = new TH2F ("hPt_NPileUpVertTracks_TimeCut2","pT of cluster vs N pile-up Tracks vertex, -25 < tof < 75 ns",
                                           nptbins,ptmin,ptmax, 20,0,20 );
		fhPtNPileUpTrkVtxTimeCut2->SetYTitle("# vertex ");
		fhPtNPileUpTrkVtxTimeCut2->SetXTitle("p_{T} (GeV/c)");
		outputContainer->Add(fhPtNPileUpTrkVtxTimeCut2);
    
  }
  
  //Keep neutral meson selection histograms if requiered
  //Setting done in AliNeutralMesonSelection
  
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
Int_t AliAnaPi0EbE::GetMCIndex(const Int_t tag)
{
  
  // Assign mc index depending on MC bit set
  
  if       ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)  )
  {
    return kmcPi0 ;
  }//pi0
  else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)  )
  {
    return kmcEta ;
  }//eta
  else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) &&
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion) )
  {
    return kmcConversion ;
  }//conversion photon
  else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) )
  {
    return kmcPhoton ;
  }//photon   no conversion
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
void AliAnaPi0EbE::HasPairSameMCMother(AliAODPWG4Particle * photon1,
                                       AliAODPWG4Particle * photon2,
                                       Int_t & label, Int_t & tag)
{
  // Check the labels of pare in case mother was same pi0 or eta
  // Set the new AOD accordingly
  
  Int_t  label1 = photon1->GetLabel();
  Int_t  label2 = photon2->GetLabel();
  
  if(label1 < 0 || label2 < 0 ) return ;
  
  //Int_t tag1 = GetMCAnalysisUtils()->CheckOrigin(label1, GetReader());
  //Int_t tag2 = GetMCAnalysisUtils()->CheckOrigin(label2, GetReader());
  Int_t tag1 = photon1->GetTag();
  Int_t tag2 = photon2->GetTag();
  
  if(GetDebug() > 0) printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - Origin of: photon1 %d; photon2 %d \n",tag1, tag2);
  if( (GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCPi0Decay) &&
       GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCPi0Decay)    ) ||
     (GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCEtaDecay) &&
      GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCEtaDecay)    )
     )
  {
    
    //Check if pi0/eta mother is the same
    if(GetReader()->ReadStack())
    {
      if(label1>=0)
      {
        TParticle * mother1 = GetMCStack()->Particle(label1);//photon in kine tree
        label1 = mother1->GetFirstMother();
        //mother1 = GetMCStack()->Particle(label1);//pi0
      }
      if(label2>=0)
      {
        TParticle * mother2 = GetMCStack()->Particle(label2);//photon in kine tree
        label2 = mother2->GetFirstMother();
        //mother2 = GetMCStack()->Particle(label2);//pi0
      }
    } // STACK
    else if(GetReader()->ReadAODMCParticles())
    {//&& (input > -1)){
      if(label1>=0)
      {
        AliAODMCParticle * mother1 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles())->At(label1);//photon in kine tree
        label1 = mother1->GetMother();
        //mother1 = GetMCStack()->Particle(label1);//pi0
      }
      if(label2>=0)
      {
        AliAODMCParticle * mother2 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles())->At(label2);//photon in kine tree
        label2 = mother2->GetMother();
        //mother2 = GetMCStack()->Particle(label2);//pi0
      }
    }// AOD
    
    //printf("mother1 %d, mother2 %d\n",label1,label2);
    if( label1 == label2 && label1>=0 )
    {
      
      label = label1;
      
      TLorentzVector mom1 = *(photon1->Momentum());
      TLorentzVector mom2 = *(photon2->Momentum());
      
      Double_t angle = mom2.Angle(mom1.Vect());
      Double_t mass  = (mom1+mom2).M();
      Double_t epair = (mom1+mom2).E();
      
      if(GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCPi0Decay))
      {
        fhMassPairMCPi0 ->Fill(epair,mass);
        fhAnglePairMCPi0->Fill(epair,angle);
        GetMCAnalysisUtils()->SetTagBit(tag,AliMCAnalysisUtils::kMCPi0);
      }
      else
      {
        fhMassPairMCEta ->Fill(epair,mass);
        fhAnglePairMCEta->Fill(epair,angle);
        GetMCAnalysisUtils()->SetTagBit(tag,AliMCAnalysisUtils::kMCEta);
      }
      
    } // same label
  } // both from eta or pi0 decay
  
}

//____________________________________________________________________________
void AliAnaPi0EbE::Init()
{
  //Init
  //Do some checks
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn() && NewOutputAOD()){
    printf("AliAnaPi0EbE::Init() - !!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD()){
    printf("AliAnaPi0EbE::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
}

//____________________________________________________________________________
void AliAnaPi0EbE::InitParameters()
{
  //Initialize the parameters of the analysis.
  AddToHistogramsName("AnaPi0EbE_");
  
  fInputAODGammaConvName = "PhotonsCTS" ;
  fAnaType = kIMCalo ;
  fCalorimeter = "EMCAL" ;
  fMinDist  = 2.;
  fMinDist2 = 4.;
  fMinDist3 = 5.;
  
  fNLMECutMin[0] = 10.;
  fNLMECutMin[1] = 6. ;
  fNLMECutMin[2] = 6. ;
  
}

//__________________________________________________________________
void  AliAnaPi0EbE::MakeAnalysisFillAOD()
{
  //Do analysis and fill aods
  
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
}

//____________________________________________
void  AliAnaPi0EbE::MakeInvMassInCalorimeter()
{
  //Do analysis and fill aods
  //Search for the photon decay in calorimeters
  //Read photon list from AOD, produced in class AliAnaPhoton
  //Check if 2 photons have the mass of the pi0.
  
  TLorentzVector mom1;
  TLorentzVector mom2;
  TLorentzVector mom ;
  
  Int_t tag   = 0;
  Int_t label = 0;
  
  if(!GetInputAODBranch()){
    printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - No input calo photons in AOD with name branch < %s >, STOP \n",GetInputAODName().Data());
    abort();
  }
  
  //Get shower shape information of clusters
  TObjArray *clusters = 0;
  if     (fCalorimeter=="EMCAL") clusters = GetEMCALClusters();
  else if(fCalorimeter=="PHOS")  clusters = GetPHOSClusters() ;
  
  for(Int_t iphoton = 0; iphoton < GetInputAODBranch()->GetEntriesFast()-1; iphoton++){
    AliAODPWG4Particle * photon1 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(iphoton));
    
    //Vertex cut in case of mixed events
    Int_t evtIndex1 = 0 ;
    if(GetMixedEvent())
      evtIndex1 = GetMixedEvent()->EventIndexForCaloCluster(photon1->GetCaloLabel(0)) ;
    if(TMath::Abs(GetVertex(evtIndex1)[2]) > GetZvertexCut()) continue ;  //vertex cut
    mom1 = *(photon1->Momentum());
    
    //Get original cluster, to recover some information
    Int_t iclus = -1;
    AliVCluster *cluster1 = FindCluster(clusters,photon1->GetCaloLabel(0),iclus);
    
    if(!cluster1){
      printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - First cluster not found\n");
      return;
    }
    
    for(Int_t jphoton = iphoton+1; jphoton < GetInputAODBranch()->GetEntriesFast(); jphoton++)
    {
      AliAODPWG4Particle * photon2 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(jphoton));
      
      Int_t evtIndex2 = 0 ;
      if(GetMixedEvent())
        evtIndex2 = GetMixedEvent()->EventIndexForCaloCluster(photon2->GetCaloLabel(0)) ;
      
      if(GetMixedEvent() && (evtIndex1 == evtIndex2))
        continue ;
      
      if(TMath::Abs(GetVertex(evtIndex2)[2]) > GetZvertexCut()) continue ;  //vertex cut
      
      mom2 = *(photon2->Momentum());
      
      //Get original cluster, to recover some information
      Int_t iclus2;
      AliVCluster *cluster2 = FindCluster(clusters,photon2->GetCaloLabel(0),iclus2,iclus+1);
      
      if(!cluster2)
      {
        printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - Second cluster not found\n");
        return;
      }
      
      Float_t e1    = photon1->E();
      Float_t e2    = photon2->E();
      
      //Select clusters with good time window difference
      Float_t tof1  = cluster1->GetTOF()*1e9;;
      Float_t tof2  = cluster2->GetTOF()*1e9;;
      Double_t t12diff = tof1-tof2;
      fhEPairDiffTime->Fill(e1+e2,    t12diff);
      if(TMath::Abs(t12diff) > GetPairTimeCut()) continue;
      
      //Play with the MC stack if available
      if(IsDataMC()) HasPairSameMCMother(photon1, photon2, label, tag) ;
      
      // Check the invariant mass for different selection on the local maxima
      // Name of AOD method TO BE FIXED
      Int_t nMaxima1 = photon1->GetFiducialArea();
      Int_t nMaxima2 = photon2->GetFiducialArea();
      
      Double_t mass  = (mom1+mom2).M();
      Double_t epair = (mom1+mom2).E();
      
      if(nMaxima1==nMaxima2)
      {
        if     (nMaxima1==1) fhMassPairLocMax[0]->Fill(epair,mass);
        else if(nMaxima1==2) fhMassPairLocMax[1]->Fill(epair,mass);
        else                 fhMassPairLocMax[2]->Fill(epair,mass);
      }
      else if(nMaxima1==1 || nMaxima2==1)
      {
        if  (nMaxima1==2 || nMaxima2==2) fhMassPairLocMax[3]->Fill(epair,mass);
        else                             fhMassPairLocMax[4]->Fill(epair,mass);
      }
      else
        fhMassPairLocMax[5]->Fill(epair,mass);
      
      // combinations with SS axis cut and NLM cut
      if(nMaxima1 == 1 && cluster2->GetM02() > 0.3) fhMassPairLocMax[6]->Fill(epair,mass);
      if(nMaxima2 == 1 && cluster1->GetM02() > 0.3) fhMassPairLocMax[6]->Fill(epair,mass);
      if(nMaxima1 >  1 && cluster2->GetM02() < 0.3 && cluster2->GetM02()> 0.1 ) fhMassPairLocMax[7]->Fill(epair,mass);
      if(nMaxima2 >  1 && cluster1->GetM02() < 0.3 && cluster1->GetM02()> 0.1 ) fhMassPairLocMax[7]->Fill(epair,mass);
      
      //Skip events with too few or too many  NLM
      if((nMaxima1 < fNLMCutMin || nMaxima1 > fNLMCutMax) || (nMaxima2 < fNLMCutMin || nMaxima2 > fNLMCutMax)) continue ;
      
      if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - NLM of out of range: cluster1 %d, cluster2 %d \n",nMaxima1, nMaxima2);
      
      //Mass of all pairs
      fhMass->Fill(epair,(mom1+mom2).M());
      
      //Select good pair (good phi, pt cuts, aperture and invariant mass)
      if(GetNeutralMesonSelection()->SelectPair(mom1, mom2,fCalorimeter))
      {
        if(GetDebug()>1)
          printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - Selected gamma pair: pt %f, phi %f, eta%f \n",(mom1+mom2).Pt(), (mom1+mom2).Phi()*180./3.1416, (mom1+mom2).Eta());
        
        //Fill some histograms about shower shape
        if(fFillSelectClHisto && clusters && GetReader()->GetDataType()!=AliCaloTrackReader::kMC)
        {
          FillSelectedClusterHistograms(cluster1, nMaxima1, photon1->GetTag());
          FillSelectedClusterHistograms(cluster2, nMaxima2, photon2->GetTag());
        }
        
        // Tag both photons as decay
        photon1->SetTagged(kTRUE);
        photon2->SetTagged(kTRUE);
        
        fhPtDecay->Fill(photon1->Pt());
        fhEDecay ->Fill(photon1->E() );
        
        fhPtDecay->Fill(photon2->Pt());
        fhEDecay ->Fill(photon2->E() );
        
        //Create AOD for analysis
        mom = mom1+mom2;
        
        //Mass of selected pairs
        fhSelectedMass->Fill(epair,mom.M());
        
        // Fill histograms to undertand pile-up before other cuts applied
        // Remember to relax time cuts in the reader
        FillPileUpHistograms(mom.Pt(),((cluster1->GetTOF()+cluster2->GetTOF())*1e9)/2,cluster1);
        
        AliAODPWG4Particle pi0 = AliAODPWG4Particle(mom);
        
        pi0.SetIdentifiedParticleType(AliCaloPID::kPi0);
        pi0.SetDetector(photon1->GetDetector());
        
        // MC
        pi0.SetLabel(label);
        pi0.SetTag(tag);
        
        //Set the indeces of the original caloclusters
        pi0.SetCaloLabel(photon1->GetCaloLabel(0), photon2->GetCaloLabel(0));
        //pi0.SetInputFileIndex(input);
        
        AddAODParticle(pi0);
        
      }//pi0
      
    }//2n photon loop
    
  }//1st photon loop
  
  if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - End fill AODs \n");
  
}

//__________________________________________________
void  AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS()
{
  //Do analysis and fill aods
  //Search for the photon decay in calorimeters
  //Read photon list from AOD, produced in class AliAnaPhoton and AliGammaConversion
  //Check if 2 photons have the mass of the pi0.
  
  TLorentzVector mom1;
  TLorentzVector mom2;
  TLorentzVector mom ;
  Int_t tag   = 0;
  Int_t label = 0;
  Int_t evtIndex = 0;
  
  // Check calorimeter input
  if(!GetInputAODBranch()){
    printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - No input calo photons in AOD branch with name < %s > , STOP\n",GetInputAODName().Data());
    abort();
  }
  
  // Get the array with conversion photons
  TClonesArray * inputAODGammaConv = (TClonesArray *) GetReader()->GetOutputEvent()->FindListObject(fInputAODGammaConvName);
  if(!inputAODGammaConv) {
    
    inputAODGammaConv = (TClonesArray *) GetReader()->GetInputEvent()->FindListObject(fInputAODGammaConvName);
    
    if(!inputAODGammaConv) {
      printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - No input gamma conversions in AOD branch with name < %s >\n",fInputAODGammaConvName.Data());
      
      return;
    }
  }
  
  //Get shower shape information of clusters
  TObjArray *clusters = 0;
  if     (fCalorimeter=="EMCAL") clusters = GetEMCALClusters();
  else if(fCalorimeter=="PHOS")  clusters = GetPHOSClusters() ;
  
  Int_t nCTS  = inputAODGammaConv->GetEntriesFast();
  Int_t nCalo = GetInputAODBranch()->GetEntriesFast();
  if(nCTS<=0 || nCalo <=0)
  {
    if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - nCalo %d, nCTS %d, cannot loop\n",nCalo,nCTS);
    return;
  }
  
  if(GetDebug() > 1)
    printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - Number of conversion photons %d\n",nCTS);
  
  // Do the loop, first calo, second CTS
  for(Int_t iphoton = 0; iphoton < GetInputAODBranch()->GetEntriesFast(); iphoton++){
    AliAODPWG4Particle * photon1 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(iphoton));
    mom1 = *(photon1->Momentum());
    
    //Get original cluster, to recover some information
    Int_t iclus = -1;
    AliVCluster *cluster = FindCluster(clusters,photon1->GetCaloLabel(0),iclus);
    
    for(Int_t jphoton = 0; jphoton < nCTS; jphoton++){
      AliAODPWG4Particle * photon2 =  (AliAODPWG4Particle*) (inputAODGammaConv->At(jphoton));
      if(GetMixedEvent())
        evtIndex = GetMixedEvent()->EventIndexForCaloCluster(photon2->GetCaloLabel(0)) ;
      if(TMath::Abs(GetVertex(evtIndex)[2]) > GetZvertexCut()) continue ;  //vertex cut
      
      mom2 = *(photon2->Momentum());
      
      Double_t mass  = (mom1+mom2).M();
      Double_t epair = (mom1+mom2).E();
      
      Int_t nMaxima = photon1->GetFiducialArea();
      if     (nMaxima==1) fhMassPairLocMax[0]->Fill(epair,mass);
      else if(nMaxima==2) fhMassPairLocMax[1]->Fill(epair,mass);
      else                fhMassPairLocMax[2]->Fill(epair,mass);
      
      if(nMaxima < fNLMCutMin || nMaxima > fNLMCutMax) continue ;
      if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - NLM %d of out of range \n",nMaxima);
      
      //Play with the MC stack if available
      if(IsDataMC())
      {
        Int_t	label2 = photon2->GetLabel();
        if(label2 >= 0 )photon2->SetTag(GetMCAnalysisUtils()->CheckOrigin(label2, GetReader()));
        
        HasPairSameMCMother(photon1, photon2, label, tag) ;
      }
      
      //Mass of selected pairs
      fhMass->Fill(epair,(mom1+mom2).M());
      
      //Select good pair (good phi, pt cuts, aperture and invariant mass)
      if(GetNeutralMesonSelection()->SelectPair(mom1, mom2,fCalorimeter))
      {
        if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - Selected gamma pair: pt %f, phi %f, eta%f\n",(mom1+mom2).Pt(), (mom1+mom2).Phi()*180./3.1416, (mom1+mom2).Eta());
        
        //Fill some histograms about shower shape
        if(fFillSelectClHisto && cluster && GetReader()->GetDataType()!=AliCaloTrackReader::kMC)
        {
          FillSelectedClusterHistograms(cluster, nMaxima, photon1->GetTag());
        }
        
        // Tag both photons as decay
        photon1->SetTagged(kTRUE);
        photon2->SetTagged(kTRUE);
        
        fhPtDecay->Fill(photon1->Pt());
        fhEDecay ->Fill(photon1->E() );
        
        //Create AOD for analysis
        
        mom = mom1+mom2;
        
        //Mass of selected pairs
        fhSelectedMass->Fill(epair,mom.M());
        
        // Fill histograms to undertand pile-up before other cuts applied
        // Remember to relax time cuts in the reader
        if(cluster) FillPileUpHistograms(mom.Pt(),cluster->GetTOF()*1e9,cluster);
        
        AliAODPWG4Particle pi0 = AliAODPWG4Particle(mom);
        
        pi0.SetIdentifiedParticleType(AliCaloPID::kPi0);
        pi0.SetDetector(photon1->GetDetector());
        
        // MC
        pi0.SetLabel(label);
        pi0.SetTag(tag);
        
        //Set the indeces of the original tracks or caloclusters
        pi0.SetCaloLabel(photon1->GetCaloLabel(0), -1);
        pi0.SetTrackLabel(photon2->GetTrackLabel(0), photon2->GetTrackLabel(1));
        //pi0.SetInputFileIndex(input);
        
        AddAODParticle(pi0);
        
      }//pi0
    }//2n photon loop
    
  }//1st photon loop
  
  if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - End fill AODs \n");
  
}


//_________________________________________________
void  AliAnaPi0EbE::MakeShowerShapeIdentification()
{
  //Search for pi0 in fCalorimeter with shower shape analysis
  
  TObjArray * pl        = 0x0;
  AliVCaloCells * cells = 0x0;
  //Select the Calorimeter of the photon
  if      (fCalorimeter == "PHOS" )
  {
    pl    = GetPHOSClusters();
    cells = GetPHOSCells();
  }
  else if (fCalorimeter == "EMCAL")
  {
    pl    = GetEMCALClusters();
    cells = GetEMCALCells();
  }
  
  if(!pl)
  {
    Info("MakeShowerShapeIdentification","TObjArray with %s clusters is NULL!\n",fCalorimeter.Data());
    return;
  }
	
  TLorentzVector mom ;
  for(Int_t icalo = 0; icalo < pl->GetEntriesFast(); icalo++)
  {
    AliVCluster * calo = (AliVCluster*) (pl->At(icalo));
    
    Int_t evtIndex = 0 ;
    if (GetMixedEvent())
    {
      evtIndex=GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ;
    }
    
    if(TMath::Abs(GetVertex(evtIndex)[2]) > GetZvertexCut()) continue ;  //vertex cut
    
    //Get Momentum vector,
    Double_t vertex[]={0,0,0};
    if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    {
      calo->GetMomentum(mom,GetVertex(evtIndex)) ;
    }//Assume that come from vertex in straight line
    else
    {
      calo->GetMomentum(mom,vertex) ;
    }
	  
    //If too small or big pt, skip it
    if(mom.E() < GetMinEnergy() || mom.E() > GetMaxEnergy() ) continue ;
    
    //Check acceptance selection
    if(IsFiducialCutOn())
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
      if(! in ) continue ;
    }
    
    if(GetDebug() > 1)
      printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - FillAOD: Min pt cut and fiducial cut passed: pt %3.2f, phi %2.2f, eta %1.2f\n",mom.Pt(),mom.Phi(),mom.Eta());
    
    //Play with the MC stack if available
    //Check origin of the candidates
    Int_t tag	= 0 ;
    if(IsDataMC())
    {
      tag = GetMCAnalysisUtils()->CheckOrigin(calo->GetLabels(),calo->GetNLabels(),GetReader());
      //GetMCAnalysisUtils()->CheckMultipleOrigin(calo->GetLabels(),calo->GetNLabels(), GetReader(), aodpi0.GetInputFileIndex(), tag);
      if(GetDebug() > 0) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - Origin of candidate %d\n",tag);
    }
    
    //Int_t nMaxima = GetCaloUtils()->GetNumberOfLocalMaxima(calo, cells); // NLM
    
    //Check Distance to Bad channel, set bit.
    Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
    if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
    if(distBad < fMinDist){ //In bad channel (PHOS cristal size 2.2x2.2 cm)
      //FillRejectedClusterHistograms(mom,tag,nMaxima);
      continue ;
    }
    if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - FillAOD: Bad channel cut passed %4.2f\n",distBad);
    
    //If too low number of cells, skip it
    if ( calo->GetNCells() < GetCaloPID()->GetClusterSplittingMinNCells())
    {
      //FillRejectedClusterHistograms(mom,tag,nMaxima);
      continue ;
    }
    
    if(GetDebug() > 1)
      printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - FillAOD: N cells cut passed %d > %d\n",
             calo->GetNCells(), GetCaloPID()->GetClusterSplittingMinNCells());
    
    //.......................................
    // TOF cut, BE CAREFUL WITH THIS CUT
    Double_t tof = calo->GetTOF()*1e9;
    if(tof < fTimeCutMin || tof > fTimeCutMax)
    {
      //FillRejectedClusterHistograms(mom,tag,nMaxima);
      continue ;
    }

    //Check PID
    //PID selection or bit setting
    Int_t    nMaxima  = 0;
    Double_t mass     = 0, angle    = 0;
    Int_t    absId1   =-1, absId2   =-1;
    Float_t  distbad1 =-1, distbad2 =-1;
    Bool_t   fidcut1  = 0, fidcut2  = 0;
    TLorentzVector    l1, l2;

    Int_t idPartType = GetCaloPID()->GetIdentifiedParticleTypeFromClusterSplitting(calo,cells,GetCaloUtils(),
                                                                                   GetVertex(evtIndex),nMaxima,
                                                                                   mass,angle,l1,l2,absId1,absId2,
                                                                                   distbad1,distbad2,fidcut1,fidcut2) ;
    
    
    if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - PDG of identified particle %d\n",idPartType);
    
    
    // Skip events where one of the new clusters (lowest energy) is close to an EMCal border or a bad channel
    if( (fCheckSplitDistToBad) &&
       (!fidcut2 || !fidcut1 || distbad1 < fMinDist || distbad2 < fMinDist))
    {
      if(GetDebug() > 1)
        Info("MakeShowerShapeIdentification", "Dist to bad channel cl %f, cl1 %f, cl2 %f; fid cl1 %d, cl2 %d \n",
               calo->GetDistanceToBadChannel(),distbad1,distbad2, fidcut1,fidcut2);
      
      //FillRejectedClusterHistograms(mom,tag,nMaxima);
      continue ;
    }
    
    //Skip events with too few or too many  NLM
    if(nMaxima < fNLMCutMin || nMaxima > fNLMCutMax)
    {
      //FillRejectedClusterHistograms(mom,tag,nMaxima);
      continue ;
    }
    
    if(GetDebug() > 1)
      printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - NLM %d accepted \n",nMaxima);
    
    //Skip matched clusters with tracks
    if(fRejectTrackMatch && IsTrackMatched(calo, GetReader()->GetInputEvent()))
    {
      FillRejectedClusterHistograms(mom,tag,nMaxima);
      continue ;
    }

    Float_t e1 = l1.Energy();
    Float_t e2 = l2.Energy();
    TLorentzVector l12 = l1+l2;
    Float_t ptSplit = l12.Pt();
    Float_t  eSplit = e1+e2;
    
    Int_t   mcIndex   =-1;
    Int_t   noverlaps = 0;
    Float_t ptprim    = 0;
    if(IsDataMC())
    {
      mcIndex = GetMCIndex(tag);
      
      Bool_t ok      = kFALSE;
      Int_t  mcLabel = calo->GetLabel();
      
      TLorentzVector primary = GetMCAnalysisUtils()->GetMother(mcLabel,GetReader(),ok);
      
      Int_t mesonLabel = -1;
      
      if(mcIndex == kmcPi0 || mcIndex == kmcEta)
      {
        if(mcIndex == kmcPi0)
        {
          TLorentzVector grandmom = GetMCAnalysisUtils()->GetMotherWithPDG(mcLabel,111,GetReader(),ok,mesonLabel);
          if(grandmom.E() > 0 && ok) ptprim =  grandmom.Pt();
        }
        else
        {
          TLorentzVector grandmom = GetMCAnalysisUtils()->GetMotherWithPDG(mcLabel,221,GetReader(),ok,mesonLabel);
          if(grandmom.E() > 0 && ok) ptprim =  grandmom.Pt();
        }
      }
            
      const UInt_t nlabels = calo->GetNLabels();
      Int_t overpdg[nlabels];
      noverlaps = GetMCAnalysisUtils()->GetNOverlaps(calo->GetLabels(), nlabels,tag,mesonLabel,GetReader(),overpdg);
    }
    
    //mass of all clusters
    fhMass       ->Fill(mom.E() ,mass);
    fhMassPt     ->Fill(mom.Pt(),mass);
    fhMassSplitPt->Fill(ptSplit ,mass);
    
    if(IsDataMC())
    {
      fhMCMassPt[mcIndex]     ->Fill(mom.Pt(),mass);
      fhMCMassSplitPt[mcIndex]->Fill(ptSplit ,mass);
      if(mcIndex==kmcPi0)
      {
        fhMCPi0PtRecoPtPrim     ->Fill(mom.Pt(),ptprim);
        fhMCPi0SplitPtRecoPtPrim->Fill(ptSplit ,ptprim);
      }
      else if(mcIndex==kmcEta)
      {
        fhMCEtaPtRecoPtPrim     ->Fill(mom.Pt(),ptprim);
        fhMCEtaSplitPtRecoPtPrim->Fill(ptSplit ,ptprim);
      }

      if(noverlaps==0)
      {
        if(mcIndex==kmcPi0)
        {
          fhMCPi0PtRecoPtPrimNoOverlap     ->Fill(mom.Pt(),ptprim);
          fhMCPi0SplitPtRecoPtPrimNoOverlap->Fill(ptSplit ,ptprim);
        }
        else if(mcIndex==kmcEta)
        {
          fhMCEtaPtRecoPtPrimNoOverlap     ->Fill(mom.Pt(),ptprim);
          fhMCEtaSplitPtRecoPtPrimNoOverlap->Fill(ptSplit ,ptprim);
        }
        
        fhMassNoOverlap       ->Fill(mom.E() ,mass);
        fhMassPtNoOverlap     ->Fill(mom.Pt(),mass);
        fhMassSplitPtNoOverlap->Fill(ptSplit ,mass);
        
        fhMCMassPtNoOverlap[mcIndex]     ->Fill(mom.Pt(),mass);
        fhMCMassSplitPtNoOverlap[mcIndex]->Fill(ptSplit ,mass);
      }
    }
    
    // Asymmetry of all clusters
    Float_t asy =-10;
    
    if(e1+e2 > 0) asy = (e1-e2) / (e1+e2);
    fhAsymmetry->Fill(mom.E(),asy);
    
    if(IsDataMC())
    {
      fhMCEAsymmetry[mcIndex]->Fill(mom.E(),asy);
    }
    
    // If cluster does not pass pid, not pi0/eta, skip it.
    if     (GetOutputAODName().Contains("Pi0") && idPartType != AliCaloPID::kPi0)
    {
      if(GetDebug() > 1) Info("MakeShowerShapeIdentification","Cluster is not Pi0\n");
      FillRejectedClusterHistograms(mom,tag,nMaxima);
      continue ;
    }
    
    else if(GetOutputAODName().Contains("Eta") && idPartType != AliCaloPID::kEta)
    {
      if(GetDebug() > 1) Info("MakeShowerShapeIdentification","Cluster is not Eta\n");
      FillRejectedClusterHistograms(mom,tag,nMaxima);
      continue ;
    }
    
    if(GetDebug() > 1)
      Info("MakeShowerShapeIdentification","Pi0/Eta selection cuts passed: pT %3.2f, pdg %d\n",
             mom.Pt(), idPartType);
    
    //Mass and asymmetry of selected pairs
    fhSelectedAsymmetry  ->Fill(mom.E() ,asy );
    fhSelectedMass       ->Fill(mom.E() ,mass);
    fhSelectedMassPt     ->Fill(mom.Pt(),mass);
    fhSelectedMassSplitPt->Fill(ptSplit ,mass);
    
    if(IsDataMC())
    {
      if(mcIndex==kmcPi0)
      {
        fhMCPi0SelectedPtRecoPtPrim     ->Fill(mom.Pt(),ptprim);
        fhMCPi0SelectedSplitPtRecoPtPrim->Fill(ptSplit ,ptprim);
      }
      else if(mcIndex==kmcEta)
      {
        fhMCEtaSelectedPtRecoPtPrim     ->Fill(mom.Pt(),ptprim);
        fhMCEtaSelectedSplitPtRecoPtPrim->Fill(ptSplit ,ptprim);
      }
      
      if(noverlaps==0)
      {
        fhSelectedMassNoOverlap       ->Fill(mom.E() ,mass);
        fhSelectedMassPtNoOverlap     ->Fill(mom.Pt(),mass);
        fhSelectedMassSplitPtNoOverlap->Fill(ptSplit ,mass);
        
        if(mcIndex==kmcPi0)
        {
          fhMCPi0SelectedPtRecoPtPrimNoOverlap     ->Fill(mom.Pt(),ptprim);
          fhMCPi0SelectedSplitPtRecoPtPrimNoOverlap->Fill(ptSplit ,ptprim);
        }
        else if(mcIndex==kmcEta)
        {
          fhMCEtaSelectedPtRecoPtPrimNoOverlap     ->Fill(mom.Pt(),ptprim);
          fhMCEtaSelectedSplitPtRecoPtPrimNoOverlap->Fill(ptSplit ,ptprim);
        }
      }
    }
    
    fhSplitE        ->Fill( eSplit);
    fhSplitPt       ->Fill(ptSplit);
    Float_t phi = mom.Phi();
    if(phi<0) phi+=TMath::TwoPi();
    fhSplitPtPhi    ->Fill(ptSplit,phi);
    fhSplitPtEta    ->Fill(ptSplit,mom.Eta());
    fhNLocMaxSplitPt->Fill(ptSplit ,nMaxima);
    fhNLocMaxPt     ->Fill(mom.Pt(),nMaxima);
    
    //Check split-clusters with good time window difference
    Double_t tof1  = cells->GetCellTime(absId1);
    GetCaloUtils()->RecalibrateCellTime(tof1, fCalorimeter, absId1,GetReader()->GetInputEvent()->GetBunchCrossNumber());
    tof1*=1.e9;
    
    Double_t tof2  = cells->GetCellTime(absId2);
    GetCaloUtils()->RecalibrateCellTime(tof2, fCalorimeter, absId2,GetReader()->GetInputEvent()->GetBunchCrossNumber());
    tof2*=1.e9;
    
    Double_t t12diff = tof1-tof2;
    fhEPairDiffTime->Fill(e1+e2,    t12diff);
    
    if(IsDataMC())
    {
      fhMCSplitE        [mcIndex]->Fill( eSplit);
      fhMCSplitPt       [mcIndex]->Fill(ptSplit);
      fhMCSplitPtPhi    [mcIndex]->Fill(ptSplit,phi);
      fhMCSplitPtEta    [mcIndex]->Fill(ptSplit,mom.Eta());
      fhMCNLocMaxSplitPt[mcIndex]->Fill(ptSplit ,nMaxima);
      fhMCNLocMaxPt     [mcIndex]->Fill(mom.Pt(),nMaxima);
      
      fhMCSelectedMassPt     [mcIndex]->Fill(mom.Pt(),mass);
      fhMCSelectedMassSplitPt[mcIndex]->Fill(ptSplit,mass);
      if(noverlaps==0)
      {
        fhMCSelectedMassPtNoOverlap     [mcIndex]->Fill(mom.Pt(),mass);
        fhMCSelectedMassSplitPtNoOverlap[mcIndex]->Fill(ptSplit,mass);
      }
    }
    
    //-----------------------
    //Create AOD for analysis
    
    if(nMaxima == 1 && fNLMECutMin[0] > mom.E()) continue;
    if(nMaxima == 2 && fNLMECutMin[1] > mom.E()) continue;
    if(nMaxima >  2 && fNLMECutMin[2] > mom.E()) continue;
    
    AliAODPWG4Particle aodpi0 = AliAODPWG4Particle(mom);
    aodpi0.SetLabel(calo->GetLabel());
    
    //Set the indeces of the original caloclusters
    aodpi0.SetCaloLabel(calo->GetID(),-1);
    aodpi0.SetDetector(fCalorimeter);
    
    if     (distBad > fMinDist3) aodpi0.SetDistToBad(2) ;
    else if(distBad > fMinDist2) aodpi0.SetDistToBad(1) ;
    else                         aodpi0.SetDistToBad(0) ;
    
    // Check if cluster is pi0 via cluster splitting
    aodpi0.SetIdentifiedParticleType(idPartType);
    
    // Add number of local maxima to AOD, method name in AOD to be FIXED
    aodpi0.SetFiducialArea(nMaxima);
    
    aodpi0.SetTag(tag);
    
    //Fill some histograms about shower shape
    if(fFillSelectClHisto && GetReader()->GetDataType()!=AliCaloTrackReader::kMC)
    {
      FillSelectedClusterHistograms(calo, nMaxima, tag, asy);
    }
    
    // Fill histograms to undertand pile-up before other cuts applied
    // Remember to relax time cuts in the reader
    Double_t tofcluster   = calo->GetTOF()*1e9;
    Double_t tofclusterUS = TMath::Abs(tofcluster);
    
    FillPileUpHistograms(aodpi0.Pt(),tofcluster,calo);
    
    Int_t id = GetReader()->GetTriggerClusterId();
    if(fFillEMCALBCHistograms && fCalorimeter=="EMCAL" && id >=0 )
    {
      Float_t phicluster = aodpi0.Phi();
      if(phicluster < 0) phicluster+=TMath::TwoPi();
      
      if(calo->E() > 2)
      {
        if      (tofclusterUS < 25) fhEtaPhiEMCALBC0->Fill(aodpi0.Eta(), phicluster);
        else if (tofclusterUS < 75) fhEtaPhiEMCALBC1->Fill(aodpi0.Eta(), phicluster);
        else                        fhEtaPhiEMCALBCN->Fill(aodpi0.Eta(), phicluster);
      }
      
      Int_t bc = GetReader()->GetTriggerClusterBC();
      if(TMath::Abs(bc) < 6  && !GetReader()->IsBadCellTriggerEvent() && !GetReader()->IsExoticEvent() )
      {
        if(GetReader()->IsTriggerMatched())
        {
          if(calo->E() > 2) fhEtaPhiTriggerEMCALBC[bc+5]->Fill(aodpi0.Eta(), phicluster);
          fhTimeTriggerEMCALBC[bc+5]->Fill(calo->E(), tofcluster);
          if(GetReader()->IsPileUpFromSPD()) fhTimeTriggerEMCALBCPileUpSPD[bc+5]->Fill(calo->E(), tofcluster);
        }
        else
        {
          if(calo->E() > 2) fhEtaPhiTriggerEMCALBCUM[bc+5]->Fill(aodpi0.Eta(), phicluster);
          fhTimeTriggerEMCALBCUM[bc+5]->Fill(calo->E(), tofcluster);
          
          if(bc==0)
          {
            if(GetReader()->IsTriggerMatchedOpenCuts(0)) fhTimeTriggerEMCALBC0UMReMatchOpenTime   ->Fill(calo->E(), tofcluster);
            if(GetReader()->IsTriggerMatchedOpenCuts(1)) fhTimeTriggerEMCALBC0UMReMatchCheckNeigh ->Fill(calo->E(), tofcluster);
            if(GetReader()->IsTriggerMatchedOpenCuts(2)) fhTimeTriggerEMCALBC0UMReMatchBoth       ->Fill(calo->E(), tofcluster);
          }
         }
      }
      else if(TMath::Abs(bc) >= 6)
        Info("MakeShowerShapeIdentification","Trigger BC not expected = %d\n",bc);
    }
    
    //Add AOD with pi0 object to aod branch
    AddAODParticle(aodpi0);
    
  }//loop
  
  if(GetDebug() > 1) Info("MakeShowerShapeIdentification","End fill AODs \n");
  
}
//______________________________________________
void  AliAnaPi0EbE::MakeAnalysisFillHistograms()
{
  //Do analysis and fill histograms
  
  if(!GetOutputAODBranch())
  {
    AliFatal(Form("No output pi0 in AOD branch with name < %s >,STOP \n",GetOutputAODName().Data()));
  }
  //Loop on stored AOD pi0
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) Info("MakeAnalysisFillHistograms","aod branch entries %d\n", naod);
  
  Float_t cen = GetEventCentrality();
  Float_t ep  = GetEventPlaneAngle();
  
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    
    AliAODPWG4Particle* pi0 =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = pi0->GetIdentifiedParticleType();
	  
    if(IsCaloPIDOn() && pdg != AliCaloPID::kPi0) continue;
    
    //Fill pi0 histograms
    Float_t ener  = pi0->E();
    Float_t pt    = pi0->Pt();
    Float_t phi   = pi0->Phi();
    if(phi < 0) phi+=TMath::TwoPi();
    Float_t eta = pi0->Eta();
    
    fhPt     ->Fill(pt  );
    fhE      ->Fill(ener);
    
    fhEEta   ->Fill(ener,eta);
    fhEPhi   ->Fill(ener,phi);
    fhPtEta  ->Fill(pt  ,eta);
    fhPtPhi  ->Fill(pt  ,phi);
    fhEtaPhi ->Fill(eta ,phi);
    
    fhPtCentrality ->Fill(pt,cen) ;
    fhPtEventPlane ->Fill(pt,ep ) ;
    
    if(IsDataMC())
    {
      Int_t tag     = pi0->GetTag();
      Int_t mcIndex = GetMCIndex(tag);
      
      fhMCE  [mcIndex] ->Fill(ener);
      fhMCPt [mcIndex] ->Fill(pt);
      fhMCPhi[mcIndex] ->Fill(pt,phi);
      fhMCEta[mcIndex] ->Fill(pt,eta);
      
      fhMCPtCentrality[mcIndex]->Fill(pt,cen);
      
      if((mcIndex==kmcPhoton || mcIndex==kmcPi0 || mcIndex==kmcEta) && fAnaType==kSSCalo)
      {
        Float_t efracMC   = 0;
        Int_t   label     = pi0->GetLabel();
        Int_t   momlabel  = -1;
        Bool_t  ok        = kFALSE;
        
        TLorentzVector mom   = GetMCAnalysisUtils()->GetMother(label,GetReader(),ok);
        if(!ok) continue;
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0))
        {
          TLorentzVector grandmom = GetMCAnalysisUtils()->GetMotherWithPDG(label,111,GetReader(),ok,momlabel);
          if(grandmom.E() > 0 && ok)
          {
            efracMC =  grandmom.E()/ener;
            fhMCPi0PtGenRecoFraction ->Fill(pt,efracMC);
          }
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))
        {
          fhMCPi0DecayPt->Fill(pt);
          TLorentzVector grandmom = GetMCAnalysisUtils()->GetMotherWithPDG(label,111,GetReader(),ok,momlabel);
          if(grandmom.E() > 0 && ok)
          {
            efracMC =  mom.E()/grandmom.E();
            fhMCPi0DecayPtFraction ->Fill(pt,efracMC);
          }
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta))
        {
          TLorentzVector grandmom = GetMCAnalysisUtils()->GetMotherWithPDG(label,221,GetReader(),ok,momlabel);
          if(grandmom.E() > 0 && ok)
          {
            efracMC =  grandmom.E()/ener;
            fhMCEtaPtGenRecoFraction ->Fill(pt,efracMC);
          }
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay))
        {
          fhMCEtaDecayPt->Fill(pt);
          TLorentzVector grandmom = GetMCAnalysisUtils()->GetMotherWithPDG(label,221,GetReader(),ok,momlabel);
          if(grandmom.E() > 0 && ok)
          {
            efracMC =  mom.E()/grandmom.E();
            fhMCEtaDecayPtFraction ->Fill(pt,efracMC);
          }
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))
        {
          fhMCOtherDecayPt->Fill(pt);
        }
        
      }
      
    }//Histograms with MC
    
  }// aod loop
  
}

//__________________________________________________________________
void AliAnaPi0EbE::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print("");
  printf("Analysis Type = %d \n",  fAnaType) ;
  if(fAnaType == kSSCalo)
  {
    printf("Calorimeter            =     %s\n", fCalorimeter.Data()) ;
    printf("Min Distance to Bad Channel   = %2.1f\n",fMinDist);
    printf("Min Distance to Bad Channel 2 = %2.1f\n",fMinDist2);
    printf("Min Distance to Bad Channel 3 = %2.1f\n",fMinDist3);
  } 
  printf("    \n") ;
  
} 


