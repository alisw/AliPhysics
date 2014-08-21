/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes hereby granted      *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//_________________________________________________________________________
// Class for analysis of particle isolation
// Input is selected particles put in AOD branch (AliAODPWG4ParticleCorrelation)
//
//  Class created from old AliPHOSGammaJet
//  (see AliRoot versions previous Release 4-09)
//
// -- Author: Gustavo Conesa (LNF-INFN)

//-Yaxian Mao (add the possibility for different IC method with different pt range, 01/10/2010)
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include <TClonesArray.h>
#include <TList.h>
#include <TObjString.h>
#include <TH2F.h>
#include <TClass.h>
#include <TH2F.h>
#include "TParticle.h"
#include "TDatabasePDG.h"



// --- Analysis system ---
#include "AliAnaParticleIsolation.h"
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliIsolationCut.h"
#include "AliFiducialCut.h"
#include "AliMCAnalysisUtils.h"
#include "AliNeutralMesonSelection.h"
#include "AliAODMCParticle.h"
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliMCAnalysisUtils.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
// --- Detectors ---
#include "AliEMCALGeometry.h"
#include "AliPHOSGeoUtils.h"

ClassImp(AliAnaParticleIsolation)

//______________________________________________________________________________
AliAnaParticleIsolation::AliAnaParticleIsolation() :
AliAnaCaloTrackCorrBaseClass(),
fCalorimeter(""),                 fIsoDetector(""),
fReMakeIC(0),                     fMakeSeveralIC(0),
fFillPileUpHistograms(0),
fFillTMHisto(0),                  fFillSSHisto(1),
fFillUEBandSubtractHistograms(1), fFillCellHistograms(0),
fFillHighMultHistograms(0),       fFillTaggedDecayHistograms(0),
fNDecayBits(0),                   fDecayBits(),
fFillNLMHistograms(0),
fLeadingOnly(0),                  fCheckLeadingWithNeutralClusters(0),
// Several IC
fNCones(0),                       fNPtThresFrac(0),
fConeSizes(),                     fPtThresholds(),
fPtFractions(),                   fSumPtThresholds(),
// Histograms
fhEIso(0),                        fhPtIso(0),
fhPtCentralityIso(0),             fhPtEventPlaneIso(0),
fhPtNLocMaxIso(0),
fhPhiIso(0),                      fhEtaIso(0),                              fhEtaPhiIso(0),
fhEtaPhiNoIso(0),
fhENoIso(0),                      fhPtNoIso(0),                             fhPtNLocMaxNoIso(0),
fhPtInCone(0),
fhPtClusterInCone(0),             fhPtCellInCone(0),                        fhPtTrackInCone(0),
fhPtTrackInConeOtherBC(0),        fhPtTrackInConeOtherBCPileUpSPD(0),
fhPtTrackInConeBC0(0),            fhPtTrackInConeVtxBC0(0),
fhPtTrackInConeBC0PileUpSPD(0),
fhPtInConePileUp(),               fhPtInConeCent(0),
fhPerpConeSumPt(0),               fhPtInPerpCone(0),
fhEtaPhiInConeCluster(0),         fhEtaPhiCluster(0),
fhEtaPhiInConeTrack(0),           fhEtaPhiTrack(0),
fhEtaBandCluster(0),              fhPhiBandCluster(0),
fhEtaBandTrack(0),                fhPhiBandTrack(0),
fhEtaBandCell(0),                 fhPhiBandCell(0),
fhConePtLead(0),                  fhConePtLeadCluster(0),                   fhConePtLeadTrack(0),
fhConeSumPt(0),                   fhConeSumPtCellTrack(0),
fhConeSumPtCell(0),               fhConeSumPtCluster(0),                    fhConeSumPtTrack(0),
fhConeSumPtEtaBandUECluster(0),             fhConeSumPtPhiBandUECluster(0),
fhConeSumPtEtaBandUETrack(0),               fhConeSumPtPhiBandUETrack(0),
fhConeSumPtEtaBandUECell(0),                fhConeSumPtPhiBandUECell(0),
fhConeSumPtTrigEtaPhi(0),
fhConeSumPtCellTrackTrigEtaPhi(0),
fhConeSumPtEtaBandUEClusterTrigEtaPhi(0),   fhConeSumPtPhiBandUEClusterTrigEtaPhi(0),
fhConeSumPtEtaBandUETrackTrigEtaPhi(0),     fhConeSumPtPhiBandUETrackTrigEtaPhi(0),
fhConeSumPtEtaBandUECellTrigEtaPhi(0),      fhConeSumPtPhiBandUECellTrigEtaPhi(0),
fhConeSumPtEtaUESub(0),                     fhConeSumPtPhiUESub(0),
fhConeSumPtEtaUESubTrigEtaPhi(0),           fhConeSumPtPhiUESubTrigEtaPhi(0),
fhConeSumPtEtaUESubTrackCell(0),            fhConeSumPtPhiUESubTrackCell(0),
fhConeSumPtEtaUESubTrackCellTrigEtaPhi(0),  fhConeSumPtPhiUESubTrackCellTrigEtaPhi(0),
fhConeSumPtEtaUESubCluster(0),              fhConeSumPtPhiUESubCluster(0),
fhConeSumPtEtaUESubClusterTrigEtaPhi(0),    fhConeSumPtPhiUESubClusterTrigEtaPhi(0),
fhConeSumPtEtaUESubCell(0),                 fhConeSumPtPhiUESubCell(0),
fhConeSumPtEtaUESubCellTrigEtaPhi(0),       fhConeSumPtPhiUESubCellTrigEtaPhi(0),
fhConeSumPtEtaUESubTrack(0),                fhConeSumPtPhiUESubTrack(0),
fhConeSumPtEtaUESubTrackTrigEtaPhi(0),      fhConeSumPtPhiUESubTrackTrigEtaPhi(0),
fhFractionTrackOutConeEta(0),               fhFractionTrackOutConeEtaTrigEtaPhi(0),
fhFractionClusterOutConeEta(0),             fhFractionClusterOutConeEtaTrigEtaPhi(0),
fhFractionClusterOutConePhi(0),             fhFractionClusterOutConePhiTrigEtaPhi(0),
fhFractionCellOutConeEta(0),                fhFractionCellOutConeEtaTrigEtaPhi(0),
fhFractionCellOutConePhi(0),                fhFractionCellOutConePhiTrigEtaPhi(0),
fhConeSumPtClustervsTrack(0),
fhConeSumPtEtaUESubClustervsTrack(0),       fhConeSumPtPhiUESubClustervsTrack(0),
fhConeSumPtCellvsTrack(0),
fhConeSumPtEtaUESubCellvsTrack(0),          fhConeSumPtPhiUESubCellvsTrack(0),
fhEtaBandClustervsTrack(0),                 fhPhiBandClustervsTrack(0),
fhEtaBandNormClustervsTrack(0),             fhPhiBandNormClustervsTrack(0),
fhEtaBandCellvsTrack(0),                    fhPhiBandCellvsTrack(0),
fhEtaBandNormCellvsTrack(0),                fhPhiBandNormCellvsTrack(0),
fhConeSumPtSubvsConeSumPtTotPhiTrack(0),    fhConeSumPtSubNormvsConeSumPtTotPhiTrack(0),
fhConeSumPtSubvsConeSumPtTotEtaTrack(0),    fhConeSumPtSubNormvsConeSumPtTotEtaTrack(0),
fhConeSumPtSubvsConeSumPtTotPhiCluster(0),  fhConeSumPtSubNormvsConeSumPtTotPhiCluster(0),
fhConeSumPtSubvsConeSumPtTotEtaCluster(0),  fhConeSumPtSubNormvsConeSumPtTotEtaCluster(0),
fhConeSumPtSubvsConeSumPtTotPhiCell(0),     fhConeSumPtSubNormvsConeSumPtTotPhiCell(0),
fhConeSumPtSubvsConeSumPtTotEtaCell(0),     fhConeSumPtSubNormvsConeSumPtTotEtaCell(0),
fhConeSumPtVSUETracksEtaBand(0),            fhConeSumPtVSUETracksPhiBand(0),
fhConeSumPtVSUEClusterEtaBand(0),           fhConeSumPtVSUEClusterPhiBand(0),

// Number of local maxima in cluster
fhNLocMax(),
fhELambda0LocMax1(),              fhELambda1LocMax1(),
fhELambda0LocMax2(),              fhELambda1LocMax2(),
fhELambda0LocMaxN(),              fhELambda1LocMaxN(),
// PileUp
fhEIsoPileUp(),                   fhPtIsoPileUp(),
fhENoIsoPileUp(),                 fhPtNoIsoPileUp(),
fhTimeENoCut(0),                  fhTimeESPD(0),                  fhTimeESPDMulti(0),
fhTimeNPileUpVertSPD(0),          fhTimeNPileUpVertTrack(0),
fhTimeNPileUpVertContributors(0),
fhTimePileUpMainVertexZDistance(0), fhTimePileUpMainVertexZDiamond(0)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();
  
  for(Int_t i = 0; i < 5 ; i++)
  {
    fConeSizes[i]      = 0 ;
    
    for(Int_t imc = 0; imc < 9; imc++)
      fhSumPtLeadingPtMC[imc][i] = 0 ;
    
    for(Int_t j = 0; j < 5 ; j++)
    {
      fhPtThresIsolated             [i][j] = 0 ;
      fhPtFracIsolated              [i][j] = 0 ;
      fhSumPtIsolated               [i][j] = 0 ;
      
      fhEtaPhiPtThresIso            [i][j] = 0 ;
      fhEtaPhiPtThresDecayIso       [i][j] = 0 ;
      fhPtPtThresDecayIso           [i][j] = 0 ;
      
      fhEtaPhiPtFracIso             [i][j] = 0 ;
      fhEtaPhiPtFracDecayIso        [i][j] = 0 ;
      fhPtPtFracDecayIso            [i][j] = 0 ;
      fhPtPtSumDecayIso             [i][j] = 0 ;
      fhPtSumDensityIso             [i][j] = 0 ;
      fhPtSumDensityDecayIso        [i][j] = 0 ;
      fhEtaPhiSumDensityIso         [i][j] = 0 ;
      fhEtaPhiSumDensityDecayIso    [i][j] = 0 ;
      fhPtFracPtSumIso              [i][j] = 0 ;
      fhPtFracPtSumDecayIso         [i][j] = 0 ;
      fhEtaPhiFracPtSumIso          [i][j] = 0 ;
      fhEtaPhiFracPtSumDecayIso     [i][j] = 0 ;
      
      for(Int_t imc = 0; imc < 9; imc++)
      {
        fhPtThresIsolatedMC[imc][i][j] = 0 ;
        fhPtFracIsolatedMC [imc][i][j] = 0 ;
        fhSumPtIsolatedMC  [imc][i][j] = 0 ;
        
      }
    }
  }
  
  for(Int_t ibit =0; ibit< 4; ibit++)
  {
    fhPtDecayIso      [ibit] = 0;
    fhPtDecayNoIso    [ibit] = 0;
    fhEtaPhiDecayIso  [ibit] = 0;
    fhEtaPhiDecayNoIso[ibit] = 0;
    for(Int_t imc = 0; imc < 9; imc++)
    {
      fhPtDecayIsoMC  [ibit][imc]    = 0;
      fhPtDecayNoIsoMC[ibit][imc]    = 0;
    }
  }
  
  for(Int_t i = 0; i < 5 ; i++)
  {
    fPtFractions    [i] = 0 ;
    fPtThresholds   [i] = 0 ;
    fSumPtThresholds[i] = 0 ;
    
    fhSumPtLeadingPt    [i] = 0 ;
    fhPtLeadingPt       [i] = 0 ;
    fhPerpSumPtLeadingPt[i] = 0 ;
    fhPerpPtLeadingPt   [i] = 0 ;
  }
  
  for(Int_t imc = 0; imc < 9; imc++)
  {
    fhPtNoIsoMC  [imc]    = 0;
    fhPtIsoMC    [imc]    = 0;
    fhPhiIsoMC   [imc]    = 0;
    fhEtaIsoMC   [imc]    = 0;
    fhPtLambda0MC[imc][0] = 0;
    fhPtLambda0MC[imc][1] = 0;
  }
  
  for(Int_t i = 0; i < 2 ; i++)
  {
    fhTrackMatchedDEta[i] = 0 ;             fhTrackMatchedDPhi[i] = 0 ;   fhTrackMatchedDEtaDPhi  [i] = 0 ;
    fhdEdx            [i] = 0 ;             fhEOverP          [i] = 0 ;   fhTrackMatchedMCParticle[i] = 0 ;
    fhELambda0        [i] = 0 ;             fhELambda1        [i] = 0 ;   fhPtLambda0       [i] = 0 ;
    fhELambda0TRD     [i] = 0 ;             fhELambda1TRD     [i] = 0 ;   fhPtLambda0TRD    [i] = 0 ;
    
    // Number of local maxima in cluster
    fhNLocMax        [i] = 0 ;
    fhELambda0LocMax1[i] = 0 ;              fhELambda1LocMax1[i] = 0 ;
    fhELambda0LocMax2[i] = 0 ;              fhELambda1LocMax2[i] = 0 ;
    fhELambda0LocMaxN[i] = 0 ;              fhELambda1LocMaxN[i] = 0 ;
  }
  
  // Acceptance
  for(Int_t i = 0; i < 6; i++)
  {
    fhPtPrimMCiso[i] = 0;
    fhEPrimMC    [i] = 0;
    fhEtaPrimMC  [i] = 0;
    fhPhiPrimMC  [i] = 0;
  }
  
  // Pile-Up
  
  for(Int_t i = 0 ; i < 7 ; i++)
  {
    fhPtInConePileUp[i] = 0 ;
    fhEIsoPileUp    [i] = 0 ;
    fhPtIsoPileUp   [i] = 0 ;
    fhENoIsoPileUp  [i] = 0 ;
    fhPtNoIsoPileUp [i] = 0 ;
  }
  
}

//_______________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateCaloUEBand(AliAODPWG4ParticleCorrelation * pCandidate,
                                                  Float_t & etaBandPtSum, Float_t & phiBandPtSum)
{
  // Get the clusters pT or sum of pT in phi/eta bands or at 45 degrees from trigger
  
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyCharged ) return ;
  
  Float_t conesize   = GetIsolationCut()->GetConeSize();
  TLorentzVector mom ;
  
  //Select the Calorimeter
  TObjArray * pl = 0x0;
  if      (fCalorimeter == "PHOS" )
    pl    = GetPHOSClusters();
  else if (fCalorimeter == "EMCAL")
    pl    = GetEMCALClusters();
  
  if(!pl) return ;
  
  //Get vertex for cluster momentum calculation
  Double_t vertex[] = {0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    GetReader()->GetVertex(vertex);
  
  Float_t ptTrig    = pCandidate->Pt() ;
  Float_t phiTrig   = pCandidate->Phi();
  Float_t etaTrig   = pCandidate->Eta();
  
  for(Int_t icluster=0; icluster < pl->GetEntriesFast(); icluster++)
  {
    AliVCluster* cluster = (AliVCluster *) pl->At(icluster);
    
    if(!cluster)
    {
      printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Cluster not available?");
      continue;
    }
    
    //Do not count the candidate (photon or pi0) or the daughters of the candidate
    if(cluster->GetID() == pCandidate->GetCaloLabel(0) ||
       cluster->GetID() == pCandidate->GetCaloLabel(1)   ) continue ;
    
    //Remove matched clusters to tracks if Neutral and Track info is used
    if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged &&
       IsTrackMatched(cluster,GetReader()->GetInputEvent())) continue ;
    
    cluster->GetMomentum(mom,vertex) ;//Assume that come from vertex in straight line
    
    //exclude particles in cone
    Float_t rad = GetIsolationCut()->Radius(etaTrig, phiTrig, mom.Eta(), mom.Phi());
    
    // histo of eta and phi for all clusters
    fhEtaPhiCluster->Fill(mom.Eta(), mom.Phi());
    if(rad < conesize) {
    	// histos for all clusters in cone
      fhEtaPhiInConeCluster->Fill(mom.Eta(), mom.Phi());
      continue ;
    }
    //fill histogram for UE in phi band in EMCal acceptance
    if(mom.Eta() > (etaTrig-conesize) && mom.Eta()  < (etaTrig+conesize))
    {
      phiBandPtSum+=mom.Pt();
      fhPhiBandCluster->Fill(mom.Eta(),mom.Phi());
      
    }
    
    //fill histogram for UE in eta band in EMCal acceptance
    if(mom.Phi() > (phiTrig-conesize) && mom.Phi() < (phiTrig+conesize))
    {
      etaBandPtSum+=mom.Pt();
      fhEtaBandCluster->Fill(mom.Eta(),mom.Phi());
    }
  }
  
  fhConeSumPtEtaBandUECluster          ->Fill(ptTrig  ,       etaBandPtSum);
  fhConeSumPtPhiBandUECluster          ->Fill(ptTrig  ,       phiBandPtSum);
  fhConeSumPtEtaBandUEClusterTrigEtaPhi->Fill(etaTrig,phiTrig,etaBandPtSum);
  fhConeSumPtPhiBandUEClusterTrigEtaPhi->Fill(etaTrig,phiTrig,phiBandPtSum);
  
}

//________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateCaloCellUEBand(AliAODPWG4ParticleCorrelation * pCandidate,
                                                      Float_t & etaBandPtSumCells, Float_t & phiBandPtSumCells)
{
  // Get the cells amplitude or sum of amplitude in phi/eta bands or at 45 degrees from trigger
  
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyCharged ) return ;
  
  Float_t conesize = GetIsolationCut()->GetConeSize();
  
  Float_t phiTrig = pCandidate->Phi();
  if(phiTrig<0) phiTrig += TMath::TwoPi();
  Float_t etaTrig = pCandidate->Eta();
  
  if(pCandidate->GetDetector()=="EMCAL")
  {
    AliEMCALGeometry* eGeom = AliEMCALGeometry::GetInstance();
    Int_t absId = -999;
    
    if (eGeom->GetAbsCellIdFromEtaPhi(etaTrig,phiTrig,absId))
    {
      if(!eGeom->CheckAbsCellId(absId)) return ;
      
      // Get absolute (col,row) of trigger particle
      Int_t nSupMod = eGeom->GetSuperModuleNumber(absId);
      Int_t nModule = -1;
      Int_t imEta=-1, imPhi=-1;
      Int_t ieta =-1, iphi =-1;
      
      if (eGeom->GetCellIndex(absId,nSupMod,nModule,imPhi,imEta))
      {
        eGeom->GetCellPhiEtaIndexInSModule(nSupMod,nModule,imPhi,imEta,iphi,ieta);
        
        Int_t colTrig = ieta;
        if (nSupMod % 2) colTrig = AliEMCALGeoParams::fgkEMCALCols + ieta ;
        Int_t rowTrig = iphi + AliEMCALGeoParams::fgkEMCALRows*int(nSupMod/2);
        
        Int_t sqrSize = int(conesize/0.0143);
        
        AliVCaloCells * cells = GetEMCALCells();
        
        Int_t nTotalRows = AliEMCALGeoParams::fgkEMCALRows*16/3 ; // 24*(16/3) 5 full-size Sectors (2 SM) + 1 one-third Sector (2 SM)
        Int_t nTotalCols = 2*AliEMCALGeoParams::fgkEMCALCols;
        //  printf("nTotalRows %i, nTotalCols %i\n",nTotalRows,nTotalCols);
        // Loop on cells in eta band
        
	  	  Int_t irowmin = rowTrig-sqrSize;
        if(irowmin<0) irowmin=0;
        Int_t irowmax = rowTrig+sqrSize;
        if(irowmax>AliEMCALGeoParams::fgkEMCALRows) irowmax=AliEMCALGeoParams::fgkEMCALRows;
        
        
        for(Int_t irow = irowmin; irow <irowmax; irow++)
        {
          for(Int_t icol = 0; icol < nTotalCols; icol++)
          {
            Int_t inSector = int(irow/AliEMCALGeoParams::fgkEMCALRows);
            if(inSector==5) continue;
            Int_t inSupMod = -1;
            Int_t icolLoc  = -1;
            if(icol < AliEMCALGeoParams::fgkEMCALCols)
            {
              inSupMod = 2*inSector + 1;
              icolLoc  = icol;
            }
            else if(icol > AliEMCALGeoParams::fgkEMCALCols - 1)
            {
              inSupMod = 2*inSector;
              icolLoc  = icol-AliEMCALGeoParams::fgkEMCALCols;
            }
            
            Int_t irowLoc  = irow - AliEMCALGeoParams::fgkEMCALRows*inSector ;
            
            // Exclude cells in cone
            if(TMath::Abs(icol-colTrig) < sqrSize || TMath::Abs(irow-rowTrig) < sqrSize){
              continue ;
            }
            Int_t iabsId = eGeom->GetAbsCellIdFromCellIndexes(inSupMod,irowLoc,icolLoc);
            if(!eGeom->CheckAbsCellId(iabsId)) continue;
            etaBandPtSumCells += cells->GetCellAmplitude(iabsId);
            fhEtaBandCell->Fill(colTrig,rowTrig);
            
            //          printf("ETA inSupMod %i,irowLoc %i,icolLoc %i, iabsId %i, etaBandPtSumCells %f\n",nSupMod,irowLoc,icolLoc,iabsId,etaBandPtSumCells);
          }
        }
        Int_t icolmin = colTrig-sqrSize;
        if(icolmin<0) icolmin=0;
        Int_t icolmax = colTrig+sqrSize;
        if(icolmax>AliEMCALGeoParams::fgkEMCALCols) icolmax=AliEMCALGeoParams::fgkEMCALCols;
	      
        // Loop on cells in phi band
        for(Int_t icol = icolmin; icol < icolmax; icol++)
        {
          for(Int_t irow = 0; irow < nTotalRows; irow++)
          {
            Int_t inSector = int(irow/AliEMCALGeoParams::fgkEMCALRows);
            if(inSector==5) continue;
            Int_t inSupMod = -1;
            Int_t icolLoc  = -1;
            //    printf("icol %i, irow %i, inSector %i\n",icol,irow ,inSector);
            if(icol < AliEMCALGeoParams::fgkEMCALCols)
            {
              //	printf("icol < AliEMCALGeoParams::fgkEMCALCols %i\n",AliEMCALGeoParams::fgkEMCALCols );
              inSupMod = 2*inSector + 1;
              icolLoc  = icol;
            }
            else if(icol > AliEMCALGeoParams::fgkEMCALCols - 1)
            {
              //      printf("icol > AliEMCALGeoParams::fgkEMCALCols -1 %i\n",AliEMCALGeoParams::fgkEMCALCols -1 );
              inSupMod = 2*inSector;
              icolLoc  = icol-AliEMCALGeoParams::fgkEMCALCols;
            }
            
            Int_t irowLoc  = irow - AliEMCALGeoParams::fgkEMCALRows*inSector ;   // Stesso problema di sopra //
            
            // Exclude cells in cone
            if(TMath::Abs(icol-colTrig) < sqrSize) {
              //printf("TMath::Abs(icol-colTrig) %i < sqrSize %i\n",TMath::Abs(icol-colTrig) ,sqrSize);continue ;
            }
            if(TMath::Abs(irow-rowTrig) < sqrSize) {
              //printf("TMath::Abs(irow-rowTrig) %i < sqrSize %i\n",TMath::Abs(irow-rowTrig) ,sqrSize);continue ;
            }
            
            Int_t iabsId = eGeom->GetAbsCellIdFromCellIndexes(inSupMod,irowLoc,icolLoc);
            if(!eGeom->CheckAbsCellId(iabsId)) {printf("!eGeom->CheckAbsCellId(iabsId=%i) inSupMod %i irowLoc %i icolLoc %i \n",iabsId,inSupMod, irowLoc, icolLoc);continue;}
            phiBandPtSumCells += cells->GetCellAmplitude(iabsId);
            fhPhiBandCell->Fill(colTrig,rowTrig);
            //printf("inSupMod %i,irowLoc %i,icolLoc %i, iabsId %i, phiBandPtSumCells %f\n",nSupMod,irowLoc,icolLoc,iabsId,phiBandPtSumCells);
          }
        }
      }
    }
  }
  
  Float_t ptTrig = pCandidate->Pt();
  
  fhConeSumPtEtaBandUECell          ->Fill(ptTrig ,        etaBandPtSumCells);
  fhConeSumPtPhiBandUECell          ->Fill(ptTrig ,        phiBandPtSumCells);
  fhConeSumPtEtaBandUECellTrigEtaPhi->Fill(etaTrig,phiTrig,etaBandPtSumCells);
  fhConeSumPtPhiBandUECellTrigEtaPhi->Fill(etaTrig,phiTrig,phiBandPtSumCells);
  
}

//________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateTrackUEBand(AliAODPWG4ParticleCorrelation * pCandidate,
                                                   Float_t & etaBandPtSum, Float_t & phiBandPtSum)
{
  // Get the track pT or sum of pT in phi/eta bands or at 45 degrees from trigger
  
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyNeutral ) return ;
  
  Float_t conesize   = GetIsolationCut()->GetConeSize();
  
  Double_t sumptPerp= 0. ;
  Float_t ptTrig    = pCandidate->Pt() ;
  Float_t phiTrig   = pCandidate->Phi();
  Float_t etaTrig   = pCandidate->Eta();
  
  TObjArray * trackList   = GetCTSTracks() ;
  for(Int_t itrack=0; itrack < trackList->GetEntriesFast(); itrack++)
  {
    AliVTrack* track = (AliVTrack *) trackList->At(itrack);
    
    if(!track)
    {
      printf("AliAnaParticleIsolation::CalculateTrackUEBand() - Track not available?");
      continue;
    }
    
    //Do not count the candidate (pion, conversion photon) or the daughters of the candidate
    if(track->GetID() == pCandidate->GetTrackLabel(0) || track->GetID() == pCandidate->GetTrackLabel(1) ||
       track->GetID() == pCandidate->GetTrackLabel(2) || track->GetID() == pCandidate->GetTrackLabel(3)   ) continue ;
    
    // histo of eta:phi for all tracks
    fhEtaPhiTrack->Fill(track->Eta(),track->Phi());
    
    //exclude particles in cone
    Float_t rad = GetIsolationCut()->Radius(etaTrig, phiTrig, track->Eta(), track->Phi());
    if(rad < conesize) {
    	// histo of eta:phi for all tracks in cone
      fhEtaPhiInConeTrack->Fill(track->Eta(),track->Phi());
      continue ;
    }
    
    //fill histogram for UE in phi band
    if(track->Eta() > (etaTrig-conesize) && track->Eta()  < (etaTrig+conesize))
    {
      phiBandPtSum+=track->Pt();
      fhPhiBandTrack->Fill(track->Eta(),track->Phi());
    }
    
    //fill histogram for UE in eta band in EMCal acceptance
    if(track->Phi() > (phiTrig-conesize) && track->Phi() < (phiTrig+conesize))
    {
      etaBandPtSum+=track->Pt();
      fhEtaBandTrack->Fill(track->Eta(),track->Phi());
    }
    
    //fill the histograms at +-45 degrees in phi from trigger particle, perpedicular to trigger axis in phi
    Double_t dPhi = phiTrig - track->Phi() + TMath::PiOver2();
    Double_t dEta = etaTrig - track->Eta();
    Double_t arg  = dPhi*dPhi + dEta*dEta;
    if(TMath::Sqrt(arg) < conesize)
    {
      fhPtInPerpCone->Fill(ptTrig,TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py()));
      sumptPerp+=track->Pt();
    }
    
    dPhi = phiTrig - track->Phi() - TMath::PiOver2();
    arg  = dPhi*dPhi + dEta*dEta;
    if(TMath::Sqrt(arg) < conesize)
    {
      fhPtInPerpCone->Fill(ptTrig,TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py()));
      sumptPerp+=track->Pt();
    }
  }
  
  fhPerpConeSumPt                    ->Fill(ptTrig ,        sumptPerp   );
  fhConeSumPtEtaBandUETrack          ->Fill(ptTrig ,        etaBandPtSum);
  fhConeSumPtPhiBandUETrack          ->Fill(ptTrig ,        phiBandPtSum);
  fhConeSumPtEtaBandUETrackTrigEtaPhi->Fill(etaTrig,phiTrig,etaBandPtSum);
  fhConeSumPtPhiBandUETrackTrigEtaPhi->Fill(etaTrig,phiTrig,phiBandPtSum);
  
}



//_____________________________________________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateNormalizeUEBandPerUnitArea(AliAODPWG4ParticleCorrelation * pCandidate, Float_t coneptsumCluster,
                                                                  Float_t coneptsumCell,          Float_t coneptsumTrack,
                                                                  Float_t &etaBandptsumTrackNorm, Float_t &etaBandptsumClusterNorm)
{
  //normalize phi/eta band per area unit
  
  Float_t etaUEptsumTrack   = 0 ;
  Float_t phiUEptsumTrack   = 0 ;
  Float_t etaUEptsumCluster = 0 ;
  Float_t phiUEptsumCluster = 0 ;
  Float_t etaUEptsumCell    = 0 ;
  Float_t phiUEptsumCell    = 0 ;
  
  Int_t   partTypeInCone    = GetIsolationCut()->GetParticleTypeInCone();
  
  // Do the normalization
  
  Float_t conesize  = GetIsolationCut()->GetConeSize();
  Float_t coneA     = conesize*conesize*TMath::Pi(); // A = pi R^2, isolation cone area
  Float_t ptTrig    = pCandidate->Pt() ;
  Float_t phiTrig   = pCandidate->Phi();
  Float_t etaTrig   = pCandidate->Eta();
  
  
  // ------ //
  // Tracks //
  // ------ //
  Float_t phiUEptsumTrackNorm  = 0 ;
  Float_t etaUEptsumTrackNorm  = 0 ;
  Float_t coneptsumTrackSubPhi = 0 ;
  Float_t coneptsumTrackSubEta = 0 ;
  Float_t coneptsumTrackSubPhiNorm = 0 ;
  Float_t coneptsumTrackSubEtaNorm = 0 ;
  etaBandptsumTrackNorm = 0 ;
  
  if( partTypeInCone!=AliIsolationCut::kOnlyNeutral )
  {
    // Sum the pT in the phi or eta band for clusters or tracks
    CalculateTrackUEBand   (pCandidate,etaUEptsumTrack  ,phiUEptsumTrack  );// rajouter ici l'histo eta phi
    
    //Fill histos
    fhConeSumPtVSUETracksEtaBand->Fill(coneptsumTrack,etaUEptsumTrack);
    fhConeSumPtVSUETracksPhiBand->Fill(coneptsumTrack,phiUEptsumTrack);
    
    
    Float_t correctConeSumTrack    = 1;
    Float_t correctConeSumTrackPhi = 1;
    
    GetIsolationCut()->CalculateUEBandTrackNormalization(GetReader(),etaTrig, phiTrig,
                                                         phiUEptsumTrack,etaUEptsumTrack,
                                                         phiUEptsumTrackNorm,etaUEptsumTrackNorm,
                                                         correctConeSumTrack,correctConeSumTrackPhi);
    
    coneptsumTrackSubPhi = coneptsumTrack - phiUEptsumTrackNorm;
    coneptsumTrackSubEta = coneptsumTrack - etaUEptsumTrackNorm;
    
    etaBandptsumTrackNorm = etaUEptsumTrackNorm;
    
    fhConeSumPtPhiUESubTrack           ->Fill(ptTrig ,          coneptsumTrackSubPhi);
    fhConeSumPtPhiUESubTrackTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumTrackSubPhi);
    fhConeSumPtEtaUESubTrack           ->Fill(ptTrig ,          coneptsumTrackSubEta);
    fhConeSumPtEtaUESubTrackTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumTrackSubEta);
    
    fhFractionTrackOutConeEta          ->Fill(ptTrig ,         correctConeSumTrack-1);
    fhFractionTrackOutConeEtaTrigEtaPhi->Fill(etaTrig, phiTrig,correctConeSumTrack-1);
    
    if(coneptsumTrack > 0)
    {
      coneptsumTrackSubPhiNorm = coneptsumTrackSubPhi/coneptsumTrack;
    	coneptsumTrackSubEtaNorm = coneptsumTrackSubEta/coneptsumTrack;
    }
    
    fhConeSumPtSubvsConeSumPtTotPhiTrack    ->Fill(coneptsumTrack,coneptsumTrackSubPhi);
    fhConeSumPtSubNormvsConeSumPtTotPhiTrack->Fill(coneptsumTrack,coneptsumTrackSubPhiNorm);
    fhConeSumPtSubvsConeSumPtTotEtaTrack    ->Fill(coneptsumTrack,coneptsumTrackSubEta);
    fhConeSumPtSubNormvsConeSumPtTotEtaTrack->Fill(coneptsumTrack,coneptsumTrackSubEtaNorm);
    
  }
  
  // ------------------------ //
  // EMCal Clusters and cells //
  // ------------------------ //
  Float_t phiUEptsumClusterNorm  = 0 ;
  Float_t etaUEptsumClusterNorm  = 0 ;
  Float_t coneptsumClusterSubPhi = 0 ;
  Float_t coneptsumClusterSubEta = 0 ;
  Float_t coneptsumClusterSubPhiNorm = 0 ;
  Float_t coneptsumClusterSubEtaNorm = 0 ;
  Float_t phiUEptsumCellNorm     = 0 ;
  Float_t etaUEptsumCellNorm     = 0 ;
  Float_t coneptsumCellSubPhi    = 0 ;
  Float_t coneptsumCellSubEta    = 0 ;
  Float_t coneptsumCellSubPhiNorm = 0 ;
  Float_t coneptsumCellSubEtaNorm = 0 ;
  etaBandptsumClusterNorm = 0;
  
  if( partTypeInCone!=AliIsolationCut::kOnlyCharged )
  {
    
    // -------------- //
    // EMCal clusters //
    // -------------- //
    
    // Sum the pT in the phi or eta band for clusters or tracks
    CalculateCaloUEBand    (pCandidate,etaUEptsumCluster,phiUEptsumCluster);// rajouter ici l'histo eta phi
    
    //Fill histos
    fhConeSumPtVSUEClusterEtaBand->Fill(coneptsumCluster,etaUEptsumCluster);
    fhConeSumPtVSUEClusterPhiBand->Fill(coneptsumCluster,phiUEptsumCluster);
    
    
    Float_t correctConeSumClusterEta = 1;
    Float_t correctConeSumClusterPhi = 1;
    
    GetIsolationCut()->CalculateUEBandClusterNormalization(GetReader(),etaTrig, phiTrig,
                                                           phiUEptsumCluster,etaUEptsumCluster,
                                                           phiUEptsumClusterNorm,etaUEptsumClusterNorm,
                                                           correctConeSumClusterEta,correctConeSumClusterPhi);
    
    // In case that cone is out of eta and phi side, we are over correcting, not too often with the current cuts ...
    // Comment if not used
    //  Float_t coneBadCellsCoeff   =1;
    //  Float_t etaBandBadCellsCoeff=1;
    //  Float_t phiBandBadCellsCoeff=1;
    //  GetIsolationCut()->GetCoeffNormBadCell(pCandidate,   GetReader(),coneBadCellsCoeff,etaBandBadCellsCoeff,phiBandBadCellsCoeff) ;
    
    //coneptsumCluster=coneptsumCluster*coneBadCellsCoeff*correctConeSumClusterEta*correctConeSumClusterPhi;
    
    coneptsumClusterSubPhi = coneptsumCluster - phiUEptsumClusterNorm;
    coneptsumClusterSubEta = coneptsumCluster - etaUEptsumClusterNorm;
    
    etaBandptsumClusterNorm = etaUEptsumClusterNorm;
    
    fhConeSumPtPhiUESubCluster           ->Fill(ptTrig ,          coneptsumClusterSubPhi);
    fhConeSumPtPhiUESubClusterTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumClusterSubPhi);
    fhConeSumPtEtaUESubCluster           ->Fill(ptTrig ,          coneptsumClusterSubEta);
    fhConeSumPtEtaUESubClusterTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumClusterSubEta);
    
    fhFractionClusterOutConeEta          ->Fill(ptTrig ,          correctConeSumClusterEta-1);
    fhFractionClusterOutConeEtaTrigEtaPhi->Fill(etaTrig, phiTrig, correctConeSumClusterEta-1);
    fhFractionClusterOutConePhi          ->Fill(ptTrig ,          correctConeSumClusterPhi-1);
    fhFractionClusterOutConePhiTrigEtaPhi->Fill(etaTrig, phiTrig, correctConeSumClusterPhi-1);
    
    if(coneptsumCluster!=0)
    {
	    coneptsumClusterSubPhiNorm = coneptsumClusterSubPhi/coneptsumCluster;
    	coneptsumClusterSubEtaNorm = coneptsumClusterSubEta/coneptsumCluster;
    }
    
    fhConeSumPtSubvsConeSumPtTotPhiCluster    ->Fill(coneptsumCluster,coneptsumClusterSubPhi);
    fhConeSumPtSubNormvsConeSumPtTotPhiCluster->Fill(coneptsumCluster,coneptsumClusterSubPhiNorm);
    fhConeSumPtSubvsConeSumPtTotEtaCluster    ->Fill(coneptsumCluster,coneptsumClusterSubEta);
    fhConeSumPtSubNormvsConeSumPtTotEtaCluster->Fill(coneptsumCluster,coneptsumClusterSubEtaNorm);
    
    // ----------- //
    // EMCal Cells //
    // ----------- //
    
    if(fFillCellHistograms)
    {
      // Sum the pT in the phi or eta band for clusters or tracks
      CalculateCaloCellUEBand(pCandidate,etaUEptsumCell   ,phiUEptsumCell   );
      
      // Move to AliIsolationCut the calculation not the histograms??
      
      //Careful here if EMCal limits changed .. 2010 (4 SM) to 2011-12 (10 SM), for the moment consider 100 deg in phi
      Float_t emcEtaSize = 0.7*2; // TO FIX
      Float_t emcPhiSize = TMath::DegToRad()*100.; // TO FIX
      
      if(((2*conesize*emcPhiSize)-coneA)!=0)phiUEptsumCellNorm = phiUEptsumCell*(coneA / ((2*conesize*emcPhiSize)-coneA));
      if(((2*conesize*emcEtaSize)-coneA)!=0)etaUEptsumCellNorm = etaUEptsumCell*(coneA / ((2*conesize*emcEtaSize)-coneA));
      
      // Need to correct coneptsumCluster by the fraction of the cone out of the calorimeter cut acceptance!
      
      Float_t correctConeSumCellEta = 1;
      if(TMath::Abs(etaTrig)+conesize > emcEtaSize/2.)
      {
        Float_t excess = TMath::Abs(etaTrig) + conesize - emcEtaSize/2.;
        correctConeSumCellEta = GetIsolationCut()->CalculateExcessAreaFraction(excess);
        //printf("Excess EMC-Eta %2.3f, coneA %2.2f,  excessA %2.2f, angle %2.2f,factor %2.2f\n",excess,coneA, excessA, angle*TMath::RadToDeg(), correctConeSumClusterEta);
        // Need to correct phi band surface if part of the cone falls out of track cut acceptance!
        if(((2*(conesize-excess)*emcPhiSize)-(coneA-correctConeSumCellEta))!=0)phiUEptsumCellNorm = phiUEptsumCell*(coneA / ((2*(conesize-excess)*emcPhiSize)-(coneA-correctConeSumCellEta)));
      }
      
      Float_t correctConeSumCellPhi = 1;
      //printf("EMCPhiTrig %2.2f, conesize %2.2f, sum %2.2f, rest %2.2f \n",phiTrig*TMath::RadToDeg(),conesize*TMath::RadToDeg(),(phiTrig+conesize)*TMath::RadToDeg(),(phiTrig-conesize)*TMath::RadToDeg() );
      if((phiTrig+conesize > 180*TMath::DegToRad()) ||
         (phiTrig-conesize <  80*TMath::DegToRad()))
      {
        Float_t excess = 0;
        if( phiTrig+conesize > 180*TMath::DegToRad() ) excess = conesize + phiTrig - 180*TMath::DegToRad() ;
        else                                           excess = conesize - phiTrig +  80*TMath::DegToRad() ;
        
        correctConeSumCellPhi = GetIsolationCut()->CalculateExcessAreaFraction(excess);
        //printf("Excess EMC-Phi %2.3f, coneA %2.2f,  excessA %2.2f, angle %2.2f,factor %2.2f\n",excess,coneA, excessA, angle*TMath::RadToDeg(), correctConeSumClusterPhi);
        
        // Need to correct eta band surface if part of the cone falls out of track cut acceptance!
        if(((2*(conesize-excess)*emcEtaSize)-(coneA-correctConeSumCellPhi))!=0)etaUEptsumCellNorm = etaUEptsumCell*(coneA / ((2*(conesize-excess)*emcEtaSize)-(coneA-correctConeSumCellPhi)));
        
      }
      
      // In case that cone is out of eta and phi side, we are over correcting, not too often with the current cuts ...
      coneptsumCellSubPhi = coneptsumCell*correctConeSumCellEta*correctConeSumCellPhi - phiUEptsumCellNorm;
      coneptsumCellSubEta = coneptsumCell*correctConeSumCellEta*correctConeSumCellPhi - etaUEptsumCellNorm;
      
      fhConeSumPtPhiUESubCell           ->Fill(ptTrig ,          coneptsumCellSubPhi);
      fhConeSumPtPhiUESubCellTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumCellSubPhi);
      fhConeSumPtEtaUESubCell           ->Fill(ptTrig ,          coneptsumCellSubEta);
      fhConeSumPtEtaUESubCellTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumCellSubEta);
      
      fhFractionCellOutConeEta          ->Fill(ptTrig ,          correctConeSumCellEta-1);
      fhFractionCellOutConeEtaTrigEtaPhi->Fill(etaTrig, phiTrig, correctConeSumCellEta-1);
      fhFractionCellOutConePhi          ->Fill(ptTrig ,          correctConeSumCellPhi-1);
      fhFractionCellOutConePhiTrigEtaPhi->Fill(etaTrig, phiTrig, correctConeSumCellPhi-1);
      if(coneptsumCell!=0)
      {
        coneptsumCellSubPhiNorm = coneptsumCellSubPhi/coneptsumCell;
        coneptsumCellSubEtaNorm = coneptsumCellSubEta/coneptsumCell;
      }
      
      fhConeSumPtSubvsConeSumPtTotPhiCell    ->Fill(coneptsumCell,coneptsumCellSubPhi);
      fhConeSumPtSubNormvsConeSumPtTotPhiCell->Fill(coneptsumCell,coneptsumCellSubPhiNorm);
      fhConeSumPtSubvsConeSumPtTotEtaCell    ->Fill(coneptsumCell,coneptsumCellSubEta);
      fhConeSumPtSubNormvsConeSumPtTotEtaCell->Fill(coneptsumCell,coneptsumCellSubEtaNorm);
    }
  }
  
  if( partTypeInCone==AliIsolationCut::kNeutralAndCharged )
  {
    // --------------------------- //
    // Tracks and clusters in cone //
    // --------------------------- //
    
    Double_t sumPhiUESub = coneptsumClusterSubPhi + coneptsumTrackSubPhi;
    Double_t sumEtaUESub = coneptsumClusterSubEta + coneptsumTrackSubEta;
    
    fhConeSumPtPhiUESub          ->Fill(ptTrig ,          sumPhiUESub);
    fhConeSumPtPhiUESubTrigEtaPhi->Fill(etaTrig, phiTrig, sumPhiUESub);
    fhConeSumPtEtaUESub          ->Fill(ptTrig ,          sumEtaUESub);
    fhConeSumPtEtaUESubTrigEtaPhi->Fill(etaTrig, phiTrig, sumEtaUESub);
    
    fhEtaBandClustervsTrack    ->Fill(etaUEptsumCluster    ,etaUEptsumTrack    );
    fhPhiBandClustervsTrack    ->Fill(phiUEptsumCluster    ,phiUEptsumTrack    );
    fhEtaBandNormClustervsTrack->Fill(etaUEptsumClusterNorm,etaUEptsumTrackNorm);
    fhPhiBandNormClustervsTrack->Fill(phiUEptsumClusterNorm,phiUEptsumTrackNorm);
    
    fhConeSumPtEtaUESubClustervsTrack->Fill(coneptsumClusterSubEta,coneptsumTrackSubEta);
    fhConeSumPtPhiUESubClustervsTrack->Fill(coneptsumClusterSubPhi,coneptsumTrackSubPhi);
    
    // ------------------------ //
    // Tracks and cells in cone //
    // ------------------------ //
    
    if(fFillCellHistograms)
    {
      Double_t sumPhiUESubTrackCell = coneptsumCellSubPhi + coneptsumTrackSubPhi;
      Double_t sumEtaUESubTrackCell = coneptsumCellSubEta + coneptsumTrackSubEta;
      
      fhConeSumPtPhiUESubTrackCell          ->Fill(ptTrig ,          sumPhiUESubTrackCell);
      fhConeSumPtPhiUESubTrackCellTrigEtaPhi->Fill(etaTrig, phiTrig, sumPhiUESubTrackCell);
      fhConeSumPtEtaUESubTrackCell          ->Fill(ptTrig ,          sumEtaUESubTrackCell);
      fhConeSumPtEtaUESubTrackCellTrigEtaPhi->Fill(etaTrig, phiTrig, sumEtaUESubTrackCell);
      
      fhEtaBandCellvsTrack    ->Fill(etaUEptsumCell    ,etaUEptsumTrack    );
      fhPhiBandCellvsTrack    ->Fill(phiUEptsumCell    ,phiUEptsumTrack    );
      fhEtaBandNormCellvsTrack->Fill(etaUEptsumCellNorm,etaUEptsumTrackNorm);
      fhPhiBandNormCellvsTrack->Fill(phiUEptsumCellNorm,phiUEptsumTrackNorm);
      
      fhConeSumPtEtaUESubCellvsTrack->Fill(coneptsumCellSubEta,coneptsumTrackSubEta);
      fhConeSumPtPhiUESubCellvsTrack->Fill(coneptsumCellSubPhi,coneptsumTrackSubPhi);
    }
    
  }
}


//______________________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateCaloSignalInCone(AliAODPWG4ParticleCorrelation * aodParticle,
                                                        Float_t & coneptsumCluster, Float_t & coneptLeadCluster)
{
  // Get the cluster pT or sum of pT in isolation cone
  coneptLeadCluster = 0;
  coneptsumCluster  = 0;
  
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyCharged ) return ;
  
  //Recover reference arrays with clusters and tracks
  TObjArray * refclusters = aodParticle->GetObjArray(GetAODObjArrayName()+"Clusters");
  if(!refclusters) return ;
  
  Float_t ptTrig = aodParticle->Pt();
  
  //Get vertex for cluster momentum calculation
  Double_t vertex[] = {0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    GetReader()->GetVertex(vertex);
  
  TLorentzVector mom ;
  for(Int_t icalo=0; icalo < refclusters->GetEntriesFast(); icalo++)
  {
    AliVCluster* calo = (AliVCluster *) refclusters->At(icalo);
    calo->GetMomentum(mom,vertex) ;//Assume that come from vertex in straight line
    
    fhPtInCone       ->Fill(ptTrig, mom.Pt());
    fhPtClusterInCone->Fill(ptTrig, mom.Pt());
    
    if(fFillPileUpHistograms)
    {
      if(GetReader()->IsPileUpFromSPD())               fhPtInConePileUp[0]->Fill(ptTrig,mom.Pt());
      if(GetReader()->IsPileUpFromEMCal())             fhPtInConePileUp[1]->Fill(ptTrig,mom.Pt());
      if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtInConePileUp[2]->Fill(ptTrig,mom.Pt());
      if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtInConePileUp[3]->Fill(ptTrig,mom.Pt());
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtInConePileUp[4]->Fill(ptTrig,mom.Pt());
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtInConePileUp[5]->Fill(ptTrig,mom.Pt());
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtInConePileUp[6]->Fill(ptTrig,mom.Pt());
    }
    
    if(fFillHighMultHistograms) fhPtInConeCent->Fill(GetEventCentrality(),mom.Pt());
    
    coneptsumCluster+=mom.Pt();
    if(mom.Pt() > coneptLeadCluster) coneptLeadCluster = mom.Pt();
  }
  
  fhConeSumPtCluster ->Fill(ptTrig, coneptsumCluster );
  fhConePtLeadCluster->Fill(ptTrig, coneptLeadCluster);
}

//______________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateCaloCellSignalInCone(AliAODPWG4ParticleCorrelation * aodParticle,
                                                            Float_t & coneptsumCell)
{
  // Get the cell amplityde or sum of amplitudes in isolation cone
  // Mising: Remove signal cells in cone in case the trigger is a cluster!
  
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyCharged ) return ;
  
  Float_t conesize = GetIsolationCut()->GetConeSize();
  
  Float_t  ptTrig  = aodParticle->Pt();
  Float_t  phiTrig = aodParticle->Phi();
  if(phiTrig<0) phiTrig += TMath::TwoPi();
  Float_t  etaTrig = aodParticle->Eta();
  
  if(aodParticle->GetDetector()=="EMCAL")
  {
    AliEMCALGeometry* eGeom = AliEMCALGeometry::GetInstance();
    Int_t absId = -999;
    
    if (eGeom->GetAbsCellIdFromEtaPhi(etaTrig,phiTrig,absId))
    {
      if(!eGeom->CheckAbsCellId(absId)) return ;
      
      // Get absolute (col,row) of trigger particle
      Int_t nSupMod = eGeom->GetSuperModuleNumber(absId);
      Int_t nModule = -1;
      Int_t imEta=-1, imPhi=-1;
      Int_t ieta =-1, iphi =-1;
      
      if (eGeom->GetCellIndex(absId,nSupMod,nModule,imPhi,imEta))
      {
        Int_t iEta=-1, iPhi=-1;
        eGeom->GetCellPhiEtaIndexInSModule(nSupMod,nModule,imPhi,imEta,iphi,ieta);
        
        Int_t colTrig = iEta;
        if (nSupMod % 2) colTrig = AliEMCALGeoParams::fgkEMCALCols + iEta ;
        Int_t rowTrig = iPhi + AliEMCALGeoParams::fgkEMCALRows*int(nSupMod/2);
        
        Int_t sqrSize = int(conesize/0.0143);
        
        AliVCaloCells * cells = GetEMCALCells();
        
        // Loop on cells in cone
        for(Int_t irow = rowTrig-sqrSize; irow < rowTrig+sqrSize; irow++)
        {
          for(Int_t icol = colTrig-sqrSize; icol < colTrig+sqrSize; icol++)
          {
            Int_t inSector = int(irow/AliEMCALGeoParams::fgkEMCALRows);
            if(inSector==5) continue;
            
            Int_t inSupMod = -1;
            Int_t icolLoc  = -1;
            if(icol < AliEMCALGeoParams::fgkEMCALCols)
            {
              inSupMod = 2*inSector + 1;
              icolLoc  = icol;
            }
            else if(icol > AliEMCALGeoParams::fgkEMCALCols - 1)
            {
              inSupMod = 2*inSector;
              icolLoc  = icol-AliEMCALGeoParams::fgkEMCALCols;
            }
            
            Int_t irowLoc  = irow - AliEMCALGeoParams::fgkEMCALRows*inSector ;
            
            Int_t iabsId = eGeom->GetAbsCellIdFromCellIndexes(inSupMod,irowLoc,icolLoc);
            if(!eGeom->CheckAbsCellId(iabsId)) continue;
            
            fhPtCellInCone->Fill(ptTrig, cells->GetCellAmplitude(iabsId));
            coneptsumCell += cells->GetCellAmplitude(iabsId);
          }
        }
      }
    }
  }
  
  fhConeSumPtCell->Fill(ptTrig,coneptsumCell);
  
}

//___________________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateTrackSignalInCone(AliAODPWG4ParticleCorrelation * aodParticle,
                                                         Float_t & coneptsumTrack, Float_t & coneptLeadTrack)
{
  // Get the track pT or sum of pT in isolation cone
  
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyNeutral ) return ;
  
  //Recover reference arrays with clusters and tracks
  TObjArray * reftracks   = aodParticle->GetObjArray(GetAODObjArrayName()+"Tracks");
  if(!reftracks) return ;
  
  Float_t  ptTrig = aodParticle->Pt();
  Double_t bz     = GetReader()->GetInputEvent()->GetMagneticField();
  
  for(Int_t itrack=0; itrack < reftracks->GetEntriesFast(); itrack++)
  {
    AliVTrack* track = (AliVTrack *) reftracks->At(itrack);
    Float_t pTtrack = track->Pt();
    
    fhPtInCone     ->Fill(ptTrig,pTtrack);
    fhPtTrackInCone->Fill(ptTrig,pTtrack);
    
    if(fFillPileUpHistograms)
    {
      ULong_t status = track->GetStatus();
      Bool_t okTOF = ( (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ) ;
      //Double32_t tof = track->GetTOFsignal()*1e-3;
      Int_t trackBC = track->GetTOFBunchCrossing(bz);
      
      if     ( okTOF && trackBC!=0 ) fhPtTrackInConeOtherBC->Fill(ptTrig,pTtrack);
      else if( okTOF && trackBC==0 ) fhPtTrackInConeBC0    ->Fill(ptTrig,pTtrack);
      
      Int_t vtxBC = GetReader()->GetVertexBC();
      if(vtxBC == 0 || vtxBC==AliVTrack::kTOFBCNA) fhPtTrackInConeVtxBC0->Fill(ptTrig,pTtrack);
      
      if(GetReader()->IsPileUpFromSPD())             { fhPtInConePileUp[0]->Fill(ptTrig,pTtrack);
        if(okTOF && trackBC!=0 )                         fhPtTrackInConeOtherBCPileUpSPD->Fill(ptTrig,pTtrack);
        if(okTOF && trackBC==0 )                         fhPtTrackInConeBC0PileUpSPD    ->Fill(ptTrig,pTtrack); }
      if(GetReader()->IsPileUpFromEMCal())             fhPtInConePileUp[1]->Fill(ptTrig,pTtrack);
      if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtInConePileUp[2]->Fill(ptTrig,pTtrack);
      if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtInConePileUp[3]->Fill(ptTrig,pTtrack);
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtInConePileUp[4]->Fill(ptTrig,pTtrack);
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtInConePileUp[5]->Fill(ptTrig,pTtrack);
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtInConePileUp[6]->Fill(ptTrig,pTtrack);
    }
    
    if(fFillHighMultHistograms) fhPtInConeCent->Fill(GetEventCentrality(),pTtrack);
    
    coneptsumTrack+=pTtrack;
    if(pTtrack > coneptLeadTrack) coneptLeadTrack = pTtrack;
  }

  fhConeSumPtTrack ->Fill(ptTrig, coneptsumTrack );
  fhConePtLeadTrack->Fill(ptTrig, coneptLeadTrack);

}

//_________________________________________________________________
void AliAnaParticleIsolation::FillPileUpHistograms(Int_t clusterID)
{
  // Fill some histograms to understand pile-up
  
  if(clusterID < 0 )
  {
    printf("AliAnaParticleIsolation::FillPileUpHistograms(), ID of cluster = %d, not possible! ", clusterID);
    return;
  }
  
  Int_t iclus = -1;
  TObjArray* clusters = 0x0;
  if     (fCalorimeter == "EMCAL") clusters = GetEMCALClusters();
  else if(fCalorimeter == "PHOS" ) clusters = GetPHOSClusters();
  
  Float_t energy = 0;
  Float_t time   = -1000;
  
  if(clusters)
  {
    AliVCluster *cluster = FindCluster(clusters,clusterID,iclus);
    energy = cluster->E();
    time   = cluster->GetTOF()*1e9;
  }
  
  //printf("E %f, time %f\n",energy,time);
  AliVEvent * event = GetReader()->GetInputEvent();
  
  fhTimeENoCut->Fill(energy,time);
  if(GetReader()->IsPileUpFromSPD())     fhTimeESPD     ->Fill(energy,time);
  if(event->IsPileupFromSPDInMultBins()) fhTimeESPDMulti->Fill(energy,time);
  
  if(energy < 8) return; // Fill time figures for high energy clusters not too close to trigger threshold
  
  AliESDEvent* esdEv = dynamic_cast<AliESDEvent*> (event);
  AliAODEvent* aodEv = dynamic_cast<AliAODEvent*> (event);
  
  // N pile up vertices
  Int_t nVerticesSPD    = -1;
  Int_t nVerticesTracks = -1;
  
  if      (esdEv)
  {
    nVerticesSPD    = esdEv->GetNumberOfPileupVerticesSPD();
    nVerticesTracks = esdEv->GetNumberOfPileupVerticesTracks();
    
  }//ESD
  else if (aodEv)
  {
    nVerticesSPD    = aodEv->GetNumberOfPileupVerticesSPD();
    nVerticesTracks = aodEv->GetNumberOfPileupVerticesTracks();
  }//AOD
  
  fhTimeNPileUpVertSPD  ->Fill(time,nVerticesSPD);
  fhTimeNPileUpVertTrack->Fill(time,nVerticesTracks);
  
  //printf("Is SPD %d, Is SPD Multi %d, n spd %d, n track %d\n",
  //       GetReader()->IsPileUpFromSPD(),event->IsPileupFromSPDInMultBins(),nVerticesSPD,nVerticesTracks);
  
  Int_t ncont = -1;
  Float_t z1 = -1, z2 = -1;
  Float_t diamZ = -1;
  for(Int_t iVert=0; iVert<nVerticesSPD;iVert++)
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
    
  }// loop
}

//_____________________________________________________________________________________________________________________
void AliAnaParticleIsolation::FillTrackMatchingShowerShapeControlHistograms(AliAODPWG4ParticleCorrelation  *pCandidate,
                                                                            Int_t mcIndex)
{
  // Fill Track matching and Shower Shape control histograms
  if(!fFillTMHisto &&  !fFillSSHisto) return;
  
  Int_t  clusterID = pCandidate->GetCaloLabel(0) ;
  Int_t  nMaxima   = pCandidate->GetFiducialArea(); // bad name, just place holder for the moment
  Int_t  mcTag     = pCandidate->GetTag() ;
  Bool_t isolated  = pCandidate->IsIsolated();
  
  if(clusterID < 0 )
  {
    printf("AliAnaParticleIsolation::FillTrackMatchingShowerShapeControlHistograms(), ID of cluster = %d, not possible! \n", clusterID);
    return;
  }
  
  Int_t iclus = -1;
  TObjArray* clusters = 0x0;
  if     (fCalorimeter == "EMCAL") clusters = GetEMCALClusters();
  else if(fCalorimeter == "PHOS" ) clusters = GetPHOSClusters();
  
  Float_t energy = pCandidate->E();
  Float_t pt     = pCandidate->Pt();
  
  if(clusters)
  {
    AliVCluster *cluster = FindCluster(clusters,clusterID,iclus);
    
    if(fFillSSHisto)
    {
      fhELambda0 [isolated]->Fill(energy, cluster->GetM02() );
      fhPtLambda0[isolated]->Fill(pt,     cluster->GetM02() );
      fhELambda1 [isolated]->Fill(energy, cluster->GetM20() );
      
      if(IsDataMC())
      {
        if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))
          fhPtLambda0MC[kmcPhoton][isolated]->Fill(pt, cluster->GetM02());
        
        fhPtLambda0MC[mcIndex][isolated]->Fill(pt, cluster->GetM02());
      }
      
      if(fCalorimeter == "EMCAL" &&  GetFirstSMCoveredByTRD() >= 0 &&
         GetModuleNumber(cluster) >= GetFirstSMCoveredByTRD()  )
      {
        fhELambda0TRD [isolated]->Fill(energy, cluster->GetM02() );
        fhPtLambda0TRD[isolated]->Fill(pt    , cluster->GetM02() );
        fhELambda1TRD [isolated]->Fill(energy, cluster->GetM20() );
      }
      
      if(fFillNLMHistograms)
      {
        fhNLocMax[isolated]->Fill(energy,nMaxima);
        if     (nMaxima==1) { fhELambda0LocMax1[isolated]->Fill(energy,cluster->GetM02()); fhELambda1LocMax1[isolated]->Fill(energy,cluster->GetM20()); }
        else if(nMaxima==2) { fhELambda0LocMax2[isolated]->Fill(energy,cluster->GetM02()); fhELambda1LocMax2[isolated]->Fill(energy,cluster->GetM20()); }
        else                { fhELambda0LocMaxN[isolated]->Fill(energy,cluster->GetM02()); fhELambda1LocMaxN[isolated]->Fill(energy,cluster->GetM20()); }
      }
    } // SS histo fill
    
    
    if(fFillTMHisto)
    {
      Float_t dZ  = cluster->GetTrackDz();
      Float_t dR  = cluster->GetTrackDx();
      
      if(cluster->IsEMCAL() && GetCaloUtils()->IsRecalculationOfClusterTrackMatchingOn())
      {
        dR = 2000., dZ = 2000.;
        GetCaloUtils()->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dZ,dR);
      }
      
      //printf("ParticleIsolation: dPhi %f, dEta %f\n",dR,dZ);
      if(fhTrackMatchedDEta[isolated] && TMath::Abs(dR) < 999)
      {
        fhTrackMatchedDEta[isolated]->Fill(energy,dZ);
        fhTrackMatchedDPhi[isolated]->Fill(energy,dR);
        if(energy > 0.5) fhTrackMatchedDEtaDPhi[isolated]->Fill(dZ,dR);
      }
      
      // Check dEdx and E/p of matched clusters
      
      if(TMath::Abs(dZ) < 0.05 && TMath::Abs(dR) < 0.05)
      {
        
        AliVTrack *track = GetCaloUtils()->GetMatchedTrack(cluster, GetReader()->GetInputEvent());
        
        if(track)
        {
          Float_t dEdx = track->GetTPCsignal();
          fhdEdx[isolated]->Fill(cluster->E(), dEdx);
          
          Float_t eOverp = cluster->E()/track->P();
          fhEOverP[isolated]->Fill(cluster->E(),  eOverp);
        }
        //else
        //  printf("AliAnaParticleIsolation::FillTrackMatchingShowerShapeHistograms() - Residual OK but (dR, dZ)= (%2.4f,%2.4f) no track associated WHAT? \n", dR,dZ);
        
        
        if(IsDataMC())
        {
          if ( !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion)  )
          {
            if       ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)      ||
                      GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)       ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 2.5 );
            else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton)    ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 0.5 );
            else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron)  ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 1.5 );
            else                                                                                   fhTrackMatchedMCParticle[isolated]->Fill(energy, 3.5 );
            
          }
          else
          {
            if       ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)      ||
                      GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)       ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 6.5 );
            else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton)    ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 4.5 );
            else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron)  ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 5.5 );
            else                                                                                   fhTrackMatchedMCParticle[isolated]->Fill(energy, 7.5 );
          }
          
        }  // MC
        
      } // match window
      
    }// TM histos fill
    
  } // clusters array available
  
}

//______________________________________________________
TObjString *  AliAnaParticleIsolation::GetAnalysisCuts()
{
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar, buffersize,"--- AliAnaParticleIsolation ---\n") ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"Isolation Cand Detector: %s\n",fIsoDetector.Data()) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fReMakeIC =%d (Flag for reisolation during histogram filling) \n",fReMakeIC) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fMakeSeveralIC=%d (Flag for isolation with several cuts at the same time ) \n",fMakeSeveralIC) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fFillTMHisto=%d (Flag for track matching histograms) \n",fFillTMHisto) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fFillSSHisto=%d (Flag for shower shape histograms) \n",fFillSSHisto) ;
  parList+=onePar ;
  
  if(fMakeSeveralIC)
  {
    snprintf(onePar, buffersize,"fNCones =%d (Number of cone sizes) \n",fNCones) ;
    parList+=onePar ;
    snprintf(onePar, buffersize,"fNPtThresFrac=%d (Flag for isolation with several cuts at the same time ) \n",fNPtThresFrac) ;
    parList+=onePar ;
    
    for(Int_t icone = 0; icone < fNCones ; icone++)
    {
      snprintf(onePar, buffersize,"fConeSizes[%d]=%1.2f (isolation cone size) \n",icone, fConeSizes[icone]) ;
      parList+=onePar ;
    }
    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++)
    {
      snprintf(onePar, buffersize,"fPtThresholds[%d]=%1.2f (isolation pt threshold) \n",ipt, fPtThresholds[ipt]) ;
      parList+=onePar ;
    }
    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++)
    {
      snprintf(onePar, buffersize,"fPtFractions[%d]=%1.2f (isolation pt fraction threshold) \n",ipt, fPtFractions[ipt]) ;
      parList+=onePar ;
    }
    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++)
    {
      snprintf(onePar, buffersize,"fSumPtThresholds[%d]=%1.2f (isolation sum pt threshold) \n",ipt, fSumPtThresholds[ipt]) ;
      parList+=onePar ;
    }
  }
  
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in IC class.
  if(!fMakeSeveralIC)parList += GetIsolationCut()->GetICParametersList() ;
  
  return new TObjString(parList) ;
  
}

//________________________________________________________
TList *  AliAnaParticleIsolation::GetCreateOutputObjects()
{
  // Create histograms to be saved in output file and
  // store them in outputContainer
  TList * outputContainer = new TList() ;
  outputContainer->SetName("IsolatedParticleHistos") ;
  
  Int_t   nptbins  = GetHistogramRanges()->GetHistoPtBins();
  Int_t   nphibins = GetHistogramRanges()->GetHistoPhiBins();
  Int_t   netabins = GetHistogramRanges()->GetHistoEtaBins();
  Float_t ptmax    = GetHistogramRanges()->GetHistoPtMax();
  Float_t phimax   = GetHistogramRanges()->GetHistoPhiMax();
  Float_t etamax   = GetHistogramRanges()->GetHistoEtaMax();
  Float_t ptmin    = GetHistogramRanges()->GetHistoPtMin();
  Float_t phimin   = GetHistogramRanges()->GetHistoPhiMin();
  Float_t etamin   = GetHistogramRanges()->GetHistoEtaMin();
  Int_t   ssbins   = GetHistogramRanges()->GetHistoShowerShapeBins();
  Float_t ssmax    = GetHistogramRanges()->GetHistoShowerShapeMax();
  Float_t ssmin    = GetHistogramRanges()->GetHistoShowerShapeMin();
  Int_t   ntimebins= GetHistogramRanges()->GetHistoTimeBins();
  Float_t timemax  = GetHistogramRanges()->GetHistoTimeMax();
  Float_t timemin  = GetHistogramRanges()->GetHistoTimeMin();
  
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
  
  Int_t   nptsumbins    = GetHistogramRanges()->GetHistoNPtSumBins();
  Float_t ptsummax      = GetHistogramRanges()->GetHistoPtSumMax();
  Float_t ptsummin      = GetHistogramRanges()->GetHistoPtSumMin();
  Int_t   nptinconebins = GetHistogramRanges()->GetHistoNPtInConeBins();
  Float_t ptinconemax   = GetHistogramRanges()->GetHistoPtInConeMax();
  Float_t ptinconemin   = GetHistogramRanges()->GetHistoPtInConeMin();
  
  //Float_t ptthre    = GetIsolationCut()->GetPtThreshold();
  //Float_t ptsumthre = GetIsolationCut()->GetSumPtThreshold();
  //Float_t ptfrac    = GetIsolationCut()->GetPtFraction();
  Float_t r         = GetIsolationCut()->GetConeSize();
  Int_t   method    = GetIsolationCut()->GetICMethod() ;
  Int_t   particle  = GetIsolationCut()->GetParticleTypeInCone() ;
  
  TString sThreshold = "";
  if      ( method == AliIsolationCut::kSumPtIC )
  {
    sThreshold = Form(", %2.2f < #Sigma #it{p}_{T}^{in cone} < %2.2f GeV/#it{c}",
                      GetIsolationCut()->GetSumPtThreshold(), GetIsolationCut()->GetSumPtThresholdMax());
    if(GetIsolationCut()->GetSumPtThresholdMax() > 200)
      sThreshold = Form(", #Sigma #it{p}_{T}^{in cone} = %2.2f GeV/#it{c}",
                        GetIsolationCut()->GetSumPtThreshold());
  }
  else if ( method == AliIsolationCut::kPtThresIC)
  {
    sThreshold = Form(", %2.2f < #it{p}_{T}^{th} < %2.2f GeV/#it{c}",
                      GetIsolationCut()->GetPtThreshold(),GetIsolationCut()->GetPtThresholdMax());
    if(GetIsolationCut()->GetSumPtThreshold() > 200)
      sThreshold = Form(", #it{p}_{T}^{th} = %2.2f GeV/#it{c}",
                        GetIsolationCut()->GetPtThreshold());
  }
  else if ( method == AliIsolationCut::kPtFracIC)
    sThreshold = Form(", #Sigma #it{p}_{T}^{in cone}/#it{p}_{T}^{trig} = %2.2f" ,
                      GetIsolationCut()->GetPtFraction());
  
  TString sParticle = ", x^{0,#pm}";
  if      ( particle == AliIsolationCut::kOnlyNeutral )  sParticle = ", x^{0}";
  else if ( particle == AliIsolationCut::kOnlyCharged )  sParticle = ", x^{#pm}";
  
  TString parTitle = Form("#it{R} = %2.2f%s%s",GetIsolationCut()->GetConeSize(), sThreshold.Data(),sParticle.Data());
  
  TString pileUpName[] = {"SPD","EMCAL","SPDOrEMCAL","SPDAndEMCAL","SPDAndNotEMCAL","EMCALAndNotSPD","NotSPDAndNotEMCAL"} ;
  
  // MC histograms title and name
  TString mcPartType[] = { "#gamma", "#gamma_{prompt}", "#gamma_{fragmentation}",
    "#pi^{0} (merged #gamma)","#gamma_{#pi decay}",
    "#gamma_{#eta decay}","#gamma_{other decay}",
    "e^{#pm}","hadrons?"} ;
  
  TString mcPartName[] = { "Photon","PhotonPrompt","PhotonFrag",
    "Pi0","Pi0Decay","EtaDecay","OtherDecay",
    "Electron","Hadron"} ;
  
  // Primary MC histograms title and name
  TString pptype[] = { "#gamma", "#gamma_{#pi decay}","#gamma_{other decay}",
    "#gamma_{prompt}","#gamma_{fragmentation}","#gamma_{ISR}"} ;
  
  TString ppname[] = { "Photon","PhotonPi0Decay","PhotonOtherDecay",
    "PhotonPrompt","PhotonFrag","PhotonISR"} ;
  
  // Not Isolated histograms, reference histograms
  
  fhENoIso  = new TH1F("hENoIso",
                       Form("Number of not isolated leading particles vs #it{p}_{T}, %s",parTitle.Data()),
                       nptbins,ptmin,ptmax);
  fhENoIso->SetYTitle("#it{counts}");
  fhENoIso->SetXTitle("E (GeV/#it{c})");
  outputContainer->Add(fhENoIso) ;
  
  fhPtNoIso  = new TH1F("hPtNoIso",
                        Form("Number of not isolated leading particles vs #it{p}_{T}, %s",parTitle.Data()),
                        nptbins,ptmin,ptmax);
  fhPtNoIso->SetYTitle("#it{counts}");
  fhPtNoIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtNoIso) ;
  
  fhEtaPhiNoIso  = new TH2F("hEtaPhiNoIso",
                            Form("Number of not isolated leading particles #eta vs #phi, %s",parTitle.Data()),
                            netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiNoIso->SetXTitle("#eta");
  fhEtaPhiNoIso->SetYTitle("#phi");
  outputContainer->Add(fhEtaPhiNoIso) ;
  
  if(IsDataMC())
  {
    // For histograms in arrays, index in the array, corresponding to any particle origin
    
    for(Int_t imc = 0; imc < 9; imc++)
    {
      
      fhPtNoIsoMC[imc]  = new TH1F(Form("hPtNoIsoMC%s",mcPartName[imc].Data()),
                                   Form("#it{p}_{T} of NOT isolated %s, %s",mcPartType[imc].Data(),parTitle.Data()),
                                   nptbins,ptmin,ptmax);
      fhPtNoIsoMC[imc]->SetYTitle("#it{counts}");
      fhPtNoIsoMC[imc]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
      outputContainer->Add(fhPtNoIsoMC[imc]) ;
      
      fhPtIsoMC[imc]  = new TH1F(Form("hPtMC%s",mcPartName[imc].Data()),
                                 Form("#it{p}_{T} of isolated %s, %s",mcPartType[imc].Data(),parTitle.Data()),
                                 nptbins,ptmin,ptmax);
      fhPtIsoMC[imc]->SetYTitle("#it{counts}");
      fhPtIsoMC[imc]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
      outputContainer->Add(fhPtIsoMC[imc]) ;
      
      fhPhiIsoMC[imc]  = new TH2F(Form("hPhiMC%s",mcPartName[imc].Data()),
                                  Form("#phi vs #it{p}_{T} of isolated %s, %s",mcPartType[imc].Data(),parTitle.Data()),
                                  nptbins,ptmin,ptmax,nphibins,phimin,phimax);
      fhPhiIsoMC[imc]->SetYTitle("#phi");
      fhPhiIsoMC[imc]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
      outputContainer->Add(fhPhiIsoMC[imc]) ;
      
      fhEtaIsoMC[imc]  = new TH2F(Form("hEtaMC%s",mcPartName[imc].Data()),
                                  Form("#phi vs #it{p}_{T} of isolated %s, %s",mcPartType[imc].Data(),parTitle.Data()),
                                  nptbins,ptmin,ptmax,netabins,etamin,etamax);
      fhEtaIsoMC[imc]->SetYTitle("#eta");
      fhEtaIsoMC[imc]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
      outputContainer->Add(fhEtaIsoMC[imc]) ;
    }
  }
  
  // Histograms for tagged candidates as decay
  if(fFillTaggedDecayHistograms)
  {
    for(Int_t ibit = 0; ibit < fNDecayBits; ibit++)
    {
      fhPtDecayNoIso[ibit]  =
      new TH1F(Form("hPtDecayNoIso_bit%d",fDecayBits[ibit]),
               Form("Number of not isolated leading pi0 decay particles vs #it{p}_{T}, bit %d, %s",fDecayBits[ibit],parTitle.Data()),
               nptbins,ptmin,ptmax);
      fhPtDecayNoIso[ibit]->SetYTitle("#it{counts}");
      fhPtDecayNoIso[ibit]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtDecayNoIso[ibit]) ;
      
      fhEtaPhiDecayNoIso[ibit]  =
      new TH2F(Form("hEtaPhiDecayNoIso_bit%d",fDecayBits[ibit]),
               Form("Number of not isolated leading Pi0 decay particles #eta vs #phi, bit %d, %s",fDecayBits[ibit],parTitle.Data()),
               netabins,etamin,etamax,nphibins,phimin,phimax);
      fhEtaPhiDecayNoIso[ibit]->SetXTitle("#eta");
      fhEtaPhiDecayNoIso[ibit]->SetYTitle("#phi");
      outputContainer->Add(fhEtaPhiDecayNoIso[ibit]) ;
      
      if(!fMakeSeveralIC)
      {
        fhPtDecayIso[ibit]  =
        new TH1F(Form("hPtDecayIso_bit%d",fDecayBits[ibit]),
                 Form("Number of isolated #pi^{0} decay particles vs #it{p}_{T}, bit %d, %s",fDecayBits[ibit],parTitle.Data()),
                 nptbins,ptmin,ptmax);
        fhPtDecayIso[ibit]->SetYTitle("#it{counts}");
        fhPtDecayIso[ibit]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
        outputContainer->Add(fhPtDecayIso[ibit]) ;
        
        fhEtaPhiDecayIso[ibit]  =
        new TH2F(Form("hEtaPhiDecayIso_bit%d",fDecayBits[ibit]),
                 Form("Number of isolated Pi0 decay particles #eta vs #phi, bit %d, %s",fDecayBits[ibit],parTitle.Data()),
                 netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiDecayIso[ibit]->SetXTitle("#eta");
        fhEtaPhiDecayIso[ibit]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiDecayIso[ibit]) ;
      }
      
      if(IsDataMC())
      {
        for(Int_t imc = 0; imc < 9; imc++)
        {
          
          fhPtDecayNoIsoMC[ibit][imc]  =
          new TH1F(Form("hPtDecayNoIso_bit%d_MC%s",fDecayBits[ibit],mcPartName[imc].Data()),
                   Form("#it{p}_{T} of NOT isolated, decay bit %d,  %s, %s",fDecayBits[ibit],mcPartType[imc].Data(),parTitle.Data()),
                   nptbins,ptmin,ptmax);
          fhPtDecayNoIsoMC[ibit][imc]->SetYTitle("#it{counts}");
          fhPtDecayNoIsoMC[ibit][imc]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
          outputContainer->Add(fhPtDecayNoIsoMC[ibit][imc]) ;
          
          if(!fMakeSeveralIC)
          {
            fhPtDecayIsoMC[ibit][imc]  =
            new TH1F(Form("hPtDecay_bit%d_MC%s",fDecayBits[ibit],mcPartName[imc].Data()),
                     Form("#it{p}_{T} of isolated %s, decay bit %d, %s",mcPartType[imc].Data(),fDecayBits[ibit],parTitle.Data()),
                     nptbins,ptmin,ptmax);
            fhPtDecayIsoMC[ibit][imc]->SetYTitle("#it{counts}");
            fhPtDecayIsoMC[ibit][imc]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add(fhPtDecayIsoMC[ibit][imc]) ;
          }
        }// MC particle loop
      }// MC
    } // bit loop
  }// decay
  
  if(!fMakeSeveralIC)
  {
    TString isoName [] = {"NoIso",""};
    TString isoTitle[] = {"Not isolated"  ,"isolated"};
    
    fhEIso   = new TH1F("hE",
                        Form("Number of isolated particles vs E, %s",parTitle.Data()),
                        nptbins,ptmin,ptmax);
    fhEIso->SetYTitle("d#it{N} / d#it{E}");
    fhEIso->SetXTitle("#it{E} (GeV/#it{c})");
    outputContainer->Add(fhEIso) ;
    
    fhPtIso  = new TH1F("hPt",
                        Form("Number of isolated particles vs #it{p}_{T}, %s",parTitle.Data()),
                        nptbins,ptmin,ptmax);
    fhPtIso->SetYTitle("d#it{N} / #it{p}_{T}");
    fhPtIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtIso) ;
    
    fhPhiIso  = new TH2F("hPhi",
                         Form("Number of isolated particles vs #phi, %s",parTitle.Data()),
                         nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhPhiIso->SetYTitle("#phi");
    fhPhiIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPhiIso) ;
    
    fhEtaIso  = new TH2F("hEta",
                         Form("Number of isolated particles vs #eta, %s",parTitle.Data()),
                         nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhEtaIso->SetYTitle("#eta");
    fhEtaIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhEtaIso) ;
    
    fhEtaPhiIso  = new TH2F("hEtaPhiIso",
                            Form("Number of isolated particles #eta vs #phi, %s",parTitle.Data()),
                            netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiIso->SetXTitle("#eta");
    fhEtaPhiIso->SetYTitle("#phi");
    outputContainer->Add(fhEtaPhiIso) ;
    
    if(fFillHighMultHistograms)
    {
      fhPtCentralityIso  = new TH2F("hPtCentrality",
                                    Form("centrality vs #it{p}_{T} for isolated particles, %s",parTitle.Data()),
                                    nptbins,ptmin,ptmax, 100,0,100);
      fhPtCentralityIso->SetYTitle("centrality");
      fhPtCentralityIso->SetXTitle("#it{p}_{T}(GeV/#it{c})");
      outputContainer->Add(fhPtCentralityIso) ;
      
      fhPtEventPlaneIso  = new TH2F("hPtEventPlane",
                                    Form("event plane angle vs #it{p}_{T} for isolated particles, %s",parTitle.Data()),
                                    nptbins,ptmin,ptmax, 100,0,TMath::Pi());
      fhPtEventPlaneIso->SetYTitle("Event plane angle (rad)");
      fhPtEventPlaneIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtEventPlaneIso) ;
    }
    
    if(fFillNLMHistograms)
    {
      fhPtNLocMaxIso  = new TH2F("hPtNLocMax",
                                 Form("Number of isolated particles vs #it{p}_{T}, %s",parTitle.Data()),
                                 nptbins,ptmin,ptmax,10,0,10);
      fhPtNLocMaxIso->SetYTitle("#it{NLM}");
      fhPtNLocMaxIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      
      fhPtNLocMaxNoIso  = new TH2F("hPtNLocMaxNoIso",
                                   Form("Number of not isolated particles vs #it{p}_{T}, %s",parTitle.Data()),
                                   nptbins,ptmin,ptmax,10,0,10);
      fhPtNLocMaxNoIso->SetYTitle("#it{NLM}");
      fhPtNLocMaxNoIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtNLocMaxNoIso) ;
    }

    fhConePtLead  = new TH2F("hConePtLead",
                            Form("Track or Cluster  leading #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                            nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
    fhConePtLead->SetYTitle("#it{p}_{T, leading} (GeV/#it{c})");
    fhConePtLead->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
    outputContainer->Add(fhConePtLead) ;

    
    fhConeSumPt  = new TH2F("hConePtSum",
                            Form("Track and Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                            nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
    fhConeSumPt->SetYTitle("#Sigma #it{p}_{T}");
    fhConeSumPt->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPt) ;
    
    fhConeSumPtTrigEtaPhi  = new TH2F("hConePtSumTrigEtaPhi",
                                      Form("Trigger #eta vs #phi, #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                      netabins,etamin,etamax,nphibins,phimin,phimax);
    fhConeSumPtTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
    fhConeSumPtTrigEtaPhi->SetXTitle("#eta_{trigger}");
    fhConeSumPtTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
    outputContainer->Add(fhConeSumPtTrigEtaPhi) ;
    
    fhPtInCone  = new TH2F("hPtInCone",
                           Form("#it{p}_{T} of clusters and tracks in isolation cone for #it{R} =  %2.2f",r),
                           nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
    fhPtInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
    fhPtInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtInCone) ;
    
    if(fFillHighMultHistograms)
    {
      fhPtInConeCent  = new TH2F("hPtInConeCent",
                                 Form("#it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                 100,0,100,nptinconebins,ptinconemin,ptinconemax);
      fhPtInConeCent->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtInConeCent->SetXTitle("centrality");
      outputContainer->Add(fhPtInConeCent) ;
    }
    
    // Cluster only histograms
    if(GetIsolationCut()->GetParticleTypeInCone()!=AliIsolationCut::kOnlyCharged)
    {
      fhConeSumPtCluster  = new TH2F("hConePtSumCluster",
                                     Form("Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                     nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtCluster->SetYTitle("#Sigma #it{p}_{T}");
      fhConeSumPtCluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtCluster) ;
      
      fhConePtLeadCluster  = new TH2F("hConeLeadPtCluster",
                                    Form("Cluster leading in isolation cone for #it{R} =  %2.2f",r),
                                    nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConePtLeadCluster->SetYTitle("#it{p}_{T, leading} (GeV/#it{c})");
      fhConePtLeadCluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConePtLeadCluster) ;

      
      if(fFillCellHistograms)
      {
        fhConeSumPtCell  = new TH2F("hConePtSumCell",
                                    Form("Cell #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                    nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtCell->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtCell->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtCell) ;
      }
      
      if(fFillUEBandSubtractHistograms)
      {
        fhConeSumPtEtaBandUECluster  = new TH2F("hConePtSumEtaBandUECluster",
                                                "#Sigma cluster #it{p}_{T} in UE Eta Band",
                                                nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtEtaBandUECluster->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaBandUECluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaBandUECluster) ;
        
        fhConeSumPtPhiBandUECluster  = new TH2F("hConePtSumPhiBandUECluster",
                                                "#Sigma cluster #it{p}_{T} UE Phi Band",
                                                nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtPhiBandUECluster->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiBandUECluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiBandUECluster) ;
        
        fhConeSumPtEtaBandUEClusterTrigEtaPhi  = new TH2F("hConePtSumEtaBandUEClusterTrigEtaPhi",
                                                          "Trigger #eta vs #phi, #Sigma cluster #it{p}_{T} in UE Eta Band",
                                                          netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtEtaBandUEClusterTrigEtaPhi) ;
        
        fhConeSumPtPhiBandUEClusterTrigEtaPhi  = new TH2F("hConePtSumPhiBandUEClusterTrigEtaPhi",
                                                          "Trigger #eta vs #phi, #Sigma cluster #it{p}_{T} UE Phi Band",
                                                          netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtPhiBandUEClusterTrigEtaPhi) ;
        if(fFillCellHistograms)
        {
          
          fhConeSumPtEtaBandUECell  = new TH2F("hConePtSumEtaBandUECell",
                                               "#Sigma cell #it{p}_{T} in UE Eta Band",
                                               nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
          fhConeSumPtEtaBandUECell->SetYTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaBandUECell->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtEtaBandUECell) ;
          
          fhConeSumPtPhiBandUECell  = new TH2F("hConePtSumPhiBandUECell",
                                               "#Sigma cell #it{p}_{T} UE Phi Band",
                                               nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
          fhConeSumPtPhiBandUECell->SetYTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiBandUECell->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPhiBandUECell) ;
          
          fhConeSumPtEtaBandUECellTrigEtaPhi  = new TH2F("hConePtSumEtaBandUECellTrigEtaPhi",
                                                         "Trigger #eta vs #phi, #Sigma cell #it{p}_{T} in UE Eta Band",
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaBandUECellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaBandUECellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaBandUECellTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaBandUECellTrigEtaPhi) ;
          
          fhConeSumPtPhiBandUECellTrigEtaPhi  = new TH2F("hConePtSumPhiBandUECellTrigEtaPhi",
                                                         "Trigger #eta vs #phi, #Sigma cell #it{p}_{T} UE Phi Band",
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiBandUECellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiBandUECellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiBandUECellTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtPhiBandUECellTrigEtaPhi) ;
        }
        
        fhEtaBandCluster  = new TH2F("hEtaBandCluster",
                                     Form("#eta vs #phi of clusters in #eta band isolation cone for #it{R} =  %2.2f",r),
                                     netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaBandCluster->SetXTitle("#eta");
        fhEtaBandCluster->SetYTitle("#phi");
        outputContainer->Add(fhEtaBandCluster) ;
        
        fhPhiBandCluster  = new TH2F("hPhiBandCluster",
                                     Form("#eta vs #phi of clusters in #phi band isolation cone for #it{R} =  %2.2f",r),
                                     netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhPhiBandCluster->SetXTitle("#eta");
        fhPhiBandCluster->SetYTitle("#phi");
        outputContainer->Add(fhPhiBandCluster) ;
        
        fhEtaPhiInConeCluster= new TH2F("hEtaPhiInConeCluster",
                                        Form("#eta vs #phi of clusters in cone for #it{R} =  %2.2f",r),
                                        netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaPhiInConeCluster->SetXTitle("#eta");
        fhEtaPhiInConeCluster->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiInConeCluster) ;
        
        fhEtaPhiCluster= new TH2F("hEtaPhiCluster",
                                  Form("#eta vs #phi of all clusters"),
                                  netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaPhiCluster->SetXTitle("#eta");
        fhEtaPhiCluster->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiCluster) ;
        
      }
      
      fhPtClusterInCone  = new TH2F("hPtClusterInCone",
                                    Form("#it{p}_{T} of clusters in isolation cone for #it{R} =  %2.2f",r),
                                    nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtClusterInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtClusterInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtClusterInCone) ;
      
      if(fFillCellHistograms)
      {
        fhPtCellInCone  = new TH2F("hPtCellInCone",
                                   Form("#it{p}_{T} of cells in isolation cone for #it{R} =  %2.2f",r),
                                   nptbins,ptmin,ptmax,1000,0,50);
        fhPtCellInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtCellInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtCellInCone) ;
        
        fhEtaBandCell  = new TH2F("hEtaBandCell",
                                  Form("#col vs #row of cells in #eta band isolation cone for #it{R} =  %2.2f",r),
                                  96,0,95,128,0,127);
        fhEtaBandCell->SetXTitle("#col");
        fhEtaBandCell->SetYTitle("#row");
        outputContainer->Add(fhEtaBandCell) ;
        
        fhPhiBandCell  = new TH2F("hPhiBandCell",
                                  Form("#col vs #row of cells in #phi band isolation cone for #it{R} =  %2.2f",r),
                                  96,0,95,128,0,127);
        fhPhiBandCell->SetXTitle("#col");
        fhPhiBandCell->SetYTitle("#row");
        outputContainer->Add(fhPhiBandCell) ;
      }
      
      if(fFillUEBandSubtractHistograms)
      {
        fhConeSumPtEtaUESubCluster  = new TH2F("hConeSumPtEtaUESubCluster",
                                               Form("Clusters #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                               nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtEtaUESubCluster->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaUESubCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaUESubCluster) ;
        
        fhConeSumPtPhiUESubCluster  = new TH2F("hConeSumPtPhiUESubCluster",
                                               Form("Clusters #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                               nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtPhiUESubCluster->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiUESubCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiUESubCluster) ;
        
        fhConeSumPtEtaUESubClusterTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubClusterTrigEtaPhi",
                                                         Form("Trigger #eta vs #phi, Clusters #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtEtaUESubClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaUESubClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtEtaUESubClusterTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtEtaUESubClusterTrigEtaPhi) ;
        
        fhConeSumPtPhiUESubClusterTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubClusterTrigEtaPhi",
                                                         Form("Trigger #eta vs #phi, Clusters #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtPhiUESubClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiUESubClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtPhiUESubClusterTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtPhiUESubClusterTrigEtaPhi) ;
        
        if(fFillCellHistograms)
        {
          fhConeSumPtEtaUESubCell  = new TH2F("hConeSumPtEtaUESubCell",
                                              Form("Cells #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                              nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtEtaUESubCell->SetYTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaUESubCell->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtEtaUESubCell) ;
          
          fhConeSumPtPhiUESubCell  = new TH2F("hConeSumPtPhiUESubCell",
                                              Form("Cells #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                              nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtPhiUESubCell->SetYTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiUESubCell->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPhiUESubCell) ;
          
          fhConeSumPtEtaUESubCellTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubCellTrigEtaPhi",
                                                        Form("Trigger #eta vs #phi, Cells #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                        netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaUESubCellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaUESubCellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaUESubCellTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaUESubCellTrigEtaPhi) ;
          
          fhConeSumPtPhiUESubCellTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubCellTrigEtaPhi",
                                                        Form("Trigger #eta vs #phi, Cells #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                                        netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiUESubCellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiUESubCellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiUESubCellTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtPhiUESubCellTrigEtaPhi) ;
        }
        
        fhFractionClusterOutConeEta  = new TH2F("hFractionClusterOutConeEta",
                                                Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #eta acceptance",r),
                                                nptbins,ptmin,ptmax,100,0,1);
        fhFractionClusterOutConeEta->SetYTitle("#it{fraction}");
        fhFractionClusterOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
        outputContainer->Add(fhFractionClusterOutConeEta) ;
        
        fhFractionClusterOutConeEtaTrigEtaPhi  = new TH2F("hFractionClusterOutConeEtaTrigEtaPhi",
                                                          Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #eta acceptance, in trigger #eta-#phi ",r),
                                                          netabins,etamin,etamax,nphibins,phimin,phimax);
        fhFractionClusterOutConeEtaTrigEtaPhi->SetZTitle("#it{fraction}");
        fhFractionClusterOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhFractionClusterOutConeEtaTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhFractionClusterOutConeEtaTrigEtaPhi) ;
        
        fhFractionClusterOutConePhi  = new TH2F("hFractionClusterOutConePhi",
                                                Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #phi acceptance",r),
                                                nptbins,ptmin,ptmax,100,0,1);
        fhFractionClusterOutConePhi->SetYTitle("#it{fraction}");
        fhFractionClusterOutConePhi->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
        outputContainer->Add(fhFractionClusterOutConePhi) ;
        
        fhFractionClusterOutConePhiTrigEtaPhi  = new TH2F("hFractionClusterOutConePhiTrigEtaPhi",
                                                          Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #phi acceptance, in trigger #eta-#phi ",r),
                                                          netabins,etamin,etamax,nphibins,phimin,phimax);
        fhFractionClusterOutConePhiTrigEtaPhi->SetZTitle("#it{fraction}");
        fhFractionClusterOutConePhiTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhFractionClusterOutConePhiTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhFractionClusterOutConePhiTrigEtaPhi) ;
        
        fhConeSumPtSubvsConeSumPtTotPhiCluster = new TH2F("hConeSumPtSubvsConeSumPtTotPhiCluster",
                                                          Form("#Sigma #it{p}_{T} in cone after bkg sub from #phi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                          nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubvsConeSumPtTotPhiCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubvsConeSumPtTotPhiCluster->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubvsConeSumPtTotPhiCluster);
        
        fhConeSumPtSubNormvsConeSumPtTotPhiCluster = new TH2F("hConeSumPtSubNormvsConeSumPtTotPhiCluster",
                                                              Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #phi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                              nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubNormvsConeSumPtTotPhiCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubNormvsConeSumPtTotPhiCluster->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotPhiCluster);
        
        fhConeSumPtSubvsConeSumPtTotEtaCluster = new TH2F("hConeSumPtSubvsConeSumPtTotEtaCluster",
                                                          Form("#Sigma #it{p}_{T} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                          nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubvsConeSumPtTotEtaCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubvsConeSumPtTotEtaCluster->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubvsConeSumPtTotEtaCluster);
        
        fhConeSumPtSubNormvsConeSumPtTotEtaCluster = new TH2F("hConeSumPtSubNormvsConeSumPtTotEtaCluster",
                                                              Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                              nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubNormvsConeSumPtTotEtaCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubNormvsConeSumPtTotEtaCluster->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotEtaCluster);
        
        fhConeSumPtVSUEClusterEtaBand  = new TH2F("hConeSumPtVSUEClusterEtaBand",
                                                  Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in eta band for cluster (before normalization), R=%2.2f",r),
                                                  nptsumbins,ptsummin,ptsummax,2*nptsumbins,ptsummin,2*ptsummax);
        fhConeSumPtVSUEClusterEtaBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
        fhConeSumPtVSUEClusterEtaBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtVSUEClusterEtaBand);
        
        fhConeSumPtVSUEClusterPhiBand  = new TH2F("hConeSumPtVSUEClusterPhiBand",
                                                  Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in phi band for cluster (before normalization), R=%2.2f",r),
                                                  nptsumbins,ptsummin,ptsummax,8*nptsumbins,ptsummin,8*ptsummax);
        fhConeSumPtVSUEClusterPhiBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
        fhConeSumPtVSUEClusterPhiBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtVSUEClusterPhiBand);
        
        if(fFillCellHistograms)
        {
          fhFractionCellOutConeEta  = new TH2F("hFractionCellOutConeEta",
                                               Form("Fraction of the isolation cone #it{R} =  %2.2f, out of cells #eta acceptance",r),
                                               nptbins,ptmin,ptmax,100,0,1);
          fhFractionCellOutConeEta->SetYTitle("#it{fraction}");
          fhFractionCellOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
          outputContainer->Add(fhFractionCellOutConeEta) ;
          
          fhFractionCellOutConeEtaTrigEtaPhi  = new TH2F("hFractionCellOutConeEtaTrigEtaPhi",
                                                         Form("Fraction of the isolation cone #it{R} =  %2.2f, out of cells #eta acceptance, in trigger #eta-#phi ",r),
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhFractionCellOutConeEtaTrigEtaPhi->SetZTitle("#it{fraction}");
          fhFractionCellOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhFractionCellOutConeEtaTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhFractionCellOutConeEtaTrigEtaPhi) ;
          
          fhFractionCellOutConePhi  = new TH2F("hFractionCellOutConePhi",
                                               Form("Fraction of the isolation cone #it{R} =  %2.2f, out of cells #phi acceptance",r),
                                               nptbins,ptmin,ptmax,100,0,1);
          fhFractionCellOutConePhi->SetYTitle("#it{fraction}");
          fhFractionCellOutConePhi->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
          outputContainer->Add(fhFractionCellOutConePhi) ;
          
          fhFractionCellOutConePhiTrigEtaPhi  = new TH2F("hFractionCellOutConePhiTrigEtaPhi",
                                                         Form("Fraction of the isolation cone #it{R} =  %2.2f, out of cells #phi acceptance, in trigger #eta-#phi ",r),
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhFractionCellOutConePhiTrigEtaPhi->SetZTitle("#it{fraction}");
          fhFractionCellOutConePhiTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhFractionCellOutConePhiTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhFractionCellOutConePhiTrigEtaPhi) ;
          
          
          fhConeSumPtSubvsConeSumPtTotPhiCell = new TH2F("hConeSumPtSubvsConeSumPtTotPhiCell",
                                                         Form("#Sigma #it{p}_{T} in cone after bkg sub from #phi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                         nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtSubvsConeSumPtTotPhiCell->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubvsConeSumPtTotPhiCell->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubvsConeSumPtTotPhiCell);
          
          fhConeSumPtSubNormvsConeSumPtTotPhiCell = new TH2F("hConeSumPtSubNormvsConeSumPtTotPhiCell",
                                                             Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #phi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                             nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtSubNormvsConeSumPtTotPhiCell->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubNormvsConeSumPtTotPhiCell->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotPhiCell);
          
          fhConeSumPtSubvsConeSumPtTotEtaCell = new TH2F("hConeSumPtSubvsConeSumPtTotEtaCell",
                                                         Form("#Sigma #it{p}_{T} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                         nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtSubvsConeSumPtTotEtaCell->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubvsConeSumPtTotEtaCell->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubvsConeSumPtTotEtaCell);
          
          fhConeSumPtSubNormvsConeSumPtTotEtaCell = new TH2F("hConeSumPtSubNormvsConeSumPtTotEtaCell",
                                                             Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                             nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtSubNormvsConeSumPtTotEtaCell->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubNormvsConeSumPtTotEtaCell->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotEtaCell);
        }
      }
    }
    
    // Track only histograms
    if(GetIsolationCut()->GetParticleTypeInCone()!=AliIsolationCut::kOnlyNeutral)
    {
      fhConeSumPtTrack  = new TH2F("hConePtSumTrack",
                                   Form("Track #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                   nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrack->SetYTitle("#Sigma #it{p}_{T}");
      fhConeSumPtTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtTrack) ;

      fhConePtLeadTrack  = new TH2F("hConeLeadPtTrack",
                                   Form("Track leading in isolation cone for #it{R} =  %2.2f",r),
                                   nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConePtLeadTrack->SetYTitle("#it{p}_{T, leading} (GeV/#it{c})");
      fhConePtLeadTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConePtLeadTrack) ;
      
      fhPtTrackInCone  = new TH2F("hPtTrackInCone",
                                  Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f",r),
                                  nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInCone) ;
      
      
      if(fFillUEBandSubtractHistograms)
      {
        fhConeSumPtEtaBandUETrack  = new TH2F("hConePtSumEtaBandUETrack",
                                              "#Sigma track #it{p}_{T} in UE Eta Band",
                                              nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtEtaBandUETrack->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaBandUETrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaBandUETrack) ;
        
        fhConeSumPtPhiBandUETrack  = new TH2F("hConePtSumPhiBandUETrack",
                                              "#Sigma track #it{p}_{T} in UE Phi Band",
                                              nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax*8);
        fhConeSumPtPhiBandUETrack->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiBandUETrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiBandUETrack) ;
        
        
        fhConeSumPtEtaBandUETrackTrigEtaPhi  = new TH2F("hConePtSumEtaBandUETrackTrigEtaPhi",
                                                        "Trigger #eta vs #phi, #Sigma track #it{p}_{T} in UE Eta Band",
                                                        netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtEtaBandUETrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaBandUETrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtEtaBandUETrackTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtEtaBandUETrackTrigEtaPhi) ;
        
        fhConeSumPtPhiBandUETrackTrigEtaPhi  = new TH2F("hConePtSumPhiBandUETrackTrigEtaPhi",
                                                        "Trigger #eta vs #phi, #Sigma track #it{p}_{T} in UE Phi Band",
                                                        netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtPhiBandUETrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiBandUETrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtPhiBandUETrackTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtPhiBandUETrackTrigEtaPhi) ;
        
        fhEtaBandTrack  = new TH2F("hEtaBandTrack",
                                   Form("#eta vs #phi of tracks in #eta band isolation cone for #it{R} =  %2.2f",r),
                                   netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaBandTrack->SetXTitle("#eta");
        fhEtaBandTrack->SetYTitle("#phi");
        outputContainer->Add(fhEtaBandTrack) ;
        
        fhPhiBandTrack  = new TH2F("hPhiBandTrack",
                                   Form("#eta vs #phi of tracks in #phi band isolation cone for #it{R} =  %2.2f",r),
                                   netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhPhiBandTrack->SetXTitle("#eta");
        fhPhiBandTrack->SetYTitle("#phi");
        outputContainer->Add(fhPhiBandTrack) ;
        
        fhConeSumPtEtaUESubTrack  = new TH2F("hConeSumPtEtaUESubTrack",
                                             Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                             nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtEtaUESubTrack->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaUESubTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaUESubTrack) ;
        
        fhConeSumPtPhiUESubTrack  = new TH2F("hConeSumPtPhiUESubTrack",
                                             Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                             nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtPhiUESubTrack->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiUESubTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiUESubTrack) ;
        
        fhConeSumPtEtaUESubTrackTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubTrackTrigEtaPhi",
                                                       Form("Trigger #eta vs #phi, Tracks #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                       netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtEtaUESubTrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaUESubTrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtEtaUESubTrackTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtEtaUESubTrackTrigEtaPhi) ;
        
        fhConeSumPtPhiUESubTrackTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubTrackTrigEtaPhi",
                                                       Form("Trigger #eta vs #phi, Tracks #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                                       netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtPhiUESubTrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiUESubTrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtPhiUESubTrackTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtPhiUESubTrackTrigEtaPhi) ;
        
        fhFractionTrackOutConeEta  = new TH2F("hFractionTrackOutConeEta",
                                              Form("Fraction of the isolation cone #it{R} =  %2.2f, out of tracks #eta acceptance",r),
                                              nptbins,ptmin,ptmax,100,0,1);
        fhFractionTrackOutConeEta->SetYTitle("#it{fraction}");
        fhFractionTrackOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
        outputContainer->Add(fhFractionTrackOutConeEta) ;
        
        fhFractionTrackOutConeEtaTrigEtaPhi  = new TH2F("hFractionTrackOutConeEtaTrigEtaPhi",
                                                        Form("Fraction of the isolation cone #it{R} =  %2.2f, out of tracks #eta acceptance, in trigger #eta-#phi ",r),
                                                        netabins,etamin,etamax,nphibins,phimin,phimax);
        fhFractionTrackOutConeEtaTrigEtaPhi->SetZTitle("#it{fraction}");
        fhFractionTrackOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhFractionTrackOutConeEtaTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhFractionTrackOutConeEtaTrigEtaPhi) ;
        
        fhConeSumPtSubvsConeSumPtTotPhiTrack = new TH2F("hConeSumPtSubvsConeSumPtTotPhiTrack",
                                                        Form("#Sigma #it{p}_{T} in cone after bkg sub from #phi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                        nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubvsConeSumPtTotPhiTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubvsConeSumPtTotPhiTrack->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubvsConeSumPtTotPhiTrack);
        
        fhConeSumPtSubNormvsConeSumPtTotPhiTrack = new TH2F("hConeSumPtSubNormvsConeSumPtTotPhiTrack",
                                                            Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #phi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                            nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubNormvsConeSumPtTotPhiTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubNormvsConeSumPtTotPhiTrack->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotPhiTrack);
        
        fhConeSumPtSubvsConeSumPtTotEtaTrack = new TH2F("hConeSumPtSubvsConeSumPtTotEtaTrack",
                                                        Form("#Sigma #it{p}_{T} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                        nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubvsConeSumPtTotEtaTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubvsConeSumPtTotEtaTrack->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubvsConeSumPtTotEtaTrack);
        
        fhConeSumPtSubNormvsConeSumPtTotEtaTrack = new TH2F("hConeSumPtSubNormvsConeSumPtTotEtaTrack",
                                                            Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                            nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubNormvsConeSumPtTotEtaTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubNormvsConeSumPtTotEtaTrack->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotEtaTrack);
        
        
        // UE in perpendicular cone
        fhPerpConeSumPt  = new TH2F("hPerpConePtSum",
                                    Form("#Sigma #it{p}_{T} in isolation cone at #pm 45 degree phi from trigger particle, #it{R} =  %2.2f",r),
                                    nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhPerpConeSumPt->SetYTitle("#Sigma #it{p}_{T}");
        fhPerpConeSumPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPerpConeSumPt) ;
        
        fhPtInPerpCone  = new TH2F("hPtInPerpCone",
                                   Form("#it{p}_{T} in isolation cone at #pm 45 degree phi from trigger particle, #it{R} =  %2.2f",r),
                                   nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtInPerpCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtInPerpCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtInPerpCone) ;
        
        fhEtaPhiTrack= new TH2F("hEtaPhiTrack",
                                Form("#eta vs #phi of all Tracks"),
                                netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaPhiTrack->SetXTitle("#eta");
        fhEtaPhiTrack->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiTrack) ;
        
        fhEtaPhiInConeTrack= new TH2F("hEtaPhiInConeTrack",
                                      Form("#eta vs #phi of Tracks in cone for #it{R} =  %2.2f",r),
                                      netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaPhiInConeTrack->SetXTitle("#eta");
        fhEtaPhiInConeTrack->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiInConeTrack) ;
        
        fhConeSumPtVSUETracksEtaBand  = new TH2F("hConeSumPtVSUETracksEtaBand",
                                                 Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in eta band for tracks (before normalization), R=%2.2f",r),
                                                 nptsumbins,ptsummin,ptsummax,2*nptsumbins,ptsummin,2*ptsummax);
        fhConeSumPtVSUETracksEtaBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
        fhConeSumPtVSUETracksEtaBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtVSUETracksEtaBand);
        
        fhConeSumPtVSUETracksPhiBand  = new TH2F("hConeSumPtVSUETracksPhiBand",
                                                 Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in phi band for tracks (before normalization), R=%2.2f",r),
                                                 nptsumbins,ptsummin,ptsummax,8*nptsumbins,ptsummin,8*ptsummax);
        fhConeSumPtVSUETracksPhiBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
        fhConeSumPtVSUETracksPhiBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtVSUETracksPhiBand);
      }
    }
    
    if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
    {
      fhConeSumPtClustervsTrack   = new TH2F("hConePtSumClustervsTrack",
                                             Form("Track vs Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                             nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
      fhConeSumPtClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
      outputContainer->Add(fhConeSumPtClustervsTrack) ;
      
      if(fFillCellHistograms)
      {
        fhConeSumPtCellvsTrack   = new TH2F("hConePtSumCellvsTrack",
                                            Form("Track vs cell #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                            nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
        fhConeSumPtCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhConeSumPtCellvsTrack) ;
        
        fhConeSumPtCellTrack = new TH2F("hConePtSumCellTrack",
                                        Form("Track and Cell #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                        nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtCellTrack->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtCellTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtCellTrack) ;
        
        fhConeSumPtCellTrackTrigEtaPhi  = new TH2F("hConePtSumCellTrackTrigEtaPhi",
                                                   Form("Trigger #eta vs #phi, #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                                   netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtCellTrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtCellTrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtCellTrackTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtCellTrackTrigEtaPhi) ;
        
      }
      
      if(fFillUEBandSubtractHistograms)
      {
        fhConeSumPtEtaUESub  = new TH2F("hConeSumPtEtaUESub",
                                        Form("#Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                        nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtEtaUESub->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaUESub->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaUESub) ;
        
        fhConeSumPtPhiUESub  = new TH2F("hConeSumPtPhiUESub",
                                        Form("#Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                        nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtPhiUESub->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiUESub->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiUESub) ;
        
        fhConeSumPtEtaUESubTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubTrigEtaPhi",
                                                  Form("Trigger #eta vs #phi, #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                  netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtEtaUESubTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaUESubTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtEtaUESubTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtEtaUESubTrigEtaPhi) ;
        
        fhConeSumPtPhiUESubTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubTrigEtaPhi",
                                                  Form("Trigger #eta vs #phi, #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                                  netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtPhiUESubTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiUESubTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtPhiUESubTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtPhiUESubTrigEtaPhi) ;
        
        fhConeSumPtEtaUESubClustervsTrack   = new TH2F("hConePtSumEtaUESubClustervsTrack",
                                                       Form("Track vs Cluster #Sigma #it{p}_{T} UE sub eta band in isolation cone for #it{R} =  %2.2f",r),
                                                       2*nptsumbins,-ptsummax,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtEtaUESubClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhConeSumPtEtaUESubClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhConeSumPtEtaUESubClustervsTrack) ;
        
        fhConeSumPtPhiUESubClustervsTrack   = new TH2F("hConePhiUESubPtSumClustervsTrack",
                                                       Form("Track vs Cluster #Sigma #it{p}_{T} UE sub phi band in isolation cone for #it{R} =  %2.2f",r),
                                                       2*nptsumbins,-ptsummax,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtPhiUESubClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhConeSumPtPhiUESubClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhConeSumPtPhiUESubClustervsTrack) ;
        
        fhEtaBandClustervsTrack   = new TH2F("hEtaBandClustervsTrack",
                                             Form("Track vs Cluster #Sigma #it{p}_{T} in Eta band in isolation cone for #it{R} =  %2.2f",r),
                                             nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhEtaBandClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhEtaBandClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhEtaBandClustervsTrack) ;
        
        fhPhiBandClustervsTrack   = new TH2F("hPhiBandClustervsTrack",
                                             Form("Track vs Cluster #Sigma #it{p}_{T} in Phi band in isolation cone for #it{R} =  %2.2f",r),
                                             nptsumbins,ptsummin,ptsummax*4,nptsumbins,ptsummin,ptsummax*8);
        fhPhiBandClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhPhiBandClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhPhiBandClustervsTrack) ;
        
        fhEtaBandNormClustervsTrack   = new TH2F("hEtaBandNormClustervsTrack",
                                                 Form("Track vs Cluster #Sigma #it{p}_{T} in Eta band in isolation cone for #it{R} =  %2.2f",r),
                                                 nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhEtaBandNormClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhEtaBandNormClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhEtaBandNormClustervsTrack) ;
        
        fhPhiBandNormClustervsTrack   = new TH2F("hPhiBandNormClustervsTrack",
                                                 Form("Track vs Cluster #Sigma #it{p}_{T} in Phi band in isolation cone for #it{R} =  %2.2f",r),
                                                 nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhPhiBandNormClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhPhiBandNormClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhPhiBandNormClustervsTrack) ;
        
        fhConeSumPtEtaUESubClustervsTrack   = new TH2F("hConePtSumEtaUESubClustervsTrack",
                                                       Form("Track vs Cluster #Sigma #it{p}_{T} UE sub eta band in isolation cone for #it{R} =  %2.2f",r),
                                                       2*nptsumbins,-ptsummax,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtEtaUESubClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhConeSumPtEtaUESubClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhConeSumPtEtaUESubClustervsTrack) ;
        
        fhConeSumPtPhiUESubClustervsTrack   = new TH2F("hConePhiUESubPtSumClustervsTrack",
                                                       Form("Track vs Cluster #Sigma #it{p}_{T} UE sub phi band in isolation cone for #it{R} =  %2.2f",r),
                                                       2*nptsumbins,-ptsummax,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtPhiUESubClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhConeSumPtPhiUESubClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhConeSumPtPhiUESubClustervsTrack) ;
        
        if(fFillCellHistograms)
        {
          
          fhConeSumPtEtaUESubCellvsTrack   = new TH2F("hConePtSumEtaUESubCellvsTrack",
                                                      Form("Track vs Cell #Sigma #it{p}_{T} UE sub eta band in isolation cone for #it{R} =  %2.2f",r),
                                                      2*nptsumbins,-ptsummax,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtEtaUESubCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
          fhConeSumPtEtaUESubCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
          outputContainer->Add(fhConeSumPtEtaUESubCellvsTrack) ;
          
          fhConeSumPtPhiUESubCellvsTrack   = new TH2F("hConePhiUESubPtSumCellvsTrack",
                                                      Form("Track vs Cell #Sigma #it{p}_{T} UE sub phi band in isolation cone for #it{R} =  %2.2f",r),
                                                      2*nptsumbins,-ptsummax,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtPhiUESubCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
          fhConeSumPtPhiUESubCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
          outputContainer->Add(fhConeSumPtPhiUESubCellvsTrack) ;
          
          fhEtaBandCellvsTrack   = new TH2F("hEtaBandCellvsTrack",
                                            Form("Track vs Cell #Sigma #it{p}_{T} in Eta band in isolation cone for #it{R} =  %2.2f",r),
                                            nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
          fhEtaBandCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
          fhEtaBandCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
          outputContainer->Add(fhEtaBandCellvsTrack) ;
          
          fhPhiBandCellvsTrack   = new TH2F("hPhiBandCellvsTrack",
                                            Form("Track vs Cell #Sigma #it{p}_{T} in Phi band in isolation cone for #it{R} =  %2.2f",r),
                                            nptsumbins,ptsummin,ptsummax*4,nptsumbins,ptsummin,ptsummax*8);
          fhPhiBandCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
          fhPhiBandCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
          outputContainer->Add(fhPhiBandCellvsTrack) ;
          
          fhEtaBandNormCellvsTrack   = new TH2F("hEtaBandNormCellvsTrack",
                                                Form("Track vs Cell #Sigma #it{p}_{T} in Eta band in isolation cone for #it{R} =  %2.2f",r),
                                                nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
          fhEtaBandNormCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
          fhEtaBandNormCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
          outputContainer->Add(fhEtaBandNormCellvsTrack) ;
          
          fhPhiBandNormCellvsTrack   = new TH2F("hPhiBandNormCellvsTrack",
                                                Form("Track vs Cell #Sigma #it{p}_{T} in Phi band in isolation cone for #it{R} =  %2.2f",r),
                                                nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
          fhPhiBandNormCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
          fhPhiBandNormCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
          outputContainer->Add(fhPhiBandNormCellvsTrack) ;
          
          fhConeSumPtEtaUESubTrackCell  = new TH2F("hConeSumPtEtaUESubTrackCell",
                                                   Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                   nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtEtaUESubTrackCell->SetYTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaUESubTrackCell->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtEtaUESubTrackCell) ;
          
          fhConeSumPtPhiUESubTrackCell  = new TH2F("hConeSumPtPhiUESubTrackCell",
                                                   Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                                   nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtPhiUESubTrackCell->SetYTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiUESubTrackCell->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPhiUESubTrackCell) ;
          
          fhConeSumPtEtaUESubTrackCellTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubTrackCellTrigEtaPhi",
                                                             Form("Trigger #eta vs #phi, Tracks #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                             netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaUESubTrackCellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaUESubTrackCellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaUESubTrackCellTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaUESubTrackCellTrigEtaPhi) ;
          
          fhConeSumPtPhiUESubTrackCellTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubTrackCellTrigEtaPhi",
                                                             Form("Trigger #eta vs #phi, Tracks #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                                             netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiUESubTrackCellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiUESubTrackCellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiUESubTrackCellTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtPhiUESubTrackCellTrigEtaPhi) ;
        }
      }
    }
    
    for(Int_t iso = 0; iso < 2; iso++)
    {
      if(fFillTMHisto)
      {
        fhTrackMatchedDEta[iso]  = new TH2F
        (Form("hTrackMatchedDEta%s",isoName[iso].Data()),
         Form("%s - d#eta of cluster-track vs cluster energy, %s",isoTitle[iso].Data(),parTitle.Data()),
         nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
        fhTrackMatchedDEta[iso]->SetYTitle("d#eta");
        fhTrackMatchedDEta[iso]->SetXTitle("E_{cluster} (GeV)");
        
        fhTrackMatchedDPhi[iso]  = new TH2F
        (Form("hTrackMatchedDPhi%s",isoName[iso].Data()),
         Form("%s - d#phi of cluster-track vs cluster energy, %s",isoTitle[iso].Data(),parTitle.Data()),
         nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
        fhTrackMatchedDPhi[iso]->SetYTitle("d#phi (rad)");
        fhTrackMatchedDPhi[iso]->SetXTitle("E_{cluster} (GeV)");
        
        fhTrackMatchedDEtaDPhi[iso]  = new TH2F
        (Form("hTrackMatchedDEtaDPhi%s",isoName[iso].Data()),
         Form("%s - d#eta vs d#phi of cluster-track, %s",isoTitle[iso].Data(),parTitle.Data()),
         nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
        fhTrackMatchedDEtaDPhi[iso]->SetYTitle("d#phi (rad)");
        fhTrackMatchedDEtaDPhi[iso]->SetXTitle("d#eta");
        
        outputContainer->Add(fhTrackMatchedDEta[iso]) ;
        outputContainer->Add(fhTrackMatchedDPhi[iso]) ;
        outputContainer->Add(fhTrackMatchedDEtaDPhi[iso]) ;
        
        fhdEdx[iso]  = new TH2F
        (Form("hdEdx%s",isoName[iso].Data()),
         Form("%s - Matched track <d#it{E}/d#it{x}> vs cluster #it{E}, %s",isoTitle[iso].Data(),parTitle.Data()),
         nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
        fhdEdx[iso]->SetXTitle("#it{E} (GeV)");
        fhdEdx[iso]->SetYTitle("<d#it{E}/d#it{x}>");
        outputContainer->Add(fhdEdx[iso]);
        
        fhEOverP[iso]  = new TH2F
        (Form("hEOverP%s",isoName[iso].Data()),
         Form("%s - Matched track #it{E}/#it{p} vs cluster, %s",isoTitle[iso].Data(),parTitle.Data()),
         nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
        fhEOverP[iso]->SetXTitle("#it{E} (GeV)");
        fhEOverP[iso]->SetYTitle("#it{E}/#it{p}");
        outputContainer->Add(fhEOverP[iso]);
        
        if(IsDataMC())
        {
          fhTrackMatchedMCParticle[iso]  = new TH2F
          (Form("hTrackMatchedMCParticle%s",isoName[iso].Data()),
           Form("%s - Origin of particle vs cluster #it{E}, %s",isoTitle[iso].Data(),parTitle.Data()),
           nptbins,ptmin,ptmax,8,0,8);
          fhTrackMatchedMCParticle[iso]->SetXTitle("#it{E} (GeV)");
          //fhTrackMatchedMCParticle[iso]->SetYTitle("Particle type");
          
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(1 ,"Photon");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(2 ,"Electron");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(3 ,"Meson Merged");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(4 ,"Rest");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(5 ,"Conv. Photon");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(6 ,"Conv. Electron");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(7 ,"Conv. Merged");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(8 ,"Conv. Rest");
          
          outputContainer->Add(fhTrackMatchedMCParticle[iso]);
        }
      }
      
      if(fFillSSHisto)
      {
        fhELambda0[iso]  = new TH2F
        (Form("hELambda0%s",isoName[iso].Data()),
         Form("%s cluster : #it{E} vs #lambda_{0}, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhELambda0[iso]->SetYTitle("#lambda_{0}^{2}");
        fhELambda0[iso]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhELambda0[iso]) ;
        
        fhELambda1[iso]  = new TH2F
        (Form("hELambda1%s",isoName[iso].Data()),
         Form("%s cluster: #it{E} vs #lambda_{1}, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhELambda1[iso]->SetYTitle("#lambda_{1}^{2}");
        fhELambda1[iso]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhELambda1[iso]) ;
        
        fhPtLambda0[iso]  = new TH2F
        (Form("hPtLambda0%s",isoName[iso].Data()),
         Form("%s cluster : #it{p}_{T} vs #lambda_{0}, %s",isoTitle[iso].Data(), parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhPtLambda0[iso]->SetYTitle("#lambda_{0}^{2}");
        fhPtLambda0[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtLambda0[iso]) ;
        
        if(IsDataMC())
        {
          for(Int_t imc = 0; imc < 9; imc++)
          {
            fhPtLambda0MC[imc][iso]  = new TH2F(Form("hPtLambda0%s_MC%s",isoName[iso].Data(),mcPartName[imc].Data()),
                                                Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s",isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhPtLambda0MC[imc][iso]->SetYTitle("#lambda_{0}^{2}");
            fhPtLambda0MC[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add( fhPtLambda0MC[imc][iso]) ;
          }
        }
        
        if(fIsoDetector=="EMCAL" &&  GetFirstSMCoveredByTRD() >= 0)
        {
          fhPtLambda0TRD[iso]  = new TH2F
          (Form("hPtLambda0TRD%s",isoName[iso].Data()),
           Form("%s cluster: #it{p}_{T} vs #lambda_{0}, SM behind TRD, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtLambda0TRD[iso]->SetYTitle("#lambda_{0}^{2}");
          fhPtLambda0TRD[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtLambda0TRD[iso]) ;
          
          fhELambda0TRD[iso]  = new TH2F
          (Form("hELambda0TRD%s",isoName[iso].Data()),
           Form("%s cluster: #it{E} vs #lambda_{0}, SM behind TRD, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda0TRD[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0TRD[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda0TRD[iso]) ;
          
          fhELambda1TRD[iso]  = new TH2F
          (Form("hELambda1TRD%s",isoName[iso].Data()),
           Form("%s cluster: #it{E} vs #lambda_{1}, SM behind TRD, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda1TRD[iso]->SetYTitle("#lambda_{1}^{2}");
          fhELambda1TRD[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda1TRD[iso]) ;
        }
        
        if(fFillNLMHistograms)
        {
          fhNLocMax[iso] = new TH2F
          (Form("hNLocMax%s",isoName[iso].Data()),
           Form("%s - Number of local maxima in cluster, %s",isoTitle[iso].Data(),parTitle.Data()),
           nptbins,ptmin,ptmax,10,0,10);
          fhNLocMax[iso]->SetYTitle("#it{NLM}");
          fhNLocMax[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhNLocMax[iso]) ;
          
          fhELambda0LocMax1[iso]  = new TH2F
          (Form("hELambda0LocMax1%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{0}, #it{NLM}=1, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda0LocMax1[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0LocMax1[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda0LocMax1[iso]) ;
          
          fhELambda1LocMax1[iso]  = new TH2F
          (Form("hELambda1LocMax1%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{1}, #it{NLM}=1, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda1LocMax1[iso]->SetYTitle("#lambda_{1}^{2}");
          fhELambda1LocMax1[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda1LocMax1[iso]) ;
          
          fhELambda0LocMax2[iso]  = new TH2F
          (Form("hELambda0LocMax2%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{0}, #it{NLM}=2, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda0LocMax2[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0LocMax2[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda0LocMax2[iso]) ;
          
          fhELambda1LocMax2[iso]  = new TH2F
          (Form("hELambda1LocMax2%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{1}, #it{NLM}=2, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda1LocMax2[iso]->SetYTitle("#lambda_{1}^{2}");
          fhELambda1LocMax2[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda1LocMax2[iso]) ;
          
          fhELambda0LocMaxN[iso]  = new TH2F
          ( Form("hELambda0LocMaxN%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{0}, #it{NLM}>2, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda0LocMaxN[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0LocMaxN[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda0LocMaxN[iso]) ;
          
          fhELambda1LocMaxN[iso]  = new TH2F
          (Form("hELambda1LocMaxN%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{1}, #it{NLM}>2, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda1LocMaxN[iso]->SetYTitle("#lambda_{1}^{2}");
          fhELambda1LocMaxN[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda1LocMaxN[iso]) ;
        } // NLM
      } // SS histo
    } // control histograms for isolated and non isolated objects
    
    
    if(fFillPileUpHistograms)
    {
      fhPtTrackInConeOtherBC  = new TH2F("hPtTrackInConeOtherBC",
                                         Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC!=0",r),
                                         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeOtherBC->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeOtherBC->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeOtherBC) ;
      
      fhPtTrackInConeOtherBCPileUpSPD  = new TH2F("hPtTrackInConeOtherBCPileUpSPD",
                                                  Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC!=0, pile-up from SPD",r),
                                                  nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeOtherBCPileUpSPD->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeOtherBCPileUpSPD->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeOtherBCPileUpSPD) ;
      
      fhPtTrackInConeBC0  = new TH2F("hPtTrackInConeBC0",
                                     Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC==0",r),
                                     nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeBC0->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeBC0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeBC0) ;
      
      fhPtTrackInConeVtxBC0  = new TH2F("hPtTrackInConeVtxBC0",
                                        Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC==0",r),
                                        nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeVtxBC0->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeVtxBC0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeVtxBC0) ;
      
      
      fhPtTrackInConeBC0PileUpSPD  = new TH2F("hPtTrackInConeBC0PileUpSPD",
                                              Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC==0, pile-up from SPD",r),
                                              nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeBC0PileUpSPD->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeBC0PileUpSPD->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeBC0PileUpSPD) ;
      
      
      for (Int_t i = 0; i < 7 ; i++)
      {
        fhPtInConePileUp[i]  = new TH2F(Form("hPtInConePileUp%s",pileUpName[i].Data()),
                                        Form("#it{p}_{T} in isolation cone for #it{R} =  %2.2f, from pile-up (%s)",r,pileUpName[i].Data()),
                                        nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtInConePileUp[i]->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtInConePileUp[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtInConePileUp[i]) ;
      }
    }
    
    if(IsDataMC())
    {
      // For histograms in arrays, index in the array, corresponding to any particle origin
      
      for(Int_t i = 0; i < 6; i++)
      {
        fhEPrimMC[i]  = new TH1F(Form("hEPrim_MC%s",ppname[i].Data()),
                                 Form("primary photon  %s : #it{E}, %s",pptype[i].Data(),parTitle.Data()),
                                 nptbins,ptmin,ptmax);
        fhEPrimMC[i]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEPrimMC[i]) ;
        
        fhPtPrimMCiso[i]  = new TH1F(Form("hPtPrim_MCiso%s",ppname[i].Data()),
                                     Form("primary isolated photon %s : #it{p}_{T}, %s",pptype[i].Data(),parTitle.Data()),
                                     nptbins,ptmin,ptmax);
        fhPtPrimMCiso[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCiso[i]) ;
        
        fhEtaPrimMC[i]  = new TH2F(Form("hEtaPrim_MC%s",ppname[i].Data()),
                                   Form("primary photon %s : #eta vs #it{p}_{T}, %s",pptype[i].Data(),parTitle.Data()),
                                   nptbins,ptmin,ptmax,200,-2,2);
        fhEtaPrimMC[i]->SetYTitle("#eta");
        fhEtaPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhEtaPrimMC[i]) ;
        
        fhPhiPrimMC[i]  = new TH2F(Form("hPhiPrim_MC%s",ppname[i].Data()),
                                   Form("primary photon %s : #phi vs #it{p}_{T}, %s",pptype[i].Data(),parTitle.Data()),
                                   nptbins,ptmin,ptmax,200,0.,TMath::TwoPi());
        fhPhiPrimMC[i]->SetYTitle("#phi");
        fhPhiPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPhiPrimMC[i]) ;
      }
      
    }//Histos with MC
    
  }
  
  if(fMakeSeveralIC)
  {
    const Int_t buffersize = 255;
    char name[buffersize];
    char title[buffersize];
    for(Int_t icone = 0; icone<fNCones; icone++)
    {
      // sum pt in cone vs. pt leading
      snprintf(name, buffersize,"hSumPtLeadingPt_Cone_%d",icone);
      snprintf(title, buffersize,"#Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",fConeSizes[icone]);
      fhSumPtLeadingPt[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhSumPtLeadingPt[icone] ->SetYTitle("#sum_{cone}#it{p}_{T} (GeV/#it{c})");//#Sigma #it{p}_{T}
      fhSumPtLeadingPt[icone] ->SetXTitle("#it{p}_{T}^{leading} (GeV/#it{c})");
      outputContainer->Add(fhSumPtLeadingPt[icone]) ;
      
      // pt in cone vs. pt leading
      snprintf(name, buffersize,"hPtLeadingPt_Cone_%d",icone);
      snprintf(title, buffersize,"#it{p}_{T} in isolation cone for #it{R} =  %2.2f",fConeSizes[icone]);
      fhPtLeadingPt[icone]  = new TH2F(name, title,  nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtLeadingPt[icone] ->SetYTitle("#it{p}_{T}^{cone} (GeV/#it{c})");
      fhPtLeadingPt[icone] ->SetXTitle("#it{p}_{T}^{leading} (GeV/#it{c})");
      outputContainer->Add(fhPtLeadingPt[icone]) ;
      
      // sum pt in cone vs. pt leading in the forward region (for background subtraction studies)
      snprintf(name, buffersize,"hPerpSumPtLeadingPt_Cone_%d",icone);
      snprintf(title, buffersize,"#Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",fConeSizes[icone]);
      fhPerpSumPtLeadingPt[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhPerpSumPtLeadingPt[icone] ->SetYTitle("#sum_{cone}#it{p}_{T} (GeV/#it{c})");//#Sigma #it{p}_{T}
      fhPerpSumPtLeadingPt[icone] ->SetXTitle("#it{p}_{T}^{leading} (GeV/#it{c})");
      outputContainer->Add(fhPerpSumPtLeadingPt[icone]) ;
      
      // pt in cone vs. pt leading in the forward region (for background subtraction studies)
      snprintf(name, buffersize,"hPerpPtLeadingPt_Cone_%d",icone);
      snprintf(title, buffersize,"#it{p}_{T} in isolation cone for #it{R} =  %2.2f",fConeSizes[icone]);
      fhPerpPtLeadingPt[icone]  = new TH2F(name, title,  nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPerpPtLeadingPt[icone] ->SetYTitle("#it{p}_{T}^{cone} (GeV/#it{c})");
      fhPerpPtLeadingPt[icone] ->SetXTitle("#it{p}_{T}^{leading} (GeV/#it{c})");
      outputContainer->Add(fhPerpPtLeadingPt[icone]) ;
      
      if(IsDataMC())
      {
        for(Int_t imc = 0; imc < 9; imc++)
        {
          snprintf(name , buffersize,"hSumPtLeadingPt_MC%s_Cone_%d",mcPartName[imc].Data(),icone);
          snprintf(title, buffersize,"Candidate %s #it{p}_{T} vs cone #Sigma #it{p}_{T} for #it{R}=%2.2f",mcPartType[imc].Data(),fConeSizes[icone]);
          fhSumPtLeadingPtMC[imc][icone] = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
          fhSumPtLeadingPtMC[imc][icone]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhSumPtLeadingPtMC[imc][icone]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
          outputContainer->Add(fhSumPtLeadingPtMC[imc][icone]) ;
        }
      }//Histos with MC
      
      for(Int_t ipt = 0; ipt<fNPtThresFrac;ipt++)
      {
        snprintf(name, buffersize,"hPtThres_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #it{p}_{T} distribution for #it{R} =  %2.2f and #it{p}_{T}^{th} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtThresholds[ipt]);
        fhPtThresIsolated[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
        fhPtThresIsolated[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtThresIsolated[icone][ipt]) ;
        
        snprintf(name, buffersize,"hPtFrac_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #it{p}_{T} distribution for #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
        fhPtFracIsolated[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
        fhPtFracIsolated[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtFracIsolated[icone][ipt]) ;
        
        snprintf(name, buffersize,"hSumPt_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #it{p}_{T} distribution for #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhSumPtIsolated[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
        // fhSumPtIsolated[icone][ipt]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhSumPtIsolated[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhSumPtIsolated[icone][ipt]) ;
        
        snprintf(name, buffersize,"hPtSumDensity_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #it{p}_{T} distribution for density in #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhPtSumDensityIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
        //fhPtSumIsolated[icone][ipt]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhPtSumDensityIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtSumDensityIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hPtFracPtSum_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #it{p}_{T} distribution for PtFracPtSum in #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
        fhPtFracPtSumIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
        //fhPtSumIsolated[icone][ipt]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhPtFracPtSumIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtFracPtSumIso[icone][ipt]) ;
        
        // eta:phi
        snprintf(name, buffersize,"hEtaPhiPtThres_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for #it{R} =  %2.2f and #it{p}_{T}^{th} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtThresholds[ipt]);
        fhEtaPhiPtThresIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtThresIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtThresIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiPtThresIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiPtFrac_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
        fhEtaPhiPtFracIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtFracIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtFracIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiPtFracIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiPtSum_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhEtaPhiPtSumIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtSumIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtSumIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiPtSumIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiSumDensity_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for density #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhEtaPhiSumDensityIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiSumDensityIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiSumDensityIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiSumDensityIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiFracPtSum_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for FracPtSum #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
        fhEtaPhiFracPtSumIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiFracPtSumIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiFracPtSumIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiFracPtSumIso[icone][ipt]) ;
        
        if(fFillTaggedDecayHistograms)
        {
          // pt decays isolated
          snprintf(name, buffersize,"hPtThres_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #it{p}_{T} distribution for #it{R} =  %2.2f and #it{p}_{T}^{th} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtThresholds[ipt]);
          fhPtPtThresDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtPtThresDecayIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtPtThresDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hPtFrac_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #it{p}_{T} distribution for #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
          fhPtPtFracDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtPtFracDecayIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtPtFracDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #it{p}_{T} distribution for #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
          fhPtPtSumDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
          //  fhPtPtSumDecayIso[icone]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhPtPtSumDecayIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtPtSumDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hPtSumDensity_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #it{p}_{T} distribution for density in #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
          fhPtSumDensityDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
          //  fhPtPtSumDecayIso[icone]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhPtSumDensityDecayIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtSumDensityDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hPtFracPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #it{p}_{T} distribution for PtFracPtSum in #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
          fhPtFracPtSumDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
          //  fhPtPtSumDecayIso[icone]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhPtFracPtSumDecayIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtFracPtSumDecayIso[icone][ipt]) ;
          
          // eta:phi decays
          snprintf(name, buffersize,"hEtaPhiPtThres_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for #it{R} =  %2.2f and #it{p}_{T}^{th} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtThresholds[ipt]);
          fhEtaPhiPtThresDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiPtThresDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiPtThresDecayIso[icone][ipt]->SetYTitle("#phi");
          outputContainer->Add(fhEtaPhiPtThresDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hEtaPhiPtFrac_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
          fhEtaPhiPtFracDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiPtFracDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiPtFracDecayIso[icone][ipt]->SetYTitle("#phi");
          outputContainer->Add(fhEtaPhiPtFracDecayIso[icone][ipt]) ;
          
          
          snprintf(name, buffersize,"hEtaPhiPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
          fhEtaPhiPtSumDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiPtSumDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiPtSumDecayIso[icone][ipt]->SetYTitle("#phi");
          outputContainer->Add(fhEtaPhiPtSumDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hEtaPhiSumDensity_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for density #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
          fhEtaPhiSumDensityDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiSumDensityDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiSumDensityDecayIso[icone][ipt]->SetYTitle("#phi");
          outputContainer->Add(fhEtaPhiSumDensityDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hEtaPhiFracPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for FracPtSum #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
          fhEtaPhiFracPtSumDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiFracPtSumDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiFracPtSumDecayIso[icone][ipt]->SetYTitle("#phi");
          outputContainer->Add(fhEtaPhiFracPtSumDecayIso[icone][ipt]) ;
          
        }
        
        if(IsDataMC())
        {
          for(Int_t imc = 0; imc < 9; imc++)
          {
            snprintf(name , buffersize,"hPtThreshMC%s_Cone_%d_Pt%d",mcPartName[imc].Data(),icone,ipt);
            snprintf(title, buffersize,"Isolated %s #it{p}_{T} for #it{R}=%2.2f and #it{p}_{T}^{th}=%2.2f",
                     mcPartType[imc].Data(),fConeSizes[icone], fPtThresholds[ipt]);
            fhPtThresIsolatedMC[imc][icone][ipt] = new TH1F(name, title,nptbins,ptmin,ptmax);
            fhPtThresIsolatedMC[imc][icone][ipt]->SetYTitle("#it{counts}");
            fhPtThresIsolatedMC[imc][icone][ipt]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add(fhPtThresIsolatedMC[imc][icone][ipt]) ;
            
            
            snprintf(name , buffersize,"hPtFracMC%s_Cone_%d_Pt%d",mcPartName[imc].Data(),icone,ipt);
            snprintf(title, buffersize,"Isolated %s #it{p}_{T} for #it{R}=%2.2f and #Sigma #it{p}_{T}^{in cone}/#it{p}_{T}^{trig}=%2.2f",
                     mcPartType[imc].Data(),fConeSizes[icone], fPtFractions[ipt]);
            fhPtFracIsolatedMC[imc][icone][ipt] = new TH1F(name, title,nptbins,ptmin,ptmax);
            fhPtFracIsolatedMC[imc][icone][ipt]->SetYTitle("#it{counts}");
            fhPtFracIsolatedMC[imc][icone][ipt]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add(fhPtFracIsolatedMC[imc][icone][ipt]) ;
            
            snprintf(name , buffersize,"hSumPtMC%s_Cone_%d_Pt%d",mcPartName[imc].Data(),icone,ipt);
            snprintf(title, buffersize,"Isolated %s #it{p}_{T} for #it{R}=%2.2f and #Sigma #it{p}_{T}^{in cone}=%2.2f",
                     mcPartType[imc].Data(),fConeSizes[icone], fSumPtThresholds[ipt]);
            fhSumPtIsolatedMC[imc][icone][ipt] = new TH1F(name, title,nptbins,ptmin,ptmax);
            fhSumPtIsolatedMC[imc][icone][ipt]->SetYTitle("#it{counts}");
            fhSumPtIsolatedMC[imc][icone][ipt]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add(fhSumPtIsolatedMC[imc][icone][ipt]) ;
          }
        }//Histos with MC
      }//icone loop
    }//ipt loop
  }
  
  if(fFillPileUpHistograms)
  {
    for (Int_t i = 0; i < 7 ; i++)
    {
      fhEIsoPileUp[i]   = new TH1F(Form("hEPileUp%s",pileUpName[i].Data()),
                                   Form("Number of isolated particles vs E, %s, pile-up event by %s",parTitle.Data(),pileUpName[i].Data()),
                                   nptbins,ptmin,ptmax);
      fhEIsoPileUp[i]->SetYTitle("d#it{N} / d#it{E}");
      fhEIsoPileUp[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhEIsoPileUp[i]) ;
      
      fhPtIsoPileUp[i]  = new TH1F(Form("hPtPileUp%s",pileUpName[i].Data()),
                                   Form("Number of isolated particles vs #it{p}_{T}, %s, pile-up event by %s",parTitle.Data(),pileUpName[i].Data()),
                                   nptbins,ptmin,ptmax);
      fhPtIsoPileUp[i]->SetYTitle("d#it{N} / #it{p}_{T}");
      fhPtIsoPileUp[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtIsoPileUp[i]) ;
      
      fhENoIsoPileUp[i]   = new TH1F(Form("hENoIsoPileUp%s",pileUpName[i].Data()),
                                     Form("Number of not isolated particles vs E, %s, pile-up event by %s",parTitle.Data(),pileUpName[i].Data()),
                                     nptbins,ptmin,ptmax);
      fhENoIsoPileUp[i]->SetYTitle("d#it{N} / dE");
      fhENoIsoPileUp[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhENoIsoPileUp[i]) ;
      
      fhPtNoIsoPileUp[i]  = new TH1F(Form("hPtNoIsoPileUp%s",pileUpName[i].Data()),
                                     Form("Number of not isolated particles vs #it{p}_{T}, %s, pile-up event by %s",parTitle.Data(),pileUpName[i].Data()),
                                     nptbins,ptmin,ptmax);
      fhPtNoIsoPileUp[i]->SetYTitle("d#it{N} / #it{p}_{T}");
      fhPtNoIsoPileUp[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtNoIsoPileUp[i]) ;
    }
    
    fhTimeENoCut  = new TH2F ("hTimeE_NoCut","time of cluster vs E of clusters, no cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeENoCut->SetXTitle("#it{E} (GeV)");
    fhTimeENoCut->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeENoCut);
    
    fhTimeESPD  = new TH2F ("hTimeE_SPD","time of cluster vs E of clusters, SPD cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeESPD->SetXTitle("#it{E} (GeV)");
    fhTimeESPD->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeESPD);
    
    fhTimeESPDMulti  = new TH2F ("hTimeE_SPDMulti","time of cluster vs E of clusters, SPD multi cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeESPDMulti->SetXTitle("#it{E} (GeV)");
    fhTimeESPDMulti->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeESPDMulti);
    
    fhTimeNPileUpVertSPD  = new TH2F ("hTime_NPileUpVertSPD","time of cluster vs N pile-up SPD vertex", ntimebins,timemin,timemax,50,0,50);
    fhTimeNPileUpVertSPD->SetYTitle("# vertex ");
    fhTimeNPileUpVertSPD->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeNPileUpVertSPD);
    
    fhTimeNPileUpVertTrack  = new TH2F ("hTime_NPileUpVertTracks","time of cluster vs N pile-up Tracks vertex", ntimebins,timemin,timemax, 50,0,50 );
    fhTimeNPileUpVertTrack->SetYTitle("# vertex ");
    fhTimeNPileUpVertTrack->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeNPileUpVertTrack);
    
    fhTimeNPileUpVertContributors  = new TH2F ("hTime_NPileUpVertContributors","time of cluster vs N constributors to pile-up SPD vertex", ntimebins,timemin,timemax,50,0,50);
    fhTimeNPileUpVertContributors->SetYTitle("# vertex ");
    fhTimeNPileUpVertContributors->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeNPileUpVertContributors);
    
    fhTimePileUpMainVertexZDistance  = new TH2F ("hTime_PileUpMainVertexZDistance","time of cluster vs distance in Z pile-up SPD vertex - main SPD vertex",ntimebins,timemin,timemax,100,0,50);
    fhTimePileUpMainVertexZDistance->SetYTitle("distance #it{z} (cm) ");
    fhTimePileUpMainVertexZDistance->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimePileUpMainVertexZDistance);
    
    fhTimePileUpMainVertexZDiamond  = new TH2F ("hTime_PileUpMainVertexZDiamond","time of cluster vs distance in Z pile-up SPD vertex - z diamond",ntimebins,timemin,timemax,100,0,50);
    fhTimePileUpMainVertexZDiamond->SetYTitle("diamond distance #it{z} (cm) ");
    fhTimePileUpMainVertexZDiamond->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimePileUpMainVertexZDiamond);
  }
  
  return outputContainer ;
  
}

//____________________________________________________
Int_t AliAnaParticleIsolation::GetMCIndex(Int_t mcTag)
{
  // Histogram index depending on origin of candidate
  
  if(!IsDataMC()) return -1;
  
  if     (GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPrompt))
  {
    return kmcPrompt;
  }
  else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCFragmentation))
  {
    return kmcFragment;
  }
  else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0))
  {
    return kmcPi0;
  }
  else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0Decay))
  {
    return kmcPi0Decay;
  }
  else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEtaDecay))
  {
    return kmcEtaDecay;
  }
  else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCOtherDecay))
  {
    return kmcOtherDecay;
  }
  else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron))
  {
    return kmcElectron;
  }
  else // anything else
  {
    // careful can contain also other decays, to be checked.
    return kmcHadron;
  }
}

//__________________________________
void AliAnaParticleIsolation::Init()
{
  // Do some checks and init stuff
  
  // In case of several cone and thresholds analysis, open the cuts for the filling of the
  // track and cluster reference arrays in cone when done in the MakeAnalysisFillAOD().
  // The different cones, thresholds are tested for this list of tracks, clusters.
  if(fMakeSeveralIC)
  {
    printf("AliAnaParticleIsolation::Init() - Open default isolation cuts for multiple Isolation analysis\n");
    GetIsolationCut()->SetPtThreshold(100);
    GetIsolationCut()->SetPtFraction(100);
    GetIsolationCut()->SetConeSize(1);
  }
  
  if(!GetReader()->IsCTSSwitchedOn() && GetIsolationCut()->GetParticleTypeInCone()!=AliIsolationCut::kOnlyNeutral)
    AliFatal("STOP!: You want to use CTS tracks in analysis but not read!! \n!!Check the configuration file!!\n");
  
}

//____________________________________________
void AliAnaParticleIsolation::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  SetInputAODName("PWG4Particle");
  SetAODObjArrayName("IsolationCone");
  AddToHistogramsName("AnaIsolation_");
  
  fCalorimeter = "EMCAL" ;
  fIsoDetector = "EMCAL" ;
  
  fReMakeIC = kFALSE ;
  fMakeSeveralIC = kFALSE ;
  
  fLeadingOnly = kTRUE;
  fCheckLeadingWithNeutralClusters = kTRUE;
  
  fNDecayBits = 1;
  fDecayBits[0] = AliNeutralMesonSelection::kPi0;
  fDecayBits[1] = AliNeutralMesonSelection::kEta;
  fDecayBits[2] = AliNeutralMesonSelection::kPi0Side;
  fDecayBits[3] = AliNeutralMesonSelection::kEtaSide;
  
  //----------- Several IC-----------------
  fNCones             = 5 ;
  fNPtThresFrac       = 5 ;
  fConeSizes      [0] = 0.1;     fConeSizes      [1] = 0.2;   fConeSizes      [2] = 0.3; fConeSizes      [3] = 0.4;  fConeSizes      [4] = 0.5;
  fPtThresholds   [0] = 1.;      fPtThresholds   [1] = 2.;    fPtThresholds   [2] = 3.;  fPtThresholds   [3] = 4.;   fPtThresholds   [4] = 5.;
  fPtFractions    [0] = 0.05;    fPtFractions    [1] = 0.075; fPtFractions    [2] = 0.1; fPtFractions    [3] = 1.25; fPtFractions    [4] = 1.5;
  fSumPtThresholds[0] = 1.;      fSumPtThresholds[1] = 2.;    fSumPtThresholds[2] = 3.;  fSumPtThresholds[3] = 4.;   fSumPtThresholds[4] = 5.;
  
}

//_________________________________________________________________________________________
Bool_t AliAnaParticleIsolation::IsTriggerTheNearSideEventLeadingParticle(Int_t & idLeading)
{
  // Check if the what of the selected isolation candidates is leading particle in the same hemisphere
  // comparing with all the candidates, all the tracks or all the clusters.
  
  Double_t ptTrig      = GetMinPt();
  Double_t phiTrig     = 0 ;
  Int_t index          =-1 ;
  AliAODPWG4ParticleCorrelation* pLeading = 0;
  
  // Loop on stored AOD particles, find leading trigger on the selected list, with at least min pT.
  
  for(Int_t iaod = 0; iaod < GetInputAODBranch()->GetEntriesFast() ; iaod++)
  {
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    particle->SetLeadingParticle(kFALSE); // set it later
    
    // Vertex cut in case of mixing
    if(GetMixedEvent())
    {
      Int_t check = CheckMixedEventVertex(particle->GetCaloLabel(0), particle->GetTrackLabel(0));
      if(check ==  0) continue;
      if(check == -1) return kFALSE; // not sure if it is correct.
    }
    
    //check if it is low pt trigger particle
    if((particle->Pt() < GetIsolationCut()->GetPtThreshold() ||
        particle->Pt() < GetIsolationCut()->GetSumPtThreshold()) &&
       !fMakeSeveralIC)
    {
      continue ; //trigger should not come from underlying event
    }

    // find the leading particles with highest momentum
    if (particle->Pt() > ptTrig)
    {
      ptTrig   = particle->Pt() ;
      phiTrig  = particle->Phi();
      index    = iaod     ;
      pLeading = particle ;
    }
  }// finish search of leading trigger particle on the AOD branch.
  
  if(index < 0) return kFALSE;
  
  idLeading = index;
  
  //printf("AOD leading pT %2.2f, ID %d\n", pLeading->Pt(),pLeading->GetCaloLabel(0));
  
  if(phiTrig < 0 ) phiTrig += TMath::TwoPi();
  
  // Compare if it is the leading of all tracks
  
  TVector3 p3;
  for(Int_t ipr = 0;ipr < GetCTSTracks()->GetEntriesFast() ; ipr ++ )
  {
    AliVTrack * track = (AliVTrack *) (GetCTSTracks()->At(ipr)) ;
    
    if(track->GetID() == pLeading->GetTrackLabel(0) || track->GetID() == pLeading->GetTrackLabel(1) ||
       track->GetID() == pLeading->GetTrackLabel(2) || track->GetID() == pLeading->GetTrackLabel(3)   ) continue ;
    
    Double_t mom[3] = {track->Px(),track->Py(),track->Pz()};
    p3.SetXYZ(mom[0],mom[1],mom[2]);
    Float_t pt   = p3.Pt();
    Float_t phi  = p3.Phi() ;
    if(phi < 0) phi+=TMath::TwoPi();
    
    //skip this event if near side associated particle pt larger than trigger
    
    if(pt > ptTrig && TMath::Abs(phi-phiTrig) < TMath::PiOver2())  return kFALSE;
    
  }// track loop
  
  // Compare if it is leading of all calorimeter clusters
  
  if(fCheckLeadingWithNeutralClusters)
  {
    // Select the calorimeter cluster list
    TObjArray * nePl = 0x0;
    if      (pLeading->GetDetector() == "PHOS" )
      nePl = GetPHOSClusters();
    else
      nePl = GetEMCALClusters();
    
    if(!nePl) return kTRUE; // Do the selection just with the tracks if no calorimeter is available.
    
    TLorentzVector lv;
    for(Int_t ipr = 0;ipr < nePl->GetEntriesFast() ; ipr ++ )
    {
      AliVCluster * cluster = (AliVCluster *) (nePl->At(ipr)) ;
      
      if(cluster->GetID() == pLeading->GetCaloLabel(0) || cluster->GetID() == pLeading->GetCaloLabel(1) ) continue ;
      
      cluster->GetMomentum(lv,GetVertex(0));
      
      Float_t pt   = lv.Pt();
      Float_t phi  = lv.Phi() ;
      if(phi < 0) phi+=TMath::TwoPi();
      
      if(IsTrackMatched(cluster,GetReader()->GetInputEvent())) continue ; // avoid charged clusters, already covered by tracks, or cluster merging with track.
      
      // skip this event if near side associated particle pt larger than trigger
      // not really needed for calorimeter, unless DCal is included
     
      if(pt > ptTrig && TMath::Abs(phi-phiTrig) < TMath::PiOver2()) return kFALSE ;

    }// cluster loop
  } // check neutral clusters
  
  idLeading = index ;
  pLeading->SetLeadingParticle(kTRUE);
  
  if( GetDebug() > 1 )
    printf("AliAnaParticleIsolation::IsTriggerTheNearSideEventLeadingParticle() - Particle AOD with index %d is leading with pT %2.2f\n",
           idLeading, pLeading->Pt());
  
  return kTRUE;
  
}

//__________________________________________________
void  AliAnaParticleIsolation::MakeAnalysisFillAOD()
{
  // Do analysis and fill aods
  // Search for the isolated photon in fCalorimeter with GetMinPt() < pt < GetMaxPt()
  // and if the particle is leading in the near side (if requested)
  
  if(!GetInputAODBranch())
    AliFatal(Form("No input particles in AOD with name branch < %s >, STOP",GetInputAODName().Data()));
  
  if(strcmp(GetInputAODBranch()->GetClass()->GetName(), "AliAODPWG4ParticleCorrelation"))
    AliFatal(Form("Wrong type of AOD object, change AOD class name in input AOD: It should be <AliAODPWG4ParticleCorrelation> and not <%s> \n",GetInputAODBranch()->GetClass()->GetName()));
  
  Int_t n = 0, nfrac = 0;
  Bool_t  isolated  = kFALSE ;
  Float_t coneptsum = 0, coneptlead = 0;
  TObjArray * pl    = 0x0; ;
  
  //Select the calorimeter for candidate isolation with neutral particles
  if      (fCalorimeter == "PHOS" )
    pl = GetPHOSClusters();
  else if (fCalorimeter == "EMCAL")
    pl = GetEMCALClusters();
  
  //Loop on AOD branch, filled previously in AliAnaPhoton, find leading particle to do isolation only with it
  TLorentzVector mom ;
  Int_t idLeading = -1 ;
  Int_t iaod0 = 0;
  Int_t naod  = GetInputAODBranch()->GetEntriesFast();
  
  if(GetDebug() > 0)
    printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() - Input aod branch entries %d\n", naod);
  
  if(IsLeadingOnlyOn())
  {
    Bool_t leading = IsTriggerTheNearSideEventLeadingParticle(idLeading);
    if(!leading)
    {
      if(GetDebug() > 1) printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() - Not leading; End fill AODs \n");
      return;
    }
    iaod0 = idLeading  ; // first entry in particle loop
    naod  = idLeading+1; // last entry in particle loop
  }
  
  // Check isolation of list of candidate particles or leading particle
  
  for(Int_t iaod = iaod0; iaod < naod; iaod++)
  {
    AliAODPWG4ParticleCorrelation * aodinput =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    if(IsLeadingOnlyOn()) aodinput->SetLeadingParticle(kTRUE);
  
    // Check isolation only of clusters in fiducial region
    if(IsFiducialCutOn())
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(*aodinput->Momentum(),aodinput->GetDetector()) ;
      if(! in ) return ;
    }
    
    //If too small or too large pt, skip
    Float_t pt = aodinput->Pt();
    if(pt < GetMinPt() || pt > GetMaxPt() ) continue ;

    //check if it is low pt trigger particle
    if( ( pt < GetIsolationCut()->GetPtThreshold() ||  pt < GetIsolationCut()->GetSumPtThreshold() ) &&
       !fMakeSeveralIC)
    {
      continue ; //trigger should not come from underlying event
    }
    
    //After cuts, study isolation
    n=0; nfrac = 0; isolated = kFALSE; coneptsum = 0; coneptlead = 0;
    GetIsolationCut()->MakeIsolationCut(GetCTSTracks(),pl,
                                        GetReader(), GetCaloPID(),
                                        kTRUE, aodinput, GetAODObjArrayName(),
                                        n,nfrac,coneptsum,coneptlead,isolated);
    
    if(!fMakeSeveralIC) aodinput->SetIsolated(isolated);
  } // particle isolationloop
  
  if(GetDebug() > 1)
  {
    if(isolated) printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() : Particle %d IS ISOLATED \n",idLeading);
    printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() - End fill AODs \n");
  }
  
}

//_________________________________________________________
void  AliAnaParticleIsolation::MakeAnalysisFillHistograms()
{
  // Do analysis and fill histograms
  
  // In case of simulated data, fill acceptance histograms
  // Not ready for multiple case analysis.
  if(IsDataMC() && !fMakeSeveralIC) FillAcceptanceHistograms();
  
  //Loop on stored AOD
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0)
    printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Histo aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliAODPWG4ParticleCorrelation* aod =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    if(IsLeadingOnlyOn() && !aod->IsLeadingParticle()) continue; // Try to isolate only leading cluster or track
    
    // Check isolation only of clusters in fiducial region
    if(IsFiducialCutOn())
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(*aod->Momentum(),aod->GetDetector()) ;
      if(! in ) continue ;
    }
    
    Float_t pt         = aod->Pt();
    
    //If too small or too large pt, skip
    if(pt < GetMinPt() || pt > GetMaxPt() ) continue ;
    
    Int_t mcTag        = aod->GetTag() ;
    Int_t mcIndex      = GetMCIndex(mcTag);
    
    // --- In case of redoing isolation from delta AOD ----
    // Not standard case, not used since its implementation
    if(fMakeSeveralIC)
    {
      //Analysis of multiple IC at same time
      MakeSeveralICAnalysis(aod,mcIndex);
      continue;
    }
    
    // --- In case of redoing isolation multiple cuts ----
    
    if(fReMakeIC)
    {
      //In case a more strict IC is needed in the produced AOD
      Bool_t  isolated = kFALSE;
      Int_t   n = 0, nfrac = 0;
      Float_t coneptsum = 0, coneptlead = 0;
      
      //Recover reference arrays with clusters and tracks
      TObjArray * refclusters = aod->GetObjArray(GetAODObjArrayName()+"Clusters");
      TObjArray * reftracks   = aod->GetObjArray(GetAODObjArrayName()+"Tracks");
      
      GetIsolationCut()->MakeIsolationCut(reftracks,   refclusters,
                                          GetReader(), GetCaloPID(),
                                          kFALSE, aod, "",
                                          n,nfrac,coneptsum,coneptlead,isolated);
    }
    
    Bool_t  isolated   = aod->IsIsolated();
    Float_t energy     = aod->E();
    Float_t phi        = aod->Phi();
    Float_t eta        = aod->Eta();

    Int_t   decayTag = 0;
    if(fFillTaggedDecayHistograms)
    {
      decayTag = aod->GetBtag(); // temporary
      if(decayTag < 0) decayTag = 0; // temporary
    }
    
    if(GetDebug() > 0) printf(" AliAnaParticleIsolation::MakeAnalysisFillHistograms() - pt %1.1f, eta %1.1f, phi %1.1f, Isolated %d\n",
                              pt, eta, phi, isolated);
    
    //---------------------------------------------------------------
    // Fill Shower shape and track matching histograms
    //---------------------------------------------------------------
    
    FillTrackMatchingShowerShapeControlHistograms(aod,mcIndex);
    
    //---------------------------------------------------------------
    // Fill pt/sum pT distribution of particles in cone or in UE band
    //---------------------------------------------------------------
    
    Float_t coneptLeadCluster= 0;
    Float_t coneptLeadTrack  = 0;
    Float_t coneptsumCluster = 0;
    Float_t coneptsumTrack   = 0;
    Float_t coneptsumCell    = 0;
    Float_t etaBandptsumClusterNorm = 0;
    Float_t etaBandptsumTrackNorm   = 0;
    
    CalculateTrackSignalInCone   (aod,coneptsumTrack  , coneptLeadTrack  );
    CalculateCaloSignalInCone    (aod,coneptsumCluster, coneptLeadCluster);
    if(fFillCellHistograms)
      CalculateCaloCellSignalInCone(aod,coneptsumCell   );
    
    if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged)
    {
      fhConeSumPtClustervsTrack     ->Fill(coneptsumCluster,coneptsumTrack);
      if(fFillCellHistograms)
      {
        fhConeSumPtCellvsTrack        ->Fill(coneptsumCell,   coneptsumTrack);
        fhConeSumPtCellTrack          ->Fill(pt,     coneptsumTrack+coneptsumCell);
        fhConeSumPtCellTrackTrigEtaPhi->Fill(eta,phi,coneptsumTrack+coneptsumCell);
      }
    }
    
    fhConeSumPt              ->Fill(pt,     coneptsumTrack+coneptsumCluster);
    fhConeSumPtTrigEtaPhi    ->Fill(eta,phi,coneptsumTrack+coneptsumCluster);
    
    if(coneptLeadTrack > coneptLeadCluster) fhConePtLead->Fill(pt, coneptLeadTrack  );
    else                                    fhConePtLead->Fill(pt, coneptLeadCluster);
    
    if(GetDebug() > 1)
      printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Particle %d Energy Sum in Isolation Cone %2.2f\n", iaod, coneptsumTrack+coneptsumCluster);
    
    //normalize phi/eta band per area unit
    if(fFillUEBandSubtractHistograms)
      CalculateNormalizeUEBandPerUnitArea(aod, coneptsumCluster, coneptsumCell, coneptsumTrack, etaBandptsumTrackNorm, etaBandptsumClusterNorm) ;
    
    //  printf("Histograms analysis : cluster pt = %f, etaBandTrack = %f, etaBandCluster = %f, isolation = %d\n",aod->Pt(),etaBandptsumTrackNorm,etaBandptsumClusterNorm,aod->IsIsolated());
    
    
    //---------------------------------------------------------------
    // Isolated/ Non isolated histograms
    //---------------------------------------------------------------
    
    if(isolated)
    {
      if(GetDebug() > 1)
        printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Particle %d ISOLATED: fill histograms\n", iaod);
      
      fhEIso      ->Fill(energy);
      fhPtIso     ->Fill(pt);
      fhPhiIso    ->Fill(pt,phi);
      fhEtaIso    ->Fill(pt,eta);
      fhEtaPhiIso ->Fill(eta,phi);
      
      if(IsDataMC())
      {
        // For histograms in arrays, index in the array, corresponding to any particle origin
        if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))
        {
          fhPtIsoMC [kmcPhoton]->Fill(pt);
          fhPhiIsoMC[kmcPhoton]->Fill(pt,phi);
          fhEtaIsoMC[kmcPhoton]->Fill(pt,eta);
        }
        
        fhPtIsoMC [mcIndex]->Fill(pt);
        fhPhiIsoMC[mcIndex]->Fill(pt,phi);
        fhEtaIsoMC[mcIndex]->Fill(pt,eta);
      }//Histograms with MC
      
      // Candidates tagged as decay in another analysis (AliAnaPi0EbE)
      if(fFillTaggedDecayHistograms && decayTag > 0)
      {
        for(Int_t ibit = 0; ibit < fNDecayBits; ibit++)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[ibit]))
          {
            fhPtDecayIso[ibit]    ->Fill(pt);
            fhEtaPhiDecayIso[ibit]->Fill(eta,phi);
            
            if(IsDataMC())
            {
              if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))
                fhPtDecayIsoMC[ibit][kmcPhoton]->Fill(pt);
              
              fhPtDecayIsoMC[ibit][mcIndex]->Fill(pt);
            }
          } // bit ok
        } // bit loop
      } // decay histograms
      
      if(fFillNLMHistograms)
        fhPtNLocMaxIso ->Fill(pt,aod->GetFiducialArea()) ; // remember to change method name
      
      if(fFillHighMultHistograms)
      {
        fhPtCentralityIso ->Fill(pt,GetEventCentrality()) ;
        fhPtEventPlaneIso ->Fill(pt,GetEventPlaneAngle() ) ;
      }

      if(fFillPileUpHistograms)
      {
        if(GetReader()->IsPileUpFromSPD())               { fhEIsoPileUp[0] ->Fill(energy) ; fhPtIsoPileUp[0]->Fill(pt) ; }
        if(GetReader()->IsPileUpFromEMCal())             { fhEIsoPileUp[1] ->Fill(energy) ; fhPtIsoPileUp[1]->Fill(pt) ; }
        if(GetReader()->IsPileUpFromSPDOrEMCal())        { fhEIsoPileUp[2] ->Fill(energy) ; fhPtIsoPileUp[2]->Fill(pt) ; }
        if(GetReader()->IsPileUpFromSPDAndEMCal())       { fhEIsoPileUp[3] ->Fill(energy) ; fhPtIsoPileUp[3]->Fill(pt) ; }
        if(GetReader()->IsPileUpFromSPDAndNotEMCal())    { fhEIsoPileUp[4] ->Fill(energy) ; fhPtIsoPileUp[4]->Fill(pt) ; }
        if(GetReader()->IsPileUpFromEMCalAndNotSPD())    { fhEIsoPileUp[5] ->Fill(energy) ; fhPtIsoPileUp[5]->Fill(pt) ; }
        if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) { fhEIsoPileUp[6] ->Fill(energy) ; fhPtIsoPileUp[6]->Fill(pt) ; }
        
        // Fill histograms to undertand pile-up before other cuts applied
        // Remember to relax time cuts in the reader
        FillPileUpHistograms(aod->GetCaloLabel(0));
      }

    }//Isolated histograms
    else // NON isolated
    {
      if(GetDebug() > 1)
        printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Particle %d NOT ISOLATED, fill histograms\n", iaod);
      
      fhENoIso        ->Fill(energy);
      fhPtNoIso       ->Fill(pt);
      fhEtaPhiNoIso   ->Fill(eta,phi);
      
      if(IsDataMC())
      {
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
          fhPtNoIsoMC[kmcPhoton]->Fill(pt);
        
        fhPtNoIsoMC[mcIndex]->Fill(pt);
      }
      
      // Candidates tagged as decay in another analysis (AliAnaPi0EbE)
      if(fFillTaggedDecayHistograms && decayTag > 0)
      {
        for(Int_t ibit = 0; ibit < fNDecayBits; ibit++)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[ibit]))
          {
            fhPtDecayNoIso[ibit]    ->Fill(pt);
            fhEtaPhiDecayNoIso[ibit]->Fill(eta,phi);
            
            if(IsDataMC())
            {
              if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))
                fhPtDecayNoIsoMC[ibit][kmcPhoton]->Fill(pt);
              
              fhPtDecayNoIsoMC[ibit][mcIndex]->Fill(pt);
            }
          } // bit ok
        } // bit loop
      } // decay histograms
      
      if(fFillNLMHistograms)
        fhPtNLocMaxNoIso ->Fill(pt,aod->GetFiducialArea()); // remember to change method name
      
      if(fFillPileUpHistograms)
      {
        if(GetReader()->IsPileUpFromSPD())                { fhENoIsoPileUp[0] ->Fill(energy) ; fhPtNoIsoPileUp[0]->Fill(pt) ; }
        if(GetReader()->IsPileUpFromEMCal())              { fhENoIsoPileUp[1] ->Fill(energy) ; fhPtNoIsoPileUp[1]->Fill(pt) ; }
        if(GetReader()->IsPileUpFromSPDOrEMCal())         { fhENoIsoPileUp[2] ->Fill(energy) ; fhPtNoIsoPileUp[2]->Fill(pt) ; }
        if(GetReader()->IsPileUpFromSPDAndEMCal())        { fhENoIsoPileUp[3] ->Fill(energy) ; fhPtNoIsoPileUp[3]->Fill(pt) ; }
        if(GetReader()->IsPileUpFromSPDAndNotEMCal())     { fhENoIsoPileUp[4] ->Fill(energy) ; fhPtNoIsoPileUp[4]->Fill(pt) ; }
        if(GetReader()->IsPileUpFromEMCalAndNotSPD())     { fhENoIsoPileUp[5] ->Fill(energy) ; fhPtNoIsoPileUp[5]->Fill(pt) ; }
        if(GetReader()->IsPileUpFromNotSPDAndNotEMCal())  { fhENoIsoPileUp[6] ->Fill(energy) ; fhPtNoIsoPileUp[6]->Fill(pt) ; }
      }
    } // non iso
  }// aod loop

}

//______________________________________________________________________
void AliAnaParticleIsolation::FillAcceptanceHistograms()
{
  // Fill acceptance histograms if MC data is available
  // Only when particle is in the calorimeter. Rethink if CTS is used.
  
  if(GetDebug() > 0) printf("AliAnaParticleIsolation::FillAcceptanceHistograms() - Start \n");
  
  //Double_t photonY   = -100 ;
  Double_t photonE   = -1 ;
  Double_t photonPt  = -1 ;
  Double_t photonPhi =  100 ;
  Double_t photonEta = -1 ;
  
  Int_t    pdg       =  0 ;
  Int_t    status    =  0 ;
  Int_t    tag       =  0 ;
  Int_t    mcIndex   =  0 ;
  Int_t    nprim     = 0;
  
  TParticle        * primStack = 0;
  AliAODMCParticle * primAOD   = 0;
  TLorentzVector lv;
  
  // Get the ESD MC particles container
  AliStack * stack = 0;
  if( GetReader()->ReadStack() )
  {
    stack = GetMCStack();
    if(!stack ) return;
    nprim = stack->GetNtrack();
  }
  
  // Get the AOD MC particles container
  TClonesArray * mcparticles = 0;
  if( GetReader()->ReadAODMCParticles() )
  {
    mcparticles = GetReader()->GetAODMCParticles();
    if( !mcparticles ) return;
    nprim = mcparticles->GetEntriesFast();
  }
  
  for(Int_t i=0 ; i < nprim; i++)
  {
    if(GetReader()->AcceptOnlyHIJINGLabels() && !GetReader()->IsHIJINGLabel(i)) continue ;
    
    if(GetReader()->ReadStack())
    {
      primStack = stack->Particle(i) ;
      if(!primStack)
      {
        printf("AliAnaParticleIsolation::FillAcceptanceHistograms() - ESD primaries pointer not available!!\n");
        continue;
      }
      
      pdg    = primStack->GetPdgCode();
      status = primStack->GetStatusCode();
      
      if(primStack->Energy() == TMath::Abs(primStack->Pz()))  continue ; //Protection against floating point exception
      
      //printf("i %d, %s %d  %s %d \n",i, stack->Particle(i)->GetName(), stack->Particle(i)->GetPdgCode(),
      //       prim->GetName(), prim->GetPdgCode());
      
      //photonY   = 0.5*TMath::Log((prim->Energy()-prim->Pz())/(prim->Energy()+prim->Pz())) ;
      
      //Photon kinematics
      primStack->Momentum(lv);
      
    }
    else
    {
      primAOD = (AliAODMCParticle *) mcparticles->At(i);
      if(!primAOD)
      {
        printf("AliAnaParticleIsolation::FillAcceptanceHistograms() - AOD primaries pointer not available!!\n");
        continue;
      }
      
      pdg    = primAOD->GetPdgCode();
      status = primAOD->GetStatus();
      
      if(primAOD->E() == TMath::Abs(primAOD->Pz()))  continue ; //Protection against floating point exception
      
      //photonY   = 0.5*TMath::Log((prim->Energy()-prim->Pz())/(prim->Energy()+prim->Pz())) ;
      
      //Photon kinematics
      lv.SetPxPyPzE(primAOD->Px(),primAOD->Py(),primAOD->Pz(),primAOD->E());
    }
    
    // Select only photons in the final state
    if(pdg != 22 ) continue ;
    
    // Consider only final state particles, but this depends on generator,
    // status 1 is the usual one, in case of not being ok, leave the possibility
    // to not consider this.
    if(status != 1 && GetMCAnalysisUtils()->GetMCGenerator()!="" ) continue ;
    
    // If too small or too large pt, skip, same cut as for data analysis
    photonPt  = lv.Pt () ;
    
    if(photonPt < GetMinPt() || photonPt > GetMaxPt() ) continue ;
    
    photonE   = lv.E  () ;
    photonEta = lv.Eta() ;
    photonPhi = lv.Phi() ;
    
    if(photonPhi < 0) photonPhi+=TMath::TwoPi();
    
    // Check if photons hit the Calorimeter acceptance
    if(IsRealCaloAcceptanceOn() && fIsoDetector!="CTS") // defined on base class
    {
      if(GetReader()->ReadStack()          &&
         !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(fIsoDetector, primStack)) continue ;
      if(GetReader()->ReadAODMCParticles() &&
         !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(fIsoDetector, primAOD  )) continue ;
    }
    
    // Check same fidutial borders as in data analysis on top of real acceptance if real was requested.
    if(!GetFiducialCut()->IsInFiducialCut(lv,fIsoDetector)) continue ;
    
    // Get tag of this particle photon from fragmentation, decay, prompt ...
    // Set the origin of the photon.
    tag = GetMCAnalysisUtils()->CheckOrigin(i,GetReader());
    
    if(!GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
    {
      // A conversion photon from a hadron, skip this kind of photon
      // printf("AliAnaPhoton::FillAcceptanceHistograms() - not a photon, weird!\n ");
      // GetMCAnalysisUtils()->PrintMCTag(tag);
      
      continue;
    }
    
    //
    if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt) )
    {
      mcIndex = kmcPrimPrompt;
    }
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation) )
    {
      mcIndex = kmcPrimFrag ;
    }
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR))
    {
      mcIndex = kmcPrimISR;
    }
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))
    {
      mcIndex = kmcPrimPi0Decay;
    }
    else if( (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) ||
              GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay)))
    {
      mcIndex = kmcPrimOtherDecay;
    }
    else
    {
      // Other decay but from non final state particle
      mcIndex = kmcPrimOtherDecay;
    }//Other origin
    
    // ////////////////////ISO MC/////////////////////////
    Double_t sumPtInCone = 0; Double_t dR=0. ;
    TParticle        * mcisopStack = 0;
    AliAODMCParticle * mcisopAOD   = 0;
    Int_t partInConeStatus = -1, partInConeMother = -1;
    Double_t partInConePt = 0, partInConeEta = 0, partInConePhi = 0;
    Int_t npart = 0;
    for(Int_t ip = 0; ip < nprim ; ip++)
    {
      if(ip==i) continue;
      
      if( GetReader()->ReadStack() )
      {
        mcisopStack = static_cast<TParticle*>(stack->Particle(ip));
        if( !mcisopStack ) continue;
        partInConeStatus = mcisopStack->GetStatusCode();
        partInConeMother = mcisopStack->GetMother(0);
        partInConePt     = mcisopStack->Pt();
        partInConeEta    = mcisopStack->Eta();
        partInConePhi    = mcisopStack->Phi();
      }
      else
      {
        mcisopAOD   = (AliAODMCParticle *) mcparticles->At(ip);
        if( !mcisopAOD )   continue;
        partInConeStatus = mcisopAOD->GetStatus();
        partInConeMother = mcisopAOD->GetMother();
        partInConePt     = mcisopAOD->Pt();
        partInConeEta    = mcisopAOD->Eta();
        partInConePhi    = mcisopAOD->Phi();
      }
      
      // Consider only final state particles, but this depends on generator,
      // status 1 is the usual one, in case of not being ok, leave the possibility
      // to not consider this.
      if( partInConeStatus != 1 && GetMCAnalysisUtils()->GetMCGenerator()!="" ) continue ;
      
      if( partInConeMother == i ) continue;
      
      if( partInConePt < GetReader()->GetCTSPtMin() ) continue;
      // Careful!!!, cut for TPC tracks and calorimeter clusters in cone can be different
      
      // TVector3 vmcv(mcisop->Px(),mcisop->Py(), mcisop->Pz());
      // if(vmcv.Perp()>1)
      //   continue;
      
      //
      // Add here Acceptance cut???, charged particles CTS fid cut, neutral Calo real acceptance.
      //
      
      dR = GetIsolationCut()->Radius(photonEta, photonPhi, partInConeEta, partInConePhi);
      
      if(dR > GetIsolationCut()->GetConeSize())
        continue;
      
      sumPtInCone += partInConePt;
      if(partInConePt > GetIsolationCut()->GetPtThreshold() &&
         partInConePt < GetIsolationCut()->GetPtThresholdMax()) npart++;
    }
    
    ///////END ISO MC/////////////////////////
    
    // Fill the histograms, only those in the defined calorimeter acceptance
    
    fhEtaPrimMC[kmcPrimPhoton]->Fill(photonPt , photonEta) ;
    fhPhiPrimMC[kmcPrimPhoton]->Fill(photonPt , photonPhi) ;
    fhEPrimMC  [kmcPrimPhoton]->Fill(photonE) ;
    
    fhEtaPrimMC[mcIndex]->Fill(photonPt , photonEta) ;
    fhPhiPrimMC[mcIndex]->Fill(photonPt , photonPhi) ;
    fhEPrimMC  [mcIndex]->Fill(photonE ) ;
    
    // Isolated?
    Bool_t isolated = kFALSE;
    if(GetIsolationCut()->GetICMethod() == AliIsolationCut::kSumPtIC   &&
       (sumPtInCone < GetIsolationCut()->GetSumPtThreshold() ||
        sumPtInCone > GetIsolationCut()->GetSumPtThresholdMax()))
      isolated = kTRUE;
    
    if(GetIsolationCut()->GetICMethod() == AliIsolationCut::kPtThresIC &&
       npart == 0)
      isolated = kTRUE;
    
    if(isolated)
    {
      fhPtPrimMCiso [mcIndex]      ->Fill(photonPt) ;
      fhPtPrimMCiso [kmcPrimPhoton]->Fill(photonPt) ;
    }
    
  }//loop on primaries
  
  if(GetDebug() > 0) printf("AliAnaParticleIsolation::FillAcceptanceHistograms() - End \n");
  
}


//_____________________________________________________________________________________
void  AliAnaParticleIsolation::MakeSeveralICAnalysis(AliAODPWG4ParticleCorrelation* ph,
                                                     Int_t mcIndex)
{
  
  //Isolation Cut Analysis for both methods and different pt cuts and cones
  Float_t ptC   = ph->Pt();
  Float_t etaC  = ph->Eta();
  Float_t phiC  = ph->Phi();
  Int_t   tag   = ph->GetTag();

  Int_t   decayTag = 0;
  if(fFillTaggedDecayHistograms)
  {
    decayTag = ph->GetBtag(); // temporary
    if(decayTag < 0) decayTag = 0; // temporary
  }

  if(GetDebug() > 0)
    printf("AliAnaParticleIsolation::MakeSeveralICAnalysis() - Isolate pT %2.2f, decay tag %d\n",ptC, decayTag);
  
  //Keep original setting used when filling AODs, reset at end of analysis
  Float_t ptthresorg = GetIsolationCut()->GetPtThreshold();
  Float_t ptfracorg  = GetIsolationCut()->GetPtFraction();
  Float_t ptsumcorg  = GetIsolationCut()->GetSumPtThreshold();
  Float_t rorg       = GetIsolationCut()->GetConeSize();
  
  Float_t coneptsum = 0, coneptlead = 0;
  Int_t   n    [10][10];//[fNCones][fNPtThresFrac];
  Int_t   nfrac[10][10];//[fNCones][fNPtThresFrac];
  Bool_t  isolated  = kFALSE;
  
  // Fill hist with all particles before isolation criteria
  fhENoIso     ->Fill(ph->E());
  fhPtNoIso    ->Fill(ptC);
  fhEtaPhiNoIso->Fill(etaC,phiC);
  
  if(IsDataMC())
  {
    if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
      fhPtNoIsoMC[kmcPhoton]->Fill(ptC);
    
    fhPtNoIsoMC[mcIndex]->Fill(ptC);
  }
  
  // Candidates tagged as decay in another analysis (AliAnaPi0EbE)
  if(fFillTaggedDecayHistograms && decayTag > 0)
  {
    for(Int_t ibit = 0; ibit < fNDecayBits; ibit++)
    {
      if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[ibit]))
      {
        fhPtDecayNoIso[ibit]    ->Fill(ptC);
        fhEtaPhiDecayNoIso[ibit]->Fill(etaC,phiC);
        
        if(IsDataMC())
        {
          if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
            fhPtDecayNoIsoMC[ibit][kmcPhoton]->Fill(ptC);
          
          fhPtDecayNoIsoMC[ibit][mcIndex]->Fill(ptC);
        }
      } // bit ok
    } // bit loop
  } // decay histograms
  
  //Get vertex for photon momentum calculation
  Double_t vertex[] = {0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    GetReader()->GetVertex(vertex);
  
  //Loop on cone sizes
  for(Int_t icone = 0; icone<fNCones; icone++)
  {
    //Recover reference arrays with clusters and tracks
    TObjArray * refclusters = ph->GetObjArray(GetAODObjArrayName()+"Clusters");
    TObjArray * reftracks   = ph->GetObjArray(GetAODObjArrayName()+"Tracks");
    
    //If too small or too large pt, skip
    if(ptC < GetMinPt() || ptC > GetMaxPt() ) continue ;
    
    //In case a more strict IC is needed in the produced AOD
    
    isolated = kFALSE; coneptsum = 0; coneptlead = 0;
    
    GetIsolationCut()->SetSumPtThreshold(100);
    GetIsolationCut()->SetPtThreshold(100);
    GetIsolationCut()->SetPtFraction(100);
    GetIsolationCut()->SetConeSize(fConeSizes[icone]);
    
    // retreive pt tracks to fill histo vs. pt leading
    //Fill pt distribution of particles in cone
    //fhPtLeadingPt(),fhPerpSumPtLeadingPt(),fhPerpPtLeadingPt(),
    
    // Tracks in perpendicular cones
    Double_t sumptPerp = 0. ;
    TObjArray * trackList   = GetCTSTracks() ;
    for(Int_t itrack=0; itrack < trackList->GetEntriesFast(); itrack++)
    {
      AliVTrack* track = (AliVTrack *) trackList->At(itrack);
      //fill the histograms at forward range
      if(!track)
      {
        printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Track not available?");
        continue;
      }
      
      Double_t dPhi = phiC - track->Phi() + TMath::PiOver2();
      Double_t dEta = etaC - track->Eta();
      Double_t arg  = dPhi*dPhi + dEta*dEta;
      if(TMath::Sqrt(arg) < fConeSizes[icone])
      {
        fhPerpPtLeadingPt[icone]->Fill(ptC,TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py()));
        sumptPerp+=track->Pt();
      }
      
      dPhi = phiC - track->Phi() - TMath::PiOver2();
      arg  = dPhi*dPhi + dEta*dEta;
      if(TMath::Sqrt(arg) < fConeSizes[icone])
      {
        fhPerpPtLeadingPt[icone]->Fill(ptC,TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py()));
        sumptPerp+=track->Pt();
      }
    }
    
    fhPerpSumPtLeadingPt[icone]->Fill(ptC,sumptPerp);
    
    // Tracks in isolation cone, pT distribution and sum
    if(reftracks && GetIsolationCut()->GetParticleTypeInCone()!= AliIsolationCut::kOnlyNeutral)
    {
      for(Int_t itrack=0; itrack < reftracks->GetEntriesFast(); itrack++)
      {
        AliVTrack* track = (AliVTrack *) reftracks->At(itrack);
        
        Float_t rad = GetIsolationCut()->Radius(etaC, phiC, track->Eta(), track->Phi());
        
        if(rad > fConeSizes[icone]) continue ;
        
        fhPtLeadingPt[icone]->Fill(ptC, track->Pt());
        coneptsum += track->Pt();
      }
    }
    
    // Clusters in isolation cone, pT distribution and sum
    if(refclusters && GetIsolationCut()->GetParticleTypeInCone()!= AliIsolationCut::kOnlyCharged)
    {
      TLorentzVector mom ;
      for(Int_t icalo=0; icalo < refclusters->GetEntriesFast(); icalo++)
      {
        AliVCluster* calo = (AliVCluster *) refclusters->At(icalo);
        
        calo->GetMomentum(mom,vertex) ;//Assume that come from vertex in straight line
        
        Float_t rad = GetIsolationCut()->Radius(etaC, phiC, mom.Eta(), mom.Phi());
        
        if(rad > fConeSizes[icone]) continue ;
        
        fhPtLeadingPt[icone]->Fill(ptC, mom.Pt());
        coneptsum += mom.Pt();
      }
    }
    
    fhSumPtLeadingPt[icone]->Fill(ptC,coneptsum);
    
    if(IsDataMC())
    {
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
        fhSumPtLeadingPtMC[kmcPhoton][icone]->Fill(ptC,coneptsum) ;
      
      fhSumPtLeadingPtMC[mcIndex][icone]->Fill(ptC,coneptsum) ;
    }
    
    ///////////////////
    
    //Loop on pt thresholds
    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++)
    {
      n    [icone][ipt]=0;
      nfrac[icone][ipt]=0;
      GetIsolationCut()->SetPtThreshold(fPtThresholds[ipt]);
      GetIsolationCut()->SetPtFraction(fPtFractions[ipt]) ;
      GetIsolationCut()->SetSumPtThreshold(fSumPtThresholds[ipt]);
      
      GetIsolationCut()->MakeIsolationCut(reftracks, refclusters,
                                          GetReader(), GetCaloPID(),
                                          kFALSE, ph, "",
                                          n[icone][ipt],nfrac[icone][ipt],
                                          coneptsum, coneptlead, isolated);
      
      // Normal pT threshold cut
      
      if(GetDebug() > 0)
      {
        printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - cone size %1.1f, ptThres  %1.1f, sumptThresh  %1.1f\n",
               fConeSizes[icone],fPtThresholds[ipt],fSumPtThresholds[ipt]);
        printf("\t n %d, nfrac %d, coneptsum %2.2f\n",
               n[icone][ipt],nfrac[icone][ipt],coneptsum);
        
        printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - pt %1.1f, eta %1.1f, phi %1.1f\n",ptC, etaC, phiC);
      }
      
      if(n[icone][ipt] == 0)
      {
        if(GetDebug() > 0)
          printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - filling pt threshold loop\n");
        
        fhPtThresIsolated [icone][ipt]->Fill(ptC);
        fhEtaPhiPtThresIso[icone][ipt]->Fill(etaC,phiC);
        
        if( fFillTaggedDecayHistograms && decayTag > 0 && fNDecayBits > 0)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[0]))
          {
            fhPtPtThresDecayIso    [icone][ipt]->Fill(ptC);
            fhEtaPhiPtThresDecayIso[icone][ipt]->Fill(etaC,phiC);
          }
        }
        
        if(IsDataMC())
        {
          if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) )
            fhPtThresIsolatedMC[kmcPhoton][icone][ipt]->Fill(ptC) ;
          
          fhPtThresIsolatedMC[mcIndex][icone][ipt]->Fill(ptC) ;
          
        }
      }
      
      // pt in cone fraction
      if(nfrac[icone][ipt] == 0)
      {
        if(GetDebug() > 0)
          printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - filling frac loop\n");
        
        fhPtFracIsolated [icone][ipt]->Fill(ptC);
        fhEtaPhiPtFracIso[icone][ipt]->Fill(etaC,phiC);
        
        if( fFillTaggedDecayHistograms && decayTag > 0 && fNDecayBits > 0)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[0]))
          {
            fhPtPtFracDecayIso    [icone][ipt]->Fill(ptC);
            fhEtaPhiPtFracDecayIso[icone][ipt]->Fill(etaC,phiC);
          }
        }
        
        if(IsDataMC())
        {
          if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
            fhPtFracIsolatedMC[kmcPhoton][icone][ipt]->Fill(ptC) ;
          
          fhPtFracIsolatedMC[mcIndex][icone][ipt]->Fill(ptC) ;
        }
      }
      
      if(GetDebug()>0)
        printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - checking IC method : %i\n",GetIsolationCut()->GetICMethod());
      
      //Pt threshold on pt cand/ sum in cone histograms
      if(coneptsum < fSumPtThresholds[ipt])
      {
        if(GetDebug() > 0 )
          printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - filling sum loop\n");
        
        fhSumPtIsolated [icone][ipt]->Fill(ptC) ;
        fhEtaPhiPtSumIso[icone][ipt]->Fill(etaC, phiC) ;
        
        if( fFillTaggedDecayHistograms && decayTag > 0 && fNDecayBits > 0)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[0]))
          {
            fhPtPtSumDecayIso[icone][ipt]->Fill(ptC);
            fhEtaPhiPtSumDecayIso[icone][ipt]->Fill(etaC, phiC) ;
          }
        }
        
        if(IsDataMC())
        {
          if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
            fhSumPtIsolatedMC[kmcPhoton][icone][ipt]->Fill(ptC) ;
          
          fhSumPtIsolatedMC[mcIndex][icone][ipt]->Fill(ptC) ;
        }
      }
      
      // pt sum pt frac method
      //    if( ((fPtFractions[ipt]*ptC < fSumPtThresholds[ipt]) && (coneptsum < fSumPtThresholds[ipt])) || ((fPtFractions[ipt]*ptC > fSumPtThresholds[ipt]) && (coneptsum < fPtFractions[ipt]*ptC)) )
      
      if(coneptsum < fPtFractions[ipt]*ptC)
      {
        if(GetDebug() > 0)
          printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - filling PtFrac PtSum loop\n");
        
        fhPtFracPtSumIso    [icone][ipt]->Fill(ptC) ;
        fhEtaPhiFracPtSumIso[icone][ipt]->Fill(etaC,phiC) ;
        
        if( fFillTaggedDecayHistograms && decayTag > 0 && fNDecayBits > 0)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[0]))
          {
            fhPtFracPtSumDecayIso    [icone][ipt]->Fill(ptC);
            fhEtaPhiFracPtSumDecayIso[icone][ipt]->Fill(etaC,phiC);
          }
        }
      }
      
      // density method
      Float_t cellDensity = GetIsolationCut()->GetCellDensity( ph, GetReader());
      if(coneptsum < fSumPtThresholds[ipt]*cellDensity)
      {
        if(GetDebug() > 0)
          printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - filling density loop\n");
        
        fhPtSumDensityIso    [icone][ipt]->Fill(ptC) ;
        fhEtaPhiSumDensityIso[icone][ipt]->Fill(etaC,phiC) ;
        
        if( fFillTaggedDecayHistograms && decayTag > 0 && fNDecayBits > 0)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[0]))
          {
            fhPtSumDensityDecayIso    [icone][ipt]->Fill(ptC);
            fhEtaPhiSumDensityDecayIso[icone][ipt]->Fill(etaC, phiC);
          }
        }
      }
    }//pt thresh loop
    
    
  }//cone size loop
  
  //Reset original parameters for AOD analysis
  GetIsolationCut()->SetPtThreshold(ptthresorg);
  GetIsolationCut()->SetPtFraction(ptfracorg);
  GetIsolationCut()->SetSumPtThreshold(ptsumcorg);
  GetIsolationCut()->SetConeSize(rorg);
  
}

//_____________________________________________________________
void AliAnaParticleIsolation::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  
  printf("ReMake Isolation          = %d \n",  fReMakeIC) ;
  printf("Make Several Isolation    = %d \n",  fMakeSeveralIC) ;
  printf("Calorimeter for isolation = %s \n",  fCalorimeter.Data()) ;
  printf("Detector for candidate isolation = %s \n", fIsoDetector.Data()) ;
  
  if(fMakeSeveralIC)
  {
    printf("N Cone Sizes       =     %d\n", fNCones) ;
    printf("Cone Sizes          =    \n") ;
    for(Int_t i = 0; i < fNCones; i++)
      printf("  %1.2f;",  fConeSizes[i]) ;
    printf("    \n") ;
    
    printf("N pT thresholds/fractions = %d\n", fNPtThresFrac) ;
    printf(" pT thresholds         =    \n") ;
    for(Int_t i = 0; i < fNPtThresFrac; i++)
      printf("   %2.2f;",  fPtThresholds[i]) ;
    
    printf("    \n") ;
    
    printf(" pT fractions         =    \n") ;
    for(Int_t i = 0; i < fNPtThresFrac; i++)
      printf("   %2.2f;",  fPtFractions[i]) ;
    
    printf("    \n") ;
    
    printf("sum pT thresholds         =    \n") ;
    for(Int_t i = 0; i < fNPtThresFrac; i++)
      printf("   %2.2f;",  fSumPtThresholds[i]) ;
    
    
  }
  
  printf("    \n") ;
  
} 

