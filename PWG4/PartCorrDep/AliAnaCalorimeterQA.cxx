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
/* $Id: $ */

//_________________________________________________________________________
// Class to check results from simulations or reconstructed real data. 
// Fill few histograms and do some checking plots
//
//-- Author: Gustavo Conesa (INFN-LNF)
//_________________________________________________________________________


// --- ROOT system ---
//#include "Riostream.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TROOT.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TLegend.h"
#include <TObjString.h>

//---- AliRoot system ----
#include "AliAnaCalorimeterQA.h"
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliVCaloCells.h"
#include "AliFiducialCut.h"
#include "AliAODTrack.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAODMCParticle.h"
#include "AliMCAnalysisUtils.h"
#include "AliAODPid.h"
#include "AliExternalTrackParam.h"

ClassImp(AliAnaCalorimeterQA)

//________________________________________
AliAnaCalorimeterQA::AliAnaCalorimeterQA() : 
AliAnaPartCorrBaseClass(),             fCalorimeter(""), 

//Switches
fFillAllPosHisto(kFALSE),              fFillAllPosHisto2(kTRUE), 
fFillAllTH12(kFALSE),                  fFillAllTH3(kTRUE), 
fFillAllTMHisto(kTRUE),                fFillAllPi0Histo(kTRUE),                 
fCorrelate(kTRUE),                     fStudyBadClusters(kFALSE),               
fStudyClustersAsymmetry(kFALSE),       fStudyWeight(kFALSE),

//Parameters and cuts
fNModules(12),                         fNRCU(2),
fNMaxCols(48),                         fNMaxRows(24),  
fTimeCutMin(-1),                       fTimeCutMax(9999999),
fEMCALCellAmpMin(0),                   fPHOSCellAmpMin(0), 

//Histograms
fhE(0),                                fhPt(0),                                
fhPhi(0),                              fhEta(0),                               fhEtaPhiE(0),
fhECharged(0),                         fhPtCharged(0),             
fhPhiCharged(0),                       fhEtaCharged(0),                        fhEtaPhiECharged(0), 

//Invariant mass
fhIM(0 ),                              fhAsym(0), 
fhNCellsPerCluster(0),                 fhNCellsPerClusterNoCut(0), 
fhNCellsPerClusterMIP(0),              fhNCellsPerClusterMIPCharged(0), 
fhNCellsvsClusterMaxCellDiffE0(0),     fhNCellsvsClusterMaxCellDiffE2(0),      
fhNCellsvsClusterMaxCellDiffE6(0),     fhNClusters(0),    

//Timing
fhClusterTimeEnergy(0),                fhCellTimeSpreadRespectToCellMax(0),  
fhCellIdCellLargeTimeSpread(0),        fhClusterPairDiffTimeE(0),
fhClusterMaxCellCloseCellRatio(0),     fhClusterMaxCellCloseCellDiff(0), 
fhClusterMaxCellDiff(0),               fhClusterMaxCellDiffNoCut(0), 
fhClusterMaxCellDiffAverageTime(0),    fhClusterMaxCellDiffAverageNoMaxTime(0),
fhClusterMaxCellDiffWeightedTime(0),   fhClusterMaxCellDiffWeightedNoMaxTime(0), 
fhLambda0vsClusterMaxCellDiffE0(0),    fhLambda0vsClusterMaxCellDiffE2(0),     fhLambda0vsClusterMaxCellDiffE6(0),
fhLambda0(0),                          fhLambda1(0),                           fhDispersion(0),

//bad clusters
fhBadClusterEnergy(0),                 fhBadClusterTimeEnergy(0),              
fhBadClusterPairDiffTimeE(0),          fhBadCellTimeSpreadRespectToCellMax(0), 
fhBadClusterMaxCellCloseCellRatio(0),  fhBadClusterMaxCellCloseCellDiff(0),    fhBadClusterMaxCellDiff(0),
fhBadClusterMaxCellDiffAverageTime(0), fhBadClusterMaxCellDiffAverageNoMaxTime(0), 
fhBadClusterMaxCellDiffWeightedTime(0),fhBadClusterMaxCellDiffWeightedNoMaxTime(0), 

//Position
fhRNCells(0),                          fhXNCells(0),               
fhYNCells(0),                          fhZNCells(0),
fhRE(0),                               fhXE(0),                    
fhYE(0),                               fhZE(0),    
fhXYZ(0),
fhRCellE(0),                           fhXCellE(0),                
fhYCellE(0),                           fhZCellE(0),
fhXYZCell(0),
fhDeltaCellClusterRNCells(0),          fhDeltaCellClusterXNCells(0),
fhDeltaCellClusterYNCells(0),          fhDeltaCellClusterZNCells(0),
fhDeltaCellClusterRE(0),               fhDeltaCellClusterXE(0),     
fhDeltaCellClusterYE(0),               fhDeltaCellClusterZE(0),

// Cells
fhNCells(0),                           fhAmplitude(0),             
fhAmpId(0),                            fhEtaPhiAmp(0), 
fhTime(0),                             fhTimeId(0),                            fhTimeAmp(0), 
fhCaloCorrNClusters(0),                fhCaloCorrEClusters(0),     
fhCaloCorrNCells(0),                   fhCaloCorrECells(0),
fhCaloV0SCorrNClusters(0),             fhCaloV0SCorrEClusters(0),              
fhCaloV0SCorrNCells(0),                fhCaloV0SCorrECells(0),
fhCaloV0MCorrNClusters(0),             fhCaloV0MCorrEClusters(0),  
fhCaloV0MCorrNCells(0),                fhCaloV0MCorrECells(0),
fhCaloTrackMCorrNClusters(0),          fhCaloTrackMCorrEClusters(0), 
fhCaloTrackMCorrNCells(0),             fhCaloTrackMCorrECells(0),

//Super-Module dependent histgrams
fhEMod(0),                             fhAmpMod(0),                            fhTimeMod(0),  
fhNClustersMod(0),                     fhNCellsMod(0),
fhNCellsPerClusterMod(0),              fhNCellsPerClusterModNoCut(0), 

fhGridCells(0),                        fhGridCellsE(0),                        fhGridCellsTime(0), 
fhTimeAmpPerRCU(0),                    fhIMMod(0),              

// Weight studies
fhECellClusterRatio(0),                fhECellClusterLogRatio(0),                 
fhEMaxCellClusterRatio(0),             fhEMaxCellClusterLogRatio(0),                

// MC and reco
fhRecoMCE(),                           fhRecoMCPhi(),                          fhRecoMCEta(), 
fhRecoMCDeltaE(),                      fhRecoMCRatioE(),                      
fhRecoMCDeltaPhi(),                    fhRecoMCDeltaEta(),               

// MC only
fhGenMCE(),                            fhGenMCEtaPhi(),   
fhGenMCAccE(),                         fhGenMCAccEtaPhi(),   

//matched MC
fhEMVxyz(0),                           fhEMR(0),                   
fhHaVxyz(0),                           fhHaR(0),
fh1pOverE(0),                          fh1dR(0),                   
fh2EledEdx(0),                         fh2MatchdEdx(0),
fhMCEle1pOverE(0),                     fhMCEle1dR(0),                          fhMCEle2MatchdEdx(0),
fhMCChHad1pOverE(0),                   fhMCChHad1dR(0),                        fhMCChHad2MatchdEdx(0),
fhMCNeutral1pOverE(0),                 fhMCNeutral1dR(0),                      fhMCNeutral2MatchdEdx(0), fh1pOverER02(0),           
fhMCEle1pOverER02(0),                  fhMCChHad1pOverER02(0),                 fhMCNeutral1pOverER02(0)
{
  //Default Ctor
  
  //Weight studies
  for(Int_t i =0; i < 7; i++){
    fhLambda0ForW0[i] = 0;
    fhLambda1ForW0[i] = 0;
    
    for(Int_t j = 0; j < 5; j++){
      fhLambda0ForW0MC[i][j] = 0;
      fhLambda1ForW0MC[i][j] = 0;
    }
    
  }
  
  //Cluster size
  fhDeltaIEtaDeltaIPhiE0[0] = 0 ;         fhDeltaIEtaDeltaIPhiE2[0] = 0;          fhDeltaIEtaDeltaIPhiE6[0] = 0; 
  fhDeltaIEtaDeltaIPhiE0[1] = 0 ;         fhDeltaIEtaDeltaIPhiE2[1] = 0;          fhDeltaIEtaDeltaIPhiE6[1] = 0; 
  fhDeltaIA[0]              = 0 ;         fhDeltaIAL0[0]            = 0;          fhDeltaIAL1[0]            = 0;
  fhDeltaIA[1]              = 0 ;         fhDeltaIAL0[1]            = 0;          fhDeltaIAL1[1]            = 0;                         
  fhDeltaIANCells[0]        = 0 ;         fhDeltaIANCells[1]        = 0;
  fhDeltaIAMC[0]            = 0 ;         fhDeltaIAMC[1]            = 0;
  fhDeltaIAMC[2]            = 0 ;         fhDeltaIAMC[3]            = 0;
  
  // MC
  
  for(Int_t i = 0; i < 6; i++){
    
    fhRecoMCE[i][0]         = 0; fhRecoMCE[i][1]        = 0;  
    fhRecoMCPhi[i][0]       = 0; fhRecoMCPhi[i][1]      = 0;  
    fhRecoMCEta[i][0]       = 0; fhRecoMCEta[i][1]      = 0;  
    fhRecoMCDeltaE[i][0]    = 0; fhRecoMCDeltaE[i][1]   = 0;  
    fhRecoMCRatioE[i][0]    = 0; fhRecoMCRatioE[i][1]   = 0;  
    fhRecoMCDeltaPhi[i][0]  = 0; fhRecoMCDeltaPhi[i][1] = 0;  
    fhRecoMCDeltaEta[i][0]  = 0; fhRecoMCDeltaEta[i][1] = 0;  
    
  }
  
  //Initialize parameters
  InitParameters();
}

//_______________________________________________________________________________________________________________
void AliAnaCalorimeterQA::BadClusterHistograms(AliVCluster* clus, TObjArray *caloClusters, AliVCaloCells * cells, 
                                               const Int_t absIdMax, const Double_t maxCellFraction,
                                               const Double_t tmax,  Double_t timeAverages[4]
                                               )
{
  //Bad cluster histograms
  
  fhBadClusterEnergy     ->Fill(clus->E());
  Double_t tof = clus->GetTOF()*1.e9;
  fhBadClusterTimeEnergy ->Fill(clus->E(),tof);
  fhBadClusterMaxCellDiff->Fill(clus->E(),maxCellFraction);
  //Clusters in event time difference
  
  for(Int_t iclus2 = 0; iclus2 < caloClusters->GetEntriesFast(); iclus2++ ){
    
    AliVCluster* clus2 =  (AliVCluster*)caloClusters->At(iclus2);
    
    if(clus->GetID()==clus2->GetID()) continue;
    
    if(clus->GetM02() > 0.01 && clus2->GetM02() > 0.01) {
      
      Double_t tof2   = clus->GetTOF()*1.e9;      
      fhBadClusterPairDiffTimeE  ->Fill(clus->E(), (tof-tof2));
      
    }
  } // loop
  
  // Max cell compared to other cells in cluster
  if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
    fhBadClusterMaxCellDiffAverageTime      ->Fill(clus->E(),tmax-timeAverages[0]);
    fhBadClusterMaxCellDiffAverageNoMaxTime ->Fill(clus->E(),tmax-timeAverages[2]);
    fhBadClusterMaxCellDiffWeightedTime     ->Fill(clus->E(),tmax-timeAverages[1]);
    fhBadClusterMaxCellDiffWeightedNoMaxTime->Fill(clus->E(),tmax-timeAverages[3]);
  }           
  
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) {
    Int_t absId  = clus->GetCellsAbsId()[ipos]; 
    if(absId!=absIdMax){
      
      Float_t frac = cells->GetCellAmplitude(absId)/cells->GetCellAmplitude(absIdMax);
      
      fhBadClusterMaxCellCloseCellRatio->Fill(clus->E(),frac);
      fhBadClusterMaxCellCloseCellDiff ->Fill(clus->E(),cells->GetCellAmplitude(absIdMax)-cells->GetCellAmplitude(absId));
      
      if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
        Double_t time  = cells->GetCellTime(absId);
        RecalibrateCellTime(time,absId);
        
        Float_t diff = (tmax-time*1e9);
        fhBadCellTimeSpreadRespectToCellMax->Fill(clus->E(), diff);
        
      } // ESD
    }// Not max
  }//loop  
  
}

//___________________________________________________________________________________________________________
void AliAnaCalorimeterQA::CalculateAverageTime(AliVCluster *clus, AliVCaloCells* cells,  Double_t timeAverages[4])
{
  // Calculate time averages and weights
  
  // First recalculate energy in case non linearity was applied
  Float_t  energy = 0;
  Float_t  ampMax = 0, amp = 0;
  Int_t    absIdMax =-1;
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) {
    
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    amp = cells->GetCellAmplitude(id);
    RecalibrateCellAmplitude(amp,id);
    
    energy    += amp;
    
    if(amp> ampMax) {
      ampMax   = amp;
      absIdMax = id;
    }
    
  } // energy loop       
  
  // Calculate average time of cells in cluster and weighted average
  Double_t aTime = 0; Double_t aTimeNoMax = 0;
  Double_t wTime = 0; Double_t wTimeNoMax = 0;
  Float_t  wTot  = 0;
  Double_t time  = 0;
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) {
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    amp = cells->GetCellAmplitude(id);
    time  = cells->GetCellTime(id);
    
    //Recalibrate energy and time
    RecalibrateCellAmplitude(amp, id);    
    RecalibrateCellTime     (time,id);
    
    Double_t w    = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(cells->GetCellAmplitude(id),energy);
    aTime += time*1e9;
    wTime += time*1e9 * w;
    wTot  += w;
    
    if(id != absIdMax){
      aTimeNoMax += time*1e9;
      wTimeNoMax += time*1e9 * w;
    }
  }        
  
  aTime  /= clus->GetNCells();
  if(clus->GetNCells() > 1 ) 
    aTimeNoMax /= (clus->GetNCells()-1);
  
  if(wTot>0){
    wTime      /= wTot;
    wTimeNoMax /= wTot;
  }
  
  timeAverages[0] = aTime;        timeAverages[1] = wTime; 
  timeAverages[2] = aTimeNoMax;   timeAverages[3] = wTimeNoMax; 
  
}

//____________________________________________________________
void AliAnaCalorimeterQA::CellHistograms(AliVCaloCells *cells)
{
  // Plot histograms related to cells only
  
  Int_t ncells = cells->GetNumberOfCells();
  
  if(GetDebug() > 0) 
    printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - %s cell entries %d\n", fCalorimeter.Data(), ncells );    
  
  //Init arrays and used variables
  Int_t *nCellsInModule = new Int_t[fNModules];
  for(Int_t imod = 0; imod < fNModules; imod++ ) nCellsInModule[imod] = 0;
  Int_t    icol   = -1;
  Int_t    irow   = -1;
  Int_t    iRCU   = -1;
  Float_t  amp    = 0.;
  Double_t time   = 0.;
  Int_t    id     = -1;
  Float_t  recalF = 1.;  
  
  for (Int_t iCell = 0; iCell < cells->GetNumberOfCells(); iCell++) {      
    if(GetDebug() > 2)  
      printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - Cell : amp %f, absId %d \n", cells->GetAmplitude(iCell), cells->GetCellNumber(iCell));
    Int_t nModule = GetModuleNumberCellIndexes(cells->GetCellNumber(iCell),fCalorimeter, icol, irow, iRCU);
    if(GetDebug() > 2) 
      printf("\t module %d, column %d, row %d \n", nModule,icol,irow);
    
    if(nModule < fNModules) {	
      
      //Check if the cell is a bad channel
      if(GetCaloUtils()->IsBadChannelsRemovalSwitchedOn()){
        if(fCalorimeter=="EMCAL"){
          if(GetCaloUtils()->GetEMCALChannelStatus(nModule,icol,irow)) continue;
        }
        else {
          if(GetCaloUtils()->GetPHOSChannelStatus(nModule,icol,irow)) {
            printf("PHOS bad channel\n");
            continue;
          }
        }
      } // use bad channel map
      
      amp     = cells->GetAmplitude(iCell)*recalF;
      time    = cells->GetTime(iCell);
      id      = cells->GetCellNumber(iCell);
      
      // Amplitude recalibration if set
      RecalibrateCellAmplitude(amp,id);
      
      // Time recalibration if set
      RecalibrateCellTime(time,id);
      
      //Transform time to ns
      time *= 1.0e9;
      
      //Remove noisy channels, only possible in ESDs
      if(GetReader()->GetDataType() == AliCaloTrackReader::kESD){
        if(time < fTimeCutMin || time > fTimeCutMax){
          if(GetDebug() > 0 )printf("AliAnaCalorimeterQA - Remove cell with Time %f\n",time);
          continue;
        }
      }
      
      fhAmplitude->Fill(amp);
      fhAmpId    ->Fill(amp,id);
      fhAmpMod   ->Fill(amp,nModule);
      
      if ((fCalorimeter=="EMCAL" && amp > fEMCALCellAmpMin) ||
          (fCalorimeter=="PHOS"  && amp > fPHOSCellAmpMin))   {
        
        nCellsInModule[nModule]++ ;
        
        Int_t icols = icol;
        Int_t irows = irow;
        if(fCalorimeter=="EMCAL"){
          icols = (nModule % 2) ? icol + fNMaxCols : icol;				
          irows = irow + fNMaxRows * Int_t(nModule / 2);
        }
        else {
          irows = irow + fNMaxRows * fNModules;
        }
        
        fhGridCells ->Fill(icols,irows);
        fhGridCellsE->Fill(icols,irows,amp);
        
        if(GetReader()->GetDataType() == AliCaloTrackReader::kESD){
          //printf("%s: time %g\n",fCalorimeter.Data(), time);
          fhTime     ->Fill(time);
          fhTimeId   ->Fill(time,id);
          fhTimeAmp  ->Fill(amp,time);
          fhGridCellsTime->Fill(icols,irows,time);
          fhTimeMod  ->Fill(time,nModule);
          fhTimeAmpPerRCU  [nModule*fNRCU+iRCU]->Fill(amp, time);
          
        }
      }
      
      //Get Eta-Phi position of Cell
      if(fFillAllPosHisto)
      {
        if(fCalorimeter=="EMCAL" && GetCaloUtils()->IsEMCALGeoMatrixSet()){
          Float_t celleta = 0.;
          Float_t cellphi = 0.;
          GetEMCALGeometry()->EtaPhiFromIndex(id, celleta, cellphi); 
          
          fhEtaPhiAmp->Fill(celleta,cellphi,amp);
          Double_t cellpos[] = {0, 0, 0};
          GetEMCALGeometry()->GetGlobal(id, cellpos);
          fhXCellE->Fill(cellpos[0],amp)  ; 
          fhYCellE->Fill(cellpos[1],amp)  ; 
          fhZCellE->Fill(cellpos[2],amp)  ;
          Float_t rcell = TMath::Sqrt(cellpos[0]*cellpos[0]+cellpos[1]*cellpos[1]);//+cellpos[2]*cellpos[2]);
          fhRCellE->Fill(rcell,amp)  ;
          fhXYZCell->Fill(cellpos[0],cellpos[1],cellpos[2])  ;
        }//EMCAL Cells
        else if(fCalorimeter=="PHOS" && GetCaloUtils()->IsPHOSGeoMatrixSet()){
          TVector3 xyz;
          Int_t relId[4], module;
          Float_t xCell, zCell;
          
          GetPHOSGeometry()->AbsToRelNumbering(id,relId);
          module = relId[0];
          GetPHOSGeometry()->RelPosInModule(relId,xCell,zCell);
          GetPHOSGeometry()->Local2Global(module,xCell,zCell,xyz);
          Float_t rcell = TMath::Sqrt(xyz.X()*xyz.X()+xyz.Y()*xyz.Y());
          fhXCellE ->Fill(xyz.X(),amp)  ; 
          fhYCellE ->Fill(xyz.Y(),amp)  ; 
          fhZCellE ->Fill(xyz.Z(),amp)  ;
          fhRCellE ->Fill(rcell  ,amp)  ;
          fhXYZCell->Fill(xyz.X(),xyz.Y(),xyz.Z())  ;
        }//PHOS cells
      }//fill cell position histograms
      
      if     (fCalorimeter=="EMCAL" && amp > fEMCALCellAmpMin) ncells ++ ;
      else if(fCalorimeter=="PHOS"  && amp > fPHOSCellAmpMin)  ncells ++ ;
      //else  
      //  printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - no %s CELLS passed the analysis cut\n",fCalorimeter.Data());    
    }//nmodules
  }//cell loop
  
  if(ncells > 0 )fhNCells->Fill(ncells) ; //fill the cells after the cut 
  
  //Number of cells per module
  for(Int_t imod = 0; imod < fNModules; imod++ ) {
    
    if(GetDebug() > 1) 
      printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - module %d calo %s cells %d\n", imod, fCalorimeter.Data(), nCellsInModule[imod]); 
    
    fhNCellsMod->Fill(nCellsInModule[imod],imod) ;
    
  }
  
  delete [] nCellsInModule;
  
}

//__________________________________________________________________________
void AliAnaCalorimeterQA::CellInClusterPositionHistograms(AliVCluster* clus)
{
  // Fill histograms releated to cell position
  
  
  Int_t nCaloCellsPerCluster = clus->GetNCells();
  UShort_t * indexList = clus->GetCellsAbsId();
  Float_t pos[3];
  clus->GetPosition(pos);
  Float_t clEnergy = clus->E();
  
  //Loop on cluster cells
  for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) {
    
    //	printf("Index %d\n",ipos);
    Int_t absId  = indexList[ipos]; 
    
    //Get position of cell compare to cluster
    
    if(fCalorimeter=="EMCAL" && GetCaloUtils()->IsEMCALGeoMatrixSet()){
      
      Double_t cellpos[] = {0, 0, 0};
      GetEMCALGeometry()->GetGlobal(absId, cellpos);
      
      fhDeltaCellClusterXNCells->Fill(pos[0]-cellpos[0],nCaloCellsPerCluster) ; 
      fhDeltaCellClusterYNCells->Fill(pos[1]-cellpos[1],nCaloCellsPerCluster) ; 
      fhDeltaCellClusterZNCells->Fill(pos[2]-cellpos[2],nCaloCellsPerCluster) ;
      
      fhDeltaCellClusterXE->Fill(pos[0]-cellpos[0],clEnergy)  ; 
      fhDeltaCellClusterYE->Fill(pos[1]-cellpos[1],clEnergy)  ; 
      fhDeltaCellClusterZE->Fill(pos[2]-cellpos[2],clEnergy)  ; 
      
      Float_t r     = TMath::Sqrt(pos[0]    *pos[0]     + pos[1]    * pos[1]    );
      Float_t rcell = TMath::Sqrt(cellpos[0]*cellpos[0] + cellpos[1]* cellpos[1]);
      
      fhDeltaCellClusterRNCells->Fill(r-rcell, nCaloCellsPerCluster) ; 
      fhDeltaCellClusterRE     ->Fill(r-rcell, clEnergy)  ; 			
      
    }//EMCAL and its matrices are available
    else if(fCalorimeter=="PHOS" && GetCaloUtils()->IsPHOSGeoMatrixSet()){
      TVector3 xyz;
      Int_t relId[4], module;
      Float_t xCell, zCell;
      
      GetPHOSGeometry()->AbsToRelNumbering(absId,relId);
      module = relId[0];
      GetPHOSGeometry()->RelPosInModule(relId,xCell,zCell);
      GetPHOSGeometry()->Local2Global(module,xCell,zCell,xyz);
      
      fhDeltaCellClusterXNCells->Fill(pos[0]-xyz.X(),nCaloCellsPerCluster) ; 
      fhDeltaCellClusterYNCells->Fill(pos[1]-xyz.Y(),nCaloCellsPerCluster) ; 
      fhDeltaCellClusterZNCells->Fill(pos[2]-xyz.Z(),nCaloCellsPerCluster) ;
      
      fhDeltaCellClusterXE->Fill(pos[0]-xyz.X(),clEnergy)  ; 
      fhDeltaCellClusterYE->Fill(pos[1]-xyz.Y(),clEnergy)  ; 
      fhDeltaCellClusterZE->Fill(pos[2]-xyz.Z(),clEnergy)  ; 
      
      Float_t r     = TMath::Sqrt(pos[0]  * pos[0]  + pos[1]  * pos[1] );
      Float_t rcell = TMath::Sqrt(xyz.X() * xyz.X() + xyz.Y() * xyz.Y());
      
      fhDeltaCellClusterRNCells->Fill(r-rcell, nCaloCellsPerCluster) ; 
      fhDeltaCellClusterRE     ->Fill(r-rcell, clEnergy)  ; 
      
    }//PHOS and its matrices are available
  }// cluster cell loop
}

//___________________________________________________________________________________________
void AliAnaCalorimeterQA::ClusterAsymmetryHistograms(AliVCluster* clus, const Int_t absIdMax)
{
  // Study the shape of the cluster in cell units terms
  
  //No use to study clusters with less than 4 cells
  if(clus->GetNCells() <=3 ) return;
  
  Int_t dIeta = 0;
  Int_t dIphi = 0;
  
  Int_t ietaMax=-1; Int_t iphiMax = 0; Int_t rcuMax = 0;
  Int_t smMax = GetModuleNumberCellIndexes(absIdMax,fCalorimeter, ietaMax, iphiMax, rcuMax);
  
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) {
    
    Int_t absId = clus->GetCellsAbsId()[ipos];
    
    Int_t ieta=-1; Int_t iphi = 0; Int_t rcu = 0;
    Int_t sm = GetModuleNumberCellIndexes(absId,fCalorimeter, ieta, iphi, rcu);
    
    if(dIphi < TMath::Abs(iphi-iphiMax)) dIphi = TMath::Abs(iphi-iphiMax);
    
    if(smMax==sm){
      if(dIeta < TMath::Abs(ieta-ietaMax)) dIeta = TMath::Abs(ieta-ietaMax);
    }
    else {
      Int_t ietaShift    = ieta;
      Int_t ietaMaxShift = ietaMax;
      if (ieta > ietaMax)  ietaMaxShift+=48;
      else                 ietaShift   +=48;
      if(dIeta < TMath::Abs(ietaShift-ietaMaxShift)) dIeta = TMath::Abs(ietaShift-ietaMaxShift);
    }
    
    //if(TMath::Abs(clus->GetM20()) < 0.0001 && clus->GetNCells() > 3){
    //  printf("Good : E %f, mcells %d, l0 %f, l1 %f, d %f, cell max t %f, cluster TOF %f, sm %d, icol %d, irow %d; Max icol %d, irow %d \n", 
    //            clus->E(), clus->GetNCells(),clus->GetM02(), clus->GetM20(), clus->GetDispersion(),tmax, tof,sm,ieta,iphi, ietaMax, iphiMax);
    //}
    
  }// fill cell-cluster histogram loop
  
  // Was cluster matched?
  Bool_t matched = GetCaloPID()->IsTrackMatched(clus,GetCaloUtils());
  
  if     (clus->E() < 2 ) fhDeltaIEtaDeltaIPhiE0[matched]->Fill(dIeta,dIphi);
  else if(clus->E() < 6 ) fhDeltaIEtaDeltaIPhiE2[matched]->Fill(dIeta,dIphi);
  else                    fhDeltaIEtaDeltaIPhiE6[matched]->Fill(dIeta,dIphi);
  
  Float_t dIA    = 1.*(dIphi-dIeta)/(dIeta+dIphi);
  fhDeltaIA[matched]->Fill(clus->E(),dIA);
  
  if(clus->E() > 0.5){
    
    fhDeltaIAL0[matched]->Fill(clus->GetM02(),dIA);
    fhDeltaIAL1[matched]->Fill(clus->GetM20(),dIA);
    fhDeltaIANCells[matched]->Fill(clus->GetNCells(),dIA);
    
  }
  
  // Origin of  clusters
  Int_t  nLabel = clus->GetNLabels();
  Int_t* labels = clus->GetLabels();
  if(IsDataMC()){
    Int_t tag = GetMCAnalysisUtils()->CheckOrigin(labels,nLabel, GetReader(),0);
    if(   GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton) && 
       !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0)    && 
       !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta)    &&
       !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)        ){
      fhDeltaIAMC[0]->Fill(clus->E(),dIA);//Pure Photon
    }
    else if ( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCElectron) &&
             !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)  ){
      fhDeltaIAMC[1]->Fill(clus->E(),dIA);//Pure electron
    }
    else if ( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)  ){
      fhDeltaIAMC[2]->Fill(clus->E(),dIA);//Converted cluster
    }
    else if(!GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton)){ 
      fhDeltaIAMC[3]->Fill(clus->E(),dIA);//Hadrons
    }
    
  }  // MC
  
}

//___________________________________________________________________________________________________________
void AliAnaCalorimeterQA::ClusterHistograms(AliVCluster* clus,TObjArray *caloClusters, AliVCaloCells * cells, 
                                            const Int_t absIdMax, const Double_t maxCellFraction,
                                            const Double_t tmax,  Double_t timeAverages[4])
{
  //Fill CaloCluster related histograms
  
  Int_t    nCaloCellsPerCluster = clus->GetNCells();
  Int_t    nModule = GetModuleNumber(clus);
  Double_t tof     = clus->GetTOF()*1.e9;
  
  if     (clus->E() < 2.){
    fhLambda0vsClusterMaxCellDiffE0->Fill(clus->GetM02(),      maxCellFraction);
    fhNCellsvsClusterMaxCellDiffE0 ->Fill(nCaloCellsPerCluster,maxCellFraction);
  }
  else if(clus->E() < 6.){
    fhLambda0vsClusterMaxCellDiffE2->Fill(clus->GetM02(),      maxCellFraction);
    fhNCellsvsClusterMaxCellDiffE2 ->Fill(nCaloCellsPerCluster,maxCellFraction);
  }
  else{
    fhLambda0vsClusterMaxCellDiffE6->Fill(clus->GetM02(),      maxCellFraction);  
    fhNCellsvsClusterMaxCellDiffE6 ->Fill(nCaloCellsPerCluster,maxCellFraction);
  }
  
  fhNCellsPerClusterNoCut  ->Fill(clus->E(), nCaloCellsPerCluster);
  if(nModule >=0 && nModule < fNModules) fhNCellsPerClusterModNoCut[nModule]->Fill(clus->E(), nCaloCellsPerCluster);
  
  fhClusterMaxCellDiffNoCut->Fill(clus->E(),maxCellFraction);
  
  
  //Check bad clusters if requested and rejection was not on
  if(!IsGoodCluster(clus->E(),nCaloCellsPerCluster)) return;
  
  fhLambda0    ->Fill(clus->E(),clus->GetM02());
  fhLambda1    ->Fill(clus->E(),clus->GetM20());
  fhDispersion ->Fill(clus->E(),clus->GetDispersion());
  
  fhClusterMaxCellDiff->Fill(clus->E(),maxCellFraction);        
  fhClusterTimeEnergy ->Fill(clus->E(),tof);
  
  //Clusters in event time difference
  for(Int_t iclus2 = 0; iclus2 < caloClusters->GetEntriesFast(); iclus2++ ){
    
    AliVCluster* clus2 =  (AliVCluster*) caloClusters->At(iclus2);
    
    if(clus->GetID()==clus2->GetID()) continue;
    
    if(clus->GetM02() > 0.01 && clus2->GetM02() > 0.01) {
      Double_t tof2   = clus->GetTOF()*1.e9;          
      fhClusterPairDiffTimeE  ->Fill(clus->E(), tof-tof2);
    }
  }        
  
  if(nCaloCellsPerCluster > 1){
    
    // check time of cells respect to max energy cell
    
    if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
      fhClusterMaxCellDiffAverageTime      ->Fill(clus->E(),tmax-timeAverages[0]);
      fhClusterMaxCellDiffAverageNoMaxTime ->Fill(clus->E(),tmax-timeAverages[2]);
      fhClusterMaxCellDiffWeightedTime     ->Fill(clus->E(),tmax-timeAverages[1]);
      fhClusterMaxCellDiffWeightedNoMaxTime->Fill(clus->E(),tmax-timeAverages[3]);
    }
    
    for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) {
      
      Int_t absId  = clus->GetCellsAbsId()[ipos];             
      if(absId == absIdMax) continue;
      
      Float_t frac = cells->GetCellAmplitude(absId)/cells->GetCellAmplitude(absIdMax);            
      fhClusterMaxCellCloseCellRatio->Fill(clus->E(),frac);
      fhClusterMaxCellCloseCellDiff ->Fill(clus->E(),cells->GetCellAmplitude(absIdMax)-cells->GetCellAmplitude(absId));
      
      if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
        
        Double_t time  = cells->GetCellTime(absId);
        RecalibrateCellTime(time,absId);
        
        Float_t diff = (tmax-time*1.0e9);
        fhCellTimeSpreadRespectToCellMax->Fill(clus->E(), diff);
        if(TMath::Abs(TMath::Abs(diff) > 100) && clus->E() > 1 ) fhCellIdCellLargeTimeSpread->Fill(absId);
      }
      
    }// fill cell-cluster histogram loop
    
  }//check time and energy of cells respect to max energy cell if cluster of more than 1 cell
  
  
  // Get vertex for photon momentum calculation and event selection
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  
  TLorentzVector mom  ;
  clus->GetMomentum(mom,v);
  
  Float_t e   = mom.E();
  Float_t pt  = mom.Pt();
  Float_t eta = mom.Eta();
  Float_t phi = mom.Phi();
  if(phi < 0) phi +=TMath::TwoPi();
  
  if(GetDebug() > 0) {
    printf("AliAnaCalorimeterQA::ClusterHistograms() - cluster: E %2.3f, pT %2.3f, eta %2.3f, phi %2.3f \n",e,pt,eta,phi*TMath::RadToDeg());
  }
  
  fhE     ->Fill(e);	
  if(nModule >=0 && nModule < fNModules) fhEMod->Fill(e,nModule);
  if(fFillAllTH12){
    fhPt     ->Fill(pt);
    fhPhi    ->Fill(phi);
    fhEta    ->Fill(eta);
  }
  
  if(fFillAllTH3)
    fhEtaPhiE->Fill(eta,phi,e);
  
  //Cells per cluster
  fhNCellsPerCluster   ->Fill(e, nCaloCellsPerCluster);
  if((fCalorimeter=="EMCAL" && GetReader()->GetEMCALPtMin() < 0.3) ||
     (fCalorimeter=="PHOS"  && GetReader()->GetPHOSPtMin()  < 0.3)) fhNCellsPerClusterMIP->Fill(e, nCaloCellsPerCluster);
  
  //Position
  if(fFillAllPosHisto2){
    
    Float_t pos[3] ;     
    clus->GetPosition(pos);
    
    fhXE     ->Fill(pos[0],e);
    fhYE     ->Fill(pos[1],e);
    fhZE     ->Fill(pos[2],e);
    if(fFillAllTH3)
      fhXYZ    ->Fill(pos[0], pos[1],pos[2]);
    
    fhXNCells->Fill(pos[0],nCaloCellsPerCluster);
    fhYNCells->Fill(pos[1],nCaloCellsPerCluster);
    fhZNCells->Fill(pos[2],nCaloCellsPerCluster);
    Float_t rxyz = TMath::Sqrt(pos[0]*pos[0]+pos[1]*pos[1]);//+pos[2]*pos[2]);
    fhRE     ->Fill(rxyz,e);
    fhRNCells->Fill(rxyz  ,nCaloCellsPerCluster);
  }
  
  if(nModule >=0 && nModule < fNModules) fhNCellsPerClusterMod[nModule]->Fill(e, nCaloCellsPerCluster);
  
}

//_____________________________________________________________________________________________
void AliAnaCalorimeterQA::ClusterLoopHistograms(TObjArray *caloClusters, AliVCaloCells* cells)
{
  // Fill clusters related histograms
  
  TLorentzVector mom  ;
  Int_t  nLabel                = 0  ;
  Int_t *labels                = 0x0;
  Int_t  nCaloClusters         = caloClusters->GetEntriesFast() ;
  Int_t  nCaloClustersAccepted = 0  ;
  Int_t  nCaloCellsPerCluster  = 0  ;
  Bool_t matched               = kFALSE;
  Int_t  nModule               =-1  ;
  
  // Get vertex for photon momentum calculation and event selection
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  
  Int_t *nClustersInModule     = new Int_t[fNModules];
  for(Int_t imod = 0; imod < fNModules; imod++ ) nClustersInModule[imod] = 0;
  
  if(GetDebug() > 0)
    printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - In %s there are %d clusters \n", fCalorimeter.Data(), nCaloClusters);
  
  // Loop over CaloClusters
  for(Int_t iclus = 0; iclus < nCaloClusters; iclus++){
    
    if(GetDebug() > 0) printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - cluster: %d/%d, data %d \n",
                              iclus+1,nCaloClusters,GetReader()->GetDataType());
    
    AliVCluster* clus =  (AliVCluster*)caloClusters->At(iclus);
    
    // Get the fraction of the cluster energy that carries the cell with highest energy and its absId
    Float_t maxCellFraction = 0.;
    Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(cells, clus,maxCellFraction);
    
    //Cut on time of clusters
    Double_t tof = clus->GetTOF()*1.e9;
    if(tof < fTimeCutMin || tof > fTimeCutMax){ 
      if(GetDebug() > 0 )printf("AliAnaCalorimeterQA - Remove cluster with TOF %f\n",tof);
      continue;
    }    
    
    // Get cluster kinematics
    clus->GetMomentum(mom,v);
    
    // Check only certain regions
    Bool_t in = kTRUE;
    if(IsFiducialCutOn()) in =  GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
    if(!in) continue;
    
    nLabel = clus->GetNLabels();
    labels = clus->GetLabels();
    
    // Cells per cluster
    nCaloCellsPerCluster = clus->GetNCells();
    
    // Cluster mathed with track?
    matched = GetCaloPID()->IsTrackMatched(clus,GetCaloUtils());
    
    // Get some time averages
    Double_t averTime[4] = {0.,0.,0.,0.};
    CalculateAverageTime(clus, cells, averTime);
    
    //Get time of max cell
    Double_t tmax  = cells->GetCellTime(absIdMax);
    RecalibrateCellTime(tmax,absIdMax);
    tmax*=1.e9;
    
    //Check bad clusters if requested and rejection was not on
    Bool_t goodCluster = IsGoodCluster(clus->E(),nCaloCellsPerCluster);
    
    // Fill histograms related to single cluster 
    
    if(!goodCluster) 
      BadClusterHistograms(clus, caloClusters, cells, absIdMax, 
                           maxCellFraction, tmax, averTime);     
    
    ClusterHistograms(clus, caloClusters, cells, absIdMax, 
                      maxCellFraction, tmax, averTime);	
    
    if(!goodCluster) continue;
    
    nCaloClustersAccepted++;
    nModule = GetModuleNumber(clus);
    if(nModule >=0 && nModule < fNModules) {
      if     (fCalorimeter=="EMCAL" && mom.E() > 2*fEMCALCellAmpMin)  nClustersInModule[nModule]++;
      else if(fCalorimeter=="PHOS"  && mom.E() > 2*fPHOSCellAmpMin )  nClustersInModule[nModule]++;
    }  
    
    //Cluster size with respect to cell with maximum energy in cell units
    if(fStudyClustersAsymmetry) ClusterAsymmetryHistograms(clus,absIdMax);
    
    // Cluster weights
    if(fStudyWeight) WeightHistograms(clus, cells);
    
    // Cells in cluster position
    if(fFillAllPosHisto) CellInClusterPositionHistograms(clus);
    
    // Fill histograms related to single cluster, mc vs data
    Int_t  mcOK = kFALSE;
    Int_t  pdg  = -1;
    if(IsDataMC() && nLabel > 0 && labels) 
      mcOK = ClusterMCHistograms(mom, matched, labels, nLabel, pdg);
    
    // Matched clusters with tracks, also do some MC comparison, needs input from ClusterMCHistograms
    if( matched &&  fFillAllTMHisto)
      ClusterMatchedWithTrackHistograms(clus,mom,mcOK,pdg);	        
    
    // Invariant mass
    if(fFillAllPi0Histo && nCaloClusters > 1 && nCaloCellsPerCluster > 1) 
      InvariantMassHistograms(iclus, mom, nModule, caloClusters);
    
  }//cluster loop
  
  // Number of clusters histograms
  if(nCaloClustersAccepted > 0) fhNClusters->Fill(nCaloClustersAccepted);
  
  // Number of clusters per module
  for(Int_t imod = 0; imod < fNModules; imod++ ){ 
    if(GetDebug() > 1) 
      printf("AliAnaCalorimeterQA::ClusterLoopHistograms() - module %d calo %s clusters %d\n", imod, fCalorimeter.Data(), nClustersInModule[imod]); 
    fhNClustersMod->Fill(nClustersInModule[imod],imod);
  }
  
  delete [] nClustersInModule;
  
}

//_____________________________________________________________________________________________
Bool_t AliAnaCalorimeterQA::ClusterMCHistograms(const TLorentzVector mom, const Bool_t matched,
                                                const Int_t * labels, const Int_t nLabels, Int_t & pdg )
{
  
  //Fill histograms only possible when simulation
  
  if(GetDebug() > 1) {
    printf("\t Primaries: nlabels %d\n",nLabels);
    if(!nLabels || !labels) printf("\t Strange, no labels!!!\n");
  }  
  
  Float_t e   = mom.E();
  Float_t eta = mom.Eta();
  Float_t phi = mom.Phi();
  if(phi < 0) phi +=TMath::TwoPi();
  
  AliAODMCParticle * aodprimary  = 0x0;
  TParticle * primary = 0x0;
  
  //Play with the MC stack if available
  Int_t label = labels[0];
  
  if(label < 0) {
    if(GetDebug() >= 0) printf("AliAnaCalorimeterQA::ClusterHistograms() *** bad label ***:  label %d \n", label);
    return kFALSE;
  }
  
  Int_t pdg0  =-1;Int_t status = -1; Int_t iMother = -1; Int_t iParent = -1;
  Float_t vxMC= 0; Float_t vyMC = 0;	
  Float_t eMC = 0; Float_t ptMC= 0; Float_t phiMC =0; Float_t etaMC = 0;
  Int_t charge = 0;	
  
  //Check the origin.
  Int_t tag = GetMCAnalysisUtils()->CheckOrigin(labels,nLabels, GetReader(),0);
  
  if(GetReader()->ReadStack() && !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCUnknown)){ //it MC stack and known tag
    
    if( label >= GetMCStack()->GetNtrack()) {
      if(GetDebug() >= 0) printf("AliAnaCalorimeterQA::ClusterHistograms() *** large label ***:  label %d, n tracks %d \n", label, GetMCStack()->GetNtrack());
      return kFALSE;
    }
    
    primary = GetMCStack()->Particle(label);
    iMother = label;
    pdg0    = TMath::Abs(primary->GetPdgCode());
    pdg     = pdg0;
    status  = primary->GetStatusCode();
    vxMC    = primary->Vx();
    vyMC    = primary->Vy();
    iParent = primary->GetFirstMother();
    
    if(GetDebug() > 1 ) {
      printf("AliAnaCalorimeterQA::ClusterHistograms() - Cluster most contributing mother: \n");
      printf("\t Mother label %d, pdg %d, %s, status %d, parent %d \n",iMother, pdg0, primary->GetName(),status, iParent);
    }
    
    //Get final particle, no conversion products
    if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)){
      //Get the parent
      primary = GetMCStack()->Particle(iParent);
      pdg = TMath::Abs(primary->GetPdgCode());
      if(GetDebug() > 1 ) printf("AliAnaCalorimeterQA::ClusterHistograms() - Converted cluster!. Find before conversion: \n");
      while((pdg == 22 || pdg == 11) && status != 1){
        iMother = iParent;
        primary = GetMCStack()->Particle(iMother);
        status  = primary->GetStatusCode();
        iParent = primary->GetFirstMother();
        pdg     = TMath::Abs(primary->GetPdgCode());
        if(GetDebug() > 1 )printf("\t pdg %d, index %d, %s, status %d \n",pdg, iMother,  primary->GetName(),status);	
      }	
      
      if(GetDebug() > 1 ) {
        printf("AliAnaCalorimeterQA::ClusterHistograms() - Converted Cluster mother before conversion: \n");
        printf("\t Mother label %d, pdg %d, %s, status %d, parent %d \n",iMother, pdg, primary->GetName(), status, iParent);
      }
      
    }
    
    //Overlapped pi0 (or eta, there will be very few), get the meson
    if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0) || 
       GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta)){
      if(GetDebug() > 1 ) printf("AliAnaCalorimeterQA::ClusterHistograms() - Overlapped Meson decay!, Find it: \n");
      while(pdg != 111 && pdg != 221){
        iMother = iParent;
        primary = GetMCStack()->Particle(iMother);
        status  = primary->GetStatusCode();
        iParent = primary->GetFirstMother();
        pdg     = TMath::Abs(primary->GetPdgCode());
        if(GetDebug() > 1 ) printf("\t pdg %d, %s, index %d\n",pdg,  primary->GetName(),iMother);
        if(iMother==-1) {
          printf("AliAnaCalorimeterQA::ClusterHistograms() - Tagged as Overlapped photon but meson not found, why?\n");
          //break;
        }
      }
      
      if(GetDebug() > 2 ) printf("AliAnaCalorimeterQA::ClusterHistograms() - Overlapped %s decay, label %d \n", 
                                 primary->GetName(),iMother);
    }
    
    eMC    = primary->Energy();
    ptMC   = primary->Pt();
    phiMC  = primary->Phi();
    etaMC  = primary->Eta();
    pdg    = TMath::Abs(primary->GetPdgCode());
    charge = (Int_t) TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
    
  }
  else if(GetReader()->ReadAODMCParticles() && !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCUnknown)){//it MC AOD and known tag
    //Get the list of MC particles
    if(!GetReader()->GetAODMCParticles(0)) 
      AliFatal("MCParticles not available!");
    
    aodprimary = (AliAODMCParticle*) (GetReader()->GetAODMCParticles(0))->At(label);
    iMother = label;
    pdg0    = TMath::Abs(aodprimary->GetPdgCode());
    pdg     = pdg0;
    status  = aodprimary->IsPrimary();
    vxMC    = aodprimary->Xv();
    vyMC    = aodprimary->Yv();
    iParent = aodprimary->GetMother();
    
    if(GetDebug() > 1 ) {
      printf("AliAnaCalorimeterQA::ClusterHistograms() - Cluster most contributing mother: \n");
      printf("\t Mother label %d, pdg %d, Primary? %d, Physical Primary? %d, parent %d \n",
             iMother, pdg0, aodprimary->IsPrimary(), aodprimary->IsPhysicalPrimary(), iParent);
    }
    
    //Get final particle, no conversion products
    if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)){
      if(GetDebug() > 1 ) 
        printf("AliAnaCalorimeterQA::ClusterHistograms() - Converted cluster!. Find before conversion: \n");
      //Get the parent
      aodprimary = (AliAODMCParticle*)(GetReader()->GetAODMCParticles(0))->At(iParent);
      pdg = TMath::Abs(aodprimary->GetPdgCode());
      while ((pdg == 22 || pdg == 11) && !aodprimary->IsPhysicalPrimary()) {
        iMother    = iParent;
        aodprimary = (AliAODMCParticle*)(GetReader()->GetAODMCParticles(0))->At(iMother);
        status     = aodprimary->IsPrimary();
        iParent    = aodprimary->GetMother();
        pdg        = TMath::Abs(aodprimary->GetPdgCode());
        if(GetDebug() > 1 )
          printf("\t pdg %d, index %d, Primary? %d, Physical Primary? %d \n",
                 pdg, iMother, aodprimary->IsPrimary(), aodprimary->IsPhysicalPrimary());	
      }	
      
      if(GetDebug() > 1 ) {
        printf("AliAnaCalorimeterQA::ClusterHistograms() - Converted Cluster mother before conversion: \n");
        printf("\t Mother label %d, pdg %d, parent %d, Primary? %d, Physical Primary? %d \n",
               iMother, pdg, iParent, aodprimary->IsPrimary(), aodprimary->IsPhysicalPrimary());
      }
      
    }
    
    //Overlapped pi0 (or eta, there will be very few), get the meson
    if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0) || 
       GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta)){
      if(GetDebug() > 1 ) printf("AliAnaCalorimeterQA::ClusterHistograms() - Overlapped Meson decay!, Find it: PDG %d, mom %d \n",pdg, iMother);
      while(pdg != 111 && pdg != 221){
        
        iMother    = iParent;
        aodprimary = (AliAODMCParticle*)(GetReader()->GetAODMCParticles(0))->At(iMother);
        status     = aodprimary->IsPrimary();
        iParent    = aodprimary->GetMother();
        pdg        = TMath::Abs(aodprimary->GetPdgCode());
        
        if(GetDebug() > 1 ) printf("\t pdg %d, index %d\n",pdg, iMother);
        
        if(iMother==-1) {
          printf("AliAnaCalorimeterQA::ClusterHistograms() - Tagged as Overlapped photon but meson not found, why?\n");
          //break;
        }
      }	
      
      if(GetDebug() > 2 ) printf("AliAnaCalorimeterQA::ClusterHistograms() - Overlapped %s decay, label %d \n", 
                                 aodprimary->GetName(),iMother);
    }	
    
    status = aodprimary->IsPrimary();
    eMC    = aodprimary->E();
    ptMC   = aodprimary->Pt();
    phiMC  = aodprimary->Phi();
    etaMC  = aodprimary->Eta();
    pdg    = TMath::Abs(aodprimary->GetPdgCode());
    charge = aodprimary->Charge();
    
  }
  
  //Float_t vz = primary->Vz();
  Float_t rVMC = TMath::Sqrt(vxMC*vxMC + vyMC*vyMC);
  if((pdg == 22 || TMath::Abs(pdg)==11) && status!=1) {
    fhEMVxyz   ->Fill(vxMC,vyMC);//,vz);
    fhEMR      ->Fill(e,rVMC);
  }
  
  //printf("reco e %f, pt %f, phi %f, eta %f \n", e, pt, phi, eta);
  //printf("prim e %f, pt %f, phi %f, eta %f \n", eMC,ptMC,phiMC ,etaMC );
  //printf("vertex: vx %f, vy %f, vz %f, r %f \n", vxMC, vyMC, vz, r);
  
  //Overlapped pi0 (or eta, there will be very few)
  if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0)){
    fhRecoMCE  [mcPi0][matched]     ->Fill(e,eMC);	
    if(e > 0.5 && eMC > 0.5) fhRecoMCEta[mcPi0][(matched)]->Fill(eta,etaMC);	
    if(e > 0.5 && eMC > 0.5) fhRecoMCPhi[mcPi0][(matched)]->Fill(phi,phiMC);
    if(eMC > 0) fhRecoMCRatioE  [mcPi0][(matched)]->Fill(e,e/eMC);
    fhRecoMCDeltaE  [mcPi0][(matched)]->Fill(e,eMC-e);
    fhRecoMCDeltaPhi[mcPi0][(matched)]->Fill(e,phiMC-phi);
    fhRecoMCDeltaEta[mcPi0][(matched)]->Fill(e,etaMC-eta);
  }//Overlapped pizero decay
  else     if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta)){
    fhRecoMCE  [mcEta][(matched)]     ->Fill(e,eMC);	
    if(e > 0.5 && eMC > 0.5) fhRecoMCEta[mcEta][(matched)]->Fill(eta,etaMC);	
    if(e > 0.5 && eMC > 0.5) fhRecoMCPhi[mcEta][(matched)]->Fill(phi,phiMC);
    if(eMC > 0) fhRecoMCRatioE  [mcEta][(matched)]->Fill(e,e/eMC);
    fhRecoMCDeltaE  [mcEta][(matched)]->Fill(e,eMC-e);
    fhRecoMCDeltaPhi[mcEta][(matched)]->Fill(e,phiMC-phi);
    fhRecoMCDeltaEta[mcEta][(matched)]->Fill(e,etaMC-eta);
  }//Overlapped eta decay
  else if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton)){
    fhRecoMCE  [mcPhoton][(matched)]     ->Fill(e,eMC);	
    if(e > 0.5 && eMC > 0.5) fhRecoMCEta[mcPhoton][(matched)]->Fill(eta,etaMC);	
    if(e > 0.5 && eMC > 0.5) fhRecoMCPhi[mcPhoton][(matched)]->Fill(phi,phiMC);
    if(eMC > 0) fhRecoMCRatioE  [mcPhoton][(matched)]->Fill(e,e/eMC);
    fhRecoMCDeltaE  [mcPhoton][(matched)]->Fill(e,eMC-e);
    fhRecoMCDeltaPhi[mcPhoton][(matched)]->Fill(e,phiMC-phi);
    fhRecoMCDeltaEta[mcPhoton][(matched)]->Fill(e,etaMC-eta);      
  }//photon
  else if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCElectron)){
    fhRecoMCE  [mcElectron][(matched)]     ->Fill(e,eMC);	
    if(e > 0.5 && eMC > 0.5) fhRecoMCEta[mcElectron][(matched)]->Fill(eta,etaMC);	
    if(e > 0.5 && eMC > 0.5) fhRecoMCPhi[mcElectron][(matched)]->Fill(phi,phiMC);
    if(eMC > 0) fhRecoMCRatioE  [mcElectron][(matched)]->Fill(e,e/eMC);
    fhRecoMCDeltaE  [mcElectron][(matched)]->Fill(e,eMC-e);
    fhRecoMCDeltaPhi[mcElectron][(matched)]->Fill(e,phiMC-phi);
    fhRecoMCDeltaEta[mcElectron][(matched)]->Fill(e,etaMC-eta);
    fhEMVxyz   ->Fill(vxMC,vyMC);//,vz);
    fhEMR      ->Fill(e,rVMC);
  }
  else if(charge == 0){
    fhRecoMCE  [mcNeHadron][(matched)]     ->Fill(e,eMC);	
    if(e > 0.5 && eMC > 0.5) fhRecoMCEta[mcNeHadron][(matched)]->Fill(eta,etaMC);	
    if(e > 0.5 && eMC > 0.5) fhRecoMCPhi[mcNeHadron][(matched)]->Fill(phi,phiMC);
    if(eMC > 0) fhRecoMCRatioE  [mcNeHadron][(matched)]->Fill(e,e/eMC);
    fhRecoMCDeltaE  [mcNeHadron][(matched)]->Fill(e,eMC-e);
    fhRecoMCDeltaPhi[mcNeHadron][(matched)]->Fill(e,phiMC-phi);
    fhRecoMCDeltaEta[mcNeHadron][(matched)]->Fill(e,etaMC-eta);      
    fhHaVxyz     ->Fill(vxMC,vyMC);//,vz);
    fhHaR        ->Fill(e,rVMC);
  }
  else if(charge!=0){
    fhRecoMCE  [mcChHadron][(matched)]     ->Fill(e,eMC);	
    if(e > 0.5 && eMC > 0.5) fhRecoMCEta[mcChHadron][(matched)]->Fill(eta,etaMC);	
    if(e > 0.5 && eMC > 0.5) fhRecoMCPhi[mcChHadron][(matched)]->Fill(phi,phiMC);
    if(eMC > 0) fhRecoMCRatioE  [mcChHadron][(matched)]->Fill(e,e/eMC);
    fhRecoMCDeltaE  [mcChHadron][(matched)]->Fill(e,eMC-e);
    fhRecoMCDeltaPhi[mcChHadron][(matched)]->Fill(e,phiMC-phi);
    fhRecoMCDeltaEta[mcChHadron][(matched)]->Fill(e,etaMC-eta);     
    fhHaVxyz     ->Fill(vxMC,vyMC);//,vz);
    fhHaR        ->Fill(e,rVMC);
  }
  
  if(primary || aodprimary) return kTRUE ;
  else                      return kFALSE;
  
}

//________________________________________________________________________________________________
void AliAnaCalorimeterQA::ClusterMatchedWithTrackHistograms(AliVCluster *clus, TLorentzVector mom, 
                                                            const Bool_t okPrimary, const Int_t pdg)
{
  //Histograms for clusters matched with tracks
  Float_t e   = mom.E();
  Float_t pt  = mom.Pt();
  Float_t eta = mom.Eta();
  Float_t phi = mom.Phi();
  if(phi < 0) phi +=TMath::TwoPi();
  
  if(fFillAllTH12){
    fhECharged      ->Fill(e);	
    fhPtCharged     ->Fill(pt);
    fhPhiCharged    ->Fill(phi);
    fhEtaCharged    ->Fill(eta);
  }
  
  AliVTrack * track = 0x0;
  if(!strcmp("AliESDCaloCluster",Form("%s",clus->ClassName())))
  { 
    track = dynamic_cast<AliVTrack*> (GetReader()->GetInputEvent()->GetTrack(clus->GetTrackMatchedIndex())); 
  }
  else //AOD
  {
    track = dynamic_cast<AliVTrack*> (clus->GetTrackMatched(0)); 
  }
  
  if(!track) return ;
  
  if(fFillAllTMHisto){
    if(fFillAllTH3)fhEtaPhiECharged->Fill(eta,phi,e);		
    if((fCalorimeter=="EMCAL" && GetReader()->GetEMCALPtMin() < 0.3) ||
       (fCalorimeter=="PHOS"  && GetReader()->GetPHOSPtMin()  < 0.3))   fhNCellsPerClusterMIPCharged->Fill(e, clus->GetNCells());
  }
  //printf("track index %d ntracks %d\n", esd->GetNumberOfTracks());
  
  //Study the track and matched cluster if track exists.
  if(!track) return;
  Double_t emcpos[3] = {0.,0.,0.};
  Double_t emcmom[3] = {0.,0.,0.};
  Double_t radius    = 441.0; //[cm] EMCAL radius +13cm
  Double_t bfield    = 0.;
  Double_t tphi      = 0;
  Double_t teta      = 0;
  Double_t tmom      = 0;
  Double_t tpt       = 0;
  Double_t tmom2     = 0;
  Double_t tpcSignal = 0;
  Bool_t okpos = kFALSE;
  Bool_t okmom = kFALSE;
  Bool_t okout = kFALSE;
  Int_t nITS   = 0;
  Int_t nTPC   = 0;
  
  //In case of ESDs get the parameters in this way
  if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
    if (track->GetOuterParam() ) {
      okout = kTRUE;
      
      bfield = GetReader()->GetInputEvent()->GetMagneticField();
      okpos = track->GetOuterParam()->GetXYZAt(radius,bfield,emcpos);
      okmom = track->GetOuterParam()->GetPxPyPzAt(radius,bfield,emcmom);
      if(!(okpos && okmom)) return;
      
      TVector3 position(emcpos[0],emcpos[1],emcpos[2]);
      TVector3 momentum(emcmom[0],emcmom[1],emcmom[2]);
      tphi = position.Phi();
      teta = position.Eta();
      tmom = momentum.Mag();
      
      tpt       = track->Pt();
      tmom2     = track->P();
      tpcSignal = track->GetTPCsignal();
      
      nITS = track->GetNcls(0);
      nTPC = track->GetNcls(1);
    }//Outer param available 
  }// ESDs
  else if(GetReader()->GetDataType()==AliCaloTrackReader::kAOD) {
    AliAODPid* pid = (AliAODPid*) ((AliAODTrack *) track)->GetDetPid();
    if (pid) {
      okout = kTRUE;
      pid->GetEMCALPosition(emcpos);
      pid->GetEMCALMomentum(emcmom);	
      
      TVector3 position(emcpos[0],emcpos[1],emcpos[2]);
      TVector3 momentum(emcmom[0],emcmom[1],emcmom[2]);
      tphi = position.Phi();
      teta = position.Eta();
      tmom = momentum.Mag();
      
      tpt       = track->Pt();
      tmom2     = track->P();
      tpcSignal = pid->GetTPCsignal();
      
    }//pid 
  }//AODs
  
  if(okout){
    //printf("okout\n");
    Double_t deta = teta - eta;
    Double_t dphi = tphi - phi;
    if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
    Double_t dR = sqrt(dphi*dphi + deta*deta);
    
    Double_t pOverE = tmom/e;
    
    fh1pOverE->Fill(tpt, pOverE);
    if(dR < 0.02) fh1pOverER02->Fill(tpt,pOverE);
    
    fh1dR->Fill(dR);
    fh2MatchdEdx->Fill(tmom2,tpcSignal);
    
    if(IsDataMC() && okPrimary){ 
      Double_t  charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
      
      if(TMath::Abs(pdg) == 11){
        fhMCEle1pOverE->Fill(tpt,pOverE);
        fhMCEle1dR->Fill(dR);
        fhMCEle2MatchdEdx->Fill(tmom2,tpcSignal);		
        if(dR < 0.02) fhMCEle1pOverER02->Fill(tpt,pOverE);
      }
      else if(charge!=0){
        fhMCChHad1pOverE->Fill(tpt,pOverE);
        fhMCChHad1dR->Fill(dR);
        fhMCChHad2MatchdEdx->Fill(tmom2,tpcSignal);	
        if(dR < 0.02) fhMCChHad1pOverER02->Fill(tpt,pOverE);
      }
      else if(charge == 0){
        fhMCNeutral1pOverE->Fill(tpt,pOverE);
        fhMCNeutral1dR->Fill(dR);
        fhMCNeutral2MatchdEdx->Fill(tmom2,tpcSignal);	
        if(dR < 0.02) fhMCNeutral1pOverER02->Fill(tpt,pOverE);
      }
    }//DataMC
    
    if(dR < 0.02 && pOverE > 0.5 && pOverE < 1.5
       && clus->GetNCells() > 1 && nITS > 3 && nTPC > 20) {
      fh2EledEdx->Fill(tmom2,tpcSignal);
    }
  }
  else{//no ESD external param or AODPid
    
    if(GetDebug() >= 0) printf("No ESD external param or AliAODPid \n");
    
  }//No out params  
}

//___________________________________
void AliAnaCalorimeterQA::Correlate()
{
  // Correlate information from PHOS and EMCAL and with V0 and track multiplicity
  
  //Clusters 
  TObjArray * caloClustersEMCAL = GetEMCALClusters();
  TObjArray * caloClustersPHOS  = GetPHOSClusters();
  
  Int_t nclEMCAL = caloClustersEMCAL->GetEntriesFast();
  Int_t nclPHOS  = caloClustersPHOS ->GetEntriesFast();
  
  Float_t sumClusterEnergyEMCAL = 0;
  Float_t sumClusterEnergyPHOS  = 0;
  Int_t iclus = 0;
  for(iclus = 0 ; iclus <  caloClustersEMCAL->GetEntriesFast() ; iclus++) 
    sumClusterEnergyEMCAL += ((AliVCluster*)caloClustersEMCAL->At(iclus))->E();
  for(iclus = 0 ; iclus <  caloClustersPHOS->GetEntriesFast(); iclus++) 
    sumClusterEnergyPHOS += ((AliVCluster*)caloClustersPHOS->At(iclus))->E();
  
  
  //Cells
  
  AliVCaloCells * cellsEMCAL = GetEMCALCells();
  AliVCaloCells * cellsPHOS  = GetPHOSCells();
  
  Int_t ncellsEMCAL = cellsEMCAL->GetNumberOfCells();
  Int_t ncellsPHOS  = cellsPHOS ->GetNumberOfCells();
  
  Float_t sumCellEnergyEMCAL = 0;
  Float_t sumCellEnergyPHOS  = 0;
  Int_t icell = 0;
  for(icell = 0 ; icell < cellsEMCAL->GetNumberOfCells()  ; icell++) 
    sumCellEnergyEMCAL += cellsEMCAL->GetAmplitude(icell);
  for(icell = 0 ; icell <  cellsPHOS->GetNumberOfCells(); icell++) 
    sumCellEnergyPHOS += cellsPHOS->GetAmplitude(icell);
  
  
  //Fill Histograms
  fhCaloCorrNClusters->Fill(nclEMCAL,nclPHOS);
  fhCaloCorrEClusters->Fill(sumClusterEnergyEMCAL,sumClusterEnergyPHOS);
  fhCaloCorrNCells   ->Fill(ncellsEMCAL,ncellsPHOS);
  fhCaloCorrECells   ->Fill(sumCellEnergyEMCAL,sumCellEnergyPHOS);
  
  Int_t v0S = GetV0Signal(0)+GetV0Signal(1);
  Int_t v0M = GetV0Multiplicity(0)+GetV0Multiplicity(1);
  Int_t trM = GetTrackMultiplicity();
  if(fCalorimeter=="PHOS"){
    fhCaloV0MCorrNClusters   ->Fill(v0M,nclPHOS);
    fhCaloV0MCorrEClusters   ->Fill(v0M,sumClusterEnergyPHOS);
    fhCaloV0MCorrNCells      ->Fill(v0M,ncellsPHOS);
    fhCaloV0MCorrECells      ->Fill(v0M,sumCellEnergyPHOS);
    
    fhCaloV0SCorrNClusters   ->Fill(v0S,nclPHOS);
    fhCaloV0SCorrEClusters   ->Fill(v0S,sumClusterEnergyPHOS);
    fhCaloV0SCorrNCells      ->Fill(v0S,ncellsPHOS);
    fhCaloV0SCorrECells      ->Fill(v0S,sumCellEnergyPHOS);
    
    fhCaloTrackMCorrNClusters->Fill(trM,nclPHOS);
    fhCaloTrackMCorrEClusters->Fill(trM,sumClusterEnergyPHOS);    
    fhCaloTrackMCorrNCells   ->Fill(trM,ncellsPHOS);
    fhCaloTrackMCorrECells   ->Fill(trM,sumCellEnergyPHOS);
  }
  else{
    fhCaloV0MCorrNClusters   ->Fill(v0M,nclEMCAL);
    fhCaloV0MCorrEClusters   ->Fill(v0M,sumClusterEnergyEMCAL);
    fhCaloV0MCorrNCells      ->Fill(v0M,ncellsEMCAL);
    fhCaloV0MCorrECells      ->Fill(v0M,sumCellEnergyEMCAL);
    
    fhCaloV0SCorrNClusters   ->Fill(v0S,nclEMCAL);
    fhCaloV0SCorrEClusters   ->Fill(v0S,sumClusterEnergyEMCAL);
    fhCaloV0SCorrNCells      ->Fill(v0S,ncellsEMCAL);
    fhCaloV0SCorrECells      ->Fill(v0S,sumCellEnergyEMCAL);
    
    fhCaloTrackMCorrNClusters->Fill(trM,nclEMCAL);
    fhCaloTrackMCorrEClusters->Fill(trM,sumClusterEnergyEMCAL);    
    fhCaloTrackMCorrNCells   ->Fill(trM,ncellsEMCAL);
    fhCaloTrackMCorrECells   ->Fill(trM,sumCellEnergyEMCAL);
  }
  
  if(GetDebug() > 0 )
  {
    printf("AliAnaCalorimeterQA::Correlate(): \n");
    printf("\t EMCAL: N cells %d, N clusters  %d, summed E cells %f, summed E clusters %f \n",
           ncellsEMCAL,nclEMCAL, sumCellEnergyEMCAL,sumClusterEnergyEMCAL);
    printf("\t PHOS : N cells %d, N clusters  %d, summed E cells %f, summed E clusters %f \n",
           ncellsPHOS,nclPHOS,sumCellEnergyPHOS,sumClusterEnergyPHOS);
    printf("\t V0 : Signal %d, Multiplicity  %d, Track Multiplicity %d \n", v0S,v0M,trM);
  }
  
}

//__________________________________________________
TObjString * AliAnaCalorimeterQA::GetAnalysisCuts()
{  	
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaCalorimeterQA ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Time Cut : %2.2f < T < %2.2f ns  \n",fTimeCutMin, fTimeCutMax) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"PHOS Cell Amplitude > %2.2f GeV, EMCAL Cell Amplitude > %2.2f GeV  \n",fPHOSCellAmpMin, fEMCALCellAmpMin) ;
  parList+=onePar ;
  //Get parameters set in base class.
  //parList += GetBaseParametersList() ;
  
  //Get parameters set in FiducialCut class (not available yet)
  //parlist += GetFidCut()->GetFidCutParametersList() 
	
  return new TObjString(parList) ;
}

//____________________________________________________
TList * AliAnaCalorimeterQA::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("QAHistos") ; 
  
  //Histograms
  Int_t nptbins     = GetHistoPtBins(); 	        Float_t ptmax     = GetHistoPtMax();           Float_t ptmin     = GetHistoPtMin();
  Int_t nfineptbins = GetHistoFinePtBins(); 	    Float_t ptfinemax = GetHistoFinePtMax();       Float_t ptfinemin = GetHistoFinePtMin();
  Int_t nphibins    = GetHistoPhiBins();     	    Float_t phimax    = GetHistoPhiMax();          Float_t phimin    = GetHistoPhiMin();
  Int_t netabins    = GetHistoEtaBins();          Float_t etamax    = GetHistoEtaMax();          Float_t etamin    = GetHistoEtaMin();	
  Int_t nmassbins   = GetHistoMassBins();         Float_t massmax   = GetHistoMassMax(); 	       Float_t massmin   = GetHistoMassMin();
  Int_t nasymbins   = GetHistoAsymmetryBins();    Float_t asymmax   = GetHistoAsymmetryMax();    Float_t asymmin   = GetHistoAsymmetryMin();
  Int_t nPoverEbins = GetHistoPOverEBins();       Float_t pOverEmax = GetHistoPOverEMax();       Float_t pOverEmin = GetHistoPOverEMin();
  Int_t ndedxbins   = GetHistodEdxBins();         Float_t dedxmax   = GetHistodEdxMax();         Float_t dedxmin   = GetHistodEdxMin();
  Int_t ndRbins     = GetHistodRBins();           Float_t dRmax     = GetHistodRMax();           Float_t dRmin     = GetHistodRMin();
  Int_t ntimebins   = GetHistoTimeBins();         Float_t timemax   = GetHistoTimeMax();         Float_t timemin   = GetHistoTimeMin();       
  Int_t nclbins     = GetHistoNClustersBins();    Int_t   nclmax    = GetHistoNClustersMax();    Int_t   nclmin    = GetHistoNClustersMin(); 
  Int_t ncebins     = GetHistoNCellsBins();       Int_t   ncemax    = GetHistoNCellsMax();       Int_t   ncemin    = GetHistoNCellsMin(); 
  Int_t nceclbins   = GetHistoNClusterCellBins(); Int_t   nceclmax  = GetHistoNClusterCellMax(); Int_t   nceclmin  = GetHistoNClusterCellMin(); 
  Int_t nvdistbins  = GetHistoVertexDistBins();   Float_t vdistmax  = GetHistoVertexDistMax();   Float_t vdistmin  = GetHistoVertexDistMin();
  Int_t rbins       = GetHistoRBins();            Float_t rmax      = GetHistoRMax();            Float_t rmin      = GetHistoRMin(); 
  Int_t xbins       = GetHistoXBins();            Float_t xmax      = GetHistoXMax();            Float_t xmin      = GetHistoXMin(); 
  Int_t ybins       = GetHistoYBins();            Float_t ymax      = GetHistoYMax();            Float_t ymin      = GetHistoYMin(); 
  Int_t zbins       = GetHistoZBins();            Float_t zmax      = GetHistoZMax();            Float_t zmin      = GetHistoZMin(); 
  Int_t ssbins      = GetHistoShowerShapeBins();  Float_t ssmax     = GetHistoShowerShapeMax();  Float_t ssmin     = GetHistoShowerShapeMin();
  Int_t tdbins      = GetHistoDiffTimeBins() ;    Float_t tdmax     = GetHistoDiffTimeMax();     Float_t tdmin     = GetHistoDiffTimeMin();
  
  Int_t nv0sbins    = GetHistoV0SignalBins();          Int_t nv0smax = GetHistoV0SignalMax();          Int_t nv0smin = GetHistoV0SignalMin(); 
  Int_t nv0mbins    = GetHistoV0MultiplicityBins();    Int_t nv0mmax = GetHistoV0MultiplicityMax();    Int_t nv0mmin = GetHistoV0MultiplicityMin(); 
  Int_t ntrmbins    = GetHistoTrackMultiplicityBins(); Int_t ntrmmax = GetHistoTrackMultiplicityMax(); Int_t ntrmmin = GetHistoTrackMultiplicityMin(); 
  
  //EMCAL
  fNMaxCols = 48;
  fNMaxRows = 24;
  fNRCU     = 2 ;
  //PHOS
  if(fCalorimeter=="PHOS"){
    fNMaxCols = 56;
    fNMaxRows = 64;
    fNRCU     = 4 ;
  }
  
  fhE  = new TH1F ("hE","E reconstructed clusters ", nptbins*5,ptmin,ptmax*5);  
  fhE->SetXTitle("E (GeV)");
  outputContainer->Add(fhE);
  
  if(fFillAllTH12){
    fhPt  = new TH1F ("hPt","p_{T} reconstructed clusters", nptbins,ptmin,ptmax); 
    fhPt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPt);
    
    fhPhi  = new TH1F ("hPhi","#phi reconstructed clusters ",nphibins,phimin,phimax); 
    fhPhi->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhi);
    
    fhEta  = new TH1F ("hEta","#eta reconstructed clusters ",netabins,etamin,etamax); 
    fhEta->SetXTitle("#eta ");
    outputContainer->Add(fhEta);
  }
  
  if(fFillAllTH3){
    fhEtaPhiE  = new TH3F ("hEtaPhiE","#eta vs #phi vs energy, reconstructed clusters",
                           netabins,etamin,etamax,nphibins,phimin,phimax,nptbins,ptmin,ptmax); 
    fhEtaPhiE->SetXTitle("#eta ");
    fhEtaPhiE->SetYTitle("#phi (rad)");
    fhEtaPhiE->SetZTitle("E (GeV) ");
    outputContainer->Add(fhEtaPhiE);
  }
  
  fhClusterTimeEnergy  = new TH2F ("hClusterTimeEnergy","energy vs TOF, reconstructed clusters",
                                   nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
  fhClusterTimeEnergy->SetXTitle("E (GeV) ");
  fhClusterTimeEnergy->SetYTitle("TOF (ns)");
  outputContainer->Add(fhClusterTimeEnergy);
  
  fhClusterPairDiffTimeE = new TH2F("hClusterPairDiffTimeE","cluster pair time difference vs E, only good clusters",
                                    nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
  fhClusterPairDiffTimeE->SetXTitle("E_{cluster} (GeV)");
  fhClusterPairDiffTimeE->SetYTitle("#Delta t (ns)");
  outputContainer->Add(fhClusterPairDiffTimeE);  
  
  fhLambda0  = new TH2F ("hLambda0","shower shape, #lambda^{2}_{0} vs E",
                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
  fhLambda0->SetXTitle("E_{cluster}");
  fhLambda0->SetYTitle("#lambda^{2}_{0}");
  outputContainer->Add(fhLambda0); 
  
  fhLambda1  = new TH2F ("hLambda1","shower shape, #lambda^{2}_{1} vs E for bad cluster ",
                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
  fhLambda1->SetXTitle("E_{cluster}");
  fhLambda1->SetYTitle("#lambda^{2}_{1}");
  outputContainer->Add(fhLambda1); 
  
  fhDispersion  = new TH2F ("hDispersion","shower shape, Dispersion^{2} vs E for bad cluster ",
                            nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
  fhDispersion->SetXTitle("E_{cluster}");
  fhDispersion->SetYTitle("Dispersion");
  outputContainer->Add(fhDispersion);       
  
  fhClusterMaxCellCloseCellRatio  = new TH2F ("hClusterMaxCellCloseCellRatio","energy vs ratio of max cell / neighbour cell, reconstructed clusters",
                                              nptbins,ptmin,ptmax, 100,0,1.); 
  fhClusterMaxCellCloseCellRatio->SetXTitle("E_{cluster} (GeV) ");
  fhClusterMaxCellCloseCellRatio->SetYTitle("E_{cell i}/E_{cell max}");
  outputContainer->Add(fhClusterMaxCellCloseCellRatio);
  
  fhClusterMaxCellCloseCellDiff  = new TH2F ("hClusterMaxCellCloseCellDiff","energy vs ratio of max cell / neighbour cell, reconstructed clusters",
                                             nptbins,ptmin,ptmax, 500,0,100.); 
  fhClusterMaxCellCloseCellDiff->SetXTitle("E_{cluster} (GeV) ");
  fhClusterMaxCellCloseCellDiff->SetYTitle("E_{cell max}-E_{cell i} (GeV)");
  outputContainer->Add(fhClusterMaxCellCloseCellDiff);
  
  fhClusterMaxCellDiff  = new TH2F ("hClusterMaxCellDiff","energy vs difference of cluster energy - max cell energy / cluster energy, good clusters",
                                    nptbins,ptmin,ptmax, 500,0,1.); 
  fhClusterMaxCellDiff->SetXTitle("E_{cluster} (GeV) ");
  fhClusterMaxCellDiff->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  outputContainer->Add(fhClusterMaxCellDiff);  
  
  fhClusterMaxCellDiffNoCut  = new TH2F ("hClusterMaxCellDiffNoCut","energy vs difference of cluster energy - max cell energy / cluster energy",
                                         nptbins,ptmin,ptmax, 500,0,1.); 
  fhClusterMaxCellDiffNoCut->SetXTitle("E_{cluster} (GeV) ");
  fhClusterMaxCellDiffNoCut->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  outputContainer->Add(fhClusterMaxCellDiffNoCut);  
  
  fhLambda0vsClusterMaxCellDiffE0  = new TH2F ("hLambda0vsClusterMaxCellDiffE0","shower shape, #lambda^{2}_{0} vs fraction of energy carried by max cell, E < 2 GeV ",
                                               ssbins,ssmin,ssmax,500,0,1.); 
  fhLambda0vsClusterMaxCellDiffE0->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  fhLambda0vsClusterMaxCellDiffE0->SetXTitle("#lambda^{2}_{0}");
  outputContainer->Add(fhLambda0vsClusterMaxCellDiffE0); 
  
  fhLambda0vsClusterMaxCellDiffE2  = new TH2F ("hLambda0vsClusterMaxCellDiffE2","shower shape, #lambda^{2}_{0} vs fraction of energy carried by max cell, 2 < E < 6 GeV ",
                                               ssbins,ssmin,ssmax,500,0,1.); 
  fhLambda0vsClusterMaxCellDiffE2->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  fhLambda0vsClusterMaxCellDiffE2->SetXTitle("#lambda^{2}_{0}");
  outputContainer->Add(fhLambda0vsClusterMaxCellDiffE2); 
  
  fhLambda0vsClusterMaxCellDiffE6  = new TH2F ("hLambda0vsClusterMaxCellDiffE6","shower shape, #lambda^{2}_{0} vs fraction of energy carried by max cell, E > 6 ",
                                               ssbins,ssmin,ssmax,500,0,1.); 
  fhLambda0vsClusterMaxCellDiffE6->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  fhLambda0vsClusterMaxCellDiffE6->SetXTitle("#lambda^{2}_{0}");
  outputContainer->Add(fhLambda0vsClusterMaxCellDiffE6); 
  
  fhNCellsvsClusterMaxCellDiffE0  = new TH2F ("hNCellsvsClusterMaxCellDiffE0","N cells per cluster vs fraction of energy carried by max cell, E < 2 GeV ",
                                              nceclbins,nceclmin,nceclmax,500,0,1.); 
  fhNCellsvsClusterMaxCellDiffE0->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  fhNCellsvsClusterMaxCellDiffE0->SetXTitle("N cells per cluster");
  outputContainer->Add(fhNCellsvsClusterMaxCellDiffE0); 
  
  fhNCellsvsClusterMaxCellDiffE2  = new TH2F ("hNCellsvsClusterMaxCellDiffE2","N cells per cluster vs fraction of energy carried by max cell, 2 < E < 6 GeV ",
                                              nceclbins,nceclmin,nceclmax,500,0,1.); 
  fhNCellsvsClusterMaxCellDiffE2->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  fhNCellsvsClusterMaxCellDiffE2->SetXTitle("N cells per cluster");
  outputContainer->Add(fhNCellsvsClusterMaxCellDiffE2); 
  
  fhNCellsvsClusterMaxCellDiffE6  = new TH2F ("hNCellsvsClusterMaxCellDiffE6","N cells per cluster vs fraction of energy carried by max cell, E > 6 ",
                                              nceclbins,nceclmin,nceclmax,500,0,1.); 
  fhNCellsvsClusterMaxCellDiffE6->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  fhNCellsvsClusterMaxCellDiffE6->SetXTitle("N cells per cluster");
  outputContainer->Add(fhNCellsvsClusterMaxCellDiffE6); 
  
  if(fCalorimeter=="EMCAL" && !GetCaloUtils()->GetEMCALRecoUtils()->IsRejectExoticCluster() && fStudyBadClusters){
    
    fhBadClusterEnergy  = new TH1F ("hBadClusterEnergy","Bad cluster energy", nptbins,ptmin,ptmax); 
    fhBadClusterEnergy->SetXTitle("E_{cluster} (GeV) ");
    outputContainer->Add(fhBadClusterEnergy);
    
    fhBadClusterMaxCellCloseCellRatio  = new TH2F ("hBadClusterMaxCellCloseCellRatio","energy vs ratio of max cell / neighbour cell constributing cell, reconstructed bad clusters",
                                                   nptbins,ptmin,ptmax, 100,0,1.); 
    fhBadClusterMaxCellCloseCellRatio->SetXTitle("E_{cluster} (GeV) ");
    fhBadClusterMaxCellCloseCellRatio->SetYTitle("ratio");
    outputContainer->Add(fhBadClusterMaxCellCloseCellRatio);
    
    fhBadClusterMaxCellCloseCellDiff  = new TH2F ("hBadClusterMaxCellCloseCellDiff","energy vs ratio of max cell - neighbour cell constributing cell, reconstructed bad clusters",
                                                  nptbins,ptmin,ptmax, 500,0,100); 
    fhBadClusterMaxCellCloseCellDiff->SetXTitle("E_{cluster} (GeV) ");
    fhBadClusterMaxCellCloseCellDiff->SetYTitle("E_{cell max} - E_{cell i} (GeV)");
    outputContainer->Add(fhBadClusterMaxCellCloseCellDiff);    
    
    fhBadClusterMaxCellDiff  = new TH2F ("hBadClusterMaxCellDiff","energy vs difference of cluster energy - max cell energy / cluster energy for bad clusters",
                                         nptbins,ptmin,ptmax, 500,0,1.); 
    fhBadClusterMaxCellDiff->SetXTitle("E_{cluster} (GeV) ");
    fhBadClusterMaxCellDiff->SetYTitle("(E_{cluster} - E_{cell max}) / E_{cluster}");
    outputContainer->Add(fhBadClusterMaxCellDiff);
    
    fhBadClusterTimeEnergy  = new TH2F ("hBadClusterTimeEnergy","energy vs TOF of reconstructed bad clusters",
                                        nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
    fhBadClusterTimeEnergy->SetXTitle("E_{cluster} (GeV) ");
    fhBadClusterTimeEnergy->SetYTitle("TOF (ns)");
    outputContainer->Add(fhBadClusterTimeEnergy);    
    
    fhBadClusterPairDiffTimeE = new TH2F("hBadClusterPairDiffTimeE","cluster pair time difference (bad - good) vs E from bad cluster",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
    fhBadClusterPairDiffTimeE->SetXTitle("E_{bad cluster} (GeV)");
    fhBadClusterPairDiffTimeE->SetYTitle("#Delta t (ns)");
    outputContainer->Add(fhBadClusterPairDiffTimeE);    
    
    if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
      fhBadCellTimeSpreadRespectToCellMax = new TH2F ("hBadCellTimeSpreadRespectToCellMax","t_{cell max}-t_{cell i} from bad cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax); 
      fhBadCellTimeSpreadRespectToCellMax->SetXTitle("E (GeV)");
      fhBadCellTimeSpreadRespectToCellMax->SetYTitle("#Delta t_{cell max - i} (ns)");
      outputContainer->Add(fhBadCellTimeSpreadRespectToCellMax);
      
      fhBadClusterMaxCellDiffAverageTime = new TH2F ("hBadClusterMaxCellDiffAverageTime","t_{cell max}-t_{average} from bad cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax); 
      fhBadClusterMaxCellDiffAverageTime->SetXTitle("E (GeV)");
      fhBadClusterMaxCellDiffAverageTime->SetYTitle("#Delta t_{cell max - average} (ns)");
      outputContainer->Add(fhBadClusterMaxCellDiffAverageTime);
      
      fhBadClusterMaxCellDiffAverageNoMaxTime = new TH2F ("hBadClusterMaxCellDiffAverageNoMaxTime","t_{cell max}-t_{average} from bad cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax); 
      fhBadClusterMaxCellDiffAverageNoMaxTime->SetXTitle("E (GeV)");
      fhBadClusterMaxCellDiffAverageNoMaxTime->SetYTitle("#Delta t_{cell max - average} (ns)");
      outputContainer->Add(fhBadClusterMaxCellDiffAverageNoMaxTime);
      
      fhBadClusterMaxCellDiffWeightedTime = new TH2F ("hBadClusterMaxCellDiffWeightedTime","t_{cell max}-t_{weighted} from bad cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax); 
      fhBadClusterMaxCellDiffWeightedTime->SetXTitle("E (GeV)");
      fhBadClusterMaxCellDiffWeightedTime->SetYTitle("#Delta t_{cell max - weighted} (ns)");
      outputContainer->Add(fhBadClusterMaxCellDiffWeightedTime);
      
      fhBadClusterMaxCellDiffWeightedNoMaxTime = new TH2F ("hBadClusterMaxCellDiffWeightedNoMaxTime","t_{cell max}-t_{average} from bad cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax); 
      fhBadClusterMaxCellDiffWeightedNoMaxTime->SetXTitle("E (GeV)");
      fhBadClusterMaxCellDiffWeightedNoMaxTime->SetYTitle("#Delta t_{cell max - weighted} (ns)");
      outputContainer->Add(fhBadClusterMaxCellDiffWeightedNoMaxTime);      
      
    }  
    
  }
  
  // Cluster size in terms of cells
  if(fStudyClustersAsymmetry){
    fhDeltaIEtaDeltaIPhiE0[0]  = new TH2F ("hDeltaIEtaDeltaIPhiE0"," Cluster size in columns vs rows for E < 2 GeV, n cells > 3",
                                           50,0,50,50,0,50); 
    fhDeltaIEtaDeltaIPhiE0[0]->SetXTitle("#Delta Column");
    fhDeltaIEtaDeltaIPhiE0[0]->SetYTitle("#Delta Row");
    outputContainer->Add(fhDeltaIEtaDeltaIPhiE0[0]); 
    
    fhDeltaIEtaDeltaIPhiE2[0]  = new TH2F ("hDeltaIEtaDeltaIPhiE2"," Cluster size in columns vs rows for 2 <E < 6 GeV, n cells > 3",
                                           50,0,50,50,0,50); 
    fhDeltaIEtaDeltaIPhiE2[0]->SetXTitle("#Delta Column");
    fhDeltaIEtaDeltaIPhiE2[0]->SetYTitle("#Delta Row");
    outputContainer->Add(fhDeltaIEtaDeltaIPhiE2[0]); 
    
    fhDeltaIEtaDeltaIPhiE6[0]  = new TH2F ("hDeltaIEtaDeltaIPhiE6"," Cluster size in columns vs rows for E > 6 GeV, n cells > 3",
                                           50,0,50,50,0,50); 
    fhDeltaIEtaDeltaIPhiE6[0]->SetXTitle("#Delta Column");
    fhDeltaIEtaDeltaIPhiE6[0]->SetYTitle("#Delta Row");
    outputContainer->Add(fhDeltaIEtaDeltaIPhiE6[0]); 
    
    fhDeltaIA[0]  = new TH2F ("hDeltaIA"," Cluster *asymmetry* in cell units vs E",
                              nptbins,ptmin,ptmax,21,-1.05,1.05); 
    fhDeltaIA[0]->SetXTitle("E_{cluster}");
    fhDeltaIA[0]->SetYTitle("A_{cell in cluster}");
    outputContainer->Add(fhDeltaIA[0]); 
    
    fhDeltaIAL0[0]  = new TH2F ("hDeltaIAL0"," Cluster *asymmetry* in cell units vs #lambda^{2}_{0}",
                                ssbins,ssmin,ssmax,21,-1.05,1.05); 
    fhDeltaIAL0[0]->SetXTitle("#lambda^{2}_{0}");
    fhDeltaIAL0[0]->SetYTitle("A_{cell in cluster}");
    outputContainer->Add(fhDeltaIAL0[0]); 
    
    fhDeltaIAL1[0]  = new TH2F ("hDeltaIAL1"," Cluster *asymmetry* in cell units vs #lambda^{2}_{1}",
                                ssbins,ssmin,ssmax,21,-1.05,1.05); 
    fhDeltaIAL1[0]->SetXTitle("#lambda^{2}_{1}");
    fhDeltaIAL1[0]->SetYTitle("A_{cell in cluster}");
    outputContainer->Add(fhDeltaIAL1[0]); 
    
    fhDeltaIANCells[0]  = new TH2F ("hDeltaIANCells"," Cluster *asymmetry* in cell units vs N cells in cluster",
                                    nceclbins,nceclmin,nceclmax,21,-1.05,1.05); 
    fhDeltaIANCells[0]->SetXTitle("N_{cell in cluster}");
    fhDeltaIANCells[0]->SetYTitle("A_{cell in cluster}");
    outputContainer->Add(fhDeltaIANCells[0]); 
    
    
    fhDeltaIEtaDeltaIPhiE0[1]  = new TH2F ("hDeltaIEtaDeltaIPhiE0Charged"," Cluster size in columns vs rows for E < 2 GeV, n cells > 3, matched with track",
                                           50,0,50,50,0,50); 
    fhDeltaIEtaDeltaIPhiE0[1]->SetXTitle("#Delta Column");
    fhDeltaIEtaDeltaIPhiE0[1]->SetYTitle("#Delta Row");
    outputContainer->Add(fhDeltaIEtaDeltaIPhiE0[1]); 
    
    fhDeltaIEtaDeltaIPhiE2[1]  = new TH2F ("hDeltaIEtaDeltaIPhiE2Charged"," Cluster size in columns vs rows for 2 <E < 6 GeV, n cells > 3, matched with track",
                                           50,0,50,50,0,50); 
    fhDeltaIEtaDeltaIPhiE2[1]->SetXTitle("#Delta Column");
    fhDeltaIEtaDeltaIPhiE2[1]->SetYTitle("#Delta Row");
    outputContainer->Add(fhDeltaIEtaDeltaIPhiE2[1]); 
    
    fhDeltaIEtaDeltaIPhiE6[1]  = new TH2F ("hDeltaIEtaDeltaIPhiE6Charged"," Cluster size in columns vs rows for E > 6 GeV, n cells > 3, matched with track",
                                           50,0,50,50,0,50); 
    fhDeltaIEtaDeltaIPhiE6[1]->SetXTitle("#Delta Column");
    fhDeltaIEtaDeltaIPhiE6[1]->SetYTitle("#Delta Row");
    outputContainer->Add(fhDeltaIEtaDeltaIPhiE6[1]); 
    
    fhDeltaIA[1]  = new TH2F ("hDeltaIACharged"," Cluster *asymmetry* in cell units vs E, matched with track",
                              nptbins,ptmin,ptmax,21,-1.05,1.05); 
    fhDeltaIA[1]->SetXTitle("E_{cluster}");
    fhDeltaIA[1]->SetYTitle("A_{cell in cluster}");
    outputContainer->Add(fhDeltaIA[1]); 
    
    fhDeltaIAL0[1]  = new TH2F ("hDeltaIAL0Charged"," Cluster *asymmetry* in cell units vs #lambda^{2}_{0}, matched with track",
                                ssbins,ssmin,ssmax,21,-1.05,1.05); 
    fhDeltaIAL0[1]->SetXTitle("#lambda^{2}_{0}");
    fhDeltaIAL0[1]->SetYTitle("A_{cell in cluster}");
    outputContainer->Add(fhDeltaIAL0[1]); 
    
    fhDeltaIAL1[1]  = new TH2F ("hDeltaIAL1Charged"," Cluster *asymmetry* in cell units vs #lambda^{2}_{1}, matched with track",
                                ssbins,ssmin,ssmax,21,-1.05,1.05); 
    fhDeltaIAL1[1]->SetXTitle("#lambda^{2}_{1}");
    fhDeltaIAL1[1]->SetYTitle("A_{cell in cluster}");
    outputContainer->Add(fhDeltaIAL1[1]); 
    
    fhDeltaIANCells[1]  = new TH2F ("hDeltaIANCellsCharged"," Cluster *asymmetry* in cell units vs N cells in cluster, matched with track",
                                    nceclbins,nceclmin,nceclmax,21,-1.05,1.05); 
    fhDeltaIANCells[1]->SetXTitle("N_{cell in cluster}");
    fhDeltaIANCells[1]->SetYTitle("A_{cell in cluster}");
    outputContainer->Add(fhDeltaIANCells[1]); 
    
    if(IsDataMC()){
      TString particle[]={"Photon","Electron","Conversion","Hadron"};
      for (Int_t iPart = 0; iPart < 4; iPart++) {
        
        fhDeltaIAMC[iPart]  = new TH2F (Form("hDeltaIA_MC%s",particle[iPart].Data()),Form(" Cluster *asymmetry* in cell units vs E, from %s",particle[iPart].Data()),
                                        nptbins,ptmin,ptmax,21,-1.05,1.05); 
        fhDeltaIAMC[iPart]->SetXTitle("E_{cluster}");
        fhDeltaIAMC[iPart]->SetYTitle("A_{cell in cluster}");
        outputContainer->Add(fhDeltaIAMC[iPart]);     
      }
    }
  }
  
  if(fStudyWeight){
    
    fhECellClusterRatio  = new TH2F ("hECellClusterRatio"," cell energy / cluster energy vs cluster energy",
                                     nptbins,ptmin,ptmax, 100,0,1.); 
    fhECellClusterRatio->SetXTitle("E_{cluster} (GeV) ");
    fhECellClusterRatio->SetYTitle("E_{cell i}/E_{cluster}");
    outputContainer->Add(fhECellClusterRatio);
    
    fhECellClusterLogRatio  = new TH2F ("hECellClusterLogRatio"," Log(cell energy / cluster energy) vs cluster energy",
                                        nptbins,ptmin,ptmax, 100,-10,10); 
    fhECellClusterLogRatio->SetXTitle("E_{cluster} (GeV) ");
    fhECellClusterLogRatio->SetYTitle("E_{cell i}/E_{cluster}");
    outputContainer->Add(fhECellClusterLogRatio);
    
    fhEMaxCellClusterRatio  = new TH2F ("hEMaxCellClusterRatio"," max cell energy / cluster energy vs cluster energy",
                                        nptbins,ptmin,ptmax, 100,0,1.); 
    fhEMaxCellClusterRatio->SetXTitle("E_{cluster} (GeV) ");
    fhEMaxCellClusterRatio->SetYTitle("E_{max cell}/E_{cluster}");
    outputContainer->Add(fhEMaxCellClusterRatio);
    
    fhEMaxCellClusterLogRatio  = new TH2F ("hEMaxCellClusterLogRatio"," Log(max cell energy / cluster energy) vs cluster energy",
                                           nptbins,ptmin,ptmax, 100,-10,10); 
    fhEMaxCellClusterLogRatio->SetXTitle("E_{cluster} (GeV) ");
    fhEMaxCellClusterLogRatio->SetYTitle("E_{max cell}/E_{cluster}");
    outputContainer->Add(fhEMaxCellClusterLogRatio);
    
    for(Int_t iw = 0; iw < 7; iw++){
      fhLambda0ForW0[iw]  = new TH2F (Form("hLambda0ForW0%d",iw),Form("shower shape, #lambda^{2}_{0} vs E, w0 = %1.1f",3+0.5*iw),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLambda0ForW0[iw]->SetXTitle("E_{cluster}");
      fhLambda0ForW0[iw]->SetYTitle("#lambda^{2}_{0}");
      outputContainer->Add(fhLambda0ForW0[iw]); 
      
      fhLambda1ForW0[iw]  = new TH2F (Form("hLambda1ForW0%d",iw),Form("shower shape, #lambda^{2}_{1} vs E, w0 = %1.1f",3+0.5*iw),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLambda1ForW0[iw]->SetXTitle("E_{cluster}");
      fhLambda1ForW0[iw]->SetYTitle("#lambda^{2}_{1}");
      outputContainer->Add(fhLambda1ForW0[iw]); 
      
      if(IsDataMC()){
        TString mcnames[] = {"Photon", "Electron","Conversion","Pi0","Hadron"};
        for(Int_t imc = 0; imc < 5; imc++){
          fhLambda0ForW0MC[iw][imc]  = new TH2F (Form("hLambda0ForW0%d_MC%s",iw,mcnames[imc].Data()),
                                                 Form("shower shape, #lambda^{2}_{0} vs E, w0 = %1.1f, for MC %s",3+0.5*iw,mcnames[imc].Data()),
                                                 nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhLambda0ForW0MC[iw][imc]->SetXTitle("E_{cluster}");
          fhLambda0ForW0MC[iw][imc]->SetYTitle("#lambda^{2}_{0}");
          outputContainer->Add(fhLambda0ForW0MC[iw][imc]); 
          
          fhLambda1ForW0MC[iw][imc]  = new TH2F (Form("hLambda1ForW0%d_MC%s",iw,mcnames[imc].Data()),
                                                 Form("shower shape, #lambda^{2}_{1} vs E, w0 = %1.1f, for MC %s",3+0.5*iw,mcnames[imc].Data()),
                                                 nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhLambda1ForW0MC[iw][imc]->SetXTitle("E_{cluster}");
          fhLambda1ForW0MC[iw][imc]->SetYTitle("#lambda^{2}_{1}");
          outputContainer->Add(fhLambda1ForW0MC[iw][imc]); 
        }
      }
      
    } 
    
  }
  
  //Track Matching
  if(fFillAllTMHisto){
    if(fFillAllTH12){
      fhECharged  = new TH1F ("hECharged","E reconstructed clusters, matched with track", nptbins,ptmin,ptmax); 
      fhECharged->SetXTitle("E (GeV)");
      outputContainer->Add(fhECharged);
      
      fhPtCharged  = new TH1F ("hPtCharged","p_{T} reconstructed clusters, matched with track", nptbins,ptmin,ptmax); 
      fhPtCharged->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtCharged);
      
      fhPhiCharged  = new TH1F ("hPhiCharged","#phi reconstructed clusters, matched with track",nphibins,phimin,phimax); 
      fhPhiCharged->SetXTitle("#phi (rad)");
      outputContainer->Add(fhPhiCharged);
      
      fhEtaCharged  = new TH1F ("hEtaCharged","#eta reconstructed clusters, matched with track",netabins,etamin,etamax); 
      fhEtaCharged->SetXTitle("#eta ");
      outputContainer->Add(fhEtaCharged);
    }
    if(fFillAllTH3){
      fhEtaPhiECharged  = new TH3F ("hEtaPhiECharged","#eta vs #phi, reconstructed clusters, matched with track",
                                    netabins,etamin,etamax,nphibins,phimin,phimax,nptbins,ptmin,ptmax); 
      fhEtaPhiECharged->SetXTitle("#eta ");
      fhEtaPhiECharged->SetYTitle("#phi ");
      fhEtaPhiECharged->SetZTitle("E (GeV) ");
      outputContainer->Add(fhEtaPhiECharged);	
    }
    
    fh1pOverE = new TH2F("h1pOverE","TRACK matches p/E",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fh1pOverE->SetYTitle("p/E");
    fh1pOverE->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fh1pOverE);
    
    fh1dR = new TH1F("h1dR","TRACK matches dR",ndRbins,dRmin,dRmax);
    fh1dR->SetXTitle("#Delta R (rad)");
    outputContainer->Add(fh1dR) ;
    
    fh2MatchdEdx = new TH2F("h2MatchdEdx","dE/dx vs. p for all matches",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fh2MatchdEdx->SetXTitle("p (GeV/c)");
    fh2MatchdEdx->SetYTitle("<dE/dx>");
    outputContainer->Add(fh2MatchdEdx);
    
    fh2EledEdx = new TH2F("h2EledEdx","dE/dx vs. p for electrons",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fh2EledEdx->SetXTitle("p (GeV/c)");
    fh2EledEdx->SetYTitle("<dE/dx>");
    outputContainer->Add(fh2EledEdx) ;
    
    fh1pOverER02 = new TH2F("h1pOverER02","TRACK matches p/E, all",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fh1pOverER02->SetYTitle("p/E");
    fh1pOverER02->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fh1pOverER02);	
  }
  
  if(fFillAllPi0Histo){
    fhIM  = new TH2F ("hIM","Cluster pairs Invariant mass vs reconstructed pair energy, ncell > 1",nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
    fhIM->SetXTitle("p_{T, cluster pairs} (GeV) ");
    fhIM->SetYTitle("M_{cluster pairs} (GeV/c^{2})");
    outputContainer->Add(fhIM);
    
    fhAsym  = new TH2F ("hAssym","Cluster pairs Asymmetry vs reconstructed pair energy",nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax); 
    fhAsym->SetXTitle("p_{T, cluster pairs} (GeV) ");
    fhAsym->SetYTitle("Asymmetry");
    outputContainer->Add(fhAsym);	
    
  }
  
  fhNCellsPerClusterNoCut  = new TH2F ("hNCellsPerClusterNoCut","# cells per cluster vs energy vs #eta, no bad clusters cut",
                                       nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
  fhNCellsPerClusterNoCut->SetXTitle("E (GeV)");
  fhNCellsPerClusterNoCut->SetYTitle("n cells");
  outputContainer->Add(fhNCellsPerClusterNoCut);
  
  fhNCellsPerCluster  = new TH2F ("hNCellsPerCluster","# cells per cluster vs energy vs #eta",nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
  fhNCellsPerCluster->SetXTitle("E (GeV)");
  fhNCellsPerCluster->SetYTitle("n cells");
  outputContainer->Add(fhNCellsPerCluster);
  
  if((fCalorimeter=="EMCAL" && GetReader()->GetEMCALPtMin() < 0.3) ||
     (fCalorimeter=="PHOS"  && GetReader()->GetPHOSPtMin()  < 0.3)) {
    fhNCellsPerClusterMIP  = new TH2F ("hNCellsPerClusterMIP","# cells per cluster vs energy vs #eta, smaller bin for MIP search", 
                                       40,0.,2., 11,0,10); 
    fhNCellsPerClusterMIP->SetXTitle("E (GeV)");
    fhNCellsPerClusterMIP->SetYTitle("n cells");
    outputContainer->Add(fhNCellsPerClusterMIP);
    
    
    if(fFillAllTMHisto){
      fhNCellsPerClusterMIPCharged  = new TH2F ("hNCellsPerClusterMIPCharged","# cells per track-matched cluster vs energy vs #eta, smaller bin for MIP search", 
                                                40,0.,2., 11,0,10); 
      fhNCellsPerClusterMIPCharged->SetXTitle("E (GeV)");
      fhNCellsPerClusterMIPCharged->SetYTitle("n cells");
      outputContainer->Add(fhNCellsPerClusterMIPCharged);
    }
	}
  
  fhNClusters  = new TH1F ("hNClusters","# clusters", nclbins,nclmin,nclmax); 
  fhNClusters->SetXTitle("number of clusters");
  outputContainer->Add(fhNClusters);
  
  if(fFillAllPosHisto2){
    
    if(fFillAllTH3){
      fhXYZ  = new TH3F ("hXYZ","Cluster: x vs y vs z",xbins,xmin,xmax,ybins,ymin,ymax,zbins,zmin,zmax); 
      fhXYZ->SetXTitle("x (cm)");
      fhXYZ->SetYTitle("y (cm)");
      fhXYZ->SetZTitle("z (cm) ");
      outputContainer->Add(fhXYZ);  
    }
    
    fhXNCells  = new TH2F ("hXNCells","Cluster X position vs N Cells per Cluster",xbins,xmin,xmax,nceclbins,nceclmin,nceclmax); 
    fhXNCells->SetXTitle("x (cm)");
    fhXNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhXNCells);
    
    fhZNCells  = new TH2F ("hZNCells","Cluster Z position vs N Cells per Cluster",zbins,zmin,zmax,nceclbins,nceclmin,nceclmax); 
    fhZNCells->SetXTitle("z (cm)");
    fhZNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhZNCells);
    
    fhXE  = new TH2F ("hXE","Cluster X position vs cluster energy",xbins,xmin,xmax,nptbins,ptmin,ptmax); 
    fhXE->SetXTitle("x (cm)");
    fhXE->SetYTitle("E (GeV)");
    outputContainer->Add(fhXE);
    
    fhZE  = new TH2F ("hZE","Cluster Z position vs cluster energy",zbins,zmin,zmax,nptbins,ptmin,ptmax); 
    fhZE->SetXTitle("z (cm)");
    fhZE->SetYTitle("E (GeV)");
    outputContainer->Add(fhZE);    
    
    fhRNCells  = new TH2F ("hRNCells","Cluster R position vs N Cells per Cluster",rbins,rmin,rmax,nceclbins,nceclmin,nceclmax); 
    fhRNCells->SetXTitle("r = #sqrt{x^{2}+y^{2}} (cm)");
    fhRNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhRNCells);
    
    
    fhYNCells  = new TH2F ("hYNCells","Cluster Y position vs N Cells per Cluster",ybins,ymin,ymax,nceclbins,nceclmin,nceclmax); 
    fhYNCells->SetXTitle("y (cm)");
    fhYNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhYNCells);
    
    fhRE  = new TH2F ("hRE","Cluster R position vs cluster energy",rbins,rmin,rmax,nptbins,ptmin,ptmax); 
    fhRE->SetXTitle("r = #sqrt{x^{2}+y^{2}} (cm)");
    fhRE->SetYTitle("E (GeV)");
    outputContainer->Add(fhRE);
    
    fhYE  = new TH2F ("hYE","Cluster Y position vs cluster energy",ybins,ymin,ymax,nptbins,ptmin,ptmax); 
    fhYE->SetXTitle("y (cm)");
    fhYE->SetYTitle("E (GeV)");
    outputContainer->Add(fhYE);
  }
  if(fFillAllPosHisto){
    
    fhRCellE  = new TH2F ("hRCellE","Cell R position vs cell energy",rbins,rmin,rmax,nptbins,ptmin,ptmax); 
    fhRCellE->SetXTitle("r = #sqrt{x^{2}+y^{2}} (cm)");
    fhRCellE->SetYTitle("E (GeV)");
    outputContainer->Add(fhRCellE);
    
    fhXCellE  = new TH2F ("hXCellE","Cell X position vs cell energy",xbins,xmin,xmax,nptbins,ptmin,ptmax); 
    fhXCellE->SetXTitle("x (cm)");
    fhXCellE->SetYTitle("E (GeV)");
    outputContainer->Add(fhXCellE);
    
    fhYCellE  = new TH2F ("hYCellE","Cell Y position vs cell energy",ybins,ymin,ymax,nptbins,ptmin,ptmax); 
    fhYCellE->SetXTitle("y (cm)");
    fhYCellE->SetYTitle("E (GeV)");
    outputContainer->Add(fhYCellE);
    
    fhZCellE  = new TH2F ("hZCellE","Cell Z position vs cell energy",zbins,zmin,zmax,nptbins,ptmin,ptmax); 
    fhZCellE->SetXTitle("z (cm)");
    fhZCellE->SetYTitle("E (GeV)");
    outputContainer->Add(fhZCellE);
    
    fhXYZCell  = new TH3F ("hXYZCell","Cell : x vs y vs z",xbins,xmin,xmax,ybins,ymin,ymax,zbins,zmin,zmax); 
    fhXYZCell->SetXTitle("x (cm)");
    fhXYZCell->SetYTitle("y (cm)");
    fhXYZCell->SetZTitle("z (cm)");
    outputContainer->Add(fhXYZCell);
    
    
    Float_t dx = TMath::Abs(xmin)+TMath::Abs(xmax);
    Float_t dy = TMath::Abs(ymin)+TMath::Abs(ymax);
    Float_t dz = TMath::Abs(zmin)+TMath::Abs(zmax);
    Float_t dr = TMath::Abs(rmin)+TMath::Abs(rmax);
    
    fhDeltaCellClusterRNCells  = new TH2F ("hDeltaCellClusterRNCells","Cluster-Cell R position vs N Cells per Cluster",rbins*2,-dr,dr,nceclbins,nceclmin,nceclmax); 
    fhDeltaCellClusterRNCells->SetXTitle("r = #sqrt{x^{2}+y^{2}} (cm)");
    fhDeltaCellClusterRNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhDeltaCellClusterRNCells);
    
    fhDeltaCellClusterXNCells  = new TH2F ("hDeltaCellClusterXNCells","Cluster-Cell X position vs N Cells per Cluster",xbins*2,-dx,dx,nceclbins,nceclmin,nceclmax); 
    fhDeltaCellClusterXNCells->SetXTitle("x (cm)");
    fhDeltaCellClusterXNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhDeltaCellClusterXNCells);
    
    fhDeltaCellClusterYNCells  = new TH2F ("hDeltaCellClusterYNCells","Cluster-Cell Y position vs N Cells per Cluster",ybins*2,-dy,dy,nceclbins,nceclmin,nceclmax); 
    fhDeltaCellClusterYNCells->SetXTitle("y (cm)");
    fhDeltaCellClusterYNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhDeltaCellClusterYNCells);
    
    fhDeltaCellClusterZNCells  = new TH2F ("hDeltaCellClusterZNCells","Cluster-Cell Z position vs N Cells per Cluster",zbins*2,-dz,dz,nceclbins,nceclmin,nceclmax); 
    fhDeltaCellClusterZNCells->SetXTitle("z (cm)");
    fhDeltaCellClusterZNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhDeltaCellClusterZNCells);
    
    fhDeltaCellClusterRE  = new TH2F ("hDeltaCellClusterRE","Cluster-Cell R position vs cluster energy",rbins*2,-dr,dr,nptbins,ptmin,ptmax); 
    fhDeltaCellClusterRE->SetXTitle("r = #sqrt{x^{2}+y^{2}} (cm)");
    fhDeltaCellClusterRE->SetYTitle("E (GeV)");
    outputContainer->Add(fhDeltaCellClusterRE);		
    
    fhDeltaCellClusterXE  = new TH2F ("hDeltaCellClusterXE","Cluster-Cell X position vs cluster energy",xbins*2,-dx,dx,nptbins,ptmin,ptmax); 
    fhDeltaCellClusterXE->SetXTitle("x (cm)");
    fhDeltaCellClusterXE->SetYTitle("E (GeV)");
    outputContainer->Add(fhDeltaCellClusterXE);
    
    fhDeltaCellClusterYE  = new TH2F ("hDeltaCellClusterYE","Cluster-Cell Y position vs cluster energy",ybins*2,-dy,dy,nptbins,ptmin,ptmax); 
    fhDeltaCellClusterYE->SetXTitle("y (cm)");
    fhDeltaCellClusterYE->SetYTitle("E (GeV)");
    outputContainer->Add(fhDeltaCellClusterYE);
    
    fhDeltaCellClusterZE  = new TH2F ("hDeltaCellClusterZE","Cluster-Cell Z position vs cluster energy",zbins*2,-dz,dz,nptbins,ptmin,ptmax); 
    fhDeltaCellClusterZE->SetXTitle("z (cm)");
    fhDeltaCellClusterZE->SetYTitle("E (GeV)");
    outputContainer->Add(fhDeltaCellClusterZE);
    
    fhEtaPhiAmp  = new TH3F ("hEtaPhiAmp","Cell #eta vs cell #phi vs cell energy",netabins,etamin,etamax,nphibins,phimin,phimax,nptbins,ptmin,ptmax); 
    fhEtaPhiAmp->SetXTitle("#eta ");
    fhEtaPhiAmp->SetYTitle("#phi (rad)");
    fhEtaPhiAmp->SetZTitle("E (GeV) ");
    outputContainer->Add(fhEtaPhiAmp);		
    
  }
  
  //Calo cells
  fhNCells  = new TH1F ("hNCells","# cells", ncebins,ncemin,ncemax); 
  fhNCells->SetXTitle("n cells");
  outputContainer->Add(fhNCells);
  
  fhAmplitude  = new TH1F ("hAmplitude","Cell Energy", nptbins*2,ptmin,ptmax); 
  fhAmplitude->SetXTitle("Cell Energy (GeV)");
  outputContainer->Add(fhAmplitude);
  
  fhAmpId  = new TH2F ("hAmpId","Cell Energy", nfineptbins,ptfinemin,ptfinemax,fNMaxRows*fNMaxCols*fNModules,0,fNMaxRows*fNMaxCols*fNModules); 
  fhAmpId->SetXTitle("Cell Energy (GeV)");
  outputContainer->Add(fhAmpId);
  
  //Cell Time histograms, time only available in ESDs
  if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
    
    fhCellTimeSpreadRespectToCellMax = new TH2F ("hCellTimeSpreadRespectToCellMax","t_{cell max}-t_{cell i} per cluster", nptbins,ptmin,ptmax,tdbins,tdmin,tdmax); 
    fhCellTimeSpreadRespectToCellMax->SetXTitle("E (GeV)");
    fhCellTimeSpreadRespectToCellMax->SetYTitle("#Delta t_{cell max-i} (ns)");
    outputContainer->Add(fhCellTimeSpreadRespectToCellMax);
    
    fhClusterMaxCellDiffAverageTime = new TH2F ("hClusterMaxCellDiffAverageTime","t_{cell max}-t_{average} per cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax); 
    fhClusterMaxCellDiffAverageTime->SetXTitle("E (GeV)");
    fhClusterMaxCellDiffAverageTime->SetYTitle("#Delta t_{cell max - average} (ns)");
    outputContainer->Add(fhClusterMaxCellDiffAverageTime);
    
    fhClusterMaxCellDiffAverageNoMaxTime = new TH2F ("hClusterMaxCellDiffAverageNoMaxTime","t_{cell max}-t_{average} per cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax); 
    fhClusterMaxCellDiffAverageNoMaxTime->SetXTitle("E (GeV)");
    fhClusterMaxCellDiffAverageNoMaxTime->SetYTitle("#Delta t_{cell max - average} (ns)");
    outputContainer->Add(fhClusterMaxCellDiffAverageNoMaxTime);
    
    fhClusterMaxCellDiffWeightedTime = new TH2F ("hClusterMaxCellDiffWeightedTime","t_{cell max}-t_{weighted} per cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax); 
    fhClusterMaxCellDiffWeightedTime->SetXTitle("E (GeV)");
    fhClusterMaxCellDiffWeightedTime->SetYTitle("#Delta t_{cell max - weighted} (ns)");
    outputContainer->Add(fhClusterMaxCellDiffWeightedTime);
    
    fhClusterMaxCellDiffWeightedNoMaxTime = new TH2F ("hClusterMaxCellDiffWeightedNoMaxTime","t_{cell max}-t_{average} per cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax); 
    fhClusterMaxCellDiffWeightedNoMaxTime->SetXTitle("E (GeV)");
    fhClusterMaxCellDiffWeightedNoMaxTime->SetYTitle("#Delta t_{cell max - weighted} (ns)");
    outputContainer->Add(fhClusterMaxCellDiffWeightedNoMaxTime);      
    
    fhCellIdCellLargeTimeSpread= new TH1F ("hCellIdCellLargeTimeSpread","Cells with time 100 ns larger than cell max in cluster ", 
                                           fNMaxCols*fNMaxRows*fNModules,0,fNMaxCols*fNMaxRows*fNModules); 
    fhCellIdCellLargeTimeSpread->SetXTitle("Absolute Cell Id");
    outputContainer->Add(fhCellIdCellLargeTimeSpread);
    
    fhTime  = new TH1F ("hTime","Cell Time",ntimebins,timemin,timemax); 
    fhTime->SetXTitle("Cell Time (ns)");
    outputContainer->Add(fhTime);
    
    fhTimeId  = new TH2F ("hTimeId","Cell Time vs Absolute Id",
                          ntimebins,timemin,timemax,fNMaxRows*fNMaxCols*fNModules,0,fNMaxRows*fNMaxCols*fNModules); 
    fhTimeId->SetXTitle("Cell Time (ns)");
    fhTimeId->SetYTitle("Cell Absolute Id");
    outputContainer->Add(fhTimeId);
    
    fhTimeAmp  = new TH2F ("hTimeAmp","Cell Time vs Cell Energy",nptbins*2,ptmin,ptmax,ntimebins,timemin,timemax); 
    fhTimeAmp->SetYTitle("Cell Time (ns)");
    fhTimeAmp->SetXTitle("Cell Energy (GeV)");
    outputContainer->Add(fhTimeAmp);
    
  }
  
	
  if(fCorrelate){
    //PHOS vs EMCAL
    fhCaloCorrNClusters  = new TH2F ("hCaloCorrNClusters","# clusters in EMCAL vs PHOS", nclbins,nclmin,nclmax,nclbins,nclmin,nclmax); 
    fhCaloCorrNClusters->SetXTitle("number of clusters in EMCAL");
    fhCaloCorrNClusters->SetYTitle("number of clusters in PHOS");
    outputContainer->Add(fhCaloCorrNClusters);
    
    fhCaloCorrEClusters  = new TH2F ("hCaloCorrEClusters","summed energy of clusters in EMCAL vs PHOS", nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
    fhCaloCorrEClusters->SetXTitle("#Sigma E of clusters in EMCAL (GeV)");
    fhCaloCorrEClusters->SetYTitle("#Sigma E of clusters in PHOS (GeV)");
    outputContainer->Add(fhCaloCorrEClusters);
    
    fhCaloCorrNCells  = new TH2F ("hCaloCorrNCells","# Cells in EMCAL vs PHOS", ncebins,ncemin,ncemax, ncebins,ncemin,ncemax); 
    fhCaloCorrNCells->SetXTitle("number of Cells in EMCAL");
    fhCaloCorrNCells->SetYTitle("number of Cells in PHOS");
    outputContainer->Add(fhCaloCorrNCells);
    
    fhCaloCorrECells  = new TH2F ("hCaloCorrECells","summed energy of Cells in EMCAL vs PHOS", nptbins*2,ptmin,ptmax*2,nptbins*2,ptmin,ptmax*2); 
    fhCaloCorrECells->SetXTitle("#Sigma E of Cells in EMCAL (GeV)");
    fhCaloCorrECells->SetYTitle("#Sigma E of Cells in PHOS (GeV)");
    outputContainer->Add(fhCaloCorrECells);
    
    //Calorimeter VS V0 signal
    fhCaloV0SCorrNClusters  = new TH2F ("hCaloV0SNClusters",Form("# clusters in %s vs V0 signal",fCalorimeter.Data()), nv0sbins,nv0smin,nv0smax,nclbins,nclmin,nclmax); 
    fhCaloV0SCorrNClusters->SetXTitle("V0 signal");
    fhCaloV0SCorrNClusters->SetYTitle(Form("number of clusters in %s",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0SCorrNClusters);
    
    fhCaloV0SCorrEClusters  = new TH2F ("hCaloV0SEClusters",Form("summed energy of clusters in %s vs V0 signal",fCalorimeter.Data()), nv0sbins,nv0smin,nv0smax,nptbins,ptmin,ptmax); 
    fhCaloV0SCorrEClusters->SetXTitle("V0 signal");
    fhCaloV0SCorrEClusters->SetYTitle(Form("#Sigma E of clusters in %s (GeV)",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0SCorrEClusters);
    
    fhCaloV0SCorrNCells  = new TH2F ("hCaloV0SNCells",Form("# Cells in %s vs V0 signal",fCalorimeter.Data()), nv0sbins,nv0smin,nv0smax, ncebins,ncemin,ncemax); 
    fhCaloV0SCorrNCells->SetXTitle("V0 signal");
    fhCaloV0SCorrNCells->SetYTitle(Form("number of Cells in %s",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0SCorrNCells);
    
    fhCaloV0SCorrECells  = new TH2F ("hCaloV0SECells",Form("summed energy of Cells in %s vs V0 signal",fCalorimeter.Data()), nv0sbins,nv0smin,nv0smax,nptbins,ptmin,ptmax); 
    fhCaloV0SCorrECells->SetXTitle("V0 signal");
    fhCaloV0SCorrECells->SetYTitle(Form("#Sigma E of Cells in %s (GeV)",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0SCorrECells);    
    
    //Calorimeter VS V0 multiplicity
    fhCaloV0MCorrNClusters  = new TH2F ("hCaloV0MNClusters",Form("# clusters in %s vs V0 signal",fCalorimeter.Data()), nv0mbins,nv0mmin,nv0mmax,nclbins,nclmin,nclmax); 
    fhCaloV0MCorrNClusters->SetXTitle("V0 signal");
    fhCaloV0MCorrNClusters->SetYTitle(Form("number of clusters in %s",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0MCorrNClusters);
    
    fhCaloV0MCorrEClusters  = new TH2F ("hCaloV0MEClusters",Form("summed energy of clusters in %s vs V0 signal",fCalorimeter.Data()), nv0mbins,nv0mmin,nv0mmax,nptbins,ptmin,ptmax); 
    fhCaloV0MCorrEClusters->SetXTitle("V0 signal");
    fhCaloV0MCorrEClusters->SetYTitle(Form("#Sigma E of clusters in %s (GeV)",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0MCorrEClusters);
    
    fhCaloV0MCorrNCells  = new TH2F ("hCaloV0MNCells",Form("# Cells in %s vs V0 signal",fCalorimeter.Data()), nv0mbins,nv0mmin,nv0mmax, ncebins,ncemin,ncemax); 
    fhCaloV0MCorrNCells->SetXTitle("V0 signal");
    fhCaloV0MCorrNCells->SetYTitle(Form("number of Cells in %s",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0MCorrNCells);
    
    fhCaloV0MCorrECells  = new TH2F ("hCaloV0MECells",Form("summed energy of Cells in %s vs V0 signal",fCalorimeter.Data()), nv0mbins,nv0mmin,nv0mmax,nptbins,ptmin,ptmax); 
    fhCaloV0MCorrECells->SetXTitle("V0 signal");
    fhCaloV0MCorrECells->SetYTitle(Form("#Sigma E of Cells in %s (GeV)",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0MCorrECells);    
    
    //Calorimeter VS Track multiplicity
    fhCaloTrackMCorrNClusters  = new TH2F ("hCaloTrackMNClusters",Form("# clusters in %s vs # tracks",fCalorimeter.Data()), ntrmbins,ntrmmin,ntrmmax,nclbins,nclmin,nclmax); 
    fhCaloTrackMCorrNClusters->SetXTitle("# tracks");
    fhCaloTrackMCorrNClusters->SetYTitle(Form("number of clusters in %s",fCalorimeter.Data()));
    outputContainer->Add(fhCaloTrackMCorrNClusters);
    
    fhCaloTrackMCorrEClusters  = new TH2F ("hCaloTrackMEClusters",Form("summed energy of clusters in %s vs # tracks",fCalorimeter.Data()), ntrmbins,ntrmmin,ntrmmax,nptbins,ptmin,ptmax); 
    fhCaloTrackMCorrEClusters->SetXTitle("# tracks");
    fhCaloTrackMCorrEClusters->SetYTitle(Form("#Sigma E of clusters in %s (GeV)",fCalorimeter.Data()));
    outputContainer->Add(fhCaloTrackMCorrEClusters);
    
    fhCaloTrackMCorrNCells  = new TH2F ("hCaloTrackMNCells",Form("# Cells in %s vs # tracks",fCalorimeter.Data()), ntrmbins,ntrmmin,ntrmmax, ncebins,ncemin,ncemax); 
    fhCaloTrackMCorrNCells->SetXTitle("# tracks");
    fhCaloTrackMCorrNCells->SetYTitle(Form("number of Cells in %s",fCalorimeter.Data()));
    outputContainer->Add(fhCaloTrackMCorrNCells);
    
    fhCaloTrackMCorrECells  = new TH2F ("hCaloTrackMECells",Form("summed energy of Cells in %s vs # tracks",fCalorimeter.Data()), ntrmbins,ntrmmin,ntrmmax,nptbins,ptmin,ptmax); 
    fhCaloTrackMCorrECells->SetXTitle("# tracks");
    fhCaloTrackMCorrECells->SetYTitle(Form("#Sigma E of Cells in %s (GeV)",fCalorimeter.Data()));
    outputContainer->Add(fhCaloTrackMCorrECells);    
    
    
  }//correlate calorimeters
  
  //Module histograms
  
  fhEMod  = new TH2F ("hE_Mod","Cluster reconstructed Energy in each present Module",nptbins,ptmin,ptmax,fNModules,0,fNModules); 
  fhEMod->SetXTitle("E (GeV)");
  fhEMod->SetYTitle("Module");
  outputContainer->Add(fhEMod);
  
  fhAmpMod  = new TH2F ("hAmp_Mod","Cell energy in each present Module",nptbins,ptmin,ptmax,fNModules,0,fNModules); 
  fhAmpMod->SetXTitle("E (GeV)");
  fhAmpMod->SetYTitle("Module");
  outputContainer->Add(fhAmpMod);
  
  fhTimeMod  = new TH2F ("hTime_Mod","Cell time in each present Module",ntimebins,timemin,timemax,fNModules,0,fNModules); 
  fhTimeMod->SetXTitle("t (ns)");
  fhTimeMod->SetYTitle("Module");
  outputContainer->Add(fhTimeMod);
  
  fhNClustersMod  = new TH2F ("hNClusters_Mod","# clusters vs Module", nclbins,nclmin,nclmax,fNModules,0,fNModules); 
  fhNClustersMod->SetXTitle("number of clusters");
  fhNClustersMod->SetYTitle("Module");
  outputContainer->Add(fhNClustersMod);
  
  fhNCellsMod  = new TH2F ("hNCells_Mod","# cells vs Module", ncebins,ncemin,ncemax,fNModules,0,fNModules); 
  fhNCellsMod->SetXTitle("n cells");
  fhNCellsMod->SetYTitle("Module");
  outputContainer->Add(fhNCellsMod);
  
  Int_t colmaxs = fNMaxCols;
  Int_t rowmaxs = fNMaxRows;
  if(fCalorimeter=="EMCAL"){
    colmaxs=2*fNMaxCols;
    rowmaxs=Int_t(fNModules/2)*fNMaxRows;
  }
  else{
    rowmaxs=fNModules*fNMaxRows;
  }
  
  fhGridCells  = new TH2F ("hGridCells",Form("Entries in grid of cells"), 
                           colmaxs+2,-1.5,colmaxs+0.5, rowmaxs+2,-1.5,rowmaxs+0.5); 
  fhGridCells->SetYTitle("row (phi direction)");
  fhGridCells->SetXTitle("column (eta direction)");
  outputContainer->Add(fhGridCells);
  
  fhGridCellsE  = new TH2F ("hGridCellsE","Accumulated energy in grid of cells", 
                            colmaxs+2,-1.5,colmaxs+0.5, rowmaxs+2,-1.5,rowmaxs+0.5); 
  fhGridCellsE->SetYTitle("row (phi direction)");
  fhGridCellsE->SetXTitle("column (eta direction)");
  outputContainer->Add(fhGridCellsE);
  
  fhGridCellsTime  = new TH2F ("hGridCellsTime","Accumulated time in grid of cells", 
                               colmaxs+2,-1.5,colmaxs+0.5, rowmaxs+2,-1.5,rowmaxs+0.5); 
  fhGridCellsTime->SetYTitle("row (phi direction)");
  fhGridCellsTime->SetXTitle("column (eta direction)");
  outputContainer->Add(fhGridCellsTime);  
  
  fhNCellsPerClusterMod  = new TH2F*[fNModules];
  fhNCellsPerClusterModNoCut = new TH2F*[fNModules];
  fhTimeAmpPerRCU        = new TH2F*[fNModules*fNRCU];
  fhIMMod                = new TH2F*[fNModules];
  
  for(Int_t imod = 0; imod < fNModules; imod++){
    
    fhNCellsPerClusterMod[imod]  = new TH2F (Form("hNCellsPerCluster_Mod%d",imod),
                                             Form("# cells per cluster vs cluster energy in Module %d",imod), 
                                             nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
    fhNCellsPerClusterMod[imod]->SetXTitle("E (GeV)");
    fhNCellsPerClusterMod[imod]->SetYTitle("n cells");
    outputContainer->Add(fhNCellsPerClusterMod[imod]);
    
    fhNCellsPerClusterModNoCut[imod]  = new TH2F (Form("hNCellsPerClusterNoCut_Mod%d",imod),
                                                  Form("# cells per cluster vs cluster energy in Module %d, no cut",imod), 
                                                  nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
    fhNCellsPerClusterModNoCut[imod]->SetXTitle("E (GeV)");
    fhNCellsPerClusterModNoCut[imod]->SetYTitle("n cells");
    outputContainer->Add(fhNCellsPerClusterModNoCut[imod]);
    
    if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
      
      for(Int_t ircu = 0; ircu < fNRCU; ircu++){
        fhTimeAmpPerRCU[imod*fNRCU+ircu]  = new TH2F (Form("hTimeAmp_Mod%d_RCU%d",imod,ircu),
                                                      Form("Cell Energy vs Cell Time in Module %d, RCU %d ",imod,ircu), 
                                                      nptbins,ptmin,ptmax,ntimebins,timemin,timemax); 
        fhTimeAmpPerRCU[imod*fNRCU+ircu]->SetXTitle("E (GeV)");
        fhTimeAmpPerRCU[imod*fNRCU+ircu]->SetYTitle("time (ns)");
        outputContainer->Add(fhTimeAmpPerRCU[imod*fNRCU+ircu]);
        
      }
    }
    
    if(fFillAllPi0Histo){
      fhIMMod[imod]  = new TH2F (Form("hIM_Mod%d",imod),
                                 Form("Cluster pairs Invariant mass vs reconstructed pair energy in Module %d, n cell > 1",imod),
                                 nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
      fhIMMod[imod]->SetXTitle("p_{T, cluster pairs} (GeV) ");
      fhIMMod[imod]->SetYTitle("M_{cluster pairs} (GeV/c^{2})");
      outputContainer->Add(fhIMMod[imod]);
      
    }
  }
  
  //Monte Carlo Histograms
  
  TString particleName[] = { "Photon", "Pi0", "Eta", "Electron", "NeutralHadron", "ChargedHadron" };
  
  if(IsDataMC()){
    for(Int_t iPart = 0; iPart < 6; iPart++){
      for(Int_t iCh = 0; iCh < 2; iCh++){
        
        fhRecoMCRatioE[iPart][iCh]  = new TH2F (Form("hRecoMCRatioE_%s_Match%d",particleName[iPart].Data(),iCh),
                                                Form("Reco/Gen E, %s, Matched %d",particleName[iPart].Data(),iCh), 
                                                nptbins, ptmin, ptmax, 200,0,2); 
        fhRecoMCRatioE[iPart][iCh]->SetXTitle("E_{reconstructed}/E_{generated}");
        outputContainer->Add(fhRecoMCRatioE[iPart][iCh]);
        
        
        fhRecoMCDeltaE[iPart][iCh]  = new TH2F (Form("hRecoMCDeltaE_%s_Match%d",particleName[iPart].Data(),iCh),
                                                Form("MC - Reco E, %s, Matched %d",particleName[iPart].Data(),iCh), 
                                                nptbins, ptmin, ptmax, nptbins*2,-ptmax,ptmax); 
        fhRecoMCDeltaE[iPart][iCh]->SetXTitle("#Delta E (GeV)");
        outputContainer->Add(fhRecoMCDeltaE[iPart][iCh]);
        
        fhRecoMCDeltaPhi[iPart][iCh]  = new TH2F (Form("hRecoMCDeltaPhi_%s_Match%d",particleName[iPart].Data(),iCh),
                                                  Form("MC - Reco #phi, %s, Matched %d",particleName[iPart].Data(),iCh),
                                                  nptbins, ptmin, ptmax, nphibins*2,-phimax,phimax); 
        fhRecoMCDeltaPhi[iPart][iCh]->SetXTitle("#Delta #phi (rad)");
        outputContainer->Add(fhRecoMCDeltaPhi[iPart][iCh]);
        
        fhRecoMCDeltaEta[iPart][iCh]  = new TH2F (Form("hRecoMCDeltaEta_%s_Match%d",particleName[iPart].Data(),iCh),
                                                  Form("MC- Reco #eta, %s, Matched %d",particleName[iPart].Data(),iCh),
                                                  nptbins, ptmin, ptmax,netabins*2,-etamax,etamax); 
        fhRecoMCDeltaEta[iPart][iCh]->SetXTitle("#Delta #eta ");
        outputContainer->Add(fhRecoMCDeltaEta[iPart][iCh]);
        
        fhRecoMCE[iPart][iCh]  = new TH2F (Form("hRecoMCE_%s_Match%d",particleName[iPart].Data(),iCh),
                                           Form("E distribution, reconstructed vs generated, %s, Matched %d",particleName[iPart].Data(),iCh),
                                           nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
        fhRecoMCE[iPart][iCh]->SetXTitle("E_{rec} (GeV)");
        fhRecoMCE[iPart][iCh]->SetYTitle("E_{gen} (GeV)");
        outputContainer->Add(fhRecoMCE[iPart][iCh]);	  
        
        fhRecoMCPhi[iPart][iCh]  = new TH2F (Form("hRecoMCPhi_%s_Match%d",particleName[iPart].Data(),iCh),
                                             Form("#phi distribution, reconstructed vs generated, %s, Matched %d",particleName[iPart].Data(),iCh),
                                             nphibins,phimin,phimax, nphibins,phimin,phimax); 
        fhRecoMCPhi[iPart][iCh]->SetXTitle("#phi_{rec} (rad)");
        fhRecoMCPhi[iPart][iCh]->SetYTitle("#phi_{gen} (rad)");
        outputContainer->Add(fhRecoMCPhi[iPart][iCh]);
        
        fhRecoMCEta[iPart][iCh]  = new TH2F (Form("hRecoMCEta_%s_Match%d",particleName[iPart].Data(),iCh),
                                             Form("#eta distribution, reconstructed vs generated, %s, Matched %d",particleName[iPart].Data(),iCh), 
                                             netabins,etamin,etamax,netabins,etamin,etamax); 
        fhRecoMCEta[iPart][iCh]->SetXTitle("#eta_{rec} ");
        fhRecoMCEta[iPart][iCh]->SetYTitle("#eta_{gen} ");
        outputContainer->Add(fhRecoMCEta[iPart][iCh]);
      }
    }  
    
    //Pure MC
    for(Int_t iPart = 0; iPart < 4; iPart++){
      fhGenMCE[iPart]     = new TH1F(Form("hGenMCE_%s",particleName[iPart].Data()) ,
                                     Form("p_{T} of generated %s",particleName[iPart].Data()),
                                     nptbins,ptmin,ptmax);
      fhGenMCEtaPhi[iPart] = new TH2F(Form("hGenMCEtaPhi_%s",particleName[iPart].Data()),
                                      Form("Y vs #phi of generated %s",particleName[iPart].Data()),
                                      netabins,etamin,etamax,nphibins,phimin,phimax);
    	
      fhGenMCE[iPart]    ->SetXTitle("p_{T} (GeV/c)");
      fhGenMCEtaPhi[iPart]->SetXTitle("#eta");
      fhGenMCEtaPhi[iPart]->SetYTitle("#phi (rad)");
      
      outputContainer->Add(fhGenMCE[iPart]);
      outputContainer->Add(fhGenMCEtaPhi[iPart]);
      
      
      fhGenMCAccE[iPart]     = new TH1F(Form("hGenMCAccE_%s",particleName[iPart].Data()) ,
                                        Form("p_{T} of generated %s",particleName[iPart].Data()),
                                        nptbins,ptmin,ptmax);
      fhGenMCAccEtaPhi[iPart] = new TH2F(Form("hGenMCAccEtaPhi_%s",particleName[iPart].Data()),
                                         Form("Y vs #phi of generated %s",particleName[iPart].Data()),
                                         netabins,etamin,etamax,nphibins,phimin,phimax);
    	
      fhGenMCAccE[iPart]    ->SetXTitle("p_{T} (GeV/c)");
      fhGenMCAccEtaPhi[iPart]->SetXTitle("#eta");
      fhGenMCAccEtaPhi[iPart]->SetYTitle("#phi (rad)");
      
      outputContainer->Add(fhGenMCAccE[iPart]);
      outputContainer->Add(fhGenMCAccEtaPhi[iPart]);
      
    }    
    
    //Vertex of generated particles 
    
    fhEMVxyz  = new TH2F ("hEMVxyz","Production vertex of reconstructed ElectroMagnetic particles",nvdistbins,vdistmin,vdistmax,nvdistbins,vdistmin,vdistmax);//,100,0,500); 
    fhEMVxyz->SetXTitle("v_{x}");
    fhEMVxyz->SetYTitle("v_{y}");
    //fhEMVxyz->SetZTitle("v_{z}");
    outputContainer->Add(fhEMVxyz);
    
    fhHaVxyz  = new TH2F ("hHaVxyz","Production vertex of reconstructed hadrons",nvdistbins,vdistmin,vdistmax,nvdistbins,vdistmin,vdistmax);//,100,0,500); 
    fhHaVxyz->SetXTitle("v_{x}");
    fhHaVxyz->SetYTitle("v_{y}");
    //fhHaVxyz->SetZTitle("v_{z}");
    outputContainer->Add(fhHaVxyz);
    
    fhEMR  = new TH2F ("hEMR","Distance to production vertex of reconstructed ElectroMagnetic particles vs E rec",nptbins,ptmin,ptmax,nvdistbins,vdistmin,vdistmax); 
    fhEMR->SetXTitle("E (GeV)");
    fhEMR->SetYTitle("TMath::Sqrt(v_{x}^{2}+v_{y}^{2})");
    outputContainer->Add(fhEMR);
    
    fhHaR  = new TH2F ("hHaR","Distance to production vertex of reconstructed Hadrons vs E rec",nptbins,ptmin,ptmax,nvdistbins,vdistmin,vdistmax); 
    fhHaR->SetXTitle("E (GeV)");
    fhHaR->SetYTitle("TMath::Sqrt(v_{x}^{2}+v_{y}^{2})");
    outputContainer->Add(fhHaR);
    
    
    //Track Matching 
    
    fhMCEle1pOverE = new TH2F("hMCEle1pOverE","TRACK matches p/E, MC electrons",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fhMCEle1pOverE->SetYTitle("p/E");
    fhMCEle1pOverE->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMCEle1pOverE);
    
    fhMCEle1dR = new TH1F("hMCEle1dR","TRACK matches dR, MC electrons",ndRbins,dRmin,dRmax);
    fhMCEle1dR->SetXTitle("#Delta R (rad)");
    outputContainer->Add(fhMCEle1dR) ;
    
    fhMCEle2MatchdEdx = new TH2F("hMCEle2MatchdEdx","dE/dx vs. p for all matches, MC electrons",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fhMCEle2MatchdEdx->SetXTitle("p (GeV/c)");
    fhMCEle2MatchdEdx->SetYTitle("<dE/dx>");
    outputContainer->Add(fhMCEle2MatchdEdx);
    
    fhMCChHad1pOverE = new TH2F("hMCChHad1pOverE","TRACK matches p/E, MC charged hadrons",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fhMCChHad1pOverE->SetYTitle("p/E");
    fhMCChHad1pOverE->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMCChHad1pOverE);
    
    fhMCChHad1dR = new TH1F("hMCChHad1dR","TRACK matches dR, MC charged hadrons",ndRbins,dRmin,dRmax);
    fhMCChHad1dR->SetXTitle("#Delta R (rad)");
    outputContainer->Add(fhMCChHad1dR) ;
    
    fhMCChHad2MatchdEdx = new TH2F("hMCChHad2MatchdEdx","dE/dx vs. p for all matches, MC charged hadrons",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fhMCChHad2MatchdEdx->SetXTitle("p (GeV/c)");
    fhMCChHad2MatchdEdx->SetYTitle("<dE/dx>");
    outputContainer->Add(fhMCChHad2MatchdEdx);
    
    fhMCNeutral1pOverE = new TH2F("hMCNeutral1pOverE","TRACK matches p/E, MC neutrals",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fhMCNeutral1pOverE->SetYTitle("p/E");
    fhMCNeutral1pOverE->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMCNeutral1pOverE);
    
    fhMCNeutral1dR = new TH1F("hMCNeutral1dR","TRACK matches dR, MC neutrals",ndRbins,dRmin,dRmax);
    fhMCNeutral1dR->SetXTitle("#Delta R (rad)");
    outputContainer->Add(fhMCNeutral1dR) ;
    
    fhMCNeutral2MatchdEdx = new TH2F("hMCNeutral2MatchdEdx","dE/dx vs. p for all matches, MC neutrals",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fhMCNeutral2MatchdEdx->SetXTitle("p (GeV/c)");
    fhMCNeutral2MatchdEdx->SetYTitle("<dE/dx>");
    outputContainer->Add(fhMCNeutral2MatchdEdx);
    
    fhMCEle1pOverER02 = new TH2F("hMCEle1pOverER02","TRACK matches p/E, MC electrons",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fhMCEle1pOverER02->SetYTitle("p/E");
    fhMCEle1pOverER02->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMCEle1pOverER02);
    
    fhMCChHad1pOverER02 = new TH2F("hMCChHad1pOverER02","TRACK matches p/E, MC charged hadrons",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fhMCChHad1pOverER02->SetYTitle("p/E");
    fhMCChHad1pOverER02->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMCChHad1pOverER02);
    
    fhMCNeutral1pOverER02 = new TH2F("hMCNeutral1pOverER02","TRACK matches p/E, MC neutrals",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fhMCNeutral1pOverER02->SetYTitle("p/E");
    fhMCNeutral1pOverER02->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMCNeutral1pOverER02);
  }
  
  //  for(Int_t i = 0; i < outputContainer->GetEntries() ; i++)
  //    printf("i=%d, name= %s\n",i,outputContainer->At(i)->GetName());
  
  return outputContainer;
}

//_____________________________________________________________________________________________
void AliAnaCalorimeterQA::InvariantMassHistograms(const Int_t iclus,   const TLorentzVector mom, 
                                                  const Int_t nModule, TObjArray* caloClusters) 
{
  // Fill Invariant mass histograms
  
  if(GetDebug()>1) printf("AliAnaCalorimeterQA::InvariantMassHistograms() - Start \n");
  
  //Get vertex for photon momentum calculation and event selection
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  
  Int_t nModule2      = -1;
  TLorentzVector mom2 ;
  Int_t nCaloClusters = caloClusters->GetEntriesFast();
  
  for(Int_t jclus = iclus + 1 ; jclus < nCaloClusters ; jclus++) {
    AliVCluster* clus2 =  (AliVCluster*)caloClusters->At(jclus);
    
    if( clus2->GetNCells() <= 1 || !IsGoodCluster(clus2->E(),clus2->GetNCells())) continue;
    
    //Get cluster kinematics
    clus2->GetMomentum(mom2,v);
    
    //Check only certain regions
    Bool_t in2 = kTRUE;
    if(IsFiducialCutOn()) in2 =  GetFiducialCut()->IsInFiducialCut(mom2,fCalorimeter) ;
    if(!in2) continue;	
    
    //Get module of cluster
    nModule2 = GetModuleNumber(clus2);
    
    //Fill histograms
    
    //All modules
    fhIM  ->Fill((mom+mom2).Pt(),(mom+mom2).M());

    //Single module
    if(nModule == nModule2 && nModule >=0 && nModule < fNModules){
      fhIMMod[nModule]->Fill((mom+mom2).Pt(),(mom+mom2).M());
    }
    
    //Asymetry histograms
    fhAsym->Fill((mom+mom2).Pt(),TMath::Abs((mom.E()-mom2.E())/(mom.E()+mom2.E())));
    
  }// 2nd cluster loop
  
}

//______________________________
void AliAnaCalorimeterQA::Init()
{ 
  //Check if the data or settings are ok
  
  if(fCalorimeter != "PHOS" && fCalorimeter !="EMCAL")
    AliFatal(Form("Wrong calorimeter name <%s>", fCalorimeter.Data()));
  
  if(GetReader()->GetDataType()== AliCaloTrackReader::kMC)
    AliFatal("Analysis of reconstructed data, MC reader not aplicable");
  
}

//________________________________________
void AliAnaCalorimeterQA::InitParameters()
{ 
  //Initialize the parameters of the analysis.
  AddToHistogramsName("AnaCaloQA_");
  
  fCalorimeter     = "EMCAL"; //or PHOS
  fNModules        = 12; // set maximum to maximum number of EMCAL modules
  fNRCU            = 2;  // set maximum number of RCU in EMCAL per SM
  fTimeCutMin      = -1;
  fTimeCutMax      = 9999999;
  fEMCALCellAmpMin = 0.2;
  fPHOSCellAmpMin  = 0.2;
  
}

//________________________________________________________________________________________________
Bool_t AliAnaCalorimeterQA::IsGoodCluster(const Float_t energy, const Int_t nCaloCellsPerCluster)
{
  //Identify cluster as exotic or not
  
  if(fCalorimeter=="EMCAL" && !GetCaloUtils()->GetEMCALRecoUtils()->IsRejectExoticCluster() && fStudyBadClusters){
    Float_t minNCells = 1;
    if(energy > 7) minNCells = TMath::Max(1,TMath::FloorNint(1 + TMath::Log(energy - 7 )*1.7 ));
    if(nCaloCellsPerCluster <= minNCells) {        
      return kFALSE;
    }
  }
  
  return kTRUE;
}

//___________________________________________________
void AliAnaCalorimeterQA::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");
  
  printf("Select Calorimeter %s \n",fCalorimeter.Data());
  printf("Time Cut: %3.1f < TOF  < %3.1f\n", fTimeCutMin, fTimeCutMax);
  printf("EMCAL Min Amplitude   : %2.1f GeV/c\n", fEMCALCellAmpMin) ;
  printf("PHOS Min Amplitude    : %2.1f GeV/c\n", fPHOSCellAmpMin) ;
  
} 

//____________________________________________________________________________
void AliAnaCalorimeterQA::RecalibrateCellTime(Double_t & time, const Int_t id)
{
  // Recalculate time if time recalibration available 
  
  if(fCalorimeter == "EMCAL" && GetCaloUtils()->GetEMCALRecoUtils()->IsTimeRecalibrationOn()) {
    GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTime(id,GetReader()->GetInputEvent()->GetBunchCrossNumber(),time);
  }
  
}

//___________________________________________________________________________________
void AliAnaCalorimeterQA::RecalibrateCellAmplitude(Float_t & amp, const Int_t id)
{
  //Recaculate cell energy if recalibration factor
  
  Int_t icol     = -1; Int_t irow     = -1; Int_t iRCU     = -1;
  Int_t nModule  = GetModuleNumberCellIndexes(id,fCalorimeter, icol, irow, iRCU);
  
  if (GetCaloUtils()->IsRecalibrationOn()) {
    if(fCalorimeter == "PHOS") {
      amp *= GetCaloUtils()->GetPHOSChannelRecalibrationFactor(nModule,icol,irow);
    }
    else		                   {
      amp *= GetCaloUtils()->GetEMCALChannelRecalibrationFactor(nModule,icol,irow);
    }
  }
}

//_____________________________________________________
void  AliAnaCalorimeterQA::MakeAnalysisFillHistograms() 
{
  //Fill Calorimeter QA histograms
  
  //Play with the MC stack if available	
  if(IsDataMC()) MCHistograms();	
  
  //Get List with CaloClusters  
  TObjArray * caloClusters = NULL;
  if      (fCalorimeter == "PHOS")  caloClusters = GetPHOSClusters();
  else if (fCalorimeter == "EMCAL") caloClusters = GetEMCALClusters();
  else 
    AliFatal(Form("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - Wrong calorimeter name <%s>, END\n", fCalorimeter.Data()));
  
  // Do not do anything if there are no clusters
  if(caloClusters->GetEntriesFast() == 0) return;
  
  //Get List with CaloCells
  AliVCaloCells * cells = 0x0; 
  if(fCalorimeter == "PHOS") cells =  GetPHOSCells();
  else			                 cells =  GetEMCALCells();
  
  if(!caloClusters || !cells) {
    AliFatal(Form("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - No CaloClusters or CaloCells available\n"));
    return; // trick coverity
  }
  
  // Correlate Calorimeters and V0 and track Multiplicity  
  if(fCorrelate)	Correlate();
  
  // Clusters 
  ClusterLoopHistograms(caloClusters,cells);
  
  // Cells  
  CellHistograms(cells);
  
  if(GetDebug() > 0)
    printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - End \n");
  
}

//______________________________________
void AliAnaCalorimeterQA::MCHistograms()
{
  //Get the MC arrays and do some checks before filling MC histograms
  
  TLorentzVector mom  ;
  
  if(GetReader()->ReadStack()){
    
    if(!GetMCStack()) 
      AliFatal("Stack not available, is the MC handler called?\n");
    
    //Fill some pure MC histograms, only primaries.
    for(Int_t i=0 ; i<GetMCStack()->GetNprimary(); i++){//Only primary particles, for all MC transport put GetNtrack()
      TParticle *primary = GetMCStack()->Particle(i) ;
      
      if (primary->GetStatusCode() > 11) continue; //Working for PYTHIA and simple generators, check for HERWIG 
      primary->Momentum(mom);
      MCHistograms(mom,TMath::Abs(primary->GetPdgCode()));
    } //primary loop
  }
  else if(GetReader()->ReadAODMCParticles()){
    
    if(!GetReader()->GetAODMCParticles(0)) 	
      AliFatal("AODMCParticles not available!");
    
    //Fill some pure MC histograms, only primaries.
    for(Int_t i=0 ; i < (GetReader()->GetAODMCParticles(0))->GetEntriesFast(); i++){
      AliAODMCParticle *aodprimary = (AliAODMCParticle*) (GetReader()->GetAODMCParticles(0))->At(i) ;
      
      if (!aodprimary->IsPrimary()) continue; //accept all which is not MC transport generated. Don't know how to avoid partons
      
      mom.SetPxPyPzE(aodprimary->Px(),aodprimary->Py(),aodprimary->Pz(),aodprimary->E());
      MCHistograms(mom,TMath::Abs(aodprimary->GetPdgCode()));
    } //primary loop
    
  }
}

//_______________________________________________________________________________
void AliAnaCalorimeterQA::MCHistograms(const TLorentzVector mom, const Int_t pdg)
{
  //Fill pure monte carlo related histograms
	
  Float_t eMC    = mom.E();
  Float_t phiMC  = mom.Phi();
  if(phiMC < 0) 
    phiMC  += TMath::TwoPi();
  Float_t etaMC  = mom.Eta();
  
  if (TMath::Abs(etaMC) > 1) return;
  
  Bool_t in = kFALSE;
  
  //Rough stimate of acceptance for pi0, Eta and electrons
  if(fCalorimeter == "PHOS"){
    
    if(GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter)) 
      in = kTRUE ;
    if(GetDebug() > 2) printf("AliAnaCalorimeterQA::MCHistograms() - In %s fiducial cut acceptance? %d\n",fCalorimeter.Data(),in);
    
  }	   
  else if(fCalorimeter == "EMCAL" && GetCaloUtils()->IsEMCALGeoMatrixSet()){
    if(GetEMCALGeometry()){
      
      Int_t absID=0;
      GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(mom.Eta(),mom.Phi(),absID);
      
      if( absID >= 0) 
        in = kTRUE;
      
      if(GetDebug() > 2) printf("AliAnaCalorimeterQA::MCHistograms() - In %s Real acceptance? %d\n",fCalorimeter.Data(),in);
    }
    else{
      if(GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter)) 
        in = kTRUE ;
      if(GetDebug() > 2) printf("AliAnaCalorimeterQA::MCHistograms() - In %s fiducial cut acceptance? %d\n",fCalorimeter.Data(),in);
    }
  }	  
  
  if (pdg==22) {
    fhGenMCE[mcPhoton]    ->Fill(eMC);
    if(eMC > 0.5) fhGenMCEtaPhi[mcPhoton]->Fill(etaMC,phiMC);
    if(in){
      fhGenMCAccE[mcPhoton]    ->Fill(eMC);
      if(eMC > 0.5) fhGenMCAccEtaPhi[mcPhoton]->Fill(etaMC,phiMC);	
    }
  }
  else if (pdg==111) {
    fhGenMCE[mcPi0]    ->Fill(eMC);
    if(eMC > 0.5) fhGenMCEtaPhi[mcPi0]->Fill(etaMC,phiMC);
    if(in){
      fhGenMCAccE[mcPi0]    ->Fill(eMC);
      if(eMC > 0.5) fhGenMCAccEtaPhi[mcPi0]->Fill(etaMC,phiMC);	
    }
  }
  else if (pdg==221) {
    fhGenMCE[mcEta]    ->Fill(eMC);
    if(eMC > 0.5) fhGenMCEtaPhi[mcEta]->Fill(etaMC,phiMC);
    if(in){
      fhGenMCAccE[mcEta]    ->Fill(eMC);
      if(eMC > 0.5) fhGenMCAccEtaPhi[mcEta]->Fill(etaMC,phiMC);	
    }    
  }
  else if (TMath::Abs(pdg)==11) {
    fhGenMCE[mcElectron]    ->Fill(eMC);
    if(eMC > 0.5) fhGenMCEtaPhi[mcElectron]->Fill(etaMC,phiMC);
    if(in){
      fhGenMCAccE[mcElectron]    ->Fill(eMC);
      if(eMC > 0.5) fhGenMCAccEtaPhi[mcElectron]->Fill(etaMC,phiMC);	
    }
  }	
}

//_________________________________________________________________________________
void AliAnaCalorimeterQA::WeightHistograms(AliVCluster *clus, AliVCaloCells* cells)
{
  // Calculate weights
  
  // First recalculate energy in case non linearity was applied
  Float_t  energy = 0;
  Float_t  ampMax = 0;  
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) {
    
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    RecalibrateCellAmplitude(amp,id);
    
    energy    += amp;
    
    if(amp> ampMax) 
      ampMax = amp;
    
  } // energy loop       
  
  if(energy <=0 ) {
    printf("AliAnaCalorimeterQA::WeightHistograms()- Wrong calculated energy %f\n",energy);
    return;
  }
  
  fhEMaxCellClusterRatio   ->Fill(energy,ampMax/energy);
  fhEMaxCellClusterLogRatio->Fill(energy,TMath::Log(ampMax/energy));
  
  //Get the ratio and log ratio to all cells in cluster
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) {
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    RecalibrateCellAmplitude(amp,id);
    
    fhECellClusterRatio   ->Fill(energy,amp/energy);
    fhECellClusterLogRatio->Fill(energy,TMath::Log(amp/energy));
  }        
  
  //Recalculate shower shape for different W0
  if(fCalorimeter=="EMCAL"){
    
    Float_t l0org = clus->GetM02();
    Float_t l1org = clus->GetM20();
    Float_t dorg  = clus->GetDispersion();
    
    for(Int_t iw = 0; iw < 7; iw++){
      GetCaloUtils()->GetEMCALRecoUtils()->SetW0(3+iw*0.5); 
      GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), cells, clus);
      
      fhLambda0ForW0[iw]->Fill(energy,clus->GetM02());
      fhLambda1ForW0[iw]->Fill(energy,clus->GetM20());
      
      if(IsDataMC()){  
        
        Int_t tag = GetMCAnalysisUtils()->CheckOrigin(clus->GetLabels(),clus->GetNLabels(), GetReader(),0);
        
        if(   GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton) && 
           !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0)      && 
           !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta)      &&
           !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)        ){
          fhLambda0ForW0MC[iw][0]->Fill(energy,clus->GetM02());
          fhLambda1ForW0MC[iw][0]->Fill(energy,clus->GetM20());
        } // Pure Photon
        else if( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCElectron) && 
                !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)        ){
          fhLambda0ForW0MC[iw][1]->Fill(energy,clus->GetM02());
          fhLambda1ForW0MC[iw][1]->Fill(energy,clus->GetM20());
        } // Electron
        else if( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)        ){
          fhLambda0ForW0MC[iw][2]->Fill(energy,clus->GetM02());
          fhLambda1ForW0MC[iw][2]->Fill(energy,clus->GetM20());
        } // Conversion
        else if( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0) ){
          fhLambda0ForW0MC[iw][3]->Fill(energy,clus->GetM02());
          fhLambda1ForW0MC[iw][3]->Fill(energy,clus->GetM20());
        }// Pi0
        else if(!GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta) && 
                !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton) ){
          fhLambda0ForW0MC[iw][4]->Fill(energy,clus->GetM02());
          fhLambda1ForW0MC[iw][4]->Fill(energy,clus->GetM20());
        }// Hadron
        
      }// Is MC
    } // w0 loop
    
    // Set the original values back
    clus->SetM02(l0org);
    clus->SetM20(l1org);
    clus->SetDispersion(dorg);
    
  }// EMCAL
  
}



