/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *

 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <iostream>

// Root 
#include <TRefArray.h>
#include <TList.h>
#include <TH1F.h>
#include <TGeoManager.h>
#include <TFile.h>

// AliRoot
#include "AliAnalysisTaskEMCALPi0CalibSelectionV2.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliMultSelection.h"
#include "AliEMCALGeometry.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliEMCALRecoUtils.h"
#include "AliOADBContainer.h"
#include "AliDataFile.h"
#include "AliEmcalCorrectionTask.h"
#include "AliEmcalCorrectionComponent.h"
#include "AliEmcalCorrectionClusterNonLinearity.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEMCALPi0CalibSelectionV2) ;
/// \endcond

///
/// Default constructor. Arrays initialization is done here.
//______________________________________________________________________________________________
AliAnalysisTaskEMCALPi0CalibSelectionV2::AliAnalysisTaskEMCALPi0CalibSelectionV2() :
AliAnalysisTaskSE(),  
fEMCALGeo(0x0),           fLoadMatrices(0),
fEMCALGeoName("EMCAL_COMPLETE12SMV1_DCAL_8SM"),
fTriggerName("EMC"),      
fRecoUtils(NULL),
fEMCALInitialized(kFALSE),
isMC(kFALSE),
fOADBFilePath(""),        
fRecalPosition(kTRUE),
fCaloClustersArr(0x0),    fEMCALCells(0x0),
//fCuts(0x0),               
fOutputContainer(0x0),
fCheckCentrality(kFALSE), fCentralityClass("V0M"),  fCentWithEventSel(kFALSE),
fCentMin(-1),             fCentMax(10000000), 
fVertex(),                
fImportGeometryFromFile(1), fImportGeometryFilePath(""),
fEmin(0.5),               fEmax(15.),        
fL0min(0.01),             fL0max(0.5),
fOpAnglemin(0.),          fOpAnglemax(3.0),
fDTimeCut(100.),          fTimeMax(1000000),        fTimeMin(-1000000),
fAsyCut(1.),              fMinNCells(2),            fGroupNCells(0),
fSameSM(kFALSE),          
fNMaskCellColumns(11),    fMaskCellColumns(0x0),
fInvMassCutMin(110.),     fInvMassCutMax(160.),
// Histograms binning
fNbins(300),              fMinBin(0.),              fMaxBin(300.),
fNTimeBins(1000),         fMinTimeBin(0.),          fMaxTimeBin(1000.),
fNEnergybins(1000),       fMinEnergyBin(0.),        fMaxEnergyBin(100.),
// Temporal
fMomentum1(),             fMomentum2(),             fMomentum12(),
// Histograms
fHmgg(0x0),               fHmggDifferentSM(0x0), 
fHmggMaskFrame(0x0),      fHmggDifferentSMMaskFrame(0x0), 
fHOpeningAngle(0x0),      fHOpeningAngleDifferentSM(0x0),  
fHAsymmetry(0x0),         fHAsymmetryDifferentSM(0x0),  
fhNEvents(0x0),
fhCentrality(0x0),        fhCentralitySelected(0x0),
fhClusterTime(0x0) {
  
  for(Int_t iMod=0; iMod < AliEMCALGeoParams::fgkEMCALModules; iMod++) {
    for(Int_t iX=0; iX<24; iX++) {
      for(Int_t iZ=0; iZ<48; iZ++) {
        fHmpi0[iMod][iZ][iX]   = 0 ;
      }
    } 
  }
  
  fVertex[0]=fVertex[1]=fVertex[2]=-1000;
  
  fHTpi0[0]= 0 ;
  fHTpi0[1]= 0 ;
  fHTpi0[2]= 0 ;
  fHTpi0[3]= 0 ;
  
  fMaskCellColumns    = new Int_t[fNMaskCellColumns];
  fMaskCellColumns[0] = 6 ;  fMaskCellColumns[1] = 7 ;  fMaskCellColumns[2] = 8 ; 
  fMaskCellColumns[3] = 35;  fMaskCellColumns[4] = 36;  fMaskCellColumns[5] = 37; 
  fMaskCellColumns[6] = 12+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[7] = 13+AliEMCALGeoParams::fgkEMCALCols;
  fMaskCellColumns[8] = 40+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[9] = 41+AliEMCALGeoParams::fgkEMCALCols; 
  fMaskCellColumns[10]= 42+AliEMCALGeoParams::fgkEMCALCols; 

  
  for(Int_t iSM = 0; iSM < AliEMCALGeoParams::fgkEMCALModules; iSM++) {
    fHmggSM[iSM]                     = 0;
    fHmggSMMaskFrame[iSM]            = 0;
    fHOpeningAngleSM[iSM]            = 0;
    fHOpeningAnglePairSM[iSM]        = 0;
    fHAsymmetrySM[iSM]               = 0;
    fHAsymmetryPairSM[iSM]           = 0;
    fhTowerDecayPhotonHit[iSM]       = 0;
    fhTowerDecayPhotonEnergy[iSM]    = 0;
    fhTowerDecayPhotonAsymmetry[iSM] = 0;
    fhTowerDecayPhotonHitMaskFrame[iSM]= 0;
    fMatrix[iSM]                     = 0x0;
    fhClusterTimeSM[iSM]             = 0;
  }
}


/// Constructor with name as option. Arrays initialization is done here.
///
/// \param name: Name of task.
///
//______________________________________________________________________________________________
AliAnalysisTaskEMCALPi0CalibSelectionV2::AliAnalysisTaskEMCALPi0CalibSelectionV2(const char* name) :
AliAnalysisTaskSE(name),  
fEMCALGeo(0x0),           fLoadMatrices(0),
fEMCALGeoName("EMCAL_COMPLETE12SMV1_DCAL_8SM"),
fTriggerName("EMC"),      
fRecoUtils(NULL),
fEMCALInitialized(kFALSE),
isMC(kFALSE),
fOADBFilePath(""),        
fRecalPosition(kTRUE),
fCaloClustersArr(0x0),    fEMCALCells(0x0),
//fCuts(0x0),               
fOutputContainer(0x0),
fCheckCentrality(kFALSE), fCentralityClass("V0M"),  fCentWithEventSel(kFALSE),
fCentMin(-1),             fCentMax(10000000), 
fVertex(),
fImportGeometryFromFile(1), fImportGeometryFilePath(""),
fEmin(0.5),               fEmax(15.),     
fL0min(0.01),             fL0max(0.5),
fOpAnglemin(0.),          fOpAnglemax(3.0),
fDTimeCut(100.),          fTimeMax(1000000),        fTimeMin(-1000000),
fAsyCut(1.),              fMinNCells(2),            fGroupNCells(0),
fSameSM(kFALSE),          
fNMaskCellColumns(11),    fMaskCellColumns(0x0),
fInvMassCutMin(110.),     fInvMassCutMax(160.),
// Histograms binning
fNbins(300),              fMinBin(0.),              fMaxBin(300.),
fNTimeBins(1000),         fMinTimeBin(0.),          fMaxTimeBin(1000.),
fNEnergybins(1000),       fMinEnergyBin(0.),        fMaxEnergyBin(100.),
// Temporal
fMomentum1(),             fMomentum2(),             fMomentum12(),
// Histograms
fHmgg(0x0),               fHmggDifferentSM(0x0), 
fHmggMaskFrame(0x0),      fHmggDifferentSMMaskFrame(0x0), 
fHOpeningAngle(0x0),      fHOpeningAngleDifferentSM(0x0),  
fHAsymmetry(0x0),         fHAsymmetryDifferentSM(0x0),  
fhNEvents(0x0),
fhCentrality(0x0),        fhCentralitySelected(0x0),
fhClusterTime(0x0) {
  
  for(Int_t iMod=0; iMod < AliEMCALGeoParams::fgkEMCALModules; iMod++) {
    for(Int_t iX=0; iX<24; iX++) {
      for(Int_t iZ=0; iZ<48; iZ++) {
        fHmpi0[iMod][iZ][iX]   = 0 ;
      }
    } 
  }
  
  fVertex[0]=fVertex[1]=fVertex[2]=-1000;
  
  fHTpi0[0]= 0 ;
  fHTpi0[1]= 0 ;
  fHTpi0[2]= 0 ;
  fHTpi0[3]= 0 ;
  
  fMaskCellColumns    = new Int_t[fNMaskCellColumns];
  fMaskCellColumns[0] = 6 ;  fMaskCellColumns[1] = 7 ;  fMaskCellColumns[2] = 8 ; 
  fMaskCellColumns[3] = 35;  fMaskCellColumns[4] = 36;  fMaskCellColumns[5] = 37; 
  fMaskCellColumns[6] = 12+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[7] = 13+AliEMCALGeoParams::fgkEMCALCols;
  fMaskCellColumns[8] = 40+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[9] = 41+AliEMCALGeoParams::fgkEMCALCols; 
  fMaskCellColumns[10]= 42+AliEMCALGeoParams::fgkEMCALCols; 

  
  for(Int_t iSM = 0; iSM < AliEMCALGeoParams::fgkEMCALModules; iSM++) {
    fHmggSM[iSM]                     = 0;
    fHmggSMMaskFrame[iSM]            = 0;
    fHOpeningAngleSM[iSM]            = 0;
    fHOpeningAnglePairSM[iSM]        = 0;
    fHAsymmetrySM[iSM]               = 0;
    fHAsymmetryPairSM[iSM]           = 0;
    fhTowerDecayPhotonHit[iSM]       = 0;
    fhTowerDecayPhotonEnergy[iSM]    = 0;
    fhTowerDecayPhotonAsymmetry[iSM] = 0;
    fhTowerDecayPhotonHitMaskFrame[iSM]= 0;
    fMatrix[iSM]                     = 0x0;
    fhClusterTimeSM[iSM]             = 0;
  }
    
  DefineOutput(1, TList::Class());
//DefineOutput(2, TList::Class());  // will contain cuts or local params
}

///
/// Destructor.
//_____________________________________________________________________________
AliAnalysisTaskEMCALPi0CalibSelectionV2::~AliAnalysisTaskEMCALPi0CalibSelectionV2() {
  
  if(fOutputContainer) {
    fOutputContainer->Delete() ; 
    delete fOutputContainer ;
  }
  
  if(fEMCALGeo)         delete fEMCALGeo  ;
  if(fNMaskCellColumns) delete [] fMaskCellColumns;
}

///
/// Accept cluster routine
//__________________________________________________________
Bool_t AliAnalysisTaskEMCALPi0CalibSelectionV2::AcceptCluster( AliVCluster* c1){
  // Exclude bad channels
  
  Float_t e1i = c1->E();   // cluster energy before correction   
  
  if (!c1->IsEMCAL()) 
    return kFALSE;
  if (c1->GetNCells() < fMinNCells) 
    return kFALSE;
  if (e1i > fEmax) 
    return kFALSE;
  if (e1i < fEmin) 
    return kFALSE;
  if (c1->GetIsExotic())
    return kFALSE;
  if ( !isMC && ( c1->GetTOF()*1.e9 > fTimeMax || c1->GetTOF()*1.e9 < fTimeMin))
    return kFALSE;
  if (c1->GetM02() < fL0min || c1->GetM02() > fL0max)
    return kFALSE;
  
  return kTRUE;

}

///
/// Initialize EMCAL
//__________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::InitializeEMCAL( AliVEvent *event ){
    AliEmcalCorrectionTask* emcalCorrTask=0x0;

    emcalCorrTask  = (AliEmcalCorrectionTask*) AliAnalysisManager::GetAnalysisManager()->GetTopTasks()->FindObject("AliEmcalCorrectionTask");

    if( emcalCorrTask ){
      AliEmcalCorrectionComponent * emcalCorrComponent = 0x0;
      emcalCorrComponent = emcalCorrTask->GetCorrectionComponent("AliEmcalCorrectionClusterNonLinearity");
      if( emcalCorrComponent ){
        fRecoUtils        = emcalCorrComponent->GetRecoUtils();
      } else {
        emcalCorrComponent = emcalCorrTask->GetCorrectionComponent("AliEmcalCorrectionCellBadChannel");
        if( emcalCorrComponent ){
          fRecoUtils        = emcalCorrComponent->GetRecoUtils();
        }
      }
    }
    
    if (fRecoUtils) {
      fEMCALInitialized = kTRUE;
      fRecoUtils->SetNumberOfCellsFromEMCALBorder(0);
      return;
    }
}

///
/// Fill the invariant mass analysis per channel with the
/// corrected clusters, and other general histograms.
//__________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::FillHistograms() {

  Int_t absId1   = -1;
  Int_t iSupMod1 = -1;
  Int_t iphi1    = -1;
  Int_t ieta1    = -1;
  Int_t absId2   = -1;
  Int_t iSupMod2 = -1;
  Int_t iphi2    = -1;
  Int_t ieta2    = -1;
  Bool_t shared  = kFALSE;
  
  Int_t bc  = InputEvent()->GetBunchCrossNumber();
  
  for(Int_t iClu=0; iClu<fCaloClustersArr->GetEntriesFast()-1; iClu++) {
    
    AliVCluster *c1 = (AliVCluster *) fCaloClustersArr->At(iClu);
    
    if (!AcceptCluster(c1)) continue;
    
    fRecoUtils->GetMaxEnergyCell(fEMCALGeo, fEMCALCells,c1,absId1,iSupMod1,ieta1,iphi1,shared); //RECOUTILS

    if( isEMC && !(isDMC) &&  iSupMod1 > 11 ) continue;
    if( isDMC && !(isEMC) &&  iSupMod1 < 12 ) continue;
    
    c1->GetMomentum(fMomentum1,fVertex);
    
    // Check if cluster is in fidutial region, not too close to borders
    if(fEMCALGeo == NULL || fEMCALCells == NULL) continue;
    Bool_t in1 = kTRUE;
    // fRecoUtils->CheckCellFiducialRegion(fEMCALGeo, c1, fEMCALCells); //RECOUTILS
    
    // Clusters not facing frame structures
    Bool_t mask1 = MaskFrameCluster(iSupMod1, ieta1);
    
    Double_t time1 = c1->GetTOF()*1.e9;
    
    fhClusterTime            ->Fill(c1->E(),time1);
    fhClusterTimeSM[iSupMod1]->Fill(c1->E(),time1);
    
    // Combine cluster with other clusters and get the invariant mass
    for (Int_t jClu=iClu+1; jClu < fCaloClustersArr->GetEntriesFast(); jClu++) {
      AliAODCaloCluster *c2 = (AliAODCaloCluster *) fCaloClustersArr->At(jClu);
    
      if (!AcceptCluster(c2)) continue;
      
      fRecoUtils->GetMaxEnergyCell(fEMCALGeo, fEMCALCells,c2,absId2,iSupMod2,ieta2,iphi2,shared); //RECOUTILS

      if( isEMC && !(isDMC) &&  iSupMod2 > 11 ) continue;
      if( isDMC && !(isEMC) &&  iSupMod2 < 12 ) continue;
      
      c2->GetMomentum(fMomentum2,fVertex);
      
      fMomentum12 = fMomentum1+fMomentum2;
      Float_t invmass = fMomentum12.M()*1000; 
      
      //Asimetry cut      
      Float_t asym = TMath::Abs(fMomentum1.E()-fMomentum2.E())/(fMomentum1.E()+fMomentum2.E());
      if(asym > fAsyCut) continue;
      
      //Time cut
      Double_t time2 = c2->GetTOF()*1.e9;
      
      // fhClusterPairDiffTime->Fill(fMomentum12.E(),time1-time2);
      if(TMath::Abs(time1-time2) > fDTimeCut) continue;

      //Opening angle cut to reject too close clusters
      if( ((fMomentum1.Angle(fMomentum2.Vect()) < fOpAnglemin) || (fMomentum1.Angle(fMomentum2.Vect()) > fOpAnglemax)))  continue;
      
      if(invmass < fMaxBin && invmass > fMinBin ) {
        //Check if cluster is in fidutial region, not too close to borders
        Bool_t in2 = kTRUE;
        // fRecoUtils->CheckCellFiducialRegion(fEMCALGeo, c2, fEMCALCells); //RECOUTILS
        
        // Clusters not facing frame structures
        Bool_t mask2 = MaskFrameCluster(iSupMod2, ieta2);         
        
        if(in1 && in2) { //acceptance cuts
          fHmgg->Fill(invmass,fMomentum12.Pt()); 
          
          if(iSupMod1==iSupMod2) {
            fHmggSM[iSupMod1]->Fill(invmass,fMomentum12.Pt()); 
            // fhClusterPairDiffTimeSameSM[iSupMod1]->Fill(fMomentum12.E(),time1-time2);
          } else {           
            fHmggDifferentSM ->Fill(invmass,fMomentum12.Pt());
          }        
          
          if(!mask1 && !mask2) {  // Pair not facing frame
            fHmggMaskFrame->Fill(invmass,fMomentum12.Pt()); 
            
            if(iSupMod1==iSupMod2) fHmggSMMaskFrame[iSupMod1]->Fill(invmass,fMomentum12.Pt()); 
            else                   fHmggDifferentSMMaskFrame ->Fill(invmass,fMomentum12.Pt());
          }
          
          if(invmass > fInvMassCutMin && invmass < fInvMassCutMax) { //restrict to clusters really close to pi0 peak (pair in 100 < mass < 160)

            // Check time of cells in both clusters, and fill time histogram
            for(Int_t icell = 0; icell < c1->GetNCells(); icell++) {
              Int_t absID = c1->GetCellAbsId(icell);   
              fHTpi0[bc%4]->Fill(absID, fEMCALCells->GetCellTime(absID)*1.e9);  
            }
            
            for(Int_t icell = 0; icell < c2->GetNCells(); icell++) {
              Int_t absID = c2->GetCellAbsId(icell);   
              fHTpi0[bc%4]->Fill(absID, fEMCALCells->GetCellTime(absID)*1.e9);  
            }
            
            //Opening angle of 2 photons
            Float_t opangle = fMomentum1.Angle(fMomentum2.Vect())*TMath::RadToDeg();
            
            fHOpeningAngle ->Fill(opangle,fMomentum12.Pt()); 
            fHAsymmetry    ->Fill(asym,fMomentum12.Pt()); 
            
            if(iSupMod1==iSupMod2) {
              fHOpeningAngleSM[iSupMod1] ->Fill(opangle,fMomentum12.Pt());
              fHAsymmetrySM[iSupMod1]    ->Fill(asym,fMomentum12.Pt());
            } else {      
              fHOpeningAngleDifferentSM  ->Fill(opangle,fMomentum12.Pt());
              fHAsymmetryDifferentSM     ->Fill(asym,fMomentum12.Pt());
            }
            
          }
          
        }
        
        //In case of filling only channels with second cluster in same SM
        if(fSameSM && iSupMod1!=iSupMod2) continue;
        
        if (fGroupNCells == 0) {
          fHmpi0[iSupMod1][ieta1][iphi1]->Fill(invmass);
          fHmpi0[iSupMod2][ieta2][iphi2]->Fill(invmass);
          
          if(invmass > fInvMassCutMin && invmass < fInvMassCutMax) { //restrict to clusters really close to pi0 peak
            fhTowerDecayPhotonHit      [iSupMod1]->Fill(ieta1,iphi1);
            fhTowerDecayPhotonEnergy   [iSupMod1]->Fill(ieta1,iphi1,fMomentum1.E());
            fhTowerDecayPhotonAsymmetry[iSupMod1]->Fill(ieta1,iphi1,asym);
            
            fhTowerDecayPhotonHit      [iSupMod2]->Fill(ieta2,iphi2);
            fhTowerDecayPhotonEnergy   [iSupMod2]->Fill(ieta2,iphi2,fMomentum2.E());
            fhTowerDecayPhotonAsymmetry[iSupMod2]->Fill(ieta2,iphi2,asym);
            
            if(!mask1)fhTowerDecayPhotonHitMaskFrame[iSupMod1]->Fill(ieta1,iphi1);
            if(!mask2)fhTowerDecayPhotonHitMaskFrame[iSupMod2]->Fill(ieta2,iphi2);
          }
        }	else { //group cells
          //printf("Regroup N %d, eta1 %d, phi1 %d, eta2 %d, phi2 %d \n",fGroupNCells, ieta1, iphi1, ieta2, iphi2);
          for (Int_t i = -fGroupNCells; i < fGroupNCells+1; i++) {
            for (Int_t j = -fGroupNCells; j < fGroupNCells+1; j++) {              
              Int_t absId11 = fEMCALGeo->GetAbsCellIdFromCellIndexes(iSupMod1, iphi1+j, ieta1+i); 
              Int_t absId22 = fEMCALGeo->GetAbsCellIdFromCellIndexes(iSupMod2, iphi2+j, ieta2+i); 
              
              Bool_t ok1 = kFALSE;
              Bool_t ok2 = kFALSE;
              
              for(Int_t icell = 0; icell < c1->GetNCells(); icell++) {
                if(c1->GetCellsAbsId()[icell] == absId11) ok1=kTRUE;
              }
              
              for(Int_t icell = 0; icell < c2->GetNCells(); icell++) {
                if(c2->GetCellsAbsId()[icell] == absId22) ok2=kTRUE;
              }
              
              if(ok1 && (ieta1+i >= 0) && (iphi1+j >= 0) && (ieta1+i < 48) && (iphi1+j < 24)) {
                fHmpi0[iSupMod1][ieta1+i][iphi1+j]->Fill(invmass);
              }
                
              if(ok2 && (ieta2+i >= 0) && (iphi2+j >= 0) && (ieta2+i < 48) && (iphi2+j < 24)) {
                fHmpi0[iSupMod2][ieta2+i][iphi2+j]->Fill(invmass);
              }
            }
          }
        }
        
      }
    }
  } // end of loop over EMCAL clusters
  
}


///
/// Init geometry and set the geometry matrix, for the first event, skip the rest.
/// Also set once the run dependent calibrations.
//________________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::InitGeometryMatrices() {
  Int_t runnumber = InputEvent()->GetRunNumber() ;
  
  //
  // Load default geo matrices if requested
  if(fImportGeometryFromFile && !gGeoManager) {
    if(fImportGeometryFilePath=="") { // If not specified, set location depending on run number
      // "$ALICE_ROOT/EVE/alice-data/default_geo.root"
      if     (runnumber <  140000) fImportGeometryFilePath = AliDataFile::GetFileNameOADB("EMCAL/geometry_2010.root").data();
      else if(runnumber <  171000) fImportGeometryFilePath = AliDataFile::GetFileNameOADB("EMCAL/geometry_2011.root").data();
      else if(runnumber <  198000) fImportGeometryFilePath = AliDataFile::GetFileNameOADB("EMCAL/geometry_2012.root").data(); // 2012-2013
      else                         fImportGeometryFilePath = AliDataFile::GetFileNameOADB("EMCAL/geometry_2015.root").data(); // >= 2015
    }
    AliInfo(Form("Import %s",fImportGeometryFilePath.Data()));
    TGeoManager::Import(fImportGeometryFilePath) ; // default need file "geometry.root" in local dir!!!!
  }
  
  //
  if(fLoadMatrices) { // Load matrices
    AliInfo("Load user defined EMCAL geometry matrices");
    // OADB if available
    AliOADBContainer emcGeoMat("AliEMCALgeo");

    if(fOADBFilePath!="")
      emcGeoMat.InitFromFile(Form("%s/EMCALlocal2master.root",fOADBFilePath.Data()),"AliEMCALgeo");
    else
      emcGeoMat.InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALlocal2master.root").data(),"AliEMCALgeo");

    TObjArray *matEMCAL=(TObjArray*)emcGeoMat.GetObject(runnumber,"EmcalMatrices");
    
    for(Int_t mod = 0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++) { //SM loop
      
      if (!fMatrix[mod]) { // Get it from OADB
        AliDebug(1,Form("EMCAL matrices SM %d, %p",mod,((TGeoHMatrix*) matEMCAL->At(mod))));
        fMatrix[mod] = (TGeoHMatrix*) matEMCAL->At(mod) ;
      }        
      
      if(fMatrix[mod]) {  
        if(DebugLevel() > 1) 
          fMatrix[mod]->Print();
        fEMCALGeo->SetMisalMatrix(fMatrix[mod],mod) ;  
      } else if(gGeoManager) {
        AliWarning(Form("Set matrix for SM %d from gGeoManager",mod));
        fEMCALGeo->SetMisalMatrix(fEMCALGeo->GetMatrixForSuperModuleFromGeoManager(mod),mod) ;
      } else {
        AliError(Form("Alignment matrix for SM %d is not available",mod));
      }
    }
  } else if(!gGeoManager) { // Load matrices from Data
    AliInfo("Get geo matrices from data");
    //Still not implemented in AOD, just a workaround to be able to work at least with ESDs	
    if(!strcmp(InputEvent()->GetName(),"AliAODEvent")) {
      AliWarning("Use ideal geometry, values geometry matrix not kept in AODs");
    } else {// ESD
      AliDebug(1,"AliAnalysisTaskEMCALClusterize Load Misaligned matrices");
      
      for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++) {
        if(InputEvent()->GetEMCALMatrix(mod)) {
          if(DebugLevel() > 1) 
            InputEvent()->GetEMCALMatrix(mod)->Print();
          fEMCALGeo->SetMisalMatrix(InputEvent()->GetEMCALMatrix(mod),mod) ;
        } 
      }
    }
  } else if(gGeoManager) { // Load default gGeoManager matrices
    for(Int_t mod = 0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++) {
      AliWarning(Form("Set matrix for SM %d from gGeoManager",mod));
      fEMCALGeo->SetMisalMatrix(fEMCALGeo->GetMatrixForSuperModuleFromGeoManager(mod),mod) ;
    }
  }

}


///
/// Create output container, init geometry.
//___________________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::UserCreateOutputObjects() {

  if(!fEMCALGeo)fEMCALGeo =  AliEMCALGeometry::GetInstance(fEMCALGeoName) ;
  Int_t nSM = (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules();
  
  fOutputContainer = new TList();
  const Int_t buffersize = 255;
  char hname[buffersize], htitl[buffersize];
  
  fhNEvents        = new TH1I("hNEvents", "Number of analyzed events"   , 1 , 0 , 1  ) ;
  fOutputContainer->Add(fhNEvents);
  
  if ( fCheckCentrality ) {
    fhCentrality   = new TH1F ("hCentrality",
                                Form("Number of events in centrality bins, |vz|<10 cm, method <%s> ",fCentralityClass.Data()),
                                1000,-1000.,1000.) ;
    fhCentrality->SetXTitle("Centrality bin");
    fOutputContainer->Add(fhCentrality) ;  
    
    fhCentralitySelected   = new TH1F ("hCentralitySelected",
                                        Form("Number of selected events in centrality bin, |vz|<10 cm, method <%s> ",fCentralityClass.Data()),
                                        100,0.,100.) ;
    fhCentralitySelected->SetXTitle("Centrality bin");
    fOutputContainer->Add(fhCentralitySelected) ;  
  }

  fHmgg = new TH2F("hmgg","2-cluster invariant mass",fNbins,fMinBin,fMaxBin,100,0,10);
  fHmgg->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
  fHmgg->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
  fOutputContainer->Add(fHmgg);
  
  fHmggDifferentSM = new TH2F("hmggDifferentSM","2-cluster invariant mass, different SM",fNbins,fMinBin,fMaxBin,100,0,10);
  fHmggDifferentSM->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
  fHmggDifferentSM->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
  fOutputContainer->Add(fHmggDifferentSM);
  
  fHOpeningAngle = new TH2F("hopang","2-cluster opening angle",100,0.,50.,100,0,10);
  fHOpeningAngle->SetXTitle("#alpha_{#gamma #gamma}");
  fHOpeningAngle->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
  fOutputContainer->Add(fHOpeningAngle);
  
  fHOpeningAngleDifferentSM = new TH2F("hopangDifferentSM","2-cluster opening angle, different SM",100,0,50.,100,0,10);
  fHOpeningAngleDifferentSM->SetXTitle("#alpha_{#gamma #gamma}");
  fHOpeningAngleDifferentSM->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
  fOutputContainer->Add(fHOpeningAngleDifferentSM);
   
  fHAsymmetry = new TH2F("hasym","2-cluster opening angle",100,0.,1.,100,0,10);
  fHAsymmetry->SetXTitle("a");
  fHAsymmetry->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
  fOutputContainer->Add(fHAsymmetry);
  
  fHAsymmetryDifferentSM = new TH2F("hasymDifferentSM","2-cluster opening angle, different SM",100,0,1.,100,0,10);
  fHAsymmetryDifferentSM->SetXTitle("a");
  fHAsymmetryDifferentSM->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
  fOutputContainer->Add(fHAsymmetryDifferentSM);
  
  fHmggMaskFrame = new TH2F("hmggMaskFrame","2-cluster invariant mass, frame masked",fNbins,fMinBin,fMaxBin,100,0,10);
  fHmggMaskFrame->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
  fHmggMaskFrame->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
  fOutputContainer->Add(fHmggMaskFrame);
  
  fHmggDifferentSMMaskFrame = new TH2F("hmggDifferentSMMaskFrame","2-cluster invariant mass, different SM, frame masked",
                                       fNbins,fMinBin,fMaxBin,100,0,10);
  fHmggDifferentSMMaskFrame->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
  fHmggDifferentSMMaskFrame->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
  fOutputContainer->Add(fHmggDifferentSMMaskFrame);
  
  for(Int_t iSM = 0; iSM < nSM; iSM++) {
    snprintf(hname, buffersize, "hmgg_SM%d",iSM);
    snprintf(htitl, buffersize, "Two-gamma inv. mass for super mod %d",iSM);
    fHmggSM[iSM] = new TH2F(hname,htitl,fNbins,fMinBin,fMaxBin,100,0,10);
    fHmggSM[iSM]->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
    fHmggSM[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHmggSM[iSM]);
    
    snprintf(hname, buffersize, "hmgg_SM%d_MaskFrame",iSM);
    snprintf(htitl, buffersize, "Two-gamma inv. mass for super mod %d",iSM);
    fHmggSMMaskFrame[iSM] = new TH2F(hname,htitl,fNbins,fMinBin,fMaxBin,100,0,10);
    fHmggSMMaskFrame[iSM]->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
    fHmggSMMaskFrame[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHmggSMMaskFrame[iSM]);
    
    snprintf(hname, buffersize, "hopang_SM%d",iSM);
    snprintf(htitl, buffersize, "Opening angle for super mod %d",iSM);
    fHOpeningAngleSM[iSM] = new TH2F(hname,htitl,100,0.,50.,100,0,10);
    fHOpeningAngleSM[iSM]->SetXTitle("#alpha_{#gamma #gamma} (deg)");
    fHOpeningAngleSM[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHOpeningAngleSM[iSM]);
    
    snprintf(hname,buffersize, "hopang_PairSM%d",iSM);
    snprintf(htitl,buffersize, "Opening angle for SM pair: %d",iSM);
    fHOpeningAnglePairSM[iSM] = new TH2F(hname,htitl,100,0.,50.,100,0,10);
    fHOpeningAnglePairSM[iSM]->SetXTitle("#alpha_{#gamma #gamma} (deg)");
    fHOpeningAnglePairSM[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHOpeningAnglePairSM[iSM]);    
    
    snprintf(hname, buffersize, "hasym_SM%d",iSM);
    snprintf(htitl, buffersize, "Asymmetry for super mod %d",iSM);
    fHAsymmetrySM[iSM] = new TH2F(hname,htitl,100,0.,1.,100,0,10);
    fHAsymmetrySM[iSM]->SetXTitle("a");
    fHAsymmetrySM[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHAsymmetrySM[iSM]);
    
    snprintf(hname,buffersize, "hasym_PairSM%d",iSM);
    snprintf(htitl,buffersize, "Asymmetry for SM pair: %d",iSM);
    fHAsymmetryPairSM[iSM] = new TH2F(hname,htitl,100,0.,1.,100,0,10);
    fHAsymmetryPairSM[iSM]->SetXTitle("a");
    fHAsymmetryPairSM[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHAsymmetryPairSM[iSM]);    
    
    Int_t colmax = 48;
    Int_t rowmax = 24;
    
    fhTowerDecayPhotonHit[iSM]  = new TH2F (Form("hTowerDecPhotonHit_Mod%d",iSM),
                                            Form("Entries in grid of cells in Module %d",iSM), 
                                            colmax+2,-1.5,colmax+0.5, rowmax+2,-1.5,rowmax+0.5); 
    fhTowerDecayPhotonHit[iSM]->SetYTitle("row (phi direction)");
    fhTowerDecayPhotonHit[iSM]->SetXTitle("column (eta direction)");
    fOutputContainer->Add(fhTowerDecayPhotonHit[iSM]);
    
    fhTowerDecayPhotonEnergy[iSM]  = new TH2F (Form("hTowerDecPhotonEnergy_Mod%d",iSM),
                                               Form("Accumulated energy in grid of cells in Module %d",iSM), 
                                               colmax+2,-1.5,colmax+0.5, rowmax+2,-1.5,rowmax+0.5); 
    fhTowerDecayPhotonEnergy[iSM]->SetYTitle("row (phi direction)");
    fhTowerDecayPhotonEnergy[iSM]->SetXTitle("column (eta direction)");
    fOutputContainer->Add(fhTowerDecayPhotonEnergy[iSM]);
    
    fhTowerDecayPhotonAsymmetry[iSM]  = new TH2F (Form("hTowerDecPhotonAsymmetry_Mod%d",iSM),
                                                  Form("Accumulated asymmetry in grid of cells in Module %d",iSM), 
                                                  colmax+2,-1.5,colmax+0.5, rowmax+2,-1.5,rowmax+0.5); 
    fhTowerDecayPhotonAsymmetry[iSM]->SetYTitle("row (phi direction)");
    fhTowerDecayPhotonAsymmetry[iSM]->SetXTitle("column (eta direction)");
    fOutputContainer->Add(fhTowerDecayPhotonAsymmetry[iSM]);
    
    fhTowerDecayPhotonHitMaskFrame[iSM]  = new TH2F (Form("hTowerDecPhotonHit_Mod%d_MaskFrame",iSM),Form("Entries in grid of cells in Module %d",iSM), 
                                                     colmax+2,-1.5,colmax+0.5, rowmax+2,-1.5,rowmax+0.5); 
    fhTowerDecayPhotonHitMaskFrame[iSM]->SetYTitle("row (phi direction)");
    fhTowerDecayPhotonHitMaskFrame[iSM]->SetXTitle("column (eta direction)");
    fOutputContainer->Add(fhTowerDecayPhotonHitMaskFrame[iSM]);
    
    fhClusterTimeSM[iSM] = new TH2F(Form("hClusterTime_SM%d",iSM),"cluster time vs E",100,0,10, 200,-1000,1000);
    fhClusterTimeSM[iSM]->SetXTitle("E (GeV)");
    fhClusterTimeSM[iSM]->SetYTitle("t (ns)");
    fOutputContainer->Add(fhClusterTimeSM[iSM]);
  }  
  
  Int_t nchannels = nSM*AliEMCALGeoParams::fgkEMCALRows*AliEMCALGeoParams::fgkEMCALCols;
  for(Int_t ibc = 0; ibc < 4; ibc++) {
    fHTpi0[ibc] = new TH2F(Form("hTime_BC%d",ibc),Form("Time of cell clusters under pi0 peak, bunch crossing %d",ibc),
                           nchannels,0,nchannels, fNTimeBins,fMinTimeBin,fMaxTimeBin);
    fOutputContainer->Add(fHTpi0[ibc]);       
    fHTpi0[ibc]->SetYTitle("time (ns)");
    fHTpi0[ibc]->SetXTitle("abs. Id. ");
  }
  
  fhClusterTime = new TH2F("hClusterTime","cluster time vs E",100,0,10, 100,0,1000);
  fhClusterTime->SetXTitle("E (GeV)");
  fhClusterTime->SetYTitle("t (ns)");
  fOutputContainer->Add(fhClusterTime);
  
  for(Int_t iMod=0; iMod < nSM; iMod++) {
    for(Int_t iRow=0; iRow < AliEMCALGeoParams::fgkEMCALRows; iRow++) {
      for(Int_t iCol=0; iCol < AliEMCALGeoParams::fgkEMCALCols; iCol++) {
        snprintf(hname,buffersize, "%d_%d_%d",iMod,iCol,iRow);
        snprintf(htitl,buffersize, "Two-gamma inv. mass for super mod %d, cell(col,row)=(%d,%d)",iMod,iCol,iRow);
        fHmpi0[iMod][iCol][iRow] = new TH1F(hname,htitl,fNbins,fMinBin,fMaxBin);
        fHmpi0[iMod][iCol][iRow]->SetXTitle("mass (MeV/c^{2})");
        fOutputContainer->Add(fHmpi0[iMod][iCol][iRow]);
      }
    }
  }
  
  fOutputContainer->SetOwner(kTRUE);
    
  PostData(1,fOutputContainer);
}

///
/// Check if cell is in one of the regions where we have significant amount of material in front of EMCAL.
/// \return True if this cell is in one problematic region
/// \param iSM: supermodule number of the cell.
/// \param ieta: column index of the cell.
//______________________________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPi0CalibSelectionV2::MaskFrameCluster(Int_t iSM, Int_t ieta) const
{
  Int_t icol = ieta;
  if(iSM%2) icol+=48; // Impair SM, shift index [0-47] to [48-96]
  
  if (fNMaskCellColumns && fMaskCellColumns) {
    for (Int_t imask = 0; imask < fNMaskCellColumns; imask++) {
      if(icol==fMaskCellColumns[imask]) return kTRUE;
    }
  }
  
  return kFALSE;
}


///
/// Check if trigger selected
///
//____________________________________________________________________
Bool_t AliAnalysisTaskEMCALPi0CalibSelectionV2::IsTriggerSelected(AliVEvent *event){

  if(fTriggerName!="" && !isMC ) {
      AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()); 
      if (fInputHandler==NULL) return kFALSE;

      if( fInputHandler->GetEventSelection() ) {
        if( fTriggerName == "INT7" ){
          fOfflineTriggerMask = AliVEvent::kINT7;
          if( !(fOfflineTriggerMask) ) return kFALSE;
          isEMC = kTRUE; 
          isDMC = kTRUE;
        } else if ( fTriggerName == "EMC7" ) {
          fOfflineTriggerMask = AliVEvent::kEMC7;
          if( !(fOfflineTriggerMask ) ) return kFALSE;

          TString firedTrigClass = event->GetFiredTriggerClasses();
          std::cout << "Trigger: " << firedTrigClass << std::endl;

          if( firedTrigClass.Contains("CEMC7") ) {
            isEMC = kTRUE;
            if( firedTrigClass.Contains("CDMC7") ){
              isDMC = kTRUE;
            }
          } else if( firedTrigClass.Contains("CDMC7") ){
            isDMC = kTRUE;
          }
        }
        return kTRUE;
      }

  } else if (isMC){
    AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()); 
      if (fInputHandler==NULL) return kFALSE;
      if( fInputHandler->GetEventSelection() ) {
        fOfflineTriggerMask = AliVEvent::kAny;
        if( !(fOfflineTriggerMask) ) return kFALSE;
        isEMC = kTRUE;
        isDMC = kTRUE;
      }
      return kTRUE;
  } else {
    std::cout << "No trigger selected" << std::endl;
    return kFALSE;
  }

}

///
/// Main method, do the analysis per event:
/// * first, select the events;
/// * then, correct the clusters if needed;
/// * finally, fill the histograms per channel after recalibration.
//__________________________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::UserExec(Option_t* /* option */) {

  // Get the input event

  AliVEvent* event = 0;
  event = InputEvent();
  
  if(!event) {
    AliWarning("Input event not available!");
    return;
  }

  // Acccess once the geometry matrix and temperature corrections and calibration coefficients
  if(fhNEvents->GetEntries() == 1) {
    InitGeometryMatrices();
  }

  // Event selection

  isEMC=kFALSE;
  isDMC=kFALSE;

  if( !IsTriggerSelected(event) ) return;
  if( !(isEMC) && !(isDMC) ) return;

  if( !fEMCALInitialized ) {
    InitializeEMCAL( event );
    // if( !fEMCALInitialized ) return;
  }

  // Centrality selection
  
  if ( fCheckCentrality ) {
    AliMultSelection* multSelection = (AliMultSelection * ) event->FindListObject("MultSelection") ;
    if ( multSelection ) {
      Float_t cent = multSelection->GetMultiplicityPercentile(fCentralityClass, fCentWithEventSel);
      fhCentrality->Fill(cent);
      
      AliDebug(1,Form("Centrality %2.1f for class <%s>, event sel %d\n",cent,fCentralityClass.Data(),fCentWithEventSel));
      //printf("Centrality %2.1f for class <%s>, select in [%1.1f,%1.1f]\n",cent,fCentralityClass.Data(),fCentMin,fCentMax);
      
      if ( cent < fCentMin || cent >= fCentMax ) return ;
      
      fhCentralitySelected->Fill(cent);
    }
  }
  
  AliDebug(1,Form("<<< %s: Event %d >>>",event->GetName(), (Int_t)Entry()));
  
  // Get the primary vertex
    
  event->GetPrimaryVertex()->GetXYZ(fVertex) ;
  
  AliDebug(1,Form("Vertex: (%.3f,%.3f,%.3f)",fVertex[0],fVertex[1],fVertex[2]));
  
  fhNEvents->Fill(0); //Count the events to be analyzed

  //Get the list of clusters and cells
  fEMCALCells       = event->GetEMCALCells();

  fCaloClustersArr  = new TRefArray();
  event->GetEMCALClusters(fCaloClustersArr);
  
  AliDebug(1,Form("N CaloClusters: %d - N CaloCells %d",fCaloClustersArr->GetEntriesFast(), fEMCALCells->GetNumberOfCells()));

  FillHistograms();  

  delete fCaloClustersArr;
  
  PostData(1,fOutputContainer);
}

///
/// Print settings.
//_____________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::PrintInfo()
{
  printf("Cluster cuts: %2.2f < E < %2.2f GeV; number of cells > %d; Assymetry < %1.2f, pair time diff < %2.2f, %2.2f < t < %2.2f ns\n", 
         fEmin,fEmax, fMinNCells, fAsyCut, fDTimeCut,fTimeMin,fTimeMax) ;
  
  printf("Group %d cells\n", fGroupNCells) ;
    
  printf("Histograms: bins %d; energy range: %2.2f < E < %2.2f MeV\n",fNbins,fMinBin,fMaxBin) ;
  
  printf("OADB path        : %s\n",fOADBFilePath .Data());  
  printf("EMCAL Geometry name: < %s >, Load Matrices %d\n",fEMCALGeoName.Data(), fLoadMatrices) ;

  if ( fCheckCentrality )
    printf("Select centrality bin [%1.1f,%1.1f] for class <%s>, event sel. %d\n",fCentMin,fCentMax,fCentralityClass.Data(),fCentWithEventSel);
  
  if(fLoadMatrices) { 
    for(Int_t ism = 0; ism < AliEMCALGeoParams::fgkEMCALModules; ism++) {
      if(fMatrix[ism]) fMatrix[ism]->Print(); 
      }
  }
}

///
/// Set the total number of columns to be masked in the analysis
/// \param n: number of columns
//_____________________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::SetNMaskCellColumns(Int_t n)
{
    if(n > fNMaskCellColumns) {
        delete [] fMaskCellColumns ;
        fMaskCellColumns = new Int_t[n] ;
    }
    
    fNMaskCellColumns = n ;
}

///
/// Define which column must be masked and its position in the array of columns
/// \param ipos: position in the list of columns.
/// \param icol: column to be masked in EMCal mapping coordinate.
//___________________________________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::SetMaskCellColumn(Int_t ipos, Int_t icol)
{
    if(ipos < fNMaskCellColumns) fMaskCellColumns[ipos] = icol            ;
    else AliWarning("Mask column not set, position larger than allocated set size first") ;
}


///
/// Create cuts/param objects and publish to slot. Comment out for the moment.
//______________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::Terminate(Option_t*)
{
  AliDebug(1,"Not implemented");

}

