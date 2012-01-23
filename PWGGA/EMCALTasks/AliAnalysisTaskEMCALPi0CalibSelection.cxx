/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Boris Polishchuk                                               *
 * Adapted to AOD reading by Gustavo Conesa                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//---------------------------------------------------------------------------// 
//                                                                           //
// Fill histograms (one per cell) with two-cluster invariant mass            //
// using calibration coefficients of the previous iteration.                 //
// Histogram for a given cell is filled if the most energy of one cluster    //
// is deposited in this cell and the other cluster could be anywherein EMCAL.//
//                                                                           //
//---------------------------------------------------------------------------//

//#include <cstdlib>
//#include <Riostream.h>
// Root 
#include "TLorentzVector.h"
#include "TRefArray.h"
#include "TList.h"
#include "TH1F.h"
#include <TGeoManager.h>

// AliRoot
#include "AliAnalysisTaskEMCALPi0CalibSelection.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliEMCALGeometry.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliEMCALRecoUtils.h"

ClassImp(AliAnalysisTaskEMCALPi0CalibSelection)


//__________________________________________________
AliAnalysisTaskEMCALPi0CalibSelection::AliAnalysisTaskEMCALPi0CalibSelection(const char* name) :
AliAnalysisTaskSE(name),fEMCALGeo(0x0), 
fEmin(0.5),               fEmax(15.),      
fL0min(0.01),             fL0max(0.5),              fDTimeCut(100.), 
fAsyCut(1.),              fMinNCells(2),            fGroupNCells(0),
fLogWeight(4.5),          fSameSM(kFALSE),          fFilteredInput(kFALSE),
fCorrectClusters(kFALSE), fEMCALGeoName("EMCAL_COMPLETEV1"), 
fTriggerName("EMC"),      fRecoUtils(new AliEMCALRecoUtils), 
fCuts(0x0),               fLoadMatrices(0),
fNMaskCellColumns(11),    fMaskCellColumns(0x0),
fInvMassCutMin(110.),     fInvMassCutMax(160.),
//Histograms
fOutputContainer(0x0),    fNbins(300),              
fMinBin(0.),              fMaxBin(300.),   
fNTimeBins(1000),         fMinTimeBin(0.),          fMaxTimeBin(1000.),   
fHmgg(0x0),               fHmggDifferentSM(0x0), 
fHmggMaskFrame(0x0),      fHmggDifferentSMMaskFrame(0x0), 
fHOpeningAngle(0x0),      fHOpeningAngleDifferentSM(0x0),  
fHIncidentAngle(0x0),     fHIncidentAngleDifferentSM(0x0),
fHAsymmetry(0x0),         fHAsymmetryDifferentSM(0x0),  
fhNEvents(0x0),
fhClusterTime(0x0),       fhClusterPairDiffTime(0x0)
{
  //Named constructor which should be used.
  
  for(Int_t iMod=0; iMod < AliEMCALGeoParams::fgkEMCALModules; iMod++) {
    for(Int_t iX=0; iX<24; iX++) {
      for(Int_t iZ=0; iZ<48; iZ++) {
        fHmpi0[iMod][iZ][iX]   = 0 ;
      }
    } 
  }
  
  fHTpi0[0]= 0 ;
  fHTpi0[1]= 0 ;
  fHTpi0[2]= 0 ;
  fHTpi0[3]= 0 ;
  
  
  fMaskCellColumns = new Int_t[fNMaskCellColumns];
  fMaskCellColumns[0] = 6 ;  fMaskCellColumns[1] = 7 ;  fMaskCellColumns[2] = 8 ; 
  fMaskCellColumns[3] = 35;  fMaskCellColumns[4] = 36;  fMaskCellColumns[5] = 37; 
  fMaskCellColumns[6] = 12+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[7] = 13+AliEMCALGeoParams::fgkEMCALCols;
  fMaskCellColumns[8] = 40+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[9] = 41+AliEMCALGeoParams::fgkEMCALCols; 
  fMaskCellColumns[10]= 42+AliEMCALGeoParams::fgkEMCALCols; 
  
  for(Int_t iSMPair = 0; iSMPair < AliEMCALGeoParams::fgkEMCALModules/2; iSMPair++) {
    fHmggPairSameSectorSM[iSMPair]          = 0;
    fHmggPairSameSectorSMMaskFrame[iSMPair] = 0;
    fhClusterPairDiffTimeSameSector[iSMPair]= 0;
  }
  for(Int_t iSMPair = 0; iSMPair < AliEMCALGeoParams::fgkEMCALModules-2; iSMPair++){ 
    fHmggPairSameSideSM[iSMPair]            = 0;
    fHmggPairSameSideSMMaskFrame[iSMPair]   = 0;
    fhClusterPairDiffTimeSameSide[iSMPair]  = 0;
  }
  
  for(Int_t iSM = 0; iSM < AliEMCALGeoParams::fgkEMCALModules; iSM++) {
    fHmggSM[iSM]                     = 0;
    fHmggSMMaskFrame[iSM]            = 0;
    fHOpeningAngleSM[iSM]            = 0;
    fHOpeningAnglePairSM[iSM]        = 0;
    fHAsymmetrySM[iSM]               = 0;
    fHAsymmetryPairSM[iSM]           = 0;
    fHIncidentAngleSM[iSM]           = 0;
    fHIncidentAnglePairSM[iSM]       = 0;
    fhTowerDecayPhotonHit[iSM]       = 0;
    fhTowerDecayPhotonEnergy[iSM]    = 0;
    fhTowerDecayPhotonAsymmetry[iSM] = 0;
    fhTowerDecayPhotonHitMaskFrame[iSM]= 0;
    fMatrix[iSM]                     = 0x0;
    fhClusterTimeSM[iSM]             = 0;
    fhClusterPairDiffTimeSameSM[iSM] = 0;
  }
  
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());  // will contain cuts or local params
  
}

//__________________________________________________
AliAnalysisTaskEMCALPi0CalibSelection::~AliAnalysisTaskEMCALPi0CalibSelection()
{
  //Destructor.
  
  if(fOutputContainer){
    fOutputContainer->Delete() ; 
    delete fOutputContainer ;
  }
  
  if(fEMCALGeo)         delete fEMCALGeo  ;
  if(fRecoUtils)        delete fRecoUtils ;
  if(fNMaskCellColumns) delete [] fMaskCellColumns;
  
}

//_____________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::LocalInit()
{
  // Local Initialization
  
  // Create cuts/param objects and publish to slot
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  fCuts = new TList();
  
  snprintf(onePar,buffersize, "Custer cuts: %2.2f < E < %2.2f GeV; %2.2f < Lambda0_2 < %2.2f GeV; min number of cells %d; Assymetry cut %1.2f, time1-time2 < %2.2f; %3.1f < Mass < %3.1f", 
           fEmin,fEmax, fL0min, fL0max, fMinNCells, fAsyCut, fDTimeCut, fInvMassCutMin, fInvMassCutMax) ;
  fCuts->Add(new TObjString(onePar));
  snprintf(onePar,buffersize, "Group %d cells;", fGroupNCells) ;
  fCuts->Add(new TObjString(onePar));
  snprintf(onePar,buffersize, "Cluster maximal cell away from border at least %d cells;", fRecoUtils->GetNumberOfCellsFromEMCALBorder()) ;
  fCuts->Add(new TObjString(onePar));
  snprintf(onePar,buffersize, "Histograms, Mass bins %d; energy range: %2.2f < E < %2.2f GeV;",fNbins,fMinBin,fMaxBin) ;
  fCuts->Add(new TObjString(onePar));
  snprintf(onePar,buffersize, "Histograms, Time bins %d; energy range: %2.2f < E < %2.2f GeV;",fNTimeBins,fMinTimeBin,fMaxTimeBin) ;
  fCuts->Add(new TObjString(onePar));
  snprintf(onePar,buffersize, "Switchs: Remove Bad Channels? %d; Use filtered input? %d;  Correct Clusters? %d, Mass per channel same SM clusters? %d ",
           fRecoUtils->IsBadChannelsRemovalSwitchedOn(),fFilteredInput,fCorrectClusters, fSameSM) ;
  fCuts->Add(new TObjString(onePar));
  snprintf(onePar,buffersize, "EMCAL Geometry name: < %s >, Load Matrices? %d",fEMCALGeoName.Data(),fLoadMatrices) ;
  fCuts->Add(new TObjString(onePar));
  
  fCuts ->SetOwner(kTRUE);
  
  // Post Data
  PostData(2, fCuts);
  
}

//__________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::UserCreateOutputObjects()
{
  //Create output container, init geometry 
  
  fEMCALGeo =  AliEMCALGeometry::GetInstance(fEMCALGeoName) ;	
  Int_t nSM = (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules();
  
  fOutputContainer = new TList();
  const Int_t buffersize = 255;
  char hname[buffersize], htitl[buffersize];
  
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
  
  Int_t nchannels = nSM*AliEMCALGeoParams::fgkEMCALRows*AliEMCALGeoParams::fgkEMCALCols;
  for(Int_t ibc = 0; ibc < 4; ibc++){
    fHTpi0[ibc] = new TH2F(Form("hTime_BC%d",ibc),Form("Time of cell clusters under pi0 peak, bunch crossing %d",ibc),
                           nchannels,0,nchannels, fNTimeBins,fMinTimeBin,fMaxTimeBin);
    fOutputContainer->Add(fHTpi0[ibc]);       
    fHTpi0[ibc]->SetYTitle("time (ns)");
    fHTpi0[ibc]->SetXTitle("abs. Id. ");
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
  
  fHIncidentAngle = new TH2F("hinang","#gamma incident angle in SM",100,0.,20.,100,0,10);
  fHIncidentAngle->SetXTitle("#alpha_{#gamma SM center}");
  fHIncidentAngle->SetYTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fHIncidentAngle);
  
  fHIncidentAngleDifferentSM = new TH2F("hinangDifferentSM","#gamma incident angle in SM, different SM pair",100,0,20.,100,0,10);
  fHIncidentAngleDifferentSM->SetXTitle("#alpha_{#gamma - SM center}");
  fHIncidentAngleDifferentSM->SetYTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fHIncidentAngleDifferentSM);
  
  fHAsymmetry = new TH2F("hasym","2-cluster opening angle",100,0.,1.,100,0,10);
  fHAsymmetry->SetXTitle("a");
  fHAsymmetry->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
  fOutputContainer->Add(fHAsymmetry);
  
  fHAsymmetryDifferentSM = new TH2F("hasymDifferentSM","2-cluster opening angle, different SM",100,0,1.,100,0,10);
  fHAsymmetryDifferentSM->SetXTitle("a");
  fHAsymmetryDifferentSM->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
  fOutputContainer->Add(fHAsymmetryDifferentSM);
  
  
  //TString pairname[] = {"A side (0-2)", "C side (1-3)","Row 0 (0-1)", "Row 1 (2-3)"};
  
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
    
    
    if(iSM < nSM/2){
      snprintf(hname,buffersize, "hmgg_PairSameSectorSM%d",iSM);
      snprintf(htitl,buffersize, "Two-gamma inv. mass for SM pair Sector: %d",iSM);
      fHmggPairSameSectorSM[iSM] = new TH2F(hname,htitl,fNbins,fMinBin,fMaxBin,100,0,10);
      fHmggPairSameSectorSM[iSM]->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
      fHmggPairSameSectorSM[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
      fOutputContainer->Add(fHmggPairSameSectorSM[iSM]);
      
      snprintf(hname,buffersize, "hmgg_PairSameSectorSM%d_MaskFrame",iSM);
      snprintf(htitl,buffersize, "Two-gamma inv. mass for SM pair Sector: %d",iSM);
      fHmggPairSameSectorSMMaskFrame[iSM] = new TH2F(hname,htitl,fNbins,fMinBin,fMaxBin,100,0,10);
      fHmggPairSameSectorSMMaskFrame[iSM]->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
      fHmggPairSameSectorSMMaskFrame[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
      fOutputContainer->Add(fHmggPairSameSectorSMMaskFrame[iSM]);
      
      fhClusterPairDiffTimeSameSector[iSM] = new TH2F(Form("hClusterPairDiffTimeSameSector%d",iSM),
                                                      Form("cluster pair time difference vs E, Sector %d",iSM),
                                                      100,0,10, 200,-100,100);
      fhClusterPairDiffTimeSameSector[iSM]->SetXTitle("E_{pair} (GeV)");
      fhClusterPairDiffTimeSameSector[iSM]->SetYTitle("#Delta t (ns)");
      fOutputContainer->Add(fhClusterPairDiffTimeSameSector[iSM]);
      
      
    }
    
    if(iSM < nSM-2){
      snprintf(hname,buffersize, "hmgg_PairSameSideSM%d",iSM);
      snprintf(htitl,buffersize, "Two-gamma inv. mass for SM pair Sector: %d",iSM);
      fHmggPairSameSideSM[iSM] = new TH2F(hname,htitl,fNbins,fMinBin,fMaxBin,100,0,10);
      fHmggPairSameSideSM[iSM]->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
      fHmggPairSameSideSM[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
      fOutputContainer->Add(fHmggPairSameSideSM[iSM]);
      
      snprintf(hname,buffersize, "hmgg_PairSameSideSM%d_MaskFrame",iSM);
      snprintf(htitl,buffersize, "Two-gamma inv. mass for SM pair Sector: %d",iSM);
      fHmggPairSameSideSMMaskFrame[iSM] = new TH2F(hname,htitl,fNbins,fMinBin,fMaxBin,100,0,10);
      fHmggPairSameSideSMMaskFrame[iSM]->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
      fHmggPairSameSideSMMaskFrame[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
      fOutputContainer->Add(fHmggPairSameSideSMMaskFrame[iSM]);   
      
      fhClusterPairDiffTimeSameSide[iSM] = new TH2F(Form("hClusterPairDiffTimeSameSide%d",iSM),
                                                    Form("cluster pair time difference vs E,  Side %d",iSM),
                                                    100,0,10, 200,-100,100);
      fhClusterPairDiffTimeSameSide[iSM]->SetXTitle("E_{pair} (GeV)");
      fhClusterPairDiffTimeSameSide[iSM]->SetYTitle("#Delta t (ns)");
      fOutputContainer->Add(fhClusterPairDiffTimeSameSide[iSM]);
      
    }
    
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
    
    snprintf(hname, buffersize, "hinang_SM%d",iSM);
    snprintf(htitl, buffersize, "Incident angle for super mod %d",iSM);
    fHIncidentAngleSM[iSM] = new TH2F(hname,htitl,100,0.,20.,100,0,10);
    fHIncidentAngleSM[iSM]->SetXTitle("#alpha_{#gamma - SM center} (deg)");
    fHIncidentAngleSM[iSM]->SetYTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fHIncidentAngleSM[iSM]);
    
    snprintf(hname,buffersize, "hinang_PairSM%d",iSM);
    snprintf(htitl,buffersize, "Incident angle for SM pair: %d",iSM);
    fHIncidentAnglePairSM[iSM] = new TH2F(hname,htitl,100,0.,20.,100,0,10);
    fHIncidentAnglePairSM[iSM]->SetXTitle("#alpha_{#gamma - SM center} (deg)");
    fHIncidentAnglePairSM[iSM]->SetYTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fHIncidentAnglePairSM[iSM]);   
    
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
    
    fhClusterTimeSM[iSM] = new TH2F(Form("hClusterTime_SM%d",iSM),"cluster time vs E",100,0,10, 100,0,1000);
    fhClusterTimeSM[iSM]->SetXTitle("E (GeV)");
    fhClusterTimeSM[iSM]->SetYTitle("t (ns)");
    fOutputContainer->Add(fhClusterTimeSM[iSM]);
    
    fhClusterPairDiffTimeSameSM[iSM] = new TH2F(Form("hClusterPairDiffTimeSameSM%d",iSM),
                                                Form("cluster pair time difference vs E, SM %d",iSM),
                                                100,0,10, 200,-100,100);
    fhClusterPairDiffTimeSameSM[iSM]->SetXTitle("E (GeV)");
    fhClusterPairDiffTimeSameSM[iSM]->SetYTitle("#Delta t (ns)");
    fOutputContainer->Add(fhClusterPairDiffTimeSameSM[iSM]);
    
  }  
  
  fhClusterTime = new TH2F("hClusterTime","cluster time vs E",100,0,10, 100,0,1000);
  fhClusterTime->SetXTitle("E (GeV)");
  fhClusterTime->SetYTitle("t (ns)");
  fOutputContainer->Add(fhClusterTime);
  
  fhClusterPairDiffTime = new TH2F("hClusterPairDiffTime","cluster pair time difference vs E",100,0,10, 800,-400,400);
  fhClusterPairDiffTime->SetXTitle("E_{pair} (GeV)");
  fhClusterPairDiffTime->SetYTitle("#Delta t (ns)");
  fOutputContainer->Add(fhClusterPairDiffTime);
  
  
  fhNEvents        = new TH1I("hNEvents", "Number of analyzed events"   , 1 , 0 , 1  ) ;
  fOutputContainer->Add(fhNEvents);
  
  fOutputContainer->SetOwner(kTRUE);
  
  //  fCalibData = new AliEMCALCalibData();
  
  PostData(1,fOutputContainer);
  
}

//__________________________________________________
Bool_t AliAnalysisTaskEMCALPi0CalibSelection::MaskFrameCluster(const Int_t iSM,  const Int_t ieta) const {
  //Check if cell is in one of the regions where we have significant amount of material in front of EMCAL
  
  Int_t icol = ieta;
  if(iSM%2) icol+=48; // Impair SM, shift index [0-47] to [48-96]
  
  if (fNMaskCellColumns && fMaskCellColumns) {
    for (Int_t imask = 0; imask < fNMaskCellColumns; imask++) {
      if(icol==fMaskCellColumns[imask]) return kTRUE;
    }
  }
  
  return kFALSE;
  
}

//__________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::UserExec(Option_t* /* option */)
{
  //Analysis per event.
  
  if(fRecoUtils->GetParticleType()!=AliEMCALRecoUtils::kPhoton){
    printf("Wrong particle type for cluster position recalculation! = %d\n", fRecoUtils->GetParticleType());
    abort();
  }
  
  if(!(((AliESDEvent*)InputEvent())->GetFiredTriggerClasses()).Contains(fTriggerName)) {
    //printf("Reject Event %d, FiredClass %s\n",(Int_t)Entry(),(((AliESDEvent*)InputEvent())->GetFiredTriggerClasses()).Data());
    return;
  }
  
  fhNEvents->Fill(0); //Event analyzed
  
  //Get the input event
  AliVEvent* event = 0;
  if(fFilteredInput) event = AODEvent();
  else               event = InputEvent();
  
  if(!event) {
    printf("Input event not available!\n");
    return;
  }
  
  if(DebugLevel() > 1) 
    printf("AliAnalysisTaskEMCALPi0CalibSelection <<< %s: Event %d >>>\n",event->GetName(), (Int_t)Entry());
  
  
  //Get the primary vertex
  Double_t v[3];
  event->GetPrimaryVertex()->GetXYZ(v) ;
  
  if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Vertex: (%.3f,%.3f,%.3f)\n",v[0],v[1],v[2]);
  
  //Int_t runNum = aod->GetRunNumber();
  //if(DebugLevel() > 1) printf("Run number: %d\n",runNum);
  
  Int_t nSM = (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules();
  //Get the matrix with geometry information
  if(fhNEvents->GetEntries()==1){
    if(fLoadMatrices){
      printf("AliAnalysisTaskEMCALPi0CalibSelection::UserExec() - Load user defined geometry matrices\n");
      for(Int_t mod=0; mod < nSM ; mod++){
        if(fMatrix[mod]){
          if(DebugLevel() > 1) 
            fMatrix[mod]->Print();
          fEMCALGeo->SetMisalMatrix(fMatrix[mod],mod) ;   
        }
      }//SM loop
    }//Load matrices
    else if(!gGeoManager){
      printf("AliAnalysisTaskEMCALPi0CalibSelection::UserExec() - Get geo matrices from data\n");
      //Still not implemented in AOD, just a workaround to be able to work at least with ESDs	
      if(!strcmp(event->GetName(),"AliAODEvent")) {
        if(DebugLevel() > 1) 
          printf("AliAnalysisTaskEMCALPi0CalibSelection Use ideal geometry, values geometry matrix not kept in AODs.\n");
      }//AOD
      else {	
        if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Load Misaligned matrices. \n");
        AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event) ;
        if(!esd) {
          printf("AliAnalysisTaskEMCALPi0CalibSelection::UserExec() - This event does not contain ESDs?");
          return;
        }
        for(Int_t mod=0; mod < nSM; mod++){
          if(DebugLevel() > 1) 
            esd->GetEMCALMatrix(mod)->Print();
          if(esd->GetEMCALMatrix(mod)) fEMCALGeo->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
        } 
      }//ESD
    }//Load matrices from Data 
  }//first event
  
  if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Will use fLogWeight %.3f .\n",fLogWeight);
  Int_t absId1   = -1;
  Int_t iSupMod1 = -1;
  Int_t iphi1    = -1;
  Int_t ieta1    = -1;
  Int_t absId2   = -1;
  Int_t iSupMod2 = -1;
  Int_t iphi2    = -1;
  Int_t ieta2    = -1;
  Bool_t shared  = kFALSE;
  
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector p12;
  
  //Get the list of clusters
  TRefArray * caloClustersArr  = new TRefArray();
  event->GetEMCALClusters(caloClustersArr);
  const Int_t kNumberOfEMCALClusters   = caloClustersArr->GetEntries() ;
  if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection - N CaloClusters: %d \n", kNumberOfEMCALClusters);
  
  // Get EMCAL cells
  AliVCaloCells *emCells = event->GetEMCALCells();
  
  // loop over EMCAL clusters
  //----------------------------------------------------------
  // First recalibrate and recalculate energy and position
  Float_t pos[]={0,0,0};
  if(fCorrectClusters){
    for(Int_t iClu=0; iClu<kNumberOfEMCALClusters; iClu++) {
      AliVCluster *c1 = (AliVCluster *) caloClustersArr->At(iClu);
      
      Float_t e1i = c1->E();   // cluster energy before correction   
      if      (e1i < fEmin) continue;
      else if (e1i > fEmax) continue;
      else if (c1->GetNCells() < fMinNCells) continue; 
      else if (c1->GetM02() < fL0min || c1->GetM02() > fL0max) continue;
        
      if(fRecoUtils->ClusterContainsBadChannel(fEMCALGeo, c1->GetCellsAbsId(), c1->GetNCells())) continue;	
      
      if(DebugLevel() > 2)
      { 
        printf("Std  : i %d, E %f, dispersion %f, m02 %f, m20 %f\n",c1->GetID(),c1->E(),c1->GetDispersion(),c1->GetM02(),c1->GetM20());
        c1->GetPosition(pos);
        printf("Std  : i %d, x %f, y %f, z %f\n",c1->GetID(), pos[0], pos[1], pos[2]);
      }
      
      //Correct cluster energy and position if requested, and not corrected previously, by default Off
      if(fRecoUtils->IsRecalibrationOn())	{
        fRecoUtils->RecalibrateClusterEnergy(fEMCALGeo, c1, emCells);
        fRecoUtils->RecalculateClusterShowerShapeParameters(fEMCALGeo, emCells,c1);
        fRecoUtils->RecalculateClusterPID(c1);
      }
      if(DebugLevel() > 2) 
        printf("Energy: after recalibration %f; \n",c1->E());
      
      // Recalculate cluster position
      fRecoUtils->RecalculateClusterPosition(fEMCALGeo, emCells,c1);
      
      // Correct Non-Linearity
      c1->SetE(fRecoUtils->CorrectClusterEnergyLinearity(c1));
      if(DebugLevel() > 2) 
        printf("\t after linearity correction %f\n",c1->E());
      
      //In case of MC analysis, to match resolution/calibration in real data
      c1->SetE(fRecoUtils->SmearClusterEnergy(c1));
      if(DebugLevel() > 2) 
        printf("\t after smearing %f\n",c1->E());      
      
      if(DebugLevel() > 2)
      { 
        printf("Cor  : i %d, E %f, dispersion %f, m02 %f, m20 %f\n",c1->GetID(),c1->E(),c1->GetDispersion(),c1->GetM02(),c1->GetM20());
        c1->GetPosition(pos);
        printf("Cor  : i %d, x %f, y %f, z %f\n",c1->GetID(), pos[0], pos[1], pos[2]);
      }    
    }    
  }
  
  //----------------------------------------------------------
  //Now the invariant mass analysis with the corrected clusters
  Int_t bc = event->GetBunchCrossNumber();
  for(Int_t iClu=0; iClu<kNumberOfEMCALClusters-1; iClu++) {
    
    AliVCluster *c1 = (AliVCluster *) caloClustersArr->At(iClu);
    if(fRecoUtils->ClusterContainsBadChannel(fEMCALGeo, c1->GetCellsAbsId(), c1->GetNCells())) continue;	
    
    Float_t e1i = c1->E();   // cluster energy before correction   
    if      (e1i < fEmin) continue;
    else if (e1i > fEmax) continue;
    else if (!fRecoUtils->IsGoodCluster(c1,fEMCALGeo,emCells,bc));
    else if (c1->GetNCells() < fMinNCells) continue; 
    else if (c1->GetM02() < fL0min || c1->GetM02() > fL0max) continue;

    if(DebugLevel() > 2)
    { 
      printf("IMA  : i %d, E %f, dispersion %f, m02 %f, m20 %f\n",c1->GetID(),e1i,c1->GetDispersion(),c1->GetM02(),c1->GetM20());
      c1->GetPosition(pos);
      printf("IMA  : i %d, x %f, y %f, z %f\n",c1->GetID(), pos[0], pos[1], pos[2]);
    }
    
    fRecoUtils->GetMaxEnergyCell(fEMCALGeo, emCells,c1,absId1,iSupMod1,ieta1,iphi1,shared);
    c1->GetMomentum(p1,v);
    
    //Check if cluster is in fidutial region, not too close to borders
    Bool_t in1 = fRecoUtils->CheckCellFiducialRegion(fEMCALGeo, c1, emCells);
    // Clusters not facing frame structures
    Bool_t mask1 = MaskFrameCluster(iSupMod1, ieta1);
    //if(mask1) printf("Reject eta %d SM %d\n",ieta1, iSupMod1);
    
    Double_t time1 = c1->GetTOF()*1.e9;
    fhClusterTime            ->Fill(c1->E(),time1);
    fhClusterTimeSM[iSupMod1]->Fill(c1->E(),time1);
    
    // Combine cluster with other clusters and get the invariant mass
    for (Int_t jClu=iClu+1; jClu<kNumberOfEMCALClusters; jClu++) {
      AliAODCaloCluster *c2 = (AliAODCaloCluster *) caloClustersArr->At(jClu);
      

      Float_t e2i = c2->E();
      if      (e2i < fEmin) continue;
      else if (e2i > fEmax) continue;
      else if (!fRecoUtils->IsGoodCluster(c2,fEMCALGeo,emCells,bc))continue;
      else if (c2->GetNCells() < fMinNCells) continue; 
      else if (c2->GetM02() < fL0min || c2->GetM02() > fL0max) continue;

      
      fRecoUtils->GetMaxEnergyCell(fEMCALGeo, emCells,c2,absId2,iSupMod2,ieta2,iphi2,shared);
      c2->GetMomentum(p2,v);
      
      p12 = p1+p2;
      Float_t invmass = p12.M()*1000; 
      //printf("*** mass %f\n",invmass);
      
      //Asimetry cut      
      Float_t asym = TMath::Abs(p1.E()-p2.E())/(p1.E()+p2.E());
      //printf("asymmetry %f\n",asym);
      if(asym > fAsyCut) continue;
      
      //Time cut
      Double_t time2 = c2->GetTOF()*1.e9;
      fhClusterPairDiffTime->Fill(p12.E(),time1-time2);
      if(TMath::Abs(time1-time2) > fDTimeCut) continue;
      
      if(invmass < fMaxBin && invmass > fMinBin ){
        
        //Check if cluster is in fidutial region, not too close to borders
        Bool_t in2 = fRecoUtils->CheckCellFiducialRegion(fEMCALGeo, c2, emCells);
        
        // Clusters not facing frame structures
        Bool_t mask2 = MaskFrameCluster(iSupMod2, ieta2);         
        //if(mask2) printf("Reject eta %d SM %d\n",ieta2, iSupMod2);
        
        if(in1 && in2){
          
          fHmgg->Fill(invmass,p12.Pt()); 
          
          if(iSupMod1==iSupMod2) {
            fHmggSM[iSupMod1]->Fill(invmass,p12.Pt()); 
            fhClusterPairDiffTimeSameSM[iSupMod1]->Fill(p12.E(),time1-time2);
          }
          else                   
            fHmggDifferentSM ->Fill(invmass,p12.Pt());
          
          // Same sector
          Int_t j=0;
          for(Int_t i = 0; i < nSM/2; i++){
            j=2*i;
            if((iSupMod1==j && iSupMod2==j+1) || (iSupMod1==j+1 && iSupMod2==j)) {
              fHmggPairSameSectorSM[i]->Fill(invmass,p12.Pt());
              fhClusterPairDiffTimeSameSector[i]->Fill(p12.E(),time1-time2);
            } 
          }
          
          // Same side
          for(Int_t i = 0; i < nSM-2; i++){
            if((iSupMod1==i && iSupMod2==i+2) || (iSupMod1==i+2 && iSupMod2==i)) {
              fHmggPairSameSideSM[i]->Fill(invmass,p12.Pt()); 
              fhClusterPairDiffTimeSameSide[i]->Fill(p12.E(),time1-time2);
            }
          }
          
          
          if(!mask1 && !mask2){
            
            fHmggMaskFrame->Fill(invmass,p12.Pt()); 
            
            if(iSupMod1==iSupMod2) fHmggSMMaskFrame[iSupMod1]->Fill(invmass,p12.Pt()); 
            else                   fHmggDifferentSMMaskFrame ->Fill(invmass,p12.Pt());
            
            // Same sector
            j=0;
            for(Int_t i = 0; i < nSM/2; i++){
              j=2*i;
              if((iSupMod1==j && iSupMod2==j+1) || (iSupMod1==j+1 && iSupMod2==j)) fHmggPairSameSectorSMMaskFrame[i]->Fill(invmass,p12.Pt()); 
            }
            
            // Same side
            for(Int_t i = 0; i < nSM-2; i++){
              if((iSupMod1==i && iSupMod2==i+2) || (iSupMod1==i+2 && iSupMod2==i)) fHmggPairSameSideSMMaskFrame[i]->Fill(invmass,p12.Pt()); 
            }
            
          }// Pair not facing frame
          
          
          if(invmass > fInvMassCutMin && invmass < fInvMassCutMax){//restrict to clusters really close to pi0 peak
            
            
            // Check time of cells in both clusters, and fill time histogram
            for(Int_t icell = 0; icell < c1->GetNCells(); icell++){
              Int_t absID = c1->GetCellAbsId(icell);   
              fHTpi0[bc%4]->Fill(absID, emCells->GetCellTime(absID)*1.e9);  
            }
            
            for(Int_t icell = 0; icell < c2->GetNCells(); icell++){
              Int_t absID = c2->GetCellAbsId(icell);   
              fHTpi0[bc%4]->Fill(absID, emCells->GetCellTime(absID)*1.e9);  
            }
            
             //Opening angle of 2 photons
            Float_t opangle = p1.Angle(p2.Vect())*TMath::RadToDeg();
            //printf("*******>>>>>>>> In PEAK pt %f, angle %f \n",p12.Pt(),opangle);
            
            //Incident angle of each photon
            Float_t inangle1 =0., inangle2=0.;
            if(gGeoManager){
              Float_t posSM1cen[3]={0.,0.,0.};
              Float_t depth = fRecoUtils->GetDepth(p1.Energy(),fRecoUtils->GetParticleType(),iSupMod1);
              fEMCALGeo->RecalculateTowerPosition(11.5, 23.5, iSupMod1, depth, fRecoUtils->GetMisalTransShiftArray(),fRecoUtils->GetMisalRotShiftArray(),posSM1cen); 
              Float_t posSM2cen[3]={0.,0.,0.}; 
              depth = fRecoUtils->GetDepth(p2.Energy(),fRecoUtils->GetParticleType(),iSupMod2);
              fEMCALGeo->RecalculateTowerPosition(11.5, 23.5, iSupMod2, depth, fRecoUtils->GetMisalTransShiftArray(),fRecoUtils->GetMisalRotShiftArray(),posSM2cen); 
              //printf("SM1 %d pos (%2.3f,%2.3f,%2.3f) \n",iSupMod1,posSM1cen[0],posSM1cen[1],posSM1cen[2]);
              //printf("SM2 %d pos (%2.3f,%2.3f,%2.3f) \n",iSupMod2,posSM2cen[0],posSM2cen[1],posSM2cen[2]);
              
              TVector3 vecSM1cen(posSM1cen[0]-v[0],posSM1cen[1]-v[1],posSM1cen[2]-v[2]); 
              TVector3 vecSM2cen(posSM2cen[0]-v[0],posSM2cen[1]-v[1],posSM2cen[2]-v[2]); 
              inangle1 = p1.Angle(vecSM1cen)*TMath::RadToDeg();
              inangle2 = p2.Angle(vecSM2cen)*TMath::RadToDeg();
              //printf("Incident angle: cluster 1 %2.3f; cluster 2 %2.3f\n",inangle1,inangle2);
            }
            
            fHOpeningAngle ->Fill(opangle,p12.Pt()); 
            fHIncidentAngle->Fill(inangle1,p1.Pt()); 
            fHIncidentAngle->Fill(inangle2,p2.Pt()); 
            fHAsymmetry    ->Fill(asym,p12.Pt()); 
            
            if(iSupMod1==iSupMod2) {
              fHOpeningAngleSM[iSupMod1] ->Fill(opangle,p12.Pt());
              fHIncidentAngleSM[iSupMod1]->Fill(inangle1,p1.Pt());
              fHIncidentAngleSM[iSupMod1]->Fill(inangle2,p2.Pt());
              fHAsymmetrySM[iSupMod1]    ->Fill(asym,p12.Pt());
            }
            else{      
              fHOpeningAngleDifferentSM  ->Fill(opangle,p12.Pt());
              fHIncidentAngleDifferentSM ->Fill(inangle1,p1.Pt());
              fHIncidentAngleDifferentSM ->Fill(inangle2,p2.Pt());
              fHAsymmetryDifferentSM     ->Fill(asym,p12.Pt());
            }
            
            if((iSupMod1==0 && iSupMod2==2) || (iSupMod1==2 && iSupMod2==0)) {
              fHOpeningAnglePairSM[0] ->Fill(opangle,p12.Pt()); 
              fHIncidentAnglePairSM[0]->Fill(inangle1,p1.Pt());
              fHIncidentAnglePairSM[0]->Fill(inangle2,p2.Pt());
              fHAsymmetryPairSM[0]    ->Fill(asym,p12.Pt());
              
            } 
            if((iSupMod1==1 && iSupMod2==3) || (iSupMod1==3 && iSupMod2==1)) {
              fHOpeningAnglePairSM[1] ->Fill(opangle,p12.Pt()); 
              fHIncidentAnglePairSM[1]->Fill(inangle1,p1.Pt());
              fHIncidentAnglePairSM[1]->Fill(inangle2,p2.Pt());
              fHAsymmetryPairSM[1]    ->Fill(asym,p12.Pt());
              
            }
            
            if((iSupMod1==0 && iSupMod2==1) || (iSupMod1==1 && iSupMod2==0)) {
              fHOpeningAnglePairSM[2] ->Fill(opangle,p12.Pt()); 
              fHIncidentAnglePairSM[2]->Fill(inangle1,p1.Pt());
              fHIncidentAnglePairSM[2]->Fill(inangle2,p2.Pt());
              fHAsymmetryPairSM[2]    ->Fill(asym,p12.Pt());
              
              
            }
            if((iSupMod1==2 && iSupMod2==3) || (iSupMod1==3 && iSupMod2==2)) {
              fHOpeningAnglePairSM[3] ->Fill(opangle,p12.Pt()); 
              fHIncidentAnglePairSM[3]->Fill(inangle1,p1.Pt());
              fHIncidentAnglePairSM[3]->Fill(inangle2,p2.Pt());
              fHAsymmetryPairSM[3]    ->Fill(asym,p12.Pt());
            }
            
          }// pair in 100 < mass < 160
          
        }//in acceptance cuts
        
        //In case of filling only channels with second cluster in same SM
        if(fSameSM && iSupMod1!=iSupMod2) continue;
        
        if (fGroupNCells == 0){
          fHmpi0[iSupMod1][ieta1][iphi1]->Fill(invmass);
          fHmpi0[iSupMod2][ieta2][iphi2]->Fill(invmass);
          
          if(invmass > 100. && invmass < 160.){//restrict to clusters really close to pi0 peak
            fhTowerDecayPhotonHit      [iSupMod1]->Fill(ieta1,iphi1);
            fhTowerDecayPhotonEnergy   [iSupMod1]->Fill(ieta1,iphi1,p1.E());
            fhTowerDecayPhotonAsymmetry[iSupMod1]->Fill(ieta1,iphi1,asym);
            
            fhTowerDecayPhotonHit      [iSupMod2]->Fill(ieta2,iphi2);
            fhTowerDecayPhotonEnergy   [iSupMod2]->Fill(ieta2,iphi2,p2.E());
            fhTowerDecayPhotonAsymmetry[iSupMod2]->Fill(ieta2,iphi2,asym);
            
            if(!mask1)fhTowerDecayPhotonHitMaskFrame[iSupMod1]->Fill(ieta1,iphi1);
            if(!mask2)fhTowerDecayPhotonHitMaskFrame[iSupMod2]->Fill(ieta2,iphi2);
            
          }// pair in mass of pi0
        }	
        else  {
          //printf("Regroup N %d, eta1 %d, phi1 %d, eta2 %d, phi2 %d \n",fGroupNCells, ieta1, iphi1, ieta2, iphi2);
          for (Int_t i = -fGroupNCells; i < fGroupNCells+1; i++) {
            for (Int_t j = -fGroupNCells; j < fGroupNCells+1; j++) {
              //printf("\t i %d, j %d\n",i,j);
              
              Int_t absId11 = fEMCALGeo->GetAbsCellIdFromCellIndexes(iSupMod1, iphi1+j, ieta1+i); 
              Int_t absId22 = fEMCALGeo->GetAbsCellIdFromCellIndexes(iSupMod2, iphi2+j, ieta2+i); 
              Bool_t ok1 = kFALSE;
              Bool_t ok2 = kFALSE;
              for(Int_t icell = 0; icell < c1->GetNCells(); icell++){
                if(c1->GetCellsAbsId()[icell] == absId11) ok1=kTRUE;
              }
              for(Int_t icell = 0; icell < c2->GetNCells(); icell++){
                if(c2->GetCellsAbsId()[icell] == absId22) ok2=kTRUE;
              }
              
              if(ok1 && (ieta1+i >= 0) && (iphi1+j >= 0) && (ieta1+i < 48) && (iphi1+j < 24)){
                //printf("\t \t eta1+i %d, phi1+j %d\n", ieta1+i, iphi1+j);
                fHmpi0[iSupMod1][ieta1+i][iphi1+j]->Fill(invmass);
              }
              if(ok2 && (ieta2+i >= 0) && (iphi2+j >= 0) && (ieta2+i < 48) && (iphi2+j < 24)){
                //printf("\t \t eta2+i %d, phi2+j %d\n", ieta2+i, iphi2+j);
                fHmpi0[iSupMod2][ieta2+i][iphi2+j]->Fill(invmass);
              }
            }// j loop
          }//i loop
        }//group cells
        
        if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Mass in (SM%d,%d,%d) and  (SM%d,%d,%d): %.3f GeV  E1_i=%f E1_ii=%f  E2_i=%f E2_ii=%f\n",
                                    iSupMod1,iphi1,ieta1,iSupMod2,iphi2,ieta2,p12.M(),e1i,c1->E(),e2i,c2->E());
      }
      
    }
    
  } // end of loop over EMCAL clusters
  
  delete caloClustersArr;
  
  PostData(1,fOutputContainer);
  
}

//_____________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::PrintInfo(){
  
  //Print settings
  printf("Cluster cuts: %2.2f < E < %2.2f GeV; number of cells > %d; Assymetry < %1.2f, pair time diff < %2.2f\n", 
         fEmin,fEmax, fMinNCells, fAsyCut, fDTimeCut) ;
  printf("Group %d cells\n", fGroupNCells) ;
  printf("Cluster maximal cell away from border at least %d cells\n", fRecoUtils->GetNumberOfCellsFromEMCALBorder()) ;
  printf("Histograms: bins %d; energy range: %2.2f < E < %2.2f GeV\n",fNbins,fMinBin,fMaxBin) ;
  printf("Switchs:\n \t Remove Bad Channels? %d; Use filtered input? %d;  Correct Clusters? %d, \n \t Mass per channel same SM clusters? %d\n",
         fRecoUtils->IsBadChannelsRemovalSwitchedOn(),fFilteredInput,fCorrectClusters, fSameSM) ;
  printf("EMCAL Geometry name: < %s >, Load Matrices %d\n",fEMCALGeoName.Data(), fLoadMatrices) ;
  if(fLoadMatrices) {for(Int_t ism = 0; ism < AliEMCALGeoParams::fgkEMCALModules; ism++) fMatrix[ism]->Print() ; }
  
  
}

