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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  EMCAL tender, apply corrections to EMCAL clusters                        //
//  and do track matching.                                                   //
//  Author: Deepa Thomas (Utrecht University)                                // 
//  Later mods/rewrite: Jiri Kral (University of Jyvaskyla)                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TInterpreter.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliAnalysisManager.h>
#include <AliTender.h>
#include "AliOADBContainer.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliMagF.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALTenderSupply.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALClusterizer.h"
#include "AliEMCALRecParam.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALAfterBurnerUF.h"
#include "AliEMCALClusterizerNxN.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALClusterizerv2.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRecParam.h"

ClassImp(AliEMCALTenderSupply)

AliEMCALTenderSupply::AliEMCALTenderSupply() :
AliTenderSupply()
,fEMCALGeo(0x0)
,fEMCALGeoName("EMCAL_COMPLETEV1")
,fEMCALRecoUtils(0)
,fConfigName("")
,fDebugLevel(0)
,fNonLinearFunc(-1) 
,fNonLinearThreshold(-1)
,fReCalibCluster(kFALSE)
,fUpdateCell(kFALSE)  
,fCalibrateEnergy(kFALSE)
,fCalibrateTime(kFALSE)
,fDoNonLinearity(kFALSE)
,fBadCellRemove(kFALSE)
,fRejectExoticCells(kFALSE)
,fRejectExoticClusters(kFALSE)
,fClusterBadChannelCheck(kFALSE)
,fRecalClusPos(kFALSE)
,fFiducial(kFALSE) 
,fNCellsFromEMCALBorder(-1)
,fRecalDistToBadChannels(kFALSE)
,fRecalShowerShape(kFALSE)
,fInputTree(0)
,fInputFile(0)
,fFilepass(0) 
,fMass(-1)
,fStep(-1)
,fCutEtaPhiSum(kTRUE)
,fCutEtaPhiSeparate(kFALSE)
,fRcut(-1)
,fEtacut(-1)
,fPhicut(-1)
,fBasePath("")
,fReClusterize(kFALSE)
,fClusterizer(0)
,fGeomMatrixSet(kFALSE)
,fLoadGeomMatrices(kFALSE)
,fRecParam(0x0)
,fDoTrackMatch(kFALSE)
,fDoUpdateOnly(kFALSE)
,fUnfolder(0)
,fDigitsArr(0)
,fClusterArr(0)
,fMisalignSurvey(kdefault)  
,fExoticCellFraction(-1)
,fExoticCellDiffTime(-1)
,fExoticCellMinAmplitude(-1)
,fRecoParamsOCDBLoaded(kFALSE)
{
  // Default constructor.
  for(Int_t i = 0; i < 10; i++) fEMCALMatrix[i] = 0 ;
  fEMCALRecoUtils  = new AliEMCALRecoUtils;
}

//_____________________________________________________
AliEMCALTenderSupply::AliEMCALTenderSupply(const char *name, const AliTender *tender) :
AliTenderSupply(name,tender)
,fEMCALGeo(0x0)
,fEMCALGeoName("EMCAL_COMPLETEV1")
,fEMCALRecoUtils(0)
,fConfigName("")
,fDebugLevel(0)
,fNonLinearFunc(-1) 
,fNonLinearThreshold(-1)        
,fReCalibCluster(kFALSE)  
,fUpdateCell(kFALSE)  
,fCalibrateEnergy(kFALSE)
,fCalibrateTime(kFALSE)
,fDoNonLinearity(kFALSE)
,fBadCellRemove(kFALSE)
,fRejectExoticCells(kFALSE)
,fRejectExoticClusters(kFALSE)
,fClusterBadChannelCheck(kFALSE)
,fRecalClusPos(kFALSE)
,fFiducial(kFALSE) 
,fNCellsFromEMCALBorder(-1)  
,fRecalDistToBadChannels(kFALSE)  
,fRecalShowerShape(kFALSE)
,fInputTree(0)  
,fInputFile(0)
,fFilepass(0) 
,fMass(-1)
,fStep(-1)
,fCutEtaPhiSum(kTRUE)
,fCutEtaPhiSeparate(kFALSE)
,fRcut(-1)  
,fEtacut(-1)  
,fPhicut(-1)  
,fBasePath("")
,fReClusterize(kFALSE)
,fClusterizer(0)
,fGeomMatrixSet(kFALSE)
,fLoadGeomMatrices(kFALSE)
,fRecParam(0x0)
,fDoTrackMatch(kFALSE)
,fDoUpdateOnly(kFALSE)
,fUnfolder(0)
,fDigitsArr(0)
,fClusterArr(0)
,fMisalignSurvey(kdefault)  
,fExoticCellFraction(-1)
,fExoticCellDiffTime(-1)
,fExoticCellMinAmplitude(-1)
,fRecoParamsOCDBLoaded(kFALSE)
{
  // Named constructor
  
  for(Int_t i = 0; i < 10; i++) fEMCALMatrix[i] = 0 ;
  fEMCALRecoUtils  = new AliEMCALRecoUtils;
}

//_____________________________________________________
AliEMCALTenderSupply::~AliEMCALTenderSupply()
{
  //Destructor
  
  delete fEMCALRecoUtils;
  delete fRecParam;
  delete fUnfolder;
  if (!fClusterizer) {
    fDigitsArr->Clear("C");
    delete fDigitsArr; 
  } else {
    delete fClusterizer;
    fDigitsArr = 0;
  }
}

//_____________________________________________________
void AliEMCALTenderSupply::Init()
{
  // Initialise EMCAL tender.

  if (fDebugLevel>0) 
    AliWarning("Init EMCAL Tender supply"); 
  
  if (fConfigName.Length()>0 && gROOT->LoadMacro(fConfigName) >=0) {
    AliDebug(1, Form("Loading settings from macro %s", fConfigName.Data()));
    AliEMCALTenderSupply *tender = (AliEMCALTenderSupply*)gInterpreter->ProcessLine("ConfigEMCALTenderSupply()");
    fDebugLevel             = tender->fDebugLevel;
    fEMCALGeoName           = tender->fEMCALGeoName; 
    fEMCALRecoUtils         = tender->fEMCALRecoUtils; 
    fConfigName             = tender->fConfigName;
    fNonLinearFunc          = tender->fNonLinearFunc;
    fNonLinearThreshold     = tender->fNonLinearThreshold;
    fReCalibCluster         = tender->fReCalibCluster;
    fUpdateCell             = tender->fUpdateCell;
    fRecalClusPos           = tender->fRecalClusPos;
    fCalibrateEnergy        = tender->fCalibrateEnergy;
    fCalibrateTime          = tender->fCalibrateTime;
    fFiducial               = tender->fFiducial;
    fNCellsFromEMCALBorder  = tender->fNCellsFromEMCALBorder;
    fRecalDistToBadChannels = tender->fRecalDistToBadChannels;    
    fRecalShowerShape       = tender->fRecalShowerShape;
    fClusterBadChannelCheck = tender->fClusterBadChannelCheck;
    fBadCellRemove          = tender->fBadCellRemove;
    fRejectExoticCells      = tender->fRejectExoticCells;
    fRejectExoticClusters   = tender->fRejectExoticClusters;
    fMass                   = tender->fMass;
    fStep                   = tender->fStep;
    fCutEtaPhiSum           = tender->fCutEtaPhiSum;
    fCutEtaPhiSeparate      = tender->fCutEtaPhiSeparate;
    fRcut                   = tender->fRcut;
    fEtacut                 = tender->fEtacut;
    fPhicut                 = tender->fPhicut;
    fReClusterize           = tender->fReClusterize;
    fLoadGeomMatrices       = tender->fLoadGeomMatrices;
    fRecParam               = tender->fRecParam;
    fDoNonLinearity         = tender->fDoNonLinearity;
    fDoTrackMatch           = tender->fDoTrackMatch;
    fDoUpdateOnly           = tender->fDoUpdateOnly;
    fMisalignSurvey         = tender->fMisalignSurvey;
    fExoticCellFraction     = tender->fExoticCellFraction;
    fExoticCellDiffTime     = tender->fExoticCellDiffTime;
    fExoticCellMinAmplitude = tender->fExoticCellMinAmplitude;
    fRecoParamsOCDBLoaded   = tender->fRecoParamsOCDBLoaded;

    for(Int_t i = 0; i < 10; i++) 
      fEMCALMatrix[i] = tender->fEMCALMatrix[i] ;
  }
  
  if (fDebugLevel>0){
    AliInfo( "Emcal Tender settings: ======================================" ); 
    AliInfo( "------------ Switches --------------------------" ); 
    AliInfo( Form( "BadCellRemove : %d", fBadCellRemove )); 
    AliInfo( Form( "ExoticCellRemove : %d", fRejectExoticCells )); 
    AliInfo( Form( "CalibrateEnergy : %d", fCalibrateEnergy )); 
    AliInfo( Form( "CalibrateTime : %d", fCalibrateTime )); 
    AliInfo( Form( "UpdateCell : %d", fUpdateCell )); 
    AliInfo( Form( "DoUpdateOnly : %d", fDoUpdateOnly )); 
    AliInfo( Form( "Reclustering : %d", fReClusterize )); 
    AliInfo( Form( "ClusterBadChannelCheck : %d", fClusterBadChannelCheck )); 
    AliInfo( Form( "ClusterExoticChannelCheck : %d", fRejectExoticClusters )); 
    AliInfo( Form( "CellFiducialRegion : %d", fFiducial )); 
    AliInfo( Form( "ReCalibrateCluster : %d", fReCalibCluster )); 
    AliInfo( Form( "RecalculateClusPos : %d", fRecalClusPos )); 
    AliInfo( Form( "RecalShowerShape : %d", fRecalShowerShape )); 
    AliInfo( Form( "NonLinearityCorrection : %d", fDoNonLinearity )); 
    AliInfo( Form( "RecalDistBadChannel : %d", fRecalDistToBadChannels )); 
    AliInfo( Form( "TrackMatch : %d", fDoTrackMatch )); 
    AliInfo( "------------ Variables -------------------------" ); 
    AliInfo( Form( "DebugLevel : %d", fDebugLevel )); 
    AliInfo( Form( "BasePath : %s", fBasePath.Data() )); 
    AliInfo( Form( "ConfigFileName : %s", fConfigName.Data() )); 
    AliInfo( Form( "EMCALGeometryName : %s", fEMCALGeoName.Data() )); 
    AliInfo( Form( "NonLinearityFunction : %d", fNonLinearFunc )); 
    AliInfo( Form( "NonLinearityThreshold : %d", fNonLinearThreshold )); 
    AliInfo( Form( "MisalignmentMatrixSurvey : %d", fMisalignSurvey )); 
    AliInfo( Form( "NumberOfCellsFromEMCALBorder : %d", fNCellsFromEMCALBorder )); 
    AliInfo( Form( "RCut : %f", fRcut )); 
    AliInfo( Form( "Mass : %f", fMass )); 
    AliInfo( Form( "Step : %f", fStep )); 
    AliInfo( Form( "EtaCut : %f", fEtacut )); 
    AliInfo( Form( "PhiCut : %f", fPhicut )); 
    AliInfo( Form( "ExoticCellFraction : %f", fExoticCellFraction )); 
    AliInfo( Form( "ExoticCellDiffTime : %f", fExoticCellDiffTime )); 
    AliInfo( Form( "ExoticCellMinAmplitude : %f", fExoticCellMinAmplitude )); 
    AliInfo( "=============================================================" ); 
  }

  // Init geometry  
  fEMCALGeo = AliEMCALGeometry::GetInstance(fEMCALGeoName) ;

  // digits array
  fDigitsArr       = new TClonesArray("AliEMCALDigit",1000);

  // Initialising non-linearity parameters
  if( fNonLinearThreshold != -1 )
    fEMCALRecoUtils->SetNonLinearityThreshold(fNonLinearThreshold);
  if( fNonLinearFunc != -1 )
    fEMCALRecoUtils->SetNonLinearityFunction(fNonLinearFunc);

  // missalignment function
  fEMCALRecoUtils->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);

  // fiducial cut
  // do not do the eta0 fiducial cut
  if( fNCellsFromEMCALBorder != -1 )
    fEMCALRecoUtils->SetNumberOfCellsFromEMCALBorder(fNCellsFromEMCALBorder);
  fEMCALRecoUtils->SwitchOnNoFiducialBorderInEMCALEta0();
    
  // exotic cell rejection
  if( fExoticCellFraction != -1 )
    fEMCALRecoUtils->SetExoticCellFractionCut( fExoticCellFraction );
  if( fExoticCellDiffTime != -1 )
    fEMCALRecoUtils->SetExoticCellDiffTimeCut( fExoticCellDiffTime );
  if( fExoticCellMinAmplitude != -1 )
    fEMCALRecoUtils->SetExoticCellMinAmplitudeCut( fExoticCellMinAmplitude );

  // Setting track matching parameters ... mass, step size and residual cut
  if( fMass != -1 )
    fEMCALRecoUtils->SetMass(fMass);
  if( fStep != -1 )
    fEMCALRecoUtils->SetStep(fStep);
  
  // spatial cut based on separate eta/phi or common processing
  if(fCutEtaPhiSum){ 
    fEMCALRecoUtils->SwitchOnCutEtaPhiSum(); 
    if( fRcut != -1 )
      fEMCALRecoUtils->SetCutR(fRcut);
  } else if (fCutEtaPhiSeparate) {
    fEMCALRecoUtils->SwitchOnCutEtaPhiSeparate();
    if( fEtacut != -1 )
      fEMCALRecoUtils->SetCutEta(fEtacut);
    if( fPhicut != -1 )
      fEMCALRecoUtils->SetCutPhi(fPhicut);
  }
}

//_____________________________________________________
void AliEMCALTenderSupply::ProcessEvent()
{
  // Event loop.
  
  AliESDEvent *event = fTender->GetEvent();
  if (!event) {
    AliError("ESD event ptr = 0, returning");
    return;
  }
  
  // Initialising parameters once per run number
  if (fTender->RunChanged()){ 
    
    AliWarning( "Run changed, initializing parameters" );

    // get pass
    GetPass();

    // define what recalib parameters are needed for various switches
    // this is based on implementation in AliEMCALRecoUtils
    Bool_t needRecoParam   = fReClusterize;
    Bool_t needBadChannels = fBadCellRemove   | fClusterBadChannelCheck | fRecalDistToBadChannels | fReClusterize;
    Bool_t needRecalib     = fCalibrateEnergy | fReClusterize;
    Bool_t needTimecalib   = fCalibrateTime   | fReClusterize;
    Bool_t needMisalign    = fRecalClusPos    | fReClusterize;
    Bool_t needClusterizer = fReClusterize;

    // initiate reco params from OCDB or load some defaults on OCDB failure
    // will not overwrive, if those have been provided by user
    if( needRecoParam ){
      Int_t initRC = InitRecParam();
      
      if( initRC == 0 )
        AliError("Reco params load from OCDB failed! Defaults loaded.");
      if( initRC == 1 )
        AliWarning("Reco params loaded from OCDB.");
      if( initRC > 1 )
        AliWarning("Reco params not loaded from OCDB (user defined or previously failed to load).");
    }

    // Init bad channels
    if( needBadChannels ){
      Int_t fInitBC = InitBadChannels();
      if (fInitBC==0)
        AliError("InitBadChannels returned false, returning");
      if (fInitBC==1)
        AliWarning("InitBadChannels OK");
      if (fInitBC>1)
        AliWarning(Form("No external hot channel set: %d - %s", event->GetRunNumber(), fFilepass.Data()));
    }

    // init recalibration factors
    if( needRecalib ) { 
      Int_t fInitRecalib = InitRecalib();
      if (fInitRecalib==0)
        AliError("InitRecalib returned false, returning");
      if (fInitRecalib==1)
        AliWarning("InitRecalib OK");
      if (fInitRecalib>1)
        AliWarning(Form("No recalibration available: %d - %s", event->GetRunNumber(), fFilepass.Data()));
        fReCalibCluster = kFALSE;
    }
    
    // init time calibration
    if( needTimecalib ){
      Int_t initTC = InitTimeCalibration();
      if ( !initTC ) 
        AliError("InitTimeCalibration returned false, returning");
      if (initTC==1)
        AliWarning("InitTimeCalib OK");
      if( initTC > 1 )
        AliWarning(Form("No external time calibration set: %d - %s", event->GetRunNumber(), fFilepass.Data()));
    }

    // init misalignment matrix
    if( needMisalign ) { 
      if (!InitMisalignMatrix())
        AliError("InitMisalignmentMatrix returned false, returning");
      else
        AliWarning("InitMisalignMatrix OK");
    }
    
    // init clusterizer
    if( needClusterizer ) {
      if (!InitClusterization()) 
        AliError("InitClusterization returned false, returning");
      else
        AliWarning("InitClusterization OK");
    }
    
    if(fDebugLevel>1) 
      fEMCALRecoUtils->Print("");
  }
  
  // disable implied switches -------------------------------------------------
  // AliEMCALRecoUtils or clusterizer functions alredy include some
  // recalibration so based on those implied callibration te switches
  // here are disabled to avoid duplication
    
  // clusterizer does cluster energy recalibration, position recomputation
  // and shower shape
  if( fReClusterize ){
    fReCalibCluster   = kFALSE;
    fRecalClusPos     = kFALSE;
    fRecalShowerShape = kFALSE;
  }
  
  // CONFIGURE THE RECO UTILS -------------------------------------------------
  // configure the reco utils
  // this option does energy recalibration
  if( fCalibrateEnergy )
    fEMCALRecoUtils->SwitchOnRecalibration();
  else
    fEMCALRecoUtils->SwitchOffRecalibration();
  
  // allows time calibration
  if( fCalibrateTime )
    fEMCALRecoUtils->SwitchOnTimeRecalibration();
  else
    fEMCALRecoUtils->SwitchOffTimeRecalibration();

  // allows to zero bad cells
  if( fBadCellRemove )
    fEMCALRecoUtils->SwitchOnBadChannelsRemoval();
  else
    fEMCALRecoUtils->SwitchOffBadChannelsRemoval();
  
  // distance to bad channel recalibration
  if( fRecalDistToBadChannels )
    fEMCALRecoUtils->SwitchOnDistToBadChannelRecalculation();
  else
    fEMCALRecoUtils->SwitchOffDistToBadChannelRecalculation();

  // exclude exotic cells
  if( fRejectExoticCells )
    fEMCALRecoUtils->SwitchOnRejectExoticCell();
  else
    fEMCALRecoUtils->SwitchOffRejectExoticCell();
  
  // exclude clusters with exotic cells
  if( fRejectExoticClusters )
    fEMCALRecoUtils->SwitchOnRejectExoticCluster();
  else
    fEMCALRecoUtils->SwitchOffRejectExoticCluster();

  // TODO: not implemented switches
  // SwitchOnClusterEnergySmearing
  // SwitchOnRunDepCorrection

  // START PROCESSING ---------------------------------------------------------
  // Test if cells present
  AliESDCaloCells *cells= event->GetEMCALCells();
  if (cells->GetNumberOfCells()<=0) 
  {
    if(fDebugLevel>1) 
      AliWarning(Form("Number of EMCAL cells = %d, returning", cells->GetNumberOfCells()));
    return;
  }
  
  if (fDebugLevel>2)
    AliInfo(Form("Re-calibrate cluster %d\n",fReCalibCluster));

  // mark the cells not recalibrated in case of selected
  // time, energy recalibration or bad channel removal
  if( fCalibrateEnergy || fCalibrateTime || fBadCellRemove )
    fEMCALRecoUtils->ResetCellsCalibrated();
  
  // CELL RECALIBRATION -------------------------------------------------------
  // cell objects will be updated
  // the cell calibrations are also processed locally any time those are needed
  // in case that the cell objects are not to be updated here for later use
  if( fUpdateCell )
  {
    // do the update
    // includes exotic cell check in the UpdateCells function - is not provided
    // by the reco utils
    UpdateCells();

    // switch off recalibrations so those are not done multiple times
    // this is just for safety, the recalibrated flag of cell object
    // should not allow for farther processing anyways
    fEMCALRecoUtils->SwitchOffRecalibration();
    fEMCALRecoUtils->SwitchOffTimeRecalibration();  
  
    if (fDoUpdateOnly)
      return;
  }

  // RECLUSTERIZATION ---------------------------------------------------------
  if (fReClusterize)
  {
    FillDigitsArray();
    Clusterize();
    UpdateClusters();
  }
  
  // Store good clusters
  TClonesArray *clusArr = dynamic_cast<TClonesArray*>(event->FindListObject("caloClusters"));
  if (!clusArr) 
    clusArr = dynamic_cast<TClonesArray*>(event->FindListObject("CaloClusters"));
  if (!clusArr) {
    AliWarning(Form("No cluster array, number of cells in event = %d, returning", cells->GetNumberOfCells()));
    return;
  }

  // PROCESS SINGLE CLUSTER RECALIBRATION -------------------------------------
  // now go through clusters one by one and process separate correction
  // as those were defined or not
  Int_t nclusters = clusArr->GetEntriesFast();
  for (Int_t icluster=0; icluster < nclusters; ++icluster) 
  { 
    AliVCluster *clust = static_cast<AliVCluster*>(clusArr->At(icluster));
    if (!clust) 
      continue;
    if  (!clust->IsEMCAL()) 
      continue;

    // REMOVE CLUSTERS WITH BAD CELLS -----------------------------
    if( fClusterBadChannelCheck )
    {
      // careful, the the ClusterContainsBadChannel is dependent on
      // SwitchOnBadChannelsRemoval, switching it ON automatically
      // and returning to original value after processing
      Bool_t badRemoval = fEMCALRecoUtils->IsBadChannelsRemovalSwitchedOn();
      fEMCALRecoUtils->SwitchOnBadChannelsRemoval();
      
      Bool_t badResult = fEMCALRecoUtils->ClusterContainsBadChannel(fEMCALGeo, clust->GetCellsAbsId(), clust->GetNCells());

      // switch the bad channels removal back
      if( ! badRemoval )
        fEMCALRecoUtils->SwitchOffBadChannelsRemoval();
      
      if( badResult )
      {
        delete clusArr->RemoveAt(icluster);
        continue; //TODO is it really needed to remove it? Or should we flag it?
      }
    }
    
    // REMOVE EXOTIC CLUSTERS -------------------------------------
    // does process local cell recalibration energy and time without replacing
    // the global cell values, in case of no cell recalib done yet
    if( fRejectExoticClusters )
    {
      // careful, the IsExoticCluster is dependent on
      // SwitchOnRejectExoticCell, switching it ON automatically
      // and returning to original value after processing
      Bool_t exRemoval = fEMCALRecoUtils->IsRejectExoticCell();
      fEMCALRecoUtils->SwitchOnRejectExoticCell();

      // get bunch crossing
      Int_t bunchCrossNo = fTender->GetEvent()->GetBunchCrossNumber();

      Bool_t exResult = fEMCALRecoUtils->IsExoticCluster(clust, cells, bunchCrossNo );

      // switch the exotic channels removal back
      if( ! exRemoval )
        fEMCALRecoUtils->SwitchOffRejectExoticCell();
      
      if( exResult )
      {
        delete clusArr->RemoveAt(icluster);
        continue; //TODO is it really needed to remove it? Or should we flag it?
      }
    }
    
    // FIDUCIAL CUT -----------------------------------------------
    if (fFiducial)
    {
      // depends on SetNumberOfCellsFromEMCALBorder
      // SwitchOnNoFiducialBorderInEMCALEta0
      if (!fEMCALRecoUtils->CheckCellFiducialRegion(fEMCALGeo, clust, cells)){
        delete clusArr->RemoveAt(icluster);
        continue; //TODO it would be nice to store the distance
      }
    }
    
    // CLUSTER ENERGY ---------------------------------------------
    // does process local cell recalibration energy and time without replacing
    // the global cell values, in case of no cell recalib done yet
    if( fReCalibCluster ) 
      fEMCALRecoUtils->RecalibrateClusterEnergy(fEMCALGeo, clust, cells);

    // CLUSTER POSITION -------------------------------------------
    // does process local cell energy recalibration, if enabled and cells
    // not calibratied yet
    if( fRecalClusPos ) 
      fEMCALRecoUtils->RecalculateClusterPosition(fEMCALGeo, cells, clust);
    
    // SHOWER SHAPE -----------------------------------------------
    if( fRecalShowerShape )
      fEMCALRecoUtils->RecalculateClusterShowerShapeParameters(fEMCALGeo, cells, clust);  

    // NONLINEARITY -----------------------------------------------
    if( fDoNonLinearity )
    {
      Float_t correctedEnergy = fEMCALRecoUtils->CorrectClusterEnergyLinearity(clust);
      clust->SetE(correctedEnergy);
    }

    // DISTANCE TO BAD CHANNELS -----------------------------------
    if( fRecalDistToBadChannels )
      fEMCALRecoUtils->RecalculateClusterDistanceToBadChannel(fEMCALGeo, cells, clust);  
  }

  clusArr->Compress();

  if (!fDoTrackMatch)
    return;

  // TRACK MATCHING -----------------------------------------------------------
  if (!TGeoGlobalMagField::Instance()->GetField()) 
  {
    event->InitMagneticField();
  }
  
  fEMCALRecoUtils->FindMatches(event,0x0,fEMCALGeo);
  fEMCALRecoUtils->SetClusterMatchedToTrack(event);
  fEMCALRecoUtils->SetTracksMatchedToCluster(event);
}

//_____________________________________________________
Bool_t AliEMCALTenderSupply::InitMisalignMatrix()
{
  // Initialising misalignment matrices
  
  AliESDEvent *event = fTender->GetEvent();
  if (!event) 
    return kFALSE;
  
  if (fGeomMatrixSet) 
  {
    AliInfo("Misalignment matrix already set");  
    return kTRUE;
  }
  
  if (fDebugLevel>0) 
    AliInfo("Initialising misalignment matrix");  
  
  if (fLoadGeomMatrices) {
    for(Int_t mod=0; mod < fEMCALGeo->GetNumberOfSuperModules(); ++mod)
    {
      if (fEMCALMatrix[mod]){
        if(fDebugLevel > 2) 
          fEMCALMatrix[mod]->Print();
        fEMCALGeo->SetMisalMatrix(fEMCALMatrix[mod],mod);  
      }
    }
    fGeomMatrixSet = kTRUE;
    return kTRUE;
  }
  
  Int_t runGM = event->GetRunNumber();
  TObjArray *mobj = 0;

 if(fMisalignSurvey == kdefault)
 { //take default alignment corresponding to run no
    AliOADBContainer emcalgeoCont(Form("emcal"));
    emcalgeoCont.InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALlocal2master.root",Form("AliEMCALgeo"));
    mobj=(TObjArray*)emcalgeoCont.GetObject(runGM,"EmcalMatrices");
 }

 if(fMisalignSurvey == kSurveybyS)
 { //take alignment at sector level
  if (runGM <= 140000) { //2010 data
    AliOADBContainer emcalgeoCont(Form("emcal2010"));
    emcalgeoCont.InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALlocal2master.root",Form("AliEMCALgeo"));
    mobj=(TObjArray*)emcalgeoCont.GetObject(100,"survey10");
    
  } else if (runGM>140000)
  { // 2011 LHC11a pass1 data
    AliOADBContainer emcalgeoCont(Form("emcal2011"));
    emcalgeoCont.InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALlocal2master.root",Form("AliEMCALgeo"));
    mobj=(TObjArray*)emcalgeoCont.GetObject(100,"survey11byS");      
  }
 }

 if(fMisalignSurvey == kSurveybyM)
 { //take alignment at module level
  if (runGM <= 140000) { //2010 data
    AliOADBContainer emcalgeoCont(Form("emcal2010"));
    emcalgeoCont.InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALlocal2master.root",Form("AliEMCALgeo"));
    mobj=(TObjArray*)emcalgeoCont.GetObject(100,"survey10");
    
  } else if (runGM>140000) 
  { // 2011 LHC11a pass1 data
    AliOADBContainer emcalgeoCont(Form("emcal2011"));
    emcalgeoCont.InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALlocal2master.root",Form("AliEMCALgeo"));
    mobj=(TObjArray*)emcalgeoCont.GetObject(100,"survey11byM");      
  }
 }

  if(!mobj)
  {
    AliFatal("Geometry matrix array not found");
    return kFALSE;
  }
  
 for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
 {
   fEMCALMatrix[mod] = (TGeoHMatrix*) mobj->At(mod);
   fEMCALGeo->SetMisalMatrix(fEMCALMatrix[mod],mod); 
   fEMCALMatrix[mod]->Print();
 }
  
  return kTRUE;
}

//_____________________________________________________
Int_t AliEMCALTenderSupply::InitBadChannels()
{
  // Initialising bad channel maps
  AliESDEvent *event = fTender->GetEvent();
  if (!event) 
    return 0;
  
  if (fDebugLevel>0) 
    AliInfo("Initialising Bad channel map");
  
  // init default maps first
  if( !fEMCALRecoUtils->GetEMCALBadChannelStatusMapArray() )
    fEMCALRecoUtils->InitEMCALBadChannelStatusMap() ;
  
  Int_t runBC = event->GetRunNumber();
  
  AliOADBContainer *contBC = new AliOADBContainer("");
  if (fBasePath!="")
  { //if fBasePath specified in the ->SetBasePath()
    if (fDebugLevel>0) AliInfo(Form("Loading Bad Channels OADB from given path %s",fBasePath.Data()));
    
    TFile *fbad=new TFile(Form("%s/EMCALBadChannels.root",fBasePath.Data()),"read");
    if (!fbad || fbad->IsZombie())
    {
      AliFatal(Form("EMCALBadChannels.root was not found in the path provided: %s",fBasePath.Data()));
      return 0;
    }  
    
    if (fbad) delete fbad;
    
    contBC->InitFromFile(Form("%s/EMCALBadChannels.root",fBasePath.Data()),"AliEMCALBadChannels");    
  } 
  else 
  { // Else choose the one in the $ALICE_ROOT directory
    if (fDebugLevel>0) AliInfo("Loading Bad Channels OADB from $ALICE_ROOT/OADB/EMCAL");
    
    TFile *fbad=new TFile("$ALICE_ROOT/OADB/EMCAL/EMCALBadChannels.root","read");
    if (!fbad || fbad->IsZombie())
    {
      AliFatal("$ALICE_ROOT/OADB/EMCAL/EMCALBadChannels.root was not found");
      return 0;
    }  
      
    if (fbad) delete fbad;
    
    contBC->InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALBadChannels.root","AliEMCALBadChannels"); 
  }
  
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(runBC);
  if (!arrayBC)
  {
    AliError(Form("No external hot channel set for run number: %d", runBC));
    return 2; 
  }
  
  TObjArray *arrayBCpass=(TObjArray*)arrayBC->FindObject(fFilepass); // Here, it looks for a specific pass
  if (!arrayBCpass)
  {
    AliError(Form("No external hot channel set for: %d -%s", runBC,fFilepass.Data()));
    return 2; 
  }

  if (fDebugLevel>0) arrayBCpass->Print();

  Int_t sms = fEMCALGeo->GetEMCGeometry()->GetNumberOfSuperModules();
  for (Int_t i=0; i<sms; ++i) 
  {
    TH2I *h = fEMCALRecoUtils->GetEMCALChannelStatusMap(i);
    if (h)
      delete h;
    h=(TH2I*)arrayBCpass->FindObject(Form("EMCALBadChannelMap_Mod%d",i));

    if (!h) 
    {
      AliError(Form("Can not get EMCALBadChannelMap_Mod%d",i));
      continue;
    }
    h->SetDirectory(0);
    fEMCALRecoUtils->SetEMCALChannelStatusMap(i,h);
  }
  return 1;  
}

//_____________________________________________________
Int_t AliEMCALTenderSupply::InitRecalib()
{
  // Initialising recalibration factors.
  
  AliESDEvent *event = fTender->GetEvent();
  if (!event) 
    return 0;
  
  if (fDebugLevel>0) 
    AliInfo("Initialising recalibration factors");
  
  // init default maps first
  if( !fEMCALRecoUtils->GetEMCALRecalibrationFactorsArray() )
    fEMCALRecoUtils->InitEMCALRecalibrationFactors() ;

  Int_t runRC = event->GetRunNumber();
      
  AliOADBContainer *contRF=new AliOADBContainer("");
  if (fBasePath!="") 
  { //if fBasePath specified in the ->SetBasePath()
    if (fDebugLevel>0)  AliInfo(Form("Loading Recalib OADB from given path %s",fBasePath.Data()));
    
    TFile *fRecalib= new TFile(Form("%s/EMCALRecalib.root",fBasePath.Data()),"read");
    if (!fRecalib || fRecalib->IsZombie()) 
    {
      AliFatal(Form("EMCALRecalib.root not found in %s",fBasePath.Data()));
      return 0;
    }
    
    if (fRecalib) delete fRecalib;
    
    contRF->InitFromFile(Form("%s/EMCALRecalib.root",fBasePath.Data()),"AliEMCALRecalib");
  }
    else
    { // Else choose the one in the $ALICE_ROOT directory
      if (fDebugLevel>0)  AliInfo("Loading Recalib OADB from $ALICE_ROOT/OADB/EMCAL");
      
      TFile *fRecalib= new TFile("$ALICE_ROOT/OADB/EMCAL/EMCALRecalib.root","read");
      if (!fRecalib || fRecalib->IsZombie()) 
      {
        AliFatal("$ALICE_ROOT/OADB/EMCAL/EMCALRecalib.root was not found");
        return 0;
      }
      
      if (fRecalib) delete fRecalib;
      
      contRF->InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALRecalib.root","AliEMCALRecalib");     
    }

  TObjArray *recal=(TObjArray*)contRF->GetObject(runRC); //GetObject(int runnumber)
  if (!recal)
  {
    AliError(Form("No Objects for run: %d",runRC));
    return 2;
  } 

  TObjArray *recalpass=(TObjArray*)recal->FindObject(fFilepass);
  if (!recalpass)
  {
    AliError(Form("No Objects for run: %d - %s",runRC,fFilepass.Data()));
    return 2;
  }

  TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
  if (!recalib)
  {
    AliError(Form("No Recalib histos found for  %d - %s",runRC,fFilepass.Data())); 
    return 2;
  }

  if (fDebugLevel>0) recalib->Print();

  Int_t sms = fEMCALGeo->GetEMCGeometry()->GetNumberOfSuperModules();
  for (Int_t i=0; i<sms; ++i) 
  {
    TH2F *h = fEMCALRecoUtils->GetEMCALChannelRecalibrationFactors(i);
    if (h)
      delete h;
    h = (TH2F*)recalib->FindObject(Form("EMCALRecalFactors_SM%d",i));
    if (!h) 
    {
      AliError(Form("Could not load EMCALRecalFactors_SM%d",i));
      continue;
    }
    h->SetDirectory(0);
    fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(i,h);
  }
  return 1;
}

//_____________________________________________________
Int_t AliEMCALTenderSupply::InitTimeCalibration()
{
  // Initialising bad channel maps
  AliESDEvent *event = fTender->GetEvent();
  if (!event) 
    return 0;
  
  if (fDebugLevel>0) 
    AliInfo("Initialising time calibration map");
  
  // init default maps first
  if( !fEMCALRecoUtils->GetEMCALTimeRecalibrationFactorsArray() )
    fEMCALRecoUtils->InitEMCALTimeRecalibrationFactors() ;

  Int_t runBC = event->GetRunNumber();
  
  AliOADBContainer *contBC = new AliOADBContainer("");
  if (fBasePath!="")
  { //if fBasePath specified in the ->SetBasePath()
    if (fDebugLevel>0) AliInfo(Form("Loading time calibration OADB from given path %s",fBasePath.Data()));
    
    TFile *fbad=new TFile(Form("%s/EMCALTimeCalib.root",fBasePath.Data()),"read");
    if (!fbad || fbad->IsZombie())
    {
      AliFatal(Form("EMCALTimeCalib.root was not found in the path provided: %s",fBasePath.Data()));
      return 0;
    }  
    
    if (fbad) delete fbad;
    
    contBC->InitFromFile(Form("%s/EMCALTimeCalib.root",fBasePath.Data()),"AliEMCALTimeCalib");    
  } 
  else 
  { // Else choose the one in the $ALICE_ROOT directory
    if (fDebugLevel>0) AliInfo("Loading time calibration OADB from $ALICE_ROOT/OADB/EMCAL");
    
    TFile *fbad=new TFile("$ALICE_ROOT/OADB/EMCAL/EMCALTimeCalib.root","read");
    if (!fbad || fbad->IsZombie())
    {
      AliFatal("$ALICE_ROOT/OADB/EMCAL/EMCALTimeCalib.root was not found");
      return 0;
    }  
      
    if (fbad) delete fbad;
    
    contBC->InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALTimeCalib.root","AliEMCALTimeCalib"); 
  }
  
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(runBC);
  if (!arrayBC)
  {
    AliError(Form("No external time calibration set for run number: %d", runBC));
    return 2; 
  }
  
  TObjArray *arrayBCpass=(TObjArray*)arrayBC->FindObject(fFilepass); // Here, it looks for a specific pass
  if (!arrayBCpass)
  {
    AliError(Form("No external time calibration set for: %d -%s", runBC,fFilepass.Data()));
    return 2; 
  }

  if (fDebugLevel>0) arrayBCpass->Print();

  for( Int_t i = 0; i < 4; i++ )
  {
    TH1F *h = fEMCALRecoUtils->GetEMCALChannelTimeRecalibrationFactors( i );
    if( h )
      delete h;
    
    h = (TH1F*)arrayBCpass->FindObject(Form("hAllTimeAvBC%d",i));
    
    if (!h)
    {
      AliError(Form("Can not get hAllTimeAvBC%d",i));
      continue;
    }
    h->SetDirectory(0);
    fEMCALRecoUtils->SetEMCALChannelTimeRecalibrationFactors(i,h);
  }
  return 1;  
}

//_____________________________________________________
void AliEMCALTenderSupply::UpdateCells()
{
  //Remove bad cells from the cell list
  //Recalibrate energy and time cells 
  //This is required for later reclusterization

  AliESDCaloCells *cells = fTender->GetEvent()->GetEMCALCells();
  Int_t bunchCrossNo = fTender->GetEvent()->GetBunchCrossNumber();

  fEMCALRecoUtils->RecalibrateCells(cells, bunchCrossNo); 
  
  // remove exotic cells - loop through cells and zero the exotic ones
  // just like with bad cell rejection in reco utils (inside RecalibrateCells)
  if( fRejectExoticCells )
  {
    Int_t    absId  =-1;
    Bool_t   isExot = kFALSE;
  
    // loop through cells
    Int_t nEMcell  = cells->GetNumberOfCells() ;  
    for (Int_t iCell = 0; iCell < nEMcell; iCell++) 
    { 
      absId  = cells->GetCellNumber(iCell);
    
      isExot = fEMCALRecoUtils->IsExoticCell( absId, cells, bunchCrossNo ); 
      // zero if exotic
      if( isExot )
        cells->SetCell( iCell, absId, 0.0, 0.0 );
    } // cell loop
  } // reject exotic cells

  cells->Sort();
}

//_____________________________________________________
Int_t AliEMCALTenderSupply::InitRecParam()
{
  // Initalize the reconstructor parameters from OCDB
  // load some default on OCDB failure
  
  Int_t runNum = -1;
  AliCDBManager    * man   = 0x0 ;
  TObjArray        * arr   = 0x0 ;
  AliEMCALRecParam * pars  = 0x0 ;
  AliCDBEntry      * entry = 0x0;
  const AliESDRun  * run   = 0x0 ;
  TString beamType ;
  
  // clean the previous reco params, if those came from OCDB
  // we do not want to erase user provided params, do we
  if( fRecoParamsOCDBLoaded ){
    if( fRecParam != 0 ){
      delete fRecParam;
      fRecParam = 0;
    }
    // zero the OCDB loaded flag
    fRecoParamsOCDBLoaded = kFALSE;
  }
  
  // exit if reco params exist (probably shipped by the user already)
  if( fRecParam != 0 )
    return 2;
  
  if (fDebugLevel>0) 
    AliInfo("Initialize the recParam");

  // get run details
  run      = fTender->GetEvent()->GetESDRun();
  beamType = run->GetBeamType();
  runNum   = fTender->GetEvent()->GetRunNumber();

  // OCDB manager should already exist
  // and have a default storage defined (done by AliTender)
  man = AliCDBManager::Instance();

  // load the file data
  if(man)
    entry = man->Get("EMCAL/Calib/RecoParam", runNum);
  
  if( entry )
    arr = (TObjArray*)(entry->GetObject());
  
  if( arr ){
    // load given parameters based on beam type
    if( beamType == "A-A" ){
      if( fDebugLevel > 0 )
        AliInfo( "Initializing A-A reco params." );
      pars = (AliEMCALRecParam*)arr->FindObject( "High Flux - Pb+Pb" );
    }
    else{
      if( fDebugLevel > 0 )
        AliInfo( "Initializing p-p reco params." );
      pars = (AliEMCALRecParam*)arr->FindObject( "Low Flux - p+p" );
    }
    
    // set the parameters, if found    
    if( pars ){
      if( fDebugLevel > 0 )
        AliInfo( "OCDB reco params set." );
      
      fRecParam = pars;
      fRecoParamsOCDBLoaded = kTRUE;
    }
    
    arr->Clear();
    delete arr;
  }

  // set some defaults if OCDB did not succede
  if( !fRecoParamsOCDBLoaded ){
    fRecParam = new AliEMCALRecParam();
    
    if ( beamType == "A-A"){
      fRecParam->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv2);
      fRecParam->SetClusteringThreshold(0.100);
      fRecParam->SetMinECut(0.050);
    } 
    else 
    {
      fRecParam->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv1);
      fRecParam->SetClusteringThreshold(0.100);
      fRecParam->SetMinECut(0.050);
    }
  }
  
  if( fRecoParamsOCDBLoaded )
    return 1;
  else
    return 0;
}

//_____________________________________________________
Bool_t AliEMCALTenderSupply::InitClusterization()
{
  // Initialing clusterization/unfolding algorithm and set all the needed parameters.
  
  AliESDEvent *event = fTender->GetEvent();
  if (!event) 
    return kFALSE;
  
  if (fDebugLevel>0) 
    AliInfo(Form("Initialising reclustering parameters: Clusterizer type: %d",fRecParam->GetClusterizerFlag()));
  
  //---setup clusterizer
  delete fClusterizer;
  if     (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv1)
    fClusterizer = new AliEMCALClusterizerv1 (fEMCALGeo);
  else if (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv2) 
    fClusterizer = new AliEMCALClusterizerv2(fEMCALGeo);
  else if (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerNxN) 
  {
    AliEMCALClusterizerNxN *clusterizer = new AliEMCALClusterizerNxN(fEMCALGeo);
    clusterizer->SetNRowDiff(fRecParam->GetNRowDiff());
    clusterizer->SetNColDiff(fRecParam->GetNColDiff());
    fClusterizer = clusterizer;
  } 
  else 
  {
    AliFatal(Form("Clusterizer < %d > not available", fRecParam->GetClusterizerFlag()));
    return kFALSE;
  }
  
  // Set the clustering parameters
  fClusterizer->SetECAClusteringThreshold( fRecParam->GetClusteringThreshold() );
  fClusterizer->SetECALogWeight          ( fRecParam->GetW0()                  );
  fClusterizer->SetMinECut               ( fRecParam->GetMinECut()             );    
  fClusterizer->SetUnfolding             ( fRecParam->GetUnfold()              );
  fClusterizer->SetECALocalMaxCut        ( fRecParam->GetLocMaxCut()           );
  fClusterizer->SetTimeCut               ( fRecParam->GetTimeCut()             );
  fClusterizer->SetTimeMin               ( fRecParam->GetTimeMin()             );
  fClusterizer->SetTimeMax               ( fRecParam->GetTimeMax()             );
  fClusterizer->SetInputCalibrated       ( kTRUE                               );
  fClusterizer->SetJustClusters          ( kTRUE                               );  
  
  // In case of unfolding after clusterization is requested, set the corresponding parameters
  if (fRecParam->GetUnfold()) 
  {
    for (Int_t i = 0; i < 8; ++i) 
    {
      fClusterizer->SetSSPars(i, fRecParam->GetSSPars(i));
    }
    for (Int_t i = 0; i < 3; ++i)
    {
      fClusterizer->SetPar5  (i, fRecParam->GetPar5(i));
      fClusterizer->SetPar6  (i, fRecParam->GetPar6(i));
    }
    fClusterizer->InitClusterUnfolding();
  }
  
  fClusterizer->SetDigitsArr(fDigitsArr);
  fClusterizer->SetOutput(0);
  fClusterArr = const_cast<TObjArray *>(fClusterizer->GetRecPoints());
  return kTRUE;
}

//_____________________________________________________
void AliEMCALTenderSupply::FillDigitsArray()
{
  // Fill digits from cells to a TClonesArray.
  
  AliESDEvent *event = fTender->GetEvent();
  if (!event)
    return;
  
  fDigitsArr->Clear("C");
  AliESDCaloCells *cells = event->GetEMCALCells();
  Int_t ncells = cells->GetNumberOfCells();
  for (Int_t icell = 0, idigit = 0; icell < ncells; ++icell) 
  {
    Double_t cellAmplitude=0, cellTime=0;
    Short_t  cellNumber=0;

    if (cells->GetCell(icell, cellNumber, cellAmplitude, cellTime) != kTRUE)
      break;

    // Do not add if already too low (some cells set to 0 if bad channels)
    if (cellAmplitude < fRecParam->GetMinECut() ) 
      continue;

    // If requested, do not include exotic cells
   if (fEMCALRecoUtils->IsExoticCell(cellNumber,cells,event->GetBunchCrossNumber())) 
      continue;
        
    AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->New(idigit));
    digit->SetId(cellNumber);
    digit->SetTime(cellTime);
    digit->SetTimeR(cellTime);
    digit->SetIndexInList(idigit);
    digit->SetType(AliEMCALDigit::kHG);
    digit->SetAmplitude(cellAmplitude);
    idigit++;
  }
}

//_____________________________________________________
void AliEMCALTenderSupply::Clusterize()
{
  // Clusterize.
  
  fClusterizer->Digits2Clusters("");
}

//_____________________________________________________
void AliEMCALTenderSupply::UpdateClusters()
{
  // Update ESD cluster list.
  
  AliESDEvent *event = fTender->GetEvent();
  if (!event)
    return;
  
  TClonesArray *clus = dynamic_cast<TClonesArray*>(event->FindListObject("caloClusters"));
  if (!clus) 
    clus = dynamic_cast<TClonesArray*>(event->FindListObject("CaloClusters"));
  if (!clus) 
  {
    AliError(" Null pointer to calo clusters array, returning");
    return;
  }
  
  Int_t nents = clus->GetEntriesFast();
  for (Int_t i=0; i < nents; ++i) 
  {
    AliESDCaloCluster *c = static_cast<AliESDCaloCluster*>(clus->At(i));
    if (!c)
      continue;
    if (c->IsEMCAL())
      delete clus->RemoveAt(i);
  }
  
  clus->Compress();
  
  RecPoints2Clusters(clus);
  
}

//_____________________________________________________
void AliEMCALTenderSupply::RecPoints2Clusters(TClonesArray *clus)
{
  // Convert AliEMCALRecoPoints to AliESDCaloClusters.
  // Cluster energy, global position, cells and their amplitude fractions are restored.
  
  AliESDEvent *event = fTender->GetEvent();
  if (!event)
    return;

  Int_t ncls = fClusterArr->GetEntriesFast();
  for(Int_t i=0, nout=clus->GetEntriesFast(); i < ncls; ++i) 
  {
    AliEMCALRecPoint *recpoint = static_cast<AliEMCALRecPoint*>(fClusterArr->At(i));
    
    Int_t ncellsTrue = 0;
    const Int_t ncells = recpoint->GetMultiplicity();
    UShort_t   absIds[ncells];  
    Double32_t ratios[ncells];
    Int_t *dlist = recpoint->GetDigitsList();
    Float_t *elist = recpoint->GetEnergiesList();
    for (Int_t c = 0; c < ncells; ++c) 
    {
      AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->At(dlist[c]));
      absIds[ncellsTrue] = digit->GetId();
      ratios[ncellsTrue] = elist[c]/digit->GetAmplitude();
      if (ratios[ncellsTrue] < 0.001) 
        continue;
      ++ncellsTrue;
    }
    
    if (ncellsTrue < 1) 
    {
      AliWarning("Skipping cluster with no cells");
      continue;
    }
    
    // calculate new cluster position
    TVector3 gpos;
    recpoint->GetGlobalPosition(gpos);
    Float_t g[3];
    gpos.GetXYZ(g);
    
    AliESDCaloCluster *c = static_cast<AliESDCaloCluster*>(clus->New(nout++));
    c->SetID(nout-1); 
    c->SetType(AliVCluster::kEMCALClusterv1);
    c->SetE(recpoint->GetEnergy());
    c->SetPosition(g);
    c->SetNCells(ncellsTrue);
    c->SetDispersion(recpoint->GetDispersion());
    c->SetEmcCpvDistance(-1);            //not yet implemented
    c->SetChi2(-1);                      //not yet implemented
    c->SetTOF(recpoint->GetTime()) ;     //time-of-flight
    c->SetNExMax(recpoint->GetNExMax()); //number of local maxima
    Float_t elipAxis[2];
    recpoint->GetElipsAxis(elipAxis);
    c->SetM02(elipAxis[0]*elipAxis[0]) ;
    c->SetM20(elipAxis[1]*elipAxis[1]) ;
    AliESDCaloCluster *cesd = static_cast<AliESDCaloCluster*>(c);
    cesd->SetCellsAbsId(absIds);
    cesd->SetCellsAmplitudeFraction(ratios);
  }
}

//_____________________________________________________
void AliEMCALTenderSupply::GetPass()
{
  // Get passx from filename.
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  fInputTree = mgr->GetTree();
  
  if (!fInputTree) 
  {
    AliError("Pointer to tree = 0, returning");
    return;
  }
  
  fInputFile = fInputTree->GetCurrentFile();
  if (!fInputFile) {
    AliError("Null pointer input file, returning");
    return;
  }
  
  TString fname(fInputFile->GetName());
  if      (fname.Contains("pass1")) fFilepass = TString("pass1");
  else if (fname.Contains("pass2")) fFilepass = TString("pass2");
  else if (fname.Contains("pass3")) fFilepass = TString("pass3");
  else 
  {
    AliError(Form("Pass number string not found: %s", fname.Data()));
    return;            
  }
}

