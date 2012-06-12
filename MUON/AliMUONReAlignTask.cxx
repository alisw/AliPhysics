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

// $Id$

//-----------------------------------------------------------------------------
/// \class AliMUONReAlignTask
/// AliAnalysisTask to realign the MUON spectrometer.
/// The Task reads as input ESDs moves the clusters of a MUONTrack acoording
/// to the re aligned geometry taken from a misalignment file in the OCDB 
/// and refit the track. Then it produces a AliMUONClusterInfo object for each
/// cluster. The output is a TTree of AliMUONClusterInfo
///
/// \author Javier Castillo, CEA/Saclay - Irfu/SPhN
//-----------------------------------------------------------------------------

#include <fstream>

#include <TString.h>
#include <TError.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TRandom.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <Riostream.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliMagF.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliGeomManager.h"

#include "AliMpConstants.h"
#include "AliMpCDB.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpPad.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONPadInfo.h"
#include "AliMUONClusterInfo.h"
#include "AliMUONRecoParam.h"
#include "AliMUONESDInterface.h"
#include "AliMUONRefitter.h"
#include "AliMUONRecoParam.h"
#include "AliMUONVDigit.h"
#include "AliMUONVCluster.h"
#include "AliMUONTrack.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONGeometryTransformer.h"

#include "AliMUONReAlignTask.h"

using std::cout;
using std::endl;
///\cond CLASSIMP   
ClassImp(AliMUONReAlignTask)
///\endcond

//________________________________________________________________________
AliMUONReAlignTask::AliMUONReAlignTask(const char *name, const char *geofilename, const char *defaultocdb, const char *misalignocdb) 
  : AliAnalysisTask(name, ""),
    fESD(0x0),
    fClusterInfoTree(0x0),
    fClusterInfo(0x0),
    fESDInterface(0x0),
    fRefitter(0x0),
    fRecoParam(0x0),
    fGeoFilename(geofilename),
    fMisAlignOCDB(misalignocdb),
    fDefaultOCDB(defaultocdb),
    fGeoTransformer(0x0),
    fNewGeoTransformer(0x0),
    fGainStore(0x0),
    fPedStore(0x0),
    fPrintLevel(0),
    fLastRun(-1)
{
  /// Default Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes a TTree
  DefineOutput(0, TTree::Class());

  fClusterInfo = new AliMUONClusterInfo();
  fESDInterface = new AliMUONESDInterface();
  fGeoTransformer = new AliMUONGeometryTransformer();
  fNewGeoTransformer = new AliMUONGeometryTransformer();
}

//________________________________________________________________________
AliMUONReAlignTask::AliMUONReAlignTask(const AliMUONReAlignTask& obj)
  : AliAnalysisTask(obj),
    fESD(0x0),
    fClusterInfoTree(0x0),
    fClusterInfo(0x0),
    fESDInterface(0x0),
    fRefitter(0x0),
    fRecoParam(0x0),
    fGeoFilename(""),
    fMisAlignOCDB(""),
    fDefaultOCDB(""),
    fGeoTransformer(0x0),
    fNewGeoTransformer(0x0),
    fGainStore(0x0),
    fPedStore(0x0),
    fPrintLevel(0),
    fLastRun(-1)
{
  /// Copy constructor
  fESD = obj.fESD;
  fClusterInfoTree = obj.fClusterInfoTree;
  fClusterInfo = obj.fClusterInfo;
  fESDInterface = obj.fESDInterface;
  fRefitter = obj.fRefitter;
  fRecoParam = obj.fRecoParam;
  fGeoFilename = obj.fGeoFilename;
  fMisAlignOCDB = obj.fMisAlignOCDB;
  fDefaultOCDB = obj.fDefaultOCDB;
  fGeoTransformer = obj.fGeoTransformer;
  fNewGeoTransformer = obj.fNewGeoTransformer;
  fGainStore = obj.fGainStore;
  fPedStore = obj.fPedStore;
  fPrintLevel = obj.fPrintLevel;
  fLastRun = obj.fLastRun;
}

//________________________________________________________________________________
AliMUONReAlignTask& AliMUONReAlignTask::operator=(const AliMUONReAlignTask& other)
{
  /// Assignment
  if(&other == this) return *this;
  AliAnalysisTask::operator=(other);
  fESD = other.fESD;
  fClusterInfoTree = other.fClusterInfoTree;
  fClusterInfo = other.fClusterInfo;
  fESDInterface = other.fESDInterface;
  fRefitter = other.fRefitter;
  fRecoParam = other.fRecoParam;
  fGeoFilename = other.fGeoFilename;
  fMisAlignOCDB = other.fMisAlignOCDB;
  fDefaultOCDB = other.fDefaultOCDB;
  fGeoTransformer = other.fGeoTransformer;
  fNewGeoTransformer = other.fNewGeoTransformer;
  fGainStore = other.fGainStore;
  fPedStore = other.fPedStore;
  fPrintLevel = other.fPrintLevel;
  fLastRun = other.fLastRun;

  return *this;
} 


//________________________________________________________________________
AliMUONReAlignTask::~AliMUONReAlignTask() 
{ 
  /// Destructor
  if (fESDInterface) delete fESDInterface;
  if (fGeoTransformer) delete fGeoTransformer;
  if (fNewGeoTransformer) delete fNewGeoTransformer;
}

//________________________________________________________________________
void AliMUONReAlignTask::LocalInit() 
{
  /// Local initialization, called once per task on the client machine 
  /// where the analysis train is assembled
  AliMpCDB::LoadMpSegmentation();

  // prepare the refitting
  gRandom->SetSeed(0);
  Prepare(fGeoFilename.Data(),fDefaultOCDB.Data(),fMisAlignOCDB.Data());

  fRefitter = new AliMUONRefitter(fRecoParam);
  fRefitter->Connect(fESDInterface);
  
  // Original geotransformer
  fGeoTransformer->LoadGeometryData();   
  // Apply mis alignments
  AliGeomManager::ApplyAlignObjsFromCDB("MUON");    
  // New geotransformer
  fNewGeoTransformer->LoadGeometryData();   
  
}

//________________________________________________________________________
void AliMUONReAlignTask::ConnectInputData(Option_t *) 
{
  /// Connect ESD here. Called on each input data change.
  // Connect ESD here
  TTree* esdTree = dynamic_cast<TTree*> (GetInputData(0));
  if (!esdTree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());   
    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } 
    else {
      fESD = esdH->GetEvent();
    }
  }
}

//________________________________________________________________________
void AliMUONReAlignTask::CreateOutputObjects()
{
  /// Executed once on each worker (machine actually running the analysis code)
  //
  // This method has to be called INSIDE the user redefined CreateOutputObjects
  // method, before creating each object corresponding to the output containers
  // that are to be written to a file. This need to be done in general for the big output
  // objects that may not fit memory during processing. 
  OpenFile(0); 
  
  fClusterInfoTree = new TTree("clusterInfoTree","clusterInfoTree");
  fClusterInfoTree->Branch("clusterInfo", &fClusterInfo, 32000, 99);

}

//________________________________________________________________________
void AliMUONReAlignTask::Exec(Option_t *) 
{
  /// Main loop, called for each event
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  // Prepare Gain and Pedestal store
  if (!fGainStore || fLastRun!=fESD->GetRunNumber()){
    fGainStore = AliMUONCalibrationData::CreateGains(fESD->GetRunNumber());
    fLastRun = fESD->GetRunNumber();
  }
  if (!fPedStore || fLastRun!=fESD->GetRunNumber()){
    fPedStore = AliMUONCalibrationData::CreatePedestals(fESD->GetRunNumber());
    fLastRun = fESD->GetRunNumber();
  }

  AliMUONPadInfo padInfo;

  Int_t nTracks = (Int_t)fESD->GetNumberOfMuonTracks();
  if (nTracks < 1) return;
  
  // load the current event
  fESDInterface->LoadEvent(*fESD);
  AliMUONVDigitStore* digitStore = fESDInterface->GetDigits();

  Double_t lX = 0.;
  Double_t lY = 0.;
  Double_t lZ = 0.;
  Double_t gX = 0.;
  Double_t gY = 0.;
  Double_t gZ = 0.;
  // loop over cluster to modify their position
  AliMUONVCluster *cluster;
  TIter next(fESDInterface->CreateClusterIterator());
  while ((cluster = static_cast<AliMUONVCluster*>(next()))) {
    //       cout << "Original cluster" << endl;
    //       cluster->Print();
    gX = cluster->GetX();
    gY = cluster->GetY();
    gZ = cluster->GetZ();
    fGeoTransformer->Global2Local(cluster->GetDetElemId(),gX,gY,gZ,lX,lY,lZ);
    fNewGeoTransformer->Local2Global(cluster->GetDetElemId(),lX,lY,lZ,gX,gY,gZ);
    cluster->SetXYZ(gX,gY,gZ);
    //       cout << "Aligned cluster" << endl;
    //       cluster->Print();
  }

  // refit the tracks from digits
  AliMUONVTrackStore* newTrackStore = fRefitter->ReconstructFromClusters();
    
  //----------------------------------------------//
  // ------ fill new ESD and print results ------ //
  //----------------------------------------------//
  // loop over the list of ESD tracks
  TClonesArray *esdTracks = (TClonesArray*) fESD->FindListObject("MuonTracks");
  for (Int_t iTrack = 0; iTrack <  nTracks; iTrack++) {      
    // get the ESD track
    AliESDMuonTrack* esdTrack = (AliESDMuonTrack*) esdTracks->UncheckedAt(iTrack);
    
    // skip ghost tracks (leave them unchanged in the new ESD file)
    if (!esdTrack->ContainTrackerData()) continue;
      
    // get the corresponding MUON track
    AliMUONTrack* track = fESDInterface->FindTrack(esdTrack->GetUniqueID());
      
    // Find the corresponding re-fitted MUON track
    AliMUONTrack* newTrack = (AliMUONTrack*) newTrackStore->FindObject(esdTrack->GetUniqueID());
    
    // print initial and re-fitted track parameters at first cluster if any
    if (fPrintLevel>0) {
      cout<<"            ----------------track #"<<iTrack+1<<"----------------"<<endl;
      cout<<"before refit:"<<endl;
      AliMUONTrackParam *param = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->First();
      param->Print("FULL");
      if (fPrintLevel>1) param->GetCovariances().Print();
      if (!newTrack) continue;
      cout<<"after refit:"<<endl;
      param = (AliMUONTrackParam*) newTrack->GetTrackParamAtCluster()->First();
      param->Print("FULL");
      if (fPrintLevel>1) param->GetCovariances().Print();
      cout<<"            ----------------------------------------"<<endl;
    }
    
    // Cluster Info part
    UInt_t muonClusterMap = BuildClusterMap(*newTrack);
	
    // loop over clusters
    AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(newTrack->GetTrackParamAtCluster()->First());
    while (trackParam) {
      fClusterInfo->Clear("C");
      
      // fill cluster info
      cluster = trackParam->GetClusterPtr();
      fClusterInfo->SetRunId(fESD->GetRunNumber());
      fClusterInfo->SetEventId(fESD->GetEventNumberInFile());
      fClusterInfo->SetZ(cluster->GetZ());
      fClusterInfo->SetClusterId(cluster->GetUniqueID());
      fClusterInfo->SetClusterXY(cluster->GetX(), cluster->GetY());
      fClusterInfo->SetClusterXYErr(cluster->GetErrX(), cluster->GetErrY());
      fClusterInfo->SetClusterChi2(cluster->GetChi2());
      fClusterInfo->SetClusterCharge(cluster->GetCharge());
      
      // fill track info
      fClusterInfo->SetTrackId(newTrack->GetUniqueID());
      fClusterInfo->SetTrackXY(trackParam->GetNonBendingCoor(), trackParam->GetBendingCoor());
      fClusterInfo->SetTrackThetaXY(TMath::ATan(trackParam->GetNonBendingSlope()), TMath::ATan(trackParam->GetBendingSlope()));
      fClusterInfo->SetTrackP(trackParam->P());
      const TMatrixD paramCov = trackParam->GetCovariances();
      fClusterInfo->SetTrackXYErr(TMath::Sqrt(paramCov(0,0)), TMath::Sqrt(paramCov(2,2)));
      fClusterInfo->SetTrackChi2(newTrack->GetNormalizedChi2());
      fClusterInfo->SetTrackCharge((Short_t)trackParam->GetCharge());
      fClusterInfo->SetTrackNHits(newTrack->GetNClusters());
      fClusterInfo->SetTrackChamberHitMap(muonClusterMap);
      
      // fill pad info if available	  
      for (Int_t i=0; i<cluster->GetNDigits(); i++) {
	AliMUONVDigit* digit = digitStore->FindObject(cluster->GetDigitId(i));
	if (!digit) continue;
	
	// pad location
	const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->
	  GetMpSegmentation(digit->DetElemId(),AliMp::GetCathodType(digit->Cathode()));
	AliMpPad pad = seg->PadByIndices(digit->PadX(), digit->PadY());
	
	// calibration parameters
	AliMUONVCalibParam* ped =  fPedStore ? static_cast<AliMUONVCalibParam*>(fPedStore->FindObject(digit->DetElemId(), digit->ManuId())) : 0x0;
	AliMUONVCalibParam* gain = fGainStore ? static_cast<AliMUONVCalibParam*>(fGainStore->FindObject(digit->DetElemId(), digit->ManuId())) : 0x0;
	Int_t manuChannel = digit->ManuChannel();
	Int_t planeType = 0;
	if ( digit->ManuId() & AliMpConstants::ManuMask(AliMp::kNonBendingPlane)) {
	  planeType = 1;
	}
	
	// fill pad info
	padInfo.SetPadId(digit->GetUniqueID());
	padInfo.SetPadPlaneType(planeType);
	padInfo.SetPadXY(pad.GetPositionX(), pad.GetPositionY());
	padInfo.SetPadDimXY(pad.GetDimensionX(), pad.GetDimensionY());
	padInfo.SetPadCharge((Double_t)digit->Charge());
	padInfo.SetPadADC(digit->ADC());
	padInfo.SetSaturated(digit->IsSaturated());
	padInfo.SetCalibrated(digit->IsCalibrated());
	if (ped) {
	  padInfo.SetPedestal(ped->ValueAsFloatFast(manuChannel,0), ped->ValueAsFloatFast(manuChannel,1));
	} else {
	  padInfo.SetPedestal(-250.,-5.);
	}
	if (gain) {
	  padInfo.SetGain(gain->ValueAsFloatFast(manuChannel,0), gain->ValueAsFloatFast(manuChannel,1),
			  gain->ValueAsIntFast(manuChannel,2), gain->ValueAsIntFast(manuChannel,3));
	} else {
	  padInfo.SetGain(-1.,-0.1,-4095,-1);
	} 	
	fClusterInfo->AddPad(padInfo);
      }
      	  
      // fill cluster info tree
      fClusterInfoTree->Fill();
      trackParam = static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->After(trackParam));
    }
  }
  // free memory
  delete newTrackStore;
    
  
  // Post final data. Write histo list to a file with option "RECREATE"
  PostData(0,fClusterInfoTree);
}      

//________________________________________________________________________
void AliMUONReAlignTask::Terminate(const Option_t*)
{
  /// Called once per task on the client machine at the end of the analysis.

}

//-----------------------------------------------------------------------
void AliMUONReAlignTask::Prepare(const char* geoFilename, const char* defaultOCDB, const char* misAlignOCDB)
{
  /// Set the geometry, the magnetic field, the mapping and the reconstruction parameters
  
  // Import TGeo geometry (needed by AliMUONTrackExtrap::ExtrapToVertex)
  if (!gGeoManager) {
    AliGeomManager::LoadGeometry(geoFilename);
    if (!gGeoManager) {
      Error("AliMUONReAlignTask", "getting geometry from file %s failed", "generated/galice.root");
      return;
    }
  }
    
  // Load mapping
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(defaultOCDB);
  man->SetSpecificStorage("MUON/Align/Data",misAlignOCDB);
  man->Print();
  man->SetRun(0);
  if ( ! AliMpCDB::LoadDDLStore() ) {
    Error("MUONRefit","Could not access mapping from OCDB !");
    exit(-1);
  }

  // set mag field
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    printf("Loading field map...\n");
    AliGRPManager *grpMan = new AliGRPManager();
    grpMan->ReadGRPEntry();
    grpMan->SetMagField();
    delete grpMan;
  }
  // set the magnetic field for track extrapolations
  AliMUONTrackExtrap::SetField();
  
  // Load initial reconstruction parameters from OCDB
  // reconstruction parameters
  fRecoParam = AliMUONRecoParam::GetCosmicParam();
  
  // digit selection
  fRecoParam->SetPadGoodnessMask(0x400BE80);
//   TString caliboption = caliboption1;
//   if ( calib == 2 ) caliboption = caliboption2;
//   fRecoParam->SetCalibrationMode("NOGAIN"caliboption.Data());
  fRecoParam->SetCalibrationMode("NOGAIN");
  
  // chamber resolution (incuding misalignment)
  for (Int_t iCh=0; iCh<10; iCh++) {
    fRecoParam->SetDefaultNonBendingReso(iCh,0.4);
    fRecoParam->SetDefaultBendingReso(iCh,0.4);
  }
  fRecoParam->SetMaxNonBendingDistanceToTrack(10.);
  fRecoParam->SetMaxBendingDistanceToTrack(10.);
  
  // cut on (non)bending slopes
  //fRecoParam->SetMaxNonBendingSlope(0.6);
  //fRecoParam->SetMaxBendingSlope(0.6);
  
  // tracking algorithm
  //  fRecoParam->MakeMoreTrackCandidates(kTRUE);
  fRecoParam->RequestStation(0, kFALSE);
  fRecoParam->RequestStation(2, kFALSE);
//   fRecoParam->RequestStation(3, kFALSE);
//   fRecoParam->RequestStation(4, kFALSE);
  fRecoParam->SetSigmaCutForTracking(7.);
  fRecoParam->ImproveTracks(kTRUE, 7.);
  Info("MUONRefit", "\n initial recontruction parameters:");
  fRecoParam->Print("FULL");
  
  AliMUONESDInterface::ResetTracker(fRecoParam);  
}

//-----------------------------------------------------------------------
UInt_t AliMUONReAlignTask::BuildClusterMap(AliMUONTrack &track)
{
  /// Build the map of clusters in tracking chambers
  
  UInt_t muonClusterMap = 0;
  
  AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->First());
  while (trackParam) {
    
    muonClusterMap |= BIT(trackParam->GetClusterPtr()->GetChamberId());
    
    trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->After(trackParam));
  }
  
  return muonClusterMap;
  
}
