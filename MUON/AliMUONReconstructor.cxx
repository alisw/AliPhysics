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
/* $Id$ */

//-----------------------------------------------------------------------------
/// \class AliMUONReconstructor
///
/// Implementation of AliReconstructor for MUON subsystem.
///
/// The clustering mode and the associated parameters can be changed by using
/// AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLow(High)FluxParam();
/// muonRecoParam->Set...(); // see methods in AliMUONRecoParam.h for details
/// AliRecoParam::Instance()->RegisterRecoParam(muonRecoParam);
///
/// Valid modes are :
///
/// SIMPLEFIT : use the AliMUONClusterFinderSimpleFit clusterizer
///
/// SIMPLEFITV3 : SIMPLEFIT with preclustering=PRECLUSTERV3
///
/// MLEM : use AliMUONClusterFinderMLEM and AliMUONPreClusterFinder for preclustering (default)
/// MLEMV2 : MLEM with preclustering=PRECLUSTERV2
/// MLEMV3 : MLEM with preclustering=PRECLUSTERV3
///
/// PRECLUSTER : use only AliMUONPreClusterFinder. Only for debug as
/// the produced clusters do not have a position, hence the tracking will not
/// work
/// PRECLUSTERV2 : another version of the preclustering
/// PRECLUSTERV3 : yet another version of the preclustering
///
/// COG : use AliMUONClusterFinderCOG clusterizer. Not really a production
/// option either, as center-of-gravity is generally not a good estimate
/// of the cluster position...
///
/// PEAKCOG : COG cluster finder around local maxima
/// PEAKFIT : fit around local maxima with up to 3 peaks, COG otherwise
///
/// NOCLUSTERING : bypass completely the clustering stage
///
/// ------
///
/// The behavior of the MUON reconstruction can also be changed, besides
/// the usual methods found in AliReconstruction (e.g. to disable tracking)
/// by using AliReconstruction::SetOption("MUON",options)
/// where options should be a space separated string.
///
/// Valid options are :
///
/// SAVEDIGITS : if you want to save in the TreeD the *calibrated* digits
///     that are used for the clustering
///
/// DIGITSTOREV1 : use the V1 implementation of the digitstore 
/// DIGITSTOREV2R : use the V2R implementation of the digitstore 
///
/// NOLOCALRECONSTRUCTION : for debug, to disable local reconstruction (and hence
/// "recover" old behavior)
///
/// TRIGGERDISABLE : disable the treatment of MUON trigger
///
/// USEFASTDECODER : makes the digit maker class use the high performance decoder
///                  AliMUONTrackerDDLDecoder instead of AliMUONPayloadTracker.
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONReconstructor.h"

#include "AliMUONCalibrationData.h"
#include "AliMUONClusterFinderCOG.h"
#include "AliMUONClusterFinderMLEM.h"
#include "AliMUONClusterFinderSimpleFit.h"
#include "AliMUONClusterFinderPeakCOG.h"
#include "AliMUONClusterFinderPeakFit.h"
#include "AliMUONClusterStoreV1.h"
#include "AliMUONClusterStoreV2.h"
#include "AliMUONConstants.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONDigitMaker.h"
#include "AliMUONDigitStoreV1.h"
#include "AliMUONDigitStoreV2R.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONPreClusterFinder.h"
#include "AliMUONPreClusterFinderV2.h"
#include "AliMUONPreClusterFinderV3.h"
#include "AliMUONRecoParam.h"
#include "AliMUONSimpleClusterServer.h"
#include "AliMUONTracker.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONTriggerStoreV1.h"
#include "AliMUONVClusterFinder.h"
#include "AliMUONVClusterServer.h"
#include "AliMUONVTrackStore.h"

#include "AliMpArea.h"
#include "AliMpCDB.h"
#include "AliMpConstants.h"

#include "AliRecoParam.h"
#include "AliRawReader.h"
#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliLog.h"

#include <Riostream.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TTree.h>

/// \cond CLASSIMP
ClassImp(AliMUONReconstructor)
/// \endcond 

AliMUONRecoParam* AliMUONReconstructor::fgRecoParam = 0x0; // reconstruction parameters

//_____________________________________________________________________________
AliMUONReconstructor::AliMUONReconstructor() : 
AliReconstructor(),
fCrateManager(0x0),
fDigitMaker(0x0),
fTransformer(new AliMUONGeometryTransformer()),
fDigitStore(0x0),
fTriggerCircuit(0x0),
fCalibrationData(0x0),
fDigitCalibrator(0x0),
fClusterServer(0x0),
fTriggerStore(0x0),
fTrackStore(0x0)
{
  /// normal ctor

  // Load mapping
  if ( ! AliMpCDB::LoadDDLStore() ) {
    AliFatal("Could not access mapping from OCDB !");
  }
  
  // Load geometry data
  fTransformer->LoadGeometryData();
  
}

//_____________________________________________________________________________
AliMUONReconstructor::~AliMUONReconstructor()
{
  /// dtor
  delete fDigitMaker;
  delete fDigitStore;
  delete fTransformer;
  delete fCrateManager;
  delete fTriggerCircuit;
  delete fCalibrationData;
  delete fDigitCalibrator;
  delete fClusterServer;
  delete fTriggerStore;
  delete fTrackStore;
}

//_____________________________________________________________________________
const AliMUONRecoParam* AliMUONReconstructor::GetRecoParam()
{
  /// get reconstruction parameters
  
  if (!fgRecoParam) {
    
    // get reconstruction parameters from AliRecoParam if any
    TObjArray *recoParams = AliRecoParam::Instance()->GetRecoParam("MUON");
    
    if (recoParams) {
      
      fgRecoParam = (AliMUONRecoParam*) recoParams->Last();
      
    } else {
      
      // initialize reconstruction parameters if not already done
      cout<<"W-AliMUONReconstructor::GetRecoParam: Reconstruction parameters not initialized - Use default one"<<endl;
      fgRecoParam = AliMUONRecoParam::GetLowFluxParam();
      AliRecoParam::Instance()->RegisterRecoParam(fgRecoParam);
      
    }
    
  }
  
  return fgRecoParam;
}

//_____________________________________________________________________________
void
AliMUONReconstructor::Calibrate(AliMUONVDigitStore& digitStore) const
{
  /// Calibrate the digitStore
  if (!fDigitCalibrator)
  {
    CreateCalibrator();
  }
  AliCodeTimerAuto(Form("%s::Calibrate(AliMUONVDigitStore*)",fDigitCalibrator->ClassName()))
  fDigitCalibrator->Calibrate(digitStore);  
}

//_____________________________________________________________________________
void
AliMUONReconstructor::ConvertDigits(AliRawReader* rawReader, 
                                    AliMUONVDigitStore* digitStore,
                                    AliMUONVTriggerStore* triggerStore) const
{
  /// Convert raw data into digit and trigger stores
  CreateDigitMaker();
  
  AliCodeTimerStart(Form("%s::Raw2Digits(AliRawReader*,AliMUONVDigitStore*,AliMUONVTriggerStore*)",
                    fDigitMaker->ClassName()))
  fDigitMaker->Raw2Digits(rawReader,digitStore,triggerStore);
  AliCodeTimerStop(Form("%s::Raw2Digits(AliRawReader*,AliMUONVDigitStore*,AliMUONVTriggerStore*)",
                         fDigitMaker->ClassName()))
  Calibrate(*digitStore);
}

//_____________________________________________________________________________
void 
AliMUONReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
   /// convert raw data into a digit tree
  AliCodeTimerAuto("")

  Bool_t alone = ( TriggerStore() == 0 );
  
  Bool_t ok = DigitStore()->Connect(*digitsTree,alone);
  if ( TriggerStore() ) 
  {
    ok = ok && TriggerStore()->Connect(*digitsTree,kFALSE);
  }
  
  if (!ok)
  {
    AliError("Could not make branches on TreeD");
  }
  else
  {
    ConvertDigits(rawReader,DigitStore(),TriggerStore());
    AliCodeTimerStart("Fill digits")
    digitsTree->Fill();
    AliCodeTimerStop("Fill digits")
    DigitStore()->Clear();
  }
}

//_____________________________________________________________________________
AliMUONTriggerCrateStore*
AliMUONReconstructor::CrateManager() const
{
  /// Return (and create if necessary) the trigger crate store
  if (fCrateManager) return fCrateManager;
  fCrateManager = new AliMUONTriggerCrateStore;
  fCrateManager->ReadFromFile();
  return fCrateManager;
}

//_____________________________________________________________________________
void
AliMUONReconstructor::CreateDigitMaker() const
{
  /// Create (and create if necessary) the digit maker
  if (fDigitMaker) return;

  AliCodeTimerAuto("")

  TString option = GetOption();
  Bool_t enableErrorLogging = kTRUE;
  Bool_t useFastDecoder = kFALSE;
  if (option.Contains("USEFASTDECODER"))
  {
    useFastDecoder = kTRUE;
  }
  fDigitMaker = new AliMUONDigitMaker(enableErrorLogging, useFastDecoder);
  option.ToUpper();
  if ( option.Contains("SAVEDIGITS" ))
    {
      fDigitMaker->SetMakeTriggerDigits(kTRUE);
    }
}

//_____________________________________________________________________________
void 
AliMUONReconstructor::CreateTriggerCircuit() const
{
  /// Return (and create if necessary) the trigger circuit object
  if (fTriggerCircuit) return;

  AliCodeTimerAuto("")

  fTriggerCircuit = new AliMUONTriggerCircuit(fTransformer);

}

//_____________________________________________________________________________
AliTracker* 
AliMUONReconstructor::CreateTracker() const
{
  /// Create the MUONTracker object
  
  CreateTriggerCircuit();
  CreateDigitMaker();
  CreateClusterServer();

  AliMUONTracker* tracker(0x0);
  
  if ( ! AliMUONReconstructor::GetRecoParam()->CombineClusterTrackReco() )
  {
    tracker = new AliMUONTracker(0x0,
                                 *DigitStore(),
                                 fDigitMaker,
                                 fTransformer,
                                 fTriggerCircuit);
  }
  else
  {
    tracker = new AliMUONTracker(fClusterServer,
                                 *DigitStore(),
                                 fDigitMaker,
                                 fTransformer,
                                 fTriggerCircuit);
  }
  
  
  return tracker;
}

//_____________________________________________________________________________
AliMUONVClusterFinder*
AliMUONReconstructor::CreateClusterFinder(const char* clusterFinderType)
{
  /// Create a given cluster finder instance
  
  AliCodeTimerAutoGeneral("")

  AliMUONVClusterFinder* clusterFinder(0x0);
  
  TString opt(clusterFinderType);
  opt.ToUpper();
  
  if ( strstr(opt,"PRECLUSTERV2") )
  {
    clusterFinder = new AliMUONPreClusterFinderV2;
  }    
  else if ( strstr(opt,"PRECLUSTERV3") )
  {
    clusterFinder = new AliMUONPreClusterFinderV3;
  }  
  else if ( strstr(opt,"PRECLUSTER") )
  {
    clusterFinder = new AliMUONPreClusterFinder;
  }  
  else if ( strstr(opt,"PEAKCOG") )
  {
    clusterFinder = new AliMUONClusterFinderPeakCOG(kFALSE,new AliMUONPreClusterFinder);
  }
  else if ( strstr(opt,"PEAKFIT") )
  {
    clusterFinder = new AliMUONClusterFinderPeakFit(kFALSE,new AliMUONPreClusterFinder);
  }
  else if ( strstr(opt,"COG") )
  {
    clusterFinder = new AliMUONClusterFinderCOG(new AliMUONPreClusterFinder);
  }  
  else if ( strstr(opt,"SIMPLEFITV3") )
  {
    clusterFinder = new AliMUONClusterFinderSimpleFit(new AliMUONClusterFinderCOG(new AliMUONPreClusterFinderV3));
  }
  else if ( strstr(opt,"SIMPLEFIT") )
  {
    clusterFinder = new AliMUONClusterFinderSimpleFit(new AliMUONClusterFinderCOG(new AliMUONPreClusterFinder));
  }
  else if ( strstr(opt,"MLEM:DRAW") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kTRUE,new AliMUONPreClusterFinder);
  }
  else if ( strstr(opt,"MLEMV3") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kFALSE,new AliMUONPreClusterFinderV3);
  } 
  else if ( strstr(opt,"MLEMV2") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kFALSE,new AliMUONPreClusterFinderV2);
  } 
  else if ( strstr(opt,"MLEM") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kFALSE,new AliMUONPreClusterFinder);
  } 
  else
  {
    AliErrorClass(Form("clustering mode \"%s\" does not exist",opt.Data()));
    return 0x0;
  }
  
  return clusterFinder;
}

//_____________________________________________________________________________
void
AliMUONReconstructor::CreateClusterServer() const
{
  /// Create cluster server
  
  if ( fClusterServer ) return;
  
  AliCodeTimerAuto("");
    
  AliMUONVClusterFinder* clusterFinder = CreateClusterFinder(GetRecoParam()->GetClusteringMode());
  
  if ( !clusterFinder ) return;
  
  AliInfo(Form("Will use %s for clusterizing",clusterFinder->ClassName()));
  
  fClusterServer = new AliMUONSimpleClusterServer(clusterFinder,*fTransformer);
}

//_____________________________________________________________________________
void
AliMUONReconstructor::CreateCalibrator() const
{
  /// Create the calibrator
  
  AliCodeTimerAuto("")
  
  Int_t runNumber = AliCDBManager::Instance()->GetRun();

  AliInfo("Calibration will occur.");
  
  fCalibrationData = new AliMUONCalibrationData(runNumber);
  if ( !fCalibrationData->IsValid() )
  {
    AliError("Could not retrieve calibrations !");
    delete fCalibrationData;
    fCalibrationData = 0x0;
    return;
  }    
  
  // Check that we get all the calibrations we'll need
  if ( !fCalibrationData->Pedestals() ||
       !fCalibrationData->Gains() ||
       !fCalibrationData->HV() )
  {
    AliFatal("Could not access all required calibration data");
  }
  
  TString opt(GetOption());
  opt.ToUpper();
  
  if ( strstr(opt,"NOSTATUSMAP") )
  {
    AliWarning("NOSTATUSMAP is obsolete");
  }

  TString calibMode = GetRecoParam()->GetCalibrationMode();

  fDigitCalibrator = new AliMUONDigitCalibrator(*fCalibrationData,calibMode.Data());
}

//_____________________________________________________________________________
AliMUONVDigitStore*
AliMUONReconstructor::DigitStore() const
{
  /// Return (and create if necessary) the digit container
  if (!fDigitStore) 
  {
    TString sopt(GetOption());
    sopt.ToUpper();
    
    AliInfo(Form("Options=%s",sopt.Data()));
    
    if ( sopt.Contains("DIGITSTOREV1") )
    {
      fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV1");
    }
    else if ( sopt.Contains("DIGITSTOREV2R") ) 
    {
      fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV2R");
    }
    else if ( sopt.Contains("DIGITSTOREV2S") ) 
    {
      fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV2S");
    }
    
    if (!fDigitStore) fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV2R");
    
    AliInfo(Form("Will use %s to store digits during reconstruction",fDigitStore->ClassName()));
  }
  return fDigitStore;
}

//_____________________________________________________________________________
void
AliMUONReconstructor::FillTreeR(AliMUONVTriggerStore* triggerStore,
                                TTree& clustersTree) const
{
  /// Write the trigger and cluster information into TreeR
  
  AliCodeTimerAuto("")

  AliDebug(1,"");
  
  Bool_t ok(kFALSE);
  Bool_t alone(kTRUE); // is trigger the only info in TreeR ?
  
  if ( ! AliMUONReconstructor::GetRecoParam()->CombineClusterTrackReco() )
  {
    alone = kFALSE; // we'll get both tracker and trigger information in TreeR
  }
  
  if ( triggerStore ) 
  {
    ok = triggerStore->Connect(clustersTree,alone);
    if (!ok)
    {
      AliError("Could not create triggerStore branches in TreeR");
    }
  }

  AliMUONVClusterStore* clusterStore(0x0);
  
  if ( !alone )
  {
    clusterStore = new AliMUONClusterStoreV2;
    
    CreateClusterServer();
    
    TIter next(DigitStore()->CreateIterator());
    fClusterServer->UseDigits(next);

    AliMpArea area;
    
    AliDebug(1,Form("Doing full clusterization in local reconstruction using %s ",fClusterServer->ClassName()));
    
    for ( Int_t i = 0; i < AliMpConstants::NofTrackingChambers(); ++i ) 
    {
      if (AliMUONReconstructor::GetRecoParam()->UseChamber(i))
      {
        if ( i >= 6 && AliMUONReconstructor::GetRecoParam()->BypassSt45() ) continue;
        
        fClusterServer->Clusterize(i,*clusterStore,area);
      }
    }
    
    Bool_t cok = clusterStore->Connect(clustersTree,alone);
    
    if (!cok) AliError("Could not connect clusterStore to clusterTree");
    
    AliDebug(1,Form("Number of clusters found = %d",clusterStore->GetSize()));
    
    StdoutToAliDebug(1,clusterStore->Print());
  }
         
  if (ok) // at least one type of branches created successfully
  {
    clustersTree.Fill();
  }
  
  delete clusterStore;
}

//_____________________________________________________________________________
Bool_t 
AliMUONReconstructor::HasDigitConversion() const
{
  /// We *do* have digit conversion, but we might advertise it only 
  /// if we want to save the digits.
  
  TString opt(GetOption());
  opt.ToUpper();
  if ( opt.Contains("SAVEDIGITS" ) && !opt.Contains("NOLOCALRECONSTRUCTION") )
  {
    return kTRUE;
  }
  else
  {
    return kFALSE;
  }
}

//_____________________________________________________________________________
void 
AliMUONReconstructor::Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const
{
  /// This method is called by AliReconstruction if HasLocalReconstruction()==kTRUE AND
  /// HasDigitConversion()==kFALSE
  
  if ( !clustersTree ) 
  {
    AliError("clustersTree is 0x0 !");
    return;
  }
  
  ConvertDigits(rawReader,DigitStore(),TriggerStore());

  FillTreeR(TriggerStore(),*clustersTree);
}

//_____________________________________________________________________________
void 
AliMUONReconstructor::Reconstruct(TTree* digitsTree, TTree* clustersTree) const
{
  /// This method is called by AliReconstruction if HasLocalReconstruction()==kTRUE
  /// AND HasDigitConversion()==kTRUE
  
  AliCodeTimerAuto("")
  
  AliDebug(1,"");
  
  if (!digitsTree || !clustersTree) 
  {
    AliError(Form("Tree is null : digitsTree=%p clustersTree=%p",
                  digitsTree,clustersTree));
    return;
  }

  if (!fDigitStore)
  {
    fDigitStore = AliMUONVDigitStore::Create(*digitsTree);
    if (!fDigitStore)
    {
      AliError(Form("Could not get DigitStore from %s",digitsTree->GetName()));
    }
    else
    {
      AliInfo(Form("Created %s from %s",fDigitStore->ClassName(),digitsTree->GetName()));
    }
  }
  if (!fTriggerStore)
  {
    fTriggerStore = AliMUONVTriggerStore::Create(*digitsTree);
    if (!fTriggerStore)
    {
      AliError(Form("Could not get TriggerStore from %s",digitsTree->GetName()));
    }
    else
    {
      AliInfo(Form("Created %s from %s",fTriggerStore->ClassName(),digitsTree->GetName()));
    }
  }
  
  if (!fTriggerStore && !fDigitStore)
  {
    AliError("No store at all. Nothing to do.");
    return;
  }
  
  // insure we start with empty stores
  if ( fDigitStore ) 
  {
    fDigitStore->Clear(); 
    Bool_t alone = ( fTriggerStore ? kFALSE : kTRUE );
    Bool_t ok = fDigitStore->Connect(*digitsTree,alone);
    if (!ok)
    {
      AliError("Could not connect digitStore to digitsTree");
      return;
    }
  }
  if ( fTriggerStore ) 
  {
    fTriggerStore->Clear();
    Bool_t alone = ( fDigitStore ? kFALSE : kTRUE );
    Bool_t ok = fTriggerStore->Connect(*digitsTree,alone);
    if (!ok)
    {
      AliError("Could not connect triggerStore to digitsTree");
      return;
    }
  }
  
  digitsTree->GetEvent(0);
  
  if ( fDigitStore ) 
  {
    // Insure we got calibrated digits (if we reconstruct from pure simulated,
    // i.e. w/o going through raw data, this will be the case)
    TIter next(fDigitStore->CreateIterator());
    AliMUONVDigit* digit = static_cast<AliMUONVDigit*>(next());
    if (digit && !digit->IsCalibrated())
    {
      Calibrate(*fDigitStore);
    }
  }
    
  FillTreeR(fTriggerStore,*clustersTree);
}

//_____________________________________________________________________________
AliMUONVTriggerStore*
AliMUONReconstructor::TriggerStore() const
{
  /// Return (and create if necessary and allowed) the trigger container
  TString sopt(GetOption());
  sopt.ToUpper();
  
  if (sopt.Contains("TRIGGERDISABLE"))
  {
    delete fTriggerStore;
    fTriggerStore = 0x0;
  }
  else
  {
    if (!fTriggerStore)
    {
      fTriggerStore = new AliMUONTriggerStoreV1;
    }
  }
  return fTriggerStore;
}
