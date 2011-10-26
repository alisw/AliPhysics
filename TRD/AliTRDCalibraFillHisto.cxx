//**************************************************************************
// * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// *                                                                        *
// * Author: The ALICE Off-line Project.                                    *
// * Contributors are mentioned in the code where appropriate.              *
//   *                                                                        *
//   * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//    * without fee, provided that the above copyright notice appears in all   *
//    * copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************/

/* $Id$ */

/////////////////////////////////////////////////////////////////////////////////
//                                                                             
// AliTRDCalibraFillHisto                                                               
//                                                                             
// This class is for the TRD calibration of the relative gain factor, the drift velocity,
// the time 0 and the pad response function. It fills histos or vectors.        
// It can be used for the calibration per chamber but also per group of pads and eventually per pad.
// The user has to choose with the functions SetNz and SetNrphi the precision of the calibration (see AliTRDCalibraMode). 
// 2D Histograms (Histo2d) or vectors (Vector2d), then converted in Trees, will be filled
// from RAW DATA in a run or from reconstructed TRD tracks during the offline tracking 
// in the function "FollowBackProlongation" (AliTRDtracker)
// Per default the functions to fill are off.                                   
//                        
// Authors:
//   R. Bailhache (R.Bailhache@gsi.de, rbailhache@ikf.uni-frankfurt.de)
//   J. Book (jbook@ikf.uni-frankfurt.de)
//                            
//////////////////////////////////////////////////////////////////////////////////////

#include <TProfile2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TH2.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TTreeStream.h>
#include <TVectorD.h>
#include <TLinearFitter.h>

#include "AliLog.h"

#include "AliTRDCalibraFillHisto.h"
#include "AliTRDCalibraMode.h"
#include "AliTRDCalibraVector.h"
#include "AliTRDCalibraVdriftLinearFit.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "AliTRDpadPlane.h"
#include "AliTRDcluster.h"
#include "AliTRDtrackV1.h"
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliTRDgeometry.h"
#include "./Cal/AliTRDCalROC.h"
#include "./Cal/AliTRDCalPad.h"
#include "./Cal/AliTRDCalDet.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDdigitsParam.h"
#include "AliTRDSignalIndex.h"
#include "AliTRDarrayADC.h"

#include "AliTRDrawStream.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"

#ifdef ALI_DATE
#include "event.h"
#endif


ClassImp(AliTRDCalibraFillHisto)

AliTRDCalibraFillHisto* AliTRDCalibraFillHisto::fgInstance = 0;
Bool_t AliTRDCalibraFillHisto::fgTerminated = kFALSE;

//_____________singleton implementation_________________________________________________
AliTRDCalibraFillHisto *AliTRDCalibraFillHisto::Instance()
{
  //
  // Singleton implementation
  //

  if (fgTerminated != kFALSE) {
    return 0;
  }

  if (fgInstance == 0) {
    fgInstance = new AliTRDCalibraFillHisto();
  }

  return fgInstance;

}

//______________________________________________________________________________________
void AliTRDCalibraFillHisto::Terminate()
{
  //
  // Singleton implementation
  // Deletes the instance of this class
  //

  fgTerminated = kTRUE;

  if (fgInstance != 0) {
    delete fgInstance;
    fgInstance = 0;
  }

}

//______________________________________________________________________________________
AliTRDCalibraFillHisto::AliTRDCalibraFillHisto()
  :TObject()
  ,fGeo(0)
  ,fCalibDB(0)
  ,fIsHLT(kFALSE)
  ,fCH2dOn(kFALSE)
  ,fPH2dOn(kFALSE)
  ,fPRF2dOn(kFALSE)
  ,fHisto2d(kFALSE)
  ,fVector2d(kFALSE)
  ,fLinearFitterOn(kFALSE)
  ,fLinearFitterDebugOn(kFALSE)
  ,fRelativeScale(0)
  ,fThresholdClusterPRF2(15.0)
  ,fLimitChargeIntegration(kFALSE)
  ,fFillWithZero(kFALSE)
  ,fNormalizeNbOfCluster(kFALSE)
  ,fMaxCluster(0)
  ,fNbMaxCluster(0)
  ,fFirstRunGain(0)
  ,fVersionGainUsed(0)
  ,fSubVersionGainUsed(0)
  ,fFirstRunGainLocal(0)
  ,fVersionGainLocalUsed(0)
  ,fSubVersionGainLocalUsed(0)
  ,fFirstRunVdrift(0)
  ,fVersionVdriftUsed(0) 
  ,fSubVersionVdriftUsed(0)
  ,fFirstRunExB(0)
  ,fVersionExBUsed(0) 
  ,fSubVersionExBUsed(0)
  ,fCalibraMode(new AliTRDCalibraMode())
  ,fDebugStreamer(0)
  ,fDebugLevel(0)
  ,fDetectorPreviousTrack(-1)
  ,fMCMPrevious(-1)
  ,fROBPrevious(-1)
  ,fNumberClusters(1)
  ,fNumberClustersf(30)
  ,fNumberClustersProcent(0.5)
  ,fThresholdClustersDAQ(120.0)
  ,fNumberRowDAQ(2)
  ,fNumberColDAQ(4)
  ,fProcent(6.0)
  ,fDifference(17)
  ,fNumberTrack(0)
  ,fTimeMax(0)
  ,fSf(10.0)
  ,fNumberBinCharge(50)
  ,fNumberBinPRF(10)
  ,fNgroupprf(3)
  ,fAmpTotal(0x0)
  ,fPHPlace(0x0)
  ,fPHValue(0x0)
  ,fGoodTracklet(kTRUE)
  ,fLinearFitterTracklet(0x0)
  ,fEntriesCH(0x0)
  ,fEntriesLinearFitter(0x0)
  ,fCalibraVector(0x0)
  ,fPH2d(0x0)
  ,fPRF2d(0x0)
  ,fCH2d(0x0)
  ,fLinearFitterArray(540)
  ,fLinearVdriftFit(0x0)
  ,fCalDetGain(0x0)
  ,fCalROCGain(0x0)
{
  //
  // Default constructor
  //

  //
  // Init some default values
  //

  fNumberUsedCh[0]       = 0;
  fNumberUsedCh[1]       = 0;
  fNumberUsedPh[0]       = 0;
  fNumberUsedPh[1]       = 0;
  
  fGeo = new AliTRDgeometry();
  fCalibDB = AliTRDcalibDB::Instance();
}

//______________________________________________________________________________________
AliTRDCalibraFillHisto::AliTRDCalibraFillHisto(const AliTRDCalibraFillHisto &c)
  :TObject(c)
  ,fGeo(0)
  ,fCalibDB(0)
  ,fIsHLT(c.fIsHLT)
  ,fCH2dOn(c.fCH2dOn)
  ,fPH2dOn(c.fPH2dOn)
  ,fPRF2dOn(c.fPRF2dOn)
  ,fHisto2d(c.fHisto2d)
  ,fVector2d(c.fVector2d)
  ,fLinearFitterOn(c.fLinearFitterOn)
  ,fLinearFitterDebugOn(c.fLinearFitterDebugOn)
  ,fRelativeScale(c.fRelativeScale)
  ,fThresholdClusterPRF2(c.fThresholdClusterPRF2)
  ,fLimitChargeIntegration(c.fLimitChargeIntegration)
  ,fFillWithZero(c.fFillWithZero)
  ,fNormalizeNbOfCluster(c.fNormalizeNbOfCluster)
  ,fMaxCluster(c.fMaxCluster)
  ,fNbMaxCluster(c.fNbMaxCluster)
  ,fFirstRunGain(c.fFirstRunGain)
  ,fVersionGainUsed(c.fVersionGainUsed)
  ,fSubVersionGainUsed(c.fSubVersionGainUsed)
  ,fFirstRunGainLocal(c.fFirstRunGainLocal)
  ,fVersionGainLocalUsed(c.fVersionGainLocalUsed)
  ,fSubVersionGainLocalUsed(c.fSubVersionGainLocalUsed)
  ,fFirstRunVdrift(c.fFirstRunVdrift)
  ,fVersionVdriftUsed(c.fVersionVdriftUsed) 
  ,fSubVersionVdriftUsed(c.fSubVersionVdriftUsed)
  ,fFirstRunExB(c.fFirstRunExB)
  ,fVersionExBUsed(c.fVersionExBUsed) 
  ,fSubVersionExBUsed(c.fSubVersionExBUsed)
  ,fCalibraMode(0x0)
  ,fDebugStreamer(0)
  ,fDebugLevel(c.fDebugLevel)
  ,fDetectorPreviousTrack(c.fDetectorPreviousTrack)
  ,fMCMPrevious(c.fMCMPrevious)
  ,fROBPrevious(c.fROBPrevious)
  ,fNumberClusters(c.fNumberClusters)
  ,fNumberClustersf(c.fNumberClustersf)
  ,fNumberClustersProcent(c.fNumberClustersProcent)
  ,fThresholdClustersDAQ(c.fThresholdClustersDAQ)
  ,fNumberRowDAQ(c.fNumberRowDAQ)
  ,fNumberColDAQ(c.fNumberColDAQ)
  ,fProcent(c.fProcent)
  ,fDifference(c.fDifference)
  ,fNumberTrack(c.fNumberTrack)
  ,fTimeMax(c.fTimeMax)
  ,fSf(c.fSf)
  ,fNumberBinCharge(c.fNumberBinCharge)
  ,fNumberBinPRF(c.fNumberBinPRF)
  ,fNgroupprf(c.fNgroupprf)
  ,fAmpTotal(0x0)
  ,fPHPlace(0x0)
  ,fPHValue(0x0)
  ,fGoodTracklet(c.fGoodTracklet)
  ,fLinearFitterTracklet(0x0)
  ,fEntriesCH(0x0)
  ,fEntriesLinearFitter(0x0)
  ,fCalibraVector(0x0)
  ,fPH2d(0x0)
  ,fPRF2d(0x0)
  ,fCH2d(0x0)
  ,fLinearFitterArray(540)
  ,fLinearVdriftFit(0x0)
  ,fCalDetGain(0x0)
  ,fCalROCGain(0x0)
{
  //
  // Copy constructor
  //
  if(c.fCalibraMode)   fCalibraMode = new AliTRDCalibraMode(*c.fCalibraMode);
  if(c.fCalibraVector) fCalibraVector = new AliTRDCalibraVector(*c.fCalibraVector);
  if(c.fPH2d) {
    fPH2d = new TProfile2D(*c.fPH2d);
    fPH2d->SetDirectory(0);
  }
  if(c.fPRF2d) {
    fPRF2d = new TProfile2D(*c.fPRF2d);
    fPRF2d->SetDirectory(0);
  }
  if(c.fCH2d) {
    fCH2d = new TH2I(*c.fCH2d);
    fCH2d->SetDirectory(0);
  }
  if(c.fLinearVdriftFit){
    fLinearVdriftFit = new AliTRDCalibraVdriftLinearFit(*c.fLinearVdriftFit);
  }

  if(c.fCalDetGain)  fCalDetGain   = new AliTRDCalDet(*c.fCalDetGain);
  if(c.fCalROCGain)  fCalROCGain   = new AliTRDCalROC(*c.fCalROCGain);

  if (fGeo) {
    delete fGeo;
  }
  fGeo = new AliTRDgeometry();
  fCalibDB = AliTRDcalibDB::Instance();

  fNumberUsedCh[0]       = 0;
  fNumberUsedCh[1]       = 0;
  fNumberUsedPh[0]       = 0;
  fNumberUsedPh[1]       = 0;

}

//____________________________________________________________________________________
AliTRDCalibraFillHisto::~AliTRDCalibraFillHisto()
{
  //
  // AliTRDCalibraFillHisto destructor
  //

  ClearHistos();
  if ( fDebugStreamer ) delete fDebugStreamer;

  if ( fCalDetGain )  delete fCalDetGain;
  if ( fCalROCGain )  delete fCalROCGain;

  if( fLinearFitterTracklet ) { delete fLinearFitterTracklet; }
  
  delete [] fPHPlace;
  delete [] fPHValue;
  delete [] fEntriesCH;
  delete [] fEntriesLinearFitter;
  delete [] fAmpTotal;
  
  for(Int_t idet=0; idet<AliTRDgeometry::kNdet; idet++){ 
    TLinearFitter *f = (TLinearFitter*)fLinearFitterArray.At(idet);
    if(f) { delete f;}
  }
  if(fLinearVdriftFit) delete fLinearVdriftFit;
  if (fGeo) {
    delete fGeo;
  }
  
}
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::Destroy() 
{
  //
  // Delete instance 
  //

  if (fgInstance) {
    delete fgInstance;
    fgInstance = 0x0;
  }
}
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::DestroyDebugStreamer() 
{
  //
  // Delete DebugStreamer
  //

  if ( fDebugStreamer ) delete fDebugStreamer;
  fDebugStreamer = 0x0;
 
}
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::ClearHistos() 
{
  //
  // Delete the histos
  //

  if (fPH2d) {
    delete fPH2d;
    fPH2d  = 0x0;
  }
  if (fCH2d) {
    delete fCH2d;
    fCH2d  = 0x0;
  }
  if (fPRF2d) {
    delete fPRF2d;
    fPRF2d = 0x0;
  }
  
}
//////////////////////////////////////////////////////////////////////////////////
// calibration with AliTRDtrackV1: Init, Update 
//////////////////////////////////////////////////////////////////////////////////
//____________Functions for initialising the AliTRDCalibraFillHisto in the code_________
Bool_t AliTRDCalibraFillHisto::Init2Dhistos(Int_t nboftimebin)
{
  //
  // Init the histograms and stuff to be filled 
  //

  // DB Setting
  // Get cal
  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam");
    return kFALSE;
  }

  // Some parameters
  if(nboftimebin > 0) fTimeMax = nboftimebin;
  else fTimeMax = cal->GetNumberOfTimeBinsDCS();
  if(fTimeMax <= 0) fTimeMax = 30;
  printf("////////////////////////////////////////////\n");
  printf("Number of time bins in calibration component %d\n",fTimeMax);
  printf("////////////////////////////////////////////\n");
  fSf                 = parCom->GetSamplingFrequency();
  if(!fNormalizeNbOfCluster) fRelativeScale = 20.0;
  else fRelativeScale = 1.18;
  fNumberClustersf    = fTimeMax;
  fNumberClusters     = (Int_t)(fNumberClustersProcent*fTimeMax);
 
  // Init linear fitter
  if(!fLinearFitterTracklet) {
    fLinearFitterTracklet = new TLinearFitter(2,"pol1");
    fLinearFitterTracklet->StoreData(kTRUE);
  }

  // Calcul Xbins Chambd0, Chamb2
  Int_t ntotal0 = CalculateTotalNumberOfBins(0);
  Int_t ntotal1 = CalculateTotalNumberOfBins(1);
  Int_t ntotal2 = CalculateTotalNumberOfBins(2);

  // If vector method On initialised all the stuff
  if(fVector2d){   
    fCalibraVector = new AliTRDCalibraVector();
    fCalibraVector->SetNumberBinCharge(fNumberBinCharge);
    fCalibraVector->SetTimeMax(fTimeMax);
    if(fNgroupprf != 0) {
      fCalibraVector->SetNumberBinPRF(2*fNgroupprf*fNumberBinPRF);
      fCalibraVector->SetPRFRange((Float_t)(3.0*fNgroupprf));
    }
    else {
      fCalibraVector->SetNumberBinPRF(fNumberBinPRF);
      fCalibraVector->SetPRFRange(1.5);
    }
    for(Int_t k = 0; k < 3; k++){
      fCalibraVector->SetDetCha0(k,fCalibraMode->GetDetChamb0(k));
      fCalibraVector->SetDetCha2(k,fCalibraMode->GetDetChamb2(k));
    }
    fCalibraVector->SetNzNrphi(0,fCalibraMode->GetNz(0),fCalibraMode->GetNrphi(0));
    fCalibraVector->SetNzNrphi(1,fCalibraMode->GetNz(1),fCalibraMode->GetNrphi(1));
    fCalibraVector->SetNzNrphi(2,fCalibraMode->GetNz(2),fCalibraMode->GetNrphi(2));
    fCalibraVector->SetNbGroupPRF(fNgroupprf);
  }
 
  // Create the 2D histos corresponding to the pad groupCalibration mode
  if (fCH2dOn) {

    AliInfo(Form("The pad calibration mode for the relative gain calibration: Nz %d, and Nrphi %d"
                ,fCalibraMode->GetNz(0)
                ,fCalibraMode->GetNrphi(0)));
    
    // Create the 2D histo
    if (fHisto2d) {
      CreateCH2d(ntotal0);
    }
    // Variable
    fAmpTotal = new Float_t[TMath::Max(fCalibraMode->GetDetChamb2(0),fCalibraMode->GetDetChamb0(0))];
    for (Int_t k = 0; k < TMath::Max(fCalibraMode->GetDetChamb2(0),fCalibraMode->GetDetChamb0(0)); k++) {
      fAmpTotal[k] = 0.0;
    } 
    //Statistics
    fEntriesCH = new Int_t[ntotal0];
    for(Int_t k = 0; k < ntotal0; k++){
      fEntriesCH[k] = 0;
    }
    
  }
  if (fPH2dOn) {

    AliInfo(Form("The pad calibration mode for the drift velocity calibration: Nz %d, and Nrphi %d"
                ,fCalibraMode->GetNz(1)
                ,fCalibraMode->GetNrphi(1)));
    
    // Create the 2D histo
    if (fHisto2d) {
      CreatePH2d(ntotal1);
    }
    // Variable
    fPHPlace = new Short_t[fTimeMax];
    for (Int_t k = 0; k < fTimeMax; k++) {
      fPHPlace[k] = -1;
    } 
    fPHValue = new Float_t[fTimeMax];
    for (Int_t k = 0; k < fTimeMax; k++) {
      fPHValue[k] = 0.0;
    }
  }
  if (fLinearFitterOn) {
    if(fLinearFitterDebugOn) {
      fLinearFitterArray.SetName("ArrayLinearFitters");
      fEntriesLinearFitter = new Int_t[540];
      for(Int_t k = 0; k < 540; k++){
	fEntriesLinearFitter[k] = 0;
      }
    }
    fLinearVdriftFit = new AliTRDCalibraVdriftLinearFit();
    TString nameee("Ver");
    nameee += fVersionExBUsed;
    nameee += "Subver";
    nameee += fSubVersionExBUsed;
    nameee += "FirstRun";
    nameee += fFirstRunExB;
    nameee += "Nz";
    fLinearVdriftFit->SetNameCalibUsed(nameee); 
  }

  if (fPRF2dOn) {

    AliInfo(Form("The pad calibration mode for the PRF calibration: Nz %d, and Nrphi %d"
                ,fCalibraMode->GetNz(2)
                ,fCalibraMode->GetNrphi(2)));
    // Create the 2D histo
    if (fHisto2d) {
      CreatePRF2d(ntotal2);
    }
  }

  return kTRUE;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t AliTRDCalibraFillHisto::InitCalDet()
{
  //
  // Init the Gain Cal Det 
  //

  // DB Setting
  // Get cal
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("TRD/Calib/ChamberGainFactor",fFirstRunGain,fVersionGainUsed,fSubVersionGainUsed);
  if(!entry) {
    AliError("No gain det calibration entry found");
    return kFALSE;
  }
  AliTRDCalDet *calDet = (AliTRDCalDet *)entry->GetObject();
  if(!calDet) {
    AliError("No calDet gain found");
    return kFALSE;
  }
   

  if( fCalDetGain ){ 
    fCalDetGain->~AliTRDCalDet();
    new(fCalDetGain) AliTRDCalDet(*(calDet));
  }else fCalDetGain = new AliTRDCalDet(*(calDet));
  
  
  // title CH2d
  TString name("Ver");
  name += fVersionGainUsed;
  name += "Subver";
  name += fSubVersionGainUsed;
  name += "FirstRun";
  name += fFirstRunGain;
  name += "Nz";
  name += fCalibraMode->GetNz(0);
  name += "Nrphi";
  name += fCalibraMode->GetNrphi(0);

  fCH2d->SetTitle(name);  
  
  // title PH2d
  TString namee("Ver");
  namee += fVersionVdriftUsed;
  namee += "Subver";
  namee += fSubVersionVdriftUsed;
  namee += "FirstRun";
  namee += fFirstRunVdrift;
  namee += "Nz";
  namee += fCalibraMode->GetNz(1);
  namee += "Nrphi";
  namee += fCalibraMode->GetNrphi(1);
  
  fPH2d->SetTitle(namee);  

  // title AliTRDCalibraVdriftLinearFit
  TString nameee("Ver");
  nameee += fVersionExBUsed;
  nameee += "Subver";
  nameee += fSubVersionExBUsed;
  nameee += "FirstRun";
  nameee += fFirstRunExB;
  nameee += "Nz";

  
  fLinearVdriftFit->SetNameCalibUsed(nameee);  



  return kTRUE;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t AliTRDCalibraFillHisto::InitCalPad(Int_t detector)
{
  //
  // Init the Gain Cal Pad 
  //

  // DB Setting
  // Get cal
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("TRD/Calib/LocalGainFactor",AliCDBManager::Instance()->GetRun(),fVersionGainLocalUsed,fSubVersionGainLocalUsed);
  if(!entry) {
    AliError("No gain pad calibration entry found");
    return kFALSE;
  }
  AliTRDCalPad *calPad = (AliTRDCalPad *)entry->GetObject();
  if(!calPad) {
    AliError("No calPad gain found");
    return kFALSE;
  }
  AliTRDCalROC *calRoc = (AliTRDCalROC *)calPad->GetCalROC(detector);
  if(!calRoc) {
    AliError("No calRoc gain found");
    return kFALSE;
  }
  
  if( fCalROCGain ){ 
    fCalROCGain->~AliTRDCalROC();
    new(fCalROCGain) AliTRDCalROC(*(calRoc));
  }else fCalROCGain = new AliTRDCalROC(*(calRoc));
  

  
 
  
  return kTRUE;

}
//____________Offline tracking in the AliTRDtracker____________________________
Bool_t AliTRDCalibraFillHisto::UpdateHistogramsV1(const AliTRDtrackV1 *t)
{
  //
  // Use AliTRDtrackV1 for the calibration
  //

  
  const AliTRDseedV1 *tracklet = 0x0;          // tracklet per plane
  AliTRDcluster *cl      = 0x0;                // cluster attached now to the tracklet
  AliTRDcluster *cls     = 0x0;                // shared cluster attached now to the tracklet
  Bool_t         newtr   = kTRUE;              // new track
  
  // Get cal
  //  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  /*
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }
*/
  if (!fCalibDB) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

  
  ///////////////////////////
  // loop over the tracklet
  ///////////////////////////
  for(Int_t itr = 0; itr < 6; itr++){
    
    if(!(tracklet = t->GetTracklet(itr))) continue;
    if(!tracklet->IsOK()) continue;
    fNumberTrack++; 
    ResetfVariablestracklet();

    //////////////////////////////////////////
    // localisation of the tracklet and dqdl
    //////////////////////////////////////////
    Int_t layer    = tracklet->GetPlane();
    Int_t ic = 0;
    while(!(cl = tracklet->GetClusters(ic++))) continue;
    Int_t detector = cl->GetDetector();
    if (detector != fDetectorPreviousTrack) {
      // if not a new track
      if(!newtr){
	// don't use the rest of this track if in the same plane
	if (layer == GetLayer(fDetectorPreviousTrack)) {
	  //printf("bad tracklet, same layer for detector %d\n",detector);
	  break;
	}
      }
      //Localise the detector bin
      LocalisationDetectorXbins(detector);
      // Get calib objects
      if(!fIsHLT) InitCalPad(detector);	
            
      // reset
      fDetectorPreviousTrack = detector;
    }
    newtr = kFALSE;

    ////////////////////////////
    // loop over the clusters
    ////////////////////////////
    Int_t nbclusters = 0;
    for(int jc=0; jc<AliTRDseedV1::kNtb; jc++){
      if(!(cl = tracklet->GetClusters(jc))) continue;
      nbclusters++;
      
      // Store the info bis of the tracklet
      Int_t row = cl->GetPadRow();
      Int_t col = cl->GetPadCol();
      CheckGoodTrackletV1(cl);
      Int_t     group[2] = {0,0};
      if(fCH2dOn)  group[0]  = CalculateCalibrationGroup(0,row,col);
      if(fPH2dOn)  group[1]  = CalculateCalibrationGroup(1,row,col);
      // Add the charge if shared cluster
      cls = tracklet->GetClusters(jc+AliTRDseedV1::kNtb);
      //
      StoreInfoCHPHtrack(cl, tracklet->GetdQdl(jc),group,row,col,cls);
    }
    
    ////////////////////////////////////////
    // Fill the stuffs if a good tracklet
    ////////////////////////////////////////
    if (fGoodTracklet) {

      // drift velocity unables to cut bad tracklets 
      Bool_t  pass = FindP1TrackPHtrackletV1(tracklet, nbclusters);
	
      //printf("pass %d and nbclusters %d\n",pass,nbclusters);

      // Gain calibration
      if (fCH2dOn) {
	FillTheInfoOfTheTrackCH(nbclusters);
      }
	
      // PH calibration
      if (fPH2dOn) {
	FillTheInfoOfTheTrackPH();    
      }
	
      if(pass && fPRF2dOn) HandlePRFtrackletV1(tracklet,nbclusters);
		
    } // if a good tracklet
  }
  
  return kTRUE;
  
}
///////////////////////////////////////////////////////////////////////////////////
// Routine inside the update with AliTRDtrack
///////////////////////////////////////////////////////////////////////////////////
//____________Offine tracking in the AliTRDtracker_____________________________
Bool_t AliTRDCalibraFillHisto::FindP1TrackPHtrackletV1(const AliTRDseedV1 *tracklet, Int_t nbclusters)
{
  //
  // Drift velocity calibration:
  // Fit the clusters with a straight line
  // From the slope find the drift velocity
  //

  ////////////////////////////////////////////////
  //Number of points: if less than 3 return kFALSE
  /////////////////////////////////////////////////
  if(nbclusters <= 2) return kFALSE;

  ////////////
  //Variables
  ////////////
  // results of the linear fit
  Double_t dydt                       = 0.0;                                // dydt tracklet after straight line fit
  Double_t errorpar                   = 0.0;                                // error after straight line fit on dy/dt
  Double_t pointError                 = 0.0;                                // error after straight line fit 
  // pad row problemes: avoid tracklet that cross pad rows, tilting angle in the constant
  Int_t    crossrow                   = 0;                                  // if it crosses a pad row
  Int_t    rowp                       = -1;                                 // if it crosses a pad row
  Float_t  tnt                        = tracklet->GetTilt();                // tan tiltingangle
  fLinearFitterTracklet->ClearPoints();  
 
  
  ///////////////////////////////////////////
  // Take the parameters of the track
  //////////////////////////////////////////
  // take now the snp, tnp and tgl from the track
  Double_t snp = tracklet->GetSnp();             // sin dy/dx at the end of the chamber
  Double_t tnp = 0.0;                            // dy/dx at the end of the chamber 
  if( TMath::Abs(snp) <  1.){
    tnp = snp / TMath::Sqrt((1.-snp)*(1.+snp));
  } 
  Double_t tgl  = tracklet->GetTgl();           // dz/dl
  Double_t dzdx = tgl*TMath::Sqrt(1+tnp*tnp);   // dz/dx calculated from dz/dl
  // at the entrance
  //Double_t tnp = tracklet->GetYref(1);      // dy/dx at the entrance of the chamber
  //Double_t tgl = tracklet->GetZref(1);      // dz/dl at the entrance of the chamber
  //Double_t dzdx = tgl;                      //*TMath::Sqrt(1+tnp*tnp); // dz/dx from dz/dl
  // at the end with correction due to linear fit
  //Double_t tnp = tracklet->GetYfit(1);      // dy/dx at the end of the chamber after fit correction
  //Double_t tgl = tracklet->GetZfit(1);      // dz/dl at the end of the chamber after fit correction 


  ////////////////////////////
  // loop over the clusters
  ////////////////////////////
  Int_t  nbli = 0;
  AliTRDcluster *cl                   = 0x0;
  //////////////////////////////
  // Check no shared clusters
  //////////////////////////////
  for(int icc=AliTRDseedV1::kNtb; icc<AliTRDseedV1::kNclusters; icc++){
    cl = tracklet->GetClusters(icc);
    if(cl)  crossrow = 1;
  }
  //////////////////////////////////
  // Loop clusters
  //////////////////////////////////
  for(int ic=0; ic<AliTRDseedV1::kNtb; ic++){
    if(!(cl = tracklet->GetClusters(ic))) continue;
    if((fLimitChargeIntegration) && (!cl->IsInChamber())) continue;
    
    Double_t ycluster                 = cl->GetY();
    Int_t time                        = cl->GetPadTime();
    Double_t timeis                   = time/fSf;
    //See if cross two pad rows
    Int_t    row                      = cl->GetPadRow();
    if(rowp==-1) rowp                 = row;
    if(row != rowp) crossrow          = 1;

    fLinearFitterTracklet->AddPoint(&timeis,ycluster,1);
    nbli++;  

    
  }
  
  ////////////////////////////////////
  // Do the straight line fit now
  ///////////////////////////////////
  if(nbli <= 2){ 
    fLinearFitterTracklet->ClearPoints();  
    return kFALSE; 
  }
  TVectorD pars;
  fLinearFitterTracklet->Eval();
  fLinearFitterTracklet->GetParameters(pars);
  pointError  =  TMath::Sqrt(fLinearFitterTracklet->GetChisquare()/(nbli-2));
  errorpar    =  fLinearFitterTracklet->GetParError(1)*pointError;
  dydt        = pars[1]; 
  //printf("chis %f, nbli %d, pointError %f, parError %f, errorpar %f\n",fLinearFitterTracklet->GetChisquare(),nbli,pointError,fLinearFitterTracklet->GetParError(1),errorpar);
  fLinearFitterTracklet->ClearPoints();  
 
  ////////////////////////////////
  // Debug stuff
  /////////////////////////////// 


  if(fDebugLevel > 0){
    if ( !fDebugStreamer ) {
      //debug stream
      TDirectory *backup = gDirectory;
      fDebugStreamer = new TTreeSRedirector("TRDdebugCalibraFill.root");
      if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
    } 
    

    Int_t layer = GetLayer(fDetectorPreviousTrack);
           
    (* fDebugStreamer) << "FindP1TrackPHtrackletV1"<<
      //"snpright="<<snpright<<
      "nbli="<<nbli<<
      "nbclusters="<<nbclusters<<
      "detector="<<fDetectorPreviousTrack<<
      "layer="<<layer<<
      "snp="<<snp<<
      "tnp="<<tnp<<
      "tgl="<<tgl<<
      "tnt="<<tnt<<
      "dydt="<<dydt<<
      "dzdx="<<dzdx<<
      "crossrow="<<crossrow<<
      "errorpar="<<errorpar<<
      "pointError="<<pointError<<
      "\n";

  }
  
  /////////////////////////
  // Cuts quality
  ////////////////////////
  
  if(nbclusters < fNumberClusters) return kFALSE;
  if(nbclusters > fNumberClustersf) return kFALSE;
  if(pointError >= 0.3) return kFALSE;
  if(crossrow == 1) return kTRUE;
  
  ///////////////////////
  // Fill
  //////////////////////

  if(fLinearFitterOn){
    //Add to the linear fitter of the detector
    if( TMath::Abs(snp) <  1.){
      Double_t x = tnp-dzdx*tnt; 
      if(fLinearFitterDebugOn) {
	(GetLinearFitter(fDetectorPreviousTrack,kTRUE))->AddPoint(&x,dydt);
	fEntriesLinearFitter[fDetectorPreviousTrack]++;
      }
      fLinearVdriftFit->Update(fDetectorPreviousTrack,x,pars[1]);
    }
  }
  
  return kTRUE;
}
//____________Offine tracking in the AliTRDtracker_____________________________
Bool_t AliTRDCalibraFillHisto::HandlePRFtrackletV1(const AliTRDseedV1 *tracklet, Int_t nbclusters)
{
  //
  // PRF width calibration
  // Assume a Gaussian shape: determinate the position of the three pad clusters
  // Fit with a straight line
  // Take the fitted values for all the clusters (3 or 2 pad clusters)
  // Fill the PRF as function of angle of the track
  //
  //

  //printf("begin\n");
  ///////////////////////////////////////////
  // Take the parameters of the track
  //////////////////////////////////////////
  // take now the snp, tnp and tgl from the track
  Double_t snp = tracklet->GetSnp();             // sin dy/dx at the end of the chamber
  Double_t tnp = 0.0;                            // dy/dx at the end of the chamber 
  if( TMath::Abs(snp) <  1.){
    tnp = snp / TMath::Sqrt((1.-snp)*(1.+snp));
  } 
  Double_t tgl  = tracklet->GetTgl();           // dz/dl
  Double_t dzdx = tgl*TMath::Sqrt(1+tnp*tnp);   // dz/dx calculated from dz/dl
  // at the entrance
  //Double_t tnp = tracklet->GetYref(1);      // dy/dx at the entrance of the chamber
  //Double_t tgl = tracklet->GetZref(1);      // dz/dl at the entrance of the chamber
  //Double_t dzdx = tgl;                      //*TMath::Sqrt(1+tnp*tnp); // dz/dx from dz/dl
  // at the end with correction due to linear fit
  //Double_t tnp = tracklet->GetYfit(1);      // dy/dx at the end of the chamber after fit correction
  //Double_t tgl = tracklet->GetZfit(1);      // dz/dl at the end of the chamber after fit correction 

  ///////////////////////////////
  // Calculate tnp group shift
  ///////////////////////////////
  Bool_t echec   = kFALSE;
  Double_t shift = 0.0;
  //Calculate the shift in x coresponding to this tnp
  if(fNgroupprf != 0.0){
    shift      = -3.0*(fNgroupprf-1)-1.5;
    Double_t limithigh  = -0.2*(fNgroupprf-1);
    if((tnp < (-0.2*fNgroupprf)) || (tnp > (0.2*fNgroupprf))) echec = kTRUE;
    else{
      while(tnp > limithigh){
	limithigh += 0.2;
	shift += 3.0;
      }
    }
  }
  // do nothing if out of tnp range
  //printf("echec %d\n",(Int_t)echec);
  if(echec) return kFALSE;

  ///////////////////////
  // Variables
  //////////////////////

  Int_t nb3pc    = 0;              // number of three pads clusters used for fit 
  // to see the difference between the fit and the 3 pad clusters position
  Double_t padPositions[AliTRDseedV1::kNtb];
  memset(padPositions, 0, AliTRDseedV1::kNtb*sizeof(Double_t)); 
  fLinearFitterTracklet->ClearPoints();  
  
  //printf("loop clusters \n");
  ////////////////////////////
  // loop over the clusters
  ////////////////////////////
  AliTRDcluster *cl                   = 0x0;
  for(int ic=0; ic<AliTRDseedV1::kNtb; ic++){
    // reject shared clusters on pad row
    if((ic+AliTRDseedV1::kNtb) < AliTRDseedV1::kNclusters) {
      cl = tracklet->GetClusters(ic+AliTRDseedV1::kNtb);
      if(cl) continue;
    }
    cl = tracklet->GetClusters(ic);
    if(!cl) continue;
    Double_t     time  = cl->GetPadTime();
    if((time<=7) || (time>=21)) continue;
    Short_t  *signals  = cl->GetSignals(); 
    Float_t xcenter    = 0.0;    
    Bool_t  echec1      = kTRUE;   

    /////////////////////////////////////////////////////////////
    // Center 3 balanced: position with the center of the pad
    /////////////////////////////////////////////////////////////
    if ((((Float_t) signals[3]) > 0.0) && 
	(((Float_t) signals[2]) > 0.0) && 
	(((Float_t) signals[4]) > 0.0)) {
      echec1 = kFALSE;
      // Security if the denomiateur is 0 
      if ((((Float_t) (((Float_t) signals[3]) * ((Float_t) signals[3]))) / 
	   ((Float_t) (((Float_t) signals[2]) * ((Float_t) signals[4])))) != 1.0) {
	xcenter = 0.5 * (TMath::Log((Float_t) (((Float_t) signals[4]) / ((Float_t) signals[2]))))
	  / (TMath::Log(((Float_t) (((Float_t) signals[3]) * ((Float_t) signals[3]))) 
			/ ((Float_t) (((Float_t) signals[2]) * ((Float_t) signals[4])))));
      }
      else {
	echec1 = kTRUE;
      }
    }
    if(TMath::Abs(xcenter) > 0.5) echec1 = kTRUE;
    if(echec1) continue;

    ////////////////////////////////////////////////////////
    //if no echec1: calculate with the position of the pad
    // Position of the cluster
    // fill the linear fitter
    ///////////////////////////////////////////////////////
    Double_t       padPosition = xcenter +  cl->GetPadCol();
    padPositions[ic]            = padPosition;
    nb3pc++;
    fLinearFitterTracklet->AddPoint(&time, padPosition,1);


  }//clusters loop

  //printf("Fin loop clusters \n");
  //////////////////////////////
  // fit with a straight line
  /////////////////////////////
  if(nb3pc < 3){ 
    fLinearFitterTracklet->ClearPoints();  
    return kFALSE;
  }
  fLinearFitterTracklet->Eval();
  TVectorD line(2);
  fLinearFitterTracklet->GetParameters(line);
  Float_t  pointError  = -1.0;
  if( fLinearFitterTracklet->GetChisquare()>=0.0) {
  pointError  =  TMath::Sqrt( fLinearFitterTracklet->GetChisquare()/(nb3pc-2));
  }
  fLinearFitterTracklet->ClearPoints();  
 
  //printf("PRF second loop \n");
  ////////////////////////////////////////////////
  // Fill the PRF: Second loop over clusters
  //////////////////////////////////////////////
  for(int ic=0; ic<AliTRDseedV1::kNtb; ic++){
    // reject shared clusters on pad row
    cl = tracklet->GetClusters(ic+AliTRDseedV1::kNtb);
    if(((ic+AliTRDseedV1::kNtb) < AliTRDseedV1::kNclusters) && (cl)) continue;
    //
    cl = tracklet->GetClusters(ic);
    if(!cl) continue;

    Short_t  *signals      = cl->GetSignals();              // signal
    Double_t     time      = cl->GetPadTime();         // time bin
    Float_t padPosTracklet = line[0]+line[1]*time;          // reconstruct position from fit
    Float_t padPos         = cl->GetPadCol();               // middle pad
    Double_t dpad          = padPosTracklet - padPos;       // reconstruct position relative to middle pad from fit 
    Float_t ycenter        = 0.0;                           // relative center charge
    Float_t ymin           = 0.0;                           // relative left charge
    Float_t ymax           = 0.0;                           // relative right charge
  
    ////////////////////////////////////////////////////////////////
    // Calculate ycenter, ymin and ymax even for two pad clusters
    ////////////////////////////////////////////////////////////////
    if(((((Float_t) signals[3]) > 0.0) && (((Float_t) signals[2]) > 0.0)) ||
       ((((Float_t) signals[3]) > 0.0) && (((Float_t) signals[4]) > 0.0))){
      Float_t sum     = ((Float_t) signals[2]) + ((Float_t) signals[3]) + ((Float_t) signals[4]);
      if(sum > 0.0) ycenter = ((Float_t) signals[3])/ sum;
      if(sum > 0.0) ymin    = ((Float_t) signals[2])/ sum;
      if(sum > 0.0) ymax    = ((Float_t) signals[4])/ sum; 
    }
    
    /////////////////////////
    // Calibration group
    ////////////////////////
    Int_t     rowcl        = cl->GetPadRow();                           // row of cluster
    Int_t     colcl        = cl->GetPadCol();                           // col of cluster 
    Int_t     grouplocal   = CalculateCalibrationGroup(2,rowcl,colcl);  // calcul the corresponding group
    Int_t     caligroup    = fCalibraMode->GetXbins(2)+ grouplocal;     // calcul the corresponding group
    Float_t   xcl          = cl->GetY();                                // y cluster
    Float_t   qcl          = cl->GetQ();                                // charge cluster 
    Int_t     layer        = GetLayer(fDetectorPreviousTrack);          // layer 
    Int_t     stack        = GetStack(fDetectorPreviousTrack);          // stack  
    Double_t  xdiff        = dpad;                                      // reconstructed position constant
    Double_t  x            = dpad;                                      // reconstructed position moved
    Float_t   ep           = pointError;                                // error of fit
    Float_t   signal1      = (Float_t)signals[1];                       // signal at the border
    Float_t   signal3      = (Float_t)signals[3];                       // signal
    Float_t   signal2      = (Float_t)signals[2];                       // signal
    Float_t   signal4      = (Float_t)signals[4];                       // signal
    Float_t   signal5      = (Float_t)signals[5];                       // signal at the border
   


    /////////////////////
    // Debug stuff
    ////////////////////

    if(fDebugLevel > 0){
      if ( !fDebugStreamer ) {
	//debug stream
	TDirectory *backup = gDirectory;
	fDebugStreamer = new TTreeSRedirector("TRDdebugCalibraFill.root");
	if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
      }     
     
      x = xdiff;
      Int_t type=0;
      Float_t y = ycenter;
      (* fDebugStreamer) << "HandlePRFtrackletV1"<<
	"caligroup="<<caligroup<<
	"detector="<<fDetectorPreviousTrack<<
	"layer="<<layer<<
	"stack="<<stack<<
	"npoints="<<nbclusters<<
	"Np="<<nb3pc<<
	"ep="<<ep<<
	"type="<<type<<
       	"snp="<<snp<<
       	"tnp="<<tnp<<
	"tgl="<<tgl<<  
	"dzdx="<<dzdx<< 
	"padPos="<<padPos<<
	"padPosition="<<padPositions[ic]<<
	"padPosTracklet="<<padPosTracklet<<
	"x="<<x<<
	"y="<<y<<	    
	"xcl="<<xcl<<
	"qcl="<<qcl<<
	"signal1="<<signal1<<
	"signal2="<<signal2<<
	"signal3="<<signal3<<
	"signal4="<<signal4<<
	"signal5="<<signal5<<
	"time="<<time<<
	"\n";
      x=-(xdiff+1);
      y = ymin;
      type=-1;
      (* fDebugStreamer) << "HandlePRFtrackletV1"<<
	"caligroup="<<caligroup<<
	"detector="<<fDetectorPreviousTrack<<
	"layer="<<layer<<
	"stack="<<stack<<
	"npoints="<<nbclusters<<
	"Np="<<nb3pc<<
	"ep="<<ep<<
	"type="<<type<<
      	"snp="<<snp<<
      	"tnp="<<tnp<<
	"tgl="<<tgl<<  
	"dzdx="<<dzdx<< 
	"padPos="<<padPos<<
	"padPosition="<<padPositions[ic]<<
	"padPosTracklet="<<padPosTracklet<<
	"x="<<x<<
	"y="<<y<<
	"xcl="<<xcl<<
	"qcl="<<qcl<<
	"signal1="<<signal1<<
	"signal2="<<signal2<<
	"signal3="<<signal3<<
	"signal4="<<signal4<<
	"signal5="<<signal5<<
	"time="<<time<<
	"\n";
      x=1-xdiff;
      y = ymax;
      type=1;
      (* fDebugStreamer) << "HandlePRFtrackletV1"<<
	"caligroup="<<caligroup<<
	"detector="<<fDetectorPreviousTrack<<
	"layer="<<layer<<
	"stack="<<stack<<
	"npoints="<<nbclusters<<
	"Np="<<nb3pc<<
	"ep="<<ep<<
	"type="<<type<<
       	"snp="<<snp<<	
       	"tnp="<<tnp<<	
	"tgl="<<tgl<<  
	"dzdx="<<dzdx<< 
	"padPos="<<padPos<<
	"padPosition="<<padPositions[ic]<<
	"padPosTracklet="<<padPosTracklet<<
	"x="<<x<<
	"y="<<y<<
	"xcl="<<xcl<<
	"qcl="<<qcl<<
	"signal1="<<signal1<<
	"signal2="<<signal2<<
	"signal3="<<signal3<<
	"signal4="<<signal4<<
	"signal5="<<signal5<<
	"time="<<time<<
	"\n";
      
    }
    
    /////////////////////
    // Cuts quality
    /////////////////////
    if(nbclusters < fNumberClusters) continue;
    if(nbclusters > fNumberClustersf) continue;
    if(nb3pc <= 5) continue;
    if((time >= 21) || (time < 7)) continue;
    if(TMath::Abs(qcl) < 80) continue; 
    if( TMath::Abs(snp) >  1.) continue;


    ////////////////////////
    // Fill the histos
    ///////////////////////
    if (fHisto2d) {
      if(TMath::Abs(dpad) < 1.5) {
	fPRF2d->Fill(shift+dpad,(caligroup+0.5),ycenter);
	fPRF2d->Fill(shift-dpad,(caligroup+0.5),ycenter);
	//printf("place %f, ycenter %f\n",(shift+dpad),ycenter);
      }
      if((ymin > 0.0) && (TMath::Abs(dpad+1.0) < 1.5)) {
	fPRF2d->Fill(shift-(dpad+1.0),(caligroup+0.5),ymin);
	fPRF2d->Fill(shift+(dpad+1.0),(caligroup+0.5),ymin);
      }
      if((ymax > 0.0) && (TMath::Abs(dpad-1.0) < 1.5)) {
	fPRF2d->Fill(shift+1.0-dpad,(caligroup+0.5),ymax);
	fPRF2d->Fill(shift-1.0+dpad,(caligroup+0.5),ymax);
      }
    }
    // vector method
    if (fVector2d) {
      if(TMath::Abs(dpad) < 1.5) {
	fCalibraVector->UpdateVectorPRF(fDetectorPreviousTrack,grouplocal,shift+dpad,ycenter);
	fCalibraVector->UpdateVectorPRF(fDetectorPreviousTrack,grouplocal,shift-dpad,ycenter);
      }
      if((ymin > 0.0) && (TMath::Abs(dpad+1.0) < 1.5)) {
	fCalibraVector->UpdateVectorPRF(fDetectorPreviousTrack,grouplocal,shift-(dpad+1.0),ymin);
	fCalibraVector->UpdateVectorPRF(fDetectorPreviousTrack,grouplocal,shift+(dpad+1.0),ymin);
      }
      if((ymax > 0.0)  && (TMath::Abs(dpad-1.0) < 1.5)) {
	fCalibraVector->UpdateVectorPRF(fDetectorPreviousTrack,grouplocal,shift+1.0-dpad,ymax);
	fCalibraVector->UpdateVectorPRF(fDetectorPreviousTrack,grouplocal,shift-1.0+dpad,ymax);
      }
    }
  } // second loop over clusters


  return kTRUE;
}
///////////////////////////////////////////////////////////////////////////////////////
// Pad row col stuff: see if masked or not
///////////////////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::CheckGoodTrackletV1(const AliTRDcluster *cl)
{
  //
  // See if we are not near a masked pad
  //

  if(cl->IsMasked()) fGoodTracklet = kFALSE;

  
}
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::CheckGoodTrackletV0(const Int_t detector,const Int_t row,const Int_t col)
{
  //
  // See if we are not near a masked pad
  //

  if (!IsPadOn(detector, col, row)) {
    fGoodTracklet = kFALSE;
  }

  if (col > 0) {
    if (!IsPadOn(detector, col-1, row)) {
      fGoodTracklet = kFALSE;
    }
  }

  if (col < 143) {
    if (!IsPadOn(detector, col+1, row)) {
      fGoodTracklet = kFALSE;
    }
  }
  
}
//_____________________________________________________________________________
Bool_t AliTRDCalibraFillHisto::IsPadOn(Int_t detector, Int_t row, Int_t col) const
{
  //
  // Look in the choosen database if the pad is On.
  // If no the track will be "not good"
  //

  // Get the parameter object
  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }
  
  if (!cal->IsChamberGood(detector)     || 
       cal->IsChamberNoData(detector)        ||
       cal->IsPadMasked(detector,col,row)) {
    return kFALSE;
  }
  else {
    return kTRUE;
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////
// Calibration groups: calculate the number of groups, localise...
////////////////////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________________
Int_t AliTRDCalibraFillHisto::CalculateCalibrationGroup(Int_t i, Int_t row, Int_t col) const
{
  //
  // Calculate the calibration group number for i
  //
 
  // Row of the cluster and position in the pad groups
  Int_t posr = 0;
  if (fCalibraMode->GetNnZ(i) != 0) {
    posr = (Int_t) row / fCalibraMode->GetNnZ(i);
  }
 
      
  // Col of the cluster and position in the pad groups
  Int_t posc = 0;
  if (fCalibraMode->GetNnRphi(i) != 0) {
    posc = (Int_t) col / fCalibraMode->GetNnRphi(i);
  }
  
  return posc*fCalibraMode->GetNfragZ(i)+posr;
  
}
//____________________________________________________________________________________
Int_t AliTRDCalibraFillHisto::CalculateTotalNumberOfBins(Int_t i)
{
  //
  // Calculate the total number of calibration groups
  //
  
  Int_t ntotal = 0;

  // All together
  if((fCalibraMode->GetNz(i)==100) && (fCalibraMode->GetNrphi(i)==100)){
    ntotal = 1;
    AliInfo(Form("Total number of Xbins: %d for i %d",ntotal,i));
    return ntotal;
  }

  // Per Supermodule
  if((fCalibraMode->GetNz(i)==10) && (fCalibraMode->GetNrphi(i)==10)){
    ntotal = 18;
    AliInfo(Form("Total number of Xbins: %d for i %d",ntotal,i));
    return ntotal;
  }

  // More
  fCalibraMode->ModePadCalibration(2,i);
  fCalibraMode->ModePadFragmentation(0,2,0,i);
  fCalibraMode->SetDetChamb2(i);
  ntotal += 6 * 18 * fCalibraMode->GetDetChamb2(i);
  fCalibraMode->ModePadCalibration(0,i);
  fCalibraMode->ModePadFragmentation(0,0,0,i);
  fCalibraMode->SetDetChamb0(i);
  ntotal += 6 * 4 * 18 * fCalibraMode->GetDetChamb0(i);
  AliInfo(Form("Total number of Xbins: %d for i %d",ntotal,i));
  return ntotal;

}
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::SetNz(Int_t i, Short_t Nz)
{
  //
  // Set the mode of calibration group in the z direction for the parameter i
  // 

  if ((Nz >= 0) && 
      (Nz <  5)) {
    fCalibraMode->SetNz(i, Nz); 
  }
  else { 
    AliInfo("You have to choose between 0 and 4");
  }

}
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::SetNrphi(Int_t i, Short_t Nrphi)
{
  //
  // Set the mode of calibration group in the rphi direction for the parameter i
  //
 
  if ((Nrphi >= 0) && 
      (Nrphi <  7)) {
    fCalibraMode->SetNrphi(i ,Nrphi); 
  }
  else {
    AliInfo("You have to choose between 0 and 6");
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFillHisto::SetAllTogether(Int_t i)
{
  //
  // Set the mode of calibration group all together
  //
  if(fVector2d == kTRUE) {
    AliInfo("Can not work with the vector method");
    return;
  }
  fCalibraMode->SetAllTogether(i);
  
}

//_____________________________________________________________________________
void AliTRDCalibraFillHisto::SetPerSuperModule(Int_t i)
{
  //
  // Set the mode of calibration group per supermodule
  //
  if(fVector2d == kTRUE) {
    AliInfo("Can not work with the vector method");
    return;
  }
  fCalibraMode->SetPerSuperModule(i);
  
}

//____________Set the pad calibration variables for the detector_______________
Bool_t AliTRDCalibraFillHisto::LocalisationDetectorXbins(Int_t detector)
{
  //
  // For the detector calcul the first Xbins and set the number of row
  // and col pads per calibration groups, the number of calibration
  // groups in the detector.
  //
  
  // first Xbins of the detector
  if (fCH2dOn) {
    fCalibraMode->CalculXBins(detector,0);
  }
  if (fPH2dOn) {
    fCalibraMode->CalculXBins(detector,1);
  }
  if (fPRF2dOn) {
    fCalibraMode->CalculXBins(detector,2);
  }

  // fragmentation of idect
  for (Int_t i = 0; i < 3; i++) {
    fCalibraMode->ModePadCalibration((Int_t) GetStack(detector),i);
    fCalibraMode->ModePadFragmentation((Int_t) GetLayer(detector)
                       , (Int_t) GetStack(detector)
                       , (Int_t) GetSector(detector),i);
  }
  
  return kTRUE;

}
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::SetNumberGroupsPRF(Short_t numberGroupsPRF)
{
  //
  // Should be between 0 and 6
  //
 
  if ((numberGroupsPRF < 0) || (numberGroupsPRF > 6)) {
    AliInfo("The number of groups must be between 0 and 6!");
  } 
  else {
    fNgroupprf = numberGroupsPRF;
  }

} 
///////////////////////////////////////////////////////////////////////////////////////////
// Per tracklet: store or reset the info, fill the histos with the info
//////////////////////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::StoreInfoCHPHtrack(const AliTRDcluster *cl,const Double_t dqdl,const Int_t *group,const Int_t row,const Int_t col,const AliTRDcluster *cls)
{
  //
  // Store the infos in fAmpTotal, fPHPlace and fPHValue
  // Correct from the gain correction before
  // cls is shared cluster if any
  //
  
  //printf("StoreInfoCHPHtrack\n");

  // time bin of the cluster not corrected
  Int_t    time     = cl->GetPadTime();
  Float_t  charge   = TMath::Abs(cl->GetQ());  
  if(cls) {
    charge += TMath::Abs(cls->GetQ());
    //printf("AliTRDCalibraFillHisto::Add the cluster charge");
  }

  //printf("Store::time %d, amplitude %f\n",time,dqdl);
  
  //Correct for the gain coefficient used in the database for reconstruction
  Float_t correctthegain = 1.0;
  if(fIsHLT) correctthegain = fCalDetGain->GetValue(fDetectorPreviousTrack);
  else correctthegain = fCalDetGain->GetValue(fDetectorPreviousTrack)*fCalROCGain->GetValue(col,row);
  Float_t correction    = 1.0;
  Float_t normalisation = 6.67;
  // we divide with gain in AliTRDclusterizer::Transform...
  if( correctthegain > 0 ) normalisation /= correctthegain;


  // take dd/dl corrected from the angle of the track
  correction = dqdl / (normalisation);
  

  // Fill the fAmpTotal with the charge
  if (fCH2dOn) {
    if((!fLimitChargeIntegration) || (cl->IsInChamber())) {
      //printf("Store::group %d, amplitude %f\n",group[0],correction);
      fAmpTotal[(Int_t) group[0]] += correction;
    }
  }

  // Fill the fPHPlace and value
  if (fPH2dOn) {
    fPHPlace[time] = group[1];
    fPHValue[time] = charge;
  }
  
}
//____________Offine tracking in the AliTRDtracker_____________________________
void AliTRDCalibraFillHisto::ResetfVariablestracklet()
{
  //
  // Reset values per tracklet
  //

  //Reset good tracklet
  fGoodTracklet = kTRUE;

  // Reset the fPHValue
  if (fPH2dOn) {
    //Reset the fPHValue and fPHPlace
    for (Int_t k = 0; k < fTimeMax; k++) {
      fPHValue[k] = 0.0;
      fPHPlace[k] = -1;
    }
  }

  // Reset the fAmpTotal where we put value
  if (fCH2dOn) {
    for (Int_t k = 0; k < fCalibraMode->GetNfragZ(0)*fCalibraMode->GetNfragRphi(0); k++) {
      fAmpTotal[k] = 0.0;
    }
  }
}
//____________Offine tracking in the AliTRDtracker_____________________________
void AliTRDCalibraFillHisto::FillTheInfoOfTheTrackCH(Int_t nbclusters)
{
  //
  // For the offline tracking or mcm tracklets
  // This function will be called in the functions UpdateHistogram... 
  // to fill the info of a track for the relativ gain calibration
  //
	
  Int_t nb            =  0;   // Nombre de zones traversees
  Int_t fd            = -1;   // Premiere zone non nulle
  Float_t totalcharge = 0.0;  // Total charge for the supermodule histo

  //printf("CH2d nbclusters %d, fNumberClusters %d, fNumberClustersf %d\n",nbclusters,fNumberClusters,fNumberClustersf);

  if(nbclusters < fNumberClusters) return;
  if(nbclusters > fNumberClustersf) return;


  // Normalize with the number of clusters
  Double_t normalizeCst = fRelativeScale;
  if(fNormalizeNbOfCluster) normalizeCst = normalizeCst*nbclusters;

  //printf("Number of groups in one detector %d\n",fCalibraMode->GetNfragZ(0)*fCalibraMode->GetNfragRphi(0));
  
  // See if the track goes through different zones
  for (Int_t k = 0; k < fCalibraMode->GetNfragZ(0)*fCalibraMode->GetNfragRphi(0); k++) {
    //printf("fAmpTotal %f for %d\n",fAmpTotal[k],k);
    if (fAmpTotal[k] > 0.0) {
      totalcharge += fAmpTotal[k];
      nb++;
      if (nb == 1) {
        fd = k;
      }
    }
  }

  //printf("CH2d: nb %d, fd %d, calibration group %d, amplitude %f, detector %d\n",nb,fd,fCalibraMode->GetXbins(0),fAmpTotal[fd]/normalizeCst,fDetectorPreviousTrack);
    
  switch (nb)
    { 
    case 1:
      fNumberUsedCh[0]++;
      fEntriesCH[fCalibraMode->GetXbins(0)+fd]++;
      if (fHisto2d) {
	FillCH2d(fCalibraMode->GetXbins(0)+fd,fAmpTotal[fd]/normalizeCst);
	//fCH2d->Fill(fAmpTotal[fd]/normalizeCst,fCalibraMode->GetXbins(0)+fd+0.5);
      }
      if (fVector2d) {
	fCalibraVector->UpdateVectorCH(fDetectorPreviousTrack,fd,fAmpTotal[fd]/normalizeCst);
      }
      break;
    case 2:
      if ((fAmpTotal[fd]   > 0.0) && 
	  (fAmpTotal[fd+1] > 0.0)) {
	// One of the two very big
	if (fAmpTotal[fd] > fProcent*fAmpTotal[fd+1]) {
	  if (fHisto2d) {
	    FillCH2d(fCalibraMode->GetXbins(0)+fd,fAmpTotal[fd]/normalizeCst);
	    //fCH2d->Fill(fAmpTotal[fd]/normalizeCst,fCalibraMode->GetXbins(0)+fd+0.5);
	  }
	  if (fVector2d) {
	    fCalibraVector->UpdateVectorCH(fDetectorPreviousTrack,fd,fAmpTotal[fd]/normalizeCst);
	  }
	  fNumberUsedCh[1]++;
	  fEntriesCH[fCalibraMode->GetXbins(0)+fd]++;
	}
	if (fAmpTotal[fd+1] > fProcent*fAmpTotal[fd]) {
	  if (fHisto2d) {
	    FillCH2d(fCalibraMode->GetXbins(0)+fd+1,fAmpTotal[fd+1]/normalizeCst);
	    //fCH2d->Fill(fAmpTotal[fd+1]/normalizeCst,fCalibraMode->GetXbins(0)+fd+1.5);
	  }
	  if (fVector2d) {
	    fCalibraVector->UpdateVectorCH(fDetectorPreviousTrack,fd+1,fAmpTotal[fd+1]/normalizeCst);
	  }
	  fNumberUsedCh[1]++;
	  fEntriesCH[fCalibraMode->GetXbins(0)+fd+1]++;
	}
      }
      if (fCalibraMode->GetNfragZ(0) > 1) {
	if (fAmpTotal[fd] > 0.0) {
	  if ((fd+fCalibraMode->GetNfragZ(0)) < (fCalibraMode->GetNfragZ(0)*fCalibraMode->GetNfragRphi(0))) {
	    if (fAmpTotal[fd+fCalibraMode->GetNfragZ(0)] > 0.0) {
	      // One of the two very big
	      if (fAmpTotal[fd] > fProcent*fAmpTotal[fd+fCalibraMode->GetNfragZ(0)]) {
		if (fHisto2d) {
		  FillCH2d(fCalibraMode->GetXbins(0)+fd,fAmpTotal[fd]/normalizeCst);
		  //fCH2d->Fill(fAmpTotal[fd]/normalizeCst,fCalibraMode->GetXbins(0)+fd+0.5);
		}
		if (fVector2d) {
		  fCalibraVector->UpdateVectorCH(fDetectorPreviousTrack,fd,fAmpTotal[fd]/normalizeCst);
		}
		fNumberUsedCh[1]++;
		fEntriesCH[fCalibraMode->GetXbins(0)+fd]++;
	      }
	      if (fAmpTotal[fd+fCalibraMode->GetNfragZ(0)] > fProcent*fAmpTotal[fd]) {
		if (fHisto2d) {
		  FillCH2d(fCalibraMode->GetXbins(0)+fd+fCalibraMode->GetNfragZ(0),fAmpTotal[fd+fCalibraMode->GetNfragZ(0)]/normalizeCst);
		  //fCH2d->Fill(fAmpTotal[fd+fCalibraMode->GetNfragZ(0)]/normalizeCst,fCalibraMode->GetXbins(0)+fd+fCalibraMode->GetNfragZ(0)+0.5);
		}
		fNumberUsedCh[1]++;
		fEntriesCH[fCalibraMode->GetXbins(0)+fd+fCalibraMode->GetNfragZ(0)]++;
		if (fVector2d) {
		  fCalibraVector->UpdateVectorCH(fDetectorPreviousTrack,fd+fCalibraMode->GetNfragZ(0),fAmpTotal[fd+fCalibraMode->GetNfragZ(0)]/normalizeCst);
		}
	      }
	    }
	  }
	}
      }
      break;
    default: break;
  }
}
//____________Offine tracking in the AliTRDtracker_____________________________
void AliTRDCalibraFillHisto::FillTheInfoOfTheTrackPH()
{
  //
  // For the offline tracking or mcm tracklets
  // This function will be called in the functions UpdateHistogram... 
  // to fill the info of a track for the drift velocity  calibration
  //
    
  Int_t nb  =  1; // Nombre de zones traversees 1, 2 ou plus de 3
  Int_t fd1 = -1; // Premiere zone non nulle
  Int_t fd2 = -1; // Deuxieme zone non nulle
  Int_t k1  = -1; // Debut de la premiere zone
  Int_t k2  = -1; // Debut de la seconde zone
  Int_t nbclusters = 0; // number of clusters



  // See if the track goes through different zones
  for (Int_t k = 0; k < fTimeMax; k++) {
    if (fPHValue[k] > 0.0) {
      nbclusters++;
      if (fd1 == -1) {
	fd1 = fPHPlace[k];
	k1  = k;	      
      }
      if (fPHPlace[k] != fd1) {
	if (fd2 == -1) {
	  k2  = k;
	  fd2 = fPHPlace[k];
	  nb  = 2;
	}
	if (fPHPlace[k] != fd2) {
          nb = 3;
	}
      }
    }
  }

  // See if noise before and after
  if(fMaxCluster > 0) {
    if(fPHValue[0] > fMaxCluster) return;
    if(fTimeMax > fNbMaxCluster) {
      for(Int_t k = (fTimeMax-fNbMaxCluster); k < fTimeMax; k++){
	if(fPHValue[k] > fMaxCluster) return;
      }
    }
  }

  //printf("nbclusters %d, low limit %d, high limit %d\n",nbclusters,fNumberClusters,fNumberClustersf);

  if(nbclusters < fNumberClusters) return;
  if(nbclusters > fNumberClustersf) return;
  
  switch(nb)
    {
    case 1:
      fNumberUsedPh[0]++;
      for (Int_t i = 0; i < fTimeMax; i++) {
	if (fHisto2d) {
	  if(fFillWithZero) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) fPHValue[i]);
	  else {
	    if(((Float_t) fPHValue[i] > 0.0)) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) fPHValue[i]);
	      }
	  //printf("Fill the time bin %d with %f\n",i,fPHValue[i]);
	}
	if (fVector2d) {
	  if(fFillWithZero) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd1,i,fPHValue[i]);
	  else {
	    if(((Float_t) fPHValue[i] > 0.0)) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd1,i,fPHValue[i]);  
	  }
	}
      }
      break;
    case 2:
      if ((fd1 == fd2+1) || 
	  (fd2 == fd1+1)) {
	// One of the two fast all the think
	if (k2 > (k1+fDifference)) {
	  //we choose to fill the fd1 with all the values
	  fNumberUsedPh[1]++;
	  for (Int_t i = 0; i < fTimeMax; i++) {
	    if (fHisto2d) {
	      if(fFillWithZero) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) fPHValue[i]);
	      else {
		if(((Float_t) fPHValue[i] > 0.0)) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) fPHValue[i]);
		  }
	    }
	    if (fVector2d) {
	      if(fFillWithZero) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd1,i,fPHValue[i]);
	      else {
		if(((Float_t) fPHValue[i] > 0.0)) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd1,i,fPHValue[i]);
		  }
	    }
	  }
	}
	if ((k2+fDifference) < fTimeMax) {
	  //we choose to fill the fd2 with all the values
	  fNumberUsedPh[1]++;
	  for (Int_t i = 0; i < fTimeMax; i++) {
	    if (fHisto2d) {
	      if(fFillWithZero) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd2)+0.5,(Float_t) fPHValue[i]);
	      else {
		if(((Float_t) fPHValue[i] > 0.0)) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd2)+0.5,(Float_t) fPHValue[i]);
	      }
	    }
	  if (fVector2d) {
	    if(fFillWithZero) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd2,i,fPHValue[i]);
	    else {
	      if(((Float_t) fPHValue[i] > 0.0)) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd2,i,fPHValue[i]);
	    }
	  }
	  }
	}
      }
      // Two zones voisines sinon rien!
      if (fCalibraMode->GetNfragZ(1) > 1) {
	// Case 2
	if ((fd1+fCalibraMode->GetNfragZ(1)) < (fCalibraMode->GetNfragZ(1)*fCalibraMode->GetNfragRphi(1))) {
	  if (fd2 == (fd1+fCalibraMode->GetNfragZ(1))) {
	    // One of the two fast all the think
	    if (k2 > (k1+fDifference)) {
	      //we choose to fill the fd1 with all the values
	      fNumberUsedPh[1]++;
	      for (Int_t i = 0; i < fTimeMax; i++) {
		if (fHisto2d) {
		  if(fFillWithZero) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) fPHValue[i]);
		  else {
		    if(((Float_t) fPHValue[i] > 0.0)) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) fPHValue[i]);
		  }
		}
		if (fVector2d) {
		  if(fFillWithZero) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd1,i,fPHValue[i]);
		  else {
		    if(((Float_t) fPHValue[i] > 0.0)) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd1,i,fPHValue[i]);
		  }
		}
	      }
	    }
	    if ((k2+fDifference) < fTimeMax) {
	      //we choose to fill the fd2 with all the values
	      fNumberUsedPh[1]++;
	      for (Int_t i = 0; i < fTimeMax; i++) {
		if (fHisto2d) {
		  if(fFillWithZero) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd2)+0.5,(Float_t) fPHValue[i]);
		  else {
		    if(((Float_t) fPHValue[i] > 0.0)) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd2)+0.5,(Float_t) fPHValue[i]);
		  }
		}
		if (fVector2d) {
		  if(fFillWithZero) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd2,i,fPHValue[i]);
		  else {
		    if(((Float_t) fPHValue[i] > 0.0)) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd2,i,fPHValue[i]);
		  }
		}
	      }
	    }
	  }
	}
	// Two zones voisines sinon rien!
	// Case 3
	if ((fd1 - fCalibraMode->GetNfragZ(1)) >= 0) {
	  if (fd2 == (fd1 - fCalibraMode->GetNfragZ(1))) {
	    // One of the two fast all the think
	    if (k2 > (k1 + fDifference)) {
	      //we choose to fill the fd1 with all the values
	      fNumberUsedPh[1]++;
	      for (Int_t i = 0; i < fTimeMax; i++) {
		if (fHisto2d) {
		  if(fFillWithZero) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) fPHValue[i]);
		  else {
		    if(((Float_t) fPHValue[i] > 0.0)) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) fPHValue[i]);
		  }
		}
		if (fVector2d) {
		  if(fFillWithZero) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd1,i,fPHValue[i]);
		  else {
		    if(((Float_t) fPHValue[i] > 0.0)) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd1,i,fPHValue[i]);
		  }
		}
	      }
	    }
	    if ((k2+fDifference) < fTimeMax) {
	      //we choose to fill the fd2 with all the values
	      fNumberUsedPh[1]++;
	      for (Int_t i = 0; i < fTimeMax; i++) {
		if (fHisto2d) {
		  if(fFillWithZero) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd2)+0.5,(Float_t) fPHValue[i]);
		  else {
		    if(((Float_t) fPHValue[i] > 0.0)) fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd2)+0.5,(Float_t) fPHValue[i]);
		  }
		}
		if (fVector2d) {
		  if(fFillWithZero) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd2,i,fPHValue[i]);
		  else {
		    if(((Float_t) fPHValue[i] > 0.0)) fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd2,i,fPHValue[i]);
		  }
		}
	      }
	    }
	  }
	}
      }
      break;
    default: break;
    } 
}
//////////////////////////////////////////////////////////////////////////////////////////
// DAQ process functions
/////////////////////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________
Int_t AliTRDCalibraFillHisto::ProcessEventDAQ(AliRawReader *rawReader)
 { //main
  //
  // Event Processing loop - AliTRDrawStream
  // 
  // 0 timebin problem
  // 1 no input
  // 2 input
  // Same algorithm as TestBeam but different reader
  //

  AliTRDrawStream *rawStream = new AliTRDrawStream(rawReader);

  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager(kTRUE);
  digitsManager->CreateArrays();
    
  Int_t withInput = 1;
  
  Double_t phvalue[16][144][36];
  for(Int_t k = 0; k < 36; k++){
    for(Int_t j = 0; j < 16; j++){
      for(Int_t c = 0; c < 144; c++){
	phvalue[j][c][k] = 0.0;
      }
    }
  }
  
  fDetectorPreviousTrack = -1;
  fMCMPrevious           = -1;
  fROBPrevious           = -1;
  
  Int_t nbtimebin = 0;                                        
  Int_t baseline  = 10;  

  
    fTimeMax = 0;
       
    Int_t det    = 0;
    while ((det = rawStream->NextChamber(digitsManager, NULL, NULL)) >= 0) { //idetector

      if (digitsManager->GetIndexes(det)->HasEntry()) {//QA
	//	printf("there is ADC data on this chamber!\n");

	AliTRDarrayADC *digits = (AliTRDarrayADC *) digitsManager->GetDigits(det); //mod
	if (digits->HasData()) { //array
	  
	  AliTRDSignalIndex   *indexes = digitsManager->GetIndexes(det);
	  if (indexes->IsAllocated() == kFALSE) {
	    AliError("Indexes do not exist!");
	    break;
	  }
	  Int_t iRow  = 0;
	  Int_t iCol  = 0;
	  indexes->ResetCounters();
	  
	  while (indexes->NextRCIndex(iRow, iCol)) { //column,row
	    //printf(" det %d \t row %d \t col %d \t digit\n",det,iRow,iCol);
	    //while (rawStream->Next()) {
	    
	    Int_t idetector = det;                                             //  current detector
	    //Int_t imcm      = rawStream->GetMCM();                            //  current MCM
	    //Int_t irob      = rawStream->GetROB();                            //  current ROB
	    
	  
	    if((fDetectorPreviousTrack != idetector) && (fDetectorPreviousTrack != -1)) {
	      // Fill
	      withInput = TMath::Max(FillDAQ(phvalue),withInput);
	      
	      // reset
	      for(Int_t k = 0; k < 36; k++){
		for(Int_t j = 0; j < 16; j++){
		  for(Int_t c = 0; c < 144; c++){
		    phvalue[j][c][k] = 0.0;
		  }
		}
	      }
	    }
	    
	    fDetectorPreviousTrack = idetector;
	    //fMCMPrevious           = imcm;
	    //fROBPrevious           = irob;
	    
	    //	  nbtimebin              = rawStream->GetNumberOfTimeBins();              //  number of time bins read from data
	    AliTRDdigitsParam *digitParam = (AliTRDdigitsParam *)digitsManager->GetDigitsParam();
	    nbtimebin              = digitParam->GetNTimeBins(det);              //  number of time bins read from data
	    baseline               = digitParam->GetADCbaseline(det);            //  baseline
	    
	    if(nbtimebin == 0) return 0;
	    if((fTimeMax != 0) && (nbtimebin != fTimeMax)) return 0;
	    fTimeMax          = nbtimebin;
	    
	    fNumberClustersf    = fTimeMax;
	    fNumberClusters     = (Int_t)(fNumberClustersProcent*fTimeMax);
	  	  
	    
	    for(Int_t itime = 0; itime < nbtimebin; itime++) {
	      //	    phvalue[row][col][itime] = signal[itime]-baseline;
	      phvalue[iRow][iCol][itime] = (Short_t)(digits->GetData(iRow,iCol,itime) - baseline);
	      /*if(phvalue[iRow][iCol][itime] >= 20) {
		 printf("----------> phvalue[%d][%d][%d] %d  baseline %d \n",
		       iRow,
		       iCol,
		       itime,
		       (Short_t)(digits->GetData(iRow,iCol,itime)),
		       baseline);
		       }*/
	    }
	    
	  }//column,row
	  
	  // fill the last one
	  if(fDetectorPreviousTrack != -1){
	    
	    // Fill
	    withInput = TMath::Max(FillDAQ(phvalue),withInput);
	    //	    printf("\n ---> withinput %d\n\n",withInput);
	    // reset
	    for(Int_t k = 0; k < 36; k++){
	      for(Int_t j = 0; j < 16; j++){
		for(Int_t c = 0; c < 144; c++){
		  phvalue[j][c][k] = 0.0;
		}
	      }
	    }
	  }
	  
	}//array
      }//QA
      digitsManager->ClearArrays(det);
    }//idetector
    delete digitsManager;

    delete rawStream;
    return withInput;
 }//main
//_____________________________________________________________________
//////////////////////////////////////////////////////////////////////////////
// Routine inside the DAQ process
/////////////////////////////////////////////////////////////////////////////
//_______________________________________________________________________
Int_t AliTRDCalibraFillHisto::FillDAQ(Double_t phvalue[16][144][36]){

  //
  // Look for the maximum by collapsing over the time
  // Sum over four pad col and two pad row
  //

  Int_t used = 0;


  Int_t idect = fDetectorPreviousTrack;      
  //printf("Enter Detector %d\n",fDetectorPreviousTrack);
  Double_t sum[36];
  for(Int_t tb = 0; tb < 36; tb++){
    sum[tb] = 0.0;
  }

  //fGoodTracklet = kTRUE;
  //fDetectorPreviousTrack = 0;  


  ///////////////////////////
  // look for maximum
  /////////////////////////

  Int_t imaxRow = 0;
  Int_t imaxCol = 0;
  Double_t integralMax = -1;
  
  for (Int_t ir = 1; ir <= 15; ir++)
    {
      for (Int_t ic = 2; ic <= 142; ic++)
	{
	  Double_t integral = 0;		  
	  for (Int_t ishiftR = 0; ishiftR < fNumberRowDAQ; ishiftR++)
	    {
	      for (Int_t ishiftC = -fNumberColDAQ; ishiftC < fNumberColDAQ; ishiftC++)
		{
		  if (ir + ishiftR >= 1 && ir + ishiftR <= 16 &&
		      ic + ishiftC >= 1 && ic + ishiftC <= 144)
		    {

		      for(Int_t tb = 0; tb< fTimeMax; tb++){
			integral += phvalue[ir + ishiftR-1][ic + ishiftC-1][tb];
		      }// addtb
		    } //addsignal
		} //shiftC
	    } // shiftR
	  if (integralMax < integral)
	    {
	      imaxRow = ir;
	      imaxCol = ic;
	      integralMax = integral;
	      
	    } // check max integral
	} //ic
    } // ir

  //  printf("imaxRow %d, imaxCol %d, fTimeMax %d, integralMax %f\n",imaxRow,imaxCol,fTimeMax, integralMax);
  //if((imaxRow == 0) || (imaxRow >= 15) || (imaxCol <= 3) || (imaxCol >= 140)) {
  //  used=1;
  //  return used;
  // }

  if(((imaxRow + fNumberRowDAQ) > 16) || (imaxRow == 0) || ((imaxCol - fNumberColDAQ) <= 1) || ((imaxCol + fNumberColDAQ) >= 144)) {
    used=1;
    return used;
  }
  //CheckGoodTrackletV0(fDetectorPreviousTrack,imaxRow,imaxCol);
  //if(!fGoodTracklet) used = 1;;
  
  //  /////////////////////////////////////////////////////
  // sum ober 2 row and 4 pad cols for each time bins
  //  ////////////////////////////////////////////////////	  
  
  
  
  for (Int_t ishiftR = 0; ishiftR < fNumberRowDAQ; ishiftR++)
    {
      for (Int_t ishiftC = -fNumberColDAQ; ishiftC < fNumberColDAQ; ishiftC++)
	{
	  if (imaxRow + ishiftR >= 1 && imaxRow + ishiftR <= 16 &&
	      imaxCol + ishiftC >= 1 && imaxCol + ishiftC <= 144)
	    { 
	      for(Int_t it = 0; it < fTimeMax; it++){
		sum[it] += phvalue[imaxRow + ishiftR-1][imaxCol + ishiftC-1][it];
	      } 
	    }
	} // col shift
    }// row shift

  Int_t nbcl = 0;
  Double_t sumcharge = 0.0;
  for(Int_t it = 0; it < fTimeMax; it++){
    sumcharge += sum[it];
    if(sum[it] > fThresholdClustersDAQ)  nbcl++;
  }


  /////////////////////////////////////////////////////////
  // Debug
  ////////////////////////////////////////////////////////
  if(fDebugLevel > 0){
    if ( !fDebugStreamer ) {
      //debug stream
      TDirectory *backup = gDirectory;
      fDebugStreamer = new TTreeSRedirector("TRDdebugCalibraFill.root");
      if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
    }     

    Double_t amph0 = sum[0];
    Double_t amphlast = sum[fTimeMax-1];
    Double_t rms      = TMath::RMS(fTimeMax,sum);
    Int_t    goodtracklet = (Int_t) fGoodTracklet;
    for(Int_t it = 0; it < fTimeMax; it++){
      Double_t clustera = sum[it]; 

    (* fDebugStreamer) << "FillDAQa"<<
      "ampTotal="<<sumcharge<<
      "row="<<imaxRow<<
      "col="<<imaxCol<<
      "detector="<<idect<<
      "amph0="<<amph0<<
      "amphlast="<<amphlast<<
      "goodtracklet="<<goodtracklet<<
      "clustera="<<clustera<<
      "it="<<it<<
      "rms="<<rms<<
      "nbcl="<<nbcl<<
      "\n"; 
    }
  }

  ////////////////////////////////////////////////////////
  // fill
  ///////////////////////////////////////////////////////
  //printf("fNumberClusters %d, fNumberClustersf %d\n",fNumberClusters,fNumberClustersf);
  if(sum[0] > 100.0) used = 1; 
  if(nbcl < fNumberClusters) used = 1;
  if(nbcl > fNumberClustersf) used = 1;

  //if(fDetectorPreviousTrack == 15){
  //  printf("rms %f and first time bin %f\n",TMath::RMS(fTimeMax,sum),sum[0]);
  //}
  //if((TMath::RMS(fTimeMax,sum) <= 10.0) && (sum[0] > 200.0)) return 1;
  if(used == 0){
    for(Int_t it = 0; it < fTimeMax; it++){
      if(fFillWithZero) UpdateDAQ(fDetectorPreviousTrack,0,0,it,sum[it],fTimeMax); 
      else{
	if(sum[it] > 0.0) UpdateDAQ(fDetectorPreviousTrack,0,0,it,sum[it],fTimeMax); 
      } 
      //if(fFillWithZero) UpdateDAQ(0,0,0,it,sum[it],fTimeMax);
      //else{
      // if(sum[it] > 0.0) UpdateDAQ(0,0,0,it,sum[it],fTimeMax);
      //}
    }
    
   
    //((TH2I *)GetCH2d()->Fill(sumcharge/30.0,fDetectorPreviousTrack));
    used = 2;
    //printf("Pass Detector %d\n",fDetectorPreviousTrack);

  }
 
  return used;
  
}
//____________Online trackling in AliTRDtrigger________________________________
Bool_t AliTRDCalibraFillHisto::UpdateDAQ(Int_t det, Int_t /*row*/, Int_t /*col*/, Int_t timebin, Float_t signal, Int_t nbtimebins)
{
  //
  // For the DAQ
  // Fill a simple average pulse height
  //

  
  ((TProfile2D *)GetPH2d(nbtimebins,fSf))->Fill((Float_t) timebin/fSf,det+0.5,(Float_t) signal); 

 
  return kTRUE;
  
}
//____________Write_____________________________________________________
//_____________________________________________________________________
void AliTRDCalibraFillHisto::Write2d(const Char_t *filename, Bool_t append)
{
  //
  //  Write infos to file
  //
  
  //For debugging
  if ( fDebugStreamer ) {
    delete fDebugStreamer;
    fDebugStreamer = 0x0;
  }

  AliInfo(Form("Numbertrack: %d Numberusedch[0]: %d, Numberusedch[1]: %d Numberusedph[0]: %d, Numberusedph[1]: %d"
	       ,fNumberTrack
	       ,fNumberUsedCh[0]
	       ,fNumberUsedCh[1]
	       ,fNumberUsedPh[0]
	       ,fNumberUsedPh[1]));
  
  TDirectory *backup = gDirectory;
  TString option;
  
  if ( append )
    option = "update";
  else
    option = "recreate";
  
  TFile f(filename,option.Data());
  
  TStopwatch stopwatch;
  stopwatch.Start();
  if(fVector2d) {
    f.WriteTObject(fCalibraVector);
  }

  if (fCH2dOn ) {
    if (fHisto2d) {
      f.WriteTObject(fCH2d);
    }
  }
  if (fPH2dOn ) {
    if (fHisto2d) {
      f.WriteTObject(fPH2d);
    }
  }
  if (fPRF2dOn) {
    if (fHisto2d) {
	f.WriteTObject(fPRF2d);
    }
  }
  if(fLinearFitterOn){
    if(fLinearFitterDebugOn) AnalyseLinearFitter();
    f.WriteTObject(fLinearVdriftFit);
  }
   
  f.Close();
  
  if ( backup ) backup->cd();
  
  AliInfo(Form("Execution time Write2d: R:%.2fs C:%.2fs"
	       ,stopwatch.RealTime(),stopwatch.CpuTime()));
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Stats stuff
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//___________________________________________probe the histos__________________________________________________
Double_t *AliTRDCalibraFillHisto::StatH(TH2 *h, Int_t i)
{
  //
  // Check the number of stats in h, 0 is TH2I 1 is TProfile2D
  // debug mode with 2 for TH2I and 3 for TProfile2D
  // It gives a pointer to a Double_t[7] with the info following...
  // [0] : number of calibration groups with entries
  // [1] : minimal number of entries found
  // [2] : calibration group number of the min
  // [3] : maximal number of entries found
  // [4] : calibration group number of the max
  // [5] : mean number of entries found
  // [6] : mean relative error
  //

  Double_t *info = new Double_t[7];
   
  // Number of Xbins (detectors or groups of pads)
  Int_t    nbins   = h->GetNbinsY(); //number of calibration groups
  Int_t    nxbins  = h->GetNbinsX(); //number of bins per histo

  // Initialise
  Double_t nbwe = 0; //number of calibration groups with entries
  Double_t minentries = 0; //minimal number of entries found
  Double_t maxentries = 0; //maximal number of entries found
  Double_t placemin = 0; //calibration group number of the min
  Double_t placemax = -1; //calibration group number of the max
  Double_t meanstats = 0.0; //mean number of entries over the calibration group with at least ome entry
  Double_t meanrelativerror = 0.0; //mean relativ error in the TProfile2D

  Double_t counter = 0;

  //Debug
  TH1F *nbEntries = 0x0;//distribution of the number of entries
  TH1F *nbEntriesPerGroup = 0x0;//Number of entries per group
  TProfile *nbEntriesPerSp = 0x0;//Number of entries for one supermodule
    
  // Beginning of the loop over the calibration groups 
  for (Int_t idect = 0; idect < nbins; idect++) {

    TH1I *projch = (TH1I *) h->ProjectionX("projch",idect+1,idect+1,(Option_t *)"e");
    projch->SetDirectory(0);
    
    // Number of entries for this calibration group
    Double_t nentries = 0.0;
    if((i%2) == 0){
      for (Int_t k = 0; k < nxbins; k++) {
	nentries += h->GetBinContent(h->GetBin(k+1,idect+1));
      }
    }
    else{
      for (Int_t k = 0; k < nxbins; k++) {
	nentries += ((TProfile2D *)h)->GetBinEntries(h->GetBin(k+1,idect+1));
	if(h->GetBinContent(h->GetBin(k+1,idect+1)) != 0) {
	  meanrelativerror += (h->GetBinError(h->GetBin(k+1,idect+1))/(TMath::Abs(h->GetBinContent(h->GetBin(k+1,idect+1)))));
	  counter++;
	} 
      }
    }

    //Debug
    if(i > 1){
      if(nentries > 0){
	if(!((Bool_t)nbEntries)) nbEntries = new TH1F("Number of entries","Number of entries",100,(Int_t)nentries/2,nentries*2);
	nbEntries->SetDirectory(0);
	nbEntries->Fill(nentries);
	if(!((Bool_t)nbEntriesPerGroup)) nbEntriesPerGroup = new TH1F("Number of entries per group","Number of entries per group",nbins,0,nbins);
	nbEntriesPerGroup->SetDirectory(0);
	nbEntriesPerGroup->Fill(idect+0.5,nentries);
	if(!((Bool_t)nbEntriesPerSp)) nbEntriesPerSp = new TProfile("Number of entries per supermodule","Number of entries per supermodule",(Int_t)(nbins/18),0,(Int_t)(nbins/18));
	nbEntriesPerSp->SetDirectory(0);
	nbEntriesPerSp->Fill((idect%((Int_t)(nbins/18)))+0.5,nentries);
      }
    }

    //min amd max
    if(nentries > maxentries){
      maxentries = nentries;
      placemax = idect;
    }
    if(idect == 0) {
      minentries = nentries;
    }
    if(nentries < minentries){
      minentries = nentries;
      placemin = idect;
    }
    //nbwe
    if(nentries > 0) {
      nbwe++;
      meanstats += nentries;
    }
  }//calibration groups loop
  
  if(nbwe > 0) meanstats /= nbwe;
  if(counter > 0) meanrelativerror /= counter;

  AliInfo(Form("There are %f calibration groups with entries",nbwe));
  AliInfo(Form("The minimum number of entries is %f for the group %f",minentries,placemin));
  AliInfo(Form("The maximum number of entries is %f for the group %f",maxentries,placemax));
  AliInfo(Form("The mean number of entries is %f",meanstats));
  if((i%2) == 1) AliInfo(Form("The mean relative error is %f",meanrelativerror));

  info[0] = nbwe;
  info[1] = minentries;
  info[2] = placemin;
  info[3] = maxentries;
  info[4] = placemax;
  info[5] = meanstats;
  info[6] = meanrelativerror;

  if(nbEntries && nbEntriesPerSp && nbEntriesPerGroup){
    gStyle->SetPalette(1);
    gStyle->SetOptStat(1111);
    gStyle->SetPadBorderMode(0);
    gStyle->SetCanvasColor(10);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.01);
    TCanvas *stat = new TCanvas("stat","",50,50,600,800);
    stat->Divide(2,1);
    stat->cd(1);
    nbEntries->Draw("");
    stat->cd(2);
    nbEntriesPerSp->SetStats(0);
    nbEntriesPerSp->Draw("");
    TCanvas *stat1 = new TCanvas("stat1","",50,50,600,800);
    stat1->cd();
    nbEntriesPerGroup->SetStats(0);
    nbEntriesPerGroup->Draw("");
  }

  return info;

}
//____________________________________________________________________________
Double_t *AliTRDCalibraFillHisto::GetMeanMedianRMSNumberCH()
{
  //
  // Return a Int_t[4] with:
  // 0 Mean number of entries
  // 1 median of number of entries
  // 2 rms of number of entries
  // 3 number of group with entries
  //

  Double_t *stat = new Double_t[4];
  stat[3]             = 0.0;

  Int_t    nbofgroups = CalculateTotalNumberOfBins(0);
  
  Double_t *weight = new Double_t[nbofgroups];
  Double_t *nonul = new Double_t[nbofgroups];
 
  for(Int_t k = 0; k < nbofgroups; k++){
    if(fEntriesCH[k] > 0) {
      weight[k] = 1.0;
      nonul[(Int_t)stat[3]] = fEntriesCH[k];
      stat[3]++;
    }
    else weight[k] = 0.0;
  }
  stat[0]          = TMath::Mean(nbofgroups,fEntriesCH,weight); 
  stat[1]          = TMath::Median(nbofgroups,fEntriesCH,weight); 
  stat[2]          = TMath::RMS((Int_t)stat[3],nonul); 

  delete [] weight;
  delete [] nonul;

  return stat;

}
//____________________________________________________________________________
Double_t *AliTRDCalibraFillHisto::GetMeanMedianRMSNumberLinearFitter() const
{
  //
  // Return a Int_t[4] with:
  // 0 Mean number of entries
  // 1 median of number of entries
  // 2 rms of number of entries
  // 3 number of group with entries
  //

  Double_t *stat      = new Double_t[4]; 
  stat[3]             = 0.0;

  Int_t    nbofgroups = 540;
  Double_t *weight    = new Double_t[nbofgroups];
  Int_t    *nonul     = new Int_t[nbofgroups]; 

  for(Int_t k = 0; k < nbofgroups; k++){
    if(fEntriesLinearFitter[k] > 0) {
      weight[k] = 1.0;
      nonul[(Int_t) stat[3]] = fEntriesLinearFitter[k];
      stat[3]++;     
    }
    else weight[k] = 0.0;
  }
  stat[0]          = TMath::Mean(nbofgroups,fEntriesLinearFitter,weight); 
  stat[1]          = TMath::Median(nbofgroups,fEntriesLinearFitter,weight); 
  stat[2]          = TMath::RMS((Int_t)stat[3],nonul); 

  delete [] weight;
  delete [] nonul;

  return stat;

}
//////////////////////////////////////////////////////////////////////////////////////
// Create Histos
//////////////////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::CreatePRF2d(Int_t nn)
{
  //
  // Create the 2D histos: here we have 2*fNgroupprf bins in tnp of 0.2 amplitude each
  // If fNgroupprf is zero then no binning in tnp
  //

  TString name("Nz");
  name += fCalibraMode->GetNz(2);
  name += "Nrphi";
  name += fCalibraMode->GetNrphi(2);
  name += "Ngp";
  name += fNgroupprf;

  if(fNgroupprf != 0){
    
    fPRF2d = new TProfile2D("PRF2d",(const Char_t *) name
			    ,2*fNgroupprf*fNumberBinPRF,-3.0*fNgroupprf,3.0*fNgroupprf,nn,0,nn);
    fPRF2d->SetYTitle("Det/pad groups");
    fPRF2d->SetXTitle("Position x/W [pad width units]");
    fPRF2d->SetZTitle("Q_{i}/Q_{total}");
    fPRF2d->SetStats(0);
  }
  else{
    fPRF2d = new TProfile2D("PRF2d",(const Char_t *) name
			    ,fNumberBinPRF,-1.5,1.5,nn,0,nn);
    fPRF2d->SetYTitle("Det/pad groups");
    fPRF2d->SetXTitle("Position x/W [pad width units]");
    fPRF2d->SetZTitle("Q_{i}/Q_{total}");
    fPRF2d->SetStats(0);
  }

}

//_____________________________________________________________________________
void AliTRDCalibraFillHisto::CreatePH2d(Int_t nn)
{
  //
  // Create the 2D histos
  //

  TString name("Ver");
  name += fVersionVdriftUsed;
  name += "Subver";
  name += fSubVersionVdriftUsed;
  name += "FirstRun";
  name += fFirstRunVdrift;
  name += "Nz";
  name += fCalibraMode->GetNz(1);
  name += "Nrphi";
  name += fCalibraMode->GetNrphi(1);
  
  fPH2d = new TProfile2D("PH2d",(const Char_t *) name
			 ,fTimeMax,-0.5/fSf,(Float_t) (fTimeMax-0.5)/fSf
			 ,nn,0,nn);
  fPH2d->SetYTitle("Det/pad groups");
  fPH2d->SetXTitle("time [#mus]");
  fPH2d->SetZTitle("<PH> [a.u.]");
  fPH2d->SetStats(0);

}
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::CreateCH2d(Int_t nn)
{
  //
  // Create the 2D histos
  //

  TString name("Ver");
  name += fVersionGainUsed;
  name += "Subver";
  name += fSubVersionGainUsed;
  name += "FirstRun";
  name += fFirstRunGain;
  name += "Nz";
  name += fCalibraMode->GetNz(0);
  name += "Nrphi";
  name += fCalibraMode->GetNrphi(0);
  
  fCH2d = new TH2I("CH2d",(const Char_t *) name
		   ,fNumberBinCharge,0,300,nn,0,nn);
  fCH2d->SetYTitle("Det/pad groups");
  fCH2d->SetXTitle("charge deposit [a.u]");
  fCH2d->SetZTitle("counts");
  fCH2d->SetStats(0);
  fCH2d->Sumw2();

}
//////////////////////////////////////////////////////////////////////////////////
// Set relative scale
/////////////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::SetRelativeScale(Float_t RelativeScale)
{
  //
  // Set the factor that will divide the deposited charge
  // to fit in the histo range [0,300]
  //
 
  if (RelativeScale > 0.0) {
    fRelativeScale = RelativeScale;
  } 
  else {
    AliInfo("RelativeScale must be strict positif!");
  }

}
//////////////////////////////////////////////////////////////////////////////////
// Quick way to fill a histo
//////////////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________
void  AliTRDCalibraFillHisto::FillCH2d(Int_t x, Float_t y)
{
  //
  // FillCH2d: Marian style
  // 
  
  //skip simply the value out of range
  if((y>=300.0) || (y<0.0)) return;
  
  //Calcul the y place
  Int_t yplace = (Int_t) (fNumberBinCharge*y/300.0)+1;
  Int_t place = (fNumberBinCharge+2)*(x+1)+yplace;
  
  //Fill
  fCH2d->GetArray()[place]++;

}
 
//////////////////////////////////////////////////////////////////////////////////
// Geometrical functions
///////////////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________________
Int_t AliTRDCalibraFillHisto::GetLayer(Int_t d) const
{
  //
  // Reconstruct the layer number from the detector number
  //

  return ((Int_t) (d % 6));

}

//_____________________________________________________________________________
Int_t AliTRDCalibraFillHisto::GetStack(Int_t d) const
{
  //
  // Reconstruct the stack number from the detector number
  //
  const Int_t kNlayer = 6;

  return ((Int_t) (d % 30) / kNlayer);

}

//_____________________________________________________________________________
Int_t AliTRDCalibraFillHisto::GetSector(Int_t d) const
{
  //
  // Reconstruct the sector number from the detector number
  //
  Int_t fg = 30;

  return ((Int_t) (d / fg));

}
///////////////////////////////////////////////////////////////////////////////////
// Getter functions for DAQ of the CH2d and the PH2d
//////////////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________
TProfile2D* AliTRDCalibraFillHisto::GetPH2d(Int_t nbtimebin, Float_t samplefrequency)
{
    //
    // return pointer to fPH2d TProfile2D
    // create a new TProfile2D if it doesn't exist allready
    //
    if ( fPH2d )
	return fPH2d;

    // Some parameters
    fTimeMax = nbtimebin;
    fSf      = samplefrequency;
  
    CreatePH2d(540);

    return fPH2d;
}
//_____________________________________________________________________
TH2I* AliTRDCalibraFillHisto::GetCH2d()
{
    //
    // return pointer to fCH2d TH2I
    // create a new TH2I if it doesn't exist allready
    //
    if ( fCH2d )
        return fCH2d;

    CreateCH2d(540);

    return fCH2d;
}
////////////////////////////////////////////////////////////////////////////////////////////
// Drift velocity calibration
///////////////////////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________
TLinearFitter* AliTRDCalibraFillHisto::GetLinearFitter(Int_t detector, Bool_t force)
{
    //
    // return pointer to TLinearFitter Calibration
    // if force is true create a new TLinearFitter if it doesn't exist allready
    //

  if ((!force) || (fLinearFitterArray.UncheckedAt(detector))){
    return (TLinearFitter*)fLinearFitterArray.UncheckedAt(detector);
  }

  // if we are forced and TLinearFitter doesn't yet exist create it

  // new TLinearFitter
  TLinearFitter *linearfitter = new TLinearFitter(2,"pol1");
  fLinearFitterArray.AddAt(linearfitter,detector);
  return linearfitter;
}

//____________________________________________________________________________
void AliTRDCalibraFillHisto::AnalyseLinearFitter()
{
  //
  // Analyse array of linear fitter because can not be written
  // Store two arrays: one with the param the other one with the error param + number of entries
  //

  for(Int_t k = 0; k < 540; k++){
    TLinearFitter *linearfitter = GetLinearFitter(k);
    if((linearfitter!=0) && (fEntriesLinearFitter[k]>10)){
      TVectorD  *par  = new TVectorD(2);
      TVectorD   pare = TVectorD(2);
      TVectorD  *parE = new TVectorD(3);
      if((linearfitter->EvalRobust(0.8)==0)) {
	//linearfitter->Eval();
	linearfitter->GetParameters(*par);
	//linearfitter->GetErrors(pare);
	//Float_t  ppointError =  TMath::Sqrt(TMath::Abs(linearfitter->GetChisquare())/fEntriesLinearFitter[k]);
	//(*parE)[0] = pare[0]*ppointError;
	//(*parE)[1] = pare[1]*ppointError;

	(*parE)[0] = 0.0;
	(*parE)[1] = 0.0;
	(*parE)[2] = (Double_t) fEntriesLinearFitter[k];
	((TObjArray *)fLinearVdriftFit->GetPArray())->AddAt(par,k);
	((TObjArray *)fLinearVdriftFit->GetEArray())->AddAt(parE,k);
      }
    }
  }
}



