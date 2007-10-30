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
// Author:
//   R. Bailhache (R.Bailhache@gsi.de)
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

#include "AliLog.h"

#include "AliTRDCalibraFillHisto.h"
#include "AliTRDCalibraMode.h"
#include "AliTRDCalibraVector.h"
#include "AliTRDCalibraVdriftLinearFit.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "AliTRDmcmTracklet.h"
#include "AliTRDpadPlane.h"
#include "AliTRDcluster.h"
#include "AliTRDtrack.h"
#include "AliTRDRawStreamV2.h"
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliTRDgeometry.h"
#include "./Cal/AliTRDCalROC.h"
#include "./Cal/AliTRDCalDet.h"

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
  ,fMITracking(kFALSE)
  ,fMcmTracking(kFALSE)
  ,fMcmCorrectAngle(kFALSE)
  ,fCH2dOn(kFALSE)
  ,fPH2dOn(kFALSE)
  ,fPRF2dOn(kFALSE)
  ,fHisto2d(kFALSE)
  ,fVector2d(kFALSE)
  ,fLinearFitterOn(kFALSE)
  ,fLinearFitterDebugOn(kFALSE)
  ,fRelativeScale(0)
  ,fThresholdClusterPRF2(15.0)
  ,fCalibraMode(new AliTRDCalibraMode())
  ,fDebugStreamer(0)
  ,fDebugLevel(0)
  ,fDetectorAliTRDtrack(kFALSE)
  ,fDetectorPreviousTrack(-1)
  ,fMCMPrevious(-1)
  ,fROBPrevious(-1)
  ,fNumberClusters(18)
  ,fProcent(6.0)
  ,fDifference(17)
  ,fNumberTrack(0)
  ,fTimeMax(0)
  ,fSf(10.0)
  ,fNumberBinCharge(100)
  ,fNumberBinPRF(40)
  ,fNgroupprf(0)
  ,fListClusters(new TObjArray()) 
  ,fPar0(0x0)
  ,fPar1(0x0)
  ,fPar2(0x0)
  ,fPar3(0x0)
  ,fPar4(0x0)
  ,fAmpTotal(0x0)
  ,fPHPlace(0x0)
  ,fPHValue(0x0)
  ,fGoodTracklet(kTRUE)
  ,fGoodTrack(kTRUE)
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
  ,fCalDetT0(0x0)
  ,fCalROCT0(0x0)
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

}

//______________________________________________________________________________________
AliTRDCalibraFillHisto::AliTRDCalibraFillHisto(const AliTRDCalibraFillHisto &c)
  :TObject(c)
  ,fGeo(0)
  ,fMITracking(c.fMITracking)
  ,fMcmTracking(c.fMcmTracking)
  ,fMcmCorrectAngle(c.fMcmCorrectAngle)
  ,fCH2dOn(c.fCH2dOn)
  ,fPH2dOn(c.fPH2dOn)
  ,fPRF2dOn(c.fPRF2dOn)
  ,fHisto2d(c.fHisto2d)
  ,fVector2d(c.fVector2d)
  ,fLinearFitterOn(c.fLinearFitterOn)
  ,fLinearFitterDebugOn(c.fLinearFitterDebugOn)
  ,fRelativeScale(c.fRelativeScale)
  ,fThresholdClusterPRF2(c.fThresholdClusterPRF2)
  ,fCalibraMode(0x0)
  ,fDebugStreamer(0)
  ,fDebugLevel(c.fDebugLevel)
  ,fDetectorAliTRDtrack(c.fDetectorAliTRDtrack)
  ,fDetectorPreviousTrack(c.fDetectorPreviousTrack)
  ,fMCMPrevious(c.fMCMPrevious)
  ,fROBPrevious(c.fROBPrevious)
  ,fNumberClusters(c.fNumberClusters)
  ,fProcent(c.fProcent)
  ,fDifference(c.fDifference)
  ,fNumberTrack(c.fNumberTrack)
  ,fTimeMax(c.fTimeMax)
  ,fSf(c.fSf)
  ,fNumberBinCharge(c.fNumberBinCharge)
  ,fNumberBinPRF(c.fNumberBinPRF)
  ,fNgroupprf(c.fNgroupprf)
  ,fListClusters(new TObjArray())
  ,fPar0(0x0)
  ,fPar1(0x0)
  ,fPar2(0x0)
  ,fPar3(0x0)
  ,fPar4(0x0)
  ,fAmpTotal(0x0)
  ,fPHPlace(0x0)
  ,fPHValue(0x0)
  ,fGoodTracklet(c.fGoodTracklet)
  ,fGoodTrack(c.fGoodTrack)
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
  ,fCalDetT0(0x0)
  ,fCalROCT0(0x0)
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
  if(c.fCalDetT0)    fCalDetT0     = new AliTRDCalDet(*c.fCalDetT0);
  if(c.fCalROCGain)  fCalROCGain   = new AliTRDCalROC(*c.fCalROCGain);
  if(c.fCalROCT0)    fCalROCT0     = new AliTRDCalROC(*c.fCalROCT0);

  if (fGeo) {
    delete fGeo;
  }
  fGeo = new AliTRDgeometry();
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
  if ( fCalDetT0 )    delete fCalDetT0;
  if ( fCalROCGain )  delete fCalROCGain;
  if ( fCalROCT0 )    delete fCalROCT0; 

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

//____________Functions for initialising the AliTRDCalibraFillHisto in the code_________
Bool_t AliTRDCalibraFillHisto::Init2Dhistos()
{
  //
  // For the offline tracking
  // This function will be called in the function AliReconstruction::Run() 
  // Init the calibration mode (Nz, Nrphi), the 2D histograms if fHisto2d = kTRUE, 
  //

  Init2Dhistostrack();
  
  //Init the tracklet parameters
  fPar0 = new Double_t[fTimeMax];
  fPar1 = new Double_t[fTimeMax];
  fPar2 = new Double_t[fTimeMax];
  fPar3 = new Double_t[fTimeMax];
  fPar4 = new Double_t[fTimeMax];
  
  for(Int_t k = 0; k < fTimeMax; k++){
    fPar0[k] = 0.0;
    fPar1[k] = 0.0;
    fPar2[k] = 0.0;
    fPar3[k] = 0.0;
    fPar4[k] = 0.0;
  }
  return kTRUE;
}

//____________Functions for initialising the AliTRDCalibraFillHisto in the code_________
Bool_t AliTRDCalibraFillHisto::Init2Dhistostrack()
{
  //
  // For the offline tracking
  // This function will be called in the function AliReconstruction::Run() 
  // Init the calibration mode (Nz, Nrphi), the 2D histograms if fHisto2d = kTRUE, 
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
  fTimeMax = cal->GetNumberOfTimeBins();
  fSf      = parCom->GetSamplingFrequency();
  fRelativeScale = 20;
 
  //calib object from database used for reconstruction
  if(fCalDetGain) delete fCalDetGain;
  fCalDetGain  = new AliTRDCalDet(*(cal->GetGainFactorDet()));
  if(fCalDetT0)   delete fCalDetT0;
  fCalDetT0   = new AliTRDCalDet(*(cal->GetT0Det()));

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
    TString namech("Nz");
    namech += fCalibraMode->GetNz(0);
    namech += "Nrphi";
    namech += fCalibraMode->GetNrphi(0);
    fCalibraVector->SetNameCH((const char* ) namech);
    TString nameph("Nz");
    nameph += fCalibraMode->GetNz(1);
    nameph += "Nrphi";
    nameph += fCalibraMode->GetNrphi(1);
    fCalibraVector->SetNamePH((const char* ) nameph);
    TString nameprf("Nz");
    nameprf += fCalibraMode->GetNz(2);
    nameprf += "Nrphi";
    nameprf += fCalibraMode->GetNrphi(2);
    nameprf += "Ngp";
    nameprf += fNgroupprf;
    fCalibraVector->SetNamePRF((const char* ) nameprf);
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
    //fLinearFitterArray.Expand(540);
    fLinearFitterArray.SetName("ArrayLinearFitters");
    fEntriesLinearFitter = new Int_t[540];
    for(Int_t k = 0; k < 540; k++){
      fEntriesLinearFitter[k] = 0;
    }
    fLinearVdriftFit = new AliTRDCalibraVdriftLinearFit();
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

//____________Functions for filling the histos in the code_____________________

//____________Offine tracking in the AliTRDtracker_____________________________
Bool_t AliTRDCalibraFillHisto::ResetTrack()
{
  //
  // For the offline tracking
  // This function will be called in the function
  // AliTRDtracker::FollowBackPropagation() at the beginning 
  // Reset the parameter to know we have a new TRD track
  //
  
  fDetectorAliTRDtrack = kFALSE;
  //fGoodTrack           = kTRUE;
  return kTRUE;

}
//____________Offine tracking in the AliTRDtracker_____________________________
void AliTRDCalibraFillHisto::ResetfVariables()
{
  //
  // Reset values per tracklet
  //

  // Reset the list of clusters
  fListClusters->Clear();

  //Reset the tracklet parameters
  for(Int_t k = 0; k < fTimeMax; k++){
    fPar0[k] = 0.0;
    fPar1[k] = 0.0;
    fPar2[k] = 0.0;
    fPar3[k] = 0.0; 
    fPar4[k] = 0.0;
  }

  ResetfVariablestrack();
  
}
//____________Offine tracking in the AliTRDtracker_____________________________
void AliTRDCalibraFillHisto::ResetfVariablestrack()
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
//____________Offline tracking in the AliTRDtracker____________________________
Bool_t AliTRDCalibraFillHisto::UpdateHistograms(AliTRDtrack *t)
{
  //
  // For the offline tracking
  // This function will be called in the function
  // AliTRDtracker::FollowBackPropagation() in the loop over the clusters
  // of TRD tracks 
  // Fill the 2D histos or the vectors with the info of the clusters at
  // the end of a detectors if the track is "good"
  //


  
  AliTRDcluster *cl = 0x0;
  Int_t index0 = 0;
  Int_t index1 = 0;
  
  // reset if good track
  fGoodTrack = kTRUE;

  
  // loop over the clusters
  while((cl = t->GetCluster(index1))){

    // Localisation of the detector
    Int_t detector = cl->GetDetector();

   
    // Fill the infos for the previous clusters if not the same
    // detector anymore but this time it should be the same track
    if ((detector != fDetectorPreviousTrack) && 
	(index0 != index1)) {
      
      fNumberTrack++;   
         
      //printf("detector %d, fPreviousdetector %d, plane %d, planeprevious %d, index0 %d, index1 %d la\n",detector,fDetectorPreviousTrack,GetPlane(detector),GetPlane(fDetectorPreviousTrack),index0,index1);

      //If the same track, then look if the previous detector is in
      //the same plane, if yes: not a good track
      //FollowBack
      //if (fDetectorAliTRDtrack && 
      // (GetPlane(detector) <= GetPlane(fDetectorPreviousTrack))) {
      //Follow
      if ((GetPlane(detector) >= GetPlane(fDetectorPreviousTrack))) {
	fGoodTrack = kFALSE;
      }
      
      // Fill only if the track doesn't touch a masked pad or doesn't
      // appear in the middle (fGoodTrack)
      if (fGoodTrack && fGoodTracklet) {
	
	// drift velocity unables to cut bad tracklets 
	Bool_t  pass = FindP1TrackPHtrack(t,index0,index1);
	
	// Gain calibration
	if (fCH2dOn) {
	  FillTheInfoOfTheTrackCH();
	}
	
	// PH calibration
	if (fPH2dOn) {
	  FillTheInfoOfTheTrackPH();    
	}
	
	if(pass && fPRF2dOn) HandlePRFtrack(t,index0,index1);
	
	
      } // if a good track
 
      // reset stuff     
      ResetfVariablestrack();
      index0 = index1;
   
    } // Fill at the end the charge
    
    // Calcul the position of the detector and take the calib objects
    if (detector != fDetectorPreviousTrack) {
      
      //Localise the detector
      LocalisationDetectorXbins(detector);
      
      // Get cal
      AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
      if (!cal) {
	AliInfo("Could not get calibDB");
	return kFALSE;
      }
      
      // Get calib objects
      if( fCalROCGain ) delete fCalROCGain;
      fCalROCGain = new AliTRDCalROC(*(cal->GetGainFactorROC(detector)));
      if( fCalROCT0 )   delete fCalROCT0;
      fCalROCT0   = new AliTRDCalROC(*(cal->GetT0ROC(detector)));
      
    }
    
    // Reset the detectbjobsor
    fDetectorPreviousTrack = detector;

    // Store the info bis of the tracklet
    Int_t *rowcol   = CalculateRowCol(cl);
    CheckGoodTracklet(detector,rowcol);
    Int_t     group[2] = {0,0};
    if(fCH2dOn)  group[0]  = CalculateCalibrationGroup(0,rowcol);
    if(fPH2dOn)  group[1]  = CalculateCalibrationGroup(1,rowcol);
    StoreInfoCHPHtrack(cl,t,index1,group,rowcol);
         
    index1++;

  } // while on clusters

  // Fill the last plane
  if( index0 != index1 ){

    //printf("fPreviousdetector %d, planeprevious %d, index0 %d, index1 %d li\n",fDetectorPreviousTrack,GetPlane(fDetectorPreviousTrack),index0,index1);
    
    fNumberTrack++; 
    
    if (fGoodTrack && fGoodTracklet) {
      
      // drift velocity unables to cut bad tracklets 
      Bool_t  pass = FindP1TrackPHtrack(t,index0,index1);
      
      // Gain calibration
      if (fCH2dOn) {
	FillTheInfoOfTheTrackCH();
      }
      
      // PH calibration
      if (fPH2dOn) {
	FillTheInfoOfTheTrackPH();    
      }
      
      if(pass && fPRF2dOn) HandlePRFtrack(t,index0,index1);
          
    } // if a good track
    
  }

  // reset stuff     
  ResetfVariablestrack();
   
  return kTRUE;
  
}
//____________Offline tracking in the AliTRDtracker____________________________
Bool_t AliTRDCalibraFillHisto::UpdateHistograms(AliTRDcluster *cl, AliTRDtrack *t)
{
  //
  // For the offline tracking
  // This function will be called in the function
  // AliTRDtracker::FollowBackPropagation() in the loop over the clusters
  // of TRD tracks 
  // Fill the 2D histos or the vectors with the info of the clusters at
  // the end of a detectors if the track is "good"
  //

  // Localisation of the detector
  Int_t detector = cl->GetDetector();
 
  // Fill the infos for the previous clusters if not the same
  // detector anymore or if not the same track
  if (((detector != fDetectorPreviousTrack) || (!fDetectorAliTRDtrack)) && 
      (fDetectorPreviousTrack != -1)) {
    
    fNumberTrack++;   
    
    // If the same track, then look if the previous detector is in
    // the same plane, if yes: not a good track
    //FollowBack
    if (fDetectorAliTRDtrack && 
    	(GetPlane(detector) <= GetPlane(fDetectorPreviousTrack))) {
    //Follow
    //if (fDetectorAliTRDtrack && 
    //    (GetPlane(detector) >= GetPlane(fDetectorPreviousTrack))) {
      fGoodTrack = kFALSE;
    }

    // Fill only if the track doesn't touch a masked pad or doesn't
    // appear in the middle (fGoodTrack)
    if (fGoodTrack && fGoodTracklet) {

      // drift velocity unables to cut bad tracklets 
      Bool_t  pass = FindP1TrackPH();
      
      // Gain calibration
      if (fCH2dOn) {
	FillTheInfoOfTheTrackCH();
      }
      
      // PH calibration
      if (fPH2dOn) {
	FillTheInfoOfTheTrackPH();    
      }

      if(pass && fPRF2dOn) HandlePRF();

      
    } // if a good track
    
    ResetfVariables();
    if(!fDetectorAliTRDtrack) fGoodTrack = kTRUE;
    
  } // Fill at the end the charge
  
  // Calcul the position of the detector and take the calib objects
  if (detector != fDetectorPreviousTrack) {
    //Localise the detector
    LocalisationDetectorXbins(detector);
    
    // Get cal
    AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
    if (!cal) {
      AliInfo("Could not get calibDB");
      return kFALSE;
    }
    
    // Get calib objects
    if( fCalROCGain ) delete fCalROCGain;
    fCalROCGain = new AliTRDCalROC(*(cal->GetGainFactorROC(detector)));
    if( fCalROCT0 )   delete fCalROCT0;
    fCalROCT0   = new AliTRDCalROC(*(cal->GetT0ROC(detector)));
  }
 
  // Reset the detector
  fDetectorPreviousTrack = detector;
  fDetectorAliTRDtrack   = kTRUE;

  // Store the infos of the tracklets
  AliTRDcluster *kcl = new AliTRDcluster(*cl);
  fListClusters->Add((TObject *)kcl);
  Int_t time = cl->GetLocalTimeBin();
  fPar0[time] = t->GetY();
  fPar1[time] = t->GetZ();
  fPar2[time] = t->GetSnp();
  fPar3[time] = t->GetTgl();
  fPar4[time] = t->GetSigned1Pt();

  // Store the info bis of the tracklet
  Int_t *rowcol   = CalculateRowCol(cl);
  CheckGoodTracklet(detector,rowcol);
  Int_t     group[2] = {0,0};
  if(fCH2dOn)  group[0]  = CalculateCalibrationGroup(0,rowcol);
  if(fPH2dOn)  group[1]  = CalculateCalibrationGroup(1,rowcol);
  StoreInfoCHPH(cl,t,group,rowcol);
   
  return kTRUE;
  
}
//____________Online trackling in AliTRDtrigger________________________________
Bool_t AliTRDCalibraFillHisto::UpdateHistogramcm(AliTRDmcmTracklet *trk)
{
  //
  // For the tracking
  // This function will be called in the function AliTRDtrigger::TestTracklet
  // before applying the pt cut on the tracklets 
  // Fill the infos for the tracklets fTrkTest if the tracklets is "good"
  //
  
  // Localisation of the Xbins involved
  Int_t idect = trk->GetDetector();
  fDetectorPreviousTrack = idect;
  LocalisationDetectorXbins(idect);

  // Get the parameter object
  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }
   
  // Reset
  ResetfVariables();

  // Get calib objects
  if( fCalROCGain ) delete fCalROCGain;
  fCalROCGain = new AliTRDCalROC(*(cal->GetGainFactorROC(idect)));
  if( fCalROCT0 )   delete fCalROCT0;
  fCalROCT0   = new AliTRDCalROC(*(cal->GetT0ROC(idect)));
   
  // Row of the tracklet and position in the pad groups
  Int_t *rowcol  = new Int_t[2];
  rowcol[0]     = trk->GetRow();
  Int_t group[3] = {-1,-1,-1};
  
  // Eventuelle correction due to track angle in z direction
  Float_t correction = 1.0;
  if (fMcmCorrectAngle) {
    Float_t z = trk->GetRowz();
    Float_t r = trk->GetTime0();
    correction = r / TMath::Sqrt((r*r+z*z));
  }

  // Boucle sur les clusters
  // Condition on number of cluster: don't come from the middle of the detector
  if (trk->GetNclusters() >= fNumberClusters) {

    for (Int_t icl = 0; icl < trk->GetNclusters(); icl++) {

      Float_t amp[3] = { 0.0, 0.0, 0.0 };
      Int_t   time   = trk->GetClusterTime(icl);
      rowcol[1]      = trk->GetClusterCol(icl);
            
      amp[0] = trk->GetClusterADC(icl)[0] * correction;
      amp[1] = trk->GetClusterADC(icl)[1] * correction;
      amp[2] = trk->GetClusterADC(icl)[2] * correction;

      
      if ((amp[0] < 0.0) || 
          (amp[1] < 0.0) || 
          (amp[2] < 0.0)) {
        continue;
      }

      // Col of cluster and position in the pad groups
      if(fCH2dOn)  {
	group[0] = CalculateCalibrationGroup(0,rowcol);
	fAmpTotal[(Int_t) group[0]] += (Float_t) (amp[0]+amp[1]+amp[2]);
      }
      if(fPH2dOn)  {
	group[1] = CalculateCalibrationGroup(1,rowcol);
	fPHPlace[time] = group[1];
	fPHValue[time] = (Float_t) (amp[0]+amp[1]+amp[2]);
      }
      if(fPRF2dOn) group[2] = CalculateCalibrationGroup(2,rowcol);

      // See if we are not near a masked pad fGoodTracklet
      CheckGoodTracklet(idect,rowcol);
               
      // Fill PRF direct without tnp bins...only for monitoring...
      if (fPRF2dOn && fGoodTracklet) {
	
	if ((amp[0] > fThresholdClusterPRF2) && 
            (amp[1] > fThresholdClusterPRF2) && 
            (amp[2] > fThresholdClusterPRF2) && 
            ((amp[0]*amp[2]/(amp[1]*amp[1])) < 0.06)) {
	
	  // Security of the denomiateur is 0
	  if ((((Float_t) (((Float_t) amp[1]) * ((Float_t) amp[1]))) 
             / ((Float_t) (((Float_t) amp[0]) * ((Float_t) amp[2])))) != 1.0) {
	    Float_t xcenter = 0.5 * (TMath::Log(amp[2] / amp[0]))
                                  / (TMath::Log((amp[1]*amp[1]) / (amp[0]*amp[2])));
	    Float_t ycenter = amp[1] / (amp[0] + amp[1] + amp[2]);

	    if (TMath::Abs(xcenter) < 0.5) {
	      Float_t yminus = amp[0] / (amp[0]+amp[1]+amp[2]);
	      Float_t ymax   = amp[2] / (amp[0]+amp[1]+amp[2]);
	      // Fill only if it is in the drift region!
	      if (((Float_t) time / fSf) > 0.3) {
		if (fHisto2d) {
		  fPRF2d->Fill(xcenter,(fCalibraMode->GetXbins(2)+group[2]+0.5),ycenter);
		  fPRF2d->Fill(-(xcenter+1.0),(fCalibraMode->GetXbins(2)+group[2]+0.5),yminus);
		  fPRF2d->Fill((1.0-xcenter),(fCalibraMode->GetXbins(2)+group[2]+0.5),ymax);
		}
		if (fVector2d) {
		  fCalibraVector->UpdateVectorPRF(idect,group[2],xcenter,ycenter);
		  fCalibraVector->UpdateVectorPRF(idect,group[2],-(xcenter+1.0),yminus);
		  fCalibraVector->UpdateVectorPRF(idect,group[2],(1.0-xcenter),ymax);
		}
	      }//in the drift region 
	    }//in the middle
	  }//denominateur security
	}//cluster shape and thresholds
      }//good and PRF On
      
    } // Boucle clusters
    
    // Fill the charge
    if(fGoodTracklet){
      if (fCH2dOn) FillTheInfoOfTheTrackCH();
      if (fPH2dOn) FillTheInfoOfTheTrackPH();	
    }

    fNumberTrack++;
        
  } // Condition on number of clusters

  return kTRUE;
  
}
//_____________________________________________________________________________
Int_t *AliTRDCalibraFillHisto::CalculateRowCol(AliTRDcluster *cl) const
{
  //
  // Calculate the row and col number of the cluster
  //


  Int_t *rowcol = new Int_t[2];
  rowcol[0] =  0;
  rowcol[1] =  0;

  // Localisation of the detector
  Int_t detector = cl->GetDetector();
  Int_t chamber  = GetChamber(detector);
  Int_t plane    = GetPlane(detector);

  // Localisation of the cluster
  Double_t pos[3] = { 0.0, 0.0, 0.0 };
  pos[0] = ((AliCluster *)cl)->GetX();
  pos[1] = cl->GetY();
  pos[2] = cl->GetZ();

  // Position of the cluster
  AliTRDpadPlane *padplane  = fGeo->GetPadPlane(plane,chamber);
  Int_t    row              = padplane->GetPadRowNumber(pos[2]);
  //Do not take from here because it was corrected from ExB already....
  //Double_t offsetz         = padplane->GetPadRowOffset(row,pos[2]);
  //Double_t offsettilt      = padplane->GetTiltOffset(offsetz);
  //Int_t    col             = padplane->GetPadColNumber(pos[1] + offsettilt,offsetz);
  //Int_t    col             = padplane->GetPadColNumber(pos[1]+offsettilt);
  Int_t    col               = cl->GetPadCol(); 

  //return
  rowcol[0]     = row;
  rowcol[1]     = col; 
  return rowcol;
  
}
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::CheckGoodTracklet(Int_t detector, Int_t *rowcol)
{
  //
  // See if we are not near a masked pad
  //

  Int_t row = rowcol[0];
  Int_t col = rowcol[1];

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
Bool_t AliTRDCalibraFillHisto::IsPadOn(Int_t detector, Int_t col, Int_t row) const
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
  
  if (!cal->IsChamberInstalled(detector)     || 
       cal->IsChamberMasked(detector)        ||
       cal->IsPadMasked(detector,col,row)) {
    return kFALSE;
  }
  else {
    return kTRUE;
  }
  
}
//_____________________________________________________________________________
Int_t AliTRDCalibraFillHisto::CalculateCalibrationGroup(Int_t i, Int_t *rowcol) const
{
  //
  // Calculate the calibration group number for i
  //
 
  // Row of the cluster and position in the pad groups
  Int_t posr = 0;
  if (fCalibraMode->GetNnZ(i) != 0) {
    posr = (Int_t) rowcol[0] / fCalibraMode->GetNnZ(i);
  }
 
      
  // Col of the cluster and position in the pad groups
  Int_t posc = 0;
  if (fCalibraMode->GetNnRphi(i) != 0) {
    posc = (Int_t) rowcol[1] / fCalibraMode->GetNnRphi(i);
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
    fCalibraMode->ModePadCalibration((Int_t) GetChamber(detector),i);
    fCalibraMode->ModePadFragmentation((Int_t) GetPlane(detector)
                       , (Int_t) GetChamber(detector)
                       , (Int_t) GetSector(detector),i);
  }
  
  return kTRUE;

}
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::StoreInfoCHPH(AliTRDcluster *cl, AliTRDtrack *t, Int_t *group, Int_t *rowcol)
{
  //
  // Store the infos in fAmpTotal, fPHPlace and fPHValue
  //
  
  // Charge in the cluster
  Float_t  q        = TMath::Abs(cl->GetQ());
  Int_t    time     = cl->GetLocalTimeBin();

  //Correct for the gain coefficient used in the database for reconstruction
  Float_t correctthegain = fCalDetGain->GetValue(fDetectorPreviousTrack)*fCalROCGain->GetValue(rowcol[1],rowcol[0]);
  Float_t correcttheT0   = fCalDetT0->GetValue(fDetectorPreviousTrack)+fCalROCT0->GetValue(rowcol[1],rowcol[0]);
    
  // we substract correcttheT0 in AliTRDclusterizerV1::MakeClusters (line 458)
  Int_t timec            = Arrondi((Double_t)(time+correcttheT0));
  if((correcttheT0+0.5)==(int(correcttheT0+0.5))) {
    timec++;
  }
  if( timec < 0 ) return;


  // Correction due to the track angle
  Float_t correction    = 1.0;
  Float_t normalisation = 6.67;
  // we divide with gain in AliTRDclusterizerV1::Transform...
  if( correctthegain > 0 ) normalisation /= correctthegain;
  if ((q >0) && (t->GetNdedx() > 0)) {
    correction = t->GetClusterdQdl((t->GetNdedx() - 1)) / (normalisation);
  }

  // Fill the fAmpTotal with the charge
  if (fCH2dOn) {
    fAmpTotal[(Int_t) group[0]] += correction;
  }

  // Fill the fPHPlace and value
  if (fPH2dOn) {
    fPHPlace[timec] = group[1];
    fPHValue[timec] = correction;
  }
  
}
//_____________________________________________________________________________
void AliTRDCalibraFillHisto::StoreInfoCHPHtrack(AliTRDcluster *cl, AliTRDtrack *t, Int_t index, Int_t *group, Int_t *rowcol)
{
  //
  // Store the infos in fAmpTotal, fPHPlace and fPHValue
  //
  
  // Charge in the cluster
  Float_t  q        = TMath::Abs(cl->GetQ());
  Int_t    time     = cl->GetLocalTimeBin();
   
  //Correct for the gain coefficient used in the database for reconstruction
  Float_t correctthegain = fCalDetGain->GetValue(fDetectorPreviousTrack)*fCalROCGain->GetValue(rowcol[1],rowcol[0]);
  Float_t correcttheT0   = fCalDetT0->GetValue(fDetectorPreviousTrack)+fCalROCT0->GetValue(rowcol[1],rowcol[0]);
    
  // we substract correcttheT0 in AliTRDclusterizerV1::MakeClusters (line 458)
  Int_t timec            = Arrondi((Double_t)(time+correcttheT0));
  if((correcttheT0+0.5)==(int(correcttheT0+0.5))) {
    timec++;
  }
  if( timec < 0 ) return;

  // Correction due to the track angle
  Float_t correction    = 1.0;
  Float_t normalisation = 6.67;
  // we divide with gain in AliTRDclusterizerV1::Transform...
  if( correctthegain > 0 ) normalisation /= correctthegain;
  if (q >0) {
    correction = t->GetClusterdQdl(index) / (normalisation);
  }

  // Fill the fAmpTotal with the charge
  if (fCH2dOn) {
    fAmpTotal[(Int_t) group[0]] += correction;
  }

  // Fill the fPHPlace and value
  if (fPH2dOn) {
    fPHPlace[timec] = group[1];
    fPHValue[timec] = correction;
  }
  
}
//_____________________________________________________________________
Int_t AliTRDCalibraFillHisto::ProcessEventDAQ(AliTRDRawStreamV2 *rawStream, Bool_t nocheck)
{
  //
  // Event Processing loop - AliTRDRawStreamV2
  // 0 timebin problem
  // 1 no input
  // 2 input
  //
  
  Int_t withInput = 1;
  
  Int_t phvalue[21][36];
  for(Int_t k = 0; k < 36; k++){
    for(Int_t j = 0; j < 21; j++){
      phvalue[j][k] = 10;
    }
  }
  
  fDetectorPreviousTrack = -1;
  fMCMPrevious           = -1;
  fROBPrevious           = -1;
  Int_t nbtimebin = 0;                                        
  Int_t baseline  = 10;  

  // For selecting the signal
  Double_t mean[21]   = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  Int_t first[21]     = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 
  Int_t    select[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  if(!nocheck){
  
    fTimeMax = 0;
       
    while (rawStream->Next()) {
      
      Int_t idetector = rawStream->GetDet();                            //  current detector
      Int_t imcm      = rawStream->GetMCM();                            //  current MCM
      Int_t irob      = rawStream->GetROB();                            //  current ROB
      
      if(((fMCMPrevious != imcm) || (fDetectorPreviousTrack != idetector) || (fROBPrevious != irob)) && (fDetectorPreviousTrack != -1)){
	
	// take the mean values and check the first time bin
	for(Int_t j = 0; j < 21; j++){       
	  if(TMath::RMS(fTimeMax,phvalue[j]) != 0.0) mean[j] = TMath::Mean(fTimeMax,phvalue[j]);
	  else mean[j] = 0.0;
	  if(phvalue[j][0] > 200.0) first[j] = 1;
	  else first[j] = 0;
	}
	
	// select
	for(Int_t j = 1; j < 20; j++){
	  if((first[j-1] == 0) && (first[j] ==0) && (first[j+1] == 0) && (mean[j-1] > (baseline+5.0)) && (mean[j] > (baseline+10.0)) && (mean[j+1] > (baseline+5.0)) && (mean[j] >= mean[j-1]) && (mean[j] >= mean[j+1])){
	    select[j] = 1;
	  }
	  else select[j] = 0;
	}

	// fill
	for(Int_t j = 1; j < 20; j++){
	  if(select[j] == 1){
	    withInput = 2;
	    for(Int_t k = 0; k < fTimeMax; k++){
	      if((phvalue[j][k] >= phvalue[j-1][k]) && (phvalue[j][k] >= phvalue[j+1][k])){
		UpdateDAQ(fDetectorPreviousTrack,0,0,k,(phvalue[j-1][k]+phvalue[j][k]+phvalue[j+1][k]),fTimeMax);
	      }
	      else{
		if((j < 19) && (phvalue[j+1][k] >= phvalue[j][k]) && (phvalue[j+1][k] >= phvalue[j+2][k])){
		  UpdateDAQ(fDetectorPreviousTrack,0,0,k,(phvalue[j][k]+phvalue[j+1][k]+phvalue[j+2][k]),fTimeMax);
		}
		else UpdateDAQ(fDetectorPreviousTrack,0,0,k,(3*baseline),fTimeMax);
	      }	   
	    }
	  }
	}

	// reset
	for(Int_t k = 0; k < 36; k++){
	  for(Int_t j = 0; j < 21; j++){
	    phvalue[j][k] = baseline;
	  }
	}
      }

      fDetectorPreviousTrack = idetector;
      fMCMPrevious           = imcm;
      fROBPrevious           = irob;

      nbtimebin         = rawStream->GetNumberOfTimeBins();              //  number of time bins read from data
      if(nbtimebin == 0) return 0;
      if((fTimeMax != 0) && (nbtimebin != fTimeMax)) return 0;
      fTimeMax          = nbtimebin;

      //baseline          = rawStream->GetCommonAdditive();                // common additive baseline
     
      Int_t iTimeBin    = rawStream->GetTimeBin();                       //  current time bin
      Int_t *signal     = rawStream->GetSignals();                       //  current ADC signal
      Int_t col         = (rawStream->GetCol())%18;                      //  current COL MCM
           
      if((col < 0) || (col >= 21)) return 0;  
      if((imcm>=16) || (imcm < 0)) return 0;  
           
      Int_t fin     = TMath::Min(fTimeMax,(iTimeBin+3));
      Int_t n       = 0;
      for(Int_t itime = iTimeBin; itime < fin; itime++){
	if(signal[n]> (baseline+3)) phvalue[col][itime] = signal[n];
	n++;
      }
    }
    
    // fill the last one
    if(fDetectorPreviousTrack != -1){

      // take the mean values and check the first time bin
      for(Int_t j = 0; j < 21; j++){       
	if(TMath::RMS(fTimeMax,phvalue[j]) != 0.0) mean[j] = TMath::Mean(fTimeMax,phvalue[j]);
	else mean[j] = 0.0;
	  if(phvalue[j][0] > 200.0) first[j] = 1;
	  else first[j] = 0;
      }
      
      // select
      for(Int_t j = 1; j < 20; j++){
	if((first[j-1] == 0) && (first[j] ==0) && (first[j+1] == 0) && (mean[j-1] > (baseline+5.0)) && (mean[j] > (baseline+10.0)) && (mean[j+1] > (baseline+5.0)) && (mean[j] >= mean[j-1]) && (mean[j] >= mean[j+1])){
	  select[j] = 1;
	}
	else select[j] = 0;
      }
      
      // fill
      for(Int_t j = 1; j < 20; j++){
	if(select[j] == 1){
	  withInput = 2;
	  for(Int_t k = 0; k < fTimeMax; k++){
	    if((phvalue[j][k] >= phvalue[j-1][k]) && (phvalue[j][k] >= phvalue[j+1][k])){
	      UpdateDAQ(fDetectorPreviousTrack,0,0,k,(phvalue[j-1][k]+phvalue[j][k]+phvalue[j+1][k]),fTimeMax);
	    }
	    else{
	      if((j < 19) && (phvalue[j+1][k] >= phvalue[j][k]) && (phvalue[j+1][k] >= phvalue[j+2][k])){
		UpdateDAQ(fDetectorPreviousTrack,0,0,k,(phvalue[j][k]+phvalue[j+1][k]+phvalue[j+2][k]),fTimeMax);
	      }
	      else UpdateDAQ(fDetectorPreviousTrack,0,0,k,(3*baseline),fTimeMax);
	    }	   
	  }
	}
      }
      
      // reset
      for(Int_t k = 0; k < 36; k++){
	for(Int_t j = 0; j < 21; j++){
	  phvalue[j][k] = baseline;
	}
      }
    }
    
  }
  else{

    while (rawStream->Next()) {

      Int_t idetector = rawStream->GetDet();                            //  current detector
      Int_t imcm      = rawStream->GetMCM();                            //  current MCM
      Int_t irob      = rawStream->GetROB();                            //  current ROB

      if(((fMCMPrevious != imcm) || (fDetectorPreviousTrack != idetector) || (fROBPrevious != irob)) && (fDetectorPreviousTrack != -1)){

	// take the mean values and check the first time bin
	for(Int_t j = 0; j < 21; j++){       
	  if(TMath::RMS(fTimeMax,phvalue[j]) != 0.0) mean[j] = TMath::Mean(fTimeMax,phvalue[j]);
	  else mean[j] = 0.0;
	  if(phvalue[j][0] > 200.0) first[j] = 1;
	  else first[j] = 0;
	}
	
	// select
	for(Int_t j = 1; j < 20; j++){
	  if((first[j-1] == 0) && (first[j] ==0) && (first[j+1] == 0) && (mean[j-1] > (baseline+5.0)) && (mean[j] > (baseline+10.0)) && (mean[j+1] > (baseline+5.0)) && (mean[j] >= mean[j-1]) && (mean[j] >= mean[j+1])){
	    select[j] = 1;
	  }
	  else select[j] = 0;
	}
	
      // fill
	for(Int_t j = 1; j < 20; j++){
	  if(select[j] == 1){
	    withInput = 2;
	    for(Int_t k = 0; k < fTimeMax; k++){
	      if((phvalue[j][k] >= phvalue[j-1][k]) && (phvalue[j][k] >= phvalue[j+1][k])){
		UpdateDAQ(fDetectorPreviousTrack,0,0,k,(phvalue[j-1][k]+phvalue[j][k]+phvalue[j+1][k]),fTimeMax);
	      }
	      else{
		if((j < 19) && (phvalue[j+1][k] >= phvalue[j][k]) && (phvalue[j+1][k] >= phvalue[j+2][k])){
		  UpdateDAQ(fDetectorPreviousTrack,0,0,k,(phvalue[j][k]+phvalue[j+1][k]+phvalue[j+2][k]),fTimeMax);
		}
		else UpdateDAQ(fDetectorPreviousTrack,0,0,k,3*baseline,fTimeMax);
	      }	   
	    }
	  }
	}
	
	// reset
	for(Int_t k = 0; k < 36; k++){
	  for(Int_t j = 0; j < 21; j++){
	    phvalue[j][k] = baseline;
	  }
	}
      }
      
      fDetectorPreviousTrack = idetector;
      fMCMPrevious           = imcm;
      fROBPrevious           = irob;
      


      //baseline          = rawStream->GetCommonAdditive();                //  common baseline
      
      nbtimebin         = rawStream->GetNumberOfTimeBins();              //  number of time bins read from data
      Int_t iTimeBin    = rawStream->GetTimeBin();                       //  current time bin
      Int_t *signal     = rawStream->GetSignals();                       //  current ADC signal
      Int_t col         = (rawStream->GetCol())%18;                      //  current COL MCM

      Int_t fin     = TMath::Min(nbtimebin,(iTimeBin+3));
      Int_t n       = 0;
      
      if((col < 0) || (col >= 21)) return 0;  
      if((imcm>=16) || (imcm < 0)) return 0;  
      
      for(Int_t itime = iTimeBin; itime < fin; itime++){
	if(signal[n]>13) phvalue[col][itime] = signal[n];
	n++;
      }
    }
    
    // fill the last one
    if(fDetectorPreviousTrack != -1){
      
      // take the mean values and check the first time bin
      for(Int_t j = 0; j < 21; j++){       
	if(TMath::RMS(fTimeMax,phvalue[j]) != 0.0) mean[j] = TMath::Mean(fTimeMax,phvalue[j]);
	else mean[j] = 0.0;
	if(phvalue[j][0] > 200.0) first[j] = 1;
	else first[j] = 0;
      }
      
      // select
      for(Int_t j = 1; j < 20; j++){
	if((first[j-1] == 0) && (first[j] ==0) && (first[j+1] == 0) && (mean[j-1] > (baseline+5.0)) && (mean[j] > (baseline+10.0)) && (mean[j+1] > (baseline+5.0)) && (mean[j] >= mean[j-1]) && (mean[j] >= mean[j+1])){
	  select[j] = 1;
	}
	else select[j] = 0;
      }
      
      // fill
      for(Int_t j = 1; j < 20; j++){
	if(select[j] == 1){
	  withInput = 2;
	  for(Int_t k = 0; k < fTimeMax; k++){
	    if((phvalue[j][k] >= phvalue[j-1][k]) && (phvalue[j][k] >= phvalue[j+1][k])){
	      UpdateDAQ(fDetectorPreviousTrack,0,0,k,(phvalue[j-1][k]+phvalue[j][k]+phvalue[j+1][k]),fTimeMax);
	    }
	    else{
	      if((j < 19) && (phvalue[j+1][k] >= phvalue[j][k]) && (phvalue[j+1][k] >= phvalue[j+2][k])){
		UpdateDAQ(fDetectorPreviousTrack,0,0,k,(phvalue[j][k]+phvalue[j+1][k]+phvalue[j+2][k]),fTimeMax);
	      }
	      else UpdateDAQ(fDetectorPreviousTrack,0,0,k,3*baseline,fTimeMax);
	    }	   
	  }
	}
      }
      
      // reset
      for(Int_t k = 0; k < 36; k++){
	for(Int_t j = 0; j < 21; j++){
	  phvalue[j][k] = baseline;
	}
      }
    }
  }
  
  return withInput;
  
}
//_____________________________________________________________________
Int_t AliTRDCalibraFillHisto::ProcessEventDAQ(AliRawReader *rawReader, Bool_t nocheck)
{
  //
  //  Event processing loop - AliRawReader
  //


  AliTRDRawStreamV2 rawStream(rawReader);

  rawReader->Select("TRD");

  return ProcessEventDAQ(&rawStream, nocheck);
}
//_________________________________________________________________________
Int_t AliTRDCalibraFillHisto::ProcessEventDAQ(
#ifdef ALI_DATE
					       eventHeaderStruct *event,
					       Bool_t nocheck
#else
					       eventHeaderStruct* /*event*/,
					       Bool_t /*nocheck*/
	    
#endif 
				   )
{
  //
  //  process date event
  //
#ifdef ALI_DATE
    AliRawReader *rawReader = new AliRawReaderDate((void*)event);
    Int_t result=ProcessEventDAQ(rawReader, nocheck);
    delete rawReader;
    return result;
#else
    Fatal("AliTRDCalibraFillHisto", "this class was compiled without DATE");
    return 0;
#endif

}
//____________Online trackling in AliTRDtrigger________________________________
Bool_t AliTRDCalibraFillHisto::UpdateDAQ(Int_t det, Int_t /*row*/, Int_t /*col*/, Int_t timebin, Int_t signal, Int_t nbtimebins)
{
  //
  // For the DAQ
  // Fill a simple average pulse height
  //
  
  // Localisation of the Xbins involved
  //LocalisationDetectorXbins(det);

  // Row  and position in the pad groups
  //Int_t posr = 0;
  //if (fCalibraMode->GetNnZ(1) != 0) {
  //  posr = (Int_t) row / fCalibraMode->GetNnZ(1);
  //}
 
  // Col of cluster and position in the pad groups
  //Int_t posc = 0;
  //if (fCalibraMode->GetNnRphi(1) != 0) {
  //  posc = (Int_t) col / fCalibraMode->GetNnRphi(1);
  //}
  
  //fPH2d->Fill((Float_t) timebin/fSf,(fCalibraMode->GetXbins(1)+posc*fCalibraMode->GetNfragZ(1)+posr)+0.5,(Float_t) signal);   
  
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
    AnalyseLinearFitter();
    f.WriteTObject(fLinearVdriftFit);
  }
   
  f.Close();
  
  if ( backup ) backup->cd();
  
  AliInfo(Form("Execution time Write2d: R:%.2fs C:%.2fs"
	       ,stopwatch.RealTime(),stopwatch.CpuTime()));
}
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
      if((!((Bool_t)nbEntries)) && (nentries > 0)){
	nbEntries = new TH1F("Number of entries","Number of entries"
                               ,100,(Int_t)nentries/2,nentries*2);
	nbEntries->SetDirectory(0);
	nbEntriesPerGroup = new TH1F("Number of entries per group","Number of entries per group"
                               ,nbins,0,nbins);
	nbEntriesPerGroup->SetDirectory(0);
	nbEntriesPerSp = new TProfile("Number of entries per supermodule","Number of entries per supermodule"
                               ,(Int_t)(nbins/18),0,(Int_t)(nbins/18));
	nbEntriesPerSp->SetDirectory(0);
      }
      if(nbEntries){
	if(nentries > 0) nbEntries->Fill(nentries);
	nbEntriesPerGroup->Fill(idect+0.5,nentries);
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

  if(i > 1){
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

  Double_t *stat      = new Double_t[4]; 
  stat[3]             = 0.0;

  Int_t    nbofgroups = CalculateTotalNumberOfBins(0);
  Double_t *weight    = new Double_t[nbofgroups];
  Int_t    *nonul     = new Int_t[nbofgroups];
 
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

  return stat;

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
//____________Protected Functions______________________________________________
//____________Create the 2D histo to be filled online__________________________
//
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

  TString name("Nz");
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

  TString name("Nz");
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

//____________Offine tracking in the AliTRDtracker_____________________________
void AliTRDCalibraFillHisto::FillTheInfoOfTheTrackCH()
{
  //
  // For the offline tracking or mcm tracklets
  // This function will be called in the functions UpdateHistogram... 
  // to fill the info of a track for the relativ gain calibration
  //
	
  Int_t nb            =  0;   // Nombre de zones traversees
  Int_t fd            = -1;   // Premiere zone non nulle
  Float_t totalcharge = 0.0;  // Total charge for the supermodule histo

 
  
  
  // See if the track goes through different zones
  for (Int_t k = 0; k < fCalibraMode->GetNfragZ(0)*fCalibraMode->GetNfragRphi(0); k++) {
    if (fAmpTotal[k] > 0.0) {
      totalcharge += fAmpTotal[k];
      nb++;
      if (nb == 1) {
        fd = k;
      }
    }
  }

  
  switch (nb)
    { 
    case 1:
      fNumberUsedCh[0]++;
      fEntriesCH[fCalibraMode->GetXbins(0)+fd]++;
      if (fHisto2d) {
	FillCH2d(fCalibraMode->GetXbins(0)+fd,fAmpTotal[fd]/fRelativeScale);
 	//fCH2d->Fill(fAmpTotal[fd]/fRelativeScale,fCalibraMode->GetXbins(0)+fd+0.5);
      }
      if (fVector2d) {
	fCalibraVector->UpdateVectorCH(fDetectorPreviousTrack,fd,fAmpTotal[fd]/fRelativeScale);
      }
      break;
    case 2:
      if ((fAmpTotal[fd]   > 0.0) && 
	  (fAmpTotal[fd+1] > 0.0)) {
	// One of the two very big
	if (fAmpTotal[fd] > fProcent*fAmpTotal[fd+1]) {
	  if (fHisto2d) {
	    FillCH2d(fCalibraMode->GetXbins(0)+fd,fAmpTotal[fd]/fRelativeScale);
	    //fCH2d->Fill(fAmpTotal[fd]/fRelativeScale,fCalibraMode->GetXbins(0)+fd+0.5);
	  }
	  if (fVector2d) {
	    fCalibraVector->UpdateVectorCH(fDetectorPreviousTrack,fd,fAmpTotal[fd]/fRelativeScale);
	  }
	  fNumberUsedCh[1]++;
	  fEntriesCH[fCalibraMode->GetXbins(0)+fd]++;
	}
	if (fAmpTotal[fd+1] > fProcent*fAmpTotal[fd]) {
	  if (fHisto2d) {
	    FillCH2d(fCalibraMode->GetXbins(0)+fd+1,fAmpTotal[fd+1]/fRelativeScale);
	    //fCH2d->Fill(fAmpTotal[fd+1]/fRelativeScale,fCalibraMode->GetXbins(0)+fd+1.5);
	  }
	  if (fVector2d) {
	    fCalibraVector->UpdateVectorCH(fDetectorPreviousTrack,fd+1,fAmpTotal[fd+1]/fRelativeScale);
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
		  FillCH2d(fCalibraMode->GetXbins(0)+fd,fAmpTotal[fd]/fRelativeScale);
		  //fCH2d->Fill(fAmpTotal[fd]/fRelativeScale,fCalibraMode->GetXbins(0)+fd+0.5);
		}
		if (fVector2d) {
		  fCalibraVector->UpdateVectorCH(fDetectorPreviousTrack,fd,fAmpTotal[fd]/fRelativeScale);
		}
		fNumberUsedCh[1]++;
		fEntriesCH[fCalibraMode->GetXbins(0)+fd]++;
	      }
	      if (fAmpTotal[fd+fCalibraMode->GetNfragZ(0)] > fProcent*fAmpTotal[fd]) {
		if (fHisto2d) {
		  FillCH2d(fCalibraMode->GetXbins(0)+fd+fCalibraMode->GetNfragZ(0),fAmpTotal[fd+fCalibraMode->GetNfragZ(0)]/fRelativeScale);
		  //fCH2d->Fill(fAmpTotal[fd+fCalibraMode->GetNfragZ(0)]/fRelativeScale,fCalibraMode->GetXbins(0)+fd+fCalibraMode->GetNfragZ(0)+0.5);
		}
		fNumberUsedCh[1]++;
		fEntriesCH[fCalibraMode->GetXbins(0)+fd+fCalibraMode->GetNfragZ(0)]++;
		if (fVector2d) {
		  fCalibraVector->UpdateVectorCH(fDetectorPreviousTrack,fd+fCalibraMode->GetNfragZ(0),fAmpTotal[fd+fCalibraMode->GetNfragZ(0)]/fRelativeScale);
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

  // See if the track goes through different zones
  for (Int_t k = 0; k < fTimeMax; k++) {
    if (fPHValue[k] > 0.0) {
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

  
  switch(nb)
    {
    case 1:
      fNumberUsedPh[0]++;
      for (Int_t i = 0; i < fTimeMax; i++) {
	if (fHisto2d) {
	  fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) fPHValue[i]);
	}
	if (fVector2d) {
	  fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd1,i,fPHValue[i]);
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
	      fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) fPHValue[i]);
	    }
	    if (fVector2d) {
	      fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd1,i,fPHValue[i]);
	    }
	  }
	}
	if ((k2+fDifference) < fTimeMax) {
	  //we choose to fill the fd2 with all the values
	  fNumberUsedPh[1]++;
	  for (Int_t i = 0; i < fTimeMax; i++) {
	    if (fHisto2d) {
	      fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd2)+0.5,(Float_t) fPHValue[i]);
	    }
	  if (fVector2d) {
	    fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd2,i,fPHValue[i]);
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
		  fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) fPHValue[i]);
		}
		if (fVector2d) {
		  fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd1,i,fPHValue[i]);
		}
	      }
	    }
	    if ((k2+fDifference) < fTimeMax) {
	      //we choose to fill the fd2 with all the values
	      fNumberUsedPh[1]++;
	      for (Int_t i = 0; i < fTimeMax; i++) {
		if (fHisto2d) {
		  fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd2)+0.5,(Float_t) fPHValue[i]);
		}
		if (fVector2d) {
		  fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd2,i,fPHValue[i]);
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
		  fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) fPHValue[i]);
		}
		if (fVector2d) {
		  fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd1,i,fPHValue[i]);
		}
	      }
	    }
	    if ((k2+fDifference) < fTimeMax) {
	      //we choose to fill the fd2 with all the values
	      fNumberUsedPh[1]++;
	      for (Int_t i = 0; i < fTimeMax; i++) {
		if (fHisto2d) {
		  fPH2d->Fill((Float_t) i/fSf,(fCalibraMode->GetXbins(1)+fd2)+0.5,(Float_t) fPHValue[i]);
		}
		if (fVector2d) {
		  fCalibraVector->UpdateVectorPH(fDetectorPreviousTrack,fd2,i,fPHValue[i]);
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
Bool_t AliTRDCalibraFillHisto::FindP1TrackPH()
{
  //
  // For the offline tracking
  // This function will be called in the functions UpdateHistogram... 
  // to fill the find the parameter P1 of a track for the drift velocity  calibration
  //

   
  //Number of points: if less than 3 return kFALSE
  Int_t npoints = fListClusters->GetEntriesFast();
  if(npoints <= 2) return kFALSE;

  //Variables
  TLinearFitter linearFitterTracklet      = TLinearFitter(2,"pol1");        // TLinearFitter per tracklet
  Double_t snp                        = 0.0;                                // sin angle in the plan yx track
  Double_t y                          = 0.0;                                // y clusters in the middle of the chamber
  Double_t z                          = 0.0;                                // z cluster  in the middle of the chamber
  Double_t dydt                       = 0.0;                                // dydt tracklet after straight line fit
  Double_t tnp                        = 0.0;                                // tan angle in the plan xy track
  Double_t tgl                        = 0.0;                                // dz/dl and not dz/dx!  
  Double_t errorpar                   = 0.0;                                // error after straight line fit on dy/dt
  Double_t pointError                 = 0.0;                                // error after straight line fit 
  Int_t    detector                   = ((AliTRDcluster *) fListClusters->At(0))->GetDetector(); //detector
  Int_t    snpright                   = 1;                                  // if we took in the middle snp
  Int_t    crossrow                   = 0;                                  // if it crosses a pad row
  Double_t  tiltingangle              = 0;                                  // tiltingangle of the pad
  Float_t   dzdx                      = 0;                                  // dz/dx now from dz/dl
  Int_t     nbli                      = 0;                                  // number linear fitter points
  AliTRDpadPlane *padplane            = fGeo->GetPadPlane(GetPlane(detector),GetChamber(detector));

  linearFitterTracklet.StoreData(kFALSE);
  linearFitterTracklet.ClearPoints();
  
  //if more than one row
  Int_t    rowp                       = -1;                              // if it crosses a pad row

  //tiltingangle
  tiltingangle                        = padplane->GetTiltingAngle();
  Float_t  tnt                        = TMath::Tan(tiltingangle/180.*TMath::Pi()); // tan tiltingangle

  //Fill with points
  for(Int_t k = 0; k < npoints; k++){
    
    AliTRDcluster *cl                 = (AliTRDcluster *) fListClusters->At(k);
    Double_t ycluster                 = cl->GetY();
    Int_t time                        = cl->GetLocalTimeBin();
    Double_t timeis                   = time/fSf;
    //See if cross two pad rows
    Int_t    row                      = padplane->GetPadRowNumber(cl->GetZ());
    if(k==0) rowp                     = row;
    if(row != rowp) crossrow          = 1;
    //Take in the middle of the chamber
    //FollowBack
    if(time > (Int_t) 10) {
    //Follow
    //if(time < (Int_t) 11) {
      z   = cl->GetZ();
      y   = cl->GetY();  
      snp = fPar2[time];
      tgl = fPar3[time];
    }
    linearFitterTracklet.AddPoint(&timeis,ycluster,1);
    nbli++;
  }
  //FollowBack
  if(((AliTRDcluster *) fListClusters->At(0))->GetLocalTimeBin() < 10) snpright = 0;
  //Follow
  //if(((AliTRDcluster *) fListClusters->At(0))->GetLocalTimeBin() >= 11) snpright = 0;
  if(nbli <= 2) return kFALSE; 
  
  // Do the straight line fit now
  TVectorD pars;
  linearFitterTracklet.Eval();
  linearFitterTracklet.GetParameters(pars);
  pointError  =  TMath::Sqrt(linearFitterTracklet.GetChisquare()/nbli);
  errorpar    =  linearFitterTracklet.GetParError(1)*pointError;
  dydt  = pars[1]; 
 
  if( TMath::Abs(snp) <  1.){
    tnp = snp / (TMath::Sqrt(1-(snp*snp)));
  } 
  dzdx = tgl*TMath::Sqrt(1+tnp*tnp);

  if(fDebugLevel > 0){
    if ( !fDebugStreamer ) {
      //debug stream
      TDirectory *backup = gDirectory;
      fDebugStreamer = new TTreeSRedirector("TRDdebugCalibraFill.root");
      if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
    } 
    
    (* fDebugStreamer) << "VDRIFT0"<<
      "npoints="<<npoints<<
      "\n"; 
  
    
    (* fDebugStreamer) << "VDRIFT"<<
      "snpright="<<snpright<<
      "npoints="<<npoints<<
      "nbli="<<nbli<<
      "detector="<<detector<<
      "snp="<<snp<<
      "tnp="<<tnp<<
      "tgl="<<tgl<<
      "tnt="<<tnt<<
      "y="<<y<<
      "z="<<z<<
      "dydt="<<dydt<<
      "dzdx="<<dzdx<<
      "crossrow="<<crossrow<<
      "errorpar="<<errorpar<<
      "pointError="<<pointError<<
      "\n";     

  }
  
  if(npoints < fNumberClusters) return kFALSE;
  if(snpright == 0) return kFALSE;
  if(pointError >= 0.1) return kFALSE;
  if(crossrow == 1) return kFALSE;
  
  if(fLinearFitterOn){
    //Add to the linear fitter of the detector
    if( TMath::Abs(snp) <  1.){
      Double_t x = tnp-dzdx*tnt; 
      (GetLinearFitter(detector,kTRUE))->AddPoint(&x,dydt);
      if(fLinearFitterDebugOn) {
	fLinearVdriftFit->Update(detector,x,pars[1]);
      }
      fEntriesLinearFitter[detector]++;
    }
  }
  //AliInfo("End of FindP1TrackPH with success!")
  return kTRUE;

}
//____________Offine tracking in the AliTRDtracker_____________________________
Bool_t AliTRDCalibraFillHisto::HandlePRF()
{
  //
  // For the offline tracking
  // Fit the tracklet with a line and take the position as reference for the PRF
  //

  //Number of points
  Int_t npoints  = fListClusters->GetEntriesFast();                         // number of total points
  Int_t nb3pc    = 0;                                                       // number of three pads clusters used for fit 
  Int_t detector = ((AliTRDcluster *) fListClusters->At(0))->GetDetector(); // detector
 

  // To see the difference due to the fit
  Double_t *padPositions;
  padPositions = new Double_t[npoints];
  for(Int_t k = 0; k < npoints; k++){
    padPositions[k] = 0.0;
  } 


  //Find the position by a fit
  TLinearFitter fitter(2,"pol1");
  fitter.StoreData(kFALSE);
  fitter.ClearPoints();
  for(Int_t k = 0;  k < npoints; k++){
    //Take the cluster
    AliTRDcluster *cl  = (AliTRDcluster *) fListClusters->At(k);
    Short_t  *signals  = cl->GetSignals();
    Double_t     time  = cl->GetLocalTimeBin();
    //Calculate x if possible 
    Float_t xcenter    = 0.0;    
    Bool_t  echec      = kTRUE;   
    if((time<=7) || (time>=21)) continue; 
    // Center 3 balanced: position with the center of the pad
    if ((((Float_t) signals[3]) > 0.0) && 
	(((Float_t) signals[2]) > 0.0) && 
	(((Float_t) signals[4]) > 0.0)) {
      echec = kFALSE;
      // Security if the denomiateur is 0 
      if ((((Float_t) (((Float_t) signals[3]) * ((Float_t) signals[3]))) / 
	   ((Float_t) (((Float_t) signals[2]) * ((Float_t) signals[4])))) != 1.0) {
	xcenter = 0.5 * (TMath::Log((Float_t) (((Float_t) signals[4]) / ((Float_t) signals[2]))))
	  / (TMath::Log(((Float_t) (((Float_t) signals[3]) * ((Float_t) signals[3]))) 
			/ ((Float_t) (((Float_t) signals[2]) * ((Float_t) signals[4])))));
      }
      else {
	echec = kTRUE;
      }
    }
    if(TMath::Abs(xcenter) > 0.5) echec = kTRUE;
    if(echec) continue;
    //if no echec: calculate with the position of the pad
    // Position of the cluster
    Double_t       padPosition = xcenter +  cl->GetPadCol();
    padPositions[k]            = padPosition;
    nb3pc++;
    fitter.AddPoint(&time, padPosition,1);
  }//clusters loop

  //printf("nb3pc %d, npoints %d\n",nb3pc,npoints);
  if(nb3pc < 3) return kFALSE;
  fitter.Eval();
  TVectorD line(2);
  fitter.GetParameters(line);
  Float_t  pointError  = -1.0;
  pointError  =  TMath::Sqrt(fitter.GetChisquare()/nb3pc);
  

  // Now fill the PRF  
  for(Int_t k = 0;  k < npoints; k++){
    //Take the cluster
    AliTRDcluster *cl      = (AliTRDcluster *) fListClusters->At(k);
    Short_t  *signals      = cl->GetSignals();              // signal
    Double_t     time      = cl->GetLocalTimeBin();         // time bin
    Float_t padPosTracklet = line[0]+line[1]*time;          // reconstruct position from fit
    Float_t padPos         = cl->GetPadCol();               // middle pad
    Double_t dpad          = padPosTracklet - padPos;       // reconstruct position relative to middle pad from fit 
    Float_t ycenter        = 0.0;                           // relative center charge
    Float_t ymin           = 0.0;                           // relative left charge
    Float_t ymax           = 0.0;                           // relative right charge
    Double_t tgl           = fPar3[(Int_t)time];            // dz/dl and not dz/dx
    Double_t pt            = fPar4[(Int_t)time];            // pt
    Float_t  dzdx          = 0.0;                           // dzdx


    //Requiere simply two pads clusters at least
    if(((((Float_t) signals[3]) > 0.0) && (((Float_t) signals[2]) > 0.0)) ||
       ((((Float_t) signals[3]) > 0.0) && (((Float_t) signals[4]) > 0.0))){
      Float_t sum     = ((Float_t) signals[2]) + ((Float_t) signals[3]) + ((Float_t) signals[4]);
      if(sum > 0.0) ycenter = ((Float_t) signals[3])/ sum;
      if(sum > 0.0) ymin    = ((Float_t) signals[2])/ sum;
      if(sum > 0.0) ymax    = ((Float_t) signals[4])/ sum; 
    }
    
    //calibration group
    Int_t    *rowcol       = CalculateRowCol(cl);                       // calcul col and row pad of the cluster
    Int_t     grouplocal   = CalculateCalibrationGroup(2,rowcol);       // calcul the corresponding group
    Int_t     caligroup    = fCalibraMode->GetXbins(2)+ grouplocal;     // calcul the corresponding group
    Double_t  snp          = fPar2[(Int_t)time];                        // sin angle in xy plan
    Float_t   xcl          = cl->GetY();                                // y cluster
    Float_t   qcl          = cl->GetQ();                                // charge cluster 
    Int_t     plane        = GetPlane(detector);                        // plane 
    Int_t     chamber      = GetChamber(detector);                      // chamber  
    Double_t  xdiff        = dpad;                                      // reconstructed position constant
    Double_t  x            = dpad;                                      // reconstructed position moved
    Float_t   ep           = pointError;                                // error of fit
    Float_t   signal1      = (Float_t)signals[1];                       // signal at the border
    Float_t   signal3      = (Float_t)signals[3];                       // signal
    Float_t   signal2      = (Float_t)signals[2];                       // signal
    Float_t   signal4      = (Float_t)signals[4];                       // signal
    Float_t   signal5      = (Float_t)signals[5];                       // signal at the border
    Float_t tnp            = 0.0;
    if(TMath::Abs(snp) < 1.0){
      tnp = snp / (TMath::Sqrt(1-snp*snp));
      dzdx = tgl*TMath::Sqrt(1+tnp*tnp);
    }
   

    if(fDebugLevel > 0){
      if ( !fDebugStreamer ) {
	//debug stream
	TDirectory *backup = gDirectory;
	fDebugStreamer = new TTreeSRedirector("TRDdebugCalibraFill.root");
	if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
      }     
      
      (* fDebugStreamer) << "PRF0"<<
	"caligroup="<<caligroup<<
	"detector="<<detector<<
	"plane="<<plane<<
	"chamber="<<chamber<<
	"npoints="<<npoints<<
	"Np="<<nb3pc<<
	"ep="<<ep<<
	"snp="<<snp<<
	"tnp="<<tnp<<    
	"tgl="<<tgl<<  
	"dzdx="<<dzdx<<  
	"pt="<<pt<<    
	"padPos="<<padPos<<
	"padPosition="<<padPositions[k]<<
	"padPosTracklet="<<padPosTracklet<<
	"x="<<x<<
	"ycenter="<<ycenter<<
	"ymin="<<ymin<<
	"ymax="<<ymax<<
	"xcl="<<xcl<<
	"qcl="<<qcl<< 
	"signal1="<<signal1<<    
	"signal2="<<signal2<<
	"signal3="<<signal3<<
	"signal4="<<signal4<<
	"signal5="<<signal5<<
	"time="<<time<<
	"\n";     
      x = xdiff;
      Int_t type=0;
      Float_t y = ycenter;
      (* fDebugStreamer) << "PRFALL"<<
	"caligroup="<<caligroup<<
	"detector="<<detector<<
	"plane="<<plane<<
	"chamber="<<chamber<<
	"npoints="<<npoints<<
	"Np="<<nb3pc<<
	"ep="<<ep<<
	"type="<<type<<
	"snp="<<snp<<
	"tnp="<<tnp<<
	"tgl="<<tgl<<  
	"dzdx="<<dzdx<< 
	"pt="<<pt<< 
	"padPos="<<padPos<<
	"padPosition="<<padPositions[k]<<
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
      (* fDebugStreamer) << "PRFALL"<<
	"caligroup="<<caligroup<<
	"detector="<<detector<<
	"plane="<<plane<<
	"chamber="<<chamber<<
	"npoints="<<npoints<<
	"Np="<<nb3pc<<
	"ep="<<ep<<
	"type="<<type<<
	"snp="<<snp<<
	"tnp="<<tnp<<
	"tgl="<<tgl<<  
	"dzdx="<<dzdx<< 
	"pt="<<pt<< 
	"padPos="<<padPos<<
	"padPosition="<<padPositions[k]<<
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
      
      (* fDebugStreamer) << "PRFALL"<<
	"caligroup="<<caligroup<<
	"detector="<<detector<<
	"plane="<<plane<<
	"chamber="<<chamber<<
	"npoints="<<npoints<<
	"Np="<<nb3pc<<
	"ep="<<ep<<
	"type="<<type<<
	"snp="<<snp<<
	"tnp="<<tnp<<	
	"tgl="<<tgl<<  
	"dzdx="<<dzdx<< 
	"pt="<<pt<< 
	"padPos="<<padPos<<
	"padPosition="<<padPositions[k]<<
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
    
    // some cuts
    if(npoints < fNumberClusters) continue;
    if(nb3pc <= 5) continue;
    if((time >= 21) || (time < 7)) continue;
    if(TMath::Abs(snp) >= 1.0) continue;
    if(qcl < 80) continue; 
    
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
    if (fHisto2d && !echec) {
      if(TMath::Abs(dpad) < 1.5) {
	fPRF2d->Fill(shift+dpad,(caligroup+0.5),ycenter);
	fPRF2d->Fill(shift-dpad,(caligroup+0.5),ycenter);
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
    //Not equivalent anymore here!
    if (fVector2d && !echec) {
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
  }
  return kTRUE;
  
}
//____________Offine tracking in the AliTRDtracker_____________________________
Bool_t AliTRDCalibraFillHisto::FindP1TrackPHtrack(AliTRDtrack *t, Int_t index0, Int_t index1)
{
  //
  // For the offline tracking
  // This function will be called in the functions UpdateHistogram... 
  // to fill the find the parameter P1 of a track for the drift velocity  calibration
  //

    
  //Number of points: if less than 3 return kFALSE
  Int_t npoints = index1-index0;
  if(npoints <= 2) return kFALSE;

  //Variables
  TLinearFitter linearFitterTracklet  = TLinearFitter(2,"pol1");            // TLinearFitter per tracklet
  Double_t snp                        = 0.0;                                // sin angle in the plan yx track
  Double_t y                          = 0.0;                                // y clusters in the middle of the chamber
  Double_t z                          = 0.0;                                // z cluster  in the middle of the chamber
  Double_t dydt                       = 0.0;                                // dydt tracklet after straight line fit
  Double_t tnp                        = 0.0;                                // tan angle in the plan xy track
  Double_t tgl                        = 0.0;                                // dz/dl and not dz/dx!  
  Double_t errorpar                   = 0.0;                                // error after straight line fit on dy/dt
  Double_t pointError                 = 0.0;                                // error after straight line fit 
  Int_t    detector                   = ((AliTRDcluster *) t->GetCluster(index0))->GetDetector(); //detector
  //Int_t    snpright                   = 1;                                  // if we took in the middle snp
  Int_t    crossrow                   = 0;                                  // if it crosses a pad row
  Double_t  tiltingangle              = 0;                                  // tiltingangle of the pad
  Float_t   dzdx                      = 0;                                  // dz/dx now from dz/dl
  Int_t     nbli                      = 0;                                  // number linear fitter points
  AliTRDpadPlane *padplane            = fGeo->GetPadPlane(GetPlane(detector),GetChamber(detector));

  linearFitterTracklet.StoreData(kFALSE);
  linearFitterTracklet.ClearPoints();
  
  //if more than one row
  Int_t    rowp                       = -1;                              // if it crosses a pad row

  //tiltingangle
  tiltingangle                        = padplane->GetTiltingAngle();
  Float_t  tnt                        = TMath::Tan(tiltingangle/180.*TMath::Pi()); // tan tiltingangle

  //Fill with points
  for(Int_t k = 0; k < npoints; k++){
    
    AliTRDcluster *cl                 = (AliTRDcluster *) t->GetCluster(k+index0);
    Double_t ycluster                 = cl->GetY();
    Int_t time                        = cl->GetLocalTimeBin();
    Double_t timeis                   = time/fSf;
    //See if cross two pad rows
    Int_t    row                      = padplane->GetPadRowNumber(cl->GetZ());
    if(k==0) rowp                     = row;
    if(row != rowp) crossrow          = 1;
    //Take in the middle of the chamber
    //FollowBack
    //if(time > (Int_t) 10) {
    //Follow
    if(time < (Int_t) 11) {
      z   = cl->GetZ();
      y   = cl->GetY();  
    }
    linearFitterTracklet.AddPoint(&timeis,ycluster,1);
    nbli++;
  }

  // take now the snp, tnp and tgl from the track
  snp = t->GetSnpPlane(GetPlane(detector));
  tgl = t->GetTglPlane(GetPlane(detector));

  //FollowBack
  //if(((AliTRDcluster *) t->GetCluster(index0))->GetLocalTimeBin() < 10) snpright = 0;
  //Follow
  //if(((AliTRDcluster *) t->GetCluster(index0))->GetLocalTimeBin() >= 11) snpright = 0;
  if(nbli <= 2) return kFALSE; 
  
  // Do the straight line fit now
  TVectorD pars;
  linearFitterTracklet.Eval();
  linearFitterTracklet.GetParameters(pars);
  pointError  =  TMath::Sqrt(linearFitterTracklet.GetChisquare()/nbli);
  errorpar    =  linearFitterTracklet.GetParError(1)*pointError;
  dydt  = pars[1]; 
 
  if( TMath::Abs(snp) <  1.){
    tnp = snp / (TMath::Sqrt(1-(snp*snp)));
  } 
  dzdx = tgl*TMath::Sqrt(1+tnp*tnp);

  if(fDebugLevel > 0){
    if ( !fDebugStreamer ) {
      //debug stream
      TDirectory *backup = gDirectory;
      fDebugStreamer = new TTreeSRedirector("TRDdebugCalibraFill.root");
      if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
    } 
    
        
    (* fDebugStreamer) << "VDRIFT"<<
      //"snpright="<<snpright<<
      "npoints="<<npoints<<
      "nbli="<<nbli<<
      "detector="<<detector<<
      "snp="<<snp<<
      "tnp="<<tnp<<
      "tgl="<<tgl<<
      "tnt="<<tnt<<
      "y="<<y<<
      "z="<<z<<
      "dydt="<<dydt<<
      "dzdx="<<dzdx<<
      "crossrow="<<crossrow<<
      "errorpar="<<errorpar<<
      "pointError="<<pointError<<
      "\n";     

  }
  
  if(npoints < fNumberClusters) return kFALSE;
  //if(snpright == 0) return kFALSE;
  if(pointError >= 0.1) return kFALSE;
  if(crossrow == 1) return kFALSE;
  
  if(fLinearFitterOn){
    //Add to the linear fitter of the detector
    if( TMath::Abs(snp) <  1.){
      Double_t x = tnp-dzdx*tnt; 
      (GetLinearFitter(detector,kTRUE))->AddPoint(&x,dydt);
      if(fLinearFitterDebugOn) {
	fLinearVdriftFit->Update(detector,x,pars[1]);
      }
      fEntriesLinearFitter[detector]++;
    }
  }
  //AliInfo("End of FindP1TrackPH with success!")
  return kTRUE;

}
//____________Offine tracking in the AliTRDtracker_____________________________
Bool_t AliTRDCalibraFillHisto::HandlePRFtrack(AliTRDtrack *t, Int_t index0, Int_t index1)
{
  //
  // For the offline tracking
  // Fit the tracklet with a line and take the position as reference for the PRF
  //

  //Number of points
  Int_t npoints  = index1-index0;                                           // number of total points
  Int_t nb3pc    = 0;                                                       // number of three pads clusters used for fit 
  Int_t detector = ((AliTRDcluster *) t->GetCluster(index0))->GetDetector(); // detector
 

  // To see the difference due to the fit
  Double_t *padPositions;
  padPositions = new Double_t[npoints];
  for(Int_t k = 0; k < npoints; k++){
    padPositions[k] = 0.0;
  } 


  //Find the position by a fit
  TLinearFitter fitter(2,"pol1");
  fitter.StoreData(kFALSE);
  fitter.ClearPoints();
  for(Int_t k = 0;  k < npoints; k++){
    //Take the cluster
    AliTRDcluster *cl  = (AliTRDcluster *) t->GetCluster(k+index0);
    Short_t  *signals  = cl->GetSignals();
    Double_t     time  = cl->GetLocalTimeBin();
    //Calculate x if possible 
    Float_t xcenter    = 0.0;    
    Bool_t  echec      = kTRUE;   
    if((time<=7) || (time>=21)) continue; 
    // Center 3 balanced: position with the center of the pad
    if ((((Float_t) signals[3]) > 0.0) && 
	(((Float_t) signals[2]) > 0.0) && 
	(((Float_t) signals[4]) > 0.0)) {
      echec = kFALSE;
      // Security if the denomiateur is 0 
      if ((((Float_t) (((Float_t) signals[3]) * ((Float_t) signals[3]))) / 
	   ((Float_t) (((Float_t) signals[2]) * ((Float_t) signals[4])))) != 1.0) {
	xcenter = 0.5 * (TMath::Log((Float_t) (((Float_t) signals[4]) / ((Float_t) signals[2]))))
	  / (TMath::Log(((Float_t) (((Float_t) signals[3]) * ((Float_t) signals[3]))) 
			/ ((Float_t) (((Float_t) signals[2]) * ((Float_t) signals[4])))));
      }
      else {
	echec = kTRUE;
      }
    }
    if(TMath::Abs(xcenter) > 0.5) echec = kTRUE;
    if(echec) continue;
    //if no echec: calculate with the position of the pad
    // Position of the cluster
    Double_t       padPosition = xcenter +  cl->GetPadCol();
    padPositions[k]            = padPosition;
    nb3pc++;
    fitter.AddPoint(&time, padPosition,1);
  }//clusters loop

  //printf("nb3pc %d, npoints %d\n",nb3pc,npoints);
  if(nb3pc < 3) return kFALSE;
  fitter.Eval();
  TVectorD line(2);
  fitter.GetParameters(line);
  Float_t  pointError  = -1.0;
  pointError  =  TMath::Sqrt(fitter.GetChisquare()/nb3pc);
  
  // Take the tgl and snp with the track t now
  Double_t  tgl = t->GetTglPlane(GetPlane(detector)); //dz/dl and not dz/dx
  Double_t  snp = t->GetSnpPlane(GetPlane(detector)); // sin angle in xy plan
  Float_t  dzdx = 0.0;                                // dzdx
  Float_t  tnp  = 0.0;
  if(TMath::Abs(snp) < 1.0){
    tnp = snp / (TMath::Sqrt(1-snp*snp));
    dzdx = tgl*TMath::Sqrt(1+tnp*tnp);
  }


  // Now fill the PRF  
  for(Int_t k = 0;  k < npoints; k++){
    //Take the cluster
    AliTRDcluster *cl      = (AliTRDcluster *) t->GetCluster(k+index0);
    Short_t  *signals      = cl->GetSignals();              // signal
    Double_t     time      = cl->GetLocalTimeBin();         // time bin
    Float_t padPosTracklet = line[0]+line[1]*time;          // reconstruct position from fit
    Float_t padPos         = cl->GetPadCol();               // middle pad
    Double_t dpad          = padPosTracklet - padPos;       // reconstruct position relative to middle pad from fit 
    Float_t ycenter        = 0.0;                           // relative center charge
    Float_t ymin           = 0.0;                           // relative left charge
    Float_t ymax           = 0.0;                           // relative right charge
  


    //Requiere simply two pads clusters at least
    if(((((Float_t) signals[3]) > 0.0) && (((Float_t) signals[2]) > 0.0)) ||
       ((((Float_t) signals[3]) > 0.0) && (((Float_t) signals[4]) > 0.0))){
      Float_t sum     = ((Float_t) signals[2]) + ((Float_t) signals[3]) + ((Float_t) signals[4]);
      if(sum > 0.0) ycenter = ((Float_t) signals[3])/ sum;
      if(sum > 0.0) ymin    = ((Float_t) signals[2])/ sum;
      if(sum > 0.0) ymax    = ((Float_t) signals[4])/ sum; 
    }
    
    //calibration group
    Int_t    *rowcol       = CalculateRowCol(cl);                       // calcul col and row pad of the cluster
    Int_t     grouplocal   = CalculateCalibrationGroup(2,rowcol);       // calcul the corresponding group
    Int_t     caligroup    = fCalibraMode->GetXbins(2)+ grouplocal;     // calcul the corresponding group
    Float_t   xcl          = cl->GetY();                                // y cluster
    Float_t   qcl          = cl->GetQ();                                // charge cluster 
    Int_t     plane        = GetPlane(detector);                        // plane 
    Int_t     chamber      = GetChamber(detector);                      // chamber  
    Double_t  xdiff        = dpad;                                      // reconstructed position constant
    Double_t  x            = dpad;                                      // reconstructed position moved
    Float_t   ep           = pointError;                                // error of fit
    Float_t   signal1      = (Float_t)signals[1];                       // signal at the border
    Float_t   signal3      = (Float_t)signals[3];                       // signal
    Float_t   signal2      = (Float_t)signals[2];                       // signal
    Float_t   signal4      = (Float_t)signals[4];                       // signal
    Float_t   signal5      = (Float_t)signals[5];                       // signal at the border
   
   

    if(fDebugLevel > 0){
      if ( !fDebugStreamer ) {
	//debug stream
	TDirectory *backup = gDirectory;
	fDebugStreamer = new TTreeSRedirector("TRDdebugCalibraFill.root");
	if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
      }     
      
      (* fDebugStreamer) << "PRF0"<<
	"caligroup="<<caligroup<<
	"detector="<<detector<<
	"plane="<<plane<<
	"chamber="<<chamber<<
	"npoints="<<npoints<<
	"Np="<<nb3pc<<
	"ep="<<ep<<
	"snp="<<snp<<
	"tnp="<<tnp<<    
	"tgl="<<tgl<<  
	"dzdx="<<dzdx<<  
       	"padPos="<<padPos<<
	"padPosition="<<padPositions[k]<<
	"padPosTracklet="<<padPosTracklet<<
	"x="<<x<<
	"ycenter="<<ycenter<<
	"ymin="<<ymin<<
	"ymax="<<ymax<<
	"xcl="<<xcl<<
	"qcl="<<qcl<< 
	"signal1="<<signal1<<    
	"signal2="<<signal2<<
	"signal3="<<signal3<<
	"signal4="<<signal4<<
	"signal5="<<signal5<<
	"time="<<time<<
	"\n";     
      x = xdiff;
      Int_t type=0;
      Float_t y = ycenter;
      (* fDebugStreamer) << "PRFALL"<<
	"caligroup="<<caligroup<<
	"detector="<<detector<<
	"plane="<<plane<<
	"chamber="<<chamber<<
	"npoints="<<npoints<<
	"Np="<<nb3pc<<
	"ep="<<ep<<
	"type="<<type<<
	"snp="<<snp<<
	"tnp="<<tnp<<
	"tgl="<<tgl<<  
	"dzdx="<<dzdx<< 
	"padPos="<<padPos<<
	"padPosition="<<padPositions[k]<<
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
      (* fDebugStreamer) << "PRFALL"<<
	"caligroup="<<caligroup<<
	"detector="<<detector<<
	"plane="<<plane<<
	"chamber="<<chamber<<
	"npoints="<<npoints<<
	"Np="<<nb3pc<<
	"ep="<<ep<<
	"type="<<type<<
	"snp="<<snp<<
	"tnp="<<tnp<<
	"tgl="<<tgl<<  
	"dzdx="<<dzdx<< 
	"padPos="<<padPos<<
	"padPosition="<<padPositions[k]<<
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
      
      (* fDebugStreamer) << "PRFALL"<<
	"caligroup="<<caligroup<<
	"detector="<<detector<<
	"plane="<<plane<<
	"chamber="<<chamber<<
	"npoints="<<npoints<<
	"Np="<<nb3pc<<
	"ep="<<ep<<
	"type="<<type<<
	"snp="<<snp<<
	"tnp="<<tnp<<	
	"tgl="<<tgl<<  
	"dzdx="<<dzdx<< 
	"padPos="<<padPos<<
	"padPosition="<<padPositions[k]<<
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
    
    // some cuts
    if(npoints < fNumberClusters) continue;
    if(nb3pc <= 5) continue;
    if((time >= 21) || (time < 7)) continue;
    if(TMath::Abs(snp) >= 1.0) continue;
    if(qcl < 80) continue; 
    
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
    if (fHisto2d && !echec) {
      if(TMath::Abs(dpad) < 1.5) {
	fPRF2d->Fill(shift+dpad,(caligroup+0.5),ycenter);
	fPRF2d->Fill(shift-dpad,(caligroup+0.5),ycenter);
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
    //Not equivalent anymore here!
    if (fVector2d && !echec) {
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
  }
  return kTRUE;
  
}
//
//____________Some basic geometry function_____________________________________
//
//_____________________________________________________________________________
Int_t AliTRDCalibraFillHisto::GetPlane(Int_t d) const
{
  //
  // Reconstruct the plane number from the detector number
  //

  return ((Int_t) (d % 6));

}

//_____________________________________________________________________________
Int_t AliTRDCalibraFillHisto::GetChamber(Int_t d) const
{
  //
  // Reconstruct the chamber number from the detector number
  //
  Int_t fgkNplan = 6;

  return ((Int_t) (d % 30) / fgkNplan);

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
  
    /*
      AliInfo(Form("The pad calibration mode for the drift velocity calibration: Nz %d, and Nrphi %d"
      ,fCalibraMode->GetNz(1)
      ,fCalibraMode->GetNrphi(1)));
      
      // Calcul the number of Xbins
      Int_t ntotal1 = 0;
      fCalibraMode->ModePadCalibration(2,1);
      fCalibraMode->ModePadFragmentation(0,2,0,1);
      fCalibraMode->SetDetChamb2(1);
      ntotal1 += 6 * 18 * fCalibraMode->GetDetChamb2(1);
      fCalibraMode->ModePadCalibration(0,1);
      fCalibraMode->ModePadFragmentation(0,0,0,1);
      fCalibraMode->SetDetChamb0(1);
      ntotal1 += 6 * 4 * 18 * fCalibraMode->GetDetChamb0(1);
      AliInfo(Form("Total number of Xbins: %d",ntotal1));
      
      CreatePH2d(ntotal1);
    */  

    CreatePH2d(540);

    return fPH2d;
}
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
      linearfitter->Eval();
      linearfitter->GetParameters(*par);
      linearfitter->GetErrors(pare);
      Float_t  ppointError =  TMath::Sqrt(TMath::Abs(linearfitter->GetChisquare())/fEntriesLinearFitter[k]);
      (*parE)[0] = pare[0]*ppointError;
      (*parE)[1] = pare[1]*ppointError;
      (*parE)[2] = (Double_t) fEntriesLinearFitter[k];
      ((TObjArray *)fLinearVdriftFit->GetPArray())->AddAt(par,k);
      ((TObjArray *)fLinearVdriftFit->GetEArray())->AddAt(parE,k);
    }
  }
}
//________________________________________________________________________________
Int_t AliTRDCalibraFillHisto::Arrondi(Double_t x) const
{
   // Partie entiere of the (x+0.5)
   
  int i;
  if (x >= (-0.5)) {
    i = int(x + 0.5);
    //if (x + 0.5 == Float_t(i) && i & 1) i--;
  } else {
    i = int(x - 0.5);
    //if (x - 0.5 == Float_t(i) && i & 1) i++;
    if((x-0.5)==i) i++;

  }
  return i;
}
