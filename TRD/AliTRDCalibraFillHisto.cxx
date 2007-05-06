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
// The user has to choose with the functions SetNz and SetNrphi the precision of the
// calibration (see AliTRDCalibraMode). 
// 2D Histograms (Histo2d) or vectors (Vector2d), then converted in Trees, will be filled
// from RAW DATA in a run or from reconstructed TRD tracks during the offline tracking 
// in the function "FollowBackProlongation" (AliTRDtracker)
// Per default the functions to fill are off.                                   
//                        
// Author:
//   R. Bailhache (R.Bailhache@gsi.de)
//                            
//////////////////////////////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TChain.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TH2.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TROOT.h>

#include "AliLog.h"
#include "AliCDBManager.h"

#include "AliTRDCalibraFillHisto.h"
#include "AliTRDCalibraMode.h"
#include "AliTRDCalibraVector.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "AliTRDmcmTracklet.h"
#include "AliTRDpadPlane.h"
#include "AliTRDcluster.h"
#include "AliTRDtrack.h"


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
  ,fMITracking(kFALSE)
  ,fMcmTracking(kFALSE)
  ,fMcmCorrectAngle(kFALSE)
  ,fCH2dOn(kFALSE)
  ,fPH2dOn(kFALSE)
  ,fPRF2dOn(kFALSE)
  ,fHisto2d(kFALSE)
  ,fVector2d(kFALSE)
  ,fRelativeScale(0)
  ,fCountRelativeScale(0)
  ,fRelativeScaleAuto(kFALSE)
  ,fThresholdClusterPRF1(0.0)
  ,fThresholdClusterPRF2(0.0)
  ,fCenterOfflineCluster(kFALSE)
  ,fWriteName(0)
  ,fCalibraMode(0)
  ,fDetectorAliTRDtrack(kFALSE)
  ,fChamberAliTRDtrack(-1)
  ,fDetectorPreviousTrack(-1)
  ,fGoodTrack(kTRUE)
  ,fAmpTotal(0x0)
  ,fPHPlace(0x0)
  ,fPHValue(0x0)
  ,fNumberClusters(0)
  ,fProcent(0.0)
  ,fDifference(0)
  ,fNumberTrack(0)
  ,fTimeMax(0)
  ,fSf(0.0)
  ,fNumberBinCharge(0)
  ,fNumberBinPRF(0)
  ,fCalibraVector(0)
  ,fPH2d(0x0)
  ,fPRF2d(0x0)
  ,fCH2d(0x0)
{
  //
  // Default constructor
  //

  fCalibraMode = new AliTRDCalibraMode();

  // Write
  for (Int_t i = 0; i < 3; i++) {
    fWrite[i]     = kFALSE;
  }

  // Init
  Init();
  
}

//______________________________________________________________________________________
AliTRDCalibraFillHisto::AliTRDCalibraFillHisto(const AliTRDCalibraFillHisto &c)
  :TObject(c)
  ,fMITracking(kFALSE)
  ,fMcmTracking(kFALSE)
  ,fMcmCorrectAngle(kFALSE)
  ,fCH2dOn(kFALSE)
  ,fPH2dOn(kFALSE)
  ,fPRF2dOn(kFALSE)
  ,fHisto2d(kFALSE)
  ,fVector2d(kFALSE)
  ,fRelativeScale(0)
  ,fCountRelativeScale(0)
  ,fRelativeScaleAuto(kFALSE)
  ,fThresholdClusterPRF1(0.0)
  ,fThresholdClusterPRF2(0.0)
  ,fCenterOfflineCluster(kFALSE)
  ,fWriteName(0)
  ,fCalibraMode(0)
  ,fDetectorAliTRDtrack(kFALSE)
  ,fChamberAliTRDtrack(-1)
  ,fDetectorPreviousTrack(-1)
  ,fGoodTrack(kTRUE)
  ,fAmpTotal(0x0)
  ,fPHPlace(0x0)
  ,fPHValue(0x0)
  ,fNumberClusters(0)
  ,fProcent(0.0)
  ,fDifference(0)
  ,fNumberTrack(0)
  ,fTimeMax(0)
  ,fSf(0.0)
  ,fNumberBinCharge(0)
  ,fNumberBinPRF(0)
  ,fCalibraVector(0)
  ,fPH2d(0x0)
  ,fPRF2d(0x0)
  ,fCH2d(0x0) 
{
  //
  // Copy constructor
  //

}

//____________________________________________________________________________________
AliTRDCalibraFillHisto::~AliTRDCalibraFillHisto()
{
  //
  // AliTRDCalibraFillHisto destructor
  //

  ClearHistos();
  
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

//_____________________________________________________________________________
void AliTRDCalibraFillHisto::Init() 
{
  //
  // Init some default values
  //

  // How to fill the 2D
  fThresholdClusterPRF1 = 2.0;
  fThresholdClusterPRF2 = 15.0;
  
  // Store the Info
  fNumberBinCharge      = 100;
  fNumberBinPRF         = 40;
  
  // Write
  fWriteName            = "TRD.calibration.root";
  
  // Internal variables
  
  // Fill the 2D histos in the offline tracking
  fDetectorPreviousTrack = -1;
  fChamberAliTRDtrack    = -1;
  fGoodTrack             = kTRUE;

  fProcent               = 6.0;
  fDifference            = 17;
  fNumberClusters        = 18;
  fNumberTrack           = 0;
  fNumberUsedCh[0]       = 0;
  fNumberUsedCh[1]       = 0;
  fNumberUsedPh[0]       = 0;
  fNumberUsedPh[1]       = 0;
 
}

//____________Functions for initialising the AliTRDCalibraFillHisto in the code_________
Bool_t AliTRDCalibraFillHisto::Init2Dhistos()
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
  if (fRelativeScaleAuto) {
    fRelativeScale = 0;
  }
  else {
    fRelativeScale = 20;
  }

  //If vector method On initialised all the stuff
  if(fVector2d){
    fCalibraVector = new AliTRDCalibraVector();
  }


  // Create the 2D histos corresponding to the pad groupCalibration mode
  if (fCH2dOn) {

    AliInfo(Form("The pad calibration mode for the relative gain calibration: Nz %d, and Nrphi %d"
                ,fCalibraMode->GetNz(0)
                ,fCalibraMode->GetNrphi(0)));
    
    // Calcul the number of Xbins
    Int_t Ntotal0 = 0;
    fCalibraMode->ModePadCalibration(2,0);
    fCalibraMode->ModePadFragmentation(0,2,0,0);
    fCalibraMode->SetDetChamb2(0);
    Ntotal0 += 6 * 18 * fCalibraMode->GetDetChamb2(0);
    fCalibraMode->ModePadCalibration(0,0);
    fCalibraMode->ModePadFragmentation(0,0,0,0);
    fCalibraMode->SetDetChamb0(0);
    Ntotal0 += 6 * 4 * 18 * fCalibraMode->GetDetChamb0(0);
    AliInfo(Form("Total number of Xbins: %d",Ntotal0));

    // Create the 2D histo
    if (fHisto2d) {
      CreateCH2d(Ntotal0);
    }
    if (fVector2d) {
      fCalibraVector->SetNumberBinCharge(fNumberBinCharge);
    }

    // Variable
    fAmpTotal = new Float_t[TMath::Max(fCalibraMode->GetDetChamb2(0),fCalibraMode->GetDetChamb0(0))];
    for (Int_t k = 0; k < TMath::Max(fCalibraMode->GetDetChamb2(0),fCalibraMode->GetDetChamb0(0)); k++) {
      fAmpTotal[k] = 0.0;
    } 

  }

  if (fPH2dOn) {

    AliInfo(Form("The pad calibration mode for the drift velocity calibration: Nz %d, and Nrphi %d"
                ,fCalibraMode->GetNz(1)
                ,fCalibraMode->GetNrphi(1)));
    
    // Calcul the number of Xbins
    Int_t Ntotal1 = 0;
    fCalibraMode->ModePadCalibration(2,1);
    fCalibraMode->ModePadFragmentation(0,2,0,1);
    fCalibraMode->SetDetChamb2(1);
    Ntotal1 += 6 * 18 * fCalibraMode->GetDetChamb2(1);
    fCalibraMode->ModePadCalibration(0,1);
    fCalibraMode->ModePadFragmentation(0,0,0,1);
    fCalibraMode->SetDetChamb0(1);
    Ntotal1 += 6 * 4 * 18 * fCalibraMode->GetDetChamb0(1);
    AliInfo(Form("Total number of Xbins: %d",Ntotal1));

    // Create the 2D histo
    if (fHisto2d) {
      CreatePH2d(Ntotal1);
    }
    if (fVector2d) {
      fCalibraVector->SetTimeMax(fTimeMax);
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

  if (fPRF2dOn) {

    AliInfo(Form("The pad calibration mode for the PRF calibration: Nz %d, and Nrphi %d"
                ,fCalibraMode->GetNz(2)
                ,fCalibraMode->GetNrphi(2)));
    
    // Calcul the number of Xbins
    Int_t Ntotal2 = 0;
    fCalibraMode->ModePadCalibration(2,2);
    fCalibraMode->ModePadFragmentation(0,2,0,2);
    fCalibraMode->SetDetChamb2(2);
    Ntotal2 += 6 * 18 * fCalibraMode->GetDetChamb2(2);
    fCalibraMode->ModePadCalibration(0,2);
    fCalibraMode->ModePadFragmentation(0,0,0,2);
    fCalibraMode->SetDetChamb0(2);
    Ntotal2 += 6 * 4 * 18 * fCalibraMode->GetDetChamb0(2);
    AliInfo(Form("Total number of Xbins: %d",Ntotal2));

    // Create the 2D histo
    if (fHisto2d) {
      CreatePRF2d(Ntotal2);
    }
    if (fVector2d) {
      fCalibraVector->SetNumberBinPRF(fNumberBinPRF);
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

  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam");
    return kFALSE;
  }

  // Get the parameter object
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }
 
  // Localisation of the detector
  Int_t detector = cl->GetDetector();
  Int_t chamber  = GetChamber(detector);
  Int_t plane    = GetPlane(detector);

  // Fill the infos for the previous clusters if not the same
  // detector anymore or if not the same track
  if (((detector != fDetectorPreviousTrack) || (!fDetectorAliTRDtrack)) && 
      (fDetectorPreviousTrack != -1)) {
    
    fNumberTrack++;   
    
    // If the same track, then look if the previous detector is in
    // the same plane, if yes: not a good track
    if (fDetectorAliTRDtrack && 
        (GetPlane(detector) <= GetPlane(fDetectorPreviousTrack))) {
      fGoodTrack = kFALSE;
    }

    // Fill only if the track doesn't touch a masked pad or doesn't
    // appear in the middle (fGoodTrack)
    if (fGoodTrack) {

      // Gain calibration
      if (fCH2dOn) {
	FillTheInfoOfTheTrackCH();
      }
      
      // PH calibration
      if (fPH2dOn) {
	FillTheInfoOfTheTrackPH();    
      }
      
    } // if a good track
    
    ResetfVariables();
    
  } // Fill at the end the charge
  
  // Calcul the position of the detector
  if (detector != fDetectorPreviousTrack) {
    LocalisationDetectorXbins(detector);
  }

  // Reset the good track for the PRF
  Bool_t good = kTRUE;
  
  // Localisation of the cluster
  Double_t pos[3] = { 0.0, 0.0, 0.0 };
  pos[0] = cl->GetX();
  pos[1] = cl->GetY();
  pos[2] = cl->GetZ();
  Int_t    time   = cl->GetLocalTimeBin();
  
  // Reset the detector
  fDetectorPreviousTrack = detector;
  fDetectorAliTRDtrack   = kTRUE;
  
  // Position of the cluster
  AliTRDpadPlane *padplane = parCom->GetPadPlane(plane,chamber);
  Int_t    row        = padplane->GetPadRowNumber(pos[2]);
  Double_t offsetz    = padplane->GetPadRowOffset(row,pos[2]);
  Double_t offsettilt = padplane->GetTiltOffset(offsetz);
  Int_t    col        = padplane->GetPadColNumber(pos[1]+offsettilt,offsetz);
  
  // See if we are not near a masked pad
  if (!IsPadOn(detector,col,row)) {
    good       = kFALSE;
    fGoodTrack = kFALSE;
  }

  if (col > 0) {
    if (!IsPadOn(detector,col-1,row)) {
      fGoodTrack = kFALSE;
      good       = kFALSE;
    }
  }

  if (col < 143) {
    if (!IsPadOn(detector,col+1,row)) {
      fGoodTrack = kFALSE;
      good       = kFALSE;
    }
  }

  // Row of the cluster and position in the pad groups
  Int_t posr[3] = { 0, 0, 0 };
  if ((fCH2dOn)  && (fCalibraMode->GetNnZ(0) != 0)) {
    posr[0] = (Int_t) row / fCalibraMode->GetNnZ(0);
  }
  if ((fPH2dOn)  && (fCalibraMode->GetNnZ(1) != 0)) {
    posr[1] = (Int_t) row / fCalibraMode->GetNnZ(1);
  }
  if ((fPRF2dOn) && (fCalibraMode->GetNnZ(2) != 0)) {
    posr[2] = (Int_t) row / fCalibraMode->GetNnZ(2);
  }  
      
  // Col of the cluster and position in the pad groups
  Int_t posc[3] = { 0, 0, 0 };
  if ((fCH2dOn)  && (fCalibraMode->GetNnRphi(0) != 0)) {
    posc[0] = (Int_t) col / fCalibraMode->GetNnRphi(0);
  }
  if ((fPH2dOn)  && (fCalibraMode->GetNnRphi(1) != 0)) {
    posc[1] = (Int_t) col / fCalibraMode->GetNnRphi(1);
  }
  if ((fPRF2dOn) && (fCalibraMode->GetNnRphi(2) != 0)) {
    posc[2] = (Int_t) col / fCalibraMode->GetNnRphi(2);
  }

  // Charge in the cluster
  // For the moment take the abs
  Float_t  q        = TMath::Abs(cl->GetQ());
  Short_t  *signals = cl->GetSignals();
 
  // Correction due to the track angle
  Float_t correction    = 1.0;
  Float_t normalisation = 6.67;
  if ((q >0) && (t->GetNdedx() > 0)) {
    correction = t->GetClusterdQdl((t->GetNdedx() - 1)) / (normalisation);
  }

  // Fill the fAmpTotal with the charge
  if (fCH2dOn) {
    fAmpTotal[(Int_t) (posc[0]*fCalibraMode->GetNfragZ(0)+posr[0])] += correction;
  }

  // Fill the fPHPlace and value
  if (fPH2dOn) {
    fPHPlace[time] = posc[1]*fCalibraMode->GetNfragZ(1)+posr[1];
    fPHValue[time] = correction;
  }

  // Fill direct the PRF
  if ((fPRF2dOn) && (good)) {

    Float_t yminus  = 0.0;
    Float_t xcenter = 0.0;
    Float_t ycenter = 0.0;
    Float_t ymax    = 0.0;
    Bool_t  echec   = kFALSE;
    
    if ((cl->From3pad()) && (!cl->IsUsed())) { 
         
      // Center 3 balanced and cut on the cluster shape
      if ((((Float_t) signals[3]) > fThresholdClusterPRF2) && 
          (((Float_t) signals[2]) > fThresholdClusterPRF2) && 
          (((Float_t) signals[4]) > fThresholdClusterPRF2) && 
          (((Float_t) signals[1]) < fThresholdClusterPRF1) && 
          (((Float_t) signals[5]) < fThresholdClusterPRF1) && 
          ((((Float_t) signals[2])*((Float_t) signals[4])/(((Float_t) signals[3])*((Float_t) signals[3]))) < 0.06)) {


	//First calculate the position of the cluster and the y
	//echec enables to repair cases where it fails
	//
	// Col correspond to signals[3]
	if (fCenterOfflineCluster) {
          xcenter = cl->GetCenter();
	}
	else {
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
	//after having calculating the position calculate the y
	if (TMath::Abs(xcenter) < 0.5) {
	  ycenter = (Float_t) (((Float_t) signals[3]) 
                            / (((Float_t) signals[2]) + ((Float_t) signals[3]) + (((Float_t) signals[4]))));
	  yminus  = (Float_t) (((Float_t) signals[2]) 
                            / (((Float_t) signals[2]) + ((Float_t) signals[3]) + (((Float_t) signals[4]))));
	  ymax    = (Float_t) (((Float_t) signals[4]) 
                            / (((Float_t) signals[2]) + ((Float_t) signals[3]) + (((Float_t) signals[4]))));
	  //If the charge of the cluster is too far away from the corrected one cut
	  if ((TMath::Abs(((Float_t) signals[2]) + ((Float_t) signals[3]) + (((Float_t) signals[4])) - q) > 10.0)) {
            echec = kTRUE;
	  }
	}
	else {
	  echec = kTRUE;
	}
      
	//Then Fill the histo if no echec
	//
	// Fill only if it is in the drift region!
	if ((((Float_t) (((Float_t) time) / fSf)) > 0.3) && (!echec)) {
	  if (fHisto2d) {
	    fPRF2d->Fill((fCalibraMode->GetXbins(2)+posc[2]*fCalibraMode->GetNfragZ(2)+posr[2]+0.5),xcenter,ycenter);
	    fPRF2d->Fill((fCalibraMode->GetXbins(2)+posc[2]*fCalibraMode->GetNfragZ(2)+posr[2]+0.5),-(xcenter+1.0),yminus);
	    fPRF2d->Fill((fCalibraMode->GetXbins(2)+posc[2]*fCalibraMode->GetNfragZ(2)+posr[2]+0.5),1.0-xcenter,ymax);
	  }
	  if (fVector2d) {
	    fCalibraVector->UpdateVectorPRF(fCalibraMode->GetXbins(2)+posc[2]*fCalibraMode->GetNfragZ(2)+posr[2],xcenter,ycenter);
	    fCalibraVector->UpdateVectorPRF(fCalibraMode->GetXbins(2)+posc[2]*fCalibraMode->GetNfragZ(2)+posr[2],-(xcenter+1.0),yminus);
	    fCalibraVector->UpdateVectorPRF(fCalibraMode->GetXbins(2)+posc[2]*fCalibraMode->GetNfragZ(2)+posr[2],1.0-xcenter,ymax);
	  }
	} // If in the drift region
      } // center 3 balanced and cut on the cluster shape
    } // Cluster isole
  } // PRF2dOn	
  
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
  LocalisationDetectorXbins(idect);

  // Get the parameter object
  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }
   
  // Reset
  ResetfVariables();
 
  // Row of the tracklet and position in the pad groups
  Int_t row     = trk->GetRow();
  Int_t posr[3] = { 0, 0, 0 };
  if ((fCH2dOn)  && (fCalibraMode->GetNnZ(0) != 0)) {
    posr[0] = (Int_t) row / fCalibraMode->GetNnZ(0);
  }
  if ((fPH2dOn)  && (fCalibraMode->GetNnZ(1) != 0)) {
    posr[1] = (Int_t) row / fCalibraMode->GetNnZ(1);
  }
  if ((fPRF2dOn) && (fCalibraMode->GetNnZ(2) != 0)) {
    posr[2] = (Int_t) row / fCalibraMode->GetNnZ(2);
  }
 
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
      Int_t   col    = trk->GetClusterCol(icl);
            
      amp[0] = trk->GetClusterADC(icl)[0] * correction;
      amp[1] = trk->GetClusterADC(icl)[1] * correction;
      amp[2] = trk->GetClusterADC(icl)[2] * correction;

      
      if ((amp[0] < 0.0) || 
          (amp[1] < 0.0) || 
          (amp[2] < 0.0)) {
        continue;
      }

      // Col of cluster and position in the pad groups
      Int_t posc[3] = { 0, 0, 0 };
      if ((fCH2dOn)  && (fCalibraMode->GetNnRphi(0) != 0)) {
        posc[0] = (Int_t) col / fCalibraMode->GetNnRphi(0);
      }
      if ((fPH2dOn)  && (fCalibraMode->GetNnRphi(1) != 0)) {
        posc[1] = (Int_t) col / fCalibraMode->GetNnRphi(1);
      }
      if ((fPRF2dOn) && (fCalibraMode->GetNnRphi(2) != 0)) {
        posc[2] = (Int_t) col / fCalibraMode->GetNnRphi(2);
      }

      // See if we are not near a masked pad
      Bool_t good = kTRUE;
      if (!IsPadOn(idect,col,row)) {
	fGoodTrack = kFALSE;
	good       = kFALSE;
      }

      if (col >   0) {
	if (!IsPadOn(idect,col-1,row)) {
	  fGoodTrack = kFALSE;
	  good       = kFALSE;
	}
      }
      
      if (col < 143) {
	if (!IsPadOn(idect,col+1,row)) {
	  fGoodTrack = kFALSE;
	  good       = kFALSE;
	}
      }

      // Total spectrum
      if (fPH2dOn) {
        fPHPlace[time] = posc[1] * fCalibraMode->GetNfragZ(1) + posr[1];
      }

      if (fCH2dOn) {
	fAmpTotal[(Int_t) (posc[0]*fCalibraMode->GetNfragZ(0)+posr[0])] += (Float_t) (amp[0]+amp[1]+amp[2]);
      }
      if (fPH2dOn) {
	fPHValue[time] = (Float_t) (amp[0]+amp[1]+amp[2]);
      }
      
            
      // Fill PRF direct
      if (fPRF2dOn && good) {
	
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
		  fPRF2d->Fill((fCalibraMode->GetXbins(2)+posc[2]*fCalibraMode->GetNfragZ(2)+posr[2]+0.5),xcenter,ycenter);
		  fPRF2d->Fill((fCalibraMode->GetXbins(2)+posc[2]*fCalibraMode->GetNfragZ(2)+posr[2]+0.5),-(xcenter+1.0),yminus);
		  fPRF2d->Fill((fCalibraMode->GetXbins(2)+posc[2]*fCalibraMode->GetNfragZ(2)+posr[2]+0.5),(1.0-xcenter),ymax);
		}
		if (fVector2d) {
		  fCalibraVector->UpdateVectorPRF((fCalibraMode->GetXbins(2)+posc[2]*fCalibraMode->GetNfragZ(2)+posr[2]),xcenter,ycenter);
		  fCalibraVector->UpdateVectorPRF(fCalibraMode->GetXbins(2)+posc[2]*fCalibraMode->GetNfragZ(2)+posr[2],-(xcenter+1.0),yminus);
		  fCalibraVector->UpdateVectorPRF(fCalibraMode->GetXbins(2)+posc[2]*fCalibraMode->GetNfragZ(2)+posr[2],(1.0-xcenter),ymax);
		}
	      }//in the drift region 
	    }//in the middle
	  }//denominateur security
	}//cluster shape and thresholds
      }//good and PRF On
      
    } // Boucle clusters
    
    // Fill the charge
    if (fCH2dOn && fGoodTrack) {
      FillTheInfoOfTheTrackCH();
    }

    // PH calibration
    if (fPH2dOn && fGoodTrack) {
      FillTheInfoOfTheTrackPH();	
    }

    fNumberTrack++;
        
  } // Condition on number of clusters

  return kTRUE;
  
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

//____________Functions for plotting the 2D____________________________________

//_____________________________________________________________________________
void AliTRDCalibraFillHisto::Plot2d()
{
  //
  // Plot the 2D histos 
  //
 
  if (fPH2dOn) {
    TCanvas *cph2d = new TCanvas("cph2d","",50,50,600,800);
    cph2d->cd();
    fPH2d->Draw("LEGO");
  }
  if (fCH2dOn) {
    TCanvas *cch2d = new TCanvas("cch2d","",50,50,600,800);
    cch2d->cd();
    fCH2d->Draw("LEGO");
  }
  if (fPRF2dOn) {
    TCanvas *cPRF2d = new TCanvas("cPRF2d","",50,50,600,800);
    cPRF2d->cd();
    fPRF2d->Draw("LEGO");
  }

}

//____________Writing the 2D___________________________________________________

//_____________________________________________________________________________
Bool_t AliTRDCalibraFillHisto::Write2d()
{
  //
  // Write the 2D histograms or the vectors converted in trees in the file
  // "TRD.calibration.root" 
  //
  
  TFile *fout = TFile::Open(fWriteName,"RECREATE");
  // Check if the file could be opened
  if (!fout || !fout->IsOpen()) {
    AliInfo("No File found!");
    return kFALSE;
  }
  AliInfo(Form("Numbertrack: %d Numberusedch[0]: %d, Numberusedch[1]: %d Numberusedph[0]: %d, Numberusedph[1]: %d"
              ,fNumberTrack
              ,fNumberUsedCh[0]
              ,fNumberUsedCh[1]
              ,fNumberUsedPh[0]
              ,fNumberUsedPh[1]));
  
  TStopwatch stopwatch;
  stopwatch.Start();
  AliInfo("Write2d");

  if ((fCH2dOn ) && (fWrite[0])) {
    if (fHisto2d) {
      fout->WriteTObject(fCH2d);
    }
    if (fVector2d) {
      TString name("Nz");
      name += fCalibraMode->GetNz(0);
      name += "Nrphi";
      name += fCalibraMode->GetNrphi(0);
      TTree *treeCH2d = fCalibraVector->ConvertVectorCTTreeHisto(fCalibraVector->GetVectorCH(),fCalibraVector->GetPlaCH(),"treeCH2d",(const char *) name);
      fout->WriteTObject(treeCH2d);
    }
  }
  if ((fPH2dOn ) && (fWrite[1])) {
    if (fHisto2d) {
      fout->WriteTObject(fPH2d);
    }
    if (fVector2d) {
      TString name("Nz");
      name += fCalibraMode->GetNz(1);
      name += "Nrphi";
      name += fCalibraMode->GetNrphi(1);
      TTree *treePH2d = fCalibraVector->ConvertVectorPTreeHisto(fCalibraVector->GetVectorPH(),fCalibraVector->GetPlaPH(),"treePH2d",(const char *) name);
      fout->WriteTObject(treePH2d);
    }
  }
  if ((fPRF2dOn ) && (fWrite[2])) {
    if (fHisto2d) {
      fout->WriteTObject(fPRF2d);
    }
    if (fVector2d) {
      TString name("Nz");
      name += fCalibraMode->GetNz(2);
      name += "Nrphi";
      name += fCalibraMode->GetNrphi(2);
      TTree *treePRF2d = fCalibraVector->ConvertVectorPTreeHisto(fCalibraVector->GetVectorPRF(),fCalibraVector->GetPlaPRF(),"treePRF2d",(const char *) name);
      fout->WriteTObject(treePRF2d);
    }
  }
  
  fout->Close();
  
  AliInfo(Form("Execution time Write2d: R:%.2fs C:%.2fs"
	      ,stopwatch.RealTime(),stopwatch.CpuTime()));

  return kTRUE;
  
}

//____________Probe the histos_________________________________________________
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
  // [6] : mean relativ error
  //

  Double_t *info = new Double_t[7];
   
  // Number of Xbins (detectors or groups of pads)
  Int_t    nbins   = h->GetNbinsX(); //number of calibration groups
  Int_t    nybins  = h->GetNbinsY(); //number of bins per histo

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
  TH1F *NbEntries = 0x0;//distribution of the number of entries
  TH1F *NbEntriesPerGroup = 0x0;//Number of entries per group
  TProfile *NbEntriesPerSp = 0x0;//Number of entries for one supermodule
    
  // Beginning of the loop over the calibration groups 
  for (Int_t idect = 0; idect < nbins; idect++) {

    TH1I *projch = (TH1I *) h->ProjectionY("projch",idect+1,idect+1,(Option_t *)"e");
    projch->SetDirectory(0);
    
    // Number of entries for this calibration group
    Double_t nentries = 0.0;
    if((i%2) == 0){
      for (Int_t k = 0; k < nybins; k++) {
	nentries += h->GetBinContent(h->GetBin(idect+1,k+1));
      }
    }
    else{
      for (Int_t k = 0; k < nybins; k++) {
	nentries += ((TProfile2D *)h)->GetBinEntries(h->GetBin(idect+1,k+1));
	if(h->GetBinContent(h->GetBin(idect+1,k+1)) != 0) {
	  meanrelativerror += (h->GetBinError(h->GetBin(idect+1,k+1))
                            / (TMath::Abs(h->GetBinContent(h->GetBin(idect+1,k+1)))));
	  counter++;
	} 
      }
    }

    //Debug
    if(i > 1){
      if((!((Bool_t)NbEntries)) && (nentries > 0)){
	NbEntries = new TH1F("Number of entries","Number of entries"
                               ,100,(Int_t)nentries/2,nentries*2);
	NbEntries->SetDirectory(0);
	NbEntriesPerGroup = new TH1F("Number of entries per group","Number of entries per group"
                               ,nbins,0,nbins);
	NbEntriesPerGroup->SetDirectory(0);
	NbEntriesPerSp = new TProfile("Number of entries per supermodule","Number of entries per supermodule"
                               ,(Int_t)(nbins/18),0,(Int_t)(nbins/18));
	NbEntriesPerSp->SetDirectory(0);
      }
      if(NbEntries){
	if(nentries > 0) NbEntries->Fill(nentries);
	NbEntriesPerGroup->Fill(idect+0.5,nentries);
	NbEntriesPerSp->Fill((idect%((Int_t)(nbins/18)))+0.5,nentries);
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
    NbEntries->Draw("");
    stat->cd(2);
    NbEntriesPerSp->SetStats(0);
    NbEntriesPerSp->Draw("");
    TCanvas *stat1 = new TCanvas("stat1","",50,50,600,800);
    stat1->cd();
    NbEntriesPerGroup->SetStats(0);
    NbEntriesPerGroup->Draw("");
  }

  return info;

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
  // Create the 2D histos
  //

  TString name("Nz");
  name += fCalibraMode->GetNz(2);
  name += "Nrphi";
  name += fCalibraMode->GetNrphi(2);

  fPRF2d = new TProfile2D("PRF2d",(const Char_t *) name
                                 ,nn,0,nn,fNumberBinPRF,-1.5,1.5);
  fPRF2d->SetXTitle("Det/pad groups");
  fPRF2d->SetYTitle("Position x/W [pad width units]");
  fPRF2d->SetZTitle("Q_{i}/Q_{total}");
  fPRF2d->SetStats(0);

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
                               ,nn,0,nn,fTimeMax
                               ,-0.5/fSf,(Float_t) (fTimeMax-0.5)/fSf);
  fPH2d->SetXTitle("Det/pad groups");
  fPH2d->SetYTitle("time [#mus]");
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
                         ,nn,0,nn,fNumberBinCharge,0,300);
  fCH2d->SetXTitle("Det/pad groups");
  fCH2d->SetYTitle("charge deposit [a.u]");
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
	
  Int_t nb =  0; // Nombre de zones traversees
  Int_t fd = -1; // Premiere zone non nulle
  
  
  // See if the track goes through different zones
  for (Int_t k = 0; k < fCalibraMode->GetNfragZ(0)*fCalibraMode->GetNfragRphi(0); k++) {
    if (fAmpTotal[k] > 0.0) {
      nb++;
      if (nb == 1) {
        fd = k;
      }
    }
  }
 
  // If automatic scale
  if ((fCountRelativeScale < 100) && (fRelativeScaleAuto)) {
    // Take only the one zone track
    if (nb == 1) {
      fRelativeScale += fAmpTotal[fd] * 0.014 * 0.01;
      fCountRelativeScale++;
    }
  }

  // We fill the CH2d after having scale with the first 100
  if ((fCountRelativeScale >= 100) && (fRelativeScaleAuto)) {
    // Case of track with only one zone
    if (nb == 1) {
      if (fHisto2d) {
        fCH2d->Fill(fCalibraMode->GetXbins(0)+fd+0.5,fAmpTotal[fd]/fRelativeScale);
      }
      if (fVector2d) {
        fCalibraVector->UpdateVectorCH(fCalibraMode->GetXbins(0)+fd,fAmpTotal[fd]/fRelativeScale);
      }
    } // Case 1 zone
    // Case of track with two zones
    if (nb == 2) {
      // Two zones voisines sinon rien!
      if ((fAmpTotal[fd]   > 0.0) && 
          (fAmpTotal[fd+1] > 0.0)) {
	// One of the two very big
	if (fAmpTotal[fd] > fProcent*fAmpTotal[fd+1]) {
	  if (fHisto2d) {
            fCH2d->Fill(fCalibraMode->GetXbins(0)+fd+0.5,fAmpTotal[fd]/fRelativeScale);
	  }
	  if (fVector2d) {
            fCalibraVector->UpdateVectorCH(fCalibraMode->GetXbins(0)+fd,fAmpTotal[fd]/fRelativeScale);
	  }
	}
	if (fAmpTotal[fd+1] > fProcent*fAmpTotal[fd])  {
	  if (fHisto2d) {
            fCH2d->Fill(fCalibraMode->GetXbins(0)+fd+1.5,fAmpTotal[fd+1]/fRelativeScale);
	  }
	  if (fVector2d) {
            fCalibraVector->UpdateVectorCH(fCalibraMode->GetXbins(0)+fd,fAmpTotal[fd+1]/fRelativeScale);
	  }
	}
      }
    } // Case 2 zones
  }

  // Fill with no automatic scale
  if (!fRelativeScaleAuto) {
    // Case of track with only one zone
    if (nb == 1) {
      fNumberUsedCh[0]++;
      if (fHisto2d) {
        fCH2d->Fill(fCalibraMode->GetXbins(0)+fd+0.5,fAmpTotal[fd]/fRelativeScale);
      }
      if (fVector2d) {
        fCalibraVector->UpdateVectorCH(fCalibraMode->GetXbins(0)+fd,fAmpTotal[fd]/fRelativeScale);
      }
    } // Case 1 zone
    // Case of track with two zones
    if (nb == 2) {
      // Two zones voisines sinon rien!
      // Case 1
      if ((fAmpTotal[fd]   > 0.0) && 
          (fAmpTotal[fd+1] > 0.0)) {
	// One of the two very big
	if (fAmpTotal[fd] > fProcent*fAmpTotal[fd+1]) {
	  if (fHisto2d) {
            fCH2d->Fill(fCalibraMode->GetXbins(0)+fd+0.5,fAmpTotal[fd]/fRelativeScale);
	  }
	  if (fVector2d) {
            fCalibraVector->UpdateVectorCH(fCalibraMode->GetXbins(0)+fd,fAmpTotal[fd]/fRelativeScale);
	  }
	  fNumberUsedCh[1]++;
	}
	if (fAmpTotal[fd+1] > fProcent*fAmpTotal[fd]) {
	  if (fHisto2d) {
            fCH2d->Fill(fCalibraMode->GetXbins(0)+fd+1.5,fAmpTotal[fd+1]/fRelativeScale);
	  }
	  if (fVector2d) {
            fCalibraVector->UpdateVectorCH(fCalibraMode->GetXbins(0)+fd+1,fAmpTotal[fd+1]/fRelativeScale);
	  }
	  fNumberUsedCh[1]++;
	}
      }
      // Case 2
      if (fCalibraMode->GetNfragZ(0) > 1) {
	if (fAmpTotal[fd] > 0.0) {
	  if ((fd+fCalibraMode->GetNfragZ(0)) < (fCalibraMode->GetNfragZ(0)*fCalibraMode->GetNfragRphi(0))) {
	    if (fAmpTotal[fd+fCalibraMode->GetNfragZ(0)] > 0.0) {
	      // One of the two very big
	      if (fAmpTotal[fd] > fProcent*fAmpTotal[fd+fCalibraMode->GetNfragZ(0)]) {
		if (fHisto2d) {
                  fCH2d->Fill(fCalibraMode->GetXbins(0)+fd+0.5,fAmpTotal[fd]/fRelativeScale);
		}
		if (fVector2d) {
                  fCalibraVector->UpdateVectorCH(fCalibraMode->GetXbins(0)+fd,fAmpTotal[fd]/fRelativeScale);
		}
		fNumberUsedCh[1]++;
	      }
	      if (fAmpTotal[fd+fCalibraMode->GetNfragZ(0)] > fProcent*fAmpTotal[fd]) {
		if (fHisto2d) {
                  fCH2d->Fill(fCalibraMode->GetXbins(0)+fd+fCalibraMode->GetNfragZ(0)
                            + 0.5,fAmpTotal[fd+fCalibraMode->GetNfragZ(0)]/fRelativeScale);
		}
		fNumberUsedCh[1]++;
		if (fVector2d) {
                  fCalibraVector->UpdateVectorCH(fCalibraMode->GetXbins(0)+fd+fCalibraMode->GetNfragZ(0)
                                                ,fAmpTotal[fd+fCalibraMode->GetNfragZ(0)]/fRelativeScale);
		}
	      }
	    }
	  }
	}
      }
    } // Case 2 zones

  }

}

//____________Offine tracking in the AliTRDtracker_____________________________
void AliTRDCalibraFillHisto::ResetfVariables()
{
  //
  // Reset values of fAmpTotal, fPHValue and fPHPlace for
  // the updateHistogram... functions
  //

  // Reset the good track
  fGoodTrack = kTRUE;
  
  // Reset the fAmpTotal where we put value
  if (fCH2dOn) {
    for (Int_t k = 0; k < fCalibraMode->GetNfragZ(0)*fCalibraMode->GetNfragRphi(0); k++) {
      fAmpTotal[k] = 0.0;
    }
  }
  
  // Reset the fPHValue
  if (fPH2dOn) {
    for (Int_t k = 0; k < fTimeMax; k++) {
      fPHValue[k] = 0.0;
      fPHPlace[k] = -1;
    }
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
  
  // Fill 
  // Case of track with only one zone
  if (nb == 1) {
    fNumberUsedPh[0]++;
    //fd1 is the only zone
    for (Int_t i = 0; i < fTimeMax; i++) {
      if (fHisto2d) {
	fPH2d->Fill((fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
      }
      if (fVector2d) {
	fCalibraVector->UpdateVectorPH((fCalibraMode->GetXbins(1)+fd1),i,fPHValue[i]);
      }
    }
  } // Case 1 zone
  // Case of track with two zones
  if (nb == 2) {
    // Two zones voisines sinon rien!
    // Case 1
    if ((fd1 == fd2+1) || 
        (fd2 == fd1+1)) {
      // One of the two fast all the think
      if (k2 > (k1+fDifference)) {
	//we choose to fill the fd1 with all the values
	fNumberUsedPh[1]++;
	for (Int_t i = 0; i < fTimeMax; i++) {
	  if (fHisto2d) {
	    fPH2d->Fill((fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
	  }
	  if (fVector2d) {
	    fCalibraVector->UpdateVectorPH((fCalibraMode->GetXbins(1)+fd1),i,fPHValue[i]);
	  }
	}
      }
      if ((k2+fDifference) < fTimeMax) {
	//we choose to fill the fd2 with all the values
	fNumberUsedPh[1]++;
	for (Int_t i = 0; i < fTimeMax; i++) {
	  if (fHisto2d) {
	    fPH2d->Fill((fCalibraMode->GetXbins(1)+fd2)+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
	  }
	  if (fVector2d) {
	    fCalibraVector->UpdateVectorPH((fCalibraMode->GetXbins(1)+fd2),i,fPHValue[i]);
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
		fPH2d->Fill((fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
	      }
	      if (fVector2d) {
		fCalibraVector->UpdateVectorPH((fCalibraMode->GetXbins(1)+fd1),i,fPHValue[i]);
	      }
	    }
	  }
	  if ((k2+fDifference) < fTimeMax) {
	    //we choose to fill the fd2 with all the values
	    fNumberUsedPh[1]++;
	    for (Int_t i = 0; i < fTimeMax; i++) {
	      if (fHisto2d) {
		fPH2d->Fill((fCalibraMode->GetXbins(1)+fd2)+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
	      }
	      if (fVector2d) {
		fCalibraVector->UpdateVectorPH((fCalibraMode->GetXbins(1)+fd2),i,fPHValue[i]);
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
		fPH2d->Fill((fCalibraMode->GetXbins(1)+fd1)+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
	      }
	      if (fVector2d) {
		fCalibraVector->UpdateVectorPH((fCalibraMode->GetXbins(1)+fd1),i,fPHValue[i]);
	      }
	    }
	  }
	  if ((k2+fDifference) < fTimeMax) {
	    //we choose to fill the fd2 with all the values
	    fNumberUsedPh[1]++;
	    for (Int_t i = 0; i < fTimeMax; i++) {
	      if (fHisto2d) {
		fPH2d->Fill((fCalibraMode->GetXbins(1)+fd2)+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
	      }
	      if (fVector2d) {
		fCalibraVector->UpdateVectorPH((fCalibraMode->GetXbins(1)+fd2),i,fPHValue[i]);
	      }
	    }
	  }
	}
      }
    }

  } // case 2 zones

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
