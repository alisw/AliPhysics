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
// AliHLTTRDCalibra                                                               
//                                                                             
// This class is for the TRD calibration of the relative gain factor, the drift velocity,
// the time 0 and the pad response function.        
// It can be used for the calibration per chamber but also per group of pads and eventually per pad.
// The user has to choose with the functions SetNz and SetNrphi the precision of the calibration. 
//Begin_Html
/*
<br>
<CENTER>
<TABLE border=1>
<TR><TD><center>Nz</center></TD><TD><center> 0 </center></TD><TD><center> 1 </center></TD><TD><center> 2 </center></TD><TD><center> 3 </center></TD><TD><center> 4 </center></TD></TR>
<TR><TD><CENTER>group of row pads per detector</CENTER></TD><TD><CENTER>1</CENTER></TD><TD><CENTER>2</CENTER></TD><TD><CENTER>4</CENTER></TD><TD><CENTER>6(chamb2)<br> 8(others chambers)</CENTER></TD><TD><CENTER>12 (chamb2)<br> 16 (chamb0)</CENTER></TD></TR>
<TR><TD><CENTER>row pads per group</CENTER></TD><TD><CENTER>12 (chamb2)<br> 16 (chamb0)</CENTER></TD><TD><CENTER>6 (chamb2)<br> 8 (chamb0)</CENTER></TD><TD><CENTER>3 (chamb2)<br> 4 (chamb0)</CENTER></TD><TD><CENTER>2</CENTER></TD><TD><CENTER>1</CENTER></TD></TR>
<TR><TD><CENTER>~distance [cm]</CENTER></TD><TD><CENTER>106 (chamb2)<br> 130 (chamb0)</CENTER></TD><TD><CENTER>53 (chamb2)<br> 65 (chamb0)</CENTER></TD><TD><CENTER>26.5 (chamb2)<br> 32.5 (chamb0)</CENTER></TD><TD><CENTER>17 (chamb2)<br> 17 (chamb0)</CENTER></TD><TD><CENTER>9 (chamb2)<br> 9 (chamb0)</CENTER></TD></TR>
<CAPTION>In the z direction</CAPTION>
</TABLE>
</CENTER>
<CENTER>
<br>
<TABLE border=1>
<TR><TD><center>Nrphi</center></TD><TD><center> 0 </center></TD><TD><center> 1 </center></TD><TD><center> 2 </center></TD><TD><center> 3 </center></TD><TD><center> 4 </center></TD><TD><center> 5 </center></TD><TD><center> 6 </center></TD></TR>
<TR><TD><CENTER>group of col pads per detector</CENTER></TD><TD><CENTER>1</CENTER></TD><TD><CENTER>2</CENTER></TD><TD><CENTER>4</CENTER></TD><TD><CENTER>8</CENTER></TD><TD><CENTER>16</CENTER></TD><TD><center>36</center></TD><TD><center>144</center></TD></TR>
<TR><TD><CENTER>col pads per group</CENTER></TD><TD><CENTER>144</CENTER></TD><TD><CENTER>72</CENTER></TD><TD><CENTER>36</CENTER></TD><TD><CENTER>18</CENTER></TD><TD><CENTER>9</CENTER></TD><TD><center>4</center></TD><TD><center>1</center></TD></TR>
<TR><TD><CENTER>~distance [cm]</CENTER></TD><TD><CENTER>113.4</CENTER></TD><TD><CENTER>56.7</CENTER></TD><TD><CENTER>25.3</CENTER></TD><TD><CENTER>14.3</CENTER></TD><TD><CENTER>7.25</CENTER></TD><TD><center>3.2</center></TD><TD><center>0.8</center></TD></TR>
<CAPTION>In the rphi direction</CAPTION>
</TABLE>
</CENTER>
<br>
*/
//End_Html 
//
// Fill histograms or vectors
//----------------------------
//   
// 2D Histograms (Histo2d) or vectors (Vector2d), then converted in Trees, will be filled
// from RAW DATA in a run or from reconstructed TRD tracks during the offline tracking 
// in the function "FollowBackProlongation" (AliTRDtracker)
// Per default the functions to fill are off.                                   
//
//Begin_Html
/*
Example of 2D histos for the relative gain (left) and the drift velocity (right) calibration of the sector 13 <br>
<center>
<img src="./gif/2dhisto.gif" width="600" height="350"><br>
</center>
*/
//End_Html    
// Author:
//   R. Bailhache (R.Bailhache@gsi.de)
//                            
//////////////////////////////////////////////////////////////////////////////////////

#include <TH1I.h>
#include <TProfile2D.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TROOT.h>

#include "AliLog.h"
#include "AliCDBManager.h"

#include "AliHLTTRDCalibra.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "AliTRDmcmTracklet.h"


ClassImp(AliHLTTRDCalibra)

AliHLTTRDCalibra* AliHLTTRDCalibra::fgInstance = 0;
Bool_t AliHLTTRDCalibra::fgTerminated = kFALSE;

//_____________singleton implementation_________________________________________________
AliHLTTRDCalibra *AliHLTTRDCalibra::Instance()
{
  //
  // Singleton implementation
  //

  if (fgTerminated != kFALSE) {
    return 0;
  }

  if (fgInstance == 0) {
    fgInstance = new AliHLTTRDCalibra();
  }

  return fgInstance;

}

//______________________________________________________________________________________
void AliHLTTRDCalibra::Terminate()
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
AliHLTTRDCalibra::AliHLTTRDCalibra()
  :TObject()
  ,fOn(kFALSE)
  ,fRelativeScale(0)
  ,fThresholdClusterPRF2(0.0)
  ,fWriteName(0)
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
  ,fPH2d(0x0)
  ,fPRF2d(0x0)
  ,fCH2d(0x0)
{
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 3; i++) {
    fNz[i]    = 0;
    fNrphi[i] = 0;
  }

  for (Int_t k = 0; k < 3; k++) {
    fNtotal[k]    = 0;
    fDetChamb2[k] = 0;
    fDetChamb0[k] = 0;
  }

  // Init
  Init();
  
}

//______________________________________________________________________________________
AliHLTTRDCalibra::AliHLTTRDCalibra(const AliHLTTRDCalibra &c)
  :TObject(c)
  ,fOn(kFALSE)
  ,fRelativeScale(0)
  ,fThresholdClusterPRF2(0.0)
  ,fWriteName(0)
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
  ,fPH2d(0x0)
  ,fPRF2d(0x0)
  ,fCH2d(0x0)
{
  //
  // Copy constructor
  //

}

//____________________________________________________________________________________
AliHLTTRDCalibra::~AliHLTTRDCalibra()
{
  //
  // AliHLTTRDCalibra destructor
  //

  ClearHistos();
  
}

//_____________________________________________________________________________
void AliHLTTRDCalibra::Destroy() 
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
void AliHLTTRDCalibra::ClearHistos() 
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
void AliHLTTRDCalibra::Init() 
{
  //
  // Init some default values
  //

  // How to fill the 2D
  fRelativeScale        = 20.0;
  fThresholdClusterPRF2 = 3.0;
  
  // Store the Info
  fNumberBinCharge      = 100;
  fNumberBinPRF         = 20;
  
  // Write
  fWriteName            = "TRD.calibration.root";
     
  // Internal variables
  
  // Fill the 2D histos in the offline tracking
  fGoodTrack             = kTRUE;

  fProcent               = 6.0;
  fDifference            = 17;
  fNumberClusters        = 18;
  fNumberTrack           = 0;
  fNumberUsedCh[0]       = 0;
  fNumberUsedCh[1]       = 0;
  fNumberUsedPh[0]       = 0;
  fNumberUsedPh[1]       = 0;
  
  // Pad calibration
  for (Int_t i = 0; i < 3; i++) {
    fRowMin[i]    = -1;
    fRowMax[i]    = -1;
    fColMax[i]    = -1;
    fColMin[i]    = -1;
    fNnZ[i]       = -1;
    fNnRphi[i]    = -1;
    fNfragZ[i]    = -1;
    fNfragRphi[i] = -1;
    fXbins[i]     = -1;
  }
  
}
//____________Functions for initialising the AliHLTTRDCalibra in the code_________
Bool_t AliHLTTRDCalibra::Init2Dhistos()
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
 
  // Create the 2D histos corresponding to the pad groupCalibration mode
  AliInfo(Form("The pad calibration mode for the relative gain calibration: Nz %d, and Nrphi %d"
	       ,fNz[0]
	       ,fNrphi[0]));
  
  // Calcul the number of Xbins
  fNtotal[0] = 0;
  ModePadCalibration(2,0);
  ModePadFragmentation(0,2,0,0);
  fDetChamb2[0] = fNfragZ[0] * fNfragRphi[0];
  fNtotal[0] += 6 * 18 * fDetChamb2[0];
  ModePadCalibration(0,0);
  ModePadFragmentation(0,0,0,0);
  fDetChamb0[0] = fNfragZ[0] * fNfragRphi[0];
  fNtotal[0] += 6 * 4 * 18 * fDetChamb0[0];
  AliInfo(Form("Total number of Xbins: %d",fNtotal[0]));
  
  CreateCH2d(fNtotal[0]);
    
  // Variable
  fAmpTotal = new Float_t[TMath::Max(fDetChamb2[0],fDetChamb0[0])];
  for (Int_t k = 0; k < TMath::Max(fDetChamb2[0],fDetChamb0[0]); k++) {
    fAmpTotal[k] = 0.0;
  } 
  
  
  AliInfo(Form("The pad calibration mode for the drift velocity calibration: Nz %d, and Nrphi %d"
	       ,fNz[1]
	       ,fNrphi[1]));
  
  // Calcul the number of Xbins
  fNtotal[1] = 0;
  ModePadCalibration(2,1);
  ModePadFragmentation(0,2,0,1);
  fDetChamb2[1] = fNfragZ[1]*fNfragRphi[1];
  fNtotal[1] += 6 * 18 * fDetChamb2[1];
  ModePadCalibration(0,1);
  ModePadFragmentation(0,0,0,1);
  fDetChamb0[1] = fNfragZ[1] * fNfragRphi[1];
  fNtotal[1] += 6 * 4 * 18 * fDetChamb0[1];
  AliInfo(Form("Total number of Xbins: %d",fNtotal[1]));
  
  // Create the 2D histo
  CreatePH2d(fNtotal[1]);
  
  // Variable
  fPHPlace = new Short_t[fTimeMax];
  for (Int_t k = 0; k < fTimeMax; k++) {
    fPHPlace[k] = -1;
  } 
  fPHValue = new Float_t[fTimeMax];
  for (Int_t k = 0; k < fTimeMax; k++) {
    fPHValue[k] = -1.0;
  }
  
  AliInfo(Form("The pad calibration mode for the PRF calibration: Nz %d, and Nrphi %d"
	       ,fNz[2]
	       ,fNrphi[2]));
  
  // Calcul the number of Xbins
  fNtotal[2] = 0;
  ModePadCalibration(2,2);
  ModePadFragmentation(0,2,0,2);
  fDetChamb2[2] = fNfragZ[2] * fNfragRphi[2];
  fNtotal[2] += 6 * 18 * fDetChamb2[2];
  ModePadCalibration(0,2);
  ModePadFragmentation(0,0,0,2);
  fDetChamb0[2] = fNfragZ[2] * fNfragRphi[2];
  fNtotal[2] += 6 * 4 * 18 * fDetChamb0[2];
  AliInfo(Form("Total number of Xbins: %d",fNtotal[2]));
  
  // Create the 2D histo
  CreatePRF2d(fNtotal[2]);
  
  return kTRUE;
  
}

//____________Functions for filling the histos in the code_____________________

//____________Online trackling in AliTRDtrigger________________________________
Bool_t AliHLTTRDCalibra::UpdateHistogramcm(AliTRDmcmTracklet *trk)
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
  if (fNnZ[0] != 0) {
    posr[0] = (Int_t) row / fNnZ[0];
  }
  if (fNnZ[1] != 0) {
    posr[1] = (Int_t) row / fNnZ[1];
  }
  if (fNnZ[2] != 0) {
    posr[2] = (Int_t) row / fNnZ[2];
  }
 
  // Eventuelle correction due to track angle in z direction
  Float_t correction = 1.0;
  Float_t z = trk->GetRowz();
  Float_t r = trk->GetTime0();
  correction = r / TMath::Sqrt((r*r+z*z));
  
  //Count the tracks
  fNumberTrack++;

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
      if (fNnRphi[0] != 0) {
        posc[0] = (Int_t) col / fNnRphi[0];
      }
      if (fNnRphi[1] != 0) {
        posc[1] = (Int_t) col / fNnRphi[1];
      }
      if (fNnRphi[2] != 0) {
        posc[2] = (Int_t) col / fNnRphi[2];
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
      fPHPlace[time] = posc[1] * fNfragZ[1] + posr[1];
      fAmpTotal[(Int_t) (posc[0]*fNfragZ[0]+posr[0])] += (Float_t) (amp[0]+amp[1]+amp[2]);
      fPHValue[time] = (Float_t) (amp[0]+amp[1]+amp[2]);
            
      
      // Fill PRF direct
      if (good) {
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
	    if ((xcenter > -0.5) && 
                (xcenter <  0.5)) {
	      Float_t yminus = amp[0] / (amp[0]+amp[1]+amp[2]);
	      Float_t ymax   = amp[2] / (amp[0]+amp[1]+amp[2]);
	      // Fill only if it is in the drift region!
	      if (((Float_t) time / fSf) > 0.3) {
		fPRF2d->Fill((fXbins[2]+posc[2]*fNfragZ[2]+posr[2]+0.5),xcenter,ycenter);
		if (xcenter < 0.0) {
		  fPRF2d->Fill((fXbins[2]+posc[2]*fNfragZ[2]+posr[2]+0.5),-(xcenter+1.0),yminus);
		}
		if (xcenter > 0.0) {
		  fPRF2d->Fill((fXbins[2]+posc[2]*fNfragZ[2]+posr[2]+0.5),(1.0-xcenter),ymax);
		}
	      } 
	    }
	  }
	}
      }
      
    } // Boucle clusters
    
    // Fill the charge and PH
    if (fGoodTrack) {
      FillTheInfoOfTheTrackCH();
      FillTheInfoOfTheTrackPH();
    }
    
  } // Condition on number of clusters
  
  return kTRUE;
  
}

//____________Functions for seeing if the pad is really okey___________________

//_____________________________________________________________________________
Bool_t AliHLTTRDCalibra::IsPadOn(Int_t detector, Int_t col, Int_t row) const
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
//____________Writing the 2D___________________________________________________

//_____________________________________________________________________________
Bool_t AliHLTTRDCalibra::Write2d()
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

  fout->WriteTObject(fCH2d);
  fout->WriteTObject(fPH2d);
  fout->WriteTObject(fPRF2d);
  
  fout->Close();
  
  AliInfo(Form("Execution time Write2d: R:%.2fs C:%.2fs"
	      ,stopwatch.RealTime(),stopwatch.CpuTime()));

  return kTRUE;
  
}
//_____________________________________________________________________________
void AliHLTTRDCalibra::SetRelativeScale(Float_t RelativeScale)
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
void AliHLTTRDCalibra::SetNz(Int_t i, Short_t Nz)
{
  //
  // Set the mode of calibration group in the z direction for the parameter i
  // 

  if ((Nz >= 0) && 
      (Nz <  5)) {
    fNz[i] = Nz; 
  }
  else { 
    AliInfo("You have to choose between 0 and 4");
  }

}

//_____________________________________________________________________________
void AliHLTTRDCalibra::SetNrphi(Int_t i, Short_t Nrphi)
{
  //
  // Set the mode of calibration group in the rphi direction for the parameter i
  //
 
  if ((Nrphi >= 0) && 
      (Nrphi <  7)) {
    fNrphi[i] = Nrphi; 
  }
  else {
    AliInfo("You have to choose between 0 and 6");
  }

}
//____________Pad Calibration Public___________________________________________

//____________Define the number of pads per group for one detector and one calibration
void AliHLTTRDCalibra::ModePadCalibration(Int_t iChamb, Int_t i)
{
  //
  // Definition of the calibration mode
  // from Nz and Nrphi, the number of row and col pads per calibration groups are setted
  //


  fNnZ[i]    = 0;
  fNnRphi[i] = 0;
  
  if ((fNz[i] == 0) && (iChamb == 2)) {
    fNnZ[i] = 12;
  }
  if ((fNz[i] == 0) && (iChamb != 2)) {
    fNnZ[i] = 16;
  }  
  if ((fNz[i] == 1) && (iChamb == 2)) {
    fNnZ[i] = 6;
  }
  if ((fNz[i] == 1) && (iChamb != 2)) {
    fNnZ[i] = 8;
  }
  if ((fNz[i] == 2) && (iChamb == 2)) {
    fNnZ[i] = 3;
  }
  if ((fNz[i] == 2) && (iChamb != 2)) {
    fNnZ[i] = 4;
  }
  if (fNz[i] == 3) {
    fNnZ[i] = 2;
  }
  if (fNz[i] == 4) {
    fNnZ[i] = 1;
  }
   
  if (fNrphi[i] == 0) {
    fNnRphi[i] = 144;
  }
  if (fNrphi[i] == 1) {
    fNnRphi[i] = 72;
  } 
  if (fNrphi[i] == 2) {
    fNnRphi[i] = 36;
  } 
  if (fNrphi[i] == 3) {
    fNnRphi[i] = 18;
  } 
  if (fNrphi[i] == 4) {
    fNnRphi[i] = 9;
  } 
  if (fNrphi[i] == 5) {
    fNnRphi[i] = 4;
  } 
  if (fNrphi[i] == 6) {
    fNnRphi[i] = 1;
  } 

}

//____________Define the number of pad groups in one detector for one calibration
Bool_t AliHLTTRDCalibra::ModePadFragmentation(Int_t iPlane,Int_t iChamb, Int_t iSect, Int_t i)
{
  //
  // Definition of the calibration mode
  // From the number of row and col pads per calibration groups the
  // number of calibration groups are setted
  //

  fNfragZ[i]    = 0;
  fNfragRphi[i] = 0;
  
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }

  // A little geometry:
  Int_t rowMax = parCom->GetRowMax(iPlane,iChamb,iSect);
  Int_t colMax = parCom->GetColMax(iPlane);
  
  // The fragmentation
  if (fNnZ[i]    != 0) {
    fNfragZ[i]    = (Int_t) rowMax / fNnZ[i];
  }

  if (fNnRphi[i] != 0) {
    fNfragRphi[i] = (Int_t) colMax / fNnRphi[i];
  }

  return kTRUE;

}

//____________Protected Functions______________________________________________
//____________Create the 2D histo to be filled online__________________________
//

//_____________________________________________________________________________
void AliHLTTRDCalibra::CreatePRF2d(Int_t nn)
{
  //
  // Create the 2D histos
  //

  TString name("Nz");
  name += fNz[2];
  name += "Nrphi";
  name += fNrphi[2];

  fPRF2d = new TProfile2D("PRF2d",(const Char_t *) name
                                 ,nn,0,nn,fNumberBinPRF,-1.0,1.0);
  fPRF2d->SetXTitle("Det/pad groups");
  fPRF2d->SetYTitle("Position x/W [pad width units]");
  fPRF2d->SetZTitle("Q_{i}/Q_{total}");
  fPRF2d->SetStats(0);

}

//_____________________________________________________________________________
void AliHLTTRDCalibra::CreatePH2d(Int_t nn)
{
  //
  // Create the 2D histos
  //

  TString name("Nz");
  name += fNz[1];
  name += "Nrphi";
  name += fNrphi[1];

  fPH2d = new TProfile2D("PH2d",(const Char_t *) name
                               ,nn,0,nn,fTimeMax
                               ,-0.5/fSf,(Float_t) (fTimeMax-0.5)/fSf);
  fPH2d->SetXTitle("Det/pad groups");
  fPH2d->SetYTitle("time [#mus]");
  fPH2d->SetZTitle("<PH> [a.u.]");
  fPH2d->SetStats(0);

}

//_____________________________________________________________________________
void AliHLTTRDCalibra::CreateCH2d(Int_t nn)
{
  //
  // Create the 2D histos
  //

  TString name("Nz");
  name += fNz[0];
  name += "Nrphi";
  name += fNrphi[0];

  fCH2d = new TH2I("CH2d",(const Char_t *) name
                         ,nn,0,nn,fNumberBinCharge,0,300);
  fCH2d->SetXTitle("Det/pad groups");
  fCH2d->SetYTitle("charge deposit [a.u]");
  fCH2d->SetZTitle("counts");
  fCH2d->SetStats(0);
  fCH2d->Sumw2();

}

//____________Offine tracking in the AliTRDtracker_____________________________
void AliHLTTRDCalibra::FillTheInfoOfTheTrackCH()
{
  //
  // For the offline tracking or mcm tracklets
  // This function will be called in the functions UpdateHistogram... 
  // to fill the info of a track for the relativ gain calibration
  //
	
  Int_t nb =  0; // Nombre de zones traversees
  Int_t fd = -1; // Premiere zone non nulle
  
  
  // See if the track goes through different zones
  for (Int_t k = 0; k < fNfragZ[0]*fNfragRphi[0]; k++) {
    if (fAmpTotal[k] > 0.0) {
      nb++;
      if (nb == 1) {
        fd = k;
      }
    }
  }
  
  // Case of track with only one zone
  if (nb == 1) {
    fNumberUsedCh[0]++;
    fCH2d->Fill(fXbins[0]+fd+0.5,fAmpTotal[fd]/fRelativeScale);
  } // Case 1 zone
    // Case of track with two zones
  if (nb == 2) {
    // Two zones voisines sinon rien!
    // Case 1
    if ((fAmpTotal[fd]   > 0.0) && 
	(fAmpTotal[fd+1] > 0.0)) {
      // One of the two very big
      if (fAmpTotal[fd] > fProcent*fAmpTotal[fd+1]) {
	fCH2d->Fill(fXbins[0]+fd+0.5,fAmpTotal[fd]/fRelativeScale);
	fNumberUsedCh[1]++;
      }
      if (fAmpTotal[fd+1] > fProcent*fAmpTotal[fd]) {
	fCH2d->Fill(fXbins[0]+fd+1.5,fAmpTotal[fd+1]/fRelativeScale);
	fNumberUsedCh[1]++;
      }
    }
    // Case 2
    if (fNfragZ[0] > 1) {
      if (fAmpTotal[fd] > 0.0) {
	if ((fd+fNfragZ[0]) < (fNfragZ[0]*fNfragRphi[0])) {
	  if (fAmpTotal[fd+fNfragZ[0]] > 0.0) {
	    // One of the two very big
	    if (fAmpTotal[fd] > fProcent*fAmpTotal[fd+fNfragZ[0]]) {
	      fCH2d->Fill(fXbins[0]+fd+0.5,fAmpTotal[fd]/fRelativeScale);
	      fNumberUsedCh[1]++;
	    }
	    if (fAmpTotal[fd+fNfragZ[0]] > fProcent*fAmpTotal[fd]) {
	      fCH2d->Fill(fXbins[0]+fd+fNfragZ[0]+0.5,fAmpTotal[fd+fNfragZ[0]]/fRelativeScale);
	      fNumberUsedCh[1]++;
	    }
	  }
	}
      }
    }
  } // Case 2 zones
  
}

//____________Offine tracking in the AliTRDtracker_____________________________
void AliHLTTRDCalibra::ResetfVariables()
{
  //
  // Reset values of fAmpTotal, fPHValue and fPHPlace for
  // the updateHistogram... functions
  //

  // Reset the good track
  fGoodTrack = kTRUE;
  
  // Reset the fAmpTotal where we put value
  for (Int_t k = 0; k < fNfragZ[0]*fNfragRphi[0]; k++) {
    fAmpTotal[k] = 0.0;
  }
  
  // Reset the fPHValue
  for (Int_t k = 0; k < fTimeMax; k++) {
    fPHValue[k] = -1.0;
    fPHPlace[k] = -1;
  }

}

//____________Offine tracking in the AliTRDtracker_____________________________
void AliHLTTRDCalibra::FillTheInfoOfTheTrackPH()
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
    for (Int_t i = 0; i < fTimeMax; i++) {
      if (fPHValue[i] > 0.0) {
	fPH2d->Fill((fXbins[1]+fPHPlace[i])+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
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
	fNumberUsedPh[1]++;
	for (Int_t i = k1; i < k2; i++) {
	  if (fPHValue[i] > 0.0) {
	    fPH2d->Fill((fXbins[1]+fPHPlace[i])+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
	  }
	}
      }
      if ((k2+fDifference) < fTimeMax) {
	fNumberUsedPh[1]++;
	for (Int_t i = k2; i < fTimeMax; i++) {
	  if (fPHValue[i] > 0.0) {
	    fPH2d->Fill((fXbins[1]+fPHPlace[i])+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
	  }
	}
      }
    }
    // Two zones voisines sinon rien!
    if (fNfragZ[1] > 1) {
      // Case 2
      if ((fd1+fNfragZ[1]) < (fNfragZ[1]*fNfragRphi[1])) {
	if (fd2 == (fd1+fNfragZ[1])) {
	  // One of the two fast all the think
	  if (k2 > (k1+fDifference)) {
	    fNumberUsedPh[1]++;
	    for (Int_t i = k1; i < k2; i++) {
	      if (fPHValue[i] > 0.0) {
		fPH2d->Fill((fXbins[1]+fPHPlace[i])+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
	      }
	    }
	  }
	  if ((k2+fDifference) < fTimeMax) {
	    fNumberUsedPh[1]++;
	    for (Int_t i = k2; i < fTimeMax; i++) {
	      if (fPHValue[i] > 0.0) {
		fPH2d->Fill((fXbins[1]+fPHPlace[i])+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
	      }
	    }
	  }
	}
      }
      // Two zones voisines sinon rien!
      // Case 3
      if ((fd1 - fNfragZ[1]) >= 0) {
	if (fd2 == (fd1 - fNfragZ[1])) {
	  // One of the two fast all the think
	  if (k2 > (k1 + fDifference)) {
	    fNumberUsedPh[1]++;
	    for (Int_t i = k1; i < k2; i++) {
	      if (fPHValue[i] > 0.0) {
		fPH2d->Fill((fXbins[1]+fPHPlace[i])+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
	      }
	    }
	  }
	  if ((k2+fDifference) < fTimeMax) {
	    fNumberUsedPh[1]++;
	    for (Int_t i = k2; i < fTimeMax; i++) {
	      if (fPHValue[i] > 0.0) {
		fPH2d->Fill((fXbins[1]+fPHPlace[i])+0.5,(Float_t) i/fSf,(Float_t) fPHValue[i]);
	      }
	    }
	  }
	}
      }
    }

  } // case 2 zones

}

//____________Set the pad calibration variables for the detector_______________
Bool_t AliHLTTRDCalibra::LocalisationDetectorXbins(Int_t detector)
{
  //
  // For the detector calcul the first Xbins and set the number of row
  // and col pads per calibration groups, the number of calibration
  // groups in the detector.
  //
  
  // first Xbins of the detector
  CalculXBins(detector,0);
  CalculXBins(detector,1);
  CalculXBins(detector,2);
  

  // fragmentation of idect
  for (Int_t i = 0; i < 3; i++) {
    ModePadCalibration((Int_t) GetChamber(detector),i);
    ModePadFragmentation((Int_t) GetPlane(detector)
			 , (Int_t) GetChamber(detector)
			 , (Int_t) GetSector(detector),i);
  }
  
  return kTRUE;

}
//____________Pad group calibration mode_______________________________________
//

//_____________________________________________________________________________
void AliHLTTRDCalibra::ReconstructionRowPadGroup(Int_t idect, Int_t i)
{
  //
  // For the calibration group idect in a detector calculate the
  // first and last row pad and col pad.
  // The pads in the interval will have the same calibrated coefficients
  //

  Int_t posc = -1;
  Int_t posr = -1;
  fRowMin[i] = -1;
  fRowMax[i] = -1;
  fColMin[i] = -1;
  fColMax[i] = -1;
  
  if (fNfragZ[i]    != 0) {
    posc = (Int_t) idect / fNfragZ[i];
  }
  if (fNfragRphi[i] != 0) {
    posr = (Int_t) idect % fNfragZ[i];
  }
  fRowMin[i] = posr     * fNnZ[i];
  fRowMax[i] = (posr+1) * fNnZ[i];
  fColMin[i] = posc     * fNnRphi[i];
  fColMax[i] = (posc+1) * fNnRphi[i];

}

//_____________________________________________________________________________
void AliHLTTRDCalibra::CalculXBins(Int_t idect, Int_t i)
{
  //
  // For the detector idect calcul the first Xbins
  //

  fXbins[i] = 0;
 
  // In which sector?
  Int_t sector = GetSector(idect);
  fXbins[i] += sector*(6*fDetChamb2[i]+6*4*fDetChamb0[i]);
 
  // In which chamber?
  Int_t chamber = GetChamber(idect);
  Int_t kc      = 0;
  while (kc < chamber) {
    if (kc == 2) {
      fXbins[i] += 6 * fDetChamb2[i];
    }
    else {
      fXbins[i] += 6 * fDetChamb0[i];
    }
    kc ++;
  }
  
  // In which plane?
  Int_t plane = GetPlane(idect);
  if (chamber == 2) {
    fXbins[i] += plane*fDetChamb2[i];
  }
  else {
    fXbins[i] += plane*fDetChamb0[i];
  }
 
}
//
//____________Some basic geometry function_____________________________________
//

//_____________________________________________________________________________
Int_t AliHLTTRDCalibra::GetPlane(Int_t d) const
{
  //
  // Reconstruct the plane number from the detector number
  //

  return ((Int_t) (d % 6));

}

//_____________________________________________________________________________
Int_t AliHLTTRDCalibra::GetChamber(Int_t d) const
{
  //
  // Reconstruct the chamber number from the detector number
  //
  Int_t fgkNplan = 6;

  return ((Int_t) (d % 30) / fgkNplan);

}

//_____________________________________________________________________________
Int_t AliHLTTRDCalibra::GetSector(Int_t d) const
{
  //
  // Reconstruct the sector number from the detector number
  //
  Int_t fg = 30;

  return ((Int_t) (d / fg));

}
