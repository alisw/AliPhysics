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

//_________________________________________________________________________
//  A singleton. This class should be used in the analysis stage to get 
//  reconstructed objects: Digits, RecPoints, TrackSegments and RecParticles,
//  instead of directly reading them from galice.root file. 
//                  
//  MvL Feb 2006:
//  The AliEMCALLoader now holds the TClonesArray and TObjArray for reading
//  Hits, Dgits, SDigits and RecPoints. Filling is managed in the GetEvent()
//  method.
//
//  Creation/writing of files is managed by the relevant parts of the 
//  reconstruction software (AliEMCALDigitiser etx)
//
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)
//*--         Completely redesigned by Dmitri Peressounko March 2001  
//
//*-- YS June 2001 : renamed the original AliEMCALIndexToObject and make
//*--         systematic usage of TFolders without changing the interface
// 
//*-- Marco van Leeuwen, Jan 2006: complete revision to simplify reading
//*--         and fit better in general ALICE scheme
//
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

#include "TTree.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliEMCALLoader.h"
#include "AliLog.h"
#include "AliCDBLocal.h"
#include "AliCDBStorage.h"
#include "AliCDBManager.h"

ClassImp(AliEMCALLoader)
  
const TString AliEMCALLoader::fgkECARecPointsBranchName("EMCALECARP");//Name for branch with ECA Reconstructed Points
AliEMCALCalibData* AliEMCALLoader::fCalibData = 0; //calibation data
AliEMCALAlignData* AliEMCALLoader::fAlignData = 0; //alignment data

//____________________________________________________________________________ 
AliEMCALLoader::AliEMCALLoader()
{
  fDebug = 0;
  fHits = new TClonesArray("AliEMCALHit");
  fDigits = new TClonesArray("AliEMCALDigit");
  fSDigits = new TClonesArray("AliEMCALDigit");
  fRecPoints = new TObjArray();
}

//____________________________________________________________________________ 
AliEMCALLoader::AliEMCALLoader(const Char_t *detname,const Char_t *eventfoldername):
  AliLoader(detname,eventfoldername)
{
  fDebug=0;
  fHits = new TClonesArray("AliEMCALHit");
  fDigits = new TClonesArray("AliEMCALDigit");
  fSDigits = new TClonesArray("AliEMCALDigit");
  fRecPoints = new TObjArray();
}

//____________________________________________________________________________ 
AliEMCALLoader::~AliEMCALLoader()
{
  // Disconnect trees and remove arrays
  if (TreeH())
    TreeH()->SetBranchAddress(fDetectorName,0);
  if (TreeD())
    TreeD()->SetBranchAddress(fDetectorName,0);
  if (TreeS())
    TreeS()->SetBranchAddress(fDetectorName,0);
  if (TreeR())
    TreeR()->SetBranchAddress(fgkECARecPointsBranchName,0);
  delete fHits;
  delete fDigits;
  delete fSDigits;
  delete fRecPoints;
}

//____________________________________________________________________________
AliEMCALAlignData* AliEMCALLoader::AlignData()
{
  // Check if the instance of AliEMCALAlignData exists, and return it
  
  if( !(AliCDBManager::Instance()->IsDefaultStorageSet()) )
    fAlignData=0x0;
  return fAlignData;
}

//____________________________________________________________________________ 
AliEMCALCalibData* AliEMCALLoader::CalibData()
{ 
  // Check if the instance of AliEMCALCalibData exists, and return it

  if( !(AliCDBManager::Instance()->IsDefaultStorageSet()) ) 
    fCalibData=0x0;
  
  return fCalibData;
  
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::CalibrateRaw(Double_t energy, Int_t module, 
				   Int_t column, Int_t row)
{
  // Convert energy into digitized amplitude for a cell relId
  // It is a user responsilibity to open CDB and set
  // AliEMCALCalibData object by the following operators:
  // 
  // AliCDBLocal *loc = new AliCDBLocal("deCalibDB");
  // AliEMCALCalibData* clb = (AliEMCALCalibData*)AliCDBStorage::Instance()
  //    ->Get(path_to_calibdata,run_number);
  // AliEMCALGetter* gime = AliEMCALGetter::Instance("galice.root");
  // gime->SetCalibData(clb);

  if (CalibData() == 0)
    Warning("CalibrateRaw","Calibration DB is not initiated!");

  Float_t gainFactor = 0.00305;//0.0015; // width of one ADC channel in GeV
  Float_t pedestal   = 0.009;//0.005;  // pedestals

  if(CalibData()) {
    gainFactor = CalibData()->GetADCchannel (module,column,row);
    pedestal   = CalibData()->GetADCpedestal(module,column,row);
  }
  
  Int_t   amp = static_cast<Int_t>( (energy - pedestal) / gainFactor + 0.5 ) ; 
  return amp;
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::LoadHits(Option_t* opt) {
  
  Int_t status = AliLoader::LoadHits(opt);  // First call AliLoader to do all the groundwork
  
  TTree *treeH = TreeH();
  
  if (treeH) {
    treeH->SetBranchAddress(fDetectorName,&fHits);
    if (treeH->GetEntries() > 1)
      AliWarning("Multiple arrays in treeH no longer supported");
    treeH->GetEvent(0);
  }
 
  return status;
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::LoadSDigits(Option_t* opt) {
  
  Int_t status = AliLoader::LoadSDigits(opt);  // First call AliLoader to do all the groundwork
  
  TTree *treeS = TreeS();
  
  if (treeS) {
    treeS->SetBranchAddress(fDetectorName,&fSDigits);
    treeS->GetEvent(0);
  }
 
  return status;
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::LoadDigits(Option_t* opt) {
  
  Int_t status = AliLoader::LoadDigits(opt);  // First call AliLoader to do all the groundwork
  
  TTree *treeD = TreeD();
  
  if (treeD) {
    treeD->SetBranchAddress(fDetectorName,&fDigits);
    treeD->GetEvent(0);
  }
 
  return status;
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::LoadRecPoints(Option_t* opt) {
  
  Int_t status = AliLoader::LoadRecPoints(opt);  // First call AliLoader to do all the groundwork
  
  TTree *treeR = TreeR();
  if (treeR) {
    treeR->SetBranchAddress(fgkECARecPointsBranchName,&fRecPoints);
    treeR->GetEvent(0);
  }
  
  return status;
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::GetEvent() {
  
  AliLoader::GetEvent();  // First call AliLoader to do all the groundwork
  
  // Now connect and fill TClonesArray

  // Hits
   TTree *treeH = TreeH();
   
   if (treeH) {
     treeH->SetBranchAddress(fDetectorName,&fHits);
     if (treeH->GetEntries() > 1)
       AliWarning("Multiple arrays in treeH no longer supported");
     treeH->GetEvent(0);
   }
   
   
   // SDigits
   TTree *treeS = TreeS();
   if (treeS) {
     treeS->SetBranchAddress(fDetectorName,&fSDigits);
     treeS->GetEvent(0);
   }
   
   // Digits
   TTree *treeD = TreeD();
   if (treeD) {
     treeD->SetBranchAddress(fDetectorName,&fDigits);
     treeD->GetEvent(0);
   }

   // RecPoints
   TTree *treeR = TreeR();
   if (treeR) {
     treeR->SetBranchAddress(fgkECARecPointsBranchName,&fRecPoints);
     treeR->GetEvent(0);
   }
   
   return 0;
}
