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

ClassImp(AliEMCALLoader)
  
const TString AliEMCALLoader::fgkECARecPointsBranchName("EMCALECARP");//Name for branch with ECA Reconstructed Points

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
