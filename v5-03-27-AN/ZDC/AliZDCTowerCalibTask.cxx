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

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TDecompLU.h" 

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliZDCTowerCalibTask.h"
#include "AliZDCTowerCalib.h"
#include "AliCDBManager.h"
#include "AliCDBId.h"
#include "AliCDBRunRange.h"
#include "AliCDBStorage.h"

ClassImp(AliZDCTowerCalibTask)

//________________________________________________________________________
AliZDCTowerCalibTask::AliZDCTowerCalibTask() 
  : AliAnalysisTask(), 
    fESD(0),
    fAZNA(),
    fAZNC(),
    fBZNA(),
    fBZNC(),
    fADCMin(0)
{
}

//________________________________________________________________________
AliZDCTowerCalibTask::AliZDCTowerCalibTask(const char *name) 
  : AliAnalysisTask(name, ""), 
    fESD(0),
    fAZNA(),
    fAZNC(),
    fBZNA(),
    fBZNC(),
    fADCMin(0)
{
  // Default constructor
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
//   // Output slot #0 writes into a TH1 container
//   DefineOutput(0, TH1F::Class());
  
  fAZNA.ResizeTo(4,4); 
  fAZNC.ResizeTo(4,4); 
  fBZNA.ResizeTo(4); 
  fBZNC.ResizeTo(4); 
  
  fAZNA.Zero(); 
  fAZNC.Zero(); 
  fBZNA.Zero(); 
  fBZNC.Zero(); 
  fADCMin = 0; 
}

//________________________________________________________________________
void AliZDCTowerCalibTask::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  
  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) Printf("ERROR: Could not read chain from input slot 0");
  else {
    // Disable all branches and enable only the needed ones
    // The next two lines are different when data produced as AliESDEvent is read
    /*
      tree->SetBranchStatus("*", kFALSE);
      tree->SetBranchStatus("fTracks.*", kTRUE);
    */
    
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!esdH) Printf("ERROR: Could not get ESDInputHandler");
    else fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliZDCTowerCalibTask::CreateOutputObjects(){
  // Create histograms
  // Called once
}

//________________________________________________________________________
void AliZDCTowerCalibTask::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  AliESDZDC *zdcEvent = fESD->GetESDZDC();
  Double_t *ezn1 = 0, *ezn2 = 0; 
  if(fESD->GetEventType()==7){
    ezn1 = (Double_t *) zdcEvent->GetZN1TowerEnergy(); // znc
    ezn2 = (Double_t *) zdcEvent->GetZN2TowerEnergy(); // zna
  }
  Double_t ezn1sum = ezn1[1] + ezn1[2] + ezn1[3] + ezn1[4]; 
  Double_t ezn2sum = ezn2[1] + ezn2[2] + ezn2[3] + ezn2[4]; 
  // select event
  Bool_t isZNAAccepted = kTRUE; 
  if (ezn2[0] < fADCMin || ezn2sum < fADCMin) isZNAAccepted = kFALSE; 
  Bool_t isZNCAccepted = kTRUE; 
  if (ezn1[0] < fADCMin || ezn1sum < fADCMin) isZNCAccepted = kFALSE; 
  // calculate coefficient matrix and known terms 
  for (Int_t i=0; i<4; i++) { 
    for (Int_t j=0; j<4; j++) { 
      if (isZNAAccepted) fAZNA[i][j] += ezn2[ i + 1 ] * ezn2[ j + 1 ] / ezn2[0]; 
      if (isZNCAccepted) fAZNC[i][j] += ezn1[ i + 1 ] * ezn1[ j + 1 ] / ezn1[0]; 
    }
    if (isZNAAccepted) fBZNA[i] += ezn2[ i + 1 ]; 
    if (isZNCAccepted) fBZNC[i] += ezn1[ i + 1 ];     
  }
  // Post output data.
  //  PostData(0, fHistPt);
}      

//________________________________________________________________________
void AliZDCTowerCalibTask::Terminate(Option_t *) 
{
  // solve the system of linear equations giving the calibration coefficients 
  // for zna and znc

  TDecompLU luZNA(fAZNA);
  Bool_t isZNAOk;
  TVectorD coeffZNA = luZNA.Solve(fBZNA,isZNAOk);
  if (isZNAOk) { 
    printf ("ZNA: coefficient matrix:\n"); 
    fAZNA.Print(); 
    printf ("ZNA: known terms:\n"); 
    fBZNA.Print(); 
    printf ("ZNA: Fitted calibration coefficients:\n"); 
    coeffZNA.Print(); 
  }
  else { 
    printf ("Singular coefficient matrix for ZNA! Setting all coefficients to 1\n"); 
    for (Int_t i=0; i<4; i++) coeffZNA[i] = 1; 
  } 

  TDecompLU luZNC(fAZNC);
  Bool_t isZNCOk;
  TVectorD coeffZNC = luZNC.Solve(fBZNC,isZNCOk);
  if (isZNCOk) { 
    printf ("ZNC: coefficient matrix:\n"); 
    fAZNC.Print(); 
    printf ("ZNC: known terms:\n"); 
    fBZNC.Print(); 
    printf ("ZNC: Fitted calibration coefficients:\n"); 
    coeffZNC.Print(); 
  }
  else {
    printf ("Singular coefficient matrix for ZNC! Setting all coefficients to 1\n"); 
    for (Int_t i=0; i<4; i++) coeffZNC[i] = 1; 
  } 
  // write to OCDB



  AliZDCTowerCalib *towerCalib = new AliZDCTowerCalib();
  
  towerCalib->SetZN1EqualCoeff(0, 1.);
  towerCalib->SetZN2EqualCoeff(0, 1.);
  
  for(Int_t j=1; j<5; j++){  
    towerCalib->SetZN1EqualCoeff(j, coeffZNC[j+1]);
    towerCalib->SetZN2EqualCoeff(j, coeffZNA[j+1]);
  }
  towerCalib->Print("");
  
  AliCDBManager *manager = AliCDBManager::Instance();
  manager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Oppedisano");
  md->SetComment("Calibration object for ZDC written by a macro");
  md->SetObjectClassName("AliZDCTowerCalib");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  
  AliCDBId id("ZDC/Calib/TowerCalib",fESD->GetRunNumber(),AliCDBRunRange::Infinity());
  manager->GetDefaultStorage()->Put(towerCalib,id, md);
}
