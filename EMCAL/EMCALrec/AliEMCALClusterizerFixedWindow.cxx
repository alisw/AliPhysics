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

// This class derives from AliEMCALClustrerizer

// --- Root ---
#include <TMath.h> 
#include <TMinuit.h>
#include <TTree.h> 
#include <TBenchmark.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TH1I.h>

// --- AliRoot ---
#include "AliLog.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALCalibData.h"
#include "AliEMCALCalibTime.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALUnfolding.h"

#include "AliEMCALClusterizerFixedWindow.h"

ClassImp(AliEMCALClusterizerFixedWindow)

//__________________________________________________________________________________________
AliEMCALClusterizerFixedWindow::AliEMCALClusterizerFixedWindow() :
  AliEMCALClusterizer(), 
  fNphi(4), 
  fNeta(4), 
  fShiftPhi(2), 
  fShiftEta(2),
  fTRUshift(0),
  fNEtaDigitsSupMod(0),
  fNPhiDigitsSupMod(0),
  fNTRUPhi(0),
  fNTRUEta(0),
  fNEtaDigits(0),
  fNPhiDigits(0),
  fMaxShiftPhi(0),
  fMaxShiftEta(0),
  fNDigitsCluster(0),
  fNClusEtaNoShift(0),
  fNClusPhiNoShift(0),
  fNClusters(0),
  fNTotalClus(0),
  fClustersArray(0),
  fInitialized(0)
{
  // Constructor
}

//__________________________________________________________________________________________
AliEMCALClusterizerFixedWindow::AliEMCALClusterizerFixedWindow(AliEMCALGeometry* geometry) :
  AliEMCALClusterizer(geometry), 
  fNphi(4), 
  fNeta(4), 
  fShiftPhi(2), 
  fShiftEta(2),
  fTRUshift(0),
  fNEtaDigitsSupMod(0),
  fNPhiDigitsSupMod(0),
  fNTRUPhi(0),
  fNTRUEta(0),
  fNEtaDigits(0),
  fNPhiDigits(0),
  fMaxShiftPhi(0),
  fMaxShiftEta(0),
  fNDigitsCluster(0),
  fNClusEtaNoShift(0),
  fNClusPhiNoShift(0),
  fNClusters(0),
  fNTotalClus(0),
  fClustersArray(0),
  fInitialized(0)
{
  // Constructor
}

//__________________________________________________________________________________________
AliEMCALClusterizerFixedWindow::AliEMCALClusterizerFixedWindow(AliEMCALGeometry* geometry, 
                                                               AliEMCALCalibData * calib, 
                                                               AliEMCALCalibTime * calibt, 
                                                               AliCaloCalibPedestal * caloped) :
  AliEMCALClusterizer(geometry, calib, calibt, caloped),
  fNphi(4), 
  fNeta(4), 
  fShiftPhi(2), 
  fShiftEta(2),
  fTRUshift(0),
  fNEtaDigitsSupMod(0),
  fNPhiDigitsSupMod(0),
  fNTRUPhi(0),
  fNTRUEta(0),
  fNEtaDigits(0),
  fNPhiDigits(0),
  fMaxShiftPhi(0),
  fMaxShiftEta(0),
  fNDigitsCluster(0),
  fNClusEtaNoShift(0),
  fNClusPhiNoShift(0),
  fNClusters(0),
  fNTotalClus(0),
  fClustersArray(0),
  fInitialized(0)
{
  // Constructor
}

//__________________________________________________________________________________________
AliEMCALClusterizerFixedWindow::~AliEMCALClusterizerFixedWindow()
{
  // Destructor

  if (fClustersArray) {
    for (Int_t i = 0; i < fNTotalClus; i++) {
      if (fClustersArray[i]) {
	delete[] fClustersArray[i];
	fClustersArray[i] = 0;
      }
    }
    delete[] fClustersArray;
    fClustersArray = 0;
  }

}

//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::SetNphi (Int_t n) 
{
  // Set fNphi; if clusterizer already initialized gives a warning and does nothing
  
  if (fInitialized != 0)
    AliWarning("Clusterizer already initialized. Unable to change the parameters.");
  else
    fNphi = n;
}

//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::SetNeta (Int_t n) 
{
  // Set fNeta; if clusterizer already initialized gives a warning and does nothing
  
  if (fInitialized != 0)
    AliWarning("Clusterizer already initialized. Unable to change the parameters.");
  else
    fNeta = n;
}

//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::SetShiftPhi (Int_t s) 
{
  // Set fShiftPhi; if clusterizer already initialized gives a warning and does nothing
  
  if (fInitialized != 0)
    AliWarning("Clusterizer already initialized. Unable to change the parameters.");
  else
    fShiftPhi = s;
}

//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::SetShiftEta (Int_t s) 
{
  // Set fShiftEta; if clusterizer already initialized gives a warning and does nothing
  
  if (fInitialized != 0)
    AliWarning("Clusterizer already initialized. Unable to change the parameters.");
  else
    fShiftEta = s;
}

//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::SetTRUshift(Bool_t b) 
{
  // Set fTRUshift; if clusterizer already initialized gives a warning and does nothing
  
  if (fInitialized != 0)
    AliWarning("Clusterizer already initialized. Unable to change the parameters.");
  else
    fTRUshift = b;
}


//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::Digits2Clusters(Option_t * option)
{
  // Steering method to perform clusterization for the current event 
 
  static Float_t cputime = 0;
  static Float_t realtime = 0;

  if (strstr(option,"tim"))
    gBenchmark->Start("EMCALClusterizer"); 
	
  if (strstr(option,"print"))
    Print(""); 
	
  //Get calibration parameters from file or digitizer default values.
  GetCalibrationParameters();
	
  //Get dead channel map from file or digitizer default values.
  GetCaloCalibPedestal();
	
  MakeClusters();  //only the real clusters
	
  if (fToUnfold) {
    fClusterUnfolding->SetInput(fNumberOfECAClusters,fRecPoints,fDigitsArr);
    fClusterUnfolding->MakeUnfolding();
  }
	
  //Evaluate position, dispersion and other RecPoint properties for EC section 
  for (Int_t index = 0; index < fRecPoints->GetEntries(); index++) { 
    AliEMCALRecPoint * rp = dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(index));
    if (rp) {
      rp->EvalAll(fECAW0,fDigitsArr,fJustClusters);
      AliDebug(5, Form("MAX INDEX %d ", rp->GetMaximalEnergyIndex()));
      //For each rec.point set the distance to the nearest bad crystal
      if (fCaloPed)
        rp->EvalDistanceToBadChannels(fCaloPed);
    }
  }
  
  fRecPoints->Sort();
	
  for (Int_t index = 0; index < fRecPoints->GetEntries(); index++) {
    AliEMCALRecPoint *rp = dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(index));
    if (rp) {
      rp->SetIndexInList(index);
    }
    else AliFatal("RecPoint NULL!!");
  }
	
  if (fTreeR)
    fTreeR->Fill();
	
  if (strstr(option,"deb") || strstr(option,"all"))  
    PrintRecPoints(option);
	
  AliDebug(1,Form("EMCAL Clusterizer found %d Rec Points",fRecPoints->GetEntriesFast()));
	
  if (strstr(option,"tim")) {
    gBenchmark->Stop("EMCALClusterizer");
    Printf("Exec took %f CPU time (%f real time) for clusterizing", 
           gBenchmark->GetCpuTime("EMCALClusterizer")-cputime,gBenchmark->GetRealTime("EMCALClusterizer")-realtime);
    cputime = gBenchmark->GetCpuTime("EMCALClusterizer");
    realtime = gBenchmark->GetRealTime("EMCALClusterizer");
  }    
}

//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::ExecOnce()
{
  // Initialize clusterizer.
  
  fInitialized = -1;

  if (!fGeom) {
    AliError("Did not get geometry!");
    return;
  }
	
  // Defining geometry and clusterization parameter
  fNEtaDigitsSupMod = fGeom->GetNEta() * fGeom->GetNETAdiv(); // always 48?;
  fNPhiDigitsSupMod = fGeom->GetNPhi() * fGeom->GetNPHIdiv(); // always 24?;
  
  fNTRUPhi = 1;
  fNTRUEta = 1;
  
  fNEtaDigits = fNEtaDigitsSupMod * fGeom->GetNumberOfSuperModules() / fGeom->GetNPhiSuperModule();
  fNPhiDigits = fNPhiDigitsSupMod * fGeom->GetNPhiSuperModule();    
  
  if (fTRUshift){
    fNTRUPhi = fGeom->GetNPhiSuperModule() * 3;
    fNTRUEta = fGeom->GetNumberOfSuperModules() / fGeom->GetNPhiSuperModule();
    fNEtaDigits /= fNTRUEta;
    fNPhiDigits /= fNTRUPhi;
  }

  // Check if clusterizer parameter are compatible with calorimeter geometry
  if (fNEtaDigits < fNeta){
    AliError(Form("Error: fNeta = %d is greater than nEtaDigits = %d.",fNeta,fNEtaDigits));
    return;
  }
  if (fNPhiDigits < fNphi){
    AliError(Form("Error: fNphi = %d is greater than nPhiDigits = %d.",fNphi,fNPhiDigits));
    return;
  }
  if (fNEtaDigits % fShiftEta != 0){
    AliError(Form("Error: fShiftEta = %d is such that clusters cannot slide the whole calorimeter (nEtaDigits = %d).",fShiftEta,fNEtaDigits));
    return;
  }
  if (fNPhiDigits % fShiftPhi != 0){
    AliError(Form("Error: fShiftPhi = %d is such that clusters cannot slide the whole calorimeter (nPhiDigits = %d).",fShiftPhi,fNPhiDigits));
    return;
  }
  if (fNeta % fShiftEta != 0){
    AliError(Form("Error: fShiftEta = %d is not divisor of fNeta = %d.",fShiftEta,fNeta));
    return;
  }
  if (fNphi % fShiftPhi != 0){
    AliError(Form("Error: fShiftPhi = %d is not divisor of fNphi = %d).",fShiftPhi,fNphi));
    return;
  }
  
  fMaxShiftPhi = fNphi / fShiftPhi;
  fMaxShiftEta = fNeta / fShiftEta;
  
  fNClusEtaNoShift = fNEtaDigits / fNeta;
  fNClusPhiNoShift = fNPhiDigits / fNphi;

  fNClusters = fNClusEtaNoShift * fNClusPhiNoShift * fNTRUEta * fNTRUPhi;
 
  fNTotalClus = fNClusters * fMaxShiftEta * fMaxShiftPhi;

  fNDigitsCluster = fNphi * fNeta;

  if (fClustersArray) {
    for (Int_t i = 0; i < fNTotalClus; i++) {
      if (fClustersArray[i]) {
	delete[] fClustersArray[i];
	fClustersArray[i] = 0;
      }
    }
    delete[] fClustersArray;
    fClustersArray = 0;
  }

  fClustersArray = new AliEMCALDigit**[fNTotalClus];
  for (Int_t i = 0; i < fNTotalClus; i++) {
    fClustersArray[i] = new AliEMCALDigit*[fNDigitsCluster];
    for (Int_t j = 0; j < fNDigitsCluster; j++) {
      fClustersArray[i][j] = 0;
    }
  }

  AliDebug(1,Form("****ExecOnce*****\n"
		  "fNphi = %d, fNeta = %d, fShiftPhi = %d, fShiftEta = %d, fTRUshift = %d\n"
		  "fNEtaDigitsSupMod = %d, fNPhiDigitsSupMod = %d, fNTRUPhi = %d, fNTRUEta = %d, fNEtaDigits = %d, fNPhiDigits = %d\n"
		  "fMaxShiftPhi = %d, fMaxShiftEta = %d, fNDigitsCluster = %d, fNClusEtaNoShift = %d, fNClusPhiNoShift = %d\n"
		  "fNClusters = %d, fNTotalClus = %d\n",
		  fNphi,fNeta,fShiftPhi,fShiftEta,fTRUshift,
		  fNEtaDigitsSupMod,fNPhiDigitsSupMod,fNTRUPhi,fNTRUEta,fNEtaDigits,fNPhiDigits,
		  fMaxShiftPhi,fMaxShiftEta,fNDigitsCluster,fNClusEtaNoShift,fNClusPhiNoShift,
		  fNClusters,fNTotalClus));

  fInitialized = 1;
}
  
//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::MakeClusters()
{
  // Make clusters.

  fNumberOfECAClusters = 0;
  fRecPoints->Delete();

  if (fInitialized == 0)
    ExecOnce();

  if (fInitialized == -1) {
    AliError(Form("%s: error initializing the clusterizer. No clusterization will be performed.",GetName()));
    return;
  }
  
  // Set up TObjArray with pointers to digits to work on calibrated digits 
  TObjArray *digitsC = new TObjArray();
  AliEMCALDigit *digit;
  Float_t dEnergyCalibrated = 0.0, ehs = 0.0, time = 0.0;
  TIter nextdigit(fDigitsArr);
  while ((digit = dynamic_cast<AliEMCALDigit*>(nextdigit()))) { // calibrate and clean up digits
    dEnergyCalibrated =  digit->GetAmplitude();
    time              =  digit->GetTime();
    Calibrate(dEnergyCalibrated, time, digit->GetId());
    digit->SetCalibAmp(dEnergyCalibrated);
    digit->SetTime(time);
    if (dEnergyCalibrated < fMinECut || time > fTimeMax || time < fTimeMin) {
      continue;
    }
    else if (!fGeom->CheckAbsCellId(digit->GetId())) {
      continue;
    }
    else {
      ehs += dEnergyCalibrated;
      digitsC->AddLast(digit);
    }
  } 
  
  AliDebug(1,Form("MakeClusters: Number of digits %d  -> (e %f), ehs %f\n",
                  fDigitsArr->GetEntries(),fMinECut,ehs));
   
  Int_t nSupMod=0, nModule=0, nIphi=0, nIeta=0;
  Int_t iphi=0, ieta=0;  // cell eta-phi indexes in SM

  for (Int_t ishiftPhi = 0; ishiftPhi < fMaxShiftPhi; ishiftPhi++){
    Int_t nClusPhi = (fNPhiDigits - fShiftPhi * ishiftPhi) / fNphi;
    
    for (Int_t ishiftEta = 0; ishiftEta < fMaxShiftEta; ishiftEta++) {
      
      Int_t nClusEta = (fNEtaDigits - fShiftEta * ishiftEta) / fNeta; 
      
      Int_t iTotalClus = fNClusters * (ishiftPhi * fMaxShiftEta + ishiftEta);
      
      TIter nextdigitC(digitsC);
      while ((digit = dynamic_cast<AliEMCALDigit*>(nextdigitC()))) { // scan over the list of digitsC
        
        fGeom->GetCellIndex (digit->GetId(), nSupMod, nModule, nIphi, nIeta);
        fGeom->GetCellPhiEtaIndexInSModule (nSupMod, nModule, nIphi, nIeta, iphi, ieta);
        
        Int_t iphi_eff = iphi - fShiftPhi * ishiftPhi + fNPhiDigitsSupMod * (nSupMod / 2); // N supermodules along phi
        
        Int_t iTRUphi = iphi_eff / fNPhiDigits;
        
        iphi_eff -= iTRUphi * fNPhiDigits;
        
        Int_t iClusPhi = iphi_eff / fNphi; 
        
        if (iphi_eff < 0 || iClusPhi >= nClusPhi) 
          continue;
        
        Int_t ieta_eff = ieta - fShiftEta * ishiftEta + fNEtaDigitsSupMod * (nSupMod % 2); // 2 supermodules along eta
        
        Int_t iTRUeta = ieta_eff / fNEtaDigits;
        
        ieta_eff -= iTRUeta * fNEtaDigits;
        
        Int_t iClusEta = ieta_eff / fNeta; 
        
        if (ieta_eff < 0 || iClusEta >= nClusEta) 
          continue;
        
        iphi_eff += iTRUphi * fNPhiDigits;
        iClusPhi = iphi_eff / fNphi; 
        
        ieta_eff += iTRUeta * fNEtaDigits;
        iClusEta = ieta_eff / fNeta; 
        
        Int_t iCluster = iClusPhi + iClusEta * fNClusPhiNoShift * fNTRUPhi; 
        Int_t iDigit = iphi_eff % fNphi + (ieta_eff % fNeta) * fNphi;

        if (iCluster >= fNClusters){
          AliError(Form("iCluster out of range! iCluster = %d, fNClusters = %d (should never happen...)", iCluster, fNClusters));
          return;
        }
        
        iCluster += iTotalClus;

        if (iCluster >= fNTotalClus){
          AliError(Form("iCluster out of range! iCluster = %d, fNTotalClus = %d (should never happen...)", iCluster, fNTotalClus));
          return;
        }

        if (iDigit >= fNDigitsCluster){
          AliError(Form("iDigit out of range! iDigit = %d, fNDigitsCluster = %d (should never happen...)", iDigit, fNDigitsCluster));
          return;
        }

        if (fClustersArray[iCluster][iDigit] != 0){
          AliError("Digit already added! (should never happen...)");
          return;
        }
        
        fClustersArray[iCluster][iDigit] = digit;
        
      } // loop on digit
      
    } // loop on eta shift
    
  } // loop on phi shift
  
  AliEMCALRecPoint *recPoint = 0;
  Bool_t recPointOk = kFALSE;
  for (Int_t iCluster = 0; iCluster < fNTotalClus; iCluster++) {

    if (!recPoint) {
      if(fNumberOfECAClusters >= fRecPoints->GetSize()) fRecPoints->Expand(fNumberOfECAClusters*2+1);

      (*fRecPoints)[fNumberOfECAClusters] = new AliEMCALRecPoint("");
      recPoint = static_cast<AliEMCALRecPoint*>(fRecPoints->At(fNumberOfECAClusters));
    }
		
    if (recPoint) {
      recPoint->SetClusterType(AliVCluster::kEMCALClusterv1);
      recPoint->SetUniqueID(iCluster);
      fNumberOfECAClusters++;

      for (Int_t iDigit = 0; iDigit < fNDigitsCluster; iDigit++) {
        if (fClustersArray[iCluster][iDigit] == 0) continue;
        digit = fClustersArray[iCluster][iDigit];
        recPoint->AddDigit(*digit, digit->GetCalibAmp(), kFALSE); //Time or TimeR?
        fClustersArray[iCluster][iDigit] = 0;
	recPointOk = kTRUE;
      }

      if (recPointOk) { // unset the pointer so that a new rec point will be allocated in the next iteration
	recPoint = 0;
	recPointOk = kFALSE;
      }
    }
    else {
      AliError("Error allocating rec points!");
      break;
    }
  }

  if (!recPointOk) {
    fRecPoints->RemoveLast();
    fNumberOfECAClusters--;
  }

  delete digitsC;
  AliDebug(1, Form("MakeClusters: Number of digits %d  -> (e %f)\n", fDigitsArr->GetEntries(),fMinECut));
  AliDebug(1, Form("total no of clusters %d from %d digits", fNumberOfECAClusters, fDigitsArr->GetEntriesFast())); 
}
