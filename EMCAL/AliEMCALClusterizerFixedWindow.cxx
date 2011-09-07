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
#include "AliESDCaloCluster.h"
#include "AliEMCALUnfolding.h"

#include "AliEMCALClusterizerFixedWindow.h"

ClassImp(AliEMCALClusterizerFixedWindow)

//__________________________________________________________________________________________
  AliEMCALClusterizerFixedWindow::AliEMCALClusterizerFixedWindow()
    : AliEMCALClusterizer(), 
      fNphi(4), 
      fNeta(4), 
      fShiftPhi(2), 
      fShiftEta(2),
      fTRUshift(0),
      fClustersArray(0)
{
  // Constructor
}

//__________________________________________________________________________________________
AliEMCALClusterizerFixedWindow::AliEMCALClusterizerFixedWindow(AliEMCALGeometry* geometry)
  : AliEMCALClusterizer(geometry), 
    fNphi(4), 
    fNeta(4), 
    fShiftPhi(2), 
    fShiftEta(2),
    fTRUshift(0),
    fClustersArray(0)
{
  // Constructor
}

//__________________________________________________________________________________________
AliEMCALClusterizerFixedWindow::AliEMCALClusterizerFixedWindow(AliEMCALGeometry* geometry, AliEMCALCalibData * calib, AliCaloCalibPedestal * caloped)
  : AliEMCALClusterizer(geometry, calib, caloped),
    fNphi(4), 
    fNeta(4), 
    fShiftPhi(2), 
    fShiftEta(2),
    fTRUshift(0),
    fClustersArray(0)
{
  // Constructor
}

//__________________________________________________________________________________________
AliEMCALClusterizerFixedWindow::~AliEMCALClusterizerFixedWindow()
{
  // Destructor
  
  delete fClustersArray;
}

//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::SetNphi (Int_t n) 
{
  // Set fNphi; if clusterizer already initialized gives a warning and does nothing
  
  if (fClustersArray)
    AliWarning("Clusterizer already initialized. Unable to change the parameters.");
  else
    fNphi = n;
}

//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::SetNeta (Int_t n) 
{
  // Set fNeta; if clusterizer already initialized gives a warning and does nothing
  
  if (fClustersArray)
    AliWarning("Clusterizer already initialized. Unable to change the parameters.");
  else
    fNeta = n;
}

//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::SetShiftPhi (Int_t s) 
{
  // Set fShiftPhi; if clusterizer already initialized gives a warning and does nothing
  
  if (fClustersArray)
    AliWarning("Clusterizer already initialized. Unable to change the parameters.");
  else
    fShiftPhi = s;
}

//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::SetShiftEta (Int_t s) 
{
  // Set fShiftEta; if clusterizer already initialized gives a warning and does nothing
  
  if (fClustersArray)
    AliWarning("Clusterizer already initialized. Unable to change the parameters.");
  else
    fShiftEta = s;
}

//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::SetTRUshift(Bool_t b) 
{
  // Set fTRUshift; if clusterizer already initialized gives a warning and does nothing
  
  if (fClustersArray)
    AliWarning("Clusterizer already initialized. Unable to change the parameters.");
  else
    fTRUshift = b;
}


//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::Digits2Clusters(Option_t * option)
{
  // Steering method to perform clusterization for the current event 
	
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
    printf("Exec took %f seconds for Clusterizing", 
           gBenchmark->GetCpuTime("EMCALClusterizer"));
  }    
}

//__________________________________________________________________________________________
void AliEMCALClusterizerFixedWindow::MakeClusters()
{
  // Make clusters
	
  if (fGeom == 0) 
    AliFatal("Did not get geometry from EMCALLoader");
	
  fNumberOfECAClusters = 0;
  fRecPoints->Delete();
	
  Int_t nSupMod=0, nModule=0, nIphi=0, nIeta=0, iphi=0, ieta=0;
	
  // Defining geometry and clusterization parameter
  Int_t nEtaDigitsSupMod = fGeom->GetNEta() * fGeom->GetNETAdiv(); // always 48?;
  Int_t nPhiDigitsSupMod = fGeom->GetNPhi() * fGeom->GetNPHIdiv(); // always 24?;
  
  Int_t nTRUPhi = 1;
  Int_t nTRUEta = 1;
  
  Int_t nEtaDigits = nEtaDigitsSupMod * fGeom->GetNumberOfSuperModules() / fGeom->GetNPhiSuperModule();
  Int_t nPhiDigits = nPhiDigitsSupMod * fGeom->GetNPhiSuperModule();    
  
  if (fTRUshift){
    nTRUPhi = fGeom->GetNPhiSuperModule() * 3;
    nTRUEta = fGeom->GetNumberOfSuperModules() / fGeom->GetNPhiSuperModule();
    nEtaDigits /= nTRUEta;
    nPhiDigits /= nTRUPhi;
  }

  // Check if clusterizer parameter are compatible with calorimeter geometry
  if (nEtaDigits < fNeta){
    AliFatal(Form("Error: fNeta = %d is greater than nEtaDigits = %d.",fNeta,nEtaDigits));
    return;
  }
  if (nPhiDigits < fNphi){
    AliFatal(Form("Error: fNphi = %d is greater than nPhiDigits = %d.",fNphi,nPhiDigits));
    return;
  }
  if (nEtaDigits % fShiftEta != 0){
    AliFatal(Form("Error: fShiftEta = %d is such that clusters cannot slide the whole calorimeter (nEtaDigits = %d).",fShiftEta,nEtaDigits));
    return;
  }
  if (nPhiDigits % fShiftPhi != 0){
    AliFatal(Form("Error: fShiftPhi = %d is such that clusters cannot slide the whole calorimeter (nPhiDigits = %d).",fShiftPhi,nPhiDigits));
    return;
  }
  if (fNeta % fShiftEta != 0){
    AliFatal(Form("Error: fShiftEta = %d is not divisor of fNeta = %d.",fShiftEta,fNeta));
    return;
  }
  if (fNphi % fShiftPhi != 0){
    AliFatal(Form("Error: fShiftPhi = %d is not divisor of fNphi = %d).",fShiftPhi,fNphi));
    return;
  }
  
  Int_t maxiShiftPhi = fNphi / fShiftPhi;
  Int_t maxiShiftEta = fNeta / fShiftEta;
	
  Int_t nDigitsCluster = fNphi * fNeta;
  
  Int_t nClusEtaNoShift = nEtaDigits / fNeta;
  Int_t nClusPhiNoShift = nPhiDigits / fNphi;
  
  Int_t nClusters =  nClusEtaNoShift * nClusPhiNoShift * nTRUEta * nTRUPhi;
  
  Int_t nTotalClus = nClusters * maxiShiftEta * maxiShiftPhi;
  
  if (!fClustersArray) {
    fClustersArray = new AliEMCALDigit**[nTotalClus];
    for (Int_t i = 0; i < nTotalClus; i++)
    {
      fClustersArray[i] = NULL;
    }
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
    if (dEnergyCalibrated < fMinECut) {
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
   
  for (Int_t ishiftPhi = 0; ishiftPhi < maxiShiftPhi; ishiftPhi++){
    Int_t nClusPhi = (nPhiDigits - fShiftPhi * ishiftPhi) / fNphi;
    
    for (Int_t ishiftEta = 0; ishiftEta < maxiShiftEta; ishiftEta++) {
      
      Int_t nClusEta = (nEtaDigits - fShiftEta * ishiftEta) / fNeta; 
      
      Int_t iTotalClus = nClusters * (ishiftPhi * maxiShiftEta + ishiftEta);
      
      TIter nextdigitC(digitsC);
      while ((digit = dynamic_cast<AliEMCALDigit*>(nextdigitC()))) { // scan over the list of digitsC
        
        fGeom->GetCellIndex (digit->GetId(), nSupMod, nModule, nIphi, nIeta);
        fGeom->GetCellPhiEtaIndexInSModule (nSupMod, nModule, nIphi, nIeta, iphi, ieta);
        
        Int_t iphi_eff = iphi - fShiftPhi * ishiftPhi + nPhiDigitsSupMod * (nSupMod / 2); // N supermodules along phi
        
        Int_t iTRUphi = iphi_eff / nPhiDigits;
        
        iphi_eff -= iTRUphi * nPhiDigits;
        
        Int_t iClusPhi = iphi_eff / fNphi; 
        
        if (iphi_eff < 0 || iClusPhi >= nClusPhi) 
          continue;
        
        Int_t ieta_eff = ieta - fShiftEta * ishiftEta + nEtaDigitsSupMod * (nSupMod % 2); // 2 supermodules along eta
        
        Int_t iTRUeta = ieta_eff / nEtaDigits;
        
        ieta_eff -= iTRUeta * nEtaDigits;
        
        Int_t iClusEta = ieta_eff / fNeta; 
        
        if (ieta_eff < 0 || iClusEta >= nClusEta) 
          continue;
        
        iphi_eff += iTRUphi * nPhiDigits;
        iClusPhi = iphi_eff / fNphi; 
        
        ieta_eff += iTRUeta * nEtaDigits;
        iClusEta = ieta_eff / fNeta; 
        
        Int_t iCluster = iClusPhi + iClusEta * nClusPhiNoShift * nTRUPhi; 
        Int_t iDigit = iphi_eff % fNphi + (ieta_eff % fNeta) * fNphi;

        if (iCluster >= nClusters){
          AliWarning(Form("iCluster out of range! iCluster = %d, nClusters = %d", iCluster, nClusters));
          return;
        }
        
        iCluster += iTotalClus;
        
        if (fClustersArray[iCluster] == NULL){
          fNumberOfECAClusters++;
          fClustersArray[iCluster] = new AliEMCALDigit*[nDigitsCluster];
          for (Int_t i = 0; i < nDigitsCluster; i++){
            fClustersArray[iCluster][i] = NULL;
          }
        }
        
        if (fClustersArray[iCluster][iDigit] != NULL){
          AliWarning("Digit already added!");
          return;
        }
        
        fClustersArray[iCluster][iDigit] = digit;
        
      } // loop on digit
      
    } // loop on eta shift
    
  } // loop on phi shift
	
  Int_t iRecPoint = 0;
  for (Int_t iCluster = 0; iCluster < nTotalClus; iCluster++){
    
    if (fClustersArray[iCluster] == NULL) continue;
		
    (*fRecPoints)[iRecPoint] = new AliEMCALRecPoint("");
    AliEMCALRecPoint *recPoint = dynamic_cast<AliEMCALRecPoint*> (fRecPoints->At(iRecPoint));
		
    if (recPoint) {
      iRecPoint++;       
      recPoint->SetClusterType(AliVCluster::kEMCALClusterv1);
      recPoint->SetUniqueID(iCluster);

      for (Int_t iDigit = 0; iDigit < nDigitsCluster; iDigit++){
        if (fClustersArray[iCluster][iDigit] == NULL) continue;
        digit = fClustersArray[iCluster][iDigit];
        recPoint->AddDigit(*digit, digit->GetCalibAmp(), kFALSE); //Time or TimeR?
        fClustersArray[iCluster][iDigit] = NULL;
      }
    }
    
    delete[] fClustersArray[iCluster];
    fClustersArray[iCluster] = NULL;
  }
  
  AliDebug(1, Form("MakeClusters: Number of digits %d  -> (e %f)\n", fDigitsArr->GetEntries(),fMinECut));
	
  AliDebug(1, Form("total no of clusters %d from %d digits", fNumberOfECAClusters, fDigitsArr->GetEntriesFast())); 
}
