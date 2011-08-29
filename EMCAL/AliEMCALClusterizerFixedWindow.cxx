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

// This class derives from AliEMCALClustrerizer but also keeps the API of AliEMCALClusterizerv1
// Algorithm:
// 1. peek the most energetic cell
// 2. assign it as a center of the cluster and add cells surrounding it: 3x3, 5x5...
// 3. remove the cells contributing to the cluster
// 4. start from 1 for the remaining clusters
// 5. cluster splitting (not implemented yet) - use the shape analysis to resolve the energy sharing
// - for high energy clusters check the surrounding of the 3x3 clusters for extra energy 
// (merge 3x3 clusters and resolve the internal energy sharing - case for 2 clusters merged)
// Use Case:
//  root [0] AliEMCALClusterizerFixedWindow * cl = new AliEMCALClusterizerFixedWindow("galice.root")  
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//               //reads gAlice from header file "..."                      
//  root [1] cl->ExecuteTask()  
//               //finds RecPoints in all events stored in galice.root
//  root [2] cl->SetDigitsBranch("digits2") 
//               //sets another title for Digitis (input) branch
//  root [3] cl->SetRecPointsBranch("recp2")  
//               //sets another title four output branches
//  root [4] cl->SetTowerLocalMaxCut(0.03)  
//               //set clusterization parameters
//  root [5] cl->ExecuteTask("deb all time")  
//               //once more finds RecPoints options are 
//               // deb - print number of found rec points
//               // deb all - print number of found RecPoints and some their characteristics 
//               // time - print benchmarking results

#include "AliEMCALClusterizerFixedWindow.h"

// --- ROOT system ---
#include <TMath.h> 
#include <TMinuit.h>
#include <TTree.h> 
#include <TBenchmark.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TH1I.h>

// --- Standard library ---
#include <cassert>
//#include <iostream>
//#include <fstream>

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALCalibData.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALUnfolding.h"
#include "AliEMCALFixedWindowClusterInfo.h"

ClassImp(AliEMCALClusterizerFixedWindow)

//____________________________________________________________________________
AliEMCALClusterizerFixedWindow::AliEMCALClusterizerFixedWindow()
: AliEMCALClusterizer(), 
nPhi(4), 
nEta(4), 
shiftPhi(2), 
shiftEta(2),
fTRUshift(0), 
clusters_array(0),
fClustersInfo(new AliEMCALFixedWindowClusterInfo("clustersInfo"))
{
	// ctor with the indication of the file where header Tree and digits Tree are stored

}

//____________________________________________________________________________
AliEMCALClusterizerFixedWindow::AliEMCALClusterizerFixedWindow(AliEMCALGeometry* geometry)
: AliEMCALClusterizer(geometry), 
nPhi(4), 
nEta(4), 
shiftPhi(2), 
shiftEta(2), 
fTRUshift(0), 
clusters_array(0),
fClustersInfo(new AliEMCALFixedWindowClusterInfo("clustersInfo"))
{
	// ctor with the indication of the file where header Tree and digits Tree are stored
	// use this contructor to avoid usage of Init() which uses runloader
	// change needed by HLT - MP
}

//____________________________________________________________________________
AliEMCALClusterizerFixedWindow::AliEMCALClusterizerFixedWindow(AliEMCALGeometry* geometry, AliEMCALCalibData * calib, AliCaloCalibPedestal * caloped)
: AliEMCALClusterizer(geometry, calib, caloped),
nPhi(4), 
nEta(4), 
shiftPhi(2), 
shiftEta(2), 
fTRUshift(0), 
clusters_array(0),
fClustersInfo(new AliEMCALFixedWindowClusterInfo("clustersInfo"))
{
	// ctor, geometry and calibration are initialized elsewhere.
}

//____________________________________________________________________________
AliEMCALClusterizerFixedWindow::~AliEMCALClusterizerFixedWindow()
{
	// dtor
  delete fClustersInfo;
}

//____________________________________________________________________________
void AliEMCALClusterizerFixedWindow::Digits2Clusters(Option_t * option)
{
	// Steering method to perform clusterization for the current event 
	// in AliEMCALLoader
	
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
	
	//fRecPoints->Sort();
	
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

//____________________________________________________________________________
void AliEMCALClusterizerFixedWindow::MakeClusters()
{
	// Make clusters
	
	if (fGeom == 0) 
		AliFatal("Did not get geometry from EMCALLoader");
	
	fNumberOfECAClusters = 0;
	fRecPoints->Delete();
  
  if (fClustersInfo->GetLastElementId() > 0)
    fClustersInfo->Clear();
	
	Int_t nSupMod=0, nModule=0, nIphi=0, nIeta=0, iphi=0, ieta=0;
	
	// Defining geometry and clusterization parameter
	Int_t nEtaDigitsSupMod = fGeom->GetNEta() * fGeom->GetNETAdiv(); // always 48?;
	Int_t nPhiDigitsSupMod = fGeom->GetNPhi() * fGeom->GetNPHIdiv(); // always 24?;
  
  Int_t nTRUPhi = 1;
  Int_t nTRUEta = 1;
  
  Int_t nEtaDigits = nEtaDigitsSupMod * fGeom->GetNumberOfSuperModules() / fGeom->GetNPhiSuperModule();
  Int_t nPhiDigits = nPhiDigitsSupMod * fGeom->GetNPhiSuperModule();    
  
  if (fTRUshift)
  {
    nTRUPhi = fGeom->GetNPhiSuperModule() * 3;
    nTRUEta = fGeom->GetNumberOfSuperModules() / fGeom->GetNPhiSuperModule();
    nEtaDigits /= nTRUEta;
    nPhiDigits /= nTRUPhi;
  }

  // Check if clusterizer parameter are compatible with calorimeter geometry
  if (nEtaDigits < nEta)
	{
		AliFatal(Form("Error: nEta = %d is greater than nEtaDigits = %d.",nEta,nEtaDigits));
		return;
	}
	if (nPhiDigits < nPhi)
	{
		AliFatal(Form("Error: nPhi = %d is greater than nPhiDigits = %d.",nPhi,nPhiDigits));
		return;
	}
	if (nEtaDigits % shiftEta != 0)
	{
		AliFatal(Form("Error: shiftEta = %d is such that clusters cannot slide the whole calorimeter (nEtaDigits = %d).",shiftEta,nEtaDigits));
		return;
	}
	if (nPhiDigits % shiftPhi != 0)
	{
		AliFatal(Form("Error: shiftPhi = %d is such that clusters cannot slide the whole calorimeter (nPhiDigits = %d).",shiftPhi,nPhiDigits));
		return;
	}
	if (nEta % shiftEta != 0)
	{
		AliFatal(Form("Error: shiftEta = %d is not divisor of nEta = %d.",shiftEta,nEta));
		return;
	}
	if (nPhi % shiftPhi != 0)
	{
		AliFatal(Form("Error: shiftPhi = %d is not divisor of nPhi = %d).",shiftPhi,nPhi));
		return;
	}
  
  Int_t maxiShiftPhi = nPhi / shiftPhi;
  Int_t maxiShiftEta = nEta / shiftEta;
	
	Int_t nDigitsCluster = nPhi * nEta;
  
  Int_t nClusEtaNoShift = nEtaDigits / nEta;
  Int_t nClusPhiNoShift = nPhiDigits / nPhi;
  
  Int_t nClusters =  nClusEtaNoShift * nClusPhiNoShift * nTRUEta * nTRUPhi;
  
  Int_t nTotalClus = nClusters * maxiShiftEta * maxiShiftPhi;
  
  if (!clusters_array)
  {
    clusters_array = new AliEMCALDigit**[nTotalClus];
    for (Int_t i = 0; i < nTotalClus; i++)
    {
      clusters_array[i] = NULL;
    }
  }
  
  AliEMCALDigit *digit = 0;
  
  for (Int_t ishiftPhi = 0; ishiftPhi < maxiShiftPhi; ishiftPhi++)
  {
    Int_t nClusPhi = (nPhiDigits - shiftPhi * ishiftPhi) / nPhi;
    
    for (Int_t ishiftEta = 0; ishiftEta < maxiShiftEta; ishiftEta++)
    {
      
      Int_t nClusEta = (nEtaDigits - shiftEta * ishiftEta) / nEta; 
      
      Int_t iTotalClus = nClusters * (ishiftPhi * maxiShiftEta + ishiftEta);
      
      TIter nextdigit(fDigitsArr);
      
      nextdigit.Reset();
      
      while (digit = static_cast<AliEMCALDigit*>(nextdigit()))
      {
        fGeom->GetCellIndex (digit->GetId(), nSupMod, nModule, nIphi, nIeta);
        fGeom->GetCellPhiEtaIndexInSModule (nSupMod, nModule, nIphi, nIeta, iphi, ieta);
        
        Int_t iphi_eff = iphi - shiftPhi * ishiftPhi + nPhiDigitsSupMod * (nSupMod / 2); // N supermodules along phi
        
        Int_t iTRUphi = iphi_eff / nPhiDigits;
        
        iphi_eff -= iTRUphi * nPhiDigits;
        
        Int_t iClusPhi = iphi_eff / nPhi; 
        
        if (iphi_eff < 0 || iClusPhi >= nClusPhi) 
          continue;
        
        Int_t ieta_eff = ieta - shiftEta * ishiftEta + nEtaDigitsSupMod * (nSupMod % 2); // 2 supermodules along eta
        
        Int_t iTRUeta = ieta_eff / nEtaDigits;
        
        ieta_eff -= iTRUeta * nEtaDigits;
        
        Int_t iClusEta = ieta_eff / nEta; 
        
        if (ieta_eff < 0 || iClusEta >= nClusEta) 
          continue;
        
        iphi_eff += iTRUphi * nPhiDigits;
        iClusPhi = iphi_eff / nPhi; 
        
        ieta_eff += iTRUeta * nEtaDigits;
        iClusEta = ieta_eff / nEta; 
        
        Int_t iCluster = iClusPhi + iClusEta * nClusPhiNoShift * nTRUPhi; 
        Int_t iDigit = iphi_eff % nPhi + (ieta_eff % nEta) * nPhi;

        
        if (iCluster >= nClusters)
        {
          AliFatal(Form("ERROR: iCluster out of range! iCluster = %d, nClusters = %d", iCluster, nClusters));
          return;
        }
        
        iCluster += iTotalClus;
        
        if (clusters_array[iCluster] == NULL)
        {
          fNumberOfECAClusters++;
          clusters_array[iCluster] = new AliEMCALDigit*[nDigitsCluster];
          for (Int_t i = 0; i < nDigitsCluster; i++)
          {
            clusters_array[iCluster][i] = NULL;
          }
          
          fClustersInfo->Add(iCluster, -1, iClusEta, iClusPhi);
        }
        
        if (clusters_array[iCluster][iDigit] != NULL)
        {
          AliFatal("ERROR: digit already added!");
          return;
        }
        
        clusters_array[iCluster][iDigit] = digit;
        
      } // loop on digit
      
    } // loop on eta shift
    
	} // loop on phi shift
	
  Int_t iRecPoint = 0;
	for (Int_t iCluster = 0; iCluster < nTotalClus; iCluster++)
	{
    
		if (clusters_array[iCluster] == NULL) continue;
		
		(*fRecPoints)[iRecPoint] = new AliEMCALRecPoint("");
		AliEMCALRecPoint *recPoint = dynamic_cast<AliEMCALRecPoint*> (fRecPoints->At(iRecPoint));
		
		if (recPoint) 
		{
      if (fClustersInfo->ContainsIndex(iRecPoint))
        AliFatal(Form("ERROR: index present already, %d", iRecPoint));
      
      fClustersInfo->SetIndexFromId(iCluster, iRecPoint);
      
			iRecPoint++;       
			recPoint->SetClusterType(AliVCluster::kEMCALClusterv1);
			// note: this way the sharing info is lost!
			for (Int_t iDigit = 0; iDigit < nDigitsCluster; iDigit++)
			{
				if (clusters_array[iCluster][iDigit] == NULL) continue;
				digit = clusters_array[iCluster][iDigit];
        Float_t dEnergyCalibrated = digit->GetAmplitude();
        Float_t time              = digit->GetTime();
        Calibrate(dEnergyCalibrated,time,digit->GetId());
				digit->SetCalibAmp(dEnergyCalibrated);
				recPoint->AddDigit(*digit, dEnergyCalibrated, kFALSE); //Time or TimeR?
        clusters_array[iCluster][iDigit] = NULL;
			}
		}
    
    delete[] clusters_array[iCluster];
    clusters_array[iCluster] = NULL;
	}
  
	AliDebug(1,Form("MakeClusters: Number of digits %d  -> (e %f)\n",
					fDigitsArr->GetEntries(),fMinECut));
	
	AliDebug(1,Form("total no of clusters %d from %d digits",fNumberOfECAClusters,fDigitsArr->GetEntriesFast())); 
}
