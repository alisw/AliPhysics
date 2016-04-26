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

//-- Author: Yves Schutz (SUBATECH)  & Dmitri Peressounko (SUBATECH & Kurchatov Institute)
//--         Gustavo Conesa (LPSC-Grenoble), move common clusterizer functionalities to mother class
//////////////////////////////////////////////////////////////////////////////
//  Clusterization class. Performs clusterization (collects neighbouring active cells) and 
//  unfolds the clusters having several local maxima.  
//  Results are stored in TreeR#, branches EMCALTowerRP (EMC recPoints),
//  EMCALPreShoRP (CPV RecPoints) and AliEMCALClusterizer (Clusterizer with all 
//  parameters including input digits branch title, thresholds etc.)
//

// --- ROOT system ---

#include <TFile.h> 
#include <TMath.h> 
#include <TMinuit.h>
#include <TTree.h> 
#include <TBenchmark.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TList.h>
#include <TClonesArray.h>

// --- Standard library ---
#include <cassert>

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliEMCALClusterizerv2.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALCalibData.h"
#include "AliEMCALCalibTime.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALUnfolding.h"

ClassImp(AliEMCALClusterizerv2)

//____________________________________________________________________________
AliEMCALClusterizerv2::AliEMCALClusterizerv2() 
  : AliEMCALClusterizerv1(), fDoEnGradCut(1)
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
}

//____________________________________________________________________________
AliEMCALClusterizerv2::AliEMCALClusterizerv2(AliEMCALGeometry* geometry)
  : AliEMCALClusterizerv1(geometry), fDoEnGradCut(1)
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
  // use this contructor to avoid usage of Init() which uses runloader
  // change needed by HLT - MP
}

//____________________________________________________________________________
AliEMCALClusterizerv2::AliEMCALClusterizerv2(AliEMCALGeometry* geometry, 
                                             AliEMCALCalibData* calib, 
                                             AliEMCALCalibTime* calibt, 
                                             AliCaloCalibPedestal* caloped)
  : AliEMCALClusterizerv1(geometry, calib, calibt, caloped), fDoEnGradCut(1)
{
  // ctor, geometry and calibration are initialized elsewhere.
}

//____________________________________________________________________________
AliEMCALClusterizerv2::~AliEMCALClusterizerv2()
{
  // dtor
}

//____________________________________________________________________________
Int_t AliEMCALClusterizerv2::AreNeighbours(AliEMCALDigit* d1, AliEMCALDigit* d2, Bool_t& shared) const
{ 
  // Gives the neighbourness of two digits = 0 are not neighbour; continue searching 
  //                                       = 1 are neighbour
  //                                       = 2 is in different SM; continue searching 
  // In case it is in different SM, but same phi rack, check if neigbours at eta=0
  // neighbours are defined as digits having at least a common side 
  // The order of d1 and d2 is important: first (d1) should be a digit already in a cluster 
  //                                      which is compared to a digit (d2)  not yet in a cluster  
  
  if (fDoEnGradCut) {
    if (d2->GetCalibAmp()>d1->GetCalibAmp()+fECALocMaxCut){
      return 3; // energy of neighboring cell should be smaller in order to become a neighbor
    }  
  }
  return AliEMCALClusterizerv1::AreNeighbours(d1,d2,shared);
}

//____________________________________________________________________________
void AliEMCALClusterizerv2::MakeClusters()
{
  // Make list of clusters. Start from highest energy cell.
  
  if (fGeom==0) 
    AliFatal("Did not get geometry from EMCALLoader");
  
  fRecPoints->Delete();

  // set up TObjArray with pointers to digits to work on, calibrate digits 
  TObjArray digitsC;
  Double_t ehs = 0.0;
  AliEMCALDigit *digit = 0;
  TIter nextdigit(fDigitsArr);
  while ( (digit = static_cast<AliEMCALDigit*>(nextdigit())) ) {
    Float_t dEnergyCalibrated = digit->GetAmplitude();
    Float_t time              = digit->GetTime();
    Calibrate(dEnergyCalibrated, time ,digit->GetId());
    digit->SetCalibAmp(dEnergyCalibrated);
    digit->SetTime(time);
    if (dEnergyCalibrated < fMinECut || time > fTimeMax || time < fTimeMin) 
      continue;
    if (!fGeom->CheckAbsCellId(digit->GetId()))
      continue;
    ehs += dEnergyCalibrated;
    digitsC.AddLast(digit);
  }

  AliDebug(1,Form("MakeClusters: Number of digits %d  -> ehs %f (minE %f)\n",
                  fDigitsArr->GetEntries(),ehs,fMinECut));

  TIter nextdigitC(&digitsC);
  while (1) {
    Int_t   iMaxEnergyDigit = -1;
    Float_t dMaxEnergyDigit = -1;
    AliEMCALDigit *pMaxEnergyDigit = 0;
    nextdigitC.Reset();
    while ( (digit = static_cast<AliEMCALDigit *>(nextdigitC())) ) { 
      Float_t dEnergyCalibrated = digit->GetCalibAmp();
      if (dEnergyCalibrated>fECAClusteringThreshold && dEnergyCalibrated>dMaxEnergyDigit) {
        dMaxEnergyDigit = dEnergyCalibrated;
        iMaxEnergyDigit = digit->GetId();
        pMaxEnergyDigit = digit;
      }
    }
    if (iMaxEnergyDigit<0 || digitsC.GetEntries() <= 0) {
      break;
    }

    if (fNumberOfECAClusters>=fRecPoints->GetSize()) 
      fRecPoints->Expand(2*fNumberOfECAClusters+1);

    AliEMCALRecPoint *recPoint = new  AliEMCALRecPoint(""); 
    recPoint->SetClusterType(AliVCluster::kEMCALClusterv1);
    recPoint->AddDigit(*pMaxEnergyDigit, pMaxEnergyDigit->GetCalibAmp(), kFALSE);
    fRecPoints->AddAt(recPoint, fNumberOfECAClusters++);
    digitsC.Remove(pMaxEnergyDigit); 
    TObjArray clusterDigits;
    clusterDigits.AddLast(pMaxEnergyDigit);	
    TIter nextClusterDigit(&clusterDigits);
    Float_t time = pMaxEnergyDigit->GetTime(); 

    AliDebug(1,Form("MakeClusters: Max digit found id = %d, ene = %f , clus.th. = %f \n", 
                    iMaxEnergyDigit, dMaxEnergyDigit, fECAClusteringThreshold));

    while ( (digit = static_cast<AliEMCALDigit*>(nextClusterDigit())) ) { // scan over digits in cluster 
      TIter nextdigitN(&digitsC); 
      AliEMCALDigit *digitN = 0; // digi neighbor
      while ( (digitN = static_cast<AliEMCALDigit*>(nextdigitN())) ) { // scan over all digits to look for neighbours
        //Do not add digits with too different time 
        if (TMath::Abs(time - digitN->GetTime()) > fTimeCut ) 
          continue;
        Bool_t shared = kFALSE; //cluster shared by 2 SuperModules?
        if (AreNeighbours(digit, digitN, shared)==1) {
          recPoint->AddDigit(*digitN, digitN->GetCalibAmp(), shared);
          clusterDigits.AddLast(digitN); 
          digitsC.Remove(digitN); 
        } 
      }
    }
    AliDebug(2,Form("MakeClusters: %d digitd, energy %f \n", clusterDigits.GetEntries(), recPoint->GetEnergy())); 
  }
  AliDebug(1,Form("total no of clusters %d from %d digits",fNumberOfECAClusters,fDigitsArr->GetEntriesFast())); 
}

