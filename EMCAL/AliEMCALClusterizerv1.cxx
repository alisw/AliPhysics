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
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALCalibData.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALUnfolding.h"

ClassImp(AliEMCALClusterizerv1)

//____________________________________________________________________________
AliEMCALClusterizerv1::AliEMCALClusterizerv1(): AliEMCALClusterizer()
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
  
  Init();
}

//____________________________________________________________________________
AliEMCALClusterizerv1::AliEMCALClusterizerv1(AliEMCALGeometry* geometry)
  : AliEMCALClusterizer(geometry)
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
  // use this contructor to avoid usage of Init() which uses runloader
  // change needed by HLT - MP
}

//____________________________________________________________________________
AliEMCALClusterizerv1::AliEMCALClusterizerv1(AliEMCALGeometry* geometry, AliEMCALCalibData * calib, AliCaloCalibPedestal * caloped)
: AliEMCALClusterizer(geometry, calib, caloped)
{
  // ctor, geometry and calibration are initialized elsewhere.
}

//____________________________________________________________________________
  AliEMCALClusterizerv1::~AliEMCALClusterizerv1()
{
  // dtor
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::Digits2Clusters(Option_t * option)
{
  // Steering method to perform clusterization for the current event 
  // in AliEMCALLoader
  
  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALClusterizer"); 
  
  if(strstr(option,"print"))
    Print(""); 
  
  //Get calibration parameters from file or digitizer default values.
  GetCalibrationParameters();
  
  //Get dead channel map from file or digitizer default values.
  GetCaloCalibPedestal();
	
  fNumberOfECAClusters = 0;
  
  MakeClusters();  //only the real clusters
  
  if(fToUnfold){
    fClusterUnfolding->SetInput(fNumberOfECAClusters,fRecPoints,fDigitsArr);
    fClusterUnfolding->MakeUnfolding();
  }
    
  //Evaluate position, dispersion and other RecPoint properties for EC section 
  Int_t index;
  for(index = 0; index < fRecPoints->GetEntries(); index++) {
    AliEMCALRecPoint * rp = dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(index));
    if(rp){
      rp->EvalAll(fECAW0,fDigitsArr,fIsInputCalibrated);
      //For each rec.point set the distance to the nearest bad crystal
      if (fCaloPed)
        rp->EvalDistanceToBadChannels(fCaloPed);
    }
    else AliFatal("Null rec point in list!");
  }
  
  fRecPoints->Sort();
  
  for(index = 0; index < fRecPoints->GetEntries(); index++) {
    AliEMCALRecPoint * rp = dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(index));
    if(rp){
      rp->SetIndexInList(index);
      rp->Print();
    }
    else AliFatal("Null rec point in list!");
  }
  
  if (fTreeR)
    fTreeR->Fill();
  
  if(strstr(option,"deb") || strstr(option,"all"))  
    PrintRecPoints(option);
  
  AliDebug(1,Form("EMCAL Clusterizer found %d Rec Points",fRecPoints->GetEntriesFast()));
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALClusterizer");
    printf("Exec took %f seconds for Clusterizing", 
           gBenchmark->GetCpuTime("EMCALClusterizer"));
  }    
}

//____________________________________________________________________________
Int_t AliEMCALClusterizerv1::AreNeighbours(AliEMCALDigit * d1, AliEMCALDigit * d2, Bool_t & shared) const
{
  // Gives the neighbourness of two digits = 0 are not neighbour; continue searching 
  //                                       = 1 are neighbour
  //                                       = 2 is in different SM; continue searching 
  // In case it is in different SM, but same phi rack, check if neigbours at eta=0
  // neighbours are defined as digits having at least a common side 
  // The order of d1 and d2 is important: first (d1) should be a digit already in a cluster 
  //                                      which is compared to a digit (d2)  not yet in a cluster  
  
  Int_t nSupMod1=0, nModule1=0, nIphi1=0, nIeta1=0, iphi1=0, ieta1=0;
  Int_t nSupMod2=0, nModule2=0, nIphi2=0, nIeta2=0, iphi2=0, ieta2=0;
  
  shared = kFALSE;
  
  fGeom->GetCellIndex(d1->GetId(), nSupMod1,nModule1,nIphi1,nIeta1);
  fGeom->GetCellIndex(d2->GetId(), nSupMod2,nModule2,nIphi2,nIeta2);
  fGeom->GetCellPhiEtaIndexInSModule(nSupMod1,nModule1,nIphi1,nIeta1, iphi1,ieta1);
  fGeom->GetCellPhiEtaIndexInSModule(nSupMod2,nModule2,nIphi2,nIeta2, iphi2,ieta2);
  
  //If different SM, check if they are in the same phi, then consider cells close to eta=0 as neighbours; May 2010
  if (nSupMod1 != nSupMod2 ) {
    //Check if the 2 SM are in the same PHI position (0,1), (2,3), ...
    Float_t smPhi1 = fGeom->GetEMCGeometry()->GetPhiCenterOfSM(nSupMod1);
    Float_t smPhi2 = fGeom->GetEMCGeometry()->GetPhiCenterOfSM(nSupMod2);
    
    if(!TMath::AreEqualAbs(smPhi1, smPhi2, 1e-3)) return 2; //Not phi rack equal, not neighbours
    
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2
    // C Side impair SM, nSupMod%2=1; A side pair SM nSupMod%2=0
    if(nSupMod1%2) ieta1+=AliEMCALGeoParams::fgkEMCALCols;
    else           ieta2+=AliEMCALGeoParams::fgkEMCALCols;
    
    shared = kTRUE; // maybe a shared cluster, we know this later, set it for the moment.
  } //Different SM, same phi
  
  Int_t rowdiff = TMath::Abs(iphi1 - iphi2);  
  Int_t coldiff = TMath::Abs(ieta1 - ieta2);  
  
  // neighbours with at least common side; May 11, 2007
  if ((coldiff==0 && TMath::Abs(rowdiff)==1) || (rowdiff==0 && TMath::Abs(coldiff)==1)) {  
    //Diagonal?
    //if ((coldiff==0 && TMath::Abs(rowdiff==1)) || (rowdiff==0 && TMath::Abs(coldiff==1)) || (TMath::Abs(rowdiff)==1 && TMath::Abs(coldiff==1))) rv = 1;
    
    if (gDebug == 2) 
      printf("AliEMCALClusterizerv1::AreNeighbours(): id1=%d, (row %d, col %d) ; id2=%d, (row %d, col %d), shared %d \n",
	     d1->GetId(), iphi1,ieta1, d2->GetId(), iphi2,ieta2, shared);   
    return 1;
  } //Neighbours
  else {
    shared = kFALSE;
    return 2; 
  } //Not neighbours
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::MakeClusters()
{
  // Steering method to construct the clusters stored in a list of Reconstructed Points
  // A cluster is defined as a list of neighbour digits
  // Mar 03, 2007 by PAI
  
  if (fGeom==0) AliFatal("Did not get geometry from EMCALLoader");
  
  fRecPoints->Delete();
  
  // Set up TObjArray with pointers to digits to work on 
  TObjArray *digitsC = new TObjArray();
  TIter nextdigit(fDigitsArr);
  AliEMCALDigit *digit;
  while ( (digit = dynamic_cast<AliEMCALDigit*>(nextdigit())) ) {
    digitsC->AddLast(digit);
  }
  
  double e = 0.0, ehs = 0.0;
  TIter nextdigitC(digitsC);
  while ( (digit = dynamic_cast<AliEMCALDigit *>(nextdigitC())) ) { // clean up digits
    e = Calibrate(digit->GetAmplitude(), digit->GetTime(),digit->GetId());//Time or TimeR?
    if ( e < fMinECut) //|| digit->GetTimeR() > fTimeCut ) // time window of cell checked in calibrate
      digitsC->Remove(digit);
    else    
      ehs += e;
  } 
  AliDebug(1,Form("MakeClusters: Number of digits %d  -> (e %f), ehs %f\n",
                  fDigitsArr->GetEntries(),fMinECut,ehs));
  
  nextdigitC.Reset();
  
  while ( (digit = dynamic_cast<AliEMCALDigit *>(nextdigitC())) ) { // scan over the list of digitsC
    TArrayI clusterECAdigitslist(fDigitsArr->GetEntries());
    
    if(fGeom->CheckAbsCellId(digit->GetId()) && (Calibrate(digit->GetAmplitude(), digit->GetTime(),digit->GetId()) > fECAClusteringThreshold  ) ){
      // start a new Tower RecPoint
      if(fNumberOfECAClusters >= fRecPoints->GetSize()) fRecPoints->Expand(2*fNumberOfECAClusters+1);
      
      AliEMCALRecPoint *recPoint = new  AliEMCALRecPoint(""); 
      fRecPoints->AddAt(recPoint, fNumberOfECAClusters);
      recPoint = dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(fNumberOfECAClusters)); 
      if (recPoint) {
        fNumberOfECAClusters++; 
        
        recPoint->SetClusterType(AliVCluster::kEMCALClusterv1);
        
        recPoint->AddDigit(*digit, Calibrate(digit->GetAmplitude(), digit->GetTime(),digit->GetId()),kFALSE); //Time or TimeR?
        TObjArray clusterDigits;
        clusterDigits.AddLast(digit);	
        digitsC->Remove(digit); 
        
        AliDebug(1,Form("MakeClusters: OK id = %d, ene = %f , cell.th. = %f \n", digit->GetId(),
                        Calibrate(digit->GetAmplitude(),digit->GetTime(),digit->GetId()), fECAClusteringThreshold));  //Time or TimeR?
        Float_t time = digit->GetTime();//Time or TimeR?
        // Grow cluster by finding neighbours
        TIter nextClusterDigit(&clusterDigits);
        while ( (digit = dynamic_cast<AliEMCALDigit*>(nextClusterDigit())) ) { // scan over digits in cluster 
          TIter nextdigitN(digitsC); 
          AliEMCALDigit *digitN = 0; // digi neighbor
          while ( (digitN = (AliEMCALDigit *)nextdigitN()) ) { // scan over all digits to look for neighbours
            
            //Do not add digits with too different time 
            Bool_t shared = kFALSE;//cluster shared by 2 SuperModules?
            if(TMath::Abs(time - digitN->GetTime()) > fTimeCut ) continue; //Time or TimeR?
            if (AreNeighbours(digit, digitN, shared)==1) {      // call (digit,digitN) in THAT order !!!!! 
              recPoint->AddDigit(*digitN, Calibrate(digitN->GetAmplitude(), digitN->GetTime(), digitN->GetId()),shared); //Time or TimeR?
              clusterDigits.AddLast(digitN); 
              digitsC->Remove(digitN); 
            } // if(ineb==1)
          } // scan over digits
        } // scan over digits already in cluster
        
        AliDebug(2,Form("MakeClusters: %d digitd, energy %f \n", clusterDigits.GetEntries(), recPoint->GetEnergy())); 
      }//recpoint
      else AliFatal("Null recpoint in array!");
    } // If seed found
  } // while digit 
  
  delete digitsC;
  
  AliDebug(1,Form("total no of clusters %d from %d digits",fNumberOfECAClusters,fDigitsArr->GetEntriesFast())); 
}
