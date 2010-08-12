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
//  root [0] AliEMCALClusterizerNxN * cl = new AliEMCALClusterizerNxN("galice.root")  
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

// --- ROOT system ---
#include <TMath.h> 
#include <TMinuit.h>
#include <TTree.h> 
#include <TBenchmark.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TClonesArray.h>

// --- Standard library ---
#include <cassert>

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliEMCALClusterizerNxN.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALCalibData.h"
#include "AliESDCaloCluster.h"

ClassImp(AliEMCALClusterizerNxN)

Bool_t AliEMCALClusterizerNxN::fgkIsInputCalibrated = kFALSE;

//____________________________________________________________________________
AliEMCALClusterizerNxN::AliEMCALClusterizerNxN()
  : AliEMCALClusterizer()
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
}

//____________________________________________________________________________
AliEMCALClusterizerNxN::AliEMCALClusterizerNxN(AliEMCALGeometry* geometry)
  : AliEMCALClusterizer(geometry)
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
  // use this contructor to avoid usage of Init() which uses runloader
  // change needed by HLT - MP

}

//____________________________________________________________________________
AliEMCALClusterizerNxN::AliEMCALClusterizerNxN(AliEMCALGeometry* geometry, AliEMCALCalibData * calib, AliCaloCalibPedestal * caloped)
: AliEMCALClusterizer(geometry, calib, caloped)
{
	// ctor, geometry and calibration are initialized elsewhere.
				
}


//____________________________________________________________________________
  AliEMCALClusterizerNxN::~AliEMCALClusterizerNxN()
{
  // dtor
}


//____________________________________________________________________________
void AliEMCALClusterizerNxN::Digits2Clusters(Option_t * option)
{
  // Steering method to perform clusterization for the current event 
  // in AliEMCALLoader

  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALClusterizer"); 
  
  if(strstr(option,"print"))
    Print("") ; 
 
  //Get calibration parameters from file or digitizer default values.
  GetCalibrationParameters() ;

  //Get dead channel map from file or digitizer default values.
  GetCaloCalibPedestal() ;
	
  fNumberOfECAClusters = 0;

  MakeClusters() ;  //only the real clusters

  if(fToUnfold)
    MakeUnfolding() ;

  Int_t index ;

  //Evaluate position, dispersion and other RecPoint properties for EC section                      
  for(index = 0; index < fRecPoints->GetEntries(); index++) 
    {
      dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(index))->EvalAll(fECAW0,fDigitsArr) ;
      AliDebug(5, Form("MAX INDEX %d ", dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(index))->GetMaximalEnergyIndex()));
      //For each rec.point set the distance to the nearest bad crystal
      dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(index))->EvalDistanceToBadChannels(fCaloPed);
    }
  
  fRecPoints->Sort() ;

  for(index = 0; index < fRecPoints->GetEntries(); index++) 
    {
      (dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(index)))->SetIndexInList(index) ;
      (dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(index)))->Print();
    }
  
  fTreeR->Fill();
  
  if(strstr(option,"deb") || strstr(option,"all"))  
    PrintRecPoints(option) ;

  AliDebug(1,Form("EMCAL Clusterizer found %d Rec Points",fRecPoints->GetEntriesFast()));

  fRecPoints->Delete();

  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALClusterizer");
    printf("Exec took %f seconds for Clusterizing", 
	   gBenchmark->GetCpuTime("EMCALClusterizer"));
  }    
}

//____________________________________________________________________________
Int_t AliEMCALClusterizerNxN::AreNeighbours(AliEMCALDigit * d1, AliEMCALDigit * d2, Bool_t & shared) const
{
	// Gives the neighbourness of two digits = 0 are not neighbour ; continue searching 
	//                                       = 1 are neighbour
	//                                       = 2 is in different SM; continue searching 
	// In case it is in different SM, but same phi rack, check if neigbours at eta=0
	// neighbours are defined as digits having at least a common side 
	// The order of d1 and d2 is important: first (d1) should be a digit already in a cluster 
	//                                      which is compared to a digit (d2)  not yet in a cluster  
	
	static Int_t nSupMod1=0, nModule1=0, nIphi1=0, nIeta1=0, iphi1=0, ieta1=0;
	static Int_t nSupMod2=0, nModule2=0, nIphi2=0, nIeta2=0, iphi2=0, ieta2=0;
	static Int_t rowdiff=0, coldiff=0;

	shared = kFALSE;
	
	fGeom->GetCellIndex(d1->GetId(), nSupMod1,nModule1,nIphi1,nIeta1);
	fGeom->GetCellIndex(d2->GetId(), nSupMod2,nModule2,nIphi2,nIeta2);
	fGeom->GetCellPhiEtaIndexInSModule(nSupMod1,nModule1,nIphi1,nIeta1, iphi1,ieta1);
	fGeom->GetCellPhiEtaIndexInSModule(nSupMod2,nModule2,nIphi2,nIeta2, iphi2,ieta2);
	
	//If different SM, check if they are in the same phi, then consider cells close to eta=0 as neighbours; May 2010
	if(nSupMod1 != nSupMod2 ) 
	  {
	    //Check if the 2 SM are in the same PHI position (0,1), (2,3), ...
	    Float_t smPhi1 = fGeom->GetEMCGeometry()->GetPhiCenterOfSM(nSupMod1);
	    Float_t smPhi2 = fGeom->GetEMCGeometry()->GetPhiCenterOfSM(nSupMod2);
	    
	    if(!TMath::AreEqualAbs(smPhi1, smPhi2, 1e-3)) return 2; //Not phi rack equal, not neighbours
	    
	    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2
	    // C Side impair SM, nSupMod%2=1; A side pair SM nSupMod%2=0
	    if(nSupMod1%2) ieta1+=AliEMCALGeoParams::fgkEMCALCols;
	    else           ieta2+=AliEMCALGeoParams::fgkEMCALCols;
	    
	    shared = kTRUE; // maybe a shared cluster, we know this later, set it for the moment.
	    
	}//Different SM, same phi
	
	rowdiff = TMath::Abs(iphi1 - iphi2);  
	coldiff = TMath::Abs(ieta1 - ieta2) ;  

	// neighbours +-1 in col and row
	if ( TMath::Abs(coldiff) < 2 && TMath::Abs(rowdiff) < 2)
	  {
	    
	    AliDebug(9, Form("AliEMCALClusterizerNxN::AreNeighbours(): id1=%d, (row %d, col %d) ; id2=%d, (row %d, col %d), shared %d \n",
			     d1->GetId(), iphi1,ieta1, d2->GetId(), iphi2,ieta2, shared));
	    
	    return 1;
	  }//Neighbours
	else 
	  {
	    AliDebug(9, Form("NOT AliEMCALClusterizerNxN::AreNeighbours(): id1=%d, (row %d, col %d) ; id2=%d, (row %d, col %d), shared %d \n",
			     d1->GetId(), iphi1,ieta1, d2->GetId(), iphi2,ieta2, shared));
	    shared = kFALSE;
	    return 2 ; 
	  }//Not neighbours
}

//____________________________________________________________________________
void AliEMCALClusterizerNxN::MakeClusters()
{
  // Steering method to construct the clusters stored in a list of Reconstructed Points
  // A cluster is defined as a list of neighbour digits
  // Mar 03, 2007 by PAI

  if (fGeom==0) AliFatal("Did not get geometry from EMCALLoader");

  fRecPoints->Clear();

  // Set up TObjArray with pointers to digits to work on 
  //TObjArray *digitsC = new TObjArray();
  TObjArray digitsC;
  TIter nextdigit(fDigitsArr);
  AliEMCALDigit *digit = 0;
  while ( (digit = dynamic_cast<AliEMCALDigit*>(nextdigit())) ) {
    digitsC.AddLast(digit);
  }

  TIter nextdigitC(&digitsC);

  AliDebug(1,Form("MakeClusters: Number of digits %d  -> (e %f)\n",
		  fDigitsArr->GetEntries(),fMinECut));

  Bool_t bDone = kFALSE;
  while ( bDone != kTRUE )
    {
      //first sort the digits:
      Int_t iMaxEnergyDigit = -1;
      Float_t dMaxEnergyDigit = -1;
      AliEMCALDigit *pMaxEnergyDigit = 0;
      nextdigitC.Reset();
      while ( (digit = dynamic_cast<AliEMCALDigit *>(nextdigitC())) ) 
	{ // scan over the list of digitsC
	  Float_t dEnergyCalibrated = Calibrate(digit->GetAmplitude(), digit->GetTime(),digit->GetId());
	  //AliDebug(5, Form("-> Digit ENERGY: %1.5f", dEnergyCalibrated));
	  
	  //if(fGeom->CheckAbsCellId(digit->GetId()) && dEnergyCalibrated > fECAClusteringThreshold  )
	  if(fGeom->CheckAbsCellId(digit->GetId()) && dEnergyCalibrated > 0.0) // no threshold!
	    {
	      if (dEnergyCalibrated > dMaxEnergyDigit)
		{
		  dMaxEnergyDigit = dEnergyCalibrated;
		  iMaxEnergyDigit = digit->GetId();
		  pMaxEnergyDigit = digit;
		}
	    }
	}

      if (iMaxEnergyDigit < 0 || digitsC.GetEntries() <= 0) 
	{
	  bDone = kTRUE;
	  continue;
	}

      AliDebug (2, Form("Max digit found: %1.2f AbsId: %d", dMaxEnergyDigit, iMaxEnergyDigit));
      AliDebug(5, Form("Max Digit ENERGY: %1.5f", dMaxEnergyDigit));

      // keep the candidate digits in a list
      TList clusterDigitList;
      clusterDigitList.SetOwner(kFALSE);
      clusterDigitList.AddLast(pMaxEnergyDigit);	 

      Double_t clusterCandidateEnergy = dMaxEnergyDigit;

      // now loop over the resto of the digits and cluster into NxN cluster 
      // we do not actually cluster yet: we keep them in the list clusterDigitList
      nextdigitC.Reset();
      while ( (digit = dynamic_cast<AliEMCALDigit *>(nextdigitC())) ) 
	{ // scan over the list of digitsC
	  if (digit == pMaxEnergyDigit) continue;
	  Float_t dEnergyCalibrated = Calibrate(digit->GetAmplitude(), digit->GetTime(),digit->GetId());
	  AliDebug(5, Form("-> Digit ENERGY: %1.5f", dEnergyCalibrated));
	  //if(fGeom->CheckAbsCellId(digit->GetId()) && dEnergyCalibrated > fECAClusteringThreshold  )
	  if(fGeom->CheckAbsCellId(digit->GetId()) && dEnergyCalibrated > 0.0  )
	    {
	      Float_t time = pMaxEnergyDigit->GetTime(); //Time or TimeR?
	      if(TMath::Abs(time - digit->GetTime()) > fTimeCut ) continue; //Time or TimeR?
	      Bool_t shared = kFALSE; //cluster shared by 2 SuperModules?
	      if (AreNeighbours(pMaxEnergyDigit, digit, shared) == 1) // call (digit,digitN) in THAT order !!!!! 
		{      
		  clusterDigitList.AddLast(digit) ;
		  clusterCandidateEnergy += dEnergyCalibrated;
		}
	    }
	}// loop over the next digits

      // start a cluster here only if a cluster energy is larger than clustering threshold
      //if (clusterCandidateEnergy > 0.1)
      AliDebug(5, Form("Clusterization threshold is %f MeV", fECAClusteringThreshold));
      if (clusterCandidateEnergy > fECAClusteringThreshold)
	{
	  if(fNumberOfECAClusters >= fRecPoints->GetSize()) fRecPoints->Expand(2*fNumberOfECAClusters+1) ;
      
	  AliEMCALRecPoint *recPoint = new  AliEMCALRecPoint("") ; 
	  fRecPoints->AddAt(recPoint, fNumberOfECAClusters) ;
	  recPoint = dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(fNumberOfECAClusters)) ; 
	  fNumberOfECAClusters++ ;       
	  recPoint->SetClusterType(AliESDCaloCluster::kEMCALClusterv1);

	  AliDebug(9, Form("Number of cells per cluster (max is 9!): %d", clusterDigitList.GetEntries()));
	  for (Int_t idig = 0; idig < clusterDigitList.GetEntries(); idig++)
	    {

	      digit = (AliEMCALDigit*)clusterDigitList.At(idig);
	      Float_t dEnergyCalibrated = Calibrate(digit->GetAmplitude(), digit->GetTime(),digit->GetId());
	      AliDebug(5, Form(" Adding digit %d", digit->GetId()));
	      // note: this way the sharing info is lost!
	      recPoint->AddDigit(*digit, dEnergyCalibrated, kFALSE) ; //Time or TimeR?
	      digitsC.Remove(digit); 		  
	    }
	}
      else
	{
	  // we do not want to start clustering in the same spot!
	  // but in this case we may NOT reuse this seed for another cluster!
	  // need a better bookeeping?
	  digitsC.Remove(pMaxEnergyDigit);
	}

      AliDebug (2, Form("Number of digits left: %d", digitsC.GetEntries()));      
    } // while ! done 

  //delete digitsC ; //nope we use an object
  
  AliDebug(1,Form("total no of clusters %d from %d digits",fNumberOfECAClusters,fDigitsArr->GetEntriesFast())); 
}

//____________________________________________________________________________
void AliEMCALClusterizerNxN::MakeUnfolding()
{
  // Unfolds clusters using the shape of an ElectroMagnetic shower
  // Performs unfolding of all clusters
		
  if(fNumberOfECAClusters > 0){
    if (fGeom==0)
      AliFatal("Did not get geometry from EMCALLoader") ;
    Int_t nModulesToUnfold = fGeom->GetNCells();

    Int_t numberofNotUnfolded = fNumberOfECAClusters ;
    Int_t index ;
    for(index = 0 ; index < numberofNotUnfolded ; index++){

      AliEMCALRecPoint * recPoint = dynamic_cast<AliEMCALRecPoint *>( fRecPoints->At(index) ) ;

      TVector3 gpos;
      Int_t absId;
      recPoint->GetGlobalPosition(gpos);
      fGeom->GetAbsCellIdFromEtaPhi(gpos.Eta(),gpos.Phi(),absId);
      if(absId > nModulesToUnfold)
        break ;

      Int_t nMultipl = recPoint->GetMultiplicity() ;
      AliEMCALDigit ** maxAt = new AliEMCALDigit*[nMultipl] ;
      Float_t * maxAtEnergy = new Float_t[nMultipl] ;
      Int_t nMax = recPoint->GetNumberOfLocalMax(maxAt, maxAtEnergy,fECALocMaxCut,fDigitsArr) ;

      if( nMax > 1 ) {     // if cluster is very flat (no pronounced maximum) then nMax = 0
        //UnfoldCluster(recPoint, nMax, maxAt, maxAtEnergy) ;
        fRecPoints->Remove(recPoint);
        fRecPoints->Compress() ;
        index-- ;
        fNumberOfECAClusters-- ;
        numberofNotUnfolded-- ;
      }
      else{
        recPoint->SetNExMax(1) ; //Only one local maximum
      }

      delete[] maxAt ;
      delete[] maxAtEnergy ;
    }
  }
  // End of Unfolding of clusters
}

//____________________________________________________________________________
Double_t  AliEMCALClusterizerNxN::ShowerShape(Double_t x, Double_t y)
{ 
  // Shape of the shower
  // If you change this function, change also the gradient evaluation in ChiSquare()
  //
  Double_t r = sqrt(x*x+y*y);
  Double_t r133  = TMath::Power(r, 1.33) ;
  Double_t r669  = TMath::Power(r, 6.69) ;
  Double_t shape = TMath::Exp( -r133 * (1. / (1.57 + 0.0860 * r133) - 0.55 / (1 + 0.000563 * r669) ) ) ;
  return shape ;
}

//____________________________________________________________________________
void  AliEMCALClusterizerNxN::UnfoldCluster(AliEMCALRecPoint * /*iniTower*/, 
					    Int_t /*nMax*/, 
					    AliEMCALDigit ** /*maxAt*/, 
					    Float_t * /*maxAtEnergy*/)
{
  //
  // Performs the unfolding of a cluster with nMax overlapping showers 
  //
  AliWarning("Not implemented. To be.");
}

//___________________________________________________________________
void   AliEMCALClusterizerNxN::SetInputCalibrated(Bool_t val)
{
  //
  // input is calibrated - the case when we run already on ESD
  //
  AliEMCALClusterizerNxN::fgkIsInputCalibrated = val;
}
