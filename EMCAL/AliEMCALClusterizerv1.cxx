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
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new  IO (à la PHOS)
//  Mar 2007, Aleksei Pavlinov - new algoritmh of pseudo clusters
//////////////////////////////////////////////////////////////////////////////
//  Clusterization class. Performs clusterization (collects neighbouring active cells) and 
//  unfolds the clusters having several local maxima.  
//  Results are stored in TreeR#, branches EMCALTowerRP (EMC recPoints),
//  EMCALPreShoRP (CPV RecPoints) and AliEMCALClusterizer (Clusterizer with all 
//  parameters including input digits branch title, thresholds etc.)
//  This TTask is normally called from Reconstructioner, but can as well be used in 
//  standalone mode.
// Use Case:
//  root [0] AliEMCALClusterizerv1 * cl = new AliEMCALClusterizerv1("galice.root")  
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
#include <cassert>

class TROOT;
#include <TH1.h>
#include <TFile.h> 
class TFolder;
#include <TMath.h> 
#include <TMinuit.h>
#include <TTree.h> 
class TSystem; 
#include <TBenchmark.h>
#include <TBrowser.h>
#include <TROOT.h>

// --- Standard library ---


// --- AliRoot header files ---
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliESD.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCAL.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecParam.h"
#include "AliEMCALReconstructor.h"
#include "AliCDBManager.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALCalibData.h"
class AliCDBStorage;
#include "AliCDBEntry.h"

ClassImp(AliEMCALClusterizerv1)

//____________________________________________________________________________
AliEMCALClusterizerv1::AliEMCALClusterizerv1()
  : AliEMCALClusterizer(),
    fGeom(0),
    fDefaultInit(kFALSE),
    fToUnfold(kFALSE),
    fNumberOfECAClusters(0),fCalibData(0),fCaloPed(0),
    fADCchannelECA(0.),fADCpedestalECA(0.),fECAClusteringThreshold(0.),fECALocMaxCut(0.),
    fECAW0(0.),fTimeCut(1.),fTimeMin(-1.),fTimeMax(1.),fMinECut(0.)
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
  
  Init() ;
}

//____________________________________________________________________________
AliEMCALClusterizerv1::AliEMCALClusterizerv1(AliEMCALGeometry* geometry)
  : AliEMCALClusterizer(),
    fGeom(geometry),
    fDefaultInit(kFALSE),
    fToUnfold(kFALSE),
    fNumberOfECAClusters(0),fCalibData(0), fCaloPed(0),
    fADCchannelECA(0.),fADCpedestalECA(0.),fECAClusteringThreshold(0.),fECALocMaxCut(0.),
    fECAW0(0.),fTimeCut(1.),fTimeMin(-1.),fTimeMax(1.),fMinECut(0.)
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
  // use this contructor to avoid usage of Init() which uses runloader
  // change needed by HLT - MP

  // Note for the future: the use on runloader should be avoided or optional at least
  // another way is to make Init virtual and protected at least such that the deriving classes can overload
  // Init() ;
  //

  if (!fGeom)
    {
      AliFatal("Geometry not initialized.");
    }

  if(!gMinuit)
    gMinuit = new TMinuit(100) ;

}

//____________________________________________________________________________
AliEMCALClusterizerv1::AliEMCALClusterizerv1(AliEMCALGeometry* geometry, AliEMCALCalibData * calib, AliCaloCalibPedestal * caloped)
: AliEMCALClusterizer(),
fGeom(geometry),
fDefaultInit(kFALSE),
fToUnfold(kFALSE),
fNumberOfECAClusters(0),fCalibData(calib), fCaloPed(caloped),
fADCchannelECA(0.),fADCpedestalECA(0.),fECAClusteringThreshold(0.),fECALocMaxCut(0.),
fECAW0(0.),fTimeCut(1.),fTimeMin(-1.),fTimeMax(1.),fMinECut(0.)
{
	// ctor, geometry and calibration are initialized elsewhere.
	
	if (!fGeom)
		AliFatal("Geometry not initialized.");
		
	if(!gMinuit)
		gMinuit = new TMinuit(100) ;
	
}


//____________________________________________________________________________
  AliEMCALClusterizerv1::~AliEMCALClusterizerv1()
{
  // dtor
}

//____________________________________________________________________________
Float_t  AliEMCALClusterizerv1::Calibrate(const Float_t amp, const Float_t time, const Int_t absId) 
{
 
  // Convert digitized amplitude into energy.
  // Calibration parameters are taken from calibration data base for raw data,
  // or from digitizer parameters for simulated data.

  if(fCalibData){
    
    if (fGeom==0)
      AliFatal("Did not get geometry from EMCALLoader") ;
    
    Int_t iSupMod = -1;
    Int_t nModule = -1;
    Int_t nIphi   = -1;
    Int_t nIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    
    Bool_t bCell = fGeom->GetCellIndex(absId, iSupMod, nModule, nIphi, nIeta) ;
    if(!bCell) {
      fGeom->PrintGeometry();
      Error("Calibrate()"," Wrong cell id number : %i", absId);
      assert(0);
    }

    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,nModule,nIphi, nIeta,iphi,ieta);
	  
	// Check if channel is bad (dead or hot), in this case return 0.	
	// Gustavo: 15-12-09 In case of RAW data this selection is already done, but not in simulation.
	// for the moment keep it here but remember to do the selection at the sdigitizer level 
	// and remove it from here
	Int_t channelStatus = (Int_t)(fCaloPed->GetDeadMap(iSupMod))->GetBinContent(ieta,iphi);
	if(channelStatus == AliCaloCalibPedestal::kHot || channelStatus == AliCaloCalibPedestal::kDead) {
		  AliDebug(2,Form("Tower from SM %d, ieta %d, iphi %d is BAD : status %d !!!",iSupMod,ieta,iphi, channelStatus));
		  return 0;
	}
	//Check if time is too large or too small, indication of a noisy channel, remove in this case
	if(time > fTimeMax || time < fTimeMin) return 0;
	  
    fADCchannelECA  = fCalibData->GetADCchannel (iSupMod,ieta,iphi);
    fADCpedestalECA = fCalibData->GetADCpedestal(iSupMod,ieta,iphi);
  
   return -fADCpedestalECA + amp * fADCchannelECA ;        
 
  }
  else //Return energy with default parameters if calibration is not available
    return -fADCpedestalECA + amp * fADCchannelECA ; 
  
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::Digits2Clusters(Option_t * option)
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
  for(index = 0; index < fRecPoints->GetEntries(); index++) {
      dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(index))->EvalAll(fECAW0,fDigitsArr) ;
	  //For each rec.point set the distance to the nearest bad crystal
	  dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(index))->EvalDistanceToBadChannels(fCaloPed);
  }

  fRecPoints->Sort() ;

  for(index = 0; index < fRecPoints->GetEntries(); index++) {
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
Bool_t AliEMCALClusterizerv1::FindFit(AliEMCALRecPoint * recPoint, AliEMCALDigit ** maxAt, 
				      const Float_t* maxAtEnergy,
				      Int_t nPar, Float_t * fitparameters) const
{
  // Calls TMinuit to fit the energy distribution of a cluster with several maxima
  // The initial values for fitting procedure are set equal to the
  // positions of local maxima.       
  // Cluster will be fitted as a superposition of nPar/3
  // electromagnetic showers

  if (fGeom==0) AliFatal("Did not get geometry from EMCALLoader");

  gMinuit->mncler();                     // Reset Minuit's list of paramters
  gMinuit->SetPrintLevel(-1) ;           // No Printout
  gMinuit->SetFCN(AliEMCALClusterizerv1::UnfoldingChiSquare) ;
  // To set the address of the minimization function
  TList * toMinuit = new TList();
  toMinuit->AddAt(recPoint,0) ;
  toMinuit->AddAt(fDigitsArr,1) ;
  toMinuit->AddAt(fGeom,2) ;

  gMinuit->SetObjectFit(toMinuit) ;         // To tranfer pointer to UnfoldingChiSquare

  // filling initial values for fit parameters
  AliEMCALDigit * digit ;

  Int_t ierflg  = 0;
  Int_t index   = 0 ;
  Int_t nDigits = (Int_t) nPar / 3 ;

  Int_t iDigit ;

  for(iDigit = 0; iDigit < nDigits; iDigit++){
    digit = maxAt[iDigit];
    Double_t x = 0.;
    Double_t y = 0.;
    Double_t z = 0.;

    fGeom->RelPosCellInSModule(digit->GetId(), y, x, z);

    Float_t energy = maxAtEnergy[iDigit] ;

    gMinuit->mnparm(index, "x",  x, 0.1, 0, 0, ierflg) ;
    index++ ;
    if(ierflg != 0){
      Error("FindFit", "EMCAL Unfolding  Unable to set initial value for fit procedure : x = %f", x ) ;
      return kFALSE;
    }
    gMinuit->mnparm(index, "z",  z, 0.1, 0, 0, ierflg) ;
    index++ ;
    if(ierflg != 0){
      Error("FindFit", "EMCAL Unfolding  Unable to set initial value for fit procedure : z = %f", z) ;
      return kFALSE;
    }
    gMinuit->mnparm(index, "Energy",  energy , 0.05*energy, 0., 4.*energy, ierflg) ;
    index++ ;
    if(ierflg != 0){
      Error("FindFit", "EMCAL Unfolding  Unable to set initial value for fit procedure : energy = %f", energy) ;
      return kFALSE;
    }
  }

  Double_t p0 = 0.1 ; // "Tolerance" Evaluation stops when EDM = 0.0001*p0 ; 
                      // The number of function call slightly depends on it.
  //Double_t p1 = 1.0 ;
  Double_t p2 = 0.0 ;

  gMinuit->mnexcm("SET STR", &p2, 0, ierflg) ;   // force TMinuit to reduce function calls
  //  gMinuit->mnexcm("SET GRA", &p1, 1, ierflg) ;   // force TMinuit to use my gradient
  gMinuit->SetMaxIterations(5);
  gMinuit->mnexcm("SET NOW", &p2 , 0, ierflg) ;  // No Warnings
  gMinuit->mnexcm("MIGRAD", &p0, 0, ierflg) ;    // minimize

  if(ierflg == 4){  // Minimum not found
    Error("FindFit", "EMCAL Unfolding  Fit not converged, cluster abandoned " ) ;
    return kFALSE ;
  }
  for(index = 0; index < nPar; index++){
    Double_t err ;
    Double_t val ;
    gMinuit->GetParameter(index, val, err) ;    // Returns value and error of parameter index
    fitparameters[index] = val ;
  }

  delete toMinuit ;
  return kTRUE;

}

//____________________________________________________________________________
void AliEMCALClusterizerv1::GetCalibrationParameters() 
{
  // Set calibration parameters:
  // if calibration database exists, they are read from database,
  // otherwise, they are taken from digitizer.
  //
  // It is a user responsilibity to open CDB before reconstruction, 
  // for example: 
  // AliCDBStorage* storage = AliCDBManager::Instance()->GetStorage("local://CalibDB");

  //Check if calibration is stored in data base

  if(!fCalibData)
    {
      AliCDBEntry *entry = (AliCDBEntry*) 
	AliCDBManager::Instance()->Get("EMCAL/Calib/Data");
      if (entry) fCalibData =  (AliEMCALCalibData*) entry->GetObject();
    }
  
  if(!fCalibData)
    AliFatal("Calibration parameters not found in CDB!");
 
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::GetCaloCalibPedestal() 
{
	// Set calibration parameters:
	// if calibration database exists, they are read from database,
	// otherwise, they are taken from digitizer.
	//
	// It is a user responsilibity to open CDB before reconstruction, 
	// for example: 
	// AliCDBStorage* storage = AliCDBManager::Instance()->GetStorage("local://CalibDB");
	
	//Check if calibration is stored in data base
	
	if(!fCaloPed)
    {
		AliCDBEntry *entry = (AliCDBEntry*) 
		AliCDBManager::Instance()->Get("EMCAL/Calib/Pedestals");
		if (entry) fCaloPed =  (AliCaloCalibPedestal*) entry->GetObject();
    }
	
	if(!fCaloPed)
		AliFatal("Pedestal info not found in CDB!");
	
}


//____________________________________________________________________________
void AliEMCALClusterizerv1::Init()
{
  // Make all memory allocations which can not be done in default constructor.
  // Attach the Clusterizer task to the list of EMCAL tasks
  
  AliRunLoader *rl = AliRunLoader::Instance();
  if (rl->GetAliRun() && rl->GetAliRun()->GetDetector("EMCAL"))
    fGeom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
  else 
    fGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());

  AliDebug(1,Form("geom 0x%x",fGeom));

  if(!gMinuit) 
    gMinuit = new TMinuit(100) ;

}

//____________________________________________________________________________
void AliEMCALClusterizerv1::InitParameters()
{ 
  // Initializes the parameters for the Clusterizer
  fNumberOfECAClusters = 0;

  fCalibData               = 0 ;
  fCaloPed                 = 0 ;
	
  const AliEMCALRecParam* recParam = AliEMCALReconstructor::GetRecParam();
  if(!recParam) {
    AliFatal("Reconstruction parameters for EMCAL not set!");
  } else {
    fECAClusteringThreshold = recParam->GetClusteringThreshold();
    fECAW0                  = recParam->GetW0();
    fMinECut                = recParam->GetMinECut();    
    fToUnfold               = recParam->GetUnfold();
    if(fToUnfold) AliWarning("Cluster Unfolding ON. Implementing only for eta=0 case!!!"); 
    fECALocMaxCut           = recParam->GetLocMaxCut();
    fTimeCut                = recParam->GetTimeCut();
    fTimeMin                = recParam->GetTimeMin();
    fTimeMax                = recParam->GetTimeMax();

    AliDebug(1,Form("Reconstruction parameters: fECAClusteringThreshold=%.3f GeV, fECAW=%.3f, fMinECut=%.3f GeV, fToUnfold=%d, fECALocMaxCut=%.3f GeV, fTimeCut=%e s,fTimeMin=%e s,fTimeMax=%e s",
		    fECAClusteringThreshold,fECAW0,fMinECut,fToUnfold,fECALocMaxCut,fTimeCut, fTimeMin, fTimeMax));
  }

}

//____________________________________________________________________________
Int_t AliEMCALClusterizerv1::AreNeighbours(AliEMCALDigit * d1, AliEMCALDigit * d2) const
{
  // Gives the neighbourness of two digits = 0 are not neighbour ; continue searching 
  //                                       = 1 are neighbour
  //                                       = 2 is in different SM; continue searching 
  // neighbours are defined as digits having at least a common vertex 
  // The order of d1 and d2 is important: first (d1) should be a digit already in a cluster 
  //                                      which is compared to a digit (d2)  not yet in a cluster  

  static Int_t rv; 
  static Int_t nSupMod1=0, nModule1=0, nIphi1=0, nIeta1=0, iphi1=0, ieta1=0;
  static Int_t nSupMod2=0, nModule2=0, nIphi2=0, nIeta2=0, iphi2=0, ieta2=0;
  static Int_t rowdiff, coldiff;
  rv = 0 ; 

  fGeom->GetCellIndex(d1->GetId(), nSupMod1,nModule1,nIphi1,nIeta1);
  fGeom->GetCellIndex(d2->GetId(), nSupMod2,nModule2,nIphi2,nIeta2);

  // Do not aggregate cells in different SM
  if(nSupMod1 != nSupMod2 ) return  2;
  
  fGeom->GetCellPhiEtaIndexInSModule(nSupMod1,nModule1,nIphi1,nIeta1, iphi1,ieta1);
  fGeom->GetCellPhiEtaIndexInSModule(nSupMod2,nModule2,nIphi2,nIeta2, iphi2,ieta2);
  
  rowdiff = TMath::Abs(iphi1 - iphi2);  
  coldiff = TMath::Abs(ieta1 - ieta2) ;  
  
  // neighbours in same SM with at least common side; May 11, 2007
  if ((coldiff==0 && TMath::Abs(rowdiff)==1) || (rowdiff==0 && TMath::Abs(coldiff)==1)) rv = 1;  

  //Diagonal?
  //if ((coldiff==0 && TMath::Abs(rowdiff==1)) || (rowdiff==0 && TMath::Abs(coldiff==1)) || (TMath::Abs(rowdiff)==1 && TMath::Abs(coldiff==1))) rv = 1;
	
  if (gDebug == 2 && rv==1) 
  printf("AreNeighbours: id1=%d, (%d,%d) \n id2=%d, (%d,%d) \n",d1->GetId(), iphi1,ieta1, d2->GetId(), iphi2,ieta2);   
  
  return rv ; 
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::MakeClusters()
{
  // Steering method to construct the clusters stored in a list of Reconstructed Points
  // A cluster is defined as a list of neighbour digits
  // Mar 03, 2007 by PAI

  if (fGeom==0) AliFatal("Did not get geometry from EMCALLoader");

  fRecPoints->Clear();

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
  AliDebug(1,Form("MakeClusters: Number of digits %d  -> (e %f), ehs %d\n",
		  fDigitsArr->GetEntries(),fMinECut,ehs));

  nextdigitC.Reset();

  while ( (digit = dynamic_cast<AliEMCALDigit *>(nextdigitC())) ) { // scan over the list of digitsC
    TArrayI clusterECAdigitslist(fDigitsArr->GetEntries());

    if(fGeom->CheckAbsCellId(digit->GetId()) && (Calibrate(digit->GetAmplitude(), digit->GetTime(),digit->GetId()) > fECAClusteringThreshold  ) ){
      // start a new Tower RecPoint
      if(fNumberOfECAClusters >= fRecPoints->GetSize()) fRecPoints->Expand(2*fNumberOfECAClusters+1) ;

      AliEMCALRecPoint *recPoint = new  AliEMCALRecPoint("") ; 
      fRecPoints->AddAt(recPoint, fNumberOfECAClusters) ;
      recPoint = dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(fNumberOfECAClusters)) ; 
      fNumberOfECAClusters++ ; 

      recPoint->SetClusterType(AliESDCaloCluster::kEMCALClusterv1);

      recPoint->AddDigit(*digit, Calibrate(digit->GetAmplitude(), digit->GetTime(),digit->GetId())) ; //Time or TimeR?
      TObjArray clusterDigits;
      clusterDigits.AddLast(digit);	
      digitsC->Remove(digit) ; 

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
			if(TMath::Abs(time - digitN->GetTime()) > fTimeCut ) continue; //Time or TimeR?
			
			if (AreNeighbours(digit, digitN)==1) {      // call (digit,digitN) in THAT oder !!!!! 
		
				recPoint->AddDigit(*digitN, Calibrate(digitN->GetAmplitude(), digitN->GetTime(), digitN->GetId()) ) ;//Time or TimeR?
				clusterDigits.AddLast(digitN) ; 
				digitsC->Remove(digitN) ; 
			} // if(ineb==1)
		} // scan over digits
      } // scan over digits already in cluster
      if(recPoint)
        AliDebug(2,Form("MakeClusters: %d digitd, energy %f \n", clusterDigits.GetEntries(), recPoint->GetEnergy())); 
    } // If seed found
  } // while digit 

  delete digitsC ;
  
  AliDebug(1,Form("total no of clusters %d from %d digits",fNumberOfECAClusters,fDigitsArr->GetEntriesFast())); 
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::MakeUnfolding()
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
        UnfoldCluster(recPoint, nMax, maxAt, maxAtEnergy) ;
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
Double_t  AliEMCALClusterizerv1::ShowerShape(Double_t x, Double_t y)
{ 
  // Shape of the shower
  // If you change this function, change also the gradient evaluation in ChiSquare()

  Double_t r = sqrt(x*x+y*y);
  Double_t r133  = TMath::Power(r, 1.33) ;
  Double_t r669  = TMath::Power(r, 6.69) ;
  Double_t shape = TMath::Exp( -r133 * (1. / (1.57 + 0.0860 * r133) - 0.55 / (1 + 0.000563 * r669) ) ) ;
  return shape ;
}

//____________________________________________________________________________
void  AliEMCALClusterizerv1::UnfoldCluster(AliEMCALRecPoint * iniTower, 
					   Int_t nMax, 
					   AliEMCALDigit ** maxAt, 
					   Float_t * maxAtEnergy)
{
  // Performs the unfolding of a cluster with nMax overlapping showers 
  Int_t nPar = 3 * nMax ;
  Float_t * fitparameters = new Float_t[nPar] ;

  if (fGeom==0)
    AliFatal("Did not get geometry from EMCALLoader") ;

  Bool_t rv = FindFit(iniTower, maxAt, maxAtEnergy, nPar, fitparameters) ;
  if( !rv ) {
    // Fit failed, return and remove cluster
    iniTower->SetNExMax(-1) ;
    delete[] fitparameters ;
    return ;
  }

  // create unfolded rec points and fill them with new energy lists
  // First calculate energy deposited in each sell in accordance with
  // fit (without fluctuations): efit[]
  // and later correct this number in acordance with actual energy
  // deposition

  Int_t nDigits = iniTower->GetMultiplicity() ;
  Float_t * efit = new Float_t[nDigits] ;
  Double_t xDigit=0.,yDigit=0.,zDigit=0. ;
  Float_t xpar=0.,zpar=0.,epar=0.  ;

  AliEMCALDigit * digit = 0 ;
  Int_t * digitsList = iniTower->GetDigitsList() ;

  Int_t iparam ;
  Int_t iDigit ;
  for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){
    digit = dynamic_cast<AliEMCALDigit*>( fDigitsArr->At(digitsList[iDigit] ) ) ;
    fGeom->RelPosCellInSModule(digit->GetId(), yDigit, xDigit, zDigit);
    efit[iDigit] = 0;

    iparam = 0 ;
    while(iparam < nPar ){
      xpar = fitparameters[iparam] ;
      zpar = fitparameters[iparam+1] ;
      epar = fitparameters[iparam+2] ;
      iparam += 3 ;
      efit[iDigit] += epar * ShowerShape(xDigit - xpar,zDigit - zpar) ;
    }
  }


  // Now create new RecPoints and fill energy lists with efit corrected to fluctuations
  // so that energy deposited in each cell is distributed between new clusters proportionally
  // to its contribution to efit

  Float_t * energiesList = iniTower->GetEnergiesList() ;
  Float_t ratio ;

  iparam = 0 ;
  while(iparam < nPar ){
    xpar = fitparameters[iparam] ;
    zpar = fitparameters[iparam+1] ;
    epar = fitparameters[iparam+2] ;
    iparam += 3 ;

    AliEMCALRecPoint * recPoint = 0 ;

    if(fNumberOfECAClusters >= fRecPoints->GetSize())
      fRecPoints->Expand(2*fNumberOfECAClusters) ;

    (*fRecPoints)[fNumberOfECAClusters] = new AliEMCALRecPoint("") ;
    recPoint = dynamic_cast<AliEMCALRecPoint *>( fRecPoints->At(fNumberOfECAClusters) ) ;
    fNumberOfECAClusters++ ;
    recPoint->SetNExMax((Int_t)nPar/3) ;

    Float_t eDigit ;
    for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){
      digit = dynamic_cast<AliEMCALDigit*>( fDigitsArr->At( digitsList[iDigit] ) ) ;
      fGeom->RelPosCellInSModule(digit->GetId(), yDigit, xDigit, zDigit);

      ratio = epar * ShowerShape(xDigit - xpar,zDigit - zpar) / efit[iDigit] ;
      eDigit = energiesList[iDigit] * ratio ;
      recPoint->AddDigit( *digit, eDigit ) ;
    }
  }

  delete[] fitparameters ;
  delete[] efit ;

}

//_____________________________________________________________________________
void AliEMCALClusterizerv1::UnfoldingChiSquare(Int_t & nPar, Double_t * Grad,
					       Double_t & fret,
					       Double_t * x, Int_t iflag)
{
  // Calculates the Chi square for the cluster unfolding minimization
  // Number of parameters, Gradient, Chi squared, parameters, what to do

  TList * toMinuit = dynamic_cast<TList*>( gMinuit->GetObjectFit() ) ;

  AliEMCALRecPoint * recPoint = dynamic_cast<AliEMCALRecPoint*>( toMinuit->At(0) )  ;
  TClonesArray * digits = dynamic_cast<TClonesArray*>( toMinuit->At(1) )  ;
  // A bit buggy way to get an access to the geometry
  // To be revised!
  AliEMCALGeometry *geom = dynamic_cast<AliEMCALGeometry *>(toMinuit->At(2));

  Int_t * digitsList     = recPoint->GetDigitsList() ;

  Int_t nOdigits = recPoint->GetDigitsMultiplicity() ;

  Float_t * energiesList = recPoint->GetEnergiesList() ;

  fret = 0. ;
  Int_t iparam ;

  if(iflag == 2)
    for(iparam = 0 ; iparam < nPar ; iparam++)
      Grad[iparam] = 0 ; // Will evaluate gradient

  Double_t efit ;

  AliEMCALDigit * digit ;
  Int_t iDigit ;

  for( iDigit = 0 ; iDigit < nOdigits ; iDigit++) {

    digit = dynamic_cast<AliEMCALDigit*>( digits->At( digitsList[iDigit] ) );

    Double_t xDigit=0 ;
    Double_t zDigit=0 ;
    Double_t yDigit=0 ;//not used yet, assumed to be 0

    geom->RelPosCellInSModule(digit->GetId(), yDigit, xDigit, zDigit);

    if(iflag == 2){  // calculate gradient
      Int_t iParam = 0 ;
      efit = 0 ;
      while(iParam < nPar ){
        Double_t dx = (xDigit - x[iParam]) ;
        iParam++ ;
        Double_t dz = (zDigit - x[iParam]) ;
        iParam++ ;
        efit += x[iParam] * ShowerShape(dx,dz) ;
        iParam++ ;
      }
      Double_t sum = 2. * (efit - energiesList[iDigit]) / energiesList[iDigit] ; // Here we assume, that sigma = sqrt(E)
      iParam = 0 ;
      while(iParam < nPar ){
        Double_t xpar = x[iParam] ;
        Double_t zpar = x[iParam+1] ;
        Double_t epar = x[iParam+2] ;
        Double_t dr = TMath::Sqrt( (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar) );
        Double_t shape = sum * ShowerShape(xDigit - xpar,zDigit - zpar) ;
        Double_t r133 =  TMath::Power(dr, 1.33);
        Double_t r669 = TMath::Power(dr,6.69);
        Double_t deriv =-1.33 * TMath::Power(dr,0.33)*dr * ( 1.57 / ( (1.57 + 0.0860 * r133) * (1.57 + 0.0860 * r133) )
                                                             - 0.55 / (1 + 0.000563 * r669) / ( (1 + 0.000563 * r669) * (1 + 0.000563 * r669) ) ) ;

        Grad[iParam] += epar * shape * deriv * (xpar - xDigit) ;  // Derivative over x
        iParam++ ;
        Grad[iParam] += epar * shape * deriv * (zpar - zDigit) ;  // Derivative over z
        iParam++ ;
        Grad[iParam] += shape ;                                  // Derivative over energy
	iParam++ ;
      }
    }
    efit = 0;
    iparam = 0 ;


    while(iparam < nPar ){
      Double_t xpar = x[iparam] ;
      Double_t zpar = x[iparam+1] ;
      Double_t epar = x[iparam+2] ;
      iparam += 3 ;
      efit += epar * ShowerShape(xDigit - xpar,zDigit - zpar) ;
    }

    fret += (efit-energiesList[iDigit])*(efit-energiesList[iDigit])/energiesList[iDigit] ;
    // Here we assume, that sigma = sqrt(E) 
  }
}
//____________________________________________________________________________
void AliEMCALClusterizerv1::Print(Option_t * /*option*/)const
{
  // Print clusterizer parameters

  TString message("\n") ; 
  
  if( strcmp(GetName(), "") !=0 ){
    
    // Print parameters
 
    TString taskName(Version()) ;
    
    printf("--------------- "); 
    printf("%s",taskName.Data()) ; 
    printf(" "); 
    printf("Clusterizing digits: "); 
    printf("\n                       ECA Local Maximum cut    = %f", fECALocMaxCut); 
    printf("\n                       ECA Logarithmic weight   = %f", fECAW0); 
    if(fToUnfold)
      printf("\nUnfolding on\n");
    else
      printf("\nUnfolding off\n");
    
    printf("------------------------------------------------------------------"); 
  }
  else
    printf("AliEMCALClusterizerv1 not initialized ") ;
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::PrintRecPoints(Option_t * option)
{
  // Prints list of RecPoints produced at the current pass of AliEMCALClusterizer
  if(strstr(option,"deb")) {
    printf("PrintRecPoints: Clusterization result:") ; 
  
    printf("           Found %d ECA Rec Points\n ", 
	 fRecPoints->GetEntriesFast()) ; 
  }

  if(strstr(option,"all")) {
    if(strstr(option,"deb")) {
      printf("\n-----------------------------------------------------------------------\n") ;
      printf("Clusters in ECAL section\n") ;
      printf("Index    Ene(GeV) Multi Module     GX    GY   GZ  lX    lY   lZ   Dispersion Lambda 1   Lambda 2  # of prim  Primaries list\n") ;
    }
   Int_t index =0;

    for (index = 0 ; index < fRecPoints->GetEntries() ; index++) {
      AliEMCALRecPoint * rp = dynamic_cast<AliEMCALRecPoint * >(fRecPoints->At(index)) ; 
      TVector3  globalpos;  
      //rp->GetGlobalPosition(globalpos);
      TVector3  localpos;  
      rp->GetLocalPosition(localpos);
      Float_t lambda[2]; 
      rp->GetElipsAxis(lambda);
      Int_t * primaries; 
      Int_t nprimaries;
      primaries = rp->GetPrimaries(nprimaries);
      if(strstr(option,"deb")) 
      printf("\n%6d  %8.4f  %3d     %4.1f    %4.1f %4.1f  %4.1f %4.1f %4.1f    %4.1f   %4f  %4f    %2d     : ", 
	     rp->GetIndexInList(), rp->GetEnergy(), rp->GetMultiplicity(),
	     globalpos.X(), globalpos.Y(), globalpos.Z(), localpos.X(), localpos.Y(), localpos.Z(), 
	     rp->GetDispersion(), lambda[0], lambda[1], nprimaries) ; 
      if(strstr(option,"deb")){ 
        for (Int_t iprimary=0; iprimary<nprimaries; iprimary++) {
	  printf("%d ", primaries[iprimary] ) ; 
        }
      }
    }

    if(strstr(option,"deb"))
    printf("\n-----------------------------------------------------------------------\n");
  }
}

//___________________________________________________________________
void  AliEMCALClusterizerv1::PrintRecoInfo()
{
  printf(" AliEMCALClusterizerv1::PrintRecoInfo() : version %s \n", Version() );

}
