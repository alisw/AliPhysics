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

//*-- Author: Yves Schutz (SUBATECH)  & Dmitri Peressounko (SUBATECH & Kurchatov Institute)
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new  IO (à la PHOS)
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
#include "AliEMCALLoader.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCAL.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALHistoUtilities.h"
#include "AliCDBManager.h"

class AliCDBStorage;
#include "AliCDBEntry.h"

ClassImp(AliEMCALClusterizerv1)

//____________________________________________________________________________
AliEMCALClusterizerv1::AliEMCALClusterizerv1() 
  : AliEMCALClusterizer(),
    fHists(0),fPointE(0),fPointL1(0),fPointL2(0),
    fPointDis(0),fPointMult(0),fDigitAmp(0),fMaxE(0),
    fMaxL1(0),fMaxL2(0),fMaxDis(0),fGeom(0),
    fDefaultInit(kTRUE),
    fToUnfold(kFALSE),
    fNumberOfECAClusters(0),fNTowerInGroup(0),fCalibData(0),
    fADCchannelECA(0.),fADCpedestalECA(0.),fECAClusteringThreshold(0.),fECALocMaxCut(0.),
    fECAW0(0.),fRecPointsInRun(0),fTimeGate(0.),fMinECut(0.)
{
  // default ctor (to be used mainly by Streamer)
  
  InitParameters() ; 
  fGeom = AliEMCALGeometry::GetInstance();
  fGeom->GetTransformationForSM(); // Global <-> Local
}

//____________________________________________________________________________
AliEMCALClusterizerv1::AliEMCALClusterizerv1(const TString alirunFileName, const TString eventFolderName)
  : AliEMCALClusterizer(alirunFileName, eventFolderName),
    fHists(0),fPointE(0),fPointL1(0),fPointL2(0),
    fPointDis(0),fPointMult(0),fDigitAmp(0),fMaxE(0),
    fMaxL1(0),fMaxL2(0),fMaxDis(0),fGeom(0),
    fDefaultInit(kFALSE),
    fToUnfold(kFALSE),
    fNumberOfECAClusters(0),fNTowerInGroup(0),fCalibData(0),
    fADCchannelECA(0.),fADCpedestalECA(0.),fECAClusteringThreshold(0.),fECALocMaxCut(0.),
    fECAW0(0.),fRecPointsInRun(0),fTimeGate(0.),fMinECut(0.)
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
  
  InitParameters() ; 
  Init() ;
}

//____________________________________________________________________________
AliEMCALClusterizerv1::AliEMCALClusterizerv1(const AliEMCALClusterizerv1& clus)
  : AliEMCALClusterizer(clus),
    fHists(clus.fHists),
    fPointE(clus.fPointE),
    fPointL1(clus.fPointL1),
    fPointL2(clus.fPointL2),
    fPointDis(clus.fPointDis),
    fPointMult(clus.fPointMult),
    fDigitAmp(clus.fDigitAmp),
    fMaxE(clus.fMaxE),
    fMaxL1(clus.fMaxL1),
    fMaxL2(clus.fMaxL2),
    fMaxDis(clus.fMaxDis),
    fGeom(clus.fGeom),
    fDefaultInit(clus.fDefaultInit),
    fToUnfold(clus.fToUnfold),
    fNumberOfECAClusters(clus.fNumberOfECAClusters),
    fNTowerInGroup(clus.fNTowerInGroup),
    fCalibData(clus.fCalibData),
    fADCchannelECA(clus.fADCchannelECA),
    fADCpedestalECA(clus.fADCpedestalECA),
    fECAClusteringThreshold(clus.fECAClusteringThreshold),
    fECALocMaxCut(clus.fECALocMaxCut),
    fECAW0(clus.fECAW0),
    fRecPointsInRun(clus.fRecPointsInRun),
    fTimeGate(clus.fTimeGate),
    fMinECut(clus.fMinECut)
{
  //copy ctor
}

//____________________________________________________________________________
  AliEMCALClusterizerv1::~AliEMCALClusterizerv1()
{
  // dtor
}

//____________________________________________________________________________
const TString AliEMCALClusterizerv1::BranchName() const 
{ 
   return GetName();
}

//____________________________________________________________________________
Float_t  AliEMCALClusterizerv1::Calibrate(Int_t amp, Int_t AbsId) 
{
 
  // Convert digitized amplitude into energy.
  // Calibration parameters are taken from calibration data base for raw data,
  // or from digitizer parameters for simulated data.

  if(fCalibData){

    //JLK 13-Mar-2006
    //We now get geometry at a higher level
    //
    // Loader
    //    AliRunLoader *rl = AliRunLoader::GetRunLoader();
    
    // Load EMCAL Geomtry
    //    rl->LoadgAlice(); 
    //AliRun   * gAlice = rl->GetAliRun(); 
    //AliEMCAL * emcal  = (AliEMCAL*)gAlice->GetDetector("EMCAL");
    //AliEMCALGeometry * geom = emcal->GetGeometry();
    
    if (fGeom==0)
      AliFatal("Did not get geometry from EMCALLoader") ;
    
    Int_t iSupMod = -1;
    Int_t nModule  = -1;
    Int_t nIphi   = -1;
    Int_t nIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    
    Bool_t bCell = fGeom->GetCellIndex(AbsId, iSupMod, nModule, nIphi, nIeta) ;
    if(!bCell) {
      fGeom->PrintGeometry();
      Error("Calibrate()"," Wrong cell id number : %i", AbsId);
      assert(0);
    }
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,nModule,nIphi, nIeta,iphi,ieta);

    fADCchannelECA  = fCalibData->GetADCchannel (iSupMod,ieta,iphi);
    fADCpedestalECA = fCalibData->GetADCpedestal(iSupMod,ieta,iphi);
    return -fADCpedestalECA + amp * fADCchannelECA ;        
 
  }
  else //Return energy with default parameters if calibration is not available
    return -fADCpedestalECA + amp * fADCchannelECA ; 
  
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::Exec(Option_t * option)
{
  // Steering method to perform clusterization for events
  // in the range from fFirstEvent to fLastEvent.
  // This range is optionally set by SetEventRange().
  // if fLastEvent=-1 (by default), then process events until the end.

  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALClusterizer"); 
  
  if(strstr(option,"print"))
    Print("") ; 

  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));

  //Get calibration parameters from file or digitizer default values.
  GetCalibrationParameters() ;

  if (fLastEvent == -1) 
    fLastEvent = rl->GetNumberOfEvents() - 1;
  Int_t nEvents   = fLastEvent - fFirstEvent + 1;

  Int_t ievent ;
  rl->LoadDigits("EMCAL");
  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {
    rl->GetEvent(ievent);

    fNumberOfECAClusters = 0;

    if(strstr(option,"pseudo"))
      MakeClusters("pseudo") ;  //both types
    else
      MakeClusters("") ;  //only the real clusters

    if(fToUnfold)
      MakeUnfolding() ;

    WriteRecPoints() ;

    if(strstr(option,"deb") || strstr(option,"all"))  
      PrintRecPoints(option) ;

    //increment the total number of recpoints per run   
    fRecPointsInRun += emcalLoader->RecPoints()->GetEntriesFast() ;  
  }
  
  Unload();

  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALClusterizer");
    printf("Exec took %f seconds for Clusterizing %f seconds per event", 
	 gBenchmark->GetCpuTime("EMCALClusterizer"), gBenchmark->GetCpuTime("EMCALClusterizer")/nEvents );
  }

}

//____________________________________________________________________________
Bool_t AliEMCALClusterizerv1::FindFit(AliEMCALRecPoint * emcRP, AliEMCALDigit ** maxAt, Float_t * maxAtEnergy,
				    Int_t nPar, Float_t * fitparameters) const
{ 
  // Calls TMinuit to fit the energy distribution of a cluster with several maxima 
  // The initial values for fitting procedure are set equal to the positions of local maxima.
  // Cluster will be fitted as a superposition of nPar/3 electromagnetic showers

  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));
  TClonesArray *digits = emcalLoader->Digits();

  gMinuit->mncler();                     // Reset Minuit's list of paramters
  gMinuit->SetPrintLevel(-1) ;           // No Printout
  gMinuit->SetFCN(AliEMCALClusterizerv1::UnfoldingChiSquare) ;  
                                         // To set the address of the minimization function 
  TList * toMinuit = new TList();
  toMinuit->AddAt(emcRP,0) ;
  toMinuit->AddAt(digits,1) ;
  
  gMinuit->SetObjectFit(toMinuit) ;         // To tranfer pointer to UnfoldingChiSquare

  // filling initial values for fit parameters
  AliEMCALDigit * digit ;

  Int_t ierflg  = 0; 
  Int_t index   = 0 ;
  Int_t nDigits = (Int_t) nPar / 3 ;

  Int_t iDigit ;

  for(iDigit = 0; iDigit < nDigits; iDigit++){
    digit = maxAt[iDigit]; 

    Float_t x = 0.;
    Float_t z = 0.;
    //   have to be tune for TRD1; May 31,06
    //   Int_t relid[2] ;
    //   fGeom->AbsToRelNumbering(digit->GetId(), relid) ; // obsolete method
    //   fGeom->PosInAlice(relid, x, z) ;                  // obsolete method

    Float_t energy = maxAtEnergy[iDigit] ;

    gMinuit->mnparm(index, "x",  x, 0.1, 0, 0, ierflg) ;
    index++ ;   
    if(ierflg != 0){ 
      Error("FindFit", "EMCAL Unfolding  Unable to set initial value for fit procedure : x = %f",  x ) ;
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

  Double_t p0 = 0.1 ; // "Tolerance" Evaluation stops when EDM = 0.0001*p0 ; The number of function call slightly
                      //  depends on it. 
  Double_t p1 = 1.0 ;
  Double_t p2 = 0.0 ;

  gMinuit->mnexcm("SET STR", &p2, 0, ierflg) ;   // force TMinuit to reduce function calls  
  gMinuit->mnexcm("SET GRA", &p1, 1, ierflg) ;   // force TMinuit to use my gradient  
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
   if(AliCDBManager::Instance()->IsDefaultStorageSet()){
     AliCDBEntry *entry = (AliCDBEntry*) 
     AliCDBManager::Instance()->Get("EMCAL/Calib/Data");
     if (entry) fCalibData = (AliEMCALCalibData*) entry->GetObject();
   }
   if(!fCalibData)
     {
       //If calibration is not available use default parameters
       //Loader
       AliEMCALLoader *emcalLoader = 
	 dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));
       if ( !emcalLoader->Digitizer() ) 
	 emcalLoader->LoadDigitizer();
       AliEMCALDigitizer * dig = dynamic_cast<AliEMCALDigitizer*>(emcalLoader->Digitizer());
       
       fADCchannelECA   = dig->GetECAchannel() ;
       fADCpedestalECA  = dig->GetECApedestal();
     }
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::Init()
{
  // Make all memory allocations which can not be done in default constructor.
  // Attach the Clusterizer task to the list of EMCAL tasks
  
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  fGeom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
  fGeom->GetTransformationForSM(); // Global <-> Local
  AliInfo(Form("geom 0x%x",fGeom));

  if(!gMinuit) 
    gMinuit = new TMinuit(100) ;

  fHists = BookHists();
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::InitParameters()
{ 
  // Initializes the parameters for the Clusterizer
  fNumberOfECAClusters = 0;

  fNTowerInGroup = 36;  //Produces maximum of 80 pseudoclusters per event

  fECAClusteringThreshold   = 0.1;  // value obtained from Aleksei
  fECALocMaxCut = 0.03; // ??

  fECAW0    = 4.5;
  fTimeGate = 1.e-8 ; 
  fToUnfold = kFALSE ;
  fRecPointsInRun  = 0 ;
  fMinECut = 0.01; // have to be tune

  fCalibData               = 0 ;
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
  if(nSupMod1 != nSupMod2) return 2; // different SM

  fGeom->GetCellPhiEtaIndexInSModule(nSupMod1,nModule1,nIphi1,nIeta1, iphi1,ieta1);
  fGeom->GetCellPhiEtaIndexInSModule(nSupMod2,nModule2,nIphi2,nIeta2, iphi2,ieta2);

  rowdiff = TMath::Abs(iphi1 - iphi2);  
  coldiff = TMath::Abs(ieta1 - ieta2) ;  
  
  //  if (( coldiff <= 1 )  && ( rowdiff <= 1 )) rv = 1;  // neighbours with at least commom vertex
  if(coldiff + rowdiff <= 1) rv = 1;  // neighbours with at least commom side; Nov 6,2006
 
  if (gDebug == 2 && rv==1) 
  printf("AreNeighbours: neighbours=%d, id1=%d, relid1=%d,%d \n id2=%d, relid2=%d,%d \n", 
	 rv, d1->GetId(), iphi1,ieta1, d2->GetId(), iphi2,ieta2);   
  
  return rv ; 
}

//____________________________________________________________________________
Int_t AliEMCALClusterizerv1::AreInGroup(AliEMCALDigit * d1, AliEMCALDigit * d2) const
{
  // Tells whether two digits fall within the same supermodule and
  // tower grouping.  The number of towers in a group is controlled by
  // the parameter nTowersInGroup
  //    = 0 are not in same group but continue searching 
  //    = 1 same group
  //    = 2 is in different SM, quit from searching
  //    = 3 different tower group, quit from searching
  //
  // The order of d1 and d2 is important: first (d1) should be a digit 
  // already in a cluster which is compared to a digit (d2)  not yet in a cluster  

  //JLK Question: does the quit from searching assume that the digits
  //are ordered, so that once you are in a different SM, you'll not
  //find another in the list that will match?  How about my TowerGroup search?

  static Int_t rv;
  static Int_t nSupMod1=0, nModule1=0, nIphi1=0, nIeta1=0, iphi1=0, ieta1=0;
  static Int_t nSupMod2=0, nModule2=0, nIphi2=0, nIeta2=0, iphi2=0, ieta2=0;
  static Int_t towerGroup1 = -1, towerGroup2 = -1;
  rv = 0 ;

  fGeom->GetCellIndex(d1->GetId(), nSupMod1,nModule1,nIphi1,nIeta1);
  fGeom->GetCellIndex(d2->GetId(), nSupMod2,nModule2,nIphi2,nIeta2);
  if(nSupMod1 != nSupMod2) return 2; // different SM

  static Int_t nTowerInSM = fGeom->GetNCellsInSupMod()/fGeom->GetNCellsInModule();

  //figure out which tower grouping each digit belongs to
  for(int it = 0; it < nTowerInSM/fNTowerInGroup; it++) {
    if(nModule1 <= nTowerInSM - it*fNTowerInGroup) towerGroup1 = it;
    if(nModule2 <= nTowerInSM - it*fNTowerInGroup) towerGroup2 = it;
  }
  if(towerGroup1 != towerGroup2) return 3; //different Towergroup

  //same SM, same towergroup, we're happy
  if(towerGroup1 == towerGroup2 && towerGroup2 >= 0)
    rv = 1;

  if (gDebug == 2 && rv==1)
  printf("AreInGroup: neighbours=%d, id1=%d, relid1=%d,%d \n id2=%d, relid2=%d,%d \n",
         rv, d1->GetId(), iphi1,ieta1, d2->GetId(), iphi2,ieta2);

  return rv ;
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::Unload() 
{
  // Unloads the Digits and RecPoints
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));
    
  emcalLoader->UnloadDigits() ; 
  emcalLoader->UnloadRecPoints() ; 
}
 
//____________________________________________________________________________
void AliEMCALClusterizerv1::WriteRecPoints()
{

  // Creates new branches with given title
  // fills and writes into TreeR.

  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));

  TObjArray * aECARecPoints = emcalLoader->RecPoints() ; 

  TClonesArray * digits = emcalLoader->Digits() ; 
  TTree * treeR = emcalLoader->TreeR();  
  if ( treeR==0 ) {
    emcalLoader->MakeRecPointsContainer();
    treeR = emcalLoader->TreeR();  
  }
  Int_t index ;

  //Evaluate position, dispersion and other RecPoint properties for EC section
  for(index = 0; index < aECARecPoints->GetEntries(); index++)
    (dynamic_cast<AliEMCALRecPoint *>(aECARecPoints->At(index)))->EvalAll(fECAW0,digits) ;
  
  aECARecPoints->Sort() ;

  for(index = 0; index < aECARecPoints->GetEntries(); index++) {
    (dynamic_cast<AliEMCALRecPoint *>(aECARecPoints->At(index)))->SetIndexInList(index) ;
    (dynamic_cast<AliEMCALRecPoint *>(aECARecPoints->At(index)))->Print();
  }

  Int_t bufferSize = 32000 ;    
  Int_t splitlevel = 0 ; 

  //EC section branch
  
  TBranch * branchECA = 0;
  if ((branchECA = treeR->GetBranch("EMCALECARP")))
    branchECA->SetAddress(&aECARecPoints);
  else
    treeR->Branch("EMCALECARP","TObjArray",&aECARecPoints,bufferSize,splitlevel);
  treeR->Fill() ;

  emcalLoader->WriteRecPoints("OVERWRITE");

}

//____________________________________________________________________________
void AliEMCALClusterizerv1::MakeClusters(char* option)
{
  // Steering method to construct the clusters stored in a list of Reconstructed Points
  // A cluster is defined as a list of neighbour digits

  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));

  TObjArray * aECARecPoints = emcalLoader->RecPoints() ; 
    
  if (fGeom==0) 
     AliFatal("Did not get geometry from EMCALLoader");

  aECARecPoints->Clear();

  //Start with pseudoclusters, if option
  if(strstr(option,"pseudo")) { // ?? Nov 3, 2006 
    TClonesArray * digits  = emcalLoader->Digits() ; 
    TClonesArray * digitsC =  dynamic_cast<TClonesArray*>(digits->Clone());

    TIter nextdigit(digitsC) ;
    AliEMCALDigit * digit;

    AliEMCALRecPoint * recPoint = 0 ; 
    int ineb=0;
    nextdigit.Reset();

    // PseudoClusterization starts  
    while ( (digit = dynamic_cast<AliEMCALDigit *>(nextdigit())) ) { // scan over the list of digitsC
      recPoint = 0 ; 
      TArrayI clusterECAdigitslist(1000); // what is this   

      if(fGeom->CheckAbsCellId(digit->GetId()) ) { //Is this an EMCAL digit? Just maing sure...
	Int_t iDigitInECACluster = 0; // start a new Tower RecPoint
	if(fNumberOfECAClusters >= aECARecPoints->GetSize()) aECARecPoints->Expand(2*fNumberOfECAClusters+1) ;
	AliEMCALRecPoint * rp = new  AliEMCALRecPoint("") ; 
	aECARecPoints->AddAt(rp, fNumberOfECAClusters) ;
	recPoint = dynamic_cast<AliEMCALRecPoint *>(aECARecPoints->At(fNumberOfECAClusters)) ; 
	fNumberOfECAClusters++ ; 

	recPoint->SetClusterType(AliESDCaloCluster::kPseudoCluster);

	recPoint->AddDigit(*digit, digit->GetAmp()) ; 
	clusterECAdigitslist[iDigitInECACluster] = digit->GetIndexInList() ;	
	iDigitInECACluster++ ; 
	digitsC->Remove(digit) ; 
	AliDebug(1,Form("MakePseudoClusters: OK id = %d, adc = %f \n", digit->GetId(), digit->GetAmp()));
	nextdigit.Reset(); // will start from beggining
      
	AliEMCALDigit * digitN = 0; // digi neighbor
	Int_t index = 0 ;

	// Find the neighbours
	while (index < iDigitInECACluster){ // scan over digits already in cluster 
	  digit =  (AliEMCALDigit*)digits->At(clusterECAdigitslist[index]);
	  index++ ; 
	  while ( (digitN = (AliEMCALDigit *)nextdigit())) { // scan over the reduced list of digits 
	    ineb = AreInGroup(digit, digitN);       // call (digit,digitN) in THAT oder !!!!! 
	    switch (ineb ) {
              case 0 :   // not a neighbour
	       break ;
	      case 1 :   // are neighbours 
	       recPoint->AddDigit(*digitN, digitN->GetAmp() ) ;
	       clusterECAdigitslist[iDigitInECACluster] = digitN->GetIndexInList() ; 
	       iDigitInECACluster++ ; 
	       digitsC->Remove(digitN) ;
	       break ;
	      case 2 :   // different SM
	       break ;
              case 3 :   // different Tower group
	       break ;
	    } // switch
	  } // scan over the reduced list of digits
	} // scan over digits already in cluster
	nextdigit.Reset() ;  // will start from beggining
      }
    }
    if(recPoint)
      AliDebug(1,Form("MakeClusters: cl.e %f \n", recPoint->GetEnergy())); 
    //if(recPoint) cout << "cl.e " << recPoint->GetEnergy() << endl; 
    delete digitsC ;
  }

  //Now do real clusters
  TClonesArray * digits  = emcalLoader->Digits() ; 
  TClonesArray * digitsC =  dynamic_cast<TClonesArray*>(digits->Clone()); // will work with thic copy

  TIter nextdigitC(digitsC) ;
  AliEMCALDigit * digit;

  AliEMCALRecPoint * recPoint = 0 ; 
  int ineb=0;
  nextdigitC.Reset();

  double e = 0.0, ehs = 0.0;
  while ( (digit = dynamic_cast<AliEMCALDigit *>(nextdigitC())) ) { // clean up digits
    e = Calibrate(digit->GetAmp(), digit->GetId());
    AliEMCALHistoUtilities::FillH1(fHists, 10, digit->GetAmp());
    AliEMCALHistoUtilities::FillH1(fHists, 11, e);
    if(e < fMinECut ) digitsC->Remove(digit);
    else              ehs += e;
  } 
  AliDebug(1,Form("MakeClusters: Number of digits %d  -> (e %f), ehs %d\n",
		  digits->GetEntries(),fMinECut,ehs));
  //cout << " Number of digits " << digits->GetEntries() << " -> (e>" <<fMinECut <<")";
  //cout << digitsC->GetEntries()<< " ehs "<<ehs<<endl; 

  // Clusterization starts    
  //  cout << "Outer Loop" << endl;
  ineb=0;
  nextdigitC.Reset();
  recPoint = 0 ; 
  while ( (digit = dynamic_cast<AliEMCALDigit *>(nextdigitC())) ) { // scan over the list of digitsC
    recPoint = 0 ; 
    TArrayI clusterECAdigitslist(1000); // what is this   

    if(fGeom->CheckAbsCellId(digit->GetId()) && (Calibrate(digit->GetAmp(), digit->GetId()) > fECAClusteringThreshold  ) ){
      Int_t iDigitInECACluster = 0; // start a new Tower RecPoint
      if(fNumberOfECAClusters >= aECARecPoints->GetSize()) aECARecPoints->Expand(2*fNumberOfECAClusters+1) ;
      AliEMCALRecPoint * rp = new  AliEMCALRecPoint("") ; 
      aECARecPoints->AddAt(rp, fNumberOfECAClusters) ;
      recPoint = dynamic_cast<AliEMCALRecPoint *>(aECARecPoints->At(fNumberOfECAClusters)) ; 
      fNumberOfECAClusters++ ; 

      recPoint->SetClusterType(AliESDCaloCluster::kClusterv1);

      recPoint->AddDigit(*digit, Calibrate(digit->GetAmp(), digit->GetId())) ; 
      clusterECAdigitslist[iDigitInECACluster] = digit->GetIndexInList() ;	
      iDigitInECACluster++ ; 
      digitsC->Remove(digit) ; 
      AliDebug(1,Form("MakeClusters: OK id = %d, ene = %f , thre = %f \n", digit->GetId(),
      Calibrate(digit->GetAmp(),digit->GetId()), fECAClusteringThreshold));  
      nextdigitC.Reset(); // will start from beggining
      
      AliEMCALDigit * digitN = 0; // digi neighbor
      Int_t index = 0 ;

      // Find the neighbours
    NEIGHBOURS:
      while (index < iDigitInECACluster){ // scan over digits already in cluster 
	digit =  (AliEMCALDigit*)digits->At(clusterECAdigitslist[index]);
	index++ ; 
        while ( (digitN = (AliEMCALDigit *)nextdigitC())) { // scan over the reduced list of digits 

	  ineb = AreNeighbours(digit, digitN);       // call (digit,digitN) in THAT oder !!!!! 
          if(ineb==1) {
	    recPoint->AddDigit(*digitN, Calibrate(digitN->GetAmp(),digitN->GetId()) ) ;
	    clusterECAdigitslist[iDigitInECACluster] = digitN->GetIndexInList() ; 
	    iDigitInECACluster++ ; 
	    digitsC->Remove(digitN) ; 
      // Have to start from begining for clusters digits too ; Nov 3,2006
            index = 0;
            nextdigitC.Reset() ; 
            goto NEIGHBOURS;
	  } // if(ineb==1)

        } // scan over the reduced list of digits
        nextdigitC.Reset() ;  // will start from beginning
      } // scan over digits already in cluster
    }
  } // while digit 
  if(recPoint)
    AliDebug(1,Form("MakeClusters: cl.e %f \n", recPoint->GetEnergy())); 
  //if(recPoint) cout << "cl.e " << recPoint->GetEnergy() << endl; 
  delete digitsC ;

  AliDebug(1,Form("total no of clusters %d from %d digits",fNumberOfECAClusters,digits->GetEntriesFast())); 
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::MakeUnfolding() const
{
  Fatal("AliEMCALClusterizerv1::MakeUnfolding", "--> Unfolding not implemented") ;
 
}

//____________________________________________________________________________
Double_t  AliEMCALClusterizerv1::ShowerShape(Double_t r)
{ 
  // Shape of the shower (see EMCAL TDR)
  // If you change this function, change also the gradient evaluation in ChiSquare()

  Double_t r4    = r*r*r*r ;
  Double_t r295  = TMath::Power(r, 2.95) ;
  Double_t shape = TMath::Exp( -r4 * (1. / (2.32 + 0.26 * r4) + 0.0316 / (1 + 0.0652 * r295) ) ) ;
  return shape ;
}

//____________________________________________________________________________
void  AliEMCALClusterizerv1::UnfoldCluster(AliEMCALRecPoint * /*iniTower*/, 
					   Int_t /*nMax*/, 
					   AliEMCALDigit ** /*maxAt*/, 
					   Float_t * /*maxAtEnergy*/) const
{
  // Performs the unfolding of a cluster with nMax overlapping showers 
  
  Fatal("UnfoldCluster", "--> Unfolding not implemented") ;

}

//_____________________________________________________________________________
void AliEMCALClusterizerv1::UnfoldingChiSquare(Int_t & /*nPar*/, Double_t * /*Grad*/,
					       Double_t & /*fret*/,
					       Double_t * /*x*/, Int_t /*iflag*/)
{
  // Calculates the Chi square for the cluster unfolding minimization
  // Number of parameters, Gradient, Chi squared, parameters, what to do
  
  ::Fatal("UnfoldingChiSquare","Unfolding not implemented") ;
}
//____________________________________________________________________________
void AliEMCALClusterizerv1::Print(Option_t * /*option*/)const
{
  // Print clusterizer parameters

  TString message("\n") ; 
  
  if( strcmp(GetName(), "") !=0 ){
    
    // Print parameters
 
    TString taskName(GetName()) ;
    taskName.ReplaceAll(Version(), "") ;
    
    printf("--------------- "); 
    printf(taskName.Data()) ; 
    printf(" "); 
    printf(GetTitle()) ; 
    printf("-----------\n");  
    printf("Clusterizing digits from the file: "); 
    printf(taskName.Data());  
    printf("\n                           Branch: "); 
    printf(GetName()); 
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
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));
  TObjArray * aECARecPoints = emcalLoader->RecPoints() ; 
    
  if(strstr(option,"deb")) {
    printf("PrintRecPoints: Clusterization result:") ; 
  
    printf("event # %d\n", emcalLoader->GetRunLoader()->GetEventNumber() ) ;
    printf("           Found %d ECA Rec Points\n ", 
	 aECARecPoints->GetEntriesFast()) ; 
  }

  fRecPointsInRun +=  aECARecPoints->GetEntriesFast() ; 
  
  if(strstr(option,"all")) {
    if(strstr(option,"deb")) {
      printf("\n-----------------------------------------------------------------------\n") ;
      printf("Clusters in ECAL section\n") ;
      printf("Index    Ene(GeV) Multi Module     GX    GY   GZ  lX    lY   lZ   Dispersion Lambda 1   Lambda 2  # of prim  Primaries list\n") ;
    }
   Int_t index =0;
   Float_t maxE=0; 
   Float_t maxL1=0; 
   Float_t maxL2=0; 
   Float_t maxDis=0; 

    AliEMCALHistoUtilities::FillH1(fHists, 12, double(aECARecPoints->GetEntries()));

    for (index = 0 ; index < aECARecPoints->GetEntries() ; index++) {
      AliEMCALRecPoint * rp = dynamic_cast<AliEMCALRecPoint * >(aECARecPoints->At(index)) ; 
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
  /////////////
      if(rp->GetEnergy()>maxE){
	      maxE=rp->GetEnergy();
	      maxL1=lambda[0];
	      maxL2=lambda[1];
	      maxDis=rp->GetDispersion();
      }
      fPointE->Fill(rp->GetEnergy());
      fPointL1->Fill(lambda[0]);
      fPointL2->Fill(lambda[1]);
      fPointDis->Fill(rp->GetDispersion());
      fPointMult->Fill(rp->GetMultiplicity());
      ///////////// 
      if(strstr(option,"deb")){ 
        for (Int_t iprimary=0; iprimary<nprimaries; iprimary++) {
	  printf("%d ", primaries[iprimary] ) ; 
        }
      }
    }

      fMaxE->Fill(maxE);
      fMaxL1->Fill(maxL1);
      fMaxL2->Fill(maxL2);
      fMaxDis->Fill(maxDis);

    if(strstr(option,"deb"))
    printf("\n-----------------------------------------------------------------------\n");
  }
}
TList* AliEMCALClusterizerv1::BookHists()
{
  //set up histograms for monitoring clusterizer performance

  gROOT->cd();

	fPointE = new TH1F("00_pointE","point energy", 2000, 0.0, 150.);
	fPointL1 = new TH1F("01_pointL1","point L1", 1000, 0.0, 3.);
	fPointL2 = new TH1F("02_pointL2","point L2", 1000, 0.0, 3.);
	fPointDis = new TH1F("03_pointDisp","point dispersion", 1000, 0.0, 10.);
	fPointMult = new TH1F("04_pointMult","#cell in point(cluster)", 101, -0.5, 100.5);
	fDigitAmp = new TH1F("05_digitAmp","Digit Amplitude", 2000, 0.0, 5000.);
	fMaxE = new TH1F("06_maxE","Max point energy", 2000, 0.0, 150.);
	fMaxL1 = new TH1F("07_maxL1","Largest (first) of eigenvalue of covariance matrix", 1000, 0.0, 3.);
	fMaxL2 = new TH1F("08_maxL2","Smalest (second) of eigenvalue of covariace matrix", 1000, 0.0, 3.);
	fMaxDis = new TH1F("09_maxDis","Point dispersion", 1000, 0.0, 10.); // 9
	//
        new TH1F("10_adcOfDigits","adc of digits(threshold control)", 1001, -0.5, 1000.5);   // 10
        new TH1F("11_energyOfDigits","energy of digits(threshold control)", 1000, 0.0, 1.);  // 11
        new TH1F("12_numberOfPoints","number of points(clusters)", 101, -0.5, 100.5);        // 12

  return AliEMCALHistoUtilities::MoveHistsToList("EmcalClusterizerv1ControlHists", kFALSE);
}

void AliEMCALClusterizerv1::SaveHists(const char *fn)
{
  AliEMCALHistoUtilities::SaveListOfHists(fHists, fn, kTRUE);
}

void  AliEMCALClusterizerv1::PrintRecoInfo()
{
  printf(" AliEMCALClusterizerv1::PrintRecoInfo() : version %s \n", Version() );
  TH1F *h = (TH1F*)fHists->At(12);
  if(h) {
    printf(" ## Multiplicity of RecPoints ## \n");
    for(int i=1; i<=h->GetNbinsX(); i++) {
      int nbin = int((*h)[i]);
      int mult = int(h->GetBinCenter(i));
      if(nbin > 0) printf(" %i : %5.5i %6.3f %% \n", mult, nbin, 100.*nbin/h->GetEntries()); 
    }    
  }
}

void AliEMCALClusterizerv1::DrawLambdasHists()
{
  if(fMaxL1) {
    fMaxL1->Draw();
    if(fMaxL2) fMaxL2->Draw("same");
    if(fMaxDis) {
      fMaxDis->Draw("same");
    }
  }
}

void AliEMCALClusterizerv1::Browse(TBrowser* b)
{
  if(fHists) b->Add(fHists);
  if(fGeom)  b->Add(fGeom);
  TTask::Browse(b);
}
