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
//  Base class for the clusterization algorithm (pure abstract)
//*--
//*-- Author: Yves Schutz  SUBATECH 
// 
//   Clusterization mother class. Contains common methods/data members of different 
//   clusterizers. GCB 2010
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TClonesArray.h"
#include "TTree.h"
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
#include <cassert>

// --- AliRoot header files ---
#include "AliEMCALClusterizer.h"
#include "AliEMCALReconstructor.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliLog.h"
#include "AliEMCAL.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALRecParam.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecParam.h"
#include "AliCDBManager.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALCalibData.h"
class AliCDBStorage;
#include "AliCDBEntry.h"

ClassImp(AliEMCALClusterizer)

//____________________________________________________________________________
AliEMCALClusterizer::AliEMCALClusterizer():
  fDigitsArr(NULL),
  fTreeR(NULL),
  fRecPoints(NULL),
  fGeom(NULL),
  fCalibData(NULL), 
  fCaloPed(NULL),
  fADCchannelECA(0.),fADCpedestalECA(0.),
  fTimeMin(-1.),fTimeMax(1.),fTimeCut(1.),
  fDefaultInit(kFALSE),fToUnfold(kFALSE),
  fNumberOfECAClusters(0), fECAClusteringThreshold(0.),
  fECALocMaxCut(0.),fECAW0(0.),fMinECut(0.)
{
  // ctor
  
  Init();
  
}

//____________________________________________________________________________
AliEMCALClusterizer::AliEMCALClusterizer(AliEMCALGeometry* geometry): 
  fDigitsArr(NULL),
  fTreeR(NULL),
  fRecPoints(NULL),
  fGeom(geometry),
  fCalibData(NULL), 
  fCaloPed(NULL),
  fADCchannelECA(0.),fADCpedestalECA(0.),
  fTimeMin(-1.),fTimeMax(1.),fTimeCut(1.),
  fDefaultInit(kFALSE),fToUnfold(kFALSE),
  fNumberOfECAClusters(0), fECAClusteringThreshold(0.),
  fECALocMaxCut(0.),fECAW0(0.),fMinECut(0.)
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
  
}

//____________________________________________________________________________
AliEMCALClusterizer::AliEMCALClusterizer(AliEMCALGeometry* geometry, AliEMCALCalibData * calib, AliCaloCalibPedestal * caloped): 
  fDigitsArr(NULL),
  fTreeR(NULL),
  fRecPoints(NULL),
  fGeom(geometry),
  fCalibData(calib),
  fCaloPed(caloped),
  fADCchannelECA(0.),fADCpedestalECA(0.),
  fTimeMin(-1.),fTimeMax(1.),fTimeCut(1.),
  fDefaultInit(kFALSE),fToUnfold(kFALSE),
  fNumberOfECAClusters(0), fECAClusteringThreshold(0.),
  fECALocMaxCut(0.),fECAW0(0.),fMinECut(0.)
{
	// ctor, geometry and calibration are initialized elsewhere.
	
	if (!fGeom)
		AliFatal("Geometry not initialized.");
  
}


//____________________________________________________________________________
AliEMCALClusterizer::~AliEMCALClusterizer()
{
  // dtor
  //Already deleted in AliEMCALReconstructor.

//   if(fGeom)      delete fGeom;
//   if(fCalibData) delete fCalibData;
//   if(fCaloPed)   delete fCaloPed;

//   if (fDigitsArr) {
//     fDigitsArr->Clear("C");
//     delete fDigitsArr;
//   }
//   if (fRecPoints) {
//     fRecPoints->Delete();
//     delete fRecPoints;
//  }
}

//____________________________________________________________________________
Float_t  AliEMCALClusterizer::Calibrate(const Float_t amp, const Float_t time, const Int_t absId) 
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
void AliEMCALClusterizer::GetCalibrationParameters() 
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
void AliEMCALClusterizer::GetCaloCalibPedestal() 
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
void AliEMCALClusterizer::Init()
{
  // Make all memory allocations which can not be done in default constructor.
  // Attach the Clusterizer task to the list of EMCAL tasks
  
  AliRunLoader *rl = AliRunLoader::Instance();
  if (rl->GetAliRun() && rl->GetAliRun()->GetDetector("EMCAL"))
    fGeom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
  else 
    fGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
  
  AliDebug(1,Form("geom %p",fGeom));
  
  if(!gMinuit) 
    gMinuit = new TMinuit(100) ;
  
}

//____________________________________________________________________________
void AliEMCALClusterizer::InitParameters()
{ 
  // Initializes the parameters for the Clusterizer
  fNumberOfECAClusters = 0 ;
  fCalibData           = 0 ;
  fCaloPed             = 0 ;
	
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
void AliEMCALClusterizer::Print(Option_t * /*option*/)const
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
    printf("AliEMCALClusterizer not initialized ") ;
}

//____________________________________________________________________________
void AliEMCALClusterizer::PrintRecPoints(Option_t * option)
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
    Int_t index; 
    for (index =  0 ; index < fRecPoints->GetEntries() ; index++) {
      AliEMCALRecPoint * rp = dynamic_cast<AliEMCALRecPoint * >(fRecPoints->At(index)) ; 
      TVector3  globalpos;  
      //rp->GetGlobalPosition(globalpos);
      TVector3  localpos;  
      rp->GetLocalPosition(localpos);
      Float_t lambda[2]; 
      rp->GetElipsAxis(lambda);
      
      Int_t nprimaries=0;
      Int_t * primaries = rp->GetPrimaries(nprimaries);
      
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
void  AliEMCALClusterizer::PrintRecoInfo()
{
  //print reco version
  printf(" AliEMCALClusterizer::PrintRecoInfo() : version %s \n", Version() );
  
}

//____________________________________________________________________________
void AliEMCALClusterizer::SetInput(TTree *digitsTree)
{
  // Read the digits from the input tree
  TBranch *branch = digitsTree->GetBranch("EMCAL");
  if (!branch) { 
    AliError("can't get the branch with the EMCAL digits !");
    return;
  }
  if (!fDigitsArr)
    fDigitsArr = new TClonesArray("AliEMCALDigit",100);
  branch->SetAddress(&fDigitsArr);
  branch->GetEntry(0);
}

//____________________________________________________________________________
void AliEMCALClusterizer::SetOutput(TTree *clustersTree)
{
  // Read the digits from the input tree
  fTreeR = clustersTree;
  
  AliDebug(9, "Making array for EMCAL clusters");
  fRecPoints = new TObjArray(100) ;
  Int_t split = 0;
  Int_t bufsize = 32000;
  fTreeR->Branch("EMCALECARP", "TObjArray", &fRecPoints, bufsize, split);
}



