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

// --- ROOT system ---
#include <TTree.h>
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
#include "AliEMCALCalibTime.h"
class AliCDBStorage;
#include "AliCDBEntry.h"

/// \cond CLASSIMP
ClassImp(AliEMCALClusterizer) ;
/// \endcond

///
/// Default Constructor
//____________________________________________________________________________
AliEMCALClusterizer::AliEMCALClusterizer():
  fIsInputCalibrated(kFALSE),
  fJustClusters(kFALSE),
  fDigitsArr(NULL),
  fTreeR(NULL),
  fRecPoints(NULL),
  fGeom(NULL),
  fCalibData(NULL), 
  fCalibTime(NULL), 
  fCaloPed(NULL),
  fADCchannelECA(0.),fADCpedestalECA(0.), fTimeECA(0.),
  fTimeMin(-1.),fTimeMax(1.),fTimeCut(1.), fTimeCalibration(0),
  fDefaultInit(kFALSE),fToUnfold(kFALSE),
  fNumberOfECAClusters(0), fECAClusteringThreshold(0.),
  fECALocMaxCut(0.),fECAW0(0.),fMinECut(0.),fRejectBelowThreshold(0),
  fClusterUnfolding(NULL)
{  
  Init();
}

///
/// Constructor with the indication of the file where header Tree and digits Tree are stored.
/// Use this contructor to avoid usage of Init() which uses runloader.
/// Change needed by HLT - MP.
/// Note for the future: the use on runloader should be avoided or optional at least.
/// Another way is to make Init virtual and protected at least 
/// such that the deriving classes can overload Init();
/// 
/// \param geometry: EMCal geometry pointer
//____________________________________________________________________________
AliEMCALClusterizer::AliEMCALClusterizer(AliEMCALGeometry* geometry): 
  fIsInputCalibrated(kFALSE),
  fJustClusters(kFALSE),
  fDigitsArr(NULL),
  fTreeR(NULL),
  fRecPoints(NULL),
  fGeom(geometry),
  fCalibData(NULL), 
  fCalibTime(NULL), 
  fCaloPed(NULL),
  fADCchannelECA(0.),fADCpedestalECA(0.), fTimeECA(0.),
  fTimeMin(-1.),fTimeMax(1.),fTimeCut(1.), fTimeCalibration(0),
  fDefaultInit(kFALSE),fToUnfold(kFALSE),
  fNumberOfECAClusters(0), fECAClusteringThreshold(0.),
  fECALocMaxCut(0.),fECAW0(0.),fMinECut(0.),fRejectBelowThreshold(0),
  fClusterUnfolding(NULL)
{  
  if (!fGeom)
  {
    AliFatal("Geometry not initialized.");
  }
  
  Int_t i=0;
  for (i = 0; i < 8; i++)
    fSSPars[i] = 0.;
  
  for (i = 0; i < 3; i++)
  {
    fPar5[i] = 0.;
    fPar6[i] = 0.;
  }
}

///
/// Constructor, geometry and calibration are initialized elsewhere.
///
/// \param geometry: EMCal geometry pointer
/// \param calib: EMCal energy calibration container
/// \param calibt: EMCal time calibration container
/// \param caloped: EMCal bad map container
//____________________________________________________________________________
AliEMCALClusterizer::AliEMCALClusterizer(AliEMCALGeometry *geometry, 
                                         AliEMCALCalibData *calib, 
                                         AliEMCALCalibTime *calibt, 
                                         AliCaloCalibPedestal *caloped): 
  fIsInputCalibrated(kFALSE),
  fJustClusters(kFALSE),
  fDigitsArr(NULL),
  fTreeR(NULL),
  fRecPoints(NULL),
  fGeom(geometry),
  fCalibData(calib),
  fCalibTime(calibt),
  fCaloPed(caloped),
  fADCchannelECA(0.),fADCpedestalECA(0.), fTimeECA(0.),
  fTimeMin(-1.),fTimeMax(1.),fTimeCut(1.), fTimeCalibration(0),
  fDefaultInit(kFALSE),fToUnfold(kFALSE),
  fNumberOfECAClusters(0), fECAClusteringThreshold(0.),
  fECALocMaxCut(0.),fECAW0(0.),fMinECut(0.),fRejectBelowThreshold(0),
  fClusterUnfolding(NULL)
{  
  if (!fGeom)
    AliFatal("Geometry not initialized.");
  
  Int_t i=0;
  for (i = 0; i < 8; i++)
    fSSPars[i] = 0.;
  
  for (i = 0; i < 3; i++) 
  {
    fPar5[i] = 0.;
    fPar6[i] = 0.;
  }
}

///
/// Destructtor
///
/// Already deleted in AliEMCALReconstructor.
//____________________________________________________________________________
AliEMCALClusterizer::~AliEMCALClusterizer()
{
  if(fClusterUnfolding) delete fClusterUnfolding;

  // make sure we delete the rec points array
  DeleteRecPoints();

  //Delete digits array
  DeleteDigits();
}

///
/// Free the cluster array
//____________________________________________________________________________
void AliEMCALClusterizer::DeleteRecPoints()
{
  if (fRecPoints) 
    {
      AliDebug(2, "Deleting fRecPoints.");
      fRecPoints->Delete();
      delete fRecPoints;
      fRecPoints = 0;
    }
}

///
/// Free the digits array
//____________________________________________________________________________
void AliEMCALClusterizer::DeleteDigits()
{
  if (fDigitsArr) 
    {
      AliDebug(2, "Deleting fDigitsArr.");
      fDigitsArr->Clear("C");
      delete fDigitsArr;
      fDigitsArr = 0;
    }
}

///
/// Convert digitized amplitude into energy, calibrate time
/// Calibration parameters are taken from OCDB : OCDB/EMCAL/Calib/Data
/// Return energy with default parameters if calibration is not available
/// 
/// \param amp: cell energy input and calibrated output
/// \param time: cell time input and calibrated output
/// \param absId: cell absolute ID number
//____________________________________________________________________________
void AliEMCALClusterizer::Calibrate(Float_t & amp, Float_t & time, const Int_t absId) 
{
  if (!fCalibData && !fCaloPed) 
  {
    if (fIsInputCalibrated == kTRUE)
    {
      AliDebug(10, Form("Input already calibrated!"));
      return ;
    }
    else{
      AliFatal("OCDB calibration and bad map parameters are not available");
      return;
    }   
  }
  
  if (fGeom==0)
    AliFatal("Did not get geometry from EMCALLoader") ;
  
  Int_t iSupMod = -1;
  Int_t nModule = -1;
  Int_t nIphi   = -1;
  Int_t nIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  
  Bool_t bCell = fGeom->GetCellIndex(absId, iSupMod, nModule, nIphi, nIeta) ;
  if(!bCell)
  {
    fGeom->PrintGeometry();
    AliError(Form("Wrong cell id number : %i", absId));
    //assert(0); // GCB: This aborts reconstruction of raw simulations 
    //where simulation had more SM than default geometry, 
    //change to return 0, to avoid aborting good generations.
    amp  = 0;
    time = 0;
    return ;
  }
  
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,nModule,nIphi, nIeta,iphi,ieta);
  
  // Check if channel is bad (dead or hot), in this case return 0.	
  // Gustavo: 15-12-09 In case of RAW data this selection is already done, but not in simulation.
  // for the moment keep it here but remember to do the selection at the sdigitizer level 
  // and remove it from here
  if (fCaloPed) 
  {
    Int_t channelStatus = (Int_t)(fCaloPed->GetDeadMap(iSupMod))->GetBinContent(ieta,iphi);
    if(channelStatus == AliCaloCalibPedestal::kHot || channelStatus == AliCaloCalibPedestal::kDead) {
      AliDebug(2,Form("Tower from SM %d, ieta %d, iphi %d is BAD : status %d !!!",iSupMod,ieta,iphi, channelStatus));
      amp  = 0 ;
      time = 0 ;
      return ;
    }
  }
  
  if (fIsInputCalibrated || !fCalibData)
  {
    AliDebug(10, Form("Input already calibrated!"));
    return ;
  }
  
  // Energy calibration
  fADCchannelECA  = fCalibData->GetADCchannel (iSupMod,ieta,iphi);
  fADCpedestalECA = fCalibData->GetADCpedestal(iSupMod,ieta,iphi);
  
  amp   = amp * fADCchannelECA - fADCpedestalECA ;  
  
  // Time calibration
  if(fCalibTime && fTimeCalibration)
  {
    Int_t bc = 0; // Get somehow the bunch crossing number
    
    fTimeECA        = fCalibTime->GetTimeChannel(iSupMod,ieta,iphi, bc) * 1e-9; // ns to s
    
    time -= fTimeECA ;
  }
}

///
/// Set calibration parameters:
/// If calibration database exists, they are read from database,
/// otherwise, they are taken from digitizer.
/// It is a user responsilibity to open CDB before reconstruction, 
/// for example: 
/// AliCDBStorage* storage = AliCDBManager::Instance()->GetStorage("local://CalibDB");
//____________________________________________________________________________
void AliEMCALClusterizer::GetCalibrationParameters() 
{
  if (fIsInputCalibrated)
    return;
  
  // Check if energy calibration is stored in data base
  if(!fCalibData)
  {
    AliCDBEntry *entry = (AliCDBEntry*) 
    AliCDBManager::Instance()->Get("EMCAL/Calib/Data");
    if (entry) fCalibData =  (AliEMCALCalibData*) entry->GetObject();
  }
  
  if(!fCalibData)
    AliFatal("Energy calibration parameters not found in CDB!");
  
  if(fTimeCalibration)
  {
    // Check if energy calibration is stored in data base
    if(!fCalibTime)
    {
      AliCDBEntry *entry = (AliCDBEntry*) 
      AliCDBManager::Instance()->Get("EMCAL/Calib/Time");
      if (entry) fCalibTime =  (AliEMCALCalibTime*) entry->GetObject();
    }
    
    if(!fCalibTime)
      AliFatal("Time calibration parameters not found in CDB!");
  }
}

///
/// Set calibration parameters:
/// if calibration database exists, they are read from database,
/// otherwise, they are taken from digitizer.
///
/// It is a user responsilibity to open CDB before reconstruction, 
/// for example: 
/// AliCDBStorage* storage = AliCDBManager::Instance()->GetStorage("local://CalibDB");
//____________________________________________________________________________
void AliEMCALClusterizer::GetCaloCalibPedestal() 
{
  if (fIsInputCalibrated)
    return;
  
  // Check if calibration is stored in data base
  if(!fCaloPed)
  {
    AliCDBEntry *entry = (AliCDBEntry*) 
    AliCDBManager::Instance()->Get("EMCAL/Calib/Pedestals");
    if (entry) fCaloPed =  (AliCaloCalibPedestal*) entry->GetObject();
  }
  
  if(!fCaloPed)
    AliFatal("Pedestal info not found in CDB!");
}

///
/// Make all memory allocations which can not be done in default constructor.
/// Attach the Clusterizer task to the list of EMCAL tasks
//____________________________________________________________________________
void AliEMCALClusterizer::Init()
{  
  AliRunLoader *rl = AliRunLoader::Instance();
  if (rl->GetAliRun())
  {
    AliEMCAL* emcal = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"));
    if(emcal)fGeom = emcal->GetGeometry();
  }
  
  if(!fGeom)
    fGeom = AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
  
  
  AliDebug(1,Form("geom %p",fGeom));
  
  if(!gMinuit) 
    gMinuit = new TMinuit(100) ;
  
  Int_t i=0;
  for (i = 0; i < 8; i++)
    fSSPars[i] = 0.;
  
  for (i = 0; i < 3; i++) 
  {
    fPar5[i] = 0.;
    fPar6[i] = 0.;
  }
}

///
/// Initializes the parameters for the Clusterizer from AliEMCALReconstructor::GetRecParam().
//____________________________________________________________________________
void AliEMCALClusterizer::InitParameters()
{ 
  return InitParameters(AliEMCALReconstructor::GetRecParam());
}

///
/// Initializes the parameters for the Clusterizer
///
/// \param recParam: Object with reco params
//____________________________________________________________________________
void AliEMCALClusterizer::InitParameters(const AliEMCALRecParam* recParam)
{ 
  fNumberOfECAClusters = 0 ;
  fCalibData           = 0 ;
  fCalibTime           = 0 ;
  fCaloPed             = 0 ;
  
  if(!recParam) 
    AliFatal("Reconstruction parameters for EMCAL not set!");
  
  
  fECAClusteringThreshold = recParam->GetClusteringThreshold();
  fECAW0                  = recParam->GetW0();
  fMinECut                = recParam->GetMinECut();    
  fRejectBelowThreshold   = recParam->GetRejectBelowThreshold();
  fToUnfold               = recParam->GetUnfold();
  fECALocMaxCut           = recParam->GetLocMaxCut();
  fTimeCut                = recParam->GetTimeCut();
  fTimeMin                = recParam->GetTimeMin();
  fTimeMax                = recParam->GetTimeMax();
  fTimeCalibration        = recParam->IsTimeCalibrationOn();
  
  //For NxN
  SetNRowDiff(recParam->GetNRowDiff());
  SetNColDiff(recParam->GetNColDiff());
  
  AliDebug(1,Form("Reconstruction parameters: fECAClusteringThreshold=%.3f GeV, fECAW=%.3f, fMinECut=%.3f GeV, "
                  "fToUnfold=%d, fECALocMaxCut=%.3f GeV, "
                  "fTimeCut=%e s,fTimeMin=%e s,fTimeMax=%e s,"
                  "fTimeCalibration %d, fRejectBelowThreshold=%d",
                  fECAClusteringThreshold,fECAW0,fMinECut,fToUnfold,fECALocMaxCut,
                  fTimeCut, fTimeMin, fTimeMax, fTimeCalibration,fRejectBelowThreshold));
  
  if (fToUnfold) 
  {
    Int_t i=0;
    for (i = 0; i < 8; i++) 
    {
      fSSPars[i] = recParam->GetSSPars(i);
    } //end of loop over parameters
    
    for (i = 0; i < 3; i++) 
    {
      fPar5[i] = recParam->GetPar5(i);
      fPar6[i] = recParam->GetPar6(i);
    } //end of loop over parameters
    
    InitClusterUnfolding();
    
    for (i = 0; i < 8; i++) 
    {
      AliDebug(1,Form("unfolding shower shape parameters: fSSPars=%f \n",fSSPars[i]));
    }
    
    for (i = 0; i < 3; i++) 
    {
      AliDebug(1,Form("unfolding parameter 5: fPar5=%f \n",fPar5[i]));
      AliDebug(1,Form("unfolding parameter 6: fPar6=%f \n",fPar6[i]));
    }
  } // to unfold
}

///
/// Print clusterizer parameters
//____________________________________________________________________________
void AliEMCALClusterizer::Print(Option_t * /*option*/) const
{  
  TString message("\n") ; 
  
  if (strcmp(GetName(),"") == 0)
  {
    printf("AliEMCALClusterizer not initialized\n");
    return;
  }
  
  // Print parameters
  TString taskName(Version()) ;
  
  printf("--------------- "); 
  printf("%s",taskName.Data()) ; 
  printf(" "); 
  printf("Clusterizing digits: "); 
  printf("\n                       ECA Local Maximum cut    = %f", fECALocMaxCut); 
  printf("\n                       ECA Logarithmic weight   = %f", fECAW0); 
  if (fToUnfold)
  {
    printf("\nUnfolding on\n");
    printf("Unfolding parameters: fSSpars: \n");
    Int_t i=0;
    for (i = 0; i < 8; i++) 
    {
      printf("fSSPars[%d] = %f \n", i, fSSPars[i]);
    }
    
    printf("Unfolding parameter 5 and 6: fPar5 and fPar6: \n");
    
    for (i = 0; i < 3; i++) 
    {
      printf("fPar5[%d] = %f \n", i, fPar5[i]);
      printf("fPar6[%d] = %f \n", i, fPar6[i]);
    }
  }
  else
    printf("\nUnfolding off\n");
  
  printf("------------------------------------------------------------------"); 
}

///
/// Prints list of RecPoints produced at the current pass of AliEMCALClusterizer
///
/// \param option: flag the parameters to print
//____________________________________________________________________________
void AliEMCALClusterizer::PrintRecPoints(Option_t * option)
{
  if (strstr(option,"deb"))
  {
    printf("PrintRecPoints: Clusterization result:") ; 
    printf("           Found %d ECA Rec Points\n ", 
           fRecPoints->GetEntriesFast()) ; 
  }
  
  if (strstr(option,"all")) 
  {
    if (strstr(option,"deb")) 
    {
      printf("\n-----------------------------------------------------------------------\n") ;
      printf("Clusters in ECAL section\n") ;
      printf("Index    Ene(GeV) Multi Module     GX    GY   GZ  lX    lY   lZ   Dispersion Lambda 1   Lambda 2  # of prim  Primaries list\n") ;
    }
    
    Int_t index; 
    for (index =  0 ; index < fRecPoints->GetEntries() ; index++) 
    {
      AliEMCALRecPoint * rp = dynamic_cast<AliEMCALRecPoint * >(fRecPoints->At(index)) ; 
      if (!rp) 
        continue;
      
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
      
      if(strstr(option,"deb"))
      { 
        for (Int_t iprimary=0; iprimary<nprimaries; iprimary++) 
        {
          printf("%d ", primaries[iprimary] ) ; 
        }
      }
    }
    
    if(strstr(option,"deb"))
      printf("\n-----------------------------------------------------------------------\n");
  }
}

///
/// Print reco version
//___________________________________________________________________
void  AliEMCALClusterizer::PrintRecoInfo()
{
  printf(" AliEMCALClusterizer::PrintRecoInfo() : version %s \n", Version() );
}

///
/// Read the digits from the input tree
//____________________________________________________________________________
void AliEMCALClusterizer::SetInput(TTree *digitsTree)
{
  TBranch *branch = digitsTree->GetBranch("EMCAL");
  if (!branch)
  {
    AliError("can't get the branch with the EMCAL digits !");
    return;
  }
  
  if (!fDigitsArr) fDigitsArr = new TClonesArray("AliEMCALDigit",100);
  else             fDigitsArr->Clear("C"); // avoid leak
  
  branch->SetAddress(&fDigitsArr);
  branch->GetEntry(0);
}

///
/// Set the tree with output clusters
//____________________________________________________________________________
void AliEMCALClusterizer::SetOutput(TTree *clustersTree)
{
  AliDebug(9, "Making array for EMCAL clusters");
  
  fRecPoints = new TObjArray(1000);
  
  if (clustersTree) 
  {
    fTreeR = clustersTree;
    Int_t split = 0;
    Int_t bufsize = 32000;
    fTreeR->Branch("EMCALECARP", "TObjArray", &fRecPoints, bufsize, split);
  }
}

///
/// Flag to indicate that input is calibrated - the case when we run already on ESD
//___________________________________________________________________
void AliEMCALClusterizer::SetInputCalibrated(Bool_t val)
{
  fIsInputCalibrated = val;
}

///
/// Flag to indicate that we are running on ESDs, when calling 
/// rp->EvalAll(fECAW0,fDigitsArr,fJustClusters); in derived classes
//___________________________________________________________________
void AliEMCALClusterizer::SetJustClusters(Bool_t val)
{
  fJustClusters = val;
}
