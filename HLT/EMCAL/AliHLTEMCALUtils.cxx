#include "AliHLTEMCALUtils.h"

#include "TClass.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TUnixSystem.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TArrayS.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliRawReaderMemory.h"

#include "AliEMCALRecParam.h"
#include "AliEMCALReconstructor.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALClusterizerv1.h" 
#include "AliEMCALRawUtils.h"      

#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALDigit.h"
#include "AliESDEvent.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALPID.h"

#include <fstream>
#include <iostream>
using namespace std;

ClassImp(AliHLTEMCALUtils);

Int_t                  AliHLTEMCALUtils::fgDebug        = 0;
AliEMCALGeometry*      AliHLTEMCALUtils::fgGeom         = NULL;
AliEMCALClusterizerv1* AliHLTEMCALUtils::fgClusterizer  = NULL;
AliEMCALRecParam*      AliHLTEMCALUtils::fgRecParam     = NULL;
AliEMCALRawUtils*      AliHLTEMCALUtils::fgRawUtils     = NULL;
TFile*                 AliHLTEMCALUtils::fgGeometryFile = NULL;
TGeoManager*           AliHLTEMCALUtils::fgGeoManager   = NULL;

TTree*                 AliHLTEMCALUtils::fgClustersTreeTmp = NULL;
TTree*                 AliHLTEMCALUtils::fgDigitsTreeTmp   = NULL;  
TClonesArray*          AliHLTEMCALUtils::fgDigitsArrTmp    = NULL;  
TBranch*               AliHLTEMCALUtils::fgDigitsBranchTmp = NULL;  


//____________________________________________________________________________
AliHLTEMCALUtils::AliHLTEMCALUtils()
  : TObject()
{
  //
  // Default constructor
  //
  ;
}

//____________________________________________________________________________
AliHLTEMCALUtils::~AliHLTEMCALUtils()
{
  //
  // Destructor
  //
  ;
}

//____________________________________________________________________________
void AliHLTEMCALUtils::Cleanup()
{
  DeleteStaticMembers();
  DeleteReconstructionTrees();
}

//____________________________________________________________________________
void AliHLTEMCALUtils::DeleteStaticMembers()
{
  if (fgRecParam != AliEMCALReconstructor::GetRecParam())
    {
      delete fgRecParam;
      fgRecParam = NULL;
    }

  delete fgRawUtils;
  fgRawUtils = NULL;

  DeleteGeoManager();

  delete fgGeom;
  fgGeom = NULL;
      
  delete fgClusterizer;
  fgClusterizer = NULL;
}

//____________________________________________________________________________
AliHLTEMCALUtils::AliHLTEMCALUtils(const AliHLTEMCALUtils & /*t*/)  
  : TObject()
{
  //
  // copy ctor not to be used
  //
  AliFatal("May not use.");  
}
  
//____________________________________________________________________________
AliHLTEMCALUtils& AliHLTEMCALUtils::operator = (const AliHLTEMCALUtils & /*t*/)  
{
  //
  // assignement operator not to be used
  //
  AliFatal("May not use.") ;
  return *this ; 
}

//____________________________________________________________________________
void AliHLTEMCALUtils::InitRecParam()
{
  //
  // Please check the AliEMCALReconstructor for comparison
  // Check if the instance of AliEMCALRecParam exists, 
  // if not, get it from OCDB if available, otherwise create a default one
  //

  fgRecParam = (AliEMCALRecParam*) AliEMCALReconstructor::GetRecParam();
  if (fgRecParam)
    return;
  
 if (!fgRecParam  && (AliCDBManager::Instance()->IsDefaultStorageSet())) 
   {
     AliCDBEntry *entry = (AliCDBEntry*) 
       AliCDBManager::Instance()->Get("EMCAL/Config/RecParam");
     if (entry) fgRecParam = (AliEMCALRecParam*) entry->GetObject();
   }
 
  if(!fgRecParam)
    {
      AliWarningClass("The Reconstruction parameters initialized to default.");
      fgRecParam = new AliEMCALRecParam;
    }

  if (!fgRecParam)
    {
      AliErrorClass("Unable to init the reco params. Something is really wrong. Memory?");
    }
  
  AliEMCALReconstructor::SetRecParam(fgRecParam);
}

//____________________________________________________________________________
AliEMCALRecParam* AliHLTEMCALUtils::GetRecParam()
{
  //
  // Init the parameters and reuse the Reconstructor
  //

  AliHLTEMCALUtils::InitRecParam();
  return fgRecParam;
}

//____________________________________________________________________________
AliEMCALRawUtils*      AliHLTEMCALUtils::GetRawUtils(AliEMCALGeometry *pGeometry)
{
  //
  // Init EMCAL raw utils
  // Use the geometry or try to figure out the default
  // The OCDB must be initialized. 
  // Otherwise the default contructor of the raw utils will complain.
  //

  if (fgRawUtils == NULL)
    {
      if (AliCDBManager::Instance()->IsDefaultStorageSet())
	{
	  if (pGeometry)
	    {
	      fgRawUtils = new AliEMCALRawUtils(pGeometry);
	    }
	  else
	    {
	      fgRawUtils = new AliEMCALRawUtils(AliHLTEMCALUtils::GetGeometry());
	    }

	  AliInfoClass("Raw Utils initialized.");
	}
      else
	{
	  AliErrorClass("OCDB not initialized. Unable to init raw utils.");
	}
    }

  return fgRawUtils;
}

//____________________________________________________________________________
AliEMCALClusterizerv1* AliHLTEMCALUtils::GetClusterizer()
{
  //
  // Init EMCAL clusterizer
  //

  if (fgClusterizer == NULL)
    {
      AliEMCALGeometry *geom = GetGeometry();
      fgClusterizer = new AliEMCALClusterizerv1(geom);
      AliInfoClass("ClusterizerV1 initialized.");
   }

  return fgClusterizer;
}

//____________________________________________________________________________
AliEMCALGeometry*      AliHLTEMCALUtils::GetGeometry()
{
  //
  // Init EMCAL geometry
  //

  if (fgGeom == NULL)
    {
      AliInfoClass(Form("Using default geometry"));
      fgGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
    }
  return fgGeom;
}

//____________________________________________________________________________
void AliHLTEMCALUtils::DeleteGeoManager()
{
  //
  // Delete the geom manager and the geom TFile
  //

  if (fgGeometryFile)
    {
      fgGeometryFile->Close();
      delete fgGeometryFile;
      fgGeometryFile = NULL;      
    }
  
  if (fgGeoManager)
    {
      delete fgGeoManager;
      fgGeoManager = NULL;
    }
}

//____________________________________________________________________________
Bool_t AliHLTEMCALUtils::LoadGeoManagerFromFile(const char *fname)
{
  //
  // Open the geom file and get the geom manager
  //

  Bool_t bReturn = kFALSE;

  DeleteGeoManager();

  fgGeometryFile = TFile::Open(fname);
  if (fgGeometryFile)
    {
      fgGeoManager = (TGeoManager *)fgGeometryFile->Get("Geometry");
    }

  if (!fgGeoManager)
    {
      DeleteGeoManager();
      AliErrorClass(Form("Unable to retrieve TGeoManager <Geometry> from %s", fname));
    }

  bReturn = kTRUE;

  return bReturn;
}

//____________________________________________________________________________
Bool_t AliHLTEMCALUtils::InitFakeCDB(const char *cdbpath, Int_t runid)
{
  //
  // Init fake OCDB - use with care.
  // For dummy usage only!
  // Make sure there is no conflict with the existing ocdb instance.
  //

  Bool_t bReturn = kFALSE;

  TString sPath = cdbpath;

  if (sPath.Length() <= 0)
    {
      sPath  = "local://";
      sPath += gSystem->Getenv("ALICE_ROOT");
      AliInfoClass(Form("Setting cdbpath to %s", sPath.Data()));
    }

  AliCDBManager* pCDB = AliCDBManager::Instance();
  if (!pCDB)
    {
      AliErrorClass("Unable to get the CDBManager instance.");      
    }
  else
    {
      pCDB->SetRun(runid); // THIS HAS TO BE RETRIEVED !!!
      pCDB->SetDefaultStorage(sPath.Data());
      AliInfoClass(Form("Fake CDB initialized with path %s", sPath.Data()));
      bReturn = kTRUE;
    }
  
  return bReturn;
}

//____________________________________________________________________________
Bool_t AliHLTEMCALUtils::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree, Option_t* sOption)
{
  // This method is a reusage the corresponding one of AliEMCALReconstructor
  // Conversion from raw data to EMCAL digits.
  // Works on a single-event basis

  // Make sure raw utils are initialized (they depend on geometry and OCDB) !

  // We Assume the rawReader, digitsTree and sOption are valid!

  // Here we assume the raw reader is setup properly, so we do not need extra reset!
  // rawReader->Reset() ; 

  Bool_t bReturn = kFALSE;
  if (!fgRawUtils)
    {
      // Try to initialize - this is potentially dangerous
      // We will use default geometry - most probably
      // The OCDB must exist! otherwise we will fail.

      if (GetRawUtils() == NULL)
	{
	  // This is no go.
	  AliErrorClass("Failed. No RawUtils.");
	  return bReturn;
	}
    }

  if (!fgRecParam)
    {
      // Again we try to init the RecParams to default
      // One should better do it properly before!
      // The ocdb has to be there! for the recparams init
      
      if (GetRecParam() == NULL)
	{
	  // This is no go.
	  AliErrorClass("Failed. No RecParams.");
	  return bReturn;
	}      
    }

  TClonesArray *digitsArr = new TClonesArray("AliEMCALDigit",200);
  Int_t bufsize = 32000;
  digitsTree->Branch("EMCAL", &digitsArr, bufsize);

  fgRawUtils->SetOption(sOption);

  fgRawUtils->SetRawFormatHighLowGainFactor(fgRecParam->GetHighLowGainFactor());
  fgRawUtils->SetRawFormatOrder(fgRecParam->GetOrderParameter());
  fgRawUtils->SetRawFormatTau(fgRecParam->GetTau());
  fgRawUtils->SetNoiseThreshold(fgRecParam->GetNoiseThreshold());
  fgRawUtils->SetNPedSamples(fgRecParam->GetNPedSamples());

  fgRawUtils->Raw2Digits(rawReader,digitsArr);

  digitsTree->Fill();
  digitsArr->Delete();
  delete digitsArr;

  bReturn = kTRUE;
  return bReturn;
}

//____________________________________________________________________________
TTree *AliHLTEMCALUtils::GetDigitsTree()
{
  // 
  // Create the digits tree and the digits arrays
  // 

  if (fgDigitsTreeTmp == NULL)
    {
      fgDigitsTreeTmp = new TTree("EMCALdigits", "EMCAL digits default");      
    }

  if (fgDigitsArrTmp == NULL)
    {
      fgDigitsArrTmp = new TClonesArray("AliEMCALDigit",200);
    }
  
  if (fgDigitsBranchTmp == NULL)
    {
      Int_t bufsize = 32000;
      fgDigitsBranchTmp = fgDigitsTreeTmp->Branch("EMCAL", &fgDigitsArrTmp, bufsize);  
    }

  fgDigitsBranchTmp->SetAddress(&fgDigitsArrTmp);

  return fgDigitsTreeTmp;
}

//____________________________________________________________________________
TTree *AliHLTEMCALUtils::GetClustersTree()
{
  // 
  // Create the clusters tree if necessary
  // 
  
  if (fgClustersTreeTmp == NULL)
    {
      fgClustersTreeTmp = new TTree("EMCALclusters", "EMCAL clusters");
    }
  
  return fgClustersTreeTmp;
}

//____________________________________________________________________________
void AliHLTEMCALUtils::ResetReconstructionTrees()
{
  // 
  // Reset the content of the reconstruction trees (digits and clusters)
  //

  TTree *tmp = 0;
  tmp = GetDigitsTree();
  if (tmp)
    {
      fgDigitsArrTmp->Delete();
      tmp->Reset();
    }

  tmp = GetClustersTree();
  if (tmp)
    {
      tmp->Reset();
    }
}

//____________________________________________________________________________
void AliHLTEMCALUtils::DeleteReconstructionTrees()
{
  // 
  // Delete the reconstruction trees (digits and clusters)
  //

  TTree *tmp = 0;
  tmp = GetDigitsTree();
  if (tmp)
    {
      fgDigitsArrTmp->Delete();
      delete fgClustersTreeTmp;
      fgClustersTreeTmp = NULL;
    }

  tmp = GetClustersTree();
  if (tmp)
    {
      delete fgDigitsTreeTmp;
      fgDigitsTreeTmp = 0;
    }
}

//____________________________________________________________________________
Bool_t AliHLTEMCALUtils::Raw2Clusters(AliRawReader* rawReader, TTree* clustersTree, Option_t* sDigitsOption)
{
  //
  // Here we go from raw to clusters (not realeasing digits on the way)
  // This method combines ConvertDigits and Reconstruct of AliEMCALReconstructor
  //
  
  // Here we assume the raw reader is setup properly, so we do not need extra reset!
  // rawReader->Reset() ; 

  Bool_t bReturn = kFALSE;
  if (!fgRawUtils)
    {
      // Try to initialize - this is potentially dangerous
      // We will use default geometry - most probably
      // The OCDB must exist! otherwise we will fail.

      if (GetRawUtils() == NULL)
	{
	  // This is no go.
	  AliErrorClass("Failed. No RawUtils.");
	  return bReturn;
	}
    }

  if (!fgRecParam)
    {
      // Again we try to init the RecParams to default
      // One should better do it properly before!
      // The ocdb has to be there! for the recparams init
      
      if (GetRecParam() == NULL)
	{
	  // This is no go.
	  AliErrorClass("Failed. No RecParams.");
	  return bReturn;
	}      
    }

  if (GetDigitsTree() == NULL)
    {
      // This is no go.
      AliErrorClass("Failed. No Digits Tree.");
      return bReturn;	  
    }
  
  if (GetClusterizer() == NULL)
    {
      // This is no go.
      AliErrorClass("Failed. No clusterizer.");
      return bReturn;	  
    }

  // raw->digits
  fgRawUtils->SetOption(sDigitsOption);
  fgRawUtils->SetRawFormatHighLowGainFactor(fgRecParam->GetHighLowGainFactor());
  fgRawUtils->SetRawFormatOrder(fgRecParam->GetOrderParameter());
  fgRawUtils->SetRawFormatTau(fgRecParam->GetTau());
  fgRawUtils->SetNoiseThreshold(fgRecParam->GetNoiseThreshold());
  fgRawUtils->SetNPedSamples(fgRecParam->GetNPedSamples());

  fgDigitsArrTmp->Delete();
  fgRawUtils->Raw2Digits(rawReader, fgDigitsArrTmp);
  fgDigitsTreeTmp->Fill();

  // digits->clusters
  fgClusterizer->SetOutput(clustersTree);

  if (fgDigitsArrTmp->GetEntries() > 0) 
    {
      fgClusterizer->SetInput(fgDigitsTreeTmp);      

      if (fgDebug > 0)
	{
	  fgClusterizer->Digits2Clusters("deb all") ;
	}
      else
	{
	  fgClusterizer->Digits2Clusters("");
	}

      fgClusterizer->Clear();      
    }

  bReturn = kTRUE;
  return bReturn;  
}

//____________________________________________________________________________
Bool_t AliHLTEMCALUtils::RawBuffer2Clusters(UChar_t *buffer, ULong_t buffSize, 
					    Int_t eqID, 
					    TTree* clustersTree, Option_t* sDigitsOption)
{
  //
  // Method provided for convenience
  // Take the raw (DDL) mem buffer and put the clusters into the tree
  //

  Bool_t bReturn = kFALSE;

  AliRawReaderMemory rawReader;
  rawReader.Reset();
  rawReader.SetEquipmentID(eqID);
  rawReader.SetMemory(buffer, buffSize);

  bReturn = AliHLTEMCALUtils::Raw2Clusters(&rawReader, clustersTree, sDigitsOption);

  return bReturn;  
}

//____________________________________________________________________________
TTree* AliHLTEMCALUtils::RawBuffer2Clusters(UChar_t *buffer, ULong_t buffSize, 
					    Int_t eqID, 
					    Option_t* sDigitsOption)
{
  //
  // Method provided for convenience
  // Take the raw (DDL) mem buffer and put the clusters into the tree
  //

  Bool_t bReturn = kFALSE;

  AliRawReaderMemory rawReader;
  rawReader.Reset();
  rawReader.SetMemory(buffer, buffSize);
  rawReader.SetEquipmentID(eqID);

  TTree *clustersTree = GetClustersTree();

  bReturn = AliHLTEMCALUtils::Raw2Clusters(&rawReader, clustersTree, sDigitsOption);
  if (bReturn == kFALSE)
    {
      clustersTree->Delete();
      delete clustersTree;
      clustersTree = 0;
    }

  return clustersTree;  
}

//____________________________________________________________________________
UChar_t *AliHLTEMCALUtils::FileToMemory(const char *fname, UChar_t *outpbuffer, ULong_t &buffSize)
{
  //
  // Method provided for convenience
  // Take the raw (DDL) file and copy it to a mem buffer
  //

  buffSize = 0;
  outpbuffer = NULL;
  
  ifstream inputFile(fname, ifstream::in);
  if (!inputFile.good())
    {
      AliErrorClass(Form("Unable to open file", fname));
      return outpbuffer;
    }

  inputFile.seekg (0, ios::end);
  buffSize = inputFile.tellg();
  inputFile.seekg (0, ios::beg);

  void *buff = calloc(buffSize, sizeof(char));  

  if (!buff)
    {
      AliErrorClass(Form("Unable to allocate memory: %d bytes", buffSize)); 
      return outpbuffer;      
    }

  inputFile.read((char*)buff, buffSize);
  outpbuffer = (UChar_t*)buff;

  AliInfoClass(Form("File read: %s. Allocated %d bytes at 0x%x", fname, buffSize, outpbuffer)); 
  return outpbuffer;
}

//____________________________________________________________________________
Bool_t AliHLTEMCALUtils::FillESD(TTree* digitsTree, TTree* clustersTree, AliESDEvent* esd)
{
  //
  // Fill the ESD just like the AliEMCALReconstructor
  //
  // This is offline code taken directly from the AliEMCALReconstructor
  // Factorization of the function in the AliEMCALReconstructor should allow calling the code directly

  // Trigger part on raw data has been removed.

  //########################################
  //##############Fill CaloCells############
  //########################################

  if (digitsTree == 0)
    {
      AliErrorClass("Digits tree is NULL!");
      return kFALSE;
    }

  if (clustersTree == 0)
    {
      AliErrorClass("Clusters tree is NULL!");
      return kFALSE;
    }

  if (esd == 0)
    {
      AliErrorClass("Pointer to ESD not set!");
      return kFALSE;
    }

  TClonesArray *digits = new TClonesArray("AliEMCALDigit",1000);
  TBranch *branchdig = digitsTree->GetBranch("EMCAL");
  if (!branchdig) 
    { 
      AliErrorClass("Can't get the branch with the PHOS digits !");
      return kFALSE;
    }

  branchdig->SetAddress(&digits);
  digitsTree->GetEvent(0);
  Int_t nDigits = digits->GetEntries(), idignew = 0 ;

  AliESDCaloCells &emcCells = *(esd->GetEMCALCells());
  emcCells.CreateContainer(nDigits);
  emcCells.SetType(AliESDCaloCells::kEMCALCell);
  for (Int_t idig = 0 ; idig < nDigits ; idig++) 
    {
      const AliEMCALDigit * dig = (const AliEMCALDigit*)digits->At(idig);
      if(dig->GetAmp() > 0 )
	{
	  emcCells.SetCell(idignew,dig->GetId(),dig->GetAmp(), dig->GetTime());   
	  idignew++;
	}
    }
  emcCells.SetNumberOfCells(idignew);
  emcCells.Sort();
  
  //------------------------------------------------------------
  //-----------------CLUSTERS-----------------------------
  //------------------------------------------------------------

  TObjArray *clusters = new TObjArray(100);
  TBranch *branch = clustersTree->GetBranch("EMCALECARP");
  branch->SetAddress(&clusters);
  clustersTree->GetEvent(0);

  Int_t nClusters = clusters->GetEntries(),  nClustersNew=0;

  esd->SetFirstEMCALCluster(esd->GetNumberOfCaloClusters()); // Put after Phos clusters 

  //######################################################
  //#######################TRACK MATCHING###############
  //######################################################
  //Fill list of integers, each one is index of track to which the cluster belongs.

  // step 1 - initialize array of matched track indexes
  Int_t *matchedTrack = new Int_t[nClusters];
  for (Int_t iclus = 0; iclus < nClusters; iclus++)
    matchedTrack[iclus] = -1;  // neg. index --> no matched track
  
  // step 2, change the flag for all matched clusters found in tracks
  Int_t iemcalMatch = -1;
  Int_t endtpc = esd->GetNumberOfTracks();
  for (Int_t itrack = 0; itrack < endtpc; itrack++) 
    {
      AliESDtrack * track = esd->GetTrack(itrack) ; // retrieve track
      iemcalMatch = track->GetEMCALcluster();
      if(iemcalMatch >= 0) matchedTrack[iemcalMatch] = itrack;
    } 
  
  //########################################
  //##############Fill CaloClusters#############
  //########################################
  esd->SetNumberOfEMCALClusters(nClusters);
  for (Int_t iClust = 0 ; iClust < nClusters ; iClust++) 
    {
      const AliEMCALRecPoint * clust = (const AliEMCALRecPoint*)clusters->At(iClust);
      //if(clust->GetClusterType()== AliESDCaloCluster::kEMCALClusterv1) nRP++; else nPC++;

      // Get information from EMCAL reconstruction points
      Float_t xyz[3];
      TVector3 gpos;
      clust->GetGlobalPosition(gpos);
      for (Int_t ixyz=0; ixyz<3; ixyz++)
	xyz[ixyz] = gpos[ixyz];
      Float_t elipAxis[2];
      clust->GetElipsAxis(elipAxis);
      //Create digits lists
      Int_t cellMult = clust->GetMultiplicity();
      //TArrayS digiList(digitMult);
      Float_t *amplFloat = clust->GetEnergiesList();
      Int_t   *digitInts = clust->GetAbsId();
      TArrayS absIdList(cellMult);
      //Uncomment when unfolding is done
      //TArrayD fracList(cellMult);
      
      Int_t newCellMult = 0;
      for (Int_t iCell=0; iCell<cellMult; iCell++) 
	{
	  if (amplFloat[iCell] > 0) {
	    absIdList[newCellMult] = (UShort_t)(digitInts[iCell]);
	    //Uncomment when unfolding is done
	    //fracList[newCellMult] = amplFloat[iCell]/emcCells.GetCellAmplitude(digitInts[iCell]);
	    newCellMult++;
	  }
	}

      absIdList.Set(newCellMult);
      //Uncomment when unfolding is done
      //fracList.Set(newCellMult);

      if(newCellMult > 0) { // accept cluster if it has some digit
	nClustersNew++;
	//Primaries
	Int_t  parentMult  = 0;
	Int_t *parentList =  clust->GetParents(parentMult);
	// fills the ESDCaloCluster
	AliESDCaloCluster * ec = new AliESDCaloCluster() ;
	ec->SetClusterType(AliESDCaloCluster::kEMCALClusterv1);
	ec->SetPosition(xyz);
	ec->SetE(clust->GetEnergy());
	ec->SetNCells(newCellMult);
	//Change type of list from short to ushort
	UShort_t *newAbsIdList  = new UShort_t[newCellMult];
	//Uncomment when unfolding is done
	//Double_t *newFracList  = new Double_t[newCellMult];
	for(Int_t i = 0; i < newCellMult ; i++) 
	  {
	    newAbsIdList[i]=absIdList[i];
	    //Uncomment when unfolding is done
	    //newFracList[i]=fracList[i];
	  }
	
	ec->SetCellsAbsId(newAbsIdList);
	//Uncomment when unfolding is done
	//ec->SetCellsAmplitudeFraction(newFracList);
	ec->SetClusterDisp(clust->GetDispersion());
	ec->SetClusterChi2(-1); //not yet implemented
	ec->SetM02(elipAxis[0]*elipAxis[0]) ;
	ec->SetM20(elipAxis[1]*elipAxis[1]) ;
	ec->SetTOF(clust->GetTime()) ; //time-of-fligh
	ec->SetNExMax(clust->GetNExMax());          //number of local maxima
	TArrayI arrayTrackMatched(1);// Only one track, temporal solution.
	arrayTrackMatched[0]= matchedTrack[iClust];
	ec->AddTracksMatched(arrayTrackMatched);
	
	TArrayI arrayParents(parentMult,parentList);
	ec->AddLabels(arrayParents);
	
	// add the cluster to the esd object
	esd->AddCaloCluster(ec);
	delete ec;
	delete [] newAbsIdList ;
	//delete [] newFracList ;
      }
    } // cycle on clusters
  
  delete [] matchedTrack;
  
  esd->SetNumberOfEMCALClusters(nClustersNew);
  //if(nClustersNew != nClusters)
  //printf(" ##### nClusters %i -> new %i ##### \n", nClusters, nClustersNew );
  
  //Fill ESDCaloCluster with PID weights
  AliEMCALPID *pid = new AliEMCALPID;
  //pid->SetPrintInfo(kTRUE);
  pid->SetReconstructor(kTRUE);
  pid->RunPID(esd);
  delete pid;

  delete digits;
  delete clusters;

  // printf(" ## AliEMCALReconstructor::FillESD() is ended : ncl %i -> %i ### \n ",nClusters, nClustersNew); 
  return kTRUE;
}
