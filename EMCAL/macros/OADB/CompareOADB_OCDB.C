///
/// \file CompareOADB_OCDB.C
/// \ingroup EMCAL_OADB
/// \brief Compare contents in EMCal OADB and OCDB
///
/// Compare the contents of the calibration, bad channels, etc in OCDB and OADB
/// You need connexion to the grid to run this macro
/// The OADB file can be in a place different than the default in aliroot
///
/// \author Gustavo Conesa Balbastre, <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-CNRS
///

#if !defined(__CINT__)
#include <TString.h>
#include <TH2F.h>
#include <TObjArray.h>
#include <TGrid.h>

#include "AliCaloCalibPedestal.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALCalibData.h"
#include "AliEMCALCalibTime.h"

#include "AliCDBManager.h"
#include "AliOADBContainer.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#endif

// Common parameters in different methods:
Int_t   fRun      ; ///< reference run number
TString fPathOADB ; ///< path to AODB file
TString fPass     ; ///< pass name
Bool_t  fPrintAll ; ///< debug level

AliEMCALGeometry * fGeom    ; ///< EMCAL geometry
AliCDBStorage    * fStorage ; ///< OCDB

///
/// Check energy recalibration.
//___________________________
void CheckEnergyCalibration()
{  
  AliOADBContainer *contRF=new AliOADBContainer("");
  contRF->InitFromFile(Form("%s/EMCALRecalib.root",fPathOADB.Data()),"AliEMCALRecalib");
  
  const Int_t nSM = fGeom->GetNumberOfSuperModules();
  
  // Get the OCDB object
  
  AliEMCALCalibData* cparam = (AliEMCALCalibData*) (fStorage->Get("EMCAL/Calib/Data", fRun)->GetObject());

//  // Access directly the OCDB file and not the latest version
//  TFile * f = TFile::Open("alien:///alice/data/2011/OCDB/EMCAL/Calib/Data/Run144484_999999999_v2_s0.root","READ");
//  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
//  AliEMCALCalibData* cparam2 = (AliEMCALCalibData*)  cdb->GetObject();

  // Get the OADB object
  
  TH2F *h[12];
  
  TObjArray *recal=(TObjArray*)contRF->GetObject(fRun);
  
  if(!recal)
  {
    printf("Energy recalibration OADB not found  3\n");
    return;
  }
  
  TObjArray *recalpass=(TObjArray*)recal->FindObject(fPass);
  
  if(!recalpass)
  {
    printf("Energy recalibration OADB not found 2\n");
    return;
  }
  
  TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
  
  if(!recalib)
  {
    printf("Energy recalibration OADB not found 1\n");
    return;
  }
  
  for (Int_t i=0; i < nSM; ++i)
  {
    h[i] = (TH2F*)recalib->FindObject(Form("EMCALRecalFactors_SM%d",i));
    if (!h[i])
    {
      printf("Could not load EMCALRecalFactors_SM%d\n",i);
      continue;
    }
  }
  
  // Do the comparison
  Float_t paramOCDB = -1;
  Float_t paramOADB = -1;
  Int_t nDiff = 0;
//  Int_t nDiff2 = 0;
  Int_t iCol = -1, iRow = -1, iSM =-1, iMod = -1,iIphi =-1,iIeta = -1;
  for(Int_t i=0;i < nSM*24*48; i++)
  {
    //printf("AbsID %d\n",i);
    fGeom->GetCellIndex(i,iSM,iMod,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSM,iMod, iIphi, iIeta,iRow,iCol);
    
    Float_t paramOCDB = -1;
    if(cparam) paramOCDB = cparam->GetADCchannel(iSM,iCol,iRow);

//    Float_t paramOCDB2 = -1;
//    if(cparam) paramOCDB2 = cparam2->GetADCchannel(iSM,iCol,iRow);
    
    Float_t paramOADB = -1;
    if(h[iSM]) paramOADB = h[iSM]->GetBinContent(iCol,iRow);
    paramOADB*=0.0162; // Transformation into OCDB parameter, it will not work for all channels
    
    if    (fPrintAll)
    printf("Recalibration parameter: ID %d, col %d, row %d, sm %d  OCDB %1.4f, OADB %1.4f\n",
           i,iCol,iRow,iSM,paramOCDB, paramOADB);
    else if(paramOADB > 0)
    {
      if ( paramOCDB/paramOADB > 1.01 || paramOCDB/paramOADB < 0.99 )
      {
        printf("DIFFERENT recalibration parameter: ID %d, col %d, row %d, sm %d  OCDB %1.4f, OADB %1.4f, ratio OCDB/OADB %1.4f\n",//, old OCDB param %1.4f\n",
               i,iCol,iRow,iSM,paramOCDB, paramOADB,paramOCDB/paramOADB);//,paramOCDB2);
        nDiff++;
      }
      else if(paramOADB <= 0)
      {
        printf("DIFFERENT recalibration parameter: ID %d, col %d, row %d, sm %d  OCDB %1.4f, OADB %1.4f\n",//, old OCDB param %1.4f\n",
               i,iCol,iRow,iSM,paramOCDB, paramOADB);//,paramOCDB2);
        nDiff++;
      }
      
//      if(TMath::Abs(paramOCDB2-0.0162)> 0.0001)
//      {
//        printf("\t Different initial calib! %1.4f\n",paramOCDB2);
//        nDiff2++;
//      }
    }
  }
  
  if(!fPrintAll) printf("Total number of different channels %d / %d\n",nDiff,nSM*24*48);//, origin %d\n",nDiff,nSM*24*48,nDiff2);
}

///
/// Check energy online calibration
//___________________________
void CheckEnergyOnlineCalibration()
{  
  AliOADBContainer *contRF=new AliOADBContainer("");
  contRF->InitFromFile(Form("%s/EMCALCalibOnlineRef.root",fPathOADB.Data()),"AliEMCALRecalib");
  
  const Int_t nSM = fGeom->GetNumberOfSuperModules();
  
  // Get the OCDB object
  // Access directly the OCDB file and not the latest version
  TFile * f = 0;
  if     (fRun < 140000)
    f = TFile::Open("alien:///alice/data/2010/OCDB/EMCAL/Calib/Data/Run113461_999999999_v3_s0.root","READ");
  else if(fRun < 171000)
    f = TFile::Open("alien:///alice/data/2011/OCDB/EMCAL/Calib/Data/Run144484_999999999_v3_s0.root","READ");
  else if(fRun < 198000)
    f = TFile::Open("alien:///alice/data/2012/OCDB/EMCAL/Calib/Data/Run177115_999999999_v2_s0.root","READ");
  else {
    printf("Run not available\n");
    return;
  }
  
  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
  AliEMCALCalibData* cparam = (AliEMCALCalibData*)  cdb->GetObject();
  
  // Get the OADB object
  
  TH2F *h[12];
  
  TObjArray *cal=(TObjArray*)contRF->GetObject(fRun);
  
  if(!cal)
  {
    printf("Energy online calibration OADB not found  3\n");
    return;
  }

  for (Int_t i=0; i < nSM; ++i)
  {
    h[i] = (TH2F*)cal->FindObject(Form("EMCALRecalFactors_SM%d",i));
    if (!h[i])
    {
      printf("Could not load EMCALRecalFactors_SM%d\n",i);
      continue;
    }
  }
  
  // Do the comparison
  Float_t paramOCDB = -1;
  Float_t paramOADB = -1;
  Int_t nDiff = 0;
  //  Int_t nDiff2 = 0;
  Int_t iCol = -1, iRow = -1, iSM =-1, iMod = -1,iIphi =-1,iIeta = -1;
  for(Int_t i=0;i < nSM*24*48; i++)
  {
    //printf("AbsID %d\n",i);
    fGeom->GetCellIndex(i,iSM,iMod,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSM,iMod, iIphi, iIeta,iRow,iCol);
    
    Float_t paramOCDB = -1;
    if(cparam) paramOCDB = cparam->GetADCchannel(iSM,iCol,iRow);
    
    Float_t paramOADB = -1;
    if(h[iSM]) paramOADB = h[iSM]->GetBinContent(iCol,iRow);
    
    if    (fPrintAll)
    printf("Online Calib: ID %d, col %d, row %d, sm %d  OCDB %1.4f, OADB %1.4f\n",
           i,iCol,iRow,iSM,paramOCDB, paramOADB);
    
    if (TMath::Abs(paramOCDB-paramOADB) > 0.001)//  || (TMath::Abs(paramOCDB-0.0162)> 0.001))
    {
      printf("Online Calib DIFFERENT: ID %d, col %d, row %d, sm %d  OCDB %1.4f, OADB %1.4f, diff OCDB-OADB %1.4f\n",
             i,iCol,iRow,iSM,paramOCDB, paramOADB,paramOCDB-paramOADB);
      nDiff++;
    }
    
    
  }
  
  if(!fPrintAll) printf("Total number of different channels %d / %d\n",nDiff,nSM*24*48);//, origin %d\n",nDiff,nSM*24*48,nDiff2);
}

///
/// Get the OCDB bad channels and compare to the OADB ones.
//____________________
void CheckBadChannels()
{  
  //const Int_t nSM = static_const (fGeom->GetNumberOfSuperModules());
  const Int_t nSM = static_cast<const Int_t> (fGeom->GetNumberOfSuperModules());

  // Access OCDB histograms
  AliCaloCalibPedestal* caloped  = (AliCaloCalibPedestal*) (fStorage->Get("EMCAL/Calib/Pedestals", fRun)->GetObject());
    
  // Access directly the OCDB file and not the latest version
  //TFile * f = TFile::Open("alien:///alice/data/2011/OCDB/EMCAL/Calib/Pedestals/Run145954_146856_v3_s0.root","READ");
  //AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
  //AliCaloCalibPedestal * caloped = (AliCaloCalibPedestal *) cdb->GetObject();  
  
  TObjArray map = caloped->GetDeadMap();

  // Access OADB histograms
  TH2I *hbm[12];
  
  AliOADBContainer *contBC=new AliOADBContainer("");
  
  contBC->InitFromFile(Form("%s/EMCALBadChannels.root",fPathOADB.Data()),"AliEMCALBadChannels"); 
  
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(fRun);
  
  if(!arrayBC)
  {
    printf("--- Bad map not available \n");
    return;
  }
  
  for (Int_t i=0; i < nSM; ++i)
  {
    hbm[i] = (TH2I*) arrayBC->FindObject(Form("EMCALBadChannelMap_Mod%d",i));
    
    if (!hbm[i])
    {
      printf("Can not get EMCALBadChannelMap_Mod%d\n",i);
      continue;
    }
  }
  
  Int_t badMapOCDB = -1;
  Int_t badMapOADB = -1;
  Int_t iCol = -1, iRow = -1, iSM =-1, iMod = -1,iIphi =-1,iIeta = -1;    
  for(Int_t i = 0; i < nSM*24*48; i++)
  {
    //printf("AbsID %d\n",i);
    fGeom->GetCellIndex(i,iSM,iMod,iIphi,iIeta); 
    fGeom->GetCellPhiEtaIndexInSModule(iSM,iMod, iIphi, iIeta,iRow,iCol);
    if(map.At(iSM))badMapOCDB = ((TH2F*)map.At(iSM))->GetBinContent(iCol,iRow);
    else badMapOCDB  = -1;
    if(hbm[iSM])   badMapOADB = hbm[iSM]->GetBinContent(iCol,iRow);
     else badMapOADB = -1;

    if    (fPrintAll && badMapOCDB > 0)
      printf("STATUS: ID %d, col %d, row %d, sm %d  OCDB %d, OADB %d\n",
             i,iCol,iRow,iSM,badMapOCDB, badMapOADB);
    else if(badMapOCDB!=badMapOADB)
      printf("DIFFERENT STATUS: ID %d, col %d, row %d, sm %d  OCDB %d, OADB %d\n",
             i,iCol,iRow,iSM,badMapOCDB, badMapOADB);

  }
}

///
/// Get the OCDB time calibration shifts and compare to the OADB ones.
//________________________
void CheckTimeCalibration()
{
  // Get the OCDB object

  AliEMCALCalibTime* cparam = (AliEMCALCalibTime*) (fStorage->Get("EMCAL/Calib/Time", fRun)->GetObject());

  // Access directly the OCDB file and not the latest version
//  TFile * f = TFile::Open("2010/OCDB/EMCAL/Calib/Time/Run122195_126437_v0_s1.root","READ");
//
//  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
//  AliEMCALCalibTime* cparam = (AliEMCALCalibTime*)  cdb->GetObject();

  if(!cparam)
  {
    printf("OCDB not available\n");
    return;
  }

  // Access OADB
  AliOADBContainer *contTRF=new AliOADBContainer("");
  contTRF->InitFromFile(Form("%s/EMCALTimeCalib.root",fPathOADB.Data()),"AliEMCALTimeCalib");
  
  TObjArray *trecal=(TObjArray*)contTRF->GetObject(fRun); 
  
  if(!trecal)
  {
    return;
    printf("\t array not available\n");
  }
  
  TObjArray *trecalpass=(TObjArray*)trecal->FindObject(fPass);
  
  if(!trecalpass)
  {
    return;
    printf("\t pass not available\n");
  }
  
  Int_t iTower  = -1;
  Int_t iIphi   = -1, iIeta   = -1;
  Int_t iSupMod = -1, iSupModMax = -1;
  Int_t iphi = -1, ieta =-1;
  
  for (Int_t ibc = 0; ibc < 4; ++ibc) 
  {
    TH1F *h = (TH1F*)trecalpass->FindObject(Form("hAllTimeAvBC%d",ibc));
    
    printf("\t BC%d -  %p\n",ibc,h);
    
    if (!h) 
    {
      printf("Could not load hAllTimeAvBC%d\n",ibc);
      continue;
    }
    
    // Fill params here
        
    const Int_t nSM = static_cast<const Int_t> (fGeom->GetNumberOfSuperModules());

    for(Int_t absId = 0; absId < nSM*24*48; absId++)
    {
      fGeom->GetCellIndex(absId,iSupMod,iTower,iIphi,iIeta); 
      fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);    
      
      Float_t paramOADB = h->GetBinContent(absId);
            
      Float_t paramOCDB = cparam->GetTimeChannel(iSupMod,ieta,iphi,ibc);
      
      
      if    (fPrintAll && TMath::Abs(paramOCDB-paramOADB) < 0.001)
        printf("OK: ID %d, col %d, row %d, sm %d  OCDB %2.3f, OADB %2.3f\n",
               absId,ieta,iphi,iSupMod,paramOCDB,paramOADB);
      else if(TMath::Abs(paramOCDB-paramOADB) > 0.001)
        printf("DIFFERENT VALUES: ID %d, col %d, row %d, sm %d  OCDB %2.3f, OADB %2.3f\n",
               absId,ieta,iphi,iSupMod,paramOCDB,paramOADB);
    }
    
    // end
  } // bunch crossing loop
}

///
/// Main method
/// Executes the for the selected run and parameter
///
/// \param run: run number
/// \param pathOADB: path location of OADB file
/// \param checkObject: Type of comparison, 1-bad channels, 2-energy, 3-energy online, 4-Time
/// \param pass: for some objects, pass string is needed
/// \param printAll: option to print messages
///
void CompareOADB_OCDB(Int_t run = 196000, TString pathOADB = "$ALICE_PHYSICS/OADB/EMCAL",
                      Int_t checkObject = 1, TString pass = "pass2", Bool_t printAll = kFALSE)
{  
  fRun      = run;
  fPathOADB = pathOADB;
  fPass     = pass;
  fPrintAll = printAll;
  
  printf("*** Settings: Run %d; object %d; pass %s; \n"
         "              path %s\n",
           fRun,checkObject,fPass.Data(),fPathOADB.Data() );
  
  gSystem->Load("libOADB");
  TGrid::Connect("alien://");
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  //man->SetDefaultStorage("local://2013/OCDB");
  man->SetRun(fRun);
  
  fStorage = man->GetDefaultStorage();
  
  // Instantiate EMCAL geometry for the first time
  
  fGeom = 0;
  if     (fRun < 140000) fGeom = AliEMCALGeometry::GetInstance("EMCAL_FIRSTYEARV1"); // 2010
  else if(fRun < 198000) fGeom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");  // 2011-2012-2013
  else                   fGeom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM"); // Run2
  
  if     (checkObject == 0) CheckBadChannels            ();
  else if(checkObject == 1) CheckEnergyCalibration      ();
  else if(checkObject == 2) CheckEnergyOnlineCalibration();
  if     (checkObject == 3) CheckTimeCalibration        ();
  else printf("non existing object option\n");
  
  printf("*** Comparisons ended *** \n");
}  

