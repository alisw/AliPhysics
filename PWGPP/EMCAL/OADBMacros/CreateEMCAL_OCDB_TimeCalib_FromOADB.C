///
/// \file CreateEMCAL_OCDB_TimeCalib_FromOADB.C
/// \ingroup EMCAL_OADB
/// \brief Create time calibration OCDB from OADB
///
/// Compare the contents of the calibration, bad channels, etc in OCDB and AODB
/// You need connexion to the grid to run this macro
/// The OADB file can be in a place different than the default in aliroot
///
/// \author Gustavo Conesa Balbastre, <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-CNRS
///

//#if !defined(__CINT__)
#include <TH1F.h>
#include <TSystem.h>
#include <TString.h>
#include <TObjArray.h>

#include <Riostream.h>

#include "AliEMCALGeometry.h"
#include "AliEMCALCalibTime.h"
#include "AliOADBContainer.h"

#include "AliCDBManager.h"
#include "AliOADBContainer.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
//#endif

AliEMCALGeometry * geom    = 0; ///< Geometry pointer
AliOADBContainer * contTRF = 0; ///< OADB container

///
/// Create OADB file for a given period
///
/// \param runmin: minimum run interval
/// \param runmax: maximum run interval
/// \param period: period name added to comment
/// \param year: year to set SM number and file location
///
void CreatePeriod(Int_t runmin, Int_t runmax, TString period, Int_t year)
{
  printf("*** Create OCDB for period %d-%d ***\n",runmin,runmax);
 
  AliEMCALCalibTime *calibti=new AliEMCALCalibTime("EMCAL");
  
  TObjArray *trecal=(TObjArray*)contTRF->GetObject(runmin); 
  
  if(!trecal)
  {
    return;
    printf("\t array not available\n");
  }
  
  TObjArray *trecalpass=(TObjArray*)trecal->FindObject("pass3");
  
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
    
    Int_t nCells = h->GetNbinsX();
    //printf("N bins %d\n",nCells);
    
    Int_t nSM = 20;
    if(year < 2014) nSM = 10;
    
    for(Int_t absId = 0; absId < nCells; absId++)
    {
      geom->GetCellIndex(absId,iSupMod,iTower,iIphi,iIeta); 
      geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);    
      
      if(iSupMod >= nSM) continue;
      
      Float_t content = h->GetBinContent(absId);
      //printf("ID %d, col %d, row %d, sm %d  content %d\n",
      //       absId,ieta,iphi,iSupMod,content);
      calibti->SetTimeChannel(iSupMod,ieta,iphi,ibc,content);
    }
    
    // end
    
  } // bunch crossing loop
    
    
  AliCDBMetaData md;
  md.SetComment(Form("Time calibration parameters for period %s, runs %d-%d",period.Data(),runmin,runmax));
  //md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Gustavo Conesa");
  
  AliCDBId id("EMCAL/Calib/Time",runmin,runmax); // create in EMCAL/Calib/Time dbFolder
  
  AliCDBManager* man = AliCDBManager::Instance();
  AliCDBStorage* loc = man->GetStorage(Form("local://%d/OCDB",year));
  loc->Put(calibti, id, &md);
  
}

/// Main method
/// Create OCDB for different periods at the same time.
/// 
void CreateEMCAL_OCDB_TimeCalib_FromOADB()
{  
  // Instantiate EMCAL geometry for the first time
  geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM"); // Run2
  
  // OADB common stuff
  gSystem->Load("libOADB");
  
  contTRF=new AliOADBContainer("");
  
  TString oadbFilePathEMCAL = "$ALICE_PHYSICS/OADB/EMCAL" ;
  
  contTRF->InitFromFile(Form("%s/EMCALTimeCalib.root",oadbFilePathEMCAL.Data()),"AliEMCALTimeCalib");
  
  CreatePeriod(122195,126437,"LHC10d, pp, 7 TeV",2010); 
  CreatePeriod(136851,140000,"LHC10h, Pb-Pb, 2.76 TeV",2010); 
  CreatePeriod(144871,146459,"LHC11a, pp, 7 TeV",2011);
  CreatePeriod(146686,146860,"LHC11a, pp, 2.76 TeV",2011);
  CreatePeriod(148531,150629,"LHC11b, pp, 7 TeV",2011);
  CreatePeriod(151636,155384,"LHC11c, pp, 7 TeV",2011);
  CreatePeriod(156477,159635,"LHC11d, pp, 7 TeV",2011);
  CreatePeriod(160670,162740,"LHC11e, pp, 7 TeV",2011);
  CreatePeriod(162933,165746,"LHC11f, pp, 7 TeV",2011);
  CreatePeriod(166529,170593,"LHC11h, Pb-Pb, 2.76 TeV",2011);
  CreatePeriod(176326,193766,"LHC12a-i",2012);
  CreatePeriod(195344,197692,"LHC13b-g",2013);
}
