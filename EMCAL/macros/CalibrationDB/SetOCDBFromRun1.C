///
/// \file SetOCDBFromRun1.C
/// \ingroup EMCAL_CalibDB
/// \brief Set OCDB for EMCAL energy calibration starting from Run1
///
/// Create OCDB file for a given period from already existing OCDB file
/// include the online calibration parameters
///
/// \author Gustavo Conesa Balbastre, <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-CNRS
///

#if !defined(__CINT__)
#include <TSystem.h>
#include <TString.h>
#include <TGrid.h>

#include <Riostream.h>

#include "AliEMCALGeometry.h"
#include "AliEMCALCalibData.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#endif

///
/// Main method
///
/// \param year: year to set geometry and run range
/// \param printAll: verbosity checks 
///
void SetOCDBFromRun1(Int_t year = 2010, Bool_t printAll = kFALSE)
{  
  TGrid::Connect("alien://");
  
  Int_t run = 182325; //2012
  if(year == 2010) run = 134908;
  if(year == 2011) run = 159582;
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  man->SetRun(run);
  AliCDBStorage *storage = man->GetDefaultStorage();
  
  // Instantiate EMCAL geometry for the first time
  AliEMCALGeometry * geom;
  if     (year == 2010) geom = AliEMCALGeometry::GetInstance("EMCAL_FIRSTYEARV1"); // 2010
  else                  geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");  // 2011-2012-2013
//else                  geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM"); // Run2

  const Int_t nSM = geom->GetNumberOfSuperModules();
  
  // Get the final OCDB object
  
  AliEMCALCalibData* cparam = (AliEMCALCalibData*) (storage->Get("EMCAL/Calib/Data", run)->GetObject());

  // Access OCDB file with the first version of the calibration
  TString        first = "Run177115_999999999_v2_s0.root";
  if(year==2010) first = "Run113461_999999999_v3_s0.root";
  if(year==2011) first = "Run144484_999999999_v3_s0.root";
  
  TFile * f = TFile::Open(Form("alien:///alice/data/%d/OCDB/EMCAL/Calib/Data/%s",year,first.Data()),"READ");
  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
  AliEMCALCalibData* cparam1 = (AliEMCALCalibData*)  cdb->GetObject();
  
  // New OCDB container
  AliEMCALCalibData *cparamnew=new AliEMCALCalibData("EMCAL");
  
  // Do the comparison
  Float_t param  = -1;
  Float_t param1 = -1;
  Int_t iCol = -1, iRow = -1, iSM =-1, iMod = -1,iIphi =-1,iIeta = -1;
  for(Int_t i=0;i < nSM*24*48; i++)
  {
    //printf("AbsID %d\n",i);
    geom->GetCellIndex(i,iSM,iMod,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSM,iMod, iIphi, iIeta,iRow,iCol);
    
    Float_t param = -1;
    if( cparam  ) param  = cparam ->GetADCchannel(iSM,iCol,iRow);
    
    Float_t param1 = -1;
    if( cparam1 ) param1 = cparam1->GetADCchannel(iSM,iCol,iRow);
    
    if    (printAll)
      printf("ID %d, col %d, row %d, sm %d  final %1.4f, first %1.4f\n",
             i,iCol,iRow,iSM,param, param1);
    cparamnew->SetADCchannel      (iSM,iCol,iRow,param );
    cparamnew->SetADCchannelOnline(iSM,iCol,iRow,param1);
  }
  
  // Create OCDB File
  AliCDBMetaData md;
  md.SetComment("Calibration after calibration with pi0, store also first online calibration");
  md.SetBeamPeriod(0);
  md.SetResponsible("Gustavo Conesa");
  md.SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  
  // Careful, select here the first run where this calibration is valid
  Int_t            firstRun = 172439; // 2012-13
  if(year == 2010) firstRun = 113461;
  if(year == 2011) firstRun = 144484;
  
  AliCDBId id("EMCAL/Calib/Data",firstRun,AliCDBRunRange::Infinity()); // create in EMCAL/Calib/Data DBFolder
  
  AliCDBManager* man2 = AliCDBManager::Instance();
  AliCDBStorage* loc = man2->GetStorage(Form("local://%d",year));
  loc->Put(cparamnew, id, &md);
}


