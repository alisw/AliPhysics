///
/// \file  TestMatrices.C
/// \ingroup EMCAL_AlignDB
/// \brief Print alignment matrices
///
/// Macro to print the aligment matrices stored in the OCDB depending on year
/// These parameters are used during reconstruction and simulation
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__)
#include <Riostream.h>

#include "AliLog.h"
#include "AliAlignObjParams.h"
#include "AliEMCALGeometry.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif

///
/// Main method, access local or grid file, print its content.
///
/// \param year: bool 0, year 2010; 1, year 2011
///
void TestMatrices(Int_t year = 2011)
{
  AliLog::SetClassDebugLevel("AliCDBManager",1);

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
      
  AliEMCALGeometry *geo;

  if      (year == 2011)
  {
    man->SetSpecificStorage("EMCAL/Align/Data",             "alien://folder=/alice/data/2011/OCDB");
    man->SetRun(146805);
    geo = AliEMCALGeometry::GetInstanceFromRunNumber(146805);
  }
  else if (year == 2010)
  {
    man->SetSpecificStorage("EMCAL/Align/Data",             "alien://folder=/alice/data/2010/OCDB");
    man->SetRun(137366);
    geo = AliEMCALGeometry::GetInstanceFromRunNumber(137366);
  }
  else if (year == 2012 || year == 2013)
  {
    man->SetSpecificStorage("EMCAL/Align/Data",             "alien://folder=/alice/data/2012/OCDB");
    man->SetRun(190000);
    geo = AliEMCALGeometry::GetInstanceFromRunNumber(190000);
  }
  else
  {
    man->SetSpecificStorage("EMCAL/Align/Data",             "alien://folder=/alice/data/2015/OCDB");
    man->SetRun(220000);
    geo = AliEMCALGeometry::GetInstanceFromRunNumber(220000);
  }
  
  AliGeomManager::LoadGeometry();
  AliGeomManager::ApplyAlignObjsFromCDB("EMCAL");
  
  for (Int_t i=0;i<(geo->GetEMCGeometry())->GetNumberOfSuperModules();++i)
    geo->GetMatrixForSuperModule(i)->Print();

  AliEMCALEMCGeometry *emc = geo->GetEMCGeometry();
  
  Double_t phimin = emc->GetArm1PhiMin();
  Double_t phimax = emc->GetArm1PhiMax();
  
  cout << phimin << " " << phimax << endl;
  
  emc->PrintGeometry();
}
