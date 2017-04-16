
///
/// \file  TestGeometry.C
/// \ingroup EMCAL_AlignDB
/// \brief Print geometry info from simulation files
///
/// Print geometry info from simulation files
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__)
#include <Riostream.h>
#include <TGeoManager.h>

#include "AliRun.h"
#include "AliEMCALLoader.h"
#include "AliEMCALGeometry.h"
#endif

///
/// Main method.
///
void TestEMCALGeometry()
{
  // Getting EMCAL Detector and Geometry.
  AliRunLoader *rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),"read");
 
  if (rl == 0x0)
      cout<<"Can not instatiate the Run Loader"<<endl;

  rl->LoadgAlice();//Needed to get geometry

  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>
    (rl->GetDetectorLoader("EMCAL"));

  //AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance();
  //AliEMCALGeometry *geom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
  TGeoManager::Import("geometry.root");

  gGeoManager->CheckOverlaps();
  gGeoManager->PrintOverlaps();

  /*
  AliRun * alirun   = rl->GetAliRun(); // Needed to get Geometry
  AliEMCAL * emcal  = (AliEMCAL*)alirun->GetDetector("EMCAL");
  AliEMCALGeometry * geom = emcal->GetGeometry();
  
  if (geom==0)
    cout<<"Did not get geometry from EMCALLoader"<<endl;

  geom->PrintGeometry();
  */

}

