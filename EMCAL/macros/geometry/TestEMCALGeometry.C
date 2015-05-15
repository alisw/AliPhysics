

// Test Macro, shows how to load Digits and Geometry, and how can we get 
// some of the parameters and variables.
// Author: Gustavo Conesa

void TestEMCALGeometry()
{
   
  // Getting EMCAL Detector and Geometry.
  
  AliRunLoader *rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),
			  "read");
 
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

