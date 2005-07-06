///////////////////////////////////////////////////////
// Macro to create the galice.root file for         ///
// beam test analysis (using AliITSBeamTestDigitizer //
// class). Before using AliITSBeamTestDigitizer it   // 
// is necessary to run this macro to create the      //
// galice.root file.                                 //
// Options: Nov04 for integrated                     // 
//          ITS beam test November 2004              //
//                                                   //
//          Aug04 for SDD beam test August 2004      //
// E. Crescio: crescio@to.infn.it                    //
///////////////////////////////////////////////////////

void AliITSBeamTest(Char_t* opt="Nov04"){

  TString choice(opt);
  Bool_t aug04 = choice.Contains("Aug04");
  Bool_t nov04 = choice.Contains("Nov04");

  AliRunLoader* rl = AliRunLoader::Open("galice.root",
		    AliConfig::GetDefaultEventFolderName(),"recreate");

  gAlice->SetRunLoader(rl);  
  rl->SetEventFolderName();

  AliITS* bt;
  if(nov04){
    if(gGeoManager) delete gGeoManager;
    gGeoManager = new TGeoManager("ITSGeometry","ITS Simulation Geometry Manager");
    TGeoMaterial *vacmat = new TGeoMaterial("Vacuum",0,0,0);
    TGeoMedium   *vacmed = new TGeoMedium("Vacuum_med",1,vacmat);
    TGeoVolume *aLICE = gGeoManager->MakeBox("ALICE",vacmed,100.,100.,200.);
    gGeoManager->SetTopVolume(aLICE);  
    bt = new AliITSvBeamTestITS04("ITS","ITS beam test");
    bt->CreateGeometry();
    bt->Init();
  }
  if(aug04){
    bt = new AliITSvSDD03("ITS",2004);
    gSystem->Load("libgeant321");
    new TGeant3("C++ Interface to Geant3");
    if(strcmp(gMC->GetName(),"TGeant3")) {
	Fatal("Init","TGeant3 should be instantiated in advance");
	return;
    } 
    bt->CreateMaterials();
    bt->CreateGeometry();
    bt->Init();
  }
  gAlice->AddModule(bt);
  bt->SetDefaults();   
  rl->AddLoader(bt);
  // rl->MakeTree("E");  
  //rl->WriteHeader("OVERWRITE");
  rl->WriteRunLoader("OVERWRITE");
  rl->WriteAliRun("OVERWRITE"); 

}
