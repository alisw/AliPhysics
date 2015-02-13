/**************************************************************************
  // Create TClonesArray of zero misalignment objects for PMD
  //
  // Macro to randomly displace the 4 sectors of the PMD
  // in each plane. Each sector (to be misaligned) 
  // of PMD houses the following :
  // (a) 6 modules of preshower plane
  // (b) 6 modules of veto plane
  // (c) The FEE boards on back plane of each module
  // (d) 6 modules of convertor plates
  // The clustering is done module-wise
  // The actual amount displacement will be provided
  // by the survey data and has to be converted into
  // displacement in x,y,z,theta, phi and psi 
   
  // Now specify the path of the module to be misaligned
  // as followed in the PMD geant
   
     _____________
    |    |        |
    | 1  |   3    |
    |    |________|
    |____|___|    |
    |        | 2  |
    |   4    |    |
    |________|____|
    
    // Misalignment Matrix is expected to be
    // same for sectors 1 and 4 
    // and for the sectors 2 and 3
    // As these will be mounted on the same
    // Steel plate 
 * sjena@cern.ch
 * Mon Nov 22 19:54:27 CET 2010
 *         

OCDB/PMD/Align/Data
            
 **************************************************************************/

void MakePMDZeroMisAlignment(){
  const char* macroname = "MakePMDZeroMisAlignment.C";

  //Create a TClonesArray of Align Object to store displacement Angles
  TClonesArray *array = new TClonesArray("AliAlignObjParams",10);
  TClonesArray &alobj = *array;
  
  Int_t iIndex=0; //  let all modules have index=0 in a layer with no LUT
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);
  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  Int_t i, j=0;

  for(i=1; i<=4; i++){
    TString snSector(Form("PMD/Sector%d",i));
    new(alobj[j++]) AliAlignObjParams(snSector.Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // Create a File to store the alignment data
    const char* filename = "PMDzeroMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"PMDAlignObjs","kSingleKey");
    f.Close();
  }else{
  // save in CDB storage
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    Info(macroname,"Saving alignment objects in CDB storage %s",
	 Storage.Data());
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Satyajit Jena");
    md->SetComment("Zero misalignment for PMD");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("PMD/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }
  array->Delete();

}
