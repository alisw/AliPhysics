void MakeMUONResMisAlignment(Bool_t volpaths = true,
                             Bool_t transforms = true, 
                             Bool_t svmaps = true)
// Macro for generating the geometry data files:
// (volpath.dat, transform.dat, svmap.dat)
// and local CDB storage with zero, residual and full misalignment
//
// The generated files do not replace the existing ones
// but have different names (with extension ".out").
//
//  Author: I. Hrivnacova, IPN Orsay
//
{
  // Initialize
  gAlice->Init("$ALICE_ROOT/MUON/Config.C");
  cout << "Init done " << endl;

  // Get MUON detector
  AliMUON* muon = (AliMUON*)gAlice->GetModule("MUON");
  if (!muon) {
    AliFatal("MUON detector not defined.");
    return 0;
  }  

  // Get geometry builder
  AliMUONGeometryBuilder* builder = muon ->GetGeometryBuilder();
  
  if (volpaths) {
    cout << "Generating volpath file ..." << endl;
    builder->GetTransformer()->WriteVolumePaths("volpath.dat.out");
  }  

  if (transforms) {
    cout << "Generating transformation file ..." << endl;
    builder->GetTransformer()->WriteTransformations("transform.dat.out");
  }  

  if (svmaps) {
    cout << "Generating svmaps file ..." << endl;
    builder->WriteSVMaps();
  }  

  cout << "Generating residual misalignment data in  MUON/ResMisAlignCDB/Data..." << endl;
  
  AliMUONGeometryMisAligner misAligner(0.0, 0.004, 0.0, 0.003, 0.0, 0.0023);
  AliMUONGeometryTransformer* newTransform 
    = misAligner.MisAlign(builder->GetTransformer(), true);
  TClonesArray* array = newTransform->GetMisAlignmentData();

  if(!gSystem->Getenv("$TOCDB")){
    // Create a File to store the alignment data
    TFile f("MUONresidualMisalignment.root","RECREATE");
    if(!f) {cerr<<"cannot open file for output\n";}
    
    f.cd();
    f.WriteObject(array,"MUONAlignObjs ","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdbManager = AliCDBManager::Instance();
    AliCDBStorage* storage = cdbManager->GetStorage(Storage);
    AliCDBMetaData* cdbData = new AliCDBMetaData();
    cdbData->SetResponsible("Dimuon Offline project");
    cdbData->SetComment("MUON alignment objects with residual misalignment");
    cdbData->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("MUON/Align/Data", 0, 9999999); 
    storage->Put(array, id, cdbData);
  }
  
  delete newTransform;
}   

