void ApplyAlignment(const char* fileMA)
{
  // apply ITSU alignment from the file
  //
  if (!gGeoManager) {
    if (!gSystem->AccessPathName("geometry.root")) {
      printf("Loading geometry.root from current directory\n");
      AliGeomManager::LoadGeometry("geometry.root"); //load geom from default CDB storage      
    }
    else {
      printf("No geometry in memory and not geometry.root in current directory\n");
      return;
    }
  }
  else {
    if (gGeoManager->IsLocked()) {
      printf("There is geometry in memory but it is locked");
      return;
    }
  }
  //
  TFile* fl = TFile::Open(fileMA);
  if (!fl) {
    printf("Failed to open misalignments file %s\n",fileMA);
  }
  TClonesArray* arr = (TClonesArray*)fl->Get("ITSUAlignObjs");
  if (!arr) {
    AliCDBEntry* cdbe = (AliCDBEntry*) fl->Get("AliCDBEntry");
    if (!cdbe) {
      printf("File %s does not contain recognizable misalignment\n",fileMA);
      return;
    }
    arr = (TClonesArray*)cdbe->GetObject();
  }
  //
  if (!arr->IsA()==TClonesArray::Class()) {
    printf("The object in %s is not TClonesArray\n",fileMA);
    return;
  }
  printf("Applying misalignment from %s geometry in memory\n",fileMA);
  AliGeomManager::ApplyAlignObjsToGeom(*arr);
  gGeoManager->LockGeometry();
  //
}
