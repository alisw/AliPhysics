void
CheckOverlaps(const char* file="geometry.root", 
	      Bool_t align=kFALSE, Bool_t sample=kTRUE)
{
  TObjArray* checked = new TObjArray();
  
  AliGeomManager::LoadGeometry(file);
  if (align)
    AliGeomManager::ApplyAlignObjsToGeom("FMDfullMisalignment.root", 
					 "FMDAlignment");
  TObjArray*        l = gGeoManager->GetListOfPhysicalNodes();
  TIter             next(l);
  TGeoPhysicalNode* pn = 0;
  TGeoVolume*       v  = 0;
  while ((pn = static_cast<TGeoPhysicalNode*>(next()))) { 
    pn->cd();
    gGeoManager->CdUp();
    v = gGeoManager->GetCurrentVolume();
    if (checked->FindObject(v)) continue;
    
    std::cout << "Checking " << v->GetName() << std::endl;
    v->CheckOverlaps(0.01);
    Int_t n = gGeoManager->GetListOfOverlaps()->GetEntriesFast();
    if (n) { 
      gGeoManager->GetListOfOverlaps()->ls();
    }
    checked->Add(v);
    
    if (!sample) continue;

    // gGeoManager->ClearOverlaps();
    gGeoManager->SetCheckingOverlaps();
    TGeoNode*    start = gGeoManager->GetCurrentNode();
    TGeoVolume*  vol   = start->GetVolume();
    TGeoIterator gnext(vol);
    TGeoNode*    node;
    TString      path;
    while ((node = gnext())) {
      gnext.GetPath(path);
      // std::cout << " Checking: " <<  path.Data() << std::endl;
      node->GetVolume()->CheckOverlaps(0.01,"s");
    }
    gGeoManager->SetCheckingOverlaps(kFALSE);
    n = gGeoManager->GetListOfOverlaps()->GetEntriesFast();
    if (n) {
      gGeoManager->GetListOfOverlaps()->ls();
      pn->Print();
    }
  }
}

