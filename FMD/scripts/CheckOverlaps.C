void
CheckOverlaps(Bool_t align=kTRUE, Bool_t sample=kTRUE)
{
  AliGeomManager::LoadGeometry("geometry.root");
  if (align)
    AliGeomManager::ApplyAlignObjsToGeom("FMDfullMisalignment.root", 
					 "FMDAlignment");
  TObjArray*        l = gGeoManager->GetListOfPhysicalNodes();
  TIter             next(l);
  TGeoPhysicalNode* pn = 0;
  TGeoVolume*       v  = 0;
  while ((pn = static_cast<TGeoPhysicalNode*>(next()))) { 
    pn->cd();
    v = gGeoManager->GetCurrentVolume();
    std::cout << "Checking " << v->GetName() << std::endl;
    v->CheckOverlaps(0.01);
    if (gGeoManager->GetListOfOverlaps()->GetEntriesFast()) 
      gGeoManager->GetListOfOverlaps()->ls();
    
    if (!sample) continue;

    gGeoManager->ClearOverlaps();
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
    if (gGeoManager->GetListOfOverlaps()->GetEntriesFast()) {
      gGeoManager->GetListOfOverlaps()->ls();
      pn->Print();
    }
  }
}

