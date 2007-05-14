TList medias;

void
ConvertGeom()
{
  
  std::ostream output("foo.C");
  TGeoVolume*  top = gGeoManager->GetTopVolume();
  TGeoIterator next(top);
  TGeoNode*    node = 0;
  
  // Iterate through all nodes, and write out the mediums used. 
  while ((node = next())) {
    TGeoMedium* med = node->GetMedium();
    if (medias->Find(med)) {
      Info("ConvertGeom", "Already has medium %s", med->GetName());
      continue;
    }
    medias->Add(med);
    output << "  {\n"
	   << "    Double_t p[] = { " << std::flush;
    for (Int_t i = 0; i < 10; i++) {
      if (i != 0) output << ", ";
      output << med->GetPar(i);
    }
    output << std::endl;
  }
}

    
    
