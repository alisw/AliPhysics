
void MakePMDAlignmentObjs()
{
  // To create PMD alignment objects using the class AliPMDMisAligner
  AliGeomManager::LoadGeometry("geometry.root");
  AliPMDMisAligner *mA = new AliPMDMisAligner();
  mA->SetMisalType("residual"); // or "ideal" or "full"
  TClonesArray* array = mA->MakeAlObjsArray();
  array->Print();


  return;
}
