/**************************************************************************
 * To Create PMD Alignment ObJect  using the class AliPMDMisAligner
 * sjena@cern.ch
 * Mon Nov 22 19:54:27 CET 2010
 *                     
 **************************************************************************/

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
