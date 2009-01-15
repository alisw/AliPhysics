void
TestSurveyToAlignObjs()
{
  AliGeomManager::LoadGeometry("geometry.root");

  AliFMDSurveyToAlignObjs convert;
  convert.LoadSurveyFromLocalFile("$ALICE_ROOT/FMD/Survey_XXX_FMD.txt");
  convert.Run();
}

