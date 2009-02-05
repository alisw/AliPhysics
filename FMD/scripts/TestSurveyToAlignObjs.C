void
TestSurveyToAlignObjs(Int_t det=1, Bool_t cdbStore=false)
{
  if (det < 1 || det > 2) { 
    Error("TestSurveyToAlignObjs", "Invalid detector %d (must be 1,2, or 3)", 
	  det);
    return;
  }

  const char* files[] = { 
    "$ALICE_ROOT/FMD/Survey_943928_FMD.txt", 
    "$ALICE_ROOT/FMD/Survey_976326_FMD.txt", 
    0 
  };
  
  
  AliGeomManager::LoadGeometry("geometry.root");

  AliFMDSurveyToAlignObjs convert;
  if (!convert.LoadSurveyFromLocalFile(files[det-1])) { 
    Error("TestSurveyToAlignObjs", "Failed to load %s", files[det-1]);
    return;
  }
  convert.Run();
  convert.CreateAlignObjs();
  convert.GetAlignObjArray()->Print();

  if (!cdbStore) 
    convert.StoreAlignObjToFile("FMD_Survey.root", "FMD");
  else 
    convert.StoreAlignObjToCDB("FMD/Align/Data", "FMD");
}

