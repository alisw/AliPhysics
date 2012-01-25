int runCreatePhysicsSelection()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWGLFtotEt");
  gInterpreter->GenerateDictionary("std::map<int, AliPhysicsSelection*>", "AliPhysicsSelection.h;map")  ;
  gInterpreter->GenerateDictionary("std::pair<int, AliPhysicsSelection*>", "AliPhysicsSelection.h;utility");
  
  gROOT->LoadMacro("CreatePhysicsSelection.C+g");
  
  Init();
  //GetPhysicsSelection();
  
  //for(Int_t i = 0; i < 10; i++)
  {
    UpdatePhysicsSelection(1, "physicsSelections.root");
  }
  
}
