void initGeom()
{
    // Macro replacing the gAlice->Init() to initialize the geometry
    // Suggested by Raffaele Grosso <Raffaele.Grosso@cern.ch>
    
    AliCDBManager* man = AliCDBManager::Instance();
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    man->SetRun(0);
    
    gROOT->LoadMacro("Config.C");
    gInterpreter->ProcessLine(gAlice->GetConfigFunction());
    gAlice->GetMCApp()->Init();
}
