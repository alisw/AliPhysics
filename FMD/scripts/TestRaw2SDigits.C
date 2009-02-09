void
TestRaw2SDigits()
{
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(0);
  AliLog::SetModuleDebugLevel("FMD", 2);
  AliSimulation sim;
  sim.SetConfigFile("$ALICE_ROOT/FMD/Config.C");
  sim.ConvertRaw2SDigits("raw.root", 0);
}
