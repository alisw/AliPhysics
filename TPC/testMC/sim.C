void sim(char * configPath, char * tpcDBpath, Int_t nevents){
  //
  //
  //
  //AliCDBManager * man = AliCDBManager::Instance();
  //man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  //man->SetRun(0);
  //man->SetSpecificStorage("TPC/*/*",tpcDBpath);

  AliSimulation sim;
  sim.SetConfigFile(configPath);
  //sim.SetSpecificStorage("TPC/*/*",tpcDBpath);
  sim.Run(nevents);


}
