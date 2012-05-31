void WriteOCDB(Int_t firstRun=0, Int_t lastRun=AliCDBRunRange::Infinity())
{
  //Write mapping to OCDB.
  //Run this macro from the same directory where mapping files *.data resides.
  //Author: Boris Polishchuk.

  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();

  AliCDBId id("PHOS/Calib/Mapping",firstRun,lastRun);
  AliCDBMetaData md;

  //Create mapping object.
  const char* mapFiles[20] = {
    "Mod0RCU0.data",
    "Mod0RCU1.data",
    "Mod0RCU2.data",
    "Mod0RCU3.data",
    "Mod1RCU0.data",
    "Mod1RCU1.data",
    "Mod1RCU2.data",
    "Mod1RCU3.data",
    "Mod2RCU0.data",
    "Mod2RCU1.data",
    "Mod2RCU2.data",
    "Mod2RCU3.data",
    "Mod3RCU0.data",
    "Mod3RCU1.data",
    "Mod3RCU2.data",
    "Mod3RCU3.data",
    "Mod4RCU0.data",
    "Mod4RCU1.data",
    "Mod4RCU2.data",
    "Mod4RCU3.data"
  };

  TString path = "./";
  
  path += "Mod";
  TString path2;
  TString path3;
  Int_t iMap = 0;
  AliAltroMapping* mapping;
  
  TObjArray objMap(20);
  
  for(Int_t iMod = 0; iMod < 5; iMod++) {
    path2 = path;
    path2 += iMod;
    path2 += "RCU";

    for(Int_t iRCU=0; iRCU<4; iRCU++) {
      path3 = path2;
      path3 += iRCU;
      path3 += ".data";
      AliAltroMapping* mapping = new AliCaloAltroMapping(path3.Data());
      objMap.AddAt(mapping,iMap);
      iMap++;
    }
  }

  //Put mapping object into OCDB: $ALICE_ROOT/OCDB/PHOS/Calib/Mapping/
  storage->Put(&objMap,id,&md);
}
