void CreateNoiseCut()
{

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");	
  AliPMDNoiseCut *ncut = new AliPMDNoiseCut();


  for(Int_t imod = 0; imod < 48; imod++)
    {
      Float_t cut = 4.;
      ncut->SetNoiseCut(imod,cut);
    }

  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Basanta Nandi");
  md->SetComment("Noise cut for different modules of PMD");

  AliCDBId id("PMD/Calib/NoiseCut",0,AliCDBRunRange::Infinity());

  man->GetDefaultStorage()->Put(ncut,id, md);

}
