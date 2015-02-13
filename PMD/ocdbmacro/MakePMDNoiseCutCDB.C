/**************************************************************************
 * To Create PMD Noise map into the OCDB object
 * sjena@cern.ch
 * Mon Nov 22 19:54:27 CET 2010
 *                     
 **************************************************************************/

void MakePMDNoiseCutCDB() {

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");	
  AliPMDNoiseCut *ncut = new AliPMDNoiseCut();


  for(Int_t imod = 0; imod < 48; imod++)
    {
      Float_t cut = 4.;
      ncut->SetNoiseCut(imod,cut);
    }

  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Satyajit Jena");
  md->SetComment("Noise cut for different modules of PMD");

  AliCDBId id("PMD/Calib/NoiseCut",0,AliCDBRunRange::Infinity());

  man->GetDefaultStorage()->Put(ncut,id, md);

}
