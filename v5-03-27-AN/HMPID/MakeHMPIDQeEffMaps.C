void MakeHMPIDQeEffMaps()
{
  //
  // Create HMPID Measured Quantum Efficiency Maps in OCDB
  // QE measurement/scan done by Antonello 
  // Data are extracted from the "Photo Cathode Production excel files
  // Graphs contain: x and y position on PC surface and the normalized current
  //
  // Position of the photo cathodes (real name) in the HMPID modules
  // RICH0 , PC0-real name PC77 ||   PC1-real name PC74 || PC2-real name PC79 || PC3-real name PC70 || PC4-real name PC72 || PC5-real name PC48 
  // RICH1 , PC0-real name PC61 ||   PC1-real name PC65 || PC2-real name PC55 || PC3-real name PC54 || PC4-real name PC73 || PC5-real name PC62
  // RICH2 , PC0-real name PC60 ||   PC1-real name PC37 || PC2-real name PC59 || PC3-real name PC56 || PC4-real name PC38 || PC5-real name PC40
  // RICH3 , PC0-real name PC42 ||   PC1-real name PC41 || PC2-real name PC44 || PC3-real name PC43 || PC4-real name PC46 || PC5-real name PC45 
  // RICH4 , PC0-real name PC57v2 || PC1-real name PC66 || PC2-real name PC67 || PC3-real name PC68 || PC4-real name PC64 || PC5-real name PC63 
  // RICH5 , PC0-real name PC53 ||   PC1-real name PC47 || PC2-real name PC51 || PC3-real name PC49 || PC4-real name PC52 || PC5-real name PC50 
  // RICH6 , PC0-real name PC71 ||   PC1-real name PC84 || PC2-real name PC82 || PC3-real name PC83 || PC4-real name PC75 || PC5-real name PC81 
  //
 
  TGraph2D *hmpQeEffGraph[7][6];

  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  TObjArray *hmpQeEffMapArray = new TObjArray();

  for(Int_t ich=0;ich<7;ich++){
   for(Int_t ipc=0;ipc<6;ipc++){
      hmpQeEffGraph[ich][ipc]=new TGraph2D(Form("HMPIDQeMapping/HmpidQeMapMod%dPc%d.txt",ich,ipc),"%lg %lg %lg","");
      hmpQeEffGraph[ich][ipc]->SetName(Form("HmpidQeMapMod%dPc%d",ich,ipc));
      hmpQeEffMapArray->AddLast(hmpQeEffGraph[ich][ipc]);
      }
   }
 
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Domenico DiBari");
  md->SetComment("Quantum Efficiany Maps for HMPID simulation");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("HMPID/Calib/QeMap",0,AliCDBRunRange::Infinity());
  man->GetDefaultStorage()->Put(hmpQeEffMapArray,id, md);
  
  return;
}
