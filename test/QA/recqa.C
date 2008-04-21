void recqa()
{
  const char * kYear = "08" ; 
  AliGeomManager::LoadGeometry("geometry.root");	
  AliCDBManager * man = AliCDBManager::Instance();
  //man->SetDefaultStorage("alien://Folder=/alice/simulation/2007/PDC07_v4-09-Rev-00/Ideal/CDB/");
  man->SetDefaultStorage("local://$ALICE_ROOT");
  TString detectors("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  AliQA::SetQARefStorage(Form("%s%s/", AliQA::GetQARefDefaultStorage(), kYear)) ;
  AliQA::SetQARefDataDirName("Sim") ; //Data, Pedestals, BlackEvent, .....
   
  AliQADataMakerSteer qas ; 
  qas.Run(detectors.Data(), AliQA::kRECPOINTS);
  //qas.Reset() ;
  qas.Run(detectors.Data(), AliQA::kESDS, kTRUE);
  
}
