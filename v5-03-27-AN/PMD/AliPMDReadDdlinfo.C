void AliPMDReadDdlinfo()
{
  TFile * f = TFile::Open("/Users/basanta/ALISOFT/PMD/VarInit/OCDB/PMD/Calib/Ddlinfo/Run0_999999999_v0_s0.root");
  //  TFile * f = TFile::Open("$ALICE_ROOT/OCDB/PMD/Calib/Mapping/Run0_999999999_v0_s0.root");

  f->ls();
  if (!AliCDBEntry)
    {
      printf("Something is wrong ************ \n");
    }
  else if(AliCDBEntry)
    {
      AliCDBEntry->PrintId(); 
      AliCDBEntry->PrintMetaData();
      AliPMDddlinfoData * mapdata = 0;
      mapdata = (AliPMDddlinfoData*)AliCDBEntry->GetObject();
      mapdata->Print(); 
    }
}
