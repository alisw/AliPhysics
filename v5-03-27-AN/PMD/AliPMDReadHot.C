void AliPMDReadHot()
{
  TFile * f = TFile::Open("Run0_0_v0_s3.root");
  f->ls();
  if (!AliCDBEntry)
    {
      printf("Something is wrong ************ \n");
    }
  else if(AliCDBEntry)
    {
      AliCDBEntry->PrintId(); 
      AliCDBEntry->PrintMetaData();
      AliPMDHotData * hot = 0;
      hot = (AliPMDHotData*)AliCDBEntry->GetObject();
      hot->Print(); 
    }
}
