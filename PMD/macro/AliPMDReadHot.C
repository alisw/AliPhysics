void AliPMDReadHot(const Char_t *filename) {
  TFile * f = TFile::Open(filename);
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
