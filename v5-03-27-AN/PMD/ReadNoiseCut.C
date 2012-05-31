void ReadNoiseCut()
{
  TFile * f = TFile::Open("$ALICE_ROOT/OCDB/PMD/Calib/NoiseCut/Run0_999999999_v0_s0.root");

  f->ls();
  if (!AliCDBEntry)
    {
      printf("Something is wrong ************ \n");
    }
  else if(AliCDBEntry)
    {
      AliCDBEntry->PrintId(); 
      AliCDBEntry->PrintMetaData();
      AliPMDNoiseCut * ncut = 0;
      ncut = (AliPMDNoiseCut*)AliCDBEntry->GetObject();
      ncut->Print(); 
    }
}
