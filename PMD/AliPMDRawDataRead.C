//
void AliPMDRawDataRead()
{
  // To read PMD raw data

  Int_t ievt = 0;
  AliRawReaderFile reader(ievt);
  AliPMDRawStream stream(&reader);

  //  FILE *fpw2 = fopen("rawread2.txt","w");

  while(stream.Next())
    {
      printf("%d %d %d %d\n",stream.GetModule(),stream.GetRow(),
	     stream.GetColumn(),stream.GetSignal());

      Int_t det = stream.GetDetector();
      Int_t smn = stream.GetSMN();
      Int_t mcm = stream.GetMCM();
      Int_t chno = stream.GetChannel();
      Int_t row = stream.GetRow();
      Int_t col = stream.GetColumn();
      Int_t sig = stream.GetSignal();

      //  fprintf(fpw2,"%d %d %d %d %d %d %d\n",det, smn, row, col,
      //      mcm, chno, sig);

    }

}
