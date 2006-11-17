//
void AliPMDRootDataRead()
{
  // To read PMD raw root data and fetch the adc value for each cell

  TObjArray pmdddlcont;

  Int_t ievt = 2;

  Bool_t junk;

  AliRawReaderRoot reader("raw_6b_61.root",ievt);
  //reader.NextEvent();
  //reader.NextEvent();



  /*
  reader.ReadHeader();
  cout << "LDC ID =       " << reader.GetLDCId()       << endl;
  cout << "Equipment ID = " << reader.GetEquipmentId() << endl;
  cout << "Data Size =    " << reader.GetDataSize()    << endl;
  */

  AliPMDRawStream stream(&reader);

  Int_t indexDDL = 0;


  for (Int_t iddl = 0; iddl < 6; iddl++)
    {
      
      reader.Select("PMD", iddl, iddl);

      junk = stream.DdlData(&pmdddlcont);

      Int_t ientries = pmdddlcont.GetEntries();
      for (Int_t ient = 0; ient < ientries; ient++)
	{
	  AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont.UncheckedAt(ient);
	  
	  Int_t det = pmdddl->GetDetector();
	  Int_t smn = pmdddl->GetSMN();
	  //Int_t mcm = pmdddl->GetMCM();
	  //Int_t chno = pmdddl->GetChannel();
	  Int_t row = pmdddl->GetRow();
	  Int_t col = pmdddl->GetColumn();
	  Int_t sig = pmdddl->GetSignal();

	  // cout << row << " " << col << " " << sig << endl;

	}
      pmdddlcont.Clear();

    }


}
