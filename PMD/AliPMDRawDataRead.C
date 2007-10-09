// To read PMD raw data

void AliPMDRawDataRead(Int_t iEvent)
{
  TObjArray pmdddlcont;
  
  for(Int_t ievt = 0; ievt < iEvent; ievt++)
    {
      AliRawReaderFile reader(ievt);
      AliPMDRawStream stream(&reader);
      
      Int_t indexDDL = 0;
      for (Int_t iddl = 0; iddl < 6; iddl++)
	{
	  
	  reader.Select("PMD", iddl, iddl);
	  //cout << reader.GetDataSize() << endl;
	  stream.DdlData(indexDDL,&pmdddlcont);
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
	      cout<<sig<<endl;
	      
	    }
	  pmdddlcont.Clear();
	  
	}

    }
}