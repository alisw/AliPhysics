// To read PMD raw data

void AliPMDRawDataRead(Int_t iEvent)
{
  TObjArray pmdddlcont;
  
  TH2F *h2 = new TH2F("h2"," ",100,-100.,100.,100,-100.,100.);
  Float_t xx, yy;
  Int_t xpad, ypad;

  AliPMDUtility cc;

  for(Int_t ievt = 0; ievt < iEvent; ievt++)
    {
      AliRawReaderFile reader(ievt);
      AliPMDRawStream stream(&reader);
      
      Int_t iddl = -1;
      while ((iddl = stream.DdlData(&pmdddlcont)) >=0)
	{
	    //cout << " inside the macro DDLNO = " << iddl << endl; 
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
	      //cout<<sig<<endl;
	      if(smn <12)
	      {
		  xpad = col;
		  ypad = row;
	      }
	      else if(smn >=12 && smn < 24)
	      {
		  xpad = row;
		  ypad = col;
	      }
   
	      if (det == 1)
	      {
		  // Draw only for PRE plane
		  cc.RectGeomCellPos(smn,xpad,ypad,xx,yy);

		  h2->Fill(xx,yy);
	      }



	    }
	  pmdddlcont.Clear();
	  
	}

    }

  h2->Draw();
}
