// To read PMD raw root data and fetch the adc value for each cell
void AliPMDRootDataRead(Int_t NEVT = 10)
{
  TObjArray pmdddlcont;

  gBenchmark->Start("");
  
  Bool_t junk;

  Int_t   xpad, ypad;
  Float_t xx, yy;

  AliPMDUtility cc;
  TH2F *h2 = new TH2F("h2","Y vs. X",200,-100.,100.,200,-100.,100.);


  for(Int_t ievt=3; ievt < NEVT; ievt++)
  {
  
  AliRawReaderRoot reader("08000033646024.10.root",ievt);
  // reader.NextEvent();
  cout<<" Processing Event No  : "<<ievt<<endl;
  /*
    reader.ReadHeader();
    cout << "LDC ID =       " << reader.GetLDCId()       << endl;
    cout << "Equipment ID = " << reader.GetEquipmentId() << endl;
    cout << "Data Size =    " << reader.GetDataSize()    << endl;
  */


  AliPMDRawStream stream(&reader);

  Int_t iddl = -1;
  while((iddl = stream.DdlData(&pmdddlcont)) >= 0)
  {
      Int_t ientries = pmdddlcont.GetEntries();

      //cout << "iddl = " << iddl << " ientries = " << ientries << endl;

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

//	  cout << iddl<<"  "<<row << " " << col << " " << sig << endl;


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
	  
	  if (det == 0)
	  {
	      cc.RectGeomCellPos(smn,xpad,ypad,xx,yy);
	      h2->Fill(xx,yy); 
	  }
	  
	}

      pmdddlcont.Clear();
    }

  }
  h2->Draw();

  gBenchmark->Show("");


}
