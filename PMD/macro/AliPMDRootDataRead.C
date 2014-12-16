// To read PMD raw root data and fetch the adc value for each cell
void AliPMDRootDataRead(Char_t *file="rawfile.root",Int_t runNr = 0,
			Int_t NEVT = 10)
{
  TObjArray pmdddlcont;

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(runNr);

  gBenchmark->Start("");
  gStyle->SetOptStat(0);

  Int_t   xpad, ypad;
  Float_t xx, yy;

  AliPMDUtility cc;
  TH2F *h2 = new TH2F("h2","Y vs. X",200,-100.,100.,200,-100.,100.);

  TH1F *h1 = new TH1F("h1","",200,0.,200.);


  for(Int_t ievt=495; ievt < NEVT; ievt++)
    {
  
      //AliRawReaderRoot reader("09000095029017.10.root",ievt);
      AliRawReaderRoot reader(file,ievt);

      cout<<" Processing Event No  : "<<ievt<<endl;

      AliPMDRawStream stream(&reader);
      
      Int_t iddl = -1;
      
      while((iddl = stream.DdlData(&pmdddlcont)) >= 0)
	{

	  Int_t ientries = pmdddlcont.GetEntries();
	  
	  cout << "iddl = " << iddl << " ientries = " << ientries << endl;
	  
	  for (Int_t ient = 0; ient < ientries; ient++)
	    {
	      AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont.UncheckedAt(ient);
	  
	      Int_t det = pmdddl->GetDetector();
	      Int_t smn = pmdddl->GetSMN();
	      Int_t mcm = pmdddl->GetMCM();
	      Int_t pbus = pmdddl->GetPatchBusId();
	      //Int_t chno = pmdddl->GetChannel();
	      
	      Int_t row = pmdddl->GetRow();
	      Int_t col = pmdddl->GetColumn();
	      Int_t sig = pmdddl->GetSignal();

	      //if (iddl == 5)
	      //printf("%d %d %d %d %d %d\n",iddl,det,smn,pbus,row,col);
	      
	      if(mcm == 0) continue;
	      

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
		  h2->Fill(xx,yy,sig); 
		  h1->Fill(sig);
		}
	      
	    }
	  
      pmdddlcont.Delete();
	}
      
    }//event loop
  h2->Draw();
  cc.SetWriteModule(1);  // set here 0 not to print the module no
  cc.DrawPMDModule(0);   // 0 : for preshower, 1:for CPV plane
  gBenchmark->Show("");
}
