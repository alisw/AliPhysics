// To read PMD raw root data and fetch the adc value for each cell
void AliPMDRootDataRead(Int_t NEVT = 10)
{
  TObjArray pmdddlcont;

  gBenchmark->Start("");
  gStyle->SetOptStat(0);

  TFile *pedfile = new TFile("PMD_PED.root");

  if(!pedfile)
    {
      printf("ERROR --- NO PEDESTAL (PMD_PED.root) FILE IS FOUND IN THE CURRENT DIRECTORY--- STOP GAIN DA\n");
      return -3;
    }


  Float_t fPedMeanRMS[2][24][48][96];

  for(Int_t i = 0; i < 2; i++)
  {
      for(Int_t j = 0; j < 24; j++)
      {
	  for(Int_t k = 0; k < 48; k++)
	  {
	      for(Int_t l = 0; l < 96; l++)
	      {
		  fPedMeanRMS[i][j][k][l] = 0.;
	      }
	  }
      }
  }

  Int_t det, sm, row, col;
  Float_t mean, rms;

  TTree *ped =(TTree*)pedfile->Get("ped");

  ped->SetBranchAddress("det",&det);
  ped->SetBranchAddress("sm",&sm);
  ped->SetBranchAddress("row",&row);
  ped->SetBranchAddress("col",&col);
  ped->SetBranchAddress("mean",&mean);
  ped->SetBranchAddress("rms",&rms);

  Int_t nentries = (Int_t)ped->GetEntries();

  for (Int_t ient = 0; ient < nentries; ient++)
  {
      ped->GetEntry(ient);
      fPedMeanRMS[det][sm][row][col] = mean + 10.*rms;
      //printf("Mean= %f, RMS= %f, PedMeanRMS=%f\n",mean,rms,fPedMeanRMS[det][sm][row][col]);

    }


  pedfile->Close();
  delete pedfile;
  pedfile = 0x0;



  Int_t   xpad, ypad;
  Float_t xx, yy;

  AliPMDUtility cc;
  TH2F *h2 = new TH2F("h2","Y vs. X",200,-100.,100.,200,-100.,100.);

  TH1F *h1 = new TH1F("h1","",200,0.,200.);


  for(Int_t ievt=3; ievt < NEVT; ievt++)
  {
  
  AliRawReaderRoot reader("08000042184000.10.root",ievt);

  cout<<" Processing Event No  : "<<ievt<<endl;

  AliPMDRawStream stream(&reader);

  Int_t iddl = -1;

  while((iddl = stream.DdlData(&pmdddlcont)) >= 0)
  {

      Int_t ientries = pmdddlcont.GetEntries();

      //cout << "iddl = " << iddl << " ientries = " << ientries << endl;

      for (Int_t ient = 0; ient < ientries; ient++)
	{
	  AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont.UncheckedAt(ient);
	  
	  det = pmdddl->GetDetector();
	  Int_t smn = pmdddl->GetSMN();
	  Int_t mcm = pmdddl->GetMCM();
	  Int_t pbus = pmdddl->GetPatchBusId();
	  //Int_t chno = pmdddl->GetChannel();

	  row = pmdddl->GetRow();
	  col = pmdddl->GetColumn();
	  Int_t sig = pmdddl->GetSignal();

	  if(mcm == 0) continue;

	  Int_t sig1 = sig - (Int_t)fPedMeanRMS[det][smn][row][col];


	  if (sig1 > 0)
	      sig = sig1;
	  else
	      sig = 0;

	  if (ievt == 159 && sig > 0) cout << row << " " << col << " " << sig << endl;


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


  h2->Draw();

  

  }



  cc.SetWriteModule(1);
  cc.DrawPMDModule();
  gBenchmark->Show("");


}
