//
// script with functions to use AliJetSpectrumUnfolding class
//


void loadlibs(){
  // load all the needed libs to run wit root only

  gSystem->Load("libTree");
  gSystem->Load("libPhysics");
  gSystem->Load("libHist");
  gSystem->Load("libVMC");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libAOD");
  gSystem->Load("libESD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libJETAN");
  gSystem->Load("libPWG4JetTasks");


}



void draw(const char* fileName = "unfolded_pwg4spec.root", const char* folder = "unfolding", Bool_t proj = kTRUE)
{

  loadlibs();

  AliJetSpectrumUnfolding* jetSpec = new AliJetSpectrumUnfolding(folder, folder);

  TFile::Open(fileName);
  jetSpec->LoadHistograms();
  
  
  if (proj)
  {
    canvas1 = new TCanvas("Response Map Projection", "Response Map Projection", 500, 500);
    canvas1->Divide(2);
  
    Int_t style = 1;
    const Int_t NRGBs = 5;
    const Int_t NCont = 500;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    gStyle->SetPalette(style);

    canvas1->cd(1);
    gPad->SetLogz();
    h2 = (jetSpec->GetCorrelation())->Projection(2,3);
    h2->SetXTitle("z^{lp}_{gen}");
    h2->SetYTitle("z^{lp}_{rec}");    
    h2->DrawCopy("colz");
  
    canvas1->cd(2);
    gPad->SetLogz();  
    h3 = (jetSpec->GetCorrelation())->Projection(0,1);
    h3->SetXTitle("E^{jet}_{gen} [GeV]");
    h3->SetYTitle("E^{jet}_{rec} [GeV]");    
    h3->DrawCopy("colz");
  }
  jetSpec->DrawComparison("Draw_unfolded_pwg4spec", jetSpec->GetGenSpectrum());

  return;
}

//________________________________________________________________________________________________________________
void unfold(const char* fileNameGen = "gen_pwg4spec.root", const char* folder = "unfolding", const char* fileNameRec = "rec_pwg4spec.root", const char* fileNameUnf = "unfolded_pwg4spec.root")
{
  // function to load jet spectra from the output file of the task AliAnalysisTaskJetSpectrum
  // to do the unfolding

  loadlibs();

  AliJetSpectrumUnfolding* jetSpec = new AliJetSpectrumUnfolding(folder, folder);

  TFile::Open(fileNameGen);
  jetSpec->LoadHistograms();

  TFile::Open(fileNameRec);
  TH2F* hist = (TH2F*) gFile->Get("unfolding/fRecSpectrum");
  jetSpec->SetRecSpectrum(hist);
  hist = (TH2F*) gFile->Get("unfolding/fGenSpectrum");
  if(hist->GetEntries()>0)jetSpec->SetGenSpectrum(hist);

  jetSpec->ApplyBayesianMethod(0.3, 20, 0, 0);
  // last parameter = calculateErrors  <- this method to calculate the errors takes a lot of time
   
  TFile* file = TFile::Open(fileNameUnf,"RECREATE");
  jetSpec->SaveHistograms();
  file->Close();
}

//___________________________________________________________________________
void buildSpectra(Int_t caseNo, const char* inFile, const char* outFile)
{
  // function to test 2D Bayesian unfolding with other spectra
  // build from a function

  loadlibs();



  AliJetSpectrumUnfolding* jetSpec = new AliJetSpectrumUnfolding("unfolding", "unfolding");

  TFile::Open(inFile);
  jetSpec->LoadHistograms("unfolding");
  
  TCanvas *c1 = new TCanvas();TH2 *h2 = (jetSpec->GetCorrelation())->Projection(0,3);
  h2->DrawCopy("colz");
  c1->Update();
  if(getchar()=='q')return;
  

  switch (caseNo)
  {
    case 1: func = new TF2("func", "501-x-y"); break;
    case 2: func = new TF2("func", "1000 * 1/(y+x+1)"); break;
    case 3: func = new TF2("func", "1./(y*pow(x,5.7))"); break;
    case 4: func = new TF2("func", "exp(-0.1*x - 0.1*y)"); break;
    case 5: func = new TF2("func", "x*x + (y**3)/100."); break;
    case 6: func = new TF2("func", "1000*y*exp(-0.1*x)"); break;
    case 7: func = new TF2("func", "exp(-((x-100.)/(0.3*x))**2 - ((y-0.6)/(0.8*y))**2)"); break;
    default: return;
  }

  //new TCanvas; func->Draw();

  jetSpec->SetGenRecFromFunc(func);
              
  TFile* file = TFile::Open(outFile,"RECREATE");
  jetSpec->SaveHistograms();
  
  h2 = (jetSpec->GetCorrelation())->Projection(0,3);
  h2->DrawCopy("colz");
  c1->Update();
  file->Close();

  //new TCanvas; jetSpec->GetRecSpectrum()->DrawCopy("colz");
  //new TCanvas; jetSpec->GetGenSpectrum()->DrawCopy("colz");
}

//_____________________________________________________________________________________________
void buildResponseMap(const char* fileName = "responseMap.root")
{
  // function to build a Response Map with a gaussian distribution
  loadlibs();

  AliJetSpectrumUnfolding* jetSpec = new AliJetSpectrumUnfolding("unfolding", "unfolding");

  TF2* func = new TF2("func", "exp(-((x-[0])/[1])**2 - ((y-[2])/[3])**2)");
  
  bool bPerfect = false;
  
  Double_t var[4];
  Float_t sigmax, sigmay;
  for (Int_t tx=1; tx<=jetSpec->GetCorrelation()->GetAxis(0)->GetNbins(); tx++)
    for (Int_t ty=1; ty<=jetSpec->GetCorrelation()->GetAxis(2)->GetNbins(); ty++)
    {
      var[0] = jetSpec->GetCorrelation()->GetAxis(0)->GetBinCenter(tx);
      var[2] = jetSpec->GetCorrelation()->GetAxis(2)->GetBinCenter(ty);
      sigmax = 0.2*var[0];
      sigmay = 0.2*var[2];
      func->SetParameters(var[0],sigmax,var[2],sigmay);
      for (Int_t mx=1; mx<=jetSpec->GetCorrelation()->GetAxis(1)->GetNbins(); mx++)
        for (Int_t my=1; my<=jetSpec->GetCorrelation()->GetAxis(3)->GetNbins(); my++)
        {
          var[1] = jetSpec->GetCorrelation()->GetAxis(1)->GetBinCenter(mx);
	  var[3] = jetSpec->GetCorrelation()->GetAxis(3)->GetBinCenter(my);


	  if(bPerfect){
	    if (var[1]==var[0] && var[3]==var[2])
	      jetSpec->GetCorrelation()->Fill(var,1);
	  }
	  else {
	    // cut at  sigma
	    if (TMath::Abs(var[1]-var[0]) < sigmax || TMath::Abs(var[3]-var[2]) < sigmay)
	      jetSpec->GetCorrelation()->Fill(var,func->Eval(var[1],var[3]));
	  }
        }
    }


  TFile* file = TFile::Open(fileName,"RECREATE");
  jetSpec->SaveHistograms();
  file->Close();
}

//_____________________________________________________________________________
void StatisticalUncertainties(const char* fileNameGen = "gen_pwg4spec.root", const char* folder = "unfolding", const char* fileNameRec = "rec_pwg4spec.root", const char* fileNameOut = "Uncertainties2DBayesUnf.root")
{
  // This function gives the statistical uncertainties due to the 2D Bayesian Unfolding

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libJETAN");
  gSystem->Load("libPWG4JetTasks");

  AliJetSpectrumUnfolding* jetSpec = new AliJetSpectrumUnfolding(folder, folder);

  TFile::Open(fileNameRec);
  jetSpec->LoadHistograms();

  TFile::Open(fileNameGen);
  TH2F* hist = (TH2F*) gFile->Get("unfolding/fGenSpectrum");
  jetSpec->SetGenSpectrum(hist);
    
  // create sigma histogram
  TH2F* sigma = (TH2F*)(jetSpec->GetGenSpectrum())->Clone("sigma");
  sigma->Reset();
    
  THnSparseF* hcorr = (THnSparseF*)(jetSpec->GetCorrelation())->Clone();
  TH2F*       htrue = (TH2F*)(jetSpec->GetGenSpectrum())->Clone();  
  TH2F*       hmeas = (TH2F*)(jetSpec->GetRecSpectrum())->Clone();    
  TH2F*       hunfo = (TH2F*)(jetSpec->GetUnfSpectrum())->Clone();      
  
  Int_t nIterations = 1000;
  Float_t binContent;
  for(Int_t i=0; i<nIterations; i++)
  {
    printf("iteration = %d\n", i);
    // reset histograms
    jetSpec->GetRecSpectrum()->Reset();
    jetSpec->GetGenSpectrum()->Reset();
    jetSpec->GetCorrelation()->Reset();
    jetSpec->GetUnfSpectrum()->Reset();
  
    THnSparseF* tmpcorr = (THnSparseF*)hcorr->Clone("tmpcorr"); 
    TH2F*       tmptrue = (TH2F*)htrue->Clone("tmptrue");  
    
    jetSpec->SetGenSpectrum(tmptrue);
    jetSpec->SetCorrelation(tmpcorr);
  
    // randomize reconstructed distribution
    for (Int_t me=1; me<=hmeas->GetNbinsX(); me++)
      for (Int_t mz=1; mz<=hmeas->GetNbinsY(); mz++)
      {
        binContent = hmeas->GetBinContent(me,mz);
        if (binContent>0)
        {
          TF1* poisson = new TF1("poisson", "TMath::Poisson(x,[0])",binContent*0.25, binContent*1.25);
          poisson->SetParameters(binContent,0.);
          binContent = poisson->GetRandom();
          delete poisson;
        }  
        jetSpec->GetRecSpectrum()->SetBinContent(me,mz, binContent);
      } 
        
    // unfold
    jetSpec->ApplyBayesianMethod(0.2, 20, 0, 0);
    
    // calculate sigma^2
    for (Int_t te=1; te<=sigma->GetNbinsX(); te++)
      for (Int_t tz=1; tz<=sigma->GetNbinsY(); tz++)
      {
        if (htrue->GetBinContent(te,tz)!=0)
        {
          binContent = (jetSpec->GetUnfSpectrum()->GetBinContent(te,tz) -
                        htrue->GetBinContent(te,tz))/htrue->GetBinContent(te,tz);
          binContent *= binContent;
          sigma->SetBinContent(te,tz, binContent + sigma->GetBinContent(te,tz));
        }  
      } 
  }
 
  // calculate sigma   
  for (Int_t te=1; te<=sigma->GetNbinsX(); te++)
    for (Int_t tz=1; tz<=sigma->GetNbinsY(); tz++)
    {
      binContent = sigma->GetBinContent(te,tz);
      binContent = TMath::Sqrt(binContent/(Float_t)nIterations);
      sigma->SetBinContent(te,tz, binContent);
    } 
          
  const Int_t NRGBs = 5;
  const Int_t NCont = 500;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);   
  gStyle->SetPalette(1);  

  new TCanvas();
  gPad->SetLogz();
  sigma->SetTitle("#sigma((U_{R} - U)/U)");
  sigma->SetXTitle("E^{jet} [GeV]");
  sigma->SetYTitle("z^{lp}");  
  sigma->DrawCopy();
  
  TFile* file = TFile::Open(fileNameOut,"RECREATE");
  sigma->Write();
  file->Close();   
}

//_______________________________________________________________________________________________________________
void FillSpecFromFiles(const char* fileNameReal = "histos_pwg4spec.root",const char* fileNameSim = "histos_pwg4spec.root")
{
  // This functions avoids problems when the number of bins or the bin limits
  // in the histograms of the AliJetSpectrumUnfolding and AliAnalysisTaskJetSpectrum classes
  // are different.

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libJETAN");
  gSystem->Load("libPWG4JetTasks");

  file = new TFile(fileNameSim,"read");
  tlist = dynamic_cast<TList*> (file->Get("pwg4spec"));

  THnSparseF *fhCorrelation  = 0;
  for(int ic = 0;ic<3;++ic){
    THnSparseF *hTmp = (THnSparseF*)(tlist->FindObject(Form("fhnCorrelation_%d",ic)));    
    if(fhCorrelation==0)fhCorrelation =  (THnSparseF*)hTmp->Clone("fhnCorrelation");
    else fhCorrelation->Add(hTmp);
  }
  TH2F *fhEGenZGen = (TH2F*)(tlist->FindObject("fh2EGenZGen"));

  AliJetSpectrumUnfolding *jetSpec = new AliJetSpectrumUnfolding("unfolding","unfolding");

  // generated jets (true distribution)
  for (Int_t te=1; te<=fhEGenZGen->GetNbinsX(); te++)
    for (Int_t tz=1; tz<=fhEGenZGen->GetNbinsY(); tz++)
    {
       Float_t ej = fhEGenZGen->GetXaxis()->GetBinCenter(te);
       Float_t  z = fhEGenZGen->GetYaxis()->GetBinCenter(tz);
       Int_t bine = jetSpec->GetGenSpectrum()->GetXaxis()->FindBin(ej);
       Int_t binz = jetSpec->GetGenSpectrum()->GetYaxis()->FindBin(z);
       Float_t cont = jetSpec->GetGenSpectrum()->GetBinContent(bine,binz);
       Float_t err  = jetSpec->GetGenSpectrum()->GetBinError(bine,binz);
       jetSpec->GetGenSpectrum()->SetBinContent(bine, binz, cont + fhEGenZGen->GetBinContent(te, tz));
       jetSpec->GetGenSpectrum()->SetBinError(bine, binz, err + fhEGenZGen->GetBinError(te, tz));
    }


  Printf("Bins %.0E",jetSpec->GetCorrelation()->GetNbins());
  for (Int_t idx=1; idx<=fhCorrelation->GetNbins(); idx++)
  {
    //printf("idx: %d\n",idx);
    Double_t var[4];
    Int_t bin[4];
    Float_t BinContent = fhCorrelation->GetBinContent(idx, bin);
    var[0] = fhCorrelation->GetAxis(0)->GetBinCenter(bin[0]);
    var[1] = fhCorrelation->GetAxis(1)->GetBinCenter(bin[1]);
    var[2] = fhCorrelation->GetAxis(2)->GetBinCenter(bin[2]);
    var[3] = fhCorrelation->GetAxis(3)->GetBinCenter(bin[3]);
    bin[0] = jetSpec->GetCorrelation()->GetAxis(0)->FindBin(var[0]);
    bin[1] = jetSpec->GetCorrelation()->GetAxis(1)->FindBin(var[1]);    
    bin[2] = jetSpec->GetCorrelation()->GetAxis(2)->FindBin(var[2]);
    bin[3] = jetSpec->GetCorrelation()->GetAxis(3)->FindBin(var[3]);        
    Float_t cont = jetSpec->GetCorrelation()->GetBinContent(bin);
    Float_t err  = jetSpec->GetCorrelation()->GetBinError(bin);
    jetSpec->GetCorrelation()->SetBinContent(bin, cont + fhCorrelation->GetBinContent(idx));
    jetSpec->GetCorrelation()->SetBinError(bin, err + fhCorrelation->GetBinError(idx));
  }
  // need number of entries for correct drawing
  jetSpec->GetCorrelation()->SetEntries(fhCorrelation->GetEntries());


  file = TFile::Open("gen_pwg4spec.root", "RECREATE");
  jetSpec->SaveHistograms();
  Printf("Bins %.0E",jetSpec->GetCorrelation()->GetNbins());
  file->Close();
  jetSpec->GetGenSpectrum()->Reset();
  printf("true distribution has been set\n");

  // reconstructed jets (measured distribution)


  // Response function
  jetSpec->GetCorrelation()->Reset();
  jetSpec->GetCorrelation()->Sumw2();
  jetSpec->GetGenSpectrum()->Reset();
  jetSpec->GetRecSpectrum()->Reset();


  file = new TFile(fileNameReal,"read");
  tlist = dynamic_cast<TList*> (file->Get("pwg4spec"));

  
  TH2F *fhERecZRec = (TH2F*)(tlist->FindObject("fh2ERecZRec"));  

  for (Int_t me=1; me<=fhERecZRec->GetNbinsX(); me++)
    for (Int_t mz=1; mz<=fhERecZRec->GetNbinsY(); mz++)
    {
       Float_t erec = fhERecZRec->GetXaxis()->GetBinCenter(me);
       Float_t   zm = fhERecZRec->GetYaxis()->GetBinCenter(mz);
       Int_t bine   = jetSpec->GetRecSpectrum()->GetXaxis()->FindBin(erec);
       Int_t binz   = jetSpec->GetRecSpectrum()->GetYaxis()->FindBin(zm);
       Float_t cont = jetSpec->GetRecSpectrum()->GetBinContent(bine, binz);
       Float_t err  = jetSpec->GetRecSpectrum()->GetBinError(bine, binz);
       jetSpec->GetRecSpectrum()->SetBinContent(bine, binz, cont + fhERecZRec->GetBinContent(me, mz));
       jetSpec->GetRecSpectrum()->SetBinError(bine, binz, err + fhERecZRec->GetBinError(me, mz));
    }


  // for control again, but now from the rec file
  // generated jets (true distribution)
  TH2F *fhEGenZGen = (TH2F*)(tlist->FindObject("fh2EGenZGen"));  
  for (Int_t te=1; te<=fhEGenZGen->GetNbinsX(); te++)
    for (Int_t tz=1; tz<=fhEGenZGen->GetNbinsY(); tz++)
    {
       Float_t ej = fhEGenZGen->GetXaxis()->GetBinCenter(te);
       Float_t  z = fhEGenZGen->GetYaxis()->GetBinCenter(tz);
       Int_t bine = jetSpec->GetGenSpectrum()->GetXaxis()->FindBin(ej);
       Int_t binz = jetSpec->GetGenSpectrum()->GetYaxis()->FindBin(z);
       Float_t cont = jetSpec->GetGenSpectrum()->GetBinContent(bine,binz);
       Float_t err  = jetSpec->GetGenSpectrum()->GetBinError(bine,binz);
       jetSpec->GetGenSpectrum()->SetBinContent(bine, binz, cont + fhEGenZGen->GetBinContent(te, tz));
       jetSpec->GetGenSpectrum()->SetBinError(bine, binz, err + fhEGenZGen->GetBinError(te, tz));
    }


  file = TFile::Open("rec_pwg4spec.root", "RECREATE");
  jetSpec->SaveHistograms();
  file->Close();

  printf("reconstructed distribution has been set\n");    
  printf("response map has been set\n");
  
}

void correct(){
  // simple steering to correct a given distribution;
  loadlibs();
  // rec and gen
  //  FillSpecFromFiles("pwg4spec_15-50_all.root","pwg4spec_allpt.root");
  FillSpecFromFiles("pwg4spec_allpt.root","pwg4spec_allpt.root");

  char name[100];
  sprintf(name, "unfolded_pwg4spec.root");

  unfold("gen_pwg4spec.root", "unfolding", "rec_pwg4spec.root", name);
  //draw(name, "unfolding", 1); 

}

void mergeJetAnaOutput(){
  // This is used to merge the analysis-output from different 
  // data samples/pt_hard bins
  // in case the eventweigth was set to xsection/ntrials already, this
  // is not needed. Both methods only work in case we do not mix different 
  // pt_hard bins, and do not have overlapping bins

  const Int_t nBins = 2;
  // LHC08q jetjet100: Mean = 1.42483e-03, RMS = 6.642e-05
  // LHC08r jetjet50: Mean = 2.44068e-02, RMS = 1.144e-03
  // LHC08v jetjet15-50: Mean = 2.168291 , RMS = 7.119e-02
  // const Float_t xsection[nBins] = {2.168291,2.44068e-02};
  Float_t xsection[nBins] = {0,0};
  Float_t nTrials[nBins] = {0,0};
  Float_t sf[nBins] = {0,0};

  const char *cFile[nBins] = {"pwg4spec_0000000-1000000_LHC08v_jetjet15-50.root",
			      "pwg4spec_0000000-1000000_LHC08r_jetjet50.root"};


  TList *lIn[nBins];
  TFile *fIn[nBins];
  for(int ib = 0;ib < nBins;++ib){
    fIn[ib] = TFile::Open(cFile[ib]);
    lIn[ib] = (TList*)fIn[ib]->Get("pwg4spec");
    TH1* hTrials = (TH1F*)lIn[ib]->FindObject("fh1PtHard_Trials");
    TProfile* hXsec = (TProfile*)lIn[ib]->FindObject("fh1Xsec");
    xsection[ib] = hXsec->GetBinContent(1);
    nTrials[ib] = hTrials->Integral();
    sf[ib] = xsection[ib]/ nTrials[ib];
  }

  TFile *fOut = new TFile("pwg4spec_allpt.root","RECREATE");
  TList *lOut = new TList();
  lOut->SetName(lIn[0]->GetName());
  // for the start scale all...
  for(int ie = 0; ie < lIn[0]->GetEntries();++ie){
    TH1 *h1Add = 0;
    THnSparse *hnAdd = 0;
    Printf("%d: %s",ie, lIn[0]->At(ie)->GetName());
    for(int ib = 0;ib < nBins;++ib){
      // dynamic cast does not work with cint
      TObject *h = lIn[ib]->At(ie);
      if(h->InheritsFrom("TH1")){
	Printf("ib %d",ib);
	TH1 *h1 = (TH1*)h;
	if(ib==0){
	  h1Add = (TH1*)h1->Clone(h1->GetName());
	  h1Add->Scale(sf[ib]);
	}
	else{
	  h1Add->Add(h1,sf[ib]);
	}
      }
      else if(h->InheritsFrom("THnSparse")){
	Printf("ib %d",ib);
	THnSparse *hn = (THnSparse*)h;
	if(ib==0){
	  hnAdd = (THnSparse*)hn->Clone(hn->GetName());
	  hnAdd->Scale(sf[ib]);
	}
	else{
	  hnAdd->Add(hn,sf[ib]);
	}
      }
      

    }// ib
    if(h1Add)lOut->Add(h1Add);
    else if(hnAdd)lOut->Add(hnAdd);
  }
  fOut->cd();
  lOut->Write(lOut->GetName(),TObject::kSingleKey);
  fOut->Close();



}
