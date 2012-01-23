void FitCDFLocal() {

  ///////////////////////////////////////////////////////////////////
  //
  // Example macro to read local N-Tuples of JPSI 
  // and bkg candidates and perform log-likelihood 
  // minimization using the minimization handler class
  //
  // Origin: C. Di Giglio
  //
  ///////////////////////////////////////////////////////////////////

  Bool_t useParFiles=kFALSE;
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/LoadLibraries.C");
  LoadLibraries(useParFiles);

  TH1F* hCsiMCPithya = new TH1F();
  TH1F* hCsiMCevtgen = new TH1F();

  TString datadir=".";
  TString useFile = datadir.Data();
  useFile.Append("/CmpBdecayKine.root"); // file with MC templates for J/psi (<-B) X distribution (in different pT bins)
  TFile *fCmpBdecayKine = new TFile(useFile);
  hCsiMCPithya = (TH1F*)fCmpBdecayKine->Get("PsproperPithya46GeV"); // Pithya case: 4 - 6 GeV
  hCsiMCPithya->Scale(1./hCsiMCPithya->Integral());
  hCsiMCevtgen = (TH1F*)fCmpBdecayKine->Get("PsproperEvtGen46GeV"); // EvtGen case: 4 - 6 GeV
  hCsiMCevtgen->Scale(1./hCsiMCevtgen->Integral());

  Double_t* x=0x0; Double_t* m=0x0; Int_t n=0;
  AliAnalysisBtoJPSItoEle* aBtoJPSItoEle =new AliAnalysisBtoJPSItoEle();
  Double_t paramInputValues [16] = /*  
                                     paramInputValues[0]   ----> fRadius;
                                     paramInputValues[1]   ----> fTheta;
                                     paramInputValues[2]   ----> fPhi;
                                     paramInputValues[3]   ----> fOneOvLamPlus;
                                     paramInputValues[4]   ----> fOneOvLamMinus;
                                     paramInputValues[5]   ----> fOneOvLamSym;
                                     paramInputValues[6]   ----> fMSlope;
                                     paramInputValues[7]   ----> fB;
                                     paramInputValues[8]   ----> fFsig;
                                     paramInputValues[9]   ----> fMmean;
                                     paramInputValues[10]  ----> fNexp;
                                     paramInputValues[11]  ----> fSigma;
                                     paramInputValues[12]  ----> fAlpha;
                                     paramInputValues[13]  ----> fNorm;
                                     paramInputValues[14]  ----> fSigmaResol;
                                     paramInputValues[15]  ----> fNResol;
                                   */
                                   { 0.5,TMath::Pi()/4.,
                                     TMath::Pi()/4.,
                                     0.00210,
                                     0.00480,
                                     0.,
                                     0.,
                                     0.227272,
                                     0.3981,
                                     3.09568,
                                     0.9671,
                                     0.0239,
                                     0.6599,
                                     0.04587,
                                     55.4, 
                                     0.033 }; // Starting values for parameters 

  TFile *f= new TFile("CdfFit_OneYear.root"); // file with N-tuples for one year collected statistics (in different pT bins)
  //TNtuple *nt=(TNtuple*)f->Get("fNtupleJPSI12GeV");
  //TNtuple *nt=(TNtuple*)f->Get("fNtupleJPSI23GeV");
  //TNtuple *nt=(TNtuple*)f->Get("fNtupleJPSI34GeV");
  TNtuple *nt=(TNtuple*)f->Get("fNtupleJPSI46GeV"); // N-Tuples in 4 - 6 Gev
  nt->ls();

  aBtoJPSItoEle->ReadCandidates(nt,x,m,n); // read N-Tuples
  printf("+++\n+++  
          Number of total candidates (prim J/psi + secondary J/psi + bkg) ---> %d candidates 
          \n+++\n",n);

  aBtoJPSItoEle->SetFitHandler(x,m,n); // Set the fit handler with given values of x, m, # of candidates 

  aBtoJPSItoEle->CloneMCtemplate(hCsiMCPithya);    // clone MC template and copy internally
                                                   // in this way any model can be setted from outside
  //aBtoJPSItoEle->CloneMCtemplate(hCsiMCEvtGen);  // clone MC template and copy internally
                                                   // in this way any model can be setted from outside

  aBtoJPSItoEle->SetCsiMC(); // Pass the MC template to the CDF fit function

  AliBtoJPSItoEleCDFfitHandler* fithandler = aBtoJPSItoEle->GetCDFFitHandler(); // Get the fit handler

  //
  // Set some fit options through the handler class
  //
  fithandler->SetErrorDef(0.5); // tells Minuit that the error interval is the one in which
                                // the function differs from the minimum for less than setted value

  fithandler->SetPrintStatus(kTRUE);

  fithandler->SetParamStartValues(paramInputValues);
  fithandler->SetResolutionConstants(); 
  fithandler->SetCrystalBallFunction(kTRUE);
  fithandler->SetMassWndLow(0.6);
  fithandler->SetMassWndHigh(0.4);
  //fithandler->FixParam(5,kTRUE); // Fix some parameters to their input values
  //fithandler->FixParam(14,kTRUE);
  //fithandler->FixParam(15,kTRUE);
  fithandler->SetMassWndLow(0.341696); // ----> M_low = 2.75 GeV
  fithandler->SetMassWndHigh(0.308304); // ----> M_high = 3.4 GeV

  aBtoJPSItoEle->DoMinimization();

  f->Close();
}
