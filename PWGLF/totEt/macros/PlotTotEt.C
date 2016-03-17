//Awful naming but this is also for plotting data
void PrintCent(int i);
float etaacc = 1.0;
//void PlotSimPbPb(char *filename = "rootFiles/LHC11a4_bis/Et.ESD.new.sim.LHC11a4_bis.root", Bool_t sim = kFALSE, Bool_t sysErr = kTRUE, char *corrfilename = "rootFiles/corrections/corrections.LHC11a4_bis.PbPb.ForData.root", Bool_t isEMCal = kTRUE, Bool_t isTPC = kTRUE){
void PlotTotEt(char *filename = "rootFiles/LHC15oPass1/Et.ESD.sim.LHC15oPass1.Run244918.root", Bool_t sim = kFALSE, Bool_t sysErr = kTRUE, char *corrfilename = "rootFiles/corrections/corrections.LHC15k1.PbPb.244918.ForData.2015.root", Bool_t isEMCal = kTRUE, Bool_t isTPC = kTRUE, Int_t nCB = 20){
// void PlotTotEt(char *filename = "rootFiles/LHC10h/Et.ESD.sim.LHC10h.Run139465.root", Bool_t sim = kFALSE, Bool_t sysErr = kTRUE, char *corrfilename = "rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.ForData.root", Bool_t isEMCal = kTRUE, Bool_t isTPC = kFALSE, Int_t nCB = 20){
  const Int_t nCentBins = nCB;
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
//   gROOT->LoadMacro("macros/loadCompiledLibraries.C");
//   loadCompiledLibraries();
//   TFile *corrfile = new TFile(corrfilename);
//   AliAnalysisHadEtCorrections *corrections = NULL;
//   if(isEMCal) corrections= (AliAnalysisHadEtCorrections *) corrfile->Get("hadCorrectionEMCAL");
//   else{ corrections= (AliAnalysisHadEtCorrections *) corrfile->Get("hadCorrectionPHOS");}
//   if(!corrections){cerr<<"Error could not find corrections in "<<corrfilename<<endl; corrfile->ls(); return;}

  gROOT->LoadMacro("macros/loadLocalHadLibraries.C");
  loadLocalHadLibraries();
  TFile *corrfile = new TFile(corrfilename);
  AliAnalysisHadEtCorrectionsLocal *corrections = NULL;
  if(isEMCal) corrections= (AliAnalysisHadEtCorrectionsLocal *) corrfile->Get("hadCorrectionEMCAL");
  else{ corrections= (AliAnalysisHadEtCorrectionsLocal *) corrfile->Get("hadCorrectionPHOS");}
  if(!corrections){cerr<<"Error could not find corrections in "<<corrfilename<<endl; corrfile->ls(); return;}


  etaacc = 2*corrections->GetEtaCut();
  cout<<"etaacc corr "<<etaacc<<endl;

  TFile *file = new TFile(filename);
  int nbins = nCentBins-2;
  Float_t npart[nCentBins] = {382.8,329.7,281.1,238.6,201.8,
		       169.5,141.6,116.2,94.09,75.45,
		       59.40,45.75,34.48,25.07,17.76,
		       12.49,8.810,6.150,4.351,3.060};
  Float_t nparterr[nCentBins]={3.1,4.6,4.8,4.2,3.8,
			4.2,3.5,3.1,3.6,2.6,
			2.0,1.5,1.6,1.1,1.4,
			1.0,0.37,0.21,0.18,0.087};
  Float_t npartnch[] = {382.8,329.7,260.5,186.4,128.9,85.0,52.8,30.0,15.8};
  Float_t nch[] = {1601,1294,966,649,426,261,149,76,35};
  Float_t nchErr[] = {60,49,37,23,15,9,6,4,2};
  Float_t etHadSim[nCentBins] = new Float_t[nCentBins];
  Float_t etTotSim[nCentBins] = new Float_t[nCentBins];
  Float_t etHadRec[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRec[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecCombinedEmcal[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecCombinedPhos[nCentBins] = new Float_t[nCentBins];
  Float_t etHadRecITS[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecITS[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecCombinedEmcalITS[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecCombinedPhosITS[nCentBins] = new Float_t[nCentBins];
  Float_t etNchRec[9] = new Float_t[9];
  Float_t etHadSimErr[nCentBins] = new Float_t[nCentBins];
  Float_t etTotSimErr[nCentBins] = new Float_t[nCentBins];
  Float_t etHadRecErr[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecErr[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecCombinedPhosErr[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecCombinedEmcalErr[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecCombinedPhosStatErr[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecCombinedEmcalStatErr[nCentBins] = new Float_t[nCentBins];
  Float_t etHadRecITSErr[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecITSErr[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecCombinedEmcalITSErr[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecCombinedPhosITSErr[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecCombinedEmcalITSStatErr[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecCombinedPhosITSStatErr[nCentBins] = new Float_t[nCentBins];
  Float_t etNchRecErr[9] = new Float_t[9];
  Float_t hadsyserr[nCentBins] = new Float_t[nCentBins];
  Float_t totsyserr[nCentBins] = new Float_t[nCentBins];
  Float_t totsyserrNoPID[nCentBins] = new Float_t[nCentBins];
  Float_t totsyserrNoPIDNegEta[nCentBins] = new Float_t[nCentBins];
  Float_t totsyserrNoPIDPosEta[nCentBins] = new Float_t[nCentBins];
  Float_t totsyserrNoPIDLimitedPhi[nCentBins] = new Float_t[nCentBins];

  Float_t etTotRecNoPID[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecITSNoPID[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecErrNoPID[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecITSErrNoPID[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecNoPIDNegEta[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecITSNoPIDNegEta[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecErrNoPIDNegEta[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecITSErrNoPIDNegEta[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecNoPIDPosEta[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecITSNoPIDPosEta[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecErrNoPIDPosEta[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecITSErrNoPIDPosEta[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecNoPIDLimitedPhi[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecITSNoPIDLimitedPhi[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecErrNoPIDLimitedPhi[nCentBins] = new Float_t[nCentBins];
  Float_t etTotRecITSErrNoPIDLimitedPhi[nCentBins] = new Float_t[nCentBins];

  TObjArray histoHadReco(nCentBins);
  TObjArray histoHadRecoITS(nCentBins);
  TObjArray histoTotReco(nCentBins);
  TObjArray histoTotRecoITS(nCentBins);
  TObjArray histoTotRecoNoPID(nCentBins);
  TObjArray histoTotRecoNoPIDITS(nCentBins);
  TObjArray histoTotRecoNoPIDNegEta(nCentBins);
  TObjArray histoTotRecoNoPIDNegEtaITS(nCentBins);
  TObjArray histoTotRecoNoPIDPosEta(nCentBins);
  TObjArray histoTotRecoNoPIDPosEtaITS(nCentBins);
  TObjArray histoTotRecoNoPIDLimitedPhi(nCentBins);
  TObjArray histoTotRecoNoPIDLimitedPhiITS(nCentBins);
  Int_t colors[] = {TColor::kRed,TColor::kOrange+7,TColor::kOrange+1,TColor::kOrange-2,TColor::kYellow+1,TColor::kSpring+9,TColor::kSpring-5,TColor::kGreen+1,TColor::kGreen-2,TColor::kTeal+5,TColor::kTeal-5,TColor::kCyan-7,TColor::kCyan-2,TColor::kAzure+2,TColor::kAzure+7,TColor::kAzure-1,TColor::kAzure-6,TColor::kBlue+3,1,1,1};

  for(int i=0;i<=nbins;i++){
    etHadSim[i]=0;
    etTotSim[i]=0;
    etHadRec[i]=0;
    etTotRec[i]=0;
    etHadRecITS[i]=0;
    etTotRecITS[i]=0;
    etHadSimErr[i]=0;
    etTotSimErr[i]=0;
    etHadRecErr[i]=0;
    etTotRecErr[i]=0;
    etHadRecITSErr[i]=0;
    etTotRecITSErr[i]=0;
    histoHadReco[i] =  out2->FindObject(Form("RecoHadEtFullAcceptanceTPCCB%i",i));
    ((TH1D*)histoHadReco[i])->GetXaxis()->SetRange(2, ((TH1D*)histoHadReco[i])->GetNbinsX() );
    etHadRec[i]=((TH1D*)histoHadReco[i])->GetMean();
    etHadRecErr[i]=((TH1D*)histoHadReco[i])->GetMeanError();
    histoHadRecoITS[i] =  out2->FindObject(Form("RecoHadEtFullAcceptanceITSCB%i",i));
    ((TH1D*)histoHadRecoITS[i])->GetXaxis()->SetRange(2, ((TH1D*)histoHadRecoITS[i])->GetNbinsX() );
    etHadRecITS[i]=((TH1D*)histoHadRecoITS[i])->GetMean();
    etHadRecITSErr[i]=((TH1D*)histoHadRecoITS[i])->GetMeanError();

    histoTotReco[i] =  out2->FindObject(Form("RecoTotEtFullAcceptanceTPCCB%i",i));
    ((TH1D*)histoTotReco[i])->GetXaxis()->SetRange(2, ((TH1D*)histoTotReco[i])->GetNbinsX() );
    etTotRec[i]=((TH1D*)histoTotReco[i])->GetMean();
    etTotRecErr[i]=((TH1D*)histoTotReco[i])->GetMeanError();
    histoTotRecoITS[i] =  out2->FindObject(Form("RecoTotEtFullAcceptanceITSCB%i",i));
    ((TH1D*)histoTotRecoITS[i])->GetXaxis()->SetRange(2, ((TH1D*)histoTotRecoITS[i])->GetNbinsX() );
    etTotRecITS[i]=((TH1D*)histoTotRecoITS[i])->GetMean();
    etTotRecITSErr[i]=((TH1D*)histoTotRecoITS[i])->GetMeanError();//had

    histoTotRecoNoPID[i] =  out2->FindObject(Form("RecoTotEtFullAcceptanceTPCNoPIDCB%i",i));
    ((TH1D*)histoTotRecoNoPID[i])->GetXaxis()->SetRange(2, ((TH1D*)histoTotRecoNoPID[i])->GetNbinsX() );
    etTotRecNoPID[i]=((TH1D*)histoTotRecoNoPID[i])->GetMean();
    etTotRecErrNoPID[i]=((TH1D*)histoTotRecoNoPID[i])->GetMeanError();
    cout<<"i "<<i<< " ET "<<etTotRecNoPID[i]<<" +/- "<<etTotRecErrNoPID[i]<<" entries "<<((TH1D*)histoTotRecoNoPID[i])->GetEntries()<<" name "<<Form("RecoTotEtFullAcceptanceTPCNoPIDCB%i",i)<<endl;
    histoTotRecoNoPIDITS[i] =  out2->FindObject(Form("RecoTotEtFullAcceptanceITSNoPIDCB%i",i));
    ((TH1D*)histoTotRecoNoPIDITS[i])->GetXaxis()->SetRange(2, ((TH1D*)histoTotRecoNoPIDITS[i])->GetNbinsX() );
    etTotRecITSNoPID[i]=((TH1D*)histoTotRecoNoPIDITS[i])->GetMean();
    etTotRecITSErrNoPID[i]=((TH1D*)histoTotRecoNoPIDITS[i])->GetMeanError();//Reco

    histoTotRecoNoPIDNegEta[i] =  out2->FindObject(Form("RecoTotEtFullAcceptanceTPCNoPIDNegEtaCB%i",i));
    ((TH1D*)histoTotRecoNoPIDNegEta[i])->GetXaxis()->SetRange(2, ((TH1D*)histoTotRecoNoPIDNegEta[i])->GetNbinsX() );
    etTotRecNoPIDNegEta[i]=((TH1D*)histoTotRecoNoPIDNegEta[i])->GetMean();
    etTotRecErrNoPIDNegEta[i]=((TH1D*)histoTotRecoNoPIDNegEta[i])->GetMeanError();

    histoTotRecoNoPIDPosEta[i] =  out2->FindObject(Form("RecoTotEtFullAcceptanceTPCNoPIDPosEtaCB%i",i));
    ((TH1D*)histoTotRecoNoPIDPosEta[i])->GetXaxis()->SetRange(2, ((TH1D*)histoTotRecoNoPIDPosEta[i])->GetNbinsX() );
    etTotRecNoPIDPosEta[i]=((TH1D*)histoTotRecoNoPIDPosEta[i])->GetMean();
    etTotRecErrNoPIDPosEta[i]=((TH1D*)histoTotRecoNoPIDPosEta[i])->GetMeanError();


    histoTotRecoNoPIDLimitedPhi[i] =  out2->FindObject(Form("RecoTotEtFullAcceptanceTPCNoPIDLimitedPhiCB%i",i));
    ((TH1D*)histoTotRecoNoPIDLimitedPhi[i])->GetXaxis()->SetRange(2, ((TH1D*)histoTotRecoNoPIDLimitedPhi[i])->GetNbinsX() );
    etTotRecNoPIDLimitedPhi[i]=((TH1D*)histoTotRecoNoPIDLimitedPhi[i])->GetMean();
    etTotRecErrNoPIDLimitedPhi[i]=((TH1D*)histoTotRecoNoPIDLimitedPhi[i])->GetMeanError();


  }
  TCanvas *c1 = new TCanvas("c1","Reconstructed Et",1200,1000);
  c1->SetTopMargin(0.04);
  c1->SetRightMargin(0.04);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.159091);
  c1->SetLogy();
  c1->Divide(2,2);
  c1->cd(1)->SetLogy();
  TH1D *histoHadRecoAll = out2->FindObject("RecoHadEtFullAcceptanceTPC");
  histoHadRecoAll->Rebin(2);
  histoHadRecoAll->Draw();
  for(int i=0;i<=nbins;i++){
    ((TH1D*)histoHadReco[i])->Rebin(2);
    ((TH1D*)histoHadReco[i])->SetLineColor(colors[i]);
    ((TH1D*)histoHadReco[i])->Draw("same");
  }
  c1->cd(2)->SetLogy();
  TH1D *histoHadRecoITSAll = out2->FindObject("RecoHadEtFullAcceptanceITS");
  histoHadRecoITSAll->Rebin(2);
  histoHadRecoITSAll->Draw();
  for(int i=0;i<=nbins;i++){
    ((TH1D*)histoHadRecoITS[i])->Rebin(2);
    ((TH1D*)histoHadRecoITS[i])->SetLineColor(colors[i]);
    ((TH1D*)histoHadRecoITS[i])->Draw("same");
  }
  c1->cd(3)->SetLogy();
  TH1D *histoTotRecoAll = out2->FindObject("RecoTotEtFullAcceptanceTPC");
  histoTotRecoAll->Rebin(2);
  histoTotRecoAll->Draw();
  for(int i=0;i<=nbins;i++){
    ((TH1D*)histoTotReco[i])->Rebin(2);
    ((TH1D*)histoTotReco[i])->SetLineColor(colors[i]);
    ((TH1D*)histoTotReco[i])->Draw("same");
  }
  c1->cd(4)->SetLogy();
  TH1D *histoTotRecoITSAll = out2->FindObject("RecoTotEtFullAcceptanceITS");
  histoTotRecoITSAll->Rebin(2);
  histoTotRecoITSAll->Draw();
  for(int i=0;i<=nbins;i++){
    ((TH1D*)histoTotRecoITS[i])->Rebin(2);
    ((TH1D*)histoTotRecoITS[i])->SetLineColor(colors[i]);
    ((TH1D*)histoTotRecoITS[i])->Draw("same");
  }


  for(int i=0;i<=nbins;i++){
    cout<<"npart "<<npart[i]<<" reco "<<etTotRec[i]/etaacc<<" +/- "<<etTotRecErr[i]/etaacc<<" reco ITS "<<etTotRecITS[i]/etaacc<<" +/- "<<etTotRecITSErr[i]/etaacc;
    cout<<" reco "<<etHadRec[i]/etaacc<<" +/- "<<etHadRecErr[i]/etaacc<<" reco ITS "<<etHadRecITS[i]/etaacc<<" +/- "<<etHadRecITSErr[i]/etaacc;
    cout<<endl;
  }

  Float_t scale = 2.0/etaacc;//factor of 1.4 gives dN/deta
  TGraphErrors *graphHadEtReco = new TGraphErrors(18);
  TGraphAsymmErrors *graphHadEtRecoErrors = new TGraphAsymmErrors(18);

  ofstream myfileHad;
  ofstream myfileTot;
  ofstream myfileTotNoPID;
  ofstream myfileTotNoPIDNegEta;
  ofstream myfileTotNoPIDPosEta;
  ofstream myfileTotNoPIDLimitedPhi;
  TString texfilenameHad = "HadEtFromTracks.dat";
  TString texfilenameTot = "TotEtFromTracks.dat";
  TString texfilenameTotNoPID = "TotEtFromTracksNoPID.dat";
  TString texfilenameTotNoPIDNegEta = "TotEtFromTracksNoPIDNegEta.dat";
  TString texfilenameTotNoPIDPosEta = "TotEtFromTracksNoPIDPosEta.dat";
  TString texfilenameTotNoPIDLimitedPhi = "TotEtFromTracksNoPIDLimitedPhi.dat";
  myfileHad.open (texfilenameHad.Data());
  myfileTot.open (texfilenameTot.Data());
  myfileTotNoPID.open (texfilenameTotNoPID.Data());
  myfileTotNoPIDNegEta.open (texfilenameTotNoPIDNegEta.Data());
  myfileTotNoPIDPosEta.open (texfilenameTotNoPIDPosEta.Data());
  myfileTotNoPIDLimitedPhi.open (texfilenameTotNoPIDLimitedPhi.Data());


  cout<<"Hadronic and total et errors"<<endl;
  for(int i=0;i<nbins;i++){
    float valueHadronic = etHadRec[i]/npart[i]*scale;
    float value = valueHadronic;
    graphHadEtReco->SetPoint(i,npart[i],valueHadronic);
    graphHadEtReco->SetPointError(i,0.0,scale*etHadRecErr[i]/npart[i]);
    graphHadEtRecoErrors->SetPoint(i,npart[i],valueHadronic);
    graphHadEtRecoErrors->SetPointEXlow(i,0.0);
    graphHadEtRecoErrors->SetPointEXhigh(i,0.0);
    float low = corrections->GetSystematicErrorBound(etHadRec[i],kTRUE,kTRUE, isTPC);
    float high = corrections->GetSystematicErrorBound(etHadRec[i],kFALSE,kTRUE, isTPC);
    float lowtot = corrections->GetSystematicErrorBound(etTotRec[i],kTRUE,kFALSE, isTPC);
    float hightot = corrections->GetSystematicErrorBound(etTotRec[i],kFALSE,kFALSE, isTPC);
    float lowtotnopid = corrections->GetSystematicErrorBound(etTotRecNoPID[i],kTRUE,kFALSE, isTPC,kFALSE);
    float hightotnopid = corrections->GetSystematicErrorBound(etTotRecNoPID[i],kFALSE,kFALSE, isTPC,kFALSE);
    float lowtotnopidnegeta = corrections->GetSystematicErrorBound(etTotRecNoPIDNegEta[i],kTRUE,kFALSE, isTPC,kFALSE);
    float hightotnopidnegeta = corrections->GetSystematicErrorBound(etTotRecNoPIDNegEta[i],kFALSE,kFALSE, isTPC,kFALSE);
    float lowtotnopidposeta = corrections->GetSystematicErrorBound(etTotRecNoPIDPosEta[i],kTRUE,kFALSE, isTPC,kFALSE);
    float hightotnopidposeta = corrections->GetSystematicErrorBound(etTotRecNoPIDPosEta[i],kFALSE,kFALSE, isTPC,kFALSE);
    float lowtotnopidlimitedphi = corrections->GetSystematicErrorBound(etTotRecNoPIDLimitedPhi[i],kTRUE,kFALSE, isTPC,kFALSE);
    float hightotnopidlimitedphi = corrections->GetSystematicErrorBound(etTotRecNoPIDLimitedPhi[i],kFALSE,kFALSE, isTPC,kFALSE);
    float errFracNpart = nparterr[i]/npart[i];
    float syserrorLow = scale/npart[i] * low; 
    float syserrorHigh = scale/npart[i] * high;
    float totalErrLow = 0;
    if(value>0) totalErrLow = value*TMath::Sqrt(TMath::Power((value - syserrorLow)/value,2) + TMath::Power(errFracNpart,2));
    float totalErrHigh = 0;
    if(value>0) totalErrHigh = value*TMath::Sqrt(TMath::Power((value - syserrorHigh)/value,2) + TMath::Power(errFracNpart,2));
    graphHadEtRecoErrors->SetPointEYlow(i,totalErrLow);
    graphHadEtRecoErrors->SetPointEYhigh(i,totalErrHigh);

    hadsyserr[i] = TMath::Abs(high-etHadRec[i]);
    if(TMath::Abs(low-etHadRec[i])>hadsyserr[i]) hadsyserr[i] = TMath::Abs(low-etHadRec[i]);
    totsyserr[i] = TMath::Abs(hightot-etTotRec[i]);
    if(TMath::Abs(lowtot-etTotRec[i])>totsyserr[i]) totsyserr[i] = TMath::Abs(lowtot-etTotRec[i]);
    totsyserrNoPID[i] = TMath::Abs(hightotnopid-etTotRecNoPID[i]);
    if(TMath::Abs(lowtotnopid-etTotRecNoPID[i])>totsyserrNoPID[i]) totsyserrNoPID[i] = TMath::Abs(lowtotnopid-etTotRecNoPID[i]);
    totsyserrNoPIDNegEta[i] = TMath::Abs(hightotnopidnegeta-etTotRecNoPIDNegEta[i]);
    if(TMath::Abs(lowtotnopidnegeta-etTotRecNoPIDNegEta[i])>totsyserrNoPIDNegEta[i]) totsyserrNoPIDNegEta[i] = TMath::Abs(lowtotnopidnegeta-etTotRecNoPIDNegEta[i]);
    totsyserrNoPIDPosEta[i] = TMath::Abs(hightotnopidposeta-etTotRecNoPIDPosEta[i]);
    if(TMath::Abs(lowtotnopidposeta-etTotRecNoPIDPosEta[i])>totsyserrNoPIDPosEta[i]) totsyserrNoPIDPosEta[i] = TMath::Abs(lowtotnopidposeta-etTotRecNoPIDPosEta[i]);
    totsyserrNoPIDLimitedPhi[i] = TMath::Abs(hightotnopidlimitedphi-etTotRecNoPIDLimitedPhi[i]);
    if(TMath::Abs(lowtotnopidlimitedphi-etTotRecNoPIDLimitedPhi[i])>totsyserrNoPIDLimitedPhi[i]) totsyserrNoPIDLimitedPhi[i] = TMath::Abs(lowtotnopidlimitedphi-etTotRecNoPID[i]);
    
    float myscale = 1/etaacc;
    if(i<nbins){
      PrintCent(i);
      myfileHad<<Form("%i\t%i\t%2.1f\t%2.1f\t%2.1f \n",i*5,i*5+5,myscale*etHadRec[i],myscale*etHadRecErr[i],myscale*hadsyserr[i]);
      myfileTot<<Form("%i\t%i\t%2.1f\t%2.1f\t%2.1f \n",i*5,i*5+5,myscale*etTotRec[i],myscale*etTotRecErr[i],myscale*totsyserr[i]);
      myfileTotNoPID<<Form("%i\t%i\t%2.1f\t%2.1f\t%2.1f \n",i*5,i*5+5,myscale*etTotRecNoPID[i],myscale*etTotRecErrNoPID[i],myscale*totsyserrNoPID[i]);
      myfileTotNoPIDNegEta<<Form("%i\t%i\t%2.1f\t%2.1f\t%2.1f \n",i*5,i*5+5,myscale*etTotRecNoPIDNegEta[i],myscale*etTotRecErrNoPIDNegEta[i],myscale*totsyserrNoPIDNegEta[i]);
      myfileTotNoPIDPosEta<<Form("%i\t%i\t%2.1f\t%2.1f\t%2.1f \n",i*5,i*5+5,myscale*etTotRecNoPIDPosEta[i],myscale*etTotRecErrNoPIDPosEta[i],myscale*totsyserrNoPIDPosEta[i]);
      myfileTotNoPIDLimitedPhi<<Form("%i\t%i\t%2.1f\t%2.1f\t%2.1f \n",i*5,i*5+5,myscale*etTotRecNoPIDLimitedPhi[i],myscale*etTotRecErrNoPIDLimitedPhi[i],myscale*totsyserrNoPIDLimitedPhi[i]);
      cout<<Form("%2.1f & ",npart[i]);
      cout<<Form("%2.1f $\\pm$ %2.1f $\\pm$ %2.1f & ",myscale*etHadRec[i],myscale*etHadRecErr[i],myscale*hadsyserr[i]);
      cout<<Form("%2.1f $\\pm$ %2.1f $\\pm$ %2.1f",myscale*etTotRec[i],myscale*etTotRecErr[i],myscale*totsyserr[i]);
      cout<<" \\\\ "<<endl;
    }
  }
  myfileHad.close();
  myfileTot.close();
  myfileTotNoPID.close();
  myfileTotNoPIDNegEta.close();
  myfileTotNoPIDPosEta.close();
  myfileTotNoPIDLimitedPhi.close();

  graphHadEtReco->GetYaxis()->SetTitle("dE_{T}^{had}/d#eta/0.5 N_{part}");
  graphHadEtReco->GetXaxis()->SetTitle("N_{part}");
  graphHadEtReco->SetMarkerStyle(20);
  TGraphErrors *graphTotEtReco = new TGraphErrors(18);
  TGraphAsymmErrors *graphTotEtRecoErrors = new TGraphAsymmErrors(18);
  for(int i=0;i<nbins;i++){
    float value = scale*etTotRec[i]/npart[i];
    graphTotEtReco->SetPoint(i,npart[i],value);
    graphTotEtReco->SetPointError(i,0.0,scale*etTotRecErr[i]/npart[i]);
    graphTotEtRecoErrors->SetPoint(i,npart[i],scale*etTotRec[i]/npart[i]);
    graphTotEtRecoErrors->SetPointEXlow(i,0.0);
    graphTotEtRecoErrors->SetPointEXhigh(i,0.0);
    //et, low bound, is hadronic, is TPC
    float syserrorLow = scale/npart[i] * corrections->GetSystematicErrorBound(etTotRec[i],kTRUE,kFALSE, isTPC);
    float syserrorHigh = scale/npart[i] * corrections->GetSystematicErrorBound(etTotRec[i],kFALSE,kFALSE, isTPC);
    float totalErrLow = 0; 
    if(value>0) totalErrLow = value*TMath::Sqrt(TMath::Power((value - syserrorLow)/value,2) + TMath::Power(errFracNpart,2));
    float totalErrHigh = 0; 
    if(value>0) totalErrHigh = totalErrHigh = value*TMath::Sqrt(TMath::Power((value - syserrorHigh)/value,2) + TMath::Power(errFracNpart,2));
    //cout<<"sys errors "<<i<<" "<<value<<" high "<<syserrorHigh<<" low "<<syserrorLow<<endl;
    graphTotEtRecoErrors->SetPointEYlow(i,totalErrLow);
    graphTotEtRecoErrors->SetPointEYhigh(i,totalErrHigh);
  }
  graphTotEtReco->GetYaxis()->SetTitle("dE_{T}/d#eta/0.5 N_{part}");
  graphTotEtReco->GetXaxis()->SetTitle("N_{part}");
  graphTotEtReco->SetMarkerStyle(20);



  TCanvas *c5 = new TCanvas("c5","Reconstructed Ethad vs npart",700,500);
  c5->SetTopMargin(0.04);
  c5->SetRightMargin(0.04);
  c5->SetBorderSize(0);
  c5->SetFillColor(0);
  c5->SetFillColor(0);
  c5->SetBorderMode(0);
  c5->SetFrameFillColor(0);
  c5->SetFrameBorderMode(0);
  c5->SetLeftMargin(0.159091);
  //c5->SetLogy();
  graphHadEtReco->Draw("AP");
  graphHadEtRecoErrors->Draw("[]");

  TCanvas *c6 = new TCanvas("c6","Reconstructed Et vs npart",700,500);
  c6->SetTopMargin(0.04);
  c6->SetRightMargin(0.04);
  c6->SetBorderSize(0);
  c6->SetFillColor(0);
  c6->SetFillColor(0);
  c6->SetBorderMode(0);
  c6->SetFrameFillColor(0);
  c6->SetFrameBorderMode(0);
  c6->SetLeftMargin(0.159091);
  //c6->SetLogy();
  graphTotEtReco->Draw("AP");
  graphTotEtRecoErrors->Draw("[]");


  //create latex summary tables
  ofstream myfileHadEtSummaryTable;
  TString texfilenameHadEtSummaryTable = "HadEtSummaryTable.tex";
  myfileHadEtSummaryTable.open (texfilenameHadEtSummaryTable.Data());
  Float_t fnotid = corrections->GetNotIDConstCorrectionTPC();
  Float_t fnotidlow = corrections->GetNotIDConstCorrectionTPCLowBound();
  myfileHadEtSummaryTable<<"\\fptcut & "<<Form("%2.3f \$\\pm\$ %2.3f & %2.1f \\\% \\\\",1.0/corrections->GetpTCutCorrectionTPC(),1.0/corrections->GetpTCutCorrectionTPCLowBound()-1.0/corrections->GetpTCutCorrectionTPC(),(corrections->GetpTCutCorrectionTPC()-corrections->GetpTCutCorrectionTPCLowBound())/corrections->GetpTCutCorrectionTPC()*100.0)<<endl;
  myfileHadEtSummaryTable<<"\\fneutral & "<<Form("%2.3f \$\\pm\$ %2.3f & %2.1f \\\% \\\\",1.0/corrections->GetNeutralCorrection(),1.0/corrections->GetNeutralCorrectionLowBound()-1.0/corrections->GetNeutralCorrection(),(corrections->GetNeutralCorrection()-corrections->GetNeutralCorrectionLowBound())/corrections->GetNeutralCorrection()*100.0)<<endl;
  myfileHadEtSummaryTable<<"\\ftotal & "<<Form("%2.3f \$\\pm\$ %2.3f & %2.1f \\\% \\\\",1.0/corrections->GetNotHadronicCorrection(),1.0/corrections->GetNotHadronicCorrectionLowBound()-1.0/corrections->GetNotHadronicCorrection(),(corrections->GetNotHadronicCorrection()-corrections->GetNotHadronicCorrectionLowBound())/corrections->GetNotHadronicCorrection()*100.0)<<endl;
  myfileHadEtSummaryTable<<"\\fnotid & "<<Form("%2.3f \$\\pm\$ %2.3f & %2.1f \\\% \\\\",1.0/fnotid, 1.0/fnotidlow - 1.0/fnotid,(fnotid - fnotidlow)/fnotid*100.0)<<endl;
  myfileHadEtSummaryTable<<"\\fbkgd & 1.8\\% & 0.8\\% \\\\"<<endl;//hard coded because it's not easy to put in
  myfileHadEtSummaryTable<<"\\eff & 40\\% & 3\\% \\\\"<<endl;//hard coded because it's not easy to put in
  myfileHadEtSummaryTable.close();

  ofstream myfile4;
  myfile4.open ("TotEtFromTracks.tex");
  ofstream myfile5;
  myfile5.open ("HadEtFromTracks.tex");
  ofstream myfile6;
  myfile6.open ("TotEtFromTracksOnly.tex");
  for(int i=0;i<20;i++){
    if(etHadRec[i]>1e-3){
      myfile5 <<Form("%i-%i\\% & %2.1f\t%2.1f \$\\pm\$ %2.1f \\\\\n",i*5,i*5+5,myscale*etHadRec[i],myscale*etHadRecErr[i],myscale*hadsyserr[i]);
      myfile6 <<Form("%i-%i\\%  & %2.1f \$\\pm\$ %2.1f \$\\pm\$ %2.1f \\\\ \n",i*5,i*5+5,myscale*etTotRec[i],myscale*etTotRecErr[i],myscale*totsyserr[i]);
      myfile4 <<Form("%i-%i\\%  & %2.1f \$\\pm\$ %2.1f \$\\pm\$ %2.1f &",i*5,i*5+5,myscale*etTotRec[i],myscale*etTotRecErr[i],myscale*totsyserr[i]);
      if(etTotRecCombinedEmcal[i]>1e-3 && i<10){
	myfile4 <<Form("%2.1f \$\\pm\$ %2.1f \$\\pm\$ %2.1f &",etTotRecCombinedEmcal[i],etTotRecCombinedEmcalStatErr[i],etTotRecCombinedEmcalErr[i]);
      }
      else{myfile4<<" -- & ";}
      if(etTotRecCombinedPhos[i]>1e-3 && i<6){
	myfile4 <<Form("%2.1f \$\\pm\$ %2.1f \$\\pm\$ %2.1f \\\\ \n",etTotRecCombinedPhos[i],etTotRecCombinedPhosStatErr[i],etTotRecCombinedPhosErr[i]);
      }
      else{myfile4<<" -- \\\\ \n ";}
    }
  }

  myfile6.close();
  myfile5.close();
  myfile4.close();

}




void PrintCent(int i){
  switch(i){
  case 0:
    cout<<"0-5\\% &";
    break;
  case 1:
    cout<<"5-10\\% &";
    break;
  case 2:
    cout<<"10-15\\% &";
    break;
  case 3:
    cout<<"15-20\\% &";
    break;
  case 4:
    cout<<"20-25\\% &";
    break;
  case 5:
    cout<<"25-30\\% &";
    break;
  case 6:
    cout<<"30-35\\% &";
    break;
  case 7:
    cout<<"35-40\\% &";
    break;
  case 8:
    cout<<"40-45\\% &";
    break;
  case 9:
    cout<<"45-50\\% &";
    break;
  case 10:
    cout<<"50-55\\% &";
    break;
  case 11:
    cout<<"55-60\\% &";
    break;
  case 12:
    cout<<"60-65\\% &";
    break;
  case 13:
    cout<<"65-70\\% &";
    break;
  case 14:
    cout<<"70-75\\% &";
    break;
  case 15:
    cout<<"75-80\\% &";
    break;
  case 16:
    cout<<"80-85\\% &";
    break;
  case 17:
    cout<<"85-90\\% &";
    break;
  case 18:
    cout<<"90-95\\% &";
    break;
  case 19:
    cout<<"95-100\\% &";
    break;
  }
}
