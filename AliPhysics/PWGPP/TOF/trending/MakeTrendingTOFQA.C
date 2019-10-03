/*
  fbellini@cern.ch - last update on 21/02/2014
  Macro to run the TOF QA trending by accessing the std QA output, 
  to be mainly used with the automatic scripts to fill the QA repository.
  Launch with 
  aliroot -l -b -q "MakeTrendingTOFQA.C(\"${fullpath}/QAresults.root\", ${run}, ...) 
  The macro produces a file containing the tree of trending variables and the main plots.
  A feature that displays the plots in canvases must be enable when needed.
*/

Int_t MakeTrendingTOFQA(TString qafilename,       //full path of the QA output; set IsOnGrid to prepend "alien://"
			Int_t runNumber,          // run number
			Bool_t isMC=kFALSE,       //MC flag, to disable meaningless checks
			Bool_t canvasE = kFALSE,  //enable display plots on canvas and save png
			Bool_t IsOnGrid = kFALSE, //set to kTRUE to access files on the grid
			TString ocdbStorage = "raw://") //set the default ocdb storage
{   
  // macro to generate tree with TOF QA trending variables
  // access qa PWGPP output files  
  if (!qafilename) {
    Printf("Error - Invalid input file");
    return 1;
  }

  /*set graphic style*/
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTitleBorderSize(0)  ;
  gStyle->SetTitleFont(42);
  gStyle->SetTextFont(42);
  gStyle->SetStatColor(kWhite); 
  gStyle->SetStatBorderSize(1);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(10);

  char defaultQAoutput[30]="QAresults.root";
  char * treePostFileName="trending.root";
  
  if (IsOnGrid) TGrid::Connect("alien://");
  TFile * fin = TFile::Open(qafilename,"r");
  if (!fin) {
    Printf("ERROR: QA output not found. Exiting...\n");
    return -1;
  } else {
    Printf("INFO: QA output file %s open. \n",fin->GetName());
  }
  
  //access histograms lists
  char tofQAdirName[15]="TOF_Performance";
  char genListName[15]="cGeneralTOFqa";
  char t0ListName[15]="cTimeZeroTOFqa";
  char pidListName[15]="cPIDTOFqa"; 
  char posListName[15]="cPositiveTOFqa";
  char negListName[15]="cNegativeTOFqa";
  
  TDirectoryFile * tofQAdir=(TDirectoryFile*)fin->Get(tofQAdirName);
  if (!tofQAdir) {
    Printf("ERROR: TOF QA directory not present in input file.\n");
    return -1;
  }
  TList * generalList=(TList*)tofQAdir->Get(genListName);
  TList  *timeZeroList=(TList*)tofQAdir->Get(t0ListName);
  TList  *pidList=(TList*)tofQAdir->Get(pidListName);
  TList  *posList=(TList*)tofQAdir->Get(posListName);
  TList  *negList=(TList*)tofQAdir->Get(negListName);
  
  if (!generalList) Printf("WARNING: general QA histograms absent or not accessible\n");
  if (!timeZeroList) Printf("WARNING: timeZero QA histograms absent or not accessible\n");
  if (!pidList) Printf("WARNING: PID QA histograms absent or not accessible\n");
  if (!posList) Printf("WARNING: general QA histograms for positive tracks absent or not accessible\n");
  if (!negList) Printf("WARNING: general QA histograms for negative tracks absent or not accessible\n");
  
  if ( (!generalList) && (!timeZeroList) && (!pidList) ){
    printf("ERROR: no QA available \n");
    return -1;
  }
  
  Printf(":::: Getting post-analysis info for run %i",runNumber);
  TFile * trendFile = new TFile(treePostFileName,"recreate");

  Double_t avTime=-9999., peakTime=-9999., spreadTime=-9999., peakTimeErr=-9999., spreadTimeErr=-9999., negTimeRatio=-9999.,
    avRawTime=-9999., peakRawTime=-9999., spreadRawTime=-9999., peakRawTimeErr=-9999., spreadRawTimeErr=-9999., avTot=-9999., peakTot=-9999.,spreadTot=-9999.,  peakTotErr=-9999.,spreadTotErr=-9999.,
    orphansRatio=-9999., avL=-9999., negLratio=-9999.,
    effPt1=-9999., effPt2=-9999., matchEffLinFit1Gev=-9999.,matchEffLinFit1GevErr=-9999.;
  
  Double_t avPiDiffTime=-9999.,peakPiDiffTime=-9999., spreadPiDiffTime=-9999.,peakPiDiffTimeErr=-9999., spreadPiDiffTimeErr=-9999.;
  
  Double_t avT0A=-9999.,peakT0A=-9999., spreadT0A=-9999.,peakT0AErr=-9999., spreadT0AErr=-9999.;
  Double_t avT0C=-9999.,peakT0C=-9999., spreadT0C=-9999.,peakT0CErr=-9999., spreadT0CErr=-9999.;
  Double_t avT0AC=-9999.,peakT0AC=-9999., spreadT0AC=-9999.,peakT0ACErr=-9999., spreadT0ACErr=-9999.;
  Double_t avT0res=-9999.,peakT0res=-9999., spreadT0res=-9999.,peakT0resErr=-9999., spreadT0resErr=-9999.;
  Double_t avT0fillRes=-9999.;

  Float_t avMulti=0;
  Float_t fractionEventsWHits=-9999.;
  /*number of good (HW ok && efficient && !noisy) TOF channels from OCDB*/
  Double_t goodChannelRatio=0.0;

  TTree * ttree=new TTree("trending","tree of trending variables");
  ttree->Branch("run",&runNumber,"run/I");
  ttree->Branch("avMulti",&avMulti,"avMulti/F");
  ttree->Branch("fractionEventsWHits",&fractionEventsWHits,"fractionEventsWHits/F");
  ttree->Branch("goodChannelsRatio",&goodChannelRatio,"goodChannelRatio/D");
  ttree->Branch("avTime",&avTime,"avTime/D"); //mean time
  ttree->Branch("peakTime",&peakTime,"peakTime/D"); //main peak time after fit
  ttree->Branch("spreadTime",&spreadTime,"spreadTime/D"); //spread of main peak of time after fit
  ttree->Branch("peakTimeErr",&peakTimeErr,"peakTimeErr/D"); //main peak time after fit error
  ttree->Branch("spreadTimeErr",&spreadTimeErr,"spreadTimeErr/D"); //spread of main peak of time after fit error
  ttree->Branch("negTimeRatio",&negTimeRatio,"negTimeRatio/D"); //negative time ratio
  
  ttree->Branch("avRawTime",&avRawTime,"avRawTime/D"); //mean raw time
  ttree->Branch("peakRawTime",&peakRawTime,"peakRawTime/D"); //mean peak of RAW TIME after fit
  ttree->Branch("spreadRawTime",&spreadRawTime,"spreadRawTime/D"); //spread of main peak of raw time after fit
  ttree->Branch("peakRawTimeErr",&peakRawTimeErr,"peakRawTimeErr/D"); //main peak raw  time after fit error
  ttree->Branch("spreadRawTimeErr",&spreadRawTimeErr,"spreadRawTimeErr/D"); //spread of  raw main peak of time after fit error
  
  ttree->Branch("avTot",&avTot,"avTot/D"); //main peak tot
  ttree->Branch("peakTot",&peakTot,"peakTot/D"); // main peak of tot after fit
  ttree->Branch("spreadTot",&spreadTot,"spreadTot/D"); //spread of main peak of tot after fit
  ttree->Branch("peakTotErr",&peakTotErr,"peakTotErr/D"); // main peak of tot after fit
  ttree->Branch("spreadTotErr",&spreadTotErr,"spreadTotErr/D"); //spread of main peak of tot after fit
  
  ttree->Branch("orphansRatio",&orphansRatio,"orphansRatio/D"); //orphans ratio

  ttree->Branch("avL",&avL,"avL/D"); //mean track length
  ttree->Branch("negLratio",&negLratio,"negLratio/D");//ratio of tracks with track length <350 cm
  ttree->Branch("effPt1",&effPt1,"effPt1/D");//matching eff at 1 GeV/c
  ttree->Branch("effPt2",&effPt2,"effPt2/D"); //matching eff at 2 GeV/c
  ttree->Branch("matchEffLinFit1Gev",&matchEffLinFit1Gev,"matchEffLinFit1Gev/D");//matching eff fit param 
  ttree->Branch("matchEffLinFit1GevErr",&matchEffLinFit1GevErr,"matchEffLinFit1GevErr/D");////matching eff fit param error
  
  ttree->Branch("avPiDiffTime",&avPiDiffTime,"avPiDiffTime/D"); //mean t-texp
  ttree->Branch("peakPiDiffTime",&peakPiDiffTime,"peakPiDiffTime/D"); //main peak t-texp after fit
  ttree->Branch("spreadPiDiffTime",&spreadPiDiffTime,"spreadPiDiffTime/D"); //spread of main peak t-texp after fit
  ttree->Branch("peakPiDiffTimeErr",&peakPiDiffTimeErr,"peakPiDiffTimeErr/D"); //main peak t-texp after fit error
  ttree->Branch("spreadPiDiffTimeErr",&spreadPiDiffTimeErr,"spreadPiDiffTimeErr/D"); //spread of main peak of t-texp after fit error

  ttree->Branch("avT0A",&avT0A,"avT0A/D"); //main peak t0A
  ttree->Branch("peakT0A",&peakT0A,"peakT0A/D"); // main peak of t0A after fit
  ttree->Branch("spreadT0A",&spreadT0A,"spreadTot/D"); //spread of main peak of t0A after fit
  ttree->Branch("peakT0AErr",&peakT0AErr,"peakT0AErr/D"); // main peak of t0A after fit
  ttree->Branch("spreadT0AErr",&spreadT0AErr,"spreadT0AErr/D"); //spread of main peak of t0A after fit

  ttree->Branch("avT0C",&avT0C,"avT0C/D"); //main peak t0C
  ttree->Branch("peakT0C",&peakT0C,"peakT0C/D"); // main peak of t0C after fit
  ttree->Branch("spreadT0C",&spreadT0C,"spreadT0C/D"); //spread of main peak of t0C after fit
  ttree->Branch("peakT0CErr",&peakT0CErr,"peakT0CErr/D"); // main peak of t0C after fit
  ttree->Branch("spreadT0CErr",&spreadT0CErr,"spreadT0CErr/D"); //spread of main peak of t0C after fit
 
  ttree->Branch("avT0AC",&avT0AC,"avT0AC/D"); //main peak t0AC
  ttree->Branch("peakT0AC",&peakT0AC,"peakT0AC/D"); // main peak of t0AC after fit
  ttree->Branch("spreadT0AC",&spreadT0AC,"spreadT0AC/D"); //spread of main peak of t0AC after fit
  ttree->Branch("peakT0ACErr",&peakT0ACErr,"peakT0ACErr/D"); // main peak of t0AC after fit
  ttree->Branch("spreadT0ACErr",&spreadT0ACErr,"spreadT0ACErr/D"); //spread of main peak of t0AC after fit
 
  ttree->Branch("avT0res",&avT0res,"avT0res/D"); //main peak t0AC
  ttree->Branch("peakT0res",&peakT0res,"peakT0res/D"); // main peak of t0AC after fit
  ttree->Branch("spreadT0res",&spreadT0res,"spreadT0res/D"); //spread of main peak of t0AC after fit
  ttree->Branch("peakT0resErr",&peakT0resErr,"peakT0resErr/D"); // main peak of t0AC after fit
  ttree->Branch("spreadT0resErr",&spreadT0resErr,"spreadT0resErr/D"); //spread of main peak of t0AC after fit
  ttree->Branch("avT0fillRes",&avT0fillRes,"avT0fillRes/D"); //t0 fill res

  //save quantities for trending
  goodChannelRatio=(Double_t)GetGoodTOFChannelsRatio(runNumber,kFALSE,ocdbStorage);
	
  //--------------------------------- Multiplicity ----------------------------------//

  TH1F * hMulti = (TH1F*) generalList->FindObject("hTOFmatchedPerEvt");
  TH1F* hFractionEventsWhits = new TH1F("hFractionEventsWhits","hFractionEventsWhits;fraction of events with hits (%)",200,0.,100.);
  Float_t fraction=0.0;
  if (hMulti->GetEntries()>0.0) {
    fraction = ((Float_t) hMulti->GetBinContent(1))/((Float_t) hMulti->GetEntries());
    avMulti = hMulti->GetMean();
  } else fraction=0.0;
  hFractionEventsWhits->Fill(fraction*100.);
  
  //--------------------------------- T0F signal ----------------------------------//
  TH1F * hRawTime = (TH1F*)generalList->FindObject("hTOFmatchedESDrawTime");
  if ((hRawTime)&&(hRawTime->GetEntries()>0)){
    avRawTime=hRawTime->GetMean();
    if (!isMC){
      hRawTime->Fit("landau","RQ0","",200.,250.);
      if (hRawTime->GetFunction("landau")) {
	peakRawTime=(hRawTime->GetFunction("landau"))->GetParameter(1);
	spreadRawTime=(hRawTime->GetFunction("landau"))->GetParameter(2);
	peakRawTimeErr=(hRawTime->GetFunction("landau"))->GetParError(1);
	spreadRawTimeErr=(hRawTime->GetFunction("landau"))->GetParError(2);
      // 	printf("Main peak raw time (landau): mean = %f +- %f\n",peakTime,peakTimeErr );
      // printf("Main peak raw time (landau): spread = %f +- %f\n",spreadRawTime,spreadRawTimeErr );
      }
    } else {
      printf("Reminder: Raw time not available in MC simulated data.");
    }
  }
  MakeUpHisto(hRawTime, "matched tracks", 21, kGreen+2);
  hRawTime->Rebin(2);
  
  TH1F * hTime = (TH1F*)generalList->FindObject("hTOFmatchedESDtime");
  if ((hTime)&&(hTime->GetEntries()>0)) {
    avTime=hTime->GetMean();
    hTime->Fit("landau","RQ0","",0.,50.);
    if (hTime->GetFunction("landau")) {
      peakTime=(hTime->GetFunction("landau"))->GetParameter(1);
      spreadTime=(hTime->GetFunction("landau"))->GetParameter(2);
      peakTimeErr=(hTime->GetFunction("landau"))->GetParError(1);
      spreadTimeErr=(hTime->GetFunction("landau"))->GetParError(2);
      negTimeRatio=((Double_t)hTime->Integral(1,3)*100.)/((Double_t)hTime->Integral());
    // printf("Main peak time (landau): mean = %f +- %f\n",peakTime,peakTimeErr );
    // printf("Main peak time (landau): spread = %f +- %f\n",spreadTime,spreadTimeErr );
    // printf("Ratio of tracks with time<12.5 ns / total = %f\n",negTimeRatio );
    }
    MakeUpHisto(hTime, "matched tracks", 20, kBlue+2);
    hTime->Rebin(2);
    
    TLegend *lTime = new TLegend(0.7125881,0.6052519,0.979435,0.7408306,NULL,"brNDC");
    lTime->SetTextSize(0.04281433);
    lTime->AddEntry(hRawTime, "raw","L");
    lTime->AddEntry(hTime, "ESD","L"); 
    lTime->SetFillColor(kWhite);
    lTime->SetShadowColor(0);
  }
      
  TH1F * hTot = (TH1F*)generalList->FindObject("hTOFmatchedESDToT");
  if ((hTot)&&(hTot->GetEntries()>0)){
    avTot=hTot->GetMean();
    hTot->Fit("gaus","","",0.,50.);
    if (hTot->GetFunction("gaus")) {
      peakTot=(hTot->GetFunction("gaus"))->GetParameter(1);
      spreadTot=(hTot->GetFunction("gaus"))->GetParameter(2);
      peakTotErr=(hTot->GetFunction("gaus"))->GetParError(1);
      spreadTotErr=(hTot->GetFunction("gaus"))->GetParError(2);
    // printf("Main peak ToT (gaus): mean = %f +- %f\n",peakTot,peakTotErr );
    // printf("Main peak ToT (gaus): spread = %f +- %f\n",spreadTot,spreadTotErr );	
    }
  }      
  MakeUpHisto(hTot, "matched tracks", 8, kViolet-3);
  
  char orphansTxt[200];
  if (hTot->GetEntries()>1){
    orphansRatio=((Float_t) hTot->GetBinContent(1))/((Float_t) hTot->GetEntries()) ;
  }
  sprintf(orphansTxt,"orphans/matched = %4.2f%%",orphansRatio*100.);
  TPaveText *tOrphans = new TPaveText(0.38,0.63,0.88,0.7, "NDC");
  tOrphans->SetBorderSize(0);
  tOrphans->SetTextSize(0.045);
  tOrphans->SetFillColor(0); //white background
  tOrphans->SetTextAlign(12);
  tOrphans->SetTextColor(kViolet-3);
  tOrphans->AddText(orphansTxt);
  
  TH1F * hL=(TH1F*)generalList->FindObject("hTOFmatchedESDtrkLength");
  char negLengthTxt[200];
  if (hL->GetEntries()>0){
    avL=hL->GetMean();
    negLratio=(hL->Integral(1,750))/((Float_t) hL->GetEntries()) ;
  }
  MakeUpHisto(hL, "matched tracks", 1, kBlue+2);
  sprintf(negLengthTxt,"trk with L<350cm /matched = %4.2f%%", negLratio*100.);
  TPaveText *tLength = new TPaveText(0.15,0.83,0.65,0.87, "NDC");
  tLength->SetBorderSize(0);
  tLength->SetTextSize(0.04);
  tLength->SetFillColor(0); //white background
  tLength->SetTextAlign(11);
  tLength->SetTextColor(kOrange-3);
  tLength->AddText(negLengthTxt);

  //--------------------------------- residuals -------------------------------------//
  TH2F* hDxPos4profile = (TH2F*) generalList->FindObject("hTOFmatchedDxVsPtPos");
  TH2F* hDxNeg4profile = (TH2F*) generalList->FindObject("hTOFmatchedDxVsPtNeg");    
  const Int_t ybinMin = 0;
  const Int_t ybinMax = hDxPos4profile->GetYaxis()->GetNbins() ;
  TProfile * profDxPos = (TProfile*)hDxPos4profile->ProfileX("profDxPos", ybinMin,ybinMax); 
  profDxPos->SetLineWidth(2);
  TProfile * profDxNeg = (TProfile*)hDxNeg4profile->ProfileX("profDxNeg", ybinMin, ybinMax); 
  profDxNeg->SetLineWidth(2);  
 
  TH1 *profRatioPosOverNegDx = (TH1*) profDxPos->Clone();
  profRatioPosOverNegDx->SetName("profRatioPosOverNegDx");
  profRatioPosOverNegDx->Divide((TH1*) profDxNeg);
  profRatioPosOverNegDx->GetYaxis()->SetRangeUser(-5.,5.);
  profRatioPosOverNegDx->GetXaxis()->SetRangeUser(0.,2.);

  TH2F* hTOFmatchedDzVsStrip = (TH2F*)generalList->FindObject("hTOFmatchedDzVsStrip");

  // fout->cd();
  // hTOFmatchedDzVsStrip->Write();
  // hDxPos4profile->Write();
  // hDxNeg4profile->Write();
  // profDxPos->Write();
  // profDxNeg->Write();
  // profRatioPosOverNegDx->Write();
 
  //--------------------------------- matching eff ----------------------------------//
  //matching as function of pT
  TH1F * hMatchingVsPt = new TH1F("hMatchingVsPt","Matching probability vs. Pt; Pt(GeV/c); matching probability", 50, 0., 5. );
  TH1F * hDenom = (TH1F*)generalList->FindObject("hESDprimaryTrackPt"); 
  if (hDenom) {  
    hMatchingVsPt=(TH1F*) generalList->FindObject("hTOFmatchedESDPt")->Clone(); 
    hMatchingVsPt->Rebin(5);
    hDenom->Rebin(5);
    hMatchingVsPt->Divide(hDenom);
    hMatchingVsPt->GetYaxis()->SetTitle("matching efficiency");
    hMatchingVsPt->SetTitle("TOF matching efficiency as function of transverse momentum");
    hMatchingVsPt->GetYaxis()->SetRangeUser(0,1.2); 
  }

  if (hMatchingVsPt->GetEntries()>0){
    hMatchingVsPt->Fit("pol0","","",1.0,10.);
    hMatchingVsPt->Draw();
    if (hMatchingVsPt->GetFunction("pol0")){
      matchEffLinFit1Gev=(hMatchingVsPt->GetFunction("pol0"))->GetParameter(0);
      matchEffLinFit1GevErr=(hMatchingVsPt->GetFunction("pol0"))->GetParError(0);	
      //printf("Matching efficiency fit param is %f +- %f\n",matchEffLinFit1Gev,matchEffLinFit1GevErr );
    }
  } else {
    printf("WARNING: matching efficiency plot has 0 entries. Skipped!\n");
  }
  MakeUpHisto(hMatchingVsPt, "efficiency", 1, kBlue+2);

  //matching as function of eta
  TH1F * hMatchingVsEta =new TH1F("hMatchingVsEta","Matching probability vs. #\Eta; #\Eta; matching probability", 20, -1., 1.);
  hDenom->Clear();
  hDenom=(TH1F*)generalList->FindObject("hTOFprimaryESDeta"); 
  if (hDenom) {  
    hMatchingVsEta=(TH1F*) generalList->FindObject("hTOFmatchedESDeta")->Clone(); 
    hMatchingVsEta->Rebin(5);
    hDenom->Rebin(5);
    hMatchingVsEta->Divide(hDenom);
    hMatchingVsEta->GetXaxis()->SetRangeUser(-1,1);
    hMatchingVsEta->GetYaxis()->SetTitle("matching efficiency");
    hMatchingVsEta->GetYaxis()->SetRangeUser(0,1.2);
    hMatchingVsEta->SetTitle("TOF matching efficiency as function of pseudorapidity");
  }
  MakeUpHisto(hMatchingVsEta, "efficiency", 1, kBlue+2);

  
  //matching as function of phi
  TH1F * hMatchingVsPhi = new TH1F("hMatchingVsPhi","Matching probability vs. Phi; Phi(rad); matching probability", 628, 0., 6.28);
  hDenom->Clear();
  hDenom=(TH1F*)generalList->FindObject("hTOFprimaryESDphi");  
  if (hDenom) {  
    hMatchingVsPhi=(TH1F*) generalList->FindObject("hTOFmatchedESDphi")->Clone(); 
    hMatchingVsPhi->Rebin(2);
    hDenom->Rebin(2);    
    hMatchingVsPhi->Divide(hDenom);
    hMatchingVsPhi->GetYaxis()->SetTitle("matching efficiency");
    hMatchingVsPhi->SetTitle("TOF matching efficiency as function of phi");
    hMatchingVsPhi->GetYaxis()->SetRangeUser(0,1.2);
  }
  MakeUpHisto(hMatchingVsPhi, "efficiency", 1, kBlue+2);

  // fout->cd();
  // hMatchingVsPt->Write();
  // hMatchingVsEta->Write();
  // hMatchingVsPhi->Write();
  

  //--------------------------------- t-texp ----------------------------------//
  TH2F * hBetaP=(TH2F*)pidList->FindObject("hTOFmatchedESDpVsBeta");
  if (hBetaP) hBetaP->GetYaxis()->SetRangeUser(0.,1.2);
  
  TH1F * hMass=(TH1F*)pidList->FindObject("hTOFmatchedMass");
  MakeUpHisto(hMass, "tracks", 1, kBlue+2);
  // hMass->SetFillColor(kAzure+10);
  // hMass->SetFillStyle(1001);
  hMass->Rebin(2);
  
  //pions
  TH2F * hDiffTimeT0TOFPion1GeV=(TH2F*)pidList->FindObject("hTOFmatchedTimePion1GeV"); 
  TH1F * hPionDiff=(TH1F*)pidList->FindObject("hTOFmatchedExpTimePi"); 
  TH1F * hKaonDiff=(TH1F*)pidList->FindObject("hTOFmatchedExpTimeKa"); 
  TH1F * hProtonDiff=(TH1F*)pidList->FindObject("hTOFmatchedExpTimePro"); 
  if ((hPionDiff)&&(hPionDiff->GetEntries()>0)) {
    avPiDiffTime=hPionDiff->GetMean();
    hPionDiff->Fit("gaus","","",-1000.,500.);
    if (hPionDiff->GetFunction("gaus")){
      peakPiDiffTime=(hPionDiff->GetFunction("gaus"))->GetParameter(1);
      spreadPiDiffTime=(hPionDiff->GetFunction("gaus"))->GetParameter(2);
      peakPiDiffTimeErr=(hPionDiff->GetFunction("gaus"))->GetParError(1);
      spreadPiDiffTimeErr=(hPionDiff->GetFunction("gaus"))->GetParError(2);
      // printf("Main peak t-t_exp (gaus): mean = %f +- %f\n",peakPiDiffTime,peakPiDiffTimeErr );
      // printf("Main peak t-t_exp (gaus): spread = %f +- %f\n",spreadPiDiffTime,spreadPiDiffTimeErr );
    }
  }
  
  TH2F * hDiffTimePi=(TH2F*)pidList->FindObject("hTOFmatchedExpTimePiVsP"); 
  hDiffTimePi->GetYaxis()->SetRangeUser(-5000.,5000.);
  hDiffTimePi->SetTitle("PIONS t-t_{exp,#pi} (from tracking) vs. P");

  // TProfile * profDiffTimePi = (TProfile*)hDiffTimePi->ProfileX("profDiffTimePi", 490, 510); 
  // if (profDiffTimePi){
  //   profDiffTimePi->SetLineWidth(2);
  //   profDiffTimePi->SetLineColor(kRed+2); 
  // }

  //Kaon
  TH2F * hDiffTimeKa=(TH2F*)pidList->FindObject("hTOFmatchedExpTimeKaVsP");  
  hDiffTimeKa->SetTitle("KAONS t-t_{exp,K} (from tracking) vs. P");
  hDiffTimeKa->GetYaxis()->SetRangeUser(-5000.,5000.);
  
  // TProfile * profDiffTimeKa = (TProfile*)hDiffTimeKa->ProfileX("profDiffTimeKa", 490, 510); 
  // if (profDiffTimeKa) {
  //   profDiffTimeKa->SetLineWidth(2);
  //   profDiffTimeKa->SetLineColor(kBlue);  
  // }

  //Protons
  TH2F * hDiffTimePro=(TH2F*)pidList->FindObject("hTOFmatchedExpTimeProVsP"); 
  hDiffTimePro->SetTitle("PROTONS t-t_{exp,p} (from tracking) vs. P");
  hDiffTimePro->GetYaxis()->SetRangeUser(-5000.,5000.);
  
  // TProfile * profDiffTimePro = (TProfile*)hDiffTimePro->ProfileX("profDiffTimePro", 490, 510); 
  // if (profDiffTimePro) {
  //   profDiffTimePro->SetLineWidth(2);
  //   profDiffTimePro->SetLineColor(kGreen+2);  
  // }

  /*
  TH2F * hDiffTimePiTh=(TH2F*)pidList->FindObject("hTOFtheoreticalExpTimePiVsP");  
  hDiffTimePiTh->GetYaxis()->SetRangeUser(-2000.,2000.); 

  TProfile * profDiffTimePiTh = (TProfile*)hDiffTimePiTh->ProfileX("profDiffTimePiTh", 490, 510);
  if (profDiffTimePiTh) {
    profDiffTimePiTh->SetLineWidth(2);
    profDiffTimePiTh->SetLineColor(kRed+2);
  }

  TH2F * hDiffTimeKaTh=(TH2F*)pidList->FindObject("hTOFtheoreticalExpTimeKaVsP"); 
  hDiffTimeKaTh->GetYaxis()->SetRangeUser(-2000.,2000.);

  TProfile * profDiffTimeKaTh = (TProfile*)hDiffTimeKaTh->ProfileX("profDiffTimeKaTh", 490, 510); 
  if (profDiffTimeKaTh) {
    profDiffTimeKaTh->SetLineWidth(2);
    profDiffTimeKaTh->SetLineColor(kBlue);
  }

  TH2F * hDiffTimeProTh=(TH2F*)pidList->FindObject("hTOFtheoreticalExpTimePro"); 
  hDiffTimeProTh->GetYaxis()->SetRangeUser(-2000.,2000.);
 
  TProfile * profDiffTimeProTh = (TProfile*)hDiffTimeProTh->ProfileX("profDiffTimeProTh", 490, 510);
  if (profDiffTimeProTh) {
    profDiffTimeProTh->SetLineWidth(2);
    profDiffTimeProTh->SetLineColor(kGreen+2);
  }
  TLegend * lPid=new TLegend(0.75,0.75,0.95,0.95,"PID");
  lPid->AddEntry(profDiffTimePi,"#pi^{#pm}","l");
  lPid->AddEntry(profDiffTimeKa,"K^{#pm}","l");
  lPid->AddEntry(profDiffTimePro,"p^{#pm}","l");
  
     if (canvasE){
     TCanvas *cPidPerformanceTh= new TCanvas("cPidPerformanceTh","summary of pid performance - theoretical times",700,700);
     cPidPerformanceTh->Divide(2,2);
     cPidPerformanceTh->cd(1);
     gPad->SetLogz();
     gPad->SetGridx();
     gPad->SetGridy();
     hDiffTimePiTh->Draw("colz");
     profDiffTimePiTh->Draw("same");
     cPidPerformanceTh->cd(2);
     gPad->SetLogz();
     gPad->SetGridx();
     gPad->SetGridy();
     hDiffTimeKaTh->Draw("colz");
     profDiffTimeKaTh->Draw("same");
     cPidPerformanceTh->cd(3);
     gPad->SetLogz();
     gPad->SetGridx();
     gPad->SetGridy();
     hDiffTimeProTh->Draw("colz");
     profDiffTimeProTh->Draw("same");
     }
    fout->cd();
    hPionDiff->Write();
    hKaonDiff->Write();
    hProtonDiff->Write();
    hBetaP->Write();
    hMass->Write();
    hDiffTimeT0TOFPion1GeV->Write();
    hDiffTimePi->Write();
    profDiffTimePi->Write();
    hDiffTimeKa->Write();
    profDiffTimeKa->Write();
    hDiffTimePro->Write();
    profDiffTimePro->Write();
    //lPid->Draw();
    hDiffTimePiTh->Write();
    profDiffTimePiTh->Write();
    hDiffTimeKaTh->Write();
    profDiffTimeKaTh->Write();
    hDiffTimeProTh->Write();
    profDiffTimeProTh->Write();
  */

  TH2F * hSigmaPi=(TH2F*)pidList->FindObject("hTOFExpSigmaPi"); 
  hSigmaPi->GetYaxis()->SetRangeUser(-5.,5.);
  TProfile * profSigmaPi = (TProfile*)hSigmaPi->ProfileX("profSigmaPi"); 
  profSigmaPi->SetLineWidth(2);
  profSigmaPi->SetLineColor(kRed+2); 

  TH2F * hSigmaKa=(TH2F*)pidList->FindObject("hTOFExpSigmaKa"); 
  hSigmaKa->GetYaxis()->SetRangeUser(-5.,5.);
  TProfile * profSigmaKa = (TProfile*)hSigmaKa->ProfileX("profSigmaKa"); 
  profSigmaKa->SetLineWidth(2);
  profSigmaKa->SetLineColor(kBlue);  

  TH2F * hSigmaPro=(TH2F*)pidList->FindObject("hTOFExpSigmaPro"); 
  hSigmaPro->GetYaxis()->SetRangeUser(-5.,5.);
  TProfile * profSigmaPro = (TProfile*)hSigmaPro->ProfileX("profSigmaPro"); 
  profSigmaPro->SetLineWidth(2);
  profSigmaPro->SetLineColor(kGreen+2);  

  // fout->cd();
  // hSigmaPi->Write();
  // profSigmaPi->Write();
  // hSigmaKa->Write();
  // profSigmaKa->Write();
  // hSigmaPro->Write();
  // profSigmaPro->Write();
  
  //--------------------------------- T0 detector ----------------------------------//
	
  TH1F*hT0A=(TH1F*)timeZeroList->FindObject("hEventT0DetA");
  if ((hT0A)&&(hT0A->GetEntries()>0)) {
    avT0A = hT0A->GetMean();
    hT0A->Fit("gaus","RQ0", "", -1000., 1000.);
    if (hT0A->GetFunction("gaus")) {
      peakT0A=(hT0A->GetFunction("gaus"))->GetParameter(1);
      spreadT0A=(hT0A->GetFunction("gaus"))->GetParameter(2);
      peakT0AErr=(hT0A->GetFunction("gaus"))->GetParError(1);
      spreadT0AErr=(hT0A->GetFunction("gaus"))->GetParError(2);	
      // printf("Main peak T0A(gaus): mean = %f +- %f\n",peakT0A,peakT0AErr );
      // printf("Main peak T0A (gaus): spread = %f +- %f\n",spreadT0A,spreadT0AErr );
      //add integral of main peak over total
    }
  }
  MakeUpHisto(hT0A, "events", 8, kBlue);
  hT0A->Rebin(2);

  TH1F*hT0C=(TH1F*)timeZeroList->FindObject("hEventT0DetC");
  if ((hT0C)&&(hT0C->GetEntries()>0)) {
    avT0C=hT0C->GetMean();
    hT0C->Fit("gaus","RQ0","", -1000., 1000.);
    if (hT0C->GetFunction("gaus")) {
      peakT0C=(hT0C->GetFunction("gaus"))->GetParameter(1);
      spreadT0C=(hT0C->GetFunction("gaus"))->GetParameter(2);
      peakT0CErr=(hT0C->GetFunction("gaus"))->GetParError(1);
      spreadT0CErr=(hT0C->GetFunction("gaus"))->GetParError(2);	
      // printf("Main peak T0C(gaus): mean = %f +- %f\n",peakT0C,peakT0CErr );
      // printf("Main peak T0C (gaus): spread = %f +- %f\n",spreadT0C,spreadT0CErr );
      //add integral of main peak over total
    }
  }
  MakeUpHisto(hT0C, "events", 8, kGreen+1);
  hT0C->Rebin(2);
	
  TH1F*hT0AC=(TH1F*)timeZeroList->FindObject("hEventT0DetAND");
  if ((hT0AC)&&(hT0AC->GetEntries()>0)) {
    avT0AC=hT0AC->GetMean();
    hT0AC->Fit("gaus","RQ0", "",-1000., 1000.);
    if (hT0AC->GetFunction("gaus")) {
      peakT0AC=(hT0AC->GetFunction("gaus"))->GetParameter(1);
      spreadT0AC=(hT0AC->GetFunction("gaus"))->GetParameter(2);
      peakT0ACErr=(hT0AC->GetFunction("gaus"))->GetParError(1);
      spreadT0ACErr=(hT0AC->GetFunction("gaus"))->GetParError(2);	
      // printf("Main peak T0AC(gaus): mean = %f +- %f\n",peakT0AC,peakT0ACErr );
      // printf("Main peak T0AC (gaus): spread = %f +- %f\n",spreadT0AC,spreadT0ACErr );	 
    }
  }
  MakeUpHisto(hT0AC, "events", 8, kRed+1);
  hT0AC->Rebin(2);
  
  TLegend *lT0 = new TLegend(0.7125881,0.6052519,0.979435,0.7408306,NULL,"brNDC");
  lT0->SetTextSize(0.041);
  lT0->AddEntry(hT0AC, "T0 A&C","L");
  lT0->AddEntry(hT0A, "T0 A","L"); 
  lT0->AddEntry(hT0C, "T0 C","L");
  lT0->SetFillColor(kWhite);
  lT0->SetShadowColor(0);
  
  TH1F*hT0res=(TH1F*)timeZeroList->FindObject("hT0DetRes");
  if ((hT0res)&&(hT0res->GetEntries()>0)) {
    avT0res=hT0res->GetMean();
    hT0res->Fit("gaus");
    if (hT0res->GetFunction("gaus")) {
      peakT0res=(hT0res->GetFunction("gaus"))->GetParameter(1);
      spreadT0res=(hT0res->GetFunction("gaus"))->GetParameter(2);
      peakT0resErr=(hT0res->GetFunction("gaus"))->GetParError(1);
      spreadT0resErr=(hT0res->GetFunction("gaus"))->GetParError(2);	
      // printf("Main peak T0res(gaus): mean = %f +- %f\n",peakT0res,peakT0resErr );
      // printf("Main peak T0res (gaus): spread = %f +- %f\n",spreadT0res,spreadT0resErr );	 
    //add integral of main peak over total
    }
  }
  TH1F*hT0fillRes=(TH1F*)timeZeroList->FindObject("hT0fillRes");
  if ((hT0fillRes)&&(hT0fillRes->GetEntries()>0)) {
    avT0fillRes=hT0fillRes->GetMean();
  }
  
  // fout->cd();
  // hT0AC->Write();
  // hT0A->Write();
  // hT0C->Write();
  // hT0res->Write();
  // hT0fillRes->Write();
  // lT0->Write();
 	
  //Fill tree and save to file
  ttree->Fill();
  printf("==============  Saving trending quantities in tree for run %i ===============\n",runNumber);
  trendFile->cd();
  ttree->Write();
  trendFile->Close();
  
  if (canvasE){
    // TString plotDir(Form("Plots_run%d",runNumber));
    // gSystem->Exec(Form("mkdir %s",plotDir.Data()));
    TString plotDir(".");

    TCanvas *cTrackProperties= new TCanvas("cTrackProperties","summary of matched tracks properties", 1200, 500);
    cTrackProperties->Divide(3,1);
    cTrackProperties->cd(1);
    gPad->SetLogy();
    hTime->Draw("");
    hRawTime ->Draw("same");
    lTime->Draw();  
    cTrackProperties->cd(2);
    gPad->SetLogy();
    hTot->Draw("");
    tOrphans->Draw(); 
    cTrackProperties->cd(3);
    gPad->SetLogy();
    hL->Draw("");
    tLength->Draw(); 

    TCanvas *cResiduals= new TCanvas("residuals","residuals", 900,500);
    cResiduals->Divide(2,1);
    cResiduals->cd(1);
    gPad->SetLogz();
    hDxPos4profile->GetYaxis()->SetRangeUser(-5.,5.);
    hDxPos4profile->Draw("colz");
    profDxPos->SetLineColor(kRed);
    profDxPos ->Draw("same");
    cResiduals->cd(2);
    gPad->SetLogz();
    hDxNeg4profile->GetYaxis()->SetRangeUser(-5.,5.); 
    hDxNeg4profile->Draw("colz");
    profDxNeg->SetLineColor(kBlue);
    profDxNeg->Draw("same"); 

    TCanvas* cProfile = new TCanvas("cProfile","cProfile",50,50, 750,550);
    cProfile->cd();
    gPad->SetLogz();
    hTOFmatchedDzVsStrip->Draw("colz");
    Int_t binmin = hTOFmatchedDzVsStrip->GetYaxis()->FindBin(-3);
    Int_t binmax = hTOFmatchedDzVsStrip->GetYaxis()->FindBin(3);
    TProfile* hDzProfile = (TProfile*)hTOFmatchedDzVsStrip->ProfileX("hDzProfile",binmin, binmax);
    hDzProfile->SetLineWidth(3);
    hDzProfile->Draw("same");
    
    TCanvas *cMatchingPerformance= new TCanvas("cMatchingPerformance","summary of matching performance",1200,500);
    cMatchingPerformance->Divide(3,1);
    cMatchingPerformance->cd(1);
    hMatchingVsPt->Draw();
    cMatchingPerformance->cd(2);
    hMatchingVsEta->Draw();
    cMatchingPerformance->cd(3);
    hMatchingVsPhi->Draw();
    
    TCanvas *cPidPerformance= new TCanvas("cPidPerformance","summary of pid performance", 900,500);
    cPidPerformance->Divide(2,1);
    cPidPerformance->cd(1);
    gPad->SetLogz();
    hBetaP->Draw("colz");   
    cPidPerformance->cd(2);
    gPad->SetLogy();
    hMass->Draw("HIST ");

    TCanvas *cPidPerformance2= new TCanvas("cPidPerformance2","summary of pid performance - expected times", 1200, 500);
    cPidPerformance2->Divide(3,1);
    cPidPerformance2->cd(1);
    gPad->SetLogz();
    hDiffTimePi->Draw("colz");
    //profDiffTimePi->Draw("same");
    cPidPerformance2->cd(2);
    gPad->SetLogz();
    hDiffTimeKa->Draw("colz");
    //profDiffTimeKa->Draw("same");
    cPidPerformance2->cd(3);
    gPad->SetLogz();
    hDiffTimePro->Draw("colz");
    //profDiffTimePro->Draw("same");
  
    TLegend * lSigmaPid=new TLegend(0.75,0.75,0.95,0.95,"#sigma_{PID}");
    lSigmaPid->AddEntry(profSigmaPi,"#pi^{#pm}","l");
    lSigmaPid->AddEntry(profSigmaKa,"K^{#pm}","l");
    lSigmaPid->AddEntry(profSigmaPro,"p^{#pm}","l");	  
    TCanvas *cPidPerformance3= new TCanvas("cPidPerformance3","summary of pid performance - sigmas",1200,500);
    cPidPerformance3->Divide(3,1);
    cPidPerformance3->cd(1);
    gPad->SetLogz();
    hSigmaPi->Draw("colz");
    profSigmaPi->Draw("same");
    cPidPerformance3->cd(2);
    gPad->SetLogz();
    hSigmaKa->Draw("colz");
    profSigmaKa->Draw("same");
    cPidPerformance3->cd(3);
    gPad->SetLogz();
    hSigmaPro->Draw("colz");
    profSigmaPro->Draw("same");
  
    TCanvas *cT0detector= new TCanvas("cT0detector","T0 detector",800,600);
    cT0detector->Divide(2,1);
    cT0detector->cd(1);
    gPad->SetGridx();
    hT0AC->Draw("");
    hT0AC->SetTitle("timeZero measured by T0 detector");
    hT0A ->Draw("same");
    hT0C ->Draw("same");
    lT0->Draw();  
    cT0detector->cd(2);
    hT0res->Draw();
    
    cPidPerformance3->Print(Form("%s/%i_PID_sigmas.png", plotDir.Data(), runNumber));
    cPidPerformance->Print(Form("%s/%i_PID.png",plotDir.Data(), runNumber));
    //cPidPerformanceTh->Print(Form("%s/PID_theoreticalTimes.png",plotDir.Data()));
    cPidPerformance2->Print(Form("%s/%i_PID_ExpTimes.png",plotDir.Data(), runNumber));
    cMatchingPerformance->Print(Form("%s/%i_Matching.png",plotDir.Data(), runNumber));
    cTrackProperties->Print(Form("%s/%i_TrackProperties.png",plotDir.Data(), runNumber));
    cResiduals->Print(Form("%s/%i_Residuals.png",plotDir.Data(), runNumber));
    cProfile->Print(Form("%s/%i_ProfileDZvsStripNumber.png",plotDir.Data(), runNumber));
    cT0detector->Print(Form("%s/%i_T0Detector.png",plotDir.Data(), runNumber));
  }
 
  return  0;
}


//----------------------------------------------------------
Double_t GetGoodTOFChannelsRatio(Int_t run = -1, Bool_t saveMap = kFALSE, TString OCDBstorage = "raw://")
{
  /*
    It retrieves from OCDB the number of good (= efficient && not noisy && HW ok) TOF channels.
    Optionally is saves the channel map
  */
  if (run<=0) {
    printf("MakeTrendingTOFqa.C - ERROR in CheckCalibStatus(): invalid run number. Please set a run number.\n"); 
    return 0.0;
  }
  
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(OCDBstorage.Data());
  cdb->SetRun(run);
  
  AliCDBEntry *cdbe = cdb->Get("TOF/Calib/Status");
  if (!cdbe) {
    printf("MakeTrendingTOFqa.C - ERROR in CheckCalibStatus(): OCDB entry not available. Please, try again.\n");
    return 0.0;
  }  

  AliTOFChannelOnlineStatusArray *array = (AliTOFChannelOnlineStatusArray *)cdbe->GetObject();
  TH2F *hOkMap = new TH2F("hOkMap", "Ok map (!noisy & !problematic & efficient);sector;strip", 72, 0., 18., 91, 0., 91.);

  AliTOFcalibHisto calibHisto;
  calibHisto.LoadCalibHisto();
  AliTOFcalib calib;
  calib.Init();
  Int_t sector, sectorStrip, padx, fea;
  Float_t hitmapx, hitmapy;
  for (Int_t i = 0; i <  array->GetSize(); i++) {
    sector = calibHisto.GetCalibMap(AliTOFcalibHisto::kSector, i);
    sectorStrip = calibHisto.GetCalibMap(AliTOFcalibHisto::kSectorStrip, i);
    padx = calibHisto.GetCalibMap(AliTOFcalibHisto::kPadX, i);
    fea = padx / 12;
    hitmapx = sector + ((Double_t)(3 - fea) + 0.5) / 4.;
    hitmapy = sectorStrip;
    if ( !(array->GetNoiseStatus(i) == AliTOFChannelOnlineStatusArray::kTOFNoiseBad)   &&
	 (calib.IsChannelEnabled(i,kTRUE,kTRUE)))
      hOkMap->Fill(hitmapx,hitmapy);
  }
  Int_t nOk=(Int_t) hOkMap->GetEntries();
  Double_t ratioOk=nOk/152928.;
  if (saveMap) hOkMap->SaveAs(Form("run%i_OKChannelsMap.root",run));
  cout << "###    Run " << run << ": TOF channels ok = " << nOk << "/ total 152928 channels = " << ratioOk*100. << "% of whole TOF" << endl;
  return ratioOk;
}

//----------------------------------------------------------
void MakeUpHisto(TH1* histo, TString titleY, Int_t marker=20, Color_t color=kBlue+2)
{
  if (!histo) return;
  histo->SetMarkerStyle(marker);
  histo->SetMarkerSize(0.7);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->SetFillColor(kWhite);
  histo->SetFillStyle(0); 
  histo->GetYaxis()->SetTitle(titleY.Data());
  histo->GetYaxis()->SetTitleOffset(1.35);
  histo->GetXaxis()->SetLabelSize(0.03);
  
  return;
}
