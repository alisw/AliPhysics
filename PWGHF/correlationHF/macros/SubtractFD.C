TString localcode=gSystem->ExpandPathName("$ALICE_PHYSICS/../src/PWGHF/correlationHF/macros/");
TString strsystem="none";
Bool_t fitcodeIsLoaded=kFALSE;
void SubtractMCclosureModulation(TH1D *h, Double_t ptD, Double_t ptTrmin, Double_t ptTrmax);
void SetSystemStringForTemplateFDnames(TString str){
  strsystem=str;
}
void SetLocalCodeDir(TString dirloc){
  localcode=dirloc;
}
//______________________________________________________________________________________________
double v2modulation(double* x, double* p)
// function for v2 modulation
{
	double mod;
	mod = p[0] *(1.+ 2*p[1]*p[2]*TMath::Cos(2*x[0]));
	return mod;
}
//______________________________________________________________________________________________
TF1 * v2function(Double_t pedestal, Double_t v2D, Double_t v2had, Int_t colour =1){
    
    TString name = "v2modulfunc";
    name += Form("ped%.2f_D%.2f_Had%.2f",pedestal,v2D,v2had);
    
    TF1* fv2modulation = new TF1(name.Data(), v2modulation, -0.5*TMath::Pi(), 1.5 *TMath::Pi(),3);
    fv2modulation->SetParameter(0,pedestal);
    fv2modulation->SetParameter(1,v2D);
    fv2modulation->SetParameter(2,v2had);
    fv2modulation->SetLineColor(colour);
    
    return fv2modulation;
}



//______________________________________________________________________________________________
TH1D * subtractpedestal(TH1D *histo,Double_t pedestal){
    TString nameoutput = histo->GetName();
    nameoutput += "_subtr_";
    nameoutput += "pedestal";
    
    
    TH1D * outputhisto = (TH1D*)histo->Clone(nameoutput.Data());
    outputhisto->Reset();
    Double_t value = 0;
    
    for(Int_t iBin = 1; iBin <= histo->GetNbinsX();iBin++){
        
        
        value = histo->GetBinContent(iBin);
        value -= pedestal;
        outputhisto->SetBinContent(iBin,value);
        
        outputhisto->SetBinError(iBin,histo->GetBinError(iBin));
        
    }
    
    //outputhisto->SetMarkerColor(color);
    //outputhisto->SetLineColor(color);
    //outputhisto->SetMarkerStyle(19+color);
    
    return outputhisto;
}





//______________________________________________________________________________________________
void ModulatePedestalWithv2 (TH1D * histo, TF1 *v2){
    TString nameoutput = histo->GetName();
    nameoutput += "_subtr_";
    nameoutput += v2->GetName();
    
    //cout << "crahs 1 " << endl;
    
    TH1D * v2modulation = (TH1D*)histo->Clone(nameoutput.Data());
    v2modulation->Reset();
    // cout << "crahs 2 " << endl;
    Double_t value = 0;
    
    for(Int_t iBin = 1; iBin <= histo->GetNbinsX();iBin++){
        
        
        //value = 0; //histo->GetBinContent(iBin);
        value = v2->Eval(histo->GetBinCenter(iBin));
        
        v2modulation->SetBinContent(iBin,value);
        v2modulation->SetBinError(iBin,histo->GetBinError(iBin));
    }
    
    histo->Add(v2modulation);
    //histo = outputhisto;
    
   // return outputhisto;
}
//______________________________________________________________________________________________
void GetTemplateFromFit(TH1D *h,TH1D *hOut,TString strCanv="cFit",Int_t methodFD=0, Double_t v2D = 0, Double_t v2had = 0){

    // modulation should be fine now... try to cross check and plot the outputs :)
    
    cout << "Processing histo " << hOut->GetName() << endl;
    cout << "Applying v2 values: D " << v2D << ", had " << v2had << endl;
    
    // Fit Template from B
    TCanvas *cFit=new TCanvas(strCanv.Data(),strCanv.Data(),700,700);
    cFit->cd();
    h->Draw();
  
    if(methodFD==0){
      TF1 *fitFunction=FitPlotsShort(h,2,0,0);
      for(Int_t j=1;j<=h->GetNbinsX();j++){
	hOut->SetBinContent(j,fitFunction->Eval(hOut->GetBinCenter(j)));
	hOut->SetBinError(j,0);
      }
      hOut->SetLineColor(kViolet);
      hOut->Draw("same");
      
    }
    else if(methodFD==1){
      TCanvas * test = new TCanvas("test","test",0,0,1000,800);
      test->Divide(2,1);
      
      Double_t nsybc, ensybc,asybc, easybc,ebase;
      Double_t nsybcMC, ensybcMC,asybcMC, easybcMC,ebaseMC;
      cFit->cd();
      TF1 *fitFunctionMC=FitPlots(h,2,-2,3,nsybcMC, ensybcMC,asybcMC, easybcMC,kFALSE);
      Double_t baseMC=GetBaseline(ebaseMC);//fitFunctionMC->GetParameter(0);
      
      TCanvas *cFitData=new TCanvas(Form("Data%s",strCanv.Data()),"cGetTemplatefromFit_FitData",700,700);
      cFitData->cd();
    
      TH1D *hDataForFit=(TH1D*)hOut->Clone("hDataForFit");
      hDataForFit->Draw("same");
      
      test->cd(1);
      hDataForFit->Draw("ep");

      TF1 *fitFunctionData=FitPlots(hDataForFit,1,5,3,nsybc, ensybc,asybc, easybc,kFALSE);
      Double_t baseData=GetBaseline(ebase);//fitFunctionData->GetParameter(0);
      
      // generate v2 function
      TF1* v2mod = v2function(baseData, v2D, v2had,2);
      

      for(Int_t j=1;j<=h->GetNbinsX();j++){
	//hOut->SetBinContent(j,fitFunctionMC->Eval(hOut->GetBinCenter(j))+baseData-baseMC);
        hOut->SetBinContent(j,fitFunctionMC->Eval(hOut->GetBinCenter(j))-baseMC);
	hOut->SetBinError(j,0);
      }
      
      // TH1D *clone = (TH1D*)hOut->Clone("clone");
      
      //  TH1D * clone2 = subtractpedestal(clone,-1.*baseData);
      
      // clone2->SetMarkerColor(2);
      ModulatePedestalWithv2(hOut,v2mod);
      // test->cd(2);
      // hOut->Draw("ep");
      // clone2->Draw("sameep");
      // v2mod->Draw("same");
      
      cFit->cd();
      hOut->SetLineColor(kViolet);
      hOut->DrawClone("same");
    
  }

  return;
}

//______________________________________________________________________________________________
void SubtractFDexploitingClassDzero(Double_t ptmin, Double_t ptmax, Double_t ptassoc, Double_t ptassocMax,
                                    TString strfileout="output.root",Int_t methodSubtr=2,Double_t purity=0.963,
                                    Int_t rebin=1,TString correlationDataFile="1Dcorr.root",
                                    TString spectraMacroOutput="HFPtSpectrum_Nb.root",
                                    TString strdirTempl="temp/",
                                    Int_t system =0, Double_t v2D = 0, Double_t v2Had = 0,Int_t systoption = 3, Int_t oldnames=1, Bool_t subtrMCclos=kFALSE){
    
    TFile *fDataCorr=TFile::Open(correlationDataFile.Data(),"READ");
    TCanvas *cData;
    TH1D *hData;
    if(oldnames) {
		cData=(TCanvas*)fDataCorr->Get("c3");
		hData=(TH1D*)cData->FindObject("hsubtract_norm");
	} else {
		cData=(TCanvas*)fDataCorr->Get(Form("cFinal_%1.1fto%1.1f",ptassoc,ptassocMax));
		hData=(TH1D*)cData->FindObject("h1D_SubtrNorm");
	} 

    TFile *fSpectrum=TFile::Open(spectraMacroOutput.Data(),"READ");
    TString namefpr="gFcConservative";
    if(system==2) namefpr="gFPromptCombined";
    TGraphAsymmErrors *gr=(TGraphAsymmErrors*)fSpectrum->Get(namefpr.Data());

    const Int_t Dzero = 0;
    SubtractFDexploitingClass(Dzero,hData,gr,ptmin,ptmax,ptassoc,ptassocMax,strdirTempl,strfileout,purity,methodSubtr,1,system,v2D,v2Had,systoption,subtrMCclos);
    fDataCorr->Close();
    fSpectrum->Close();
}


//______________________________________________________________________________________________
void SubtractFDexploitingClassDstar(Double_t ptmin, Double_t ptmax, Double_t ptassoc, Double_t ptassocMax,
                                    TString strfileout="output.root",Int_t methodSubtr=2,Double_t purity=0.963,
                                    Int_t rebin=1,TString correlationDataFile="1Dcorr.root",
                                    TString spectraMacroOutput="HFPtSpectrum_Nb.root",
                                    TString strdirTempl="temp/",
                                    Int_t system =0, Double_t v2D = 0, Double_t v2Had = 0,Int_t systoption = 3, Int_t oldnames=1, Bool_t subtrMCclos=kFALSE){
    
    TFile *fDataCorr=TFile::Open(correlationDataFile.Data(),"READ");
    TCanvas *cData;
    TH1D *hData;
    if(oldnames) {
		cData=(TCanvas*)fDataCorr->Get("FinalDphiCorrelationsCanvas");
		hData=(TH1D*)cData->FindObject("DMesonHadronDPhi");
	} else {
		cData=(TCanvas*)fDataCorr->Get(Form("cFinal_%1.1fto%1.1f",ptassoc,ptassocMax));
		hData=(TH1D*)cData->FindObject("h1D_SubtrNorm");
	} 
    
    TFile *fSpectrum=TFile::Open(spectraMacroOutput.Data(),"READ");
    TString namefpr="gFcConservative";
    if(system==2) namefpr="gFPromptCombined";
    TGraphAsymmErrors *gr=(TGraphAsymmErrors*)fSpectrum->Get(namefpr.Data());
    
    const Int_t Dstar = 1;
    SubtractFDexploitingClass(Dstar,hData,gr,ptmin,ptmax,ptassoc,ptassocMax,strdirTempl,strfileout,purity,methodSubtr,1,system,v2D,v2Had,systoption,subtrMCclos);
    fDataCorr->Close();
    fSpectrum->Close();    
}

//______________________________________________________________________________________________
void SubtractFDexploitingClassDplus(Double_t ptmin, Double_t ptmax, Double_t ptassoc, Double_t ptassocMax,
                                    TString strfileout="output.root",Int_t methodSubtr=2,Double_t purity=0.963,
                                    Int_t rebin=1,TString correlationDataFile="1Dcorr.root",
                                    TString spectraMacroOutput="HFPtSpectrum_Nb.root",
                                    TString strdirTempl="temp/",
                                    Int_t system =0, Double_t v2D = 0, Double_t v2Had = 0,Int_t systoption = 3, Int_t oldnames=1, Bool_t subtrMCclos=kFALSE){
    
    TFile *fDataCorr=TFile::Open(correlationDataFile.Data(),"READ");
    TCanvas *cData;
    TH1D *hData;
    if(oldnames) {
		cData=(TCanvas*)fDataCorr->Get("corrCorrected1D");
		hData=(TH1D*)cData->FindObject("ME_Corrected_1D");
	} else {
		cData=(TCanvas*)fDataCorr->Get(Form("cFinal_%1.1fto%1.1f",ptassoc,ptassocMax));
		hData=(TH1D*)cData->FindObject("h1D_SubtrNorm");
	} 
    
    TFile *fSpectrum=TFile::Open(spectraMacroOutput.Data(),"READ");
    TString namefpr="gFcConservative";
    if(system==2) namefpr="gFPromptCombined";
    TGraphAsymmErrors *gr=(TGraphAsymmErrors*)fSpectrum->Get(namefpr.Data());

    const Int_t Dplus = 2;
    SubtractFDexploitingClass(Dplus,hData,gr,ptmin,ptmax,ptassoc,ptassocMax,strdirTempl,strfileout,purity,methodSubtr,1,system,v2D,v2Had,systoption,subtrMCclos);
    fDataCorr->Close();
    fSpectrum->Close();
}

//______________________________________________________________________________________________
void SubtractFDexploitingClassDzerov2Modulations(Double_t ptmin, Double_t ptmax, Double_t ptassoc,Double_t ptassocmax,
                                                 TString strfileout="output.root",Int_t methodSubtr=2,Double_t purity=0.965,
                                                 Int_t rebin=1,TString correlationDataFile="1Dcorr.root",
                                                 TString spectraMacroOutput="HFPtSpectrum_Nb.root", TString strdirTempl="temp",
                                                 Int_t system = 1, Int_t systoption=3, Int_t oldnames=1, Bool_t subtrMCclos=kFALSE){
    
    /*Parameter
     (Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],
     outputfilename.Data(),2,purities[ihadpt],
     1,inputcorrelation.Data(),inputfc.Data(),
     templatedir.Data(),1,systmode); // collsyst is always set to 1!
     */
    
    TFile *fDataCorr=TFile::Open(correlationDataFile.Data(),"READ");
    TCanvas *cData;
    TH1D *hData;
    if(oldnames) {
		cData=(TCanvas*)fDataCorr->Get("c3");
		hData=(TH1D*)cData->FindObject("hsubtract_norm");
	} else {
		cData=(TCanvas*)fDataCorr->Get(Form("cFinal_%1.1fto%1.1f",ptassoc,ptassocmax));
		hData=(TH1D*)cData->FindObject("h1D_SubtrNorm");
	} 
    
    TFile *fSpectrum=TFile::Open(spectraMacroOutput.Data(),"READ");
    TString namefpr="gFcConservative";
    if(system==2) namefpr="gFPromptCombined";
    TGraphAsymmErrors *gr=(TGraphAsymmErrors*)fSpectrum->Get(namefpr.Data());
    
    const Int_t Dzero = 0;
    SubtractFDexploitingClassv2Modulations(Dzero,hData,gr,ptmin,ptmax,ptassoc,ptassocmax,strdirTempl,strfileout,purity,methodSubtr,1,system,systoption,subtrMCclos);
}


//______________________________________________________________________________________________
void SubtractFDexploitingClassDstarv2Modulations(Double_t ptmin,Double_t ptmax,Double_t ptassoc,Double_t ptassocmax,
                                                 TString strfileout="output.root",Int_t methodSubtr=2,Double_t purity=0.965,
                                                 Int_t rebin=1,TString correlationDataFile="1Dcorr.root",
                                                 TString spectraMacroOutput="HFPtSpectrum_Nb.root", TString strdirTempl="temp",
                                                 Int_t system = 1, Int_t systoption=3, Int_t oldnames=1, Bool_t subtrMCclos=kFALSE){
    
    /*Parameter
     (Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],
     outputfilename.Data(),2,purities[ihadpt],
     1,inputcorrelation.Data(),inputfc.Data(),
     templatedir.Data(),1,systmode); // collsyst is always set to 1!
     */
    
    TFile *fDataCorr=TFile::Open(correlationDataFile.Data(),"READ");
    TCanvas *cData;
    TH1D *hData;
    if(oldnames) {
		cData=(TCanvas*)fDataCorr->Get("FinalDphiCorrelationsCanvas");
		hData=(TH1D*)cData->FindObject("DMesonHadronDPhi");
	} else {
		cData=(TCanvas*)fDataCorr->Get(Form("cFinal_%1.1fto%1.1f",ptassoc,ptassocmax));
		hData=(TH1D*)cData->FindObject("h1D_SubtrNorm");
	} 
    
    TFile *fSpectrum=TFile::Open(spectraMacroOutput.Data(),"READ");
    TString namefpr="gFcConservative";
    if(system==2) namefpr="gFPromptCombined";
    TGraphAsymmErrors *gr=(TGraphAsymmErrors*)fSpectrum->Get(namefpr.Data());
    
    const Int_t Dstar = 1;
    SubtractFDexploitingClassv2Modulations(Dstar,hData,gr,ptmin,ptmax,ptassoc,ptassocmax,strdirTempl,strfileout,purity,methodSubtr,1,system,systoption,subtrMCclos);
    fDataCorr->Close();
    fSpectrum->Close();
}

//______________________________________________________________________________________________
void SubtractFDexploitingClassDplusv2Modulations(Double_t ptmin,Double_t ptmax,Double_t ptassoc,Double_t ptassocmax,
                                                 TString strfileout="output.root",Int_t methodSubtr=2,Double_t purity=0.965,
                                                 Int_t rebin=1,TString correlationDataFile="1Dcorr.root",
                                                 TString spectraMacroOutput="HFPtSpectrum_Nb.root", TString strdirTempl="temp",
                                                 Int_t system = 1, Int_t systoption=3, Int_t oldnames=1, Bool_t subtrMCclos=kFALSE){
    
    /*Parameter 
     (Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],
     outputfilename.Data(),2,purities[ihadpt],
     1,inputcorrelation.Data(),inputfc.Data(),
     templatedir.Data(),1,systmode); // collsyst is always set to 1!
    */
    
    TFile *fDataCorr=TFile::Open(correlationDataFile.Data(),"READ");
    TCanvas *cData;
    TH1D *hData;
    if(oldnames) {
		cData=(TCanvas*)fDataCorr->Get("corrCorrected1D");
		hData=(TH1D*)cData->FindObject("ME_Corrected_1D");
	} else {
		cData=(TCanvas*)fDataCorr->Get(Form("cFinal_%1.1fto%1.1f",ptassoc,ptassocmax));
		hData=(TH1D*)cData->FindObject("h1D_SubtrNorm");
	} 
    
    TFile *fSpectrum=TFile::Open(spectraMacroOutput.Data(),"READ");
    TString namefpr="gFcConservative";
    if(system==2) namefpr="gFPromptCombined";
    TGraphAsymmErrors *gr=(TGraphAsymmErrors*)fSpectrum->Get(namefpr.Data());

    const Int_t Dplus = 2;
    SubtractFDexploitingClassv2Modulations(Dplus,hData,gr,ptmin,ptmax,ptassoc,ptassocmax,strdirTempl,strfileout,purity,methodSubtr,1,system,systoption,subtrMCclos);
    fDataCorr->Close();
    fSpectrum->Close();
}

//______________________________________________________________________________________________
void SubtractFDexploitingClass(Int_t meson,TH1D *hData,TGraphAsymmErrors *grFprompt,Double_t ptmin,Double_t ptmax,Double_t ptassoc,Double_t ptassocMax,TString strdirTempl,TString strfileout="FDoutput",Double_t purity=0.963,Int_t methodSubtr=2,Int_t rebin=1, Int_t system =0 /*0 is pp, 1 is pPb*/, Double_t v2D = 0, Double_t v2Had = 0, Int_t systoption = 3, Bool_t subtrMCclos = kFALSE){
    
    if(system ==0){
        if(TMath::Abs(v2D)>0.00000001 || TMath::Abs(v2Had)>0.00000001){
            cout << "Sorry, it is pp - let's keep assuming that in this system v2 is still 0 :) " << endl;
            cout << "Resetting v2 D meson and v2 hadron to 0" << endl;
        v2D = 0; v2Had = 0;
        }
    }
    cout << " " << endl;
    if(meson==0)cout << "---Dzero---  " << endl;
    else if(meson==1)cout << "---Dstar---  " << endl;
    else if(meson==2)cout << "---Dplus---  " << endl;
    else cout << "---No Meson---  " << endl;
    cout << "======== Inserted v2 modulation " << endl;
    cout << "for D meson = " << v2D << endl;
    cout << "for hadron = " << v2Had << endl;
    cout << " " << endl;
    
    if(!fitcodeIsLoaded){
      gROOT->LoadMacro(Form("%s/FitPlots.C",localcode.Data()));
      fitcodeIsLoaded=kTRUE;
    }
    hData->Sumw2();
    hData->Rebin(rebin);
    hData->Scale(purity*1./(Double_t)rebin);

    //Apply modulation from MC closure test (as correction on data points! pPb2016 preliminary approach!)
    if(subtrMCclos) SubtractMCclosureModulation(hData,0.5*(ptmin+ptmax),ptassoc,ptassocMax);

    TString strmes;
    if(meson==0){
        strmes="Dzero";
    }
    else if(meson==1){
        strmes="Dstar";//"Dstar";
    }
    else if(meson==2){
        strmes="Dplus";
    }
    else return;

    strfileout += Form("_v2D%.2f_v2had%.2f.root",v2Had,v2D);
    cout << "Saving outptut as " << strfileout << endl;
    cout << "checkmate 1 in SubtractFDexploitingClass" << endl;
   
    TFile *fOut=new TFile(strfileout.Data(),"RECREATE");
    fOut->cd();
    
    AliHFCorrelationFDsubtraction *fdSubtracter=new AliHFCorrelationFDsubtraction();
    fdSubtracter->SetUncorrectedHistogram(hData);
    fdSubtracter->SetDptRange(ptmin,ptmax);
    fdSubtracter->SetMethod(1);
    fdSubtracter->SetNTemplates(3);
    fdSubtracter->SetSystOption(systoption);// default one
    fdSubtracter->Init();
    fdSubtracter->SetFpromptGraphFc(grFprompt);
    
    
    TString templPromptPerugia0="", templPromptPerugia2010="", templPromptPerugia2011="";
    TString templPerugia0="", templPerugia2010="",templPerugia2011="";
    
    TString strsystnametempl;
    if(system ==0){
      if(strsystem.EqualTo("none")){
	strsystnametempl="pp";
      }
      else {
	strsystnametempl=strsystem;
      }
        cout << "Templates for on pp... Loading Perugia0, Perugia2010, Perugia2011" << endl;
        templPromptPerugia0=Form("%s/%sCorrelationPlotsPerugia0Pt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);
        templPerugia0=Form("%s/%sCorrelationPlotsPerugia0Pt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);

        templPromptPerugia2010=Form("%s/%sCorrelationPlotsPerugia2010Pt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);
        templPerugia2010=Form("%s/%sCorrelationPlotsPerugia2010Pt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);

        templPromptPerugia2011=Form("%s/%sCorrelationPlotsPerugia2011Pt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);
        templPerugia2011=Form("%s/%sCorrelationPlotsPerugia2011Pt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);
    }
    
    else if(system ==1){
      if(strsystem.EqualTo("none")){
	strsystnametempl="pPb";
      }
      else {
	strsystnametempl=strsystem;
      }
        cout << "Templates for pPb... Loading Perugia0, Perugia2010, Perugia2011" << endl;
        templPromptPerugia0=Form("%s/%sCorrelationPlotsPerugia0wBoostPt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);
        templPerugia0=Form("%s/%sCorrelationPlotsPerugia0wBoostPt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);

        templPromptPerugia2010=Form("%s/%sCorrelationPlotsPerugia2010wBoostPt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);
        templPerugia2010=Form("%s/%sCorrelationPlotsPerugia2010wBoostPt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);

        templPromptPerugia2011=Form("%s/%sCorrelationPlotsPerugia2011wBoostPt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);
        templPerugia2011=Form("%s/%sCorrelationPlotsPerugia2011wBoostPt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);
    }

     else if(system ==2){
      if(strsystem.EqualTo("none")){
  strsystnametempl="pPb";
      }
      else {
  strsystnametempl=strsystem;
      }
        cout << "Templates for pPb... Loading Perugia2011, Perugia2010, PYTHIA8" << endl;
        templPromptPerugia0=Form("%s/%sCorrelationPlotsPerugia2011wBoostPt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);
        templPerugia0=Form("%s/%sCorrelationPlotsPerugia2011wBoostPt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);

        templPromptPerugia2010=Form("%s/%sCorrelationPlotsPerugia2010wBoostPt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);
        templPerugia2010=Form("%s/%sCorrelationPlotsPerugia2010wBoostPt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);

        templPromptPerugia2011=Form("%s/%sCorrelationPlotsPYTHIA8wBoostPt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);
        templPerugia2011=Form("%s/%sCorrelationPlotsPYTHIA8wBoostPt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocMax);
    }   
    
    TH1D **hFDtempl=new TH1D*[3];
    TH1D **hPromptTempl=new TH1D*[3];
    TH1D **hFDtemplFile=new TH1D*[3];
    
    //1. Perugia0
    TFile *fFDtemplCorr=TFile::Open(templPerugia0.Data(),"READ");
    fOut->cd();
    hFDtemplFile[0]=(TH1D*)(fFDtemplCorr->Get("hCorrDeltaPhi")->Clone());
    hFDtemplFile[0]->Draw();
    TFile *fFDtemplPromptCorr=TFile::Open(templPromptPerugia0.Data(),"READ");
    fOut->cd();
    hPromptTempl[0]=(TH1D*)(fFDtemplPromptCorr->Get("hCorrDeltaPhi")->Clone());
    hPromptTempl[0]->SetTitle("PromptPerugia0");
    hPromptTempl[0]->SetName("PromptPerugia0");
    fFDtemplCorr->Close();
    fFDtemplPromptCorr->Close();

    //2. Perugia 2010
    TFile *fFDtemplCorr=TFile::Open(templPerugia2010.Data(),"READ");
    fOut->cd();
    hFDtemplFile[1]=(TH1D*)(fFDtemplCorr->Get("hCorrDeltaPhi")->Clone());
    TFile *fFDtemplPromptCorr=TFile::Open(templPromptPerugia2010.Data(),"READ");
    fOut->cd();
    hPromptTempl[1]=(TH1D*)(fFDtemplPromptCorr->Get("hCorrDeltaPhi")->Clone());
    hPromptTempl[1]->SetTitle("PromptPerugia2010");
    hPromptTempl[1]->SetName("PromptPerugia2010");
    fFDtemplCorr->Close();
    fFDtemplPromptCorr->Close();

    //3. Perugia 2011
    TFile *fFDtemplCorr=TFile::Open(templPerugia2011.Data(),"READ");
    fOut->cd();
    hFDtemplFile[2]=(TH1D*)(fFDtemplCorr->Get("hCorrDeltaPhi")->Clone());
    TFile *fFDtemplPromptCorr=TFile::Open(templPromptPerugia2011.Data(),"READ");
    fOut->cd();
    hPromptTempl[2]=(TH1D*)(fFDtemplPromptCorr->Get("hCorrDeltaPhi")->Clone());
    hPromptTempl[2]->SetTitle("PromptPerugia2011");
    hPromptTempl[2]->SetName("PromptPerugia2011");
    fFDtemplCorr->Close();
    fFDtemplPromptCorr->Close();


    Printf("The pointer to the objects are: %p , %p ",    hPromptTempl[0],    hFDtemplFile[0]);
    hPromptTempl[0]->GetTitle();
    hFDtemplFile[0]->GetTitle();
    Printf("And they are not empty");

    cout << "crash here 1" << endl;
    cout<<"Macro path works?"<< localcode.Data() <<endl;
    Printf(localcode.Data());
    TH1D * cloner = NULL;
    
    if(methodSubtr==0)
        hFDtempl[0]=(TH1D*)hFDtemplFile[0]->Clone("hTemplDfromB");
    else if(methodSubtr==1){
        hFDtempl[0]=(TH1D*)hData->Clone("hTemplDfromBFitFunc");
        cloner = (TH1D*)hFDtempl[0]->Clone("clonetest");
        GetTemplateFromFit(hFDtemplFile[0],hFDtempl[0],"cFitFD",0,v2D,v2Had);
    }
    else if(TMath::Abs(methodSubtr)==2){
        hFDtempl[0]=(TH1D*)hData->Clone("hTemplDfromBFitFunc");
        cloner = (TH1D*)hFDtempl[0]->Clone("clonetest");
        GetTemplateFromFit(hFDtemplFile[0],hFDtempl[0],"cFitFD",1,v2D,v2Had);
    }
    if(!hFDtempl[0])return;
    fdSubtracter->AddTemplateHisto(hFDtempl[0]);
    cloner->SetMarkerColor(2);
    
    TCanvas * check = new TCanvas("check","check",0,0,1000,800);
    check->cd();
    cloner->Draw("ep");
    cout << "test 2" << endl;
    hFDtempl[0]->Draw("ep");
    
   
    //2. Perugia 2010
    if(methodSubtr==0)
        hFDtempl[1]=(TH1D*)hFDtemplFile[1]->Clone("hTemplDfromB");
    else if(methodSubtr==1){
        hFDtempl[1]=(TH1D*)hData->Clone("hTemplDfromBFitFunc");
        GetTemplateFromFit(hFDtemplFile[1],hFDtempl[1],"cFitFD",0,v2D,v2Had);
    }
    else if(TMath::Abs(methodSubtr)==2){// Increase the baseline of MC template in order to match the baseline observed in data
        hFDtempl[1]=(TH1D*)hData->Clone("hTemplDfromBFitFunc");
        GetTemplateFromFit(hFDtemplFile[1],hFDtempl[1],"cFitFD",1,v2D,v2Had);
    }
    if(!hFDtempl[1])return;
    cout << "test 1" << endl;
    fdSubtracter->AddTemplateHisto(hFDtempl[1]);
    
  

    //3. Perugia 2011
    if(methodSubtr==0)
        hFDtempl[2]=(TH1D*)hFDtemplFile[2]->Clone("hTemplDfromB");
    else if(methodSubtr==1){
        hFDtempl[2]=(TH1D*)hData->Clone("hTemplDfromBFitFunc");
        GetTemplateFromFit(hFDtemplFile[2],hFDtempl[2],"cFitFD",0,v2D,v2Had);
    }
    else if(TMath::Abs(methodSubtr)==2){// Increase the baseline of MC template in order to match the baseline observed in data
        hFDtempl[2]=(TH1D*)hData->Clone("hTemplDfromBFitFunc");
        GetTemplateFromFit(hFDtemplFile[2],hFDtempl[2],"cFitFD",1,v2D,v2Had);
    }
    if(!hFDtempl[2])return;
    fdSubtracter->AddTemplateHisto(hFDtempl[2]);
    
    TH1D * template0 = fdSubtracter->GetTemplate(0);
    TH1D * template1 = fdSubtracter->GetTemplate(1);
    TH1D * template2 = fdSubtracter->GetTemplate(2);
    
    // NOW CALCULATE THE ENVELOPE
    fdSubtracter->CalculateEnvelope();
    
    TCanvas *cHistEnvelope=new TCanvas("cHistEnvelope","cHistEnvelope",700,700);
    cHistEnvelope->cd();
    TH1D *hEnvMin=fdSubtracter->GetHistoEnvelopeMin();
    TH1D *hEnvMax=fdSubtracter->GetHistoEnvelopeMax();
    hEnvMin->Draw();
    hEnvMax->Draw("same");
    
    TCanvas *cEnvelope=new TCanvas("cEnvelope","cEnvelope",700,700);
    cEnvelope->cd();
    TGraphAsymmErrors *grEnv=fdSubtracter->GetGraphEnvelope();
    grEnv->Draw("ap");
    
    TCanvas *cHistEnvelopeRatio=new TCanvas("cHistEnvelopeRatio","cHistEnvelopeRatio",700,700);
    cHistEnvelopeRatio->cd();
    TH1D *hEnvMinRatio=fdSubtracter->GetHistoEnvelopeRatioMin();
    TH1D *hEnvMaxRatio=fdSubtracter->GetHistoEnvelopeRatioMax();
    hEnvMinRatio->Draw();
    hEnvMaxRatio->Draw("same");
    
    TCanvas *cGrEnvelopeRatio=new TCanvas("cGrEnvelopeRatio","cGrEnvelopeRatio",700,700);
    cGrEnvelopeRatio->cd();
    TGraphAsymmErrors *grEnvRatio=fdSubtracter->GetGraphEnvelopeRatio();
    grEnvRatio->Draw("ap");
    
    TCanvas *cfinalPlot=new TCanvas("cFinal","cFinal",700,700);
    cfinalPlot->cd();
    TH1D *hFinal=fdSubtracter->GetCentralSubtractedPlot();

    AliHFDhadronCorrSystUnc *oUnc=new AliHFDhadronCorrSystUnc();
    oUnc->SetName("SystematicUncertainty");
    if(system == 0)oUnc->InitStandardUncertaintiesPP2010(meson,0.5*(ptmin+ptmax),ptassoc); // load pp uncertainties
    if(system == 1)oUnc->InitStandardUncertaintiesPPb2013(meson,0.5*(ptmin+ptmax),ptassoc); // load p-Pb uncertainties
    if(system == 2)oUnc->InitStandardUncertaintiesPPb2016(meson,0.5*(ptmin+ptmax),ptassoc,ptassocMax); // load p-Pb uncertainties
    oUnc->SetHistoBeautyFDmin(fdSubtracter->GetHistoRelSystUncMin(),"",kTRUE);
    oUnc->SetHistoBeautyFDmax(fdSubtracter->GetHistoRelSystUncMax(),"",kTRUE);
    
    TCanvas *canvFinal=oUnc->BuildSystUncertaintyPlotVsDeltaPhi(hFinal,1);
    canvFinal->cd();
    hPromptTempl[0]->Draw("same");
    hPromptTempl[0]->SetLineColor(kBlue);
    hPromptTempl[1]->Draw("same");
    hPromptTempl[1]->SetLineColor(kGreen);
    hPromptTempl[2]->Draw("same");
    hPromptTempl[2]->SetLineColor(kOrange);
    
    
    TCanvas *cErrorNonFlatOnly=new TCanvas("cErrorNonFlatOnly","cErrorNonFlatOnly",700,700);
    cErrorNonFlatOnly->cd();
    TGraphAsymmErrors *grNonFlat=oUnc->GetTotNonFlatUncGraph();
    hFinal->Draw();
    grNonFlat->Draw("E2");
    
     
    oUnc->Write();
    canvFinal->Write();
    hFinal->Write();
    cfinalPlot->Write();
    hEnvMin->Write();
    hEnvMax->Write();
    template0->Write();
    template1->Write();
    template2->Write();
    fOut->Close();
    cout << "checkmate 2" << endl;
    
}

//______________________________________________________________________________________________
void SubtractFDexploitingClassv2Modulations(const Int_t meson,TH1D *hData,TGraphAsymmErrors *grFprompt,Double_t ptmin,Double_t ptmax,Double_t ptassoc,Double_t ptassocmax,TString strdirTempl,TString strfileout="FDoutput",Double_t purity=0.963,Int_t methodSubtr=2,Int_t rebin=1, Int_t system =0, Int_t systoption =3, Bool_t subtrMCclos = kFALSE){
    
    cout << "This might take a while - go and get a coffee, you deserved it :) " << endl;
    gSystem->Sleep(2000);
    
    if(system ==0){
        cout << "Sorry, it is pp : Going Back: Exiting...." << endl;
        return;
    }
    
    Double_t v2hadmin, v2hadmax, v2Dmin, v2Dmax;
    v2Dmin = 0.05; v2Dmax = 0.13;
    
    if(ptassoc < 0.6){v2hadmin = 0.05; v2hadmax = 0.07;}
    if(ptassoc > 0.6){v2hadmin = 0.08; v2hadmax = 0.13;}
    cout << " " << endl;
    cout << "======== Inserted v2 modulation " << endl;
    if(meson==0)cout << "---Dzero---  " << endl;
    else if(meson==1)cout << "---Dstar---  " << endl;
    else if(meson==2)cout << "---Dplus---  " << endl;
    else cout << "---Dplus---  " << endl;

    cout << "for D meson = " << v2Dmin << "," << v2Dmax <<  endl;
    cout << "for hadron = " << v2hadmin << "," << v2hadmax << endl;
    cout << " " << endl;
    
    if(!fitcodeIsLoaded){
      gROOT->LoadMacro(Form("%s/FitPlots.C",localcode.Data()));
    }
    hData->Sumw2();
    hData->Rebin(rebin);
    hData->Scale(purity*1./(Double_t)rebin);

    //Apply modulation from MC closure test (as correction on data points! pPb2016 preliminary approach!)
    if(subtrMCclos) SubtractMCclosureModulation(hData,0.5*(ptmin+ptmax),ptassoc,ptassocmax);

    TString strmes;
    if(meson==0)strmes="Dzero";
    else if(meson==1)strmes="Dstar";//"D0";//"Dstar";
    else if(meson==2)strmes="Dplus";
    else return;
    
    TFile *fOut=new TFile(strfileout.Data(),"RECREATE");
    fOut->cd();

    AliHFCorrelationFDsubtraction *fdSubtracter=new AliHFCorrelationFDsubtraction();
    fdSubtracter->SetUncorrectedHistogram(hData);
    fdSubtracter->SetDptRange(ptmin,ptmax);
    fdSubtracter->SetMethod(1);
    fdSubtracter->SetNTemplates(15);
    fdSubtracter->SetSystOption(systoption);// default one
    fdSubtracter->Init();
    fdSubtracter->SetFpromptGraphFc(grFprompt);
    

    TString templPromptPerugia0="", templPromptPerugia2010="", templPromptPerugia2011="";
    TString templPerugia0="", templPerugia2010="", templPerugia2011="";
    TString strsystnametempl;

    if(system ==1){
        cout << "Templates for pPb... Loading Perugia0, Perugia2010, Perugia2011" << endl;
      if(strsystem.EqualTo("none")){
	strsystnametempl="pPb";
      }
      else {
	strsystnametempl=strsystem;
      }
        templPromptPerugia0=Form("%s/%sCorrelationPlotsPerugia0wBoostPt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocmax);
        templPromptPerugia2010=Form("%s/%sCorrelationPlotsPerugia2010wBoostPt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocmax);
        templPromptPerugia2011=Form("%s/%sCorrelationPlotsPerugia2011wBoostPt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocmax);
        templPerugia0=Form("%s/%sCorrelationPlotsPerugia0wBoostPt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocmax);
        templPerugia2010=Form("%s/%sCorrelationPlotsPerugia2010wBoostPt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocmax);
        templPerugia2011=Form("%s/%sCorrelationPlotsPerugia2011wBoostPt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocmax);
    }
    if(system ==2){
        cout << "Templates for pPb... Loading PYTHIA8, Perugia2011, Perugia2010" << endl;
      if(strsystem.EqualTo("none")){
        strsystnametempl="pPb";
      }
      else {
        strsystnametempl=strsystem;
      }
        templPromptPerugia0=Form("%s/%sCorrelationPlotsPYTHIA8wBoostPt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocmax);
        templPromptPerugia2010=Form("%s/%sCorrelationPlotsPerugia2011wBoostPt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocmax);
        templPromptPerugia2011=Form("%s/%sCorrelationPlotsPerugia2010wBoostPt%sfromC%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocmax);
        templPerugia0=Form("%s/%sCorrelationPlotsPYTHIA8wBoostPt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocmax);
        templPerugia2010=Form("%s/%sCorrelationPlotsPerugia2011wBoostPt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocmax);
        templPerugia2011=Form("%s/%sCorrelationPlotsPerugia2010wBoostPt%sfromB%.0fTo%.0f_ptAssall%.1fto%.1f_DeltaEta10.root",strdirTempl.Data(),strsystnametempl.Data(),strmes.Data(),ptmin,ptmax,ptassoc,ptassocmax);
    }    
    else cout<< "Not pPb - Where is your files ?? " << endl;

    TH1D **hFDtempl=new TH1D*[15];
    TH1D **hPromptTempl=new TH1D*[15];
    TH1D **hFDtemplFile=new TH1D*[3];
    
    //1. Load Perugia0
    TFile *fFDtemplCorr=TFile::Open(templPerugia0.Data(),"READ");
    fOut->cd();
    hFDtemplFile[0]=(TH1D*)(fFDtemplCorr->Get("hCorrDeltaPhi")->Clone());
    hFDtemplFile[0]->Draw();
    TFile *fFDtemplPromptCorr=TFile::Open(templPromptPerugia0.Data(),"READ");
    fOut->cd();
    hPromptTempl[0]=(TH1D*)(fFDtemplPromptCorr->Get("hCorrDeltaPhi")->Clone());
    hPromptTempl[0]->SetTitle("PromptPerugia0");
    hPromptTempl[0]->SetName("PromptPerugia0");
    fFDtemplCorr->Close();
    fFDtemplPromptCorr->Close();
    if(methodSubtr==0) {
        hFDtempl[0]=(TH1D*)hFDtemplFile[0]->Clone("hTemplDfromB0"); // check those numbers
        hFDtempl[1]=(TH1D*)hFDtemplFile[0]->Clone("hTemplDfromB1");
        hFDtempl[2]=(TH1D*)hFDtemplFile[0]->Clone("hTemplDfromB2");
        hFDtempl[3]=(TH1D*)hFDtemplFile[0]->Clone("hTemplDfromB3");
        hFDtempl[4]=(TH1D*)hFDtemplFile[0]->Clone("hTemplDfromB4");
    }
    else if(methodSubtr==1 || TMath::Abs(methodSubtr)==2){
        hFDtempl[0]=(TH1D*)hData->Clone("hTemplDfromBFitFunc0");
        hFDtempl[1]=(TH1D*)hData->Clone("hTemplDfromBFitFunc1");
        hFDtempl[2]=(TH1D*)hData->Clone("hTemplDfromBFitFunc2");
        hFDtempl[3]=(TH1D*)hData->Clone("hTemplDfromBFitFunc3");
        hFDtempl[4]=(TH1D*)hData->Clone("hTemplDfromBFitFunc4");
    }
    
    
    //load perugia 2010
    TFile *fFDtemplCorr=TFile::Open(templPerugia2010.Data(),"READ");
    fOut->cd();
    hFDtemplFile[1]=(TH1D*)(fFDtemplCorr->Get("hCorrDeltaPhi")->Clone());
    TFile *fFDtemplPromptCorr=TFile::Open(templPromptPerugia2010.Data(),"READ");
    fOut->cd();
    hPromptTempl[1]=(TH1D*)(fFDtemplPromptCorr->Get("hCorrDeltaPhi")->Clone());
    hPromptTempl[1]->SetTitle("PromptPerugia2010");
    hPromptTempl[1]->SetName("PromptPerugia2010");
    fFDtemplCorr->Close();
    fFDtemplPromptCorr->Close();
    if(methodSubtr==0) {
        hFDtempl[5]=(TH1D*)hFDtemplFile[1]->Clone("hTemplDfromB5"); // check those numbers
        hFDtempl[6]=(TH1D*)hFDtemplFile[1]->Clone("hTemplDfromB6");
        hFDtempl[7]=(TH1D*)hFDtemplFile[1]->Clone("hTemplDfromB7");
        hFDtempl[8]=(TH1D*)hFDtemplFile[1]->Clone("hTemplDfromB8");
        hFDtempl[9]=(TH1D*)hFDtemplFile[1]->Clone("hTemplDfromB9");
    }
    else if(methodSubtr==1 || TMath::Abs(methodSubtr)==2){
        hFDtempl[5]=(TH1D*)hData->Clone("hTemplDfromBFitFunc5");
        hFDtempl[6]=(TH1D*)hData->Clone("hTemplDfromBFitFunc6");
        hFDtempl[7]=(TH1D*)hData->Clone("hTemplDfromBFitFunc7");
        hFDtempl[8]=(TH1D*)hData->Clone("hTemplDfromBFitFunc8");
        hFDtempl[9]=(TH1D*)hData->Clone("hTemplDfromBFitFunc9");
    }
    
    
    //load perugia 2011
    TFile *fFDtemplCorr=TFile::Open(templPerugia2011.Data(),"READ");
    fOut->cd();
    hFDtemplFile[2]=(TH1D*)(fFDtemplCorr->Get("hCorrDeltaPhi")->Clone());
    TFile *fFDtemplPromptCorr=TFile::Open(templPromptPerugia2011.Data(),"READ");
    fOut->cd();
    hPromptTempl[2]=(TH1D*)(fFDtemplPromptCorr->Get("hCorrDeltaPhi")->Clone());
    hPromptTempl[2]->SetTitle("PromptPerugia2011");
    hPromptTempl[2]->SetName("PromptPerugia2011");
    fFDtemplCorr->Close();
    fFDtemplPromptCorr->Close();
    if(methodSubtr==0){
        hFDtempl[10]=(TH1D*)hFDtemplFile[2]->Clone("hTemplDfromB10"); // check those numbers
        hFDtempl[11]=(TH1D*)hFDtemplFile[2]->Clone("hTemplDfromB11");
        hFDtempl[12]=(TH1D*)hFDtemplFile[2]->Clone("hTemplDfromB12");
        hFDtempl[13]=(TH1D*)hFDtemplFile[2]->Clone("hTemplDfromB13");
        hFDtempl[14]=(TH1D*)hFDtemplFile[2]->Clone("hTemplDfromB14");
    }
    else if(methodSubtr==1 || TMath::Abs(methodSubtr)==2){
        hFDtempl[10]=(TH1D*)hData->Clone("hTemplDfromBFitFunc10");
        hFDtempl[11]=(TH1D*)hData->Clone("hTemplDfromBFitFunc11");
        hFDtempl[12]=(TH1D*)hData->Clone("hTemplDfromBFitFunc12");
        hFDtempl[13]=(TH1D*)hData->Clone("hTemplDfromBFitFunc13");
        hFDtempl[14]=(TH1D*)hData->Clone("hTemplDfromBFitFunc14");
    }
    
    
    
    
    Int_t mode = -1.;
    if(methodSubtr==1) mode = 0;
    if(methodSubtr==2) mode = 1;
    for(Int_t k = 0; k<3; k++){
        // looping on different generators
      GetTemplateFromFit(hFDtemplFile[k],hFDtempl[k*5],"cFitFD",mode,0,0); // no v2
      GetTemplateFromFit(hFDtemplFile[k],hFDtempl[k*5+1],"cFitFD",mode,v2hadmax,v2Dmax); // maximum values of v2
      GetTemplateFromFit(hFDtemplFile[k],hFDtempl[k*5+2],"cFitFD",mode,v2hadmin,v2Dmax); //  values of v2
      GetTemplateFromFit(hFDtemplFile[k],hFDtempl[k*5+3],"cFitFD",mode,v2hadmax,v2Dmin); //  values of v2
      GetTemplateFromFit(hFDtemplFile[k],hFDtempl[k*5+4],"cFitFD",mode,v2hadmin,v2Dmin); //  values of v2

    }
 
 
    //adding templates to subtractor
    for(Int_t i = 0; i<15; i++){
      fdSubtracter->AddTemplateHisto(hFDtempl[i]);
    }

    
    // NOW CALCULATE THE ENVELOPE
    fdSubtracter->CalculateEnvelope();
    
    
    TCanvas *cHistEnvelope=new TCanvas("cHistEnvelope","cHistEnvelope",700,700);
    cHistEnvelope->cd();
    TH1D *hEnvMin=fdSubtracter->GetHistoEnvelopeMin();
    TH1D *hEnvMax=fdSubtracter->GetHistoEnvelopeMax();
    hEnvMin->Draw();
    hEnvMax->Draw("same");
    
    TCanvas *cEnvelope=new TCanvas("cEnvelope","cEnvelope",700,700);
    cEnvelope->cd();
    TGraphAsymmErrors *grEnv=fdSubtracter->GetGraphEnvelope();
    grEnv->Draw("ap");
    
    
    TCanvas *cHistEnvelopeRatio=new TCanvas("cHistEnvelopeRatio","cHistEnvelopeRatio",700,700);
    cHistEnvelopeRatio->cd();
    TH1D *hEnvMinRatio=fdSubtracter->GetHistoEnvelopeRatioMin();
    TH1D *hEnvMaxRatio=fdSubtracter->GetHistoEnvelopeRatioMax();
    hEnvMinRatio->Draw();
    hEnvMaxRatio->Draw("same");
    
    
    TCanvas *cGrEnvelopeRatio=new TCanvas("cGrEnvelopeRatio","cGrEnvelopeRatio",700,700);
    cGrEnvelopeRatio->cd();
    TGraphAsymmErrors *grEnvRatio=fdSubtracter->GetGraphEnvelopeRatio();
    grEnvRatio->Draw("ap");
    
    
    TCanvas *cfinalPlot=new TCanvas("cFinal","cFinal",700,700);
    cfinalPlot->cd();
    TH1D *hFinal=fdSubtracter->GetCentralSubtractedPlot();
    
    AliHFDhadronCorrSystUnc *oUnc=new AliHFDhadronCorrSystUnc();
    oUnc->SetName("SystematicUncertainty");
    if(system ==0) oUnc->InitStandardUncertainties2010(meson,0.5*(ptmin+ptmax),ptassoc); // load pp uncertainties
    if(system ==1) oUnc->InitStandardUncertaintiesPPb2013(meson,0.5*(ptmin+ptmax),ptassoc); // load p-Pb uncertainties
    if(system ==2) oUnc->InitStandardUncertaintiesPPb2016(meson,0.5*(ptmin+ptmax),ptassoc,ptassocmax); // load p-Pb uncertainties
    oUnc->SetHistoBeautyFDmin(fdSubtracter->GetHistoRelSystUncMin(),"",kTRUE);
    oUnc->SetHistoBeautyFDmax(fdSubtracter->GetHistoRelSystUncMax(),"",kTRUE);
    
    TCanvas *canvFinal=oUnc->BuildSystUncertaintyPlotVsDeltaPhi(hFinal,1);
    canvFinal->cd();
    hPromptTempl[0]->Draw("same");
    hPromptTempl[0]->SetLineColor(kBlue);
    hPromptTempl[1]->Draw("same");
    hPromptTempl[1]->SetLineColor(kGreen);
    hPromptTempl[2]->Draw("same");
    hPromptTempl[2]->SetLineColor(kOrange);
    
    
    TCanvas *cErrorNonFlatOnly=new TCanvas("cErrorNonFlatOnly","cErrorNonFlatOnly",700,700);
    cErrorNonFlatOnly->cd();
    TGraphAsymmErrors *grNonFlat=oUnc->GetTotNonFlatUncGraph();
    hFinal->Draw();
    grNonFlat->Draw("E2");
    
    
 //   strfileout += Form("_v2D%.2f_v2had%.2f.root",v2Had,v2D);
    
    cout << "Saving outptut as " << strfileout << endl;
    
    
    cout << "check :) 1" << endl;
    oUnc->Write();
    canvFinal->Write();
    hFinal->Write();
    cfinalPlot->Write();
    hEnvMin->Write();
    hEnvMax->Write();
    for(Int_t k =0; k<15; k++){
        (fdSubtracter->GetTemplate(k))->Write();
    }
   
    cout << "check :) 2" << endl;
    //fdSubtracter->Write();
    fOut->Close();
     cout << "check :) 3" << endl;
   
}

//______________________________________________________________________________________________
void OpenOutputFileAndDraw(TString strfile,Double_t ptminD,Double_t ptmaxD,TString strMeson="D^{0}",Double_t ptminAss=0.3, Double_t pTmaxAss=1.0, Double_t deltaeta=1, Int_t system=0,Double_t max=10,Bool_t rangeLabel=kFALSE){
  gStyle->SetOptStat(0000);
  if(!fitcodeIsLoaded){
    gROOT->LoadMacro(Form("%s/FitPlots.C",localcode.Data()));
  }

  TFile *f=TFile::Open(strfile.Data(),"READ");
  AliHFDhadronCorrSystUnc *syst=(AliHFDhadronCorrSystUnc*)f->Get("SystematicUncertainty");
  TH1D *hUncCorrMin=syst->GetHistoTotFlatMin();
  TH1D *hUncCorrMax=syst->GetHistoTotFlatMax();
  TH1D *hFDsub=(TH1D*)f->Get("hDataCorrectedTempl0CentrFprompt");
  TGraphAsymmErrors *gr=syst->GetTotNonFlatUncGraph();
  gStyle->SetOptStat(0000);
  
  TCanvas *cDraw=new TCanvas("cDraw","cDraw",800,800);
  cDraw->cd();
  cDraw->SetLeftMargin(0.15);
  cDraw->SetRightMargin(0.05);
  cDraw->SetTicks();

  hFDsub->SetLineColor(kBlack);
  hFDsub->SetMarkerColor(kBlack);
  hFDsub->SetXTitle("#Delta#varphi (rad)");
  hFDsub->SetYTitle(Form("#frac{1}{#it{N}_{%s}}#frac{d#it{N}^{assoc}}{d#Delta#varphi} (rad^{-1})",strMeson.Data()));
  hFDsub->GetYaxis()->SetTitleOffset(1.3);
  hFDsub->GetYaxis()->SetRangeUser(0,max);
  hFDsub->SetTitle("");
  gr->SetLineColor(kBlack);
  gr->SetMarkerColor(kBlack);
  gr->SetFillStyle(0);
  hFDsub->Draw();
  gr->Draw("E2");

  TLatex *tlTitle=new TLatex(0.18,0.85,Form("#bf{%s meson-charged particle azimuthal correlations}",strMeson.Data()));
  tlTitle->SetNDC();
  tlTitle->Draw();
  tlTitle->SetTextSize(0.03);

  TLatex *tSystem=new TLatex(0.18,0.80,"#bf{pp, #sqrt{#it{s}} = 7 TeV, L_{int} = 5 nb^{-1}}");
  if(system==1) tSystem->SetTitle("#bf{p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, L_{int} = 50 #mub^{-1}}");
  tSystem->SetNDC();
  tSystem->SetTextSize(0.03);
  tSystem->Draw();
  
  TLatex *tptD=new TLatex(0.18,0.75,Form("#bf{%.0f < #it{p}_{T}^{%s} < %.0f GeV/c,  #it{p}_{T}^{assoc} = %.1f - %.1f  GeV/c}",ptminD,ptmaxD,strMeson.Data(),ptminAss,pTmaxAss));
  tptD->SetNDC();
  tptD->SetTextSize(0.03);
  tptD->Draw();

  TLatex *tDEta=new TLatex(0.18,0.70,Form("#bf{|#Delta#eta| < %.1f}", deltaeta));
  if(rangeLabel) tDEta->SetTitle(Form("#bf{|y_{%s}|<0.5, |#eta_{assoc}|<0.8, |#Delta#eta|<%.1f, corrected by Event Mixing}",strMeson.Data(),deltaeta));
  tDEta->SetNDC();
  tDEta->SetTextSize(0.03);
  tDEta->Draw();

  //TLatex *tlTitle=new TLatex(0.60,0.64,"ALICE Preliminary");
  //tlTitle->SetNDC();
  //tlTitle->Draw();
  //tlTitle->SetTextFont(42);
  //tlTitle->SetTextSize(0.04);

  TLatex *tUncertainty;
  if(TMath::Abs(hUncCorrMin->GetBinContent(1)-hUncCorrMax->GetBinContent(1))<0.001)tUncertainty=new TLatex(0.2,0.6,Form("#bf{%.0f#% scale uncertainty}",hUncCorrMin->GetBinContent(1)*100.));
  else tUncertainty=new TLatex(0.18,0.62,Form("#bf{{}^{+%.0f%s}_{-%.0f%s} scale uncertainty}","%","%",TMath::Abs(hUncCorrMax->GetBinContent(1))*100.,TMath::Abs(hUncCorrMin->GetBinContent(1)*100.)));
  tUncertainty->SetNDC();
  tUncertainty->SetTextSize(0.04);
  tUncertainty->Draw();

  
  TCanvas *cVaryHisto=new TCanvas("cVaryHisto","cVaryHisto",700,700);
  cVaryHisto->cd();
  hFDsub->DrawCopy();
  TH1D *hVaryUp=syst->GetVariedHisto(hFDsub,gr,1);
    if(strMeson == "D^{0}"){
        
        hVaryUp->SetMarkerColor(kRed);
        hVaryUp->SetMarkerStyle(20);
        
    }else if(strMeson == "D^{+}"){
        
        hVaryUp->SetMarkerColor(kGreen+3);
        hVaryUp->SetMarkerStyle(21);
        
    }else if(strMeson == "D^{*}"){
        
        cout << "WAs D+ " << endl;
        hVaryUp->SetMarkerColor(kAzure-2);
        hVaryUp->SetMarkerStyle(22);
        
    }
    
  TH1D *hVaryDown=syst->GetVariedHisto(hFDsub,gr,0);
    if(strMeson == "D^{0}"){
        
        hVaryDown->SetMarkerColor(kRed);
        hVaryDown->SetMarkerStyle(20);
        
    }else if(strMeson == "D^{+}"){
        
        hVaryDown->SetMarkerColor(kGreen+3);
        hVaryDown->SetMarkerStyle(21);
        
    }else if(strMeson == "D^{*}"){
        
        hVaryDown->SetMarkerColor(kAzure-2);
        hVaryDown->SetMarkerStyle(22);
        
    }
    
  hVaryUp->Draw("same");
  hVaryDown->Draw("same");

  TString strfileout="CanvaAndVariedHisto";
  if(strMeson.Contains("D^{0}"))strfileout.Append(Form("DzeroPt%.0fto%.0fassocPt%.1f.root",ptminD,ptmaxD,ptminAss));
  else if(strMeson.Contains("D^{*+}"))strfileout.Append(Form("DstarPt%.0fto%.0fassocPt%.1f.root",ptminD,ptmaxD,ptminAss));
  else if(strMeson.Contains("D^{+}"))strfileout.Append(Form("DplusPt%.0fto%.0fassocPt%.1fto.root",ptminD,ptmaxD,ptminAss));

  cDraw->Update();
  TFile *fout=new TFile(strfileout.Data(),"RECREATE");
  strfileout.ReplaceAll(".root",".png");
  fout->cd();
  cDraw->Write();
  cDraw->Print(strfileout.Data());
  hVaryUp->Write();
  hVaryDown->Write();
  fout->Close();
}


//______________________________________________________________________________________________
void OpenOutputFileAndDrawReflect(TString strfile,Double_t ptminD,Double_t ptmaxD,Int_t meson=AliHFCorrelationUtils::kDzero,Double_t ptminAss=0.3,Double_t ptmaxAss=-99,Double_t deltaeta=1, Int_t system=0,Int_t reflect=0,Double_t max=10,Int_t rebin=-1, Bool_t rangeLabel=kFALSE) {
  gStyle->SetOptStat(0000);
  TString strMeson="D^{0}";
  if(meson==AliHFCorrelationUtils::kDzero)strMeson="D^{0}";
  else if(meson==AliHFCorrelationUtils::kDstar)strMeson="D^{*+}";
  else if(meson==AliHFCorrelationUtils::kDplus)strMeson="D^{+}";
  else {
    Printf("Wrong meson selected");
    return;
  }
  TFile *f=TFile::Open(strfile.Data(),"READ");
  AliHFDhadronCorrSystUnc *so=(AliHFDhadronCorrSystUnc*)f->Get("SystematicUncertainty");
  so->SetName("InputSystematicUncertainty");


  TH1D *h=(TH1D*)f->Get("hDataCorrectedTempl0CentrFprompt");

  TH1D *hFDm=so->GetHistoFDmin();
  TH1D *hFDM=so->GetHistoFDmax();
  TH1D *hFDsub,*hFDm2,*hFDM2;
  if(reflect==1){
    hFDsub=AliHFCorrelationUtils::ReflectHisto(h,0.5);
    hFDm2=AliHFCorrelationUtils::ReflectHisto(hFDm,0.5);
    hFDM2=AliHFCorrelationUtils::ReflectHisto(hFDM,0.5);
  }
  else{
    hFDsub=h;
    hFDM2=hFDM;
    hFDm2=hFDm;
  }
  AliHFDhadronCorrSystUnc *syst=new AliHFDhadronCorrSystUnc();
  syst->SetName("SystematicUncertainty");
  if(rebin>0){
    hFDsub->Rebin(rebin);
    hFDsub->Scale(1./rebin);

    hFDm2->Rebin(rebin);
    hFDm2->Scale(1./rebin);
    hFDM2->Rebin(rebin);
    hFDM2->Scale(1./rebin);

  }
  syst->SetHistoTemplate(hFDsub);
  if(system==2)syst->InitStandardUncertaintiesPPb2016(meson,(ptmaxD+ptminD)/2.,ptminAss,ptmaxAss);
  if(system==AliHFCorrelationUtils::kpPb)syst->InitStandardUncertaintiesPPb2013(meson,(ptmaxD+ptminD)/2.,ptminAss);
  else if(system==AliHFCorrelationUtils::kpp)syst->InitStandardUncertaintiesPP2010(meson,(ptmaxD+ptminD)/2.,ptminAss);
  syst->SetHistoBeautyFDmin(hFDm2,"",kTRUE);
  syst->SetHistoBeautyFDmax(hFDM2,"",kTRUE);
  
  syst->BuildSystUncertaintyPlotVsDeltaPhi(hFDsub,1);

  TH1D *hUncCorrMin=syst->GetHistoTotFlatMin();
  TH1D *hUncCorrMax=syst->GetHistoTotFlatMax();
  //  TH1D *hFDsub=(TH1D*)f->Get("hDataCorrectedTempl0CentrFprompt");
  TGraphAsymmErrors *gr=syst->GetTotNonFlatUncGraph();
  gStyle->SetOptStat(0000);
  
  TCanvas *cDraw=new TCanvas("cDraw","cDraw",800,800);
  cDraw->cd();
  cDraw->SetLeftMargin(0.15);
  cDraw->SetRightMargin(0.05);
  cDraw->SetTicks();

  hFDsub->SetLineColor(kBlack);
  hFDsub->SetMarkerColor(kBlack);
  hFDsub->SetXTitle("#Delta#varphi (rad)");
  hFDsub->SetYTitle(Form("#frac{1}{#it{N}_{%s}}#frac{d#it{N}^{assoc}}{d#Delta#varphi} (rad^{-1})",strMeson.Data()));
  hFDsub->GetYaxis()->SetTitleOffset(1.3);
  hFDsub->GetYaxis()->SetRangeUser(0,max);
  hFDsub->SetTitle("");
  gr->SetLineColor(kBlack);
  gr->SetMarkerColor(kBlack);
  gr->SetFillStyle(0);
  hFDsub->SetMarkerColor(kRed);
  hFDsub->SetMarkerStyle(20);
  hFDsub->Draw();
  gr->Draw("E2");

  TLatex *tlTitle=new TLatex(0.18,0.85,Form("#bf{%s meson-charged particle azimuthal correlations}",strMeson.Data()));
  tlTitle->SetNDC();
  tlTitle->Draw();
  tlTitle->SetTextSize(0.03);

  TLatex *tSystem=new TLatex(0.18,0.80,"#bf{pp, #sqrt{#it{s}} = 7 TeV, L_{int} = 5 nb^{-1}}");
  if(system==1 || system==2) tSystem->SetTitle("#bf{p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, L_{int} = 50 #mub^{-1}}");
  tSystem->SetNDC();
  tSystem->SetTextSize(0.03);
  tSystem->Draw();
  
  TLatex *tptD;
  if(ptmaxAss<0)tptD=new TLatex(0.18,0.75,Form("#bf{%.0f < #it{p}_{T}^{%s} < %.0f GeV/c,  #it{p}_{T}^{assoc} > %.1f GeV/c}",ptminD,ptmaxD,strMeson.Data(),ptminAss));
  else tptD=new TLatex(0.18,0.75,Form("#bf{%.0f < #it{p}_{T}^{%s} < %.0f GeV/c, %.1f< #it{p}_{T}^{assoc} < %.1f GeV/c}",ptminD,ptmaxD,strMeson.Data(),ptminAss,ptmaxAss));
  tptD->SetNDC();
  tptD->SetTextSize(0.03);
  tptD->Draw();

  TLatex *tDEta=new TLatex(0.18,0.70,Form("#bf{|#Delta#eta|<%.1f}",deltaeta));
  if(rangeLabel) tDEta->SetTitle(Form("#bf{|y_{%s}|<0.5, |#eta_{assoc}|<0.8, |#Delta#eta|<%.1f, corrected by Event Mixing}",strMeson.Data(),deltaeta));
  tDEta->SetNDC();
  tDEta->SetTextSize(0.03);
  tDEta->Draw();

  TLatex *tlTitle=new TLatex(0.60,0.64,"ALICE Preliminary");
  tlTitle->SetNDC();
  tlTitle->Draw();
  tlTitle->SetTextFont(42);
  tlTitle->SetTextSize(0.04);

  TLatex *tUncertainty;
  if(TMath::Abs(hUncCorrMin->GetBinContent(1)-hUncCorrMax->GetBinContent(1))<0.001)tUncertainty=new TLatex(0.2,0.6,Form("#bf{%.0f#% scale uncertainty}",hUncCorrMin->GetBinContent(1)*100.));
  else tUncertainty=new TLatex(0.18,0.62,Form("#bf{{}^{+%.0f%s}_{-%.0f%s} scale uncertainty}","%","%",TMath::Abs(hUncCorrMax->GetBinContent(1))*100.,TMath::Abs(hUncCorrMin->GetBinContent(1)*100.)));
  tUncertainty->SetNDC();
  tUncertainty->SetTextSize(0.04);
  tUncertainty->Draw();

  
  TCanvas *cVaryHisto=new TCanvas("cVaryHisto","cVaryHisto",700,700);
  cVaryHisto->cd();
  hFDsub->DrawCopy();
  TH1D *hVaryUp=syst->GetVariedHisto(hFDsub,gr,1);
    if(strMeson == "D^{0}"){
        
        hVaryUp->SetMarkerColor(kRed);
        hVaryUp->SetMarkerStyle(20);
        
    }else if(strMeson == "D^{+}"){
        
        hVaryUp->SetMarkerColor(kGreen+3);
        hVaryUp->SetMarkerStyle(21);
        
    }else if(strMeson == "D^{*}"){
        
        cout << "WAs D+ " << endl;
        hVaryUp->SetMarkerColor(kAzure-2);
        hVaryUp->SetMarkerStyle(22);
        
    }
    
  TH1D *hVaryDown=syst->GetVariedHisto(hFDsub,gr,0);
    if(strMeson == "D^{0}"){
        
        hVaryDown->SetMarkerColor(kRed);
        hVaryDown->SetMarkerStyle(20);
        
    }else if(strMeson == "D^{+}"){
        
        hVaryDown->SetMarkerColor(kGreen+3);
        hVaryDown->SetMarkerStyle(21);
        
    }else if(strMeson == "D^{*}"){
        
        hVaryDown->SetMarkerColor(kAzure-2);
        hVaryDown->SetMarkerStyle(22);
        
    }
    hVaryUp->Draw("same");
  hVaryDown->Draw("same");

  TString strfileout="CanvaAndVariedHisto";
  if(system==AliHFCorrelationUtils::kpp)strfileout.Append("pp");
  else if(system==AliHFCorrelationUtils::kpPb || system==2)strfileout.Append("pPb");

  if(ptmaxAss<0){
    if(strMeson.Contains("D^{0}"))strfileout.Append(Form("DzeroPt%.0fto%.0fassocPt%.1fto99.0.root",ptminD,ptmaxD,ptminAss));
    else if(strMeson.Contains("D^{*+}"))strfileout.Append(Form("DstarPt%.0fto%.0fassocPt%.1fto99.0.root",ptminD,ptmaxD,ptminAss));
    else if(strMeson.Contains("D^{+}"))strfileout.Append(Form("DplusPt%.0fto%.0fassocPt%.1fto99.0.root",ptminD,ptmaxD,ptminAss));
  }
  else{
    if(strMeson.Contains("D^{0}"))strfileout.Append(Form("DzeroPt%.0fto%.0fassocPt%.1fto%.1f.root",ptminD,ptmaxD,ptminAss,ptmaxAss));
    else if(strMeson.Contains("D^{*+}"))strfileout.Append(Form("DstarPt%.0fto%.0fassocPt%.1fto%.1f.root",ptminD,ptmaxD,ptminAss,ptmaxAss));
    else if(strMeson.Contains("D^{+}"))strfileout.Append(Form("DplusPt%.0fto%.0fassocPt%.1fto%.1f.root",ptminD,ptmaxD,ptminAss,ptmaxAss));
  }

  cDraw->Update();
  TFile *fout=new TFile(strfileout.Data(),"RECREATE");
  strfileout.ReplaceAll(".root",".png");
  fout->cd();
  cDraw->Write();  
  cDraw->SaveAs(strfileout.Data());
  strfileout.ReplaceAll(".png",".eps");
  cDraw->SaveAs(strfileout.Data());
  hVaryUp->Write();
  hVaryDown->Write();
  hFDsub->Write();
  syst->Write();
  fout->Close();
}

//______________________________________________________________________________________________
void SubtractMCclosureModulation(TH1D *h, Double_t ptD, Double_t ptTrmin, Double_t ptTrmax) {

   Int_t bin1L = h->GetXaxis()->FindBin(0.1);
   Int_t bin1R = h->GetXaxis()->FindBin(-0.1);
   Int_t bin2L = h->GetXaxis()->FindBin(0.3);
   Int_t bin2R = h->GetXaxis()->FindBin(-0.3);
   Int_t bin3L = h->GetXaxis()->FindBin(0.5);
   Int_t bin3R = h->GetXaxis()->FindBin(-0.5);
   Int_t bin4L = h->GetXaxis()->FindBin(0.7);
   Int_t bin4R = h->GetXaxis()->FindBin(-0.7);
   Int_t bin5L = h->GetXaxis()->FindBin(0.9);
   Int_t bin5R = h->GetXaxis()->FindBin(-0.9);

   Double_t mod[5] = {0.,0.,0.,0.,0.};
   AliHFCorrelationUtils::GetMCClosureModulation(ptD,ptTrmin,ptTrmax,mod);
   if(mod[0]==-999) printf("Error in retrieving MC closure modulation! Results will be not corrected by it\n");

   printf("Applying modulations of %1.4f, %1.4f, %1.4f, %1.4f, %1.4f\n",mod[0],mod[1],mod[2],mod[3],mod[4]);
   printf("Bin1L old data %1.5f\n",h->GetBinContent(bin1L));
   printf("Bin1R old data %1.5f\n",h->GetBinContent(bin1R));
   printf("Bin2L old data %1.5f\n",h->GetBinContent(bin2L));
   printf("Bin2R old data %1.5f\n",h->GetBinContent(bin2R));
   printf("Bin3L old data %1.5f\n",h->GetBinContent(bin3L));
   printf("Bin3R old data %1.5f\n",h->GetBinContent(bin3R));
   printf("Bin4L old data %1.5f\n",h->GetBinContent(bin4L));
   printf("Bin4R old data %1.5f\n",h->GetBinContent(bin4R));
   printf("Bin5L old data %1.5f\n",h->GetBinContent(bin5L));
   printf("Bin5R old data %1.5f\n",h->GetBinContent(bin5R));
   
   //Note that the modulation is applyied, bin-by-bin, for the five three bins next to deltaPhi=0, before template correction, and after purity correction!
   //The modulation is SUBTRACTED from data (so, if negative, as in the shoulders, it's actually added to data)

   h->SetBinContent(bin1L,h->GetBinContent(bin1L)*(mod[0]));
   h->SetBinContent(bin1R,h->GetBinContent(bin1R)*(mod[0]));
   h->SetBinContent(bin2L,h->GetBinContent(bin2L)*(mod[1]));
   h->SetBinContent(bin2R,h->GetBinContent(bin2R)*(mod[1]));
   h->SetBinContent(bin3L,h->GetBinContent(bin3L)*(mod[2]));
   h->SetBinContent(bin3R,h->GetBinContent(bin3R)*(mod[2]));
   h->SetBinContent(bin4L,h->GetBinContent(bin4L)*(mod[3]));
   h->SetBinContent(bin4R,h->GetBinContent(bin4R)*(mod[3]));
   h->SetBinContent(bin5L,h->GetBinContent(bin5L)*(mod[4]));
   h->SetBinContent(bin5R,h->GetBinContent(bin5R)*(mod[4]));	
   h->SetBinError(bin1L,h->GetBinError(bin1L)*(mod[0]));
   h->SetBinError(bin1R,h->GetBinError(bin1R)*(mod[0]));
   h->SetBinError(bin2L,h->GetBinError(bin2L)*(mod[1]));
   h->SetBinError(bin2R,h->GetBinError(bin2R)*(mod[1]));
   h->SetBinError(bin3L,h->GetBinError(bin3L)*(mod[2]));
   h->SetBinError(bin3R,h->GetBinError(bin3R)*(mod[2]));
   h->SetBinError(bin4L,h->GetBinError(bin4L)*(mod[3]));
   h->SetBinError(bin4R,h->GetBinError(bin4R)*(mod[3]));
   h->SetBinError(bin5L,h->GetBinError(bin5L)*(mod[4]));
   h->SetBinError(bin5R,h->GetBinError(bin5R)*(mod[4]));

   printf("Bin1L new data %1.5f\n",h->GetBinContent(bin1L));
   printf("Bin1R new data %1.5f\n",h->GetBinContent(bin1R));
   printf("Bin2L new data %1.5f\n",h->GetBinContent(bin2L));
   printf("Bin2R new data %1.5f\n",h->GetBinContent(bin2R));
   printf("Bin3L new data %1.5f\n",h->GetBinContent(bin3L));
   printf("Bin3R new data %1.5f\n",h->GetBinContent(bin3R));
   printf("Bin4L new data %1.5f\n",h->GetBinContent(bin4L));
   printf("Bin4R new data %1.5f\n",h->GetBinContent(bin4R));
   printf("Bin5L new data %1.5f\n",h->GetBinContent(bin5L));
   printf("Bin5R new data %1.5f\n",h->GetBinContent(bin5R));

   return;
}
