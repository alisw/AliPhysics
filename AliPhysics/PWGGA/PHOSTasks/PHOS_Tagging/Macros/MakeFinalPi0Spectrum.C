void MakeFinalPi0Spectrum(const Int_t cen=6)
{
  TFile fout("PHOS_pi0_pPb_MB.root","update") ;
  TFile *fin  = new TFile("raw_MBmix.root"); 
  TFile *fin2  = new TFile("raw_MBmix2.root"); 
 
  const Int_t nPID=1 ;
  char cPID[12][12] ;
/*  
  sprintf(cPID[0],"Emin3_All") ;
  sprintf(cPID[1],"Emin3_Disp");
  sprintf(cPID[2],"Emin3_CPV") ;
  sprintf(cPID[3],"Emin3_Both"); 
  */
  sprintf(cPID[0],"All") ;
  sprintf(cPID[1],"Disp");
  sprintf(cPID[2],"CPV") ;
  sprintf(cPID[3],"Both"); 
  

  TH1D * nr1[nPID] ;
  TH1D * nr1int[nPID] ;
  TH1D * nr2[nPID] ;
  TH1D * nr2int[nPID] ;
  TH1D * nr1B[nPID] ;
  TH1D * nr1intB[nPID] ;
  TH1D * nr2B[nPID] ;
  TH1D * nr2intB[nPID] ;

  char key[55];


  TFile * feff = TFile::Open("Efficiency/PHOS_pi0_eff.root") ;
  TF1 * eff[nPID] ;
  TH1D * heff[nPID] ;
  TF1 * effInt[nPID] ;
  TH1D * heffInt[nPID] ;
  TH1D * ha[nPID] ; //efficiency with errors
  for(Int_t iPID=0; iPID<nPID; iPID++){
    sprintf(key,"yeild1_GS_Emin3_%s_cen%d",cPID[iPID],cen) ;
    nr1[iPID]=(TH1D*)fin->Get(key) ; 
    nr1B[iPID]=(TH1D*)fin2->Get(key) ; 
    sprintf(key,"yeild1_int_GS_Emin3_%s_cen%d",cPID[iPID],cen) ;
    nr1int[iPID]=(TH1D*)fin->Get(key) ; 
    nr1intB[iPID]=(TH1D*)fin2->Get(key) ; 
    sprintf(key,"yeild2_CB_Emin3_%s_cen%d",cPID[iPID],cen) ;
    nr2[iPID]=(TH1D*)fin->Get(key) ; 
    nr2B[iPID]=(TH1D*)fin2->Get(key) ; 
    sprintf(key,"yeild2_int_CB_Emin3_%s_cen%d",cPID[iPID],cen) ;
    nr2int[iPID]=(TH1D*)fin->Get(key) ; 
    nr2intB[iPID]=(TH1D*)fin2->Get(key) ; 
  
    //Divide by bin width etc: dy pt dpt dphi
    for(Int_t i=1;i<= nr1[iPID]->GetNbinsX();i++){

      Double_t pt  =TMath::TwoPi()*nr1[iPID]->GetXaxis()->GetBinCenter(i);
      nr1[iPID]->SetBinContent(i,nr1[iPID]->GetBinContent(i)/nr1[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      nr1[iPID]->SetBinError(i,nr1[iPID]->GetBinError(i)/nr1[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      nr1int[iPID]->SetBinContent(i,nr1int[iPID]->GetBinContent(i)/nr1int[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      nr1int[iPID]->SetBinError(i,nr1int[iPID]->GetBinError(i)/nr1int[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;

      nr2[iPID]->SetBinContent(i,nr2[iPID]->GetBinContent(i)/nr2[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      nr2[iPID]->SetBinError(i,nr2[iPID]->GetBinError(i)/nr2[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      nr2int[iPID]->SetBinContent(i,nr2int[iPID]->GetBinContent(i)/nr2int[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      nr2int[iPID]->SetBinError(i,nr2int[iPID]->GetBinError(i)/nr2int[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      
      nr1B[iPID]->SetBinContent(i,nr1B[iPID]->GetBinContent(i)/nr1B[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      nr1B[iPID]->SetBinError(i,nr1B[iPID]->GetBinError(i)/nr1B[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      nr1intB[iPID]->SetBinContent(i,nr1intB[iPID]->GetBinContent(i)/nr1intB[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      nr1intB[iPID]->SetBinError(i,nr1intB[iPID]->GetBinError(i)/nr1intB[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;

      nr2B[iPID]->SetBinContent(i,nr2B[iPID]->GetBinContent(i)/nr2B[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      nr2B[iPID]->SetBinError(i,nr2B[iPID]->GetBinError(i)/nr2B[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      nr2intB[iPID]->SetBinContent(i,nr2intB[iPID]->GetBinContent(i)/nr2intB[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      nr2intB[iPID]->SetBinError(i,nr2intB[iPID]->GetBinError(i)/nr2intB[iPID]->GetXaxis()->GetBinWidth(i)/pt) ;
      
    }

    //Correct for efficiency
    eff[iPID] = (TF1*)feff->Get(Form("eff_Pi0_Gaus_pPb_Emin3_%s_cen%d",cPID[iPID],6)) ; 
    eff[iPID]->SetRange(0.,70.) ;
    heff[iPID]=(TH1D*)feff->Get(Form("yeild_GS_Emin3_%s_cen%d",cPID[iPID],6)) ;
    nr1[iPID]->Divide(heff[iPID]) ;
    nr1B[iPID]->Divide(heff[iPID]) ;
//     nr1[iPID]->Divide(eff[iPID]) ;
//     nr1B[iPID]->Divide(eff[iPID]) ;
// 
    effInt[iPID] = (TF1*)feff->Get(Form("eff_int_Pi0_Gaus_pPb_Emin3_%s_cen%d",cPID[iPID],6)) ;
    effInt[iPID]->SetRange(0.,70.) ;
    heffInt[iPID]=(TH1D*)feff->Get(Form("yeild_int_GS_Emin3_%s_cen%d",cPID[iPID],6)) ;
    nr1int[iPID]->Divide(heffInt[iPID]) ;
    nr1intB[iPID]->Divide(heffInt[iPID]) ;
//     nr1int[iPID]->Divide(effInt[iPID]) ;
//     nr1intB[iPID]->Divide(effInt[iPID]) ;
    
    heff[iPID]=(TH1D*)feff->Get(Form("yeild_CB_Emin3_%s_cen%d",cPID[iPID],6)) ;
    nr2[iPID]->Divide(heff[iPID]) ;
    nr2B[iPID]->Divide(heff[iPID]) ;
    heffInt[iPID]=(TH1D*)feff->Get(Form("yeild_int_CB_Emin3_%s_cen%d",cPID[iPID],6)) ;
    nr2int[iPID]->Divide(heffInt[iPID]) ;
    nr2intB[iPID]->Divide(heffInt[iPID]) ;

  }
  

  //First for each PID evaluate spectrum (mean of pol2,3,int2,3) 
  //and sys error related to Bg subtraction
  //Stat err is the smallest error in a given bin in nr1,nr2,nr1int,nr2int
  TH1D * hStatRY[nPID] ;
  TH1D * hSysRY[nPID] ;
  for(Int_t iPID=0; iPID<nPID; iPID++){
    hStatRY[iPID] = (TH1D*)nr1[iPID]->Clone(Form("hStatRawYield%d",iPID)) ;
    hSysRY[iPID] = (TH1D*)nr1[iPID]->Clone(Form("hSysRawYield%d",iPID)) ;
    hStatRY[iPID]->Reset() ;
    hSysRY[iPID]->Reset() ;
    for(Int_t i=1;i<= hStatRY[iPID]->GetNbinsX();i++){
      Double_t w=0.;
      Double_t mean=0. ;
      Double_t rms=0. ;
      Double_t stat=nr1[iPID]->GetBinError(i) ;
      Double_t wi=nr1[iPID]->GetBinError(i) ;
      if(wi>0.){
	if(nr1[iPID]->GetBinError(i)<stat) stat=nr1[iPID]->GetBinError(i) ;
	mean+=nr1[iPID]->GetBinContent(i)/wi/wi ;
	rms+=nr1[iPID]->GetBinContent(i)*nr1[iPID]->GetBinContent(i)/wi/wi ;
	w+=1./wi/wi ;
      }
      wi=nr2[iPID]->GetBinError(i) ;
      if(wi>0.){
	if(nr2[iPID]->GetBinError(i)<stat) stat=nr2[iPID]->GetBinError(i) ;
	mean+=nr2[iPID]->GetBinContent(i)/wi/wi ;
	rms+=nr2[iPID]->GetBinContent(i)*nr2[iPID]->GetBinContent(i)/wi/wi ;
	w+=1./wi/wi ;
      }
      wi=nr1int[iPID]->GetBinError(i) ;
      if(wi>0.){
	if(nr1int[iPID]->GetBinError(i)<stat) stat=nr1int[iPID]->GetBinError(i) ;
	mean+=nr1int[iPID]->GetBinContent(i)/wi/wi ;
	rms+=nr1int[iPID]->GetBinContent(i)*nr1int[iPID]->GetBinContent(i)/wi/wi ;
	w+=1./wi/wi ;
      }
      wi=nr2int[iPID]->GetBinError(i) ;
      if(wi>0.){
	if(nr2int[iPID]->GetBinError(i)<stat) stat=nr2int[iPID]->GetBinError(i) ;
	mean+=nr2int[iPID]->GetBinContent(i)/wi/wi ;
	rms+=nr2int[iPID]->GetBinContent(i)*nr2int[iPID]->GetBinContent(i)/wi/wi ;
	w+=1./wi/wi ;
      }
      if(w>0.){
	mean=mean/w ;
	rms=rms/w-mean*mean ;
        if(TMath::Abs(rms)<1.e-8*mean*mean)rms=0.;
        hStatRY[iPID]->SetBinContent(i,mean) ;
        hStatRY[iPID]->SetBinError(i,stat) ;
        hSysRY[iPID]->SetBinContent(i,mean) ;
        hSysRY[iPID]->SetBinError(i,TMath::Sqrt(rms)) ;
      }	
    }
  }  
  //Estimate Bg
  TH1D * hSysBg = (TH1D*)hStatRY[0]   ->Clone(Form("hPi0_PbPb_cen%d_SystBg",cen)) ;
  hSysBg->Reset() ;
  for(Int_t i=1;i<= hStatRY[0]->GetNbinsX();i++){
      Double_t w=0.;
      Double_t mean=0. ;
      Double_t wi=nr1[0]->GetBinError(i) ;
      if(wi>0.){
	mean+=(nr1[0]->GetBinContent(i)-nr1B[0]->GetBinContent(i))*(nr1[0]->GetBinContent(i)-nr1B[0]->GetBinContent(i))/wi/wi ;
	w+=1./wi/wi ;
      }
      wi=nr2[0]->GetBinError(i) ;
      if(wi>0.){
	mean+=(nr2[0]->GetBinContent(i)-nr2B[0]->GetBinContent(i))*(nr2[0]->GetBinContent(i)-nr2B[0]->GetBinContent(i))/wi/wi ;
	w+=1./wi/wi ;
      }
      wi=nr1int[0]->GetBinError(i) ;
      if(wi>0.){
	mean+=(nr1int[0]->GetBinContent(i)-nr1intB[0]->GetBinContent(i))*(nr1int[0]->GetBinContent(i)-nr1intB[0]->GetBinContent(i))/wi/wi ;
	w+=1./wi/wi ;
      }
      wi=nr2int[0]->GetBinError(i) ;
      if(wi>0.){
	mean+=(nr2int[0]->GetBinContent(i)-nr2intB[0]->GetBinContent(i))*(nr2int[0]->GetBinContent(i)-nr2intB[0]->GetBinContent(i))/wi/wi ;
	w+=1./wi/wi ;
      }
      if(w>0.){
	mean=mean/w ;
        if(hStatRY[0]->GetBinContent(i))
//         hSysBg->SetBinContent(i,TMath::Sqrt(mean)/hStatRY[0]->GetBinContent(i)) ;
        hSysBg->SetBinContent(i,hStatRY[0]->GetBinContent(i)) ;
        hSysBg->SetBinError(i,TMath::Sqrt(mean)) ;
      }	

  }  
//  hSysBg->Draw() ; return ;
  
  
  
//hStatRY[3]->Draw() ;return ;

  //Final spectrum  
  TH1D * hStat = (TH1D*)nr1[0]->Clone(Form("hPi0_PbPb_cen%d_Stat",cen)) ;
  //And full sys error
  TH1D * hSys  = (TH1D*)hStat   ->Clone(Form("hPi0_PbPb_cen%d_Syst",cen)) ;
  //Sys errors appearing in Raa
  TH1D * hSysRaa = (TH1D*)hStat   ->Clone(Form("hPi0_PbPb_cen%d_SystRaa",cen)) ;
  //Sys errors appearing in gdir
  TH1D * hSysA = (TH1D*)hStat   ->Clone(Form("hPi0_PbPb_cen%d_SystA",cen)) ;
  TH1D * hSysB = (TH1D*)hStat   ->Clone(Form("hPi0_PbPb_cen%d_SystB",cen)) ;
  TH1D * hSysC = (TH1D*)hStat   ->Clone(Form("hPi0_PbPb_cen%d_SystC",cen)) ;
  //Uncorrelated sys errors 
  TH1D * hSysUC = (TH1D*)hStat   ->Clone(Form("hPi0_PbPb_cen%d_SystUnc",cen)) ;
 
  TH1D * hSysSignal = (TH1D*)hStat   ->Clone(Form("hPi0_PbPb_cen%d_SystSignal",cen)) ;
  TH1D * hSysAccept = (TH1D*)hStat   ->Clone(Form("hPi0_PbPb_cen%d_SystAccept",cen)) ;
  TH1D * hSysMaterial = (TH1D*)hStat   ->Clone(Form("hPi0_PbPb_cen%d_SystMaterial",cen)) ;
  TH1D * hSysNonlin = (TH1D*)hStat   ->Clone(Form("hPi0_PbPb_cen%d_SystNonLin",cen)) ;
  TH1D * hSysEscale = (TH1D*)hStat   ->Clone(Form("hPi0_PbPb_cen%d_SystEcale",cen)) ;
  

  TString centext[10] = {"0-20%", "20-40%","40-60%","60-80%","80-100%","60-80%","0-100%", "0-10%", "40-80%", "0-40%"};
  hStat->SetTitle(Form("#pi^{0} prod.yield, %s, stat.error"       ,centext[cen].Data()));
  hSys ->SetTitle(Form("#pi^{0} prod.yield, %s, syst.error"       ,centext[cen].Data()));

  //Sys error related to raw yield extraction
  TH1D * hSysRawYield  = (TH1D*)hStat->Clone(Form("hSysRawYield%d_Syst",cen)) ;  
  //Sys error related to PID
    
  hSys->Reset() ;
  hStat->Reset() ;
  for(Int_t i=1; i<=hStat->GetNbinsX(); i++){
    Double_t w=0. ;
    Double_t mean=0. ;
    Double_t stat=hStatRY[0]->GetBinError(i) ;
    Double_t sysRawYield = 0. ;
    for(Int_t iPID=0; iPID<nPID; iPID++){
       Double_t wiStat = hStatRY[iPID]->GetBinError(i) ;
       Double_t wiSys = hSysRY[iPID]->GetBinError(i) ;
       if(stat>wiStat) stat=wiStat ;
       Double_t wi2 = wiStat*wiStat+wiSys*wiSys ;
       if(wi2>0){
         mean+=hStatRY[iPID]->GetBinContent(i)/wi2 ;
         w+=1./wi2 ;
       }
       sysRawYield+=wiSys ; //THis will be arithmetic average
    }
    if(w>0){
      hStat->SetBinContent(i,mean/w) ;
      hStat->SetBinError(i,stat) ;
//      
      hSysRawYield->SetBinContent(i,mean/w) ;
      hSysRawYield->SetBinError(i,sysRawYield/nPID) ;
    }
  }
/*  
  TH1D * box=new TH1D("box","aaa",200,0.,20.) ;
  box->Draw() ;
  TLegend * l = new TLegend(0.,0.7,0.9,0.9) ;
  for(Int_t iPID=0; iPID<nPID; iPID++){
    hStatRY[iPID]->Divide(hStat) ;
    hStatRY[iPID]->SetMarkerStyle(20+iPID) ;
    hStatRY[iPID]->SetMarkerColor(1+iPID) ;
    hStatRY[iPID]->Draw("same") ;
    l->AddEntry(hStatRY[iPID],cPID[iPID],"p") ;
  }
  l->Draw() ;
  return ;
  //hStat->Draw() ; return ;
*/
  //Add other sys errors   2.05921e+00,6.77850e+00
  TF1 * globalE  = new TF1("globalE","TMath::Abs(1.-((1.+x/1.005/7.14/0.163)/(1+x*1.005/7.14/0.163))^7.14)",1.,30.) ;  
     
  TF1 * permod = new TF1("permod","0.0",0.,30.) ;
  TF1 * conv = new TF1("conversion","0.035",0.,30.) ;
  TF1 * conv2 = new TF1("conversion2","0.017",0.,30.) ;
  TF1 * accept= new TF1("accept","0.01",0.,30.) ;
  TF1 * pileup= new TF1("pileup","0.01",0.,30.) ;
  TF1 * nonlin= new TF1("nl","0.015+7.38*exp(-x/0.24)",0.,30.) ;
  permod->SetLineColor(52) ;
  permod->SetLineStyle(8) ;

  
  TF1 * fEffErr = new TF1("effErr","0.02",1.,20.) ; ;
  
  
 
// Boris feed-down
  TFile * ffeed = new TFile("feed_down.root") ;
  TF1 * ftmp = (TF1*)ffeed->Get("fun1") ;
  TF1 * feeddown = new TF1("feeddown",p4,0.9,20.,5) ;
  feeddown->SetParameters(ftmp->GetParameters()) ;
//   feeddown->Draw() ; return ;
  
  TF1 * feeddownErr= new TF1("fd","0.02",0.,30.) ;
  
//Print Errors
Double_t pt1=1.5 ;
Double_t pt2=7.5;
Int_t i1 = hStat->GetXaxis()->FindBin(pt1) ;
Int_t i2 = hStat->GetXaxis()->FindBin(pt2) ;
printf("Raw Yield &  $\\pm$%f    & $\\pm$%f  \n",100.*hSysRawYield->GetBinError(i1)/hSysRawYield->GetBinContent(i1),
       100.*hSysRawYield->GetBinError(i2)/hSysRawYield->GetBinContent(i2)) ;
printf("Backgournd &  $\\pm$%f    & $\\pm$%f  \n",100.*hSysBg->GetBinError(i1)/hSysBg->GetBinContent(i1),
       100.*hSysBg->GetBinError(i2)/hSysBg->GetBinContent(i2)) ;
printf("Global E &  $\\pm$%f    & $\\pm$%f  \n",100.*globalE->Eval(pt1),
       100.*globalE->Eval(pt2));
printf("Non-linearity &  $\\pm$%f    & $\\pm$ %f \n ",100.*nonlin->Eval(pt1),
       100.*nonlin->Eval(pt2));
  
  //Add errors in sys errors
  for(Int_t i=1;i<=hSys->GetNbinsX();i++){
    Double_t pt=hStat->GetXaxis()->GetBinCenter(i) ;
    Double_t mean= hStat->GetBinContent(i) ;
    Double_t a=hSysRawYield->GetBinError(i) ;  //RawYiel
    a=a*a+hSysBg->GetBinError(i)*hSysBg->GetBinError(i) ;
    a=TMath::Sqrt(a) ;
    Double_t b=0; //hSysPID->GetBinError(i) ;  //PID
    Double_t effErr = 0.02; //Estimated comparing single pi0 and LHC13b2_efix_p4 simulations
    Double_t tot= mean*mean*(
			     nonlin->Eval(pt)      *nonlin->Eval(pt) +
			     conv->Eval(pt)        *conv->Eval(pt) +
			     accept->Eval(pt)      *accept->Eval(pt) +
			     permod->Eval(pt)      *permod->Eval(pt) +
			     pileup->Eval(pt)      *pileup->Eval(pt) +
			     globalE->Eval(pt)     *globalE->Eval(pt) +
			     feeddownErr->Eval(pt)     *feeddownErr->Eval(pt) +
			     effErr*effErr
			     ); 
    hSys ->SetBinContent(i,mean) ;
    hSys ->SetBinError(i,TMath::Sqrt(tot+a*a+b*b)) ;
    //Now for Raa
    Double_t tot2= mean*mean*(
			     pileup->Eval(pt)      *pileup->Eval(pt) +
			     accept->Eval(pt)      *accept->Eval(pt) +
			     feeddownErr->Eval(pt)     *feeddownErr->Eval(pt) +
			     effErr*effErr
			     ); 
    hSysRaa ->SetBinContent(i,mean) ;
    hSysRaa ->SetBinError(i,TMath::Sqrt(tot2+a*a+b*b)) ;
    hSysA ->SetBinContent(i,mean) ;
    hSysA ->SetBinError(i,TMath::Sqrt(a*a)) ;
    hSysB ->SetBinContent(i,mean) ;
    hSysB ->SetBinError(i,TMath::Sqrt(mean*mean*(effErr*effErr + 
              nonlin->Eval(pt)*nonlin->Eval(pt) + 
              globalE->Eval(pt)*globalE->Eval(pt)+
              feeddownErr->Eval(pt)*feeddownErr->Eval(pt)  ))) ;
    hSysC ->SetBinContent(i,mean) ;
    hSysC ->SetBinError(i,TMath::Sqrt(mean*mean*(
			     accept->Eval(pt)      *accept->Eval(pt) +
			     permod->Eval(pt)      *permod->Eval(pt) +
			     pileup->Eval(pt)      *pileup->Eval(pt) 
			     ))) ;
    hSysSignal->SetBinContent(i,mean) ;                         
    hSysSignal->SetBinError(i,a) ;
//     hSysBg->SetBinContent(i,mean) ; //already filled                         
//     hSysBg->SetBinError(i,a) ;
    hSysAccept->SetBinContent(i,mean) ; 
    hSysAccept->SetBinError(i,mean*TMath::Sqrt(effErr*effErr + accept->Eval(pt)*accept->Eval(pt))) ;
    hSysMaterial->SetBinContent(i,mean) ;
    hSysMaterial->SetBinError(i,mean*conv->Eval(pt)) ;
    hSysNonlin->SetBinContent(i,mean) ;
    hSysNonlin->SetBinError(i,mean*nonlin->Eval(pt)) ;
    hSysEscale->SetBinContent(i,mean) ;
    hSysEscale->SetBinError(i,mean*globalE->Eval(pt)) ;                             
  }
  
  
  
  
printf("Total &  $\\pm$%f    & $\\pm$%f \n",100.*hSys->GetBinError(i1)/hSys->GetBinContent(i1),
       100.*hSys->GetBinError(i2)/hSys->GetBinContent(i2)) ;
  
  
  
  //SetRange
  Double_t ptmin=1.;
  Double_t ptmax=20.99;
  for(Int_t i=1;i<=hSys->GetNbinsX();i++){
    Double_t pt=hSys->GetXaxis()->GetBinCenter(i) ;
    if(pt<ptmin || pt>ptmax){   
      hStat ->SetBinContent(i,0.) ;
      hStat ->SetBinError(i,0.) ;
      hSys ->SetBinContent(i,0.) ;
      hSys ->SetBinError(i,0.) ;
      hSysRaa->SetBinContent(i,0.) ;
      hSysRaa->SetBinError(i,0.) ;
      hSysA->SetBinContent(i,0.) ;
      hSysA->SetBinError(i,0.) ;
      hSysB->SetBinContent(i,0.) ;
      hSysB->SetBinError(i,0.) ;
      hSysC->SetBinContent(i,0.) ;
      hSysC->SetBinError(i,0.) ;
      
      hSysSignal->SetBinContent(i,0.) ;                         
      hSysSignal->SetBinError(i,0.) ;
      hSysBg->SetBinContent(i,0.) ; //already filled                         
      hSysBg->SetBinError(i,0.) ;
      hSysAccept->SetBinContent(i,0.) ; 
      hSysAccept->SetBinError(i,0.) ;
      hSysMaterial->SetBinContent(i,0.) ;
      hSysMaterial->SetBinError(i,0.) ;
      hSysNonlin->SetBinContent(i,0.) ;
      hSysNonlin->SetBinError(i,0.) ;
      hSysEscale->SetBinContent(i,0.) ;
      hSysEscale->SetBinError(i,0.) ;                             
      
    }
  }  

  
  //Draw sys errors
  TH1D * hRelSys = (TH1D*)hSys->Clone("RelSys") ;
  for(Int_t i=1;i<=hSys->GetNbinsX();i++){
    Double_t mean= hSys->GetBinContent(i) ;
    Double_t a=hSys->GetBinError(i) ;
    Double_t bs=hSysRawYield->GetBinError(i) ;
    if(mean>0){
      hRelSys->SetBinContent(i,a/mean) ;
      hRelSys->SetBinError(i,0.) ;
      hSysRawYield->SetBinContent(i,bs/mean) ;
      hSysRawYield->SetBinError(i,0.) ;
    }
    else{
      hRelSys->SetBinContent(i,0.) ;
      hSysRawYield->SetBinContent(i,0.) ;
      hSysRawYield->SetBinError(i,0.) ;
    }
  }
//  Smooth(hSysPID) ;

  
  TCanvas * c = new TCanvas("SysErrors","Systematic errors") ;
  c->cd() ;
  gPad->SetLogx() ;
  hRelSys->SetMinimum(0.) ;
  hRelSys->SetMaximum(0.3) ;
  hRelSys->SetLineColor(2) ;
  hRelSys->SetLineWidth(2) ;
  hSysRawYield->SetLineWidth(2) ;
  hSysRawYield->SetLineColor(kBlue) ;
  feeddownErr->SetLineColor(kOrange+9) ;
  feeddownErr->SetLineStyle(3) ;
  
  fEffErr->SetLineStyle(1) ;
  fEffErr->SetLineColor(kAzure) ;
  fEffErr->SetMarkerColor(kAzure) ;
//   fEffErr->SetMarkerStyle(33) ;
  
  globalE->SetLineColor(kGreen+2) ;
  nonlin->SetLineColor(kAzure+10) ;
  nonlin->SetLineStyle(9) ;
  conv->SetLineColor(6) ;
  accept->SetLineColor(kOrange) ;
  hRelSys->SetXTitle("p_{t} (GeV/c)") ;
  hRelSys->SetYTitle("Rel.error") ;
  hRelSys->GetYaxis()->SetTitleOffset(1.2) ;
  hRelSys->GetXaxis()->SetRangeUser(0.8,20.) ;
  hRelSys->SetMinimum(0.006) ;
  hRelSys->SetMaximum(1.) ;
  hRelSys->Draw("h") ;
  hSysRawYield->Draw("hsame") ;
  globalE->Draw("same") ; 
  nonlin ->Draw("same") ;
  conv   ->Draw("same") ;
  permod   ->Draw("same") ;
  accept ->Draw("same") ;
   fEffErr ->Draw("same") ;
  feeddownErr->Draw("same") ;
  TLegend * l = new TLegend(0.3,0.5,0.6,0.9) ;
  l->AddEntry(hSysRawYield,"Raw extraction","l") ;
  l->AddEntry(conv,"Conversion","l") ;
  l->AddEntry(nonlin,"Non-linearity","l") ;
  l->AddEntry(accept,"Acceptance","l");
  l->AddEntry(fEffErr,"Efficiency","l");
  l->AddEntry(feeddownErr,"Feed-down correction","l");
  l->AddEntry(globalE,"Global E scale","l") ;
  l->AddEntry(hRelSys,"Total error","l") ;
  l->Draw() ;

  
  TCanvas * csp = new TCanvas("FinalSpectrum","Final Spectrum") ;
  csp->SetLogy();
  
  hStat->SetXTitle("p_{T}, GeV/c");
  hStat->SetYTitle("1/N_{ev}d^{3}N/dyd^{2}p_{t} (GeV^{-2}c^{2})");
  hSys ->SetXTitle("p_{T}, GeV/c");
  hSys ->SetYTitle("1/N_{ev}d^{3}N/dyd^{2}p_{t} (GeV^{-2}c^{2})");

  hSys ->SetFillColor(kBlue-10) ;
  hSys ->SetFillStyle(0);
  hSys ->Draw("E2") ;
  hStat->SetMarkerStyle(20) ;
  hStat->SetMarkerColor(4) ;
  hStat->SetLineColor(4) ;
  hStat->Draw("same") ;

 
 // Feed-down correction
  for(Int_t i=1; i<=hStat->GetNbinsX();i++){
    Double_t x = hStat->GetXaxis()->GetBinCenter(i) ;
    Double_t correction = 1.-feeddown->Eval(x) ;
    hStat->SetBinContent(i,correction*hStat->GetBinContent(i)) ;
    hStat->SetBinError(i,correction*hStat->GetBinError(i)) ;
    hSys->SetBinContent(i,correction*hSys->GetBinContent(i)) ;
    hSys->SetBinError(i,correction*hSys->GetBinError(i)) ;
    hSysRaa->SetBinContent(i,correction*hSys->GetBinContent(i)) ;
    hSysRaa->SetBinError(i,correction*hSys->GetBinError(i)) ;
    hSysA->SetBinContent(i,correction*hSysA->GetBinContent(i)) ;
    hSysA->SetBinError(i,correction*hSysA->GetBinError(i)) ;
    hSysB->SetBinContent(i,correction*hSysB->GetBinContent(i)) ;
    hSysB->SetBinError(i,correction*hSysB->GetBinError(i)) ;
    hSysC->SetBinContent(i,correction*hSysC->GetBinContent(i)) ;
    hSysC->SetBinError(i,correction*hSysC->GetBinError(i)) ;
    
      hSysSignal->SetBinContent(i,correction*hSysSignal->GetBinContent(i)) ;                         
      hSysSignal->SetBinError(i,correction*hSysSignal->GetBinError(i)) ;
      hSysBg->SetBinContent(i,correction*hSysBg->GetBinContent(i)) ; //already filled                         
      hSysBg->SetBinError(i,correction*hSysBg->GetBinError(i)) ;
      hSysAccept->SetBinContent(i,correction*hSysAccept->GetBinContent(i)) ; 
      hSysAccept->SetBinError(i,correction*hSysAccept->GetBinError(i)) ;
      hSysMaterial->SetBinContent(i,correction*hSysMaterial->GetBinContent(i)) ;
      hSysMaterial->SetBinError(i,correction*hSysMaterial->GetBinError(i)) ;
      hSysNonlin->SetBinContent(i,correction*hSysNonlin->GetBinContent(i)) ;
      hSysNonlin->SetBinError(i,correction*hSysNonlin->GetBinError(i)) ;
      hSysEscale->SetBinContent(i,correction*hSysEscale->GetBinContent(i)) ;
      hSysEscale->SetBinError(i,correction*hSysEscale->GetBinError(i)) ;                             
    
 }  


  TH1D * hStatNoBW = (TH1D*)hStat->Clone(Form("hPi0_PbPb_cen%d_NoBW_Stat",cen)) ;
  TH1D * hSysNoBW = (TH1D*)hSys->Clone(Form("hPi0_PbPb_cen%d_NoBW_Syst",cen)) ;
  TH1D * hSysUCNoBW = (TH1D*)hSys->Clone(Form("hPi0_PbPb_cen%d_NoBW_SystUnc",cen)) ;
  
  //Apply bin width correction  

//   TH1D * hBWcorr = BinWidthCorrection(hStat) ;
//   hSys ->Divide(hBWcorr) ;
//   hSysRaa->Divide(hBWcorr) ;
//   hStat->Divide(hBWcorr) ;
// //  hSysG->Divide(hBWcorr) ;
//   hSysUC->Divide(hBWcorr) ;



//  hStat->Draw() ;
  
fout.cd() ;  
  hSys ->Write(0,TObject::kOverwrite) ;
  hSysRaa->Write(0,TObject::kOverwrite) ;
  hSysA->Write(0,TObject::kOverwrite) ;
  hSysB->Write(0,TObject::kOverwrite) ;
  hSysC->Write(0,TObject::kOverwrite) ;
  hStat->Write(0,TObject::kOverwrite) ;
  
  hSysSignal->Write(0,TObject::kOverwrite) ;            
  hSysBg->Write(0,TObject::kOverwrite) ;                     
  hSysAccept->Write(0,TObject::kOverwrite) ;
  hSysMaterial->Write(0,TObject::kOverwrite) ;
  hSysNonlin->Write(0,TObject::kOverwrite) ;
  hSysEscale->Write(0,TObject::kOverwrite) ;
  
  fout.Close() ;
  
}

//-----------------------------------------------------------------------------
TH1D * BinWidthCorrection(TH1D * h){
  //We apply bin width a-la PHENIX 
  //Use Tsalis fit to calculate shift in y direction
 
 TF1 * fit = new TF1("Tsalis","[0]*(1.+(sqrt(x*x+0.135*0.135)-0.135)/([2]*[1]))^-[2]*([3]+[4]*x+[5]*x*x+[6]*x*x*x)",0.5,25.) ;
  fit->SetParameters(100.,0.2,11.,4.41717,0.833991,-0.418400,0.0566629) ;

 //  TF1 * fit = new TF1("hag","[0]*(([2]*[1]+sqrt(x*x+0.135*0.135))/([2]*[1]+0.135))^-[2]",0.5,25.) ;
//  TF1 * fit = new TF1("hag","[0]*([1]+sqrt(x*x+0.135*0.135))^-[2]+[3]*exp(-x/[4])",0.5,25.) ;
//  fit->SetParameters(10.,0.2,8.,1.,0.5) ;
  TCanvas * corr = new TCanvas("BWcorr","Bin width correction") ;
  Int_t col[6]={kRed,kOrange,kMagenta,kGreen,kCyan,kBlue} ;
  TH1D * hcorr[20] ;
  char key[55] ;
  Double_t rMax=10 ;
  Int_t iteration=0 ;
  TH1D * htmp = (TH1D*)h->Clone("tmp") ;
  
  while(iteration<6){
    printf(" Iteration %d: rMax=%f \n",iteration, rMax) ;
    htmp->Fit(fit,"N","",1.,20.) ;
    sprintf(key,"Ineration%d",iteration) ;
    hcorr[iteration]=(TH1D*)h->Clone(key);
    rMax= 0; 
    for(Int_t i=1;i<=h->GetNbinsX();i++){
      Double_t a=h->GetXaxis()->GetBinLowEdge(i) ;
      Double_t b=h->GetXaxis()->GetBinUpEdge(i) ;
      Double_t r=fit->Integral(a,b)/(b-a)/fit->Eval(0.5*(a+b)) ;
      hcorr[iteration]->SetBinContent(i,r) ;
      hcorr[iteration]->SetBinError(i,0.) ;
      if(rMax<r)rMax=r ;
    }
    delete htmp ;
    htmp = (TH1D*)h->Clone("tmp") ;
    htmp->Divide(hcorr[iteration]) ;
    corr->cd() ;
    hcorr[iteration]->SetLineColor(col[iteration]);
    if(iteration==0)
      hcorr[iteration]->Draw() ;
    else
      hcorr[iteration]->Draw("same") ;
    corr->Update() ;
    iteration++ ;
  } 

  return hcorr[5] ;
}
void Smooth(TH1D * h){
  
 const Int_t start=10;
 Int_t n=h->GetNbinsX();
 for(Int_t i=start+1; i<h->GetNbinsX();i++){
  Double_t av=h->GetBinContent(i-1)+h->GetBinContent(i)+h->GetBinContent(i+1) ;
  av/=3. ;
  h->SetBinContent(i,av) ;
 }
 Double_t av=h->GetBinContent(n-1)+2.*h->GetBinContent(i) ;
  av/=3. ;
  h->SetBinContent(n,av) ;
  
  
}
//Parameterization from Boris
Double_t p4(Double_t * x, Double_t * par)
{
   Double_t z=TMath::Min(7.,x[0]);
    
    Double_t f = par[0] + par[1]*z + par[2]*z*z + par[3]*z*z*z + par[4]*z*z*z*z;
    return f;
}
