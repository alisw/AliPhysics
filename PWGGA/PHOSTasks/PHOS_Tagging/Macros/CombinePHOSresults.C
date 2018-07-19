void CombinePHOSresults(){
 
  char centext[10] ;
  sprintf(centext,"00-100") ;

  TFile * fout = new TFile("PHOS_pPb_19092017.root","recreate") ;  
    
  TFile *fpi0 = new TFile("PHOS_pi0_pPb_MB.root") ;

  //Data mass, width
  TFile * fMass = new TFile("raw_MBmix.root") ;
  
  
  //MC: Efficiency
  TFile * fileeff = TFile::Open("Efficiency/PHOS_pi0_eff.root") ;
 
  //MC: Mass
  TFile * fMCMass = new TFile("Efficiency/mass_MB.root") ;
  

    //pi0s
    TH1D * hpi0Stat= (TH1D*)fpi0->Get(Form("hPi0_PbPb_cen6_Stat")) ;  
    hpi0Stat->SetName(Form("hPHOS_pi0_pPb_cen%s_Stat",centext)) ;
    hpi0Stat->SetTitle(Form("pi0 yield, centrality %s%%, stat. errors",centext)) ;
    TH1D * hpi0Syst= (TH1D*)fpi0->Get(Form("hPi0_PbPb_cen6_Syst")) ;
    hpi0Syst->SetName(Form("hPHOS_pi0_pPb_cen%s_SystTotal",centext)) ;
    hpi0Syst->SetTitle(Form("pi0 yield, centrality %s%%, total sys. errors for Double Ratio",centext)) ;
    TH1D * hpi0SystA = (TH1D*)fpi0->Get(Form("hPi0_PbPb_cen6_SystA")) ;
    hpi0SystA->SetName(Form("hPHOS_pi0_pPb_cen%s_SystA",centext)) ;
    hpi0SystA->SetTitle(Form("pi0 yield (fitted), centrality %s%%, type A sys. errors",centext)) ;
    TH1D * hpi0SystB = (TH1D*)fpi0->Get(Form("hPi0_PbPb_cen6_SystB")) ;
    hpi0SystB->SetName(Form("hPHOS_pi0_pPb_cen%s_SystB",centext)) ;
    hpi0SystB->SetTitle(Form("pi0 yield (fitted), centrality %s%%, type B sys. errors (eff)",centext)) ;
    TH1D * hpi0SystC = (TH1D*)fpi0->Get(Form("hPi0_PbPb_cen6_SystC")) ;
    hpi0SystC->SetName(Form("hPHOS_pi0_pPb_cen%s_SystC",centext)) ;
    hpi0SystC->SetTitle(Form("pi0 yield, centrality %s%%, type C sys. errors",centext)) ;
    
    TH1D * hpi0SystSig = (TH1D*)fpi0->Get(Form("hPi0_PbPb_cen6_SystSignal")) ;
    hpi0SystSig->SetName(Form("hPHOS_pi0_pPb_cen%s_SystSignalExtr",centext)) ;
    hpi0SystSig->SetTitle(Form("pi0 yield, centrality %s%%, signal extraction sys. errors",centext)) ;
    
    TH1D * hpi0SystBg = (TH1D*)fpi0->Get(Form("hPi0_PbPb_cen6_SystBg")) ;
    hpi0SystBg->SetName(Form("hPHOS_pi0_pPb_cen%s_SystBg",centext)) ;
    hpi0SystBg->SetTitle(Form("pi0 yield, centrality %s%%, bg. subtrction sys. errors",centext)) ;

    
    TH1D * hpi0SystAcc = (TH1D*)fpi0->Get(Form("hPi0_PbPb_cen6_SystAccept")) ;
    hpi0SystAcc->SetName(Form("hPHOS_pi0_pPb_cen%s_SystAccept",centext)) ;
    hpi0SystAcc->SetTitle(Form("pi0 yield, centrality %s%%, acceptance sys. errors",centext)) ;
   
    TH1D * hpi0SystMat = (TH1D*)fpi0->Get(Form("hPi0_PbPb_cen6_SystMaterial")) ;
    hpi0SystMat->SetName(Form("hPHOS_pi0_pPb_cen%s_SystMaterial",centext)) ;
    hpi0SystMat->SetTitle(Form("pi0 yield, centrality %s%%, material sys. errors",centext)) ;

    TH1D * hpi0SystNL = (TH1D*)fpi0->Get(Form("hPi0_PbPb_cen6_SystNonLin")) ;
    hpi0SystNL->SetName(Form("hPHOS_pi0_pPb_cen%s_SystNonLin",centext)) ;
    hpi0SystNL->SetTitle(Form("pi0 yield, centrality %s%%, nonlinearity sys. errors",centext)) ;

   TH1D * hpi0SystE = (TH1D*)fpi0->Get(Form("hPi0_PbPb_cen6_SystEcale")) ;
    hpi0SystE->SetName(Form("hPHOS_pi0_pPb_cen%s_SystEscale",centext)) ;
    hpi0SystE->SetTitle(Form("pi0 yield, centrality %s%%, E scale sys. errors",centext)) ;


    
    //Mass, Data
    TH1D * hdataMass = (TH1D*)fMass->Get("mass1_GS_Emin3_All_cen6") ;
    hdataMass->SetName("mass_data_GS_All") ;
    hdataMass->SetTitle("Peak position in data, Gaus fit, no PID") ;
    
    TH1D * hdataWidth = (TH1D*)fMass->Get("width1_GS_Emin3_All_cen6") ;
    hdataWidth->SetName("width_data_GS_All") ;
    hdataWidth->SetTitle("Peak width in data, Gaus fit, no PID") ;
    
    //MC
    TH1D * hMCMass = (TH1D*)fMCMass->Get("mass1_GS_Emin3_All_cen6") ;
    hMCMass->SetName("mass_MC_GS_All") ;
    hMCMass->SetTitle("Peak position in MC, Gaus fit, no PID") ;
    
    TH1D * hMCWidth = (TH1D*)fMCMass->Get("width1_GS_Emin3_All_cen6") ;
    hMCWidth->SetName("width_MC_GS_All") ;
    hMCWidth->SetTitle("Peak width in MC, Gaus fit, no PID") ;
    
    
    //Efficiency
    TF1 * feff = (TF1*)fileeff->Get(Form("eff_Pi0_Gaus_pPb_Emin3_All_cen6")) ; 
    feff->SetRange(1.,20.) ;
    feff->SetName("efficiency_GS_All") ;
    TH1D * heff=(TH1D*)fileeff->Get(Form("yeild_GS_Emin3_All_cen6")) ;
    heff->SetName("hefficiency_GS_All") ;
    heff->SetTitle("pi0 reconstruction efficiency, Gaus fit, no PID") ;
    
    
    //Invariant mass distributions
    TFile *fmgg = new TFile("mgg.root") ;
    TH1D * hMggReal = (TH1D*)fmgg->Get("real") ;
    hMggReal->SetName("hMggReal") ;
    hMggReal->SetTitle("Invariant mass distribution, Real") ;
    
    TH1D * hMggSignal = (TH1D*)fmgg->Get("Signal") ;
    hMggSignal->SetName("hMggSignal") ;
    hMggSignal->SetTitle("Invariant mass distribution, Signal") ;
    
    TH1D * hMggMixed = (TH1D*)fmgg->Get("mixed") ;
    hMggMixed->SetName("hMggBg") ;
    hMggMixed->SetTitle("Invariant mass distribution, background from mixed events") ;
    
    
    TH1D * hMggMixedNP = (TH1D*)fmgg->Get("BgPolCorrection") ;
    hMggMixedNP->SetName("hMggBgNoPolCorrection") ;
    hMggMixedNP->SetTitle("Invariant mass distribution, background from mixed events, no polinomial correction") ;

    
    TFile * ffeed = new TFile("feed_down.root") ;
    TF1 * ftmp = (TF1*)ffeed->Get("fun1") ;
    TF1 * feeddown = new TF1("FeeddownCorrection",p4,0.9,20.,5) ;
    feeddown->SetParameters(ftmp->GetParameters()) ;

    
    Int_t nbin=33 ;
    Double_t xa[34] ={0.8,1.0,1.2,1.4,1.6, 1.8,2.0,2.2,2.4,2.6, 2.8,3.0,3.2,3.4,3.6, 3.8,4.0,4.5,5.0,5.5, 6.,7.,8.,10.,12.,16.,20.,22.,  24.,26.,28.,30.,35.,  40.};
    TH1F * hSystTotal= new TH1F("hSystTotal","Summary of total sys err.",nbin,xa);
    for(Int_t i=1; i<=nbin; i++){
      if(hpi0Syst->GetBinContent(i)>0.)
        hSystTotal->SetBinContent(i,hpi0Syst->GetBinError(i)/hpi0Syst->GetBinContent(i)) ;  
        
    }
    
    
    //Writing
    TDirectory *dir = fout->mkdir(Form("PHOS_pPb_502_Centrality_%s",centext),
                                  Form("PHOS pi0 spectrum pPb 5.02 TeV, centrality %s%%",centext));
    dir->cd();
 
    //pi0s
    hpi0Stat->Write(0,TObject::kOverwrite) ;  
    hpi0Syst->Write(0,TObject::kOverwrite) ;      
    hpi0SystA->Write(0,TObject::kOverwrite) ;  
    hpi0SystB->Write(0,TObject::kOverwrite) ;  
    hpi0SystC->Write(0,TObject::kOverwrite) ;  
    
    hpi0SystSig->Write(0,TObject::kOverwrite) ; 
    hpi0SystBg->Write(0,TObject::kOverwrite) ; 
    hpi0SystAcc->Write(0,TObject::kOverwrite) ; 
    hpi0SystMat->Write(0,TObject::kOverwrite) ; 
    hpi0SystNL->Write(0,TObject::kOverwrite) ; 
    hpi0SystE->Write(0,TObject::kOverwrite) ; 
    
    
    hdataMass->Write(0,TObject::kOverwrite) ;  
    hdataWidth->Write(0,TObject::kOverwrite) ;  
    hMCMass->Write(0,TObject::kOverwrite) ;  
    hMCWidth->Write(0,TObject::kOverwrite) ;  
    feff->Write(0,TObject::kOverwrite) ;  
    heff->Write(0,TObject::kOverwrite) ;  
    feeddown->Write(0,TObject::kOverwrite) ; 
 
    hMggReal->Write(0,TObject::kOverwrite) ;
    hMggSignal->Write(0,TObject::kOverwrite) ;
    hMggMixed->Write(0,TObject::kOverwrite) ;
    hMggMixedNP->Write(0,TObject::kOverwrite) ;
    
    //a-la Tsubasa
    fout->cd() ;
    hSystTotal->Write(0,TObject::kOverwrite) ;
    hpi0Stat->Write("hCor_stat",TObject::kOverwrite) ;  
    hpi0Syst->Write("hCor_syst",TObject::kOverwrite) ;     
    hdataMass->Write("Gmean_Real",TObject::kOverwrite) ;  
    hMCMass->Write("Gmean_MC",TObject::kOverwrite) ;  
    hdataWidth->Write("Gsigma_Real",TObject::kOverwrite) ;  
    hMCWidth->Write("Gsigma_MC",TObject::kOverwrite) ;  
   
    
    
    
  
   
}

//Parameterization from Boris
Double_t p4(Double_t * x, Double_t * par)
{
   Double_t z=TMath::Min(7.,x[0]);
    
    Double_t f = par[0] + par[1]*z + par[2]*z*z + par[3]*z*z*z + par[4]*z*z*z*z;
    return f;
}
