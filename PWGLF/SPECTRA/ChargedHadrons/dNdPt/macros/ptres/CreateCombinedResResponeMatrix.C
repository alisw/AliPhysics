{
   const char* outfile = "combinedptresponse.root";
   
   Int_t nperbin = 1000000;
   Double_t resfactor = 1.0; 
   Double_t extrasmear = 0.0;
   Double_t meanshifta = 0.0;
   Double_t sigmashifta = 0.0; 
   Double_t meanshiftc = 0.0;
   Double_t sigmashiftc = 0.0;
   
   TH2D* ptresponse = new TH2D("ptresponse","ptresponse",4000,0,60,4000,0,60);
   ptresponse->GetXaxis()->SetTitle("pT,gen");
   ptresponse->GetYaxis()->SetTitle("pT,rec");
   
   TH2D* ptresponsemc = new TH2D("ptresponsemc","ptresponsemc",4000,0,60,4000,0,60);
   ptresponsemc->GetXaxis()->SetTitle("pT,gen");
   ptresponsemc->GetYaxis()->SetTitle("pT,rec");
  
   TH2D* ptresponsedata = new TH2D("ptresponsedata","ptresponsedata",4000,0,60,4000,0,60);
   ptresponsedata->GetXaxis()->SetTitle("pT,gen");
   ptresponsedata->GetYaxis()->SetTitle("pT,rec");   
  
    const char* mcfilename = "/home/mknichel/charged-particle-spectra/train_output/LF_pp_MC/1101_LHC17l3b/AnalysisResults.root";
    const char* datafilename = "/home/mknichel/charged-particle-spectra/train_output/LF_pp/1227_LHC17pq_pass1/AnalysisResults.root";

    TFile * fmc = TFile::Open(mcfilename,"READ");
    TFile * fdata = TFile::Open(datafilename,"READ");

    THnSparseD* hMCCov   = (THnSparseD*)((TList*)((TDirectoryFile*)fmc->Get("TaskPtResStudy"))->Get("TaskPtResStudy"))->FindObject("fHistPtResCov");
    THnSparseD* hDataCov     = (THnSparseD*)((TList*)((TDirectoryFile*)fdata->Get("TaskPtResStudy"))->Get("TaskPtResStudy"))->FindObject("fHistPtResCov");

    THnSparseD* hMCRes   = (THnSparseD*)((TList*)((TDirectoryFile*)fmc->Get("TaskPtResStudy"))->Get("TaskPtResStudy"))->FindObject("fHistPtRes");
    THnSparseD* hDataRes     = (THnSparseD*)((TList*)((TDirectoryFile*)fdata->Get("TaskPtResStudy"))->Get("TaskPtResStudy"))->FindObject("fHistPtRes");

    THnSparseD* hMCResMC   = (THnSparseD*)((TList*)((TDirectoryFile*)fmc->Get("TaskPtResStudy"))->Get("TaskPtResStudy"))->FindObject("fHistPtResMC");
    
    TH2D* mcptresponse = (TH2D*) hMCResMC->Projection(0,1)->Clone("mcptresponse");
   
    TH2D* cov = (TH2D*) hDataCov->Projection(1,0)->Clone("cov");
    
    TH1D* scaledcov = new TH1D("scaledcov","scaledcov",1000,0,0.1);
    
    //create combined histogram pos+neg
    for (int i=1; i<=200; i++) {Double_t pt1 = cov->GetXaxis()->GetBinCenter(i); if (pt1>0) continue; for (int y=1; y<=1000; y++) { cov->Fill(-pt1,cov->GetYaxis()->GetBinCenter(y),cov->GetBinContent(i,y)); } }
    cov->GetXaxis()->SetRangeUser(0,0.19);
    TH1D* prof = (TH1D*) cov->ProfileX();
    
    TF1* fit = new TF1("linfit","[0]+x*[1]+x*x*[2]",0,1);
    TF1* normres = new TF1("normres","x*[0]+x*x*[1]",0,1);    
    
    prof->Fit(fit);prof->Fit(fit);prof->Fit(fit);
    normres->SetParameters(fit->GetParameter(1),fit->GetParameter(2));    

    for (int i=1; i<=200; i++) {Double_t pt1 = cov->GetXaxis()->GetBinCenter(i); if ((pt1<0) || (pt1>0.2)) continue; for (int y=1; y<=1000; y++) { Double_t s1pt = cov->GetYaxis()->GetBinCenter(y); scaledcov->Fill(s1pt-normres->Eval(pt1),cov->GetBinContent(i,y));  } }

   delete gRandom;
   gRandom = new TRandom3(0);
   printf("hier!\n");
   for(int x=1; x <= mcptresponse->GetNbinsX(); x++) {
       for (int y=1; y <= mcptresponse->GetNbinsY(); y++) {
           ptresponsemc->SetBinContent(x,y,mcptresponse->GetBinContent(x,y));
       }
   }
   
   
   for (int i=1; i<=ptresponsedata->GetNbinsX(); i++) {
       Double_t pt = ptresponsedata->GetXaxis()->GetBinCenter(i); 
       if (pt<5) continue;
       Double_t pt1 = 1./pt;
       for (int j=0; j<nperbin; j++) {
            Double_t sspt1 = scaledcov->GetRandom();
            Double_t spt1 = sspt1 + normres->Eval(pt1);
            spt1 *= resfactor;
            if (extrasmear!=0) spt1 = TMath::Sqrt(spt1*spt1+extrasmear*extrasmear);
            int side = gRandom->Integer(2);
            int charge = gRandom->Integer(2);
            Double_t dpt1 = 0;        
            if (side)    { dpt1 = gRandom->Gaus(meanshifta,sigmashifta); }
                    else { dpt1 = gRandom->Gaus(meanshiftc,sigmashiftc); }            
            if (charge) { pt1 += dpt1; }
               else { pt1 -= dpt1; }      
        Double_t pt1s = gRandom->Gaus(pt1,spt1);
        Double_t pts = 1./pt1s; 
        ptresponsedata->Fill(pt,pts);
       }
   }
       
    
//normalize to 1  in each row
for (int i=1; i<=ptresponsedata->GetNbinsX(); i++) {
    Double_t integral = ptresponsedata->Integral(i,i,1,ptresponsedata->GetNbinsY());
    if (integral<1) continue;
    for (int j=0; j<=ptresponsedata->GetNbinsY(); j++) { 
        ptresponsedata->SetBinContent(i,j,ptresponsedata->GetBinContent(i,j)/integral);
    }
}    

//normalize to 1  in each row
for (int i=1; i<=ptresponsemc->GetNbinsX(); i++) {
    Double_t integral = ptresponsemc->Integral(i,i,1,ptresponsemc->GetNbinsY());
    if (integral<1) continue;
    for (int j=0; j<=ptresponsemc->GetNbinsY(); j++) { 
        ptresponsemc->SetBinContent(i,j,ptresponsemc->GetBinContent(i,j)/integral);
    }
}    

//combine data and mc
for (int i=1; i<=ptresponse->GetNbinsX(); i++) {        
    for (int j=0; j<=ptresponse->GetNbinsY(); j++) { 
        Double_t ptgen =  ptresponse->GetXaxis()->GetBinCenter(i); 
        Double_t mc    = ptresponsemc->GetBinContent(i,j);
        Double_t data  = ptresponsedata->GetBinContent(i,j);
        Double_t comb  = 0;
        if (ptgen<5) { comb = mc; }
        else if (ptgen>9) { comb = data; }       
        else { comb = 0.25*(ptgen-5)*data + 0.25*(9-ptgen)*mc; }
        ptresponse->SetBinContent(i,j,comb);        
        ptresponse->SetBinError(i,j,0);
    }   
} 

//normalize again, just in case something went wrong, set errors to zero
for (int i=1; i<=ptresponse->GetNbinsX(); i++) {
    Double_t integral = ptresponse->Integral(i,i,1,ptresponse->GetNbinsY());
    if (integral == 0) continue;
    for (int j=0; j<=ptresponse->GetNbinsY(); j++) { 
        ptresponse->SetBinContent(i,j,ptresponse->GetBinContent(i,j)/integral);
        ptresponse->SetBinError(i,j,0);
    }   
}        
  
TFile* fout = TFile::Open(outfile,"RECREATE");
ptresponse->Write();
fout->Close();
ptresponse->Draw("COLZ");    
    
}