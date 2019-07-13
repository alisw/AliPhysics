{
   const char* infile = "combinedptresponse.root";
   const char* fitfile = "spectrafits.root";
   
   const char* outfile = "ptresponsematrix_LHC17pq_pass1.root";
   
   
   TFile* fin = TFile::Open(infile,"READ");
   
   TH2D* ptresponse = (TH2D*) fin->Get("ptresponse");
   
   
   TFile* ffit = TFile::Open(fitfile,"READ");
   
   TF1* fit_lomult = (TF1*) ffit->Get("fitsmallnch");
   TF1* fit_himult = (TF1*) ffit->Get("fitlargench");
   TH2D* mario = (TH2D*) ffit->Get("marioresposnematrix");
   
   // these are rec vs. gen
   TH2D* ptresp_himult =(TH2D*) mario->Clone("ptresp_himult");
   TH2D* ptresp_lomult =(TH2D*) mario->Clone("ptresp_lomult");
   
   //loop over gen pt in step of the ptresponse
   for (int i=1; i<=ptresponse->GetXaxis()->GetNbins(); ++i) {
       Double_t ptgen = ptresponse->GetXaxis()->GetBinCenter(i);
       Double_t pt1 =  ptresponse->GetXaxis()->GetBinLowEdge(i);
       Double_t pt2 = pt1 + ptresponse->GetXaxis()->GetBinWidth(i);
       Double_t w_lomult = fit_lomult->Integral(pt1,pt2);
       Double_t w_himult = fit_himult->Integral(pt1,pt2);
       // loop over rec pt
       for (int j=1; j<=ptresponse->GetYaxis()->GetNbins(); ++j) {    
           Double_t ptrec = ptresponse->GetYaxis()->GetBinCenter(j);
           Double_t w = ptresponse->GetBinContent(i,j);
           if (w == 0) continue;
           ptresp_lomult->Fill(ptrec,ptgen,w_lomult*w);
           ptresp_himult->Fill(ptrec,ptgen,w_himult*w);       
       }
   }

   TH2D* ptresp_himult_norm =(TH2D*) ptresp_himult->Clone("ptresp_himult_norm");
   TH2D* ptresp_lomult_norm =(TH2D*) ptresp_lomult->Clone("ptresp_lomult_norm");   
   
//normalize to 1  in each ptgen
for (int i=1; i<=ptresp_himult_norm->GetNbinsY(); i++) {
    Double_t integral = ptresp_himult_norm->Integral(1,ptresp_himult_norm->GetNbinsX(),i,i);
    if (integral<=0) continue;
    for (int j=0; j<=ptresp_himult_norm->GetNbinsX(); j++) { 
        ptresp_himult_norm->SetBinContent(j,i,ptresp_himult_norm->GetBinContent(j,i)/integral);
    }
}
   
//normalize to 1  in each ptgen
for (int i=1; i<=ptresp_lomult_norm->GetNbinsY(); i++) {
    Double_t integral = ptresp_lomult_norm->Integral(1,ptresp_lomult_norm->GetNbinsX(),i,i);
    if (integral<=0) continue;
    for (int j=0; j<=ptresp_lomult_norm->GetNbinsX(); j++) { 
        ptresp_lomult_norm->SetBinContent(j,i,ptresp_lomult_norm->GetBinContent(j,i)/integral);
    }
}
   TFile* fout = TFile::Open(outfile,"RECREATE");
ptresp_lomult->Write();
ptresp_himult->Write();
ptresp_himult_norm->Write();
ptresp_lomult_norm->Write();
fout->Close();
ptresponse->Draw("COLZ");    
    
    
}