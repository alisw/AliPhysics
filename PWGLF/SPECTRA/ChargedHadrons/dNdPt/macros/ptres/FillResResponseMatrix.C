{
   const char* infile = "scaledcov.root";
   const char* outfile = "ptresponse.root";
   
   Int_t nperbin = 1000000;
   Double_t resfactor = 1.0; 
   Double_t extrasmear = 0.0;
   Double_t meanshifta = 0.0;
   Double_t sigmashifta = 0.0; 
   Double_t meanshiftc = 0.0;
   Double_t sigmashiftc = 0.0;
    
   TH2D* ptresponse = new TH2D("ptresponse","ptresponse",4000,0,200,4000,0,200);
   ptresponse->GetXaxis()->SetTitle("pT,gen");
   ptresponse->GetYaxis()->SetTitle("pT,rec");
 
   TFile* fc = TFile::Open(infile,"READ");
   TH1F* sc = (TH1F*) fc->Get("scaledcov");
   TF1*  nr = (TF1*)  fc->Get("normres"); 
   
   delete gRandom;
   gRandom = new TRandom3(0);
   
   for (int i=1; i<=ptresponse->GetNbinsX(); i++) {
       Double_t pt = ptresponse->GetXaxis()->GetBinCenter(i);
       if (pt<7) continue;
       Double_t pt1 = 1./pt;
       for (int j=0; j<nperbin; j++) {
            Double_t sspt1 = sc->GetRandom();
            Double_t spt1 = sspt1 + nr->Eval(pt1);
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
        ptresponse->Fill(pt,pts);
       }
   }
    
//normalize to 1  in each row
for (int i=1; i<=ptresponse->GetNbinsX(); i++) {
    Double_t integral = ptresponse->Integral(i,i,1,ptresponse->GetNbinsY());
    if (integral<1) continue;
    for (int j=0; j<=ptresponse->GetNbinsY(); j++) { 
        ptresponse->SetBinContent(i,j,ptresponse->GetBinContent(i,j)/integral);
    }
}
    
  
TFile* fout = TFile::Open(outfile,"RECREATE");
ptresponse->Write();
fout->Close();
ptresponse->Draw("COLZ");
}
