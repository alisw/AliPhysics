/***************************************************
processJETriggerQA:
To procees JE Triggere QA wagon's output

****************************************************/



const Float_t ptmin =  0. ; //lower cutoff of jet pt spectrum

void processJETriggerQA_V2(TString strFileIn    = "AnalysisResults.root", 
			   TString suftype      = "eps",
			   Float_t jetR         = 0.2, 
			   Float_t minTrkPT     = 0.15, 
			   Float_t minClusterET = 0.3, 
			   Int_t run            = 0,
			   TString trigsuffix   = "",
			   const char* outfile  = "JETriggerQA_outfile.root"
			   ){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TString prefix = "fig_je_TriggerQA_";

  TFile * f1 = TFile::Open(strFileIn.Data());

  if(!trigsuffix=="" || !trigsuffix=="EJE" || !trigsuffix=="EGA")
  {cout << "Unknown trigger suffix. It should be either empty, EJE or EGA." << endl; return;}

  //Load histogram list
  TString folder = Form("TriggerQA_Jet_AKTFullR%02d0_PicoTracks_pT%04d_CaloClustersCorr_ET%04d_pt_scheme_Jet_AKTChargedR%02d0_PicoTracks_pT%04d_CaloClustersCorr_ET%04d_pt_scheme_TC%s/TriggerQA_Jet_AKTFullR%02d0_PicoTracks_pT%04d_CaloClustersCorr_ET%04d_pt_scheme_Jet_AKTChargedR%02d0_PicoTracks_pT%04d_CaloClustersCorr_ET%04d_pt_scheme_TC%s", 
       TMath::Nint(jetR*10), TMath::Nint(minTrkPT*1000), TMath::Nint(minClusterET*1000),
       TMath::Nint(jetR*10), TMath::Nint(minTrkPT*1000), TMath::Nint(minClusterET*1000),trigsuffix.Data(),
       TMath::Nint(jetR*10), TMath::Nint(minTrkPT*1000), TMath::Nint(minClusterET*1000),  
       TMath::Nint(jetR*10), TMath::Nint(minTrkPT*1000), TMath::Nint(minClusterET*1000),trigsuffix.Data());
  TList *histList = 0x0;
  histList =  (TList*)f1->Get(folder);
  if(histList==0){
    cout << "Could not find " << folder << " in " << strFileIn << endl;
    return;
  }

  //---------------------------------------------------------------------------------------------------
  //       jet histograms
  //---------------------------------------------------------------------------------------------------

  const Int_t kJetType = 2;
  TString suffix[kJetType] = {"Charged","Full"};
  TH3F *h3PtEtaPhiJet[kJetType];
  TH1F *hPtJet[kJetType];
  TH2F *hEtaPhiJet[kJetType];
  TH2F *hRhoCent[kJetType];

  TH1F *hNEventSel = histList->FindObject("fhNEvents");
  Float_t nEvents  = hNEventSel->GetBinContent(2);

  TProfile *hTriggerbit = histList->FindObject("fhTriggerbit");
  UInt_t triggerbit = hTriggerbit->GetBinContent(1);

  for(Int_t itype = 0; itype < kJetType; itype++){

     h3PtEtaPhiJet[itype] = (TH3F*)  histList->FindObject(Form("fh3PtEtaPhiJet%s",suffix[itype].Data()));
     if(! h3PtEtaPhiJet[itype]) continue;

     //jet pt spectra
     Int_t binMin = 1;
     if(ptmin>0.) binMin = h3PtEtaPhiJet[itype]->GetXaxis()->FindBin(ptmin+0.00001);
     h3PtEtaPhiJet[itype]->GetXaxis()->SetRange(binMin, h3PtEtaPhiJet[itype]->GetNbinsX());

     hPtJet[itype] = (TH1F*) h3PtEtaPhiJet[itype]->Project3D("x");
     hPtJet[itype]->Scale(1./nEvents,"width");
     SetHist((TH1F*) hPtJet[itype],"p_{T,corr}^{jet} (GeV)","1/N_{evt} dN/dp_{T,corr}^{jet} (GeV^{-1})");
     hPtJet[itype]->SetName(Form("hPtJet%s",suffix[itype].Data())); 

     //eta versus phi
     hEtaPhiJet[itype] = (TH2F*) h3PtEtaPhiJet[itype]->Project3D("yz");
     SetHist((TH1F*) hEtaPhiJet[itype],"#varphi^{jet} (rad)","#eta^{jet}");
     hEtaPhiJet[itype]->SetName(Form("hEtaPhiJet%s",suffix[itype].Data()));
  }

  for(Int_t itype = 0; itype < kJetType; itype++){
  //rho versus centrality
    hRhoCent[itype] = (TH2F*) histList->FindObject(Form("fHistRhovsCent%s",suffix[itype].Data()));
     if(!hRhoCent[itype]) continue;
     SetHist((TH1F*) hRhoCent[itype], "Centrality (%)", Form("%s#rho (GeV/c*rad^{-1})", (suffix[itype]=="Full")?"s":"" ));
     hRhoCent[itype]->SetName(Form("hRhoCent%s",suffix[itype].Data()));
  }

  //______________
  //Draw histograms

  TCanvas *c[100];
  TH1F *frame[100];
  Int_t nCan = 0;
  TLegend *leg; 

  TString longtrigname = "";
  UInt_t emcalPhysSelbit = 1<<31;
  if(triggerbit==AliVEvent::kAny) longtrigname = "kAny";
  else if(triggerbit==AliVEvent::kAnyINT) longtrigname = "kAnyINT";
  else if(triggerbit==emcalPhysSelbit) longtrigname = "EmcalPhysicsSelectionTask";
  else {
    if(triggerbit & AliVEvent::kCentral) longtrigname += "kCentral";
    if(triggerbit & AliVEvent::kSemiCentral) longtrigname += "kSemiCentral";
    if(triggerbit & AliVEvent::kMB)  longtrigname += "kMB";
  }

  if(trigsuffix=="EJE")  longtrigname += "AndkEMCEJE";
  if(trigsuffix=="EGA")  longtrigname += "AndkEMCEGA";
 
  for(Int_t itype = 0; itype < kJetType; itype++){ //loop over charged and full jets

     //draw pt spectrum
     if(!hPtJet[itype]) continue; 
     c[nCan] = new TCanvas(Form("c%d",nCan),Form("c%d: Pt %s jets",nCan,suffix[itype].Data()),600,450);
     SetCanvas((TCanvas*) c[nCan]);
     c[nCan]->SetLogy();

     frame[nCan] = gPad->DrawFrame(hPtJet[itype]->GetBinLowEdge(1),  
                                   1e-7,
                                   hPtJet[itype]->GetBinLowEdge(hPtJet[itype]->GetNbinsX()+1), 
                                   hPtJet[itype]->GetBinContent(hPtJet[itype]->GetMaximumBin())*2.);

     SetHist((TH1F*) frame[nCan],hPtJet[itype]->GetXaxis()->GetTitle(),hPtJet[itype]->GetYaxis()->GetTitle());
 

     hPtJet[itype]->DrawCopy("same");
 
     leg = new TLegend(0.15,0.5,0.88,0.88);
     SetLeg(leg);
     
     TString txt =  Form("%s jets AKT R=%.1f",suffix[itype].Data(),jetR);
     if(run>0) txt += Form(" run:%d",run); 
     leg->AddEntry((TObject*) 0, txt.Data(),"");
     leg->AddEntry((TObject*) 0, Form("Trigger: %s",longtrigname.Data()),"");
     leg->AddEntry((TObject*) 0, Form("p_{T,trk}> %d MeV",TMath::Nint(minTrkPT*1000)),"");
     if(itype==1) leg->AddEntry((TObject*) 0, Form("E_{T}>%d MeV",TMath::Nint(minClusterET*1000)),""); 
     leg->AddEntry((TObject*) 0, Form("#it{N}_{events} = %.0f",nEvents),"");
     leg->Draw();

     c[nCan]->SaveAs(Form("%s_Pt_AKT%02d_pT%04d_ET%04d_Run%d_Trigger%s_%s.%s",prefix.Data(),TMath::Nint(jetR*10),
			  TMath::Nint(minTrkPT*1000),TMath::Nint(minClusterET*1000), run, longtrigname.Data(), suffix[itype].Data(), suftype.Data()));

     nCan++;
  }//end of loop over charged and full jets
  //_________________
  //draw eta versus phi
  for(Int_t itype = 0; itype < kJetType; itype++){ //loop over charged and full jets
     if(!hEtaPhiJet[itype]) continue; 
     c[nCan] = new TCanvas(Form("c%d",nCan),Form("c%d: eta-phi %s jets",nCan,suffix[itype].Data()),600,450);
     SetCanvas((TCanvas*) c[nCan]);
     c[nCan]->SetRightMargin(0.15);

     frame[nCan] = gPad->DrawFrame(hEtaPhiJet[itype]->GetXaxis()->GetBinLowEdge(1),  
                                   hEtaPhiJet[itype]->GetYaxis()->GetBinLowEdge(1),
                                   hEtaPhiJet[itype]->GetXaxis()->GetBinLowEdge(hEtaPhiJet[itype]->GetNbinsX()),
                                   hEtaPhiJet[itype]->GetYaxis()->GetBinLowEdge(hEtaPhiJet[itype]->GetNbinsY()));

     SetHist((TH1F*) frame[nCan],hEtaPhiJet[itype]->GetXaxis()->GetTitle(),hEtaPhiJet[itype]->GetYaxis()->GetTitle());
 

     hEtaPhiJet[itype]->DrawCopy("same,colz");
 
     leg = new TLegend(0.15,0.5,0.88,0.88);
     SetLeg(leg);

     TString txt =  Form("%s jets AKT R=%.1f",suffix[itype].Data(),jetR);
     if(run>0) txt += Form(" run:%d",run); 
     leg->AddEntry((TObject*) 0, txt.Data(),"");
     leg->AddEntry((TObject*) 0, Form("Trigger: %s",longtrigname.Data()),"");
     leg->AddEntry((TObject*) 0, Form("p_{T,trk}> %d MeV",TMath::Nint(minTrkPT*1000)),"");
     if(itype==1) leg->AddEntry((TObject*) 0, Form("E_{T}>%d MeV",TMath::Nint(minClusterET*1000)),"");
     leg->AddEntry((TObject*) 0, Form("#it{N}_{events} = %.0f",nEvents),"");
     leg->AddEntry((TObject*) 0, Form("p_{T,corr}^{jet} > %.1f GeV", ptmin),"");
     leg->Draw();

     c[nCan]->SaveAs(Form("%s_EtaPhi_AKT%02d_pT%04d_ET%04d_Run%d_Trigger%s_%s.%s",prefix.Data(),TMath::Nint(jetR*10),
			  TMath::Nint(minTrkPT*1000),TMath::Nint(minClusterET*1000), run, longtrigname.Data(), suffix[itype].Data(), suftype.Data()));

     nCan++;
  }//end of the loop over charged and full jets
  //_________________
  //draw rho versus centrality
  for(Int_t itype = 0; itype < kJetType; itype++){ //loop over charged and full jets
     if(hRhoCent[itype]) {
       c[nCan] = new TCanvas(Form("c%d",nCan),Form("c%d: Rho-Cent",nCan),600,450);
       SetCanvas((TCanvas*) c[nCan]);
       c[nCan]->SetRightMargin(0.15);
       c[nCan]->SetLogz();
  
       frame[nCan] = gPad->DrawFrame(hRhoCent[itype]->GetXaxis()->GetBinLowEdge(1),
                                     hRhoCent[itype]->GetYaxis()->GetBinLowEdge(1),
                                     hRhoCent[itype]->GetXaxis()->GetBinLowEdge(hRhoCent[itype]->GetNbinsX()),
                                     hRhoCent[itype]->GetYaxis()->GetBinLowEdge(hRhoCent[itype]->GetNbinsY()));
  
       SetHist((TH1F*) frame[nCan],hRhoCent[itype]->GetXaxis()->GetTitle(),hRhoCent[itype]->GetYaxis()->GetTitle());
  
       hRhoCent[itype]->DrawCopy("colz");
  
       leg = new TLegend(0.15,0.5,0.88,0.88);
       SetLeg(leg);
  
       TString txt =  Form("%s jets AKT R=%.1f",suffix[itype].Data(),jetR);
       if(run>0) txt += Form(" run:%d",run); 
  
       leg->AddEntry((TObject*) 0, txt.Data(),"");
       leg->AddEntry((TObject*) 0, Form("Trigger: %s",longtrigname.Data()),"");
       leg->AddEntry((TObject*) 0, Form("#it{N}_{events} = %.0f",nEvents),"");
       leg->Draw();
  
       c[nCan]->SaveAs(Form("%s_RhoCent_AKT%02d_pT%04d_ET%04d_Run%d_Trigger%s_%s.%s",prefix.Data(),TMath::Nint(jetR*10),
			    TMath::Nint(minTrkPT*1000),TMath::Nint(minClusterET*1000), run, longtrigname.Data(), suffix[itype].Data(), suftype.Data()));
       nCan++; 
    }
  }//end of the loop over charged and full jets

  //---------------------------------------------------------------------------------------------------
  //                       WRITE OUTPUT TO ROOT FILE
  //---------------------------------------------------------------------------------------------------


  /* Standalone output
    TFile *histOut = new TFile(Form("%s_AKT%02d_pT%04d_ET%04d_jetPtMin%.1f_Run%d.root",prefix.Data(),
    TMath::Nint(jetR*10), TMath::Nint(minTrkPT*1000), TMath::Nint(minClusterET*1000),ptmin,run),"RECREATE");
    
    for(Int_t itype = 0; itype < kJetType; itype++){ //loop over charged and full jets
     
    if(hPtJet[itype])     hPtJet[itype]->Write();
    if(hEtaPhiJet[itype]) hEtaPhiJet[itype]->Write();
    if(hRhoCent[itype])   hRhoCent[itype]->Write();
    }
    
    histOut->Close();
  
  */

  // Common output - 
  // Added by sjena

  TFile *fout = TFile::Open(outfile,"UPDATE");
  fout->ls();
  
  TDirectoryFile *cdd = NULL;
  cdd = (TDirectoryFile*)fout->Get("JE");
  if(!cdd) {
    Printf("Warning: JE <dir> doesn't exist, creating a new one");
    cdd = (TDirectoryFile*)fout->mkdir("JE");
  }
  cdd->cd();
  cdd->ls();
  
  for(Int_t itype = 0; itype < kJetType; itype++){ //loop over charged and full jets
    
    if(hPtJet[itype])     hPtJet[itype]->Write(Form("%s%d_%s",prefix.Data(), itype, hPtJet[itype]->GetName()));
    if(hEtaPhiJet[itype]) hEtaPhiJet[itype]->Write(Form("%s%d_%s",prefix.Data(), itype, hEtaPhiJet[itype]->GetName()));
    if(hRhoCent[itype])   hRhoCent[itype]->Write(Form("%s%d_%s",prefix.Data(), itype, hRhoCent[itype]->GetName()));
  }
  
  fout->Close();



}
//__________________________________________________________

void SetHist(TH1* h,TString titx, TString tity){

   h->GetXaxis()->SetTitle(titx.Data());
   h->GetYaxis()->SetTitle(tity.Data());
   h->GetXaxis()->SetTitleSize(0.06);
   h->GetYaxis()->SetTitleSize(0.06);
   h->GetYaxis()->SetTitleOffset(1.);
   h->GetXaxis()->SetTitleOffset(1.);
   h->SetLineWidth(3);

}
//_____________________________________________________________________

void SetCanvas(TCanvas* c){
   c->SetLeftMargin(0.15);
   c->SetBottomMargin(0.15);
   c->SetRightMargin(0.05);
   c->SetTopMargin(0.05);
   c->SetTickx();
   c->SetTicky();
}
//_____________________________________________________________________

void SetLeg(TLegend* le){
   le->SetFillColor(10);
   le->SetBorderSize(0);
   le->SetFillStyle(0);
   le->SetTextSize(0.05);
}
