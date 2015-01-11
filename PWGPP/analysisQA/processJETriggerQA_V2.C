/***************************************************
processJETriggerQA:
To procees JE Triggere QA wagon's output

****************************************************/



const Float_t ptmin =  0. ; //lower cutoff of jet pt spectrum

void processJETriggerQA_V2(TString strFileIn    = "AnalysisResults.root", 
                           TString strFileIn2   = "",
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
  Bool_t drawComp = kFALSE;
  if(strFileIn2!="") {
    TFile * f2 = TFile::Open(strFileIn2.Data());
    drawComp = kTRUE;
  }

  if(!trigsuffix=="" || !trigsuffix=="EJE" || !trigsuffix=="EGA")
  {cout << "Unknown trigger suffix. It should be either empty, EJE or EGA." << endl; return;}

  //Load histogram list
  TString folder = Form("TriggerQA_Jet_AKTFullR%02d0_PicoTracks_pT%04d_CaloClustersCorr_ET%04d_pt_scheme_Jet_AKTChargedR%02d0_PicoTracks_pT%04d_CaloClustersCorr_ET%04d_pt_scheme_TC%s/TriggerQA_Jet_AKTFullR%02d0_PicoTracks_pT%04d_CaloClustersCorr_ET%04d_pt_scheme_Jet_AKTChargedR%02d0_PicoTracks_pT%04d_CaloClustersCorr_ET%04d_pt_scheme_TC%s", 
       TMath::Nint(jetR*10), TMath::Nint(minTrkPT*1000), TMath::Nint(minClusterET*1000),
       TMath::Nint(jetR*10), TMath::Nint(minTrkPT*1000), TMath::Nint(minClusterET*1000),trigsuffix.Data(),
       TMath::Nint(jetR*10), TMath::Nint(minTrkPT*1000), TMath::Nint(minClusterET*1000),  
       TMath::Nint(jetR*10), TMath::Nint(minTrkPT*1000), TMath::Nint(minClusterET*1000),trigsuffix.Data());
  TList *histList[2];
  histList[0] =  (TList*)f1->Get(folder);
  if(histList[0]==0){
    cout << "Could not find " << folder << " in " << strFileIn << endl;
    return;
  }
  //Load histogram list for comparison
  if(drawComp) {
    TList *histList[1] = 0x0;
    histList[1] =  (TList*)f2->Get(folder);
    if(histList==0){
      cout << "Could not find " << folder << " in " << strFileIn2 << endl;
      return;
    }
  }

  //---------------------------------------------------------------------------------------------------
  //       jet histograms
  //---------------------------------------------------------------------------------------------------

  const Int_t kJetType = 2;
  TString suffix[kJetType] = {"Charged","Full"};
  TString suffix2[2] = {"","_old"};
  TH3F *h3PtEtaPhiJet[kJetType][2];
  TH1F *hPtJet[kJetType][2];
  TH2F *hEtaPhiJet[kJetType][2];
  TH2F *hRhoCent[kJetType][2];
  //Number [0] is for the original data, [1] is for the comparison graphs.

  TH1F *hNEventSel[2]; Float_t nEvents[2];
  hNEventSel[0] = (TH1F*) histList[0]->FindObject("fhNEvents");
  Float_t nEvents[0]  = hNEventSel[0]->GetBinContent(2);
  if(drawComp) {
    TH1F *hNEventSel[1] = (TH1F*) histList[1]->FindObject("fhNEvents");
    Float_t nEvents[1]  = hNEventSel[1]->GetBinContent(2);
  }

  TProfile *hTriggerbit = histList[0]->FindObject("fhTriggerbit");
  UInt_t triggerbit = hTriggerbit->GetBinContent(1);

  for(Int_t itype = 0; itype < kJetType; itype++){
    for(Int_t comp = 0; comp < 2; comp++){
      if((comp==1)&&(drawComp==kFALSE)) continue;
      h3PtEtaPhiJet[itype][comp] = (TH3F*)  histList[comp]->FindObject(Form("fh3PtEtaPhiJet%s",suffix[itype].Data()));
      if(! h3PtEtaPhiJet[itype][comp]) continue;
 
      //jet pt spectra
      Int_t binMin = 1;

      if(ptmin>0.) binMin = h3PtEtaPhiJet[itype][comp]->GetXaxis()->FindBin(ptmin+0.00001);
      h3PtEtaPhiJet[itype][comp]->GetXaxis()->SetRange(binMin, h3PtEtaPhiJet[itype][comp]->GetNbinsX());
 
      hPtJet[itype][comp] = (TH1F*) h3PtEtaPhiJet[itype][comp]->Project3D("x");
      hPtJet[itype][comp]->Scale(1./nEvents[comp],"width");
      SetHist((TH1F*) hPtJet[itype][comp],"p_{T,corr}^{jet} [GeV/c]","1/N_{evt} dN/dp_{T,corr}^{jet} [(GeV/c)^{-1}]");
      hPtJet[itype][comp]->SetName(Form("hPtJet%s%s",suffix[itype].Data(),suffix2[comp].Data())); 
 
      //eta versus phi
      hEtaPhiJet[itype][comp] = (TH2F*) h3PtEtaPhiJet[itype][comp]->Project3D("yz");
      SetHist((TH1F*) hEtaPhiJet[itype][comp],"#varphi^{jet} (rad)","#eta^{jet}");
      hEtaPhiJet[itype][comp]->SetName(Form("hEtaPhiJet%s%s",suffix[itype].Data(),suffix2[comp].Data()));
    }
  }

  for(Int_t itype = 0; itype < kJetType; itype++){
    for(Int_t comp = 0; comp < 2; comp++){
      if((comp==1)&&(drawComp==kFALSE)) continue;
      //rho versus centrality
      hRhoCent[itype][comp] = (TH2F*) histList[comp]->FindObject(Form("fHistRhovsCent%s",suffix[itype].Data()));
      if(!hRhoCent[itype][comp]) continue;
      SetHist((TH1F*) hRhoCent[itype][comp], "Centrality (%)", Form("%s#rho [GeV/c*rad^{-1}]", (suffix[itype]=="Full")?"s":"" ));
      hRhoCent[itype][comp]->SetName(Form("hRhoCent%s%s",suffix[itype].Data(),suffix2[comp].Data()));
    }
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

  if(trigsuffix=="EJE")  longtrigname += " and kEMCEJE";
  if(trigsuffix=="EGA")  longtrigname += " and kEMCEGA";
 
  for(Int_t itype = 0; itype < kJetType; itype++){ //loop over charged and full jets

     //draw pt spectrum
     if(!hPtJet[itype][0]) continue; 
     c[nCan] = new TCanvas(Form("c%d",nCan),Form("c%d: Pt %s jets",nCan,suffix[itype].Data()),600,450);
     SetCanvas((TCanvas*) c[nCan]);
     c[nCan]->SetLogy();

     frame[nCan] = gPad->DrawFrame(hPtJet[itype][0]->GetBinLowEdge(1),  
                                   1e-7,
                                   hPtJet[itype][0]->GetBinLowEdge(hPtJet[itype][0]->GetNbinsX()+1), 
                                   hPtJet[itype][0]->GetBinContent(hPtJet[itype][0]->GetMaximumBin())*2.);

     SetHist((TH1F*) frame[nCan],hPtJet[itype][0]->GetXaxis()->GetTitle(),hPtJet[itype][0]->GetYaxis()->GetTitle());

     if(drawComp) {
       hPtJet[itype][1]->SetLineColor(2);
       hPtJet[itype][1]->SetLineWidth(5);
       hPtJet[itype][1]->DrawCopy("same");
     }
     hPtJet[itype][0]->DrawCopy("same");
 
     leg = new TLegend(0.15,0.5,0.88,0.88);
     SetLeg(leg);
     
     TString txt =  Form("%s jets AKT R=%.1f",suffix[itype].Data(),jetR);
     if(run>0) txt += Form(" run:%d",run); 
     leg->AddEntry((TObject*) 0, txt.Data(),"");
     leg->AddEntry((TObject*) 0, Form("Trigger: %s",longtrigname.Data()),"");
     leg->AddEntry((TObject*) 0, Form("p_{T,trk}> %d MeV/c",TMath::Nint(minTrkPT*1000)),"");
     if(itype==1) leg->AddEntry((TObject*) 0, Form("E_{T}>%d MeV/c",TMath::Nint(minClusterET*1000)),"");
     if(drawComp){
       leg->AddEntry(hPtJet[itype][0], "Recent data", "l");
       leg->AddEntry(hPtJet[itype][1], "Previous data", "l");
     }
     leg->AddEntry((TObject*) 0, Form("#it{N}_{events} = %.0f",nEvents),"");
     leg->Draw();

     c[nCan]->SaveAs(Form("%s_Pt_AKT%02d_pT%04d_ET%04d_Run%d_Trigger%s_%s.%s",prefix.Data(),TMath::Nint(jetR*10),
			  TMath::Nint(minTrkPT*1000),TMath::Nint(minClusterET*1000), run, longtrigname.Data(), suffix[itype].Data(), suftype.Data()));

     nCan++;
  }//end of loop over charged and full jets
  //_________________
  //draw eta versus phi
  for(Int_t itype = 0; itype < kJetType; itype++){ //loop over charged and full jets
    if(!hEtaPhiJet[itype][0]) continue;
    c[nCan] = new TCanvas(Form("c%d",nCan),Form("c%d: eta-phi %s jets",nCan,suffix[itype].Data()),600,450);
    if(drawComp) c[nCan]->Divide(2,1);

    for(Int_t comp = 0; comp < 2; comp++) {
      if((comp==1)&&(!drawComp)) continue;
      c[nCan]->cd(comp+1);
      SetCanvas(gPad);
      gPad->SetRightMargin(0.14);
//      c[nCan]->SetRightMargin(0.15);
  
      frame[nCan] = gPad->DrawFrame(hEtaPhiJet[itype][comp]->GetXaxis()->GetBinLowEdge(1),  
                                     hEtaPhiJet[itype][comp]->GetYaxis()->GetBinLowEdge(1),
                                     hEtaPhiJet[itype][comp]->GetXaxis()->GetBinLowEdge(hEtaPhiJet[itype][comp]->GetNbinsX()),
                                     hEtaPhiJet[itype][comp]->GetYaxis()->GetBinLowEdge(hEtaPhiJet[itype][comp]->GetNbinsY()));
  
      SetHist((TH1F*) frame[nCan],hEtaPhiJet[itype][comp]->GetXaxis()->GetTitle(),hEtaPhiJet[itype][comp]->GetYaxis()->GetTitle());
  
      hEtaPhiJet[itype][comp]->DrawCopy("colz");
  
      leg = new TLegend(0.15,0.5,0.88,0.88);
      SetLeg(leg);
  
      TString txt =  Form("%s jets AKT R=%.1f",suffix[itype].Data(),jetR);
      if(run>0) txt += Form(" run:%d",run); 
      leg->AddEntry((TObject*) 0, txt.Data(),"");
      if(drawComp) leg->AddEntry((TObject*) 0, comp ? "Recent production" : "Previous production", "");
      leg->AddEntry((TObject*) 0, Form("Trigger: %s",longtrigname.Data()),"");
      leg->AddEntry((TObject*) 0, Form("p_{T,trk}> %d MeV/c",TMath::Nint(minTrkPT*1000)),"");
      if(itype==1) leg->AddEntry((TObject*) 0, Form("E_{T}>%d MeV/c",TMath::Nint(minClusterET*1000)),"");
      leg->AddEntry((TObject*) 0, Form("#it{N}_{events} = %.0f",nEvents),"");
      leg->AddEntry((TObject*) 0, Form("p_{T,corr}^{jet} > %.1f GeV/c", ptmin),"");
      leg->Draw();
    }

    c[nCan]->SaveAs(Form("%s_EtaPhi_AKT%02d_pT%04d_ET%04d_Run%d_Trigger%s_%s.%s",prefix.Data(),TMath::Nint(jetR*10),
			  TMath::Nint(minTrkPT*1000),TMath::Nint(minClusterET*1000), run, longtrigname.Data(), suffix[itype].Data(), suftype.Data()));

    nCan++;
  }//end of the loop over charged and full jets
  //_________________
  //draw rho versus centrality
  for(Int_t itype = 0; itype < kJetType; itype++){ //loop over charged and full jets
     if(!hRhoCent[itype][0]) continue;
     c[nCan] = new TCanvas(Form("c%d",nCan),Form("c%d: Rho-Cent",nCan),600,450);
     if(drawComp) c[nCan]->Divide(2,1);

     for(Int_t comp = 0; comp < 2; comp++) {
       if((comp==1)&&(!drawComp)) continue;
       c[nCan]->cd(comp+1);
       SetCanvas(gPad);
       gPad->SetRightMargin(0.15);
       c[nCan]->SetLogz();
  
       frame[nCan] = gPad->DrawFrame(hRhoCent[itype][comp]->GetXaxis()->GetBinLowEdge(1),
                                     hRhoCent[itype][comp]->GetYaxis()->GetBinLowEdge(1),
                                     hRhoCent[itype][comp]->GetXaxis()->GetBinLowEdge(hRhoCent[itype][comp]->GetNbinsX()),
                                     hRhoCent[itype][comp]->GetYaxis()->GetBinLowEdge(hRhoCent[itype][comp]->GetNbinsY()));
  
       SetHist((TH1F*) frame[nCan],hRhoCent[itype][comp]->GetXaxis()->GetTitle(),hRhoCent[itype][comp]->GetYaxis()->GetTitle());

       hRhoCent[itype][comp]->DrawCopy("colz");
  
       leg = new TLegend(-0.05,0.5,0.88,0.88);
       SetLeg(leg);
  
       TString txt =  Form("%s jets AKT R=%.1f",suffix[itype].Data(),jetR);
       if(run>0) txt += Form(" run:%d",run); 
  
       leg->AddEntry((TObject*) 0, txt.Data(),"");
       if(drawComp) leg->AddEntry((TObject*) 0, comp ? "Recent production" : "Previous production", "");
       leg->AddEntry((TObject*) 0, Form("Trigger: %s",longtrigname.Data()),"");
       leg->AddEntry((TObject*) 0, Form("#it{N}_{events} = %.0f",nEvents),"");
       leg->Draw();
    }

    c[nCan]->SaveAs(Form("%s_RhoCent_AKT%02d_pT%04d_ET%04d_Run%d_Trigger%s_%s.%s",prefix.Data(),TMath::Nint(jetR*10),
                         TMath::Nint(minTrkPT*1000),TMath::Nint(minClusterET*1000), run, longtrigname.Data(), suffix[itype].Data(), suftype.Data()));
    nCan++;
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

    if(hPtJet[itype][0])     hPtJet[itype][0]->Write(Form("%s%d_%s",prefix.Data(), itype, hPtJet[itype][0]->GetName()));
    if(drawComp) if(hPtJet[itype][1])     hPtJet[itype][1]->Write(Form("%s%d_%s",prefix.Data(), itype, hPtJet[itype][1]->GetName()));
    if(hEtaPhiJet[itype][0]) hEtaPhiJet[itype][0]->Write(Form("%s%d_%s",prefix.Data(), itype, hEtaPhiJet[itype][0]->GetName()));
    if(drawComp) if(hEtaPhiJet[itype][1]) hEtaPhiJet[itype][1]->Write(Form("%s%d_%s",prefix.Data(), itype, hEtaPhiJet[itype][1]->GetName()));
    if(hRhoCent[itype][0])   hRhoCent[itype][0]->Write(Form("%s%d_%s",prefix.Data(), itype, hRhoCent[itype][0]->GetName()));
    if(drawComp) if(hRhoCent[itype][1])   hRhoCent[itype][1]->Write(Form("%s%d_%s",prefix.Data(), itype, hRhoCent[itype][1]->GetName()));
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

void SetCanvas(TVirtualPad* c){
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
