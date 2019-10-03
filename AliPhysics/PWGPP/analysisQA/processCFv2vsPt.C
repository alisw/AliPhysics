/***************************************************************
 processCFv2vsPt.C
 Post Processing of Flow QA task in Analysis QA train 
 Clone copy of v2vsPt_0d.C from Ante

 Modification Done // sjena
 To save file and putting ani unique naming convention

***************************************************************/

// Name of the merged, large statistics file obtained with the merging macros:
TString mergedFileName = "AnalysisResults.root";

// Final output file name holding correct final results for large statistics sample:
TString outputFileName = "AnalysisResults.root";

void processCFv2vsPt(const char *filePath = "AnalysisResults.root",
		     TString suffix="eps", 
		     const char *outfile="CFv2vsPt_output.root") {
 // Load flow libraries:
 gSystem->Load("libPWGflowBase");

 // Let's go:
 TString mergedFileFullPathName(filePath);
 TFile *mergedFile = NULL;
 if(gSystem->AccessPathName(mergedFileFullPathName.Data(),kFileExists))
 {
  cout<<endl;
  cout<<" WARNING: Couldn't find a file: "<<mergedFileFullPathName.Data()<<endl;
  cout<<endl;
  exit(0);
 } else
   {
    // Create temporarily copy of <mergedFileName> if neccessary:
    if(!(mergedFileName == outputFileName))
    {
     TSystemFile *fileTemp = new TSystemFile(mergedFileFullPathName.Data(),".");
     fileTemp->Copy("mergedAnalysisResultsTemp.root");
     delete fileTemp;
    }
    // Access <mergedFileName>:
    mergedFile = TFile::Open(mergedFileFullPathName.Data(),"UPDATE");
  }

 // Accessing the file and updating it:  
 TList* mergedFileKeys = mergedFile->GetListOfKeys();
 for(Int_t i=0; i<mergedFileKeys->GetEntries(); i++)
 {
  TDirectory* directory = dynamic_cast<TDirectory*>(mergedFile->Get(mergedFileKeys->At(i)->GetName()));
  if (!directory) continue;
  TList* listTemp = directory->GetListOfKeys();
  if(!listTemp) continue;
  for (Int_t icent=0; icent<listTemp->GetEntries(); icent++)
  {
   TList* list = dynamic_cast<TList*>(directory->Get(listTemp->At(icent)->GetName()));
   if (!list) continue;
   // QC:
   if(TString(list->GetName()).Contains("QC"))
   {
    AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
    qc->GetOutputHistograms(list);
    Bool_t bApplyCorrectionForNUA = kFALSE; // apply correction for non-uniform acceptance
    qc->GetIntFlowFlags()->SetBinContent(3,(Int_t)bApplyCorrectionForNUA);
    qc->Finish();
    directory->Add(list,kTRUE);
    directory->Write(directory->GetName(),TObject::kSingleKey+TObject::kWriteDelete);
   } // if(TString(list->GetName()).Contains("QC"))
  } // for (Int_t icent=0; icent<listTemp->GetEntries(); icent++)
 } // for(Int_t i=0; i<mergedFileKeys->GetEntries(); i++)

 // Close the final output file:
 delete mergedFile;

 // Changing the final name, if neccessary:
 if(!(mergedFileName == outputFileName))
 {
  TSystemFile *outputFileFinal = new TSystemFile(mergedFileName.Data(),".");
  outputFileFinal->Rename(outputFileName.Data());
  delete outputFileFinal;
  TSystemFile *mergedFileFinal = new TSystemFile("mergedAnalysisResultsTemp.root",".");
  mergedFileFinal->Rename(mergedFileName.Data());
  delete mergedFileFinal;
 } // end of if(!(mergedFileName == outputFileName))

 // Finally, make some nice plots...:
 // Access "v2 vs. pT" histogram from the specified output file: 
 TFile *f = TFile::Open(filePath,"READ");
 TDirectoryFile *df = (TDirectoryFile*)f->FindObjectAny("outputQCanalysisTPCstandalone");
 if(!df){df = (TDirectoryFile*)f->FindObjectAny("outputQCanalysis");} // for toy MC
 if(!df){cout<<"df is NULL pointer!!!!"<<endl; exit(0);}
 TList *l = NULL;
 TList *listTemp = df->GetListOfKeys();
 if(listTemp && listTemp->GetEntries() == 1)
 {
  TString listName = listTemp->At(0)->GetName(); 
  df->GetObject(listName.Data(),l);
 } 
 TList *diffFlowList = dynamic_cast<TList*> l->FindObject("Differential Flow");
 TList *diffFlowResultsList = dynamic_cast<TList*> diffFlowList->FindObject("Results");
 TList *diffFlowDifFlowList = dynamic_cast<TList*> diffFlowResultsList->FindObject("Differential flow (POI, p_{T})");

 // v2{2} vs. pT:
 TH1D *v22 = dynamic_cast<TH1D*> diffFlowDifFlowList->FindObject("fDiffFlow, POI, p_{T}, v'{2}");
 if(!v22){cout<<"v22 is NULL"<<endl; exit(0);}
 v22->SetMarkerStyle(kFullCircle);
 v22->SetMarkerColor(kBlue);
 v22->SetLineColor(kBlue);

 // v2{4} vs. pT:
 TH1D *v24 = dynamic_cast<TH1D*> diffFlowDifFlowList->FindObject("fDiffFlow, POI, p_{T}, v'{4}");
 if(!v24){cout<<"v24 is NULL"<<endl; exit(0);}
 v24->SetMarkerStyle(kFullSquare);
 v24->SetMarkerColor(kRed);
 v24->SetLineColor(kRed);

 // Legend:
 TLegend *legend = new TLegend(0.15,0.67,0.37,0.87);
 legend->SetFillStyle(0); // white legend background
 legend->SetTextSize(0.04);
 legend->SetBorderSize(0.0);
 //legend->SetTextFont(22);
 legend->SetMargin(0.1);  
 //legend->AddEntry(inputValues_g,"theoretical/input values","p");
 legend->AddEntry(v22,"v_{2}{2}","p");
 legend->AddEntry(v24,"v_{2}{4}","p");
 
 // Style histogram:
 TH1D *hist = new TH1D("","",10,0.,10.);
 hist->SetStats(kFALSE);
 hist->SetTitle("");
 hist->GetXaxis()->SetTitle("p_{T}");
 Double_t maxFlow_y = v22->GetMaximum(0.75); // Hardwiring that max_value is <= 0.75 for aesthetic reasons
 if(v24->GetMaximum(0.75)>maxFlow_y){maxFlow_y = v24->GetMaximum(0.75);}
 hist->GetYaxis()->SetRangeUser(0.0,maxFlow_y);

 // Final drawing:
 // v2 vs. pT:
 TCanvas *c1 = new TCanvas("c1","v2 vs. pT"); 
 hist->Draw();
 legend->Draw("same");
 v22->Draw("psame");
 v24->Draw("psame");

 c1->SaveAs(Form("fig_cf_flow_fDiffFlow.%s",suffix.Data()));

 //===========================================================================================  

 // Bias from non-uniform azimuthal acceptance (NUA):
 TList *intFlowList = dynamic_cast<TList*> l->FindObject("Integrated Flow");
 intFlowList = dynamic_cast<TList*> intFlowList->FindObject("Results");
 TH1D *fIntFlowDetectorBias = dynamic_cast<TH1D*> intFlowList->FindObject("fIntFlowDetectorBias");
 fIntFlowDetectorBias->SetStats(kFALSE);
 fIntFlowDetectorBias->SetTitle("Quantifying bias from non-uniform acceptance");
 fIntFlowDetectorBias->SetFillColor(kBlue-10);
 fIntFlowDetectorBias->GetXaxis()->SetRangeUser(0.,2.);
 fIntFlowDetectorBias->GetYaxis()->SetRangeUser(0.8044,1.0544);
 fIntFlowDetectorBias->SetBinError(1,0.); // stat. errors need to be reimplemented
 fIntFlowDetectorBias->SetBinError(2,0.); // stat. errors need to be reimplemented

 // Final drawing:
 TCanvas *c2 = new TCanvas("c2","NUA"); 
 fIntFlowDetectorBias->Draw();

 c2->SaveAs(Form("fig_cf_flow_IntFlowDetectorBias.%s",suffix.Data()));

 //===========================================================================================  

 // Differential (vs. pT) cumulants:
 if(!diffFlowResultsList){cout<<"diffFlowResultsList is NULL"<<endl; exit(0);}
 TList *diffFlowDifCumulantsList = dynamic_cast<TList*> diffFlowResultsList->FindObject("Differential Q-cumulants (POI, p_{T})");
 if(!diffFlowDifCumulantsList){cout<<"diffFlowDifCumulantsList is NULL"<<endl; exit(0);}

 // d2{2} vs. pT:
 TH1D *d22 = dynamic_cast<TH1D*> diffFlowDifCumulantsList->FindObject("fDiffFlowCumulants, POI, p_{T}, QC{2'}");
 if(!d22){cout<<"d22 is NULL"<<endl; exit(0);}
 d22->SetMarkerStyle(kFullCircle);
 d22->SetMarkerColor(kBlue);
 d22->SetLineColor(kBlue);
 d22->Scale(10); // Rescale this histogram by 10^2

 // d2{4} vs. pT:
 TH1D *d24 = dynamic_cast<TH1D*> diffFlowDifCumulantsList->FindObject("fDiffFlowCumulants, POI, p_{T}, QC{4'}");
 if(!d24){cout<<"d24 is NULL"<<endl; exit(0);}
 d24->SetMarkerStyle(kFullSquare);
 d24->SetMarkerColor(kRed);
 d24->SetLineColor(kRed);
 d24->Scale(1000); // Rescale this histogram by 10^4

 // Legend:
 TLegend *legend_c = new TLegend(0.15,0.67,0.37,0.87);
 legend_c->SetFillStyle(0); // white legend background
 legend_c->SetTextSize(0.04);
 legend_c->SetBorderSize(0.0);
 //legend_c->SetTextFont(22);
 legend_c->SetMargin(0.1);  
 legend_c->AddEntry(d22,"d_{2}{2} #times 10","p");
 legend_c->AddEntry(d24,"d_{2}{4} #times 10^{3}","p");
 
 // Style histogram:
 TH1D *hist_c = new TH1D("","",10,0.,10.);
 hist_c->SetStats(kFALSE);
 hist_c->SetTitle("");
 hist_c->GetXaxis()->SetTitle("p_{T}");
 Double_t maxCumulant_y = d22->GetMaximum(0.75); // Hardwiring that max_value is <= 0.75 for aesthetic reasons
 if(d24->GetMaximum(0.75)>maxCumulant_y){maxCumulant_y = d24->GetMaximum(0.75);}
 Double_t minCumulant_y = d22->GetMinimum(-0.75); // Hardwiring that max_value is <= 0.75 for aesthetic reasons
 if(d24->GetMinimum(-0.75)<minCumulant_y){minCumulant_y = d24->GetMinimum(-0.75);}
 hist_c->GetYaxis()->SetRangeUser(minCumulant_y,maxCumulant_y);

 // Final drawing:
 // d2 vs. pT:
 TCanvas *c3 = new TCanvas("c3","differential Q-cumulants vs. pT"); 
 hist_c->Draw();
 legend_c->Draw("same");
 d22->Draw("psame");
 d24->Draw("psame");

 c3->SaveAs(Form("fig_cf_flow_fDiffFlowCumulants.%s",suffix.Data()));

 //f->Close();

 cout<<endl;
  
 //Added by sjena
 TFile *fout = TFile::Open(outfile,"UPDATE");
 fout->ls();
 
 TDirectoryFile *cdd = NULL;
 cdd = (TDirectoryFile*)fout->Get("CF");
 if(!cdd) {
    Printf("Warning: CF <dir> doesn't exist, creating a new one");
    cdd = (TDirectoryFile*)fout->mkdir("CF");
 }
 cdd->cd();
 cdd->ls();
 
 
 v22->SetName(Form("fig_cf_flow_22", v22->GetName()));
 // v22->SetName(Form("fig_cf_flow_fDiffFlow22"));
 v22->Write();
 
 v24->SetName(Form("fig_cf_flow_24", v24->GetName()));
 // v24->SetName(Form("fig_cf_flow_fDiffFlow24"));
 v24->Write();
 
 fIntFlowDetectorBias->SetName(Form("fig_cf_flow_IntFlowDetectorBias", fIntFlowDetectorBias->GetName()));
 fIntFlowDetectorBias->Write();
 
 fout->cd();
 fout->Close();
 

 
} // void void v2vsPt_0d()

