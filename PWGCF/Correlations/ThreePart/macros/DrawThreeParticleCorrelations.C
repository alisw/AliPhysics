// #include "AliAnalysisTaskCorrelation3p.h"
#include "TFile.h"
#include "TObject.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TNamed.h"
// #include <ZTrees.h>
//Run with aliroot  -b -l ../run-single-task.C'("compile", "AliAnalysisTaskCorrelation3p.cxx")' DrawThreeParticleCorrelations.C
void Correct(const char* histname, TDirectory* SameDir, TDirectory* METADir, TDirectory* METriggerdir, TFile* resultsfile)
{
  TDirectory* ResultsSame=resultsfile->GetDirectory("Same");
  TDirectory * ResultsMinusMETA=resultsfile->GetDirectory("ResultMinusMETA");
  TDirectory * ResultsMinusMETrigger = resultsfile->GetDirectory("ResultMinusMETrigger");
  TDirectory * ResultsMinusBoth = resultsfile->GetDirectory("ResultsMinusBoth");
  TDirectory * Compare = resultsfile->GetDirectory("Compare");

  if(!ResultsSame) ResultsSame=resultsfile->mkdir("Same");
  if(!ResultsMinusMETA) ResultsMinusMETA=resultsfile->mkdir("ResultMinusMETA");
  if(!ResultsMinusMETrigger) ResultsMinusMETrigger=resultsfile->mkdir("ResultMinusMETrigger");
  if(!ResultsMinusBoth) ResultsMinusBoth=resultsfile->mkdir("ResultsMinusBoth");
  if(!Compare) Compare = resultsfile->mkdir("Compare");
  
  TH2D * SameHist = dynamic_cast<TH2D*>(SameDir->Get(histname)->Clone(Form("%s%s",histname,"Same")));
  TH2D * SameMinusMETA = dynamic_cast<TH2D*>(SameHist->Clone(Form("%s%s",histname,"SameMinusMETA")));
  SameMinusMETA->SetTitle(Form("%s%s",SameMinusMETA->GetTitle()," with the corrected Mixed event with trigger and first associated from same event substracted."));
  TH2D * SameMinusMETrigger = dynamic_cast<TH2D*>(SameHist->Clone(Form("%s%s",histname,"SameMinusMETrigger")));
  SameMinusMETrigger->SetTitle(Form("%s%s",SameMinusMETrigger->GetTitle()," with the corrected Mixed event with both associated from same event substracted."));
  TH2D * SameMinusBoth = dynamic_cast<TH2D*>(SameHist->Clone(Form("%s%s",histname,"SameMinusBoth")));
  SameMinusBoth->SetTitle(Form("%s%s",SameMinusBoth->GetTitle()," with both mixed event correlations substracted."));
  TH2D * METAhist = dynamic_cast<TH2D*>(METADir->Get(histname)->Clone(Form("%s%s",histname,"META")));
  TH2D * METriggerhist = dynamic_cast<TH2D*>(METriggerdir->Get(histname)->Clone(Form("%s%s",histname,"METrigger")));

  SameMinusMETA->Add(METAhist,-1.0);
  SameMinusMETrigger->Add(METriggerhist,-1.0);
  SameMinusBoth->Add(METAhist,-1.0);
  SameMinusBoth->Add(METriggerhist,-1.0);
  TCanvas * c1 = new TCanvas(Form("%s%s",histname,"Canvas"));
  c1->Divide(2,2);
  c1->cd(1);
  SameHist->Draw("surf2");
  ResultsSame->cd();
  SameHist->Write(histname);
  c1->cd(2);
  SameMinusMETA->Draw("surf2");
  ResultsMinusMETA->cd(); 
  SameMinusMETA->Write(histname);
  c1->cd(3);
  SameMinusMETrigger->Draw("surf2");
  ResultsMinusMETrigger->cd();
  SameMinusMETrigger->Write(histname);
  c1->cd(4);
  SameMinusBoth->Draw("surf2"); 
  ResultsMinusBoth->cd();
  SameMinusBoth->Write(histname);
  Compare->cd();
  c1->Write(histname);

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1);
  SameMinusBoth->Draw("surf2"); 
  c1->cd(2);
  SameMinusBoth->Draw("surf3"); 
  c1->cd(3);
  Int_t binpihn = SameMinusBoth->GetYaxis()->FindBin(0+0.2);
  Int_t binpiln = SameMinusBoth->GetYaxis()->FindBin(0-0.2);
  TH1D* projY = SameMinusBoth->ProjectionX(Form("%s%s",SameMinusBoth->GetName(),"_nearside"),binpiln,binpihn);
  projY->SetStats(0);
  projY->GetYaxis()->SetRangeUser(0., 1.1*projY->GetBinContent(projY->GetMaximumBin()));
  projY->SetTitle("Integral over the near side peak with #Delta#eta_{12} = 0#pm 0.2");
  projY->GetYaxis()->SetTitle(SameMinusBoth->GetZaxis()->GetTitle());
  projY->SetTitleSize(0.04,"xy");
  projY->SetTitleOffset(1.05,"xy");
  projY->Draw("E");
  c1->cd(4);
  Int_t binpih = SameMinusBoth->GetYaxis()->FindBin(TMath::Pi()+0.2);
  Int_t binpil = SameMinusBoth->GetYaxis()->FindBin(TMath::Pi()-0.2);
  TH1D* projX = SameMinusBoth->ProjectionX(Form("%s%s",SameMinusBoth->GetName(),"_px"),binpil,binpih);
  projX->SetStats(0);
  projX->SetTitle("Integral over a slice of the away side peak around with #Delta#eta_{12} = #pi#pm 0.2");
  projX->GetYaxis()->SetTitle(SameMinusBoth->GetZaxis()->GetTitle());
  projX->SetTitleSize(0.04,"xy");
  projX->SetTitleOffset(1.05,"xy");
  projX->Draw("E");

  ResultsMinusBoth->cd();
  c1->Write(Form("%s%s",histname,"Canvas"));
  delete c1;
}

void Correct(const char* histname1,const char* histname2,const char* histname3,const char* histname4, TDirectory* SameDir, TDirectory* METADir, TDirectory* METriggerdir, TFile* resultsfile, bool ispp=false, bool isPbPb=false, Double_t * METAscale=NULL, Double_t * METriggerscale =NULL)
{
  TDirectory* ResultsSame=resultsfile->GetDirectory("Same");
  TDirectory * ResultsMinusMETA=resultsfile->GetDirectory("ResultMinusMETA");
  TDirectory * ResultsMinusMETrigger = resultsfile->GetDirectory("ResultMinusMETrigger");
  TDirectory * ResultsMinusBoth = resultsfile->GetDirectory("ResultsMinusBoth");
  TDirectory * Compare = resultsfile->GetDirectory("Compare");
  if(!ResultsSame) ResultsSame=resultsfile->mkdir("Same");
  if(!ResultsMinusMETA) ResultsMinusMETA=resultsfile->mkdir("ResultMinusMETA");
  if(!ResultsMinusMETrigger) ResultsMinusMETrigger=resultsfile->mkdir("ResultMinusMETrigger");
  if(!ResultsMinusBoth) ResultsMinusBoth=resultsfile->mkdir("ResultsMinusBoth");
  if(!Compare) Compare = resultsfile->mkdir("Compare");
  
  TH2D * SameHist1 = dynamic_cast<TH2D*>(SameDir->Get(histname1)->Clone(Form("%s%s",histname1,"Same")));
  TH2D * SameMinusMETA1 = dynamic_cast<TH2D*>(SameHist1->Clone(Form("%s%s",histname1,"SameMinusMETA")));
  SameMinusMETA1->SetTitle(Form("%s%s",SameMinusMETA1->GetTitle()," with the corrected Mixed event with trigger and first associated from same event substracted."));
  TH2D * SameMinusMETrigger1 = dynamic_cast<TH2D*>(SameHist1->Clone(Form("%s%s",histname1,"SameMinusMETrigger")));
  SameMinusMETrigger1->SetTitle(Form("%s%s",SameMinusMETrigger1->GetTitle()," with the corrected Mixed event with both associated from same event substracted."));
  TH2D * SameMinusBoth1 = dynamic_cast<TH2D*>(SameHist1->Clone(Form("%s%s",histname1,"SameMinusBoth")));
  SameMinusBoth1->SetTitle(Form("%s%s",SameMinusBoth1->GetTitle()," with both mixed event correlations substracted."));
  TH2D * METAhist1 = dynamic_cast<TH2D*>(METADir->Get(histname1)->Clone(Form("%s%s",histname1,"META")));
  TH2D * METriggerhist1 = dynamic_cast<TH2D*>(METriggerdir->Get(histname1)->Clone(Form("%s%s",histname1,"METrigger")));
  if(ispp){
    METAscale = SameHist1->Integral(SameHist1->GetXaxis()->FindBin(0),SameHist1->GetXaxis()->FindBin(0),1,SameHist1->GetYaxis()->FindBin(-1));
    METAscale +=SameHist1->Integral(1,SameHist1->GetXaxis()->FindBin(-1),SameHist1->GetYaxis()->FindBin(0),SameHist1->GetYaxis()->FindBin(0));
    Double_t Divideby = METAhist1->Integral(METAhist1->GetXaxis()->FindBin(0),METAhist1->GetXaxis()->FindBin(0),1,METAhist1->GetYaxis()->FindBin(-1));
    Divideby +=METAhist1->Integral(1,METAhist1->GetXaxis()->FindBin(-1),METAhist1->GetYaxis()->FindBin(0),METAhist1->GetYaxis()->FindBin(0));
    if(Divideby >0) METAscale = METAscale/Divideby;
    else METAscale = 0.0;
    METAhist1->Scale(METAscale);
    METriggerscale = SameHist1->Integral(SameHist1->GetXaxis()->FindBin(1.0),SameHist1->GetXaxis()->FindBin(2.0),SameHist1->GetYaxis()->FindBin(1.0),SameHist1->GetYaxis()->FindBin(2.0));
    Double_t Divideby2 = METriggerhist1->Integral(METriggerhist1->GetXaxis()->FindBin(1.0),METriggerhist1->GetXaxis()->FindBin(2.0),METriggerhist1->GetYaxis()->FindBin(1.0),METriggerhist1->GetYaxis()->FindBin(2.0));
    if(Divideby2 >0) METriggerscale = METriggerscale/Divideby2;
    else METriggerscale = 0.0;
    METriggerhist1->Scale(METriggerscale);
  }
  
  
  TH2D * SameHist2 = dynamic_cast<TH2D*>(SameDir->Get(histname2)->Clone(Form("%s%s",histname2,"Same")));
  TH2D * SameMinusMETA2 = dynamic_cast<TH2D*>(SameHist2->Clone(Form("%s%s",histname2,"SameMinusMETA")));
  SameMinusMETA2->SetTitle(Form("%s%s",SameMinusMETA2->GetTitle()," with the corrected Mixed event with trigger and first associated from same event substracted."));
  TH2D * SameMinusMETrigger2 = dynamic_cast<TH2D*>(SameHist2->Clone(Form("%s%s",histname2,"SameMinusMETrigger")));
  SameMinusMETrigger2->SetTitle(Form("%s%s",SameMinusMETrigger2->GetTitle()," with the corrected Mixed event with both associated from same event substracted."));
  TH2D * SameMinusBoth2 = dynamic_cast<TH2D*>(SameHist2->Clone(Form("%s%s",histname2,"SameMinusBoth")));
  SameMinusBoth2->SetTitle(Form("%s%s",SameMinusBoth2->GetTitle()," with both mixed event correlations substracted."));
  TH2D * METAhist2 = dynamic_cast<TH2D*>(METADir->Get(histname2)->Clone(Form("%s%s",histname2,"META")));
  TH2D * METriggerhist2 = dynamic_cast<TH2D*>(METriggerdir->Get(histname2)->Clone(Form("%s%s",histname2,"METrigger")));
  
  TH2D * SameHist3 = dynamic_cast<TH2D*>(SameDir->Get(histname3)->Clone(Form("%s%s",histname3,"Same")));
  TH2D * SameMinusMETA3 = dynamic_cast<TH2D*>(SameHist3->Clone(Form("%s%s",histname3,"SameMinusMETA")));
  SameMinusMETA3->SetTitle(Form("%s%s",SameMinusMETA3->GetTitle()," with the corrected Mixed event with trigger and first associated from same event substracted."));
  TH2D * SameMinusMETrigger3 = dynamic_cast<TH2D*>(SameHist3->Clone(Form("%s%s",histname3,"SameMinusMETrigger")));
  SameMinusMETrigger3->SetTitle(Form("%s%s",SameMinusMETrigger3->GetTitle()," with the corrected Mixed event with both associated from same event substracted."));
  TH2D * SameMinusBoth3 = dynamic_cast<TH2D*>(SameHist3->Clone(Form("%s%s",histname3,"SameMinusBoth")));
  SameMinusBoth3->SetTitle(Form("%s%s",SameMinusBoth3->GetTitle()," with both mixed event correlations substracted."));
  TH2D * METAhist3 = dynamic_cast<TH2D*>(METADir->Get(histname3)->Clone(Form("%s%s",histname3,"META")));
  TH2D * METriggerhist3 = dynamic_cast<TH2D*>(METriggerdir->Get(histname3)->Clone(Form("%s%s",histname3,"METrigger")));
  
  TH2D * SameHist4 = dynamic_cast<TH2D*>(SameDir->Get(histname4)->Clone(Form("%s%s",histname4,"Same")));
  TH2D * SameMinusMETA4 = dynamic_cast<TH2D*>(SameHist4->Clone(Form("%s%s",histname4,"SameMinusMETA")));
  SameMinusMETA4->SetTitle(Form("%s%s",SameMinusMETA4->GetTitle()," with the corrected Mixed event with trigger and first associated from same event substracted."));
  TH2D * SameMinusMETrigger4 = dynamic_cast<TH2D*>(SameHist4->Clone(Form("%s%s",histname4,"SameMinusMETrigger")));
  SameMinusMETrigger4->SetTitle(Form("%s%s",SameMinusMETrigger4->GetTitle()," with the corrected Mixed event with both associated from same event substracted."));
  TH2D * SameMinusBoth4 = dynamic_cast<TH2D*>(SameHist4->Clone(Form("%s%s",histname4,"SameMinusBoth")));
  SameMinusBoth4->SetTitle(Form("%s%s",SameMinusBoth4->GetTitle()," with both mixed event correlations substracted."));
  TH2D * METAhist4 = dynamic_cast<TH2D*>(METADir->Get(histname4)->Clone(Form("%s%s",histname4,"META")));
  TH2D * METriggerhist4 = dynamic_cast<TH2D*>(METriggerdir->Get(histname4)->Clone(Form("%s%s",histname4,"METrigger")));
  
  SameMinusMETA1->Add(METAhist1,-1.0);
  SameMinusMETrigger1->Add(METriggerhist1,-1.0);
  SameMinusBoth1->Add(METAhist1,-1.0);
  SameMinusBoth1->Add(METriggerhist1,-1.0);
  TCanvas * c1 = new TCanvas("shosw");
  c1->Divide(2,2);
  c1->cd(1);
  SameHist1->Draw("surf2");
  ResultsSame->cd();
  SameHist1->Write(histname1);
  c1->cd(2);
  SameMinusMETA1->Draw("surf2");
  ResultsMinusMETA->cd(); 
  SameMinusMETA1->Write(histname1);
  c1->cd(3);
  SameMinusMETrigger1->Draw("surf2");
  ResultsMinusMETrigger->cd();
  SameMinusMETrigger1->Write(histname1);
  c1->cd(4);
  SameMinusBoth1->Draw("surf2"); 
  ResultsMinusBoth->cd();
  SameMinusBoth1->Write(histname1);
  Compare->cd();
  c1->Write(histname1);
  
  c1->Clear();
  SameMinusMETA2->Add(METAhist2,-1.0);
  SameMinusMETrigger2->Add(METriggerhist2,-1.0);
  SameMinusBoth2->Add(METAhist2,-1.0);
  SameMinusBoth2->Add(METriggerhist2,-1.0);
  c1->Divide(2,2);
  c1->cd(1);
  SameHist2->Draw("surf2");
  ResultsSame->cd();
  SameHist2->Write(histname2);
  c1->cd(2);
  SameMinusMETA2->Draw("surf2");
  ResultsMinusMETA->cd(); 
  SameMinusMETA2->Write(histname2);
  c1->cd(3);
  SameMinusMETrigger2->Draw("surf2");
  ResultsMinusMETrigger->cd();
  SameMinusMETrigger2->Write(histname2);
  c1->cd(4);
  SameMinusBoth2->Draw("surf2"); 
  ResultsMinusBoth->cd();
  SameMinusBoth2->Write(histname2);
  Compare->cd();
  c1->Write(histname2);
  
  c1->Clear();
  SameMinusMETA3->Add(METAhist3,-1.0);
  SameMinusMETrigger3->Add(METriggerhist3,-1.0);
  SameMinusBoth3->Add(METAhist3,-1.0);
  SameMinusBoth3->Add(METriggerhist3,-1.0);
  c1->Divide(2,2);
  c1->cd(1);
  SameHist3->Draw("surf2");
  ResultsSame->cd();
  SameHist3->Write(histname3);
  c1->cd(2);
  SameMinusMETA3->Draw("surf2");
  ResultsMinusMETA->cd(); 
  SameMinusMETA3->Write(histname3);
  c1->cd(3);
  SameMinusMETrigger3->Draw("surf2");
  ResultsMinusMETrigger->cd();
  SameMinusMETrigger3->Write(histname3);
  c1->cd(4);
  SameMinusBoth3->Draw("surf2"); 
  ResultsMinusBoth->cd();
  SameMinusBoth3->Write(histname3);
  Compare->cd();
  c1->Write(histname3);
  
  c1->Clear();
  SameMinusMETA4->Add(METAhist4,-1.0);
  SameMinusMETrigger4->Add(METriggerhist4,-1.0);
  SameMinusBoth4->Add(METAhist4,-1.0);
  SameMinusBoth4->Add(METriggerhist4,-1.0);
  c1->Divide(2,2);
  c1->cd(1);
  SameHist4->Draw("surf2");
  ResultsSame->cd();
  SameHist4->Write(histname4);
  c1->cd(2);
  SameMinusMETA4->Draw("surf2");
  ResultsMinusMETA->cd(); 
  SameMinusMETA4->Write(histname4);
  c1->cd(4);
  SameMinusMETrigger4->Draw("surf2");
  ResultsMinusMETrigger->cd();
  SameMinusMETrigger4->Write(histname4);
  c1->cd(4);
  SameMinusBoth4->Draw("surf2"); 
  ResultsMinusBoth->cd();
  SameMinusBoth4->Write(histname4);
  Compare->cd();
  c1->Write(histname4);
  
  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1);
  SameMinusBoth1->Draw("surf2"); 
  c1->cd(2);
  SameMinusBoth2->Draw("surf2"); 
  c1->cd(3);
  SameMinusBoth3->Draw("surf2"); 
  c1->cd(4);
  SameMinusBoth4->Draw("surf2"); 
  
  ResultsMinusBoth->cd();
  c1->Write(Form("%s%s",histname1,"Canvas"));
  delete c1;
}

void DrawThreeParticleCorrelations(const char* options="")
{
  bool bdrawbins=true;
  bool bdrawbinstats=false;
  bool makefile =false;
  bool makefileall =false;
  bool makefilemixed =false;
  bool makefiletriggermixed =false;
  bool makefilegen =false;
  bool removeextracor =false;
  bool ispp = false;
  bool isPbPb = false;
  TString underdirectory("");
  TString delimiter(" ");
  TStringToken token(options, delimiter);
  while (token.NextToken()) {
    const char* key=0;
    TString argument=token;
    if (argument.CompareTo("-h")==0 ||argument.CompareTo("--help")==0) {
      AliInfo(Form("Draw options:"
		    "\n\t  mergefirst   - merges the Multiplicity bins and Vz bins before correcting with mixed events (default merges after mixing)."
		    "\n\t  drawbins     - each multiplicity-Vz bin is drawn alone."
		    "\n\t binstats      - draw bin statistics."
		    "\n\t mixed		- draw mixed."
		    "\n\t triggermixed  - draw trigger mixed."
		    ));
      return;
      }
    key="binstats";
    if (argument.CompareTo(key)==0) {
      bdrawbins=false;
      bdrawbinstats=true;
      underdirectory.Clear();
      underdirectory.Append("/binstats");
      continue;
    }
    key="makefile";
    if (argument.CompareTo(key)==0) {
      bdrawbins=false;
      makefile=true;
      continue;
    }    
    key="makefilemixedall";
    if (argument.CompareTo(key)==0) {
      bdrawbins=false;
      makefileall=true;
      continue;
    }    
    key="makefilemixed";
    if (argument.CompareTo(key)==0) {
      bdrawbins=false;
      makefilemixed=true;
      continue;
    }    
    key="makefiletriggermixed";
    if (argument.CompareTo(key)==0) {
      bdrawbins=false;
      makefiletriggermixed=true;
      continue;
    }    
    key="makefilegen";
    if (argument.CompareTo(key)==0) {
      bdrawbins=false;
      makefilegen=true;
      continue;
    }    
    key="correct";
    if(argument.CompareTo(key)==0){
      bdrawbins=false;
      removeextracor=true;
      continue;
    }
    key="pp";
    if(argument.CompareTo(key)==0){
      ispp=true;
      continue;
    }
    key="PbPb";
    if(argument.CompareTo(key)==0){
      isPbPb=true;
      continue;
    }
  }
  if(removeextracor){
    TFile * rfile = TFile::Open("results.root","READ");
    TFile * corrected = new TFile("corrected.root","RECREATE");
    if(!rfile){cout << "Run with makefile, makefilemixed, makefiletriggermixed or makefilegen first." <<endl;return;}
    TDirectory * sameventdiv = rfile->GetDirectory("divided");
    TDirectory * METAdiv = rfile->GetDirectory("META/divided");
    TDirectory * METriggerdiv = rfile->GetDirectory("METrigger/divided");
    Double_t * METAscale, METriggerscale;
    if(!ispp&&!isPbPb){
      Correct("DPhi_1_DPHI_2_from3D","DPhi_1_DPHI_2_near_from3D","DPhi_1_DPHI_2_mid_from3D","DPhi_1_DPHI_2_far_from3D",sameventdiv,METAdiv,METriggerdiv,corrected);
      Correct("DPhi_1_DEta_12_from3D",sameventdiv,METAdiv,METriggerdiv,corrected);
      Correct("DPhi_1_DEta_12_SameSide_from3D",sameventdiv,METAdiv,METriggerdiv,corrected);
    }
    else{
      Correct("DPhi_1_DPHI_2_from3D","DPhi_1_DPHI_2_near_from3D","DPhi_1_DPHI_2_mid_from3D","DPhi_1_DPHI_2_far_from3D",sameventdiv,METAdiv,METriggerdiv,corrected,ispp,isPbPb,METAscale,METriggerscale);
      Correct("DPhi_1_DEta_12_from3D",sameventdiv,METAdiv,METriggerdiv,corrected);
      Correct("DPhi_1_DEta_12_SameSide_from3D",sameventdiv,METAdiv,METriggerdiv,corrected);      
      cout << METAscale<<endl;
    }

	
//     TCanvas * c1 = new TCanvas("show");
//     c1->Divide(2,2);
//     sameventdiv->cd();
// //     sameventdiv->ls();
//     TH2D * DPHIDPHIs	 = dynamic_cast<TH2D*>(sameventdiv->Get("DPhi_1_DPHI_2_from2D")->Clone("DPhi_1_DPHI_2_from2Dsame"));
//     TH2D * DPHIDPHIscor1 = dynamic_cast<TH2D*>(DPHIDPHIs->Clone("DPhi_1_DPHI_2_from2Dcor1"));
//     TH2D * DPHIDPHIscor2 = dynamic_cast<TH2D*>(DPHIDPHIs->Clone("DPhi_1_DPHI_2_from2Dcor2"));
//     TH2D * DPHIDPHIscor3 = dynamic_cast<TH2D*>(DPHIDPHIs->Clone("DPhi_1_DPHI_2_from2Dcor3"));
// 
//     c1->cd(1);
//     DPHIDPHIs->Draw("colz");
//     METAdiv->cd();
//     TH2D * DPHIDPHIMETA;
//     METAdiv->GetObject("DPhi_1_DPHI_2_from2D",DPHIDPHIMETA);
//     METriggerdiv->cd();
//     TH2D * DPHIDPHIMETrigger;
//     METriggerdiv->GetObject("DPhi_1_DPHI_2_from2D",DPHIDPHIMETrigger);
//     
//     DPHIDPHIscor1->Add(DPHIDPHIMETA,-1.0);
//     DPHIDPHIscor3->Add(DPHIDPHIMETA,-1.0);
//     c1->cd(2);
//     DPHIDPHIscor1->Draw("colz");
//     DPHIDPHIscor2->Add(DPHIDPHIMETrigger,-1.0);
//     DPHIDPHIscor3->Add(DPHIDPHIMETrigger,-1.0);
//     c1->cd(3);
//     DPHIDPHIscor2->Draw("colz");
//     c1->cd(4);
//     DPHIDPHIscor3->Draw("colz");
    rfile->Close();corrected->Close();
    return;
  }
  
  TFile * ffile = TFile::Open("AnalysisResults.root", "READ");
  if(!ffile) return;  
  if(makefilegen){
    TDirectory * dir = ffile->GetDirectory("ThreePartGenerator");
    TList* fList;
    dir->GetObject("ThreePartGeneratorCoutput1",fList);
    AliCorrelation3p* generated;
    if(!(makefilemixed||makefiletriggermixed)) {generated = dynamic_cast<AliCorrelation3p*>fList->FindObject("tracktrigger_correlation_4_8");generated->MakeResultsFile("");}
    if(makefilemixed){generated = dynamic_cast<AliCorrelation3p*>fList->FindObject("tracktrigger_correlation_4_8META");generated->MakeResultsFile("META");}
    if(makefiletriggermixed)  {generated = dynamic_cast<AliCorrelation3p*>fList->FindObject("tracktrigger_correlation_4_8METrigger");generated->MakeResultsFile("METrigger");}

//     generated->MakeResultsFile("");
    return;
  }
  
  TStringToken triggers("Tracks pi0s"," ");
  while(triggers.NextToken()) {
    for (int i =0;i<2;i++){
      TString Directory(Form("ThreePart%s",triggers.Data()));
      if(i == 1) Directory.Append("PbPb");
      TDirectory * dir = ffile->GetDirectory(Directory.Data());
      if(!dir) continue;//directory does not exist.
      TList* fList;
      dir->GetObject(Form("%s%s",Directory.Data(),"Coutput1"),fList);
      if(!fList) continue;//list does not exist.
      if(bdrawbinstats){
	TH1D* nass = (TH1D*)fList->FindObject("NAssociated");
	TH1D* ntrig = (TH1D*)fList->FindObject("Ntriggers");
	TH1D* nasstrig =(TH1D*)fList->FindObject("NAssociatedETriggered");
	TCanvas* c = new TCanvas("numberoftriggers");
	c->Divide(2,2);
	c->cd(1);
	gPad->SetLogy();
// 	if(i==0)ntrig->GetXaxis()->SetRange(0,15);
	ntrig->GetXaxis()->SetTitle("N_{triggers}");
	ntrig->SetStats(0);
	ntrig->Draw("E");
	double eventwithoutatrigger = ntrig->GetBinContent(1);
	int eventwithoutatriggerint = ntrig->GetBinContent(1);
	double eventswithatrigger = ntrig->Integral(2,ntrig->GetNbinsX());
	int eventswithatriggerint = ntrig->Integral(2,ntrig->GetNbinsX());
	double ratio = eventswithatrigger/(eventwithoutatrigger+eventwithoutatrigger)*100;
	c->cd(2);
	gPad->SetLogy();
// 	if(i ==0 ) nass->GetXaxis()->SetRange(0,25);
// 	if(i ==1 ) nass->GetXaxis()->SetRange(0,100);
	nass->GetXaxis()->SetTitle("N_{associated}");
	nass->SetStats(0);
	nass->Draw("E");
	c->cd(3);
	gPad->SetLogy();
// 	if(i==0)nasstrig->GetXaxis()->SetRange(0,25);
// 	if(i==1)nasstrig->GetXaxis()->SetRange(0,100);
	nasstrig->GetXaxis()->SetTitle("N_{associated}");
	nasstrig->SetTitle("Number of Associated in Events that contain at least one trigger.");
	nasstrig->SetStats(0);
	nasstrig->Draw("E");
	double eventswithatriggerandtwoassociated = nasstrig->Integral(3,nasstrig->GetNbinsX());
	c->cd(4);
	TPaveText* text = new TPaveText(.05,.1,.95,.8,"nb");
	text->SetLineColor(0);
	text->SetFillStyle(0);
	text->AddText(Form("#bullet %4.1f %% of all events contain at least one trigger.",ratio));
	text->AddText(Form("#bullet %.1e events with at least one trigger.",eventswithatrigger));
	text->AddText(Form("#bullet %.1e events with\n at least one trigger and at least two associated.",eventswithatriggerandtwoassociated));
	text->Draw();
	c->Print("Triggersperevent.eps");
	TH3D* eventsbeforeselection = (TH3D*)fList->FindObject("Eventbeforeselection");
	eventsbeforeselection->GetYaxis()->SetTitle("Multiplicity");
	eventsbeforeselection->GetXaxis()->SetTitle("ZVertex");
	eventsbeforeselection->GetZaxis()->SetTitle("Centrality");
	eventsbeforeselection->SetTitleSize(0.045,"xyz");
	eventsbeforeselection->SetTitleOffset(1.3,"xyz");
	TH3D* eventsafterselection  = (TH3D*)fList->FindObject("Eventafterselection");
	eventsafterselection->GetYaxis()->SetTitle("Multiplicity");
	eventsafterselection->GetXaxis()->SetTitle("ZVertex");
	eventsafterselection->GetZaxis()->SetTitle("Centrality");
	eventsafterselection->SetTitleSize(0.045,"xyz");
	eventsafterselection->SetTitleOffset(1.3,"xyz");
	TH2D* multvertexbeforeselection = (TH2D*) eventsbeforeselection->Project3D("yx");
	multvertexbeforeselection->SetTitle("multiplicity vs z vertex before event selection.");
	TH2D* centvertexbeforeselection = (TH2D*) eventsbeforeselection->Project3D("zx");
	centvertexbeforeselection->SetTitle("centrality vs z vertex before event selection.");
	TH1D* multbeforeselection = (TH1D*) eventsbeforeselection->Project3D("x");
	multbeforeselection->SetTitle("multiplicity before event selection");
	TH1D* vertexbeforeselection = (TH1D*)  eventsbeforeselection->Project3D("y");
	vertexbeforeselection->SetTitle("vertex before event selection");
	TH1D* centralitybeforeselection = (TH1D*)  eventsbeforeselection->Project3D("z");
	centralitybeforeselection->SetTitle("centrality before event selection");
	TH2D* multvertexafterselection = (TH2D*) eventsafterselection->Project3D("yx");
	multvertexafterselection->SetTitle("multiplicity vs z vertex after event selection.");
	TH2D* centvertexafterselection = (TH2D*) eventsafterselection->Project3D("zx");
	centvertexafterselection->SetTitle("centrality vs z vertex after event selection.");
	TH1D* multafterselection = (TH1D*) eventsafterselection->Project3D("y");
	multafterselection->SetTitle("multiplicity after event selection");
	if(i==0)multafterselection->GetXaxis()->SetRange(1,multafterselection->GetXaxis()->FindBin(200));
// 	if(i==1)multafterselection->GetYaxis()->SetRange(0,multafterselection->GetXaxis()->FindBin(500));
	TH1D* vertexafterselection = (TH1D*)  eventsafterselection->Project3D("x");
	vertexafterselection->SetTitle("vertex after event selection");
	TH1D* centralityafterselection = (TH1D*)  eventsafterselection->Project3D("z");
	centralityafterselection->SetTitle("centrality after event selection");
	c->Clear();
	c->Divide(2,2);
	c->cd(1);
	if(i==0)multvertexbeforeselection->GetYaxis()->SetRange(1,multvertexafterselection->GetYaxis()->FindBin(200));
// 	if(i==1)multvertexbeforeselection->GetYaxis()->SetRange(0,multvertexafterselection->GetYaxis()->FindBin(500));
	multvertexbeforeselection->Draw("surf3");
	multvertexbeforeselection->SetStats(0);
	gPad->SetTheta(30); // default is 30
	gPad->SetPhi(150); // default is 30
	gPad->Update();
	c->cd(2);
	if(i==0)multvertexafterselection->GetYaxis()->SetRange(1,multvertexafterselection->GetYaxis()->FindBin(200));
// 	if(i==1)multvertexafterselection->GetXaxis()->SetRange(0,multafterselection->GetXaxis()->FindBin(500));
	multvertexafterselection->Draw("surf3");
	multvertexafterselection->SetStats(0);
	gPad->SetTheta(30); // default is 30
	gPad->SetPhi(150); // default is 30
	gPad->Update();
	c->cd(3);
	multafterselection->SetStats(0);
	multafterselection->SetTitleOffset(1.1,"xy");
	multafterselection->Draw("E");
	gPad->SetLogy();
	c->cd(4);
	vertexafterselection->SetStats(0);
	vertexafterselection->SetTitleOffset(1.1,"xy");
	vertexafterselection->Draw("E");
	c->Print("MultVertex.eps");
	if(i==1){
	  c->Clear();
	  c->Divide(2,2);
	  c->cd(1);
	  centvertexbeforeselection->Draw("surf3");
	  centvertexbeforeselection->SetStats(0);
	  centvertexbeforeselection->GetYaxis()->SetRange(1,centvertexbeforeselection->GetNbinsX());

	  gPad->SetTheta(30); // default is 30
	  gPad->SetPhi(150); // default is 30
	  gPad->Update();
	  c->cd(2);
	  centvertexafterselection->Draw("surf3");
	  centvertexafterselection->SetStats(0);
	  gPad->SetTheta(30); // default is 30
	  gPad->SetPhi(150); // default is 30	  
	  gPad->Update();
	  c->cd(3);
	  centralityafterselection->SetStats(0);
	  centralityafterselection->SetTitleOffset(1.1,"xy");
	  centralityafterselection->Draw("E");
// 	  gPad->SetLogy();
	  c->cd(4);
	  vertexafterselection->SetStats(0);
	  vertexafterselection->SetTitleOffset(1.1,"xy");
	  vertexafterselection->Draw("E");
	  c->Print("CentVertex.eps");
	}
	delete c;
      }
      //Extract the appropriate Correlation object
      TString corobj(Form("%strigger_correlation_4_8",triggers.Data()));
      TString corobj2(Form("%strigger_correlation_8_16",triggers.Data()));

      if(triggers.CompareTo("Tracks")==0){
	corobj.ReplaceAll("Tracks","track");
	corobj2.ReplaceAll("Tracks","track");

	AliCorrelation3p* signaltrack = dynamic_cast<AliCorrelation3p*>(fList->FindObject(corobj.Data()));
	TString save = TString(corobj.Data());
	corobj.Append("META");
	AliCorrelation3p* mixedtrack = dynamic_cast<AliCorrelation3p*>(fList->FindObject(corobj.Data()));
	corobj = TString(save.Data());
	corobj.Append("METrigger");	
	AliCorrelation3p* mixedtriggertrack = dynamic_cast<AliCorrelation3p*>(fList->FindObject(corobj.Data()));
	corobj = TString(save.Data());
	corobj.Append("MEAllfull");	
	AliCorrelation3p* mixedtrackall = dynamic_cast<AliCorrelation3p*>(fList->FindObject(corobj.Data()));
	if(!signaltrack){
	  signaltrack = dynamic_cast<AliCorrelation3p*>(fList->FindObject(corobj2.Data()));
	  if(!signaltrack) continue;
	  mixedtrack= dynamic_cast<AliCorrelation3p*>(fList->FindObject(Form("%s%s",corobj2.Data(),"META")));
	  mixedtrackall= dynamic_cast<AliCorrelation3p*>(fList->FindObject(Form("%s%s",corobj2.Data(),"MEAllfull")));

	  mixedtriggertrack= dynamic_cast<AliCorrelation3p*>(fList->FindObject(Form("%s%s",corobj2.Data(),"METrigger")));
	}

	if(makefile&&0) signaltrack->MakeResultsFile("");

	if(makefile&&1) signaltrack->MakeResultsFile("");
	if(makefilemixed&&0) mixedtrack->MakeResultsFile("META");
	if(makefileall) mixedtrackall->MakeResultsFile("MEALL");
	if(makefilemixed&&1) mixedtrack->MakeResultsFile("META");
	if(makefiletriggermixed&&0) mixedtriggertrack->MakeResultsFile("METrigger");
	if(makefiletriggermixed&&1) mixedtriggertrack->MakeResultsFile("METrigger");	
      }
//       if(triggers.CompareTo("pi0s")==0){
// 	AliCorrelation3ppi0* signalpi0 = fList->FindObject(corobj.Data());
// 	if(!signalpi0) continue;
// 	corobj.Append("ME");
// 	AliCorrelation3ppi0* mixedpi0 = fList->FindObject(corobj.Data());
// 	if(!mixedpi0)continue;
// 	corobj.Append("Trigger");
// 	AliCorrelation3ppi0* triggermixedpi0 = fList->FindObject(corobj.Data());
// 	if(!triggermixedpi0)continue;
//       }
    }
  }
}

//Prepare the local directories needed:
//   gSystem->MakeDirectory("imgs");
      //prepare the specific folders
//       gSystem->MakeDirectory(Form("imgs/%s",Directory.Data()));
//       if(underdirectory.CompareTo("")!=0){Directory.Append(underdirectory);gSystem->MakeDirectory(Form("imgs/%s",Directory.Data()));}
//       if((underdirectory.CompareTo("/rawmixed")!=0)&&(underdirectory.CompareTo("/binstats")!=0))){gSystem->MakeDirectory(Form("imgs/%s/corrected",Directory.Data()));}