/* 
 Comparison macro for D*, D0, D+
 in good style
 You can choose 2 plots to be shown,
 as arguments of the main method
 by: Fabio (fabio.colamaria@cern.ch)
*/

TString inputdirectory = "";
TCanvas **Canvas;
TString *filenames;

void SetInputDirectory(TString strdir){
  inputdirectory=strdir;
}

//main function
void DoNiceSpecieComparisonPlot(TString ptD1 ="5to8", TString ptThr1 ="0.3to1.0", TString sys1 = "pp", TString ptD2 ="8to16", TString ptThr2 ="1.0to99.0", TString sys2 = "pPb") {
  gStyle->SetOptStat(0000);
gStyle->SetOptFit(000);

  TString ptD[2] = {ptD1.Data(),ptD2.Data()};
  TString ptThr[2] = {ptThr1.Data(),ptThr2.Data()};
  TString sys[2] = {sys1.Data(),sys2.Data()};

  LoadFileNamesToCompare(ptThr[0],sys[0],ptThr[1],sys[1]);

  Canvas = new TCanvas *[2];
  Canvas[0] = GetCanvas(0,"CanvasCorrelation");   // get p-Pb
  Canvas[0]->SetName("canvaInput0");
  Canvas[1] = GetCanvas(1,"CanvasCorrelation");   // get p-Pb
  Canvas[1]->SetName("canvaInput1");

  TCanvas *cOut = new TCanvas("cOut","Nice Style D-meson comparison",150,150,1500,750);
  cOut->Divide(2,1);
  TPad *pad[2];

  // load the selected pads in the canvas
  for(int i=0; i<=1; i++) {
    cOut->cd(i+1);
    Int_t suffPad = 0;
    if(ptD[i]=="3to5") suffPad=1;  else if(ptD[i]=="5to8") suffPad=2;  else if(ptD[i]=="8to16") suffPad=3;
    if(sys[i]=="pPb") suffPad-=1;
    if(suffPad==0) {printf("Error! Uncorrect pT_D range selected!"); return;}
    pad[i] = (TPad*)Canvas[i]->FindObject(Form("CanvasCorrelation_%d",suffPad));
    pad[i]->SetPad(0,0,1,1);
    pad[i]->Draw();
  }

  printf("Canvas and pads correctly opened and loaded");

  //adjust the graphic layout
  RestylePlot(pad[0],0);
  RestylePlot(pad[1],1);

  // saving the canvases in .root and .png
  TString directname="Output_SngCav_Comparison";
  SaveCanvas(cOut, directname);
  printf("... Plot correctly produced and saved!");

}

//_______________________________________________________________________
void LoadFileNamesToCompare(TString ptThr1 ="0.3to1.0", TString sys1 = "pp", TString ptThr2 ="1.0to99.0", TString sys2 = "pPb"){
 
  filenames = new TString[2];
 
  filenames[0] = Form("%s/Comparison_DHCorrelations_assopT%s_%s.root",inputdirectory.Data(),ptThr1.Data(),sys1.Data());
  filenames[1] = Form("%s/Comparison_DHCorrelations_assopT%s_%s.root",inputdirectory.Data(),ptThr2.Data(),sys2.Data());

}  

//_______________________________________________________________________
TCanvas * GetCanvas(Int_t i, TString canvasname = "cDraw"){

  TString path = filenames[i];
  //cout <<"file #" <<i<<": Reading File from path --->> " << path << endl;
  TFile * file = TFile::Open(path.Data(),"READ");
  TCanvas * c = (TCanvas*)file->Get(canvasname.Data());  
  return c;  
}

//_______________________________________________________________________
void RestylePlot(TPad * pad, Int_t pos){

  TList *lc=pad->GetListOfPrimitives();
  Int_t entries=lc->GetEntries();
  lc->ls();
  Int_t syst=1;
  Int_t nextMeson=0;

  pad->SetTickx();
  pad->SetTicky();
  pad->SetMargin(0.18,0.00,0.10,0.05);
  pad->SetFillStyle(0);

  Double_t maxY[3] = {0.,0.,0.};

  TLegend *leg = new TLegend(0.35,0.48,0.92,0.72);
  leg->SetLineColor(kWhite);

  TLegend *legSuperimp = new TLegend(0.35,0.48,0.92,0.72);
  legSuperimp->SetLineColor(kWhite);
  legSuperimp->SetFillStyle(0);

  TString dZeroUnc = "";
  TString dPlusUnc = "";
  TString dStarUnc = "";    

  TH1D *h1;
  TH1D *h2;
  TH1D *h3;
  
  TH1D *h1Superimp;
  TH1D *h2Superimp;
  TH1D *h3Superimp;

  for(Int_t jl=0;jl<entries;jl++){
    TObject *obj=(TObject*)lc->At(jl);
    TString strName=obj->ClassName();

    printf("Primitive %d is: %s\n",jl,strName.Data());
    if(strName.Contains("TFrame")) {      TFrame *tf=(TFrame*)obj; tf->SetX1(0.2); continue;}

    if(strName.Contains("TLatex")) {
      TLatex *tl=(TLatex*)obj;
      TString str=tl->GetTitle();
      str.ReplaceAll("#bf","#font[42]");
      if(str.Contains("#it{p}_{T}^{D^{0}}")) {tl->SetTitle(""); continue;}

      if(str.Contains("Comparison")) tl->SetTitle("");

      if(str.Contains("#Delta#eta")) tl->SetTitle("");

      if(str.Contains("scale uncertainty") && tl->GetTextColor()==kBlack) tl->SetTitle("");

      if(str.Contains("D^{}  meson-charged")) {
   	/*if(pos==0) {str.ReplaceAll("D^{}  meson","D meson"); tl->SetX(tl->GetX()+0.01); tl->SetTextSize(0.035);}*/
	/*else*/ tl->SetTitle("");
      }

      if(str.Contains("D^{0}"))str.ReplaceAll("D^{0}","D^{0 }");
      if(str.Contains("D^{+}"))str.ReplaceAll("D^{+}","D^{+ }");
      if(str.Contains("TeV")) {
	if(str.Contains("pp")) {
	  syst=0;  //needed to set correctly the y range later
	  tl->SetTitle("pp, #sqrt{#it{s}} = 7 TeV");
        } else tl->SetTitle("p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
	tl->SetX(tl->GetX()+0.04); 
	/*if(pos!=0)*/ tl->SetY(tl->GetY()+0.08); 
	tl->SetTextSize(0.045);
      }

      if(str.Contains("GeV")) {
	tl->SetX(tl->GetX()+0.04);
	/*if(pos!=0)*/ tl->SetY(tl->GetY()+0.065); 
	str.ReplaceAll("GeV/c","GeV/#it{c}");
	if(pos==0) str.ReplaceAll("^{D}","^{D}_{cms}");
	tl->SetTitle(str.Data());
	tl->SetTextSize(0.045);
      }

      if(str.Contains("assoc")) {
	tl->SetX(tl->GetX());
	/*if(pos!=0)*/ tl->SetY(tl->GetY()-0.065); 
	str.ReplaceAll("GeV/c","GeV/#it{c}");
	str.ReplaceAll("1.0","1");
	tl->SetTitle(str.Data());
	tl->SetTextSize(0.045);
      }

      if(str.Contains("scale uncertainty") && tl->GetTextColor()==kRed) {
	tl->SetX(0.5);
	tl->SetY(0.67);
	str.ReplaceAll("#cbar",",");
	str.ReplaceAll("+","#plus"); 
	str.ReplaceAll("-","#minus"); 
	tl->SetTitle("");
	tl->SetTextSize(0.042);
	dZeroUnc=str.Data();
	
      }

      if(str.Contains("scale uncertainty") && tl->GetTextColor()==kAzure-2) {
	tl->SetX(0.5);
	tl->SetY(0.595); 
	str.ReplaceAll("#cbar",","); 
	str.ReplaceAll("+","#plus"); 
	str.ReplaceAll("-","#minus"); 
	tl->SetTitle("");
	tl->SetTextSize(0.042);
	dStarUnc=str.Data();
      }

      if(str.Contains("scale uncertainty") && tl->GetTextColor()==kGreen+3) {
	tl->SetX(0.5);
	tl->SetY(0.52); 
	str.ReplaceAll("#cbar",","); 
	str.ReplaceAll("+","#plus"); 
	str.ReplaceAll("-","#minus"); 
	tl->SetTitle("");
	tl->SetTextSize(0.042);
	dPlusUnc=str.Data();	
      }

    } //end of TLatex 'if'

    if(strName.Contains("TH1")) {
      TH1D *hist = (TH1D*)lc->At(jl); 
      if(hist->GetMarkerColor()==kRed) {
	maxY[0] = hist->GetBinContent(hist->GetMaximumBin());
	hist->SetMarkerSize(1.8);
	hist->GetXaxis()->SetTitle("#Delta#varphi (rad)");
	hist->GetXaxis()->CenterTitle(kTRUE);
	hist->GetXaxis()->SetTitleSize(0.05);
	hist->GetXaxis()->SetTitleOffset(1.00);
	hist->GetXaxis()->SetLabelSize(0.04);
	hist->GetYaxis()->SetTitle("#frac{1}{#it{N}_{D}} #frac{d#it{N}^{assoc}}{d#Delta#varphi} (rad^{-1})");
	hist->GetYaxis()->CenterTitle(kTRUE);
	hist->GetYaxis()->SetTitleSize(0.05);
	hist->GetYaxis()->SetTitleOffset(1.14);
	hist->GetYaxis()->SetLabelSize(0.04);
	hist->SetLineColor(kRed);	
	hist->SetLineWidth(1);	
        if(dZeroUnc.CompareTo("")) {
          TLegendEntry* l1 = leg->AddEntry(hist,Form("%s",dZeroUnc.Data()),"elp");
          l1->SetTextColor(kRed);
        } 
	h1 = (TH1D*)hist->Clone();
        h1Superimp = (TH1D*)hist->Clone();
        h1Superimp->SetMarkerStyle(kCircle);
        h1Superimp->SetMarkerColor(kRed+1);
printf("GetBinContent %f, %f\n",hist->GetBinContent(2),h1Superimp->GetBinContent(2));
        if(dZeroUnc.CompareTo("")) {
          TLegendEntry* l1sup = legSuperimp->AddEntry(h1Superimp,"","elp");
          l1sup->SetTextColor(kRed);          	  
	}
        nextMeson=1;  
      }
      if(hist->GetMarkerColor()==kAzure-2) {
	maxY[1] = hist->GetBinContent(hist->GetMaximumBin());
	hist->SetMarkerSize(2.5);	
	hist->SetMarkerStyle(33);
	hist->SetLineColor(kAzure-2);
	hist->SetLineWidth(1);		
	h2 = (TH1D*)hist->Clone();
        h2Superimp = (TH1D*)hist->Clone();
        h2Superimp->SetMarkerStyle(27);
        h2Superimp->SetMarkerColor(kAzure+3);
	nextMeson=2;
        TLegendEntry* l2 = leg->AddEntry(hist,Form("%s",dStarUnc.Data()),"elp");
        l2->SetTextColor(kAzure-2);	
        TLegendEntry* l2sup = legSuperimp->AddEntry(h2Superimp,"","elp");
        l2sup->SetTextColor(kAzure-2);    	
      }
      if(hist->GetMarkerColor()==kGreen+3) {
	maxY[2] = hist->GetBinContent(hist->GetMaximumBin());
	hist->SetMarkerSize(1.8);
	hist->SetLineColor(kGreen+3);
	hist->SetLineWidth(1);	
	h3 = (TH1D*)hist->Clone();	
        h3Superimp = (TH1D*)hist->Clone();
        h3Superimp->SetMarkerStyle(25);
        h3Superimp->SetMarkerColor(kGreen+4);
	nextMeson=3;
        TLegendEntry* l3 = leg->AddEntry(hist,Form("%s",dPlusUnc.Data()),"elp");
        l3->SetTextColor(kGreen+3);		
        TLegendEntry* l3sup = legSuperimp->AddEntry(h3Superimp,"","elp");
        l3sup->SetTextColor(kGreen+3);	
      }

    } //end of TH1D 'if'

    if(strName.Contains("TGraphAsymmErrors")) {
      TGraphAsymmErrors *errBox = (TGraphAsymmErrors*)lc->At(jl);  
      if(nextMeson==1) {
	errBox->SetLineColor(kRed);
	errBox->SetLineWidth(1);
	nextMeson=0;
      }
      if(nextMeson==2) {
	errBox->SetLineColor(kAzure-2);
	errBox->SetLineWidth(1);
	nextMeson=0;
      }
      if(nextMeson==3) {
	errBox->SetLineColor(kGreen+3);
	errBox->SetLineWidth(1);
	nextMeson=0;
      }

    } //end of TGraphAsymmErrors 'if'


  } //end of cycle on primitives

  pad->cd();
  if(pos==0 || pos==1){ //write to both panels (IRC comment)
	TLatex *tlAlice=new TLatex(0.82,0.875,Form("ALICE"));
	tlAlice->SetNDC();
	tlAlice->Draw();
	tlAlice->SetTextSize(0.050);
      }  
/*
      if(syst==0) {
	  TLatex *tly=new TLatex(0.19,0.74,Form("#bf{|#Delta#eta|<1.0}, #bf{|#it{y}^{D}|<0.5}"));
	  tly->SetNDC();
	  tly->Draw();
	  tly->SetTextSize(0.035);
      }
      else {
	  TLatex *tly=new TLatex(0.19,0.74,Form("#bf{|#Delta#eta|<1.0}, #bf{-0.96<#it{y}^{D}_{cms}<0.04}"));
	  tly->SetNDC();
	  tly->Draw();
	  tly->SetTextSize(0.035);
      }
*/
  Double_t maxFinalY = TMath::Max(maxY[0],TMath::Max(maxY[1],maxY[2])) * 1.9;
  maxFinalY = (TMath::Floor(maxFinalY*2)+1)/2.;

  /*if(pos==0)*/ leg->Draw("same");
  legSuperimp->Draw("same");

  TH1D *histFin = (TH1D*)pad->FindObject("hDataCorrectedTempl0CentrFpromptReflected");
  histFin->SetMaximum(maxFinalY);
  histFin->GetYaxis()->SetTitleOffset(histFin->GetYaxis()->GetTitleOffset()+0.3);

  h1->Draw("same");
  h1Superimp->Draw("same");
  h3->Draw("same");
  h3Superimp->Draw("same");
  h2->Draw("same");
  h2Superimp->Draw("same");

}  

//_______________________________________________________________________
void SaveCanvas(TCanvas * c, TString directory){
  
  TString outputDir = "";
  TString nameoutput = "Comparison_DHCorrelations_NiceStyle";  
   
  c->SaveAs(Form("%s/%s.root",directory.Data(),nameoutput.Data()));
  c->SaveAs(Form("%s/%s.png",directory.Data(),nameoutput.Data()));
  c->SaveAs(Form("%s/%s.eps",directory.Data(),nameoutput.Data()));
  c->SaveAs(Form("%s/%s.pdf",directory.Data(),nameoutput.Data()));
  
}




