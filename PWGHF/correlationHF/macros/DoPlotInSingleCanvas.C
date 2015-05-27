//##################################################################
//
//  Macro managing the plotting of single mesons and average results
//  for different kinematic ranges in single canvases
//    Jitendra Kumar, jitendra.kumar@cern.ch
// 
//##################################################################

TString inputdirectory = "";
Int_t nhistos=6; //pPb !! IMPORTANT to change -->9 for pp (this is done already in the LoadFile method below)
TString *filenames; 
TCanvas **Canvas;
TCanvas *coutput;

void SetInputDirectory(TString strdir){
  inputdirectory=strdir;
}
void DoAllPlotInSingleCanvasStandardPaths(Int_t avopt){// standard path --> average plots are in directory = inputdirectory/Averages
    
    DoPlotInSingleCanvas("Dzero","pp");
    DoPlotInSingleCanvas("Dstar","pp");
    DoPlotInSingleCanvas("Dplus","pp");
    DoPlotInSingleCanvas("Dzero","pPb");
    DoPlotInSingleCanvas("Dstar","pPb");
    DoPlotInSingleCanvas("Dplus","pPb");

    inputdirectory.Append("/Averages");
    if(avopt==0||avopt==2){
      DoPlotInSingleCanvas("WeightedAverageDzeroDstarDplus","pp");
      DoPlotInSingleCanvas("WeightedAverageDzeroDstarDplus","pPb");
    }
    if(avopt==1||avopt==2){
      DoPlotInSingleCanvas("ArithmeticAverageDzeroDstarDplus","pp");
      DoPlotInSingleCanvas("ArithmeticAverageDzeroDstarDplus","pPb");
    }
    
}

void Reset(){
  delete [] filenames;
  delete filenames;

  //  for(Int_t k=0;k<nhistos;k++){
    //    if(Canvas[k])delete Canvas[k];
  //  }
  delete [] Canvas;
  delete coutput;
}

//_______________________________________________________________________
void DoPlotInSingleCanvas(TString particle ="Dstar",TString system = "pp"){
    
  gStyle->SetOptStat(0);
  Reset();  
  LoadFileNames(particle,system);
  cout <<" -------------------> Started plotting in Single Canvas" << endl;  
  

  Canvas = new TCanvas*[nhistos];
  Canvas[0] = GetCanvas(0,"cDraw");   // get p-Pb/pp histo canvas

  TList *l0=Canvas[0]->GetListOfPrimitives();
  TH1D *hDraw;
  Int_t entries=l0->GetEntries();
  for(Int_t jl=0;jl<entries;jl++){
    TObject *obj=l0->At(jl);
      TString strName=obj->ClassName();
      if(strName.Contains("TH1")){
	//Printf("Found the histo");
	hDraw=(TH1D*)obj->Clone("hDraw");
	hDraw->Reset(0);

	hDraw->SetXTitle("#Delta#varphi  (rad)");
	hDraw->GetXaxis()->SetTitleFont(42);
	hDraw->GetXaxis()->SetLabelFont(42);
	hDraw->GetXaxis()->SetTitleOffset(1.1);
	hDraw->GetXaxis()->SetTitleSize(0.042);

    
	if(particle=="Dstar")hDraw->SetYTitle("#frac{1}{#it{N}_{D^{*+}}}#frac{d#it{N}^{assoc}}{d#Delta#varphi} (rad^{-1})");
	hDraw->GetYaxis()->SetLabelFont(42);
	hDraw->GetYaxis()->SetTitleFont(42);
	hDraw->GetYaxis()->SetTitleOffset(1.3);
	hDraw->GetYaxis()->SetTitleSize(0.042);	
      }
    }
    
    Int_t CanvasSizeX, CanvasSizeY;
    if(system=="pp")CanvasSizeX =1030, CanvasSizeY = 1000;
    else if(system=="pPb")CanvasSizeX =1260, CanvasSizeY = 800;
    else cout << "Check System input !! " <<endl;
    coutput=new TCanvas("coutput","coutput",CanvasSizeX,CanvasSizeY);
    if(system=="pp")coutput->Divide(3,3);
    else if(system=="pPb")coutput->Divide(3,2);
    else cout << "Check System input !! " <<endl;
    
    coutput->Draw();
    coutput->SetLeftMargin(Canvas[0]->GetLeftMargin()-0.08);
    coutput->SetRightMargin(Canvas[0]->GetRightMargin()+0.08);
    //coutput->SetFrameBorderMode(0); // new to May 27
    // loop on histos
    for(Int_t k=0; k<nhistos; k++){
      
      if(k>0)Canvas[k] = GetCanvas(k,"cDraw");   // get p-Pb
      TList *lc=Canvas[k]->GetListOfPrimitives();
      Canvas[k]->SetName(Form("canvaInput%d",k));
      //Canvas[k]->Draw();
      coutput->cd(k+1);
      if(k>0) entries=lc->GetEntries();
      
      TPad *pd=(TPad*)coutput->cd(k+1);
      pd->SetTickx();
      pd->SetTicky();
      //pd->SetMargin(0,0,0,0);
      pd->SetLeftMargin(Canvas[k]->GetLeftMargin());
      //pd->SetRightMargin(Canvas[k]->GetRightMargin()+0.08);
      
      TH1D *hDrL=(TH1D*)hDraw->DrawCopy();

      if(particle.Contains("Average") && k==0){
	TLatex *tlTitle=new TLatex(0.18,0.30,Form("#bf{Average D^{0}, D^{+}, D^{*+}}"));
	tlTitle->SetNDC();
	tlTitle->Draw();
	tlTitle->SetTextSize(0.0425);
      }
      
      for(Int_t jl=0;jl<entries;jl++){
	TObject *obj=lc->At(jl);
	TString strName=obj->ClassName();
	
	if(strName.Contains("TFrame")){
	  continue;}
	
	if(strName.Contains("TPad")){
	  continue;}
	
	if(strName.Contains("TH1")){
	  TH1D *hcur=(TH1D*)obj;
	  if(system=="pp")hcur->SetMarkerSize(0.7);
	  else if(system=="pPb")hcur->SetMarkerSize(0.9);
	  else Printf("Check your system name");
	  
	  if(particle=="Dzero")hcur->SetMarkerStyle(20);
	  else if(particle=="Dplus")hcur->SetMarkerStyle(21);
	  else if(particle=="Dstar")hcur->SetMarkerStyle(22);
	  else Printf("Check your particle name");
	  
	  if(particle=="Dzero")hcur->SetMarkerColor(kRed);
	  else if(particle=="Dplus")hcur->SetMarkerColor(kGreen+3);
	  else if(particle=="Dstar")hcur->SetMarkerColor(kAzure-2);
	  else Printf("Check your particle name");
        
      hcur->SetLineWidth(2);
	  Double_t val=hcur->GetMaximum();
	  hDrL->SetMaximum(val);
	}
	

	// Style for first histogram 

	if(k==0){
	  
	  if(strName.Contains("TLatex")){
            TLatex *tl=(TLatex*)obj;
            TString str=tl->GetTitle();
            
	    // Title 
            if(str.Contains("charged")){
	      if(particle=="Avg") {
		tl->SetTitle("#bf{D meson - charged particle azimuthal correlations}");
	      }
	      tl->SetX(tl->GetX()+0.010);
	      tl->SetY(tl->GetY()-0.015);
	      if(system=="pp")tl->SetTextSize(0.0360);
	      else if(system=="pPb")tl->SetTextSize(0.0360);
	      else Printf("Check your particle name");
	      
            }
            
	    // Dplus and Assso pT range
            if(str.Contains("{p}_{T}^{assoc}")){
	      tl->SetX(tl->GetX()+0.005);
	      if(system=="pp"){
		tl->SetTextSize(0.045);
		tl->SetY(tl->GetY()-0.066);
	      }
	      else if(system=="pPb"){
		tl->SetTextSize(0.043);
		if(particle=="Dstar")tl->SetY(tl->GetY()+0.016);
		else if(particle=="Dplus") tl->SetY(tl->GetY()+0.030);
		else if(particle=="Dzero") tl->SetY(tl->GetY()+0.030);
		else Printf("Check your particle name");
		
	      }
	      else Printf("Check your particle name");
            }
            
	    
	    // ETA
            if(str.Contains("eta")){
	      if(system=="pp"){
		tl->SetX(tl->GetX()+0.005);
		tl->SetY(tl->GetY()-0.098); 
		tl->SetTextSize(0.045);
	      }
	      else if(system=="pPb"){
		tl->SetX(tl->GetX()+0.005);
		tl->SetTextSize(0.043);
		if(particle=="Dstar")tl->SetY(tl->GetY()-0.018);
		else if(particle=="Dplus") tl->SetY(tl->GetY()-0.018);
		else if(particle=="Dzero") tl->SetY(tl->GetY()-0.018);
		else Printf("Check your particle name");
	      }
	      else Printf("Check your particle name"); 
            }
	    
	    // SCALE Uncertainity 
            if(str.Contains("scale uncertainty")){
	      if(system=="pp"){
                tl->SetX(tl->GetX()+0.23);
		if(particle=="Dstar") tl->SetY(tl->GetY()-0.06);
		else if(particle=="Dplus")  tl->SetY(tl->GetY()-0.04);
		else if(particle=="Dzero")  tl->SetY(tl->GetY()-0.04);
                tl->SetTextSize(0.053);
	      }
              
	      else if(system=="pPb"){
		tl->SetX(tl->GetX()+0.25);
		tl->SetY(tl->GetY()+0.04);
		tl->SetTextSize(0.047);
	      }
	      else Printf("Check your particle name");              
            }
            
	    
	    if(str.Contains("ALICE Preliminary")){
	      if(system=="pp"){
		tl->SetX(tl->GetX()-0.28);
		tl->SetY(tl->GetY()-0.467);
		tl->SetTextSize(0.058);
	      }
	      else if(system=="pPb"){
		tl->SetX(tl->GetX()-0.26);
		tl->SetY(tl->GetY()-0.37);
		tl->SetTextSize(0.052);
	      }
	      else Printf("Check your particle name");
	      
            }
          
          
          if(str.Contains("TeV")){
              TString strTitle=tl->GetTitle();

              if(system=="pp"){
              tl->SetX(tl->GetX()+0.010);
              tl->SetY(tl->GetY()-0.029);
	          tl->SetTextSize(0.045);
                  
              Printf("%s Line Written here is <--- ",strTitle.Data());
              strTitle.ReplaceAll("pp,","pp, ");
              strTitle.ReplaceAll("L_{int}","#it{L}_{int} ");
              if(particle=="Dstar")  strTitle.ReplaceAll("#sqrt{s","#sqrt{#it{s} ");
              else if(particle=="Dzero")  strTitle.ReplaceAll("#sqrt{s","#sqrt{#it{s} ");
                  
              if(particle=="Dstar")  strTitle.ReplaceAll("=5.02 TeV"," = 5.02 TeV");

              //else if(particle=="Dzero")  strTitle.ReplaceAll("#sqrt{s","#it{s} ");

                tl->SetTitle(strTitle.Data());
               //#bf{pp, #sqrt{#it{s}} = 7 TeV, L_{int} = 5 nb^{-1}}

              }
              else if(system=="pPb"){
		tl->SetX(tl->GetX()+0.03);
		tl->SetY(tl->GetY()-0.62);
		tl->SetTextSize(0.042);
        Printf("%s Line Written here is <--",strTitle.Data());
        strTitle.ReplaceAll("L_{int}"," #it{L}_{int}");
        strTitle.ReplaceAll("p-Pb,","p-Pb,  ");
        if(particle=="Dstar")  strTitle.ReplaceAll("#sqrt{s","#sqrt{#it{s}");
        else if(particle=="Dzero")  strTitle.ReplaceAll("#sqrt{s","#sqrt{#it{s}");

        if(particle=="Dstar")  strTitle.ReplaceAll("=7 TeV"," = 7 TeV");

        tl->SetTitle(strTitle.Data());
    };
              //else Printf("Check your particle name");
          }
          
          
          
	  }
	}
        
        
	
	if(k>0){
	  if(strName.Contains("TLatex")){
	    TLatex *tl=(TLatex*)obj;
	    TString str1=tl->GetTitle();
	    if(str1.Contains("charged")){
	      continue;
	    }
	    TString str2=tl->GetTitle();
	    if(str2.Contains("TeV")){
	      continue;
	    }
	    TString str3=tl->GetTitle();
	    if(str3.Contains("ALICE")){
	      continue;
	    }
            
	    TString str4=tl->GetTitle();
	    if(str4.Contains("|#Delta#eta|")){
	      continue;
	    }
            
	    TString str5=tl->GetTitle();
	    if(str4.Contains("{p}_{T}^{assoc}")){
	      
	      if(system=="pp") {
		tl->SetX(tl->GetX()+0.03);
		tl->SetY(tl->GetY()+0.020);
		tl->SetTextSize(0.044);
	      }
	      else if(system=="pPb") {
		tl->SetX(tl->GetX()+0.03);
		tl->SetY(tl->GetY()+0.050);
		tl->SetTextSize(0.042);
	      }
	      
	    }
            
	    TString str6=tl->GetTitle();
	    if(str6.Contains("scale uncertainty")){
	      if(system=="pp") {
		tl->SetX(tl->GetX()+0.16);
		tl->SetY(tl->GetY()-0.02);
		tl->SetTextSize(0.060);
	      }
	      else if(system=="pPb") {
                tl->SetX(tl->GetX()+0.23);
                tl->SetY(tl->GetY()+0.06);
                tl->SetTextSize(0.047);
	      }
	      
	    }
	    
	  }
	}
	
	coutput->cd(k+1);
	
	if(strName.Contains("TGraph")){
	  //Printf("There is a TGraph");
	  TGraphAsymmErrors *gr=(TGraphAsymmErrors*)obj;
	  //   gr->SetMarkerStyle();
	  //   gr->SetMarkerColor(0);
	  if(system=="pp")gr->SetMarkerSize(0.7);
	  else if(system=="pPb")gr->SetMarkerSize(0.9);
	  
	  if(particle == "Dzero")gr->SetLineColor(kRed);
	  else if(particle == "Dplus")gr->SetLineColor(kGreen+3);
	  else if(particle == "Dstar")gr->SetLineColor(kAzure-2);

      //gr->SetMarkerStyle(25);

	  gr->SetLineWidth(2);
	  gr->Draw("E2");
	  
	  TGraphAsymmErrors *gr2=(TGraphAsymmErrors*)gr->Clone("grHelp");
	  //gr2->SetMarkerColor(kBlack);
	  //gr2->SetLineColor(kBlack);
	  //gr2->SetMarkerSize(0.3);
	  //gr2->Draw("p");
	}
	else  obj->Draw("same");
        
      }
    }
    
    
    coutput->Draw();
    
    // saving the canvases in .root and .png
    TString ptoutput="", direcname="";
    ptoutput += "CanvasOf";
    direcname += "Output_Plots/";
    direcname += particle;
    
    // changes name of final canvas/histo
    if(particle == "Dzero") ptoutput += "_DzeroHadronCorrelations";
    if(particle == "Dstar") ptoutput += "_DstarHadronCorrelations";
    if(particle == "Dplus") ptoutput += "_DplusHadronCorrelations";
    if(particle.Contains("Average")) ptoutput += particle.Data();
    // changes name of final canvas/histo
    if(system == "pp") ptoutput += "_pp";
    if(system == "pPb") ptoutput += "_pPb";
    
    SaveCanvas(coutput, direcname, ptoutput);
    Printf("...  .. .. Done !");

    return;
    
    
}


////_______________________________________________________________________
void LoadFileNames(TString Part ="Dplus",TString Sys = "pPb"){

  Int_t skip3to5pPb=0;
  if(Sys.EqualTo("pp")){
    nhistos=9;
  }
  else if(Sys.EqualTo("pPb")){
    nhistos=6; 
    skip3to5pPb=3;
  }
  filenames=new TString[nhistos];


  if(!(Part.Contains("Average"))){
      
      Printf(" -------------------> Adding %s in %s coll",Part.Data(),Sys.Data());
      if(Sys.EqualTo("pp")){
	filenames[0] = Form("%s/CanvaAndVariedHisto%s%sPt3to5assocPt0.3to1.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
	filenames[1] = Form("%s/CanvaAndVariedHisto%s%sPt3to5assocPt0.3to99.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
	filenames[2] = Form("%s/CanvaAndVariedHisto%s%sPt3to5assocPt1.0to99.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      }
      filenames[3-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s%sPt5to8assocPt0.3to1.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      filenames[4-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s%sPt5to8assocPt0.3to99.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      filenames[5-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s%sPt5to8assocPt1.0to99.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      filenames[6-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s%sPt8to16assocPt0.3to1.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      filenames[7-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s%sPt8to16assocPt0.3to99.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      filenames[8-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s%sPt8to16assocPt1.0to99.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      
  }  
  else {
      if(Sys.EqualTo("pp")){
	filenames[0] = Form("%s/CanvaAndVariedHisto%s_%s_Pt3to5assocPt0.3to1.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
	filenames[1] = Form("%s/CanvaAndVariedHisto%s_%s_Pt3to5assocPt0.3to99.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
	filenames[2] = Form("%s/CanvaAndVariedHisto%s_%s_Pt3to5assocPt1.0to99.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
      }
      filenames[3-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s_%s_Pt5to8assocPt0.3to1.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
      filenames[4-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s_%s_Pt5to8assocPt0.3to99.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
      filenames[5-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s_%s_Pt5to8assocPt1.0to99.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
      filenames[6-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s_%s_Pt8to16assocPt0.3to1.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
      filenames[7-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s_%s_Pt8to16assocPt0.3to99.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
      filenames[8-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s_%s_Pt8to16assocPt1.0to99.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());  
  }    
  
}

//_______________________________________________________________________
TCanvas * GetCanvas(Int_t i, TString canvasname = "cDraw"){
  
  //load the histogram with
  TString path = filenames[i];
  cout << "Reading File from path: " << path << endl;  
  TFile * file = TFile::Open(path.Data(),"READ");
  TCanvas * c = (TCanvas*)file->Get(canvasname.Data());
  
  return c;
  
}


//_______________________________________________________________________
void SaveCanvas(TCanvas * c, TString directory, TString name){
    
    TString outputDir = "";//Plots/15_May/";
    //    gSystem->Exec(Form("mkdir %s",outputDir.Data()));
    
    if(directory != ""){outputDir += directory;
        TString exec = "mkdir -p ";
        exec += outputDir;
        cout << exec << endl;
        gSystem->Exec(exec.Data());
    }
    
    
    TString plotsout = "";//"Canvas_pT_05_";
    plotsout += name;
    c->SaveAs(Form("%s/%s.root",outputDir.Data(),plotsout.Data()));
    c->SaveAs(Form("%s/%s.png",outputDir.Data(),plotsout.Data()));
    c->SaveAs(Form("%s/%s.eps",outputDir.Data(),plotsout.Data()));

}





