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

	hDraw->SetXTitle("#Delta#varphi (rad)");
	hDraw->GetXaxis()->SetTitleFont(43);
	hDraw->GetXaxis()->SetLabelSize(18);
	hDraw->GetXaxis()->SetLabelFont(43);
	hDraw->GetXaxis()->SetTitleOffset(2.5);
	hDraw->GetXaxis()->SetTitleSize(25);

	if(particle.Contains("Average")) hDraw->SetYTitle("#frac{1}{N_{D}}#frac{dN^{assoc}}{d#Delta#varphi} (rad^{-1})");
	hDraw->GetYaxis()->SetLabelFont(43);
	hDraw->GetYaxis()->SetLabelSize(18);
	hDraw->GetYaxis()->SetTitleFont(43);
	hDraw->GetYaxis()->SetTitleOffset(3.6);
	hDraw->GetYaxis()->SetTitleSize(25);	
      }
    }
    
    Int_t CanvasSizeX, CanvasSizeY;
    if(system=="pp")CanvasSizeX =1030, CanvasSizeY = 1000;
    else if(system=="pPb")CanvasSizeX =1200, CanvasSizeY = 800;
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
      pd->SetLeftMargin(Canvas[k]->GetLeftMargin());
      pd->SetFillStyle(0);
    
      TH1D *hDrL=(TH1D*)hDraw->DrawCopy();
      if(k>=3) {
	hDrL->SetYTitle("");
	hDrL->GetYaxis()->SetLabelSize(0);
      }
      if(k%3!=2) {
	hDrL->SetXTitle("");
	hDrL->GetXaxis()->SetLabelSize(0);
      }

      if(particle.Contains("Average") && k==0){
	TLatex *tlTitle=new TLatex(0.28,0.22,Form("#bf{Average D^{0}, D^{+}, D^{*+}}"));
	tlTitle->SetNDC();
	tlTitle->Draw();
	tlTitle->SetTextFont(63);
	tlTitle->SetTextSize(22);
	if(system=="pPb") {tlTitle->SetX(0.3); tlTitle->SetY(0.19);}

	if(system=="pp") {
	  TLatex *tly=new TLatex(0.45,0.65,Form("#bf{|#Delta#eta|<1.0}, #bf{|#it{y}^{D}|<0.5}"));
	  tly->SetNDC();
	  tly->Draw();
	  tly->SetTextFont(63);
	  tly->SetTextSize(19);

    	/*  TLatex *tlAlice=new TLatex(0.69,0.27,Form("#bf{ALICE}"));
	  tlAlice->SetNDC();
	  tlAlice->Draw();
	  tlAlice->SetTextSize(0.75);*/
	}
	else if(system=="pPb"){
	  TLatex *tly=new TLatex(0.37,0.675,Form("#bf{|#Delta#eta|<1.0}, #bf{-0.96<#it{y}^{D}_{cms}<0.04}"));
	  tly->SetNDC();
	  tly->Draw();
	  tly->SetTextFont(63);
	  tly->SetTextSize(19);

    	/*  TLatex *tlAlice=new TLatex(0.66,0.33,Form("#bf{ALICE}"));
	  tlAlice->SetNDC();
	  tlAlice->Draw();
	  tlAlice->SetTextSize(0.75);*/
	}
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
	  if(system=="pp")hcur->SetMarkerSize(0.6);
	  else if(system=="pPb")hcur->SetMarkerSize(0.6);
	  else Printf("Check your system name");
	  
	  if(particle=="Dzero")hcur->SetMarkerStyle(20);
	  else if(particle=="Dplus")hcur->SetMarkerStyle(21);
	  else if(particle=="Dstar")hcur->SetMarkerStyle(22);
	  else Printf("Check your particle name");
	  
	  if(particle=="Dzero")hcur->SetMarkerColor(kRed);
	  else if(particle=="Dplus")hcur->SetMarkerColor(kGreen+3);
	  else if(particle=="Dstar")hcur->SetMarkerColor(kAzure-2);
	  else Printf("Check your particle name");
        
          hcur->SetLineWidth(1.5);
	  Double_t val=hcur->GetMaximum();
	  hDrL->SetMaximum(val-0.1);
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
	      tl->SetTextFont(63);
	      if(system=="pp")tl->SetTextSize(0.0360);
	      else if(system=="pPb")tl->SetTextSize(0.0360);
	      else Printf("Check your system name");
	      
            }
            
	    // D meson pT range
            if(str.Contains("{p}_{T}^{D}")){
	      tl->SetX(tl->GetX()+0.005);
	      if(system=="pp"){
	  	tl->SetTextFont(63);
		tl->SetTextSize(19);
		if(particle.Contains("Average")) {tl->SetX(0.45); tl->SetY(0.85);}
	      }
	      else if(system=="pPb"){
	  	tl->SetTextFont(63);
		tl->SetTextSize(19);
		if(particle.Contains("Average")) {tl->SetX(0.45); tl->SetY(0.855);}
		else Printf("Check your particle name");
		
	      }
	      else Printf("Check your particle name");
            }

	    // Dplus and Assso pT range
            if(str.Contains("{p}_{T}^{assoc}")){
	      tl->SetX(tl->GetX()+0.005);
  	  /*    if(str.Contains("99")) {
	        if(str.Contains("0.3")) tl->SetTitle("#it{p}_{T}^{assoc}>0.3 GeV/#it{c}");
	        if(str.Contains("1")) tl->SetTitle("#it{p}_{T}^{assoc}>1 GeV/#it{c}");
	      }*/
	      if(system=="pp"){
	  	tl->SetTextFont(63);
		tl->SetTextSize(19);
		tl->SetY(tl->GetY()-0.066);
		if(particle.Contains("Average")) {tl->SetX(0.45); tl->SetY(0.75);}
	      }
	      else if(system=="pPb"){
	  	tl->SetTextFont(63);
		tl->SetTextSize(19);
		if(particle=="Dstar")tl->SetY(tl->GetY()+0.016);
		else if(particle=="Dplus") tl->SetY(tl->GetY()+0.030);
		else if(particle=="Dzero") tl->SetY(tl->GetY()+0.030);
		else if(particle.Contains("Average")) {tl->SetX(0.45); tl->SetY(0.77);}
		else Printf("Check your particle name");
		
	      }
	      else Printf("Check your system name");
            }
            
	    
	    // ETA
            if(str.Contains("eta")){
	      if(system=="pp"){
		tl->SetX(tl->GetX()+0.005);
		tl->SetY(tl->GetY()-0.11); 
	  	tl->SetTextFont(63);
		tl->SetTextSize(19);
		if(particle.Contains("Average")) tl->SetTitle(""); //{tl->SetX(0.2); tl->SetY(0.74);}
	      }
	      else if(system=="pPb"){
		tl->SetX(tl->GetX()+0.005);
	  	tl->SetTextFont(63);
		tl->SetTextSize(19);
		if(particle=="Dstar")tl->SetY(tl->GetY()-0.018);
		else if(particle=="Dplus") tl->SetY(tl->GetY()-0.018);
		else if(particle=="Dzero") tl->SetY(tl->GetY()-0.018);
		else if(particle.Contains("Average")) tl->SetTitle(""); //{tl->SetX(0.2); tl->SetY(0.74);}
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
		else if(particle.Contains("Average"))  {tl->SetX(0.35); tl->SetY(0.53);}
	  	tl->SetTextFont(63);
                tl->SetTextSize(22);
	      }
              
	      else if(system=="pPb"){
		tl->SetX(tl->GetX()+0.25);
		tl->SetY(tl->GetY()+0.04);
	  	tl->SetTextFont(63);
		tl->SetTextSize(22);
		if(particle.Contains("Average"))  {tl->SetX(0.35); tl->SetY(0.575);}
	      }
	      else Printf("Check your particle name");              
            }
         
          if(str.Contains("TeV")){
              TString strTitle=tl->GetTitle();

              if(system=="pp"){
              tl->SetX(tl->GetX()+0.010);
              tl->SetY(tl->GetY()-0.029);
	      tl->SetTextFont(63);
	      tl->SetTextSize(20);
                  
              Printf("%s Line Written here is <--- ",strTitle.Data());
              strTitle.ReplaceAll("pp,","ALICE, pp,");
              strTitle.ReplaceAll(", L_{int} = 5 nb^{-1}","");
              if(particle=="Dstar")  strTitle.ReplaceAll("#sqrt{s","#sqrt{#it{s} ");
              else if(particle=="Dzero")  strTitle.ReplaceAll("#sqrt{s","#sqrt{#it{s} ");
              else if(particle=="Dstar")  strTitle.ReplaceAll("=5.02 TeV"," = 5.02 TeV");
              else if(particle.Contains("Average"))  {tl->SetX(0.28); tl->SetY(0.12);}

              //else if(particle=="Dzero")  strTitle.ReplaceAll("#sqrt{s","#it{s} ");

                tl->SetTitle(strTitle.Data());
               //#bf{pp, #sqrt{#it{s}} = 7 TeV, L_{int} = 5 nb^{-1}}

              }
              else if(system=="pPb"){
		tl->SetX(tl->GetX()+0.03);
		tl->SetY(tl->GetY()-0.62);
	  	tl->SetTextFont(63);
		tl->SetTextSize(20);
        	Printf("%s Line Written here is <--",strTitle.Data());
       	 	strTitle.ReplaceAll(", L_{int} = 50 #mub^{-1}","");
        	strTitle.ReplaceAll("p-Pb,","ALICE, p-Pb,");
        	if(particle=="Dstar")  strTitle.ReplaceAll("#sqrt{s","#sqrt{#it{s}");
        	else if(particle=="Dzero")  strTitle.ReplaceAll("#sqrt{s","#sqrt{#it{s}");
       		else if(particle=="Dstar")  strTitle.ReplaceAll("=7 TeV"," = 7 TeV");
              	else if(particle.Contains("Average"))  {tl->SetX(0.3); tl->SetY(0.11);}

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
	    if(str5.Contains("{p}_{T}^{assoc}")) {
  	 /*     if(str5.Contains("99")) {
	        if(str5.Contains("0.3")) tl->SetTitle("#it{p}_{T}^{assoc}>0.3 GeV/#it{c}");
	        if(str5.Contains("1")) tl->SetTitle("#it{p}_{T}^{assoc}>1 GeV/#it{c}");
	      }*/
	      if(system=="pp") {
	  	tl->SetTextFont(63);
		tl->SetTextSize(19);
		if(particle.Contains("Average")) {tl->SetX(0.45); tl->SetY(0.71);}
                if(k>=3) tl->SetX(tl->GetX()-0.16);
                if(k%3==1) tl->SetY(tl->GetY()+0.05);
                if(k%3==2) tl->SetY(tl->GetY()+0.08);
	      }
	      else if(system=="pPb") {
	  	tl->SetTextFont(63);
		tl->SetTextSize(19);
		if(particle.Contains("Average")) {tl->SetX(0.45); tl->SetY(0.71);}
                if(k>=3) tl->SetX(tl->GetX()-0.2);
                if(k%3==1) tl->SetY(tl->GetY()+0.05);
                if(k%3==2) tl->SetY(tl->GetY()+0.08);
	      }
	      
	    }
            
	    TString str6=tl->GetTitle();
	    if(str6.Contains("scale uncertainty")){
	      if(system=="pp") {
	  	tl->SetTextFont(63);
		tl->SetTextSize(22);
		if(particle.Contains("Average"))  {tl->SetX(0.35); tl->SetY(0.58);}
                if(k>=3) tl->SetX(tl->GetX()-0.16);
                if(k%3==1) tl->SetY(tl->GetY()+0.05);
                if(k%3==2) tl->SetY(tl->GetY()+0.08);
	      }
	      else if(system=="pPb") {
	  	tl->SetTextFont(63);
                tl->SetTextSize(22);
		if(particle.Contains("Average"))  {tl->SetX(0.35); tl->SetY(0.58);}
                if(k>=3) {tl->SetX(tl->GetX()-0.18); tl->SetY(tl->GetY()+0.015);}
                if(k%3==1) tl->SetY(tl->GetY()+0.05);
                if(k%3==2) tl->SetY(tl->GetY()+0.08);
	      }
	      
	    }

	    TString str7=tl->GetTitle();
	    if(str7.Contains("{p}_{T}^{D}")){
	      
	      if(system=="pp") {
	  	tl->SetTextFont(63);
		tl->SetTextSize(19);
		if(particle.Contains("Average")) {tl->SetX(0.45); tl->SetY(0.81);}
                if(k>=3) tl->SetX(tl->GetX()-0.16);
                if(k%3==1) tl->SetY(tl->GetY()+0.05);
                if(k%3==2) tl->SetY(tl->GetY()+0.08);
	      }
	      else if(system=="pPb") {
	  	tl->SetTextFont(63);
		tl->SetTextSize(19);
		if(particle.Contains("Average")) {tl->SetX(0.45); tl->SetY(0.81);}
                if(k>=3) tl->SetX(tl->GetX()-0.2);
                if(k%3==1) tl->SetY(tl->GetY()+0.05);
                if(k%3==2) tl->SetY(tl->GetY()+0.08);
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
	  if(system=="pp")gr->SetMarkerSize(0.6);
	  else if(system=="pPb")gr->SetMarkerSize(0.6);
	  
	  if(particle == "Dzero")gr->SetLineColor(kRed);
	  else if(particle == "Dplus")gr->SetLineColor(kGreen+3);
	  else if(particle == "Dstar")gr->SetLineColor(kAzure-2);

      //gr->SetMarkerStyle(25);

	  gr->SetLineWidth(1.5);
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
	filenames[0] = Form("%s/CanvaAndVariedHisto%s%sPt3to5assocPt0.3to99.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
	filenames[1] = Form("%s/CanvaAndVariedHisto%s%sPt3to5assocPt0.3to1.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
	filenames[2] = Form("%s/CanvaAndVariedHisto%s%sPt3to5assocPt1.0to99.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      }
      filenames[3-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s%sPt5to8assocPt0.3to99.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      filenames[4-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s%sPt5to8assocPt0.3to1.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      filenames[5-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s%sPt5to8assocPt1.0to99.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      filenames[6-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s%sPt8to16assocPt0.3to1.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      filenames[7-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s%sPt8to16assocPt0.3to99.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      filenames[8-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s%sPt8to16assocPt1.0to99.0.root",inputdirectory.Data(),Sys.Data(),Part.Data());
      
  }  
  else {
      if(Sys.EqualTo("pp")){
	filenames[0] = Form("%s/CanvaAndVariedHisto%s_%s_Pt3to5assocPt0.3to99.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
	filenames[1] = Form("%s/CanvaAndVariedHisto%s_%s_Pt3to5assocPt0.3to1.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
	filenames[2] = Form("%s/CanvaAndVariedHisto%s_%s_Pt3to5assocPt1.0to99.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
      }
      filenames[3-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s_%s_Pt5to8assocPt0.3to99.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
      filenames[4-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s_%s_Pt5to8assocPt0.3to1.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
      filenames[5-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s_%s_Pt5to8assocPt1.0to99.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
      filenames[6-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s_%s_Pt8to16assocPt0.3to99.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
      filenames[7-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%s_%s_Pt8to16assocPt0.3to1.0.root",inputdirectory.Data(),Part.Data(),Sys.Data());
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
    c->SaveAs(Form("%s/%s.pdf",outputDir.Data(),plotsout.Data()));
    c->SaveAs(Form("%s/%s.eps",outputDir.Data(),plotsout.Data()));

}





