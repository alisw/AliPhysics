/* 
 Comparoson macro for D*, D0, D+
 SYSTEM:  p-Pb at 5.02TeV
 Dmeson pT : High and mid pT 
 Associated pT : 0.3-1.0 GeV/c
 by: Jitendra (jitendra.kumar@cern.ch)
*/
//using namespace std::

TString inputdirectory = "";
Int_t nhistos=6; //pPb !! IMPORTANT to change -->9 for pp
void SetInputDirectory(TString strdir){
  inputdirectory=strdir;
}

TString *filenames;
TCanvas **Canvas;
TCanvas *CanvasCorrelation;

void Reset(){
  delete [] filenames;
  delete filenames;

  //  for(Int_t k=0;k<nhistos;k++){
    //    if(Canvas[k])delete Canvas[k];
  //  }
  delete [] Canvas;
  delete CanvasCorrelation;
}

void DoPlotComparedot3to1pPb(){
  DoPlotSingle1Canvas("0.3to1.0","pPb");
}

//main function
void DoPlotSingle1Canvas(TString TrackpTthr ="0.3to1.0",TString system = "pPb"){
  
  gStyle->SetOptStat(0);
  Int_t CanvasSizeX, CanvasSizeY;
  if(system=="pp")CanvasSizeX =1030, CanvasSizeY = 1000;
  else if(system=="pPb")CanvasSizeX =1200, CanvasSizeY = 600;
  else cout << "Check System input !! " <<endl;
  CanvasCorrelation=new TCanvas("CanvasCorrelation","CanvasCorrelation",CanvasSizeX,CanvasSizeY);
  CanvasCorrelation->SetLeftMargin(0.18409396);
  //CanvasCorrelation->SetRightMargin(0.05872483);
  //CanvasCorrelation->SetTopMargin(0.07678883);
  //CanvasCorrelation->SetBottomMargin(0.1239092);
  
  if(system=="pp")CanvasCorrelation->Divide(3,3);
  else if(system=="pPb")CanvasCorrelation->Divide(2,1);
  else cout << "Check System input !! " <<endl;
  //cout->Draw();
  //cout->SetLeftMargin(Canvas[0]->GetLeftMargin()-0.08);
  //cout->SetRightMargin(Canvas[0]->GetRightMargin()+0.08);
  //cout->SetFrameBorderMode(0);
  
  TString CanvasName = "";

  LoadFileNamesToCompare(TrackpTthr,system);
  cout << " " << endl << "Adding file.." << endl;
  for(Int_t k = 0; k<nhistos; k++)cout << k+1 <<" ---------> # " << filenames[k].Data()<<endl;
  
  Int_t entries;
  Canvas = new TCanvas *[nhistos];
  Int_t kTemp = 0;
  for(Int_t k=0; k<nhistos; k++){
    
    gPad->SetLeftMargin(0.15009396);
    gPad->SetRightMargin(-0.12409396);
    
    if(k<3)kTemp=1;
    else if(3<=k<6)kTemp=2;
    else if(6<=k)kTemp=3;
    
    Canvas[k] = GetCanvas(k,"cDraw");   // get p-Pb
    Canvas[k]->SetName(Form("canvaInput%d",k));
    
    TList *lc=Canvas[k]->GetListOfPrimitives();
    entries=lc->GetEntries();
    
    for(Int_t jl=0;jl<entries;jl++){
      TObject *obj=lc->At(jl);
      TString strName=obj->ClassName();
      if(strName.Contains("TFrame"))continue;
      if(strName.Contains("TPad"))continue;

      if(strName.Contains("TH1")){
	TH1D *hcur=(TH1D*)obj;
	if(system=="pp")hcur->SetMarkerSize(0.7);
	else if(system=="pPb")hcur->SetMarkerSize(0.9);
	else Printf("Check your system name");
      }
      
      if(strName.Contains("TLatex")){
	TLatex *tl=(TLatex*)obj;
	TString str=tl->GetTitle();
	if(str.Contains("D^{0}")){
	  hcur->SetMarkerColor(kRed);
	  hcur->SetMarkerStyle(20);
	}else if(str.Contains("D^{+}")){
	  hcur->SetMarkerColor(kGreen+3);
	  hcur->SetMarkerStyle(21);
	}else if(str.Contains("D^{*+}")){
	  // hcur->SetMarkerColor(5);
	  hcur->SetMarkerColor(kAzure-2);
	  hcur->SetMarkerStyle(22);
	}
      }
      
      if(strName.Contains("TGraph")){
	//Printf("There is a TGraph");
	TGraphAsymmErrors *gr=(TGraphAsymmErrors*)obj;
	//   gr->SetMarkerStyle();
	//   gr->SetMarkerColor(0);
	if(system=="pp")gr->SetMarkerSize(0.7);
	else if(system=="pPb")gr->SetMarkerSize(0.9);
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
    
    CanvasCorrelation->cd(kTemp);
    
    if(k==0 || k==3 || k==6){  
      hcur->GetYaxis()->SetRangeUser(0,12);
      hcur->GetYaxis()->SetTitle("#frac{1}{#it{N}_{D^{}}}#frac{d#it{N}^{assoc}}{d#Delta#varphi} (rad^{-1})");
      hcur->Draw("ep");
      
      //dot3 to 1 = 0.15, 0.21, 027
      TLatex *tlTitleDplus=new TLatex(0.48,0.15,"#bf{D^{+}#cbar{}^{+20%}_{-18%} scale uncertainty}");
      tlTitleDplus->SetNDC();
      tlTitleDplus->SetTextColor(kGreen+3);
      tlTitleDplus->SetTextSize(0.033);
      tlTitleDplus->Draw();
      
      TLatex *tlTitleDStar=new TLatex(0.48,0.21,"#bf{D^{*+}#cbar{}^{+16%}_{-13%} scale uncertainty}");
      tlTitleDStar->SetNDC();
      tlTitleDStar->SetTextColor(kAzure-2);
      tlTitleDStar->SetTextSize(0.033);
      tlTitleDStar->Draw();
      
      TLatex *tlTitleDZero=new TLatex(0.48,0.27,"#bf{D^{0}#cbar{}^{+16%}_{-13%} scale uncertainty}");
      tlTitleDZero->SetNDC();
      tlTitleDZero->SetTextColor(kRed);
      tlTitleDZero->SetTextSize(0.033);
      tlTitleDZero->Draw();
      
      if(kTemp==1)TLatex *tlTitle00=new TLatex(0.18,0.75,"#bf{5 < #it{p}_{T}^{D^{0}} < 8 GeV/c}, #bf{-0.96 < #it{y}^{D}_{cms} < 0.04}");
      if(kTemp==2)TLatex *tlTitle00=new TLatex(0.18,0.75,"#bf{8 < #it{p}_{T}^{D^{0}} < 16 GeV/c}, #bf{-0.96 < #it{y}^{D}_{cms} < 0.04}");
      
      tlTitle00->SetNDC();
      //tlTitle0->SetTextColor(kRed);
      tlTitle00->SetTextSize(0.033);
      tlTitle00->Draw();
      
      for(Int_t jl=0;jl<entries;jl++){      
	TObject *obj=lc->At(jl);
	TString strName=obj->ClassName();
	TString str=obj->GetTitle();
	
	if(strName.Contains("TLatex")){
	  TLatex *tl=(TLatex*)obj;
	  
	  if(str.Contains("-charged")){
	    TString strTitle=tl->GetTitle();
	    Printf("\n%s <-- Default txt(1)",strTitle.Data());
	    strTitle.ReplaceAll("D^{0}","D^{} ");
	    Printf("%s <-- Replaced with txt(1)",strTitle.Data());
	    tl->SetTitle(strTitle.Data());
	  }
          
	  if(str.Contains("ALICE")){
	    TString strTitle2=tl->GetTitle();
	    Printf("\n%s <-- Default txt(2)",strTitle2.Data());
	    strTitle2.ReplaceAll("ALICE Preliminary","");
	    Printf("%s <-- Replaced with txt(2)",strTitle2.Data());
	    tl->SetTitle(strTitle2.Data());
	  }
          
	  if(str.Contains("pp")){
	    TString strTitle3=tl->GetTitle();
	    Printf("\n%s <-- Default txt(3)",strTitle3.Data());
	    strTitle3.ReplaceAll("pp, #sqrt{s}=7 TeV, L_{int} = 5 nb^{-1}","p-Pb, #sqrt{s}=5.02 TeV, L_{int} = 50 #mub^{-1}");
	    Printf("%s <-- Replaced with txt(3)",strTitle3.Data());
	    tl->SetTitle(strTitle3.Data());
	  }
          
	  if(str.Contains("assoc")){
	    TString strTitle5=tl->GetTitle();
	    Printf("\n%s <-- Default txt(4)",strTitle5.Data());
	    if(kTemp==1)strTitle5.ReplaceAll(strTitle5.Data(),"#bf{0.3 < #it{p}_{T}^{assoc} < 1 GeV/#it{c}}, #bf{|#Delta#eta| < 1.0}");
	    if(kTemp==2)strTitle5.ReplaceAll(strTitle5.Data(),"#bf{0.3 < #it{p}_{T}^{assoc} < 1 GeV/#it{c}}, #bf{|#Delta#eta| < 1.0}");
	    Printf("%s <-- Replaced with txt(4)",strTitle5.Data());
	    tl->SetTitle(strTitle5.Data());
	  }
          
	  // SCALE Uncertainity
	  if(str.Contains("scale uncertainty")){
	    TString strTitle4=tl->GetTitle();
	    strTitle4.ReplaceAll("{}^{+16%}_{-13%} scale uncertainty","D^{0}, D^{+}, D^{*+} Comparison");
	    tl->SetTitle(strTitle4.Data());
	    tl->SetX(tl->GetX()+0.10); // +0.10 for dot3to1
	    tl->SetY(tl->GetY()+0.01);
	    tl->SetTextSize(0.023);
	    
	    TString strTitle41=tl->GetTitle();
	    cout << "============> "<<tl->GetTextSize() << ",  "  << tl->GetY()<< endl;
	    if(system=="pp"){
	      tl->SetX(tl->GetX()+0.23);
	      tl->SetTextSize(0.033);
	    }
	    else if(system=="pPb"){
	      tl->SetX(tl->GetX()+0.25);
	      tl->SetY(tl->GetY()+0.04);
	      tl->SetTextSize(0.037);
	    }
	    cout << "----------> " << strTitle41.Data() << endl;  
	  }
          tl->Draw();
	}
        
	if(strName.Contains("TGraph")){
	  //Printf("There is a TGraph");
	  TGraphAsymmErrors *gr=(TGraphAsymmErrors*)obj;
	  //   gr->SetMarkerStyle();
	  //   gr->SetMarkerColor(0);
	  if(system=="pp")gr->SetMarkerSize(0.7);
	  else if(system=="pPb")gr->SetMarkerSize(0.9);
	  gr->SetLineWidth(2);
	  gr->Draw("E2");
          
	  TGraphAsymmErrors *gr2=(TGraphAsymmErrors*)gr->Clone("grHelp");
	  //gr2->SetMarkerColor(kBlack);
	  //gr2->SetLineColor(kBlack);
	  //gr2->SetMarkerSize(0.3);
	  //gr2->Draw("p");
	}
	//else  obj->Draw("same"); 
      }      
    }
    
    if(k!=0 || k!=3 || k!=6){

      hcur->Draw("sameep");      
      for(Int_t jl=0;jl<entries;jl++){
	
	TObject *obj=lc->At(jl);
	TString strName=obj->ClassName();
        
	if(strName.Contains("TGraph")){
	  //Printf("There is a TGraph");
	  TGraphAsymmErrors *gr=(TGraphAsymmErrors*)obj;
	  //   gr->SetMarkerStyle();
	  //   gr->SetMarkerColor(0);
	  if(system=="pp")gr->SetMarkerSize(0.7);
	  else if(system=="pPb")gr->SetMarkerSize(0.9);
	  gr->SetLineWidth(2);
	  gr->Draw("E2");
       
	  TGraphAsymmErrors *gr2=(TGraphAsymmErrors*)gr->Clone("grHelp");
	  //gr2->SetMarkerColor(kBlack);
	  //gr2->SetLineColor(kBlack);
	  //gr2->SetMarkerSize(0.3);
	  //gr2->Draw("p");
	}
	//else  obj->Draw("same");
      }
    }
   
    CanvasCorrelation->Update();
  }
    
    // saving the canvases in .root and .png
    TString ptoutput="", direcname="";
    direcname += "Output_SngCav_Comparison";
    ptoutput += "Comparison_DHCorrelations_assopT";
    ptoutput += TrackpTthr.Data();
    if(system == "pp") ptoutput += "_pp";
    if(system == "pPb") ptoutput += "_pPb";
    
    SaveCanvas(CanvasCorrelation, direcname, ptoutput);
    Printf("...  .. .. Done !");
    
  return;  
}


//_______________________________________________________________________
void LoadFileNamesToCompare(TString TrackpT ="0.3to1.0",TString Sys = "pPb"){
  
  Int_t skip3to5pPb=0;
  if(Sys.EqualTo("pp")){
    nhistos=9;
  }
  else if(Sys.EqualTo("pPb")){
    nhistos=6; 
    skip3to5pPb=3;
  }
  filenames=new TString[nhistos];
  Printf(" -------------------> Adding %s in %s coll",TrackpT.Data(),Sys.Data());
  if(Sys.EqualTo("pp")){
    filenames[0] = Form("%s/CanvaAndVariedHisto%sDzeroPt3to5assocPt%s.root",inputdirectory.Data(),Sys.Data(),TrackpT.Data());
    filenames[1] = Form("%s/CanvaAndVariedHisto%sDplusPt3to5assocPt%s.root",inputdirectory.Data(),Sys.Data(),TrackpT.Data());
    filenames[2] = Form("%s/CanvaAndVariedHisto%sDstarPt3to5assocPt%s.root",inputdirectory.Data(),Sys.Data(),TrackpT.Data());
  }
  filenames[3-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%sDzeroPt5to8assocPt%s.root",inputdirectory.Data(),Sys.Data(),TrackpT.Data());
  filenames[4-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%sDplusPt5to8assocPt%s.root",inputdirectory.Data(),Sys.Data(),TrackpT.Data());
  filenames[5-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%sDstarPt5to8assocPt%s.root",inputdirectory.Data(),Sys.Data(),TrackpT.Data());
  filenames[6-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%sDzeroPt8to16assocPt%s.root",inputdirectory.Data(),Sys.Data(),TrackpT.Data());
  filenames[7-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%sDplusPt8to16assocPt%s.root",inputdirectory.Data(),Sys.Data(),TrackpT.Data());
  filenames[8-skip3to5pPb] = Form("%s/CanvaAndVariedHisto%sDstarPt8to16assocPt%s.root",inputdirectory.Data(),Sys.Data(),TrackpT.Data());
  
}  

  

//_______________________________________________________________________
TCanvas * GetCanvas(Int_t i, TString canvasname = "cDraw"){
  
  //load the histogram with
  TString path = filenames[i];
  //cout <<"file #" <<i<<": Reading File from path --->> " << path << endl;
  TFile * file = TFile::Open(path.Data(),"READ");
  TCanvas * c = (TCanvas*)file->Get(canvasname.Data());
  return c;
}


//_______________________________________________________________________
void SaveCanvas(TCanvas * c, TString directory, TString name){
  
  TString outputDir = "";//Plots/15_May/";
  
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
