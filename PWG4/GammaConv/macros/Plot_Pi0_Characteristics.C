#include <fstream>
#include <Riostream.h>
/*
 *
 *
 *
 */

extern TRandom *gRandom;
extern TBenchmark *gBenchmark;
extern TSystem *gSystem;

void Plot_Pi0_Characteristics(const char *inputRootFile = "Pi0Characteristics",const char *path = "./Output/"){	

  gROOT->Reset();	
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);
  TString filename = Form("%s%s.root",path,inputRootFile);	
  TFile f(filename.Data());  
  
  TList *histograms = f.GetListOfKeys(); // get the list of directories in the file
  
  TString cutSelectionArray[15];
  Int_t cutsAdded=0;
  //Per pt bin histograms:
  TH1F* YieldPerBin[32];
  TH1F* StoBPerBin[32];
  TH1F* SignificancePerBin[32];
  TH1F* MassPerBin[32];
  TH1F* FWHMPerBin[32];

  Float_t lowBinLimits[32];
  Float_t highBinLimits[32];

  TH1F * Raw_Yield;
  TH1F * StoB;
  TH1F * Significance;
  TH1F * Mass;
  TH1F * FWHM;



  Int_t colorCounterRawYield=1;
  Int_t histogramCounterRawYield=2;
  TCanvas *canvasRawYield = new  TCanvas("canvasRawYield","",200,10,600,600);
  canvasRawYield->SetFillColor(0);
  canvasRawYield->SetLogy();
  TLegend *legRawYield= new TLegend(0.7,0.7,1.,1.);
  legRawYield->SetFillColor(0);

  Int_t colorCounterStoB=1;
  Int_t histogramCounterStoB=2;
  TCanvas *canvasStoB = new  TCanvas("canvasStoB","",200,10,600,600);
  canvasStoB->SetFillColor(0);
  canvasStoB->SetLogy();
  TLegend *legStoB= new TLegend(0.7,0.7,1.,1.);
  legStoB->SetFillColor(0);

  Int_t colorCounterSignificance=1;
  Int_t histogramCounterSignificance=2;
  TCanvas *canvasSignificance = new  TCanvas("canvasSignificance","",200,10,600,600);
  canvasSignificance->SetFillColor(0);
  canvasSignificance->SetLogy();
  TLegend *legSignificance= new TLegend(0.7,0.7,1.,1.);
  legSignificance->SetFillColor(0);

  Int_t colorCounterMass=1;
  Int_t histogramCounterMass=2;
  TCanvas *canvasMass = new  TCanvas("canvasMass","",200,10,600,600);
  canvasMass->SetFillColor(0);
  //  canvasMass->SetLogy();
  TLegend *legMass= new TLegend(0.7,0.7,1.,1.);
  legMass->SetFillColor(0);

  Int_t colorCounterFWHM=1;
  Int_t histogramCounterFWHM=2;
  TCanvas *canvasFWHM = new  TCanvas("canvasFWHM","",200,10,600,600);
  canvasFWHM->SetFillColor(0);
  //  canvasFWHM->SetLogy();
  TLegend *legFWHM= new TLegend(0.7,0.7,1.,1.);
  legFWHM->SetFillColor(0);


  for(Int_t entFile=0;entFile<histograms->GetEntries();entFile++){

    TString histoname = histograms->At(entFile)->GetName();
    if(histoname.Contains("Raw_Yield")){
      canvasRawYield->cd();
      cout<<"Histogram contains Raw_Yield: "<<histoname.Data()<<endl;
      if(Raw_Yield==NULL){
	cout<<"Setting first histogram"<<endl;
	Raw_Yield = (TH1F*)f.Get(histoname.Data());
	cout<<"Histoname: "<<histoname.Data()<<endl;
	TString cutValue= histoname(histoname.Index("Pi0_")+4,histoname.Length());
	cutSelectionArray[cutsAdded]= cutValue;
	cutsAdded++;
	cout<<"Cut Value is: "<<cutValue.Data()<<endl;
	for(Int_t bin=1;bin<Raw_Yield->GetNbinsX();bin++){
	  TString perBinHistoName = Form("Raw Yield %f <pt <%f",Raw_Yield->GetBinLowEdge(bin),Raw_Yield->GetBinLowEdge(bin)+Raw_Yield->GetBinWidth(bin));
	  lowBinLimits[bin]=Raw_Yield->GetBinLowEdge(bin);
	  highBinLimits[bin]=Raw_Yield->GetBinLowEdge(bin)+Raw_Yield->GetBinWidth(bin);
	  YieldPerBin[bin]= new TH1F(perBinHistoName.Data(),perBinHistoName.Data(),32,0,32);
	  YieldPerBin[bin]->SetBinContent(histogramCounterRawYield,Raw_Yield->GetBinContent(bin));
	  YieldPerBin[bin]->SetBinError(histogramCounterRawYield,Raw_Yield->GetBinError(bin));
	  YieldPerBin[bin]->GetXaxis()->SetBinLabel(histogramCounterRawYield,cutValue.Data());
	}
	cout<<"Raw_Yield: "<<Raw_Yield<<endl;
	legRawYield->AddEntry(Raw_Yield,cutValue.Data(),"l");	     
	//       	Raw_Yield->DrawCopy();
       	Raw_Yield->Draw();
	canvasRawYield->Update();
      }
      else{
	colorCounterRawYield++;
	TH1F* histogram = (TH1F*)f.Get(histoname.Data());
	cout<<"histogram number of bins: "<<histogram->GetNbinsX()<<endl;
	TString cutValue= histoname(histoname.Index("Pi0_")+4,histoname.Length());
	cutSelectionArray[cutsAdded]= cutValue;
	cutsAdded++;
	for(Int_t bin=1;bin<histogram->GetNbinsX();bin++){
	  //	   YieldPerBin[bin]= new TH1F(Form("yield bin %d",bin),Form("yield bin %d",bin),32,0,32);
	  YieldPerBin[bin]->SetBinContent(histogramCounterRawYield,histogram->GetBinContent(bin));
	  YieldPerBin[bin]->SetBinError(histogramCounterRawYield,histogram->GetBinError(bin));
	  YieldPerBin[bin]->GetXaxis()->SetBinLabel(histogramCounterRawYield,cutValue.Data());
	}

	histogram->SetLineColor(colorCounterRawYield);
	legRawYield->AddEntry(histogram,cutValue.Data(),"l");
	legRawYield->Draw();
	//	histogram->DrawCopy("same");
	histogram->Draw("same");
	canvasRawYield->Update();
      }
      histogramCounterRawYield+=2;
      
    }

    if(histoname.Contains("SB_Pi0")){
      canvasStoB->cd();
      cout<<"Histogram contains SB: "<<histoname.Data()<<endl;
      if(StoB==NULL){
	cout<<"Setting first histogram"<<endl;
	StoB = (TH1F*)f.Get(histoname.Data());
	cout<<"Histoname: "<<histoname.Data()<<endl;
	TString cutValue= histoname(histoname.Index("Pi0_")+4,histoname.Length());
	cout<<"Cut Value is: "<<cutValue.Data()<<endl;
	for(Int_t bin=1;bin<Raw_Yield->GetNbinsX();bin++){
	  TString perBinHistoName = Form("S/B %f <pt <%f",StoB->GetBinLowEdge(bin),StoB->GetBinLowEdge(bin)+StoB->GetBinWidth(bin));
	  StoBPerBin[bin]= new TH1F(perBinHistoName.Data(),perBinHistoName.Data(),32,0,32);
	  StoBPerBin[bin]->SetBinContent(histogramCounterStoB,StoB->GetBinContent(bin));
	  StoBPerBin[bin]->SetBinError(histogramCounterStoB,StoB->GetBinError(bin));
	  StoBPerBin[bin]->GetXaxis()->SetBinLabel(histogramCounterStoB,cutValue.Data());
	}
	legStoB->AddEntry(StoB,cutValue.Data(),"l");	     
	StoB->DrawCopy();
	canvasStoB->Update();
      }
      else{
	colorCounterStoB++;
	TH1F* histogram = (TH1F*)f.Get(histoname.Data());
	cout<<"histogram number of bins: "<<histogram->GetNbinsX()<<endl;
	TString cutValue= histoname(histoname.Index("Pi0_")+4,histoname.Length());
	for(Int_t bin=1;bin<histogram->GetNbinsX();bin++){
	  //	   YieldPerBin[bin]= new TH1F(Form("yield bin %d",bin),Form("yield bin %d",bin),32,0,32);
	  StoBPerBin[bin]->SetBinContent(histogramCounterStoB,histogram->GetBinContent(bin));
	  StoBPerBin[bin]->SetBinError(histogramCounterStoB,histogram->GetBinError(bin));
	  StoBPerBin[bin]->GetXaxis()->SetBinLabel(histogramCounterStoB,cutValue.Data());
	}

	histogram->SetLineColor(colorCounterStoB);
	legStoB->AddEntry(histogram,cutValue.Data(),"l");
	legStoB->Draw();
	histogram->DrawCopy("same");
	canvasStoB->Update();
      }
      histogramCounterStoB+=2;
    }

    if(histoname.Contains("Significance_Pi0")){
      canvasSignificance->cd();
      cout<<"Histogram contains Significance: "<<histoname.Data()<<endl;
      if(Significance==NULL){
	cout<<"Setting first histogram"<<endl;
	Significance = (TH1F*)f.Get(histoname.Data());
	cout<<"Histoname: "<<histoname.Data()<<endl;
	TString cutValue= histoname(histoname.Index("Pi0_")+4,histoname.Length());
	cout<<"Cut Value is: "<<cutValue.Data()<<endl;
	for(Int_t bin=1;bin<Raw_Yield->GetNbinsX();bin++){
	  TString perBinHistoName = Form("Significance %f <pt <%f",Significance->GetBinLowEdge(bin),Significance->GetBinLowEdge(bin)+Significance->GetBinWidth(bin));
	  SignificancePerBin[bin]= new TH1F(perBinHistoName.Data(),perBinHistoName.Data(),32,0,32);
	  SignificancePerBin[bin]->SetBinContent(histogramCounterSignificance,Significance->GetBinContent(bin));
	  SignificancePerBin[bin]->SetBinError(histogramCounterSignificance,Significance->GetBinError(bin));
	  SignificancePerBin[bin]->GetXaxis()->SetBinLabel(histogramCounterSignificance,cutValue.Data());
	}
	legSignificance->AddEntry(Significance,cutValue.Data(),"l");	     
	Significance->DrawCopy();
	canvasSignificance->Update();
	//	Raw_Yield = new TH1F("Raw_Yield","Raw_Yield",histogram->GetN)
      }
      else{
	colorCounterSignificance++;
	TH1F* histogram = (TH1F*)f.Get(histoname.Data());
	cout<<"histogram number of bins: "<<histogram->GetNbinsX()<<endl;
	TString cutValue= histoname(histoname.Index("Pi0_")+4,histoname.Length());
	for(Int_t bin=1;bin<histogram->GetNbinsX();bin++){
	  //	   YieldPerBin[bin]= new TH1F(Form("yield bin %d",bin),Form("yield bin %d",bin),32,0,32);
	  SignificancePerBin[bin]->SetBinContent(histogramCounterSignificance,histogram->GetBinContent(bin));
	  SignificancePerBin[bin]->SetBinError(histogramCounterSignificance,histogram->GetBinError(bin));
	  SignificancePerBin[bin]->GetXaxis()->SetBinLabel(histogramCounterSignificance,cutValue.Data());
	}

	histogram->SetLineColor(colorCounterSignificance);
	legSignificance->AddEntry(histogram,cutValue.Data(),"l");
	legSignificance->Draw();
	histogram->DrawCopy("same");
	canvasSignificance->Update();
      }
      histogramCounterSignificance+=2;
    }

    if(histoname.Contains("Mass_Pi0")){
      canvasMass->cd();
      cout<<"Histogram contains Mass: "<<histoname.Data()<<endl;
      if(Mass==NULL){
	cout<<"Setting first histogram"<<endl;
	Mass = (TH1F*)f.Get(histoname.Data());
	Mass->SetMinimum(0.1);
	Mass->SetMaximum(.16);
	cout<<"Histoname: "<<histoname.Data()<<endl;
	TString cutValue= histoname(histoname.Index("Pi0_")+4,histoname.Length());
	cout<<"Cut Value is: "<<cutValue.Data()<<endl;
	for(Int_t bin=1;bin<Raw_Yield->GetNbinsX();bin++){
	  TString perBinHistoName = Form("Mass %f <pt <%f",Mass->GetBinLowEdge(bin),Mass->GetBinLowEdge(bin)+Mass->GetBinWidth(bin));
	  MassPerBin[bin]= new TH1F(perBinHistoName.Data(),perBinHistoName.Data(),32,0,32);
	  MassPerBin[bin]->SetBinContent(histogramCounterMass,Mass->GetBinContent(bin));
	  MassPerBin[bin]->SetBinError(histogramCounterMass,Mass->GetBinError(bin));
	  MassPerBin[bin]->GetXaxis()->SetBinLabel(histogramCounterMass,cutValue.Data());
	}
	legMass->AddEntry(Mass,cutValue.Data(),"l");	     
	Mass->DrawCopy();
	canvasMass->Update();
	//	Raw_Yield = new TH1F("Raw_Yield","Raw_Yield",histogram->GetN)
      }
      else{
	colorCounterMass++;
	TH1F* histogram = (TH1F*)f.Get(histoname.Data());
	cout<<"histogram number of bins: "<<histogram->GetNbinsX()<<endl;
	TString cutValue= histoname(histoname.Index("Pi0_")+4,histoname.Length());
	for(Int_t bin=1;bin<histogram->GetNbinsX();bin++){
	  //	   YieldPerBin[bin]= new TH1F(Form("yield bin %d",bin),Form("yield bin %d",bin),32,0,32);
	  MassPerBin[bin]->SetBinContent(histogramCounterMass,histogram->GetBinContent(bin));
	  MassPerBin[bin]->SetBinError(histogramCounterMass,histogram->GetBinError(bin));
	  MassPerBin[bin]->GetXaxis()->SetBinLabel(histogramCounterMass,cutValue.Data());
	}

	histogram->SetLineColor(colorCounterMass);
	legMass->AddEntry(histogram,cutValue.Data(),"l");
	legMass->Draw();
	histogram->DrawCopy("same");
	canvasMass->Update();
      }
      histogramCounterMass+=2;
    }

    if(histoname.Contains("FWHM_Pi0")){
      canvasFWHM->cd();
      cout<<"Histogram contains FWHM: "<<histoname.Data()<<endl;
      if(FWHM==NULL){
	cout<<"Setting first histogram"<<endl;
	FWHM = (TH1F*)f.Get(histoname.Data());
	cout<<"Histoname: "<<histoname.Data()<<endl;
	TString cutValue= histoname(histoname.Index("Pi0_")+4,histoname.Length());
	cout<<"Cut Value is: "<<cutValue.Data()<<endl;
	for(Int_t bin=1;bin<Raw_Yield->GetNbinsX();bin++){
	  TString perBinHistoName = Form("FWHM %f <pt <%f",FWHM->GetBinLowEdge(bin),FWHM->GetBinLowEdge(bin)+FWHM->GetBinWidth(bin));
	  FWHMPerBin[bin]= new TH1F(perBinHistoName.Data(),perBinHistoName.Data(),32,0,32);
	  FWHMPerBin[bin]->SetBinContent(histogramCounterFWHM,FWHM->GetBinContent(bin));
	  FWHMPerBin[bin]->SetBinError(histogramCounterFWHM,FWHM->GetBinError(bin));
	  FWHMPerBin[bin]->GetXaxis()->SetBinLabel(histogramCounterFWHM,cutValue.Data());
	}
	legFWHM->AddEntry(FWHM,cutValue.Data(),"l");	     
	FWHM->DrawCopy();
	canvasFWHM->Update();
	//	Raw_Yield = new TH1F("Raw_Yield","Raw_Yield",histogram->GetN)
      }
      else{
	colorCounterFWHM++;
	TH1F* histogram = (TH1F*)f.Get(histoname.Data());
	cout<<"histogram number of bins: "<<histogram->GetNbinsX()<<endl;
	TString cutValue= histoname(histoname.Index("Pi0_")+4,histoname.Length());
	for(Int_t bin=1;bin<histogram->GetNbinsX();bin++){
	  //	   YieldPerBin[bin]= new TH1F(Form("yield bin %d",bin),Form("yield bin %d",bin),32,0,32);
	  FWHMPerBin[bin]->SetBinContent(histogramCounterFWHM,histogram->GetBinContent(bin));
	  FWHMPerBin[bin]->SetBinError(histogramCounterFWHM,histogram->GetBinError(bin));
	  FWHMPerBin[bin]->GetXaxis()->SetBinLabel(histogramCounterFWHM,cutValue.Data());
	}

	histogram->SetLineColor(colorCounterFWHM);
	legFWHM->AddEntry(histogram,cutValue.Data(),"l");
	legFWHM->Draw();
	histogram->DrawCopy("same");
	canvasFWHM->Update();
      }
      histogramCounterFWHM+=2;
    }

  }//end of for loop over histograms 



  TPostScript *ps_characteristics;
  ps_characteristics = new TPostScript(Form("%sPi0Characteristics.ps",path),111);	
  ps_characteristics->NewPage();

 
  canvasRawYield->Update();
  canvasRawYield->Draw();
  ps_characteristics->NewPage();
  canvasSignificance->Update();
  canvasSignificance->Draw();
  ps_characteristics->NewPage();
  canvasStoB->Update();
  canvasStoB->Draw();
  ps_characteristics->NewPage();
  canvasMass->Update();
  canvasMass->Draw();
  ps_characteristics->NewPage();
  canvasFWHM->Update();
  canvasFWHM->Draw();
  ps_characteristics->NewPage();



 
  /*  
  
  TCanvas *canvas1 = new TCanvas("Integrated1","",10,10,700,1000);
  TPad *pad1 = new TPad("pad1","",0.,0.,1.,1.,0);
  pad1->SetFillColor(0);
  pad1->GetFrame()->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->Divide(1,2);
  pad1->Draw();
  
  pad1->cd(1);
  Raw_Yield->Draw();
  
  pad1->cd(2);
  StoB->Draw();
  ps_characteristics->NewPage();
  
  
  TCanvas *canvas2 = new TCanvas("Integrated2","",10,10,700,1000);
  TPad *pad2 = new TPad("pad2","",0.,0.,1.,1.,0);
  pad2->SetFillColor(0);
  pad2->GetFrame()->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->Divide(1,2);
  pad2->Draw();
  
  pad2->cd(1);
  Significance->Draw();
  
  pad2->cd(2);
  Mass->Draw();
  
  ps_characteristics->NewPage();   
  
  TCanvas *canvas3 = new TCanvas("Integrated3","",10,10,700,1000);
  TPad *pad3 = new TPad("pad3","",0.,0.,1.,1.,0);
  pad3->SetFillColor(0);
  pad3->GetFrame()->SetFillColor(0);
  pad3->SetBorderMode(0);
  pad3->Divide(1,2);
  pad3->Draw();
  
  pad3->cd(1);
  FWHM->Draw();
  
  pad3->cd(2);
  //    Mass->Draw();
  */
  ps_characteristics->NewPage();
  
  for(Int_t bin=2;bin<Raw_Yield->GetNbinsX();bin++){
    TString canvasname= Form("bin %d",bin);
    TCanvas *binC = new TCanvas(canvasname.Data(),"",10,10,700,1000);
    TString padname= Form("pad %d",bin);
    //    TPad *pad = new TPad(padname.Data(),"",0.05,0.05,0.95,0.95,0);
    TPad *pad = new TPad(padname.Data(),"",0.,0.,1.,1.,0);
    pad->SetFillColor(0);
    pad->GetFrame()->SetFillColor(0);
    pad->SetBorderMode(0);
    pad->Divide(1,5);
    pad->Draw();

    pad->cd(1);
    YieldPerBin[bin]->Draw();

    pad->cd(2);
    StoBPerBin[bin]->Draw();
    
    pad->cd(3);
    SignificancePerBin[bin]->Draw();
 
    pad->cd(4);
    MassPerBin[bin]->Draw();

    pad->cd(5);
    FWHMPerBin[bin]->Draw();

    binC->Update();
    binC->Close();

    ps_characteristics->NewPage();
    binC->Close();
  }

  Int_t rebinValue=4;

  for(Int_t cuts=0;cuts<cutsAdded;cuts++){
    cout<<"CUT: "<<cutSelectionArray[cuts].Data()<<endl;
    //    for(Int_t bin=2;bin<15;bin++){
      TCanvas *canvasTest = new  TCanvas("canvastest","",200,10,600,600);
      TPad *pad = new TPad(padname.Data(),"",0.,0.,1.,1.,0);
      pad->SetFillColor(0);
      pad->GetFrame()->SetFillColor(0);
      pad->SetBorderMode(0);
      pad->Divide(4,4);
      pad->Draw();
    for(Int_t bin=2;bin<18;bin++){
      pad->cd(bin-1);
      TString namet= Form("Mapping_Reco_InvMass_in_Pt_Bin%s%02d",cutSelectionArray[cuts].Data(),bin);
      cout<<"Getting histogram: "<<namet.Data()<<endl;
      TH1F * signalt = (TH1F*)f.Get(namet.Data());
      signalt->Rebin(rebinValue);
      TString titlet= Form("Inv_Mass_cut%s_pt[%f,%f]",cutSelectionArray[cuts].Data(),lowBinLimits[bin],highBinLimits[bin]);
      signalt->SetTitle(titlet.Data());
      signalt->Sumw2();
      signalt->SetAxisRange(0.,0.7);
      signalt->Draw();

      TString nameb= Form("Mapping_Back_InvMass_in_Pt_Bin%s%02d",cutSelectionArray[cuts].Data(),bin);
      cout<<"Getting histogram: "<<nameb.Data()<<endl;
      TH1F * signalb = (TH1F*)f.Get(nameb.Data());
      signalb->Rebin(rebinValue);
      TString titleb= Form("Inv_Mass_cut%s_pt[%f,%f]",cutSelectionArray[cuts].Data(),lowBinLimits[bin],highBinLimits[bin]);
      signalb->SetTitle(titleb.Data());
      signalb->SetAxisRange(0.,0.7);
      signalb->SetLineColor(4);
      signalb->Draw("same");
      canvasTest->Update();

      //      ps_characteristics->NewPage();
    }
    ps_characteristics->NewPage();

    TCanvas *canvasTestDiff = new  TCanvas("canvastestdiff","",200,10,600,600);
    TPad *padD = new TPad(padname.Data(),"",0.,0.,1.,1.,0);
    canvasTestDiff->SetFillColor(0);
    padD->SetFillColor(0);
    padD->GetFrame()->SetFillColor(0);
    padD->SetBorderMode(0);
    padD->Divide(4,4);
    padD->Draw();
    for(Int_t bin=2;bin<18;bin++){
      padD->cd(bin-1);
      TString name= Form("Mapping_Signal_InvMass_in_Pt_Bin%s%02d",cutSelectionArray[cuts].Data(),bin);
      cout<<"Getting histogram: "<<name.Data()<<endl;
      TH1F * signal = (TH1F*)f.Get(name.Data());
      signal->Rebin(rebinValue);
      TString title= Form("Signal_Inv_Mass_cut%s_pt[%f,%f]",cutSelectionArray[cuts].Data(),lowBinLimits[bin],highBinLimits[bin]);
      signal->SetTitle(title.Data());
      signal->SetAxisRange(0.,0.7);
      signal->Draw();
      
      canvasTestDiff->Update();

      //      ps_characteristics->NewPage();
    }
    ps_characteristics->NewPage();
  }

  TCanvas *canvasTest = new  TCanvas("canvastest","",200,10,600,600);
  TPad *pad = new TPad(padname.Data(),"",0.,0.,1.,1.,0);
  pad->SetFillColor(0);
  pad->GetFrame()->SetFillColor(0);
  pad->SetBorderMode(0);
  pad->Divide(3,3);
  pad->Draw();
  for(Int_t bin=0;bin<cutsAdded;bin++){
    cout<<"CUT: "<<cutSelectionArray[bin].Data()<<endl;
    pad->cd(bin+1);
    pad->cd(bin+1)->SetLogz(1);
    TString namet= Form("ESD_Mother_InvMass_%s",cutSelectionArray[bin].Data());
    cout<<"Getting histogram: "<<namet.Data()<<endl;
    TH1F * massAll = (TH1F*)f.Get(namet.Data());
    
    TString titlet= Form("CutId%s",cutSelectionArray[bin].Data());
    massAll->SetTitle(titlet.Data());
    massAll->Draw();
    canvasTest->Update();
    canvasTest->Print("massAll.gif");
  }
  ps_characteristics->NewPage();

    //    for(Int_t bin=2;bin<15;bin++){
  TCanvas *canvasTest = new  TCanvas("canvastest","",200,10,600,600);
  TPad *pad = new TPad(padname.Data(),"",0.,0.,1.,1.,0);
  pad->SetFillColor(0);
  pad->GetFrame()->SetFillColor(0);
  pad->SetBorderMode(0);
  pad->Divide(3,3);
  pad->Draw();
  for(Int_t bin=0;bin<cutsAdded;bin++){
    cout<<"CUT: "<<cutSelectionArray[bin].Data()<<endl;
    pad->cd(bin+1);
    pad->cd(bin+1)->SetLogx(1);
    pad->cd(bin+1)->SetLogz(1);
    TString namet= Form("ESD_ConvGamma_E_dEdxP_%s",cutSelectionArray[bin].Data());
    cout<<"Getting histogram: "<<namet.Data()<<endl;
    TH1F * dedxp = (TH1F*)f.Get(namet.Data());
    
    TString titlet= Form("CutId%s",cutSelectionArray[bin].Data());
    dedxp->SetTitle(titlet.Data());
    dedxp->Draw("col2");
    canvasTest->Update();
    canvasTest->Print("dedxp.gif");
  }
  ps_characteristics->NewPage();
  TCanvas *canvasTest = new  TCanvas("canvastest","",200,10,600,600);
  TPad *pad = new TPad(padname.Data(),"",0.,0.,1.,1.,0);
  pad->SetFillColor(0);
  pad->GetFrame()->SetFillColor(0);
  pad->SetBorderMode(0);
  pad->Divide(3,3);
  pad->Draw();
  for(Int_t bin=0;bin<cutsAdded;bin++){
    cout<<"CUT: "<<cutSelectionArray[bin].Data()<<endl;
    pad->cd(bin+1);
    pad->cd(bin+1)->SetLogz(1);
    TString namet= Form("ESD_ConvGamma_alfa_qt_%s",cutSelectionArray[bin].Data());
    cout<<"Getting histogram: "<<namet.Data()<<endl;
    TH1F * armen = (TH1F*)f.Get(namet.Data());
    
    TString titlet= Form("CutId%s",cutSelectionArray[bin].Data());
    armen->SetTitle(titlet.Data());
    armen->Draw("col2");
    canvasTest->Update();
    canvasTest->Print("armen.gif");
  }

  ps_characteristics->Close();
}
