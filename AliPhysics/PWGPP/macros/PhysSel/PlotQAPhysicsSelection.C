//using namespace std;

#include <iostream>

#include "QALHC11c_QA69.h"
#include "DefinePlots.h"
#include "TGraphErrors.h"


const char * GetLocalFileName(Int_t run, const char * suffix, const char * path);
TH1F * GetEmpty(const char * name, Int_t nfile);
TGraphErrors * GetGraphRej(TGraphErrors * gr, TList * rejRunList, const char * reason, Float_t &mean, Bool_t doDraw) ;
double GetMedian(double* arr, int n);
double meanMed(double* vec, int np, double nsigmaCut, int &nrej, int *rejList);
TGraphErrors * graph[kNGraphs];

//---------------------------------------------------------
TH1F *fAccOverAll[3];
TH1F *fBGOverAll[3];
TH1F *fV0BGOverAll[3];
TH1F *fV0ABGOverAll[3];
TH1F *fV0CBGOverAll[3];
TH1F *fAfterOverBefore[3];
TH1F *fV0AOverV0C[3];
TH1F *fF02OverAll[3];
TH1F *fF01OverAll[3];

void PlotQAPhysicsSelection() {


  for(Int_t i=0;i<3;++i){
    fAccOverAll[i]=new TH1F(Form("fAccOverAll_%d",i),Form("fAccOverAll_%d",i),100,0,1);
    fBGOverAll[i]=new TH1F(Form("fBGOverAll_%d",i),Form("fBGOverAll_%d",i),50,0,0.5);
    fV0BGOverAll[i]=new TH1F(Form("fV0BGOverAll_%d",i),Form("fV0BGOverAll_%d",i),50,0,0.1);
    fV0ABGOverAll[i]=new TH1F(Form("fV0ABGOverAll_%d",i),Form("fV0ABGOverAll_%d",i),50,0,0.1);
    fV0CBGOverAll[i]=new TH1F(Form("fV0CBGOverAll_%d",i),Form("fV0CBGOverAll_%d",i),50,0,0.1);
    fAfterOverBefore[i]=new TH1F(Form("fAfterOverBefore_%d",i),Form("fAfterOverBefore_%d",i),100,0.5,1);
    fV0AOverV0C[i]=new TH1F(Form("fV0AOverV0C_%d",i),Form("fV0AOverV0C_%d",i),100,0.7,1.2);
    fF01OverAll[i]=new TH1F(Form("fF01OverAll_%d",i),Form("fF01OverAll_%d",i),200,0.6,1.1);
    fF02OverAll[i]=new TH1F(Form("fF02OverAll_%d",i),Form("fF02OverAll_%d",i),200,0.6,1.1);
  }


  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // Book graphs
  TGraphErrors * graph[kNGraphs] = {0};
  for(Int_t igraph = 0; igraph < kNGraphs; igraph++){
    graph[igraph] = new TGraphErrors;
    graph[igraph]->SetName(gnames[igraph]);
    graph[igraph]->GetXaxis()->SetTitle("run");
    graph[igraph]->GetYaxis()->SetTitle(ylabels[igraph]);
    graph[igraph]->SetMarkerStyle(20);
  }

  //loading libraries
  loadlibs();

  TList * listEmptyRuns    = new TList(); // keep track of empty runs
  TList * listRejectedRuns = new TList(); // keep track of rejected runes

  // Count ncycles
  Int_t cycle=0;
  while(QAcycle[++cycle]>0) {
    
    cout << "." ;
  }  
  cout << "  Ncycles " << cycle <<endl;
  cout << QAcycle[0] ;

  // connect to grid
  if(!localMode) TGrid::Connect("alien://");  

  // do not use scientific notation for run number
  TGaxis::SetMaxDigits(7)  ;




  
  // loop over all files
  Int_t ifile =-1;
  Int_t ifileGood = 0;
  Int_t ifileNotEmpty  = 0;
  while (runs[++ifile] > 0) {
    
    Long_t *id,*size,*flags,*mt;
    



    // check if QA file is available and open it
    // try different QA train outputs
    TString file ;
    TFile *fr=0;
    TFile *fc=0; // centrality, only in local mode for the time being
    for(Int_t c=0;c<cycle;c++){ 
      if(fr)break; 
      if(output.Contains("AOD"))
	file.Form("alien://%s/%09d/%s%3.3d/event_stat.root",location.Data(),runs[ifile],
		  output.Data(),QAcycle[c] );
      else
	file.Form("alien://%s/%09d/%s%d/event_stat.root",location.Data(),runs[ifile],
		  output.Data(),QAcycle[c] );
      if(localMode){ 
	if(!GetLocalFileName(runs[ifile], localSuffix, localPath))continue;
	file = GetLocalFileName(runs[ifile], localSuffix, localPath);
	if(gSystem->GetPathInfo(file,id,size,flags,mt)){Printf("not found");continue;}
      }
  
      Printf("\nBegin of reading: %s", file.Data());
      fr=TFile::Open(file);     

    }
 


    // set value to -1 by default:
    for(Int_t igraph = 0; igraph < kNGraphs; igraph++){
      graph[igraph]->SetPoint(ifile, runs[ifile], -1);
    }

    // If the file is not available, continue
    if(!fr){
      Printf("File %d is not available.\n",runs[ifile]);
      listEmptyRuns->Add(new TObjString(Form("File not Existing [%d]",runs[ifile])));
      continue;
    }


    // get stats and fill graphs
    TH2F * hStats = (TH2F*) fr->Get("fHistStatistics"); 
    if(!localMode) {
      gSystem->Exec(Form("alien_cp %s %s",file.Data(), GetLocalFileName(runs[ifile], localSuffix, localPath)));
      cout << Form("alien_cp %s %s",file.Data(), GetLocalFileName(runs[ifile], localSuffix, localPath)) <<endl;
    }

 

    Int_t rowC0SMH   = -1;
    Int_t rowCMBAC   = -1;
    Int_t rowCMBACS2 = -1;
    Int_t nbiny = hStats->GetNbinsY();
    for(Int_t ibiny = 1; ibiny <= nbiny; ibiny++){
      TString label = hStats->GetYaxis()->GetBinLabel(ibiny);
      if(label.Contains(trigger3))   rowC0SMH   = ibiny; 
      if(label.Contains(trigger1))   rowCMBAC   = ibiny; 
      if(label.Contains(trigger2)) rowCMBACS2 = ibiny; 
    }

    //Number of events in the selected trigger class
    Float_t C0SMH   = hStats->GetBinContent(AliPhysicsSelection::kStatTriggerClass,rowC0SMH);    
    Float_t CMBAC   = hStats->GetBinContent(AliPhysicsSelection::kStatTriggerClass,rowCMBAC);
    Float_t CMBACS2 = hStats->GetBinContent(AliPhysicsSelection::kStatTriggerClass,rowCMBACS2);


    //Number of events after physics selection
    Float_t C0SMH_APS   = hStats->GetBinContent(AliPhysicsSelection::kStatOffline,rowC0SMH);    
    Float_t CMBAC_APS   = hStats->GetBinContent(AliPhysicsSelection::kStatOffline,rowCMBAC);
    Float_t CMBACS2_APS = hStats->GetBinContent(AliPhysicsSelection::kStatOffline,rowCMBACS2);


    //Fraction of Events rejected as background because out of the correlation tracklets vs clusters
    Float_t C0SMHBG   = hStats->GetBinContent(AliPhysicsSelection::kStatBG,rowC0SMH);    
    Float_t CMBACBG   = hStats->GetBinContent(AliPhysicsSelection::kStatBG,rowCMBAC);
    Float_t CMBACS2BG = hStats->GetBinContent(AliPhysicsSelection::kStatBG,rowCMBACS2);

 


    //Events with a signal from V0A in the collision time window
    Float_t C0SMHV0A   = hStats->GetBinContent(AliPhysicsSelection::kStatV0A,rowC0SMH);    
    Float_t CMBACV0A   = hStats->GetBinContent(AliPhysicsSelection::kStatV0A,rowCMBAC);
    Float_t CMBACS2V0A = hStats->GetBinContent(AliPhysicsSelection::kStatV0A,rowCMBACS2);
    
    //Events with a signal from V0C in the collision time window
    Float_t C0SMHV0C   = hStats->GetBinContent(AliPhysicsSelection::kStatV0C,rowC0SMH);    
    Float_t CMBACV0C   = hStats->GetBinContent(AliPhysicsSelection::kStatV0C,rowCMBAC);
    Float_t CMBACS2V0C = hStats->GetBinContent(AliPhysicsSelection::kStatV0C,rowCMBACS2);

    if(C0SMHV0C>0)graph[kGraphV0AOverV0CSMH]->SetPoint(ifile,runs[ifile], C0SMHV0A/C0SMHV0C);
    if(CMBACV0C>0)graph[kGraphV0AOverV0CAC]->SetPoint(ifile,runs[ifile], CMBACV0A/CMBACV0C);
    if(CMBACS2V0C>0)graph[kGraphV0AOverV0CAC2]->SetPoint(ifile,runs[ifile], CMBACS2V0A/CMBACS2V0C);

 
    //Events flagged as BG from V0 side A
    Float_t C0SMHV0ABG   = hStats->GetBinContent(AliPhysicsSelection::kStatV0ABG,rowC0SMH);    
    Float_t CMBACV0ABG   = hStats->GetBinContent(AliPhysicsSelection::kStatV0ABG,rowCMBAC);
    Float_t CMBACS2V0ABG = hStats->GetBinContent(AliPhysicsSelection::kStatV0ABG,rowCMBACS2);
    
    //Events flagged as BG from V0 side C
    Float_t C0SMHV0CBG   = hStats->GetBinContent(AliPhysicsSelection::kStatV0CBG,rowC0SMH);    
    Float_t CMBACV0CBG   = hStats->GetBinContent(AliPhysicsSelection::kStatV0CBG,rowCMBAC);
    Float_t CMBACS2V0CBG = hStats->GetBinContent(AliPhysicsSelection::kStatV0CBG,rowCMBACS2);


    //Number of events with more than 1 chip hit in the pixels, computed offline
    Float_t C0SMHF01   = hStats->GetBinContent(AliPhysicsSelection::kStatFO1,rowC0SMH);    
    Float_t CMBACF01   = hStats->GetBinContent(AliPhysicsSelection::kStatFO1,rowCMBAC);
    Float_t CMBACS2F01 = hStats->GetBinContent(AliPhysicsSelection::kStatFO1,rowCMBACS2);


    //Number of events with more than 2 chip hit in the pixels, computed offline
    Float_t C0SMHF02   = hStats->GetBinContent(AliPhysicsSelection::kStatFO2,rowC0SMH);    
    Float_t CMBACF02   = hStats->GetBinContent(AliPhysicsSelection::kStatFO2,rowCMBAC);
    Float_t CMBACS2F02 = hStats->GetBinContent(AliPhysicsSelection::kStatFO2,rowCMBACS2);


    //Accepted by V0
    Float_t C0SMHV0   = hStats->GetBinContent(AliPhysicsSelection::kStatV0,rowC0SMH);    
    Float_t CMBACV0   = hStats->GetBinContent(AliPhysicsSelection::kStatV0,rowCMBAC);
    Float_t CMBACS2V0 = hStats->GetBinContent(AliPhysicsSelection::kStatV0,rowCMBACS2);



    //Events flagged as BG from V0 and the A or C side
    Float_t C0SMHV0BG   = hStats->GetBinContent(AliPhysicsSelection::kStatV0ABG,rowC0SMH)  +hStats->GetBinContent(AliPhysicsSelection::kStatV0CBG,rowC0SMH)  ;
    Float_t CMBACV0BG   = hStats->GetBinContent(AliPhysicsSelection::kStatV0ABG,rowCMBAC)  +hStats->GetBinContent(AliPhysicsSelection::kStatV0CBG,rowCMBAC)  ;
    Float_t CMBACS2V0BG = hStats->GetBinContent(AliPhysicsSelection::kStatV0ABG,rowCMBACS2)+hStats->GetBinContent(AliPhysicsSelection::kStatV0CBG,rowCMBACS2);


   //Events passing the ZDC time cut on the correlation between the sum and the difference of the timing (rejects the "debunched" events)
    
    Float_t C0SMHZDC   = hStats->GetBinContent(AliPhysicsSelection::kStatZDCTime,rowC0SMH);
    Float_t CMBACZDC   = hStats->GetBinContent(AliPhysicsSelection::kStatZDCTime,rowCMBAC);
    Float_t CMBACS2ZDC = hStats->GetBinContent(AliPhysicsSelection::kStatZDCTime,rowCMBACS2);
                                                                    

    //Accepted events                                               kStatAccepted
    Float_t C0SMHACC   = hStats->GetBinContent(AliPhysicsSelection::kStatAccepted,rowC0SMH);
    Float_t CMBACACC   = hStats->GetBinContent(AliPhysicsSelection::kStatAccepted,rowCMBAC);
    Float_t CMBACS2ACC = hStats->GetBinContent(AliPhysicsSelection::kStatAccepted,rowCMBACS2);


    cout<<"C0SMHBG="<<C0SMHBG<<"   C0SMHACC="<<C0SMHACC<<endl;





    cout << runs[ifile] << " " << CMBACS2 << " " << CMBACS2V0BG << " " << " " << CMBACS2ACC << endl;
    
  
   

   
    

    if(CMBAC>0) {
      graph[kGraphBGOverAllAC]->SetPoint(ifile,runs[ifile], CMBACBG / CMBAC);
      graph[kGraphV0BGOverAllAC]->SetPoint(ifile,runs[ifile], CMBACV0BG / CMBAC);
      graph[kGraphV0ABGOverAllAC]->SetPoint(ifile,runs[ifile], CMBACV0ABG / CMBAC);//
      graph[kGraphV0CBGOverAllAC]->SetPoint(ifile,runs[ifile], CMBACV0CBG / CMBAC);//
      graph[kGraphACCOverAllAC]    ->SetPoint(ifile,runs[ifile],CMBACACC/CMBAC);
      graph[kGraphNevACratioAfter]    ->SetPoint(ifile,runs[ifile],CMBAC_APS/CMBAC);
      //F0
      graph[kGraphFO1OverAllAC]    ->SetPoint(ifile,runs[ifile],CMBACF01/CMBAC);
      graph[kGraphFO2OverAllAC]    ->SetPoint(ifile,runs[ifile],CMBACF02/CMBAC);



    }
    if(CMBACS2>0){
      graph[kGraphBGOverAllAC2]->SetPoint(ifile,runs[ifile], CMBACS2BG / CMBACS2);
      graph[kGraphV0BGOverAllAC2]->SetPoint(ifile,runs[ifile], CMBACS2V0BG / CMBACS2);
      graph[kGraphV0ABGOverAllAC2]->SetPoint(ifile,runs[ifile], CMBACS2V0ABG / CMBACS2);//
      graph[kGraphV0CBGOverAllAC2]->SetPoint(ifile,runs[ifile], CMBACS2V0CBG / CMBACS2);//
      graph[kGraphACCOverAllACS2]->SetPoint(ifile,runs[ifile],CMBACS2ACC/CMBACS2);
      graph[kGraphNevAC2ratioAfter]->SetPoint(ifile,runs[ifile],CMBACS2_APS/CMBACS2);

      //FO
      graph[kGraphFO1OverAllAC2]    ->SetPoint(ifile,runs[ifile],CMBACS2F01/CMBACS2);
      graph[kGraphFO2OverAllAC2]    ->SetPoint(ifile,runs[ifile],CMBACS2F02/CMBACS2);



    }
    if(C0SMH>0){
      graph[kGraphBGOverAllSMH]->SetPoint(ifile,runs[ifile], C0SMHBG / C0SMH);
      graph[kGraphV0BGOverAllSMH]->SetPoint(ifile,runs[ifile], C0SMHV0BG / C0SMH);
      graph[kGraphV0ABGOverAllSMH]->SetPoint(ifile,runs[ifile],  C0SMHV0ABG/ C0SMH);
      graph[kGraphV0CBGOverAllSMH]->SetPoint(ifile,runs[ifile],  C0SMHV0CBG/ C0SMH);
      graph[kGraphACCOverAllSMH]   ->SetPoint(ifile,runs[ifile],C0SMHACC/C0SMH);
      graph[kGraphNevSMHratioAfter]   ->SetPoint(ifile,runs[ifile],C0SMH_APS/C0SMH);

      graph[kGraphFO1OverAllSMH]    ->SetPoint(ifile,runs[ifile], C0SMHF01/C0SMH);
      graph[kGraphFO2OverAllSMH]    ->SetPoint(ifile,runs[ifile], C0SMHF02/C0SMH);


    } 
    
    
    
    

    graph[kGraphNevSMH]->SetPoint(ifile,runs[ifile],C0SMHACC);
    graph[kGraphNevAC ]->SetPoint(ifile,runs[ifile],CMBACACC);
    graph[kGraphNevAC2]->SetPoint(ifile,runs[ifile],CMBACS2ACC);
    graph[kGraphNevSMHBefore]->SetPoint(ifile,runs[ifile],C0SMH);
    graph[kGraphNevACBefore ]->SetPoint(ifile,runs[ifile],CMBAC);
    graph[kGraphNevAC2Before]->SetPoint(ifile,runs[ifile],CMBACS2);
    graph[kGraphNevSMHAfter]->SetPoint(ifile,runs[ifile],C0SMH_APS);
    graph[kGraphNevACAfter ]->SetPoint(ifile,runs[ifile],CMBAC_APS);
    graph[kGraphNevAC2After]->SetPoint(ifile,runs[ifile],CMBACS2_APS);

 

    ifileGood++;
    if(C0SMH>0 || CMBACS2>0 || CMBAC>0) ifileNotEmpty++;
    else  listEmptyRuns->Add(new TObjString(Form("No events [%d]",runs[ifile])));
 

  }


 


  // Set bin labels with run number and save graphs
  // also prepare a table with number of events
  AliLatexTable tableEvts (7,"c|ccc|ccc");
  //tableEvts.InsertCustomRow("& \\multicolumn{3}{c}{Before Phys. Sel.} & \\multicolumn{3}{c}{After Phys. Sel.}\\\\");
  tableEvts.InsertCustomRow("Run & CINT1B & CEMC1B & CSH1B &  CINT1B & CEMC1B & CSH1B \\\\"); 
 



  fou = new TFile("qaLHC10apass2.root","recreate");
  fou->cd();

  for(Int_t igraph = 0; igraph < kNGraphs; igraph++){
    //    TAxis * ax = graph[igraph]->GetHistogram()->GetXaxis();
    graph[igraph]->GetXaxis()->SetTitle("run");
    graph[igraph]->GetYaxis()->SetTitle(ylabels[igraph]);

    graph[igraph]->GetXaxis()->SetLabelSize(0.06);
    graph[igraph]->GetXaxis()->SetTitleSize(0.07);
    graph[igraph]->GetXaxis()->SetTitleOffset(0.5);
    graph[igraph]->GetYaxis()->SetLabelSize(0.08);
    graph[igraph]->GetYaxis()->SetTitleSize(0.06);
    graph[igraph]->GetYaxis()->SetTitleOffset(0.7);
    
    
    graph[igraph]->GetYaxis()->CenterTitle(0);
    graph[igraph]->GetYaxis()->SetLabelFont(42);
    graph[igraph]->GetYaxis()->SetTitleFont(42);
    graph[igraph]->GetYaxis()->SetNoExponent(kTRUE);
    graph[igraph]->GetXaxis()->SetLabelFont(42);
    graph[igraph]->GetXaxis()->SetTitleFont(52);

    graph[igraph]->SetMarkerStyle(30);
    graph[igraph]->SetMarkerColor(4);


    graph[igraph]->SetMinimum(0);
    graph[igraph]->SetName(gnames[igraph]);
    graph[igraph]->SetTitle(ylabels[igraph]);
    //    graph[igraph]->SetMaximum(1);
    graph[igraph]->Write();
  }

 
 
  fou->Close();


  
  Int_t npoint = graph[kGraphNevAC2] -> GetN();
  for(Int_t ipoint = 0; ipoint < npoint; ipoint++){      
    tableEvts.SetNextCol(TMath::Nint(graph[kGraphNevAC2]->GetX()[ipoint]));
    tableEvts.SetNextCol(TMath::Nint(graph[kGraphNevACBefore]->GetY()[ipoint]));
    tableEvts.SetNextCol(TMath::Nint(graph[kGraphNevAC2Before]->GetY()[ipoint]));
    tableEvts.SetNextCol(TMath::Nint(graph[kGraphNevSMHBefore]->GetY()[ipoint]));
    tableEvts.SetNextCol(TMath::Nint(graph[kGraphNevAC]->GetY()[ipoint]));
    tableEvts.SetNextCol(TMath::Nint(graph[kGraphNevAC2]->GetY()[ipoint]));
    tableEvts.SetNextCol(TMath::Nint(graph[kGraphNevSMH]->GetY()[ipoint]));
    //tableEvts.SetNextCol(TMath::Nint(graph[kGraphNev90]->GetY()[ipoint]));
    tableEvts.InsertRow();
    //   ax->SetBinLabel(ax->FindBin(ipoint), Form("%d", runs[ipoint-1]));
  }    




  Float_t meanDummy;
  

  /*
    THREE TRIGGERS IN THE SAME CANVAS
    
  */
  c4a = new TCanvas("GraphACCOverAll","GraphACCOverAll",1000,700);
  

  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();


  c4a->cd();
  pad1->Draw();


  pad1->cd();


  graph[kGraphACCOverAllAC]->Draw("PA");
  graph[kGraphACCOverAllACRej] = GetGraphRej(graph[kGraphACCOverAllAC] , listRejectedRuns, "Acc/All AC" ,  meanDummy, 1);

  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fAccOverAll[0]->Draw();


  pad2->cd();


  graph[kGraphACCOverAllACS2]->Draw("PA");
  graph[kGraphACCOverAllACS2Rej] = GetGraphRej(graph[kGraphACCOverAllACS2] , listRejectedRuns, "Acc/All [ACS2]" ,  meanDummy, 1);

  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  //  GetHisto(graph[kGraphACCOverAllACS2],fAccOverAll[1]);
  fAccOverAll[1]->Draw();



  pad3->cd();


  graph[kGraphACCOverAllSMH]->Draw("PA");
  graph[kGraphACCOverAllSMHRej] = GetGraphRej(graph[kGraphACCOverAllSMH] , listRejectedRuns, "Acc/All SMH" ,  meanDummy, 1);

  //fAccOverAll[i]

  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  //  GetHisto(graph[kGraphACCOverAllSMH],fAccOverAll[2]);
  fAccOverAll[2]->Draw();
  fAccOverAll[2]->GetXaxis()->SetTitle("Accepted / All");

  /*
  for(Int_t ipoint = 0; ipoint < graph[kGraphACCOverAllAC]->GetN(); ipoint++){
    fhistTest->Fill(graph[kGraphACCOverAllAC]->GetY()[ipoint]);
  }

  pad1->cd(4);
  fhistTest->Sumw2();
  fhistTest->SetMarkerStyle(25);
  fhistTest->GetXaxis()->SetTitle(ylabels[kGraphACCOverAllAC]);
  fhistTest->GetYaxis()->SetTitle("Entries");
  fhistTest->Draw();
  */
  //pad1->cd(2);





  c4a->Update();
  gSystem->ProcessEvents();
  c4a->SaveAs(Form("picturesLHC11hAOD50/c4a_%s.png",c4a->GetName()));



  c5a = new TCanvas("GraphNev","GraphNev",1000,700);
  
  TPad *pad1 = new TPad("pad1",
			"The pad with the function",0.01,0.01,0.99,0.94,0);
  pad1->Draw();


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();




  pad1->Divide(1,3);
  pad1->cd(1);
  graph[kGraphNevAC]->Draw("PA");

  pad1->cd(2);
  graph[kGraphNevAC2]->Draw("PA");


  pad1->cd(3);
  graph[kGraphNevSMH]->Draw("PA");
 


  c5a->Update();
  gSystem->ProcessEvents();
  c5a->SaveAs(Form("picturesLHC11hAOD50/c5a_%s.png",c5a->GetName()));





  c7a = new TCanvas("GraphBGOverAll","GraphACCOverAll",1000,700);
  

  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();


  c7a->cd();
  pad1->Draw();


  pad1->cd();

  graph[kGraphBGOverAllAC]->Draw("PA");
  graph[kGraphBGOverAllAC] = GetGraphRej(graph[kGraphBGOverAllAC] , listRejectedRuns, "BG/All AC" ,  meanDummy, 1);

  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fBGOverAll[0]->Draw();
 






  pad2->cd();

  graph[kGraphBGOverAllAC2]->Draw("PA");
  graph[kGraphBGOverAllAC2] = GetGraphRej(graph[kGraphBGOverAllAC2] , listRejectedRuns, "BG/All AC2" ,  meanDummy, 1);


  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fBGOverAll[1]->Draw();



  pad3->cd();

 
  graph[kGraphBGOverAllSMH]->Draw("PA");
  graph[kGraphBGOverAllSMHRej] = GetGraphRej(graph[kGraphBGOverAllSMH] , listRejectedRuns, "BG/All SMH" ,  meanDummy, 1);


  pad3a->cd();
  pad3a->cd()->SetGridx(1);

  fBGOverAll[2]->GetXaxis()->SetTitle("BG / All");
  fBGOverAll[2]->Draw();




  c7a->Update();
  gSystem->ProcessEvents();
  c7a->SaveAs(Form("picturesLHC11hAOD50/c7a_%s.png",c7a->GetName()));



  c7b = new TCanvas("GraphV0BGOverAll","GraphV0BGOverAll",1000,700);
    TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();


  c7b->cd();



  pad1->cd();

  graph[kGraphV0BGOverAllAC]->Draw("PA");
  graph[kGraphV0BGOverAllAC] = GetGraphRej(graph[kGraphV0BGOverAllAC] , listRejectedRuns, "V0BG/All AC" ,  meanDummy, 1);

  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fV0BGOverAll[0]->Draw();


  pad2->cd();
  graph[kGraphV0BGOverAllAC2]->Draw("PA");
  graph[kGraphV0BGOverAllAC2] = GetGraphRej(graph[kGraphV0BGOverAllAC2] , listRejectedRuns, "V0BG/All AC2" ,  meanDummy, 1);

  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fV0BGOverAll[1]->Draw();

  pad3->cd();
  graph[kGraphV0BGOverAllSMH]->Draw("PA");
  graph[kGraphV0BGOverAllSMHRej] = GetGraphRej(graph[kGraphV0BGOverAllSMH] , listRejectedRuns, "V0BG/All SMH" ,  meanDummy, 1);

  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fV0BGOverAll[2]->GetXaxis()->SetTitle("V0BG / All");
  fV0BGOverAll[2]->Draw();



  c7b->Update();
  gSystem->ProcessEvents();
  c7b->SaveAs(Form("picturesLHC11hAOD50/c7b_%s.png",c7b->GetName()));










  c7c = new TCanvas("GraphV0A_BGOverAll","GraphV0A_BGOverAll",1000,700);
   TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();


  c7c->cd();

  pad1->cd(); 


  graph[kGraphV0ABGOverAllAC]->Draw("PA");
  graph[kGraphV0ABGOverAllAC] = GetGraphRej(graph[kGraphV0ABGOverAllAC] , listRejectedRuns, "V0A_BG/All AC" ,  meanDummy, 1);

  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fV0ABGOverAll[0]->Draw();



  pad2->cd(); 
  graph[kGraphV0ABGOverAllAC2]->Draw("PA");
  graph[kGraphV0ABGOverAllAC2] = GetGraphRej(graph[kGraphV0ABGOverAllAC2] , listRejectedRuns, "V0A_BG/All AC2" ,  meanDummy, 1);

  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fV0ABGOverAll[1]->Draw();

  pad3->cd(); 
  graph[kGraphV0ABGOverAllSMH]->Draw("PA");
  graph[kGraphV0ABGOverAllSMHRej] = GetGraphRej(graph[kGraphV0ABGOverAllSMH] , listRejectedRuns, "V0A_BG/All SMH" ,  meanDummy, 1);


  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fV0ABGOverAll[2]->Draw();
  fV0ABGOverAll[2]->GetXaxis()->SetTitle("V0ABG / All");


  c7c->Update();
  gSystem->ProcessEvents();
  c7c->SaveAs(Form("picturesLHC11hAOD50/c7c_%s.png",c7c->GetName()));





  c7d = new TCanvas("GraphV0C_BGOverAll","GraphV0C_BGOverAll",1000,700);
  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();


  c7d->cd();

  pad1->cd();
  graph[kGraphV0CBGOverAllAC]->Draw("PA");
  graph[kGraphV0CBGOverAllAC] = GetGraphRej(graph[kGraphV0CBGOverAllAC] , listRejectedRuns, "V0C_BG/All AC" ,  meanDummy, 1);

  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fV0CBGOverAll[0]->Draw();


  pad2->cd();
  graph[kGraphV0CBGOverAllAC2]->Draw("PA");
  graph[kGraphV0CBGOverAllAC2] = GetGraphRej(graph[kGraphV0CBGOverAllAC2] , listRejectedRuns, "V0C_BG/All AC2" ,  meanDummy, 1);

  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fV0CBGOverAll[1]->Draw();


	   pad3->cd();
  graph[kGraphV0CBGOverAllSMH]->Draw("PA");
  graph[kGraphV0CBGOverAllSMHRej] = GetGraphRej(graph[kGraphV0CBGOverAllSMH] , listRejectedRuns, "V0A_BG/All SMH" ,  meanDummy, 1);


  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fV0CBGOverAll[2]->Draw();
  fV0CBGOverAll[2]->GetXaxis()->SetTitle("V0CBG / All");



  c7d->Update();
  gSystem->ProcessEvents();
  c7d->SaveAs(Form("picturesLHC11hAOD50/c7d_%s.png",c7d->GetName()));





  c8a = new TCanvas("GraphAfterOverBefore","GraphAfterOverBefore",1000,700);
  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();
  
  pad1->cd();
  
  graph[kGraphNevACratioAfter]->Draw("PA");
  graph[kGraphNevACratioAfterRej] = GetGraphRej(graph[kGraphNevACratioAfter] , listRejectedRuns, "V0BG/All AC" ,  meanDummy, 1);
  
  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fAfterOverBefore[0]->Draw();


  pad2->cd();
  graph[kGraphNevAC2ratioAfter]->Draw("PA");
  graph[kGraphNevAC2ratioAfterRej] = GetGraphRej(graph[kGraphNevAC2ratioAfter] , listRejectedRuns, "V0BG/All AC2" ,  meanDummy, 1);

  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fAfterOverBefore[1]->Draw();

  pad3->cd();
  graph[kGraphNevSMHratioAfter]->Draw("PA");
  graph[kGraphNevSMHratioAfterRej] = GetGraphRej(graph[kGraphNevSMHratioAfter] , listRejectedRuns, "V0BG/All SMH" ,  meanDummy, 1);

  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fAfterOverBefore[2]->Draw();
  fAfterOverBefore[2]->GetXaxis()->SetTitle("After_PS / Before_PS");


  c8a->Update();
  gSystem->ProcessEvents();
  c8a->SaveAs(Form("picturesLHC11hAOD50/c8a_%s.png",c8a->GetName()));





  c8b = new TCanvas("GraphNevBefore","GraphNevBefore",1000,700);
  
  TPad *pad1 = new TPad("pad1",
			"The pad with the function",0.01,0.01,0.99,0.94,0);
  pad1->Draw();


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();




  pad1->Divide(1,3);

  pad1->cd(1);
  graph[kGraphNevACBefore]->Draw("PA");

  pad1->cd(2);
  graph[kGraphNevAC2Before]->Draw("PA");

  pad1->cd(3);
  graph[kGraphNevSMHBefore]->Draw("PA");



  c8b->Update();
  gSystem->ProcessEvents();
  c8b->SaveAs(Form("picturesLHC11hAOD50/c8b_%s.png",c8b->GetName()));





  c8c = new TCanvas("GraphNevAfter","GraphNevAfter",1000,700);
  
  TPad *pad1 = new TPad("pad1",
			"The pad with the function",0.01,0.01,0.99,0.94,0);
  pad1->Draw();


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();




  pad1->Divide(1,3);

  pad1->cd(1);
  graph[kGraphNevACAfter]->Draw("PA");


  pad1->cd(2);
  graph[kGraphNevAC2After]->Draw("PA");


  pad1->cd(3);

  graph[kGraphNevSMHAfter]->Draw("PA");



  c8c->Update();
  gSystem->ProcessEvents();
  c8c->SaveAs(Form("picturesLHC11hAOD50/c8c_%s.png",c8c->GetName()));








  c9a = new TCanvas("GraphV0AOverV0C","GraphV0AOverV0C",1000,700);
  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();

  pad1->cd();
  graph[kGraphV0AOverV0CAC]->Draw("PA");
  graph[kGraphV0AOverV0CACRej] = GetGraphRej(graph[kGraphV0AOverV0CAC] , listRejectedRuns, "Nev_V0A/Nev_V0C AC" ,  meanDummy, 1);

  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fV0AOverV0C[0]->Draw();


  pad2->cd();
  graph[kGraphV0AOverV0CAC2]->Draw("PA");
  graph[kGraphV0AOverV0CAC2Rej] = GetGraphRej(graph[kGraphV0AOverV0CAC2] , listRejectedRuns, "Nev_V0A/Nev_V0C AC2" ,  meanDummy, 1);

  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fV0AOverV0C[1]->Draw();


  pad3->cd();
  graph[kGraphV0AOverV0CSMH]->Draw("PA");
  graph[kGraphV0AOverV0CSMHRej] = GetGraphRej(graph[kGraphV0AOverV0CSMH] , listRejectedRuns, "Nev_V0A/Nev_V0C SMH" ,  meanDummy, 1);

  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fV0AOverV0C[2]->Draw();
  fV0AOverV0C[2]->GetXaxis()->SetTitle("V0A / VOC");



  c9a->Update();
  gSystem->ProcessEvents();
  c9a->SaveAs(Form("picturesLHC11hAOD50/c9a_%s.png",c9a->GetName()));



  c10a = new TCanvas("GraphFO1OverAll","GraphFO1OverAll",1000,700);
  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();  


  pad1->cd();
  graph[kGraphFO1OverAllAC]->Draw("PA");
  graph[kGraphFO1OverAllACRej] = GetGraphRej(graph[kGraphFO1OverAllAC] , listRejectedRuns, "Nev_FO1/All AC" ,  meanDummy, 1);

  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fF01OverAll[0]->Draw();


  pad2->cd();
  graph[kGraphFO1OverAllAC2]->Draw("PA");
  graph[kGraphFO1OverAllAC2Rej] = GetGraphRej(graph[kGraphFO1OverAllAC2] , listRejectedRuns, "Nev_FO1/All AC2" ,  meanDummy, 1);

  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fF01OverAll[1]->Draw();


  pad3->cd();
  graph[kGraphFO1OverAllSMH]->Draw("PA");
  graph[kGraphFO1OverAllSMHRej] = GetGraphRej(graph[kGraphFO1OverAllSMH] , listRejectedRuns, "Nev_FO1/All SMH" ,  meanDummy, 1);

  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fF01OverAll[2]->Draw();
  fF01OverAll[2]->GetXaxis()->SetTitle("F01 / All");




  c10a->Update();
  gSystem->ProcessEvents();
  c10a->SaveAs(Form("picturesLHC11hAOD50/c10a_%s.png",c10a->GetName()));




  c10b = new TCanvas("GraphFO2OverAll","GraphFO2OverAll",1000,700);
  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();  


  pad1->cd();
  graph[kGraphFO2OverAllAC]->Draw("PA");
  graph[kGraphFO2OverAllACRej] = GetGraphRej(graph[kGraphFO2OverAllAC] , listRejectedRuns, "Nev_FO2/All AC" ,  meanDummy, 1);

  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fF02OverAll[0]->Draw();
  fF02OverAll[0]->GetXaxis()->SetTitle("F02 / All");


  pad2->cd();
  graph[kGraphFO2OverAllAC2]->Draw("PA");
  graph[kGraphFO2OverAllAC2Rej] = GetGraphRej(graph[kGraphFO2OverAllAC2] , listRejectedRuns, "Nev_FO2/All AC2" ,  meanDummy, 1);

  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fF02OverAll[1]->Draw();
  fF02OverAll[1]->GetXaxis()->SetTitle("F02 / All");

  pad3->cd();
  graph[kGraphFO2OverAllSMH]->Draw("PA");
  graph[kGraphFO2OverAllSMHRej] = GetGraphRej(graph[kGraphFO2OverAllSMH] , listRejectedRuns, "Nev_FO2/All SMH" ,  meanDummy, 1);

  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fF02OverAll[2]->Draw();
  fF02OverAll[2]->GetXaxis()->SetTitle("F02 / All");





  c10b->Update();
  gSystem->ProcessEvents();
  c10b->SaveAs(Form("picturesLHC11hAOD50/c10b_%s.png",c10b->GetName()));



  cout << "Files Statistics" << endl;
  cout << " Total     [" << ifile         << "]" << endl;
  cout << " Available [" << ifileGood     << "]" << endl;  
  cout << " Not Empty [" << ifileNotEmpty << "]" <<  endl;

  cout << "" << endl;
  cout << "All Events" << endl <<  endl;
  tableEvts.PrintTable("CSV");

  cout << "Empty or missing files" << endl<< endl;  
  listEmptyRuns->Sort();
  listEmptyRuns->Print();

  cout << "Suspicious Runs" << endl << endl;
  listRejectedRuns->Sort();
  listRejectedRuns->Print();


 

  

}

TH1F * GetEmpty(const char * name, Int_t nfile) {
  TH1F * he01 = new TH1F(TString("hempty")+name, "hempty", nfile, -0.5, nfile-0.5);
  for(Int_t ilab = 0; ilab < nfile; ilab++){
    he01->GetXaxis()->SetBinLabel(ilab+1, Form("%d", runs[ilab]));
  }
  he01->SetMinimum(0);
  he01->SetMaximum(1);
  he01->SetXTitle("run");
  he01->SetYTitle(name);
  return he01;
}


const char * GetLocalFileName(Int_t run, const char * suffix, const char * path) {
  // returns the filename of the local copy of the event_stat file
  static TString name;
  //  name.Form("%s/event_stat_%s_%d.root", path, suffix, run);
  name.Form("%s/event_stat_%d.root", path, run);
  return name.Data();

}

double meanMed(double* vec, int np, double nsigmaCut, int &nrej, int *rejList)
{
  // compute the mean of the array "vec" rejecting the outliers
  // if rejlist array is provided, fill indices of rejected values
  // This method works assuming that the "good" points are nearly normally 
  // distrubted around the mean (i.e. there is no significant trend)
  //
  // create copy of the vector
  double *vec1 = new Double_t[np];
  memcpy(vec1, vec, np*sizeof(double));
  //
  // get median of the vector as a robust initial estimate of the mean
  cout << "Points " << np << endl;
  double md = GetMedian(vec1,np);
  //
  // compute squared differences to median, 
  for (int i=0;i<np;i++) {
    vec1[i] = TMath::Abs(vec1[i] - md);
  }
  //
  // compute median squared difference for robust estimate of the sigma
  double sigmd = GetMedian(vec1,np);
  //
  printf("Median Value: %+e | Median abs residual: %e\n",md,sigmd);
  //
  sigmd /= 0.6745; // go to sigma assuming normal distribution
  printf("Estimate of sigma: %+e\n",sigmd);
  // if the estimate of the sigma is null, do not look for outliers

  cout<<"md="<<md<<endl;

  if(!sigmd) return md;

 


  // compute mean rejecting more than nsigmaCut deviations
  double mean = 0;
  int npok = 0;
  for (int i=0;i<np;i++) {
    double dev =  TMath::Abs(vec[i]-md)/sigmd;
    if ( dev < nsigmaCut ) {
      mean += vec[i];
      npok++;
    }
    else {
      printf("Reject#%d (Y=%+e) : deviation: %e\n",i,vec[i],dev);
      if (rejList) rejList[nrej] = i; 
      nrej++;
    }
  }
  //
  delete vec1;
  return npok ? mean/npok : 0;
  //


  cout<<"ending meanMed"<<endl;

}

double GetMedian(double* arr, int n)
{
  // get median by Wirths algroithm
  int i=0,j=0,l=0,m=0;
  double x;
  l=0 ; m=n-1;
  int k = (n&1) ? n/2 : n/2-1;
  while (l<m) {
    x=arr[k] ;
    i=l ;
    j=m ;
    do {
      //      cout << arr[i] << " " << arr[j] << i<<","<<j << " " << arr[k]  << " " << k << " " << x<< endl;
      
      while (arr[i]<x) i++ ;
      while (x<arr[j]) j-- ;
      if (i<=j) { // swap elements	  
	//	cout << "Swapping" << endl;	  
	double t = arr[i];
	arr[i] = arr[j];
	arr[j] = t;
	i++ ; j-- ;
      }      
    } while (i<=j) ;
    if (j<k) l=i ;
    if (k<i) m=j ;
  }
  return arr[k] ;
}

// double GetMedian(double* arr, int n)
// {
//   // get median by Wirths algroithm
//   int i,j,l,m;
//   double x;
//   l=0 ; m=n-1;
//   int k = (n&1) ? n/2 : n/2-1;
//   while (l<m) {
//     x=arr[k] ;
//     i=l ;
//     j=m ;
//     do {
//       //      cout << i << " " << j << endl;
      
//       //skip runs which were not set (value is -1)
//       while (arr[i]<x) i++ ;
//       while (x<arr[j]) j-- ;
//       if(i<=j) { // swap elements
// 	// if((TMath::Abs(arr[i]+1)< 0.0001) && ((TMath::Abs(arr[j]+1)<0.0001))){
// 	//   i++ ; j-- ;
// 	// } else {
// 	  double t = arr[i];
// 	  arr[i] = arr[j];
// 	  arr[j] = t;
// 	  i++ ; j-- ;
// 	  //	}
//       }      
//     } while (i<=j) ;
//     if (j<k) l=i ;
//     if (k<i) m=j ;
    
//   }
//   return arr[k] ;
// }

TGraphErrors * GetGraphRej(TGraphErrors * gr, TList * rejRunList, const char * reason, Float_t &mean, Bool_t doDraw) {

  //Returns a new graph containing only rejected points

  const Double_t nsigmaCut = 3;

  int *rejList = new Int_t[gr->GetN()];
  int nrej = 0;

  Double_t * array = new Double_t[gr->GetN()];
  Int_t * correspondenceFullArray = new Int_t[gr->GetN()];
  Int_t npoint = 0; 
  //exclude from the array all the -1
  for(Int_t ipoint = 0; ipoint < gr->GetN(); ipoint++){
    if (TMath::Abs(gr->GetY()[ipoint]+1)>0.001) { // skip points for which there is no file (==-1)
      array[npoint] = gr->GetY()[ipoint];
      correspondenceFullArray[npoint] = ipoint;
      npoint++;
    } 
  }

  cout<<"calling meanMed"<<endl;
  mean = meanMed(array,npoint,nsigmaCut, nrej, rejList);

  cout<<"nrej="<<nrej<<"  mean="<<mean<<endl;
  
 
    
  TGraphErrors *grrej = new TGraphErrors(nrej);
  for (int i=0;i<nrej;i++) {
    grrej->SetPoint(i, gr->GetX()[correspondenceFullArray[rejList[i]]], gr->GetY()[correspondenceFullArray[rejList[i]]]);
    grrej->SetPointError(i, gr->GetEX()[correspondenceFullArray[rejList[i]]], gr->GetEY()[correspondenceFullArray[rejList[i]]]);
    if(!knownProblems.Contains(Form("%d",(int)gr->GetX()[correspondenceFullArray[rejList[i]]])))
      rejRunList->Add(new TObjString(Form("[%d], %s",  (int)gr->GetX()[correspondenceFullArray[rejList[i]]], reason)));
  }
  grrej->SetMarkerColor(kRed);
  grrej->SetMarkerStyle(gr->GetMarkerStyle());
  


  delete rejList;

  if(doDraw) {
    Float_t minXDraw = gr->GetXaxis()->GetBinLowEdge(1);
    Float_t maxXDraw = gr->GetXaxis()->GetBinLowEdge(gr->GetXaxis()->GetNbins());

    grrej->Draw("P");
    TLine * l = new TLine (minXDraw, mean,maxXDraw, mean);
    l->SetLineStyle(kDashed);
    l->Draw();


  }

  return grrej;



} 

void loadlibs()
{
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libMinuit");
  gSystem->Load("libPWG2spectra");
  gSystem->Load("libPWG0base"); 
}
