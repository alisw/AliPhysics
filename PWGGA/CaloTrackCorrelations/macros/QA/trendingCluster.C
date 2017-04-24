/// \file trendingCluster.C
/// \ingroup CaloTrackCorrMacrosQA
/// \brief Do EMCal QA  trending plots
///
/// Example macro to check QA outputs from the histograms itself
/// Executed with Root 
///
/// The input file needs: 
/// 1.  the root output from AliAnaCalorimeterQA (or QA train) run by run ie runnumber.root 
/// placed in the directory period/pass/
/// 2.  the run list (runlist.txt) of these run output placed in the directory period/pass/ runlist mean : index runnumber//
/// Trigger 
/// 2 trigger options " MB"  : CaloQA_default output from train
///                   " EMC" : CaloQA_EMC7
///
/// Configured by options of EMCAL  checker with:
/// a. cell multiplicity
/// b. cluster multiplicity/event
/// c. Cells per Cluster 
/// d. mean cluster energy 
///
/// more checker values could be added based on request (more canvas)
///
/// \author Yaxian Mao, Wuhan
/// \author Marie Germain,  SUBATECH


/* $Id:  $ */
//--------------------------------------------------

//
// Author: Yaxian Mao
// Modified M. Germain
//

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3D.h>
#include <Riostream.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGrid.h>
#include <TFileMerger.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGridCollection.h>
#include <TROOT.h>
#include <TGridResult.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <stdio.h>

#include <string> 
#include <fstream>
#include <iostream>
#include <sstream>  // Required for stringstreams

using namespace std;

void  trendingCluster(TString fCalorimeter = "EMCAL", TString period = "LHC11h", TString pass = "pass1_HLT", const Int_t n = 10, TString fTrigger = "MB"){
  
  FILE * pFile;
  
  TString file = "";
  if (fTrigger=="EMC") file = "/scratch/alicehp2/germain/QA/"+period+"/"+ pass + "/runlist.txt" ;
  else                 file = "/scratch/alicehp2/germain/QA/"+period+"/"+ pass + "/runlistMB.txt" ;
  
  //cout<<file<<endl;
  pFile = fopen(file.Data(), "r"); //open the text file where include the run list and correct run index
  
  cout<<file<<endl;
  cout << " fcalo: " << fCalorimeter << "; period: " << period << "; pass: " << pass << "  trigger "<<fTrigger<<  endl; 
  char outfilename [100]     ;
  
  Int_t index[500];
  Int_t p;
  Int_t q;
  Int_t ncols;
  Int_t nlines = 0 ;
  Int_t RunId[500] ;
  Double_t x[500] ;
  Double_t xrun[500] ;
  //TString label; 
  
  while (1){
    ncols = fscanf(pFile,"%d  %d ",&p,&q);
    if (ncols< 0) break;
    x[nlines]=p;
    index[nlines]=p;
    RunId[nlines]=q;
    xrun[nlines]=1.*q;
    nlines++;
  }
  fclose(pFile);
  const Int_t nRun = nlines ;
  
  Double_t xe[nRun] ;
  Double_t nEvent[nRun] ;
  //  Double_t TimeMean[nRun] ;
  //  Double_t TimeRMS[nRun] ;
  Double_t CellMean[nRun] ;
  Double_t CellRMS[nRun] ;
  Double_t ClusterMean[nRun] ;
  Double_t ClusterRMS[nRun] ;
  Double_t EtotalMean[nRun] ;//total energy deposited per event
  Double_t EtotalRMS[nRun] ;
  
  Double_t CellPerClusterMean[nRun] ;
  Double_t CellPerClusterRMS[nRun] ;
  Double_t ECell1Mean[nRun] ;//total energy deposited per event without 1 cell clusters
  Double_t ECell1RMS[nRun] ;
  
  TFile * f ;
  TDirectoryFile * dir;
  TList * outputList;
  
  TH1F * fhNEvents;
  TH1F * fhE;
  TH1F * fhNClusters;
  TH1F * fhNCells;
  TH2F * fhIM ;
  TH3D * fhNCellsPerCluster ;
  TH2F * fhTimeAmp ;
  
  TH1F * NCells[n];
  TH1F * NClusters[n];
  TH2F * NCellsPerCluster[n];
  TH1F * E[n];
  
  TGraphErrors * AverNcellsSM[n];
  TGraphErrors * AverNclustersSM[n];
  TGraphErrors * AverNcellsPerClusterSM[n];
  TGraphErrors * AverESM[n];
  TGraphErrors * AverTimeSM[n];
  TGraphErrors * AverMggSM[n];
  TGraphErrors * AverEcell1SM[n];
  
  Double_t CellMeanSM[n][nRun] ;
  Double_t CellRMSSM[n][nRun] ;
  Double_t ClusterMeanSM[n][nRun] ;
  Double_t ClusterRMSSM[n][nRun] ;
  Double_t EtotalMeanSM[n][nRun] ;//total energy deposited per event
  Double_t EtotalRMSSM[n][nRun] ;
  Double_t CellPerClusterMeanSM[n][nRun] ;
  Double_t CellPerClusterRMSSM[n][nRun] ;
  Double_t ECell1MeanSM[n][nRun] ;//total energy deposited per event without 1 cell clusters
  Double_t ECell1RMSSM[n][nRun] ;
  
  
  TString namefile = "/scratch/alicehp2/germain/QA/"+period+"/"+pass+"/"+ fCalorimeter + period + pass + fTrigger+"data.txt";
  
  fstream QAData(namefile, ios::out); //write the QA check values at the end
  
  cout << " namefile " << namefile << endl;
  cout << " nRun " << nRun <<  " index(nRun) " << index[nRun-1]  << endl;
  
  for(Int_t i = 0 ; i < nRun ; i++){
    
    xe[i] = 0.5 ; 
    
    TString name = "/scratch/alicehp2/germain/QA/"+period +"/"+ pass + "/";
    
    name += RunId[i] ;
    name += ".root";
    f = TFile::Open(name.Data(),"read") ;
    if (!f) continue; 
    
    //  cout << " i = " << i << " file opend Output" <<  RunId[i]<< endl;
    
    if(fTrigger=="EMC"){        dir = (TDirectoryFile *)f->Get("CaloQA_EMC7");
      outputList = (TList*)dir->Get("CaloQA_EMC7");
    }
    else{
      dir = (TDirectoryFile *)f->Get("CaloQA_default");
      outputList = (TList*)dir->Get("CaloQA_default");
    }
    
    //define the averages for checking Histograms
    //
    
    fhNEvents =(TH1F *)outputList->FindObject("hNEvents");
    nEvent[i]=fhNEvents->GetEntries();
    cout <<  " Run " << RunId[i]<< " nevent " << nEvent[i]<< endl;
    if( nEvent[i] == 0) continue ;
    if( nEvent[i] <  2) continue ;
    fhE = (TH1F *)outputList->FindObject(fCalorimeter+"_hE");
    
    Double_t energy = 0. ;
    
    for(Int_t ibin = fhE->FindBin(0.3) ; ibin <fhE->FindBin(50.) ; ibin++){ //Starting from 0.3eV
      energy+=fhE->GetBinCenter(ibin)*fhE->GetBinContent(ibin);
    }
    EtotalMean[i]=energy/fhE->Integral(fhE->FindBin(0.3), fhE->FindBin(50.)) ;
    EtotalRMS[i]=fhE->GetMeanError();
    
    
    
    //for single module check
    for(Int_t ism = 0 ; ism < n ; ism++){
      TString nameNCell = Form("%s_hNCells_Mod%d",fCalorimeter.Data(),ism);
      
      TString nameNCluster = Form("%s_hNClusters_Mod%d",fCalorimeter.Data(),ism);
      TString nameNCellPerCluster = Form("%s_hNCellsPerCluster_Mod%d",fCalorimeter.Data(),ism);
      TString nameE = Form("%s_hE_Mod%d",fCalorimeter.Data(),ism);
      
      NCells[ism] = (TH1F*)outputList->FindObject(nameNCell.Data());
      NClusters[ism] = (TH1F*)outputList->FindObject(nameNCluster.Data());
      NCellsPerCluster[ism] = (TH2F*)outputList->FindObject(nameNCellPerCluster.Data());
      E[ism] = (TH1F*)outputList->FindObject(nameE.Data());
      CellMeanSM[ism][i]=NCells[ism]->GetMean();
      CellRMSSM[ism][i]=NCells[ism]->GetMeanError();
      ClusterMeanSM[ism][i]=NClusters[ism]->GetMean();
      ClusterRMSSM[ism][i]=NClusters[ism]->GetMeanError();
      CellPerClusterMeanSM[ism][i]=NCellsPerCluster[ism]->GetMean(2);
      CellPerClusterRMSSM[ism][i]=NCellsPerCluster[ism]->GetMeanError(2);
      //   cout<<"SM = "<<ism<<" Mean : = "<<CellPerClusterMeanSM[ism][i]<<endl ;
      ECell1MeanSM[ism][i] =NCellsPerCluster[ism]->ProjectionX("",2,300,"")->Integral(5,100)/(nEvent[i]);
      ECell1RMSSM[ism][i] =NCellsPerCluster[ism]->ProjectionX("",2,300,"")->GetMeanError();
      Double_t energySM = 0. ;
      for(Int_t ibin = E[ism]->FindBin(0.3) ; ibin <E[ism]->FindBin(50.) ; ibin++){ //starting from 0.3GeV
        energySM+=E[ism]->GetBinCenter(ibin)*(E[ism]->GetBinContent(ibin));
      }
      EtotalMeanSM[ism][i]=energySM/(E[ism]->Integral(E[ism]->FindBin(0.3),E[ism]->FindBin(50.)));
      
      EtotalRMSSM[ism][i]=E[ism]->GetMeanError();
      
      if(ism==0) {
        fhNCells = (TH1F*)NCells[ism]->Clone("NCells");
        fhNClusters = (TH1F*)NClusters[ism]->Clone("NClusters");  
      }
      else {
        fhNCells->Add(NCells[ism],1);
        fhNClusters->Add(NClusters[ism],1);
      }
    } //per SM loop
    ClusterMean[i]=fhNClusters->GetMean();
    ClusterRMS[i]=fhNClusters->GetMeanError();  
    CellMean[i]=fhNCells->GetMean();
    CellRMS[i]=fhNCells->GetMeanError();  
    outputList->Clear() ; 
    dir->Close();
    f->Close();
    f=NULL;
    dir=NULL;
    outputList=NULL;
    
    //if you want to write all the QA check values at the end, otherwise just comment out below
    // becareful with different detectors as the check output are different since different modules/SM for different detetcor
    
    QAData <<i+1<<"   "<< RunId[i] <<"    "<< nEvent[i]        
    //   <<"   "<< RunId[i]<<"   "<<CellMean[i]
    //   <<" "<< CellMeanSM[0][i] <<"   "<< CellMeanSM[1][i]
    //     <<"   "<<CellMeanSM[2][i] 
    //     <<"   "<< RunId[i] <<" "<<ClusterMean[i] 
    //     <<" "<< ClusterMeanSM[0][i]<<" "<< ClusterMeanSM[1][i]
    //     <<" "<< ClusterMeanSM[2][i]
    //     <<"   "<< RunId[i]<<"   "<< CellPerClusterMean[i]
    //     <<" "<< CellPerClusterMeanSM[0][i]<<" "<< CellPerClusterMeanSM[1][i]
    //     <<" "<< CellPerClusterMeanSM[2][i]
    //     <<"   "<< RunId[i]<<"   "<< EtotalMean[i]
    //     <<" "<< EtotalMeanSM[0][i] <<" "<< EtotalMeanSM[1][i]
    //     <<" "<< EtotalMeanSM[2][i] 
    
    <<"\n" ; 
    
  } // end loop on nrun
  
  
  QAData.close();
  
  TString base = "/scratch/alicehp2/germain/QA/";
  base += period ;
  base += "/";
  base += pass ;
  base += "/";
  
  base += fTrigger;
  TString ClusterAverages ; ClusterAverages = base +  "ClusterAverages.gif";
  TString Entries  ;  Entries = base + "Nentries.gif";
  TString ClusterAveragesEnergy ; ClusterAveragesEnergy = base +  "ClusterAveragesEnergy.gif";
  TString ClusterAveragesEntries ; ClusterAveragesEntries = base +  "ClusterAveragesEntries.gif";
  TString ClusterAveragesCells ; ClusterAveragesCells = base +  "ClusterAveragescells.gif";
  
  
  cout << "c11 nEvents" << endl;
  //just for the canvas defination
  
  cout << " index(0)" << index[nRun-1] << " index(20) " << index[20] <<endl;
  TH1F * dummy = new TH1F("dummy", "dummy", nRun, 0., nRun+0.5); 
  //  TH1F * dummy = new TH1F("dummy", "dummy", index[nRun-1], 0., index[nRun-1]+0.5); 
  dummy->SetTitle("") ; 
  dummy->SetStats(kFALSE) ; 
  dummy->SetAxisRange(0., nRun, "X") ; 
  
  for(Int_t i = 0 ; i < nRun ; i++){
    TString label = " ";
    label+=RunId[i];
    cout <<" run "<< RunId[i] << " label " <<label << endl;
    
    dummy->GetXaxis()->SetBinLabel(i+1,label.Data());
    dummy->GetXaxis()->LabelsOption("v");
  }
  
  
  
  //number of events passes  physics selection for each run
  TCanvas  * c11 = new TCanvas("nEvents", "nEvents", 1000, 500);
  c11->SetFillColor(0);
  c11->SetBorderSize(0);
  c11->SetFrameBorderMode(0); 
  gStyle->SetOptStat(0);  
  c11->SetLogy();
  c11->SetGrid();
  // dummy->GetXaxis()->SetTitle("RUN");    
  dummy->GetXaxis()->SetTitleOffset(0.05);
  dummy->GetYaxis()->SetTitle("N_{events}"); 
  dummy->SetMinimum(1.) ;  //should addjust based on the statistics 
  dummy->SetMaximum(1.e6) ; //should addjust based on the statistics 
  dummy->Draw();
  TGraph * nEvents = new TGraph(nRun, x, nEvent);
  nEvents->SetMarkerStyle(20);
  nEvents->SetMarkerColor(1);
  nEvents->SetLineColor(2);
  nEvents->Draw("same PL") ;
  c11->Update();
  
  
  if (fTrigger=="MB")sprintf(outfilename,"nEventMB.gif");
  if (fTrigger=="EMC")sprintf(outfilename,"nEventEMC.gif");
  
  
  c11->SaveAs(Entries);
  
  
  cout << "c1 Aver NCell" << endl;
  
  TCanvas  * c1 = new TCanvas("AverNCell", "AverNCell", 600, 600);
  // c1->SetLogy();
  c1->SetFillColor(0);
  c1->SetBorderSize(0);
  c1->SetFrameBorderMode(0); 
  gStyle->SetOptStat(0);  
  TH1F * h1 = (TH1F*)dummy->Clone("");
  h1->GetXaxis()->SetTitle("RUN Index");    
  h1->GetYaxis()->SetTitle("<N_{cells}>");  
  h1->SetMinimum(0.) ; 
  h1->SetMaximum(10.) ; 
  if(fCalorimeter=="EMCAL") h1->SetMaximum(5.) ;
  h1->Draw();
  
  TGraphErrors * AverNcells = new TGraphErrors(nRun, x, CellMean, xe, CellRMS);
  
  AverNcells->SetMarkerColor(1);
  AverNcells->SetMarkerStyle(20);
  AverNcells->Draw("same P") ;  
  for(Int_t ism = 0 ; ism < n ; ism++){
    AverNcellsSM[ism] = new TGraphErrors(nRun, x, CellMeanSM[ism], xe, CellRMSSM[ism]);
    AverNcellsSM[ism]->SetMarkerColor(ism+2);
    AverNcellsSM[ism]->SetMarkerStyle(21+ism);
    AverNcellsSM[ism]->Draw("same P");
  }
  
  
  
  TLegend  * l1 = new TLegend(0.4, 0.6, 0.75, 0.85);
  l1->SetFillColor(0);
  l1->SetBorderSize(0);
  l1->SetTextSize(0.02);
  l1->AddEntry(AverNcells, "<# of cells>", "") ;
  l1->AddEntry(AverNcells, Form("det = %s",fCalorimeter.Data()), "") ;
  l1->AddEntry(AverNcells,"average", "p");
  for(Int_t ism = 0 ; ism < n ; ism++){
    TString projname = Form("SM_ %d",ism);
    l1->AddEntry(AverNcellsSM[ism],projname.Data(), "p");
  }
  l1->Draw("same");
  c1->Update();
  
  
  TCanvas  * c200 = new TCanvas("ClusterAverages", "ClusterAverages", 1000, 500);
  c200->SetFillColor(0);
  c200->SetBorderSize(0);
  c200->SetFrameBorderMode(0); 
  c200->SetGrid();
  
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02);
  gPad->SetGrid();
  
  TH1F * h2 = (TH1F*)dummy->Clone("");
  
  h2->GetYaxis()->SetTitle("<N_{clusters}>/event");  
  
  if (fTrigger=="EMC") h2->SetMaximum(70) ;
  else                 h2->SetMaximum(50) ;
  h2->SetMinimum(10) ; // for Pb Pb
  h2->Draw();
  
  TGraphErrors * AverNclusters = new TGraphErrors(nRun, x, ClusterMean, xe, ClusterRMS);
  AverNclusters->SetMarkerStyle(20);
  AverNclusters->SetMarkerColor(1);
  AverNclusters->Draw("same P") ;  
  
  for(Int_t ism = 0 ; ism < n ; ism++){
    AverNclustersSM[ism] = new TGraphErrors(nRun, x, ClusterMeanSM[ism], xe, ClusterRMSSM[ism]);
    AverNclustersSM[ism]->SetMarkerColor(ism+2);
    AverNclustersSM[ism]->SetMarkerStyle(21+ism);
    
    AverNclustersSM[ism]->Draw("same P");  
  }
  
  
  TLegend  * l200 = new TLegend(0.4, 0.6, 0.75, 0.85);
  l200->SetFillColor(0);
  l200->SetBorderSize(0);
  l200->SetTextSize(0.02);
  l200->AddEntry(AverNclusters, "<# of clusters>", "") ;
  l200->AddEntry(AverNclusters, Form("det = %s",fCalorimeter.Data()), "") ;
  l200->AddEntry(AverNclusters,"average", "p");
  
  for(Int_t ism = 0 ; ism < n ; ism++){
    TString projname = Form("SM_ %d",ism);
    l200->AddEntry(AverNclustersSM[ism],projname.Data(), "p");
  }
  
  l200->Draw("same");
  c200->Update();
  c200->SaveAs(ClusterAveragesEntries);
  
  TCanvas  * c201 = new TCanvas("ClusterAveragesEnergy", "ClusterAveragesEnergy", 1000, 500);
  c201->SetFillColor(0);
  c201->SetBorderSize(0);
  c201->SetFrameBorderMode(0); 
  c201->SetGrid();
  
  
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02);
  gPad->SetGrid();
  
  
  TH1F * h3 = (TH1F*)dummy->Clone("");
  
  h3->GetYaxis()->SetTitle("<E> (GeV)");  
  h3->SetMinimum(0.2) ; 
  h3->SetMaximum(1.) ; 
  h3->Draw();
  
  TGraphErrors * AverE = new TGraphErrors(nRun, x, EtotalMean, xe, EtotalRMS);
  AverE->SetMarkerStyle(20);
  AverE->SetMarkerColor(1);
  AverE->Draw("same P");
  
  for(Int_t ism = 0 ; ism < n ; ism++){
    AverESM[ism] = new TGraphErrors(nRun, x, EtotalMeanSM[ism], xe, EtotalRMSSM[ism]);
    AverESM[ism]->SetMarkerColor(ism+2);
    AverESM[ism]->SetMarkerStyle(21+ism);
    AverESM[ism]->Draw("same P");
    
  }
  
  TLegend  * l3 = new TLegend(0.4, 0.6, 0.75, 0.85);
  l3->SetFillColor(0);
  l3->SetBorderSize(0);
  l3->SetTextSize(0.02);
  l3->AddEntry(AverE, "<E>", "") ;
  l3->AddEntry(AverE, Form("det = %s",fCalorimeter.Data()), "") ;
  l3->AddEntry(AverE,"average", "p");
  
  for(Int_t ism = 0 ; ism < n ; ism++){
    TString projname = Form("SM_ %d",ism);
    l3->AddEntry(AverESM[ism],projname.Data(), "p");
  }
  
  l3->Draw("same");
  
  c201->SaveAs(ClusterAveragesEnergy);
  
  
  TCanvas  * c202 = new TCanvas("ClusterAveragesCells", "ClusterAveragesCells", 1000, 500);
  c202->SetFillColor(0);
  c202->SetBorderSize(0);
  c202->SetFrameBorderMode(0); 
  c202->SetGrid();
  
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02); 
  
  gPad->SetGrid();
  
  
  TH1F * h4 = (TH1F*)dummy->Clone("");
  
  h4->GetYaxis()->SetTitle("<N_{CellsPerCluster}>");  
  h4->SetMinimum(1.5) ; 
  h4->SetMaximum(5.5) ; 
  h4->Draw();
  
  
  TGraphErrors * AverCellPerCluster = new TGraphErrors(nRun, x, CellPerClusterMean, xe, CellPerClusterRMS);
  AverCellPerCluster->SetMarkerStyle(20);
  AverCellPerCluster->SetMarkerColor(1);
  
  for(Int_t ism = 0 ; ism < n ; ism++){
    AverNcellsPerClusterSM[ism] = new TGraphErrors(nRun, x, CellPerClusterMeanSM[ism], xe, CellPerClusterRMSSM[ism]);
    AverNcellsPerClusterSM[ism]->SetMarkerColor(ism+2);
    AverNcellsPerClusterSM[ism]->SetMarkerStyle(21+ism);
    AverNcellsPerClusterSM[ism]->Draw("same P");
    
  }
  
  TLegend  * l4 = new TLegend(0.4, 0.6, 0.75, 0.85);
  l4->SetFillColor(0);
  l4->SetBorderSize(0);
  l4->SetTextSize(0.02);
  l4->AddEntry(AverCellPerCluster, "<# of cells per cluster>", "") ;
  l4->AddEntry(AverCellPerCluster, Form("det = %s",fCalorimeter.Data()), "") ;
  l4->AddEntry(AverCellPerCluster,"average", "p");
  
  for(Int_t ism = 0 ; ism < n ; ism++){
    TString projname = Form("SM_ %d",ism);
    l4->AddEntry(AverNcellsPerClusterSM[ism],projname.Data(), "p");
  }
  
  l4->Draw("same");
  
  c202->SaveAs(ClusterAveragesCells);
  
}
