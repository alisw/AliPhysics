// $Id$
/*
 * Drawing macro for reading the THnSparse output of the 
 * HLT/QA/tasks/AliAnalysisTaskHLTCentralBarrel.cxx task.
 * 
 * The cuts are user defined in lines 156-158 as arguments of the 
 * function cutStsudies(...).
 * 
 * The input file contains information about the run number
 * and date that will appear in the canvas title.
 * 
 * The cuts are turned into strings and incorporated in the
 * name of the output file, which contains canvases with event
 * and track properties for HLT and offline.
 *  
 * Since the run information is available, there will be a new
 * folder created and the canvas ROOT files will be saved in there.
 *
 * The macro should be compiled before running:
 *
 * root[0] .L drawTHnSparse.C++
 * root[1] drawTHnSparse("../HLT-OFFLINE-CentralBarrel-comparison.root")
 *
 * @ingroup alihlt_qa
 * @author Kalliopi.Kanaki@ift.uib.no 
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TList.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"

#include <iostream>
//#include <sstream>
#include <cstdlib>

//using std::stringstream;
using std::endl;
#endif


//---------- forward declerations ---------------//

TString cutStudies( TCanvas* can1, TCanvas* can2, TCanvas* can3, TString folder,
//void cutStudies( TCanvas* can1, TCanvas* can2, TCanvas* can3, TString folder,
                    THnSparse* htrackHLT, THnSparse* htrackOFF, TText* hText,
            	    double minEta,   double maxEta,
	    	    int minTrackMult,int maxTrackMult,
            	    double minPt,    double maxPt,
	    	    double minDCAr,  double maxDCAr,
	    	    double minDCAz,  double maxDCAz,
	    	    int minTPCclus,  int maxTPCclus,
	    	    int minITSclus,  int maxITSclus,
		    int vertexStatus
                  );
void printStats(TH1D *hlt, TH1D *off);
void defineYaxisMax(TH1D *hlt, TH1D *off);
void printLegend(TLegend *l, TH1D *hlt, TH1D *off);
void plotAid(TCanvas* can, THnSparse* hHLT, THnSparse* hOFF, TText* hText, TH1D *hlt, TH1D *off, TLegend *l, int size);	   
void plot2D(TCanvas* can, THnSparse* h,
            double minEta,    double maxEta,
            double minPt,     double maxPt,
	    double minDCAr,   double maxDCAr,
	    double minDCAz,   double maxDCAz,
	    int minTPCclus,   int maxTPCclus,
	    int minITSclus,   int maxITSclus, 
 	    int minTrackMult, int maxTrackMult
           );
void plotEventQuantities(TCanvas* can, THnSparse* heventHLT, THnSparse* heventOFF, TText* hText);
TString fix1DTitle(const char* c);
TString fix2DTitle(const char* c1, const char* c2);
TString cutsToString( double minEta,	double maxEta,
 	              double minPt,	double maxPt,
	   	      double minDCAr,	double maxDCAr,
	   	      double minDCAz,	double maxDCAz,
	   	      int minTPCclus,	int maxTPCclus,
	   	      int minITSclus,	int maxITSclus, 
	              int minTrackMult, int maxTrackMult,
		      int VS
	            );
TString itoa(int i);
void plotTrackQuantities( TCanvas* can, THnSparse* htrackHLT, THnSparse* htrackOFF, TText* hText,
           	    	  double minEta,    double maxEta,
 	                  double minPt,     double maxPt,
	   	    	  double minDCAr,   double maxDCAr,
	   	    	  double minDCAz,   double maxDCAz,
	   	    	  int minTPCclus,   int maxTPCclus,
	   	    	  int minITSclus,   int maxITSclus, 
	                  int minTrackMult, int maxTrackMult,
			  int VS
                        );
			
vector<TString> outputNames;

//------------------------------------------------------------------//		
			
void drawTHnSparse(TString inputFile){
 
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(10);
  TH1::AddDirectory(kFALSE);

  TFile *file = TFile::Open(inputFile);
  if(!file){
    printf("Error: No file %s in folder.\n", inputFile.Data());
    return;
  }

  TList *list = static_cast<TList*>(file->Get("esd_thnsparse"));
  if(!list){
    printf("Error: No List contained in file %s.\n", inputFile.Data());
    return;
  }

  THnSparseF *heventHLT = (THnSparseF*)list->FindObject("fEventHLT"); 
  if(!heventHLT){
      printf("Error: There is no HLT THnSparse object in file %s\n", inputFile.Data());
  }
  THnSparseF *heventOFF = (THnSparseF*)list->FindObject("fEventOFF");  
  if(!heventOFF){
      printf("Error: There is no OFF THnSparse object in file %s\n", inputFile.Data());
  } 
  THnSparseF *htrackHLT = (THnSparseF*)list->FindObject("fTrackHLT");
  if(!htrackHLT){
      printf("Error: No HLT THnSparse object found\n");
      return;
  } 
  THnSparseF *htrackOFF = (THnSparseF*)list->FindObject("fTrackOFF");  
  if(!htrackOFF){
      printf("Error: No OFF THnSparse object found\n");
      return;
  }
      
  TText *hText = (TText*)list->FindObject("text");
  if(!hText) printf("No hText\n");
  
  TString t = "event properties for ";
  t+=hText->GetTitle();
  
  TString folder = hText->GetTitle();
  folder.ReplaceAll(" ",""); 
  folder.Remove(0,3);
  folder.Remove(6,15);
   
  gSystem->Exec("mkdir "+folder);
  
  TCanvas *can0 = new TCanvas("can0",t,                            900,600); can0->Divide(3,2);  
  TCanvas *can1 = new TCanvas("can1","track properties",          1100,900); can1->Divide(4,3);
  TCanvas *can2 = new TCanvas("can2","2-D HLT track correlations",1200,800); can2->Divide(4,2);
  TCanvas *can3 = new TCanvas("can3","2-D OFF track correlations",1200,800); can3->Divide(4,2);

  plotEventQuantities(can0,heventHLT,heventOFF,hText);
  can0->SaveAs(folder+"/event_properties.root");
  can0->SaveAs(folder+"/event_properties.png");
    
  TString s = "";                                                       // eta   mult      pt     DCAr    DCAz    TPCclus ITSclus  vertexStatus
  s = cutStudies(can1, can2, can3, folder, htrackHLT, htrackOFF, hText, -2, 2, 0, 20000, 0, 200, -80, 80, -80, 80,  0, 200, 0, 6, 2); outputNames.push_back(s); 
  s = cutStudies(can1, can2, can3, folder, htrackHLT, htrackOFF, hText, -2, 2, 0, 20000, 0, 200, -80, 80, -80, 80,  0, 200, 2, 6, 2); outputNames.push_back(s); 
  s = cutStudies(can1, can2, can3, folder, htrackHLT, htrackOFF, hText, -2, 2, 0, 20000, 0, 200, -80, 80, -80, 80,  0, 200, 4, 6, 2); outputNames.push_back(s); 
 
  TCanvas *ov = new TCanvas("ov","",1100,900);
  ov->Divide(4,3);

  TCanvas *ca[outputNames.size()]; 
  TFile   *ff[outputNames.size()]; 
  TPad    *pad[12]; 
  TH1D    *g[outputNames.size()];
  
  for(int j=1; j<12; j++){ // not 13, last pad is empty (TODO)        
    for(UInt_t i=0; i<outputNames.size(); i++){  
           
        ff[i] = TFile::Open(outputNames[i].Data());   
        if(!ff[i] || ff[i]->IsZombie()){
           printf("Non-existent, corrupted or zombie file %s\n", outputNames[i].Data());
           return;
        } 
        ca[i]  = (TCanvas*)ff[i]->GetObjectUnchecked("can1");		    
	if(!ca[i]){
	   printf("Empty canvas in file %s.\n", outputNames[i].Data());
	   continue;
	}	
        pad[j] = (TPad*)ca[i]->GetListOfPrimitives()->FindObject(Form("can1_%d",j));         	
        if(!pad[j]){
           printf("Empty pad in canvas %s.\n", ca[i]->GetName());
           continue;	     
        }
        g[i] =(TH1D*)pad[j]->FindObject(Form("fTrackHLT_proj_%d",j-1));
	if(!g[i]){
	   printf("Empty histogram for i=%d, file %s.\n", i, outputNames[i].Data());
	   continue;
	}
        
        ov->cd(j);	
        if(i==0) g[i]->Draw();
        else { 
	  g[i]->SetLineColor(i+1); 
	  defineYaxisMax(g[0], g[i]);
	  g[i]->Draw("sames");
	}					 
    }
  }  
  file->Close();    
}

TString cutStudies( TCanvas* can1, TCanvas* can2, TCanvas* can3, TString folder,
//void cutStudies( TCanvas* can1, TCanvas* can2, TCanvas* can3, TString folder,
                    THnSparse* htrackHLT, THnSparse* htrackOFF, TText* hText,
            	    double minEta,   double maxEta,
	    	    int minTrackMult,int maxTrackMult,
            	    double minPt,    double maxPt,
	    	    double minDCAr,  double maxDCAr,
	    	    double minDCAz,  double maxDCAz,
	    	    int minTPCclus,  int maxTPCclus,
	    	    int minITSclus,  int maxITSclus,
		    int vertexStatus
                  )
{
  plotTrackQuantities(can1, htrackHLT, htrackOFF, hText, 
                      minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, minITSclus, maxITSclus, minTrackMult, maxTrackMult, vertexStatus);
  
  plot2D(can2, htrackHLT, minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, minITSclus, maxITSclus, minTrackMult, maxTrackMult);
  
  plot2D(can3, htrackOFF, minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, minITSclus, maxITSclus, minTrackMult, maxTrackMult);
  
  TString cuts = cutsToString(minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, 
                              minITSclus, maxITSclus, minTrackMult, maxTrackMult, vertexStatus);

  can1->SaveAs(folder+"/track_properties_"+cuts+".root");
  can1->SaveAs(folder+"/track_properties_"+cuts+".png");  
  can2->SaveAs(folder+"/HLT_2D_track_correlations_"+cuts+".root");
  can2->SaveAs(folder+"/HLT_2D_track_correlations_"+cuts+".png");
  can3->SaveAs(folder+"/OFF_2D_track_correlations_"+cuts+".root");
  can3->SaveAs(folder+"/OFF_2D_track_correlations_"+cuts+".png");
  
  return folder+"/track_properties_"+cuts+".root";
}

void printStats(TH1D *hlt, TH1D *off){  
  gPad->Update();
  TPaveStats *st1 = (TPaveStats*)hlt->FindObject("stats");
  st1->SetLineColor(0);

  gPad->Update();
  TPaveStats *st2 = (TPaveStats*)off->FindObject("stats");
  st2->SetY2NDC(st1->GetY1NDC()-0.05);
  st2->SetY1NDC(st2->GetY2NDC()-TMath::Abs(st1->GetY1NDC()-st1->GetY2NDC()));
  st2->SetLineColor(0);
  st2->SetTextColor(off->GetLineColor());
  st2->SetFillStyle(0);
  st2->Draw();
}

void defineYaxisMax(TH1D *hlt, TH1D *off){ 
  if(hlt->GetMaximum() > off->GetMaximum()) off->SetMaximum(1.1*hlt->GetMaximum());
  else hlt->SetMaximum(1.1*off->GetMaximum());
}

void printLegend(TLegend *l, TH1D *hlt, TH1D *off){  
  l->SetFillColor(10); 
  l->SetLineColor(10);
  l->AddEntry(hlt, "HLT", "l");
  l->AddEntry(off, "OFF", "l");
  l->Draw("same");
}

//====================== for 1D distributions ===============================//

void plotAid(TCanvas* can, THnSparse* hHLT, THnSparse* hOFF, TText* hText, TH1D *hlt, TH1D *off, TLegend *l, int size,
             double minEta,    double maxEta,
             double minPt,     double maxPt,
	     double minDCAr,   double maxDCAr,
	     double minDCAz,   double maxDCAz,
	     int minTPCclus,   int maxTPCclus,
	     int minITSclus,   int maxITSclus, 
 	     int minTrackMult, int maxTrackMult,
	     int VS
           )
{     
   hHLT->GetAxis(0)->SetRangeUser(minPt,maxPt);
   hHLT->GetAxis(1)->SetRangeUser(minTPCclus,maxTPCclus);
   hHLT->GetAxis(3)->SetRangeUser(minEta, maxEta);
   hHLT->GetAxis(5)->SetRangeUser(minDCAr, maxDCAr);
   hHLT->GetAxis(6)->SetRangeUser(minDCAz, maxDCAz);
   hHLT->GetAxis(10)->SetRangeUser(minITSclus, maxITSclus);
   hHLT->GetAxis(11)->SetRangeUser(minTrackMult, maxTrackMult);
   if(VS!=2) hHLT->GetAxis(12)->SetRangeUser(VS,VS);
   
   hOFF->GetAxis(0)->SetRangeUser(minPt,maxPt);
   hOFF->GetAxis(1)->SetRangeUser(minTPCclus,maxTPCclus);
   hOFF->GetAxis(3)->SetRangeUser(minEta, maxEta);
   hOFF->GetAxis(5)->SetRangeUser(minDCAr, maxDCAr);
   hOFF->GetAxis(6)->SetRangeUser(minDCAz, maxDCAz);
   hOFF->GetAxis(10)->SetRangeUser(minITSclus, maxITSclus);
   hOFF->GetAxis(11)->SetRangeUser(minTrackMult, maxTrackMult);
   if(VS!=2) hOFF->GetAxis(12)->SetRangeUser(VS,VS);
       
   for(int i=0; i<size; i++){
  
      hlt = hHLT->Projection(i); if(!hlt){ printf("plotAid: empty HLT histogram\n"); continue; }
      off = hOFF->Projection(i); if(!off){ printf("plotAid: empty OFF histogram\n"); continue; }
      
      hlt->SetTitle(fix1DTitle(hHLT->Projection(i)->GetTitle())); 
   
      TString s = hlt->GetTitle();      
      if(s.Contains("p_")){ 
         s+=" (GeV/c)";
	 hlt->SetXTitle(s);     
      }
      else if(s.Contains("theta") || s.Contains("phi")){
         s+=" (rad)"; 
	 hlt->SetXTitle(s);
      }
      else if(s.Contains("DCA")){
         s+=" (cm)";
	 hlt->SetXTitle(s);
      }
        
      defineYaxisMax(hlt, off);
      off->SetLineColor(2);
      //off->SetLineStyle(2);
     
      can->cd(i+1);
      hlt->Draw();
      off->Draw("sames");
      printStats(hlt, off);
      
      if(i==0){
         printLegend(l,hlt,off);
      	 TPaveText *pave = new TPaveText(50,500,140,2100);
         pave->SetFillColor(kWhite);
         pave->SetLineColor(kWhite);
	 pave->SetShadowColor(kWhite);
         //TText *t1 = 
	 pave->AddText(itoa(minEta)+" < #eta < "+itoa(maxEta));
	 pave->AddText(itoa(minPt)+" < p_{T} (GeV/c) < "+itoa(maxPt));
	 pave->AddText(itoa(minTrackMult)+" < track mult < "+itoa(maxTrackMult));
	 pave->AddText(itoa(minDCAr)+" < DCAr (cm) < "+itoa(maxDCAr));
	 pave->AddText(itoa(minDCAz)+" < DCAz (cm) < "+itoa(maxDCAz));
	 pave->AddText(itoa(minTPCclus)+" < TPC clusters/track < "+itoa(maxTPCclus));
	 pave->AddText(itoa(minITSclus)+" < ITS clusters/track < "+itoa(maxITSclus));  
	 if(VS!=2) pave->AddText("vertex status: "+itoa(VS));
    	 
    	 pave->Draw();
    	 can->Update();
      }      
   } 
}

//====================== for 2D distributions ===============================//

void plot2D(TCanvas* can, THnSparse* h,
            double minEta,    double maxEta,
            double minPt,     double maxPt,
	    double minDCAr,   double maxDCAr,
	    double minDCAz,   double maxDCAz,
	    int minTPCclus,   int maxTPCclus,
	    int minITSclus,   int maxITSclus, 
 	    int minTrackMult, int maxTrackMult
           )
{
  h->GetAxis(0)->SetRangeUser(minPt,maxPt);
  h->GetAxis(1)->SetRangeUser(minTPCclus,maxTPCclus);
  h->GetAxis(3)->SetRangeUser(minEta, maxEta);
  h->GetAxis(5)->SetRangeUser(minDCAr, maxDCAr);
  h->GetAxis(6)->SetRangeUser(minDCAz, maxDCAz);
  h->GetAxis(10)->SetRangeUser(minITSclus, maxITSclus);
  h->GetAxis(11)->SetRangeUser(minTrackMult, maxTrackMult);
  
  can->cd(1);    
  TH2D *ht = h->Projection(1,0);
  ht->SetTitle(fix2DTitle(h->Projection(1)->GetTitle(), h->Projection(0)->GetTitle()));

  TString s = "";
  s = fix1DTitle(h->Projection(0)->GetTitle())+" (GeV/c)";
  ht->SetXTitle(s);
  ht->Draw("colz");
  
  can->cd(2);
  ht = h->Projection(1,3);
  ht->SetTitle(fix2DTitle(h->Projection(1)->GetTitle(), h->Projection(3)->GetTitle()));
  ht->Draw("colz");
  
  can->cd(3);
  ht = h->Projection(1,5);
  ht->SetTitle(fix2DTitle(h->Projection(1)->GetTitle(), h->Projection(5)->GetTitle()));
  s = fix1DTitle(h->Projection(5)->GetTitle())+" (cm)";
  ht->SetXTitle(s);
  ht->Draw("colz");
  
  can->cd(4);
  ht = h->Projection(1,6);
  ht->SetTitle(fix2DTitle(h->Projection(1)->GetTitle(), h->Projection(6)->GetTitle()));
  s = fix1DTitle(h->Projection(6)->GetTitle())+" (cm)";
  ht->SetXTitle(s);
  ht->Draw("colz");
  
  can->cd(5);
  ht = h->Projection(6,0);
  ht->SetTitle(fix2DTitle(h->Projection(6)->GetTitle(), h->Projection(0)->GetTitle()));
  s = fix1DTitle(h->Projection(0)->GetTitle())+" (GeV/c)";
  ht->SetXTitle(s);
  s = fix1DTitle(h->Projection(6)->GetTitle())+" (cm)";
  ht->SetYTitle(s);
  ht->Draw("colz");
  
  can->cd(6);
  ht = h->Projection(6,3);
  ht->SetTitle(fix2DTitle(h->Projection(6)->GetTitle(), h->Projection(3)->GetTitle()));
  s = fix1DTitle(h->Projection(6)->GetTitle())+" (cm)";
  ht->SetYTitle(s);
  ht->Draw("colz");
  
  can->cd(7);
  ht = h->Projection(3,0);
  ht->SetTitle(fix2DTitle(h->Projection(3)->GetTitle(), h->Projection(0)->GetTitle()));
  s = fix1DTitle(h->Projection(0)->GetTitle())+" (GeV/c)";
  ht->SetXTitle(s);
  ht->Draw("colz");
  
  can->cd(8);
  ht = h->Projection(3,4);
  ht->SetTitle(fix2DTitle(h->Projection(3)->GetTitle(), h->Projection(4)->GetTitle()));
  s = fix1DTitle(h->Projection(4)->GetTitle())+" (rad)";
  ht->SetXTitle(s);
  ht->Draw("colz"); 
}

void plotAid(TCanvas* can, THnSparse* hHLT, THnSparse* hOFF, TText* hText, TH1D *hlt, TH1D *off, TLegend *l, int size){
 
  for(int i=0; i<size; i++){         
      hlt = hHLT->Projection(i);
      off = hOFF->Projection(i); 
      hlt->SetTitle(fix1DTitle(hHLT->Projection(i)->GetTitle()));      
      TString s = hlt->GetTitle();      
      if(s.Contains("primary")){ 
         s+=" (cm)";
	 hlt->SetXTitle(s);     
      }
     
      defineYaxisMax(hlt, off);
      off->SetLineColor(2);
     
      can->cd(i+1);
      hlt->Draw();
      off->Draw("sames");
      printStats(hlt, off);
      
      if(i==0) printLegend(l,hlt,off);
   } 
}

void plotEventQuantities(TCanvas* can, THnSparse* heventHLT, THnSparse* heventOFF, TText* hText){

  TH1D *hlt = NULL;
  TH1D *off = NULL;
 
  TLegend *leg1 = new TLegend(0.6,0.6,0.8,0.8);
 
  plotAid(can, heventHLT, heventOFF, hText, hlt, off, leg1, 6);  
  return;
}

void plotTrackQuantities( TCanvas* can, THnSparse* htrackHLT, THnSparse* htrackOFF, TText* hText,
           	    	  double minEta,    double maxEta,
 	                  double minPt,     double maxPt,
	   	    	  double minDCAr,   double maxDCAr,
	   	    	  double minDCAz,   double maxDCAz,
	   	    	  int minTPCclus,   int maxTPCclus,
	   	    	  int minITSclus,   int maxITSclus, 
	                  int minTrackMult, int maxTrackMult,
			  int VS
                        )
{
  TH1D *hlt = NULL;
  TH1D *off = NULL;
 
  TLegend *leg1 = new TLegend(0.6,0.6,0.8,0.8);
  plotAid(can, htrackHLT, htrackOFF, hText, hlt, off, leg1, 11, 
          minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, minITSclus, maxITSclus, minTrackMult, maxTrackMult, VS);  
  return;
}

TString fix1DTitle(const char* c){
  TString tmp = c;
  tmp.ReplaceAll("projection ","");
  return tmp;   
}

TString fix2DTitle(const char* c1, const char* c2){
  TString tmp = fix1DTitle(c1)+" vs."+fix1DTitle(c2);
  return tmp;
}

TString cutsToString( double minEta,	double maxEta,
 	              double minPt,	double maxPt,
	   	      double minDCAr,	double maxDCAr,
	   	      double minDCAz,	double maxDCAz,
	   	      int minTPCclus,	int maxTPCclus,
	   	      int minITSclus,	int maxITSclus, 
	              int minTrackMult, int maxTrackMult,
		      int VS
	            )
{
  TString cuts = "";
  cuts += "eta"+itoa(minEta)+"_"+itoa(maxEta)+"_";
  cuts += "Pt"+itoa(minPt)+"_"+itoa(maxPt)+"_";
  cuts += "TM"+itoa(minTrackMult)+"_"+itoa(maxTrackMult)+"_";
  cuts += "DCAr"+itoa(minDCAr)+"_"+itoa(maxDCAr)+"_";
  cuts += "DCAz"+itoa(minDCAz)+"_"+itoa(maxDCAz)+"_";
  cuts += "TPCclus"+itoa(minTPCclus)+"_"+itoa(maxTPCclus)+"_";
  cuts += "ITSclus"+itoa(minITSclus)+"_"+itoa(maxITSclus)+"_";
  
  if(VS==2) cuts.Chop(); 
  else cuts += "VS"+itoa(VS); 
  return cuts;
}

TString itoa(int i){
  //stringstream si;
  //si << i;
  //return si.str();
  TString x; x += i; return x;
}

TString itoa(double i){
  //stringstream si;
  //si << i;
  //return si.str();
  TString x; x += i; return x;
}
