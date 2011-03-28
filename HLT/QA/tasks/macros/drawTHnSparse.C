// $Id$
/*
 * Drawing macro for reading the THnSparse output of the 
 * HLT/QA/tasks/AliAnalysisTaskHLTCentralBarrel.cxx task.
 * 
 * The cuts are user defined in lines 188-190 as arguments of the 
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
 * root[1] drawTHnSparse("HLT-OFFLINE-CentralBarrel-comparison.root")
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
#include <cstdlib>
using std::endl;
#endif

//---------- forward declerations ---------------//

TString cutStudies( TCanvas* can1, TCanvas* can2, TCanvas* can3, TString folder,
                    THnSparse* htrackHLT, THnSparse* htrackOFF, TText* hText,
            	    double minEta,   double maxEta,
	    	    int minTrackMult,int maxTrackMult,
            	    double minPt,    double maxPt,
	    	    double minDCAr,  double maxDCAr,
	    	    double minDCAz,  double maxDCAz,
	    	    int minTPCclus,  int maxTPCclus,
	    	    int minITSclus,  int maxITSclus,
		    int vs,          float vz,
		    float minCent,   float maxCent
                  );
void printStats(TH1D *hlt, TH1D *off);
void defineYaxisMax(TH1D *hlt, TH1D *off);
void printLegend(TLegend *l, TH1D *hlt, TH1D *off);
void plotAid(TCanvas* can, THnSparse* hHLT, THnSparse* hOFF, TText* hText, TH1D *hlt, TH1D *off, TLegend *l);	   
void plot2D( TCanvas* can,     THnSparse* h,
             double minEta,    double maxEta,
             double minPt,     double maxPt,
	     double minDCAr,   double maxDCAr,
	     double minDCAz,   double maxDCAz,
	     int minTPCclus,   int maxTPCclus,
	     int minITSclus,   int maxITSclus, 
 	     int minTrackMult, int maxTrackMult, 
	     int vs,           float vz,
	     float minCent,    float maxCent
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
		      int vs,           float vz
		      //float minCent,    float maxCent
	            );
void plotTrackQuantities( TCanvas* can,     THnSparse* htrackHLT, THnSparse* htrackOFF, TText* hText,
           	    	  double minEta,    double maxEta,
 	                  double minPt,     double maxPt,
	   	    	  double minDCAr,   double maxDCAr,
	   	    	  double minDCAz,   double maxDCAz,
	   	    	  int minTPCclus,   int maxTPCclus,
	   	    	  int minITSclus,   int maxITSclus, 
	                  int minTrackMult, int maxTrackMult,
			  int vs,           float vz,
			  float minCent,    float maxCent
                        );
			
vector<TString> outputNames;

//------------------------------------------------------------------//		
			
void drawTHnSparse(TString inputFile="HLT-OFFLINE-CentralBarrel-comparison.root"){
 
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(10);
  gStyle->SetTitleX(gStyle->GetPadLeftMargin());
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
      return;
  }
  THnSparseF *heventOFF = (THnSparseF*)list->FindObject("fEventOFF");  
  if(!heventOFF){
      printf("Error: There is no OFF THnSparse object in file %s\n", inputFile.Data());
      return;
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
  
  
  TString folder = hText->GetTitle();
  folder.ReplaceAll(" ",""); 
  folder.Remove(0,3);
  folder.Remove(6,15);
   
  gSystem->Exec("mkdir "+folder);
  TString canvasTitle ="";
  
  if(heventHLT->GetEntries()>0 || heventOFF->GetEntries()>0){ 
     // create and fill the canvas only if the respective THnSparse with event properties has been filled     
     canvasTitle = "event properties for ";
     canvasTitle+=hText->GetTitle();     
     TCanvas *can0 = new TCanvas("can0",canvasTitle, 900,600); 
     can0->Divide(4,2); 
     plotEventQuantities(can0,heventHLT,heventOFF,hText);
     can0->SaveAs(folder+"/event_properties.root");
     can0->SaveAs(folder+"/event_properties.png");
  }
  
  if(htrackHLT->GetEntries()>0 || htrackOFF->GetEntries()>0){
     // create and fill the canvas only if the respective THnSparse with track properties has been filled
     canvasTitle = "track properties for ";
     canvasTitle+=hText->GetTitle();     
     TCanvas *can1 = new TCanvas("can1",canvasTitle,1100,900);                  can1->Divide(3,3); // the first 9 variables filled in the THnSparse
     
     canvasTitle = "2-D HLT track distributions for ";
     canvasTitle+=hText->GetTitle();     
     TCanvas *can2 = new TCanvas("can2",canvasTitle,1200,800); can2->Divide(4,2);
   
     canvasTitle = "2-D OFF track distributions for ";
     canvasTitle+=hText->GetTitle();     
     TCanvas *can3 = new TCanvas("can3",canvasTitle,1200,800); can3->Divide(4,2);
     
     int counter=0; 
     // counts how many times the function cutStudies() is called
     // if more than once, then it creates and fills the canvas with the overlapping hlt distributions for the various sets of cuts
     
     TString s = "";                                                       // eta   mult    pt     DCAr     DCAz    TPCclus  ITSclus  vertexStatus   vertexZ    
     s = cutStudies(can1, can2, can3, folder, htrackHLT, htrackOFF, hText, -1, 1, 0, 2000, 0,  5, -10, 10, -10, 10,  0, 200, 0, 6,       2,	     10, 0, 100); outputNames.push_back(s); counter++;
     	     
     if(counter>=2){     
        
	canvasTitle = "overlaid HLT track distributions for ";
        canvasTitle+=hText->GetTitle();
        TCanvas *ov = new TCanvas("ov",canvasTitle,1100,900);
   	ov->Divide(3,3);

	TCanvas *ca; TFile *ff;	TPad *pad; 
   	TH1D *g[outputNames.size()];
   	
   	for(int j=1; j<10; j++){// loop over the pads of the canvas "ov" with dimension 3x3     
   	  for(UInt_t i=0; i<outputNames.size(); i++){   // loop over the files with different sets of cuts
   		 
	      ff = TFile::Open(outputNames[i].Data());   
   	      if(!ff || ff->IsZombie()){
   		 printf("Non-existent, corrupted or zombie file %s\n", outputNames[i].Data());
   		 return;
   	      } 
   	      ca  = (TCanvas*)ff->GetObjectUnchecked("can1");		  
   	      if(!ca){
   		 printf("Empty canvas in file %s.\n", outputNames[i].Data());
   		 continue;
   	      }       
   	      pad = (TPad*)ca->GetListOfPrimitives()->FindObject(Form("can1_%d",j));	      
   	      if(!pad){
   	         printf("Empty pad in canvas %s.\n", ca->GetName());
   		 continue;	   
   	      }
   	      g[i] =(TH1D*)pad->FindObject(Form("fTrackHLT_proj_%d",j-1));
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
   	      ff->Close();				       
   	  }
   	}  
     }
  } 
  file->Close();  
}

// ============== main function for filling the track properties, 1D and 2D ================ //

TString cutStudies( TCanvas* can1, TCanvas* can2, TCanvas* can3, TString folder,
                    THnSparse* htrackHLT, THnSparse* htrackOFF, TText* hText,
            	    double minEta,   double maxEta,
	    	    int minTrackMult,int maxTrackMult,
            	    double minPt,    double maxPt,
	    	    double minDCAr,  double maxDCAr,
	    	    double minDCAz,  double maxDCAz,
	    	    int minTPCclus,  int maxTPCclus,
	    	    int minITSclus,  int maxITSclus,
		    int vs,          float vz,
		    float minCent,   float maxCent
                  ){
  plotTrackQuantities(can1, htrackHLT, htrackOFF, hText, 
                      minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, minITSclus, maxITSclus, minTrackMult, maxTrackMult, vs, vz,
		      minCent, maxCent);
  
  plot2D(can2, htrackHLT, minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, minITSclus, maxITSclus, minTrackMult, maxTrackMult, vs, vz,
         minCent, maxCent);  
  plot2D(can3, htrackOFF, minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, minITSclus, maxITSclus, minTrackMult, maxTrackMult, vs, vz,
         minCent, maxCent);
  
  TString cuts = cutsToString(minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, 
                              minITSclus, maxITSclus, minTrackMult, maxTrackMult, vs, vz);

  can1->SaveAs(folder+"/track_properties_"+cuts+".root");
  can1->SaveAs(folder+"/track_properties_"+cuts+".png");  
  can2->SaveAs(folder+"/HLT_2D_track_correlations_"+cuts+".root");
  can2->SaveAs(folder+"/HLT_2D_track_correlations_"+cuts+".png");
  can3->SaveAs(folder+"/OFF_2D_track_correlations_"+cuts+".root");
  can3->SaveAs(folder+"/OFF_2D_track_correlations_"+cuts+".png");
  
  return folder+"/track_properties_"+cuts+".root";
}

//====================== for 1D track distributions ===============================//

void plotAid(TCanvas* can, THnSparse* hHLT, THnSparse* hOFF, TText* /*hText*/, TH1D *hlt, TH1D *off, TLegend *l, int size,
             double minEta,    double maxEta,
             double minPt,     double maxPt,
	     double minDCAr,   double maxDCAr,
	     double minDCAz,   double maxDCAz,
	     int minTPCclus,   int maxTPCclus,
	     int minITSclus,   int maxITSclus, 
 	     int minTrackMult, int maxTrackMult,
	     int vs,           float vz,
	     float minCent,    float maxCent
           )
{     
   hHLT->GetAxis(0)->SetRangeUser(minPt,maxPt);
   hHLT->GetAxis(1)->SetRangeUser(minTPCclus,maxTPCclus);
   hHLT->GetAxis(3)->SetRangeUser(minEta, maxEta);
   hHLT->GetAxis(5)->SetRangeUser(minDCAr, maxDCAr);
   hHLT->GetAxis(6)->SetRangeUser(minDCAz, maxDCAz);
   hHLT->GetAxis(8)->SetRangeUser(minITSclus, maxITSclus);
   hHLT->GetAxis(9)->SetRangeUser(minTrackMult, maxTrackMult);
   if(vs!=2) hHLT->GetAxis(10)->SetRangeUser(vs,vs);
   hHLT->GetAxis(11)->SetRangeUser(-TMath::Abs(vz), TMath::Abs(vz));
   if(hHLT->GetNdimensions()==13) hHLT->GetAxis(12)->SetRangeUser(minCent, maxCent);
   
   hOFF->GetAxis(0)->SetRangeUser(minPt,maxPt);
   hOFF->GetAxis(1)->SetRangeUser(minTPCclus,maxTPCclus);
   hOFF->GetAxis(3)->SetRangeUser(minEta, maxEta);
   hOFF->GetAxis(5)->SetRangeUser(minDCAr, maxDCAr);
   hOFF->GetAxis(6)->SetRangeUser(minDCAz, maxDCAz);
   hOFF->GetAxis(8)->SetRangeUser(minITSclus, maxITSclus);
   hOFF->GetAxis(9)->SetRangeUser(minTrackMult, maxTrackMult);
   if(vs!=2) hOFF->GetAxis(10)->SetRangeUser(vs,vs);
   hOFF->GetAxis(11)->SetRangeUser(-TMath::Abs(vz), TMath::Abs(vz));
   if(hOFF->GetNdimensions()==13) hOFF->GetAxis(12)->SetRangeUser(minCent, maxCent);
       
   for(int i=0; i<size; i++){
  
      hlt = hHLT->Projection(i); if(!hlt){ printf("plotAid for track properties: empty HLT histogram\n"); continue; }
      off = hOFF->Projection(i); if(!off){ printf("plotAid for track properties: empty OFF histogram\n"); continue; }
      
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
    
      can->cd(i+1);
      hlt->Draw();
      off->Draw("sames");
      printStats(hlt, off);
      
      if(i==0){
         printLegend(l,hlt,off);
      	 TPaveText *pave = new TPaveText(25,5600,140,36400);
         pave->SetFillColor(kWhite);
         pave->SetLineColor(kWhite);
	 pave->SetShadowColor(kWhite);
	 s=""; s+=minEta; s+=" < eta < "; s+=maxEta; pave->AddText(s);
	 s=""; s+=minPt; s+=" < pt (GeV/c) < "; s+=maxPt; pave->AddText(s);
	 s=""; s+=minTrackMult; s+=" < track mult < "; s+=maxTrackMult; pave->AddText(s);
	 s=""; s+=minDCAr; s+=" < DCAr (cm) < "; s+=maxDCAr; pave->AddText(s);
	 s=""; s+=minDCAz; s+=" < DCAz (cm) < "; s+=maxDCAz; pave->AddText(s);
	 s=""; s+=minTPCclus; s+=" < TPC clusters/track < "; s+=maxTPCclus; pave->AddText(s);
	 s=""; s+=minITSclus; s+=" < ITS clusters/track < "; s+=maxITSclus; pave->AddText(s);
	 if(vs!=2) { s=""; s+="vertex status "; s+=vs; pave->AddText(s); }
    	 pave->Draw();
    	 can->Update();
      }      
   } 
}

//====================== for 2D track distributions ===============================//

void plot2D(TCanvas* can,     THnSparse* h,
            double minEta,    double maxEta,
            double minPt,     double maxPt,
	    double minDCAr,   double maxDCAr,
	    double minDCAz,   double maxDCAz,
	    int minTPCclus,   int maxTPCclus,
	    int minITSclus,   int maxITSclus, 
 	    int minTrackMult, int maxTrackMult,
	    int vs,           float vz,
	    float minCent,    float maxCent
           )
{
  h->GetAxis(0)->SetRangeUser(minPt,maxPt);
  h->GetAxis(1)->SetRangeUser(minTPCclus,maxTPCclus);
  h->GetAxis(3)->SetRangeUser(minEta, maxEta);
  h->GetAxis(5)->SetRangeUser(minDCAr, maxDCAr);
  h->GetAxis(6)->SetRangeUser(minDCAz, maxDCAz);
  h->GetAxis(8)->SetRangeUser(minITSclus, maxITSclus);
  h->GetAxis(9)->SetRangeUser(minTrackMult, maxTrackMult);
  if(vs!=2) h->GetAxis(10)->SetRangeUser(vs,vs);
  h->GetAxis(11)->SetRangeUser(-TMath::Abs(vz), TMath::Abs(vz));
  if(h->GetNdimensions()==13) h->GetAxis(12)->SetRangeUser(minCent, maxCent);
  
  can->cd(1);    
  TH2D *ht = h->Projection(1,0); // TPC clusters/track vs. pt
  ht->SetTitle(fix2DTitle(h->Projection(1)->GetTitle(), h->Projection(0)->GetTitle()));

  TString s = "";
  s = fix1DTitle(h->Projection(0)->GetTitle())+" (GeV/c)";
  ht->SetXTitle(s);
  ht->Draw("colz");
  
  can->cd(2);
  ht = h->Projection(1,3); // TPC clusters/track vs. eta
  ht->SetTitle(fix2DTitle(h->Projection(1)->GetTitle(), h->Projection(3)->GetTitle()));
  ht->Draw("colz");
  
  can->cd(3);
  ht = h->Projection(1,5); // TPC clusters/track vs. DCAr
  ht->SetTitle(fix2DTitle(h->Projection(1)->GetTitle(), h->Projection(5)->GetTitle()));
  s = fix1DTitle(h->Projection(5)->GetTitle())+" (cm)";
  ht->SetXTitle(s);
  ht->Draw("colz");
  
  can->cd(4);
  ht = h->Projection(1,6); // TPC clusters/track vs. DCAz
  ht->SetTitle(fix2DTitle(h->Projection(1)->GetTitle(), h->Projection(6)->GetTitle()));
  s = fix1DTitle(h->Projection(6)->GetTitle())+" (cm)";
  ht->SetXTitle(s);
  ht->Draw("colz");
  
  can->cd(5);
  ht = h->Projection(5,0); // DCAr vs. pt
  ht->SetTitle(fix2DTitle(h->Projection(5)->GetTitle(), h->Projection(0)->GetTitle()));
  s = fix1DTitle(h->Projection(0)->GetTitle())+" (GeV/c)";
  ht->SetXTitle(s);
  s = fix1DTitle(h->Projection(5)->GetTitle())+" (cm)";
  ht->SetYTitle(s);
  ht->Draw("colz");
  
  can->cd(6);
  ht = h->Projection(6,3); // DCAz vs. pt
  ht->SetTitle(fix2DTitle(h->Projection(6)->GetTitle(), h->Projection(3)->GetTitle()));
  s = fix1DTitle(h->Projection(6)->GetTitle())+" (cm)";
  ht->SetYTitle(s);
  ht->Draw("colz");
  
  can->cd(7);
  ht = h->Projection(3,0); // eta vs. pt
  ht->SetTitle(fix2DTitle(h->Projection(3)->GetTitle(), h->Projection(0)->GetTitle()));
  s = fix1DTitle(h->Projection(0)->GetTitle())+" (GeV/c)";
  ht->SetXTitle(s);
  ht->Draw("colz");
  
  can->cd(8);
  ht = h->Projection(3,4); // eta vs. phi
  ht->SetTitle(fix2DTitle(h->Projection(3)->GetTitle(), h->Projection(4)->GetTitle()));
  s = fix1DTitle(h->Projection(4)->GetTitle())+" (rad)";
  ht->SetXTitle(s);
  ht->Draw("colz"); 
}

// ================== for event properties ==================== //

void plotAid(TCanvas* can, THnSparse* hHLT, THnSparse* hOFF, TText* /*hText*/, TH1D *hlt, TH1D *off, TLegend *l){
      
  int size=0;
  if(hHLT->GetNdimensions()==6 || hOFF->GetNdimensions()==6) size = 6;
  else size = 7; 
  for(int i=0; i<size; i++){ 
  
     hHLT->GetAxis(5)->SetRangeUser(1,1); // select events with found primary vertex       
     hOFF->GetAxis(5)->SetRangeUser(1,1);
     
     hlt = hHLT->Projection(i); if(!hlt){ printf("plotAid for event properties: empty HLT histogram\n"); continue; }
     off = hOFF->Projection(i); if(!off){ printf("plotAid for event properties: empty OFF histogram\n"); continue; }
     hlt->SetTitle(fix1DTitle(hHLT->Projection(i)->GetTitle()));      
     off->SetTitle(fix1DTitle(hOFF->Projection(i)->GetTitle()));      
     TString s = hlt->GetTitle();      
     if(s.Contains("primary")){ 
  	s+=" (cm)";
        hlt->SetXTitle(s); 
        off->SetXTitle(s);    
     }
    
     defineYaxisMax(hlt, off);
     off->SetLineColor(2);
    
     can->cd(i+1);
     hlt->Draw();
     off->Draw("sames");
     printStats(hlt, off);
     
     if(i==0) printLegend(l,hlt,off);
  } 
  can->cd(7);
  TH2D *h = hHLT->Projection(3,4);
  TString s1 = fix2DTitle(hHLT->Projection(3)->GetTitle(), hHLT->Projection(4)->GetTitle()); s1+=" (HLT)";
  h->SetTitle(s1);
  h->Draw("colz");
 
  can->cd(8);
  TH2D *o = hOFF->Projection(3,4);
  TString s2 = fix2DTitle(hOFF->Projection(3)->GetTitle(), hOFF->Projection(4)->GetTitle()); s2+=" (OFF)";
  o->SetTitle(s2);
  o->Draw("colz");
}

void plotEventQuantities(TCanvas* can, THnSparse* heventHLT, THnSparse* heventOFF, TText* hText){

  TH1D *hlt = NULL;
  TH1D *off = NULL;
 
  TLegend *leg1 = new TLegend(0.6,0.6,0.8,0.8);
 
  plotAid(can, heventHLT, heventOFF, hText, hlt, off, leg1);  
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
			  int vs,           float vz,
			  float minCent,    float maxCent
                        )
{
  TH1D *hlt = NULL;
  TH1D *off = NULL;
 
  TLegend *leg1 = new TLegend(0.6,0.6,0.8,0.8);
    
  plotAid(can, htrackHLT, htrackOFF, hText, hlt, off, leg1, 9, 
          minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, minITSclus, maxITSclus, minTrackMult, maxTrackMult, vs, vz, minCent, maxCent);  
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
		      int vs,           float vz
//		      float minCent,    float maxCent
	            )
{
  TString cuts = "";
  char s[300]; sprintf(s, "eta%2g_%2g_Pt%2g_%2g_TM%d_%d_DCAr%2g_%2g_DCAz%2g_%2g_TPCclus%d_%d_ITSclus%d_%d_Zvertex%2g", 
                           minEta,maxEta,minPt,maxPt,minTrackMult,maxTrackMult,minDCAr,maxDCAr,minDCAz,maxDCAz,minTPCclus,maxTPCclus,minITSclus,maxITSclus,vz);
  cuts = s; cuts.ReplaceAll(" ","");
 
  if(vs!=2){
    char v[10]; 
    sprintf(v, "_VS%d",vs);
    cuts+=v;
  }
  return cuts;
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
