// $Id$
/*
 * Drawing macro for reading the THnSparse output of the 
 * HLT/QA/tasks/AliAnalysisTaskHLTCentralBarrel.cxx task.
 * 
 * The cuts are user defined around lines 124-128 as members
 * of an array of structures. Now there are a few predefined sets of cuts.
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
 * In case more than 1 sets of cuts have been used, the macro creates a canvas
 * and overlays the different histograms both for HLT and OFF, if the latter
 * is available.
 * 
 * If the user wants more flexibility adding the TLegend that explains the cuts,
 * the macro HLT/QA/tasks/macros/overlayPlots.C should be used. For the current one
 * the legends of the histograms are somewhat difficult to print at the moment.
 * The text keeps appearing in the wrong location (Kelly, 17.05.2011).
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
#include "TH3D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
using std::endl;
#endif

//---------- forward declerations ---------------//

struct cuts { float minEta; float maxEta; float minPt; float maxPt; float minDCAr; float maxDCAr; float minDCAz; float maxDCAz;
	      int minTPCclus; int maxTPCclus; int minITSclus; int maxITSclus; int vertexStatus; float vertexZ; int minV0cent; int maxV0cent;	 
            };
vector<TString> outputNames;
void cutStudies( TString folder,  THnSparse* htrackHLT, THnSparse* htrackOFF, TText* hText, cuts cut);
void printStats(TH1D *hlt, TH1D *off);
void defineYaxisMax(TH1D *h1, TH1D *h2);
void printLegend(TH1D *hlt, TH1D *off);
void plot2D( THnSparse* h, TText *hText, TString folder, cuts cut );
void plotEventQuantities(THnSparse* heventHLT, THnSparse* heventOFF, TText* hText, TString folder,TString fBeamType);
TString fix1DTitle(const char* c);
TString fix2DTitle(const char* c1, const char* c2);
TString cutsToString( cuts cut );
void plotTrackQuantities( THnSparse* htrackHLT, THnSparse* htrackOFF, TText* hText, TString folder, cuts cut);
void printCuts(cuts cut);	

//------------------------------------------------------------------//		

void drawTHnSparse(TString beamType="p-p", TString inputFile="HLT-OFFLINE-CentralBarrel-comparison.root"){
 
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat("emr");
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

  TString folder = "CentralBarrelTask_";
  folder += hText->GetTitle();
  folder.ReplaceAll(" ",""); 
  folder.ReplaceAll(",","_");    
  gSystem->Exec("mkdir "+folder); // create a folder whose name contains run number and date of run
  
  if(heventHLT->GetEntries()>0 || heventOFF->GetEntries()>0) plotEventQuantities(heventHLT,heventOFF,hText,folder,beamType);
  else {
    if(heventHLT->GetEntries()==0) printf("\nThe HLT event THnSparse contains 0 entries\n");
    if(heventOFF->GetEntries()==0) printf("\nThe OFF event THnSparse contains 0 entries\n");
  }   
   
  cuts p[] = {                                                                                                                                                                               
         // eta    pt      DCAr     DCAz    TPCclus  ITSclus  vtxStatus  |vtxZ|   V0cent                                                                                                
       {-1, 1,  0,   10, -10, 10, -10, 10,  0, 200,   0, 10,      2,       10,   0, 100},                                                                                                    
       {-1, 1,  0.30,10,  -7,  7,  -7,  7,  0, 200,   0, 10,      2,       10,   0, 100},                                                                                                    
       {-1, 1,  0.60,10,  -7,  7,  -7,  7,  0, 200,   0, 10,      2,       10,   0, 100},                                                                                                    
       {-1, 1,  0.90,10,  -7,  7,  -7,  7,  0, 200,   0, 10,      2,       10,   0, 100}                                                                                                     
  };                                                                                                                                                                                         
  const int nCutSets = sizeof(p)/sizeof(cuts);
 
  for(int i=0; i<nCutSets; ++i) cutStudies(folder, htrackHLT, htrackOFF, hText, p[i]); 
  
  // the TString vector outputNames is filled at the end of the cutStudies() call with the name of the output name for every cut.  
  // If there is more than 1 set of cuts, then there will be 2 new canvases created, HLT and OFF respectively with the overlaid track properties
     
  if(outputNames.size()>=2){          
     TString tmp = "overlaid HLT track distributions for ";
     tmp += hText->GetTitle();
     TCanvas *ovHLT = new TCanvas("ovHLT",tmp,1200,800);
     ovHLT->Divide(4,2);
     
     tmp = "overlaid OFF track distributions for ";
     tmp += hText->GetTitle();
     TCanvas *ovOFF = new TCanvas("ovOFF",tmp,1200,800);
     ovOFF->Divide(4,2);

     TCanvas *ca; TFile *ff; TPad *pad; 
     TH1D *hlt[outputNames.size()];
     TH1D *off[outputNames.size()];
          
     for(int j=1; j<9; j++){ // loop over the pads of the canvas "ov" with dimension 3x3     
       for(UInt_t i=0; i<outputNames.size(); i++){ // loop over the files with different sets of cuts
    	      
     	   ff = TFile::Open(outputNames[i].Data());   
     	   if(!ff || ff->IsZombie()){
     	      printf("Non-existent, corrupted or zombie file %s\n", outputNames[i].Data());
     	      continue;
     	   } 
     	   ca  = (TCanvas*)ff->GetObjectUnchecked("can3");	       
     	   if(!ca){
     	      printf("Empty canvas in file %s\n", outputNames[i].Data());
     	      continue;
     	   }	   
     	   pad = (TPad*)ca->GetListOfPrimitives()->FindObject(Form("can3_%d",j));	   
     	   if(!pad){
     	      printf("Empty pad in canvas %s\n", ca->GetName());
     	      continue; 	
     	   }

     	   hlt[i] =(TH1D*)pad->FindObject(Form("fTrackHLT_proj_%d",j-1)); 	   	   	   
     	   if(!hlt[i]){
     	      printf("Empty HLT histogram for cuts i=%d, file %s\n", i, outputNames[i].Data());
     	      continue;
     	   }
	   
	   ovHLT->cd(j); 	       	  
   	   if(i==0){
	      hlt[i]->Draw();
	   }
     	   else { 
     	     hlt[i]->SetLineColor(i+1); 
     	     defineYaxisMax(hlt[0], hlt[i]); 
     	     hlt[i]->Draw("sames");
     	   }
     	   if(i>0) printStats(hlt[i-1], hlt[i]);
	   
  	   off[i] =(TH1D*)pad->FindObject(Form("fTrackOFF_proj_%d",j-1));
	   if(!off[i]){
     	      printf("Empty OFF histogram for cuts i=%d, file %s\n", i, outputNames[i].Data());
     	      continue;
     	   }
	   
	   ovOFF->cd(j);	   
     	   if(i==0){
	      off[i]->SetLineColor(kBlack); 
	      off[i]->Draw();
              TPaveStats *st = (TPaveStats*)off[i]->FindObject("stats"); 
              st->SetTextColor(kBlack);
	      ovOFF->Update();
	   }
     	   else { 
     	     off[i]->SetLineColor(i+1); 
     	     defineYaxisMax(off[0], off[i]); 
     	     off[i]->Draw("sames");
     	   }
	   if(i>0) printStats(off[i-1], off[i]);
	   
     	   ff->Close(); 				    
       } // end of loop over files 
     } // end of loop over canvas pads
     file->Close();  
     ovHLT->Update(); 
     ovHLT->SaveAs(folder+"/overlaid_HLT_track_cuts.root");
     ovHLT->Print(folder+"/overlaid_HLT_track_cuts.png");  
     ovOFF->Update();
     ovOFF->SaveAs(folder+"/overlaid_OFF_track_cuts.root");
     ovOFF->Print(folder+"/overlaid_OFF_track_cuts.png");  
     //delete ovHLT;
     //delete ovOFF;
  } // end if for counter>=2
  return;
}

// ============== main function for filling the track properties, 1D and 2D ================ //

void cutStudies( TString folder, THnSparse* htrackHLT, THnSparse* htrackOFF, TText* hText, cuts cut){
 
  printCuts(cut);
  plotTrackQuantities(htrackHLT, htrackOFF, hText, folder, cut); 
  if(htrackHLT->GetEntries()>0) plot2D(htrackHLT, hText, folder, cut);  
  if(htrackOFF->GetEntries()>0) plot2D(htrackOFF, hText, folder, cut);
   
  TString strcuts = cutsToString(cut);  
  outputNames.push_back(folder+"/track_properties_"+strcuts+".root");
  return;
}

void plotEventQuantities(THnSparse* heventHLT, THnSparse* heventOFF, TText* hText, TString folder, TString fBeamType){

  TString tmp = "vertex event properties for ";
  tmp += hText->GetTitle();
  TCanvas *can1 = new TCanvas("can1",tmp,1200,700);
  can1->Divide(3,2); 
  
  heventHLT->GetAxis(8)->SetRangeUser(1,1); // select events with existing primary vertex	
  heventOFF->GetAxis(8)->SetRangeUser(1,1);
  TH1D *hlt = NULL; 
  TH1D *off = NULL; 
  
  for(int i=0; i<6; i++){ // loop for HLT/OFF primary and SPD vertex xyz
      can1->cd(i+1);
      hlt = heventHLT->Projection(i); if(!hlt){ printf("plotEventQuantities: empty HLT histogram, projection %d\n",i); continue; }
      off = heventOFF->Projection(i); if(!off){ printf("plotEventQuantities: empty OFF histogram, projection %d\n",i); continue; }
      off->SetLineColor(2); 
      hlt->SetTitle(fix1DTitle(heventHLT->Projection(i)->GetTitle())); 
      off->SetTitle(fix1DTitle(heventOFF->Projection(i)->GetTitle()));       
      TString s = hlt->GetTitle();	
      if(s.Contains("primary")){ 
  	 s+=" (cm)";
  	 hlt->SetXTitle(s); 
  	 off->SetXTitle(s);    
      }
      if(off->GetEntries()>0) defineYaxisMax(hlt, off);
      
      if(hlt->GetEntries()>0) hlt->Draw();
      if(off->GetEntries()>0) off->Draw("sames");
      if(hlt->GetEntries()>0 && off->GetEntries()>0) printStats(hlt, off);
      if(off->GetEntries()>0 && i==0 ) printLegend(hlt,off);
  }
  
  tmp = "general event properties for ";
  tmp += hText->GetTitle(); 
  TCanvas *can2 = new TCanvas("can2",tmp,1200,700);
  can2->Divide(3,2);
  
  can2->cd(1); // track multiplicity
  if(fBeamType.Contains("p-p")){
    heventHLT->GetAxis(7)->SetRangeUser(0,1000);
    heventOFF->GetAxis(7)->SetRangeUser(0,1000);
  }
  hlt = heventHLT->Projection(7); 
  off = heventOFF->Projection(7); 
  off->SetLineColor(2); 
  hlt->SetTitle(fix1DTitle(heventHLT->Projection(7)->GetTitle())); 
  off->SetTitle(fix1DTitle(heventOFF->Projection(7)->GetTitle()));	 
  if(off->GetEntries()>0) defineYaxisMax(hlt, off);
  if(hlt->GetEntries()>0) hlt->Draw();
  if(off->GetEntries()>0){ off->Draw("sames"); printLegend(hlt,off); }
  if(hlt->GetEntries()>0 && off->GetEntries()>0) printStats(hlt, off);

  can2->cd(2); // number of contributors
  if(fBeamType.Contains("p-p")){
    heventHLT->GetAxis(6)->SetRangeUser(0,1000);
    heventOFF->GetAxis(6)->SetRangeUser(0,1000);
  }
  hlt = heventHLT->Projection(6); 
  off = heventOFF->Projection(6); 
  off->SetLineColor(2); 
  hlt->SetTitle(fix1DTitle(heventHLT->Projection(6)->GetTitle())); 
  off->SetTitle(fix1DTitle(heventOFF->Projection(6)->GetTitle()));	 
  if(off->GetEntries()>0) defineYaxisMax(hlt, off);
  if(hlt->GetEntries()>0) hlt->Draw();
  if(off->GetEntries()>0) off->Draw("sames");
  if(hlt->GetEntries()>0 && off->GetEntries()>0) printStats(hlt, off);
   
  can2->cd(3); // vertex status
  hlt = heventHLT->Projection(8); 
  off = heventOFF->Projection(8); 
  off->SetLineColor(2); 
  hlt->SetTitle(fix1DTitle(heventHLT->Projection(8)->GetTitle())); 
  off->SetTitle(fix1DTitle(heventOFF->Projection(8)->GetTitle()));	 
  if(off->GetEntries()>0) defineYaxisMax(hlt, off);
  if(hlt->GetEntries()>0) hlt->Draw();
  if(off->GetEntries()>0) off->Draw("sames");
  if(hlt->GetEntries()>0 && off->GetEntries()>0) printStats(hlt, off); 
  
  can2->cd(4); // # of contributors vs. track multiplicity for HLT
  TH2D *h = heventHLT->Projection(6,7);
  TString s1 = fix2DTitle(heventHLT->Projection(6)->GetTitle(), heventHLT->Projection(7)->GetTitle()); s1+=" (HLT)";
  h->SetTitle(s1);
  h->Draw("colz");
  
  can2->cd(5); // # of contributors vs. track multiplicity for OFF
  TH2D *o = heventOFF->Projection(6,7);
  s1 = fix2DTitle(heventOFF->Projection(6)->GetTitle(), heventOFF->Projection(7)->GetTitle()); s1+=" (OFF)";
  o->SetTitle(s1);
  o->Draw("colz");
   
  if(heventHLT->GetNdimensions()==10){
     can2->cd(6);
     off = heventOFF->Projection(9); // V0 centrality, taken from the offline ESD
     off->SetTitle(fix1DTitle(heventOFF->Projection(9)->GetTitle()));
     off->SetLineColor(2);
     off->Draw();
  }  
  
  can1->SaveAs(folder+"/vertex_event_properties.root");
  can1->SaveAs(folder+"/vertex_event_properties.png");
  can2->SaveAs(folder+"/general_event_properties.root");
  can2->SaveAs(folder+"/general_event_properties.png");
  //delete can1;
  //delete can2;
    
  return;
}

void plotTrackQuantities( THnSparse* htrackHLT, THnSparse* htrackOFF, TText* hText, TString folder, cuts cut ){   
  
  htrackHLT->GetAxis(0)->SetRangeUser(cut.minPt,cut.maxPt);
  htrackHLT->GetAxis(1)->SetRangeUser(cut.minTPCclus,cut.maxTPCclus);
  htrackHLT->GetAxis(2)->SetRangeUser(cut.minEta, cut.maxEta);
  htrackHLT->GetAxis(4)->SetRangeUser(cut.minDCAr, cut.maxDCAr);
  htrackHLT->GetAxis(5)->SetRangeUser(cut.minDCAz, cut.maxDCAz);
  htrackHLT->GetAxis(7)->SetRangeUser(cut.minITSclus, cut.maxITSclus);
  if(cut.vertexStatus!=2) htrackHLT->GetAxis(8)->SetRangeUser(cut.vertexStatus, cut.vertexStatus);
  htrackHLT->GetAxis(9)->SetRangeUser(-TMath::Abs(cut.vertexZ), TMath::Abs(cut.vertexZ));
  if(htrackHLT->GetNdimensions()==11) htrackHLT->GetAxis(10)->SetRangeUser(cut.minV0cent, cut.maxV0cent);
  
  htrackOFF->GetAxis(0)->SetRangeUser(cut.minPt,cut.maxPt);
  htrackOFF->GetAxis(1)->SetRangeUser(cut.minTPCclus,cut.maxTPCclus);
  htrackOFF->GetAxis(2)->SetRangeUser(cut.minEta, cut.maxEta);
  htrackOFF->GetAxis(4)->SetRangeUser(cut.minDCAr, cut.maxDCAr);
  htrackOFF->GetAxis(5)->SetRangeUser(cut.minDCAz, cut.maxDCAz);
  htrackOFF->GetAxis(7)->SetRangeUser(cut.minITSclus, cut.maxITSclus);
  if(cut.vertexStatus!=2) htrackOFF->GetAxis(8)->SetRangeUser(cut.vertexStatus,cut.vertexStatus);
  htrackOFF->GetAxis(9)->SetRangeUser(-TMath::Abs(cut.vertexZ), TMath::Abs(cut.vertexZ));
  if(htrackOFF->GetNdimensions()==11) htrackOFF->GetAxis(10)->SetRangeUser(cut.minV0cent, cut.maxV0cent);
  
  TString tmp = "";
  tmp += "track properties for ";
  tmp+=hText->GetTitle();     
  TCanvas *can3 = new TCanvas("can3",tmp,1600,1000); 
  can3->Divide(4,2); 
  // the first 9 track related variables filled in the THnSparse
  //  0	    1     2    3     4      5      6     7   
  // pt  TPCcl  eta   phi   DCAr  DCAz charge  ITScl 
  
  TH1D *hlt = NULL; 
  TH1D *off = NULL;

  TLegend *leg = new TLegend(0.25,0.2,0.85,0.85);

  for(int i=0; i<8; i++){  
    hlt = htrackHLT->Projection(i); if(!hlt){ printf("plotTrackQuantities: empty HLT histogram for projection %d\n",i); continue; }
    off = htrackOFF->Projection(i); if(!off){ printf("plotTrackQuantities: empty OFF histogram for projection %d\n",i); continue; }     
    hlt->SetTitle(fix1DTitle(htrackHLT->Projection(i)->GetTitle())); 
    off->SetTitle(fix1DTitle(htrackOFF->Projection(i)->GetTitle())); // is necessary in cases where only offline data is available
  
    TString s = hlt->GetTitle();      
    if(s.Contains("p_")){ 
      s+=" (GeV/c)";
      hlt->SetXTitle(s);     
    }
    else if(s.Contains("phi")){
      s+=" (rad)"; 
      hlt->SetXTitle(s);
    }
    else if(s.Contains("DCA")){
      s+=" (cm)";
      hlt->SetXTitle(s);
    }
       
    if(off->GetEntries()>0) defineYaxisMax(hlt, off);
    off->SetLineColor(2);
   
    can3->cd(i+1);
    if(hlt->GetEntries()>0) hlt->Draw();
    if(off->GetEntries()>0) off->Draw("sames");
    if(hlt->GetEntries()>0 && off->GetEntries()>0) printStats(hlt, off);
    if(off->GetEntries()>0 && i==0 ) printLegend(hlt,off);
     
    if(hlt->GetEntries()>0 && i==0){
      leg->SetFillColor(kWhite);
      leg->SetLineColor(kWhite);
      leg->SetShadowColor(kWhite);
      s=""; s+=cut.minEta; s+="<#eta< "; s+=cut.maxEta;leg->AddEntry((TObject*)0, s, "");  
      s=""; s+=cut.minPt; s.Resize(4); s+="<p_{T} (GeV/c)<"; s+=cut.maxPt;  leg->AddEntry((TObject*)0, s, ""); 
      s=""; s+=cut.minDCAr; s+="<DCAr (cm)<"; s+=cut.maxDCAr; leg->AddEntry((TObject*)0, s, "");  
      s=""; s+=cut.minDCAz; s+="<DCAz (cm)<"; s+=cut.maxDCAz; leg->AddEntry((TObject*)0, s, "");  
      if(cut.minTPCclus!=0 && cut.maxTPCclus!=200){ s=""; s+=cut.minTPCclus; s+="<TPC cls/tr<"; s+=cut.maxTPCclus; leg->AddEntry((TObject*)0, s, "");}
      if(cut.minITSclus!=0 && cut.maxITSclus!=10){ s=""; s+=cut.minITSclus; s+="<ITS cls/tr<"; s+=cut.maxITSclus; leg->AddEntry((TObject*)0, s, "");}
      s=""; s+="|vertexZ| (cm)<"; s+=TMath::Abs(cut.vertexZ); leg->AddEntry((TObject*)0, s, ""); 
      if(cut.vertexStatus!=2) { s=""; s+="vertex status "; s+=cut.vertexStatus; leg->AddEntry((TObject*)0, s, ""); }
      if(htrackHLT->GetNdimensions()==11) { s=""; s+=cut.minV0cent; s+=" < V0 centr < "; s+=cut.maxV0cent; leg->AddEntry((TObject*)0, s, "");  }
      
      leg->Draw();
      
      can3->Update();
	}      
  } 
  TString strcuts = cutsToString(cut);
  
  can3->SaveAs(folder+"/track_properties_"+strcuts+".root");
  can3->SaveAs(folder+"/track_properties_"+strcuts+".png");
  delete can3;
  delete leg;

  return;
}

//====================== for 2D track distributions ===============================//

void plot2D( THnSparse* h, TText *hText, TString folder, cuts cut ){

  h->GetAxis(0)->SetRangeUser(cut.minPt,cut.maxPt);
  h->GetAxis(1)->SetRangeUser(cut.minTPCclus,cut.maxTPCclus);
  h->GetAxis(2)->SetRangeUser(cut.minEta, cut.maxEta);
  h->GetAxis(4)->SetRangeUser(cut.minDCAr, cut.maxDCAr);
  h->GetAxis(5)->SetRangeUser(cut.minDCAz, cut.maxDCAz);
  h->GetAxis(7)->SetRangeUser(cut.minITSclus,cut. maxITSclus);
  if(cut.vertexStatus!=2) h->GetAxis(8)->SetRangeUser(cut.vertexStatus,cut.vertexStatus);
  h->GetAxis(9)->SetRangeUser(-TMath::Abs(cut.vertexZ), TMath::Abs(cut.vertexZ));
  if(h->GetNdimensions()==11) h->GetAxis(10)->SetRangeUser(cut.minV0cent, cut.maxV0cent);
  
  TString name = h->GetName();
  TString tmp = "";
  if(name.Contains("HLT")) tmp = "HLT 2-D track distributions for ";
  else tmp = "OFF 2-D track distributions for ";

  tmp += hText->GetTitle();
  TCanvas *can4 = new TCanvas("can4",tmp,1200,800);
  can4->Divide(4,2);

  can4->cd(1);    
  TH2D *ht = h->Projection(1,0); // TPC clusters/track vs. pt
  ht->SetTitle(fix2DTitle(h->Projection(1)->GetTitle(), h->Projection(0)->GetTitle()));

  TString s = "";
  s = fix1DTitle(h->Projection(0)->GetTitle())+" (GeV/c)";
  ht->SetXTitle(s);
  ht->Draw("colz");
  
  can4->cd(2);
  ht = h->Projection(1,2); // TPC clusters/track vs. eta
  ht->SetTitle(fix2DTitle(h->Projection(1)->GetTitle(), h->Projection(2)->GetTitle()));
  ht->Draw("colz");
  
  can4->cd(3);
  ht = h->Projection(1,4); // TPC clusters/track vs. DCAr
  ht->SetTitle(fix2DTitle(h->Projection(1)->GetTitle(), h->Projection(4)->GetTitle()));
  s = fix1DTitle(h->Projection(4)->GetTitle())+" (cm)";
  ht->SetXTitle(s);
  ht->Draw("colz");
  
  can4->cd(4);
  ht = h->Projection(1,5); // TPC clusters/track vs. DCAz
  ht->SetTitle(fix2DTitle(h->Projection(1)->GetTitle(), h->Projection(5)->GetTitle()));
  s = fix1DTitle(h->Projection(5)->GetTitle())+" (cm)";
  ht->SetXTitle(s);
  ht->Draw("colz");
  
  can4->cd(5);
  ht = h->Projection(4,0); // DCAr vs. pt
  ht->SetTitle(fix2DTitle(h->Projection(4)->GetTitle(), h->Projection(0)->GetTitle()));
  s = fix1DTitle(h->Projection(0)->GetTitle())+" (GeV/c)";
  ht->SetXTitle(s);
  s = fix1DTitle(h->Projection(4)->GetTitle())+" (cm)";
  ht->SetYTitle(s);
  ht->Draw("colz");
  
  can4->cd(6);
  ht = h->Projection(5,3); // DCAz vs. pt
  ht->SetTitle(fix2DTitle(h->Projection(5)->GetTitle(), h->Projection(3)->GetTitle()));
  s = fix1DTitle(h->Projection(5)->GetTitle())+" (cm)";
  ht->SetYTitle(s);
  ht->Draw("colz");
  
  can4->cd(7);
  ht = h->Projection(2,0); // eta vs. pt
  ht->SetTitle(fix2DTitle(h->Projection(2)->GetTitle(), h->Projection(0)->GetTitle()));
  s = fix1DTitle(h->Projection(0)->GetTitle())+" (GeV/c)";
  ht->SetXTitle(s);
  ht->Draw("colz");
  
  can4->cd(8);
  ht = h->Projection(2,3); // eta vs. phi
  ht->SetTitle(fix2DTitle(h->Projection(2)->GetTitle(), h->Projection(3)->GetTitle()));
  s = fix1DTitle(h->Projection(3)->GetTitle())+" (rad)";
  ht->SetXTitle(s);
  ht->Draw("colz");
 
  TString strcuts = cutsToString(cut);
  if(name.Contains("HLT")){
     can4->SaveAs(folder+"/HLT_2D_track_correlations_"+strcuts+".root"); 
     can4->SaveAs(folder+"/HLT_2D_track_correlations_"+strcuts+".png"); 
  } else {
     can4->SaveAs(folder+"/OFF_2D_track_correlations_"+strcuts+".root"); 
     can4->SaveAs(folder+"/OFF_2D_track_correlations_"+strcuts+".png"); 
  }
  //delete can4;
  return;
}

TString cutsToString( cuts cut ){
 
  TString strcuts = "";
  char s[300]; sprintf(s, "eta%2g_%2g_Pt%2g_%2g_DCAr%2g_%2g_DCAz%2g_%2g_TPCclus%d_%d_ITSclus%d_%d_Zvertex%2g_cent%d_%d", 
                           cut.minEta,cut.maxEta,cut.minPt,cut.maxPt,cut.minDCAr,cut.maxDCAr,cut.minDCAz,cut.maxDCAz,cut.minTPCclus,cut.maxTPCclus,cut.minITSclus,
			   cut.maxITSclus,cut.vertexZ,cut.minV0cent,cut.maxV0cent);
  strcuts = s; strcuts.ReplaceAll(" ","");
 
  if(cut.vertexStatus!=2){
    char v[10]; 
    sprintf(v, "_VS%d",cut.vertexStatus);
    strcuts+=v;
  }
  return strcuts;
}

void printStats(TH1D *hlt, TH1D *off){  
  gPad->Update();
  TPaveStats *st1 = (TPaveStats*)hlt->FindObject("stats");
  st1->SetLineColor(0);
  //st1->SetTextSize(7);
  //st1->SetTextFont(8);

  gPad->Update();
  TPaveStats *st2 = (TPaveStats*)off->FindObject("stats");
  st2->SetY2NDC(st1->GetY1NDC()-0.05);
  st2->SetY1NDC(st2->GetY2NDC()-TMath::Abs(st1->GetY1NDC()-st1->GetY2NDC()));
  st2->SetLineColor(0);
  st2->SetTextColor(off->GetLineColor());
  st2->SetFillStyle(0);
  //st2->SetTextSize(7);
  //st2->SetTextFont(8);
  st2->Draw();
  return;
}

void defineYaxisMax(TH1D *h1, TH1D *h2){ 
  //Y axis
  if(h1->GetMaximum() > h2->GetMaximum()) h2->SetMaximum(1.1*h1->GetMaximum());
  else h1->SetMaximum(1.1*h2->GetMaximum());
  
  h1->SetMinimum(0);
  h2->SetMinimum(0);
 
  // X axis  
  double xmin, xmax;  
  if(h1->GetBinLowEdge(1) > h2->GetBinLowEdge(1)) xmin = h1->GetBinLowEdge(1);
  else xmin = h2->GetBinLowEdge(1);
  if(h1->GetBinLowEdge(h1->GetNbinsX()+1) > h2->GetBinLowEdge(h1->GetNbinsX()+1)) xmax = h1->GetBinLowEdge(h1->GetNbinsX()+1);
  else xmax = h2->GetBinLowEdge(h2->GetNbinsX()+1);
  
  h2->SetAxisRange(xmin, xmax, "X");  
  return;
}

void printLegend(TH1D *hlt, TH1D *off){  
  TLegend *l = new TLegend(0.68,0.16,0.88,0.36);
  l->SetFillColor(10); 
  l->SetLineColor(10);
  l->AddEntry(hlt, "HLT", "l");
  l->AddEntry(off, "OFF", "l");
  l->Draw("same");
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

void printCuts(cuts cut){
  printf("\n %2g<eta<%2g, %2g<pt (GeV/c)<%2g, %2g<DCAr (cm)<%2g, %2g<DCAz (cm)<%2g, %d<TPCclus<%d, %d<ITSclus<%d, vertex status %d, |vertexZ| (cm)<%2g, %d<centr (%)<%d \n\n",cut.minEta, cut.maxEta, cut.minPt, cut.maxPt, cut.minDCAr, cut.maxDCAr, cut.minDCAz, cut.maxDCAz, cut.minTPCclus, cut.maxTPCclus, cut.minITSclus, cut.maxITSclus, cut.vertexStatus, cut.vertexZ, cut.minV0cent, cut.maxV0cent );
}
