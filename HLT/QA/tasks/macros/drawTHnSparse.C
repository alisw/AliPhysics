// $Id$
/*
 * Drawing macro for reading the THnSparse output of the 
 * HLT/QA/tasks/AliAnalysisTaskHLTCentralBarrel.cxx task.
 * 
 * The cuts are user defined as arguments of the function.
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
 * @ingroup alihlt_qa
 * @author Kalliopi.Kanaki@ift.uib.no 
 */

void drawTHnSparse( TString inputFile,
            	    double minEta=-2,   double maxEta=2,
	    	    int minTrackMult=0, int maxTrackMult=20000,
            	    double minPt=0,     double maxPt=200,
	    	    double minDCAr=-80, double maxDCAr=80,
	    	    double minDCAz=-80, double maxDCAz=80,
	    	    int minTPCclus=0,   int maxTPCclus=200,
	    	    int minITSclus=1,   int maxITSclus=6,
		    int vertexStatus=2
            	  )
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(10);
  TH1::AddDirectory(kFALSE);

  TFile *file = TFile::Open(inputFile);
  if(!file){
    printf("Error: No file %s in folder.\n", inputFile);
    return;
  }

  TList *list = static_cast<TList*>(file->Get("esd_thnsparse"));
  if(!list){
    printf("Error: No List contained in file %s.\n", inputFile);
    return;
  }

  THnSparseF *heventHLT = (THnSparse*)list->FindObject("fEventHLT"); 
  if(!heventHLT){
      printf("Error: There is no HLT THnSparse object in file %s\n", inputFile);
  }
  THnSparseF *heventOFF = (THnSparse*)list->FindObject("fEventOFF");  
  if(!heventOFF){
      printf("Error: There is no OFF THnSparse object in file %s\n", inputFile);
  } 
  THnSparseF *htrackHLT = (THnSparse*)list->FindObject("fTrackHLT");
  if(!htrackHLT){
      printf("Error: No HLT THnSparse object found\n");
      return;
  } 
  THnSparseF *htrackOFF = (THnSparse*)list->FindObject("fTrackOFF");  
  if(!htrackOFF){
      printf("Error: No OFF THnSparse object found\n");
      return;
  }
    
  
  TText *hText = list->FindObject("text");
  if(!hText) printf("No hText\n");
  
  TString t = "event properties for ";
  t+=hText->GetTitle();
   
  TCanvas *can0 = new TCanvas("can0",t,900,600); 
  can0->Divide(3,2);
  plotEventQuantities(can0,heventOFF,heventHLT,hText);
   
  TCanvas *can1 = new TCanvas("can1","track properties",1100,900); 
  can1->Divide(4,3);
  plotTrackQuantities(can1, htrackOFF, htrackHLT, hText, 
                      minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, minITSclus, maxITSclus, minTrackMult, maxTrackMult, vertexStatus);
  
  TCanvas *can2 = new TCanvas("can2","2-D HLT track correlations",1200,800); 
  can2->Divide(4,2);
  plot2D(can2, htrackHLT, minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, minITSclus, maxITSclus, minTrackMult, maxTrackMult);
  
  TCanvas *can3 = new TCanvas("can3","2-D OFF track correlations",1200,800); 
  can3->Divide(4,2);
  plot2D(can3, htrackOFF, minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, minITSclus, maxITSclus, minTrackMult, maxTrackMult);
  
  TString folder = hText->GetTitle();
  folder.ReplaceAll(" ",""); 
  folder.Remove(0,3);
  folder.Remove(6,15);
   
  gSystem->Exec("mkdir "+folder);
  can0->SaveAs(folder+"/event_properties.root");
  
 
  stringstream sMinEta, sMaxEta;
  sMinEta << minEta; sMaxEta << maxEta;
 
  stringstream sMinTM, sMaxTM;
  sMinTM << minTrackMult; sMaxTM << maxTrackMult;
 
  stringstream sMinPt, sMaxPt;
  sMinPt << minPt; sMaxPt << maxPt;
 
  stringstream sMinDCAr, sMaxDCAr;
  sMinDCAr << minDCAr; sMaxDCAr << maxDCAr;

  stringstream sMinDCAz, sMaxDCAz;
  sMinDCAz << minDCAz; sMaxDCAz << maxDCAz;

  stringstream sMinTPCclus, sMaxTPCclus;
  sMinTPCclus << minTPCclus; sMaxTPCclus << maxTPCclus;
 
  stringstream sMinITSclus, sMaxITSclus;
  sMinITSclus << minITSclus; sMaxITSclus << maxITSclus;
  
  stringstream sVS; sVS << vertexStatus;
  
  TString trackName = "track_properties_";
  trackName += "eta"+sMinEta.str()+"_"+sMaxEta.str()+"_";
  trackName += "Pt"+sMinPt.str()+"_"+sMaxPt.str()+"_";
  trackName += "TM"+sMinTM.str()+"_"+sMaxTM.str()+"_";
  trackName += "DCAr"+sMinDCAr.str()+"_"+sMaxDCAr.str()+"_";
  trackName += "DCAz"+sMinDCAz.str()+"_"+sMaxDCAz.str()+"_";
  trackName += "TPCclus"+sMinTPCclus.str()+"_"+sMaxTPCclus.str()+"_";
  trackName += "ITSclus"+sMinITSclus.str()+"_"+sMaxITSclus.str()+"_";
  if(vertexStatus==2){
     trackName.Chop(); trackName += ".root";
  }
  else { 
     trackName += "VS"+sVS.str();   
     trackName += ".root";
  }
  can1->SaveAs(folder+"/"+trackName);
  
  can2->SaveAs(folder+"/HLT_2D_track_correlations.root");
  can3->SaveAs(folder+"/OFF_2D_track_correlations.root");

  return;
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

void plotAid(TCanvas* can, THnSparse* hOFF, THnSparse* hHLT, TText* hText, TH1D *hlt, TH1D *off, TLegend *l, int size,
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
  
      hlt = hHLT->Projection(i);
      off = hOFF->Projection(i); 
      hlt->SetTitle(hHLT->Projection(i)->GetTitle()); 
      defineYaxisMax(hlt, off);
      off->SetLineColor(2);
     
      can->cd(i+1);
      hlt->Draw();
      off->Draw("sames");
      printStats(hlt, off);
   } 
   printLegend(l,hlt,off);  
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
  h->Projection(1,0)->Draw("colz");	      
  can->cd(2);
  h->Projection(1,3)->Draw("colz");	   
  can->cd(3);
  h->Projection(1,5)->Draw("colz"); 
  can->cd(4);
  h->Projection(1,6)->Draw("colz"); 
  
  can->cd(5);
  h->Projection(6,0)->Draw("colz");	      
  can->cd(6);
  h->Projection(6,3)->Draw("colz");	   
  can->cd(7);
  h->Projection(3,0)->Draw("colz"); 
  can->cd(8);
  h->Projection(3,4)->Draw("colz"); 
 
}

void plotAid(TCanvas* can, THnSparse* hOFF, THnSparse* hHLT, TText* hText, TH1D *hlt, TH1D *off, TLegend *l, int size){
 
  for(int i=0; i<size; i++){         
      hlt = hHLT->Projection(i);
      off = hOFF->Projection(i); 
      hlt->SetTitle(hHLT->Projection(i)->GetTitle()); 
      defineYaxisMax(hlt, off);
      off->SetLineColor(2);
     
      can->cd(i+1);
      hlt->Draw();
      off->Draw("sames");
      printStats(hlt, off);
   } 
   printLegend(l,hlt,off);  
}

void plotEventQuantities(TCanvas* can, THnSparse* heventOFF, THnSparse* heventHLT, TText* hText){

  TH1D *hlt = NULL;
  TH1D *off = NULL;
 
  TLegend *leg1 = new TLegend(0.6,0.6,0.8,0.8);
 
  plotAid(can, heventOFF, heventHLT, hText, hlt, off, leg1, 6);  
  return;
}

void plotTrackQuantities( TCanvas* can, THnSparse* htrackOFF, THnSparse* htrackHLT, TText* hText,
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
  plotAid(can, htrackOFF, htrackHLT, hText, hlt, off, leg1, 11, 
          minEta, maxEta, minPt, maxPt, minDCAr, maxDCAr, minDCAz, maxDCAz, minTPCclus, maxTPCclus, minITSclus, maxITSclus, minTrackMult, maxTrackMult, VS);  
return;
}
