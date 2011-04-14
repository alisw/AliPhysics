// $Id$
/*
 * Macro to overlay the histograms produced by  
 * HLT/QA/tasks/macros/drawTHnSparse.C
 * 
 * It assumes a file where the input is specified in 
 * the following format:
 * 
 * number of files
 * file1 legend1
 * file2 legend2
 * ...
 * ...
 * 
 * @ingroup alihlt_qa
 * @author Kalliopi.Kanaki@ift.uib.no 
 */
void overlayPlots(const char* option="HLT"/* or "OFF" */, string fi="files.txt"){
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetTitleX(gStyle->GetPadLeftMargin());
 
  char filenames[100];
  sprintf(filenames,"%s",fi.c_str());
  ifstream in(filenames);
  if(!in){
    printf("File %s does not exist", fi.Data());
    break;
  }
  string c;
  TString f;
  int nr_textfile = 0;

  in>>nr_textfile;
  if(!in.good()) break;
  printf("Number of files: %d\n", nr_textfile);

  const int nr_files = nr_textfile;
  TString file[nr_files];
  string cutnames[nr_files];

  nr_textfile=0;
  while(nr_textfile < nr_files){
    in >> f >> c;
    if(!in.good()) break;
    file[nr_textfile] = f;
    cutnames[nr_textfile] = c; 
    Printf("\nfile %d : %s", nr_textfile, f.Data());
    nr_textfile++;

  }
  in.close();
  
  TCanvas *ca;
  TFile   *ff; 
  TPad    *pad; 
  TH1D    *g[nr_files];
  
  TCanvas *d = new TCanvas("d",Form("%s cut studies",option),1100,900);
  d->Divide(3,3);
  //d->Divide(3,2);
  
  TLegend *l = new TLegend(0.6,0.6,0.8,0.8);
  l->SetFillColor(10);
  l->SetLineColor(10);
  
  char cut[100];  
   
  //for(int j=1; j<7; j++){ 
  for(int j=1; j<10; j++){ 
     for(int i=0; i<nr_files; i++){ 
            
        ff = TFile::Open(file[i].Data());   
        if(!ff || ff->IsZombie()){
           printf("Non-existent, corrupted or zombie file %s\n", file[i].Data());
           return;
        } 
        //ca  = (TCanvas*)ff->GetObjectUnchecked("can0");		    
        ca  = (TCanvas*)ff->GetObjectUnchecked("can3");		    
	if(!ca){
	   printf("Empty canvas in file %s.\n", file[i].Data());
	   continue;
	}	
        //pad = (TPad*)ca->GetListOfPrimitives()->FindObject(Form("can0_%d",j));         	
        pad = (TPad*)ca->GetListOfPrimitives()->FindObject(Form("can3_%d",j));         	
        if(!pad){
           printf("Empty pad in canvas %s.\n", ca->GetName());
           continue;	     
        }
        //g[i] =(TH1D*)pad->FindObject(Form("fEvent%s_proj_%d",option,j-1));
        g[i] =(TH1D*)pad->FindObject(Form("fTrack%s_proj_%d",option,j-1));
	if(!g[i]){
	   printf("Empty histogram for i=%d, file %s.\n", i, file[i].Data());
	   continue;
	}
        
        d->cd(j);
        if(i==0){
	  g[i]->SetLineColor(kBlack); 
	  TPaveStats *st = (TPaveStats*)g[i]->FindObject("stats"); 
	  st->SetTextColor(kBlack);
	  g[i]->Draw();
	}
        else { 
	  g[i]->SetLineColor(i+1); 
	  defineYaxisMax(g[0], g[i]);
	  g[i]->Draw("sames");
	}					 
        if(i>0) printStats(g[i-1], g[i]);
        
        ff->Close();
        sprintf( cut,"%s",cutnames[i].c_str() );
	if(j==1) l->AddEntry(g[i],cut,"l");	    	
	else continue;
    }
    if(j==1) l->Draw("same");
  }
  return;
}

void printStats(TH1D *h1, TH1D *h2){  
  gPad->Update();
  TPaveStats *st1 = (TPaveStats*)h1->FindObject("stats");
  st1->SetLineColor(0);

  gPad->Update();
  TPaveStats *st2 = (TPaveStats*)h2->FindObject("stats");
  st2->SetY2NDC(st1->GetY1NDC()-0.05);
  st2->SetY1NDC(st2->GetY2NDC()-TMath::Abs(st1->GetY1NDC()-st1->GetY2NDC()));
  st2->SetLineColor(0);
  st2->SetTextColor(h2->GetLineColor());
  st2->SetFillStyle(0);
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
