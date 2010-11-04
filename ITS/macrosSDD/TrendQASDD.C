#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TFileMerger.h>
#include <TAlienFile.h>
//#include <TExec.h>
#include <TSystem.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <Riostream.h>
#include <TObjArray.h>
#include <TClass.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#endif
void TrendQASDD( Char_t filelist[150]="LHC10b.txt", Char_t filename[20]="FileQAtrend",Char_t pass[7]="pass1",Int_t year=2010, Char_t period[7]="LHC10b", Bool_t kUseOriginalFile=kFALSE , Char_t especiename[50]="LowMultiplicity")
{
  Int_t  m;
  vector<Int_t> n;
  //Char_t FileName[200];//="/home/msicilia/Desktop/createlists";
  char FileName[150];
  Char_t FileName1[150];
  sprintf(FileName1,"%s",filelist);
  FILE * pFile;
  pFile = fopen (FileName1,"r");
  while (fscanf (pFile,"%d\n",&m)==1){
    
    if(kUseOriginalFile)sprintf(FileName,"run%d/%s/File.QA.%d.%s.%s.Run.%d.root",m,pass,year,period,pass,m);
    else sprintf(FileName, "run%d/%s/%s.root", m,pass,filename);
    if(gSystem->Exec(Form("ls %s >/dev/null 2>&1", FileName))==0){
      TFile mergedfile(FileName);
      if(kUseOriginalFile){
      if(mergedfile.GetKey("ITS")==0x0){
	printf("Run %d, In this run ITS QA has not been executed.-- Exit file\n",m);
	continue;
      }
      else {n.push_back(m);} //this is to grow the vector
      }else {n.push_back(m);}
    }
    else{printf("Run %d,: no file %s present -- Continue\n",m, FileName);continue;}
  }
  fclose (pFile);
  
  Char_t filepath[200];
  TH1F *histocharge;
  TH1F *histochargetrend[2];
  
  char name[100];
  char title[100];
  
  for(Int_t layer=3;layer<5;layer++)
    {
      sprintf(name,"histochargetrend_%d", layer);
      sprintf(title,"Trend of the Charge of the Layer %d", layer);
      histochargetrend[layer-3]=new TH1F(name,title,n.size(),0.,(float)n.size());
    }
  
  gStyle->SetOptStat(0);	
  TCanvas *canvas[6];
  
  for(Int_t i=0;i<3;i++)
    {
      sprintf(name,"canvas%i",i);
      canvas[i]=new TCanvas(name,name);
      canvas[i]->Divide(2,1);
    }
  
  TLegend *legend =new TLegend(0.83,0.8,0.97,0.2);
  TLegend *legend2=new TLegend(0.83,0.8,0.97,0.2);

  Char_t runnmbr[20];
  Float_t fmax=0;
  Float_t fmaxold=0;
  Float_t fmaxmargin=0;

  Float_t fmaxtime=0;
  Float_t fmaxoldtime=0;
  Float_t fmaxmargintime=0;

  Int_t j=0;
  Int_t hh=0;
  //Int_t kk=0;
  char buffer [10];

  printf("%d \n",n.size());
  
  for(Int_t layer=3;layer<5;layer++){
    for (Int_t irun=0; irun<(Int_t)(n.size());irun++)
      {
	if(kUseOriginalFile)sprintf(FileName,"run%d/%s/File.QA.%d.%s.%s.Run.%d.root",n[irun],pass,year,period,pass,n[irun]);
	else sprintf(FileName, "run%d/%s/%s.root", n[irun],pass,filename);
	
	TFile mergedfile(FileName);
	if(kUseOriginalFile){
	if(mergedfile.GetKey("ITS")==0x0){
	  printf("Run %d, In this run ITS QA has not been executed.-- Exit file\n",n[irun]);
	  continue;
	}
	}
	if(kUseOriginalFile)sprintf(filepath,"ITS/RecPoints/%s/%s_SDDLay%dTotCh",especiename,especiename,layer);
	else sprintf(filepath,"%s_SDDLay%dTotCh",especiename,layer);
	histocharge=(TH1F*)(mergedfile.Get(filepath));	 
	
	if(histocharge)
	  {
	    sprintf(buffer, "%d", n[irun]);
	    histochargetrend[layer-3]->Fill(buffer,histocharge->GetMean());
	    histochargetrend[layer-3]->SetMinimum(0.);
	    histochargetrend[layer-3]->SetMaximum(250.);
	    histochargetrend[layer-3]->GetXaxis()->SetTitle("# Run");
	    histochargetrend[layer-3]->GetYaxis()->SetTitle("keV");
	    fmax=histocharge->GetMaximum();
	    if (j==0)	{ fmaxold=fmax;}
	    if (j!=0) { if(fmaxold<fmax) { fmaxold=fmax;} }
	    j++;
	  }//end if histocharge
	
	histocharge=NULL;
	
	//drift time
	
	if(kUseOriginalFile)sprintf(filepath,"ITS/RecPoints/%s/%s_SDDdrifttime_Layer%d",especiename,especiename,layer);
	else sprintf(filepath,"%s_SDDdrifttime_Layer%d",especiename,layer);
	histocharge=(TH1F*)(mergedfile.Get(filepath));	 
	
	if(histocharge)
	  {
	    sprintf(buffer, "%d", n[irun]);	
	    fmaxtime=histocharge->GetMaximum();
	    if (j==0)	{ fmaxoldtime=fmaxtime;}
	    if (j!=0) { if(fmaxoldtime<fmaxtime) { fmaxoldtime=fmaxtime;} }
	    j++;
	  }//end if histocharge
	
	histocharge=NULL;
	
      }//end for run
    fmaxmargin=fmaxold;
    fmaxmargintime=1.023*fmaxoldtime;
    printf("Max Y Range is %4.1f time %4.1f\n",fmaxmargin,fmaxmargintime);
    
    canvas[0]->SetFillColor(46);
    canvas[0]->cd(layer-2)->SetFillColor(46);
    histochargetrend[layer-3]->SetMarkerStyle(20);
    histochargetrend[layer-3]->SetMarkerSize(1);
    histochargetrend[layer-3]->SetMarkerColor(kGreen+4);
    histochargetrend[layer-3]->SetFillColor(42);
    canvas[0]->cd(layer-2)->SetFrameFillColor(kAzure-9);
    histochargetrend[layer-3]->SetBarWidth(0.9);
    histochargetrend[layer-3]->SetBarOffset(0.05);
    histochargetrend[layer-3]->DrawCopy("P");
 
    j=0;
    hh=0;
    Int_t color=0;
    for (Int_t irun=0; irun< (Int_t)(n.size());irun++)
      {
	
	if(kUseOriginalFile) sprintf(FileName, "run%d/%s/File.QA.%d.%s.%s.Run.%d.root", n[irun],pass,year,period,pass,n[irun]);
	else sprintf(FileName, "run%d/%s/%s.root", n[irun],pass,filename);	   
	
	TFile mergedfile(FileName);
	if(kUseOriginalFile){
	if(mergedfile.GetKey("ITS")==0x0){
	  printf("run %d In this run ITS QA has not been executed.--- Exit file\n",n[irun]);
	  continue;
	}
	}	
	if(kUseOriginalFile)sprintf(filepath,"ITS/RecPoints/%s/%s_SDDLay%dTotCh",especiename,especiename,layer);
	else sprintf(filepath,"%s_SDDLay%dTotCh",especiename,layer);	
	sprintf(runnmbr,"run%d",n[irun]);
	printf("%s  %s",filepath,runnmbr);	
	histocharge=(TH1F*)(mergedfile.Get(filepath));	
	if(histocharge)
	  {	    
	    color=2+j;
	    if (color>=12 && color<=27) { color=j+28; }
	    histocharge->SetLineColor(color);
	    if (j==0)
	      {	
		canvas[1]->cd(layer-2);

		canvas[1]->cd(layer-2)->SetFrameFillColor(kGray+3);
		histocharge->SetMaximum(fmaxmargin);
		histocharge->GetYaxis()->SetTitle("");
		histocharge->DrawNormalized();
		histocharge->SetMaximum(fmaxmargin);
		canvas[1]->cd(layer-2)->Update();
	      }	   
	    if (j!=0)
	      {
	
		canvas[1]->cd(layer-2);
		histocharge->DrawNormalized("same");
		canvas[1]->cd(layer-2)->Update();
	      }	    
	    if(layer==3)legend->AddEntry(histocharge,runnmbr,"l");
	    canvas[1]->cd(layer-2);
	    legend->Draw();	    
	    canvas[1]->Update();
	    j++;
	    printf("...Found\n");
	    
	  }//end if histocharge
	else{printf("...Not Found....\n");}
	histocharge=NULL;
	
	    //======================================================== Drift Time	
	if(kUseOriginalFile)sprintf(filepath,"ITS/RecPoints/%s/%s_SDDdrifttime_Layer%d",especiename,especiename,layer);
	else sprintf(filepath,"%s_SDDdrifttime_Layer%d",especiename,layer);	
	sprintf(runnmbr,"run %d",n[irun]);
	printf("%s  %s",filepath,runnmbr);	
	histocharge=(TH1F*)(mergedfile.Get(filepath));	
	if(histocharge)
	  {	    
	    printf("...Found\n");
	    histocharge->SetLineColor(color);
	    if (hh==0)
	      {	
		canvas[2]->cd(layer-2);
		canvas[2]->cd(layer-2)->SetFrameFillColor(kGray+3);
		histocharge->SetMaximum(fmaxmargintime);
		histocharge->GetYaxis()->SetTitle("");
		histocharge->DrawNormalized();
		histocharge->SetMaximum(fmaxmargintime);
	      }	   
	    if (hh!=0)
	      {
		canvas[2]->cd(layer-2);
		histocharge->DrawNormalized("same");
	      }	    
	    if(layer==3)legend2->AddEntry(histocharge,runnmbr,"l");
	    canvas[2]->cd(layer-2);
	    legend2->Draw();	    
	    canvas[2]->Update();
	    hh++;
	  }//end if histocharge
	else{printf("...Not Found....\n");}
	histocharge=NULL;
      }//enf for 
      }//end for layer
  
  TFile trendfile(Form("SDDQAtrend%s%s.root",period,pass),"recreate");  
  trendfile.cd();
  for(Int_t ican=0;ican<3;ican++)canvas[ican]->Write();
  trendfile.Close();
  Char_t psfile[50];
  sprintf(psfile,"SDDtrend%s%s.ps",period,pass);
  canvas[0]->Print(Form("%s[",psfile));
  for(Int_t ifile=0;ifile<3;ifile++){canvas[ifile]->Print(psfile);}

  canvas[2]->Print(Form("%s]",psfile));

  
  delete histocharge;
	



}//end macro
