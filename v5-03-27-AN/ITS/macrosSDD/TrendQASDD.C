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
#include <TROOT.h>
#endif
void TrendQASDD( Char_t filelist[150]="LHC10b.txt", Char_t filename[20]="FileQAtrend",Char_t pass[7]="pass1",Int_t year=2010, Char_t period[7]="LHC10b", Bool_t kDOSuperimpose=kFALSE,Bool_t kUseOriginalFile=kFALSE , Char_t especiename[50]="LowMultiplicity")
{
  gROOT->SetStyle("Plain");
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

  char name[100];
  char title[100];

  TH1F *histocharge;
  TH1F *histochargerp;

  sprintf(name,"histoRAWeventtrend");
  sprintf(title,"Events used for the QA analysis - RAW");
  TH1F *histoRAWeventtrend=new TH1F(name,title,n.size(),0.,(float)n.size());
  histoRAWeventtrend->GetXaxis()->SetTitle("# Run");
  histoRAWeventtrend->GetXaxis()->SetTicks("-");
  histoRAWeventtrend->GetXaxis()->SetLabelSize(0.02);
  histoRAWeventtrend->GetXaxis()->SetTitleSize(0.02);
  histoRAWeventtrend->GetXaxis()->SetTitleOffset(-2.);
  histoRAWeventtrend->GetXaxis()->SetTickLength(0.01);
  histoRAWeventtrend->GetYaxis()->SetTitle("Events");
  histoRAWeventtrend->GetYaxis()->SetLabelSize(0.02);
  histoRAWeventtrend->GetYaxis()->SetTitleSize(0.02);
  histoRAWeventtrend->GetYaxis()->SetTitleOffset(1.5);
  histoRAWeventtrend->SetMarkerStyle(20);
  histoRAWeventtrend->SetMarkerColor(1);

  sprintf(name,"histoRPeventtrend");
  sprintf(title,"Events used for the QA analysis - RP");
  TH1F *histoRPeventtrend=new TH1F(name,title,n.size(),0.,(float)n.size());
  histoRPeventtrend->GetXaxis()->SetTitle("# Run");
  histoRPeventtrend->GetXaxis()->SetTicks("-");
  histoRPeventtrend->GetXaxis()->SetLabelSize(0.02);
  histoRPeventtrend->GetXaxis()->SetTitleSize(0.02);
  histoRPeventtrend->GetXaxis()->SetTitleOffset(-2.);
  histoRPeventtrend->GetXaxis()->SetTickLength(0.01);
  histoRPeventtrend->GetYaxis()->SetTitle("Events");
  histoRPeventtrend->GetYaxis()->SetLabelSize(0.02);
  histoRPeventtrend->GetYaxis()->SetTitleSize(0.02);
  histoRPeventtrend->GetYaxis()->SetTitleOffset(1.5);
  histoRPeventtrend->SetMarkerStyle(24);
  histoRPeventtrend->SetMarkerColor(4);


  TH1F *histoRAWnormevents[3];
  TH1F *histoRPnormevents[3];
  
  TH1F *histoRAWfilledmodules[3];
  TH1F *histoRPfilledmodules[3];

  TH1F *histoRAWfilleddr[3];
  TH1F *histoRPfilleddr[3];

  TH1F *histoRAWactivemodules[3];
  TH1F *histoRPactivemodules[3];

  TH1F *histoRAWactivedr[3];
  TH1F *histoRPactivedr[3];

  TH1F *histoRAWoverthmodules[3];
  TH1F *histoRPoverthmodules[3];

  TH1F *histoRAWoverthdr[3];
  TH1F *histoRPoverthdr[3];

  TH1F *histochargetrend[2];
  TH1F *historadiustrend[2];
  

  //index on the histograms

  Int_t event=1;
  Int_t normevtRAW[3]={2,3,4};
  Int_t normevtRP[3]={2,3,4};

  Int_t filledmodRAW[3]={6,20,34};
  Int_t filledmodRP[3] ={6,20,36};

  Int_t filleddrRAW[3]={8,22,36};
  Int_t filleddrRP[3] ={8,22,38};

  Int_t activemodRAW[3]={5,19,33};
  Int_t activemodRP[3] ={5,19,35};

  Int_t activedrRAW[3]={7,21,35};
  Int_t activedrRP[3] ={7,21,37};

  Int_t overthmodRAW[3]={17,31,45};
  Int_t overthmodRP[3] ={17,33,49};

  Int_t overthdrRAW[3]={18,32,46};
  Int_t overthdrRP[3] ={18,34,50};

  Int_t chargeindex[2]={31,47};
  Int_t radiusindex[2]={32,48};

  const TString histoname[3]={"Total","Layer 3","Layer 4"};
  const TString histolegend[3]={"All","Layer 3","Layer 4"};
  const TString histolegend2[2]={"Filled","Active"};
  const TString histolegend3[2]={"Raw","RecPoints"};

  const TString drawoptionevents[3]={"P","PSAME","PSAME"};
  const TString drawoptionfill[3]={"P","P","PSAME"};
  const TString drawoptionactive[3]={"PSAME","PSAME","PSAME"};

  TLegend *legend =new TLegend(0.83,0.8,0.97,0.2);
  TLegend *legend2=new TLegend(0.83,0.8,0.97,0.2);
  TLegend *legendevents=new TLegend(0.83,01.,0.97,0.9);
  TLegend *legendchecks =new TLegend(0.83,1.,0.97,0.9);
  TLegend *legendchecks2 =new TLegend(0.83,1.,0.97,0.9);
  TLegend *legendchecks3 =new TLegend(0.83,1.,0.97,0.9);
  TLegend *legendchecks4 =new TLegend(0.83,1.,0.97,0.9);

  Char_t runnmbr[20];
  Float_t fmax[2];
  for(Int_t i=0;i<2;i++)fmax[i]=0;
  Float_t fmaxold=0;
  Float_t fmaxmargin[2];
  for(Int_t i=0;i<2;i++)fmaxmargin[i]=0;

  Float_t fmaxtime=0;
  Float_t fmaxoldtime=0;
  Float_t fmaxmargintime[2];
  for(Int_t i=0;i<2;i++)fmaxmargintime[i]=0;

  Int_t j=0;
  Int_t hh=0;

  for(Int_t ii=0;ii<3;ii++)

    {
      sprintf(name,"histoRAWnormevent%d", ii);
      sprintf(title,"Trend of the Norm Events  %s RAW",histoname[ii].Data());
      histoRAWnormevents[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRAWnormevents[ii]->SetStats(0);
      histoRAWnormevents[ii]->SetMinimum(0.);
      histoRAWnormevents[ii]->GetXaxis()->SetTitle("# Run");
      histoRAWnormevents[ii]->GetXaxis()->SetTicks("-");
      histoRAWnormevents[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRAWnormevents[ii]->GetXaxis()->SetTitleSize(0.02);
      histoRAWnormevents[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRAWnormevents[ii]->GetXaxis()->SetTickLength(0.01);
      histoRAWnormevents[ii]->GetYaxis()->SetTitle("entries/(goodch*totevents)");
      histoRAWnormevents[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRAWnormevents[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRAWnormevents[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRAWnormevents[ii]->SetMarkerStyle(20+ii);
      histoRAWnormevents[ii]->SetMarkerColor(1+ii);

      sprintf(name,"histoRPnormevent%d", ii);
      sprintf(title,"Trend of the Norm Events  %s RP", histoname[ii].Data());
      histoRPnormevents[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRPnormevents[ii]->SetStats(0);
      histoRPnormevents[ii]->SetMinimum(0.);
      histoRPnormevents[ii]->GetXaxis()->SetTitle("# Run");
      histoRPnormevents[ii]->GetXaxis()->SetTicks("-");
      histoRPnormevents[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRPnormevents[ii]->GetXaxis()->SetTickLength(0.01);
      histoRPnormevents[ii]->GetXaxis()->SetTitleSize(0.02);
      histoRPnormevents[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRPnormevents[ii]->GetYaxis()->SetTitle("entries/(goodch*totevents)");
      histoRPnormevents[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRPnormevents[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRPnormevents[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRPnormevents[ii]->SetMarkerStyle(20+ii);
      histoRPnormevents[ii]->SetMarkerColor(1+ii);


      sprintf(name,"histoRAWfilledmodules%d", ii);
      sprintf(title,"Trend of the filled modules %s RAW", histoname[ii].Data());
      histoRAWfilledmodules[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRAWfilledmodules[ii]->SetStats(0);
      histoRAWfilledmodules[ii]->SetMaximum(272.);
      histoRAWfilledmodules[ii]->SetMinimum(-1.);
      histoRAWfilledmodules[ii]->GetXaxis()->SetTitle("# Run");
      histoRAWfilledmodules[ii]->GetXaxis()->SetTicks("-");
      histoRAWfilledmodules[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRAWfilledmodules[ii]->GetXaxis()->SetTickLength(0.01);
      histoRAWfilledmodules[ii]->GetXaxis()->SetTitleSize(0.02);
      histoRAWfilledmodules[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRAWfilledmodules[ii]->GetYaxis()->SetTitle("modules");
      histoRAWfilledmodules[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRAWfilledmodules[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRAWfilledmodules[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRAWfilledmodules[ii]->SetMarkerStyle(20+ii);
      histoRAWfilledmodules[ii]->SetMarkerColor(1+ii);

      sprintf(name,"histoRPfilledmodules%d", ii);
      sprintf(title,"Trend of the filled modules %s RP",  histoname[ii].Data());
      histoRPfilledmodules[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRPfilledmodules[ii]->GetXaxis()->SetTitle("# Run");
      histoRPfilledmodules[ii]->SetStats(0);
      histoRPfilledmodules[ii]->SetMaximum(272.);
      histoRPfilledmodules[ii]->SetMinimum(-1.);
      histoRPfilledmodules[ii]->GetXaxis()->SetTicks("-");
      histoRPfilledmodules[ii]->GetYaxis()->SetTitle("modules");
      histoRPfilledmodules[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRPfilledmodules[ii]->GetXaxis()->SetTickLength(0.01);
      histoRPfilledmodules[ii]->GetXaxis()->SetTitleSize(0.02);
      histoRPfilledmodules[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRPfilledmodules[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRPfilledmodules[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRPfilledmodules[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRPfilledmodules[ii]->SetMarkerStyle(20+ii);
      histoRPfilledmodules[ii]->SetMarkerColor(1+ii);

      sprintf(name,"histoRAWfilleddr%d", ii);
      sprintf(title,"Trend of the filled drift regions %s RAW",  histoname[ii].Data());
      histoRAWfilleddr[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRAWfilleddr[ii]->GetXaxis()->SetTitle("# Run");
      histoRAWfilleddr[ii]->SetStats(0);
      histoRAWfilleddr[ii]->SetMaximum(272.);
      histoRAWfilleddr[ii]->SetMinimum(-1.);
      histoRAWfilleddr[ii]->GetXaxis()->SetTicks("-");
      histoRAWfilleddr[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRAWfilleddr[ii]->GetXaxis()->SetTickLength(0.01);
      histoRAWfilleddr[ii]->GetXaxis()->SetTitleSize(0.02);
      histoRAWfilleddr[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRAWfilleddr[ii]->GetYaxis()->SetTitle("Drift Regions");
      histoRAWfilleddr[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRAWfilleddr[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRAWfilleddr[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRAWfilleddr[ii]->SetMarkerStyle(20+ii);
      histoRAWfilleddr[ii]->SetMarkerColor(1+ii);

      sprintf(name,"histoRPfilleddr%d", ii);
      sprintf(title,"Trend of the filled drift regions %s RP", histoname[ii].Data());
      histoRPfilleddr[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRPfilleddr[ii]->SetStats(0);
      histoRPfilleddr[ii]->SetMaximum(272.);
      histoRPfilleddr[ii]->SetMinimum(-1.);
      histoRPfilleddr[ii]->GetXaxis()->SetTitle("# Run");
      histoRPfilleddr[ii]->GetXaxis()->SetTicks("-");
      histoRPfilleddr[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRPfilleddr[ii]->GetXaxis()->SetTickLength(0.01);
      histoRPfilleddr[ii]->GetXaxis()->SetTitleSize(0.02);
      histoRPfilleddr[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRPfilleddr[ii]->GetYaxis()->SetTitle("Drift Regions");
      histoRPfilleddr[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRPfilleddr[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRPfilleddr[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRPfilleddr[ii]->SetMarkerStyle(20+ii);
      histoRPfilleddr[ii]->SetMarkerColor(1+ii);


      sprintf(name,"histoRAWactivemodules%d", ii);
      sprintf(title,"Trend of the active module %s RAW",  histoname[ii].Data());
      histoRAWactivemodules[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRAWactivemodules[ii]->SetStats(0);
      histoRAWactivemodules[ii]->GetXaxis()->SetTitle("# Run");
      histoRAWactivemodules[ii]->GetXaxis()->SetTicks("-");
      histoRAWactivemodules[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRAWactivemodules[ii]->GetXaxis()->SetTickLength(0.01);
      histoRAWactivemodules[ii]->GetXaxis()->SetTitleSize(0.02);
      histoRAWactivemodules[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRAWactivemodules[ii]->GetYaxis()->SetTitle("modules");
      histoRAWactivemodules[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRAWactivemodules[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRAWactivemodules[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRAWactivemodules[ii]->SetMarkerStyle(24+ii);
      histoRAWactivemodules[ii]->SetMarkerColor(4);



      sprintf(name,"histoRPactivemodules%d", ii);
      sprintf(title,"Trend of the active modules %s RP",  histoname[ii].Data());
      histoRPactivemodules[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRPactivemodules[ii]->SetStats(0);
      histoRPactivemodules[ii]->GetXaxis()->SetTitle("# Run");
      histoRPactivemodules[ii]->GetXaxis()->SetTicks("-");
      histoRPactivemodules[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRPactivemodules[ii]->GetXaxis()->SetTickLength(0.01);
      histoRPactivemodules[ii]->GetXaxis()->SetTitleSize(0.02);
      histoRPactivemodules[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRPactivemodules[ii]->GetYaxis()->SetTitle("modules");
      histoRPactivemodules[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRPactivemodules[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRPactivemodules[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRPactivemodules[ii]->SetMarkerStyle(24+ii);
      histoRPactivemodules[ii]->SetMarkerColor(4);


      sprintf(name,"histoRAWactivedr%d", ii);
      sprintf(title,"Trend of the active drift region %s RAW", histoname[ii].Data());
      histoRAWactivedr[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRAWactivedr[ii]->SetStats(0);
      histoRAWactivedr[ii]->GetXaxis()->SetTitle("# Run");
      histoRAWactivedr[ii]->GetXaxis()->SetTicks("-");
      histoRAWactivedr[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRAWactivedr[ii]->GetXaxis()->SetTickLength(0.01);
      histoRAWactivedr[ii]->GetXaxis()->SetTitleSize(0.02);
      histoRAWactivedr[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRAWactivedr[ii]->GetYaxis()->SetTitle("Drift Regions");
      histoRAWactivedr[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRAWactivedr[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRAWactivedr[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRAWactivedr[ii]->SetMarkerStyle(24+ii);
      histoRAWactivedr[ii]->SetMarkerColor(4);

      sprintf(name,"histoRPactivedr%d", ii);
      sprintf(title,"Trend of the active drift region %s RP",  histoname[ii].Data());
      histoRPactivedr[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRPactivedr[ii]->SetStats(0);
      histoRPactivedr[ii]->GetXaxis()->SetTitle("# Run");
      histoRPactivedr[ii]->GetXaxis()->SetTicks("-");
      histoRPactivedr[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRPactivedr[ii]->GetXaxis()->SetTickLength(0.01);
      histoRPactivedr[ii]->GetXaxis()->SetTitleSize(0.02);
      histoRPactivedr[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRPactivedr[ii]->GetYaxis()->SetTitle("Drift Regions");
      histoRPactivedr[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRPactivedr[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRPactivedr[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRPactivedr[ii]->SetMarkerStyle(24+ii);
      histoRPactivedr[ii]->SetMarkerColor(4);

      sprintf(name,"histoRAWoverthmodules%d", ii);
      sprintf(title,"Trend of the over th modules %s RAW", histoname[ii].Data());
      histoRAWoverthmodules[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRAWoverthmodules[ii]->SetStats(0);
      histoRAWoverthmodules[ii]->GetXaxis()->SetTitle("# Run");
      histoRAWoverthmodules[ii]->GetXaxis()->SetTicks("-");
      histoRAWoverthmodules[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRAWoverthmodules[ii]->GetXaxis()->SetTickLength(0.01);
      histoRAWoverthmodules[ii]->GetXaxis()->SetTitleSize(0.02);
      histoRAWoverthmodules[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRAWoverthmodules[ii]->GetYaxis()->SetTitle("modules");
      histoRAWoverthmodules[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRAWoverthmodules[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRAWoverthmodules[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRAWoverthmodules[ii]->SetMarkerStyle(20+ii);
      histoRAWoverthmodules[ii]->SetMarkerColor(1+ii);

      sprintf(name,"histoRPoverthmodules%d", ii);
      sprintf(title,"Trend of the over th modules %s RP",  histoname[ii].Data());
      histoRPoverthmodules[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRPoverthmodules[ii]->SetStats(0);
      histoRPoverthmodules[ii]->GetXaxis()->SetTitle("# Run");
      histoRPoverthmodules[ii]->GetXaxis()->SetTicks("-");
      histoRPoverthmodules[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRPoverthmodules[ii]->GetXaxis()->SetTickLength(0.01);
      histoRPoverthmodules[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRPoverthmodules[ii]->GetYaxis()->SetTitle("modules");
      histoRPoverthmodules[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRPoverthmodules[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRPoverthmodules[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRPoverthmodules[ii]->SetMarkerStyle(20+ii);
      histoRPoverthmodules[ii]->SetMarkerColor(1+ii);


      sprintf(name,"histoRAWoverthdr%d", ii);
      sprintf(title,"Trend of the over th drift regions %s RAW",  histoname[ii].Data());
      histoRAWoverthdr[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRAWoverthdr[ii]->SetStats(0);
      histoRAWoverthdr[ii]->GetXaxis()->SetTitle("# Run");
      histoRAWoverthdr[ii]->GetXaxis()->SetTicks("-");
      histoRAWoverthdr[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRAWoverthdr[ii]->GetXaxis()->SetTickLength(0.01);
      histoRAWoverthdr[ii]->GetXaxis()->SetTitleSize(0.02);
      histoRAWoverthdr[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRAWoverthdr[ii]->GetYaxis()->SetTitle("Drift Regions");
      histoRAWoverthdr[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRAWoverthdr[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRAWoverthdr[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRAWoverthdr[ii]->SetMarkerStyle(20+ii);
      histoRAWoverthdr[ii]->SetMarkerColor(1+ii);


      sprintf(name,"histoRPoverthdr%d", ii);
      sprintf(title,"Trend of the over th drift regions %s RP",  histoname[ii].Data());
      histoRPoverthdr[ii]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histoRPoverthdr[ii]->SetStats(0);
      histoRPoverthdr[ii]->GetXaxis()->SetTitle("# Run");
      histoRPoverthdr[ii]->GetXaxis()->SetTicks("-");
      histoRPoverthdr[ii]->GetXaxis()->SetLabelSize(0.02);
      histoRPoverthdr[ii]->GetXaxis()->SetTickLength(0.01);
      histoRPoverthdr[ii]->GetXaxis()->SetTitleSize(0.02);
      histoRPoverthdr[ii]->GetXaxis()->SetTitleOffset(-2.);
      histoRPoverthdr[ii]->GetYaxis()->SetTitle("Drift Regions");
      histoRPoverthdr[ii]->GetYaxis()->SetLabelSize(0.02);
      histoRPoverthdr[ii]->GetYaxis()->SetTitleSize(0.02);
      histoRPoverthdr[ii]->GetYaxis()->SetTitleOffset(1.5);
      histoRPoverthdr[ii]->SetMarkerStyle(20+ii);
      histoRPoverthdr[ii]->SetMarkerColor(1+ii);

    }

  for(Int_t layer=3;layer<5;layer++)
    {
      sprintf(name,"histochargetrend_%d", layer);
      sprintf(title,"Trend of the Charge of the Layer %d", layer);
      histochargetrend[layer-3]=new TH1F(name,title,n.size(),0.,(float)n.size());
      histochargetrend[layer-3]->SetStats(0);
      histochargetrend[layer-3]->GetXaxis()->SetTicks("-");
      histochargetrend[layer-3]->GetXaxis()->SetLabelSize(0.02);
      histochargetrend[layer-3]->GetXaxis()->SetTickLength(0.01);
      histochargetrend[layer-3]->GetYaxis()->SetTitleSize(0.02);
      histochargetrend[layer-3]->GetXaxis()->SetTitleOffset(-2.);
      histochargetrend[layer-3]->GetXaxis()->SetLabelSize(0.02);
      histochargetrend[layer-3]->GetYaxis()->SetTitleSize(0.02);
      histochargetrend[layer-3]->GetYaxis()->SetTitleOffset(1.5);
      histochargetrend[layer-3]->GetXaxis()->SetTitle("# Run");
      histochargetrend[layer-3]->GetYaxis()->SetTitle("keV");
      histochargetrend[layer-3]->SetMinimum(0.);
      histochargetrend[layer-3]->SetMaximum(250.);
      histochargetrend[layer-3]->SetMarkerStyle(20);
      histochargetrend[layer-3]->SetMarkerSize(1);
      histochargetrend[layer-3]->SetMarkerColor(kGreen+4);
      histochargetrend[layer-3]->SetFillColor(42);


      sprintf(name,"historadiustrend_%d", layer);
      sprintf(title,"Trend of the radius of the Layer %d", layer);
      historadiustrend[layer-3]=new TH1F(name,title,n.size(),0.,(float)n.size());
      historadiustrend[layer-3]->SetStats(0);
      historadiustrend[layer-3]->GetXaxis()->SetTicks("-");
      historadiustrend[layer-3]->GetXaxis()->SetLabelSize(0.02);
      historadiustrend[layer-3]->GetXaxis()->SetTickLength(0.01);
      historadiustrend[layer-3]->GetXaxis()->SetTitleSize(0.02);
      historadiustrend[layer-3]->GetXaxis()->SetTitleOffset(-2.);
      historadiustrend[layer-3]->GetYaxis()->SetLabelSize(0.02);
      historadiustrend[layer-3]->GetYaxis()->SetTitleSize(0.02);
      historadiustrend[layer-3]->GetYaxis()->SetTitleOffset(1.5);
      historadiustrend[layer-3]->GetXaxis()->SetTitle("# Run");
      historadiustrend[layer-3]->GetYaxis()->SetTitle("r [cm]");
      historadiustrend[layer-3]->SetMinimum(0.);
      historadiustrend[layer-3]->SetMaximum(26.);
      historadiustrend[layer-3]->SetMarkerStyle(20);
      historadiustrend[layer-3]->SetMarkerSize(1);
      historadiustrend[layer-3]->SetMarkerColor(kGreen+4);
      historadiustrend[layer-3]->SetFillColor(42);
    }
  

  gStyle->SetOptStat(0);
  	
  TCanvas *canvas[12];
  
  for(Int_t i=0;i<12;i++)
    {
      sprintf(name,"canvas%i",i);
      //printf("--> %s \n",name);
      
      canvas[i]=new TCanvas(name,name);

    }
  
  char buffer [10];
  for (Int_t irun=0; irun<(Int_t)(n.size());irun++)
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
	sprintf(buffer, "%d", n[irun]);
	if(kUseOriginalFile)sprintf(filepath,"ITS/Raws/%s/Expert/%s_SDDRawDataCheck",especiename,especiename);
	else sprintf(filepath,"%s_SDDRawDataCheck",especiename);	
	sprintf(runnmbr,"run%d",n[irun]);
	//printf("%s  %s\n",filepath,runnmbr);	
	histocharge=(TH1F*)(mergedfile.Get(filepath));
	if(histocharge)
	  {
	    histoRAWeventtrend->Fill(buffer,(histocharge->GetBinContent(event)));
	  }
	
	if(kUseOriginalFile)sprintf(filepath,"ITS/RecPoints/%s/Expert/%s_SDDRecPointDataCheck",especiename,especiename);
	else sprintf(filepath,"%s_SDDRecPointCheck",especiename);	
	sprintf(runnmbr,"run%d",n[irun]);
	//printf("%s  %s\n",filepath,runnmbr);	
	histochargerp=(TH1F*)(mergedfile.Get(filepath));	
	if(histochargerp)
	  {
	    histoRPeventtrend->Fill(buffer,histochargerp->GetBinContent(event));
	  }
	for(Int_t iindex=0;iindex<3;iindex++)
	  {
	    if(histocharge)
	      {
		histoRAWnormevents[iindex]->Fill(buffer,histocharge->GetBinContent(normevtRAW[iindex]));  
		histoRAWfilledmodules[iindex]->Fill(buffer,histocharge->GetBinContent(filledmodRAW[iindex]));	    
		histoRAWfilleddr[iindex]->Fill(buffer,histocharge->GetBinContent(filleddrRAW[iindex]));
		histoRAWactivemodules[iindex]->Fill(buffer,histocharge->GetBinContent(activemodRAW[iindex]));	    
		histoRAWactivedr[iindex]->Fill(buffer,histocharge->GetBinContent(activedrRAW[iindex]));	    
		histoRAWoverthmodules[iindex]->Fill(buffer,histocharge->GetBinContent(overthmodRAW[iindex]));	  
		histoRAWoverthdr[iindex]->Fill(buffer,histocharge->GetBinContent(overthdrRAW[iindex]));
	      }
	    
	    if(histochargerp)
	      {
		histoRPnormevents[iindex]->Fill(buffer,histochargerp->GetBinContent(normevtRP[iindex]));
		histoRPfilledmodules[iindex]->Fill(buffer,histochargerp->GetBinContent(filledmodRP[iindex]));
		histoRPfilleddr[iindex]->Fill(buffer,histochargerp->GetBinContent(filleddrRP[iindex]));
		histoRPactivemodules[iindex]->Fill(buffer,histochargerp->GetBinContent(activemodRP[iindex]));
		histoRPactivedr[iindex]->Fill(buffer,histochargerp->GetBinContent(activedrRP[iindex]));
		histoRPoverthmodules[iindex]->Fill(buffer,histochargerp->GetBinContent(overthmodRP[iindex]));
		histoRPoverthdr[iindex]->Fill(buffer,histochargerp->GetBinContent(overthdrRP[iindex]));
		if(iindex<2){
		  histochargetrend[iindex]->Fill(buffer,histochargerp->GetBinContent(chargeindex[iindex]));
		  histochargetrend[iindex]->SetBinError(irun+1,histochargerp->GetBinError(chargeindex[iindex]));
		  historadiustrend[iindex]->Fill(buffer,histochargerp->GetBinContent(radiusindex[iindex]));
		  historadiustrend[iindex]->SetBinError(irun+1,histochargerp->GetBinError(radiusindex[iindex]));
		}

	      }
	  }
	histocharge=NULL;
	histochargerp=NULL;
    }//end loop on run
  
  canvas[3]->cd();
  //  histoRAWeventtrend->GetXaxis()->LabelsOption("v");
  histoRAWeventtrend->Draw("E0");
  //  histoRPeventtrend->GetXaxis()->LabelsOption("v");
  histoRPeventtrend->Draw("E0SAME");
  legendevents->AddEntry(histoRAWeventtrend,histolegend3[0].Data(),"P");
  legendevents->AddEntry(histoRPeventtrend,histolegend3[1].Data(),"P");
  legendevents->Draw();
  canvas[3]->Update();
  
  canvas[4]->Divide(2,1);  
  canvas[4]->cd(1);
  Char_t legendtext[50];
  for(Int_t i=0;i<3;i++)
    {
      //      histoRAWnormevents[i]->GetXaxis()->LabelsOption("v");
      histoRAWnormevents[i]->Draw(drawoptionevents[i].Data());
      //      sprintf(legendtext,"%s %s",histolegend);
      legendchecks->AddEntry(histoRAWnormevents[i],histolegend[i].Data(),"P");
    }
  legendchecks->Draw();
  canvas[4]->cd(2);

  for(Int_t i=0;i<3;i++)
    {
      histoRPnormevents[i]->Draw(drawoptionevents[i].Data());
      //legendchecks->AddEntry(histoRPnormevent,histolegend[i].Data(),"P");

    }
  legendchecks->Draw("same");
  canvas[4]->Update();
  
  canvas[5]->Divide(2,1);  
  canvas[5]->cd(1);
  for(Int_t i=1;i<3;i++)
    {
      //histoRAWfilledmodules[i]->GetXaxis()->LabelsOption("v");
      histoRAWfilledmodules[i]->Draw(drawoptionfill[i].Data());
      sprintf(legendtext,"%s %s",histolegend[i].Data(),histolegend2[0].Data());
      legendchecks2->AddEntry(histoRAWfilledmodules[i],legendtext,"P");
      //histoRAWactivemodules[i]->GetXaxis()->LabelsOption("v");
      histoRAWactivemodules[i]->Draw(drawoptionactive[i].Data());
      sprintf(legendtext,"%s %s",histolegend[i].Data(),histolegend2[1].Data());
      legendchecks2->AddEntry(histoRAWactivemodules[i],legendtext,"P");
    }
  legendchecks2->Draw("same");
  canvas[5]->cd(2);
  for(Int_t i=1;i<3;i++)
    {
      //histoRPfilledmodules[i]->GetXaxis()->LabelsOption("v");
      histoRPfilledmodules[i]->Draw(drawoptionfill[i].Data());
      //histoRPactivemodules[i]->GetXaxis()->LabelsOption("v");
      histoRPactivemodules[i]->Draw(drawoptionactive[i].Data());
    }
  legendchecks2->Draw("same");
  canvas[5]->Update();

  canvas[6]->Divide(2,1);  
  canvas[6]->cd(1);
  for(Int_t i=1;i<3;i++)
    {
      //histoRAWfilleddr[i]->GetXaxis()->LabelsOption("v");
      histoRAWfilleddr[i]->Draw(drawoptionfill[i].Data());
      //histoRAWactivedr[i]->GetXaxis()->LabelsOption("v");
      histoRAWactivedr[i]->Draw(drawoptionactive[i].Data());
    }
  legendchecks2->Draw("same");
  canvas[6]->cd(2);
  for(Int_t i=1;i<3;i++)
    {
      //histoRPfilleddr[i]->GetXaxis()->LabelsOption("v");
      histoRPfilleddr[i]->Draw(drawoptionfill[i].Data());
      //histoRPactivedr[i]->GetXaxis()->LabelsOption("v");
      histoRPactivedr[i]->Draw(drawoptionactive[i].Data());
    }
  legendchecks2->Draw("same");
  canvas[6]->Update();



  ///////////////////////////////////////////////


  
  canvas[10]->Divide(2,1);  
  canvas[10]->cd(1);
  
  //histoRAWfilledmodules[i]->GetXaxis()->LabelsOption("v");
  histoRAWfilledmodules[0]->Draw(drawoptionfill[0].Data());
  sprintf(legendtext,"%s %s",histolegend[0].Data(),histolegend2[0].Data());
  legendchecks3->AddEntry(histoRAWfilledmodules[0],legendtext,"P");
  //histoRAWactivemodules[i]->GetXaxis()->LabelsOption("v");
  histoRAWactivemodules[0]->Draw(drawoptionactive[0].Data());
  sprintf(legendtext,"%s %s",histolegend[0].Data(),histolegend2[1].Data());
  legendchecks3->AddEntry(histoRAWactivemodules[0],legendtext,"P");
  legendchecks3->Draw("same");
  canvas[10]->cd(2);
  
  //histoRPfilledmodules[i]->GetXaxis()->LabelsOption("v");
  histoRPfilledmodules[0]->Draw(drawoptionfill[0].Data());
  //histoRPactivemodules[i]->GetXaxis()->LabelsOption("v");
  histoRPactivemodules[0]->Draw(drawoptionactive[0].Data());
  legendchecks3->Draw("same");
  canvas[10]->Update();
  
  canvas[11]->Divide(2,1);  
  canvas[11]->cd(1);
  
  //histoRAWfilleddr[i]->GetXaxis()->LabelsOption("v");
  histoRAWfilleddr[0]->Draw(drawoptionfill[0].Data());
  //histoRAWactivedr[i]->GetXaxis()->LabelsOption("v");
  histoRAWactivedr[0]->Draw(drawoptionactive[0].Data());
  legendchecks3->Draw("same");
  canvas[11]->cd(2);
  
  //histoRPfilleddr[i]->GetXaxis()->LabelsOption("v");
  histoRPfilleddr[0]->Draw(drawoptionfill[0].Data());
  //histoRPactivedr[i]->GetXaxis()->LabelsOption("v");
  histoRPactivedr[0]->Draw(drawoptionactive[0].Data());
  legendchecks3->Draw("same");
  canvas[11]->Update();


  ///////////////////////////////////////////////

  canvas[8]->Divide(2,1);  
  canvas[8]->cd(1);
  for(Int_t i=0;i<3;i++)
    {
      //histoRAWoverthmodules[i]->GetXaxis()->LabelsOption("v");
      histoRAWoverthmodules[i]->Draw(drawoptionfill[i].Data());
      //histoRAWoverthmodules[i]->Draw(drawoptionactive[i].Data());
    }
  canvas[8]->cd(2);
  for(Int_t i=0;i<3;i++)
    {
      //histoRPoverthmodules[i]->GetXaxis()->LabelsOption("v");
      histoRPoverthmodules[i]->Draw(drawoptionfill[i].Data());
      //histoRPoverthmodules[i]->Draw(drawoptionactive[i].Data());
    }
  canvas[8]->Update();

  canvas[9]->Divide(2,1);  
  canvas[9]->cd(1);
  for(Int_t i=0;i<3;i++)
    {
      histoRAWoverthdr[i]->GetXaxis()->LabelsOption("v");
      histoRAWoverthdr[i]->Draw(drawoptionfill[i].Data());
      //histoRAWoverthdr[i]->Draw(drawoptionactive[i].Data());
    }
  canvas[9]->cd(2);
  for(Int_t i=0;i<3;i++)
    {
      //histoRPoverthdr[i]->GetXaxis()->LabelsOption("v");
      histoRPoverthdr[i]->Draw(drawoptionfill[i].Data());
      //histoRPoverthdr[i]->Draw(drawoptionactive[i].Data());
    }
  canvas[9]->Update();
  
  //Int_t kk=0;


  printf("Total runs: %d \n",(int)(n.size()));
  canvas[0]->Divide(2,1);

  for(Int_t layer=3;layer<5;layer++){
    canvas[0]->SetFillColor(46);
    canvas[0]->cd(layer-2)->SetFillColor(46);
    //histochargetrend[layer-3]->GetXaxis()->LabelsOption("v");
    histochargetrend[layer-3]->SetMarkerStyle(20);
    histochargetrend[layer-3]->SetMarkerSize(1);
    histochargetrend[layer-3]->SetMarkerColor(kGreen+4);
    histochargetrend[layer-3]->SetFillColor(42);
    canvas[0]->cd(layer-2)->SetFrameFillColor(kAzure-9);
    //    histochargetrend[layer-3]->GetXaxis()->LabelsOption("v");
    histochargetrend[layer-3]->SetBarWidth(0.9);
    histochargetrend[layer-3]->SetBarOffset(0.05);
    histochargetrend[layer-3]->DrawCopy("P");
    canvas[0]->Update();
 
    canvas[7]->cd();
    canvas[7]->SetFillColor(46);
    //canvas[7]->SetFillColor(46);
    //historadiustrend[layer-3]->GetXaxis()->LabelsOption("v");
    historadiustrend[layer-3]->SetMarkerStyle(20);
    historadiustrend[layer-3]->SetMarkerSize(1);
    historadiustrend[layer-3]->SetMarkerColor(8+layer-3);
    legendchecks4->AddEntry(historadiustrend[layer-3],histolegend[layer-2].Data(),"P");
    //historadiustrend[layer-3]->SetFillColor(42);
    canvas[7]->SetFrameFillColor(kAzure-9);
    historadiustrend[layer-3]->SetBarWidth(0.9);
    historadiustrend[layer-3]->SetBarOffset(0.05);
    historadiustrend[layer-3]->DrawCopy(drawoptionevents[layer-3].Data());
    canvas[7]->Update();
  }
  legendchecks4->Draw("same");
  canvas[7]->Update();
  canvas[1]->Divide(2,1);
  canvas[2]->Divide(2,1);
  //canvas[7]->Divide(2,1);
  
  if(kDOSuperimpose==kTRUE)
    {
      Int_t jj=0;
      for (Int_t irun=0; irun<(Int_t)(n.size());irun++)
	{
	  for(Int_t ilayer=3;ilayer<5;ilayer++)
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
	    if(kUseOriginalFile)sprintf(filepath,"ITS/RecPoints/%s/%s_SDDLay%dTotCh",especiename,especiename,ilayer);
	    else sprintf(filepath,"%s_SDDLay%dTotCh",especiename,ilayer);
	    histocharge=(TH1F*)(mergedfile.Get(filepath));	 
	    
	    if(histocharge)
	      {
		sprintf(buffer, "%d", n[irun]);
		
		fmax[ilayer-3]=histocharge->GetMaximum();
		if (jj==0)	{ fmaxold=fmax[ilayer-3];}
		if (jj!=0) { if(fmaxold<fmax[ilayer-3]) { fmaxold=fmax[ilayer-3];} }
		//j++;
	      }//end if histocharge
	    
	    histocharge=NULL;
	    
	    //drift time
	    
	    if(kUseOriginalFile)sprintf(filepath,"ITS/RecPoints/%s/%s_SDDdrifttime_Layer%d",especiename,especiename,ilayer);
	    else sprintf(filepath,"%s_SDDdrifttime_Layer%d",especiename,ilayer);
	    histocharge=(TH1F*)(mergedfile.Get(filepath));	 
	    
	    if(histocharge)
	      {
		sprintf(buffer, "%d", n[irun]);	
		fmaxtime=histocharge->GetMaximum();
		if (jj==0)	{ fmaxoldtime=fmaxtime;}
		if (jj!=0) { if(fmaxoldtime<fmaxtime) { fmaxoldtime=fmaxtime;} }
		//j++;
	      }//end if histocharge
	    jj++;
	    histocharge=NULL;
	    
	  }//end for run
	}//end layer

    fmaxmargin[0]=fmaxold;
    fmaxmargintime[0]=1.023*fmaxoldtime;
  

    fmaxmargin[1]=fmaxmargin[0];
    fmaxmargintime[1]=fmaxmargintime[0];


  for(Int_t layer=3;layer<5;layer++){
  
    printf("Max Y Range is %4.1f time %4.1f\n",fmaxmargin[layer-3],fmaxmargintime[layer-3]);  
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
		histocharge->SetMaximum(fmaxmargin[layer-3]);
		histocharge->GetYaxis()->SetTitle("");
		histocharge->DrawNormalized();
		canvas[1]->cd(layer-2)->Update();
		histocharge->GetYaxis()->SetRangeUser(0.,fmaxmargin[layer-3]);
		canvas[1]->cd(layer-2)->Update();
	      }	   
	    if (j!=0)
	      {
		//printf("new j %d\n",j);
		canvas[1]->cd(layer-2);
		histocharge->DrawNormalized("same");
		canvas[1]->cd(layer-2)->Update();
	      }	    
	    if(layer==3)legend->AddEntry(histocharge,runnmbr,"l");
	    canvas[1]->cd(layer-2);
	    legend->Draw("same");	    
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
		histocharge->SetMaximum(fmaxmargintime[layer-3]);
		histocharge->GetYaxis()->SetTitle("");
		histocharge->DrawNormalized();
	      }	   
	    if (hh!=0)
	      {
		canvas[2]->cd(layer-2);
		histocharge->DrawNormalized("same");
	      }	    
	    if(layer==3)legend2->AddEntry(histocharge,runnmbr,"l");
	    canvas[2]->cd(layer-2);
	    legend2->Draw("same");	    
	    canvas[2]->Update();
	    hh++;
	  }//end if histocharge
	else{printf("...Not Found....\n");}
	histocharge=NULL;
      }//enf for 
      }//end for layer
    }
  
  TFile trendfile(Form("SDDQAtrend%s%s.root",period,pass),"recreate");  
  trendfile.cd();
  for(Int_t ican=0;ican<10;ican++)canvas[ican]->Write();
  trendfile.Close();
  Char_t psfile[50];
  sprintf(psfile,"SDDtrend%s%s.ps",period,pass);
  canvas[0]->Print(Form("%s[",psfile));
  for(Int_t ifile=0;ifile<12;ifile++){canvas[ifile]->Print(psfile);}

  canvas[11]->Print(Form("%s]",psfile));
  
  
  delete histocharge;
  delete histochargerp;



}//end macro
