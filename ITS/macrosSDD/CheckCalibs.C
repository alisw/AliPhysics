#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TGLayout.h>
#include <TGWindow.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#endif

// Macro to display the output of SDD calibration runs
// -> copies the ASCII files with the DA output from the LDC
// -> plot various histograms and graphs
// Origin: F. Prino (prino@to.infn.it)

class CheckCalibInterface : public TGMainFrame {

private:
  TGCompositeFrame    *fHor1;
  TGCompositeFrame    *fHor2;
  TGCompositeFrame    *fHor3;
  TGCompositeFrame    *fHor4;
  TGTextButton        *fCopyFiles;
  TGTextButton        *fShowPedest;
  TGTextButton        *fShowPulser;
  TGTextButton        *fShowInject;
  TGTextButton        *fShowOneMod;
  TGTextButton        *fExit;
  TGGroupFrame        *fGframeALL;
  TGGroupFrame        *fGframeSING;
  TGGroupFrame        *fGframe1;
  TGGroupFrame        *fGframe2;
  TGNumberEntry       *fDDL;
  TGNumberEntry       *fChannel;
  Int_t fNumDDL;
  Int_t fNumChannel;

public:
  CheckCalibInterface(const TGWindow *p, UInt_t w, UInt_t h);
  virtual ~CheckCalibInterface();
  void DoSetlabel();

  static void CopyFiles();
  static void ShowPedestal();
  static void ShowPulser();
  static void ShowInjector();
  static void ShowSingleModule(Int_t iddl, Int_t ichan);
  static void ClearAll();

  ClassDef(CheckCalibInterface, 0)
};
                          
CheckCalibInterface::CheckCalibInterface(const TGWindow *p, UInt_t w, UInt_t h)
   : TGMainFrame(p, w, h)
{

  fHor1 = new TGHorizontalFrame(this);
  fHor2 = new TGHorizontalFrame(this);
  fHor3 = new TGHorizontalFrame(this);
  fHor4 = new TGHorizontalFrame(this);

  TGLayoutHints *lh=new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);
  TGLayoutHints *lhc=new TGLayoutHints(kLHintsCenterY | kLHintsExpandY, 5, 5, 5, 5);
  TGLayoutHints *lh1=new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5);

  fCopyFiles = new TGTextButton(fHor1, "&Copy Files", "CheckCalibInterface::CopyFiles()");
  fHor1->AddFrame(fCopyFiles, lh1);

  fGframeALL = new TGGroupFrame(fHor2,"All Modules",kHorizontalFrame);
  fShowPedest = new TGTextButton(fGframeALL, "&Show Pedestal","CheckCalibInterface::ShowPedestal()");
  fShowPulser = new TGTextButton(fGframeALL, "&Show Pulser","CheckCalibInterface::ShowPulser()");
  fShowInject = new TGTextButton(fGframeALL, "&Show Injector","CheckCalibInterface::ShowInjector()");    
  fGframeALL->AddFrame(fShowPedest, lh1);
  fGframeALL->AddFrame(fShowPulser, lh1);
  fGframeALL->AddFrame(fShowInject, lh1);
  fHor2->AddFrame(fGframeALL, lh1);

  fGframeSING = new TGGroupFrame(fHor3,"Single Module",kHorizontalFrame);
  fGframe1 = new TGGroupFrame(fGframeSING, "DDL");
  fDDL = new TGNumberEntry(fGframe1, 0, 2,-1, TGNumberFormat::kNESInteger,
			   TGNumberFormat::kNEANonNegative, 
			   TGNumberFormat::kNELLimitMinMax,
			   0, 23);
  fGframe2 = new TGGroupFrame(fGframeSING, "Channel");
  fChannel = new TGNumberEntry(fGframe2, 0, 2,-1, TGNumberFormat::kNESInteger,
			       TGNumberFormat::kNEANonNegative, 
			       TGNumberFormat::kNELLimitMinMax,
			       0, 11);
  fGframe1->AddFrame(fDDL,lh);
  fGframe2->AddFrame(fChannel,lh);
  fGframeSING->AddFrame(fGframe1,lh);
  fGframeSING->AddFrame(fGframe2,lh);
  
  fDDL->Connect("ValueSet(Long_t)", "CheckCalibInterface", this, "DoSetlabel()");
  (fDDL->GetNumberEntry())->Connect("ReturnPressed()", "CheckCalibInterface", this, "DoSetlabel()");
  
  fChannel->Connect("ValueSet(Long_t)", "CheckCalibInterface", this, "DoSetlabel()");
  (fChannel->GetNumberEntry())->Connect("ReturnPressed()", "CheckCalibInterface", this, "DoSetlabel()");
  fHor3->AddFrame(fGframeSING, lh1);
  
  fShowOneMod = new TGTextButton(fGframeSING, "&Show Selected Module","CheckCalibInterface::ShowSingleModule(0,0)");
  fGframeSING->AddFrame(fShowOneMod, lhc);

  fExit = new TGTextButton(fHor4, "&Exit", "gApplication->Terminate(0)");
  fHor4->AddFrame(fExit, lh1);

  AddFrame(fHor1,lh1);
  AddFrame(fHor2,lh1);
  AddFrame(fHor3,lh1);
  AddFrame(fHor4,lh1);

  SetCleanup(kDeepCleanup);
  SetWindowName("Main Control");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();
}

CheckCalibInterface::~CheckCalibInterface()
{
   // Destructor.
   
   Cleanup();
}


void CheckCalibInterface::CopyFiles(){
  //
  TString command;
  TString ldcName[6]={"aldaqpc083","aldaqpc084","aldaqpc041",
		      "aldaqpc082","aldaqpc085","aldaqpc086"};
  gSystem->Exec("rm -rf calibFiles");
  gSystem->Exec("mkdir calibFiles");
  for(Int_t iLDC=0; iLDC<6; iLDC++){
    Int_t firstDDL=iLDC*4;
    Int_t lastDDL=iLDC*4+3;
    command.Form("scp %s:/dateSite/ldc-SDD-%02d-%02d-0/work/SDDbase_step2_LDC.tar calibFiles/.",ldcName[iLDC].Data(),firstDDL,lastDDL);
    gSystem->Exec(command.Data());
    command.Form("scp %s:/dateSite/ldc-SDD-%02d-%02d-0/work/SDDbase_LDC.tar calibFiles/.",ldcName[iLDC].Data(),firstDDL,lastDDL);
    gSystem->Exec(command.Data());
    command.Form("scp %s:/dateSite/ldc-SDD-%02d-%02d-0/work/SDDinj_LDC.tar calibFiles/.",ldcName[iLDC].Data(),firstDDL,lastDDL);
    gSystem->Exec(command.Data());
    gSystem->Exec("cd calibFiles; tar xvf SDDbase_step2_LDC.tar; tar xvf SDDbase_LDC.tar; tar xvf SDDinj_LDC.tar; cd ..");
  }
  printf("-------------- DONE ---------------\n");

  return;
}

void CheckCalibInterface::ShowPedestal(){
  ClearAll();

  TH2F* hbadchan=new TH2F("hbadchan","Number of bad channels",24,-0.25,11.75,24,-0.5,23.5);
  TH2F* hrawnoisemod=new TH2F("hrawnoisemod","",288,-0.5,287.5,50,0.,10.);
  TH1F* hrawnoise=new TH1F("hrawnoise","",100,0.,10.);
  TH2F* hcornoisemod=new TH2F("hcornoisemod","",288,-0.5,287.5,50,0.,10.);
  TH1F* hcornoise=new TH1F("hcornoise","",100,0.,10.);
  TH2F* hbasemod=new TH2F("hbasemod","",288,-0.5,287.5,50,0.,150.);
  TH1F* hbase=new TH1F("hbase","",100,0.,150.);
  Int_t retfscf;
  TString inpFileName;
  Float_t baseline,rawnoise,cmn,corn;
  Int_t isgoodan,i,basmin,basoff;
  Int_t th,tl;
  Int_t nGoodAnodes=0;
  for(Int_t iddl=0;iddl<24;iddl++){
    for(Int_t imod=0;imod<12;imod++){
      for(Int_t isid=0;isid<2;isid++){
	inpFileName.Form("./calibFiles/SDDbase_step2_ddl%02dc%02d_sid%d.data",iddl,imod,isid);
	FILE* basFil = fopen(inpFileName.Data(),"read");
	Int_t sideId=imod*2+isid;
	Int_t modId=iddl*12+imod;
	if (basFil == 0) hbadchan->SetBinContent(sideId+1,iddl+1,256);
	else{
	  retfscf=fscanf(basFil,"%d\n",&th);
	  retfscf=fscanf(basFil,"%d\n",&tl);
	  if(th==255 && tl==20) continue;
	  for(Int_t ian=0;ian<256;ian++){
	    retfscf=fscanf(basFil,"%d %d %f %d %d %f %f %f\n",&i,&isgoodan,&baseline,&basmin,&basoff,&rawnoise,&cmn,&corn);
	    hrawnoisemod->Fill(modId,rawnoise);
	    hrawnoise->Fill(rawnoise);
	    hcornoisemod->Fill(modId,corn);
	    hcornoise->Fill(corn);
	    hbasemod->Fill(modId,baseline);
	    hbase->Fill(baseline);
	    if(!isgoodan){
	      hbadchan->SetBinContent(sideId+1,iddl+1, 1+hbadchan->GetBinContent(sideId+1,iddl+1));
	    }else{
	      nGoodAnodes++;
	    }
	  }
	  fclose(basFil);
	}
      }
    }
  }
  hrawnoisemod->SetStats(0);
  hcornoisemod->SetStats(0);
  hbasemod->SetStats(0);
  hbadchan->SetStats(0);
  gStyle->SetPalette(1);
  TString txtCountGood;
  Float_t fracGood=100.*(Float_t)nGoodAnodes/(520.*260.);
  txtCountGood.Form("Number of GoodAnodes = %d (%5.1f\%)\n",nGoodAnodes,fracGood);

  TCanvas* c0=new TCanvas("c0","Bad Channels",800,900);
  c0->SetRightMargin(0.14);
  c0->SetBottomMargin(0.2);
  hbadchan->Draw("colz"); 
  hbadchan->GetXaxis()->SetTitle("Channel");
  hbadchan->GetYaxis()->SetTitle("DDL");
  hbadchan->GetXaxis()->SetTickLength(0);
  hbadchan->GetYaxis()->SetTickLength(0);
  c0->cd();
  TLine** linv=new TLine*[12];
  for(Int_t i=0;i<12;i++){
    linv[i]=new TLine(i+0.75,-0.5,i+0.75,23.5);
    linv[i]->SetLineColor(kGray+1);
    linv[i]->Draw();
  }
  TLine** linh=new TLine*[24];
  for(Int_t i=0;i<24;i++){
    linh[i]=new TLine(-0.25,i+0.5,11.75,i+0.5);
    linh[i]->SetLineColor(kGray+1);
    linh[i]->Draw();
  }
  TLatex* tg=new TLatex(0.1,0.05,txtCountGood.Data());
  tg->SetNDC();
  tg->SetTextColor(4);
  tg->SetTextSize(0.04);
  tg->Draw();
  c0->Update();
 
  TCanvas* c1=new TCanvas("c1","Baseline",1200,700);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hbase->Draw();
  hbase->GetXaxis()->SetTitle("Baseline (ADC)");
  c1->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hbasemod->Draw("colz");
  hbasemod->GetXaxis()->SetTitle("Module index (=DDL*12+Channel)");
  hbasemod->GetYaxis()->SetTitle("Baseline (ADC)");
  hbasemod->GetYaxis()->SetTitleOffset(1.35);

  TCanvas* c2=new TCanvas("c2","Raw Noise",1200,700);
  c2->Divide(2,1);
  c2->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hrawnoise->Draw();
  hrawnoise->GetXaxis()->SetTitle("Raw Noise (ADC)");
  c2->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hrawnoisemod->Draw("colz");
  hrawnoisemod->GetXaxis()->SetTitle("Module index (=DDL*12+Channel)");
  hrawnoisemod->GetYaxis()->SetTitle("Raw Noise (ADC)");
  hrawnoisemod->GetYaxis()->SetTitleOffset(1.35);

  TCanvas* c3=new TCanvas("c3","Corr Noise",1200,700);
  c3->Divide(2,1);
  c3->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hcornoise->Draw();
  hcornoise->GetXaxis()->SetTitle("Raw Noise (ADC)");
  c3->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hcornoisemod->Draw("colz");
  hcornoisemod->GetXaxis()->SetTitle("Module index (=DDL*12+Channel)");
  hcornoisemod->GetYaxis()->SetTitle("Raw Noise (ADC)");
  hcornoisemod->GetYaxis()->SetTitleOffset(1.35);  
}

void CheckCalibInterface::ShowPulser(){
  ClearAll();

  TH2F* hbadchan=new TH2F("hbadchan","Number of bad channels",24,-0.25,11.75,24,-0.5,23.5);
  TH2F* hrawnoisemod=new TH2F("hrawnoisemod","",288,-0.5,287.5,50,0.,10.);
  TH1F* hrawnoise=new TH1F("hrawnoise","",100,0.,10.);
  TH2F* hcornoisemod=new TH2F("hcornoisemod","",288,-0.5,287.5,50,0.,10.);
  TH1F* hcornoise=new TH1F("hcornoise","",100,0.,10.);
  TH2F* hbasemod=new TH2F("hbasemod","",288,-0.5,287.5,50,0.,150.);
  TH1F* hbase=new TH1F("hbase","",100,0.,150.);
  TH2F* hgainmod=new TH2F("hgainmod","",288,-0.5,287.5,50,0.,5.);
  TH1F* hgain=new TH1F("hgain","",100,0.,5.);

  Int_t retfscf;
  TString inpFileName;
  Float_t baseline,rawnoise,cmn,corn,gain;
  Int_t isgoodan,i,basmin,basoff;
  Int_t th,tl;
  Int_t iii,jjj,kkk;
  Int_t nGoodAnodes=0;

  for(Int_t iddl=0;iddl<24;iddl++){
    for(Int_t imod=0;imod<12;imod++){
      for(Int_t isid=0;isid<2;isid++){
	inpFileName.Form("./calibFiles/SDDbase_ddl%02dc%02d_sid%d.data",iddl,imod,isid);
	FILE* basFil = fopen(inpFileName.Data(),"read");
	Int_t sideId=imod*2+isid;
	Int_t modId=iddl*12+imod;
	if (basFil == 0) hbadchan->SetBinContent(sideId+1,iddl+1,256);
	else{
	  retfscf=fscanf(basFil,"%d %d %d\n",&iii,&jjj,&kkk);
	  if(kkk==0) continue;
	  retfscf=fscanf(basFil,"%d\n",&th);
	  retfscf=fscanf(basFil,"%d\n",&tl);
	  if(th==255 && tl==20) continue;
	  for(Int_t ian=0;ian<256;ian++){
	    retfscf=fscanf(basFil,"%d %d %f %d %d %f %f %f %f\n",&i,&isgoodan,&baseline,&basmin,&basoff,&rawnoise,&cmn,&corn,&gain);
	    hrawnoisemod->Fill(modId,rawnoise);
	    hrawnoise->Fill(rawnoise);
	    hcornoisemod->Fill(modId,corn);
	    hcornoise->Fill(corn);
	    hbasemod->Fill(modId,baseline);
	    hbase->Fill(baseline);
	    hgainmod->Fill(modId,gain);
	    hgain->Fill(gain);
	    if(!isgoodan){
	      hbadchan->SetBinContent(sideId+1,iddl+1, 1+hbadchan->GetBinContent(sideId+1,iddl+1));
	    }else{
	      nGoodAnodes++;
	    }
	  }
	  fclose(basFil);
	}
      }
    }
  }
  hrawnoisemod->SetStats(0);
  hcornoisemod->SetStats(0);
  hbasemod->SetStats(0);
  hgainmod->SetStats(0);
  hbadchan->SetStats(0);
  gStyle->SetPalette(1);
  TString txtCountGood;
  Float_t fracGood=100.*(Float_t)nGoodAnodes/(520.*260.);
  txtCountGood.Form("Number of GoodAnodes = %d (%5.1f\%)\n",nGoodAnodes,fracGood);

  TCanvas* c0=new TCanvas("c0","Bad Channels",800,900);
  c0->SetRightMargin(0.14);
  c0->SetBottomMargin(0.2);
  hbadchan->Draw("colz"); 
  hbadchan->GetXaxis()->SetTitle("Channel");
  hbadchan->GetYaxis()->SetTitle("DDL");
  hbadchan->GetXaxis()->SetTickLength(0);
  hbadchan->GetYaxis()->SetTickLength(0);
  c0->cd();
  TLine** linv=new TLine*[12];
  for(Int_t i=0;i<12;i++){
    linv[i]=new TLine(i+0.75,-0.5,i+0.75,23.5);
    linv[i]->SetLineColor(kGray+1);
    linv[i]->Draw();
  }
  TLine** linh=new TLine*[24];
  for(Int_t i=0;i<24;i++){
    linh[i]=new TLine(-0.25,i+0.5,11.75,i+0.5);
    linh[i]->SetLineColor(kGray+1);
    linh[i]->Draw();
  }
  TLatex* tg=new TLatex(0.1,0.05,txtCountGood.Data());
  tg->SetNDC();
  tg->SetTextColor(4);
  tg->SetTextSize(0.04);
  tg->Draw();
  c0->Update();
 
  TCanvas* c1=new TCanvas("c1","Baseline",1200,700);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hbase->Draw();
  hbase->GetXaxis()->SetTitle("Baseline (ADC)");
  c1->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hbasemod->Draw("colz");
  hbasemod->GetXaxis()->SetTitle("Module index (=DDL*12+Channel)");
  hbasemod->GetYaxis()->SetTitle("Baseline (ADC)");
  hbasemod->GetYaxis()->SetTitleOffset(1.35);

  TCanvas* c2=new TCanvas("c2","Raw Noise",1200,700);
  c2->Divide(2,1);
  c2->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hrawnoise->Draw();
  hrawnoise->GetXaxis()->SetTitle("Raw Noise (ADC)");
  c2->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hrawnoisemod->Draw("colz");
  hrawnoisemod->GetXaxis()->SetTitle("Module index (=DDL*12+Channel)");
  hrawnoisemod->GetYaxis()->SetTitle("Raw Noise (ADC)");
  hrawnoisemod->GetYaxis()->SetTitleOffset(1.35);

  TCanvas* c3=new TCanvas("c3","Corr Noise",1200,700);
  c3->Divide(2,1);
  c3->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hcornoise->Draw();
  hcornoise->GetXaxis()->SetTitle("Raw Noise (ADC)");
  c3->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hcornoisemod->Draw("colz");
  hcornoisemod->GetXaxis()->SetTitle("Module index (=DDL*12+Channel)");
  hcornoisemod->GetYaxis()->SetTitle("Raw Noise (ADC)");
  hcornoisemod->GetYaxis()->SetTitleOffset(1.35);

  TCanvas* c4=new TCanvas("c4","Gain",1200,700);
  c4->Divide(2,1);
  c4->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hgain->Draw();
  hgain->GetXaxis()->SetTitle("Gain (ADC/DAC)");
  c4->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.14);
  hgainmod->Draw("colz");
  hgainmod->GetXaxis()->SetTitle("Module index (=DDL*12+Channel)");
  hgainmod->GetYaxis()->SetTitle("Gain (ADC/DAC)");
  hgainmod->GetYaxis()->SetTitleOffset(1.35);

  return;
}

void CheckCalibInterface::ShowInjector(){
  ClearAll();

  TH2F* hinjstatus=new TH2F("hinjstatus","Injector Status",24,-0.25,11.75,24,-0.5,23.5);
  TGraph *vvsmod0=new TGraph(0);
  TGraph *vvsmod1=new TGraph(0);
  TGraph *poldegvsmod0=new TGraph(0); 
  TGraph *poldegvsmod1=new TGraph(0); 
  TGraph *anmaxvsmod0=new TGraph(0); 
  TGraph *anmaxvsmod1=new TGraph(0); 
  TGraph *dvcevsmod0=new TGraph(0);
  TGraph *dvcevsmod1=new TGraph(0);
  TGraph *dveevsmod0=new TGraph(0);
  TGraph *dveevsmod1=new TGraph(0);
  vvsmod0->SetTitle("Drift Speed vs. mod. number");
  vvsmod1->SetTitle("Drift Speed vs. mod. number");
  poldegvsmod0->SetTitle("Degree of poly fit vs. mod. number");
  poldegvsmod1->SetTitle("Degree of poly fit vs. mod. number");
  anmaxvsmod0->SetTitle("Anode with max. vdrift vs. mod. number");
  anmaxvsmod1->SetTitle("Anode with max. vdrift vs. mod. number");
  dvcevsmod0->SetTitle("Delta Vdrift 128-0 vs. mod. number");
  dvcevsmod1->SetTitle("Delta Vdrift 128-0 vs. mod. number");
  dveevsmod0->SetTitle("Delta Vdrift 256-0 vs. mod. number");
  dveevsmod1->SetTitle("Delta Vdrift 256-0 vs. mod. number");

  TF1* fPoly=new TF1("fPoly","pol3",0.,256.);

  Int_t evNumb,polDeg; 
  UInt_t timeStamp,statusInj;
  Int_t retfscf;
  TString inpFileName;
  Float_t auxP;
  Int_t iGoodInj=0;
  for(Int_t iddl=0;iddl<24;iddl++){
    for(Int_t imod=0;imod<12;imod++){
      for(Int_t isid=0;isid<2;isid++){
	inpFileName.Form("./calibFiles/SDDinj_ddl%02dc%02d_sid%d.data",iddl,imod,isid);
	FILE* injFil = fopen(inpFileName.Data(),"read");
	Int_t sideId=imod*2+isid;
	Int_t modId=iddl*12+imod;
	if (injFil == 0){ 
	  hinjstatus->SetBinContent(sideId+1,iddl+1,0);
	}else{
	  Bool_t firstEvent=kTRUE;
	  retfscf=fscanf(injFil,"%d",&polDeg);
	  while (!feof(injFil)){
	    retfscf=fscanf(injFil,"%d %u ",&evNumb,&timeStamp);
	    if(evNumb==-99){
	      statusInj=timeStamp;
	      Int_t n7=(statusInj&(0x1F<<25))>>25;
	      Int_t n6=(statusInj&(0x1F<<20))>>20;
	      Int_t n5=(statusInj&(0x1F<<15))>>15;
	      Int_t n4=(statusInj&(0x1F<<5))>>10;
	      Int_t n3=(statusInj&(0x1F<<5))>>5;
	      Int_t n2=statusInj&0x1F;
	      Float_t aveStatus=(7.*n7+6.*n6+5.*n5+4.*n4+3.*n3+2.*n2)/32;
	      hinjstatus->SetBinContent(sideId+1,iddl+1,aveStatus);
	    }else{
	      if(feof(injFil)) break;
	      for(Int_t ic=0;ic<4;ic++){ 
		retfscf=fscanf(injFil,"%f ",&auxP);
		fPoly->SetParameter(ic,auxP);		
	      }
	      if(firstEvent && polDeg>0){
		firstEvent=kFALSE;
		++iGoodInj;
		if(isid==0){
		  vvsmod0->SetPoint(vvsmod0->GetN(),(Float_t)modId,fPoly->Eval(128));
		  poldegvsmod0->SetPoint(poldegvsmod0->GetN(),(Float_t)modId,polDeg);
		  anmaxvsmod0->SetPoint(anmaxvsmod0->GetN(),(Float_t)modId,fPoly->GetMaximumX(0.,256.));
		  dvcevsmod0->SetPoint(dvcevsmod0->GetN(),(Float_t)modId,fPoly->Eval(128)-fPoly->Eval(0));
		  dveevsmod0->SetPoint(dveevsmod0->GetN(),(Float_t)modId,fPoly->Eval(256)-fPoly->Eval(0));
		}else{
		  vvsmod1->SetPoint(vvsmod1->GetN(),(Float_t)modId,fPoly->Eval(128));
		  poldegvsmod1->SetPoint(poldegvsmod1->GetN(),(Float_t)modId,polDeg);
		  anmaxvsmod1->SetPoint(anmaxvsmod1->GetN(),(Float_t)modId,fPoly->GetMaximumX(0.,256.));
		  dvcevsmod1->SetPoint(dvcevsmod1->GetN(),(Float_t)modId,fPoly->Eval(128)-fPoly->Eval(0));
		  dveevsmod1->SetPoint(dveevsmod1->GetN(),(Float_t)modId,fPoly->Eval(256)-fPoly->Eval(0));
		}
	      }
	    }
	  }	
	  fclose(injFil);
	}
      }
    }
  }
  delete fPoly;

  TString countmods;
  countmods.Form("Number of half-modules with drift speed from injectors = %d",iGoodInj);
  gStyle->SetPalette(59);
  hinjstatus->SetStats(0);
  hinjstatus->SetMinimum(-0.01);
  hinjstatus->SetMaximum(7.);
  TCanvas* c0=new TCanvas("c0","Injector status",800,900);
  c0->SetRightMargin(0.14);
  c0->SetBottomMargin(0.2);
  hinjstatus->Draw("colz"); 
  hinjstatus->GetXaxis()->SetTitle("Channel");
  hinjstatus->GetYaxis()->SetTitle("DDL");
  hinjstatus->GetXaxis()->SetTickLength(0);
  hinjstatus->GetYaxis()->SetTickLength(0);
  c0->cd();
  TLine** linv=new TLine*[12];
  for(Int_t i=0;i<12;i++){
    linv[i]=new TLine(i+0.75,-0.5,i+0.75,23.5);
    linv[i]->SetLineColor(kGray+1);
    linv[i]->Draw();
  }
  TLine** linh=new TLine*[24];
  for(Int_t i=0;i<24;i++){
    linh[i]=new TLine(-0.25,i+0.5,11.75,i+0.5);
    linh[i]->SetLineColor(kGray+1);
    linh[i]->Draw();
  }
  TLatex* t3=new TLatex(0.1,0.05,countmods.Data());
  t3->SetNDC();
  t3->SetTextColor(4);
  t3->SetTextSize(0.03);
  t3->Draw();
  c0->Update();

  TLatex* tleft=new TLatex(0.2,0.82,"Side 0");
  tleft->SetNDC();
  tleft->SetTextColor(1);
  TLatex* tright=new TLatex(0.2,0.75,"Side 1");
  tright->SetNDC();
  tright->SetTextColor(2);

  TCanvas* c1;
  c1=new TCanvas("c1","Vdrift vs. mod",1000,700);  
  vvsmod0->SetMarkerStyle(20);
  vvsmod0->Draw("AP");
  vvsmod0->GetXaxis()->SetLimits(-1,290);
  vvsmod0->GetXaxis()->SetTitle("Module index (=DDL*12+Channel)");
  vvsmod0->GetYaxis()->SetTitle("Vdrift (#mum/ns)");
  vvsmod1->SetMarkerStyle(21);
  vvsmod1->SetMarkerColor(2);
  vvsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();
  TLine* ltop=new TLine(12*8-0.5,vvsmod0->GetYaxis()->GetXmin(),12*8-0.5,vvsmod0->GetYaxis()->GetXmax());
  ltop->SetLineColor(4);
  ltop->Draw();
  TLatex* ttop=new TLatex(12*3.5,vvsmod0->GetYaxis()->GetXmin()+0.05,"TOP");
  ttop->SetTextColor(4);
  ttop->Draw();
  TLine* lmed=new TLine(12*16-0.5,vvsmod0->GetYaxis()->GetXmin(),12*16-0.5,vvsmod0->GetYaxis()->GetXmax());
  lmed->SetLineColor(4);
  lmed->Draw();
  TLatex* tmed=new TLatex(12*11.5,vvsmod0->GetYaxis()->GetXmin()+0.05,"MED");
  tmed->SetTextColor(4);
  tmed->Draw();
  TLatex* tbot=new TLatex(12*19.5,vvsmod0->GetYaxis()->GetXmin()+0.05,"BOT");
  tbot->SetTextColor(4);
  tbot->Draw();


  TCanvas* c2=new TCanvas("c2","Params vs. mod",900,900);
  c2->Divide(2,2);
  
  c2->cd(1);
  gPad->SetLeftMargin(0.14);
  poldegvsmod0->SetMarkerStyle(20);
  poldegvsmod0->Draw("AP");
  poldegvsmod0->GetXaxis()->SetLimits(-1,290);
  poldegvsmod0->GetXaxis()->SetTitle("Module index (=DDL*12+Channel)");
  poldegvsmod0->GetYaxis()->SetTitle("Degree of Polynomial fit");
  poldegvsmod0->GetYaxis()->SetTitleOffset(1.4);
  poldegvsmod1->SetMarkerStyle(21);
  poldegvsmod1->SetMarkerColor(2);
  poldegvsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();
  c2->cd(2);
  gPad->SetLeftMargin(0.14);
  anmaxvsmod0->SetMarkerStyle(20);
  anmaxvsmod0->Draw("AP");
  anmaxvsmod0->GetXaxis()->SetLimits(-1,290);
  anmaxvsmod0->GetXaxis()->SetTitle("Module index (=DDL*12+Channel)");
  anmaxvsmod0->GetYaxis()->SetTitle("Anode with max. drift speed");
  anmaxvsmod0->GetYaxis()->SetTitleOffset(1.4);
  anmaxvsmod1->SetMarkerStyle(21);
  anmaxvsmod1->SetMarkerColor(2);
  anmaxvsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();
  c2->cd(3);
  gPad->SetLeftMargin(0.14);
  dvcevsmod0->SetMarkerStyle(20);
  dvcevsmod0->Draw("AP");
  dvcevsmod0->GetXaxis()->SetLimits(-1,290);
  dvcevsmod0->GetXaxis()->SetTitle("Module index (=DDL*12+Channel)");
  dvcevsmod0->GetYaxis()->SetTitle("vdrift(anode128)-vdrift(anode0)");
  dvcevsmod0->GetYaxis()->SetTitleOffset(1.4);
  dvcevsmod1->SetMarkerStyle(21);
  dvcevsmod1->SetMarkerColor(2);
  dvcevsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();
  c2->cd(4);
  gPad->SetLeftMargin(0.14);
  dveevsmod0->SetMarkerStyle(20);
  dveevsmod0->Draw("AP");
  dveevsmod0->GetXaxis()->SetLimits(-1,290);
  dveevsmod0->GetYaxis()->SetTitleOffset(1.4);
  dveevsmod0->GetXaxis()->SetTitle("Module index (=DDL*12+Channel)");
  dveevsmod0->GetYaxis()->SetTitle("vdrift(anode256)-vdrift(anode0)");
  dveevsmod1->SetMarkerStyle(21);
  dveevsmod1->SetMarkerColor(2);
  dveevsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();

  return;
}

void CheckCalibInterface::ShowSingleModule(Int_t iddl, Int_t ichan){
  //
  ClearAll();
  TString inpFileName1,inpFileName2,inpFileName3;
  Int_t retfscf;
  Float_t baseline,rawnoise,cmn,corn,gain;
  Int_t isgoodan,i,basmin,basoff;
  Int_t th,tl;
  Int_t iii,jjj,kkk;
  Int_t evNumb,polDeg; 
  UInt_t timeStamp,statusInj;
  Float_t auxP;
  TF1* fPoly=new TF1("fPoly","pol3",0.,256.);

  TGraph* gbasel=new TGraph(0);
  TGraph* gbaser=new TGraph(0);
  TGraph* grawnl=new TGraph(0);
  TGraph* grawnr=new TGraph(0);
  TGraph* gcornl=new TGraph(0);
  TGraph* gcornr=new TGraph(0);
  TGraph* ggainl=new TGraph(0);
  TGraph* ggainr=new TGraph(0);
  TGraph* gstpdl=new TGraph(0);
  TGraph* gstpdr=new TGraph(0);
  TGraph* gstpul=new TGraph(0);
  TGraph* gstpur=new TGraph(0);
  TGraph* gdrspl=new TGraph(0);
  TGraph* gdrspr=new TGraph(0);
  gbasel->SetTitle("Baseline Left");
  gbaser->SetTitle("Baseline Right");
  grawnl->SetTitle("Noise Left");
  grawnr->SetTitle("Noise Right");
  gcornl->SetTitle("Noise Left");
  gcornr->SetTitle("Noise Right");
  ggainl->SetTitle("Gain Left");
  ggainr->SetTitle("Gain Right");
  gstpdl->SetTitle("Status Left");
  gstpdr->SetTitle("Status Right");
  gstpul->SetTitle("Status Left");
  gstpur->SetTitle("Status Right");
  gdrspl->SetTitle("Drift Speed Left");
  gdrspr->SetTitle("Drift Speed Right");

  for(Int_t isid=0;isid<2;isid++){
    inpFileName1.Form("./calibFiles/SDDbase_step2_ddl%02dc%02d_sid%d.data",iddl,ichan,isid);
    inpFileName2.Form("./calibFiles/SDDbase_ddl%02dc%02d_sid%d.data",iddl,ichan,isid);
    inpFileName3.Form("./calibFiles/SDDinj_ddl%02dc%02d_sid%d.data",iddl,ichan,isid);
  
  
    FILE* basFil = fopen(inpFileName1.Data(),"read");
    if (basFil != 0){
      retfscf=fscanf(basFil,"%d\n",&th);
      retfscf=fscanf(basFil,"%d\n",&tl);
      if(th==255 && tl==20) continue;
      for(Int_t ian=0;ian<256;ian++){
	retfscf=fscanf(basFil,"%d %d %f %d %d %f %f %f\n",&i,&isgoodan,&baseline,&basmin,&basoff,&rawnoise,&cmn,&corn);
	if(isid==0){
	  gbasel->SetPoint(gbasel->GetN(),(Float_t)i,baseline);
	  grawnl->SetPoint(grawnl->GetN(),(Float_t)i,rawnoise);
	  gcornl->SetPoint(gcornl->GetN(),(Float_t)i,corn);
	  gstpdl->SetPoint(gstpdl->GetN(),(Float_t)i,(Float_t)isgoodan);
	}else{
	  gbaser->SetPoint(gbaser->GetN(),(Float_t)i,baseline);
	  grawnr->SetPoint(grawnr->GetN(),(Float_t)i,rawnoise);
	  gcornr->SetPoint(gcornr->GetN(),(Float_t)i,corn);
	  gstpdr->SetPoint(gstpdr->GetN(),(Float_t)i,(Float_t)isgoodan);
	}
      }
      fclose(basFil);
    }
	
    FILE* pulFil = fopen(inpFileName2.Data(),"read");
    if (pulFil != 0){
      retfscf=fscanf(pulFil,"%d %d %d\n",&iii,&jjj,&kkk);
      retfscf=fscanf(pulFil,"%d\n",&th);
      retfscf=fscanf(pulFil,"%d\n",&tl);
      if(th==255 && tl==20) continue;
      for(Int_t ian=0;ian<256;ian++){
	retfscf=fscanf(pulFil,"%d %d %f %d %d %f %f %f %f\n",&i,&isgoodan,&baseline,&basmin,&basoff,&rawnoise,&cmn,&corn,&gain);
	if(isid==0){
	  ggainl->SetPoint(ggainl->GetN(),(Float_t)i,gain);
	  gstpul->SetPoint(gstpul->GetN(),(Float_t)i,(Float_t)isgoodan);
	}else{
	  ggainr->SetPoint(ggainr->GetN(),(Float_t)i,gain);
	  gstpur->SetPoint(gstpur->GetN(),(Float_t)i,(Float_t)isgoodan);
	}
      }
      fclose(pulFil);
    }
    
    FILE* injFil = fopen(inpFileName3.Data(),"read");
    Bool_t firstEvent=kTRUE;
    if (injFil != 0){
      retfscf=fscanf(injFil,"%d",&polDeg);
      while (!feof(injFil)){
	retfscf=fscanf(injFil,"%d %u ",&evNumb,&timeStamp);
	if(evNumb==-99){
	  statusInj=timeStamp;
	}else{
	  if(feof(injFil)) break;
	  for(Int_t ic=0;ic<4;ic++){ 
	    retfscf=fscanf(injFil,"%f ",&auxP);
	    fPoly->SetParameter(ic,auxP);
	  }	  
	}
	if(firstEvent==kTRUE && polDeg>0){
	  firstEvent=kFALSE;
	  for(Int_t ian=0; ian<256; ian+=8){
	    if(isid==0) gdrspl->SetPoint(gdrspl->GetN(),(Float_t)ian,fPoly->Eval(ian));
	    else gdrspr->SetPoint(gdrspr->GetN(),(Float_t)ian,fPoly->Eval(ian));
	    
	  }
	}
      }
    }
  }


  TCanvas * c0=new TCanvas("c0","Baselines",900,600);
  c0->Divide(2,1);
  c0->cd(1);
  gPad->SetLeftMargin(0.14);
  gbasel->SetMarkerStyle(7);
  gbasel->Draw("AP");
  gbasel->GetXaxis()->SetTitle("Anode");
  gbasel->GetYaxis()->SetTitle("Baseline (ADC)");
  gbasel->GetYaxis()->SetTitleOffset(1.35);
  c0->cd(2);
  gPad->SetLeftMargin(0.14);
  gbaser->SetMarkerStyle(7);
  gbaser->Draw("AP");
  gbaser->GetXaxis()->SetTitle("Anode");
  gbaser->GetYaxis()->SetTitle("Baseline (ADC)");
  gbaser->GetYaxis()->SetTitleOffset(1.35);

  TCanvas * c1=new TCanvas("c1","Noise",900,600);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLeftMargin(0.14);
  TLatex* t1=new TLatex(0.6,0.8,"Raw");
  t1->SetNDC();
  TLatex* t2=new TLatex(0.6,0.75,"Corrected");
  t2->SetTextColor(2);
  t2->SetNDC();
  grawnl->SetMarkerStyle(7);
  grawnl->Draw("AP");
  grawnl->GetXaxis()->SetTitle("Anode");
  grawnl->GetYaxis()->SetTitle("Noise (ADC)");
  grawnl->GetYaxis()->SetTitleOffset(1.35);
  gcornl->SetMarkerStyle(7);
  gcornl->SetMarkerColor(2);
  gcornl->Draw("SAMEP");
  t1->Draw();
  t2->Draw();
  c1->cd(2);
  gPad->SetLeftMargin(0.14);
  grawnr->SetMarkerStyle(7);
  grawnr->Draw("AP");
  grawnr->GetXaxis()->SetTitle("Anode");
  grawnr->GetYaxis()->SetTitle("Noise (ADC)");
  grawnr->GetYaxis()->SetTitleOffset(1.35);
  gcornr->SetMarkerStyle(7);
  gcornr->SetMarkerColor(2);
  gcornr->Draw("SAMEP");
  t1->Draw();
  t2->Draw();
 
  TCanvas * c2=new TCanvas("c2","Gain",900,600);
  c2->Divide(2,1);
  c2->cd(1);
  gPad->SetLeftMargin(0.14);
  ggainl->SetMarkerStyle(7);
  ggainl->Draw("AP");
  ggainl->GetXaxis()->SetTitle("Anode");
  ggainl->GetYaxis()->SetTitle("Gain (ADC/DAC)");
  ggainl->GetYaxis()->SetTitleOffset(1.35);
  c2->cd(2);
  gPad->SetLeftMargin(0.14);
  ggainr->SetMarkerStyle(7);
  ggainr->Draw("AP");
  ggainr->GetXaxis()->SetTitle("Anode");
  ggainr->GetYaxis()->SetTitle("Gain (ADC/DAC)");
  ggainr->GetYaxis()->SetTitleOffset(1.35);

  TCanvas * c3=new TCanvas("c3","Anode Status",900,600);
  c3->Divide(2,1);
  c3->cd(1);
  gPad->SetLeftMargin(0.14);
  TLatex* t3=new TLatex(0.6,0.85,"Pedestal");
  t3->SetNDC();
  TLatex* t4=new TLatex(0.6,0.8,"Pulser");
  t4->SetTextColor(4);
  t4->SetNDC();
  gstpdl->SetMarkerStyle(7);
  gstpdl->Draw("AP");
  gstpdl->SetMinimum(-0.001);
  gstpdl->SetMaximum(1.2);
  gstpdl->GetXaxis()->SetTitle("Anode");
  gstpdl->GetYaxis()->SetTitle("Status (1=OK, 0=BAD)");
  gstpdl->GetYaxis()->SetTitleOffset(1.35);
  gstpul->SetMarkerStyle(7);
  gstpul->SetMarkerColor(4);
  gstpul->Draw("SAMEP");
  t3->Draw();
  t4->Draw();
  c3->cd(2);
  gPad->SetLeftMargin(0.14);
  gstpdr->SetMarkerStyle(7);
  gstpdr->Draw("AP");
  gstpdr->SetMinimum(-0.001);
  gstpdr->SetMaximum(1.2);
  gstpdr->GetXaxis()->SetTitle("Anode");
  gstpdr->GetYaxis()->SetTitle("Status (1=OK, 0=BAD)");
  gstpdr->GetYaxis()->SetTitleOffset(1.35);
  gstpur->SetMarkerStyle(7);
  gstpur->SetMarkerColor(4);
  gstpur->Draw("SAMEP");
  t3->Draw();
  t4->Draw();

  TCanvas * c4=new TCanvas("c4","Drift Speed",900,600);
  c4->Divide(2,1);
  c4->cd(1);
  gPad->SetLeftMargin(0.14);
  gdrspl->SetMarkerStyle(7);
  gdrspl->Draw("APL");
  gdrspl->GetXaxis()->SetTitle("Anode");
  gdrspl->GetYaxis()->SetTitle("Drift Speed (#mum/ns)");
  gdrspl->GetYaxis()->SetTitleOffset(1.35);
  c4->cd(2);
  gPad->SetLeftMargin(0.14);
  gdrspr->SetMarkerStyle(7);
  gdrspr->Draw("APL");
  gdrspr->GetXaxis()->SetTitle("Anode");
  gdrspr->GetYaxis()->SetTitle("Drift Speed (#mum/ns)");
  gdrspr->GetYaxis()->SetTitleOffset(1.35);

  return;
}

void CheckCalibInterface::DoSetlabel()
{
   // Slot method connected to the ValueSet(Long_t) signal.
   // It displays the value set in TGNumberEntry widget.
  fNumDDL=fDDL->GetNumberEntry()->GetIntNumber();
  fNumChannel=fChannel->GetNumberEntry()->GetIntNumber();
  fShowOneMod->SetCommand(Form("CheckCalibInterface::ShowSingleModule(%d,%d)",
			       fNumDDL,fNumChannel));
}

void CheckCalibInterface::ClearAll(){

  gROOT->Clear();
  if(gROOT->FindObject("c0")) delete gROOT->FindObject("c0");
  if(gROOT->FindObject("c1")) delete gROOT->FindObject("c1");
  if(gROOT->FindObject("c2")) delete gROOT->FindObject("c2");
  if(gROOT->FindObject("c3")) delete gROOT->FindObject("c3");
  if(gROOT->FindObject("c4")) delete gROOT->FindObject("c4");
  if(gROOT->FindObject("hbadchan")) delete gROOT->FindObject("hbadchan");
  if(gROOT->FindObject("hrawnoisemod")) delete gROOT->FindObject("hrawnoisemod");
  if(gROOT->FindObject("hrawnoise")) delete gROOT->FindObject("hrawnoise");
  if(gROOT->FindObject("hcornoisemod")) delete gROOT->FindObject("hcornoisemod");
  if(gROOT->FindObject("hcornoise")) delete gROOT->FindObject("hcornoise");
  if(gROOT->FindObject("hbasemod")) delete gROOT->FindObject("hbasemod");
  if(gROOT->FindObject("hbase")) delete gROOT->FindObject("hbase");
  if(gROOT->FindObject("hgainmod")) delete gROOT->FindObject("hgainmod");
  if(gROOT->FindObject("hgain")) delete gROOT->FindObject("hgainnoise");

}

void CheckCalibs()
{
   new CheckCalibInterface(gClient->GetRoot(), 400, 400); 
}
