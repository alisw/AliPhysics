#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <fstream>
#include <TFile.h>
#include <TString.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TClass.h>
#include <TDatabasePDG.h>
#include <TParameter.h>
#include <AliCounterCollection.h>
#include <AliRDHFCuts.h>

#include "MakePIDqaReport.C"

#endif

/* $Id$ */ 

TString *periodsname;

//read the file and take list and stat

Bool_t ReadFile(TList* &list,TList* &hstat, TString listname,TString partname,TString path="./",TString filename=/*"AnalysisResults.root"*/"PWG3histograms.root", TString dirname="PWG3_D2H_QA");
Bool_t ReadFileMore(TList* &list,TList* &hstat, AliRDHFCuts* &cutobj, TString listname,TString partname,TString path="./",TString filename=/*"AnalysisResults.root"*/"PWG3histograms.root", TString dirname="PWG3_D2H_QA");
void SuperimposeBBToTPCSignal(Int_t period /*0=LHC10bc, 1=LHC10d, 2=LHC10h*/,TCanvas* cpid, Int_t set);
void TPCBetheBloch(Int_t set);
Bool_t ReadFilesForCompilation(TString inputtextfile, TList**& lists, Int_t&numb, TString*& legend);
void CompilationEvSelection(Int_t n,TList** lists, TString* legend);
void CompilationTrackSelection(Int_t n,TList** lists, TString* legend);

Bool_t ReadFile(TList* &list,TList* &lstat, TString listname,TString partname,TString path,TString filename, TString dirname){

  TString lstatname="nEntriesQA", cutobjname="";
  filename.Prepend(path);
  listname+=partname;
  lstatname+=partname;

  TFile* f=new TFile(filename.Data());
  if(!f->IsOpen()){
    cout<<filename.Data()<<" not found"<<endl;
    return kFALSE;
  }
  TDirectoryFile* dir=(TDirectoryFile*)f->Get(dirname);
  if(!dir){
    cout<<dirname.Data()<<" not found in "<<filename.Data()<<endl;
    f->ls();
    return kFALSE;
  }

  list=(TList*)dir->Get(listname);
  if(!list){
    cout<<"List "<<listname.Data()<<" not found"<<endl;
    dir->ls();
    return kFALSE;
  }

  lstat=(TList*)dir->Get(lstatname);
  if(!lstat){
    cout<<lstatname.Data()<<" not found"<<endl;
    return kFALSE;
  }

  return kTRUE;
}

Bool_t ReadFileMore(TList* &list,TList* &lstat, AliRDHFCuts* &cutobj, TString listname,TString partname,TString path,TString filename,TString dirname){

  TString lstatname="nEntriesQA", cutobjname="";
  filename.Prepend(path);
  listname+=partname;
  lstatname+=partname;

  if(partname.Contains("Dplus")) cutobjname="AnalysisCuts";//"DplustoKpipiCutsStandard";
  else{
    if(partname.Contains("D0")) cutobjname="D0toKpiCutsStandard";//"D0toKpiCuts"
    else{
      if(partname.Contains("Dstar")) cutobjname="DStartoKpipiCuts";
      else{
	if(partname.Contains("Ds")) cutobjname="DstoKKpiCuts";
	else{
	  if(partname.Contains("D04")) cutobjname="D0toKpipipiCuts";
	  else{
	    if(partname.Contains("Lc")) cutobjname="LctopKpiAnalysisCuts";
	  }
	}
      }
    }
  }

  TFile* f=new TFile(filename.Data());
  if(!f->IsOpen()){
    cout<<filename.Data()<<" not found"<<endl;
    return kFALSE;
  }
  TDirectoryFile* dir=(TDirectoryFile*)f->Get(dirname);
  if(!dir){
    cout<<dirname.Data()<<" not found  in "<<filename.Data()<<endl;
    return kFALSE;
  }

  list=(TList*)dir->Get(listname);
  if(!list){
    cout<<"List "<<listname.Data()<<" not found"<<endl;
    dir->ls();
    return kFALSE;
  }

  lstat=(TList*)dir->Get(lstatname);
  if(!lstat){
    cout<<lstatname.Data()<<" not found"<<endl;
    return kFALSE;
  }

  cutobj=(AliRDHFCuts*)dir->Get(cutobjname);
  if(!cutobj){
    cout<<cutobjname.Data()<<" not found"<<endl;
    return kFALSE;
  }

  return kTRUE;
}

//draw "track related" histograms (list "outputTrack")
void DrawOutputTrack(TString partname="D00100",TString textleg="",TString path="./", Bool_t superimpose=kFALSE, TString suffixdir="",TString filename="AnalysisResults.root", Bool_t normint=kTRUE /*normalize at integral, or at nevents*/){
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetPalette(1);

  TString listname="outputTrack",name1="",name2="",path2="",filename2="PWG3histograms.root",legtext1="",legtext2="",suffixdir2="";
  TString tmp="y";
  name1=suffixdir;

  if(superimpose){
    cout<<"Enter the names:\n>";
    cin>>name1;
    cout<<">";
    cin>>name2;
    cout<<"Are they in the same output file? (y/n)"<<endl;
    cin>>tmp;
    if(tmp=="n"){
      cout<<"Path: \n";
      cout<<">";
      cin>>path2;
      cout<<"Filename: "<<endl;
      cout<<">";
      cin>>filename2;
      cout<<"Suffix in the dir name, if any (otherwhise write no)\n>";
      cin>>suffixdir2;
      if(suffixdir2=="no") suffixdir2="";
      cout<<"Text in the legend:\n 1)";
      cin>>legtext1;
      cout<<"2)";
      cin>>legtext2;
    }
    
  }

  Float_t nevents=0,nevents2=0;
  TList* list = NULL;
  TList * lstat = NULL;
  TH1F * hEntr = NULL;
  TH1F* hHasFilterBit = NULL;

  TString dirname="PWG3_D2H_QA";
  //  dirname+=suffixdir;
  Bool_t isRead=ReadFile(list,lstat,listname,Form("%s%s",partname.Data(),name1.Data()),path,filename,Form("%s%s",dirname.Data(),suffixdir.Data()));
  if(!isRead) return;
  if(!list || !lstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }
  hEntr = (TH1F*) lstat->FindObject("hNentries");
  hHasFilterBit = (TH1F*) lstat->FindObject("HasSelBit");

  nevents=hEntr->GetBinContent(1);
  TPaveText *pvtxt=new TPaveText(0.6,0.6,0.9,0.9,"NDC");
  pvtxt->SetBorderSize(0);
  pvtxt->SetFillStyle(0);
  pvtxt->AddText(legtext1);

  TList* llist = NULL;
  TList* llstat = NULL;
  TH1F* hEntr1 = NULL;
  TH1F* hHasFilterBit1 = NULL;

  if(superimpose){
    isRead=ReadFile(llist,llstat,listname,Form("%s%s",partname.Data(),name2.Data()),path2,filename2,Form("%s%s",dirname.Data(),suffixdir2.Data()));
    if(!isRead) return;
    if(!llist || !llstat){
      cout<<":-( null pointers..."<<endl;
      return;
    }
    hEntr1 = (TH1F*) llstat->FindObject("hNentries");
    hHasFilterBit1 = (TH1F*) llstat->FindObject("HasSelBit");

    nevents2=hEntr1->GetBinContent(1);
    TText *redtext=pvtxt->AddText(legtext2);
    redtext->SetTextColor(kRed);
    hEntr1->Scale(hEntr->Integral()/hEntr1->Integral());
    hHasFilterBit1->Scale(hEntr->Integral()/hEntr1->Integral());
  }

  for(Int_t i=0;i<list->GetEntries();i++){
    TH1F* h=(TH1F*)list->At(i);
    TH1F* hh=0x0;
    TH1F* hr=0x0;
    if(superimpose){
      hh=(TH1F*)llist->At(i);
      hr=(TH1F*)hh->Clone(Form("%s_ratio",hh->GetName()));
      if(hh && TMath::Abs(hh->Integral()) > 1e-6)hh->Scale(h->Integral()/hh->Integral());
    }
    if(!h || (superimpose && !hh)){
      cout<<"Histogram "<<i<<" not found"<<endl;
      continue;
    }
    if(superimpose){
      if(normint) hh->Scale(h->Integral()/hh->Integral());
      else hh->Scale(nevents/nevents2);
      hEntr1->SetLineColor(kRed);
      hh->SetLineColor(kRed);
      hr->Divide(h);
    }

    TCanvas* c=new TCanvas(Form("c%s",h->GetName()),h->GetName());
    c->cd();
    c->SetGrid();
    TString hname=h->GetName();
    if(!hname.Contains("nCls")){
      c->SetLogy();
      if(hname.Contains("Layer")){
	for(Int_t ibin=1;ibin<=h->GetNbinsX();ibin++){
	  h->GetXaxis()->SetBinLabel(ibin+1,Form("%d",ibin));
	}
	h->GetXaxis()->SetLabelSize(0.06);
	h->GetXaxis()->SetRangeUser(0,6); //comment to see ntracks!
      }
      //h->SetMinimum(1);
      if(hname.Contains("vs")||hname.Contains("dauphi")|| hname.Contains("dauphi_filt")){
	h->Draw("colz");
      }else h->Draw();

      if(superimpose) 
	{
	  hh->Draw("sames");
	  TCanvas* c2=new TCanvas(Form("c2%s",h->GetName()),h->GetName());
	  c2->cd();
	  c2->SetGrid();
	  hr->Draw();
	  c2->SaveAs(Form("%s%s%s%sRatio.png",c->GetName(),name1.Data(),name2.Data(),textleg.Data()));

	}
    } else {
      h->Draw("htext0");
      if(superimpose)hh->Draw("htext0sames");
    }
    c->cd();
    pvtxt->Draw();
    c->SaveAs(Form("%s%s%s%s.png",c->GetName(),name1.Data(),name2.Data(),textleg.Data()));
    c->SaveAs(Form("%s%s%s%s.eps",c->GetName(),name1.Data(),name2.Data(),textleg.Data()));
  }
  

  TCanvas* cstFilter=new TCanvas("cstFilter","Stat_HasFilterBit");
  cstFilter->SetGridy();
  cstFilter->cd();
  hHasFilterBit->Draw("htext0");

  if(superimpose) {
    hHasFilterBit1->Draw("htext0sames");
    pvtxt->Draw();
  }
  cstFilter->SaveAs(Form("%s%s.png",hHasFilterBit->GetName(),textleg.Data()));
  cstFilter->SaveAs(Form("%s%s.eps",hHasFilterBit->GetName(),textleg.Data()));



  TCanvas* cst=new TCanvas("cst","Stat");
  cst->SetGridy();
  cst->cd();
  hEntr->Draw("htext0");

  if(superimpose) {
    hEntr1->Draw("htext0sames");
    pvtxt->Draw();
  }
  cst->SaveAs(Form("%s%s.png",hEntr->GetName(),textleg.Data()));
  cst->SaveAs(Form("%s%s.eps",hEntr->GetName(),textleg.Data()));

  TH1F* hd0fb4=(TH1F*)list->FindObject("hd0TracksFilterBit4");
  TH1F* hd0SPD1=(TH1F*)list->FindObject("hd0TracksSPDin");
  TH1F* hd0SPDany=(TH1F*)list->FindObject("hd0TracksSPDany");
  TH1F* hd0TPCITScuts=(TH1F*)list->FindObject("hd0TracksTPCITSSPDany");
  TH1F* hhd0fb4=0x0; TH1F* hhd0SPD1=0x0; TH1F* hhd0SPDany=0x0; TH1F* hhd0TPCITScuts=0x0;
  if(superimpose){
    hhd0fb4=(TH1F*)llist->FindObject("hd0TracksFilterBit4");
    hhd0SPD1=(TH1F*)llist->FindObject("hd0TracksSPDin");
    hhd0SPDany=(TH1F*)llist->FindObject("hd0TracksSPDany");
    hhd0TPCITScuts=(TH1F*)llist->FindObject("hd0TracksTPCITSSPDany");

  }
  if(hd0fb4 && hd0SPD1 && hd0SPDany && hd0TPCITScuts){
    TCanvas* ctrsel=new TCanvas("ctrsel","Track Sel");
    ctrsel->SetLogy();
    hd0SPD1->SetLineColor(kCyan+3);
    hd0SPD1->Draw();
    ctrsel->Update();
    TPaveStats *st1=(TPaveStats*)hd0SPD1->GetListOfFunctions()->FindObject("stats");
    st1->SetTextColor(kCyan+3);
    st1->SetY1NDC(0.71);
    st1->SetY2NDC(0.9);
    hd0SPDany->SetLineColor(4);
    hd0SPDany->Draw("sames");
    ctrsel->Update();
    TPaveStats *st2=(TPaveStats*)hd0SPDany->GetListOfFunctions()->FindObject("stats");
    st2->SetY1NDC(0.51);
    st2->SetY2NDC(0.7);
    st2->SetTextColor(4);
    hd0fb4->SetLineColor(2);
    hd0fb4->Draw("sames");
    ctrsel->Update();
    TPaveStats *st3=(TPaveStats*)hd0fb4->GetListOfFunctions()->FindObject("stats");
    st3->SetY1NDC(0.31);
    st3->SetY2NDC(0.5);
    st3->SetTextColor(2);
    hd0TPCITScuts->SetLineColor(kGreen+1);
    hd0TPCITScuts->Draw("sames");
    ctrsel->Update();
    TPaveStats *st4=(TPaveStats*)hd0TPCITScuts->GetListOfFunctions()->FindObject("stats");
    st4->SetY1NDC(0.71);
    st4->SetY2NDC(0.9);
    st4->SetX1NDC(0.55);
    st4->SetX2NDC(0.75);
    st4->SetTextColor(kGreen+1);
    ctrsel->Modified();
    TLegend* leg=new TLegend(0.15,0.5,0.45,0.78);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hd0SPD1,"kITSrefit+SPD inner","L");
    leg->AddEntry(hd0SPDany,"kITSrefit+SPD any","L");
    leg->AddEntry(hd0TPCITScuts,"TPC+ITS cuts+SPD any","L");
    leg->AddEntry(hd0fb4,"Filter Bit 4","L");
    leg->Draw();
    
    if(superimpose && hhd0fb4 && hhd0SPD1 && hhd0SPDany && hhd0TPCITScuts){
      hhd0SPD1->SetLineColor(kCyan+3);
      hhd0SPD1->SetStats(0);
      hhd0SPD1->SetLineStyle(3);
      hhd0SPDany->SetLineColor(4);
      hhd0SPDany->SetStats(0);
      hhd0SPDany->SetLineStyle(3);
      hhd0TPCITScuts->SetLineColor(kGreen+1);
      hhd0TPCITScuts->SetStats(0);
      hhd0TPCITScuts->SetLineStyle(3);
      hhd0fb4->SetLineColor(2);
      hhd0fb4->SetStats(0);
      hhd0fb4->SetLineStyle(3);
      ctrsel->cd();
      hhd0SPD1->Draw("sames");
      hhd0SPDany->Draw("sames");
      hhd0TPCITScuts->Draw("sames");
      hhd0fb4->Draw("sames");

    }
    ctrsel->SaveAs("ImpactParameterTrackSel.eps");
    ctrsel->SaveAs("ImpactParameterTrackSel.png");
    
  }
}

//draw "pid related" histograms (list "outputPID")
//period=-999 to draw the pull instead of the cut
void DrawOutputPID(TString partname="D00100", Int_t mode=0/*0=with pull, 1=with nsigma*/,TString textleg="",TString path="./",TString suffixdir="", TString filename="AnalysisResults.root"){
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  Int_t period=2 ,set=0;
  if(mode==1){
    cout<<"Choose period: \n-LHC10h -> 2;\n-LHC10de -> 1;\n-LHC10bc -> 0"<<endl;
    cin>>period;
    if(period>0){
      cout<<"Choose set: "<<endl;
      if(period==2) cout<<"-pass1 -> 0;\n-pass2 -> 1"<<endl;
      cin>>set;
    }
  }

  TString listname="outputPid";
  TString dirname="PWG3_D2H_QA";
  dirname+=suffixdir;

  TList* list = NULL;
  TList * lstat = NULL;
  //needed only for mode 1
  AliRDHFCuts* cutobj;
  AliAODPidHF* aodpid;
  Double_t nsigmaTOF=0;
  Double_t nsigmaTPC[3]={},plimTPC[2]={};

  if(mode==1){
    Bool_t isRead=ReadFileMore(list,lstat,cutobj,listname,partname+suffixdir,path,filename,dirname);
    if(!isRead) return;
    if(!list || !lstat){
      cout<<":-( null pointers..."<<endl;
      return;
    }
    aodpid=(AliAODPidHF*)cutobj->GetPidHF();
    if(!aodpid){
      cout<<"PidHF object not found! cannot get the nsigma values"<<endl;
      return;
    }
    nsigmaTOF=aodpid->GetSigma(3);
  
    nsigmaTPC[0]=aodpid->GetSigma(0);
    nsigmaTPC[1]=aodpid->GetSigma(1);
    nsigmaTPC[2]=aodpid->GetSigma(2);
    aodpid->GetPLimit(plimTPC);

  }else{
    Bool_t isRead=ReadFile(list,lstat,listname,partname+suffixdir,path,filename,dirname);
    if(!isRead) return;
    if(!list || !lstat){
      cout<<":-( null pointers..."<<endl;
      return;
    }
  }


  TPaveText *txtsigmaTOF=new TPaveText(0.1,0.65,0.5,0.9,"NDC");
  txtsigmaTOF->SetBorderSize(0);
  txtsigmaTOF->SetFillStyle(0);
  txtsigmaTOF->AddText(Form("nsigmacut from cutobj = %.1f",nsigmaTOF));
  TLine lTOF;
  lTOF.SetLineColor(kMagenta+1);
  lTOF.SetLineStyle(2);
  lTOF.SetLineWidth(3);

  TPaveText *txtsigmaTPC=new TPaveText(0.3,0.6,0.6,0.9,"NDC");
  txtsigmaTPC->SetBorderSize(0);
  txtsigmaTPC->SetFillStyle(0);
  txtsigmaTPC->AddText("nsigmacut from cutobj \n");
  txtsigmaTPC->AddText(Form("p < %.1f : %.1f \n",plimTPC[0],nsigmaTPC[0]));
  txtsigmaTPC->AddText(Form("%.1f < p < %.1f : %.1f \n",plimTPC[0],plimTPC[1],nsigmaTPC[1]));
  txtsigmaTPC->AddText(Form("p > %.1f : %.1f \n",plimTPC[1],nsigmaTPC[2]));
  TLine lTPC;
  lTPC.SetLineColor(kMagenta+1);
  lTPC.SetLineStyle(2);
  lTPC.SetLineWidth(3);

  // TCanvas *ctest=new TCanvas("text","Test text");
  // ctest->cd();
  // txtsigmaTPC->Draw();
  // txtsigmaTOF->Draw();


  for(Int_t i=0;i<list->GetEntries();i++){
    TClass* objtype=list->At(i)->IsA();
    TString tpname=objtype->GetName();

    if(tpname=="TH1F"){
      TH1F* h=(TH1F*)list->At(i);

      if(!h){
	cout<<"Histogram "<<i<<" not found"<<endl;
	continue;
      }
      //h->Scale(1./h->Integral("width"));
      TCanvas* c=new TCanvas(Form("c%s",h->GetName()),h->GetName());
      c->SetLogz();
      c->cd();
      h->Draw();
    
      //write
      c->SaveAs(Form("%s%s.png",h->GetName(),textleg.Data()));
      c->SaveAs(Form("%s%s.eps",h->GetName(),textleg.Data()));
      TFile* fout=new TFile(Form("%s.root",h->GetName()),"recreate");
      fout->cd();
      c->Write();
    }
  
    if(tpname=="TH2F"){
      TH2F* h=(TH2F*)list->At(i);
      
      if(!h){
	cout<<"Histogram "<<i<<" not found"<<endl;
	continue;
      }
      TString hname=h->GetName();
      h->Sumw2();
      if(h->Integral("width")==0) {cout<<"Empty histogram, skip\n"; continue;}
      h->Scale(1./h->Integral("width"));

      Double_t maxzaxis=h->GetBinContent(h->GetMaximumBin());
      Double_t minzaxis=h->GetBinContent(h->GetMinimumBin());
      printf("Minimum = %f, maximum = %f\n",minzaxis,maxzaxis);
      TH2F* hallzrange=(TH2F*)h->Clone(Form("%swholez",hname.Data()));
      hallzrange->SetAxisRange(1e-07,maxzaxis,"Z");
      //hallzrange->SetAxisRange(minzaxis,maxzaxis,"Z");

      TCanvas* cwholez=new TCanvas(Form("c%swholez",hname.Data()),Form("%s down to lowest z",hname.Data()));
      cwholez->SetLogz();
      hallzrange->Draw("colz");
      cwholez->SaveAs(Form("%swholez.png",h->GetName()));
      cwholez->SaveAs(Form("%swholez.eps",h->GetName()));

      if(hname.Contains("hTOFtimeKaonHyptime")){
	TCanvas* cz=new TCanvas(Form("c%szoom",hname.Data()),Form("%szoom",hname.Data()));
	cz->SetLogz();
	TH2F* hz=(TH2F*)h->Clone(Form("%sz",hname.Data()));
	hz->Draw("colz");
	hz->SetAxisRange(-1500,1500,"Y");
	hz->SetAxisRange(0.,5.,"X");
	//write
	cz->SaveAs(Form("%szoom.png",h->GetName()));
	cz->SaveAs(Form("%szoom.eps",h->GetName()));
      }

      TCanvas* c=new TCanvas(Form("c%s",hname.Data()),hname.Data());
      c->SetLogz();
      //c->SetLogx();
      TCanvas* c2=new TCanvas(Form("c2%s",hname.Data()),hname.Data());
      c2->SetLogz();

      c->cd();
      h->DrawClone("colz");

      if (hname.Contains("Sig") || hname.Contains("sigma"))h->SetAxisRange(-5,5,"Y");
      c2->cd();
      //if (hname.Contains("TOFtime"))h->SetAxisRange(-1500,1500,"Y");
      h->SetAxisRange(0.,5.,"X");
     
      h->Draw("colz");
     
      //TCanvas *test=new TCanvas("test","test");
      if(mode==0){
      	 AddFit(h);	
      }else{ //mode 1

	if(hname.Contains("TOFsigma")) {

	  c->cd();
	  txtsigmaTOF->Draw();
	  lTOF.DrawLine(.2,nsigmaTOF,20,nsigmaTOF);
	  lTOF.DrawLine(.2,-1*nsigmaTOF,4.,-1*nsigmaTOF);

	}
      

	if(hname.Contains("TPCsigma")){

	  c->cd();
	  txtsigmaTPC->Draw();
	  lTPC.DrawLine(0.,nsigmaTPC[0],plimTPC[0],nsigmaTPC[0]);
	  lTPC.DrawLine(plimTPC[0],nsigmaTPC[1],plimTPC[1],nsigmaTPC[1]);
	  lTPC.DrawLine(plimTPC[1],nsigmaTPC[2],4,nsigmaTPC[2]);
	  lTPC.DrawLine(0.,-1*nsigmaTPC[0],plimTPC[0],-1*nsigmaTPC[0]);
	  lTPC.DrawLine(plimTPC[0],-1*nsigmaTPC[1],plimTPC[1],-1*nsigmaTPC[1]);
	  lTPC.DrawLine(plimTPC[1],-1*nsigmaTPC[2],4,-1*nsigmaTPC[2]);
	}

	if(hname.Contains("TPCsigvsp")){
	  SuperimposeBBToTPCSignal(period,c,set);
	}
      }
	
      //write
      c->SaveAs(Form("%s%d.png",h->GetName(),mode));
      c->SaveAs(Form("%s%d.eps",h->GetName(),mode));
      c2->SaveAs(Form("%s2%d.png",h->GetName(),mode));
      c2->SaveAs(Form("%s2%d.eps",h->GetName(),mode));

      TFile* fout=new TFile(Form("%s%d.root",h->GetName(),mode),"recreate");
      fout->cd();
      c->Write();
      c2->Write();
    }
  }
}

void SuperimposeBBToTPCSignal(Int_t period /*0=LHC10bc, 1=LHC10d, 2=LHC10h*/,TCanvas* cpid,Int_t set /*see below*/){

  TFile* fBethe=new TFile("BetheBlochTPC.root");
  if(!fBethe->IsOpen()){
    TPCBetheBloch(set);
    fBethe=new TFile("BetheBlochTPC.root");
  }
  const Int_t npart=4;
  TString partnames[npart]={"Kaon","Pion","Electron","Proton"};
  for(Int_t ipart=0;ipart<npart;ipart++){
    TString grname=Form("%sP%d",partnames[ipart].Data(),period);
    TGraph* gr=(TGraph*)fBethe->Get(grname);
    cpid->cd();
    gr->SetLineColor(1);
    gr->SetLineWidth(2);
    gr->Draw("L");
  }

  //cpid->SaveAs(Form("%sBB.png",hname.Data()));
}

//draw and save Bethe Bloch from TPC in different periods
void TPCBetheBloch(Int_t set){
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(0);

  AliTPCPIDResponse *tpcResp=new AliTPCPIDResponse();

  const Int_t npart=4;
  //Double_t masses[npart]={TDatabasePDG::Instance()->GetParticle(321)->Mass()/*Kaon*/,TDatabasePDG::Instance()->GetParticle(211)->Mass()/*Pion*/,TDatabasePDG::Instance()->GetParticle(11)->Mass()/*Electron*/,TDatabasePDG::Instance()->GetParticle(2212)->Mass()/*Proton*/};
  TString partnames[npart]={"Kaon","Pion","Electron","Proton"};
  //printf("%s = %.4f,%s = %.4f,%s = %.4f\n",partnames[0].Data(),masses[0],partnames[1].Data(),masses[1],partnames[2].Data(),masses[2]);
  TCanvas *cBethe=new TCanvas("cBethe","Bethe Bloch K pi e p");
  Int_t nperiods=4; //LHC10b+c, LHC10d, LHC10h, MC
  Double_t alephParameters[5]={};
  Int_t nsets=1/*LHC10bc*/+2/*LHC10de*/+2/*LHC10h*/+3/*MC*/;

  periodsname=new TString[nsets];
  cout<<"Creating the file of the Bethe Bloch"<<endl;
  TFile* fout=new TFile("BetheBlochTPC.root","recreate");

  for(Int_t iperiod=0;iperiod<nperiods;iperiod++){
    cout<<"Period "<<iperiod<<" : ";
    if(iperiod==0){ //LHC10bc
      
      alephParameters[0] = 0.0283086/0.97;
      alephParameters[1] = 2.63394e+01;
      alephParameters[2] = 5.04114e-11;
      alephParameters[3] = 2.12543e+00;
      alephParameters[4] = 4.88663e+00;
      periodsname[0]="dataLHC10bc";  
    }
    if(iperiod==1){ //LHC10de,low energy
      if(set==0){   
	alephParameters[0] = 1.63246/50.;
	alephParameters[1] = 2.20028e+01;
	alephParameters[2] = TMath::Exp(-2.48879e+01);
	alephParameters[3] = 2.39804e+00;
	alephParameters[4] = 5.12090e+00;
	periodsname[1]="dataLHC10deold"; 
      }
      if(set==1){
	alephParameters[0] = 1.34490e+00/50.;
	alephParameters[1] =  2.69455e+01;
	alephParameters[2] =  TMath::Exp(-2.97552e+01);
	alephParameters[3] = 2.35339e+00;
	alephParameters[4] = 5.98079e+00;
	periodsname[2]="dataLHC10denew";
      }
    }

    if(iperiod==2){//LHC10h
      if(set==0){//pass1 
	alephParameters[0]=1.25202/50.;
	alephParameters[1]=2.74992e+01;
	alephParameters[2]=TMath::Exp(-3.31517e+01);
	alephParameters[3]=2.46246;
	alephParameters[4]=6.78938;
	periodsname[3]="dataLHC10hpass1";
      }
      if (set==1){//pass2 (AOD044)
	alephParameters[0]=1.25202/50.;
	alephParameters[1]=2.74992e+01;
	alephParameters[2]=TMath::Exp(-3.31517e+01);
	alephParameters[3]=2.46246;
	alephParameters[4]=6.78938;
	periodsname[4]="dataLHC10hpass2";
      }
    }
    if(iperiod==3){ //MC
      if(set==0){
	alephParameters[0] = 2.15898e+00/50.;
	alephParameters[1] = 1.75295e+01;
	alephParameters[2] = 3.40030e-09;
	alephParameters[3] = 1.96178e+00;
	alephParameters[4] = 3.91720e+00;
	periodsname[5]="MCold";
      }
      if(set==1){ //new
	alephParameters[0] = 1.44405/50;
	alephParameters[1] = 2.35409e+01;
	alephParameters[2] = TMath::Exp(-2.90330e+01);
	alephParameters[3] = 2.10681;
	alephParameters[4] = 4.62254;
	periodsname[6]="MCnew";
      }

      if(set==2){ //new BB from Francesco
	alephParameters[0] = 0.029021;
	alephParameters[1] = 25.4181;
	alephParameters[2] = 4.66596e-08;
	alephParameters[3] = 1.90008;
	alephParameters[4] = 4.63783;
	periodsname[7]="MCBBFrancesco";
      }

      if(set==3){ //low energy 2011
	alephParameters[0] = 0.0207667;
	alephParameters[1] = 29.9936;
	alephParameters[2] = 3.87866e-11;
	alephParameters[3] = 2.17291;
	alephParameters[4] = 7.1623;
	//periodsname[8]="MClowen2011";
      }


    }
    //cout<<periodsname[iperiod]<<endl;
    tpcResp->SetBetheBlochParameters(alephParameters[0],alephParameters[1],alephParameters[2],alephParameters[3],alephParameters[4]);
    cout<<"here"<<endl;
    for(Int_t ipart=0;ipart<npart;ipart++){

      const Int_t n=1000;
      Double_t p[n],bethe[n];

      for(Int_t k=0;k<n;k++){ //loop on the momentum steps
	p[k]=0.0001+k*4./n; //limits 0.-4. GeV/c
	//cout<<p[k]<<"\t";
	//bethe[k]=-tpcResp->Bethe(p[k]/masses[ipart]);
	AliPID::EParticleType ptype=AliPID::kKaon;
	if(ipart==1) ptype=AliPID::kPion;
	if(ipart==2) ptype=AliPID::kElectron;
	if(ipart==3) ptype=AliPID::kProton;
	bethe[k]=tpcResp->GetExpectedSignal(p[k],ptype);
      }
      //cout<<endl;
      TGraph *gr=new TGraph(n,p,bethe);
      gr->SetName(Form("%sP%d",partnames[ipart].Data(),iperiod));
      gr->SetTitle(Form("%sP%d;p (GeV/c);",partnames[ipart].Data(),iperiod));
      gr->SetLineColor(ipart+1);
      gr->SetMarkerColor(ipart+1);
      gr->GetYaxis()->SetRangeUser(35,100);
      cBethe->cd();
      if(iperiod==0 && ipart==0)gr->DrawClone("AL");
      else gr->DrawClone("L");

      fout->cd();
      gr->Write();
    }

  }
  TParameter<int> sett;
  sett.SetVal(set);
  fout->cd();
  sett.Write();

  fout->Close();
}

void DrawOutputCentrality(TString partname="D0",TString textleg="",TString path="./", Bool_t superimpose=kFALSE,TString suffixdir="",TString filename=/*"AnalysisResults.root"*/"PWG3histograms.root"){
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetPalette(1);

  TString listname="outputCentrCheck",partname2="",path2="",suffixdir2="",filename2="PWG3histograms.root";

  // if(superimpose){
  //   cout<<"Enter the names:\n>";
  //   cin>>name1;
  //   cout<<">";
  //   cin>>name2;
  // }
  // TString listname="outputTrack",name1="",name2="";
  TString tmp="y";

  if(superimpose){
    cout<<"##Second file\n";
    cout<<"Enter the name:\n";
    cout<<">";
    cin>>partname2;
    cout<<"Are they in the same output file? (y/n)"<<endl;
    cin>>tmp;
    if(tmp=="n"){
      cout<<"Path: \n";
      cout<<">";
      cin>>path2;
      cout<<"Dir name:\n";
      cout<<">";
      cin>>suffixdir2;
      cout<<"Filename: "<<endl;
      cout<<">";
      cin>>filename2;
    }
    
  }
  // Int_t nhist=1;
  // TString *name=0x0;
  // if(superimpose){
  //   cout<<"Number of histogram to superimpose: ";
  //   cin>>nhist;
  //   name=new TString[nhist];
  //   for (Int_t j=0;j<nhist;j++){
  //     cout<<">";
  //     cin>>name[j];
  //   }
  // }

  TList* list = NULL;
  TList * lstat = NULL;
  TH1F* hEntr = NULL;

  TString dirname="PWG3_D2H_QA",dirname2=dirname;
  dirname+=suffixdir;
  dirname2+=suffixdir2;
  Bool_t isRead=ReadFile(list,lstat,listname,partname.Data(),path,filename,dirname);
  if(!isRead) return;
  if(!list || !lstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }

  TPaveText *pvtxt=new TPaveText(0.6,0.6,0.9,0.9,"NDC");
  pvtxt->SetBorderSize(0);
  pvtxt->SetFillStyle(0);
  pvtxt->AddText(partname);

  TList* llist = NULL;
  TList* llstat = NULL;
  TH1F* hEntr1 = NULL;

  if(superimpose){
    isRead=ReadFile(llist,llstat,listname,partname2.Data(),path2,filename2,dirname2);
    if(!isRead) return;
    if(!llist || !llstat){
      cout<<":-( null pointers..."<<endl;
      return;
    }
    TText *redtext=pvtxt->AddText(partname2);
    redtext->SetTextColor(kRed);

  }


  TCanvas* cst=new TCanvas("cst","Stat");
  cst->SetGridy();
  cst->cd();
  hEntr=(TH1F*)lstat->FindObject("hNentries");
  hEntr1=(TH1F*)llstat->FindObject("hNentries");
  Int_t nevents=hEntr->GetBinContent(1);
  hEntr->Draw("htext0");
  cst->SaveAs(Form("%s%s.png",hEntr->GetName(),textleg.Data()));
  cst->SaveAs(Form("%s%s.eps",hEntr->GetName(),textleg.Data()));
  Int_t nevents080=1,nnevents080=1;

  //TCanvas *spare=new TCanvas("sparecv","Spare");

  for(Int_t i=0;i<list->GetEntries();i++){

    TClass* objtype=list->At(i)->IsA();
    TString tpname=objtype->GetName();

    if(tpname=="TH1F"){

      TH1F* h=(TH1F*)list->At(i);
      TH1F* hh=0x0;
      if(superimpose){
	hh=(TH1F*)llist->At(i);
      }
      if(!h || (superimpose && !hh)){
	cout<<"Histogram "<<i<<" not found"<<endl;
	continue;
      }
      if(superimpose){
	hEntr1->SetLineColor(kRed);
	hh->SetLineColor(kRed);
      }

      TCanvas* c=new TCanvas(Form("c%s",h->GetName()),h->GetName());
      TPaveText *pvtxt2=new TPaveText(0.6,0.6,0.9,0.9,"NDC");
      pvtxt2->SetBorderSize(0);
      pvtxt2->SetFillStyle(0);

      c->cd();
      c->SetGrid();
      c->SetLogy();
      Int_t entries=h->Integral();
      pvtxt2->AddText(Form("%.1f %s of the events",(Double_t)entries/(Double_t)nevents*100,"%"));
      h->Draw();
      if(superimpose) {
	hh->Draw("sames");
	pvtxt->Draw();
      }
      pvtxt2->Draw();
      c->SaveAs(Form("%s%s.pdf",c->GetName(),textleg.Data()));
      c->SaveAs(Form("%s%s.eps",c->GetName(),textleg.Data()));
    }
    if(tpname=="TH2F"){
      TH2F* hhh=(TH2F*)list->At(i);
      if(!hhh){
	cout<<"Histogram "<<i<<" not found"<<endl;
	continue;
      }
      TCanvas* c=new TCanvas(Form("c%s",hhh->GetName()),hhh->GetName());
      TPaveText *pvtxt3=new TPaveText(0.6,0.6,0.9,0.9,"NDC");
      pvtxt3->SetBorderSize(0);
      pvtxt3->SetFillStyle(0);

      c->cd();
      c->SetGrid();
      Int_t entries=hhh->Integral();
      pvtxt3->AddText(Form("%.1f %s of the events",(Double_t)entries/(Double_t)nevents*100,"%"));
      hhh->Draw("colz");
      c->SetLogz();
      pvtxt3->Draw();
      c->SaveAs(Form("%s%s.pdf",c->GetName(),textleg.Data()));
      c->SaveAs(Form("%s%s.eps",c->GetName(),textleg.Data()));
    }
  }
  
  
  listname="countersCentrality";

  isRead=ReadFile(list,lstat,listname,partname.Data(),path,filename,dirname);
  if(!isRead) return;
  if(!list || !lstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }


  if(superimpose){
    isRead=ReadFile(llist,llstat,listname,partname2.Data(),path2,filename2,dirname2);
    if(!isRead) return;
    if(!llist || !llstat){
      cout<<":-( null pointers..."<<endl;
      return;
    }
    TText *redtext=pvtxt->AddText(partname2);
    redtext->SetTextColor(kRed);

  }

  TH1F* hallcntr=0x0;
  TH1F* hhallcntr=0x0;
  cout<<"normalizing to 0-80% as a check"<<endl;
  Int_t ncentr=10;//check this
  TH1F* h020=0x0;
  TH1F* h2080=0x0;
  TH1F* hh020=0x0;
  TH1F* hh2080=0x0;

  TCanvas *cvnocnt=new TCanvas("cvnocnt","No Centrality estimation",800,400);
  cvnocnt->Divide(2,1);
  TCanvas *ccent=0x0;

  for(Int_t i=0;i<list->GetEntries();i++){
    AliCounterCollection* coll=(AliCounterCollection*)list->At(i);
    AliCounterCollection* colle=0x0;
    if(superimpose) colle=(AliCounterCollection*)llist->At(i);
    coll->SortRubric("run");//sort by run number

    h020=0x0;
    h2080=0x0;
    hh020=0x0;
    hh2080=0x0;
    
    hallcntr=0x0; 
    hhallcntr=0x0; 

    TH1F* hbad=(TH1F*)coll->Get("run",Form("centralityclass:-990_-980"));
    cvnocnt->cd(i+1);
    if(hbad) hbad->Draw();

    ccent=new TCanvas(Form("ccent%s",coll->GetName()),Form("Centrality vs Run (%s)",coll->GetName()),1400,800);
    ccent->SetTicky();
    ccent->Divide(4,2);
    
    TH1F* hh=0x0;

    for(Int_t ic=0;ic<8/*ncentr*/;ic++){ //normalizing to 0-80% as a check

      TH1F* h=(TH1F*)coll->Get("run",Form("centralityclass:%d_%d",ic*10,ic*10+10));
      h->SetName(Form("h%d%d",i,ic));
      if(!hallcntr) {
	hallcntr=(TH1F*)h->Clone("hallcntr");
	hallcntr->Sumw2();
      } else {
	hallcntr->Add(h);
      }
      
      nevents080+=h->Integral();

      if(superimpose){
	hh=(TH1F*)colle->Get("run",Form("centralityclass:%d_%d",ic*10,ic*10+10));
	hh->SetName(Form("hh%d%d",i,ic));
	if(!hhallcntr) {
	  hhallcntr=(TH1F*)hh->Clone("hhallcntr");
	  hhallcntr->Sumw2();
	}else hhallcntr->Add(hh);

	nnevents080+=hh->Integral();
	
      }
    }

    for(Int_t ic=0;ic<ncentr;ic++){

      TH1F* h=(TH1F*)coll->Get("run",Form("centralityclass:%d_%d",ic*10,ic*10+10));
      h->SetName(Form("h%d%d",i,ic));
      h->Sumw2();
      
      if(ic>=0 && ic<=1){ //0-20
	if(!h020) {
	  h020=(TH1F*)h->Clone(Form("h020%s",coll->GetName()));
	  h020->SetTitle(Form("Centrality 0-20 %s",coll->GetName()));
	  if(superimpose){
	    hh020=(TH1F*)hh->Clone(Form("hh020%s",coll->GetName()));
	    hh020->SetTitle(Form("Centrality 0-20 %s",coll->GetName()));
	  }
	}
	else {
	  h020->Add(h);
	  if(superimpose)hh020->Add(hh);
	}
      }
      if(ic>=2 && ic<=7){ //20-80
	if(!h2080) {
	  h2080=(TH1F*)h->Clone(Form("h2080%s",coll->GetName()));
	  h2080->SetTitle(Form("Centrality 20-80 %s",coll->GetName()));
	  if(superimpose){
	    hh2080=(TH1F*)hh->Clone(Form("hh2080%s",coll->GetName()));
	    hh2080->SetTitle(Form("Centrality 20-80 %s",coll->GetName()));
	  }
	}
	else {
	  h2080->Add(h);
	  if(superimpose)hh2080->Add(hh);
	}
	
      }

      h->Divide(hallcntr);

      if(ic<8){
	ccent->cd(ic+1);
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetTitleOffset(1.5);
	h->SetMinimum(0);
	//h->GetYaxis()->SetRangeUser(0.,0.15);
	h->DrawClone();
      }
      /*
	if(ic==0&&i==0){
	spare->cd();
	h->Draw();
	}
      */
      // ccent->cd(1);
      // h->SetLineColor(ic+1);
      // if(ic==0)h->DrawClone();
      // else h->DrawClone("sames");
    }
    h020->Divide(hallcntr);
    if(superimpose){
      hh020->Divide(hhallcntr);
      hh020->SetLineColor(2);
      hh020->SetMarkerColor(2);
    }

    /*//draw 0-20 and 20-80 in the multi pad canvas (increase divisions before uncommenting)
    ccent->cd(ncentr+1);
    h020->DrawClone();
    if(superimpose){
      hh020->DrawClone("sames");
    }
    */
    TCanvas* cv020=new TCanvas(Form("cv020-%d",i),"0-20% vs run number",1400,600);
    cv020->cd();
    cv020->SetTicky();
    h020->GetYaxis()->SetRangeUser(0.,1.);
    h020->DrawClone();
    if(superimpose)hh020->DrawClone("sames");
    cv020->SaveAs(Form("cv020-%d.pdf",i));
    cv020->SaveAs(Form("cv020-%d.eps",i));

    h2080->Divide(hallcntr);
    if(superimpose) {
      hh2080->Divide(hhallcntr);
      hh2080->SetLineColor(2);
      hh2080->SetMarkerColor(2);
    }

    /*
    ccent->cd(ncentr+2);
    h2080->DrawClone();
   
    if(superimpose){
      hh2080->DrawClone("sames");
    }
    */
    TCanvas* cv2080=new TCanvas(Form("cv2080-%d",i),"20-80% vs run number",1400,600);
    cv2080->cd();
    cv2080->SetTicky();
    h2080->GetYaxis()->SetRangeUser(0.,1.);
    h2080->DrawClone();
    if(superimpose)hh2080->DrawClone("sames");
    cv2080->SaveAs(Form("cv2080-%d.pdf",i));
    cv2080->SaveAs(Form("cv2080-%d.eps",i));

    ccent->SaveAs(Form("%s%s.pdf",ccent->GetName(),textleg.Data()));
    ccent->SaveAs(Form("%s%s.eps",ccent->GetName(),textleg.Data()));
  }
  
}

void DrawProjections(TString partname="D0",TString h2dname="hMultvsPercentile",Int_t groupnbins=5,Float_t fitmin=15,Float_t fitmax=50,TString direction="X",TString path="./",TString suffixdir="", TString filename="AnalysisResults.root", TString fitfunc="pol0"/*option "nofit" does not fit*/){
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetPalette(1);

  TString listname="outputCentrCheck";
  TString dirname="PWG3_D2H_QA";
  dirname+=suffixdir;

  TList* list = NULL;
  TList * lstat = NULL;
  TH1F* hEntr = NULL;

  Bool_t isRead=ReadFile(list,lstat,listname,partname,path,filename,dirname);
  if(!isRead) return;
  if(!list || !lstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }
  hEntr = (TH1F*) lstat->FindObject("hNentries");
  Double_t nevents=hEntr->GetBinContent(5); //ev good vertex

  TH2F* h2=(TH2F*)list->FindObject(h2dname);
  if(!h2){
    cout<<h2dname.Data()<<" not found"<<endl;
    list->ls();
    return;
  }
  TCanvas* cv2d=new TCanvas("cv2d",h2->GetName());
  cv2d->cd();
  cv2d->SetLogz();
  cv2d->SetGrid();
  h2->Draw("colz");
  TPaveText *pvst=new TPaveText(0.6,0.2,0.9,0.7,"NDC");
  pvst->SetBorderSize(0);
  pvst->SetFillStyle(0);
  pvst->AddText("Bin -> Cont/nEvVtx");

  //nsteps=group bins in the Y(X) direction if projecting on the X(Y) direction
  Int_t nsteps=0;

  if(direction=="X") nsteps=h2->GetNbinsY()/groupnbins;
  if(direction=="Y") nsteps=h2->GetNbinsX()/groupnbins;
  cout<<"Grouping bins by " <<groupnbins<<" I obtaine "<<nsteps<<" projections"<<endl;

  TCanvas *cvpj=new TCanvas(Form("cvpj%s%s",direction.Data(),h2dname.Data()),Form("cvpj%s",direction.Data()),1200,800);
  cvpj->Divide((Int_t)(nsteps/3)+1,3);
  TFile* fout=new TFile(Form("proj%s%s.root",direction.Data(),h2dname.Data()), "recreate");
  //Float_t maxx[nsteps];
  //Float_t maxx[12]={9000,9000,6000,4000,2000,1400,800,500,200,100,40,25};
  Double_t integralpernev[nsteps];

  Double_t minx=0,maxx=0;
  if(direction=="X"){
    minx=h2->GetYaxis()->GetXmin();
    maxx=h2->GetYaxis()->GetXmax();
  }
  if(direction=="Y"){
    minx=h2->GetXaxis()->GetXmin();
    maxx=h2->GetXaxis()->GetXmax();
  }
  printf("Plotting from %.1f to %.1f\n",minx,maxx);
  TCanvas *cintegral=new TCanvas("cintegral","Integral of each projection");
  TH1F* hint=new TH1F("hint","Integral of each projection;Centrality (%);Entries",nsteps,minx,maxx);
  Double_t minint=999999999,maxint=0;

  for(Int_t i=0;i<nsteps;i++){
    TH1F* h=0x0;
    // if(direction=="X")h=(TH1F*)h2->ProjectionX(Form("px%d",i),i+kbins,i+2*kbins);
    // if(direction=="Y")h=(TH1F*)h2->ProjectionY(Form("py%d",i),i+kbins,i+2*kbins);
    if(direction=="X")h=(TH1F*)h2->ProjectionX(Form("px%d",i),groupnbins*i+1,groupnbins*(i+1));
    if(direction=="Y")h=(TH1F*)h2->ProjectionY(Form("py%d",i),groupnbins*i+1,groupnbins*(i+1));
    Double_t projint=h->Integral();
    cout<<"Integral of projection "<<i<<" = "<<projint<<endl;
    hint->SetBinContent(i+1,projint);
    hint->SetBinError(i+1,TMath::Sqrt(projint));

    if(projint<1e-7) continue;
    if(minint>projint) minint=projint;
    if(projint>maxint) maxint=projint;
    integralpernev[i]=h->Integral()/nevents;

    TPaveText *pvtxt=new TPaveText(0.6,0.6,0.9,0.9,"NDC");
    pvtxt->SetBorderSize(0);
    pvtxt->SetFillStyle(0);
    pvtxt->AddText(Form("%.0f - %.0f",h2->GetYaxis()->GetBinLowEdge(groupnbins*i+1),h2->GetYaxis()->GetBinLowEdge(groupnbins*(i+1))));
    pvst->AddText(Form("%.0f - %.0f -> %.2f",h2->GetYaxis()->GetBinLowEdge(groupnbins*i+1),h2->GetYaxis()->GetBinLowEdge((groupnbins*(i+1))),integralpernev[i]));

    cvpj->cd(i+1);
    //h->GetXaxis()->SetRangeUser(0,maxx[i]);
    h->Draw();
    pvtxt->Draw();
    fout->cd();
    h->Write();
  }
  cvpj->SaveAs(Form("cvpj%s%s.pdf",direction.Data(),h2dname.Data()));
  cvpj->SaveAs(Form("cvpj%s%s.eps",direction.Data(),h2dname.Data()));

  cv2d->cd();
  pvst->Draw();
  cv2d->SaveAs(Form("%s.pdf",h2->GetName()));
  cv2d->SaveAs(Form("%s.eps",h2->GetName()));

  cintegral->cd();
  hint->SetMarkerStyle(20);
  hint->Draw("PE");
  if(!fitfunc.Contains("nofit")){
    hint->Fit(fitfunc.Data(),"RL","PE",fitmin,fitmax);
    TF1* fpolfit=hint->GetFunction(fitfunc.Data());
    TPaveText *txtvar=new TPaveText(0.3,0.1,0.9,0.4,"NDC");
    txtvar->SetBorderSize(0);
    txtvar->SetFillStyle(0);
    //txtvar->AddText(Form("Full spread %.0f- %.0f",maxint,minint));
    txtvar->AddText(Form("Fit in %.1f-%.1f; ",fitmin,fitmax));
    for(Int_t ipar=0;ipar<fpolfit->GetNpar();ipar++){
      txtvar->AddText(Form("par%d = %.0f, ",ipar, fpolfit->GetParameter(ipar)));
    }
    txtvar->AddText(Form("#tilde{#chi}^{2} = %.2f",fpolfit->GetChisquare()/fpolfit->GetNDF()));
    txtvar->AddText(Form("bin width = %.1f %s",hint->GetBinWidth(3),"%"));
    txtvar->Draw();
  }
  fout->cd();
  hint->Write();
  cintegral->SaveAs(Form("%s.pdf",hint->GetName()));
  cintegral->SaveAs(Form("%s.eps",hint->GetName()));
}

void DrawEventSelection(TString partname="D0", TString path="./",TString suffixdir="",TString filename="AnalysisResults.root"){
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  TString listname="outputEvSel";
  TString dirname="PWG3_D2H_QA";
  dirname+=suffixdir;

  TList* list= NULL;
  TList * lstat = NULL;
  TH1F* hEntr= NULL;
  

  Bool_t isRead=ReadFile(list,lstat,listname,partname,path,filename,dirname);
  if(!isRead) return;
  if(!list || !lstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }
  hEntr =(TH1F *)lstat->FindObject("hNentries");
  //Double_t neventsgv=hstat->Integral(5,5); //ev good vertex

  for(Int_t i=0;i<list->GetEntries();i++){

    TClass* objtype=list->At(i)->IsA();
    TString tpname=objtype->GetName();

    if(tpname=="TH1F"){
      TH1F* htmp=(TH1F*)list->At(i);
      TCanvas* c=new TCanvas(Form("c%s",htmp->GetName()),Form("c%s",htmp->GetName()));
      c->cd();
      htmp->Draw();
      c->SaveAs(Form("%s.pdf",htmp->GetName()));
      c->SaveAs(Form("%s.eps",htmp->GetName()));
    }

    if(tpname=="TH2F"){
      TH2F* htmp=(TH2F*)list->At(i);
      TCanvas* c=new TCanvas(Form("c%s",htmp->GetName()),Form("c%s",htmp->GetName()),1200,800);
      c->cd();
      htmp->SetMarkerSize(1.3);
      htmp->Draw("colzhtext45");
      c->SaveAs(Form("%s.pdf",htmp->GetName()));
      c->SaveAs(Form("%s.eps",htmp->GetName()));
    }
  }

  AliCounterCollection* coll=(AliCounterCollection*)list->FindObject("trigCounter");
  if(!coll) {
    cout<<"Trigger counter not found"<<endl;
    return;
  }
  
  coll->SortRubric("run");//sort by run number

  TString collname=coll->GetName();

  TString keywords=coll->GetKeyWords("triggertype");

  Int_t nkeys=keywords.CountChar(',')+1;

  TString *words = new TString[nkeys];
  for(Int_t k=0;k<nkeys;k++) words[k]="";
  printf("Keywords: ");
  Int_t count=0;
  for(Int_t l=0;l<keywords.Length();l++){
    if(keywords[l] != ',') words[count]+=keywords[l];
    else {
      printf("%s ",words[count].Data());
      count++;
    }
  }
  cout<<endl;

  TH1D** htrig=new TH1D*[nkeys]; //each trigger type in one histogram of counts vs run
  TH1D** htrignormAll=new TH1D*[nkeys]; //normalized to the counts in kAll
  TH1D** htrignormAny=new TH1D*[nkeys]; //normalized to the counts in kAny
  TCanvas* ctrigfractionAll=new TCanvas("cvtrigfrac","Fraction of given trigger type vs run (kAll)",1400,800);
  TCanvas* ctrigfractionAny=new TCanvas("cvtrigfracAll","Fraction of given trigger type vs run (kAny)",1400,800);
  TLegend* legtr=new TLegend(0.15,0.5,0.35,0.8);
  legtr->SetBorderSize(0);
  legtr->SetFillStyle(0);

  TH1D* htrigkAny=new TH1D();
  htrigkAny = coll->Get("run",Form("triggertype:%s","ANY"));
  htrigkAny->SetName(Form("h%s","ANY"));
  htrigkAny->SetTitle("Trigger type;RUN; counts");
  htrigkAny->Sumw2();
  TH1D* htrigkAll=new TH1D();
  htrigkAll = coll->Get("run",Form("triggertype:%s","ALL"));
  htrigkAll->SetName(Form("h%s","ALL"));
  htrigkAll->SetTitle("Trigger type;RUN; counts");
  htrigkAll->Sumw2();

  for(Int_t k=0;k<nkeys;k++){
    htrig[k]=coll->Get("run",Form("triggertype:%s",words[k].Data()));
    htrig[k]->SetName(Form("h%s",words[k].Data()));
    htrig[k]->SetTitle("Trigger type;RUN; counts");
    htrig[k]->SetMarkerColor(k+1);
    htrig[k]->SetMarkerStyle(k+20);
    htrig[k]->Sumw2();
    legtr->AddEntry(htrig[k],Form("%s",words[k].Data()),"P");
    //drawings
    //1) counts of a given trigger over counts in kAll
    htrignormAll[k]=(TH1D*)htrig[k]->Clone(Form("h%snormAll",words[k].Data()));
    htrignormAll[k]->SetTitle("Trigger type over ALL trigger;RUN; counts/countsALL");
    htrignormAll[k]->Divide(htrig[k],htrigkAll,1.,1.,"B");
    htrignormAll[k]->GetXaxis()->SetRangeUser(0,1.1);

    ctrigfractionAll->cd();
    if(k>0)htrignormAll[k]->Draw("PEsames");
    else htrignormAll[k]->Draw("PE");

    //1) counts of a given trigger over counts in kAny
    htrignormAny[k]=(TH1D*)htrig[k]->Clone(Form("h%snormAny",words[k].Data()));
    htrignormAny[k]->SetTitle("Trigger type over ANY trigger;RUN; counts/countsANY");
    htrignormAny[k]->Divide(htrig[k],htrigkAny,1.,1.,"B");
    htrignormAny[k]->GetXaxis()->SetRangeUser(0,1.1);
    
    ctrigfractionAny->cd();
    if(k>0)htrignormAny[k]->Draw("PEsames");
    else htrignormAny[k]->Draw("PE");
  }

  ctrigfractionAll->cd();
  legtr->Draw();
  ctrigfractionAll->SaveAs("TrgFractionOverALL.pdf");
  ctrigfractionAll->SaveAs("TrgFractionOverALL.eps");
  
  ctrigfractionAny->cd();
  legtr->Draw();
  ctrigfractionAny->SaveAs("TrgFractionOverANY.pdf");
  ctrigfractionAny->SaveAs("TrgFractionOverANY.eps");

  delete [] words;
  delete htrigkAll; delete htrigkAny;

  //draw number of events per trigger
  TH2F* hTrigMul=(TH2F*)list->FindObject("hTrigMul");

  TH1F *hnEvPerTrig=(TH1F*)hTrigMul->ProjectionX("hnEvPerTrig");
  hnEvPerTrig->SetTitle("Number of events per trigger");
  TCanvas* cnev=new TCanvas("cnev","Number of events",1400,800);
  cnev->cd();
  hnEvPerTrig->Draw("htext0");
  cnev->SaveAs(Form("%s.pdf",cnev->GetName()));
  cnev->SaveAs(Form("%s.eps",cnev->GetName()));

  //draw number of events selected per trigger
  TH2F* hTrigMulSel=(TH2F*)list->FindObject("hTrigMulSel");

  TH1F *hnEvPerTrigSel=(TH1F*)hTrigMulSel->ProjectionX("hnEvPerTrigSel");
  hnEvPerTrigSel->SetTitle("Number of events selected per trigger");
  TCanvas* cnevs=new TCanvas("cnevs","Number of events selected",1400,800);
  cnevs->cd();
  hnEvPerTrigSel->Draw("htext0");
  cnevs->SaveAs(Form("%s%s.pdf",cnevs->GetName(),suffixdir.Data()));
  cnevs->SaveAs(Form("%s%s.eps",cnevs->GetName(),suffixdir.Data()));

  TH1F* hnEvSelPerc=(TH1F*)hnEvPerTrigSel->Clone("hnEvSelFrac");
  hnEvSelPerc->SetTitle("Fraction of event selected per trigger");
  hnEvSelPerc->Divide(hnEvPerTrig);
  TCanvas* cnevsp=new TCanvas("cnevsp","Fraction of events selected",1400,800);
  cnevsp->cd();
  hnEvSelPerc->Draw("htext0");
  cnevsp->SaveAs(Form("%s%s.pdf",cnevsp->GetName(),suffixdir.Data()));
  cnevsp->SaveAs(Form("%s%s.eps",cnevsp->GetName(),suffixdir.Data()));
}

void WriteTextFileForCompilation(TString inputtextfile, Int_t analysistype){
   //This method writes a text file with the name given in input containing the path of the QA output file, the directories and the lists of histograms to be read to perform a comparison of QA outputs
   // analysistype = 1 is the OutputTrack
   // analysistype = 2 is the OutputEventSelection
   // other outputs can be implemented

   // The paths and trigger type (or more in general the suffix of the directories to be read have to be set in this method)
   // The legend is set either equal to the trigger type if they are different, or it should be set manually
   // Run this methos before CompilationOfTriggeredEvents unless the text file has been already produced
   
  Printf("WARNING! Did you customize the parameters in the macro?");
  const Int_t n=8; //customize this
  TString prepath="/data/Work/HFQA/LHC12QATrainOutputs/";//"./";
  TString pathname[n]={"LHC12a/Train165/","LHC12b/Train166/","LHC12d/Train168/","LHC12e/Train170/","LHC12f/Train171/","LHC12g/Train172/"};

  TString dirname[n];
  TString listname[n];
  ofstream myfile;
  myfile.open(inputtextfile);
  myfile<<n<<endl;

  TString filename="AnalysisResults.root";
  TString basedirname="PWG3_D2H_QA";
  TString trgnames[n];//={"INT7","EMC7","EMCJET7","EMCGAMMA7","INT8","EMC8","EMCJET8","EMCGAMMA8"}; //customize this
  TString baselistname="", partname="D00100";
  TString legend[n]={"LHC12a165","LHC12b166","LHC12d168","LHC12e170","LHC12f171","LHC12g172"}; //Customize this if want it to be different from trgnames
  if(analysistype==1) baselistname="outputTrack";
  if(analysistype==2) baselistname="outputEvSel";
  //... others could be added

  //All in the same directory on your computer
  for(Int_t i=0;i<n;i++){
    //set
    pathname[i]=prepath+pathname[i];//"./";
    trgnames[i]="EMC7";
    dirname[i]=basedirname+trgnames[i];
    listname[i]=baselistname+partname+trgnames[i];
    if(legend[i]=="") legend[i]=trgnames[i];
    
    //Printf("Trigger name is %s",trgnames[i].Data());
    //Printf("Legend is %s",legend[i].Data());
    
    //write text file
    myfile<<endl;
    myfile<<endl;
    myfile<<pathname[i].Data()<<filename.Data()<<endl;
    myfile<<dirname[i].Data()<<endl;
    myfile<<listname[i].Data()<<endl;
    myfile<<legend[i].Data()<<endl;
      
  }
  myfile.close();

}

Bool_t ReadFilesForCompilation(TString inputtextfile, TList**& lists, Int_t& numb, TString*& legend){

   //Method to read the QA files indicated in the text file in input (to be written with WriteTextFileForCompilation() method)
  ifstream myfile;
  myfile.open(inputtextfile);
  if(!myfile.is_open()) Printf("No files found, did you make it correctly?");
  Int_t n; 
  myfile >> n;
  lists=new TList*[n];
  legend= new TString[n];
  numb=n;

  Int_t k=0;

  while(myfile){
   
    TString filename;
    myfile>>filename;
    //Printf("Filename = %s",filename.Data());
    TFile *fin=new TFile(filename);
    if(!fin->IsOpen()){
      Printf("File %s not found",filename.Data());
      return kFALSE;
    }
   
    TString dirname;
    myfile>>dirname;
    //Printf("Dirname = %s",dirname.Data());
    TDirectoryFile* dir=(TDirectoryFile*)fin->Get(dirname);
    if(!dir){
      Printf("Directory %s not found in %s",dirname.Data(),filename.Data());
      fin->ls();
      return kFALSE;
    }
   
    TString listname;
    myfile>>listname;
    //Printf("Listname = %s",listname.Data());
    lists[k]=(TList*)dir->Get(listname);
    if(!lists[k]){
      Printf("List %s not found in %s",listname.Data(),dirname.Data());
      dir->ls();
      return kFALSE;
    }
    
    TString temp;
    myfile>>temp;
    legend[k]=temp;
    Printf("Legend = %s",legend[k].Data());
    if(!legend[k]){
      Printf("Legend %s not found",filename.Data());
      return kFALSE;
    }
    

    k++;
    if(k==n)return kTRUE;// Needed as the while loop does not always realise the end is reached. 
  }
  return kTRUE;
}

void CompilationOfTriggeredEvents(TString inputtextfile,Int_t analysistype){
   //This method draws the histograms superimposed

   // analysistype = 1 is the OutputTrack
   // analysistype = 2 is the OutputEventSelection
   // other outputs can be implemented

  TList **lists=0x0;
  Int_t n;
  TString* legend=0x0;
  Bool_t okread=ReadFilesForCompilation(inputtextfile,lists,n,legend);
  if(!okread){
    Printf("Did you write %s with  WriteTextFileForCompilation(%s)?",inputtextfile.Data(),inputtextfile.Data());
    return;
  }

  if(analysistype==2){CompilationEvSelection(n,lists,legend);}
  if(analysistype==1){CompilationTrackSelection(n,lists,legend);}
}

void CompilationEvSelection(Int_t n,TList** lists, TString* legend){
   // Specific method for EventSelection output (2)

  TCanvas *cnevs=new TCanvas("cnevs","Number of events selected",1400,800);
  TCanvas *cnevsp=new TCanvas("cnevsp","Fraction of events selected",1400,800);
  TH1F** hnEvPerTrig=new TH1F*[n];
  TH1F** hnEvPerTrigSel=new TH1F*[n];
  TH1F** hnEvSelPerc=new TH1F*[n];
  
  TLegend* leg=new TLegend(0.15,0.5,0.45,0.78);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  Bool_t first=kTRUE;
  Int_t nentriesl=lists[0]->GetEntries();
  TCanvas** c=new TCanvas*[nentriesl];
  Int_t nth1f=0;
  
  for(Int_t i=0;i<n;i++){
   
    TH2F* hTrigMul=(TH2F*)lists[i]->FindObject("hTrigMul");

    hnEvPerTrig[i]=(TH1F*)hTrigMul->ProjectionX(Form("hnEvPerTrig%d",i));
    hnEvPerTrig[i]->SetTitle("Number of events selected per trigger");
   

    TH2F* hTrigMulSel=(TH2F*)lists[i]->FindObject("hTrigMulSel");
    
    hnEvPerTrigSel[i]=(TH1F*)hTrigMulSel->ProjectionX(Form("hnEvPerTrigSel%d",i));
    hnEvPerTrigSel[i]->SetTitle("Number of events selected per trigger");
    hnEvPerTrigSel[i]->SetLineColor(i+1);
    hnEvPerTrigSel[i]->SetLineWidth(2);

    cnevs->cd();
    if(i==0) hnEvPerTrigSel[i]->Draw();
    else  hnEvPerTrigSel[i]->Draw("sames");

    hnEvSelPerc[i]=(TH1F*)hnEvPerTrigSel[i]->Clone(Form("hnEvSelFrac%d",i));
    hnEvSelPerc[i]->SetTitle("Fraction of event selected per trigger");
    hnEvSelPerc[i]->SetLineColor(i+1);
    hnEvSelPerc[i]->SetLineWidth(2);
    hnEvSelPerc[i]->Divide(hnEvPerTrig[i]);
    cnevsp->cd();
    if(i==0) hnEvSelPerc[i]->Draw("htext0");
    else  hnEvSelPerc[i]->Draw("htext0sames");
    nth1f=0; //restart counting per each file 
    
    //fill legend
    leg->AddEntry(hnEvSelPerc[i],legend[i],"L"); 
    
    for(Int_t j=0; j<nentriesl; j++){
       TClass* objtype=lists[i]->At(j)->IsA();
       TString tpname=objtype->GetName();
       
       if (tpname=="TH1F"){
       	  TH1F* htmp=(TH1F*)lists[i]->At(j);
       	  
       	  htmp->SetLineColor(1+i);
       	  if(htmp->Integral()>0) htmp->Scale(1./htmp->Integral());

       	  if(first) {
       	     c[nth1f]=new TCanvas(Form("c%s",htmp->GetName()),Form("c%s",htmp->GetName()));       	     
       	     c[nth1f]->cd();
       	     htmp->Draw();
       	       
       	  }
       	  else {
       	     c[nth1f]->cd();
       	     htmp->Draw("sames");
       	  }
       	  nth1f++;
       }
    }

    first=kFALSE;
 
  }
  
  for(Int_t j=0;j<nth1f;j++){
     c[j]->cd();
     leg->Draw();
     c[j]->SaveAs(Form("%s.pdf",c[j]->GetName()));
     c[j]->SaveAs(Form("%s.eps",c[j]->GetName()));
  }
  
  cnevs->cd();
  leg->Draw();
  cnevs->SaveAs(Form("%s.eps",cnevs->GetName()));
  cnevs->SaveAs(Form("%s.pdf",cnevs->GetName()));
  
  cnevsp->cd();
  leg->Draw();
  cnevsp->SaveAs(Form("%s.eps",cnevsp->GetName()));
  cnevsp->SaveAs(Form("%s.pdf",cnevsp->GetName()));
  TCanvas *test=new TCanvas();
  test->cd(); leg->Draw();
}


void CompilationTrackSelection(Int_t n,TList** lists,TString* legend)
{
   // Specific method for Track output (1)

  for(Int_t i=0;i<lists[0]->GetEntries();i++){
    TObjArray* plotseth=new TObjArray();
    TObjArray* plotsethr=new TObjArray();
    
    Bool_t is2d=false;
    is2d = lists[0]->At(i)->InheritsFrom("TH2");
    if(is2d) continue;
    
    for(Int_t j=0; j<n; j++){
      TH1F* temph=(TH1F*)lists[j]->At(i);
      TH1F * temphr=0x0;
      TH1F * ploto=0x0; //0th plot to divide by
      
      // plot ratios and scale all plots to integral
      temphr= (TH1F*)temph->Clone(Form("%s_ratio",temph->GetName()));
      if(j==0)ploto = temphr;
      else ploto = (TH1F*)plotseth->At(0);
      temphr->Divide(ploto);
      if(temphr->Integral()>0){temphr-> Scale(1./temphr->Integral());}// take out for nonscaled ratios
      if(temph->Integral()>0){ temph-> Scale(1./temph->Integral());}
      
      plotseth->AddLast(new TH1F(*temph));
      plotsethr->AddLast(new TH1F(*temphr));
    }

    
    TH1F *h = (TH1F*)plotseth->At(0);
    TCanvas* c=new TCanvas(Form("c%s",h->GetName()),h->GetName());
    c->cd();
    c->SetGrid();
    TString hname=h->GetName();
    
    if(!hname.Contains("nCls")){
      c->SetLogy();
      if(hname.Contains("Layer")){
	for(Int_t ibin=1;ibin<=h->GetNbinsX();ibin++){
	  h->GetXaxis()->SetBinLabel(ibin+1,Form("%d",ibin));
	}
	h->GetXaxis()->SetLabelSize(0.06);
	h->GetXaxis()->SetRangeUser(0,6); //comment to see ntracks!
      }
      //h->SetMinimum(1);
      h->Draw();
      
      TLegend* leg=new TLegend(0.15,0.5,0.45,0.78);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h,legend[0],"L");
      
      for(Int_t j=1; j<n; j++){
	TH1F * h2 = (TH1F*)plotseth->At(j);
	h2->SetLineColor(1+j);
	h2->Draw("sames");
	leg->AddEntry(h2,legend[j],"L");
      }
      
      leg->Draw();
      
      TCanvas* c2=new TCanvas(Form("c2%s",h->GetName()),h->GetName());
      c2->cd();
      c2->SetGrid();
      TH1F *hr = (TH1F*)plotsethr->At(1);
      hr->SetLineColor(2);
      hr->Draw();
      TLegend* leg2=new TLegend(0.15,0.5,0.45,0.78);
      leg2->SetFillStyle(0);
      leg2->SetBorderSize(0);
      
      TString ratioleg;
      ratioleg+=legend[1];
      ratioleg+="/";
      ratioleg+=legend[0];
    
      leg2->AddEntry(hr,ratioleg,"L");
      
      for(Int_t j=2; j<n; j++){
	TH1F * hr2 = (TH1F*)plotsethr->At(j);
	hr2->SetLineColor(1+j);
	hr2->Draw("sames");
      
	TString ratioleg2;
	ratioleg2+=legend[j];
	ratioleg2+="/";
	ratioleg2+=legend[0];

	leg2->AddEntry(hr2,ratioleg2,"L");
      }
      leg2->Draw();
      c2->SaveAs(Form("%s%s%sRatio.pdf",c->GetName(),legend[0].Data(), legend[n-1].Data()));
      c2->SaveAs(Form("%s%s%sRatio.eps",c->GetName(),legend[0].Data(), legend[n-1].Data()));
    } 
    else {
      h->Draw("htext0");
    
      TLegend* leg=new TLegend(0.15,0.5,0.45,0.78);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h,legend[0],"L");
      
      for(Int_t j=1; j<n; j++){
	TH1F * h2 = (TH1F*)plotseth->At(j);
	h2->SetLineColor(1+j);
	leg->AddEntry(h2,legend[j],"L");
	h2->Draw("htext0sames");
      }
      
      leg->Draw();
    }
    c->cd();
    c->SaveAs(Form("%s%s%s.eps",c->GetName(),legend[0].Data(), legend[n-1].Data()));
    c->SaveAs(Form("%s%s%s.pdf",c->GetName(),legend[0].Data(), legend[n-1].Data()));		
  }
  
  Int_t check=0;
  Int_t maxfa=1;
  if(n<=10)maxfa=n;
  else cout<<"Warning: To many files for combinationplot, show only original"<<endl;
  TLegend* leg3=new TLegend(0.15,0.5,0.45,0.78);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  TCanvas* ctrsel=new TCanvas("ctrsel","Track Sel");
  ctrsel->SetLogy();
  for(Int_t i=0; i<maxfa; i++){
    TH1F* hd0fb4=(TH1F*)lists[i]->FindObject("hd0TracksFilterBit4");
    TH1F* hd0SPD1=(TH1F*)lists[i]->FindObject("hd0TracksSPDin");
    TH1F* hd0SPDany=(TH1F*)lists[i]->FindObject("hd0TracksSPDany");
    TH1F* hd0TPCITScuts=(TH1F*)lists[i]->FindObject("hd0TracksTPCITSSPDany");
    if(hd0fb4 && hd0SPD1 && hd0SPDany && hd0TPCITScuts){
      if(i==0){check=1;}
      else{if(check==0)return;}
      ctrsel->cd();
      hd0SPD1->SetLineColor(kCyan+3);
      hd0SPDany->SetLineColor(4);
      hd0TPCITScuts->SetLineColor(kGreen+1);
      hd0fb4->SetLineColor(2);
      if(i==0){
	hd0SPD1->Draw();
	ctrsel->Update();
	TPaveStats *st1=(TPaveStats*)hd0SPD1->GetListOfFunctions()->FindObject("stats");
	st1->SetTextColor(kCyan+3);
	st1->SetY1NDC(0.71);
	st1->SetY2NDC(0.9);
	hd0SPDany->Draw("sames");
	ctrsel->Update();
	TPaveStats *st2=(TPaveStats*)hd0SPDany->GetListOfFunctions()->FindObject("stats");
	st2->SetY1NDC(0.51);
	st2->SetY2NDC(0.7);
	st2->SetTextColor(4);
	hd0fb4->Draw("sames");
	ctrsel->Update();
	TPaveStats *st3=(TPaveStats*)hd0fb4->GetListOfFunctions()->FindObject("stats");
	st3->SetY1NDC(0.31);
	st3->SetY2NDC(0.5);
	st3->SetTextColor(2);
	hd0TPCITScuts->Draw("sames");
	ctrsel->Update();
	TPaveStats *st4=(TPaveStats*)hd0TPCITScuts->GetListOfFunctions()->FindObject("stats");
       st4->SetY1NDC(0.71);
       st4->SetY2NDC(0.9);
       st4->SetX1NDC(0.55);
       st4->SetX2NDC(0.75);
       st4->SetTextColor(kGreen+1);
       ctrsel->Modified();
       leg3->AddEntry(hd0SPD1,"kITSrefit+SPD inner","L");
       leg3->AddEntry(hd0SPDany,"kITSrefit+SPD any","L");
       leg3->AddEntry(hd0TPCITScuts,"TPC+ITS cuts+SPD any","L");
       leg3->AddEntry(hd0fb4,"Filter Bit 4","L");
       leg3->AddEntry(hd0SPD1, legend[i], "L");
      }
     else{	
       hd0SPD1->SetStats(0);
       hd0SPD1->SetLineStyle(i+1);
       hd0SPD1->Draw("sames"); 
       hd0SPDany->SetStats(0);
       hd0SPDany->SetLineStyle(i+1);
       hd0TPCITScuts->SetStats(0);
       hd0TPCITScuts->SetLineStyle(i+1);
       hd0fb4->SetStats(0);
       hd0fb4->SetLineStyle(i+1);
       ctrsel->cd();
       hd0SPD1->Draw("sames");
       hd0SPDany->Draw("sames");
       hd0TPCITScuts->Draw("sames");
       hd0fb4->Draw("sames");
       leg3->AddEntry(hd0SPD1, legend[i], "L");				
     }
     
    }
  }
  leg3->Draw();
  ctrsel->SaveAs("ImpactParameterTrackSel.eps");
  ctrsel->SaveAs("ImpactParameterTrackSel.pdf");
}

