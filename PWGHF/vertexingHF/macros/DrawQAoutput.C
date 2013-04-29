#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
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
#endif

TString *periodsname;

//read the file and take list and stat

Bool_t ReadFile(TList* &list,TH1F* &hstat, TString listname,TString partname,TString path="./",TString filename=/*"AnalysisResults.root"*/"PWG3histograms.root", TString dirname="PWG3_D2H_QA");
Bool_t ReadFileMore(TList* &list,TH1F* &hstat, AliRDHFCuts* &cutobj, TString listname,TString partname,TString path="./",TString filename=/*"AnalysisResults.root"*/"PWG3histograms.root", TString dirname="PWG3_D2H_QA");
void SuperimposeBBToTPCSignal(Int_t period /*0=LHC10bc, 1=LHC10d, 2=LHC10h*/,TCanvas* cpid, Int_t set);
void TPCBetheBloch(Int_t set);

Bool_t ReadFile(TList* &list,TH1F* &hstat, TString listname,TString partname,TString path,TString filename, TString dirname){

  TString hstatname="nEntriesQA", cutobjname="";
  filename.Prepend(path);
  listname+=partname;
  hstatname+=partname;

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

  hstat=(TH1F*)dir->Get(hstatname);
  if(!hstat){
    cout<<hstatname.Data()<<" not found"<<endl;
    return kFALSE;
  }

  return kTRUE;
}

Bool_t ReadFileMore(TList* &list,TH1F* &hstat, AliRDHFCuts* &cutobj, TString listname,TString partname,TString path,TString filename,TString dirname){

  TString hstatname="nEntriesQA", cutobjname="";
  filename.Prepend(path);
  listname+=partname;
  hstatname+=partname;

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
	    else{
	      if(partname.Contains("LcToV0x")) cutobjname="LctoV0AnalysisCuts";
	    }
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

  hstat=(TH1F*)dir->Get(hstatname);
  if(!hstat){
    cout<<hstatname.Data()<<" not found"<<endl;
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
void DrawOutputTrack(TString partname="D0",TString textleg="",TString path="./", Bool_t superimpose=kFALSE, TString suffixdir="",TString filename=/*"AnalysisResults.root"*/"PWG3histograms.root"){
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetPalette(1);

  TString listname="outputTrack",name1="",name2="",path2="",filename2="PWG3histograms.root";
  TString tmp="y";

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
    }
    
  }

  TList* list;
  TH1F * hstat;
  TString dirname="PWG3_D2H_QA";
  dirname+=suffixdir;
  Bool_t isRead=ReadFile(list,hstat,listname,Form("%s%s",partname.Data(),name1.Data()),path,filename,dirname);
  if(!isRead) return;
  if(!list || !hstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }
  TPaveText *pvtxt=new TPaveText(0.6,0.6,0.9,0.9,"NDC");
  pvtxt->SetBorderSize(0);
  pvtxt->SetFillStyle(0);
  pvtxt->AddText(name1);

  TList* llist;
  TH1F* hhstat;
  if(superimpose){
    isRead=ReadFile(llist,hhstat,listname,Form("%s%s",partname.Data(),name2.Data()),path2,filename2,dirname);
    if(!isRead) return;
    if(!llist || !hhstat){
      cout<<":-( null pointers..."<<endl;
      return;
    }
    TText *redtext=pvtxt->AddText(name2);
    redtext->SetTextColor(kRed);
    hhstat->Scale(hstat->Integral()/hhstat->Integral());

  }

  for(Int_t i=0;i<list->GetEntries();i++){
    TH1F* h=(TH1F*)list->At(i);
    TH1F* hh=0x0;
    TH1F* hr=0x0;
    if(superimpose){
      hh=(TH1F*)llist->At(i);
      hr=(TH1F*)hh->Clone(Form("%s_ratio",hh->GetName()));
      hh->Scale(h->Integral()/hh->Integral());
    }
    if(!h || (superimpose && !hh)){
      cout<<"Histogram "<<i<<" not found"<<endl;
      continue;
    }
    if(superimpose){
      hh->Scale(h->Integral()/hh->Integral());
      hhstat->SetLineColor(kRed);
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
      h->Draw();
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
  
  TCanvas* cst=new TCanvas("cst","Stat");
  cst->SetGridy();
  cst->cd();
  hstat->Draw("htext0");
  if(superimpose) {
    hhstat->Draw("htext0sames");
    pvtxt->Draw();
  }
  cst->SaveAs(Form("%s%s.png",hstat->GetName(),textleg.Data()));
  cst->SaveAs(Form("%s%s.eps",hstat->GetName(),textleg.Data()));

  TH1F* hd0fb4=(TH1F*)list->FindObject("hd0TracksFilterBit4");
  TH1F* hd0SPD1=(TH1F*)list->FindObject("hd0TracksSPDin");
  TH1F* hd0SPDany=(TH1F*)list->FindObject("hd0TracksSPDany");
  TH1F* hd0TPCITScuts=(TH1F*)list->FindObject("hd0TracksTPCITSSPDany");
  if(hd0fb4 && hd0SPD1 && hd0SPDany && hd0TPCITScuts){
    TCanvas* ctrsel=new TCanvas("ctrsel","Track Sel");
    ctrsel->SetLogy();
    hd0SPD1->Draw();
    ctrsel->Update();
    TPaveStats *st1=(TPaveStats*)hd0SPD1->GetListOfFunctions()->FindObject("stats");
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
    
    ctrsel->SaveAs("ImpactParameterTrackSel.eps");
    ctrsel->SaveAs("ImpactParameterTrackSel.png");
    
  }
}

//draw "pid related" histograms (list "outputPID")
//period=-999 to draw the pull instead of the cut
void DrawOutputPID(TString partname="D0", Int_t mode=0/*0=with pull, 1=with nsigma*/,TString textleg="",TString path="./",TString suffixdir="", TString filename="AnalysisResults.root"){
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

  TList* list;
  TH1F * hstat;
  //needed only for mode 1
  AliRDHFCuts* cutobj;
  AliAODPidHF* aodpid;
  Double_t nsigmaTOF=0;
  Double_t nsigmaTPC[3]={},plimTPC[2]={};

  if(mode==1){
    Bool_t isRead=ReadFileMore(list,hstat,cutobj,listname,partname,path,filename,dirname);
    if(!isRead) return;
    if(!list || !hstat){
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
    Bool_t isRead=ReadFile(list,hstat,listname,partname,path,filename,dirname);
    if(!isRead) return;
    if(!list || !hstat){
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
	//mean and pull, code from Jens Wiechula
	TF1 fg("fg","gaus",-2.,2.); // fit range +- 2 sigma
	TLine l;
	TObjArray arr;

	//h->Draw("colz");
	fg.SetParameters(1,0,1);
	h->FitSlicesY(&fg,0,-1,0,"NQR",&arr);

	TH1 *hM=(TH1*)arr.At(1);
	hM->SetMarkerStyle(20);
	hM->SetMarkerSize(.5);
	hM->DrawClone("sames");

	TH1 *hS=(TH1*)arr.At(2);
	hS->SetMarkerStyle(20);
	hS->SetMarkerSize(.5);
	hS->SetMarkerColor(kRed);
	hS->SetLineColor(kRed);
	hS->DrawClone("same");

	l.SetLineColor(kBlack);
	l.DrawLine(.2,0,20,0);
	l.SetLineColor(kRed);
	l.DrawLine(.2,1,20,1);
	
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

  TList* list;
  TH1F * hstat;

  TString dirname="PWG3_D2H_QA",dirname2=dirname;
  dirname+=suffixdir;
  dirname2+=suffixdir2;
  Bool_t isRead=ReadFile(list,hstat,listname,partname.Data(),path,filename,dirname);
  if(!isRead) return;
  if(!list || !hstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }

  TPaveText *pvtxt=new TPaveText(0.6,0.6,0.9,0.9,"NDC");
  pvtxt->SetBorderSize(0);
  pvtxt->SetFillStyle(0);
  pvtxt->AddText(partname);

  TList* llist;
  TH1F* hhstat;
  if(superimpose){
    isRead=ReadFile(llist,hhstat,listname,partname2.Data(),path2,filename2,dirname2);
    if(!isRead) return;
    if(!llist || !hhstat){
      cout<<":-( null pointers..."<<endl;
      return;
    }
    TText *redtext=pvtxt->AddText(partname2);
    redtext->SetTextColor(kRed);

  }


  TCanvas* cst=new TCanvas("cst","Stat");
  cst->SetGridy();
  cst->cd();
  Int_t nevents=hstat->GetBinContent(1);
  hstat->Draw("htext0");
  cst->SaveAs(Form("%s%s.png",hstat->GetName(),textleg.Data()));
  cst->SaveAs(Form("%s%s.eps",hstat->GetName(),textleg.Data()));
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
	hhstat->SetLineColor(kRed);
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
      TH2F* h=(TH2F*)list->At(i);
      if(!h){
	cout<<"Histogram "<<i<<" not found"<<endl;
	continue;
      }
      TCanvas* c=new TCanvas(Form("c%s",h->GetName()),h->GetName());
      TPaveText *pvtxt3=new TPaveText(0.6,0.6,0.9,0.9,"NDC");
      pvtxt3->SetBorderSize(0);
      pvtxt3->SetFillStyle(0);

      c->cd();
      c->SetGrid();
      Int_t entries=h->Integral();
      pvtxt3->AddText(Form("%.1f %s of the events",(Double_t)entries/(Double_t)nevents*100,"%"));
      h->Draw("colz");
      c->SetLogz();
      pvtxt3->Draw();
      c->SaveAs(Form("%s%s.pdf",c->GetName(),textleg.Data()));
      c->SaveAs(Form("%s%s.eps",c->GetName(),textleg.Data()));
    }
  }
  
  
  listname="countersCentrality";

  isRead=ReadFile(list,hstat,listname,partname.Data(),path,filename,dirname);
  if(!isRead) return;
  if(!list || !hstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }


  if(superimpose){
    isRead=ReadFile(llist,hhstat,listname,partname2.Data(),path2,filename2,dirname2);
    if(!isRead) return;
    if(!llist || !hhstat){
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

  TList* list;
  TH1F * hstat;

  Bool_t isRead=ReadFile(list,hstat,listname,partname,path,filename,dirname);
  if(!isRead) return;
  if(!list || !hstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }
  Double_t nevents=hstat->GetBinContent(5); //ev good vertex

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

  TList* list;
  TH1F * hstat;

  Bool_t isRead=ReadFile(list,hstat,listname,partname,path,filename,dirname);
  if(!isRead) return;
  if(!list || !hstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }
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
  TH1D** htrignorm=new TH1D*[nkeys]; //normalized to the counts in kAny
  TCanvas* ctrigfraction=new TCanvas("cvtrigfrac","Fraction of given trigger type vs run",1400,800);
  TLegend* legtr=new TLegend(0.15,0.5,0.35,0.8);
  legtr->SetBorderSize(0);
  legtr->SetFillStyle(0);
  for(Int_t k=0;k<nkeys;k++){
    htrig[k]=coll->Get("run",Form("triggertype:%s",words[k].Data()));
    htrig[k]->SetName(Form("h%s",words[k].Data()));
    htrig[k]->SetTitle("Trigger type;RUN; counts");
    htrig[k]->SetMarkerColor(k+1);
    htrig[k]->SetMarkerStyle(k+20);
    htrig[k]->Sumw2();
    legtr->AddEntry(htrig[k],Form("%s",words[k].Data()),"P");
    //drawings
    //1) counts of a given trigger over counts in kAny
    htrignorm[k]=(TH1D*)htrig[k]->Clone(Form("h%snormAny",words[k].Data()));
    htrignorm[k]->SetTitle("Trigger type over ANY trigger;RUN; counts/countsANY");
    htrignorm[k]->Divide(htrig[k],htrig[0],1.,1.,"B");
    htrignorm[k]->GetXaxis()->SetRangeUser(0,1.1);
    
    ctrigfraction->cd();
    if(k>0)htrignorm[k]->Draw("PEsames");
    else htrignorm[k]->Draw("PE");
  } 
  
  ctrigfraction->cd();
  legtr->Draw();
  ctrigfraction->SaveAs("TrgFractionOverANY.pdf");
  ctrigfraction->SaveAs("TrgFractionOverANY.eps");

  delete [] words;

}
