///==========================================================================
///
///    macro to plot centrality bin values 
///==========================================================================
///
#if !defined(__CINT__) || defined(__MAKECINT__)

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <Riostream.h>
#include <TClassTable.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLine.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TGraphErrors.h>

#endif

void getCentrality(TH1 *histNch, Float_t ff);

//double centPercent[]={5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100};
double centPercent[]={5.,10.,20.,40.,60.,80.,100.};
//double centPercent[]={20.,40.,60.,80.,100.};
//double centPercent[]={5.,10.,20.,40.,60.,100.};
//double centPercent[]={20.,40.,60.,100.};
//double centPercent[]={60,100.};
const Int_t nbins  = sizeof(centPercent)/sizeof(double); 

TArrayF* binUp = new TArrayF(nbins);
TArrayF* Multbin = new TArrayF(nbins);

TH1F *h= new TH1F("h","h;Centrality [%]; Multiplicity per N_{ancestor}/<m>",nbins,0,100);

void plotGlauberCenVars(Float_t eff=1., const Char_t* file="ZNA_ntuple_195483.root")
{
  TFile *f = TFile::Open(file);
  TNtuple* ntuple = dynamic_cast<TNtuple*> (f->Get("gnt"));
   
  TGraphErrors *gNpart=new TGraphErrors(0);
  gNpart->SetName("gNpart"); 
  TGraphErrors *gNcoll=new TGraphErrors(0);
  gNcoll->SetName("gNcoll"); 
  TGraphErrors *gtAA=new TGraphErrors(0);
  gtAA->SetName("gtAA"); 
  
  /*TFile *ffd = TFile::Open("hZNAcalibRUN195483.root");
  TH1F * hd = dynamic_cast<TH1F*> (ffd->Get(("hZNA")));
  hd->Sumw2();*/
  //
  TFile *ff = TFile::Open("ZNA_fit_195483.root");
  TH1F * hd = dynamic_cast<TH1F*> (ff->Get(("hZNA")));
  hd->Sumw2();
  TH1F * hg = dynamic_cast<TH1F*> (ff->Get(("hZNA_GLAU")));
  hd->SetMarkerColor(kBlue+3);
  hd->SetMarkerSize(1.);
  hd->SetLineColor(kBlue+3);
  hd->SetLineWidth(2);
  hd->SetMarkerStyle(20);
  hd->SetLineWidth(2);
//  hg->Scale(1./hd->GetEntries());
//  hd->Scale(1./hd->GetEntries());
  hd->SetMinimum(1.e-01);
  hd->SetXTitle("E_{ZNA} (TeV)");
  hg->SetLineColor(kPink-2);
  hg->SetLineWidth(2);
  
  TH1F* hist = (TH1F*) hg->Clone("hist");

  //---------------------------------------------------
  getCentrality(hist, eff);
  //---------------------------------------------------

  TCanvas* canvas = new TCanvas("canvas","Multiplicity",200,200,600,600);
  canvas->cd();
  canvas->SetLogy();
  hd->Draw("pe");
  //hd->GetXaxis()->SetRangeUser(0.,130.);
  hd->SetMinimum(0.01);
  hg->Draw("SAME");

  float low = 0;
  float high = hist->GetNbinsX();
  for(int i=0; i<binUp->GetSize(); i++){
      low = binUp->At(i);
      hist->GetXaxis()->SetRange(low+1, high);
      hist->SetFillColor((i%2==0)?0:kAzure+6);
      hist->SetLineColor((i%2==0)?0:kAzure+6);
      printf(" bin %d  low %f  high %f\n",i,low,high);
      hist->DrawCopy("h same");
      high=low;
  }
  hd->Draw("esame");
  hg->Draw("SAME");
  canvas->Print("plotGlauber.gif");   
  
  TCanvas* canvas2 = new TCanvas("canvas2","NPart");
  canvas2->cd();
  canvas2->SetLogy();
  TH1F *hist2 = new TH1F("hist2","N_{part}",35,0.,35);
  ntuple->Project("hist2","fNpart");
  //hist2->SetStats(0);
  hist2->SetTitle("");
  hist2->GetXaxis()->SetTitle("NPart");
  hist2->GetXaxis()->SetTitleSize(0.05);
  hist2->GetXaxis()->SetLabelSize(0.04);
  hist2->GetXaxis()->SetTitleOffset(1.2);
  hist2->GetYaxis()->SetTitle("");
  hist2->GetYaxis()->SetTitleOffset(1.3);
  hist2->GetYaxis()->SetTitleSize(0.05);
  hist2->GetYaxis()->SetLabelSize(0.04);
  hist2->DrawCopy();
  
  float lownp=0;
  float highnp=5000;
  TH1F *htemp10[nbins];
  printf("\n ***** N_part \n");
  for(int i=0; i<Multbin->GetSize(); i++){
      lownp = Multbin->At(i);
      char cuts[120];
      char histname[20];
      sprintf(cuts,"Etot>%f && Etot<=%f",lownp,highnp);
      sprintf(histname,"htemp10[%i]",i);
      htemp10[i] = new TH1F(histname,"N_{part}",35,0.,35);
      //printf(" cut: %s\n", cuts);
      ntuple->Project(histname,"fNpart",cuts);
      htemp10[i]->SetLineColor(i+1);
      htemp10[i]->Draw("same");
      cout  << i << " | " << lownp << " | " << highnp << " | " << setprecision(3) << 
      htemp10[i]->GetMean() << " | " << htemp10[i]->GetRMS() << " | " << endl;
      gNpart->SetPoint(i,Float_t(i),htemp10[i]->GetMean());
      gNpart->SetPointError(i,0,htemp10[i]->GetRMS());
      highnp = lownp;
  }
  cout << endl;
  
  TCanvas* canvas3 = new TCanvas("canvas3","NColl");
  canvas3->SetLogy();
  TH1F *hist3 = new TH1F("hist3","N_{coll}",35,0.,35);
  ntuple->Project("hist3","fNcoll");
  //hist3->SetStats(0);
  hist3->SetTitle("");
  hist3->GetXaxis()->SetTitle("NColl");
  hist3->GetXaxis()->SetTitleSize(0.05);
  hist3->GetXaxis()->SetLabelSize(0.04);
  hist3->GetXaxis()->SetTitleOffset(1.2);
  hist3->GetXaxis()->SetTitle("");
  hist3->GetXaxis()->SetTitleOffset(1.3);
  hist3->GetXaxis()->SetTitleSize(0.05);
  hist3->GetXaxis()->SetLabelSize(0.04);
  hist3->DrawCopy();
  
  float lownc = 0;
  float highnc = 5000;
  TH1F *htemp11[nbins];
  printf("\n ***** N_coll \n");
  for(int i=0; i<Multbin->GetSize(); i++){
      lownc = Multbin->At(i);
      char cuts[120];
      char histname[20];
      sprintf(cuts,"Etot>%f && Etot<=%f",lownc,highnc);
      sprintf(histname,"htemp11[%i]",i);
      htemp11[i] = new TH1F(histname,"N_{coll}",35,0.,35.);
      ntuple->Project(histname,"fNcoll",cuts);
      htemp11[i]->SetLineColor(i+1);
      htemp11[i]->Draw("same");
      cout << setprecision(3) << htemp11[i]->GetMean() << " | " << htemp11[i]->GetRMS() << " | " << endl;
      gNcoll->SetPoint(i,Float_t(i),htemp11[i]->GetMean());
      gNcoll->SetPointError(i,0,htemp11[i]->GetRMS());
      highnc = lownc;
  }
  cout << endl;
  
  TCanvas* canvas4 = new TCanvas("canvas4","Impact Parameter");
  canvas4->cd();
  TH1F *hist4 = new TH1F("hist4","b",100,0.,16.);
  ntuple->Project("hist4","fB");
  //hist4->SetStats(0);
  hist4->SetTitle("");
  hist4->GetXaxis()->SetTitle("b");
  hist4->GetXaxis()->SetTitleSize(0.05);
  hist4->GetXaxis()->SetLabelSize(0.04);
  hist4->GetXaxis()->SetTitleOffset(1.2);
  hist4->GetYaxis()->SetTitle("");
  hist4->GetYaxis()->SetTitleOffset(1.3);
  hist4->GetYaxis()->SetTitleSize(0.05);
  hist4->GetYaxis()->SetLabelSize(0.04);
  hist4->DrawCopy();
  
  float lowb = 0;
  float highb = 5000;
  TH1F *htemp12[nbins];
  printf("\n ***** b \n");
  for(int i=0; i<Multbin->GetSize(); i++){
      lowb = Multbin->At(i);
      char cuts[100];
      char histname[25];
      sprintf(cuts,"Etot>%f && Etot<=%f",lowb,highb);
      sprintf(histname,"htemp12[%i]",i);
      htemp12[i] = new TH1F(histname,"b",100,0.,16.);
      //printf(" cut: %s\n", cuts);
      ntuple->Project(histname,"fB",cuts);
      htemp12[i]->SetLineColor(i+1);
      htemp12[i]->DrawCopy("same");
      cout << i << " | " << lowb << " | " << highb << " | " << setprecision(3) << htemp12[i]->GetMean() << " | " << htemp12[i]->GetRMS() << " | " << endl;
      highb = lowb;
  }
  
  TCanvas* canvas5 = new TCanvas("canvas5","Taa");
  canvas5->SetLogy();
  TH1F *hist5 = new TH1F("hist5","T_{AA}",100,0.,0.5);
  ntuple->Project("hist5","fTaa");
  //hist5->SetStats(0);
  hist5->SetTitle("");
  hist5->GetXaxis()->SetTitle("tAA");
  hist5->GetXaxis()->SetTitleSize(0.05);
  hist5->GetXaxis()->SetLabelSize(0.04);
  hist5->GetXaxis()->SetTitleOffset(1.2);
  hist5->GetYaxis()->SetTitle("");
  hist5->GetYaxis()->SetTitleOffset(1.3);
  hist5->GetYaxis()->SetTitleSize(0.05);
  hist5->GetYaxis()->SetLabelSize(0.04);
  hist5->DrawCopy();
  
  float lowtaa = 0;
  float hightaa = 5000;
  TH1F *htemp13[nbins];
  printf("\n ***** T_AA \n");
  for (int i=0; i<Multbin->GetSize(); i++){
      lowtaa = Multbin->At(i);
      char cuts[100];
      char histname[100];
      sprintf(cuts,"Etot>%f && Etot<%f",lowtaa,hightaa);
      //printf(" cut: %s\n", cuts);
      sprintf(histname,"htemp13[%i]",i);
      htemp13[i] = new TH1F(histname,"b",100,0.,0.5);
      ntuple->Project(histname,"fTaa",cuts);
      htemp13[i]->SetLineColor(i+1);
      htemp13[i]->DrawCopy("same");
      cout << setprecision(3) << htemp13[i]->GetMean() << " | " << htemp13[i]->GetRMS() << " | " << endl;  
      gtAA->SetPoint(i,Float_t(i),htemp13[i]->GetMean());
      gtAA->SetPointError(i,0,htemp13[i]->GetRMS());
      hightaa = lowtaa;
  }

  /*TCanvas* canvas6 = new TCanvas("canvas6","Mean Mult");
  canvas6->SetLogy();
  //ntuple->Draw("ntot/Npart/23.867>>hmultperanc");
  ntuple->Draw("Etot/(0.801*Npart+(1-0.801)*Ncoll)/3.9>>hmultperanc");
  TH1F* hmultperanc = (TH1F*)gPad->GetPrimitive("hmultperanc");
  TH1F* hist6 = (TH1F*)hmultperanc->Clone();
  hist6->SetStats(0);
  hist6->SetTitle("");
  hist6->GetXaxis()->SetTitle("Mult/NPart");
  hist6->GetXaxis()->SetTitleSize(0.05);
  hist6->GetXaxis()->SetLabelSize(0.04);
  hist6->GetXaxis()->SetTitleOffset(1.);
  hist6->GetYaxis()->SetTitle("");
  hist6->GetYaxis()->SetTitleOffset(1.);
  hist6->GetYaxis()->SetTitleSize(0.05);
  hist6->GetYaxis()->SetLabelSize(0.04);
  hist6->DrawCopy();
  
  low=0;
  high=50000;
  for (int i=0; i<Multbin->GetSize(); i++)
    {
      low=Multbin->At(i);
      char cuts[100];
      char histtitle1[100];
      char histtitle[100];
      sprintf(cuts,"Etot>%i && Etot<%i",low,high);
      //sprintf(histtitle,"ntot/Npart/23.867>>htemp%i(100)",80+i);
      sprintf(histtitle,"Etot/(0.801*Npart+(1-0.801)*Ncoll)/3.9>>htemp%i(100)",80+i);
      sprintf(histtitle1,"htemp%i",80+i);
      ntuple->Draw(histtitle,cuts,"same");
      TH1F* htemp14 = (TH1F*)gPad->GetPrimitive(histtitle1);
      htemp14 = (TH1F*)gPad->GetPrimitive(histtitle1);
      htemp14->SetFillColor(i+1);
      htemp14->DrawCopy("same");
      //cout  << i << " | " << low << " | " << high << " | " << setprecision(3) << htemp14->GetMean() << " | " << htemp14->GetRMS() << " | " << endl;
      high=low;
      
      h->SetBinContent(i+1,htemp14->GetMean());
      h->SetBinError(i+1,0.01*htemp14->GetRMS());  
    }
  
  TCanvas* c7 = new TCanvas("c7","c7");
  c7->cd();
  h->Draw();
  */
  
  //TFile *outrootfile = new TFile("OutGlauber_Hijing.root","RECREATE");
  TFile *outrootfile = new TFile("test.root","RECREATE");
  outrootfile->cd();
  //h->Write();
  gNpart->Write();
  gNcoll->Write();
  gtAA->Write();
  hd->Write();
  outrootfile->Close();


}


void getCentrality(TH1 *histNch, Float_t ff)
// histNch - histo of multiplicity distribution (better with binsize=1)
// ff fraction of accepted events. All losses are assumed to occur in most
// peripheral bin
{
  
  double sum = histNch->Integral(); 
  int nbinsx = int (histNch->GetNbinsX());
  printf(" ZNA histo has %d bins in range: %1.2f-%1.2f TeV -> integral %f entries %f\n\n",
        nbinsx, histNch->GetBinLowEdge(1), 
	histNch->GetBinLowEdge(nbinsx)+histNch->GetBinWidth(nbinsx),
	sum, histNch->GetEntries());
  double accumulo=0., frac=0.;
  int ic=0;
  for(int ib=nbinsx; ib>0; ib--){
    Double_t content = histNch->GetBinContent(ib);
    accumulo += content;
    frac = accumulo/sum*100.*ff;
    //if(content>0.) printf(" bin %d  x %f content %f cumulative %f fraction %f \n",
      //     ib, histNch->GetBinCenter(ib), content, accumulo, frac);
    if(frac >= centPercent[ic]){
      binUp->SetAt(ib, ic);
      Multbin->SetAt(histNch->GetBinLowEdge(ib), ic);
      cout<<"   -> centrality = "<<centPercent[ic]<<"   ZN >= "<< histNch->GetBinLowEdge(ib) <<endl;
      ic++;
    }
    if(ic==nbins) break;
  }
  printf(" \n float multCent[%i] = {",nbins);
  // cout <<" \n float multCent[nbins] = {";
  
  for(int ich=nbins-1; ich>-1; ich--){
    cout<< histNch->GetBinLowEdge(binUp->At(ich));
    if(ich!=0) cout<<", ";
  }
  cout<<"};\n"<<endl;
}
