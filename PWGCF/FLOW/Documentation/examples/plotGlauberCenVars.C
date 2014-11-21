///==========================================================================
///
///    macro to plot centrality bin values 
///==========================================================================
///
#include <cstdlib>

const Int_t nbins=11;
double centPercent[nbins]={5.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100};
TArrayI* binUp = new TArrayI(nbins);
TArrayI* Multbin = new TArrayI(nbins);

void plotGlauberCenVars(const Char_t* file="glau_pbpb_ntuple.root")
{
 TFile f(file);
 TNtuple* ntuple = f.Get("nt_Pb_Pb");

 TCanvas* canvas = new TCanvas("Multiplicity","Multiplicity");
 canvas->SetLogy();
 ntuple->Draw("dNdEtaNBD>>htemp(500)");
 TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
 TH1F* hist = (TH1F*)htemp->Clone();
 hist->SetTitle("");
 hist->SetStats(0);
 TAxis *axis = hist->GetXaxis();
 axis->SetTitle("Multiplicity");
 //axis->CenterTitle(kTRUE);
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 axis->SetTitleOffset(1.2);
 axis = hist->GetYaxis();
 axis->SetTitle("");
 axis->SetTitleOffset(1.3);
 //axis->CenterTitle(kTRUE);
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 hist->DrawCopy();

 //---------------------------------------------------
 getCentrality(hist);
 //---------------------------------------------------

 int low=0;
 int high=hist->GetNbinsX();
 //printf("size: %i\n",binUp->GetSize());
 for (int i=0; i<binUp->GetSize(); i++)
 {
   low=binUp->At(i);
   //printf("%i, %i\n",low,high);
   hist->GetXaxis()->SetRange(low+1,high);
   hist->SetFillColor((i%2==0)?0:kGray);
   hist->DrawCopy("same");
   high=low;
 }

 TCanvas* canvas2 = new TCanvas("N Wounded","N Wounded");
 canvas2->SetLogy();
 ntuple->Draw("Npart>>htemp2(200)");
 TH1F* htemp2 = (TH1F*)gPad->GetPrimitive("htemp2");
 TH1F* hist2 = (TH1F*)htemp2->Clone();
 hist2->SetStats(0);
 hist2->SetTitle("");
 axis = hist2->GetXaxis();
 axis->SetTitle("N Wounded");
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 axis->SetTitleOffset(1.2);
 axis = hist2->GetYaxis();
 axis->SetTitle("");
 axis->SetTitleOffset(1.3);
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 hist2->DrawCopy();


 int low=0;
 int high=5000;
 //printf("size: %i\n",Multbin->GetSize());
 TH1F* htemp10 = new TH1F();
 for (int i=0; i<Multbin->GetSize(); i++)
 {
   low=Multbin->At(i);
   //printf("%i, %i\n",low,high);
   char cuts[100];
   char histtitle1[100];
   char histtitle[100];
   sprintf(cuts,"dNdEtaNBD>%i && dNdEtaNBD<%i",low,high);
   //printf("dNdEtaNBD>%i && dNdEtaNBD<%i\n",low,high);
   sprintf(histtitle,"Npart>>htemp%i(200)",20+i);
   sprintf(histtitle1,"htemp%i",20+i);
   ntuple->Draw(histtitle,cuts,"same");
   htemp10 = (TH1F*)gPad->GetPrimitive(histtitle1);
   htemp10->SetLineColor(i+1);
   htemp10->DrawCopy("same");
   high=low;
 }

 TCanvas* canvas3 = new TCanvas("N Binary","N Binary");
 canvas3->SetLogy();
 ntuple->Draw("Ncoll>>htemp3(500)");
 TH1F* htemp3 = (TH1F*)gPad->GetPrimitive("htemp3");
 TH1F* hist3 = (TH1F*)htemp3->Clone();
 hist3->SetStats(0);
 hist3->SetTitle("");
 axis = hist3->GetXaxis();
 axis->SetTitle("N Binary");
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 axis->SetTitleOffset(1.2);
 axis = hist3->GetYaxis();
 axis->SetTitle("");
 axis->SetTitleOffset(1.3);
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 hist3->DrawCopy();

 int low=0;
 int high=5000;
 for (int i=0; i<Multbin->GetSize(); i++)
 {
   low=Multbin->At(i);
   char cuts[100];
   char histtitle[100];
   char histtitle1[100];
   sprintf(cuts,"dNdEtaNBD>%i && dNdEtaNBD<%i",low,high);
   sprintf(histtitle,"Ncoll>>htemp%i(500)",40+i);
   sprintf(histtitle1,"htemp%i",40+i);
   ntuple->Draw(histtitle,cuts,"same");
   TH1F* htemp11 = (TH1F*)gPad->GetPrimitive(histtitle1);
   htemp11->SetLineColor(i+1);
   htemp11->DrawCopy("same");
   high=low;
 }

 TCanvas* canvas4 = new TCanvas("Impact Parameter","Impact Parameter");
 ntuple->Draw("B>>htemp4(100)");
 TH1F* htemp4 = (TH1F*)gPad->GetPrimitive("htemp4");
 TH1F* hist4 = (TH1F*)htemp4->Clone();
 hist4->SetStats(0);
 hist4->SetTitle("");
 axis = hist4->GetXaxis();
 axis->SetTitle("b");
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 axis->SetTitleOffset(1.2);
 axis = hist4->GetYaxis();
 axis->SetTitle("");
 axis->SetTitleOffset(1.3);
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 hist4->DrawCopy();

 int low=0;
 int high=5000;
 for (int i=0; i<Multbin->GetSize(); i++)
 {
   low=Multbin->At(i);
   char cuts[100];
   char histtitle[100];
   char histtitle1[100];
   sprintf(cuts,"dNdEtaNBD>%i && dNdEtaNBD<%i",low,high);
   sprintf(histtitle,"B>>htemp%i(100)",50+i);
   sprintf(histtitle1,"htemp%i",50+i);
   ntuple->Draw(histtitle,cuts,"same");
   TH1F* htemp12 = (TH1F*)gPad->GetPrimitive(histtitle1);
   htemp12->SetLineColor(i+1);
   htemp12->DrawCopy("same");
   high=low;
 }

 TCanvas* canvas5 = new TCanvas("Eccentricity","Eccentricity");
 canvas5->SetLogy();
 ntuple->Draw("VarE>>htemp5(100)");
 TH1F* htemp5 = (TH1F*)gPad->GetPrimitive("htemp5");
 TH1F* hist5 = (TH1F*)htemp5->Clone();
 hist5->SetStats(0);
 hist5->SetTitle("");
 axis = hist5->GetXaxis();
 axis->SetTitle("standard eccentricity");
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 axis->SetTitleOffset(1.2);
 axis = hist5->GetYaxis();
 axis->SetTitle("");
 axis->SetTitleOffset(1.3);
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 hist5->DrawCopy();

 int low=0;
 int high=5000;
 for (int i=0; i<Multbin->GetSize(); i++)
 {
   low=Multbin->At(i);
   char cuts[100];
   char histtitle[100];
   char histtitle1[100];
   sprintf(cuts,"dNdEtaNBD>%i && dNdEtaNBD<%i",low,high);
   sprintf(histtitle,"VarE>>htemp%i(100)",60+i);
   sprintf(histtitle1,"htemp%i",60+i);
   ntuple->Draw(histtitle,cuts,"same");
   TH1F* htemp13 = (TH1F*)gPad->GetPrimitive(histtitle1);
   htemp13->SetLineColor(i+1);
   htemp13->DrawCopy("same");
   high=low;
 }

 TCanvas* canvas6 = new TCanvas("EccentricityPart","EccentricityPart");
 canvas6->SetLogy();
 ntuple->Draw("VarEPart>>htemp6(100)");
 TH1F* htemp6 = (TH1F*)gPad->GetPrimitive("htemp6");
 TH1F* hist6 = (TH1F*)htemp6->Clone();
 hist6->SetStats(0);
 hist6->SetTitle("");
 axis = hist6->GetXaxis();
 axis->SetTitle("participant eccentricity");
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 axis->SetTitleOffset(1.2);
 axis = hist6->GetYaxis();
 axis->SetTitle("");
 axis->SetTitleOffset(1.3);
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 hist6->DrawCopy();

 int low=0;
 int high=5000;
 for (int i=0; i<Multbin->GetSize(); i++)
 {
   low=Multbin->At(i);
   char cuts[100];
   char histtitle[100];
   char histtitle1[100];
   sprintf(cuts,"dNdEtaNBD>%i && dNdEtaNBD<%i",low,high);
   sprintf(histtitle,"VarEPart>>htemp%i(100)",70+i);
   sprintf(histtitle1,"htemp%i",70+i);
   ntuple->Draw(histtitle,cuts,"same");
   TH1F* htemp14 = (TH1F*)gPad->GetPrimitive(histtitle1);
   htemp14->SetLineColor(i+1);
   htemp14->DrawCopy("same");
   high=low;
 }

}


void getCentrality(TH1 *histNch, Float_t ff=1.0)
// histNch - histo of multiplicity distribution (better with binsize=1)
// ff fraction of accepted events. All losses are assumed to occur in most
// peripheral bin
{

 //double sum= histNch->GetEntries() - histNch->GetBinContent(1);
 double sum= histNch->Integral(); 
 int nbinsx=histNch->GetNbinsX();
 double frac=0.;
 int ic=0;
 for (int ib=nbinsx;ib>0;ib--){
   frac += histNch->GetBinContent(ib)/sum*100.*ff;
   if(frac > centPercent[ic]){
     binUp->SetAt(ib,ic);
     Multbin->SetAt(histNch->GetBinCenter(ib),ic);
     cout<<" centrality="<<centPercent[ic]<<"   mult <="<< histNch->GetBinCenter(ib) <<endl;
     ic++;
   }
   if(ic==nbins) break;
 }
 printf(" \n float multCent[%i] = {",nbins);
 // cout <<" \n float multCent[nbins] = {";

 for (int ic=nbins-1; ic>-1; ic--){
   cout<< histNch->GetBinCenter(binUp->At(ic));
   if (ic!=0) cout<<", ";
 }
 cout<<"};\n"<<endl;
}
