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

void plotCen2()
{
  // TFile f("GlauberMC_PbPb_ntuple.root");
 TFile f("dataNch.root");
 // TNtuple* ntuple = f.Get("nt_Pb_Pb");

 TCanvas* canvas = new TCanvas("Multiplicity","Multiplicity");
 canvas->SetLogy();
 // ntuple->Draw("dNdEtaNBD>>htemp(500)");
 // TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
 // TH1F* hist = (TH1F*)htemp->Clone();
 TH1F* hist = (TH1F*)hm2->Clone();
 double sum = hist->Integral();
 double width = hist->GetBinWidth(10);
 hist->Scale(1./sum/width);

 hist->SetTitle("");
 hist->SetStats(0);
 TAxis *axis = hist->GetXaxis();
 axis->SetTitle("Multiplicity");
 //axis->CenterTitle(kTRUE);
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 axis->SetTitleOffset(1.2);
 axis = hist->GetYaxis();
 axis->SetTitle("1/N_{ev} dN_{ev}/dn_{ch}");
 axis->SetTitleOffset(1.3);
 //axis->CenterTitle(kTRUE);
 axis->SetTitleSize(0.05);
 axis->SetLabelSize(0.04);
 hist->DrawCopy();
 

 //---------------------------------------------------
 // getCentrality(hist);
 //---------------------------------------------------

 int low=0;
 int high=hist->GetNbinsX();
 // float multCent[9] = {37.5, 82.5, 167.5, 297.5, 502.5, 782.5, 1157.5, 1692.5, 2057.5};
 float multCent[10] = {37.5, 82.5, 167.5, 297.5, 502.5, 782.5, 1157.5, 1692.5, 2057.5, 3000};
 int binCent[10];
 binCent[9]=high;
 int ic=0;
 for(int ib=1;ib<high;ib++){
   if(hist->GetBinCenter(ib) > multCent[ic]){
     binCent[ic]=ib;
     ic++;
   }
 }
 
 //printf("size: %i\n",binUp->GetSize());
 // for (int i=0; i<binUp->GetSize(); i++)
 for (int i=0; i<9; i++){
   hist->GetXaxis()->SetRange(binCent[i]+1,binCent[i+1]);
   hist->SetFillColor((i%2==0)?0:kGray);
   hist->DrawCopy("same");
 }  
 TLatex la;
  la.SetTextAlign(12);
  la.SetTextSize(0.045);
  la.SetIndiceSize(0.5);
  la.SetTextAngle(90);
  la.DrawLatex(2200,0.000005,"0 - 5%" );
  la.DrawLatex(1900,0.000005,"5% - 10%" );
  la.DrawLatex(1400,0.000005,"10% - 20%" );
  la.DrawLatex(950,0.000005,"20% - 30%" );
  la.DrawLatex(650,0.000005,"30% - 40%" );
  la.DrawLatex(400,0.000005,"40% - 50%" );
  la.DrawLatex(230,0.000005,"50% - 60%" );


 canvas->Print("mult_dist.png");
 canvas->Print("mult_dist.pdf");

}

