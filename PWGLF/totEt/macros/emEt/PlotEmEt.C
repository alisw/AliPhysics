
#ifndef __CINT__
#include "TGraphErrors.h"
#include "TTree.h"
#include "TList.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include <iostream>
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#endif

void PlotEmEt(TString filename = "Et.ESD.realPbPb.PHOS.root", TString det = "Emcal", bool sim = false, double scaleFactor = 1, double pionScale = 1.0)
{

    Double_t x[10], xerr[10],  y_rel_error[10] ,yerr[10], xpion[10], xpionerr[10], ypion[10], ypionerr[10], ysim[10], yrec[10];

    TCanvas *cents = new TCanvas("cents", "cents", 1920, 1000);
    cents->Divide(3,3);
    cents->cd(1);

    int n = 0;
    x[n] = xpion[n] = 15.8;
    ypion[n] = 0;
    xerr[n] = xpionerr[n] = 0.6;
    n++;

    x[n] = xpion[n] = 30.0;
    ypion[n] = 16.7/25./(xpion[n]/2)*scaleFactor*pionScale;
    xerr[n] = xpionerr[n] = 1.3;
    n++;

    x[n] = xpion[n] = 52.8;
    ypion[n] = 32.9/25./(xpion[n]/2)*scaleFactor*pionScale;
    xerr[n] = xpionerr[n] = 2;
    n++;

    x[n] = xpion[n] = 85;
    ypion[n] = 58.5/25./(xpion[n]/2)*scaleFactor*pionScale;
    xerr[n] = xpionerr[n] = 2.6;
    n++;

    x[n] = xpion[n] = 128.9;
    ypion[n] = 96.1/25./(xpion[n]/2)*scaleFactor*pionScale;
    xerr[n] = xpionerr[n] = 3.3;
    n++;

    x[n] = xpion[n] = 186.4;
    ypion[n] = 147.9/25./(xpion[n]/2)*scaleFactor*pionScale;
    xerr[n] = xpionerr[n] = 3.9;
    n++;

    x[n] = xpion[n] = 260.5;
    ypion[n] = 220.3/25./(xpion[n]/2)*scaleFactor*pionScale;
    xerr[n] = xpionerr[n] = 4.4;
    n++;

    xpion[n] = 329.7;
    ypion[n] = 293.3/25./(xpion[n]/2)*scaleFactor*pionScale;
    xpionerr[n] = 6;
    x[n] = 356;
    xerr[n] = 6;
    n++;

    xpion[n] = 382.8;
    xpionerr[n] = 6;
    ypion[n] = 362.1/25./(xpion[n]/2)*scaleFactor*pionScale;


    TFile *f = TFile::Open(filename, "READ");
    if (!f)
    {
        std::cerr << "Could not open file: " << filename << std::endl;
    }

    TList *l = dynamic_cast<TList*>(f->Get("out1"));

    if (!l)
    {
        std::cerr << "Could not get object list from: " << filename << std::endl;
    }
    TString treename = "fEventSummaryTree"+det+"Rec";
    TTree *t = dynamic_cast<TTree*>(l->FindObject(treename.Data()));

    if (!t)
    {
        std::cerr << "Could not get tree from: " << filename << std::endl;
	return;
    }
    TLegend *centLegends[10];
    for(int i = 0; i < n; i++)
    {
        TString cut = "fCentClass==" + TString::Format("%d", n-i-1);
        std::cout << cut << std::endl;
        cents->cd(i+1);
        cents->cd(i+1)->SetLogy(1);
        t->Draw("fTotNeutralEt", cut);
	t->GetHistogram()->SetStats(0);
	t->GetHistogram()->SetTitle(TString::Format("E_{T}, Centrality class: %d", n-i-1 ));
	cents->cd(i+1)->Update();
	centLegends[i] = new TLegend(0.55,0.7,0.9,0.8);
    
	centLegends[i]->SetFillColor(kWhite);
	centLegends[i]->AddEntry(t->GetHistogram(),"Reconstructed E_{T} distribution", "l");
        std::cout << t->GetHistogram()->GetMean() << " - " << x[i]*0.5 << std::endl;
	yrec[i] = t->GetHistogram()->GetMean()/(x[i]*0.5)*scaleFactor;
	yerr[i] = 0;//t->GetHistogram()->GetRMS();
	if(!sim)
	{
	  if(ypion[i])
	  {
	    y_rel_error[i] = TMath::Abs(yrec[i] - ypion[i])/ypion[i];
	  }
	  else
	  {
	    y_rel_error[i] = 0;
	  }
	}
	else
	{
	  y_rel_error[i] = t->GetHistogram()->GetMean()/(x[i]*0.5);
	}

    }

    if(sim)
    {
        t = dynamic_cast<TTree*>(l->FindObject("fEventSummaryTreePhosMC"));

        if (!t)
        {
            std::cerr << "Could not get sim tree from: " << filename << std::endl;
	    return;
        }

        std::cout << n << std::endl;
        for(int i = 0; i < n; i++)
        {
            TString cut = "fCentClass==" + TString::Format("%d", n-i-1);
	    cents->cd(i+1);
	    t->SetLineColor(kRed);
	    
            t->Draw("fTotNeutralEt", cut,"SAME");
	    centLegends[i]->AddEntry(t->GetHistogram(),"True E_{T} distribution", "l");
	    centLegends[i]->Draw();
	    ysim[i] = t->GetHistogram()->GetMean()/(x[i]*0.5)*scaleFactor;

            yerr[i] = 0;//t->GetHistogram()->GetRMS();

	    y_rel_error[i] = ( (y_rel_error[i]-(t->GetHistogram()->GetMean()/(x[i]*0.5)) )/(t->GetHistogram()->GetMean()/(x[i]*0.5)));
        }
    }
    if(0)
    {
        y_rel_error[0] = 1.9*scaleFactor;
        yerr[0] = 0.0;

        y_rel_error[1] = 2.3*scaleFactor;
        yerr[1] = 0;

        y_rel_error[2] = 2.7*scaleFactor;
        yerr[2] = 0;

        y_rel_error[3] = 3.1*scaleFactor;
        yerr[3] = 0;

        y_rel_error[4] = 3.5*scaleFactor;
        yerr[4] = 0.0;

        y_rel_error[5] = 3.7*scaleFactor;
        yerr[5] = 0.0;

        y_rel_error[6] = 4.07*scaleFactor;
        yerr[6] =0;

        y_rel_error[7] = 4.61*scaleFactor;
        yerr[7] = 0;

    }

    TCanvas *c2 = new TCanvas("c2", "dE_{T}/d#eta#frac{1}{0.5*N_{part}} [GeV]",1920, 1000);
   
    TGraphErrors *gr = new TGraphErrors(n-2,x+2,yrec+2,xerr,yerr);
    //TGraph *gr = new TGraph(n,x2,y2);
    gr->SetTitle("");
    gr->GetYaxis()->SetTitle("dE_{T}/d#eta#frac{1}{0.5*N_{part}} [GeV]");
    //gr->GetYaxis()->SetTitle("Relative error in reconstructed vs true EM E_{T} ");

    gr->GetXaxis()->SetTitle("N_{part}");

    gr->GetXaxis()->SetRangeUser(0, 400);

    gr->GetYaxis()->SetRangeUser(0, 0.15*scaleFactor);//0.15*scaleFactor);
    gr->SetMarkerStyle(20);
    gr->Draw("AP");

    TLegend *leg = new TLegend(0.1,0.7,0.55,0.9);
    leg->SetFillColor(kWhite);
    leg->AddEntry(gr,"Reconstructed EM E_{T}", "lp");
    if(!sim) 
    {
      TGraphErrors *gr2 = new TGraphErrors(n+1-2,xpion+2,ypion+2,xpionerr+2,ypionerr+2);
      gr2->SetTitle("");

      gr2->GetXaxis()->SetTitle("N_{part}");

      gr2->GetXaxis()->SetRangeUser(0, 400);
      gr2->GetYaxis()->SetRangeUser(0,0.15*scaleFactor);
      gr2->SetMarkerStyle(20);
      gr2->SetMarkerColor(kRed);
      gr2->Draw("P");
      leg->AddEntry(gr2,"Reconstructed charged pion E_{T} (scaled by 0.5)", "lp");
    }
    else
    {
      TGraphErrors *gr2 = new TGraphErrors(n+1-2,x+2,ysim+2,xerr+2,yerr+2);
      gr2->SetTitle("");

      gr2->GetXaxis()->SetTitle("N_{part}");

      gr2->GetYaxis()->SetRangeUser(0,0.15*scaleFactor);
      gr2->SetMarkerStyle(20);
      gr2->SetMarkerColor(kRed);
      gr2->Draw("P");
      leg->AddEntry(gr2,"True E_{T}", "lp");
    }
      
    
    leg->Draw();
    
    
    TString title;
    TString ytitle;
    TString xtitle = "N_{part}";
    if(sim)
    {
	title = "Relative Discrepancy for reconstructed vs true";
	ytitle = "#frac{E_{T,rec}-E_{T,true}}{E_{T,true}}";
    }
    else
    {
      title = "Relative Discrepancy for EM E_{T} vs EM #pi^{+/-}";
      ytitle = "#frac{E_{T,rec}-E_{T,#pi^{+/-}}}{E_{T,#pi^{+/-}}}";
    }
     
    TCanvas *c3 = new TCanvas("c3", title,1920, 1000);
    TGraphErrors *gr3 = 0;
    if(sim)
    {
      gr3 = new TGraphErrors(n-2,x+2,y_rel_error+2,xerr+2,yerr+2);
    }
    else
    {
      gr3 = new TGraphErrors(n-3,x+2,y_rel_error+2,xerr+2,yerr+2);
    }
    gr3->SetTitle(title);
    gr3->GetYaxis()->SetTitle(ytitle);
    gr3->GetXaxis()->SetTitle(xtitle);
    gr3->GetXaxis()->SetRangeUser(0, 400);
    if(sim) gr3->GetYaxis()->SetRangeUser(-0.2,0.2);
    else gr3->GetYaxis()->SetRangeUser(-0.3,0.3);
    gr3->SetMarkerStyle(20);
    gr3->Draw("AP");
    
    TString name = "dn_deta_npart_scale_factor_" + TString::Format("%.2f", scaleFactor);
    if(sim) name += "_sim.png";
    else name += "_real.png";
    c2->Print(name, "png");

    name = "rel_error_";
    if(sim) name += "_sim.png";
    else name += "_real.png";
    c3->Print(name, "png");

    if(sim) name = "tot_neutral_et_sim.png";
    else name = "tot_neutral_et_real.png";
    cents->Print(name, "png");



}
