#include <TStyle.h>
#include <TColor.h>
#include <TLegend.h>
#include <iostream>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TH1D.h>
#include <complex>
#include <TComplex.h>
#include <iostream>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TBox.h>
#include <TMultiGraph.h>
#include <Riostream.h>
#include <TPad.h>

using namespace std;

Double_t threshold = 0.139570*2.;
double PhaseSpaceFactor(double x, double T, double pT){
 return ( x / ( sqrt( pow(x,2) + pow(pT,2) ) ) ) * exp( -sqrt( pow(x,2) + pow(pT,2) )/T );
}

double Width(double x, double M0, double Gam0, int Spin){
 double w = pow( (x*x - threshold*threshold)/(M0*M0 - threshold*threshold),0.5+(double)Spin )*Gam0*M0/x;
 return w;
}

double breitWigner(double x, double Amp, double M0, double Gam0, int Spin){
 double br = Amp*x*M0*Width(x, M0, Gam0, Spin);
 br /= ( pow(M0*M0-x*x,2) + M0*M0*pow(Width(x, M0, Gam0, Spin),2) );
 return br;
}
double background(double x, double ind, double b1, double b2){
 double bg = pow(x-threshold,ind)*exp(b1*x + b2*x*x);
 return bg;
}
double breitWigner_f0(double *x, double* par){
 return breitWigner(x[0], par[1], par[0], par[2], 0);
}
double sfit(double *x, double* par){

 double bgA = par[0];
 double bgind = par[1];
 double bgb1 = par[2];
 double bgb2 = par[3];

 double f0m = par[4];
 double f0A = par[5];
 double f0g = par[6];
 int f0s = 0;

 double f2m = par[7];
 double f2A = par[8];
 double f2g = par[9];
 int f2s = 2;

 double rhom = par[10];
 double rhoA = par[11];
 double rhog = par[12];
 int rhos = 1;

 double T = par[13];
 double pT = par[14];

 double breitWigner_f0  = breitWigner(x[0],f0A,f0m,f0g,f0s);
 double breitWigner_f2  = breitWigner(x[0],f2A,f2m,f2g,f2s);
 double breitWigner_rho  = breitWigner(x[0],rhoA,rhom,rhog,rhos);

 double bgfunc = bgA*background(x[0],bgind,bgb1,bgb2);

// return ( breitWigner_f0 + breitWigner_f2 + breitWigner_rho )*PhaseSpaceFactor(x[0],T,pT) + bgfunc;
 return breitWigner_f0 + breitWigner_f2 + breitWigner_rho + bgfunc;
}
double LevyTsallisRho(double *x, double *par){
 double mass = 0.750;
 return x[0] * par[0] *
         ( ( par[1]-1 )*( par[1]-2 ) ) /
         ( par[1]*par[2] * (par[1]*par[2] + mass*( par[1]-2 ) ) ) *
         pow( ( 1 + ( sqrt( mass*mass + x[0]*x[0] ) - mass ) / par[1] / par[2] ), -par[1] );
}
double LevyTsallisF0(double *x, double *par){
 double mass = 0.980;
 return x[0] * par[0] *
         ( ( par[1]-1 )*( par[1]-2 ) ) /
         ( par[1]*par[2] * (par[1]*par[2] + mass*( par[1]-2 ) ) ) *
         pow( ( 1 + ( sqrt( mass*mass + x[0]*x[0] ) - mass ) / par[1] / par[2] ), -par[1] );
}
double LevyTsallisF2(double *x, double *par){
 double mass = 1.2755;
 return x[0] * par[0] *
         ( ( par[1]-1 )*( par[1]-2 ) ) /
         ( par[1]*par[2] * (par[1]*par[2] + mass*( par[1]-2 ) ) ) *
         pow( ( 1 + ( sqrt( mass*mass + x[0]*x[0] ) - mass ) / par[1] / par[2] ), -par[1] );
}
double LevyTsallisRho1stMom(double *x, double *par){
 double mass = 0.750;
 return x[0] * par[0] *
         ( ( par[1]-1 )*( par[1]-2 ) ) /
         ( par[1]*par[2] * (par[1]*par[2] + mass*( par[1]-2 ) ) ) *
         pow( ( 1 + ( sqrt( mass*mass + x[0]*x[0] ) - mass ) / par[1] / par[2] ), -par[1] ) *
         x[0];
}
double LevyTsallisF01stMom(double *x, double *par){
 double mass = 0.980;
 return x[0] * par[0] *
         ( ( par[1]-1 )*( par[1]-2 ) ) /
         ( par[1]*par[2] * (par[1]*par[2] + mass*( par[1]-2 ) ) ) *
         pow( ( 1 + ( sqrt( mass*mass + x[0]*x[0] ) - mass ) / par[1] / par[2] ), -par[1] ) *
         x[0];

}
double LevyTsallisF21stMom(double *x, double *par){
 double mass = 1.2755;
 return x[0] * par[0] *
         ( ( par[1]-1 )*( par[1]-2 ) ) /
         ( par[1]*par[2] * (par[1]*par[2] + mass*( par[1]-2 ) ) ) *
         pow( ( 1 + ( sqrt( mass*mass + x[0]*x[0] ) - mass ) / par[1] / par[2] ), -par[1] ) *
         x[0];
}
void GetSignals_13TeV(){
 gStyle->SetTitleFont(22,"y");
 gStyle->SetTitleFont(22,"z");
 gStyle->SetLabelFont(22,"x");
 gStyle->SetLabelFont(22,"y");
 gStyle->SetLabelFont(22,"z");

// TFile* fin = new TFile("./ResultFiles/AnalysisResults_f0f2LHCAODSysZSysTrkSysPID.root","read");
 TFile* fin = new TFile("./ResultFiles/AnalysisResults.root","read");
 TFile* fmc = new TFile("./ResultFiles/AnalysisResults_f0f2_Correction_13TeV.root","read");
 TFile* fmcTrig = new TFile("./ResultFiles/AnalysisResults_f0f2_TrigEff_13TeV.root","read");

 auto output = fin->Get("output");
 auto output_mc = fmc->Get("output");
 auto output_trig = fmcTrig->Get("output");

 THnSparse* evtSel = (THnSparse*)output->FindObject("EvtSelector");
 THnSparse* SP_trig = (THnSparse*)output_trig->FindObject("TrigEffMult");
 THnSparse* hs = (THnSparse*)output->FindObject("hInvMass");

 const int nbins_mult = 6;
 double multmin[nbins_mult] = {
	0, 5, 10, 20, 50, 0 };
 double multmax[nbins_mult] = {
	5, 10, 20, 50, 100, 100 };

 const int nbins_pt = 8;
 
 double ptmin[nbins_pt] = {
//         0.0, 0.3, 0.6, 1.0, 1.5,
//         2.0, 2.5, 3.0, 3.5, 4.0,
//         5.0, 6.0, 7.0, 8.0, 0.0};
	0.3, 1.0, 2.0, 3.0, 4.0,
	5.0, 6.0, 0.3 };
 double ptmax[nbins_pt] = {
//	 0.3, 0.6, 1.0, 1.5,
//         2.0, 2.5, 3.0, 3.5, 4.0,
//         5.0, 6.0, 7.0, 8.0, 13.0, 13.0};
	1.0, 2.0, 3.0, 4.0, 5.0,
	6.0, 8.0, 8.0 };

/*
 double ptmin[nbins_pt] = {
	0.0, 0.6, 1.5, 2.5, 3.5,
	5.0, 7.0, 9.0, 0.0};
 double ptmax[nbins_pt] = {
	0.6, 1.5, 2.5, 3.5, 5.0,
	7.0, 9.0, 13.0, 13.0 };
*/
 //event class definition


 TH1D* hTrigEff[nbins_mult];
 double TrigEff[nbins_mult];
 for(int i=0;i<nbins_mult;i++){
        SP_trig->GetAxis(0)->SetRangeUser( multmin[i], multmax[i] );
        hTrigEff[i] = (TH1D*)SP_trig->Projection(1);
        TrigEff[i] = hTrigEff[i]->GetMean();
 }
//Trigger efficiency

 TH1D* hMult;
 double MultEvt[nbins_mult];
 evtSel->GetAxis(0)->SetRangeUser( -10, 10 );
 hMult = (TH1D*)evtSel->Projection(1);
 for(int l=0;l<nbins_mult;l++){
        MultEvt[l] = 0.0;
        for(int i=0;i<hMult->GetNbinsX();i++){
                if( hMult->GetBinCenter(i+1) > multmin[l] &&
                    hMult->GetBinCenter(i+1) < multmax[l] ){
                        MultEvt[l] += hMult->GetBinContent(i+1);
                }
        }
 }
//event number

 TH1D* hMassNP[nbins_mult][nbins_pt];
 TH1D* hMassPP[nbins_mult][nbins_pt];
 TH1D* hMassNN[nbins_mult][nbins_pt];
 TH1D* hMassLS[nbins_mult][nbins_pt];
 TH1D* hMassSigLS[nbins_mult][nbins_pt];

 TLegend* CosmeticLegend_invmass = new TLegend(0.6,0.7,0.85,0.89);
 CosmeticLegend_invmass->SetLineWidth(0.0);
 CosmeticLegend_invmass->SetFillColorAlpha(0,0);

 TCanvas* ca = new TCanvas("ca","ca",1200,800);
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(0); 
 hs->GetAxis(5)->SetRangeUser( 1, 1);
 hs->GetAxis(1)->SetRangeUser( -10, 10);
 for(int i=0;i<nbins_mult;i++){
	hs->GetAxis(2)->SetRangeUser( multmin[i], multmax[i] );
	for(int j=0;j<nbins_pt;j++){
		hs->GetAxis(3)->SetRangeUser( ptmin[j], ptmax[j]-0.001 );
		hs->GetAxis(0)->SetRange(1,1);
		hMassNP[i][j] = (TH1D*)hs->Projection(4);
		hMassNP[i][j]->SetTitle(Form(" +-, %.0lf-%.0lf %%, %.1lf < #font[12]{p}_{#font[22]{T}} < %.1lf",multmin[i],multmax[i],ptmin[j],ptmax[j] ));
		hMassNP[i][j]->GetXaxis()->SetTitle("M_{#pi#pi} (GeV/#font[12]{c^{2}})");


		hs->GetAxis(0)->SetRange(2,2);
		hMassNN[i][j] = (TH1D*)hs->Projection(4);

		hs->GetAxis(0)->SetRange(3,3);
		hMassPP[i][j] = (TH1D*)hs->Projection(4);

		hMassLS[i][j] = (TH1D*)hMassNN[i][j]->Clone();	
		for(int k=0;k<hMassNN[i][j]->GetNbinsX();k++){
			hMassLS[i][j]->SetBinContent( k+1,
				2.0*sqrt( hMassNN[i][j]->GetBinContent(k+1)*hMassPP[i][j]->GetBinContent(k+1) ) );
			hMassLS[i][j]->SetBinError( k+1,
				sqrt( hMassNN[i][j]->GetBinContent(k+1)*hMassPP[i][j]->GetBinContent(k+1) )*
				sqrt( 1.0/hMassNN[i][j]->GetBinContent(k+1)+1.0/hMassPP[i][j]->GetBinContent(k+1) ) );
		}

		ca->cd();

		hMassNP[i][j]->SetMinimum(0.001);
		hMassNP[i][j]->GetXaxis()->SetRangeUser(0.3,1.8);
		hMassNP[i][j]->SetMarkerStyle(20);
		hMassNP[i][j]->SetMarkerColor(1);
		hMassNP[i][j]->SetLineColor(1);
		hMassNP[i][j]->Draw();

		hMassLS[i][j]->SetMarkerStyle(21);
		hMassLS[i][j]->SetMarkerColor(2);
		hMassLS[i][j]->SetLineColor(2);
		hMassLS[i][j]->Draw("same");

                CosmeticLegend_invmass->SetHeader("ALICE, pp, #sqrt{s} = 13 TeV","C");
                CosmeticLegend_invmass->AddEntry((TObject*)0,Form("%.0lf-%.0lf %%, %.1lf < #font[12]{p}_{#font[22]{T}} < %.1lf (GeV/#font[12]{c})",multmin[i],multmax[i],ptmin[j],ptmax[j]),"" );
                CosmeticLegend_invmass->AddEntry((TObject*)0,"#font[22]{f}_{0}(980) #rightarrow #pi^{+}#pi^{-} (|y|<0.5)","");
                CosmeticLegend_invmass->AddEntry(hMassNP[i][j],"Unlike-Sign Pairs","P");
                CosmeticLegend_invmass->AddEntry(hMassLS[i][j],"Like-Sign Pairs","P");
                CosmeticLegend_invmass->Draw();

		ca->SaveAs(Form("./GetSignals_13TeV/LS_%d_%d.pdf",i,j));
		CosmeticLegend_invmass->Clear();

		hMassSigLS[i][j] = (TH1D*)hMassNP[i][j]->Clone();
		hMassSigLS[i][j]->Add( hMassLS[i][j], -1.0 );
		hMassSigLS[i][j]->SetTitle(Form(" (+-) - 2#sqrt{++#times--} Mult[%.1f,%.1f], pT[%.1f,%.1f]",multmin[i],multmax[i],ptmin[j],ptmax[j]));
		hMassSigLS[i][j]->GetYaxis()->SetTitle("Number of Counts / 5 (MeV/#font[12]{c^{2}})");
		hMassSigLS[i][j]->GetXaxis()->SetRangeUser(0.3,1.8);
		hMassSigLS[i][j]->SetMarkerStyle(24);
		hMassSigLS[i][j]->SetMarkerColor(kBlack);
		hMassSigLS[i][j]->SetLineColor(kBlack);

	
	}
 }
//Like-sign Subtraction

 const int Npar = 13;
 const int Ncomp = 4;
 int MagIndex[Ncomp] = {0,5,8,11};
 char CompName[Ncomp][100]={"BG","F0","F2","#rho"};

 double *FitPar[nbins_mult][nbins_pt];
 double *FitParErr[nbins_mult][nbins_pt];
 double ReducedChi2[nbins_mult][nbins_pt];

 char ParName[Npar][1000]={
        "BGAmp","BGSlope","BG1st","BG2nd",
        "f0Mass","f0Amp","f0W",
        "f2Mass","f2Amp","f2W",
        "rhoMass","rhoAmp","rhoW" };

 TF1* fitBin[nbins_mult][nbins_pt];
 TF1* fitBinComp[nbins_mult][nbins_pt][Ncomp];
 TFitResultPtr fitResults[nbins_mult][nbins_pt];

 TF1* fbreitWigner_f0[nbins_mult][nbins_pt];

 for(int i=0;i<nbins_mult;i++){
	for(int j=0;j<nbins_pt;j++){
		fitBin[i][j] = new TF1("f1",sfit,0.7,1.7,Npar);
		fbreitWigner_f0[i][j] = new TF1("f1",breitWigner_f0,0.7,1.7,3);
		for(int k=0;k<Ncomp;k++){
			fitBinComp[i][j][k] = new TF1("f1",sfit,0.7,1.7,Npar);
		}
	}
 }

 double BGAmpMin = 1e2;
// double BGAmpMax = 1.1e7 * MultEvt[nbins_mult-1] / 4.0e8;
 double BGAmpMax = 2e7;

 double BGSlopeMin = -1.5;
 double BGSlopeMax = 2.5;

 double BGInd1Min = -5.0;
 double BGInd1Max = -1.0;

 double BGInd2Min = -1.5;
 double BGInd2Max = 1.5;

 double f0MassMinRange = 0.95;
 double f0MassMaxRange = 0.985;

 double f0AmpMin = 0.1;
 double f0AmpMax = 2e4;

 double f0WMin = 0.025;
 double f0WMax = 0.1;

 double f2MassMinRange = 1.2751-0.0012*5.0;
 double f2MassMaxRange = 1.2751+0.0012*5.0;

 double f2AmpMin = 0.01;
 double f2AmpMax = 2e4;

 double f2WMin = 0.1851-0.0024*0.02;
 double f2WMax = 0.1851+0.0029*0.02;

 double rhoMassMin = 0.74;
 double rhoMassMax = 0.75;

 double rhoAmpMin = 1.0;
 double rhoAmpMax = 1e7;

 double rhoWMin = 0.1462-0.001;
 double rhoWMax = 0.1462+0.001;


 double FitParMin[Npar] = {
        BGAmpMin, BGSlopeMin, BGInd1Min, BGInd2Min, f0MassMinRange,
        f0AmpMin, f0WMin, f2MassMinRange, f2AmpMin, f2WMin,
        rhoMassMin, rhoAmpMin, rhoWMin };

 double FitParMax[Npar] = {
        BGAmpMax, BGSlopeMax, BGInd1Max, BGInd2Max, f0MassMaxRange,
        f0AmpMax, f0WMax, f2MassMaxRange, f2AmpMax, f2WMax,
        rhoMassMax, rhoAmpMax, rhoWMax };


 TCanvas* cb = new TCanvas("cb","cb",1300,800);
 cb->DrawFrame(0,0,1,1);
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(0);
 cb->SetRightMargin(0.4);

 TLatex* latex_rho = new TLatex();
 TLatex* latex_f0 = new TLatex();
 TLatex* latex_f2 = new TLatex();
 TLatex* latex_evt = new TLatex();
 TLatex* latex_bg = new TLatex();

 latex_rho->SetTextColor(6);
 latex_f0->SetTextColor(3);
 latex_f2->SetTextColor(4);

 latex_rho->SetTextFont(22);
 latex_f0->SetTextFont(22);
 latex_f2->SetTextFont(22);
 latex_evt->SetTextFont(22);
 latex_bg->SetTextFont(22);

 latex_evt->SetTextSize(0.03);
 latex_rho->SetTextSize(0.03);
 latex_f0->SetTextSize(0.03);
 latex_f2->SetTextSize(0.03);
 latex_bg->SetTextSize(0.03);

 for(int i=0;i<nbins_mult;i++){
        for(int j=0;j<nbins_pt;j++){

  		fitBin[i][j]->SetParLimits(0,hMassSigLS[i][j]->GetMaximum()*1.0,BGAmpMax);
                if( j>4 && (i!=3) ) fitBin[i][j]->SetParLimits(1,0.2,5.0);
                if( j==6 && i==1 ) fitBin[i][j]->SetParLimits(1,1.5,2.5);
                if( j==6 && i==2 ) fitBin[i][j]->SetParLimits(1,1.0,1.2);
                if( i==3 && j >3 ) fitBin[i][j]->SetParLimits(1,1.2,2.0);
                if( i==5 && j==5 ) fitBin[i][j]->SetParLimits(1,1.0,1.1);
                if( i==2 && j==3 ) fitBin[i][j]->SetParLimits(1,2.3,3.5);
                fitBin[i][j]->SetParLimits(2,BGInd1Min,BGInd1Max);
                fitBin[i][j]->SetParLimits(3,BGInd2Min,BGInd2Max);
                if( i==2 && j==6 ) fitBin[i][j]->SetParLimits(3,0.3,1.0);
                if( i==2 && j==3 ) fitBin[i][j]->SetParLimits(3,-0.3,0.0);

		fitBin[i][j]->SetParLimits(4,0.965,0.985);
		fitBin[i][j]->SetParLimits(5,f0AmpMin,f0AmpMax);
		fitBin[i][j]->FixParameter(6, 0.040 );

		fitBin[i][j]->SetParLimits(7,f2MassMinRange,f2MassMaxRange);
		fitBin[i][j]->SetParLimits(8,f2AmpMin,f2AmpMax);
                fitBin[i][j]->SetParLimits(9,f2WMin,f2WMax);

                fitBin[i][j]->SetParLimits(10,rhoMassMin,rhoMassMax);
                fitBin[i][j]->SetParLimits(11,rhoAmpMin,rhoAmpMax);
                fitBin[i][j]->SetParLimits(12,rhoWMin,rhoWMax);

		fitResults[i][j] = (TFitResultPtr)hMassSigLS[i][j]->Fit(fitBin[i][j],"s","",0.75,1.7);


                fitBin[i][j]->SetParLimits(0,fitBin[i][j]->GetParameter(0)*0.9, fitBin[i][j]->GetParameter(0)*1.1 );
                fitBin[i][j]->FixParameter(1, fitBin[i][j]->GetParameter(1) );
                fitBin[i][j]->FixParameter(2, fitBin[i][j]->GetParameter(2) );
                fitBin[i][j]->FixParameter(3, fitBin[i][j]->GetParameter(3) );

		fitBin[i][j]->SetParLimits(4,f0MassMinRange,f0MassMaxRange);
		fitBin[i][j]->SetParLimits(5,f0AmpMin,f0AmpMax);
		fitBin[i][j]->SetParLimits(6,f0WMin,f0WMax);

		fitBin[i][j]->SetParLimits(7,f2MassMinRange,f2MassMaxRange);
		fitBin[i][j]->SetParLimits(8,f2AmpMin,f2AmpMax);
		fitBin[i][j]->SetParLimits(9,f2WMin,f2WMax);

		fitBin[i][j]->FixParameter(10, fitBin[i][j]->GetParameter(10) );
		fitBin[i][j]->SetParLimits(11,rhoAmpMin,rhoAmpMax);
		fitBin[i][j]->FixParameter(12, fitBin[i][j]->GetParameter(12) );

                fitResults[i][j] = (TFitResultPtr)hMassSigLS[i][j]->Fit(fitBin[i][j],"s","",0.8,1.7);

		cb->cd();                
		cb->SetRightMargin(0.4);

		hMassSigLS[i][j]->GetXaxis()->SetRangeUser(0.7,1.76);
		hMassSigLS[i][j]->GetFunction("f1")->SetLineColor(kBlack);
		hMassSigLS[i][j]->SetMinimum(0.01);
		hMassSigLS[i][j]->Draw();

		FitPar[i][j] = (double*)fitBin[i][j]->GetParameters();
                FitParErr[i][j] = (double*)fitBin[i][j]->GetParErrors();

                if( fitBin[i][j]->GetNDF() != 0 )
                ReducedChi2[i][j] = fitBin[i][j]->GetChisquare() / fitBin[i][j]->GetNDF();
                else{  ReducedChi2[i][j]= 1e4; }

                for(int k=0;k<Ncomp;k++){
                	fitBinComp[i][j][k]->SetParameters( fitBin[i][j]->GetParameters() );
                        for(int l=0;l<Ncomp;l++){
                        	if( k==l ) continue;
                                fitBinComp[i][j][k]->SetParameter( MagIndex[l], 0.0 );
                        }
                        fitBinComp[i][j][k]->SetLineColor(k+2 + (k+2)/5);
                        fitBinComp[i][j][k]->Draw("same");
		}


		fbreitWigner_f0[i][j]->SetParameter(0, fitBin[i][j]->GetParameter(4) );
		fbreitWigner_f0[i][j]->SetParameter(1, fitBin[i][j]->GetParameter(5) );
		fbreitWigner_f0[i][j]->SetParameter(2, fitBin[i][j]->GetParameter(6) );

		fbreitWigner_f0[i][j]->SetParError(0, fitBin[i][j]->GetParError(4) );
		fbreitWigner_f0[i][j]->SetParError(1, fitBin[i][j]->GetParError(5) );
		fbreitWigner_f0[i][j]->SetParError(2, fitBin[i][j]->GetParError(6) );



                cb->cd();
                latex_evt->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*1.05,
                	Form("#splitline{Event Class}{%.1lf < #font[12]{p}_{#font[22]{T}}/(GeV/#font[12]{c}) <%.1lf,\t %.0lf-%.0lf %%}",
                        ptmin[j],ptmax[j],multmin[i],multmax[i]) );
                latex_evt->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.98,
                        Form("#chi^{2} / NDF = %.2lf/%d = %.2lf",fitBin[i][j]->GetChisquare(),fitBin[i][j]->GetNDF(),fitBin[i][j]->GetChisquare()/fitBin[i][j]->GetNDF()));

                latex_rho->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.92,
                	Form("#rho^{0} Mass : %.3lf #pm %.3lf (GeV/#font[12]{c^{2}})",fitBin[i][j]->GetParameter(10),fitBin[i][j]->GetParError(10)));
                latex_rho->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.87,
                        Form("#rho^{0} Width : %.3lf #pm %.3lf (GeV/#font[12]{c^{2}}) (constrained)",fitBin[i][j]->GetParameter(12),fitBin[i][j]->GetParError(12)));
                latex_rho->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.82,
                        Form("#rho^{0} Amplitude : %.3lf #pm %.3lf ",fitBin[i][j]->GetParameter(11),fitBin[i][j]->GetParError(11)));

                latex_f0->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.77,
                	Form("#font[22]{f}_{0} Mass : %.3lf #pm %.3lf (GeV/#font[12]{c^{2}})",fitBin[i][j]->GetParameter(4),fitBin[i][j]->GetParError(4)));
                latex_f0->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.72,
                        Form("#font[22]{f}_{0} Width : %.3lf #pm %.3lf (GeV/#font[12]{c^{2}}) (constrained)",fitBin[i][j]->GetParameter(6),fitBin[i][j]->GetParError(6)));
                latex_f0->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.67,
                        Form("#font[22]{f}_{0} Amplitude : %.3lf #pm %.3lf ",fitBin[i][j]->GetParameter(5),fitBin[i][j]->GetParError(5)));

                latex_f2->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.62,
                        Form("#font[22]{f}_{2} Mass : %.3lf #pm %.3lf (GeV/#font[12]{c^{2}}) (constrained)",fitBin[i][j]->GetParameter(7),fitBin[i][j]->GetParError(7)));
                latex_f2->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.57,
                        Form("#font[22]{f}_{2} Width : %.3lf #pm %.3lf (GeV/#font[12]{c^{2}}) (constrained)",fitBin[i][j]->GetParameter(9),fitBin[i][j]->GetParError(9)));
                latex_f2->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.52,
                        Form("#font[22]{f}_{2} Amplitude : %.3lf #pm %.3lf ",fitBin[i][j]->GetParameter(8),fitBin[i][j]->GetParError(8)));

                latex_bg->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.42,
                        Form("#splitline{#color[2]{Bkg. Par.} #color[1]{(}#color[2]{A}#color[1]{ (m-2m_{#pi})^{slp}#font[22]{exp}(c_{1}m+c_{2}m^{2}) } )}{= %.3lf (fixed)}",fitBin[i][j]->GetParameter(0) ));
                latex_bg->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.33,
                        Form("#splitline{#color[2]{Bkg. Par.} #color[1]{(A}#color[1]{ (m-2m_{#pi})^{#color[2]{slp}}#font[22]{exp}(c_{1}m+c_{2}m^{2}) } )}{= %.3lf (fixed)}",fitBin[i][j]->GetParameter(1) ));
                latex_bg->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.24,
                        Form("#splitline{#color[2]{Bkg. Par.} #color[1]{(A}#color[1]{ (m-2m_{#pi})^{slp}#font[22]{exp}(#color[2]{c_{1}}m+c_{2}m^{2}) } )}{= %.3lf (fixed)}",fitBin[i][j]->GetParameter(2) ));
                latex_bg->DrawLatex(1.8,hMassSigLS[i][j]->GetMaximum()*0.15,
                        Form("#splitline{#color[2]{Bkg. Par.} #color[1]{(A}#color[1]{ (m-2m_{#pi})^{slp}#font[22]{exp}(c_{1}m+#color[2]{c_{2}}m^{2}) } )}{= %.3lf (fixed)}",fitBin[i][j]->GetParameter(3) ));


                cb->SaveAs(Form("./GetSignals_13TeV/FitResults_%d_%d.pdf",i,j));
	}
 }
//Fit procedure

 TGraphErrors* gFitPar[Npar][nbins_mult];
 TGraphErrors* gFitParErr[Npar][nbins_mult];
 TGraphErrors* gChi2[nbins_mult];

 double ParX[nbins_pt];
 double ParY[Npar][nbins_pt];
 double ParErrX[nbins_pt];
 double ParErrY[Npar][nbins_pt];

 int badcount[Npar];
 double BadParX[Npar][nbins_pt];
 double BadParY[Npar][nbins_pt];
 double BadParErrX[Npar][nbins_pt];
 double BadParErrY[Npar][nbins_pt];

 for(int i=0;i<nbins_mult;i++){
	for(int j=0;j<Npar;j++){
		badcount[j] = 0;
	}
	for(int j=0;j<nbins_pt-1;j++){
		ParX[j] = (ptmin[j]+ptmax[j])/2.0;
		ParErrX[j] = (ptmax[j]-ptmin[j])/2.0;

		for(int k=0;k<Npar;k++){
			ParY[k][j] = FitPar[i][j][k];
			ParErrY[k][j] = FitParErr[i][j][k];
			if( ReducedChi2[i][j] > 3.0 ){
				BadParX[k][badcount[k]] = ParX[j];
				BadParErrX[k][badcount[k]] = ParErrX[j];
				BadParY[k][badcount[k]] = ParY[k][j];
				BadParErrY[k][badcount[k]] = ParErrY[k][j];
				badcount[k]++;
			}
		}
	}
	for(int k=0;k<Npar;k++){
		gFitPar[k][i] = new TGraphErrors( nbins_pt-1, ParX, ParY[k], ParErrX, ParErrY[k] );
		gFitParErr[k][i] = new TGraphErrors( badcount[k], BadParX[k], BadParY[k], BadParErrX[k], BadParErrY[k] );

		gFitPar[k][i]->SetMaximum( FitParMax[k]*1.05 );
		gFitPar[k][i]->SetMinimum( FitParMin[k]*0.95 );
	}

	for(int k=0;k<Npar;k++){
		ca->cd();
		if(k==6){
//			gStyle->SetOptFit(1);
//			gFitPar[k][i]->Fit("pol0");
		}
		gFitPar[k][i]->SetMarkerStyle(20);
		gFitPar[k][i]->SetMarkerColor(1);
		gFitPar[k][i]->SetLineColor(1);
		gFitPar[k][i]->Draw("AP");

		gFitParErr[k][i]->SetMarkerStyle(22);
		gFitParErr[k][i]->SetMarkerColor(2);
		gFitParErr[k][i]->SetLineColor(2);
		gFitParErr[k][i]->SetMarkerSize(2.0);
		gFitParErr[k][i]->Draw("P");

		ca->SaveAs(Form("./GetSignals_13TeV/%s_%d.pdf",ParName[k],i) );
	} 
 }

 const int nparticle = 3;
 TH1D* pTGenSpectra[nparticle][nbins_mult];
 TH1D* pTAcpSpectra[nparticle][nbins_mult];
 TH1D* pTCorrectionSpectra[nparticle][nbins_mult];
 THnSparse* hMCgen;
 THnSparse* hMCacp;

 char MCGenName[nparticle][1000] = {
         "hRhoGenParticle",
         "hF0GenParticle",
         "hF2GenParticle" };
 char MCAcpName[nparticle][1000] = {
         "hRhoTrueParticleADDPID",
         "hF0TrueParticleADDPID",
         "hF2TrueParticleADDPID" };

 double NewXbins[ nbins_pt ];
 for(int i=0;i<nbins_pt-1;i++){
	NewXbins[i] = ptmin[i];
 } NewXbins[nbins_pt-1] = ptmax[nbins_pt-1];

 double CorrectionFactor[nparticle][nbins_mult][nbins_pt];
 double CorrectionFactorErr[nparticle][nbins_mult][nbins_pt];

 for(int i=0;i<nparticle;i++){
	hMCgen = (THnSparse*)output_mc->FindObject(MCGenName[i]);
	hMCacp = (THnSparse*)output_mc->FindObject(MCAcpName[i]);
	for(int j=0;j<nbins_mult;j++){
		hMCgen->GetAxis(1)->SetRangeUser( multmin[j], multmax[j] );
		hMCacp->GetAxis(1)->SetRangeUser( multmin[j], multmax[j] );

		pTGenSpectra[i][j] = (TH1D*)hMCgen->Projection(2);
		pTAcpSpectra[i][j] = (TH1D*)hMCacp->Projection(2);

		pTGenSpectra[i][j] = (TH1D*)pTGenSpectra[i][j]->Rebin( nbins_pt-1, "", NewXbins );
		pTAcpSpectra[i][j] = (TH1D*)pTAcpSpectra[i][j]->Rebin( nbins_pt-1, "", NewXbins );

		pTCorrectionSpectra[i][j] = (TH1D*)pTAcpSpectra[i][j]->Clone(0);
		pTCorrectionSpectra[i][j]->Divide( pTGenSpectra[i][j] );

		for(int k=0;k<pTCorrectionSpectra[i][j]->GetNbinsX();k++){
			CorrectionFactor[i][j][k] = pTCorrectionSpectra[i][j]->GetBinContent(k+1);
			CorrectionFactorErr[i][j][k] = pTCorrectionSpectra[i][j]->GetBinError(k+1);
		}
		CorrectionFactor[i][j][pTCorrectionSpectra[i][j]->GetNbinsX()] =
                	( pTAcpSpectra[i][j]->GetEntries() )/( pTGenSpectra[i][j]->GetEntries() );
                CorrectionFactorErr[i][j][pTCorrectionSpectra[i][j]->GetNbinsX()] =
                	( pTAcpSpectra[i][j]->GetEntries() )/( pTGenSpectra[i][j]->GetEntries() )*
                        sqrt( 1.0/pTAcpSpectra[i][j]->GetEntries() + 1.0/pTGenSpectra[i][j]->GetEntries() );
	}
 }
//tracking efficiency and geometric acceptance

 double Integration[nbins_mult][nbins_pt];
 double IntegrationErr[nbins_mult][nbins_pt];

 double InvYieldF0[nbins_mult][nbins_pt];
 double InvYieldF0Err[nbins_mult][nbins_pt]; 

 double InvYieldF2[nbins_mult][nbins_pt];
 double InvYieldF2Err[nbins_mult][nbins_pt];


 double Br = 0.46;
 double BrErr = 0.06;

 for(int i=0;i<nbins_mult;i++){
	for(int j=0;j<nbins_pt;j++){
		Integration[i][j] = fitBinComp[i][j][1]->Integral(
			fitBinComp[i][j][1]->GetParameter(4)-
			fitBinComp[i][j][1]->GetParameter(6)*2.0,
			fitBinComp[i][j][1]->GetParameter(4)+
			fitBinComp[i][j][1]->GetParameter(6)*2.0 );
		Integration[i][j] /= hMassSigLS[i][j]->GetBinWidth(1);

		IntegrationErr[i][j] = fitBinComp[i][j][1]->IntegralError(
			fitBinComp[i][j][1]->GetParameter(4)-
			fitBinComp[i][j][1]->GetParameter(6)*2.0,
			fitBinComp[i][j][1]->GetParameter(4)+
			fitBinComp[i][j][1]->GetParameter(6)*2.0,
			fitBinComp[i][j][1]->GetParameters(),
			fitResults[i][j]->GetCovarianceMatrix().GetMatrixArray() );
		IntegrationErr[i][j] /= hMassSigLS[i][j]->GetBinWidth(1);

	
		InvYieldF0[i][j] = (1.0/MultEvt[i]) * (1.0/(ptmax[j]-ptmin[j])) * Integration[i][j] * TrigEff[i] / (CorrectionFactor[1][i][j]*Br);
		InvYieldF0Err[i][j] = sqrt( 
			pow( (1.0/MultEvt[i]) * (1.0/(ptmax[j]-ptmin[j])) * IntegrationErr[i][j] * TrigEff[i] / (CorrectionFactor[1][i][j]*Br),2 ) +
			pow( (1.0/MultEvt[i]) * (1.0/(ptmax[j]-ptmin[j])) * Integration[i][j] * TrigEff[i] * CorrectionFactorErr[1][i][j] / (pow(CorrectionFactor[1][i][j],2)*Br),2 ) +
			pow( (1.0/MultEvt[i]) * (1.0/(ptmax[j]-ptmin[j])) * Integration[i][j] * TrigEff[i] * BrErr / (CorrectionFactor[1][i][j]*Br*Br),2 ) ) ;

	}
 }
//Yield Extraction




 TH1D* hF0Yield[nbins_mult];
 for(int i=0;i<nbins_mult;i++){
	hF0Yield[i] = new TH1D("","",nbins_pt-1,0,100);
	hF0Yield[i] = (TH1D*)hF0Yield[i]->Rebin( nbins_pt-1, Form("Yield of #font[22]{f}_{0}, %.0lf-%.0lf %%",multmin[i],multmax[i]), NewXbins );
	hF0Yield[i]->SetTitle( Form("Yield of #font[22]{f}_{0}, %.0lf-%.0lf %%",multmin[i],multmax[i]) );
        hF0Yield[i]->SetName(Form("hy_%d",i));
	for(int j=0;j<nbins_pt-1;j++){
		hF0Yield[i]->SetBinContent(j+1, InvYieldF0[i][j] );
		hF0Yield[i]->SetBinError(j+1, InvYieldF0Err[i][j] );
	}
	hF0Yield[i]->SetMarkerColor(kBlack);
	hF0Yield[i]->SetLineColor(kBlack);
	hF0Yield[i]->SetMarkerStyle(20);
	hF0Yield[i]->GetXaxis()->SetTitle("#font[12]{p}_{#font[22]{T}} (GeV/#font[12]{c})");
	hF0Yield[i]->GetYaxis()->SetTitle("(#varepsilon_{trig}/N_{evt}) (d^{2}N/d#font[12]{p}_{#font[22]{T}}dy) (1/Correction) (1/BR)");

	ca->cd();
	ca->SetLogy();
	hF0Yield[i]->Draw();
	ca->SaveAs(Form("./GetSignals_13TeV/f0Yield%d.pdf",i));
 }

 ca->cd();
 gPad->SetLogy(0);
 TLegend* legf0width = new TLegend(0.7,0.7,0.9,0.9);
 legf0width->SetNColumns(2);

 for(int i=0;i<nbins_mult;i++){
	gFitPar[6][i]->SetMaximum(0.11);
	gFitPar[6][i]->SetMinimum(0.00);
	gFitPar[6][i]->SetLineColor(i+1);
	gFitPar[6][i]->SetMarkerColor(i+1);
	gFitPar[6][i]->SetMarkerStyle(20+i);
	gFitPar[6][i]->GetXaxis()->SetTitle("#font[12]{p}_{#font[22]{T}} (GeV/#font[12]{c})");
	gFitPar[6][i]->GetYaxis()->SetTitle("f0 Width (GeV/#font[12]{c}^{2})");
	if(i==0) gFitPar[6][i]->Draw("AP");
	gFitPar[6][i]->Draw("P");
	legf0width->AddEntry( gFitPar[6][i], Form("%.0lf-%.0lf %%",multmin[i],multmax[i] ),"lp");
 }
 legf0width->Draw("same");

 TLine* f0WMinLine = new TLine(gFitPar[6][0]->GetX()[0] - gFitPar[6][0]->GetEX()[0], f0WMin,
	gFitPar[6][0]->GetX()[ gFitPar[6][0]->GetN()-1 ] + gFitPar[6][0]->GetEX()[ gFitPar[6][0]->GetN()-1 ], f0WMin);
 TLine* f0WMaxLine = new TLine(gFitPar[6][0]->GetX()[0] - gFitPar[6][0]->GetEX()[0], f0WMax,
        gFitPar[6][0]->GetX()[ gFitPar[6][0]->GetN()-1 ] + gFitPar[6][0]->GetEX()[ gFitPar[6][0]->GetN()-1 ], f0WMax);

 f0WMinLine->SetLineStyle(2);
 f0WMinLine->Draw("same");
 f0WMaxLine->SetLineStyle(2);
 f0WMaxLine->Draw("same");

 ca->SaveAs(Form("./GetSignals_13TeV/f0WidthAll.pdf"));

 legf0width->Clear();
 for(int i=0;i<nbins_mult;i++){
        gFitPar[4][i]->SetMaximum(1.00);
        gFitPar[4][i]->SetMinimum(0.91);
        gFitPar[4][i]->SetLineColor(i+1);
        gFitPar[4][i]->SetMarkerColor(i+1);
        gFitPar[4][i]->SetMarkerStyle(20+i);
        gFitPar[4][i]->GetXaxis()->SetTitle("#font[12]{p}_{#font[22]{T}} (GeV/#font[12]{c})");
        gFitPar[4][i]->GetYaxis()->SetTitle("f0 Mass (GeV/#font[12]{c}^{2})");
        if(i==0) gFitPar[4][i]->Draw("AP");
        gFitPar[4][i]->Draw("P");
        legf0width->AddEntry( gFitPar[4][i], Form("%.0lf-%.0lf %%",multmin[i],multmax[i] ),"lp");
 }
 legf0width->Draw("same");

 TLine* f0MMinLine = new TLine(gFitPar[6][0]->GetX()[0] - gFitPar[6][0]->GetEX()[0], f0MassMinRange,
        gFitPar[6][0]->GetX()[ gFitPar[6][0]->GetN()-1 ] + gFitPar[6][0]->GetEX()[ gFitPar[6][0]->GetN()-1 ], f0MassMinRange);
 TLine* f0MMaxLine = new TLine(gFitPar[6][0]->GetX()[0] - gFitPar[6][0]->GetEX()[0], f0MassMaxRange,
        gFitPar[6][0]->GetX()[ gFitPar[6][0]->GetN()-1 ] + gFitPar[6][0]->GetEX()[ gFitPar[6][0]->GetN()-1 ], f0MassMaxRange);

 f0MMinLine->SetLineStyle(2);
 f0MMinLine->Draw("same");
 f0MMaxLine->SetLineStyle(2);
 f0MMaxLine->Draw("same");

 ca->SaveAs(Form("./GetSignals_13TeV/f0MassAll.pdf"));

 double XdatForComp[13] = {
        0.15384615384615374, 0.4556213017751478, 0.8106508875739644, 1.2485207100591715, 1.7573964497041423,
        2.2544378698224854, 2.751479289940828, 3.248520710059172, 3.757396449704142, 4.502958579881657,
        5.497041420118343, 6.502958579881657, 7.497041420118343 };
 double YdatForComp[13] = {
        0.013151371960322117, 0.030136976430556763, 0.02716962231219719, 0.0125799367222851, 0.006957348048080166,
        0.004022551451063787, 0.0018081511281811573, 0.00089928792213658385, 0.00039928792213658385, 0.0002598908543435751,
        0.00008067749620268826, 0.00004022551451063779, 0.00002128012999008446 };

 TGraph* gF0IncYieldAt5TeV = new TGraph( 13, XdatForComp, YdatForComp );

 ca->cd();
 gPad->SetLogy();
 hF0Yield[nbins_mult-1]->SetMaximum(0.1);
 hF0Yield[nbins_mult-1]->Draw("EP");
 gF0IncYieldAt5TeV->SetMarkerColor(kRed);
 gF0IncYieldAt5TeV->SetLineColor(kRed);
 gF0IncYieldAt5TeV->SetMarkerStyle(21);
 gF0IncYieldAt5TeV->Draw("P");

 ca->SaveAs("./GetSignals_13TeV/YieldComp.pdf");

 
 TFile* fout = new TFile("./ResultFigures/SigOut_13TeV.root","recreate");
 for(int i=0;i<nbins_mult;i++){
	hF0Yield[i]->Write();
	gFitPar[6][i]->Write();
 }
 gF0IncYieldAt5TeV->Write();
 }
