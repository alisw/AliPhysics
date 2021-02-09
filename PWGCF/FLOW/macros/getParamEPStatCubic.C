#include <TMinuit.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLegend.h>
#include <TFile.h>

#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <vector>
using namespace std;

//#include "gsl/gsl_sf_hyperg.h"

double v2i,v2ei,v4i,v4ei,v6i,v6ei,v8i,v8ei;

//=====================
Double_t eps22EP(Double_t alpha, Double_t eps0)
{

    Double_t hypergf1 = ROOT::Math::hyperg(1.5, 1., alpha+2., eps0*eps0);
    
    //Double_t hypergf1 = gsl_sf_hyperg_2F1(1.5, 1., alpha+2., eps0*eps0);
    
  Double_t qc2 = 1. - alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1;
  Double_t eps22 = TMath::Sqrt(qc2);

  return eps22;

}


Double_t eps24EP(Double_t alpha, Double_t eps0)
{

    Double_t hypergf1 = ROOT::Math::hyperg(1.5, 1., alpha+2., eps0*eps0);
    Double_t hypergf2 = ROOT::Math::hyperg(2.5, 2., alpha+3., eps0*eps0);
    
    //Double_t hypergf1 = gsl_sf_hyperg_2F1(1.5, 1., alpha+2., eps0*eps0);
    //Double_t hypergf2 = gsl_sf_hyperg_2F1(2.5, 2., alpha+3., eps0*eps0);
    
    
  Double_t qc4 = 1. - 2.*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1 + 2.*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1 - alpha/(alpha+2.)*(1.-eps0*eps0)*(1.-eps0*eps0)*hypergf2;
  Double_t eps24 = TMath::Power(qc4, 1./4.);

  return eps24;

}


Double_t eps26EP(Double_t alpha, Double_t eps0)
{

    Double_t hypergf1 = ROOT::Math::hyperg(1.5, 1., alpha+2., eps0*eps0);
    Double_t hypergf2 = ROOT::Math::hyperg(2.5, 2., alpha+3., eps0*eps0);
    Double_t hypergf3 = ROOT::Math::hyperg(3.5, 3., alpha+4., eps0*eps0);
    
    /*
    Double_t hypergf1 = gsl_sf_hyperg_2F1(1.5, 1., alpha+2., eps0*eps0);
    Double_t hypergf2 = gsl_sf_hyperg_2F1(2.5, 2., alpha+3., eps0*eps0);
    Double_t hypergf3 = gsl_sf_hyperg_2F1(3.5, 3., alpha+4., eps0*eps0);
    */
    
    
  Double_t qc6 = 1. + 9./2.*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1 - 3.*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1 + 3.*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1*(3./4.*alpha/(alpha+2.)*(1.-eps0*eps0)*(1.-eps0*eps0)*hypergf2 - 1.) - 3./2.*alpha/(alpha+2.)*(1.-eps0*eps0)*(1.-eps0*eps0)*hypergf2 - 1./4.*alpha/(alpha+3.)*(1.-eps0*eps0)*(1.-eps0*eps0)*(1.-eps0*eps0)*hypergf3;
    
  Double_t eps26 = TMath::Power(qc6, 1./6.);

  return eps26;

}


Double_t eps28EP(Double_t alpha, Double_t eps0)
{
    
    Double_t hypergf1 = ROOT::Math::hyperg(1.5, 1., alpha+2., eps0*eps0);
    Double_t hypergf2 = ROOT::Math::hyperg(2.5, 2., alpha+3., eps0*eps0);
    Double_t hypergf3 = ROOT::Math::hyperg(3.5, 3., alpha+4., eps0*eps0);
    Double_t hypergf4 = ROOT::Math::hyperg(4.5, 4., alpha+5., eps0*eps0);
    
    /*
    Double_t hypergf1 = gsl_sf_hyperg_2F1(1.5, 1., alpha+2., eps0*eps0);
    Double_t hypergf2 = gsl_sf_hyperg_2F1(2.5, 2., alpha+3., eps0*eps0);
    Double_t hypergf3 = gsl_sf_hyperg_2F1(3.5, 3., alpha+4., eps0*eps0);
    Double_t hypergf4 = gsl_sf_hyperg_2F1(4.5, 4., alpha+5., eps0*eps0);
    */
    
    Double_t term2 = -288.*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1;
    Double_t term3= 144.*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1;
    Double_t term4 = -66.*alpha/(alpha+2.)*(1.-eps0*eps0)*(1.-eps0*eps0)*hypergf2;
    Double_t term5 = 18.*alpha/(alpha+2.)*(1.-eps0*eps0)*(1.-eps0*eps0)*hypergf2*alpha/(alpha+2.)*(1.-eps0*eps0)*(1.-eps0*eps0)*hypergf2;
    Double_t term6 = -24.*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1*(-11. + 6.*alpha/(alpha+2.)*(1.-eps0*eps0)*(1.-eps0*eps0)*hypergf2);
    Double_t term7 = -12.*alpha/(alpha+3.)*(1.-eps0*eps0)*(1.-eps0*eps0)*(1.-eps0*eps0)*hypergf3;
    Double_t term8 = 4.*alpha/(alpha+1.)*(1.-eps0*eps0)*hypergf1*(-33. + 42.*alpha/(alpha+2.)*(1.-eps0*eps0)*(1.-eps0*eps0)*hypergf2 + 4.*alpha/(alpha+3.)*(1.-eps0*eps0)*(1.-eps0*eps0)*(1.-eps0*eps0)*hypergf3);
    Double_t term9 = -alpha/(alpha+4.)*(1.-eps0*eps0)*(1.-eps0*eps0)*(1.-eps0*eps0)*(1.-eps0*eps0)*hypergf4;
    
    Double_t qc8 = (33. + term2 + term3+ term4 + term5 + term6 + term7 + term8 + term9)/33.;
    
    Double_t eps28 = TMath::Power(qc8, 1./8.);
    
    return eps28;
    
    
}


//==================================



//==============================================================================
void fcn(Int_t &npar, Double_t *gin, Double_t &chisq, Double_t *par, Int_t iflag) 
{ 
  double eps0    = par[0];
  double alpha   = par[1];
  double kapa    = par[2];
  double kapapr  = par[3];

  //double eps0 = v/kapa;
  //double alpha = sigm/kapa;

  double v22ep = kapa*eps22EP(alpha, eps0) + kapapr*kapa*eps22EP(alpha, eps0)*eps22EP(alpha, eps0)*eps22EP(alpha, eps0);
  double v24ep = kapa*eps24EP(alpha, eps0) + kapapr*kapa*eps24EP(alpha, eps0)*eps24EP(alpha, eps0)*eps24EP(alpha, eps0);
  double v26ep = kapa*eps26EP(alpha, eps0) + kapapr*kapa*eps26EP(alpha, eps0)*eps26EP(alpha, eps0)*eps26EP(alpha, eps0);
  double v28ep = kapa*eps28EP(alpha, eps0) + kapapr*kapa*eps28EP(alpha, eps0)*eps28EP(alpha, eps0)*eps28EP(alpha, eps0);

  double d1 = (v22ep - v2i)/v2ei;  
  double d2 = (v24ep - v4i)/v4ei;  
  double d3 = (v26ep - v6i)/v6ei;
  double d4 = (v28ep - v8i)/v8ei;

  chisq = d1*d1+d2*d2+d3*d3+d4*d4;
   //   cout <<"chisq="<< chisq << endl;
}
//================================================================

//          fitting 

//===========================================================================

void meanvEP246(double v2, double v2e, 
		double v4, double v4e,
		double v6, double v6e,
        double v8, double v8e,
		Double_t& eps0, Double_t& eps0Err,
		Double_t& alpha, Double_t& alphaErr,
		Double_t& kapa, Double_t& kapaErr,
        Double_t& kapapr, Double_t& kapaprErr)
{

  v2i  = v2;
  v2ei = v2e;
  v4i  = v4;
  v4ei = v4e;
  v6i  = v6;
  v6ei = v6e;
  v8i  = v8;
  v8ei = v8e;

  cout << "v2="<<v2<<"  v2e=" <<v2e<< endl;
  cout << "v4="<<v4<<"  v4e=" <<v4e<< endl;
  cout << "v6="<<v6<<"  v6e=" <<v6e<< endl;
  cout << "v8="<<v8<<"  v8e=" <<v8e<< endl;
     //     minuit->mnexcm("MIGRAD", arglist ,1,ierflg);
  if(v2 > v4 ) {
    cout<<"STOP"<<endl;
  }

  //initialize TMinuit with a maximum of 3 params
  TMinuit *minuit = new TMinuit(4);
  minuit->SetFCN(fcn);

  //   minuit->mnexcm("SET ERR", arglist ,1,ierflg);
  // Set starting values and step sizes for parameters
  Double_t arglist[10];
  Int_t ierflg = 0;
  // Now ready for minimization step
  arglist[0] = 2000;
  arglist[1] = 1;

  
  Double_t step[4] = {0.01 , 0.01 , 0.01, 0.01};
  Double_t lowLim[4] = {0, 0, 0, 0};
  Double_t upLim[4] = {1, 500, 1, 0.3};//kpr can be 1

    Double_t vstart[4] = {0.3, 30, 0.5, 0.1};
    //Double_t vstart[4] = {0.3, 30, 0.4, 0.1};// for kpr range 0-1 (in lowLim/upLim above)
    
  minuit->mnparm(0, "eps0", vstart[0], step[0], lowLim[0], upLim[0],ierflg);
  minuit->mnparm(1, "alpha", vstart[1], step[1], lowLim[1], upLim[1],ierflg);
  minuit->mnparm(2, "kapa", vstart[2], step[2], lowLim[2], upLim[2],ierflg);
    minuit->mnparm(3, "kapapr", vstart[3], step[3], lowLim[3], upLim[3],ierflg);
  
  minuit->SetPrintLevel(0);
  
  minuit->mnexcm("MIGRAD", arglist , 2, ierflg);

  
  Int_t nfits = 0;
  TString status = minuit->fCstatu.Data();

  while ((!status.Contains("CONVERGED")) && (nfits < 10)){

    Double_t vstartn[4] = {0.05*(nfits+1), 10.*(nfits+1), 0.05*(nfits+1), 0.1};

    minuit->mnparm(0, "eps0", vstartn[0], step[0], lowLim[0], upLim[0],ierflg);
    minuit->mnparm(1, "alpha", vstartn[1], step[1], lowLim[1], upLim[1],ierflg);
    minuit->mnparm(2, "kapa", vstartn[2], step[2], lowLim[2], upLim[2],ierflg);
      minuit->mnparm(3, "kapapr", vstartn[3], step[3], lowLim[3], upLim[3],ierflg);
    
    minuit->SetPrintLevel(0);

    minuit->mnexcm("MIGRAD", arglist , 2, ierflg);

    status = minuit->fCstatu.Data();
    nfits++;
  }


  //  minuit->Migrad();
  minuit->GetParameter(0, eps0, eps0Err);
  minuit->GetParameter(1, alpha, alphaErr);
  minuit->GetParameter(2, kapa, kapaErr);
    minuit->GetParameter(3, kapapr, kapaprErr);

    cout << " eps0 = " <<eps0<<" +/- "<<eps0Err<<endl;
    cout << " alpha = "<<alpha<<" +/- "<<alphaErr<<endl;
    cout << " kapa = " <<kapa<<" +/- "<<kapaErr<<endl;
    cout << " kapapr = " <<kapapr<<" +/- "<<kapaprErr<<endl;

  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  //  minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //  minuit->mnprin(3,amin);
}


void getParamEP(Int_t optDrawEP = 1)
{

    /*
    gSystem->Load("/usr/local/lib/libgsl.dylib");
    gSystem->Load("/usr/local/lib/libgslcblas.dylib");
    gSystem->SetIncludePath("-I/usr/local/include");
    */
    
    
  //for root
  gSystem->Load("libMathCore.so");
  gSystem->Load("libMathMore.so"); 
  gSystem->AddIncludePath("-I$ROOTSYS/include"); 


    ifstream inFile;
    inFile.open("v22468Run2JacopoStat.txt");


  double cent, v2sp, v2sper, v24, v24er, v26, v26er, v28, v28er;
  vector<double> vcent, vv22, vv22er, vv24, vv24er, vv26, vv26er, vv28, vv28er;
  
  while (!inFile.eof()) {
    inFile>>cent>>v2sp>>v2sper>>v24>>v24er>>v26>>v26er>>v28>>v28er;
    vcent.push_back(cent);
    vv22.push_back(v2sp);
    vv22er.push_back(v2sper);
    vv24.push_back(v24);
    vv24er.push_back(v24er);
    vv26.push_back(v26);
    vv26er.push_back(v26er);
    vv28.push_back(v28);
    vv28er.push_back(v28er);

  }

  inFile.close();

    
    Int_t nPoints = vcent.size();
    
    cout<<"Number of points = "<<nPoints<<endl;
  
    TGraphErrors* grEps0 = new TGraphErrors(nPoints-1);
    grEps0->SetName("grEps0");
    TGraphErrors* grAlpha = new TGraphErrors(nPoints-1);
    grAlpha->SetName("grAlpha");
    TGraphErrors* grKappa = new TGraphErrors(nPoints-1);
    grKappa->SetName("grKappa");
    TGraphErrors* grKappaPr = new TGraphErrors(nPoints-1);
    grKappaPr->SetName("grKappaPr");
   
    
    
    TGraphErrors* grv22 = new TGraphErrors(nPoints-1);
    TGraphErrors* grv24 = new TGraphErrors(nPoints-1);
    TGraphErrors* grv26 = new TGraphErrors(nPoints-1);
    TGraphErrors* grv28 = new TGraphErrors(nPoints-1);
        
    TGraphErrors* grv22EP = new TGraphErrors(nPoints-1);
    TGraphErrors* grv24EP = new TGraphErrors(nPoints-1);
    TGraphErrors* grv26EP = new TGraphErrors(nPoints-1);
    TGraphErrors* grv28EP = new TGraphErrors(nPoints-1);
    
    
    TGraphErrors* grv22Rat = new TGraphErrors(nPoints-1);
    grv22Rat->SetName("grv22Rat");
    TGraphErrors* grv24Rat = new TGraphErrors(nPoints-1);
    grv24Rat->SetName("grv24Rat");
    TGraphErrors* grv26Rat = new TGraphErrors(nPoints-1);
    grv26Rat->SetName("grv26Rat");
    TGraphErrors* grv28Rat = new TGraphErrors(nPoints-1);
    grv28Rat->SetName("grv28Rat");

    
    for (Int_t i = 0; i < nPoints-1; i++){
    
        cout<<"     "<<endl;
        cout<<"******"<<endl;
        cout<<"!!!!!!!!! Fitting centrality "<<vcent[i]<<endl;
        cout<<"******"<<endl;

        Double_t Eps0, Eps0Err, Alpha, AlphaErr, Kappa, KappaErr, Kappapr, KappaprErr;
        meanvEP246(vv22[i], vv22er[i], vv24[i], vv24er[i], vv26[i], vv26er[i], vv28[i], vv28er[i], Eps0, Eps0Err, Alpha, AlphaErr, Kappa, KappaErr, Kappapr, KappaprErr);

        grEps0->SetPoint(i, vcent[i], Eps0);
        grEps0->SetPointError(i, 0, Eps0Err);

        grAlpha->SetPoint(i, vcent[i], Alpha);
        grAlpha->SetPointError(i, 0, AlphaErr);

        grKappa->SetPoint(i, vcent[i], Kappa);
        grKappa->SetPointError(i, 0, KappaErr);
        
        grKappaPr->SetPoint(i, vcent[i], Kappapr);
        grKappaPr->SetPointError(i, 0, KappaprErr);

        
        
        if (optDrawEP){
        
            Double_t V22ep = Kappa*eps22EP(Alpha, Eps0) + Kappa*Kappapr*eps22EP(Alpha, Eps0)*eps22EP(Alpha, Eps0)*eps22EP(Alpha, Eps0);
            Double_t V24ep = Kappa*eps24EP(Alpha, Eps0) + Kappa*Kappapr*eps24EP(Alpha, Eps0)*eps24EP(Alpha, Eps0)*eps24EP(Alpha, Eps0);
            Double_t V26ep = Kappa*eps26EP(Alpha, Eps0) + Kappa*Kappapr*eps26EP(Alpha, Eps0)*eps26EP(Alpha, Eps0)*eps26EP(Alpha, Eps0);
            Double_t V28ep = Kappa*eps28EP(Alpha, Eps0) + Kappa*Kappapr*eps28EP(Alpha, Eps0)*eps28EP(Alpha, Eps0)*eps28EP(Alpha, Eps0);
        
        
            Double_t V22epErr = KappaErr*eps22EP(AlphaErr, Eps0Err) + KappaErr*KappaprErr*eps22EP(AlphaErr, Eps0Err)*eps22EP(AlphaErr, Eps0Err)*eps22EP(AlphaErr, Eps0Err);
            Double_t V24epErr = KappaErr*eps24EP(AlphaErr, Eps0Err) + KappaErr*KappaprErr*eps24EP(AlphaErr, Eps0Err)*eps24EP(AlphaErr, Eps0Err)*eps24EP(AlphaErr, Eps0Err);
            Double_t V26epErr = KappaErr*eps26EP(AlphaErr, Eps0Err) + KappaErr*KappaprErr*eps26EP(AlphaErr, Eps0Err)*eps26EP(AlphaErr, Eps0Err)*eps26EP(AlphaErr, Eps0Err);
            Double_t V28epErr = KappaErr*eps28EP(AlphaErr, Eps0Err) + KappaErr*KappaprErr*eps28EP(AlphaErr, Eps0Err)*eps28EP(AlphaErr, Eps0Err)*eps28EP(AlphaErr, Eps0Err);
        
      
            grv22->SetPoint(i, vcent[i], vv22[i]);
            grv22->SetPointError(i, 0, vv22er[i]);
      
            grv24->SetPoint(i, vcent[i], vv24[i]);
            grv24->SetPointError(i, 0, vv24er[i]);
      
            grv26->SetPoint(i, vcent[i], vv26[i]);
            grv26->SetPointError(i, 0, vv26er[i]);
        
            grv28->SetPoint(i, vcent[i], vv28[i]);
            grv28->SetPointError(i, 0, vv28er[i]);
        
      
            grv22EP->SetPoint(i, vcent[i], V22ep);
            grv22EP->SetPointError(i, 0, V22epErr);
      
            grv24EP->SetPoint(i, vcent[i], V24ep);
            grv24EP->SetPointError(i, 0, V24epErr);
      
            grv26EP->SetPoint(i, vcent[i], V26ep);
            grv26EP->SetPointError(i, 0, V26epErr);
        
            grv28EP->SetPoint(i, vcent[i], V28ep);
            grv28EP->SetPointError(i, 0, V28epErr);
        
        
        
            Double_t ratV22 = V22ep/vv22[i];
            //Double_t errRatV22 = V22ep/vv22[i]*TMath::Sqrt(TMath::Abs(V22epErr*V22epErr/V22ep/V22ep + vv22er[i]*vv22er[i]/vv22[i]/vv22[i]));
            Double_t errRatV22 = TMath::Sqrt(TMath::Abs(vv22er[i]*vv22er[i] - V22epErr*V22epErr));
        
            Double_t ratV24 = V24ep/vv24[i];
            //Double_t errRatV24 = V24ep/vv24[i]*TMath::Sqrt(TMath::Abs(V24epErr*V24epErr/V24ep/V24ep + vv24er[i]*vv24er[i]/vv24[i]/vv24[i]));
            Double_t errRatV24 = TMath::Sqrt(TMath::Abs(vv24er[i]*vv24er[i] - V24epErr*V24epErr));
        
            Double_t ratV26 = V26ep/vv26[i];
            //Double_t errRatV26 = V26ep/vv26[i]*TMath::Sqrt(TMath::Abs(V26epErr*V26epErr/V26ep/V26ep + vv26er[i]*vv26er[i]/vv26[i]/vv26[i]));
            Double_t errRatV26 = TMath::Sqrt(TMath::Abs(vv26er[i]*vv26er[i] - V26epErr*V26epErr));
        
            Double_t ratV28 = V28ep/vv28[i];
            //Double_t errRatV28 = V28ep/vv28[i]*TMath::Sqrt(TMath::Abs(V28epErr*V28epErr/V28ep/V28ep + vv28er[i]*vv28er[i]/vv28[i]/vv28[i]));
            Double_t errRatV28 = TMath::Sqrt(TMath::Abs(vv28er[i]*vv28er[i] - V28epErr*V28epErr));
        
        
            grv22Rat->SetPoint(i, vcent[i], ratV22);
            grv22Rat->SetPointError(i, 0, errRatV22);
        
            grv24Rat->SetPoint(i, vcent[i], ratV24);
            grv24Rat->SetPointError(i, 0, errRatV24);
        
            grv26Rat->SetPoint(i, vcent[i], ratV26);
            grv26Rat->SetPointError(i, 0, errRatV26);
        
            grv28Rat->SetPoint(i, vcent[i], ratV28);
            grv28Rat->SetPointError(i, 0, errRatV28);
        
        
        }

  }


  TCanvas* cEps0 = new TCanvas("cEps0", "cEps0");
  cEps0->cd();

  TH1D* hdumeps = new TH1D("hdumeps", "; centrality percentile; #epsilon_{0}", 100, 0, 100);
  hdumeps->SetMaximum(1);
  hdumeps->SetMinimum(0);
  hdumeps->Draw();

  grEps0->SetMarkerStyle(20);
  grEps0->SetMarkerColor(1);
  grEps0->SetLineColor(1);
  grEps0->Draw("sameP");


  TCanvas* cAlph = new TCanvas("cAlph", "cAlph");
  cAlph->cd();

  TH1D* hdumalp = new TH1D("hdumalp", "; centrality percentile; #alpha", 100, 0, 100);
  hdumalp->SetMaximum(180);
  hdumalp->SetMinimum(0);
  hdumalp->Draw();

  grAlpha->SetMarkerStyle(20);
  grAlpha->SetMarkerColor(1);
  grAlpha->SetLineColor(1);
  grAlpha->Draw("sameP");
  


  TCanvas* cKap = new TCanvas("cKap", "cKap");
  cKap->cd();

  TH1D* hdumkap = new TH1D("hdumkap", "; centrality percentile; #Kappa", 100, 0, 100);
  hdumkap->SetMaximum(1);
  hdumkap->SetMinimum(0);
  hdumkap->Draw();

  grKappa->SetMarkerStyle(20);
  grKappa->SetMarkerColor(1);
  grKappa->SetLineColor(1);
  grKappa->Draw("sameP");
    
    
    
    TCanvas* cKappr = new TCanvas("cKappr", "cKappr");
    cKappr->cd();
    
    TH1D* hdumkappr = new TH1D("hdumkappr", "; centrality percentile; #Kappa^{'}", 100, 0, 100);
    hdumkappr->SetMaximum(1);
    hdumkappr->SetMinimum(0);
    hdumkappr->Draw();
    
    grKappaPr->SetMarkerStyle(20);
    grKappaPr->SetMarkerColor(1);
    grKappaPr->SetLineColor(1);
    grKappaPr->Draw("sameP");
    
    
    
   if (optDrawEP){
       
       TCanvas* cv2 = new TCanvas("cv2", "cv2");
       cv2->cd();
       
       TH1D* hdumv2 = new TH1D("hdumv2", "; centrality percentile; v_{2}", 100, 0, 100);
       hdumv2->SetMaximum(0.15);
       hdumv2->SetMinimum(0);
       hdumv2->Draw();
       
       grv22->SetLineColor(1);
       grv22->SetMarkerColor(1);
       grv22->SetMarkerStyle(20);
       grv22->Draw("Psame");
       
       grv22EP->SetLineColor(2);
       grv22EP->SetMarkerColor(2);
       grv22EP->SetMarkerStyle(24);
       grv22EP->Draw("Psame");
       
       
       grv24->SetLineColor(1);
       grv24->SetMarkerColor(1);
       grv24->SetMarkerStyle(21);
       grv24->Draw("Psame");
       
       grv24EP->SetLineColor(2);
       grv24EP->SetMarkerColor(2);
       grv24EP->SetMarkerStyle(25);
       grv24EP->Draw("Psame");
       
       
       grv26->SetLineColor(kGreen+2);
       grv26->SetMarkerColor(kGreen+2);
       grv26->SetMarkerStyle(22);
       grv26->Draw("Psame");
       
       grv26EP->SetLineColor(4);
       grv26EP->SetMarkerColor(4);
       grv26EP->SetMarkerStyle(26);
       grv26EP->Draw("Psame");
       
       
       grv28->SetLineColor(kMagenta+2);
       grv28->SetMarkerColor(kMagenta+2);
       grv28->SetMarkerStyle(23);
       grv28->Draw("Psame");
       
       grv28EP->SetLineColor(kCyan+2);
       grv28EP->SetMarkerColor(kCyan+2);
       grv28EP->SetMarkerStyle(32);
       grv28EP->Draw("Psame");
       
       
       TLegend* lg = new TLegend(0.12, 0.53, 0.28, 0.88);
       lg->SetFillColor(0);
       lg->SetBorderSize(0);
       lg->AddEntry(grv22, "v_{2}{2, |#Delta#eta|>2} (Data)", "LP");
       lg->AddEntry(grv24, "v_{2}{4} (Data)", "LP");
       lg->AddEntry(grv26, "v_{2}{6} (Data)", "LP");
       lg->AddEntry(grv28, "v_{2}{8} (Data)", "LP");
       lg->AddEntry(grv22EP, "v_{2}{2, |#Delta#eta|>2} (Fit)", "LP");
       lg->AddEntry(grv24EP, "v_{2}{4} (Fit)", "LP");
       lg->AddEntry(grv26EP, "v_{2}{6} (Fit)", "LP");
       lg->AddEntry(grv28EP, "v_{2}{8} (Fit)", "LP");
       lg->Draw();
       
       //cv2->SaveAs("v2_comp_data_fit.png");
       
       
       
       TCanvas* cv2r = new TCanvas("cv2r", "cv2r");
       cv2r->SetGridy();
       cv2r->cd();
       
       TH1D* hdumv2r = new TH1D("hdumv2r", "; centrality percentile; v_{2} (fit/data)", 100, 0, 100);
       hdumv2r->SetMaximum(1.03);
       hdumv2r->SetMinimum(0.95);
       //hdumv2r->SetMaximum(1.03);
       //hdumv2r->SetMinimum(0.96);
       hdumv2r->Draw();
       
       grv22Rat->SetLineColor(1);
       grv22Rat->SetMarkerColor(1);
       grv22Rat->SetMarkerStyle(20);
       grv22Rat->Draw("Psame");
       
       grv24Rat->SetLineColor(2);
       grv24Rat->SetMarkerColor(2);
       grv24Rat->SetMarkerStyle(25);
       grv24Rat->Draw("Psame");
       
       grv26Rat->SetLineColor(4);
       grv26Rat->SetMarkerColor(4);
       grv26Rat->SetMarkerStyle(26);
       grv26Rat->Draw("Psame");
       
       grv28Rat->SetLineColor(kMagenta+2);
       grv28Rat->SetMarkerColor(kMagenta);
       grv28Rat->SetMarkerStyle(32);
       grv28Rat->Draw("Psame");
       
       
       TLegend* lgr = new TLegend(0.3, 0.68, 0.5, 0.88);
       lgr->SetFillColor(0);
       lgr->SetBorderSize(0);
       lgr->AddEntry(grv22Rat, "v_{2}{2, |#Delta#eta|>2}", "LP");
       lgr->AddEntry(grv24Rat, "v_{2}{4}", "LP");
       lgr->AddEntry(grv26Rat, "v_{2}{6}", "LP");
       lgr->AddEntry(grv28Rat, "v_{2}{8}", "LP");
       lgr->Draw();
       
       //cv2r->SaveAs("v2_comp_data_fit_rat.png");
       
    }
    

    //TFile* out = new TFile("param_fitEP_run2_JacopoStatCubicKpr01.root", "RECREATE");
    TFile* out = new TFile("param_fitEP_run2_JacopoStatCubic.root", "RECREATE");
    grEps0->Write();
    grAlpha->Write();
    grKappa->Write();
    grKappaPr->Write();
    if (optDrawEP){
        grv22Rat->Write();
        grv24Rat->Write();
        grv26Rat->Write();
        grv28Rat->Write();
    }
    out->Close();
    
    

}
