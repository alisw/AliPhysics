/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

//  -------------------------------------------------------
//  pt spectra extrapolation in the 0. - 0.2 region using
//  Boltzmann-Gibbs Blast Wave model or Tsallis Blast
//  Wave model for azimuthal isotropic  expansion in
//  highly central collisions analysis
//  author: Cristian Andrei
//          acristian@niham.nipne.ro
//  ----------------------------------------------------------


#include "TH1D.h"
#include "TFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
// #include "TString.h"
// #include "TPaveStats.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TF1.h"

#include "TVirtualFitter.h"
#include "Math/WrappedFunction.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"

#include "AliCFContainer.h"
#include "AliCFDataGrid.h"
#include "AliCFEffGrid.h"

#include "AliAnalysisCentralExtrapolate.h"


ClassImp(AliAnalysisCentralExtrapolate)

//________________________________________________________________________
AliAnalysisCentralExtrapolate::AliAnalysisCentralExtrapolate(const char *name) 
  : TObject()
  ,fPartType()
  ,fInputList(0)
  ,fResultsList(0x0)

{
// Constructor

//  printf("AliAnalysisCentralExtrapolate::AliAnalysisCentralExtrapolate(const char *name)\n");
	fInputList = new TList();
	fResultsList = new TList();
	
}
  
AliAnalysisCentralExtrapolate::~AliAnalysisCentralExtrapolate() {
	
	if(fInputList) delete fInputList;
	if(fResultsList) delete fResultsList;

}

//************************************************************************************
//____________________________________________________________________________________
// using global variables to avoid the limitations of ROOT::Math::WrappedMultiFunction<>
static Double_t mass = 0.;
static Double_t pt = 0.;
static Double_t T = 0.;
static Double_t betaS = 0.;
static Double_t q = 0.;
static Double_t mt = 0.;

void AliAnalysisCentralExtrapolate::ApplyEff(){
// applies the efficiency map to the pt spectra

    TH1::SetDefaultSumw2();
//     Int_t stepGen = 0;
    Int_t stepRec = 1;

	AliCFContainer *data = 0;
	
	
	TFile *file1 = new TFile("$ALICE_ROOT/PWG2/data/AliAnalysisCentralEfficiency.root", "read");
    
    AliCFEffGrid *eff = 0;
    
    if(fPartType.Contains("kPi")){
		data = dynamic_cast<AliCFContainer*>(fInputList->FindObject("TaskCentral_CFCont_Pi"));
		eff = dynamic_cast<AliCFEffGrid*>(file1->Get("eff_Pi"));
		mass= 0.13957;//(GeV - pions)
// 		ccorrdata =new TCanvas("Pions ccorrdata","Pions - corrected data",0,0,600,800);
		printf("\n\n**************************************\n");
		printf("\tRunning for pions!\n");
	}
	else if(fPartType.Contains("kK")){
		data = dynamic_cast<AliCFContainer*>(fInputList->FindObject("TaskCentral_CFCont_K"));
		eff = dynamic_cast<AliCFEffGrid*>(file1->Get("eff_K"));
		mass= 0.49368;//(GeV - kaons)
// 		ccorrdata =new TCanvas("Kaons ccorrdata","Kaons - corrected data",0,0,600,800);
		printf("\n\n**************************************\n");
 		printf("\tRunning for kaons!\n");
	}
	else if(fPartType.Contains("kProton")){
		data = dynamic_cast<AliCFContainer*>(fInputList->FindObject("TaskCentral_CFCont_P"));
		eff = dynamic_cast<AliCFEffGrid*>(file1->Get("eff_P"));
		mass = 0.93827;//(GeV - protons)
// 		ccorrdata =new TCanvas("Protons ccorrdata","Protons - corrected data",0,0,600,800);
		printf("\n\n**************************************\n");
 		printf("\tRunning for protons!\n");
	}
	else printf("Unsupported particle type!\n");

    if(!data){
		printf("Unable to get CFContainer! \n");
		return;
	}
	
	if(!eff){
		printf("No Eff Grid found! \n");
		return;
    }
	
// 	TCanvas *ccorrdata = new TCanvas();
// 	ccorrdata->Divide(1,2);
// 	ccorrdata->cd(1);
// 	ccorrdata->cd(1)->SetLogy();

 
    AliCFDataGrid *corrdata = new AliCFDataGrid("corrdata","corrected data",*data);

//correct selection step "reconstructed"
    corrdata->SetMeasured(stepRec); //set data to be corrected
    corrdata->ApplyEffCorrection(*eff);//apply the correction for efficiency

//     TH1D *hPtMC = data->ShowProjection(0,0); //MC distribution ShowProjection(ivar, istep)
//     hPtMC->SetMarkerStyle(20);
//     hPtMC->SetMarkerColor(kGreen-3);
//     hPtMC->GetXaxis()->SetTitle("p_{T}(GeV/c)");
//     hPtMC->GetYaxis()->SetTitle("#frac{dN}{p_{T}dp_{T}}");
//     hPtMC->Draw("p e1");
// 
//     TH1D *hPtESDI = corrdata->GetData()->Project(0); //uncorrected ESD  Project(ivar)
//     hPtESDI->SetMarkerStyle(25);
//     hPtESDI->SetMarkerColor(kBlue);
//     hPtESDI->GetXaxis()->SetTitle("Pt");
//     hPtESDI->Draw("p e1 same");

    TH1D *hPtESDCorr = corrdata->Project(0); //corrected data
    hPtESDCorr->SetMarkerStyle(26);
    hPtESDCorr->SetMarkerColor(kRed);        //ESD corrected 
	hPtESDCorr->SetName("hPtESDCorr");
// 	hPtESDCorr->Draw("p e1 same");

// 	TLegend *leg2 = new TLegend(0.40,0.65,0.87,0.8);
//     leg2->AddEntry(hPtMC," generated","p");
//     leg2->AddEntry(hPtESDI," reconstructed (not corrected)","p");    
//     leg2->AddEntry(hPtESDCorr," reconstructed & corrected","p");    
// 	leg2->Draw();


// 	ccorrdata->cd(2);
// 	ccorrdata->cd(2)->SetLogy();
// 	
// 	TH1D *efficiency = eff->Project(0); //the efficiency vs pt (ipt = 0)
// 	efficiency->GetXaxis()->SetTitle("p_{T}(GeV/c)");
// 	efficiency->Draw();


	TH1F *hNoEvt = dynamic_cast<TH1F*>(fInputList->FindObject("TaskCentral_NoEvt"));
	if(!hNoEvt){
		printf("Unable to get the number of events! \n");
		return;
	}

	
    Int_t noEvt = (Int_t)(hNoEvt->GetEntries());

    printf("\n** No of processed events = %i **\n",noEvt);

    Double_t scale = 1.0/noEvt;

    TH1D *hPtESDCorrNorm = (TH1D*)hPtESDCorr->Clone("hPtESDCorrNorm");
    hPtESDCorrNorm->Scale(scale);

	fResultsList->SetName(fPartType);
	fResultsList->Add(hPtESDCorr);
	fResultsList->Add(hPtESDCorrNorm);
	
}


//*********************************************************************************
//_________________________________________________________________________________
// fit Tsallis


Double_t Integ1(const Double_t *x){
//the function under the integral - TBW

    const Double_t rMax = 11.0;//(fm)

    Double_t betaR = betaS*(x[0]/rMax);


    Double_t rho = TMath::ATanH(betaR);

    
    Double_t temp = 1.0+((q-1.0)*(mt*TMath::CosH(x[2])*TMath::CosH(rho) - pt*TMath::SinH(rho)*TMath::Cos(x[1]))/T);

    Double_t power = -1.0/(q-1.0);

    Double_t temp1 = pow(temp, power);

    Double_t f = TMath::CosH(x[1])*x[0]*temp1;

    return f;
}


// TF1 requires the function to have the ( )( Double_t *, Double_t *) signature 
Double_t Tsallis(Double_t* const x, Double_t* const par){ 
//computes the triple integral in the TBW formula

    pt = x[0]; 

    T = par[0];
    betaS = par[1];
    q = par[2];

    mt = sqrt(mass*mass+ pt*pt);
	

   ROOT::Math::WrappedMultiFunction<> f1(Integ1,3);

    ROOT::Math::IntegratorMultiDim ig(f1, ROOT::Math::IntegrationMultiDim::kPLAIN);


    Double_t xmin[3] = {0.0, -TMath::Pi(), -0.5}; //rmin, phi_min, ymin
    Double_t xmax[3] = {11.0, TMath::Pi(), 0.5}; //rmax, phi_max, ymax

    Double_t rez = mt*ig.Integral(xmin,xmax);

    return rez;
}


void AliAnalysisCentralExtrapolate::TsallisFit(){
// fits and extrpolates the pt spectrum using TBW model

	printf("---------------------------------------------\n");
	printf("*****   Tsallis Blast Wave Fit    *****\n");

	printf("AliAnalysisCentralExtrapolate::Tsallis mass = %g\n", mass);


    TStopwatch timer; 

    TH1::SetDefaultSumw2();
    TH1::AddDirectory(0);

    timer.Start(); 

	TH1D *ptSpectra  = dynamic_cast<TH1D*>(fResultsList->FindObject("hPtESDCorrNorm"));
	if(!ptSpectra){
		printf("TsallisFit: Can't get the normalized spectrum\n");
		return;	
	}

	if(!ptSpectra->GetEntries()){
		printf("TsallisFit: The fit data is empty!\n");
		return;
	}

    TVirtualFitter::SetDefaultFitter("Minuit2");//options Minuit, Minuit2, Fumili, Fumili2

    TF1 *ptFit = new TF1("ptFit",Tsallis,0.2,4.0,3);

    gStyle->SetOptFit(1112);

//     ptFit->SetParName(0,"T");
//     ptFit->SetParName(1,"betaS");
//     ptFit->SetParName(2,"q");

    ptFit->SetParameters(0.1,0.5,1.05); //GeV
    ptFit->SetParLimits(0,0.0,0.5);//GeV
    ptFit->SetParLimits(1,0.0,1.0);
    ptFit->SetParLimits(2,1.00001,1.5);

	ptSpectra->Fit("ptFit","Q R B ME N");

    timer.Stop(); 
    timer.Print();

//  ----------- print the fit result ----------
    Double_t chi2 = ptFit->GetChisquare();
    Double_t ndf = ptFit->GetNDF();
    Double_t p1 = ptFit->GetParameter(0);
    Double_t param1Err = ptFit->GetParError(0);
    Double_t p2 = ptFit->GetParameter(1);
    Double_t param2Err = ptFit->GetParError(1);
    Double_t p3 = ptFit->GetParameter(2);
    Double_t param3Err = ptFit->GetParError(2);

    printf("chi2/ndf = %f/%f\n", chi2,ndf);
    printf("p1 = %f\tparam1Err = %f\n", p1, param1Err);
    printf("p2 = %f\tparam2Err = %f\n", p2, param2Err);
    printf("p3 = %f\tparam3Err = %f\n", p3, param3Err);


// create a new, "extended" histogram
    TH1F *hPtExtTsallis = new TH1F("PtExtTsallis","Pt Corr Norm Ext",25,0.0,5.0);
    
    Double_t bin, binerr, test;
    
    for(Int_t i=0; i<ptSpectra->GetNbinsX()+1;i++){
	
		bin = ptSpectra->GetBinContent(i);
		binerr = ptSpectra->GetBinError(i);
		test = ptSpectra->GetBinLowEdge(i);

		hPtExtTsallis->SetBinContent(i, bin);
		hPtExtTsallis->SetBinError(i, binerr);
    }
    
    Double_t eval;
    eval = ptFit->Eval(0.1);
    printf("extrapolated pt value = %g \n", eval);
    printf("****************************************************\n\n");
	
    hPtExtTsallis->SetBinContent(1, eval);

    hPtExtTsallis->SetMarkerStyle(24);
    hPtExtTsallis->SetMarkerColor(kRed);        //ESD corrected 
    hPtExtTsallis->GetXaxis()->SetTitle("Pt");


	fResultsList->Add(hPtExtTsallis);

}

//*********************************************************************************
//_________________________________________________________________________________
// fit Boltzmann

Double_t func(Double_t r){ 
//the function under the integral - BGBW

    const Double_t rMax = 11.0;//(fm)

    Double_t betaR = betaS*(r/rMax);


    Double_t rho = TMath::ATanH(betaR);


	Double_t mt = sqrt(mass*mass + pt*pt);// !!! par[0]= pt !!!!!

    Double_t argI0 = (pt*TMath::SinH(rho))/T; //T = par[1]
	if(argI0>700.0) return 0.0; // !! floating point exception protection

    Double_t argK1 = (mt*TMath::CosH(rho))/T;

	
    Double_t i0 = TMath::BesselI0(argI0);


    Double_t k1 = TMath::BesselK1(argK1);


    Double_t f = r*mt*i0*k1;
	
    return f;
}


Double_t Boltzmann(Double_t* const x, Double_t* const par){
//computes the integral in the BGBW formula

	pt = x[0];
    T = par[0];
	
    betaS = par[1];


    ROOT::Math::WrappedFunction<> f1(func);

    ROOT::Math::Integrator ig(f1, ROOT::Math::IntegrationOneDim::kGAUSS,1.E-12,1.E-12,10000);

    Double_t rez = ig.Integral(0.0,11.0);

    return rez;
}


void AliAnalysisCentralExtrapolate::BoltzmannFit(){
//fits and extrapoates the pt spectrum using the BGBW model

	printf("---------------------------------------------------\n");
	printf("*** Boltzmann-Gibbs Blast Wave Fit  ***\n");

	printf("AliAnalysisCentralExtrapolate::Boltzmann mass = %g\n", mass);


    TStopwatch timer; 

    TH1::SetDefaultSumw2();
    TH1::AddDirectory(0);

    timer.Start(); 


	TH1D *ptSpectra  = dynamic_cast<TH1D*>(fResultsList->FindObject("hPtESDCorrNorm"));
	if(!ptSpectra){
		printf("BoltzmannFit: Can't get the normalized spectrum\n");
		return;	
	}

	printf("pt spectra get entries: %f\n",ptSpectra->GetEntries());

	if(!ptSpectra->GetEntries()){
		printf("BoltzmannFit: The fit data is empty!\n");
		return;
	}

    gStyle->SetOptStat("neRM");

    TVirtualFitter::SetDefaultFitter("Minuit2");

    TF1 *ptFit = new TF1("ptFit",Boltzmann,0.2,4.0,2); 

    gStyle->SetOptFit(1112);

//     ptFit->SetParName(0,"T");
//     ptFit->SetParName(1,"betaS");

    ptFit->SetParameters(0.1,0.5); //GeV
    ptFit->SetParLimits(0,0.001,0.5);//GeV
    ptFit->SetParLimits(1,0.0,1.0);

	ptSpectra->Fit("ptFit","Q R B ME I N");

    timer.Stop();
    timer.Print();

//  ----------- print the fit results ----------
    Double_t chi2 = ptFit->GetChisquare();
    Double_t ndf = ptFit->GetNDF();
    Double_t p1 = ptFit->GetParameter(0);
    Double_t param1Err = ptFit->GetParError(0);
    Double_t p2 = ptFit->GetParameter(1);
    Double_t param2Err = ptFit->GetParError(1);

    printf("chi2/ndf = %f/%f = %f\n", chi2,ndf,chi2/ndf);
    printf("p1 = %f (GeV)\tparam1Err = %f\n", p1, param1Err);
    printf("p2 = %f\tparam2Err = %f\n", p2, param2Err);


    TH1F *hPtExtBoltzmann = new TH1F("PtExtBoltzmann","Pt Corr Norm Ext",25,0.0,5.0);

    Double_t bin, binerr, test;

    for(Int_t i=0; i<ptSpectra->GetNbinsX()+1;i++){
	
		bin = ptSpectra->GetBinContent(i);
		binerr = ptSpectra->GetBinError(i);
		test = ptSpectra->GetBinLowEdge(i);
	
		hPtExtBoltzmann->SetBinContent(i, bin);
		hPtExtBoltzmann->SetBinError(i, binerr);
    }

    Double_t eval;
    eval = ptFit->Eval(0.1);
    printf("extrapolated pt value = %g \n", eval);
    printf("********************************************\n\n");

    hPtExtBoltzmann->SetBinContent(1, eval);

    hPtExtBoltzmann->SetMarkerStyle(24);
    hPtExtBoltzmann->SetMarkerColor(kRed);        //ESD corrected 
    hPtExtBoltzmann->GetXaxis()->SetTitle("Pt");


	fResultsList->Add(hPtExtBoltzmann);

}


