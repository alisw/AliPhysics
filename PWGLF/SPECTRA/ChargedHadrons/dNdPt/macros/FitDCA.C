#include "TROOT.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TFractionFitter.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "AlidNdPtCutAnalysis.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooGlobalFunc.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooHistFunc.h"
#include "RooRealSumPdf.h"
using namespace RooFit ;

/*
  Macro to perform templates fit of DCA distributions.
  It is supposed to analyse output root files from task PWGLF/SPECTRA/ChargedHadrons/dNdPt/AlidNdPtCutAnalysis.cxx
  For use with other tasks/other analysis: there are two alternative functions for the fits (performRooFit, performFit) that can be copied and used standalone, 
  see explanation inside the functions to know the input needed.
  
  
   Usage of this macro:
   1) compilation 
   gSystem->SetIncludePath("-I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/../src/PWGLF/SPECTRA/ChargedHadrons/dNdPt/ -I$ALICE_PHYSICS/../src/OADB -I$ALICE_PHYSICS/../src/PWGUD/base/")
   .L Fit.C+

   2) Example of usage - run with the default parameters: 
   Fit("test",2,80,100)
   it will run the two-template fit selecting the centrality range 80-100
 
   See function Fit for the explanation of arguments
   See function Examples for examples of usage 
   Set the global variables to change the fit behavior

   Origin: federica.sozzi@cern.ch
*/

// --------------------------------
// switches to change fit behavior
// ---------------------------------
Bool_t rebin = false; // Set to true to rebin the histograms
Bool_t useInFrac = true;// Set to true to use integrals of MC histograms as initial conditions for the fit parameters

Bool_t rooFit = true;//choose rooFit minimizer (default); if set to false, use TFractionFitter

// Bins and histograms for the pt projections:
 const Int_t Nbins = 8;
 Double_t binsPt[Nbins+1] = {0.1,0.5,1,1.5,2,3,4,5,10};
//const Int_t Nbins = 3;
//Double_t binsPt[Nbins+1] = {0.1,0.5,1,1.5};


//Max. number of templates foreseen :
const Int_t maxTemp = 5 ; 
//number of templates in the fit - it is changed in the Fit function, depending on the type of fit chosen via input parameter of Fit
Int_t nTemp = 3;
// Constrain on the parameters:  Used only in case of TFractionFitter (RooFit uses recursive parameters)
Double_t minV[maxTemp] = {0.7, 0.0 , 0, 0, 0};
Double_t maxV[maxTemp] = {1, 0.3, 0.2,0.2,0.05};
//..and their initial values (these are overwritten  if useInFrac is set to true)
Double_t val[maxTemp] = {0.75, 0.15 , 0.05,0.04,0.01};
	
//vector to store the fraction of the MC true
Double_t frac[maxTemp] = {0,0,0,0,0}; 
Double_t fracCorr[maxTemp] = {0,0,0,0,0}; 
//vector to store the fraction of MC templates really used
Double_t fracTemp[maxTemp] = {0,0,0,0,0}; 
Double_t fracTempCorr[maxTemp] = {0,0,0,0,0}; 


//these are relative weights for Lambda and K0, for the templates type "4".
//It seems that they do not have such a big influence on the result.... 
//Float_t weightsL_K0[Nbins] ={0.35,0.55,0.85}; //these numbers are the ratios L+antiL/K0 in pp- should be rechecked
//Float_t weightsL_K0[Nbins] ={0.92,0.30,0.50}; //these numbers are ROUGH estimation for the ratios L+antiL/K0 in PbPb peripheral 60-80 - should be rechecked
Float_t weightsL_K0[Nbins] ={1,0.20,0.4}; //these numbers are ROUGH estimation for the ratios L+antiL/K0 in PbPb central 0-5- should be rechecked


//File names for data and MC
TString fileMC, fileRD;

//  cutMode is the set of cuts contained in the files, it can be :
// "cutMode7000" : TPC only
// "cutMode8000" : TPC+ITS
TString cutMode="cutMode8000";

//Marker colors for the different templates
Int_t color[maxTemp]={2,3,4,6,7};
//histograms taken from the input files (data and MC) 
TH1D *hDCAy_pt, *hDCAyMCPrim_pt, *hDCAyMCSecDecays_pt, *hDCAyMCSecMaterial_pt, *hDCAyMCSecDecaysK0s_pt,*hDCAyMCSecDecaysLambda_pt;

//histograms for the templates and the data used in the fit
TH1D *htemplates[maxTemp], *hdata;
//histograms containing the adjusted MC template  from the fit (TFractionFitter only)
TH1D *htemplatesPred[maxTemp];

//function used to define the DCA range for different pT (x), case of TPC + ITS cuts - cutMode8000
TF1 *fDCA = new TF1("fDCA", "0.0182+0.035/x",0.15,200);


void Fit(TString outputname="", const Int_t nfrac = 3, const Float_t cenMin=0, const Float_t cenMax=-1, TString fileRD= "/lustre/nyx/alice/users/fsozzi/train/V012.PbPb2015/2016-04-12_1615.23621/mergedPeriods/PbPb/5.023ATeV/LHC15o.pass_lowint_firstphys/AnalysisResults.root", TString fileMC="/lustre/nyx/alice/users/fsozzi/train/V012.MC_PbPb2015/2016-04-12_1609.23621/mergedPeriods/MC_PbPb/5.023ATeV/LHC15k1a1/AnalysisResults.root" );

TH1D *performFit(std::vector<Double_t> &res, std::vector<Double_t> &errRes, Float_t &csq);
TH1D *performRooFit(std::vector<Double_t> &res, std::vector<Double_t> &errRes, Float_t &csq);

TH3D* GetDataHistoFromFile( TString fileRD, const Float_t cenMax, TString cutMode);
TH3D* GetEtaPtMCPrimHistoFromFile(TString fileMC, const Float_t cenMax,TString cutMode);
TH3D* GetEtaPtMCSecDecaysHistoFromFile(TString fileMC, const Float_t cenMax,TString cutMode);
TH3D* GetEtaPtMCSecMaterialHistoFromFile(TString fileMC, const Float_t cenMax,TString cutMode); 
TH3D* GetEtaPtMCSecDecaysK0sHistoFromFile(TString fileMC, const Float_t cenMax,TString cutMode); 
TH3D* GetEtaPtMCSecDecaysLambdaHistoFromFile(TString fileMC, const Float_t cenMax,TString cutMode); 
 
void SetCutMode(TString cut){cutMode=cut;}
TString GetCutMode(){return cutMode;}
void SetTemplatesNumber(Int_t num){nTemp=num;}
Int_t GetTemplatesNumber(){return nTemp;}
void SetData(TH1D *h){hdata=h;}
void buildTemplates(Int_t fitType, Int_t bin);


void Examples()
{

	// This function is not intended to be called:
	// it contains examples of lines to be copied in a ROOT interactive shell
	// to perform different analysis.
	// Modify the arguments in the "Fit()" call to perform different type of fits, arguments are explained in the Fit() function below.

	
	// copy the lines correspoding to your analysis
	//
	
  	// ------------------------------------------
	//           Analysis on pp data 
	// ------------------------------------------
    //NB original files for pp are on hera, copied on nyx
	// Fit pp, 2 templates, 8000 cuts ( without golden chi2 cut)
	fileRD = "/lustre/nyx/alice/users/fsozzi/train/V012.pp/2016-04-11_1705.23620/mergedPeriods/pp/5TeV/LHC15n.pass1/AnalysisResults.root";//full statistics 8000 cuts 2nd version, without golden chi2 cut
	fileMC = "/lustre/nyx/alice/users/fsozzi/train/V012.MC_pp/2016-04-11_2214.23620/mergedPeriods/MC_pp/5TeV/LHC15l1b2/AnalysisResults.root";//full statistics 8000 cuts 2nd version, without golden chi2 cut
//	fileRD = "/lustre/nyx/alice/users/fsozzi/train/V012.pp/2016-04-11_1705.23620/mergedRuns/pp/5TeV/LHC15n.pass1/244351.ana/AnalysisResults.root";
//	fileMC = "/lustre/nyx/alice/users/fsozzi/train/V012.MC_pp/2016-04-11_2214.23620/mergedRuns/MC_pp/5TeV/LHC15l1a2//AnalysisResults.root";
	Fit("",2,0,-1,fileRD,fileMC);


    // Fit pp, 2 templates, 7000 cuts
	fileRD = "/lustre/nyx/alice/users/fsozzi/train/V012.pp/2016-06-01_1418.23713/mergedPeriods/pp/5TeV/LHC15n.pass1/AnalysisResults.root";//full statistics 7000 cuts
	fileMC = "/lustre/nyx/alice/users/fsozzi/train/V012.MC_pp/2016-06-01_1805.23713/mergedPeriods/MC_pp/5TeV/LHC15l1a2/AnalysisResults.root";//7000 cuts
	SetCutMode("cutMode7000");
	Fit("",2,0,-1,fileRD,fileMC);

  	// ------------------------------------------
	//           Analysis on PbPb data 
	// ------------------------------------------
	
   
 	// Fit PbPb, 2 templates, 8000 cuts (without golden chi2 cut), all centrality ranges
	//fileRD = "/lustre/nyx/alice/users/fsozzi/train/V012.PbPb2015/2016-04-12_1615.23621/mergedPeriods/PbPb/5.023ATeV/LHC15o.pass_lowint_firstphys/AnalysisResults.root";//pass1_1stphysics
	fileRD = "/lustre/nyx/alice/users/fsozzi/train/V012.PbPb2015/2016-05-27_1554.23708/mergedPeriods/LHC15o.pass2/AnalysisResults.root";//pass2_lowIR - full statistics 8000 cuts 2nd version, without golden chi2 cut
	fileMC = "/lustre/nyx/alice/users/fsozzi/train/V012.MC_PbPb2015/2016-04-12_1609.23621/mergedPeriods/MC_PbPb/5.023ATeV/LHC15k1a1/AnalysisResults.root";
  	Fit("",2,0,-1,fileRD,fileMC);

	//Fit PbPb, 2 templates, 7000 cuts,  all centrality	
	//fileRD = "/lustre/nyx/alice/users/fsozzi/train/V012.PbPb2015/2016-03-23_2326.23570/mergedPeriods/PbPb/5.023ATeV/LHC15o.pass_lowint_firstphys/AnalysisResults.root";//pass1_1stphysics
	fileRD = "/lustre/nyx/alice/users/fsozzi/train/V012.PbPb2015/2016-05-26_1757.23708/mergedPeriods/LHC15o.pass2/AnalysisResults.root";//all run, centrality analysis, pass2!!
	fileMC = "/lustre/nyx/alice/users/fsozzi/train/V012.MC_PbPb2015/2016-03-23_2330.23570/mergedPeriods/MC_PbPb/5.023ATeV/LHC15k1a1/AnalysisResults.root";//full statistics 7000 cuts, centrality analysis
	SetCutMode("cutMode7000");
	Fit("",3,0,-1,fileRD,fileMC);



	
}

void buildTemplates(Int_t fitType, Int_t bin){
	
	// this functions builds the different possible templates
	// To use it, call the main function Fit setting the argument nfrac to one of these numbers
	
	switch(fitType){
	case 1: // 2 templates: (Primaries and sec from material) + sec from decay
		htemplates[0] = hDCAyMCPrim_pt;
		htemplates[0] ->Add(hDCAyMCSecMaterial_pt);  htemplates[0]->SetName("Prim+Sec material");
		htemplates[1] = hDCAyMCSecDecays_pt;         htemplates[1]->SetName("Sec from decay");
		break;
	case 2:  // 2 templates: Primaries + (sec from decay and material)
		htemplates[0] = hDCAyMCPrim_pt;                htemplates[0]->SetName("Primaries");
		htemplates[1] = hDCAyMCSecDecays_pt;
		htemplates[1] ->Add(hDCAyMCSecMaterial_pt);    htemplates[1]->SetName("Secondaries");
		break;
	case 3:  // 3 templates: primaries + secondaries + sec from material 
		htemplates[0]=hDCAyMCPrim_pt;        htemplates[0]->SetName("Primaries");
		htemplates[1]=hDCAyMCSecDecays_pt;   htemplates[1]->SetName("Sec from decay");
		htemplates[2]=hDCAyMCSecMaterial_pt; htemplates[2]->SetName("Sec from mat");
		break;
	case 4:  // 4 templates: prim + (k0 and L) + sec from mat+ (sec from decay!=k0,L)  
		cout<<"!!!!!!!!!!!!!!!!!!!! "<<endl;
		cout<< "please check that you selected the correct weight factor for your sample, in variable weightsL_K0"<<endl;
		cout<<"!!!!!!!!!!!!!!!!!!!! "<<endl;
		usleep(3000000);//sleep for 3 seconds to read the warning message :)

		htemplates[0] = hDCAyMCPrim_pt;                 htemplates[0]->SetName("Primaries");
		htemplates[1] = (TH1D *)hDCAyMCSecDecaysK0s_pt->Clone();
		htemplates[1] ->Add(hDCAyMCSecDecaysLambda_pt, (hDCAyMCSecDecaysK0s_pt->GetEntries()/hDCAyMCSecDecaysLambda_pt->GetEntries()) *weightsL_K0[bin]); htemplates[1]->SetName("K0+w*#Lambda");
		htemplates[2] = hDCAyMCSecMaterial_pt;          htemplates[2]->SetName("Sec from mat");
		htemplates[3] = hDCAyMCSecDecays_pt;
		htemplates[3] ->Add(hDCAyMCSecDecaysK0s_pt,-1);
		htemplates[3] ->Add(hDCAyMCSecDecaysLambda_pt,-1);                htemplates[3]->SetName("Sec != K0,#Lambda"); 
		break;
	case 5:  // 5 templates
		htemplates[0] = hDCAyMCPrim_pt;            htemplates[0]->SetName("Primaries");
		htemplates[1] = hDCAyMCSecDecaysK0s_pt;    htemplates[1]->SetName("K0");
		htemplates[2] = hDCAyMCSecDecaysLambda_pt; htemplates[2]->SetName("#Lambda");
		htemplates[3] = hDCAyMCSecDecays_pt;
		htemplates[3] ->Add(htemplates[2],-1);     
		htemplates[3] ->Add(htemplates[1],-1);     htemplates[3]->SetName("Sec != K0,#Lambda");
		htemplates[4] = hDCAyMCSecMaterial_pt;     htemplates[4]->SetName("Sec from mat");
		break;
	case 6:  // 3 templates: (prim+ sec from mat) + (k0 and L) + (sec from decay!=k0,L)  
		htemplates[0] = hDCAyMCPrim_pt;                
		htemplates[0] ->Add(hDCAyMCSecMaterial_pt);  htemplates[0]->SetName("Prim+Sec material");
		htemplates[1] = (TH1D *)hDCAyMCSecDecaysK0s_pt->Clone();
		htemplates[1] ->Add(hDCAyMCSecDecaysLambda_pt, (hDCAyMCSecDecaysK0s_pt->GetEntries()/hDCAyMCSecDecaysLambda_pt->GetEntries()) *weightsL_K0[bin]); htemplates[1]->SetName("K0+w*#Lambda");
		htemplates[2] = hDCAyMCSecDecays_pt;
		htemplates[2] ->Add(hDCAyMCSecDecaysK0s_pt,-1);
		htemplates[2] ->Add(hDCAyMCSecDecaysLambda_pt,-1);                htemplates[2]->SetName("Sec != K0,#Lambda"); 
		break;

		//You can add your own definition of templates here, adding a "case".
		//Please be aware that dependening on your definition,
		//you could have to change some lines in the Fit function also, to draw exceptions (e.g. proper normalizazion of the fractions)
		
	default: //like case 3 
		htemplates[0]=hDCAyMCPrim_pt;htemplates[0]->SetName("Primaries");
		htemplates[1]=hDCAyMCSecDecays_pt;htemplates[1]->SetName("Sec from decay");
		htemplates[2]=hDCAyMCSecMaterial_pt;htemplates[2]->SetName("Sec from mat");
		break;
	}

	//define a color for each template
	for (Int_t i = 0; i < nTemp; i++) {
		htemplates[i]->SetLineColor(color[i]);
	}

		
}


void Fit(TString outputname, const Int_t nfrac, const Float_t cenMin, const Float_t cenMax, TString fileRD, TString fileMC)
{

	//Arguments : 
    //	outputname: suffix for the canvases name 
	//  nfrac is used to set the number and type of templates in the fit. It should be a number between 1 and 6 - see function buildTemplates to choose the fit  you want

	//   
    //  cenMin, cenMax: set if you want to select a centrality range (range goes from 0 to 100)
	// fileRD and fileMC are the name of the input files, respectively real data and MC


	
	//check that nfrac makes sense
	if(nfrac<1 || nfrac>6)
	{
		cout<<"nfrac should be between 1 and 6 "<<endl;
		return;
	}
    //	set the templates for the fit according to the type of fit chosen
	if(nfrac<3) SetTemplatesNumber(2);
	if(nfrac==3||nfrac==6) SetTemplatesNumber(3);
	if(nfrac==4) SetTemplatesNumber(4);
	if(nfrac==5) SetTemplatesNumber(5);

    //Some style stuff
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetLabelSize(0.05,"x");
	gStyle->SetLabelSize(0.05,"Z");
	gStyle->SetTitleSize(0.05,"X");
	gStyle->SetTitleSize(0.05,"Y");
	gStyle->SetOptStat(false);

	//histogram that will contain the ratio between the fit results and the data
	TH1D *ratio;
	
	//graph containing the fractions (result, initial conditions,  fraction from MC, result and fraction from MC in DCA range)
	TGraphAsymmErrors *gr[nTemp]; TGraph *grIn[nTemp]; TGraph *grFr[nTemp]; TGraph *grCorr[nTemp]; TGraph *grFrCorr[nTemp]; 
	for (Int_t i = 0; i < nTemp; i++) {
		gr[i]=new TGraphAsymmErrors(Nbins);
		gr[i]->SetMarkerStyle(20+i);
		gr[i]->SetMarkerColor(color[i]);
		gr[i]->SetLineColor(color[i]);
		grIn[i]=new TGraph(Nbins);
		grIn[i]->SetMarkerStyle(24+i);
		grIn[i]->SetMarkerColor(color[i]);
		grIn[i]->SetLineColor(color[i]);

		grFr[i]=new TGraph(Nbins);
		grFr[i]->SetMarkerStyle(31);
		grFr[i]->SetMarkerColor(color[i]);
		grFr[i]->SetLineColor(color[i]);

		grFrCorr[i]=new TGraph(Nbins);
		grFrCorr[i]->SetMarkerStyle(31);
		grFrCorr[i]->SetMarkerColor(color[i]);
		grFrCorr[i]->SetLineColor(color[i]);

		grCorr[i]=new TGraph(Nbins);
		grCorr[i]->SetMarkerStyle(31);
		grCorr[i]->SetMarkerColor(color[i]);
		grCorr[i]->SetLineColor(color[i]);
	}


	
// get histo from data 
	TH3D *fDCAyEtaPt = GetDataHistoFromFile(fileRD,cenMax,cutMode);

// get histo from MC 
	TH3D *fDCAyEtaPtMCPrim = GetEtaPtMCPrimHistoFromFile( fileMC,  cenMax, cutMode);
	TH3D *fDCAyEtaPtMCSecDecays = GetEtaPtMCSecDecaysHistoFromFile( fileMC,  cenMax, cutMode);
	TH3D *fDCAyEtaPtMCSecMaterial = GetEtaPtMCSecMaterialHistoFromFile( fileMC,  cenMax, cutMode);

	TH3D *fDCAyEtaPtMCSecDecaysK0s = GetEtaPtMCSecDecaysK0sHistoFromFile( fileMC,  cenMax, cutMode);
	TH3D *fDCAyEtaPtMCSecDecaysLambda = GetEtaPtMCSecDecaysLambdaHistoFromFile( fileMC,  cenMax, cutMode);
	
	//prepare canvases
	Int_t nraws;
	if(Nbins%3==0) nraws = Nbins/3;
	else nraws = 1+Nbins/3;
	//canvas for templates
	TCanvas *can = new TCanvas("can", "can",13,33,983,250*nraws);
	can->Divide(3,nraws);
	//canvas for ratios
	TCanvas *canratio = new TCanvas("canratio", "canratio",13,33,983,250*nraws);
	canratio->Divide(3,nraws);

    // variables to store the fit results 
	std::vector<Double_t> res, errRes;


	Float_t scalingValue[Nbins], scalingValueDCA[Nbins]; //these scaling factors are needed to properly evaluate the fraction of primaries, in case the template uses primaries + sec from material (nfrac==1or6)
	 
	//translate centrality value in  bin number:
	Float_t binCenMin=	fDCAyEtaPt->GetYaxis()->FindBin(cenMin);
	Float_t binCenMax=	fDCAyEtaPt->GetYaxis()->FindBin(cenMax-0.000001);
	cout<<cenMin<<" "<<cenMax<<endl;
	cout<<binCenMin<<" "<<binCenMax<<endl;
	
    //Loop over pt bins
	for (Int_t i=0; i<Nbins; i++){

		//translate pt value in pt bin number:
		Float_t binMin=	fDCAyEtaPt->GetZaxis()->FindBin(binsPt[i]);
		Float_t binMax=	fDCAyEtaPt->GetZaxis()->FindBin(binsPt[i+1]-0.000001);
		fDCAyEtaPt->GetZaxis()->SetRange(binMin,binMax);
		Float_t meanPt = fDCAyEtaPt->GetMean(3); 
   
		cout<<"Performing fit num  "<<i<<", pt bin ["<<binsPt[i]<<" "<<binsPt[i+1]<<" ]"<<endl;
		cout<<"Mean pt is : "<<meanPt<<endl;
		cout<<"corresponding to bin num."<< binMin<<" "<<binMax<<endl;
		cout<<"-----------------------------------"<<endl;

		
	
        //make projections to the DCA variable for each pt (and centrality) bin
		//these are the histograms used for the fit
		hDCAy_pt = (TH1D*)fDCAyEtaPt->ProjectionX(Form("hDCAy_pt%d",i+1),binCenMin,binCenMax,binMin, binMax);

		hDCAyMCPrim_pt = (TH1D*)fDCAyEtaPtMCPrim->ProjectionX(Form("hDCAyMCPrim_pt%d",i+1),binCenMin,binCenMax,binMin, binMax);

		hDCAyMCSecDecays_pt = (TH1D*)fDCAyEtaPtMCSecDecays->ProjectionX(Form("hDCAyMCSecDecays_pt%d",i+1),binCenMin,binCenMax,binMin, binMax);

		hDCAyMCSecMaterial_pt = (TH1D*)fDCAyEtaPtMCSecMaterial->ProjectionX(Form("hDCAyMCSecMaterial_pt%d",i+1),binCenMin,binCenMax,binMin, binMax);

		hDCAyMCSecDecaysK0s_pt = (TH1D*)fDCAyEtaPtMCSecDecaysK0s->ProjectionX(Form("hDCAyMCSecDecaysK0s_pt%d",i+1),binCenMin,binCenMax,binMin, binMax);

		hDCAyMCSecDecaysLambda_pt = (TH1D*)fDCAyEtaPtMCSecDecaysLambda->ProjectionX(Form("hDCAyMCSecLambda_pt%d",i+1),binCenMin,binCenMax,binMin, binMax);

		
		//Rebin in case option was set to true
		if(rebin){
			hDCAy_pt->Rebin(2);
			hDCAyMCPrim_pt->Rebin(2);
			hDCAyMCSecDecays_pt->Rebin(2);
			hDCAyMCSecMaterial_pt->Rebin(2);
			hDCAyMCSecDecaysK0s_pt->Rebin(2);
			hDCAyMCSecDecaysLambda_pt->Rebin(2);
		}

		
	
		
		// Evaluate the DCA range used in the analysis using the formula and the pt mean value in the bin - this is used only in case the "TPC+ITS"cuts are used, namely "8000cuts"
		Float_t rangeDCAmin=1,rangeDCAmax=100;
		if(GetCutMode().Contains("8000")){
		 rangeDCAmin = hDCAyMCPrim_pt->GetXaxis()->FindBin(-fDCA->Eval(meanPt)+0.000001);
		 rangeDCAmax = hDCAyMCPrim_pt->GetXaxis()->FindBin(fDCA->Eval(meanPt)-0.000001);
		}
		else //get the range of the histogram
		{
			rangeDCAmax = hDCAyMCPrim_pt->GetNbinsX();
		}
		//	cout<<rangeDCAmin<<" "<<rangeDCAmax<<endl;
		
		//--> evaluate the fractions in the MC -- NB these are not used anymore - 
		cout<<"Statistics of the histograms: "<<hDCAy_pt->Integral()<<" "<<hDCAyMCPrim_pt->Integral()<<" "<<hDCAyMCSecDecays_pt->Integral()<<" "<<hDCAyMCSecMaterial_pt->Integral()<<" "<<endl;
		Double_t tot = hDCAyMCPrim_pt->Integral()+hDCAyMCSecDecays_pt->Integral()+hDCAyMCSecMaterial_pt->Integral();
				
		Float_t totDCArange = hDCAyMCPrim_pt->Integral(rangeDCAmin,rangeDCAmax)+hDCAyMCSecDecays_pt->Integral(rangeDCAmin,rangeDCAmax)+hDCAyMCSecMaterial_pt->Integral(rangeDCAmin,rangeDCAmax);

		//these are the fractions from pure MC
		frac[0] = hDCAyMCPrim_pt->Integral()/tot;// this contains the sum of primaries and Material if nfrac == 1
		frac[1] = hDCAyMCSecDecays_pt->Integral()/tot; // this contains the sum of MCDecay and Material if nfrac == 2
		frac[2] = hDCAyMCSecMaterial_pt->Integral()/tot; 

		//these are the fractions from pure MC in the smaller DCA range 
		fracCorr[0] = hDCAyMCPrim_pt->Integral(rangeDCAmin,rangeDCAmax)/totDCArange;// this contains the sum of primaries and Material if nfrac == 1
		fracCorr[1] = hDCAyMCSecDecays_pt->Integral(rangeDCAmin,rangeDCAmax)/totDCArange; // this contains the sum of MCDecay and Material if nfrac == 2
		fracCorr[2] = hDCAyMCSecMaterial_pt->Integral(rangeDCAmin,rangeDCAmax)/totDCArange; 

		//<-----------------
		

		//in case of fits 1 and 6, the primaries are summed with secondaries from material, therefore a scaling factor to evaluate the actual primaries is needed
		if(nfrac==1||nfrac==6){
			
			scalingValue[i]=hDCAyMCPrim_pt->Integral()/(hDCAyMCPrim_pt->Integral()+hDCAyMCSecMaterial_pt->Integral());
			scalingValueDCA[i]=hDCAyMCPrim_pt->Integral(rangeDCAmin,rangeDCAmax)/(hDCAyMCPrim_pt->Integral(rangeDCAmin,rangeDCAmax)+hDCAyMCSecMaterial_pt->Integral(rangeDCAmin,rangeDCAmax));

		}


		//set the histogram for the data
		SetData(hDCAy_pt);
		//Build the MC templates from the histograms
		buildTemplates(nfrac,i);
		 
		//evaluate the fractions from the templates, they can be used as initial value in the fit  
		tot=0;
		totDCArange=0;
		for (Int_t k = 0; k<nTemp ; k++){
			tot+=htemplates[k]->Integral();
			totDCArange+=htemplates[k]->Integral(rangeDCAmin,rangeDCAmax);
		}
		for (Int_t k = 0; k<nTemp ; k++){
			fracTemp[k]=htemplates[k]->Integral()/tot;
			fracTempCorr[k]=htemplates[k]->Integral(rangeDCAmin,rangeDCAmax)/totDCArange;
		}
		if(useInFrac){	  
			for (Int_t k = 0; k<nTemp ; k++){
				val[k]=fracTemp[k];
			}
		}

		//perform the fit - res and errRes will contain the fractions and their errors; result is the fitted model
		Float_t csq=0;
		TH1D *result=0;
		if(rooFit)
			result = performRooFit(res, errRes, csq);
		else
			result = performFit(res, errRes, csq);

		
		// for (Int_t ind = 0; ind < nTemp; ind++) {
		// 	cout<<res[ind]<<" "<<errRes[ind]<<" "<<csq<<endl;
		// }

		
//prepare graphs with the results
		for (Int_t ind = 0; ind < nTemp; ind++) {
			gr[ind]->SetPoint(i,meanPt,res[ind]);
			gr[ind]->SetPointError(i,-binsPt[i]+meanPt,binsPt[i+1]-meanPt,errRes[ind]*0.5,errRes[ind]*0.5);
			
			//and those with the initial conditions,  
			grIn[ind]->SetPoint(i,meanPt,val[ind]);

            //the MC fractions,
			grFr[ind]->SetPoint(i,meanPt,fracTemp[ind]);

			// the MC fractions in the smaller DCA range
			grFrCorr[ind]->SetPoint(i,meanPt,fracTempCorr[ind]);
		}

		// Evaluate the fraction in a DCA subrange :
		// - Evaluate the normalization for the fraction as: Sum_j (frac_j * I_j(R) / I_j)
		Float_t norm=0;
		for (Int_t ind = 0; ind < nTemp; ind++) {
		
			norm += res[ind]*htemplates[ind]->Integral(rangeDCAmin,rangeDCAmax)/htemplates[ind]->Integral();
		}
		// - Fill the graph with the corrected values evaluated as : (frac_j * I_j(R) / I_j)/Sum_j (frac_j * I_j(R) / I_j)
		for (Int_t ind = 0; ind < nTemp; ind++) {
			grCorr[ind]->SetPoint(i,meanPt,(1/norm)*res[ind]*htemplates[ind]->Integral(rangeDCAmin,rangeDCAmax)/htemplates[ind]->Integral());
		}


//Draw distributions
		can->cd(i+1);
		gPad->SetLogy();
		hdata->SetMinimum(1);//minimum set to 1 for display reason in log scale
		hdata->SetTitle("");
		hdata->GetXaxis()->SetLabelSize(0.05);
		hdata->GetXaxis()->SetTitleSize(0.05);
		hdata->Draw("");
		if (result!= 0x0)
		{
			result->SetMarkerStyle(8);
			result->SetMarkerSize(0.4);
			//	result->Draw("epsame");

		}
		for (Int_t ind = 0; ind < nTemp; ind++) {
			htemplates[ind]->Draw("same");
	// 		htemplatesPred[ind]->SetLineStyle(2);
	// htemplatesPred[ind]->Draw("same");
		}

		TLatex tex(-0.1,1,"");
		tex.SetTextFont(32);
		tex.SetTextSize(0.08);
		tex.DrawLatex(hdata->GetBinLowEdge(2),hdata->GetMaximum()*0.1,TString::Format("%g<p_{T}<%g GeV/c",binsPt[i],binsPt[i+1]).Data());
  
		// Draw the ratio of results and histogram
		canratio->cd(i+1); 
		gPad->SetGridy();

		// To display errors
		// hdata[i]->Sumw2();
		// result->Sumw2();

		TH1D *histProj=new TH1D("proj","proj",100,0.8,1.2);
		if (result!= 0x0){
			//	cout<<"........"<<endl;
			ratio = (TH1D *)result->Clone("");
			ratio ->SetTitle(" ");
			ratio ->Divide(hdata);
			ratio ->GetXaxis()->SetLabelSize(0.05);
			ratio ->GetXaxis()->SetTitleSize(0.05);
			ratio ->GetYaxis()->SetRangeUser(0.9,1.1);
			ratio ->Print();
			//	ratio ->Draw();
			for (Int_t k=1; k < ratio ->GetNbinsX(); k++) {
				histProj->Fill(ratio ->GetBinContent(k));
				cout<<ratio ->GetBinContent(k)<<endl;
			}
			histProj->Draw();
			histProj->Draw();
			tex.DrawLatex(histProj->GetBinLowEdge(0.85),1.08,TString::Format("%g", histProj->GetRMS()).Data());	
			
			//tex.DrawLatex(ratio->GetBinLowEdge(2),1.08,TString::Format("%g<p_{T}<%g GeV/c",binsPt[i],binsPt[i+1]).Data());
			//tex.DrawLatex(ratio->GetBinLowEdge(2),0.92,TString::Format("#chi^{2}=%g",csq).Data());
		}
	}
	canratio->Print(TString::Format("ratios%s.gif",outputname.Data()));
	can->Print(TString::Format("histograms%s.gif",outputname.Data()));
	can->Print(TString::Format("histograms%s.eps",outputname.Data()));
	canratio->Print(TString::Format("ratios%s.pdf",outputname.Data()));
	can->Print(TString::Format("histograms%s.pdf",outputname.Data()));

	//Draw the results
	TCanvas *can2 = new TCanvas("can2","can2");
	can2->cd();can2->SetLogx();
	gr[0]->SetTitle(" ");
	gr[0]->GetXaxis()->SetRangeUser(0.1,binsPt[Nbins]);
	gr[0]->GetYaxis()->SetRangeUser(0,1);
	gr[0]->GetYaxis()->SetTitle("fractions");
	gr[0]->GetXaxis()->SetTitle("p_{T} (GeV/c)         ");
	gr[0]->Draw("ap");

	TLine *l = new TLine(0,0,0,0);
	TLegend *leg = new TLegend(0.5,0.4,0.68,0.6);
	leg->SetHeader("Color Legend");
	for (Int_t ind = 0; ind < nTemp; ind++) {
		gr[ind]->Draw("p");
		grIn[ind]->Draw("p");
		grFr[ind]->Draw("p");
		l= new TLine(binsPt[0],minV[ind],binsPt[Nbins],minV[ind]);
		l->SetLineStyle(2);
		l->SetLineColor(color[ind]);
		if(!rooFit)	l->Draw("same");
		l = new TLine(binsPt[0],maxV[ind],binsPt[Nbins],maxV[ind]);
		l->SetLineStyle(2);
		l->SetLineColor(color[ind]);
		if(!rooFit)	l->Draw("same");

		leg->AddEntry(gr[ind],htemplates[ind]->GetName(),"lep");
	}


	leg->DrawClone();
	leg = new TLegend(0.68,0.4,0.9,0.6);
	leg->SetHeader("Type legend");
	leg->AddEntry(gr[nTemp-1],"Fit value","lep");
	leg->AddEntry(grIn[nTemp-1]," Initial value ","p");
	leg->AddEntry(grFr[nTemp-1]," Fraction from MC ","p");
	leg->AddEntry(l,"Boundaries ","l");
	leg->Draw();


	can2->Print(TString::Format("results%s.gif",outputname.Data()));
	can2->Print(TString::Format("results%s.pdf",outputname.Data()));


	//Final printing - tables with <pt> and fractions
	cout<<"MC fractions "<<setprecision(3)<<endl;
	for (Int_t k = 0; k < Nbins; k++) {
		cout<<grFr[0]->GetX()[k]<<"\t";
		for (Int_t i = 0; i < nTemp; i++) {
			cout<<grFr[i]->GetY()[k]<<"\t";
		}
		cout<<" "<<endl;
	}

	cout<<" MC fractions  in smaller DCA range "<<setprecision(3)<<endl;
	for (Int_t k = 0; k < Nbins; k++) {
		cout<<grFrCorr[0]->GetX()[k]<<"\t";
		for (Int_t i = 0; i < nTemp; i++) {
			cout<<grFrCorr[i]->GetY()[k]<<"\t";
		}
		cout<<" "<<endl;
	}
	cout<<"Results "<<setprecision(3)<<endl;
	for (Int_t k = 0; k < Nbins; k++) {
		cout<<gr[0]->GetX()[k]<<"\t";
		//	cout<<"XXXX \t";
		for (Int_t i = 0; i < nTemp; i++) {
			//		cout<<gr[i]->GetY()[k]<<"\pm"<<gr[i]->GetEY()[k]<<"\t";
			cout<<gr[i]->GetY()[k]<<"\t";
		}
		cout<<" "<<endl;
	}
// cout<<"Results "<<setprecision(3)<<endl;
// 		//cout<<gr[0]->GetX()[k]<<"\t";
// 		for (Int_t i = 0; i < nTemp; i++) {
// 			cout<<"XXXX \t";
// 			for (Int_t k = 0; k < Nbins; k++) {
// 	//		cout<<gr[i]->GetY()[k]<<"\pm"<<gr[i]->GetEY()[k]<<"\t";
// 				cout<<gr[i]->GetY()[k]<<"\t";
// 		}
// 		cout<<" "<<endl;
// 	}

	cout<<" Results in smaller DCA range "<<setprecision(3)<<endl;
	for (Int_t k = 0; k < Nbins; k++) {
		cout<<"  "<<grCorr[0]->GetX()[k]<<"\t";
		for (Int_t i = 0; i < nTemp; i++) {
			//		cout<<gr[i]->GetY()[k]<<"\pm"<<gr[i]->GetEY()[k]<<"\t";
			cout<<grCorr[i]->GetY()[k]<<"\t";
	}
		cout<<" "<<endl;
	}
	cout<<" Results in smaller DCA range "<<setprecision(3)<<endl;
	for (Int_t k = 0; k < Nbins; k++) {
		cout<<"XXXX "<<grCorr[0]->GetX()[k]<<"\t";
		for (Int_t i = 0; i < nTemp; i++) {
			//		cout<<gr[i]->GetY()[k]<<"\pm"<<gr[i]->GetEY()[k]<<"\t";
			cout<<grCorr[i]->GetY()[k]<<"\t";
			cout<<grFrCorr[i]->GetY()[k]<<"\t";
	}
		cout<<" "<<endl;
	}	
	if(nfrac==1||nfrac==6){
		cout<<" Results for template 1 scaled to obtain the real primaries "<<setprecision(3)<<endl;
		cout<<" all range "<<setprecision(3)<<endl;
		for (Int_t k = 0; k < Nbins; k++) {
			cout<<scalingValue[k]*gr[0]->GetY()[k]<<"\t";
			cout<<" "<<endl;
		}
		cout<<" DCA range "<<setprecision(3)<<endl;
		for (Int_t k = 0; k < Nbins; k++) {
			cout<<scalingValueDCA[k]*grCorr[0]->GetY()[k]<<"\t";
			cout<<" "<<endl;
		}
	}
}


TH3D* GetDataHistoFromFile( TString fileRD, const Float_t cenMax, TString cutMode){
	TFile *file = TFile::Open(fileRD);
	
	if(!file)
	{
		cout<< "GetDataHistoFromFile did not find the file " <<endl;
		return 0;
	}
	file->cd("dNdPtCutAnalysis");

// get object
	AlidNdPtCutAnalysis *obj = (AlidNdPtCutAnalysis *)gDirectory->Get(cutMode)->FindObject("dNdPtCutAnalysisObj");
	if(!obj){cout<<"file not containing object dNdPtCutAnalysisObj"<<endl; return 0 ; }
	TH3D *fDCAyEtaPt;
	if(cenMax==-1)//no centrality range is selected, therefore take the original histo with eta 
		fDCAyEtaPt = (TH3D*)obj->GetDCAyEtaPt();
	else // otherwise take the one in bins of centrality
	 	fDCAyEtaPt = (TH3D*)obj->GetDCAyCenPt();


	file->Close();
	return (fDCAyEtaPt);
}


TH3D* GetEtaPtMCPrimHistoFromFile(TString fileMC, const Float_t cenMax,TString cutMode){ 
//
// open mc file 
	TFile *file1 = TFile::Open(fileMC);
	if(!file1) {
		cout<< "GetEtaPtMCPrimHistoFromFile did not find the file " <<endl;
		return 0;
	}
	file1->cd("dNdPtCutAnalysis");

	AlidNdPtCutAnalysis *objmc = (AlidNdPtCutAnalysis *)gDirectory->Get(cutMode)->FindObject("dNdPtCutAnalysisObj");
	if(!objmc){cout<<"file not containing object dNdPtCutAnalysisObj"<<endl; return 0 ; }
	TH3D *fDCAyEtaPtMCPrim;
	if(cenMax==-1)//no centrality range is selected, therefore take the original histo with eta 
	{
		fDCAyEtaPtMCPrim = (TH3D*)objmc->GetDCAyEtaPtMCPrim();
	}
	else // otherwise take the one in bins of centrality
	{
	 	fDCAyEtaPtMCPrim = (TH3D*)objmc->GetDCAyCenPtMCPrim();
	}
	file1->Close();
	return fDCAyEtaPtMCPrim;
}

TH3D* GetEtaPtMCSecDecaysHistoFromFile(TString fileMC, const Float_t cenMax,TString cutMode){ 
//
// open mc file 
	TFile *file1 = TFile::Open(fileMC);
	if(!file1)
	{
		cout<< "GetEtaPtMCSecDecaysHistoFromFile did not find the file " <<endl;
		return 0;
	}
	file1->cd("dNdPtCutAnalysis");

	AlidNdPtCutAnalysis *objmc = (AlidNdPtCutAnalysis *)gDirectory->Get(cutMode)->FindObject("dNdPtCutAnalysisObj");
	if(!objmc){cout<<"file not containing object dNdPtCutAnalysisObj"<<endl; return 0 ; }
	TH3D *fDCAyEtaPtMCSecDecays;
	if(cenMax==-1)//no centrality range is selected, therefore take the original histo with eta 
	{
		fDCAyEtaPtMCSecDecays = (TH3D*)objmc->GetDCAyEtaPtMCSecDecays();
	}
	else // otherwise take the one in bins of centrality
	{
	 	fDCAyEtaPtMCSecDecays = (TH3D*)objmc->GetDCAyCenPtMCSecDecays();
	}
	file1->Close();
	return fDCAyEtaPtMCSecDecays;
}

TH3D* GetEtaPtMCSecMaterialHistoFromFile(TString fileMC, const Float_t cenMax,TString cutMode){ 
//
// open mc file 
	TFile *file1 = TFile::Open(fileMC);
	if(!file1) {
		cout<< "GetEtaPtMCSecMaterialHistoFromFile did not find the file " <<endl;
		return 0;
	}
	file1->cd("dNdPtCutAnalysis");

	AlidNdPtCutAnalysis *objmc = (AlidNdPtCutAnalysis *)gDirectory->Get(cutMode)->FindObject("dNdPtCutAnalysisObj");
	TH3D *fDCAyEtaPtMCSecMaterial;
	if(cenMax==-1)//no centrality range is selected, therefore take the original histo with eta 
	{
		fDCAyEtaPtMCSecMaterial = (TH3D*)objmc->GetDCAyEtaPtMCSecMaterial();
	}
	else // otherwise take the one in bins of centrality
	{
		fDCAyEtaPtMCSecMaterial = (TH3D*)objmc->GetDCAyCenPtMCSecMaterial();
	}
	file1->Close();
	return fDCAyEtaPtMCSecMaterial;
}

TH3D* GetEtaPtMCSecDecaysK0sHistoFromFile(TString fileMC, const Float_t cenMax,TString cutMode){ 
//
// open mc file
	
	TFile *file1 = TFile::Open(fileMC);
	if(!file1) {
		cout<< "GetEtaPtMCSecDecaysK0sHistoFromFile did not find the file " <<endl;
		return 0;
	}
	file1->cd("dNdPtCutAnalysis");

	AlidNdPtCutAnalysis *objmc = (AlidNdPtCutAnalysis *)gDirectory->Get(cutMode)->FindObject("dNdPtCutAnalysisObj");
	if(!objmc){cout<<"file not containing object dNdPtCutAnalysisObj"<<endl; return 0 ; }
	TH3D *fDCAyEtaPtMCSecDecaysK0s;
	if(cenMax==-1)//no centrality range is selected, therefore take the original histo with eta 
	{
		fDCAyEtaPtMCSecDecaysK0s = (TH3D*)objmc->GetDCAyEtaPtMCSecDecaysK0s();
	}
	else // otherwise take the one in bins of centrality
	{
		fDCAyEtaPtMCSecDecaysK0s = (TH3D*)objmc->GetDCAyCenPtMCSecDecaysK0s();
	}
	file1->Close();
	return fDCAyEtaPtMCSecDecaysK0s;
}


TH3D* GetEtaPtMCSecDecaysLambdaHistoFromFile(TString fileMC, const Float_t cenMax,TString cutMode){ 
//
// open mc file 
	TFile *file1 = TFile::Open(fileMC);
	if(!file1) {
		cout<< "GetEtaPtMCSecDecaysLambdaHistoFromFile did not find the file " <<endl;
		return 0;
	}
	file1->cd("dNdPtCutAnalysis");

	AlidNdPtCutAnalysis *objmc = (AlidNdPtCutAnalysis *)gDirectory->Get(cutMode)->FindObject("dNdPtCutAnalysisObj");
	if(!objmc){cout<<"file not containing object dNdPtCutAnalysisObj"<<endl; return 0 ; }
	TH3D *fDCAyEtaPtMCSecDecaysLambda;
	if(cenMax==-1)//no centrality range is selected, therefore take the original histo with eta 
	{
		fDCAyEtaPtMCSecDecaysLambda = (TH3D*)objmc->GetDCAyEtaPtMCSecDecaysLambda();
	}
	else // otherwise take the one in bins of centrality
	{
		fDCAyEtaPtMCSecDecaysLambda = (TH3D*)objmc->GetDCAyCenPtMCSecDecaysLambda();
	}
	file1->Close();
	return fDCAyEtaPtMCSecDecaysLambda;
}


TH1D* performFit(std::vector<Double_t> &res, std::vector<Double_t> &errRes, Float_t &csq) {

	// this function contains the template fit with TFractionFitter
	// it needs the global variables htemplates and hdata to be filled, the vector for the limits (minV, maxV) and for the initial values (val) of the parameters
	// Fill the vectors with the results and their errors (res, errRes)
	// return the fitted histogram
	
			
// MC histograms are put in this array
	TObjArray *mc = new TObjArray(nTemp);
	for (Int_t ind = 0; ind < nTemp; ind++) {
		mc->Add(htemplates[ind]);
	}
		 
// Fit data with MC templates
	TFractionFitter* fit = new TFractionFitter(hdata, mc);
	//fit->SetRangeX(50,250); // fit range in bins

	//exclude bins
	// for (j = 44; j < 56; j++) {
	// 	fit->ExcludeBin(j);
	// }

	//constrains and initial conditions on the parameters
	for (Int_t ind = 0; ind < nTemp; ind++) {
		fit->Constrain(1+ind, minV[ind], maxV[ind]);
	}
	TVirtualFitter* vFit = fit->GetFitter();
//		vFit->SetPrecision(0.1);//this seems not to have any effect ... why ? 
	for (Int_t ind = 0; ind < nTemp; ind++) {
		vFit->SetParameter(ind,htemplates[ind]->GetName(),val[ind],0.001,minV[ind], maxV[ind]);
	}

	//First fit
	Int_t status = fit->Fit(); 
	cout << "fit status: " << status << endl;
	TH1D* result =0;
	if (status == 0)  {
		result = (TH1D*) fit->GetPlot();
	}
	else {
		//second fit, to improve the cases with "status 4" 
		status = fit->Fit(); 
		cout << "fit status: " << status << endl;
		if (status == 0)  {
			result = (TH1D*) fit->GetPlot();
		}
	}
	// fit->UnConstrain(1);
	// fit->UnConstrain(2);
	// fit->UnConstrain(3);
	// 		status = fit->Fit(); 

	//get the results and put them in the variables f_i,e_i
	//Float_t csq=0;
	csq = (Float_t)fit->GetChisquare()/(Float_t)fit->GetNDF();
	cout<< " csq: "<<csq<<endl;
	for (Int_t ind = 0; ind < nTemp; ind++) {
		res.push_back(1);errRes.push_back(1);
		fit->GetResult(ind,res[ind],errRes[ind]);
		cout << " res =" <<  res[ind]<< " err=" << errRes[ind] << endl; 
	}

	//MINOS ERRORS ANALYSIS, in order to evaluate correct errors (due to the constraints)
	//same results as the parabolic errors
	// Double_t arglist[2];
	// arglist[0] = 500;
	// arglist[1] = 1e-4;
	// //	vFit->ExecuteCommand("MIGRAD",arglist,2);
	// vFit->ExecuteCommand("MINOS",arglist,0);
  
	for (Int_t ind = 0; ind < nTemp; ind++) {
			
		// Get also the "true" MC distribution, namely the templates used
		// Careful: always check that these histograms are reasonably similar to the initial ones - this is not always the case!
		htemplatesPred[ind]=(TH1D *)fit->GetMCPrediction(ind);
	}

	delete fit;
	delete mc;
	return result;

}

TH1D* performRooFit(std::vector<Double_t> &res, std::vector<Double_t> &errRes,  Float_t &csq) 
{
   	// this function contains the template fit RooFIT
	// it needs the global variables htemplates and hdata to be filled and the vector of the initial values (val)
	// Fill the vectors with the results and their errors (res, errRes)
	// return the fitted histogram
	//
	// This function uses "rigid templates" (templates histograms are used as they are - not good in case of low statistics);
     // Todo: implement a smoothing or something similar 

	//declare dca variable - limits taken from the actual histograms
	RooRealVar x("x","x",hdata->GetXaxis()->GetXmin(),hdata->GetXaxis()->GetXmax());

 	// Import data and MC histograms in RooFit variables
	RooDataHist rooFitData( hdata->GetName(),hdata->GetTitle(),x,hdata );
	std::vector<RooDataHist> dataHistTemplates;
	for (Int_t ind = 0; ind < nTemp; ind++) {
		dataHistTemplates.push_back(RooDataHist(TString::Format("h%d",ind),TString::Format("h%d",ind),x,htemplates[ind]));
	}


	//build the same for PDF objects
	std::vector<RooHistPdf> histPDFTemplates;
	for (Int_t ind = 0; ind < nTemp; ind++) {
		histPDFTemplates.push_back(RooHistPdf(TString::Format("hpdf%d",ind),TString::Format("hpdf%d",ind),x,dataHistTemplates[ind]));
	}

	//and put them in a container (NB: different loops because putting everything in one loop create weird objects - not clear why)
	RooArgList templatesPDFList;
	for (Int_t ind = 0; ind < nTemp; ind++) {
		templatesPDFList.add(histPDFTemplates[ind]);
	}

	
	//constrains and initial conditions on the parameters
	std::vector<RooRealVar> coefficients;
	for (Int_t ind = 0; ind < nTemp; ind++) {
	
		coefficients.push_back(RooRealVar(TString::Format("v%d",ind),TString::Format("v%d",ind), val[ind],0,1)) ;
	}

// coefficient list  must be shorter of 1 unit wrt the template list (the last coefficient is constrained to be 1 - Sum (others )
	RooArgList fractionList;
	for (Int_t ind = 0; ind < nTemp-1; ind++) {
		fractionList.add(coefficients[ind]);
	}
	

	//construct the recursive (last option set to kTRUE) sum of the single PDFs 
	RooAddPdf model2("model2","model2",templatesPDFList,fractionList,kTRUE) ;
	//fit our model to the data histogram and store the results (option Save is needed)
	RooFitResult *resultsModel2 = model2.fitTo(rooFitData,Save(1));


	
    // retrieve absolute fractions resulting from the fit
    //
	// Formula for the absolute fraction from the recursive factors of the model
	// v0, v1*(1-v0), v2*(1-v1)*(1-v0) and so on,
	// apart for the last that depends only on the other components: (1-v_{n-1})*...*(1-v0)
	std::vector<RooFormulaVar> absoluteFraction;
	RooArgSet listResults;
	TString formula; 
	for (Int_t ind = 0; ind < nTemp; ind++) {
		if (ind==nTemp-1) formula.Form("1");
		else {
			listResults.add(coefficients[ind]);
			formula.Form("v%d",ind) ;
		}
		for (Int_t k = 0; k < ind ; k++) {
		 	formula.Append(TString::Format("*(1-v%d)",k)); 
		}
//		cout<<formula.Data()<<endl; 
		absoluteFraction.push_back(RooFormulaVar(TString::Format("af%d",ind),TString::Format("af%d",ind),formula.Data(),listResults));
		//insert the results in the proper vector
		res.insert(res.begin()+ind,absoluteFraction[ind].getVal()) ;
 		errRes.insert(errRes.begin()+ind,absoluteFraction[ind].getPropagatedError(*resultsModel2)) ;
	// cout<<	absoluteFraction[ind].getVal()<<endl;
	// 	cout<<	absoluteFraction[ind].getPropagatedError(*resultsModel2)<<endl;
	}
	
    //take the model as an histogram
	TH1D* result = (TH1D*)model2.createHistogram("x");
	//model is a PDF normalized to 1, multiply for the entries of the data histogram
	//NB: do not call here sumw2, because one wants poissonian errors AFTER the rescaling
	result ->Scale(hdata->Integral()); 
	
	RooArgSet *paramList=model2.getParameters(rooFitData);
	paramList->Print("v") ;

	TCanvas *ctemp = new TCanvas();
	ctemp->cd();
	RooPlot* xframe = x.frame() ; 
	rooFitData.plotOn(xframe, Name("data"), MarkerColor(kBlack));
	model2.plotOn(xframe, Name("model2"), MarkerColor(kBlue));

	for (Int_t ind = 0; ind < nTemp; ind++) {
		model2.plotOn(xframe, Components(TString::Format("hpdf%d",ind)), LineColor(color[ind]), LineStyle(kDashed));
	}
	xframe->Draw();
	csq = xframe->chiSquare("model2", "data", nTemp-1);
	ctemp->Delete();
	for (Int_t ind = 0; ind < nTemp; ind++) {
		cout<<res[ind]<<" "<<errRes[ind]<<" "<<csq<<endl;
	}
	return (result);

}

