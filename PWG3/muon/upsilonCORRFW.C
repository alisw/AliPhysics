#include "TH1.h"
#include "TH3.h"
#include "TF1.h"
#include "TF3.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TArrow.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TSpectrum.h"

// RooFit
#include "RooGenericPdf.h"
#include "RooFitResult.h"
#include "RooRealVar.h"


///////////////////////////
// Default Options (Global)
///////////////////////////

//________________
// container check

// display container contents
Bool_t IsShowContOn = kFALSE;
// check bin error
Bool_t IsCheckBinError = kFALSE;

//_______________________
// User setting variables

// enumerate variables
enum { NEVENT, Y, PT, IMASS, TRIG, PTMU, PMU, TRIGSIDE, RABSMU, CHARGE, ETAMU, THABSMU, VZMU, DCAMU }; 

// set efficiency matrix bin 
const Char_t *var0name="y";
const Char_t *var1name="pt";
const Int_t var0bin = 5;
const Int_t var1bin = 7;
Double_t var0lim[var0bin+1] = {-4.0,-3.7,-3.4,-3.1,-2.8,-2.5};
Double_t var1lim[var1bin+1] = {0,2,4,6,8,10,15,20};
const Double_t binnorm = (Double_t)5/2;	// for the last two pt bin

//______________
// Analysis cuts

// Use pt trigger cut
Bool_t IsHpt = kTRUE;

//_______________
// fit parameters 
Bool_t UseRooFit = kFALSE;
if(UseRooFit) { using namespace RooFit; }

// rebinning
Bool_t rebinned = kTRUE;
// 3 gauss for UPSILON family
Char_t *resname = "Upsilon";
Double_t fitrange[2] = {7, 12};		// fit range min, max
Double_t usrrange[2] = {7, 12}; 	// set range user
// setting fit paramter : single exponential(par[2]) + 3 gaussian(par[9])
// param[11] = { A, B, N_Y(1S), mean_Y(1S), sigma_Y(1S), N_Y(2S), mean_Y(2S), sigma_Y(2S), N_Y(3S), mean_Y(3S), sigma_Y(3S)
Double_t param[11] = {3000,-0.3,6.6,9.46,0.230,1.8,10.023,0.230,1,10.355,0.230};
// fix parameter switch
Bool_t IsFixSigma = kTRUE;

////////////////////////////
// UPSILON ANALYSIS MACRO //
////////////////////////////

//______________________________________________________

Bool_t upsilonCORRFW(const char *dataname,								// input file name
										 const Bool_t IsCreateEffMat=kFALSE,	// switch for efficiency matrix creation
										 const Bool_t IsExtSig=kFALSE,					// switch for fitting
										 const Bool_t IsDoCorrection=kFALSE		// switch for correction
										 )
{
// print out startup messages
	printf("******************************************************\n");
	printf("     Starting Correction Framework for Upsilon...\n");
	printf("******************************************************\n");
	printf("\nChecking options...\n");
	printf("Create new efficiency matrix...%s\n",IsCreateEffMat ? "yes":"no");
	printf("Extract Signal...%s\n",IsExtSig ? "yes":"no");
	printf("Do correction...%s\n",IsDoCorrection ? "yes":"no");
	printf("Checking done...\n");

// output files
	char *file_effmat = "";	// efficiency matrix
	char *file_datmat = ""; // data matrix
	char *file_result = ""; // correction result 

// call a fuction to construct efficiency matrix
	if(IsCreateEffMat) file_effmat = createEffMat(dataname);
	else {	// use default efficiency matrix
		//file_effmat = "effCont_mat_residual_aod_50k_pdccut.root";
		file_effmat = "efficiency.container.CFMuonResUpsilon.local.AOD.MC.root";
		printf("Default efficiency matrix: %s\n",file_effmat);
	}

// call a function to perform fit to extract #nignal
	if(IsExtSig) file_datmat = extSignal(dataname);
	else {
		file_datmat = "data.container.CFMuonResUpsilon.PDC09.local.AOD.MC.root";
		printf("Default data container: %s\n",file_datmat);
	}

// call a function to perform a correction on data
	if(IsDoCorrection) file_result = doCorrection(dataname, file_effmat, file_datmat);

// print out closing messages
	printf("******************************************************\n");
	printf("    Finishing Correction Framework for Upsilon...\n");
	printf("******************************************************\n");
	
// display results
	printf("\n here are results!!\n");
	printf("efficiency container is %s\n",file_effmat);
	printf("data container is %s\n",file_datmat);
	printf("correction result file is %s\n",file_result);

	return kTRUE;
}
  
//////////////////////////////
// Create Efficiency matrix //
//////////////////////////////

char* createEffMat(const char *dataname)
{

// print out starting messages
	printf("\nNow starting efficiency matrix creation...\n");

	LoadLib();

// open input container
	printf("Opening container from %s...",dataname);
  TFile *file = new TFile(dataname);	
	AliCFContainer *cont = (AliCFContainer*) (file->Get("container"));
	printf("done\n");

// print out container step and variables
	//cont->Print();
// get N step
	Int_t nstep = cont->GetNStep();
// get N variables
	Int_t nvar = cont->GetNVar();
	printf("Nstep: %d\n", nstep);
	printf("Nvar: %d\n", nvar);
// casting constant for use of array
	const Int_t stepnum = nstep;
	const Int_t varnum = nvar;	
// set variable: bin, min and max
	Double_t varbin[varnum]={ 5,  15, 100, 300, 40, 100, 100, 4, 300, 6,  100, 200, 100, 100 };
	Double_t varmin[varnum]={ 0,  -4,   0,   0,  0,   0,   0, 0,  15, 0, -4.5, 170, -50,   0 };
	Double_t varmax[varnum]={ 5,-2.5, 100, 300,  4, 100, 100, 4,  90, 6, -2.0, 180,  50, 100 };
// set variable epsilon
	Double_t varEpsilon[varnum];
	for (Int_t ie=0; ie<varnum; ie++) varEpsilon[ie] = (varmax[ie]-varmin[ie])/(100*varbin[ie]);
// Display container contents
	if(IsShowContOn) {
		TCanvas *csc[stepnum];

		for(Int_t i=0; i<2; i++) {
			csc[i] = new TCanvas(Form("csc%d",i),Form("container(%d)",i),0,0,1000,600);
			csc[i]->Divide(nvar/2,2);
			for(Int_t j=0; j<nvar; j++) {
				csc[i]->cd(j+1); cont->ShowProjection(j,i)->Draw();
			}
		}
	}
	
//___________
// GridSparse 

// GridSparse from container
	AliCFGridSparse *genSpar = (AliCFGridSparse*)cont->GetGrid(0);		// GEN
	AliCFGridSparse *recSpar = (AliCFGridSparse*)cont->GetGrid(1);		// REC 

// variables for efficiency matrix
	const Int_t nStep=2; 													// number of steps: MC and REC
	const Int_t nVar=2;														// number of variables: var0 and var1
	const Int_t var0dim=var0bin;									// number of bin in var0 dimension of efficiency matrix
	const Int_t var1dim=var1bin; 									// number of bin in var1 dimension of efficiency matrix
	const Int_t binMat[nVar]={var0dim,var1dim};		// summary array

// create new container for efficiency 
	AliCFContainer *effCont = new AliCFContainer("effCont","Upsilon efficiency matrix",nStep,nVar,binMat);
	effCont->SetBinLimits(0,var0lim);
	effCont->SetBinLimits(1,var1lim);
// create new gridsparse for efficiency
// gridsparse of generated events
	AliCFGridSparse *genGrid = new AliCFGridSparse("genGrid","Generated Upsilon matrix",nVar,binMat);
	genGrid->SetBinLimits(0,var0lim);
	genGrid->SetBinLimits(1,var1lim);
// gridsparse of reconstructed events
	AliCFGridSparse *recGrid = new AliCFGridSparse("recGrid","Reconstructed Upsilon matrix",nVar,binMat);
	recGrid->SetBinLimits(0,var0lim);
	recGrid->SetBinLimits(1,var1lim);
	
//____________________
// Loop on matrix bins

// cut on single muon
// theta cut
	varmin[THABSMU] = 171+varEpsilon[THABSMU];
	varmax[THABSMU] = 178-varEpsilon[THABSMU];
// eta cut
	varmin[ETAMU] = -4.0+varEpsilon[ETAMU];
	varmax[ETAMU] = -2.5-varEpsilon[ETAMU];
// rabs cut
	//varmin[RABSMU] = 17.6+varEpsilon[RABSMU];
	//varmin[RABSMU] = 80.0-varEpsilon[RABSMU];

// canvas for fitting
	TCanvas *cfit = new TCanvas("cfit","cfit",0,0,1600,1125);
	cfit->Divide(7,5);

	printf("Loop on matrix bins...");
	for(Int_t bin0=0; bin0<var0dim; bin0++) {			// loop on rapidity (var0dim)
		varmin[Y] = var0lim[bin0]+varEpsilon[Y];
		varmax[Y] = var0lim[bin0+1]-varEpsilon[Y];
		printf("%s range: %.2f < %s < %.2f\n",var0name,varmin[Y],var0name,varmax[Y]);

		for(Int_t bin1=0; bin1<var1dim; bin1++) {		// loop on pt (var1dim)
			varmin[PT] = var1lim[bin1]+varEpsilon[PT];
			varmax[PT] = var1lim[bin1+1]-varEpsilon[PT];
			printf("%s range: %.2f < %s < %.2f\n",var1name,varmin[PT],var1name,varmax[PT]);

			TH1D* gmass = (TH1D*)genSpar->Slice(IMASS,-1,-1,varmin,varmax);
			TH1D* rmass = (TH1D*)recSpar->Slice(IMASS,-1,-1,varmin,varmax);

			cfit->cd((bin0*7)+bin1+1);
			if(rebinned) rmass->Rebin(2);
			// fit reconstructed resonances to get the #_signal
			// declare array to retreive fit results
			Double_t par[4];
			Double_t *parErr;
			// landau (x) gauss : #_signal,mean,sigma,width
			TF1 *func = new TF1("func",gaussXlandau,8,10,4);
			func->SetParameters(30,9.46,0.09,0.004);
			// simple gauss
			//TF1 *func = new TF1("func",normGauss,8,10,3);
			//func->SetParameters(30,9.46,0.09);
			func->SetLineColor(kBlue);
			// fit with likelihood method
			rmass->Fit(func,"0L");
			// retreive fit parameters
			func->GetParameters(&par[0]);
			// fit second times
			func->SetParameters(par);
			rmass->Fit(func,"RL+");
			// retreive fit parameters
			func->GetParameters(&par[0]);
			parErr = func->GetParErrors();
			rmass->Draw("E");
			cfit->Update();

			Double_t genValue = bin1<5 ? gmass->Integral():gmass->Integral()/binnorm;
			Double_t genValueErr = TMath::Sqrt(genValue);
			Double_t recValue = bin1<5 ? par[0]/0.1:(par[0]/0.1)/binnorm;
			Double_t recValueErr = bin1<5 ? parErr[0]/0.1:(parErr[0]/0.1)/binnorm;

			printf("genValue = %.0f +/- %.0f\n",genValue,genValueErr);
			printf("recValue = %.0f +/- %.0f\n",recValue,recValueErr);

			Int_t binCell[nVar]={bin0+1,bin1+1};
			genGrid->SetElement(binCell,genValue);
			genGrid->SetElementError(binCell,genValueErr);
			recGrid->SetElement(binCell,recValue);
			recGrid->SetElementError(binCell,recValueErr);
		}
	}
	printf("done\n");
	printf("saving fit plots as %s...",Form("recfitplots.%s",dataname));
	cfit->SaveAs(Form("recfitplots.%s",dataname));
	printf("done\n");


	printf("fill efficiency container...");
//__________________________
// fill efficiency container
	effCont->SetGrid(0,genGrid);
	effCont->SetGrid(1,recGrid);
	TH1D *gvar0 = effCont->ShowProjection(0,0);
	TH1D *gvar1 = effCont->ShowProjection(1,0);
	TH1D *rvar0 = effCont->ShowProjection(0,1);
	TH1D *rvar1 = effCont->ShowProjection(1,1);
	TCanvas *cvars = new TCanvas("cvars","variables",0,0,800,800);
	cvars->Divide(2,2);
	cvars->cd(1); gvar0->Draw();
	cvars->cd(2); gvar1->Draw();
	cvars->cd(3); rvar0->Draw();
	cvars->cd(4); rvar1->Draw();
	printf("done\n");
// save efficiency container 
	char *outfile = Form("efficiency.%s",dataname);
	printf("saving efficiency container as %s...",outfile);
	effCont->Save(outfile);
	printf("done\n");

//_________________________
// create efficiency matrix
	printf("creating efficiency matrix...");
	AliCFEffGrid *matEff = new AliCFEffGrid("matEff",Form("%s efficiency matrix",resname),*effCont);
	printf("done\n");
// calcualte efficiency
	printf("calculating efficiency...");
	matEff->CalculateEfficiency(1,0);			// REC over MC
	printf("done\n");

//________________________
// display efficiency plot
	TCanvas *ceff = new TCanvas("ceff",Form("%s efficiency",resname),0,0,1200,400);
	ceff->Divide(3,1);
// efficiency over var0
	ceff->cd(1);
	TH1D *effvar0 = matEff->Project(0);
	effvar0->SetName(Form("eff_%s",var0name));
	effvar0->SetMinimum(0); effvar0->SetMaximum(1);
	effvar0->SetTitle(Form("%s efficiency(%s)",resname,var0name));
	effvar0->GetXaxis()->SetTitle(var0name);
	effvar0->GetYaxis()->SetTitle("Efficiency");
	effvar0->Draw("error");
// efficiency over var1
	ceff->cd(2);
	TH1D *effvar1 = matEff->Project(1);
	effvar1->SetName(Form("eff_%s",var1name));
	effvar1->SetMinimum(0); effvar1->SetMaximum(1);
	effvar1->SetTitle(Form("%s efficiency(%s)",resname,var1name));
	effvar1->GetXaxis()->SetTitle(var1name);
	effvar1->GetYaxis()->SetTitle("Efficiency");
	effvar1->Draw("error");
// 2-dimensional efficiency plot (var0, var1)
	ceff->cd(3);
	TH1D *eff2D = matEff->Project(0,1);
	eff2D->SetName("eff_2D");
	eff2D->SetMinimum(0); eff2D->SetMaximum(1);
	eff2D->SetTitle(Form("%s efficiency(%s,%s)",resname,var0name,var1name));
	eff2D->GetXaxis()->SetTitle(var0name);
	eff2D->GetYaxis()->SetTitle(var1name);
	eff2D->GetZaxis()->SetTitle("Efficiency");
	eff2D->Draw("colz");

	printf("\nEfficiency matrix created...\n");
	return outfile;
}

////////////////////
// Extract Signal //
////////////////////

char* extSignal(const char *dataname)
{

// print out starting messages
	printf("\nNow extracting #signal by fitting...\n");

	LoadLib();

// open input container
	printf("Opening container from %s...",dataname);
  TFile *file = new TFile(dataname);	
	AliCFContainer *cont = (AliCFContainer*) (file->Get("container"));
	printf("done\n");

// check dataname : PDC, LHC or Unknown (1, 2, 0)
	Int_t dFlag = 0;
	if(strstr(dataname,"PDC09")) { printf("this is PDC09 data...\n"); dFlag = 1; }
	else if(strstr(dataname,"LHC")) { printf("this is LHC data...\n"); dFlag = 2; }
	else { printf("unknown data...\n"); dFlag = 0; }

// print out container step and variables
// get N step
	Int_t nstep = cont->GetNStep();
// get N variables
	Int_t nvar = cont->GetNVar();
	printf("Nstep: %d\n", nstep);
	printf("Nvar: %d\n", nvar);
// casting constant for use of array
	const Int_t stepnum = nstep;
	const Int_t varnum = nvar;	
// set variable: bin, min and max
	Double_t varbin[varnum]={ 5,  15, 100, 300, 40, 100, 100, 4, 300, 6,  100, 200, 100, 100 };
	Double_t varmin[varnum]={ 0,  -4,   0,   0,  0,   0,   0, 0,  15, 0, -4.5, 170, -50,   0 };
	Double_t varmax[varnum]={ 5,-2.5, 100, 300,  4, 100, 100, 4,  90, 6, -2.0, 180,  50, 100 };
// set variable epsilon
	Double_t varEpsilon[varnum];
	for (Int_t ie=0; ie<varnum; ie++) varEpsilon[ie] = (varmax[ie]-varmin[ie])/(100*varbin[ie]);
// Display container contents
	if(IsShowContOn) {
		TCanvas *csc[stepnum];

		for(Int_t i=0; i<nstep; i++) {
			csc[i] = new TCanvas(Form("csc%d",i),Form("container(%d)",i),0,0,1000,600);
			csc[i]->Divide(nvar/2,2);
			for(Int_t j=0; j<nvar; j++) {
				csc[i]->cd(j+1); cont->ShowProjection(j,i)->Draw();
			}
		}

/*
		TH1D *hCheckMass = cont->ShowProjection(IMASS,1);
		Int_t nbins = hCheckMass->GetXaxis()->GetNbins();
		printf("nbins = %d\n", nbins);
		for(i=1;i<=nbins;i++) printf("bin %d: %f - %f\n",i, hCheckMass->GetBinError(i), TMath::Sqrt(hCheckMass->GetBinContent(i)));
		TCanvas *ctemp = new TCanvas("ctemp","ctemp",0,0,500,500);
		hCheckMass->Draw();
		*/
	}
	
//___________
// GridSparse 

// GridSparse from container
	AliCFGridSparse *genSpar = (AliCFGridSparse*)cont->GetGrid(0);		// GEN for PDC09
	AliCFGridSparse *cintSpar = (AliCFGridSparse*)cont->GetGrid(1);		// REC for PDC09 || CINT && !CMUS
	AliCFGridSparse *cmusSpar = (AliCFGridSparse*)cont->GetGrid(4);		// CMUS only

// variables for data matrix
	const Int_t nStep=2; 													// number of steps: REC(CINT && !CMUS) and CMUS only
	const Int_t nVar=2;														// number of variables: var0 and var1
	const Int_t var0dim=var0bin;									// number of bin in var0 dimension of data matrix
	const Int_t var1dim=var1bin; 									// number of bin in var1 dimension of data matrix
	const Int_t binMat[nVar]={var0dim,var1dim};		// summary array

// create new container for data
	AliCFContainer *dataCont = new AliCFContainer("dataCont","Upsilon data matrix",nStep,nVar,binMat);
	dataCont->SetBinLimits(0,var0lim);
	dataCont->SetBinLimits(1,var1lim);
// create new gridsparse for data
// gridsparse of REC (or CINT && !CMUS) events
	AliCFGridSparse *cintGrid = new AliCFGridSparse("cintGrid","Upsilon data matrix(CINT)",nVar,binMat);
	cintGrid->SetBinLimits(0,var0lim);
	cintGrid->SetBinLimits(1,var1lim);
// gridsparse of CMUS only events
	AliCFGridSparse *cmusGrid = new AliCFGridSparse("cmusGrid","Upsilon matrix(CMUS)",nVar,binMat);
	cmusGrid->SetBinLimits(0,var0lim);
	cmusGrid->SetBinLimits(1,var1lim);
	
//_____________________
// Loops on matrix bins 
	
	if(dFlag == 2) { // LHC data
	// event statistics
		for(Int_t i=0; i<5; i++) {
			varmin[NEVENT] = i+varEpsilon[NEVENT];
			varmax[NEVENT] = i+1-varEpsilon[NEVENT];
			Int_t nEvtInt = ((TH1D*)cintSpar->Slice(IMASS,-1,-1,varmin,varmax))->Integral();
			Int_t nEvtMus = ((TH1D*)cmusSpar->Slice(IMASS,-1,-1,varmin,varmax))->Integral();
			printf("#Event (%d) = %d(CINT1 = %d, CMUS1 = %d)\n",i,nEvtInt+nEvtMus,nEvtInt,nEvtMus);
		}
	
	// common cut on events
	// Rabs cut
		varmin[RABSMU] = 17.6+varEpsilon[RABSMU]; varmax[RABSMU] = 80.0-varEpsilon[RABSMU];
	// acceptance cut on dimuon
		varmin[Y] = -4+varEpsilon[Y];	varmax[Y] = -2.5-varEpsilon[Y];
	
	// canvas for display
		//TCanvas *cmassi = new TCanvas("cmassi","mass distribution(CINT)",0,0,900,600);
		//cmassi->Divide(3,2);
		//TCanvas *cmassm = new TCanvas("cmassm","mass distribution(CMUS)",0,0,900,600);
		//cmassm->Divide(3,2);
	
		TH1D *hmassi[3], *hmassm[3];
	
		// physics selection (0: off, 1: on)
		for(Int_t i=0;i<2;i++) {
			// track-trigger matching
			for(Int_t j=0;j<3;j++) {
				if(IsHpt) {
					switch (j) {
						case 0:
							varmin[TRIG] = 0+varEpsilon[TRIG];
							break;
						case 1:
							varmin[TRIG] = 3+varEpsilon[TRIG];
							break;
						case 2:
							varmin[TRIG] = 3.3+varEpsilon[TRIG];
							break;
					}
				}
				else {
					switch (j) {
						case 0:
							varmin[TRIG] = 0+varEpsilon[TRIG];
							break;
						case 1:
							varmin[TRIG] = 2+varEpsilon[TRIG];
							break;
						case 2:
							varmin[TRIG] = 2.2+varEpsilon[TRIG];
							break;
					}
				}
				// (CINT && !CMUS)
				varmin[NEVENT] = i+1+varEpsilon[NEVENT]; varmax[NEVENT] = i+2-varEpsilon[NEVENT];
				hmassi[j] = (TH1D*)cintSpar->Slice(IMASS,-1,-1,varmin,varmax);
				hmassi[j]->GetXaxis()->SetRangeUser(usrrange[0],usrrange[1]);
				//printf("#Event (CINT1B-%dmatch,hpt) = %d\n",j,hmassi[j]->Integral());
				//cmassi->cd(i*3+j+1);
				//hmassi[j]->Draw();
	
				// (CMUS only)
				varmin[NEVENT] = i+3+varEpsilon[NEVENT]; varmax[NEVENT] = i+4-varEpsilon[NEVENT];
				hmassm[j] = (TH1D*)cmusSpar->Slice(IMASS,-1,-1,varmin,varmax);
				if(rebinned) hmassm[j]->Rebin(4);
				hmassm[j]->GetXaxis()->SetRangeUser(usrrange[0],usrrange[1]);
				//printf("#Event (CMUS1B-%dmatch,hpt) = %d\n",j,hmassm[j]->Integral());
				//cmassm->cd(i*3+j+1);
				//hmassm[j]->Draw();
			}
		}
	
		// draw N matching mass plot and fit
		TCanvas *c1[3];
		TPaveText *info;
	
		Char_t *ptext[6];
	
		for(Int_t i=0; i<3; i++) {
			
			Double_t par[15] = {0, };
			Double_t *parErr = 0;;
			Double_t chi2=0, ndf=0, mean=0, mean_err=0, sigma=0, sigma_err=0, imin=0, imax=0, norm=0, norm_err=0, nsig=0, nsig_err=0,nbkg=0;
	
			//cout << "Bin cont. = " << hmassm[i]->GetBinContent(20) << " +- " << hmassm[i]->GetBinError(20) << endl;
			//cout << "sqrt(Bin cont.) = " << TMath::Sqrt(hmassm[i]->GetBinContent(20)) << endl;
	
			TF1 *bkg = new TF1("bkg",expo,fitrange[0],fitrange[1],2);
 			bkg->SetParameters(param[0],param[1]);
			hmassm[i]->Fit(bkg,"NR");
			bkg->GetParameters(&param[0]);
	
			TF1 *gl = new TF1("gl",global,fitrange[0],fitrange[1],11);
			gl->SetParameters(&param[0]);
	
			// fix yield 2S and 3S normalized from CMS result
			gl->FixParameter(5,1.8); gl->FixParameter(8,1);
			// fix mean from PDG
			gl->FixParameter(6,param[6]); gl->FixParameter(9,param[9]);
			// fix sigma from simulation with current alignment(v2)
			if(IsFixSigma) {gl->FixParameter(4,param[4]); gl->FixParameter(7,param[7]); gl->FixParameter(10,param[10]);}
			else {
				if(i!=2) {gl->SetParLimits(4,param[4]-0.2,param[4]+0.2);gl->FixParameter(7,param[7]); gl->FixParameter(10,param[10]);}
				else {gl->SetParLimits(4,param[4]-0.2,param[4]+0.2);gl->FixParameter(7,param[7]); gl->FixParameter(10,param[10]);}
			}
	
			gl->SetLineColor(kBlue);
			gl->SetLineWidth(2);
			hmassm[i]->Fit(gl,"NR");
			gl->GetParameters(&param[0]);
			parErr = gl->GetParErrors();
	
 	 		bkg->SetParameters(param[0],param[1]);
			bkg->SetParErrors(parErr);
			bkg->SetLineColor(kRed);
			bkg->SetLineStyle(2);
			bkg->SetLineWidth(2);
	
			TF1 *sig = new TF1("sig",normGauss3,fitrange[0],fitrange[1],9);
			sig->SetParameters(param[2],param[3],param[4],param[5],param[6],param[7],param[8],param[9],param[10]);
			sig->SetParErrors(&parErr[2]);
			sig->SetLineColor(kRed);
			sig->SetLineWidth(2);
	
			TF1 *sig1 = new TF1("sig1",normGauss,fitrange[0],fitrange[1],3);
			sig1->SetParameters(param[2],param[3],param[4]);
			sig1->SetParErrors(&parErr[2]);
			sig1->SetLineColor(kMagenta);
			sig1->SetLineWidth(2);
	
			TF1 *sig2 = new TF1("sig2",normGauss,fitrange[0],fitrange[1],3);
			sig2->SetParameters(param[5],param[6],param[7]);
			sig2->SetParErrors(&parErr[5]);
			sig2->SetLineColor(kViolet);
			sig2->SetLineWidth(2);
	
			TF1 *sig3 = new TF1("sig3",normGauss,fitrange[0],fitrange[1],3);
			sig3->SetParameters(param[8],param[9],param[10]);
			sig3->SetParErrors(&parErr[8]);
			sig3->SetLineColor(kPink);
			sig3->SetLineWidth(2);
	
			c1[i] = new TCanvas(Form("c1%d",i),Form("c1%d",i),0,0,800,800);
	
			c1[i]->cd(1); //c1[i]->cd(1)->SetLogy();
			hmassm[i]->SetMinimum(0); 
			hmassm[i]->DrawCopy();
			sig->Draw("same");
			//bkg_2->Draw("same");
			bkg->Draw("same");
			gl->Draw("same");
			sig1->Draw("same");
			sig2->Draw("same");
			sig3->Draw("same");
			
			chi2 = gl->GetChisquare();
			ndf = gl->GetNDF();
			mean = sig->GetParameter(1);
			mean_err = sig->GetParError(1);
			sigma = sig->GetParameter(2);
			sigma_err = sig->GetParError(2);
			imin = mean-2*sigma;
			imax = mean+2*sigma;
			nsig = sig->Integral(imin, imax);
			nsig_err = TMath::Sqrt(nsig);
			//nsig = sig->GetParameter(0);
			//nsig_err = sig->GetParError(0);
			if(rebinned) {
				norm = sig->GetParameter(0)/0.2;
				norm_err = sig->GetParError(0)/0.2;
				nsig = sig1->Integral(imin,imax)/0.2;
				nsig_err = TMath::Sqrt(nsig);
				nbkg = bkg->Integral(imin, imax)/0.2;
			}
			else {
				norm = sig->GetParameter(0)/0.1;
				norm_err = sig->GetParError(0)/0.1;
				nsig = sig1->Integral(imin,imax)/0.1;
				nsig_err = TMath::Sqrt(nsig);
				nbkg = bkg->Integral(imin, imax)/0.1;
			}
			
			TPaveText *info = new TPaveText(0.58,0.48,0.96,0.70,"brNDC");
			info->InsertText(Form("N_{%s} = %.0f #pm %.0f",resname, norm, norm_err));
	 		info->InsertText(Form("Mass = %.3f #pm %.3f GeV/c^{2}", mean, mean_err));
	 		if(IsFixSigma) info->InsertText(Form("#sigma_{%s} = %.0f (fixed) MeV/c^{2}",resname,sigma*1000));
	 		else info->InsertText(Form("#sigma_{%s} = %.0f #pm %.0f MeV/c^{2}",resname,sigma*1000, sigma_err*1000));
			info->InsertText(Form("S/B (2#sigma)= %.3f", nsig/nbkg));
			info->InsertText(Form("Significance (2#sigma) = %.3f" , nsig/TMath::Sqrt(nsig+nbkg)));
	
	 		printf("%d #geq trigger matching(%s)\n",i,IsHpt ? "Hpt" : "Lpt");
	 		printf("#Chi^{2}/ndf = %.3f\n", chi2/ndf);
			printf("N_Upsilon = %.0f +- %.0f\n",norm, norm_err);
			printf("S = %f, B = %f, S/B = %f, sqrt(S+B) = %f, S/sqrt(S+B) = %f\n", nsig, nbkg, nsig/nbkg, TMath::Sqrt(nsig+nbkg), nsig/TMath::Sqrt(nsig+nbkg));
			
			info->SetFillColor(0);
			info->Draw("same");
	
	  	TLegend *leg = new TLegend(0.5928571,0.7442857,0.96,0.9542857,NULL,"brNDC");
 			leg->SetBorderSize(0);
 			leg->SetTextFont(62);
 			leg->SetLineColor(1);
 			leg->SetLineStyle(1);
	  	leg->SetLineWidth(1);
 	  	leg->SetFillColor(0);
 	  	leg->SetFillStyle(1001);
	
   		TLegendEntry *entry=leg->AddEntry(hmassm[i],"Data","lpf");
  		entry->SetFillStyle(1001);
 			entry->SetLineColor(1);
	 		entry->SetLineStyle(1);
 			entry->SetLineWidth(1);
   		entry->SetMarkerColor(1);
   		entry->SetMarkerStyle(8);
   		entry->SetMarkerSize(1);
	
   		entry=leg->AddEntry(sig,"Signal only","lpf");
   		entry->SetFillColor(19);
   		entry->SetLineColor(kRed);
   		entry->SetLineStyle(1);
   		entry->SetLineWidth(2);
   		entry->SetMarkerColor(1);
   		entry->SetMarkerStyle(1);
   		entry->SetMarkerSize(1);
				
   		entry=leg->AddEntry(bkg,"Background only","lpf");
   		entry->SetFillColor(19);
   		entry->SetLineColor(kRed);
   		entry->SetLineStyle(2);
   		entry->SetLineWidth(2);
   		entry->SetMarkerColor(1);
   		entry->SetMarkerStyle(1);
   		entry->SetMarkerSize(1);

   		entry=leg->AddEntry(gl,"Signal + Background","lpf");
   			entry->SetFillColor(19);
   		entry->SetLineColor(kBlue);
   		entry->SetLineStyle(1);
   		entry->SetLineWidth(2);
   		entry->SetMarkerColor(1);
   		entry->SetMarkerStyle(1);
   		entry->SetMarkerSize(1);
   		leg->Draw();
		} 
	}	// dFlag ==2 (LHC data)

	if(dFlag == 1) {	// PDC09 
	// cut on single muon
	// theta cut
		//varmin[THABSMU] = 171+varEpsilon[THABSMU];
		//varmax[THABSMU] = 178-varEpsilon[THABSMU];
	// eta cut
		//varmin[ETAMU] = -4.0+varEpsilon[ETAMU];
		//varmax[ETAMU] = -2.5-varEpsilon[ETAMU];

	// canvas for fitting
		TCanvas *cfit = new TCanvas("cfit","cfit",0,0,1600,1125);
		cfit->Divide(7,5);


		printf("Loop on matrix bins...");
		for(Int_t bin0=0; bin0<var0dim; bin0++) {			// loop on rapidity (var0dim)
			varmin[Y] = var0lim[bin0]+varEpsilon[Y];
			varmax[Y] = var0lim[bin0+1]-varEpsilon[Y];
			printf("%s range: %.2f < %s < %.2f\n",var0name,varmin[Y],var0name,varmax[Y]);
	
			for(Int_t bin1=0; bin1<var1dim; bin1++) {		// loop on pt (var1dim)
				varmin[PT] = var1lim[bin1]+varEpsilon[PT];
				varmax[PT] = var1lim[bin1+1]-varEpsilon[PT];
				printf("%s range: %.2f < %s < %.2f\n",var1name,varmin[PT],var1name,varmax[PT]);

				cfit->cd((bin0*7)+bin1+1);
	
				TH1D* gmass = (TH1D*)genSpar->Slice(IMASS,-1,-1,varmin,varmax);
				TH1D* rmass = (TH1D*)cintSpar->Slice(IMASS,-1,-1,varmin,varmax);
				rmass->Rebin(2);
				gmass->GetXaxis()->SetRangeUser(usrrange[0],usrrange[1]);
				rmass->GetXaxis()->SetRangeUser(usrrange[0],usrrange[1]);

				Double_t genValue, genValueErr, recValue, recValueErr;
				if(!UseRooFit) {
				// fit reconstructed resonances to get the #_signal
				// declare array to retreive fit results
				Double_t param_pdc[12] = {1,1,1,9.46,0.1,0.04,1,10.023,0.1,1,10.38,0.1};
				Double_t *parErr;

				// fit functions
				TF1 *func_sig[3], *func_bkg, *func_gl;

				//___________
				// background
				func_bkg = new TF1("func_bkg",expo,fitrange[0],fitrange[1],2);
				func_bkg->SetParameters(param_pdc[0],param_pdc[1]);
				func_bkg->SetLineColor(kRed);
				func_bkg->SetLineStyle(2);
				func_bkg->SetLineWidth(2);
				rmass->Fit(func_bkg,"NR");
				func_bkg->GetParameters(&param_pdc[0]);

				//_______
				// signal
				// Y(2S)+Y(3S)
				func_sig[1] = new TF1("func_sig2",normGauss2,9.80,10.60,6);
				func_sig[1]->SetParameters(param_pdc[6],param_pdc[7],param_pdc[8],param_pdc[9],param_pdc[10],param_pdc[11]);
				func_sig[1]->SetLineColor(kMagenta);
				func_sig[1]->SetLineWidth(2);
				rmass->Fit(func_sig[1],"NR");
				func_sig[1]->GetParameters(&param_pdc[6]);

/*
				func_sig[0] = new TF1("func_sig1",gaussXlandau,9.20,9.70,4);
				func_sig[0]->SetParameters(param_pdc[2],param_pdc[3],param_pdc[4],param_pdc[5]);
				func_sig[0]->SetLineColor(kPink);
				func_sig[0]->SetLineWidth(2);
				rmass->Fit(func_sig[0],"NR");
				func_sig[0]->GetParameters(&param_pdc[2]);

				func_sig[2] = new TF1("func_sig3",normGauss,10.32,10.57,3);
				func_sig[2]->SetParameters(param_pdc[9],param_pdc[10],param_pdc[11]);
				func_sig[2]->SetLineColor(kViolet);
				func_sig[2]->SetLineWidth(2);
				rmass->Fit(func_sig[2],"0");
				func_sig[2]->SetParameters(param_pdc[9],param_pdc[10],param_pdc[11]);
				rmass->Fit(func_sig[2],"R");
				func_sig[2]->GetParameters(&param_pdc[9]);
				*/

				// Y(1S)+Y(2S)+Y(3S)+Background
				func_gl = new TF1("func_gl",global_pdc,fitrange[0],fitrange[1],12);
				func_gl->SetParameters(param_pdc);
				func_gl->SetParLimits(3,9.30,9.60);
				//func_gl->SetParLimits(4,0.05,0.15);
				func_gl->SetParLimits(7,9.95,10.08);
				//func_gl->SetParLimits(8,0.05,0.15);
				func_gl->SetParLimits(10,10.32,10.57);
				//func_gl->SetParLimits(11,0.05,0.15);
				func_gl->SetLineColor(kBlue);
				func_gl->SetLineWidth(2);
				rmass->Fit(func_gl,"NR");
				func_gl->GetParameters(&param_pdc[0]);
				parErr = func_gl->GetParErrors();

				rmass->Draw("E");
				func_bkg->Draw("same");

				func_gl->Draw("same");
				cfit->Update();
	
				genValue = bin1<5 ? gmass->Integral():gmass->Integral()/binnorm;
				genValueErr = TMath::Sqrt(genValue);
				recValue = bin1<5 ? param_pdc[2]/0.1:(param_pdc[2]/0.1)/binnorm;
				recValueErr = bin1<5 ? parErr[2]/0.1:(parErr[2]/0.1)/binnorm;
				}
				else {		// RooFit
					// Define RooFit variable of the x axis value
					RooRealVar mass("mass","Mass [GeV/c^{2}]",fitrange[0],fitrange[1]);

					// Convert TH1D class to RooFit
					RooDataHist hst("hst","hst",RooArgList(mass),rmass);

					// Define RooFit variable of signal
					RooRealVar Nsig0("Nsig0","Nsig0",0,1000);
					RooRealVar Nsig1("Nsig1","Nsig1",0,1000);
					RooRealVar Nsig2("Nsig2","Nsig2",0,1000);
					RooRealVar Nbg("Nbg","Nbg",-10,5000);

					// Define RooFit variable of gaussian for Y(1S)
					RooRealVar mG0("mG0","mG0",9.30,9.60);
					RooRealVar sG0("sG0","sG0",0.05,0.2);
					//RooGenericPdf genG("genG","genG","Nsig/(sG*TMath::Sqrt(2*TMath::Pi())) * TMath::Exp(-((mass-mG)*(mass-mG))/(2*(sG*sG)))",RooArgSet(mass,Nsig,sG,mG));
					RooGaussian gauss0("gauss0","gauss0",mass,mG0,sG0);

					// Define RooFit variable of gaussian for Y(2S)
					RooRealVar mG1("mG1","mG1",9.80,10.10);
					RooRealVar sG1("sG1","sG1",0.05,0.2);
					RooGaussian gauss1("gauss1","gauss1",mass,mG1,sG1);

					// Define RooFit variable of gaussian for Y(3S)
					RooRealVar mG2("mG2","mG2",10.30,10.60);
					RooRealVar sG2("sG2","sG2",0.05,0.2);
					RooGaussian gauss2("gauss2","gauss2",mass,mG2,sG2);

					// Define RooFit variable of exponential for background
					RooRealVar c("c","c",-20,20);
					RooExponential exp("exp","exp",mass,c);
					
					RooAddPdf glpdf("glpdf","glpdf",RooArgList(gauss0,gauss1,gauss2,exp),RooArgList(Nsig0,Nsig1,Nsig2,Nbg));

					RooFitResult *fitres = glpdf.fitTo(hst,"mhe0");

					// Get a frame from RooRealVar
					RooPlot *frame = mass.frame(Title("#Upsilon mass plot"));
					
					// Get Chi2, mean, sigma and Nsig
					printf("chi2 = %6.4f\n",frame->chiSquare());
					printf("mean  (1S) = %6.2f +/- %4.2f [GeV]\n",mG0.getVal(),mG0.getError());
					printf("sigma (1S) = %6.2f +/- %4.2f [GeV]\n",sG0.getVal(),sG0.getError());
					printf("Nsig  (1S) = %6.0f +/- %4.0f\n",Nsig0.getVal(),Nsig0.getError());
					genValue = bin1<5 ? gmass->Integral():gmass->Integral()/binnorm;
					genValueErr = TMath::Sqrt(genValue);
					recValue = bin1<5 ? Nsig0.getVal():Nsig0.getVal()/binnorm;
					recValueErr = bin1<5 ? Nsig0.getError():Nsig0.getError()/binnorm;

					// Add RooDataHist in the frame
					hst.plotOn(frame);
					glpdf.plotOn(frame,Components(gauss0),LineColor(kRed),LineStyle(kDashed));
					glpdf.plotOn(frame);

					frame->Draw();
				}

				printf("genValue = %.0f +/- %.0f\n",genValue,genValueErr);
				printf("recValue = %.0f +/- %.0f\n",recValue,recValueErr);
	
				// fill grid
				Int_t binCell[nVar]={bin0+1,bin1+1};
				cintGrid->SetElement(binCell,genValue);
				cintGrid->SetElementError(binCell,genValueErr);
				cmusGrid->SetElement(binCell,recValue);
				cmusGrid->SetElementError(binCell,recValueErr);
			}
		}
	}

// set data container
	dataCont->SetGrid(0, cintGrid);
	dataCont->SetGrid(1, cmusGrid);
	
	TCanvas *cdat = new TCanvas("cdat",Form("%s data container",resname),0,0,1200,400);
	cdat->Divide(3,1);
	cdat->cd(1); dataCont->Project(1,0)->Draw();
	cdat->cd(2); dataCont->Project(1,1)->Draw();
	cdat->cd(3); dataCont->Project(1,0,1)->Draw("colz");

// save data container 
	char *outfile = Form("data.%s",dataname);
	printf("saving data container as %s...",outfile);
	dataCont->Save(outfile);
	printf("done\n");

	return outfile;
}	// extSignal 

///////////////////
// Do Correction //
///////////////////

char* doCorrection(char *dataname, char* file_effmat, char* file_datmat)
{
	
// print out starting messages
	printf("\nNow extracting #signal by fitting...\n");

	LoadLib();

	TFile *fEffCont = new TFile(file_effmat);
	TFile *fDataCont = new TFile(file_datmat);

	AliCFContainer *effCont = (AliCFContainer*)fEffCont->Get("effCont");
	AliCFContainer *dataCont = (AliCFContainer*)fDataCont->Get("dataCont");

	AliCFEffGrid *matEff = new AliCFEffGrid("matEff",Form("%s efficiency matrix",resname),*effCont);
	AliCFDataGrid *matGen = new AliCFDataGrid("matGen",Form("%s gen matrix",resname), *dataCont,0);
	AliCFDataGrid *matData = new AliCFDataGrid("matData",Form("%s data matrix",resname),*dataCont,1);

	// calculate efficiency
	matEff->CalculateEfficiency(1,0);

//________________________
// display efficiency plot
	TCanvas *ceff = new TCanvas("ceff",Form("%s efficiency",resname),0,0,1200,400);
	ceff->Divide(3,1);
// efficiency over var0
	ceff->cd(1);
	TH1D *effvar0 = matEff->Project(0);
	effvar0->SetName(Form("eff_%s",var0name));
	effvar0->SetMinimum(0); effvar0->SetMaximum(1);
	effvar0->SetTitle(Form("%s efficiency(%s)",resname,var0name));
	effvar0->GetXaxis()->SetTitle(var0name);
	effvar0->GetYaxis()->SetTitle("Efficiency");
	effvar0->Draw("error");
// efficiency over var1
	ceff->cd(2);
	TH1D *effvar1 = matEff->Project(1);
	effvar1->SetName(Form("eff_%s",var1name));
	effvar1->SetMinimum(0); effvar1->SetMaximum(1);
	effvar1->SetTitle(Form("%s efficiency(%s)",resname,var1name));
	effvar1->GetXaxis()->SetTitle(var1name);
	effvar1->GetYaxis()->SetTitle("Efficiency");
	effvar1->Draw("error");
// 2-dimensional efficiency plot (var0, var1)
	ceff->cd(3);
	TH1D *eff2D = matEff->Project(0,1);
	eff2D->SetName("eff_2D");
	eff2D->SetMinimum(0); eff2D->SetMaximum(1);
	eff2D->SetTitle(Form("%s efficiency(%s,%s)",resname,var0name,var1name));
	eff2D->GetXaxis()->SetTitle(var0name);
	eff2D->GetYaxis()->SetTitle(var1name);
	eff2D->GetZaxis()->SetTitle("Efficiency");
	eff2D->Draw("colz");
//_________________
// display gen plot
	TCanvas *ccor = new TCanvas("ccor",Form("%s correction",resname),0,0,1200,400);
	ccor->Divide(3,1);
// cora over var0
	ccor->cd(1);
	TH1D *genvar0 = matGen->Project(0);
	genvar0->SetName(Form("cor_%s",var0name));
	genvar0->SetMinimum(0);
	genvar0->SetLineColor(kRed);
	genvar0->SetMarkerStyle(22);
	genvar0->SetMarkerColor(kRed);
	genvar0->SetTitle(Form("%s correction(%s)",resname,var0name));
	genvar0->GetXaxis()->SetTitle(var0name);
	genvar0->GetYaxis()->SetTitle("dN/dy");
	genvar0->Draw();
// cora over var1
	ccor->cd(2);
	TH1D *genvar1 = matGen->Project(1);
	genvar1->SetName(Form("cor_%s",var1name));
	genvar1->SetMinimum(0);
	genvar1->SetLineColor(kRed);
	genvar1->SetMarkerStyle(22);
	genvar1->SetMarkerColor(kRed);
	genvar1->SetTitle(Form("%s correction(%s)",resname,var1name));
	genvar1->GetXaxis()->SetTitle(var1name);
	genvar1->GetYaxis()->SetTitle("dN/dpt");
	genvar1->Draw();

//__________________
// display data plot
// data over var0
	ccor->cd(1);
	TH1D *datvar0 = matData->Project(0);
	datvar0->SetName(Form("dat_%s",var0name));
	datvar0->SetMinimum(0);
	datvar0->SetMarkerStyle(22);
	datvar0->SetTitle(Form("%s before correction(%s)",resname,var0name));
	datvar0->GetXaxis()->SetTitle(var0name);
	datvar0->GetYaxis()->SetTitle("dN/dy");
	datvar0->Draw("same");
// data over var1
	ccor->cd(2);
	TH1D *datvar1 = matData->Project(1);
	datvar1->SetName(Form("dat_%s",var1name));
	datvar1->SetMinimum(0);
	datvar1->SetMarkerStyle(22);
	datvar1->SetTitle(Form("%s before correction(%s)",resname,var1name));
	datvar1->GetXaxis()->SetTitle(var1name);
	datvar1->GetYaxis()->SetTitle("dN/dpt");
	datvar1->Draw("same");

//___________________
// perform correction
	matData->ApplyEffCorrection(*matEff);

//________________________
// display correction plot
// correction over var0
	ccor->cd(1);
	TH1D *corvar0 = matData->Project(0);
	corvar0->SetName(Form("cor_%s",var0name));
	corvar0->SetMinimum(0);
	corvar0->SetLineColor(kBlue);
	corvar0->SetMarkerStyle(22);
	corvar0->SetMarkerColor(kBlue);
	corvar0->SetTitle(Form("%s correction(%s)",resname,var0name));
	corvar0->GetXaxis()->SetTitle(var0name);
	corvar0->GetYaxis()->SetTitle("dN/dy");
	corvar0->Draw("same");
// correction over var1
	ccor->cd(2);
	TH1D *corvar1 = matData->Project(1);
	corvar1->SetName(Form("cor_%s",var1name));
	corvar1->SetMinimum(0);
	corvar1->SetLineColor(kBlue);
	corvar1->SetMarkerStyle(22);
	corvar1->SetMarkerColor(kBlue);
	corvar1->SetTitle(Form("%s correction(%s)",resname,var1name));
	corvar1->GetXaxis()->SetTitle(var1name);
	corvar1->GetYaxis()->SetTitle("dN/dpt");
	corvar1->Draw("same");
// 2D correction over var0 and var1
	ccor->cd(3);
	TH2D *cor2D = matData->Project(0,1);
	cor2D->SetName("cor_2D");
	cor2D->SetMinimum(0);
	cor2D->SetTitle(Form("%s correction(%s,%s)",resname,var0name,var1name));
	cor2D->GetXaxis()->SetTitle(var0name);
	cor2D->GetYaxis()->SetTitle(var1name);
	cor2D->GetZaxis()->SetTitle("dN/dpty");
	cor2D->Draw("colz");


	printf("\nCorrection matrix created...\n");
	
// save correction result
	char *outfile = Form("correction-result.%s",dataname);
	printf("saving correction result as %s...",outfile);
	TFile file(outfile,"RECREATE");
	file.cd();
	corvar0->Write();
	corvar1->Write();
	cor2D->Write();
	file.Write();
	file.Close();

	printf("done\n");

	return outfile;
}

//////////////////////////
// DEFINE FIT FUNCTIONS //
//////////////////////////

// single exponential
Double_t expo(Double_t *x, Double_t *par) {
	return par[0] * TMath::Exp((par[1]*x[0]));
}

// double exponential with exclusive region - UPSILON
Double_t dexpo(Double_t *x, Double_t *par) {
	// exclusive region
  if (x[0] > 9.0 && x[0] < 10.0) {
   	TF1::RejectPoint();
   	return 0;
 	}
	return expo(x, &par[0]) + expo(x, &par[2]);
}

// double exponential
Double_t dexpo2(Double_t *x, Double_t *par) {
	return expo(x, &par[0]) + expo(x, &par[2]);
}

// signal : gauss
// fit each resonance
Double_t normGauss(Double_t *x, Double_t *par) {
	return par[0] / (par[2]*TMath::Sqrt(2*TMath::Pi())) * TMath::Exp(-((x[0]-par[1])*(x[0]-par[1]))/(2*(par[2]*par[2])));
}

// fit two resonance
Double_t normGauss2(Double_t *x, Double_t *par) {
	return normGauss(x, &par[0]) + normGauss(x,&par[3]);
}

// fit three resonances
Double_t normGauss3(Double_t *x, Double_t *par) {
	return normGauss(x, &par[0]) + normGauss(x, &par[3]) + normGauss(x, &par[6]);
}

// convolution of gauss and landau
Double_t gaussXlandau(Double_t *var, Double_t *par) {
  // Gaussian convoluted with an inverted Landau function:
  // TMath::Gaus(x,mean,sigma)
  // TMath::Landau(x,mean,width)

  Double_t N = par[0];  // normalization factor = number of entries
  Double_t mG = par[1]; // mean of the gaussian
  Double_t sG = par[2]; // sigma of the gaussian
  Double_t wL = par[3]; // width of the Landau
  Double_t x = var[0];  // variable
  // Numerical constant
  Double_t sqrt2pi = 2.506628275;  // Sqrt{2*Pi}
  Double_t mL0 = -0.22278298;      // maximum Landau location for mean=0
  // Range of integral convolution
  Double_t Nsigma = 5;
  Double_t Sigma = sG + wL;
  //Double_t xMin = mG - Nsigma*Sigma;
  //Double_t xMax = mG + Nsigma*Sigma;
  Double_t xMin = 0;
  Double_t xMax = 2*mG;
  // Integral convolution of a Gaussian with an inverted Landau
  Int_t Nstep = 2000;
  Double_t step = (xMax-xMin) / Nstep;
  Double_t sum = 0;
  Double_t xx;
  for (Double_t i=1.0; i<Nstep; i++) {
    xx = xMin + (i-0.5)*step;
    sum += TMath::Gaus(xx,mG,sG) * TMath::Landau(-(x-xx),0,wL);
  }
  // Normalized result
  return N*sum*step/(sqrt2pi*sG*wL);
}

// global fit : UPSILON Family w/ sigle exponential
Double_t global(Double_t *x, Double_t *par) {
	return expo(x,&par[0]) + normGauss3(x, &par[2]);
}

// global fit for PDC09 : UPSILON Family w/ sigle exponential
Double_t global_pdc(Double_t *x, Double_t *par) {
	return expo(x,&par[0]) + gaussXlandau(x, &par[2]) + normGauss2(x, &par[6]);
}

// Crystal ball fit
Double_t CrystalBall(Double_t *x, Double_t *par) {
	//Crystal ball function for signal
	//parameters are 0:alpha,1:n,2:mean,3:sigma,4:number of expected events
	Double_t m=x[0];
	Double_t alpha=par[0];
	Double_t n=par[1];
	Double_t m0=par[2];
	Double_t sigma=par[3];
	Double_t N=par[4];
	
	Double_t t = (m-m0)/sigma;
	if(alpha<0) t = -t;

	Double_t absAlpha = TMath::Abs((Double_t)alpha);

	if(t >= -absAlpha) {
		return N/(sigma*TMath::Sqrt(2*TMath::Pi()))*exp(-0.5*t*t);
		//return N*exp(-0.5*t*t);
	}
	else {
		Double_t intn=0;
		Double_t fracn=modf(n, &intn);
		if(n>10) return 0;
		if((n/absAlpha<0) && TMath::Abs(fracn)>0) return 0; // TMath::Power(x,y): negative x and non-integer y not supported!!
		Double_t a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
		Double_t b = n/absAlpha - absAlpha;

		return N/(sigma*TMath::Sqrt(2*TMath::Pi()))*(a/TMath::Power(b - t, n));
		//return N*(a/TMath::Power(b - t, n));
	}
}

// Loading libraries
void LoadLib() {
	printf("Loading libraries...");
// loading libraries
	gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include  -I$ROOTSYS/include"); 
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD") ;
	gSystem->Load("libANALYSIS.so");
	gSystem->Load("libANALYSISalice.so");
	gSystem->Load("libCORRFW.so");

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(0);
	printf("done\n");
}
