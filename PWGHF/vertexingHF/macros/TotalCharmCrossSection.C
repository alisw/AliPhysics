//Author: Jeremy Wilkinson, jwilkinson@physi.uni-heidelberg.de
//Routine to extrapolate D meson cross sections using FONLL calculations
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>
#include <iostream>
#include <TH1.h>
#include <TH1D.h>
#include <TKey.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TPave.h>
#include <TPaveText.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <fstream>

void TotalCharmCrossSection() {

	// FONLL predictions from _pdf.rooT files
	TFile *D0_05 = TFile::Open("D07TeV_FONLL_y0.5_pdf.root");
	TFile *Dst_05 = TFile::Open("Dstar7TeV_FONLL_y0.5_pdf.root");
	TFile *Dplus_05 = TFile::Open("Dplus7TeV_FONLL_y0.5_pdf.root");

	TFile *D0_ally = TFile::Open("D07TeV_FONLL_y12_pdf.root");
	TFile *Dst_ally = TFile::Open("Dstar7TeV_FONLL_y12_pdf.root");
	TFile *Dplus_ally = TFile::Open("Dplus7TeV_FONLL_y12_pdf.root");
	
	//Spectra for error propagation
	TFile *specD0 = TFile::Open("HFPtSpectrum_D0Kpi_combinedFD_rebinnedth_150311_newsigma_7TeV.root");
	TFile *specDst = TFile::Open("HFPtSpectrum_DstarD0pi_combineFD_rebinnedth_newnorm_120311.root");
	TFile *specDplus = TFile::Open("HFPtSpectrum_DplusKpipi_combinedFD_150311_newsigma_7TeV.root");

	TH1F *hD0stat = (TH1F*)specD0->Get("histoSigmaCorr");
	TH1F *hDststat = (TH1F*)specDst->Get("histoSigmaCorr");
	TH1F *hDplusstat = (TH1F*)specDplus->Get("histoSigmaCorr");


	//Initialising the histo arrays
	TH1F **histosD0_ally, **histosD0_05, **histosDst_ally, **histosDst_05, **histosDplus_ally, **histosDplus_05;
	
	histosD0_ally = new TH1F*[11];
	histosD0_05 = new TH1F*[11];
	histosDst_ally = new TH1F*[11];
	histosDst_05 = new TH1F*[11];
	histosDplus_ally = new TH1F*[11];
	histosDplus_05 = new TH1F*[11];


//Getting all the necessary histos for extrapolation ratios

	histosD0_ally[0] = (TH1F*)D0_ally->Get("FONLL_central_prediction");
	histosD0_ally[1] = (TH1F*)D0_ally->Get("FONLL_fr55_prediction");
	histosD0_ally[2] = (TH1F*)D0_ally->Get("FONLL_fr22_prediction");
	histosD0_ally[3] = (TH1F*)D0_ally->Get("FONLL_fr21_prediction");
	histosD0_ally[4] = (TH1F*)D0_ally->Get("FONLL_fr12_prediction");
	histosD0_ally[5] = (TH1F*)D0_ally->Get("FONLL_fr15_prediction");
	histosD0_ally[6] = (TH1F*)D0_ally->Get("FONLL_fr51_prediction");
	histosD0_ally[7] = (TH1F*)D0_ally->Get("FONLL_minmass_prediction");
	histosD0_ally[8] = (TH1F*)D0_ally->Get("FONLL_maxmass_prediction");
	histosD0_ally[9] = (TH1F*)D0_ally->Get("FONLL_minpdf_prediction");
	histosD0_ally[10] = (TH1F*)D0_ally->Get("FONLL_maxpdf_prediction");
	
	histosD0_05[0] = (TH1F*)D0_05->Get("FONLL_central_prediction");
	histosD0_05[1] = (TH1F*)D0_05->Get("FONLL_fr55_prediction");
	histosD0_05[2] = (TH1F*)D0_05->Get("FONLL_fr22_prediction");
	histosD0_05[3] = (TH1F*)D0_05->Get("FONLL_fr21_prediction");
	histosD0_05[4] = (TH1F*)D0_05->Get("FONLL_fr12_prediction");
	histosD0_05[5] = (TH1F*)D0_05->Get("FONLL_fr15_prediction");
	histosD0_05[6] = (TH1F*)D0_05->Get("FONLL_fr51_prediction");
	histosD0_05[7] = (TH1F*)D0_05->Get("FONLL_minmass_prediction");
	histosD0_05[8] = (TH1F*)D0_05->Get("FONLL_maxmass_prediction");
	histosD0_05[9] = (TH1F*)D0_05->Get("FONLL_minpdf_prediction");
	histosD0_05[10] = (TH1F*)D0_05->Get("FONLL_maxpdf_prediction");
	
	histosDst_ally[0] = (TH1F*)Dst_ally->Get("FONLL_central_prediction");
	histosDst_ally[1] = (TH1F*)Dst_ally->Get("FONLL_fr55_prediction");
	histosDst_ally[2] = (TH1F*)Dst_ally->Get("FONLL_fr22_prediction");
	histosDst_ally[3] = (TH1F*)Dst_ally->Get("FONLL_fr21_prediction");
	histosDst_ally[4] = (TH1F*)Dst_ally->Get("FONLL_fr12_prediction");
	histosDst_ally[5] = (TH1F*)Dst_ally->Get("FONLL_fr15_prediction");
	histosDst_ally[6] = (TH1F*)Dst_ally->Get("FONLL_fr51_prediction");
	histosDst_ally[7] = (TH1F*)Dst_ally->Get("FONLL_minmass_prediction");
	histosDst_ally[8] = (TH1F*)Dst_ally->Get("FONLL_maxmass_prediction");
	histosDst_ally[9] = (TH1F*)Dst_ally->Get("FONLL_minpdf_prediction");
	histosDst_ally[10] = (TH1F*)Dst_ally->Get("FONLL_maxpdf_prediction");

	
	histosDst_05[0] = (TH1F*)Dst_05->Get("FONLL_central_prediction");
	histosDst_05[1] = (TH1F*)Dst_05->Get("FONLL_fr55_prediction");
	histosDst_05[2] = (TH1F*)Dst_05->Get("FONLL_fr22_prediction");
	histosDst_05[3] = (TH1F*)Dst_05->Get("FONLL_fr21_prediction");
	histosDst_05[4] = (TH1F*)Dst_05->Get("FONLL_fr12_prediction");
	histosDst_05[5] = (TH1F*)Dst_05->Get("FONLL_fr15_prediction");
	histosDst_05[6] = (TH1F*)Dst_05->Get("FONLL_fr51_prediction");
	histosDst_05[7] = (TH1F*)Dst_05->Get("FONLL_minmass_prediction");
	histosDst_05[8] = (TH1F*)Dst_05->Get("FONLL_maxmass_prediction");
	histosDst_05[9] = (TH1F*)Dst_05->Get("FONLL_minpdf_prediction");
	histosDst_05[10] = (TH1F*)Dst_05->Get("FONLL_maxpdf_prediction");
	
	histosDplus_ally[0] = (TH1F*)Dplus_ally->Get("FONLL_central_prediction");
	histosDplus_ally[1] = (TH1F*)Dplus_ally->Get("FONLL_fr55_prediction");
	histosDplus_ally[2] = (TH1F*)Dplus_ally->Get("FONLL_fr22_prediction");
	histosDplus_ally[3] = (TH1F*)Dplus_ally->Get("FONLL_fr21_prediction");
	histosDplus_ally[4] = (TH1F*)Dplus_ally->Get("FONLL_fr12_prediction");
	histosDplus_ally[5] = (TH1F*)Dplus_ally->Get("FONLL_fr15_prediction");
	histosDplus_ally[6] = (TH1F*)Dplus_ally->Get("FONLL_fr51_prediction");
	histosDplus_ally[7] = (TH1F*)Dplus_ally->Get("FONLL_minmass_prediction");
	histosDplus_ally[8] = (TH1F*)Dplus_ally->Get("FONLL_maxmass_prediction");
	histosDplus_ally[9] = (TH1F*)Dplus_ally->Get("FONLL_minpdf_prediction");
	histosDplus_ally[10] = (TH1F*)Dplus_ally->Get("FONLL_maxpdf_prediction");
	
	histosDplus_05[0] = (TH1F*)Dplus_05->Get("FONLL_central_prediction");
	histosDplus_05[1] = (TH1F*)Dplus_05->Get("FONLL_fr55_prediction");
	histosDplus_05[2] = (TH1F*)Dplus_05->Get("FONLL_fr22_prediction");
	histosDplus_05[3] = (TH1F*)Dplus_05->Get("FONLL_fr21_prediction");
	histosDplus_05[4] = (TH1F*)Dplus_05->Get("FONLL_fr12_prediction");
	histosDplus_05[5] = (TH1F*)Dplus_05->Get("FONLL_fr15_prediction");
	histosDplus_05[6] = (TH1F*)Dplus_05->Get("FONLL_fr51_prediction");
	histosDplus_05[7] = (TH1F*)Dplus_05->Get("FONLL_minmass_prediction");
	histosDplus_05[8] = (TH1F*)Dplus_05->Get("FONLL_maxmass_prediction");
	histosDplus_05[9] = (TH1F*)Dplus_05->Get("FONLL_minpdf_prediction");
	histosDplus_05[10] = (TH1F*)Dplus_05->Get("FONLL_maxpdf_prediction");

//Output file for extrap. ratios

	ofstream outD0("outputD0.txt"), outDst("outputDst.txt"), outDplus("outputDplus.txt");

//Initialise doubles for max, min, central for each meson
Double_t ratscalmaxD0 = 0, ratscalminD0 = 9999, ratscalmaxDst = 0, ratscalminDst = 9999, ratscalmaxDplus = 0, ratscalminDplus = 9999,ratmassmaxD0 = 0, ratmassminD0 = 9999, ratmassmaxDst = 0, ratmassminDst = 9999, ratmassmaxDplus = 0, ratmassminDplus = 9999, ratpdfmaxD0 = 0, ratpdfminD0 = 9999, ratpdfmaxDst = 0, ratpdfminDst = 9999, ratpdfmaxDplus = 0, ratpdfminDplus = 9999, ratD0cent, ratDstcent, ratDpluscent,
ratD0[11], ratDst[11], ratDplus[11];


	//Central value is with central parameters
	ratD0cent = histosD0_ally[0]->Integral("width")/histosD0_05[0]->Integral(21,120,"width");
	ratDstcent = histosDst_ally[0]->Integral("width")/histosDst_05[0]->Integral(21,120,"width");
	ratDpluscent = histosDplus_ally[0]->Integral("width")/histosDplus_05[0]->Integral(21,120,"width");
//Headers
outD0 << "Name\tpT\ty\toverall\n";
outDst << "Name\tpT\ty\toverall\n";
outDplus << "Name\tpT\ty\toverall\n";


//Determining maximum and minimum ratios, and writing all ratios to file for later analysis 
	for (Int_t i = 0; i<=10; i++) {

		ratD0[i] = histosD0_ally[i]->Integral("width")/histosD0_05[i]->Integral(21,120,"width");

							//Scale variation for histos between i=1 and i=6
			if (ratD0[i] > ratscalmaxD0 && 0 <= i <= 6) {ratscalmaxD0 = ratD0[i];}
			if (ratD0[i] < ratscalminD0 && 0 <= i <= 6) {ratscalminD0 = ratD0[i];}
			
							//Mass variation for 0, 7, 8
			if (ratD0[i] > ratmassmaxD0 && (i == 0 || i == 7 || i == 8)) {ratmassmaxD0 = ratD0[i];}
			if (ratD0[i] < ratmassminD0 && (i == 0 || i == 7 || i == 8)) {ratmassminD0 = ratD0[i];}

							//PDF uncertainty is 0, 9 and 10
			if (ratD0[i] > ratpdfmaxD0 && (i == 0 || i == 9 || i == 10)) {ratpdfmaxD0 = ratD0[i];}
			if (ratD0[i] < ratpdfminD0 && (i == 0 || i == 9 || i == 10)) {ratpdfminD0 = ratD0[i];}

		outD0 << histosD0_05[i]->GetName() << '\t'
		<< histosD0_05[i]->Integral("width")/histosD0_05[i]->Integral(21,120,"width") << '\t'
		<< histosD0_ally[i]->Integral("width")/histosD0_05[i]->Integral("width") << '\t'
		<< ratD0[i] << endl;
		
		ratDst[i] = histosDst_ally[i]->Integral("width")/histosDst_05[i]->Integral(21,120,"width");
			if (ratDst[i] > ratscalmaxDst && 0 <= i <= 6) {ratscalmaxDst = ratDst[i];}
			if (ratDst[i] < ratscalminDst && 0 <= i <= 6) {ratscalminDst = ratDst[i];}
			
			if (ratDst[i] > ratmassmaxDst && (i == 0 || i == 7 || i == 8)) {ratmassmaxDst = ratDst[i];}
			if (ratDst[i] < ratmassminDst && (i == 0 || i == 7 || i == 8)) {ratmassminDst = ratDst[i];}
			
			if (ratDst[i] > ratpdfmaxDst && (i == 0 || i == 9 || i == 10)) {ratpdfmaxDst = ratDst[i];}
			if (ratDst[i] < ratpdfminDst && (i == 0 || i == 9 || i == 10)) {ratpdfminDst = ratDst[i];}

		outDst << histosDst_05[i]->GetName() << '\t'
		<< histosDst_05[i]->Integral("width")/histosDst_05[i]->Integral(21,120,"width") << '\t'
		<< histosDst_ally[i]->Integral("width")/histosDst_05[i]->Integral("width") << '\t'
		<< ratDst[i] << endl;


		ratDplus[i] = histosDplus_ally[i]->Integral("width")/histosDplus_05[i]->Integral(21,120,"width");
			if (ratDplus[i] > ratscalmaxDplus && 0 <= i <= 6) {ratscalmaxDplus = ratDplus[i];}
			if (ratDplus[i] < ratscalminDplus && 0 <= i <= 6) {ratscalminDplus = ratDplus[i];}
			
			if (ratDplus[i] > ratmassmaxDplus && (i == 0 || i == 7 || i == 8)) {ratmassmaxDplus = ratDplus[i];}
			if (ratDplus[i] < ratmassminDplus && (i == 0 || i == 7 || i == 8)) {ratmassminDplus = ratDplus[i];}

			if (ratDplus[i] > ratpdfmaxDplus && (i == 0 || i == 9 || i == 10)) {ratpdfmaxDplus = ratDplus[i];}
			if (ratDplus[i] < ratpdfminDplus && (i == 0 || i == 9 || i == 10)) {ratpdfminDplus = ratDplus[i];}

		outDplus << histosDplus_05[i]->GetName() << '\t'
		<< histosDplus_05[i]->Integral("width")/histosDplus_05[i]->Integral(21,120,"width") << '\t'
		<< histosDplus_ally[i]->Integral("width")/histosDplus_05[i]->Integral("width") << '\t'
		<< ratDplus[i] << endl;

	
	
	}
	//Maximum and minimum error on centrals:
	Double_t ratscalD0UpErr, ratscalD0LowErr, ratscalDstUpErr, ratscalDstLowErr, ratscalDplusUpErr, ratscalDplusLowErr, ratmassD0UpErr, ratmassD0LowErr, ratmassDstUpErr, ratmassDstLowErr, ratmassDplusUpErr, ratmassDplusLowErr, ratpdfD0UpErr, ratpdfD0LowErr, ratpdfDstUpErr, ratpdfDstLowErr, ratpdfDplusUpErr, ratpdfDplusLowErr, ratD0UpErr, ratD0LowErr, ratDstUpErr, ratDstLowErr, ratDplusUpErr, ratDplusLowErr;

	ratscalD0UpErr = ratscalmaxD0 - ratD0cent;
	ratscalD0LowErr = ratD0cent - ratscalminD0;
	ratscalDstUpErr = ratscalmaxDst - ratDstcent;
	ratscalDstLowErr = ratDstcent - ratscalminDst;
	ratscalDplusUpErr = ratscalmaxDplus - ratDpluscent;
	ratscalDplusLowErr = ratDpluscent - ratscalminDplus;

	ratmassD0UpErr = ratmassmaxD0 - ratD0cent;
	ratmassD0LowErr = ratD0cent - ratmassminD0;
	ratmassDstUpErr = ratmassmaxDst - ratDstcent;
	ratmassDstLowErr = ratDstcent - ratmassminDst;
	ratmassDplusUpErr = ratmassmaxDplus - ratDpluscent;
	ratmassDplusLowErr = ratDpluscent - ratmassminDplus;

	ratpdfD0UpErr = ratpdfmaxD0 - ratD0cent;
	ratpdfD0LowErr = ratD0cent - ratpdfminD0;
	ratpdfDstUpErr = ratpdfmaxDst - ratDstcent;
	ratpdfDstLowErr = ratDstcent - ratpdfminDst;
	ratpdfDplusUpErr = ratpdfmaxDplus - ratDpluscent;
	ratpdfDplusLowErr = ratDpluscent - ratpdfminDplus;

	//Overall uncertainties are mass uncertainty and scale uncertainty added in quadrature
	ratD0UpErr = sqrt(ratscalD0UpErr*ratscalD0UpErr+ratmassD0UpErr*ratmassD0UpErr+ratpdfD0UpErr*ratpdfD0UpErr);
	ratD0LowErr = sqrt(ratscalD0LowErr*ratscalD0LowErr+ratmassD0LowErr*ratmassD0LowErr+ratpdfD0LowErr*ratpdfD0LowErr);
	ratDstUpErr = sqrt(ratscalDstUpErr*ratscalDstUpErr+ratmassDstUpErr*ratmassDstUpErr+ratpdfDstUpErr*ratpdfDstUpErr);
	ratDstLowErr = sqrt(ratscalDstLowErr*ratscalDstLowErr+ratmassDstLowErr*ratmassDstLowErr+ratpdfDstLowErr*ratpdfDstLowErr);
	ratDplusUpErr = sqrt(ratscalDplusUpErr*ratscalDplusUpErr+ratmassDplusUpErr*ratmassDplusUpErr+ratpdfDplusUpErr*ratpdfDplusUpErr);
	ratDplusLowErr = sqrt(ratscalDplusLowErr*ratscalDplusLowErr+ratmassDplusLowErr*ratmassDplusLowErr+ratpdfDplusLowErr*ratpdfDplusLowErr);
	

	//Output ratios to file for future reference:

	ofstream outratios("ratios.txt");
	outratios << "D0\ncentral\t" << ratD0cent << "\nscales:\t+" << ratscalD0UpErr << "\t-" << ratscalD0LowErr << "\nmass:\t+" << ratmassD0UpErr << "\t-" << ratmassD0LowErr << "\nPDF:\t+" << ratpdfD0UpErr << "\t-" << ratpdfD0LowErr << "\nOverall:\t+" << ratD0UpErr << "\t-" << ratD0LowErr
<<	"\n\nD*\ncentral\t" << ratDstcent << "\nscales:\t+" << ratscalDstUpErr << "\t-" << ratscalDstLowErr << "\nmass:\t+" << ratmassDstUpErr << "\t-" << ratmassDstLowErr << "\nPDF:\t+" << ratpdfDstUpErr << "\t-" << ratpdfDstLowErr << "\nOverall:\t+" << ratDstUpErr << "\t-" << ratDstLowErr
<<	"\n\nD+\ncentral\t" << ratDpluscent << "\nscales:\t+" << ratscalDplusUpErr << "\t-" << ratscalDplusLowErr << "\nmass:\t+" << ratmassDplusUpErr << "\t-" << ratmassDplusLowErr << "\nPDF:\t+" << ratpdfDplusUpErr << "\t-" << ratpdfDplusLowErr << "\nOverall:\t+" << ratDplusUpErr << "\t-" << ratDplusLowErr << '\n';
	//Visible spectra:
	Double_t D0stati[6], Dststati[6], Dplusstati[6];

	Double_t D0stat =0, Dststat=0, Dplusstat=0, D0statsq=0, Dststatsq=0, Dplusstatsq=0, D0systhi=0, D0systlo=0, Dstsysthi=0, Dstsystlo=0, Dplussysthi=0, Dplussystlo=0;
	for (Int_t i = 4; i<=9; i++) {
		D0stati[i-4] = hD0stat->GetBinError(i)/1e6/3.89e-2;
		Dststati[i-4] = hDststat->GetBinError(i-1)/1e6/3.89e-2/67.7e-2;
		Dplusstati[i-4] = hDplusstat->GetBinError(i-3)/1e6/0.0922;
		
		//Statistical added in quadrature (bin widths included)
		D0statsq += pow(D0stati[i-4]*hD0stat->GetBinWidth(i),2);
		Dststatsq += pow(Dststati[i-4]*hDststat->GetBinWidth(i-1),2);
		Dplusstatsq += pow(Dplusstati[i-4]*hDplusstat->GetBinWidth(i-3),2);
			
	}


	//Rooting to get actual values from delta-sigma-squared
	D0stat = sqrt(D0statsq); Dststat = sqrt(Dststatsq); Dplusstat = sqrt(Dplusstatsq);

	//Corrected systematics are not in a ROOT file, input here manually
	D0systhi = 967844/1e6/3.89e-2;
	D0systlo = 1849380/1e6/3.89e-2;
	Dstsysthi = 232479/1e6/3.89e-2/67.7e-2;
	Dstsystlo = 379356/1e6/3.89e-2/67.7e-2;
	Dplussysthi = 1710680/1e6/0.0922;
	Dplussystlo = 3402420/1e6/0.0922;

	
	//Visible D cross sections from data file, "width" switch incorporates bin widths
	Double_t D0vis = hD0stat->Integral("width")/1e6/3.89e-2,
		Dstvis = hDststat->Integral("width")/1e6/67.7e-2/3.89e-2,
		Dplusvis = hDplusstat->Integral("width")/1e6/0.0922;
	
	//Branching ratio errors as percentage:
	Double_t ErrBrD0Perc = 0.0002785/0.0216673, ErrBrDstPerc = 9.29153e-5/0.0062678, ErrBrDplusPerc = 0.0004746/0.0208372;
	
	//Luminosity errors as percentages

	Double_t ErrLumD0Perc=0.07, ErrLumDstPerc=0.07, ErrLumDplusPerc=0.07;

	//Initialising doubles for final results:
	Double_t crossD0, crossD0stat, crossD0systhi, crossD0systlo, crossD0lum, crossD0br, crossD0extrhi, crossD0extrlo, crossDst, crossDststat, crossDstsysthi, crossDstsystlo, crossDstlum, crossDstbr, crossDstextrhi, crossDstextrlo, crossDplus, crossDplusstat, crossDplussysthi, crossDplussystlo, crossDpluslum, crossDplusbr, crossDplusextrhi, crossDplusextrlo;

	//Also initialising for final ccbar cross sec:
	Double_t ccbarD0, ccbarD0stat, ccbarD0systhi, ccbarD0systlo, ccbarD0lum, ccbarD0extrhi, ccbarD0extrlo, ccbarDst, ccbarDststat, ccbarDstsysthi, ccbarDstsystlo, ccbarDstlum, ccbarDstextrhi, ccbarDstextrlo, ccbarDplus, ccbarDplusstat, ccbarDplussysthi, ccbarDplussystlo, ccbarDpluslum, ccbarDplusextrhi, ccbarDplusextrlo;

	//Cross section is visible multiplied by central ratio
	crossD0 = D0vis * ratD0cent;
	crossDst = Dstvis * ratDstcent;
	crossDplus = Dplusvis * ratDpluscent;
	
	//Stat on each is statitstical of visible multiplied by central ratio
	crossD0stat = D0stat * ratD0cent;
	crossDststat = Dststat * ratDstcent;
	crossDplusstat = Dplusstat * ratDpluscent;
	
	//Systematic the same, but upper and lower separately
	crossD0systhi = D0systhi * ratD0cent;
	crossD0systlo = D0systlo * ratD0cent;
	crossDstsysthi = Dstsysthi * ratDstcent;
	crossDstsystlo = Dstsystlo * ratDstcent;
	crossDplussysthi = Dplussysthi * ratDpluscent;
	crossDplussystlo = Dplussystlo * ratDpluscent;
	
	//Luminosity errors are percentage multiplied by central result
	crossD0lum = ErrLumD0Perc * crossD0;
	crossDstlum = ErrLumDstPerc * crossDst;
	crossDpluslum = ErrLumDplusPerc * crossDplus;
	
	//Similarly for branching ratios
	crossD0br = ErrBrD0Perc * crossD0;
	crossDstbr = ErrBrDstPerc * crossDst;
	crossDplusbr = ErrBrDplusPerc * crossDplus;

	//Extrapolation hi/lo errors are overall hi/lo ratio errors multiplied by visible
	crossD0extrhi = ratD0UpErr * D0vis;
	crossD0extrlo = ratD0LowErr * D0vis;
	crossDstextrhi = ratDstUpErr * Dstvis;
	crossDstextrlo = ratDstLowErr * Dstvis;
	crossDplusextrhi = ratDplusUpErr * Dplusvis;
	crossDplusextrlo = ratDplusLowErr * Dplusvis;

	//ccbar cross section: divide result by fragmentation ratio. This is calculated the same way for all uncertainties
	//Fragmentations:
	Double_t D0frag = 0.557, Dstfrag = 0.238, Dplusfrag = 0.226;

	//D0
	ccbarD0 = crossD0/D0frag;
	ccbarD0stat = crossD0stat/D0frag;
	ccbarD0systhi = crossD0systhi/D0frag;
	ccbarD0systlo = crossD0systlo/D0frag;
	ccbarD0lum = crossD0lum/D0frag;
	ccbarD0extrhi = crossD0extrhi/D0frag;
	ccbarD0extrlo = crossD0extrlo/D0frag;

	//D*
	ccbarDst = crossDst/Dstfrag;
	ccbarDststat = crossDststat/Dstfrag;
	ccbarDstsysthi = crossDstsysthi/Dstfrag;
	ccbarDstsystlo = crossDstsystlo/Dstfrag;
	ccbarDstlum = crossDstlum/Dstfrag;
	ccbarDstextrhi = crossDstextrhi/Dstfrag;
	ccbarDstextrlo = crossDstextrlo/Dstfrag;

	//D+
	ccbarDplus = crossDplus/Dplusfrag;
	ccbarDplusstat = crossDplusstat/Dplusfrag;
	ccbarDplussysthi = crossDplussysthi/Dplusfrag;
	ccbarDplussystlo = crossDplussystlo/Dplusfrag;
	ccbarDpluslum = crossDpluslum/Dplusfrag;
	ccbarDplusextrhi = crossDplusextrhi/Dplusfrag;
	ccbarDplusextrlo = crossDplusextrlo/Dplusfrag;

	
	//To combine: Add D0 + D+, divide by total frag. ratio of the two. Take weighted average of this value with D* using 1/(stat)^2 as weight.


	//Initialising variables:
	Double_t totalfrag, D0Dplusccbar, D0Dplusccbarstat, D0Dplusccbarsysthi, D0Dplusccbarsystlo, D0Dplusccbarlum, D0Dplusccbarextrhi, D0Dplusccbarextrlo;

	totalfrag = D0frag + Dplusfrag;

	//Propagating this method through:
	D0Dplusccbar = (crossD0 + crossDplus)/totalfrag;
	//Errors are quadratic sum of D0 and D+ errors, divided by total frag. ratio
	D0Dplusccbarstat = sqrt(crossD0stat * crossD0stat + crossDplusstat * crossDplusstat)/totalfrag;
	D0Dplusccbarsysthi = sqrt(crossD0systhi * crossD0systhi + crossDplussysthi * crossDplussysthi)/totalfrag;
	D0Dplusccbarsystlo = sqrt(crossD0systlo * crossD0systlo + crossDplussystlo * crossDplussystlo)/totalfrag;
	D0Dplusccbarlum = sqrt(crossD0lum * crossD0lum + crossDpluslum * crossDpluslum)/totalfrag;
	D0Dplusccbarextrhi = sqrt(crossD0extrhi * crossD0extrhi + crossDplusextrhi * crossDplusextrhi)/totalfrag;
	D0Dplusccbarextrlo = sqrt(crossD0extrlo * crossD0extrlo + crossDplusextrlo * crossDplusextrlo)/totalfrag;
	
cout << "Just to test, this is D0Dplusccbar:   " << D0Dplusccbar << "\n\n\n";

	
	//Weighted average:
	Double_t weightD0Dplus, weightDst, totalweight, overallccbar, overallccbarstat, overallccbarsysthi, overallccbarsystlo, 	overallccbarlum, overallccbarextrhi, overallccbarextrlo;
	

	//Weights are 1/sigma(stat)^2 (absolute value of error, not relative errors)
	weightD0Dplus = 1/(D0Dplusccbarstat*D0Dplusccbarstat);
	weightDst = 1/(ccbarDststat*ccbarDststat);

	totalweight = weightD0Dplus + weightDst;

	//Statistical error is square root of 1/(total weight):	
	overallccbarstat = 1/sqrt(totalweight);
	
	//Weighted averages of everything else:
	overallccbar = (D0Dplusccbar * weightD0Dplus + ccbarDst * weightDst)/totalweight;
	overallccbarsysthi = (D0Dplusccbarsysthi * weightD0Dplus + ccbarDstsysthi * weightDst)/totalweight;
	overallccbarsystlo = (D0Dplusccbarsystlo * weightD0Dplus + ccbarDstsystlo * weightDst)/totalweight;
	overallccbarlum =    (D0Dplusccbarlum * weightD0Dplus + ccbarDstlum * weightDst)/totalweight;
	overallccbarextrhi = (D0Dplusccbarextrhi * weightD0Dplus + ccbarDstextrhi * weightDst)/totalweight;
	overallccbarextrlo = (D0Dplusccbarextrlo * weightD0Dplus + ccbarDstextrlo * weightDst)/totalweight;


	//Outputs all cross sections to screen
	cout << "D0:\nVisible =\t" << D0vis << " ± " << D0stat << " (stat.) + " << D0systhi << " - " << D0systlo << " (syst.) microbarn\nExtrapolated =\t" << crossD0 << " ± " << crossD0stat << " (stat.) +" << crossD0systhi <<" -" << crossD0systlo << " (syst.) ± " << crossD0lum << " (lum) ± " << crossD0br << " (br.) +" << crossD0extrhi << " -" << crossD0extrlo << " (extr.) microbarn\nTotal sigma =\t" << ccbarD0 << " ± " << ccbarD0stat << " (stat.) +" << ccbarD0systhi <<" -" << ccbarD0systlo << " (syst.) ± " << ccbarD0lum << " (lum) +" << ccbarD0extrhi << " -" << ccbarD0extrlo << " (extr.) microbarn\n"
	<< "Dst:\nVisible =\t" << Dstvis << " ± " << Dststat << " (stat.) + " << Dstsysthi << " - " << Dstsystlo << " (syst.) microbarn\nExtrapolated =\t" << crossDst << " ± " << crossDststat << " (stat.) +" << crossDstsysthi <<" -" << crossDstsystlo << " (syst.) ± " << crossDstlum << " (lum) ± " << crossDstbr << " (br.) +" << crossDstextrhi << " -" << crossDstextrlo << " (extr.) microbarn\nTotal sigma =\t" << ccbarDst << " ± " << ccbarDststat << " (stat.) +" << ccbarDstsysthi <<" -" << ccbarDstsystlo << " (syst.) ± " << ccbarDstlum << " (lum) +" << ccbarDstextrhi << " -" << ccbarDstextrlo << " (extr.) microbarn\n"
	<< "Dplus:\nVisible =\t" << Dplusvis << " ± " << Dplusstat << " (stat.) + " << Dplussysthi << " - " << Dplussystlo << " (syst.) microbarn\nExtrapolated =\t" << crossDplus << " ± " << crossDplusstat << " (stat.) +" << crossDplussysthi <<" -" << crossDplussystlo << " (syst.) ± " << crossDpluslum << " (lum) ± " << crossDplusbr << " (br.) +" << crossDplusextrhi << " -" << crossDplusextrlo << " (extr.) microbarn\nTotal sigma =\t" << ccbarDplus << " ± " << ccbarDplusstat << " (stat.) +" << ccbarDplussysthi <<" -" << ccbarDplussystlo << " (syst.) ± " << ccbarDpluslum << " (lum) +" << ccbarDplusextrhi << " -" << ccbarDplusextrlo << " (extr.) microbarn\n\nWeighted average =\t" << overallccbar << " ± " << overallccbarstat << " (stat.) + " << overallccbarsysthi << " - " << overallccbarsystlo << " (syst.) ± " << overallccbarlum << " (lum.) + " << overallccbarextrhi << " - " << overallccbarextrlo << " (extr.) µb\n";




}
