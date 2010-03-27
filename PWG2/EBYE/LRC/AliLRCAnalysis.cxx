/**************************************************************************
 * Author: Panos Christakoglou.                                           *
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
//-------------------------------------------------------------------------
//    Description: 
//    This class include into LRC library for Long-Range Correlation analysis
//    it is base class for NN, PtN, PtPt
//    implements base methods for thees classes
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------

/* $Id$ */

//-------------------------------------------------------------------------
//         LRC library for Long-Range Correlation analysis
//
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------

#include "AliLRCAnalysis.h"
#include "Riostream.h"
#include "TFile.h"
#include "AliLRCFit.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"

class gStyle;
class math;

ClassImp(AliLRCAnalysis) 

/******************************************************
 * AliLRCAnalysis class
 ******************************************************/

AliLRCAnalysis::AliLRCAnalysis(): fPrAbs(new TH1D()), fPrRel(new TH1D()), fPrf(new TH1D()), fPrb(new TH1D()), fileHist(new TFile()), fdptb(.0), fEntries(0), fSx((char*)" "), fSy((char*)" "), fxFitMin(.0), fxFitMax(.0), farel(.0), fbrel(.0), farelError(0.), fbrelError(0.), fXi2rel(.0), faabs(.0), fbabs(.0), faabsError(0.), fbabsError(0.), fXi2abs(.0){
//Empty constructor
}

AliLRCAnalysis::AliLRCAnalysis(const AliLRCAnalysis& a):fPrAbs(a.fPrAbs), fPrRel(a.fPrRel), fPrf(a.fPrf), fPrb(a.fPrb), fileHist(a.fileHist), fdptb(a.fdptb), fEntries(a.fEntries), fSx(a.fSx), fSy(a.fSy), fxFitMin(a.fxFitMin), fxFitMax(a.fxFitMax), farel(a.farel), fbrel(a.fbrel), farelError(a.farelError), fbrelError(a.fbrelError), fXi2rel(a.fXi2rel), faabs(a.faabs), fbabs(a.fbabs),  faabsError(a.faabsError), fbabsError(a.fbabsError), fXi2abs(a.fXi2abs){
//Constructor
}

AliLRCAnalysis& AliLRCAnalysis::operator= (const AliLRCAnalysis& a){
//Operator =
	if(this!=&a){
		fPrAbs = a.fPrAbs;
		fPrRel = a.fPrRel;
		fPrf = a.fPrf;
		fPrb = a.fPrb;
		fileHist = a.fileHist;
		fSx = a.fSx;
		fSy = a.fSy;
		fdptb = a.fdptb;
		fEntries = a.fEntries;
		fxFitMin = a.fxFitMin;
		fxFitMax = a.fxFitMax;
 		farel = a.farel;
		fbrel = a.fbrel;
 		farelError = a.farelError;
 		fbrelError = a.fbrelError;
		fXi2rel = a.fXi2rel;
		faabs = a.faabs;
		fbabs = a.fbabs; 
		faabsError = a.faabsError;
		fbabsError = a.fbabsError;
		fXi2abs = a.fXi2abs;
	}
	return *this;
}


AliLRCAnalysis::~AliLRCAnalysis() {
//Destructor
    delete fPrAbs;
    delete fPrRel;
    delete fileHist;
	delete fPrf;
	delete fPrb;
}

double AliLRCAnalysis::HI2(TH1D * const h, double a, double b, double xmin, double xmax) const {
//hi square calculation of approximation 1d hist with ax+b between xmin and xmax
    int trueN = 0;
    double fhi2=0;
    double num;
    double fN=h->GetNbinsX();
    double f = h->GetXaxis()->GetXmin();
    double fdf = (h->GetXaxis()->GetXmax()-f)/fN;
    int fNmin = int((xmin-f)/fdf)+1;
    int fNmax = int((xmax-f)/fdf)+1;
    for(int i=fNmin; i<=fNmax; i++) {
        double fw = h->GetBinError(i);
        if(fw){
            num = b*(i*fdf-fdf/2+f)+a-h->GetBinContent(i);
            fhi2 = fhi2 + (num*num) / (fw*fw);
            trueN++;
        }
    }
	//cout << "trueN " << trueN << endl;
    return trueN > 2 ? fhi2/(trueN-2.) : -1;
}

double AliLRCAnalysis::HI2(TH1D * const h, double a, double b) const {
//hi square calculation of approximation 1d hist with ax+b
    int trueN = 0;
    double fhi2=0;
    double num;
    int fN=h->GetNbinsX();
    double f = h->GetXaxis()->GetXmin();
    double fdf = (h->GetXaxis()->GetXmax()-f)/fN;
    int fNmin = 1;
    int fNmax = fN;
    for(int i=fNmin; i<=fNmax; i++) {
        double fw = h->GetBinError(i);
        if(fw){
            num = b*(i*fdf-fdf/2+f)+a-h->GetBinContent(i);
            fhi2 = fhi2 + (num*num) / (fw*fw);
            trueN++;
        }
    }
	//cout << "trueN2 " << trueN << endl;
    return trueN > 2 ? fhi2/(trueN-2.) : -1;
}


void AliLRCAnalysis::CreateHist(char *name, char *nameAbs, char *nameRel, char *atitleF, char *atitleB, char *rtitleF, char *rtitleB, TH2D* sourceHist) {
//Create absolute and relation var histogramm
    TProfile* profX = (TProfile*) sourceHist->ProfileX(name, 1, sourceHist->GetNbinsY());
    fEntries = (int) sourceHist->GetEntries();
    fPrAbs = new TH1D(nameAbs, profX->GetTitle(), profX->GetXaxis()->GetNbins(), profX->GetXaxis()->GetXmin(), profX->GetXaxis()->GetXmax());
    fPrAbs->SetOption("E");
    fPrRel = new TH1D(nameRel, profX->GetTitle(), profX->GetXaxis()->GetNbins(), profX->GetXaxis()->GetXmin()/sourceHist->ProjectionX()->GetMean(), profX->GetXaxis()->GetXmax()/sourceHist->ProjectionX()->GetMean());
    fPrRel->SetOption("E");
    fPrAbs->GetXaxis()->SetTitle(atitleF);
    fPrAbs->GetYaxis()->SetTitle(atitleB);
    fPrRel->GetXaxis()->SetTitle(rtitleF);
    fPrRel->GetYaxis()->SetTitle(rtitleB);
    fPrf = (TH1D*) sourceHist->ProjectionX();
    fPrb = (TH1D*) sourceHist->ProjectionY();
    fPrf->GetXaxis()->SetTitle(atitleF);
    fPrf->GetYaxis()->SetTitle("Tracks");
    fPrb->GetXaxis()->SetTitle(atitleB);
    fPrb->GetYaxis()->SetTitle("Tracks");
    fSx = atitleF;
    fSy = atitleB;
    double mnf = fPrf->GetMean();
    fxFitMin = mnf-2*TMath::Sqrt(mnf);
    fxFitMax = mnf+2*TMath::Sqrt(mnf);
    //delete profX;

}

void AliLRCAnalysis::SetBinsRange(int binMin, int binMax){
//Set the bin range
	TH1D* h=fPrf;
    Int_t n=h->GetNbinsX();
    fxFitMin = h->GetXaxis()->GetXmin();
    fxFitMax = h->GetXaxis()->GetXmax();
    double df = (fxFitMax-fxFitMin)/n;
    for(int i=1; i<=n; i++)
	if(h->GetBinContent(i) != 0){
		fxFitMin +=  (i + binMin) * df;
		break;
	}
    for(int i=1; i<=n; i++)
	if(h->GetBinContent(n-i) != 0){
		fxFitMax -= (i + binMax) * df;
		break;
	}
}

bool AliLRCAnalysis::SetFitRange(double xMin, double xMax){
//Set the fit range
	if(xMax < xMin){
		return false;
	}
	this->fxFitMin = xMin;
	this->fxFitMax = xMax;
	return true;
}
void AliLRCAnalysis::SetFullFitRange(){
//Set fitting on full range
	TH1D* h=fPrf;
    fxFitMin = h->GetXaxis()->GetXmin();
    //fxFitMax = h->GetXaxis()->GetXmax();
	
	for ( int binI = h->GetNbinsX(); binI > 0; binI-- )
	{
		if ( h->GetBinContent(binI) != 0 )
		{
			fxFitMax = h->GetBinLowEdge(binI+1) ;
			break;
		}	
	}
}


void AliLRCAnalysis::SetXmin(double xMin){
	fxFitMin = xMin;
}

void AliLRCAnalysis::SetXmax(double xMax){
	fxFitMax = xMax;
}

double AliLRCAnalysis::GetArel() const {
	return farel;
}

double AliLRCAnalysis::GetBrel() const {
	return fbrel;
}

double AliLRCAnalysis::GetArelError() const {
	return farelError;
}

double AliLRCAnalysis::GetBrelError() const {
	return fbrelError;
}

double AliLRCAnalysis::GetXi2rel() const {
	return fXi2rel;
}

double AliLRCAnalysis::GetAabs() const {
	return faabs;
}

double AliLRCAnalysis::GetBabs() const {
	return fbabs;
}

double AliLRCAnalysis::GetAabsError() const {
	return faabsError;
}

double AliLRCAnalysis::GetBabsError() const {
	return fbabsError;
}
double AliLRCAnalysis::GetXi2abs() const {
	return fXi2abs;
}

void AliLRCAnalysis::Calculate(){
//Calculate all
	double mnf;
	fPrAbs->SetStats(0);
    	fPrAbs->Fit("pol1", "0", "", fxFitMin, fxFitMax);
	TF1 *fitt1 = fPrAbs->GetFunction("pol1");
	faabs = fitt1->GetParameter(0);
	fbabs = fitt1->GetParameter(1);
	faabsError = fitt1->GetParError(0);
	fbabsError = fitt1->GetParError(1);	
	fXi2abs = HI2(fPrAbs, fitt1->GetParameter(0), fitt1->GetParameter(1), fxFitMin, fxFitMax);
	
	mnf = fPrf->GetMean(); 
	fPrRel->SetStats(0);
	AliLRCFit *fit1 = new AliLRCFit(fPrRel, fxFitMin/mnf, fxFitMax/mnf);
	TF1 *f1 = new TF1("f1", "[0] + [1]*x", 0, fPrRel->GetXaxis()->GetXmax());
	f1->SetParameter(0,fit1->Geta());
	f1->SetParameter(1,fit1->Getb());
	fPrRel->Fit("f1", "0", "", fxFitMin/mnf, fxFitMax/mnf);
	farel = fit1->Geta();
	fbrel = fit1->Getb();
	farelError = fit1->Getda();
	fbrelError = fit1->Getdb();
	fXi2rel = fit1->Gethi2();
}

void AliLRCAnalysis::DrawAbs() {
//Draw abs var hist with ALL info
	int * mas = new int [N_PL_FLAGS];
	for ( int i = 0; i < N_PL_FLAGS; i++ )
		mas[i] = 1;
	DrawAbsPure( mas, 1 );
	delete []mas;
}

void AliLRCAnalysis::DrawAbs( int * mas ) {
//Draw abs var hist with REQUESTED BY ARRAY info
	DrawAbsPure( mas, 1 );
}

void AliLRCAnalysis::DrawAbsPure( const int * const mDrawArray, bool drawPaveLabel ) {
// Draw abs var histrogram
    //double mnf;
    double y1, y2, x1, x2;
    Int_t i, n;
    char str[50];
    //mnf = fPrf->GetMean();
    fPrAbs->SetStats(0);
    fPrAbs->Fit("pol1", "", "", fxFitMin, fxFitMax);
    //fPrAbs->Fit("pol1");
    TF1 *fit1 = fPrAbs->GetFunction("pol1");
    y1=fPrAbs->GetBinContent(1)-fPrAbs->GetBinError(1);
    y2=fPrAbs->GetBinContent(1)+fPrAbs->GetBinError(1);
    
    n=fPrAbs->GetNbinsX();
    for(i=2; i<=n; i++){
        if(fPrAbs->GetBinContent(i)-fPrAbs->GetBinError(i)<y1)
            y1=fPrAbs->GetBinContent(i)-fPrAbs->GetBinError(i);
        if(fPrAbs->GetBinContent(i)+fPrAbs->GetBinError(i)>y2)
            y2=fPrAbs->GetBinContent(i)+fPrAbs->GetBinError(i);
    }
 	fPrAbs->DrawCopy();
	
	x1 = fPrAbs->GetXaxis()->GetXmin();
    x2 = fPrAbs->GetXaxis()->GetXmax();
    

	if ( drawPaveLabel )
	{	
		int nDatas = 0;
		for ( int j = 0; j < 9; j++)
			if ( mDrawArray[j] ) nDatas++;
		double aXshift = (x2-x1)/7;
		double aYshift = (y2-y1)/20;

		TPaveText *pt1 = new TPaveText(x1+(x2-x1)/2 + aXshift, y1+aYshift, x2-(x2-x1)/6 + aXshift, y1+(y2-y1)/3*2/9*nDatas + aYshift);
	
		sprintf(str, "Entries = %i", fEntries);
		if ( mDrawArray[0] ) pt1->AddText(str);
		sprintf(str, "a = %g #pm %g", GetRoundWithError(fit1->GetParameter(0), fit1->GetParError(0)), GetRoundWithPrecision(fit1->GetParError(0), 2)); //fit1->GetParameter(0), fit1->GetParError(0));
		if ( mDrawArray[1] ) pt1->AddText(str);
		sprintf(str, "b = %g #pm %g", GetRoundWithError(fit1->GetParameter(1), fit1->GetParError(1)), GetRoundWithPrecision(fit1->GetParError(1), 2)); //fit1->GetParameter(1), fit1->GetParError(1));
		if ( mDrawArray[2] ) pt1->AddText(str);
		sprintf(str, "#hat{#chi}^{2} = #chi^{2}/(n-2) = %g", GetRoundWithPrecision(HI2(fPrAbs, fit1->GetParameter(0), fit1->GetParameter(1), fxFitMin, fxFitMax), 3));
		if ( mDrawArray[3] ) pt1->AddText(str);
		sprintf(str, "<%s> = %g " , fSx, GetRoundWithPrecision(fPrf->GetMean(), 3));
		if ( mDrawArray[4] ) pt1->AddText(str);
		
		sprintf(str, "<%s> = %g", fSy,  GetRoundWithPrecision(fPrb->GetMean(),3));
		if ( mDrawArray[5] ) pt1->AddText(str);
		
		sprintf(str, "<<%s>> = %g " , fSx, GetRoundWithPrecision(fPrf->GetRMS(), 3));
		if ( mDrawArray[6] ) pt1->AddText(str);
		sprintf(str, "<<%s>> = %g", fSy,  GetRoundWithPrecision(fPrb->GetRMS(), 3));
		if ( mDrawArray[7] ) pt1->AddText(str);
		
		if ( fdptb ) {
			sprintf(str, "d%s = %g", fSy,  GetRoundWithPrecision(fdptb, 3));
			if ( mDrawArray[8] ) pt1->AddText(str);
		}
		
		pt1->SetTextAlign(12);
		pt1->SetTextFont(42);
		pt1->SetFillColor(4000);
		//pt1->SetFillStyle(4100);
		pt1->SetShadowColor(4000);
		pt1->SetBorderSize(0);
		
		pt1->DrawClone("same");
	}
}

void AliLRCAnalysis::DrawRel() {
//Draw rel var hist with ALL info
	int * mas = new int [N_PL_FLAGS];
	for ( int i = 0; i < N_PL_FLAGS; i++ )
		mas[i] = 1;
	DrawRelPure( mas, 1 );
	delete []mas;
}

void AliLRCAnalysis::DrawRel( int * mas ) {
//Draw rel var hist with REQUESTED BY ARRAY info
	DrawRelPure( mas, 1 );
}

void AliLRCAnalysis::DrawRelPure( const int * const mDrawArray, bool drawPaveLabel ) {
// Draw rel var histogram
	double mnf;
	double y1, y2, x1, x2;
	Int_t i, n;
	char str[50];
		
	mnf = fPrf->GetMean();  
		
	fPrRel->SetStats(0);
	AliLRCFit *fit1 = new AliLRCFit(fPrRel, fxFitMin/mnf, fxFitMax/mnf);
	TF1 *f1 = new TF1("f1", "[0] + [1]*x", 0, fPrRel->GetXaxis()->GetXmax());
	f1->SetParameter(0,fit1->Geta());
	f1->SetParameter(1,fit1->Getb());
	fPrRel->Fit("f1", "", "", fxFitMin/mnf, fxFitMax/mnf);
    y1=fPrRel->GetBinContent(1)-fPrRel->GetBinError(1);
    y2=fPrRel->GetBinContent(1)+fPrRel->GetBinError(1);
    n=fPrRel->GetNbinsX();
    for(i=2; i<=n; i++) {
        if(fPrRel->GetBinContent(i)-fPrRel->GetBinError(i)<y1)
            y1=fPrRel->GetBinContent(i)-fPrRel->GetBinError(i);
        if(fPrRel->GetBinContent(i)+fPrRel->GetBinError(i)>y2)
            y2=fPrRel->GetBinContent(i)+fPrRel->GetBinError(i);
    }
	fPrRel->DrawCopy();

    x1 = fPrRel->GetXaxis()->GetXmin();
    x2 = fPrRel->GetXaxis()->GetXmax();
	
	if ( drawPaveLabel )
	{	
		int nDatas = 0;
		for ( int j = 0; j < 9; j++)
			if ( mDrawArray[j] ) nDatas++;
		double aXshift = (x2-x1)/7;
		double aYshift = (y2-y1)/20;
		
		TPaveText *pt1 = new TPaveText(x1+(x2-x1)/2 + aXshift, y1+aYshift, x2-(x2-x1)/6 + aXshift, y1+(y2-y1)/3*2/9*nDatas + aYshift);
		sprintf(str, "Entries = %i", fEntries);
		if ( mDrawArray[0] ) pt1->AddText(str);
		sprintf(str, "a = %g #pm %g", GetRoundWithError(fit1->Geta(), fit1->Getda()), GetRoundWithPrecision(fit1->Getda(), 2)); //fit1->GetParameter(0), fit1->GetParError(0));
		if ( mDrawArray[1] ) pt1->AddText(str);
		sprintf(str, "b = %g #pm %g", GetRoundWithError(fit1->Getb(), fit1->Getdb()), GetRoundWithPrecision(fit1->Getdb(), 2)); //fit1->GetParameter(1), fit1->GetParError(1));
		if ( mDrawArray[2] ) pt1->AddText(str);
		sprintf(str, "#hat{#chi}^{2} = #chi^{2}/(n-2) = %g", GetRoundWithPrecision(fit1->Gethi2(), 3));
		if ( mDrawArray[3] ) pt1->AddText(str);
		sprintf(str, "<%s> = %g " , fSx, GetRoundWithPrecision(fPrf->GetMean(), 3));
		if ( mDrawArray[4] ) pt1->AddText(str);
		
		sprintf(str, "<%s> = %g", fSy,  GetRoundWithPrecision(fPrb->GetMean(),3));
		if ( mDrawArray[5] ) pt1->AddText(str);
		
		sprintf(str, "<<%s>> = %g " , fSx, GetRoundWithPrecision(fPrf->GetRMS(), 3));
		if ( mDrawArray[6] ) pt1->AddText(str);
		sprintf(str, "<<%s>> = %g", fSy,  GetRoundWithPrecision(fPrb->GetRMS(), 3));
		if ( mDrawArray[7] ) pt1->AddText(str);
		
		if ( fdptb ) {
			sprintf(str, "d%s = %g", fSy,  GetRoundWithPrecision(fdptb, 3));
			if ( mDrawArray[8] ) pt1->AddText(str);
		}
	
		pt1->SetTextAlign(12);
		pt1->SetTextFont(42);
		pt1->SetFillColor(4000);
		//pt1->SetFillStyle(4100);
		pt1->SetShadowColor(4000);
		pt1->SetBorderSize(0);
		
		pt1->DrawClone();//"s(0,0)");
	}
}

void AliLRCAnalysis::SetGraphics() const {
// Set root graph style
	TStyle tempSt;
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    gStyle->SetStatColor(0);
    gStyle->SetPadColor(0);
}

double AliLRCAnalysis::Integral(TH2D* const source, Int_t nbin) const {
// calculate the integrall for x bin and y bins of 2d histogramm
    double sum = 0;
    for(Int_t i = 1; i<=source->GetNbinsY(); i++) {
        sum += source->GetYaxis()->GetBinCenter(i)*source->GetBinContent(nbin, i);
    }
    return sum;
}

void AliLRCAnalysis::SetErrors(TH2D* source, const char *name){
// Calculate errors for NN
	TProfile* profX = (TProfile*) source->ProfileX(name, 1, source->GetNbinsY());
	for(int i = 0; i < profX->GetNbinsX(); i++)		
	{
		fPrAbs->SetBinContent(i, profX->GetBinContent(i));
		if(fPrf->GetBinContent(i)!=0)
		{
		   fPrAbs->SetBinError(i,TMath::Sqrt(profX->GetBinContent(i)/fPrf->GetBinContent(i)));
		}
		fPrRel->SetBinContent(i, fPrAbs->GetBinContent(i)/fPrb->GetMean());
		fPrRel->SetBinError(i,fPrAbs->GetBinError(i)/fPrb->GetMean());	
	}	
}

void AliLRCAnalysis::SetErrors(TH2D* source, const char *name, double ptd, TH2D* nb){
//Calculate arrors for ptn and ptpt
	TProfile* profX = (TProfile*) source->ProfileX(name, 1, source->GetNbinsY());
	fdptb = ptd;
	double pt;
	for(int i = 0; i < profX->GetNbinsX(); i++)		
	{
		fPrAbs->SetBinContent(i, profX->GetBinContent(i));
		if(fPrf->GetBinContent(i)!=0)
		{
			  pt = profX->GetBinContent(i);
			  fPrAbs->SetBinError(i,ptd*TMath::Sqrt(Integral(nb,i))/fPrf->GetBinContent(i));
		}
		fPrRel->SetBinContent(i, fPrAbs->GetBinContent(i)/fPrb->GetMean());
		fPrRel->SetBinError(i,fPrAbs->GetBinError(i)/fPrb->GetMean());	
	}	

}

void AliLRCAnalysis::SetErrors(TH2D* source, const char *name, double ptd, TProfile* nb){
//Calculate arrors for ptn and ptpt
	TProfile* profX = (TProfile*) source->ProfileX(name, 1, source->GetNbinsY());
	fdptb = ptd;
	double pt;
	for(int i = 0; i < profX->GetNbinsX(); i++)		
	{
		fPrAbs->SetBinContent(i, profX->GetBinContent(i));
		if(fPrf->GetBinContent(i)!=0)
		{
			  pt = profX->GetBinContent(i);
			  fPrAbs->SetBinError(i,ptd*TMath::Sqrt(nb->GetBinContent(i)/fPrf->GetBinContent(i)));
		}
		fPrRel->SetBinContent(i, fPrAbs->GetBinContent(i)/fPrb->GetMean());
		fPrRel->SetBinError(i,fPrAbs->GetBinError(i)/fPrb->GetMean());	
	}	

}

double AliLRCAnalysis::GetRoundWithError( double value, double error ){
//Rounding error and value with DEFAULT precision
	return GetRoundValueErrorPrecision( value, error, 2 );
}
double AliLRCAnalysis::GetRoundWithError( double value, double error, int pres ){
//Rounding error and value with REQUESTED precision
	return GetRoundValueErrorPrecision( value, error, pres );
}
double AliLRCAnalysis::GetRoundWithPrecision( double value, int pres ){
//Rounding error and value with requested precision
	return GetRoundValueErrorPrecision( value, 0, pres );
}

double AliLRCAnalysis::GetRoundValueErrorPrecision( double value, double error, int pres ) const {
	//Rounding error and value with requested precision
	//value == value, error == error
	//if single argument(value=0) - calculate without errors:

	int i = 0;
	double order = 1;
	bool noError = false;
	pres -= 1;
	
	cout << "Before rounding: " << value << " " << error << endl;

	if ( !error)
	{
		error = value;
		noError = true;
	}
	while ( ((int)error)%10 == 0 || i == 10 )
	{
		i++;
		error*=10;
		order*=10;
	}

	for ( int j = 0; j < pres; j++ )
	{
		error*=10;
		order*=10;
		i++;
	}
	int adding = 0;
	if ( ((int)(error*10))%10 > 4 && ((int)(error*10))%10 != 9 ) //trouble: if we round 19 to 20 - zero disappeares!
		adding = 1;
	error = (double)((int)error + adding)/order;
	
	if ( noError )
	{
		cout << "After rounding: " << error << endl;
		return error;
	}
	else
	{
		for ( int j = 0; j < i; j++ )
			value*=10;
		
		adding = 0;
		if ( ((int)(value*10))%10 > 4 && ((int)(value*10))%10 != 9 )
			adding = 1;
		value = (double)((int)value + adding)/order;
		
		cout << "After rounding: " << value << " " << error << endl;
		return value; //taking into account ERROR
	}
}

