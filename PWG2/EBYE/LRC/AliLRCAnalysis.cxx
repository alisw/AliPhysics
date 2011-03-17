/**************************************************************************
 * Author: Andrey Ivanov.                                           *
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
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch ,
//    Andrey Ivanov (SPbSU-CERN), Igor Altsebeev (SPbSU-CERN) 
//-------------------------------------------------------------------------

#include "AliLRCAnalysis.h"
#include "TFile.h"
#include "AliLRCFit.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPaveText.h"
#include "TF1.h"
#include "math.h"
#include "TStyle.h"
class gStyle;


ClassImp(AliLRCAnalysis) 

/******************************************************
 * AliLRCAnalysis class
 ******************************************************/

AliLRCAnalysis::AliLRCAnalysis():TObject(), fPrAbs(new TH1D()), fPrRel(new TH1D()), fPrf(new TH1D()), fPrb(new TH1D()), fileHist(new TFile()), fdptb(.0), fEntries(0), fSx((char*)" "), fSy((char*)" "), fxFitMin(.0), fxFitMax(.0), fNsigma(2.), farel(.0), fbrel(.0), farelError(.0), fbrelError(.0), fXi2rel(.0), faabs(.0), fbabs(.0), faabsError(.0), fbabsError(.0), fXi2abs(.0),fFitMethod(0){
//Empty constructor
}

AliLRCAnalysis::AliLRCAnalysis(const AliLRCAnalysis& a):TObject(a),fPrAbs(a.fPrAbs), fPrRel(a.fPrRel), fPrf(a.fPrf), fPrb(a.fPrb), fileHist(a.fileHist), fdptb(a.fdptb), fEntries(a.fEntries), fSx(a.fSx), fSy(a.fSy), fxFitMin(a.fxFitMin), fxFitMax(a.fxFitMax), fNsigma(a.fNsigma), farel(a.farel), fbrel(a.fbrel), farelError(a.farelError), fbrelError(a.fbrelError), fXi2rel(a.fXi2rel), faabs(a.faabs), fbabs(a.fbabs),faabsError(a.faabsError), fbabsError(a.fbabsError), fXi2abs(a.fXi2abs),fFitMethod(a.fFitMethod){
//Constructor
}

AliLRCAnalysis& AliLRCAnalysis::operator= (const AliLRCAnalysis& a){
//Operator =
	if(this!=&a){
		TObject::operator= (a);
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
		fNsigma = a.fNsigma;
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
		fFitMethod = a.fFitMethod;
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
    double rms = fPrf->GetRMS();
    fxFitMin = mnf-fNsigma*rms;//sqrt(rms);
    fxFitMax = mnf+fNsigma*rms;//sqrt(rms);//mnf+2*sqrt(rms);
    //delete profX;

}

void AliLRCAnalysis::SetBinsRange(int binMin, int binMax){
//Set the bin range
	TH1D* h=fPrf;
    int n=h->GetNbinsX();
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
bool AliLRCAnalysis::SetFitRangeMin(double xMin){
//Set the fit range min
	fxFitMin = xMin;
	return true;
}

bool AliLRCAnalysis::SetFitRangeMax(double xMax){
//Set the fit range min
	//this->fxFitMax = xMax;

	TH1D* h=fPrf;
	double maxBorder = 0.;

	for ( int binI = h->GetNbinsX(); binI > 0; binI-- )
	{
		if ( h->GetBinContent(binI) != 0 )
		{
			maxBorder = h->GetBinLowEdge(binI+1) ;
			break;
		}	
	}
	fxFitMax = TMath::Max( fxFitMin, TMath::Min( xMax, maxBorder ) );

	return true;
}

void AliLRCAnalysis::SetFullFitRange(){
//Set fitting on full range
	TH1D* h=fPrf;
    fxFitMin = h->GetXaxis()->GetXmin();
	
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
void AliLRCAnalysis::SetNsigma(double nSigma){
//Sets fiting range in forward value RMSs
	fNsigma = nSigma;
	double mnf = fPrf->GetMean();
	double rms = fPrf->GetRMS();
	fxFitMin = mnf-fNsigma*rms;
	fxFitMax = mnf+fNsigma*rms;

//	cout << "### Fit params: N sigma is " << fNsigma 
//		<< ", mean = " << mnf << ", rms = " << rms << ", minFitX = " << fxFitMin << ", maxFitX = " << fxFitMax << endl;
}

void AliLRCAnalysis::SetFitMethod(int id){
//Choose fit method
	fFitMethod = id;
}
double AliLRCAnalysis::GetFitXmin() const
{
	return fxFitMin;
}

double AliLRCAnalysis::GetFitXmax() const 
{
	return fxFitMax;
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
	double mnf = fPrf->GetMean(); 

	// new fit for abs
	AliLRCFit *fit1 = new AliLRCFit(fPrRel, fxFitMin/mnf, fxFitMax/mnf, 1.);
	AliLRCFit *fit2 = new AliLRCFit(fPrAbs, fxFitMin, fxFitMax, 0.);
	farel = fit1->Geta();
	fbrel = fit1->Getb();

	faabs = fit2->Geta();
	fbabs = fit2->Getb();

	if ( fFitMethod==0 )
	{
		faabsError = fit2->Getda();
		fbabsError = fit2->Getdb();

		farelError = fit1->Getda();
		fbrelError = fit1->Getdb();
	}
	if ( fFitMethod==1 )
	{
		faabsError = fit2->Getda1();
		fbabsError = fit2->Getdb1();

		farelError = fit1->Getda1();
		fbrelError = fit1->Getdb1();
	}

	fXi2rel = fit1->Gethi2();
	fXi2abs = fit2->Gethi2();
}

void AliLRCAnalysis::DrawAbs() {
//Draw abs var hist with ALL info
	int * mas = new int [fgkPlotFlags];
	for ( int i = 0; i < fgkPlotFlags; i++ )
		mas[i] = 1;
	DrawAbsPure( mas, 1 );
	delete []mas;
}

void AliLRCAnalysis::DrawAbs( const int * const mas ) {
//Draw abs var hist with REQUESTED BY ARRAY info
	DrawAbsPure( mas, 1 );
}

void AliLRCAnalysis::DrawAbsPure( const int * const mDrawArray, bool drawPaveLabel ) {
	Calculate();
// Draw abs var histrogram
	DrawHist( mDrawArray, drawPaveLabel, faabs, fbabs, faabsError, fbabsError, fPrAbs, 0 );
}

void AliLRCAnalysis::DrawRel() {
//Draw rel var hist with ALL info
	int * mas = new int [fgkPlotFlags];
	for ( int i = 0; i < fgkPlotFlags; i++ )
		mas[i] = 1;
	DrawRelPure( mas, 1 );
	delete []mas;
}

void AliLRCAnalysis::DrawRel( const int * const mas ) {
//Draw rel var hist with REQUESTED BY ARRAY info
	DrawRelPure( mas, 1 );
}

void AliLRCAnalysis::DrawRelPure( const int * const mDrawArray, bool drawPaveLabel ) {
	Calculate();
// Draw rel var histogram
	DrawHist( mDrawArray, drawPaveLabel, farel, fbrel, farelError, fbrelError, fPrRel, 1 );
}
void AliLRCAnalysis::DrawHist( const int * const mDrawArray, bool drawPaveLabel, double aCoef, double bCoef,
		double aCoefError, double bCoefError, TH1D* profToDraw, int histType ) {
// Method called by DrawRelPure or DrawAbsPure to draw corresponding LRC  histo 
	double mnf = fPrf->GetMean();
	double y1, y2, x1, x2;
	Int_t i, n;
	char str[50];
    //mnf = fPrf->GetMean();
	profToDraw->SetStats(0);
    //profToDraw->Fit("pol1", "", "", fxFitMin, fxFitMax);
    //profToDraw->Fit("pol1");
    //TF1 *fit1 = profToDraw->GetFunction("pol1");

	//draw fit line
	TF1 *f1 = 0x0;
	if ( histType == 0 )	{
		f1 = new TF1( "f1", "[0]+[1]*x", fxFitMin, fxFitMax);
		f1->SetLineColor(kRed);
		f1->SetParameters(aCoef,bCoef); 
	}
	else if ( histType == 1 )	{
		f1 = new TF1( "f1", "[0]+[1]*(x-1)", fxFitMin/mnf, fxFitMax/mnf);
		f1->SetLineColor(kGreen);
		f1->SetParameters(aCoef,bCoef); 
	}
	//cout << " set draw params: a=" << aCoef << " b=" << bCoef << endl;

	y1=profToDraw->GetBinContent(1)-profToDraw->GetBinError(1);
	y2=profToDraw->GetBinContent(1)+profToDraw->GetBinError(1);
	n=profToDraw->GetNbinsX();

	for(i=2; i<=n; i++)
	{
        	if(profToDraw->GetBinContent(i)-profToDraw->GetBinError(i)<y1)
            		y1=profToDraw->GetBinContent(i)-profToDraw->GetBinError(i);
        	if(profToDraw->GetBinContent(i)+profToDraw->GetBinError(i)>y2)
            		y2=profToDraw->GetBinContent(i)+profToDraw->GetBinError(i);
	}

	profToDraw->DrawCopy();
	f1->DrawCopy("same");
	
	x1 = profToDraw->GetXaxis()->GetXmin();
	x2 = profToDraw->GetXaxis()->GetXmax();

	if ( drawPaveLabel )
	{	
		int nDatas = 0;
		for ( int j = 0; j < 9; j++)
			if ( mDrawArray[j] ) nDatas++;
		double aXshift = (x2-x1)/7;
		double aYshift = (y2-y1)/20;

		TPaveText *pt1 = new TPaveText(x1+(x2-x1)/2 + aXshift, y1+aYshift, x2-(x2-x1)/6 + aXshift, y1+(y2-y1)/3*2/9*nDatas + aYshift);
		snprintf(str,50, "Events = %i", fEntries);
		//sprintf(str, "Events = %i", fEntries);
		if ( mDrawArray[0] ) pt1->AddText(str);
		snprintf(str, 50,"a = %g #pm %g", GetRoundWithError( aCoef, aCoefError ), GetRoundWithPrecision(aCoefError, 2)); //fit1->GetParameter(0), fit1->GetParError(0));
		//sprintf(str, "a = %g #pm %g", GetRoundWithError( aCoef, aCoefError ), GetRoundWithPrecision(aCoefError, 2)); //fit1->GetParameter(0), fit1->GetParError(0));
		if ( mDrawArray[1] ) pt1->AddText(str);
		snprintf(str, 50,"b = %g #pm %g", GetRoundWithError( bCoef, bCoefError ), GetRoundWithPrecision(bCoefError, 2)); //fit1->GetParameter(1), fit1->GetParError(1));
		//sprintf(str, "b = %g #pm %g", GetRoundWithError( bCoef, bCoefError ), GetRoundWithPrecision(bCoefError, 2)); //fit1->GetParameter(1), fit1->GetParError(1));
		if ( mDrawArray[2] ) pt1->AddText(str);
		snprintf(str, 50,"#hat{#chi}^{2} = #chi^{2}/(n-2) = %g", GetRoundWithPrecision(fXi2abs, 3));
		//sprintf(str, "#hat{#chi}^{2} = #chi^{2}/(n-2) = %g", GetRoundWithPrecision(fXi2abs, 3));
		if ( mDrawArray[3] ) pt1->AddText(str);
		snprintf(str,50, "<%s> = %g " , fSx, GetRoundWithPrecision(fPrf->GetMean(), 3));
		//sprintf(str, "<%s> = %g " , fSx, GetRoundWithPrecision(fPrf->GetMean(), 3));
		if ( mDrawArray[4] ) pt1->AddText(str);
		
		snprintf(str,50, "<%s> = %g", fSy,  GetRoundWithPrecision(fPrb->GetMean(),3));
		//sprintf(str, "<%s> = %g", fSy,  GetRoundWithPrecision(fPrb->GetMean(),3));
		if ( mDrawArray[5] ) pt1->AddText(str);
		
		snprintf(str, 50,"<<%s>> = %g " , fSx, GetRoundWithPrecision(fPrf->GetRMS(), 3));
		//sprintf(str, "<<%s>> = %g " , fSx, GetRoundWithPrecision(fPrf->GetRMS(), 3));
		if ( mDrawArray[6] ) pt1->AddText(str);
		snprintf(str, 50,"<<%s>> = %g", fSy,  GetRoundWithPrecision(fPrb->GetRMS(), 3));
		//sprintf(str, "<<%s>> = %g", fSy,  GetRoundWithPrecision(fPrb->GetRMS(), 3));
		if ( mDrawArray[7] ) pt1->AddText(str);
		
		if ( fdptb ) {
			snprintf(str,50, "d%s = %g", fSy,  GetRoundWithPrecision(fdptb, 3));
			//sprintf(str, "d%s = %g", fSy,  GetRoundWithPrecision(fdptb, 3));
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
		   fPrAbs->SetBinError(i,sqrt(profX->GetBinContent(i)/fPrf->GetBinContent(i)));
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
			  fPrAbs->SetBinError(i,ptd*sqrt(Integral(nb,i))/fPrf->GetBinContent(i));
		}
		fPrRel->SetBinContent(i, fPrAbs->GetBinContent(i)/fPrb->GetMean());
		fPrRel->SetBinError(i,fPrAbs->GetBinError(i)/fPrb->GetMean());	
	}	

}

void AliLRCAnalysis::SetErrors(TH2D* source, const char *name, double ptd, const TProfile* nb){
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
			  fPrAbs->SetBinError(i,ptd*sqrt(nb->GetBinContent(i)/fPrf->GetBinContent(i)));
		}
		fPrRel->SetBinContent(i, fPrAbs->GetBinContent(i)/fPrb->GetMean());
		fPrRel->SetBinError(i,fPrAbs->GetBinError(i)/fPrb->GetMean());	
	}	

}

double AliLRCAnalysis::GetRoundWithError( double value, double error ) const {
//Rounding error and value with DEFAULT precision
	return GetRoundValueErrorPrecision( value, error, 2 );
}
double AliLRCAnalysis::GetRoundWithError( double value, double error, int pres ) const {
//Rounding error and value with REQUESTED precision
	return GetRoundValueErrorPrecision( value, error, pres );
}
double AliLRCAnalysis::GetRoundWithPrecision( double value, int pres ) const {
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
	
	//cout << "Before rounding: " << value << " " << error << endl;

	if ( !error)
	{
		error = value;
		noError = true;
	}
	while ( ((int)error)%10 == 0 && i < 6 )
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
		//cout << "After rounding: " << error << endl;
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
		
		//cout << "After rounding: " << value << " " << error << endl;
		return value; //taking into account ERROR
	}
}

TH1D* AliLRCAnalysis::GetAbsHisto() const
{
//Returns final histo in absolute variables
return fPrAbs;
}
TH1D* AliLRCAnalysis::GetRelHisto() const
{
//Returns final histo in relative variables
return fPrRel;
}

TH1D* AliLRCAnalysis::GetForwardValueDist() const
{
//Returns destribution of value used in forward window
return fPrf;
}
TH1D* AliLRCAnalysis::GetBackwardValueDist() const
{
//Returns destribution of value used in backward window
return fPrb;
}


