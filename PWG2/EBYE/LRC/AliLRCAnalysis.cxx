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

ClassImp(AliLRCAnalysis) 

/******************************************************
 * AliLRCAnalysis class
 ******************************************************/

AliLRCAnalysis::AliLRCAnalysis(): fPrAbs(new TH1D()), fPrRel(new TH1D()), fPrf(new TH1D()), fPrb(new TH1D()), fileHist(new TFile()), fdptb(.0), fEntries(0), fSx((char*)" "), fSy((char*)" "), fxFitMin(.0), fxFitMax(.0), fa_rel(.0), fb_rel(.0), fXi2_rel(.0), fa_abs(.0), fb_abs(.0), fXi2_abs(.0){
//Empty constructor
}

AliLRCAnalysis::AliLRCAnalysis(const AliLRCAnalysis& a):fPrAbs(a.fPrAbs), fPrRel(a.fPrRel), fPrf(a.fPrf), fPrb(a.fPrb), fileHist(a.fileHist), fdptb(a.fdptb), fEntries(a.fEntries), fSx(a.fSx), fSy(a.fSy), fxFitMin(a.fxFitMin), fxFitMax(a.fxFitMax), fa_rel(a.fa_rel), fb_rel(a.fb_rel), fXi2_rel(a.fXi2_rel), fa_abs(a.fa_abs), fb_abs(a.fb_abs), fXi2_abs(a.fXi2_abs){
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
 		fa_rel = a.fa_rel;
		fb_rel = a.fb_rel;
		fXi2_rel = a.fXi2_rel;
		fa_abs = a.fa_abs;
		fb_abs = a.fb_abs; 
		fXi2_abs = a.fXi2_abs;
	}
	return *this;
}


AliLRCAnalysis::~AliLRCAnalysis() {
//Destructor
    delete fPrAbs;
    delete fPrRel;
    delete fileHist;
}

double AliLRCAnalysis::HI2(TH1D *h, double a, double b, double xmin, double xmax) const {
//hi square calculation of approximation 1d hist with ax+b between xmin and xmax
    double trueN = 0;
    double fhi2=0;
    double num;
    double fN=h->GetNbinsX();
    double f = h->GetXaxis()->GetXmin();
    double fdf = (h->GetXaxis()->GetXmax()-f)/fN;
    int fNmin = int((xmin-f)/fdf)+1;
    int fNmax = int((xmax-f)/fdf)+1;
    for(int i=fNmin; i<=fNmax; i++) {
        double fw = h->GetBinError(i);
        if(fw!=0){
            num = b*(i*fdf-fdf/2+f)+a-h->GetBinContent(i);
            fhi2 = fhi2 + (num*num) / (fw*fw);
            trueN++;
        }
    }
    return fhi2/(trueN-2);
}

double AliLRCAnalysis::HI2(TH1D *h, double a, double b) const {
//hi square calculation of approximation 1d hist with ax+b
    double trueN = 0;
    double fhi2=0;
    double num;
    int fN=h->GetNbinsX();
    double f = h->GetXaxis()->GetXmin();
    double fdf = (h->GetXaxis()->GetXmax()-f)/fN;
    int fNmin = 1;
    int fNmax = fN;
    for(int i=fNmin; i<=fNmax; i++) {
        double fw = h->GetBinError(i);
        if(fw!=0){
            num = b*(i*fdf-fdf/2+f)+a-h->GetBinContent(i);
            fhi2 = fhi2 + (num*num) / (fw*fw);
            trueN++;
        }
    }
    return fhi2/(trueN-2);
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
    fxFitMin = mnf-2*sqrt(mnf);
    fxFitMax = mnf+2*sqrt(mnf);

}

void AliLRCAnalysis::SetBinsRange(int binMin, int binMax){
	TH1D* h=fPrf;
    double N=h->GetNbinsX();
    fxFitMin = h->GetXaxis()->GetXmin();
    fxFitMax = h->GetXaxis()->GetXmax();
    double df = (fxFitMax-fxFitMin)/N;
    for(int i=1; i<=N; i++)
	if(h->GetBinContent(i) != 0){
		fxFitMin +=  (i + binMin) * df;
		break;
	}
    for(int i=1; i<=N; i++)
	if(h->GetBinContent(N-i) != 0){
		fxFitMax -= (i + binMax) * df;
		break;
	}
}

bool AliLRCAnalysis::SetFitRange(double xMin, double xMax){
//������������� ������� �������������. ���������� false, ���� xMin ������ xMax. ��������� �� ���������� �� ����, � ���� � ���� ��� �����.((
	if(xMax < xMin){
		return false;
	}
	this->fxFitMin = xMin;
	this->fxFitMax = xMax;
	return true;
}

void AliLRCAnalysis::SetXmin(double xMin){
	fxFitMin = xMin;
}

void AliLRCAnalysis::SetXmax(double xMax){
	fxFitMax = xMax;
}

double AliLRCAnalysis::GetArel(){
	return fa_rel;
}

double AliLRCAnalysis::GetBrel(){
	return fb_rel;
}

double AliLRCAnalysis::GetXi2rel(){
	return fXi2_rel;
}

double AliLRCAnalysis::GetAabs(){
	return fa_abs;
}

double AliLRCAnalysis::GetBabs(){
	return fb_abs;
}

double AliLRCAnalysis::GetXi2abs(){
	return fXi2_abs;
}

void AliLRCAnalysis::Calculate(){
	double mnf;
	fPrAbs->SetStats(0);
    	fPrAbs->Fit("pol1", "0", "", fxFitMin, fxFitMax);
	TF1 *fitt1 = fPrAbs->GetFunction("pol1");
	fa_abs = fitt1->GetParameter(0);
	fb_abs = fitt1->GetParameter(1);
	fXi2_abs = HI2(fPrAbs, fitt1->GetParameter(0), fitt1->GetParameter(1), fxFitMin, fxFitMax);
	
	mnf = fPrf->GetMean(); 
	fPrRel->SetStats(0);
	AliLRCFit *fit1 = new AliLRCFit(fPrRel, fxFitMin/mnf, fxFitMax/mnf);
	TF1 *f1 = new TF1("f1", "[0] + [1]*x", 0, fPrRel->GetXaxis()->GetXmax());
	f1->SetParameter(0,fit1->Geta());
	f1->SetParameter(1,fit1->Getb());
	fPrRel->Fit("f1", "0", "", fxFitMin/mnf, fxFitMax/mnf);
	fa_rel = fit1->Geta();
	fb_rel = fit1->Getb();
	fXi2_rel = fit1->Gethi2();
}

void AliLRCAnalysis::DrawAbs() {
// Draw abs var histrogramm
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
    x1 = fPrAbs->GetXaxis()->GetXmin();
    x2 = fPrAbs->GetXaxis()->GetXmax();
    TPaveText *pt1 = new TPaveText(x1+(x2-x1)/4, y1+(y2-y1)/2, x2-(x2-x1)/4, y2);
    pt1->SetTextSize(0.03);
    
    sprintf(str, "Entries = %i", fEntries);
    pt1->AddText(str);
    sprintf(str, "a = %f #pm %f", fit1->GetParameter(0), fit1->GetParError(0));
    pt1->AddText(str);
    sprintf(str, "b = %f #pm %f", fit1->GetParameter(1), fit1->GetParError(1));
    pt1->AddText(str);
    sprintf(str, "#hat{#chi}^{2} = #chi^{2}/(n-2) = %f", HI2(fPrAbs, fit1->GetParameter(0), fit1->GetParameter(1), fxFitMin, fxFitMax));
    pt1->AddText(str);
    sprintf(str, "<%s> = %f " , fSx, fPrf->GetMean());
    pt1->AddText(str);
    sprintf(str, "<%s> = %f", fSy,  fPrb->GetMean());
    pt1->AddText(str);
    
    sprintf(str, "<<%s>> = %f " , fSx, fPrf->GetRMS());
    pt1->AddText(str);
    sprintf(str, "<<%s>> = %f", fSy,  fPrb->GetRMS());
    pt1->AddText(str);
    if(fdptb){
        sprintf(str, "d%s = %f", fSy,  fdptb);
        pt1->AddText(str);
    }
    
    
    pt1->SetTextAlign(12);
    pt1->SetTextFont(42);
    pt1->SetFillColor(0);
    fPrAbs->DrawCopy();
	pt1->DrawClone();
}

void AliLRCAnalysis::DrawRel() {
// Draw rel var histogramm
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
    x1 = fPrRel->GetXaxis()->GetXmin();
    x2 = fPrRel->GetXaxis()->GetXmax();
    TPaveText *pt1 = new TPaveText(x1+(x2-x1)/4, y1+(y2-y1)/2, x2-(x2-x1)/4, y2);
    pt1->SetTextSize(0.03);
   sprintf(str, "Entries = %i", fEntries);
   pt1->AddText(str);
    sprintf(str, "a = %f #pm %f", fit1->Geta(), fit1->Getda());
    pt1->AddText(str);
    sprintf(str, "b = %f #pm %f", fit1->Getb(), fit1->Getdb());
    pt1->AddText(str);
    sprintf(str, "#hat{#chi}^{2} = #chi^{2}/(n-2) = %f", fit1->Gethi2());
    pt1->AddText(str);
    sprintf(str, "<%s> = %f " , fSx, fPrf->GetMean());
    pt1->AddText(str);
    sprintf(str, "<%s> = %f", fSy, fPrb->GetMean());
    pt1->AddText(str);
    sprintf(str, "<<%s>> = %f " , fSx, fPrf->GetRMS());
    pt1->AddText(str);
    sprintf(str, "<<%s>> = %f", fSy, fPrb->GetRMS());
    pt1->AddText(str);
    if(fdptb){
        sprintf(str, "d%s = %f", fSy,  fdptb);
        pt1->AddText(str);
    }
    pt1->SetTextAlign(12);
    pt1->SetTextFont(42);
    pt1->SetFillColor(0);
    fPrRel->DrawCopy();
    pt1->DrawClone();
}

void AliLRCAnalysis::SetGraphics() const {
// Set root graph style
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    gStyle->SetStatColor(0);
    gStyle->SetPadColor(0);
}

double AliLRCAnalysis::Integral(TH2D* source, Int_t nbin) const {
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
			  fPrAbs->SetBinError(i,ptd*sqrt(nb->GetBinContent(i)/fPrf->GetBinContent(i)));
		}
		fPrRel->SetBinContent(i, fPrAbs->GetBinContent(i)/fPrb->GetMean());
		fPrRel->SetBinError(i,fPrAbs->GetBinError(i)/fPrb->GetMean());	
	}	

}


