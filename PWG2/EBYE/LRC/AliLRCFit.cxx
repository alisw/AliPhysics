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
//    it makes fit of the 1d histogramm
//    calculates ax+b coefficients with error and hi square
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------

/* $Id$ */

#include "AliLRCFit.h"


ClassImp(AliLRCFit) 


AliLRCFit::AliLRCFit():fN (0), fTrueN(0), fNmin (0), fNmax(0), fS1(.0), fSz(.0), fSfz(.0), fSf(.0), fSf2(.0), fhi2(.0), fw(.0), fz(.0), f(.0), fnum(.0), fdf(.0), fdelta(.0), fa(.0), fb(.0), fda(.0), fdb(.0), fda1(.0), fdb1(.0), fxmin(.0), fxmax(.0){
//Empty constructor
}

AliLRCFit::AliLRCFit(TH1D* h):fN (0), fTrueN(0), fNmin (0), fNmax(0), fS1(.0), fSz(.0), fSfz(.0), fSf(.0), fSf2(.0), fhi2(.0), fw(.0), fz(.0), f(.0), fnum(.0), fdf(.0), fdelta(.0), fa(.0), fb(.0), fda(.0), fdb(.0), fda1(.0), fdb1(.0), fxmin(.0), fxmax(.0){
    //Constructor make fit of 1d histogramm
    fxmin = h->GetXaxis()->GetXmin();
    fxmax = h->GetXaxis()->GetXmax();
    fN=h->GetNbinsX();
    fNmin=1;
    fNmax=fN;
    fnum = h->GetXaxis()->GetXmin();
    fdf = (h->GetXaxis()->GetXmax()-fnum)/fN;
    for(int i=1; i<=fN; i++){
        f=i*fdf-fdf/2+fnum;
        fw = h->GetBinError(i);
        if(fw!=0){
            fz = h->GetBinContent(i);
            fS1 = fS1 + 1/(fw*fw);
            fSz = fSz + fz/(fw*fw);
            fSfz = fSfz + f*fz/(fw*fw);
            fSf = fSf + f/(fw*fw);
            fSf2 = fSf2 + (f*f)/(fw*fw);
        }
    }
    fdelta = fS1*fSf2 - fSf*fSf;
    fb = (fS1*fSfz - fSf*fSz)/fdelta;
    fa = ((fSf2-fSf)*fSz - (fSf-fS1)*fSfz)/fdelta;
    fda = sqrt((fS1+fSf2-2*fSf)/fdelta);
    fdb = sqrt(fS1/fdelta);
    fdb1 = 0;
    fda1 = 0;
    f = h->GetXaxis()->GetXmin();
    for(int i=1; i<=fN; i++){
        f=i*fdf-fdf/2+fnum;
        fw = h->GetBinError(i);
        if(fw!=0){
            fz = h->GetBinContent(i);
            fdb1 = fdb1 + ((fS1*f - fSf)*(fS1*f - fSf)/(fw*fw)) * ((fz-fa-fb*(f-1))*(fz-fa-fb*(f-1))/(fw*fw));
            fda1 = fda1 + (((fSf2-fSf)-(fSf-fS1)*f)*((fSf2-fSf)-(fSf-fS1)*f)/(fw*fw)) * ((fz-fa+fb-fb*f)*(fz-fa+fb-fb*f)/(fw*fw));
            fhi2 = fhi2 + ((fz-(fa-fb)-fb*f) * (fz-(fa-fb)-fb*f)) / (fw*fw);
            fTrueN++;
        }
    }
    fdb1 = sqrt(fdb1/(fdelta*fdelta));
    fda1 = sqrt(fda1/(fdelta*fdelta));
    fhi2 = fhi2 / (fTrueN-2);
}

AliLRCFit::AliLRCFit(TH1D *h, double xmin, double xmax):fN (0), fTrueN(0), fNmin (0), fNmax(0), fS1(.0), fSz(.0), fSfz(.0), fSf(.0), fSf2(.0), fhi2(.0), fw(.0), fz(.0), f(.0), fnum(.0), fdf(.0), fdelta(.0), fa(.0), fb(.0), fda(.0), fdb(.0), fda1(.0), fdb1(.0), fxmin(.0), fxmax(.0){
     //Constructor make fit of 1d histogramm between xmin and xmax
    fxmin = xmin;
    fxmax = xmax;
    fN=h->GetNbinsX();
    fnum = h->GetXaxis()->GetXmin();
    fdf = (h->GetXaxis()->GetXmax()-fnum)/fN;
    fNmin=int((xmin-f)/fdf)+1;
    fNmax=int((xmax-f)/fdf)+1;
    for(int i=fNmin; i<=fNmax; i++){
        f=i*fdf-fdf/2+fnum;
        fw = h->GetBinError(i);
        if(fw!=0){
            fz = h->GetBinContent(i);
            fS1 = fS1 + 1/(fw*fw);
            fSz = fSz + fz/(fw*fw);
            fSfz = fSfz + f*fz/(fw*fw);
            fSf = fSf + f/(fw*fw);
            fSf2 = fSf2 + (f*f)/(fw*fw);
        }
    }
    fdelta = fS1*fSf2 - fSf*fSf;
    fb = (fS1*fSfz - fSf*fSz)/fdelta;
    fa = ((fSf2-fSf)*fSz - (fSf-fS1)*fSfz)/fdelta;
    fda = sqrt((fS1+fSf2-2*fSf)/fdelta);
    fdb = sqrt(fS1/fdelta);
    fdb1 = 0;
    fda1 = 0;
    for(int i=fNmin; i<=fNmax; i++){
        f=i*fdf-fdf/2+fnum;
        fw = h->GetBinError(i);
        if(fw!=0){
            fz = h->GetBinContent(i);
            fdb1 = fdb1 + ((fS1*f - fSf)*(fS1*f - fSf)/(fw*fw)) * ((fz-fa-fb*(f-1))*(fz-fa-fb*(f-1))/(fw*fw));
            fda1 = fda1 + (((fSf2-fSf)-(fSf-fS1)*f)*((fSf2-fSf)-(fSf-fS1)*f)/(fw*fw)) * ((fz-fa+fb-fb*f)*(fz-fa+fb-fb*f)/(fw*fw));
            fhi2 = fhi2 + ((fz-(fa-fb)-fb*f) * (fz-(fa-fb)-fb*f)) / (fw*fw);
            fTrueN++;
        }
    }
    fdb1 = sqrt(fdb1/(fdelta*fdelta));
    fda1 = sqrt(fda1/(fdelta*fdelta));
    fhi2 = fhi2 / (fTrueN-2);
}

AliLRCFit::~AliLRCFit() {
}

double AliLRCFit::Geta() const {return fa;}
double AliLRCFit::Getb() const {return fb;}
double AliLRCFit::Getda() const {return fda;}
double AliLRCFit::Getdb() const {return fdb;}
double AliLRCFit::Getda1() const {return fda1;}
double AliLRCFit::Getdb1() const {return fdb1;}
double AliLRCFit::Gethi2() const {return fhi2;}
double AliLRCFit::Getf() const {return f;}
double AliLRCFit::Getxmin() const {return fxmin;}
double AliLRCFit::Getxmax() const {return fxmax;}
int AliLRCFit::GetN() const {return fN;}




