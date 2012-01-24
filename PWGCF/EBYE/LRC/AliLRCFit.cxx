/**************************************************************************
 * Author: Andrey Ivanov.                                                 *
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
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch,
//    Andrey Ivanov (SPbSU-CERN), Igor Altsebeev (SPbSU-CERN) 
//-------------------------------------------------------------------------

/* $Id$ */

#include "AliLRCFit.h"
#include "TH1D.h"
#include "math.h"

class TH1D;

ClassImp(AliLRCFit) 


AliLRCFit::AliLRCFit():TObject(),fN (0), fTrueN(0), fNmin (0), fNmax(0), fS1(.0), fSz(.0), fSfz(.0), fSf(.0), fSf2(.0), fhi2(.0), fw(.0), fz(.0), f(.0), fnum(.0), fdf(.0), fdelta(.0), fa(.0), fb(.0), fda(.0), fdb(.0), fda1(.0), fdb1(.0), fxmin(.0), fxmax(.0){
//Empty constructor
}

AliLRCFit::AliLRCFit(TH1D * const h,double xShift):TObject(),fN (0), fTrueN(0), fNmin (0), fNmax(0), fS1(.0), fSz(.0), fSfz(.0), fSf(.0), fSf2(.0), fhi2(.0), fw(.0), fz(.0), f(.0), fnum(.0), fdf(.0), fdelta(.0), fa(.0), fb(.0), fda(.0), fdb(.0), fda1(.0), fdb1(.0), fxmin(.0), fxmax(.0){
    //Constructor make fit of 1d histogramm
AliLRCFit(h,h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),xShift);
}

AliLRCFit::AliLRCFit(TH1D * const h, double xmin, double xmax,double xShift):TObject(),fN (0), fTrueN(0), fNmin (0), fNmax(0), fS1(.0), fSz(.0), fSfz(.0), fSf(.0), fSf2(.0), fhi2(.0), fw(.0), fz(.0), f(.0), fnum(.0), fdf(.0), fdelta(.0), fa(.0), fb(.0), fda(.0), fdb(.0), fda1(.0), fdb1(.0), fxmin(.0), fxmax(.0){
     //Constructor make fit of 1d histogramm between xmin and xmax
	fNmin=h->GetXaxis()->FindBin(xmin);
	fNmax=h->GetXaxis()->FindBin(xmax);
	//double xShift=0;
	for(int i=fNmin; i<=fNmax; i++){
	fw = h->GetBinError(i);
	f= h->GetBinCenter(i)-xShift;
	fz= h->GetBinContent(i);
	if(fw){
		fTrueN++;
		fS1  += 1.0/(fw*fw);
		fSz  += fz/(fw*fw);
		fSf  += f/(fw*fw);
		fSfz += f*fz/(fw*fw);
		fSf2 += (f*f)/(fw*fw);
		}
	}
	if(fTrueN<2)return;
        fdelta = fS1*fSf2 - fSf*fSf;
        fa = (fSf2*fSz - fSf*fSfz)/fdelta;
	fb = (fS1*fSfz - fSf*fSz)/fdelta;
 	fda = sqrt(fSf2/fdelta);
 	fdb = sqrt(fS1/fdelta);
        fdb1 = 0;
        fda1 = 0;
	fhi2 = 0;
        for(int i=fNmin; i<=fNmax; i++){
	fw = h->GetBinError(i);
	f= h->GetBinCenter(i)-xShift;
	fz= h->GetBinContent(i);
	if(fw){
		double fDletaZ2=fz-fa-fb*f;
		fDletaZ2*=fDletaZ2;
		fhi2 += fDletaZ2/(fw*fw);
		fda1 += ( fDletaZ2/(fdelta*fdelta) ) * ((fSf2 - fSf*f)/ (fw*fw))*((fSf2 - fSf*f)/ (fw*fw));
		fdb1 += ( fDletaZ2/(fdelta*fdelta) ) * ((fS1*f - fSf)/ (fw*fw))*((fS1*f - fSf)/ (fw*fw));
		}
	}
	fda1 = sqrt( fda1 );
	fdb1 = sqrt( fdb1 );
	fhi2 = fTrueN > 2 ? fhi2 / (fTrueN-2) : -1;
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
double AliLRCFit::GetFitRange() const
{
//Returns range between xmin and xmax
return (fxmax-fxmin);
}



