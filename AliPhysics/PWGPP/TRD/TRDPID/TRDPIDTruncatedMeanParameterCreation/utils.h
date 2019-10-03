#ifndef UTILS_H 
#define UTILS_H 

#include <TH3.h>
#include <TH1.h>
#include <TAxis.h>


 //======== calculate beta*gamma ======================
Double_t betaGamma(const Double_t p, const Double_t m)
{
    return p/m;
}


//======== calculate pseudorapitidy ==================
Double_t pseudorapidity(const Double_t theta)
{
    Double_t tmpAngle = theta+TMath::Pi()/2;    // global theta is needed
    Double_t tmp = TMath::Tan(tmpAngle/2);
    
    return -TMath::Log(tmp);
}


//======== calculate momentum mean ===================
Double_t momentumMean(const Double_t *p)
{
    //
    // returns mean momentum of 'filled' layers
    //
    
    Double_t sum = 0;
    Double_t norm = 0;
    
    for (Int_t i = 0; i < 6; i++) {
        if (p[i] > 0) {
            norm++;
            sum += p[i];
        }
    }
    
   return sum/norm;

    
    //
    // returns momentum of first filled layer
    //
    /*
    for (Int_t i = 0; i < 6; i++) {
        if (p[i] > 0) return p[i];
    }
    */
}


//======== return maxVal of array ====================
Double_t maxVal(const Double_t *arr)
{
    Double_t maxVal = arr[0];
    
    for (Int_t i = 1; i < (int) (sizeof(arr)/sizeof(arr[0])); i++) {
        if (arr[i] > maxVal) {
            maxVal = arr[i];
        }
    }
        
    return maxVal;
}


//======== return minVal of array ====================
Double_t minVal(const Double_t *arr)
{
    Double_t minVal = arr[0];
    
    for (Int_t i = 1; i < (int) (sizeof(arr)/sizeof(arr[0])); i++) {
        if (arr[i] < minVal) {
            minVal = arr[i];
        }
    }
    
    return minVal;
}


//======== calculate intCharge mean ==================
Double_t intChargeMean(const Double_t *p)
{
    Double_t sum = 0;
    Double_t norm = 0;
    
    for (Int_t i = 0; i < 48; i++) {
        if (i != 0 || i % 8 != 0) continue;
        if (p[i] > 0) {
            norm++;
            sum += p[i];
        }
    }
    
    return sum/norm;
}


//======== get mpv =================================== (or basically any bin content from 1D histogram)
Double_t expMPV(TH1D *hist, const Double_t bg)
{
    return hist->GetBinContent(hist->FindBin(bg));
}


				       

//======== particle momentum/betaGamma cuts ========================

bool electronMomentumCut(Double_t *p,  Float_t nSigma)
{
  // cout << "MomentumCuts" << MomentumCuts << endl;
  if (MomentumCuts) //turn false; // cuts switched off
    {
      cout << "Wieso kTrue" << endl;
      if ((momentumMean(p) > 4 && nSigma  < 0) || (momentumMean(p) > 4 && nSigma > 2)) return true; // original comment: really helpful (no reasons given)
    }
  return false; // no cut
}


bool pionMomentumCut(Double_t *p)
{
  if (MomentumCuts) // return false; // cuts switched off
    {
      //if (momentumMean(p) > 0.8 && momentumMean(p) < 2.5) return true;
      if (betaGamma(momentumMean(p), mPion) >10) return true; // no reasons known (seems strange)
      if (betaGamma(momentumMean(p), mPion) <4) return true;
    }
  return false; // no cut
}


bool protonMomentumCut(Double_t *p)
{
  if (MomentumCuts)  // cuts switched on
    {
      if (momentumMean(p) > 1.2) return true; // no reasons known
      if (betaGamma(momentumMean(p), mProton) <0.9) return true;
    }
  return false; // no cut
}







//======== theta/eta cut =================================
bool thetaCut(const Double_t theta)
{
    // cut on eta boundaries of stacks
    Double_t eta = pseudorapidity(theta);
    
    if ((eta > 0.851)
        || (eta < 0.536 && eta > 0.527)
        || (eta < 0.157 && eta > 0.145)
        || (eta < -0.145 && eta > -0.157)
        || (eta < -0.527 && eta > -0.536)
        || (eta < -0.851)) return true;
    else return false;
    
    /*
    // "arbitrary" fiducial cut -> delta_eta = - delta_theta / sin(theta) => reasonable change delta_eta = +/- 0.1 to the inside
    Double_t eta = pseudorapidity(theta);
    
    if ((eta > 0.841)
        || (eta < 0.546 && eta > 0.517)
        || (eta < 0.167 && eta > 0.135)
        || (eta < -0.135 && eta > -0.167)
        || (eta < -0.517 && eta > -0.546)
        || (eta < -0.841)) return true;
    else return false;
    */
}


//======== local theta cut ===========================
bool localThetaCut(const Double_t *theta, const Int_t nch, Bool_t fTheta)
{
    //Double_t eta[nch];
    Double_t eta;
    
    if(fTheta) {
        eta = (Double_t)pseudorapidity(theta[0]);
    } else {
        eta = (Double_t)theta[0];
    }
    
    bool cut = false;
    
    //for (Int_t i = 0; i < nch; i++) {         // eta[i]
        if ((eta > 0.841)
            || (eta < 0.546     && eta > 0.517)
            || (eta < 0.167     && eta > 0.135)
            || (eta < -0.135    && eta > -0.167)
            || (eta < -0.517    && eta > -0.546)
            || (eta < -0.841)) {
            cut = true;
            //break;
        }
    //}
    
    return cut;
}


//======== get center of supermodule =================
Double_t centerOfSupermodule(const Int_t numberOfSuperModule)
{
    // returns phi in center of supermodule in degrees
    return 10 + numberOfSuperModule*20;
}


//======== cut on SM =================================
bool supermoduleCut(const Int_t * supermoduleArray, const Double_t globalPhi)
{
    bool cut = true;
    
    for (Int_t i = 0; i < (int) (sizeof(supermoduleArray)/sizeof(supermoduleArray[0])); i++) {
        if ((globalPhi*180/TMath::Pi() > centerOfSupermodule(supermoduleArray[i])-10) ||
            (globalPhi*180/TMath::Pi() < centerOfSupermodule(supermoduleArray[i])+10)) {
            cut = false;
            break;
        }
    }
    
    return cut;
}


//======== y cut =====================================
bool yCut(const Double_t *y, const Double_t *yMax)
{
    bool cut = false;
    
    for (Int_t layer = 0; layer < 6; layer++){
        if (TMath::Abs(y[layer]) > yMax[layer]) {
            cut = true;
            break;
        }
        else continue;
    }
    return cut;
}


//======== nSigma cut ================================ (w/ exclusion cut)
/*bool nSigmaCut(const Float_t *nSigmas, const Double_t *nSigmaCuts, const Double_t *nSigmaExclusionCuts, const Int_t pId) {
    bool cut = false;
    
    switch (pId) {
        case 11:
            if ((nSigmas[0] < nSigmaCuts[0]) || (nSigmas[0] > nSigmaCuts[1])) cut = true;
            break;
        case 211:
            if ((nSigmas[1] < nSigmaCuts[2]) || (nSigmas[1] > nSigmaCuts[3])) cut = true;
            break;
        case 2212:
            if (((      nSigmas[0] > nSigmaExclusionCuts[0]) && (nSigmas[0] < nSigmaExclusionCuts[1]))
                || ((   nSigmas[1] > nSigmaExclusionCuts[2]) && (nSigmas[1] < nSigmaExclusionCuts[3]))) cut = true;
            else if ((  nSigmas[2] < nSigmaCuts[4]) || (nSigmas[2] > nSigmaCuts[5])) cut = true;
            break;
    }
    
    return cut;
}*/

//======== nSigma cut ================================ (w/o exclusion cut)
bool nSigmaCut(const Float_t *nSigmas, const Double_t *nSigmaCuts, const Int_t pId) {
    bool cut = false;
    
    switch (pId) {
        case 11:
            if ((nSigmas[0] < nSigmaCuts[0]) || (nSigmas[0] > nSigmaCuts[1])) cut = true;
            break;
        case 211:
            if ((nSigmas[1] < nSigmaCuts[2]) || (nSigmas[1] > nSigmaCuts[3])) cut = true;
            break;
        case 2212:
            if ((nSigmas[2] < nSigmaCuts[4]) || (nSigmas[2] > nSigmaCuts[5])) cut = true;
            break;
    }
    
    return cut;
}


//======== missing slices cut ========================
bool missingSlicesCut(Double_t *slices, Int_t nch) {
    bool cut;

    Int_t emptySlices = 0;
    
    for (Int_t i = 0; i < 48; i++) {
        if (slices[i] <= 0) emptySlices += 1;
    }
    
    if (emptySlices > (6 - nch)*8) {
        cut = true;
    } else {
        cut = false;
    }

    return cut;
}


//======== eta asymmetry correction ==================
//
// OLD: better use 2D correction -> etaCorrFactor
//
Double_t EtaCorrFactorV1(const TF1* mpvDist, const Double_t eta)
{
    Double_t mpvDist_0      = mpvDist->Eval(0);
    Double_t mpvDist_eta    = mpvDist->Eval(eta);
    
    //
    // 1) c(eta) = f(0)/f(eta)      => dE/dx' = dE/dx * c(eta)
    //
    if (mpvDist_eta > 0) {
        Double_t c = mpvDist_0/mpvDist_eta;
        return c;
    } else {
        return 1;
    }
    
    //
    // 2) c(eta) = f(eta) - f(0)    => dE/dx' = dE/dx - c(eta)
    //
    //Double_t c = mpvDist_eta - mpvDist_0;
    //return dEdx - c;
}

//======== eta asymmetry correction ==================
Double_t EtaCorrFactor(const Double_t bg, const Double_t eta, TH2D* etaMap)
{
    Double_t corrFactor = 1;
    
    if (eta >= -0.9 && eta <= 0.9) {
        if (bg >= 0.3 && bg <= 1e4) {
            Double_t tmp = etaMap->GetBinContent(etaMap->FindBin(eta, bg));
            
            if (tmp != 0) corrFactor = tmp;
	    else corrFactor=0;
        }
    }
    
    return corrFactor;
}

//======== eta asymmetry correction for new TPC variable TPC tgl  ==================
Double_t EtaCorrFactorTPCtgl(const Double_t bg, const Double_t eta, TH2D* etaMap)
{
    Double_t corrFactor = 1;
    
    //  if (eta >= -1. && eta <= 1.) { // just changed variable range
    if (bg >= 0.3 && bg <= 1e4) {
      Double_t tmp = etaMap->GetBinContent(etaMap->FindBin(eta, bg));
      if (tmp > 1E-5) corrFactor = tmp;
      else corrFactor=0;
    }
    //  }
    
    return corrFactor;
}

Double_t NclsCorrFactor(const Double_t bg, const Int_t Ncls, TH2D* NclsMap)
{
  Double_t corrFactor = 1.;

  // if (Ncls >=68 && Ncls <= 134) {
    if (bg >= 0.3 && bg <= 1e4) {
      Double_t tmp = NclsMap->GetBinContent(NclsMap->FindBin(Ncls, bg));
      if (tmp != 0) corrFactor = tmp;
      else corrFactor=0;
    }
    // }
  return corrFactor;
}

Double_t CentCorrFactor(const Double_t bg, const Double_t cent, TH2D* CentMap)
{
  Double_t corrFactor = 1.;

  
  if (bg >= 0.3 && bg <= 1e4) {
    Double_t tmp = CentMap->GetBinContent(CentMap->FindBin(cent, bg));
    if (tmp != 0) corrFactor = tmp;
    else corrFactor=0;
  }
  
  return corrFactor;
}

//======== Gaus fit function ========================= (Xianguo)
Double_t RawGaus(const Double_t *xx, const Double_t *par)
{
    return par[0]*TMath::Gaus(xx[0], par[1], TMath::Abs(par[2]), kTRUE);
}


//======== Resolution fit function =================== (Xianguo)
Double_t ResFunc(const Double_t *xx, const Double_t *par)
{
    //return par[0] + par[1] * TMath::Power(xx[0], par[2]);
    //return par[0] + TMath::Abs(par[1]) * TMath::Power(xx[0], par[2]);

    // test
    Double_t p0 = TMath::Abs(par[0]);
    Double_t p1 = TMath::Abs(par[1]);
    
    return TMath::Sqrt(p0*p0 + p1*p1/xx[0]);
}


//======== array of axis bins ========================  (Xianguo)
Double_t * GetAxisArray(TAxis * aa)
{
    // define array of doubles
    TArrayD * xs = new TArrayD(*(aa->GetXbins()));
    
    if (xs->GetSize() == 0) {                           // GetSize() returns number of elements
        const Int_t nbin = aa->GetNbins();
        xs->Set(nbin+1);                                // set number of elements
        for (Int_t ii = 0; ii <= nbin; ii++) {
            xs->SetAt(aa->GetBinUpEdge(ii), ii);        // GetBinUpEdge() returns upper edge of bin
        }
    }
    
    Double_t * bins = xs->GetArray();
    if (!bins) {
        printf("GetAxisArray bins null!! %d\n", xs->GetSize());
        exit(1);
    }
    
    return bins;
}


//======== return chi^2 ==============================  (Xianguo)
Double_t GetFCNChisquare(const Double_t pars[])
{
    Double_t f = 0;
    
    for (Int_t idata = 0; idata < fgFitNData; idata++) {
        const Double_t evaly    = fgFitFunc(&(fgFitX[idata]), pars);
        const Double_t dy       = fgFitY[idata] - evaly;
        
        Double_t dfdx = 0;
        
        //calculated df/dx if necessary
        if (fgFitEx[idata] != 0) {
            const Double_t derix0 = fgFitX[idata] - fgFitEx[idata];
            const Double_t derix1 = fgFitX[idata] + fgFitEx[idata];
            const Double_t deriy0 = fgFitFunc(&derix0, pars);
            const Double_t deriy1 = fgFitFunc(&derix1, pars);
            dfdx = (deriy1-deriy0)/(derix1-derix0);
        }
        f += dy*dy/(fgFitEy[idata]*fgFitEy[idata] + dfdx*fgFitEx[idata]*dfdx*fgFitEx[idata]);
    }
    
    return f;
}


//======== used by Minuit to minimize chi^2 ==========  (Xianguo)
void MinFCNChisquare(Int_t &npar, Double_t *gin, Double_t &f, Double_t *pars, Int_t iflag)
{
    f = GetFCNChisquare(pars);
}


//======== FitKernel =================================  (Xianguo)
Int_t FitKernel(const Int_t npar, Double_t pars[], Double_t errs[], Double_t chi[])
{
    // initialize TMinuit
    TMinuit * mnt = new TMinuit(npar);
    mnt->SetFCN(MinFCNChisquare);
    
    const Double_t errdef = 1;
    mnt->SetErrorDef(errdef);

    mnt->SetPrintLevel(-1); // set output level

    // Choose different Fit Strategy
    Double_t StrArg[1]={2};
    Int_t ErrFlag=-999;
    //if (npar !=8) mnt->mnexcm("SET STR", StrArg, 1, ErrFlag );
    


    Int_t nfree = 0;
    for (Int_t ipar = 0; ipar < npar; ipar++) {
        const Double_t start    = pars[ipar];
        Double_t ss       = TMath::Abs(pars[ipar])  > EPSILON ? TMath::Abs(pars[ipar])*1e-1 : 1e-3;   // old 6e-3
	//important if limits are set -- reasonable for p0 of gaussian fits / but not for the mean value
        const Double_t ll       = 0;
        const Double_t hh       = 0;
        
	if (npar !=8 && ipar ==1) {
	  // ss = TMath::Abs(pars[2])*1e-3; // set stepsize for gaussian mean proportional to RMS
	}

	if (npar == 8) {
	  //	  ss = TMath::Abs(pars[ipar])*1e-3;
	}
        // try parameter limits: -> if fit problematic: fix one parameter to range! (works for ALEPH + TR case, therefore npar == 8)
	/*  if (npar == 8) {
            if (ipar == 0) {
                mnt->DefineParameter(ipar, Form("p%d",ipar), start, ss, 0, 1.5);
            } else {
                mnt->DefineParameter(ipar, Form("p%d",ipar), start, ss, ll,hh);
		}
	     if (ipar == 2) {
	        mnt->DefineParameter(ipar, Form("p%d",ipar), start, ss, 0, 8.);
            } else {
                mnt->DefineParameter(ipar, Form("p%d",ipar), start, ss, ll,hh);
		}
	    
        } else {
            mnt->DefineParameter(ipar, Form("p%d",ipar), start, ss, ll,hh);
	    }*/

	/*	if (npar != 8 && ipar == 1) {
	  mnt->DefineParameter(ipar, Form("p%d", ipar), start, ss, pars[1]-pars[2], pars[1]+pars[2]); // restrict gaussian mean to RMS deviation from mean
	}
        else {
	  mnt->DefineParameter(ipar, Form("p%d",ipar), start, ss, ll,hh);
	  }*/

	mnt->DefineParameter(ipar, Form("p%d",ipar), start, ss, ll,hh);
       
	nfree++;
    }

    if (npar == 8) {
      //for (int j=0; j<7; j++) { mnt->FixParameter(j); }
      //nfree-=7;
    }


    //cout << "initial CHisquare" << GetFCNChisquare(pars) << endl;
    //cout << fMaxIter << " " << fTol << endl;

    // first fit using MIGRAD
    Double_t arg[2] = {Double_t(fMaxIter), fTol};
    Int_t migradflag = -999;
    mnt->mnexcm("MIGRAD", arg, 2, migradflag);
    
    //get par right after MIGRAD before HESSE modifies them, initial pars not changed
    Double_t migp[10];
    Double_t dummy = -999;
    for(Int_t ipar = 0; ipar < npar; ipar++){
        mnt->GetParameter(ipar, migp[ipar], dummy);
    }
    
    // second fit using HESSE    
    Int_t hesseflag = -999;
    mnt->mnexcm("HESSE", arg, 1, hesseflag);
    
    const Int_t kfail = migradflag+hesseflag;
    if (kfail) {
        printf("FitKernel fitting error !! migradflag = %d, hesseflag = %d: keep raw par and set err to -999!\n", migradflag, hesseflag);
        if (errs) {
            for (Int_t ipar = 0; ipar < npar; ipar++) {
                errs[ipar] = -999;
            }
        }
        if (chi) {
            chi[0] = chi[1] = -999;
        }
    }
    else {
        for (Int_t ipar = 0; ipar < npar; ipar++) {
            if (errs) {
                mnt->GetParameter(ipar, dummy, errs[ipar]);
            }
            
            pars[ipar] = migp[ipar];
        }
        if (chi) {
            chi[0] = GetFCNChisquare(pars);
            if (nfree != mnt->GetNumFreePars()) {
                printf("TMinuitFitKernel nfree!=mnt->GetNumFreePars() %d %d\n", nfree, mnt->GetNumFreePars());
                exit(1);
            }
            chi[1] = -nfree;        //need to add nx outside FitKernel
        }
    }
    
    delete mnt;
    
    return kfail;
}


//======== binned fit function using FitKernel =======  (Xianguo)
Int_t BinnedFit(const Int_t nx, const Double_t xdata[], const Double_t ydata[], const Double_t *xerr, const Double_t yerr[], FFunc ffunc, const Int_t npar, Double_t pars[], Double_t errs[], Double_t chi[])
{
    if (nx <= 0) {
        printf("ChisquareFit error nx <= 0 %d\n", nx);
	//  exit(1);
       	return 0;
    }
    
    fgFitNData  = nx;
    fgFitX      = new Double_t[nx];
    fgFitY      = new Double_t[nx];
    fgFitEx     = new Double_t[nx];
    fgFitEy     = new Double_t[nx];
    
    for (Int_t ix = 0; ix < fgFitNData; ix++) {
        fgFitX[ix]  = xdata[ix];
        fgFitY[ix]  = ydata[ix];
	//cout << fgFitX[ix] << "  " << fgFitY[ix]<< endl;
        fgFitEx[ix] = xerr ? xerr[ix] : 0;
        fgFitEy[ix] = yerr[ix];
    }
    
    fgFitFunc = ffunc;
    
    const Int_t kfail = FitKernel(npar, pars, errs, chi);
    
   if (chi) {
        chi[1] += nx;
        printf("\tChisquareFit Chi2/NDOF %.1f/%.0f = %.1f\n\n", chi[0], chi[1], chi[1] ? chi[0]/chi[1] : -999);
    }
    
    
    delete fgFitX;
    delete fgFitY;
    delete fgFitEx;
    delete fgFitEy;
    
    return kfail;
}


//======== chisquare fit using BinnedFit =============  (Xianguo)
//Int_t ChisquareFit(const TH1 *hh, FFunc ffunc, const Int_t npar, Double_t pars[], Double_t errs[], Double_t chi[], const Bool_t kxerr)
Int_t ChisquareFit(const TH1 *hh, FFunc ffunc, const Int_t npar, Double_t pars[], Double_t errs[], Double_t chi[], const Double_t thres, const Bool_t kxerr)
{
    // cout << thres << endl;
    const Int_t nmax = hh->GetNbinsX();
    
    Double_t xdata[nmax];
    Double_t ydata[nmax];
    Double_t xerr[nmax];
    Double_t yerr[nmax];
    
    Int_t nx = 0;
    // cout << hh->GetXaxis()->GetFirst() << " " << hh->GetXaxis()->GetLast() << endl;
    for(Int_t ii = hh->GetXaxis()->GetFirst(); ii <= hh->GetXaxis()->GetLast(); ii++){
      const Double_t ey = hh->GetBinError(ii);
      const Double_t yy = hh->GetBinContent(ii);
      if (ey < EPSILON) {
	if (yy > EPSILON) {
	  printf("ChisquareFit error! %d %e %e\n", ii, ey, yy);
	  exit(1);
	  //	return kFALSE;
	}
	//	cout << "continue" << endl;
	continue;
      }
        
      // cout << "yy " << yy << " thres  " << thres << "nx " << nx << endl;
      if (yy >= thres) {
	const Double_t xx = hh->GetBinCenter(ii);
	const Double_t ex = hh->GetBinWidth(ii)/TMath::Sqrt(12);        // assuming uniform distribution in bin
        
	xdata[nx]   = xx;
	ydata[nx]   = yy;
	xerr[nx]    = ex;
	yerr[nx]    = ey;
	nx++;
      }
      // else cout << yy << " " << thres << endl;
      //cout << "nx" << nx << endl;
    }
    
    return BinnedFit(nx, xdata, ydata, kxerr ? xerr : 0x0, yerr, ffunc, npar, pars, errs, chi);
}


//======== logarithmic binning =======================  (Xianguo)
void BinLogX(TAxis *axis)
{
    // Method for the correct logarithmic binning of histograms
    // copied and modified from AliTPCcalibBase
    
    const Int_t bins    = axis->GetNbins();
    const Double_t from = axis->GetXmin();
    const Double_t to   = axis->GetXmax();
    
    if (from<EPSILON){
        //printf("BinLogX warning xmin < epsilon! nothing done, axis not set. %e\n", from);  exit(1);
        return;
    }
    
    Double_t * new_bins = new Double_t[bins + 1];
    
    new_bins[0] = from;
    const Double_t factor = pow(to/from, 1./bins);
    
    for (int i = 1; i <= bins; i++) {
        new_bins[i] = factor * new_bins[i-1];
    }
    
    axis->Set(bins, new_bins);
    delete [] new_bins;
}


//======== return fit ================================  (Xianguo)
TH1D * GetHfit(const TString hname, FFunc func, const Double_t par[], const Double_t xmin, const Double_t xmax, const Bool_t kbinlogx = kFALSE)
{
    TH1D * h1 = new TH1D(hname, "", 1000, xmin, xmax);
    
    if (kbinlogx) BinLogX(h1->GetXaxis());
    
    for(Int_t ii = 1; ii <= h1->GetNbinsX(); ii++){
        const Double_t xx = h1->GetBinCenter(ii);
        const Double_t yy = func(&xx, par);
        h1->SetBinContent(ii, yy);
    }
    
    return h1;
}


//======== fit slices ================================  (Xianguo)
void FitSlicesY(const TH2D *hInput, TH1D *&hnor, TH1D *&hmpv, TH1D *&hwid, TH1D *&hres, TH1D *&hchi, const Int_t thres, const Double_t thresRMS, TList *ll)
{
  // bin numbers
  const Int_t xFirst  = hInput->GetXaxis()->GetFirst();
  const Int_t xLast   = hInput->GetXaxis()->GetLast();
  const Int_t yFirst  = hInput->GetYaxis()->GetFirst();
  const Int_t yLast   = hInput->GetYaxis()->GetLast();
    
  // number of bins
  const Int_t nBinsX = hInput->GetNbinsX();
  const Int_t nBinsY = hInput->GetNbinsY();
    
  // lower (xMin)/upper (xMax) edges of first/last bin
  const Double_t xMin = hInput->GetXaxis()->GetXmin();
  const Double_t xMax = hInput->GetXaxis()->GetXmax();
  const Double_t yMin = hInput->GetYaxis()->GetXmin();
  const Double_t yMax = hInput->GetYaxis()->GetXmax();
   
  // define histograms
  hnor = new TH1D(Form("%s_amp",hInput->GetName()), "", nBinsX, xMin, xMax);
  hmpv = new TH1D(Form("%s_mpv",hInput->GetName()), "", nBinsX, xMin, xMax);
  hwid = new TH1D(Form("%s_wid",hInput->GetName()), "", nBinsX, xMin, xMax);
  hres = new TH1D(Form("%s_res",hInput->GetName()), "", nBinsX, xMin, xMax);
  hchi = new TH1D(Form("%s_chi",hInput->GetName()), "", nBinsX, xMin, xMax);
    
  // array of bins of hInput->GetXaxis()
 const Double_t *hInputXBins = GetAxisArray(hInput->GetXaxis());
    
  // initialize axis with variable bins
  hnor->GetXaxis()->Set(nBinsX, hInputXBins);
  hmpv->GetXaxis()->Set(nBinsX, hInputXBins);
  hwid->GetXaxis()->Set(nBinsX, hInputXBins);
  hres->GetXaxis()->Set(nBinsX, hInputXBins);
  hchi->GetXaxis()->Set(nBinsX, hInputXBins);
    
  for(Int_t ix = xFirst; ix <= xLast; ix++){
    //cout << ix << endl << ix << endl << ix << endl;


    // temporary slices histogram
    TH1D *htmp = new TH1D(Form("%s_%d", hInput->GetName(), ix), "", nBinsY, yMin, yMax);
        
    //checked, ok
    const Double_t *hInputYBins = GetAxisArray(hInput->GetYaxis());
    //hInput->GetYaxis()->GetXbins()->GetArray();
        
    // if (hInputYBins != 0) set xAxis bins
    if (hInputYBins){
      htmp->GetXaxis()->Set(nBinsY, hInputYBins);
    }
        
    // declare normalized histogram bin content
    Double_t ntot = 0;
      
    for (Int_t iy = yFirst; iy <= yLast; iy++){
      const Double_t be = hInput->GetBinError(ix, iy);
      const Double_t bc = hInput->GetBinContent(ix, iy);
            
      if (be < EPSILON) {
	if(bc > EPSILON) {
	  Printf("FitSlicesY error %d %d %e %e\n", ix, iy, be, bc);
	  exit(1);
	}
	continue;
      }
            
      htmp->SetBinContent(iy, bc);
      htmp->SetBinError(iy, be);
        
      ntot += (bc/be)*(bc/be); // (N/Sqrt(N))^2=Ni
      
        
      //if(be) printf("test %d %d : %f %f %f\n", ix, iy, bc, be, pow(bc/be,2));
    }
        
    // fill _amp histogram
    hnor->SetBinContent(ix,ntot);
    hnor->SetBinError(ix, 0);

    // discard slices w/ nEntries below thres
    if (htmp->Integral() < thres) {
      // cout << "Slice below threshold" << htmp->Integral()  << " " <<thres<< endl; 
      delete htmp;
      continue;
    }
        
    // discard slices w/ RMS above thres
    if (htmp->GetRMS() > thresRMS) {
      //cout << "Slice RMS above threshold" << htmp->GetRMS() << " " << thresRMS << endl;
      delete htmp;
      continue;
    }
        

    // xiangu
    //if (htmp->GetEntries() < thres || htmp->GetRMS() < EPSILON) {
    //    delete htmp;
    //    continue;
    //}
        



    //test htmp->Draw();
    //Double_t pars[10] = {htmp->Integral(0, htmp->GetNbinsX()+1)*htmp->GetBinWidth(1), htmp->GetMean(), htmp->GetRMS()};
    Double_t pars[10] = {htmp->GetMaximum()*(2.5*htmp->GetRMS()), htmp->GetMean(), htmp->GetRMS()}; // More reasonable to devide p0 by approx 2.5
    Double_t errs[10], chi[10];
        
    TH1D * tmpslicefit = 0x0;
        
    // fit histogram slices with gaus
    Int_t MaxBin=htmp->GetMaximumBin();
    Double_t AveragefitThres=htmp->GetBinContent(MaxBin-1)+htmp->GetBinContent(MaxBin)+htmp->GetBinContent(MaxBin+1);


    const Double_t fitThres = AveragefitThres/3.*GaussFrac; // htmp->GetMaximum()*GaussFrac;
    //for (int k=0; k<10; k++) { cout << errs[k] << endl; }

    const Int_t kfail = ChisquareFit(htmp, RawGaus, 3, pars, errs, chi, fitThres, kFALSE);
      
    if (kfail) {
      printf("FitSlicesY ChisquareFit fail! %s %d n %.0f, direct continue\n", htmp->GetName(), kfail, htmp->Integral(0, htmp->GetNbinsX()+1));
      delete htmp;
      continue;
    }
        
    // get fitted slices histogram
    tmpslicefit = GetHfit(Form("%sfit", htmp->GetName()), RawGaus, pars, htmp->GetXaxis()->GetXmin(), htmp->GetXaxis()->GetXmax());
        
    // discard "bad" fits
    if (tmpslicefit->GetRMS() > thresRMS) {
      cout << "Slice RMS above threshold" << tmpslicefit->GetRMS() << " " << thresRMS << endl;
      delete tmpslicefit;
      continue;
    }


    // additional condition for bad fits (hesse has no error message for bad errors if migrad fails. (Florian)
    if (chi[0]/chi[1]>maxChi) {
      //cout << "Bad Fit / Chi Square" << chi[0]/chi[1] << endl;
      //ChiVetoCount++;
      //delete tmpslicefit;
      //continue;
    }
        
    // set par[2] to abs value for width histogram
    pars[2] = TMath::Abs(pars[2]);
      
    // mpv
    hmpv->SetBinContent(ix, pars[1]);
    hmpv->SetBinError(ix, errs[1]);
        
    // width
    hwid->SetBinContent(ix, pars[2]);
    hwid->SetBinError(ix, errs[2]);
        
    // resolution (width/mpv)
    hres->SetBinContent(ix, fabs(pars[1]) > EPSILON ? pars[2]/fabs(pars[1]) : 0);  // if(fabs(pars[1]) > EPSILON) {pars[2]/fabs(pars[1])} else {0}
    hres->SetBinError(ix, fabs(pars[1]) > EPSILON ? errs[2]/fabs(pars[1]) : 0);    // not exact but for MPV>>Width
        
    // reduced chi
    hchi->SetBinContent(ix, chi[1] >= 1 ? chi[0]/chi[1] : 0);
    hchi->SetBinError(ix, 0);
        
    // add slices histograms and fits to list
    if (ll) {
      ll->Add(htmp);
      if(tmpslicefit) ll->Add(tmpslicefit);
    }
    else {
      delete htmp;
      delete tmpslicefit;
    }
  }
    
  TH1 *hhs[] = {hnor, hmpv, hwid, hres, hchi};
  const TString yt[] = {"N", "MPV", "#sigma", "#sigma/MPV", "#chi^{2}/NDOF"};
  const Int_t nh = sizeof(hhs)/sizeof(TH1*);
    
  for(Int_t ii = 0; ii < nh; ii++){
    hhs[ii]->SetYTitle(Form("%s of %s", yt[ii].Data(), hInput->GetYaxis()->GetTitle()));
    hhs[ii]->SetXTitle(hInput->GetXaxis()->GetTitle());
    hhs[ii]->GetYaxis()->SetTitleOffset(hInput->GetYaxis()->GetTitleOffset());
    hhs[ii]->SetTitle(hInput->GetTitle());
  }

  cout << ChiVetoCount << "Slices have been discarded" << endl;
}


//======== SliceAverage  ================================  (FHerrmann)
void SliceAverage(const TH2D *hInput,  TH1D *&hmpv,  TList *ll)
{
  // bin numbers
  const Int_t xFirst  = hInput->GetXaxis()->GetFirst();
  const Int_t xLast   = hInput->GetXaxis()->GetLast();
  const Int_t yFirst  = hInput->GetYaxis()->GetFirst();
  const Int_t yLast   = hInput->GetYaxis()->GetLast();
    
  // number of bins
  const Int_t nBinsX = hInput->GetNbinsX();
  const Int_t nBinsY = hInput->GetNbinsY();
    
  // lower (xMin)/upper (xMax) edges of first/last bin
  const Double_t xMin = hInput->GetXaxis()->GetXmin();
  const Double_t xMax = hInput->GetXaxis()->GetXmax();
  const Double_t yMin = hInput->GetYaxis()->GetXmin();
  const Double_t yMax = hInput->GetYaxis()->GetXmax();
   
  // define histograms
  hmpv = new TH1D(Form("%s_mpv",hInput->GetName()), "", nBinsX, xMin, xMax);
    
  // array of bins of hInput->GetXaxis()
  const Double_t *hInputXBins = GetAxisArray(hInput->GetXaxis());
 
  // initialize axis with variable bins
  hmpv->GetXaxis()->Set(nBinsX, hInputXBins);
    
  for(Int_t ix = xFirst; ix <= xLast; ix++){
    //cout << ix << endl << ix << endl << ix << endl;


    // temporary slices histogram
    TH1D *htmp = new TH1D(Form("%s_%d", hInput->GetName(), ix), "", nBinsY, yMin, yMax);
        
    //checked, ok
    const Double_t *hInputYBins = GetAxisArray(hInput->GetYaxis());
    //hInput->GetYaxis()->GetXbins()->GetArray();
        
    // if (hInputYBins != 0) set xAxis bins
    if (hInputYBins){
      htmp->GetXaxis()->Set(nBinsY, hInputYBins);
    }
        
    // declare normalized histogram bin content
    Double_t ntot = 0;
      
    for (Int_t iy = yFirst; iy <= yLast; iy++){
      const Double_t be = hInput->GetBinError(ix, iy);
      const Double_t bc = hInput->GetBinContent(ix, iy);
            
      if (be < EPSILON) {
	if(bc > EPSILON) {
	  Printf("FitSlicesY error %d %d %e %e\n", ix, iy, be, bc);
	  exit(1);
	}
	continue;
      }
            
      htmp->SetBinContent(iy, bc);
      htmp->SetBinError(iy, be);
        
      ntot += (bc/be)*(bc/be); // (N/Sqrt(N))^2=Ni
      
        
      //if(be) printf("test %d %d : %f %f %f\n", ix, iy, bc, be, pow(bc/be,2));
    }
      
    // mpv
    hmpv->SetBinContent(ix, htmp->GetMean());
    hmpv->SetBinError(ix, 0);
 
  }
    

 


  TH1 *hhs[] = {hmpv};
  const TString yt[] = {"MPV"};
  const Int_t nh = sizeof(hhs)/sizeof(TH1*);
    
  for(Int_t ii = 0; ii < nh; ii++){
    hhs[ii]->SetYTitle(Form("%s of %s", yt[ii].Data(), hInput->GetYaxis()->GetTitle()));
    hhs[ii]->SetXTitle(hInput->GetXaxis()->GetTitle());
    hhs[ii]->GetYaxis()->SetTitleOffset(hInput->GetYaxis()->GetTitleOffset());
    hhs[ii]->SetTitle(hInput->GetTitle());
  }

  cout << ChiVetoCount << "Slices have been discarded" << endl;
}


//======== normalized 2D histogram ===================  (Xianguo)
TH2D* NormalHist(const TH2D *hraw, const Double_t thres, const Bool_t kmax)
{
  TH2D *hh = (TH2D*)hraw->Clone(Form("%snor", hraw->GetName()));
  hh->Scale(0);
    
  const Int_t xFirst  = hh->GetXaxis()->GetFirst();
  const Int_t xLast   = hh->GetXaxis()->GetLast();
  const Int_t yFirst  = hh->GetYaxis()->GetFirst();
  const Int_t yLast   = hh->GetYaxis()->GetLast();
    
  Double_t hmax = -1e10;
  Double_t hmin = 1e10;
  Double_t nent = 0;
    
  for (Int_t ix = xFirst; ix <= xLast; ix++) {
    //if option "e" is specified, the errors are computed. if option "o" original axis range of the taget axes will be kept, but only bins inside the selected range will be filled.
    TH1D * sliceh = hraw->ProjectionY(Form("tmpnormalhist%sx%d", hh->GetName(), ix), ix, ix, "oe");
    const Double_t tot = sliceh->GetEntries();
        
    TH1D * pdfh = 0x0;
        
    if (tot > EPSILON) {
      nent += tot;
      
      Double_t imax = -999;
            
      imax = sliceh->GetBinContent(sliceh->GetMaximumBin());
            
      for (Int_t iy = yFirst; iy <= yLast; iy++) {
	const Double_t cont = kmax ? sliceh->GetBinContent(iy)/imax : pdfh->GetBinContent(iy);          // pdfh is empty...
	const Double_t ierr = kmax ? sliceh->GetBinError(iy)/imax   : pdfh->GetBinError(iy);
        
	if (tot > thres && cont > 0) {
	  hh->SetBinContent(  ix, iy, cont);
	  hh->SetBinError(    ix, iy, ierr);
          
	  if(cont > hmax) hmax = cont;
	  if(cont < hmin) hmin = cont;
	}
      }
    }
    
    delete pdfh;
    delete sliceh;
  }
  
  hh->SetEntries(nent);
  hh->SetMinimum(0.99*hmin);
  hh->SetMaximum(1.1*hmax);
    
  return hh;
}


//======== TR part of parametrization ================  (Xianguo)
Double_t MeanTR(const Double_t * xx,  const Double_t * par)
{
  //
  //  ALEPH+LOGISTIC parametrization for dEdx+TR, in unit of MIP
  //  npar = 4
  //
    
  //  cout << "TR" << endl;
  // for (int i=0; i<4; i++) { cout << par[i] << "\t";}


  const Double_t bg       = xx[0];
  const Double_t gamma    = sqrt(1 + bg*bg);
    
  const Double_t p0 = TMath::Abs(par[1]);
  const Double_t p1 = TMath::Abs(par[2]);
  const Double_t p2 = TMath::Abs(par[3]);
    
  const Double_t zz       = TMath::Log(gamma);
  const Double_t tryield  = p0/(1 + TMath::Exp(-p1 * (zz - p2)));
    
  return par[0] + tryield;
}


//======== ALEPH parametrization =====================  (Xianguo)
Double_t ALEPH(const Double_t * xx,  const Double_t * par)
{
  //
  //ALEPH parametrization for dEdx
  //npar = 5
  //
  
  //cout << "Aleph" << endl;
  //for (int i=0; i<5; i++) { cout << par[i] << "\t";}


  const Double_t bg   = xx[0];
  const Double_t beta = bg/TMath::Sqrt(1. + bg*bg);
    
  const Double_t p0 = TMath::Abs(par[0]);
  const Double_t p1 = TMath::Abs(par[1]);
    
  //after redefining p2 (<0) -> ln P2
  //const Double_t p2 = par[2];
    
  //restore back
  const Double_t p2 = TMath::Abs(par[2]);
  
  const Double_t p3 = TMath::Abs(par[3]);
  const Double_t p4 = TMath::Abs(par[4]);
  
  const Double_t aa = TMath::Power(beta, p3);
    
  //numerically very bad when p2~1e-15, bg~1e3, p4~5.
  const Double_t bb = TMath::Log( p2 + TMath::Power(1./bg, p4) );
    
  //+1 to avoid the singularity when bg<1 and p
  // More reasonable to devide 4>>1
  //very^inf important! with this hesse NOTPOS->OK and no nan or inf in migrad
  //i.e. numerically very stable
  //the reason is when bg<1 ln blows up very fast with p4>>1 and nan/inf appears in migrad search
  //the fitting will adjust the parameters as if bg is not shifted, the fitted curve looks fine!!
  //-- 2012 Aug 8
  //----> fail for 10h, not used, restore back!
  //const Double_t lbg = TMath::Log(bg);
  
  //redefine p2(<0) -> ln P2
  //corresponds to alephFit.C fAlephOPt = 11
  //const Double_t bb = p2 + TMath::Log( 1 + TMath::Exp(-p2-p4*lbg) );
    
  //printf("test----- %f %f -- %f %f %f %f %f --- %f %f %f\n", bg, beta, p0, p1, p2, p3, p4, p0/aa, aa, bb);
    
  return (p1-aa-bb)*p0/aa;
}


//======== dEdx part of parametrization ==============  (Xianguo)
Double_t MeandEdx(const Double_t * xx,  const Double_t * par)
{
    return ALEPH(xx, par);
}


//======== dEdx + TR parametrization =================  (Xianguo)
Double_t MeandEdxTR(const Double_t * xx, const Double_t * pin)
{
    //
    //  ALEPH+LOGISTIC parametrization for dEdx+TR, in unit of MIP
    //  npar = 8 = 3+5
    //
    
    Double_t ptr[4] = {0};
    for (int ii = 0; ii < 3; ii++) {
        ptr[ii+1] = pin[ii];
    }

    return MeanTR(xx, ptr) + MeandEdx(xx, &(pin[3]));
    
    // dE/dx only:
    //return MeandEdx(xx, &(pin[3]));
}


// ============ general Track and layer cuts ==========

bool trackAndLayerCuts(Int_t trdNCh, Int_t trdNCls, Int_t pId) {
  if (trdNCh < nchMin  || trdNCh > nchMax) return true; // layer cut
  // if (tracks.trdNch() != nch) return true; // old version, but problems for fitting each slice separately
  // cout << trdNCh << endl;
  bool cut = false;
  switch (pId) {
  case 11:
    if (trdNCls/trdNCh < ncls[0])  cut = true;
    break;
  case 211:
    if (trdNCls/trdNCh < ncls[1])  cut = true;
    break;
  case 2212:
     if (trdNCls/trdNCh < ncls[2])  cut = true;
    break;
  }

 
  //  if (trdNCls/trdNCh < ncls) return true; // general quality track cuts

  return cut;
}

bool runCut(Int_t run) {
  if (RunPeriod=="LHC13bc") {
  }
  
  if (RunPeriod=="LHC15o")
    {
      if (run==246870 || run== 246855  || run== 246433  || run== 246392  || run== 246391  || run== 246148  || run== 245963  || run== 245793  || run== 245785  || run== 245775  || run== 245766  || run== 245759  || run== 245752 || run==  245738 || run== 245731 || run== 245729) 
	{ // bad quality runs
	  return true;
	};
      if (run==246152 ||run==246984 || run==246758||run==246048||run==246807||run==246181||run==246087||run==246272 ||run==246757 ||run==246434 ||run==246151 ||run==246042 ||run==246488 ||run==246804 ||run==246948 ||run==246982) // HIR >5.5 Part1
	{
	  return false;
	}

      if(run==246846 ||run==246928 ||run==246115 ||run==246945 ||run==246180 ||run==246805||run==246217 ||run==246851 ||run==245683 ||run==246847 ||run==246750 ||run==246271 ||run==246487 ||run==246431 ||run==246844 ||run==246037 ||run==246178 ||run==246980 ||run==246845 ||run==246428 ||run==246424 ||run==245923 ||run==246113 ||run==246036) //HIR > 5.5 kHz Part2
	{
	  return false;
	}
      else
	{
	  return true;
	}
    }
  
  if (RunPeriod=="LHC15n")
    {
      // cout << "LHC15n" << run << endl;
      if (run==244377 || run ==  244364 || run ==  244359 || run ==  244355|| run ==  244351 || run ==  244343 || run ==  244340) // LIR
      
	// if( run == 244628 || run ==  244627|| run== 244626 || run == 244619 || run== 244618 || run == 244617 || run== 244542|| run== 244540||run == 244531||run == 244484|| run ==244483|| run ==244482||run == 244481||run == 244480) 
	{
	  //	cout << run << endl;
	  return false;
	}
      else
	{
	  return true;
	}
  }

  

  return false;
}

// at the moment deactivated because Eta and Y are not available on ESD level. 
bool geometricCuts(Double_t* trdEtaLocal, Int_t trdNCh, Double_t* trdY) {
  // if (localThetaCut(trdEtaLocal, trdNCh, kFALSE)) return true;
  // if (yCut(trdY, (Double_t*)yMaxLayer)) return true;

  return false;
}


void copyFile(const char * Source, const char * Destination) {
  std::ifstream  src(Source, std::ios::binary);
  std::ofstream  dst(Destination,   std::ios::binary);
  dst << src.rdbuf();
}

	
/*      
// For correction Maps old
bool checkSliceContent(TH3D * Hist, Double_t XValue, Double_t YValue) {
  TH1D * tmpHist;
  Int_t XBin, YBin;
  XBin = Hist->GetXaxis()->FindBin(XValue);
  YBin = Hist->GetYaxis()->FindBin(YValue);
  tmpHist = Hist->ProjectionZ("tmpHist", XBin, XBin, YBin, YBin);
  if (tmpHist->Integral() > 100) return false;
  else return true;
}
*/

// For correction Maps old
bool checkSliceContent(TH2D * Hist, Double_t XValue, Double_t YValue) {
  Int_t XBin, YBin;
  XBin = Hist->GetXaxis()->FindBin(XValue);
  YBin = Hist->GetYaxis()->FindBin(YValue);
  if (Hist->GetBinContent(XBin, YBin) > 100000) return false;
  else return true;
}



#endif
