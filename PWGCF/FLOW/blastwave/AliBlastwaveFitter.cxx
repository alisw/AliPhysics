#include<stdio.h>

#include"AliBlastwaveFitter.h"

#include"TMath.h"
#include"TH1.h"
#include"TGraphErrors.h"
#include"TGraphAsymmErrors.h"

ClassImp(AliBlastwaveFitter);

AliBlastwaveFitter::AliBlastwaveFitter(const char *name) :
    TNamed(name,name),
    fNfunction(0),
    fMinuit(new TMinuit()),
    fMinos(kFALSE)
{
}
//------------------------------------------------------------------------------
AliBlastwaveFitter::AliBlastwaveFitter() :
    TNamed("BlastwaveFitter","BlastwaveFitter"),
    fNfunction(0),
    fMinuit(new TMinuit()),
    fMinos(kFALSE)
{  
}
//------------------------------------------------------------------------------
AliBlastwaveFitter::~AliBlastwaveFitter(){
}
//------------------------------------------------------------------------------
Int_t AliBlastwaveFitter::CheckAvailability(){
    Int_t something = 2;
    Int_t consistency = 0;
    Int_t nparameters = -1;
    for(Int_t i=0;i < fNfunction;i++){
	if(!fFunc[i]){
	    fSpectraFlag[i] = 0;
	    fV2Flag[i] = 0;
	    continue;
	}
	if(fFunc[i]->GetSpectraFit() && fFunc[i]->GetSpectrumObj()){
	    TString classe(fFunc[i]->GetSpectrumObj()->ClassName());
	    if(classe.Contains("TH1")) fSpectraFlag[i] = 1;
	    else if(classe.Contains("TGraphE")) fSpectraFlag[i] = 2;
	    else fSpectraFlag[i] = 0;
	}
	else fSpectraFlag[i] = 0;
	if(fFunc[i]->GetV2Fit() && fFunc[i]->GetV2Obj()){
	    TString classe(fFunc[i]->GetV2Obj()->ClassName());
	    if(classe.Contains("TH1")) fV2Flag[i] = 1;
	    else if(classe.Contains("TGraphErr")) fV2Flag[i] = 2;
	    else if(classe.Contains("TGraphAsymmErr")) fV2Flag[i] = 3;
	    else fV2Flag[i] = 0;
	}
	else fV2Flag[i] = 0;

	if(fSpectraFlag[i] || fV2Flag[i]) something = 0;
	else continue;

	Int_t npar = fFunc[i]->GetNpar();
	if(nparameters == -1) nparameters = npar;
	else if(npar != nparameters) consistency = 1;
    }

    return (consistency + something);
}
//------------------------------------------------------------------------------
Int_t AliBlastwaveFitter::PrepareToFit(){
    Int_t check = CheckAvailability();

    Int_t npar=0;

    if(check%2){
	printf("Some problems with number of parameters are found\n");
	return check;
    }
    if((check/2)%2){
	printf("Nothing to fit\n");
	return check;
    }
    printf("Summary of fit functions and histos\n");

    Bool_t kFlowAvailable=kFALSE;

    for(Int_t i=0;i < fNfunction;i++){
	if(!fFunc[i]) continue;
	printf("------------------------------------------------------------------\n");
	printf("%2i) SpectraFlag = %i -- V2Flag = %i for \"%s\" (class:%s)  with mass = %6.3f GeV/c^2\n",i,fSpectraFlag[i],fV2Flag[i],fFunc[i]->GetName(),fFunc[i]->ClassName(),fFunc[i]->GetMass());
	if(fSpectraFlag[i]){
	    printf("   data spectra:\"%s\" (class:%s)\n",fFunc[i]->GetSpectrumObj()->GetName(),fFunc[i]->GetSpectrumObj()->ClassName());
	    npar = fFunc[i]->GetNpar();
	}
	if(fV2Flag[i]){
	    printf("   data v2     :\"%s\" (class:%s)\n",fFunc[i]->GetV2Obj()->GetName(),fFunc[i]->GetV2Obj()->ClassName());
	    npar = fFunc[i]->GetNpar();
	    kFlowAvailable=kTRUE;	
	}
    }
    printf("------------------------------------------------------------------\n");

    fgNparReal = npar;

    fgNDGF = 0;
    fgNfunctionCurrent = fNfunction;
    for(Int_t i=0;i < fNfunction;i++){
	fgFuncC[i] = fFunc[i];
	fgSpectraFlagC[i] = fSpectraFlag[i];
	fgV2FlagC[i] = fV2Flag[i];
    }
    // prepare minuit
    if(fgNparReal >= 25){
	fMinuit->BuildArrays(fgNparReal);
    }
    fMinuit->SetFCN(AliBlastwaveFitter::FCN);
    Double_t arglist[3];
    Int_t ierflg = 0;
    arglist[0] = 1;
    fMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

    // par limits of the medium
    for(Int_t i=0;i < npar;i++){
      fMinuit->mnparm(i/*#par*/, fFunc[0]->GetParName(i)/*name*/, fFunc[0]->GetParStart(i)/*startvalue*/, fFunc[0]->GetParStep(i)/*step*/, fFunc[0]->GetParMin(i)/*min*/,fFunc[0]->GetParMax(i)/*max*/, ierflg);
    }
    
    if(!kFlowAvailable){ // if there is no flow histos
	for(Int_t i=0;i < fNfunction;i++){
	    if(fSpectraFlag[i]){
		fFunc[i]->SwitchOffFlow(fMinuit);
		i = fNfunction;
	    }
	}
    }

    // set strategy
    arglist[0] = 1; // !!!
    fMinuit->mnexcm("SET STRATEGY", arglist, 1, ierflg);
  
    return 0;
}
//------------------------------------------------------------------------------
Int_t AliBlastwaveFitter::Fit(){
    Int_t npar = fFunc[0]->GetNpar();
    /* start MIGRAD minimization */
    Double_t arglist[20];
    Int_t ierflg = 0;
    arglist[0] = 500000;
    arglist[1] = 1.;
    fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    /* set strategy */
    arglist[0] = 2;
    fMinuit->mnexcm("SET STRATEGY", arglist, 1, ierflg);
    
    /* start MIGRAD minimization */
    arglist[0] = 500000;
    arglist[1] = 1.;
    fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    /* start IMPROVE minimization */
    arglist[0] = 500000;
    fMinuit->mnexcm("IMPROVE", arglist, 1, ierflg);
    
    /* start MINOS */
    if(fMinos){
      arglist[0] = 500000;
      for(Int_t i=0;i < npar;i++){
  	arglist[i+1] = i;
      }
      fMinuit->mnexcm("MINOS", arglist, npar+1, ierflg);
    }

    /* print results */
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    fMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    fMinuit->mnprin(4, amin);

    for(Int_t i=0;i < fgNfunctionCurrent;i++){
	fgFuncC[i]->Terminate();
    }

    return 0;
}
//------------------------------------------------------------------------------
void AliBlastwaveFitter::FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    Double_t chi = 0., val, vale, pull,pt;
    Int_t iNorm = 0,nparFunc;

    if(gin[0]==1){}

//    printf("%i) %f %f %f %f %f %f %f\n",npar,par[0],par[1],par[2],par[3],par[4],par[5],par[6]);


    fgNDGF=0;  
  
    for(Int_t i=0;i < fgNfunctionCurrent;i++){
	if(!fgFuncC[i]) continue;
	
	// Set parameter in the current function
  	nparFunc = fgFuncC[i]->GetNpar();
	for(Int_t j=0;j<nparFunc;j++)
	    fgFuncC[i]->SetParameter(j,par[j]);

	// Set normalization in the current function
	if(fgSpectraFlagC[i]){// && nparFunc+iNorm < npar){
	    fgFuncC[i]->SetNormalization();//par[nparFunc+iNorm]);
	    iNorm++;
	}

	if(fgSpectraFlagC[i]){
	    // spectra comparison
	    TH1 *hForFit = (TH1 *) fgFuncC[i]->GetSpectrumObj();
	    TGraphErrors *gForFit = (TGraphErrors *) fgFuncC[i]->GetSpectrumObj();

	    if(fgSpectraFlagC[i]==1){
		for (Int_t ibin = 1; ibin <= hForFit->GetNbinsX(); ibin++) {
		    pt = hForFit->GetBinCenter(ibin);
		    if(pt < fgFuncC[i]->GetMaxPt() && pt > fgFuncC[i]->GetMinPt()){
			val = hForFit->GetBinContent(ibin);
			vale = hForFit->GetBinError(ibin);
			if(vale>0 && TMath::Abs(val) > vale+0.00001){
			  pull = (val - fgFuncC[i]->EvalYield(pt))/vale;
			  chi += pull * pull;
			  fgNDGF++;
			}
		    }
		} 
	    }
	    else if(fgSpectraFlagC[i]==2){
		for (Int_t ibin = 0; ibin < gForFit->GetN(); ibin++) {
		    pt = gForFit->GetX()[ibin];
		    if(pt < fgFuncC[i]->GetMaxPt() && pt > fgFuncC[i]->GetMinPt()){
			val = gForFit->GetY()[ibin];
			vale = gForFit->GetEY()[ibin];
			if(vale>0){
			  pull = (val - fgFuncC[i]->EvalYield(pt))/vale;
			  chi += pull * pull;
			  fgNDGF++;
			}
		    }		    
		}		
	    }
	}
	if(fgV2FlagC[i]){
	    // v2 comparison
	    TH1 *hForFit = (TH1 *) fgFuncC[i]->GetV2Obj();
	    TGraphErrors *gForFit = (TGraphErrors *) fgFuncC[i]->GetV2Obj();
	    TGraphAsymmErrors *gForFitA = (TGraphAsymmErrors *) fgFuncC[i]->GetV2Obj();

	    if(fgV2FlagC[i]==1){
		for (Int_t ibin = 1; ibin <= hForFit->GetNbinsX(); ibin++) {
		    pt = hForFit->GetBinCenter(ibin);
		    if(pt < fgFuncC[i]->GetMaxPt() && pt > fgFuncC[i]->GetMinPt()){
			val = hForFit->GetBinContent(ibin);
			vale = hForFit->GetBinError(ibin);
			if(vale>0 && TMath::Abs(val) > vale+0.00001){
			  pull = (val - fgFuncC[i]->EvalV2(pt))/vale;
			  chi += pull * pull;
			  fgNDGF++;
			}
		    }
		} 
	    }
	    else if(fgV2FlagC[i]==2){
		for (Int_t ibin = 0; ibin < gForFit->GetN(); ibin++) {
		    pt = gForFit->GetX()[ibin];
		    if(pt < fgFuncC[i]->GetMaxPt() && pt > fgFuncC[i]->GetMinPt()){
			val = gForFit->GetY()[ibin];
			vale = gForFit->GetEY()[ibin];
			if(vale>0){
			  pull = (val - fgFuncC[i]->EvalV2(pt))/vale;
			  chi += pull * pull;
			  fgNDGF++;
			}
		    }		    
		}		
	    }
	    else if(fgV2FlagC[i]==3){
		for (Int_t ibin = 0; ibin < gForFitA->GetN(); ibin++) {
		    pt = gForFit->GetX()[ibin];
		    if(pt < fgFuncC[i]->GetMaxPt() && pt > fgFuncC[i]->GetMinPt()){
			val = gForFit->GetY()[ibin];
			if(val - fgFuncC[i]->EvalV2(pt) > 0) vale = gForFit->GetEYlow()[ibin];
			else vale = gForFit->GetEYhigh()[ibin];
			if(vale>0){
			  pull = (val - fgFuncC[i]->EvalV2(pt))/vale;
			  chi += pull * pull;
			  fgNDGF++;
			}
		    }		    
		}		
	    }
	}
    }

    fgNDGF -= npar + iNorm;

    if(fgNDGF > 0) fgChi2 = chi/fgNDGF;

    f = chi;
    iflag = 0;
}
//------------------------------------------------------------------------------
TGraph*  AliBlastwaveFitter::DoContour(Int_t np,Int_t ip1,Int_t ip2,Float_t nsigma){

  TGraph *gCorr=NULL;
  
  Double_t par1,par2,err;
  fMinuit->GetParameter(ip1, par1,err);
  fMinuit->GetParameter(ip2, par2,err);

  Int_t trial = 0;
  Double_t scale = 1;
  while((! gCorr || gCorr->GetN() < np/4+5) && trial < 10){
    fMinuit->SetErrorDef(nsigma*nsigma*scale); //3 sigma contour
    gCorr = (TGraph *) fMinuit->Contour(np,ip1,ip2);
    scale *= 0.5;
    trial++;
  }
  scale *= 2;
  scale = TMath::Sqrt(scale);

  if(scale < 0.7 && gCorr){
    Double_t *x = gCorr->GetX();
    Double_t *y = gCorr->GetY();
    for(Int_t i=0;i < gCorr->GetN();i++){
      gCorr->SetPoint(i,par1+(x[i]-par1)/scale,par2+(y[i]-par2)/scale);
    }
  }

  if(trial > 1) printf("Contour realized with %3.1f sigma instead of %3.1f and then rescaled\n",nsigma*TMath::Sqrt(scale),nsigma);

  if(gCorr) gCorr->SetMarkerStyle(20);
  
  return gCorr;

}
 //------------------------------------------------------------------------------
TGraph* AliBlastwaveFitter::DoContourBetaT(Int_t np,Int_t iBoostOrBeta,Int_t iT,Float_t nsigma){


  TGraph *gCorr=NULL;
  
  Double_t par1,par2,err;
  fMinuit->GetParameter(iBoostOrBeta, par1,err);
  fMinuit->GetParameter(iT, par2,err);

  const Int_t maxnpar = 20;
  Double_t par[maxnpar];

  Int_t npar = fMinuit->GetNumPars();
  for(Int_t i=0;i < TMath::Min(npar,maxnpar);i++){
    fMinuit->GetParameter(i, par[i],err);
  }

  Int_t trial = 0;
  Double_t scale = 1;
  while((! gCorr || gCorr->GetN() < np/4+5) && trial < 10){
    fMinuit->SetErrorDef(nsigma*nsigma*scale); //3 sigma contour
    gCorr = (TGraph *) fMinuit->Contour(np,iBoostOrBeta,iT);
    scale *= 0.5;
    trial++;
  }
  scale *= 2;
  scale = TMath::Sqrt(scale);

  if(gCorr){
    Double_t *x = gCorr->GetX();
    Double_t *y = gCorr->GetY();
    for(Int_t i=0;i < gCorr->GetN();i++){
      Float_t meanboost = par1+(x[i]-par1)/scale;
      par[iBoostOrBeta] = meanboost;
      par[iT] = y[i];
      Float_t meanbeta = fgFuncC[0]->GetMeanBeta(par);
      gCorr->SetPoint(i,meanbeta,par2+(y[i]-par2)/scale);
    }
  }

  if(trial > 1) printf("Contour realized with %3.1f sigma instead of %3.1f and then rescaled\n",nsigma*TMath::Sqrt(scale),nsigma);

  if(gCorr) gCorr->SetMarkerStyle(20);
  
  return gCorr;

}
 //------------------------------------------------------------------------------
TGraph* AliBlastwaveFitter::ConvertContourFromBoostToBeta(TGraph *g,Int_t iBoostOrBeta,Int_t iT){

  if(!g) return NULL;

  TGraph *gCorr=new TGraph(g->GetN());
  
  const Int_t maxnpar = 20;
  Double_t par[maxnpar],err;

  Int_t npar = fMinuit->GetNumPars();
  for(Int_t i=0;i < TMath::Min(npar,maxnpar);i++){
    fMinuit->GetParameter(i, par[i],err);
  }

  if(gCorr){
    Double_t *x = g->GetX();
    Double_t *y = g->GetY();
    for(Int_t i=0;i < gCorr->GetN();i++){
      par[iBoostOrBeta] = x[i];
      par[iT] = y[i];
      Float_t meanbeta = fgFuncC[0]->GetMeanBeta(par);
      gCorr->SetPoint(i,meanbeta,y[i]);
    }
  }

  if(gCorr) gCorr->SetMarkerStyle(20);
  
  return gCorr;

}
