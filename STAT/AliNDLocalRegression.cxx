/**************************************************************************
 * Copyright(c) 2006-07, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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
//                Implementation of the AliNDLocalRegression class
//-------------------------------------------------------------------------

/*
  Related task: https://alice.its.cern.ch/jira/browse/ATO-193

  Algorithm secription - see: 
  Kernel_smoother: Local polynomial regression   
  http://en.wikipedia.org/w/index.php?title=Kernel_smoother&oldid=627785784
  

  Formally, the local polynomial regression is computed by solving a weighted least square problem.
  Weights are provided as a width of the gausian kernel. 
  Local fit parameters are computed on the grid defined by axis set defiend by THn.
  For example use please check UnitTest:
  .L $ALICE_ROOT/../src/STAT/test/AliNDLocalRegressionTest.C+
  //
  Init:
  AliNDLocalRegression *pfitNDIdeal=0; 
  pfitNDIdeal->SetHistogram((THn*)(hN->Clone()));
  pfitNDIdeal->SetCuts(3,0.8,1);                  // outlier rejection setting see /AliNDLocalRegressionTest.C:UnitTestGaussNoisePlusOutliers() for motivation
  pfitNDIdeal->MakeFit(treeIn, "val:err", "xyz0:xyz1","Entry$%2==1", "0.05:0.05","2:2",0.001);

  Usage: 

  Double_t xyz[2]={1,2}
  pfitNDIdeal->Eval(xyz);

  In TFormulas:
  pfitNDGaus0->AddVisualCorrection(pfitNDGaus0,2); 
  pfitNDGaus1->AddVisualCorrection(pfitNDGaus0,3);
  treeIn->Draw("(AliNDLocalRegression::GetCorrND(3,xyz0,xyz1)-AliNDLocalRegression::GetCorrND(2,xyz0,xyz1))/sqrt(AliNDLocalRegression::GetCorrNDError(3,xyz0,xyz1)**2+AliNDLocalRegression::GetCorrNDError(2,xyz0,xyz1)**2)>>pullsGaus01(200,-20,20)","","");


  To do:
     1.) Statistical error of the local interpolation ignores Gaussian kernel weights 
         errors are overestimated - find a proper mathematical formula to estimate statistical error of estimator
     2.) Implent regularization for smoothing  - requesting approximate smoothnes in values and derivative
     

  author: marian.ivanov@cern.ch
*/

#include <TVectorD.h>

#include "AliNDLocalRegression.h"
#include "AliLog.h"

#include "THn.h"
#include "TObjString.h"
#include "TTreeStream.h"
#include "AliMathBase.h"
#include "TMatrixD.h"
#include "TRobustEstimator.h"
#include "AliMathBase.h"
#include "TStopwatch.h"

ClassImp(AliNDLocalRegression)

TObjArray *AliNDLocalRegression::fgVisualCorrection=0;
// instance of correction for visualization
Int_t AliNDLocalRegression::fgVerboseLevel=1000;

AliNDLocalRegression::AliNDLocalRegression():
  TNamed(),           
  fHistPoints(0),            // ND histogram defining regression granularity
  fRobustFractionLTS(0),           //   fraction of data used for the robust mean and robust rms estimator (LTS https://en.wikipedia.org/wiki/Least_trimmed_squares)
  fRobustRMSLTSCut(0),           //  cut on the robust RMS  |value-localmean|<fRobustRMSLTSCut*localRMS
  fCutType(0),                    //  type of the cut 0- no cut 1-cut localmean=median, 2-cut localmen=rosbut mean 
  fInputTree(0),             // input tree - object is not owner
  fStreamer(0),              // optional streamer 
  fFormulaVal(0),            // value:err  definition formula
  fSelection(0),             // point selector formula
  fFormulaVar(0),            //: separated variable   definition formula
  fKernelWidthFormula(0),    //: separated  - kernel width for the regression
  fPolDimensionFormula(0),   //: separated  - polynom for the regression
  fNParameters(0),           // number of local paramters to fit
  fLocalFitParam(0),         // local fit parameters 
  fLocalFitQuality(0),         // local fit quality
  fLocalFitCovar(0),          // local fit covariance matrix
  fLocalRobustStat(0),         // local robust statistic
  fBinIndex(0),                  //[fNParameters] working arrays current bin index
  fBinCenter(0),                 //[fNParameters] working current local variables - bin center
  fBinDelta(0),                  //[fNParameters] working current local variables - bin delta
  fBinWidth(0),                  //[fNParameters] working current local variables - bin delta
  fUseBinNorm(kFALSE)            //  switch make polynom  in units of bins (kTRUE)  or  in natural units (kFALSE)
{
  if (!fgVisualCorrection) fgVisualCorrection= new TObjArray;
}

AliNDLocalRegression::AliNDLocalRegression(const char* name, const char* title):
  TNamed(name,title),				
 fHistPoints(0),            // ND histogram defining regression granularity
  fRobustFractionLTS(0),           //   fraction of data used for the robust mean and robust rms estimator (LTS https://en.wikipedia.org/wiki/Least_trimmed_squares)
  fRobustRMSLTSCut(0),           //  cut on the robust RMS  |value-localmean|<fRobustRMSLTSCut*localRMS
  fCutType(0),                    //  type of the cut 0- no cut 1-cut localmean=median, 2-cut localmen=rosbut mean 
  fInputTree(0),             // input tree - object is not owner
  fStreamer(0),              // optional streamer 
  fFormulaVal(0),            // value:err  definition formula
  fSelection(0),             // point selector formula
  fFormulaVar(0),            //: separated variable   definition formula
  fKernelWidthFormula(0),    //: separated  - kernel width for the regression
  fPolDimensionFormula(0),   //: separated  - polynom for the regression
  fNParameters(0),           // number of local paramters to fit
  fLocalFitParam(0),         // local fit parameters 
  fLocalFitQuality(0),         // local fit quality
  fLocalFitCovar(0),          // local fit covariance matrix
  fLocalRobustStat(0),         // local robust statistic
  fBinIndex(0),                  //[fNParameters] working arrays current bin index
  fBinCenter(0),                 //[fNParameters] working current local variables - bin center
  fBinDelta(0),                  //[fNParameters] working current local variables - bin delta
  fBinWidth(0),                  //[fNParameters] working current local variables - bin delta
  fUseBinNorm(kFALSE)            //  switch make polynom  in units of bins (kTRUE)  or  in natural units (kFALSE)
{
}

AliNDLocalRegression::~AliNDLocalRegression(){
  //
  // destructor
  //
  if (fHistPoints) delete fHistPoints; 
  if (fStreamer)   delete fStreamer;
  //  fInputTree(0),             //! input tree - object is not owner
  delete fFormulaVal;            // value:err  definition formula
  delete fSelection;             // point selector formula
  delete fFormulaVar;            //: separated variable   definition formula
  delete fKernelWidthFormula;    //: separated  - kernel width for the regression
  delete fPolDimensionFormula;   //: separated  - polynom for the regression
  //
  if (fLocalFitParam) fLocalFitParam->Delete();
  if (fLocalFitQuality) fLocalFitQuality->Delete();
  if (fLocalFitCovar) fLocalFitCovar->Delete();
  delete fLocalFitParam;         // local fit parameters 
  delete fLocalFitQuality;       // local fit quality
  delete fLocalFitCovar;         // local fit covariance matrix
  //
  delete fLocalRobustStat;
  //
  delete[] fBinIndex;
  delete[] fBinCenter;
  delete[] fBinDelta;
}

void AliNDLocalRegression::SetHistogram(THn* histo ){
  //
  // Setup the local regression ayout according THn hitogram binning
  //
  if (fHistPoints!=0){
    AliError("Hostogram initialized\n");
    return ;
  }
  fHistPoints=histo;
  fLocalFitParam = new TObjArray(fHistPoints->GetNbins());
  fLocalFitParam->SetOwner(kTRUE);
  fLocalFitQuality = new TObjArray(fHistPoints->GetNbins());
  fLocalFitQuality->SetOwner(kTRUE);
  fLocalFitCovar = new TObjArray(fHistPoints->GetNbins());
  fLocalFitCovar->SetOwner(kTRUE);
  
  //
  // Check histogram
  //
  Int_t ndim = histo->GetNdimensions();
  Bool_t isOK=kTRUE;
  for (Int_t idim=0; idim<ndim; idim++){
    TAxis * axis = histo->GetAxis(idim);
    if (axis->GetNbins()<2) {
      AliError(TString::Format("Invalid binning nbins<2 %d",  axis->GetNbins()).Data());
    }
    if (axis->GetXmin()>=axis->GetXmax()) {
      AliError(TString::Format("Invalid range <%f,%f", axis->GetXmin(),axis->GetXmax()).Data());
    }    
  }


}
void  AliNDLocalRegression::SetCuts(Double_t nSigma, Double_t robustFraction, Int_t estimator){
  //
  //
  //
  fRobustFractionLTS=robustFraction;    //  fraction of data used for the robust mean and robust rms estimator (LTS https://en.wikipedia.org/wiki/Least_trimmed_squares)
  fRobustRMSLTSCut=nSigma;              //  cut on the robust RMS  |value-localmean|<fRobustRMSLTSCut*localRMS
  fCutType=estimator;                   //  type of the cut 0- no cut 1-cut localmean=median, 2-cut localmen=rosbut mean 

}

Bool_t    AliNDLocalRegression::CleanCovariance(){
  //
  //  Clean covariance matrix if not needed anymore
  //
  if (fLocalFitCovar) delete fLocalFitCovar;
  fLocalFitCovar=0;
};


Bool_t AliNDLocalRegression::MakeFit(TTree * tree , const char* formulaVal, const char * formulaVar, const char*selection, const char * formulaKernel, const char * dimensionFormula, Double_t weightCut, Int_t entries, Bool_t useBinNorm){
  //
  //  Make a local fit in grid as specified by the input THn histogram
  //  Histogram has to be set before invocation of method
  //
  //  Output:
  //    array of fit parameters and covariance matrices  
  //  
  //  Input Parameters:
  //   tree        - input tree
  //   formulaVal  - : separated variable:error string
  //   formulaVar  - : separate varaible list
  //   selection   - selection (cut) for TTreeDraw
  //   kernelWidth - : separated list of width of kernel for local fitting
  //   dimenstionFormula - dummy for the moment
  //
  //Algorithm:
  //   1.) Check consistency of input data
  //
  //   2.) Cache input data from tree to the array of vector TVectorD
  //
  //   3.) Calculate robust local mean and robust local RMS in case outlier removal algorithm specified
  // 
  //   4.) Make local fit
  //
  //  const Double_t kEpsilon=1e-6;
  const Int_t kMaxDim=100;
  Int_t binRange[kMaxDim]={0};
  Double_t nchi2Cut=-2*TMath::Log(weightCut); // transform probability to nsigma cut 
  if (fHistPoints==NULL){
    AliError("ND histogram not initialized\n");
    return kFALSE;
  }
  if (tree==NULL || tree->GetEntries()==0){
    AliError("Empty tree\n");
    return kFALSE;
  }
  if (formulaVar==NULL || formulaVar==0) {
    AliError("Empty variable list\n");
    return kFALSE;
  }
  if (formulaKernel==NULL) {
    AliError("Kernel width not specified\n");
    return kFALSE;
  }
  fUseBinNorm=useBinNorm;
  //  
  fInputTree= tree;  // should be better TRef?
  fFormulaVal           = new TObjString(formulaVal);
  fFormulaVar           = new TObjString(formulaVar);
  fSelection            = new TObjString(selection);
  fKernelWidthFormula   = new TObjString(formulaKernel);
  fPolDimensionFormula  = new TObjString(dimensionFormula);
  TObjArray * arrayFormulaVar=fFormulaVar->String().Tokenize(":");
  Int_t nvarFormula = arrayFormulaVar->GetEntries();
  if (nvarFormula!=fHistPoints->GetNdimensions()){
    AliError("Histogram/points mismatch\n");
    return kFALSE;
  }
  TObjArray * arrayKernel=fKernelWidthFormula->String().Tokenize(":");
  Int_t nwidthFormula = arrayKernel->GetEntries();
  if (nvarFormula!=nwidthFormula){
    delete arrayKernel;
    delete arrayFormulaVar;
    AliError("Variable/Kernel mismath\n");
    return kFALSE;
  }
  fNParameters=nvarFormula;
  //
  // 2.) Load input data
  //
  //
  Int_t entriesVal = tree->Draw(formulaVal,selection,"goffpara",entries);
  if (entriesVal==0) {
    AliError(TString::Format("Empty point list\t%s\t%s\n",formulaVal,selection).Data());
    return kFALSE; 
  }
  if (tree->GetVal(0)==NULL || (tree->GetVal(1)==NULL)){
    AliError(TString::Format("Wrong selection\t%s\t%s\n",formulaVar,selection).Data());
    return kFALSE; 
  }
  TVectorD values(entriesVal,tree->GetVal(0));
  TVectorD errors(entriesVal,tree->GetVal(1));
  Double_t *pvalues = values.GetMatrixArray();
  Double_t *perrors = errors.GetMatrixArray();


  // 2.b) variables
  TObjArray pointArray(fNParameters);
  Int_t entriesVar = tree->Draw(formulaVar,selection,"goffpara",entries);
  if (entriesVal!=entriesVar) {
    AliError(TString::Format("Wrong selection\t%s\t%s\n",formulaVar,selection).Data());
    return kFALSE; 
  }
  for (Int_t ipar=0; ipar<fNParameters; ipar++) pointArray.AddAt(new TVectorD(entriesVar,tree->GetVal(ipar)),ipar);
  // 2.c) kernel array 
  TObjArray kernelArrayI2(fNParameters);
  Int_t entriesKernel = tree->Draw(formulaKernel,selection,"goffpara",entries);
  for (Int_t ipar=0; ipar<fNParameters; ipar++) {
    TVectorD* vdI2 = new TVectorD(entriesVar,tree->GetVal(ipar));
    for (int k=entriesVar;k--;) { // to speed up, precalculate inverse squared
      double kv = (*vdI2)[k];
      if (TMath::Abs(kv)<1e-12) AliFatalClassF("Kernel width=%f for entry %d of par:%d",kv,k,ipar);
      (*vdI2)[k] = 1./(kv*kv);
    }
    kernelArrayI2.AddAt(vdI2,ipar);
  }
  //
  Double_t *pvecVar[kMaxDim]={0};
  Double_t *pvecKernelI2[kMaxDim]={0};
  for (Int_t idim=0; idim<pointArray.GetEntries(); idim++){
    pvecVar[idim]=((TVectorD*)(pointArray.At(idim)))->GetMatrixArray();
    pvecKernelI2[idim]=((TVectorD*)(kernelArrayI2.At(idim)))->GetMatrixArray();
    binRange[idim]=fHistPoints->GetAxis(idim)->GetNbins();
  }
  //
  //
  //
  Int_t nbins = fHistPoints->GetNbins();
  fBinIndex   = new Int_t[fHistPoints->GetNdimensions()];
  fBinCenter  = new Double_t[fHistPoints->GetNdimensions()];
  fBinDelta   = new Double_t[fHistPoints->GetNdimensions()];
  fBinWidth   = new Double_t[fHistPoints->GetNdimensions()];

  //
  // 3.) 
  //
  if (fCutType>0 && fRobustRMSLTSCut>0){
    MakeRobustStatistic(values, errors,  pointArray, kernelArrayI2, weightCut, fRobustFractionLTS);
  }
  //

  // 4.) Make local fits
  //

  Double_t *binHypFit  = new Double_t[2*fHistPoints->GetNdimensions()];
  //
  TLinearFitter fitter(1+2*fNParameters,TString::Format("hyp%d",2*fNParameters).Data());
  for (Int_t ibin=0; ibin<nbins; ibin++) {
    fHistPoints->GetBinContent(ibin,fBinIndex); 
    Bool_t isUnderFlowBin=kFALSE;
    Bool_t isOverFlowBin=kFALSE;
    for (Int_t idim=0; idim<fNParameters; idim++) {      
      if (fBinIndex[idim]==0) isUnderFlowBin=kTRUE;
      if (fBinIndex[idim]>binRange[idim]) isOverFlowBin=kTRUE;      
      fBinCenter[idim]=fHistPoints->GetAxis(idim)->GetBinCenter(fBinIndex[idim]);
      fBinWidth[idim]=fHistPoints->GetAxis(idim)->GetBinWidth(fBinIndex[idim]);
    }
    if (isUnderFlowBin || isOverFlowBin) continue;
    fitter.ClearPoints();
    // add fit points    
    for (Int_t ipoint=0; ipoint<entriesVal; ipoint++){
      Double_t sumChi2=0;
      if (fCutType>0 && fRobustRMSLTSCut>0){
	Double_t localRMS=(*fLocalRobustStat)(ibin,2);
	Double_t localMean=(*fLocalRobustStat)(ibin,1);
	Double_t localMedian=(*fLocalRobustStat)(ibin,0);
	if (fCutType==1){
	  if (TMath::Abs(pvalues[ipoint]-localMedian)>fRobustRMSLTSCut*localRMS) continue;
	}
	if (fCutType==2){
	  if (TMath::Abs(pvalues[ipoint]-localMean)>fRobustRMSLTSCut*localRMS) continue;
	}
      }
      for (Int_t idim=0; idim<fNParameters; idim++){
	//TVectorD &vecVar=*((TVectorD*)(pointArray.UncheckedAt(idim)));
	//TVectorD &vecKernel=*((TVectorD*)(kernelArray.UncheckedAt(idim)));
	fBinDelta[idim]=pvecVar[idim][ipoint]-fBinCenter[idim];       	
	sumChi2+= (fBinDelta[idim]*fBinDelta[idim]) * pvecKernelI2[idim][ipoint];
	if (sumChi2>nchi2Cut) break;//continue;
	if (fUseBinNorm){
	  binHypFit[2*idim]=fBinDelta[idim]/fBinWidth[idim];
	  binHypFit[2*idim+1]=binHypFit[2*idim]*binHypFit[2*idim];
	}else{
	  binHypFit[2*idim]=fBinDelta[idim];
	  binHypFit[2*idim+1]=fBinDelta[idim]*fBinDelta[idim];
	}
      }      
      if (sumChi2>nchi2Cut) continue;
      //      Double_t weight=TMath::Exp(-sumChi2*0.5);
      //      fitter.AddPoint(binHypFit,pvalues[ipoint], perrors[ipoint]/weight);
      Double_t weightI=TMath::Exp(sumChi2*0.5);
      fitter.AddPoint(binHypFit,pvalues[ipoint], perrors[ipoint]*weightI);
    }
    TVectorD * fitParam=new TVectorD(fNParameters*2+1);
    TVectorD * fitQuality=new TVectorD(3);
    TMatrixD * fitCovar=new TMatrixD(fNParameters*2+1,fNParameters*2+1);
    Double_t normRMS=0;
    Int_t nBinPoints=fitter.GetNpoints();
    Bool_t fitOK=kFALSE;
    (*fitQuality)[0]=0;
    (*fitQuality)[1]=0;
    (*fitQuality)[2]=0;

    if (fitter.GetNpoints()>fNParameters*2+2){
      fitOK = (fitter.Eval()==0);
      if (fitOK){
	normRMS=fitter.GetChisquare()/(fitter.GetNpoints()-fitter.GetNumberFreeParameters());
	fitter.GetParameters(*fitParam);
	fitter.GetCovarianceMatrix(*fitCovar);
	(*fitQuality)[0]=nBinPoints;
	(*fitQuality)[1]=normRMS;    
	(*fitQuality)[2]=ibin;    	
	fLocalFitParam->AddAt(fitParam,ibin);
	fLocalFitQuality->AddAt(fitQuality,ibin);
	fLocalFitCovar->AddAt(fitCovar,ibin);      
      }
    }
    if (fStreamer){
      TVectorD pfBinCenter(fNParameters, fBinCenter);
      Double_t median=0,mean=0,rms=0;
      if (fLocalRobustStat){
	median=(*fLocalRobustStat)(ibin,0);
	mean=(*fLocalRobustStat)(ibin,1);
	rms=(*fLocalRobustStat)(ibin,2);
      }
      (*fStreamer)<<"localFit"<<
	"ibin="<<ibin<<                // bin index
	"fitOK="<<fitOK<< 
	"localMedian="<<median<<
	"localMean="<<mean<<
	"localRMS="<<rms<<
	"nBinPoints="<<nBinPoints<<    // center of the bin
	"binCenter.="<<&pfBinCenter<<  // 
	"normRMS="<<normRMS<<          
	"fitParam.="<<fitParam<<
	"fitCovar.="<<fitCovar<<
	"fitOK="<<fitOK<<
	"\n";
    }
    if (!fitOK) { // avoid memory leak for failed fits
      delete fitParam;
      delete fitQuality;
      delete fitCovar;
    }
  }

  return kTRUE;
}


Bool_t  AliNDLocalRegression::MakeRobustStatistic(TVectorD &values,TVectorD &errors,  TObjArray &pointArray,  TObjArray &kernelArrayI2, Double_t weightCut, Double_t robustFraction){
  //
  // Calculate robust statistic information
  // use raw array to make faster calcualtion
  const Int_t kMaxDim=100;
  Double_t *pvalues=values.GetMatrixArray();
  Double_t *pvecVar[kMaxDim]={0};
  Double_t *pvecKernelI2[kMaxDim]={0};
  for (Int_t idim=0; idim<pointArray.GetEntries(); idim++){
    pvecVar[idim]=((TVectorD*)(pointArray.At(idim)))->GetMatrixArray();
    pvecKernelI2[idim]=((TVectorD*)(kernelArrayI2.At(idim)))->GetMatrixArray();
    
  }


  Double_t nchi2Cut=-2*TMath::Log(weightCut); // transform probability to nsigma cut 
  if (robustFraction>1) robustFraction=1;

  Int_t nbins = fHistPoints->GetNbins();    // 
  Int_t npoints= values.GetNrows();         // number of points for fit
  if (fLocalRobustStat){
    delete fLocalRobustStat;
  }
  fLocalRobustStat=new TMatrixD(nbins,3);

  TVectorD valueLocal(npoints);
  for (Int_t ibin=0; ibin<nbins; ibin++){
    fHistPoints->GetBinContent(ibin,fBinIndex); // 
    for (Int_t idim=0; idim<fNParameters; idim++){
      fBinCenter[idim]=fHistPoints->GetAxis(idim)->GetBinCenter(fBinIndex[idim]);
      fBinWidth[idim]=fHistPoints->GetAxis(idim)->GetBinWidth(fBinIndex[idim]);
    }
    Int_t indexLocal=0;
    for (Int_t ipoint=0; ipoint<npoints; ipoint++){
      Double_t sumChi2=0;
      for (Int_t idim=0; idim<fNParameters; idim++){
	//TVectorD &vecVar=*((TVectorD*)(pointArray.UncheckedAt(idim)));
	//TVectorD &vecKernel=*((TVectorD*)(kernelArray.UncheckedAt(idim)));
	fBinDelta[idim]=pvecVar[idim][ipoint]-fBinCenter[idim];       	
	sumChi2+= (fBinDelta[idim]*fBinDelta[idim]) * pvecKernelI2[idim][ipoint];
	if (sumChi2>nchi2Cut) break; //continue;
      }      
      if (sumChi2>nchi2Cut) continue;
      valueLocal[indexLocal]=pvalues[ipoint];
      indexLocal++;
    }
    Double_t median=0,meanX=0, rmsX=0;
    if (indexLocal*robustFraction-1>3){
      median=TMath::Median(indexLocal,valueLocal.GetMatrixArray());
      AliMathBase::EvaluateUni(indexLocal,valueLocal.GetMatrixArray(), meanX,rmsX, indexLocal*robustFraction-1);
    }
    (*fLocalRobustStat)(ibin,0)=median;
    (*fLocalRobustStat)(ibin,1)=meanX;
    (*fLocalRobustStat)(ibin,2)=rmsX;
  }
}



Double_t AliNDLocalRegression::Eval(Double_t *point ){
  //
  //
  // 
  const Double_t almost0=0.00000001;
  // backward compatibility
  if(!fBinWidth){
    fBinWidth   = new Double_t[fHistPoints->GetNdimensions()];
  }

  for (Int_t iDim=0; iDim<fNParameters; iDim++){
    if (point[iDim]<= fHistPoints->GetAxis(iDim)->GetXmin())   point[iDim]=fHistPoints->GetAxis(iDim)->GetXmin()+almost0*fHistPoints->GetAxis(iDim)->GetBinWidth(0);
    if (point[iDim]>= fHistPoints->GetAxis(iDim)->GetXmax())   point[iDim]=fHistPoints->GetAxis(iDim)->GetXmax()-almost0*fHistPoints->GetAxis(iDim)->GetBinWidth(0);
  }

  Int_t ibin = fHistPoints->GetBin(point);
  Bool_t rangeOK=kTRUE;
  if (ibin>=fLocalFitParam->GetEntriesFast() ){
    rangeOK=kFALSE;
  }else{
    if (fLocalFitParam->UncheckedAt(ibin)==NULL) {
      rangeOK=kFALSE; 
    }
  }
  if (!rangeOK) return 0;

  fHistPoints->GetBinContent(ibin,fBinIndex); 
  for (Int_t idim=0; idim<fNParameters; idim++){
    fBinCenter[idim]=fHistPoints->GetAxis(idim)->GetBinCenter(fBinIndex[idim]);
    fBinWidth[idim]=fHistPoints->GetAxis(idim)->GetBinWidth(fBinIndex[idim]);
  } 
  TVectorD &vecParam = *((TVectorD*)fLocalFitParam->At(ibin));
  Double_t value=vecParam[0];
  if (!rangeOK) return value;
  if (fUseBinNorm) {
    for (Int_t ipar=0; ipar<fNParameters; ipar++){
      Double_t delta=(point[ipar]-fBinCenter[ipar])/fBinWidth[ipar];
      value+=(vecParam[1+2*ipar]+vecParam[1+2*ipar+1]*delta)*delta;
    }
  } else {
    for (Int_t ipar=0; ipar<fNParameters; ipar++){
      Double_t delta=(point[ipar]-fBinCenter[ipar]);
      value+=(vecParam[1+2*ipar]+vecParam[1+2*ipar+1]*delta)*delta;
    }
  }
  return value;
}

Double_t AliNDLocalRegression::EvalError(Double_t *point ){
  //
  //
  // 
  if (fLocalFitCovar==NULL) {
    ::Error("AliNDLocalRegression::EvalError", "Covariance matrix not available");
    return 0 ;
  }
  for (Int_t iDim=0; iDim<fNParameters; iDim++){
    if (point[iDim]< fHistPoints->GetAxis(iDim)->GetXmin())   point[iDim]=fHistPoints->GetAxis(iDim)->GetXmin();
    if (point[iDim]> fHistPoints->GetAxis(iDim)->GetXmax())   point[iDim]=fHistPoints->GetAxis(iDim)->GetXmax();
  }

  Int_t ibin = fHistPoints->GetBin(point); 
  if (fLocalFitParam->At(ibin)==NULL) return 0;
  fHistPoints->GetBinContent(ibin,fBinIndex); 
  for (Int_t idim=0; idim<fNParameters; idim++){
    fBinCenter[idim]=fHistPoints->GetAxis(idim)->GetBinCenter(fBinIndex[idim]);
  } 
  TMatrixD &vecCovar = *((TMatrixD*)fLocalFitCovar->At(ibin));
  //TVectorD &vecQuality = *((TVectorD*)fLocalFitQuality->At(ibin));
  Double_t value=TMath::Sqrt(vecCovar(0,0));  // fill covariance to be used 
  return value;
}

Bool_t AliNDLocalRegression::Derivative(Double_t *point, Double_t *d)
{
  // fill d by partial derivatives
  //
  const Double_t almost0=0.00000001;
  for (Int_t iDim=0; iDim<fNParameters; iDim++){
    const TAxis* ax = fHistPoints->GetAxis(iDim);
    if      (point[iDim]<=ax->GetXmin()) point[iDim] = ax->GetXmin()+almost0*ax->GetBinWidth(0);
    else if (point[iDim]>=ax->GetXmax()) point[iDim] = ax->GetXmax()-almost0*ax->GetBinWidth(0);
  }

  Int_t ibin = fHistPoints->GetBin(point);
  if (ibin>=fLocalFitParam->GetEntriesFast() ||
      !fLocalFitParam->UncheckedAt(ibin) ) return kFALSE;
  //
  fHistPoints->GetBinContent(ibin,fBinIndex); 
  for (Int_t idim=0; idim<fNParameters; idim++){
    const TAxis* ax = fHistPoints->GetAxis(idim);
    fBinCenter[idim] = ax->GetBinCenter(fBinIndex[idim]);
    fBinWidth[idim] = ax->GetBinWidth(fBinIndex[idim]);
  } 
  TVectorD &vecParam = *((TVectorD*)fLocalFitParam->At(ibin));

  if (fUseBinNorm){
    for (Int_t ipar=0; ipar<fNParameters; ipar++){
      Double_t delta=(point[ipar]-fBinCenter[ipar])/fBinWidth[ipar];
      d[ipar] = vecParam[1+2*ipar] + 2.0*vecParam[1+2*ipar+1]*delta;
    }
  }else{
    for (Int_t ipar=0; ipar<fNParameters; ipar++){
      Double_t delta=(point[ipar]-fBinCenter[ipar]);
      d[ipar] = vecParam[1+2*ipar] + 2.0*vecParam[1+2*ipar+1]*delta;
    }
  }
  return kTRUE;
}

Bool_t AliNDLocalRegression::EvalAndDerivative(Double_t *point, Double_t &val, Double_t *d)
{
  // fill d by partial derivatives and calculate value in one go
  //
  const Double_t almost0=0.00000001;
  for (Int_t iDim=0; iDim<fNParameters; iDim++){
    const TAxis* ax = fHistPoints->GetAxis(iDim);
    if      (point[iDim]<=ax->GetXmin()) point[iDim] = ax->GetXmin()+almost0*ax->GetBinWidth(0);
    else if (point[iDim]>=ax->GetXmax()) point[iDim] = ax->GetXmax()-almost0*ax->GetBinWidth(0);
  }
  val = 0;
  Int_t ibin = fHistPoints->GetBin(point);
  if (ibin>=fLocalFitParam->GetEntriesFast() ||
      !fLocalFitParam->UncheckedAt(ibin) ) return kFALSE;
  //
  fHistPoints->GetBinContent(ibin,fBinIndex); 
  for (Int_t idim=0; idim<fNParameters; idim++){
    const TAxis* ax = fHistPoints->GetAxis(idim);
    fBinCenter[idim] = ax->GetBinCenter(fBinIndex[idim]);
    fBinWidth[idim] = ax->GetBinWidth(fBinIndex[idim]);
  } 
  TVectorD &vecParam = *((TVectorD*)fLocalFitParam->At(ibin));
  
  val = vecParam[0];
  if (fUseBinNorm){
    for (Int_t ipar=0; ipar<fNParameters; ipar++){
      Double_t delta=(point[ipar]-fBinCenter[ipar])/fBinWidth[ipar];
      d[ipar] = vecParam[1+2*ipar] + 2.0*vecParam[1+2*ipar+1]*delta;
      val += (vecParam[1+2*ipar]+vecParam[1+2*ipar+1]*delta)*delta;
    }
  } else {
    for (Int_t ipar=0; ipar<fNParameters; ipar++){
      Double_t delta=(point[ipar]-fBinCenter[ipar]);
      d[ipar] = vecParam[1+2*ipar] + 2.0*vecParam[1+2*ipar+1]*delta;
      val += (vecParam[1+2*ipar]+vecParam[1+2*ipar+1]*delta)*delta;
    }
  }
  return kTRUE;
}


Int_t  AliNDLocalRegression::GetVisualCorrectionIndex(const char *corName){
  //
  return TMath::Hash(corName)%1000000;
}

    
void AliNDLocalRegression::AddVisualCorrection(AliNDLocalRegression* corr, Int_t position){
  /// make correction available for visualization using
  /// TFormula, TFX and TTree::Draw
  /// important in order to check corrections and also compute dervied variables
  /// e.g correction partial derivatives
  ///
  /// NOTE - class is not owner of correction
  if (position==0) {
    position=GetVisualCorrectionIndex(corr->GetName());
  }
  
  if (!fgVisualCorrection) fgVisualCorrection=new TObjArray(1000000);
  if (position>=fgVisualCorrection->GetEntriesFast())
    fgVisualCorrection->Expand((position+10)*2);
  if (fgVisualCorrection->At(position)!=NULL){
    ::Error("AliNDLocalRegression::AddVisualCorrection","Correction %d already defined Old: %s New: %s",position,fgVisualCorrection->At(position)->GetName(), corr->GetName());
  }
  fgVisualCorrection->AddAt(corr, position);
}

AliNDLocalRegression* AliNDLocalRegression::GetVisualCorrection(Int_t position) {
  /// Get visula correction registered at index=position  
  return fgVisualCorrection? (AliNDLocalRegression*)fgVisualCorrection->At(position):0;
}

Double_t AliNDLocalRegression::GetCorrND(Double_t index, Double_t par0){
  //
  //
  AliNDLocalRegression *corr = (AliNDLocalRegression*)fgVisualCorrection->At(index);
  if (!corr) return 0;
  return corr->Eval(&par0);
}

Double_t AliNDLocalRegression::GetCorrNDError(Double_t index, Double_t par0){
  //
  //
  AliNDLocalRegression *corr = (AliNDLocalRegression*)fgVisualCorrection->At(index);
  if (!corr) return 0;
  return corr->EvalError(&par0);
}

Double_t AliNDLocalRegression::GetCorrND(Double_t index, Double_t par0, Double_t par1){
  //
  //
  AliNDLocalRegression *corr = (AliNDLocalRegression*)fgVisualCorrection->At(index);
  if (!corr) return 0;
  Double_t par[2]={par0,par1};
  return corr->Eval(par);
}
Double_t AliNDLocalRegression::GetCorrNDError(Double_t index, Double_t par0, Double_t par1){
  //
  //
  AliNDLocalRegression *corr = (AliNDLocalRegression*)fgVisualCorrection->At(index);
  if (!corr) return 0;
  Double_t par[2]={par0,par1};
  return corr->EvalError(par);
}

Double_t AliNDLocalRegression::GetCorrND(Double_t index, Double_t par0, Double_t par1, Double_t par2){
  //
  //
  AliNDLocalRegression *corr = (AliNDLocalRegression*)fgVisualCorrection->At(index);
  if (!corr) return 0;
  Double_t par[3]={par0,par1,par2};
  return corr->Eval(par);
}

Double_t AliNDLocalRegression::GetCorrNDError(Double_t index, Double_t par0, Double_t par1, Double_t par2){
  //
  //
  AliNDLocalRegression *corr = (AliNDLocalRegression*)fgVisualCorrection->At(index);
  if (!corr) return 0;
  Double_t par[3]={par0,par1,par2};
  return corr->EvalError(par);
}



Double_t AliNDLocalRegression::GetCorrND(Double_t index, Double_t par0, Double_t par1, Double_t par2, Double_t par3){
  //
  //
  AliNDLocalRegression *corr = (AliNDLocalRegression*)fgVisualCorrection->At(index);
  if (!corr) return 0;
  Double_t par[4]={par0,par1,par2,par3};
  return corr->Eval(par);
}

Double_t AliNDLocalRegression::GetCorrNDError(Double_t index, Double_t par0, Double_t par1, Double_t par2, Double_t par3){
  //
  //
  AliNDLocalRegression *corr = (AliNDLocalRegression*)fgVisualCorrection->At(index);
  if (!corr) return 0;
  Double_t par[4]={par0,par1,par2,par3};
  return corr->EvalError(par);
}




Bool_t AliNDLocalRegression::AddWeekConstrainsAtBoundaries(Int_t nDims, Int_t *indexes, Double_t *relWeight, TTreeSRedirector* pcstream, Bool_t useCommon){
   //
  // Adding week constrain AtBoundaries
  //
  //  Technique similar to "Kalman update" of measurement used at boundaries - https://en.wikipedia.org/wiki/Kalman_filter
  // 
  // 1.) Make backup of original parameters
  // 2.) Book Kalman matrices
  // 3.) Loop over all measurements bins and update mesurements -adding boundary measurements as additional measurement
  //     relWeight vector specify relative weight of such measurement  (err_i=sigma_i*refWeight_i) - not yet implemented
  // 4.) replace original parameters with constrained parameters
  //     procedure can be repeated 
  /*
    Input parameters example:
    nDims=2;
    Int_t indexes[2]={0,1};
    Double_t relWeight0[6]={1,1,1,1,1,1};
    Double_t relWeight1[6]={1,1,10,1,1,10};
    pcstream=new TTreeSRedirector("constrainStream.root","recreate");
    
    AliNDLocalRegression * regression0 = ( AliNDLocalRegression *)AliNDLocalRegression::GetVisualCorrections()->FindObject("pfitNDGaus0");
    AliNDLocalRegression * regression1 = ( AliNDLocalRegression *)AliNDLocalRegression::GetVisualCorrections()->FindObject("pfitNDGaus1");

    regressionUpdate0 = (AliNDLocalRegression *)regression0->Clone();
    regressionUpdate1 = (AliNDLocalRegression *)regression1->Clone();
    AddWeekConstrainsAtBoundaries( regressionUpdate0, nDims, indexes,relWeight0, pcstream);
    AddWeekConstrainsAtBoundaries( regressionUpdate0, nDims, indexes,relWeight0, pcstream);
    AddWeekConstrainsAtBoundaries( regressionUpdate0, nDims, indexes,relWeight0, pcstream);
    AddWeekConstrainsAtBoundaries( regressionUpdate0, nDims, indexes,relWeight0, pcstream);
    AddWeekConstrainsAtBoundaries( regressionUpdate1, nDims, indexes,relWeight1, pcstream);
    AddWeekConstrainsAtBoundaries( regressionUpdate1, nDims, indexes,relWeight1, pcstream);
    AddWeekConstrainsAtBoundaries( regressionUpdate1, nDims, indexes,relWeight1, pcstream);
    AddWeekConstrainsAtBoundaries( regressionUpdate1, nDims, indexes,relWeight1, pcstream);

    regressionUpdate0->SetName("pfitNDGaus0_Updated");
    regressionUpdate1->SetName("pfitNDGaus1_Updated");
    AliNDLocalRegression::AddVisualCorrection(regressionUpdate0);
    AliNDLocalRegression::AddVisualCorrection(regressionUpdate1);
    treeIn->SetAlias( regressionUpdate0->GetName(), TString::Format("AliNDLocalRegression::GetCorrND(%d,xyz0,xyz1+0)", regressionUpdate0->GetVisualCorrectionIndex()).Data());
     treeIn->SetAlias( regressionUpdate1->GetName(), TString::Format("AliNDLocalRegression::GetCorrND(%d,xyz0,xyz1+0)", regressionUpdate1->GetVisualCorrectionIndex()).Data());
    delete pcstream;


    TFile *f = TFile::Open("constrainStream.root")
   */
  const Double_t kScale=0.5;
  const Double_t singularity_tolerance = 1e-200;
  if (fLocalFitCovar==NULL) {
    ::Error("AliNDLocalRegression::EvalError", "Covariance matrix not available");
    return 0 ;
  }
  //
  // 1.)  Make backup of original parameters
  //
  const TObjArray *vecParamOrig    = fLocalFitParam;
  const TObjArray *vecCovarOrig    = fLocalFitCovar;
  TObjArray *vecParamUpdated = new TObjArray(fHistPoints->GetNbins());
  TObjArray *vecCovarUpdated = new TObjArray(fHistPoints->GetNbins());
  // 
  // 2.) Book local varaibles and Kalman matrices
  //  
  Int_t nParams= 1+2*fNParameters;
  Int_t nMeas= nDims*6; // update each dimension specified 2 ends 2 measurements (value and first derivative)
  
  TMatrixD matWeight(nParams,nParams);       // weight matrix for side param
  TMatrixD matCovarSide(nParams,nParams);    // reweighted covariance matrix for side parameters

  TMatrixD vecXk(nParams,1);           // X vector - parameter of the local fit at bin
  TMatrixD covXk(nParams,nParams);     // X covariance 
  TMatrixD matHk(nMeas,nParams);       // vector to mesurement (values at boundary of bin)
  TMatrixD measR(nMeas,nMeas);         // measurement error at boundary as provided by bin in local neigberhood 
  TMatrixD vecZk(nMeas,1);             // measurement at boundary 
  //
  TMatrixD measRBin(nMeas,nMeas);              // measurement error bin
  TMatrixD vecZkBin(nMeas,1);                  // measurement bin
  TMatrixD matrixTransformBin(nMeas, nParams);  // vector to measurement to calculate error matrix current bin
  //
  TMatrixD vecZkSide(3,1);                // measurement side
  TMatrixD matrixTransformSide(3,nParams);// vector to measurement to calculate error matrix side bin

  //
  TMatrixD vecYk(nMeas,1);          // Innovation or measurement residual
  TMatrixD matHkT(nParams,nMeas);
  TMatrixD matSk(nMeas,nMeas);    // Innovation (or residual) covariance
  TMatrixD matKk(nParams,nMeas);    // Optimal Kalman gain
  TMatrixD mat1(nParams,nParams);     // update covariance matrix
  TMatrixD covXk2(nParams,nParams);   // 
  TMatrixD covOut(nParams,nParams);   //
  mat1.UnitMatrix();
  //
  // 3.) Loop over all measurements bins and update mesurements -adding boundary measurements as additional measurement
  //     relWeight vector specify relative weight of such measurement  (err_i=sigma_i*refWeight_i
  const THn* his = GetHistogram();
  Int_t binIndex[999]={0};
  Int_t binIndexSide[999]={0};
  Int_t nbinsAxis[999]={0};
  Double_t binCenter[999]={0};
  Double_t binWidth[999]={0};

  if (relWeight!=NULL) for (Int_t iParam=0; iParam<nParams; iParam++){
    Int_t index=0;
    if (iParam<3)  index=iParam;
    if (iParam>=3) {
      Int_t dim=(iParam-3)/2;
      Int_t deriv=1+(iParam-3)%2;
      index=3*dim+deriv;
    }
    matWeight(iParam,iParam)=relWeight[index];
  }


  for (Int_t iDim=0; iDim<nDims; iDim++){nbinsAxis[iDim]=his->GetAxis(iDim)->GetNbins();}  
  //  Int_t nBins=fHistPoints->GetNbins();
  Int_t nBins=fLocalFitParam->GetSize();
  for (Int_t iBin=0; iBin<nBins; iBin++){   // loop over bins
    if (iBin%fgVerboseLevel==0) printf("%d\n",iBin);
    //
    his->GetBinContent(iBin,binIndex);
    for (Int_t iDim=0; iDim<nDims; iDim++) { // fill common info for bin of interest
      binCenter[iDim]= his->GetAxis(iDim)->GetBinCenter(binIndex[iDim]);
      binWidth[iDim] = his->GetAxis(iDim)->GetBinWidth(binIndex[iDim]);
    }
    if (fLocalFitParam->UncheckedAt(iBin)==NULL) continue;
    Double_t *vecParam0 = ((TVectorD*)(fLocalFitParam->UncheckedAt(iBin)))->GetMatrixArray();
    TMatrixD   matParam0(nParams,1, vecParam0);
    TMatrixD & matCovar0=*(((TMatrixD*)(fLocalFitCovar->UncheckedAt(iBin))));
    measR.Zero();
    vecZk.Zero();
    measRBin.Zero();
    vecZkBin.Zero();    
    matrixTransformBin.Zero();
    covXk=matCovar0;
    vecXk=matParam0;
    //
    //  neiborhood loop
    Int_t constCounter=0;
    Int_t constCounterDim[100]={0};
    for (Int_t iDim=0; iDim<nDims; iDim++){         // loop in n dim
      constCounterDim[iDim]=0;  // number of constraints per dimension
      for (Int_t iSide=-1; iSide<=1; iSide+=2){     // left right loop
	for (Int_t jDim=0; jDim<nDims; jDim++) binIndexSide[jDim]= binIndex[jDim];
	vecZkSide.Zero();
	matrixTransformSide.Zero();
	//
	binIndexSide[iDim]+=iSide;      
	if (binIndexSide[iDim]<0) binIndexSide[iDim]=0;
	if (binIndexSide[iDim]>his->GetAxis(iDim)->GetNbins())  binIndexSide[iDim]=his->GetAxis(iDim)->GetNbins();
	Bool_t isConst=binIndexSide[iDim]>0 &&binIndexSide[iDim]<=his->GetAxis(iDim)->GetNbins() && (fLocalFitParam)->UncheckedAt(his->GetBin(binIndexSide))!=NULL; 
	Int_t binSide=his->GetBin(binIndexSide);
	if (binSide>=nBins ) binIndexSide[iDim]=binIndex[iDim];
	if ((fLocalFitParam)->UncheckedAt(his->GetBin(binIndexSide))==NULL) binIndexSide[iDim]=binIndex[iDim];
	if (isConst)  {
	  constCounter++;
	  constCounterDim[iDim]++;
	}
	Double_t localCenter=his->GetAxis(iDim)->GetBinCenter(binIndex[iDim]);
	Double_t sideCenter= his->GetAxis(iDim)->GetBinCenter(binIndexSide[iDim]);
	Double_t position=   (iSide<0) ? his->GetAxis(iDim)->GetBinLowEdge(binIndex[iDim]) :  his->GetAxis(iDim)->GetBinUpEdge(binIndex[iDim]);
	Double_t* vecParamSide  = ((TVectorD*)(fLocalFitParam)->UncheckedAt(his->GetBin(binIndexSide)))->GetMatrixArray();
	TMatrixD   matParamSide(nParams,1, vecParamSide);
	if (relWeight==NULL){
	  matCovarSide=*((TMatrixD*)(fLocalFitCovar->UncheckedAt(his->GetBin(binIndexSide))));
	}
	if (relWeight!=NULL){
	  matCovarSide=TMatrixD( matWeight,TMatrixD::kMult,*((TMatrixD*)(fLocalFitCovar->UncheckedAt(his->GetBin(binIndexSide)))));
	  matCovarSide*=matWeight;
	}
	if (!isConst) matCovarSide*=1000;
	//
	Double_t deltaLocal=(position-localCenter);
	Double_t deltaSide=(position-sideCenter);
	if (fUseBinNorm){
	   deltaLocal/=fBinWidth[iDim];
	   deltaSide/=fBinWidth[iDim];
	}
	//
	matrixTransformSide(0,0)=1;        matrixTransformSide(0,1+2*iDim)=deltaSide;      matrixTransformSide(0,1+2*iDim+1)=deltaSide*deltaSide;
	matrixTransformSide(1,1+2*iDim)=1;   matrixTransformSide(1,1+2*iDim+1)=2*deltaSide;
	matrixTransformSide(2,1+2*iDim+1)=2;
	//
	Int_t iMeas0=6*iDim+3*(iSide+1)/2;
	matrixTransformBin(iMeas0+0,0)=1;        matrixTransformBin(iMeas0+0,1+2*iDim)=deltaLocal;      matrixTransformBin(iMeas0+0,1+2*iDim+1)=deltaSide*deltaLocal;
	matrixTransformBin(iMeas0+1,1+2*iDim)=1;   matrixTransformBin(iMeas0+1,1+2*iDim+1)=2*deltaLocal;
	matrixTransformBin(iMeas0+2,1+2*iDim+1)=2;
	//
	for (Int_t iconst=0; iconst<3; iconst++){
	  Int_t iMeas=iMeas0+iconst;
	  Double_t localMeasurement=0;
	  Double_t sideMeasurement=0;
	  if (iconst==0){ // measurement - derivative 0
	    localMeasurement=vecParam0[0]+deltaLocal*(vecParam0[1+2*iDim]+vecParam0[2+2*iDim]*deltaLocal);
	    sideMeasurement=vecParamSide[0]+deltaSide*(vecParamSide[1+2*iDim]+vecParamSide[2+2*iDim]*deltaSide);
	  }
	  if (iconst==1){ // measurement -derivative 1
	    localMeasurement=(vecParam0[1+2*iDim]+2*vecParam0[2+2*iDim]*deltaLocal);
	    sideMeasurement=(vecParamSide[1+2*iDim]+2*vecParamSide[2+2*iDim]*deltaSide);
	  }
	  if (iconst==2){
	    localMeasurement=2*vecParam0[2+2*iDim];
	    sideMeasurement=2*vecParamSide[2+2*iDim];
	  }
	  vecZkSide(iconst,0)=sideMeasurement;
	  vecZk(iMeas,0)=sideMeasurement;
	  if (!isConst) vecZk(iMeas,0)=localMeasurement;
	  vecZkBin(iMeas,0)=localMeasurement;
	}
	TMatrixD measRSide0(matrixTransformSide,TMatrixD::kMult,matCovarSide);   //     (iconst,iconst)  = (iconst,nParam)*(nParams,nParams)*(nParams,iconst
	TMatrixD matrixTransformSideT(TMatrixD::kTransposed ,matrixTransformSide);
	TMatrixD measRSide(measRSide0,TMatrixD::kMult,matrixTransformSideT);
	// update measutement Covariance matrix for given side
	for (Int_t iconst=0; iconst<3; iconst++)
	  for (Int_t jconst=0; jconst<3; jconst++){
	    measR(iMeas0+iconst,iMeas0+jconst)=measRSide(iconst,jconst);
	  }
	if (pcstream){
	  TMatrixD vecZkSideCheck(matrixTransformSide,TMatrixD::kMult,matParamSide);   //     (iconst,1)       = (iConst,nParam)*(nParams,1)	
	  //
	  (*pcstream)<<"checkSide"<<  // check agreement in 1D
	    "iBin="<<iBin<<
	    "iDim="<<iDim<<
	    "iSide="<<iSide<<
	    "vecZkSide.="<<&vecZkSide<<
	    "vecZkSideCheck.="<<&vecZkSideCheck<<
	    "measRSide.="<<&measRSide<<	  
	    "vecZk.="<<&vecZk<<
	    "vecZkBin.="<<&vecZkBin<<	    
	    "\n";
	}	
      }
    }
    //
    //
    TMatrixD measRBin0(matrixTransformBin,TMatrixD::kMult,matCovar0);   //     (iconst,iconst)  = (iconst,nParam)*(nParams,nParams)*(nParams,iconst
    TMatrixD matrixTransformBinT(TMatrixD::kTransposed ,matrixTransformBin);
    TMatrixD measRBin(measRBin0,TMatrixD::kMult,matrixTransformBinT);
    //
    // make Kalman Update of state vector with side mesurement
    //
    matHk=matrixTransformBin;
    matHkT= matrixTransformBinT;
    //
    vecYk = vecZk-matHk*vecXk;                 // Innovation or measurement residual
    if (useCommon) vecYk*=0.5;                 // in case we are using middle point use only half of delta
    matSk = (matHk*(covXk*matHkT))+measR;      // Innovation (or residual) covariance
    Double_t determinant=0;
    Int_t constCounter2=0;
    for (Int_t kDim=0; kDim<nDims; kDim++) if (constCounterDim[kDim]>0) constCounter2++;
    if (constCounter2==nDims){
      determinant= matSk.Determinant();
    }
    if (TMath::Abs(determinant)<singularity_tolerance ) {
      vecParamUpdated->AddAt(new TVectorD(*((TVectorD*)(fLocalFitParam->UncheckedAt(iBin)))),iBin); 
      vecCovarUpdated->AddAt(new TMatrixD(*((TMatrixD*)(fLocalFitCovar->UncheckedAt(iBin)))),iBin);
      AliDebug(1,TString::Format("Update matrix not possible matSk.Determinant() too small, skipping bin %d",iBin).Data());
    }else{
      matSk.Invert();
      matKk = (covXk*matHkT)*matSk;              //  Optimal Kalman gain
      vecXk += matKk*vecYk;                      //  updated vector 
      covXk2 = (mat1-(matKk*matHk));
      covOut =  covXk2*covXk; 
      //
      vecParamUpdated->AddAt(new TVectorD(nParams,vecXk.GetMatrixArray()), iBin); 
      vecCovarUpdated->AddAt(new TMatrixD(covOut), iBin); 
    }
    if (pcstream){
      TMatrixD vecZkBinCheck(matrixTransformBin,TMatrixD::kMult,matParam0); 
      TVectorD vecPos(nDims,binCenter);
      TVectorD *vecXk0= (TVectorD*)(fLocalFitParam->UncheckedAt(iBin));
      TMatrixD vecYkUpdated=(vecZk-matHk*vecXk);
      //
       (*pcstream)<<"checkBin"<<       // check agreement in all sides
	 "iBin="<<iBin<<               // bin index
	 "vecPos.="<<&vecPos<<         // bin position
	 "determinant="<<determinant<<  // determinant
	 "constCounter="<<constCounter<<
	 "constCounter2="<<constCounter<< 
	 //
	 "vecXk0.="<<vecXk0<<          // original parameter vector
	 "vecXk.="<<&vecXk<<           // parameter vector at bin after update
	 "covXk.="<<&covXk<<           // covaraince matrix before update
	 "covOut.="<<&covOut<<           // covaraince matrix after update
	 "vecZk.="<<&vecZk<<           // measurement vector - values according side measurement
	 "vecZkBin.="<<&vecZkBin<<     // expected vector according parameters for bin
	 "vecZkBinCheck.="<<&vecZkBinCheck<<   // expected vector according parameters at bin centers - crosscheck tracsrormation matrix
	 "measRBin.="<<&measRBin<<     // expected error of extrapolation
	 "measR.="<<&measR<<           // error of the side measurement
	 // tmporary data
	 "vecYk.="<<&vecYk<<           // delta vector (nparams)
	 "matSk.="<<&matSk<<           // inovation covariance (nMeas,nMeas)
	 "matKk.="<<&matKk<<          // optimal Kalman gain  (nParams,nMeas) 
	 "covXk2.="<<&covXk2<<	 
	 //
	 "vecYkUpdated.="<<&vecYkUpdated<< // diff after kalman update
	 "\n";
    }
  } 
  //
  // 4.) replace original parameters with constrained parameters
  //
  fLocalFitParam= vecParamUpdated;
  fLocalFitCovar= vecCovarUpdated;  
  return 0;
}

void AliNDLocalRegression::DumpToTree(Int_t nDiv,  TTreeStream & stream){
  //
  //
  //
  const Int_t kMaxDim=100;
  Int_t nBins=fHistPoints->GetNbins();
  
  TVectorD binLowEdge(fNParameters);
  TVectorD binLocal(fNParameters);
  TVectorF binIndexF(fNParameters);
  Double_t *pbinLowEdge= binLowEdge.GetMatrixArray();
  //
  for (Int_t iBin=0; iBin<nBins; iBin++){   // loop over bins
    if (iBin%fgVerboseLevel==0) printf("%d\n",iBin);
    fHistPoints->GetBinContent(iBin,fBinIndex);
    for (Int_t iDim=0; iDim<fNParameters; iDim++) { // fill common info for bin of interest
      binIndexF[iDim]=fBinIndex[iDim];
      fBinCenter[iDim]= fHistPoints->GetAxis(iDim)->GetBinCenter(fBinIndex[iDim]);
      binLowEdge[iDim]= fHistPoints->GetAxis(iDim)->GetBinLowEdge(fBinIndex[iDim]);
      fBinDelta[iDim] = fHistPoints->GetAxis(iDim)->GetBinWidth(fBinIndex[iDim]);
    }
    
    for (Int_t iDim=0; iDim<fNParameters; iDim++){
      for (Int_t jDim=0; jDim<fNParameters; jDim++) binLocal[jDim]=fBinCenter[jDim];
      for (Int_t iDiv=0; iDiv<nDiv; iDiv++){
	binLocal[iDim]=binLowEdge[iDim]+(fBinDelta[iDim]*(iDiv-1.))/Double_t(nDiv);
	Double_t value2= Eval(binLocal.GetMatrixArray());
	binLocal[iDim]=binLowEdge[iDim]+(fBinDelta[iDim]*Double_t(iDiv))/Double_t(nDiv);
	Double_t value=Eval(binLocal.GetMatrixArray());
	Double_t derivative=nDiv*(value-value2)/fBinDelta[iDim];
	stream<<
	  "pos.="<<&binLocal<<          // position vector
	  "iBin="<<iBin<<               // bin index
	  "iDim="<<iDim<<               // scan bin index
	  "iDiv="<<iDiv<<               // division index
	  "binIndexF.="<<&binIndexF<<   // bin index
	  "value="<<value<<             // value at 
	  "derivative="<<derivative<<   // numerical derivative
	  "\n";
      }
    }
  }
}

///   Interpolate graph using local regression with kernel width
///   \param gr          - input graph  in procedure it is assumed graph is sorted 
///   \param evalTime    - X(time) to evaluat graph
///   \param kernelWidth - gaussian kernel width - see ( http://en.wikipedia.org/w/index.php?title=Kernel_smoother&oldid=627785784)
///   \param sigmaCut    - sigma cut to reduce CPU consuption only points in +-sigmaCut*kernelWidth used
///   \param evalLog     - in case kTRUE polynomial fit of the log(yi) values used (e.g for spectra measurement (exponential  or power law)  
///   \param pol - degree of polynom used, in case number of points too small (e.g on edges) smaller polynomial used 
///       - 2 deggress of freedom requested  (1 point -> pol0  , 2 points -> pol0,  3 points ->pol1)
///   \param param and covar - if specified full parameters and covaraince matrix filled


Double_t AliNDLocalRegression::EvalGraphKernel(TGraph * gr, Double_t evalTime, Double_t kernelWidth, Double_t sigmaCut, Bool_t evalLog, Int_t pol, TVectorD *param, TMatrixD *covar){
     
  const Int_t    kMinEntries=4;
  const Double_t kMaxExp=100; 
  Int_t npoints=gr->GetN();
  Int_t index0= TMath::BinarySearch(npoints, gr->GetX(), evalTime-sigmaCut*kernelWidth);
  Int_t index1= TMath::BinarySearch(npoints, gr->GetX(), evalTime+sigmaCut*kernelWidth);
  if (index1-index0 < kMinEntries) {
    Int_t index= TMath::BinarySearch(npoints, gr->GetX(), evalTime);
    index0=index-kMinEntries/2;
    index1=index+kMinEntries/2;
  }
  if (index0<0) { index1-=index0; index0=0;}
  if (index1>=npoints) {index0-=index1-npoints;   index1=npoints-1;}
  if (index0<0) index0=0;

  Int_t lpol=TMath::Min(pol, TMath::Max(index1-index0-2,0));
  TLinearFitter fitter(lpol+1, TString::Format("pol%d",lpol));
  Double_t mkernel2=1./(kernelWidth*kernelWidth);
  for (Int_t ipoint=index0; ipoint<=index1; ipoint++){
    Double_t x=gr->GetX()[ipoint]-evalTime;
    Double_t y=gr->GetY()[ipoint];
    Double_t errorW=TMath::Exp(TMath::Min(x*x*mkernel2,kMaxExp));
    if (evalLog==kFALSE){
      fitter.AddPoint(&x, gr->GetY()[ipoint], errorW);
    }else{
      if (y>0) fitter.AddPoint(&x, TMath::Log(y) , errorW);
    }
  }
  Int_t hasFailed=(fitter.GetNpoints()>kMinEntries)? fitter.Eval():1;
  if (hasFailed) return 0;
  
  if (param){
    fitter.GetParameters(*param);
  }
  if (covar) fitter.GetCovarianceMatrix(*covar);  
  return (evalLog==0) ? fitter.GetParameter(0) : TMath::Exp(fitter.GetParameter(0));
}


