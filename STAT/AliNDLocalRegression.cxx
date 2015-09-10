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

ClassImp(AliNDLocalRegression)

TObjArray *AliNDLocalRegression::fgVisualCorrection=0;
// instance of correction for visualization


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
  fBinIndex(0),                  //[fNParameters] working arrays current bin index
  fBinCenter(0),                 //[fNParameters] working current local variables - bin center
  fBinDelta(0),                  //[fNParameters] working current local variables - bin delta
  fLocalRobustStat(0)         // local robust statistic
{
  if (!fgVisualCorrection) fgVisualCorrection= new TObjArray;
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
  delete fLocalFitParam;         // local fit parameters 
  delete fLocalFitQuality;       // local fit quality
  delete fLocalFitCovar;         // local fit covariance matrix

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
  fLocalFitQuality = new TObjArray(fHistPoints->GetNbins());
  fLocalFitCovar = new TObjArray(fHistPoints->GetNbins());
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



Bool_t AliNDLocalRegression::MakeFit(TTree * tree , const char* formulaVal, const char * formulaVar, const char*selection, const char * formulaKernel, const char * dimensionFormula, Double_t weightCut, Int_t entries){
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
  // 2.b) variables
  TObjArray pointArray(fNParameters);
  Int_t entriesVar = tree->Draw(formulaVar,selection,"goffpara",entries);
  if (entriesVal!=entriesVar) {
    AliError(TString::Format("Wrong selection\t%s\t%s\n",formulaVar,selection).Data());
    return kFALSE; 
  }
  for (Int_t ipar=0; ipar<fNParameters; ipar++) pointArray.AddAt(new TVectorD(entriesVar,tree->GetVal(ipar)),ipar);
  // 2.c) kernel array 
  TObjArray kernelArray(fNParameters);
  Int_t entriesKernel = tree->Draw(formulaKernel,selection,"goffpara",entries);
  for (Int_t ipar=0; ipar<fNParameters; ipar++) kernelArray.AddAt(new TVectorD(entriesVar,tree->GetVal(ipar)),ipar);
  //
  //
  //
  Int_t nbins = fHistPoints->GetNbins();
  fBinIndex   = new Int_t[fHistPoints->GetNdimensions()];
  fBinCenter  = new Double_t[fHistPoints->GetNdimensions()];
  fBinDelta   = new Double_t[fHistPoints->GetNdimensions()];

  //
  // 3.) 
  //
  if (fCutType>0 && fRobustRMSLTSCut>0){
    MakeRobustStatistic(values, errors,  pointArray, kernelArray, weightCut, fRobustFractionLTS);
  }
  //
  // 4.) Make local fits
  //
  Double_t *binHypFit  = new Double_t[2*fHistPoints->GetNdimensions()];
  //
  TLinearFitter fitter(1+2*fNParameters,TString::Format("hyp%d",2*fNParameters).Data());
  for (Int_t ibin=0; ibin<nbins; ibin++){
    fHistPoints->GetBinContent(ibin,fBinIndex); // 
    for (Int_t idim=0; idim<fNParameters; idim++){
      fBinCenter[idim]=fHistPoints->GetAxis(idim)->GetBinCenter(fBinIndex[idim]);
    }
    fitter.ClearPoints();
    // add fit points    
    for (Int_t ipoint=0; ipoint<entriesVal; ipoint++){
      Double_t weight=1;
      if (fCutType>0 && fRobustRMSLTSCut>0){
	Double_t localRMS=(*fLocalRobustStat)(ibin,2);
	Double_t localMean=(*fLocalRobustStat)(ibin,1);
	Double_t localMedian=(*fLocalRobustStat)(ibin,0);
	if (fCutType==1){
	  if (TMath::Abs(values[ipoint]-localMedian)>fRobustRMSLTSCut*localRMS) continue;
	}
	if (fCutType==2){
	  if (TMath::Abs(values[ipoint]-localMean)>fRobustRMSLTSCut*localRMS) continue;
	}
      }
      for (Int_t idim=0; idim<fNParameters; idim++){
	TVectorD &vecVar=*((TVectorD*)(pointArray.UncheckedAt(idim)));
	TVectorD &vecKernel=*((TVectorD*)(kernelArray.UncheckedAt(idim)));
	fBinDelta[idim]=vecVar[ipoint]-fBinCenter[idim];       	
	weight*=TMath::Gaus(fBinDelta[idim],0,vecKernel[ipoint]);
	if (weight<weightCut) continue;
	binHypFit[2*idim]=fBinDelta[idim];
	binHypFit[2*idim+1]=fBinDelta[idim]*fBinDelta[idim];
      }      
      if (weight<weightCut) continue;
      fitter.AddPoint(binHypFit,values[ipoint], errors[ipoint]/weight);
    }
    TVectorD * fitParam=new TVectorD(fNParameters*2+1);
    TVectorD * fitQuality=new TVectorD(3);
    TMatrixD * fitCovar=new TMatrixD(fNParameters*2+1,fNParameters*2+1);
    Double_t normRMS=0;
    Int_t nBinPoints=fitter.GetNpoints();
    if (fitter.GetNpoints()>fNParameters*2+2){
      fitter.Eval();
      normRMS=fitter.GetChisquare()/(fitter.GetNpoints()-fitter.GetNumberFreeParameters());
      fitter.GetParameters(*fitParam);
      fitter.GetCovarianceMatrix(*fitCovar);
      (*fitQuality)[0]=nBinPoints;
      (*fitQuality)[1]=normRMS;
    }
    fLocalFitParam->AddAt(fitParam,ibin);
    fLocalFitQuality->AddAt(fitQuality,ibin);
    fLocalFitCovar->AddAt(fitCovar,ibin);
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
	"localMedian="<<median<<
	"localMean="<<mean<<
	"localRMS="<<rms<<
	"nBinPoints="<<nBinPoints<<    // center of the bin
	"binCenter.="<<&pfBinCenter<<  // 
	"normRMS="<<normRMS<<          
	"fitParam.="<<fitParam<<
	"fitCovar.="<<fitCovar<<
	"\n";
    }
  }
  return kTRUE;
}


Bool_t  AliNDLocalRegression::MakeRobustStatistic(TVectorD &values,TVectorD &errors,  TObjArray &pointArray,  TObjArray &kernelArray, Double_t weightCut, Double_t robustFraction){
  //
  // Calculate robust statistic information
  //
  TRobustEstimator e;

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
    }
    Int_t indexLocal=0;
    for (Int_t ipoint=0; ipoint<npoints; ipoint++){
      Double_t weight=1;
      for (Int_t idim=0; idim<fNParameters; idim++){
	TVectorD &vecVar=*((TVectorD*)(pointArray.UncheckedAt(idim)));
	TVectorD &vecKernel=*((TVectorD*)(kernelArray.UncheckedAt(idim)));
	fBinDelta[idim]=vecVar[ipoint]-fBinCenter[idim];       	
	weight*=TMath::Gaus(fBinDelta[idim],0,vecKernel[ipoint]);
	if (weight<weightCut) continue;	
      }      
      if (weight<weightCut) continue;
      valueLocal[indexLocal]=values[ipoint];
      indexLocal++;
    }
    Double_t median=0,meanX=0, rmsX=0;
    if (indexLocal*robustFraction-1>3){
      median=TMath::Median(indexLocal,valueLocal.GetMatrixArray());
      e.EvaluateUni(indexLocal,valueLocal.GetMatrixArray(), meanX,rmsX, indexLocal*robustFraction-1);
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
  Int_t ibin = fHistPoints->GetBin(point);
  fHistPoints->GetBinContent(ibin,fBinIndex); 
  for (Int_t idim=0; idim<fNParameters; idim++){
    fBinCenter[idim]=fHistPoints->GetAxis(idim)->GetBinCenter(fBinIndex[idim]);
  } 
  TVectorD &vecParam = *((TVectorD*)fLocalFitParam->At(ibin));
  Double_t value=vecParam[0];
  for (Int_t ipar=0; ipar<fNParameters; ipar++){
    Double_t delta=point[ipar]-fBinCenter[ipar];
    value+=vecParam[1+2*ipar]*delta;
    value+=vecParam[1+2*ipar+1]*delta*delta;
  }
  return value;
}

Double_t AliNDLocalRegression::EvalError(Double_t *point ){
  //
  //
  // 
  Int_t ibin = fHistPoints->GetBin(point);
  fHistPoints->GetBinContent(ibin,fBinIndex); 
  for (Int_t idim=0; idim<fNParameters; idim++){
    fBinCenter[idim]=fHistPoints->GetAxis(idim)->GetBinCenter(fBinIndex[idim]);
  } 
  TMatrixD &vecCovar = *((TMatrixD*)fLocalFitCovar->At(ibin));
  //TVectorD &vecQuality = *((TVectorD*)fLocalFitQuality->At(ibin));
  Double_t value=TMath::Sqrt(vecCovar(0,0));  // fill covariance to be used 
  return value;
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


