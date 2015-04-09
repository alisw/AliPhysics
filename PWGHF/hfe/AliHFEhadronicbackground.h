#ifndef ALIHFEHADRONICBACKGROUND_H
#define ALIHFEHADRONICBACKGROUND_H

#include <TObject.h>
#include <TGraph.h>
#include <TMatrixDSym.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliPID.h"
#include "AliESDpid.h"
#include "TSpline.h"
#include "TPRegexp.h"


class AliHFEtpcPIDqa;
class TF1;
class TH1D;
class TH1;
class TH2D;
class TObjArray;
class TLegend;



class AliHFEhadronicbackground : public TObject{
  public:
    AliHFEhadronicbackground(char * run); 
    ~AliHFEhadronicbackground();
    
    void SetupElectronLandau(void);
    
    void SetUseModLandauForElectrons(Bool_t elLandau = true){useLandauForElectrons=elLandau;}
    void SetelectronExponentialParameter(Double_t parameter = true){electronExponentialParameter=parameter;}
    void SetFitRange(Double_t lower, Double_t upper){FitRangeLow=lower; FitRangeUp=upper;}

    void Process(TH2D* hSigmaTPC, Double_t pmin, Double_t pmax);

    void SetLowerCutTPC(TF1 *lcut) { fLcutTPC = lcut; }
    void SetUpperCutTPC(TF1 *ucut) { fUcutTPC = ucut; }

    void SetElectronWidth(Double_t Width);
    void SetElectronCenter(Double_t Center);

    Double_t KaonFitFunction(Double_t * dEdx, Double_t *par);
    Double_t TotalFitFunction(Double_t * dEdx, Double_t *par);
    Double_t LandauExp(Double_t * x, Double_t * par);
    Double_t ElectronLandauExp(Double_t * x, Double_t * par);
    
    void FitSlice(TH1 *hslice, Double_t momentum);
   
    void DrawDiagrams(void);

    char* GetOutputFileName(void){return fOutputFileName;}
    void SetOutputFileName(const char * newFileName);
    void SetUseGaussians(Int_t gauss=1){useGaussians=gauss;}
    void SetUseLandau(Int_t landau=1){useLandau=landau;}

    
    TF1 * returnElectronFit(int bin){return fElFits[bin-1];}
    TF1 * returnPionFit(int bin){return fPiFits[bin-1];}
    TF1 * returnKaonFit(int bin){return fKaFits[bin-1];} 
    TH1D * returnContamination(){return fContamination;}
    
    AliESDpid* SetupPIDSplinesTPC(const TString &periodin, Int_t pass=2, Bool_t isMC=kFALSE, Bool_t draw=kFALSE, Bool_t mipMom=kTRUE);
    
    TH1D * elNrDiag;

    void WriteOutputToFile(bool write=true){WriteToFile=write;}
    
  private:
    AliESDpid*PIDobject;
    bool WriteToFile;
    Bool_t useLandauForElectrons;
    Double_t electronExponentialParameter;
    double FitRangeLow;
    double FitRangeUp;

    TH1D * PiFitHist;  // for Kaon Fit
    Double_t PiHistCenter;
    TH2D * V0PionsHist;

    Int_t useGaussians;
    Int_t useLandau;
    
    Double_t * eleLanParas;

    char * fOutputFileName;

    TF1 *fLcutTPC;
    TF1 *fUcutTPC;

    TH1D ** slices; // Array of the momentum bin slices of the TPC signal
    TH1D **difference; // Absolute of the difference between signal (slices) and fit
    TH1D **ratio; // Ratio of signal and fit
    TF1 **fElFits; // fit of the electrons
    TF1 **fPiFits; // fit of the pions
    TF1 **fKaFits;
    TF1 **fAllCombined;
    
    TF1 * tempCombined;

    TH1D *fLambdaPion;
    TH1D *fCenterPi;
    TH1D *fCenterKa;
    TH1D *fCenterEl;
    TH1D *fWidthPi;
    TH1D *fWidthKa;
    TH1D *fWidthEl;
    TH1D *fAmpliPi;
    TH1D *fAmpliEl;
    TH1D *fAmpliKa;
    
    double * momentumArray;

    TH1D * fContamination;
    TH1D * fEfficiency;
    TH1D * fTRDEfficiency;

    TLegend ** FitLegend;

    Bool_t GetParametersFromFit;
    Bool_t ranCalculation;
    Bool_t fElectronWidthIsFixed;
    Bool_t fElectronCenterIsFixed;
    Double_t fFixedElectronWidth;
    Double_t fFixedElectronCenter;
    Int_t nSlices;
    Double_t * electronWidth;
    Double_t * electronCenter;
    Double_t * PionWidth;
    Double_t * PionCenter;
    Int_t nParameters;  // Number of parameters
    Double_t parMin, // Momentum at first parameter
             parStep; // momentum difference between parameter.


    ClassDef(AliHFEhadronicbackground, 1);
};
#endif

