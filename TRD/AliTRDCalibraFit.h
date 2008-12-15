#ifndef ALITRDCALIBRAFIT_H
#define ALITRDCALIBRAFIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for the HLT parameters                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#  include <TObject.h>
#endif

class TTree;
class TProfile2D;
class TGraphErrors;
class TObjArray;
class TH1I;
class TH1;
class TH1F;
class TH2I;
class TH2F;
class TF1;
class TTreeSRedirector;

class AliLog;

class AliTRDCalibraMode;
class AliTRDCalibraVector;
class AliTRDCalibraVdriftLinearFit;
class AliTRDCalDet;
class AliTRDCalROC;
class AliTRDgeometry;

class AliTRDCalibraFit : public TObject {

 public: 

  // Instance
  static AliTRDCalibraFit *Instance();
  static void Terminate();
  static void Destroy();
  void DestroyDebugStreamer();

  AliTRDCalibraFit(const AliTRDCalibraFit &c);
  AliTRDCalibraFit &operator=(const AliTRDCalibraFit &) { return *this; }

       // Function for integration range of the charge 
       void     RangeChargeIntegration(Float_t vdrift, Float_t t0, Int_t &begin, Int_t &peak, Int_t &end);
  
       // Functions fit for CH
       Bool_t   AnalyseCH(TH2I *ch);
       Bool_t   AnalyseCH(AliTRDCalibraVector *calvect);
       
       // Functions fit for PH       
       Bool_t   AnalysePH(TProfile2D *ph);
       Bool_t   AnalysePH(AliTRDCalibraVector *calvect);
       
       // Functions fit for PRF
       Bool_t   AnalysePRF(TProfile2D *prf);
       Bool_t   AnalysePRF(AliTRDCalibraVector *calvect);
       
       Bool_t   AnalysePRFMarianFit(TProfile2D *prf);
       Bool_t   AnalysePRFMarianFit(AliTRDCalibraVector *calvect);
       
       // Functions fit for vdrift/lorentzangle
       Bool_t   AnalyseLinearFitters(AliTRDCalibraVdriftLinearFit *calivdli);
       
       // Pad Calibration
       Bool_t   SetModeCalibration(const char *name, Int_t i);
       
       //Reset Function
       void     ResetVectorFit();
       
       // Some functions
       Double_t *CalculPolynomeLagrange2(Double_t *x, Double_t *y) const;
       Double_t *CalculPolynomeLagrange3(Double_t *x, Double_t *y) const;
       Double_t *CalculPolynomeLagrange4(Double_t *x, Double_t *y) const;
       
       // Fill the database
       void         PutMeanValueOtherVectorFit(Int_t ofwhat = 1, Bool_t perdetector = kFALSE);
       void         PutMeanValueOtherVectorFit2(Int_t ofwhat = 1, Bool_t perdetector = kFALSE);
       AliTRDCalDet *CreateDetObjectVdrift(TObjArray *vectorFit, Bool_t perdetector = kFALSE);
       AliTRDCalDet *CreateDetObjectGain(TObjArray *vectorFit, Bool_t meanOtherBefore=kTRUE, Double_t scaleFitFactor = 0.02431, Bool_t perdetector = kTRUE);
       AliTRDCalDet *CreateDetObjectT0(TObjArray *vectorFit, Bool_t perdetector = kFALSE);
       AliTRDCalDet *CreateDetObjectLorentzAngle(TObjArray *vectorFit);
       
       TObject      *CreatePadObjectGain(TObjArray *vectorFit = 0, Double_t scaleFitFactor = 1.0, AliTRDCalDet *detobject = 0);
       TObject      *CreatePadObjectVdrift(TObjArray *vectorFit = 0, AliTRDCalDet *detobject = 0);
       TObject      *CreatePadObjectT0(TObjArray *vectorFit = 0, AliTRDCalDet *detobject = 0);
       TObject      *CreatePadObjectPRF(TObjArray *vectorFit);
       
       // Outliers stats
       AliTRDCalDet *MakeOutliersStatDet(TObjArray *vectorFit, const char *name, Double_t &mean);
       TObject      *MakeOutliersStatPad(TObjArray *vectorFit, const char *name, Double_t &mean);
       
       
       // Correct the error
       TH1F   *CorrectTheError(TGraphErrors *hist);
       
       //
       // Set or Get the variables
       //
       
       // Fit
       void     ChooseMethod(Short_t method)                              { fMethod = method;               }
       void     SetBeginFitCharge(Float_t beginFitCharge);   
       void     SetPeriodeFitPH(Int_t periodeFitPH);   
       void     SetTakeTheMaxPH()                                         { fTakeTheMaxPH   = kTRUE;        }
       void     SetT0Shift0(Float_t t0Shift0); 
       void     SetT0Shift1(Float_t t0Shift1); 
       void     SetRangeFitPRF(Float_t rangeFitPRF);     
       void     SetAccCDB()                                               { fAccCDB         = kTRUE;        }
       void     SetMinEntries(Int_t minEntries);                    
       void     SetRebin(Short_t rebin);
       
       Int_t    GetPeriodeFitPH() const                                   { return fFitPHPeriode;           }
       Bool_t   GetTakeTheMaxPH() const                                   { return fTakeTheMaxPH;           }
       Float_t  GetT0Shift0() const                                       { return fT0Shift0;               }
       Float_t  GetT0Shift1() const                                       { return fT0Shift1;               }
       Float_t  GetRangeFitPRF() const                                    { return fRangeFitPRF;            }
       Bool_t   GetAccCDB() const                                         { return fAccCDB;                 }
       Int_t    GetMinEntries() const                                     { return fMinEntries;             }
       Short_t  GetRebin() const                                          { return fRebin;                  }
       
       // Statistics
       Int_t    GetNumberFit() const                                      { return fNumberFit;              }
       Int_t    GetNumberFitSuccess() const                               { return fNumberFitSuccess;       }
       Int_t    GetNumberEnt() const                                      { return fNumberEnt;              }
       Double_t GetStatisticMean() const                                  { return fStatisticMean;          }
       
       
       // Debug
       void     SetDebugLevel(Short_t level)                              { fDebugLevel = level;            }
       void     SetDet(Int_t iLayer, Int_t iStack, Int_t iSector)         { fDet[0]  = iLayer; 
                                                                            fDet[1]  = iStack; 
                                                                            fDet[2]  = iSector;             }
       void     SetFitVoir(Int_t fitVoir)                                 { fFitVoir = fitVoir;             }
       // Magnetic field  
       void     SetMagneticField(Float_t magneticfield)                   { fMagneticField = magneticfield; }
       
       // Get the scale factor
       Double_t GetScaleFitFactor() const                                 { return fScaleFitFactor;         }
       
       // Vector Fit getter
       TObjArray  GetVectorFit() const                                    { return fVectorFit;              }
       TObjArray  GetVectorFit2() const                                   { return fVectorFit2;             }
       
       // AliTRDCalibraMode
       AliTRDCalibraMode *GetCalibraMode()                                { return fCalibraMode;            }
       
 protected:
       
       // Geometry
       AliTRDgeometry  *fGeo;               //! The TRD geometry
       
       
       Int_t        fNumberOfBinsExpected;  // Number of bins expected  
       
       // Fit
       Short_t      fMethod;                // Method
       Float_t      fBeginFitCharge;        // The fit begins at mean/fBeginFitCharge for the gain calibration
       Int_t        fFitPHPeriode;          // Periode of the fit PH
       Bool_t       fTakeTheMaxPH;          // Take the Max for the T0 determination
       Float_t      fT0Shift0;              // T0 Shift with the maximum positive slope
       Float_t      fT0Shift1;              // T0 Shift with the maximum of the amplification region
       Float_t      fRangeFitPRF;           // The fit range for the PRF is -fRangeFitPRF +fRangeFitPRF
       Bool_t       fAccCDB;                // If there is a calibration database to be compared with....
       Int_t        fMinEntries;            // Min Entries to fit the histo
       Short_t      fRebin;                 // If you want to rebin the histo for the gain calibration 
       
       // Statistics      
       Int_t        fNumberFit;             // To know how many pad groups have been fitted
       Int_t        fNumberFitSuccess;      // To know how many pad groups have been fitted successfully
       Int_t        fNumberEnt;             // To know how many pad groups have entries in the histo
       Double_t     fStatisticMean;         // To know the mean statistic of the histos
       
       
       // Debug Modes
       TTreeSRedirector   *fDebugStreamer;         //!Debug streamer
       Short_t     fDebugLevel;            // Flag for debugging
       Int_t       fDet[3];                // Detector  visualised (layer,stack,sector) si debugging == 3 or 4
       Int_t       fFitVoir;               // Fit visualised si debugging == 2
       
       // Magnetic field lorentz angle
       Float_t     fMagneticField;        // Magnetic field lorentz angle
       
       // Calibra objects
       
       AliTRDCalibraMode *fCalibraMode;  // The calibration mode
       
       // Current values of the coefficients found and ect...
       Float_t  fCurrentCoef[2];         // Current coefs  
       Float_t  fCurrentCoefE;           // Current coefs error 
       Float_t  fCurrentCoef2[2];        // Current coefs  
       Float_t  fCurrentCoefE2;          // Current coefs error   
       Float_t  fPhd[3];                 // Begin AR and DR
       Int_t    fDect1;                  // First calibration group that will be called to be maybe fitted
       Int_t    fDect2;                  // Last calibration group that will be called to be maybe fitted
       Double_t fScaleFitFactor;         // Scale factor of the fit results for the gain
       Int_t    fEntriesCurrent;         // Entries in the current histo
       Int_t    fCountDet;               // Current detector (or first in the group)
       Int_t    fCount;                  // When the next detector comes
       Int_t    fNbDet;                  // Number of detector in the group
       
       // Current calib object
       AliTRDCalDet *fCalDet;            // Current calib object
       AliTRDCalROC *fCalROC;            // Current calib object
       AliTRDCalDet *fCalDet2;           // Current calib object
       AliTRDCalROC *fCalROC2;           // Current calib object
       
       // Current values detector
       
       Float_t *fCurrentCoefDetector;     // Current values for the detector 
       Float_t *fCurrentCoefDetector2;    // Current values for the detector   
       
       class AliTRDFitInfo : public TObject {
	 
       public:
	 
	 AliTRDFitInfo()
	   :TObject()
	   ,fCoef(0x0)
	   ,fDetector(-1)                                   { }    
	 AliTRDFitInfo(const AliTRDFitInfo &i) 
	   :TObject(i)
	   ,fCoef(0x0)
	   ,fDetector(-1)                                   { }
	 AliTRDFitInfo &operator=(const AliTRDFitInfo&)     { return *this;            }
	 virtual ~AliTRDFitInfo()                           { if(fCoef) { delete [] fCoef;} }
	 
	 void      SetCoef(Float_t *coef)                   { fCoef = coef;            }
	 void      SetDetector(Int_t detector)              { fDetector = detector;    }
	 
	 Float_t  *GetCoef()                                { return fCoef;            }
	 Int_t     GetDetector() const                      { return fDetector;        }
	 
       protected:
	 
	 Float_t  *fCoef;                        // Relative coefficient for each group of the detector
	 Int_t     fDetector;                    // Detector number
    
       };
       
       TObjArray       fVectorFit;            // Vectors to fit
       TObjArray       fVectorFit2;           // Vectors to fit
       
       //
       // A lot of internal functions......
       //
       
       // Init AliTRDCalibraFit
       Bool_t   InitFit(Int_t nbins, Int_t i);
       Bool_t   InitFitCH();
       Bool_t   InitFitPH();
       Bool_t   InitFitPRF();
       Bool_t   InitFitLinearFitter();
       
       // Not enough Statistics
       Bool_t   NotEnoughStatisticCH(Int_t idect);
       Bool_t   NotEnoughStatisticPH(Int_t idect,Double_t nentries);
       Bool_t   NotEnoughStatisticPRF(Int_t idect);
       Bool_t   NotEnoughStatisticLinearFitter();
       
       // Fill Infos Fit
       Bool_t   FillInfosFitCH(Int_t idect);
       Bool_t   FillInfosFitPH(Int_t idect,Double_t nentries);
       Bool_t   FillInfosFitPRF(Int_t idect);
       Bool_t   FillInfosFitLinearFitter();
       
       void     FillFillCH(Int_t idect);
       void     FillFillPH(Int_t idect,Double_t nentries);
       void     FillFillPRF(Int_t idect);
       void     FillFillLinearFitter();
       
       Bool_t   FillVectorFit();
       Bool_t   FillVectorFit2();
       
       // Functions... 
       void     InitfCountDetAndfCount(Int_t i);
       void     CalculNumberOfBinsExpected(Int_t i);
       void     CalculDect1Dect2(Int_t i);
       void     UpdatefCountDetAndfCount(Int_t idect, Int_t i);
       void     ReconstructFitRowMinRowMax(Int_t idect, Int_t i);
       Bool_t   CheckFitVoir();
       void     NormierungCharge();
       Bool_t   SetNrphiFromTObject(const char *name, Int_t i);
       Bool_t   SetNzFromTObject(const char *name, Int_t i);
       Int_t    GetNumberOfGroupsPRF(const char* nametitle);
       
       // Calculate the mean coefs from the database
       Bool_t   CalculVdriftCoefMean();
       Bool_t   CalculChargeCoefMean(Bool_t vrai);
       Bool_t   CalculPRFCoefMean();
       Bool_t   CalculT0CoefMean();
       Bool_t   CalculVdriftLorentzCoef();
       Float_t  GetPRFDefault(Int_t layer) const;
       void     SetCalROC(Int_t i);
       
       // Fit methods
       void     FitBisCH(TH1 *projch, Double_t mean);
       void     FitCH(TH1 *projch, Double_t mean);
       void     FitMeanW(TH1 *projch, Double_t nentries);
       void     FitMeanWSm(TH1 *projch, Float_t sumAll);
       void     FitMean(TH1 *projch, Double_t nentries, Double_t mean);
       void     FitPH(TH1 *projPH, Int_t idect);
       void     FitPRF(TH1 *projPRF);
       void     RmsPRF(TH1 *projPRF);
       Bool_t   FitPRFGausMI(Double_t *arraye,Double_t *arraym,Double_t *arrayme,Int_t nBins,Float_t xMin,Float_t xMax);
       Double_t FitGausMI(Double_t *arraye,Double_t *arraym,Double_t *arrayme,Int_t nBins, Float_t xMin,Float_t xMax,TVectorD *param, Bool_t kError= kTRUE);
       void     FitPente(TH1 *projPH);
       void     FitLagrangePoly(TH1* projPH);
       void     FitTnpRange(Double_t *arraye,Double_t *arraym,Double_t *arrayme,Int_t nbg,Int_t nybins);
       TH1I    *ReBin(TH1I *hist) const;
       TH1F    *ReBin(TH1F *hist) const;
       
       // Some basic geometry function
       virtual Int_t    GetLayer(Int_t d) const;
       virtual Int_t    GetStack(Int_t d) const;
       virtual Int_t    GetSector(Int_t d) const;
       
       // Instance of this class and so on
       static  AliTRDCalibraFit   *fgInstance;                     // Instance
       static  Bool_t             fgTerminated;                    // If terminated
       
       
 private:
       
       static  Double_t PH(Double_t *x, Double_t *par);
       static  Double_t AsymmGauss(Double_t *x, Double_t *par);
       static  Double_t FuncLandauGaus(Double_t *x, Double_t *par);
       static  Double_t LanGauFun(Double_t *x, Double_t *par);
       TF1     *LanGauFit(TH1 *his, Double_t *fitrange, Double_t *startvalues
			  , Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams
			  , Double_t *fiterrors, Double_t *chiSqr, Int_t *ndf) const;
       Int_t    LanGauPro(Double_t *params, Double_t &maxx, Double_t &fwhm);
       static  Double_t GausConstant(Double_t *x, Double_t *par); 
       
       // This is a singleton, contructor is private!
       AliTRDCalibraFit();
       virtual ~AliTRDCalibraFit();
       
       
  ClassDef(AliTRDCalibraFit,2)                 // TRD Calibration class
	 
};
  
#endif


