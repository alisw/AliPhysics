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

class AliLog;
class AliTRDCalibraMode;
class AliTRDCalibraVector;
class AliTRDCalDet;

class AliTRDCalibraFit : public TObject {

 public: 

  // Instance
  static AliTRDCalibraFit *Instance();
  static void Terminate();
  static void Destroy();

  AliTRDCalibraFit(const AliTRDCalibraFit &c);
  AliTRDCalibraFit &operator=(const AliTRDCalibraFit &) { return *this; }

  // Functions fit online
          Bool_t   FitCHOnline(TH2I *ch);
	  Bool_t   FitCHOnline();
	  Bool_t   FitCHOnline(TTree *tree);
          Bool_t   FitPHOnline(TProfile2D *ph);
	  Bool_t   FitPHOnline();
	  Bool_t   FitPHOnline(TTree *tree);
          Bool_t   FitPRFOnline(TProfile2D *prf);
	  Bool_t   FitPRFOnline();
	  Bool_t   FitPRFOnline(TTree *tree);
  
  // Pad Calibration
	  Bool_t   SetModeCalibrationFromTObject(TObject *object, Int_t i);

  // Fill the database
  TObject         *CreatePadObjectTree(TTree *tree);
  TObject         *CreatePadObjectTree(TTree *tree, Int_t i, AliTRDCalDet *detobject);
  AliTRDCalDet    *CreateDetObjectTree(TTree *tree, Int_t i);

  // Correct the error
  TH1F            *CorrectTheError(TGraphErrors *hist);
    
  //
  // Set or Get the variables
  //

  // Write
          void     SetWriteCoef(Int_t i)                                     { fWriteCoef[i]  = kTRUE;         }
          void     SetWriteNameCoef(TString writeNameCoef)                   { fWriteNameCoef = writeNameCoef; }
           
          Bool_t   GetWriteCoef(Int_t i) const                               { return fWriteCoef[i];           }
          TString  GetWriteNameCoef() const                                  { return fWriteNameCoef;          }
          
  // Fit
          void     SetFitPHOn()                                              { fFitPHOn        = kTRUE;        }
	  void     SetFitPol2On()                                            { fFitPol2On      = kTRUE;        }
	  void     SetFitLagrPolOn()                                         { fFitLagrPolOn   = kTRUE;        }
	  void     SetTakeTheMaxPH()                                         { fTakeTheMaxPH   = kTRUE;        }
          void     SetPeriodeFitPH(Int_t periodeFitPH); 
	  void     SetFitPHNDB(Int_t fitPHNDB);  
          void     SetBeginFitCharge(Float_t beginFitCharge);     
          void     SetT0Shift(Float_t t0Shift); 
          void     SetRangeFitPRF(Float_t rangeFitPRF);     
	  void     SetFitPRFOn()                                             { fFitPRFOn       = kTRUE;        }
	  void     SetRMSPRFOn()                                             { fRMSPRFOn       = kTRUE;        }
	  void     SetFitPRFNDB(Int_t fitPRFNDB);    
          void     SetMeanChargeOn()                                         { fMeanChargeOn   = kTRUE;        }
	  void     SetFitChargeBisOn()                                       { fFitChargeBisOn = kTRUE;        }
	  void     SetFitChargeOn()                                          { fFitChargeOn    = kTRUE;        }
	  void     SetFitMeanWOn()                                           { fFitMeanWOn     = kTRUE;        }
	  void     SetFitChargeNDB(Int_t fitChargeNDB);
	  void     SetAccCDB()                                               { fAccCDB         = kTRUE;        }
          void     SetMinEntries(Int_t minEntries)                           { fMinEntries     = minEntries;   }
	  void     SetRebin(Short_t rebin);
  
          Bool_t   GetFitPHOn() const                                        { return fFitPHOn;                }
	  Bool_t   GetFitPol2On() const                                      { return fFitPol2On;              }
	  Bool_t   GetFitLagrPolOn() const                                   { return fFitLagrPolOn;           }
	  Bool_t   GetTakeTheMaxPH() const                                   { return fTakeTheMaxPH;           }
          Int_t    GetPeriodeFitPH() const                                   { return fFitPHPeriode;           }
	  Int_t    GetFitPHNDB() const                                       { return fFitPHNDB;               }
          Float_t  GetBeginFitCharge() const                                 { return fBeginFitCharge;         }
          Float_t  GetT0Shift() const                                        { return fT0Shift;                }
          Float_t  GetRangeFitPRF() const                                    { return fRangeFitPRF;            }
	  Bool_t   GetFitPRFOn() const                                       { return fFitPRFOn;               }
	  Bool_t   GetRMSPRFOn() const                                       { return fRMSPRFOn;               }
	  Int_t    GetFitPRFNDB() const                                      { return fFitPRFNDB;              }
          Bool_t   GetMeanChargeOn() const                                   { return fMeanChargeOn;           }
	  Bool_t   GetFitChargeBisOn() const                                 { return fFitChargeBisOn;         }
	  Bool_t   GetFitChargeOn() const                                    { return fFitChargeOn;            }
	  Bool_t   GetFitMeanWOn() const                                     { return fFitMeanWOn;             }
	  Int_t    GetFitChargeNDB() const                                   { return fFitChargeNDB;           }
	  Bool_t   GetAccCDB() const                                         { return fAccCDB;                 }
          Int_t    GetMinEntries() const                                     { return fMinEntries;             }
	  Short_t  GetRebin() const                                          { return fRebin;                  }

  // Statistics
          Int_t    GetNumberFit() const                                      { return fNumberFit;              }
	  Int_t    GetNumberFitSuccess() const                               { return fNumberFitSuccess;       }
	  Int_t    GetNumberEnt() const                                      { return fNumberEnt;              }
          Double_t GetStatisticMean() const                                  { return fStatisticMean;          }
	
  
  // Debug
          void     SetDebug(Short_t debug)                                   { fDebug   = debug;               }
          void     SetDet(Int_t iPlane, Int_t iChamb, Int_t iSect)           { fDet[0]  = iPlane; 
                                                                               fDet[1]  = iChamb; 
                                                                               fDet[2]  = iSect;               }
          void     SetFitVoir(Int_t fitVoir)                                 { fFitVoir = fitVoir;             }
  
          Short_t  GetDebug() const                                          { return fDebug;                  }
          Int_t    GetDet(Int_t i) const                                     { return fDet[i];                 }
          Int_t    GetFitVoir() const                                        { return fFitVoir;                }

  // calibration mode
	  void     SetCalibraMode(AliTRDCalibraMode *calibramode)            { fCalibraMode = calibramode;     }
AliTRDCalibraMode  *GetCalibraMode() const                                   { return fCalibraMode;            }

 
  // Getter for the coefficient trees 
          TTree   *GetPRF() const                                            { return fPRF;                    }
          TTree   *GetGain() const                                           { return fGain;                   }
          TTree   *GetT0() const                                             { return fT0;                     }
          TTree   *GetVdrift() const                                         { return fVdrift;                 }

  // Vector method
               void  SetCalibraVector(AliTRDCalibraVector *calibraVector)    { fCalibraVector = calibraVector; }
AliTRDCalibraVector *GetCalibraVector() const                                { return fCalibraVector;          }
	  

 private:
  
  static  Double_t PH(Double_t *x, Double_t *par);
  static  Double_t AsymmGauss(Double_t *x, Double_t *par);
  static  Double_t FuncLandauGaus(Double_t *x, Double_t *par);
  static  Double_t LanGauFun(Double_t *x, Double_t *par);
          TF1     *LanGauFit(TH1 *his, Double_t *fitrange, Double_t *startvalues
                           , Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams
                           , Double_t *fiterrors, Double_t *chiSqr, Int_t *ndf);
          Int_t    LanGauPro(Double_t *params, Double_t &maxx, Double_t &fwhm);
	  Double_t *CalculPolynomeLagrange2(Double_t *x, Double_t *y);
	  Double_t *CalculPolynomeLagrange3(Double_t *x, Double_t *y);
	  Double_t *CalculPolynomeLagrange4(Double_t *x, Double_t *y);
  static  Double_t GausConstant(Double_t *x, Double_t *par); 
  
  // This is a singleton, contructor is private!
  AliTRDCalibraFit();
  virtual ~AliTRDCalibraFit();

 protected:

  // Write
          Bool_t   fWriteCoef[3];           // Do you want to write the result in a file?
          TString  fWriteNameCoef;          // Where the coef Det are written
         
  // Fit
          Bool_t   fFitPHOn;                // The fit PH On (0)
	  Bool_t   fFitPol2On;              // The fit Pol2 On (1)
	  Bool_t   fFitLagrPolOn;           // The fit LagrPol On (3)
	  Bool_t   fTakeTheMaxPH;           // Take the Max for the T0 determination
          Int_t    fFitPHPeriode;           // Periode of the fit PH
	  Int_t    fFitPHNDB;               // To choose which method will be used to fill the database for the PH 
          Float_t  fBeginFitCharge;         // The fit begins at mean/fBeginFitCharge for the gain calibration
	  Float_t  fT0Shift;                // T0 Shift with the actual method
          Float_t  fRangeFitPRF;            // The fit range for the PRF is -fRangeFitPRF +fRangeFitPRF
	  Bool_t   fFitPRFOn;               // The fit PRF Gaussian On (0)
	  Bool_t   fRMSPRFOn;               // The RMS PRF On (2)
	  Int_t    fFitPRFNDB;              // To choose which method will be used to fill the database for the PRF
          Bool_t   fMeanChargeOn;           // Mean Charge on (1)
	  Bool_t   fFitChargeBisOn;         // For an other fit function (convolution, more time consuming)(2)
	  Bool_t   fFitChargeOn;            // For the first fit function (sum of Gaus and Landau) (0)
	  Bool_t   fFitMeanWOn;             // For the Marian Mean W method (4)
	  Int_t    fFitChargeNDB;           // To choose which method will be used to fill the database for the CH 
	  Bool_t   fAccCDB;                 // If there is a calibration database to be compared with....
	  Int_t    fMinEntries;             // Min Entries to fit the histo
	  Short_t  fRebin;                  // If you want to rebin the histo for the gain calibration
         
  // Statistics      
          Int_t    fNumberFit;              // To know how many pad groups have been fitted
	  Int_t    fNumberFitSuccess;       // To know how many pad groups have been fitted successfully
	  Int_t    fNumberEnt;              // To know how many pad groups have entries in the histo
          Double_t fStatisticMean;          // To know the mean statistic of the histos

  // Debug Mode
          Short_t  fDebug;                  // For debugging 0 rien, 1 errors, 2 one fit alone, 3 one detector, 4 one detector with errors
          Int_t    fDet[3];                 // Detector  visualised (plane,chamb,sect) si debugging == 3 or 4
          Int_t    fFitVoir;                // Fit visualised si debugging == 2
  
  // Calibration mode

	  AliTRDCalibraMode *fCalibraMode;  // The calibration mode

  // The coefficients trees

          TTree   *fPRF;                    // Tree of the sigma of PRD
          TTree   *fGain;                   // Tree of the gain factor
          TTree   *fT0;                     // Tree of the time0
          TTree   *fVdrift;                 // Tree of the drift velocity

  // "Pointer" of the branch of the tree
          Int_t    fVdriftDetector;         // Branch of Vdrift
          Float_t *fVdriftPad;              // Branch of Vdrift
          Int_t    fT0Detector;             // Branch of t0
          Float_t *fT0Pad;                  // Branch of t0
          Int_t    fPRFDetector;            // Branch of PRF
          Float_t *fPRFPad;                 // Branch of PRF
          Float_t *fCoefCH;                 // Branch relative gain
 
  //
  // For debugging 
  //

  // To build the graph with the errors of the fits
	  Double_t        *fCoefCharge[5];   // Coefs resulting from the fit for the gain
	  Double_t        *fCoefChargeE[4];  // Error of the found coefs for the gain
	  Double_t        *fCoefVdrift[4];   // Coefs resulting from the fit for the drift velocity
	  Double_t        *fCoefVdriftE[3];  // Error of the found coefs for the drift velocity
	  Double_t        *fCoefT0[4];       // Coefs resulting from the fit for the drift velocity
	  Double_t        *fCoefT0E[3];      // Error of the found coefs for the drift velocity
	  Double_t        *fCoefPRF[3];      // Coefs resulting from the fit for the PRF
	  Double_t        *fCoefPRFE[2];     // Error of the found coefs for the PRF
	  TH2F    *fCoefChargeDB[4];         // Visualisation of the coef of the detecteur fDet for the gain
          TH2F    *fCoefVdriftDB[3];         // Visualisation of the coef of the detecteur fDet for the drift velocity
          TH2F    *fCoefT0DB[3];             // Visualisation of the coef of the detecteur fDet for time 0
          TH2F    *fCoefPRFDB[2];            // Visualisation of the coef of the detecteur fDet for the pad response function

  // Variables in the loop for the coef or more general
          Float_t  fChargeCoef[5];          // 4 Marian Mean W, 3 database value, 0 fit, 1 mean, 2 fit time consuming   
          Float_t  fVdriftCoef[4];          // 3 lagrangepoly, 2 database value, 1 slope method, 0 fit
          Float_t  fPRFCoef[3];             // 2 Rms, 1 database, 0 fit 
          Float_t  fT0Coef[4];              // 3 lagrangepoly, 2 database, 1 slope method, 0 fit
          Float_t  fPhd[3];                 // Begin AR and DR
	  Int_t    fDect1[3];               // First calibration group that will be called to be maybe fitted
          Int_t    fDect2[3];               // Last calibration group that will be called to be maybe fitted
          Double_t fScaleFitFactor;         // Scale factor of the fit results for the gain
	  Int_t    fEntriesCurrent;         // Entries in the current histo
	  Int_t    fCountDet[3];            // Current detector
          Int_t    fCount[3];               // When the next detector comes
 
  // Vector method

	  AliTRDCalibraVector *fCalibraVector; // The vector object

	  class AliTRDFitCHInfo : public TObject {

	  public:
	    
	    AliTRDFitCHInfo()
	      :TObject()
	      ,fCoef(0x0)
	      ,fDetector(-1)                                   { }    
	    AliTRDFitCHInfo(const AliTRDFitCHInfo &i) 
	      :TObject(i)
	      ,fCoef(0x0)
	      ,fDetector(-1)                                   { }
	    AliTRDFitCHInfo &operator=(const AliTRDFitCHInfo&) { return *this;            }
	    virtual ~AliTRDFitCHInfo()                         { }
	    
	    void      SetCoef(Float_t *coef)                   { fCoef     = coef;        }
	    void      SetDetector(Int_t detector)              { fDetector = detector;    }
	    
	    Float_t  *GetCoef() const                          { return fCoef;            }
	    Int_t     GetDetector() const                      { return fDetector;        }
	    
	  protected:
	    
	    Float_t  *fCoef;                        // Relative gain coefficient for each group of the detector
	    Int_t     fDetector;                    // Detector number
	    
	  };

	  TObjArray       *fVectorFitCH;            // Vectors to fit
  
  //
  // A lot of internal functions......
  //

  // Init AliTRDCalibraFit
          void     Init();
    
  //
  // Fit
  //  

  // Create histos if fDebug == 1 or fDebug >=3
          void     CreateFitHistoPHDB(Int_t rowMax, Int_t colMax);
          void     CreateFitHistoT0DB(Int_t rowMax, Int_t colMax);
          void     CreateFitHistoCHDB(Int_t rowMax, Int_t colMax);
          void     CreateFitHistoPRFDB(Int_t rowMax, Int_t colMax);
          void     InitArrayFitCH();
          void     InitArrayFitPH();
          void     InitArrayFitT0();
          void     InitArrayFitPRF();
  
  // CHFit functions
          Bool_t   FillVectorFitCH(Int_t countdet);
          Bool_t   InitFit(Int_t nbins, Int_t i);
          void     InitfCountDetAndfCount(Int_t i);
          void     UpdatefCountDetAndfCount(Int_t idect, Int_t i);
          void     ReconstructFitRowMinRowMax(Int_t idect, Int_t i);
          Bool_t   NotEnoughStatistic(Int_t idect, Int_t i);
          Bool_t   FillInfosFit(Int_t idect, Int_t i);
          Bool_t   WriteFitInfos(Int_t i);
          void     NormierungCharge();

  // Fill histos DB from the Coef histos 
          void     FillCoefChargeDB();
          void     FillCoefVdriftDB();
          void     FillCoefT0DB();
          void     FillCoefPRFDB();

  // Plot histos CoefPRF Coef....
          void     PlotWritePH();
          void     PlotWriteT0();
          void     PlotWriteCH();
          void     PlotWritePRF();
  
  // Plot histos DB
          void     PlotPHDB();
          void     PlotT0DB();
          void     PlotCHDB();
          void     PlotPRFDB();
  
  // Write the DB histos
          void     WritePHDB(TFile *fout);
          void     WriteT0DB(TFile *fout);
          void     WriteCHDB(TFile *fout);
          void     WritePRFDB(TFile *fout);

  // Calculate the mean coefs from the database
          Bool_t   CalculVdriftCoefMean(Int_t fect, Int_t idect);
          Bool_t   CalculChargeCoefMean(Int_t fect, Int_t idect, Bool_t vrai);
          Bool_t   CalculPRFCoefMean(Int_t fect, Int_t idect);
          Bool_t   CalculT0CoefMean(Int_t fect, Int_t idect);
          Float_t  GetPRFDefault(Int_t plane) const;

 // Fit methods
          void     FitBisCH(TH1 *projch, Int_t idect);
          void     FitCH(TH1 *projch, Int_t idect);
	  void     FitMeanW(TH1 *projch, Int_t idect);
	  void     FitMean(TH1 *projch, Int_t idect, Double_t nentries);
          void     FitPH(TH1 *projPH, Int_t idect);
          void     FitPRF(TH1 *projPRF, Int_t idect);
	  void     RmsPRF(TH1 *projPRF, Int_t idect);
          void     FitPente(TH1 *projPH, Int_t idect);
	  void     FitLagrangePoly(TH1* projPH, Int_t idect);
          TH1I    *ReBin(TH1I *hist) const;
          TH1F    *ReBin(TH1F *hist) const;
  
  // Clear
	  void     ClearTree();

  // Some basic geometry function
  virtual Int_t    GetPlane(Int_t d) const;
  virtual Int_t    GetChamber(Int_t d) const;
  virtual Int_t    GetSector(Int_t d) const;
  
  // Init, Fill and Reset the variables to default value tree Gain, PRF, Vdrift and T0
          void     InitTreePH();
          void     FillTreeVdrift(Int_t countdet);
          void     InitTreeT0();
          void     FillTreeT0(Int_t countdet);
          void     InitTreePRF();
          void     FillTreePRF(Int_t countdet);
          void     ConvertVectorFitCHTree();

  // Instance of this class and so on
  static  AliTRDCalibraFit *fgInstance;        // Instance
  static  Bool_t   fgTerminated;               // If terminated
    
  ClassDef(AliTRDCalibraFit,1)                 // TRD Calibration class

};
  
#endif


