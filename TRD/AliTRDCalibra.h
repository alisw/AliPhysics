#ifndef ALITRDCALIBRA_H
#define ALITRDCALIBRA_H
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
class AliTRDcluster;
class AliTRDtrack;
class AliTRDmcmTracklet;
class AliTRDCalDet;

class AliTRDCalibra : public TObject {

 public: 

  // Instance
  static AliTRDCalibra *Instance();
  static void Terminate();
  static void Destroy();

  AliTRDCalibra(const AliTRDCalibra &c);
  AliTRDCalibra &operator=(const AliTRDCalibra &) { return *this; }

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
 
  // Functions for initialising the AliTRDCalibra in the code
          Bool_t   Init2Dhistos();

  // Functions for filling the histos in the code
          Bool_t   ResetTrack();
          Bool_t   UpdateHistograms(AliTRDcluster *cl, AliTRDtrack *t);
          Bool_t   UpdateHistogramcm(AliTRDmcmTracklet *trk);
 
  // Is Pad on
          Bool_t   IsPadOn(Int_t detector, Int_t col, Int_t row) const;

  // Functions for plotting the 2D
          void     Plot2d();

  // Functions for writting the 2D
          Bool_t   Write2d();

  // Function fill 2D for the moment out of the code
          Bool_t   Create2DDiSimOnline(Int_t iev1, Int_t iev2);
          Bool_t   Create2DRaDaOnline(Int_t iev1, Int_t iev2);

  // Pad Calibration
          Bool_t   ModePadFragmentation(Int_t iPlane,Int_t iChamb, Int_t iSect, Int_t i);
          void     ModePadCalibration(Int_t iChamb, Int_t i);
          Bool_t   SetModeCalibrationFromTObject(TObject *object, Int_t i);

  // Fill the database
  TObject         *CreatePadObjectTree(TTree *tree);
  TObject         *CreatePadObjectTree(TTree *tree, Int_t i, AliTRDCalDet *detobject);
  AliTRDCalDet    *CreateDetObjectTree(TTree *tree, Int_t i);

  // Correct the error
  TH1F            *CorrectTheError(TGraphErrors *hist);
  TGraphErrors    *AddProfiles(TGraphErrors *hist1, TGraphErrors *hist2) const ;

  // Add two trees
  TTree           *Sum2Trees(const Char_t *filename1, const Char_t *filename2, const Char_t *variablecali);
  
  //
  // Set of Get the variables
  //

  // Choice to fill or not the 2D
          void     SetMITracking(Bool_t mitracking = kTRUE)                  { fMITracking      = mitracking;  }
          void     SetMcmTracking(Bool_t mcmtracking = kTRUE)                { fMcmTracking     = mcmtracking;  }
          void     SetMcmCorrectAngle()                                      { fMcmCorrectAngle = kTRUE; }
          void     SetPH2dOn()                                               { fPH2dOn          = kTRUE; }
          void     SetCH2dOn()                                               { fCH2dOn          = kTRUE; }
          void     SetPRF2dOn()                                              { fPRF2dOn         = kTRUE; }
          void     SetHisto2d()                                              { fHisto2d         = kTRUE; }
          void     SetVector2d()                                             { fVector2d        = kTRUE; }
  
          Bool_t   GetMITracking() const                                     { return fMITracking;       }
          Bool_t   GetMcmTracking() const                                    { return fMcmTracking;      }
          Bool_t   GetMcmCorrectAngle() const                                { return fMcmCorrectAngle;  }
          Bool_t   GetPH2dOn() const                                         { return fPH2dOn;           }
          Bool_t   GetCH2dOn() const                                         { return fCH2dOn;           }
          Bool_t   GetPRF2dOn() const                                        { return fPRF2dOn;          }
          Bool_t   GetHisto2d() const                                        { return fHisto2d;          }
          Bool_t   GetVector2d() const                                       { return fVector2d;         }
  TH2I            *GetCH2d() const                                           { return fCH2d;             }
  TProfile2D      *GetPH2d() const                                           { return fPH2d;             }
  TProfile2D      *GetPRF2d() const                                          { return fPRF2d;            }
  
  // How to fill the 2D
          void     SetRelativeScaleAuto()                                    { fRelativeScaleAuto    = kTRUE;                      }
          void     SetRelativeScale(Float_t relativeScale);                      
          void     SetThresholdDigit(Int_t digitthreshold)                   { fThresholdDigit       = digitthreshold;             }
          void     SetThresholdClusterPRF1(Float_t thresholdClusterPRF1)     { fThresholdClusterPRF1 = thresholdClusterPRF1;       }
          void     SetThresholdClusterPRF2(Float_t thresholdClusterPRF2)     { fThresholdClusterPRF2 = thresholdClusterPRF2;       }
          void     SetCenterOfflineCluster()                                 { fCenterOfflineCluster = kTRUE;                      }
          void     SetTraMaxPad()                                            { fTraMaxPad            = kTRUE;                      }
          void     SetNz(Int_t i, Short_t nz);
          void     SetNrphi(Int_t i, Short_t nrphi);
          void     SetProcent(Float_t procent)                               { fProcent              = procent;                    }
          void     SetDifference(Short_t difference)                         { fDifference           = difference;                 }
          void     SetNumberClusters(Short_t numberClusters)                 { fNumberClusters       = numberClusters;             }
          void     SetNumberBinCharge(Short_t numberBinCharge)               { fNumberBinCharge      = numberBinCharge;            }
          void     SetNumberBinPRF(Short_t numberBinPRF)                     { fNumberBinPRF         = numberBinPRF;               }
  
          Float_t  GetRelativeScale() const                                  { return fRelativeScale;        }
          Bool_t   GetRelativeScaleAuto() const                              { return fRelativeScaleAuto;    }
          Int_t    GetThresholdDigit() const                                 { return fThresholdDigit;       }
          Float_t  GetThresholdClusterPRF1() const                           { return fThresholdClusterPRF1; }
          Float_t  GetThresholdClusterPRF2() const                           { return fThresholdClusterPRF2; }
          Bool_t   GetTraMaxPad()const                                       { return fTraMaxPad;            }
          Short_t  GetNz(Int_t i) const                                      { return fNz[i];                }
          Short_t  GetNrphi(Int_t i) const                                   { return fNrphi[i];             }
          Float_t  GetProcent() const                                        { return fProcent;              }
          Short_t  GetDifference() const                                     { return fDifference;           }
          Short_t  GetNumberClusters() const                                 { return fNumberClusters;       }
          Short_t  GetNumberBinCharge() const                                { return fNumberBinCharge;      }
          Short_t  GetNumberBinPRF() const                                   { return fNumberBinPRF;         }
  
  // Write
          void     SetWriteCoef(Int_t i)                                     { fWriteCoef[i]  = kTRUE;         }
          void     SetWriteNameCoef(TString writeNameCoef)                   { fWriteNameCoef = writeNameCoef; }
          void     SetWrite(Int_t i)                                         { fWrite[i]      = kTRUE;         }
          void     SetWriteName(TString writeName)                           { fWriteName     = writeName;     }
  
          Bool_t   GetWriteCoef(Int_t i) const                               { return fWriteCoef[i];           }
          TString  GetWriteNameCoef() const                                  { return fWriteNameCoef;          }
          Bool_t   GetWrite(Int_t i) const                                   { return fWrite[i];               }
          TString  GetWriteName() const                                      { return fWriteName;              }
  
  // Fit
          void     SetFitPHOn()                                              { fFitPHOn        = kTRUE;}
          void     SetPeriodeFitPH(Int_t periodeFitPH);   
          void     SetBeginFitCharge(Float_t beginFitCharge);     
          void     SetT0Shift(Float_t t0Shift); 
          void     SetRangeFitPRF(Float_t rangeFitPRF);       
          void     SetMeanChargeOn()                                         { fMeanChargeOn   = kTRUE;      }
          void     SetAccCDB()                                               { fAccCDB         = kTRUE;      }
          void     SetFitChargeBisOn()                                       { fFitChargeBisOn = kTRUE;      }
          void     SetMinEntries(Int_t minEntries)                           { fMinEntries     = minEntries; }
  
          Bool_t   GetFitPHOn() const                                        { return fFitPHOn;        }
          Int_t    GetPeriodeFitPH() const                                   { return fFitPHPeriode;   }
          Float_t  GetBeginFitCharge() const                                 { return fBeginFitCharge; }
          Float_t  GetT0Shift() const                                        { return fT0Shift;        }
          Float_t  GetRangeFitPRF() const                                    { return fRangeFitPRF;    }
          Bool_t   GetMeanChargeOn() const                                   { return fMeanChargeOn;   }
          Bool_t   GetAccCDB() const                                         { return fAccCDB;         }
          Bool_t   GetFitChargeBisOn() const                                 { return fFitChargeBisOn; }
          Int_t    GetMinEntries() const                                     { return fMinEntries;     }
          Int_t    GetNumberFit() const                                      { return fNumberFit;      }
          Double_t GetStatisticMean() const                                  { return fStatisticMean;  }
  
  // Debug
          void     SetDebug(Short_t debug)                                   { fDebug   = debug;   }
          void     SetDet(Int_t iPlane, Int_t iChamb, Int_t iSect)           { fDet[0]  = iPlane; 
                                                                               fDet[1]  = iChamb; 
                                                                               fDet[2]  = iSect;   }
          void     SetFitVoir(Int_t fitVoir)                                 { fFitVoir = fitVoir; }
  
          Short_t  GetDebug() const                                          { return fDebug;      }
          Int_t    GetDet(Int_t i) const                                     { return fDet[i];     }
          Int_t    GetFitVoir() const                                        { return fFitVoir;    }

  //
  // Internal variables to be sure!
  //
  
  // Pad calibration
          Short_t  GetNnz(Int_t i) const                                     { return fNnZ[i];       }
          Short_t  GetNnrphi(Int_t i) const                                  { return fNnRphi[i];    }
          Short_t  GetNfragz(Int_t i) const                                  { return fNfragZ[i];    }
          Short_t  GetNfragrphi(Int_t i) const                               { return fNfragRphi[i]; }
          Short_t  GetDetChamb0(Int_t i) const                               { return fDetChamb0[i]; }
          Short_t  GetDetChamb2(Int_t i) const                               { return fDetChamb2[i]; }
    
          void     SetRebin(Short_t rebin);
          Short_t  GetRebin() const                                          { return fRebin;        }

  // Getter for the coefficient trees 
          TTree   *GetPRF() const                                            { return fPRF;          }
          TTree   *GetGain() const                                           { return fGain;         }
          TTree   *GetT0() const                                             { return fT0;           }
          TTree   *GetVdrift() const                                         { return fVdrift;       }

 private:
  
  static  Double_t PH(Double_t *x, Double_t *par);
  static  Double_t AsymmGauss(Double_t *x, Double_t *par);
  static  Double_t FuncLandauGaus(Double_t *x, Double_t *par);
  static  Double_t LanGauFun(Double_t *x, Double_t *par);
          TF1     *LanGauFit(TH1 *his, Double_t *fitrange, Double_t *startvalues
                           , Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams
                           , Double_t *fiterrors, Double_t *chiSqr, Int_t *ndf);
          Int_t    LanGauPro(Double_t *params, Double_t &maxx, Double_t &fwhm); 
  
  // This is a singleton, contructor is private!
  AliTRDCalibra();
  virtual ~AliTRDCalibra();

 protected:

  // Choice to fill or not the 2D
          Bool_t   fMITracking;             // Chose to fill the 2D histos or vectors during the offline MI tracking
          Bool_t   fMcmTracking;            // Chose to fill the 2D histos or vectors during the tracking with tracklets
          Bool_t   fMcmCorrectAngle;        // Apply correction due to the mcmtrackletangle in the z direction (only) assuming  from vertex
          Bool_t   fCH2dOn;                 // Chose to fill the 2D histos or vectors for the relative gain calibration 
          Bool_t   fPH2dOn;                 // Chose to fill the 2D histos or vectors for the drift velocity and T0
          Bool_t   fPRF2dOn;                // Chose to fill the 2D histos or vectors for the pad response function calibration
          Bool_t   fHisto2d;                // Chose to fill the 2D histos
          Bool_t   fVector2d;               // Chose to fill vectors

  // How to fill the 2D
          Float_t  fRelativeScale;          // Scale of the deposited charge
          Int_t    fCountRelativeScale;     // fCountRelativeScale first data used for the scaling
          Bool_t   fRelativeScaleAuto;      // Scaling with the first fCountRelativeScale objects
          Int_t    fThresholdDigit;         // Threshold on RawData
          Float_t  fThresholdClusterPRF1;   // Threshold on cluster pad signals for PRF peripherique
          Float_t  fThresholdClusterPRF2;   // Threshold on cluster pad signals for PRF peripherique
          Bool_t   fCenterOfflineCluster;   // Choose to use the offline determination of the center of the cluster
          Bool_t   fTraMaxPad;              // Take the Max Pad for the gain calibration and PH
          Short_t  fNz[3];                  // Mode of calibration 
          Short_t  fNrphi[3];               // Mode of calibration 
          Int_t    fNtotal[3];              // Total number of Xbins

  // Write
          Bool_t   fWriteCoef[3];           // Do you want to write the result in a file?
          TString  fWriteNameCoef;          // Where the coef Det are written
          Bool_t   fWrite[3];               // Do you want to write the 2D histo or vectors converted in a tree
          TString  fWriteName;              // Where the 2D or trees are written
  
  // Fit
          Bool_t   fFitPHOn;                // The fit PH On
          Int_t    fFitPHPeriode;           // Periode of the fit PH
          Float_t  fBeginFitCharge;         // The fit begins at mean/fBeginFitCharge for the gain calibration
          Float_t  fRangeFitPRF;            // The fit range for the PRF is -fRangeFitPRF +fRangeFitPRF
          Bool_t   fMeanChargeOn;           // Mean Charge on
          Bool_t   fFitChargeBisOn;         // For an other fit function (convolution and not sum, more time consuming)
          Float_t  fT0Shift;                // T0 Shift with the actual method
          Bool_t   fAccCDB;                 // If there is a calibration database to be compared with....
          Int_t    fNumberFit;              // To know how many pad groups have been fitted
          Double_t fStatisticMean;          // To know the mean statistic of the histos

  // Debug Mode
          Short_t  fDebug;                  // For debugging 0 rien, 1 errors, 2 one fit alone, 3 one detector, 4 one detector with errors
          Int_t    fDet[3];                 // Detector  visualised (plane,chamb,sect) si debugging == 3 or 4
          Int_t    fFitVoir;                // Fit visualised si debugging == 2
  
  //
  // Internal variables
  //

  // Storage of coef
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
  
  // Fill the 2D histos in the offline tracking
          Bool_t   fDetectorAliTRDtrack;    // Change of track
          Int_t    fChamberAliTRDtrack;     // Change of chamber
          Int_t    fDetectorPreviousTrack;  // Change of detector
          Bool_t   fGoodTrack;              // If goes through a kaputt pad
          Float_t *fAmpTotal;               // Energy deposited in the calibration group by the track
          Short_t *fPHPlace;                // Calibration group of PH
          Float_t *fPHValue;                // PH
          Short_t  fNumberClusters;         // Minimum number of clusters in the tracklets
          Float_t  fProcent;                // Limit to take the info of the most important calibration group if the track goes through 2 groups (CH)
          Short_t  fDifference;             // Limit to take the info of the most important calibration group if the track goes through 2 groups (CH)
          Int_t    fNumberTrack;            // How many tracks could be used (Debug for the moment)
          Int_t    fNumberUsedCh[2];        // How many tracks have been really used for the gain (0, strict; 1 with fProcent)
          Int_t    fNumberUsedPh[2];        // How many tracks have been really used for the drift velocity (0, strict; 1 with fDifference)

  //
  // For debugging 
  //

  // Histograms to store the coef
          TH1F    *fCoefCharge[4];          // Replica des 2D but in coefs resulting from the fit for the gain
          TH1F    *fCoefVdrift[3];          // Replica des 2D but in coefs resulting from the fit for the drift velocity
          TH1F    *fCoefPRF[2];             // Replica des 2D but in coefs resulting from the fit for the pad response function
          TH1F    *fCoefT0[3];              // Replica des 2D but in coefs resulting from the fit for time 0
          TH1F    *fDeltaCharge[3];         // Replica des 2D but in errors for each detector resulting from the fit for the gain
          TH1F    *fDeltaVdrift[2];         // Replica des 2D but in errors for each detector resulting from the fit for the drift velocity
          TH1F    *fDeltaT0[2];             // Replica des 2D but in errors for each detector resulting from the fit for time 0
          TH1F    *fDeltaPRF;               // Replica des 2D but in errors for each detector resulting from the fit for the pad response function
          TH1I    *fErrorCharge[3];         // Replica des 2D but in errors resulting from the fit for the gain
          TH1I    *fErrorVdrift[2];         // Replica des 2D but in errors resulting from the fit for the drift velocity
          TH1I    *fErrorT0[2];             // Replica des 2D but in errors resulting from the fit for time 0
          TH1I    *fErrorPRF;               // Replica des 2D but in errors resulting from the fit for the pad response function
          TH2F    *fCoefChargeDB[3];        // Visualisation of the coef of the detecteur fDet for the gain
          TH2F    *fCoefVdriftDB[2];        // Visualisation of the coef of the detecteur fDet for the drift velocity
          TH2F    *fCoefT0DB[2];            // Visualisation of the coef of the detecteur fDet for time 0
          TH2F    *fCoefPRFDB;              // Visualisation of the coef of the detecteur fDet for the pad response function

  // Variables in the loop for the coef or more general
          Float_t  fChargeCoef[4];          // 3 database value, 0 fit, 1 mean, 2 fit time consuming   
          Float_t  fVdriftCoef[3];          // 2 database value, 1 slope method, 0 fit
          Float_t  fPRFCoef[2];             // 1 database, 0 fit 
          Float_t  fT0Coef[3];              // 3 database, 1 slope method, 0 fit
          Float_t  fPhd[3];                 // Begin AR and DR
          Int_t    fTimeMax;                // Number of time bins
          Float_t  fSf;                     // Sampling frequence
          Int_t    fDect1[3];               // First calibration group that will be called to be maybe fitted
          Int_t    fDect2[3];               // Last calibration group that will be called to be maybe fitted
          Double_t fScaleFitFactor;         // Scale factor of the fit results for the gain
          Int_t    fMinEntries;             // Min Entries to fit the histo
          Int_t    fEntriesCurrent;         // Entries in the current histo
          Int_t    fCountDet[3];            // Current detector
          Int_t    fCount[3];               // When the next detector comes
          Float_t  fL3P0;                   // Parameter to be pass from the default fit of CH histo to the optional one
          Float_t  fL3P2;                   // Parameter to be pass from the default fit of CH histo to the optional one
          Float_t  fG3P2;                   // Parameter to be pass from the default fit of CH histo to the optional one

  // Pad Calibration
          Short_t  fNnZ[3];                 // Number of pad rows in a group
          Short_t  fNnRphi[3];              // Number of pad cols in a group
          Short_t  fNfragZ[3];              // Number of  pad row group
          Short_t  fNfragRphi[3];           // Number of pad col group
          Short_t  fRowMin[3];              // Limits of the group in pad row
          Short_t  fRowMax[3];              // Limits of the group in pad row
          Short_t  fColMin[3];              // Limits of the group in pad col
          Short_t  fColMax[3];              // Limits of the group in pad col
          Int_t    fXbins[3];               // First Xbins of the detector
          Short_t  fDetChamb0[3];           // Number of XBins for chamber != 2
          Short_t  fDetChamb2[3];           // Number of Xbins fir chamber 2
  
  // Methode  Alexandru store info
  class AliTRDPlace : public TObject {

   public: 
    
    AliTRDPlace()
      :TObject()
      ,fPlace(0x0)                                     { }
    AliTRDPlace(const AliTRDPlace &i)
      :TObject(i)
      ,fPlace(0x0)                                     { }
    AliTRDPlace &operator=(const AliTRDPlace&)         { return *this;            } 
    virtual ~AliTRDPlace()                             { }
    
    void      SetPlace(Int_t place)                    { fPlace = place;          }
    Int_t     GetPlace() const                         { return  fPlace;          }
    
  protected:
    
    Int_t     fPlace;                       // Place of the calibration group

  };
  
  class AliTRDCTInfo : public TObject {

   public: 
    
    AliTRDCTInfo()
      :TObject()
      ,fEntries(0x0)                                   { }
    AliTRDCTInfo(const AliTRDCTInfo &i)
      :TObject(i)
      ,fEntries(0x0)                                   { }
    AliTRDCTInfo &operator=(const AliTRDCTInfo&)       { return *this;            } 
    virtual ~AliTRDCTInfo()                            { }
    
    void      SetEntries(UShort_t *entries)            { fEntries = entries;      }

    UShort_t *GetEntries() const                       { return fEntries;         }
    
  protected:
    
    UShort_t *fEntries;                     // Current number of entries for each bin of CH

  };

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
  
  class AliTRDPInfo : public TObject {
  public:
    
    AliTRDPInfo()
      :TObject()
      ,fSum(0x0) 
      ,fSumSquare(0x0)
      ,fEntries(0x0)                                   { }
    AliTRDPInfo(const AliTRDPInfo &i)
      :TObject(i)
      ,fSum(0x0) 
      ,fSumSquare(0x0)
      ,fEntries(0x0)                                   { } 
    AliTRDPInfo &operator=(const AliTRDPInfo&)         { return *this;            }
    virtual ~AliTRDPInfo()                             { }
    
    void      SetSum(Float_t *sum)                     { fSum       = sum;        }
    void      SetSumSquare(Float_t *sumSquare)         { fSumSquare = sumSquare;  }
    void      SetEntries(UShort_t *entries)            { fEntries   = entries;    }
    
    Float_t  *GetSum() const                           { return fSum;             }
    Float_t  *GetSumSquare() const                     { return fSumSquare;       }
    UShort_t *GetEntries() const                       { return fEntries;         }
    
  protected:

    Float_t  *fSum;                         // Current mean for each bin of the average pulse height
    Float_t  *fSumSquare;                   // Current mean of square values for each bin of the average pulse height
    UShort_t *fEntries;                     // Current number of entries for each bin of the average pulse height

  }; 
   
  // PH
  // fTimeMax will define the size of fcharge
  TObjArray       *fVectorPH;               // Vectors to fill
  TObjArray       *fPlaPH;                  // Vectors to fill
  // CH
          Short_t  fNumberBinCharge;        // Number of bins for the gain factor
  TObjArray       *fVectorCH;               // Vectors to fill
  TObjArray       *fPlaCH;                  // Vectors to fill
  // FitCH
  TObjArray       *fVectorFitCH;            // Vectors to fit
  // PRF
  Short_t          fNumberBinPRF;           // Number of bin for the PRF
  TObjArray       *fVectorPRF;              // Vectors to fill
  TObjArray       *fPlaPRF;                 // Vectors to fill
 
  // Histograms to store the info from the digits, from the tracklets or from the tracks
  TProfile2D      *fPH2d;                   // 2D average pulse height
  TProfile2D      *fPRF2d;                  // 2D PRF
  TH2I            *fCH2d;                   // 2D deposited charge 
          Short_t  fRebin;                  // If you want to rebin the histo for the gain calibration
  
  //
  // A lot of internal functions......
  //

  // Init AliTRDCalibra
          void     Init();

  // Create the 2D histo to be filled Online
          void     CreateCH2d(Int_t nn);
          void     CreatePH2d(Int_t nn);
          void     CreatePRF2d(Int_t nn);  
  
  // Fill the 2D
          void     FillTheInfoOfTheTrackPH();
          void     FillTheInfoOfTheTrackCH();
          void     ResetfVariables();
          Bool_t   LocalisationDetectorXbins(Int_t detector);
 
  // Plot the 2D
          void     PlotCH2d();
          void     PlotPH2d();
          void     PlotPRF2d();
    
  //
  // Fit
  //  

  // Create histos if fDebug == 1 or fDebug >=3
          void     CreateFitHistoPHDB(Int_t rowMax, Int_t colMax);
          void     CreateFitHistoT0DB(Int_t rowMax, Int_t colMax);
          void     CreateFitHistoCHDB(Int_t rowMax, Int_t colMax);
          void     CreateFitHistoPRFDB(Int_t rowMax, Int_t colMax);
          void     CreateFitHistoCH(Int_t nbins, Double_t low, Double_t high);
          void     CreateFitHistoPH(Int_t nbins, Double_t low, Double_t high);
          void     CreateFitHistoT0(Int_t nbins, Double_t low, Double_t high);
          void     CreateFitHistoPRF(Int_t nbins, Double_t low, Double_t high);
  
  // CHFit functions
          Bool_t   FillVectorFitCH(Int_t countdet);
          Bool_t   InitFit(Int_t nbins, Double_t lowedge, Double_t upedge, Int_t i);
          void     InitfCountDetAndfCount(Int_t i);
          void     UpdatefCountDetAndfCount(Int_t idect, Int_t i);
          void     ReconstructFitRowMinRowMax(Int_t idect, Int_t i);
          Bool_t   NotEnoughStatistic(Int_t idect, Int_t i);
          Bool_t   FillInfosFit(Int_t idect, Int_t i);
          Bool_t   WriteFitInfos(Int_t i);
          void     NormierungCharge();

  // Fill histos Errors from the delta histos
          void     ErrorPH();
          void     ErrorT0();
          void     ErrorCH();
          void     ErrorPRF();

  // Fill histos DB from the Coef histos 
          void     FillCoefChargeDB();
          void     FillCoefVdriftDB();
          void     FillCoefT0DB();
          void     FillCoefPRFDB();

  // Plot histos CoefPRF Coef....
          void     PlotPH();
          void     PlotT0();
          void     PlotCH();
          void     PlotPRF();
  
  // Plot histos DB
          void     PlotPHDB();
          void     PlotT0DB();
          void     PlotCHDB();
          void     PlotPRFDB();

  // Write the Coef, delta and error histos
          void     WritePH(TFile *fout);
          void     WriteT0(TFile *fout);
          void     WriteCH(TFile *fout);
          void     WritePRF(TFile *fout);
    
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

  // Pad group calibration mode
          void     ReconstructionRowPadGroup(Int_t idect, Int_t i);
          void     CalculXBins(Int_t idect, Int_t i);  

  // Convertion vector, tree, histos....
          Int_t    SearchInVector(Int_t group, Int_t i) const;
          Int_t    SearchInTreeVector (TObjArray *vectorplace, Int_t group) const;
          Int_t    SearchBin (Float_t value, Int_t i) const;
          Bool_t   UpdateVectorCH(Int_t group, Float_t value);
          Bool_t   UpdateVectorPRF(Int_t group, Float_t x, Float_t y);
          Bool_t   UpdateVectorPH(Int_t group, Int_t time, Float_t value);
          Bool_t   UpdateVectorT0(Int_t group, Int_t time);
  TGraphErrors    *ConvertVectorPHisto(AliTRDPInfo *pInfo, const Char_t *name) const;
          TH1F    *ConvertVectorCTHisto(AliTRDCTInfo *cTInfo, const Char_t *name) const;
          TTree   *ConvertVectorCTTreeHisto(TObjArray *vVectorCT, TObjArray *pPlaCT, const Char_t *name, const Char_t *nametitle) const;
          TTree   *ConvertVectorPTreeHisto(TObjArray *vVectorP, TObjArray *pPlaP, const Char_t *name, const Char_t *nametitle) const;
  TObjArray       *ConvertTreeVector(TTree *tree) const ;
          Bool_t   MergeVectorCT(TObjArray *vVectorCT2, TObjArray *pPlaCT2);
          Bool_t   MergeVectorP(TObjArray *vVectorP2, TObjArray *pPlaP2, Int_t i);

  // Fit methods
          void     FitBisCH(TH1 *projch, Int_t idect);
          void     FitCH(TH1 *projch, Int_t idect);
          void     FitPH(TH1 *projPH, Int_t idect);
          void     FitPRF(TH1 *projPRF, Int_t idect);
          void     FitPente(TH1 *projPH, Int_t idect);
          TH1I    *ReBin(TH1I *hist) const;
          TH1F    *ReBin(TH1F *hist) const;
  
  // Clear
          void     ClearHistos();
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
  static  AliTRDCalibra *fgInstance;        // Instance
  static  Bool_t   fgTerminated;            // If terminated
    
  ClassDef(AliTRDCalibra, 1)   // TRD Calibration class

};
  
#endif


