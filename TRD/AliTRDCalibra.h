#ifndef ALITRDCALIBRA_H
#define ALITRDCALIBRA_H

#include <Riostream.h>
#include <vector>

#include <TProfile2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLine.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TObject.h>
#include <TLatex.h>
#include <TChain.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TObject.h>
#include <TStopwatch.h>
#include <TList.h>
#include <TCollection.h>

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliRawReader.h"
#include "../RAW/AliRawReaderFile.h"

#include "AliTRD.h"
#include "AliTRDgeometry.h"
#include "AliTRDSimParam.h"
#include "AliTRDCommonParam.h"
#include "AliTRDcalibDB.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDdigit.h"
#include "AliTRDcluster.h"
#include "AliTRDtrack.h"
#include "AliTRDpadPlane.h"
#include "AliTRDrawData.h"
#include "AliTRDmcmTracklet.h"
#include "../TRD/Cal/AliTRDCalDet.h"
#include "../TRD/Cal/AliTRDCalDet.h"
#include "../TRD/Cal/AliTRDCalPad.h"
#include "../TRD/Cal/AliTRDCalROC.h"

class TPInfo;
class TCTInfo;
class TFitCHInfo;

class AliTRDCalibra : public TObject {
 public: 

   // Instance
  static AliTRDCalibra* Instance();
  static void Terminate();
  static void Destroy();

  AliTRDCalibra(const AliTRDCalibra &c);
  AliTRDCalibra &operator=(const AliTRDCalibra &) { return *this; }


  // Functions fit online
  Bool_t FitCHOnline(TH2I *ch);
  Bool_t FitCHOnline();
  Bool_t FitCHOnline(TTree *tree);
  Bool_t FitPHOnline(TProfile2D *PH);
  Bool_t FitPHOnline();
  Bool_t FitPHOnline(TTree *tree);
  Bool_t FitPRFOnline(TProfile2D *Prf);
  Bool_t FitPRFOnline();
  Bool_t FitPRFOnline(TTree *tree);

 
  // Functions for initialising the AliTRDCalibra in the code
  Bool_t Init2Dhistos();
  
  
  // Functions for filling the histos in the code
  Bool_t Resettrack();
  Bool_t UpdateHistograms(AliTRDcluster *cl, AliTRDtrack *t);
  Bool_t UpdateHistogramcm(AliTRDmcmTracklet *fTrk);
 

  // Is Pad on
  Bool_t IsPadOn(Int_t detector, Int_t col, Int_t row);

  // Functions for plotting the 2D
  void   Plot2d();

  // Functions for writting the 2D
  Bool_t Write2d();
 

  // Function fill 2D for the moment out of the code
  Bool_t Create2DDiSimOnline(Int_t iev1, Int_t iev2);
  Bool_t Create2DRaDaOnline(Int_t iev1, Int_t iev2);
  
  
  // Pad Calibration
  // From Nnz and Nnrphi calculates fNfragz and fNfragrphi
  Bool_t ModePadFragmentation(Int_t iPlane,Int_t iChamb, Int_t iSect, Int_t i);
  // From Nz and Nrphi calculates Nnz and Nnrphi
  void   ModePadCalibration(Int_t iChamb, Int_t i);


  // Fill the database
  TObject*      CreatePadObjectTree(TTree *tree);
  TObject*      CreatePadObjectTree(TTree *tree, Int_t i, AliTRDCalDet *detobject);
  AliTRDCalDet* CreateDetObjectTree(TTree *tree, Int_t i);

  // Correct the error
  TH1F*         CorrectTheError(TGraphErrors *hist);
  TGraphErrors* AddProfiles(TGraphErrors *hist1, TGraphErrors *hist2);

  // Add two trees
  TTree*        Sum2Trees(const char* filename1, const char* filename2, const char *variablecali);
  

  // Set of Get the variables*******************************************************************************

  // Choice to fill or not the 2D
  void SetMItracking()                                           { fMItracking      = kTRUE; }
  void Setmcmtracking()                                          { fmcmtracking     = kTRUE; }
  void Setmcmcorrectangle()                                      { fmcmcorrectangle = kTRUE; }
  void SetPH2dOn()                                               { fPH2dOn          = kTRUE; }
  void SetCH2dOn()                                               { fCH2dOn          = kTRUE; }
  void SetPRF2dOn()                                              { fPRF2dOn         = kTRUE; }
  void SetHisto2d()                                              { fHisto2d         = kTRUE; }
  void SetVector2d()                                             { fVector2d        = kTRUE; }
  
  Bool_t GetMItracking() const                                   { return fMItracking;       }
  Bool_t Getmcmtracking() const                                  { return fmcmtracking;      }
  Bool_t Getmcmcorrectangle() const                              { return fmcmcorrectangle;  }
  Bool_t GetPH2dOn() const                                       { return fPH2dOn;           }
  Bool_t GetCH2dOn() const                                       { return fCH2dOn;           }
  Bool_t GetPRF2dOn() const                                      { return fPRF2dOn;          }
  Bool_t GetHisto2d() const                                      { return fHisto2d;          }
  Bool_t GetVector2d() const                                     { return fVector2d;         }
  TH2I       *GetCH2d() const                                    { return fCH2d;             }
  TProfile2D *GetPH2d() const                                    { return fPH2d;             }
  TProfile2D *GetPRF2d() const                                   { return fPRF2d;            }
  
  
  
  // How to fill the 2D
  void SetRelativeScaleAuto()                                    { fRelativeScaleAuto    = kTRUE;                      }
  void SetRelativeScale(Float_t RelativeScale)                   { if (RelativeScale > 0.0) {
                                                                     fRelativeScale = RelativeScale;
                                                                   } 
                                                                   else {
                                                                     AliInfo("RelativeScale must be strict positif!");
                                                                   }                                                   }
  void SetThresholdDigit(Int_t digitthreshold)                   { fThresholddigit       = digitthreshold;             }
  void SetThresholdClusterPRF1(Float_t ThresholdClusterPRF1)     { fThresholdClusterPRF1 = ThresholdClusterPRF1;       }
  void SetThresholdClusterPRF2(Float_t ThresholdClusterPRF2)     { fThresholdClusterPRF2 = ThresholdClusterPRF2;       }
  void SetCenterOfflineCluster()                                 { fCenterOfflineCluster = kTRUE;                      }
  void SetTraMaxPad()                                            { fTraMaxPad            = kTRUE;                      }
  void SetNz(Int_t i, Short_t Nz)                                { if ((Nz >= 0) && (Nz < 5)) {
                                                                     fNz[i] = Nz; 
                                                                   }
                                                                   else { 
                                                                     AliInfo("You have to choose between 0 and 4");
                                                                   }                                                   }
  void SetNrphi(Int_t i, Short_t Nrphi)                          { if ((Nrphi >= 0) && (Nrphi < 7)) {
                                                                     fNrphi[i] = Nrphi; 
                                                                   }
                                                                   else {
                                                                     AliInfo("You have to choose between 0 and 6");
                                                                   }                                                   }
  void SetProcent(Float_t procent)                               { fprocent              = procent;                    }
  void SetDifference(Short_t difference)                         { fdifference           = difference;                 }
  void SetNumberClusters(Short_t NumberClusters)                 { fNumberClusters       = NumberClusters;             }
  void SetNumberBinCharge(Short_t NumberBinCharge)               { fNumberBinCharge      = NumberBinCharge;            }
  void SetNumberBinPRF(Short_t NumberBinPRF)                     { fNumberBinPRF         = NumberBinPRF;               }
  
  Float_t GetRelativeScale() const                               { return fRelativeScale;        }
  Bool_t  GetRelativeScaleAuto() const                           { return fRelativeScaleAuto;    }
  Int_t   GetThresholdDigit() const                              { return fThresholddigit;       }
  Float_t GetThresholdClusterPRF1() const                        { return fThresholdClusterPRF1; }
  Float_t GetThresholdClusterPRF2() const                        { return fThresholdClusterPRF2; }
  Bool_t  GetTraMaxPad()const                                    { return fTraMaxPad;            }
  Short_t GetNz(Int_t i) const                                   { return fNz[i];                }
  Short_t GetNrphi(Int_t i) const                                { return fNrphi[i];             }
  Float_t GetProcent() const                                     { return fprocent;              }
  Short_t GetDifference() const                                  { return fdifference;           }
  Short_t GetNumberClusters() const                              { return fNumberClusters;       }
  Short_t GetNumberBinCharge() const                             { return fNumberBinCharge;      }
  Short_t GetNumberBinPRF() const                                { return fNumberBinPRF;         }
  
  
  // Write
  void SetWriteCoef(Int_t i)                                     { fWriteCoef[i]  = kTRUE;         }
  void SetWriteNameCoef(TString WriteNameCoef)                   { fWriteNameCoef = WriteNameCoef; }
  void SetWrite(Int_t i)                                         { fWrite[i]      = kTRUE;         }
  void SetWriteName(TString WriteName)                           { fWriteName     = WriteName;     }
  
  Bool_t  GetWriteCoef(Int_t i) const                            { return fWriteCoef[i];           }
  TString GetWriteNameCoef() const                               { return fWriteNameCoef;          }
  Bool_t  GetWrite(Int_t i) const                                { return fWrite[i];               }
  TString GetWriteName() const                                   { return fWriteName;              }
  
  
  // Fit
  void SetFitPHOn()                                              { fFitPHOn        = kTRUE;}
  void SetPeriodeFitPH(Int_t PeriodeFitPH)                       { if (PeriodeFitPH   > 0) {
                                                                     fFitPHPeriode   = PeriodeFitPH; 
                                                                   }
                                                                   else { 
                                                                     AliInfo("PeriodeFitPH must be higher than 0!");
                                                                   }                                                    }
  void SetBeginFitCharge(Float_t BeginFitCharge)                 { if (BeginFitCharge > 0) {
                                                                     fBeginFitCharge = BeginFitCharge; 
                                                                   }
                                                                   else {
                                                                     AliInfo("BeginFitCharge must be strict positif!");
                                                                   }                                                    }
  void SetT0Shift(Float_t T0Shift)                               { if (T0Shift        > 0) {
                                                                     fT0Shift        = T0Shift; 
                                                                   } 
                                                                   else {
                                                                     AliInfo("T0Shift must be strict positif!");
                                                                   }                                                    }
  void SetRangeFitPRF(Float_t RangeFitPRF)                       { if ((RangeFitPRF   >   0) && 
                                                                       (RangeFitPRF  <= 1.0)) {
                                                                     fRangeFitPRF    = RangeFitPRF;
                                                                   } 
                                                                   else {
                                                                     AliInfo("RangeFitPRF must be between 0 and 1.0");
								   }                                                    }                
  void SetMeanChargeOn()                                         { fMeanChargeOn   = kTRUE;}
  void SetFitChargeBisOn()                                       { fFitChargeBisOn = kTRUE;}
  void SetMinEntries(Int_t MinEntries)                           { fMinEntries     = MinEntries;}
  
  Bool_t  GetFitPHOn() const                                     { return fFitPHOn;        }
  Int_t   GetPeriodeFitPH() const                                { return fFitPHPeriode;   }
  Float_t GetBeginFitCharge() const                              { return fBeginFitCharge; }
  Float_t GetT0Shift() const                                     { return fT0Shift;        }
  Float_t GetRangeFitPRF() const                                 { return fRangeFitPRF;    }
  Bool_t  GetMeanChargeOn() const                                { return fMeanChargeOn;   }
  Bool_t  GetFitChargeBisOn() const                              { return fFitChargeBisOn; }
  Int_t   GetMinEntries() const                                  { return fMinEntries;     }
  
  // Debug  
  void SetDebug(Short_t Debug)                                   { fDebug   = Debug;   }
  void SetDet(Int_t iPlane, Int_t iChamb, Int_t iSect)           { fDet[0]  = iPlane; 
                                                                   fDet[1]  = iChamb; 
                                                                   fDet[2]  = iSect;   }
  void SetFitVoir(Int_t FitVoir)                                 { fFitVoir = FitVoir; }
  
  Short_t GetDebug() const                                       { return fDebug;      }
  Int_t   GetDet(Int_t i) const                                  { return fDet[i];     }
  Int_t   GetFitVoir() const                                     { return fFitVoir;    }
  
  // Internal variables to be sure!*****************************************
  
  // Pad calibration
  Short_t GetNnz(Int_t i) const                                  { return fNnz[i];       }
  Short_t GetNnrphi(Int_t i) const                               { return fNnrphi[i];    }
  Short_t GetNfragz(Int_t i) const                               { return fNfragz[i];    }
  Short_t GetNfragrphi(Int_t i) const                            { return fNfragrphi[i]; }
  
  
  void    SetRebin(Short_t Rebin)                                { if (Rebin > 0) {
                                                                     fRebin = Rebin; 
                                                                     AliInfo("You have to be sure that fRebin divides fNumberBinCharge used!");
                                                                   } 
                                                                   else {
                                                                     AliInfo("You have to choose a positiv value!");
								   }                     }
  
  Short_t GetRebin() const                                       { return fRebin;        }

  // Getter for the coefficient trees***********************************

  TTree   *GetPRF() const                                        { return fPRF;          }
  TTree   *GetGain() const                                       { return fGain;         }
  TTree   *GetT0() const                                         { return fT0;           }
  TTree   *GetVdrift() const                                     { return fVdrift;       }
  
  protected:

  // Instance of this class and so on
  static   AliTRDCalibra* fgInstance;   // Instance
  static   Bool_t fgTerminated;

  
  // Choice to fill or not the 2D
  Bool_t   fMItracking;                 // Chose to fill the 2D histos or vectors during the offline MI tracking
  Bool_t   fmcmtracking;                // Chose to fill the 2D histos or vectors during the tracking with tracklets
  Bool_t   fmcmcorrectangle;            // Apply correction due to the mcmtrackletangle in the z direction (only) assuming  from vertex
  Bool_t   fCH2dOn;                     // Chose to fill the 2D histos or vectors for the relative gain calibration 
  Bool_t   fPH2dOn;                     // Chose to fill the 2D histos or vectors for the drift velocity and T0
  Bool_t   fPRF2dOn;                    // Chose to fill the 2D histos or vectors for the pad response function calibration
  Bool_t   fHisto2d;                    // Chose to fill the 2D histos
  Bool_t   fVector2d;                   // Chose to fill vectors
  

  // How to fill the 2D
  Float_t  fRelativeScale;              // Scale of the deposited charge
  Int_t    fCountRelativeScale;         // fCountRelativeScale first data used for the scaling
  Bool_t   fRelativeScaleAuto;          // Scaling with the first fCountRelativeScale objects
  Int_t    fThresholddigit;             // Threshold on RawData
  Float_t  fThresholdClusterPRF1;       // Threshold on cluster pad signals for PRF peripherique
  Float_t  fThresholdClusterPRF2;       // Threshold on cluster pad signals for PRF peripherique
  Bool_t   fCenterOfflineCluster;       // Choose to use the offline determination of the center of the cluster
  Bool_t   fTraMaxPad;                  // Take the Max Pad for the gain calibration and PH
  Short_t  fNz[3];                      // Mode of calibration 
  Short_t  fNrphi[3];                   // Mode of calibration 
  Int_t    fNtotal[3];                  // Total number of Xbins

  
  // Write
  Bool_t   fWriteCoef[3];               // Do you want to write the result in a file?
  TString  fWriteNameCoef;              // Where the coef Det are written
  Bool_t   fWrite[3];                   // Do you want to write the 2D histo or vectors converted in a tree
  TString  fWriteName;                  // Where the 2D or trees are written
  
  
  // Fit
  Bool_t   fFitPHOn;                    // The fit PH On
  Int_t    fFitPHPeriode;               // Periode of the fit PH
  Float_t  fBeginFitCharge;             // The fit begins at mean/fBeginFitCharge for the gain calibration
  Float_t  fRangeFitPRF;                // The fit range for the PRF is -fRangeFitPRF +fRangeFitPRF
  Bool_t   fMeanChargeOn;               // Mean Charge on
  Bool_t   fFitChargeBisOn;             // For an other fit function (convolution and not sum, more time consuming)
  Float_t  fT0Shift;                    // T0 Shift with the actual method
  
  
  // Debug Mode
  Short_t  fDebug;                      // For debugging 0 rien, 1 errors, 2 one fit alone, 3 one detector, 4 one detector with errors
  Int_t    fDet[3];                     // Detector  visualised (plane,chamb,sect) si debugging == 3 or 4
  Int_t    fFitVoir;                    // Fit visualised si debugging == 2
  
  
  // Internal variables*************************************************
  
  // Storage of coef
  TTree   *fPRF;                        // Tree of the sigma of PRD
  TTree   *fGain;                       // Tree of the gain factor
  TTree   *fT0;                         // Tree of the time0
  TTree   *fVdrift;                     // Tree of the drift velocity

  // "Pointer" of the branch of the tree
  Int_t    fVdriftDetector;             // Branch of Vdrift
  Float_t *fVdriftPad;                  // Branch of Vdrift
  Int_t    fT0Detector;                 // Branch of t0
  Float_t *fT0Pad;                      // Branch of t0
  Int_t    fPRFDetector;                // Branch of PRF
  Float_t *fPRFPad;                     // Branch of PRF
  Float_t *fcoefCH;                     // Branch relative gain
  
  // Fill the 2D histos in the offline tracking
  Bool_t   fDetectorAliTRDtrack;        // Change of track
  Int_t    fChamberAliTRDtrack;         // Change of chamber
  Int_t    fDetectorprevioustrack;      // Change of detector
  Bool_t   fGoodTrack;                  // If goes through a kaputt pad
  Float_t *famptotal;                   // Energy deposited in the calibration group by the track
  Short_t *fPHplace;                    // Calibration group of PH
  Float_t *fPHvalue;                    // PH
  Short_t  fNumberClusters;             // Minimum number of clusters in the tracklets
  Float_t  fprocent;                    // Limit to take the info of the most important calibration group if the track goes through 2 groups (CH)
  Short_t  fdifference;                 // Limit to take the info of the most important calibration group if the track goes through 2 groups (CH)
  Int_t    fNumbertrack;                // How many tracks could be used (Debug for the moment)
  Int_t    fNumberusedch[2];            // How many tracks have been really used for the gain (0, strict; 1 with fprocent)
  Int_t    fNumberusedph[2];            // How many tracks have been really used for the drift velocity (0, strict; 1 with fdifference)


  // For debugging 

  // Histograms to store the coef
  TH1F    *fCoefCharge[4];              // Replica des 2D but in coefs resulting from the fit for the gain
  TH1F    *fCoefVdrift[3];              // Replica des 2D but in coefs resulting from the fit for the drift velocity
  TH1F    *fCoefPRF[2];                 // Replica des 2D but in coefs resulting from the fit for the pad response function
  TH1F    *fCoefT0[3];                  // Replica des 2D but in coefs resulting from the fit for time 0
  TH1F    *fDeltaCharge[3];             // Replica des 2D but in errors for each detector resulting from the fit for the gain
  TH1F    *fDeltaVdrift[2];             // Replica des 2D but in errors for each detector resulting from the fit for the drift velocity
  TH1F    *fDeltaT0[2];                 // Replica des 2D but in errors for each detector resulting from the fit for time 0
  TH1F    *fDeltaPRF;                   // Replica des 2D but in errors for each detector resulting from the fit for the pad response function
  TH1I    *fErrorCharge[3];             // Replica des 2D but in errors resulting from the fit for the gain
  TH1I    *fErrorVdrift[2];             // Replica des 2D but in errors resulting from the fit for the drift velocity
  TH1I    *fErrorT0[2];                 // Replica des 2D but in errors resulting from the fit for time 0
  TH1I    *fErrorPRF;                   // Replica des 2D but in errors resulting from the fit for the pad response function
  TH2F    *fCoefChargeDB[3];            // Visualisation of the coef of the detecteur fDet for the gain
  TH2F    *fCoefVdriftDB[2];            // Visualisation of the coef of the detecteur fDet for the drift velocity
  TH2F    *fCoefT0DB[2];                // Visualisation of the coef of the detecteur fDet for time 0
  TH2F    *fCoefPRFDB;                  // Visualisation of the coef of the detecteur fDet for the pad response function
  

  // Variables in the loop for the coef or more general
  Float_t  fChargeCoef[4];              // 3 database value, 0 fit, 1 mean, 2 fit time consuming   
  Float_t  fVdriftCoef[3];              // 2 database value, 1 slope method, 0 fit
  Float_t  fPRFCoef[2];                 // 1 database, 0 fit 
  Float_t  fT0Coef[3];                  // 3 database, 1 slope method, 0 fit
  Float_t  fPhd[3];                     // Begin AR and DR
  Int_t    fTimeMax;                    // Number of time bins
  Float_t  fSf;                         // Sampling frequence
  Int_t    fdect1[3];                   // First calibration group that will be called to be maybe fitted
  Int_t    fdect2[3];                   // Last calibration group that will be called to be maybe fitted
  Double_t fScalefitfactor;             // Scale factor of the fit results for the gain
  Int_t    fMinEntries;                 // Min Entries to fit the histo
  Int_t    fEntriesCurrent;             // Entries in the current histo
  Int_t    fcountdet[3];                // Current detector
  Int_t    fcount[3];                   // When the next detector comes
  Float_t  fl3P0;                       // ????
  Float_t  fl3P2;                       // ????
  Float_t  fg3P2;                       // ????

  
  // Pad Calibration
  Short_t  fNnz[3];                     // Number of pad rows in a group
  Short_t  fNnrphi[3];                  // Number of pad cols in a group
  Short_t  fNfragz[3];                  // Number of  pad row group
  Short_t  fNfragrphi[3];               // Number of pad col group
  Short_t  frowmin[3];                  // Limits of the group in pad row
  Short_t  frowmax[3];                  // Limits of the group in pad row
  Short_t  fcolmin[3];                  // Limits of the group in pad col
  Short_t  fcolmax[3];                  // Limits of the group in pad col
  Int_t    fXbins[3];                   // First Xbins of the detector
  Short_t  fdetChamb0[3];               // Number of XBins for chamber != 2
  Short_t  fdetChamb2[3];               // Number of Xbins fir chamber 2
  
  
  // Methode  Alexandru store info
  class TCTInfo : public TObject {
  public: 

    TCTInfo()
   :TObject()
   ,fentries(0x0)                            { }

    TCTInfo(const TCTInfo &i)
   :TObject(i)
   ,fentries(0x0)                            { }

    TCTInfo &operator=(const TCTInfo&)       { return *this; }

    UShort_t *fentries;       // ????

  };

  class TFitCHInfo : public TObject {
  public:

    TFitCHInfo()
   :TObject()
   ,fcoef(0x0)
   ,fDetector(-1)                            { }

    TFitCHInfo(const TFitCHInfo &i) 
   :TObject(i)
   ,fcoef(0x0)
   ,fDetector(-1)                            { }

    TFitCHInfo &operator=(const TFitCHInfo&) { return *this; }

    Float_t *fcoef;           // ????
    Int_t    fDetector;       // ????

  };
  
  class TPInfo : public TObject {
  public:

    TPInfo()
   :TObject()
   ,fsum(0x0) 
   ,fsumsquare(0x0)
   ,fentries(0x0)                            { }

    TPInfo(const TPInfo &i)
   :TObject(i)
   ,fsum(0x0) 
   ,fsumsquare(0x0)
   ,fentries(0x0)                            { }

    TPInfo &operator=(const TPInfo&)         { return *this; }

    Float_t  *fsum;           // ????
    Float_t  *fsumsquare;     // ????
    UShort_t *fentries;       // ????

  }; 
  
  
  // PH
  // fTimeMax will define the size of fcharge
  std::vector<TPInfo *>     fVectorPH;    // Vectors to fill
  std::vector<Int_t>        fPlaPH;       // Vectors to fill
  // CH
  Short_t fNumberBinCharge;               // Number of bins for the gain factor
  std::vector<TCTInfo *>    fVectorCH;    // Vectors to fill
  std::vector<Int_t>        fPlaCH;       // Vectors to fill
  // FitCH
  std::vector<TFitCHInfo *> fVectorFitCH; // Vectors to fit
  // PRF
  Short_t fNumberBinPRF;                  // Number of bin for the PRF
  std::vector<TPInfo *>     fVectorPRF;   // Vectors to fill
  std::vector<Int_t>        fPlaPRF;      // Vectors to fill
 
  
  // Histograms to store the info from the digits, from the tracklets or from the tracks
  TProfile2D *fPH2d;                      // 2D average pulse height
  TProfile2D *fPRF2d;                     // 2D PRF
  TH2I       *fCH2d;                      // 2D deposited charge 
 
  Short_t fRebin;                         // If you want to rebin the histo for the gain calibration
  

  // A lot of internal functions......*************************************************

  // Init AliTRDCalibra
  void Init();


  // Create the 2D histo to be filled Online
  void CreateCH2d(Int_t Nn);
  void CreatePH2d(Int_t Nn);
  void CreatePRF2d(Int_t Nn);
  
  
  // Fill the 2D from mcmtracklet from the TRD.Track.root file and not yet in the code
  void   FillTheInfoOfTheTrackPH();
  void   FillTheInfoOfTheTrackCH();
  void   Resetfvariables();
  Bool_t LocalisationdetectorXbins(Int_t detector);
  
 
  // Plot the 2D
  void PlotCH2d();
  void PlotPH2d();
  void PlotPRF2d();
    

  // Fit**************************************************
  
  // Create histos if fDebug == 1 or fDebug >=3
  void CreateFitHistoPHDB(Int_t rowMax, Int_t colMax);
  void CreateFitHistoT0DB(Int_t rowMax, Int_t colMax);
  void CreateFitHistoCHDB(Int_t rowMax, Int_t colMax);
  void CreateFitHistoPRFDB(Int_t rowMax, Int_t colMax);
  void CreateFitHistoCH(Int_t Nbins, Double_t Low, Double_t High);
  void CreateFitHistoPH(Int_t Nbins, Double_t Low, Double_t High);
  void CreateFitHistoT0(Int_t Nbins, Double_t Low, Double_t High);
  void CreateFitHistoPRF(Int_t Nbins, Double_t Low, Double_t High);
  
  // CHFit functions
  Bool_t FillVectorFitCH(Int_t countdet);
  Bool_t InitFit(Int_t Nbins, Double_t lowedge, Double_t upedge, Int_t i);
  void   Initfcountdetandfcount(Int_t i);
  void   Updatefcountdetandfcount(Int_t idect, Int_t i);
  void   Reconstructfitrowminrowmax(Int_t idect, Int_t i);
  Bool_t NotEnoughStatistic(Int_t idect, Int_t i);
  Bool_t FillInfosFit(Int_t idect, Int_t i);
  Bool_t WriteFitInfos(Int_t i);
  void   NormierungCharge();


  // Fill histos Errors from the delta histos
  void ErrorPH();
  void ErrorT0();
  void ErrorCH();
  void ErrorPRF();
  

  // Fill histos DB from the Coef histos 
  void FillCoefChargeDB();
  void FillCoefVdriftDB();
  void FillCoefT0DB();
  void FillCoefPRFDB();

  
  // Plot histos CoefPRF Coef....
  void PlotPH();
  void PlotT0();
  void PlotCH();
  void PlotPRF();
  

  // Plot histos DB
  void PlotPHDB();
  void PlotT0DB();
  void PlotCHDB();
  void PlotPRFDB();


  // Write the Coef, delta and error histos
  void WritePH(TFile *fout);
  void WriteT0(TFile *fout);
  void WriteCH(TFile *fout);
  void WritePRF(TFile *fout);
  
  
  //Write the DB histos
  void WritePHDB(TFile *fout);
  void WriteT0DB(TFile *fout);
  void WriteCHDB(TFile *fout);
  void WritePRFDB(TFile *fout);
  
  
  // Calculate the mean coefs from the database
  Bool_t CalculVdriftCoefMean(Int_t Dect, Int_t idect);
  Bool_t CalculChargeCoefMean(Int_t Dect, Int_t idect, Bool_t vrai);
  Bool_t CalculPRFCoefMean(Int_t Dect, Int_t idect);
  Bool_t CalculT0CoefMean(Int_t Dect, Int_t idect);
  
  
  // Pad group calibration mode
  void ReconstructionRowPadGroup(Int_t idect, Int_t i);
  void CalculXBins(Int_t idect, Int_t i);  


  // Convertion vector, tree, histos....
  Int_t  SearchInVector(Int_t group, Int_t i);
  Int_t  SearchInTreeVector(std::vector<Int_t> vectorplace, Int_t group);
  Int_t  SearchBin(Float_t value, Int_t i);
  Bool_t UpdateVectorCH(Int_t group, Float_t value);
  Bool_t UpdateVectorPRF(Int_t group, Float_t x, Float_t y);
  Bool_t UpdateVectorPH(Int_t group, Int_t time, Float_t value);
  Bool_t UpdateVectorT0(Int_t group, Int_t time);
  TGraphErrors      *ConvertVectorPHisto(TPInfo *fPInfo, const char* name);
  TH1F              *ConvertVectorCTHisto(TCTInfo *fCTInfo, const char* name);
  TTree             *ConvertVectorCTTreeHisto(std::vector<TCTInfo *> VectorCT, std::vector<Int_t> PlaCT, const char* name);
  TTree             *ConvertVectorPTreeHisto(std::vector<TPInfo *> VectorP, std::vector<Int_t> PlaP, const char* name);
  std::vector<Int_t> ConvertTreeVector(TTree *tree);
  Bool_t MergeVectorCT(std::vector<TCTInfo* > VectorCT2, std::vector<Int_t> PlaCT2);
  Bool_t MergeVectorP(std::vector<TPInfo* > VectorP2, std::vector<Int_t> PlaP2, Int_t i);

  
  // Fit methods
  void  FitBisCH(TH1* projch, Int_t idect);
  void  FitCH(TH1* projch, Int_t idect);
  void  FitPH(TH1* projPH, Int_t idect);
  void  FitPRF(TH1* projPRF, Int_t idect);
  void  FitPente(TH1* projPH, Int_t idect);
  TH1I *ReBin(TH1I* hist);
  TH1F *ReBin(TH1F* hist);
  
  
  // Clear
  void ClearHistos();
  void ClearTree();
  void ClearFile();


  // Some basic geometry function
  virtual Int_t GetPlane(Int_t d) const;
  virtual Int_t GetChamber(Int_t d) const;
  virtual Int_t GetSector(Int_t d) const;
  
  // Init, Fill and Reset the variables to default value tree Gain, PRF, Vdrift and T0
  void InittreePH();
  void FilltreeVdrift(Int_t countdet);
  void InittreeT0();
  void FilltreeT0(Int_t countdet);
  void InittreePRF();
  void FilltreePRF(Int_t countdet);
  void ConvertVectorFitCHTree();
  
  
  private:
  
  static Double_t PH(Double_t* x, Double_t* par);
  static Double_t AsymmGauss(Double_t* x, Double_t* par);
  static Double_t funcLandauGaus(Double_t* x, Double_t* par);
  static Double_t langaufun(Double_t *x, Double_t *par);
  TF1 *langaufit(TH1 *his, Double_t *fitrange, Double_t *startvalues
               , Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams
               , Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF);
  Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM); 
  
  // This is a singleton, contructor is private!
  AliTRDCalibra();
  virtual ~AliTRDCalibra();
  
  ClassDef(AliTRDCalibra, 1)   // TRD Calibration class

};
  
#endif


