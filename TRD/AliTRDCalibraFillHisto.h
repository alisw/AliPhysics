#ifndef ALITRDCALIBRAFILLHISTO_H
#define ALITRDCALIBRAFILLHISTO_H
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
#ifndef ROOT_TLinearFitter
#  include <TLinearFitter.h> 
#endif
#ifndef ROOT_TProfile2D
#  include <TProfile2D.h> 
#endif
#ifndef ROOT_TH2I
#  include <TH2I.h> 
#endif

class TProfile2D;
class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH2;
class TLinearFitter;
class TTreeSRedirector;

class AliLog;
class AliRawReader;

class AliTRDCalibraMode;
class AliTRDCalibraVector;
class AliTRDCalibraVdriftLinearFit;
class AliTRDRawStreamV2;
class AliTRDcluster;
class AliTRDtrack;
class AliTRDmcmTracklet;
class AliTRDgeometry;
class AliTRDCalDet;
class AliTRDCalROC;


struct eventHeaderStruct;

class AliTRDCalibraFillHisto : public TObject {

 public: 

  // Instance
  static AliTRDCalibraFillHisto *Instance();
  static void Terminate();
  static void Destroy();

  AliTRDCalibraFillHisto(const AliTRDCalibraFillHisto &c);
  AliTRDCalibraFillHisto &operator=(const AliTRDCalibraFillHisto &) { return *this; }

  // Functions for initialising the AliTRDCalibraFillHisto in the code
          Bool_t   Init2Dhistos();
	  Bool_t   Init2Dhistostrack();

  // Functions for filling the histos in the code
          Bool_t   ResetTrack();
	  Bool_t   UpdateHistograms(AliTRDtrack *t);
          Bool_t   UpdateHistograms(AliTRDcluster *cl, AliTRDtrack *t);
          Bool_t   UpdateHistogramcm(AliTRDmcmTracklet *trk);
 
  // Process events DAQ
	  Int_t   ProcessEventDAQ(AliTRDRawStreamV2 *rawStream, Bool_t nocheck = kFALSE);
	  Int_t   ProcessEventDAQ(AliRawReader *rawReader, Bool_t nocheck = kFALSE);
	  Int_t   ProcessEventDAQ(eventHeaderStruct *event, Bool_t nocheck = kFALSE);

	  Bool_t   UpdateDAQ(Int_t det, Int_t /*row*/, Int_t /*col*/, Int_t timebin, Int_t signal, Int_t nbtimebins);

  // Is Pad on
          Bool_t   IsPadOn(Int_t detector, Int_t col, Int_t row) const;

  // Functions for write
	  void     Write2d(const Char_t *filename = "TRD.calibration.root", Bool_t append = kFALSE);

  //For the statistics
	  Double_t *StatH(TH2 *ch, Int_t i);
	  Double_t *GetMeanMedianRMSNumberCH();
	  Double_t *GetMeanMedianRMSNumberLinearFitter() const;
       
     	 
  //
  // Set of Get the variables
  //

  // Choice to fill or not the 2D
          void     SetMITracking(Bool_t mitracking = kTRUE)                  { fMITracking      = mitracking;        }
          void     SetMcmTracking(Bool_t mcmtracking = kTRUE)                { fMcmTracking     = mcmtracking;       }
          void     SetMcmCorrectAngle(Bool_t mcmcorrectangle = kTRUE)        { fMcmCorrectAngle = mcmcorrectangle;   }
          void     SetPH2dOn(Bool_t ph2don = kTRUE)                          { fPH2dOn          = ph2don;            }
          void     SetCH2dOn(Bool_t ch2don = kTRUE)                          { fCH2dOn          = ch2don;            }
          void     SetPRF2dOn(Bool_t prf2don = kTRUE)                        { fPRF2dOn         = prf2don;           }
          void     SetHisto2d(Bool_t histo2d = kTRUE)                        { fHisto2d         = histo2d;           }
          void     SetVector2d(Bool_t vector2d = kTRUE)                      { fVector2d        = vector2d;          }
	  void     SetLinearFitterOn(Bool_t linearfitteron = kTRUE)          { fLinearFitterOn      = linearfitteron;}
	  void     SetLinearFitterDebugOn(Bool_t debug = kTRUE)              { fLinearFitterDebugOn = debug;         }
	 	  
  
          Bool_t   GetMITracking() const                                     { return fMITracking;             }
          Bool_t   GetMcmTracking() const                                    { return fMcmTracking;            }
          Bool_t   GetMcmCorrectAngle() const                                { return fMcmCorrectAngle;        }
          Bool_t   GetPH2dOn() const                                         { return fPH2dOn;                 }
          Bool_t   GetCH2dOn() const                                         { return fCH2dOn;                 }
          Bool_t   GetPRF2dOn() const                                        { return fPRF2dOn;                }
          Bool_t   GetHisto2d() const                                        { return fHisto2d;                }
          Bool_t   GetVector2d() const                                       { return fVector2d;               }
          Bool_t   GetLinearFitterOn() const                                 { return fLinearFitterOn;         }
	  Bool_t   GetLinearFitterDebugOn() const                            { return fLinearFitterDebugOn;    }
  
  TH2I            *GetCH2d() const                                           { return fCH2d;                   }
  TProfile2D      *GetPH2d(Int_t nbtimebin=24, Float_t samplefrequency= 10.0);
  TProfile2D      *GetPRF2d() const                                          { return fPRF2d;                  } 
  TObjArray        GetLinearFitterArray() const                              { return fLinearFitterArray;      }
  TLinearFitter   *GetLinearFitter(Int_t detector, Bool_t force=kFALSE);
  AliTRDCalibraVdriftLinearFit *GetVdriftLinearFit() const                   { return fLinearVdriftFit;        }
  
 
  // How to fill the 2D
          void     SetRelativeScale(Float_t relativeScale);                      
          void     SetThresholdClusterPRF2(Float_t thresholdClusterPRF2)     { fThresholdClusterPRF2 = thresholdClusterPRF2; }
	  void     SetNz(Int_t i, Short_t nz);
          void     SetNrphi(Int_t i, Short_t nrphi);
          void     SetProcent(Float_t procent)                               { fProcent              = procent;              }
          void     SetDifference(Short_t difference)                         { fDifference           = difference;           }
          void     SetNumberClusters(Short_t numberClusters)                 { fNumberClusters       = numberClusters;       }
          void     SetNumberBinCharge(Short_t numberBinCharge)               { fNumberBinCharge      = numberBinCharge;      }
          void     SetNumberBinPRF(Short_t numberBinPRF)                     { fNumberBinPRF         = numberBinPRF;         }
	  void     SetNumberGroupsPRF(Short_t numberGroupsPRF);
  
          Float_t  GetRelativeScale() const                                  { return fRelativeScale;          }
          Float_t  GetThresholdClusterPRF2() const                           { return fThresholdClusterPRF2;   }
	  Float_t  GetProcent() const                                        { return fProcent;                }
          Short_t  GetDifference() const                                     { return fDifference;             }
          Short_t  GetNumberClusters() const                                 { return fNumberClusters;         }
          Short_t  GetNumberBinCharge() const                                { return fNumberBinCharge;        }
          Short_t  GetNumberBinPRF() const                                   { return fNumberBinPRF;           }
	  Short_t  GetNumberGroupsPRF() const                                { return fNgroupprf;              }
	  Int_t    *GetEntriesLinearFitter() const                           { return fEntriesLinearFitter;    }

 // Debug
          void     SetDebugLevel(Short_t level)                              { fDebugLevel = level;           }

  // Vector method
AliTRDCalibraVector *GetCalibraVector() const                                { return fCalibraVector;          }   
  
 protected:

  // Geometry
  AliTRDgeometry  *fGeo;                    //! The TRD geometry

  // Choice to fill or not the 2D
          Bool_t   fMITracking;             // Chose to fill the 2D histos or vectors during the offline MI tracking
          Bool_t   fMcmTracking;            // Chose to fill the 2D histos or vectors during the tracking with tracklets
          Bool_t   fMcmCorrectAngle;        // Apply correction due to the mcmtrackletangle in the z direction (only) assuming  from vertex
          Bool_t   fCH2dOn;                 // Chose to fill the 2D histos or vectors for the relative gain calibration 
          Bool_t   fPH2dOn;                 // Chose to fill the 2D histos or vectors for the drift velocity and T0
          Bool_t   fPRF2dOn;                // Chose to fill the 2D histos or vectors for the pad response function calibration
          Bool_t   fHisto2d;                // Chose to fill the 2D histos
          Bool_t   fVector2d;               // Chose to fill vectors
	  Bool_t   fLinearFitterOn;         // Method with linear fit for drift velocity
	  Bool_t   fLinearFitterDebugOn;    // Method with linear fit for drift velocity

  // How to fill the 2D
          Float_t  fRelativeScale;          // Scale of the deposited charge
          Float_t  fThresholdClusterPRF2;   // Threshold on cluster pad signals
  // Calibration mode
	  AliTRDCalibraMode *fCalibraMode;  // Calibration mode

  //For debugging
	  TTreeSRedirector          *fDebugStreamer;                 //!Debug streamer
          Short_t     fDebugLevel;                                   // Flag for debugging
  //
  // Internal variables
  //

  // Fill the 2D histos in the offline tracking
          Bool_t   fDetectorAliTRDtrack;    // Change of track
	  Int_t    fDetectorPreviousTrack;  // Change of detector
	  Int_t    fMCMPrevious;            // Change of MCM
	  Int_t    fROBPrevious;            // Change of ROB
	  Short_t  fNumberClusters;         // Minimum number of clusters in the tracklets
          Float_t  fProcent;                // Limit to take the info of the most important calibration group if the track goes through 2 groups (CH)
          Short_t  fDifference;             // Limit to take the info of the most important calibration group if the track goes through 2 groups (CH)
          Int_t    fNumberTrack;            // How many tracks could be used (Debug for the moment)
          Int_t    fNumberUsedCh[2];        // How many tracks have been really used for the gain (0, strict; 1 with fProcent)
          Int_t    fNumberUsedPh[2];        // How many tracks have been really used for the drift velocity (0, strict; 1 with fDifference)
	  Int_t    fTimeMax;                // Number of time bins
          Float_t  fSf;                     // Sampling frequence
	  Short_t  fNumberBinCharge;        // Number of bins for the gain factor
	  Short_t  fNumberBinPRF;           // Number of bin for the PRF
	  Short_t  fNgroupprf;              // Number of groups in tnp bins for PRF /2.0

  // Variables per tracklet
	  TObjArray     *fListClusters;              // List of clusters
	  Double_t      *fPar0;                      // List of track parameter fP[0]
	  Double_t      *fPar1;                      // List of track parameter fP[1]
	  Double_t      *fPar2;                      // List of track parameter fP[2]
	  Double_t      *fPar3;                      // List of track parameter fP[3]
	  Double_t      *fPar4;                      // List of track paarmeter fP[4]
	  Float_t       *fAmpTotal;                  // Energy deposited in the calibration group by the track
          Short_t       *fPHPlace;                   // Calibration group of PH
          Float_t       *fPHValue;                   // PH
	  Bool_t         fGoodTracklet;              // Good tracklet
 // Variables per track
	  Bool_t         fGoodTrack;                 // no return

 //Statistics
	  Int_t         *fEntriesCH;                 // Number of entries CH
	  Int_t         *fEntriesLinearFitter;       // Number of entries LinearFitter


  //
  // Vector method
  //
  	  
	  AliTRDCalibraVector *fCalibraVector; // The vector object
 
 
  // Histograms to store the info from the digits, from the tracklets or from the tracks
	  TProfile2D      *fPH2d;                         // 2D average pulse height
	  TProfile2D      *fPRF2d;                        // 2D PRF
	  TH2I            *fCH2d;                         // 2D deposited charge
	  TObjArray       fLinearFitterArray;             // TObjArray of Linear Fitters for the detectors 
	  AliTRDCalibraVdriftLinearFit *fLinearVdriftFit; // Info Linear Fit
	  
 // Current calib object: to correct for the database used
	  AliTRDCalDet *fCalDetGain;                      // Current calib object gain
	  AliTRDCalROC *fCalROCGain;                      // Current calib object gain
	  AliTRDCalDet *fCalDetT0;                        // Current calib object T0
	  AliTRDCalROC *fCalROCT0;                        // Current calib object T0

 
  //
  // A lot of internal functions......
  //
  // Create the 2D histo to be filled Online
          void     CreateCH2d(Int_t nn);
          void     CreatePH2d(Int_t nn);
          void     CreatePRF2d(Int_t nn);
  
  // Fill the 2D
          void     FillTheInfoOfTheTrackPH();
          void     FillTheInfoOfTheTrackCH();
	  void     FillCH2d(Int_t x, Float_t y);
	  Bool_t   FindP1TrackPH();
	  Bool_t   FindP1TrackPHtrack(AliTRDtrack *t, Int_t index0, Int_t index1);
          void     ResetfVariables();
	  void     ResetfVariablestrack();
          Bool_t   LocalisationDetectorXbins(Int_t detector);
	  Int_t   *CalculateRowCol(AliTRDcluster *cl) const;
	  void     CheckGoodTracklet(Int_t detector, Int_t *rowcol);
	  Int_t    CalculateCalibrationGroup(Int_t i, Int_t *rowcol) const;
	  Int_t    CalculateTotalNumberOfBins(Int_t i);
	  void     StoreInfoCHPH(AliTRDcluster *cl, AliTRDtrack *t, Int_t *group, Int_t *rowcol);
	  void     StoreInfoCHPHtrack(AliTRDcluster *cl, AliTRDtrack *t, Int_t index, Int_t *group, Int_t *rowcol);
	  Bool_t   HandlePRF();
	  Bool_t   HandlePRFtrack(AliTRDtrack *t, Int_t index0, Int_t index1);
 
  // LinearFitter
	  void     AnalyseLinearFitter();
	  
  // Clear
          void     ClearHistos();
      
  // Some basic geometry function
  virtual Int_t    GetPlane(Int_t d) const;
  virtual Int_t    GetChamber(Int_t d) const;
  virtual Int_t    GetSector(Int_t d) const;
	  
          Int_t    Arrondi(Double_t x)const;
 

  // Instance of this class and so on
  static  AliTRDCalibraFillHisto *fgInstance;                // Instance
  static  Bool_t   fgTerminated;                             // If terminated

 private:
  
  // This is a singleton, contructor is private!
  AliTRDCalibraFillHisto();
  virtual ~AliTRDCalibraFillHisto(); 
    
  ClassDef(AliTRDCalibraFillHisto,4)                         // TRD Calibration class

};
  
#endif


