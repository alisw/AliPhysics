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

class TTree;
class TProfile2D;
class TGraphErrors;
class TGraph;
class TObjArray;
class TH1F;
class TH2I;
class TH2;

class AliLog;
class AliTRDCalibraMode;
class AliTRDCalibraVector;

class AliTRDcluster;
class AliTRDtrack;
class AliTRDmcmTracklet;

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

  //For the statistics
	  Double_t *StatH(TH2 *ch, Int_t i);
     	 
  //
  // Set of Get the variables
  //

  // Choice to fill or not the 2D
          void     SetMITracking(Bool_t mitracking = kTRUE)                  { fMITracking      = mitracking;  }
          void     SetMcmTracking(Bool_t mcmtracking = kTRUE)                { fMcmTracking     = mcmtracking; }
          void     SetMcmCorrectAngle()                                      { fMcmCorrectAngle = kTRUE;       }
          void     SetPH2dOn()                                               { fPH2dOn          = kTRUE;       }
          void     SetCH2dOn()                                               { fCH2dOn          = kTRUE;       }
          void     SetPRF2dOn()                                              { fPRF2dOn         = kTRUE;       }
          void     SetHisto2d()                                              { fHisto2d         = kTRUE;       }
          void     SetVector2d()                                             { fVector2d        = kTRUE;       }
  
          Bool_t   GetMITracking() const                                     { return fMITracking;             }
          Bool_t   GetMcmTracking() const                                    { return fMcmTracking;            }
          Bool_t   GetMcmCorrectAngle() const                                { return fMcmCorrectAngle;        }
          Bool_t   GetPH2dOn() const                                         { return fPH2dOn;                 }
          Bool_t   GetCH2dOn() const                                         { return fCH2dOn;                 }
          Bool_t   GetPRF2dOn() const                                        { return fPRF2dOn;                }
          Bool_t   GetHisto2d() const                                        { return fHisto2d;                }
          Bool_t   GetVector2d() const                                       { return fVector2d;               }
  TH2I            *GetCH2d() const                                           { return fCH2d;                   }
  TProfile2D      *GetPH2d() const                                           { return fPH2d;                   }
  TProfile2D      *GetPRF2d() const                                          { return fPRF2d;                  }
  
  // How to fill the 2D
          void     SetRelativeScaleAuto()                                    { fRelativeScaleAuto    = kTRUE;                }
          void     SetRelativeScale(Float_t relativeScale);                      
	  void     SetThresholdClusterPRF1(Float_t thresholdClusterPRF1)     { fThresholdClusterPRF1 = thresholdClusterPRF1; }
          void     SetThresholdClusterPRF2(Float_t thresholdClusterPRF2)     { fThresholdClusterPRF2 = thresholdClusterPRF2; }
          void     SetCenterOfflineCluster()                                 { fCenterOfflineCluster = kTRUE;                }
          void     SetNz(Int_t i, Short_t nz);
          void     SetNrphi(Int_t i, Short_t nrphi);
          void     SetProcent(Float_t procent)                               { fProcent              = procent;              }
          void     SetDifference(Short_t difference)                         { fDifference           = difference;           }
          void     SetNumberClusters(Short_t numberClusters)                 { fNumberClusters       = numberClusters;       }
          void     SetNumberBinCharge(Short_t numberBinCharge)               { fNumberBinCharge      = numberBinCharge;      }
          void     SetNumberBinPRF(Short_t numberBinPRF)                     { fNumberBinPRF         = numberBinPRF;         }
  
          Float_t  GetRelativeScale() const                                  { return fRelativeScale;          }
          Bool_t   GetRelativeScaleAuto() const                              { return fRelativeScaleAuto;      }
	  Float_t  GetThresholdClusterPRF1() const                           { return fThresholdClusterPRF1;   }
          Float_t  GetThresholdClusterPRF2() const                           { return fThresholdClusterPRF2;   }
	  Float_t  GetProcent() const                                        { return fProcent;                }
          Short_t  GetDifference() const                                     { return fDifference;             }
          Short_t  GetNumberClusters() const                                 { return fNumberClusters;         }
          Short_t  GetNumberBinCharge() const                                { return fNumberBinCharge;        }
          Short_t  GetNumberBinPRF() const                                   { return fNumberBinPRF;           }
  
  // Write
	  void     SetWrite(Int_t i)                                         { fWrite[i]      = kTRUE;         }
          void     SetWriteName(TString writeName)                           { fWriteName     = writeName;     }
  
	  Bool_t   GetWrite(Int_t i) const                                   { return fWrite[i];               }
          TString  GetWriteName() const                                      { return fWriteName;              }

 //  Calibration mode
AliTRDCalibraMode  *GetCalibraMode() const                                   { return fCalibraMode;            }

// Vector method
AliTRDCalibraVector *GetCalibraVector() const                                { return fCalibraVector;          }   
  
 private:
   
  // This is a singleton, contructor is private!
  AliTRDCalibraFillHisto();
  virtual ~AliTRDCalibraFillHisto();

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
	  Float_t  fThresholdClusterPRF1;   // Threshold on cluster pad signals for PRF peripherique
          Float_t  fThresholdClusterPRF2;   // Threshold on cluster pad signals for PRF peripherique
          Bool_t   fCenterOfflineCluster;   // Choose to use the offline determination of the center of the cluster
	
  // Write
	  Bool_t   fWrite[3];               // Do you want to write the 2D histo or vectors converted in a tree
          TString  fWriteName;              // Where the 2D or trees are written

	  // Calibration mode
	  AliTRDCalibraMode *fCalibraMode;  // Calibration mode

  //
  // Internal variables
  //

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
	  Int_t    fTimeMax;                // Number of time bins
          Float_t  fSf;                     // Sampling frequence
	  Short_t  fNumberBinCharge;        // Number of bins for the gain factor
	  Short_t  fNumberBinPRF;           // Number of bin for the PRF

  //
  // Vector method
  //
  
	  
	  AliTRDCalibraVector *fCalibraVector; // The vector object
 
 
  // Histograms to store the info from the digits, from the tracklets or from the tracks
  TProfile2D      *fPH2d;                   // 2D average pulse height
  TProfile2D      *fPRF2d;                  // 2D PRF
  TH2I            *fCH2d;                   // 2D deposited charge 
          
  //
  // A lot of internal functions......
  //

  // Init AliTRDCalibraFillHisto
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
 
  // Clear
          void     ClearHistos();
      
  // Some basic geometry function
  virtual Int_t    GetPlane(Int_t d) const;
  virtual Int_t    GetChamber(Int_t d) const;
  virtual Int_t    GetSector(Int_t d) const;
 

  // Instance of this class and so on
  static  AliTRDCalibraFillHisto *fgInstance;        // Instance
  static  Bool_t   fgTerminated;                     // If terminated
    
  ClassDef(AliTRDCalibraFillHisto,1)                 // TRD Calibration class

};
  
#endif


