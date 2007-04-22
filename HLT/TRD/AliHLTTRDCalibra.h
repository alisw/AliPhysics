#ifndef ALIHLTTRDCALIBRA_H
#define ALIHLTTRDCALIBRA_H
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

#include "TProfile2D.h"
#include "TH2I.h"

class TProfile2D;
class TH2I;

class AliLog;

class AliTRDmcmTracklet;

class AliHLTTRDCalibra : public TObject {

 public: 

  // Instance
  static AliHLTTRDCalibra *Instance();
  static void Terminate();
  static void Destroy();

  AliHLTTRDCalibra(const AliHLTTRDCalibra &c);
  AliHLTTRDCalibra &operator=(const AliHLTTRDCalibra &) { return *this; }

  // Functions for initialising the AliHLTTRDCalibra in the code
          Bool_t   Init2Dhistos();

  // Functions for filling the histos in the code
	  Bool_t   UpdateHistogramcm(AliTRDmcmTracklet *trk);
 
  // Is Pad on
          Bool_t   IsPadOn(Int_t detector, Int_t col, Int_t row) const;

  // Functions for writting the 2D
          Bool_t   Write2d();

  // Pad Calibration
          Bool_t   ModePadFragmentation(Int_t iPlane,Int_t iChamb, Int_t iSect, Int_t i);
          void     ModePadCalibration(Int_t iChamb, Int_t i);
          
  //
  // Set of Get the variables
  //

  // Choice to fill or not the 2D
	  void     SetOn()                                                   { fOn          = kTRUE;           }
	  Bool_t   GetOn() const                                             { return fOn;                     }
          
  //Get the 2D histos      
	  TH2I            *GetCH2d() const                                           { return fCH2d;                   }
	  TProfile2D      *GetPH2d() const                                           { return fPH2d;                   }
	  TProfile2D      *GetPRF2d() const                                          { return fPRF2d;                  }
  
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
  
	  Float_t  GetThresholdClusterPRF2() const                           { return fThresholdClusterPRF2;   }
	  Short_t  GetNz(Int_t i) const                                      { return fNz[i];                  }
          Short_t  GetNrphi(Int_t i) const                                   { return fNrphi[i];               }
          Float_t  GetProcent() const                                        { return fProcent;                }
          Short_t  GetDifference() const                                     { return fDifference;             }
          Short_t  GetNumberClusters() const                                 { return fNumberClusters;         }
          Short_t  GetNumberBinCharge() const                                { return fNumberBinCharge;        }
          Short_t  GetNumberBinPRF() const                                   { return fNumberBinPRF;           }
  
  // Write
	  void     SetWriteName(TString writeName)                           { fWriteName     = writeName;     }
	  TString  GetWriteName() const                                      { return fWriteName;              }
          
  //
  // Internal variables to be sure!
  //
  
  // Pad calibration
          Short_t  GetNnz(Int_t i) const                                     { return fNnZ[i];                 }
          Short_t  GetNnrphi(Int_t i) const                                  { return fNnRphi[i];              }
          Short_t  GetNfragz(Int_t i) const                                  { return fNfragZ[i];              }
          Short_t  GetNfragrphi(Int_t i) const                               { return fNfragRphi[i];           }
          Short_t  GetDetChamb0(Int_t i) const                               { return fDetChamb0[i];           }
          Short_t  GetDetChamb2(Int_t i) const                               { return fDetChamb2[i];           }

 private:
   
  // This is a singleton, contructor is private!
	  AliHLTTRDCalibra();
	  virtual ~AliHLTTRDCalibra();

 protected:

  // Choice to fill or not the 2D
          Bool_t   fOn;                     // Chose to fill the 2D histos 
            
  // How to fill the 2D
	  Float_t  fRelativeScale;          // Scale of the deposited charge
          Float_t  fThresholdClusterPRF2;   // Threshold on cluster pad signals for PRF peripherique
	  Short_t  fNz[3];                  // Mode of calibration 
          Short_t  fNrphi[3];               // Mode of calibration 
          Int_t    fNtotal[3];              // Total number of Xbins

  // Write
	  TString  fWriteName;              // Where the 2D histos are written
           
  //
  // Internal variables
  //
 
  // Fill the 2D histos in the offline tracking
	  Bool_t   fGoodTrack;              // If goes through a kaputt pad
          Float_t *fAmpTotal;               //! Energy deposited in the calibration group by the track
          Short_t *fPHPlace;                //! Calibration group of PH
          Float_t *fPHValue;                //! PH
          Short_t  fNumberClusters;         // Minimum number of clusters in the tracklets
          Float_t  fProcent;                // Limit to take the info of the most important calibration group if the track goes through 2 groups (CH)
          Short_t  fDifference;             // Limit to take the info of the most important calibration group if the track goes through 2 groups (CH)
          Int_t    fNumberTrack;            // How many tracks could be used
          Int_t    fNumberUsedCh[2];        // How many tracks have been really used for the gain (0, strict; 1 with fProcent)
          Int_t    fNumberUsedPh[2];        // How many tracks have been really used for the drift velocity (0, strict; 1 with fDifference)

  //
  // For debugging 
  //

 // Variables in the loop for the coef or more general
	  Int_t    fTimeMax;                // Number of time bins
          Float_t  fSf;                     // Sampling frequence
	                     
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
          Short_t  fDetChamb2[3];           // Number of Xbins for chamber 2
  
 
  // CH
          Short_t  fNumberBinCharge;        // Number of bins for the gain factor
 
  // PRF
	  Short_t          fNumberBinPRF;   // Number of bin for the PRF
  
  // Histograms to store the info from the digits, from the tracklets or from the tracks
	  TProfile2D      *fPH2d;           // 2D average pulse height
	  TProfile2D      *fPRF2d;          // 2D PRF
	  TH2I            *fCH2d;          // 2D deposited charge 
    
  //
  // A lot of internal functions......
  //

  // Init AliHLTTRDCalibra
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
 
  // Pad group calibration mode
          void     ReconstructionRowPadGroup(Int_t idect, Int_t i);
          void     CalculXBins(Int_t idect, Int_t i);  

 
  // Clear
          void     ClearHistos();
         
  // Some basic geometry function
  virtual Int_t    GetPlane(Int_t d) const;
  virtual Int_t    GetChamber(Int_t d) const;
  virtual Int_t    GetSector(Int_t d) const;
  
  // Instance of this class and so on
  static  AliHLTTRDCalibra *fgInstance;        // Instance
  static  Bool_t   fgTerminated;            // If terminated
    
  ClassDef(AliHLTTRDCalibra,2)                 // TRD Calibration class

};
  
#endif


