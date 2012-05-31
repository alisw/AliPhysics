#ifndef ALITRDCALIBRAVECTOR_H
#define ALITRDCALIBRAVECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#  include <TObject.h>
#endif

class TGraphErrors;
class TH1F;
class TObjArray;

class AliTRDPhInfo;
class AliTRDEntriesInfo;
class AliTRDPrfInfo;

class AliTRDCalibraVector : public TObject {

 public: 

  AliTRDCalibraVector();
  AliTRDCalibraVector(const AliTRDCalibraVector &c);
  virtual ~AliTRDCalibraVector();
  virtual Long64_t Merge(const TCollection* list);

  AliTRDCalibraVector& operator = (const  AliTRDCalibraVector &source);

  // Init
  void           TestInit(Int_t i, Int_t detmax);

  // Fill
  Bool_t         UpdateVectorCH(Int_t det, Int_t group, Float_t value);
  Bool_t         UpdateVectorPRF(Int_t det, Int_t group, Float_t x, Float_t y);
  Bool_t         UpdateVectorPH(Int_t det, Int_t group, Int_t time, Float_t value);


  Bool_t         FillVectorCH(Int_t det, Int_t group, Int_t bin, Int_t entries);
  Bool_t         FillVectorPRF(Int_t det, Int_t group, Int_t bin, Int_t entries, Float_t mean, Float_t square);
  Bool_t         FillVectorPH(Int_t det, Int_t group, Int_t bin, Int_t entries, Float_t mean, Float_t square);


  // Add
  Bool_t         Add(AliTRDCalibraVector *calvector);

  AliTRDCalibraVector *AddStatsPerDetectorCH();
  AliTRDCalibraVector *AddStatsPerDetectorPH();
  AliTRDCalibraVector *AddStatsPerDetectorPRF();

  
  // Fit
  TGraphErrors  *ConvertVectorPHTGraphErrors(Int_t det, Int_t group, const Char_t *name);
  TGraphErrors  *ConvertVectorPRFTGraphErrors(Int_t det, Int_t group, const Char_t *name);
  TH1F          *CorrectTheError(const TGraphErrors *hist, Int_t &nbEntries);
  TH1F          *ConvertVectorCHHisto(Int_t det, Int_t group, const Char_t *name);

  // Find
  Int_t          SearchBin(Float_t value, Int_t i) const;  
  Bool_t         FindTheMaxEntries(Int_t i, Int_t &detectormax, Int_t &groupmax);

  //
  // Set and Get methods
  //

  void           SetNumberBinCharge(Short_t numberbincharge)   { fNumberBinCharge = numberbincharge;            } 
  void           SetNumberBinPRF(Short_t numberbinprf)         { fNumberBinPRF    = numberbinprf;               } 
  void           SetTimeMax(Int_t timemax)                     { fTimeMax         = timemax;                    } 
  void           SetPRFRange(Float_t prfrange)                 { fPRFRange        = prfrange;                   }  
  void           SetDetCha0(Int_t i, Short_t total)            { fDetCha0[i]      = total;                      } 
  void           SetDetCha2(Int_t i, Short_t total)            { fDetCha2[i]      = total;                      } 
  void           SetNzNrphi(Int_t i, Int_t nz, Int_t nrphi);    
  void           SetNbGroupPRF(Int_t nbGroup)                  { fNbGroupPRF = (UChar_t) nbGroup;               } 

  Short_t        GetNumberBinCharge()const                     { return fNumberBinCharge;                       }
  Short_t        GetNumberBinPRF()const                        { return fNumberBinPRF;                          }
  Int_t          GetTimeMax()const                             { return fTimeMax;                               } 
  Float_t        GetPRFRange()const                            { return fPRFRange;                              } 
  Short_t        GetDetCha0(Int_t i) const                     { return fDetCha0[i];                            }
  Short_t        GetDetCha2(Int_t i) const                     { return fDetCha2[i];                            }
  TString        GetNamePH() const;    
  TString        GetNamePRF() const;   
  TString        GetNameCH() const;  
  Int_t          GetNz(Int_t i) const;
  Int_t          GetNrphi(Int_t i) const;
  Int_t          GetNbGroupPRF() const                         { return (Int_t)fNbGroupPRF;                     }

  Int_t          GetTotalNumberOfBinsInDetector(Int_t det, Int_t i, Int_t nbBin) const;

  TObject               *GetPHEntries(Int_t det,Bool_t force = kFALSE);
  TObject               *GetPHMean(Int_t det,Bool_t force = kFALSE);
  TObject               *GetPHSquares(Int_t det,Bool_t force = kFALSE);
  
  TObject               *GetPRFEntries(Int_t det,Bool_t force = kFALSE);
  TObject               *GetPRFMean(Int_t det,Bool_t force = kFALSE);
  TObject               *GetPRFSquares(Int_t det,Bool_t force = kFALSE);
  
  TObject               *GetCHEntries(Int_t det,Bool_t force = kFALSE);

   
  
 protected:
  
  // Current data
  
  AliTRDEntriesInfo     *fPHEntries[540];                //  PH entries
  AliTRDPhInfo          *fPHMean[540];                   //  PH Mean
  AliTRDPhInfo          *fPHSquares[540];                //  PH Squares
  
  AliTRDEntriesInfo     *fPRFEntries[540];               //  PRF entries
  AliTRDPrfInfo         *fPRFMean[540];                  //  PRF Mean
  AliTRDPrfInfo         *fPRFSquares[540];               //  PRF Squares
  
  AliTRDEntriesInfo     *fCHEntries[540];                //  CH entries
  
  UChar_t      fModeCH;                         //  Calibration mode
  UChar_t      fModePH;                         //  Calibration mode
  UChar_t      fModePRF;                        //  Calibration mode
  UChar_t      fNbGroupPRF;                     //  Nb of group PRD
  
  
  Int_t            fDetectorPH;                //!  Current detector
  Int_t            fDetectorCH;                //!  Current detector
  Int_t            fDetectorPRF;               //!  Current detector
  Bool_t           fStopFillCH;                //!  To know if we stop to fill
  TH1F            *fHisto;                     //!  Histo to be fitted
  TGraphErrors    *fGraph;                     //!  TGraphError
  AliTRDCalibraVector *fCalVector;             //!  AliTRDCalibraVector
  
  // Size of the infos
  
  Short_t          fNumberBinCharge;          // Number of bins for the gain factor
  Short_t          fNumberBinPRF;             // Number of bin for the PRF
  Int_t            fTimeMax;                  // Number of time bins
  Float_t          fPRFRange;                 // Range PRF
  Short_t          fDetCha0[3];               // Number of XBins for chamber != 2
  Short_t          fDetCha2[3];               // Number of Xbins for chamber 2
  
  // Some functions
  
  AliTRDEntriesInfo  *GetEntriesPH(Int_t det, AliTRDEntriesInfo** array, Bool_t force);
  AliTRDPhInfo       *GetMeanSquaresPH(Int_t det, AliTRDPhInfo** array, Bool_t force);
  
  AliTRDEntriesInfo  *GetEntriesPRF(Int_t det, AliTRDEntriesInfo** array, Bool_t force);
  AliTRDPrfInfo      *GetMeanSquaresPRF(Int_t det, AliTRDPrfInfo** array, Bool_t force);
  
  AliTRDEntriesInfo  *GetEntriesCH(Int_t det, AliTRDEntriesInfo** array, Bool_t force);

  ClassDef(AliTRDCalibraVector,1)                   // TRD Calibration class
    
    };
    
#endif

