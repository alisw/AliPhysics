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

class AliTRDarrayF;
class AliTRDarrayI;

class AliTRDCalibraVector : public TObject {

 public: 

  AliTRDCalibraVector();
  AliTRDCalibraVector(const AliTRDCalibraVector &c);
  virtual ~AliTRDCalibraVector();

  AliTRDCalibraVector& operator = (const  AliTRDCalibraVector &source);

  // Fill
  Int_t          SearchBin(Float_t value, Int_t i) const;  
  Bool_t         UpdateVectorCH(Int_t det, Int_t group, Float_t value);
  Bool_t         UpdateVectorPRF(Int_t det, Int_t group, Float_t x, Float_t y);
  Bool_t         UpdateVectorPH(Int_t det, Int_t group, Int_t time, Float_t value);

  // Add
  Bool_t         Add(AliTRDCalibraVector *calvector);
  
  // Fit
  TGraphErrors  *ConvertVectorPHTGraphErrors(Int_t det, Int_t group, const Char_t *name);
  TGraphErrors  *ConvertVectorPRFTGraphErrors(Int_t det, Int_t group, const Char_t *name);
  TH1F          *ConvertVectorCHHisto(Int_t det, Int_t group, const Char_t *name);

  //
  // Set and Get methods
  //

  void           SetNumberBinCharge(Short_t numberbincharge)   { fNumberBinCharge = numberbincharge;            } 
  void           SetNumberBinPRF(Short_t numberbinprf)         { fNumberBinPRF    = numberbinprf;               } 
  void           SetTimeMax(Int_t timemax)                     { fTimeMax         = timemax;                    } 
  void           SetPRFRange(Float_t prfrange)                 { fPRFRange        = prfrange;                   }  
  void           SetDetCha0(Int_t i, Short_t total)            { fDetCha0[i]      = total;                      } 
  void           SetDetCha2(Int_t i, Short_t total)            { fDetCha2[i]      = total;                      } 
  void           SetNamePH(const char* name)                   { fVectorPHEntries.SetName(name);                } 
  void           SetNamePRF(const char* name)                  { fVectorPRFEntries.SetName(name);               } 
  void           SetNameCH(const char* name)                   { fVectorCHEntries.SetName(name);                } 

  Short_t        GetNumberBinCharge()const                     { return fNumberBinCharge;                       }
  Short_t        GetNumberBinPRF()const                        { return fNumberBinPRF;                          }
  Int_t          GetTimeMax()const                             { return fTimeMax;                               } 
  Float_t        GetPRFRange()const                            { return fPRFRange;                              } 
  Short_t        GetDetCha0(Int_t i) const                     { return fDetCha0[i];                            }
  Short_t        GetDetCha2(Int_t i) const                     { return fDetCha2[i];                            }
  const char*    GetNamePH()                                   { return fVectorPHEntries.GetName();             }
  const char*    GetNamePRF()                                  { return fVectorPRFEntries.GetName();            }
  const char*    GetNameCH()                                   { return fVectorCHEntries.GetName();             }

  AliTRDarrayI  *GetPHEntries(Int_t det,Bool_t force = kFALSE);
  AliTRDarrayF  *GetPHMean(Int_t det,Bool_t force = kFALSE);
  AliTRDarrayF  *GetPHSquares(Int_t det,Bool_t force = kFALSE);

  AliTRDarrayI  *GetPRFEntries(Int_t det,Bool_t force = kFALSE);
  AliTRDarrayF  *GetPRFMean(Int_t det,Bool_t force = kFALSE);
  AliTRDarrayF  *GetPRFSquares(Int_t det,Bool_t force = kFALSE);

  AliTRDarrayI  *GetCHEntries(Int_t det,Bool_t force = kFALSE);

 protected:

  // Current data

         AliTRDarrayI     *fPHEntries;                //  Current AliTRDArrayI PH entries
         AliTRDarrayF     *fPHMean;                   //  Current AliTRDArrayF PH Mean
         AliTRDarrayF     *fPHSquares;                //  Current AliTRDArrayF PH Squares

	 AliTRDarrayI     *fPRFEntries;               //  Current AliTRDArrayI PH entries
         AliTRDarrayF     *fPRFMean;                  //  Current AliTRDArrayF PH Mean
         AliTRDarrayF     *fPRFSquares;               //  Current AliTRDArrayF PH Squares

	 AliTRDarrayI     *fCHEntries;                //  Current AliTRDArrayI PH entries

	 Int_t            fDetectorPH;                //  Current detector
	 Int_t            fDetectorCH;                //  Current detector
	 Int_t            fDetectorPRF;               //  Current detector
       	 
  // Arrays of these objects 
	  
	  TObjArray       fVectorPHMean;              // array of AliTRDarrayF of mean
	  TObjArray       fVectorPHSquares;           // array of AliTRDarrayF of squares
	  TObjArray       fVectorPHEntries;           // array of AliTRDarrayI of entries
	  TObjArray       fVectorCHEntries;           // array of AliTRDarrayI of entries
	  TObjArray       fVectorPRFMean;             // array of AliTRDarrayF of mean
	  TObjArray       fVectorPRFSquares;          // array of AliTRDarrayF of squares
	  TObjArray       fVectorPRFEntries;          // array of AliTRDarrayI of entries
	 
  // Size of the infos

	  Short_t          fNumberBinCharge;          // Number of bins for the gain factor
	  Short_t          fNumberBinPRF;             // Number of bin for the PRF
	  Int_t            fTimeMax;                  // Number of time bins
	  Float_t          fPRFRange;                 // Range PRF
	  Short_t          fDetCha0[3];               // Number of XBins for chamber != 2
          Short_t          fDetCha2[3];               // Number of Xbins for chamber 2

  // Some functions


	  AliTRDarrayI  *GetEntriesPH(Int_t det, TObjArray *array, Bool_t force);
	  AliTRDarrayF  *GetMeanSquaresPH(Int_t det, TObjArray *array, Bool_t force);

	  AliTRDarrayI  *GetEntriesPRF(Int_t det, TObjArray *array, Bool_t force);
	  AliTRDarrayF  *GetMeanSquaresPRF(Int_t det, TObjArray *array, Bool_t force);
	  
	  AliTRDarrayI  *GetEntriesCH(Int_t det, TObjArray *array, Bool_t force);

  // Some basic geometry function
          virtual Int_t    GetChamber(Int_t d) const;
  
  ClassDef(AliTRDCalibraVector,1)                   // TRD Calibration class

};
  
#endif


