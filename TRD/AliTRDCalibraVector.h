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

class TArrayF;
class TArrayI;

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
  void           SetNamePH(const char* name)                   { fNamePH = name;  } 
  void           SetNamePRF(const char* name)                  { fNamePRF = name; } 
  void           SetNameCH(const char* name)                   { fNameCH = name;  } 

  Short_t        GetNumberBinCharge()const                     { return fNumberBinCharge;                       }
  Short_t        GetNumberBinPRF()const                        { return fNumberBinPRF;                          }
  Int_t          GetTimeMax()const                             { return fTimeMax;                               } 
  Float_t        GetPRFRange()const                            { return fPRFRange;                              } 
  Short_t        GetDetCha0(Int_t i) const                     { return fDetCha0[i];                            }
  Short_t        GetDetCha2(Int_t i) const                     { return fDetCha2[i];                            }
  const char*    GetNamePH()                                   { return fNamePH;             }
  const char*    GetNamePRF()                                  { return fNamePRF;            }
  const char*    GetNameCH()                                   { return fNameCH;             }


 
  TArrayI  *GetPHEntries(Int_t det,Bool_t force = kFALSE);
  TArrayF  *GetPHMean(Int_t det,Bool_t force = kFALSE);
  TArrayF  *GetPHSquares(Int_t det,Bool_t force = kFALSE);

  TArrayI  *GetPRFEntries(Int_t det,Bool_t force = kFALSE);
  TArrayF  *GetPRFMean(Int_t det,Bool_t force = kFALSE);
  TArrayF  *GetPRFSquares(Int_t det,Bool_t force = kFALSE);

  TArrayI  *GetCHEntries(Int_t det,Bool_t force = kFALSE);

 protected:

  // Current data

         TArrayI     *fPHEntries[540];                //  Current AliTRDArrayI PH entries
         TArrayF     *fPHMean[540];                   //  Current AliTRDArrayF PH Mean
         TArrayF     *fPHSquares[540];                //  Current AliTRDArrayF PH Squares

	 TArrayI     *fPRFEntries[540];               //  Current AliTRDArrayI PH entries
         TArrayF     *fPRFMean[540];                  //  Current AliTRDArrayF PH Mean
         TArrayF     *fPRFSquares[540];               //  Current AliTRDArrayF PH Squares

	 TArrayI     *fCHEntries[540];                //  Current AliTRDArrayI PH entries

	 const char *fNameCH;                         //  Name for calibration mode
	 const char *fNamePH;                         //  Name for calibration mode
	 const char *fNamePRF;                        //  Name for calibration mode

	 Int_t            fDetectorPH;                //  Current detector
	 Int_t            fDetectorCH;                //  Current detector
	 Int_t            fDetectorPRF;               //  Current detector
       	 
  // Size of the infos

	  Short_t          fNumberBinCharge;          // Number of bins for the gain factor
	  Short_t          fNumberBinPRF;             // Number of bin for the PRF
	  Int_t            fTimeMax;                  // Number of time bins
	  Float_t          fPRFRange;                 // Range PRF
	  Short_t          fDetCha0[3];               // Number of XBins for chamber != 2
          Short_t          fDetCha2[3];               // Number of Xbins for chamber 2

  // Some functions


	  TArrayI  *GetEntriesPH(Int_t det, TArrayI** array, Bool_t force);
	  TArrayF  *GetMeanSquaresPH(Int_t det, TArrayF** array, Bool_t force);

	  TArrayI  *GetEntriesPRF(Int_t det, TArrayI** array, Bool_t force);
	  TArrayF  *GetMeanSquaresPRF(Int_t det, TArrayF** array, Bool_t force);
	  
	  TArrayI  *GetEntriesCH(Int_t det, TArrayI** array, Bool_t force);

  // Some basic geometry function
          virtual Int_t    GetStack(Int_t d) const;
  
  ClassDef(AliTRDCalibraVector,1)                   // TRD Calibration class

};
  
#endif


