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
class TTree;
class TObjArray;

class AliTRDCalibraVector : public TObject {

 public: 

  AliTRDCalibraVector();
  AliTRDCalibraVector(const AliTRDCalibraVector &c);
  virtual ~AliTRDCalibraVector();
  AliTRDCalibraVector &operator=(const AliTRDCalibraVector &) { return *this; }

  Int_t          SearchBin(Float_t value, Int_t i) const;  
  Int_t          SearchInVector(Int_t group, Int_t i) const;
  Int_t          SearchInTreeVector (TObjArray *vectorplace, Int_t group) const;
  Bool_t         UpdateVectorCH(Int_t group, Float_t value);
  Bool_t         UpdateVectorPRF(Int_t group, Float_t x, Float_t y);
  Bool_t         UpdateVectorPH(Int_t group, Int_t time, Float_t value);
  TH1F          *ConvertVectorCTHisto(Int_t place, const Char_t *name) const;
  TTree         *ConvertVectorCTTreeHisto(TObjArray *vVectorCT, TObjArray *pPlaCT
                                        , const Char_t *name, const Char_t *nametitle) const;
  TGraphErrors  *ConvertVectorPHisto(Int_t place, const Char_t *name) const;
  TTree         *ConvertVectorPTreeHisto(TObjArray *vVectorP, TObjArray *pPlaP
                                       , const Char_t *name, const Char_t *nametitle) const;
  TObjArray     *ConvertTreeVector(TTree *tree) const ;
  Bool_t         MergeVectorCT(TObjArray *vVectorCT2, TObjArray *pPlaCT2);
  Bool_t         MergeVectorP(TObjArray *vVectorP2, TObjArray *pPlaP2, Int_t i);
  //Add two trees
  TTree         *Sum2Trees(const Char_t *filename1, const Char_t *filename2, const Char_t *variablecali);
  //Correct the errors
  TGraphErrors  *AddProfiles(TGraphErrors *hist1, TGraphErrors *hist2) const ;

  //
  // Set and Get methods
  //

  void           SetNumberBinCharge(Short_t numberbincharge)   { fNumberBinCharge = numberbincharge;            } 
  void           SetNumberBinPRF(Short_t numberbinprf)         { fNumberBinPRF    = numberbinprf;               } 
  void           SetTimeMax(Int_t timemax)                     { fTimeMax         = timemax;                    }  
  void           SetVectorPH(TObjArray *vectorPH)              { fVectorPH        = vectorPH;                   }
  void           SetPlaPH(TObjArray *plaPH)                    { fPlaPH           = plaPH;                      }
  void           SetVectorCH(TObjArray *vectorCH)              { fVectorCH        = vectorCH;                   }
  void           SetPlaCH(TObjArray *plaCH)                    { fPlaCH           = plaCH;                      }
  void           SetVectorPRF(TObjArray *vectorPRF)            { fVectorPRF       = vectorPRF;                  }
  void           SetPlaPRF(TObjArray *plaPRF)                  { fPlaPRF          = plaPRF;                     }

  Short_t        GetNumberBinCharge()const                     { return fNumberBinCharge;                       }
  Short_t        GetNumberBinPRF()const                        { return fNumberBinPRF;                          }
  Int_t          GetTimeMax()const                             { return fTimeMax;                               } 
  TObjArray*     GetVectorPH()const                            { return fVectorPH;                              } 
  TObjArray*     GetPlaPH()const                               { return fPlaPH;                                 } 
  TObjArray*     GetVectorCH()const                            { return fVectorCH;                              } 
  TObjArray*     GetPlaCH()const                               { return fPlaCH;                                 } 
  TObjArray*     GetVectorPRF()const                           { return fVectorPRF;                             } 
  TObjArray*     GetPlaPRF()const                              { return fPlaPRF;                                } 
 
 protected:

 // Objects for storing the infos per calibration group
	  
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

  // Arrays of these objects 
	  
	  TObjArray       *fVectorPH;               // Vectors to fill
	  TObjArray       *fPlaPH;                  // Vectors to fill
	  
	  TObjArray       *fVectorCH;               // Vectors to fill
	  TObjArray       *fPlaCH;                  // Vectors to fill
	  
	  TObjArray       *fVectorPRF;              // Vectors to fill
	  TObjArray       *fPlaPRF;                 // Vectors to fill

  // Size of the infos

	  Short_t          fNumberBinCharge;        // Number of bins for the gain factor
	  Short_t          fNumberBinPRF;           // Number of bin for the PRF
	  Int_t            fTimeMax;                // Number of time bins

  // Some functions

	  TGraphErrors  *ConvertVectorPHistoI(AliTRDPInfo *pInfo, const Char_t *name) const;
	  TH1F          *ConvertVectorCTHistoI(AliTRDCTInfo *cTInfo, const Char_t *name) const;

  ClassDef(AliTRDCalibraVector,1)                   // TRD Calibration class

};
  
#endif


