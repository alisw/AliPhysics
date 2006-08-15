#ifndef AliTRDRECPARAM_H
#define AliTRDRECPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class containing constant reconstruction parameters                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/* $Id$ */

#include "TObject.h"

class AliTRDRecParam : public TObject {

 public:

  static  AliTRDRecParam *Instance();
  static  void Terminate();
  
  enum { kNplan =   6
       , kNcham =   5
       , kNsect =  18
       , kNdet  = 540 };
    
  AliTRDRecParam(const AliTRDRecParam &p);   
  AliTRDRecParam &operator=(const AliTRDRecParam &p);
 
  virtual void     Copy(TObject &p) const;
  
  virtual void     SetLUT(Int_t lutOn = 1)                        { fLUTOn         = lutOn;  };
  virtual void     SetClusMaxThresh(Float_t thresh)               { fClusMaxThresh = thresh; };
  virtual void     SetClusSigThresh(Float_t thresh)               { fClusSigThresh = thresh; };
          void     SetTailCancelation(Int_t tcOn = 1)             { fTCOn          = tcOn;   };
          void     SetNexponential(Int_t nexp)                    { fTCnexp        = nexp;   };
 
          Bool_t   LUTOn() const                                  { return fLUTOn;           };
  virtual Float_t  GetClusMaxThresh() const                       { return fClusMaxThresh;   };
  virtual Float_t  GetClusSigThresh() const                       { return fClusSigThresh;   };
          Bool_t   TCOn() const                                   { return fTCOn;            };
          Int_t    GetTCnexp() const                              { return fTCnexp;          };
    
  virtual Double_t LUTposition(Int_t iplane, Double_t ampL, Double_t ampC, Double_t ampR) const;
  
 protected:

  static  AliTRDRecParam *fgInstance; //  Instance of this class (singleton implementation)
  static  Bool_t   fgTerminated;      //  Defines if this class has already been terminated and
                                      //  therefore does not return instances in GetInstance anymore
    
          void     Init();
  virtual void     FillLUT();
  
          // Clusterization parameter
          Float_t  fClusMaxThresh;    //  Threshold value for cluster maximum
          Float_t  fClusSigThresh;    //  Threshold value for cluster signal
    
          Int_t    fLUTOn;            //  Switch for the lookup table method
          Int_t    fLUTbin;           //  Number of bins of the LUT
          Float_t *fLUT;              //! The lookup table
  
          Int_t    fTCOn;             //  Switch for the tail cancelation
          Int_t    fTCnexp;           //  Number of exponentials, digital filter
  
 private:

  // This is a singleton, constructor is private!  
  AliTRDRecParam();
  ~AliTRDRecParam();
  
  ClassDef(AliTRDRecParam,1)          //  The reconstruction parameter

};

#endif
