#ifndef ALITRDCALIBRAMODE_H
#define ALITRDCALIBRAMODE_H
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

class AliTRDgeometry;

class AliTRDCalibraMode : public TObject {

 public: 

  AliTRDCalibraMode();
  AliTRDCalibraMode(const AliTRDCalibraMode &c);
  virtual ~AliTRDCalibraMode();
  AliTRDCalibraMode &operator=(const AliTRDCalibraMode &) { return *this; }

   
  Bool_t   ModePadFragmentation(Int_t iLayer,Int_t iStack, Int_t iSector, Int_t i);
  void     ModePadCalibration(Int_t iStack, Int_t i);
  void     ReconstructionRowPadGroup(Int_t idect, Int_t i);
  void     CalculXBins(Int_t idect, Int_t i);
  void     ResetMinMax(Int_t i);
  
          
  //
  // Set and Get methods
  //
  
  //Set
          void     SetNz(Int_t i, Short_t nz);
          void     SetNrphi(Int_t i, Short_t nrphi);
	  void     SetRowMin(Int_t i, Short_t rowmin)                        { fRowMin[i] = rowmin;            }
	  void     SetRowMax(Int_t i, Short_t rowmax)                        { fRowMax[i] = rowmax;            }
	  void     SetColMin(Int_t i, Short_t colmin)                        { fColMin[i] = colmin;            }
	  void     SetColMax(Int_t i, Short_t colmax)                        { fColMax[i] = colmax;            } 
	  void     SetDetChamb0(Int_t i);
	  void     SetDetChamb2(Int_t i);
	  
	  void     SetPerSuperModule(Int_t i);
	  void     SetAllTogether(Int_t i);


  //Get
	  Short_t  GetNz(Int_t i) const                                      { return fNz[i];                  }
          Short_t  GetNrphi(Int_t i) const                                   { return fNrphi[i];               }
          Short_t  GetNnZ(Int_t i) const                                     { return fNnZ[i];                 }
          Short_t  GetNnRphi(Int_t i) const                                  { return fNnRphi[i];              }
          Short_t  GetNfragZ(Int_t i) const                                  { return fNfragZ[i];              }
          Short_t  GetNfragRphi(Int_t i) const                               { return fNfragRphi[i];           }
	  Short_t  GetRowMin(Int_t i) const                                  { return fRowMin[i];              }
	  Short_t  GetRowMax(Int_t i) const                                  { return fRowMax[i];              }
	  Short_t  GetColMin(Int_t i) const                                  { return fColMin[i];              }
	  Short_t  GetColMax(Int_t i) const                                  { return fColMax[i];              }
	  Short_t  GetXbins(Int_t i) const                                   { return fXbins[i];               }
          Short_t  GetDetChamb0(Int_t i) const                               { return fDetChamb0[i];           }
          Short_t  GetDetChamb2(Int_t i) const                               { return fDetChamb2[i];           }
    
     
 protected:

  // Geometry
  AliTRDgeometry  *fGeo;                    //! The TRD geometry

          Short_t  fNz[3];                  // Mode of calibration 
          Short_t  fNrphi[3];               // Mode of calibration 
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

  // Some basic geometry function
  virtual Int_t    GetLayer(Int_t d) const;
  virtual Int_t    GetStack(Int_t d) const;
  virtual Int_t    GetSector(Int_t d) const;
 
  ClassDef(AliTRDCalibraMode,2)             // TRD Calibration class

};
#endif


