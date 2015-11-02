#ifndef AliITSMFTCalibrationPix_H
#define AliITSMFTCalibrationPix_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliITSMFTCalibrationPix.h 53595 2011-12-14 16:13:16Z masera $ */
//#include "TRandom.h"
#include "TObject.h"
#include "TString.h"
#include "TArrayS.h"

////////////////////////////////////////////////////
//                                                //
// ITS response class for pixels                  //
////////////////////////////////////////////////////

class AliITSMFTCalibrationPix :  public TObject {
 public:
    AliITSMFTCalibrationPix(); // default constructor
    AliITSMFTCalibrationPix(Short_t nChips,Short_t nColPerChip,Short_t nRow); // default constructor
    AliITSMFTCalibrationPix(const AliITSMFTCalibrationPix &src);
    //    
    virtual ~AliITSMFTCalibrationPix() {;} // destructror

    virtual void   ClearBad();

    virtual Int_t  GetNrBad() const;
    virtual Int_t  GetNrBadInChip(Int_t chip) const;
    virtual Int_t  GetNrBadInColumn(Int_t col) const;

    virtual Int_t  GetBadColAt(Int_t index) const;
    virtual Int_t  GetBadRowAt(Int_t index) const;
    virtual void   GetBadPixel(Int_t index, Int_t &row, Int_t &col) const;
    virtual void   GetBadPixelSingle(Int_t i,UInt_t &row, UInt_t &col) const;

    virtual Int_t  GetNrBadSingle() const {return fNrBadSingle;}
    virtual void   SetNrBadSingle(Int_t nr) {fNrBadSingle=nr;} // used to be called SetNrBad, but misleading
    virtual void   SetBadList(TArrayS badlist) {fBadChannels=badlist;}
    virtual void   SetNrBad(Int_t /*nr*/); // Use SetNrBadSingle!!!

    // Get data type
    virtual const char  *DataType() const {return fDataType.Data();}
    // Type of data - real or simulated
    virtual void    SetDataType(const char *data="simulated") {fDataType=data;}

    virtual Bool_t IsBad() const;
    virtual Bool_t IsChipBad(Int_t chip) const;
    virtual Bool_t IsColumnBad(Int_t col) const;
    virtual Bool_t IsPixelBad(Int_t col, Int_t row) const;

    virtual void   SetChipBad(Int_t chip);
    virtual void   UnSetChipBad(Int_t chip);
    Bool_t         IsChipMarkedBad(Int_t c) const {return (fBadChips&(0x1<<c))!=0;}

    virtual void   AddBad(Int_t col, Int_t row);

    virtual Int_t  GetChipIndexFromCol(Int_t col) const;
    //    virtual Int_t  GetChipFromChipIndex(Int_t index) const;

    virtual void   GiveCompressParam(Int_t *)   const {NotImplemented("GiveCompressParam");}
    virtual void   SetDetParam(Double_t *)            {NotImplemented("SetDetParam");}
    virtual void   GetDetParam(Double_t *)      const {NotImplemented("GetDetParam");}
    virtual void   SetNDetParam(Int_t /*n*/)          {NotImplemented("SetNDetParam");}
    virtual Int_t  NDetParam()                  const {NotImplemented("NDetParam"); return 0;}
    virtual void   SetSigmaSpread(Double_t, Double_t) {NotImplemented("SetSigmaSpread");}
    virtual void   SigmaSpread(Double_t & /* p1 */,Double_t & /* p2 */) const {NotImplemented("SigmaSpread");}
    //
    void    SetColRowData(Short_t nchip, Short_t ncolperchip, Short_t nrow);
    Int_t   GetNCol()         const {return fNCol;}
    Int_t   GetNRow()         const {return fNRow;}
    Int_t   GetNColPerChip()  const {return fNColPerChip;}
    Int_t   GetNChips()       const {return fNChips;}
    //
 protected:
    void NotImplemented(const char *method) const {if(gDebug>0)
         Warning(method,"This method is not implemented for this sub-class");}

    TString  fDataType;   // data type - real or simulated

    Short_t  fNChips;          // n of chips
    Short_t  fNColPerChip;     // n of columns per chip
    Short_t  fNCol;            // number of columns
    Short_t  fNRow;            // number of rows
    Int_t    fNrBadSingle;     // Nr of SINGLE bad pixels
    UInt_t   fBadChips;        // bit pattern of completely dead chips?
    TArrayS  fBadChannels;     // Array with bad channels info (col0,row0,col1...rowN) N = fNrBadSingle
    //
    AliITSMFTCalibrationPix& operator=(const AliITSMFTCalibrationPix& source);
    //
    ClassDef(AliITSMFTCalibrationPix,2) // pixels response
};

#endif
