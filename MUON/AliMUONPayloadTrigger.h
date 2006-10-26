#ifndef ALIMUONPAYLOADTRIGGER_H
#define ALIMUONPAYLOADTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONPayloadTrigger
/// \brief Class for decoding trigger rawdata 
///
/// \author Christian Finck
///
///////////////////////////////////////////////////////////////////////////////
///
/// This class decode the payload for trigger raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliMUONDDLTrigger;
class AliMUONLocalStruct;
class AliMUONRegHeader;

class AliMUONPayloadTrigger: public TObject {
  public :
    AliMUONPayloadTrigger();
    virtual ~AliMUONPayloadTrigger();

    Bool_t Decode(UInt_t *buffer);
    void   ResetDDL();

    Int_t GetMaxReg() const {return fMaxReg;}
    Int_t GetMaxLoc() const {return fMaxLoc;}


    void SetMaxReg(Int_t reg);
    void SetMaxLoc(Int_t loc);

    AliMUONLocalStruct*     GetLocalInfo()  const {return fLocalStruct;}
    AliMUONRegHeader*       GetRegHeader()  const {return fRegHeader;}
    AliMUONDDLTrigger*      GetDDLTrigger() const {return fDDLTrigger;}

  private :

    Int_t fMaxReg;        ///< maximum number of regional cards in DATE file
    Int_t fMaxLoc;        ///< maximum number of local cards in DATE file

    AliMUONDDLTrigger*       fDDLTrigger;   //!< pointer for DDL structure
    AliMUONRegHeader*        fRegHeader;    //!< pointer for regional structure
    AliMUONLocalStruct*      fLocalStruct;  //!< pointer to local structure


    AliMUONPayloadTrigger(const AliMUONPayloadTrigger& stream);
    AliMUONPayloadTrigger& operator = (const AliMUONPayloadTrigger& stream);

    ClassDef(AliMUONPayloadTrigger, 1)    // base class for reading MUON trigger rawdata
};

#endif
