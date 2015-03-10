#ifndef ALIMUONPAYLOADTRIGGER_H
#define ALIMUONPAYLOADTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONPayloadTrigger
/// \brief Class for decoding the payload for trigger raw data 
///
//  Author Christian Finck

#include <TObject.h>

class AliMUONDDLTrigger;
class AliMUONLocalStruct;
class AliMUONRegHeader;
class AliMUONLogger;
class AliMUONLogger;

class AliMUONPayloadTrigger: public TObject {
  public :
    AliMUONPayloadTrigger();
    virtual ~AliMUONPayloadTrigger();

    Bool_t Decode(UInt_t *buffer, Bool_t scalerEvent = kFALSE);
    void   ResetDDL();

    /// Return maximum number of regional cards in DATE file
    Int_t GetMaxReg() const {return fMaxReg;}
    /// Return maximum number of local cards in DATE file
    Int_t GetMaxLoc() const {return fMaxLoc;}


    void SetMaxLoc(Int_t loc);

    /// Return pointer to local structure
    AliMUONLocalStruct*     GetLocalInfo()  const {return fLocalStruct;}
    /// Return pointer for regional structure
    AliMUONRegHeader*       GetRegHeader()  const {return fRegHeader;}
    /// Return pointer for DDL structure
    AliMUONDDLTrigger*      GetDDLTrigger() const {return fDDLTrigger;}

    /// Get number of end of DARC word errors
    Int_t   GetDarcEoWErrors() const {return fDarcEoWErrors;}
    /// Get number of end of Global word errors
    Int_t   GetGlobalEoWErrors() const {return fGlobalEoWErrors;}
    /// Get number of end of regional word errors
    Int_t   GetRegEoWErrors() const {return fRegEoWErrors;}
    /// Get number of end of local word errors
    Int_t   GetLocalEoWErrors() const {return fLocalEoWErrors;}

   /// Get Error logger
    AliMUONLogger* GetErrorLogger() const {return fLog;}
    
    /// set warnings flag
    void DisableWarnings() {fWarnings = kFALSE;} 
    
  private :
    /// Not implemented
    AliMUONPayloadTrigger(const AliMUONPayloadTrigger& stream);
    /// Not implemented
    AliMUONPayloadTrigger& operator = (const AliMUONPayloadTrigger& stream);

    void   AddErrorMessage(const Char_t* msg);
    void   SetMaxReg(Int_t reg);

    Int_t fMaxReg;        ///< maximum number of regional cards in DATE file
    Int_t fMaxLoc;        ///< maximum number of local cards in DATE file

    AliMUONDDLTrigger*       fDDLTrigger;   //!<! pointer for DDL structure
    AliMUONRegHeader*        fRegHeader;    //!<! pointer for regional structure
    AliMUONLocalStruct*      fLocalStruct;  //!<! pointer to local structure

    AliMUONLogger* fLog;                    //!<! Map of errors msg;
    Int_t   fDarcEoWErrors;                 //!<! number of end of DARC word errors;
    Int_t   fGlobalEoWErrors;               //!<! number of end of global word errors;
    Int_t   fRegEoWErrors;                  //!<! number of end of regional word errors;
    Int_t   fLocalEoWErrors;                //!<! number of end of local word errors;
    Bool_t  fWarnings;                      //!<! flag to enable/disable warnings
    Bool_t  fNofRegSet;                     //!<! true if number of regional boards is set from outside

    ClassDef(AliMUONPayloadTrigger, 3)    // base class for reading MUON trigger rawdata
};

#endif
