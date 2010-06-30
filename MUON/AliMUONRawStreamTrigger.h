#ifndef ALIMUONRAWSTREAMTRIGGER_H
#define ALIMUONRAWSTREAMTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONRawStreamTrigger
/// \brief Class for reading MUON raw digits
///
//  Author: Christian Finck

#include <TObject.h>
#include "AliMUONPayloadTrigger.h"
#include "AliMUONVRawStreamTrigger.h"
class TArrayS;

class AliRawReader;
class AliMUONDDLTrigger;
class AliMUONDarcHeader;
class AliMUONRegkHeader;
class AliMUONLocalStruct;

class AliMUONRawStreamTrigger: public AliMUONVRawStreamTrigger {
  public :
    AliMUONRawStreamTrigger();
    AliMUONRawStreamTrigger(AliRawReader* rawReader);
    virtual ~AliMUONRawStreamTrigger();

    /// Initialize iterator
    void First();

   /// Returns current DDL object during iteration
    AliMUONDDLTrigger* CurrentDDL() const { return fCurrentDDL; }
    
    /// Returns current DarcHeader object during iteration
    AliMUONDarcHeader* CurrentDarcHeader() const { return fCurrentDarcHeader; }
    
    /// Returns current RegHeader object during iteration
    AliMUONRegHeader* CurrentRegHeader() const { return fCurrentRegHeader; }
    
    /// Returns current LocalStruct object during iteration
    AliMUONLocalStruct* CurrentLocalStruct() const { return fCurrentLocalStruct; }
    
    /// Advance one step in the iteration. Returns false if finished.
    virtual Bool_t Next(UChar_t& id,   UChar_t& dec,     Bool_t& trigY, 
			UChar_t& yPos, UChar_t& sXDev,   UChar_t& xDev,
			UChar_t& xPos, Bool_t& triggerY, Bool_t& triggerX,
			TArrayS& xPattern, TArrayS& yPattern);


    virtual Bool_t   NextDDL();

    /// Return maximum number of DDLs
    Int_t GetMaxDDL() const {return fgkMaxDDL;}
    /// Return maximum number of regional cards in DATE file
    Int_t GetMaxReg() const {return fPayload->GetMaxReg();}
    /// Return maximum number of local cards in DATE file
    Int_t GetMaxLoc() const {return fPayload->GetMaxLoc();}

    //void SetMaxReg(Int_t reg);
    void SetMaxLoc(Int_t loc);

    /// Return pointer for DDL structure
    AliMUONDDLTrigger* GetDDLTrigger() const {return fPayload->GetDDLTrigger();}

    /// Return number of DDL
    Int_t GetDDL() const {return fDDL - 1;}

    /// Return pointer for payload
    AliMUONPayloadTrigger*  GetPayLoad()    const {return fPayload;}

    /// Whether the iteration is finished or not
    Bool_t IsDone() const;

    /// add error message into error logger
    void AddErrorMessage();

    /// Disable Warnings
    void DisableWarnings() {fPayload->DisableWarnings();}
    
    /// error numbers
    enum rawStreamTriggerError {
      kDarcEoWErr   = 6, ///< end of Darc word error 
      kGlobalEoWErr = 7, ///< end of Global word error
      kRegEoWErr    = 8, ///< end of Regional word error 
      kLocalEoWErr  = 9  ///< end of local word error

    };

  private :
    /// Not implemented
    AliMUONRawStreamTrigger(const AliMUONRawStreamTrigger& stream);
    /// Not implemented
    AliMUONRawStreamTrigger& operator = (const AliMUONRawStreamTrigger& stream);

    Bool_t GetNextDDL();
    Bool_t GetNextRegHeader();
    Bool_t GetNextLocalStruct();

 private:

    AliMUONPayloadTrigger* fPayload; ///< pointer to payload decoder
    AliMUONDDLTrigger* fCurrentDDL;          //!< for iterator: current ddl ptr
    Int_t fCurrentDDLIndex;                  //!< for iterator: current ddl index
    AliMUONDarcHeader* fCurrentDarcHeader;   //!< for iterator: current darc ptr
    AliMUONRegHeader* fCurrentRegHeader;     //!< for iterator: current reg ptr
    Int_t fCurrentRegHeaderIndex;            //!< for iterator: current reg index    
    AliMUONLocalStruct* fCurrentLocalStruct; //!< for iterator: current local ptr
    Int_t fCurrentLocalStructIndex;          //!< for iterator: current local index    
    Bool_t fLocalStructRead;                 //!< flag for read out local structure
    Int_t fDDL;                              //!< number of DDL    


    Bool_t fNextDDL;      ///< flag for next DDL to be read

    static const Int_t  fgkMaxDDL;       ///< maximum number of DDLs

    ClassDef(AliMUONRawStreamTrigger, 4)    // base class for reading MUON trigger rawdata
};

#endif
