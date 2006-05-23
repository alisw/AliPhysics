#ifndef ALIMUONRAWSTREAMTRIGGER_H
#define ALIMUONRAWSTREAMTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONRawStreamTrigger
/// \brief Class for reading MUON raw digits
///
/// \author Christian Finck
///
///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to MUON digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TClonesArray.h>

class AliRawReader;
class AliMUONDDLTrigger;
class AliMUONLocalStruct;
class AliMUONRegHeader;

class AliMUONRawStreamTrigger: public TObject {
  public :
    AliMUONRawStreamTrigger();
    AliMUONRawStreamTrigger(AliRawReader* rawReader);
    AliMUONRawStreamTrigger(const AliMUONRawStreamTrigger& stream);
    AliMUONRawStreamTrigger& operator = (const AliMUONRawStreamTrigger& stream);
    virtual ~AliMUONRawStreamTrigger();

    virtual Bool_t   Next();
    virtual Bool_t   NextDDL();
    virtual void     ResetDDL();

    Int_t GetMaxDDL() const {return fMaxDDL;}
    Int_t GetMaxReg() const {return fMaxReg;}
    Int_t GetMaxLoc() const {return fMaxLoc;}


    void SetMaxDDL(Int_t ddl);
    void SetMaxReg(Int_t reg);
    void SetMaxLoc(Int_t loc);


    void SetReader(AliRawReader* rawReader) {fRawReader = rawReader;}

    AliMUONLocalStruct*     GetLocalInfo()  const {return fLocalStruct;}
    AliMUONDDLTrigger*      GetDDLTrigger() const {return fDDLTrigger;}
    Int_t                   GetDDL()        const {return fDDL - 1;}

  protected :

    AliRawReader*    fRawReader;    ///< object for reading the raw data
 
    Int_t  fDDL;          ///< number of DDL
    Int_t  fSubEntries;   ///< entries of buspatch structure
    Bool_t fNextDDL;      ///< flag for next DDL to be read

    Int_t fMaxDDL;        ///< maximum number of DDL in DATE file
    Int_t fMaxReg;        ///< maximum number of regional cards in DATE file
    Int_t fMaxLoc;        ///< maximum number of local cards in DATE file

    AliMUONDDLTrigger*       fDDLTrigger;   //!< pointer for DDL structure
    AliMUONRegHeader*        fRegHeader;    //!< pointer for regional structure
    AliMUONLocalStruct*      fLocalStruct;  //!< pointer to local structure

    ClassDef(AliMUONRawStreamTrigger, 2)    // base class for reading MUON trigger rawdata
};

#endif
