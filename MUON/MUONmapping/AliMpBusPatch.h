/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: $ 

/// \ingroup management
/// \class AliMpBusPatch
/// \brief The class defines the properties of BusPatch
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MP_BUS_PATCH_H
#define ALI_MP_BUS_PATCH_H

#include <TObject.h>
#include <TString.h>

#include "AliMpArrayI.h"

class AliMpBusPatch : public  TObject {

  public:
    AliMpBusPatch(Int_t id, Int_t deId, Int_t ddlId);
    AliMpBusPatch(TRootIOCtor* /*ioCtor*/);
    virtual ~AliMpBusPatch();

    // static methods
    static Int_t GetGlobalBusID(Int_t localID, Int_t ddlID);
    static Int_t GetLocalBusID(Int_t globalID, Int_t ddlID);

    // methods 
    Bool_t AddManu(Int_t manuId);
    Bool_t SetNofManusPerModule(Int_t manuNumber = 0);
    void   SetTranslatorLabel(TString label);
    void   SetCableLabel(TString label); 
    void   SetCableLength(Float_t length);
    void   SetFrtId(Int_t id);
    void   RevertReadout();
    void   ResetReadout();
  
    // get methods
    Int_t  GetId() const;
    Int_t  GetDEId() const;
    Int_t  GetDdlId() const;
    Int_t  GetFrtId() const;
    Int_t  GetNofManus() const;
    Int_t  GetManuId(Int_t index) const;
    Bool_t HasManu(Int_t manuId) const;
    
    Int_t  GetNofPatchModules() const;
    Int_t  GetNofManusPerModule(Int_t patchModule) const;
    
    Float_t  GetCableLength() const;
    TString  GetCableLabel() const;
    TString  GetTranslatorLabel() const;
  TString GetFRTPosition() const;
  
  virtual void Print(Option_t* opt="") const;
  
  private:
    /// Not implemented
    AliMpBusPatch();
    /// Not implemented
    AliMpBusPatch(const AliMpBusPatch& rhs);
    /// Not implemented
    AliMpBusPatch& operator=(const AliMpBusPatch& rhs);

    // static data members	
    static const Int_t  fgkOffset; ///< Offset for conversion global/local ID  

    // data members	
    Int_t        fId;     ///< Identifier (unique)
    Int_t        fDEId;   ///< Detection element to which this bus patch is connected
    Int_t        fDdlId;  ///< DDL to which this bus patch is connected
    AliMpArrayI  fManus;  ///< Manu Ids connected to this bus patch
    AliMpArrayI  fNofManusPerModule; ///< Nof Manus per patch modules (PCBs)
    Float_t      fCableLength;       ///< length of the buspatch cable
    TString      fCableLabel;        ///< label of the buspatch cable
    TString      fTranslatorLabel;   ///< label of the translator board
    Int_t        fFrtId;               ///< FRT Ids connected to this bus patch

  ClassDef(AliMpBusPatch,3)  // The class collectiong electronics properties of DDL
};

// inline functions

/// Return the unique Id
inline Int_t AliMpBusPatch::GetId() const
{  return fId; }

/// Return the Detection element Id
inline Int_t AliMpBusPatch::GetDEId() const
{  return fDEId; }

/// Return the Ddl Id
inline Int_t AliMpBusPatch::GetDdlId() const
{  return fDdlId; }

/// Return the FRT Id
inline Int_t AliMpBusPatch::GetFrtId() const
{  return fFrtId; }

/// Return length of buspatch
inline Float_t  AliMpBusPatch::GetCableLength() const
{ return fCableLength; }

/// Set FRT id for buspatch
inline void  AliMpBusPatch::SetFrtId(Int_t id)
{ fFrtId = id; }

/// Set length of buspatch
inline void  AliMpBusPatch::SetCableLength(Float_t length)
{ fCableLength = length; }

/// Return label of buspatch
inline TString  AliMpBusPatch::GetCableLabel() const
{ return fCableLabel; }

/// Set label of buspatch
inline void  AliMpBusPatch::SetCableLabel(TString label)
{ fCableLabel = label; }

/// Return label of translator
inline TString  AliMpBusPatch::GetTranslatorLabel() const
{ return fCableLabel; }

/// Set label of translator
inline void  AliMpBusPatch::SetTranslatorLabel(TString label)
{ fTranslatorLabel = label; }


#endif //ALI_BUS_PATCH_H














