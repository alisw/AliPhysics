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

    // get methods
    Int_t  GetId() const;
    Int_t  GetDEId() const;
    Int_t  GetDdlId() const;
    Int_t  GetNofManus() const;
    Int_t  GetManuId(Int_t index) const;
    Bool_t HasManu(Int_t manuId) const;
    

  private:
    AliMpBusPatch();
    AliMpBusPatch(const AliMpBusPatch& rhs);
    AliMpBusPatch& operator=(const AliMpBusPatch& rhs);

    // static data members	
    static const Int_t  fgkOffset; //< Offset for conversion global/local ID  

    // data members	
    Int_t        fId;    ///< Identifier (unique)
    Int_t        fDEId;  ///< Detection element to which this bus patch is connected
    Int_t        fDdlId; ///< DDL to which this bus patch is connected
    AliMpArrayI  fManus; ///< Manu Ids connected to this bus patch
     
  ClassDef(AliMpBusPatch,1)  // The class collectiong electronics properties of DDL
};

// inline functions

/// Return the unique Id
inline Int_t AliMpBusPatch::GetId() const
{  return fId; }

/// Return the Detection element Id
inline Int_t AliMpBusPatch::GetDEId() const
{  return fDEId; }

/// Return the Ddl  Id
inline Int_t AliMpBusPatch::GetDdlId() const
{  return fDdlId; }

#endif //ALI_BUS_PATCH_H














