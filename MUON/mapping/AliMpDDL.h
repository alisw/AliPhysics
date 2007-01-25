/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpDDL.h,v 1.6 2006/05/24 13:58:16 ivana Exp $ 

/// \ingroup management
/// \class AliMpDDL
/// \brief The class defined electronics properties of DDL
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MP_DDL_H
#define ALI_MP_DDL_H

#include <TObject.h>

#include "AliMpArrayI.h"

class AliMpDDL : public  TObject {

  public:
    AliMpDDL(Int_t id);
    AliMpDDL(TRootIOCtor* /*ioCtor*/);
    virtual ~AliMpDDL();

    // methods 
    Bool_t AddDE(Int_t detElemId);
    void   FillBusPatchIds();

    // get methods
    Int_t  GetId() const;
    
    // DEs
    Int_t  GetNofDEs() const;
    Int_t  GetDEId(Int_t index) const;
    Bool_t HasDEId(Int_t detElemId) const;
    
    // Bus patches
    Int_t  GetNofBusPatches() const;
    Int_t  GetBusPatchId(Int_t index) const;
    Bool_t HasBusPatchId(Int_t busPatchId) const;
    
    // Dsp info
    Int_t  GetMaxDsp() const;
    void   GetBusPerDsp(Int_t* iBusPerDSP) const; 

  private:
    AliMpDDL();
    AliMpDDL(const AliMpDDL& rhs);
    AliMpDDL& operator=(const AliMpDDL& rhs);

    // data members	
    Int_t       fId;          ///< Identifier (unique)
    AliMpArrayI fDEIds;       ///< Detection element Ids connected to this DDL
    AliMpArrayI fBusPatchIds; ///< Bus patch Ids connected to this DDL
    
     
  ClassDef(AliMpDDL,1)  // The class collectiong electronics properties of DDL
};

// inline functions

/// Return the unique Id
inline Int_t AliMpDDL::GetId() const
{  return fId; }

#endif //ALI_MP_MANAGER_H















