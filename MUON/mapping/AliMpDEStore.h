/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpDEStore.h,v 1.6 2006/05/24 13:58:16 ivana Exp $ 

/// \ingroup management
/// \class AliMpDEStore
/// \brief The container class for detection element objects
///
/// \author Ivana Hrivnacova, IPN Orsay;
///         Laurent Aphecetche, Christian Finck, SUBATECH Nantes

#ifndef ALI_MP_DE_STORE_H
#define ALI_MP_DE_STORE_H

#include <TObject.h>
#include <TArrayI.h>

#include "AliMpExMap.h"
#include "AliMpPlaneType.h"
#include "AliMpStationType.h"
#include "AliMpIntPair.h"

class AliMpDetElement;

class AliMpDEStore : public  TObject {

  friend class AliMpDEIterator;

  public:
    AliMpDEStore(TRootIOCtor* /*ioCtor*/);
    virtual ~AliMpDEStore();
    
    // static access method
    static AliMpDEStore* Instance(); 
    
    // methods
    AliMpDetElement* GetDetElement(Int_t detElemId, Bool_t warn = true) const;
    AliMpIntPair     GetDetElemIdManu(Int_t manuSerial) const;
    
  private:
    AliMpDEStore();
    AliMpDEStore(const AliMpDEStore& rhs);
    AliMpDEStore& operator=(const AliMpDEStore& rhs);

    // methods
    Bool_t IsPlaneType(const TString& planeTypeName);
    AliMp::PlaneType   PlaneType(const TString& planeTypeName);
    AliMp::StationType StationType(const TString& stationTypeName);

    Bool_t ReadManuToSerialNbs(AliMpDetElement* detElement, 
                       AliMp::StationType stationType);
    Bool_t ReadDENames(AliMp::StationType stationType);
    void   FillDEs();

    // static data members	
    static AliMpDEStore* fgInstance;       ///< Singleton instance
    static const char    fgkCommentPrefix; ///< Comment prefix in DE names file

    // data members	
    AliMpExMap fDetElements; ///< Map between DE Ids and DE objects
      
  ClassDef(AliMpDEStore,1)  // The manager class for definition of detection element types
};

#endif //ALI_MP_MANAGER_H















