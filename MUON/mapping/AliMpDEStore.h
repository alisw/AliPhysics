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
#include "AliMpStation12Type.h"
#include "AliMpIntPair.h"

class AliMpDetElement;
class AliMpDataStreams;
class TString;

class AliMpDEStore : public  TObject {

  friend class AliMpDEIterator;

  public:
    AliMpDEStore(TRootIOCtor* ioCtor);
    virtual ~AliMpDEStore();
    
    // static access method
    static AliMpDEStore* Instance(Bool_t warn = true); 
    static AliMpDEStore* ReadData(const AliMpDataStreams& dataStreams,
                                  Bool_t warn = true);
    
    // methods
    AliMpDetElement* GetDetElement(Int_t detElemId, Bool_t warn = true) const;
    AliMpDetElement* GetDetElement(const TString& detName, Bool_t warn = true) const;
    
  private:
    AliMpDEStore(const AliMpDataStreams& dataStreams);
    /// Not implemented
    AliMpDEStore();
    /// Not implemented
    AliMpDEStore(const AliMpDEStore& rhs);
    /// Not implemented
    AliMpDEStore& operator=(const AliMpDEStore& rhs);

    // methods
    Bool_t IsPlaneType(const TString& planeTypeName);
 
    Bool_t ReadDENames(AliMp::StationType stationType, 
                       AliMq::Station12Type station12Type = AliMq::kNotSt12);
    void   FillDEs();

    // static data members	
    static AliMpDEStore* fgInstance;       ///< Singleton instance
    static const char    fgkCommentPrefix; ///< Comment prefix in DE names file

    // data members	
    const AliMpDataStreams&  fkDataStreams; //!< Data streams
    AliMpExMap  fDetElements; ///< Map between DE Ids and DE objects
      
  ClassDef(AliMpDEStore,1)  // The manager class for definition of detection element types
};

#endif //ALI_MP_MANAGER_H















