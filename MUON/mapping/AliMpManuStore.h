/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpManuStore.h,v 1.6 2006/05/24 13:58:16 ivana Exp $ 

/// \ingroup management
/// \class AliMpManuStore
/// \brief The container class for manu serial numbers
///
/// \author Ivana Hrivnacova, IPN Orsay; Christian Finck, SUBATECH Nantes

#ifndef ALI_MP_MANU_STORE_H
#define ALI_MP_MANU_STORE_H

#include <TObject.h>

#include "AliMpPlaneType.h"
#include "AliMpStationType.h"
#include "AliMpIntPair.h"

#include <TString.h>
#include <TExMap.h>

class AliMpDetElement;
class AliMpDataStreams;
class TString;

class AliMpManuStore : public  TObject {

  public:
    AliMpManuStore(TRootIOCtor* ioCtor);
    virtual ~AliMpManuStore();
    
    // static access method
    static AliMpManuStore* Instance(Bool_t warn = true); 
    static AliMpManuStore* ReadData(const AliMpDataStreams& dataStreams,
                                    Bool_t warn = true);
                                    
    static void SetWarnIfDoublon(Bool_t warn);                             
    
    
    // methods
    Bool_t  AddManu(Int_t detElemId, Int_t manuId, Int_t serialNb);

    Int_t  NofManus() const;
    Int_t  NofManus(Int_t detElemId) const;

    Int_t  GetManuSerial(Int_t detElemId, Int_t manuId) const;
    AliMpIntPair  GetDetElemIdManu(Int_t manuSerial) const;

  private:
    AliMpManuStore(const AliMpDataStreams& dataStreams);
     /// Not implemented
    AliMpManuStore();
     /// Not implemented
    AliMpManuStore(const AliMpManuStore& rhs);
    /// Not implemented
    AliMpManuStore& operator=(const AliMpManuStore& rhs);
 
    // methods
    Bool_t ReadData(const AliMpDetElement* detElement, Int_t& nofManus);
    Bool_t ReadManuSerial();
    
    // not yet in use methods
    void   ReplaceManu(Int_t detElemId, Int_t manuId, Int_t serialNb);
    Bool_t WriteData(const TString& outDir = "data_run_out");

    // static data members	
    static AliMpManuStore* fgInstance;      ///< Singleton instance
    static Bool_t          fgWarnIfDoublon; ///< Option to warn about doublons

    // data members	
    const AliMpDataStreams& fkDataStreams; //!< Data streams
    mutable TExMap fManuToSerialNbs; ///< Map from manuId to serial #   
    mutable TExMap fSerialNbToManus; ///< Map manu serial # to manuId
    mutable TExMap fNofManusInDE;    ///< Number of manus with serial nbs in DE
    Int_t          fNofManus;        ///< Total number of manus
      
  ClassDef(AliMpManuStore,1)  // The manager class for definition of detection element types
};

// inline functions

inline void AliMpManuStore::SetWarnIfDoublon(Bool_t warn) 
{ 
/// Set option to warn if the same serial number is present for more manus

  fgWarnIfDoublon = warn; 
}                                   
    
    
#endif //ALI_MP_MANU_STORE_H















