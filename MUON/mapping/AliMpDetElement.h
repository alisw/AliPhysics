/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpDetElement.h,v 1.6 2006/05/24 13:58:16 ivana Exp $ 

/// \ingroup management
/// \class AliMpDetElement
/// \brief The class defines the electronics properties of detection element
///
/// \author Ivana Hrivnacova, IPN Orsay;
///         Laurent Aphecetche, Ch. Finck, Subatech Nantes

#ifndef ALI_MP_DET_ELEMENT_H
#define ALI_MP_DET_ELEMENT_H

#include <TObject.h>
#include <TArrayI.h>
#include <TExMap.h>

#include "AliMpArrayI.h"
#include "AliMpStationType.h"
#include "AliMpPlaneType.h"
#include "AliMpCathodType.h"

class AliMpVSegmentation;

class AliMpDetElement : public  TObject {

  public:
    AliMpDetElement(Int_t id, const TString& name,
                    const TString& segType, AliMp::PlaneType planeType);
    AliMpDetElement(TRootIOCtor* /*ioCtor*/);
    virtual ~AliMpDetElement();

    // static methods
    static char GetNameSeparator(); 
    
    // methods
    Bool_t AddBusPatch(Int_t busPatchId); 
    void   AddManuSerial(Int_t manuId, Int_t serialNb); 
    void   SetDdlId(Int_t ddlId);

    // get methods
    Int_t   GetId() const;
    Int_t   GetDdlId() const;
    TString GetDEName() const;
    TString GetSegType() const;
    TString GetSegName(AliMp::CathodType cath) const;

    AliMp::PlaneType   GetPlaneType(AliMp::CathodType cath) const;
    AliMp::CathodType  GetCathodType(AliMp::PlaneType planeType) const;
    AliMp::StationType GetStationType() const;
    
    Int_t  GetNofBusPatches() const;
    Int_t  GetBusPatchId(Int_t index) const;
    Bool_t HasBusPatchId(Int_t busPatchId) const;
    
    Int_t  GetNofManus() const;    
    Int_t  GetManuSerialFromId(Int_t manuId) const;
    Int_t  GetManuIdFromSerial(Int_t serialNb) const;

  private:
    AliMpDetElement();
    AliMpDetElement(const AliMpDetElement& rhs);
    AliMpDetElement& operator=(const AliMpDetElement& rhs);

    // static data members	
    static const char  fgkNameSeparator; ///< Separator character used in DE names

    // data members	
    Int_t          fId;         ///< Identifier (unique)
    Int_t          fDdlId;      ///< DDL Id to which this DE is connected
    TString        fName;       ///< Name unique
    TString        fSegType;    ///< Segmentation type name
    AliMp::PlaneType fPlaneType;  ///< Plane type on cathod0
    //AliMpExMap     fBusPatches; ///< Bus patches connected to this detection element
    AliMpArrayI    fBusPatchIds;  ///< Bus patches connected to this detection element
    mutable TExMap fManuToSerialNbs; //< Map from manuId to serial #   
    mutable TExMap fSerialNbToManus; //< Map manu serial # to manuId
     
  ClassDef(AliMpDetElement,1)  // The manager class for definition of detection element types
};

// inline function

/// Return the name separator
inline  char AliMpDetElement::GetNameSeparator()
{ return fgkNameSeparator; }  

/// Set DDL Id
inline  void AliMpDetElement::SetDdlId(Int_t ddlId)
{ fDdlId = ddlId; }

/// Return Id
inline  Int_t   AliMpDetElement::GetId() const
{ return fId; }

/// Return DDL Id
inline  Int_t   AliMpDetElement::GetDdlId() const
{ return fDdlId; }

/// Return name
inline  TString AliMpDetElement::GetDEName() const
{ return fName; }

/// Return segmentation type name
inline  TString AliMpDetElement::GetSegType() const
{ return fSegType; }

#endif //ALI_MP_MANAGER_H















