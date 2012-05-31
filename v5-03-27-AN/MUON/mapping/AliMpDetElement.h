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
#include "AliMpPlaneType.h"
#include "AliMpCathodType.h"
#include "AliMpStationType.h"
#include "AliMpStation12Type.h"

#ifndef ALI_MP_EX_MAP_H
#  include "AliMpExMap.h"
#endif

class AliMpVSegmentation;
class AliMpArrayI;

class AliMpDetElement : public  TObject {

  public:  
    AliMpDetElement(Int_t id, const TString& name,
                    const TString& segType, AliMp::PlaneType planeType);
    AliMpDetElement(TRootIOCtor* ioCtor);
    virtual ~AliMpDetElement();

    // static methods
    static char GetNameSeparator(); 
    
    // methods
    Bool_t AddBusPatch(Int_t busPatchId); 
    void   AddManu(Int_t manuId);
    void   SetDdlId(Int_t ddlId);

    // get methods
    Int_t   GetId() const;
    Int_t   GetDdlId() const;
    TString GetDEName() const;
    TString GetSegType() const;
    TString GetSegName(AliMp::CathodType cath) const;

    AliMp::PlaneType     GetPlaneType(AliMp::CathodType cath) const;
    AliMp::CathodType    GetCathodType(AliMp::PlaneType planeType) const;
    AliMp::CathodType    GetCathodTypeFromManuId(Int_t manuId) const;
    AliMp::StationType   GetStationType() const;
    AliMq::Station12Type GetStation12Type() const;
    
    Int_t  GetNofBusPatches() const;
    Int_t  GetBusPatchId(Int_t index) const;
    Bool_t HasBusPatchId(Int_t busPatchId) const;

    Int_t  NofManus() const;
    Int_t  NofChannelsInManu(Int_t manuId) const;
    Bool_t IsExistingChannel(Int_t manuId, Int_t manuChannel) const;
    Bool_t IsConnectedChannel(Int_t manuId, Int_t manuChannel) const;
    
    const AliMpArrayI* ManusForHV(Int_t hvIndex) const;
    
           /// Return the number of channels in this detection element    
    Int_t NofChannels() const { return fNofChannels; }
    
  private:
    /// Not implemented
    AliMpDetElement();
    /// Not implemented
    AliMpDetElement(const AliMpDetElement& rhs);
    /// Not implemented
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
    
    mutable TExMap fManuList;  ///< map of manus
    mutable TExMap fTrackerChannels; ///< list of connected pads (tracker only)
    
    AliMpExMap fHVmanus; ///< map of HV->manu
    
    Int_t fNofChannels; ///< number of channels in this detection element
    
  ClassDef(AliMpDetElement,4)  // The manager class for definition of detection element types
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















