/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifType.h,v 1.11 2006/05/24 13:58:18 ivana Exp $

/// \ingroup motif
/// \class AliMpMotifType
/// \brief Class that defines the motif properties.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_TYPE_H
#define ALI_MP_MOTIF_TYPE_H

#include <TObject.h>

#include "AliMpEncodePair.h"

#ifndef ROOT_TObjArray
#  include <TObjArray.h>
#endif

#include <TString.h>

class AliMpVPadIterator;
class AliMpConnection;

class AliMpMotifType : public TObject
  {
  public:
    AliMpMotifType(const TString &id);
    AliMpMotifType(const AliMpMotifType& rhs);
    AliMpMotifType& operator=(const AliMpMotifType& rhs);
    AliMpMotifType(TRootIOCtor* ioCtor);
    virtual ~AliMpMotifType();
    
    TObject* Clone(const char* newname="") const;
    
    virtual AliMpVPadIterator* CreateIterator() const;
    
    // find methods
    AliMpConnection *FindConnectionByPadNum(Int_t padNum) const;
    AliMpConnection *FindConnectionByLocalIndices(
                         MpPair_t localIndices) const;
    AliMpConnection *FindConnectionByLocalIndices(
                         Int_t localIx, Int_t localIy) const;
    AliMpConnection *FindConnectionByGassiNum(Int_t gassiNum) const;
    AliMpConnection *FindConnectionByKaptonNum(Int_t kaptonNum) const;
    AliMpConnection *FindConnectionByBergNum(Int_t bergNum) const;
    
    MpPair_t FindLocalIndicesByPadNum(Int_t padNum) const;
    MpPair_t FindLocalIndicesByGassiNum(Int_t gassiNum) const;
    MpPair_t FindLocalIndicesByKaptonNum(Int_t kaptonNum) const;
    MpPair_t FindLocalIndicesByBergNum(Int_t bergNum) const;
    MpPair_t FindLocalIndicesByConnection(
                         const AliMpConnection* connection) const;
    
    // set methods
    void SetNofPads(Int_t nofPadsX, Int_t nofPadY);
    
    // get methods
    /// Return unique motif ID
    TString  GetID() const        {return fID;}
    /// Return number of pads in x direction
    Int_t    GetNofPadsX() const  {return fNofPadsX;}
    /// Return number of pads in y direction
    Int_t    GetNofPadsY() const  {return fNofPadsY;}
    
    Int_t    GetNofPads() const   {return fNofPads;}
    
    // Other methods
    Bool_t AddConnection(AliMpConnection* connection);

    virtual void Print(Option_t *option="") const;
    
    Int_t   PadNum(const TString &padName) const;
    
    TString PadName(Int_t padNum) const;

    Bool_t HasPadByLocalIndices(MpPair_t localIndices) const;
    Bool_t HasPadByLocalIndices(Int_t localIx, Int_t localIy) const;

    Bool_t HasPadByManuChannel(Int_t manuChannel) const;

    Bool_t HasPadByGassiNum(Int_t gassiNum) const { return HasPadByManuChannel(gassiNum); }
    
    Bool_t IsFull() const;
    
    Bool_t Save(const char* motifName) const;
    Bool_t Save() const;
    
  private:
    /// Not implemented
    AliMpMotifType();

    // methods
    void Copy(TObject& o) const;

    // static data members
    static const Int_t  fgkPadNumForA; ///< the pad number for the pad "A"
    
    // data members
    TString   fID;              ///< unique motif ID
    Int_t     fNofPadsX;        ///< number of pads in x direction
    Int_t     fNofPadsY;        ///< number of pads in y direction
    Int_t     fNofPads;    ///< total number of pads (= the number of non-void entries in the arrays below)
    Int_t     fMaxNofPads; ///< max number of pads we can hold
    TObjArray fConnectionsByLocalIndices; ///< array [ix + 64*iy ] -> AliMpConnection*
    TObjArray fConnectionsByManuChannel;  ///< array [manuChannel] -> AliMpConnection*
    
    ClassDef(AliMpMotifType,2)  // Motif type
  };

// inline functions

/// Return true if the motif conatins all pads
inline Bool_t AliMpMotifType::IsFull() const 
{ return GetNofPads() == fNofPadsX*fNofPadsY; }

#endif //ALI_MP_MOTIF_TYPE_H


