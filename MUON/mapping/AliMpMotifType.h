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

#include "AliMpContainers.h"

#include "AliMpIntPair.h"
#ifdef WITH_ROOT
#include "AliMpExMap.h"
#endif

#include <TString.h>

#ifdef WITH_STL
#include <map>
#endif

class AliMpConnection;
class AliMpVPadIterator;

class AliMpMotifType : public TObject
{
  public:
#ifdef WITH_STL
    typedef std::map< AliMpIntPair, AliMpConnection* > ConnectionMap;
    typedef ConnectionMap::const_iterator     ConnectionMapCIterator;
#endif    
#ifdef WITH_ROOT
    typedef AliMpExMap ConnectionMap;
#endif    

  public:
    AliMpMotifType(const TString &id);
    AliMpMotifType();
    virtual ~AliMpMotifType();

    virtual AliMpVPadIterator* CreateIterator() const;

    // find methods
    AliMpConnection *FindConnectionByPadNum(Int_t padNum) const;
    AliMpConnection *FindConnectionByLocalIndices(
                               const AliMpIntPair& localIndices) const;
    AliMpConnection *FindConnectionByGassiNum(Int_t gassiNum) const;
    AliMpConnection *FindConnectionByKaptonNum(Int_t kaptonNum) const;
    AliMpConnection *FindConnectionByBergNum(Int_t bergNum) const;

    AliMpIntPair FindLocalIndicesByPadNum(Int_t padNum) const;
    AliMpIntPair FindLocalIndicesByGassiNum(Int_t gassiNum) const;
    AliMpIntPair FindLocalIndicesByKaptonNum(Int_t kaptonNum) const;
    AliMpIntPair FindLocalIndicesByBergNum(Int_t bergNum) const;
    AliMpIntPair FindLocalIndicesByConnection(
                               const AliMpConnection* connection) const;

    // set methods
    void SetNofPads(Int_t nofPadsX, Int_t nofPadY);
    void SetVerboseLevel(Int_t level){fVerboseLevel=level;}
    
    // get methods
    TString  GetID() const        {return fID;}
    Int_t    GetNofPadsX() const  {return fNofPadsX;}
    Int_t    GetNofPadsY() const  {return fNofPadsY;}
    Int_t    GetNofPads() const;
    
    // Other methods
    void AddConnection(const AliMpIntPair &localIndices, 
                       AliMpConnection* connection);
    virtual void Print(Option_t *option="") const;
    Int_t   PadNum(const TString &padName) const;
    TString PadName(Int_t padNum) const;
    Bool_t  HasPad(const AliMpIntPair& localIndices) const;
    Bool_t  IsFull() const;

  private:
    // static data members
    static const Int_t  fgkPadNumForA; ///< the pad number for the pad "A"
  
    // data members
    TString   fID;              ///< unique motif ID
    Int_t     fNofPadsX;        ///< number of pads in x direction
    Int_t     fNofPadsY;        ///< number of pads in y direction
    Int_t     fVerboseLevel;    ///< verbose level
    ConnectionMap fConnections; ///< Map (ix,iy) of connections
    
  ClassDef(AliMpMotifType,1)  // Motif type
};

// inline functions

inline Bool_t AliMpMotifType::IsFull() const 
{ return GetNofPads() == fNofPadsX*fNofPadsY; }

#endif //ALI_MP_MOTIF_TYPE_H

