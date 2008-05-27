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
    /// Connection map type
    typedef std::map< AliMpIntPair, AliMpConnection* > ConnectionMap;
    /// Connection map iterator type
    typedef ConnectionMap::const_iterator     ConnectionMapCIterator;
#endif    
#ifdef WITH_ROOT
    /// Connection map type
    typedef AliMpExMap ConnectionMap;
#endif    

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
    
    // get methods
             /// Return unique motif ID
    TString  GetID() const        {return fID;}
             /// Return number of pads in x direction
    Int_t    GetNofPadsX() const  {return fNofPadsX;}
             /// Return number of pads in y direction
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
    
    Bool_t Save(const char* motifName) const;
    Bool_t Save() const;

  private:
      void Copy(TObject& o) const;
    
  private:
    /// Not implemented
    AliMpMotifType();
    // static data members
    static const Int_t  fgkPadNumForA; ///< the pad number for the pad "A"
  
    // data members
    TString   fID;              ///< unique motif ID
    Int_t     fNofPadsX;        ///< number of pads in x direction
    Int_t     fNofPadsY;        ///< number of pads in y direction
    ConnectionMap fConnections; ///< Map (ix,iy) of connections
    
  ClassDef(AliMpMotifType,1)  // Motif type
};

// inline functions

/// Return true if the motif conatins all pads
inline Bool_t AliMpMotifType::IsFull() const 
{ return GetNofPads() == fNofPadsX*fNofPadsY; }

#endif //ALI_MP_MOTIF_TYPE_H

