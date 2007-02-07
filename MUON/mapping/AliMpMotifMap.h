/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifMap.h,v 1.14 2006/05/24 13:58:18 ivana Exp $

/// \ingroup motif
/// \class AliMpMotifMap
/// \brief Motif map containers
///
/// The class defines:
/// - map of motif objects to their string IDs
/// - map of motif type objects to their string IDs
/// - map of motif position objects to their string IDs
/// - map of motif position objects to their global indices
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_MAP_H
#define ALI_MP_MOTIF_MAP_H

#include <TObject.h>

#include "AliMpContainers.h"

#ifdef WITH_STL
#include "AliMpIntPair.h"
#endif

#ifdef WITH_ROOT
#include "AliMpExMap.h"
#endif

#ifdef WITH_STL
#include <map>
#include <iterator>
#endif

class AliMpVMotif;
class AliMpMotifType;
class AliMpMotifPosition;
class AliMpMotifMap;

class TArrayI;
class TString;
class TVector2;

class AliMpMotifMap : public TObject
{
  public:
#ifdef WITH_STL
    typedef std::map<TString, AliMpVMotif*> MotifMap;
    typedef MotifMap::const_iterator        MotifMapIterator;
    typedef std::map<TString, AliMpMotifType*> MotifTypeMap;
    typedef MotifTypeMap::const_iterator       MotifTypeMapIterator;
    typedef std::map<Int_t, AliMpMotifPosition*>  MotiPositionMap;
    typedef MotiPositionMap::const_iterator       MotifPositionMapIterator;
    typedef std::map<AliMpIntPair, AliMpMotifPosition*> MotifPositionMap2;
    typedef MotifPositionMap2::const_iterator           MotifPositionMap2Iterator;
#endif    
#ifdef WITH_ROOT
    typedef AliMpExMap MotifMap;
    typedef AliMpExMap MotifTypeMap;
    typedef AliMpExMap MotifPositionMap;
    typedef AliMpExMap MotifPositionMap2;
#endif    

  public:
    AliMpMotifMap(Bool_t /*standardConstructor*/);
    AliMpMotifMap();
    virtual ~AliMpMotifMap();
    
    // methods
    Bool_t  AddMotif(AliMpVMotif* motif, Bool_t warn = true);
    Bool_t  AddMotifType(AliMpMotifType* motifType, Bool_t warn = true);
    Bool_t  AddMotifPosition(AliMpMotifPosition* motifType, Bool_t warn = true);
    void   FillMotifPositionMap2();
    virtual void Print(const char* option = "ALL") const;
    void   PrintGlobalIndices(const char* fileName) const;
    void   UpdateGlobalIndices(const char* fileName);
   
    // find methods
    AliMpVMotif*  FindMotif(const TString& motifID) const;
    AliMpVMotif*  FindMotif(const TString& motifID, const TString& motifTypeID, 
                            const TVector2& padDimensions) const;
    AliMpMotifType*      FindMotifType(const TString& motifTypeID) const;
    AliMpMotifPosition*  FindMotifPosition(Int_t motifPositionID) const;

    /// Find all motifPositionsIDs (=electronicCardNumbers) handled by this map
    void    GetAllMotifPositionsIDs(TArrayI& enc) const;
    UInt_t  GetNofMotifPositions() const;
    AliMpMotifPosition* GetMotifPosition(UInt_t index) const;

    /// Calculate total number of pads defined in the map
    Int_t CalculateNofPads() const;
     
  private:
    // methods
    void  PrintMotif(const AliMpVMotif* motif) const;
    void  PrintMotifType(const AliMpMotifType* motifType) const;
    void  PrintMotifPosition(const AliMpMotifPosition* motifPosition) const;
    void  PrintMotifPosition2(const AliMpMotifPosition* motifPosition) const;
    void  PrintMotifs() const;
    void  PrintMotifTypes() const;
    void  PrintMotifPositions() const;
    void  PrintMotifPositions2() const;
 
    // data members
    MotifMap           fMotifs;         ///< motifs map
    MotifTypeMap       fMotifTypes;     ///< motifs types map
#ifdef WITH_STL
    std::map<Int_t, AliMpMotifPosition*> fMotifPositions; ///< motif positions map by Id
#endif
#ifdef WITH_ROOT
    MotifPositionMap   fMotifPositions; ///< motifs positions map
#endif
    MotifPositionMap2  fMotifPositions2;///< motifs positions map

  ClassDef(AliMpMotifMap,1)  // motif map
};

#endif //ALI_MP_MOTIF_MAP_H
