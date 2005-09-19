/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifTypes.h,v 1.7 2005/08/26 15:43:36 ivana Exp $

/// \ingroup motif
/// AliMpMotifTypes
/// Sytem dependent types definitions for motif category.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_TYPES_H
#define ALI_MP_MOTIF_TYPES_H

#include "AliMpContainers.h"

#ifdef WITH_STL
  #include <map>
  #include <vector>
#endif

#ifdef WITH_ROOT
  #include <TExMap.h>
  #include <TObjArray.h>
#endif

#include <TVector2.h>
#include <TString.h>

#include "AliMpIntPair.h"

class AliMpConnection;
class AliMpVMotif;
class AliMpMotifType;
class AliMpMotifPosition;

#ifdef WITH_STL
#ifdef __HP_aCC
  typedef map<AliMpIntPair, AliMpConnection*> ConnectionMap_t;
  typedef ConnectionMap_t::const_iterator     ConnectionMapCIterator;
  typedef map<TString, AliMpVMotif*> MotifMap;
  typedef MotifMap::const_iterator   MotifMapIterator;
  typedef map<TString, AliMpMotifType*>  MotifTypeMap;
  typedef MotifTypeMap::const_iterator   MotifTypeMapIterator;
  typedef map<Int_t, AliMpMotifPosition*>  MotiPositionMap;
  typedef MotiPositionMap::const_iterator  MotifPositionMapIterator;
  typedef map<AliMpIntPair, AliMpMotifPosition*> MotifPositionMap2;
  typedef MotifPositionMap2::const_iterator      MotifPositionMap2Iterator;
  typedef map<string,pair<Int_t,Int_t> > PadMapType;
  typedef PadMapType::iterator PadMapTypeIterator;
  typedef vector<TVector2> DimensionsMap;
#else
  typedef std::map< AliMpIntPair, AliMpConnection* > ConnectionMap_t;
  typedef ConnectionMap_t::const_iterator     ConnectionMapCIterator;
  typedef std::map<TString, AliMpVMotif*> MotifMap;
  typedef MotifMap::const_iterator        MotifMapIterator;
  typedef std::map<TString, AliMpMotifType*> MotifTypeMap;
  typedef MotifTypeMap::const_iterator       MotifTypeMapIterator;
  typedef std::map<Int_t, AliMpMotifPosition*>  MotiPositionMap;
  typedef MotiPositionMap::const_iterator       MotifPositionMapIterator;
  typedef std::map<AliMpIntPair, AliMpMotifPosition*> MotifPositionMap2;
  typedef MotifPositionMap2::const_iterator           MotifPositionMap2Iterator;
  typedef std::map<std::string, std::pair<Int_t,Int_t> > PadMapType;
  typedef PadMapType::iterator PadMapTypeIterator;
  typedef std::vector< TVector2 > DimensionsMap;
#endif
#endif

#ifdef WITH_ROOT
  typedef TExMap     ConnectionMap_t;
  typedef TExMapIter ConnectionMapCIterator;
  typedef TExMap     MotifMap;
  typedef TExMapIter MotifMapIterator;
  typedef TExMap     MotifTypeMap;
  typedef TExMapIter MotifTypeMapIterator;
  typedef TExMap     MotifPositionMap;
  typedef TExMapIter MotifPositionMapIterator;
  typedef TExMap     MotifPositionMap2;
  typedef TExMapIter MotifPositionMap2Iterator;
  typedef TExMap     PadMapType;
  typedef TExMapIter PadMapTypeIterator;
  typedef TObjArray  DimensionsMap;
#endif

#endif //ALI_MP_MOTIF_TYPES_H
