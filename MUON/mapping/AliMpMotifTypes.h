// $Id$
// Category: motif
//
// AliMpMotifTypes
// ---------------
// Sytem dependent types definitions for motif category.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_TYPES_H
#define ALI_MP_MOTIF_TYPES_H

#include <map>
#include <vector>

#include <TVector2.h>
#include <TString.h>

#include "AliMpIntPair.h"

class AliMpConnection;
class AliMpVMotif;
class AliMpMotifType;
class AliMpMotifPosition;

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
  typedef std::vector< TVector2 > DimensionsMap;
#endif

#endif //ALI_MP_MOTIF_TYPES_H
