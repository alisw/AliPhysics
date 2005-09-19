/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSectorTypes.h,v 1.7 2005/08/26 15:43:36 ivana Exp $

/// \ingroup sector
/// AliMpSectorTypes
/// System dependent types definitions for sector category.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_SECTOR_TYPES_H
#define ALI_MP_SECTOR_TYPES_H

#include "AliMpContainers.h"

#ifdef WITH_STL
  #include <vector>
  #include <map>
  #include <set>
#endif

#ifdef WITH_ROOT
  #include <TArrayI.h>
  #include <TObjArray.h>
  #include <TList.h>
  #include <TExMap.h>
#endif

#include <string>

#include <TVector2.h>

#include "AliMpPad.h"

class AliMpVPadRowSegment;
class AliMpPadRow;
class AliMpVMotif;
class AliMpVRowSegment;
class AliMpSubZone;
class AliMpZone;
class AliMpRow;

#ifdef WITH_STL
#ifdef __HP_aCC
  typedef vector<Int_t> MotifPositionIdVector;
  typedef vector<AliMpPad> PadVector;
  typedef vector<AliMpPadRow*>  PadRowVector;
  typedef vector<AliMpVMotif*>  MotifVector;
  typedef vector<AliMpVRowSegment*>  RowSegmentVector;
  typedef vector<AliMpVPadRowSegment*>  PadRowSegmentVector;
  typedef vector<AliMpSubZone*> SubZoneVector;
  typedef vector<AliMpRow*> RowVector;
  typedef vector<AliMpZone*> ZoneVector;
  typedef map<Int_t, TVector2>  PadDimensionsMap;
  typedef PadDimensionsMap::const_iterator  PadDimensionsMapCIterator;
  typedef set<AliMpPad> PadSet;
  typedef PadSet::const_iterator PadSetIterator;
#else
  typedef std::vector<Int_t> MotifPositionIdVector;
  typedef std::vector<AliMpPad> PadVector;
  typedef std::vector<AliMpPadRow*>  PadRowVector;
  typedef std::vector<AliMpVMotif*>  MotifVector;
  typedef std::vector<AliMpVRowSegment*>  RowSegmentVector;
  typedef std::vector<AliMpVPadRowSegment*>  PadRowSegmentVector;
  typedef std::vector<AliMpSubZone*>  SubZoneVector;
  typedef std::vector<AliMpRow*> RowVector;
  typedef std::vector<AliMpZone*> ZoneVector;
  typedef std::map<Int_t, TVector2>  PadDimensionsMap;
  typedef PadDimensionsMap::const_iterator  PadDimensionsMapCIterator;
  typedef std::set<AliMpPad> PadSet;
  typedef PadSet::const_iterator PadSetIterator;
#endif
#endif

#ifdef WITH_ROOT
#ifndef __HP_aCC
  using std::string;
#endif
  typedef TArrayI    MotifPositionIdVector;
  typedef TObjArray  PadVector;
  typedef TObjArray  PadRowVector;
  typedef TObjArray  MotifVector;
  typedef TList      RowSegmentVector;
  typedef TObjArray  PadRowSegmentVector;
  typedef TObjArray  SubZoneVector;
  typedef TObjArray  RowVector;
  typedef TObjArray  ZoneVector;
  typedef TExMap     PadDimensionsMap;
  typedef TExMapIter PadDimensionsMapCIterator;
  typedef TObjArray  PadSet;
#endif

#endif //ALI_MP_SECTOR_TYPES_H
