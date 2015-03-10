/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSt345Reader.h,v 1.6 2006/05/23 13:07:47 ivana Exp $

/// \ingroup slat
/// \class AliMpSt345Reader
/// \brief Read slat and pcb ASCII files
/// 
/// \author Laurent Aphecetche

#ifndef ALI_MP_ST345_READER_H
#define ALI_MP_ST345_READER_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ROOT_Tmap
#  include "TMap.h"
#endif

#ifndef ROOT_TExMap
#  include "TExMap.h"
#endif

#ifndef ALI_MP_PLANE_TYPE_H
#  include "AliMpPlaneType.h"
#endif

class AliMpSlatMotifMap;
class AliMpSlat;
class AliMpPCB;
class AliMpDataStreams;

class AliMpSt345Reader : public TObject
{
 public:
  AliMpSt345Reader(AliMpSlatMotifMap* motifMap);
  virtual ~AliMpSt345Reader();

  AliMpSlat* ReadSlat(const AliMpDataStreams& dataStreams,
                      const char* slatType, AliMp::PlaneType planeType);

  AliMpPCB* ReadPCB(const AliMpDataStreams& dataStreams,
                    const char* pcbType);

private:
  /// Not implemented
  AliMpSt345Reader();
  /// Not implemented
  AliMpSt345Reader(const AliMpSt345Reader& rhs);
  /// Not implemented
  AliMpSt345Reader& operator=(const AliMpSt345Reader& rhs);

  AliMpSlatMotifMap* fMotifMap; //!<! storage for motifTypes and motifs...
  
  ClassDef(AliMpSt345Reader,0) // Reader for slat stations mapping files 
};

#endif
