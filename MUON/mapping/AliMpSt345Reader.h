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
class TObjString;

class AliMpSt345Reader : public TObject
{
 public:
  AliMpSt345Reader(AliMpSlatMotifMap&);
  virtual ~AliMpSt345Reader();

  AliMpSlat* ReadSlat(const char* slatType, AliMpPlaneType planeType);

  AliMpPCB* ReadPCB(const char* pcbType);

private:
    AliMpSlatMotifMap& fMotifMap; //!< storage for motifTypes and motifs...
  
  ClassDef(AliMpSt345Reader,0) // Reader for slat stations mapping files 
};

#endif
