/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpManuGeo.h,v 1.5 2006/05/24 13:58:16 ivana Exp $

/// \ingroup management
/// \class AliMpManuGeo
/// \brief Class that manages the maps manuId<>manuSerial#<>DE 
///
/// \author Ch. Finck; Subatech Nantes

#ifndef ALI_MP_MANUGEO_H
#define ALI_MP_MANUGEO_H

#include <TObject.h>

#include <TExMap.h>
#include <TArrayI.h>

class AliMpIntPair;

class AliMpManuGeo : public TObject
{

 public:

  AliMpManuGeo();
  virtual ~AliMpManuGeo();
   

  // methods
  void ReadGeomManuFiles();
  void ReadGeomManuFile(Int_t idDE);

  AliMpIntPair  GetDetElemManu(Int_t manuSerial);
  Int_t         GetManuSerial(AliMpIntPair& pair);

 private:
  AliMpManuGeo(const AliMpManuGeo& src);
  AliMpManuGeo& operator = (const AliMpManuGeo& src) ;


  TExMap fDeManuToSerialNb; //!< Map from (idDE, manuId) to manu serial #   
  TExMap fSerialNbToDeManu; //!< Map manu serial # to (idDE, manuId)


  ClassDef(AliMpManuGeo,1) //utility class for the motif type
};


#endif 
