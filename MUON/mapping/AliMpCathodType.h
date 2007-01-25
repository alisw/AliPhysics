/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpCathodType.h,v 1.8 2006/05/24 13:58:07 ivana Exp $

/// \ingroup basic
/// \enum AliMpCathodType
/// Enumeration for refering to cath0 and cath1.
///
/// \author Ivana Hrivnacova; IPN Orsay
 
#ifndef ALI_MP_CATHOD_TYPE_H
#define ALI_MP_CATHOD_TYPE_H

#include <Rtypes.h>
#include <TString.h>

namespace AliMp {

  enum CathodType
  {
    kCath0, ///< cathod 0
    kCath1  ///< cathod 1
  };

  /// Convert integer number in enum;
  AliMp::CathodType GetCathodType(Int_t cathodNumber);

  /// Return name for given cathodType
  TString CathodTypeName(AliMp::CathodType cathodType);

  /// Return the other cathod type
  AliMp::CathodType OtherCathodType(AliMp::CathodType cathodType);

} 
      
#endif //ALI_MP_CATHOD_TYPE_H
