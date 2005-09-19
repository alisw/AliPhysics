/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$ 
// $MpId: AliMpHelper.h,v 1.4 2005/09/19 19:01:09 ivana Exp $ 

/// \ingroup slat
/// \class AliMpHelper
/// \brief Helper class to parse slat mapping ascii files.
/// 
/// Author: Laurent Aphecetche

#ifndef ALI_MP_HELPER_H
#define ALI_MP_HELPER_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class TArrayI;
class TString;

class AliMpHelper : public TObject
{
 public:
  AliMpHelper();
  virtual ~AliMpHelper();

  static void DecodeName(const char* manus, char sep, TArrayI& theList);                       

  static void GetRange(const char* str, Int_t& begin, Int_t& end, 
		                   Int_t& incr, Int_t& n);

  static TString Normalize(const char* line);
                           
  ClassDef(AliMpHelper,1) // Helper for parsing slat stations mapping files 
};

#endif
