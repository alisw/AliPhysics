#ifndef ALI_MUON_ST1_TYPES_H
#define ALI_MUON_ST1_TYPES_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// AliMUONSt1Types
// ---------------
// System dependent types definitions for MUON Station1.
//
// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//

#include <string>
#include <map>
#include <vector>
#include <fstream>

#include "AliMUONSt1Containers.h"
#include "AliMUONSt1SpecialMotif.h"

class TString;
class AliMUONSt1ResponseParameter;

#ifdef __HP_aCC
  typedef vector<Int_t>  IntVector;
  typedef map<Int_t , AliMUONSt1SpecialMotif> TSpecialMap;
  typedef map<string,AliMUONSt1ResponseParameter*> TParamsMap;
  typedef map<string,TList*> TListMap;
#else
  using std::string;
  using std::vector;
  using std::multimap;
  using std::pair;

  typedef std::vector<Int_t>  IntVector;
  typedef std::map<Int_t , AliMUONSt1SpecialMotif> TSpecialMap;
  typedef std::map<std::string,AliMUONSt1ResponseParameter*> TParamsMap;
  typedef std::map<std::string,TList*> TListMap;
#endif

#ifdef ST1_WITH_STL
  #include <map>
  #ifdef __HP_aCC
    using std::map;
  #endif
#endif

#ifdef ST1_WITH_ROOT
  #include "TExMap.h"
#endif

#ifdef ST1_WITH_STL
    typedef map<Int_t , AliMUONSt1SpecialMotif> SpecialMap;
#endif
#ifdef ST1_WITH_ROOT
    typedef  TExMap  SpecialMap;
#endif

#endif //ALI_MUON_ST1_TYPES_H
