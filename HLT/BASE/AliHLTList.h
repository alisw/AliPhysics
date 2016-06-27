/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id: AliHLTList.h 64008 2016-08-25 13:09:59Z mkrzewic $ */

//-------------------------------------------------------------------------
//                          Class AliHLTList
//essentially a TList, becomes owner after deserializing
// Origin: Mikolaj Krzewicki, mkrzewic@cern.ch
//-------------------------------------------------------------------------

#ifndef ALIHLTLIST_H
#define ALIHLTLIST_H
#include "TList.h"
class AliHLTList : public TList {
public:
AliHLTList():TList() {}
virtual ~AliHLTList() {}
ClassDef(AliHLTList,1)
};
#endif

