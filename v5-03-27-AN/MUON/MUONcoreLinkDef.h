/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \file MUONcoreLinkDef.h
/// \brief The CINT link definitions for \ref core 

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class  AliMpIntPair+;
#pragma link C++ class  AliMpExMap-;
#pragma link C++ class  AliMpExMapIterator+;
#pragma link C++ class  AliMpArrayI+;
#pragma link C++ class  AliMpStringObjMap+;

#pragma link C++ class AliMUONObjectPair+;
#pragma link C++ class AliMUONStringIntMap+;
#pragma link C++ class AliMUON2DMap+;
#pragma link C++ class AliMUON2DMapIterator+;
#pragma link C++ class AliMUON2DMapIteratorByI+;
#pragma link C++ class AliMUON1DArray+;
#pragma link C++ class AliMUON1DMap+;
#pragma link C++ class AliMUONCheckItem+;
#pragma link C++ class AliMUONVStore+;
#pragma link C++ class AliMUONTreeManager+;
#pragma link C++ class AliMUONLogger+;

#pragma link C++ function operator-(const AliMpIntPair& ,const AliMpIntPair& );
#pragma link C++ function operator+(const AliMpIntPair& ,const AliMpIntPair& );
#pragma link C++ function operator<<(ostream& ,const AliMpIntPair& );

#endif


