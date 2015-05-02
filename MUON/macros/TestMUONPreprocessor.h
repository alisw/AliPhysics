/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \file TestMUONPreprocessor.h
/// \brief The definitions for the macro TestMUONPreprocessor.C

#ifndef TESTMUONPREPROCESSOR_H
#define TESTMUONPREPROCESSOR_H

class TMap;

/// Create a fake DCS alias map
TMap* CreateDCSAliasMap(const char* inputCDB, Int_t runNumber);

#endif
