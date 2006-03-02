#ifndef MUONCDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// By Laurent Aphecetche

#include "Rtypes.h"

class TList;
class AliMUONV2DStore;

void generateCalibrations(const char* cdbpath, Bool_t defaultValues = kTRUE);

TList* manuList(Bool_t reset=kFALSE);

void plotCDB(const char* calibType="MUON/Calib/Pedestals");

AliMUONV2DStore* readCDB(const char* calibType="MUON/Calib/Pedestals");

void testMakeStores(Int_t readLoop=10);

#endif
