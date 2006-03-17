#ifndef MUONCDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// By Laurent Aphecetche

#include "Rtypes.h"

class TList;
class AliMUONV1DStore;
class AliMUONV2DStore;

static const char* CDBPath = "local://$ALICE_ROOT/";

void generateTrigger(const char* cdbpath=CDBPath);

void generateCalibrations(const char* cdbpath=CDBPath, Bool_t defaultValues = kTRUE);

TList* manuList(Bool_t reset=kFALSE);

void plotCDB(const char* calibType="MUON/Calib/Pedestals");

AliMUONV2DStore* read2D(const char* calibType="MUON/Calib/Pedestals");
AliMUONV1DStore* read1D(const char* calibType="MUON/Calib/LocalBoardMasks");

void testMakeStores(Int_t readLoop=10);

#endif
