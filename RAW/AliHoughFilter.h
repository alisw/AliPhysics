#ifndef ALIHOUGHFILTER_H
#define ALIHOUGHFILTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// high level filter algorithm for TPC using a hough transformation
///

#include "AliFilter.h"


class AliHoughFilter: public AliFilter {
public:
  AliHoughFilter();

  virtual Bool_t       Filter(AliRawEvent* event, AliESD* esd);

private:
  ClassDef(AliHoughFilter, 0)   // TPC hough filter
};

#endif
