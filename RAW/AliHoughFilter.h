#ifndef ALIHOUGHFILTER_H
#define ALIHOUGHFILTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// high level filter algorithm for TPC using a hough transformation
///

#include "AliFilter.h"

class TTree;

class AliITSgeom;


class AliHoughFilter: public AliFilter {
public:
  AliHoughFilter();

  virtual Bool_t       Filter(AliRawEvent* event, AliESDEvent* esd);

  void                 RunITSclusterer(AliRawEvent* event, TTree *treeClusters);
  void                 RunITSvertexer(AliESDEvent* esd, TTree *treeClusters);
  void                 RunTPCtracking(AliRawEvent* event, AliESDEvent* esd);
  void                 RunITStracking(AliESDEvent* esd, TTree *treeClusters);
private:
  AliHoughFilter(const AliHoughFilter&);
  AliHoughFilter &operator=(const AliHoughFilter&);

  Float_t fPtmin;        //Low limit on Pt

  AliITSgeom *fITSgeom;  //Pointer to the ITS geometry

  ClassDef(AliHoughFilter, 0)   // TPC hough filter
};

#endif
