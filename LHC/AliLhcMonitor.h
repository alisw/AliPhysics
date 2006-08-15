#ifndef ALILHCMONITOR_H
#define ALILHCMONITOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TObject.h>
class AliLhcMonitor
{
 public:
    AliLhcMonitor();
    virtual ~AliLhcMonitor(){;}
    virtual void  SetMonitor(Int_t n) = 0;
    virtual void  Record()            = 0;
    virtual void  DrawPlots()         = 0;

 protected:
    int fNmax;
    ClassDef(AliLhcMonitor,1) // LHC Monitor Base Class
};

#endif
