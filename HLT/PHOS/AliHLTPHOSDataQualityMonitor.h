//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSDQM_H
#define ALIHLTPHOSDQM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


class AliHLTPHOSDataQualityMonitor
{
  AliHLTPHOSDataQualityMonitor();
  virtual ~AliHLTPHOSDataQualityMonitor();
  AliHLTPHOSDataQualityMonitor(const AliHLTPHOSDataQualityMonitor & );
  AliHLTPHOSDataQualityMonitor & operator = (const AliHLTPHOSDataQualityMonitor &)
   {
      return *this;
   };

};

#endif
