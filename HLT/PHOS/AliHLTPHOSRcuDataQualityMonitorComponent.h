//-*- Mode: C++ -*-
// $Id$

#ifndef AliHLTPHOSRCUDATAQUALITYMONITORCOMPONENT_H
#define AliHLTPHOSRCUDATAQUALITYMONITORCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */   

#include "AliHLTPHOSProcessor.h"


class AliHLTPHOSRcuDataQualityMonitorComponent:public AliHLTPHOSProcessor
{
  AliHLTPHOSRcuDataQualityMonitorComponent();
  virtual ~AliHLTPHOSRcuDataQualityMonitorComponent();
  AliHLTPHOSRcuDataQualityMonitorComponent(const  AliHLTPHOSRcuDataQualityMonitorComponent & );
  AliHLTPHOSRcuDataQualityMonitorComponent & operator = (const AliHLTPHOSRcuDataQualityMonitorComponent &)
   {
      return *this;
   };
  //  virtual int DoInit(int argc =0, const char** argv  = 0);
  //  virtual int Deinit();
  //  virtual int DoDeinit();
 };




#endif
