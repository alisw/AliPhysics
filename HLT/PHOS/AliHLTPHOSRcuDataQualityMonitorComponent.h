#ifndef AliHLTPHOSRCUDATAQUALITYMONITORCOMPONENT_H
#define AliHLTPHOSRCUDATAQUALITYMONITORCOMPONENT_H

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
