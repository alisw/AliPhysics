#ifndef ALIHLTPHOSDQM_H
#define ALIHLTPHOSDQM_H

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
