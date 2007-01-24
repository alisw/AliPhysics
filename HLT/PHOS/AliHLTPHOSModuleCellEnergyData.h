#ifndef ALIHLTPHOSNODULECELLENERGYDATA_H
#define ALIHLTPHOSNODULECELLENERGYDATA_H

struct AliHLTPHOSModuleCellEnergyData
{
  AliHLTUInt8_t fModuleID;
  unsigned long cellEnergies[64][56][2];
};


#endif
