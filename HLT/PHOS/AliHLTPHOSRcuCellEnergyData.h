#ifndef ALIHLTPHOSRCUCELLENERGYDATA_H
#define ALIHLTPHOSRCUCELLENERGYDATA_H

struct AliHLTPHOSRcuCellEnergyData
{
  AliHLTUInt8_t fRcuX;
  AliHLTUInt8_t fRcuY;
  AliHLTUInt8_t fModuleID;
  
  unsigned long cellEnergies[32][28][2];
  
};


#endif
