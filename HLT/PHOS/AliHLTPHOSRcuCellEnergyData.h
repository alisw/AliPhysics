#ifndef ALIHLTPHOSRCUCELLENERGYDATA_H
#define ALIHLTPHOSRCUCELLENERGYDATA_H

struct AliHLTPHOSRcuCellEnergyData
{
  AliHLTUInt16_t fRcuX;
  AliHLTUInt16_t fRcuZ;
  AliHLTUInt16_t fModuleID;
  
  Double_t fCellEnergies[32][28][2];
  
};


#endif
