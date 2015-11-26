#ifndef AliPHOSCpv3GConnection_h
#define AliPHOSCpv3GConnection_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Connection tables of pads within one 3Gassiplex card.
// 3Gassiplex card is a matrix of 8x6 pads.
// This class sets a correspondence between the absolute pad ID (0..47) and 
// (x,z) position of this pad within a 3Gassiplex card
// Author: Sergey Evdokimov, IHEP Protvino. Oct-2014

#include "Rtypes.h"

class AliPHOSCpv3GConnection
{
public:
  AliPHOSCpv3GConnection();
  Int_t Pad2X(Int_t pad)         {return pad2x[pad];  }
  Int_t Pad2Y(Int_t pad)         {return pad2y[pad];  }
  Int_t XY2Pad(Int_t x, Int_t y) {return xy2pad[x][y];}
private:
  Int_t pad2x[48];    //array of 48 elements containing 3gassiplex x coordinate for every channel
  Int_t pad2y[48];    //array of 48 elements containing 3gassiplex x coordinate for every channel
  Int_t xy2pad[8][6]; //2D array containing channels for every pad in 3gassiplex card
};

#endif //AliPHOSCpv3GConnection_h
