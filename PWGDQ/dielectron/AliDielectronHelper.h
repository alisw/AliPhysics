#ifndef ALIDIELECTRONHELPER_H
#define ALIDIELECTRONHELPER_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
// Dielectron helpers                                                                    //
//                                                                                       //
//                                                                                       //
// Authors:                                                                              //
//   Jens Wiechula <Jens.Wiechula@cern.ch>                                               //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////


#include <TVectorDfwd.h>

class AliKFParticle;
class AliVEvent;
class AliMCEvent;

namespace AliDielectronHelper
{



TVectorD* MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);
TVectorD* MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);
TVectorD* MakeArbitraryBinning(const char* bins);

Int_t GetNch(const AliMCEvent *ev=0x0, Double_t eta=0.9);
Int_t GetNacc(const AliVEvent *ev=0x0);
Int_t GetNaccTrcklts(const AliVEvent *ev=0x0);

void RotateKFParticle(AliKFParticle * kfParticle,Double_t angle, const AliVEvent * const ev=0x0);


}

#endif
