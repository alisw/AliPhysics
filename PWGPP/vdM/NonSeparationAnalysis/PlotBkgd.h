// -*- C++ -*-

#ifndef _PLOT_BKGD_H_
#define _PLOT_BKGD_H_

class TCanvas;
class TTree;
class TCut;

#include "AliXMLEngine.h"

void PlotBkgd(TCanvas *c1, Int_t fillNumber, TTree *TE, const AliXMLEngine::Node& n, const TCut& vtxCuts);

#endif // _PLOT_BKGD_H_
