#ifndef PLOTTING_H
#define PLOTTING_H

#include <TColor.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TStyle.h>

const int kNcolors = 7;
const Color_t kColor[kNcolors] = {kRed,kOrange,kGreen+1,kAzure+1,kBlue,kViolet+2,kBlack};

void SetHistStyle(TH1* h, int id = 0, int marker = 20, int linew = 1, int fillstyle = 0) {
  Color_t cc = kColor[id%kNcolors];
  h->SetMarkerColor(cc);
  h->SetLineColor(cc);
  h->SetFillStyle(fillstyle);
  h->SetLineWidth(linew);
  h->SetMarkerStyle(marker);
}

#endif
