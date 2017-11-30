// -*- C++ -*-
// $Id$

#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TCut.h>

#include "AliVdMData.h"

#include "PlotBkgd.h"

void PlotBkgdLumiRegion()
{
  AliVdMData d(AliVdMData::GetFileName("5568/5568.xml"));

  const TString vtxFileName = "5568/5568_vtx.root";
  const TCut    vtxCuts     = "ntrksTRKnc>=11 && chi2/ntrksTRKnc<2";
  const TCut    bkgdCuts    = "timeAD.A>55.7 && timeAD.A<57.8 && timeAD.C>64.3 && timeAD.C<66 && timeV0.C>1.8 && timeV0.C<3.5 && timeV0.A>9.8 && timeV0.A<11.4";
  const Int_t   bcidSel     = -1; // all BCIDs
  const TString pn          = "pdf/5568/LumiRegionBkgd_5568.pdf";

  TFile *f = TFile::Open(AliVdMData::GetFileName(vtxFileName));
  TTree *TE = (TTree*)gDirectory->Get("Vertex_Performance/cOutputVtxESD");

  TCanvas *c1 = new TCanvas;
  c1->SaveAs(pn+"[");
  std::for_each(d.GetScansBegin(), d.GetScansEnd(), [&](const AliXMLEngine::Node& n) {
      PlotBkgd(c1, 5568, TE, n, vtxCuts);
      c1->SaveAs(pn);
      PlotBkgd(c1, 5568, TE, n, vtxCuts*bkgdCuts);
      c1->SaveAs(pn);
    });
  c1->SaveAs(pn+"]");

  delete c1;
  delete TE;
  f->Close();

  Printf("#scans: %ld", std::distance(d.GetScansBegin(), d.GetScansEnd()));
}
