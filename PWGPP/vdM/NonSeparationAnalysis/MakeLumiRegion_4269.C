// -*- C++ -*-
// $Id$

#include <sstream>

#include "AliLuminousRegionFit.h"
#include "AliVdMData.h"

void MakeLumiRegion()
{
  AliVdMData d(AliVdMData::GetFileName("4937/4937.xml"));

  const TString vtxFileName = "4269/4269_vtx.root";
  const TCut    vtxCuts     = "ntrksTRKnc>=11 && chi2/ntrksTRKnc<2";
  const TCut    bkgdCuts    = "timeAD.A>55.7 && timeAD.A<57.8 && timeAD.C>64.3 && timeAD.C<66 && timeV0.C>1.8 && timeV0.C<3.5 && timeV0.A>9.8 && timeV0.A<11.4";
  const Int_t   bcidSel     = -1; // all BCIDs

  std::for_each(d.GetScansBegin(), d.GetScansEnd(), [&](const AliXMLEngine::Node& n) {
      std::istringstream iss(n.GetData());
      AliLuminousRegionFit f(d.GetFillNumber(),
                             AliVdMData::GetFileName(vtxFileName),
                             iss);
      f.DoFit(AliVdMData::GetScanName(n),
              AliVdMData::GetScanType(n),
              std::stoi(n.GetAttr("offset_um").GetData()),
              vtxCuts*bkgdCuts,
              bcidSel);
    });

  Printf("#scans: %ld", std::distance(d.GetScansBegin(), d.GetScansEnd()));
}
