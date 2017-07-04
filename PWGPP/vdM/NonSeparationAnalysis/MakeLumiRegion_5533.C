// -*- C++ -*-

#include <sstream>

#include "AliLuminousRegionFit.h"
#include "AliVdMData.h"

void MakeLumiRegion()
{
  AliVdMData d(AliVdMData::GetFileName("5533/5533.xml"));

  const TString vtxFileName = "5533/5533_vtx.root";
  const TCut    vtxCuts     = "ntrksTRKnc>=11 && chi2/ntrksTRKnc<2";
  const Int_t   bcidSel     = -1; // all BCIDs

  std::for_each(d.GetScansBegin(), d.GetScansEnd(), [&](const AliXMLEngine::Node& n) {
      std::istringstream iss(n.GetData());
      AliLuminousRegionFit f(d.GetFillNumber(),
                             AliVdMData::GetFileName(vtxFileName),
                             iss);
      f.DoFit(AliVdMData::GetScanName(n),
              AliVdMData::GetScanType(n),
              std::stod(n.GetAttr("offset_um").GetData()),
              vtxCuts,
              bcidSel);
    });

  Printf("#scans: %ld", std::distance(d.GetScansBegin(), d.GetScansEnd()));
}
