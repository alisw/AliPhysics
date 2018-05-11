// -*- C++ -*-

#ifndef PWGPP_VDM_PROC_H
#define PWGPP_VDM_PROC_H

#include <vector>
#include <string>

#include "AliVdMMetaData.h"
#include "AliVdMScanData.h"

void proc(const AliVdMMetaData&, AliVdMScanData&, const std::vector<std::string>& triggerNames, Bool_t computeBkgd=kTRUE);

#endif // PWGPP_VDM_PROC_H
