// -*- C++ -*-

#ifndef PWGOPP_VDM_PROC_BKGD_H
#define PWGOPP_VDM_PROC_BKGD_H

#include "AliVdMMetaData.h"
#include "AliVdMScanData.h"

extern void proc_bkgd(const AliVdMMetaData& vdmMetaData,
                      AliVdMScanData& allData,
                      const char* className);

#endif //PWGOPP_VDM_PROC_BKGD_H
