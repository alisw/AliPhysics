// -*- C++ -*-

#ifndef PWGOPP_VDM_PROC_PIELUP_H
#define PWGOPP_VDM_PROC_PIELUP_H

#include <vector>
#include "AliVdMMetaData.h"
#include "AliVdMScanData.h"

extern void proc_pileup(const AliVdMMetaData& vdmMetaData,
                        AliVdMScanData& allData,
                        const char* classAC,
                        const char* classAnotC,
                        const char* classCnotA,
                        const std::vector<Double_t>& par0,
                        Bool_t subtractBkgd=kTRUE);

#endif //PWGOPP_VDM_PROC_PIELUP_H
