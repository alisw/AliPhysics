// -*- mode: C++ -*-
//
#ifndef ALIFMDCALIBFWD_H
#define ALIFMDCALIBFWD_H
#ifndef ALIFMDUSHORTMAP_H
# include <AliFMDUShortMap.h>
#endif
#ifndef ALIFMDBOOLMAP_H
# include <AliFMDBoolMap.h>
#endif
typedef AliFMDUShortMap AliFMDCalibZeroSuppression;
typedef AliFMDBoolMap   AliFMDCalibDeadMap;
class AliFMDCalibPedestal;
class AliFMDCalibGain;
class AliFMDCalibSampleRate;
class AliFMDCalibStripRange;
class AliFMDAltroMapping;

#endif
//
// EOF
//
