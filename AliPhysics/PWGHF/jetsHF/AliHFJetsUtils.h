#ifndef ALIHFJETUTILS
#define ALIHFJETUTILS

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// ******************************************
// Manager class for HF jet analysis utils
// Author: ycorrale@cern.ch
// *******************************************

#include <vector>
#include <map>

#include <Rtypes.h>

using std::vector;
using std::map;
using std::pair;
using std::make_pair;

class AliESDtrack;

typedef vector <pair<Int_t, AliESDtrack*>> vctr_pair_int_esdTrk;
typedef vector <pair <Double_t, Int_t>>   vctr_pair_dbl_int;
typedef map    <Int_t, Bool_t>            map_int_bool;

#ifndef ALILOG_H
#include "AliLog.h"
#endif

#ifndef __MSG_COLOR__
#define __MSG_COLOR__

#define COLOR_ERROR      "\033[22;31m"   //RED
#define COLOR_WARNING    "\033[22;31;1m" //ORANGE
#define COLOR_DEBUG      "\033[22;34m"   //BLUE
#define COLOR_INFO       "\033[22;32;1m" //GREEN
#define COLOR_DEFAULT    "\033[m"        //Default

#define MSGERROR(msg)   (COLOR_ERROR   msg COLOR_DEFAULT)
#define MSGWARNING(msg) (COLOR_WARNING msg COLOR_DEFAULT)
#define MSGDEBUG(msg)   (COLOR_DEBUG   msg COLOR_DEFAULT)
#define MSGINFO(msg)    (COLOR_INFO    msg COLOR_DEFAULT)

#endif //endif __MSG_COLOR__

#endif //endif ALIHFJETUTILS
