/**************************************************************************
 * Copyright(c) 2008-2010, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

/////////////////////////////////////////////////////////////////
// Author: A. Mastroserio                                      //
// This class is used as a container for the FO scan           //
//                                                             //
/////////////////////////////////////////////////////////////////


#include <TObjArray.h>
#include "AliITSOnlineSPDfoChipConfig.h"

ClassImp(AliITSOnlineSPDfoChipConfig)

AliITSOnlineSPDfoChipConfig::AliITSOnlineSPDfoChipConfig():
 TObject(),
 fChipConfigRow(0),
 fChipConfigCol(0),
 fCounter(0),
 fMatrixId(-1)
{
//constructor  
}
//__________________________

AliITSOnlineSPDfoChipConfig::AliITSOnlineSPDfoChipConfig(Short_t measure[4]):
 TObject(),
 fChipConfigRow(measure[1]),
 fChipConfigCol(measure[2]),
 fCounter(measure[3]),
 fMatrixId(measure[0])
{
//constructor
}
//
AliITSOnlineSPDfoChipConfig::AliITSOnlineSPDfoChipConfig(const AliITSOnlineSPDfoChipConfig &p):
TObject(p),
fChipConfigRow(p.fChipConfigRow),
fChipConfigCol(p.fChipConfigCol),
fCounter(p.fCounter),
fMatrixId(p.fMatrixId)
{
  //
  //copy constructor
  //
}    
