/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Author:  Per Thomas Hille  <perthi@fys.uio.no>                 *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include  "AliHLTPHOSDDLPackedFileWriter.h"


//_________________________________________________________________________________________________
AliHLTPHOSDDLPackedFileWriter::AliHLTPHOSDDLPackedFileWriter()
{

}


//_________________________________________________________________________________________________
AliHLTPHOSDDLPackedFileWriter::~AliHLTPHOSDDLPackedFileWriter()
{

}


//_________________________________________________________________________________________________
const int 
AliHLTPHOSDDLPackedFileWriter::WriteFile(const AliHLTComponentEventData& /*evtData*/, 
			const AliHLTComponentBlockData* /*blocks*/, AliHLTComponentTriggerData& /*trigData*/, int /*evntCnt*/) const
{
  return 0;
}
