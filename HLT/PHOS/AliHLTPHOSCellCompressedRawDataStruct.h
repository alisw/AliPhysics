//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSCELLCOMPRESSEDRAWDATA_H
#define ALIHLTPHOSCELLCOMPRESSEDRAWDATA_H

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

struct  AliHLTPHOSCellCompressedRawData
{
  AliHLTPHOSCellCompressedRawData();
  virtual ~AliHLTPHOSCellCompressedRawData();
  
  AliHLTUInt16_t fSize;                    // Total size of data from this channel
  AliHLTUInt16_t fDataQuality;             // Quality
  AliHLTUInt8_t  fAlgorithm;               // The algorithm used to calculate energy and TOF
  AliHLTUInt8_t  fIsAttachedRawData;       // wether or not to attach the raw data
  AliHLTUInt8_t *fRawData;                 // the raw data if attached

};

#endif
