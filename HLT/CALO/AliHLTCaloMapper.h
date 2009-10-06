//-*- Mode: C++ -*-
// $Id: AliHLTCaloMapper.h 34264 2009-08-14 18:29:23Z odjuvsla $

#ifndef ALIHLTCALOMAPPER_H
#define ALIHLTCALOMAPPER_H


/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2006                                       *
 *                                                                        * 
 * Author: Per Thomas Hille perthi@fys.uio.no for the ALICE DCS Project.  *
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

//#include "AliHLTCaloBase.h"

#include "AliHLTCaloConstants.h"

using namespace CaloHLTConst;

#include "Rtypes.h"
#include "AliHLTLogging.h"

//#include "AliPHOSEMCAGeometry.h"

class AliHLTCaloCoordinate;

class AliHLTCaloMapper : public AliHLTLogging
//class AliHLTCaloMapper 
{
 public:
  AliHLTCaloMapper();
  virtual ~AliHLTCaloMapper();
  
  virtual void InitAltroMapping() = 0; 
  virtual void InitDDLSpecificationMapping() = 0;
  
  bool GetIsInitializedMapping();
  char* GetFilePath();
  const virtual int GetDDLFromSpec( const AliHLTUInt32_t spec )=0;
  const virtual int GetChannelID(const AliHLTUInt32_t spec, const Int_t hadd) = 0;
  void ChannelId2Coordinate(const UShort_t channelId,    AliHLTCaloCoordinate &channelCoord);
  //  UShort_t GetChannelID(Int_t specification, Int_t hwAddress);
  //  static void GetChannelCoord(UShort_t channelId, UShort_t* channelCoord);
  //  static void GetLocalCoord(UShort_t channelId, Float_t* localCoord);
  //  int GetDDLFromSpec(Int_t specification);

  struct fAltromap{ 
    int fZRow; // Coordinate in Z direction (beam direction) relatve too one RCU
    int fXCol; // Coordinate in X direction (perpendicular too beam direction an parallell to ground) relatve too one RCU
    int fGain; // Gain (high gain = 1, low gain = 0)
  };
  
  struct fDDLSpecificationMap{ 
    UInt_t fRcuX; // Coordinate in Z direction (beam direction) relatve too one RCU
    UInt_t fRcuZ; // Coordinate in X direction (perpendicular too beam direction an parallell to ground) relatve too one RCU
    UInt_t fRcuXOffset;
    UInt_t fRcuZOffset;
    int fModId; 
  };

  fAltromap *fHw2geomapPtr; //pointer to structure holding information about geometrical address 


  char fFilepath[1024];

protected:
  bool fIsInitializedMapping;
  fDDLSpecificationMap* fSpecificationMapPtr;


 private:
  AliHLTCaloMapper(const AliHLTCaloMapper & );
  AliHLTCaloMapper & operator = (const AliHLTCaloMapper &);
  //AliPHOSEMCAGeometry fPHOSGeometry;
  //static const Float_t fCellStep = 2.2;//2.*fPHOSGeometry.GetAirCellHalfSize()[0];
};

#endif
