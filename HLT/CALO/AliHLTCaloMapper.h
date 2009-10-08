//-*- Mode: C++ -*-
// $Id: AliHLTCALOMapper.h 34264 2009-08-14 18:29:23Z odjuvsla $

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


#include "AliHLTCaloConstants.h"
#include "Rtypes.h"
#include "AliHLTLogging.h"

using namespace CaloHLTConst;


class AliHLTCaloCoordinate;

class AliHLTCaloMapper : public AliHLTLogging
{
public:

  //  AliHLTCaloMapper(const int speficication);

  AliHLTCaloMapper(const unsigned long speficication);
 
  virtual ~AliHLTCaloMapper();

  // virtual void InitAltroMapping( const int specification ) = 0; 
  virtual void InitAltroMapping( const unsigned long specification ) = 0; 

  virtual void InitDDLSpecificationMapping() = 0;
  bool GetIsInitializedMapping();
  char* GetFilePath();
  //const virtual int GetDDLFromSpec( const AliHLTUInt32_t spec );
  // const virtual int GetChannelID(const AliHLTUInt32_t spec, const Int_t hadd);

  //  static int GetDDLFromSpec( const AliHLTUInt32_t spec );
  int GetDDLFromSpec( const AliHLTUInt32_t spec );

  
  static unsigned long GetSpecFromDDLIndex( const int ddlindex );

 


  UInt_t  GetChannelID(const AliHLTUInt32_t spec, const Int_t hadd);
  static void ChannelId2Coordinate(const UShort_t channelId,    AliHLTCaloCoordinate &channelCoord);
  
  // static unsigned long GetSpecFromDDLIndex( const int ddlindex );
  // static int GetChannelID(const AliHLTUInt32_t spec, const Int_t hadd);
  // static void ChannelId2Coordinate(const UShort_t channelId,    AliHLTCaloCoordinate &channelCoord);
  

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
  int fSpecification;
  bool fIsInitializedMapping;
  fDDLSpecificationMap* fSpecificationMapPtr;
  
 private:
  AliHLTCaloMapper();
  AliHLTCaloMapper(const AliHLTCaloMapper & );
  AliHLTCaloMapper & operator = (const AliHLTCaloMapper &);
  
};

#endif
