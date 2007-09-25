#ifndef ALIHLTMUONHITRECONSTRUCTOR_H
#define ALIHLTMUONHITRECONSTRUCTOR_H
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Indranil Das <indra.das@saha.ac.in>                                  *
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

///////////////////////////////////////////////
//Author : Indranil Das, SINP, INDIA
//         Sukalyan Chattopadhyay, SINP, INDIA
//         
//Email :  indra.das@saha.ac.in
//         sukalyan.chattopadhyay@saha.ac.in 
///////////////////////////////////////////////

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>

#include <TString.h>

#include <AliHLTLogging.h>

#include "AliRawReader.h"
#include "AliRawReaderFile.h"
#include "AliRawReaderRoot.h"
#include "AliRawReaderDate.h"

#if __GNUC__ < 3
#define std
#endif


typedef std::map<int,int> BusToDetElem;
typedef std::map<int,int> BusToDDL;

extern "C" struct AliHLTMUONRecHitStruct;


class AliHLTMUONHitReconstructor : public AliHLTLogging
{
public:

  struct DHLTLut{
    int fIdManuChannel;  // Manu channel address
    int fIX,fIY;  // The X,Y number of the pad.
    float fRealX,fRealY,fRealZ;  // The real coordinate of the pad.
    int fPlane,fPcbZone;  // The plane and PCB zone ID numbers.
  };

  struct DHLTPad{
    int fDetElemId;  // The detector element ID of the pad.
    int fBuspatchId;  // The bus patch ID number for the pad.
    int fIdManuChannel;  // The Manu channel address.
    int fIX,fIY;  // The X,Y number of the pad.
    float fRealX,fRealY,fRealZ;   // The real coordinate of the pad.
    int fPlane,fPcbZone;   // The plane and PCB zone ID numbers.
    int fCharge;  // The charge measured on the pad.
  };

  AliHLTMUONHitReconstructor();
  virtual ~AliHLTMUONHitReconstructor(void);

  bool LoadLookUpTable(DHLTLut* lookUpTableData, int lookUpTableId);
  bool SetBusToDetMap(BusToDetElem busToDetElem);
  bool SetBusToDDLMap(BusToDDL busToDDL);
  
  bool Run(
  		const AliHLTUInt32_t* rawData,
  		AliHLTUInt32_t rawDataSize,
  		AliHLTMUONRecHitStruct* recHit,
  		AliHLTUInt32_t& nofHit
  	);
  void SetDCCut(int dcCut) {fDCCut = dcCut;}
  void SetDebugLevel(int debugLevel) {fDebugLevel = debugLevel;}
  int GetDebugLevel() const {return fDebugLevel;}
  int GetDEId(int iBusPatch) {return fBusToDetElem[iBusPatch] ;}           
    
  int GetLutLine(int iDDL) const;

  static int GetkDetectorId() { return fgkDetectorId; }
  static int GetkDDLOffSet() { return fgkDDLOffSet; }
  static int GetkNofDDL() { return fgkNofDDL; }
  static int GetkDDLHeaderSize() { return fgkDDLHeaderSize; }
  
private: 
  static const int fgkDetectorId ;            // DDL Offset
  static const int fgkDDLOffSet ;             // DDL Offset
  static const int fgkNofDDL ;                // Number of DDL 
  static const int fgkDDLHeaderSize  ;        // DDL header size
 
protected:
  AliHLTMUONHitReconstructor(const AliHLTMUONHitReconstructor& rhs); // copy constructor
  AliHLTMUONHitReconstructor& operator=(const AliHLTMUONHitReconstructor& rhs); // assignment operator
  
private:

  static const int fgkEvenLutSize ;          // Size of the LookupTable with event DDLID
  static const int fgkOddLutSize ;           // Size of the LookupTable with odd DDLID
  static const int fgkLutLine[2];            // nof Line in LookupTable    

  static const int fgkMinIdManuChannel[2];   // Minimum value of idManuChannel in LookupTable  
  static const int fgkMaxIdManuChannel[2];   // Maximum value of idManuChannel in LookupTable  
  static const float fgkHalfPadSize[3];      // pad halflength for the pcb zones  
  
  int fkBlockHeaderSize;                     // Block header size
  int fkDspHeaderSize;                       // DSP header size
  int fkBuspatchHeaderSize;                  // buspatch header size
  
  int fDCCut;                                // DC Cut value

  DHLTPad* fPadData;                         // pointer to the array containing the information of each padhits
  DHLTLut* fLookUpTableData;                 // pointer to the array of Lookuptable data
  
  AliHLTMUONRecHitStruct* fRecPoints;       // Reconstructed hits
  AliHLTUInt32_t *fRecPointsCount;                      // nof reconstructed hit  
  AliHLTUInt32_t fMaxRecPointsCount;                    // max nof reconstructed hit  
   
  int fCentralCountB, fCentralCountNB;        // centeral hits 
  int fIdOffSet,fDDLId;                       // DDLId and DDL id offset
  int fDigitPerDDL;                           // Total nof Digits perDDL 
  
  int *fDetManuChannelIdList;                          // pointer to an array of idManuChannel
  int *fCentralChargeB,*fCentralChargeNB;              // pointer to an array of central hit
  float *fRecX,*fRecY;                                 // pointer to an array of reconstructed hit
  float *fAvgChargeX,*fAvgChargeY;                     // average charge on central pad found using CG method
  int fGetIdTotalData[336][80][2] ;                    // an array of idManuChannel with argumrnt of centralX,centralY and  planeType
  int fNofFiredDetElem,fMaxFiredPerDetElem[13];        // counter for detector elements that are fired 
  int fDebugLevel;
  BusToDetElem fBusToDetElem;             // Mapping between bus address and detector element ID.
  BusToDetElem fBusToDDL;                 // Mapping between bus address and DDL.


  bool ReadDDL(const AliHLTUInt32_t* rawData, AliHLTUInt32_t rawDataSize);
  void FindCentralHits(int minPadId, int maxPadId);
  bool FindRecHits() ;
  void RecXRecY();
  bool MergeRecHits();
  void Clear();
};

#endif // ALIHLTMUONHITRECONSTRUCTOR_H
