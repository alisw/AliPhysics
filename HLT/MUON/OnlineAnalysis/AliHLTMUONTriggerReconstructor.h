#ifndef AliHLTMUONTRIGGERRECONSTRUCTOR_H
#define AliHLTMUONTRIGGERRECONSTRUCTOR_H
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

/**********************************************************************
 Created on : 16/05/2007
 Purpose    : This class is supposed to read the trigger DDL files and 
              give the output AliHLTMUONTriggerRecordStruct
 Author     : Indranil Das, HEP Division, SINP
 Email      : indra.das@saha.ac.in | indra.ehep@gmail.com
**********************************************************************/

#include <vector>
#include <TArrayS.h>

#include <AliHLTLogging.h>
#include "AliHLTMUONTriggerRecordsBlockStruct.h"
#include "AliHLTMUONHitReconstructor.h"
#include "AliMUONTriggerCrateStore.h"

class AliHLTMUONTriggerReconstructor : public AliHLTLogging
{

 public:

  AliHLTMUONTriggerReconstructor();
  virtual ~AliHLTMUONTriggerReconstructor();

  bool LoadLookUpTable(AliHLTMUONHitReconstructor::DHLTLut* lookUpTableData, int lookUpTableId);

  //bool Run(int iEvent, int iDDL, AliHLTMUONTriggerRecordStruct trigRecord, int *nofTrigRec); // for Reading using rawreader
  bool Run(int *rawData, int *rawDataSize, AliHLTMUONTriggerRecordStruct trigRecord[], int *nofTrigRec);

  int GetLutLine(){return fgkLutLine ;}

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

  AliHLTMUONTriggerReconstructor(const AliHLTMUONTriggerReconstructor& rhs); // copy constructor
  AliHLTMUONTriggerReconstructor& operator=(const AliHLTMUONTriggerReconstructor& rhs); // assignment operator
  
 private:

  static const int fgkEvenLutSize ;           // Size of the LookupTable with event DDLID
  static const int fgkOddLutSize ;            // Size of the LookupTable with odd DDLID
  static const int fgkLutLine;                // nof Line in LookupTable    

  static const int fgkMinIdManuChannel[2];    // Minimum value of idManuChannel in LookupTable  
  static const int fgkMaxIdManuChannel[2];    // Maximum value of idManuChannel in LookupTable  
  static const float fgkHalfPadSizeXB[3];       // pad halflength for the pcb zones  
  static const float fgkHalfPadSizeYNB[2];       // pad halflength for the pcb zones  

  static const int fgkDetElem;                // nof Detection element per DDL    

  
  AliHLTMUONHitReconstructor::DHLTPad* fPadData;                          // pointer to the array containing the information of each padhits
  AliHLTMUONHitReconstructor::DHLTLut* fLookUpTableData;                  // pointer to the array of Lookuptable data
  
  AliHLTMUONTriggerRecordStruct *fRecPoints;    // Reconstructed hits
  int *fRecPointsCount;                       // nof reconstructed hit  
  int fMaxRecPointsCount;                    // max nof reconstructed hit  
  int fDigitPerDDL;                                    // Total nof Digits perDDL 


  int fGetIdTotalData[104][64][2] ;           // an array of idManuChannel with argumrnt of centralX,centralY and  planeType
  int fNofFiredDetElem,*fMaxFiredPerDetElem;        // counter for detector elements that are fired 
  int *fDetManuChannelIdList;                          // pointer to an array of idManuChannel
  int *fCentralChargeB,*fCentralChargeNB;              // pointer to an array of central hit


  int fDDLId ;
  int fIdOffSet ;

  AliMUONTriggerCrateStore* fCrateManager;

  bool MergeTrigHits(int minPadId, int maxPadId);
  bool FindTrigHits() ;

  bool ReadDDL(int *rawData, int *rawDataSize);
  bool Pattern2Pad(int nBoard, TArrayS* xyPattern, vector<AliHLTMUONHitReconstructor::DHLTPad>& padList);

};

#endif // AliHLTMUONTRIGGERRECONSTRUCTOR_H
