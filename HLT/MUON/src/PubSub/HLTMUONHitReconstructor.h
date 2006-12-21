///////////////////////////////////////////////
//Author : Indranil Das, SINP, INDIA
//         Sukalyan Chattopadhyay, SINP, INDIA
//         
//
//Email :  indra.das@saha.ac.in
//         sukalyan.chattopadhyay@saha.ac.in 
///////////////////////////////////////////////

#ifndef HLTMUONHITRECONSTRUCTOR_H
#define HLTMUONHITRECONSTRUCTOR_H

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>

// and hand control to the PubSub logging level
//#define DEBUG


using namespace std; 

struct DHLTLut{
  int fIdManuChannel;
  int fIX,fIY;
  float fRealX,fRealY,fRealZ;
  int fPlane,fPcbZone;
};

struct DHLTPad{
  int fDetElemId;
  int fBuspatchId;
  int fIdManuChannel;
  int fIX,fIY;
  float fRealX,fRealY,fRealZ;
  int fPlane,fPcbZone;
  int fCharge;
};

struct DHLTRecPoint{
  float X,Y,Z;
  int DetElemId;
};

typedef  map<int,int> BusToDetElem ;

class HLTMUONHitReconstructor 
{
public:
  HLTMUONHitReconstructor();
  virtual ~HLTMUONHitReconstructor(void);

  bool LoadLookUpTable(DHLTLut* lookUpTableData, int lookUpTableId);
  bool SetBusToDetMap(BusToDetElem busToDetElem);
  
  //bool Init();
  bool Run(const int* rawData, int rawDataSize, DHLTRecPoint recHit[], int *nofHit);

  void SetDCCut(int dcCut) {fDCCut = dcCut;}
  void SetDebugLevel(int debugLevel) {fDebugLevel = debugLevel;}
  int GetDebugLevel(){return fDebugLevel;}
  int GetDEId(int iBusPatch) {return fBusToDetElem[iBusPatch] ;}           
    
  const int GetLutLine(int iDDL);
  
  static const int fgkDetectorId ;            // DDL Offset
  static const int fgkDDLOffSet ;             // DDL Offset
  static const int fgkNofDDL ;                // Number of DDL 
  static const int fgkDDLHeaderSize  ;                      // DDL header size  
protected:
  HLTMUONHitReconstructor(const HLTMUONHitReconstructor& rhs); // copy constructor
  HLTMUONHitReconstructor& operator=(const HLTMUONHitReconstructor& rhs); // assignment operator
private:
  

  static const int fgkEvenLutSize ;           // Size of the LookupTable with event DDLID
  static const int fgkOddLutSize ;            // Size of the LookupTable with odd DDLID
  static const int fgkLutLine[2];             // nof Line in LookupTable    

  static const int fgkMinIdManuChannel[2];    // Minimum value of idManuChannel in LookupTable  
  static const int fgkMaxIdManuChannel[2];    // Maximum value of idManuChannel in LookupTable  
  static const float fgkHalfPadSize[3];       // pad halflength for the pcb zones  
  

  int fkBlockHeaderSize ;                     // Block header size
  int fkDspHeaderSize  ;                      // DSP header size
  int fkBuspatchHeaderSize ;                  // buspatch header size
  
  int fDCCut;                                 // DC Cut value

  DHLTPad* fPadData;                          // pointer to the array containing the information of each padhits
  DHLTLut* fLookUpTableData;                      // pointer to the array of Lookuptable data
  
  DHLTRecPoint *fRecPoints;                   // Reconstructed hits
  int *fRecPointsCount;                       // nof reconstructed hit  
  int fMaxRecPointsCount;                    // max nof reconstructed hit  
   
  int fCentralCountB, fCentralCountNB;                 // centeral hits 
  int fIdOffSet,fDDLId;                       // DDLId and DDL id offset
  int fDigitPerDDL;                                    // Total nof Digits perDDL 
  
  int *fDetManuChannelIdList;                          // pointer to an array of idManuChannel
  int *fCentralChargeB,*fCentralChargeNB;              // pointer to an array of central hit
  float *fRecX,*fRecY;                                 // pointer to an array of reconstructed hit
  float *fAvgChargeX,*fAvgChargeY;                     // average charge on central pad found using CG method
  int fGetIdTotalData[336][80][2] ;                    // an array of idManuChannel with argumrnt of centralX,centralY and  planeType
  int fNofFiredDetElem,fMaxFiredPerDetElem[13];        // counter for detector elements that are fired 
  int fDebugLevel;
  BusToDetElem fBusToDetElem;

  bool ReadDDL(const int* rawData, int rawDataSize);
  void FindCentralHits(int minPadId, int maxPadId);
  bool FindRecHits() ;
  void RecXRecY();
  bool MergeRecHits();
  void CleanUp();

};


#endif // HLTMUONHITRECONSTRUCTOR_H
