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

////////////////////////////////////////////////////////////////////////////////
//
//  This class provides storage for event and track information which 
//  are used for same-event as well as mixed-event analyses in AliXiStar 
//
//  authors: Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////


#include "AliXiStarPbPbEventCollection.h"

AliXiStarPbPbTrackStruct::AliXiStarPbPbTrackStruct():
  fStatus(0),
  fFilterMap(0),
  fID(0),
  fPhi(0),
  fPt(0),
  fMom(0),
  fP(),
  fCharge(0),
  fEta(0),
  fMass(0),
  fDCAXY(0),
  fDCAZ(0),
  fDCA(0),
  fX(),
  fCov(),
  fNSigmaPi(0),
  fNSigmaK(0),
  fNSigmaPr(0),
  fLabel(0),
  fNclusTPC(0)
{
  //Default constructor
}
//_____________________________________________________________________________
AliXiStarPbPbTrackStruct::AliXiStarPbPbTrackStruct(const AliXiStarPbPbTrackStruct &obj) 
  : fStatus(obj.fStatus),
    fFilterMap(obj.fFilterMap),
    fID(obj.fID),
    fPhi(obj.fPhi),
    fPt(obj.fPt),
    fMom(obj.fMom),
    fP(),
    fCharge(obj.fCharge),
    fEta(obj.fEta),
    fMass(obj.fMass),
    fDCAXY(obj.fDCAXY),
    fDCAZ(obj.fDCAZ),
    fDCA(obj.fDCA),
    fX(),
    fCov(),
    fNSigmaPi(obj.fNSigmaPi),
    fNSigmaK(obj.fNSigmaK),
    fNSigmaPr(obj.fNSigmaPr),
    fLabel(obj.fLabel),
    fNclusTPC(obj.fNclusTPC)
{
  // copy constructor
}

//_____________________________________________________________________________
AliXiStarPbPbTrackStruct &AliXiStarPbPbTrackStruct::operator=(const AliXiStarPbPbTrackStruct &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;
  
  fStatus = obj.fStatus;
  fFilterMap = obj.fFilterMap;
  fID = obj.fID;
  fPhi = obj.fPhi;
  fPt = obj.fPt;
  fMom = obj.fMom;
  fP[0] = obj.fP[0];
  fP[1] = obj.fP[1];
  fP[2] = obj.fP[2];
  fCharge = obj.fCharge;
  fEta = obj.fEta;
  fMass = obj.fMass;
  fDCAXY = obj.fDCAXY;
  fDCAZ = obj.fDCAZ;
  fDCA = obj.fDCA;
  fX[0] = obj.fX[0];
  fX[1] = obj.fX[1];
  fX[2] = obj.fX[2];
  fCov[0] = obj.fCov[0]; fCov[1] = obj.fCov[1]; fCov[2] = obj.fCov[2];
  fCov[3] = obj.fCov[3]; fCov[4] = obj.fCov[4]; fCov[5] = obj.fCov[5];
  fCov[6] = obj.fCov[6]; fCov[7] = obj.fCov[7]; fCov[8] = obj.fCov[8];
  fCov[9] = obj.fCov[9]; fCov[10] = obj.fCov[10]; fCov[11] = obj.fCov[11];
  fCov[12] = obj.fCov[12]; fCov[13] = obj.fCov[13]; fCov[14] = obj.fCov[14];
  fCov[15] = obj.fCov[15]; fCov[16] = obj.fCov[16]; fCov[17] = obj.fCov[17];
  fCov[18] = obj.fCov[18]; fCov[19] = obj.fCov[19]; fCov[20] = obj.fCov[20];
  fNSigmaPi = obj.fNSigmaPi;
  fNSigmaK = obj.fNSigmaK;
  fNSigmaPr = obj.fNSigmaPr;
  fLabel = obj.fLabel;
  fNclusTPC = obj.fNclusTPC;

  return (*this);
}

//_____________________________________________________________________________
AliXiStarPbPbTrackStruct::~AliXiStarPbPbTrackStruct()
{
  // Destructor
}
//_____________________________________________________________________________
AliXiStarPbPbEventStruct::AliXiStarPbPbEventStruct():
  fNTracks(0),
  fTracks(0)
{
  //Default constructor
}
//_____________________________________________________________________________
AliXiStarPbPbEventStruct::AliXiStarPbPbEventStruct(const AliXiStarPbPbEventStruct &obj) 
  :  fNTracks(obj.fNTracks),
     fTracks(obj.fTracks)
{
  // copy constructor
}
//_____________________________________________________________________________
AliXiStarPbPbEventStruct &AliXiStarPbPbEventStruct::operator=(const AliXiStarPbPbEventStruct &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;
  
  fNTracks = obj.fNTracks;
  fTracks = obj.fTracks;
  
  return (*this);
}
//_____________________________________________________________________________
AliXiStarPbPbEventStruct::~AliXiStarPbPbEventStruct()
{
  // Destructor
  if(fTracks) delete fTracks;
}
//_____________________________________________________________________________
AliXiStarPbPbEventCollection::AliXiStarPbPbEventCollection():
  fFIFO(0),
  fEvtStr(0)
{
//Default constructor
}
//_____________________________________________________________________________
AliXiStarPbPbEventCollection::AliXiStarPbPbEventCollection(short a):
  fFIFO(0),
  fEvtStr(0x0)
{
  // main constructor
  SetBuffSize(a);
  
  fEvtStr = new AliXiStarPbPbEventStruct[fFIFO];  //allocate pointer array
  for(Int_t ii = 0; ii < fFIFO; ii++){   //Initialize to NULL
    (fEvtStr + ii)->fTracks = NULL;
    (fEvtStr + ii)->fNTracks = 0;
    (fEvtStr + ii)->fTracks = new AliXiStarPbPbTrackStruct[20000];
  }
  
}
//_____________________________________________________________________________
AliXiStarPbPbEventCollection::AliXiStarPbPbEventCollection(const AliXiStarPbPbEventCollection &obj) 
  :  fFIFO(obj.fFIFO),
     fEvtStr(obj.fEvtStr)
{
  // copy constructor
}
//_____________________________________________________________________________
AliXiStarPbPbEventCollection &AliXiStarPbPbEventCollection::operator=(const AliXiStarPbPbEventCollection &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;
  
  fFIFO = obj.fFIFO;
  fEvtStr = obj.fEvtStr;
  
  return (*this);
}

//_____________________________________________________________________________
AliXiStarPbPbEventCollection::~AliXiStarPbPbEventCollection(){

    for(Int_t i = 0; i < fFIFO; i++){

	if((fEvtStr + i)->fTracks != NULL){
	  delete [] (fEvtStr + i)->fTracks;
	}	
	
    }
    
    delete [] fEvtStr;
    //remove histos from heap

}
//_____________________________________________________________________________
void AliXiStarPbPbEventCollection::FIFOShift(){ //Shift elements in FIFO by one and clear last element in FIFO 
  
  for(UShort_t i=fFIFO-1 ; i > 0; i--){
    for(Int_t j=0; j<(fEvtStr + i-1)->fNTracks; j++) (fEvtStr + i)->fTracks[j] = (fEvtStr + i-1)->fTracks[j];
    (fEvtStr + i)->fNTracks = (fEvtStr + i-1)->fNTracks;
  }
  
  (fEvtStr)->fNTracks=0;
 
}
