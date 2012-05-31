// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <Rtypes.h>
#include "TEveElement.h"
#include "TPolyMarker3D.h"
#include "AliPID.h"

// --------------------------------------------------------------------
// --                             ITS                                --
// --------------------------------------------------------------------

// --------------------------------------------------------------------
struct AliHLTITSSpacePointData {
  Float_t fY;         // Y coordinate in local coordinates  GetY()
  Float_t fZ;         // Z coordinate in local coordinates  GetZ()
  Float_t fSigmaY2;   // error  of the clusters             GetSigmaY2()
  Float_t fSigmaZ2;   // error  of the clusters             GetSigmaZ2()
  Float_t fSigmaYZ;   // error  of the clusters             GetSigmaYZ()
  Float_t fQ;         // Q of cluster (in ADC counts)       GetQ()
  Int_t fNy;          //number of digits in Y direction     GetNy()
  Int_t fNz;          //number of digits in Z direction     GetNz()
  Int_t fLayer;       // layer number                       GetLayer()
  Int_t fTracks[3];   // MC label                           GetLabel(i)
  Int_t fIndex;       // Detector Index                     GetDetectorIndex()
};

// --------------------------------------------------------------------
struct AliHLTITSClusterData {
  Int_t fSpacePointCnt;
  AliHLTITSSpacePointData fSpacePoints[1];
};

// --------------------------------------------------------------------
// --                             TPC                                --
// --------------------------------------------------------------------

// --------------------------------------------------------------------
struct AliHLTTPCSpacePointData {
  Float_t fX;       // == fPadRow in local system
  Float_t fY;
  Float_t fZ;
  UInt_t fID;       //contains slice patch and number
  UChar_t fPadRow;
  Float_t fSigmaY2; //error (former width) of the clusters
  Float_t fSigmaZ2; //error (former width) of the clusters
  UInt_t fCharge;
  UInt_t fQMax;     // QMax of cluster
  Bool_t fUsed;     // only used in AliHLTTPCDisplay
  Int_t fTrackN;    // only used in AliHLTTPCDisplay
};

// --------------------------------------------------------------------
struct AliHLTTPCClusterData {
  Int_t fSpacePointCnt;
  AliHLTTPCSpacePointData fSpacePoints[1];
};

// --------------------------------------------------------------------
// --                            PHOS                                --
// --------------------------------------------------------------------

// struct AliHLTCaloClusterDataStruct
// {

//   UInt_t fNCells;                                //COMMENT
//   Float_t fGlobalPos[3];                      //COMMENT
//   Float_t fEnergy;                            //COMMENT
//   Float_t fTOF;                               //COMMENT
//   Float_t fDispersion;                        //COMMENT
//   Float_t fFitQuality;                        //COMMENT
//   Float_t fM20;                               //COMMENT
//   Float_t fM02;                               //COMMENT
//   Float_t fEmcCpvDistance;                    //COMMENT
//   Float_t fDistToBadChannel;                  //COMMENT
//   Float_t fPID[AliPID::kSPECIESN];            //COMMENT
//   Int_t fID;                                     //COMMENT
//   UChar_t fNExMax;                               //COMMENT 
//   Char_t fClusterType;                           //COMMENT
//   Float_t fDistanceToBadChannel;              //COMMENT
//   UShort_t fCellsAbsId;                      //COMMENT
//   Float_t fCellsAmpFraction;              //COMMENT
// };

// struct AliHLTCaloClusterHeaderStruct
// {
//   Short_t fNClusters;
// };

// --------------------------------------------------------------------
// struct AliHLTPHOSDigitDataStruct {
//   /** The x coordinate */
//   Float_t fX;

//   /** The x coordinate */
//   Float_t fZ;

//   /** The module number */
//   Int_t fModule;

//   /** The amplitude in ADC counts */
//   Float_t fAmplitude;

//   /** The time in sample count */ 
//   Float_t fTime;

//   /* The energy in GeV */
//   Float_t fEnergy;

//   /** The gain */
//   Int_t fGain;
  
//   /** The crazyness */
//   Int_t fCrazyness; 

//   /**  The baseline */
//   Float_t fBaseline;

// };

// --------------------------------------------------------------------
// --                            MUON                                --
// --------------------------------------------------------------------

// --------------------------------------------------------------------
struct AliHLTMUONTrackerRawData {
  UInt_t fDataId;
  UShort_t fADC;
};

// --------------------------------------------------------------------
struct AliHLTMUONTrackerMappingData {
  Int_t fDetElemId;
  Int_t fIX,fIY;
  Int_t fCath;
  Float_t fRealX,fRealY,fRealZ;
  Int_t fADC;
  Float_t fPed;
  Float_t fSigma;
};

// --------------------------------------------------------------------
struct AliHLTMUONTriggerPointData {
  Int_t fDetElemId;
  Float_t fX,fY,fZ ;
};

// --------------------------------------------------------------------
struct AliHLTMUONTriggerMappingData {
  AliHLTMUONTriggerPointData fLut[8][16][4][2][16];
};

// --------------------------------------------------------------------
// --                            STRUCT                              --
// --------------------------------------------------------------------

// --------------------------------------------------------------------
void hlt_structs() {}
