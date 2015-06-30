#ifndef AliMFTConstants_H
#define AliMFTConstants_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Constants for the Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include <TObject.h>

class AliMFTConstants : public TObject {

public:
  // Geometry
  static const Int_t kNDisks = 5;
  static const Int_t kNPlanePerDisk = 2;
  static const Double_t kSensorLength;
  static const Double_t kSensorHeight;
  static const Double_t kSensorActiveHeight;
  static const Double_t kSensorActiveWidth;
  static const Double_t kSensorInterspace;
  static const Double_t kSensorSideOffset;
  static const Double_t kSensorTopOffset;
  static const Double_t kLadderOffsetToEnd;
  static const Double_t kSensorThickness;
  static const Double_t kFlexThickness; 
  static const Double_t kXPixelPitch;
  static const Double_t kYPixelPitch;
  static const Double_t kSensorMargin;

  static const Int_t kNPixelX = 1024;
  static const Int_t kNPixelY = 512;

  
  
  
  
  
  
  
  static const Int_t fNMaxPlanes = 20;

  static const Int_t fNMaxDigitsPerCluster = 50;  ///< max number of digits per cluster
  static const Double_t fCutForAvailableDigits;   ///<
  static const Double_t fCutForAttachingDigits;   ///<

  static const Int_t fNMaxMCTracksPerCluster = 10;   ///< max number of MC tracks sharing the same MFT cluster
  static const Int_t fNMaxMCTracksPerDigit = 3;      ///< max number of MC tracks sharing the same MFT digit

  static const Double_t fElossPerElectron;

  // superposition between the active elements tasselling the MFT planes, for having a full acceptance coverage even in case of 10 degrees inclined tracks
  static const Double_t fActiveSuperposition;  ///<
                                                
  static const Double_t fHeightActive;   ///< height of the active elements
  static const Double_t fHeightReadout;  ///< height of the readout elements attached to the active ones
	 
  // minimum border size between the end of the support plane and the sensors: fHeightReadout + 0.3
  static const Double_t fSupportExtMargin;  ///<

  static const Int_t fNMaxDetElemPerPlane = 1000;  ///<

  static const Double_t fRadLengthSi;    ///< expressed in cm

  static const Double_t fWidthChip;      ///< expressed in cm

  static const Double_t fPrecisionPointOfClosestApproach;  ///< precision (along z) for the research of the point of closest approach for a dimuon

  static const Double_t fZEvalKinem;     // z coordinate at which the kinematics is evaluated for the ESD and AOD tracks

  static const Double_t fXVertexTolerance;   // tolerance on the vertex for the first extrapolation of MUON tracks to I.P.
  static const Double_t fYVertexTolerance;   // tolerance on the vertex for the first extrapolation of MUON tracks to I.P.

  static const Double_t fPrimaryVertexResX;   // expected res. in Pb-Pb for the prim vtx from ITS+MFT combined vertexing (should not be used, depends on contributors)
  static const Double_t fPrimaryVertexResY;   // expected res. in Pb-Pb for the prim vtx from ITS+MFT combined vertexing (should not be used, depends on contributors)
  static const Double_t fPrimaryVertexResZ;   // expected res. in Pb-Pb for the prim vtx from ITS+MFT combined vertexing (should not be used, depends on contributors)

  static const Double_t fMisalignmentMagnitude;   // Expected misalignment magnitude (for MC, waiting for OCDB)

  static const Int_t fNMaxMuonsForPCA = 10;
  static const Int_t fNMaxPileUpEvents = 5;         // Max events to consider for pile-up studies
  static const Int_t fLabelOffsetMC = 10000000;     // Offset to be added to MC labels of tracks from merged underlying and pile-up events
  static const Int_t fNMaxLaddersPerPlane = 20;
  static const Int_t fNMaxChipsPerLadder = 5;
  static const Int_t fNMFTHalves = 2;
  static const Double_t fChipWidth;
  static const Double_t fMinDistanceLadderFromSupportRMin;
  static const Double_t fChipThickness;
  static const Double_t fChipInterspace; // Offset between two adjacent chip on a ladder
  static const Double_t fChipSideOffset; // Side Offset between the ladder edge and the chip edge
  static const Double_t fChipTopOffset; // Top Offset between the ladder edge and the chip edge


protected:

  AliMFTConstants() : TObject() {}
  virtual ~AliMFTConstants(){}

  ClassDef(AliMFTConstants, 4)    // MFT global constants 

};
	
#endif

