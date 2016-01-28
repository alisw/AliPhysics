#ifndef AliMFTConstants_H
#define AliMFTConstants_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTbase
/// \class AliMFTConstants
/// \brief Constants for the Muon Forward Tracker
///
/// Distance Unit is [cm]
///
/// \author Antonio Uras <antonio.uras@cern.ch>

#include <TObject.h>

class AliMFTConstants : public TObject {

public:
  // Geometry
  static const Int_t kNDisks = 5;             ///< \brief Number of Disk
  static const Double_t kSensorLength;        ///< \brief CMOS Sensor Length
  static const Double_t kSensorHeight;        ///< \brief CMOS Sensor Height
  static const Double_t kSensorActiveHeight;  ///< \brief CMOS Sensor Active height
  static const Double_t kSensorActiveWidth;   ///< \brief CMOS Sensor Active width
  static const Double_t kSensorThickness;     ///< \brief CMOS Sensor Thickness
  static const Double_t kXPixelPitch;         ///< \brief Pixel pitch along X
  static const Double_t kYPixelPitch;         ///< \brief Pixel pitch along Y
  static const Int_t kNPixelX = 1024;         ///< \brief Number of Pixel along X
  static const Int_t kNPixelY = 512;          ///< \brief Number of Pixel along Y
  static const Double_t kSensorMargin;        ///< \brief Inactive margin around active area

  static const Double_t kSensorInterspace;    ///< \brief Interspace between 2 sensors on a ladder
  static const Double_t kSensorSideOffset;    ///< \brief Offset of sensor compare to ladder edge (close to the beam pipe)
  static const Double_t kSensorTopOffset;     ///< \brief Offset of sensor compare to ladder top edge
  static const Double_t kLadderOffsetToEnd;   ///< \brief Offset of sensor compare to ladder connector edge
  static const Double_t kFlexThickness;       ///< \brief Flex Thickness

  
	
	
	/// Return Disk thickness in X0
	static Double_t DiskThicknessInX0(Int_t Id) {return (Id >= 0 && Id < kNDisks) ? fgDiskThicknessInX0[Id] : 0.;}
	
	/// Return Plane Z position
	static Double_t DefaultPlaneZ(Int_t Id) {return (Id >= 0 && Id < kNDisks*2) ? fgPlaneZPos[Id]+(-(Id%2)*2-1)*0.0025 : 0.;} // Move to the middle of the CMOS sensor in Z direction

  
  
  
  
  
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
	static Double_t  fgDiskThicknessInX0[kNDisks]; ///< default disk thickness in X0 for reconstruction
	static Double_t  fgPlaneZPos[2*kNDisks]; ///< default Plane Z position for reconstruction

  /// \cond CLASSIMP
  ClassDef(AliMFTConstants, 4);
  /// \endcond

};
	
#endif

