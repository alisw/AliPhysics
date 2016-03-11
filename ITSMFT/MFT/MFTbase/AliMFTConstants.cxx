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

//====================================================================================================================================================
/// \class AliMFTConstants
//      Constants for the Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliMFTConstants.h"

/// \cond CLASSIMP
ClassImp(AliMFTConstants);
/// \endcond

// Geometry Related Constants

const Double_t AliMFTConstants::kSensorLength=3.; //[cm]
const Double_t AliMFTConstants::kSensorHeight=1.5; //[cm]
const Double_t AliMFTConstants::kXPixelPitch=29.250e-4; // 29.15 micron // TODO : Check that
const Double_t AliMFTConstants::kYPixelPitch=26.880e-4; // 26.88 micron // TODO : Check that
const Double_t AliMFTConstants::kSensorMargin=29.120e-4; // 29.12 micron // TODO : Check that

const Double_t AliMFTConstants::kSensorActiveWidth  = kNPixelX * kXPixelPitch; //[cm]
const Double_t AliMFTConstants::kSensorActiveHeight = kNPixelY * kYPixelPitch; //[cm]

const Double_t AliMFTConstants::kSensorInterspace = 0.01; //[cm]  Offset between two adjacent chip on a ladder
const Double_t AliMFTConstants::kSensorSideOffset=0.04; // [cm] Side Offset between the ladder edge and the chip edge
const Double_t AliMFTConstants::kSensorTopOffset=0.04; // [cm] Top Offset between the ladder edge and the chip edge
const Double_t AliMFTConstants::kLadderOffsetToEnd=3.; // [cm] Offset between the last Chip of the ladder and the end of the ladder toward the DAQ connector
const Double_t AliMFTConstants::kSensorThickness=50.e-4; // 50 micron


const Double_t AliMFTConstants::kFlexThickness=165.e-4; // 100 micron // TODO : Change that

// Defaults parameters for track reconstruction
Double_t AliMFTConstants::fgDiskThicknessInX0[AliMFTConstants::kNDisks] = {0.008, 0.008, 0.008, 0.008, 0.008};
Double_t AliMFTConstants::fgPlaneZPos[2*AliMFTConstants::kNDisks] = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -68.0, -69.4, -76.1, -77.5};


const Double_t AliMFTConstants::fCutForAvailableDigits = 5.;
const Double_t AliMFTConstants::fCutForAttachingDigits = 1.;

const Double_t AliMFTConstants::fElossPerElectron = 3.62e-09;

const Double_t AliMFTConstants::fActiveSuperposition = 0.05;
                                 
const Double_t AliMFTConstants::fHeightActive = 1.3;
const Double_t AliMFTConstants::fHeightReadout = 0.2;

const Double_t AliMFTConstants::fSupportExtMargin = fHeightReadout + 0.3;

const Double_t AliMFTConstants::fRadLengthSi = 9.37;

const Double_t AliMFTConstants::fWidthChip = 1.0;

const Double_t AliMFTConstants::fPrecisionPointOfClosestApproach = 10.e-4;  // 10 micron

const Double_t AliMFTConstants::fZEvalKinem = 0.;

const Double_t AliMFTConstants::fXVertexTolerance = 500.e-4;    // 500 micron
const Double_t AliMFTConstants::fYVertexTolerance = 500.e-4;    // 500 micron

const Double_t AliMFTConstants::fPrimaryVertexResX = 5.e-4;   // 5 micron
const Double_t AliMFTConstants::fPrimaryVertexResY = 5.e-4;   // 5 micron
const Double_t AliMFTConstants::fPrimaryVertexResZ = 5.e-4;   // 5 micron

const Double_t AliMFTConstants::fMisalignmentMagnitude = 15.e-4;    // 15 micron

const Double_t AliMFTConstants::fChipWidth = 3.; // 3 cm ???
const Double_t AliMFTConstants::fChipThickness=500.e-4; // 50 micron
const Double_t AliMFTConstants::fMinDistanceLadderFromSupportRMin = 0.1; // 1mm ???

const Double_t AliMFTConstants::fChipInterspace=500.e-4; // 50um // Offset between two adjacent chip on a ladder
const Double_t AliMFTConstants::fChipSideOffset=500.e-4; // Side Offset between the ladder edge and the chip edge
const Double_t AliMFTConstants::fChipTopOffset=500.e-4; // Top Offset between the ladder edge and the chip edge

//====================================================================================================================================================
