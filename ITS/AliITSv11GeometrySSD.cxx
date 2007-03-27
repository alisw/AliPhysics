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


//*************************************************************************
// SSD geometry, based on ROOT geometrical modeler
//
// Enrico Cattaruzza                                    ecattar@ts.infn.it
//*************************************************************************
#include "TMath.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoMatrix.h"
#include <TGeoManager.h>
#include "AliITSv11GeometrySSD.h"
#include "TVector3.h"
#include "TGeoArb8.h"
#include "TList.h"
#include "TGeoMatrix.h"
#include "TGeoCompositeShape.h"
#include "TGeoTube.h"
#include "TGeoBBox.h"
#include "TGeoTrd1.h"
#include <iostream>
using namespace std;
/////////////////////////////////////////////////////////////////////////////////
//Parameters for SSD Geometry
/////////////////////////////////////////////////////////////////////////////////
// Layer5 (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDLay5LadderLength      = 950.7;
const Int_t AliITSv11GeometrySSD::fgkSSDLay5SensorsNumber        =  22;
const Int_t AliITSv11GeometrySSD::fgkSSDLay5LadderNumber         =  34;
const Double_t AliITSv11GeometrySSD::fgkSSDLay5RadiusMin         = 378.0;
const Double_t AliITSv11GeometrySSD::fgkSSDLay5RadiusMax         = 384.0;
const Double_t AliITSv11GeometrySSD::fgkLay5CenterITSPosition    = 467.85;
/////////////////////////////////////////////////////////////////////////////////
// Layer6 (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDLay6LadderLength      = 1068.0;
const Int_t AliITSv11GeometrySSD::fgkSSDLay6SensorsNumber        =   25;
const Int_t AliITSv11GeometrySSD::fgkSSDLay6LadderNumber         =   38;
const Double_t AliITSv11GeometrySSD::fgkSSDLay6RadiusMin         =  428.0;
const Double_t AliITSv11GeometrySSD::fgkSSDLay6RadiusMax         =  434.0;
const Double_t AliITSv11GeometrySSD::fgkLay6CenterITSPosition    = 526.50;
/////////////////////////////////////////////////////////////////////////////////
// SSD Chips and Hybrid (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Int_t AliITSv11GeometrySSD::fgkSSDChipNumber               =   6;
const Double_t AliITSv11GeometrySSD::fgkSSDChipLength            =  11.100; 
const Double_t AliITSv11GeometrySSD::fgkSSDChipWidth             =   3.850;
const Double_t AliITSv11GeometrySSD::fgkSSDChipHeight            =   0.180;
const Double_t AliITSv11GeometrySSD::fgkSSDChipSeparationLength  =   1.000;
const Double_t AliITSv11GeometrySSD::fgkSSDChipGlueLength        =   fgkSSDChipLength;
const Double_t AliITSv11GeometrySSD::fgkSSDChipGlueWidth         =   fgkSSDChipWidth;
const Double_t AliITSv11GeometrySSD::fgkSSDChipGlueHeight        =   0.030;
/////////////////////////////////////////////////////////////////////////////////
// Stiffener (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerLength       =  73.000;
const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerWidth        =   6.500;
//const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerHeight       =   3.315;
const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerHeight       =   0.315;
const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerToChipDist   =   2.500;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor0603Length   =   1.600;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor0603Width    =   0.870;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor0603Height   =   0.800;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor1812Length   =   4.600;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor1812Width    =   3.400;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor1812Height   =   1.400;
const Double_t AliITSv11GeometrySSD::fgkSSDWireLength            =  30.000;
const Double_t AliITSv11GeometrySSD::fgkSSDWireRadius            =   0.185;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorPosition[2]  = {44.32, 0.33};
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorSeparation   = 0.44;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorLength       = 2.16;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorWidth        = 3.60;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorHeight       = 
													  0.25*fgkSSDStiffenerHeight;
/////////////////////////////////////////////////////////////////////////////////
// Cooling Block (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockLength    =   3.000;
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockWidth     =   4.000;
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockHeight[3] = 
													 {1.950, 0.240, 0.300};
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockHoleRadius[2] = 
														    {1.000, 0.120};
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockHoleLength[2] = 
															{1.900, 0.400};
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockHoleCenter    =  
																	 1.500;
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockHoleHeight    =  
																	 0.300;
/////////////////////////////////////////////////////////////////////////////////
// SSD Sensor (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const char*  AliITSv11GeometrySSD::fgkSSDSensitiveVolName       = 
                                                          "SSDSensorSensitiveVol";
const Double_t AliITSv11GeometrySSD::fgkSSDSensorLength          =  42.000;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorHeight          =   0.300;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorWidth           =  75.000;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorOverlap         = 
	   											   fgkSSDSensorLength-39.1;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorInsensitiveLength      = 1.; 
const Double_t AliITSv11GeometrySSD::fgkSSDSensorInsensitiveWidth       = 1.;
/////////////////////////////////////////////////////////////////////////////////
// Flex (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDFlexFullLength       =  106.000; 
const Double_t AliITSv11GeometrySSD::fgkSSDFlexLength[4]        = 
			{0.5*(fgkSSDStiffenerLength+fgkSSDChipNumber*fgkSSDChipLength
						+(fgkSSDChipNumber-1)*fgkSSDChipSeparationLength),
			 0.5*(fgkSSDStiffenerLength+fgkSSDChipNumber*fgkSSDChipLength
						+(fgkSSDChipNumber-1)*fgkSSDChipSeparationLength)-4.000,
						  9.500, 10.000};
const Double_t AliITSv11GeometrySSD::fgkSSDFlexWidth[2]         = 
														 {  9.340,  5.380};
const Double_t AliITSv11GeometrySSD::fgkSSDFlexHeight[2]        =
														 {  0.030,  0.020};      
const Double_t AliITSv11GeometrySSD::fgkSSDFlexAngle            =   30.000;
const Double_t AliITSv11GeometrySSD::fgkSSDFlexHoleLength       =    1.430;
const Double_t AliITSv11GeometrySSD::fgkSSDFlexHoleWidth        =    3.000;
const Double_t AliITSv11GeometrySSD::fgkSSDEndFlexCompLength[6] = 
										   {3.30,4.12,4.22,1.70,0.75,7.18};
const Double_t AliITSv11GeometrySSD:: fgkSSDEndFlexCompWidth[3] =
													   {15.03,23.48,12.28};
/////////////////////////////////////////////////////////////////////////////////
// SSD Ladder Cable (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDLadderCableWidth     =   23.5;
const Double_t AliITSv11GeometrySSD::fgkSSDEndLadderCableLength =   50.000; /////to be modified
/////////////////////////////////////////////////////////////////////////////////
// SSD Module (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDModuleStiffenerPosition[2]  = 
														   { 1.000, 3.900};
const Double_t AliITSv11GeometrySSD::fgkSSDModuleSensorSupportDistance =  
																	45.600;
const Double_t AliITSv11GeometrySSD::fgkSSDModuleCoolingBlockToSensor  =  
																	 5.075;
/////////////////////////////////////////////////////////////////////////////////
// Sensor Support (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportLength		   = 
																	   5.800;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportWidth          =  
																	   2.000;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportHeight[2]      =
														     { 4.620, 5.180};
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportThickness[2] = 
														     { 0.450, 0.450};
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportPosition       =  
									  0.5*(fgkSSDModuleSensorSupportDistance
							       +    fgkSSDSensorSideSupportThickness[0])
								   -           fgkSSDSensorSideSupportLength;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportLength	   =  
									   								   5.250;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportWidth        =
																       1.680;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportHeight[2]    =
	   {fgkSSDSensorSideSupportHeight[0]+fgkSSDSensorSideSupportThickness[0],
	   fgkSSDSensorSideSupportHeight[1]+fgkSSDSensorSideSupportThickness[1]};
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportThickness[2] =
   {fgkSSDSensorSideSupportThickness[0],fgkSSDSensorSideSupportThickness[1]};
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportPosition     = 
																      19.000;
/////////////////////////////////////////////////////////////////////////////////
// Chip Cables (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDChipCablesLength[2]   = 
						 {73.12/fgkSSDChipNumber,fgkSSDChipLength+2.*0.19};
const Double_t AliITSv11GeometrySSD::fgkSSDChipCablesHeight[4]   = 
						{  0.014,  0.010, fgkSSDModuleCoolingBlockToSensor
  -  (fgkSSDSensorSideSupportHeight[1]-fgkSSDSensorSideSupportHeight[0])
  -   fgkSSDCoolingBlockHoleCenter-fgkSSDStiffenerHeight
  -   fgkSSDChipHeight-fgkSSDSensorHeight,
      fgkSSDModuleCoolingBlockToSensor
  -   fgkSSDCoolingBlockHoleCenter-fgkSSDStiffenerHeight
  -   fgkSSDChipHeight-fgkSSDSensorHeight};
const Double_t AliITSv11GeometrySSD::fgkSSDChipCablesWidth[3]    = 
												 { 11.000,  0.800,  0.600};
/////////////////////////////////////////////////////////////////////////////////
// Carbon Fiber Junction Parameters (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberJunctionLength          = 
																	   3.820;
//const Double_t AliITSv11GeometrySSD::fgkCarbonFiberJunctionLength          = 
//																	   3.780;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberJunctionWidth           = 
									  fgkSSDSensorLength-fgkSSDSensorOverlap;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberJunctionEdge[2]         = 
															 {  0.86,  0.30};
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberJunctionAngle[2]        = 
															 { 30.00, 90.00};
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberJunctionToSensorSupport = 
																	    1.78;
/////////////////////////////////////////////////////////////////////////////////
//Carbon Fiber Parameters (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberTriangleLength          = 
  fgkSSDModuleSensorSupportDistance-2.*fgkCarbonFiberJunctionToSensorSupport;  
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberTriangleAngle           = 
																	   60.00;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberSupportTopEdgeDist[2]   = 
														   {  0.751,  0.482};
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberSupportEdgeLength       =  
																	   1.630;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberSupportWidth            = 
																	   0.950;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberSupportXAxisLength      = 
	     fgkCarbonFiberTriangleLength-0.5*fgkCarbonFiberSupportTopEdgeDist[1]
	             / TMath::Cos(fgkCarbonFiberTriangleAngle*TMath::DegToRad());
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberSupportYAxisLength      = 
				  0.5*(fgkCarbonFiberJunctionWidth-fgkCarbonFiberSupportWidth)
			 - fgkCarbonFiberSupportTopEdgeDist[0]-fgkCarbonFiberSupportWidth;
/////////////////////////////////////////////////////////////////////////////////
// Carbon Fiber Lower Support Parameters (lengths are in mm)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberLowerSupportWidth             
																	  =  0.950;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberLowerSupportLowerLenght       
																	  =  1.600;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberLowerSupportHeight            
																	  =  0.830;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberLowerSupportVolumeSeparation  
											  = 0.5*fgkCarbonFiberSupportWidth;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberLowerSupportTransverseWidth   
								 = fgkCarbonFiberJunctionWidth
                                 - 2.*(fgkCarbonFiberLowerSupportWidth
								 + fgkCarbonFiberLowerSupportVolumeSeparation);
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberLowerSupportVolumePosition[2] 
								 = {fgkCarbonFiberLowerSupportWidth
								 +  fgkCarbonFiberLowerSupportVolumeSeparation,
									fgkCarbonFiberLowerSupportWidth
								 +  fgkCarbonFiberLowerSupportVolumeSeparation							
								 +  fgkCarbonFiberLowerSupportTransverseWidth};
/////////////////////////////////////////////////////////////////////////////////
// End Ladder Carbon Fiber Junction Parameters (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkEndLadderCarbonFiberLowerJunctionLength[2] = 
  {0.5*(fgkSSDLay5LadderLength
	  -fgkSSDLay5SensorsNumber*fgkCarbonFiberJunctionWidth
	  -fgkCarbonFiberLowerSupportWidth),
       0.5*(fgkSSDLay5LadderLength
	  -fgkSSDLay5SensorsNumber*fgkCarbonFiberJunctionWidth
	  +fgkCarbonFiberLowerSupportWidth)};
const Double_t AliITSv11GeometrySSD::fgkEndLadderCarbonFiberUpperJunctionLength[2] = 
		{fgkEndLadderCarbonFiberLowerJunctionLength[0]-20.4,
		 fgkEndLadderCarbonFiberLowerJunctionLength[1]-20.6};
const Double_t AliITSv11GeometrySSD::fgkEndLadderMountingBlockPosition[2] = 
						   {fgkEndLadderCarbonFiberLowerJunctionLength[0]-16.50,
						   fgkEndLadderCarbonFiberLowerJunctionLength[1]-31.50};
/////////////////////////////////////////////////////////////////////////////////
// Cooling Tube Support (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportRmax          =  1.45;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportRmin          
											  = fgkSSDCoolingBlockHoleRadius[0];
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportLength        =  8.55;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportHeight        =  0.85;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportWidth         =  2.00;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportSeparation    = 
                                      fgkSSDSensorLength-2.*fgkSSDSensorOverlap;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportToCarbonFiber = 11.70;
/////////////////////////////////////////////////////////////////////////////////
// Cooling Tube (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeRmax       = 
													  fgkCoolingTubeSupportRmin;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeRmin       =  0.96;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeLength     = 
													fgkCarbonFiberJunctionWidth;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSeparation = 
					fgkSSDModuleSensorSupportDistance+fgkSSDCoolingBlockLength;
//const Double_t AliITSv11GeometrySSD_ct::fgkCoolingTubeLength               =  39.1;
/////////////////////////////////////////////////////////////////////////////////
// SSD Mounting Block Parameters (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockLength[3]            = 
															{ 60.0, 42.2, 34.0};
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHeight[4]            =
													  {  4.0,  8.0,  5.0,  0.2};
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockWidth                =   
																		   20.0;
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHoleTrapezoidAngle   =   
																		   40.0;
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHoleTrapezoidHeight  = 
	           0.30*(fgkSSDMountingBlockHeight[1]-fgkSSDMountingBlockHeight[2]);
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHoleTrapezoidUpBasis =    
																			2.5;
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHoleTubeLength[2]    = 
																  { 56.0, 12.0}; 
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHoleTubeWidth[2]     = 
																  {  5.0,  2.9}; 
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHoleRadius           = 
																			1.0;
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockScrewHoleEdge        =   
																			6.0;
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockScrewHoleHeigth      =  
																			4.0;
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockScrewHoleRadius[2]   =
									{  1.5,fgkSSDMountingBlockScrewHoleEdge/6.};
/////////////////////////////////////////////////////////////////////////////////
AliITSv11GeometrySSD::AliITSv11GeometrySSD(){           
  ///////////////////////////////// 
  // Standard Default Constructor
  ///////////////////////////////// 
  // Initializing display colors  
  ////////////////////////////// 
  fColorCarbonFiber =  4;
  fColorRyton       =  5;
  fColorPhynox      =  5;
  fColorSilicon     =  3;
  fColorAl          =  7;
  fColorKapton      =  6;
  fColorPolyhamide  =  5;
  fColorStiffener   =  9;
  fColorEpoxy       = 30;
  CreateMaterials();
}
/////////////////////////////////////////////////////////////////////////////////
// Setting the transformation Matrices
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetSSDSensorSupportCombiTransMatrix(){
  ////////////////////////////////////////////////////////
  //Translation Parameters SSDSensorSupportAssembly:
  ////////////////////////////////////////////////////////
  const Double_t SSDSensorSupportTransX[3] = {-0.5*fgkSSDSensorSideSupportWidth,
					                           0.5*fgkSSDSensorSideSupportWidth,
					                   0.5*fgkSSDSensorCenterSupportThickness[0]
									  -    fgkSSDSensorCenterSupportPosition}; 
  const Double_t SSDSensorSupportTransY[3] = 
									   {0.5*fgkSSDSensorSideSupportThickness[0],
					                   -0.5*fgkSSDSensorSideSupportThickness[0]
						               -fgkSSDModuleSensorSupportDistance,
					                    0.5*fgkSSDSensorCenterSupportWidth
									   -0.5*fgkSSDModuleSensorSupportDistance}; 
  const Double_t SSDSensorSupportTransZ[3] = {0.,0.,
										fgkSSDSensorCenterSupportThickness[0]}; 
  ////////////////////////////////////////////////////////
  //Rotational Parameters SSDSensorSupportAssembly:
  ////////////////////////////////////////////////////////  
  const Double_t SSDSensorSupportRotPhi[3]   = {   0., 180., 270.};
  const Double_t SSDSensorSupportRotTheta[3] = {  90.,  90.,  90.};
  const Double_t SSDSensorSupportRotPsi[3]   = {- 90.,- 90.,- 90.};
  ////////////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of SSDSensorSupportAssembly:
  ////////////////////////////////////////////////////////////////
  char SSDSensorSupportCombiTransName[40];
  char SSDSensorSupportRotName[40];
  TGeoCombiTrans *SSDSensorSupportLocalMatrix[fgkSSDSensorSupportCombiTransNumber];
  for(Int_t i=0; i<fgkSSDSensorSupportCombiTransNumber; i++){ 
		sprintf(SSDSensorSupportCombiTransName,"SSDSensorSupportCombiTrans%i",i);
		sprintf(SSDSensorSupportRotName,"SSDSensorSupportRot%i",i);
		SSDSensorSupportLocalMatrix[i] =
				 new TGeoCombiTrans(SSDSensorSupportCombiTransName,
									SSDSensorSupportTransX[i],
									SSDSensorSupportTransY[i],
									SSDSensorSupportTransZ[i],
									new TGeoRotation(SSDSensorSupportRotName,
													  SSDSensorSupportRotPhi[i],
													SSDSensorSupportRotTheta[i],
												    SSDSensorSupportRotPsi[i]));
           SSDSensorSupportCombiTransMatrix[i] = SSDSensorSupportLocalMatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetSSDModuleCombiTransMatrix(Double_t SSDChipCablesHeigth){
/////////////////////////////////////////////////////////////////////////////////
//Translation Parameters SSDModuleAssembly:
////////////////////////////////////////////////////////
  const Double_t SSDModuleTransX[7] = {0.5*fgkSSDStiffenerLength,
           				         0.5*fgkSSDChipLength+0.5*(fgkSSDStiffenerLength
				       - (fgkSSDChipNumber*fgkSSDChipLength+(fgkSSDChipNumber-1)
					   * fgkSSDChipSeparationLength)),
				       - fgkSSDModuleStiffenerPosition[0]+0.5*fgkSSDSensorWidth,
				   0.5*(fgkSSDStiffenerLength-fgkSSDChipNumber*fgkSSDChipLength
					   -       (fgkSSDChipNumber-1)*fgkSSDChipSeparationLength),
				0.5*fgkSSDStiffenerLength-0.5*fgkSSDModuleSensorSupportDistance
			   - fgkSSDCoolingBlockLength,
				0.5*(fgkSSDStiffenerLength+fgkSSDChipCablesLength[1]
			   + (fgkSSDChipNumber-1)*fgkSSDChipCablesLength[0]),
				0.5*(fgkSSDStiffenerLength+fgkSSDChipNumber*fgkSSDChipLength
			   +(fgkSSDChipNumber-1)*fgkSSDChipSeparationLength)}; 
  const Double_t SSDModuleTransY[7] = {0.5*fgkSSDStiffenerWidth,
				       0.5*fgkSSDChipWidth+(fgkSSDStiffenerWidth
					 - fgkSSDStiffenerToChipDist-fgkSSDChipWidth),
				     - fgkSSDModuleStiffenerPosition[1]+0.5*fgkSSDSensorLength,
				       fgkSSDStiffenerWidth,
				       0.,
				       fgkSSDSensorLength-2.*fgkSSDModuleStiffenerPosition[1]
					 - fgkSSDStiffenerWidth
				     + fgkSSDStiffenerToChipDist+fgkSSDChipWidth,
				       fgkSSDSensorLength-2.*fgkSSDModuleStiffenerPosition[1]
					 - fgkSSDStiffenerWidth}; 
  const Double_t SSDModuleTransZ[7] = {0.5*fgkSSDStiffenerHeight,
				     - 0.5*fgkSSDChipHeight,
				       0.5*fgkSSDSensorHeight-fgkSSDStiffenerHeight-fgkSSDChipHeight
					 - SSDChipCablesHeigth,
				     - 0.5*fgkSSDFlexHeight[0],
				       fgkSSDStiffenerHeight,
				     - fgkSSDChipHeight,
				     - 0.5*fgkSSDChipHeight}; 
  ////////////////////////////////////////////////////////
  //Rotational Parameters SSDModuleAssembly:
  ////////////////////////////////////////////////////////  
  const Double_t SSDModuleRotPhi[7]   = {   0.,   0.,  90.,   0.,   0.,  90., 180.};
  const Double_t SSDModuleRotTheta[7] = {   0.,   0.,   0.,   0.,   0.,   0.,   0.};
  const Double_t SSDModuleRotPsi[7]   = {   0.,   0.,   0.,   0.,   0.,   0.,   0.};
  ////////////////////////////////////////////////////////  
  //Name of CombiTrans Transformation of SSDModuleAssembly:
  ////////////////////////////////////////////////////////  
  const char* SSDModuleCombiTransName[7] = {"SSDStiffenerCombiTrans",
											     "SSDChipCombiTrans",
											   "SSDSensorCombiTrans",
											    "SSDFlex0CombiTrans",
										 "SSDCoolingBlockCombiTrans",
										   "SSDChipCablesCombiTrans",
										        "SSDFlex1CombiTrans"};
  const char* SSDModuleRotName[7] = {"SSDStiffenerRotName",
										  "SSDChipRotName",
										"SSDSensorRotName",
										 "SSDFlex0RotName",
								  "SSDCoolingBlockRotName",
									"SSDChipCablesRotName",
										"SSDFlex1RotName"};
  TGeoCombiTrans *SSDModuleLocalMatrix[fgkSSDModuleCombiTransNumber];
  for(Int_t i=0; i<fgkSSDModuleCombiTransNumber; i++){ 
					SSDModuleLocalMatrix[i] =
					new TGeoCombiTrans(SSDModuleCombiTransName[i],
									   SSDModuleTransX[i],
									   SSDModuleTransY[i],
									   SSDModuleTransZ[i],
									   new TGeoRotation(SSDModuleRotName[i],
														SSDModuleRotPhi[i],
												        SSDModuleRotTheta[i],
												        SSDModuleRotPsi[i]));
    SSDModuleCombiTransMatrix[i] = SSDModuleLocalMatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetCarbonFiberJunctionCombiTransMatrix(){
/////////////////////////////////////////////////////////////////////////////////
  //Translation Parameters CarbonFiberJunction:
  ////////////////////////////////////////////////////////
  const Double_t CarbonFiberJunctionTransX[3] = 
				{  0.0,fgkCarbonFiberTriangleLength,fgkCarbonFiberTriangleLength
				 * TMath::Cos(fgkCarbonFiberTriangleAngle*TMath::DegToRad())};
  const Double_t CarbonFiberJunctionTransY[3] = 
									  {  0.0, 0.0,fgkCarbonFiberTriangleLength
				 *   TMath::Sin(fgkCarbonFiberTriangleAngle*TMath::DegToRad())};
  const Double_t CarbonFiberJunctionTransZ[3] = {  0.0,  0.0,  0.0};
  ////////////////////////////////////////////////////////
  //Rotational Parameters CarbonFiberJunction:
  ////////////////////////////////////////////////////////
  const Double_t CarbonFiberJunctionRotPhi[3]   = {   0., 120., 240.};
  const Double_t CarbonFiberJunctionRotTheta[3] = {   0.,   0.,   0.};
  const Double_t CarbonFiberJunctionRotPsi[3]   = {   0.,   0.,   0.};
  ///////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of CarbonFiberJunction:
  ///////////////////////////////////////////////////////////
  char CarbonFiberJunctionCombiTransName[40];
  char CarbonFiberJunctionRotName[40];
  TGeoCombiTrans *CarbonFiberJunctionLocalMatrix[fgkCarbonFiberJunctionCombiTransNumber];
  for(Int_t i=0; i<fgkCarbonFiberJunctionCombiTransNumber; i++) {
		sprintf(CarbonFiberJunctionCombiTransName,"CarbonFiberJunctionCombiTrans%i",i);
		sprintf(CarbonFiberJunctionRotName,"CarbonFiberJunctionRot%i",i);
		CarbonFiberJunctionLocalMatrix[i] =
					new TGeoCombiTrans(CarbonFiberJunctionCombiTransName,
									   CarbonFiberJunctionTransX[i],
									   CarbonFiberJunctionTransY[i],
									   CarbonFiberJunctionTransZ[i],
								new TGeoRotation(CarbonFiberJunctionRotName,
												CarbonFiberJunctionRotPhi[i],
												CarbonFiberJunctionRotTheta[i],
												CarbonFiberJunctionRotPsi[i]));
    CarbonFiberJunctionCombiTransMatrix[i] = CarbonFiberJunctionLocalMatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetEndLadderCarbonFiberJunctionCombiTransMatrix(Int_t i){
/////////////////////////////////////////////////////////////////////////////////
  //Translation Parameters EndLadderCarbonFiberJunction:
  ////////////////////////////////////////////////////////
  const Double_t EndLadderCarbonFiberJunctionTransX[3] = 
				{  0.0,fgkCarbonFiberTriangleLength,fgkCarbonFiberTriangleLength
		*            TMath::Cos(fgkCarbonFiberTriangleAngle*TMath::DegToRad())};
  const Double_t EndLadderCarbonFiberJunctionTransY[3] = 
										{  0.0, 0.0,fgkCarbonFiberTriangleLength
		*            TMath::Sin(fgkCarbonFiberTriangleAngle*TMath::DegToRad())};
  const Double_t EndLadderCarbonFiberJunctionTransZ[3] = {  0.0,  0.0,  
						   0.5*(fgkEndLadderCarbonFiberLowerJunctionLength[i]
							   -fgkEndLadderCarbonFiberUpperJunctionLength[i])};
  ////////////////////////////////////////////////////////
  //Rotational Parameters EndLadderCarbonFiberJunction:
  ////////////////////////////////////////////////////////
  const Double_t EndLadderCarbonFiberJunctionRotPhi[3]   = {   0., 120., 240.};
  const Double_t EndLadderCarbonFiberJunctionRotTheta[3] = {   0.,   0.,   0.};
  const Double_t EndLadderCarbonFiberJunctionRotPsi[3]   = {   0.,   0.,   0.};
  ///////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of CarbonFiberJunction:
  ///////////////////////////////////////////////////////////
  char EndLadderCarbonFiberJunctionCombiTransName[40];
  char EndLadderCarbonFiberJunctionRotName[40];
  TGeoCombiTrans *EndLadderCarbonFiberJunctionLocalMatrix[fgkEndLadderCarbonFiberJunctionCombiTransNumber];
  for(Int_t i=0; i<fgkEndLadderCarbonFiberJunctionCombiTransNumber; i++) {
	sprintf(EndLadderCarbonFiberJunctionCombiTransName,"EndLadderCarbonFiberJunctionCombiTrans%i",i);
	sprintf(EndLadderCarbonFiberJunctionRotName,"EndLadderCarbonFiberJunctionRot%i",i);
	EndLadderCarbonFiberJunctionLocalMatrix[i] =
	new TGeoCombiTrans(EndLadderCarbonFiberJunctionCombiTransName,
										  EndLadderCarbonFiberJunctionTransX[i],
										  EndLadderCarbonFiberJunctionTransY[i],
										  EndLadderCarbonFiberJunctionTransZ[i],
						new TGeoRotation(EndLadderCarbonFiberJunctionRotName,
									  	  EndLadderCarbonFiberJunctionRotPhi[i],
										EndLadderCarbonFiberJunctionRotTheta[i],
									    EndLadderCarbonFiberJunctionRotPsi[i]));
    EndLadderCarbonFiberJunctionCombiTransMatrix[i] = 
									 EndLadderCarbonFiberJunctionLocalMatrix[i];
  }
}
////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetCarbonFiberAssemblyCombiTransMatrix(){
/////////////////////////////////////////////////////////////////////////////////
  //Translation Parameters CarbonFiberAssembly:
  ////////////////////////////////////////////////////////
  const Double_t CarbonFiberAssemblyTransX[3] = {  0.0,  0.0,  0.0};
  const Double_t CarbonFiberAssemblyTransY[3] = 
										{  0.5*fgkCarbonFiberJunctionWidth, 0.0, 
					fgkCarbonFiberJunctionWidth-fgkCarbonFiberLowerSupportWidth
			   -    fgkCarbonFiberLowerSupportVolumePosition[0]
			   -    fgkCarbonFiberLowerSupportVolumePosition[1]};
  const Double_t CarbonFiberAssemblyTransZ[3] = 
						  {  0.0,  0.0,-  0.5*fgkCarbonFiberLowerSupportHeight};
  ////////////////////////////////////////////////////////
  //Rotational Parameters CarbonFiberAssembly:
  ////////////////////////////////////////////////////////
  const Double_t CarbonFiberAssemblyRotPhi[3]   = {   0.,  90.,   0.};
  const Double_t CarbonFiberAssemblyRotTheta[3] = {  90.,
											-fgkCarbonFiberTriangleAngle,   0.};
  const Double_t CarbonFiberAssemblyRotPsi[3]   = {   0.,- 90.,   0.};
  ///////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of CarbonFiberAssembly:
  ///////////////////////////////////////////////////////////
  char CarbonFiberAssemblyCombiTransName[30];
  char CarbonFiberAssemblyRotName[30];
  TGeoCombiTrans *CarbonFiberAssemblyLocalMatrix[fgkCarbonFiberAssemblyCombiTransNumber];
  for(Int_t i=0; i<fgkCarbonFiberAssemblyCombiTransNumber; i++) {
	sprintf(CarbonFiberAssemblyCombiTransName,"CarbonFiberAssemblyCombiTrans%i",i);
	sprintf(CarbonFiberAssemblyRotName,"CarbonFiberAssemblyRot%i",i);
	CarbonFiberAssemblyLocalMatrix[i] =
						new TGeoCombiTrans(CarbonFiberAssemblyCombiTransName,
										   CarbonFiberAssemblyTransX[i],
										   CarbonFiberAssemblyTransY[i],
										   CarbonFiberAssemblyTransZ[i],
						 new TGeoRotation(CarbonFiberAssemblyRotName,
										   CarbonFiberAssemblyRotPhi[i],
										   CarbonFiberAssemblyRotTheta[i],
										   CarbonFiberAssemblyRotPsi[i]));
    CarbonFiberAssemblyCombiTransMatrix[i] = CarbonFiberAssemblyLocalMatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetCoolingTubeSupportCombiTransMatrix(){
  ////////////////////////////////////////////////////////
  //Translation Parameters CoolingTubeSupport:
  ////////////////////////////////////////////////////////
  Double_t Phi = TMath::ASin(0.5*fgkCoolingTubeSupportHeight
													/fgkCoolingTubeSupportRmax);
  const Double_t CoolingTubeSupportTransX[2] = 
						  {0.,2.*fgkCoolingTubeSupportRmax*TMath::Cos(Phi)
						+  2.*(fgkCoolingTubeSupportLength
						-  fgkCoolingTubeSupportRmax*(1.+TMath::Cos(Phi)))
						+  fgkCarbonFiberTriangleLength
						-  2.0*fgkCarbonFiberJunctionLength};
  const Double_t CoolingTubeSupportTransY[2] = {  0.0,  0.0};
  const Double_t CoolingTubeSupportTransZ[2] = {  0.0,  0.0};
  ////////////////////////////////////////////////////////
  //Rotational Parameters CoolingTubeSupport:
  ////////////////////////////////////////////////////////
  const Double_t CoolingTubeSupportRotPhi[2]   = {   0., 180.};
  const Double_t CoolingTubeSupportRotTheta[2] = {   0.,   0.};
  const Double_t CoolingTubeSupportRotPsi[2]   = {   0.,   0.};
  ///////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of CarbonFiberJunction:
  ///////////////////////////////////////////////////////////
  char CoolingTubeSupportCombiTransName[40];
  char CoolingTubeSupportRotName[40];
  TGeoCombiTrans *CoolingTubeSupportLocalMatrix[fgkCoolingTubeSupportCombiTransNumber];
  for(Int_t i=0; i<fgkCoolingTubeSupportCombiTransNumber; i++) {
	sprintf(CoolingTubeSupportCombiTransName,"CoolingTubeSupportCombiTrans%i",i);
	sprintf(CoolingTubeSupportRotName,"CoolingTubeSupportRot%i",i);
	CoolingTubeSupportLocalMatrix[i] =
			new TGeoCombiTrans(CoolingTubeSupportCombiTransName,
							   CoolingTubeSupportTransX[i],
							   CoolingTubeSupportTransY[i],
							   CoolingTubeSupportTransZ[i],
							   new TGeoRotation(CoolingTubeSupportRotName,
												CoolingTubeSupportRotPhi[i],
												CoolingTubeSupportRotTheta[i],
												CoolingTubeSupportRotPsi[i]));
    CoolingTubeSupportCombiTransMatrix[i] = CoolingTubeSupportLocalMatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetCoolingTubeCombiTransMatrix(){
  ////////////////////////////////////////////////////////
  //Translation Parameters CoolingTube:
  ////////////////////////////////////////////////////////
  const Double_t CoolingTubeTransX[2] = {  0.,  fgkCoolingTubeSeparation};
  const Double_t CoolingTubeTransY[2] = {  fgkCoolingTubeLength/2.0, fgkCoolingTubeLength/2.0};
  const Double_t CoolingTubeTransZ[2] = {  0.0,  0.};
  ////////////////////////////////////////////////////////
  //Rotational Parameters CoolingTube:
  ////////////////////////////////////////////////////////
  const Double_t CoolingTubeRotPhi[2]   = {   0.,   0.};
  const Double_t CoolingTubeRotTheta[2] = {  90.,  90.};
  const Double_t CoolingTubeRotPsi[2]   = {   0.,   0.};
  ///////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of CarbonFiberJunction:
  ///////////////////////////////////////////////////////////
  const char* CoolingTubeCombiTransName[fgkCoolingTubeCombiTransNumber] = 
							{"CoolingTubeCombiTrans0","CoolingTubeCombiTrans1"};
  const char* CoolingTubeRorationName[fgkCoolingTubeCombiTransNumber] = 
								{"CoolingTubeRotation0","CoolingTubeRotation1"};
  TGeoCombiTrans *CoolingTubeLocalMatrix[fgkCoolingTubeCombiTransNumber];
  for(Int_t i=0; i<fgkCoolingTubeCombiTransNumber; i++) {
	CoolingTubeLocalMatrix[i] =
			new TGeoCombiTrans(CoolingTubeCombiTransName[i],
							   CoolingTubeTransX[i],
							   CoolingTubeTransY[i],
							   CoolingTubeTransZ[i], 
							   new TGeoRotation(CoolingTubeRorationName[i],
							                    CoolingTubeRotPhi[i],
							                    CoolingTubeRotTheta[i],
									    CoolingTubeRotPsi[i]    ) );
    CoolingTubeTransMatrix[i] = CoolingTubeLocalMatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetLadderSegmentCombiTransMatrix(){
/////////////////////////////////////////////////////////////////////////////////
  //Translation Parameters LadderSegment:
  ////////////////////////////////////////////////////////
	const Double_t LadderSegmentTransX[fgkLadderSegmentCombiTransNumber] = {  0.,
	 -  0.5*(fgkSSDSensorWidth-fgkCarbonFiberTriangleLength),
			 fgkCarbonFiberTriangleLength+fgkCarbonFiberJunctionToSensorSupport,
			 fgkCarbonFiberJunctionLength-(fgkCoolingTubeSupportLength
	 -       fgkCoolingTubeSupportRmax),
		0.5*(fgkCarbonFiberTriangleLength-fgkCoolingTubeSeparation)}; 
	const Double_t LadderSegmentTransY[fgkLadderSegmentCombiTransNumber] = {  0.,
	 -      (2.*fgkSSDSensorLength-fgkSSDSensorOverlap)+
			 fgkSSDModuleStiffenerPosition[1]+fgkSSDStiffenerWidth
	 +		 0.5*fgkSSDFlexHoleLength+2.*fgkCarbonFiberJunctionWidth
	 -		 0.5*(fgkCarbonFiberLowerSupportWidth+fgkSSDSensorCenterSupportLength
	 -            fgkSSDSensorCenterSupportThickness[0]),
				  fgkCarbonFiberJunctionWidth-0.5*(fgkCarbonFiberLowerSupportWidth
	 +            fgkSSDSensorCenterSupportLength
	 -            fgkSSDSensorCenterSupportThickness[0])
	 -			  fgkSSDSensorCenterSupportPosition,
				  fgkCarbonFiberJunctionWidth-fgkCarbonFiberLowerSupportWidth
	 -			  fgkCoolingTubeSupportToCarbonFiber,
												 0.0};
	const Double_t LadderSegmentTransZ[fgkLadderSegmentCombiTransNumber] = {  0.,
	 -        (fgkSSDModuleCoolingBlockToSensor+0.5*fgkCoolingTubeSupportHeight
	 -         fgkSSDSensorHeight-fgkSSDChipCablesHeight[3]-fgkSSDChipHeight),
																			 0.,
	 -     0.5*fgkCoolingTubeSupportHeight,
	 -     0.5*fgkCoolingTubeSupportHeight};
//////////////////////////////////////////////////
  //Rotational Parameters LadderSegment:
  ////////////////////////////////////////////////////////
  const Double_t LadderSegmentRotPhi[fgkLadderSegmentCombiTransNumber]   =
													  {   0.,   0.,- 90.,   0.,  0.};
  const Double_t LadderSegmentRotTheta[fgkLadderSegmentCombiTransNumber] = 
													  {   0.,   0.,   0.,  90.,  0.};
  const Double_t LadderSegmentRotPsi[fgkLadderSegmentCombiTransNumber]   = 
													  {   0.,   0.,   0.,   0.,  0.};
  //////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of LadderSegment:
  //////////////////////////////////////////////////////
  char LadderSegmentCombiTransName[40];
  char LadderSegmentRotName[40];
  TGeoCombiTrans *LadderSegmentLocalMatrix[fgkLadderSegmentCombiTransNumber];
  for(Int_t i=0; i<fgkLadderSegmentCombiTransNumber; i++) {
		sprintf(LadderSegmentCombiTransName,"LadderSegmentCombiTrans%i",i);
		sprintf(LadderSegmentRotName,"LadderSegmentRot%i",i);
		LadderSegmentLocalMatrix[i] =
					new TGeoCombiTrans(LadderSegmentCombiTransName,
									   LadderSegmentTransX[i],
									   LadderSegmentTransY[i],
									   LadderSegmentTransZ[i],
									   new TGeoRotation(LadderSegmentRotName,
														LadderSegmentRotPhi[i],
														LadderSegmentRotTheta[i],
														LadderSegmentRotPsi[i]));
    LadderSegmentCombiTransMatrix[i] = LadderSegmentLocalMatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetEndLadderSegmentCombiTransMatrix(Int_t i){
/////////////////////////////////////////////////////////////////////////////////
  //Translation Parameters EndLadderSegment:
  ////////////////////////////////////////////////////////
  const Double_t EndLadderSegmentTransX[4] = {0.0,
											  0.0,
	         -  0.25*(fgkSSDMountingBlockLength[0]
			 +	fgkSSDMountingBlockLength[1])
			 +  0.5*fgkCarbonFiberTriangleLength,
											  0.0}; 
  const Double_t EndLadderSegmentTransY[4] = 
							 {0.5*fgkEndLadderCarbonFiberLowerJunctionLength[i],
					          i==0 ? 0. : fgkCarbonFiberLowerSupportWidth,
					                 fgkEndLadderMountingBlockPosition[i],
					           (1-i)*(fgkEndLadderMountingBlockPosition[i]
										 +  0.5*fgkSSDMountingBlockWidth)}; 
  const Double_t EndLadderSegmentTransZ[4] = {0.0,
											  0.0,
										-  fgkSSDMountingBlockHeight[1]
										+  0.5*fgkSSDMountingBlockHeight[0],
									    -  0.5*fgkCarbonFiberLowerSupportHeight}; 
  ////////////////////////////////////////////////////////
  //Rotational Parameters EndLadderSegment:
  ////////////////////////////////////////////////////////  
  const Double_t EndLadderSegmentRotPhi[4]   = {   0.,  90.,   0.,   0.};
  const Double_t EndLadderSegmentRotTheta[4] = {  90.,-fgkCarbonFiberTriangleAngle,
												   0.,   0.};
  const Double_t EndLadderSegmentRotPsi[4]   = {   0.,- 90.,   0.,   0.};
  ////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of EndLadderSegment:
  ////////////////////////////////////////////////////////
  char EndLadderSegmentCombiTransName[30];
  char EndLadderSegmentRotName[30];
  TGeoCombiTrans *EndLadderSegmentLocalMatrix[fgkEndLadderSegmentCombiTransNumber];
  for(Int_t i=0; i<fgkEndLadderSegmentCombiTransNumber; i++){ 
  		sprintf(EndLadderSegmentCombiTransName,"EndLadderSegmentCombiTrans%i",i);
		sprintf(EndLadderSegmentRotName,"EndLadderSegmentRot%i",i);
		EndLadderSegmentLocalMatrix[i] =
				new TGeoCombiTrans(EndLadderSegmentCombiTransName,
								   EndLadderSegmentTransX[i],
								   EndLadderSegmentTransY[i],
								   EndLadderSegmentTransZ[i],
								   new TGeoRotation(EndLadderSegmentRotName,
													EndLadderSegmentRotPhi[i],
													EndLadderSegmentRotTheta[i],
													EndLadderSegmentRotPsi[i]));
    EndLadderSegmentCombiTransMatrix[i] = EndLadderSegmentLocalMatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetLadderCableCombiTransMatrix(Int_t iLayer){
////////////////////////////////////////////////////////////////////////////////
// Translation Parameters for LadderCable
/////////////////////////////////////////
  Double_t SSDFlexRadius = fgkSSDStiffenerHeight+2*fgkSSDFlexHeight[0]
					     + fgkSSDFlexHeight[1];  
  Double_t SSDFlexRadiusMax = (fgkSSDFlexLength[3]-fgkSSDFlexLength[2])
							/  TMath::Tan(fgkSSDFlexAngle*TMath::DegToRad());
  Double_t SSDLadderCableTransX[3];
  SSDLadderCableTransX[0] = (SSDFlexRadiusMax-fgkSSDFlexHeight[1]-SSDFlexRadius)
						  *   TMath::Sin(2.*fgkSSDFlexAngle*TMath::DegToRad())
						  *	  TMath::Cos(2.*fgkSSDFlexAngle*TMath::DegToRad());
  SSDLadderCableTransX[1] = ((SSDFlexRadiusMax-fgkSSDFlexHeight[1]-SSDFlexRadius)
						  -     SSDLadderCableTransX[0]
						  /     TMath::Sin(2.*fgkSSDFlexAngle*TMath::DegToRad()))
						  *     TMath::Cos(fgkSSDFlexAngle*TMath::DegToRad());						
  SSDLadderCableTransX[2] = (fgkSSDFlexFullLength-2.*fgkSSDFlexAngle
						  *	  TMath::DegToRad()*SSDFlexRadiusMax
						  -     fgkSSDFlexLength[2]-TMath::Pi()
						  *	  fgkSSDStiffenerHeight-fgkSSDFlexLength[0]
						  -	  fgkSSDLadderCableWidth)
						  *	  TMath::Cos(2.*fgkSSDFlexAngle*TMath::DegToRad());
  Double_t SSDLadderCableTransZ[3] = {SSDLadderCableTransX[0]
						  *	TMath::Tan(2.*fgkSSDFlexAngle*TMath::DegToRad()),
							SSDLadderCableTransX[1]
						  *	TMath::Tan(fgkSSDFlexAngle*TMath::DegToRad()),
							SSDLadderCableTransX[2]
						  *	TMath::Tan(2.*fgkSSDFlexAngle*TMath::DegToRad())};	
  TGeoRotation* LocalLadderCableRot[3];	
  LocalLadderCableRot[0] = new TGeoRotation("LocalLadderCableRot0",90.,0.,0.);
  LocalLadderCableRot[1] = new TGeoRotation("LocalLadderCableRot1",90.,60.,-90.);
  LocalLadderCableRot[2] = new TGeoRotation("LocalLadderCableRot2");
  LocalLadderCableRot[2]->SetRotation((*LocalLadderCableRot[1])
						 *			  (*LocalLadderCableRot[0]));
////////////////////////////////////////////
// LocalLadderCableCombiTransMatrix
////////////////////////////////////////////
  const Int_t LocalLadderSideCablesNumber = 2;
  const Int_t LocalLadderCombiTransNumber = 5;
  TGeoCombiTrans** LocalLadderCableCombiTransMatrix[LocalLadderSideCablesNumber];
  for(Int_t i=0; i<LocalLadderSideCablesNumber; i++) 
	LocalLadderCableCombiTransMatrix[i] = 
							   new TGeoCombiTrans*[LocalLadderCombiTransNumber];
///////////////////////////////////////////
// Left Side Ladder Cables Transformations
///////////////////////////////////////////
  LocalLadderCableCombiTransMatrix[0][0]  =
						new TGeoCombiTrans(-0.5*fgkCarbonFiberTriangleLength,
						0.,0.,NULL);
  LocalLadderCableCombiTransMatrix[0][1] = LadderSegmentCombiTransMatrix[1];
  LocalLadderCableCombiTransMatrix[0][2] = 
						new TGeoCombiTrans(fgkSSDModuleStiffenerPosition[0],
										   fgkSSDModuleStiffenerPosition[1],0.,0);
  LocalLadderCableCombiTransMatrix[0][3] = SSDModuleCombiTransMatrix[6];
  LocalLadderCableCombiTransMatrix[0][4] = 
						new TGeoCombiTrans(-SSDLadderCableTransX[0]
						-     SSDLadderCableTransX[1]-SSDLadderCableTransX[2]
						+     fgkSSDFlexLength[0]-fgkSSDFlexLength[2],
							  0.,
							  0.5*fgkSSDFlexHeight[0]+2.*(fgkSSDFlexHeight[0]
						+	  fgkSSDFlexHeight[1])+2.*fgkSSDStiffenerHeight
						+     SSDLadderCableTransZ[0]-SSDLadderCableTransZ[1]
						+	  SSDLadderCableTransZ[2],LocalLadderCableRot[2]);
///////////////////////////////////////////
// Rigth Side Ladder Cables Transformations
///////////////////////////////////////////
 for(Int_t i=0; i<LocalLadderCombiTransNumber; i++)
   LocalLadderCableCombiTransMatrix[1][i] = 
			(i!=3 ? LocalLadderCableCombiTransMatrix[0][i]:
					SSDModuleCombiTransMatrix[3]); 	
///////////////////////////////////////////
// Setting LadderCableCombiTransMatrix
///////////////////////////////////////////
Int_t BeamAxisTrans[3] = {0,0,0};
switch(iLayer){
case 5: 
	BeamAxisTrans[0] = fgkSSDLay5SensorsNumber/2;
	BeamAxisTrans[1] = BeamAxisTrans[0]+1;
	BeamAxisTrans[2] = BeamAxisTrans[0]-1;
	break;
case 6:
	BeamAxisTrans[0] = (fgkSSDLay6SensorsNumber-1)/2;
	BeamAxisTrans[1] = BeamAxisTrans[0]+1;
	BeamAxisTrans[2] = BeamAxisTrans[0];
	break;
}
 TGeoHMatrix* LocalLadderCableHMatrix[LocalLadderSideCablesNumber];
 TGeoRotation* LadderCableRot;
 TGeoTranslation* LadderCableTrans;
 for(Int_t i=0; i<LocalLadderSideCablesNumber; i++){
	LocalLadderCableHMatrix[i] = new TGeoHMatrix();
	for(Int_t j=0; j<LocalLadderCombiTransNumber; j++){
		  LocalLadderCableHMatrix[i]->MultiplyLeft(
		  LocalLadderCableCombiTransMatrix[i][LocalLadderCombiTransNumber-j-1]);
	}
	LadderCableRot = new TGeoRotation();
	LadderCableRot->SetMatrix(LocalLadderCableHMatrix[i]->GetRotationMatrix());
    LadderCableTrans = new TGeoTranslation();
    Double_t* LadderCableTransVector = LocalLadderCableHMatrix[i]->GetTranslation();
    LadderCableTrans->SetTranslation(LadderCableTransVector[0],
									 LadderCableTransVector[1]
					+                (i==0 ? BeamAxisTrans[0] : 0.)
					*				 fgkCarbonFiberJunctionWidth,
									 LadderCableTransVector[2]);	 
	LadderCableCombiTransMatrix[i] = new TGeoCombiTrans(*LadderCableTrans,
												   *LadderCableRot);
	} 
	LadderCableCombiTransMatrix[2] = 
					AddTranslationToCombiTrans(LadderCableCombiTransMatrix[1],0.,
					BeamAxisTrans[1]*fgkCarbonFiberJunctionWidth,0.);
	LadderCableCombiTransMatrix[3] = 
					AddTranslationToCombiTrans(LadderCableCombiTransMatrix[0],0.,
					BeamAxisTrans[2]*fgkCarbonFiberJunctionWidth,0.);
 }
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDSensorSupportShape(Double_t length, 
							Double_t height,Double_t width,Double_t* thickness){
////////////////////////////////////////////////////////////////////////////////
  const Int_t VertexNumber = 4;
  const Int_t ShapesNumber = 2;
  Double_t Width[2] = {width,width};
  TVector3** VertexPosition[ShapesNumber];
  for(Int_t i = 0; i<ShapesNumber; i++) VertexPosition[i] = 
													new TVector3*[VertexNumber];
  //First Shape Vertex Positioning
  VertexPosition[0][0] = new TVector3();
  VertexPosition[0][1] = new TVector3(height);
  VertexPosition[0][2] = new TVector3(thickness[0]);
  VertexPosition[0][3] = new TVector3(*VertexPosition[0][1]);
  //Second Shape Vertex Positioning
  VertexPosition[1][0] = new TVector3(*VertexPosition[0][0]);
  VertexPosition[1][1] = new TVector3(length);
  VertexPosition[1][2] = new TVector3(thickness[1]);
  VertexPosition[1][3] = new TVector3(VertexPosition[1][1]->X());
  char* SSDSensorSupportShapeName[ShapesNumber] = 
							{"SSDSensorSupportShape1","SSDSensorSupportShape2"};
  TGeoArb8* SSDSensorSupportShape[ShapesNumber];
  for(Int_t i = 0; i< ShapesNumber; i++) SSDSensorSupportShape[i] = 
					   GetArbShape(VertexPosition[i],Width,i==0 ? 
													 thickness[1]: thickness[0],
												  SSDSensorSupportShapeName[i]);
  /////////////////////////////////////
  //Setting Translations and Rotations: 
  /////////////////////////////////////
  TGeoRotation* SSDSensorSupportShapeRot[2];
  SSDSensorSupportShapeRot[0] = 
					 new TGeoRotation("SSDSensorSupportShapeRot1",180.,0.,0.);
  SSDSensorSupportShapeRot[1] = 
					 new TGeoRotation("SSDSensorSupportShapeRot2",90.,90.,-90.);
  TGeoTranslation* SSDSensorSupportShapeTrans = 
					 new TGeoTranslation("SSDSensorSupportShapeTrans",0.,0.,
															  0.5*thickness[0]);
  TGeoCombiTrans* SSDSensorSupportCombiTrans = 
	  new TGeoCombiTrans("SSDSensorSupportCombiTrans",0.5*thickness[0],width,0.,
						  new TGeoRotation((*SSDSensorSupportShapeRot[1])
						*                  (*SSDSensorSupportShapeRot[0])));
  TGeoVolume* SSDSensorSupportCompVolume = 
						  new TGeoVolumeAssembly("SSDSensorSupportCompVolume");
  SSDSensorSupportCompVolume->AddNode(new TGeoVolume("SSDSensorSupportVolume1",
							SSDSensorSupportShape[0],fgkSSDSensorSupportMedium),1,
							SSDSensorSupportShapeTrans);
  SSDSensorSupportCompVolume->AddNode(new TGeoVolume("SSDSensorSupportVolume2",
							SSDSensorSupportShape[1],fgkSSDSensorSupportMedium),1,
							SSDSensorSupportCombiTrans);
  return SSDSensorSupportCompVolume;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDSensorSupport(Int_t VolumeKind, Int_t n){
////////////////////////////////////////////////////////////////////////////////
  TGeoVolume* SSDSensorSupport;
  Double_t SideSupportThickness[2] = {fgkSSDSensorSideSupportThickness[0],
										   fgkSSDSensorSideSupportThickness[1]};
  VolumeKind == 0 ? SSDSensorSupport = GetSSDSensorSupportShape(
								fgkSSDSensorSideSupportLength,
								fgkSSDSensorSideSupportHeight[(n==0 ? 0 : 1)],
								fgkSSDSensorSideSupportWidth,
								SideSupportThickness) :
    SSDSensorSupport = GetSSDSensorSupportShape(fgkSSDSensorCenterSupportLength,
						        fgkSSDSensorCenterSupportHeight[(n==0 ? 0 : 1)],
							    fgkSSDSensorCenterSupportWidth,
							    SideSupportThickness);
  return SSDSensorSupport;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDSensorSupportAssembly(Int_t n){
////////////////////////////////////////////////////////////////////////////////
  TGeoVolume* SSDSensorSupportAssembly = 
							 new TGeoVolumeAssembly("SSDSensorSupportAssembly");
  const Int_t VolumeNumber = 2;
  TGeoVolume* SSDSensorSupport[VolumeNumber];
  for(Int_t i=0; i<VolumeNumber; i++) SSDSensorSupport[i] = 
														GetSSDSensorSupport(i,n);
  SetSSDSensorSupportCombiTransMatrix();
  for(Int_t i=0; i<fgkSSDSensorSupportCombiTransNumber; i++) 
	SSDSensorSupportAssembly->AddNode((i<2 ? SSDSensorSupport[0]:
											 SSDSensorSupport[1]),
										 i+1,SSDSensorSupportCombiTransMatrix[i]);
  return SSDSensorSupportAssembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDModule(Int_t iChipCablesHeight){
////////////////////////////////////////////////////////////////////////////////
  TGeoVolume* SSDModuleVolume[fgkSSDModuleCombiTransNumber-1];
  SSDModuleVolume[0] = GetSSDStiffenerAssembly();
  SSDModuleVolume[1] = GetSSDChipAssembly();
  SSDModuleVolume[2] = GetSSDSensor();
  SSDModuleVolume[3] = GetSSDFlexAssembly();
  SSDModuleVolume[4] = GetSSDCoolingBlockAssembly();
  SSDModuleVolume[5] = GetSSDChipCablesAssembly(fgkSSDChipCablesHeight[iChipCablesHeight+2]);
  SetSSDModuleCombiTransMatrix(fgkSSDChipCablesHeight[iChipCablesHeight+2]);
  TGeoCombiTrans* SSDModuleGlobalCombiTrans = 
								   new TGeoCombiTrans("SSDModuleGlobalCombiTrans",
									   fgkSSDModuleStiffenerPosition[0],
									   fgkSSDModuleStiffenerPosition[1],0.,NULL);
  TGeoHMatrix* SSDModuleHMatrix[fgkSSDModuleCombiTransNumber];
  TGeoVolume* SSDModule = new TGeoVolumeAssembly("SSDModule");
  for(Int_t i=0; i<fgkSSDModuleCombiTransNumber; i++){ 
	SSDModuleHMatrix[i] = new TGeoHMatrix((*SSDModuleGlobalCombiTrans)
						*				  (*SSDModuleCombiTransMatrix[i]));
	SSDModule->AddNode(i==fgkSSDModuleCombiTransNumber-1 ? SSDModuleVolume[3] : 
						  SSDModuleVolume[i],i!=fgkSSDModuleCombiTransNumber-1?1:2,
						  SSDModuleHMatrix[i]);
  }
  return SSDModule;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDSensor(){
////////////////////////////////////////////////////////////////////////////////
  Double_t fgkSSDSensorSensitiveLength = fgkSSDSensorLength-2.*fgkSSDSensorInsensitiveLength;
  Double_t fgkSSDSensorSensitiveWidth  = fgkSSDSensorWidth-2.*fgkSSDSensorInsensitiveWidth;
  TGeoBBox* SSDSensorSensitiveShape = new TGeoBBox("SSDSensorSensitiveShape",
                                                0.5*fgkSSDSensorSensitiveLength,
                                                0.5*fgkSSDSensorSensitiveWidth,
                                                0.5*fgkSSDSensorHeight);
  TGeoVolume* SSDSensorSensitive = 
      new TGeoVolume(fgkSSDSensitiveVolName,SSDSensorSensitiveShape,fgkSSDSensorMedium);
  SSDSensorSensitive->SetLineColor(fColorSilicon);
  TGeoBBox* SSDSensorInsensitiveShape[2];
  SSDSensorInsensitiveShape[0] = new TGeoBBox("SSDSensorInSensitiveShape1",
                                                0.5*fgkSSDSensorLength,
                                                0.5*fgkSSDSensorInsensitiveWidth,
                                                0.5*fgkSSDSensorHeight);
  SSDSensorInsensitiveShape[1] = new TGeoBBox("SSDSensorInSensitiveShape2",
                                                0.5*fgkSSDSensorInsensitiveWidth,
                                                0.5*fgkSSDSensorSensitiveWidth,
                                                0.5*fgkSSDSensorHeight);
  const char* SSDSensorInsensitiveName[2] = {"SSDSensorInsensitive1",
                                             "SSDSensorInsensitive2"};
  TGeoVolume* SSDSensorInsensitive[2];
  for(Int_t i=0; i<2; i++){ SSDSensorInsensitive[i] = 
      new TGeoVolume(SSDSensorInsensitiveName[i],SSDSensorInsensitiveShape[i],
                     fgkSSDSensorMedium);
      SSDSensorInsensitive[i]->SetLineColor(fColorCarbonFiber);
  }
  TGeoVolume* SSDSensorInsensitiveVol = 
                              new TGeoVolumeAssembly("SSDSensorInsensitiveVol");
  for(Int_t i=0; i<4; i++) 
            SSDSensorInsensitiveVol->AddNode(i%2==0 ? SSDSensorInsensitive[0]:
            SSDSensorInsensitive[1],i<2?1:2,
                  new TGeoTranslation(0.5*(1.-TMath::Power(-1.,i))*(i==1? 1.:-1.)
      *     (SSDSensorSensitiveShape->GetDX()+SSDSensorInsensitiveShape[1]->GetDX()),
                        0.5*(1.+TMath::Power(-1.,i))*(i==0?-1.: 1.)
      *     (SSDSensorSensitiveShape->GetDY()+SSDSensorInsensitiveShape[0]->GetDY()),0.));            
  TGeoVolume* SSDSensor = new TGeoVolumeAssembly("SSDSensor");
  SSDSensor->AddNode(SSDSensorSensitive,1),SSDSensor->AddNode(SSDSensorInsensitiveVol,1);
  return SSDSensor;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDChipAssembly(){
////////////////////////////////////////////////////////////////////////////////
  const Int_t SSDChipRowNumber = 2;
  TGeoBBox* SSDChipCompShape[2];
  SSDChipCompShape[0] =  new TGeoBBox("SSDChipCompShape",
										0.5*fgkSSDChipLength,
										0.5*fgkSSDChipWidth,
										0.5*(fgkSSDChipHeight-fgkSSDChipGlueHeight));
  SSDChipCompShape[1] =  new TGeoBBox("SSDChipGlueCompShape",
										0.5*fgkSSDChipLength,
										0.5*fgkSSDChipWidth,
										0.5*fgkSSDChipGlueHeight);
  TGeoVolume* SSDChipComp[2];
  SSDChipComp[0] = new TGeoVolume("SSDChipComp",SSDChipCompShape[0],fgkSSDChipMedium);
  SSDChipComp[1] = new TGeoVolume("SSDChipGlueComp",SSDChipCompShape[1],
								  fgkSSDChipGlueMedium);
  SSDChipComp[0]->SetLineColor(fColorSilicon);  
  SSDChipComp[1]->SetLineColor(fColorEpoxy);
  TGeoTranslation* SSDChipCompTrans[2];
  SSDChipCompTrans[0] = new TGeoTranslation(0.,0.,-SSDChipCompShape[1]->GetDZ());
  SSDChipCompTrans[1] = new TGeoTranslation(0.,0.,SSDChipCompShape[0]->GetDZ());
  TGeoVolume* SSDChip = new TGeoVolumeAssembly("SSDChip");  
  for(Int_t i=0; i<2; i++) SSDChip->AddNode(SSDChipComp[i],1,SSDChipCompTrans[i]);
  Double_t SSDChipSeparation[2] = {fgkSSDChipLength+fgkSSDChipSeparationLength,
						  fgkSSDSensorLength-2.*fgkSSDModuleStiffenerPosition[1]
				   -  2.*(fgkSSDStiffenerWidth-fgkSSDStiffenerToChipDist
				   -  0.5*fgkSSDChipWidth)};
  TGeoVolume* SSDChipAssembly = new TGeoVolumeAssembly("SSDChipAssembly"); 
  for(Int_t i=0; i<SSDChipRowNumber; i++)
    for(Int_t j=0; j<fgkSSDChipNumber; j++) 
		SSDChipAssembly->AddNode(SSDChip,fgkSSDChipNumber*i+j+1,
		new TGeoTranslation(j*SSDChipSeparation[0],i*SSDChipSeparation[1],0.));
  return SSDChipAssembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDStiffenerAssembly(){
////////////////////////////////////////////////////////////////////////////////
  const Int_t SSDStiffenerNumber = 2;
  Double_t SSDStiffenerSeparation = fgkSSDSensorLength
								  - 2.*fgkSSDModuleStiffenerPosition[1]
								  -    fgkSSDStiffenerWidth;
  TGeoVolume* SSDStiffener = new TGeoVolumeAssembly("SSDStiffener");
////////////////////////////
// Stiffener Volumes
///////////////////////////
  const Int_t StiffenerBoxNumber = 6;
  TGeoBBox* SSDStiffenerBoxShapes[StiffenerBoxNumber];
  SSDStiffenerBoxShapes[0] = new TGeoBBox("SSDStiffenerBoxShape1",
										  0.5* fgkSSDStiffenerLength,
										  0.5* fgkSSDStiffenerWidth,
										  0.5*(fgkSSDStiffenerHeight
						   -                   fgkSSDConnectorHeight));
  SSDStiffenerBoxShapes[1] = new TGeoBBox("SSDStiffenerBoxShape2",
										  0.5*(fgkSSDConnectorPosition[0]
						   -              2.0* fgkSSDConnectorLength
						   -				   fgkSSDConnectorSeparation),
										  0.5* fgkSSDStiffenerWidth,
										  0.5* fgkSSDConnectorHeight);
  SSDStiffenerBoxShapes[2] = new TGeoBBox("SSDStiffenerBoxShape3",
										  0.5*(fgkSSDConnectorSeparation
						   +              2.*  fgkSSDConnectorLength),
										  0.5* fgkSSDConnectorPosition[1],
										  0.5* fgkSSDConnectorHeight);
  SSDStiffenerBoxShapes[3] = new TGeoBBox("SSDStiffenerBoxShape4",
											SSDStiffenerBoxShapes[2]->GetDX(),
										  0.5*(fgkSSDStiffenerWidth
						   -                   fgkSSDConnectorPosition[1]
						   -                   fgkSSDConnectorWidth),
										  0.5* fgkSSDConnectorHeight);
   SSDStiffenerBoxShapes[4] = new TGeoBBox("SSDStiffenerBoxShape5",
										  0.5* fgkSSDConnectorSeparation,
										  0.5* fgkSSDConnectorWidth,
										  0.5* fgkSSDConnectorHeight);
   SSDStiffenerBoxShapes[5] = new TGeoBBox("SSDStiffenerBoxShape6",
										  0.5*(fgkSSDStiffenerLength
							 -				   fgkSSDConnectorPosition[0]),
										  0.5* fgkSSDStiffenerWidth,
										  0.5* fgkSSDConnectorHeight);
  TGeoVolume* SSDStiffenerBox[StiffenerBoxNumber];
  char SSDStiffenerBoxName[30];
  for(Int_t i=0; i<StiffenerBoxNumber; i++){ 
	sprintf(SSDStiffenerBoxName,"SSDStiffenerBox%i",i+1);
    SSDStiffenerBox[i] = new TGeoVolume(SSDStiffenerBoxName,SSDStiffenerBoxShapes[i],
									  fgkSSDStiffenerMedium);
	SSDStiffenerBox[i]->SetLineColor(fColorStiffener);
  }
////////////////////////////
// Connector 
///////////////////////////
  TGeoBBox* SSDConnectorShape =  new TGeoBBox("SSDConnectorShape",
											 0.5*fgkSSDConnectorLength,
											 0.5*fgkSSDConnectorWidth,
											 0.5*fgkSSDConnectorHeight);
  TGeoVolume* SSDConnector    = new TGeoVolume("SSDConnector",SSDConnectorShape,
												fgkSSDStiffenerConnectorMedium); 
  SSDConnector->SetLineColor(fColorAl);
  const Int_t SSDConnectorNumber = 2;
  TGeoTranslation* SSDConnectorTrans[SSDConnectorNumber];
  SSDConnectorTrans[0] = new TGeoTranslation("SSDConnectorTrans1",
        -  SSDStiffenerBoxShapes[0]->GetDX()+fgkSSDConnectorPosition[0]
		-  fgkSSDConnectorSeparation-1.5*fgkSSDConnectorLength,
		   SSDStiffenerBoxShapes[0]->GetDY()-fgkSSDConnectorPosition[1]
		-  SSDConnectorShape->GetDY(),
		   SSDStiffenerBoxShapes[0]->GetDZ()+SSDConnectorShape->GetDZ());	
  SSDConnectorTrans[1] = new TGeoTranslation("SSDConnectorTrans2",
        -  SSDStiffenerBoxShapes[0]->GetDX()+fgkSSDConnectorPosition[0]
        -  0.5*fgkSSDConnectorLength,
           SSDStiffenerBoxShapes[0]->GetDY()-fgkSSDConnectorPosition[1]
		-  SSDConnectorShape->GetDY(),
		   SSDStiffenerBoxShapes[0]->GetDZ()+SSDConnectorShape->GetDZ());
  for(Int_t i=0; i<SSDConnectorNumber; i++)
			SSDStiffener->AddNode(SSDConnector,i+1,SSDConnectorTrans[i]);	
//////////////////////////////////////
// TGeoTranslation for Stiffener Boxes
//////////////////////////////////////
  TGeoTranslation* SSDStiffenerBoxTrans[StiffenerBoxNumber];
  SSDStiffenerBoxTrans[0] = new TGeoTranslation("SSDStiffenerBoxTrans1",0.,0.,0.);
  SSDStiffenerBoxTrans[1] = new TGeoTranslation("SSDStiffenerBoxTrans2",
		 - (SSDStiffenerBoxShapes[0]->GetDX()-SSDStiffenerBoxShapes[1]->GetDX()),
		    0.,
			SSDStiffenerBoxShapes[0]->GetDZ()+SSDStiffenerBoxShapes[1]->GetDZ());
  SSDStiffenerBoxTrans[2] = new TGeoTranslation("SSDStiffenerBoxTrans3",
         - SSDStiffenerBoxShapes[0]->GetDX()-SSDStiffenerBoxShapes[2]->GetDX()
		 + fgkSSDConnectorPosition[0],
		   SSDStiffenerBoxShapes[0]->GetDY()-SSDStiffenerBoxShapes[2]->GetDY(),
		   SSDStiffenerBoxShapes[0]->GetDZ()+SSDStiffenerBoxShapes[2]->GetDZ());
  SSDStiffenerBoxTrans[3] = new TGeoTranslation("SSDStiffenerBoxTrans4",
         - SSDStiffenerBoxShapes[0]->GetDX()-SSDStiffenerBoxShapes[3]->GetDX()
		 + fgkSSDConnectorPosition[0],
		 - SSDStiffenerBoxShapes[0]->GetDY()+SSDStiffenerBoxShapes[3]->GetDY(),
		   SSDStiffenerBoxShapes[0]->GetDZ()+SSDStiffenerBoxShapes[3]->GetDZ());
  SSDStiffenerBoxTrans[4] = new TGeoTranslation("SSDStiffenerBoxTrans5",
		 - SSDStiffenerBoxShapes[0]->GetDX()+fgkSSDConnectorPosition[0]
		 - 0.5*fgkSSDConnectorSeparation-2.*SSDConnectorShape->GetDX(),
		   SSDStiffenerBoxShapes[0]->GetDY()-fgkSSDConnectorPosition[1]
		 - SSDConnectorShape->GetDY(),
		   SSDStiffenerBoxShapes[0]->GetDZ()+SSDConnectorShape->GetDZ());
  SSDStiffenerBoxTrans[5] = new TGeoTranslation("SSDStiffenerBoxTrans6",
         - SSDStiffenerBoxShapes[0]->GetDX()+fgkSSDConnectorPosition[0]
		 + SSDStiffenerBoxShapes[5]->GetDX(),
		   0.,
		   SSDStiffenerBoxShapes[0]->GetDZ()+SSDStiffenerBoxShapes[5]->GetDZ());
  for(Int_t i=0; i<StiffenerBoxNumber; i++) 
		SSDStiffener->AddNode(SSDStiffenerBox[i],1,SSDStiffenerBoxTrans[i]);	
  TGeoCombiTrans* SSDStiffenerCombiTrans[SSDStiffenerNumber];
  char SSDStiffenerCombiTransName[30];
    for(Int_t i=0; i<SSDStiffenerNumber; i++){ 
	sprintf(SSDStiffenerCombiTransName,"SSDStiffenerCombiTrans%i",i+1);
    SSDStiffenerCombiTrans[i] = new TGeoCombiTrans(SSDStiffenerCombiTransName,
			0.,i*SSDStiffenerSeparation,0.,new TGeoRotation("",180*(1-i),0.,0.));
  }
////////////////////////////
// Capacitor 0603-2200 nF
///////////////////////////
  const Int_t Capacitor0603Number = 5;
  TGeoBBox* Capacitor0603Shape =  new TGeoBBox("Capacitor0603Shape",
											 0.5*fgkSSDCapacitor0603Length,
											 0.5*fgkSSDCapacitor0603Width,
											 0.5*fgkSSDCapacitor0603Height);
  TGeoVolume* Capacitor0603 = new TGeoVolume("Capacitor0603",Capacitor0603Shape,
                                             fgkSSDStiffener0603CapacitorMedium); 
  Capacitor0603->SetLineColor(fColorAl);
////////////////////////////
// Capacitor 1812-330 nF
///////////////////////////
  TGeoBBox* Capacitor1812Shape =  new TGeoBBox("Capacitor1812Shape",
											 0.5*fgkSSDCapacitor1812Length,
											 0.5*fgkSSDCapacitor1812Width,
											 0.5*fgkSSDCapacitor1812Height);
  TGeoVolume* Capacitor1812 = new TGeoVolume("Capacitor1812",Capacitor1812Shape,
                                             fgkSSDStiffener1812CapacitorMedium); 
  Capacitor1812->SetLineColor(fColorAl);
  TGeoTranslation* Capacitor1812Trans = new TGeoTranslation("Capacitor1812Trans",
			  0.,
			  0.5*fgkSSDStiffenerWidth+SSDStiffenerSeparation
	       -  Capacitor1812Shape->GetDY()-fgkSSDConnectorPosition[1],
			  SSDStiffenerBoxShapes[0]->GetDZ()+fgkSSDConnectorHeight
		   +  0.5*fgkSSDCapacitor1812Height);
////////////////////////////
//Hybrid Wire
////////////////////////////
  Double_t WireX = 2.*(fgkSSDConnectorPosition[0]-0.5*fgkSSDStiffenerLength
				 - 0.5*fgkSSDConnectorLength)-fgkSSDConnectorLength
				 - fgkSSDConnectorSeparation;
  Double_t WireY = SSDStiffenerSeparation+fgkSSDStiffenerWidth
				 - 2.*fgkSSDConnectorPosition[1]-fgkSSDConnectorWidth;
  Double_t SSDWireRadius = TMath::Sqrt(TMath::Power(WireX,2.)
				         + TMath::Power(WireY,2));
  Double_t WireAngle = TMath::ATan(WireX/WireY);
  TGeoTube *HybridWireShape = new TGeoTube("HybridWireShape", 0., 
						fgkSSDWireRadius, 0.5*SSDWireRadius);
  TGeoVolume* HybridWire = new TGeoVolume("HybridWire",HybridWireShape,
                                             fgkSSDStiffenerHybridWireMedium); 
  HybridWire->SetLineColor(fColorPhynox);
  TGeoCombiTrans* HybridWireCombiTrans[2];
  HybridWireCombiTrans[0] = new TGeoCombiTrans("HybridWireCombiTrans1",
                   0.5*fgkSSDStiffenerLength-fgkSSDConnectorPosition[0]
				 + 1.5*fgkSSDConnectorLength+fgkSSDConnectorSeparation,
                   0.5*SSDWireRadius-0.5*fgkSSDStiffenerWidth
				 + fgkSSDConnectorPosition[1]+0.5*fgkSSDConnectorWidth,
				   SSDStiffenerBoxShapes[0]->GetDZ()+fgkSSDConnectorHeight
				 + fgkSSDWireRadius,
                   new TGeoRotation("HybridWireRot1",0.,90.,0.));
  HybridWireCombiTrans[1] = new TGeoCombiTrans("HybridWireCombiTrans2",
				   0.,
				 - 0.5*fgkSSDConnectorWidth+fgkSSDWireRadius,
				   0.,	
                   new TGeoRotation("HybridWireRot2",
				 -                  WireAngle*TMath::RadToDeg(),0.,0.));
  TGeoHMatrix* HybridWireMatrix = new TGeoHMatrix();
  HybridWireMatrix->MultiplyLeft(HybridWireCombiTrans[0]);
  HybridWireMatrix->MultiplyLeft(HybridWireCombiTrans[1]);
////////////////////////////
// Stiffener Assembly
///////////////////////////
  TGeoVolume* SSDStiffenerAssembly = 
								new TGeoVolumeAssembly("SSDStiffenerAssembly");
  SSDStiffenerAssembly->AddNode(HybridWire,1,HybridWireMatrix);
  for(Int_t i=0; i<SSDStiffenerNumber; i++) {
	SSDStiffenerAssembly->AddNode(SSDStiffener,i+1,SSDStiffenerCombiTrans[i]);
	for(Int_t j=1; j<Capacitor0603Number+1; j++){
    SSDStiffenerAssembly->AddNode(Capacitor0603,Capacitor0603Number*i+j,new TGeoTranslation("",(j-3.
	)/6*fgkSSDStiffenerLength,
					i*SSDStiffenerSeparation+
					0.5*((i==0? 1:-1)*fgkSSDStiffenerWidth
					+(i==0? -1:+1)*fgkSSDCapacitor0603Width),
					-0.5*(fgkSSDStiffenerHeight+fgkSSDCapacitor0603Height)));
	}
	if(i==1) SSDStiffenerAssembly->AddNode(Capacitor1812,1,Capacitor1812Trans);
}
  return SSDStiffenerAssembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDChipCables(Double_t SSDChipCablesHeigth, 
																	char* side){
////////////////////////////////////////////////////////////////////////////////
  const Int_t SSDChipCablesLayNumber = 2;
  Int_t SSDChipCablesColor[2] = {fColorAl,fColorPolyhamide};
  Double_t SSDChipCablesRadius[2];
  SSDChipCablesRadius[0] = 0.25*(SSDChipCablesHeigth
						 -  fgkSSDChipCablesHeight[0]
						 -  fgkSSDChipCablesHeight[1]);
  SSDChipCablesRadius[1] = SSDChipCablesRadius[0]-fgkSSDChipCablesHeight[0];
  Double_t SSDChipCablesPieceLength = 0.5*(fgkSSDChipCablesWidth[0]
								    - 2.*TMath::Pi()*SSDChipCablesRadius[0]
									- SSDChipCablesRadius[0]
									- fgkSSDChipCablesWidth[1]
									- fgkSSDChipCablesWidth[2]
									- (side=="Right" ? 0 : 
									  fgkSSDModuleStiffenerPosition[1]
									+ TMath::Pi()*(0.5*fgkSSDSensorHeight
									+ fgkSSDChipCablesHeight[0]
									+ fgkSSDChipCablesHeight[1])));
  //////////////////////////
  //Box and Tube Seg Shapes
  //////////////////////////
  char* SSDChipCablesBoxShapeName[2*SSDChipCablesLayNumber] = 
					  {"SSDChipCablesBoxLay0Shape0","SSDChipCablesBoxLay0Shape1",
					   "SSDChipCablesBoxLay1Shape0","SSDChipCablesBoxLay1Shape1"};
  char* SSDChipCablesTubeSegShapeName[2*SSDChipCablesLayNumber] = 
			  {"SSDChipCablesTubeSegLay0Shape0","SSDChipCablesTubeSegLay0Shape1",
			   "SSDChipCablesTubeSegLay1Shape0","SSDChipCablesTubeSegLay1Shape1"};
  TGeoBBox** SSDChipCablesBoxShape[SSDChipCablesLayNumber];
  TGeoTubeSeg** SSDChipCablesTubeSegShape[SSDChipCablesLayNumber];
  for(Int_t i=0; i<SSDChipCablesLayNumber; i++){
    SSDChipCablesBoxShape[i]        = new TGeoBBox*[2];
    SSDChipCablesTubeSegShape[i]    = new TGeoTubeSeg*[2+(side=="Right" ? 0 : 1)];
    SSDChipCablesBoxShape[i][0]     = new TGeoBBox(SSDChipCablesBoxShapeName[2*i],
												   0.5*SSDChipCablesPieceLength,
												   0.5*fgkSSDChipCablesLength[1],
												   0.5*fgkSSDChipCablesHeight[i]);
    SSDChipCablesBoxShape[i][1]     = new TGeoBBox(SSDChipCablesBoxShapeName[2*i+1],
						   0.5*(SSDChipCablesPieceLength+SSDChipCablesRadius[0]
					 + (side=="Right" ? 0. : fgkSSDModuleStiffenerPosition[1])),		
				    0.5*fgkSSDChipCablesLength[1],0.5*fgkSSDChipCablesHeight[i]);
    SSDChipCablesTubeSegShape[i][0] = 
					   new TGeoTubeSeg(SSDChipCablesTubeSegShapeName[2*i],						
						   SSDChipCablesRadius[1]-i*fgkSSDChipCablesHeight[1],
						   SSDChipCablesRadius[i],0.5*fgkSSDChipCablesLength[1],
						   0.,180.);
    SSDChipCablesTubeSegShape[i][1] = 
					   new TGeoTubeSeg(SSDChipCablesTubeSegShapeName[2*i+1],
					         SSDChipCablesRadius[0]+i*fgkSSDChipCablesHeight[0],
							 SSDChipCablesRadius[0]+fgkSSDChipCablesHeight[0]
					 +		 i*fgkSSDChipCablesHeight[1],
							 0.5*fgkSSDChipCablesLength[1],0.,180.);
    if(side!="Right") SSDChipCablesTubeSegShape[i][2] = 
					 new TGeoTubeSeg(SSDChipCablesTubeSegShapeName[2*i],
						 0.5*fgkSSDSensorHeight+(1-i)*fgkSSDChipCablesHeight[1],
						 0.5*fgkSSDSensorHeight+(1-i)*fgkSSDChipCablesHeight[0]
					 +   fgkSSDChipCablesHeight[1],
						 0.5*fgkSSDChipCablesLength[1],0.,180.);
  }
  //////////////////////////
  //Box under Chip
  //////////////////////////
  char SSDUnderChipCablesBoxShapeName[30];
  char SSDUnderChipCablesBoxName[30];	
  char SSDUnderChipCablesBoxTransName[30];	
  TGeoBBox* SSDUnderChipCablesBoxShape[SSDChipCablesLayNumber];
  TGeoVolume* SSDUnderChipCablesBox[SSDChipCablesLayNumber]; 
  TGeoTranslation* SSDUnderChipCablesBoxTrans[SSDChipCablesLayNumber];
  for(Int_t i=0; i<SSDChipCablesLayNumber; i++){ 
		sprintf(SSDUnderChipCablesBoxShapeName,"SSDUnderChipCablesBoxShape%i",i+1);
		sprintf(SSDUnderChipCablesBoxName,"SSDUnderChipCablesBox%i",i+1);
		sprintf(SSDUnderChipCablesBoxTransName,"SSDUnderChipCablesBoxTrans%i",i+1);
		SSDUnderChipCablesBoxShape[i] = 
								   new TGeoBBox(SSDUnderChipCablesBoxShapeName,
								   0.5*fgkSSDChipWidth,
								   0.5*fgkSSDChipCablesLength[1],
								   0.5*fgkSSDChipCablesHeight[i]);
		SSDUnderChipCablesBox[i] = new TGeoVolume(SSDUnderChipCablesBoxName,
									SSDUnderChipCablesBoxShape[i],
									(i==0?fgkSSDAlTraceChipCableMedium:fgkSSDKaptonChipCableMedium));		
        SSDUnderChipCablesBox[i]->SetLineColor(SSDChipCablesColor[i]);
		SSDUnderChipCablesBoxTrans[i] = 
						new TGeoTranslation(SSDUnderChipCablesBoxTransName,
						(side=="Right"?-1.:1.)*0.5*fgkSSDChipWidth,
						0.5*(fgkSSDChipCablesLength[0]-fgkSSDChipCablesLength[1])
						+0.5*fgkSSDChipCablesLength[1],
						(i==0?1.:-1.)*0.5*fgkSSDChipCablesHeight[1-i]);
  }
  //////////////////
  //Trapezoid Shapes
  //////////////////
  const Int_t SSDChipCablesVertexNumber = 2;
  const Int_t SSDChipCablesTrapezoidNumber = 2;
  TVector3** SSDChipCablesTrapezoidVertex[SSDChipCablesVertexNumber];
  for(Int_t i = 0; i< SSDChipCablesTrapezoidNumber; i++) 
	 SSDChipCablesTrapezoidVertex[i] = new TVector3*[SSDChipCablesVertexNumber];
  //First Shape Vertex Positioning
  SSDChipCablesTrapezoidVertex[0][0] = new TVector3();
  SSDChipCablesTrapezoidVertex[0][1] = 
		new TVector3(0.5*(fgkSSDChipCablesLength[0]-fgkSSDChipCablesLength[1]));
  //Second Shape Vertex Positioning
  SSDChipCablesTrapezoidVertex[1][0] = 
							  new TVector3(*SSDChipCablesTrapezoidVertex[0][0]);
  SSDChipCablesTrapezoidVertex[1][1] = 
							  new TVector3(*SSDChipCablesTrapezoidVertex[0][1]);
  //Setting the names of shapes and volumes
  char* SSDChipCablesTrapezoidBoxShapeName[SSDChipCablesTrapezoidNumber] = 
		  {"SSDChipCablesTrapezoidBoxShape1","SSDChipCablesTrapezoidBoxShape2"};
  char* SSDChipCablesTrapezoidShapeName[SSDChipCablesTrapezoidNumber] = 
		  {"SSDChipCablesTrapezoidShape1","SSDChipCablesTrapezoidShape2"};
  char* SSDChipCablesTrapezoidBoxName[SSDChipCablesTrapezoidNumber] = 
		  {"SSDChipCablesTrapezoidBox1","SSDChipCablesTrapezoidBox2"};
  char* SSDChipCablesTrapezoidName[SSDChipCablesTrapezoidNumber] = 
		  {"SSDChipCablesTrapezoid1","SSDChipCablesTrapezoid2"};
  char* SSDChipCablesTrapezoidAssemblyName[SSDChipCablesTrapezoidNumber] = 
		  {"SSDChipCablesTrapezoidAssembly1","SSDChipCablesTrapezoidAssembly2"};
  //Setting the Shapes
  TGeoBBox* SSDChipCablesTrapezoidBoxShape[SSDChipCablesTrapezoidNumber]; 
  TGeoArb8* SSDChipCablesTrapezoidShape[SSDChipCablesTrapezoidNumber];
  //Setting the Volumes
  TGeoVolume* SSDChipCablesTrapezoidBox[SSDChipCablesTrapezoidNumber];
  TGeoVolume* SSDChipCablesTrapezoid[SSDChipCablesTrapezoidNumber];
  TGeoVolume* SSDChipCablesTrapezoidAssembly[SSDChipCablesTrapezoidNumber]; 
  Double_t SSDChipCablesTrapezoidWidth[SSDChipCablesVertexNumber] = 
   {fgkSSDChipCablesWidth[1]+fgkSSDChipCablesWidth[2],fgkSSDChipCablesWidth[1]};
  for(Int_t i=0; i<SSDChipCablesTrapezoidNumber; i++){
    SSDChipCablesTrapezoidBoxShape[i] = 
					new TGeoBBox(SSDChipCablesTrapezoidBoxShapeName[i],
						0.5*fgkSSDChipCablesLength[1],
					    0.5*(fgkSSDChipCablesWidth[1]+fgkSSDChipCablesWidth[2]),
						0.5*fgkSSDChipCablesHeight[i]);
    SSDChipCablesTrapezoidShape[i] = 
							  GetTrapezoidShape(SSDChipCablesTrapezoidVertex[i],
							  SSDChipCablesTrapezoidWidth,
							  fgkSSDChipCablesHeight[i],
							  SSDChipCablesTrapezoidShapeName[i]);
    SSDChipCablesTrapezoidBox[i] = 
						new TGeoVolume(SSDChipCablesTrapezoidBoxName[i],
									   SSDChipCablesTrapezoidBoxShape[i],
									   (i==0?fgkSSDAlTraceChipCableMedium:fgkSSDKaptonChipCableMedium));
    SSDChipCablesTrapezoid[i] = new TGeoVolume(SSDChipCablesTrapezoidName[i],
											   SSDChipCablesTrapezoidShape[i],
											   (i==0?fgkSSDAlTraceChipCableMedium:fgkSSDKaptonChipCableMedium));
    SSDChipCablesTrapezoidBox[i]->SetLineColor(SSDChipCablesColor[i]);
    SSDChipCablesTrapezoid[i]->SetLineColor(SSDChipCablesColor[i]);
    SSDChipCablesTrapezoidAssembly[i] = 
				new TGeoVolumeAssembly(SSDChipCablesTrapezoidAssemblyName[i]);
    SSDChipCablesTrapezoidAssembly[i]->AddNode(SSDChipCablesTrapezoidBox[i],1,
				new TGeoTranslation(0.5*fgkSSDChipCablesLength[1],
				   0.5*(fgkSSDChipCablesWidth[1]+fgkSSDChipCablesWidth[2]),0.));
    SSDChipCablesTrapezoidAssembly[i]->AddNode(SSDChipCablesTrapezoid[i],0,
			new TGeoCombiTrans(0.,0.,0.,new TGeoRotation("",90.,180.,-90.)));
    SSDChipCablesTrapezoidAssembly[i]->AddNode(SSDChipCablesTrapezoid[i],1,
			new TGeoTranslation(fgkSSDChipCablesLength[1],0.,0.));
  }
  /////////////////////////////
  //Box and Tube Seg CombiTrans
  /////////////////////////////
  TGeoTranslation* SSDChipCablesBoxTrans[2*SSDChipCablesLayNumber];
  SSDChipCablesBoxTrans[0] = 
					new TGeoTranslation("SSDChipCablesLay1Box1Trans",0.,0.,0.);
  SSDChipCablesBoxTrans[1] = 
					new TGeoTranslation("SSDChipCablesLay1Box2Trans",
										 SSDChipCablesBoxShape[0][1]->GetDX()
						   -             0.5*SSDChipCablesPieceLength,
                       0.0,
						   -             2.*SSDChipCablesRadius[0]
						   -             fgkSSDChipCablesHeight[0]);
  SSDChipCablesBoxTrans[2] = new TGeoTranslation("SSDChipCablesLay2Box1Trans",
										 0.0,
										 0.0,
										 0.5*(fgkSSDChipCablesHeight[0]
						   +			 fgkSSDChipCablesHeight[1]));
  SSDChipCablesBoxTrans[3] = 
							 new TGeoTranslation("SSDChipCablesLay2Box2Trans",
										 SSDChipCablesBoxShape[1][1]->GetDX()
						   -			 0.5*SSDChipCablesPieceLength,
										 0.0,
						   -			 2.*SSDChipCablesRadius[0]
						   -			 0.5*fgkSSDChipCablesHeight[1]
						   -			 1.5*fgkSSDChipCablesHeight[0]);
  TGeoRotation* SSDChipCablesRot[3];
  SSDChipCablesRot[0] = new TGeoRotation("SSDChipCablesRot1",0.,90.,0.);
  SSDChipCablesRot[1] = new TGeoRotation("SSDChipCablesRot2",90.,90.,-90.);
  SSDChipCablesRot[2] = new TGeoRotation("SSDChipCablesRot3",90.,-90.,-90.);
  TGeoCombiTrans* SSDChipCablesTubeSegCombiTrans[2*(SSDChipCablesLayNumber+
													  (side=="Right" ? 0 : 1))];
  SSDChipCablesTubeSegCombiTrans[0] = 
				new TGeoCombiTrans("SSDChipCablesLay1TubeSeg1CombiTrans",
				0.5*SSDChipCablesPieceLength,
				0.0,
				SSDChipCablesRadius[0]
			-   0.5*fgkSSDChipCablesHeight[0],
				new TGeoRotation((*SSDChipCablesRot[1])*(*SSDChipCablesRot[0])));
  SSDChipCablesTubeSegCombiTrans[1] = 
				new TGeoCombiTrans("SSDChipCablesLay1TubeSeg2CombiTrans",
			-   0.5*SSDChipCablesPieceLength,
				0.0,
			-   SSDChipCablesRadius[0]-0.5*fgkSSDChipCablesHeight[0],
				new TGeoRotation((*SSDChipCablesRot[2])*(*SSDChipCablesRot[0])));
  SSDChipCablesTubeSegCombiTrans[2] = 
  new TGeoCombiTrans("SSDChipCablesLay2TubeSeg1CombiTrans",
				0.5*SSDChipCablesPieceLength,
				0.0,
				SSDChipCablesRadius[0]-0.5*fgkSSDChipCablesHeight[0],
				new TGeoRotation((*SSDChipCablesRot[1])*(*SSDChipCablesRot[0])));
  SSDChipCablesTubeSegCombiTrans[3] = 
				new TGeoCombiTrans("SSDChipCablesLay2TubeSeg2CombiTrans",
			-	0.5*SSDChipCablesPieceLength,
				0.0,
			-	SSDChipCablesRadius[0]+0.5*fgkSSDChipCablesHeight[0]
			-   fgkSSDChipCablesHeight[0],
				new TGeoRotation((*SSDChipCablesRot[2])*(*SSDChipCablesRot[0])));
  SSDChipCablesTubeSegCombiTrans[4] = 
				new TGeoCombiTrans("SSDChipCablesLay1TubeSeg4CombiTrans",
				0.5*SSDChipCablesPieceLength+SSDChipCablesRadius[0]
			+   fgkSSDModuleStiffenerPosition[1],
				0.0,
			-	2.*SSDChipCablesRadius[0]-0.5*fgkSSDChipCablesHeight[0]
			-   (0.5*fgkSSDSensorHeight+fgkSSDChipCablesHeight[0]
			+	fgkSSDChipCablesHeight[1]),
			new TGeoRotation((*SSDChipCablesRot[1])*(*SSDChipCablesRot[0])));
  SSDChipCablesTubeSegCombiTrans[5] = 
			new TGeoCombiTrans("SSDChipCablesLay2TubeSeg5CombiTrans",
				0.5*SSDChipCablesPieceLength+SSDChipCablesRadius[0]
			+	fgkSSDModuleStiffenerPosition[1],
				0.0,
			-	2.*SSDChipCablesRadius[0]-0.5*fgkSSDChipCablesHeight[0]
			-	(0.5*fgkSSDSensorHeight+fgkSSDChipCablesHeight[0]
			+	fgkSSDChipCablesHeight[1]),
			new TGeoRotation((*SSDChipCablesRot[1])*(*SSDChipCablesRot[0])));
  TGeoCombiTrans* SSDChipCablesTrapezoidCombiTrans[SSDChipCablesLayNumber];
  SSDChipCablesTrapezoidCombiTrans[0] = (side=="Right" ? 
			new TGeoCombiTrans("SSDChipCableLay1TrapezoidRightCombiTrans",
				0.5*SSDChipCablesPieceLength+SSDChipCablesTrapezoidWidth[0]
			+	SSDChipCablesRadius[0],
			-	0.5*fgkSSDChipCablesLength[1],
			-	fgkSSDChipCablesHeight[0]-2.*SSDChipCablesRadius[0],
			new TGeoRotation("",90.,0.,0.)) :
			new TGeoCombiTrans("SSDChipCableLay1TrapezoidLeftCombiTrans",
			-	2.*(fgkSSDChipCablesWidth[1]+fgkSSDChipCablesWidth[2])
			+	0.5*SSDChipCablesPieceLength+SSDChipCablesTrapezoidWidth[0]
			+	SSDChipCablesRadius[0]+fgkSSDModuleStiffenerPosition[1],
				0.5*fgkSSDChipCablesLength[1],
			-	2.*(fgkSSDChipCablesHeight[0]+fgkSSDChipCablesHeight[1])
			-	2.*SSDChipCablesRadius[0]-fgkSSDSensorHeight,
			new TGeoRotation("",-90.,0.,0.)));
  SSDChipCablesTrapezoidCombiTrans[1] = (side=="Right" ? 
			new TGeoCombiTrans("SSDChipCableLay2TrapezoidRightCombiTrans",
				0.5*SSDChipCablesPieceLength+SSDChipCablesTrapezoidWidth[0]
			+	SSDChipCablesRadius[0],
			-	0.5*fgkSSDChipCablesLength[1],
			-	0.5*(fgkSSDChipCablesHeight[0]+fgkSSDChipCablesHeight[1])
			-	fgkSSDChipCablesHeight[0]-2.*SSDChipCablesRadius[0],
				new TGeoRotation("",90.,0.,0.)) :
				new TGeoCombiTrans("SSDChipCableLay2TrapezoidLeftCombiTrans",
			-	2.*(fgkSSDChipCablesWidth[1]+fgkSSDChipCablesWidth[2])
			+	0.5*SSDChipCablesPieceLength+SSDChipCablesTrapezoidWidth[0]
			+	SSDChipCablesRadius[0]+fgkSSDModuleStiffenerPosition[1],
				0.5*fgkSSDChipCablesLength[1],-0.5*(fgkSSDChipCablesHeight[0]
			+	fgkSSDChipCablesHeight[1])-fgkSSDChipCablesHeight[1]
			-	fgkSSDChipCablesHeight[0]-2.*SSDChipCablesRadius[0]
			-	fgkSSDSensorHeight,new TGeoRotation("",-90.,0.,0.)));  
  //////////////////////////
  //Box and Tube Seg Volumes
  //////////////////////////
  char* SSDChipCablesBoxName[2*SSDChipCablesLayNumber] = 
							 {"SSDChipCablesLay1Box1","SSDChipCablesLay1Box2",
							  "SSDChipCablesLay2Box1","SSDChipCablesLay2Box2"};
  char* SSDChipRightCablesTubeSegName[2*SSDChipCablesLayNumber] = 
			  {"SSDChipRightCablesLay1TubeSeg1","SSDChipRightCablesLay1TubeSeg2",
			   "SSDChipRightCablesLay2TubeSeg1","SSDChipRightCablesLay2TubeSeg2"};
  char* SSDChipLeftCablesTubeSegName[2*SSDChipCablesLayNumber] = 
			  {"SSDChipLeftCablesLay1TubeSeg1","SSDChipLeftCablesLay1TubeSeg2",
			   "SSDChipLeftCablesLay2TubeSeg1","SSDChipLeftCablesLay2TubeSeg2"};
  char* SSDChipCablesLayAssemblyName[SSDChipCablesLayNumber] = 
			  {"SSDChipCablesLay1","SSDChipCablesLay2"};
  TGeoVolume** SSDChipCablesBox[SSDChipCablesLayNumber];
  TGeoVolume** SSDChipCablesTubeSeg[SSDChipCablesLayNumber];
  TGeoVolume* SSDChipCablesLayAssembly[SSDChipCablesLayNumber];
  for(Int_t i=0; i<SSDChipCablesLayNumber; i++){
    TGeoMedium* SSDChipCablesLayMed = 
            (i==0?fgkSSDAlTraceChipCableMedium:fgkSSDKaptonChipCableMedium);
    SSDChipCablesBox[i] = new TGeoVolume*[2];
    SSDChipCablesTubeSeg[i] = new TGeoVolume*[2+(side=="Right" ? 0 : 1)];
    SSDChipCablesBox[i][0] = new TGeoVolume(SSDChipCablesBoxName[2*i],
					   SSDChipCablesBoxShape[i][0],SSDChipCablesLayMed);
    SSDChipCablesBox[i][1] = new TGeoVolume(SSDChipCablesBoxName[2*i+1],
					   SSDChipCablesBoxShape[i][1],SSDChipCablesLayMed);
    SSDChipCablesTubeSeg[i][0] = new TGeoVolume(SSDChipRightCablesTubeSegName[2*i],
					   SSDChipCablesTubeSegShape[i][0],SSDChipCablesLayMed);
    SSDChipCablesTubeSeg[i][1] = new TGeoVolume(SSDChipRightCablesTubeSegName[2*i+1],
					   SSDChipCablesTubeSegShape[i][1],SSDChipCablesLayMed);
    SSDChipCablesBox[i][0]->SetLineColor(SSDChipCablesColor[i]);
    SSDChipCablesBox[i][1]->SetLineColor(SSDChipCablesColor[i]);
    SSDChipCablesTubeSeg[i][0]->SetLineColor(SSDChipCablesColor[i]);
    SSDChipCablesTubeSeg[i][1]->SetLineColor(SSDChipCablesColor[i]);
    SSDChipCablesLayAssembly[i] = new TGeoVolumeAssembly(SSDChipCablesLayAssemblyName[i]);
    SSDChipCablesLayAssembly[i]->AddNode(SSDChipCablesBox[i][0],1,
										 SSDChipCablesBoxTrans[2*i]);
    SSDChipCablesLayAssembly[i]->AddNode(SSDChipCablesBox[i][1],1,
										 SSDChipCablesBoxTrans[2*i+1]);
    SSDChipCablesLayAssembly[i]->AddNode(SSDChipCablesTubeSeg[i][0],1,
										 SSDChipCablesTubeSegCombiTrans[2*i]);
    SSDChipCablesLayAssembly[i]->AddNode(SSDChipCablesTubeSeg[i][1],1,
										 SSDChipCablesTubeSegCombiTrans[2*i+1]);
    if(side!="Right"){
      SSDChipCablesTubeSeg[i][2] = new TGeoVolume(SSDChipLeftCablesTubeSegName[2*i],
												  SSDChipCablesTubeSegShape[i][2],
												  SSDChipCablesLayMed);
      SSDChipCablesTubeSeg[i][2]->SetLineColor(SSDChipCablesColor[i]);
      SSDChipCablesLayAssembly[i]->AddNode(SSDChipCablesTubeSeg[i][2],1,
										   SSDChipCablesTubeSegCombiTrans[4+i]);
    }
    SSDChipCablesLayAssembly[i]->AddNode(SSDChipCablesTrapezoidAssembly[i],1,
										 SSDChipCablesTrapezoidCombiTrans[i]);
  }
  TGeoCombiTrans* SSDChipCablesCombiTrans[SSDChipCablesLayNumber];
  SSDChipCablesCombiTrans[0] = new TGeoCombiTrans("SSDChipCablesCombiTrans1",
					   (side=="Right" ? -1 : 1)*0.5*SSDChipCablesPieceLength,
						0.5*fgkSSDChipCablesLength[0],
					-	(2.*SSDChipCablesRadius[0]-fgkSSDChipCablesHeight[0]
					-	0.5*fgkSSDChipCablesHeight[1]),
						new TGeoRotation("",(side=="Right" ? 0 : 1)*180.,0.,0.));
  SSDChipCablesCombiTrans[1] = new TGeoCombiTrans("SSDChipCablesCombiTrans2",
						(side=="Right" ? -1 : 1)*0.5*SSDChipCablesPieceLength,
						0.5*fgkSSDChipCablesLength[0],
					-	(2.*SSDChipCablesRadius[0]-fgkSSDChipCablesHeight[0]
					-	0.5*fgkSSDChipCablesHeight[1]),
						new TGeoRotation("",(side=="Right" ? 0 : 1)*180.,0.,0.));
  TGeoVolume* SSDChipCablesAssembly = 
						new TGeoVolumeAssembly("SSDChipCables");
  for(Int_t i=0; i<SSDChipCablesLayNumber; i++){ 
		SSDChipCablesAssembly->AddNode(SSDChipCablesLayAssembly[i],1,
													SSDChipCablesCombiTrans[i]);
		SSDChipCablesAssembly->AddNode(SSDUnderChipCablesBox[i],1,SSDUnderChipCablesBoxTrans[i]);
  }
  return SSDChipCablesAssembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDChipCablesAssembly(Double_t SSDChipCablesHeigth){
////////////////////////////////////////////////////////////////////////////////
  const Int_t ChipCablesNumber = 2;
  Double_t ChipCablesTransVector = fgkSSDSensorLength
								 - 2.*fgkSSDModuleStiffenerPosition[1]
								 - 2.*(fgkSSDStiffenerWidth
								 - fgkSSDStiffenerToChipDist-fgkSSDChipWidth);
  char* SSDChipCablesName[ChipCablesNumber] = {"Right","Left"};
  TGeoVolume* SSDChipCables[ChipCablesNumber];  
  TGeoVolume* SSDChipCablesAssembly = 
					 new TGeoVolumeAssembly("SSDChipCablesAssembly");
  for(Int_t i=0; i<ChipCablesNumber; i++) SSDChipCables[i] = 
					 GetSSDChipCables(SSDChipCablesHeigth,SSDChipCablesName[i]);
  for(Int_t i=0; i<ChipCablesNumber; i++)
    for(Int_t j=0; j<fgkSSDChipNumber; j++)
      SSDChipCablesAssembly->AddNode(SSDChipCables[i],fgkSSDChipNumber*i+j+1,
			new TGeoTranslation(-(SSDChipCablesName[i]=="Left" ? 1. : 0.)
		*	ChipCablesTransVector,(j-0.5)*fgkSSDChipCablesLength[0]
		+	0.5*fgkSSDChipCablesLength[1],0.));
  return SSDChipCablesAssembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDFlex(Double_t SSDFlexRadius, Double_t SSDFlexHeigth){
////////////////////////////////////////////////////////////////////////////////  
  const Int_t SSDFlexVolumeNumber = 3;
  TGeoVolume* SSDFlexVolume[SSDFlexVolumeNumber];
  ////////////////////////
  // Setting Display Color
  ////////////////////////
  Int_t SSDFlexColor;
  SSDFlexColor = (SSDFlexHeigth == fgkSSDFlexHeight[0] ? fColorAl: fColorPolyhamide);
  TGeoMedium* SSDFlexMed = (SSDFlexHeigth == fgkSSDFlexHeight[0] ? fgkSSDAlTraceFlexMedium :
                            fgkSSDKaptonFlexMedium);
  ////////////////////////
  // SSDFlexTrapezoidShape
  ////////////////////////
  const Int_t SSDFlexVertexNumber = 2;
  Double_t SSDFlexWidth[SSDFlexVertexNumber] = {fgkSSDFlexWidth[1],
												fgkSSDFlexWidth[0]};
  TVector3* SSDFlexVertex[SSDFlexVertexNumber];
  SSDFlexVertex[0] = new TVector3();
  SSDFlexVertex[1] = new TVector3(fgkSSDFlexLength[0]-fgkSSDFlexLength[1]);
  TGeoArb8* SSDFlexTrapezoidShape = GetTrapezoidShape(SSDFlexVertex,
													  SSDFlexWidth,SSDFlexHeigth,
													  "SSDFlexTrapezoidShape");
  SSDFlexVolume[0] = new TGeoVolume("SSDFlexTrapezoid",SSDFlexTrapezoidShape,SSDFlexMed);
  SSDFlexVolume[0]->SetLineColor(SSDFlexColor);
  /////////////////////////
  //SSDFlexTubeSeg Assembly
  /////////////////////////
  const Int_t SSDFlexTubeSegNumber = 2;
  TGeoTubeSeg* SSDFlexTubeSegShape[SSDFlexTubeSegNumber];
  Double_t SSDFlexRadiusMax = (fgkSSDFlexLength[3]-fgkSSDFlexLength[2])
							/  TMath::Tan(fgkSSDFlexAngle*TMath::DegToRad());
  SSDFlexTubeSegShape[0] = new TGeoTubeSeg("SSDFlexTubeSegShape1",
						   SSDFlexRadius,SSDFlexRadius+SSDFlexHeigth,
						   0.5*fgkSSDFlexWidth[0],0.,180.);
  SSDFlexTubeSegShape[1] = new TGeoTubeSeg("SSDFlexTubeSegShape2",
						   SSDFlexRadiusMax-SSDFlexRadius-SSDFlexHeigth,
						   SSDFlexRadiusMax-SSDFlexRadius,0.5*fgkSSDFlexWidth[0],
						   0.,2.*fgkSSDFlexAngle);
  TGeoRotation** SSDFlexTubSegRot[SSDFlexTubeSegNumber];
  for(Int_t i = 0; i<SSDFlexTubeSegNumber; i++) 
									  SSDFlexTubSegRot[i] = new TGeoRotation*[2]; 
  SSDFlexTubSegRot[0][0] = new TGeoRotation("SSDFlexTubeSeg1Rot1", 0., 90.,  0.);
  SSDFlexTubSegRot[0][1] = new TGeoRotation("SSDFlexTubeSeg1Rot2",90., 90.,-90.);
  SSDFlexTubSegRot[1][0] = new TGeoRotation("SSDFlexTubeSeg2Rot1", 0.,-90.,  0.);
  SSDFlexTubSegRot[1][1] = new TGeoRotation("SSDFlexTubeSeg2Rot2",90., 90.,-90.);
  TGeoCombiTrans* SSDFlexTubeSegCombiTrans[SSDFlexTubeSegNumber];
  SSDFlexTubeSegCombiTrans[0] = new TGeoCombiTrans("SSDFlexTubeSegCombiTrans1",
								fgkSSDFlexLength[0],0.5*fgkSSDFlexWidth[0],
								SSDFlexRadius+0.5*SSDFlexHeigth,
								new TGeoRotation((*SSDFlexTubSegRot[0][1])
							*	(*SSDFlexTubSegRot[0][0])));
  SSDFlexTubeSegCombiTrans[1] = new TGeoCombiTrans("SSDFlexTubeSegCombiTrans2",
								fgkSSDFlexLength[0]-fgkSSDFlexLength[2],
								0.5*fgkSSDFlexWidth[0],
								SSDFlexRadiusMax+0.5*SSDFlexHeigth+SSDFlexRadius,
								new TGeoRotation((*SSDFlexTubSegRot[1][1])
							*	(*SSDFlexTubSegRot[1][0])));
  SSDFlexVolume[1] = new TGeoVolumeAssembly("SSDFlexTubeSegAssembly");
  TGeoVolume* SSDFlexTubeSeg[SSDFlexTubeSegNumber];
  char SSDFlexTubeSegName[30];
  for(Int_t i=0; i<SSDFlexTubeSegNumber; i++){ 
		sprintf(SSDFlexTubeSegName,"SSDFlexTubeSeg%i",i+1);
		SSDFlexTubeSeg[i] = new TGeoVolume(SSDFlexTubeSegName,SSDFlexTubeSegShape[i],
                                     SSDFlexMed);
		SSDFlexTubeSeg[i]->SetLineColor(SSDFlexColor);
        SSDFlexVolume[1]->AddNode(SSDFlexTubeSeg[i],1,SSDFlexTubeSegCombiTrans[i]);
  }
  ///////////
  //Box Shape 
  ///////////
  const Int_t SSDFlexBoxNumber = 7;
  Double_t SSDFlexBoxLength[SSDFlexBoxNumber];
  SSDFlexBoxLength[0] = 0.5*(fgkSSDChipNumber
					  *	fgkSSDChipLength+(fgkSSDChipNumber-1)
					  *	fgkSSDChipSeparationLength
					  - fgkSSDModuleSensorSupportDistance-fgkSSDFlexHoleLength)
					  - (fgkSSDFlexLength[0]-fgkSSDFlexLength[1]);
  SSDFlexBoxLength[1] = fgkSSDModuleSensorSupportDistance+fgkSSDFlexHoleLength;
  SSDFlexBoxLength[2] = 0.5*(fgkSSDModuleSensorSupportDistance
					  -	fgkSSDFlexHoleLength-fgkSSDFlexHoleWidth); 	
  SSDFlexBoxLength[3] = fgkSSDFlexHoleWidth;	
  SSDFlexBoxLength[4] = fgkSSDFlexLength[1]-SSDFlexBoxLength[0]
					  -	SSDFlexBoxLength[1];
  SSDFlexBoxLength[5] = fgkSSDFlexLength[2];	
  SSDFlexBoxLength[6] = fgkSSDFlexFullLength-2.*fgkSSDFlexAngle
					  *	TMath::DegToRad()*SSDFlexRadiusMax
					  - fgkSSDFlexLength[2]-TMath::Pi()
					  *	fgkSSDStiffenerHeight-fgkSSDFlexLength[0];	
  Double_t SSDFlexBoxWidth[SSDFlexBoxNumber];
  SSDFlexBoxWidth[0] = fgkSSDFlexWidth[0];
  SSDFlexBoxWidth[1] = fgkSSDFlexWidth[0]-fgkSSDFlexHoleWidth;
  SSDFlexBoxWidth[2] = fgkSSDFlexHoleWidth;
  SSDFlexBoxWidth[3] = SSDFlexBoxWidth[2]-fgkSSDFlexHoleLength;
  SSDFlexBoxWidth[4] = fgkSSDFlexWidth[0];
  SSDFlexBoxWidth[5] = fgkSSDFlexWidth[0];
  SSDFlexBoxWidth[6] = fgkSSDFlexWidth[0];
  TGeoBBox* SSDFlexBoxShape[SSDFlexBoxNumber+1];
  for(Int_t i=0; i<SSDFlexBoxNumber+1; i++) SSDFlexBoxShape[i] = 
						(i!= SSDFlexBoxNumber ? new TGeoBBox("SSDFlexBoxShape",
								0.5*SSDFlexBoxLength[i],
								0.5*SSDFlexBoxWidth[i],0.5*SSDFlexHeigth) : 
								SSDFlexBoxShape[2]);
  //////////////////////////////
  //SSDFlex Box Shape CombiTrans 
  //////////////////////////////
  TGeoCombiTrans* SSDFlexBoxCombiTrans[SSDFlexBoxNumber+1];
  SSDFlexBoxCombiTrans[0] = new TGeoCombiTrans("SSDFlexBoxCombiTrans0",
								SSDFlexVertex[1]->X()+0.5*SSDFlexBoxLength[0],
								0.5*fgkSSDFlexWidth[0],0.,0);
  SSDFlexBoxCombiTrans[1] = new TGeoCombiTrans("SSDFlexBoxCombiTrans1",
								SSDFlexVertex[1]->X()+SSDFlexBoxLength[0]
							+	0.5*SSDFlexBoxLength[1],
								fgkSSDFlexHoleWidth+0.5*SSDFlexBoxWidth[1],0.,0);
  SSDFlexBoxCombiTrans[2] = new TGeoCombiTrans("SSDFlexBoxCombiTrans2",
								SSDFlexVertex[1]->X()+SSDFlexBoxLength[0]
							+	fgkSSDFlexHoleLength+0.5*SSDFlexBoxLength[2],
								0.5*SSDFlexBoxWidth[2],0.,0);
  SSDFlexBoxCombiTrans[3] = new TGeoCombiTrans("SSDFlexBoxCombiTrans3",
								SSDFlexVertex[1]->X()+SSDFlexBoxLength[0]
							+	fgkSSDFlexHoleLength+SSDFlexBoxLength[2]
							+	0.5*fgkSSDFlexHoleWidth,
								fgkSSDFlexHoleLength+0.5*SSDFlexBoxWidth[3],0.,0);
  SSDFlexBoxCombiTrans[4] = new TGeoCombiTrans("SSDFlexBoxCombiTrans4",
								SSDFlexVertex[1]->X()+SSDFlexBoxLength[0]
							+	SSDFlexBoxLength[1]+0.5*SSDFlexBoxLength[4],
								0.5*fgkSSDFlexWidth[0],0.,0);
  SSDFlexBoxCombiTrans[5] = new TGeoCombiTrans("SSDFlexBoxCombiTrans5",
							-	0.5*fgkSSDFlexLength[2]+fgkSSDFlexLength[0],
								0.5*fgkSSDFlexWidth[0],
								2.*SSDFlexRadius+SSDFlexHeigth,0);
  SSDFlexBoxCombiTrans[6] = new TGeoCombiTrans("SSDFlexBoxCombiTrans6",
							-	SSDFlexBoxShape[6]->GetDX()
							+	SSDFlexBoxShape[6]->GetDX()
							*	TMath::Cos(2.*fgkSSDFlexAngle*TMath::DegToRad())
							+	fgkSSDFlexLength[0]-fgkSSDFlexLength[2]
							-	(SSDFlexRadiusMax-SSDFlexRadius-0.5*SSDFlexHeigth)
							*	TMath::Cos(fgkSSDFlexAngle*TMath::DegToRad()),
								0.5*fgkSSDFlexWidth[0],SSDFlexBoxShape[6]->GetDX()
								*TMath::Sin(2.*fgkSSDFlexAngle*TMath::DegToRad())
							+	SSDFlexHeigth+2.*SSDFlexRadius+(SSDFlexRadiusMax
							-	SSDFlexRadius-0.5*SSDFlexHeigth)
							*	TMath::Sin(fgkSSDFlexAngle*TMath::DegToRad()),
								new TGeoRotation("",90.,2.*fgkSSDFlexAngle,-90.));
  SSDFlexBoxCombiTrans[7] = new TGeoCombiTrans("SSDFlexBoxCombiTrans7",
								SSDFlexVertex[1]->X()+SSDFlexBoxLength[0]
							+	fgkSSDFlexHoleLength+1.5*SSDFlexBoxLength[2]
							+	SSDFlexBoxLength[3],
								0.5*SSDFlexBoxWidth[2],0.,0);
  ////////////////////////////
  //SSDFlex Box Shape Assembly 
  ////////////////////////////
  SSDFlexVolume[2] = new TGeoVolumeAssembly("SSDFlexBoxAssembly");
  TGeoVolume* SSDFlexBox[SSDFlexBoxNumber+1];
  TGeoVolume* SSDEndFlex = GetSSDEndFlex(SSDFlexBoxLength[6],SSDFlexHeigth);
  TGeoHMatrix* SSDEndFlexHMatrix = new TGeoHMatrix();
  TGeoRotation* SSDEndFlexRot= new TGeoRotation("SSDEndFlexRot",180.,0.,0);
  SSDEndFlexHMatrix->MultiplyLeft(SSDEndFlexRot);
  SSDEndFlexHMatrix->MultiplyLeft(SSDFlexBoxCombiTrans[6]);
  char SSDFlexBoxName[30];
  for(Int_t i=0; i<SSDFlexBoxNumber+1; i++){
	sprintf(SSDFlexBoxName,"SSDFlexBox%i",i!=SSDFlexBoxNumber?i+1:7);
	if(i==6){SSDFlexVolume[2]->AddNode(SSDEndFlex,1,SSDEndFlexHMatrix);}
	else{
    SSDFlexBox[i] = new TGeoVolume(SSDFlexBoxName,SSDFlexBoxShape[i],
                                   SSDFlexMed);
	SSDFlexBox[i]->SetLineColor(SSDFlexColor);
	SSDFlexVolume[2]->AddNode(SSDFlexBox[i],1,SSDFlexBoxCombiTrans[i]);}
 }
  //////////////////////
  //SSDFlex Construction
  //////////////////////
  TGeoVolume* SSDFlex = new TGeoVolumeAssembly("SSDFlex");
  for(Int_t i =0; i<SSDFlexVolumeNumber; i++) SSDFlex->AddNode(SSDFlexVolume[i],1);
  return SSDFlex;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDEndFlex(Double_t SSDEndFlexLength, 
														Double_t SSDFlexHeigth){
  /////////////////////////////////////////
  // Setting Display Color, Media and Index
  /////////////////////////////////////////
  Int_t SSDFlexColor;
  SSDFlexColor = (SSDFlexHeigth == fgkSSDFlexHeight[0] ? fColorAl: fColorPolyhamide);
  TGeoMedium* SSDFlexMed = (SSDFlexHeigth == fgkSSDFlexHeight[0] ? fgkSSDAlTraceFlexMedium :
                            fgkSSDKaptonFlexMedium);
  ////////////////////////
  const Int_t SSDEndFlexBoxNumber = 5;
  TGeoBBox* SSDEndFlexBBoxShape[SSDEndFlexBoxNumber];
  SSDEndFlexBBoxShape[0] = new TGeoBBox("SSDFlexBoxShape1",
								   0.5*SSDEndFlexLength,0.5*fgkSSDFlexWidth[0],
								   0.5*SSDFlexHeigth);
  SSDEndFlexBBoxShape[1] = new TGeoBBox("SSDFlexBoxShape2",
				    0.5*fgkSSDEndFlexCompLength[1],
					0.5*(fgkSSDEndFlexCompWidth[0]-fgkSSDFlexWidth[0])/2,
					0.5*SSDFlexHeigth);
  SSDEndFlexBBoxShape[2] = new TGeoBBox("SSDFlexBoxShape3",
				    0.5*fgkSSDEndFlexCompLength[2],
					0.5*(fgkSSDEndFlexCompWidth[1]-fgkSSDFlexWidth[0])/2,
					0.5*SSDFlexHeigth);
  SSDEndFlexBBoxShape[3] = new TGeoBBox("SSDFlexBoxShape4",
				    0.5*fgkSSDEndFlexCompLength[3],
					0.5*(fgkSSDEndFlexCompWidth[0]-fgkSSDFlexWidth[0])/2,
					0.5*SSDFlexHeigth);
  SSDEndFlexBBoxShape[4] = new TGeoBBox("SSDFlexBoxShape5",
				    0.5*(fgkSSDEndFlexCompLength[4]+fgkSSDEndFlexCompLength[5]),
					0.25*(fgkSSDEndFlexCompWidth[2]-fgkSSDFlexWidth[0])/2,
					0.5*SSDFlexHeigth);
  TGeoVolume* SSDEndFlexBBox[SSDEndFlexBoxNumber];  
  char SSDEndFlexBBoxName[30];
  for(Int_t i=0; i<SSDEndFlexBoxNumber; i++){
	sprintf(SSDEndFlexBBoxName,"SSDEndFlexBBox%i",i+1);
	SSDEndFlexBBox[i] = new TGeoVolume(SSDEndFlexBBoxName,
                     SSDEndFlexBBoxShape[i],
                     SSDFlexMed);
	SSDEndFlexBBox[i]->SetLineColor(SSDFlexColor);
  }
  TGeoVolume* SSDEndFlex = new TGeoVolumeAssembly("SSDEndFlex");
  Double_t PartialSumLength = 0.;
  for(Int_t i=0; i<SSDEndFlexBoxNumber+1; i++) PartialSumLength += fgkSSDEndFlexCompLength[i];
  Double_t ReferenceLength = SSDEndFlexLength-PartialSumLength;
  SSDEndFlex->AddNode(SSDEndFlexBBox[0],1);
  SSDEndFlex->AddNode(SSDEndFlexBBox[1],1,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+ReferenceLength+fgkSSDEndFlexCompLength[0]
					+  0.5*fgkSSDEndFlexCompLength[1],
					   0.5*fgkSSDFlexWidth[0]+SSDEndFlexBBoxShape[1]->GetDY(),
					   0.));
  SSDEndFlex->AddNode(SSDEndFlexBBox[1],2,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+ReferenceLength+fgkSSDEndFlexCompLength[0]
					+  0.5*fgkSSDEndFlexCompLength[1],
					-  0.5*fgkSSDFlexWidth[0]-SSDEndFlexBBoxShape[1]->GetDY(),
					   0.));
  SSDEndFlex->AddNode(SSDEndFlexBBox[2],1,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+ReferenceLength+fgkSSDEndFlexCompLength[0]
					+  fgkSSDEndFlexCompLength[1]+0.5*fgkSSDEndFlexCompLength[2],
					+  0.5*fgkSSDFlexWidth[0]+SSDEndFlexBBoxShape[2]->GetDY(),
					   0.));
  SSDEndFlex->AddNode(SSDEndFlexBBox[2],2,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+ReferenceLength+fgkSSDEndFlexCompLength[0]
					+  fgkSSDEndFlexCompLength[1]+0.5*fgkSSDEndFlexCompLength[2],
					-  0.5*fgkSSDFlexWidth[0]-SSDEndFlexBBoxShape[2]->GetDY(),
					   0.));
  SSDEndFlex->AddNode(SSDEndFlexBBox[3],1,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+ReferenceLength+fgkSSDEndFlexCompLength[0]
					+  fgkSSDEndFlexCompLength[1]+fgkSSDEndFlexCompLength[2]
					+  0.5*fgkSSDEndFlexCompLength[3],
					+  0.5*fgkSSDFlexWidth[0]+SSDEndFlexBBoxShape[3]->GetDY(),
					   0.));
  SSDEndFlex->AddNode(SSDEndFlexBBox[3],2,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+ReferenceLength+fgkSSDEndFlexCompLength[0]
					+  fgkSSDEndFlexCompLength[1]+fgkSSDEndFlexCompLength[2]
					+  0.5*fgkSSDEndFlexCompLength[3],
					-  0.5*fgkSSDFlexWidth[0]-SSDEndFlexBBoxShape[3]->GetDY(),
					   0.));
  SSDEndFlex->AddNode(SSDEndFlexBBox[4],1,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+ReferenceLength+fgkSSDEndFlexCompLength[0]
					+  fgkSSDEndFlexCompLength[1]+fgkSSDEndFlexCompLength[2]
					+  fgkSSDEndFlexCompLength[3]+0.5*(fgkSSDEndFlexCompLength[4]+fgkSSDEndFlexCompLength[5]),
					+  0.5*fgkSSDFlexWidth[0]+SSDEndFlexBBoxShape[4]->GetDY(),
					   0.));
  SSDEndFlex->AddNode(SSDEndFlexBBox[4],2,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+ReferenceLength+fgkSSDEndFlexCompLength[0]
					+  fgkSSDEndFlexCompLength[1]+fgkSSDEndFlexCompLength[2]
					+  fgkSSDEndFlexCompLength[3]+0.5*(fgkSSDEndFlexCompLength[4]
					+  fgkSSDEndFlexCompLength[5]),
					-  0.5*fgkSSDFlexWidth[0]-SSDEndFlexBBoxShape[4]->GetDY(),
					   0.));
  return SSDEndFlex;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDFlexAssembly(){
////////////////////////////////////////////////////////////////////////////////
  TGeoVolume* SSDFlexAssembly = new TGeoVolumeAssembly("SSDFlexAssembly");
  const Int_t SSDFlexLayerNumber = 4;
  Double_t SSDFlexHeight[SSDFlexLayerNumber];
  Double_t SSDFlexRadius[SSDFlexLayerNumber];
  TGeoTranslation* SSDFlexTrans[SSDFlexLayerNumber];
  for(Int_t i=0; i<SSDFlexLayerNumber; i++){ 
    SSDFlexHeight[i] = (i%2==0 ? fgkSSDFlexHeight[0] : fgkSSDFlexHeight[1]);
    SSDFlexRadius[i] = (i==0 ? fgkSSDStiffenerHeight : SSDFlexRadius[i-1]
					 +					   SSDFlexHeight[i-1]);
    SSDFlexTrans[i]  = new TGeoTranslation(0.,0.,-0.5*i*(SSDFlexHeight[0]
					 +					   SSDFlexHeight[1])); 
    SSDFlexAssembly->AddNode(GetSSDFlex(SSDFlexRadius[i],SSDFlexHeight[i]),i+1,
										   SSDFlexTrans[i]);   
  }
  return SSDFlexAssembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDCoolingBlockAssembly(){
////////////////////////////////////////////////////////////////////////////////
  const Int_t SSDCoolingBlockTransNumber = 2;
  Double_t SSDCoolingBlockTransVector[SSDCoolingBlockTransNumber] = 
					{fgkSSDModuleSensorSupportDistance+fgkSSDCoolingBlockLength,
					 fgkSSDSensorLength-2.*fgkSSDModuleStiffenerPosition[1]
				-	 fgkSSDCoolingBlockWidth};
  TGeoVolume* SSDCoolingBlock = GetSSDCoolingBlock();
  TGeoVolume* SSDCoolingBlockAssembly = 
							  new TGeoVolumeAssembly("SSDCoolingBlockAssembly");
  for(Int_t i=0; i<SSDCoolingBlockTransNumber; i++)
    for(Int_t j=0; j<SSDCoolingBlockTransNumber; j++) 
		SSDCoolingBlockAssembly->AddNode(SSDCoolingBlock,
						  SSDCoolingBlockTransNumber*i+j+1,
						  new TGeoTranslation(i*SSDCoolingBlockTransVector[0],
						  j*SSDCoolingBlockTransVector[1],0.));
  return SSDCoolingBlockAssembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDCoolingBlock(){
/////////////////////////////////////////////////////////////////////////////////
  // Center Cooling Block Hole
  ////////////////////////////
  Double_t CoolingBlockHoleAngle = TMath::ACos(0.5*fgkSSDCoolingBlockHoleLength[0]
							/fgkSSDCoolingBlockHoleRadius[0])*TMath::RadToDeg();
  Double_t CoolingBlockHoleWidth = fgkSSDCoolingBlockWidth;
  /*TGeoTubeSeg* CoolingBlockHoleShape = */ new TGeoTubeSeg("CoolingBlockHoleShape",
						0.,
						fgkSSDCoolingBlockHoleRadius[0],
						0.5*CoolingBlockHoleWidth,
						180.-CoolingBlockHoleAngle,360.+CoolingBlockHoleAngle);
  TVector3* CoolingBlockHoleVertex[3];
  CoolingBlockHoleVertex[0] = new TVector3();
  CoolingBlockHoleVertex[1] = new TVector3(fgkSSDCoolingBlockHoleRadius[0]
					*	TMath::Cos((90.-CoolingBlockHoleAngle)*TMath::DegToRad()),
						fgkSSDCoolingBlockHoleRadius[0]
					*	TMath::Sin((90.-CoolingBlockHoleAngle)*TMath::DegToRad()));
  CoolingBlockHoleVertex[2] = new TVector3(CoolingBlockHoleVertex[1]->X(),
					-	CoolingBlockHoleVertex[1]->Y());
  /* TGeoArb8* CoolingBlockTriangleHoleShape = */
						GetTriangleShape(CoolingBlockHoleVertex,
						CoolingBlockHoleWidth,"CoolingBlockTriangleHoleShape");
  TGeoRotation* CoolingBlockHoleRot = 
							  new TGeoRotation("CoolingBlockHoleRot",90,0.,0.);
  CoolingBlockHoleRot->RegisterYourself();
  /* TGeoCompositeShape* CoolingTubeHoleShape  = */ 
                              new TGeoCompositeShape("CoolingTubeHoleShape",
							  "CoolingBlockTriangleHoleShape:CoolingBlockHoleRot+"
							  "CoolingBlockHoleShape");
  ///////////////////////////
  // Cooling Block Trapezoids
  ///////////////////////////
  const Int_t VertexNumber = 4;
  const Int_t TrapezoidNumber = 2;
  TVector3** CoolingBlockTrapezoidVertex[TrapezoidNumber];
  for(Int_t i = 0; i<TrapezoidNumber; i++) CoolingBlockTrapezoidVertex[i] = 
						new TVector3*[VertexNumber]; 
  Double_t CoolingBlockComponentHeight    = fgkSSDCoolingBlockHeight[0]
				    -	fgkSSDCoolingBlockHoleCenter
					-	fgkSSDCoolingBlockHoleRadius[0]
					*	TMath::Sin(CoolingBlockHoleAngle*TMath::DegToRad());
  Double_t CoolingBlockTrapezoidLength[TrapezoidNumber] = 
					{	fgkSSDCoolingBlockLength,
						0.5*(fgkSSDCoolingBlockLength-2.
					*	(fgkSSDCoolingBlockHoleLength[1]
					-	fgkSSDCoolingBlockHoleRadius[1])
					-	fgkSSDCoolingBlockHoleLength[0])}; 
  Double_t CoolingBlockTrapezoidHeigth[TrapezoidNumber] = 
					{	fgkSSDCoolingBlockHeight[0]-CoolingBlockComponentHeight
					-	fgkSSDCoolingBlockHeight[1]-fgkSSDCoolingBlockHeight[2]
					-	fgkSSDCoolingBlockHoleRadius[1],
						CoolingBlockComponentHeight};
  Double_t CoolingBlockTrapezoidWidth[TrapezoidNumber]  = 
						{fgkSSDCoolingBlockWidth,fgkSSDCoolingBlockWidth};
  //////////////////////////
  //Vertex Positioning Shape 
  //////////////////////////
  CoolingBlockTrapezoidVertex[0][0] = new TVector3();
  CoolingBlockTrapezoidVertex[0][1] = new TVector3(CoolingBlockTrapezoidLength[0]);
  CoolingBlockTrapezoidVertex[0][2] = new TVector3(
						0.5*(CoolingBlockTrapezoidVertex[0][1]->X()
					-	2.*CoolingBlockTrapezoidLength[1]
					-	fgkSSDCoolingBlockHoleLength[0]));
  CoolingBlockTrapezoidVertex[0][3] = 
						new TVector3(CoolingBlockTrapezoidVertex[0][1]->X()
					-	CoolingBlockTrapezoidVertex[0][2]->X());
  CoolingBlockTrapezoidVertex[1][0] = new TVector3(); 
  CoolingBlockTrapezoidVertex[1][1] = new TVector3(CoolingBlockTrapezoidLength[1]);
  CoolingBlockTrapezoidVertex[1][2] = 
						new TVector3(CoolingBlockTrapezoidHeigth[1]
					/				 CoolingBlockTrapezoidHeigth[0]
					*	CoolingBlockTrapezoidVertex[0][2]->X());
  CoolingBlockTrapezoidVertex[1][3] = 
						new TVector3(CoolingBlockTrapezoidVertex[1][1]->X());
  char* CoolingBlockTrapezoidShapeName[TrapezoidNumber] = 
					{"CoolingBlockTrapezoidShape0","CoolingBlockTrapezoidShape1"};
  TGeoArb8* CoolingBlockTrapezoidShape[TrapezoidNumber];
  for(Int_t i = 0; i< TrapezoidNumber; i++) CoolingBlockTrapezoidShape[i] = 
						GetArbShape(CoolingBlockTrapezoidVertex[i],
						CoolingBlockTrapezoidWidth,
						CoolingBlockTrapezoidHeigth[i],
						CoolingBlockTrapezoidShapeName[i]);
  TGeoTranslation* CoolingBlockTrapezoidTrans = 
						new TGeoTranslation("CoolingBlockTrapezoidTrans",
						CoolingBlockTrapezoidVertex[0][2]->X(),
						0.0,
						0.5*(CoolingBlockTrapezoidHeigth[0]
					+	CoolingBlockTrapezoidHeigth[1]));
  CoolingBlockTrapezoidTrans->RegisterYourself();
  TGeoCombiTrans* CoolingBlockTrapezoidCombiTrans = 
						new TGeoCombiTrans("CoolingBlockTrapezoidCombiTrans",
						CoolingBlockTrapezoidVertex[0][3]->X(),
						fgkSSDCoolingBlockWidth,
						0.5*(CoolingBlockTrapezoidHeigth[0]
					+	CoolingBlockTrapezoidHeigth[1]),
						new TGeoRotation("",180.,0.,0.));
  CoolingBlockTrapezoidCombiTrans->RegisterYourself();
  /* TGeoCompositeShape* CoolingBlockTrapezoidCompositeShape = */ 
	new TGeoCompositeShape("CoolingBlockTrapezoidCompositeShape",
	"CoolingBlockTrapezoidShape0+CoolingBlockTrapezoidShape1:CoolingBlockTrapezoidTrans+"
	"CoolingBlockTrapezoidShape1:CoolingBlockTrapezoidCombiTrans"); 
  /////////////////////////////
  // Cooling Block Boxes Shapes
  /////////////////////////////
  const Int_t BoxNumber = 3;
  TGeoBBox* CoolingBlockBoxShape[BoxNumber];
  CoolingBlockBoxShape[0] = new TGeoBBox("CoolingBlockBoxShape0",
						0.5*fgkSSDCoolingBlockLength,
						0.5*fgkSSDCoolingBlockWidth,
						0.5*fgkSSDCoolingBlockHoleRadius[1]);
  CoolingBlockBoxShape[1] = new TGeoBBox("CoolingBlockBoxShape1",
						0.5*(fgkSSDCoolingBlockLength
					-	2.*fgkSSDCoolingBlockHoleLength[1]),
						0.5*fgkSSDCoolingBlockWidth,
						0.5*fgkSSDCoolingBlockHeight[2]);
  CoolingBlockBoxShape[2] = new TGeoBBox("CoolingBlockBoxShape2",
						0.5*fgkSSDCoolingBlockLength,
						0.5*fgkSSDCoolingBlockWidth,
						0.5*fgkSSDCoolingBlockHeight[1]);
  TGeoTranslation* CoolingBlockBoxTrans[BoxNumber-1];
  CoolingBlockBoxTrans[0] = new TGeoTranslation("CoolingBlockBoxTrans0",0.,0.,
						0.5*(fgkSSDCoolingBlockHeight[1]
					+	fgkSSDCoolingBlockHoleRadius[1])
					+	fgkSSDCoolingBlockHeight[2]);
  CoolingBlockBoxTrans[1] = new TGeoTranslation("CoolingBlockBoxTrans1",
						0.0,
						0.0,
						0.5*(fgkSSDCoolingBlockHeight[1]
					+	fgkSSDCoolingBlockHeight[2]));
  for(Int_t i=0; i<BoxNumber-1; i++) CoolingBlockBoxTrans[i]->RegisterYourself();
  /* TGeoCompositeShape* CoolingBlockBoxCompositeShape = */
	new TGeoCompositeShape("CoolingBlockBoxCompositeShape",
						   "CoolingBlockBoxShape0:CoolingBlockBoxTrans0+"
	 "CoolingBlockBoxShape1:CoolingBlockBoxTrans1+CoolingBlockBoxShape2");
  ///////////////////////
  // Cooling Block Shape
  //////////////////////
  TGeoCombiTrans* CoolingTubeHoleShapeCombiTrans = 
						new TGeoCombiTrans("CoolingTubeHoleShapeCombiTrans",
						0.5*fgkSSDCoolingBlockLength,
						0.5*fgkSSDCoolingBlockWidth,
						fgkSSDCoolingBlockHoleCenter,
						new TGeoRotation("",0.,90.,0.));
  CoolingTubeHoleShapeCombiTrans->RegisterYourself();
  TGeoTranslation* CoolingBlockTrapezoidCompositeShapeTrans = 
						new TGeoTranslation("CoolingBlockTrapezoidCompositeShapeTrans",
						0.0,
						0.0,
						0.5*CoolingBlockTrapezoidHeigth[0]+fgkSSDCoolingBlockHeight[1]+
						fgkSSDCoolingBlockHeight[2]+fgkSSDCoolingBlockHoleRadius[1]);
  CoolingBlockTrapezoidCompositeShapeTrans->RegisterYourself();
  TGeoTranslation* CoolingBlockBoxCompositeShapeTrans = 
						new TGeoTranslation("CoolingBlockBoxCompositeShapeTrans",
						0.5*fgkSSDCoolingBlockLength,
						0.5*fgkSSDCoolingBlockWidth,
						0.5*fgkSSDCoolingBlockHeight[1]);
  CoolingBlockBoxCompositeShapeTrans->RegisterYourself();
  TGeoCompositeShape* SSDCoolingBlockShape = 
		new TGeoCompositeShape("SSDCoolingBlockShape",	
		"CoolingBlockBoxCompositeShape:CoolingBlockBoxCompositeShapeTrans+"
		"CoolingBlockTrapezoidCompositeShape:CoolingBlockTrapezoidCompositeShapeTrans-"
		"CoolingTubeHoleShape:CoolingTubeHoleShapeCombiTrans");
  TGeoVolume* SSDCoolingBlock = new TGeoVolume("SSDCoolingBlock",
		SSDCoolingBlockShape,fgkSSDAlCoolBlockMedium);
  return SSDCoolingBlock;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCarbonFiberJunction(Double_t width){
////////////////////////////////////////////////////////////////////////////////
  const Int_t VertexNumber = 4;
  TVector3* Vertex[VertexNumber];
  Vertex[0] = new TVector3();
  Vertex[1] = new TVector3(fgkCarbonFiberJunctionLength,0.);
  Vertex[2] = new TVector3(fgkCarbonFiberJunctionLength
	        -	  fgkCarbonFiberJunctionEdge[1]
			*	  TMath::Cos(fgkCarbonFiberJunctionAngle[1]*TMath::DegToRad()),
				  fgkCarbonFiberJunctionEdge[1]*TMath::Sin(fgkCarbonFiberJunctionAngle[1]
			*     TMath::DegToRad()));
  Vertex[3] = new TVector3(fgkCarbonFiberJunctionEdge[0]
			*	  TMath::Cos(fgkCarbonFiberJunctionAngle[0]*TMath::DegToRad()),
				  fgkCarbonFiberJunctionEdge[0]
			*	  TMath::Sin(fgkCarbonFiberJunctionAngle[0]*TMath::DegToRad()));
  TGeoArb8* CarbonFiberJunctionShapePiece = 
						new TGeoArb8("CarbonFiberJunctionShapePiece",0.5*width);
  //////////////////////////////////
  //Setting the vertices in TGeoArb8
  //////////////////////////////////
  for(Int_t i = 0; i<2*VertexNumber; i++)
	CarbonFiberJunctionShapePiece->SetVertex(i,
							Vertex[(i < VertexNumber ? i: i-VertexNumber)]->X(),
							Vertex[(i < VertexNumber ? i : i-VertexNumber)]->Y());
  TGeoRotation* CarbonFiberJunctionRot = 
						new TGeoRotation("CarbonFiberJunctionRot",
										  180.,
										  180.,
										  180-2.*fgkCarbonFiberJunctionAngle[0]); 
  TGeoVolume* CarbonFiberJunctionPiece = 
						new TGeoVolume("CarbonFiberJunctionPiece",
						CarbonFiberJunctionShapePiece,fgkSSDCarbonFiberMedium);
  TGeoVolume* CarbonFiberJunction = 
						new TGeoVolumeAssembly("CarbonFiberJunction");
  CarbonFiberJunctionPiece->SetLineColor(fColorCarbonFiber);
  CarbonFiberJunction->AddNode(CarbonFiberJunctionPiece,1);
  CarbonFiberJunction->AddNode(CarbonFiberJunctionPiece,2,CarbonFiberJunctionRot);
  return CarbonFiberJunction;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCarbonFiberJunctionAssembly(){
////////////////////////////////////////////////////////////////////////////////
  SetCarbonFiberJunctionCombiTransMatrix();
  TGeoVolume* CarbonFiberJunctionAssembly = 
						  new TGeoVolumeAssembly("CarbonFiberJunctionAssembly");
  TGeoVolume* CarbonFiberJunction = 
						  GetCarbonFiberJunction(fgkCarbonFiberJunctionWidth);
  for(Int_t i=0; i<fgkCarbonFiberJunctionCombiTransNumber;i++) 
    CarbonFiberJunctionAssembly->AddNode(CarbonFiberJunction,i+1,
										 CarbonFiberJunctionCombiTransMatrix[i]);
  return CarbonFiberJunctionAssembly;
}
/////////////////////////////////////////////////////////////////////////////////
TList* AliITSv11GeometrySSD::GetLadderCableSegment(Double_t SSDEndLadderCableLength){
/////////////////////////////////////////////////////////////////////////////////
const Int_t LadderCableSegmentNumber = 2;
/////////////////////////////////////////
// LadderSegmentBBox Volume
/////////////////////////////////////////
  TGeoBBox* LadderCableSegmentBBoxShape[LadderCableSegmentNumber];
  const char* LadderCableSegmentBBoxShapeName[LadderCableSegmentNumber] = 
				{"LadderCableSegmentBBoxShape1","LadderCableSegmentBBoxShape2"};
  for(Int_t i=0; i<LadderCableSegmentNumber; i++) LadderCableSegmentBBoxShape[i] = 
						  new TGeoBBox(LadderCableSegmentBBoxShapeName[i],
									   0.5*fgkSSDFlexWidth[0],
									   0.5*fgkSSDLadderCableWidth,
									   0.5*fgkSSDFlexHeight[i]); 
  const char* LadderCableSegmentBBoxName[LadderCableSegmentNumber] = 
						  {"LadderCableSegmentBBox1","LadderCableSegmentBBox2"};
  TGeoVolume* LadderCableSegmentBBox[LadderCableSegmentNumber];
  for(Int_t i=0; i<LadderCableSegmentNumber; i++){ 
			LadderCableSegmentBBox[i] =
						  new TGeoVolume(LadderCableSegmentBBoxName[i],
										 LadderCableSegmentBBoxShape[i],
										 (i==0?fgkSSDAlTraceLadderCableMedium:
            fgkSSDKaptonLadderCableMedium));
			LadderCableSegmentBBox[i]->SetLineColor(i==0 ? fColorAl : 
														   fColorPolyhamide);
  }
  TGeoTranslation* LadderCableSegmentBBoxTrans[LadderCableSegmentNumber];										  
  LadderCableSegmentBBoxTrans[0] = 
						   new TGeoTranslation("LadderCableSegmentBBoxTrans1",
											   0.5*fgkSSDFlexWidth[0],
											   0.5*fgkSSDLadderCableWidth,
											   0.5*fgkSSDFlexHeight[0]);
  LadderCableSegmentBBoxTrans[1] = 
						   new TGeoTranslation("LadderCableSegmentBBoxTrans2",
											   0.5*fgkSSDFlexWidth[0],
											   0.5*fgkSSDLadderCableWidth,
											   fgkSSDFlexHeight[0]
											   +0.5*fgkSSDFlexHeight[1]);
  TGeoVolume* LadderCableSegmentBBoxAssembly = 
						   new TGeoVolumeAssembly("LadderCableSegmentBBoxAssembly"); 
  for(Int_t i=0; i<LadderCableSegmentNumber; i++)  
		LadderCableSegmentBBoxAssembly->AddNode(LadderCableSegmentBBox[i],1,
											    LadderCableSegmentBBoxTrans[i]);
/////////////////////////////////////////
// LadderSegmentArb8 Volume
/////////////////////////////////////////
  const Int_t VertexNumber = 4;
  TVector3** LadderCableSegmentVertexPosition[LadderCableSegmentNumber];
  for(Int_t i = 0; i<LadderCableSegmentNumber; i++) LadderCableSegmentVertexPosition[i] = 
												  new TVector3*[VertexNumber];
//Shape Vertex Positioning
  for(Int_t i=0; i<LadderCableSegmentNumber; i++){
	LadderCableSegmentVertexPosition[i][0] = new TVector3(0.,i*fgkSSDFlexHeight[0]);
	LadderCableSegmentVertexPosition[i][1] = new TVector3(fgkSSDLadderCableWidth,
														  i*fgkSSDFlexHeight[0]);
	LadderCableSegmentVertexPosition[i][2] = new TVector3(0.,fgkSSDFlexHeight[0]
										   +			     fgkSSDFlexHeight[1]
										   +			  i*fgkSSDFlexHeight[0]);
	LadderCableSegmentVertexPosition[i][3] = 
						   new TVector3(LadderCableSegmentVertexPosition[i][1]->X(),
										LadderCableSegmentVertexPosition[i][2]->Y());
  }
  Double_t LadderCableSegmentWidth[2][2] = {{fgkSSDFlexHeight[0],fgkSSDFlexHeight[0]},
							     		    {fgkSSDFlexHeight[1],fgkSSDFlexHeight[1]}};	
  char* LadderCableSegmentArbShapeName[LadderCableSegmentNumber] = 
					{"LadderCableSegmentArbShape1","LadderCableSegmentArbShape2"};
  TGeoArb8* LadderCableSegmentArbShape[LadderCableSegmentNumber];
  for(Int_t i = 0; i< LadderCableSegmentNumber; i++) LadderCableSegmentArbShape[i] = 
					GetArbShape(LadderCableSegmentVertexPosition[i],
								LadderCableSegmentWidth[i],
								fgkCarbonFiberJunctionWidth-fgkSSDFlexWidth[0],
								LadderCableSegmentArbShapeName[i]);
  const char* LadderCableSegmentArbName[LadderCableSegmentNumber] = 
						  {"LadderCableSegmentArb1","LadderCableSegmentArb2"};
  TGeoVolume* LadderCableSegmentArb[LadderCableSegmentNumber];
  for(Int_t i=0; i<LadderCableSegmentNumber; i++){
			 LadderCableSegmentArb[i] =
						   new TGeoVolume(LadderCableSegmentArbName[i],
										  LadderCableSegmentArbShape[i],
										  (i==0?fgkSSDAlTraceLadderCableMedium:
            fgkSSDKaptonLadderCableMedium)); 
			 LadderCableSegmentArb[i]->SetLineColor(i==0 ? fColorAl : 
														   fColorPolyhamide);
}
  TGeoRotation* LadderCableSegmentArbRot[LadderCableSegmentNumber];
  LadderCableSegmentArbRot[0] = new TGeoRotation("LadderCableSegmentArbRot1",
												 90.,90,-90.);	 
  LadderCableSegmentArbRot[1] = new TGeoRotation("LadderCableSegmentArbRot2",
												  0.,90.,0.);	 
  TGeoCombiTrans* LadderCableSegmentArbCombiTrans =  
						   new TGeoCombiTrans("LadderCableSegmentArbCombiTrans",
							   0.5*(fgkCarbonFiberJunctionWidth-fgkSSDFlexWidth[0])
							 + fgkSSDFlexWidth[0],0.,0.,
						   new TGeoRotation((*LadderCableSegmentArbRot[1])
						     *(*LadderCableSegmentArbRot[0])));
  TGeoVolume* LadderCableSegmentArbAssembly = 
						   new TGeoVolumeAssembly("LadderCableSegmentArbAssembly"); 
  for(Int_t i=0; i<LadderCableSegmentNumber; i++)
  LadderCableSegmentArbAssembly->AddNode(LadderCableSegmentArb[i],1,
										   LadderCableSegmentArbCombiTrans);
/////////////////////////////////////////
// End Ladder Cable Volume
/////////////////////////////////////////
  TGeoBBox* LadderEndCableSegmentBBoxShape[LadderCableSegmentNumber];
  const char* LadderEndCableSegmentBBoxShapeName[LadderCableSegmentNumber] = 
				{"LadderEndCableSegmentBBoxShape1","LadderEndCableSegmentBBoxShape2"};
  for(Int_t i=0; i<LadderCableSegmentNumber; i++) LadderEndCableSegmentBBoxShape[i] = 
						  new TGeoBBox(LadderEndCableSegmentBBoxShapeName[i],
									   0.5*SSDEndLadderCableLength,
									   0.5*fgkSSDLadderCableWidth,
									   0.5*fgkSSDFlexHeight[i]);
  const char* LadderEndCableSegmentBBoxName[LadderCableSegmentNumber] = 
						  {"LadderEndCableSegmentBBox1","LadderEndCableSegmentBBox2"};
  TGeoVolume* LadderEndCableSegmentBBox[LadderCableSegmentNumber];
  for(Int_t i=0; i<LadderCableSegmentNumber; i++){ 
			LadderEndCableSegmentBBox[i] =
						  new TGeoVolume(LadderEndCableSegmentBBoxName[i],
										 LadderEndCableSegmentBBoxShape[i],
										 (i==0?fgkSSDAlTraceLadderCableMedium:
            fgkSSDKaptonLadderCableMedium));
			LadderEndCableSegmentBBox[i]->SetLineColor(i==0 ? fColorAl : 
														   fColorPolyhamide);
  }
  TGeoTranslation* LadderEndCableSegmentBBoxTrans[LadderCableSegmentNumber];										  
  LadderEndCableSegmentBBoxTrans[0] = 
						   new TGeoTranslation("LadderEndCableSegmentBBoxTrans0",
											   0.5*SSDEndLadderCableLength,
											   0.5*fgkSSDLadderCableWidth,
											   0.5*fgkSSDFlexHeight[0]);
  LadderEndCableSegmentBBoxTrans[1] = 
						   new TGeoTranslation("LadderEndCableSegmentBBoxTrans1",
											   0.5*SSDEndLadderCableLength,
											   0.5*fgkSSDLadderCableWidth,
											   fgkSSDFlexHeight[0]
											   +0.5*fgkSSDFlexHeight[1]);
  TGeoVolume* LadderEndCableSegmentBBoxAssembly = 
						   new TGeoVolumeAssembly("LadderEndCableSegmentBBoxAssembly"); 
  for(Int_t i=0; i<LadderCableSegmentNumber; i++)  
		LadderEndCableSegmentBBoxAssembly->AddNode(LadderEndCableSegmentBBox[i],1,
											    LadderEndCableSegmentBBoxTrans[i]);
/////////////////////////////////////////
  TList* LadderCableSegmentList = new TList();
  LadderCableSegmentList->Add(LadderCableSegmentBBoxAssembly);
  LadderCableSegmentList->Add(LadderCableSegmentArbAssembly);
  LadderCableSegmentList->Add(LadderEndCableSegmentBBoxAssembly);
  return LadderCableSegmentList;
  }
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetLadderCable(Int_t n, Double_t SSDEndLadderCableLength){
/////////////////////////////////////////////////////////////////////////////////
  TList* LadderCableSegmentList = GetLadderCableSegment(SSDEndLadderCableLength);
  TGeoVolume* LadderCable = new TGeoVolumeAssembly("LadderCable"); 
  for(Int_t i=0; i<n; i++){
	TGeoTranslation* LadderCableTrans = new TGeoTranslation(
							i*(fgkCarbonFiberJunctionWidth),
							fgkSSDLadderCableWidth-fgkSSDFlexWidth[0],
							i*(fgkSSDFlexHeight[0]+fgkSSDFlexHeight[1]));
    LadderCable->AddNode((TGeoVolume*)LadderCableSegmentList->At(0),i+1,LadderCableTrans);  
	if(i!=n-1) LadderCable->AddNode((TGeoVolume*)LadderCableSegmentList->At(1),i+1,LadderCableTrans);  
  }
  TGeoTranslation* EndLadderCableTrans = new TGeoTranslation("EndLadderCableTrans",
					  (n-1)*fgkCarbonFiberJunctionWidth+fgkSSDFlexWidth[0],
								 fgkSSDLadderCableWidth-fgkSSDFlexWidth[0],
					  (n-1)*(fgkSSDFlexHeight[0]+fgkSSDFlexHeight[1]));
  LadderCable->AddNode((TGeoVolume*)LadderCableSegmentList->At(2),1,EndLadderCableTrans);
  return LadderCable;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetLadderCableAssembly(Int_t n, Double_t SSDEndLadderCableLength){
/////////////////////////////////////////////////////////////////////////////////
  TGeoVolume* LadderCableAssembly = new TGeoVolumeAssembly("LadderCableAssembly");
  char LadderCableTransName[30];
  for(Int_t i=0; i<n; i++){ 
	sprintf(LadderCableTransName,"LadderCableTrans%i",i+1);
    LadderCableAssembly->AddNode(GetLadderCable(n-i,SSDEndLadderCableLength),i+1,
	new TGeoTranslation(LadderCableTransName,i*fgkCarbonFiberJunctionWidth,0.,0.));
  }
  return LadderCableAssembly;
}
/////////////////////////////////////////////////////////////////////////////////
TList* AliITSv11GeometrySSD::GetLadderCableAssemblyList(Int_t n, Double_t SSDEndLadderCableLength){
/////////////////////////////////////////////////////////////////////////////////
  const Int_t LadderCableAssemblyNumber = 2;
  TGeoVolume* LadderCableAssembly = GetLadderCableAssembly(n,SSDEndLadderCableLength);
  TGeoVolume* LadderCable[LadderCableAssemblyNumber];
  char LadderCableAssemblyName[30];
  TList* LadderCableAssemblyList = new TList();
  for(Int_t i=0; i<LadderCableAssemblyNumber; i++){ 
	sprintf(LadderCableAssemblyName,"LadderCableAssembly%i",i+1);
	LadderCable[i] = new TGeoVolumeAssembly(LadderCableAssemblyName);
	LadderCable[i]->AddNode(LadderCableAssembly,i+1,i==0 ? NULL :
					 new TGeoCombiTrans((n-1)
					 *	 fgkCarbonFiberJunctionWidth+fgkSSDFlexWidth[0],
					     2.*fgkSSDLadderCableWidth+0.5*fgkSSDFlexWidth[0],
											0.,new TGeoRotation("",180,0.,0.)));
	LadderCableAssemblyList->Add(LadderCable[i]);
}
  return LadderCableAssemblyList;
}
/////////////////////////////////////////////////////////////////////////////////
TList* AliITSv11GeometrySSD::GetEndLadderCarbonFiberJunctionAssembly(){
////////////////////////////////////////////////////////////////////////////////
  const Int_t EndLabberCarbonFiberJunctionNumber = 2;
  TGeoVolume* EndLadderCarbonFiberJunctionAssembly[EndLabberCarbonFiberJunctionNumber];
  EndLadderCarbonFiberJunctionAssembly[0] = 
				new TGeoVolumeAssembly("EndLadderCarbonFiberJunctionAssembly1");
  EndLadderCarbonFiberJunctionAssembly[1] = 
				new TGeoVolumeAssembly("EndLadderCarbonFiberJunctionAssembly2");
  TGeoVolume** EndLadderCarbonFiberJunction[EndLabberCarbonFiberJunctionNumber];
  for(Int_t i=0; i<EndLabberCarbonFiberJunctionNumber; i++) 
						   EndLadderCarbonFiberJunction[i] = new TGeoVolume*[2];
  for(Int_t i=0; i<EndLabberCarbonFiberJunctionNumber; i++){
    EndLadderCarbonFiberJunction[i][0] = 
		  GetCarbonFiberJunction(fgkEndLadderCarbonFiberLowerJunctionLength[i]);
    EndLadderCarbonFiberJunction[i][1] = 
		  GetCarbonFiberJunction(fgkEndLadderCarbonFiberUpperJunctionLength[i]);
  }
  TList* EndLadderCarbonFiberJunctionList = new TList();
  for(Int_t i=0; i<EndLabberCarbonFiberJunctionNumber; i++){
    SetEndLadderCarbonFiberJunctionCombiTransMatrix(i);
    for(Int_t j=0; j<fgkCarbonFiberJunctionCombiTransNumber; j++)
      EndLadderCarbonFiberJunctionAssembly[i]->AddNode(j==2 ? 
						 EndLadderCarbonFiberJunction[i][1] : 
						 EndLadderCarbonFiberJunction[i][0],
						 j+1,EndLadderCarbonFiberJunctionCombiTransMatrix[j]);
    EndLadderCarbonFiberJunctionList->Add(EndLadderCarbonFiberJunctionAssembly[i]);
  }
  return EndLadderCarbonFiberJunctionList;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCarbonFiberSupport(){
////////////////////////////////////////////////////////////////////////////////
  const Int_t VertexNumber = 4;
  const Int_t ShapesNumber = 2;
  TVector3** VertexPosition[ShapesNumber];
  for(Int_t i=0; i<ShapesNumber; i++) VertexPosition[i] = new TVector3*[VertexNumber];
  Double_t CarbonFiberSupportXAxisEdgeProj = 
		fgkCarbonFiberSupportEdgeLength*TMath::Cos(fgkCarbonFiberJunctionAngle[0]
	*	TMath::DegToRad());
  Double_t Theta = TMath::ATan(fgkCarbonFiberSupportYAxisLength
				 /			   fgkCarbonFiberSupportXAxisLength);
  /////////////////////
  //Vertex Positioning
  ////////////////////
  VertexPosition[0][0] = new TVector3();
  VertexPosition[0][1] = new TVector3(fgkCarbonFiberSupportXAxisLength,
									  fgkCarbonFiberSupportYAxisLength);
  VertexPosition[0][2] = new TVector3(CarbonFiberSupportXAxisEdgeProj,
									  CarbonFiberSupportXAxisEdgeProj
					   *			  TMath::Tan(Theta));
  VertexPosition[0][3] = new TVector3(fgkCarbonFiberSupportXAxisLength
					   -			  CarbonFiberSupportXAxisEdgeProj,
									  fgkCarbonFiberSupportYAxisLength
					   -			  VertexPosition[0][2]->Y());
  ////////////////////////////////////////////////////
  //Setting the parameters for Isometry Transformation
  ////////////////////////////////////////////////////
  Double_t SymmetryPlanePosition = (fgkCarbonFiberSupportYAxisLength
								 +	fgkCarbonFiberSupportTopEdgeDist[0]
								 +	fgkCarbonFiberSupportWidth);
  Double_t* param = new Double_t[4]; 
  param[0] = 0., param[1] = 1., param[2] = 0., param[3] = -SymmetryPlanePosition;
  for(Int_t j=0; j<VertexNumber; j++) VertexPosition[1][j] = 
				  new TVector3((GetReflection(VertexPosition[0][j],param))->X(),
							  (GetReflection(VertexPosition[0][j],param))->Y());
  char* CarbonFiberSupportShapeName[ShapesNumber] = 
						{"CarbonFiberSupportShape1","CarbonFiberSupportShape2"};
  TGeoArb8* CarbonFiberSupportShape[ShapesNumber]; 
  Double_t Width[2] = {fgkCarbonFiberSupportWidth,fgkCarbonFiberSupportWidth};
  Double_t CarbonFiberSupportHeight = 
	  CarbonFiberSupportXAxisEdgeProj*TMath::Tan(fgkCarbonFiberJunctionAngle[0]
	  *TMath::DegToRad());
  for(Int_t i = 0; i< ShapesNumber; i++) CarbonFiberSupportShape[i] = 
					GetArbShape(VertexPosition[i],Width,CarbonFiberSupportHeight,
									CarbonFiberSupportShapeName[i],i==0 ? 1: -1);
  /////////////////////////////////////
  //Setting Translations and Rotations: 
  /////////////////////////////////////
  TGeoTranslation* CarbonFiberSupportTrans = 
						  new TGeoTranslation("CarbonFiberSupportTrans",
											  0.0,0.0,0.5*CarbonFiberSupportHeight);
  CarbonFiberSupportTrans->RegisterYourself();
  TGeoRotation* CarbonFiberCompShapeRot[2];
  CarbonFiberCompShapeRot[0] = new TGeoRotation("CarbonFiberCompShapeRot1",
											  0.0,180.0,0.0);
  CarbonFiberCompShapeRot[1] = new TGeoRotation("CarbonFiberCompShapeRot2",
										  90.,-fgkCarbonFiberTriangleAngle,-90.);
  Double_t TransVector[3] = {fgkCarbonFiberTriangleLength
						  *  TMath::Cos(fgkCarbonFiberTriangleAngle
						  *	 TMath::DegToRad()),0.,-fgkCarbonFiberTriangleLength
						  *	 TMath::Sin(fgkCarbonFiberTriangleAngle
						  *	 TMath::DegToRad())};
  TGeoCombiTrans* CarbonFiberSupportCombiTrans = 
							   new TGeoCombiTrans("CarbonFiberSupportCombiTrans",
							   TransVector[0],2.*SymmetryPlanePosition
						  +	   TransVector[1],TransVector[2],
							   new TGeoRotation((*CarbonFiberCompShapeRot[1])
						  *	   (*CarbonFiberCompShapeRot[0])));
  CarbonFiberSupportCombiTrans->RegisterYourself();
////////////////////////////////////////////////////////////////////////////////
  TGeoCompositeShape* CarbonFiberSupportCompShape = 
							new TGeoCompositeShape("CarbonFiberSupportCompShape",
							"CarbonFiberSupportShape1:CarbonFiberSupportTrans+"
							"CarbonFiberSupportShape2:CarbonFiberSupportTrans");
  TGeoVolume* CarbonFiberSupport = new TGeoVolume("CarbonFiberSupport",
						   CarbonFiberSupportCompShape,fgkSSDCarbonFiberMedium);
  CarbonFiberSupport->SetLineColor(fColorCarbonFiber);
  TGeoVolume* CarbonFiberSupportAssembly = 
						    new TGeoVolumeAssembly("CarbonFiberSupportAssembly");
  CarbonFiberSupportAssembly->AddNode(CarbonFiberSupport,1);
  CarbonFiberSupportAssembly->AddNode(CarbonFiberSupport,2,
									  CarbonFiberSupportCombiTrans);
  return CarbonFiberSupportAssembly;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCarbonFiberLowerSupport(Int_t ikind, Bool_t EndLadder){
////////////////////////////////////////////////////////////////////////////////
  const Int_t VertexNumber = 4;
  const Int_t ShapesNumber = 2;
  Double_t Width[2] = {fgkCarbonFiberLowerSupportWidth,
								fgkCarbonFiberLowerSupportWidth};
  TVector3** VertexPosition[ShapesNumber];
  for(Int_t i = 0; i<ShapesNumber; i++) VertexPosition[i] = 
						 new TVector3*[VertexNumber];
  //First Shape Vertex Positioning
  VertexPosition[0][0] = new TVector3(fgkCarbonFiberLowerSupportLowerLenght);
  VertexPosition[0][1] = new TVector3(fgkCarbonFiberTriangleLength
					   -		fgkCarbonFiberLowerSupportLowerLenght);
  VertexPosition[0][2] = new TVector3();
  VertexPosition[0][3] = new TVector3(fgkCarbonFiberTriangleLength);
  //Second Shape Vertex Positioning
  Double_t Theta = TMath::ATan((fgkCarbonFiberLowerSupportVolumePosition[1]
				 -				fgkCarbonFiberLowerSupportVolumePosition[0])
				 /				fgkCarbonFiberTriangleLength);
  VertexPosition[1][0] = new TVector3(VertexPosition[0][0]->X(),
								VertexPosition[0][0]->X()*TMath::Tan(Theta)
				 +				fgkCarbonFiberLowerSupportVolumePosition[0]);
  VertexPosition[1][1] = new TVector3(VertexPosition[0][1]->X(),
								VertexPosition[0][1]->X()*TMath::Tan(Theta)
				 +				fgkCarbonFiberLowerSupportVolumePosition[0]);
  VertexPosition[1][2] = new TVector3(0.,fgkCarbonFiberLowerSupportVolumePosition[0]);
  VertexPosition[1][3] = new TVector3(fgkCarbonFiberTriangleLength,
								fgkCarbonFiberLowerSupportVolumePosition[1]);
  char* CarbonFiberLowerSupportName[ShapesNumber] = 
			  {"CarbonFiberLowerSupportShape1","CarbonFiberLowerSupportShape2"};
  TGeoArb8* CarbonFiberLowerSupportShape[ShapesNumber];
  for(Int_t i = 0; i< ShapesNumber; i++) CarbonFiberLowerSupportShape[i] = 
								GetArbShape(VertexPosition[i],Width,
											fgkCarbonFiberLowerSupportHeight,
											CarbonFiberLowerSupportName[i]);
  ///////////////////////////////////////////////////////
  TGeoTranslation* CarbonFiberLowerSupportTrans[ShapesNumber];
  CarbonFiberLowerSupportTrans[0] = 
						new TGeoTranslation("CarbonFiberLowerSupportTrans1",
						0.0,
						VertexPosition[1][3]->Y()+VertexPosition[1][2]->Y(),
						0.0);
  CarbonFiberLowerSupportTrans[1] = 
						new TGeoTranslation("CarbonFiberLowerSupportTrans2",
						0.0,
				-		VertexPosition[1][3]->Y()-VertexPosition[1][2]->Y(),
						0.0);
  for(Int_t i = 0; i< ShapesNumber; i++) 
						CarbonFiberLowerSupportTrans[i]->RegisterYourself(); 
  ///////////////////////////////////////////////////////
  TGeoCompositeShape* CarbonFiberLowerSupportCompShape; 
  if(EndLadder==false)
    CarbonFiberLowerSupportCompShape = 
				new TGeoCompositeShape("CarbonFiberLowerSupportCompShape",
				"CarbonFiberLowerSupportShape2+"
				"CarbonFiberLowerSupportShape1:CarbonFiberLowerSupportTrans1");
  else
    if(ikind==0)
      CarbonFiberLowerSupportCompShape = 
						  (TGeoCompositeShape*)CarbonFiberLowerSupportShape[0];
    else
      CarbonFiberLowerSupportCompShape = 
	new TGeoCompositeShape("CarbonFiberLowerSupportCompShape",
				 "CarbonFiberLowerSupportShape1+"
				 "CarbonFiberLowerSupportShape1:CarbonFiberLowerSupportTrans1"); 
  TGeoVolume* CarbonFiberLowerSupport = new TGeoVolume("CarbonFiberLowerSupport",
					  CarbonFiberLowerSupportCompShape,fgkSSDCarbonFiberMedium);
  CarbonFiberLowerSupport->SetLineColor(fColorCarbonFiber);
  return CarbonFiberLowerSupport;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCarbonFiberAssemblySupport(){
////////////////////////////////////////////////////////////////////////////////
  SetCarbonFiberAssemblyCombiTransMatrix();
  TGeoVolume* CarbonFiberAssemblySupport = 
						new TGeoVolumeAssembly("CarbonFiberAssembly");
  TGeoVolume* CarbonFiberAssemblyVolumes[fgkCarbonFiberAssemblyCombiTransNumber];
  CarbonFiberAssemblyVolumes[0] = GetCarbonFiberJunctionAssembly();
  CarbonFiberAssemblyVolumes[1] = GetCarbonFiberSupport();
  CarbonFiberAssemblyVolumes[2] = GetCarbonFiberLowerSupport();
  for(Int_t i=0; i<fgkCarbonFiberAssemblyCombiTransNumber;i++) 
    CarbonFiberAssemblySupport->AddNode(CarbonFiberAssemblyVolumes[i],1,
						CarbonFiberAssemblyCombiTransMatrix[i]);
  return CarbonFiberAssemblySupport;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCoolingTubeSupport(){
////////////////////////////////////////////////////////////////////////////////
  const Int_t VertexNumber = 3;
  Double_t Phi = TMath::ASin(0.5*fgkCoolingTubeSupportHeight
			   /			 fgkCoolingTubeSupportRmax)*TMath::RadToDeg();
  /* TGeoTubeSeg* CoolingTubeSegShape =  */ 
						new TGeoTubeSeg("CoolingTubeSegShape",0.0,
							fgkCoolingTubeSupportRmax,
							0.5*fgkCoolingTubeSupportWidth,Phi,
							360-Phi);
  /* TGeoTube* CoolingTubeHoleShape =  */
						new TGeoTube("CoolingTubeHoleShape",0.0,
							fgkCoolingTubeSupportRmin,
							0.5*fgkCoolingTubeSupportWidth);
  TVector3* VertexPosition[VertexNumber];
  ///////////////////////////
  //Shape Vertex Positioning
  ///////////////////////////
  VertexPosition[0] = new TVector3();
  VertexPosition[1] = new TVector3(fgkCoolingTubeSupportRmax
					*		TMath::Cos(Phi*TMath::DegToRad()),
							fgkCoolingTubeSupportRmax
					*		TMath::Sin(Phi*TMath::DegToRad()));
  VertexPosition[2] = new TVector3(VertexPosition[1]->X(),
					-			   VertexPosition[1]->Y());
  /* TGeoArb8* CoolingTubeTriangleShape = */ GetTriangleShape(VertexPosition,
									   fgkCoolingTubeSupportWidth,
									   "CoolingTubeTriangleShape");
  Double_t* BoxOrigin = new Double_t[3];
  Double_t BoxLength = fgkCoolingTubeSupportLength-fgkCoolingTubeSupportRmax
					 - VertexPosition[1]->X();
  BoxOrigin[0] = VertexPosition[1]->X()+0.5*BoxLength, BoxOrigin[1] = BoxOrigin[2] = 0.;
  /* TGeoBBox* CoolingTubeBoxShape = */
						new TGeoBBox("CoolingTubeBoxShape",0.5*BoxLength,
							0.5*fgkCoolingTubeSupportHeight,
							0.5*fgkCoolingTubeSupportWidth,BoxOrigin);
  TGeoCompositeShape* CoolingTubeSupportShape = 
						new TGeoCompositeShape("CoolingTubeSupportShape",
						 "(CoolingTubeSegShape+CoolingTubeTriangleShape"
						 "+CoolingTubeBoxShape)-CoolingTubeHoleShape"); 
  TGeoVolume* CoolingTubeSupport = new TGeoVolume("CoolingTubeSupport",
						 CoolingTubeSupportShape,fgkSSDTubeHolderMedium);
  return CoolingTubeSupport;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCoolingTubeSupportAssembly(){
////////////////////////////////////////////////////////////////////////////////
  TGeoVolume* CoolingTubeSupportAssembly = 
						   new TGeoVolumeAssembly("CoolingTubeSupportAssembly");
  TGeoVolume* CoolingTubeSupport = GetCoolingTubeSupport();
  SetCoolingTubeSupportCombiTransMatrix();
  for(Int_t i=0; i<fgkCoolingTubeSupportCombiTransNumber;i++) 
    CoolingTubeSupportAssembly->AddNode(CoolingTubeSupport,i+1,
										 CoolingTubeSupportCombiTransMatrix[i]);
  return CoolingTubeSupportAssembly;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCoolingTube(){
////////////////////////////////////////////////////////////////////////////////
  TGeoVolume* CoolingTubeAssembly = new TGeoVolumeAssembly("CoolingTubeAssembly");
  TGeoTube *CoolingTubeShape = new TGeoTube("CoolingTubeShape", fgkCoolingTubeRmin, 
						fgkCoolingTubeRmax, fgkCoolingTubeLength/2.0);
  TGeoVolume* CoolingTube = new TGeoVolume("CoolingTube",
						 CoolingTubeShape,fgkSSDCoolingTubePhynox);
  TGeoTube *CoolingTubeInteriorShape = new TGeoTube("CoolingTubeInteriorShape", 
						0, fgkCoolingTubeRmin, 
						fgkCoolingTubeLength/2.0);
  TGeoVolume *CoolingTubeInterior = new TGeoVolume("CoolingTubeInterior",
						   CoolingTubeInteriorShape,fgkSSDCoolingTubeWater);
  CoolingTubeAssembly->AddNode(CoolingTube,1);
  CoolingTubeAssembly->AddNode(CoolingTubeInterior,2);
  return CoolingTubeAssembly;
}  
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCoolingTubeAssembly(){
////////////////////////////////////////////////////////////////////////////////
  TGeoVolume* CoolingTubeAssembly =   new TGeoVolumeAssembly("CoolingTubeAssembly");
  TGeoVolume* CoolingTube = GetCoolingTube();
  SetCoolingTubeCombiTransMatrix();
  for(Int_t i=0; i<fgkCoolingTubeCombiTransNumber;i++) 
    CoolingTubeAssembly->AddNode(CoolingTube,i+1,CoolingTubeTransMatrix[i]);
  return CoolingTubeAssembly;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetLadderSegment(Int_t iChipCablesHeight){
////////////////////////////////////////////////////////////////////////////////
    TGeoVolume* LadderSegment = new TGeoVolumeAssembly("LadderSegment");
    TGeoVolume* LadderSegmentVolumes[fgkLadderSegmentCombiTransNumber];
    LadderSegmentVolumes[0] = GetCarbonFiberAssemblySupport();
    LadderSegmentVolumes[1] = GetSSDModule(iChipCablesHeight);
    LadderSegmentVolumes[2] = GetSSDSensorSupportAssembly(iChipCablesHeight);
    LadderSegmentVolumes[3] = GetCoolingTubeSupportAssembly();
	LadderSegmentVolumes[4] = GetCoolingTubeAssembly(); 
    SetLadderSegmentCombiTransMatrix();
    for(Int_t i=0; i<fgkLadderSegmentCombiTransNumber; i++) 
      LadderSegment->AddNode(LadderSegmentVolumes[i],1,
							 LadderSegmentCombiTransMatrix[i]);
    return LadderSegment;
}
////////////////////////////////////////////////////////////////////////////////
TList* AliITSv11GeometrySSD::GetEndLadderSegment(){
////////////////////////////////////////////////////////////////////////////////
  const Int_t EndLadderSegmentNumber = 2;
  TList* EndLadderCarbonFiberJunctionList = GetEndLadderCarbonFiberJunctionAssembly();
  TGeoVolume* EndLadderSegment[EndLadderSegmentNumber];
  EndLadderSegment[0] = new TGeoVolumeAssembly("EndLadderSegment1");
  EndLadderSegment[1] = new TGeoVolumeAssembly("EndLadderSegment2");
  TGeoVolume** LadderSegmentVolumes[EndLadderSegmentNumber];
  const Int_t LadderSegmentVolumeNumber = 4;
  for(Int_t i=0; i<EndLadderSegmentNumber; i++) LadderSegmentVolumes[i] = 
						new TGeoVolume*[LadderSegmentVolumeNumber];
  LadderSegmentVolumes[0][0] = (TGeoVolume*)EndLadderCarbonFiberJunctionList->At(0);
  LadderSegmentVolumes[0][1] = GetCarbonFiberSupport();
  LadderSegmentVolumes[0][2] = GetSSDMountingBlock();
  LadderSegmentVolumes[0][3] = GetCarbonFiberLowerSupport(0,true);
  LadderSegmentVolumes[1][0] = (TGeoVolume*)EndLadderCarbonFiberJunctionList->At(1);
  LadderSegmentVolumes[1][1] = LadderSegmentVolumes[0][1];
  LadderSegmentVolumes[1][2] = LadderSegmentVolumes[0][2];
  LadderSegmentVolumes[1][3] = GetCarbonFiberLowerSupport(1,true);
  TList* EndLadderSegmentList = new TList();
  for(Int_t i=0; i<EndLadderSegmentNumber; i++){
    SetEndLadderSegmentCombiTransMatrix(i);
    for(Int_t j=0; j<LadderSegmentVolumeNumber; j++)
      EndLadderSegment[i]->AddNode(LadderSegmentVolumes[i][j],1,
								   EndLadderSegmentCombiTransMatrix[j]);
    EndLadderSegmentList->Add(EndLadderSegment[i]);
  }
  return EndLadderSegmentList;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDMountingBlock(){
////////////////////////////////////////////////////////////////////////////////
  // Mounting Block Boxes Shapes
  ///////////////////////////////////////
  const Int_t MountingBlockBoxNumber = 3;
  TGeoBBox* MountingBlockBoxShape[MountingBlockBoxNumber];
  MountingBlockBoxShape[0] = new TGeoBBox("MountingBlockBoxShape0",
							0.25*(fgkSSDMountingBlockLength[0]
						-	fgkSSDMountingBlockLength[1]),
							0.5*fgkSSDMountingBlockWidth,
							0.5*fgkSSDMountingBlockHeight[0]);
  MountingBlockBoxShape[1] = new TGeoBBox("MountingBlockBoxShape1",
							0.25*(fgkSSDMountingBlockLength[1]
						-	fgkSSDMountingBlockLength[2]),
							0.5*fgkSSDMountingBlockWidth,
							0.5*(fgkSSDMountingBlockHeight[1]
						-	fgkSSDMountingBlockHeight[3]));
  MountingBlockBoxShape[2] = new TGeoBBox("MountingBlockBoxShape2",
							0.5*fgkSSDMountingBlockLength[2],
							0.5*fgkSSDMountingBlockWidth,
							0.5*(fgkSSDMountingBlockHeight[2]
						-	fgkSSDMountingBlockHeight[3]));
  TGeoTranslation* MountingBlockBoxTrans[MountingBlockBoxNumber+2];
  MountingBlockBoxTrans[0] = new TGeoTranslation("MountingBlockBoxTrans0",0.,0.,0.);
  MountingBlockBoxTrans[1] = new TGeoTranslation("MountingBlockBoxTrans1",
							MountingBlockBoxShape[0]->GetDX()
						+	MountingBlockBoxShape[1]->GetDX(),
							0.0,
							MountingBlockBoxShape[1]->GetDZ()
						-	MountingBlockBoxShape[0]->GetDZ()
						+	fgkSSDMountingBlockHeight[3]);
  MountingBlockBoxTrans[2] = new TGeoTranslation("MountingBlockBoxTrans2",
							MountingBlockBoxShape[0]->GetDX()
						+	2.*MountingBlockBoxShape[1]->GetDX()
						+	MountingBlockBoxShape[2]->GetDX(),
							0.0,
							MountingBlockBoxShape[2]->GetDZ()
						-	MountingBlockBoxShape[0]->GetDZ()
						+	fgkSSDMountingBlockHeight[3]);
  MountingBlockBoxTrans[3] = new TGeoTranslation("MountingBlockBoxTrans3",
							MountingBlockBoxShape[0]->GetDX()
						+	MountingBlockBoxShape[1]->GetDX()
						+	2.*(MountingBlockBoxShape[1]->GetDX()
						+	MountingBlockBoxShape[2]->GetDX()),
							0.0,
							MountingBlockBoxShape[1]->GetDZ()
						-	MountingBlockBoxShape[0]->GetDZ()
						+	fgkSSDMountingBlockHeight[3]);
  MountingBlockBoxTrans[4] = new TGeoTranslation("MountingBlockBoxTrans4",
							2.*(MountingBlockBoxShape[0]->GetDX()
						+	2.*MountingBlockBoxShape[1]->GetDX()
						+   MountingBlockBoxShape[2]->GetDX()),
							0.0,
							0.0);
  for(Int_t i=0; i<MountingBlockBoxNumber+2; i++) 
									MountingBlockBoxTrans[i]->RegisterYourself();
  ///////////////////////////////////////
  // Mounting Block Trapezoid Hole Shapes
  ///////////////////////////////////////
  const Int_t HoleTrapezoidVertexNumber = 4;
  TVector3* HoleTrapezoidVertex[HoleTrapezoidVertexNumber];
  HoleTrapezoidVertex[0] = new TVector3();
  HoleTrapezoidVertex[1] = new TVector3(fgkSSDMountingBlockHoleTrapezoidHeight);
  HoleTrapezoidVertex[2] = new TVector3(*HoleTrapezoidVertex[0]);
  HoleTrapezoidVertex[3] = new TVector3(*HoleTrapezoidVertex[1]);
  Double_t HoleTrapezoidWidth[2] = {fgkSSDMountingBlockHoleTrapezoidUpBasis
						+	2.*MountingBlockBoxShape[1]->GetDX()
						*	TMath::Tan(fgkSSDMountingBlockHoleTrapezoidAngle
						*	TMath::DegToRad()),
							fgkSSDMountingBlockHoleTrapezoidUpBasis}; 
  /* TGeoArb8* HoleTrapezoidShape =*/ GetArbShape(HoleTrapezoidVertex,
							HoleTrapezoidWidth,
							2.*MountingBlockBoxShape[1]->GetDX(),
							"HoleTrapezoidShape");
  TGeoRotation* HoleTrapezoidShapeRot[2];
  HoleTrapezoidShapeRot[0] = new TGeoRotation("HoleTrapezoidShapeRot0",
							90.,-90.,-90.);
  HoleTrapezoidShapeRot[1] = new TGeoRotation("HoleTrapezoidShapeRot1",
							-180.,0.,0.);
  TGeoCombiTrans* HoleTrapezoidShapeCombiTrans = 
							new TGeoCombiTrans("HoleTrapezoidShapeCombiTrans",
							MountingBlockBoxShape[0]->GetDX()
						+	3.*MountingBlockBoxShape[1]->GetDX()
						+	2.*MountingBlockBoxShape[2]->GetDX(),
							0.5*fgkSSDMountingBlockWidth,
						-	fgkSSDMountingBlockHoleTrapezoidHeight
						+	2.*MountingBlockBoxShape[1]->GetDZ()
						-	MountingBlockBoxShape[0]->GetDZ()
						+	fgkSSDMountingBlockHeight[3],
							new TGeoRotation((*HoleTrapezoidShapeRot[1])
						*	(*HoleTrapezoidShapeRot[0])));
  HoleTrapezoidShapeCombiTrans->RegisterYourself();
  ///////////////////////////////////
  // Mounting Block Screw Hole Shapes
  ///////////////////////////////////
  const Int_t MountingBlockTubeNumber = 2;
  TGeoTube* MountingBlockTubeShape[MountingBlockTubeNumber];
  MountingBlockTubeShape[0] = new TGeoTube("MountingBlockTubeShape0",0.0,
							fgkSSDMountingBlockHoleRadius,
							MountingBlockBoxShape[0]->GetDZ());
  MountingBlockTubeShape[1] = new TGeoTube("MountingBlockTubeShape1",0.0,
							fgkSSDMountingBlockHoleRadius,
							MountingBlockBoxShape[2]->GetDZ());
  TGeoTranslation* MountingBlockTubeTrans[2*MountingBlockTubeNumber];
  MountingBlockTubeTrans[0] = new TGeoTranslation("MountingBlockTubeTrans0",
						-	0.5*(fgkSSDMountingBlockLength[0]
						-	fgkSSDMountingBlockHoleTubeLength[0]),
							0.5*fgkSSDMountingBlockWidth
						-	fgkSSDMountingBlockHoleTubeWidth[0],0.);
  MountingBlockTubeTrans[1] = new TGeoTranslation("MountingBlockTubeTrans1",
						-	0.5*(fgkSSDMountingBlockLength[0]
						-	fgkSSDMountingBlockHoleTubeLength[0])
						+	fgkSSDMountingBlockHoleTubeLength[0],
						-	0.5*fgkSSDMountingBlockWidth
						+	fgkSSDMountingBlockHoleTubeWidth[0],
							0.);
  MountingBlockTubeTrans[2] = new TGeoTranslation("MountingBlockTubeTrans2",
						-	MountingBlockBoxShape[0]->GetDX()
						+	0.5*fgkSSDMountingBlockLength[0]
						-	fgkSSDMountingBlockHoleTubeLength[1],
							0.5*fgkSSDMountingBlockWidth
						-	fgkSSDMountingBlockHoleTubeWidth[0],
							MountingBlockBoxShape[2]->GetDZ()
						-	MountingBlockBoxShape[0]->GetDZ()
						+	fgkSSDMountingBlockHeight[3]);
  MountingBlockTubeTrans[3] = new TGeoTranslation("MountingBlockTubeTrans3",
						-	MountingBlockBoxShape[0]->GetDX()
						+	0.5*fgkSSDMountingBlockLength[0],
						-	0.5*fgkSSDMountingBlockWidth
						+	fgkSSDMountingBlockHoleTubeWidth[1],
							MountingBlockBoxShape[2]->GetDZ()
						-	MountingBlockBoxShape[0]->GetDZ()
						+	fgkSSDMountingBlockHeight[3]);
  for(Int_t i=0; i<2*MountingBlockTubeNumber; i++) 
						MountingBlockTubeTrans[i]->RegisterYourself();
  /* TGeoCompositeShape* MountingBlockMainShape = */
						new TGeoCompositeShape("MountingBlockMainShape",
						"MountingBlockBoxShape0:MountingBlockBoxTrans0+"
						"MountingBlockBoxShape1:MountingBlockBoxTrans1+"
						"MountingBlockBoxShape2:MountingBlockBoxTrans2+"
						"MountingBlockBoxShape1:MountingBlockBoxTrans3+"
						"MountingBlockBoxShape0:MountingBlockBoxTrans4");
  ////////////////////////////////////////////
  // Mounting Block Screw Composite Hole Shapes
  ////////////////////////////////////////////
  const Int_t MountingBlockHoleTubeSegNumber = 4;
  /* TGeoTubeSeg* MountingBlockHoleTubeSegShape =  */ 
						new TGeoTubeSeg("MountingBlockHoleTubeSegShape",
						0.0,
						fgkSSDMountingBlockScrewHoleRadius[0],
						0.5*fgkSSDMountingBlockScrewHoleHeigth,-90.,180.);
  TGeoCombiTrans* MountingBlockHoleTubeSegCombiTrans[MountingBlockHoleTubeSegNumber];
  char* MountingBlockHoleTubeSegCombiTransName[MountingBlockHoleTubeSegNumber] = 
					{	"MountingBlockHoleTubeSegCombiTrans0",
						"MountingBlockHoleTubeSegCombiTrans1",
						"MountingBlockHoleTubeSegCombiTrans2",
						"MountingBlockHoleTubeSegCombiTrans3"};
  for(Int_t i=0; i<MountingBlockHoleTubeSegNumber; i++){
    MountingBlockHoleTubeSegCombiTrans[i] =
      new TGeoCombiTrans(MountingBlockHoleTubeSegCombiTransName[i],
						0.5*fgkSSDMountingBlockScrewHoleEdge*TMath::Sqrt(2)
					*	TMath::Cos(45*(2*i+1)*TMath::DegToRad()),
						0.5*fgkSSDMountingBlockScrewHoleEdge*TMath::Sqrt(2)
					*	TMath::Sin(45*(2*i+1)*TMath::DegToRad()),
						0.0,
						new TGeoRotation("",90.*i,0.,0.));
    MountingBlockHoleTubeSegCombiTrans[i]->RegisterYourself();
  }
  TGeoBBox* MountingBlockHoleBoxShape = 
						new TGeoBBox("MountingBlockHoleBoxShape",
						0.5*fgkSSDMountingBlockScrewHoleEdge,
						0.5*fgkSSDMountingBlockScrewHoleEdge,
						0.5*fgkSSDMountingBlockScrewHoleHeigth);
  TGeoCompositeShape* MountingBlockScrewHole[2];
  MountingBlockScrewHole[0] = 
			new TGeoCompositeShape("MountingBlockScrewHole0",
			"MountingBlockHoleTubeSegShape:MountingBlockHoleTubeSegCombiTrans0+"
			"MountingBlockHoleTubeSegShape:MountingBlockHoleTubeSegCombiTrans1+"
			"MountingBlockHoleTubeSegShape:MountingBlockHoleTubeSegCombiTrans2+"
			"MountingBlockHoleTubeSegShape:MountingBlockHoleTubeSegCombiTrans3+"
			"MountingBlockHoleBoxShape");
  /* TGeoTubeSeg* MountingBlockLowerHoleTubeSegShape = */
			new TGeoTubeSeg("MountingBlockLowerHoleTubeSegShape",
						0.0,
						fgkSSDMountingBlockScrewHoleRadius[1],
						0.5*(fgkSSDMountingBlockHoleTubeWidth[0]
					-	fgkSSDMountingBlockScrewHoleHeigth
					-	fgkSSDMountingBlockHeight[3]),0.,90.); 
  TGeoCombiTrans* MountingBlockLowerHoleTubeSegCombiTrans[MountingBlockHoleTubeSegNumber];
  char* MountingBlockLowerHoleTubeSegCombiTransName[MountingBlockHoleTubeSegNumber] =
					{	"MountingBlockLowerHoleTubeSegCombiTrans0",
						"MountingBlockLowerHoleTubeSegCombiTrans1",
						"MountingBlockLowerHoleTubeSegCombiTrans2",
						"MountingBlockLowerHoleTubeSegCombiTrans3"};
  for(Int_t i=0; i<MountingBlockHoleTubeSegNumber; i++){
    MountingBlockLowerHoleTubeSegCombiTrans[i] =
			new TGeoCombiTrans(MountingBlockLowerHoleTubeSegCombiTransName[i],
						0.5*(fgkSSDMountingBlockScrewHoleEdge
					-	2.*fgkSSDMountingBlockScrewHoleRadius[1])
					*	TMath::Sqrt(2)*TMath::Cos(45*(2*i+1)*TMath::DegToRad()),
						0.5*(fgkSSDMountingBlockScrewHoleEdge
					-	2.0*fgkSSDMountingBlockScrewHoleRadius[1])
					*	TMath::Sqrt(2)*TMath::Sin(45*(2*i+1)*TMath::DegToRad()),0.,
						new TGeoRotation("",90.*i,0.,0.));
					MountingBlockLowerHoleTubeSegCombiTrans[i]->RegisterYourself();
  }
  Double_t fgkSSDMountingBlockLowerScrewHoleEdge = fgkSSDMountingBlockScrewHoleEdge
					-	2.*fgkSSDMountingBlockScrewHoleRadius[1];
  TGeoBBox* MountingBlockLowerHoleBoxShape[2];
  MountingBlockLowerHoleBoxShape[0] = 
			new TGeoBBox("MountingBlockLowerHoleBoxShape0",
						0.5*fgkSSDMountingBlockLowerScrewHoleEdge,
						0.5*fgkSSDMountingBlockLowerScrewHoleEdge,
						0.5*(fgkSSDMountingBlockHoleTubeWidth[0]
					-	fgkSSDMountingBlockScrewHoleHeigth
					-	fgkSSDMountingBlockHeight[3]));
  MountingBlockLowerHoleBoxShape[1] = 
			new TGeoBBox("MountingBlockLowerHoleBoxShape1",
						0.5*fgkSSDMountingBlockLowerScrewHoleEdge,
						0.5*fgkSSDMountingBlockScrewHoleRadius[1],
						0.5*(fgkSSDMountingBlockHoleTubeWidth[0]
					-	fgkSSDMountingBlockScrewHoleHeigth
					-	fgkSSDMountingBlockHeight[3]));
  TGeoCombiTrans* MountingBlockLowerHoleBoxCombiTrans[MountingBlockHoleTubeSegNumber];
  char* MountingBlockLowerHoleBoxCombiTransName[MountingBlockHoleTubeSegNumber] = 
					{	"MountingBlockLowerHoleBoxCombiTrans0",
						"MountingBlockLowerHoleBoxCombiTrans1",
						"MountingBlockLowerHoleBoxCombiTrans2",
						"MountingBlockLowerHoleBoxCombiTrans3"};
  for(Int_t i=0; i<MountingBlockHoleTubeSegNumber; i++){
    MountingBlockLowerHoleBoxCombiTrans[i] =
			new TGeoCombiTrans(MountingBlockLowerHoleBoxCombiTransName[i],
						0.5*(fgkSSDMountingBlockLowerScrewHoleEdge
					+	fgkSSDMountingBlockScrewHoleRadius[1])
					*	TMath::Cos(90*(i+1)*TMath::DegToRad()),
						0.5*(fgkSSDMountingBlockLowerScrewHoleEdge
					+	fgkSSDMountingBlockScrewHoleRadius[1])
					*	TMath::Sin(90*(i+1)*TMath::DegToRad()),0.,
						new TGeoRotation("",90.*i,0.,0.));
    MountingBlockLowerHoleBoxCombiTrans[i]->RegisterYourself();
  }
  MountingBlockScrewHole[1] = new TGeoCompositeShape("MountingBlockScrewHole1",
	"MountingBlockLowerHoleTubeSegShape:MountingBlockLowerHoleTubeSegCombiTrans0+"
	"MountingBlockLowerHoleTubeSegShape:MountingBlockLowerHoleTubeSegCombiTrans1+"
	"MountingBlockLowerHoleTubeSegShape:MountingBlockLowerHoleTubeSegCombiTrans2+"
	"MountingBlockLowerHoleTubeSegShape:MountingBlockLowerHoleTubeSegCombiTrans3+"
	"MountingBlockLowerHoleBoxShape0+"
	"MountingBlockLowerHoleBoxShape1:MountingBlockLowerHoleBoxCombiTrans0+"
	"MountingBlockLowerHoleBoxShape1:MountingBlockLowerHoleBoxCombiTrans1+"
	"MountingBlockLowerHoleBoxShape1:MountingBlockLowerHoleBoxCombiTrans2+"
	"MountingBlockLowerHoleBoxShape1:MountingBlockLowerHoleBoxCombiTrans3");
  TGeoTranslation* MountingBlockScrewHole1Trans = 
			new TGeoTranslation("MountingBlockScrewHole1Trans",0.,0.,
					-	MountingBlockLowerHoleBoxShape[0]->GetDZ()
					-	MountingBlockHoleBoxShape->GetDZ());
  MountingBlockScrewHole1Trans->RegisterYourself();
  /* TGeoCompositeShape* MountingBlockHole = */ 
			new TGeoCompositeShape("MountingBlockHole",
	"MountingBlockScrewHole0+MountingBlockScrewHole1:MountingBlockScrewHole1Trans");
  TGeoTranslation* MountingBlockHoleTrans = new TGeoTranslation("MountingBlockHoleTrans",
						0.5*fgkSSDMountingBlockLength[0]
					-	MountingBlockBoxShape[0]->GetDZ(),
						0.0,
						2.*MountingBlockBoxShape[2]->GetDZ()
					-	MountingBlockBoxShape[0]->GetDZ()
					+	fgkSSDMountingBlockHeight[3]
					-	MountingBlockHoleBoxShape->GetDZ());
  MountingBlockHoleTrans->RegisterYourself();
  TGeoCompositeShape* MountingBlockShape = new TGeoCompositeShape("MountingBlockShape",
			"MountingBlockMainShape-(MountingBlockTubeShape0:MountingBlockTubeTrans0+"
			"MountingBlockTubeShape0:MountingBlockTubeTrans1+"
			"MountingBlockTubeShape1:MountingBlockTubeTrans2+"
			"MountingBlockTubeShape1:MountingBlockTubeTrans3+"
			"HoleTrapezoidShape:HoleTrapezoidShapeCombiTrans+"
			"MountingBlockHole:MountingBlockHoleTrans)");
  TGeoVolume* SSDMountingBlock = new TGeoVolume("SSDMountingBlock",
			MountingBlockShape,fgkSSDMountingBlockMedium);
  return SSDMountingBlock;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetLadder(Int_t iLayer){
  TGeoVolume* Ladder = new TGeoVolumeAssembly(iLayer==5 ? "ITSssdLay5Ladd" 
														: "ITSssdLay6Ladd");
  TGeoVolume* LadderSegment[2];
  LadderSegment[0] = GetLadderSegment(0);
  LadderSegment[1] = GetLadderSegment(1);
  TList* EndLadderSegmentList = GetEndLadderSegment();
  Double_t BeamAxisTransVector = fgkCarbonFiberJunctionWidth;
  Int_t SSDLaySensorsNumber = (iLayer==5 ? 
						fgkSSDLay5SensorsNumber : 
						fgkSSDLay6SensorsNumber);
  for(Int_t i=0; i<SSDLaySensorsNumber; i++) Ladder->AddNode(i%2==0 ? 
						LadderSegment[iLayer==5 ? 0 : 1] : 
						LadderSegment[iLayer==5 ? 1 : 0],
						SSDLaySensorsNumber-i,new TGeoTranslation("",-0.5*fgkCarbonFiberTriangleLength,
						BeamAxisTransVector*i,0.));
  Ladder->AddNode((TGeoVolume*)EndLadderSegmentList->At(0),1,
						new TGeoTranslation("",-0.5*fgkCarbonFiberTriangleLength,
						fgkCarbonFiberJunctionWidth*SSDLaySensorsNumber,0.));
  Ladder->AddNode((TGeoVolume*)EndLadderSegmentList->At(1),1,
						new TGeoCombiTrans("",0.5*fgkCarbonFiberTriangleLength,
						0.,0.,new TGeoRotation("",180.,0.,0.)));
/////////////////////////////////////////////////////////////////////////////						
/// Placing Ladder Cables
/////////////////////////////////////////////////////////////////////////////		
  SetLadderCableCombiTransMatrix(iLayer);
  Int_t SideCableNumber[2] = {0,0};
  switch(iLayer){
	case 5: 
		SideCableNumber[0] = fgkSSDLay5SensorsNumber/2+1; 
		SideCableNumber[1] = SideCableNumber[0]-2;
		break;
    case 6:
		SideCableNumber[0] = (fgkSSDLay6SensorsNumber-1)/2+1;
		SideCableNumber[1] = SideCableNumber[0]-1;
		break;
  }
  const Double_t* CarbonFiberToModulePosition = 
							 LadderSegmentCombiTransMatrix[1]->GetTranslation();
  Double_t SSDEndLadderCableLength[4];
  SSDEndLadderCableLength[0] = CarbonFiberToModulePosition[1]
							 + fgkSSDSensorLength
							 - fgkSSDModuleStiffenerPosition[1]
							 - fgkSSDStiffenerWidth 
							 - fgkSSDFlexWidth[0]
							 + fgkEndLadderCarbonFiberLowerJunctionLength[1];
  SSDEndLadderCableLength[1] = CarbonFiberToModulePosition[1]
							 + fgkSSDModuleStiffenerPosition[1]
							 + fgkSSDStiffenerWidth
							 + fgkEndLadderCarbonFiberLowerJunctionLength[1];
  SSDEndLadderCableLength[2] = SSDEndLadderCableLength[1]
							 - fgkEndLadderCarbonFiberLowerJunctionLength[1]
							 + fgkEndLadderCarbonFiberLowerJunctionLength[0];
  SSDEndLadderCableLength[3] = fgkCarbonFiberJunctionWidth-(fgkSSDSensorLength
							 + CarbonFiberToModulePosition[1]
							 - fgkSSDModuleStiffenerPosition[1]
							 - fgkSSDStiffenerWidth)
							 + fgkEndLadderCarbonFiberLowerJunctionLength[1];
  TList* LadderCableAssemblyList[4];
  const Int_t EndLadderCablesNumber = 4;
  for(Int_t i=0; i<EndLadderCablesNumber; i++){
	  LadderCableAssemblyList[i] = 
	  GetLadderCableAssemblyList(SideCableNumber[i<2?0:1],
								   SSDEndLadderCableLength[i]);
	  Ladder->AddNode((TGeoVolume*)LadderCableAssemblyList[i]->At(i%2==0?0:1),
									i<2?1:2,LadderCableCombiTransMatrix[i]);
  }
  return Ladder;
}								  
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetLayer(Int_t iLayer){
  TGeoVolume* Layer = new TGeoVolumeAssembly(iLayer==5 ? "ITSssdLayer5" 
													   : "ITSssdLayer6");
  TGeoVolume* Ladder = GetLadder(iLayer);
  /////////////////////////////////////////////////////
  // Setting the CombiTransformation to pass ITS center 
  /////////////////////////////////////////////////////
  Double_t ITSCenterTransZ = iLayer==5 ? fgkSSDLay5LadderLength
                           -             fgkLay5CenterITSPosition: 
                                         fgkSSDLay6LadderLength
                           -             fgkLay6CenterITSPosition;
  Double_t ITSSensorYTrans = fgkSSDModuleCoolingBlockToSensor
                           + 0.5*        fgkCoolingTubeSupportHeight
                           -(iLayer==5 ? fgkSSDSensorSideSupportHeight[1]
                           -             fgkSSDSensorSideSupportHeight[0]: 0.);
  TGeoRotation* ITSCenterRot[3];
  ITSCenterRot[0] = new TGeoRotation("ITSCenterRot0",90.,180.,-90.);
  ITSCenterRot[1] = new TGeoRotation("ITSCenterRot1",0.,90.,0.);
  ITSCenterRot[2] = new TGeoRotation((*ITSCenterRot[1])*(*ITSCenterRot[0]));
  TGeoCombiTrans* ITSCenterCombiTrans = new TGeoCombiTrans("ITSCenterCombiTrans",0.,
                                        fgkSSDMountingBlockHeight[1]+ITSSensorYTrans,
                                        fgkEndLadderCarbonFiberLowerJunctionLength[1]
                                      - ITSCenterTransZ,ITSCenterRot[2]);
  /////////////////////////////////////////////////////
  // Setting the Ladder into the Layer 
  /////////////////////////////////////////////////////
  Double_t LayLadderAnglePosition = 360./(iLayer==5 ? fgkSSDLay5LadderNumber : 
                                          fgkSSDLay6LadderNumber);
  Int_t LayLadderNumber        = iLayer==5 ? fgkSSDLay5LadderNumber : 
                                             fgkSSDLay6LadderNumber; 
  Double_t LayRadiusMin        = iLayer==5 ? fgkSSDLay5RadiusMin    : 
                                             fgkSSDLay6RadiusMin;
  Double_t LayRadiusMax        = iLayer==5 ? fgkSSDLay5RadiusMax    : 
                                             fgkSSDLay6RadiusMax;
  TGeoCombiTrans* LadderCombiTrans[LayLadderNumber];
  TGeoHMatrix* LadderHMatrix[LayLadderNumber];
  Double_t LayerRadius = 0.;
  char LadderCombiTransName[30], LadderRotName[30];
  for(Int_t i=0; i<LayLadderNumber;i++){
      sprintf(LadderCombiTransName,"LadderCombiTrans%i",i);
      sprintf(LadderRotName,"LaddeRot%i",i);
      switch(iLayer){
            case 5: 
                  LayerRadius = (i%2==0 ? LayRadiusMin: LayRadiusMax);
                  break;
            case 6:
                  LayerRadius = (i%2==0 ? LayRadiusMax : LayRadiusMin);
                  break;
      }
      LadderCombiTrans[i] = new TGeoCombiTrans(LadderCombiTransName,
                            LayerRadius *	TMath::Cos((i+1)
                          * LayLadderAnglePosition*TMath::DegToRad()),
                            LayerRadius *	TMath::Sin((i+1)
                          * LayLadderAnglePosition*TMath::DegToRad()),0.,
                            new TGeoRotation(LadderRotName,(i+1)
                          * LayLadderAnglePosition-90,0.,0.));
    LadderHMatrix[i] = new TGeoHMatrix((*LadderCombiTrans[i])
                          * (*ITSCenterCombiTrans));
    Layer->AddNode(Ladder,i+1,LadderHMatrix[i]);          
  }   
  return Layer;
}
////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::Layer5(TGeoVolume* moth){
////////////////////////////////////////////////////////////////////////////////
// Insert the layer 5 in the mother volume. 
////////////////////////////////////////////////////////////////////////////////
	if(! moth){
		cout << "Error::AliITSv11GeometrySSD: Can't insert layer5, mother is null!"<< endl;
		return;
	}
	fMotherVol = moth;
	moth->AddNode(GetLayer(5),1,0);
}
////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::Layer6(TGeoVolume* moth){
////////////////////////////////////////////////////////////////////////////////
// Insert the layer 6 in the mother volume. 
////////////////////////////////////////////////////////////////////////////////
	if(! moth){
		cout << "Error::AliITSv11GeometrySSD: Can't insert layer6, mother is null!"<< endl;
		return;
	}
	fMotherVol = moth;
	moth->AddNode(GetLayer(6),1,0);
}
////////////////////////////////////////////////////////////////////////////////
TGeoArb8* AliITSv11GeometrySSD::GetTrapezoidShape(TVector3* vertexpos[], 
							Double_t* width, Double_t height, char* shapename){
////////////////////////////////////////////////////////////////////////////////
  const Int_t VertexNumber = 4;
  const Int_t TransVectNumber = 2;
  TVector3* Vertex[VertexNumber];
  TVector3* TransVector[2];
  for(Int_t i=0; i<TransVectNumber; i++) 
									TransVector[i] = new TVector3(0.,width[i]);
////////////////////////////////////////////////////////////////////////////////
  //Setting the vertices
////////////////////////////////////////////////////////////////////////////////
  Vertex[0] = new TVector3(*vertexpos[0]);
  Vertex[1] = new TVector3(*vertexpos[1]);
  Vertex[2] = new TVector3(*Vertex[1]+*TransVector[1]);
  Vertex[3] = new TVector3(*Vertex[0]+*TransVector[0]);
  TGeoArb8* TrapezoidShape = new TGeoArb8(shapename,0.5*height);
  for(Int_t i=0; i<2*VertexNumber; i++) 
	TrapezoidShape->SetVertex(i,Vertex[(i<VertexNumber ? i : i-VertexNumber)]->X(),
								Vertex[(i<VertexNumber ? i : i-VertexNumber)]->Y()); 
  return TrapezoidShape;
}
////////////////////////////////////////////////////////////////////////////////
TGeoArb8* AliITSv11GeometrySSD::GetArbShape(TVector3* vertexpos[], Double_t* width, 
									Double_t height, char* shapename, Int_t isign){
  const Int_t VertexNumber = 8;
  const Int_t TransVectNumber = 2;
  TVector3* Vertex[VertexNumber];
  TVector3* TransVector[2];
  for(Int_t i=0; i<TransVectNumber; i++) TransVector[i] = new TVector3(0.,width[i]);
////////////////////////////////////////////////////////////////////////////////
  //Setting the vertices for TGeoArb8
////////////////////////////////////////////////////////////////////////////////
  Vertex[0] = new TVector3(*vertexpos[0]);
  Vertex[1] = new TVector3(*vertexpos[1]);
  Vertex[2] = new TVector3(*Vertex[1]+isign*(*TransVector[0]));
  Vertex[3] = new TVector3(*Vertex[0]+isign*(*TransVector[0]));
  Vertex[4] = new TVector3(*vertexpos[2]);
  Vertex[5] = new TVector3(*vertexpos[3]);
  Vertex[6] = new TVector3(*Vertex[5]+isign*(*TransVector[1]));
  Vertex[7] = new TVector3(*Vertex[4]+isign*(*TransVector[1]));
////////////////////////////////////////////////////////////////////////////////
  TGeoArb8* ArbShape = new TGeoArb8(shapename,0.5*height);
  for(Int_t i = 0; i<VertexNumber;i++) 
							ArbShape->SetVertex(i,Vertex[i]->X(),Vertex[i]->Y());
  return ArbShape;
} 
////////////////////////////////////////////////////////////////////////////////
TGeoArb8* AliITSv11GeometrySSD::GetTriangleShape(TVector3* vertexpos[], 
											  Double_t height, char* shapename){
////////////////////////////////////////////////////////////////////////////////
  const Int_t VertexNumber = 4;
  TVector3* Vertex[VertexNumber];
//////////////////////////////////////
  //Setting the vertices for TGeoArb8
  ////////////////////////////////////
  for(Int_t i = 0; i<VertexNumber; i++)  
  Vertex[i] = new TVector3(i!=VertexNumber-1?*vertexpos[i]:*Vertex[VertexNumber-1-i]);
  TGeoArb8* TriangleShape = new TGeoArb8(shapename,0.5*height);
  for(Int_t i = 0; i<2*VertexNumber; i++) 
	TriangleShape->SetVertex(i,Vertex[(i < VertexNumber ? i: i-VertexNumber)]->X(),
							   Vertex[(i < VertexNumber ? i : i-VertexNumber)]->Y());
  return TriangleShape;
}
////////////////////////////////////////////////////////////////////////////////
TVector3* AliITSv11GeometrySSD::GetReflection(TVector3* Vector,Double_t* Param){
////////////////////////////////////////////////////////////////////////////////
  TVector3* n = new TVector3(Param[0],Param[1],Param[2]);
  Double_t D = ((*Vector)*(*n)+Param[3])/n->Mag2();
  TVector3* ReflectedVector = new TVector3(*Vector-2*D*(*n));
  return ReflectedVector;
}
////////////////////////////////////////////////////////////////////////////////
TGeoCombiTrans* AliITSv11GeometrySSD::AddTranslationToCombiTrans(TGeoCombiTrans* ct,
                                                       Double_t dx,
                                                       Double_t dy,
                                                       Double_t dz) const{
  TGeoCombiTrans* CombiTrans = new TGeoCombiTrans(*ct);
  // Add a dx,dy,dz translation to the initial TGeoCombiTrans
  const Double_t *vect = CombiTrans->GetTranslation();
  Double_t newVect[3] = {vect[0]+dx, vect[1]+dy, vect[2]+dz};
  CombiTrans->SetTranslation(newVect);
  return CombiTrans;
}
////////////////////////////////////////////////////////////////////////////////
//To be interfaced with AliRoot Materials file
////////////////////////////////////////////////////////////////////////////////
TGeoMedium* AliITSv11GeometrySSD::GetMedium(const char* mediumName) {
  char ch[30];
  sprintf(ch, "ITS_%s",mediumName);
  TGeoMedium* medium =  gGeoManager->GetMedium(ch);
  if (! medium)
    printf("Error(AliITSv11GeometrySSD)::medium %s not found !\n", mediumName);
  return medium;
}
////////////////////////////////////////////////////////////////////////////////
//To be interfaced with AliRoot Materials file
////////////////////////////////////////////////////////////////////////////////
/*void AliITSv11GeometrySSD::CreateMaterials(){
///////////////////////////////////
// This part has to be modified
///////////////////////////////////
//  Float_t ifield =  2.000;
//  Float_t fieldm = 10.000;
//    Int_t   ifield = gAlice->Field()->Integ();
//    Float_t fieldm = gAlice->Field()->Max();
  ///////////////////////////////////
  // Silicon for Sensor
  /////////////////////////////////// 
  Int_t SINumber = 0;
  TGeoMaterial* SiMaterial = new TGeoMaterial("SiMaterial");
  TGeoMedium* Silicon = new TGeoMedium("Silicon",SINumber,SiMaterial);
  fgkSSDSensorMedium = Silicon;
  ///////////////////////////////////
  // Silicon Mixture for Sensor
  /////////////////////////////////// 
  Int_t SiMixtureNumber = 1;
  TGeoMaterial* SiMixtureMaterial = new TGeoMaterial("SiMixtureMaterial");
  TGeoMedium* SiliconMixture = new TGeoMedium("SiliconMixture",SiMixtureNumber,SiMixtureMaterial);
  fgkSSDChipMedium = SiliconMixture;
  ///////////////////////////////////
  // Stiffener Components Materials
  /////////////////////////////////// 
  Int_t K1100Number = 2;
  TGeoMaterial* K1100Material = new TGeoMaterial("K1100Material");
  TGeoMedium* K1100 = new TGeoMedium("K1100",K1100Number,K1100Material);
  fgkSSDStiffenerMedium = K1100;
  ///////////////////////////  
  // Stiffener Connectors 
  ///////////////////////////  
  Int_t SnMaterialNumber = 3;
  TGeoMaterial* SnMaterial = new TGeoMaterial("SnMaterial");
  TGeoMedium* SnMedium = new TGeoMedium("SnMedium",SnMaterialNumber,
                                                SnMaterial);
  fgkSSDStiffenerConnectorMedium = SnMedium;
  ////////////////////////////////  
  // Stiffener 0603-1812 Capacitor
  ////////////////////////////////  
  Int_t Al2O3Number = 4;
  TGeoMaterial* Al2O3Material = new TGeoMaterial("Al2O3Material");
  TGeoMedium* Al2O3Medium = new TGeoMedium("Al2O3Medium",
                                                Al2O3Number,
                                                Al2O3Material);
  fgkSSDStiffener0603CapacitorMedium = Al2O3Medium;
  fgkSSDStiffener1812CapacitorMedium = Al2O3Medium;
  ///////////////////////////  
  // Stiffener Hybrid Wire 
  ///////////////////////////  
  Int_t CuHybridWireNumber = 5;
  TGeoMaterial* CuHybridWireMaterial = new TGeoMaterial("CuHybridWireMaterial");
  TGeoMedium* CuHybridWireMedium = new TGeoMedium("CuHybridWireMedium",
                                                CuHybridWireNumber,
                                                CuHybridWireMaterial);
  fgkSSDStiffenerHybridWireMedium = CuHybridWireMedium;
  ///////////////////////////  
  // Al for Cooling Block
  ///////////////////////////  
  Int_t AlCoolBlockNumber = 6;
  TGeoMaterial* AlCoolBlockMaterial = new TGeoMaterial("AlCoolBlockMaterial");
  TGeoMedium* AlCoolBlockMedium = new TGeoMedium("AlCoolBlockMedium",
                                                AlCoolBlockNumber,
                                                AlCoolBlockMaterial);
  fgkSSDAlCoolBlockMedium = AlCoolBlockMedium;
  //////////////////////////////////////////////////////  
  // Kapton and Al for Chip Cable Flex and Ladder Cables
  //////////////////////////////////////////////////////  
  Int_t KaptonBlockNumber = 7;
  TGeoMaterial* KaptonMaterial = new TGeoMaterial("KaptonMaterial");
  TGeoMedium* KaptonMedium = new TGeoMedium("KaptonMedium",
                                                KaptonBlockNumber,
                                                KaptonMaterial);
  Int_t AlTraceBlockNumber = 8;
  TGeoMaterial* AlTraceMaterial = new TGeoMaterial("AlTraceMaterial");
  TGeoMedium* AlTraceMedium = new TGeoMedium("AlTraceMedium",
                                                AlTraceBlockNumber,
                                                AlTraceMaterial);
  fgkSSDKaptonChipCableMedium = KaptonMedium;
  fgkSSDAlTraceChipCableMedium = AlTraceMedium;
  fgkSSDKaptonFlexMedium = KaptonMedium;
  fgkSSDAlTraceFlexMedium = AlTraceMedium;
  fgkSSDKaptonLadderCableMedium = KaptonMedium;
  fgkSSDAlTraceLadderCableMedium = AlTraceMedium;
  /////////////////////////////////////////////////////////////////  
  // M55J for Carbon Fiber, CarbonFiber Lower Support and Junction
  //////////////////////////////////////////////////////////////////  
  Int_t M55JNumber = 9;
  TGeoMaterial* M55JMaterial = new TGeoMaterial("M55JMaterial");
  TGeoMedium* M55JMedium = new TGeoMedium("M55JMedium",M55JNumber,
                                           M55JMaterial);
  fgkSSDCarbonFiberMedium = M55JMedium;
  fgkSSDMountingBlockMedium = M55JMedium;
  /////////////////////////////////////////////////////////////////  
  // G10 for Detector Leg, TubeHolder
  //////////////////////////////////////////////////////////////////  
  Int_t G10Number = 10;
  TGeoMaterial* G10Material = new TGeoMaterial("G10Material");
  TGeoMedium* G10Medium = new TGeoMedium("G10Medium",G10Number,
                                           G10Material);
  fgkSSDTubeHolderMedium = G10Medium;
  fgkSSDSensorSupportMedium = G10Medium;
  /////////////////////////////////////////////////////////////////  
  // Water and Phynox for Cooling Tube
  //////////////////////////////////////////////////////////////////  
  Int_t WaterNumber = 11;
  TGeoMaterial* WaterMaterial = new TGeoMaterial("WaterMaterial");
  TGeoMedium* WaterMedium = new TGeoMedium("WaterMedium",WaterNumber,
                                           WaterMaterial);
  fgkSSDCoolingTubeWater = WaterMedium;
  Int_t PhynoxNumber = 12;
  TGeoMaterial* PhynoxMaterial = new TGeoMaterial("PhynoxMaterial");
  TGeoMedium* PhynoxMedium = new TGeoMedium("PhynoxMedium",PhynoxNumber,
                                           PhynoxMaterial);
  fgkSSDCoolingTubePhynox = PhynoxMedium;
}
*/
void AliITSv11GeometrySSD::CreateMaterials(){
///////////////////////////////////
// This part has to be modified
///////////////////////////////////
  ///////////////////////////////////
  // Silicon for Sensor
  /////////////////////////////////// 
  fgkSSDSensorMedium = GetMedium("Si");
  ///////////////////////////////////
  // Silicon Mixture for Sensor
  /////////////////////////////////// 
  fgkSSDChipMedium = GetMedium("SPD SI CHIP$");
  fgkSSDChipGlueMedium = GetMedium("EPOXY$");
  ///////////////////////////////////
  // Stiffener Components Materials
  /////////////////////////////////// 
  fgkSSDStiffenerMedium = GetMedium("ITSsddCarbonM55J");
  ///////////////////////////  
  // Stiffener Connectors 
  ///////////////////////////  
  fgkSSDStiffenerConnectorMedium = GetMedium("COPPER");
  ////////////////////////////////  
  // Stiffener 0603-1812 Capacitor
  ////////////////////////////////  
  fgkSSDStiffener0603CapacitorMedium = GetMedium("SDD ruby sph. Al2O3");
  fgkSSDStiffener1812CapacitorMedium = GetMedium("SDD ruby sph. Al2O3");
  ///////////////////////////  
  // Stiffener Hybrid Wire 
  ///////////////////////////  
  fgkSSDStiffenerHybridWireMedium = GetMedium("COPPER");
  ///////////////////////////  
  // Al for Cooling Block
  ///////////////////////////  
  fgkSSDAlCoolBlockMedium = GetMedium("ITSal");
  //////////////////////////////////////////////////////  
  // Kapton and Al for Chip Cable Flex and Ladder Cables
  //////////////////////////////////////////////////////  
  fgkSSDKaptonChipCableMedium = GetMedium("KAPTONH(POLYCH2)");
  fgkSSDAlTraceChipCableMedium = GetMedium("ITSal");
  fgkSSDKaptonFlexMedium = GetMedium("KAPTONH(POLYCH2)");
  fgkSSDAlTraceFlexMedium = GetMedium("ITSal");
  fgkSSDKaptonLadderCableMedium = GetMedium("KAPTONH(POLYCH2)");
  fgkSSDAlTraceLadderCableMedium = GetMedium("ITSal");
  /////////////////////////////////////////////////////////////////  
  // M55J for Carbon Fiber, CarbonFiber Lower Support and Junction
  //////////////////////////////////////////////////////////////////  
  fgkSSDCarbonFiberMedium = GetMedium("GEN C (M55J)$");
  /////////////////////////////////////////////////////////////////  
  // G10 for Detector Leg, TubeHolder
  //////////////////////////////////////////////////////////////////  
  fgkSSDTubeHolderMedium = GetMedium("G10FR4$");
  fgkSSDSensorSupportMedium = GetMedium("G10FR4$");
  fgkSSDMountingBlockMedium = GetMedium("G10FR4$");
  fgkSSDMountingBlockMedium = GetMedium("G10FR4$");
  /////////////////////////////////////////////////////////////////  
  // Water and Phynox for Cooling Tube
  //////////////////////////////////////////////////////////////////  
  fgkSSDCoolingTubeWater = GetMedium("WATER");
  fgkSSDCoolingTubePhynox = GetMedium("INOX$");
}
/////////////////////////////////////////////////////////////////////
