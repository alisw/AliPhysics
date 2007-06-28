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
#include "TGeoMatrix.h"
#include <TGeoManager.h>
#include "TVector3.h"
#include "TGeoArb8.h"
#include "TList.h"
#include "TGeoMatrix.h"
#include "TGeoCompositeShape.h"
#include "TGeoTube.h"
#include "TGeoBBox.h"
#include "TGeoXtru.h"
#include "AliITSv11GeometrySSD.h"
/////////////////////////////////////////////////////////////////////////////////
// Names of the Sensitive Volumes of Layer 5 and Layer 6
/////////////////////////////////////////////////////////////////////////////////
const char* AliITSv11GeometrySSD::fgSDDsensitiveVolName5 = "ITSsddSensitivL5";
const char* AliITSv11GeometrySSD::fgSDDsensitiveVolName6 = "ITSsddSensitivL6";
/////////////////////////////////////////////////////////////////////////////////
//Parameters for SSD Geometry
/////////////////////////////////////////////////////////////////////////////////
// Layer5 (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDLay5LadderLength      = 950.7*fgkmm;
const Int_t AliITSv11GeometrySSD::fgkSSDLay5SensorsNumber        =  22;
const Int_t AliITSv11GeometrySSD::fgkSSDLay5LadderNumber         =  34;
const Double_t AliITSv11GeometrySSD::fgkSSDLay5RadiusMin         = 378.0*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDLay5RadiusMax         = 384.0*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkLay5CenterITSPosition    = 467.85*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// Layer6 (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDLay6LadderLength      = 1068.0*fgkmm;
const Int_t AliITSv11GeometrySSD::fgkSSDLay6SensorsNumber        =   25;
const Int_t AliITSv11GeometrySSD::fgkSSDLay6LadderNumber         =   38;
const Double_t AliITSv11GeometrySSD::fgkSSDLay6RadiusMin         =  428.0*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDLay6RadiusMax         =  434.0*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkLay6CenterITSPosition    = 526.50*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// SSD Chips and Hybrid (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Int_t AliITSv11GeometrySSD::fgkSSDChipNumber               =   6;
const Double_t AliITSv11GeometrySSD::fgkSSDChipLength            =  11.100*fgkmm; 
const Double_t AliITSv11GeometrySSD::fgkSSDChipWidth             =   3.850*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDChipHeight            =   0.180*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDChipSeparationLength  =   1.000*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDChipGlueLength     = fgkSSDChipLength;
const Double_t AliITSv11GeometrySSD::fgkSSDChipGlueWidth      =  fgkSSDChipWidth;
const Double_t AliITSv11GeometrySSD::fgkSSDChipGlueHeight        =   0.030*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// Stiffener (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerLength       =  73.000*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerWidth        =   6.500*fgkmm;
//const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerHeight       =   3.315;
const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerHeight       =   0.315*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerToChipDist   =   2.500*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor0603Length   = 1.600*fgkmm; 
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor0603Width    =   0.870*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor0603Height   =   0.800*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor1812Length   =   4.600*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor1812Width    =   3.400*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor1812Height   =   0.700*fgkmm; // multiplied by 0.5  
const Double_t AliITSv11GeometrySSD::fgkSSDWireLength            =  30.000*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDWireRadius            =   0.185*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorPosition[2]  = 
													   {44.32*fgkmm, 0.33*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorSeparation   =	  0.44*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorLength       =	  2.16*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorWidth        =	  3.60*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorHeight       = 
													  0.25*fgkSSDStiffenerHeight;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorAlHeight     =	 0.030*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorNiHeight     =   0.002*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// Cooling Block (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockLength    =   3.000*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockWidth     =   4.000*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockHeight[3] =  
										 {1.950*fgkmm, 0.240*fgkmm, 0.300*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockHoleRadius[2] = 
													  {1.000*fgkmm, 0.120*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockHoleLength[2] = 
													  {1.900*fgkmm, 0.400*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockHoleCenter    =  
																	 1.500*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCoolingBlockHoleHeight    =  
																	 0.300*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// SSD Sensor (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const char*  AliITSv11GeometrySSD::fgkSSDSensitiveVolName       = 
														 "SSDSensorSensitiveVol";
const Double_t AliITSv11GeometrySSD::fgkSSDSensorLength          =  42.000*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorHeight          =   0.300*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorWidth           =  75.000*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorOverlap         = 
	   											   fgkSSDSensorLength-39.1*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorInsensitiveLength    = 1.*fgkmm; 
const Double_t AliITSv11GeometrySSD::fgkSSDSensorInsensitiveWidth     = 1.*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// Flex (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDFlexFullLength       =  106.000*fgkmm; 
const Double_t AliITSv11GeometrySSD::fgkSSDFlexLength[4]        = 
			{0.5 * (fgkSSDStiffenerLength+fgkSSDChipNumber*fgkSSDChipLength
				 + (fgkSSDChipNumber-1)*fgkSSDChipSeparationLength),
			 0.5 * (fgkSSDStiffenerLength+fgkSSDChipNumber*fgkSSDChipLength
				 + (fgkSSDChipNumber-1)*fgkSSDChipSeparationLength)
									   - 4.000*fgkmm, 9.500*fgkmm, 10.000*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDFlexWidth[2]         = 
												   {  9.340*fgkmm,  5.380*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDFlexHeight[2]        =
												   {  0.030*fgkmm,  0.020*fgkmm};      
const Double_t AliITSv11GeometrySSD::fgkSSDFlexAngle            =   30.000;
const Double_t AliITSv11GeometrySSD::fgkSSDFlexHoleLength       =    1.430*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDFlexHoleWidth        =    3.000*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDEndFlexCompLength[6] = 
			 {3.30*fgkmm,4.12*fgkmm,4.22*fgkmm,1.70*fgkmm,0.75*fgkmm,7.18*fgkmm};
const Double_t AliITSv11GeometrySSD:: fgkSSDEndFlexCompWidth[3] =
										   {15.03*fgkmm,23.48*fgkmm,12.28*fgkmm};
/////////////////////////////////////////////////////////////////////////////////
// SSD Ladder Cable (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDLadderCableWidth     =     23.5*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// SSD Module (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDModuleStiffenerPosition[2]  = 
													 { 1.000*fgkmm, 3.900*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDModuleSensorSupportDistance =  
																	45.600*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDModuleCoolingBlockToSensor  =  
																	 5.075*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// Sensor Support (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportLength		       = 
																	 5.800*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportWidth          =  
																	 2.000*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportHeight[2]      =
												     { 4.620*fgkmm, 5.180*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportThickness[2]   = 
													 { 0.450*fgkmm, 0.450*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportPosition       
								= 0.5 * (fgkSSDModuleSensorSupportDistance
							    +  fgkSSDSensorSideSupportThickness[0])
								-  fgkSSDSensorSideSupportLength;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportLength	   =  
									   							    5.250*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportWidth        =
																	1.680*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportHeight[2]    
								  = {fgkSSDSensorSideSupportHeight[0]
								  +  fgkSSDSensorSideSupportThickness[0],
									 fgkSSDSensorSideSupportHeight[1]
								  +  fgkSSDSensorSideSupportThickness[1]};
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportThickness[2] 
								  =  {fgkSSDSensorSideSupportThickness[0],
									  fgkSSDSensorSideSupportThickness[1]};
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportPosition     = 
																   19.000*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// Chip Cables (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDChipCablesLength[2]   = 
				  {73.12/fgkSSDChipNumber*fgkmm,fgkSSDChipLength+2.*0.19*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDChipCablesHeight[4]   = 
				  {  0.014*fgkmm,  0.010*fgkmm, fgkSSDModuleCoolingBlockToSensor
								  -  (fgkSSDSensorSideSupportHeight[1]
								  -   fgkSSDSensorSideSupportHeight[0])
								  -   fgkSSDCoolingBlockHoleCenter
								  -   fgkSSDStiffenerHeight
								  -   fgkSSDChipHeight-fgkSSDSensorHeight,
									  fgkSSDModuleCoolingBlockToSensor
								  -   fgkSSDCoolingBlockHoleCenter
								  -	  fgkSSDStiffenerHeight
								  -   fgkSSDChipHeight-fgkSSDSensorHeight};
const Double_t AliITSv11GeometrySSD::fgkSSDChipCablesWidth[3]    = 
		                            { 11.000*fgkmm,  0.800*fgkmm,  0.600*fgkmm};
/////////////////////////////////////////////////////////////////////////////////
// Carbon Fiber Junction Parameters (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberJunctionLength          = 
																	3.820*fgkmm;
//const Double_t AliITSv11GeometrySSD::fgkCarbonFiberJunctionLength          = 
//																	   3.780;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberJunctionWidth           = 
										 fgkSSDSensorLength-fgkSSDSensorOverlap;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberJunctionEdge[2]         = 
													{  0.86*fgkmm,  0.30*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberJunctionAngle[2]        = 
																{ 30.00, 90.00};
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberJunctionToSensorSupport = 
																	 1.78*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
//Carbon Fiber Parameters (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberTriangleLength          
								   = fgkSSDModuleSensorSupportDistance
								   - 2. * fgkCarbonFiberJunctionToSensorSupport;  
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberTriangleAngle = 60.00;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberSupportTopEdgeDist[2]   = 
												  {  0.751*fgkmm,  0.482*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberSupportEdgeLength  = 
																	1.630*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberSupportWidth =   0.950*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberSupportXAxisLength      
									= fgkCarbonFiberTriangleLength
									- 0.5*fgkCarbonFiberSupportTopEdgeDist[1]
									/ TMath::Cos(fgkCarbonFiberTriangleAngle
									* TMath::DegToRad());
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberSupportYAxisLength      
									= 0.5*(fgkCarbonFiberJunctionWidth
									- fgkCarbonFiberSupportWidth)
									- fgkCarbonFiberSupportTopEdgeDist[0]
									- fgkCarbonFiberSupportWidth;
/////////////////////////////////////////////////////////////////////////////////
// Carbon Fiber Lower Support Parameters (lengths are in mm)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberLowerSupportWidth             
																	  =  0.950*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberLowerSupportLowerLenght       
																	  =  1.600*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberLowerSupportHeight            
																	  =  0.830*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberLowerSupportVolumeSeparation  
											  = 0.5*fgkCarbonFiberSupportWidth;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberLowerSupportTransverseWidth   
									= fgkCarbonFiberJunctionWidth
									- 2. * (fgkCarbonFiberLowerSupportWidth
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
							   {0.5 * (fgkSSDLay5LadderLength
									-  fgkSSDLay5SensorsNumber
									*  fgkCarbonFiberJunctionWidth
									-  fgkCarbonFiberLowerSupportWidth),
								0.5 * (fgkSSDLay5LadderLength
									-  fgkSSDLay5SensorsNumber
									*  fgkCarbonFiberJunctionWidth
									+  fgkCarbonFiberLowerSupportWidth)};
const Double_t AliITSv11GeometrySSD::fgkEndLadderCarbonFiberUpperJunctionLength[2] = 
						{fgkEndLadderCarbonFiberLowerJunctionLength[0]-20.4*fgkmm,
						 fgkEndLadderCarbonFiberLowerJunctionLength[1]-20.6*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkEndLadderMountingBlockPosition[2] = 
						{fgkEndLadderCarbonFiberLowerJunctionLength[0]-16.50*fgkmm,
						 fgkEndLadderCarbonFiberLowerJunctionLength[1]-31.50*fgkmm};
/////////////////////////////////////////////////////////////////////////////////
// Cooling Tube Support (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportRmax      =  1.45*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportRmin          
											  = fgkSSDCoolingBlockHoleRadius[0];
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportLength    =  8.55*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportHeight    =  0.85*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportWidth     =  2.00*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportSeparation    = 
                                        fgkSSDSensorLength-2.*fgkSSDSensorOverlap;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportToCarbonFiber = 
																	  11.70*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// Cooling Tube (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeRmax = fgkCoolingTubeSupportRmin;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeRmin =  0.96*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeLength = 
													fgkCarbonFiberJunctionWidth;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSeparation = 
									 fgkSSDModuleSensorSupportDistance
								  +	 fgkSSDCoolingBlockLength;
//const Double_t AliITSv11GeometrySSD_ct::fgkCoolingTubeLength               =  39.1;
/////////////////////////////////////////////////////////////////////////////////
// SSD Mounting Block Parameters (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockLength[3]            = 
										   { 60.0*fgkmm, 42.2*fgkmm, 34.0*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHeight[4]            =
							   {  4.0*fgkmm,  8.0*fgkmm,  5.0*fgkmm,  0.2*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockWidth                =   
																	  20.0*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHoleTrapezoidAngle   =   
																		    40.0;
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHoleTrapezoidHeight  = 
	           0.30*(fgkSSDMountingBlockHeight[1]-fgkSSDMountingBlockHeight[2]);
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHoleTrapezoidUpBasis =    
																	  2.5*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHoleTubeLength[2]    = 
													  { 56.0*fgkmm, 12.0*fgkmm}; 
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHoleTubeWidth[2]     = 
												      {  5.0*fgkmm,  2.9*fgkmm}; 
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockHoleRadius           = 
																	  1.0*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockScrewHoleEdge        =   
																	  6.0*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockScrewHoleHeigth      =  
																	  4.0*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDMountingBlockScrewHoleRadius[2]   =
							  {  1.5*fgkmm,fgkSSDMountingBlockScrewHoleEdge/6.};
/////////////////////////////////////////////////////////////////////////////////
ClassImp(AliITSv11GeometrySSD)
/////////////////////////////////////////////////////////////////////////////////
AliITSv11GeometrySSD::AliITSv11GeometrySSD(): 
  AliITSv11Geometry(),
  fSSDChipMedium(),
  fSSDChipGlueMedium(),
  fSSDStiffenerMedium(),
  fSSDStiffenerConnectorMedium(),
  fSSDStiffener0603CapacitorMedium(),
  fSSDStiffener1812CapacitorMedium(),
  fSSDStiffenerHybridWireMedium(),
  fSSDKaptonFlexMedium(),
  fSSDAlTraceFlexMedium(),
  fSSDAlTraceLadderCableMedium(),
  fSSDKaptonLadderCableMedium(),
  fSSDKaptonChipCableMedium(),
  fSSDAlTraceChipCableMedium(),
  fSSDAlCoolBlockMedium(),
  fSSDSensorMedium(),
  fSSDSensorSupportMedium(),
  fSSDCarbonFiberMedium(),
  fSSDTubeHolderMedium(),
  fSSDCoolingTubeWater(),
  fSSDCoolingTubePhynox(),
  fSSDMountingBlockMedium(),
  fSSDAir(),
  fCreateMaterials(kFALSE),
  fTransformationMatrices(kFALSE),
  fBasicObjects(kFALSE),
  fcarbonfiberjunction(),
  fcoolingtubesupport(),
  fhybridmatrix(),
  fssdcoolingblocksystem(),
  fcoolingblocksystematrix(),
  fssdstiffenerflex(),
  fssdendflex(),
  fendladdermountingblock(),
  fSSDSensor5(),
  fSSDSensor6(),
  fSSDLayer5(),	
  fSSDLayer6(),
  fMotherVol(),
  fColorCarbonFiber(4),
  fColorRyton(5),
  fColorPhynox(14),
  fColorSilicon(3),
  fColorAl(38),
  fColorKapton(6),
  fColorPolyhamide(5),
  fColorStiffener(9),
  fColorEpoxy(30),
  fColorWater(7),
  fColorG10(41)
{
  ////////////////////////
  // Standard constructor
  ////////////////////////
}
/////////////////////////////////////////////////////////////////////////////////
AliITSv11GeometrySSD::AliITSv11GeometrySSD(const AliITSv11GeometrySSD &s):
  AliITSv11Geometry(s.GetDebug()),
  fSSDChipMedium(s.fSSDChipMedium),
  fSSDChipGlueMedium(s.fSSDChipGlueMedium),
  fSSDStiffenerMedium(s.fSSDStiffenerMedium),
  fSSDStiffenerConnectorMedium(s.fSSDStiffenerConnectorMedium),
  fSSDStiffener0603CapacitorMedium(s.fSSDStiffener0603CapacitorMedium),
  fSSDStiffener1812CapacitorMedium(s.fSSDStiffener1812CapacitorMedium),
  fSSDStiffenerHybridWireMedium(s.fSSDStiffenerHybridWireMedium),
  fSSDKaptonFlexMedium(s.fSSDKaptonFlexMedium),
  fSSDAlTraceFlexMedium(s.fSSDAlTraceFlexMedium),
  fSSDAlTraceLadderCableMedium(s.fSSDAlTraceLadderCableMedium),
  fSSDKaptonLadderCableMedium(s.fSSDKaptonLadderCableMedium),
  fSSDKaptonChipCableMedium(s.fSSDKaptonChipCableMedium),
  fSSDAlTraceChipCableMedium(s.fSSDAlTraceChipCableMedium),
  fSSDAlCoolBlockMedium(s.fSSDAlCoolBlockMedium),
  fSSDSensorMedium(s.fSSDSensorMedium),
  fSSDSensorSupportMedium(s.fSSDSensorSupportMedium),
  fSSDCarbonFiberMedium(s.fSSDCarbonFiberMedium),
  fSSDTubeHolderMedium(s.fSSDTubeHolderMedium),
  fSSDCoolingTubeWater(s.fSSDCoolingTubeWater),
  fSSDCoolingTubePhynox(s.fSSDCoolingTubePhynox),
  fSSDMountingBlockMedium(s.fSSDMountingBlockMedium),
  fSSDAir(s.fSSDAir),
  fCreateMaterials(s.fCreateMaterials),
  fTransformationMatrices(s.fTransformationMatrices),
  fBasicObjects(s.fBasicObjects),
  fcarbonfiberjunction(s.fcarbonfiberjunction),
  fcoolingtubesupport(s.fcoolingtubesupport),
  fhybridmatrix(s.fhybridmatrix),
  fssdcoolingblocksystem(s.fssdcoolingblocksystem),
  fcoolingblocksystematrix(s.fcoolingblocksystematrix),
  fssdstiffenerflex(s.fssdstiffenerflex),
  fssdendflex(s.fssdendflex),
  fendladdermountingblock(s.fendladdermountingblock),
  fSSDSensor5(s.fSSDSensor5),
  fSSDSensor6(s.fSSDSensor6),
  fSSDLayer5(s.fSSDLayer5),	
  fSSDLayer6(s.fSSDLayer6),
  fMotherVol(s.fMotherVol),
  fColorCarbonFiber(s.fColorCarbonFiber),
  fColorRyton(s.fColorRyton),
  fColorPhynox(s.fColorPhynox),
  fColorSilicon(s.fColorSilicon),
  fColorAl(s.fColorAl),
  fColorKapton(s.fColorKapton),
  fColorPolyhamide(s.fColorPolyhamide),
  fColorStiffener(s.fColorStiffener),
  fColorEpoxy(s.fColorEpoxy),
  fColorWater(s.fColorWater),
  fColorG10(s.fColorG10)
{
  ////////////////////////
  // Copy Constructor
  ////////////////////////
}
/////////////////////////////////////////////////////////////////////////////////
AliITSv11GeometrySSD& AliITSv11GeometrySSD::
operator=(const AliITSv11GeometrySSD &s){
  ////////////////////////
  // Assignment operator
  ////////////////////////
  this->~AliITSv11GeometrySSD();
  new(this) AliITSv11GeometrySSD(s); 
  return *this;
/*	
  if(&s == this) return *this;
  fMotherVol = s.fMotherVol;
  return *this;
 */
}
///////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::CreateTransformationMatrices(){
  ///////////////////////////////////////////////////////////////////////  
  // Method generating the trasformation matrix for the whole SSD Geometry   
  ///////////////////////////////////////////////////////////////////////  
  // Setting some variables for Carbon Fiber Supportmatrix creation
  //////////////////////////////////////////////////////////////////////
  Double_t carbonfibersupportxaxisEdgeproj = fgkCarbonFiberSupportEdgeLength
										   * CosD(fgkCarbonFiberJunctionAngle[0]);
  Double_t symmetryplaneposition = (fgkCarbonFiberSupportYAxisLength
								 +	fgkCarbonFiberSupportTopEdgeDist[0]
								 +	fgkCarbonFiberSupportWidth);
  Double_t carbonfibersupportheight = carbonfibersupportxaxisEdgeproj
									* TanD(fgkCarbonFiberJunctionAngle[0]);
  TGeoRotation* carbonfiberot[3];
  for(Int_t i=0; i<3; i++) carbonfiberot[i] = new TGeoRotation();
  carbonfiberot[0]->SetAngles(0.0,180.0,0.0);
  carbonfiberot[1]->SetAngles(90.,-fgkCarbonFiberTriangleAngle,-90.);
  carbonfiberot[2]->SetRotation((*carbonfiberot[1])*(*carbonfiberot[0]));
  Double_t transvector[3] = {fgkCarbonFiberTriangleLength
						  *  CosD(fgkCarbonFiberTriangleAngle),0.,
						  -  fgkCarbonFiberTriangleLength
						  *	 SinD(fgkCarbonFiberTriangleAngle)};
  ///////////////////////////////////////////
  //Setting Local Translations and Rotations: 
  ///////////////////////////////////////////
  TGeoCombiTrans* localcarbonfibersupportmatrix[3]; 
  localcarbonfibersupportmatrix[0] = new TGeoCombiTrans(0.0,0.0,
									 0.5*carbonfibersupportheight,NULL);	
  localcarbonfibersupportmatrix[1] = new TGeoCombiTrans(transvector[0],
									 2.*symmetryplaneposition+transvector[1],
									 transvector[2], carbonfiberot[2]);
  localcarbonfibersupportmatrix[2] = new TGeoCombiTrans(*carbonfiberot[1]);
  /////////////////////////////////////////////////////////////
  // Carbon Fiber Support Transformations
  /////////////////////////////////////////////////////////////
  const Int_t kcarbonfibersupportmatrixnumber[2] = {2,3};
  for(Int_t i=0; i<fgkcarbonfibersupportnumber; i++){ 
		fcarbonfibersupportmatrix[i] = new TGeoHMatrix();
		for(Int_t j=0; j<kcarbonfibersupportmatrixnumber[i]; j++)
			fcarbonfibersupportmatrix[i]->MultiplyLeft(localcarbonfibersupportmatrix[i==0?2*j:j]);
  }		
  /////////////////////////////////////////////////////////////
  // Carbon Fiber Junction Transformation
  /////////////////////////////////////////////////////////////
  const Int_t kcarbonfiberjunctionmatrixnumber = 2;
  TGeoCombiTrans** localcarbonfiberjunctionmatrix[fgkcarbonfiberjunctionumber];
  TGeoRotation** localcarbonfiberjunctionrot[fgkcarbonfiberjunctionumber];
  TGeoTranslation** localcarbonfiberjunctiontrans[fgkcarbonfiberjunctionumber];
  for(Int_t i=0; i<fgkcarbonfiberjunctionumber; i++){ 
	localcarbonfiberjunctionmatrix[i] = 
						new TGeoCombiTrans*[kcarbonfiberjunctionmatrixnumber];
	localcarbonfiberjunctionrot[i] = 
						new TGeoRotation*[kcarbonfiberjunctionmatrixnumber];
	localcarbonfiberjunctiontrans[i] = 
						new TGeoTranslation*[kcarbonfiberjunctionmatrixnumber];
  }
  ///////////////////////
  // Setting Translations
  ///////////////////////
  localcarbonfiberjunctiontrans[0][0] = new TGeoTranslation(0.,0.,0.);
  localcarbonfiberjunctiontrans[1][0] = 
				 new TGeoTranslation(fgkCarbonFiberTriangleLength,0.0,0.0);
  localcarbonfiberjunctiontrans[2][0] = 
				 new TGeoTranslation(fgkCarbonFiberTriangleLength
				 * TMath::Cos(fgkCarbonFiberTriangleAngle*TMath::DegToRad()),
				 fgkCarbonFiberTriangleLength
				 * TMath::Sin(fgkCarbonFiberTriangleAngle*TMath::DegToRad()),0.0);
  localcarbonfiberjunctiontrans[0][1] = 
				 new TGeoTranslation(0.0,0.5*fgkCarbonFiberJunctionWidth,0.0);
  localcarbonfiberjunctiontrans[1][1] = 
				 new TGeoTranslation(*localcarbonfiberjunctiontrans[0][1]);
  localcarbonfiberjunctiontrans[2][1] = 
				 new TGeoTranslation(*localcarbonfiberjunctiontrans[0][1]);
  ////////////////////
  // Setting Rotations
  ////////////////////
  for(Int_t i=0; i<fgkcarbonfiberjunctionumber; i++)
		for(Int_t j=0; j<kcarbonfiberjunctionmatrixnumber; j++)
			localcarbonfiberjunctionrot[i][j] = new TGeoRotation();
  for(Int_t i=0; i<fgkcarbonfiberjunctionumber; i++)
	localcarbonfiberjunctionrot[i][0]->SetAngles(120.*i,0.0,0.0);
  localcarbonfiberjunctionrot[0][1]->SetAngles(0.0,90.0,0.0);
  localcarbonfiberjunctionrot[1][1]->SetRotation(*localcarbonfiberjunctionrot[0][1]);
  localcarbonfiberjunctionrot[2][1]->SetRotation(*localcarbonfiberjunctionrot[0][1]);
  ////////////////////////////////////////
  // Setting Carbon Fiber Junction matrix 
  ////////////////////////////////////////
  for(Int_t i=0; i<fgkcarbonfiberjunctionumber; i++){
		fcarbonfiberjunctionmatrix[i] = new TGeoHMatrix();
		for(Int_t j=0; j<kcarbonfiberjunctionmatrixnumber; j++){
			localcarbonfiberjunctionmatrix[i][j] = 
			new TGeoCombiTrans(*localcarbonfiberjunctiontrans[i][j],
							   *localcarbonfiberjunctionrot[i][j]);
		    fcarbonfiberjunctionmatrix[i]->MultiplyLeft(localcarbonfiberjunctionmatrix[i][j]);
	 }
  }
  /////////////////////////////////////////////////////////////
  // Carbon Fiber Lower Support Transformations
  /////////////////////////////////////////////////////////////
  TGeoTranslation* localcarbonfiberlowersupportrans[2];
  localcarbonfiberlowersupportrans[0] = new TGeoTranslation(0.0,
									 fgkCarbonFiberLowerSupportVolumePosition[1]
								+    fgkCarbonFiberLowerSupportVolumePosition[0],
									 0.0);
  localcarbonfiberlowersupportrans[1] = new TGeoTranslation(0.0,
									 fgkCarbonFiberJunctionWidth
								-    fgkCarbonFiberLowerSupportWidth
								-    fgkCarbonFiberLowerSupportVolumePosition[0]
								-    fgkCarbonFiberLowerSupportVolumePosition[1],
								-    0.5*fgkCarbonFiberLowerSupportHeight);
   localcarbonfiberlowersupportrans[0]->Add(localcarbonfiberlowersupportrans[1]);
   fcarbonfiberlowersupportrans[0] = 
						new TGeoTranslation(*localcarbonfiberlowersupportrans[0]);
   fcarbonfiberlowersupportrans[1] = 
						new TGeoTranslation(*localcarbonfiberlowersupportrans[1]);
  /////////////////////////////////////////////////////////////
  // SSD Sensor Support Transformations
  /////////////////////////////////////////////////////////////
  const Int_t kssdsensorsupportmatrixnumber = 3;
  TGeoCombiTrans** localssdsensorsupportmatrix[fgkssdsensorsupportnumber];
  TGeoRotation** localssdsensorsupportrot[fgkssdsensorsupportnumber];
  TGeoTranslation** localssdsensorsupportrans[fgkssdsensorsupportnumber];
  for(Int_t i=0; i<fgkssdsensorsupportnumber; i++){ 
	localssdsensorsupportmatrix[i] = 
						new TGeoCombiTrans*[kssdsensorsupportmatrixnumber];
	localssdsensorsupportrot[i] = 
						new TGeoRotation*[kssdsensorsupportmatrixnumber];
	localssdsensorsupportrans[i] = 
						new TGeoTranslation*[kssdsensorsupportmatrixnumber];
  }
  ///////////////////////
  // Setting Translations
  ///////////////////////
  localssdsensorsupportrans[0][0] = new TGeoTranslation(0.0,
									  0.5*fgkSSDSensorSideSupportWidth,
									  0.0);
  localssdsensorsupportrans[1][0] = 
						 new TGeoTranslation(*localssdsensorsupportrans[0][0]);
  localssdsensorsupportrans[2][0] = 
						 new TGeoTranslation(*localssdsensorsupportrans[0][0]);
  localssdsensorsupportrans[0][1] = 
						 new TGeoTranslation(-0.5*fgkSSDSensorSideSupportWidth,
										0.5*fgkSSDSensorSideSupportThickness[0],
										0.0);
  localssdsensorsupportrans[1][1] = 
						 new TGeoTranslation(0.5*fgkSSDSensorSideSupportWidth,
									-   0.5*fgkSSDSensorSideSupportThickness[0]
								    -   fgkSSDModuleSensorSupportDistance,
										0.0);
  localssdsensorsupportrans[2][1] = 
						 new TGeoTranslation(0.5*fgkSSDSensorCenterSupportThickness[0]
									-    fgkSSDSensorCenterSupportPosition,
										 0.5*fgkSSDSensorCenterSupportWidth
									-    0.5*fgkSSDModuleSensorSupportDistance,
										 fgkSSDSensorCenterSupportThickness[0]);
  localssdsensorsupportrans[0][2] = 
						 new TGeoTranslation(fgkCarbonFiberTriangleLength
									+    fgkCarbonFiberJunctionToSensorSupport,
										 fgkCarbonFiberJunctionWidth
								    -    0.5*(fgkCarbonFiberLowerSupportWidth
									+    fgkSSDSensorCenterSupportLength
									-    fgkSSDSensorCenterSupportThickness[0])
									-    fgkSSDSensorCenterSupportPosition,
									     0.0);
  localssdsensorsupportrans[1][2] = 
						new TGeoTranslation(*localssdsensorsupportrans[0][2]);
  localssdsensorsupportrans[2][2] = 
						new TGeoTranslation(*localssdsensorsupportrans[0][2]);
  ////////////////////
  // Setting Rotations
  ////////////////////
  for(Int_t i=0; i<fgkssdsensorsupportnumber; i++)
		for(Int_t j=0; j<kssdsensorsupportmatrixnumber; j++)
			localssdsensorsupportrot[i][j] = new TGeoRotation();
  for(Int_t i=0; i<fgkssdsensorsupportnumber; i++){
	localssdsensorsupportrot[i][0]->SetAngles(0.0,90.0,0.0);
	localssdsensorsupportrot[i][2]->SetAngles(-90.0,0.0,0.0);
  }
  localssdsensorsupportrot[0][1]->SetAngles(0.0,90.0,-90.0);
  localssdsensorsupportrot[1][1]->SetAngles(180.0,90.0,-90.0);
  localssdsensorsupportrot[2][1]->SetAngles(270.0,90.0,-90.0);
  ////////////////////////////////////////
  // SSD Sensor Support matrix 
  ////////////////////////////////////////
  for(Int_t i=0; i<fgkssdsensorsupportnumber; i++){
		fssdsensorsupportmatrix[i] = new TGeoHMatrix();
		for(Int_t j=0; j<kssdsensorsupportmatrixnumber; j++){
			localssdsensorsupportmatrix[i][j] = 
			new TGeoCombiTrans(*localssdsensorsupportrans[i][j],
							   *localssdsensorsupportrot[i][j]);
		    fssdsensorsupportmatrix[i]->MultiplyLeft(localssdsensorsupportmatrix[i][j]);
	 }
  }
  /////////////////////////////////////////////////////////////
  // SSD Cooling Tube Support Transformations
  /////////////////////////////////////////////////////////////
  const Int_t kcoolingtubesupportmatrixnumber = 2;
  TGeoCombiTrans* localcoolingtubesupportmatrix[kcoolingtubesupportmatrixnumber];
  TGeoTranslation* localcoolingtubesupportrans[kcoolingtubesupportmatrixnumber];
  TGeoRotation* localcoolingtubesupportrot[kcoolingtubesupportmatrixnumber];
  Double_t phi = TMath::ASin(0.5*fgkCoolingTubeSupportHeight
													/fgkCoolingTubeSupportRmax);
  localcoolingtubesupportrans[0] = 
			new TGeoTranslation(2.*fgkCoolingTubeSupportRmax*TMath::Cos(phi)
						+  2.*(fgkCoolingTubeSupportLength
						-  fgkCoolingTubeSupportRmax*(1.+TMath::Cos(phi)))
						+  fgkCarbonFiberTriangleLength
						-  2.0*fgkCarbonFiberJunctionLength,0.0,0.0);
  localcoolingtubesupportrans[1] = 
			new TGeoTranslation(fgkCarbonFiberJunctionLength
					- (fgkCoolingTubeSupportLength-fgkCoolingTubeSupportRmax),
					- (2.0*fgkSSDSensorLength-fgkSSDSensorOverlap)+
						   fgkSSDModuleStiffenerPosition[1]+fgkSSDStiffenerWidth
                    +  0.5*fgkSSDFlexHoleLength+2.*fgkCarbonFiberJunctionWidth
                    -  0.5*(fgkCarbonFiberLowerSupportWidth
					+	   fgkSSDSensorCenterSupportLength
                    -      fgkSSDSensorCenterSupportThickness[0])
					+  0.5*fgkSSDSensorLength,
					-  0.5*fgkCoolingTubeSupportHeight);  
  for(Int_t i=0; i<kcoolingtubesupportmatrixnumber; i++)     
  localcoolingtubesupportrot[i] = new TGeoRotation();
  localcoolingtubesupportrot[0]->SetAngles(180.0,0.0,0.0);
  localcoolingtubesupportrot[1]->SetAngles(0.0,90.0,0.0);
  for(Int_t i=0; i<kcoolingtubesupportmatrixnumber; i++)
	localcoolingtubesupportmatrix[i] = 
		new TGeoCombiTrans(*localcoolingtubesupportrans[i],
						   *localcoolingtubesupportrot[i]);
  fcoolingtubesupportmatrix[0] = new TGeoHMatrix(*localcoolingtubesupportmatrix[1]);
  fcoolingtubesupportmatrix[1] = new TGeoHMatrix((*localcoolingtubesupportmatrix[1])*
								(*localcoolingtubesupportmatrix[0]));
  /////////////////////////////////////////////////////////////
  // SSD Cooling Tube Transformations
  /////////////////////////////////////////////////////////////
  TGeoRotation* localcoolingtuberot = new TGeoRotation();	
  localcoolingtuberot->SetAngles(0.,90.,0.);
  TGeoTranslation** localcoolingtubetrans[4];
  TVector3** localcoolingtubevect[4];
  for(Int_t i=0; i<4; i++){
	localcoolingtubevect[i] = new TVector3*[2];
	localcoolingtubetrans[i] = new TGeoTranslation*[2];
	fcoolingtubematrix[i] = new TGeoHMatrix*[2];
  }
  localcoolingtubevect[0][0] = new TVector3(-0.5*(fgkCoolingTubeSeparation
						  -fgkCarbonFiberTriangleLength),
						  - (2.0*fgkSSDSensorLength-fgkSSDSensorOverlap)+
								 fgkSSDModuleStiffenerPosition[1]+fgkSSDStiffenerWidth
						  +	 0.5*fgkSSDFlexHoleLength+2.*fgkCarbonFiberJunctionWidth
						  -	 0.5*(fgkCarbonFiberLowerSupportWidth
						  +      fgkSSDSensorCenterSupportLength
						  -      fgkSSDSensorCenterSupportThickness[0])+
							 0.5*fgkSSDSensorLength-0.25*(fgkSSDSensorLength
						  -  2.0*fgkSSDModuleStiffenerPosition[1]
						  -	 2.0*fgkSSDCoolingBlockWidth-fgkCoolingTubeSupportWidth)
						  -  0.5*fgkCoolingTubeSupportWidth,
						  -  0.5*fgkCoolingTubeSupportHeight);	
  localcoolingtubevect[0][1] = new TVector3(localcoolingtubevect[0][0]->X(),
							localcoolingtubevect[0][0]->Y()+0.5*(fgkSSDSensorLength
						  -  2.0*fgkSSDModuleStiffenerPosition[1]
						  -	 2.0*fgkSSDCoolingBlockWidth-fgkCoolingTubeSupportWidth)
						  +  fgkCoolingTubeSupportWidth,
						  localcoolingtubevect[0][0]->Z());	
  localcoolingtubevect[1][0] = new TVector3(-localcoolingtubevect[0][0]->X()
							 +				 fgkCarbonFiberTriangleLength,
											 localcoolingtubevect[0][0]->Y(),
											 localcoolingtubevect[0][0]->Z());
  localcoolingtubevect[1][1] = new TVector3(-localcoolingtubevect[0][1]->X()
							 +				 fgkCarbonFiberTriangleLength,
											 localcoolingtubevect[0][1]->Y(),
											 localcoolingtubevect[0][1]->Z());
  localcoolingtubevect[2][0] = new TVector3(-0.5*(fgkCoolingTubeSeparation
						  -	fgkCarbonFiberTriangleLength),
						  - (2.0*fgkSSDSensorLength-fgkSSDSensorOverlap)+
								 fgkSSDModuleStiffenerPosition[1]+fgkSSDStiffenerWidth
						  +	 0.5*fgkSSDFlexHoleLength+2.*fgkCarbonFiberJunctionWidth
						  -	 0.5*(fgkCarbonFiberLowerSupportWidth
						  +      fgkSSDSensorCenterSupportLength
						  -      fgkSSDSensorCenterSupportThickness[0])
						  +  fgkSSDModuleStiffenerPosition[1]
						  -  0.5*(fgkSSDModuleStiffenerPosition[1]-fgkSSDSensorOverlap),
						  -  0.5*fgkCoolingTubeSupportHeight);	
  localcoolingtubevect[2][1] = new TVector3(-localcoolingtubevect[2][0]->X()
							 +				 fgkCarbonFiberTriangleLength,
											 localcoolingtubevect[2][0]->Y(),
											 localcoolingtubevect[2][0]->Z());	
  localcoolingtubevect[3][0] = new TVector3(-0.5*(fgkCoolingTubeSeparation
						  -	fgkCarbonFiberTriangleLength),
						  - (2.0*fgkSSDSensorLength-fgkSSDSensorOverlap)+
								 fgkSSDModuleStiffenerPosition[1]+fgkSSDStiffenerWidth
						  +	 0.5*fgkSSDFlexHoleLength+2.*fgkCarbonFiberJunctionWidth
						  -	 0.5*(fgkCarbonFiberLowerSupportWidth
						  +      fgkSSDSensorCenterSupportLength
						  -      fgkSSDSensorCenterSupportThickness[0])
						  +      fgkSSDSensorLength
						  -	 0.5*fgkSSDModuleStiffenerPosition[1],
						  -  0.5*fgkCoolingTubeSupportHeight);	
  localcoolingtubevect[3][1] = new TVector3(-localcoolingtubevect[3][0]->X()
						  + fgkCarbonFiberTriangleLength,
							localcoolingtubevect[3][0]->Y(),
						  - 0.5*fgkCoolingTubeSupportHeight);	
  for(Int_t i=0; i<4; i++) 
	for(Int_t j=0; j<2; j++){
		localcoolingtubetrans[i][j] = 
			new TGeoTranslation(localcoolingtubevect[i][j]->X(),
								localcoolingtubevect[i][j]->Y(),
								localcoolingtubevect[i][j]->Z());
		fcoolingtubematrix[i][j] = new TGeoHMatrix((*localcoolingtubetrans[i][j])
							  *					(*localcoolingtuberot));
	}
  /////////////////////////////////////////////////////////////
  // SSD Hybrid Components Transformations
  /////////////////////////////////////////////////////////////
  const Int_t khybridmatrixnumber = 3;
  TGeoTranslation* localhybridtrans[khybridmatrixnumber];
  localhybridtrans[0] = new TGeoTranslation(0.5*fgkSSDStiffenerLength,
                                            0.5*fgkSSDStiffenerWidth,
                                            0.5*fgkSSDStiffenerHeight);
  localhybridtrans[1] = new TGeoTranslation(fgkSSDModuleStiffenerPosition[0],
                                            fgkSSDModuleStiffenerPosition[1],0.0);

  localhybridtrans[2] = new TGeoTranslation(
                      -  0.5*(fgkSSDSensorWidth-fgkCarbonFiberTriangleLength),
                      -      (2.*fgkSSDSensorLength-fgkSSDSensorOverlap)+
                              fgkSSDModuleStiffenerPosition[1]+fgkSSDStiffenerWidth
                      +		0.5*fgkSSDFlexHoleLength+2.*fgkCarbonFiberJunctionWidth
                      -		0.5*(fgkCarbonFiberLowerSupportWidth+fgkSSDSensorCenterSupportLength
                      -       fgkSSDSensorCenterSupportThickness[0]),
                      -      (fgkSSDModuleCoolingBlockToSensor+0.5*fgkCoolingTubeSupportHeight
                      -       fgkSSDSensorHeight-fgkSSDChipCablesHeight[3]-fgkSSDChipHeight)); 
  fhybridmatrix = new TGeoHMatrix();
  for(Int_t i=0; i<khybridmatrixnumber; i++) fhybridmatrix->MultiplyLeft(localhybridtrans[i]);
  /////////////////////////////////////////////////////////////
  // SSD Cooling Block Transformations
  /////////////////////////////////////////////////////////////
  const Int_t kcoolingblockmatrixnumber = 4;    
  TGeoTranslation* localcoolingblocktrans[kcoolingblockmatrixnumber];
  localcoolingblocktrans[0] = new TGeoTranslation(-2.*(fgkCoolingTubeSupportRmax
                            -  fgkCoolingTubeSupportRmin),0.0,
                               0.5*(fgkSSDCoolingBlockHoleCenter+fgkCoolingTubeRmax));
  localcoolingblocktrans[1] = new TGeoTranslation(0.5*fgkSSDStiffenerLength
                            -  0.5*fgkSSDModuleSensorSupportDistance-fgkSSDCoolingBlockLength,
                               0.0,fgkSSDStiffenerHeight);
  localcoolingblocktrans[2] = new TGeoTranslation(*localhybridtrans[1]);
  localcoolingblocktrans[3] = new TGeoTranslation(*localhybridtrans[2]);
  fcoolingblocksystematrix = new TGeoHMatrix();
  for(Int_t i=0; i<kcoolingblockmatrixnumber; i++)
      fcoolingblocksystematrix->MultiplyLeft(localcoolingblocktrans[i]);
  /////////////////////////////////////////////////////////////
  // SSD Stiffener Flex Transformations
  /////////////////////////////////////////////////////////////
  const Int_t klocalflexmatrixnumber = 4;
  TGeoCombiTrans** localflexmatrix[fgkflexnumber];
  for(Int_t i=0; i<fgkflexnumber; i++)    
      localflexmatrix[i] = new TGeoCombiTrans*[klocalflexmatrixnumber];
  for(Int_t i=0; i<fgkflexnumber; i++)
      for(Int_t j =0; j<klocalflexmatrixnumber; j++) 
            localflexmatrix[i][j] = new TGeoCombiTrans();
  Double_t ssdstiffenerseparation = fgkSSDSensorLength
								  - 2.*fgkSSDModuleStiffenerPosition[1]
								  -    fgkSSDStiffenerWidth;
  localflexmatrix[0][0]->SetTranslation(-fgkSSDFlexLength[0]
                                        +0.5*fgkSSDStiffenerLength,
                                         0.5*fgkSSDStiffenerWidth,
                                        -0.5*fgkSSDStiffenerHeight
                                        -0.5*fgkSSDFlexHeight[0]);
  localflexmatrix[1][0]->SetTranslation(-(fgkSSDStiffenerLength-fgkSSDFlexLength[0])
                                        +0.5*fgkSSDStiffenerLength,ssdstiffenerseparation
                                        -0.5*fgkSSDStiffenerWidth,
                                        -0.5*fgkSSDStiffenerHeight
                                        -0.5*fgkSSDFlexHeight[0]);
  TGeoRotation* localflexrot = new TGeoRotation();
  localflexrot->SetAngles(180.,0.,0.);    
  localflexmatrix[1][0]->SetRotation(localflexrot);
  for(Int_t i=0; i<fgkflexnumber; i++)
      for(Int_t j =1; j<klocalflexmatrixnumber; j++) 
            localflexmatrix[i][j]->SetTranslation(*localhybridtrans[j-1]);
  for(Int_t i=0; i<fgkflexnumber; i++){
      fstiffenerflexmatrix[i] = new TGeoHMatrix();
      for(Int_t j =0; j<klocalflexmatrixnumber; j++)   
            fstiffenerflexmatrix[i]->MultiplyLeft(localflexmatrix[i][j]);
  }
  /////////////////////////////////////////////////////////////
  // SSD End Flex Transformations
  /////////////////////////////////////////////////////////////
  TGeoRotation* localendflexrot = new TGeoRotation();
  localendflexrot->SetAngles(0.0,90.0,0.0);
  TGeoCombiTrans* localendflexmatrix = new TGeoCombiTrans();
  Double_t ssdflexradiusmax = (fgkSSDFlexLength[3]-fgkSSDFlexLength[2])
                            /  TMath::Tan(fgkSSDFlexAngle*TMath::DegToRad());
  Double_t ssdflexboxlength = fgkSSDFlexFullLength-2.*fgkSSDFlexAngle
                            * TMath::DegToRad()*ssdflexradiusmax
					                       - fgkSSDFlexLength[2]-TMath::Pi()
					                       * fgkSSDStiffenerHeight-fgkSSDFlexLength[0];
  Double_t trans = ssdflexboxlength*CosD(2.*fgkSSDFlexAngle)
                            + (ssdflexradiusmax-fgkSSDStiffenerHeight)*SinD(2.*fgkSSDFlexAngle)
                            +      fgkSSDFlexLength[2];
  localendflexmatrix->SetTranslation(fgkSSDFlexLength[0]-trans,
                              0.5*fgkSSDFlexWidth[0],
                              2.*fgkSSDStiffenerHeight
                            + 0.5*fgkSSDFlexHeight[0]);      
  localendflexmatrix->SetRotation(localendflexrot);
  for(Int_t i=0; i<fgkflexnumber; i++) 
      fendflexmatrix[i] = new TGeoHMatrix((*fstiffenerflexmatrix[i])*(*localendflexmatrix));
  /////////////////////////////////////////////////////////////
  // End Ladder Carbon Fiber Junction
  /////////////////////////////////////////////////////////////
  TGeoCombiTrans** localendladdercarbonfiberjunctionmatrix[fgkendlabbercarbonfiberjunctionumber];
  TGeoRotation** localendladdercarbonfiberjunctionrot[fgkendlabbercarbonfiberjunctionumber];    
  TGeoTranslation** localendladdercarbonfiberjunctiontrans[fgkendlabbercarbonfiberjunctionumber];    
  for(Int_t i=0; i<fgkendlabbercarbonfiberjunctionumber; i++){
      localendladdercarbonfiberjunctionmatrix[i] 
            = new TGeoCombiTrans*[fgkendladdercabonfiberjunctionmatrixnumber];
      localendladdercarbonfiberjunctionrot[i] 
            = new TGeoRotation*[fgkendladdercabonfiberjunctionmatrixnumber];
      localendladdercarbonfiberjunctiontrans[i] 
            = new TGeoTranslation*[fgkendladdercabonfiberjunctionmatrixnumber];
      fendladdercarbonfiberjunctionmatrix[i]
            = new TGeoHMatrix*[fgkendladdercabonfiberjunctionmatrixnumber];
  }
  for(Int_t i=0; i<fgkendlabbercarbonfiberjunctionumber; i++)    
      for(Int_t j=0; j<fgkendladdercabonfiberjunctionmatrixnumber; j++){
            localendladdercarbonfiberjunctionrot[i][j] = new TGeoRotation();
            localendladdercarbonfiberjunctiontrans[i][j] = new TGeoTranslation();
      }
  for(Int_t i=0; i<fgkendlabbercarbonfiberjunctionumber; i++)     
      for(Int_t j=0; j<fgkendladdercabonfiberjunctionmatrixnumber; j++)
          localendladdercarbonfiberjunctionrot[i][j]->SetAngles(120.*j,0.,0.);
  for(Int_t i=0; i<fgkendlabbercarbonfiberjunctionumber; i++){
      localendladdercarbonfiberjunctiontrans[i][1]->SetTranslation(fgkCarbonFiberTriangleLength,
                              0.0,0.0);
      localendladdercarbonfiberjunctiontrans[i][2]->SetTranslation(fgkCarbonFiberTriangleLength
		*                     CosD(fgkCarbonFiberTriangleAngle),fgkCarbonFiberTriangleLength
		*                     SinD(fgkCarbonFiberTriangleAngle),
                        0.5*(fgkEndLadderCarbonFiberLowerJunctionLength[i]
  -                          fgkEndLadderCarbonFiberUpperJunctionLength[i]));
  }
  TGeoCombiTrans* localendladdercarbonfiberjunctionglobalmatrix[fgkendlabbercarbonfiberjunctionumber];
  TGeoRotation* localendladdercarbonfiberjunctionglobalrot[fgkendlabbercarbonfiberjunctionumber];
  TGeoTranslation* localendladdercarbonfiberjunctionglobaltrans[fgkendlabbercarbonfiberjunctionumber];
  for(Int_t i=0; i<fgkendlabbercarbonfiberjunctionumber; i++){
      localendladdercarbonfiberjunctionglobalrot[i] = new TGeoRotation();
      localendladdercarbonfiberjunctionglobaltrans[i] = new TGeoTranslation();
      localendladdercarbonfiberjunctionglobalrot[i]->SetAngles(0.0,90.0,0.0);
      localendladdercarbonfiberjunctionglobaltrans[i]->SetTranslation(0.0,
            0.5*fgkEndLadderCarbonFiberLowerJunctionLength[i],0.0);
      localendladdercarbonfiberjunctionglobalmatrix[i] = 
            new TGeoCombiTrans(*localendladdercarbonfiberjunctionglobaltrans[i],
                               *localendladdercarbonfiberjunctionglobalrot[i]);
  }
  for(Int_t i=0; i<fgkendlabbercarbonfiberjunctionumber; i++)    
      for(Int_t j=0; j<fgkendladdercabonfiberjunctionmatrixnumber; j++){
            localendladdercarbonfiberjunctionmatrix[i][j] = 
                  new TGeoCombiTrans(*localendladdercarbonfiberjunctiontrans[i][j],
                                     *localendladdercarbonfiberjunctionrot[i][j]);
           fendladdercarbonfiberjunctionmatrix[i][j] =
            new TGeoHMatrix((*localendladdercarbonfiberjunctionglobalmatrix[i])
            *               (*localendladdercarbonfiberjunctionmatrix[i][j])); 
      }  
  /////////////////////////////////////////////////////////////
  // End Ladder Carbon Fiber Support
  /////////////////////////////////////////////////////////////
  TGeoTranslation* localendladdercarbonfibertrans[fgkendladdercarbonfibermatrixnumber];
  for(Int_t i=0; i<fgkendladdercarbonfibermatrixnumber; i++){
      localendladdercarbonfibertrans[i] = new TGeoTranslation();
      localendladdercarbonfibertrans[i]->SetTranslation(0.0,
            i==0 ? 0.0 :fgkCarbonFiberLowerSupportWidth,0.0);
      fendladdercarbonfibermatrix[i] = new TGeoHMatrix*[fgkcarbonfibersupportnumber];
  }
  for(Int_t i=0; i<fgkendladdercarbonfibermatrixnumber; i++)
      for(Int_t j=0; j<fgkcarbonfibersupportnumber; j++)
            fendladdercarbonfibermatrix[i][j] = 
            new TGeoHMatrix((*localendladdercarbonfibertrans[i])
            *(*fcarbonfibersupportmatrix[j]));
  /////////////////////////////////////////////////////////////
  // End Ladder SSD Mounting Block
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i<fgkendladdermountingblocknumber; i++)
      fendladdermountingblocktrans[i] = 
            new TGeoTranslation(-  0.25*(fgkSSDMountingBlockLength[0]
                                +	 fgkSSDMountingBlockLength[1])
                                +  0.5*fgkCarbonFiberTriangleLength,
                                fgkEndLadderMountingBlockPosition[i],
                                -  fgkSSDMountingBlockHeight[1]
                                +  0.5*fgkSSDMountingBlockHeight[0]);
  /////////////////////////////////////////////////////////////
  // End Ladder Carbon Fiber Lower Support
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i<fgkendladderlowersuppnumber; i++)
      fendladderlowersupptrans[i] = 
            new TGeoTranslation(0.0,(1-i)*(fgkEndLadderMountingBlockPosition[i]
                        +  0.5*fgkSSDMountingBlockWidth),
                        -  0.5*fgkCarbonFiberLowerSupportHeight);
  fendladderlowersupptrans[2] = new TGeoTranslation(0.0,
									 fgkCarbonFiberLowerSupportVolumePosition[1]
								+    fgkCarbonFiberLowerSupportVolumePosition[0],
									 0.0);
  fendladderlowersupptrans[2]->Add(fendladderlowersupptrans[1]);
 /////////////////////////////////////////////////////////////
  // Matrix for positioning Ladder into mother volume
  /////////////////////////////////////////////////////////////
  TGeoHMatrix** ladderglobalmatrix[fgkladdernumber];
  for(Int_t i=0; i<fgkladdernumber; i++) 
	ladderglobalmatrix[i] = new TGeoHMatrix*[fgkladdernumber];
  TGeoRotation* localladdermotherrot = new TGeoRotation();
  localladdermotherrot->SetAngles(0.,90.,0.);  
  TGeoTranslation* localladdermothertrans[fgkladdernumber];
  TGeoCombiTrans* localladdermothercombitrans[fgkladdernumber];
  for(Int_t i=0; i<fgkladdernumber; i++){
	localladdermothertrans[i] = new TGeoTranslation(0.,
							  - fgkEndLadderCarbonFiberLowerJunctionLength[1]
							  + fgkEndLadderCarbonFiberLowerJunctionLength[0]
							  + (i==0?fgkSSDLay5SensorsNumber:fgkSSDLay6SensorsNumber)
							  * fgkCarbonFiberJunctionWidth,0.);
	localladdermothercombitrans[i] = new TGeoCombiTrans(*localladdermothertrans[i],
														*localladdermotherrot);
	ladderglobalmatrix[0][i] = new TGeoHMatrix(*localladdermothercombitrans[i]);
	ladderglobalmatrix[1][i] = new TGeoHMatrix(ladderglobalmatrix[0][i]->Inverse());
  }
  /////////////////////////////////////////////////////////////
  // Ladder Cables Matrices
  /////////////////////////////////////////////////////////////
  Double_t ssdflexradius = fgkSSDStiffenerHeight+2*fgkSSDFlexHeight[0]
					     + fgkSSDFlexHeight[1];  
  Double_t ssdladdercabletransx[3];
  ssdladdercabletransx[0] = (ssdflexradiusmax-fgkSSDFlexHeight[1]-ssdflexradius)
						  *   SinD(2.*fgkSSDFlexAngle)
						  *	  CosD(2.*fgkSSDFlexAngle);
  ssdladdercabletransx[1] = ((ssdflexradiusmax-fgkSSDFlexHeight[1]-ssdflexradius)
						  -     ssdladdercabletransx[0]
						  /     SinD(2.*fgkSSDFlexAngle))
						  *     CosD(fgkSSDFlexAngle);						
  ssdladdercabletransx[2] = (fgkSSDFlexFullLength-2.*fgkSSDFlexAngle
						  *	  TMath::DegToRad()*ssdflexradiusmax
						  -     fgkSSDFlexLength[2]-TMath::Pi()
						  *	  fgkSSDStiffenerHeight-fgkSSDFlexLength[0]
						  -	  fgkSSDLadderCableWidth)
						  *	  CosD(2.*fgkSSDFlexAngle);
  Double_t ssdladdercabletransz[3] = {ssdladdercabletransx[0]
						  *	TanD(2.*fgkSSDFlexAngle),
							ssdladdercabletransx[1]
						  *	TanD(fgkSSDFlexAngle),
							ssdladdercabletransx[2]
						  *	TanD(2.*fgkSSDFlexAngle)};	
  TGeoRotation* localladdercablerot[3];	
  for(Int_t i=0; i<3; i++) localladdercablerot[i] = new TGeoRotation();
  localladdercablerot[0]->SetAngles(90.,0.,0.);
  localladdercablerot[1]->SetAngles(90.,60.,-90.);
  localladdercablerot[2]->SetRotation((*localladdercablerot[1])
						 *			  (*localladdercablerot[0]));
  ////////////////////////////////////////////
  // LocalLadderCableCombiTransMatrix
  ////////////////////////////////////////////
  const Int_t klocalladdersidecablesnumber = 2;
  const Int_t klocalladdercombitransnumber = 5;
  TGeoCombiTrans** localladdercablecombitransmatrix[klocalladdersidecablesnumber];
  for(Int_t i=0; i<klocalladdersidecablesnumber; i++) 
	 localladdercablecombitransmatrix[i] = 
							   new TGeoCombiTrans*[klocalladdercombitransnumber];
  ///////////////////////////////////////////
  // Left Side Ladder Cables Transformations
  ///////////////////////////////////////////
  localladdercablecombitransmatrix[0][0]  =
						new TGeoCombiTrans(-0.5*fgkCarbonFiberTriangleLength,
						0.,0.,NULL);
  localladdercablecombitransmatrix[0][1] = 
	new TGeoCombiTrans(-0.5*(fgkSSDSensorWidth-fgkCarbonFiberTriangleLength),
					   - (2.*fgkSSDSensorLength-fgkSSDSensorOverlap)+
						 fgkSSDModuleStiffenerPosition[1]+fgkSSDStiffenerWidth
					   + 0.5*fgkSSDFlexHoleLength+2.*fgkCarbonFiberJunctionWidth
					   - 0.5*(fgkCarbonFiberLowerSupportWidth
					   + fgkSSDSensorCenterSupportLength
					   - fgkSSDSensorCenterSupportThickness[0]),
					   - (fgkSSDModuleCoolingBlockToSensor
					   + 0.5*fgkCoolingTubeSupportHeight
					   - fgkSSDSensorHeight-fgkSSDChipCablesHeight[3]
					   - fgkSSDChipHeight),NULL);
  localladdercablecombitransmatrix[0][2] = 
						new TGeoCombiTrans(fgkSSDModuleStiffenerPosition[0],
										   fgkSSDModuleStiffenerPosition[1],0.,0);
  localladdercablecombitransmatrix[0][3] = new TGeoCombiTrans(
					0.5*(fgkSSDStiffenerLength+fgkSSDChipNumber*fgkSSDChipLength
				   +(fgkSSDChipNumber-1)*fgkSSDChipSeparationLength),
				   fgkSSDSensorLength-2.*fgkSSDModuleStiffenerPosition[1]
					 - fgkSSDStiffenerWidth,- 0.5*fgkSSDChipHeight,
												new TGeoRotation("",180.,0.,0.));
  localladdercablecombitransmatrix[0][4] = 
						new TGeoCombiTrans(-ssdladdercabletransx[0]
						-     ssdladdercabletransx[1]-ssdladdercabletransx[2]
						+     fgkSSDFlexLength[0]-fgkSSDFlexLength[2],
							  0.,
							  0.5*fgkSSDFlexHeight[0]+2.*(fgkSSDFlexHeight[0]
						+	  fgkSSDFlexHeight[1])+2.*fgkSSDStiffenerHeight
						+     ssdladdercabletransz[0]-ssdladdercabletransz[1]
						+	  ssdladdercabletransz[2],localladdercablerot[2]);
  ///////////////////////////////////////////
  // Rigth Side Ladder Cables Transformations
  ///////////////////////////////////////////
  TGeoCombiTrans* localladdercablessdmodulematrix = 
	new TGeoCombiTrans(0.5*(fgkSSDStiffenerLength-fgkSSDChipNumber*fgkSSDChipLength
								  - (fgkSSDChipNumber-1)*fgkSSDChipSeparationLength),
									 fgkSSDStiffenerWidth,
								  - 0.5*fgkSSDFlexHeight[0],NULL);
  for(Int_t i=0; i<klocalladdercombitransnumber; i++)
   localladdercablecombitransmatrix[1][i] = 
			(i!=3 ? new TGeoCombiTrans(*localladdercablecombitransmatrix[0][i]):
					new TGeoCombiTrans(*localladdercablessdmodulematrix)); 	
  ///////////////////////////////////////////
  // Setting LadderCableHMatrix
  ///////////////////////////////////////////
  Int_t beamaxistrans[2][3];
  beamaxistrans[0][0] = fgkSSDLay5SensorsNumber/2; 
  beamaxistrans[0][1] = beamaxistrans[0][0]+1;
  beamaxistrans[0][2] = beamaxistrans[0][0]-1;
  beamaxistrans[1][0] = (fgkSSDLay6SensorsNumber-1)/2;
  beamaxistrans[1][1] = beamaxistrans[1][0]+1;
  beamaxistrans[1][2] = beamaxistrans[1][0];
  TGeoHMatrix** localladdercablehmatrix[fgkladdercablesnumber];
  TGeoRotation* laddercablerot = new TGeoRotation();
  TGeoTranslation* laddercabletrans = new TGeoTranslation();
  TGeoCombiTrans* laddercablecombitrans = new TGeoCombiTrans();
  Double_t* laddercabletransvector;	
  for(Int_t i=0; i<fgkladdercablesnumber; i++){ 
	localladdercablehmatrix[i] = new TGeoHMatrix*[klocalladdersidecablesnumber];
	fladdercablematrix[i] = new TGeoHMatrix*[fgkladdercablematrixnumber];
  }
  for(Int_t i=0; i<fgkladdercablesnumber; i++){
	for(Int_t j=0; j<klocalladdersidecablesnumber; j++){
		localladdercablehmatrix[i][j] = new TGeoHMatrix();
		for(Int_t k=0; k<klocalladdercombitransnumber; k++){
			localladdercablehmatrix[i][j]->MultiplyLeft(
			localladdercablecombitransmatrix[j][klocalladdercombitransnumber-k-1]);
        }
		laddercablerot->SetMatrix(localladdercablehmatrix[i][j]->GetRotationMatrix());
		laddercabletransvector = localladdercablehmatrix[i][j]->GetTranslation();
		laddercabletrans->SetTranslation(laddercabletransvector[0],
									 laddercabletransvector[1]
					+                (j==0 ? beamaxistrans[i][0] : 0.)
					*				 fgkCarbonFiberJunctionWidth,
									 laddercabletransvector[2]);
		laddercablecombitrans->SetRotation(*laddercablerot);
		laddercablecombitrans->SetTranslation(*laddercabletrans);	
		fladdercablematrix[i][j] = new TGeoHMatrix(*laddercablecombitrans);
	}
    fladdercablematrix[i][2] = 
					AddTranslationToHMatrix(fladdercablematrix[i][1],0.,
					beamaxistrans[i][1]*fgkCarbonFiberJunctionWidth,0.);
	fladdercablematrix[i][3] = 
					AddTranslationToHMatrix(fladdercablematrix[i][0],0.,
					beamaxistrans[i][2]*fgkCarbonFiberJunctionWidth,0.);
  }
  for(Int_t i=0; i<fgkladdercablesnumber; i++)
	for(Int_t j=0; j<klocalladdercombitransnumber-1; j++)
		fladdercablematrix[i][j]->MultiplyLeft(ladderglobalmatrix[1][i]);
  ///////////////////////////////////////////
  // Setting Ladder HMatrix
  ///////////////////////////////////////////
  Int_t ssdlaysensorsnumber[fgkladdernumber] = {fgkSSDLay5SensorsNumber,
												fgkSSDLay6SensorsNumber};
  for(Int_t i=0; i<fgkladdernumber; i++){
	fladdermatrix[i] = new TGeoHMatrix*[ssdlaysensorsnumber[i]];
	for(Int_t j=0; j<ssdlaysensorsnumber[i]; j++){
		fladdermatrix[i][j] = new TGeoHMatrix();
		fladdermatrix[i][j]->SetDx(-0.5*fgkCarbonFiberTriangleLength);
		fladdermatrix[i][j]->SetDy(fgkCarbonFiberJunctionWidth*j);
		fladdermatrix[i][j]->MultiplyLeft(ladderglobalmatrix[1][i]);
	}
  }
  ///////////////////////////////////////////
  // Setting SSD Sensor Matrix 
  ///////////////////////////////////////////
  TGeoCombiTrans* localssdsensorcombitrans[2];
  TGeoRotation* localssdsensorrot = new TGeoRotation();	
  localssdsensorrot->SetAngles(0.,90.,0.);	
  TGeoTranslation* localssdsensortrans[2];
  for(Int_t i=0; i<2; i++) localssdsensortrans[i] = new TGeoTranslation();
  localssdsensortrans[0]->SetTranslation(0.5*fgkCarbonFiberTriangleLength,
					  -		(2.*fgkSSDSensorLength-fgkSSDSensorOverlap)+
                              fgkSSDModuleStiffenerPosition[1]+fgkSSDStiffenerWidth
                      +		0.5*fgkSSDFlexHoleLength+2.*fgkCarbonFiberJunctionWidth
                      -		0.5*(fgkCarbonFiberLowerSupportWidth+fgkSSDSensorCenterSupportLength
                      -       fgkSSDSensorCenterSupportThickness[0])+0.5*fgkSSDSensorLength,
							0.5*fgkSSDSensorHeight-0.5*fgkCoolingTubeSupportHeight
					  -		fgkSSDModuleCoolingBlockToSensor+(fgkSSDSensorSideSupportHeight[1]
					  -		fgkSSDSensorSideSupportHeight[0]));
  localssdsensortrans[1]->SetTranslation(0.5*fgkCarbonFiberTriangleLength,
					  -	   (2.*fgkSSDSensorLength-fgkSSDSensorOverlap)+
                              fgkSSDModuleStiffenerPosition[1]+fgkSSDStiffenerWidth
                      +		0.5*fgkSSDFlexHoleLength+2.*fgkCarbonFiberJunctionWidth
                      -		0.5*(fgkCarbonFiberLowerSupportWidth+fgkSSDSensorCenterSupportLength
                      -       fgkSSDSensorCenterSupportThickness[0])+0.5*fgkSSDSensorLength,
						    0.5*fgkSSDSensorHeight-0.5*fgkCoolingTubeSupportHeight
							-fgkSSDModuleCoolingBlockToSensor);
  for(Int_t i=0; i<2; i++) 
	localssdsensorcombitrans[i] = new TGeoCombiTrans(*localssdsensortrans[i],
													 *localssdsensorrot);	
    for(Int_t i=0; i<fgkladdernumber; i++){
	fssdsensormatrix[i] = new TGeoHMatrix*[ssdlaysensorsnumber[i]];
	for(Int_t j=0; j<ssdlaysensorsnumber[i]; j++){
		switch(i){
			case 0: //Ladder of Layer5  
			fssdsensormatrix[i][j] = new TGeoHMatrix((*fladdermatrix[i][j])
								   * ((j%2==0 ? *localssdsensorcombitrans[0] :
												*localssdsensorcombitrans[1])));
			break;
			case 1: //Ladder of Layer6 
			fssdsensormatrix[i][j] = new TGeoHMatrix((*fladdermatrix[i][j])
								   * ((j%2==0 ? *localssdsensorcombitrans[1] :
												*localssdsensorcombitrans[0])));
		break;
		}
	  }
  }	
  //////////////////////////
  // Setting SSD End Ladder  
  //////////////////////////
  for(Int_t i=0; i<2; i++) fendladdersegmentmatrix[i] = new TGeoHMatrix*[2];
  for(Int_t i=0; i<2; i++){
	fendladdersegmentmatrix[0][i] = new TGeoHMatrix();
	fendladdersegmentmatrix[0][i]->SetDx(-0.5*fgkCarbonFiberTriangleLength);
	fendladdersegmentmatrix[0][i]->SetDy(fgkCarbonFiberJunctionWidth*ssdlaysensorsnumber[i]);
	fendladdersegmentmatrix[0][i]->MultiplyLeft(ladderglobalmatrix[1][i]);
	fendladdersegmentmatrix[1][i] = new TGeoHMatrix();
	fendladdersegmentmatrix[1][i]->SetDx(-0.5*fgkCarbonFiberTriangleLength);
	fendladdersegmentmatrix[1][i]->RotateZ(180.0);
	fendladdersegmentmatrix[1][i]->MultiplyLeft(ladderglobalmatrix[1][i]);
   }
  /////////////////////////////////////////////////////
  // Setting the CombiTransformation to pass ITS center 
  /////////////////////////////////////////////////////
  Double_t itscentertransz[fgklayernumber];
  itscentertransz[0] = fgkSSDLay5LadderLength
					 - fgkLay5CenterITSPosition;
  itscentertransz[1] = fgkSSDLay6LadderLength
					 - fgkLay6CenterITSPosition;
  Double_t itssensortransy = fgkSSDModuleCoolingBlockToSensor
						   + 0.5*fgkCoolingTubeSupportHeight;
  TGeoRotation* itscenterrot[3];
  for(Int_t i=0; i<fgklayernumber; i++) itscenterrot[i] = new TGeoRotation();
  itscenterrot[0]->SetAngles(90.,180.,-90.);
  itscenterrot[1]->SetAngles(0.,90.,0.);
  itscenterrot[2] = new TGeoRotation((*itscenterrot[1])*(*itscenterrot[0]));
  TGeoCombiTrans* itscentercombitrans[fgklayernumber];
  for(Int_t i=0; i<fgklayernumber; i++) 
	itscentercombitrans[i] = new TGeoCombiTrans(0.,
							 itssensortransy,
							 fgkEndLadderCarbonFiberLowerJunctionLength[1]
						   - itscentertransz[i],itscenterrot[2]);
  TGeoRotation** locallayerrot[fgklayernumber];
  TGeoTranslation** locallayertrans[fgklayernumber];	
  TGeoCombiTrans** locallayercombitrans[fgklayernumber];
  TGeoTranslation* localbeamaxistrans[fgklayernumber];
  localbeamaxistrans[0] = new TGeoTranslation(0.,0.,0.5*fgkSSDLay5LadderLength
					 - fgkLay5CenterITSPosition);
  localbeamaxistrans[1] = new TGeoTranslation(0.,0.,0.5*fgkSSDLay6LadderLength
					 - fgkLay6CenterITSPosition);
  const Int_t kssdlayladdernumber[fgklayernumber] = 
			{fgkSSDLay5LadderNumber,fgkSSDLay6LadderNumber};
  for(Int_t i=0; i<fgklayernumber; i++){
    locallayerrot[i] = new TGeoRotation*[kssdlayladdernumber[i]];
    locallayertrans[i] = new TGeoTranslation*[kssdlayladdernumber[i]];
	locallayercombitrans[i] = new TGeoCombiTrans*[kssdlayladdernumber[i]];
	flayermatrix[i] = new TGeoHMatrix*[kssdlayladdernumber[i]];
  }
  Double_t layerladderangleposition[fgklayernumber] = 
		{360./fgkSSDLay5LadderNumber,360./fgkSSDLay6LadderNumber};
  Double_t layerradius = 0.;
  for(Int_t i=0; i<fgklayernumber; i++){	
	for(Int_t j=0; j<kssdlayladdernumber[i]; j++){
		switch(i){
			case 0: //Ladder of Layer5  
			layerradius = (j%2==0 ? fgkSSDLay5RadiusMin: fgkSSDLay5RadiusMax);
			break;
			case 1: //Ladder of Layer6 
			layerradius = (j%2==0 ? fgkSSDLay6RadiusMin: fgkSSDLay6RadiusMax);
		break;
		}
		locallayerrot[i][j] = new TGeoRotation();
		locallayertrans[i][j] = new TGeoTranslation();
		locallayerrot[i][j]->SetAngles(j*layerladderangleposition[i],0.,0.);
		locallayertrans[i][j]->SetTranslation(layerradius 
							  *	CosD(90.0+j*layerladderangleposition[i]),
							    layerradius 
							  * SinD(90.0+j*layerladderangleposition[i]),0.);
		locallayercombitrans[i][j] = new TGeoCombiTrans(*locallayertrans[i][j],
									 *locallayerrot[i][j]);
		flayermatrix[i][j] = new TGeoHMatrix((*locallayercombitrans[i][j])*(*itscentercombitrans[i]));
		flayermatrix[i][j]->Multiply(ladderglobalmatrix[0][i]);
		flayermatrix[i][j]->MultiplyLeft(localbeamaxistrans[i]);
	}
  }
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i< kcarbonfibersupportmatrixnumber[1]; i++){
	delete carbonfiberot[i];
	delete localcarbonfibersupportmatrix[i];
  }
  for(Int_t i=0; i< fgkcarbonfiberjunctionumber; i++){
     for(Int_t j=0; j< kcarbonfiberjunctionmatrixnumber; j++){
       delete localcarbonfiberjunctionmatrix[i][j];
	   delete localcarbonfiberjunctionrot[i][j];
	   delete localcarbonfiberjunctiontrans[i][j];
	   }
       delete [] localcarbonfiberjunctionmatrix[i];
       delete [] localcarbonfiberjunctionrot[i];
       delete [] localcarbonfiberjunctiontrans[i];
  }
  for(Int_t i=0; i<fgkcarbonfiberlowersupportnumber; i++) 
	   delete localcarbonfiberlowersupportrans[i];
  for(Int_t i=0; i< fgkssdsensorsupportnumber; i++){
     for(Int_t j=0; j< kssdsensorsupportmatrixnumber; j++){
       delete localssdsensorsupportmatrix[i][j];
	   delete localssdsensorsupportrot[i][j];
	   delete localssdsensorsupportrans[i][j];
	   }
       delete [] localssdsensorsupportmatrix[i];
       delete [] localssdsensorsupportrot[i];
       delete [] localssdsensorsupportrans[i];
  }
  for(Int_t i=0; i<kcoolingtubesupportmatrixnumber; i++){
	delete localcoolingtubesupportmatrix[i];
	delete localcoolingtubesupportrot[i];
	delete localcoolingtubesupportrans[i];
  }
  for(Int_t i=0; i<4; i++){
	for(Int_t j=0; j<2; j++){
		delete localcoolingtubevect[i][j];
		delete localcoolingtubetrans[i][j];
	}
	delete [] localcoolingtubevect[i];
	delete [] localcoolingtubetrans[i];
  }
 for(Int_t i=0; i<khybridmatrixnumber; i++) delete localhybridtrans[i];
 for(Int_t i=0; i<kcoolingblockmatrixnumber; i++) delete localcoolingblocktrans[i];
 for(Int_t i=0; i<fgkflexnumber; i++){
      for(Int_t j=1; j<klocalflexmatrixnumber; j++) 
            delete localflexmatrix[i][j];
      delete [] localflexmatrix[i];
 }
 delete localflexrot;
 delete localendflexrot;
 delete localendflexmatrix;
 for(Int_t i=0; i<fgkladdernumber; i++){ 
	delete localladdermothertrans[i];
	delete localladdermothercombitrans[i];
  }
 delete localladdermotherrot;
 for(Int_t i=0; i<fgkendlabbercarbonfiberjunctionumber; i++){    
      for(Int_t j=0; j<fgkendladdercabonfiberjunctionmatrixnumber; j++){
            delete localendladdercarbonfiberjunctionmatrix[i][j];
            delete localendladdercarbonfiberjunctionrot[i][j];
            delete localendladdercarbonfiberjunctiontrans[i][j];
      }
      delete [] localendladdercarbonfiberjunctionmatrix[i];
      delete [] localendladdercarbonfiberjunctionrot[i];
      delete [] localendladdercarbonfiberjunctiontrans[i];
      delete localendladdercarbonfiberjunctionglobalrot[i];
      delete localendladdercarbonfiberjunctionglobaltrans[i];
      delete localendladdercarbonfiberjunctionglobalmatrix[i];
 }
 for(Int_t i=0; i<fgkendladdercarbonfibermatrixnumber; i++)
      delete localendladdercarbonfibertrans[i];
  for(Int_t i=0; i<3; i++) delete localladdercablerot[i];
  for(Int_t i=0; i<klocalladdersidecablesnumber; i++){
	for(Int_t j=0; j<klocalladdercombitransnumber; j++)
		delete localladdercablecombitransmatrix[i][j];
		delete []localladdercablecombitransmatrix[i];
  }
  for(Int_t i=0; i<fgkladdercablesnumber; i++){
	for(Int_t j=0; j<klocalladdersidecablesnumber; j++)
		delete localladdercablehmatrix[i][j];
	delete []localladdercablehmatrix[i];
  }
  delete laddercablerot;
  delete laddercabletrans;
  delete laddercablecombitrans;
  delete localladdercablessdmodulematrix;
  delete localssdsensorrot;	
  for(Int_t i=0; i<2; i++){
	delete localssdsensortrans[i];
	delete localssdsensorcombitrans[i];
  }
  for(Int_t i=0; i<fgklayernumber; i++){
	for(Int_t j=0; j<kssdlayladdernumber[i]; j++){
		delete locallayerrot[i][j];
		delete locallayertrans[i][j];
		delete locallayercombitrans[i][j];
    }
	delete [] locallayerrot[i];
	delete [] locallayertrans[i];
	delete [] locallayercombitrans[i];
	delete localbeamaxistrans[i];
  }
  for(Int_t i=0; i<3; i++) delete itscenterrot[i];
  for(Int_t i=0; i<fgkladdernumber; i++){
	for(Int_t j=0; j<fgkladdernumber; j++)
		delete ladderglobalmatrix[i][j];
	delete [] ladderglobalmatrix[i];
  }
  /////////////////////////////////////////////////////////////
  fTransformationMatrices = kTRUE;	
}
///////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::CreateBasicObjects(){
  /////////////////////////////////////////////////////////////  
  // Method generating the Objects of SSD Geometry    
  /////////////////////////////////////////////////////////////
  // SSD Sensor
  ///////////////////////////////////
  SetSSDSensor();
  /////////////////////////////////////////////////////////////  
  // Carbon Fiber Support    
  /////////////////////////////////////////////////////////////  
  TList* carbonfibersupportlist = GetCarbonFiberSupportList();  
  for(Int_t i=0; i<fgkcarbonfibersupportnumber; i++) 
      fcarbonfibersupport[i] = (TGeoVolume*)carbonfibersupportlist->At(i);
  /////////////////////////////////////////////////////////////
  // Carbon Fiber Junction 
  /////////////////////////////////////////////////////////////
  fcarbonfiberjunction = GetCarbonFiberJunction(fgkCarbonFiberJunctionWidth);
  /////////////////////////////////////////////////////////////
  // Carbon Fiber Lower Support
  /////////////////////////////////////////////////////////////
  TList* carbonfiberlowersupportlist = GetCarbonFiberLowerSupportList();
  for(Int_t i=0; i<fgkcarbonfiberlowersupportnumber; i++)
 	fcarbonfiberlowersupport[i] = (TGeoVolume*)carbonfiberlowersupportlist->At(i);
  /////////////////////////////
  // SSD Sensor Support
  /////////////////////////////
  for(Int_t i=0; i<fgkvolumekind; i++) fssdsensorsupport[i] = 
										new TGeoVolume*[fgkssdsensorsupportnumber]; 
  Double_t sidesupporthickness[2] = {fgkSSDSensorSideSupportThickness[0],
									 fgkSSDSensorSideSupportThickness[1]};
  for(Int_t i=0; i<fgkssdsensorsupportnumber-1; i++){
	fssdsensorsupport[0][i] = GetSSDSensorSupport(fgkSSDSensorSideSupportLength,
											   fgkSSDSensorSideSupportHeight[i],
											   fgkSSDSensorSideSupportWidth,
											   sidesupporthickness);  
	fssdsensorsupport[1][i] = GetSSDSensorSupport(fgkSSDSensorCenterSupportLength,
											   fgkSSDSensorCenterSupportHeight[i],
											   fgkSSDSensorCenterSupportWidth,
											   sidesupporthickness);
  }
  /////////////////////////////////////////////////////////////
  // SSD Cooling Tube Support
  /////////////////////////////////////////////////////////////
  Int_t edgesnumber = 16;
  fcoolingtubesupport = GetCoolingTubeSupport(edgesnumber);	  
  /////////////////////////////////////////////////////////////
  // SSD Hybrid
  /////////////////////////////////////////////////////////////
  TList* ssdhybridcomponentslist = GetSSDHybridParts();
  for(Int_t i=0; i<fgkhybridcompnumber; i++) 
	fssdhybridcomponent[i] = (TGeoVolume*)ssdhybridcomponentslist->At(i);
  /////////////////////////////////////////////////////////////
  // SSD Cooling Block System
  /////////////////////////////////////////////////////////////
  fssdcoolingblocksystem = GetCoolingBlockSystem();
   /////////////////////////////////////////////////////////////
  // SSD Cooling Tube
  /////////////////////////////////////////////////////////////
  TList* coolingtubelist = GetCoolingTubeList();	
  for(Int_t i=0; i<fgkcoolingtubenumber; i++)	
	fcoolingtube[i] = (TGeoVolume*)coolingtubelist->At(i);
  /////////////////////////////////////////////////////////////
  // SSD Flex  
  /////////////////////////////////////////////////////////////
  fssdstiffenerflex = GetSSDStiffenerFlex();
  fssdendflex = GetSSDEndFlex();
  ///////////////////////////////////
  // End Ladder Carbon Fiber Junction
  ///////////////////////////////////
  for(Int_t i=0; i<fgkendlabbercarbonfiberjunctionumber; i++) 
						   fendladdercarbonfiberjunction[i] = 
						   new TGeoVolume*[fgkendlabbercarbonfiberjunctionumber];
  for(Int_t i=0; i<fgkendlabbercarbonfiberjunctionumber; i++){
    fendladdercarbonfiberjunction[i][0] = 
		  GetCarbonFiberJunction(fgkEndLadderCarbonFiberLowerJunctionLength[i]);
    fendladdercarbonfiberjunction[i][1] = 
		  GetCarbonFiberJunction(fgkEndLadderCarbonFiberUpperJunctionLength[i]);
  }
  ///////////////////////////////////
  // End Ladder Mounting Block
  ///////////////////////////////////
  fendladdermountingblock = GetSSDMountingBlock();
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  delete carbonfibersupportlist;
  delete carbonfiberlowersupportlist;
  delete ssdhybridcomponentslist;
  /////////////////////////////////////////////////////////////
  fBasicObjects = kTRUE;
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetSSDSensor(){
  ////////////////////////////////////////////////////////////////
  // Method generating SSD Sensors: it sets the private variables
  // fSSDSensor5, fSSDSensor6  
  ////////////////////////////////////////////////////////////////
  Double_t ssdsensitivelength = fgkSSDSensorLength-2.*fgkSSDSensorInsensitiveLength;
  Double_t ssdsensitivewidth  = fgkSSDSensorWidth-2.*fgkSSDSensorInsensitiveWidth;
  TGeoBBox* ssdsensorsensitiveshape = new TGeoBBox("SSDSensorSensitiveShape",
                                                0.5*ssdsensitivewidth,
                                                0.5*fgkSSDSensorHeight,
                                                0.5*ssdsensitivelength);
  TGeoVolume* ssdsensorsensitiveLay5 = 
	new TGeoVolume(fgSDDsensitiveVolName5,ssdsensorsensitiveshape,fSSDSensorMedium);
  TGeoVolume* ssdsensorsensitiveLay6 = 
	new TGeoVolume(fgSDDsensitiveVolName6,ssdsensorsensitiveshape,fSSDSensorMedium);
  ssdsensorsensitiveLay5->SetLineColor(fColorSilicon);
  ssdsensorsensitiveLay6->SetLineColor(fColorSilicon);
  TGeoBBox* ssdsensorinsensitiveshape[2];
  ssdsensorinsensitiveshape[0] = new TGeoBBox("SSDSensorInsensitiveShape1",
                                                0.5*fgkSSDSensorInsensitiveWidth,
                                                0.5*fgkSSDSensorHeight,
                                                0.5*fgkSSDSensorLength);
  ssdsensorinsensitiveshape[1] = new TGeoBBox("SSDSensorInsensitiveShape2",
                                                0.5*ssdsensitivewidth,
                                                0.5*fgkSSDSensorHeight,
                                                0.5*fgkSSDSensorInsensitiveWidth);
  const char* ssdsensorinsensitivename[2] = {"SSDSensorInsensitive1",
                                             "SSDSensorInsensitive2"};
  TGeoVolume* ssdsensorinsensitive[2];
  for(Int_t i=0; i<2; i++){ ssdsensorinsensitive[i] = 
      new TGeoVolume(ssdsensorinsensitivename[i],ssdsensorinsensitiveshape[i],
                     fSSDSensorMedium);
      ssdsensorinsensitive[i]->SetLineColor(fColorCarbonFiber);
  }
  /////////////////////////////////////////////////////////////
  // Virtual Volume containing SSD Sensor  
  /////////////////////////////////////////////////////////////
  TGeoBBox* virtualSSDSensorShape = new TGeoBBox("SSDSensorShape",
											     0.5*fgkSSDSensorWidth,
											     0.5*fgkSSDSensorHeight,
											     0.5*fgkSSDSensorLength);
  fSSDSensor5 = new TGeoVolume("ITSsddSensor5",virtualSSDSensorShape,
										 fSSDAir);	
  fSSDSensor6 = new TGeoVolume("ITSsddSensor6",virtualSSDSensorShape,
										 fSSDAir);	
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i<4; i++){ 
            fSSDSensor5->AddNode(i%2==0 ? ssdsensorinsensitive[0]:
            ssdsensorinsensitive[1],i<2?1:2,
			new TGeoTranslation(
			 0.5*(1.+TMath::Power(-1.,i))*(i==0?-1.: 1.)
      *   (ssdsensorsensitiveshape->GetDX()+ssdsensorinsensitiveshape[0]->GetDX()),0.,			
			0.5*(1.-TMath::Power(-1.,i))*(i==1? 1.:-1.)
      *   (ssdsensorsensitiveshape->GetDZ()+ssdsensorinsensitiveshape[1]->GetDZ())));    
            fSSDSensor6->AddNode(i%2==0 ? ssdsensorinsensitive[0]:
            ssdsensorinsensitive[1],i<2?1:2,
			new TGeoTranslation(
			 0.5*(1.+TMath::Power(-1.,i))*(i==0?-1.: 1.)
      *   (ssdsensorsensitiveshape->GetDX()+ssdsensorinsensitiveshape[0]->GetDX()),0.,			
			0.5*(1.-TMath::Power(-1.,i))*(i==1? 1.:-1.)
      *   (ssdsensorsensitiveshape->GetDZ()+ssdsensorinsensitiveshape[1]->GetDZ())));    
  }
    fSSDSensor5->AddNode(ssdsensorsensitiveLay5,1);
    fSSDSensor6->AddNode(ssdsensorsensitiveLay6,1);
}
///////////////////////////////////////////////////////////////////////////////
TList* AliITSv11GeometrySSD::GetCarbonFiberSupportList(){
  /////////////////////////////////////////////////////////////  
  // Method generating the Carbon Fiber Support   
  /////////////////////////////////////////////////////////////  
  const Int_t kvertexnumber = 4;
  const Int_t kshapesnumber = 2;
  TVector3** vertexposition[kshapesnumber];
  for(Int_t i=0; i<kshapesnumber; i++) vertexposition[i] = new TVector3*[kvertexnumber];
  Double_t carbonfibersupportxaxisEdgeproj = 
		fgkCarbonFiberSupportEdgeLength*TMath::Cos(fgkCarbonFiberJunctionAngle[0]
	*	TMath::DegToRad());
  Double_t theta = TMath::ATan(fgkCarbonFiberSupportYAxisLength
				 /			   fgkCarbonFiberSupportXAxisLength);
  /////////////////////
  //Vertex Positioning
  ////////////////////
  vertexposition[0][0] = new TVector3();
  vertexposition[0][1] = new TVector3(fgkCarbonFiberSupportXAxisLength,
									  fgkCarbonFiberSupportYAxisLength);
  vertexposition[0][2] = new TVector3(carbonfibersupportxaxisEdgeproj,
									  carbonfibersupportxaxisEdgeproj
					   *			  TMath::Tan(theta));
  vertexposition[0][3] = new TVector3(fgkCarbonFiberSupportXAxisLength
					   -			  carbonfibersupportxaxisEdgeproj,
									  fgkCarbonFiberSupportYAxisLength
					   -			  vertexposition[0][2]->Y());
  ////////////////////////////////////////////////////
  //Setting the parameters for Isometry Transformation
  ////////////////////////////////////////////////////
  Double_t symmetryplaneposition = (fgkCarbonFiberSupportYAxisLength
								 +	fgkCarbonFiberSupportTopEdgeDist[0]
								 +	fgkCarbonFiberSupportWidth);
  Double_t* param = new Double_t[4]; 
  param[0] = 0., param[1] = 1., param[2] = 0., param[3] = -symmetryplaneposition;
  for(Int_t j=0; j<kvertexnumber; j++) vertexposition[1][j] = 
				  new TVector3((GetReflection(vertexposition[0][j],param))->X(),
							  (GetReflection(vertexposition[0][j],param))->Y());
  char* carbonfibersupportshapename[kshapesnumber] = 
						{"CarbonFiberSupportShape1","CarbonFiberSupportShape2"};
  char* carbonfibersupportname[kshapesnumber] = 
						{"CarbonFiberSupport1","CarbonFiberSupport2"};
  TGeoArb8* carbonfibersupportshape[kshapesnumber]; 
  TGeoVolume* carbonfibersupport[kshapesnumber];
  TList* carbonfibersupportlist = new TList();
  Double_t width[2] = {fgkCarbonFiberSupportWidth,fgkCarbonFiberSupportWidth};
  Double_t carbonfibersupportheight = 
	  carbonfibersupportxaxisEdgeproj*TMath::Tan(fgkCarbonFiberJunctionAngle[0]
	  *TMath::DegToRad());
  for(Int_t i = 0; i< kshapesnumber; i++){
   carbonfibersupportshape[i] = 
					GetArbShape(vertexposition[i],width,carbonfibersupportheight,
								carbonfibersupportshapename[i],i==0 ? 1: -1);
   carbonfibersupport[i] = new TGeoVolume(carbonfibersupportname[i],
						   carbonfibersupportshape[i],fSSDCarbonFiberMedium);
   carbonfibersupport[i]->SetLineColor(fColorCarbonFiber);
   carbonfibersupportlist->Add(carbonfibersupport[i]);	
   }
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i< kshapesnumber; i++){
     for(Int_t j=0; j< kvertexnumber; j++)
	   delete vertexposition[i][j];
       delete [] vertexposition[i];
  }
  delete [] param;
  /////////////////////////////////////////////////////////////
   return carbonfibersupportlist;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCarbonFiberJunction(Double_t width){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Carbon Fiber Junction
  /////////////////////////////////////////////////////////////
  const Int_t kvertexnumber = 6;
  TGeoXtru* carbonfiberjunctionshape = new TGeoXtru(2);
  Double_t reflectionparam[3] = {TMath::Tan(fgkCarbonFiberJunctionAngle[0]
							  *  TMath::DegToRad()),-1.,0.};
  TVector3* vertex[kvertexnumber];
  vertex[0] = new TVector3();
  vertex[3] = new TVector3(fgkCarbonFiberJunctionEdge[0]
			*			  TMath::Cos(fgkCarbonFiberJunctionAngle[0]
			*			  TMath::DegToRad()),
						  fgkCarbonFiberJunctionEdge[0]
			*			  TMath::Sin(fgkCarbonFiberJunctionAngle[0]
			*			  TMath::DegToRad()));
  vertex[4] = new TVector3(fgkCarbonFiberJunctionLength,
						   fgkCarbonFiberJunctionEdge[1]);
  vertex[5] = new TVector3(fgkCarbonFiberJunctionLength); 
  vertex[1] = GetReflection(vertex[5],reflectionparam);	
  vertex[2] = GetReflection(vertex[4],reflectionparam);	
  Double_t xvertexpoints[6], yvertexpoints[6];
  for(Int_t i=0; i<kvertexnumber; i++) 
	  xvertexpoints[i] = vertex[i]->X(), yvertexpoints[i] = vertex[i]->Y();
  carbonfiberjunctionshape->DefinePolygon(kvertexnumber,xvertexpoints,yvertexpoints);
  carbonfiberjunctionshape->DefineSection(0,-0.5*width);
  carbonfiberjunctionshape->DefineSection(1,0.5*width);
  TGeoVolume* carbonfiberjunction = new TGeoVolume("CarbonFiberJunction",
								carbonfiberjunctionshape,fSSDCarbonFiberMedium);
  carbonfiberjunction->SetLineColor(fColorCarbonFiber);
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  for (Int_t i=0; i<kvertexnumber; i++) delete vertex[i];
  ///////////////////////////////////////////////////////////// 
  return carbonfiberjunction;
}
////////////////////////////////////////////////////////////////////////////////
TList* AliITSv11GeometrySSD::GetCarbonFiberLowerSupportList(){
  /////////////////////////////////////////////////////////////
  // Method generating the Carbon Fiber Lower Support   
  /////////////////////////////////////////////////////////////  
  const Int_t kvertexnumber = 4;
  const Int_t kshapesnumber = 2;
  Double_t width[2] = {fgkCarbonFiberLowerSupportWidth,
								fgkCarbonFiberLowerSupportWidth};
  TVector3** vertexposition[kshapesnumber];
  for(Int_t i = 0; i<kshapesnumber; i++) vertexposition[i] = 
						 new TVector3*[kvertexnumber];
  //First Shape Vertex Positioning
  vertexposition[0][0] = new TVector3(fgkCarbonFiberLowerSupportLowerLenght);
  vertexposition[0][1] = new TVector3(fgkCarbonFiberTriangleLength
					   -		fgkCarbonFiberLowerSupportLowerLenght);
  vertexposition[0][2] = new TVector3();
  vertexposition[0][3] = new TVector3(fgkCarbonFiberTriangleLength);
  //Second Shape Vertex Positioning
  Double_t theta = TMath::ATan((fgkCarbonFiberLowerSupportVolumePosition[1]
				 -				fgkCarbonFiberLowerSupportVolumePosition[0])
				 /				fgkCarbonFiberTriangleLength);
  vertexposition[1][0] = new TVector3(vertexposition[0][0]->X(),
								vertexposition[0][0]->X()*TMath::Tan(theta)
				 +				fgkCarbonFiberLowerSupportVolumePosition[0]);
  vertexposition[1][1] = new TVector3(vertexposition[0][1]->X(),
								vertexposition[0][1]->X()*TMath::Tan(theta)
				 +				fgkCarbonFiberLowerSupportVolumePosition[0]);
  vertexposition[1][2] = new TVector3(0.,fgkCarbonFiberLowerSupportVolumePosition[0]);
  vertexposition[1][3] = new TVector3(fgkCarbonFiberTriangleLength,
								fgkCarbonFiberLowerSupportVolumePosition[1]);
  char* carbonfiberlowersupportshapename[kshapesnumber] = 
			  {"CarbonFiberLowerSupportShape1","CarbonFiberLowerSupportShape2"};
  char* carbonfiberlowersupportname[kshapesnumber] = 
			  {"CarbonFiberLowerSupport1","CarbonFiberLowerSupport2"};
  TGeoArb8* carbonfiberlowersupportshape[kshapesnumber];
  TGeoVolume* carbonfiberlowersupport[kshapesnumber];
  TList* carbonfiberlowersupportlist = new TList();
  for(Int_t i = 0; i< kshapesnumber; i++){ 
	carbonfiberlowersupportshape[i] = 
								GetArbShape(vertexposition[i],width,
											fgkCarbonFiberLowerSupportHeight,
											carbonfiberlowersupportshapename[i]);
    carbonfiberlowersupport[i] = new TGeoVolume(carbonfiberlowersupportname[i],
						carbonfiberlowersupportshape[i],fSSDCarbonFiberMedium);
	carbonfiberlowersupport[i]->SetLineColor(fColorCarbonFiber);
    carbonfiberlowersupportlist->Add(carbonfiberlowersupport[i]);
  }
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i< kshapesnumber; i++){
     for(Int_t j=0; j< kvertexnumber; j++)
	   delete vertexposition[i][j];
       delete [] vertexposition[i];
  }
  /////////////////////////////////////////////////////////////
  return carbonfiberlowersupportlist;
}
///////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDSensorSupport(Double_t length, Double_t height, 
								 Double_t width, Double_t* thickness)const{
  /////////////////////////////////////////////////////////////
  // Method generating the Sensor Support   
  /////////////////////////////////////////////////////////////  
	const Int_t kvertexnumber = 6;
	TGeoXtru* ssdsensorsupportshape = new TGeoXtru(2);	
    TVector3* vertexposition[kvertexnumber];
	vertexposition[0] = new TVector3();	
	vertexposition[1] = new TVector3(0.0,length);	
	vertexposition[2] = new TVector3(thickness[1],vertexposition[1]->Y());	
	vertexposition[3] = new TVector3(vertexposition[2]->X(),thickness[0]);	
	vertexposition[4] = new TVector3(height,vertexposition[3]->Y());	
	vertexposition[5] = new TVector3(vertexposition[4]->X());	
	Double_t xvertexpoints[6], yvertexpoints[6];
	for(Int_t i=0; i<kvertexnumber; i++) 
		xvertexpoints[i] = vertexposition[i]->X(), 
		yvertexpoints[i] = vertexposition[i]->Y();
    ssdsensorsupportshape->DefinePolygon(kvertexnumber,xvertexpoints,yvertexpoints);
    ssdsensorsupportshape->DefineSection(0,-0.5*width);
    ssdsensorsupportshape->DefineSection(1,0.5*width);
    TGeoVolume* ssdsensorsupport = new TGeoVolume("SSDSensorSupport",
								 ssdsensorsupportshape,fSSDSensorSupportMedium);
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
	for (Int_t i=0; i<kvertexnumber; i++)
		delete vertexposition[i];
  /////////////////////////////////////////////////////////////
    return ssdsensorsupport;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCoolingTubeSupport(Int_t nedges){
  /////////////////////////////////////////////////////////////
  // Method generating the Cooling Tube Support
  /////////////////////////////////////////////////////////////
  if(nedges%2!=0) nedges--;	
  const Int_t kvertexnumber = nedges+5;
  Double_t phi = TMath::ASin(0.5*fgkCoolingTubeSupportHeight
			   /			 fgkCoolingTubeSupportRmax)*TMath::RadToDeg();
  Double_t angle = 90.+phi;
  Double_t psi = 90.-phi;
  ///////////////////////////////////////
  // Vertex Positioning for TGeoXTru
  ///////////////////////////////////////
  TVector3** vertexposition = new TVector3*[kvertexnumber];
  vertexposition[0] = new TVector3(fgkCoolingTubeSupportRmin*CosD(angle),
								   fgkCoolingTubeSupportRmin*SinD(angle));
  vertexposition[1] = new TVector3(fgkCoolingTubeSupportRmax*CosD(angle),
								   fgkCoolingTubeSupportRmax*SinD(angle));
  vertexposition[2] = new TVector3(vertexposition[1]->X(),
								   fgkCoolingTubeSupportRmax);
  vertexposition[3] = new TVector3(-vertexposition[1]->X(),
								   fgkCoolingTubeSupportRmax);
  vertexposition[4] = new TVector3(-vertexposition[1]->X(),
								    vertexposition[1]->Y());
  for(Int_t i=0; i<nedges; i++)
	vertexposition[i+5] = 
		new TVector3(fgkCoolingTubeSupportRmin*CosD(psi+i*(2.*phi/nedges)),
					 fgkCoolingTubeSupportRmin*SinD(psi+i*(2.*phi/nedges)));
  ///////////////////////////////////////////////////////////////////////
  // TGeoXTru Volume definition for Cooling Tube Support Arc Part
  ///////////////////////////////////////////////////////////////////////
  TGeoXtru* coolingtubesupportarcshape = new TGeoXtru(2);	
  Double_t* xvertexpoints = new Double_t[kvertexnumber];	
  Double_t* yvertexpoints = new Double_t[kvertexnumber];	
  for(Int_t i=0; i<kvertexnumber; i++){
	xvertexpoints[i] = vertexposition[i]->X();
	yvertexpoints[i] = vertexposition[i]->Y();
  } 
  coolingtubesupportarcshape->DefinePolygon(kvertexnumber,xvertexpoints,
											yvertexpoints);
  coolingtubesupportarcshape->DefineSection(0,-0.5*fgkCoolingTubeSupportWidth);
  coolingtubesupportarcshape->DefineSection(1,0.5*fgkCoolingTubeSupportWidth);
  TGeoVolume* coolingtubesupportarc = new TGeoVolume("CoolingTubeSupportArc",
								          coolingtubesupportarcshape,
										  fSSDTubeHolderMedium);
  coolingtubesupportarc->SetLineColor(fColorG10);
  //////////////////////////////////////////////////////////////////////////
  // TGeoTubeSeg Volume definition for Cooling Tube Support TGeoTubeSeg Part
  //////////////////////////////////////////////////////////////////////////
  TGeoTubeSeg* coolingtubesupportsegshape = 
							new TGeoTubeSeg(fgkCoolingTubeSupportRmin,
											fgkCoolingTubeSupportRmax,
											0.5*fgkCoolingTubeSupportWidth,
											phi,360-phi);
  TGeoVolume* coolingtubesupportseg = new TGeoVolume("CoolingTubeSupportSeg",
											coolingtubesupportsegshape,
											fSSDTubeHolderMedium);
  coolingtubesupportseg->SetLineColor(fColorG10);
  //////////////////////////////////////////////////////////////////////////
  // TGeoBBox Volume definition for Cooling Tube Support Box Part
  //////////////////////////////////////////////////////////////////////////
  Double_t* boxorigin = new Double_t[3];
  Double_t boxlength = fgkCoolingTubeSupportLength-2.*fgkCoolingTubeSupportRmax;
  boxorigin[0] = fgkCoolingTubeSupportRmax+0.5*boxlength, boxorigin[1] = boxorigin[2] = 0.;
  TGeoBBox* coolingtubesupportboxshape = new TGeoBBox(0.5*boxlength,
										 0.5*fgkCoolingTubeSupportHeight,
										 0.5*fgkCoolingTubeSupportWidth,boxorigin);
  TGeoVolume* coolingtubesupportbox = new TGeoVolume("CoolingTubeSupportBox",
                               coolingtubesupportboxshape,fSSDTubeHolderMedium);
  coolingtubesupportbox->SetLineColor(fColorG10);
  //////////////////////////////////////////////////////////////////////////
  // Cooling Tube for Cooling Tube Support 
  //////////////////////////////////////////////////////////////////////////
  TGeoXtru* coolingtubearcshape[2];
  coolingtubearcshape[0] = new TGeoXtru(2);	
  Double_t* xvert = new Double_t[nedges+2];
  Double_t* yvert = new Double_t[nedges+2];
  Double_t ratio = fgkCoolingTubeRmin/fgkCoolingTubeSupportRmin;
  ////////////////////////////////////////
  // Positioning the vertices for TGeoXTru
  ////////////////////////////////////////
  xvert[0] = 0., yvert[0] = 0.;
  xvert[1] = vertexposition[0]->X()*ratio,  yvert[1] = vertexposition[0]->Y()*ratio;
  for(Int_t i=0; i< nedges; i++)
		xvert[i+2] = vertexposition[kvertexnumber-i-1]->X()*ratio,
		yvert[i+2] = vertexposition[kvertexnumber-i-1]->Y()*ratio;
  ////////////////////////////////////////
  // Defining TGeoXTru PolyGone
  ////////////////////////////////////////
  coolingtubearcshape[0]->DefinePolygon(nedges+2,xvert,yvert);
  coolingtubearcshape[0]->DefineSection(0,-0.5*fgkCoolingTubeSupportWidth);
  coolingtubearcshape[0]->DefineSection(1,0.5*fgkCoolingTubeSupportWidth);
  coolingtubearcshape[1] = GetArcShape(2.*phi,fgkCoolingTubeRmin,
		fgkCoolingTubeRmax,nedges,fgkCoolingTubeSupportWidth);
  TGeoVolume* coolingtubearc[2];
  coolingtubearc[0] = new TGeoVolume("CoolingTubeWaterArcPart",
								  coolingtubearcshape[0],fSSDCoolingTubeWater);
  coolingtubearc[1] = new TGeoVolume("CoolingTubePhynoxArcPart",
								  coolingtubearcshape[1],fSSDCoolingTubePhynox);
  coolingtubearc[0]->SetLineColor(fColorWater);
  coolingtubearc[1]->SetLineColor(fColorPhynox);
  ////////////////////////////////////////////
  // Defining TGeoTubeSeg Part of Cooling Tube
  ////////////////////////////////////////////
  TGeoTubeSeg* coolingtubesegshape[2];
  coolingtubesegshape[0] = new TGeoTubeSeg(fgkCoolingTubeRmin,fgkCoolingTubeRmax,
							0.5*fgkCoolingTubeSupportWidth,phi,360-phi);
  coolingtubesegshape[1] = new TGeoTubeSeg(0.,fgkCoolingTubeRmin,
							0.5*fgkCoolingTubeSupportWidth,phi,360-phi);
  TGeoVolume* coolingtubeseg[2];
  coolingtubeseg[0] = new TGeoVolume("CoolingTubePhynoxPart",
								 coolingtubesegshape[0],fSSDCoolingTubePhynox);
  coolingtubeseg[1] = new TGeoVolume("CoolingTubeWaterPart",
								 coolingtubesegshape[1],fSSDCoolingTubeWater);
  coolingtubeseg[0]->SetLineColor(fColorPhynox);
  coolingtubeseg[1]->SetLineColor(fColorWater);
  /////////////////////////////////////////////////////////////
  // Virtual Volume containing Cooling Tube Support  
  /////////////////////////////////////////////////////////////
  TGeoXtru* virtualCoolingTubeSupportShape = new TGeoXtru(2);
  const Int_t kvirtualvertexnumber = 8;
  TVector3* virtualvertex[kvirtualvertexnumber];
   ////////////////////////////////////////
  // Positioning the vertices for TGeoXTru
  ////////////////////////////////////////
  virtualvertex[0] = new TVector3(-fgkCoolingTubeSupportRmax,-fgkCoolingTubeSupportRmax); 
  virtualvertex[1] = new TVector3(virtualvertex[0]->X(),-virtualvertex[0]->Y());
  virtualvertex[2] = new TVector3(-virtualvertex[0]->X(),virtualvertex[1]->Y());
  virtualvertex[3] = new TVector3(virtualvertex[2]->X(),0.5*fgkCoolingTubeSupportHeight);
  virtualvertex[4] = new TVector3(virtualvertex[3]->X()+boxlength,virtualvertex[3]->Y());
  virtualvertex[5] = new TVector3(virtualvertex[4]->X(),-virtualvertex[4]->Y());
  virtualvertex[6] = new TVector3(virtualvertex[3]->X(),-virtualvertex[3]->Y());
  virtualvertex[7] = new TVector3(virtualvertex[2]->X(),-virtualvertex[2]->Y());
  Double_t xmothervertex[kvirtualvertexnumber], ymothervertex[kvirtualvertexnumber];
  for(Int_t i=0; i< kvirtualvertexnumber; i++)
	xmothervertex[i] = virtualvertex[i]->X(),
	ymothervertex[i] = virtualvertex[i]->Y();
  ////////////////////////////////////////
  // Defining TGeoXTru PolyGone
  ////////////////////////////////////////
  virtualCoolingTubeSupportShape->DefinePolygon(kvirtualvertexnumber,xmothervertex,
																	 ymothervertex);
  virtualCoolingTubeSupportShape->DefineSection(0,-0.5*fgkCoolingTubeSupportWidth);
  virtualCoolingTubeSupportShape->DefineSection(1,0.5*fgkCoolingTubeSupportWidth);
  TGeoVolume* virtualcoolingtubesupport = new TGeoVolume("CoolingTubeSupport",
								 virtualCoolingTubeSupportShape,fSSDAir);
  ////////////////////////////////////////
  // Positioning Volumes in Virtual Volume
  ////////////////////////////////////////
  TGeoRotation* coolingtubesupportrot = new TGeoRotation(); 
  coolingtubesupportrot->SetAngles(-90.0,0.0,0.0);
  virtualcoolingtubesupport->AddNode(coolingtubesupportarc,1,coolingtubesupportrot);
  virtualcoolingtubesupport->AddNode(coolingtubesupportbox,1);
  virtualcoolingtubesupport->AddNode(coolingtubesupportseg,1);
  virtualcoolingtubesupport->AddNode(coolingtubearc[0],1,coolingtubesupportrot);
  virtualcoolingtubesupport->AddNode(coolingtubearc[1],1,coolingtubesupportrot);
  virtualcoolingtubesupport->AddNode(coolingtubeseg[0],1);
  virtualcoolingtubesupport->AddNode(coolingtubeseg[1],1);
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  delete [] vertexposition;
  delete xvertexpoints;
  delete yvertexpoints;
  delete xvert;
  delete yvert;
  for(Int_t i=0; i< kvirtualvertexnumber; i++)
	delete virtualvertex[i];
  /////////////////////////////////////////////////////////////
	return virtualcoolingtubesupport;
}
/////////////////////////////////////////////////////////////////////////////////
TList* AliITSv11GeometrySSD::GetSSDHybridParts(){
  /////////////////////////////////////////////////////////////
  // Method generating List containing SSD Hybrid Components   
  /////////////////////////////////////////////////////////////
  TList* ssdhybridlist = new TList();
  const Int_t kssdstiffenernumber = 2;
  Double_t ssdstiffenerseparation = fgkSSDSensorLength
								  - 2.*fgkSSDModuleStiffenerPosition[1]
								  -    fgkSSDStiffenerWidth;
  Double_t ssdchipcablesradius[kssdstiffenernumber];
  for(Int_t i=0; i<kssdstiffenernumber; i++)
	  ssdchipcablesradius[i] = 0.25*(fgkSSDChipCablesHeight[i+2]
			       -  fgkSSDChipCablesHeight[0]
			       -  fgkSSDChipCablesHeight[1]);
  /////////////////////////////////////////////////////////////
  // Mother Volumes Containers 
  /////////////////////////////////////////////////////////////
  const Int_t kmothernumber = 2;
  const Int_t kmothervertexnumber = 12;
  Double_t xmothervertex[kmothernumber][kmothervertexnumber]; 
  Double_t ymothervertex[kmothernumber][kmothervertexnumber]; 
  ///////////////////////
  // Setting the vertices 
  ///////////////////////
  xmothervertex[0][0]  = -0.5*fgkSSDStiffenerLength;
  xmothervertex[0][1]  = xmothervertex[0][0]; 
  xmothervertex[0][2]  = fgkSSDFlexLength[0]-0.5*fgkSSDStiffenerLength;
  xmothervertex[0][3]  = xmothervertex[0][2];
  xmothervertex[0][4]  = xmothervertex[0][0];
  xmothervertex[0][5]  = xmothervertex[0][4];
  xmothervertex[0][6]  = -xmothervertex[0][0];
  xmothervertex[0][7]  = xmothervertex[0][6];
  xmothervertex[0][8]  = -xmothervertex[0][2];
  xmothervertex[0][9]  = xmothervertex[0][8];
  xmothervertex[0][10] = xmothervertex[0][7];
  xmothervertex[0][11] = xmothervertex[0][10];
  for(Int_t i=0; i<kmothervertexnumber; i++) xmothervertex[1][i] = xmothervertex[0][i];
  for(Int_t i = 0; i<kmothernumber; i++){
      ymothervertex[i][0]  = -(fgkSSDChipWidth-(0.5*fgkSSDStiffenerWidth-fgkSSDStiffenerToChipDist)
                           + ssdchipcablesradius[i]+fgkSSDChipCablesWidth[1]+fgkSSDChipCablesWidth[2]);
      ymothervertex[i][1]  = ssdstiffenerseparation-0.5*fgkSSDStiffenerWidth-fgkSSDFlexWidth[0];
      ymothervertex[i][2]  = ymothervertex[i][1];
      ymothervertex[i][3]  = ssdstiffenerseparation-0.5*fgkSSDStiffenerWidth;
      ymothervertex[i][4]  = ymothervertex[i][3];
      ymothervertex[i][5]  = ymothervertex[i][4]+0.5*fgkSSDStiffenerWidth-ymothervertex[i][0];
      ymothervertex[i][6]  = ymothervertex[i][5];
      ymothervertex[i][7]  = 0.5*fgkSSDStiffenerWidth+fgkSSDFlexWidth[0];
      ymothervertex[i][8]  = ymothervertex[i][7];
      ymothervertex[i][9]  = 0.5*fgkSSDStiffenerWidth;
      ymothervertex[i][10] = ymothervertex[i][9];
      ymothervertex[i][11] = ymothervertex[i][0];
  }
  TGeoXtru* ssdhybridmothershape[kmothernumber];
  TGeoVolume* ssdhybridmother[kmothernumber];
  const char* ssdhybridmothername[kmothernumber] = {"SSDHybridMother1","SSDHybridMother2"};
  for(Int_t i=0; i<kmothernumber; i++){
      ssdhybridmothershape[i] = new TGeoXtru(2);
      ssdhybridmothershape[i]->DefinePolygon(kmothervertexnumber,xmothervertex[i],
                                          ymothervertex[i]);
      ssdhybridmothershape[i]->DefineSection(0,-0.5*fgkSSDStiffenerHeight-fgkSSDChipHeight      
                                               -fgkSSDChipCablesHeight[i+2]);
      ssdhybridmothershape[i]->DefineSection(1, 0.5*fgkSSDStiffenerHeight);
      ssdhybridmother[i] = new TGeoVolume(ssdhybridmothername[i],ssdhybridmothershape[i],
                                          fSSDAir);
   }   
  /////////////////////////////////////////////////////////////
  // SSD Stiffener   
  /////////////////////////////////////////////////////////////
  TGeoBBox* ssdstiffenershape = new TGeoBBox("SSDStiffenerShape",
                                             0.5*fgkSSDStiffenerLength,
                                             0.5*fgkSSDStiffenerWidth,
                                             0.5*fgkSSDStiffenerHeight);
  TGeoVolume* ssdstiffener = new TGeoVolume("SSDStiffener",ssdstiffenershape,
                                            fSSDStiffenerMedium);  
  ssdstiffener->SetLineColor(fColorStiffener); 
  TGeoTranslation* ssdstiffenertrans[kssdstiffenernumber];
  for(Int_t i=0; i<kssdstiffenernumber; i++) 
      ssdstiffenertrans[i] = new TGeoTranslation(0.,i*ssdstiffenerseparation,0.);
  /////////////////////////////////////////////////////////////
  // SSD Chip System	
  /////////////////////////////////////////////////////////////
  TList* ssdchipsystemlist = GetSSDChipSystem(); 
  Double_t ssdchipseparation = fgkSSDSensorLength
                             - 2.*fgkSSDModuleStiffenerPosition[1]
                             - 2.*(fgkSSDStiffenerWidth-fgkSSDStiffenerToChipDist
                             - 0.5*fgkSSDChipWidth);
  Double_t ssdchipsystemlength = (fgkSSDChipNumber-1)*(fgkSSDChipLength 
			       +  fgkSSDChipSeparationLength)+fgkSSDChipCablesLength[1];
  TGeoTranslation* ssdchipsystemtrans = new TGeoTranslation(0.5*fgkSSDChipCablesLength[1]
                                      - 0.5*ssdchipsystemlength,
                                        0.5*(ssdstiffenerseparation-ssdchipseparation),
                                      - 0.5*(fgkSSDChipHeight+fgkSSDStiffenerHeight)); 	
////////////////////////////
// Capacitor 0603-2200 nF
///////////////////////////
  const Int_t knapacitor0603number = 5;
  TGeoBBox* capacitor0603shape =  new TGeoBBox("Capacitor0603Shape",
											 0.5*fgkSSDCapacitor0603Length,
											 0.5*fgkSSDCapacitor0603Width,
											 0.5*fgkSSDCapacitor0603Height);
  TGeoVolume* capacitor0603 = new TGeoVolume("Capacitor0603",capacitor0603shape,
                                             fSSDStiffener0603CapacitorMedium); 
  capacitor0603->SetLineColor(fColorAl);
  for(Int_t i=0; i<kmothernumber; i++){
      for(Int_t j=0; j<kssdstiffenernumber; j++){
            ssdhybridmother[i]->AddNode(ssdstiffener,j+1,ssdstiffenertrans[j]);
            for(Int_t k=1; k<knapacitor0603number+1; k++){
                  ssdhybridmother[i]->AddNode(capacitor0603,knapacitor0603number*j+k,
                        new TGeoTranslation((k-3.)/6*fgkSSDStiffenerLength,
                                             j*ssdstiffenerseparation
                        +                    0.5*((j==0? 1:-1)*fgkSSDStiffenerWidth
                        +                    (j==0? -1:+1)*fgkSSDCapacitor0603Width),
                        -                    0.5*(fgkSSDStiffenerHeight+fgkSSDCapacitor0603Height)));
            }
      } 
      ssdhybridmother[i]->AddNode((TGeoVolume*)ssdchipsystemlist->At(i),i+1,ssdchipsystemtrans);
      ssdhybridlist->Add(ssdhybridmother[i]);
  }    
/////////////////////////////////////////////////////////////
// Mother Volume Containing Capacitor Part 
/////////////////////////////////////////////////////////////
  const Int_t kcapacitormothernumber = 8;
  Double_t xcapacitorvertex[kcapacitormothernumber];
  Double_t ycapacitorvertex[kcapacitormothernumber];  
  ///////////////////////
  // Setting the vertices 
  ///////////////////////
  xcapacitorvertex[0] = -fgkSSDConnectorPosition[0]+ssdstiffenershape->GetDX();    
  xcapacitorvertex[1] = xcapacitorvertex[0];   
  xcapacitorvertex[2] = 0.5*fgkSSDFlexHoleWidth;   
  xcapacitorvertex[3] = xcapacitorvertex[2];   
  xcapacitorvertex[4] = xcapacitorvertex[0];   
  xcapacitorvertex[5] = xcapacitorvertex[0];   
  xcapacitorvertex[6] = -xcapacitorvertex[0];   
  xcapacitorvertex[7] = xcapacitorvertex[6];   
  ycapacitorvertex[0] = -0.5*fgkSSDStiffenerWidth;    
  ycapacitorvertex[1] = ssdstiffenerseparation-0.5*fgkSSDStiffenerWidth-fgkSSDFlexHoleLength;   
  ycapacitorvertex[2] = ycapacitorvertex[1];   
  ycapacitorvertex[3] = ycapacitorvertex[2]+fgkSSDFlexHoleLength;   
  ycapacitorvertex[4] = ycapacitorvertex[3];   
  ycapacitorvertex[5] = ycapacitorvertex[4]+fgkSSDStiffenerWidth;   
  ycapacitorvertex[6] = ycapacitorvertex[5];   
  ycapacitorvertex[7] = ycapacitorvertex[0];   
  TGeoXtru* ssdhybridcapacitormothershape = new TGeoXtru(2);
  ssdhybridcapacitormothershape->DefinePolygon(kcapacitormothernumber,xcapacitorvertex,
                                              ycapacitorvertex);
  ssdhybridcapacitormothershape->DefineSection(0,0.5*fgkSSDStiffenerHeight);
  ssdhybridcapacitormothershape->DefineSection(1, 0.5*fgkSSDStiffenerHeight+fgkSSDCapacitor1812Height);
  TGeoVolume* ssdhybridcapacitormother = new TGeoVolume("SSDHybridCapacitorMother",ssdhybridcapacitormothershape,
                                          fSSDAir);
////////////////////////////
// Connector 
///////////////////////////
  const Int_t kssdconnectornumber = 2;
  TGeoBBox* ssdconnectorshape[kssdconnectornumber];
  Double_t ssdAlconnectororigin[3] = {0.0,0.0,0.5*(fgkSSDStiffenerHeight+fgkSSDConnectorAlHeight)};    
  Double_t ssdNiconnectororigin[3] = {0.0,0.0,0.5*(fgkSSDStiffenerHeight+fgkSSDConnectorNiHeight)
                                   +  fgkSSDConnectorAlHeight};  
  const char* ssdconnectorname[kssdconnectornumber] = {"SSDConnectorAl","SSDConnectorNi"};
  TGeoVolume* ssdconnector[kssdconnectornumber];
  for(Int_t i=0; i<kssdconnectornumber; i++){
      ssdconnectorshape[i] = new TGeoBBox(0.5*fgkSSDConnectorLength,
                                          0.5*fgkSSDConnectorWidth,
                                          0.5*((1-i)*fgkSSDConnectorAlHeight
                           +              i*fgkSSDConnectorNiHeight),
                             i==0? ssdAlconnectororigin : ssdNiconnectororigin);
      ssdconnector[i] = new TGeoVolume(ssdconnectorname[i],ssdconnectorshape[i],
                                       i==0 ? fSSDAlTraceFlexMedium 
                                            : fSSDStiffenerConnectorMedium);      
      ssdconnector[i]->SetLineColor(i==0 ? fColorAl : fColorPhynox);
  }
  TGeoTranslation* ssdconnectortrans[2*kssdconnectornumber];
  ssdconnectortrans[0] = new TGeoTranslation(-ssdstiffenershape->GetDX()
                       +  fgkSSDConnectorPosition[0]
                       -  fgkSSDConnectorSeparation
                       -  1.5*fgkSSDConnectorLength,
                          ssdstiffenerseparation+ssdstiffenershape->GetDY()
                       -  fgkSSDConnectorPosition[1]
                       -  ssdconnectorshape[0]->GetDY(),0.0);	
  ssdconnectortrans[1] = new TGeoTranslation(
                       -  ssdstiffenershape->GetDX()
                       +  fgkSSDConnectorPosition[0]
                       -  0.5*fgkSSDConnectorLength,
                          ssdstiffenerseparation+ssdstiffenershape->GetDY()
                       -  fgkSSDConnectorPosition[1]
                       -  ssdconnectorshape[0]->GetDY(),0.0);
  ssdconnectortrans[2] = new TGeoTranslation(+ssdstiffenershape->GetDX()
                       -  fgkSSDConnectorPosition[0]
                       +  fgkSSDConnectorSeparation
                       +  1.5*fgkSSDConnectorLength,
                          -(ssdstiffenershape->GetDY()
                       -  fgkSSDConnectorPosition[1]
                       -  ssdconnectorshape[0]->GetDY()),0.0);	
  ssdconnectortrans[3] = new TGeoTranslation(+ssdstiffenershape->GetDX()
                       -  fgkSSDConnectorPosition[0]
                       +  0.5*fgkSSDConnectorLength,
                          -(ssdstiffenershape->GetDY()
                       -  fgkSSDConnectorPosition[1]
                       -  ssdconnectorshape[0]->GetDY()),0.0);
  for(Int_t i=0; i<2*kssdconnectornumber; i++)
      for(Int_t j=0; j<kssdconnectornumber; j++)
            ssdhybridcapacitormother->AddNode(ssdconnector[j],i+1,ssdconnectortrans[i]);      
////////////////////////////
// Capacitor 1812-330 nF
/////////////////////////// 
  Double_t ssdcapacitor1812origin[3] = {0.0,0.0,0.5*(fgkSSDStiffenerHeight+fgkSSDCapacitor1812Height)};    
  TGeoBBox* capacitor1812shape =  new TGeoBBox("Capacitor1812Shape",
											 0.5*fgkSSDCapacitor1812Length,
											 0.5*fgkSSDCapacitor1812Width,
											 0.5*fgkSSDCapacitor1812Height,
            ssdcapacitor1812origin);
  TGeoVolume* capacitor1812 = new TGeoVolume("Capacitor1812",capacitor1812shape,
                                             fSSDStiffener1812CapacitorMedium); 
  capacitor1812->SetLineColor(fColorAl);
  TGeoTranslation* capacitor1812trans = new TGeoTranslation(0.0,
                                        0.5*fgkSSDStiffenerWidth+ssdstiffenerseparation
                                      - capacitor1812shape->GetDY()-fgkSSDConnectorPosition[1],0.0);
  ssdhybridcapacitormother->AddNode(capacitor1812,1,capacitor1812trans);
////////////////////////////
//Hybrid Wire
////////////////////////////
  Double_t wirex = 2.*(fgkSSDConnectorPosition[0]-0.5*fgkSSDStiffenerLength
				 - 0.5*fgkSSDConnectorLength)-fgkSSDConnectorLength
				 - fgkSSDConnectorSeparation;
  Double_t wirey = ssdstiffenerseparation+fgkSSDStiffenerWidth
				 - 2.*fgkSSDConnectorPosition[1]-fgkSSDConnectorWidth;
  Double_t ssdwireradius = TMath::Sqrt(TMath::Power(wirex,2.)
				         + TMath::Power(wirey,2));
  Double_t wireangle = TMath::ATan(wirex/wirey);
  TGeoTube *hybridwireshape = new TGeoTube("HybridWireShape", 0., 
						fgkSSDWireRadius, 0.5*ssdwireradius);
  TGeoVolume* hybridwire = new TGeoVolume("HybridWire",hybridwireshape,
                                             fSSDStiffenerHybridWireMedium); 
  hybridwire->SetLineColor(fColorPhynox);
  TGeoCombiTrans* hybridwirecombitrans[2];
  hybridwirecombitrans[0] = new TGeoCombiTrans("HybridWireCombiTrans1",
                   0.5*fgkSSDStiffenerLength-fgkSSDConnectorPosition[0]
				 + 1.5*fgkSSDConnectorLength+fgkSSDConnectorSeparation,
                   0.5*ssdwireradius-0.5*fgkSSDStiffenerWidth
				 + fgkSSDConnectorPosition[1]+0.5*fgkSSDConnectorWidth,
				   ssdstiffenershape->GetDZ()
				 + fgkSSDWireRadius+fgkSSDConnectorAlHeight+fgkSSDConnectorNiHeight,
                   new TGeoRotation("HybridWireRot1",0.,90.,0.));
  hybridwirecombitrans[1] = new TGeoCombiTrans("HybridWireCombiTrans2",
                            0.0,
                          - 0.5*fgkSSDConnectorWidth+fgkSSDWireRadius,
                            0.0,	
                            new TGeoRotation("HybridWireRot2",
                          - wireangle*TMath::RadToDeg(),0.,0.));
  TGeoHMatrix* hybridwirematrix = new TGeoHMatrix();
  hybridwirematrix->MultiplyLeft(hybridwirecombitrans[0]);
  hybridwirematrix->MultiplyLeft(hybridwirecombitrans[1]);
  ssdhybridcapacitormother->AddNode(hybridwire,1,hybridwirematrix);
  ssdhybridlist->Add(ssdhybridcapacitormother);
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  delete hybridwirecombitrans[0];
  delete hybridwirecombitrans[1];
  delete ssdchipsystemlist;
  return ssdhybridlist;
  /////////////////////////////////////////////////////////////
}
///////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCoolingBlockSystem(){
  /////////////////////////////////////////////////////////////
  // SSD Cooling Block System
  /////////////////////////////////////////////////////////////
  // SSD Cooling Block and Cooling Tube Transformations
  /////////////////////////////////////////////////////////////
  TGeoRotation* localcoolingblockrot = new TGeoRotation();
  localcoolingblockrot->SetAngles(0.,90.,0.);
  TGeoCombiTrans* localcoolingblockmatrix = 
	new TGeoCombiTrans(0.,0.5*fgkSSDCoolingBlockWidth,0.,localcoolingblockrot);
  TGeoTranslation* localcoolingblocktrans;  
  TVector3* coolingblocktransvector;
  coolingblocktransvector = new TVector3(fgkSSDModuleSensorSupportDistance
								+ fgkSSDCoolingBlockLength,
								  fgkSSDSensorLength
								- 2.*fgkSSDModuleStiffenerPosition[1]
								- fgkSSDCoolingBlockWidth);
  const Int_t kcoolingblocktransnumber = 2;
  const Int_t kcoolingblocknumber = 4;
  TGeoHMatrix* coolingblockmatrix[kcoolingblocknumber];
  TGeoHMatrix* coolingtubematrix[kcoolingblocknumber];
  TGeoRotation* localcoolingtuberot = new TGeoRotation();
  localcoolingtuberot->SetAngles(0.0,90.0,0.0);
  TGeoTranslation* localcoolingtubetrans = new TGeoTranslation();
  localcoolingtubetrans->SetTranslation(0.5*fgkSSDCoolingBlockLength,
										0.5*fgkSSDCoolingBlockWidth,
											fgkSSDCoolingBlockHoleCenter);
  TGeoCombiTrans* localcoolingtubematrix = new TGeoCombiTrans(*localcoolingtubetrans,
															  *localcoolingtuberot);
  Double_t coolingtubedistance = fgkCoolingTubeSupportRmax-fgkCoolingTubeSupportRmin;
  for(Int_t i=0; i<kcoolingblocktransnumber; i++){
	  for(Int_t j=0; j<kcoolingblocktransnumber; j++){
		localcoolingblocktrans= 
		     new TGeoTranslation(i*coolingblocktransvector->X()+2*coolingtubedistance,
  								 j*coolingblocktransvector->Y(),
								 - 0.5*(fgkSSDCoolingBlockHoleCenter
							     + fgkCoolingTubeRmax));
		coolingblockmatrix[2*i+j] = new TGeoHMatrix((*localcoolingblocktrans)
								 *(*localcoolingblockmatrix));
		coolingtubematrix[2*i+j] = new TGeoHMatrix((*localcoolingblocktrans)
								 *(*localcoolingtubematrix));
	}
  }
  /////////////////////////////////////////////////////////////
  // Virtual Volume containing CoolingBlock System   
  /////////////////////////////////////////////////////////////
  TGeoXtru* coolingsystemothershape = new TGeoXtru(2);
  const Int_t kmothervertexnumber = 16;  
  Double_t xmothervertex[kmothervertexnumber];
  Double_t ymothervertex[kmothervertexnumber];
  ///////////////////////
  // Setting the vertices 
  ///////////////////////fgkCoolingTubeSupportRmax
  xmothervertex[0] = 0.0,ymothervertex[0] = 0.0;
  xmothervertex[1] = xmothervertex[0], ymothervertex[1] = coolingblocktransvector->Y()
				   + fgkSSDCoolingBlockWidth;
  xmothervertex[2] = coolingblocktransvector->X()
				   + fgkSSDCoolingBlockLength
				   + 4*coolingtubedistance;
  ymothervertex[2] = ymothervertex[1];
  xmothervertex[3] = xmothervertex[2], ymothervertex[3] = ymothervertex[0];
  xmothervertex[4] = xmothervertex[3]-2.*coolingtubedistance-fgkSSDCoolingBlockLength;
  ymothervertex[4] = ymothervertex[0];
  xmothervertex[5] = xmothervertex[4], ymothervertex[5] = fgkSSDCoolingBlockWidth;
  xmothervertex[6] = xmothervertex[3]-coolingtubedistance; 
  ymothervertex[6] = ymothervertex[5]; 
  xmothervertex[7] = xmothervertex[6], ymothervertex[7] = ymothervertex[2]
				   - fgkSSDCoolingBlockWidth; 
  xmothervertex[8] = xmothervertex[5], ymothervertex[8] = ymothervertex[7];
  xmothervertex[9] = xmothervertex[8], ymothervertex[9] = ymothervertex[2]
				   - coolingtubedistance;
  xmothervertex[10] = fgkSSDCoolingBlockLength+2.*coolingtubedistance;
  ymothervertex[10] = ymothervertex[9];
  xmothervertex[11] = xmothervertex[10], ymothervertex[11] = ymothervertex[8];
  xmothervertex[12] = coolingtubedistance, ymothervertex[12] = ymothervertex[11];
  xmothervertex[13] = xmothervertex[12], ymothervertex[13] = fgkSSDCoolingBlockWidth;
  xmothervertex[14] = 2.*coolingtubedistance+fgkSSDCoolingBlockLength;
  ymothervertex[14] = ymothervertex[13];
  xmothervertex[15] = xmothervertex[14], ymothervertex[15] = ymothervertex[0];
  //////////////////////////////////////////////////////////
  coolingsystemothershape->DefinePolygon(kmothervertexnumber,
									xmothervertex,ymothervertex);
  coolingsystemothershape->DefineSection(0,-0.5*(fgkSSDCoolingBlockHoleCenter
											   + fgkCoolingTubeRmax));
  coolingsystemothershape->DefineSection(1, 0.5*(fgkSSDCoolingBlockHoleCenter
											   + fgkCoolingTubeRmax));
  TGeoVolume* coolingsystemother = new TGeoVolume("CoolingBlockSystem",
							  coolingsystemothershape,fSSDAir);
  /////////////////////////////////////////////////////////////
  // SSD Cooling Tube Part 
  /////////////////////////////////////////////////////////////
  TGeoTube* coolingtubeshape[fgkcoolingtubenumber];
  coolingtubeshape[0] = new TGeoTube(fgkCoolingTubeRmin,fgkCoolingTubeRmax,
										 0.5*fgkSSDCoolingBlockWidth); 
  coolingtubeshape[1] = new TGeoTube(0.0,fgkCoolingTubeRmin,
									 0.5*fgkSSDCoolingBlockWidth);
  TGeoVolume* coolingtube[fgkcoolingtubenumber];
  coolingtube[0] = new TGeoVolume("OuterCoolingTube",coolingtubeshape[0],
									fSSDCoolingTubePhynox);
  coolingtube[1] = new TGeoVolume("InnerCoolingTube",coolingtubeshape[1],
									fSSDCoolingTubeWater);
  coolingtube[0]->SetLineColor(fColorPhynox);
  coolingtube[1]->SetLineColor(fColorWater);
  TGeoVolume* ssdcoolingblock = GetSSDCoolingBlock(30);
  /////////////////////////////////////////////////////////////
  // Adding Cooling block to mother volume
  /////////////////////////////////////////////////////////////
   for(Int_t i=0; i<kcoolingblocknumber; i++){ 
	coolingsystemother->AddNode(ssdcoolingblock,i+1,coolingblockmatrix[i]);
	coolingsystemother->AddNode(coolingtube[0],i+1,coolingtubematrix[i]);
	coolingsystemother->AddNode(coolingtube[1],i+1,coolingtubematrix[i]);
  }
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
	delete coolingblocktransvector;
    delete localcoolingblocktrans;
	delete localcoolingblockrot;
	delete localcoolingblockmatrix;
	delete localcoolingtubetrans;
	delete localcoolingtuberot;
  /////////////////////////////////////////////////////////////
  // Checking overlaps	
  /////////////////////////////////////////////////////////////
	//coolingsystemother->CheckOverlaps(0.01);
  /////////////////////////////////////////////////////////////
	return coolingsystemother;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDStiffenerFlex()const{
  /////////////////////////////////////////////////////////////
  // SSD Flex
  /////////////////////////////////////////////////////////////
  const Int_t kssdflexlayernumber = 2;
  TGeoXtru* ssdflexshape[kssdflexlayernumber];
  for(Int_t i=0; i<kssdflexlayernumber; i++) ssdflexshape[i] = new TGeoXtru(2);
  const Int_t kmothervertexnumber = 17; 
  Double_t xmothervertex[kmothervertexnumber];
  Double_t ymothervertex[kmothervertexnumber];
  /////////////////////////////////////////////
  // Auxiliary variables for vertex positioning
  /////////////////////////////////////////////
  const Int_t kssdflexboxnumber = 5;
  Double_t ssdflexboxlength[kssdflexboxnumber];
  ssdflexboxlength[0] = 0.5*(fgkSSDChipNumber
					  *	fgkSSDChipLength+(fgkSSDChipNumber-1)
					  *	fgkSSDChipSeparationLength
					  - fgkSSDModuleSensorSupportDistance-fgkSSDFlexHoleLength)
					  - (fgkSSDFlexLength[0]-fgkSSDFlexLength[1]);
  ssdflexboxlength[1] = fgkSSDModuleSensorSupportDistance+fgkSSDFlexHoleLength;
  ssdflexboxlength[2] = 0.5*(fgkSSDModuleSensorSupportDistance
					  -	fgkSSDFlexHoleLength-fgkSSDFlexHoleWidth); 	
  ssdflexboxlength[3] = fgkSSDFlexHoleWidth;	
  ssdflexboxlength[4] = fgkSSDFlexLength[1]-ssdflexboxlength[0]
					  -	ssdflexboxlength[1];
  Double_t ssdflexboxwidth[kssdflexboxnumber];
  ssdflexboxwidth[0] = fgkSSDFlexWidth[0];
  ssdflexboxwidth[1] = fgkSSDFlexWidth[0]-fgkSSDFlexHoleWidth;
  ssdflexboxwidth[2] = fgkSSDFlexHoleWidth;
  ssdflexboxwidth[3] = ssdflexboxwidth[2]-fgkSSDFlexHoleLength;
  ssdflexboxwidth[4] = fgkSSDFlexWidth[0];
  ///////////////////////
  // Setting the vertices 
  ///////////////////////
  xmothervertex[0]  = 0.0;
  xmothervertex[1]  = xmothervertex[0];
  xmothervertex[2]  = fgkSSDFlexLength[0]-fgkSSDFlexLength[1];
  xmothervertex[3]  = xmothervertex[2]+ssdflexboxlength[0]+ssdflexboxlength[1]
					+ ssdflexboxlength[4];
  xmothervertex[4]  = xmothervertex[3];
  xmothervertex[5]  = xmothervertex[4]-ssdflexboxlength[4];
  xmothervertex[6]  = xmothervertex[5];
  xmothervertex[7]  = xmothervertex[6]-fgkSSDFlexHoleLength;
  xmothervertex[8]  = xmothervertex[7];
  xmothervertex[9]  = xmothervertex[8]-ssdflexboxlength[2];
  xmothervertex[10] = xmothervertex[9]; 
  xmothervertex[11] = xmothervertex[10]-ssdflexboxlength[3];
  xmothervertex[12] = xmothervertex[11];
  xmothervertex[13] = xmothervertex[12]-ssdflexboxlength[2];
  xmothervertex[14] = xmothervertex[13];
  xmothervertex[15] = xmothervertex[14]-fgkSSDFlexHoleLength;
  xmothervertex[16] = xmothervertex[15];
  ymothervertex[0]  = 0.0;
  ymothervertex[1]  = fgkSSDFlexWidth[1];
  ymothervertex[2]  = fgkSSDFlexWidth[0];
  ymothervertex[3]  = ymothervertex[2];
  ymothervertex[4]  = ymothervertex[0];
  ymothervertex[5]  = ymothervertex[4];
  ymothervertex[6]  = ssdflexboxwidth[2];
  ymothervertex[7]  = ymothervertex[6];
  ymothervertex[8]  = ymothervertex[0];
  ymothervertex[9]  = ymothervertex[8];
  ymothervertex[10] = ssdflexboxwidth[2]-ssdflexboxwidth[3];
  ymothervertex[11] = ymothervertex[10];
  ymothervertex[12] = ymothervertex[0];
  ymothervertex[13] = ymothervertex[12];
  ymothervertex[14] = ymothervertex[7];
  ymothervertex[15] = ymothervertex[14];
  ymothervertex[16] = ymothervertex[0];
  /////////////////////////////////////////////////////////////
  // First Mother Volume containing SSDFlex
  /////////////////////////////////////////////////////////////
  TGeoXtru* ssdflexmothershape = new TGeoXtru(2);
  ssdflexmothershape->DefinePolygon(kmothervertexnumber,xmothervertex,
								    ymothervertex);
  ssdflexmothershape->DefineSection(0,-1.5*fgkSSDFlexHeight[0]-2.*fgkSSDFlexHeight[1]);
  ssdflexmothershape->DefineSection(1, 0.5*fgkSSDFlexHeight[0]);
  TGeoVolume* ssdflexmother = new TGeoVolume("SSDFlexMother",ssdflexmothershape,
											 fSSDAir);
  /////////////////////////////////////////////////////////////
  // SSDFlex Layer Shapes
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i<kssdflexlayernumber; i++){
	ssdflexshape[i]->DefinePolygon(kmothervertexnumber,xmothervertex,
								   ymothervertex);
    ssdflexshape[i]->DefineSection(0,-0.5*fgkSSDFlexHeight[i]);
	ssdflexshape[i]->DefineSection(1, 0.5*fgkSSDFlexHeight[i]);
  }
  /////////////////////////////////////
  // Setting Layers into Mother Volume
  /////////////////////////////////////
  Int_t ssdflexcolor[kssdflexlayernumber] = {fColorAl,fColorPolyhamide};
  TGeoMedium* ssdflexmed[kssdflexlayernumber] = {fSSDAlTraceFlexMedium,
												 fSSDKaptonFlexMedium};
  const char* ssdflexname[2*kssdflexlayernumber] = {"AlFlexLay1","KaptonFlexLay1",
													"AlFlexLay2","KaptonFlexLay2"};
  TGeoVolume* ssdflex[2*kssdflexlayernumber];
  TGeoTranslation* ssdflextrans[2*kssdflexlayernumber];
  for(Int_t i=0; i<2*kssdflexlayernumber; i++){
	ssdflex[i] = new TGeoVolume(ssdflexname[i],
								i%2==0 ? ssdflexshape[0] : ssdflexshape[1],
								i%2==0 ? ssdflexmed[0]   : ssdflexmed[1]);
	ssdflex[i]->SetLineColor(i%2==0 ? ssdflexcolor[0] : ssdflexcolor[1]);
    ssdflextrans[i]  = new TGeoTranslation(0.,0.,-0.5*i*(fgkSSDFlexHeight[0]
					 +					   fgkSSDFlexHeight[1])); 
    ssdflexmother->AddNode(ssdflex[i],1,ssdflextrans[i]);
  }
  //ssdflexmother->CheckOverlaps(0.01);
  return ssdflexmother;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDEndFlex(){
  /////////////////////////////////////////////////////////////
  // Method generating SSD End Flex   
  /////////////////////////////////////////
  Double_t ssdflexradiusmax = (fgkSSDFlexLength[3]-fgkSSDFlexLength[2])
							/  TMath::Tan(fgkSSDFlexAngle*TMath::DegToRad());
  Double_t ssdflexboxlength = fgkSSDFlexFullLength-2.*fgkSSDFlexAngle
					        * TMath::DegToRad()*ssdflexradiusmax
					        - fgkSSDFlexLength[2]-TMath::Pi()
					        * fgkSSDStiffenerHeight-fgkSSDFlexLength[0];
  const Int_t knedges = 20;  
  const Int_t karcnumber = 2;
  TVector3* vertexposition[karcnumber*(knedges+1)];
  Double_t deltangle[karcnumber] = {2.*fgkSSDFlexAngle/knedges,180.0/knedges}; 
  Double_t angle[karcnumber] = {90.-2.*fgkSSDFlexAngle,180.0};	
  Double_t radius[karcnumber] = {ssdflexradiusmax-fgkSSDStiffenerHeight,fgkSSDStiffenerHeight};
  Double_t referenceangle[karcnumber] = {-180.0*TMath::DegToRad(),
										 - 90.0*TMath::DegToRad()};
  TVector3* referencetrans[karcnumber];
  referencetrans[0] = new TVector3(ssdflexboxlength*CosD(2.*fgkSSDFlexAngle)
					+			   radius[0]*SinD(2.*fgkSSDFlexAngle),
								   radius[0]);
  referencetrans[1] = new TVector3(referencetrans[0]->X()
					+              fgkSSDFlexLength[2],
     -              fgkSSDStiffenerHeight);
for(Int_t i=0; i<karcnumber; i++){
	for(Int_t j=0; j<knedges+1; j++){
		vertexposition[j+i*(knedges+1)] = new TVector3(radius[i]*CosD(angle[i]),
										               radius[i]*SinD(angle[i]));
		angle[i] +=  deltangle[i]*(1.0-2.0*i);
	} 	
  }
  ///////////////////////
  // Setting the vertices 
  ///////////////////////
  const Int_t kendflexlayernumber = 4;
  const Int_t kendflexvertexnumber = 4*(knedges+1)+2;
  TVector3** vertex[kendflexlayernumber];
  for(Int_t i=0; i<kendflexlayernumber; i++) 
					vertex[i] = new TVector3*[kendflexvertexnumber];
  TVector3* transvector[kendflexlayernumber+1];
  TVector3* deltatransvector = new TVector3();	
  for(Int_t i=0; i<kendflexlayernumber+1; i++) transvector[i] = new TVector3();	
  transvector[0]->SetXYZ(0.0,ssdflexboxlength*SinD(2.*fgkSSDFlexAngle)
				+		 radius[0]*(1.0-CosD(2.*fgkSSDFlexAngle)),0.0);
  for(Int_t i=1; i<kendflexlayernumber+1; i++){ 	
	deltatransvector->SetXYZ((i%2!=0?fgkSSDFlexHeight[0]:fgkSSDFlexHeight[1])
					*		  CosD(fgkSSDFlexAngle),
							  (i%2!=0?fgkSSDFlexHeight[0]:fgkSSDFlexHeight[1])
					*         SinD(fgkSSDFlexAngle),0.0);   
	*transvector[i] = *transvector[i-1]+*deltatransvector;
  }
  Double_t ratioradius[karcnumber][kendflexlayernumber+1];
  ratioradius[0][0] = 1., ratioradius[1][0] = 1.;
  for(Int_t i=0; i<karcnumber; i++){
	for(Int_t j=1; j<kendflexlayernumber+1; j++){
		ratioradius[i][j] = ratioradius[i][j-1]-TMath::Power(-1.0,i)
						  * (j%2!=0?fgkSSDFlexHeight[0]:fgkSSDFlexHeight[1])
					      /radius[i];
	}
  }
  for(Int_t i=0; i<kendflexlayernumber; i++){
	vertex[i][0] = new TVector3(transvector[i]->X(),transvector[i]->Y());
	vertex[i][1] = new TVector3(transvector[i+1]->X(),transvector[i+1]->Y());
	for(Int_t j=0; j<karcnumber*(knedges+1); j++){
		if(j<(knedges+1)){
			vertex[i][j+2] = new TVector3(vertexposition[j]->X()*ratioradius[0][i+1],
										  vertexposition[j]->Y()*ratioradius[0][i+1]);
			vertex[i][j+2]->RotateZ(referenceangle[0]);
			*vertex[i][j+2] += *referencetrans[0];
			vertex[i][4*(knedges+1)-j+1] = 
							 new TVector3(vertexposition[j]->X()*ratioradius[0][i],
										  vertexposition[j]->Y()*ratioradius[0][i]);
			vertex[i][4*(knedges+1)-j+1]->RotateZ(referenceangle[0]);
			*vertex[i][4*(knedges+1)-j+1] += *referencetrans[0];
		}
		else{
		
			vertex[i][j+2] = new TVector3(vertexposition[j]->X()*ratioradius[1][i+1],
										  vertexposition[j]->Y()*ratioradius[1][i+1]);
			vertex[i][j+2]->RotateZ(referenceangle[1]);
			*vertex[i][j+2] += *referencetrans[1];
			vertex[i][4*(knedges+1)-j+1] = 
							 new TVector3(vertexposition[j]->X()*ratioradius[1][i],
										  vertexposition[j]->Y()*ratioradius[1][i]);
			vertex[i][4*(knedges+1)-j+1]->RotateZ(referenceangle[1]);
			*vertex[i][4*(knedges+1)-j+1] += *referencetrans[1];
	   }
	}
  }
  /////////////////////////////////////////////////////////////
  // First Mother Volume containing SSDEndFlex
  /////////////////////////////////////////////////////////////
  TGeoXtru* ssdendflexmothershape = new TGeoXtru(2);
  Double_t xmothervertex[kendflexvertexnumber];
  Double_t ymothervertex[kendflexvertexnumber];
  xmothervertex[0] = vertex[0][0]->X();	
  ymothervertex[0] = vertex[0][0]->Y();
  for(Int_t i=1; i<kendflexvertexnumber; i++){
	if(i<2*(knedges+1)+2){
		xmothervertex[i] = vertex[3][i]->X();
		ymothervertex[i] = vertex[3][i]->Y();
	}
	else{
		xmothervertex[i] = vertex[0][i]->X();
		ymothervertex[i] = vertex[0][i]->Y();
	}
  }
  ssdendflexmothershape->DefinePolygon(kendflexvertexnumber,
									   xmothervertex,ymothervertex);
  ssdendflexmothershape->DefineSection(0,-0.5*fgkSSDFlexWidth[0]);
  ssdendflexmothershape->DefineSection(1, 0.5*fgkSSDFlexWidth[0]);
  TGeoVolume* ssdendflexmother = new TGeoVolume("SSDEndFlexMother",
								 ssdendflexmothershape,fSSDAir);	
  //////////////////////////////////////
  // End Flex TGeoXtru Layer Definition 
  //////////////////////////////////////
  TGeoXtru* ssdendflexshape[kendflexlayernumber];
  TGeoVolume* ssdendflex[kendflexlayernumber];
  for(Int_t i=0; i<kendflexlayernumber; i++) ssdendflexshape[i] = new TGeoXtru(2);
  Double_t xvertex[kendflexlayernumber][kendflexvertexnumber];
  Double_t yvertex[kendflexlayernumber][kendflexvertexnumber];
  Int_t ssdendflexcolor[kendflexlayernumber] = {fColorAl,fColorPolyhamide};
  TGeoMedium* ssdendflexmed[kendflexlayernumber] = {fSSDAlTraceFlexMedium,
													fSSDKaptonFlexMedium};
  const char* ssdendflexname[kendflexlayernumber] = {"AlEndFlexLay1","KaptonEndFlexLay1",
													 "AlEndFlexLay2","KaptonEndFlexLay2"};
  for(Int_t i=0; i<kendflexlayernumber; i++){
	for(Int_t j=0; j<4*(knedges+1)+2; j++){
		xvertex[i][j] = vertex[i][j]->X();
		yvertex[i][j] = vertex[i][j]->Y();
	}
  ssdendflexshape[i]->DefinePolygon(kendflexvertexnumber,xvertex[i],yvertex[i]);
  ssdendflexshape[i]->DefineSection(0,-0.5*fgkSSDFlexWidth[0]);
  ssdendflexshape[i]->DefineSection(1, 0.5*fgkSSDFlexWidth[0]);
  ssdendflex[i] = new TGeoVolume(ssdendflexname[i],ssdendflexshape[i],
								 i%2==0 ? ssdendflexmed[0] : ssdendflexmed[1]);
  ssdendflex[i]->SetLineColor(i%2==0 ? ssdendflexcolor[0] : ssdendflexcolor[1]);
  ssdendflexmother->AddNode(ssdendflex[i],1);
  }
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i<karcnumber*(knedges+1); i++) delete vertexposition[i];
  for(Int_t i=0; i<karcnumber; i++) delete referencetrans[i];
  for(Int_t i=0; i<kendflexlayernumber; i++){
	for(Int_t j=0; j<kendflexvertexnumber; j++) delete vertex[i][j];
	delete [] vertex[i];
  }
  for(Int_t i=0; i<kendflexlayernumber+1; i++) delete transvector[i];	
  delete deltatransvector;
  /////////////////////////////////////////////////////////////
  //ssdendflexmother->CheckOverlaps(0.01);
  return ssdendflexmother;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDMountingBlock(){
  /////////////////////////////////////////////////////////////
  // Method generating the Mounting Block   
  /////////////////////////////////////////////////////////////  
  // Mounting Block Boxes Shapes
  ///////////////////////////////////////
  const Int_t kmountingblockboxnumber = 3;
  TGeoBBox* mountingblockboxshape[kmountingblockboxnumber];
  mountingblockboxshape[0] = new TGeoBBox("MountingBlockBoxShape0",
							0.25*(fgkSSDMountingBlockLength[0]
						-	fgkSSDMountingBlockLength[1]),
							0.5*fgkSSDMountingBlockWidth,
							0.5*fgkSSDMountingBlockHeight[0]);
  mountingblockboxshape[1] = new TGeoBBox("MountingBlockBoxShape1",
							0.25*(fgkSSDMountingBlockLength[1]
						-	fgkSSDMountingBlockLength[2]),
							0.5*fgkSSDMountingBlockWidth,
							0.5*(fgkSSDMountingBlockHeight[1]
						-	fgkSSDMountingBlockHeight[3]));
  mountingblockboxshape[2] = new TGeoBBox("MountingBlockBoxShape2",
							0.5*fgkSSDMountingBlockLength[2],
							0.5*fgkSSDMountingBlockWidth,
							0.5*(fgkSSDMountingBlockHeight[2]
						-	fgkSSDMountingBlockHeight[3]));
  TGeoTranslation* mountingblockboxtrans[kmountingblockboxnumber+2];
  mountingblockboxtrans[0] = new TGeoTranslation("MountingBlockBoxTrans0",0.,0.,0.);
  mountingblockboxtrans[1] = new TGeoTranslation("MountingBlockBoxTrans1",
							mountingblockboxshape[0]->GetDX()
						+	mountingblockboxshape[1]->GetDX(),
							0.0,
							mountingblockboxshape[1]->GetDZ()
						-	mountingblockboxshape[0]->GetDZ()
						+	fgkSSDMountingBlockHeight[3]);
  mountingblockboxtrans[2] = new TGeoTranslation("MountingBlockBoxTrans2",
							mountingblockboxshape[0]->GetDX()
						+	2.*mountingblockboxshape[1]->GetDX()
						+	mountingblockboxshape[2]->GetDX(),
							0.0,
							mountingblockboxshape[2]->GetDZ()
						-	mountingblockboxshape[0]->GetDZ()
						+	fgkSSDMountingBlockHeight[3]);
  mountingblockboxtrans[3] = new TGeoTranslation("MountingBlockBoxTrans3",
							mountingblockboxshape[0]->GetDX()
						+	mountingblockboxshape[1]->GetDX()
						+	2.*(mountingblockboxshape[1]->GetDX()
						+	mountingblockboxshape[2]->GetDX()),
							0.0,
							mountingblockboxshape[1]->GetDZ()
						-	mountingblockboxshape[0]->GetDZ()
						+	fgkSSDMountingBlockHeight[3]);
  mountingblockboxtrans[4] = new TGeoTranslation("MountingBlockBoxTrans4",
							2.*(mountingblockboxshape[0]->GetDX()
						+	2.*mountingblockboxshape[1]->GetDX()
						+   mountingblockboxshape[2]->GetDX()),
							0.0,
							0.0);
  for(Int_t i=0; i<kmountingblockboxnumber+2; i++) 
									mountingblockboxtrans[i]->RegisterYourself();
  ///////////////////////////////////////
  // Mounting Block Trapezoid Hole Shapes
  ///////////////////////////////////////
  const Int_t kholetrapezoidvertexnumber = 4;
  TVector3* holetrapezoidvertex[kholetrapezoidvertexnumber];
  holetrapezoidvertex[0] = new TVector3();
  holetrapezoidvertex[1] = new TVector3(fgkSSDMountingBlockHoleTrapezoidHeight);
  holetrapezoidvertex[2] = new TVector3(*holetrapezoidvertex[0]);
  holetrapezoidvertex[3] = new TVector3(*holetrapezoidvertex[1]);
  Double_t holetrapezoidwidth[2] = {fgkSSDMountingBlockHoleTrapezoidUpBasis
						+	2.*mountingblockboxshape[1]->GetDX()
						*	TMath::Tan(fgkSSDMountingBlockHoleTrapezoidAngle
						*	TMath::DegToRad()),
							fgkSSDMountingBlockHoleTrapezoidUpBasis}; 
  GetArbShape(holetrapezoidvertex,
							holetrapezoidwidth,
							2.*mountingblockboxshape[1]->GetDX(),
							"HoleTrapezoidShape");
  TGeoRotation* holetrapezoidshaperot[2];
  holetrapezoidshaperot[0] = new TGeoRotation("HoleTrapezoidShapeRot0",
							90.,-90.,-90.);
  holetrapezoidshaperot[1] = new TGeoRotation("HoleTrapezoidShapeRot1",
							-180.,0.,0.);
  TGeoCombiTrans* holetrapezoidshapecombitrans = 
							new TGeoCombiTrans("HoleTrapezoidShapeCombiTrans",
							mountingblockboxshape[0]->GetDX()
						+	3.*mountingblockboxshape[1]->GetDX()
						+	2.*mountingblockboxshape[2]->GetDX(),
							0.5*fgkSSDMountingBlockWidth,
						-	fgkSSDMountingBlockHoleTrapezoidHeight
						+	2.*mountingblockboxshape[1]->GetDZ()
						-	mountingblockboxshape[0]->GetDZ()
						+	fgkSSDMountingBlockHeight[3],
							new TGeoRotation((*holetrapezoidshaperot[1])
						*	(*holetrapezoidshaperot[0])));
  holetrapezoidshapecombitrans->RegisterYourself();
  ///////////////////////////////////
  // Mounting Block Screw Hole Shapes
  ///////////////////////////////////
  const Int_t kmountingblocktubenumber = 2;
  TGeoTube* mountingblocktubeshape[kmountingblocktubenumber];
  mountingblocktubeshape[0] = new TGeoTube("MountingBlockTubeShape0",0.0,
							fgkSSDMountingBlockHoleRadius,
							mountingblockboxshape[0]->GetDZ());
  mountingblocktubeshape[1] = new TGeoTube("MountingBlockTubeShape1",0.0,
							fgkSSDMountingBlockHoleRadius,
							mountingblockboxshape[2]->GetDZ());
  TGeoTranslation* mountingblocktubetrans[2*kmountingblocktubenumber];
  mountingblocktubetrans[0] = new TGeoTranslation("MountingBlockTubeTrans0",
						-	0.5*(fgkSSDMountingBlockLength[0]
						-	fgkSSDMountingBlockHoleTubeLength[0]),
							0.5*fgkSSDMountingBlockWidth
						-	fgkSSDMountingBlockHoleTubeWidth[0],0.);
  mountingblocktubetrans[1] = new TGeoTranslation("MountingBlockTubeTrans1",
						-	0.5*(fgkSSDMountingBlockLength[0]
						-	fgkSSDMountingBlockHoleTubeLength[0])
						+	fgkSSDMountingBlockHoleTubeLength[0],
						-	0.5*fgkSSDMountingBlockWidth
						+	fgkSSDMountingBlockHoleTubeWidth[0],
							0.);
  mountingblocktubetrans[2] = new TGeoTranslation("MountingBlockTubeTrans2",
						-	mountingblockboxshape[0]->GetDX()
						+	0.5*fgkSSDMountingBlockLength[0]
						-	fgkSSDMountingBlockHoleTubeLength[1],
							0.5*fgkSSDMountingBlockWidth
						-	fgkSSDMountingBlockHoleTubeWidth[0],
							mountingblockboxshape[2]->GetDZ()
						-	mountingblockboxshape[0]->GetDZ()
						+	fgkSSDMountingBlockHeight[3]);
  mountingblocktubetrans[3] = new TGeoTranslation("MountingBlockTubeTrans3",
						-	mountingblockboxshape[0]->GetDX()
						+	0.5*fgkSSDMountingBlockLength[0],
						-	0.5*fgkSSDMountingBlockWidth
						+	fgkSSDMountingBlockHoleTubeWidth[1],
							mountingblockboxshape[2]->GetDZ()
						-	mountingblockboxshape[0]->GetDZ()
						+	fgkSSDMountingBlockHeight[3]);
  for(Int_t i=0; i<2*kmountingblocktubenumber; i++) 
						mountingblocktubetrans[i]->RegisterYourself();
						new TGeoCompositeShape("MountingBlockMainShape",
						"MountingBlockBoxShape0:MountingBlockBoxTrans0+"
						"MountingBlockBoxShape1:MountingBlockBoxTrans1+"
						"MountingBlockBoxShape2:MountingBlockBoxTrans2+"
						"MountingBlockBoxShape1:MountingBlockBoxTrans3+"
						"MountingBlockBoxShape0:MountingBlockBoxTrans4");
  ////////////////////////////////////////////
  // Mounting Block Screw Composite Hole Shapes
  ////////////////////////////////////////////
  const Int_t kmountingblockholetubesegnumber = 4;
						new TGeoTubeSeg("MountingBlockHoleTubeSegShape",
						0.0,
						fgkSSDMountingBlockScrewHoleRadius[0],
						0.5*fgkSSDMountingBlockScrewHoleHeigth,-90.,180.);
  TGeoCombiTrans* mountingblockholetubesegcombitrans[kmountingblockholetubesegnumber];
  char* mountingblockholetubesegcombitransname[kmountingblockholetubesegnumber] = 
					{	"MountingBlockHoleTubeSegCombiTrans0",
						"MountingBlockHoleTubeSegCombiTrans1",
						"MountingBlockHoleTubeSegCombiTrans2",
						"MountingBlockHoleTubeSegCombiTrans3"};
  for(Int_t i=0; i<kmountingblockholetubesegnumber; i++){
    mountingblockholetubesegcombitrans[i] =
      new TGeoCombiTrans(mountingblockholetubesegcombitransname[i],
						0.5*fgkSSDMountingBlockScrewHoleEdge*TMath::Sqrt(2)
					*	TMath::Cos(45*(2*i+1)*TMath::DegToRad()),
						0.5*fgkSSDMountingBlockScrewHoleEdge*TMath::Sqrt(2)
					*	TMath::Sin(45*(2*i+1)*TMath::DegToRad()),
						0.0,
						new TGeoRotation("",90.*i,0.,0.));
    mountingblockholetubesegcombitrans[i]->RegisterYourself();
  }
  TGeoBBox* mountingblockholeboxshape = 
						new TGeoBBox("MountingBlockHoleBoxShape",
						0.5*fgkSSDMountingBlockScrewHoleEdge,
						0.5*fgkSSDMountingBlockScrewHoleEdge,
						0.5*fgkSSDMountingBlockScrewHoleHeigth);
  TGeoCompositeShape* mountingblockscrewhole[2];
  mountingblockscrewhole[0] = 
			new TGeoCompositeShape("MountingBlockScrewHole0",
			"MountingBlockHoleTubeSegShape:MountingBlockHoleTubeSegCombiTrans0+"
			"MountingBlockHoleTubeSegShape:MountingBlockHoleTubeSegCombiTrans1+"
			"MountingBlockHoleTubeSegShape:MountingBlockHoleTubeSegCombiTrans2+"
			"MountingBlockHoleTubeSegShape:MountingBlockHoleTubeSegCombiTrans3+"
			"MountingBlockHoleBoxShape");
			new TGeoTubeSeg("MountingBlockLowerHoleTubeSegShape",
						0.0,
						fgkSSDMountingBlockScrewHoleRadius[1],
						0.5*(fgkSSDMountingBlockHoleTubeWidth[0]
					-	fgkSSDMountingBlockScrewHoleHeigth
					-	fgkSSDMountingBlockHeight[3]),0.,90.); 
  TGeoCombiTrans* mountingblocklowerholetubesegcombitrans[kmountingblockholetubesegnumber];
  char* mountingblocklowerholetubesegcombitransname[kmountingblockholetubesegnumber] =
					{	"MountingBlockLowerHoleTubeSegCombiTrans0",
						"MountingBlockLowerHoleTubeSegCombiTrans1",
						"MountingBlockLowerHoleTubeSegCombiTrans2",
						"MountingBlockLowerHoleTubeSegCombiTrans3"};
  for(Int_t i=0; i<kmountingblockholetubesegnumber; i++){
    mountingblocklowerholetubesegcombitrans[i] =
			new TGeoCombiTrans(mountingblocklowerholetubesegcombitransname[i],
						0.5*(fgkSSDMountingBlockScrewHoleEdge
					-	2.*fgkSSDMountingBlockScrewHoleRadius[1])
					*	TMath::Sqrt(2)*TMath::Cos(45*(2*i+1)*TMath::DegToRad()),
						0.5*(fgkSSDMountingBlockScrewHoleEdge
					-	2.0*fgkSSDMountingBlockScrewHoleRadius[1])
					*	TMath::Sqrt(2)*TMath::Sin(45*(2*i+1)*TMath::DegToRad()),0.,
						new TGeoRotation("",90.*i,0.,0.));
					mountingblocklowerholetubesegcombitrans[i]->RegisterYourself();
  }
  Double_t fgkSSDMountingBlockLowerScrewHoleEdge = fgkSSDMountingBlockScrewHoleEdge
					-	2.*fgkSSDMountingBlockScrewHoleRadius[1];
  TGeoBBox* mountingblocklowerholeboxshape[2];
  mountingblocklowerholeboxshape[0] = 
			new TGeoBBox("MountingBlockLowerHoleBoxShape0",
						0.5*fgkSSDMountingBlockLowerScrewHoleEdge,
						0.5*fgkSSDMountingBlockLowerScrewHoleEdge,
						0.5*(fgkSSDMountingBlockHoleTubeWidth[0]
					-	fgkSSDMountingBlockScrewHoleHeigth
					-	fgkSSDMountingBlockHeight[3]));
  mountingblocklowerholeboxshape[1] = 
			new TGeoBBox("MountingBlockLowerHoleBoxShape1",
						0.5*fgkSSDMountingBlockLowerScrewHoleEdge,
						0.5*fgkSSDMountingBlockScrewHoleRadius[1],
						0.5*(fgkSSDMountingBlockHoleTubeWidth[0]
					-	fgkSSDMountingBlockScrewHoleHeigth
					-	fgkSSDMountingBlockHeight[3]));
  TGeoCombiTrans* mountingblocklowerholeBoxcombitrans[kmountingblockholetubesegnumber];
  char* mountingBlockLowerHoleBoxCombiTransName[kmountingblockholetubesegnumber] = 
					{	"MountingBlockLowerHoleBoxCombiTrans0",
						"MountingBlockLowerHoleBoxCombiTrans1",
						"MountingBlockLowerHoleBoxCombiTrans2",
						"MountingBlockLowerHoleBoxCombiTrans3"};
  for(Int_t i=0; i<kmountingblockholetubesegnumber; i++){
    mountingblocklowerholeBoxcombitrans[i] =
			new TGeoCombiTrans(mountingBlockLowerHoleBoxCombiTransName[i],
						0.5*(fgkSSDMountingBlockLowerScrewHoleEdge
					+	fgkSSDMountingBlockScrewHoleRadius[1])
					*	TMath::Cos(90*(i+1)*TMath::DegToRad()),
						0.5*(fgkSSDMountingBlockLowerScrewHoleEdge
					+	fgkSSDMountingBlockScrewHoleRadius[1])
					*	TMath::Sin(90*(i+1)*TMath::DegToRad()),0.,
						new TGeoRotation("",90.*i,0.,0.));
    mountingblocklowerholeBoxcombitrans[i]->RegisterYourself();
  }
  mountingblockscrewhole[1] = new TGeoCompositeShape("MountingBlockScrewHole1",
	"MountingBlockLowerHoleTubeSegShape:MountingBlockLowerHoleTubeSegCombiTrans0+"
	"MountingBlockLowerHoleTubeSegShape:MountingBlockLowerHoleTubeSegCombiTrans1+"
	"MountingBlockLowerHoleTubeSegShape:MountingBlockLowerHoleTubeSegCombiTrans2+"
	"MountingBlockLowerHoleTubeSegShape:MountingBlockLowerHoleTubeSegCombiTrans3+"
	"MountingBlockLowerHoleBoxShape0+"
	"MountingBlockLowerHoleBoxShape1:MountingBlockLowerHoleBoxCombiTrans0+"
	"MountingBlockLowerHoleBoxShape1:MountingBlockLowerHoleBoxCombiTrans1+"
	"MountingBlockLowerHoleBoxShape1:MountingBlockLowerHoleBoxCombiTrans2+"
	"MountingBlockLowerHoleBoxShape1:MountingBlockLowerHoleBoxCombiTrans3");
  TGeoTranslation* mountingblockscrewhole1trans = 
			new TGeoTranslation("MountingBlockScrewHole1Trans",0.,0.,
					-	mountingblocklowerholeboxshape[0]->GetDZ()
					-	mountingblockholeboxshape->GetDZ());
  mountingblockscrewhole1trans->RegisterYourself();
			new TGeoCompositeShape("MountingBlockHole",
	"MountingBlockScrewHole0+MountingBlockScrewHole1:MountingBlockScrewHole1Trans");
  TGeoTranslation* mountingblockholetrans = new TGeoTranslation("MountingBlockHoleTrans",
						0.5*fgkSSDMountingBlockLength[0]
					-	mountingblockboxshape[0]->GetDZ(),
						0.0,
						2.*mountingblockboxshape[2]->GetDZ()
					-	mountingblockboxshape[0]->GetDZ()
					+	fgkSSDMountingBlockHeight[3]
					-	mountingblockholeboxshape->GetDZ());
  mountingblockholetrans->RegisterYourself();
  TGeoCompositeShape* mountingblockshape = new TGeoCompositeShape("MountingBlockShape",
			"MountingBlockMainShape-(MountingBlockTubeShape0:MountingBlockTubeTrans0+"
			"MountingBlockTubeShape0:MountingBlockTubeTrans1+"
			"MountingBlockTubeShape1:MountingBlockTubeTrans2+"
			"MountingBlockTubeShape1:MountingBlockTubeTrans3+"
			"HoleTrapezoidShape:HoleTrapezoidShapeCombiTrans+"
			"MountingBlockHole:MountingBlockHoleTrans)");
  TGeoVolume* ssdmountingblock = new TGeoVolume("SSDMountingBlock",
			mountingblockshape,fSSDMountingBlockMedium);
  return ssdmountingblock;
}
///////////////////////////////////////////////////////////////////////////////
TList* AliITSv11GeometrySSD::GetCoolingTubeList()const{
  /////////////////////////////////////////////////////////////
  // Method generating the Cooling Tube 
  /////////////////////////////////////////////////////////////  
   TGeoTube** coolingtubeshape[fgkcoolingtubenumber];
   for(Int_t i=0; i<fgkcoolingtubenumber; i++) coolingtubeshape[i] = 
												new	TGeoTube*[2];
   coolingtubeshape[0][0] = new TGeoTube(fgkCoolingTubeRmin,fgkCoolingTubeRmax,
					  0.25 * (fgkSSDSensorLength-2.*fgkSSDModuleStiffenerPosition[1]
						  -	2.*fgkSSDCoolingBlockWidth-fgkCoolingTubeSupportWidth));
   coolingtubeshape[0][1] = new TGeoTube(0.0,fgkCoolingTubeRmin,
										 coolingtubeshape[0][0]->GetDz());
   coolingtubeshape[1][0] = new TGeoTube(fgkCoolingTubeRmin,fgkCoolingTubeRmax,
										 0.5*(fgkSSDModuleStiffenerPosition[1]
						  -					  fgkSSDSensorOverlap));
   coolingtubeshape[1][1] = new TGeoTube(0.0,fgkCoolingTubeRmin,
										 coolingtubeshape[1][0]->GetDz());
   coolingtubeshape[2][0] = new TGeoTube(fgkCoolingTubeRmin,fgkCoolingTubeRmax,
										 0.5*fgkSSDModuleStiffenerPosition[1]);
   coolingtubeshape[2][1] = new TGeoTube(0.0,fgkCoolingTubeRmin,
										 coolingtubeshape[2][0]->GetDz());
   TGeoVolume** coolingtube[fgkcoolingtubenumber];
   for(Int_t i=0; i<fgkcoolingtubenumber; i++) coolingtube[i] = 
											 new TGeoVolume*[2];
   coolingtube[0][0] = new TGeoVolume("OuterCoolingTube1",coolingtubeshape[0][0],
									  fSSDCoolingTubePhynox);
   coolingtube[0][1] = new TGeoVolume("InnerCoolingTube1",coolingtubeshape[0][1],
									  fSSDCoolingTubeWater);
   coolingtube[1][0] = new TGeoVolume("OuterCoolingTube2",coolingtubeshape[1][0],
									  fSSDCoolingTubePhynox);
   coolingtube[1][1] = new TGeoVolume("InnerCoolingTube2",coolingtubeshape[1][1],
									  fSSDCoolingTubeWater);
   coolingtube[2][0] = new TGeoVolume("OuterCoolingTube3",coolingtubeshape[2][0],
									  fSSDCoolingTubePhynox);
   coolingtube[2][1] = new TGeoVolume("InnerCoolingTube3",coolingtubeshape[2][1],
									  fSSDCoolingTubeWater);
   for(Int_t i=0; i<fgkcoolingtubenumber; i++){
	coolingtube[i][0]->SetLineColor(fColorPhynox);
	coolingtube[i][1]->SetLineColor(fColorWater);
   }
  /////////////////////////////////////////////////////////////
  // Virtual Volume containing Cooling Tubes
  /////////////////////////////////////////////////////////////
  TGeoTube* virtualcoolingtubeshape[fgkcoolingtubenumber];
  for(Int_t i=0; i<fgkcoolingtubenumber; i++)
  virtualcoolingtubeshape[i] = new TGeoTube(coolingtubeshape[i][1]->GetRmin(),
											coolingtubeshape[i][0]->GetRmax(),
											coolingtubeshape[i][0]->GetDz());
  TGeoVolume* virtualcoolingtube[fgkcoolingtubenumber];
  virtualcoolingtube[0] = new TGeoVolume("CoolingTube1",virtualcoolingtubeshape[0],
									  fSSDAir);
  virtualcoolingtube[1] = new TGeoVolume("CoolingTube2",virtualcoolingtubeshape[1],
									  fSSDAir);
  virtualcoolingtube[2] = new TGeoVolume("CoolingTube3",virtualcoolingtubeshape[2],
									  fSSDAir);
  TList* coolingtubelist = new TList();
  for(Int_t i=0; i<fgkcoolingtubenumber; i++){
	virtualcoolingtube[i]->AddNode(coolingtube[i][0],1);
	virtualcoolingtube[i]->AddNode(coolingtube[i][1],1);
    coolingtubelist->Add(virtualcoolingtube[i]);
  }
  return coolingtubelist;
}
///////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDCoolingBlock(Int_t nedges){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Cooling Block    
  /////////////////////////////////////////////////////////////
  const Int_t kvertexnumber = 8;
  ///////////////////////////////////////
  // Vertex Positioning for TGeoXTru
  ///////////////////////////////////////
  TVector3** vertexposition = new TVector3*[2*kvertexnumber+nedges+1];
  vertexposition[0] = new TVector3(0.0,0.0);
  vertexposition[1] = new TVector3(0.0,fgkSSDCoolingBlockHeight[1]);
  vertexposition[2] = new TVector3(fgkSSDCoolingBlockHoleLength[1],
					  vertexposition[1]->Y());
  vertexposition[3] = new TVector3(vertexposition[2]->X(),
					  vertexposition[2]->Y()+fgkSSDCoolingBlockHeight[2]);
  vertexposition[4] = new TVector3(vertexposition[1]->X(),vertexposition[3]->Y());
  vertexposition[5] = new TVector3(vertexposition[4]->X(),
					+ vertexposition[3]->Y()+fgkSSDCoolingBlockHoleRadius[1]);
  vertexposition[6] = new TVector3(Xfrom2Points(vertexposition[5]->X(),
					  vertexposition[5]->Y(),0.5*(fgkSSDCoolingBlockLength
					- fgkSSDCoolingBlockHoleLength[0]
					- 4.*fgkSSDCoolingBlockHoleRadius[1]),
					  fgkSSDCoolingBlockHeight[0]
					- fgkSSDCoolingBlockHoleRadius[1],
					  fgkSSDCoolingBlockHeight[0]),fgkSSDCoolingBlockHeight[0]);
  vertexposition[7] = new TVector3(0.5*(fgkSSDCoolingBlockLength
					- fgkSSDCoolingBlockHoleLength[0]),
					  vertexposition[6]->Y());
  Double_t alpha = TMath::ACos(0.5*fgkSSDCoolingBlockHoleLength[0]
			   / fgkSSDCoolingBlockHoleRadius[0])*TMath::RadToDeg();
  Double_t phi = 180.-alpha;
  Double_t psi = 180.+2.*alpha;
  Double_t deltapsi = psi/nedges;
  Double_t radius = fgkSSDCoolingBlockHoleRadius[0]/CosD(0.5*deltapsi);
  TVector3* transvector = new TVector3(0.5*fgkSSDCoolingBlockLength,
						  fgkSSDCoolingBlockHoleCenter);
  for(Int_t i=0; i<nedges+1; i++){
	vertexposition[kvertexnumber+i] = new TVector3(radius*CosD(phi+i*deltapsi),
											       radius*SinD(phi+i*deltapsi));
   *vertexposition[kvertexnumber+i] += (*transvector);
  }
  Double_t param[4] = {1.0,0.0,0.0,-0.5*fgkSSDCoolingBlockLength};  
  for(Int_t i=0; i<kvertexnumber; i++)
    vertexposition[kvertexnumber+nedges+1+i] = 
						GetReflection(vertexposition[kvertexnumber-1-i],param);
  ///////////////////////////////////////////////////////////////////////
  // TGeoXTru Volume definition for Cooling Tube Support Arc Part
  ///////////////////////////////////////////////////////////////////////
  TGeoXtru* ssdcoolingblockshape = new TGeoXtru(2);	
  Double_t* xvertexpoints = new Double_t[2*kvertexnumber+nedges+1]; 
  Double_t* yvertexpoints = new Double_t[2*kvertexnumber+nedges+1];
  for(Int_t i=0; i<2*kvertexnumber+nedges+1; i++){
	xvertexpoints[i] = vertexposition[i]->X();
	yvertexpoints[i] = vertexposition[i]->Y();
  } 
  ssdcoolingblockshape->DefinePolygon(2*kvertexnumber+nedges+1,xvertexpoints,
											yvertexpoints);
  ssdcoolingblockshape->DefineSection(0,-0.5*fgkSSDCoolingBlockWidth);
  ssdcoolingblockshape->DefineSection(1,0.5*fgkSSDCoolingBlockWidth);
  TGeoVolume* ssdcoolingblock = new TGeoVolume("SSDCoolingBlock",
								          ssdcoolingblockshape,
										  fSSDAlCoolBlockMedium);
  ssdcoolingblock->SetLineColor(fColorAl);
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  delete [] vertexposition;
  delete xvertexpoints;
  delete yvertexpoints;
  /////////////////////////////////////////////////////////////
  return ssdcoolingblock;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDChipCables(Double_t SSDChipCablesHeigth, Int_t nedges){
  ///////////////////////////////////////////////////////
  const Int_t kssdchipcablesnumber    = 2;
  const Int_t kssdchipcableslaynumber = 2;
  const Int_t kvertexnumber			  = 4*(nedges+1)+4;
  Int_t ssdchipcablescolor[kssdchipcableslaynumber] = {fColorAl,fColorPolyhamide};
  Double_t ssdchipcablesradius[kssdchipcableslaynumber];
  ssdchipcablesradius[0] = 0.25*(SSDChipCablesHeigth
						 -  fgkSSDChipCablesHeight[0]
						 -  fgkSSDChipCablesHeight[1]);
  ssdchipcablesradius[1] = ssdchipcablesradius[0]-fgkSSDChipCablesHeight[0];
  Double_t ssdchipcablespiecelength[kssdchipcablesnumber];
  ssdchipcablespiecelength[0] = 0.5*(fgkSSDChipCablesWidth[0]
							  - 2.*TMath::Pi()*ssdchipcablesradius[0]
							  - ssdchipcablesradius[0]
							  - fgkSSDChipCablesWidth[1]
							  - fgkSSDChipCablesWidth[2]);
  ssdchipcablespiecelength[1] = ssdchipcablespiecelength[0]
							  - 0.5*(fgkSSDModuleStiffenerPosition[1]							  
							  +      fgkSSDChipCablesHeight[1]
							  +      fgkSSDSensorHeight);
  ///////////////////////////////////////////////////////
  // Vertex Positioning for TGeoXTrue Layer 1 and Layer 2
  ///////////////////////////////////////////////////////
  TVector3** vertexposition[kssdchipcableslaynumber];
  for(Int_t i=0; i<kssdchipcableslaynumber; i++) vertexposition[i] = 
												  new TVector3*[4*(nedges+1)+4];
  Double_t ratio[4];
  ratio[0] = ssdchipcablesradius[1]/ssdchipcablesradius[0]; 
  ratio[1] = (ssdchipcablesradius[1]-fgkSSDChipCablesHeight[1])
		   /  ssdchipcablesradius[0]; 
  ratio[2] = (ssdchipcablesradius[0]+fgkSSDChipCablesHeight[0])
		   /  ssdchipcablesradius[0];
  ratio[3] = (ssdchipcablesradius[0]+fgkSSDChipCablesHeight[0]
		   +  fgkSSDChipCablesHeight[1])
		   /  ssdchipcablesradius[0];
  Double_t phi = 180.;
  Double_t deltaphi = 180./nedges;
  Double_t angle = 0.0;
  Double_t **xvertexpoints = new Double_t*[kssdchipcableslaynumber];
  Double_t **yvertexpoints = new Double_t*[kssdchipcableslaynumber];
  for(Int_t i=0; i<kssdchipcableslaynumber;i++){
	xvertexpoints[i] = new Double_t[kvertexnumber];
	yvertexpoints[i] = new Double_t[kvertexnumber];
  }  
  TVector3* vertex = new TVector3();
  TVector3* transvector[kssdchipcableslaynumber];
  transvector[0] = new TVector3(fgkSSDChipWidth,
								SSDChipCablesHeigth-ssdchipcablesradius[0]);
  transvector[1] = new TVector3();
  TGeoXtru* ssdchipcableshape[kssdchipcableslaynumber*kssdchipcableslaynumber];
  TGeoVolume* ssdchipcable[kssdchipcableslaynumber*kssdchipcableslaynumber];
  const char* ssdchipcablename[kssdchipcableslaynumber*kssdchipcableslaynumber] = 
		{"SSDChipcableAllay1Left","SSDChipcableKaptonlay2Left",
		 "SSDChipcableAllay1Right","SSDChipcableKaptonlay2Right"};
  for(Int_t k=0; k<kssdchipcablesnumber; k++){
	transvector[1]->SetX(fgkSSDChipWidth-ssdchipcablespiecelength[k]);
	transvector[1]->SetY(ssdchipcablesradius[0]
				 +		 fgkSSDChipCablesHeight[0]
				 +		 fgkSSDChipCablesHeight[1]);  
	for(Int_t i=0; i<kssdchipcableslaynumber; i++){
		vertexposition[i][0] = new TVector3(0.,SSDChipCablesHeigth
							 - fgkSSDChipCablesHeight[0]-i*fgkSSDChipCablesHeight[1]);
		vertexposition[i][1] = new TVector3(0.,SSDChipCablesHeigth
							 - i*fgkSSDChipCablesHeight[0]);
		vertexposition[i][2*(nedges+1)+2] = 
					new TVector3(fgkSSDChipWidth+ssdchipcablesradius[0]
				+				 fgkSSDChipCablesWidth[1]
				+				 fgkSSDChipCablesWidth[2],
								((1.-i)*fgkSSDChipCablesHeight[i]
				+				 fgkSSDChipCablesHeight[1]));
        vertexposition[i][2*(nedges+1)+3] = 
					new TVector3(vertexposition[i][2*(nedges+1)+2]->X(),
								 vertexposition[i][2*(nedges+1)+2]->Y()
				-				 fgkSSDChipCablesHeight[i]);
	    for(Int_t j=0; j<nedges+1; j++){ 		
		    angle = 0.5*phi+TMath::Power(-1,i+1)*j*deltaphi;
			vertex->SetX(ssdchipcablesradius[0]*CosD(angle));
			vertex->SetY(ssdchipcablesradius[0]*SinD(angle));
			vertexposition[0][(nedges+1)*i+j+2] = 
						new TVector3(*vertex+*transvector[i]);
			vertexposition[1][(nedges+1)*i+j+2] = 
						new TVector3(vertex->X()*ratio[2*i]+transvector[i]->X(),
									 vertex->Y()*ratio[2*i]+transvector[i]->Y());
			vertexposition[0][(4-i)*(nedges+1)+4-j-1] = 
						new TVector3(*vertexposition[1][(nedges+1)*i+j+2]);
			vertexposition[1][(4-i)*(nedges+1)+4-j-1] = 
						new TVector3(vertex->X()*ratio[2*i+1]
							+			 transvector[i]->X(),
									   	 vertex->Y()*ratio[2*i+1]
							+	     	 transvector[i]->Y());
		}
	}
	for(Int_t i=0; i<kssdchipcableslaynumber; i++){
		for(Int_t j=0; j<kvertexnumber; j++){ 	
			xvertexpoints[i][j] = vertexposition[i][j]->X();
			yvertexpoints[i][j] = vertexposition[i][j]->Y();
		}
		ssdchipcableshape[kssdchipcablesnumber*k+i] = new TGeoXtru(2);
		ssdchipcableshape[kssdchipcablesnumber*k+i]->DefinePolygon(kvertexnumber,
										xvertexpoints[i],yvertexpoints[i]);
		ssdchipcableshape[kssdchipcablesnumber*k+i]->DefineSection(0,-0.5*fgkSSDChipCablesLength[1]);
		ssdchipcableshape[kssdchipcablesnumber*k+i]->DefineSection(1,+0.5*fgkSSDChipCablesLength[1]);
		ssdchipcable[kssdchipcablesnumber*k+i] = 
				new TGeoVolume(ssdchipcablename[kssdchipcablesnumber*k+i],
							   ssdchipcableshape[kssdchipcablesnumber*k+i],
							  (kssdchipcablesnumber*k+i)%2==0?
							   fSSDAlTraceChipCableMedium:fSSDKaptonChipCableMedium);
		ssdchipcable[kssdchipcablesnumber*k+i]->SetLineColor(ssdchipcablescolor[i]);
	}
	for(Int_t i=0; i<kssdchipcableslaynumber; i++)
		for(Int_t j=0; j<kvertexnumber; j++) delete vertexposition[i][j];
  }
  /////////////////////////////////////////////////////////////
  // Mother Volume definition 
  /////////////////////////////////////////////////////////////
  Double_t ssdchipseparation = fgkSSDSensorLength
							 - 2.*fgkSSDModuleStiffenerPosition[1]
							 - 2.*(fgkSSDStiffenerWidth-fgkSSDStiffenerToChipDist
							 - 0.5*fgkSSDChipWidth)-fgkSSDChipWidth;
  Double_t boxorigin[3] = {-0.5*ssdchipseparation,0.,0.5*SSDChipCablesHeigth}; 
  Double_t dx = ssdchipseparation+2.*(fgkSSDChipWidth+ssdchipcablesradius[0]
							  +fgkSSDChipCablesWidth[1]
							  +fgkSSDChipCablesWidth[2]);
  Double_t dy = fgkSSDChipCablesLength[1];
  Double_t dz = SSDChipCablesHeigth;
  TGeoBBox* ssdchipcablesmotherbox = new TGeoBBox(0.5*dx,0.5*dy,0.5*dz,boxorigin);
  TGeoVolume* ssdchipcablesmother = new TGeoVolume("SSDChipCablesMother",
			  ssdchipcablesmotherbox,fSSDAir);
  /////////////////////////////////////////////////////////////
  // Rotation and Translation Definition for positioning 
  /////////////////////////////////////////////////////////////
  TGeoRotation* ssdchipcablesrot[5];
  ssdchipcablesrot[0] = new TGeoRotation("",90.,180.,-90);
  ssdchipcablesrot[1] = new TGeoRotation("",0.0,90.0,0.0);
  ssdchipcablesrot[2] = new TGeoRotation((*ssdchipcablesrot[1])*(*ssdchipcablesrot[0]));
  ssdchipcablesrot[3] = new TGeoRotation("",180.,0.0,0.0);
  ssdchipcablesrot[4] = new TGeoRotation((*ssdchipcablesrot[3])*(*ssdchipcablesrot[2]));
  TGeoCombiTrans* ssdchipcablescombitrans = new TGeoCombiTrans(-ssdchipseparation,
														0.,0.,ssdchipcablesrot[2]);
  ssdchipcablesmother->AddNode(ssdchipcable[0],1,ssdchipcablesrot[4]);
  ssdchipcablesmother->AddNode(ssdchipcable[1],1,ssdchipcablesrot[4]);
  ssdchipcablesmother->AddNode(ssdchipcable[2],1,ssdchipcablescombitrans);
  ssdchipcablesmother->AddNode(ssdchipcable[3],1,ssdchipcablescombitrans);
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i<kssdchipcableslaynumber;i++){
	delete [] xvertexpoints[i];
	delete [] yvertexpoints[i];
  }
  for(Int_t i=0; i<kssdchipcableslaynumber; i++) delete [] vertexposition[i];
  for(Int_t i=0; i<kssdchipcableslaynumber; i++) delete transvector[i];
  delete vertex; 
  delete ssdchipcablesrot[0];
  delete ssdchipcablesrot[1];
  delete ssdchipcablesrot[3];
  /////////////////////////////////////////////////////////////
  return ssdchipcablesmother;
}
///////////////////////////////////////////////////////////////////////////////
TList* AliITSv11GeometrySSD::GetSSDChipSystem(){
  /////////////////////////////////////////////////////////////
  // SSD Chip Assembly
  /////////////////////////////////////////////////////////////
  TGeoVolume* ssdchipassembly = GetSSDChips();
  TList* ssdchipsystemlist = new TList();
  const Int_t knedges = 20;
  const Int_t kchipsystemnumber = 2;
  /////////////////////////////////////////////////////////////
  // Mother Volume containing SSDChipSystem
  /////////////////////////////////////////////////////////////
  TGeoXtru* chipsystemothershape[kchipsystemnumber];
  for(Int_t i=0; i<kchipsystemnumber; i++) chipsystemothershape[i] = new TGeoXtru(2);
  const Int_t kmothervertexnumber = 12;  
  Double_t xmothervertex[kchipsystemnumber][kmothervertexnumber];
  Double_t ymothervertex[kchipsystemnumber][kmothervertexnumber];
  Double_t ssdchipcablesradius[kchipsystemnumber];
  Double_t ssdchipseparation = fgkSSDSensorLength
			     - 2.*fgkSSDModuleStiffenerPosition[1]
			     - 2.*(fgkSSDStiffenerWidth
			     - fgkSSDStiffenerToChipDist-0.5*fgkSSDChipWidth);
  for(Int_t i=0; i<kchipsystemnumber; i++)
	ssdchipcablesradius[i] = 0.25*(fgkSSDChipCablesHeight[i+2]
						   -  fgkSSDChipCablesHeight[0]
						   -  fgkSSDChipCablesHeight[1]);
  ///////////////////////
  // Setting the vertices 
  ///////////////////////
  xmothervertex[0][0]  = -0.5*fgkSSDChipCablesLength[1];  
  xmothervertex[0][1]  = xmothervertex[0][0];  
  xmothervertex[0][2]  = (fgkSSDChipNumber-1)*(fgkSSDChipLength
					   + fgkSSDChipSeparationLength)+0.5*fgkSSDChipCablesLength[1];  
  xmothervertex[0][3]  = xmothervertex[0][2];  
  xmothervertex[0][4]  = 0.5*fgkSSDChipCablesLength[1];  
  xmothervertex[0][5]  = xmothervertex[0][4];  
  xmothervertex[0][6]  = xmothervertex[0][2]-0.5*fgkSSDChipCablesLength[1];  
  xmothervertex[0][7]  = xmothervertex[0][6]; 
  xmothervertex[0][8]  = 0.0;  
  xmothervertex[0][9]  = xmothervertex[0][8];  
  xmothervertex[0][10] = xmothervertex[0][4];  
  xmothervertex[0][11] = xmothervertex[0][10];  
  for(Int_t i=0; i<kmothervertexnumber; i++) 
	xmothervertex[1][i] = xmothervertex[0][i]; 
  for(Int_t i=0; i<kchipsystemnumber; i++){
	ymothervertex[i][0]  = -0.5*fgkSSDChipWidth-ssdchipcablesradius[i]
						 - fgkSSDChipCablesWidth[1]-fgkSSDChipCablesWidth[2];
	ymothervertex[i][1]  = ssdchipseparation-ymothervertex[i][0];
	ymothervertex[i][2]  = ymothervertex[i][1];
	ymothervertex[i][3]  = ymothervertex[i][0];
	ymothervertex[i][4]  = ymothervertex[i][0];
	ymothervertex[i][5]  = 0.5*fgkSSDChipWidth;
	ymothervertex[i][6]  = ymothervertex[i][5];
	ymothervertex[i][7]  = ssdchipseparation-0.5*fgkSSDChipWidth;
	ymothervertex[i][8]  = ymothervertex[i][7];
	ymothervertex[i][9]  = ymothervertex[i][5];
	ymothervertex[i][10] = ymothervertex[i][5];
	ymothervertex[i][11] = ymothervertex[i][4];
  }
  //////////////////////////////////////////////////////////
  TGeoVolume* chipsystemother[kchipsystemnumber];
  const char* chipsytemothername[kchipsystemnumber] = 
					{"SSDChipSytemother1","SSDChipSytemother2"};
  for(Int_t i=0; i<kchipsystemnumber; i++){
    chipsystemothershape[i]->DefinePolygon(kmothervertexnumber,
									xmothervertex[i],ymothervertex[i]);
    chipsystemothershape[i]->DefineSection(0,-fgkSSDChipCablesHeight[i+2]
										  -0.5*fgkSSDChipHeight);
    chipsystemothershape[i]->DefineSection(1,0.5*fgkSSDChipHeight);
    chipsystemother[i] = new TGeoVolume(chipsytemothername[i],
							  chipsystemothershape[i],fSSDAir);
  }
  /////////////////////////////////////////////////////////////
  // SSD Chip Cables
  /////////////////////////////////////////////////////////////
  TGeoVolume* ssdchipcables[kchipsystemnumber];
  TGeoRotation** ssdchipcablesrot[kchipsystemnumber];
  TGeoTranslation** ssdchipcablestrans[kchipsystemnumber];
  TGeoCombiTrans** ssdchipcablescombitrans[kchipsystemnumber];
  //////////////////
  for(Int_t i=0; i<kchipsystemnumber; i++){
		ssdchipcables[i] = 
		GetSSDChipCables(fgkSSDChipCablesHeight[i+2],knedges);
		ssdchipcablestrans[i] = new TGeoTranslation*[fgkSSDChipNumber];
		ssdchipcablesrot[i] = new TGeoRotation*[fgkSSDChipNumber];
		ssdchipcablescombitrans[i] = new TGeoCombiTrans*[fgkSSDChipNumber];
  }
  for(Int_t i=0; i<kchipsystemnumber; i++){
	for(Int_t j=0; j<fgkSSDChipNumber; j++){
		ssdchipcablestrans[i][j] = new TGeoTranslation();
		ssdchipcablesrot[i][j] = new TGeoRotation();
		ssdchipcablescombitrans[i][j] = new TGeoCombiTrans();
		ssdchipcablesrot[i][j]->SetAngles(-90.0,0.0,0.0);
		ssdchipcablestrans[i][j]->SetTranslation(j*(fgkSSDChipLength
						  +                fgkSSDChipSeparationLength),
											0.5*fgkSSDChipWidth,
						  -					0.5*fgkSSDChipHeight
						  -					fgkSSDChipCablesHeight[i+2]);
		ssdchipcablescombitrans[i][j]->SetRotation(*ssdchipcablesrot[i][j]);
		ssdchipcablescombitrans[i][j]->SetTranslation(*ssdchipcablestrans[i][j]);
		chipsystemother[i]->AddNode(ssdchipcables[i],j+1,ssdchipcablescombitrans[i][j]);

	}
	chipsystemother[i]->AddNode(ssdchipassembly,i+1);
       	ssdchipsystemlist->Add(chipsystemother[i]);	
  }
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i<kchipsystemnumber; i++){
	for(Int_t j=0; j<fgkSSDChipNumber; j++){
		delete ssdchipcablesrot[i][j];
		delete ssdchipcablestrans[i][j];
	}
	delete ssdchipcablesrot[i];
	delete ssdchipcablestrans[i];
  }
  /////////////////////////////////////////////////////////////
  return ssdchipsystemlist;
}
///////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDChips() const{
  /////////////////////////////////////////////////////////////
  // SSD Chip Assembly Generation    
  /////////////////////////////////////////////////////////////
  const Int_t kssdchiprownumber = 2;
  TGeoBBox* ssdchipcompshape[2];
  ssdchipcompshape[0] =  new TGeoBBox("SSDChipCompShape",
										0.5*fgkSSDChipLength,
										0.5*fgkSSDChipWidth,
										0.5*(fgkSSDChipHeight-fgkSSDChipGlueHeight));
  ssdchipcompshape[1] =  new TGeoBBox("SSDChipGlueCompShape",
										0.5*fgkSSDChipLength,
										0.5*fgkSSDChipWidth,
										0.5*fgkSSDChipGlueHeight);
  TGeoVolume* ssdchipcomp[2];
  ssdchipcomp[0] = new TGeoVolume("SSDChipComp",ssdchipcompshape[0],fSSDChipMedium);
  ssdchipcomp[1] = new TGeoVolume("SSDChipGlueComp",ssdchipcompshape[1],
								  fSSDChipGlueMedium);
  ssdchipcomp[0]->SetLineColor(fColorSilicon);  
  ssdchipcomp[1]->SetLineColor(fColorEpoxy);
  TGeoTranslation* ssdchipcomptrans[2];
  ssdchipcomptrans[0] = new TGeoTranslation(0.,0.,-ssdchipcompshape[1]->GetDZ());
  ssdchipcomptrans[1] = new TGeoTranslation(0.,0.,ssdchipcompshape[0]->GetDZ());
  /////////////////////////////////////////////////////////////
  // Virtual Volume containing SSDChip   
  /////////////////////////////////////////////////////////////
  TGeoBBox* ssdvirtualchipshape = new TGeoBBox("SSDChipShape",0.5*fgkSSDChipLength,
						  							         0.5*fgkSSDChipWidth,
													         0.5*fgkSSDChipHeight);
  TGeoVolume* ssdchip = new TGeoVolume("SSDChip",ssdvirtualchipshape,fSSDAir);  
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i<2; i++) ssdchip->AddNode(ssdchipcomp[i],1,ssdchipcomptrans[i]);
  Double_t ssdchipseparation[2] = {fgkSSDChipLength+fgkSSDChipSeparationLength,
						  fgkSSDSensorLength-2.*fgkSSDModuleStiffenerPosition[1]
				   -  2.*(fgkSSDStiffenerWidth-fgkSSDStiffenerToChipDist
				   -  0.5*fgkSSDChipWidth)};
  /////////////////////////////////////////////////////////////
  // Virtual Volume containing SSDChipAssembly   
  /////////////////////////////////////////////////////////////
  TGeoXtru* ssdchipmothershape = new TGeoXtru(2);
  const Int_t kssdmothervertexnumber = 2*fgkSSDChipNumber;  
  Double_t xmothervertex[kssdmothervertexnumber];
  Double_t ymothervertex[kssdmothervertexnumber];
  ///////////////////////
  // Setting the vertices 
  ///////////////////////
  xmothervertex[0] = -0.5*fgkSSDChipLength,ymothervertex[0] = -0.5*fgkSSDChipWidth;
  xmothervertex[1] = xmothervertex[0], ymothervertex[1] = ssdchipseparation[1]
				   - ymothervertex[0];
  xmothervertex[2] = (fgkSSDChipNumber-1)*ssdchipseparation[0]-xmothervertex[0];
  ymothervertex[2] = ymothervertex[1];
  xmothervertex[3] = xmothervertex[2], ymothervertex[3] = ymothervertex[0];
  xmothervertex[4] = ssdchipseparation[0]+xmothervertex[0];
  ymothervertex[4] = ymothervertex[0];
  xmothervertex[5] = xmothervertex[4], ymothervertex[5] = -ymothervertex[4];
  xmothervertex[6] = (fgkSSDChipNumber-1)*ssdchipseparation[0]
				   + (0.5*fgkSSDChipLength-fgkSSDChipWidth);
  ymothervertex[6] = ymothervertex[5];
  xmothervertex[7] = xmothervertex[6], ymothervertex[7] = ymothervertex[2]
				   - fgkSSDChipWidth;
  xmothervertex[8] = -0.5*fgkSSDChipLength+fgkSSDChipWidth;
  ymothervertex[8] = ymothervertex[7];
  xmothervertex[9] = -0.5*fgkSSDChipLength+fgkSSDChipWidth;
  ymothervertex[9] = ymothervertex[6];
  xmothervertex[10] = -xmothervertex[0], ymothervertex[10] = ymothervertex[9];
  xmothervertex[11] = xmothervertex[10], ymothervertex[11] = ymothervertex[0];
  //////////////////////////////////////////////////////////
  ssdchipmothershape->DefinePolygon(kssdmothervertexnumber,
									xmothervertex,ymothervertex);
  ssdchipmothershape->DefineSection(0,-0.5*fgkSSDChipHeight);
  ssdchipmothershape->DefineSection(1, 0.5*fgkSSDChipHeight);
  TGeoVolume* ssdchipmother = new TGeoVolume("SSDChipContainer",
							  ssdchipmothershape,fSSDAir);
   /////////////////////////////////////////////////////////////
  for(Int_t i=0; i<kssdchiprownumber; i++)
    for(Int_t j=0; j<fgkSSDChipNumber; j++) 
		ssdchipmother->AddNode(ssdchip,fgkSSDChipNumber*i+j+1,
		new TGeoTranslation(j*ssdchipseparation[0],i*ssdchipseparation[1],0.));
  return ssdchipmother;
}
/////////////////////////////////////////////////////////////////////////////////
TList* AliITSv11GeometrySSD::GetLadderCableSegment(Double_t ssdendladdercablelength){
  /////////////////////////////////////////////////////////////
  // Method returning a List containing pointers to Ladder Cable Volumes    
  /////////////////////////////////////////////////////////////
  const Int_t kladdercablesegmentnumber = 2;
  /////////////////////////////////////////
  // LadderSegmentBBox Volume
  /////////////////////////////////////////
  TGeoBBox* laddercablesegmentbboxshape[kladdercablesegmentnumber];
  const char* laddercablesegmentbboxshapename[kladdercablesegmentnumber] = 
				{"LadderCableSegmentBBoxShape1","LadderCableSegmentBBoxShape2"};
  for(Int_t i=0; i<kladdercablesegmentnumber; i++) laddercablesegmentbboxshape[i] = 
						  new TGeoBBox(laddercablesegmentbboxshapename[i],
									   0.5*fgkSSDFlexWidth[0],
									   0.5*fgkSSDLadderCableWidth,
									   0.5*fgkSSDFlexHeight[i]); 
  const char* laddercablesegmentbboxname[kladdercablesegmentnumber] = 
						  {"LadderCableSegmentBBox1","LadderCableSegmentBBox2"};
  TGeoVolume* laddercablesegmentbbox[kladdercablesegmentnumber];
  for(Int_t i=0; i<kladdercablesegmentnumber; i++){ 
			laddercablesegmentbbox[i] =
						  new TGeoVolume(laddercablesegmentbboxname[i],
										 laddercablesegmentbboxshape[i],
										 (i==0?fSSDAlTraceLadderCableMedium:
            fSSDKaptonLadderCableMedium));
			laddercablesegmentbbox[i]->SetLineColor(i==0 ? fColorAl : 
														   fColorPolyhamide);
  }
  TGeoTranslation* laddercablesegmentbboxtrans[kladdercablesegmentnumber];										  
  laddercablesegmentbboxtrans[0] = 
						   new TGeoTranslation("LadderCableSegmentBBoxTrans1",
											   0.5*fgkSSDFlexWidth[0],
											   0.5*fgkSSDLadderCableWidth,
											   0.5*fgkSSDFlexHeight[0]);
  laddercablesegmentbboxtrans[1] = 
						   new TGeoTranslation("LadderCableSegmentBBoxTrans2",
											   0.5*fgkSSDFlexWidth[0],
											   0.5*fgkSSDLadderCableWidth,
											   fgkSSDFlexHeight[0]
											   +0.5*fgkSSDFlexHeight[1]);
  TGeoVolume* laddercablesegmentbboxassembly = 
						   new TGeoVolumeAssembly("LadderCableSegmentBBoxAssembly"); 
  for(Int_t i=0; i<kladdercablesegmentnumber; i++)  
		laddercablesegmentbboxassembly->AddNode(laddercablesegmentbbox[i],1,
											    laddercablesegmentbboxtrans[i]);
/////////////////////////////////////////
// LadderSegmentArb8 Volume
/////////////////////////////////////////
  const Int_t kvertexnumber = 4;
  TVector3** laddercablesegmentvertexposition[kladdercablesegmentnumber];
  for(Int_t i = 0; i<kladdercablesegmentnumber; i++) laddercablesegmentvertexposition[i] = 
												  new TVector3*[kvertexnumber];
//Shape Vertex Positioning
  for(Int_t i=0; i<kladdercablesegmentnumber; i++){
	laddercablesegmentvertexposition[i][0] = new TVector3(0.,i*fgkSSDFlexHeight[0]);
	laddercablesegmentvertexposition[i][1] = new TVector3(fgkSSDLadderCableWidth,
														  i*fgkSSDFlexHeight[0]);
	laddercablesegmentvertexposition[i][2] = new TVector3(0.,fgkSSDFlexHeight[0]
										   +			     fgkSSDFlexHeight[1]
										   +			  i*fgkSSDFlexHeight[0]);
	laddercablesegmentvertexposition[i][3] = 
						   new TVector3(laddercablesegmentvertexposition[i][1]->X(),
										laddercablesegmentvertexposition[i][2]->Y());
  }
  Double_t laddercablesegmentwidth[2][2] = {{fgkSSDFlexHeight[0],fgkSSDFlexHeight[0]},
							     		    {fgkSSDFlexHeight[1],fgkSSDFlexHeight[1]}};	
  char* laddercablesegmentarbshapename[kladdercablesegmentnumber] = 
					{"LadderCableSegmentArbShape1","LadderCableSegmentArbShape2"};
  TGeoArb8* laddercablesegmentarbshape[kladdercablesegmentnumber];
  for(Int_t i = 0; i< kladdercablesegmentnumber; i++) laddercablesegmentarbshape[i] = 
					GetArbShape(laddercablesegmentvertexposition[i],
								laddercablesegmentwidth[i],
								fgkCarbonFiberJunctionWidth-fgkSSDFlexWidth[0],
								laddercablesegmentarbshapename[i]);
  const char* laddercablesegmentarbname[kladdercablesegmentnumber] = 
						  {"LadderCableSegmentArb1","LadderCableSegmentArb2"};
  TGeoVolume* laddercablesegmentarb[kladdercablesegmentnumber];
  for(Int_t i=0; i<kladdercablesegmentnumber; i++){
			 laddercablesegmentarb[i] =
						   new TGeoVolume(laddercablesegmentarbname[i],
										  laddercablesegmentarbshape[i],
										  (i==0?fSSDAlTraceLadderCableMedium:
            fSSDKaptonLadderCableMedium)); 
			 laddercablesegmentarb[i]->SetLineColor(i==0 ? fColorAl : 
														   fColorPolyhamide);
}
  TGeoRotation* laddercablesegmentarbrot[kladdercablesegmentnumber];
  laddercablesegmentarbrot[0] = new TGeoRotation("LadderCableSegmentArbRot1",
												 90.,90,-90.);	 
  laddercablesegmentarbrot[1] = new TGeoRotation("LadderCableSegmentArbRot2",
												  0.,90.,0.);	 
  TGeoCombiTrans* laddercablesegmentarbcombitrans =  
						   new TGeoCombiTrans("LadderCableSegmentArbCombiTrans",
							   0.5*(fgkCarbonFiberJunctionWidth-fgkSSDFlexWidth[0])
							 + fgkSSDFlexWidth[0],0.,0.,
						   new TGeoRotation((*laddercablesegmentarbrot[1])
						     *(*laddercablesegmentarbrot[0])));
  TGeoVolume* laddercablesegmentarbassembly = 
						   new TGeoVolumeAssembly("LadderCableSegmentArbAssembly"); 
  for(Int_t i=0; i<kladdercablesegmentnumber; i++)
  laddercablesegmentarbassembly->AddNode(laddercablesegmentarb[i],1,
										   laddercablesegmentarbcombitrans);
/////////////////////////////////////////
// End Ladder Cable Volume
/////////////////////////////////////////
  TGeoBBox* ladderendcablesegmentbboxshape[kladdercablesegmentnumber];
  const char* ladderendcablesegmentbboxshapename[kladdercablesegmentnumber] = 
				{"LadderEndCableSegmentBBoxShape1","LadderEndCableSegmentBBoxShape2"};
  for(Int_t i=0; i<kladdercablesegmentnumber; i++) ladderendcablesegmentbboxshape[i] = 
						  new TGeoBBox(ladderendcablesegmentbboxshapename[i],
									   0.5*ssdendladdercablelength,
									   0.5*fgkSSDLadderCableWidth,
									   0.5*fgkSSDFlexHeight[i]);
  const char* ladderendcablesegmentbboxname[kladdercablesegmentnumber] = 
						  {"LadderEndCableSegmentBBox1","LadderEndCableSegmentBBox2"};
  TGeoVolume* ladderendcablesegmentbbox[kladdercablesegmentnumber];
  for(Int_t i=0; i<kladdercablesegmentnumber; i++){ 
			ladderendcablesegmentbbox[i] =
						  new TGeoVolume(ladderendcablesegmentbboxname[i],
										 ladderendcablesegmentbboxshape[i],
										 (i==0?fSSDAlTraceLadderCableMedium:
            fSSDKaptonLadderCableMedium));
			ladderendcablesegmentbbox[i]->SetLineColor(i==0 ? fColorAl : 
														   fColorPolyhamide);
  }
  TGeoTranslation* ladderendcablesegmentbboxtrans[kladdercablesegmentnumber];										  
  ladderendcablesegmentbboxtrans[0] = 
						   new TGeoTranslation("LadderEndCableSegmentBBoxTrans0",
											   0.5*ssdendladdercablelength,
											   0.5*fgkSSDLadderCableWidth,
											   0.5*fgkSSDFlexHeight[0]);
  ladderendcablesegmentbboxtrans[1] = 
						   new TGeoTranslation("LadderEndCableSegmentBBoxTrans1",
											   0.5*ssdendladdercablelength,
											   0.5*fgkSSDLadderCableWidth,
											   fgkSSDFlexHeight[0]
											   +0.5*fgkSSDFlexHeight[1]);
  TGeoVolume* ladderendcablesegmentbboxassembly = 
						   new TGeoVolumeAssembly("LadderEndCableSegmentBBoxAssembly"); 
  for(Int_t i=0; i<kladdercablesegmentnumber; i++)  
		ladderendcablesegmentbboxassembly->AddNode(ladderendcablesegmentbbox[i],1,
											    ladderendcablesegmentbboxtrans[i]);
/////////////////////////////////////////
  TList* laddercablesegmentlist = new TList();
  laddercablesegmentlist->Add(laddercablesegmentbboxassembly);
  laddercablesegmentlist->Add(laddercablesegmentarbassembly);
  laddercablesegmentlist->Add(ladderendcablesegmentbboxassembly);
  return laddercablesegmentlist;
  }
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetLadderCable(Int_t n, Double_t ssdendladdercablelength){
  /////////////////////////////////////////////////////////////
  // Method generating Ladder Cable Volumes Assemblies    
  /////////////////////////////////////////////////////////////
  TList* laddercablesegmentlist = GetLadderCableSegment(ssdendladdercablelength);
  TGeoVolume* laddercable = new TGeoVolumeAssembly("LadderCable"); 
  for(Int_t i=0; i<n; i++){
	 TGeoTranslation* laddercabletrans = new TGeoTranslation(
							i*(fgkCarbonFiberJunctionWidth),
							fgkSSDLadderCableWidth-fgkSSDFlexWidth[0],
							i*(fgkSSDFlexHeight[0]+fgkSSDFlexHeight[1]));
    laddercable->AddNode((TGeoVolume*)laddercablesegmentlist->At(0),i+1,laddercabletrans);  
	if(i!=n-1) laddercable->AddNode((TGeoVolume*)laddercablesegmentlist->At(1),i+1,laddercabletrans);  
  }
  TGeoTranslation* endladdercabletrans = new TGeoTranslation("EndLadderCableTrans",
					  (n-1)*fgkCarbonFiberJunctionWidth+fgkSSDFlexWidth[0],
								 fgkSSDLadderCableWidth-fgkSSDFlexWidth[0],
					  (n-1)*(fgkSSDFlexHeight[0]+fgkSSDFlexHeight[1]));
  laddercable->AddNode((TGeoVolume*)laddercablesegmentlist->At(2),1,endladdercabletrans);
  return laddercable;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetLadderCableAssembly(Int_t n, Double_t ssdendladdercablelength){
  /////////////////////////////////////////////////////////////
  // Method generating Ladder Cable Volumes Assembly   
  /////////////////////////////////////////////////////////////
  TGeoVolume* laddercableassembly = new TGeoVolumeAssembly("LadderCableAssembly");
  char laddercabletransname[30];
  for(Int_t i=0; i<n; i++){ 
	sprintf(laddercabletransname,"LadderCableTrans%i",i+1);
    laddercableassembly->AddNode(GetLadderCable(n-i,ssdendladdercablelength),i+1,
	new TGeoTranslation(laddercabletransname,i*fgkCarbonFiberJunctionWidth,0.,0.));
  }
  return laddercableassembly;
}
/////////////////////////////////////////////////////////////////////////////////
TList* AliITSv11GeometrySSD::GetLadderCableAssemblyList(Int_t n, Double_t ssdendladdercablelength){
  /////////////////////////////////////////////////////////////
  // Method generating Ladder Cable List Assemblies  
  /////////////////////////////////////////////////////////////  
  const Int_t kladdercableassemblynumber = 2;
  TGeoVolume* laddercableassembly = GetLadderCableAssembly(n,ssdendladdercablelength);
  TGeoVolume* ladderCable[kladdercableassemblynumber];
  char laddercableassemblyname[30];
  TList* laddercableassemblylist = new TList();
  for(Int_t i=0; i<kladdercableassemblynumber; i++){ 
	sprintf(laddercableassemblyname,"LadderCableAssembly%i",i+1);
	ladderCable[i] = new TGeoVolumeAssembly(laddercableassemblyname);
	ladderCable[i]->AddNode(laddercableassembly,i+1,i==0 ? NULL :
					 new TGeoCombiTrans((n-1)
					 *	 fgkCarbonFiberJunctionWidth+fgkSSDFlexWidth[0],
					     2.*fgkSSDLadderCableWidth+0.5*fgkSSDFlexWidth[0],
											0.,new TGeoRotation("",180,0.,0.)));
	laddercableassemblylist->Add(ladderCable[i]);
}
  return laddercableassemblylist;
}
///////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetLadderSegment(){
  /////////////////////////////////////////////////////////////
  // Method Generating Ladder Segment Array
  /////////////////////////////////////////////////////////////
  fladdersegment[0] = new TGeoVolumeAssembly("LadderSegment1");	
  fladdersegment[1] = new TGeoVolumeAssembly("LadderSegment2");	
  if(!fCreateMaterials) CreateMaterials();
  if(!fTransformationMatrices) CreateTransformationMatrices();
  if(!fBasicObjects) CreateBasicObjects();
  for(Int_t i=0; i<fgkladdersegmentnumber; i++){
  // Placing Carbon Fiber Support	
	for(Int_t j=0; j<fgkcarbonfibersupportnumber; j++){ 
		fladdersegment[i]->AddNode(fcarbonfibersupport[0],j+1,
											fcarbonfibersupportmatrix[j]);	
		fladdersegment[i]->AddNode(fcarbonfibersupport[1],j+1,
											fcarbonfibersupportmatrix[j]);
  }
  // Placing Carbon Fiber Junction
    for(Int_t j=0; j<fgkcarbonfiberjunctionumber; j++)
        fladdersegment[i]->AddNode(fcarbonfiberjunction,j+1,
								   fcarbonfiberjunctionmatrix[j]);
  // Placing Carbon Fiber Lower Support
	for(Int_t j=0; j<fgkcarbonfiberlowersupportnumber; j++)
		fladdersegment[i]->AddNode(fcarbonfiberlowersupport[j],j+1,
			           			   fcarbonfiberlowersupportrans[j]);	
  // Placing SSD Sensor Support
    for(Int_t j=0; j<fgkssdsensorsupportnumber; j++) 
	fladdersegment[i]->AddNode(j<2 ? fssdsensorsupport[0][i] :
								     fssdsensorsupport[1][i],
							   j+1,fssdsensorsupportmatrix[j]);
  // Placing SSD Cooling Tube Support 
	for(Int_t j=0; j<fgkcoolingtubesupportnumber; j++)
		fladdersegment[i]->AddNode(fcoolingtubesupport,j+1,
								   fcoolingtubesupportmatrix[j]);
  // Placing SSD Cooling Tube  
	for(Int_t j=0; j<2; j++)
		for(Int_t k=0; k<2; k++){
		fladdersegment[i]->AddNode(fcoolingtube[0],2*j+k+1,fcoolingtubematrix[j][k]);
		fladdersegment[i]->AddNode(fcoolingtube[j+1],k+1,fcoolingtubematrix[2+j][k]);
		}
  // Placing SSD Hybrid
    switch(i){
	case 0: 
		fladdersegment[i]->AddNode(fssdhybridcomponent[0],1,fhybridmatrix);
		fladdersegment[i]->AddNode(fssdhybridcomponent[2],1,fhybridmatrix);
		break;
    case 1:
		fladdersegment[i]->AddNode(fssdhybridcomponent[1],1,fhybridmatrix);
		fladdersegment[i]->AddNode(fssdhybridcomponent[2],1,fhybridmatrix);
		break;
	}
	// Placing Cooling Block System
    fladdersegment[i]->AddNode(fssdcoolingblocksystem,1,fcoolingblocksystematrix);
	// Placing SSD Flex
	for(Int_t j=0; j<fgkflexnumber; j++){
      fladdersegment[i]->AddNode(fssdstiffenerflex,j+1,fstiffenerflexmatrix[j]);
      fladdersegment[i]->AddNode(fssdendflex,j+1,fendflexmatrix[j]);
	}
   }
}
///////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetEndLadderSegment(){
  /////////////////////////////////////////////////////////////
  // Method Generating End Ladder
  /////////////////////////////////////////////////////////////
  // End Ladder Carbon Fiber Junction 
  /////////////////////////////////////////////////////////////
  fendladdersegment[0] = new TGeoVolumeAssembly("EndLadder1");
  fendladdersegment[1] = new TGeoVolumeAssembly("EndLadder2");
  if(!fCreateMaterials) CreateMaterials();
  if(!fTransformationMatrices) CreateTransformationMatrices();
  if(!fBasicObjects) CreateBasicObjects();
  for(Int_t i=0; i<fgkendlabbercarbonfiberjunctionumber; i++){
	for(Int_t j=0; j<fgkendladdercabonfiberjunctionmatrixnumber; j++)
		fendladdersegment[i]->AddNode(j==2 ? 
							fendladdercarbonfiberjunction[i][1] : 
							fendladdercarbonfiberjunction[i][0],
							j+1,fendladdercarbonfiberjunctionmatrix[i][j]);
  }
  /////////////////////////////////////////////////////////////
  // End Ladder Carbon Fiber Support 
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i<fgkcarbonfibersupportnumber; i++)
      for(Int_t j=0; j<fgkendladdercarbonfibermatrixnumber; j++){
		  fendladdersegment[i]->AddNode(fcarbonfibersupport[0],j+1,
                  fendladdercarbonfibermatrix[i][j]);	
          fendladdersegment[i]->AddNode(fcarbonfibersupport[1],j+1,
                  fendladdercarbonfibermatrix[i][j]);	
      }
  /////////////////////////////////////////////////////////////
  // End Ladder Mounting Block
  /////////////////////////////////////////////////////////////
 // for(Int_t i=0; i<fgkendladdermountingblocknumber; i++) 
 //      fendladdersegment[i]->AddNode(fendladdermountingblock,1,
 //									 fendladdermountingblocktrans[i]);
  /////////////////////////////////////////////////////////////
  // End Ladder Lower Supports
  /////////////////////////////////////////////////////////////
  fendladdersegment[0]->AddNode(fcarbonfiberlowersupport[0],1,
								fendladderlowersupptrans[0]);
  fendladdersegment[1]->AddNode(fcarbonfiberlowersupport[0],2,
								fendladderlowersupptrans[1]);
  fendladdersegment[1]->AddNode(fcarbonfiberlowersupport[0],3,
								fendladderlowersupptrans[2]);
  //fendladdersegment[0]->CheckOverlaps(0.01);
  //fendladdersegment[1]->CheckOverlaps(0.01);
}
///////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetLadder(){
  /////////////////////////////////////////////////////////////
  // Method Generating Ladder of Layer 5 and 6
  /////////////////////////////////////////////////////////////  
  Int_t ssdlaysensorsnumber[fgkladdernumber] = {fgkSSDLay5SensorsNumber,
												fgkSSDLay6SensorsNumber};
  /////////////////////////////////////////////////////////////////////////////						
  /// Generating Ladder Mother Volume Containing Ladder 
  /////////////////////////////////////////////////////////////////////////////		 
  TGeoXtru* laddershape[fgkladdernumber];	
  for(Int_t i=0; i<fgkladdernumber; i++) laddershape[i] = new TGeoXtru(2);
  const Int_t kmothervertexnumber = 8;  
  Double_t xmothervertex[fgkladdernumber][kmothervertexnumber];
  Double_t ymothervertex[fgkladdernumber][kmothervertexnumber];
  ///////////////////////
  // Setting the vertices 
  ///////////////////////
  Double_t laddercablethickness = (fgkSSDLay6SensorsNumber+2)
								* (fgkSSDFlexHeight[0]+fgkSSDFlexHeight[1]);
  xmothervertex[0][0] = -0.5*fgkSSDSensorWidth;
  ymothervertex[0][0] = -0.5*fgkCoolingTubeSupportHeight-fgkSSDModuleCoolingBlockToSensor;
  xmothervertex[0][1] = xmothervertex[0][0];
  ymothervertex[0][1] = 0.0;
  xmothervertex[0][2] = - 0.5*fgkSSDModuleSensorSupportDistance-2.*fgkCoolingTubeSupportRmax
						- laddercablethickness/SinD(2.*fgkSSDFlexAngle);
  ymothervertex[0][2] = ymothervertex[0][1];
  ymothervertex[0][3] = 0.5*fgkCarbonFiberTriangleLength*TanD(2.*fgkSSDFlexAngle);
  xmothervertex[0][3] = xmothervertex[0][2]+ymothervertex[0][3]/TanD(2.*fgkSSDFlexAngle);		
  xmothervertex[0][4] = -xmothervertex[0][3];
  ymothervertex[0][4] = ymothervertex[0][3];
  xmothervertex[0][5] = -xmothervertex[0][2];
  ymothervertex[0][5] = ymothervertex[0][2];
  xmothervertex[0][6] = -xmothervertex[0][1];
  ymothervertex[0][6] = ymothervertex[0][1];
  xmothervertex[0][7] = -xmothervertex[0][0];
  ymothervertex[0][7] = ymothervertex[0][0];
  for(Int_t i=0; i<kmothervertexnumber; i++){
	xmothervertex[1][i] = xmothervertex[0][i];
	ymothervertex[1][i] = ymothervertex[0][i];
  }
  const char* laddername[fgkladdernumber] = {"ITSssdLay5Ladd","ITSssdLay6Ladd"};
  for(Int_t i=0; i<fgkladdernumber; i++){
	laddershape[i]->DefinePolygon(kmothervertexnumber,xmothervertex[i],
								    ymothervertex[i]);
    laddershape[i]->DefineSection(0,-fgkEndLadderCarbonFiberLowerJunctionLength[1]);
    laddershape[i]->DefineSection(1,ssdlaysensorsnumber[i]*fgkCarbonFiberJunctionWidth
											+fgkEndLadderCarbonFiberLowerJunctionLength[0]);
    fladder[i] = new TGeoVolume(laddername[i],laddershape[i],fSSDAir);
 }
///////////////////////////////////////////////////////////////////////////
 if(!fCreateMaterials) CreateMaterials();
 if(!fTransformationMatrices) CreateTransformationMatrices();
 if(!fBasicObjects) CreateBasicObjects();
 SetLadderSegment(); 
 SetEndLadderSegment();
  for(Int_t i=0; i<fgkladdernumber; i++){
	for(Int_t j=0; j<ssdlaysensorsnumber[i]; j++){
	//////////////////////////						
	/// Placing Ladder Segment
	//////////////////////////		
		fladder[i]->AddNode(j%2==0 ? fladdersegment[i==0 ? 0 : 1] :
								     fladdersegment[i==0 ? 1 : 0],
									 ssdlaysensorsnumber[i]-j-1,fladdermatrix[i][j]);
	//////////////////////////						
	/// Placing SSD Sensor
	//////////////////////////		
		fladder[i]->AddNode(i==0?fSSDSensor5:fSSDSensor6,ssdlaysensorsnumber[i]-j-1,
							fssdsensormatrix[i][j]);
	}
	///////////////////////////////						
	/// Placing End Ladder Segment
	///////////////////////////////		
    fladder[i]->AddNode(fendladdersegment[0],1,fendladdersegmentmatrix[0][i]);
	fladder[i]->AddNode(fendladdersegment[1],1,fendladdersegmentmatrix[1][i]);
   }
/////////////////////////////////////////////////////////////////////////////						
/// Placing Ladder Cables
/////////////////////////////////////////////////////////////////////////////		
  Int_t sidecablenumber[2][2];
  sidecablenumber[0][0] = fgkSSDLay5SensorsNumber/2+1; 
  sidecablenumber[0][1] = sidecablenumber[0][0]-2;
  sidecablenumber[1][0] = (fgkSSDLay6SensorsNumber-1)/2+1;
  sidecablenumber[1][1] = sidecablenumber[1][0]-1;
  Double_t carbonfibertomoduleposition[3];
  carbonfibertomoduleposition[0] = -0.5*(fgkSSDSensorWidth-fgkCarbonFiberTriangleLength);
  carbonfibertomoduleposition[1] = - (2.*fgkSSDSensorLength-fgkSSDSensorOverlap)+
			 fgkSSDModuleStiffenerPosition[1]+fgkSSDStiffenerWidth
	 +		 0.5*fgkSSDFlexHoleLength+2.*fgkCarbonFiberJunctionWidth
	 -		 0.5*(fgkCarbonFiberLowerSupportWidth+fgkSSDSensorCenterSupportLength
	 -            fgkSSDSensorCenterSupportThickness[0]);
  carbonfibertomoduleposition[2] = - (fgkSSDModuleCoolingBlockToSensor
								 +   0.5*fgkCoolingTubeSupportHeight
	 -         fgkSSDSensorHeight-fgkSSDChipCablesHeight[3]-fgkSSDChipHeight);	
  const Double_t kendladdercablecorrection = 1.72*fgkmm; //this has to be tuned
  Double_t ssdendladdercablelength[4];
  ssdendladdercablelength[0] = carbonfibertomoduleposition[1]
							 + fgkSSDSensorLength
							 - fgkSSDModuleStiffenerPosition[1]
							 - fgkSSDStiffenerWidth 
							 - fgkSSDFlexWidth[0]
							 + fgkEndLadderCarbonFiberLowerJunctionLength[1]-0.000001*kendladdercablecorrection;
  ssdendladdercablelength[1] = carbonfibertomoduleposition[1]
							 + fgkSSDModuleStiffenerPosition[1]
							 + fgkSSDStiffenerWidth
							 + fgkEndLadderCarbonFiberLowerJunctionLength[1]-0.000001*kendladdercablecorrection;
  ssdendladdercablelength[2] = ssdendladdercablelength[1]
							 - fgkEndLadderCarbonFiberLowerJunctionLength[1]
							 + fgkEndLadderCarbonFiberLowerJunctionLength[0]
							 - kendladdercablecorrection;
  ssdendladdercablelength[3] = fgkCarbonFiberJunctionWidth-(fgkSSDSensorLength
							 + carbonfibertomoduleposition[1]
							 - fgkSSDModuleStiffenerPosition[1]
							 - fgkSSDStiffenerWidth)
							 + fgkEndLadderCarbonFiberLowerJunctionLength[0]-0.000001*kendladdercablecorrection;
  TList* laddercableassemblylist[4];
  const Int_t kendladdercablesnumber = 4;
  for(Int_t i=0; i<fgkladdercablesnumber; i++)
	for(Int_t j=0; j<kendladdercablesnumber; j++){
		laddercableassemblylist[j] = 
		GetLadderCableAssemblyList(sidecablenumber[i][j<2?0:1],
								   ssdendladdercablelength[j]);
	    fladder[i]->AddNode((TGeoVolume*)laddercableassemblylist[j]->At(j%2==0?0:1),
									j<2?1:2,fladdercablematrix[i][j]);
  }
  //fladder[0]->CheckOverlaps(0.01);
  //fladder[1]->CheckOverlaps(0.01);
}
////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetLayer(){
////////////////////////////////////////////////////////////////////////////////
  // Creating Ladder of Layer 5 and Layer 6
  /////////////////////////////////////////////////////////////
  if(!fCreateMaterials) CreateMaterials();
  if(!fTransformationMatrices) CreateTransformationMatrices();
  if(!fBasicObjects) CreateBasicObjects();
  SetLadder(); // Generating the ladder of Layer5 and Layer6
  const Int_t kssdlayladdernumber[fgklayernumber] = 
			{fgkSSDLay5LadderNumber,fgkSSDLay6LadderNumber};
  /////////////////////////////////////////////////////////////
  // Generating mother volumes for Layer5 and Layer6
  /////////////////////////////////////////////////////////////
  TGeoXtru* ssdladdermothershape = (TGeoXtru*)fladder[0]->GetShape();
  TGeoTube* ssdlayershape[fgklayernumber];
  ssdlayershape[0] = new TGeoTube(fgkSSDLay5RadiusMin,fgkSSDLay5RadiusMax
				   -              ssdladdermothershape->GetY(0)
				   +			  TMath::Sqrt(TMath::Power(ssdladdermothershape->GetY(4),2.)
				   +			  TMath::Power(ssdladdermothershape->GetX(4),2.)),
								  0.5*fgkSSDLay5LadderLength);
  ssdlayershape[1] = new TGeoTube(fgkSSDLay6RadiusMin,fgkSSDLay6RadiusMax
				   -			  ssdladdermothershape->GetY(0)
				   +			  TMath::Sqrt(TMath::Power(ssdladdermothershape->GetY(4),2.)
				   +			  TMath::Power(ssdladdermothershape->GetX(4),2.)),
								  0.5*fgkSSDLay6LadderLength);
  fSSDLayer5 = new TGeoVolume("ITSssdLayer5",ssdlayershape[0],fSSDAir);
  fSSDLayer6 = new TGeoVolume("ITSssdLayer6",ssdlayershape[1],fSSDAir);
  /////////////////////////////////////////////////////////////
  Int_t *ladderindex[fgklayernumber];
  Int_t index[fgklayernumber] = {8,9};
  for(Int_t i=0; i<fgklayernumber; i++) ladderindex[i] = new Int_t[kssdlayladdernumber[i]];
  for(Int_t i=0; i<fgklayernumber; i++)	
	for(Int_t j=0; j<kssdlayladdernumber[i]; j++){
		ladderindex[i][j] = ((j>=0)&&(j<=kssdlayladdernumber[i]-index[i]-1)) ? 
							  j+index[i] : j+index[i]-kssdlayladdernumber[i]; 
		i ==0 ? fSSDLayer5->AddNode(fladder[0],ladderindex[i][j],flayermatrix[i][j]) : 
		        fSSDLayer6->AddNode(fladder[1],ladderindex[i][j],flayermatrix[i][j]);
	}
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
	for(Int_t i=0; i<fgklayernumber; i++) delete ladderindex[i];
}
////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::Layer5(TGeoVolume* moth){
  /////////////////////////////////////////////////////////////
  // Insert the layer 5 in the mother volume. 
  /////////////////////////////////////////////////////////////
  if (! moth) {
    printf("Error::AliITSv11GeometrySSD: Can't insert layer5, mother is null!\n");
    return;
  };
  if(!fSSDLayer5) SetLayer();
  fMotherVol = moth;
  TGeoTranslation* centerITSlayer5trans = new TGeoTranslation(0.,0.,-0.5*fgkSSDLay5LadderLength
										+ fgkLay5CenterITSPosition);
  moth->AddNode(fSSDLayer5,1,centerITSlayer5trans);
 }
////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::Layer6(TGeoVolume* moth){
  /////////////////////////////////////////////////////////////
  // Insert the layer 6 in the mother volume. 
  /////////////////////////////////////////////////////////////
  if (! moth) {
    printf("Error::AliITSv11GeometrySSD: Can't insert layer6, mother is null!\n");
    return;
  };
  if(!fSSDLayer6) SetLayer();
  fMotherVol = moth;
  TGeoTranslation* centerITSlayer6trans = new TGeoTranslation(0.,0.,-0.5*fgkSSDLay6LadderLength
										+ fgkLay6CenterITSPosition);
  moth->AddNode(fSSDLayer6,1,centerITSlayer6trans);
 }
 ////////////////////////////////////////////////////////////////////////////////
TGeoArb8* AliITSv11GeometrySSD::GetArbShape(TVector3* vertexpos[], Double_t* width, 
									Double_t height, char* shapename, Int_t isign) const{
  /////////////////////////////////////////////////////////////
  // Method generating an Arb shape 
  /////////////////////////////////////////////////////////////
  const Int_t kvertexnumber = 8;
  const Int_t ktransvectnumber = 2;
  TVector3* vertex[kvertexnumber];
  TVector3* transvector[2];
  for(Int_t i=0; i<ktransvectnumber; i++) transvector[i] = new TVector3(0.,width[i]);
  /////////////////////////////////////////////////////////////
  //Setting the vertices for TGeoArb8
  /////////////////////////////////////////////////////////////
  vertex[0] = new TVector3(*vertexpos[0]);
  vertex[1] = new TVector3(*vertexpos[1]);
  vertex[2] = new TVector3(*vertex[1]+isign*(*transvector[0]));
  vertex[3] = new TVector3(*vertex[0]+isign*(*transvector[0]));
  vertex[4] = new TVector3(*vertexpos[2]);
  vertex[5] = new TVector3(*vertexpos[3]);
  vertex[6] = new TVector3(*vertex[5]+isign*(*transvector[1]));
  vertex[7] = new TVector3(*vertex[4]+isign*(*transvector[1]));
  /////////////////////////////////////////////////////////////
  TGeoArb8* arbshape = new TGeoArb8(shapename,0.5*height);
  for(Int_t i = 0; i<kvertexnumber;i++) 
							arbshape->SetVertex(i,vertex[i]->X(),vertex[i]->Y());
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i< kvertexnumber; i++) delete vertex[i];  
  for(Int_t i=0; i< ktransvectnumber; i++) delete transvector[i];  
  /////////////////////////////////////////////////////////////
  return arbshape;
} 
///////////////////////////////////////////////////////////////////////////////
TGeoXtru* AliITSv11GeometrySSD::GetArcShape(Double_t phi, Double_t rmin, 
								Double_t rmax, Int_t nedges, Double_t height){
  /////////////////////////////////////////////////////////////
  // Method generating Arc shape 
  /////////////////////////////////////////////////////////////
	const Int_t kvertexnumber = 2*nedges+2;
	TGeoXtru* arcshape = new TGeoXtru(2);	
	TVector3** vertexposition[2];
	for(Int_t i=0; i<2; i++) vertexposition[i] = new TVector3*[nedges+1];
	Double_t angle = 0.;
    for(Int_t i=0; i<nedges+1; i++){ 
		angle = 90.+0.5*phi-i*(phi/nedges);
		vertexposition[0][i] = new TVector3(rmin*CosD(angle),rmin*SinD(angle));
		vertexposition[1][i] = new TVector3(rmax*CosD(angle),rmax*SinD(angle));
	}
	Double_t *xvertexpoints = new Double_t[kvertexnumber];
	Double_t *yvertexpoints = new Double_t[kvertexnumber];
	for(Int_t i=0; i<kvertexnumber; i++){ 
		if(i==0){ xvertexpoints[i] = vertexposition[0][i]->X(),
				  yvertexpoints[i] = vertexposition[0][i]->Y();	
		}
		else if(i>=1&&i<nedges+2)
		{
			xvertexpoints[i] = vertexposition[1][i-1]->X(); 
			yvertexpoints[i] = vertexposition[1][i-1]->Y(); 
		}
        else
		{
			xvertexpoints[i] = vertexposition[0][kvertexnumber-i]->X(); 
			yvertexpoints[i] = vertexposition[0][kvertexnumber-i]->Y(); 
		}
    }
  arcshape->DefinePolygon(kvertexnumber,xvertexpoints,yvertexpoints);
  arcshape->DefineSection(0,-0.5*height);
  arcshape->DefineSection(1,0.5*height);
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  for(Int_t i=0; i<2; i++){
	for(Int_t j=0; j<nedges+1; j++)
		delete vertexposition[i][j];
	delete [] vertexposition[i];
  }
  delete [] xvertexpoints;
  delete [] yvertexpoints;
  /////////////////////////////////////////////////////////////
	return arcshape;
}
////////////////////////////////////////////////////////////////////////////////
TVector3* AliITSv11GeometrySSD::GetReflection(TVector3* vector,Double_t* param) const{
  /////////////////////////////////////////////////////////////
  // Given an axis specified by param, it gives the reflection of the point
  // respect to the axis
  /////////////////////////////////////////////////////////////
  TVector3* n = new TVector3(param[0],param[1],param[2]);
  Double_t d = ((*vector)*(*n)+param[3])/n->Mag2();
  TVector3* reflectedvector = new TVector3(*vector-2*d*(*n));
  /////////////////////////////////////////////////////////////
  // Deallocating memory
  /////////////////////////////////////////////////////////////
  delete n;
  /////////////////////////////////////////////////////////////
  return reflectedvector;
}
////////////////////////////////////////////////////////////////////////////////
TGeoHMatrix* AliITSv11GeometrySSD::AddTranslationToHMatrix(TGeoHMatrix* ct,
                                                       Double_t dx,
                                                       Double_t dy,
                                                       Double_t dz) const{
  /////////////////////////////////////////////////////////////
  // Add a dx,dy,dz translation to the initial TGeoCombiTrans
  /////////////////////////////////////////////////////////////
  TGeoHMatrix* hmatrix = new TGeoHMatrix(*ct);
  const Double_t *vect = hmatrix->GetTranslation();
  Double_t newvect[3] = {vect[0]+dx, vect[1]+dy, vect[2]+dz};
  hmatrix->SetTranslation(newvect);
  TGeoHMatrix* matrix = new TGeoHMatrix(*hmatrix);
  delete hmatrix;
  return matrix;
}
////////////////////////////////////////////////////////////////////////////////
TGeoMedium* AliITSv11GeometrySSD::GetMedium(const char* mediumName) {
  /////////////////////////////////////////////////////////////
  // Method returning the Medium type 
  /////////////////////////////////////////////////////////////
  char ch[30];
  sprintf(ch, "ITS_%s",mediumName);
  TGeoMedium* medium =  gGeoManager->GetMedium(ch);
  if (! medium)
    printf("Error(AliITSv11GeometrySSD)::medium %s not found !\n", mediumName);
  return medium;
}
////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::CreateMaterials(){
///////////////////////////////////
// This part has to be modified
///////////////////////////////////
  ///////////////////////////////////
  // Silicon for Sensor
  /////////////////////////////////// 
  fSSDSensorMedium = GetMedium("SI$");
  ///////////////////////////////////
  // Silicon Mixture for Sensor
  /////////////////////////////////// 
  fSSDChipMedium = GetMedium("SPD SI CHIP$");
  fSSDChipGlueMedium = GetMedium("EPOXY$");
  ///////////////////////////////////
  // Stiffener Components Materials
  /////////////////////////////////// 
  fSSDStiffenerMedium = GetMedium("SDD C AL (M55J)$");
  ///////////////////////////  
  // Stiffener Connectors 
  ///////////////////////////  
  fSSDStiffenerConnectorMedium = GetMedium("COPPER$");
  ////////////////////////////////  
  // Stiffener 0603-1812 Capacitor
  ////////////////////////////////  
  fSSDStiffener0603CapacitorMedium = GetMedium("SDD ruby sph. Al2O3$");
  fSSDStiffener1812CapacitorMedium = GetMedium("SDD ruby sph. Al2O3$");
  ///////////////////////////  
  // Stiffener Hybrid Wire 
  ///////////////////////////  
  fSSDStiffenerHybridWireMedium = GetMedium("COPPER$");
  ///////////////////////////  
  // Al for Cooling Block
  ///////////////////////////  
  fSSDAlCoolBlockMedium = GetMedium("AL$");
  //////////////////////////////////////////////////////  
  // Kapton and Al for Chip Cable Flex and Ladder Cables
  //////////////////////////////////////////////////////  
  fSSDKaptonChipCableMedium = GetMedium("KAPTONH(POLYCH2)$");
  fSSDAlTraceChipCableMedium = GetMedium("AL$");
  fSSDKaptonFlexMedium = GetMedium("KAPTONH(POLYCH2)$");
  fSSDAlTraceFlexMedium = GetMedium("AL$");
  fSSDKaptonLadderCableMedium = GetMedium("KAPTONH(POLYCH2)$");
  fSSDAlTraceLadderCableMedium = GetMedium("AL$");
  /////////////////////////////////////////////////////////////////  
  // M55J for Carbon Fiber, CarbonFiber Lower Support and Junction
  //////////////////////////////////////////////////////////////////  
  fSSDCarbonFiberMedium = GetMedium("GEN C (M55J)$");
  /////////////////////////////////////////////////////////////////  
  // G10 for Detector Leg, TubeHolder
  //////////////////////////////////////////////////////////////////  
  fSSDTubeHolderMedium = GetMedium("G10FR4$");
  fSSDSensorSupportMedium = GetMedium("G10FR4$");
  fSSDMountingBlockMedium = GetMedium("G10FR4$");
  fSSDMountingBlockMedium = GetMedium("G10FR4$");
  /////////////////////////////////////////////////////////////////  
  // Water and Phynox for Cooling Tube
  //////////////////////////////////////////////////////////////////  
  fSSDCoolingTubeWater = GetMedium("WATER$");
  fSSDCoolingTubePhynox = GetMedium("INOX$");
  /////////////////////////////////////////////////////////////////////
  fSSDAir = GetMedium("SDD AIR$");
  fCreateMaterials = kTRUE;
}
/////////////////////////////////////////////////////////////////////
