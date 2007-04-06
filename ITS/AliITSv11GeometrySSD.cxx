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
#include "AliITSv11GeometrySSD.h"
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
const Double_t AliITSv11GeometrySSD::fgkSSDChipGlueLength        =   fgkSSDChipLength;
const Double_t AliITSv11GeometrySSD::fgkSSDChipGlueWidth         =   fgkSSDChipWidth;
const Double_t AliITSv11GeometrySSD::fgkSSDChipGlueHeight        =   0.030*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// Stiffener (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerLength       =  73.000*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerWidth        =   6.500*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerHeight       =   0.315*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDStiffenerToChipDist   =   2.500*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor0603Length   =   1.600*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor0603Width    =   0.870*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor0603Height   =   0.800*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor1812Length   =   4.600*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor1812Width    =   3.400*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDCapacitor1812Height   =   1.400*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDWireLength            =  30.000*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDWireRadius            =   0.185*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorPosition[2]  = {44.32*fgkmm, 
                                                                     0.33*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorSeparation   = 0.44*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorLength       = 2.16*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorWidth        = 3.60*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDConnectorHeight       = 
                                                   0.25*fgkSSDStiffenerHeight;
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
const Double_t AliITSv11GeometrySSD::fgkSSDSensorInsensitiveLength      = 1.*fgkmm; 
const Double_t AliITSv11GeometrySSD::fgkSSDSensorInsensitiveWidth       = 1.*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// Flex (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDFlexFullLength       =  106.000*fgkmm; 
const Double_t AliITSv11GeometrySSD::fgkSSDFlexLength[4]        = 
			{0.5*(fgkSSDStiffenerLength+fgkSSDChipNumber*fgkSSDChipLength
						+(fgkSSDChipNumber-1)*fgkSSDChipSeparationLength),
			 0.5*(fgkSSDStiffenerLength+fgkSSDChipNumber*fgkSSDChipLength
						+(fgkSSDChipNumber-1)*fgkSSDChipSeparationLength)-4.000*fgkmm,
						  9.500*fgkmm, 10.000*fgkmm};
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
const Double_t AliITSv11GeometrySSD::fgkSSDLadderCableWidth     =   23.5*fgkmm;
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
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportLength		   = 
																	   5.800*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportWidth          =  
																	   2.000*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportHeight[2]      =
														     { 4.620*fgkmm, 5.180*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportThickness[2] = 
														     { 0.450*fgkmm, 0.450*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDSensorSideSupportPosition       =  
									  0.5*(fgkSSDModuleSensorSupportDistance
							       +    fgkSSDSensorSideSupportThickness[0])
								   -           fgkSSDSensorSideSupportLength;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportLength	   =  
									   								   5.250*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportWidth        =
																       1.680*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportHeight[2]    =
	   {fgkSSDSensorSideSupportHeight[0]+fgkSSDSensorSideSupportThickness[0],
	   fgkSSDSensorSideSupportHeight[1]+fgkSSDSensorSideSupportThickness[1]};
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportThickness[2] =
   {fgkSSDSensorSideSupportThickness[0],fgkSSDSensorSideSupportThickness[1]};
const Double_t AliITSv11GeometrySSD::fgkSSDSensorCenterSupportPosition     = 
																      19.000*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// Chip Cables (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkSSDChipCablesLength[2]   = 
						 {73.12/fgkSSDChipNumber*fgkmm,fgkSSDChipLength+2.*0.19*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkSSDChipCablesHeight[4]   = 
						{  0.014*fgkmm,  0.010*fgkmm, fgkSSDModuleCoolingBlockToSensor
  -  (fgkSSDSensorSideSupportHeight[1]-fgkSSDSensorSideSupportHeight[0])
  -   fgkSSDCoolingBlockHoleCenter-fgkSSDStiffenerHeight
  -   fgkSSDChipHeight-fgkSSDSensorHeight,
      fgkSSDModuleCoolingBlockToSensor
  -   fgkSSDCoolingBlockHoleCenter-fgkSSDStiffenerHeight
  -   fgkSSDChipHeight-fgkSSDSensorHeight};
const Double_t AliITSv11GeometrySSD::fgkSSDChipCablesWidth[3]    = 
												 { 11.000*fgkmm,  0.800*fgkmm,  0.600*fgkmm};
/////////////////////////////////////////////////////////////////////////////////
// Carbon Fiber Junction Parameters (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberJunctionLength          = 
																	   3.820*fgkmm;
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
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberTriangleLength          = 
  fgkSSDModuleSensorSupportDistance-2.*fgkCarbonFiberJunctionToSensorSupport;  
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberTriangleAngle           = 
																	   60.00;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberSupportTopEdgeDist[2]   = 
														   {  0.751*fgkmm,  0.482*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberSupportEdgeLength       =  
																	   1.630*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberSupportWidth            = 
																	   0.950*fgkmm;
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
																	  =  0.950*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberLowerSupportLowerLenght       
																	  =  1.600*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCarbonFiberLowerSupportHeight            
																	  =  0.830*fgkmm;
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
		{fgkEndLadderCarbonFiberLowerJunctionLength[0]-20.4*fgkmm,
		 fgkEndLadderCarbonFiberLowerJunctionLength[1]-20.6*fgkmm};
const Double_t AliITSv11GeometrySSD::fgkEndLadderMountingBlockPosition[2] = 
						   {fgkEndLadderCarbonFiberLowerJunctionLength[0]-16.50*fgkmm,
						   fgkEndLadderCarbonFiberLowerJunctionLength[1]-31.50*fgkmm};
/////////////////////////////////////////////////////////////////////////////////
// Cooling Tube Support (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportRmax          =  1.45*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportRmin          
											  = fgkSSDCoolingBlockHoleRadius[0];
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportLength        =  8.55*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportHeight        =  0.85*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportWidth         =  2.00*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportSeparation    = 
                                      fgkSSDSensorLength-2.*fgkSSDSensorOverlap;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSupportToCarbonFiber = 11.70*fgkmm;
/////////////////////////////////////////////////////////////////////////////////
// Cooling Tube (lengths are in mm and angles in degrees)
/////////////////////////////////////////////////////////////////////////////////
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeRmax       = 
													  fgkCoolingTubeSupportRmin;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeRmin       =  0.96*fgkmm;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeLength     = 
													fgkCarbonFiberJunctionWidth;
const Double_t AliITSv11GeometrySSD::fgkCoolingTubeSeparation = 
					fgkSSDModuleSensorSupportDistance+fgkSSDCoolingBlockLength;
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
  fMotherVol(0),
  fColorCarbonFiber(4),
  fColorRyton(5),
  fColorPhynox(5),
  fColorSilicon(3),
  fColorAl(7),
  fColorKapton(6),
  fColorPolyhamide(5),
  fColorStiffener(9),
  fColorEpoxy(30)
{
  ////////////////////////
  // Standard constructor
  ////////////////////////
  CreateMaterials();
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
  fMotherVol(s.fMotherVol),
  fColorCarbonFiber(s.fColorCarbonFiber),
  fColorRyton(s.fColorRyton),
  fColorPhynox(s.fColorPhynox),
  fColorSilicon(s.fColorSilicon),
  fColorAl(s.fColorAl),
  fColorKapton(s.fColorKapton),
  fColorPolyhamide(s.fColorPolyhamide),
  fColorStiffener(s.fColorStiffener),
  fColorEpoxy(s.fColorEpoxy)
{
  ////////////////////////
  // Copy Constructor
  ////////////////////////
  CreateMaterials();
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
/////////////////////////////////////////////////////////////////////////////////
// Setting the transformation Matrices
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetSSDSensorSupportCombiTransMatrix(){
  /////////////////////////////////////////////////////////////
  // Method generating CombiTrans Matrix for SSD Sensor Support   
  /////////////////////////////////////////////////////////////
  //Translation Parameters SSDSensorSupportAssembly:
  ////////////////////////////////////////////////////////
  const Double_t kssdsensorsupporttransx[3] = {-0.5*fgkSSDSensorSideSupportWidth,
					                           0.5*fgkSSDSensorSideSupportWidth,
					                   0.5*fgkSSDSensorCenterSupportThickness[0]
									  -    fgkSSDSensorCenterSupportPosition}; 
  const Double_t kssdsensorsupporttransy[3] = 
									   {0.5*fgkSSDSensorSideSupportThickness[0],
					                   -0.5*fgkSSDSensorSideSupportThickness[0]
						               -fgkSSDModuleSensorSupportDistance,
					                    0.5*fgkSSDSensorCenterSupportWidth
									   -0.5*fgkSSDModuleSensorSupportDistance}; 
  const Double_t kssdsensorsupporttransz[3] = {0.,0.,
										fgkSSDSensorCenterSupportThickness[0]}; 
  ////////////////////////////////////////////////////////
  //Rotational Parameters SSDSensorSupportAssembly:
  ////////////////////////////////////////////////////////  
  const Double_t kssdsensorsupportrotphi[3]   = {   0., 180., 270.};
  const Double_t kssdsensorsupportrottheta[3] = {  90.,  90.,  90.};
  const Double_t kssdsensorsupportrotpsi[3]   = {- 90.,- 90.,- 90.};
  ////////////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of SSDSensorSupportAssembly:
  ////////////////////////////////////////////////////////////////
  char ssdsensorsupportcombitransname[40];
  char ssdsensorsupportrotname[40];
  TGeoCombiTrans *ssdsensorsupportlocalmatrix[fgkSSDSensorSupportCombiTransNumber];
  for(Int_t i=0; i<fgkSSDSensorSupportCombiTransNumber; i++){ 
		sprintf(ssdsensorsupportcombitransname,"SSDSensorSupportCombiTrans%i",i);
		sprintf(ssdsensorsupportrotname,"SSDSensorSupportRot%i",i);
		ssdsensorsupportlocalmatrix[i] =
				 new TGeoCombiTrans(ssdsensorsupportcombitransname,
									kssdsensorsupporttransx[i],
									kssdsensorsupporttransy[i],
									kssdsensorsupporttransz[i],
									new TGeoRotation(ssdsensorsupportrotname,
													   kssdsensorsupportrotphi[i],
													   kssdsensorsupportrottheta[i],
												    kssdsensorsupportrotpsi[i]));
           fSSDSensorSupportCombiTransMatrix[i] = ssdsensorsupportlocalmatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetSSDModuleCombiTransMatrix(Double_t SSDChipCablesHeigth){
  /////////////////////////////////////////////////////////////
  // Method generating CombiTrans Matrix for SSD Module   
  /////////////////////////////////////////////////////////////
  //Translation Parameters SSDModuleAssembly:
  ////////////////////////////////////////////////////////
  const Double_t kssdmoduletransx[7] = {0.5*fgkSSDStiffenerLength,
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
  const Double_t kssdmoduletransy[7] = {0.5*fgkSSDStiffenerWidth,
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
  const Double_t kssdmoduletransz[7] = {0.5*fgkSSDStiffenerHeight,
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
  const Double_t kssdmodulerotphi[7]   = {   0.,   0.,  90.,   0.,   0.,  90., 180.};
  const Double_t kssdmodulerottheta[7] = {   0.,   0.,   0.,   0.,   0.,   0.,   0.};
  const Double_t kssdmodulerotpsi[7]   = {   0.,   0.,   0.,   0.,   0.,   0.,   0.};
  ////////////////////////////////////////////////////////  
  //Name of CombiTrans Transformation of SSDModuleAssembly:
  ////////////////////////////////////////////////////////  
  const char* ssdmodulecombitransname[7] = {"SSDStiffenerCombiTrans",
											     "SSDChipCombiTrans",
											   "SSDSensorCombiTrans",
											    "SSDFlex0CombiTrans",
										 "SSDCoolingBlockCombiTrans",
										   "SSDChipCablesCombiTrans",
										        "SSDFlex1CombiTrans"};
  const char* ssdmodulerotname[7] = {"SSDStiffenerRotName",
										  "SSDChipRotName",
										"SSDSensorRotName",
										 "SSDFlex0RotName",
								  "SSDCoolingBlockRotName",
									"SSDChipCablesRotName",
										"SSDFlex1RotName"};
  TGeoCombiTrans *ssdmodulelocalmatrix[fgkSSDModuleCombiTransNumber];
  for(Int_t i=0; i<fgkSSDModuleCombiTransNumber; i++){ 
					ssdmodulelocalmatrix[i] =
					new TGeoCombiTrans(ssdmodulecombitransname[i],
									   kssdmoduletransx[i],
									   kssdmoduletransy[i],
									   kssdmoduletransz[i],
									   new TGeoRotation(ssdmodulerotname[i],
														      kssdmodulerotphi[i],
												        kssdmodulerottheta[i],
												        kssdmodulerotpsi[i]));
    fSSDModuleCombiTransMatrix[i] = ssdmodulelocalmatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetCarbonFiberJunctionCombiTransMatrix(){
  /////////////////////////////////////////////////////////////
  // Method generating CombiTrans Matrix for Carbon Fiber Junction   
  /////////////////////////////////////////////////////////////
  //Translation Parameters CarbonFiberJunction:
  ////////////////////////////////////////////////////////
  const Double_t kcarbonfiberjunctiontransx[3] = 
				{  0.0,fgkCarbonFiberTriangleLength,fgkCarbonFiberTriangleLength
				 * TMath::Cos(fgkCarbonFiberTriangleAngle*TMath::DegToRad())};
  const Double_t kcarbonfiberjunctiontransy[3] = 
									  {  0.0, 0.0,fgkCarbonFiberTriangleLength
				 *   TMath::Sin(fgkCarbonFiberTriangleAngle*TMath::DegToRad())};
  const Double_t kcarbonfiberjunctiontransz[3] = {  0.0,  0.0,  0.0};
  ////////////////////////////////////////////////////////
  //Rotational Parameters CarbonFiberJunction:
  ////////////////////////////////////////////////////////
  const Double_t kcarbonfiberjunctionrotphi[3]   = {   0., 120., 240.};
  const Double_t kcarbonfiberjunctionrottheta[3] = {   0.,   0.,   0.};
  const Double_t kcarbonfiberjunctionrotpsi[3]   = {   0.,   0.,   0.};
  ///////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of CarbonFiberJunction:
  ///////////////////////////////////////////////////////////
  char carbonfiberjunctioncombitransname[40];
  char carbonfiberjunctionrotname[40];
  TGeoCombiTrans *carbonfiberjunctionlocalmatrix[fgkCarbonFiberJunctionCombiTransNumber];
  for(Int_t i=0; i<fgkCarbonFiberJunctionCombiTransNumber; i++) {
		sprintf(carbonfiberjunctioncombitransname,"CarbonFiberJunctionCombiTrans%i",i);
		sprintf(carbonfiberjunctionrotname,"CarbonFiberJunctionRot%i",i);
		carbonfiberjunctionlocalmatrix[i] =
					new TGeoCombiTrans(carbonfiberjunctioncombitransname,
									   kcarbonfiberjunctiontransx[i],
									   kcarbonfiberjunctiontransy[i],
									   kcarbonfiberjunctiontransz[i],
								new TGeoRotation(carbonfiberjunctionrotname,
												kcarbonfiberjunctionrotphi[i],
												kcarbonfiberjunctionrottheta[i],
												kcarbonfiberjunctionrotpsi[i]));
    fCarbonFiberJunctionCombiTransMatrix[i] = carbonfiberjunctionlocalmatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetEndLadderCarbonFiberJunctionCombiTransMatrix(Int_t i){
  /////////////////////////////////////////////////////////////
  // Method generating CombiTrans Matrix for End Ladder Carbon Fiber Junction   
  /////////////////////////////////////////////////////////////
  //Translation Parameters EndLadderCarbonFiberJunction:
  ////////////////////////////////////////////////////////
  const Double_t kendladdercarbonfiberjunctiontransx[3] = 
				{  0.0,fgkCarbonFiberTriangleLength,fgkCarbonFiberTriangleLength
		*            TMath::Cos(fgkCarbonFiberTriangleAngle*TMath::DegToRad())};
  const Double_t kendladdercarbonfiberjunctiontransy[3] = 
										{  0.0, 0.0,fgkCarbonFiberTriangleLength
		*            TMath::Sin(fgkCarbonFiberTriangleAngle*TMath::DegToRad())};
  const Double_t kendladdercarbonfiberjunctiontransz[3] = {  0.0,  0.0,  
						   0.5*(fgkEndLadderCarbonFiberLowerJunctionLength[i]
							   -fgkEndLadderCarbonFiberUpperJunctionLength[i])};
  ////////////////////////////////////////////////////////
  //Rotational Parameters EndLadderCarbonFiberJunction:
  ////////////////////////////////////////////////////////
  const Double_t kendladdercarbonfiberjunctionrotphi[3]   = {   0., 120., 240.};
  const Double_t kendladdercarbonfiberjunctionrottheta[3] = {   0.,   0.,   0.};
  const Double_t kendladdercarbonfiberjunctionrotpsi[3]   = {   0.,   0.,   0.};
  ///////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of CarbonFiberJunction:
  ///////////////////////////////////////////////////////////
  char endladdercarbonfiberjunctioncombitransname[40];
  char endladdercarbonfiberjunctionrotname[40];
  TGeoCombiTrans *endladdercarbonfiberjunctionlocalmatrix[fgkEndLadderCarbonFiberJunctionCombiTransNumber];
  for(Int_t i=0; i<fgkEndLadderCarbonFiberJunctionCombiTransNumber; i++) {
	sprintf(endladdercarbonfiberjunctioncombitransname,"EndLadderCarbonFiberJunctionCombiTrans%i",i);
	sprintf(endladdercarbonfiberjunctionrotname,"EndLadderCarbonFiberJunctionRot%i",i);
	endladdercarbonfiberjunctionlocalmatrix[i] =
	new TGeoCombiTrans(endladdercarbonfiberjunctioncombitransname,
										  kendladdercarbonfiberjunctiontransx[i],
										  kendladdercarbonfiberjunctiontransy[i],
										  kendladdercarbonfiberjunctiontransz[i],
						new TGeoRotation(endladdercarbonfiberjunctionrotname,
									  	  kendladdercarbonfiberjunctionrotphi[i],
										    kendladdercarbonfiberjunctionrottheta[i],
									     kendladdercarbonfiberjunctionrotpsi[i]));
    fEndLadderCarbonFiberJunctionCombiTransMatrix[i] = 
									 endladdercarbonfiberjunctionlocalmatrix[i];
  }
}
////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetCarbonFiberAssemblyCombiTransMatrix(){
  /////////////////////////////////////////////////////////////
  // Method generating CombiTrans Matrix for Carbon fiber Assembly   
  /////////////////////////////////////////////////////////////
  //Translation Parameters CarbonFiberAssembly:
  ////////////////////////////////////////////////////////
  const Double_t kcarbonfiberassemblytransx[3] = {  0.0,  0.0,  0.0};
  const Double_t kcarbonfiberassemblytransy[3] = 
										{  0.5*fgkCarbonFiberJunctionWidth, 0.0, 
					fgkCarbonFiberJunctionWidth-fgkCarbonFiberLowerSupportWidth
			   -    fgkCarbonFiberLowerSupportVolumePosition[0]
			   -    fgkCarbonFiberLowerSupportVolumePosition[1]};
  const Double_t kcarbonfiberassemblytransz[3] = 
						  {  0.0,  0.0,-  0.5*fgkCarbonFiberLowerSupportHeight};
  ////////////////////////////////////////////////////////
  //Rotational Parameters CarbonFiberAssembly:
  ////////////////////////////////////////////////////////
  const Double_t kcarbonfiberassemblyrotphi[3]   = {   0.,  90.,   0.};
  const Double_t kcarbonfiberassemblyrottheta[3] = {  90.,
											-fgkCarbonFiberTriangleAngle,   0.};
  const Double_t kcarbonfiberassemblyrotpsi[3]   = {   0.,- 90.,   0.};
  ///////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of CarbonFiberAssembly:
  ///////////////////////////////////////////////////////////
  char carbonfiberassemblycombitransname[30];
  char carbonfiberassemblyrotname[30];
  TGeoCombiTrans *carbonfiberassemblylocalmatrix[fgkCarbonFiberAssemblyCombiTransNumber];
  for(Int_t i=0; i<fgkCarbonFiberAssemblyCombiTransNumber; i++) {
	sprintf(carbonfiberassemblycombitransname,"CarbonFiberAssemblyCombiTrans%i",i);
	sprintf(carbonfiberassemblyrotname,"CarbonFiberAssemblyRot%i",i);
	carbonfiberassemblylocalmatrix[i] =
						new TGeoCombiTrans(carbonfiberassemblycombitransname,
										   kcarbonfiberassemblytransx[i],
										   kcarbonfiberassemblytransy[i],
										   kcarbonfiberassemblytransz[i],
						 new TGeoRotation(carbonfiberassemblyrotname,
										   kcarbonfiberassemblyrotphi[i],
										   kcarbonfiberassemblyrottheta[i],
										   kcarbonfiberassemblyrotpsi[i]));
    fCarbonFiberAssemblyCombiTransMatrix[i] = carbonfiberassemblylocalmatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetCoolingTubeSupportCombiTransMatrix(){
  /////////////////////////////////////////////////////////////
  // Method generating CombiTrans Matrix for Cooling Tube Support   
  /////////////////////////////////////////////////////////////
  //Translation Parameters CoolingTubeSupport:
  ////////////////////////////////////////////////////////
  Double_t phi = TMath::ASin(0.5*fgkCoolingTubeSupportHeight
													/fgkCoolingTubeSupportRmax);
  const Double_t kcoolingtubesupporttransx[2] = 
						  {0.,2.*fgkCoolingTubeSupportRmax*TMath::Cos(phi)
						+  2.*(fgkCoolingTubeSupportLength
						-  fgkCoolingTubeSupportRmax*(1.+TMath::Cos(phi)))
						+  fgkCarbonFiberTriangleLength
						-  2.0*fgkCarbonFiberJunctionLength};
  const Double_t kcoolingtubesupporttransy[2] = {  0.0,  0.0};
  const Double_t kcoolingtubesupporttransz[2] = {  0.0,  0.0};
  ////////////////////////////////////////////////////////
  //Rotational Parameters CoolingTubeSupport:
  ////////////////////////////////////////////////////////
  const Double_t kcoolingsubesupportrotphi[2]   = {   0., 180.};
  const Double_t kcoolingsubesupportrottheta[2] = {   0.,   0.};
  const Double_t kcoolingsubesupportrotpsi[2]   = {   0.,   0.};
  ///////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of CarbonFiberJunction:
  ///////////////////////////////////////////////////////////
  char coolingtubesupportcombitransname[40];
  char coolingtubesupportrotname[40];
  TGeoCombiTrans *coolingtubesupportlocalmatrix[fgkCoolingTubeSupportCombiTransNumber];
  for(Int_t i=0; i<fgkCoolingTubeSupportCombiTransNumber; i++) {
	sprintf(coolingtubesupportcombitransname,"CoolingTubeSupportCombiTrans%i",i);
	sprintf(coolingtubesupportrotname,"CoolingTubeSupportRot%i",i);
	coolingtubesupportlocalmatrix[i] =
			new TGeoCombiTrans(coolingtubesupportcombitransname,
							   kcoolingtubesupporttransx[i],
							   kcoolingtubesupporttransy[i],
							   kcoolingtubesupporttransz[i],
							   new TGeoRotation(coolingtubesupportrotname,
												kcoolingsubesupportrotphi[i],
												kcoolingsubesupportrottheta[i],
												kcoolingsubesupportrotpsi[i]));
    fCoolingTubeSupportCombiTransMatrix[i] = coolingtubesupportlocalmatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetCoolingTubeCombiTransMatrix(){
  /////////////////////////////////////////////////////////////
  // Method generating CombiTrans Matrix for Cooling Tube  
  /////////////////////////////////////////////////////////////
  //Translation Parameters CoolingTube:
  ////////////////////////////////////////////////////////
  const Double_t kcoolingtubetransx[2] = {  0.,  fgkCoolingTubeSeparation};
  const Double_t kcoolingtubetransy[2] = {  fgkCoolingTubeLength/2.0, fgkCoolingTubeLength/2.0};
  const Double_t kcoolingtubetransz[2] = {  0.0,  0.};
  ////////////////////////////////////////////////////////
  //Rotational Parameters CoolingTube:
  ////////////////////////////////////////////////////////
  const Double_t kcoolingtuberotphi[2]   = {   0.,   0.};
  const Double_t kcoolingtuberottheta[2] = {  90.,  90.};
  const Double_t kcoolingtuberotpsi[2]   = {   0.,   0.};
  ///////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of CarbonFiberJunction:
  ///////////////////////////////////////////////////////////
  const char* coolingtubecombitransname[fgkCoolingTubeCombiTransNumber] = 
							{"CoolingTubeCombiTrans0","CoolingTubeCombiTrans1"};
  const char* coolingtuberotationname[fgkCoolingTubeCombiTransNumber] = 
								{"CoolingTubeRotation0","CoolingTubeRotation1"};
  TGeoCombiTrans *coolingtubelocalmatrix[fgkCoolingTubeCombiTransNumber];
  for(Int_t i=0; i<fgkCoolingTubeCombiTransNumber; i++) {
	 coolingtubelocalmatrix[i] =
			new TGeoCombiTrans(coolingtubecombitransname[i],
							   kcoolingtubetransx[i],
							   kcoolingtubetransy[i],
							   kcoolingtubetransz[i], 
							   new TGeoRotation(coolingtuberotationname[i],
							                    kcoolingtuberotphi[i],
							                    kcoolingtuberottheta[i],
                           kcoolingtuberotpsi[i]) );
    fCoolingTubeTransMatrix[i] = coolingtubelocalmatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetLadderSegmentCombiTransMatrix(){
  /////////////////////////////////////////////////////////////
  // Method generating CombiTrans Matrix for SSD Ladder Segment   
  /////////////////////////////////////////////////////////////
  //Translation Parameters LadderSegment:
  ////////////////////////////////////////////////////////
	const Double_t kladdersegmenttransx[fgkLadderSegmentCombiTransNumber] = {  0.,
	 -  0.5*(fgkSSDSensorWidth-fgkCarbonFiberTriangleLength),
			 fgkCarbonFiberTriangleLength+fgkCarbonFiberJunctionToSensorSupport,
			 fgkCarbonFiberJunctionLength-(fgkCoolingTubeSupportLength
	 -       fgkCoolingTubeSupportRmax),
		0.5*(fgkCarbonFiberTriangleLength-fgkCoolingTubeSeparation)}; 
	const Double_t kladdersegmenttransy[fgkLadderSegmentCombiTransNumber] = {  0.,
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
	const Double_t kladdersegmenttransz[fgkLadderSegmentCombiTransNumber] = {  0.,
	 -        (fgkSSDModuleCoolingBlockToSensor+0.5*fgkCoolingTubeSupportHeight
	 -         fgkSSDSensorHeight-fgkSSDChipCablesHeight[3]-fgkSSDChipHeight),
																			 0.,
	 -     0.5*fgkCoolingTubeSupportHeight,
	 -     0.5*fgkCoolingTubeSupportHeight};
//////////////////////////////////////////////////
  //Rotational Parameters LadderSegment:
  ////////////////////////////////////////////////////////
  const Double_t kladdersegmentrotphi[fgkLadderSegmentCombiTransNumber]   =
													  {   0.,   0.,- 90.,   0.,  0.};
  const Double_t kladdersegmentrottheta[fgkLadderSegmentCombiTransNumber] = 
													  {   0.,   0.,   0.,  90.,  0.};
  const Double_t kladdersegmentrotpsi[fgkLadderSegmentCombiTransNumber]   = 
													  {   0.,   0.,   0.,   0.,  0.};
  //////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of LadderSegment:
  //////////////////////////////////////////////////////
  char laddersegmentcombitransname[40];
  char laddersegmentrotname[40];
  TGeoCombiTrans *laddersegmentlocalmatrix[fgkLadderSegmentCombiTransNumber];
  for(Int_t i=0; i<fgkLadderSegmentCombiTransNumber; i++) {
		sprintf(laddersegmentcombitransname,"LadderSegmentCombiTrans%i",i);
		sprintf(laddersegmentrotname,"LadderSegmentRot%i",i);
		laddersegmentlocalmatrix[i] =
					new TGeoCombiTrans(laddersegmentcombitransname,
									   kladdersegmenttransx[i],
									   kladdersegmenttransy[i],
									   kladdersegmenttransz[i],
									   new TGeoRotation(laddersegmentrotname,
														kladdersegmentrotphi[i],
														kladdersegmentrottheta[i],
														kladdersegmentrotpsi[i]));
    fLadderSegmentCombiTransMatrix[i] = laddersegmentlocalmatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetEndLadderSegmentCombiTransMatrix(Int_t i){
  /////////////////////////////////////////////////////////////
  // Method generating CombiTrans Matrix for SSD End Ladder Segment   
  /////////////////////////////////////////////////////////////
  //Translation Parameters EndLadderSegment:
  ////////////////////////////////////////////////////////
  const Double_t kendladdersegmenttransx[fgkEndLadderSegmentCombiTransNumber] =
            {0.0,
											  0.0,
          -  0.25*(fgkSSDMountingBlockLength[0]
			       +	 fgkSSDMountingBlockLength[1])
          +  0.5*fgkCarbonFiberTriangleLength,
											  0.0}; 
  const Double_t kendladdersegmenttransy[fgkEndLadderSegmentCombiTransNumber] = 
							 {0.5*fgkEndLadderCarbonFiberLowerJunctionLength[i],
            i==0 ? 0. : fgkCarbonFiberLowerSupportWidth,
                        fgkEndLadderMountingBlockPosition[i],
					           (1-i)*(fgkEndLadderMountingBlockPosition[i]
										 +  0.5*fgkSSDMountingBlockWidth)}; 
  const Double_t kendladdersegmenttransz[fgkEndLadderSegmentCombiTransNumber] = 
            {0.0,
											  0.0,
										-  fgkSSDMountingBlockHeight[1]
										+  0.5*fgkSSDMountingBlockHeight[0],
          -  0.5*fgkCarbonFiberLowerSupportHeight}; 
  ////////////////////////////////////////////////////////
  //Rotational Parameters EndLadderSegment:
  ////////////////////////////////////////////////////////  
  const Double_t kendladdersegmentrotphi[fgkEndLadderSegmentCombiTransNumber] =  
            {   0.,  90.,   0.,   0.};
  const Double_t kendladdersegmentrottheta[fgkEndLadderSegmentCombiTransNumber] = 
            {  90.,-fgkCarbonFiberTriangleAngle, 0.,   0.};
  const Double_t kendladdersegmentrotpsi[fgkEndLadderSegmentCombiTransNumber] = 
            {   0.,- 90.,   0.,   0.};
  ////////////////////////////////////////////////////////
  //Name of CombiTrans Transformation of EndLadderSegment:
  ////////////////////////////////////////////////////////
  char endladdersegmentcombitransname[30];
  char endladdersegmentrotname[30];
  TGeoCombiTrans *endladdersegmentlocalmatrix[fgkEndLadderSegmentCombiTransNumber];
  for(Int_t i=0; i<fgkEndLadderSegmentCombiTransNumber; i++){ 
  		sprintf(endladdersegmentcombitransname,"EndLadderSegmentCombiTrans%i",i);
		sprintf(endladdersegmentrotname,"EndLadderSegmentRot%i",i);
		endladdersegmentlocalmatrix[i] =
				new TGeoCombiTrans(endladdersegmentcombitransname,
								   kendladdersegmenttransx[i],
								   kendladdersegmenttransy[i],
								   kendladdersegmenttransz[i],
								   new TGeoRotation(endladdersegmentrotname,
													kendladdersegmentrotphi[i],
													kendladdersegmentrottheta[i],
													kendladdersegmentrotpsi[i]));
    fEndLadderSegmentCombiTransMatrix[i] = endladdersegmentlocalmatrix[i];
  }
}
/////////////////////////////////////////////////////////////////////////////////
void AliITSv11GeometrySSD::SetLadderCableCombiTransMatrix(Int_t iLayer){
  /////////////////////////////////////////////////////////////
  // Method generating CombiTrans Matrix for SSD Ladder Cable   
  /////////////////////////////////////////////////////////////
  // Translation Parameters for LadderCable
  /////////////////////////////////////////
  Double_t ssdflexradius = fgkSSDStiffenerHeight+2*fgkSSDFlexHeight[0]
					     + fgkSSDFlexHeight[1];  
  Double_t ssdflexradiusmax = (fgkSSDFlexLength[3]-fgkSSDFlexLength[2])
							/  TMath::Tan(fgkSSDFlexAngle*TMath::DegToRad());
  Double_t ssdladdercabletransx[3];
  ssdladdercabletransx[0] = (ssdflexradiusmax-fgkSSDFlexHeight[1]-ssdflexradius)
						  *   TMath::Sin(2.*fgkSSDFlexAngle*TMath::DegToRad())
						  *	  TMath::Cos(2.*fgkSSDFlexAngle*TMath::DegToRad());
  ssdladdercabletransx[1] = ((ssdflexradiusmax-fgkSSDFlexHeight[1]-ssdflexradius)
						  -     ssdladdercabletransx[0]
						  /     TMath::Sin(2.*fgkSSDFlexAngle*TMath::DegToRad()))
						  *     TMath::Cos(fgkSSDFlexAngle*TMath::DegToRad());						
  ssdladdercabletransx[2] = (fgkSSDFlexFullLength-2.*fgkSSDFlexAngle
						  *	  TMath::DegToRad()*ssdflexradiusmax
						  -     fgkSSDFlexLength[2]-TMath::Pi()
						  *	  fgkSSDStiffenerHeight-fgkSSDFlexLength[0]
						  -	  fgkSSDLadderCableWidth)
						  *	  TMath::Cos(2.*fgkSSDFlexAngle*TMath::DegToRad());
  Double_t ssdladdercabletransz[3] = {ssdladdercabletransx[0]
						  *	TMath::Tan(2.*fgkSSDFlexAngle*TMath::DegToRad()),
							ssdladdercabletransx[1]
						  *	TMath::Tan(fgkSSDFlexAngle*TMath::DegToRad()),
							ssdladdercabletransx[2]
						  *	TMath::Tan(2.*fgkSSDFlexAngle*TMath::DegToRad())};	
  TGeoRotation* localladdercablerot[3];	
  localladdercablerot[0] = new TGeoRotation("LocalLadderCableRot0",90.,0.,0.);
  localladdercablerot[1] = new TGeoRotation("LocalLadderCableRot1",90.,60.,-90.);
  localladdercablerot[2] = new TGeoRotation("LocalLadderCableRot2");
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
  localladdercablecombitransmatrix[0][1] = fLadderSegmentCombiTransMatrix[1];
  localladdercablecombitransmatrix[0][2] = 
						new TGeoCombiTrans(fgkSSDModuleStiffenerPosition[0],
										   fgkSSDModuleStiffenerPosition[1],0.,0);
  localladdercablecombitransmatrix[0][3] = fSSDModuleCombiTransMatrix[6];
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
  for(Int_t i=0; i<klocalladdercombitransnumber; i++)
   localladdercablecombitransmatrix[1][i] = 
			(i!=3 ? localladdercablecombitransmatrix[0][i]:
					fSSDModuleCombiTransMatrix[3]); 	
  ///////////////////////////////////////////
  // Setting LadderCableCombiTransMatrix
  ///////////////////////////////////////////
  Int_t beamaxistrans[3] = {0,0,0};
  switch(iLayer){
  case 5: 
	  beamaxistrans[0] = fgkSSDLay5SensorsNumber/2;
	  beamaxistrans[1] = beamaxistrans[0]+1;
	  beamaxistrans[2] = beamaxistrans[0]-1;
  break;
  case 6:
	  beamaxistrans[0] = (fgkSSDLay6SensorsNumber-1)/2;
  	beamaxistrans[1] = beamaxistrans[0]+1;
  	beamaxistrans[2] = beamaxistrans[0];
	  break;
  }
  TGeoHMatrix* localladdercablehmatrix[klocalladdersidecablesnumber];
  TGeoRotation* laddercablerot;
  TGeoTranslation* laddercabletrans;
  for(Int_t i=0; i<klocalladdersidecablesnumber; i++){
	 localladdercablehmatrix[i] = new TGeoHMatrix();
	 for(Int_t j=0; j<klocalladdercombitransnumber; j++){
		   localladdercablehmatrix[i]->MultiplyLeft(
		   localladdercablecombitransmatrix[i][klocalladdercombitransnumber-j-1]);
	 }
	 laddercablerot = new TGeoRotation();
	 laddercablerot->SetMatrix(localladdercablehmatrix[i]->GetRotationMatrix());
  laddercabletrans = new TGeoTranslation();
  Double_t* laddercabletransvector = localladdercablehmatrix[i]->GetTranslation();
  laddercabletrans->SetTranslation(laddercabletransvector[0],
									 laddercabletransvector[1]
					+                (i==0 ? beamaxistrans[0] : 0.)
					*				 fgkCarbonFiberJunctionWidth,
									 laddercabletransvector[2]);	 
	 fLadderCableCombiTransMatrix[i] = new TGeoCombiTrans(*laddercabletrans,
												   *laddercablerot);
	} 
	fLadderCableCombiTransMatrix[2] = 
					AddTranslationToCombiTrans(fLadderCableCombiTransMatrix[1],0.,
					beamaxistrans[1]*fgkCarbonFiberJunctionWidth,0.);
	fLadderCableCombiTransMatrix[3] = 
					AddTranslationToCombiTrans(fLadderCableCombiTransMatrix[0],0.,
					beamaxistrans[2]*fgkCarbonFiberJunctionWidth,0.);
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDSensorSupportShape(Double_t length, 
							Double_t height,Double_t width,Double_t* thickness){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Sensor Support Shape   
  /////////////////////////////////////////////////////////////
  const Int_t kvertexnumber = 4;
  const Int_t kshapesnumber = 2;
  Double_t shapewidth[2] = {width,width};
  TVector3** vertexposition[kshapesnumber];
  for(Int_t i = 0; i<kshapesnumber; i++) vertexposition[i] = 
													new TVector3*[kvertexnumber];
  //First Shape Vertex Positioning
  vertexposition[0][0] = new TVector3();
  vertexposition[0][1] = new TVector3(height);
  vertexposition[0][2] = new TVector3(thickness[0]);
  vertexposition[0][3] = new TVector3(*vertexposition[0][1]);
  //Second Shape Vertex Positioning
  vertexposition[1][0] = new TVector3(*vertexposition[0][0]);
  vertexposition[1][1] = new TVector3(length);
  vertexposition[1][2] = new TVector3(thickness[1]);
  vertexposition[1][3] = new TVector3(vertexposition[1][1]->X());
  char* ssdsensorsupportshapename[kshapesnumber] = 
							{"SSDSensorSupportShape1","SSDSensorSupportShape2"};
  TGeoArb8* lSSDSensorSupportShape[kshapesnumber];
  for(Int_t i = 0; i< kshapesnumber; i++) lSSDSensorSupportShape[i] = 
					   GetArbShape(vertexposition[i],shapewidth,i==0 ? 
													 thickness[1]: thickness[0],
												  ssdsensorsupportshapename[i]);
  /////////////////////////////////////
  //Setting Translations and Rotations: 
  /////////////////////////////////////
  TGeoRotation* ssdsensorsupportshaperot[2];
  ssdsensorsupportshaperot[0] = 
					 new TGeoRotation("SSDSensorSupportShapeRot1",180.,0.,0.);
  ssdsensorsupportshaperot[1] = 
					 new TGeoRotation("SSDSensorSupportShapeRot2",90.,90.,-90.);
  TGeoTranslation* ssdsensorsupportshapetrans = 
					 new TGeoTranslation("SSDSensorSupportShapeTrans",0.,0.,
															  0.5*thickness[0]);
  TGeoCombiTrans* ssdsensorsupportcombitrans = 
	  new TGeoCombiTrans("SSDSensorSupportCombiTrans",0.5*thickness[0],width,0.,
						  new TGeoRotation((*ssdsensorsupportshaperot[1])
						*                  (*ssdsensorsupportshaperot[0])));
  TGeoVolume* ssdsensorsupportcompvolume = 
						  new TGeoVolumeAssembly("SSDSensorSupportCompVolume");
  ssdsensorsupportcompvolume->AddNode(new TGeoVolume("SSDSensorSupportVolume1",
							lSSDSensorSupportShape[0],fSSDSensorSupportMedium),1,
							ssdsensorsupportshapetrans);
  ssdsensorsupportcompvolume->AddNode(new TGeoVolume("SSDSensorSupportVolume2",
							lSSDSensorSupportShape[1],fSSDSensorSupportMedium),1,
							ssdsensorsupportcombitrans);
  return ssdsensorsupportcompvolume;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDSensorSupport(Int_t VolumeKind, Int_t n){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Sensor Support    
  /////////////////////////////////////////////////////////////
  TGeoVolume* ssdsensorsupport;
  Double_t sidesupporthickness[2] = {fgkSSDSensorSideSupportThickness[0],
										   fgkSSDSensorSideSupportThickness[1]};
  VolumeKind == 0 ? ssdsensorsupport = GetSSDSensorSupportShape(
								fgkSSDSensorSideSupportLength,
								fgkSSDSensorSideSupportHeight[(n==0 ? 0 : 1)],
								fgkSSDSensorSideSupportWidth,
								sidesupporthickness) :
    ssdsensorsupport = GetSSDSensorSupportShape(fgkSSDSensorCenterSupportLength,
						        fgkSSDSensorCenterSupportHeight[(n==0 ? 0 : 1)],
							    fgkSSDSensorCenterSupportWidth,
							    sidesupporthickness);
  return ssdsensorsupport;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDSensorSupportAssembly(Int_t n){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Sensor Support Assembly   
  /////////////////////////////////////////////////////////////
  TGeoVolume* ssdsensorsupportassembly = 
							 new TGeoVolumeAssembly("SSDSensorSupportAssembly");
  const Int_t kvolumenumber = 2;
  TGeoVolume* ssdsensorsupport[kvolumenumber];
  for(Int_t i=0; i<kvolumenumber; i++) ssdsensorsupport[i] = 
														GetSSDSensorSupport(i,n);
  SetSSDSensorSupportCombiTransMatrix();
  for(Int_t i=0; i<fgkSSDSensorSupportCombiTransNumber; i++) 
	ssdsensorsupportassembly->AddNode((i<2 ? ssdsensorsupport[0]:
											 ssdsensorsupport[1]),
										 i+1,fSSDSensorSupportCombiTransMatrix[i]);
  return ssdsensorsupportassembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDModule(Int_t iChipCablesHeight){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Sensor Module    
  /////////////////////////////////////////////////////////////
  TGeoVolume* ssdmodulevolume[fgkSSDModuleCombiTransNumber-1];
  ssdmodulevolume[0] = GetSSDStiffenerAssembly();
  ssdmodulevolume[1] = GetSSDChipAssembly();
  ssdmodulevolume[2] = GetSSDSensor();
  ssdmodulevolume[3] = GetSSDFlexAssembly();
  ssdmodulevolume[4] = GetSSDCoolingBlockAssembly();
  ssdmodulevolume[5] = GetSSDChipCablesAssembly(fgkSSDChipCablesHeight[iChipCablesHeight+2]);
  SetSSDModuleCombiTransMatrix(fgkSSDChipCablesHeight[iChipCablesHeight+2]);
  TGeoCombiTrans* ssdmoduleglobalcombitrans = 
								   new TGeoCombiTrans("SSDModuleGlobalCombiTrans",
									   fgkSSDModuleStiffenerPosition[0],
									   fgkSSDModuleStiffenerPosition[1],0.,NULL);
  TGeoHMatrix* ssdmodulehmatrix[fgkSSDModuleCombiTransNumber];
  TGeoVolume* ssdmodule = new TGeoVolumeAssembly("SSDModule");
  for(Int_t i=0; i<fgkSSDModuleCombiTransNumber; i++){ 
	ssdmodulehmatrix[i] = new TGeoHMatrix((*ssdmoduleglobalcombitrans)
						*				  (*fSSDModuleCombiTransMatrix[i]));
	ssdmodule->AddNode(i==fgkSSDModuleCombiTransNumber-1 ? ssdmodulevolume[3] : 
						  ssdmodulevolume[i],i!=fgkSSDModuleCombiTransNumber-1?1:2,
						  ssdmodulehmatrix[i]);
  }
  return ssdmodule;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDSensor(){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Sensor  
  /////////////////////////////////////////////////////////////
  Double_t ssdsensitivelength = fgkSSDSensorLength-2.*fgkSSDSensorInsensitiveLength;
  Double_t ssdsensitivewidth  = fgkSSDSensorWidth-2.*fgkSSDSensorInsensitiveWidth;
  TGeoBBox* ssdsensorsensitiveshape = new TGeoBBox("SSDSensorSensitiveShape",
                                                0.5*ssdsensitivelength,
                                                0.5*ssdsensitivewidth,
                                                0.5*fgkSSDSensorHeight);
  TGeoVolume* ssdsensorsensitive = 
      new TGeoVolume(fgkSSDSensitiveVolName,ssdsensorsensitiveshape,fSSDSensorMedium);
  ssdsensorsensitive->SetLineColor(fColorSilicon);
  TGeoBBox* ssdsensorinsensitiveshape[2];
  ssdsensorinsensitiveshape[0] = new TGeoBBox("SSDSensorInsensitiveShape1",
                                                0.5*fgkSSDSensorLength,
                                                0.5*fgkSSDSensorInsensitiveWidth,
                                                0.5*fgkSSDSensorHeight);
  ssdsensorinsensitiveshape[1] = new TGeoBBox("SSDSensorInsensitiveShape2",
                                                0.5*fgkSSDSensorInsensitiveWidth,
                                                0.5*ssdsensitivewidth,
                                                0.5*fgkSSDSensorHeight);
  const char* ssdsensorinsensitivename[2] = {"SSDSensorInsensitive1",
                                             "SSDSensorInsensitive2"};
  TGeoVolume* ssdsensorinsensitive[2];
  for(Int_t i=0; i<2; i++){ ssdsensorinsensitive[i] = 
      new TGeoVolume(ssdsensorinsensitivename[i],ssdsensorinsensitiveshape[i],
                     fSSDSensorMedium);
      ssdsensorinsensitive[i]->SetLineColor(fColorCarbonFiber);
  }
  TGeoVolume* ssdsensorinsensitivevol = 
                              new TGeoVolumeAssembly("SSDSensorInsensitiveVol");
  for(Int_t i=0; i<4; i++) 
            ssdsensorinsensitivevol->AddNode(i%2==0 ? ssdsensorinsensitive[0]:
            ssdsensorinsensitive[1],i<2?1:2,
                  new TGeoTranslation(0.5*(1.-TMath::Power(-1.,i))*(i==1? 1.:-1.)
      *     (ssdsensorsensitiveshape->GetDX()+ssdsensorinsensitiveshape[1]->GetDX()),
                        0.5*(1.+TMath::Power(-1.,i))*(i==0?-1.: 1.)
      *     (ssdsensorsensitiveshape->GetDY()+ssdsensorinsensitiveshape[0]->GetDY()),0.));            
  TGeoVolume* ssdsensor = new TGeoVolumeAssembly("SSDSensor");
  ssdsensor->AddNode(ssdsensorsensitive,1),ssdsensor->AddNode(ssdsensorinsensitivevol,1);
  return ssdsensor;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDChipAssembly() const{
  /////////////////////////////////////////////////////////////
  // Method generating SSD Chip Assembly    
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
  TGeoVolume* ssdchip = new TGeoVolumeAssembly("SSDChip");  
  for(Int_t i=0; i<2; i++) ssdchip->AddNode(ssdchipcomp[i],1,ssdchipcomptrans[i]);
  Double_t ssdchipseparation[2] = {fgkSSDChipLength+fgkSSDChipSeparationLength,
						  fgkSSDSensorLength-2.*fgkSSDModuleStiffenerPosition[1]
				   -  2.*(fgkSSDStiffenerWidth-fgkSSDStiffenerToChipDist
				   -  0.5*fgkSSDChipWidth)};
  TGeoVolume* ssdchipassembly = new TGeoVolumeAssembly("SSDChipAssembly"); 
  for(Int_t i=0; i<kssdchiprownumber; i++)
    for(Int_t j=0; j<fgkSSDChipNumber; j++) 
		ssdchipassembly->AddNode(ssdchip,fgkSSDChipNumber*i+j+1,
		new TGeoTranslation(j*ssdchipseparation[0],i*ssdchipseparation[1],0.));
  return ssdchipassembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDStiffenerAssembly(){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Stiffener Assembly    
  /////////////////////////////////////////////////////////////
  const Int_t kssdstiffenernumber = 2;
  Double_t ssdstiffenerseparation = fgkSSDSensorLength
								  - 2.*fgkSSDModuleStiffenerPosition[1]
								  -    fgkSSDStiffenerWidth;
  TGeoVolume* ssdstiffener = new TGeoVolumeAssembly("SSDStiffener");
////////////////////////////
// Stiffener Volumes
///////////////////////////
  const Int_t kstiffenerboxnumber = 6;
  TGeoBBox* ssdstiffenerboxshapes[kstiffenerboxnumber];
  ssdstiffenerboxshapes[0] = new TGeoBBox("SSDStiffenerBoxShape1",
										  0.5* fgkSSDStiffenerLength,
										  0.5* fgkSSDStiffenerWidth,
										  0.5*(fgkSSDStiffenerHeight
						   -                   fgkSSDConnectorHeight));
  ssdstiffenerboxshapes[1] = new TGeoBBox("SSDStiffenerBoxShape2",
										  0.5*(fgkSSDConnectorPosition[0]
						   -              2.0* fgkSSDConnectorLength
						   -				   fgkSSDConnectorSeparation),
										  0.5* fgkSSDStiffenerWidth,
										  0.5* fgkSSDConnectorHeight);
  ssdstiffenerboxshapes[2] = new TGeoBBox("SSDStiffenerBoxShape3",
										  0.5*(fgkSSDConnectorSeparation
						   +              2.*  fgkSSDConnectorLength),
										  0.5* fgkSSDConnectorPosition[1],
										  0.5* fgkSSDConnectorHeight);
  ssdstiffenerboxshapes[3] = new TGeoBBox("SSDStiffenerBoxShape4",
											ssdstiffenerboxshapes[2]->GetDX(),
										  0.5*(fgkSSDStiffenerWidth
						   -                   fgkSSDConnectorPosition[1]
						   -                   fgkSSDConnectorWidth),
										  0.5* fgkSSDConnectorHeight);
   ssdstiffenerboxshapes[4] = new TGeoBBox("SSDStiffenerBoxShape5",
										  0.5* fgkSSDConnectorSeparation,
										  0.5* fgkSSDConnectorWidth,
										  0.5* fgkSSDConnectorHeight);
   ssdstiffenerboxshapes[5] = new TGeoBBox("SSDStiffenerBoxShape6",
										  0.5*(fgkSSDStiffenerLength
							 -				   fgkSSDConnectorPosition[0]),
										  0.5* fgkSSDStiffenerWidth,
										  0.5* fgkSSDConnectorHeight);
  TGeoVolume* ssdstiffenerbox[kstiffenerboxnumber];
  char ssdtiffenerboxname[30];
  for(Int_t i=0; i<kstiffenerboxnumber; i++){ 
	 sprintf(ssdtiffenerboxname,"SSDStiffenerBox%i",i+1);
    ssdstiffenerbox[i] = new TGeoVolume(ssdtiffenerboxname,ssdstiffenerboxshapes[i],
									  fSSDStiffenerMedium);
	ssdstiffenerbox[i]->SetLineColor(fColorStiffener);
  }
////////////////////////////
// Connector 
///////////////////////////
  TGeoBBox* ssdconnectorshape =  new TGeoBBox("SSDConnectorShape",
											 0.5*fgkSSDConnectorLength,
											 0.5*fgkSSDConnectorWidth,
											 0.5*fgkSSDConnectorHeight);
  TGeoVolume* ssdconnector    = new TGeoVolume("SSDConnector",ssdconnectorshape,
												fSSDStiffenerConnectorMedium); 
  ssdconnector->SetLineColor(fColorAl);
  const Int_t kssdconnectornumber = 2;
  TGeoTranslation* ssdconnectortrans[kssdconnectornumber];
  ssdconnectortrans[0] = new TGeoTranslation("SSDConnectorTrans1",
        -  ssdstiffenerboxshapes[0]->GetDX()+fgkSSDConnectorPosition[0]
		-  fgkSSDConnectorSeparation-1.5*fgkSSDConnectorLength,
		   ssdstiffenerboxshapes[0]->GetDY()-fgkSSDConnectorPosition[1]
		-  ssdconnectorshape->GetDY(),
		   ssdstiffenerboxshapes[0]->GetDZ()+ssdconnectorshape->GetDZ());	
  ssdconnectortrans[1] = new TGeoTranslation("SSDConnectorTrans2",
        -  ssdstiffenerboxshapes[0]->GetDX()+fgkSSDConnectorPosition[0]
        -  0.5*fgkSSDConnectorLength,
           ssdstiffenerboxshapes[0]->GetDY()-fgkSSDConnectorPosition[1]
		-  ssdconnectorshape->GetDY(),
		   ssdstiffenerboxshapes[0]->GetDZ()+ssdconnectorshape->GetDZ());
  for(Int_t i=0; i<kssdconnectornumber; i++)
			ssdstiffener->AddNode(ssdconnector,i+1,ssdconnectortrans[i]);	
//////////////////////////////////////
// TGeoTranslation for Stiffener Boxes
//////////////////////////////////////
  TGeoTranslation* ssdstiffenerboxtrans[kstiffenerboxnumber];
  ssdstiffenerboxtrans[0] = new TGeoTranslation("SSDStiffenerBoxTrans1",0.,0.,0.);
  ssdstiffenerboxtrans[1] = new TGeoTranslation("SSDStiffenerBoxTrans2",
		 - (ssdstiffenerboxshapes[0]->GetDX()-ssdstiffenerboxshapes[1]->GetDX()),
		    0.,
			ssdstiffenerboxshapes[0]->GetDZ()+ssdstiffenerboxshapes[1]->GetDZ());
  ssdstiffenerboxtrans[2] = new TGeoTranslation("SSDStiffenerBoxTrans3",
         - ssdstiffenerboxshapes[0]->GetDX()-ssdstiffenerboxshapes[2]->GetDX()
		 + fgkSSDConnectorPosition[0],
		   ssdstiffenerboxshapes[0]->GetDY()-ssdstiffenerboxshapes[2]->GetDY(),
		   ssdstiffenerboxshapes[0]->GetDZ()+ssdstiffenerboxshapes[2]->GetDZ());
  ssdstiffenerboxtrans[3] = new TGeoTranslation("SSDStiffenerBoxTrans4",
         - ssdstiffenerboxshapes[0]->GetDX()-ssdstiffenerboxshapes[3]->GetDX()
		 + fgkSSDConnectorPosition[0],
		 - ssdstiffenerboxshapes[0]->GetDY()+ssdstiffenerboxshapes[3]->GetDY(),
		   ssdstiffenerboxshapes[0]->GetDZ()+ssdstiffenerboxshapes[3]->GetDZ());
  ssdstiffenerboxtrans[4] = new TGeoTranslation("SSDStiffenerBoxTrans5",
		 - ssdstiffenerboxshapes[0]->GetDX()+fgkSSDConnectorPosition[0]
		 - 0.5*fgkSSDConnectorSeparation-2.*ssdconnectorshape->GetDX(),
		   ssdstiffenerboxshapes[0]->GetDY()-fgkSSDConnectorPosition[1]
		 - ssdconnectorshape->GetDY(),
		   ssdstiffenerboxshapes[0]->GetDZ()+ssdconnectorshape->GetDZ());
  ssdstiffenerboxtrans[5] = new TGeoTranslation("SSDStiffenerBoxTrans6",
         - ssdstiffenerboxshapes[0]->GetDX()+fgkSSDConnectorPosition[0]
		 + ssdstiffenerboxshapes[5]->GetDX(),
		   0.,
		   ssdstiffenerboxshapes[0]->GetDZ()+ssdstiffenerboxshapes[5]->GetDZ());
  for(Int_t i=0; i<kstiffenerboxnumber; i++) 
		ssdstiffener->AddNode(ssdstiffenerbox[i],1,ssdstiffenerboxtrans[i]);	
  TGeoCombiTrans* ssdstiffenercombitrans[kssdstiffenernumber];
  char ssdstiffenercombitransname[30];
    for(Int_t i=0; i<kssdstiffenernumber; i++){ 
	sprintf(ssdstiffenercombitransname,"SSDStiffenerCombiTrans%i",i+1);
    ssdstiffenercombitrans[i] = new TGeoCombiTrans(ssdstiffenercombitransname,
			0.,i*ssdstiffenerseparation,0.,new TGeoRotation("",180*(1-i),0.,0.));
  }
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
////////////////////////////
// Capacitor 1812-330 nF
///////////////////////////
  TGeoBBox* capacitor1812shape =  new TGeoBBox("Capacitor1812Shape",
											 0.5*fgkSSDCapacitor1812Length,
											 0.5*fgkSSDCapacitor1812Width,
											 0.5*fgkSSDCapacitor1812Height);
  TGeoVolume* capacitor1812 = new TGeoVolume("Capacitor1812",capacitor1812shape,
                                             fSSDStiffener1812CapacitorMedium); 
  capacitor1812->SetLineColor(fColorAl);
  TGeoTranslation* capacitor1812trans = new TGeoTranslation("Capacitor1812Trans",
			  0.,
			  0.5*fgkSSDStiffenerWidth+ssdstiffenerseparation
	       -  capacitor1812shape->GetDY()-fgkSSDConnectorPosition[1],
			  ssdstiffenerboxshapes[0]->GetDZ()+fgkSSDConnectorHeight
		   +  0.5*fgkSSDCapacitor1812Height);
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
				   ssdstiffenerboxshapes[0]->GetDZ()+fgkSSDConnectorHeight
				 + fgkSSDWireRadius,
                   new TGeoRotation("HybridWireRot1",0.,90.,0.));
  hybridwirecombitrans[1] = new TGeoCombiTrans("HybridWireCombiTrans2",
				   0.,
				 - 0.5*fgkSSDConnectorWidth+fgkSSDWireRadius,
				   0.,	
                   new TGeoRotation("HybridWireRot2",
				 -                  wireangle*TMath::RadToDeg(),0.,0.));
  TGeoHMatrix* hybridwirematrix = new TGeoHMatrix();
  hybridwirematrix->MultiplyLeft(hybridwirecombitrans[0]);
  hybridwirematrix->MultiplyLeft(hybridwirecombitrans[1]);
////////////////////////////
// Stiffener Assembly
///////////////////////////
  TGeoVolume* ssdstiffenerassembly = 
								new TGeoVolumeAssembly("SSDStiffenerAssembly");
  ssdstiffenerassembly->AddNode(hybridwire,1,hybridwirematrix);
  for(Int_t i=0; i<kssdstiffenernumber; i++) {
	ssdstiffenerassembly->AddNode(ssdstiffener,i+1,ssdstiffenercombitrans[i]);
	for(Int_t j=1; j<knapacitor0603number+1; j++){
    ssdstiffenerassembly->AddNode(capacitor0603,knapacitor0603number*i+j,new TGeoTranslation("",(j-3.
	)/6*fgkSSDStiffenerLength,
					i*ssdstiffenerseparation+
					0.5*((i==0? 1:-1)*fgkSSDStiffenerWidth
					+(i==0? -1:+1)*fgkSSDCapacitor0603Width),
					-0.5*(fgkSSDStiffenerHeight+fgkSSDCapacitor0603Height)));
	}
	if(i==1) ssdstiffenerassembly->AddNode(capacitor1812,1,capacitor1812trans);
}
  return ssdstiffenerassembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDChipCables(Double_t SSDChipCablesHeigth, 
																	char* side){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Chip Cables    
  /////////////////////////////////////////////////////////////
  const Int_t kssdchipcableslaynumber = 2;
  Int_t ssdchipcablescolor[2] = {fColorAl,fColorPolyhamide};
  Double_t ssdchipcablesradius[2];
  ssdchipcablesradius[0] = 0.25*(SSDChipCablesHeigth
						 -  fgkSSDChipCablesHeight[0]
						 -  fgkSSDChipCablesHeight[1]);
  ssdchipcablesradius[1] = ssdchipcablesradius[0]-fgkSSDChipCablesHeight[0];
  Double_t ssdchipcablespiecelength = 0.5*(fgkSSDChipCablesWidth[0]
								    - 2.*TMath::Pi()*ssdchipcablesradius[0]
									- ssdchipcablesradius[0]
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
  char* ssdchipcablesboxshapename[2*kssdchipcableslaynumber] = 
					  {"SSDChipCablesBoxLay0Shape0","SSDChipCablesBoxLay0Shape1",
					   "SSDChipCablesBoxLay1Shape0","SSDChipCablesBoxLay1Shape1"};
  char* ssdchipcablestubesegshapename[2*kssdchipcableslaynumber] = 
			  {"SSDChipCablesTubeSegLay0Shape0","SSDChipCablesTubeSegLay0Shape1",
			   "SSDChipCablesTubeSegLay1Shape0","SSDChipCablesTubeSegLay1Shape1"};
  TGeoBBox** ssdchipcablesboxshape[kssdchipcableslaynumber];
  TGeoTubeSeg** ssdchipcablestubesegshape[kssdchipcableslaynumber];
  for(Int_t i=0; i<kssdchipcableslaynumber; i++){
    ssdchipcablesboxshape[i]        = new TGeoBBox*[2];
    ssdchipcablestubesegshape[i]    = new TGeoTubeSeg*[2+(side=="Right" ? 0 : 1)];
    ssdchipcablesboxshape[i][0]     = new TGeoBBox(ssdchipcablesboxshapename[2*i],
												   0.5*ssdchipcablespiecelength,
												   0.5*fgkSSDChipCablesLength[1],
												   0.5*fgkSSDChipCablesHeight[i]);
    ssdchipcablesboxshape[i][1]     = new TGeoBBox(ssdchipcablesboxshapename[2*i+1],
						   0.5*(ssdchipcablespiecelength+ssdchipcablesradius[0]
					 + (side=="Right" ? 0. : fgkSSDModuleStiffenerPosition[1])),		
				    0.5*fgkSSDChipCablesLength[1],0.5*fgkSSDChipCablesHeight[i]);
    ssdchipcablestubesegshape[i][0] = 
					   new TGeoTubeSeg(ssdchipcablestubesegshapename[2*i],						
						   ssdchipcablesradius[1]-i*fgkSSDChipCablesHeight[1],
						   ssdchipcablesradius[i],0.5*fgkSSDChipCablesLength[1],
						   0.,180.);
    ssdchipcablestubesegshape[i][1] = 
					   new TGeoTubeSeg(ssdchipcablestubesegshapename[2*i+1],
					         ssdchipcablesradius[0]+i*fgkSSDChipCablesHeight[0],
							 ssdchipcablesradius[0]+fgkSSDChipCablesHeight[0]
					 +		 i*fgkSSDChipCablesHeight[1],
							 0.5*fgkSSDChipCablesLength[1],0.,180.);
    if(side!="Right") ssdchipcablestubesegshape[i][2] = 
					 new TGeoTubeSeg(ssdchipcablestubesegshapename[2*i],
						 0.5*fgkSSDSensorHeight+(1-i)*fgkSSDChipCablesHeight[1],
						 0.5*fgkSSDSensorHeight+(1-i)*fgkSSDChipCablesHeight[0]
					 +   fgkSSDChipCablesHeight[1],
						 0.5*fgkSSDChipCablesLength[1],0.,180.);
  }
  //////////////////////////
  //Box under Chip
  //////////////////////////
  char ssdunderchipcablesboxshapename[30];
  char ssdunderchipcablesboxname[30];	
  char ssdunderchipcablesboxtransname[30];	
  TGeoBBox* ssdunderchipcablesboxshape[kssdchipcableslaynumber];
  TGeoVolume* ssdunderchipcablesbox[kssdchipcableslaynumber]; 
  TGeoTranslation* ssdunderchipcablesboxtrans[kssdchipcableslaynumber];
  for(Int_t i=0; i<kssdchipcableslaynumber; i++){ 
		sprintf(ssdunderchipcablesboxshapename,"SSDUnderChipCablesBoxShape%i",i+1);
		sprintf(ssdunderchipcablesboxname,"SSDUnderChipCablesBox%i",i+1);
		sprintf(ssdunderchipcablesboxtransname,"SSDUnderChipCablesBoxTrans%i",i+1);
		ssdunderchipcablesboxshape[i] = 
								   new TGeoBBox(ssdunderchipcablesboxshapename,
								   0.5*fgkSSDChipWidth,
								   0.5*fgkSSDChipCablesLength[1],
								   0.5*fgkSSDChipCablesHeight[i]);
		ssdunderchipcablesbox[i] = new TGeoVolume(ssdunderchipcablesboxname,
									ssdunderchipcablesboxshape[i],
									(i==0?fSSDAlTraceChipCableMedium:fSSDKaptonChipCableMedium));		
        ssdunderchipcablesbox[i]->SetLineColor(ssdchipcablescolor[i]);
		ssdunderchipcablesboxtrans[i] = 
						new TGeoTranslation(ssdunderchipcablesboxtransname,
						(side=="Right"?-1.:1.)*0.5*fgkSSDChipWidth,
						0.5*(fgkSSDChipCablesLength[0]-fgkSSDChipCablesLength[1])
						+0.5*fgkSSDChipCablesLength[1],
						(i==0?1.:-1.)*0.5*fgkSSDChipCablesHeight[1-i]);
  }
  //////////////////
  //Trapezoid Shapes
  //////////////////
  const Int_t kssdchipcablesvertexnumber = 2;
  const Int_t kssdchipcablestrapezoidnumber = 2;
  TVector3** ssdchipcablestrapezoidvertex[kssdchipcablesvertexnumber];
  for(Int_t i = 0; i< kssdchipcablestrapezoidnumber; i++) 
	 ssdchipcablestrapezoidvertex[i] = new TVector3*[kssdchipcablesvertexnumber];
  //First Shape Vertex Positioning
  ssdchipcablestrapezoidvertex[0][0] = new TVector3();
  ssdchipcablestrapezoidvertex[0][1] = 
		new TVector3(0.5*(fgkSSDChipCablesLength[0]-fgkSSDChipCablesLength[1]));
  //Second Shape Vertex Positioning
  ssdchipcablestrapezoidvertex[1][0] = 
							  new TVector3(*ssdchipcablestrapezoidvertex[0][0]);
  ssdchipcablestrapezoidvertex[1][1] = 
							  new TVector3(*ssdchipcablestrapezoidvertex[0][1]);
  //Setting the names of shapes and volumes
  char* ssdchipcablestrapezoidboxshapename[kssdchipcablestrapezoidnumber] = 
		  {"SSDChipCablesTrapezoidBoxShape1","SSDChipCablesTrapezoidBoxShape2"};
  char* ssdchipcablestrapezoidshapename[kssdchipcablestrapezoidnumber] = 
		  {"SSDChipCablesTrapezoidShape1","SSDChipCablesTrapezoidShape2"};
  char* ssdchipcablestrapezoidboxname[kssdchipcablestrapezoidnumber] = 
		  {"SSDChipCablesTrapezoidBox1","SSDChipCablesTrapezoidBox2"};
  char* ssdhipcablestrapezoidname[kssdchipcablestrapezoidnumber] = 
		  {"SSDChipCablesTrapezoid1","SSDChipCablesTrapezoid2"};
  char* ssdchipcablestrapezoidassemblyname[kssdchipcablestrapezoidnumber] = 
		  {"SSDChipCablesTrapezoidAssembly1","SSDChipCablesTrapezoidAssembly2"};
  //Setting the Shapes
  TGeoBBox* ssdchipcablestrapezoidboxshape[kssdchipcablestrapezoidnumber]; 
  TGeoArb8* ssdchipcablestrapezoidshape[kssdchipcablestrapezoidnumber];
  //Setting the Volumes
  TGeoVolume* ssdchipcablestrapezoidbox[kssdchipcablestrapezoidnumber];
  TGeoVolume* ssdchipcablestrapezoid[kssdchipcablestrapezoidnumber];
  TGeoVolume* ssdchipcablestrapezoidassembly[kssdchipcablestrapezoidnumber]; 
  Double_t ssdchipcablestrapezoidwidth[kssdchipcablesvertexnumber] = 
   {fgkSSDChipCablesWidth[1]+fgkSSDChipCablesWidth[2],fgkSSDChipCablesWidth[1]};
  for(Int_t i=0; i<kssdchipcablestrapezoidnumber; i++){
    ssdchipcablestrapezoidboxshape[i] = 
					new TGeoBBox(ssdchipcablestrapezoidboxshapename[i],
						0.5*fgkSSDChipCablesLength[1],
					    0.5*(fgkSSDChipCablesWidth[1]+fgkSSDChipCablesWidth[2]),
						0.5*fgkSSDChipCablesHeight[i]);
    ssdchipcablestrapezoidshape[i] = 
							  GetTrapezoidShape(ssdchipcablestrapezoidvertex[i],
							  ssdchipcablestrapezoidwidth,
							  fgkSSDChipCablesHeight[i],
							  ssdchipcablestrapezoidshapename[i]);
    ssdchipcablestrapezoidbox[i] = 
						new TGeoVolume(ssdchipcablestrapezoidboxname[i],
									   ssdchipcablestrapezoidboxshape[i],
									   (i==0?fSSDAlTraceChipCableMedium:fSSDKaptonChipCableMedium));
    ssdchipcablestrapezoid[i] = new TGeoVolume(ssdhipcablestrapezoidname[i],
											   ssdchipcablestrapezoidshape[i],
											   (i==0?fSSDAlTraceChipCableMedium:fSSDKaptonChipCableMedium));
    ssdchipcablestrapezoidbox[i]->SetLineColor(ssdchipcablescolor[i]);
    ssdchipcablestrapezoid[i]->SetLineColor(ssdchipcablescolor[i]);
    ssdchipcablestrapezoidassembly[i] = 
				new TGeoVolumeAssembly(ssdchipcablestrapezoidassemblyname[i]);
    ssdchipcablestrapezoidassembly[i]->AddNode(ssdchipcablestrapezoidbox[i],1,
				new TGeoTranslation(0.5*fgkSSDChipCablesLength[1],
				   0.5*(fgkSSDChipCablesWidth[1]+fgkSSDChipCablesWidth[2]),0.));
    ssdchipcablestrapezoidassembly[i]->AddNode(ssdchipcablestrapezoid[i],0,
			new TGeoCombiTrans(0.,0.,0.,new TGeoRotation("",90.,180.,-90.)));
    ssdchipcablestrapezoidassembly[i]->AddNode(ssdchipcablestrapezoid[i],1,
			new TGeoTranslation(fgkSSDChipCablesLength[1],0.,0.));
  }
  /////////////////////////////
  //Box and Tube Seg CombiTrans
  /////////////////////////////
  TGeoTranslation* ssdchipcablesboxtrans[2*kssdchipcableslaynumber];
  ssdchipcablesboxtrans[0] = 
					new TGeoTranslation("SSDChipCablesLay1Box1Trans",0.,0.,0.);
  ssdchipcablesboxtrans[1] = 
					new TGeoTranslation("SSDChipCablesLay1Box2Trans",
										 ssdchipcablesboxshape[0][1]->GetDX()
						   -             0.5*ssdchipcablespiecelength,
                       0.0,
						   -             2.*ssdchipcablesradius[0]
						   -             fgkSSDChipCablesHeight[0]);
  ssdchipcablesboxtrans[2] = new TGeoTranslation("SSDChipCablesLay2Box1Trans",
										 0.0,
										 0.0,
										 0.5*(fgkSSDChipCablesHeight[0]
						   +			 fgkSSDChipCablesHeight[1]));
  ssdchipcablesboxtrans[3] = 
							 new TGeoTranslation("SSDChipCablesLay2Box2Trans",
										 ssdchipcablesboxshape[1][1]->GetDX()
						   -			 0.5*ssdchipcablespiecelength,
										 0.0,
						   -			 2.*ssdchipcablesradius[0]
						   -			 0.5*fgkSSDChipCablesHeight[1]
						   -			 1.5*fgkSSDChipCablesHeight[0]);
  TGeoRotation* ssdchipcablesrot[3];
  ssdchipcablesrot[0] = new TGeoRotation("SSDChipCablesRot1",0.,90.,0.);
  ssdchipcablesrot[1] = new TGeoRotation("SSDChipCablesRot2",90.,90.,-90.);
  ssdchipcablesrot[2] = new TGeoRotation("SSDChipCablesRot3",90.,-90.,-90.);
  TGeoCombiTrans* ssdchipcablestubesegcombitrans[2*(kssdchipcableslaynumber+1)];    
//  TGeoCombiTrans* SSDChipCablesTubeSegCombiTrans[2*(SSDChipCablesLayNumber+
//													  (side=="Right" ? 0 : 1))];
  ssdchipcablestubesegcombitrans[0] = 
				new TGeoCombiTrans("SSDChipCablesLay1TubeSeg1CombiTrans",
				0.5*ssdchipcablespiecelength,
				0.0,
				ssdchipcablesradius[0]
			-   0.5*fgkSSDChipCablesHeight[0],
				new TGeoRotation((*ssdchipcablesrot[1])*(*ssdchipcablesrot[0])));
  ssdchipcablestubesegcombitrans[1] = 
				new TGeoCombiTrans("SSDChipCablesLay1TubeSeg2CombiTrans",
			-   0.5*ssdchipcablespiecelength,
				0.0,
			-   ssdchipcablesradius[0]-0.5*fgkSSDChipCablesHeight[0],
				new TGeoRotation((*ssdchipcablesrot[2])*(*ssdchipcablesrot[0])));
  ssdchipcablestubesegcombitrans[2] = 
  new TGeoCombiTrans("SSDChipCablesLay2TubeSeg1CombiTrans",
				0.5*ssdchipcablespiecelength,
				0.0,
				ssdchipcablesradius[0]-0.5*fgkSSDChipCablesHeight[0],
				new TGeoRotation((*ssdchipcablesrot[1])*(*ssdchipcablesrot[0])));
  ssdchipcablestubesegcombitrans[3] = 
				new TGeoCombiTrans("SSDChipCablesLay2TubeSeg2CombiTrans",
			-	0.5*ssdchipcablespiecelength,
				0.0,
			-	ssdchipcablesradius[0]+0.5*fgkSSDChipCablesHeight[0]
			-   fgkSSDChipCablesHeight[0],
				new TGeoRotation((*ssdchipcablesrot[2])*(*ssdchipcablesrot[0])));
  ssdchipcablestubesegcombitrans[4] = 
				new TGeoCombiTrans("SSDChipCablesLay1TubeSeg4CombiTrans",
				0.5*ssdchipcablespiecelength+ssdchipcablesradius[0]
			+   fgkSSDModuleStiffenerPosition[1],
				0.0,
			-	2.*ssdchipcablesradius[0]-0.5*fgkSSDChipCablesHeight[0]
			-   (0.5*fgkSSDSensorHeight+fgkSSDChipCablesHeight[0]
			+	fgkSSDChipCablesHeight[1]),
			new TGeoRotation((*ssdchipcablesrot[1])*(*ssdchipcablesrot[0])));
  ssdchipcablestubesegcombitrans[5] = 
			new TGeoCombiTrans("SSDChipCablesLay2TubeSeg5CombiTrans",
				0.5*ssdchipcablespiecelength+ssdchipcablesradius[0]
			+	fgkSSDModuleStiffenerPosition[1],
				0.0,
			-	2.*ssdchipcablesradius[0]-0.5*fgkSSDChipCablesHeight[0]
			-	(0.5*fgkSSDSensorHeight+fgkSSDChipCablesHeight[0]
			+	fgkSSDChipCablesHeight[1]),
			new TGeoRotation((*ssdchipcablesrot[1])*(*ssdchipcablesrot[0])));
  TGeoCombiTrans* ssdchipcablestrapezoidcombitrans[kssdchipcableslaynumber];
  ssdchipcablestrapezoidcombitrans[0] = (side=="Right" ? 
			new TGeoCombiTrans("SSDChipCableLay1TrapezoidRightCombiTrans",
				0.5*ssdchipcablespiecelength+ssdchipcablestrapezoidwidth[0]
			+	ssdchipcablesradius[0],
			-	0.5*fgkSSDChipCablesLength[1],
			-	fgkSSDChipCablesHeight[0]-2.*ssdchipcablesradius[0],
			new TGeoRotation("",90.,0.,0.)) :
			new TGeoCombiTrans("SSDChipCableLay1TrapezoidLeftCombiTrans",
			-	2.*(fgkSSDChipCablesWidth[1]+fgkSSDChipCablesWidth[2])
			+	0.5*ssdchipcablespiecelength+ssdchipcablestrapezoidwidth[0]
			+	ssdchipcablesradius[0]+fgkSSDModuleStiffenerPosition[1],
				0.5*fgkSSDChipCablesLength[1],
			-	2.*(fgkSSDChipCablesHeight[0]+fgkSSDChipCablesHeight[1])
			-	2.*ssdchipcablesradius[0]-fgkSSDSensorHeight,
			new TGeoRotation("",-90.,0.,0.)));
  ssdchipcablestrapezoidcombitrans[1] = (side=="Right" ? 
			new TGeoCombiTrans("SSDChipCableLay2TrapezoidRightCombiTrans",
				0.5*ssdchipcablespiecelength+ssdchipcablestrapezoidwidth[0]
			+	ssdchipcablesradius[0],
			-	0.5*fgkSSDChipCablesLength[1],
			-	0.5*(fgkSSDChipCablesHeight[0]+fgkSSDChipCablesHeight[1])
			-	fgkSSDChipCablesHeight[0]-2.*ssdchipcablesradius[0],
				new TGeoRotation("",90.,0.,0.)) :
				new TGeoCombiTrans("SSDChipCableLay2TrapezoidLeftCombiTrans",
			-	2.*(fgkSSDChipCablesWidth[1]+fgkSSDChipCablesWidth[2])
			+	0.5*ssdchipcablespiecelength+ssdchipcablestrapezoidwidth[0]
			+	ssdchipcablesradius[0]+fgkSSDModuleStiffenerPosition[1],
				0.5*fgkSSDChipCablesLength[1],-0.5*(fgkSSDChipCablesHeight[0]
			+	fgkSSDChipCablesHeight[1])-fgkSSDChipCablesHeight[1]
			-	fgkSSDChipCablesHeight[0]-2.*ssdchipcablesradius[0]
			-	fgkSSDSensorHeight,new TGeoRotation("",-90.,0.,0.)));  
  //////////////////////////
  //Box and Tube Seg Volumes
  //////////////////////////
  char* ssdchipcablesboxname[2*kssdchipcableslaynumber] = 
							 {"SSDChipCablesLay1Box1","SSDChipCablesLay1Box2",
							  "SSDChipCablesLay2Box1","SSDChipCablesLay2Box2"};
  char* ssdchiprightcablestubesegname[2*kssdchipcableslaynumber] = 
			  {"SSDChipRightCablesLay1TubeSeg1","SSDChipRightCablesLay1TubeSeg2",
			   "SSDChipRightCablesLay2TubeSeg1","SSDChipRightCablesLay2TubeSeg2"};
  char* ssdchipLeftcablestubesegname[2*kssdchipcableslaynumber] = 
			  {"SSDChipLeftCablesLay1TubeSeg1","SSDChipLeftCablesLay1TubeSeg2",
			   "SSDChipLeftCablesLay2TubeSeg1","SSDChipLeftCablesLay2TubeSeg2"};
  char* ssdchipcableslayassemblyname[kssdchipcableslaynumber] = 
			  {"SSDChipCablesLay1","SSDChipCablesLay2"};
  TGeoVolume** ssdchipcablesbox[kssdchipcableslaynumber];
  TGeoVolume** ssdchipcablestubeseg[kssdchipcableslaynumber];
  TGeoVolume* ssdchipcableslayassembly[kssdchipcableslaynumber];
  for(Int_t i=0; i<kssdchipcableslaynumber; i++){
    TGeoMedium* ssdchipcableslaymed = 
            (i==0?fSSDAlTraceChipCableMedium:fSSDKaptonChipCableMedium);
    ssdchipcablesbox[i] = new TGeoVolume*[2];
    ssdchipcablestubeseg[i] = new TGeoVolume*[2+(side=="Right" ? 0 : 1)];
    ssdchipcablesbox[i][0] = new TGeoVolume(ssdchipcablesboxname[2*i],
					   ssdchipcablesboxshape[i][0],ssdchipcableslaymed);
    ssdchipcablesbox[i][1] = new TGeoVolume(ssdchipcablesboxname[2*i+1],
					   ssdchipcablesboxshape[i][1],ssdchipcableslaymed);
    ssdchipcablestubeseg[i][0] = new TGeoVolume(ssdchiprightcablestubesegname[2*i],
					   ssdchipcablestubesegshape[i][0],ssdchipcableslaymed);
    ssdchipcablestubeseg[i][1] = new TGeoVolume(ssdchiprightcablestubesegname[2*i+1],
					   ssdchipcablestubesegshape[i][1],ssdchipcableslaymed);
    ssdchipcablesbox[i][0]->SetLineColor(ssdchipcablescolor[i]);
    ssdchipcablesbox[i][1]->SetLineColor(ssdchipcablescolor[i]);
    ssdchipcablestubeseg[i][0]->SetLineColor(ssdchipcablescolor[i]);
    ssdchipcablestubeseg[i][1]->SetLineColor(ssdchipcablescolor[i]);
    ssdchipcableslayassembly[i] = new TGeoVolumeAssembly(ssdchipcableslayassemblyname[i]);
    ssdchipcableslayassembly[i]->AddNode(ssdchipcablesbox[i][0],1,
										 ssdchipcablesboxtrans[2*i]);
    ssdchipcableslayassembly[i]->AddNode(ssdchipcablesbox[i][1],1,
										 ssdchipcablesboxtrans[2*i+1]);
    ssdchipcableslayassembly[i]->AddNode(ssdchipcablestubeseg[i][0],1,
										 ssdchipcablestubesegcombitrans[2*i]);
    ssdchipcableslayassembly[i]->AddNode(ssdchipcablestubeseg[i][1],1,
										 ssdchipcablestubesegcombitrans[2*i+1]);
    if(side!="Right"){
      ssdchipcablestubeseg[i][2] = new TGeoVolume(ssdchipLeftcablestubesegname[2*i],
												  ssdchipcablestubesegshape[i][2],
												  ssdchipcableslaymed);
      ssdchipcablestubeseg[i][2]->SetLineColor(ssdchipcablescolor[i]);
      ssdchipcableslayassembly[i]->AddNode(ssdchipcablestubeseg[i][2],1,
										   ssdchipcablestubesegcombitrans[4+i]);
    }
    ssdchipcableslayassembly[i]->AddNode(ssdchipcablestrapezoidassembly[i],1,
										 ssdchipcablestrapezoidcombitrans[i]);
  }
  TGeoCombiTrans* ssdchipcablescombitrans[kssdchipcableslaynumber];
  ssdchipcablescombitrans[0] = new TGeoCombiTrans("SSDChipCablesCombiTrans1",
					   (side=="Right" ? -1 : 1)*0.5*ssdchipcablespiecelength,
						0.5*fgkSSDChipCablesLength[0],
					-	(2.*ssdchipcablesradius[0]-fgkSSDChipCablesHeight[0]
					-	0.5*fgkSSDChipCablesHeight[1]),
						new TGeoRotation("",(side=="Right" ? 0 : 1)*180.,0.,0.));
  ssdchipcablescombitrans[1] = new TGeoCombiTrans("SSDChipCablesCombiTrans2",
						(side=="Right" ? -1 : 1)*0.5*ssdchipcablespiecelength,
						0.5*fgkSSDChipCablesLength[0],
					-	(2.*ssdchipcablesradius[0]-fgkSSDChipCablesHeight[0]
					-	0.5*fgkSSDChipCablesHeight[1]),
						new TGeoRotation("",(side=="Right" ? 0 : 1)*180.,0.,0.));
  TGeoVolume* ssdchipcablesassembly = 
						new TGeoVolumeAssembly("SSDChipCables");
  for(Int_t i=0; i<kssdchipcableslaynumber; i++){ 
		ssdchipcablesassembly->AddNode(ssdchipcableslayassembly[i],1,
													ssdchipcablescombitrans[i]);
		ssdchipcablesassembly->AddNode(ssdunderchipcablesbox[i],1,ssdunderchipcablesboxtrans[i]);
  }
  return ssdchipcablesassembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDChipCablesAssembly(Double_t SSDChipCablesHeigth){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Chip Cables Assembly    
  /////////////////////////////////////////////////////////////
  const Int_t kchipcablesnumber = 2;
  Double_t chipcablestransvector = fgkSSDSensorLength
								 - 2.*fgkSSDModuleStiffenerPosition[1]
								 - 2.*(fgkSSDStiffenerWidth
								 - fgkSSDStiffenerToChipDist-fgkSSDChipWidth);
  char* ssdchipcablesname[kchipcablesnumber] = {"Right","Left"};
  TGeoVolume* ssdchipcables[kchipcablesnumber];  
  TGeoVolume* ssdchipcablesassembly = 
					 new TGeoVolumeAssembly("SSDChipCablesAssembly");
  for(Int_t i=0; i<kchipcablesnumber; i++) ssdchipcables[i] = 
					 GetSSDChipCables(SSDChipCablesHeigth,ssdchipcablesname[i]);
  for(Int_t i=0; i<kchipcablesnumber; i++)
    for(Int_t j=0; j<fgkSSDChipNumber; j++)
      ssdchipcablesassembly->AddNode(ssdchipcables[i],fgkSSDChipNumber*i+j+1,
			new TGeoTranslation(-(ssdchipcablesname[i]=="Left" ? 1. : 0.)
		*	chipcablestransvector,(j-0.5)*fgkSSDChipCablesLength[0]
		+	0.5*fgkSSDChipCablesLength[1],0.));
  return ssdchipcablesassembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDFlex(Double_t ssdflexradius, Double_t SSDFlexHeigth){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Flex    
  /////////////////////////////////////////////////////////////
  const Int_t kssdflexvolumenumber = 3;
  TGeoVolume* ssdflexvolume[kssdflexvolumenumber];
  ////////////////////////
  // Setting Display Color
  ////////////////////////
  Int_t ssdflexcolor;
  ssdflexcolor = (SSDFlexHeigth == fgkSSDFlexHeight[0] ? fColorAl: fColorPolyhamide);
  TGeoMedium* ssdflexmed = (SSDFlexHeigth == fgkSSDFlexHeight[0] ? fSSDAlTraceFlexMedium :
                            fSSDKaptonFlexMedium);
  ////////////////////////
  // SSDFlexTrapezoidShape
  ////////////////////////
  const Int_t kssdflexvertexnumber = 2;
  Double_t ssdflexwidth[kssdflexvertexnumber] = {fgkSSDFlexWidth[1],
												fgkSSDFlexWidth[0]};
  TVector3* ssdflexvertex[kssdflexvertexnumber];
  ssdflexvertex[0] = new TVector3();
  ssdflexvertex[1] = new TVector3(fgkSSDFlexLength[0]-fgkSSDFlexLength[1]);
  TGeoArb8* ssdflextrapezoidshape = GetTrapezoidShape(ssdflexvertex,
													  ssdflexwidth,SSDFlexHeigth,
													  "SSDFlexTrapezoidShape");
  ssdflexvolume[0] = new TGeoVolume("SSDFlexTrapezoid",ssdflextrapezoidshape,ssdflexmed);
  ssdflexvolume[0]->SetLineColor(ssdflexcolor);
  /////////////////////////
  //SSDFlexTubeSeg Assembly
  /////////////////////////
  const Int_t kssdflextubesegnumber = 2;
  TGeoTubeSeg* ssdflextubesegshape[kssdflextubesegnumber];
  Double_t ssdflexradiusmax = (fgkSSDFlexLength[3]-fgkSSDFlexLength[2])
							/  TMath::Tan(fgkSSDFlexAngle*TMath::DegToRad());
  ssdflextubesegshape[0] = new TGeoTubeSeg("SSDFlexTubeSegShape1",
						   ssdflexradius,ssdflexradius+SSDFlexHeigth,
						   0.5*fgkSSDFlexWidth[0],0.,180.);
  ssdflextubesegshape[1] = new TGeoTubeSeg("SSDFlexTubeSegShape2",
						   ssdflexradiusmax-ssdflexradius-SSDFlexHeigth,
						   ssdflexradiusmax-ssdflexradius,0.5*fgkSSDFlexWidth[0],
						   0.,2.*fgkSSDFlexAngle);
  TGeoRotation** ssdflextubsegrot[kssdflextubesegnumber];
  for(Int_t i = 0; i<kssdflextubesegnumber; i++) 
									  ssdflextubsegrot[i] = new TGeoRotation*[2]; 
  ssdflextubsegrot[0][0] = new TGeoRotation("SSDFlexTubeSeg1Rot1", 0., 90.,  0.);
  ssdflextubsegrot[0][1] = new TGeoRotation("SSDFlexTubeSeg1Rot2",90., 90.,-90.);
  ssdflextubsegrot[1][0] = new TGeoRotation("SSDFlexTubeSeg2Rot1", 0.,-90.,  0.);
  ssdflextubsegrot[1][1] = new TGeoRotation("SSDFlexTubeSeg2Rot2",90., 90.,-90.);
  TGeoCombiTrans* ssdflextubesegcombitrans[kssdflextubesegnumber];
  ssdflextubesegcombitrans[0] = new TGeoCombiTrans("SSDFlexTubeSegCombiTrans1",
								fgkSSDFlexLength[0],0.5*fgkSSDFlexWidth[0],
								ssdflexradius+0.5*SSDFlexHeigth,
								new TGeoRotation((*ssdflextubsegrot[0][1])
							*	(*ssdflextubsegrot[0][0])));
  ssdflextubesegcombitrans[1] = new TGeoCombiTrans("SSDFlexTubeSegCombiTrans2",
								fgkSSDFlexLength[0]-fgkSSDFlexLength[2],
								0.5*fgkSSDFlexWidth[0],
								ssdflexradiusmax+0.5*SSDFlexHeigth+ssdflexradius,
								new TGeoRotation((*ssdflextubsegrot[1][1])
							*	(*ssdflextubsegrot[1][0])));
  ssdflexvolume[1] = new TGeoVolumeAssembly("SSDFlexTubeSegAssembly");
  TGeoVolume* ssdflextubeseg[kssdflextubesegnumber];
  char ssdflextubesegname[30];
  for(Int_t i=0; i<kssdflextubesegnumber; i++){ 
		sprintf(ssdflextubesegname,"SSDFlexTubeSeg%i",i+1);
		ssdflextubeseg[i] = new TGeoVolume(ssdflextubesegname,ssdflextubesegshape[i],
                                     ssdflexmed);
		ssdflextubeseg[i]->SetLineColor(ssdflexcolor);
        ssdflexvolume[1]->AddNode(ssdflextubeseg[i],1,ssdflextubesegcombitrans[i]);
  }
  ///////////
  //Box Shape 
  ///////////
  const Int_t kssdflexboxnumber = 7;
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
  ssdflexboxlength[5] = fgkSSDFlexLength[2];	
  ssdflexboxlength[6] = fgkSSDFlexFullLength-2.*fgkSSDFlexAngle
					  *	TMath::DegToRad()*ssdflexradiusmax
					  - fgkSSDFlexLength[2]-TMath::Pi()
					  *	fgkSSDStiffenerHeight-fgkSSDFlexLength[0];	
  Double_t ssdflexboxwidth[kssdflexboxnumber];
  ssdflexboxwidth[0] = fgkSSDFlexWidth[0];
  ssdflexboxwidth[1] = fgkSSDFlexWidth[0]-fgkSSDFlexHoleWidth;
  ssdflexboxwidth[2] = fgkSSDFlexHoleWidth;
  ssdflexboxwidth[3] = ssdflexboxwidth[2]-fgkSSDFlexHoleLength;
  ssdflexboxwidth[4] = fgkSSDFlexWidth[0];
  ssdflexboxwidth[5] = fgkSSDFlexWidth[0];
  ssdflexboxwidth[6] = fgkSSDFlexWidth[0];
  TGeoBBox* ssdflexboxshape[kssdflexboxnumber+1];
  for(Int_t i=0; i<kssdflexboxnumber+1; i++) ssdflexboxshape[i] = 
						(i!= kssdflexboxnumber ? new TGeoBBox("SSDFlexBoxShape",
								0.5*ssdflexboxlength[i],
								0.5*ssdflexboxwidth[i],0.5*SSDFlexHeigth) : 
								ssdflexboxshape[2]);
  //////////////////////////////
  //SSDFlex Box Shape CombiTrans 
  //////////////////////////////
  TGeoCombiTrans* ssdflexboxcombitrans[kssdflexboxnumber+1];
  ssdflexboxcombitrans[0] = new TGeoCombiTrans("SSDFlexBoxCombiTrans0",
								ssdflexvertex[1]->X()+0.5*ssdflexboxlength[0],
								0.5*fgkSSDFlexWidth[0],0.,0);
  ssdflexboxcombitrans[1] = new TGeoCombiTrans("SSDFlexBoxCombiTrans1",
								ssdflexvertex[1]->X()+ssdflexboxlength[0]
							+	0.5*ssdflexboxlength[1],
								fgkSSDFlexHoleWidth+0.5*ssdflexboxwidth[1],0.,0);
  ssdflexboxcombitrans[2] = new TGeoCombiTrans("SSDFlexBoxCombiTrans2",
								ssdflexvertex[1]->X()+ssdflexboxlength[0]
							+	fgkSSDFlexHoleLength+0.5*ssdflexboxlength[2],
								0.5*ssdflexboxwidth[2],0.,0);
  ssdflexboxcombitrans[3] = new TGeoCombiTrans("SSDFlexBoxCombiTrans3",
								ssdflexvertex[1]->X()+ssdflexboxlength[0]
							+	fgkSSDFlexHoleLength+ssdflexboxlength[2]
							+	0.5*fgkSSDFlexHoleWidth,
								fgkSSDFlexHoleLength+0.5*ssdflexboxwidth[3],0.,0);
  ssdflexboxcombitrans[4] = new TGeoCombiTrans("SSDFlexBoxCombiTrans4",
								ssdflexvertex[1]->X()+ssdflexboxlength[0]
							+	ssdflexboxlength[1]+0.5*ssdflexboxlength[4],
								0.5*fgkSSDFlexWidth[0],0.,0);
  ssdflexboxcombitrans[5] = new TGeoCombiTrans("SSDFlexBoxCombiTrans5",
							-	0.5*fgkSSDFlexLength[2]+fgkSSDFlexLength[0],
								0.5*fgkSSDFlexWidth[0],
								2.*ssdflexradius+SSDFlexHeigth,0);
  ssdflexboxcombitrans[6] = new TGeoCombiTrans("SSDFlexBoxCombiTrans6",
							-	ssdflexboxshape[6]->GetDX()
							+	ssdflexboxshape[6]->GetDX()
							*	TMath::Cos(2.*fgkSSDFlexAngle*TMath::DegToRad())
							+	fgkSSDFlexLength[0]-fgkSSDFlexLength[2]
							-	(ssdflexradiusmax-ssdflexradius-0.5*SSDFlexHeigth)
							*	TMath::Cos(fgkSSDFlexAngle*TMath::DegToRad()),
								0.5*fgkSSDFlexWidth[0],ssdflexboxshape[6]->GetDX()
								*TMath::Sin(2.*fgkSSDFlexAngle*TMath::DegToRad())
							+	SSDFlexHeigth+2.*ssdflexradius+(ssdflexradiusmax
							-	ssdflexradius-0.5*SSDFlexHeigth)
							*	TMath::Sin(fgkSSDFlexAngle*TMath::DegToRad()),
								new TGeoRotation("",90.,2.*fgkSSDFlexAngle,-90.));
  ssdflexboxcombitrans[7] = new TGeoCombiTrans("SSDFlexBoxCombiTrans7",
								ssdflexvertex[1]->X()+ssdflexboxlength[0]
							+	fgkSSDFlexHoleLength+1.5*ssdflexboxlength[2]
							+	ssdflexboxlength[3],
								0.5*ssdflexboxwidth[2],0.,0);
  ////////////////////////////
  //SSDFlex Box Shape Assembly 
  ////////////////////////////
  ssdflexvolume[2] = new TGeoVolumeAssembly("SSDFlexBoxAssembly");
  TGeoVolume* ssdflexbox[kssdflexboxnumber+1];
  TGeoVolume* ssdendflex = GetSSDEndFlex(ssdflexboxlength[6],SSDFlexHeigth);
  TGeoHMatrix* ssdendflexhmatrix = new TGeoHMatrix();
  TGeoRotation* ssdendflexrot = new TGeoRotation("SSDEndFlexRot",180.,0.,0);
  ssdendflexhmatrix->MultiplyLeft(ssdendflexrot);
  ssdendflexhmatrix->MultiplyLeft(ssdflexboxcombitrans[6]);
  char ssdflexboxname[30];
  for(Int_t i=0; i<kssdflexboxnumber+1; i++){
	sprintf(ssdflexboxname,"SSDFlexBox%i",i!=kssdflexboxnumber?i+1:7);
	if(i==6){ssdflexvolume[2]->AddNode(ssdendflex,1,ssdendflexhmatrix);}
	else{
    ssdflexbox[i] = new TGeoVolume(ssdflexboxname,ssdflexboxshape[i],
                                   ssdflexmed);
	ssdflexbox[i]->SetLineColor(ssdflexcolor);
	ssdflexvolume[2]->AddNode(ssdflexbox[i],1,ssdflexboxcombitrans[i]);}
 }
  //////////////////////
  //SSDFlex Construction
  //////////////////////
  TGeoVolume* ssdflex = new TGeoVolumeAssembly("SSDFlex");
  for(Int_t i =0; i<kssdflexvolumenumber; i++) ssdflex->AddNode(ssdflexvolume[i],1);
  return ssdflex;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDEndFlex(Double_t SSDEndFlexLength, 
														Double_t SSDFlexHeigth){
  /////////////////////////////////////////////////////////////
  // Method generating SSD End Flex   
  /////////////////////////////////////////
  // Setting Display Color, Media and Index
  /////////////////////////////////////////
  Int_t ssdflexcolor;
  ssdflexcolor = (SSDFlexHeigth == fgkSSDFlexHeight[0] ? fColorAl: fColorPolyhamide);
  TGeoMedium* ssdflexmed = (SSDFlexHeigth == fgkSSDFlexHeight[0] ? fSSDAlTraceFlexMedium :
                            fSSDKaptonFlexMedium);
  ////////////////////////
  const Int_t kssdendflexboxnumber = 5;
  TGeoBBox* ssdendflexbboxshape[kssdendflexboxnumber];
  ssdendflexbboxshape[0] = new TGeoBBox("SSDFlexBoxShape1",
								   0.5*SSDEndFlexLength,0.5*fgkSSDFlexWidth[0],
								   0.5*SSDFlexHeigth);
  ssdendflexbboxshape[1] = new TGeoBBox("SSDFlexBoxShape2",
				    0.5*fgkSSDEndFlexCompLength[1],
					0.5*(fgkSSDEndFlexCompWidth[0]-fgkSSDFlexWidth[0])/2,
					0.5*SSDFlexHeigth);
  ssdendflexbboxshape[2] = new TGeoBBox("SSDFlexBoxShape3",
				    0.5*fgkSSDEndFlexCompLength[2],
					0.5*(fgkSSDEndFlexCompWidth[1]-fgkSSDFlexWidth[0])/2,
					0.5*SSDFlexHeigth);
  ssdendflexbboxshape[3] = new TGeoBBox("SSDFlexBoxShape4",
				    0.5*fgkSSDEndFlexCompLength[3],
					0.5*(fgkSSDEndFlexCompWidth[0]-fgkSSDFlexWidth[0])/2,
					0.5*SSDFlexHeigth);
  ssdendflexbboxshape[4] = new TGeoBBox("SSDFlexBoxShape5",
				    0.5*(fgkSSDEndFlexCompLength[4]+fgkSSDEndFlexCompLength[5]),
					0.25*(fgkSSDEndFlexCompWidth[2]-fgkSSDFlexWidth[0])/2,
					0.5*SSDFlexHeigth);
  TGeoVolume* ssdendflexbbox[kssdendflexboxnumber];  
  char ssdendflexbboxname[30];
  for(Int_t i=0; i<kssdendflexboxnumber; i++){
	sprintf(ssdendflexbboxname,"SSDEndFlexBBox%i",i+1);
	ssdendflexbbox[i] = new TGeoVolume(ssdendflexbboxname,
                     ssdendflexbboxshape[i],
                     ssdflexmed);
	ssdendflexbbox[i]->SetLineColor(ssdflexcolor);
  }
  TGeoVolume* ssdendflex = new TGeoVolumeAssembly("SSDEndFlex");
  Double_t partialsumlength = 0.;
  for(Int_t i=0; i<kssdendflexboxnumber+1; i++) partialsumlength += fgkSSDEndFlexCompLength[i];
  Double_t referencelength = SSDEndFlexLength-partialsumlength;
  ssdendflex->AddNode(ssdendflexbbox[0],1);
  ssdendflex->AddNode(ssdendflexbbox[1],1,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+referencelength+fgkSSDEndFlexCompLength[0]
					+  0.5*fgkSSDEndFlexCompLength[1],
					   0.5*fgkSSDFlexWidth[0]+ssdendflexbboxshape[1]->GetDY(),
					   0.));
  ssdendflex->AddNode(ssdendflexbbox[1],2,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+referencelength+fgkSSDEndFlexCompLength[0]
					+  0.5*fgkSSDEndFlexCompLength[1],
					-  0.5*fgkSSDFlexWidth[0]-ssdendflexbboxshape[1]->GetDY(),
					   0.));
  ssdendflex->AddNode(ssdendflexbbox[2],1,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+referencelength+fgkSSDEndFlexCompLength[0]
					+  fgkSSDEndFlexCompLength[1]+0.5*fgkSSDEndFlexCompLength[2],
					+  0.5*fgkSSDFlexWidth[0]+ssdendflexbboxshape[2]->GetDY(),
					   0.));
  ssdendflex->AddNode(ssdendflexbbox[2],2,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+referencelength+fgkSSDEndFlexCompLength[0]
					+  fgkSSDEndFlexCompLength[1]+0.5*fgkSSDEndFlexCompLength[2],
					-  0.5*fgkSSDFlexWidth[0]-ssdendflexbboxshape[2]->GetDY(),
					   0.));
  ssdendflex->AddNode(ssdendflexbbox[3],1,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+referencelength+fgkSSDEndFlexCompLength[0]
					+  fgkSSDEndFlexCompLength[1]+fgkSSDEndFlexCompLength[2]
					+  0.5*fgkSSDEndFlexCompLength[3],
					+  0.5*fgkSSDFlexWidth[0]+ssdendflexbboxshape[3]->GetDY(),
					   0.));
  ssdendflex->AddNode(ssdendflexbbox[3],2,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+referencelength+fgkSSDEndFlexCompLength[0]
					+  fgkSSDEndFlexCompLength[1]+fgkSSDEndFlexCompLength[2]
					+  0.5*fgkSSDEndFlexCompLength[3],
					-  0.5*fgkSSDFlexWidth[0]-ssdendflexbboxshape[3]->GetDY(),
					   0.));
  ssdendflex->AddNode(ssdendflexbbox[4],1,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+referencelength+fgkSSDEndFlexCompLength[0]
					+  fgkSSDEndFlexCompLength[1]+fgkSSDEndFlexCompLength[2]
					+  fgkSSDEndFlexCompLength[3]+0.5*(fgkSSDEndFlexCompLength[4]+fgkSSDEndFlexCompLength[5]),
					+  0.5*fgkSSDFlexWidth[0]+ssdendflexbboxshape[4]->GetDY(),
					   0.));
  ssdendflex->AddNode(ssdendflexbbox[4],2,new TGeoTranslation(
					-  0.5*SSDEndFlexLength+referencelength+fgkSSDEndFlexCompLength[0]
					+  fgkSSDEndFlexCompLength[1]+fgkSSDEndFlexCompLength[2]
					+  fgkSSDEndFlexCompLength[3]+0.5*(fgkSSDEndFlexCompLength[4]
					+  fgkSSDEndFlexCompLength[5]),
					-  0.5*fgkSSDFlexWidth[0]-ssdendflexbboxshape[4]->GetDY(),
					   0.));
  return ssdendflex;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDFlexAssembly(){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Flex Assembly    
  /////////////////////////////////////////////////////////////
  TGeoVolume* ssdflexassembly = new TGeoVolumeAssembly("SSDFlexAssembly");
  const Int_t kssdflexlayernumber = 4;
  Double_t ssdflexheight[kssdflexlayernumber];
  Double_t ssdflexradius[kssdflexlayernumber];
  TGeoTranslation* ssdflextrans[kssdflexlayernumber];
  for(Int_t i=0; i<kssdflexlayernumber; i++){ 
    ssdflexheight[i] = (i%2==0 ? fgkSSDFlexHeight[0] : fgkSSDFlexHeight[1]);
    ssdflexradius[i] = (i==0 ? fgkSSDStiffenerHeight : ssdflexradius[i-1]
					 +					   ssdflexheight[i-1]);
    ssdflextrans[i]  = new TGeoTranslation(0.,0.,-0.5*i*(ssdflexheight[0]
					 +					   ssdflexheight[1])); 
    ssdflexassembly->AddNode(GetSSDFlex(ssdflexradius[i],ssdflexheight[i]),i+1,
										   ssdflextrans[i]);   
  }
  return ssdflexassembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDCoolingBlockAssembly(){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Cooling Block Assembly    
  /////////////////////////////////////////////////////////////
  const Int_t kssdcoolingblocktransnumber = 2;
  Double_t ssdcoolingblocktransvector[kssdcoolingblocktransnumber] = 
					{fgkSSDModuleSensorSupportDistance+fgkSSDCoolingBlockLength,
					 fgkSSDSensorLength-2.*fgkSSDModuleStiffenerPosition[1]
				-	 fgkSSDCoolingBlockWidth};
  TGeoVolume* ssdcoolingblock = GetSSDCoolingBlock();
  TGeoVolume* ssdcoolingblockassembly = 
							  new TGeoVolumeAssembly("SSDCoolingBlockAssembly");
  for(Int_t i=0; i<kssdcoolingblocktransnumber; i++)
    for(Int_t j=0; j<kssdcoolingblocktransnumber; j++) 
		ssdcoolingblockassembly->AddNode(ssdcoolingblock,
						  kssdcoolingblocktransnumber*i+j+1,
						  new TGeoTranslation(i*ssdcoolingblocktransvector[0],
						  j*ssdcoolingblocktransvector[1],0.));
  return ssdcoolingblockassembly;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetSSDCoolingBlock(){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Cooling Block    
  /////////////////////////////////////////////////////////////
  // Center Cooling Block Hole
  ////////////////////////////
  Double_t coolingblockholeangle = TMath::ACos(0.5*fgkSSDCoolingBlockHoleLength[0]
							/fgkSSDCoolingBlockHoleRadius[0])*TMath::RadToDeg();
  Double_t coolingblockholewidth = fgkSSDCoolingBlockWidth;
  new TGeoTubeSeg("CoolingBlockHoleShape",
						0.,
						fgkSSDCoolingBlockHoleRadius[0],
						0.5*coolingblockholewidth,
						180.-coolingblockholeangle,360.+coolingblockholeangle);
  TVector3* coolingblockholevertex[3];
  coolingblockholevertex[0] = new TVector3();
  coolingblockholevertex[1] = new TVector3(fgkSSDCoolingBlockHoleRadius[0]
					*	TMath::Cos((90.-coolingblockholeangle)*TMath::DegToRad()),
						fgkSSDCoolingBlockHoleRadius[0]
					*	TMath::Sin((90.-coolingblockholeangle)*TMath::DegToRad()));
  coolingblockholevertex[2] = new TVector3(coolingblockholevertex[1]->X(),
					-	coolingblockholevertex[1]->Y());
						GetTriangleShape(coolingblockholevertex,
						coolingblockholewidth,"CoolingBlockTriangleHoleShape");
  TGeoRotation* coolingblockholerot = 
							  new TGeoRotation("CoolingBlockHoleRot",90,0.,0.);
  coolingblockholerot->RegisterYourself();
    new TGeoCompositeShape("CoolingTubeHoleShape",
							  "CoolingBlockTriangleHoleShape:CoolingBlockHoleRot+"
							  "CoolingBlockHoleShape");
  ///////////////////////////
  // Cooling Block Trapezoids
  ///////////////////////////
  const Int_t kvertexnumber = 4;
  const Int_t ktrapezoidnumber = 2;
  TVector3** coolingblocktrapezoidvertex[ktrapezoidnumber];
  for(Int_t i = 0; i<ktrapezoidnumber; i++) coolingblocktrapezoidvertex[i] = 
						new TVector3*[kvertexnumber]; 
  Double_t coolingblockcomponentheight = fgkSSDCoolingBlockHeight[0]
				    -	fgkSSDCoolingBlockHoleCenter
					-	fgkSSDCoolingBlockHoleRadius[0]
					*	TMath::Sin(coolingblockholeangle*TMath::DegToRad());
  Double_t coolingblocktrapezoidlength[ktrapezoidnumber] = 
					{	fgkSSDCoolingBlockLength,
						0.5*(fgkSSDCoolingBlockLength-2.
					*	(fgkSSDCoolingBlockHoleLength[1]
					-	fgkSSDCoolingBlockHoleRadius[1])
					-	fgkSSDCoolingBlockHoleLength[0])}; 
  Double_t coolingblocktrapezoidheigth[ktrapezoidnumber] = 
					{	fgkSSDCoolingBlockHeight[0]-coolingblockcomponentheight
					-	fgkSSDCoolingBlockHeight[1]-fgkSSDCoolingBlockHeight[2]
					-	fgkSSDCoolingBlockHoleRadius[1],
						coolingblockcomponentheight};
  Double_t coolingblocktrapezoidwidth[ktrapezoidnumber]  = 
						{fgkSSDCoolingBlockWidth,fgkSSDCoolingBlockWidth};
  //////////////////////////
  //Vertex Positioning Shape 
  //////////////////////////
  coolingblocktrapezoidvertex[0][0] = new TVector3();
  coolingblocktrapezoidvertex[0][1] = new TVector3(coolingblocktrapezoidlength[0]);
  coolingblocktrapezoidvertex[0][2] = new TVector3(
						0.5*(coolingblocktrapezoidvertex[0][1]->X()
					-	2.*coolingblocktrapezoidlength[1]
					-	fgkSSDCoolingBlockHoleLength[0]));
  coolingblocktrapezoidvertex[0][3] = 
						new TVector3(coolingblocktrapezoidvertex[0][1]->X()
					-	coolingblocktrapezoidvertex[0][2]->X());
  coolingblocktrapezoidvertex[1][0] = new TVector3(); 
  coolingblocktrapezoidvertex[1][1] = new TVector3(coolingblocktrapezoidlength[1]);
  coolingblocktrapezoidvertex[1][2] = 
						new TVector3(coolingblocktrapezoidheigth[1]
					/				 coolingblocktrapezoidheigth[0]
					*	coolingblocktrapezoidvertex[0][2]->X());
  coolingblocktrapezoidvertex[1][3] = 
						new TVector3(coolingblocktrapezoidvertex[1][1]->X());
  char* coolingblocktrapezoidshapename[ktrapezoidnumber] = 
					{"CoolingBlockTrapezoidShape0","CoolingBlockTrapezoidShape1"};
  TGeoArb8* coolingblocktrapezoidshape[ktrapezoidnumber];
  for(Int_t i = 0; i< ktrapezoidnumber; i++) coolingblocktrapezoidshape[i] = 
						GetArbShape(coolingblocktrapezoidvertex[i],
						coolingblocktrapezoidwidth,
						coolingblocktrapezoidheigth[i],
						coolingblocktrapezoidshapename[i]);
  TGeoTranslation* coolingblocktrapezoidtrans = 
						new TGeoTranslation("CoolingBlockTrapezoidTrans",
						coolingblocktrapezoidvertex[0][2]->X(),
						0.0,
						0.5*(coolingblocktrapezoidheigth[0]
					+	coolingblocktrapezoidheigth[1]));
  coolingblocktrapezoidtrans->RegisterYourself();
  TGeoCombiTrans* coolingblocktrapezoidcombitrans = 
						new TGeoCombiTrans("CoolingBlockTrapezoidCombiTrans",
						coolingblocktrapezoidvertex[0][3]->X(),
						fgkSSDCoolingBlockWidth,
						0.5*(coolingblocktrapezoidheigth[0]
					+	coolingblocktrapezoidheigth[1]),
						new TGeoRotation("",180.,0.,0.));
  coolingblocktrapezoidcombitrans->RegisterYourself();
	new TGeoCompositeShape("CoolingBlockTrapezoidCompositeShape",
	"CoolingBlockTrapezoidShape0+CoolingBlockTrapezoidShape1:CoolingBlockTrapezoidTrans+"
	"CoolingBlockTrapezoidShape1:CoolingBlockTrapezoidCombiTrans"); 
  /////////////////////////////
  // Cooling Block Boxes Shapes
  /////////////////////////////
  const Int_t kboxnumber = 3;
  TGeoBBox* coolingblockboxshape[kboxnumber];
  coolingblockboxshape[0] = new TGeoBBox("CoolingBlockBoxShape0",
						0.5*fgkSSDCoolingBlockLength,
						0.5*fgkSSDCoolingBlockWidth,
						0.5*fgkSSDCoolingBlockHoleRadius[1]);
  coolingblockboxshape[1] = new TGeoBBox("CoolingBlockBoxShape1",
						0.5*(fgkSSDCoolingBlockLength
					-	2.*fgkSSDCoolingBlockHoleLength[1]),
						0.5*fgkSSDCoolingBlockWidth,
						0.5*fgkSSDCoolingBlockHeight[2]);
  coolingblockboxshape[2] = new TGeoBBox("CoolingBlockBoxShape2",
						0.5*fgkSSDCoolingBlockLength,
						0.5*fgkSSDCoolingBlockWidth,
						0.5*fgkSSDCoolingBlockHeight[1]);
  TGeoTranslation* coolingblockboxtrans[kboxnumber-1];
  coolingblockboxtrans[0] = new TGeoTranslation("CoolingBlockBoxTrans0",0.,0.,
						0.5*(fgkSSDCoolingBlockHeight[1]
					+	fgkSSDCoolingBlockHoleRadius[1])
					+	fgkSSDCoolingBlockHeight[2]);
  coolingblockboxtrans[1] = new TGeoTranslation("CoolingBlockBoxTrans1",
						0.0,
						0.0,
						0.5*(fgkSSDCoolingBlockHeight[1]
					+	fgkSSDCoolingBlockHeight[2]));
  for(Int_t i=0; i<kboxnumber-1; i++) coolingblockboxtrans[i]->RegisterYourself();
	new TGeoCompositeShape("CoolingBlockBoxCompositeShape",
						   "CoolingBlockBoxShape0:CoolingBlockBoxTrans0+"
	 "CoolingBlockBoxShape1:CoolingBlockBoxTrans1+CoolingBlockBoxShape2");
  ///////////////////////
  // Cooling Block Shape
  //////////////////////
  TGeoCombiTrans* coolingtubeholeshapecombitrans = 
						new TGeoCombiTrans("CoolingTubeHoleShapeCombiTrans",
						0.5*fgkSSDCoolingBlockLength,
						0.5*fgkSSDCoolingBlockWidth,
						fgkSSDCoolingBlockHoleCenter,
						new TGeoRotation("",0.,90.,0.));
  coolingtubeholeshapecombitrans->RegisterYourself();
  TGeoTranslation* coolingblocktrapezoidcompositeshapetrans = 
						new TGeoTranslation("CoolingBlockTrapezoidCompositeShapeTrans",
						0.0,
						0.0,
						0.5*coolingblocktrapezoidheigth[0]+fgkSSDCoolingBlockHeight[1]+
						fgkSSDCoolingBlockHeight[2]+fgkSSDCoolingBlockHoleRadius[1]);
  coolingblocktrapezoidcompositeshapetrans->RegisterYourself();
  TGeoTranslation* coolingblockboxcompositeshapetrans = 
						new TGeoTranslation("CoolingBlockBoxCompositeShapeTrans",
						0.5*fgkSSDCoolingBlockLength,
						0.5*fgkSSDCoolingBlockWidth,
						0.5*fgkSSDCoolingBlockHeight[1]);
  coolingblockboxcompositeshapetrans->RegisterYourself();
  TGeoCompositeShape* ssdoolingblockshape = 
		new TGeoCompositeShape("SSDCoolingBlockShape",	
		"CoolingBlockBoxCompositeShape:CoolingBlockBoxCompositeShapeTrans+"
		"CoolingBlockTrapezoidCompositeShape:CoolingBlockTrapezoidCompositeShapeTrans-"
		"CoolingTubeHoleShape:CoolingTubeHoleShapeCombiTrans");
  TGeoVolume* ssdcoolingblock = new TGeoVolume("SSDCoolingBlock",
		ssdoolingblockshape,fSSDAlCoolBlockMedium);
  return ssdcoolingblock;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCarbonFiberJunction(Double_t width){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Carbon Fiber Junction
  /////////////////////////////////////////////////////////////
  const Int_t kvertexnumber = 4;
  TVector3* vertex[kvertexnumber];
  vertex[0] = new TVector3();
  vertex[1] = new TVector3(fgkCarbonFiberJunctionLength,0.);
  vertex[2] = new TVector3(fgkCarbonFiberJunctionLength
	        -	  fgkCarbonFiberJunctionEdge[1]
			*	  TMath::Cos(fgkCarbonFiberJunctionAngle[1]*TMath::DegToRad()),
				  fgkCarbonFiberJunctionEdge[1]*TMath::Sin(fgkCarbonFiberJunctionAngle[1]
			*     TMath::DegToRad()));
  vertex[3] = new TVector3(fgkCarbonFiberJunctionEdge[0]
			*	  TMath::Cos(fgkCarbonFiberJunctionAngle[0]*TMath::DegToRad()),
				  fgkCarbonFiberJunctionEdge[0]
			*	  TMath::Sin(fgkCarbonFiberJunctionAngle[0]*TMath::DegToRad()));
  TGeoArb8* carbonfiberjunctionshapepiece = 
						new TGeoArb8("CarbonFiberJunctionShapePiece",0.5*width);
  //////////////////////////////////
  //Setting the vertices in TGeoArb8
  //////////////////////////////////
  for(Int_t i = 0; i<2*kvertexnumber; i++)
	carbonfiberjunctionshapepiece->SetVertex(i,
							vertex[(i < kvertexnumber ? i: i-kvertexnumber)]->X(),
							vertex[(i < kvertexnumber ? i : i-kvertexnumber)]->Y());
  TGeoRotation* carbonfiberjunctionrot = 
						new TGeoRotation("CarbonFiberJunctionRot",
										  180.,
										  180.,
										  180-2.*fgkCarbonFiberJunctionAngle[0]); 
  TGeoVolume* carbonfiberjunctionpiece = 
						new TGeoVolume("CarbonFiberJunctionPiece",
						carbonfiberjunctionshapepiece,fSSDCarbonFiberMedium);
  TGeoVolume* carbonfiberjunction = 
						new TGeoVolumeAssembly("CarbonFiberJunction");
  carbonfiberjunctionpiece->SetLineColor(fColorCarbonFiber);
  carbonfiberjunction->AddNode(carbonfiberjunctionpiece,1);
  carbonfiberjunction->AddNode(carbonfiberjunctionpiece,2,carbonfiberjunctionrot);
  return carbonfiberjunction;
}
/////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCarbonFiberJunctionAssembly(){
  /////////////////////////////////////////////////////////////
  // Method generating SSD Carbon Fiber Junction Assembly    
  /////////////////////////////////////////////////////////////
  SetCarbonFiberJunctionCombiTransMatrix();
  TGeoVolume* carbonfiberjunctionassembly = 
						  new TGeoVolumeAssembly("CarbonFiberJunctionAssembly");
  TGeoVolume* carbonfiberjunction = 
						  GetCarbonFiberJunction(fgkCarbonFiberJunctionWidth);
  for(Int_t i=0; i<fgkCarbonFiberJunctionCombiTransNumber;i++) 
    carbonfiberjunctionassembly->AddNode(carbonfiberjunction,i+1,
										 fCarbonFiberJunctionCombiTransMatrix[i]);
  return carbonfiberjunctionassembly;
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
/////////////////////////////////////////////////////////////////////////////////
TList* AliITSv11GeometrySSD::GetEndLadderCarbonFiberJunctionAssembly(){
  /////////////////////////////////////////////////////////////
  // Method generating the End Ladder Carbon Fiber Junction Assembly   
  /////////////////////////////////////////////////////////////  
  const Int_t kendlabbercarbonfiberjunctionumber = 2;
  TGeoVolume* endladdercarbonfiberjunctionassembly[kendlabbercarbonfiberjunctionumber];
  endladdercarbonfiberjunctionassembly[0] = 
				new TGeoVolumeAssembly("EndLadderCarbonFiberJunctionAssembly1");
  endladdercarbonfiberjunctionassembly[1] = 
				new TGeoVolumeAssembly("EndLadderCarbonFiberJunctionAssembly2");
  TGeoVolume** endladdercarbonfiberjunction[kendlabbercarbonfiberjunctionumber];
  for(Int_t i=0; i<kendlabbercarbonfiberjunctionumber; i++) 
						   endladdercarbonfiberjunction[i] = new TGeoVolume*[2];
  for(Int_t i=0; i<kendlabbercarbonfiberjunctionumber; i++){
    endladdercarbonfiberjunction[i][0] = 
		  GetCarbonFiberJunction(fgkEndLadderCarbonFiberLowerJunctionLength[i]);
    endladdercarbonfiberjunction[i][1] = 
		  GetCarbonFiberJunction(fgkEndLadderCarbonFiberUpperJunctionLength[i]);
  }
  TList* endladdercarbonfiberjunctionlist = new TList();
  for(Int_t i=0; i<kendlabbercarbonfiberjunctionumber; i++){
    SetEndLadderCarbonFiberJunctionCombiTransMatrix(i);
    for(Int_t j=0; j<fgkCarbonFiberJunctionCombiTransNumber; j++)
      endladdercarbonfiberjunctionassembly[i]->AddNode(j==2 ? 
						 endladdercarbonfiberjunction[i][1] : 
						 endladdercarbonfiberjunction[i][0],
						 j+1,fEndLadderCarbonFiberJunctionCombiTransMatrix[j]);
    endladdercarbonfiberjunctionlist->Add(endladdercarbonfiberjunctionassembly[i]);
  }
  return endladdercarbonfiberjunctionlist;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCarbonFiberSupport(){
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
  TGeoArb8* carbonfibersupportshape[kshapesnumber]; 
  Double_t width[2] = {fgkCarbonFiberSupportWidth,fgkCarbonFiberSupportWidth};
  Double_t carbonfibersupportheight = 
	  carbonfibersupportxaxisEdgeproj*TMath::Tan(fgkCarbonFiberJunctionAngle[0]
	  *TMath::DegToRad());
  for(Int_t i = 0; i< kshapesnumber; i++) carbonfibersupportshape[i] = 
					GetArbShape(vertexposition[i],width,carbonfibersupportheight,
									carbonfibersupportshapename[i],i==0 ? 1: -1);
  /////////////////////////////////////
  //Setting Translations and Rotations: 
  /////////////////////////////////////
  TGeoTranslation* carbonfibersupporttrans = 
						  new TGeoTranslation("CarbonFiberSupportTrans",
											  0.0,0.0,0.5*carbonfibersupportheight);
  carbonfibersupporttrans->RegisterYourself();
  TGeoRotation* carbonfibercompshaperot[2];
  carbonfibercompshaperot[0] = new TGeoRotation("CarbonFiberCompShapeRot1",
											  0.0,180.0,0.0);
  carbonfibercompshaperot[1] = new TGeoRotation("CarbonFiberCompShapeRot2",
										  90.,-fgkCarbonFiberTriangleAngle,-90.);
  Double_t transvector[3] = {fgkCarbonFiberTriangleLength
						  *  TMath::Cos(fgkCarbonFiberTriangleAngle
						  *	 TMath::DegToRad()),0.,-fgkCarbonFiberTriangleLength
						  *	 TMath::Sin(fgkCarbonFiberTriangleAngle
						  *	 TMath::DegToRad())};
  TGeoCombiTrans* carbonfibersupportcombitrans = 
							   new TGeoCombiTrans("CarbonFiberSupportCombiTrans",
							   transvector[0],2.*symmetryplaneposition
						  +	   transvector[1],transvector[2],
							   new TGeoRotation((*carbonfibercompshaperot[1])
						  *	   (*carbonfibercompshaperot[0])));
  carbonfibersupportcombitrans->RegisterYourself();
////////////////////////////////////////////////////////////////////////////////
  TGeoCompositeShape* carbonfibersupportcompshape = 
							new TGeoCompositeShape("CarbonFiberSupportCompShape",
							"CarbonFiberSupportShape1:CarbonFiberSupportTrans+"
							"CarbonFiberSupportShape2:CarbonFiberSupportTrans");
  TGeoVolume* carbonfibersupport = new TGeoVolume("CarbonFiberSupport",
						   carbonfibersupportcompshape,fSSDCarbonFiberMedium);
  carbonfibersupport->SetLineColor(fColorCarbonFiber);
  TGeoVolume* carbonfibersupportassembly = 
						    new TGeoVolumeAssembly("CarbonFiberSupportAssembly");
  carbonfibersupportassembly->AddNode(carbonfibersupport,1);
  carbonfibersupportassembly->AddNode(carbonfibersupport,2,
									  carbonfibersupportcombitrans);
  return carbonfibersupportassembly;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCarbonFiberLowerSupport(Int_t ikind, Bool_t EndLadder){
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
  char* carbonfiberlowersupportname[kshapesnumber] = 
			  {"CarbonFiberLowerSupportShape1","CarbonFiberLowerSupportShape2"};
  TGeoArb8* carbonfiberlowersupportshape[kshapesnumber];
  for(Int_t i = 0; i< kshapesnumber; i++) carbonfiberlowersupportshape[i] = 
								GetArbShape(vertexposition[i],width,
											fgkCarbonFiberLowerSupportHeight,
											carbonfiberlowersupportname[i]);
  ///////////////////////////////////////////////////////
  TGeoTranslation* carbonfiberlowersupporttrans[kshapesnumber];
  carbonfiberlowersupporttrans[0] = 
						new TGeoTranslation("CarbonFiberLowerSupportTrans1",
						0.0,
						vertexposition[1][3]->Y()+vertexposition[1][2]->Y(),
						0.0);
  carbonfiberlowersupporttrans[1] = 
						new TGeoTranslation("CarbonFiberLowerSupportTrans2",
						0.0,
				-		vertexposition[1][3]->Y()-vertexposition[1][2]->Y(),
						0.0);
  for(Int_t i = 0; i< kshapesnumber; i++) 
						carbonfiberlowersupporttrans[i]->RegisterYourself(); 
  ///////////////////////////////////////////////////////
  TGeoCompositeShape* carbonfiberlowersupportcompshape; 
  if(EndLadder==false)
    carbonfiberlowersupportcompshape = 
				new TGeoCompositeShape("CarbonFiberLowerSupportCompShape",
				"CarbonFiberLowerSupportShape2+"
				"CarbonFiberLowerSupportShape1:CarbonFiberLowerSupportTrans1");
  else
    if(ikind==0)
      carbonfiberlowersupportcompshape = 
						  (TGeoCompositeShape*)carbonfiberlowersupportshape[0];
    else
      carbonfiberlowersupportcompshape = 
	new TGeoCompositeShape("CarbonFiberLowerSupportCompShape",
				 "CarbonFiberLowerSupportShape1+"
				 "CarbonFiberLowerSupportShape1:CarbonFiberLowerSupportTrans1"); 
  TGeoVolume* carbonfiberlowersupport = new TGeoVolume("CarbonFiberLowerSupport",
					  carbonfiberlowersupportcompshape,fSSDCarbonFiberMedium);
  carbonfiberlowersupport->SetLineColor(fColorCarbonFiber);
  return carbonfiberlowersupport;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCarbonFiberAssemblySupport(){
  /////////////////////////////////////////////////////////////
  // Method generating the Carbon Fiber Assembly Support   
  /////////////////////////////////////////////////////////////  
  SetCarbonFiberAssemblyCombiTransMatrix();
  TGeoVolume* carbonfiberassemblysupport = 
						new TGeoVolumeAssembly("CarbonFiberAssembly");
  TGeoVolume* carbonfiberassemblyvolumes[fgkCarbonFiberAssemblyCombiTransNumber];
  carbonfiberassemblyvolumes[0] = GetCarbonFiberJunctionAssembly();
  carbonfiberassemblyvolumes[1] = GetCarbonFiberSupport();
  carbonfiberassemblyvolumes[2] = GetCarbonFiberLowerSupport();
  for(Int_t i=0; i<fgkCarbonFiberAssemblyCombiTransNumber;i++) 
    carbonfiberassemblysupport->AddNode(carbonfiberassemblyvolumes[i],1,
						fCarbonFiberAssemblyCombiTransMatrix[i]);
  return carbonfiberassemblysupport;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCoolingTubeSupport(){
  /////////////////////////////////////////////////////////////
  // Method generating the Cooling Tube Support
  /////////////////////////////////////////////////////////////  
  const Int_t kvertexnumber = 3;
  Double_t phi = TMath::ASin(0.5*fgkCoolingTubeSupportHeight
			   /			 fgkCoolingTubeSupportRmax)*TMath::RadToDeg();
						new TGeoTubeSeg("CoolingTubeSegShape",0.0,
							fgkCoolingTubeSupportRmax,
							0.5*fgkCoolingTubeSupportWidth,phi,
							360-phi);
						new TGeoTube("CoolingTubeHoleShape",0.0,
							fgkCoolingTubeSupportRmin,
							0.5*fgkCoolingTubeSupportWidth);
  TVector3* vertexposition[kvertexnumber];
  ///////////////////////////
  //Shape Vertex Positioning
  ///////////////////////////
  vertexposition[0] = new TVector3();
  vertexposition[1] = new TVector3(fgkCoolingTubeSupportRmax
					*		TMath::Cos(phi*TMath::DegToRad()),
							fgkCoolingTubeSupportRmax
					*		TMath::Sin(phi*TMath::DegToRad()));
  vertexposition[2] = new TVector3(vertexposition[1]->X(),
					-			   vertexposition[1]->Y());
  GetTriangleShape(vertexposition,
									   fgkCoolingTubeSupportWidth,
									   "CoolingTubeTriangleShape");
  Double_t* boxorigin = new Double_t[3];
  Double_t boxlength = fgkCoolingTubeSupportLength-fgkCoolingTubeSupportRmax
					 - vertexposition[1]->X();
  boxorigin[0] = vertexposition[1]->X()+0.5*boxlength, boxorigin[1] = boxorigin[2] = 0.;
						new TGeoBBox("CoolingTubeBoxShape",0.5*boxlength,
							0.5*fgkCoolingTubeSupportHeight,
							0.5*fgkCoolingTubeSupportWidth,boxorigin);
  TGeoCompositeShape* coolingtubesupportshape = 
						new TGeoCompositeShape("CoolingTubeSupportShape",
						 "(CoolingTubeSegShape+CoolingTubeTriangleShape"
						 "+CoolingTubeBoxShape)-CoolingTubeHoleShape"); 
  TGeoVolume* coolingtubesupport = new TGeoVolume("CoolingTubeSupport",
						 coolingtubesupportshape,fSSDTubeHolderMedium);
  return coolingtubesupport;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCoolingTubeSupportAssembly(){
  /////////////////////////////////////////////////////////////
  // Method generating the Cooling Tube Support Assembly
  /////////////////////////////////////////////////////////////  
  TGeoVolume* coolingtubesupportassembly = 
						   new TGeoVolumeAssembly("CoolingTubeSupportAssembly");
  TGeoVolume* coolingtubesupport = GetCoolingTubeSupport();
  SetCoolingTubeSupportCombiTransMatrix();
  for(Int_t i=0; i<fgkCoolingTubeSupportCombiTransNumber;i++) 
    coolingtubesupportassembly->AddNode(coolingtubesupport,i+1,
										 fCoolingTubeSupportCombiTransMatrix[i]);
  return coolingtubesupportassembly;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCoolingTube() const{
  /////////////////////////////////////////////////////////////
  // Method generating the Cooling Tube 
  /////////////////////////////////////////////////////////////  
  TGeoVolume* coolingtubeassembly = new TGeoVolumeAssembly("CoolingTubeAssembly");
  TGeoTube *coolingtubeshape = new TGeoTube("CoolingTubeShape", fgkCoolingTubeRmin, 
						fgkCoolingTubeRmax, fgkCoolingTubeLength/2.0);
  TGeoVolume* coolingtube = new TGeoVolume("CoolingTube",
						 coolingtubeshape,fSSDCoolingTubePhynox);
  TGeoTube *coolingtubeinteriorshape = new TGeoTube("CoolingTubeInteriorShape", 
						0, fgkCoolingTubeRmin, 
						fgkCoolingTubeLength/2.0);
  TGeoVolume *coolingtubeinterior = new TGeoVolume("CoolingTubeInterior",
						   coolingtubeinteriorshape,fSSDCoolingTubeWater);
  coolingtubeassembly->AddNode(coolingtube,1);
  coolingtubeassembly->AddNode(coolingtubeinterior,2);
  return coolingtubeassembly;
}  
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetCoolingTubeAssembly(){
  /////////////////////////////////////////////////////////////
  // Method generating the Cooling Tube Assembly
  /////////////////////////////////////////////////////////////  
  TGeoVolume* coolingtubeassembly =   new TGeoVolumeAssembly("CoolingTubeAssembly");
  TGeoVolume* coolingtube = GetCoolingTube();
  SetCoolingTubeCombiTransMatrix();
  for(Int_t i=0; i<fgkCoolingTubeCombiTransNumber;i++) 
    coolingtubeassembly->AddNode(coolingtube,i+1,fCoolingTubeTransMatrix[i]);
  return coolingtubeassembly;
}
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetLadderSegment(Int_t iChipCablesHeight){
  /////////////////////////////////////////////////////////////
  // Method generating the basic Ladder Segment element which will be replicated   
  /////////////////////////////////////////////////////////////  
  TGeoVolume*laddersegment = new TGeoVolumeAssembly("LadderSegment");
  TGeoVolume* laddersegmentvolumes[fgkLadderSegmentCombiTransNumber];
  laddersegmentvolumes[0] = GetCarbonFiberAssemblySupport();
  laddersegmentvolumes[1] = GetSSDModule(iChipCablesHeight);
  laddersegmentvolumes[2] = GetSSDSensorSupportAssembly(iChipCablesHeight);
  laddersegmentvolumes[3] = GetCoolingTubeSupportAssembly();
	 laddersegmentvolumes[4] = GetCoolingTubeAssembly(); 
  SetLadderSegmentCombiTransMatrix();
  for(Int_t i=0; i<fgkLadderSegmentCombiTransNumber; i++) 
   laddersegment->AddNode(laddersegmentvolumes[i],1,
							 fLadderSegmentCombiTransMatrix[i]);
  return laddersegment;
 }
////////////////////////////////////////////////////////////////////////////////
TList* AliITSv11GeometrySSD::GetEndLadderSegment(){
  /////////////////////////////////////////////////////////////
  // Method generating the Terminal Ladder Segment  
  /////////////////////////////////////////////////////////////  
  const Int_t kendladdersegmentnumber = 2;
  TList* endladdercarbonfiberjunctionlist = GetEndLadderCarbonFiberJunctionAssembly();
  TGeoVolume* endladdersegment[kendladdersegmentnumber];
  endladdersegment[0] = new TGeoVolumeAssembly("EndLadderSegment1");
  endladdersegment[1] = new TGeoVolumeAssembly("EndLadderSegment2");
  TGeoVolume** laddersegmentvolumes[kendladdersegmentnumber];
  const Int_t kladdersegmentvolumenumber = 4;
  for(Int_t i=0; i<kendladdersegmentnumber; i++) laddersegmentvolumes[i] = 
						new TGeoVolume*[kladdersegmentvolumenumber];
  laddersegmentvolumes[0][0] = (TGeoVolume*)endladdercarbonfiberjunctionlist->At(0);
  laddersegmentvolumes[0][1] = GetCarbonFiberSupport();
  laddersegmentvolumes[0][2] = GetSSDMountingBlock();
  laddersegmentvolumes[0][3] = GetCarbonFiberLowerSupport(0,true);
  laddersegmentvolumes[1][0] = (TGeoVolume*)endladdercarbonfiberjunctionlist->At(1);
  laddersegmentvolumes[1][1] = laddersegmentvolumes[0][1];
  laddersegmentvolumes[1][2] = laddersegmentvolumes[0][2];
  laddersegmentvolumes[1][3] = GetCarbonFiberLowerSupport(1,true);
  TList* endladdersegmentlist = new TList();
  for(Int_t i=0; i<kendladdersegmentnumber; i++){
    SetEndLadderSegmentCombiTransMatrix(i);
    for(Int_t j=0; j<kladdersegmentvolumenumber; j++)
      endladdersegment[i]->AddNode(laddersegmentvolumes[i][j],1,
								   fEndLadderSegmentCombiTransMatrix[j]);
    endladdersegmentlist->Add(endladdersegment[i]);
  }
  return endladdersegmentlist;
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
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetLadder(Int_t iLayer){
  /////////////////////////////////////////////////////////////
  // Method generating the Layer5 or Layer6 Ladder  
  /////////////////////////////////////////////////////////////  
  TGeoVolume* ladder = new TGeoVolumeAssembly(iLayer==5 ? "ITSssdLay5Ladd" 
														: "ITSssdLay6Ladd");
  TGeoVolume* laddersegment[2];
  laddersegment[0] = GetLadderSegment(0);
  laddersegment[1] = GetLadderSegment(1);
  TList* endladdersegmentlist = GetEndLadderSegment();
  Double_t beamaxistransvector = fgkCarbonFiberJunctionWidth;
  Int_t ssdlaysensorsnumber = (iLayer==5 ? 
						fgkSSDLay5SensorsNumber : 
						fgkSSDLay6SensorsNumber);
  for(Int_t i=0; i<ssdlaysensorsnumber; i++) ladder->AddNode(i%2==0 ? 
						laddersegment[iLayer==5 ? 0 : 1] : 
						laddersegment[iLayer==5 ? 1 : 0],
						ssdlaysensorsnumber-i-1,new TGeoTranslation(-0.5*fgkCarbonFiberTriangleLength,
						beamaxistransvector*i,0.));
  ladder->AddNode((TGeoVolume*)endladdersegmentlist->At(0),1,
						new TGeoTranslation("",-0.5*fgkCarbonFiberTriangleLength,
						fgkCarbonFiberJunctionWidth*ssdlaysensorsnumber,0.));
  ladder->AddNode((TGeoVolume*)endladdersegmentlist->At(1),1,
						new TGeoCombiTrans("",0.5*fgkCarbonFiberTriangleLength,
						0.,0.,new TGeoRotation("",180.,0.,0.)));
/////////////////////////////////////////////////////////////////////////////						
/// Placing Ladder Cables
/////////////////////////////////////////////////////////////////////////////		
  SetLadderCableCombiTransMatrix(iLayer);
  Int_t sidecablenumber[2] = {0,0};
  switch(iLayer){
	case 5: 
		sidecablenumber[0] = fgkSSDLay5SensorsNumber/2+1; 
		sidecablenumber[1] = sidecablenumber[0]-2;
		break;
    case 6:
		sidecablenumber[0] = (fgkSSDLay6SensorsNumber-1)/2+1;
		sidecablenumber[1] = sidecablenumber[0]-1;
		break;
  }
  const Double_t* carbonfibertomoduleposition = 
							 fLadderSegmentCombiTransMatrix[1]->GetTranslation();
  Double_t ssdendladdercablelength[4];
  ssdendladdercablelength[0] = carbonfibertomoduleposition[1]
							 + fgkSSDSensorLength
							 - fgkSSDModuleStiffenerPosition[1]
							 - fgkSSDStiffenerWidth 
							 - fgkSSDFlexWidth[0]
							 + fgkEndLadderCarbonFiberLowerJunctionLength[1];
  ssdendladdercablelength[1] = carbonfibertomoduleposition[1]
							 + fgkSSDModuleStiffenerPosition[1]
							 + fgkSSDStiffenerWidth
							 + fgkEndLadderCarbonFiberLowerJunctionLength[1];
  ssdendladdercablelength[2] = ssdendladdercablelength[1]
							 - fgkEndLadderCarbonFiberLowerJunctionLength[1]
							 + fgkEndLadderCarbonFiberLowerJunctionLength[0];
  ssdendladdercablelength[3] = fgkCarbonFiberJunctionWidth-(fgkSSDSensorLength
							 + carbonfibertomoduleposition[1]
							 - fgkSSDModuleStiffenerPosition[1]
							 - fgkSSDStiffenerWidth)
							 + fgkEndLadderCarbonFiberLowerJunctionLength[1];
  TList* laddercableassemblylist[4];
  const Int_t kendladdercablesnumber = 4;
  for(Int_t i=0; i<kendladdercablesnumber; i++){
	  laddercableassemblylist[i] = 
	  GetLadderCableAssemblyList(sidecablenumber[i<2?0:1],
								   ssdendladdercablelength[i]);
	  ladder->AddNode((TGeoVolume*)laddercableassemblylist[i]->At(i%2==0?0:1),
									i<2?1:2,fLadderCableCombiTransMatrix[i]);
  }
  return ladder;
}								  
////////////////////////////////////////////////////////////////////////////////
TGeoVolume* AliITSv11GeometrySSD::GetLayer(Int_t iLayer){
  /////////////////////////////////////////////////////////////
  // Method generating the Layer5 or Layer6  
  /////////////////////////////////////////////////////////////  
  TGeoVolume* layer = new TGeoVolumeAssembly(iLayer==5 ? "ITSssdLayer5" 
													   : "ITSssdLayer6");
  TGeoVolume* ladder = GetLadder(iLayer);
  /////////////////////////////////////////////////////
  // Setting the CombiTransformation to pass ITS center 
  /////////////////////////////////////////////////////
  Double_t itscentertransz = iLayer==5 ? fgkSSDLay5LadderLength
                           -             fgkLay5CenterITSPosition: 
                                         fgkSSDLay6LadderLength
                           -             fgkLay6CenterITSPosition;
  Double_t itssensorytrans = fgkSSDModuleCoolingBlockToSensor
                           + 0.5*        fgkCoolingTubeSupportHeight
                           -(iLayer==5 ? fgkSSDSensorSideSupportHeight[1]
                           -             fgkSSDSensorSideSupportHeight[0]: 0.);
  TGeoRotation* itscenterrot[3];
  itscenterrot[0] = new TGeoRotation("ITSCenterRot0",90.,180.,-90.);
  itscenterrot[1] = new TGeoRotation("ITSCenterRot1",0.,90.,0.);
  itscenterrot[2] = new TGeoRotation((*itscenterrot[1])*(*itscenterrot[0]));
  TGeoCombiTrans* itscentercombitrans = new TGeoCombiTrans("ITSCenterCombiTrans",0.,
                                        fgkSSDMountingBlockHeight[1]+itssensorytrans,
                                        fgkEndLadderCarbonFiberLowerJunctionLength[1]
                                      - itscentertransz,itscenterrot[2]);
  /////////////////////////////////////////////////////
  // Setting the Ladder into the Layer 
  /////////////////////////////////////////////////////
  TGeoCombiTrans* lay5laddercombitrans[fgkSSDLay5LadderNumber];
  TGeoCombiTrans* lay6laddercombitrans[fgkSSDLay6LadderNumber];
  TGeoHMatrix* lay5ladderhmatrix[fgkSSDLay5LadderNumber];
  TGeoHMatrix* lay6ladderhmatrix[fgkSSDLay6LadderNumber];
  if(iLayer==5){
      Double_t lay5ladderangleposition = 360./fgkSSDLay5LadderNumber;    
      char lay5laddercombitransname[30], lay5ladderrotname[30];
      for(Int_t i=0; i<fgkSSDLay5LadderNumber;i++){
      sprintf(lay5laddercombitransname,"Lay5LadderCombiTrans%i",i);
      sprintf(lay5ladderrotname,"Lay5LaddeRot%i",i);
      Double_t lay5layerradius = (i%2==0 ? fgkSSDLay5RadiusMin: fgkSSDLay5RadiusMax);
      lay5laddercombitrans[i] = new TGeoCombiTrans(lay5laddercombitransname,
                                lay5layerradius *	TMath::Cos((i+1)
                              * lay5ladderangleposition*TMath::DegToRad()),
                                lay5layerradius *	TMath::Sin((i+1)
                              * lay5ladderangleposition*TMath::DegToRad()),0.,
                                new TGeoRotation(lay5ladderrotname,(i+1)
                              * lay5ladderangleposition-90,0.,0.));
      lay5ladderhmatrix[i] = new TGeoHMatrix((*lay5laddercombitrans[i])
                           * (*itscentercombitrans));
      layer->AddNode(ladder,i,lay5ladderhmatrix[i]);            
      }
  }
  else{
      Double_t lay6ladderangleposition = 360./fgkSSDLay6LadderNumber;    
      char lay6laddercombitransname[30], lay6ladderrotname[30];
      for(Int_t i=0; i<fgkSSDLay6LadderNumber;i++){
      sprintf(lay6laddercombitransname,"Lay6LadderCombiTrans%i",i);
      sprintf(lay6ladderrotname,"Lay6LaddeRot%i",i);
      Double_t lay6layerradius = (i%2==0 ? fgkSSDLay6RadiusMin: fgkSSDLay6RadiusMax);
      lay6laddercombitrans[i] = new TGeoCombiTrans(lay6laddercombitransname,
                                lay6layerradius *	TMath::Cos((i+1)
                              * lay6ladderangleposition*TMath::DegToRad()),
                                lay6layerradius *	TMath::Sin((i+1)
                              * lay6ladderangleposition*TMath::DegToRad()),0.,
                                new TGeoRotation(lay6ladderrotname,(i+1)
                              * lay6ladderangleposition-90,0.,0.));
      lay6ladderhmatrix[i] = new TGeoHMatrix((*lay6laddercombitrans[i])
                           * (*itscentercombitrans));
      layer->AddNode(ladder,i,lay6ladderhmatrix[i]);            
      }
  }
  return layer;
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
	 fMotherVol = moth;
	 moth->AddNode(GetLayer(5),1,0);
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
	 fMotherVol = moth;
	 moth->AddNode(GetLayer(6),1,0);
 }
////////////////////////////////////////////////////////////////////////////////
TGeoArb8* AliITSv11GeometrySSD::GetTrapezoidShape(TVector3* vertexpos[], 
							Double_t* width, Double_t height, char* shapename) const{
  /////////////////////////////////////////////////////////////
  // Method generating a trapezoid shape 
  /////////////////////////////////////////////////////////////
  const Int_t kvertexnumber = 4;
  const Int_t ktransvectnumber = 2;
  TVector3* vertex[kvertexnumber];
  TVector3* transvector[2];
  for(Int_t i=0; i<ktransvectnumber; i++) 
									transvector[i] = new TVector3(0.,width[i]);
  /////////////////////////////////////////////////////////////
  //Setting the vertices
  /////////////////////////////////////////////////////////////
  vertex[0] = new TVector3(*vertexpos[0]);
  vertex[1] = new TVector3(*vertexpos[1]);
  vertex[2] = new TVector3(*vertex[1]+*transvector[1]);
  vertex[3] = new TVector3(*vertex[0]+*transvector[0]);
  TGeoArb8* trapezoidshape = new TGeoArb8(shapename,0.5*height);
  for(Int_t i=0; i<2*kvertexnumber; i++) 
	trapezoidshape->SetVertex(i,vertex[(i<kvertexnumber ? i : i-kvertexnumber)]->X(),
								vertex[(i<kvertexnumber ? i : i-kvertexnumber)]->Y()); 
  return trapezoidshape;
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
  return arbshape;
} 
////////////////////////////////////////////////////////////////////////////////
TGeoArb8* AliITSv11GeometrySSD::GetTriangleShape(TVector3* vertexpos[], 
											  Double_t height, char* shapename) const{
  /////////////////////////////////////////////////////////////
  // Method generating a triangle shape 
  /////////////////////////////////////////////////////////////
  const Int_t kvertexnumber = 4;
  TVector3* vertex[kvertexnumber];
//////////////////////////////////////
  //Setting the vertices for TGeoArb8
  ////////////////////////////////////
  for(Int_t i = 0; i<kvertexnumber; i++)  
  vertex[i] = new TVector3(i!=kvertexnumber-1?*vertexpos[i]:*vertex[kvertexnumber-1-i]);
  TGeoArb8* triangleshape = new TGeoArb8(shapename,0.5*height);
  for(Int_t i = 0; i<2*kvertexnumber; i++) 
	triangleshape->SetVertex(i,vertex[(i < kvertexnumber ? i: i-kvertexnumber)]->X(),
							   vertex[(i < kvertexnumber ? i : i-kvertexnumber)]->Y());
  return triangleshape;
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
  return reflectedvector;
}
////////////////////////////////////////////////////////////////////////////////
TGeoCombiTrans* AliITSv11GeometrySSD::AddTranslationToCombiTrans(TGeoCombiTrans* ct,
                                                       Double_t dx,
                                                       Double_t dy,
                                                       Double_t dz) const{
  /////////////////////////////////////////////////////////////
  // Add a dx,dy,dz translation to the initial TGeoCombiTrans
  /////////////////////////////////////////////////////////////
  TGeoCombiTrans* combiTrans = new TGeoCombiTrans(*ct);
  const Double_t *vect = combiTrans->GetTranslation();
  Double_t newvect[3] = {vect[0]+dx, vect[1]+dy, vect[2]+dz};
  combiTrans->SetTranslation(newvect);
  return combiTrans;
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
  fSSDSensorMedium = Silicon;
  ///////////////////////////////////
  // Silicon Mixture for Sensor
  /////////////////////////////////// 
  Int_t SiMixtureNumber = 1;
  TGeoMaterial* SiMixtureMaterial = new TGeoMaterial("SiMixtureMaterial");
  TGeoMedium* SiliconMixture = new TGeoMedium("SiliconMixture",SiMixtureNumber,SiMixtureMaterial);
  fSSDChipMedium = SiliconMixture;
  ///////////////////////////////////
  // Stiffener Components Materials
  /////////////////////////////////// 
  Int_t K1100Number = 2;
  TGeoMaterial* K1100Material = new TGeoMaterial("K1100Material");
  TGeoMedium* K1100 = new TGeoMedium("K1100",K1100Number,K1100Material);
  fSSDStiffenerMedium = K1100;
  ///////////////////////////  
  // Stiffener Connectors 
  ///////////////////////////  
  Int_t SnMaterialNumber = 3;
  TGeoMaterial* SnMaterial = new TGeoMaterial("SnMaterial");
  TGeoMedium* SnMedium = new TGeoMedium("SnMedium",SnMaterialNumber,
                                                SnMaterial);
  fSSDStiffenerConnectorMedium = SnMedium;
  ////////////////////////////////  
  // Stiffener 0603-1812 Capacitor
  ////////////////////////////////  
  Int_t Al2O3Number = 4;
  TGeoMaterial* Al2O3Material = new TGeoMaterial("Al2O3Material");
  TGeoMedium* Al2O3Medium = new TGeoMedium("Al2O3Medium",
                                                Al2O3Number,
                                                Al2O3Material);
  fSSDStiffener0603CapacitorMedium = Al2O3Medium;
  fSSDStiffener1812CapacitorMedium = Al2O3Medium;
  ///////////////////////////  
  // Stiffener Hybrid Wire 
  ///////////////////////////  
  Int_t CuHybridWireNumber = 5;
  TGeoMaterial* CuHybridWireMaterial = new TGeoMaterial("CuHybridWireMaterial");
  TGeoMedium* CuHybridWireMedium = new TGeoMedium("CuHybridWireMedium",
                                                CuHybridWireNumber,
                                                CuHybridWireMaterial);
  fSSDStiffenerHybridWireMedium = CuHybridWireMedium;
  ///////////////////////////  
  // Al for Cooling Block
  ///////////////////////////  
  Int_t AlCoolBlockNumber = 6;
  TGeoMaterial* AlCoolBlockMaterial = new TGeoMaterial("AlCoolBlockMaterial");
  TGeoMedium* AlCoolBlockMedium = new TGeoMedium("AlCoolBlockMedium",
                                                AlCoolBlockNumber,
                                                AlCoolBlockMaterial);
  fSSDAlCoolBlockMedium = AlCoolBlockMedium;
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
  fSSDKaptonChipCableMedium = KaptonMedium;
  fSSDAlTraceChipCableMedium = AlTraceMedium;
  fSSDKaptonFlexMedium = KaptonMedium;
  fSSDAlTraceFlexMedium = AlTraceMedium;
  fSSDKaptonLadderCableMedium = KaptonMedium;
  fSSDAlTraceLadderCableMedium = AlTraceMedium;
  /////////////////////////////////////////////////////////////////  
  // M55J for Carbon Fiber, CarbonFiber Lower Support and Junction
  //////////////////////////////////////////////////////////////////  
  Int_t M55JNumber = 9;
  TGeoMaterial* M55JMaterial = new TGeoMaterial("M55JMaterial");
  TGeoMedium* M55JMedium = new TGeoMedium("M55JMedium",M55JNumber,
                                           M55JMaterial);
  fSSDCarbonFiberMedium = M55JMedium;
  fSSDMountingBlockMedium = M55JMedium;
  /////////////////////////////////////////////////////////////////  
  // G10 for Detector Leg, TubeHolder
  //////////////////////////////////////////////////////////////////  
  Int_t G10Number = 10;
  TGeoMaterial* G10Material = new TGeoMaterial("G10Material");
  TGeoMedium* G10Medium = new TGeoMedium("G10Medium",G10Number,
                                           G10Material);
  fSSDTubeHolderMedium = G10Medium;
  fSSDSensorSupportMedium = G10Medium;
  /////////////////////////////////////////////////////////////////  
  // Water and Phynox for Cooling Tube
  //////////////////////////////////////////////////////////////////  
  Int_t WaterNumber = 11;
  TGeoMaterial* WaterMaterial = new TGeoMaterial("WaterMaterial");
  TGeoMedium* WaterMedium = new TGeoMedium("WaterMedium",WaterNumber,
                                           WaterMaterial);
  fSSDCoolingTubeWater = WaterMedium;
  Int_t PhynoxNumber = 12;
  TGeoMaterial* PhynoxMaterial = new TGeoMaterial("PhynoxMaterial");
  TGeoMedium* PhynoxMedium = new TGeoMedium("PhynoxMedium",PhynoxNumber,
                                           PhynoxMaterial);
  fSSDCoolingTubePhynox = PhynoxMedium;
}
*/
void AliITSv11GeometrySSD::CreateMaterials(){
///////////////////////////////////
// This part has to be modified
///////////////////////////////////
  ///////////////////////////////////
  // Silicon for Sensor
  /////////////////////////////////// 
  fSSDSensorMedium = GetMedium("Si");
  ///////////////////////////////////
  // Silicon Mixture for Sensor
  /////////////////////////////////// 
  fSSDChipMedium = GetMedium("SPD SI CHIP$");
  fSSDChipGlueMedium = GetMedium("EPOXY$");
  ///////////////////////////////////
  // Stiffener Components Materials
  /////////////////////////////////// 
  fSSDStiffenerMedium = GetMedium("ITSsddCarbonM55J");
  ///////////////////////////  
  // Stiffener Connectors 
  ///////////////////////////  
  fSSDStiffenerConnectorMedium = GetMedium("COPPER");
  ////////////////////////////////  
  // Stiffener 0603-1812 Capacitor
  ////////////////////////////////  
  fSSDStiffener0603CapacitorMedium = GetMedium("SDD ruby sph. Al2O3");
  fSSDStiffener1812CapacitorMedium = GetMedium("SDD ruby sph. Al2O3");
  ///////////////////////////  
  // Stiffener Hybrid Wire 
  ///////////////////////////  
  fSSDStiffenerHybridWireMedium = GetMedium("COPPER");
  ///////////////////////////  
  // Al for Cooling Block
  ///////////////////////////  
  fSSDAlCoolBlockMedium = GetMedium("ITSal");
  //////////////////////////////////////////////////////  
  // Kapton and Al for Chip Cable Flex and Ladder Cables
  //////////////////////////////////////////////////////  
  fSSDKaptonChipCableMedium = GetMedium("KAPTONH(POLYCH2)");
  fSSDAlTraceChipCableMedium = GetMedium("ITSal");
  fSSDKaptonFlexMedium = GetMedium("KAPTONH(POLYCH2)");
  fSSDAlTraceFlexMedium = GetMedium("ITSal");
  fSSDKaptonLadderCableMedium = GetMedium("KAPTONH(POLYCH2)");
  fSSDAlTraceLadderCableMedium = GetMedium("ITSal");
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
  fSSDCoolingTubeWater = GetMedium("WATER");
  fSSDCoolingTubePhynox = GetMedium("INOX$");
}
/////////////////////////////////////////////////////////////////////

