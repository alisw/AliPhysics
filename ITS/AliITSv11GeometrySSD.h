#ifndef ALIITSV11GEOMETRYSSD_H
#define ALIITSV11GEOMETRYSSD_H
//*************************************************************************
// class AliITSv11GeometrySSD
// Enrico Cattaruzza                                     ecattar@ts.infn.it
//*************************************************************************
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */
class TGeoVolume;
class TGeoShape;
class TGeoCombiTrans;
class TGeoMedium;
class TGeoCompositeShape;
class TGeoMatrix;
class TVector3;
class TGeoArb8;
class TList;
class TGeoTube;
class TGeoHMatrix;
class TGeoTranslation;
class TGeoRotation;
class TGeoXtru;
class TGeoVolumeAssembly;
#include "AliITSv11Geometry.h"
class AliITSv11GeometrySSD: public AliITSv11Geometry {
public:
  AliITSv11GeometrySSD();
  virtual ~AliITSv11GeometrySSD(){};
  /////////////////////////////////////////////////////////////////////////
  // Public methods
  /////////////////////////////////////////////////////////////////////////
  void CreateMaterials();	  // Method setting the materials 
  TGeoMedium* GetMedium(const char* mediumName);   // It returns the Medium
  const char*   GetSenstiveVolumeName5() const {return fgSSDsensitiveVolName5;};
  // it returns the Sensitive Volume of Layer 5
  const char*   GetSenstiveVolumeName6() const {return fgSSDsensitiveVolName6;};
  // it returns the Sensitive Volume of Layer 6
  TGeoVolumeAssembly* GetLadderSegment(Int_t i) const {return fladdersegment[i];}; // Get Ladder Segment
  TGeoVolumeAssembly* GetEndLadderSegment(Int_t i) const {return fendladdersegment[i];}; // Get End Ladder Segment 
  TGeoVolume* GetLadder(Int_t i) const {return fladder[i];}; // Get Ladder
//  TGeoVolumeAssembly* GetLadder(Int_t i) {return fladder[i];}; // Get Ladder
  TGeoVolumeAssembly* GetLayer(Int_t i)const {return i==5? fSSDLayer5 : fSSDLayer6;}; // Get Layer
  TGeoVolume** GetEndCapAssembly();     // End Cap Assembly
  void SetLadderSegment();				// Set Ladder Elementary Segment 
  void SetEndLadderSegment();			// Set End Ladder Segment
  void SetLadder();						// Set Ladder
  void SetLayer();						// Set Layer
  void SetSSDCone();                    // Set SSD Cone
  TGeoVolume* SetSSDCables();           // Set SSD Cables
  void Layer5(TGeoVolume* moth);        // Setting Layer 5 into mother volume
  void Layer6(TGeoVolume* moth);        // Setting Layer 6 into mother volume
  void LadderSupportLayer5(TGeoVolume* moth); // Setting Ladder Support of Layer 5
  void LadderSupportLayer6(TGeoVolume* moth); // Setting Ladder Support of Layer 6
  void EndCapSupportSystemLayer5(TGeoVolume* moth); // Setting End Cap Support + End Cap Assembly Layer 5
  void EndCapSupportSystemLayer6(TGeoVolume* moth); // Setting End Cap Support + End Cap Assembly Layer 6
  void SSDCone(TGeoVolume* moth); // Setting SSD Cone;
  void SSDCables(TGeoVolume* moth); // Setting SSD Cables;
private:
  AliITSv11GeometrySSD(const AliITSv11GeometrySSD &source);
  AliITSv11GeometrySSD& operator=(const AliITSv11GeometrySSD &source);

  /////////////////////////////////////////////////////////////////////////////////
  // Names of the Sensitive Volumes of Layer 5 and Layer 6
  /////////////////////////////////////////////////////////////////////////////////
  static const char* fgSSDsensitiveVolName5;       // sens. vol. name for lay. 5
  static const char* fgSSDsensitiveVolName6;       // sens. vol. name for lay. 6
  /////////////////////////////////////////////////////////////////////////////////
  // Variable for Vertical Disalignement of Modules
  /////////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDModuleVerticalDisalignment;  // Vertical Disalignement of Volume
  static const Double_t fgkSSDModuleSideDisalignment;  // Vertical Disalignement of Volume
  static const Double_t fgkSSDLadderVerticalDisalignment;  // Extra space at ladder support for disalignment
  static const Double_t fgkSSDTolerance;  // SSD Tolerance
  /////////////////////////////////////////////////////////////////////////
  // Layer5 (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDLay5LadderLength;    // Ladder Layer5 Length
  static const Int_t fgkSSDLay5SensorsNumber;      // Ladder Layer5 Sensor Number
  static const Int_t fgkSSDLay5LadderNumber;       // Ladder Layer5 Number
  static const Double_t fgkSSDLay5RadiusMin;       // Ladder Layer5 Min Radius
  static const Double_t fgkSSDLay5RadiusMax;       // Ladder Layer5 Max Radius
  static const Double_t fgkLay5CenterITSPosition;  // ITS center position respect
                                                   // to Ladder Layer5
  /////////////////////////////////////////////////////////////////////////
  // Layer6 (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDLay6LadderLength;    // Ladder Layer6 Length
  static const Int_t fgkSSDLay6SensorsNumber;      // Ladder Layer6 Sensor Number
  static const Int_t fgkSSDLay6LadderNumber;       // Ladder Layer6 Number
  static const Double_t fgkSSDLay6RadiusMin;       // Ladder Layer6 Min Radius
  static const Double_t fgkSSDLay6RadiusMax;       // Ladder Layer6 Max Radius 
  static const Double_t fgkLay6CenterITSPosition;  // ITS center position respect
                                                   // to Ladder Layer6
  /////////////////////////////////////////////////////////////////////////
  // SSD Chips and Hybrid
  /////////////////////////////////////////////////////////////////////////
  static const Int_t fgkSSDChipNumber;             // SSD Module Chip Number
  static const Double_t fgkSSDChipLength;          // SSD Module Chip Length
  static const Double_t fgkSSDChipWidth;           // SSD Module Chip Width
  static const Double_t fgkSSDChipHeight;          // SSD Module Chip Height
  static const Double_t fgkSSDChipSeparationLength;// SSD Module Distance between Chips
  static const Double_t fgkSSDChipGlueLength;      // SSD Module Chip Glue Layer Length
  static const Double_t fgkSSDChipGlueWidth;       // SSD Module Chip Glue Layer Width 
  static const Double_t fgkSSDChipGlueHeight;      // SSD Module Chip Glue Layer Height
  /////////////////////////////////////////////////////////////////////////
  // Stiffener Components
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDStiffenerLength;     // SSD Module Stiffener Length
  static const Double_t fgkSSDStiffenerWidth;      // SSD Module Stiffener Width
  static const Double_t fgkSSDStiffenerHeight;     // SSD Module Stiffener Height
  static const Double_t fgkSSDStiffenerToChipDist; // SSD Module Stiffener position respect 
                                                   // to sensor Edge
  static const Double_t fgkSSDCapacitor0603Length; // SSD Stiffener Capacitor 0603 Length
  static const Double_t fgkSSDCapacitor0603Width;  // SSD Stiffener Capacitor 0603 Width
  static const Double_t fgkSSDCapacitor0603Height; // SSD Stiffener Capacitor 0603 Height
  static const Double_t fgkSSDCapacitor0603CapLength; // SSD Stiffener Capacitor 1812 Cap Length 
  static const Double_t fgkSSDCapacitor1812Length; // SSD Stiffener Capacitor 1812 Length 
  static const Double_t fgkSSDCapacitor1812Width;  // SSD Stiffener Capacitor 1812 Width
  static const Double_t fgkSSDCapacitor1812Height; // SSD Stiffener Capacitor 1812 Height
  static const Double_t fgkSSDCapacitor1812CapLength; // SSD Stiffener Capacitor 1812 Cap Length 
  static const Double_t fgkSSDWireLength;          // SSD Stiffener Wire Length
  static const Double_t fgkSSDWireRadius;          // SSD Stiffener Wire Radius
  static const Double_t fgkSSDConnectorPosition[2];// SSD Connector Position respect to Stiffener
  static const Double_t fgkSSDConnectorSeparation; // SSD Connector separation distance
  static const Double_t fgkSSDConnectorLength;     // SSD Stiffener Connector Length
  static const Double_t fgkSSDConnectorWidth;      // SSD Stiffener Connector Width
  static const Double_t fgkSSDConnectorHeight;     // SSD Stiffener Connector Height
  static const Double_t fgkSSDConnectorAlHeight;     // SSD Stiffener Connector Al Height
  static const Double_t fgkSSDConnectorNiHeight;     // SSD Stiffener Connector Ni Height
  static const Double_t fgkSSDConnectorSnHeight;     // SSD Stiffener Connector Sn Height
  /////////////////////////////////////////////////////////////////////////
  // Flex
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDFlexFullLength;      // SSD Flex Full Length
  static const Double_t fgkSSDFlexLength[4];       // SSD Flex Components Length
  static const Double_t fgkSSDFlexWidth[2];        // SSD Flex Components Width
  static const Double_t fgkSSDFlexHeight[2];       // SSD Flex Layers Height
  static const Double_t fgkSSDFlexAngle;           // SSD Flex Angle 
  static const Double_t fgkSSDFlexHoleLength;      // SSD Flex Hole Length
  static const Double_t fgkSSDFlexHoleWidth;       // SSD Flex Hole Width
  static const Double_t fgkSSDEndFlexCompLength[6];// SSD End-Flex Components Length
  static const Double_t fgkSSDEndFlexCompWidth[3]; // SSD End-Flex Components Width
  /////////////////////////////////////////////////////////////////////////////////
  // SSD Ladder Cable 
  /////////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDLadderCableWidth;    // SSD Ladder Cable Width
  static const Double_t fgkSSDLadderCableHeight[2];  // SSD Ladder Cable Height (thickness)
  /////////////////////////////////////////////////////////////////////////
  // SSD Module Components 
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDModuleStiffenerPosition[2]; // SSD Module Stiffener position 
                                                          // respect to Sensor Edge
  static const Double_t fgkSSDModuleSensorSupportDistance;// SSD Module Sensor Support Position 
                                                          // respect to Sensor Edge 
  static const Double_t fgkSSDModuleCoolingBlockToSensor; // SSD Cooling Block Position 
                                                          // respect to sensor
  /////////////////////////////////////////////////////////////////////////
  // Chip Cables
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDChipCablesLength[2]; // SSD Chip Cables Components Length
  static const Double_t fgkSSDChipCablesHeight[4]; // SSD Chip Cables Components Height   
  static const Double_t fgkSSDChipCablesWidth[3];  // SSD Chip Cables Components Width
  /////////////////////////////////////////////////////////////////////////
  // Cooling Block
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDCoolingBlockLength;       // SSD Cooling Block Length
  static const Double_t fgkSSDCoolingBlockWidth;        // SSD Cooling Block Width   
  static const Double_t fgkSSDCoolingBlockHeight[3];    // SSD Cooling Block Heiht
  static const Double_t fgkSSDCoolingBlockHoleRadius[2];// SSD Cooling Block Hole Radius
  static const Double_t fgkSSDCoolingBlockHoleLength[2];// SSD Cooling Block Hole Length 
  static const Double_t fgkSSDCoolingBlockHoleCenter;   // SSD Cooling Block Hole Ceneter Position
  static const Double_t fgkSSDCoolingBlockHoleHeight;   // SSD Cooling Block Hole Height
  /////////////////////////////////////////////////////////////////////////
  // SSD Sensor 
  /////////////////////////////////////////////////////////////////////////
  static const char* fgkSSDSensitiveVolName;           // SSD Name of the Sensitive Part of the Sensor
  static const Double_t fgkSSDSensorLength;            // SSD Sensor Length              
  static const Double_t fgkSSDSensorHeight;            // SSD Sensor Height
  static const Double_t fgkSSDSensorWidth;             // SSD Sensor Width
  static const Double_t fgkSSDSensorOverlap;           // SSD Sensor Beam Axis Overlap
  static const Double_t fgkSSDSensorInsensitiveLength; // SSD Insensitive Part Length
  static const Double_t fgkSSDSensorInsensitiveWidth;  // SSD Insensitive Part Width
  /////////////////////////////////////////////////////////////////////////
  // SSD Sensor Support 
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDSensorSideSupportLength;        // SSD Side Sensor Support Length
  static const Double_t fgkSSDSensorSideSupportWidth;         // SSD Side Sensor Support Width
  static const Double_t fgkSSDSensorSideSupportHeight[2];     // SSD Side Sensor Support Height
  static const Double_t fgkSSDSensorSideSupportThickness[2];  // SSD Side Sensor Support Thickness 
  static const Double_t fgkSSDSensorSideSupportPosition;      // SSD Side Sensor Support Position 
  static const Double_t fgkSSDSensorCenterSupportLength;      // SSD Center Sensor Support Length
  static const Double_t fgkSSDSensorCenterSupportWidth;       // SSD Center Sensor Support Width
  static const Double_t fgkSSDSensorCenterSupportHeight[2];   // SSD Center Sensor Support Height
  static const Double_t fgkSSDSensorCenterSupportThickness[2];// SSD Center Sensor Support Thickness
  static const Double_t fgkSSDSensorCenterSupportPosition;    // SSD Center Sensor Support Position
  static const Int_t fgkSSDSensorSupportCombiTransNumber = 3; // Number of TGeoCombiTrans 
                                                              // for positioning volumes in Sensor Support Assembly       
  /////////////////////////////////////////////////////////////////////////
  //Parameters for Carbon Fiber 
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkCarbonFiberTriangleLength;            // Carbon Fiber Triangle Length 
  static const Double_t fgkCarbonFiberTriangleAngle;             // Carbon Fiber Triangle Angle
  static const Double_t fgkCarbonFiberSupportTopEdgeDist[2];     // Carbon Fiber Support Top Edge Components
  static const Double_t fgkCarbonFiberSupportEdgeLength;         // Carbon Fiber Support Edge Lenght
  static const Double_t fgkCarbonFiberSupportWidth;              // Carbon Fiber Support Edge Width
  static const Double_t fgkCarbonFiberSupportXAxisLength;        // Carbon Fiber Support X Axis Lenght
  static const Double_t fgkCarbonFiberSupportYAxisLength;        // Carbon Fiber Support Y Axis Lenght
  static const Int_t fgkCarbonFiberAssemblyCombiTransNumber = 3; // Number of TGeoCombiTrans 
                                                                 // for positioning volumes in Carbon Fiber Assembly 
  //////////////////////////////////////////////////////////////////////////////
  // Carbon Fiber Junction Parameters
  //////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkCarbonFiberJunctionLength;            // Carbon Fiber Junction Length             
  static const Double_t fgkCarbonFiberJunctionWidth;             // Carbon Fiber Junction Width 
  static const Double_t fgkCarbonFiberJunctionEdge[2];           // Carbon Fiber Junction Edge Length  
  static const Double_t fgkCarbonFiberJunctionAngle[2];          // Carbon Fiber Junction Angle 
  static const Double_t fgkCarbonFiberJunctionToSensorSupport;   // Carbon Fiber Junction position respect to sensor
  /////////////////////////////////////////////////////////////////////////
  //Parameters for Carbon Fiber Lower Support (lengths are in mm)
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkLowerSupportToSensorZ;                   // Distance from lower support to sensor center
  static const Double_t fgkCarbonFiberLowerSupportWidth;            // Lower Support of Carbon Fiber Width
  static const Double_t fgkCarbonFiberLowerSupportLowerLenght;      // Lower Support of Carbon Fiber Length
  static const Double_t fgkCarbonFiberLowerSupportHeight;           // Lower Support of Carbon Fiber Height
  static const Double_t fgkCarbonFiberLowerSupportTransverseWidth;  // Lower Support of Carbon Fiber Transverse separation
  static const Double_t fgkCarbonFiberLowerSupportVolumeSeparation; // Distance between Lower Supports of Carbon Fiber 
  static const Double_t fgkCarbonFiberLowerSupportVolumePosition[2];// Carbon fiber lower Support Position  
  /////////////////////////////////////////////////////////////////////////
  // End Ladder Carbon Fiber Lower Junction Parameters (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkEndLadderCarbonFiberLowerJunctionLength[2];   // End Ladder Carbon Fiber Lower Up Support length 
  static const Double_t fgkEndLadderCarbonFiberUpperJunctionLength[2];   // End Ladder Carbon Fiber Lower Down Support length 
  static const Double_t fgkEndLadderMountingBlockPosition[2];            // End Ladder Mounting Block Position 
  static const Double_t fgkendladdercoolingsupportdistance[3];			 // End Ladder Cooling Support Position
  /////////////////////////////////////////////////////////////////////////
  // Cooling Tube Support (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkCoolingTubeSupportRmax;          // Cooling Tube Support Max Radius
  static const Double_t fgkCoolingTubeSupportRmin;          // Cooling Tube Support Min Radius
  static const Double_t fgkCoolingTubeSupportLength;        // Cooling Tube Support Length
  static const Double_t fgkCoolingTubeSupportHeight;        // Cooling Tube Support Height
  static const Double_t fgkCoolingTubeSupportWidth;         // Cooling Tube Support Width
  static const Double_t fgkCoolingTubeSupportSeparation;    // Cooling Tube Support Separation
  static const Double_t fgkCoolingTubeSupportToCarbonFiber; // Cooling Tube Support position respect to Carbon Fiber  
  /////////////////////////////////////////////////////////////////////////////////
  // Cooling Tube (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkCoolingTubeRmax;       // Cooling Tube Max Radius 
  static const Double_t fgkCoolingTubeRmin;       // Cooling Tube Min Radius
  static const Double_t fgkCoolingTubeLength;     // Cooling Tube Length  
  static const Double_t fgkCoolingTubeSeparation; // Cooling Tube Separation
  static const Double_t fgkMountingBlockToSensorSupport; // Distance between Mounting block and Side Sensor Support	
  /////////////////////////////////////////////////////////////////////////
  // SSD Mounting Block Parameters (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDMountingBlockLength[3];  // SSD Mounting Block Components Lengths
  static const Double_t fgkSSDMountingBlockHeight[4];  // SSD Mounting Block Components Heights
  static const Double_t fgkSSDMountingBlockWidth;      // SSD Mounting Block Width
  static const Double_t fgkSSDMountingBlockHoleTrapezoidAngle;  // SSD Mounting Block Hole Trapezoid Angle
  static const Double_t fgkSSDMountingBlockHoleTrapezoidHeight; // SSD Mounting Block Hole Trapezoid Height
  static const Double_t fgkSSDMountingBlockHoleTrapezoidUpBasis;// SSD Mounting Block Hole Trapezoid Up Basis Length
  static const Double_t fgkSSDMountingBlockHoleTubeLength[2];   // SSD Mounting Block Hole Tube Lengths   
  static const Double_t fgkSSDMountingBlockHoleTubeWidth[2];    // SSD Mounting Block Hole Tube Width   
  static const Double_t fgkSSDMountingBlockHoleRadius;          // SSD Mounting Block Hole radius  
  static const Double_t fgkSSDMountingBlockScrewHoleEdge;       // SSD Mounting Block Screw Hole Edge  
  static const Double_t fgkSSDMountingBlockScrewHoleHeight;     // SSD Mounting Block Screw Hole Height  
  static const Double_t fgkSSDMountingBlockScrewHoleRadius[2];  // SSD Mounting Block Screw Hole Radii
  /////////////////////////////////////////////////////////////////////////
  // SSD Mounting Block Clip Parameters (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkMountingBlockClipLength;          // SSD Mounting Block Clip Length  
  static const Double_t fgkMountingBlockClipThickness;       // SSD Mounting Block Clip Thickness 
  static const Double_t fgkMountingBlockClibScrewRadius;     // SSD Mounting Block Clip Radius 
  static const Double_t fgkMountingBlockClibScrewPosition;  // SSD Mounting Block Clip Screw Position
  static const Double_t fgkMountingBlockClibWidth;          // SSD Mounting Block Clip 
  /////////////////////////////////////////////////////////////////////////////////
  // SSD Mounting Block Support Parameters (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkMountingBlockSupportWidth[2]; // SSD Mounting Block Support Width
  static const Double_t fgkMountingBlockSupportDownHeight; // SSD Mounting Block Support Down Height
  static const Double_t fgkMountingBlockSupportRadius[2];  // SSD Mounting Block Support Radius
  static const Double_t fgkMountingBlockSupportUpHeight[2]; // SSD Mounting Block Support Height
  static const Double_t fgkLadderSupportHeight;            // SSD Ladder Support Width
  static const Double_t fgkLadderSupportRingLay5Position;  // SSD Ladder Support Ring Position Layer5 respect to ITS center
  static const Double_t fgkLadderSupportRingLay6Position;  // SSD Ladder Support Ring Position Layer6 respect to ITS center
  /////////////////////////////////////////////////////////////////////////////////
  // SSD End Cap Cover Plate Parameters (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkEndCapCoverPlateSmallHoleRadius; // End Cap Cover Plate Hole Small Radious 
  static const Double_t fgkEndCapCoverPlateBigHoleRadius;   // End Cap Cover Plate Hole Big Radious
  static const Double_t fgkEndCapCoverPlateThickness;       // End Cap Cover Plate Thickness
  static const Double_t fgkEndCapCoverPlateSmallHoleSeparation[3]; // End Cap Cover Plate Hole Separation
  static const Double_t fgkEndCapCoverPlateLength[6];       // End Cap Cover Plate Length
  static const Double_t fgkEndCapCoverPlateWidth[3];        // End Cap Cover Plate Width
  static const Double_t fgkEndCapCoverPlateScrewRadiusMin;  // End Cap Cover Plate Screw Radius Min
  static const Double_t fgkEndCapCoverPlateScrewRadiusMax;  // End Cap Cover Plate Screw Radius Max
  static const Double_t fgkEndCapCoverPlateClipLength;      // End Cap Cover Plate Clip Length
  static const Double_t fgkEndCapCoverPlateClipWidth;       // End Cap Cover Plate Clip Width
  static const Double_t fgkEndCapCoverPlateDownClipLength;  // End Cap Cover Plate Down Clip Length
  static const Double_t fgkEndCapCoverPlateDownClipWidth;   // End Cap Cover Plate Down Clip Width
  /////////////////////////////////////////////////////////////////////////////////
  // SSD End Cap Cooling Tube Parameters (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkEndCapCoolingTubeAxialRadius[2]; // End Cap Cooling Tube Axial Radius
  static const Double_t fgkEndCapCoolingTubeRadiusMin; // End Cap Cooling Tube Min Radius
  static const Double_t fgkEndCapCoolingTubeRadiusMax; // End Cap Cooling Tube Max Radius
  static const Double_t fgkEndCapCoolingTubeAngle[5];  // End Cap Cooling Tube Angle
  static const Double_t fgkEndCapCoolingTubeLength[5]; // End Cap Cooling Tube Length
  static const Double_t fgkEndCapCoolingTubeToCoverSide; // End Cap Cooling Tube Position respect to CoverSide
  /////////////////////////////////////////////////////////////////////////////////
  // SSD End Cap Cover Side Parameters (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkEndCapSideCoverLength[5]; // End Cap Cover Side Length
  static const Double_t fgkEndCapSideCoverWidth[7]; //  End Cap Cover Side Width
  static const Double_t fgkEndCapSideCoverThickness; // End Cap Cover Side Thickness
  /////////////////////////////////////////////////////////////////////////////////
  // SSD End Cap Cards Parameters (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////////////
  static const Double_t  fgkEndCapCardElectBoardBackLength[3]; // End Cap Card Electronic Board Back Length 
  static const Double_t  fgkEndCapCardElectBoardBackWidth[2];  // End Cap Card Electronic Board Back Width 
  static const Double_t  fgkEndCapCardElectBoardBackThickness; // End Cap Card Electronic Board Back Thickness 
  static const Double_t  fgkEndCapCardElectBoardLength;        // End Cap Card Electronic Board Length
  static const Double_t  fgkEndCapCardElectBoardLayerWidth[2]; // End Cap Card Electronic Board Layer Width
  static const Double_t  fgkEndCapCardElectBoardLayerThickness;// End Cap Card Electronic Board Layer Thickness 
  static const Double_t  fgkEndCapCardJMDConnectorThickness;   // End Cap Card JMD Connector Thickness
  static const Double_t  fgkEndCapCardJMDConnectorLength[2];   // End Cap Card JMD Connector Length
  static const Double_t  fgkEndCapCardJMDConnectorWidth[2];    // End Cap Card JMD Connector Width
  static const Double_t  fgkEndCapCardJMDConnectorToLayer;     // End Cap Card JMD Connector to Layer Distance 
  static const Double_t  fgkEndCapCardCableConnectorLength[3]; // End Cap Card Cable Connector Length
  static const Double_t  fgkEndCapCardCableConnectorWidth[2];  // End Cap Card Cable Connector Width
  static const Double_t  fgkEndCapCardCableConnectorThickness; // End Cap Card Cable Connector Thickness
  static const Double_t  fgkEndCapCardCableConnectorDistance;  // End Cap Card Cable Connector Distance
  static const Double_t  fgkEndCapCardCableConnectorToLayer;   // End Cap Card Cable Connector To Layer Distance
  static const Double_t  fgkEndCapStripConnectionLength;       // End Cap Strip Connection Length
  static const Double_t  fgkEndCapStripConnectionThickness;    // End Cap Strip Connection Thickness
  static const Double_t  fgkEndCapStripConnectionWidth;        // End Cap Strip Connection Width
  static const Double_t  fgkEndCapInterfaceCardBLength[7];     // End Cap Interface CardB Length
  static const Double_t  fgkEndCapInterfaceCardBWidth[5];      // End Cap Interface CardB Width
  static const Double_t  fgkEndCapInterfaceCardBThickness;     // End Cap Interface CardB Thickness
  static const Double_t  fgkEndCapInterfaceElectBoardCardBThickness; // End Cap Interface Elect Board CardB Thickness 
  static const Double_t  fgkEndCapInterfaceCardBJMDConnectorSeparation; // End Cap Interface CardB JMD Connector Separation
  static const Double_t  fgkEndCapStiffenerLength;             // End Cap Stiffener Length
  static const Double_t  fgkEndCapStiffenerWidth;			   // End Cap Stiffener Width
  static const Double_t  fgkEndCapStiffenerThickness;          // End Cap Stiffener Thickness
  static const Double_t  fgkEndCapEffectiveCableRadiusMin;     // End Cap Effective Cable Radius Min
  static const Double_t  fgkEndCapEffectiveCableRadiusMax;     // End Cap Effective Cable Radius Max
  /////////////////////////////////////////////////////////////////////////////////
  // SSD End Cap SupportLayer5/6 Side Parameters (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkEndCapLay5SupportLength; // End Cap Layer5 Support Length
  static const Double_t fgkEndCapLay5SupportMiddleRadius; // End Cap Layer5 Support Middle Radius
  static const Double_t fgkEndCapLay5SupportLowRadius; // End Cap Layer5 Support Low Radius
  static const Double_t fgkEndCapLay5SupportHighWidth; // End Cap Layer5 High Radius
  static const Double_t fgkEndCapLay5SupportLowWidth; // End Cap Layer5 Low Width
  static const Double_t fgkEndCapSupportLength[2]; // End Cap Layer5/6 Support Length
  static const Double_t fgkEndCapSupportMiddleRadius[2]; // End Cap Layer5/6 Support Middle Radius
  static const Double_t fgkEndCapSupportLowRadius[2]; // End Cap Layer5/6 Support Low Radius
  static const Double_t fgkEndCapSupportHighWidth; // End Cap Layer5/6 High Radius
  static const Double_t fgkEndCapSupportLowWidth[2]; // End Cap Layer5/6 Low Width  
  static const Double_t fgkEndCapSupportCenterLay5ITSPosition; // End Cap Support Center ITS Position Layer 5
  static const Double_t fgkEndCapSupportCenterLay5Position; // End Cap Support Position Respect Z Axis Origin Layer 5 
  static const Double_t fgkEndCapSupportCenterLay6ITSPosition; // End Cap Support Center ITS Position Layer 6
  static const Double_t fgkEndCapSupportCenterLay6Position; // End Cap Support Position Respect Z Axis Origin Layer 6 
  /////////////////////////////////////////////////////////////////////////////////
  // SSD End Cap Kapton Foil Parameters (lengths are in mm and angles in degrees)
  ////////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkEndCapKaptonFoilThickness; // End Cap Kapton Foil Thickness
  static const Double_t fgkEndCapKaptonFoilLength;    // End Cap Kapton Foil Length
  static const Double_t fgkEndCapKaptonFoilWidth ;    // End Cap Kapton Foil Width
  /////////////////////////////////////////////////////////////////////////
  // SSD Cone
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDLowerPConeRadius; // SSD Cone Lower Radius
  static const Double_t fgkSSDPConeAngle;        // SSD Cone Angle
  static const Double_t fgkSSDPConeZLength[2];   // SSD Cone ZLength
  static const Double_t fgkSSDPConeLittleHoleRadius; // SSD Cone Little Hole Radius
  static const Double_t fgkSSDPConeLittleHoleLength; // SSD Cone Little Hole Length
  static const Double_t fgkSSDConeMiddleRadius; // SSD Cone Middle Radius 
  static const Double_t fgkSSDPConeMiddleLength; // SSD Cone Middle Length
  static const Double_t fgkSSDPConeMiddleWidth; // SSD Cone Middle Width
  static const Double_t fgkSSDPConeUpRadius;  // SSD Cone Up Radius
  static const Double_t fgkSSDPConeUpMaxRadius; // SSD Cone Up Max Radius
  static const Double_t fgkSSDPConeUpMiddleRadius; // SSD Cone Up Middle Radius
  static const Double_t fgkSSDPConeDownRadius; // SSD Cone Down Radius
  static const Double_t fgkSSDPConeTrapezoidAngle; // SSD Cone Trapezoid Angle
  static const Double_t fgkSSDPConeTrapezoidBasis; // SSD Cone Trapezoid Basis
  static const Double_t fgkSSDPConeExternalRadius; // SSD Cone External Radius 
  static const Double_t fgkSSDPConeRadiusWidth; // SSD Cone Radius Width
  static const Double_t fgkSSDPConeLength; // SSD Cone Length
  static const Double_t fgkSSDCentralSupportLength; //SSD Central Support Length
  static const Double_t fgkSSDCentralSupportRadius; // SSD Central Support Radius 
  static const Double_t fgkSSDCentralSupportWidth; // SSD Central Support Width
  static const Double_t fgkSSDCentralAL3SupportLength; // SSD Central Support Length
  static const Double_t fgkSSDCentralAL3SupportWidth; // SSD Central Support Width
  /////////////////////////////////////////////////////////////////////////
  // SSD Cables e Patch Panel
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDCablesLay5TubeRadiusMin; // Radius Min Cable Tube Layer 5
  static const Double_t fgkSSDCablesLay6TubeRadiusMin; // Radius Min Cable Tube Layer 6
  static const Double_t fgkSSDCablesLay5RightSideHeight;  // Width Lay 5 Cables to be fixed in order to reproduce material budget
  static const Double_t fgkSSDCablesLay6RightSideHeight;  // // Width Lay 5 Cables to be fixed in order to reproduce material budget
  static const Double_t fgkSSDCableAngle; // Angle Cable   
  static const Double_t fgkSSDCablesLay5RightSideWaterHeight;  // Width Lay 5 Water Cables to be fixed in order to reproduce material budget
  static const Double_t fgkSSDCablesPatchPanel2RB26Angle[2]; // Angle Position Patch Panel RB26
  static const Double_t fgkSSDCablesPatchPanel2RB24Angle[2]; // Angle Position Patch Panel RB24
  static const Double_t fgkSSDPatchPanel2RB26ITSDistance;    // Patch Panel RB26 Position
  static const Double_t fgkSSDPatchPanel2RB24ITSDistance;   // Patch Panel RB24 Position 
  static const Double_t fgkSSDPatchPanel2RB26Radius; // Patch Panel Radius 
  static const Double_t fgkSSDPatchPanel2RB24Radius; // Patch Panel Radius
  static const Double_t fgkSSDPatchPanelHeight; // Patch Panel Height
  static const Double_t fgkSSDCableMaterialBudgetHeight; // SSD Cable Material Budget
  /////////////////////////////////////////////////////////////////////////
  // Private methods for private members generation
  /////////////////////////////////////////////////////////////////////////

  void CreateTransformationMatrices();  // Method setting the transformation matrices
  void CreateBasicObjects();			// Method creating the basic objects of ssd geometry
  void SetSSDSensor();					// Method setting the SSD Layer 5 and 6 sensors
  TList* GetCarbonFiberSupportList();	// Method generating CarbonFiberSupport
  TGeoVolume* GetCarbonFiberJunction(Double_t width); // Method generating 
										// CarbonFiberJunction
  TList* GetCarbonFiberLowerSupportList(); 
									    // Method generating CarbonFiberLowerSupport
  TGeoVolume* GetSSDSensorSupport(Double_t length, Double_t height, 
								  Double_t width, const Double_t* thickness) const; //
										// Method generating SSDSensorSupport
  TGeoVolume* GetCoolingTubeSupport(Int_t nedges); // Method generating CoolingTubeSupport 
  TList* GetSSDHybridParts();			// Method setting Hybrid Components 
  TGeoVolume* GetCoolingBlockSystem();  // Method generating Cooling Block System
  TGeoVolume* GetSSDStiffenerFlex()const;    // Method generating StiffenerFlex
  TGeoVolume* GetSSDEndFlex();			// Method generating EndFlex
  TGeoVolume* GetSSDMountingBlock();	// Method generating Mounting Block
  TGeoVolume* GetMountingBlockClip() const;   // Method generating Mounting Block Clip
  void CreateCoolingTubes();			// Create/set cooling tubes 
  TGeoVolume* GetSSDCoolingBlock(Int_t nedges); 
									    // Method generating StiffenerFlex
  void GetSSDChipCables(TGeoVolume *&cableL, TGeoVolume *&cableR, Double_t SSDChipCablesHeigth, Int_t nedges); 
										// Method setting ChipCables
  TGeoVolume* GetSSDChip() const;     // Method generating Chips
  TList* GetLadderCableSegment(Double_t ssdendladdercablelength); 
										// Method generating LadderCableSegment
  TGeoVolume* GetLadderCable(Int_t n, Double_t ssdendladdercablelength); 
										// Method generating Ladder Cable
  TGeoVolume* GetLadderCableAssembly(Int_t n, Double_t ssdendladdercablelength); 
										// Method generating Ladder Cable Assembly
  TList* GetLadderCableAssemblyList(Int_t n, Double_t ssdendladdercablelength); 
										// Method generating Ladder Cable List
  TList* GetMountingBlockSupport(Int_t nedges); // Get Mounting Block Support
  void SetLadderSupport(Int_t nedges); // It generates the ladder support
  TGeoVolume* GetEndCapCoolingTube();  // End Cap Cooling Tube
  TGeoVolume* GetEndCapCoverPlate();   // End Cap Cover Plate
  TGeoVolume* GetEndCapSideCover() const;    // End Cap Side Cover
  TGeoVolume** GetEndCapCards() const;       // End Cap Cards
  TGeoVolume** EndCapSupport();        // End Cap Support Layer 5 and Layer 6
  void SetEndCapSupportAssembly();     // EndCap Support + End Cap Layer 5 and 6
  TGeoVolume* GetEndCapEffectiveCables(Double_t radiusmin, Double_t radiusmax, 
									   Double_t width,Int_t ncables,const char* volname); // End Cap Effective HV Cables
  TGeoXtru* GetArcShape(Double_t phi, Double_t rmin, 
					    Double_t rmax, Int_t nedges, Double_t height); 
										//Auxiliary Method for Arc Shape
  TGeoArb8* GetArbShape(TVector3 const * const vertexpos[4],const Double_t* width, 
                        Double_t height,const char* shapename,Int_t isign = 1) const;
									   // Method generating an Arb shape 
  TGeoShape* GetScrewShape(const Double_t* radius,const Int_t* edgesnumber,const Double_t* section) const;// Method Generating the Screw Shape  
  TGeoShape* GetHoleShape(Double_t radius, Int_t nedges, const Double_t *section) const;// Method Generating the Hole Shape  
  TVector3* GetReflection(const TVector3* vector,const Double_t* param) const; 
										// Given an axis specified by param,
										// it gives the reflection of the point respect to the axis
  TGeoHMatrix* AddTranslationToHMatrix(TGeoHMatrix* ct,Double_t dx,Double_t dy,
                                                       Double_t dz) const;
										// add (dx,dy,dz) translation to a initial TGeoCombiTrans
  /////////////////////////////////////////////////////////////////////////
  // Private members
  /////////////////////////////////////////////////////////////////////////
  // Materials
  /////////////////////////////////////////////////////////////////////////
  TGeoMedium* fSSDChipMedium;                    // SSD Module Chip Medium
  TGeoMedium* fSSDChipGlueMedium;                // SSD Module Chip Glue Layer Medium 
  TGeoMedium* fSSDStiffenerMedium;               // SSDStiffener Medium 
  TGeoMedium* fSSDStiffenerConnectorMedium;      // SSD Stiffener Connector Medium 
  TGeoMedium* fSSDStiffener0603CapacitorMedium;  // SSD Stiffener Capacitor 0603 Medium 
  TGeoMedium* fSSDStiffener1812CapacitorMedium;  // SSD Stiffener Capacitor 1812 Medium 
  TGeoMedium* fSSDStiffenerCapacitorCapMedium;  // SSD Stiffener Capacitor Cap Medium 
  TGeoMedium* fSSDStiffenerHybridWireMedium;     // SSD Stiffener Wire Medium  
  TGeoMedium* fSSDKaptonFlexMedium;              // SSD Flex Kapton Layer Medium    
  TGeoMedium* fSSDAlTraceFlexMedium;             // SSD Flex Al Layer Medium 
  TGeoMedium* fSSDAlTraceLadderCableMedium;      // SSD Ladder Cable Al Layer Medium
  TGeoMedium* fSSDKaptonLadderCableMedium;       // SSD Ladder Cable Kapton Layer Medium
  TGeoMedium* fSSDKaptonChipCableMedium;         // SSD Chip Cables Kapton Layer Medium 
  TGeoMedium* fSSDAlTraceChipCableMedium;        // SSD Chip Cables Al Layer Medium
  TGeoMedium* fSSDAlCoolBlockMedium;             // SSD Cooling Block Al Medium
  TGeoMedium* fSSDSensorMedium;                  // SSD Sensor Medium  
  TGeoMedium* fSSDSensorSupportMedium;                  // SSD Sensor Support Medium   
  TGeoMedium* fSSDCarbonFiberMedium;             // SSD Carbon Fiber Medium 
  TGeoMedium* fSSDTubeHolderMedium;              // Cooling Tube Support Medium
  TGeoMedium* fSSDCoolingTubeWater;              // Medium for Inner Part of Cooling Tube
  TGeoMedium* fSSDCoolingTubePhynox;             // Medium for Cooling Tube 
  TGeoMedium* fSSDSupportRingAl;                 // Medium for Support Ring
  TGeoMedium* fSSDMountingBlockMedium;           // Medium for SSD Mounting Block  
  TGeoMedium* fSSDRohaCellCone;                  // Medium for SSD Ring Cone Support
  TGeoMedium* fSSDAir;							 // SSD Air
  TGeoMedium* fSSDCopper;                        // Copper for SSD Cables
  TGeoMedium* fSSDSn;                            // Tin for SSD solderings
  /////////////////////////////////////////////////////////////////////////
  Bool_t fCreateMaterials;		  // Bool variable which verifies if materials have been created
  Bool_t fTransformationMatrices; // Bool variable which verifies if matrices have been allocated
  Bool_t fBasicObjects;          // Bool variable which verifies if basic objects have been allocated
  /////////////////////////////////////////////////////////////////////////
  // Carbon Fiber Support Matrices and Objects
  ////////////////////////////////////////////
  static const Int_t fgkcarbonfibersupportnumber = 2;				   // Support Number	
  TGeoVolume* fcarbonfibersupport[fgkcarbonfibersupportnumber];		   // Support
  TGeoHMatrix* fcarbonfibersupportmatrix[fgkcarbonfibersupportnumber]; // Support Matrix
  /////////////////////////
  // Carbon Fiber Junction
  ////////////////////////
  static const Int_t fgkcarbonfiberjunctionumber = 3;  // Carbon Fiber Number
  TGeoVolume* fcarbonfiberjunction;					   // Carbon Fiber
  TGeoHMatrix* fcarbonfiberjunctionmatrix[fgkcarbonfiberjunctionumber]; // Carbon Fiber Matrix
  /////////////////////////////
  // Carbon Fiber Lower Support
  /////////////////////////////
  static const Int_t fgkcarbonfiberlowersupportnumber = 2; // Carbon Fiber Lower Support Number
  TGeoVolume* fcarbonfiberlowersupport[fgkcarbonfiberlowersupportnumber]; // Carbon Fiber Lower Support 
  TGeoTranslation* fcarbonfiberlowersupportrans[fgkcarbonfiberlowersupportnumber];// Carbon Fiber Lower Support Translation
  /////////////////////////////
  // SSD Sensor Support
  /////////////////////////////
  static const Int_t fgkvolumekind = 2; // volumekind = 0 : side ssd support
										// volumekind = 1 : central ssd support	
  static const Int_t fgkssdsensorsupportnumber = 3; // SSD Sensor Support Number
  TGeoVolume** fssdsensorsupport[fgkvolumekind];    // SSD Sensor 
  TGeoHMatrix* fssdsensorsupportmatrix[fgkssdsensorsupportnumber]; // SSD Sensor Matrix 
  /////////////////////////////////////////////////////////////
  // SSD Cooling Tube Support
  /////////////////////////////////////////////////////////////
  static const Int_t fgkcoolingtubesupportnumber = 2; // Cooling Tube Support Number
  TGeoVolume* fcoolingtubesupport;					  // Cooling Tube Support
  TGeoHMatrix* fcoolingtubesupportmatrix[fgkcoolingtubesupportnumber]; // Cooling Tube Support Matrix 
  /////////////////////////////////////////////////////////////
  // SSD Hybrid
  /////////////////////////////////////////////////////////////
  static const Int_t fgkhybridcompnumber = 3;  // Hybrid number
  TGeoVolume* fssdhybridcomponent[fgkhybridcompnumber]; // Hybrid Components
  TGeoHMatrix* fhybridmatrix;		// Hybrid Matrix
  /////////////////////////////////////////////////////////////
  // SSD Cooling Block System
  /////////////////////////////////////////////////////////////
  static const Int_t fgkcoolingblocknumber = 4; // Cooling Block Number
  TGeoVolume* fssdcoolingblocksystem;  // Cooling Block 
  TGeoHMatrix* fcoolingblocksystematrix;  // Cooling Block Matrix 
  TGeoHMatrix* fcoolingblockmatrix[fgkcoolingblocknumber];  // Cooling System Matrix
  /////////////////////////////////////////////////////////////
  // SSD Flex  
  /////////////////////////////////////////////////////////////
  static const Int_t fgkflexnumber = 2; // Flex Number 
  TGeoVolume* fssdstiffenerflex;		// Stiffener Flex
  TGeoVolume* fssdendflex;				// End flex
  TGeoHMatrix* fstiffenerflexmatrix[fgkflexnumber]; // Stiffener Flex Matrix
  TGeoHMatrix* fendflexmatrix[fgkflexnumber];       // End Flex Matrix
  /////////////////////////////////////////
  // Cooling Tube
  /////////////////////////////////////////
  TGeoHMatrix* fcoolingtubematrix[2];  // Cooling Tube Matrix
  TGeoVolume* fcoolingtube;			// Ladder Cooling Tube 
  static const Int_t fgkendladdercoolingtubenumber = 2;		// End Ladder Cooling Tube Number  	
  TGeoVolume* fendladdercoolingtube[fgkendladdercoolingtubenumber];	// End Ladder Cooling Tube
  TGeoHMatrix* fendladdercoolingtubematrix[fgkendladdercoolingtubenumber][2];  //End ladder cooling tube matrix
  /////////////////////////////////////////
  // End Ladder Components
  /////////////////////////////////////////
  TGeoVolumeAssembly* fendladdersegment[2];  // End Ladder Segment 
  TGeoHMatrix** fendladdersegmentmatrix[2];  // End Ladder Matrix
  /////////////////////////////////////////////////////////////
  // End Ladder SSD Cooling Tube Support 
  /////////////////////////////////////////////////////////////
  TGeoHMatrix*** fendladdercoolingtubesupportmatrix; //End ladder cooling tube support matrix
  ///////////////////////////////////
  // End Ladder Carbon Fiber Junction
  ///////////////////////////////////
  static const Int_t fgkendlabbercarbonfiberjunctionumber = 2; // End Ladder Carbon fiber Junction Number
  TGeoVolume** fendladdercarbonfiberjunction[fgkendlabbercarbonfiberjunctionumber]; // End Ladder Carbon fiber Junction Volumes
  static const Int_t fgkendladdercarbonfiberjunctionmatrixnumber = 3; // End Ladder Carbon fiber Junction Matrix Number
  TGeoHMatrix** fendladdercarbonfiberjunctionmatrix[fgkendladdercarbonfiberjunctionmatrixnumber]; // End Ladder Carbon fiber Junction Matrix 
  ///////////////////////////////////
  // End Ladder Carbon Fiber Support
  ///////////////////////////////////
  static const Int_t fgkendladdercarbonfibermatrixnumber = 2; // End Ladder Carbon fiber Matrix Number
  TGeoHMatrix** fendladdercarbonfibermatrix[fgkendladdercarbonfibermatrixnumber]; // End Ladder Carbon fiber Matrix 
  ///////////////////////////////////
  // End Ladder SSD Mounting Block
  ///////////////////////////////////
  static const Int_t fgkendladdermountingblocknumber = 2; // Mounting Block Number   
  TGeoVolume* fendladdermountingblock;					  // Mounting Block
  TGeoVolume* fendladdermountingblockclip;                // Mounting Block Clip
  TGeoCombiTrans* fendladdermountingblockcombitrans[fgkendladdermountingblocknumber]; // End Ladder Mounting Block CombiTrans
  TGeoHMatrix** fendladdermountingblockclipmatrix[fgkendladdermountingblocknumber]; // End Ladder Mounting Block Clip HMatrix
  ///////////////////////////////////
  // End Ladder Lower Support
  ///////////////////////////////////
  static const Int_t fgkendladderlowersuppnumber = 2; // End Ladder Lower Support Number
  TGeoTranslation* fendladderlowersupptrans[fgkendladderlowersuppnumber+1]; // End Ladder Lower Support Translations
  /////////////////////////////////////////////////////////////////////////
  // LadderCables 
  /////////////////////////////////////////////////////////////////////////
  static const Int_t fgkladdercablematrixnumber = 4; // Number of TGeoCombiTrans 
                                                     // for positioning Ladder Cable volumes  
  static const Int_t fgkladdercablesnumber = 2;		  // Number of Ladder Cables Layers
  TGeoHMatrix **fladdercablematrix[fgkladdercablesnumber]; // TGeoCombiTrans for positioning
                                                           // Ladder Cables volumes
  ///////////////////////////////////
  // Ladder Segment
  ///////////////////////////////////
  static const Int_t fgkladdersegmentnumber = 2; // Ladder Segment Kinds Number
  TGeoVolumeAssembly* fladdersegment[fgkladdersegmentnumber]; // Ladder Segment
  ///////////////////////////////////
  // Ladder 
  ///////////////////////////////////
  static const Int_t fgkladdernumber = 2;		      // Ladder Number 
  TGeoVolume* fladder[fgkladdernumber];			      //fladder[0]: ladder of Layer 5
												      //fladder[1]: ladder of Layer 6
//  TGeoVolumeAssembly* fladder[fgkladdernumber];
  TGeoHMatrix** fladdermatrix[fgkladdernumber];       // Ladder Matrix
  ///////////////////////////////////
  // SSD Sensor
  ///////////////////////////////////
  TGeoVolume* fSSDSensor5;  // Layer 5 SSD Sensor
  TGeoVolume* fSSDSensor6;  // Layer 6 SSD Sensor
  TGeoHMatrix** fssdsensormatrix[fgkladdernumber]; // SSD Sensor Matrix
  ///////////////////////////////////
  // SSD Layer
  ///////////////////////////////////
  static const Int_t fgklayernumber = 2; // Layer Number
  TGeoVolumeAssembly* fSSDLayer5;		 // SSD Layer 5
  TGeoVolumeAssembly* fSSDLayer6;	     // SSD Layer 6
  TGeoHMatrix** flayermatrix[fgklayernumber]; // Layer Transformations
  /////////////////////////////////////////////////////////////////////////
  // Mother Volume 
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume *fMotherVol;                                    // Mother volume for ITS Layer5 and Layer6   
  TGeoVolume* GetMotherVolume() const { return fMotherVol;}; // Method returning Mother Volume
  /////////////////////////////////////////////////////////////////////////
  // Ladder Support 
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* fLay5LadderSupport[2]; // Up and Down parts of Layer5 Ladder Support
  TGeoVolume* fLay6LadderSupport[2]; // Up and Down parts of Layer6 Ladder Support
  TGeoVolumeAssembly* fLay5LadderSupportRing; // Layer5 Ladder Support Ring
  TGeoVolumeAssembly* fLay6LadderSupportRing; // Layer6 Ladder Support Ring
  /////////////////////////////////////////////////////////////////////////
  // End Cap Support + End Cap Assembly
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume** fgkEndCapSupportSystem; // End Cap Support + End Cap Assembly
  /////////////////////////////////////////////////////////////////////////
  // SSD Cone
  /////////////////////////////////////////////////////////////////////////
  TGeoVolumeAssembly*fSSDCone;  // SSD Cone  
  /////////////////////////////////////////////////////////////////////////
  // Color Display 
  /////////////////////////////////////////////////////////////////////////
  Int_t fColorCarbonFiber;    //  display colors
  Int_t fColorRyton;          //  ===
  Int_t fColorPhynox;         //  ===
  Int_t fColorSilicon;        //  ===
  Int_t fColorAl;             //  ===
  Int_t fColorNiSn;           //  ===
  Int_t fColorKapton;         //  ===
  Int_t fColorPolyhamide;     //  ===
  Int_t fColorStiffener;      //  ===
  Int_t fColorEpoxy;          //  ===
  Int_t fColorWater;		  //  ===
  Int_t fColorG10;            //  ===
ClassDef(AliITSv11GeometrySSD, 5)     // ITS v11 SSD geometry
};
#endif

