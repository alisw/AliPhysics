//*************************************************************************
// class AliITSv11GeometrySSD
// Enrico Cattaruzza                                     ecattar@ts.infn.it
//*************************************************************************
#ifndef ALIITSV11GEOMETRYSSD_H
#define ALIITSV11GEOMETRYSSD_H
class TGeoVolume;
class TGeoCombiTrans;
class TGeoMedium;
class TGeoCompositeShape;
class TGeoMatrix;
class TVector3;
class TGeoArb8;
class TList;
class TGeoTube;
class TGeoHMatrix;
#include "AliITSv11Geometry.h"
class AliITSv11GeometrySSD: public AliITSv11Geometry {
public:
  AliITSv11GeometrySSD();
  AliITSv11GeometrySSD(const AliITSv11GeometrySSD &source);
  AliITSv11GeometrySSD& operator=(const AliITSv11GeometrySSD &source);
  virtual ~AliITSv11GeometrySSD(){};
  /////////////////////////////////////////////////////////////////////////
  // Public methods
  /////////////////////////////////////////////////////////////////////////
  const char*   GetSensitiveVolumeName() const {return fgkSSDSensitiveVolName;};
  // Method returning the name of the sensitive part of detector
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDSensorSupportShape(Double_t length,Double_t height,
                                       Double_t width,Double_t* thickness);
  // Method returning the name of the sensitive part of detector
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDSensorSupport(Int_t VolumeKind,Int_t n);
  // Method returning the SSD Sensor Support     
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDSensorSupportAssembly(Int_t n);
  // Method returning the SSD Sensor Support Assembly    
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDSensor();
  // Method returning the SSD Sensor     
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDChipAssembly() const;
  // Method returning the SSD Chip Assembly    
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDChipCables(Double_t ssdchipcablesheigth,char* side);
  // Method returning the SSD Chip Cables    
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDChipCablesAssembly(Double_t ssdchipcablesheigth);
  // Method returning the SSD Chip Cables Assembly   
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDStiffenerAssembly();
  // Method returning the SSD Stiffener Assembly    
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDFlex(Double_t ssdflexradius, Double_t ssdflexheigth);
  // Method returning the SSD Flex    
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDEndFlex(Double_t ssdendflexlength,Double_t ssdflexheigth);
  // Method returning the SSD End Flex    
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDFlexAssembly();
  // Method returning the SSD Flex Assembly    
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDCoolingBlock();
  // Method returning the SSD Cooling Block     
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDCoolingBlockAssembly();
  // Method returning the SSD Cooling Block Assembly    
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDModule(Int_t ichipcablesheight);
  // Method returning the SSD Module     
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetCarbonFiberJunction(Double_t width);
  // Method returning the Carbon Fiber Junction     
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetCarbonFiberJunctionAssembly();
  // Method returning the Carbon Fiber Junction Assembly    
  /////////////////////////////////////////////////////////////////////////
  TList* GetLadderCableSegment(Double_t ssdendladdercablelength);
  // Method returning a List containing pointers to Ladder Cable Volumes    
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetLadderCable(Int_t n,Double_t ssdendladdercablelength);
  // Method generating Ladder Cable Volumes Assemblies   
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetLadderCableAssembly(Int_t n,Double_t ssdendladdercablelength);
  // Method generating Ladder Cable Volumes Assembly   
  /////////////////////////////////////////////////////////////////////////
  TList* GetLadderCableAssemblyList(Int_t n,Double_t ssdendladdercablelength);
  // Method generating Ladder Cable List Assemblies   
  /////////////////////////////////////////////////////////////////////////
  TList* GetEndLadderCarbonFiberJunctionAssembly();
  // Method generating the End Ladder Carbon Fiber Junction Assembly   
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetCarbonFiberSupport();
  // Method generating the Carbon Fiber Support   
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetCarbonFiberLowerSupport(Int_t i=0,Bool_t endladder = false);
  // Method generating the Carbon Fiber Lower Support   
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetCarbonFiberAssemblySupport();
  // Method generating the Carbon Fiber Assembly Support   
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetCoolingTubeSupport();
  // Method generating the Cooling Tube Support   
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetCoolingTubeSupportAssembly();
  // Method generating the Cooling Tube Support Assembly  
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetCoolingTube() const;
  // Method generating the Cooling Tube 
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetCoolingTubeAssembly();
  // Method generating the Cooling Tube Assembly  
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetLadderSegment(Int_t iChipCablesHeight);
  // Method generating the basic Ladder Segment element which will be replicated   
  /////////////////////////////////////////////////////////////////////////
  TList* GetEndLadderSegment();
  // Method generating the Terminal Ladder Segment  
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetSSDMountingBlock();
  // Method generating the Terminal Ladder Mounting Block  
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetLadder(Int_t iLayer);
  // Method generating the Layer5 or Layer6 Ladder  
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume* GetLayer(Int_t iLayer);
  // Method generating the Layer5 or Layer6  
  /////////////////////////////////////////////////////////////////////////
  void Layer5(TGeoVolume* mothervolume);
  // Method placing the Layer5 within mother volume  
  /////////////////////////////////////////////////////////////////////////
  void Layer6(TGeoVolume* mothervolume);
  // Method placing the Layer6 within mother volume  
  /////////////////////////////////////////////////////////////////////////
  //Auxiliary methods for shapes building
  /////////////////////////////////////////////////////////////////////////
  TVector3* GetReflection(TVector3* vector,Double_t* param) const;
  // Given an axis specified by param, it gives the reflection of the point
  // respect to the axis
  /////////////////////////////////////////////////////////////////////////
  TGeoArb8* GetArbShape(TVector3* vertexpos[],Double_t* width, 
                        Double_t height,char* shapename,Int_t isign = 1) const;
  // Method generating an arb shape 
  /////////////////////////////////////////////////////////////////////////
  TGeoArb8* GetTriangleShape(TVector3* vertexpos[],Double_t height, char* shapename) const;
  // Method generating a triangle shape 
  /////////////////////////////////////////////////////////////////////////
  TGeoArb8* GetTrapezoidShape(TVector3* vertexpos[],Double_t* width,    
                              Double_t height,char* shapename) const;
  // Method generating a trapezoid shape 
  /////////////////////////////////////////////////////////////////////////
  TGeoCombiTrans* AddTranslationToCombiTrans(TGeoCombiTrans* ct,Double_t dx,
                                             Double_t dy,Double_t dz) const;
  // add (dx,dy,dz) translation to a initial TGeoCombiTrans
  /////////////////////////////////////////////////////////////////////////
  //Auxiliary methods for material building
  /////////////////////////////////////////////////////////////////////////
  TGeoMedium* GetMedium(const char* mediumName);
  // It returns the Medium
  /////////////////////////////////////////////////////////////////////////
  void CreateMaterials();
  // it creates the matherials    
private:
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
  TGeoMedium* fSSDChipMedium;                    // SSD Module Chip Medium
  TGeoMedium* fSSDChipGlueMedium;                // SSD Module Chip Glue Layer Medium 
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
  static const Double_t fgkSSDCapacitor1812Length; // SSD Stiffener Capacitor 1812 Length 
  static const Double_t fgkSSDCapacitor1812Width;  // SSD Stiffener Capacitor 1812 Width
  static const Double_t fgkSSDCapacitor1812Height; // SSD Stiffener Capacitor 1812 Height
  static const Double_t fgkSSDWireLength;          // SSD Stiffener Wire Length
  static const Double_t fgkSSDWireRadius;          // SSD Stiffener Wire Radius
  static const Double_t fgkSSDConnectorPosition[2];// SSD Connector Position respect to Stiffener
  static const Double_t fgkSSDConnectorSeparation; // SSD Connector separation distance
  static const Double_t fgkSSDConnectorLength;     // SSD Stiffener Connector Length
  static const Double_t fgkSSDConnectorWidth;      // SSD Stiffener Connector Width
  static const Double_t fgkSSDConnectorHeight;     // SSD Stiffener Connector Height
  TGeoMedium* fSSDStiffenerMedium;               // SSDStiffener Medium 
  TGeoMedium* fSSDStiffenerConnectorMedium;      // SSD Stiffener Connector Medium 
  TGeoMedium* fSSDStiffener0603CapacitorMedium;  // SSD Stiffener Capacitor 0603 Medium 
  TGeoMedium* fSSDStiffener1812CapacitorMedium;  // SSD Stiffener Capacitor 1812 Medium 
  TGeoMedium* fSSDStiffenerHybridWireMedium;     // SSD Stiffener Wire Medium  
  /////////////////////////////////////////////////////////////////////////
  // Flex
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDFlexFullLength;      // SSD Flex Full Length
  static const Double_t fgkSSDFlexLength[4];       // SSD Flex Components Length
  static const Double_t fgkSSDFlexWidth[2];        // SSD Flex Components Width
  static const Double_t fgkSSDFlexHeight[2];       // SSD Flex Layers Heigth
  static const Double_t fgkSSDFlexAngle;           // SSD Flex Angle 
  static const Double_t fgkSSDFlexHoleLength;      // SSD Flex Hole Length
  static const Double_t fgkSSDFlexHoleWidth;       // SSD Flex Hole Width
  static const Double_t fgkSSDEndFlexCompLength[6];// SSD End-Flex Components Length
  static const Double_t fgkSSDEndFlexCompWidth[3]; // SSD End-Flex Components Width
  TGeoMedium* fSSDKaptonFlexMedium;              // SSD Flex Kapton Layer Medium    
  TGeoMedium* fSSDAlTraceFlexMedium;             // SSD Flex Al Layer Medium 
  /////////////////////////////////////////////////////////////////////////////////
  // SSD Ladder Cable 
  /////////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDLadderCableWidth;    // SSD Ladder Cable Width
  TGeoMedium* fSSDAlTraceLadderCableMedium;      // SSD Ladder Cable Al Layer Medium
  TGeoMedium* fSSDKaptonLadderCableMedium;       // SSD Ladder Cable Kapton Layer Medium
  /////////////////////////////////////////////////////////////////////////
  // SSD Module Components 
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDModuleStiffenerPosition[2]; // SSD Module Stiffener position 
                                                          // respect to Sensor Edge
  static const Double_t fgkSSDModuleSensorSupportDistance;// SSD Module Sensor Support Position 
                                                          // respect to Sensor Edge 
  static const Double_t fgkSSDModuleCoolingBlockToSensor; // SSD Cooling Block Position 
                                                          // respect to sensor
  static const Int_t fgkSSDModuleCombiTransNumber = 7;    // Number of TGeoCombiTrans 
                                                          // for positioning volumes in SSD Module 
  void SetSSDModuleCombiTransMatrix(Double_t);            // Method for generating TGeoCombiTrans for 
                                                          // for volume positioning in SSD Module
  TGeoCombiTrans *fSSDModuleCombiTransMatrix[fgkSSDModuleCombiTransNumber]; // TGeoCombiTrans roto-trans
                                                          // transformations for volume positioning
  /////////////////////////////////////////////////////////////////////////
  // Chip Cables
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDChipCablesLength[2]; // SSD Chip Cables Components Length
  static const Double_t fgkSSDChipCablesHeight[4]; // SSD Chip Cables Components Height   
  static const Double_t fgkSSDChipCablesWidth[3];  // SSD Chip Cables Components Width
  TGeoMedium* fSSDKaptonChipCableMedium;         // SSD Chip Cables Kapton Layer Medium 
  TGeoMedium* fSSDAlTraceChipCableMedium;        // SSD Chip Cables Al Layer Medium
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
  TGeoMedium* fSSDAlCoolBlockMedium;                  // SSD Cooling Block Al Medium
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
  TGeoMedium* fSSDSensorMedium;                      // SSD Sensor Medium  
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
  void SetSSDSensorSupportCombiTransMatrix();                 // Method for generating TGeoCombiTrans for 
                                                              // volume positioning in Sensor Support Assembly   
  TGeoCombiTrans *fSSDSensorSupportCombiTransMatrix[fgkSSDSensorSupportCombiTransNumber]; // TGeoCombiTrans roto-trans
                                                          // transformations for volume positioning in Sensor Support Assembly
  TGeoMedium* fSSDSensorSupportMedium;                  // SSD Sensor Support Medium   
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
  void SetCarbonFiberAssemblyCombiTransMatrix();                 // Method for generating TGeoCombiTrans for 
                                                                 // volume positioning in Carbon Fiber Assembly
  TGeoCombiTrans *fCarbonFiberAssemblyCombiTransMatrix[fgkCarbonFiberAssemblyCombiTransNumber]; // TGeoCombiTrans roto-trans
                                                                 // transformations for volume positioning in Carbon Fiber Assembly
  TGeoMedium* fSSDCarbonFiberMedium;                           // SSD Carbon Fiber Medium 
  //////////////////////////////////////////////////////////////////////////////
  // Carbon Fiber Junction Parameters
  //////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkCarbonFiberJunctionLength;            // Carbon Fiber Junction Length             
  static const Double_t fgkCarbonFiberJunctionWidth;             // Carbon Fiber Junction Width 
  static const Double_t fgkCarbonFiberJunctionEdge[2];           // Carbon Fiber Junction Edge Length  
  static const Double_t fgkCarbonFiberJunctionAngle[2];          // Carbon Fiber Junction Angle 
  static const Double_t fgkCarbonFiberJunctionToSensorSupport;   // Carbon Fiber Junction position respect to sensor
  static const Int_t fgkCarbonFiberJunctionCombiTransNumber = 3; // Number of TGeoCombiTrans 
                                                                 // for positioning volumes in Carbon Fiber Junction Assembly 
  void SetCarbonFiberJunctionCombiTransMatrix();                 // Method for generating TGeoCombiTrans for 
                                                                 // volume positioning in Carbon Fiber Junction Assembly
  TGeoCombiTrans *fCarbonFiberJunctionCombiTransMatrix[fgkCarbonFiberJunctionCombiTransNumber];// TGeoCombiTrans roto-trans
                                                                 // transformations for volume positioning in Carbon Fiber Junction Assembly
  /////////////////////////////////////////////////////////////////////////
  //Parameters for Carbon Fiber Lower Support (lengths are in mm)
  /////////////////////////////////////////////////////////////////////////
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
  static const Int_t fgkEndLadderCarbonFiberJunctionCombiTransNumber = 3;// Number of TGeoCombiTrans 
                                                                         // for positioning volumes in End LadderCarbon Fiber Junction Assembly 
  void SetEndLadderCarbonFiberJunctionCombiTransMatrix(Int_t);           // Method for generating TGeoCombiTrans for 
                                                                         // volume positioning in End Ladder Carbon Fiber Junction Assembly
  TGeoCombiTrans *fEndLadderCarbonFiberJunctionCombiTransMatrix[fgkEndLadderCarbonFiberJunctionCombiTransNumber]; // TGeoCombiTrans roto-trans
                                                                         // transformations for volume positioning in End Ladder 
                                                                         // Carbon Fiber Junction Assembly
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
  static const Int_t fgkCoolingTubeSupportCombiTransNumber = 2; // Number of TGeoCombiTrans 
                                                            // for positioning volumes in Cooling Tube Support Assembly
  void SetCoolingTubeSupportCombiTransMatrix();             // Method for generating TGeoCombiTrans for 
                                                            // volume positioning in Cooling Tube Support Assembly
  TGeoCombiTrans *fCoolingTubeSupportCombiTransMatrix[fgkCoolingTubeSupportCombiTransNumber]; // TGeoCombiTrans roto-trans
                                                            // transformations for volume positioning in  
                                                            // Cooling Tube Support Assembly
  TGeoMedium* fSSDTubeHolderMedium;                       // Cooling Tube Support Medium
  /////////////////////////////////////////////////////////////////////////////////
  // Cooling Tube (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkCoolingTubeRmax;       // Cooling Tube Max Radius 
  static const Double_t fgkCoolingTubeRmin;       // Cooling Tube Min Radius
  static const Double_t fgkCoolingTubeLength;     // Cooling Tube Length  
  static const Double_t fgkCoolingTubeSeparation; // Cooling Tube Separation
  static const Int_t    fgkCoolingTubeCombiTransNumber = 2; // Number of TGeoCombiTrans 
                                                            // for positioning volumes in Cooling Tube Assembly
  void SetCoolingTubeCombiTransMatrix();                    // Method for generating TGeoCombiTrans for 
                                                            // volume positioning in Cooling Tube Assembly
  TGeoCombiTrans *fCoolingTubeTransMatrix[fgkCoolingTubeCombiTransNumber];// TGeoCombiTrans roto-trans
                                                            // transformations for volume positioning Cooling Tube Assembly
  TGeoMedium* fSSDCoolingTubeWater;           // Medium for Inner Part of Cooling Tube
  TGeoMedium* fSSDCoolingTubePhynox;          // Medium for Cooling Tube 
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
  static const Double_t fgkSSDMountingBlockScrewHoleHeigth;     // SSD Mounting Block Screw Hole Height  
  static const Double_t fgkSSDMountingBlockScrewHoleRadius[2];  // SSD Mounting Block Screw Hole Radii
  TGeoMedium* fSSDMountingBlockMedium;                        // Medium for SSD Mounting Block  
  /////////////////////////////////////////////////////////////////////////
  // LadderSegment 
  /////////////////////////////////////////////////////////////////////////
  static const Int_t fgkLadderSegmentCombiTransNumber = 5; // Number of TGeoCombiTrans 
                                                           // for positioning volumes in Ladder Segment    
  void SetLadderSegmentCombiTransMatrix();                 // Method for generating TGeoCombiTrans for 
                                                           // volume positioning in Ladder Segment 
  TGeoCombiTrans *fLadderSegmentCombiTransMatrix[fgkLadderSegmentCombiTransNumber]; // TGeoCombiTrans roto-trans
                                                            // transformations for volume positioning in  
                                                            // in Ladder Segment
  /////////////////////////////////////////////////////////////////////////
  // End LadderSegment 
  /////////////////////////////////////////////////////////////////////////
  static const Int_t fgkEndLadderSegmentCombiTransNumber = 4; // Number of TGeoCombiTrans 
                                                           // for positioning volumes in End Ladder Segment    
  void SetEndLadderSegmentCombiTransMatrix(Int_t);         // Method for generating TGeoCombiTrans for 
                                                           // volume positioning in End Ladder Segment  
  TGeoCombiTrans *fEndLadderSegmentCombiTransMatrix[fgkEndLadderSegmentCombiTransNumber]; // TGeoCombiTrans roto-trans
                                                            // transformations for volume positioning in  
                                                            // in End Ladder Segment
  /////////////////////////////////////////////////////////////////////////
  // LadderCables 
  /////////////////////////////////////////////////////////////////////////
  static const Int_t fgkLadderCableCombiTransNumber = 4; // Number of TGeoCombiTrans 
                                                         // for positioning Ladder Cable volumes  
  void SetLadderCableCombiTransMatrix(Int_t iLayer);     // Method for generating TGeoCombiTrans 
                                                         // for positioning Ladder Cables volumes
  TGeoCombiTrans *fLadderCableCombiTransMatrix[fgkLadderCableCombiTransNumber]; // TGeoCombiTrans for positioning
                                                         // Ladder Cables volumes
  /////////////////////////////////////////////////////////////////////////
  // Mother Volume 
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume *fMotherVol;                                    // Mother volume for ITS Layer5 and Layer6   
  TGeoVolume* GetMotherVolume() const { return fMotherVol;}; // Method returning Mother Volume
  /////////////////////////////////////////////////////////////////////////
  // Color Display 
  /////////////////////////////////////////////////////////////////////////
  Int_t fColorCarbonFiber;    //  display colors
  Int_t fColorRyton;          //  ===
  Int_t fColorPhynox;         //  ===
  Int_t fColorSilicon;        //  ===
  Int_t fColorAl;             //  ===
  Int_t fColorKapton;         //  ===
  Int_t fColorPolyhamide;     //  ===
  Int_t fColorStiffener;      //  ===
  Int_t fColorEpoxy;          //  ===
  ClassDef(AliITSv11GeometrySSD, 2)     // ITS v11 SSD geometry
};
#endif
