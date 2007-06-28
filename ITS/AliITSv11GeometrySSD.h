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
class TGeoTranslation;
class TGeoXtru;
class TGeoVolumeAssembly;
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
  TGeoMedium* GetMedium(const char* mediumName);
  // It returns the Medium
  const char*   GetSenstiveVolumeName5() const {return fgSDDsensitiveVolName5;};
  // it returns the Sensitive Volume of Layer 5
  const char*   GetSenstiveVolumeName6() const {return fgSDDsensitiveVolName6;};
  // it returns the Sensitive Volume of Layer 6
  TGeoVolumeAssembly* GetLadderSegment(Int_t i){return fladdersegment[i];}; // Get Ladder Segment
  TGeoVolumeAssembly* GetEndLadderSegment(Int_t i){return fendladdersegment[i];}; // Get End Ladder Segment 
  TGeoVolume* GetLadder(Int_t i) {return fladder[i];}; // Get Ladder
  TGeoVolume* GetLayer(Int_t i)const {return i==5? fSSDLayer5 : fSSDLayer6;}; // Get Layer
  void SetLadderSegment();				// Set Ladder Elementary Segment 
  void SetEndLadderSegment();			// Set End Ladder Segment
  void SetLadder();						// Set Ladder
  void SetLayer();						// Set Layer
  void Layer5(TGeoVolume* moth);        // Setting Layer 5 into mother volume
  void Layer6(TGeoVolume* moth);        // Setting Layer 6 into mother volume
private:
  /////////////////////////////////////////////////////////////////////////////////
  // Names of the Sensitive Volumes of Layer 5 and Layer 6
  /////////////////////////////////////////////////////////////////////////////////
  static const char* fgSDDsensitiveVolName5;       // sens. vol. name for lay. 5
  static const char* fgSDDsensitiveVolName6;       // sens. vol. name for lay. 6
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
  static const Double_t fgkSSDConnectorAlHeight;     // SSD Stiffener Connector Al Height
  static const Double_t fgkSSDConnectorNiHeight;     // SSD Stiffener Connector Ni Height
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
  /////////////////////////////////////////////////////////////////////////////////
  // SSD Ladder Cable 
  /////////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDLadderCableWidth;    // SSD Ladder Cable Width
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
  /////////////////////////////////////////////////////////////////////////
  // Private methods for private members generation
  /////////////////////////////////////////////////////////////////////////
  void CreateMaterials();				// Method setting the materials 
  void CreateTransformationMatrices();  // Method setting the transformation matrices
  void CreateBasicObjects();			// Method creating the basic objects of ssd geometry
  void SetSSDSensor();					// Method setting the SSD Layer 5 and 6 sensors
  TList* GetCarbonFiberSupportList();	// Method generating CarbonFiberSupport
  TGeoVolume* GetCarbonFiberJunction(Double_t width); // Method generating 
										// CarbonFiberJunction
  TList* GetCarbonFiberLowerSupportList(); 
									    // Method generating CarbonFiberLowerSupport
  TGeoVolume* GetSSDSensorSupport(Double_t length, Double_t height, 
								  Double_t width, Double_t* thickness) const; //
										// Method generating SSDSensorSupport
  TGeoVolume* GetCoolingTubeSupport(Int_t nedges); // Method generating CoolingTubeSupport 
  TList* GetSSDHybridParts();			// Method setting Hybrid Components 
  TGeoVolume* GetCoolingBlockSystem();  // Method generating Cooling Block System
  TGeoVolume* GetSSDStiffenerFlex()const;    // Method generating StiffenerFlex
  TGeoVolume* GetSSDEndFlex();			// Method generating EndFlex
  TGeoVolume* GetSSDMountingBlock();	// Method generating Mounting Block
  TList* GetCoolingTubeList()const;			// Method generating list of Tubes
  TGeoVolume* GetSSDCoolingBlock(Int_t nedges); 
									    // Method generating StiffenerFlex
  TGeoVolume* GetSSDChipCables(Double_t SSDChipCablesHeigth, Int_t nedges); 
										// Method setting ChipCables
  TList* GetSSDChipSystem();			// Method setting Chip System
  TGeoVolume* GetSSDChips() const;     // Method generating Chips
  TList* GetLadderCableSegment(Double_t ssdendladdercablelength); 
										// Method generating LadderCableSegment
  TGeoVolume* GetLadderCable(Int_t n, Double_t ssdendladdercablelength); 
										// Method generating Ladder Cable
  TGeoVolume* GetLadderCableAssembly(Int_t n, Double_t ssdendladdercablelength); 
										// Method generating Ladder Cable Assembly
  TList* GetLadderCableAssemblyList(Int_t n, Double_t ssdendladdercablelength); 
										// Method generating Ladder Cable List
  TGeoXtru* GetArcShape(Double_t phi, Double_t rmin, 
					    Double_t rmax, Int_t nedges, Double_t height); 
										//Auxiliary Method for Arc Shape
  TGeoArb8* GetArbShape(TVector3* vertexpos[],Double_t* width, 
                        Double_t height,char* shapename,Int_t isign = 1) const;
									   // Method generating an Arb shape 
  TVector3* GetReflection(TVector3* vector,Double_t* param) const; 
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
  TGeoMedium* fSSDMountingBlockMedium;           // Medium for SSD Mounting Block  
  TGeoMedium* fSSDAir;							 // SSD Air
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
  static const Int_t fgkcoolingtubenumber = 3;				// Coling Tube Number
  TGeoHMatrix** fcoolingtubematrix[fgkcoolingtubenumber+1];  // Cooling Tube Matrix
  TGeoVolume* fcoolingtube[fgkcoolingtubenumber];	     // Cooling Tube
  /////////////////////////////////////////
  // End Ladder Components
  /////////////////////////////////////////
  TGeoVolumeAssembly* fendladdersegment[2];  // End Ladder Segment 
  TGeoHMatrix** fendladdersegmentmatrix[2];  // End Ladder Matrix
  ///////////////////////////////////
  // End Ladder Carbon Fiber Junction
  ///////////////////////////////////
  static const Int_t fgkendlabbercarbonfiberjunctionumber = 2; // End Ladder Carbon fiber Junction Number
  TGeoVolume** fendladdercarbonfiberjunction[fgkendlabbercarbonfiberjunctionumber]; // End Ladder Carbon fiber Junction Volumes
  static const Int_t fgkendladdercabonfiberjunctionmatrixnumber = 3; // End Ladder Carbon fiber Junction Matrix Number
  TGeoHMatrix** fendladdercarbonfiberjunctionmatrix[fgkendlabbercarbonfiberjunctionumber]; // End Ladder Carbon fiber Junction Matrix 
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
  TGeoTranslation* fendladdermountingblocktrans[fgkendladdermountingblocknumber]; // Mounting Block Translation
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
  TGeoVolume* fSSDLayer5;				 // SSD Layer 5
  TGeoVolume* fSSDLayer6;	             // SSD Layer 6
  TGeoHMatrix** flayermatrix[fgklayernumber]; // Layer Transformations
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
  Int_t fColorWater;		  //  ===
  Int_t fColorG10;            //  ===
  ClassDef(AliITSv11GeometrySSD, 2)     // ITS v11 SSD geometry
};
#endif
