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
class TVector3;
class TGeoArb8;
class TList;
class TGeoTube;
class TGeoHMatrix;
class AliITSv11GeometrySSD{
public:
  AliITSv11GeometrySSD();
  virtual ~AliITSv11GeometrySSD(){};
  /////////////////////////////////////////////////////////////////////////
  // Public methods
  /////////////////////////////////////////////////////////////////////////
  TGeoMedium *GetMedium(char*); // to be interfaced with AliRoot
  const char*   GetSensitiveVolumeName() const {return fgkSSDSensitiveVolName;};
  TGeoVolume* GetSSDSensorSupportShape(Double_t,Double_t,Double_t,Double_t*);
  TGeoVolume* GetSSDSensorSupport(Int_t, Int_t);
  TGeoVolume* GetSSDSensorSupportAssembly(Int_t);
  TGeoVolume* GetSSDSensor();
  TGeoVolume* GetSSDChipAssembly();
  TGeoVolume* GetSSDChipCables(Double_t,char*);
  TGeoVolume* GetSSDChipCablesAssembly(Double_t);
  TGeoVolume* GetSSDStiffenerAssembly();
  TGeoVolume* GetSSDFlex(Double_t,Double_t);
  TGeoVolume* GetSSDEndFlex(Double_t,Double_t);
  TGeoVolume* GetSSDFlexAssembly();
  TGeoVolume* GetSSDCoolingBlock();
  TGeoVolume* GetSSDCoolingBlockAssembly();
  TGeoVolume* GetSSDModule(Int_t);
  TGeoVolume* GetCarbonFiberJunction(Double_t);
  TGeoVolume* GetCarbonFiberJunctionAssembly();
  TList* GetLadderCableSegment(Double_t);
  TGeoVolume* GetLadderCable(Int_t n,Double_t);
  TGeoVolume* GetLadderCableAssembly(Int_t n,Double_t);
  TList* GetLadderCableAssemblyList(Int_t n,Double_t);
  TList* GetEndLadderCarbonFiberJunctionAssembly();
  TGeoVolume* GetCarbonFiberSupport();
  TGeoVolume* GetCarbonFiberLowerSupport(Int_t i=0, Bool_t EndLadder = false);
  TGeoVolume* GetCarbonFiberAssemblySupport();
  TGeoVolume* GetCoolingTubeSupport();
  TGeoVolume* GetCoolingTubeSupportAssembly();
  TGeoVolume* GetCoolingTube();
  TGeoVolume* GetCoolingTubeAssembly();
  TGeoVolume* GetLadderSegment(Int_t);
  TList* GetEndLadderSegment();
  TGeoVolume* GetSSDMountingBlock();
  TGeoVolume* GetLadder(Int_t);
  TGeoVolume* GetLayer(Int_t);
  void Layer5(TGeoVolume*);
  void Layer6(TGeoVolume*);
  /////////////////////////////////////////////////////////////////////////
  //Auxiliary methods for shapes building
  /////////////////////////////////////////////////////////////////////////
  TVector3* GetReflection(TVector3*,Double_t*);
  TGeoArb8* GetArbShape(TVector3* [],Double_t*,Double_t,char*,Int_t isign = 1);
  TGeoArb8* GetTriangleShape(TVector3* [],Double_t,char*);
  TGeoArb8* GetTrapezoidShape(TVector3* [],Double_t*,Double_t,char*);
  TGeoCombiTrans* AddTranslationToCombiTrans(TGeoCombiTrans*,Double_t,
							 Double_t,Double_t) const;
  /////////////////////////////////////////////////////////////////////////
  //Auxiliary methods for material building
  /////////////////////////////////////////////////////////////////////////
  TGeoMedium* AliITSv11GeometrySSD::GetMedium(const char* mediumName);
  void CreateMaterials();
private:
  /////////////////////////////////////////////////////////////////////////
  // Layer5 (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDLay5LadderLength;
  static const Int_t fgkSSDLay5SensorsNumber;
  static const Int_t fgkSSDLay5LadderNumber;
  static const Double_t fgkSSDLay5RadiusMin;
  static const Double_t fgkSSDLay5RadiusMax;
  static const Double_t fgkLay5CenterITSPosition;
  /////////////////////////////////////////////////////////////////////////
  // Layer6 (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDLay6LadderLength; //It is not used explicitely in the simulation
  static const Int_t fgkSSDLay6SensorsNumber;
  static const Int_t fgkSSDLay6LadderNumber;
  static const Double_t fgkSSDLay6RadiusMin;
  static const Double_t fgkSSDLay6RadiusMax;
  static const Double_t fgkLay6CenterITSPosition;
  /////////////////////////////////////////////////////////////////////////
  // SSD Chips and Hybrid
  /////////////////////////////////////////////////////////////////////////
  static const Int_t fgkSSDChipNumber;
  static const Double_t fgkSSDChipLength; 
  static const Double_t fgkSSDChipWidth;
  static const Double_t fgkSSDChipHeight;
  static const Double_t fgkSSDChipSeparationLength;
  static const Double_t fgkSSDChipGlueLength; 
  static const Double_t fgkSSDChipGlueWidth; 
  static const Double_t fgkSSDChipGlueHeight; 
  TGeoMedium* fgkSSDChipMedium;
  TGeoMedium* fgkSSDChipGlueMedium;
  /////////////////////////////////////////////////////////////////////////
  // Stiffener
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDStiffenerLength;
  static const Double_t fgkSSDStiffenerWidth;
  static const Double_t fgkSSDStiffenerHeight;
  static const Double_t fgkSSDStiffenerToChipDist;
  static const Double_t fgkSSDCapacitor0603Length;
  static const Double_t fgkSSDCapacitor0603Width;
  static const Double_t fgkSSDCapacitor0603Height;
  static const Double_t fgkSSDCapacitor1812Length;
  static const Double_t fgkSSDCapacitor1812Width;
  static const Double_t fgkSSDCapacitor1812Height;
  static const Double_t fgkSSDWireLength;
  static const Double_t fgkSSDWireRadius;
  static const Double_t fgkSSDConnectorPosition[2];
  static const Double_t fgkSSDConnectorSeparation;
  static const Double_t fgkSSDConnectorLength;
  static const Double_t fgkSSDConnectorWidth;
  static const Double_t fgkSSDConnectorHeight;
  TGeoMedium* fgkSSDStiffenerMedium;
  TGeoMedium* fgkSSDStiffenerConnectorMedium;
  TGeoMedium* fgkSSDStiffener0603CapacitorMedium;
  TGeoMedium* fgkSSDStiffener1812CapacitorMedium;
  TGeoMedium* fgkSSDStiffenerHybridWireMedium;
  /////////////////////////////////////////////////////////////////////////
  // Flex
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDFlexFullLength;
  static const Double_t fgkSSDFlexLength[4];
  static const Double_t fgkSSDFlexWidth[2];
  static const Double_t fgkSSDFlexHeight[2];
  static const Double_t fgkSSDFlexAngle;
  static const Double_t fgkSSDFlexHoleLength;
  static const Double_t fgkSSDFlexHoleWidth;
  static const Double_t fgkSSDEndFlexCompLength[6];
  static const Double_t fgkSSDEndFlexCompWidth[3];
  TGeoMedium* fgkSSDKaptonFlexMedium;
  TGeoMedium* fgkSSDAlTraceFlexMedium;
  /////////////////////////////////////////////////////////////////////////////////
  // SSD Ladder Cable 
  /////////////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDEndLadderCableLength;
  static const Double_t fgkSSDLadderCableWidth;
  TGeoMedium* fgkSSDAlTraceLadderCableMedium;
  TGeoMedium* fgkSSDKaptonLadderCableMedium;
  /////////////////////////////////////////////////////////////////////////
  // SSD Module
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDModuleStiffenerPosition[2];
  static const Double_t fgkSSDModuleSensorSupportDistance;
  static const Double_t fgkSSDModuleCoolingBlockToSensor;
  static const Int_t fgkSSDModuleCombiTransNumber = 7;
  void SetSSDModuleCombiTransMatrix(Double_t);
  TGeoCombiTrans *SSDModuleCombiTransMatrix[fgkSSDModuleCombiTransNumber];
  /////////////////////////////////////////////////////////////////////////
  // Chip Cables
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDChipCablesLength[2];
  static const Double_t fgkSSDChipCablesHeight[4];
  static const Double_t fgkSSDChipCablesWidth[3];
  TGeoMedium* fgkSSDKaptonChipCableMedium;
  TGeoMedium* fgkSSDAlTraceChipCableMedium;
  /////////////////////////////////////////////////////////////////////////
  // Cooling Block
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDCoolingBlockLength;
  static const Double_t fgkSSDCoolingBlockWidth;
  static const Double_t fgkSSDCoolingBlockHeight[3];
  static const Double_t fgkSSDCoolingBlockHoleRadius[2];
  static const Double_t fgkSSDCoolingBlockHoleLength[2];
  static const Double_t fgkSSDCoolingBlockHoleCenter;
  static const Double_t fgkSSDCoolingBlockHoleHeight;
  TGeoMedium* fgkSSDAlCoolBlockMedium;
  /////////////////////////////////////////////////////////////////////////
  // SSD Sensor 
  /////////////////////////////////////////////////////////////////////////
  static const char* fgkSSDSensitiveVolName;
  static const Double_t fgkSSDSensorLength;
  static const Double_t fgkSSDSensorHeight;
  static const Double_t fgkSSDSensorWidth;
  static const Double_t fgkSSDSensorOverlap;
  static const Double_t fgkSSDSensorInsensitiveLength;
  static const Double_t fgkSSDSensorInsensitiveWidth;
  TGeoMedium* fgkSSDSensorMedium;
  /////////////////////////////////////////////////////////////////////////
  // SSD Sensor Support 
  /////////////////////////////////////////////////////////////////////////
  static const Int_t fgkSSDSensorSupportCombiTransNumber = 3;
  static const Double_t fgkSSDSensorSideSupportLength;
  static const Double_t fgkSSDSensorSideSupportWidth;
  static const Double_t fgkSSDSensorSideSupportHeight[2];
  static const Double_t fgkSSDSensorSideSupportThickness[2];
  static const Double_t fgkSSDSensorSideSupportPosition;
  static const Double_t fgkSSDSensorCenterSupportLength;
  static const Double_t fgkSSDSensorCenterSupportWidth;
  static const Double_t fgkSSDSensorCenterSupportHeight[2];
  static const Double_t fgkSSDSensorCenterSupportThickness[2];
  static const Double_t fgkSSDSensorCenterSupportPosition;
  void SetSSDSensorSupportCombiTransMatrix();
  TGeoCombiTrans *SSDSensorSupportCombiTransMatrix[fgkSSDSensorSupportCombiTransNumber];
  TGeoMedium* fgkSSDSensorSupportMedium;
  /////////////////////////////////////////////////////////////////////////
  //Parameters for Carbon Fiber 
  /////////////////////////////////////////////////////////////////////////
  static const Int_t fgkCarbonFiberAssemblyCombiTransNumber = 3;
  static const Double_t fgkCarbonFiberTriangleLength;
  static const Double_t fgkCarbonFiberTriangleAngle;
  static const Double_t fgkCarbonFiberSupportTopEdgeDist[2];
  static const Double_t fgkCarbonFiberSupportEdgeLength;
  static const Double_t fgkCarbonFiberSupportWidth;
  static const Double_t fgkCarbonFiberSupportXAxisLength;
  static const Double_t fgkCarbonFiberSupportYAxisLength;
  void SetCarbonFiberAssemblyCombiTransMatrix();
  TGeoCombiTrans *CarbonFiberAssemblyCombiTransMatrix[fgkCarbonFiberAssemblyCombiTransNumber];
  TGeoMedium* fgkSSDCarbonFiberMedium;
  //////////////////////////////////////////////////////////////////////////////
  // Carbon Fiber Junction Parameters
  //////////////////////////////////////////////////////////////////////////////
  static const Int_t fgkCarbonFiberJunctionCombiTransNumber = 3;
  static const Double_t fgkCarbonFiberJunctionLength;
  static const Double_t fgkCarbonFiberJunctionWidth;
  static const Double_t fgkCarbonFiberJunctionEdge[2];
  static const Double_t fgkCarbonFiberJunctionAngle[2];
  static const Double_t fgkCarbonFiberJunctionToSensorSupport;
  void SetCarbonFiberJunctionCombiTransMatrix();
  TGeoCombiTrans *CarbonFiberJunctionCombiTransMatrix[fgkCarbonFiberJunctionCombiTransNumber];
  /////////////////////////////////////////////////////////////////////////
  //Parameters for Carbon Fiber Lower Support (lengths are in mm)
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkCarbonFiberLowerSupportWidth;
  static const Double_t fgkCarbonFiberLowerSupportLowerLenght;
  static const Double_t fgkCarbonFiberLowerSupportHeight;
  static const Double_t fgkCarbonFiberLowerSupportTransverseWidth;
  static const Double_t fgkCarbonFiberLowerSupportVolumeSeparation;
  static const Double_t fgkCarbonFiberLowerSupportVolumePosition[2];
  /////////////////////////////////////////////////////////////////////////
  // End Ladder Carbon Fiber Junction Parameters (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////
  static const Int_t fgkEndLadderCarbonFiberJunctionCombiTransNumber = 3;
  static const Double_t fgkEndLadderCarbonFiberLowerJunctionLength[2]; 
  static const Double_t fgkEndLadderCarbonFiberUpperJunctionLength[2]; 
  static const Double_t fgkEndLadderMountingBlockPosition[2];
  void SetEndLadderCarbonFiberJunctionCombiTransMatrix(Int_t);
  TGeoCombiTrans *EndLadderCarbonFiberJunctionCombiTransMatrix[fgkEndLadderCarbonFiberJunctionCombiTransNumber];
  /////////////////////////////////////////////////////////////////////////
  // Cooling Tube Support (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////
  static const Int_t fgkCoolingTubeSupportCombiTransNumber = 2;
  static const Double_t fgkCoolingTubeSupportRmax;
  static const Double_t fgkCoolingTubeSupportRmin;
  static const Double_t fgkCoolingTubeSupportLength;
  static const Double_t fgkCoolingTubeSupportHeight;
  static const Double_t fgkCoolingTubeSupportWidth;
  static const Double_t fgkCoolingTubeSupportSeparation;
  static const Double_t fgkCoolingTubeSupportToCarbonFiber;
  void SetCoolingTubeSupportCombiTransMatrix();
  TGeoCombiTrans *CoolingTubeSupportCombiTransMatrix[fgkCoolingTubeSupportCombiTransNumber];
  TGeoMedium* fgkSSDTubeHolderMedium;
  /////////////////////////////////////////////////////////////////////////////////
  // Cooling Tube (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////////////
  static const Int_t    fgkCoolingTubeCombiTransNumber = 2;
  static const Double_t fgkCoolingTubeRmax;
  static const Double_t fgkCoolingTubeRmin;
  static const Double_t fgkCoolingTubeLength;
  static const Double_t fgkCoolingTubeSeparation;
  void SetCoolingTubeCombiTransMatrix();
  TGeoCombiTrans *CoolingTubeTransMatrix[fgkCoolingTubeCombiTransNumber];
  TGeoMedium* fgkSSDCoolingTubeWater;
  TGeoMedium* fgkSSDCoolingTubePhynox;
  /////////////////////////////////////////////////////////////////////////
  // SSD Mounting Block Parameters (lengths are in mm and angles in degrees)
  /////////////////////////////////////////////////////////////////////////
  static const Double_t fgkSSDMountingBlockLength[3];  
  static const Double_t fgkSSDMountingBlockHeight[4];
  static const Double_t fgkSSDMountingBlockWidth;
  static const Double_t fgkSSDMountingBlockHoleTrapezoidAngle;
  static const Double_t fgkSSDMountingBlockHoleTrapezoidHeight;
  static const Double_t fgkSSDMountingBlockHoleTrapezoidUpBasis;
  static const Double_t fgkSSDMountingBlockHoleTubeLength[2];
  static const Double_t fgkSSDMountingBlockHoleTubeWidth[2]; 
  static const Double_t fgkSSDMountingBlockHoleRadius;
  static const Double_t fgkSSDMountingBlockScrewHoleEdge;
  static const Double_t fgkSSDMountingBlockScrewHoleHeigth;
  static const Double_t fgkSSDMountingBlockScrewHoleRadius[2];
  TGeoMedium* fgkSSDMountingBlockMedium;
  /////////////////////////////////////////////////////////////////////////
  // LadderSegment 
  /////////////////////////////////////////////////////////////////////////
  static const Int_t fgkLadderSegmentCombiTransNumber = 5;
  void SetLadderSegmentCombiTransMatrix();
  TGeoCombiTrans *LadderSegmentCombiTransMatrix[fgkLadderSegmentCombiTransNumber];
  static const Int_t fgkEndLadderSegmentCombiTransNumber = 4;
  void SetEndLadderSegmentCombiTransMatrix(Int_t);
  TGeoCombiTrans *EndLadderSegmentCombiTransMatrix[fgkEndLadderSegmentCombiTransNumber];
  /////////////////////////////////////////////////////////////////////////
  // LadderCables 
  /////////////////////////////////////////////////////////////////////////
  void SetLadderCableCombiTransMatrix(Int_t);
  TGeoCombiTrans *LadderCableCombiTransMatrix[4];
  /////////////////////////////////////////////////////////////////////////
  // Mother Volume 
  /////////////////////////////////////////////////////////////////////////
  TGeoVolume *fMotherVol;  
  TGeoVolume* GetMotherVolume() const { return fMotherVol;};
  /////////////////////////////////////////////////////////////////////////
  // Color Display 
  /////////////////////////////////////////////////////////////////////////
  Int_t fColorCarbonFiber;   
  Int_t fColorRyton;         
  Int_t fColorPhynox;        
  Int_t fColorSilicon;       
  Int_t fColorAl;   
  Int_t fColorKapton;         
  Int_t fColorPolyhamide;    
  Int_t fColorStiffener;
  Int_t fColorEpoxy;
};
#endif

