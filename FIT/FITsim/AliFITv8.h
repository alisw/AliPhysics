#ifndef FITV8_H
#define FITV8_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////
// Full geometry hits classes for detector: FIT //
///////////////////////////////////////////////////

#include "AliFIT.h"
#include "TGeoVolume.h"
#include "TGraph.h"
#include <TGeoMatrix.h>
#include <TVirtualMC.h>
#include <sstream>
#include <string>
#include <vector>

class AliFITv8 : public AliFIT {

public:
  enum constants {
    kAir = 1,
    kVac = 3,
    kCeramic = 4,
    kGlass = 6,
    kOpAir = 7,
    kAl = 15,
    kOpGlass = 16,
    kOpGlassCathode = 19,
    kSensAir = 22
  };

  AliFITv8();
  AliFITv8(const char *name, const char *title);
  AliFITv8(const AliFITv8 &o)
      : AliFIT(), fIdSens1(0), fIdSens2(0), fPMTeff(0x0) {
    ((AliFITv8 &)o).Copy(*this);
  }

  AliFITv8 &operator=(const AliFITv8 &) { return *this; }
  virtual ~AliFITv8();
  virtual void CreateGeometry();
  virtual void DefineOpticalProperties();
  virtual void AddAlignableVolumes() const;
  virtual void CreateMaterials();
  virtual void Init();
  virtual Int_t IsVersion() const { return 0; }
  Bool_t RegisterPhotoE(Double_t energy);
  void SetVZEROGeo(TGeoVolume *alice);
  void SetVZEROGeoOld(TGeoVolume *alice);
  void SetOneMCP(TGeoVolume *stl);
  virtual void StepManager();

protected:
  // T0+
  Int_t fIdSens1;   // Sensitive volume in FIT
  Int_t fIdSens2;   // Sensitive volume in FIT
  TGraph *fPMTeff;  // MCP-PMT detection efficiency
  Int_t fSenseless; // Senseless hit entry

  // V0+
  Int_t GetCellId(Int_t *vol);
  Int_t GetCellIdOld(Int_t *vol);
  Int_t nSectors, nRings;
  Int_t fIdV0Plus[16][5]; // Sensitive volumes [nSectors][nRings], if modified
                          // then update the construct in .cxx
  Int_t fCellId;          // Scintillator cell number

private:
  // Optical properties reader: e-Energy, abs-AbsorptionLength[cm], n-refractive
  // index
  Int_t ReadOptProperties(const std::string inputFilePath, Float_t **e,
                          Double_t **de, Float_t **abs, Float_t **n,
                          Float_t **qe, Int_t &kNbins) const;
  void FillOtherOptProperties(Float_t **efficAll, Float_t **rindexAir,
                              Float_t **absorAir, Float_t **rindexCathodeNext,
                              Float_t **absorbCathodeNext, Double_t **efficMet,
                              Double_t **aReflMet, const Int_t kNbins) const;
  void DeleteOptPropertiesArr(
      Float_t **e, Double_t **de, Float_t **abs, Float_t **n,
      Float_t **efficAll, Float_t **rindexAir, Float_t **absorAir,
      Float_t **rindexCathodeNext, Float_t **absorbCathodeNext,
      Double_t **efficMet,
      Double_t **aReflMet) const; // should be called at the very end of the
                                  // simulation to free the memory

  // V0+ parameters related to geometry
  Double_t fV0PlusR0, fV0PlusR1, fV0PlusR2, fV0PlusR3, fV0PlusR4, fV0PlusR5,
      fV0PlusR6;
  Double_t fV0PlusSciWd, fV0PlusFraWd, fV0PlusZposition;
  Float_t fV0PlusnMeters;

  // V0+ parameters related to light production:
  Double_t fV0PlusLightYield;       // Lightyield in BC404 (from V0A)
  Double_t fV0PlusLightAttenuation; // LightAttenuation in fibers (from V0A)
  Double_t
      fV0PlusFibToPhot; // Loss in Fibers - Photocathode Connection (from V0A)

  /* ======= Rihan: Variables and Functions from Geometry.cxx ======= */

  // Strings for volume names, etc.
  static const std::string sDetectorName; // = "FV0";   // String Initializing
                                          // is not supported by C++11 ??
  static const std::string sScintillatorName; // = "SCINT";
  static const std::string sPlasticName;      // = "PLAST";
  static const std::string sSectorName;       // = "SECTOR";
  static const std::string sCellName;         // = "CELL";
  static const std::string
      sScintillatorSectorName; // = sScintillatorName + sSectorName;
  static const std::string
      sScintillatorCellName; // = sScintillatorName + sCellName;
  static const std::string sPlasticSectorName; // = sPlasticName + sSectorName;
  static const std::string sPlasticCellName;   // = sPlasticName + sCellName;

  static const std::string sFiberName;        // = "FIBER";
  static const std::string sScrewName;        // = "SCREW";
  static const std::string sScrewHolesCSName; // = "FV0SCREWHOLES";
  static const std::string sRodName;          // = "ROD";
  static const std::string sRodHolesCSName;   // = "FV0RODHOLES";
  static const std::string sContainerName;    // = "CONTAINER";

  // static const char  sCellTypes[sNumberOfCellSectors]{'a', 'b', 'b', 'a'};

  // Cell and scintillator constants
  static const int sNumberOfCellSectors; // = 4;   // < Number of cell sectors
                                         // for one half of the detector
  static const int sNumberOfCellRings;   // = 5;   // < Number of cell rings
  // Container constants
  static const float
      sDzContainer; // = 30.0;     // <-- Depth of the metal container
  static const float sDrContainerHole; // = 4.05;     // <-- Radius of the beam
                                       // hole in the metal container
  static const float sXShiftContainerHole; // = -0.15;    // <-- x-shift of the
                                           // beam hole in the metal container
  static const float sDrMaxContainerBack;  // = 83.1;     // <-- Outer radius of
                                           // the container backplate
  static const float sDzContainerBack; // = 1.00;     // <-- Thickness of the
                                       // container backplate
  static const float sDrMinContainerFront; // = 45.7;     // <-- Inner radius of
                                           // the container frontplate
  static const float sDrMaxContainerFront; // = 83.1;     // <-- Outer radius of
                                           // the container frontplate
  static const float sDzContainerFront; // = 1.00;     // <-- Thickness of the
                                        // container frontplate
  static const float
      sDxContainerStand; // = 40.0;     // <-- Width of the container stand
  static const float sDyContainerStand;   // = 3.00;     // <-- Height of the
                                          // container stand at its center in x
  static const float sDrMinContainerCone; // = 24.3;     // <-- Inner radius at
                                          // bottom of container frontplate cone
  static const float sDzContainerCone;    // = 16.2;     // <-- Depth of the
                                          // container frontplate cone
  static const float sThicknessContainerCone; // = 0.60;     // <-- Thickness of
                                              // the container frontplate cone
  static const float
      sXYThicknessContainerCone; // = 0.975;    // <-- Radial thickness in the
                                 // xy-plane of container cone
  static const float
      sDrMinContainerOuterShield; // = 82.5;     // <-- Inner radius of outer
                                  // container shield
  static const float
      sDrMaxContainerOuterShield; // = 82.65;    // <-- Outer radius of outer
                                  // container shield
  static const float
      sDrMinContainerInnerShield; // = 4.00;     // <-- Inner radius of the
                                  // inner container shield
  static const float
      sDrMaxContainerInnerShield; // = 4.05;     // <-- Outer radius of inner
                                  // container shield
  static const float
      sDxContainerCover; // = 0.15;     // <-- Thickness of the container cover
  static const float sDxContainerStandBottom; // = 38.5;     // <-- Width of the
                                              // bottom of the container stand
  static const float
      sDyContainerStandBottom; // = 2.00;     // <-- Thickness of the bottom of
                               // the container stand

  // General geometry constants
  static const float sEpsilon; // = 0.01;     //   Used to make one spatial
                               // dimension infinitesimally larger than other
  static const float
      sDzScintillator; // = 4.00;     //   Thickness of the scintillator
  static const float
      sDzPlastic; // = 1.00;     //   Thickness of the fiber plastics
  static const float sXGlobal; // = 0.00;     //   Global x-position of the
                               // geometrical center of scintillators
  static const float sYGlobal; // = 0.00;     //   Global y-position of the
                               // geometrical center of scintillators
  static const float
      sDxHalvesSeparation; // = 0.00;     //   Separation between the left and
                           // right side of the detector
  static const float
      sDyHalvesSeparation; // = 0.00;     //   y-position of the right detector
                           // part relative to the left part
  static const float
      sDzHalvesSeparation; // = 0.00;     //   z-position of the right detector
                           // part relative to the left part
  // Local position constants
  static const float
      sXScintillator; // = sDxContainerCover;        // x-position of the right
                      // half of the scintillator.
  // FT0 starts at z=328
  static const float sZGlobal; // = 320 - sDzScintillator/2;  //   Global z-pos
                               // of geometrical center of scintillators

  static const float dxHoleCut; // = 0.20;   // <-- width of extension of hole
                                // 1,2,7 in "a"-type cell
  static const float
      sXShiftInnerRadiusScintillator; // =-0.15;   // <-- Shift of the inner
                                      // radius origin of the scintillators.
  static const float
      sDxHoleExtensionScintillator; // = 0.200;  // <-- Extension of the
                                    // scintillator holes for the metal rods
  static const float sDrHoleSmallScintillator; // = 0.265;  // <-- Radius of the
                                               // small scintillator screw hole
  static const float sDrHoleLargeScintillator; // = 0.415;  // <-- Radius of the
                                               // large scintillator screw hole
  static const float
      sDrSeparationScint;   // = 0.03+0.04; // <-- Separation b/w scint cells =
                            // paint width + half of separation
  static const float xHole; // = sDrSeparationScint + dxHoleCut; // x-position
                            // of holes 1, 2 and 7 in the "a" cell

  static const float sZScintillator; // = 0;
  static const float sXShiftScrews;  // = sXScintillator;
  static const float
      sZPlastic; // = sZScintillator + sDzScintillator / 2 + sDzPlastic / 2;
  static const float sZFiber; // = (sZPlastic + sZContainerFront) / 2; //
                              // z-position of the fiber volumes.

  static const float
      sZContainerBack; //   /// z-position of the container backplate.
  static const float
      sZContainerFront; //   /// z-position of the container frontplate.
  static const float
      sZContainerMid; //   /// z-position of the center of the container.
  static const float
      sZCone; //   /// z-position of the container frontplate cone.

  /// Screw Properties:
  static const float sZShiftScrew; // = 0;                                 //
                                   // <-- z shift of the screws. 0 means they
                                   // are aligned with the scintillator.
  static const int sNumberOfScrewTypes; // = 6; // <-- Number of the different
                                        // screw types.
  static const float
      sDrMinScrewTypes[6]; //{0.25, 0.25, 0.4, 0.4, 0.4, 0.4};     // <-- Radii
                           //of the thinner part of the screw types.
  static const float
      sDrMaxScrewTypes[6]; //{0, 0.5, 0.6, 0.6, 0.6, 0};           // <-- Radii
                           //of the thicker part of the screw types.
  static const float
      sDzMaxScrewTypes[6]; //{6.02, 13.09, 13.1, 23.1, 28.3, 5};   // <-- Length
                           //of the thinner part of the screw types.
  static const float
      sDzMinScrewTypes[6]; //{0, 6.78, 6.58, 15.98, 21.48, 0};     // <-- Length
                           //of the thicker part of the screw types.

  /// Rod Properties:
  static const float
      sZShiftRod; // = -0.05;                       // <-- z shift of the rods.
                  // 0 = aligned with tht scintillators.
  static const int sNumberOfRodTypes; // = 4;                           // <--
                                      // Number of the different screw types.
  static const float
      sDxMinRodTypes[4]; //{0.366, 0.344, 0.344, 0.344};   // <-- Width of the
                         //thinner part of the rod types.
  static const float
      sDxMaxRodTypes[4]; //{0.536, 0.566, 0.566, 0.716};   // <-- Width of the
                         //thicker part of the rod types.
  static const float
      sDyMinRodTypes[4]; //{0.5, 0.8, 0.8, 0.8};           // <-- Height of the
                         //thinner part of the rod types.
  static const float
      sDyMaxRodTypes[4]; //{0.9, 1.2, 1.2, 1.2};           // <-- Height of the
                         //thicker part of the rod types.
  static const float
      sDzMaxRodTypes[4]; //{12.5, 12.5, 22.5, 27.7};       // <-- Length of the
                         //thinner part of the rod types.
  static const float
      sDzMinRodTypes[4]; //{7.45, 7.45, 17.45, 22.65};     // <-- Length of the
                         //thicker part of the rod types.

  static const float sCellRingRadii[10]; //!

  std::vector<float>
      mRAvgRing; //! index 0 -> ring1 min, index 1 -> ring1 max and so on..
  std::vector<float>
      mRMinScintillator; //!  <-- Inner radii of scintillator rings (.at(0) ->
                         //!  ring 1, .at(4) -> ring 5)
  std::vector<float>
      mRMaxScintillator; //!  <-- Outer radii of scintillator rings (.at(0) ->
                         //!  ring 1, .at(4) -> ring 5)

  std::vector<std::string>
      mSensitiveVolumeNames; //!  <-- The names of all the sensitive volumes
  std::vector<TGeoMatrix *>
      mSectorTrans; //!  <-- Transformations of sectors (.at(0) -> sector 1)

  std::vector<std::vector<float>>
      mScrewPos; //!  <-- xyz-coordinates of all the screws
  std::vector<std::vector<float>>
      mRodPos;                 //!  <-- xyz-coordinates of all the rods
  std::vector<float> mRScrews; //!  <-- Radial distance to the screw locations
  std::vector<int> mScrewTypeIDs; //!  <-- The type ID of each screw (.at(n) ->
                                  //!  type ID of screw no. n)
  std::vector<float>
      mRScrewAndRod; //!  <-- Radii of the screw and rod positions
  std::vector<float>
      mDrMinScrews; //!  <-- Radii of the thinner part of the screws
  std::vector<float>
      mDrMaxScrews; //!  <-- Radii of the thicker part of the screws
  std::vector<float>
      mDzMaxScrews; //!  <-- Length of the thinner part of the screws
  std::vector<float>
      mDzMinScrews; //!  <-- Length of the thicker part of the screws

  std::vector<float> mDxMinRods; //! Width of the thinner part of the rods
  std::vector<float> mDxMaxRods; //! Width of the thicker part of the rods
  std::vector<float> mDyMinRods; //! Height of the thinner part of the rods
  std::vector<float> mDyMaxRods; //! Height of the thicker part of the rods
  std::vector<float> mDzMaxRods; //! Length of the thinner part of the rods
  std::vector<float> mDzMinRods; //! Length of the thicker part of the rods
  std::vector<int>
      mRodTypeIDs; //! The type ID of each rod (.at(n) -> type ID of rod no. n)

  std::vector<float>
      mRMinFiber; //!  <-- Inner radii of fiber volumes (.at(0) -> fiber 1)
  std::vector<float>
      mRMaxFiber; //!  <-- Outer radii of fiber volumes (.at(0) -> fiber 1)

  std::vector<TGeoMedium *>
      mMediumFiber; //!  Medium of the fiber volumes. at(n) -> medium of the
                    //!  n:th fiber starting from the middle.
  std::vector<TGeoMedium *> mMediumScrewTypes; //!  Medium of the screw types
  std::vector<TGeoMedium *> mMediumRodTypes;   //!  Medium of the rod types

  std::vector<int> mcIdtoRing1map;
  std::vector<int> mcIdtoRing2map;
  std::vector<int> mcIdtoRing3map;
  std::vector<int> mcIdtoRing4map;
  std::vector<int> mcIdtoRing5map;

  TGeoMatrix *mLeftTransformation;  //!  <-- Transformation for the left part of
                                    //!  the detector
  TGeoMatrix *mRightTransformation; //!  <-- Transformation for the right part
                                    //!  of the detector

  int getRingNumberFromMCcellId(int gMCcellId = 0);

  /// Helper function for creating and registering a TGeoTranslation.
  /// \param  name  The name of the translation.
  /// \param  dx    Translation dx.
  /// \param  dy    Translation dy.
  /// \param  dz    Translation dz.
  /// \return The newly created and registered TGeoTranslation.
  TGeoTranslation *createAndRegisterTrans(const std::string &name,
                                          double dx = 0, double dy = 0,
                                          double dz = 0) const;

  /// Helper function for creating and registering a TGeoRotation.
  /// \param  name  The name of the rotation.
  /// \param  dx    Translation phi.
  /// \param  dy    Translation theta.
  /// \param  dz    Translation psi.
  /// \return The newly created and registered TGeoRotation.
  TGeoRotation *createAndRegisterRot(const std::string &name, double phi = 0,
                                     double theta = 0, double psi = 0) const;

  /// Build sector assembly of specified type.
  /// \param  cellName  The type of the cells in the sector assembly.
  /// \return The sector assembly.
  TGeoVolumeAssembly *buildSectorAssembly(const std::string &cellName,
                                          Bool_t isSensitive = kFALSE) const;

  /// Build a sector of specified type and number.
  /// \param  cellType  The type of the cells in the sector.
  /// \param  iSector   The numbering of the sector.
  /// \return The sector.
  TGeoVolumeAssembly *buildSector(const std::string &cellType, int iSector,
                                  Bool_t isSensitive = kFALSE) const;

  const std::string createVolumeName(const std::string &volumeType,
                                     int number = -1) const;

  /// Create the shape for a specified screw.
  /// \param  shapeName   The name of the shape.
  /// \param  screwTypeID The number of the screw type.
  /// \param  xEpsilon    Shrinks or expands the x dimensions of the screw
  /// shape. \param  yEpsilon    Shrinks or expands the y dimensions of the
  /// screw shape. \param  zEpsilon    Shrinks or expands the z dimensions of
  /// the screw shape. \return The screw shape.
  TGeoShape *createScrewShape(const std::string &shapeName, int screwTypeID,
                              float xEpsilon = 0, float yEpsilon = 0,
                              float zEpsilon = 0) const;

  /// Create the shape for a specified rod.
  /// \param  shapeName The name of the shape.
  /// \param  rodTypeID The number of the rod type.
  /// \param  xEpsilon  Shrinks or expands the x dimensions of the rod shape.
  /// \param  yEpsilon  Shrinks or expands the y dimensions of the rod shape.
  /// \param  zEpsilon  Shrinks or expands the z dimensions of the rod shape.
  /// \return The rod shape.
  TGeoShape *createRodShape(const std::string &shapeName, int rodTypeID,
                            float xEpsilon = 0, float yEpsilon = 0,
                            float zEpsilon = 0) const;

  // Assemble the scintillator sectors.
  void
  assembleScintSectors(TGeoVolume *vFV0) const; ///   vFV0 is The FIT V0 volume.
  // Assemble the Plastic scintillators.
  void
  assemblePlasticSectors(TGeoVolume *vFV0) const; /// vFV0 is The FIT V0 volume.
  // Assemble the screwss.
  void assembleScrews(TGeoVolume *vFV0) const;
  // Assemble the rods.
  void assembleRods(TGeoVolume *vFV0) const;
  // Assemble the optical fibers.
  void assembleFibers(TGeoVolume *vFV0) const;
  /// Assemble the metal container.
  void assembleMetalContainer(TGeoVolume *vFV0) const;

  void initializeCellRingRadii();         // Initialize sector Ring Radii.
  void initializeTransformations();       // Initialize common transformations.
  void initializeSectorTransformations(); // Initialize sector transformations.
  // Intialize function for Scintillator + Plastic cells:
  void initializeCells(const std::string &cellType, const float zThickness,
                       const TGeoMedium *medium, const Bool_t isSensitive);

  void initializeScrews(); // <--- Initialize the screw volumes.
  void initializeScrewHoles();
  void initializeScrewTypeMedium();
  void initializeScrewAndRodRadii();
  void initializeScrewAndRodPositionsAndDimensions();
  void addScrewProperties(const int screwTypeID = 1, const int iRing = 0,
                          const float phi = 0.);

  void initializeRods(); // <-- Initialize the rod volumes.
  void initializeRodHoles();
  void initializeRodTypeMedium();
  void addRodProperties(const int rodTypeID = 1, const int iRing = 0);

  void initializeFiberVolumeRadii();
  void initializeFiberMedium();
  void initializeFibers();

  void
  initializeMetalContainer(); //   <-- Initialize the metal container volume.

  ClassDef(AliFITv8, 1) // Class for FIT version 6
};

#endif
