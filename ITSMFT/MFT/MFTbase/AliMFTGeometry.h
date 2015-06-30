#ifndef AliMFTGeometry_H
#define AliMFTGeometry_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved.
 *
 * See cxx source for full Copyright notice
 */
/** @file    AliMFTGeometry.h
 @author  Raphael Tieulent <raphael.tieulent@cern.ch>
 @date    June 9th 2015
 @brief   Geometry mananger for the MFT
 */
//____________________________________________________________________
//
// Muon Forward Tracker Geometry.
//
// This class is a singleton that handles the geometry parameters of
// the MFT detectors.
//
# include "AliGeometry.h"
# include "AliMFTSegmentation.h"

class AliMFTGeometryBuilder;
class TParticle;
class AliRecPoint;
class TVector3;

//__________________________________________________________________
/** @brief Singleton object of MFT geometry descriptions and parameters.
 This class is a singleton that handles the geometry parameters of
 the MFT detectors.
 
 The actual code is done by various separate classes.
 @endverbatim
 
 @ingroup MFTbase
 */
class AliMFTGeometry : public AliGeometry
{
public:
  /**
   * singleton access
   *
   * @return Singleton
   */
  static AliMFTGeometry* Instance();
  enum ObjectTypes {kHalfMFTType, kHalfDiskType, kPlaneType, kLadderType, kSensorType};
  
  virtual ~AliMFTGeometry();
  virtual void GetGlobal(const AliRecPoint * p, TVector3 & pos, TMatrixF & mat) const {};
  virtual void GetGlobal(const AliRecPoint * p, TVector3 & pos) const {};
  virtual Bool_t Impact(const TParticle * particle) const             {return kFALSE;};
  void   Build();
  
  Int_t GetObjectType(UInt_t uniqueID)  const {return ((uniqueID>>13)&0x7);};
  Int_t GetHalfMFTID(UInt_t uniqueID)   const {return ((uniqueID>>12)&0x1);};
  Int_t GetHalfDiskID(UInt_t uniqueID)  const {return ((uniqueID>>9)&0x7);};
  Int_t GetLadderID(UInt_t uniqueID)    const {return ((uniqueID>>4)&0x1F);};
  Int_t GetSensorID(UInt_t uniqueID)    const {return (uniqueID&0xF);}; 
  
  UInt_t GetObjectID(ObjectTypes type, Int_t half=0, Int_t disk=0, Int_t ladder=0, Int_t chip=0) const;
  
  /// \brief Returns TGeo ID of the volume describing the sensors
  Int_t GetSensorVolumeID()    const {return fSensorVolumeId;};
  /// \brief Set the TGeo ID of the volume describing the sensors
  void  SetSensorVolumeID(Int_t val)   { fSensorVolumeId= val;};

  /// \brief Returns pointer to the segmentation
  AliMFTSegmentation * GetSegmentation() const {return fSegmentation;};
  
  Bool_t Hit2PixelID(Double_t xHit, Double_t yHit, Double_t zHit, Int_t detElemID, Int_t &xPixel, Int_t &yPixel) const;
  void GetPixelCenter(Int_t xPixel, Int_t yPixel, Int_t detElemID, Double_t &xCenter, Double_t &yCenter, Double_t &zCenter ) const ;
  Int_t GetDiskNSensors(Int_t diskId) const;
  Int_t GetDetElemLocalID(Int_t detElem) const;

  void LoadSegmentation();

  
private:

  static AliMFTGeometry* fgInstance; // Singleton instance
  AliMFTGeometry();
  AliMFTGeometry(const char* name);
  
  AliMFTGeometryBuilder* fBuilder; ///< \brief Geometry Builder
  AliMFTSegmentation *fSegmentation; ///< \brief Segmentation of the detector
  Int_t fSensorVolumeId; ///< \brief ID of the volume describing the CMOS Sensor
  
  /// \cond CLASSIMP
  ClassDef(AliMFTGeometry, 1);
  /// \endcond

};


#endif

