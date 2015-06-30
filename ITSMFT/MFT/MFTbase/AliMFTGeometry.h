#ifndef AliMFTGeometry_H
#define AliMFTGeometry_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTbase
/// \class AliMFTGeometry
/// \brief Class Handling both Virutal Segmentation and Real Volumes
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>
/// \date June 9th, 2015

# include "AliGeometry.h"
# include "AliMFTSegmentation.h"

class AliMFTGeometryBuilder;

//__________________________________________________________________

class AliMFTGeometry : public AliGeometry
{
public:
  /// \brief Retuns MFT Geometry singleton object
  static AliMFTGeometry* Instance();
  
  enum ObjectTypes {kHalfMFTType, kHalfDiskType, kPlaneType, kLadderType, kSensorType};
  
  virtual ~AliMFTGeometry();
  virtual void GetGlobal(const AliRecPoint * p, TVector3 & pos, TMatrixF & mat) const {};
  virtual void GetGlobal(const AliRecPoint * p, TVector3 & pos) const {};
  virtual Bool_t Impact(const TParticle * particle) const             {return kFALSE;};
  void   Build();
  
  /// \brief Returns Object type based on Unique ID provided
  Int_t GetObjectType(UInt_t uniqueID)  const {return ((uniqueID>>13)&0x7);};
  /// \brief Returns Half-MFT ID based on Unique ID provided
  Int_t GetHalfMFTID(UInt_t uniqueID)   const {return ((uniqueID>>12)&0x1);};
  /// \brief Returns Half-Disk ID based on Unique ID provided
  Int_t GetHalfDiskID(UInt_t uniqueID)  const {return ((uniqueID>>9)&0x7);};
  /// \brief Returns Ladder ID based on Unique ID provided
  Int_t GetLadderID(UInt_t uniqueID)    const {return ((uniqueID>>4)&0x1F);};
  /// \brief Returns Sensor ID based on Unique ID provided
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

  static AliMFTGeometry* fgInstance; ///< \brief  Singleton instance
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

