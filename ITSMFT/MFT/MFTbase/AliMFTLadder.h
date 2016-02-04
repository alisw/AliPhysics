#ifndef AliMFTLadder_H
#define AliMFTLadder_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTbase
/// \class AliMFTLadder
/// \brief Class building the Ladder geometry
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>
/// \date June 9th, 2015

#include "TNamed.h"
#include "TGeoVolume.h"

class AliMFTLadderSegmentation;
class AliMFTFlex;
//=============================================================================================


class AliMFTLadder : public TNamed {
  
public:
  
  AliMFTLadder();
  AliMFTLadder(AliMFTLadderSegmentation *segmentation);
  
  virtual ~AliMFTLadder();
  
  TGeoVolume * CreateVolume();
  void CreateSensors();

private:

  const static Double_t kLadderDeltaY;      ///< \brief Ladder size along Y direction (height)
  const static Double_t kLadderDeltaZ;      ///< \brief Ladder size along Z direction (thickness)
  AliMFTLadderSegmentation *fSegmentation;  ///< \brief Virtual Segmentation object of the ladder
  AliMFTFlex      * fMFTFlex;               ///< \brief Flex object (\todo to be removed ?)
  TGeoVolumeAssembly * fLadderVolume;               ///< \brief Pointer to the Volume holding the ladder geometry

  
  /// \cond CLASSIMP
  ClassDef(AliMFTLadder, 1);
  /// \endcond

};

//=============================================================================================

#endif

