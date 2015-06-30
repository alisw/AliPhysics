#ifndef AliMFTChip_H
#define AliMFTChip_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTbase
/// \class AliMFTChip
/// \brief Class describing geometry of MFT CMOS MAP Chip
///
/// Units are cm and deg
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>
/// \date June 9th, 2015

#include "TNamed.h"
#include "TGeoVolume.h"

class AliMFTLadderSegmentation;
class AliMFTChipSegmentation;
//=============================================================================================


class AliMFTChip : public TNamed {
  
public:
  
  AliMFTChip();
  AliMFTChip(AliMFTChipSegmentation *segmentation, const char * ladderName);
  
  virtual ~AliMFTChip();
  
  TGeoVolume * CreateVolume();
  void GetPosition(AliMFTLadderSegmentation * ladderSeg, Int_t iChip, Double_t *pos);


private:
  
  /// \cond CLASSIMP
  ClassDef(AliMFTChip, 1);
  /// \endcond

};

//=============================================================================================

#endif

