#ifndef AliMFTGeometryBuilder_H
#define AliMFTGeometryBuilder_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTbase
/// \class AliMFTGeometryBuilder
/// \brief Class describing MFT Geometry Builder
///
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>
/// \date June 9th, 2015

#include "TNamed.h"

class AliMFTSegmentation;
class TClonesArray;
class TGeoVolumeAssembly;

//=============================================================================================


class AliMFTGeometryBuilder : public TNamed {
  
public:
  
  AliMFTGeometryBuilder();
  
  virtual ~AliMFTGeometryBuilder();
  
  void LoadSegmentation();
  
  void BuildGeometry();
  
  
protected:
  TClonesArray *fMFTHalves;                  // ! list

private:
  
  /// \cond CLASSIMP
  ClassDef(AliMFTGeometryBuilder, 1);
  /// \endcond

};

//=============================================================================================

#endif

