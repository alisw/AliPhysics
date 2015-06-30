#ifndef AliMFTTrackParam_H
#define AliMFTTrackParam_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup MFTrec
/// \class AliMFTTrackParam
/// \brief Class holding the parameter of a MFT Standalone Track
///
///
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>, IPN-Lyon
/// \date April 28th, 2015

#include "TObject.h"

//=============================================================================================

class AliMFTTrackParam : public TObject {
  
  public:
  
  AliMFTTrackParam();
  virtual ~AliMFTTrackParam();


  protected:
  
  
  /// \cond CLASSIMP
  ClassDef(AliMFTTrackParam,1);
  /// \endcond
};

//=============================================================================================

#endif
