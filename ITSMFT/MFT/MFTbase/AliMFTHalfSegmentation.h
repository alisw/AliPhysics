#ifndef AliMFTHalfSegmentation_H
#define AliMFTHalfSegmentation_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTbase
/// \class AliMFTHalfSegmentation
/// \brief Segmentation class for each half of the ALICE Muon Forward Tracker
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>
/// \date June 9th, 2015

#include "TXMLEngine.h"
#include "AliMFTSegmentation.h"
#include "AliMFTConstants.h"
#include "AliMFTVSegmentation.h"

//====================================================================================================================================================

class AliMFTHalfDiskSegmentation;

class AliMFTHalfSegmentation : public AliMFTVSegmentation {

public:
  
  AliMFTHalfSegmentation();
  AliMFTHalfSegmentation(const Char_t *initFile, const Short_t id);
  AliMFTHalfSegmentation(const AliMFTHalfSegmentation &source);

  virtual ~AliMFTHalfSegmentation();
  virtual void Clear(const Option_t* /*opt*/);
  
  Bool_t GetID() const {return (GetUniqueID()>>12);};
  
  Int_t GetNHalfDisks() const { return fMFTHalfDisks->GetEntries(); }

  AliMFTHalfDiskSegmentation* GetHalfDisk(Int_t iDisk) const { if (iDisk>=0 && iDisk<fMFTHalfDisks->GetEntries()) return (AliMFTHalfDiskSegmentation*) fMFTHalfDisks->At(iDisk); else return NULL; }
 
private:
  
  void FindHalf(TXMLEngine* xml, XMLNodePointer_t node, XMLNodePointer_t &retnode);
  void CreateHalfDisks(TXMLEngine* xml, XMLNodePointer_t node);

  TClonesArray *fMFTHalfDisks; ///< \brief Array of pointer to AliMFTHalfDiskSegmentation

  /// \cond CLASSIMP
  ClassDef(AliMFTHalfSegmentation, 1);
  /// \endcond
  
};

//====================================================================================================================================================

#endif

