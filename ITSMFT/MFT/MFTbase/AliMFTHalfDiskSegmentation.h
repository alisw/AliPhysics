#ifndef AliMFTHalfDiskSegmentation_H
#define AliMFTHalfDiskSegmentation_H 

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTbase
/// \class AliMFTHalfDiskSegmentation
/// \brief Class for the description of the structure a Half-Disk
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>
/// \date June 9th, 2015

#include "TXMLEngine.h"

#include "AliMFTLadderSegmentation.h"
#include "AliMFTVSegmentation.h"

class TClonesArray;

//====================================================================================================================================================

class AliMFTHalfDiskSegmentation : public AliMFTVSegmentation {

public:

  AliMFTHalfDiskSegmentation();
  AliMFTHalfDiskSegmentation(UInt_t uniqueID);
  AliMFTHalfDiskSegmentation(const AliMFTHalfDiskSegmentation& pt);
  
  virtual ~AliMFTHalfDiskSegmentation();

  virtual void Clear(const Option_t* /*opt*/);
  
  virtual void Print(Option_t* opt="");
  
  void CreateLadders(TXMLEngine* xml, XMLNodePointer_t node);
  
  /// \brief Get the number of Ladder on the Half-Disk really constructed 
  Int_t    GetNLaddersBuild()  const {return fLadders->GetEntriesFast();};

  /// \brief Get the number of Ladder on the Half-Disk
  Int_t    GetNLadders()  const {return fNLadders;};
  
  /// \brief Set the number of Ladder on the Half-Disk
  void    SetNLadders(Int_t val)   {fNLadders = val;};

  
  /// \brief Returns pointer to the ladder segmentation object
  /// \param iLadder Int_t : ladder number on the Half-Disk
  AliMFTLadderSegmentation* GetLadder(Int_t iLadder) { return ( (iLadder>=0 && iLadder<GetNLadders())  ? (AliMFTLadderSegmentation*) fLadders->At(iLadder) : NULL )  ; }
  
  /// \brief Returns the Z position of the half-disk
  Double_t GetZ() const {const Double_t *pos = GetTransformation()->GetTranslation(); return pos[2];};

  Int_t GetNChips();
  
private:
  
  Int_t fNLadders; ///< \brief Number of ladder holded by the half-disk

  TClonesArray *fLadders; ///< \brief Array of pointer to AliMFTLadderSegmentation
  
  /// \cond CLASSIMP
  ClassDef(AliMFTHalfDiskSegmentation, 1);
  /// \endcond

};

//====================================================================================================================================================
	
#endif

