#ifndef AliMFTLadderSegmentation_H
#define AliMFTLadderSegmentation_H 

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTbase
/// \class AliMFTLadderSegmentation
/// \brief Description of the virtual segmentation of a ladder
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>
/// \date June 9th, 2015

#include "TClonesArray.h"
#include "AliMFTVSegmentation.h"
#include "AliMFTChipSegmentation.h"

//====================================================================================================================================================

class AliMFTLadderSegmentation : public AliMFTVSegmentation {

public:

  AliMFTLadderSegmentation();
  AliMFTLadderSegmentation(UInt_t uniqueID);
  AliMFTLadderSegmentation(const AliMFTLadderSegmentation& ladder);
  AliMFTLadderSegmentation& operator=(const AliMFTLadderSegmentation& ladder);

  virtual ~AliMFTLadderSegmentation() { if(fChips){fChips->Delete(); delete fChips; fChips=NULL;} }
  virtual void Print(Option_t* opt="");
  virtual void Clear(const Option_t* /*opt*/) { if(fChips){fChips->Clear();} }
  
  AliMFTChipSegmentation* GetSensor(Int_t sensor) const ;


  void CreateSensors();
  
  /// \brief Returns number of Sensor on the ladder
  Int_t GetNSensors() const { return fNSensors; };
  /// \brief Set number of Sensor on the ladder
  void SetNSensors(Int_t val) {fNSensors = val;};
  
  AliMFTChipSegmentation* GetChip(Int_t chipNumber) const {return GetSensor(chipNumber);};

private:
  
  Int_t fNSensors;      ///< \brief Number of Sensors holded by the ladder
  TClonesArray *fChips; ///< \brief Array of pointer to AliMFTChipSegmentation

  /// \cond CLASSIMP
  ClassDef(AliMFTLadderSegmentation, 1);
  /// \endcond

};

//====================================================================================================================================================
	
#endif

