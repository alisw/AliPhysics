#ifndef ALIEMCALTRIGGERRAWDIGIT_H
#define ALIEMCALTRIGGERRAWDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <iosfwd>
#include "AliEMCALRawDigit.h" 
#include "AliEMCALTriggerTypes.h" 
#include "AliLog.h"

/// \class AliEMCALTriggerRawDigit
/// \brief EMCal trigger raw digits 
/// \ingroup EMCALbase
/// \author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
///
///  A trigger Digit containing raw samples.
class AliEMCALTriggerRawDigit : public AliEMCALRawDigit 
{

public:
	
  /// \brief Dummy constructor 
  AliEMCALTriggerRawDigit();

  /// \brief Constructor 
  /// \param ID FastOR absolute ID
  /// \param timeSamples ADC time samples
  /// \param nSamples Number of time samples
  AliEMCALTriggerRawDigit(Int_t id, Int_t timeSamples[], Int_t nSamples);
  
  /// \brief Destructor 
  virtual ~AliEMCALTriggerRawDigit();
  
  /// \brief Set trigger bit for a given trigger type
  /// \param type Trigger type (L0, L1 trigger for gamma / jet at threshold)
  /// \param mode Reconstruction mode (raw data or recalculated)
  void    SetTriggerBit(const int type, const Int_t mode) 
  { fTriggerBits = (fTriggerBits | (1 << (type + kTriggerTypeEnd * mode))) ; }
  
  /// \brief Set L0 time
  /// \param timebin L0 time
  /// \return True if L0 time is set, false if setting failed (time is larger than the maximum L0 time)
  ///
  /// L0 time is the time bin the digit contributed to a patch above threshold
  Bool_t  SetL0Time( Int_t timebin );
  
  /// \brief Get trigger bit
  /// \param type Trigger type (L0, L1 trigger for gamma / jet at threshold)
  /// \param mode Reconstruction mode (raw data or recalculated)
  /// \return 1 if bit for trigger is set, 0 if not set
  Int_t   GetTriggerBit(const TriggerType_t type, const Int_t mode) const;
  
  /// \brief Get bitmap connected to the trigger raw digits
  /// \return Bitmap with trigger bits handled by the raw digit
  Int_t   GetTriggerBits() const { return fTriggerBits ; }
  
	/// \brief Get L0 time for a given index
  /// \param[in] index Index of the time sample
  /// \param[out] time L0 time
  /// \return True if the, false if the index is out-of-range
  ///
  /// L0 time is the time bin the digit contributed to a patch above threshold
  Bool_t  GetL0Time(const Int_t index, Int_t& time) const;

  /// \brief Get L0 times
  /// \param[out ]times Output array for L0 time (size must match the number of L0 times in the raw digit)
  ///
  /// L0 time is the time bin the digit contributed to a patch above threshold
  /// \return Always true
  Bool_t  GetL0Times(Int_t times[]            ) const;

  /// Get the number of L0 times handled by this digit
  /// \return Number of L0 times (max: 10)
  Int_t   GetNL0Times(                        ) const { return fNL0Times ; }
  
  /// \brief Get L0 time sum
  /// \param time Time for which to start the time integration
  /// \return L0 time sum
  /// 
  /// Integrates the time samples starting from time over for 
  /// time samples
  Int_t   GetL0TimeSum(const Int_t time) const;
  
  /// \brief Set the L1 time sum
  /// \param ts L1 time sum
  void    SetL1TimeSum(Int_t ts) 
  { if (fL1TimeSum >= 0) AliWarning("You're overwriting digit time sum! Please check"); fL1TimeSum = ts ; }

  /// \brief Accessor for L1 time sum
  /// \return L1 time sum
  Int_t   GetL1TimeSum(        ) const { return fL1TimeSum   ; }
  
  /// \brief Set the L1 subregion
  /// \param sr L1 subregion
  void    SetL1SubRegion(Int_t sr) 
  { if (fL1SubRegion >= 0) AliWarning("You're overwriting digit subregion! Please check"); fL1SubRegion = sr ; }

  /// \brief Get the L1 subregion
  /// \return L1 subregion
  Int_t   GetL1SubRegion(      ) const { return fL1SubRegion ; }
   
  /// \brief Dump raw digit info
  ///
  /// See PrintStream for more information
  virtual void Print(const Option_t* opt) const;

  /// \brief Print digit information on the output stream
  /// \param stream Stream where to print digit information
  virtual void PrintStream(std::ostream &stream) const;
	
private: 
 
  AliEMCALTriggerRawDigit           (const AliEMCALTriggerRawDigit &cd); // Not implemented
  AliEMCALTriggerRawDigit &operator=(const AliEMCALTriggerRawDigit &cd); // Not implemented
  
  Int_t   fTriggerBits; ///< Trigger bits: Bitmap with all trigger bits (L0, L1) for raw data and recalculation mode
  Int_t   fNL0Times;    ///< Number of L0 times
  Int_t   fL0Times[10]; ///< L0 times: Time bins the patch corresponds to a patch above threshold
  
  Int_t   fL1TimeSum;   ///< L1 time sum
  Int_t   fL1SubRegion; ///< Subregion
  
  ClassDef(AliEMCALTriggerRawDigit,2) ;
};

/// \brief Output stream operator for AliEMCALTriggerRawDigit
/// \param stream Output stream used for printing the digit
/// \param dig Digit to be printed
/// \return Stream after printing 
std::ostream &operator<<(std::ostream &stream, const AliEMCALTriggerRawDigit &dig);

#endif //ALIEMCALTRIGGERRAWDIGIT_H

