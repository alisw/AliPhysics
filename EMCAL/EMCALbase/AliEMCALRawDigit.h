#ifndef ALIEMCALRAWDIGIT_H
#define ALIEMCALRAWDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <iosfwd>
#include "TObject.h" 

/// \class AliEMCALRawDigit
/// \brief EMCal raw digits object  
/// \ingroup EMCALbase
/// \author R. Guernane, LPSC-IN2P3-CNRS
///
///  A Digit containing raw samples, for trigger.
class AliEMCALRawDigit : public TObject 
{
  
public:
  /// \brief Dummy constructor
  AliEMCALRawDigit();

  /// \brief Constructor
  /// \param ID Module (Tower/FastOR) absolute ID
  /// \param timeSamples Time samples
  /// \param nSamples Number of time samples
  AliEMCALRawDigit(Int_t id, Int_t timeSamples[], Int_t nSamples);

  /// \brief Check for equalness
  /// \param rhs Digit to compare to
  /// \return True if digits have the same ID, false otherwise
  bool operator==(const AliEMCALRawDigit &rhs) const;

  /// \brief Check whether this digit has a smaller ID
  /// \param rhs Digit to compare to
  /// \return True if this digit has a smaller ID, false otherwise
  bool operator<(const AliEMCALRawDigit &rhs) const;
  
  /// \brief Destructor
  ///
  /// Delete array of time samples
  virtual ~AliEMCALRawDigit();
  
  /// \brief Clear digit 
  /// 
  /// Delete array of time samples
  void Clear(Option_t *);
  
  /// \brief Declare digit as sortable
  /// \return Always true
  Bool_t  IsSortable() const { return kTRUE;}
  
  /// \brief Comparison function
  /// \param obj Object to compare to (only AliEMCALRawDigit supported)
  /// \return 1 if id of digit is > id of other digit
  ///         -1 if id of digit is < id of other digit
  ///         0 if indices are equal
  ///
  /// Compares two digits with respect to its Id
  /// to sort according increasing Id
  Int_t   Compare(const TObject* obj) const;
  
  /// \brief Set the module ID
  /// \param id Absolute ID of the module
  void    SetId(Int_t id)           { fId        = id ; }

  /// \brief Set the digit amplitude
  /// \param amp Digit amplitude
  void    SetAmplitude(Float_t amp) { fAmplitude = amp  ; }

  /// \brief Set the digit time
  /// \param time Digit time
  void    SetTime(Float_t time)     { fTime      = time ; }
  
  /// \brief Set the raw time samples
  /// \param timeSamples Array with raw time samples
  /// \param nSamples Number of time samples
  void    SetTimeSamples(const Int_t timeSamples[], const Int_t nSamples);
  
  /// \brief Get the absolute ID of the module
  /// \return ID of the module
  Int_t   GetId()        const { return fId        ; }	

  /// \brief Get the digit amplitude
  /// \return Digit amplitude
  Float_t GetAmplitude() const { return fAmplitude ; }

  /// \brief Get the digit time
  /// \return Digit time
  Float_t GetTime()      const { return fTime      ; }

  /// \brief Get the number of time samples
  /// \return Number of time samples
  Int_t   GetNSamples()  const { return fNSamples  ; }

  /// \brief Access to the first n time samples
  /// \param[out] samples Output array for the time samples
  /// \param ns Number of time samples to be accessed
  /// \return True if the number of samples is within range, false otherwise
  Bool_t  GetSamples(Int_t samples[64], Int_t ns) const;

  /// \brief Get the raw time sample (amplitude and time) for a given indes
  /// \param[in] iSample Index of the time sample
  /// \param[out] sample Raw time sample (amplitude and time)
  /// \return True if the index is within range, false otherwise
  Bool_t  GetTimeSample(const Int_t iSample, Int_t& sample) const;

  /// \brief Access to amplitude and time of a given time sample
  /// \param in iSample Index of the time sample
  /// \param[out] timeBin Time of the time sample and if the sample was ok
  /// \param[out] amp Amplitude of the time sample
  /// \return True if the sample index is within range, false otherwise
  Bool_t  GetTimeSample(const Int_t iSample, Int_t& timeBin, Int_t& amp) const;

  /// \brief Get amplitude and time of the time sample with the maximum amplitude
  /// \param[out] amplitude Amplitude of the maximum
  /// \param[out] time Time of the maximum
  /// \return True if the digit contains at least one time sample
  Bool_t  GetMaximum(Int_t& amplitude, Int_t& time) const;
  
  /// \brief Print digit stored info
  ///
  /// Refer to PrintStream for the implementation
  virtual void Print(const Option_t* opt) const;

  /// \brief Print information of this digit to the output stream
  /// \param stream Stream where the print the digit
  virtual void PrintStream(std::ostream &stream) const;
  
protected: 
  
  AliEMCALRawDigit           (const AliEMCALRawDigit &cd); // Not implemented
  AliEMCALRawDigit &operator=(const AliEMCALRawDigit &cd); // Not implemented
  
  Int_t   fId;            ///< Module absolute ID
  
  Int_t   fNSamples;      ///< Number of time samples
  
  /// Data field for time samples 
  Int_t*  fSamples;	      //[fNSamples]
  
  Float_t fAmplitude;     ///< digit amplitude
  
  Float_t fTime;          ///< digit time 
  
  ClassDef(AliEMCALRawDigit,1)  ;
};

/// \brief Output stream operator for AliEMCALRawDigit
/// \param stream Stream where to print the digit on
/// \param dig Digit to be printed
/// \return Stream after printing the digit
std::ostream &operator<<(std::ostream &stream, const AliEMCALRawDigit &dig);

#endif // ALIEMCALRAWDIGIT_H

