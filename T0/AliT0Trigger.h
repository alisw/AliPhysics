#ifndef ALIT0TRIGGER_H
#define ALIT0TRIGGER_H

/// \ingroup sim
/// \class AliT0Trigger
/// \brief T0 trigger class
///
/////////////////////////////////////////////////
///  T0 Trigger Detector Class               //
/////////////////////////////////////////////////

#include "AliTriggerDetector.h"
class AliT0;

class AliT0Trigger : public AliTriggerDetector
{
 public:
  AliT0Trigger();  // constructor
  virtual ~AliT0Trigger(){}  // destructor
  virtual void    CreateInputs();
  virtual void    Trigger();
  
 private:
  
  AliT0 *fT0;          //!
  AliT0digit *fDigits   ; //! digits

  AliT0Trigger(const AliT0Trigger&);
  AliT0Trigger& operator=(const AliT0Trigger&);

  
  ClassDef(AliT0Trigger,1)  // T0 Trigger Detector class
};

typedef AliT0Trigger AliSTARTTrigger; // for backward compatibility

#endif








