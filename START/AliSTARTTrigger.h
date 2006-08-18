#ifndef ALISTARTTRIGGER_H
#define ALISTARTTRIGGER_H

/// \ingroup sim
/// \class AliSTARTTrigger
/// \brief START trigger class
///
/////////////////////////////////////////////////
///  START Trigger Detector Class               //
/////////////////////////////////////////////////

#include "AliTriggerDetector.h"
class AliSTART;

class AliSTARTTrigger : public AliTriggerDetector
{
 public:
  AliSTARTTrigger();  // constructor
  virtual ~AliSTARTTrigger(){}  // destructor
  virtual void    CreateInputs();
  virtual void    Trigger();
  
 private:
  
  AliSTART *fSTART;          //!
  AliSTARTdigit *fDigits   ; //! digits

  AliSTARTTrigger(const AliSTARTTrigger&);
  AliSTARTTrigger& operator=(const AliSTARTTrigger&);

  
  ClassDef(AliSTARTTrigger,1)  // START Trigger Detector class
    };
#endif








