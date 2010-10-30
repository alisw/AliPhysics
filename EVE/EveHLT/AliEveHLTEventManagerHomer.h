#ifndef ALIEVEHLTEVENTMANAGERHOMER_H
#define ALIEVEHLTEVENTMANAGERHOMER_H

class AliESDEvent;

class AliEveHLTEventManager;
#include "AliEveEventBufferHomer.h"
#include "AliEveEventBuffer.h"
class TList;

class AliEveHLTEventManagerHomer : public AliEveHLTEventManager { 

public:

  ///Constructor
  AliEveHLTEventManagerHomer();
  
  virtual ~AliEveHLTEventManagerHomer();

  void NextEvent();
  void NavigateFwd();
  void NavigateBack();

 private:


  /** copy constructor prohibited */
  AliEveHLTEventManagerHomer(const AliEveHLTEventManagerHomer&);

  /** assignment operator prohibited */
  AliEveHLTEventManagerHomer& operator=(const AliEveHLTEventManagerHomer&);
  
  AliEveEventBufferHomer * fEventBuffer;
  AliEveEventBuffer * GetEventBuffer() { return dynamic_cast<AliEveEventBuffer*>(fEventBuffer); }

  ClassDef(AliEveHLTEventManagerHomer, 0);

};

#endif
