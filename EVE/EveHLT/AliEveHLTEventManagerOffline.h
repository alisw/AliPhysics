#ifndef ALIEVEHLTEVENTMANAGEROFFLINE_H
#define ALIEVEHLTEVENTMANAGEROFFLINE_H

class AliESDEvent;

#include "AliEveHLTEventManager.h"
#include "AliEveEventBufferOffline.h"

class AliEveHLTEventManagerOffline : public AliEveHLTEventManager { 

public:

  ///Constructor
  AliEveHLTEventManagerOffline(TString filename);
  
  virtual ~AliEveHLTEventManagerOffline();

  void NextEvent();
  void NavigateFwd();
  void NavigateBack();

 private:

  ///Default constructor, private
  AliEveHLTEventManagerOffline();

  /** copy constructor prohibited */
  AliEveHLTEventManagerOffline(const AliEveHLTEventManagerOffline&);

  /** assignment operator prohibited */
  AliEveHLTEventManagerOffline& operator=(const AliEveHLTEventManagerOffline&);

  /** Process the event data */
  //Int_t ProcessEvent(AliESDEvent * event);
  AliEveEventBufferOffline * fEventBuffer;
  AliEveEventBuffer * GetEventBuffer() {return fEventBuffer;}

  ClassDef(AliEveHLTEventManagerOffline, 1);

};

#endif
