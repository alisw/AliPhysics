// Author: 2010 Svein Lindal <slindal@fys.uio.no>                        *
//         for The ALICE HLT Project.                                    *

#ifndef ALIEVEHLTEVENTMANAGERHOMER_H
#define ALIEVEHLTEVENTMANAGERHOMER_H

class AliESDEvent;

#include "AliEveHLTEventManager.h" 
#include "AliEveEventBufferHomer.h"
#include "AliEveEventBuffer.h"
class TList;
class TTimer;
class TGLOverlayButton;

class AliEveHLTEventManagerHomer : public AliEveHLTEventManager { 

public:

  ///Constructor
  AliEveHLTEventManagerHomer();
  
  virtual ~AliEveHLTEventManagerHomer();

  //Get Next Event
  void NextEvent();
  //Try to get the next event
  void TryNextEvent();
  //Get next event in buffer
  void NavigateFwd();
  //Get Previous event in buffer
  void NavigateBack();

  //Process block list
  void ProcessList(TList * blockList);

 private:


  /** copy constructor prohibited */
  AliEveHLTEventManagerHomer(const AliEveHLTEventManagerHomer&);

  /** assignment operator prohibited */
  AliEveHLTEventManagerHomer& operator=(const AliEveHLTEventManagerHomer&);
  
  AliEveEventBufferHomer * fEventBuffer; //Event buffer
  ///Get event buffer
  AliEveEventBuffer * GetEventBuffer() { return dynamic_cast<AliEveEventBuffer*>(fEventBuffer); }

  
  TTimer * fNextEventTimer;  //Timer to fetch next event
  TGLOverlayButton * fInfoButton; //Information button


  ClassDef(AliEveHLTEventManagerHomer, 0);



};

#endif
