#ifndef ALIITSV0FINDER_H
#define ALIITSV0FINDER_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------------------------
//                V0 finder on-the-fly during ITS tracking
//           Origin: Marian Ivanov, CERN, Marian.Ivanov@cern.ch
//           Extraction to a separate class: Andrea Dainese
//           Current support and development: 
//-------------------------------------------------------------------------

/* $Id$ */

class TTree;
class TTreeSRedirector;
class AliESDEvent;
class AliITStrackerMI;

//-------------------------------------------------------------------------
class AliITSV0Finder : public TObject {
public:
  AliITSV0Finder();
  //AliITSV0Finder(const AliITSV0Finder &/*v0Finder*/) {;}
  //AliITSV0Finder & operator=(const AliITSV0Finder &/*v0Finder*/) {;}
 
  virtual ~AliITSV0Finder();

  //try to find V0
  static void FindV02(AliESDEvent *event,AliITStrackerMI *tracker);  
  //try to refit  V0's
  static void RefitV02(const AliESDEvent *event,AliITStrackerMI *tracker);
  //try to update, or reject TPC  V0s
  static void UpdateTPCV0(const AliESDEvent *event,AliITStrackerMI *tracker);  

  TTreeSRedirector *GetDebugStreamer() {return fDebugStreamer;}

  static void   SetDisabled(Bool_t v=kTRUE) {fgDisabled = v;}
  static Bool_t GetDisabled()               {return fgDisabled;}

 protected:
  static Bool_t fgDisabled;     // possibilidy to disable from reconstruction in cpases

private:
  TTreeSRedirector *fDebugStreamer;      //!debug streamer
 

  AliITSV0Finder(const AliITSV0Finder& obj);
  AliITSV0Finder& operator=(const AliITSV0Finder& obj);
 
  ClassDef(AliITSV0Finder,0)   // on-the-fly V0 finder for AliITStrackerMI
};
#endif
