#ifndef ALIJETREADER_H
#define ALIJETREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
 
// Jet reader base class
// manages the reading of input for jet algorithms
// Authors: jgcn@mda.cinvestav.mx
//          Magali Estienne <magali.estienne@subatech.in2p3.fr>  
//          alexandre.shabetai@cern.ch

#include <TObject.h>

class AliJetReaderHeader;
class AliJetFillCalTrkEvent;
class AliJetCalTrkEvent;
class AliVEvent;
class AliMCEvent;

class AliJetReader : public TObject 
{
 public: 
  AliJetReader();
  virtual ~AliJetReader();

  // Getters
  AliJetCalTrkEvent*        GetCalTrkEvent()       const {return fCalTrkEvent;}
  AliJetReaderHeader*       GetReaderHeader()      const {return fReaderHeader;}
  
  // Setters
  void                      SetReaderHeader(AliJetReaderHeader* header)      {fReaderHeader = header;}  

  // Others
  void                      SetInputEvent(const TObject* esd, const TObject* aod, const AliMCEvent* mc);
  void                      InitTasks();
  Bool_t                    CreateTasks();
  Bool_t                    ExecTasks();
  Bool_t                    ProcessEvent();
  void                      WriteRHeaderToFile() const;
  void                      WriteReaderHeader();
 
 protected:
  AliJetReader(const AliJetReader& rJetReader);
  AliJetReader& operator = (const AliJetReader& rhsr);

  AliJetCalTrkEvent*        fCalTrkEvent;                     //! Pointer to calTrkEvent
  AliJetFillCalTrkEvent*    fFillEvent;                       //! Pointer to AliJetFillCalTrkEvent
  AliJetReaderHeader*       fReaderHeader;                    //  Pointer to header
  AliJetFillCalTrkEvent*    fFillEventwTrks;                  //  For charged particle task
  Int_t                     fDebug;                           //  Debug option
  AliVEvent*                fVEvent;                          //! Input event
  AliMCEvent*               fMCEvent;			      //! MC Event;
  Int_t                     fOpt;                             //  Detector config
  
  ClassDef(AliJetReader,2)                                    // jet reader class

};
 
#endif
