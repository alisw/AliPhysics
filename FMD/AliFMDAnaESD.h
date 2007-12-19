// -*- mode: C++ -*-
#ifndef ALIFMDANAESD_H
#define ALIFMDANAESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
//____________________________________________________________________
//
// Utility class for analysing ESD data. 
// This class does sharing and background correction 
//
#include <AliFMDInput.h>
class AliFMDAnaRing;
// #include "AliFMDAnaRing.h"
class TBrowser;

/** Base class for analysing FMD ESD data */
class AliFMDAnaESD : public AliFMDInput 
{
public:
  /** Constructor 
      @param lower cut. */
  AliFMDAnaESD();
  /** Destructor */
  virtual ~AliFMDAnaESD() {}
  /** Called at beginning of run 
      @return @c false on error */
  virtual Bool_t Init();
  /** Begining of event
      @param ev Event number
      @return @c false on error */
  virtual Bool_t Begin(Int_t ev);
  /** Loop over all ESD data, and call ProcessESD for each entry.
      @return  @c false on error  */
  virtual Bool_t ProcessESDs();
  /** Fill for analysis */
  virtual void Fill(Float_t phi, Float_t eta, Float_t mult) {}
  /** Called at the end of event
      @return @c false in case of errors */
  virtual Bool_t End();
  /** Called at the end of run 
      @return @c false in case of errors */
  virtual Bool_t Finish();
  /** Browse this object */
  virtual void Browse(TBrowser* b);
  /** This is a folder */ 
  virtual Bool_t IsFolder() const { return kTRUE; }
protected:
  /** Add a ring 
      @param ring Ring object */
  virtual void AddRing(AliFMDAnaRing* ring);
  /** Find a ring index
      @param det  Detector number 
      @param ring Ring id
      @return Index of ring object */
  virtual Int_t FindRing(UShort_t det, Char_t ring) const;
  /** Ring objects */ 
  AliFMDAnaRing* fRing[5];
  /** Number of events */
  ULong_t fNEvents;
  ClassDef(AliFMDAnaESD,0) // Base class for analysing FMD ESD
};

#endif
//
// EOF
//

  
