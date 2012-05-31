//-*- Mode: C++ -*-

// $Id$


/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice     
 */

/** @file   AliEveEventBufferOffline.h
    @author Svein Lindal
    @date
    @brief  Manager for HOMER in aliroot
*/

#define BUFFERSIZE 15


#ifndef ALIEVEEVENTBUFFEROFFLINE_H
#define ALIEVEEVENTBUFFEROFFLINE_H

#include "Rtypes.h"
#include "AliEveEventBuffer.h"

class TFile;
class TTree;
class AliESDEvent;
class TString;

class AliEveEventBufferOffline : public AliEveEventBuffer {


public:
  
  /** default constructor */
  AliEveEventBufferOffline(TString file);
  /** destructor */
  virtual ~AliEveEventBufferOffline();

  void ConnectToSource();
  void WriteToFile(Int_t runnumber);

private:

  //not allowed
  AliEveEventBufferOffline();

  /** copy constructor prohibited */
  AliEveEventBufferOffline(const AliEveEventBufferOffline&);

  /** assignment operator prohibited */
  AliEveEventBufferOffline& operator=(const AliEveEventBufferOffline&);

  ///Inherited from AliEveEventBuffer
  TObject * GetEventFromSource();

  TFile * fFile;  //File poineter
  Int_t fNEntries; //Number of entries
  Int_t fEventNo; //Event number
  AliESDEvent * fEvent; //Event pointer
  TTree * fTree; //TTree pointer

  ClassDef(AliEveEventBufferOffline, 0); 
};

#endif
