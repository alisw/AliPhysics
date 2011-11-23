//-*- Mode: C++ -*-

// $Id$

#ifndef ALIHLTESDTRACKCUTS_H
#define ALIHLTESDTRACKCUTS_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTESDTrackCuts.h
/// @author Jochen Thaeder <jochen@thaeder.de>
/// @brief  ESD track cuts used in the analysis of HLT data
///

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliESDtrackCuts.h"
#include "AliHLTLogging.h"

/**
 * @class AliHLTESDTrackCuts
 * This class is an HLT wrapper for ESD track cuts used in the analysis
 * 
 * It provides a method function to the 2010 standard pp cuts.
 */

class AliHLTESDTrackCuts : public AliESDtrackCuts, public AliHLTLogging {
public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** Constructor */
  AliHLTESDTrackCuts(const Char_t* name = "AliHLTESDTrackCuts", const Char_t* title = "");

  /** Destructor */
  ~AliHLTESDTrackCuts();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Selection
   * ---------------------------------------------------------------------------------
   */

  /** Selection of esd track 
   *  @param  obj  ptr to AliESDtrack
   *  @return      kTRUE if track survives the cuts
   */
  virtual Bool_t IsSelected(TObject* obj);

  /*
   * ---------------------------------------------------------------------------------
   *                           Standard Track Cut Definitions
   * ---------------------------------------------------------------------------------
   */

  /** Get standard ESD track cuts used for 2010 pp data analysis
   *  Important : Returned object has to be deleted by user !
   *
   *  !!! Be aware - this is not the final yet
   *
   *  @return new AliHLTESDTrackCuts object, to be deleted by user
   */
  static AliHLTESDTrackCuts* GetStandardTrackCuts2010pp();

  ///////////////////////////////////////////////////////////////////////////////////
  
 private:
 
  /** copy constructor prohibited */
  AliHLTESDTrackCuts(const AliHLTESDTrackCuts&);
  
  /** assignment operator prohibited */
  AliHLTESDTrackCuts& operator=(const AliHLTESDTrackCuts&);
  
  ClassDef(AliHLTESDTrackCuts, 0)
};
#endif
