//-*- Mode: C++ -*-

// $Id:  $

#ifndef ALIHLTJETREADERHEADER_H
#define ALIHLTJETREADERHEADER_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETReaderHeader.h
    @author Jochen Thaeder
    @date   
    @brief  ReaderHeader for jet finder
*/

#include "AliAnalysisCuts.h"

#include "AliJetReaderHeader.h"
#include "AliHLTLogging.h"

/**
 * @class  AliHLTJETReaderHeader
 * ReaderHeader for jet finder
 *
 * @ingroup alihlt_jet
 */

class AliHLTJETReaderHeader : public AliJetReaderHeader, public AliHLTLogging {
  
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** standard constructor */
  AliHLTJETReaderHeader();

  /** destructor */
  virtual ~AliHLTJETReaderHeader();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Setter
   * ---------------------------------------------------------------------------------
   */
  
  /** Set Analysis Cuts
   *  @param cuts ptr to AnalysisCuts
   */
  void SetAnalysisCuts( AliAnalysisCuts* cuts ) { fCuts = cuts; }

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */
  
  /** Get Analysis Cuts
   *  @return ptr to AnalysisCuts
   */
   AliAnalysisCuts* GetAnalysisCuts() { return fCuts; }


  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTJETReaderHeader (const AliHLTJETReaderHeader&);

  /** assignment operator prohibited */
  AliHLTJETReaderHeader& operator= (const AliHLTJETReaderHeader&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** Minimum pt  */
  AliAnalysisCuts           *fCuts;                   // see above

  ClassDef(AliHLTJETReaderHeader, 1)

};
#endif

