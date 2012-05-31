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

#include "AliHLTJETTrackCuts.h"
#include "AliHLTJETConeSeedCuts.h"

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
   *                                   Initialize
   * ---------------------------------------------------------------------------------
   */
  
  /** Initialize reader header for cone jet finder
   *  @return 0 on success, otherwise <0
   */
  Int_t Initialize();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Setter
   * ---------------------------------------------------------------------------------
   */
  
  /** Set Analysis Cuts
   *  @param cuts ptr to AliHLTJETTrackCuts 
   */
  void SetTrackCuts( AliHLTJETTrackCuts* cuts )  { fTrackCuts = cuts; }
 
  /** Get Analysis Cuts
   *  @param cuts ptr to AliHLTJETConeSeedCuts 
   */
  void SetSeedCuts( AliHLTJETConeSeedCuts* cuts ){ fSeedCuts = cuts; }

  /** Set the cone radius */
  void SetConeRadius( Float_t f ) { fConeRadius = f; }

  /** Set grid binning in eta */
  void SetGridEtaBinning( Float_t f ) { fGridEtaBinning = f; }

  /** Set grid binning in phi */
  void SetGridPhiBinning( Float_t f ) { fGridPhiBinning = f; }

  /** Set algorithm type */
  void SetJetAlgorithm( AliHLTJETBase::JetAlgorithmType_t a ) { fAlgorithm = a; }

  /** Set Usage of Kinematics */
  void SetUseMC( Bool_t b ) { fUseMC = b; }

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */

  /** Get Analysis Cuts
   *  @return ptr to AliHLTJETTrackCuts 
   */
  AliHLTJETTrackCuts* GetTrackCuts()     { return fTrackCuts; }
 
  /** Get Analysis Cuts
   *  @return ptr to AliHLTJETConeSeedCuts
   */
  AliHLTJETConeSeedCuts* GetSeedCuts()   { return fSeedCuts; }
  
  /** Get grid eta range */
  Float_t GetGridEtaRange()              { return fGridEtaRange; }

  /** Get grid phi range */
  Float_t GetGridPhiRange()              { return fGridPhiRange; }

  /** Get grid eta binning */
  Float_t GetGridEtaBinning()            { return fGridEtaBinning; }
  
  /** Get grid phi binning */
  Float_t GetGridPhiBinning()            { return fGridPhiBinning; }

  /** Get cone radius */
  Float_t GetConeRadius()                { return fConeRadius; }

  /** Get algorithm type */
  AliHLTJETBase::JetAlgorithmType_t GetJetAlgorithm() { return fAlgorithm; }

  /** Get Usage of Kinematics */
  Bool_t  GetUseMC()                     { return fUseMC; }

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

  /** Cuts on track selection */
  AliHLTJETTrackCuts        *fTrackCuts;            //! transient

  /** Cuts on seed selection */
  AliHLTJETConeSeedCuts     *fSeedCuts;             //! transient

  // -- Grid members

  /** Grid binning in eta */
  Float_t                    fGridEtaBinning;       // see above

  /** Grid binning in phi */
  Float_t                    fGridPhiBinning;       // see above

  // -- Grid binning 

  /** Grid eta range */
  Float_t                    fGridEtaRange;         // see above

  /** Grid phi range 
   *  abs(phimin)+phiMax + 2*coneRadius;
   */
  Float_t                    fGridPhiRange;         // see above

  // -- Algorithm members

  /** Algorithm */
  AliHLTJETBase::JetAlgorithmType_t fAlgorithm;     // see above

  /** Cone radius */
  Float_t                    fConeRadius;           // see above
  
  /** Use MC Data -- only needed for off-line */
  Bool_t                     fUseMC;                // see above 

  ClassDef(AliHLTJETReaderHeader, 1)

};
#endif

