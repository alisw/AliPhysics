//-*- Mode: C++ -*-

// $Id: AliHLTJETConeEtaPhiCell.h  $

#ifndef ALIHLTJETCONEETAPHICELL_H
#define ALIHLTJETCONEETAPHICELL_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETConeEtaPhiCell.h
    @author Jochen Thaeder
    @date   
    @brief  Cell in eta-phi space of the cone finder
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "TObjArray.h"
#include "TParticle.h"

#include "AliESDtrack.h"

#include "AliHLTLogging.h"
#include "AliHLTJETBase.h"

/**
 * @class AliHLTJETConeEtaPhiCell
 * This class contains one cell in the eta-phi space for the cone finder
 * 
 * @ingroup alihlt_jet_cone
 */

class AliHLTJETConeEtaPhiCell : public TObject, public AliHLTLogging  {

public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */
  
  /** Constructor for ESD tracks */
  AliHLTJETConeEtaPhiCell( Int_t etaIdx, Int_t phiIdx, AliESDtrack* track );

  /** Constructor for MC particles */
  AliHLTJETConeEtaPhiCell( Int_t etaIdx, Int_t phiIdx, TParticle* particle );

  /** Destructor */
  ~AliHLTJETConeEtaPhiCell();

  /** A destructor like class, called by TClonesArray->Clear("C") */
  void Clear(Option_t* option = "");

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */
 
  /** Get cell eta index */
  Int_t         GetEtaIdx()        { return fEtaIdx; }

  /** Get cell phi index */
  Int_t         GetPhiIdx()        { return fPhiIdx; }

  /** Get summed pt of cell */
  Float_t       GetPt()            { return fPt; }  

  /** Get summed eta of cell */
  Float_t       GetEta()           { return fEta; }

  /** Get summed phi of cell */
  Float_t       GetPhi()           { return fPhi; }  

  /** Get N of tracks in cell */
  Int_t         GetNTracks()       { return fNTracks; }

  /**  Get List of tracks in cell */
  TObjArray*    GetTrackList()     { return fTrackList; }

  /** Get type of tracks (TrackType_t) */
  Int_t         GetTrackType()     { return fTrackType; }


  /*
   * ---------------------------------------------------------------------------------
   *                                     Process 
   * ---------------------------------------------------------------------------------
   */

  /** Add Track to cell, using ESD track */
  void AddTrack( AliESDtrack* track );

  /** Add Track to cell, using MC particle */
  void AddTrack( TParticle* particle );
  
  ///////////////////////////////////////////////////////////////////////////////////

 private:

  /** Standard constructor prohibited */
  AliHLTJETConeEtaPhiCell();

  /** copy constructor prohibited */
  AliHLTJETConeEtaPhiCell(const AliHLTJETConeEtaPhiCell&);

  /** assignment operator prohibited */
  AliHLTJETConeEtaPhiCell& operator=(const AliHLTJETConeEtaPhiCell&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  // -- cell index parameters --

  /** cell eta index */
  Int_t      fEtaIdx;                        // see above

  /** cell phi index */
  Int_t      fPhiIdx;                        // see above

  // -- summed parameters --

  /** Summed pt of cell */
  Float_t    fPt;                            // see above

  /** Summed eta of cell */
  Float_t    fEta;                           // see above

  /** Summed phi of cell */
  Float_t    fPhi;                           // see above

  /** N of tracks in cell */
  UInt_t     fNTracks;                       // see above

  // -- Lists --

  /** TObjArray of tracks in the cell */
  TObjArray* fTrackList;                     //! transient
  
  /** Type of tracks (TrackType_t) */
  Int_t      fTrackType;                     // see above

  ClassDef(AliHLTJETConeEtaPhiCell, 1)
};
#endif

