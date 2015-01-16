//-*- Mode: C++ -*-

// $Id: AliHLTJETConeJetCandidate.h $

#ifndef ALIHLTJETCONEJETCANDIDATE_H
#define ALIHLTJETCONEJETCANDIDATE_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETConeJetCandidate.h
    @author Jochen Thaeder
    @date   
    @brief  Jet candidate of the cone finder
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "TClonesArray.h"

#include "AliHLTLogging.h"
#include "AliHLTJETBase.h"

class AliHLTJETConeEtaPhiCell;

/**
 * @class AliHLTJETConeJetCandidate
 * This class is contstructed with a seed and contains a found jet 
 * from the jet cone finder
 * Two options exist: 
 * <ul>
 *    <li>Add up the whole cell to the jet candidate<li>
 *    <li>Add up only patricles inside the coneradius to the jet candidate<li>
 * </ul>
 *
 * @ingroup alihlt_jet_cone
 */

class AliHLTJETConeJetCandidate : public TObject, public AliHLTLogging {
public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** Constructor 
   *  @param aEtaPhi   eta and phi of the seed
   *  @param aGridIdx  indeces in the grid
   *  @param useWholeCell  XXXXX
  */
  AliHLTJETConeJetCandidate( const Float_t* aEtaPhi, const Int_t* aGridIdx, 
			     Float_t coneRadius, Bool_t useWholeCell = kTRUE);

  /** Destructor */
  ~AliHLTJETConeJetCandidate();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */

  // -- Sees properties
  
  /** Get cell eta index of seed*/
  Int_t         GetSeedEtaIdx()    { return fSeedEtaIdx; }

  /** Get cell phi index of seed*/
  Int_t         GetSeedPhiIdx()    { return fSeedPhiIdx; }

  /** Get pt of seed */
  Float_t       GetSeedPt()        { return fSeedPt; }  
  
  // -- Jet properties

  /** Get eta of jet */
  Float_t       GetEta()           { return ( (!fNTracks) ? 0. : fEta / fNTracks); } 

  /** Get phi of jet  */
  Float_t       GetPhi()           { return ( (!fNTracks) ? 0. : fPhi / fNTracks); } 

  /** Get pt of jet */
  Float_t       GetPt()            { return fPt; }  

  /** Get Et of jet */
  Float_t       GetEt()            { return fPt; }  

  /*
   * ---------------------------------------------------------------------------------
   *                                     Process 
   * ---------------------------------------------------------------------------------
   */

  /** Add cell to JetCandidate 
   *  @param cell ptr to cell
   *  @return 0 on success, <0 on failure
   */
  Int_t AddCell( AliHLTJETConeEtaPhiCell* cell );



  /* XXXXXXXXXX
  void SetAll( Float_t pt, Float_t eta, Float_t phi, Int_t nTracks, Bool_t useWholeCell) {
    fPt = pt; 
    fEta = eta*nTracks; 
    fPhi = phi*nTracks; 
    fNTracks = nTracks; 
    fUseWholeCell = useWholeCell; 
  }
  */

  /*
   * ---------------------------------------------------------------------------------
   *                                Sort of JetCandidates 
   * ---------------------------------------------------------------------------------
   */

  /** Compare this class with an other instance of this class
   *  used in a TClonesArray::Sort() 
   *  @param   obj  ptr to other instance  
   *  @return  Returns 0 when equal, 1 when this is smaller 
   *  and -1 when bigger -- sorts descending
   */
  Int_t Compare( const TObject* obj) const;


  /** Defines this class as being sortable in a TClonesArray
   *  @return     always kTRUE;
   */
  Bool_t IsSortable() const  { return kTRUE; }

  ///////////////////////////////////////////////////////////////////////////////////
  
 private:
 
  /** standard constructor prohibited */
  AliHLTJETConeJetCandidate();

  /** copy constructor prohibited */
  AliHLTJETConeJetCandidate(const AliHLTJETConeJetCandidate&);

  /** assignment operator prohibited */
  AliHLTJETConeJetCandidate& operator=(const AliHLTJETConeJetCandidate&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Helper - private
   * ---------------------------------------------------------------------------------
   */

  /** Get distance squared between two points in eta-phi
   *  @param eta1    eta coordinate of first point
   *  @param phi1    phi coordinate of first point
   *  @param eta2    eta coordinate of second point
   *  @param phi2    phi coordinate of second point
   *  @return        Distance squared
   */
  Float_t GetDistance2( const Float_t eta1, const Float_t phi1, 
			const Float_t eta2, const Float_t phi2);
  
  /** Check if particle is in side the cne
   *  @param eta    eta coordinate 
   *  @param phi    phi coordinate 
   *  @return       kTRUE if it is, otherwise kFALSE
   */
  Bool_t InCone( Float_t eta, Float_t phi );

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  //  -- Seed cell index parameters

  /** Seed cell index */
  Int_t                  fSeedCellIdx;      // see above

  /** Seed cell eta index */ 
  Int_t                  fSeedEtaIdx;       // see above

  /** Seeed cell phi index */  
  Int_t                  fSeedPhiIdx;       // see above

  // -- Seed parameters

  /** seed eta */
  Float_t                fSeedEta;          // see above

  /** seed phi */
  Float_t                fSeedPhi;          // see above

  /** seed pt */
  Float_t                fSeedPt;           // see above

  // -- Summed variables

  /** Summed eta */
  Float_t                fEta;              // see above

  /** Summed phi */
  Float_t                fPhi;              // see above

  /** Summed pt */
  Float_t                fPt;               // see above

  /** Number of tracks in JetCandidate */
  UInt_t                 fNTracks;          // see above

  /** Flag if whole cell should be added, or every track */
  Bool_t                 fUseWholeCell;     // see above

  /** Cone radius squared */
  Float_t                fConeRadius2;      // see above

  ClassDef(AliHLTJETConeJetCandidate, 1)
};
#endif
