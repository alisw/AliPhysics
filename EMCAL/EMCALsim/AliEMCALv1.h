#ifndef ALIEMCALV1_H
#define ALIEMCALV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
/// \class AliEMCALv1
/// \ingroup EMCALsim
/// \brief EMCal simulation manager class v1
///
/// Implementation version v1 of EMCAL Manager class; 
/// An object of this class does not produce digits.
/// It is the one to use if you do want to produce outputs in TREEH 
///
/// This Class not stores information on all particles prior to EMCAL entry - 
/// in order to facilitate analysis. This is done by setting fIShunt =2, and 
/// flagging all parents of particles entering the EMCAL.
///
/// WARNING: Do not use for full detector simulations, use v2.
///
/// \author Sahal Yacoob (LBL /UCT)
/// \author Jennifer Klay (LBL)
/// \author Yves Schutz (IN2P3)
//_________________________________________________________________________


// --- ROOT system ---
class TClonesArray;
class TLorentzVector;
class TFile;

// --- AliRoot header files ---
#include "AliEMCALv0.h"

class AliEMCALv1 : public AliEMCALv0 
{
  
public:
  
  AliEMCALv1(void) ; 
  AliEMCALv1(const char *name, const char *title="", const Bool_t checkGeoAndRun = kTRUE) ;
  virtual ~AliEMCALv1(void) ;
  
  using AliEMCALv0::AddHit;
  virtual void  AddHit( Int_t shunt, Int_t primary, Int_t track, Int_t iparent, Float_t ienergy,
                        Int_t id, Float_t *hits, Float_t *p ) ;
  
  // Gives the version number 
  virtual       Int_t  IsVersion(void) const { return 1 ; }
  virtual const TString  Version(void) const { return TString("v1") ; }

  virtual void StepManager(void) ;
  virtual void RemapTrackHitIDs(Int_t *map);
  virtual void FinishPrimary();
  
  virtual void    SetTimeCut(Float_t tc) { fTimeCut = tc   ; }
  virtual Float_t GetTimeCut() const     { return fTimeCut ; } 
  
protected:
  
  Int_t   fCurPrimary;  ///< Current primary track
  Int_t   fCurParent;   ///< Current parent 
  Int_t   fCurTrack;    ///< Current track
  Float_t fTimeCut;     ///< Cut to remove the background from the ALICE system
  
private:
  
  AliEMCALv1              (const AliEMCALv1 & emcal);
  AliEMCALv1 & operator = (const AliEMCALv1 & /*rvalue*/);
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALv1,9) ;
  /// \endcond

};

#endif // AliEMCALV1_H
