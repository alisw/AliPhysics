#ifndef ALIEMCALV2_H
#define ALIEMCALV2_H
/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
/// \class AliEMCALv2
/// \ingroup EMCALsim
/// \brief EMCal simulation manager class v2
///
/// Implementation version v2 of EMCAL Manager class; SHASHLYK version
/// An object of this class does not produce digits.
/// It is the one to use if you do want to produce outputs in TREEH 
///
/// This Class does not store information on all particles prior to EMCAL entry - 
/// in order to facilitate analysis. This is done by setting fIShunt =2, and 
/// flagging all parents of particles entering the EMCAL.
///
/// \author Alexei Pavlinov (WSU)
/// \author M.L. Wang CCNU Wuhan & Subatech Oct-23-2012: Adapted for DCAL 
///
//_________________________________________________________________________

// --- ROOT system ---
class TBrowser;

// --- AliRoot header files ---
#include "AliEMCALv1.h"

class AliEMCALv2 : public AliEMCALv1 
{
  
public:

  AliEMCALv2(void) ; 
  AliEMCALv2(const char *name, const char *title="", const Bool_t checkGeoAndRun = kTRUE) ;
  virtual ~AliEMCALv2(void) ;

  using AliEMCALv1::AddHit;
  virtual void  AddHit( Int_t shunt, Int_t primary, Int_t track, Int_t iparent, Float_t ienergy,
                        Int_t id, Float_t *hits, Float_t *p ) ;

  virtual void StepManager(void) ;

  // Gives the version number 
  virtual       Int_t IsVersion(void) const { return 2 ; }
  virtual const TString Version(void) const { return TString("v2") ; }

 private:

  AliEMCALv2              (const AliEMCALv2 & emcal);
  AliEMCALv2 & operator = (const AliEMCALv2 & /*rvalue*/);
 
  /// \cond CLASSIMP
  ClassDef(AliEMCALv2,3) ;
  /// \endcond

};

#endif // AliEMCALV2_H
